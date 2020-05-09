import numpy as np
from aiida.parsers import Parser
from aiida.orm import Dict
#from six.moves import range
from aiida.common import OutputParsingError
from aiida.common import exceptions

# See the LICENSE.txt and AUTHORS.txt files.

# TO DO Get modules metadata from setup script.

# List of scalar values from CML to be transferred to AiiDA
#pylint: disable=invalid-name
standard_output_list = [
    'siesta:FreeE', 'siesta:E_KS', 'siesta:Ebs', 'siesta:E_Fermi', 'siesta:stot'
]  ## leave svec for later

#####################################################
# BEGINNING OF SET OF AUXILIARY FUNCTIONS           #
# They should probably be put in a separate module  #
#####################################################


def text_to_array(text, dtype):
    return np.array(text.replace("\n", "").split(), dtype=dtype)


def get_parsed_xml_doc(xml_path):

    from xml.dom import minidom

    try:
        xmldoc = minidom.parse(xml_path)
    except EOFError:
        raise OutputParsingError("Faulty Xml File")
        #xmldoc = None

    # We might want to add some extra consistency checks

    return xmldoc


def get_dict_from_xml_doc(xmldoc):

    # Scalar items

    scalar_dict = {}

    # Metadata items
    itemlist = xmldoc.getElementsByTagName('metadata')
    for item in itemlist:
        #
        # Maybe make sure that 'name' does not contain
        # forbidden characters
        #
        name = item.attributes['name'].value
        value = item.attributes['content'].value
        scalar_dict[name] = value

    # Scalar output items
    # From the last "SCF Finalization" module
    # This means that we do not record non-converged values

    itemlist = xmldoc.getElementsByTagName('module')

    scf_final = None
    for item in itemlist:
        if 'title' in list(item.attributes.keys()):
            # Get last scf finalization module
            if item.attributes['title'].value == "SCF Finalization":
                scf_final = item

    if scf_final is not None:

        # wrapped in <property> elements with a <scalar> child
        props = scf_final.getElementsByTagName('property')

        for prop in props:
            if 'dictRef' in list(prop.attributes.keys()):
                name = prop.attributes['dictRef'].value
                if name in standard_output_list:
                    data = prop.getElementsByTagName('scalar')[0]
                    value = data.childNodes[0].nodeValue
                    units = data.attributes['units'].value
                    loc_colon = units.find(':')
                    unit_name = units[loc_colon + 1:]
                    loc_colon = name.find(':')
                    reduced_name = name[loc_colon + 1:]

                    # Put units in separate entries, as in QE
                    # Use numbers (floats) instead of strings

                    scalar_dict[reduced_name] = float(value)
                    scalar_dict[reduced_name + "_units"] = unit_name

    scalar_dict['variable_geometry'] = is_variable_geometry(xmldoc)
    #
    # Sizes of orbital set (and non-zero interactions), and mesh
    #
    no_u, nnz, mesh = get_sizes_info(xmldoc)
    if no_u is not None:
        scalar_dict['no_u'] = no_u
    if nnz is not None:
        scalar_dict['nnz'] = nnz
    if mesh is not None:
        scalar_dict['mesh'] = mesh

    return scalar_dict


def is_variable_geometry(xmldoc):
    """
     Tries to guess whether the calculation involves changes in
     geometry.
     """

    itemlist = xmldoc.getElementsByTagName('module')
    for item in itemlist:
        # Check there is a step which is a "geometry optimization" one
        if 'dictRef' in list(item.attributes.keys()):
            if item.attributes['dictRef'].value == "Geom. Optim":
                return True
            
    return False


def get_sizes_info(xmldoc):
    """
     Gets the number of orbitals and non-zero interactions
     """
    no_u = None
    nnz = None
    mesh = None

    itemlist = xmldoc.getElementsByTagName('module')
    for item in itemlist:
        # Process the first "step" module, which is a "geometry" one
        if 'serial' in list(item.attributes.keys()):
            # Get properties here
            props_list = item.getElementsByTagName('property')
            for prop in props_list:
                if 'dictRef' in list(prop.attributes.keys()):
                    if prop.attributes['dictRef'].value == "siesta:no_u":
                        scalar = prop.getElementsByTagName('scalar')[0]
                        no_u = int(scalar.childNodes[0].data)
                    if prop.attributes['dictRef'].value == "siesta:nnz":
                        scalar = prop.getElementsByTagName('scalar')[0]
                        nnz = int(scalar.childNodes[0].data)
                    if prop.attributes['dictRef'].value == "siesta:ntm":
                        array = prop.getElementsByTagName('array')[0]
                        mesh = [int(s) for s in array.childNodes[0].data.split()]

    return no_u, nnz, mesh


def get_last_structure(xmldoc, input_structure):

    itemlist = xmldoc.getElementsByTagName('module')
    #
    # Use the last "geometry" module, and not the
    # "Finalization" one.

    finalmodule = None
    for item in itemlist:
        # Get a "geometry" module by the criteria:
        if 'serial' in list(item.attributes.keys()):
            if 'dictRef' in list(item.attributes.keys()):
                if item.attributes['dictRef'].value != "SCF":
                    finalmodule = item

    # In case there is no appropriate data, fall back and
    # at least return the initial structure
    # (this should not be necessary, as the initial Geometry module
    # is opened very soon)
    if finalmodule is None:
        #self.logger.warning("Returning input structure in output_structure node")
        return input_structure

    atoms = finalmodule.getElementsByTagName('atom')
    cellvectors = finalmodule.getElementsByTagName('latticeVector')

    atomlist = []

    for atm in atoms:
        kind = atm.attributes['elementType'].value
        x = atm.attributes['x3'].value
        y = atm.attributes['y3'].value
        z = atm.attributes['z3'].value
        atomlist.append([kind, [float(x), float(y), float(z)]])

    cell = []
    for latt in cellvectors:
        data = latt.childNodes[0].data.split()
        cell.append([float(s) for s in data])

    # Generally it is better to pass the input structure
    # and reset the data, since site 'names' are not handled by
    # the CML file (at least not in Siesta versions <= 4.0)

    stru = input_structure.clone()
    stru.reset_cell(cell)
    new_pos = [atom[1] for atom in atomlist]
    stru.reset_sites_positions(new_pos)

    return stru


def get_final_forces_and_stress(xmldoc):
    #
    # Extracts final forces and stress as lists of lists...
    #
    itemlist = xmldoc.getElementsByTagName('module')

    # Note: In modern versions of Siesta, forces and stresses
    # are written in the "SCF Finalization" modules at the end
    # of each geometry step.

    # Search for the last one of those modules

    scf_final = None
    for item in itemlist:
        if 'title' in list(item.attributes.keys()):
            # Get last scf finalization module
            if item.attributes['title'].value == "SCF Finalization":
                scf_final = item

    forces = None
    stress = None

    if scf_final is not None:
        props = scf_final.getElementsByTagName('property')
        for prop in props:
            if 'dictRef' in list(prop.attributes.keys()):

                if prop.attributes['dictRef'].value == 'siesta:forces':
                    mat = prop.getElementsByTagName('matrix')[0]
                    # Get flat list and reshape as list of lists
                    # using info on rows and columns in CML file
                    rows = int(mat.attributes['rows'].value)
                    cols = int(mat.attributes['columns'].value)
                    f = mat.childNodes[0].data.split()  #pylint: disable=invalid-name
                    f = [float(x) for x in f]  #pylint: disable=invalid-name
                    forces = [f[rows * i:rows * (i + 1)] for i in range(cols)]

                if prop.attributes['dictRef'].value == 'siesta:stress':
                    mat = prop.getElementsByTagName('matrix')[0]
                    # Get flat list and reshape as list of lists
                    # using info on rows and columns in CML file
                    rows = int(mat.attributes['rows'].value)
                    cols = int(mat.attributes['columns'].value)
                    s = mat.childNodes[0].data.split()  #pylint: disable=invalid-name
                    s = [float(x) for x in s]  #pylint: disable=invalid-name
                    stress = [s[rows * i:rows * (i + 1)] for i in range(cols)]

    return forces, stress


##################################
# END OF AUXILIARY FUNCTIONS SET #
##################################


class SiestaOutputParsingError(OutputParsingError):
    pass


class SiestaCMLParsingError(OutputParsingError):
    pass


class SiestaParser(Parser):
    """
    Parser for the output of Siesta.
    """

    def parse(self, **kwargs):
        """
        Receives in input a dictionary of retrieved nodes. Does all the logic here.
        """
        from aiida.engine import ExitCode

        parser_version = '1.0.1'
        parser_info = {}
        parser_info['parser_info'] = 'AiiDA Siesta Parser V. {}'.format(parser_version)

        try:
            output_folder = self.retrieved
        except exceptions.NotExistent:
            raise OutputParsingError("Folder not retrieved")

        messages_path, xml_path, json_path, bands_path = \
            self._fetch_output_files(output_folder)

        if xml_path is None:
            raise OutputParsingError("Xml file not retrieved")
        xmldoc = get_parsed_xml_doc(xml_path)
        result_dict = get_dict_from_xml_doc(xmldoc)

        if json_path is not None:
            from .json_time import get_timing_info
            global_time, timing_decomp = get_timing_info(json_path)
            if global_time is None:
                raise OutputParsingError("Cannot fully parse the time.json file")
            else:
                result_dict["global_time"] = global_time
                result_dict["timing_decomposition"] = timing_decomp

        if messages_path is None:
            # Perhaps using an old version of Siesta
            warnings_list = ['WARNING: No MESSAGES file...']
        else:
            successful, warnings_list = self.get_warnings_from_file(messages_path)  #pylint: disable=unused-variable
        result_dict["warnings"] = warnings_list

        # Add parser info dictionary
        parsed_dict = dict(list(result_dict.items()) + list(parser_info.items()))

        output_data = Dict(dict=parsed_dict)
        self.out('output_parameters', output_data)

        # If the structure has changed, save it
        in_struc = self.node.inputs.structure
        if parsed_dict['variable_geometry']:
            # Get the input structure to copy its site names, as the CML file traditionally contained only the
            # atomic symbols.
            struc = get_last_structure(xmldoc, in_struc)
            self.out('output_structure', struc)

        forces, stress = get_final_forces_and_stress(xmldoc)
        if forces is not None and stress is not None:
            from aiida.orm import ArrayData
            arraydata = ArrayData()
            arraydata.set_array('forces', np.array(forces))
            arraydata.set_array('stress', np.array(stress))
            self.out('forces_and_stress', arraydata)

        for line in parsed_dict["warnings"]:
            if u'GEOM_NOT_CONV' in line:
                return self.exit_codes.GEOM_NOT_CONV
            if u'SCF_NOT_CONV' in line:
                return self.exit_codes.SCF_NOT_CONV

        if bands_path is None:
            if "bandskpoints" in self.node.inputs:
                return self.exit_codes.BANDS_FILE_NOT_PRODUCED

        # Parse band-structure information if available
        if bands_path is not None:
            bands, coords = self.get_bands(bands_path)
            from aiida.orm import BandsData
            arraybands = BandsData()
            f = self.node.inputs.bandskpoints  #Temporary workaround due to a bug
            #f._set_reciprocal_cell()         #in KpointData (issue #2749)
            arraybands.set_kpoints(f.get_kpoints(cartesian=True))
            arraybands.labels = f.labels
            arraybands.set_bands(bands, units="eV")
            self.out('bands', arraybands)
            bandsparameters = Dict(dict={"kp_coordinates": coords})
            self.out('bands_parameters', bandsparameters)

        return ExitCode(0)

    def _fetch_output_files(self, out_folder):
        """
        Checks the output folder for standard output and standard error files, returns their absolute paths
        or "None" in case the file is not found in the remote folder.
        """
        import os

        list_of_files = out_folder._repository.list_object_names()

        messages_path = None
        xml_path = None
        json_path = None
        bands_path = None

        namexmlfile = str(self.node.get_option('prefix')) + ".xml"
        if namexmlfile in list_of_files:
            xml_path = os.path.join(out_folder._repository._get_base_folder().abspath, namexmlfile)

        if self.node.process_class._JSON_FILE in list_of_files:
            json_path = os.path.join(
                out_folder._repository._get_base_folder().abspath, self.node.process_class._JSON_FILE
            )

        if self.node.process_class._MESSAGES_FILE in list_of_files:
            messages_path = os.path.join(
                out_folder._repository._get_base_folder().abspath, self.node.process_class._MESSAGES_FILE
            )

        namebandsfile = str(self.node.get_option('prefix')) + ".bands"
        if namebandsfile in list_of_files:
            bands_path = os.path.join(out_folder._repository._get_base_folder().abspath, namebandsfile)

        return messages_path, xml_path, json_path, bands_path

    def get_warnings_from_file(self, messages_path):
        """
     Generates a list of warnings from the 'MESSAGES' file, which  contains a line per message,
     prefixed with 'INFO', 'WARNING' or 'FATAL'.

     Returns a boolean indicating success (True) or failure (False) and a list of strings.
     """
        f = open(messages_path)
        lines = f.read().split('\n')  # There will be a final '' element

        import re

        # Search for 'FATAL:' messages, log them, and return immediately
        there_are_fatals = False
        for line in lines:
            if re.match('^FATAL:.*$', line):
                self.logger.error(line)
                there_are_fatals = True

        if there_are_fatals:
            return False, lines[:-1]  # Remove last (empty) element

        # Make sure that the job did finish (and was not interrupted externally)

        normal_end = False
        for line in lines:
            if re.match('^INFO: Job completed.*$', line):
                normal_end = True

        if normal_end is False:
            lines[-1] = 'FATAL: ABNORMAL_EXTERNAL_TERMINATION'
            self.logger.error("Calculation interrupted externally")
            return False, lines

        # (Insert any other "non-success" conditions before next section)
        # (e.g.: be very picky about (some) 'WARNING:' messages)
        # Return with success flag

        return True, lines[:-1]  # Remove last (empty) element

    def get_bands(self, bands_path):
        # The parsing is different depending on whether I have Bands or Points.
        # I recognise these two situations by looking at bandskpoints.label
        # (like I did in the plugin)
        tottx = []
        f = open(bands_path)
        tottx = f.read().split()

        #ef = float(tottx[0])
        if self.node.inputs.bandskpoints.labels is None:
            #minfreq, maxfreq = float(tottx[1]), float(tottx[2])
            nbands, nspins, nkpoints = int(tottx[3]), int(tottx[4]), int(tottx[5])
            spinup = np.zeros((nkpoints, nbands))
            spindown = np.zeros((nkpoints, nbands))
            coords = np.zeros(nkpoints)
            if nspins == 2:
                block_length = nbands * 2
            else:
                block_length = nbands
            for i in range(nkpoints):
                # coords[i] = float(tottx[i * (block_length + 3) + 6])
                for j in range(nbands):
                    spinup[i, j] = (float(tottx[i * (block_length + 3) + 6 + j + 3]))
                    if nspins == 2:
                        # Probably wrong! - need to test!!!
                        spindown[i, j] = (float(tottx[i * (nbands * 2 + 3) + 6 + j + 3 + nbands]))
        else:
            #mink, maxk = float(tottx[1]), float(tottx[2])
            #minfreq, maxfreq = float(tottx[3]), float(tottx[4])
            nbands, nspins, nkpoints = int(tottx[5]), int(tottx[6]), int(tottx[7])
            spinup = np.zeros((nkpoints, nbands))
            spindown = np.zeros((nkpoints, nbands))
            coords = np.zeros(nkpoints)
            if nspins == 2:
                block_length = nbands * 2
            else:
                block_length = nbands

            for i in range(nkpoints):
                coords[i] = float(tottx[i * (block_length + 1) + 8])

                for j in range(nbands):
                    spinup[i, j] = (float(tottx[i * (block_length + 1) + 8 + j + 1]))
                    if nspins == 2:
                        # Probably wrong! - need to test!!!
                        spindown[i, j] = (float(tottx[i * (nbands * 2 + 1) + 8 + j + 1 + nbands]))
        if nspins == 2:
            bands = (spinup, spindown)
        elif nspins == 1:
            bands = spinup
        else:
            raise NotImplementedError('nspins=4: non collinear bands not implemented yet')

        return (bands, coords)
