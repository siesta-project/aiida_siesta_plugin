# -*- coding: utf-8 -*-
from __future__ import absolute_import
import numpy as np
from aiida.parsers import Parser
from aiida_siesta.calculations.siesta import SiestaCalculation
from aiida.orm import Dict
from six.moves import range
from aiida.common import OutputParsingError
from aiida.common import exceptions
# See the LICENSE.txt and AUTHORS.txt files.

# TODO Get modules metadata from setup script.

# List of scalar values from CML to be transferred to AiiDA
standard_output_list = [
    'siesta:FreeE', 'siesta:E_KS', 'siesta:Ebs', 'siesta:E_Fermi',
    'siesta:stot'
]  ## leave svec for later

#####################################################
# BEGINNING OF SET OF AUXILIARY FUNCTIONS           #
# They should probably be put in a separate module  #
#####################################################

def text_to_array(s, dtype):
    return np.array(s.replace("\n", "").split(), dtype=dtype)


def get_parsed_xml_doc(xml_path):

    from xml.dom import minidom

    try:
        xmldoc = minidom.parse(xml_path)
    except:
        raise OutputParsingError("Faulty Xml File")
        #xmldoc = None

    # We might want to add some extra consistency checks

    return xmldoc


def get_dict_from_xml_doc(xmldoc):

    # Scalar items

    scalar_dict = {}

    # Metadata items
    itemlist = xmldoc.getElementsByTagName('metadata')
    for s in itemlist:
        #
        # Maybe make sure that 'name' does not contain
        # forbidden characters
        #
        name = s.attributes['name'].value
        value = s.attributes['content'].value
        scalar_dict[name] = value

    # Scalar output items
    # From the last "SCF Finalization" module
    # This means that we do not record non-converged values

    itemlist = xmldoc.getElementsByTagName('module')

    scf_final = None
    for m in itemlist:
        if 'title' in list(m.attributes.keys()):
            # Get last scf finalization module
            if m.attributes['title'].value == "SCF Finalization":
                scf_final = m

    if scf_final is not None:

        # wrapped in <property> elements with a <scalar> child
        props = scf_final.getElementsByTagName('property')

        for s in props:
            if 'dictRef' in list(s.attributes.keys()):
                name = s.attributes['dictRef'].value
                if name in standard_output_list:
                    data = s.getElementsByTagName('scalar')[0]
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
    for m in itemlist:
        # Check the type of the first "step" module, which is a "geometry" one
        if 'serial' in list(m.attributes.keys()):
            if 'dictRef' in list(m.attributes.keys()):
                if m.attributes['dictRef'].value == "Single-Point":
                    return False
                else:
                    return True

    # If we reach this point, something is very wrong
    return False


def get_sizes_info(xmldoc):
    """
     Gets the number of orbitals and non-zero interactions
     """
    no_u = None
    nnz = None
    mesh = None

    itemlist = xmldoc.getElementsByTagName('module')
    for m in itemlist:
        # Process the first "step" module, which is a "geometry" one
        if 'serial' in list(m.attributes.keys()):
            # Get properties here
            props_list = m.getElementsByTagName('property')
            for p in props_list:
                if 'dictRef' in list(p.attributes.keys()):
                    if p.attributes['dictRef'].value == "siesta:no_u":
                        scalar = p.getElementsByTagName('scalar')[0]
                        no_u = int(scalar.childNodes[0].data)
                    if p.attributes['dictRef'].value == "siesta:nnz":
                        scalar = p.getElementsByTagName('scalar')[0]
                        nnz = int(scalar.childNodes[0].data)
                    if p.attributes['dictRef'].value == "siesta:ntm":
                        array = p.getElementsByTagName('array')[0]
                        mesh = [
                            int(s) for s in array.childNodes[0].data.split()
                        ]

    return no_u, nnz, mesh


def get_last_structure(xmldoc, input_structure):

    from aiida.plugins import DataFactory

    itemlist = xmldoc.getElementsByTagName('module')
    #
    # Use the last "geometry" module, and not the
    # "Finalization" one.

    finalmodule = None
    for m in itemlist:
        # Get a "geometry" module by the criteria:
        if 'serial' in list(m.attributes.keys()):
            if 'dictRef' in list(m.attributes.keys()):
                if m.attributes['dictRef'].value != "SCF":
                    finalmodule = m

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

    for a in atoms:
        kind = a.attributes['elementType'].value
        x = a.attributes['x3'].value
        y = a.attributes['y3'].value
        z = a.attributes['z3'].value
        atomlist.append([kind, [float(x), float(y), float(z)]])

    cell = []
    for l in cellvectors:
        data = l.childNodes[0].data.split()
        cell.append([float(s) for s in data])

    # Generally it is better to pass the input structure
    # and reset the data, since site 'names' are not handled by
    # the CML file (at least not in Siesta versions <= 4.0)
    #

    # s = input_structure.copy()
    s = input_structure.clone()
    s.reset_cell(cell)
    new_pos = [atom[1] for atom in atomlist]
    s.reset_sites_positions(new_pos)

    return s


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
    for m in itemlist:
        if 'title' in list(m.attributes.keys()):
            # Get last scf finalization module
            if m.attributes['title'].value == "SCF Finalization":
                scf_final = m

    forces = None
    stress = None

    if scf_final is not None:
        props = scf_final.getElementsByTagName('property')
        for p in props:
            if 'dictRef' in list(p.attributes.keys()):

                if p.attributes['dictRef'].value == 'siesta:forces':
                    mat = p.getElementsByTagName('matrix')[0]
                    # Get flat list and reshape as list of lists
                    # using info on rows and columns in CML file
                    rows = int(mat.attributes['rows'].value)
                    cols = int(mat.attributes['columns'].value)
                    f = mat.childNodes[0].data.split()
                    f = [float(x) for x in f]
                    forces = [f[rows * i:rows * (i + 1)] for i in range(cols)]

                if p.attributes['dictRef'].value == 'siesta:stress':
                    mat = p.getElementsByTagName('matrix')[0]
                    # Get flat list and reshape as list of lists
                    # using info on rows and columns in CML file
                    rows = int(mat.attributes['rows'].value)
                    cols = int(mat.attributes['columns'].value)
                    s = mat.childNodes[0].data.split()
                    s = [float(x) for x in s]
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
        Receives in input a dictionary of retrieved nodes.
        Does all the logic here.
        """
        from aiida.engine import ExitCode

        try:
            output_folder = self.retrieved
        except exceptions.NotExistent:
            return self.exit_codes.ERROR_NO_RETRIEVED_FOLDER

        output_path, messages_path, xml_path, json_path, bands_path = \
            self._fetch_output_files(output_folder)

        out_results = self._get_output_nodes(output_path, messages_path,
                                             xml_path, json_path, bands_path)

        for line in out_results["warnings"]:
            if u'GEOM_NOT_CONV' in line:
                return (self.exit_codes.GEOM_NOT_CONV)
            if u'SCF_NOT_CONV' in line:
                return (self.exit_codes.SCF_NOT_CONV)

        return ExitCode(0)

    def _fetch_output_files(self, out_folder):
        """
        Checks the output folder for standard output and standard error
        files, returns their absolute paths on success.

        :param retrieved: A dictionary of retrieved nodes, as obtained from the
          parser.
        """
        # from aiida.common.datastructures import calc_states
        from aiida.common.exceptions import InvalidOperation
        import os

        list_of_files = out_folder._repository.list_object_names()

        output_path = None
        messages_path = None
        xml_path = None
        json_path = None
        bands_path = None

        if self.node.get_option('output_filename') in list_of_files:
            output_path = os.path.join(
                out_folder._repository._get_base_folder().abspath,
                self.node.get_option('output_filename'))
        else:
            raise OutputParsingError("Output file not retrieved")

        if self.node.process_class._DEFAULT_XML_FILE in list_of_files:
            xml_path = os.path.join(
                out_folder._repository._get_base_folder().abspath,
                self.node.process_class._DEFAULT_XML_FILE)
        else:
            raise OutputParsingError("Xml file not retrieved")

        if self.node.process_class._DEFAULT_JSON_FILE in list_of_files:
            json_path = os.path.join(
                out_folder._repository._get_base_folder().abspath,
                self.node.process_class._DEFAULT_JSON_FILE)
#        else:
#            raise OutputParsingError("json file not retrieved")

        if self.node.process_class._DEFAULT_MESSAGES_FILE in list_of_files:
            messages_path = os.path.join(
                out_folder._repository._get_base_folder().abspath,
                self.node.process_class._DEFAULT_MESSAGES_FILE)
#        else:
#            raise OutputParsingError("message file not retrieved")

        if self.node.process_class._DEFAULT_BANDS_FILE in list_of_files:
            bands_path = os.path.join(
                out_folder._repository._get_base_folder().abspath,
                self.node.process_class._DEFAULT_BANDS_FILE)

        if bands_path is None:
            supposed_to_have_bandsfile = True
            try:
                self.node.inputs.bandskpoints
            except:
                supposed_to_have_bandsfile = False
            if supposed_to_have_bandsfile:
                raise OutputParsingError("bands file not retrieved")

        return output_path, messages_path, xml_path, json_path, bands_path



    def _get_output_nodes(self, output_path, messages_path, xml_path,
                          json_path, bands_path):
        """
        Extracts output nodes from the standard output and standard error
        files. (And XML and JSON files)
        """
        from aiida.orm import TrajectoryData
        import re

        parser_version = '1.0.1'
        parser_info = {}
        parser_info['parser_info'] = 'AiiDA Siesta Parser V. {}'.format(
            parser_version)
        parser_info['parser_warnings'] = []

        #No need for checks anymore

        xmldoc = get_parsed_xml_doc(xml_path)

        in_struc = self.node.inputs.structure

        try:
            in_settings = self.node.inputs.settings
        except exceptions.NotExistent:
            in_settings = None

        result_dict = get_dict_from_xml_doc(xmldoc)

        # Add timing information
        #if json_path is None:
        ###TODO###
        #Not sure how to implement what once was logger.info
        #Also not sure I understood the purpose of this file
        if json_path is not None:
            from .json_time import get_timing_info
            global_time, timing_decomp = get_timing_info(json_path)
            if global_time is None:
                raise OutputParsingError(
                    "Cannot fully parse the time.json file")
            else:
                result_dict["global_time"] = global_time
                result_dict["timing_decomposition"] = timing_decomp

        # Add warnings
        if messages_path is None:
            # Perhaps using an old version of Siesta
            warnings_list = ['WARNING: No MESSAGES file...']
        else:
            successful, warnings_list = self.get_warnings_from_file(
                messages_path)

        ###TODO### or maybe not
        #How do we process this successfull booelan?

        result_dict["warnings"] = warnings_list

        # Add parser info dictionary
        parsed_dict = dict(
            list(result_dict.items()) + list(parser_info.items()))

        output_data = Dict(dict=parsed_dict)

        self.out('output_parameters', output_data)

        # If the structure has changed, save it
        if is_variable_geometry(xmldoc):
            # Get the input structure to copy its site names,
            # as the CML file traditionally contained only the
            # atomic symbols.
            #
            struc = get_last_structure(xmldoc, in_struc)
            # result_list.append((self.get_linkname_outstructure(),struc))
            self.out(self.get_linkname_outstructure(), struc)

        # Save forces and stress in an ArrayData object
        forces, stress = get_final_forces_and_stress(xmldoc)

        if forces is not None and stress is not None:
            # from aiida.orm.nodes.array import ArrayData
            from aiida.orm import ArrayData
            arraydata = ArrayData()
            arraydata.set_array('forces', np.array(forces))
            arraydata.set_array('stress', np.array(stress))
            # result_list.append((self.get_linkname_fs_and_stress(),arraydata))
            self.out(self.get_linkname_fs_and_stress(), arraydata)

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
            self.out(self.get_linkname_bands(), arraybands)
            bandsparameters = Dict(dict={"kp_coordinates": coords})
            self.out(self.get_linkname_bandsparameters(), bandsparameters)

        return result_dict


    def get_warnings_from_file(self, messages_path):
        """
     Generates a list of warnings from the 'MESSAGES' file, which
     contains a line per message, prefixed with 'INFO',
     'WARNING' or 'FATAL'.

     :param messages_path:

     Returns a boolean indicating success (True) or failure (False)
     and a list of strings.
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

        # Make sure that the job did finish (and was not interrupted
        # externally)

        normal_end = False
        for line in lines:
            if re.match('^INFO: Job completed.*$', line):
                normal_end = True

        if normal_end == False:
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
        from aiida.common.exceptions import InputValidationError
        from aiida.common.exceptions import ValidationError
        tottx = []
        f = open(bands_path)
        tottx = f.read().split()

        ef = float(tottx[0])
        if self.node.inputs.bandskpoints.labels is None:
            minfreq, maxfreq = float(tottx[1]), float(tottx[2])
            nbands, nspins, nkpoints = int(tottx[3]), int(tottx[4]), int(
                tottx[5])
            spinup = np.zeros((nkpoints, nbands))
            spindown = np.zeros((nkpoints, nbands))
            coords = np.zeros(nkpoints)
            if (nspins == 2):
                block_length = nbands * 2
            else:
                block_length = nbands
            for i in range(nkpoints):
                # coords[i] = float(tottx[i * (block_length + 3) + 6])
                for j in range(nbands):
                    spinup[i, j] = (float(tottx[i * (block_length + 3) + 6 +
                                                j + 3]))
                    if (nspins == 2):
                        # Probably wrong! - need to test!!!
                        spindown[i, j] = (float(tottx[i * (nbands * 2 + 3) +
                                                      6 + j + 3 + nbands]))
        else:
            mink, maxk = float(tottx[1]), float(tottx[2])
            minfreq, maxfreq = float(tottx[3]), float(tottx[4])
            nbands, nspins, nkpoints = int(tottx[5]), int(tottx[6]), int(
                tottx[7])
            spinup = np.zeros((nkpoints, nbands))
            spindown = np.zeros((nkpoints, nbands))
            coords = np.zeros(nkpoints)
            if (nspins == 2):
                block_length = nbands * 2
            else:
                block_length = nbands

            for i in range(nkpoints):
                coords[i] = float(tottx[i * (block_length + 1) + 8])

                for j in range(nbands):
                    spinup[i, j] = (float(tottx[i * (block_length + 1) + 8 +
                                                j + 1]))
                    if (nspins == 2):
                        # Probably wrong! - need to test!!!
                        spindown[i, j] = (float(tottx[i * (nbands * 2 + 1) +
                                                      8 + j + 1 + nbands]))
        if (nspins == 2):
            bands = (spinup, spindown)
        elif (nspins == 1):
            bands = spinup
        else:
            raise NotImplementedError(
                'nspins=4: non collinear bands not implemented yet')

        return (bands, coords)

    def get_linkname_outstructure(self):
        """
        Returns the name of the link to the output_structure
        Node exists if positions or cell changed.
        """
        return 'output_structure'

    def get_linkname_fs_and_stress(self):
        """
        Returns the name of the link to the output_array
        In Siesta, Node exists to hold the final forces and stress,
        pending the implementation of trajectory data.
        """
        return 'forces_and_stress'

    def get_linkname_bands(self):
        """
        Returns the name of the link to the bands object
        In Siesta, Node exists to hold the bands,
        pending the implementation of trajectory data.
        """
        return 'bands'

    def get_linkname_bandsparameters(self):
        """
        Returns the name of the link to the bands_path.
        X-axis data for bands. Maybe should use ArrayData (db-integrity?).
        """
        return 'bands_parameters'
