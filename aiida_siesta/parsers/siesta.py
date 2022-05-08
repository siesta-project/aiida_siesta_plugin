import os
import numpy as np
from aiida.parsers import Parser
from aiida.orm import Dict
from aiida.common import OutputParsingError
from aiida.common import exceptions

# See the LICENSE.txt and AUTHORS.txt files.

#####################################################
# BEGINNING OF SET OF AUXILIARY FUNCTIONS           #
# They should probably be put in a separate module  #
#####################################################


def get_eps2(eps2_path):
    """
    Reads the eps2_path files to extract an array energy vs eps2
    """

    eps2_list = []

    with open(eps2_path, 'r') as file_h:
        for line in file_h:
            # check if the current line starts with "#"
            if line.startswith("#"):
                pass
            else:
                e_and_eps2 = [float(i) for i in line.split()]
                eps2_list.append(e_and_eps2)

    return eps2_list


def is_polarization_problem(output_path):
    """
    Check the presence of polarization errors.
    """

    thefile = open(output_path)
    lines = thefile.read().split('\n')

    for line in lines:
        if "POLARIZATION: Iteration to find the polarization" in line:
            return True

    thefile.close()

    return False


def get_min_split(output_path):
    """
    Check the presence of split_norm errors in the .out.
    If present, extract the minimum split_norm parameter. If not, return None.
    """

    thefile = open(output_path)
    lines = thefile.read().split('\n')

    min_split_norm = None
    split_norm_error = False

    for line in lines:
        if "split_norm" in line:
            split_norm_error = True
            words = line.split()

    if split_norm_error is not None:
        min_sp = words[4]
        min_split_norm = float(min_sp[:-1])

    return min_split_norm


def get_parsed_xml_doc(xml_path):

    from xml.dom import minidom

    try:
        xmldoc = minidom.parse(xml_path)
    except EOFError:
        raise OutputParsingError("Faulty Xml File")

    return xmldoc


def get_dict_from_xml_doc(xmldoc):

    # List of scalar values from CML to be transferred to AiiDA
    #pylint: disable=invalid-name
    standard_output_list = ['siesta:FreeE', 'siesta:E_KS', 'siesta:Ebs', 'siesta:E_Fermi', 'siesta:stot']

    # Scalar items
    scalar_dict = {}
    # Metadata items
    itemlist = xmldoc.getElementsByTagName('metadata')
    # Transform metagata in scalar
    for item in itemlist:
        # Maybe make sure that 'name' does not contain forbidden characters
        name = item.attributes['name'].value
        value = item.attributes['content'].value
        scalar_dict[name] = value

    # Look for "SCF Finalization" module, present if at least one scf converged
    itemlist = xmldoc.getElementsByTagName('module')
    scf_final = None
    for item in itemlist:
        if 'title' in list(item.attributes.keys()):
            # Get last scf finalization module
            if item.attributes['title'].value == "SCF Finalization":
                scf_final = item
    # In a geom_optimization run, we catch the data of last run even if the
    # geom_optimization failed

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

    #Detect if it was a geometry optimization (relax) or a single point calculation
    scalar_dict['variable_geometry'] = is_variable_geometry(xmldoc)

    # Sizes of orbital set (and non-zero interactions), and mesh
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
            ref = item.attributes['dictRef'].value
            # This is a very simple-minded approach, since
            # there might be Lua runs which are not properly
            # geometry optimizations.
            if ref in ("Geom. Optim", "LUA"):
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

    # Use the last "geometry" module, and not the "Finalization" one.
    finalmodule = None
    for item in itemlist:
        # Get a "geometry" module by the criteria:
        if 'serial' in list(item.attributes.keys()):
            if 'dictRef' in list(item.attributes.keys()):
                if item.attributes['dictRef'].value != "SCF":
                    finalmodule = item

    # In case there is no appropriate data, fall back and at least return the initial structure
    # (this should not be necessary, as the initial Geometry module is opened very soon)
    if finalmodule is None:
        return False, input_structure

    # When using floating sites, Siesta associates 'atomic positions' to them, and
    # the structure (and forces) in the XML file include these fake atoms.
    # In order to return physical structures and forces, we need to remove them.
    # Recall that the input structure is the physical one, as the floating sites
    # are specified in the 'basis' input

    number_of_real_atoms = len(input_structure.sites)

    atoms = finalmodule.getElementsByTagName('atom')
    cellvectors = finalmodule.getElementsByTagName('latticeVector')

    atomlist = []

    for atm in atoms[0:number_of_real_atoms]:
        kind = atm.attributes['elementType'].value
        x = atm.attributes['x3'].value
        y = atm.attributes['y3'].value
        z = atm.attributes['z3'].value
        atomlist.append([kind, [float(x), float(y), float(z)]])

    cell = []
    for latt in cellvectors:
        data = latt.childNodes[0].data.split()
        cell.append([float(s) for s in data])

    # Generally it is better to clone the input structure and reset the data, since site
    # 'names' are not handled by the CML file (at least not in Siesta versions <= 4.0)

    import copy
    from aiida.orm.nodes.data.structure import Site
    new_structure = input_structure.clone()
    new_structure.reset_cell(cell)
    new_structure.clear_sites()
    for i in range(number_of_real_atoms):
        new_site = Site(site=input_structure.sites[i])
        new_site.position = copy.deepcopy(atomlist[i][1])
        new_structure.append_site(new_site)

    # The most obvious alternative does not work, as the reset method below does not
    # work if the numbers do not match.
    #
    ## new_pos = [atom[1] for atom in atomlist]
    ## new_structure.reset_sites_positions(new_pos)

    return True, new_structure


def get_final_forces_and_stress(xmldoc):
    #
    # Extracts final forces and stress as lists of lists...
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


class SiestaParser(Parser):
    """
    Parser for the output of Siesta.
    """

    _version = '1.4.0'

    def parse(self, **kwargs):  # noqa: MC0001  - is mccabe too complex funct -
        """
        Receives in input a dictionary of retrieved nodes. Does all the logic here.
        """
        from aiida.engine import ExitCode

        parser_info = {}
        parser_info['parser_info'] = 'AiiDA Siesta Parser V. {}'.format(self._version)

        try:
            output_folder = self.retrieved
        except exceptions.NotExistent:
            raise OutputParsingError("Folder not retrieved")

        output_path, messages_path, xml_path, json_path, bands_path, basis_enthalpy_path, harris_en_path, eps2_path = \
            self._fetch_output_files(output_folder)

        if xml_path is None:
            raise OutputParsingError("Xml file not retrieved")
        xmldoc = get_parsed_xml_doc(xml_path)
        result_dict = get_dict_from_xml_doc(xmldoc)

        if output_path is None:
            raise OutputParsingError("output file not retrieved")

        output_dict = dict(list(result_dict.items()) + list(parser_info.items()))

        warnings_list = []

        if json_path is not None:
            from .json_time import get_timing_info
            global_time, timing_decomp = get_timing_info(json_path)
            if global_time is None:
                warnings_list.append(["Cannot fully parse the time.json file"])
            else:
                output_dict["global_time"] = global_time
                output_dict["timing_decomposition"] = timing_decomp

        if basis_enthalpy_path is not None:
            the_file = open(basis_enthalpy_path)
            bas_enthalpy = float(the_file.read().split()[0])
            the_file.close()
            output_dict["basis_enthalpy"] = bas_enthalpy
            output_dict["basis_enthalpy_units"] = "eV"
        else:
            warnings_list.append(["BASIS_ENTHALPY file not retrieved"])

        if harris_en_path is not None:
            the_file = open(harris_en_path)
            harr_enthalpy = float(the_file.read().split()[0])
            the_file.close()
            output_dict["harris_basis_enthalpy"] = harr_enthalpy
            output_dict["harris_basis_enthalpy_units"] = "eV"
        else:
            warnings_list.append(["HARRIS_BASIS_ENTHALPY file not retrieved"])

        have_errors_to_analyse = False
        if messages_path is None:
            # Perhaps using an old version of Siesta
            warnings_list.append(['WARNING: No MESSAGES file, could not check if calculation terminated correctly'])
        else:
            have_errors_to_analyse = True
            #succesful when "INFO: Job completed" is present in message files
            succesful, from_message = self._get_warnings_from_file(messages_path)
            warnings_list.append(from_message)
        output_dict["warnings"] = warnings_list

        # An output_parametrs port is always return, even if only parser's info are present
        output_data = Dict(dict=output_dict)
        self.out('output_parameters', output_data)

        #
        # When using floating sites, Siesta associates 'atomic positions' to them, and
        # the structure and forces in the XML file include these fake atoms.
        # In order to return physical structures and forces, we need to remove them.
        # Recall that the input structure is the physical one, and the floating sites
        # are specified in the 'basis' input
        #
        physical_structure = self.node.inputs.structure
        number_of_real_atoms = len(physical_structure.sites)

        # If the structure has changed, save it
        #
        if output_dict['variable_geometry']:
            in_struc = self.node.inputs.structure
            # The next function never fails. If problems arise, the initial structure is
            # returned. The input structure is also necessary because the CML file
            # traditionally contains only the atomic symbols and not the site names.
            # The returned structure does not have any floating atoms, they are removed
            # in the `get_last_structure` call.
            success, out_struc = get_last_structure(xmldoc, in_struc)
            if not success:
                self.logger.warning("Problem in parsing final structure, returning inp structure in output_structure")

            self.out('output_structure', out_struc)

        # Attempt to parse forces and stresses. In case of failure "None" is returned.
        # Therefore the function never crashes
        forces, stress = get_final_forces_and_stress(xmldoc)
        if forces is not None and stress is not None:
            from aiida.orm import ArrayData
            arraydata = ArrayData()
            arraydata.set_array('forces', np.array(forces[0:number_of_real_atoms]))
            arraydata.set_array('stress', np.array(stress))
            self.out('forces_and_stress', arraydata)

        #Attempt to parse the ion files. Files ".ion.xml" are not produced by siesta if ions file are used
        #in input (`user-basis = T`). This explains the first "if" statement. The SiestaCal input is called
        #`ions__El` (El is the element label) therefore we look for the str "ions" in any of the inputs name.
        if not any(["ions" in inp for inp in self.node.inputs]):  #pylint: disable=too-many-nested-blocks
            from aiida_siesta.data.ion import IonData
            ions = {}
            #Ions from the structure
            in_struc = self.node.inputs.structure
            for kind in in_struc.get_kind_names():
                ion_file_name = kind + ".ion.xml"
                if ion_file_name in output_folder._repository.list_object_names():
                    ion_file_path = os.path.join(output_folder._repository._get_base_folder().abspath, ion_file_name)
                    ions[kind] = IonData(ion_file_path)
                else:
                    self.logger.warning(f"no ion file retrieved for {kind}")
            #Ions from floating_sites
            if "basis" in self.node.inputs:
                basis_dict = self.node.inputs.basis.get_dict()
                if "floating_sites" in basis_dict:
                    floating_kinds = []
                    for orb in basis_dict["floating_sites"]:
                        if orb["name"] not in floating_kinds:
                            floating_kinds.append(orb["name"])
                            ion_file_name = orb["name"] + ".ion.xml"
                            if ion_file_name in output_folder._repository.list_object_names():
                                ion_file_path = os.path.join(
                                    output_folder._repository._get_base_folder().abspath, ion_file_name
                                )
                                ions[orb["name"]] = IonData(ion_file_path)
                            else:
                                self.logger.warning(f"no ion file retrieved for {orb['name']}")
            #Return the outputs
            if ions:
                self.out('ion_files', ions)

        # Error analysis
        if have_errors_to_analyse:
            # No metter if "INFO: Job completed" is present (succesfull) or not, we check for known
            # errors. They might apprear as WARNING (therefore with succesful True) or FATAL
            # (succesful False)
            for line in from_message:
                if u'split options' in line:
                    min_split = get_min_split(output_path)
                    if min_split:
                        self.logger.error("Error in split_norm option. Minimum value is {}".format(min_split))
                        return self.exit_codes.SPLIT_NORM
                if u'sys::die' in line:
                    #This is the situation when siesta dies with no specified error
                    #to be reported in "MESSAGES", unfortunately some interesting cases
                    #are treated in this way, we explore the .out file for more insights.
                    if is_polarization_problem(output_path):
                        return self.exit_codes.BASIS_POLARIZ
                if u'SCF_NOT_CONV' in line:
                    return self.exit_codes.SCF_NOT_CONV
                if u'GEOM_NOT_CONV' in line:
                    return self.exit_codes.GEOM_NOT_CONV

        #Because no known error has been found, attempt to parse bands if requested
        if bands_path is None:
            if "bandskpoints" in self.node.inputs:
                return self.exit_codes.BANDS_FILE_NOT_PRODUCED
        else:
            #bands, coords = self._get_bands(bands_path)
            try:
                bands = self._get_bands(bands_path)
            except (ValueError, IndexError):
                return self.exit_codes.BANDS_PARSE_FAIL
            from aiida.orm import BandsData
            arraybands = BandsData()
            #Reset the cell for KpointsData of bands, necessary
            #for bandskpoints without cell and if structure changed
            bkp = self.node.inputs.bandskpoints.clone()
            if output_dict['variable_geometry']:
                bkp.set_cell_from_structure(out_struc)
            else:
                bkp.set_cell_from_structure(self.node.inputs.structure)
            arraybands.set_kpointsdata(bkp)
            arraybands.set_bands(bands, units="eV")
            self.out('bands', arraybands)
            #bandsparameters = Dict(dict={"kp_coordinates": coords})
            #self.out('bands_parameters', bandsparameters)

        #Because no known error has been found, attempt to parse EPSIMG file if requested
        if eps2_path is None:
            if "optical" in self.node.inputs:
                return self.exit_codes.EPS2_FILE_NOT_PRODUCED
        else:
            eps2_list = get_eps2(eps2_path)
            optical_eps2 = ArrayData()
            optical_eps2.set_array('e_eps2', np.array(eps2_list))
            self.out('optical_eps2', optical_eps2)

        #At the very end, return a particular exit code if "INFO: Job completed"
        #was not present in the MESSAGES file, but no known error is detected.
        if have_errors_to_analyse:
            if not succesful:
                self.logger.error(
                    'The calculation finished without "INFO: Job completed", but no '
                    'error could be processed. Might be that the calculation was killed externally'
                )
                return self.exit_codes.UNEXPECTED_TERMINATION

        return ExitCode(0)

    def _fetch_output_files(self, out_folder):
        """
        Checks the output folder for standard output and standard error files, returns their absolute paths
        or "None" in case the file is not found in the remote folder.
        """

        list_of_files = out_folder._repository.list_object_names()

        output_path = None
        messages_path = None
        xml_path = None
        json_path = None
        bands_path = None
        basis_enthalpy_path = None
        harris_en_path = None
        eps2_path = None

        if self.node.get_option('output_filename') in list_of_files:
            oufil = self.node.get_option('output_filename')
            output_path = os.path.join(out_folder._repository._get_base_folder().abspath, oufil)

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

        if self.node.process_class._BASIS_ENTHALPY_FILE in list_of_files:
            basis_enthalpy_path = os.path.join(
                out_folder._repository._get_base_folder().abspath, self.node.process_class._BASIS_ENTHALPY_FILE
            )

        if self.node.process_class._HARRIS_ENTHALPY_FILE in list_of_files:
            harris_en_path = os.path.join(
                out_folder._repository._get_base_folder().abspath, self.node.process_class._HARRIS_ENTHALPY_FILE
            )

        namebandsfile = str(self.node.get_option('prefix')) + ".bands"
        if namebandsfile in list_of_files:
            bands_path = os.path.join(out_folder._repository._get_base_folder().abspath, namebandsfile)

        nameeps2file = str(self.node.get_option('prefix')) + ".EPSIMG"
        if nameeps2file in list_of_files:
            eps2_path = os.path.join(out_folder._repository._get_base_folder().abspath, nameeps2file)

        return output_path, messages_path, xml_path, json_path, bands_path, basis_enthalpy_path, harris_en_path, eps2_path  # noqa

    def _get_warnings_from_file(self, messages_path):
        """
        Generates a list of warnings from the 'MESSAGES' file, which  contains a line per message,
        prefixed with 'INFO', 'WARNING' or 'FATAL'.

        Returns a boolean indicating success (True) or failure (False) and a list of strings.
        """
        thefile = open(messages_path)
        lines = thefile.read().split('\n')  # There will be a final '' element

        import re

        # Search for 'FATAL:' and 'WARNING:' messages, log them
        normal_end = False
        for line in lines:
            if re.match('^FATAL:.*$', line):
                self.logger.error(line)
            if re.match('^WARNING:.*$', line):
                self.logger.warning(line)
            if re.match('^INFO: Job completed.*$', line):
                normal_end = True

        #If there is "INFO: Job completed" we return True and possible warnings
        #in a list (lines)
        if normal_end:
            return True, lines[:-1]

        #Otherwise we return False and possible FATAL, WARNING in a list. If
        #no FATAL, WARNING are present, this means that MESSAGES file is eampty.
        #Therefore the list is eampty but the Bool knows that something was wrong.
        return False, lines[:-1]

    def _get_bands(self, bands_path):
        # The parsing is different depending on whether I have Bands or Points.
        # I recognise these two situations by looking at bandskpoints.label
        # (like I did in the plugin)
        tottx = []
        thefile = open(bands_path)
        tottx = thefile.read().split()
        #ef = float(tottx[0])
        if self.node.inputs.bandskpoints.labels is None:
            #minfreq, maxfreq = float(tottx[1]), float(tottx[2])
            nbands, nspins, nkpoints = int(tottx[3]), int(tottx[4]), int(tottx[5])
            spinup = np.zeros((nkpoints, nbands))
            spindown = np.zeros((nkpoints, nbands))
            #coords = np.zeros(nkpoints)
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
            #coords = np.zeros(nkpoints)
            if nspins == 2:
                block_length = nbands * 2
            else:
                block_length = nbands

            for i in range(nkpoints):
                #coords[i] = float(tottx[i * (block_length + 1) + 8])
                for j in range(nbands):
                    spinup[i, j] = (float(tottx[i * (block_length + 1) + 8 + j + 1]))
                    if nspins == 2:
                        spindown[i, j] = (float(tottx[i * (nbands * 2 + 1) + 8 + j + 1 + nbands]))
        if nspins == 2:
            bands = (spinup, spindown)
        elif nspins == 1:
            bands = spinup
        else:
            raise ValueError('detected nspin > 2, something wrong')

        return bands  #, coords
