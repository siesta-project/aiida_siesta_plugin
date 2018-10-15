# -*- coding: utf-8 -*-
import os

from aiida.common.constants import elements
from aiida.common.datastructures import CalcInfo, CodeInfo
from aiida.common.exceptions import InputValidationError
from aiida.common.utils import classproperty
from aiida.orm.calculation.job import JobCalculation
from aiida.orm.data.array.kpoints import KpointsData
from aiida.orm.data.parameter import ParameterData
from aiida.orm.data.remote import RemoteData
from aiida.orm.data.structure import StructureData

from aiida_siesta.data.psf import PsfData, get_pseudos_from_structure
# Module with fdf-aware dictionary
from tkdict import FDFDict

__copyright__ = u"Copyright (c), 2015, ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE (Theory and Simulation of Materials (THEOS) and National Centre for Computational Design and Discovery of Novel Materials (NCCR MARVEL)), Switzerland and ROBERT BOSCH LLC, USA. All rights reserved."
__license__ = "MIT license, see LICENSE.txt file"
__version__ = "0.9.10"
__contributors__ = "Victor M. Garcia-Suarez, Alberto Garcia, Emanuele Bosoni"


class SiestaCalculation(JobCalculation):
    """
    Siesta calculator class for AiiDA.
    """
    _siesta_plugin_version = 'aiida-0.12.0--plugin-0.9.10'

    def _init_internal_params(self):
        super(SiestaCalculation, self)._init_internal_params()

        # Default Siesta output parser provided by AiiDA
        self._default_parser = "siesta.parser"

        # Keywords that cannot be set
        # We need to canonicalize this!

        self._aiida_blocked_keywords = ['system-name', 'system-label']

        self._aiida_blocked_keywords.append('number-of-species')
        self._aiida_blocked_keywords.append('number-of-atoms')
        #
        #
        self._aiida_blocked_keywords.append('latticeconstant')
        self._aiida_blocked_keywords.append('lattice-constant')
        self._aiida_blocked_keywords.append('atomic-coordinates-format')
        self._aiida_blocked_keywords.append('atomiccoordinatesformat')
        self._aiida_blocked_keywords.append('use-tree-timer')
        self._aiida_blocked_keywords.append('xml-write')

        # Default input and output files
        self._DEFAULT_INPUT_FILE = 'aiida.in'
        self._DEFAULT_OUTPUT_FILE = 'aiida.out'
        self._DEFAULT_XML_FILE = 'aiida.xml'
        self._DEFAULT_JSON_FILE = 'time.json'
        self._DEFAULT_MESSAGES_FILE = 'MESSAGES'
        self._DEFAULT_BANDS_FILE = 'aiida.bands'

        self._PSEUDO_SUBFOLDER = './'
        self._OUTPUT_SUBFOLDER = './'
        self._PREFIX = 'aiida'
        self._INPUT_FILE_NAME = 'aiida.fdf'
        self._OUTPUT_FILE_NAME = 'aiida.out'
        self._XML_FILE_NAME = 'aiida.xml'
        self._JSON_FILE_NAME = 'time.json'
        self._MESSAGES_FILE_NAME = 'MESSAGES'
        self._BANDS_FILE_NAME = 'aiida.bands'

        # in restarts, it will copy from the parent the following
        # (fow now, just the density matrix file)
        self._restart_copy_from = os.path.join(self._OUTPUT_SUBFOLDER, '*.DM')

        # in restarts, it will copy the previous folder in the following one
        self._restart_copy_to = self._OUTPUT_SUBFOLDER

    @classproperty
    def _use_methods(cls):
        """
        Extend the parent _use_methods with further keys.
        """
        retdict = JobCalculation._use_methods

        retdict['kpoints'] = {
            'valid_types': KpointsData,
            'additional_parameter': None,
            'linkname': 'kpoints',
            'docstring': "Use the node defining the kpoint sampling to use",
        }
        retdict["structure"] = {
            'valid_types': StructureData,
            'additional_parameter': None,
            'linkname': 'structure',
            'docstring': "Choose the input structure to use",
        }
        retdict["basis"] = {
            'valid_types': ParameterData,
            'additional_parameter': None,
            'linkname': 'basis',
            'docstring': "Choose the input basis to use",
        }
        retdict["settings"] = {
            'valid_types': ParameterData,
            'additional_parameter': None,
            'linkname': 'settings',
            'docstring': "Use an additional node for special settings",
        }
        retdict["parameters"] = {
            'valid_types': ParameterData,
            'additional_parameter': None,
            'linkname': 'parameters',
            'docstring': ("Use a node that specifies the input parameters "
                          "for the namelists"),
        }
        retdict["parent_folder"] = {
            'valid_types': RemoteData,
            'additional_parameter': None,
            'linkname': 'parent_calc_folder',
            'docstring': ("Use a remote folder as parent folder (for "
                          "restarts and similar"),
        }
        retdict["pseudo"] = {
            'valid_types': PsfData,
            'additional_parameter': "kind",
            'linkname': cls._get_linkname_pseudo,
            'docstring': ("Use a node for the PSF pseudopotential of one of "
                          "the elements in the structure. You have to pass "
                          "an additional parameter ('kind') specifying the "
                          "name of the structure kind (i.e., the name of "
                          "the species) for which you want to use this "
                          "pseudo. You can pass either a string, or a "
                          "list of strings if more than one kind uses the "
                          "same pseudo"),
        }
        retdict['bandskpoints'] = {
            'valid_types': KpointsData,
            'additional_parameter': None,
            'linkname': 'bandskpoints',
            'docstring': ("Use the node defining the kpoint sampling"
                          "to use for bands calculation"),
        }
        return retdict

    def _prepare_for_submission(self, tempfolder, inputdict):
        """
        This is the routine to be called when you want to create
        the input files and related stuff with a plugin.

        :param tempfolder: a aiida.common.folders.Folder subclass where
                           the plugin should put all its files.
        :param inputdict: a dictionary with the input nodes, as they would
                be returned by get_inputdata_dict (without the Code!)
        """

        local_copy_list = []
        remote_copy_list = []

        # Process the settings dictionary first
        # Settings can be undefined, and defaults to an empty dictionary
        settings = inputdict.pop(self.get_linkname('settings'), None)
        if settings is None:
            settings_dict = {}
        else:
            if not isinstance(settings, ParameterData):
                raise InputValidationError(
                    "settings, if specified, must be of "
                    "type ParameterData")

            # Settings converted to UPPERCASE
            # Presumably to standardize the usage and avoid
            # ambiguities
            settings_dict = _uppercase_dict(
                settings.get_dict(), dict_name='settings')

        try:
            parameters = inputdict.pop(self.get_linkname('parameters'))
        except KeyError:
            raise InputValidationError("No parameters specified for this "
                                       "calculation")
        if not isinstance(parameters, ParameterData):
            raise InputValidationError("parameters is not of type "
                                       "ParameterData")

        # Basis can be undefined, and defaults to an empty dictionary,
        # Siesta will use default parameters
        basis = inputdict.pop(self.get_linkname('basis'), None)
        if basis is None:
            input_basis = {}
        else:
            if not isinstance(basis, ParameterData):
                raise InputValidationError("basis not of type ParameterData")

#            input_basis=FDFDict(basis.get_dict())
            input_basis = basis.get_dict()

        try:
            structure = inputdict.pop(self.get_linkname('structure'))
        except KeyError:
            raise InputValidationError("No structure specified for this "
                                       "calculation")
        if not isinstance(structure, StructureData):
            raise InputValidationError(
                "structure is not of type StructureData")

        # k-points
        # It is now possible to elide the kpoints node.
        #
        # Note also that a *different* set of k-points is needed if a band
        # calculation is carried out. This should be specified somehow in
        # the 'settings' dictionary (see QE example...)

        kpoints = inputdict.pop(self.get_linkname('kpoints'), None)
        if kpoints is None:
            # Do nothing. Assume it is a gamma-point calculation
            pass
        else:
            if not isinstance(kpoints, KpointsData):
                raise InputValidationError("kpoints, if specified, must be of "
                                           "type KpointsData")

        bandskpoints = inputdict.pop(self.get_linkname('bandskpoints'), None)
        if bandskpoints is None:
            flagbands = False
        else:
            flagbands = True
            if not isinstance(bandskpoints, KpointsData):
                raise InputValidationError(
                    "kpoints for bands is not of type KpointsData")

        pseudos = {}
        # I create here a dictionary that associates each kind name to a pseudo
        for link in inputdict.keys():
            if link.startswith(self._get_linkname_pseudo_prefix()):
                kindstring = link[len(self._get_linkname_pseudo_prefix()):]
                kinds = kindstring.split('_')
                the_pseudo = inputdict.pop(link)
                if not isinstance(the_pseudo, PsfData):
                    raise InputValidationError(
                        "Pseudo for kind(s) {} is not "
                        " of type PsfData".format(",".join(kinds)))
                #
                # Note that we can associate the same pseudo object to different
                # atom kinds
                #
                for kind in kinds:
                    pseudos[kind] = the_pseudo

        parent_calc_folder = inputdict.pop(
            self.get_linkname('parent_folder'), None)
        if parent_calc_folder is not None:
            if not isinstance(parent_calc_folder, RemoteData):
                raise InputValidationError("parent_calc_folder, if specified,"
                                           "must be of type RemoteData")

        try:
            code = inputdict.pop(self.get_linkname('code'))
        except KeyError:
            raise InputValidationError(
                "No code specified for this calculation")

        # Here, there should be no more parameters...
        if inputdict:
            raise InputValidationError(
                "The following input data nodes are "
                "unrecognized: {}".format(inputdict.keys()))

        # Check structure, get species, check peudos
        kindnames = [k.name for k in structure.kinds]
        if set(kindnames) != set(pseudos.keys()):
            err_msg = ("Mismatch between the defined pseudos and the list of "
                       "kinds of the structure. Pseudos: {}; kinds: {}".format(
                           ",".join(pseudos.keys()),
                           ",".join(list(kindnames))))
            raise InputValidationError(err_msg)

        ##############################
        # END OF INITIAL INPUT CHECK #
        ##############################

        #
        # There should be a warning for duplicated (canonicalized) keys
        # in the original dictionary in the script

        input_params = FDFDict(parameters.get_dict())

        # Look for blocked keywords and
        # add the proper values to the dictionary

        for blocked_key in self._aiida_blocked_keywords:
            canonical_blocked = FDFDict.translate_key(blocked_key)
            for key in input_params:
                if key == canonical_blocked:
                    raise InputValidationError(
                        "You cannot specify explicitly the '{}' flag in the "
                        "input parameters".format(
                            input_params.get_last_key(key)))

        input_params.update({'system-name': self._PREFIX})
        input_params.update({'system-label': self._PREFIX})
        input_params.update({'use-tree-timer': 'T'})
        input_params.update({'xml-write': 'T'})

        input_params.update({'number-of-species': len(structure.kinds)})
        input_params.update({'number-of-atoms': len(structure.sites)})
        #
        # Regarding the lattice-constant parameter:
        # -- The variable "alat" is not typically kept anywhere, and
        # has already been used to define the vectors.
        # We need to specify that the units of these vectors are Ang...

        input_params.update({'lattice-constant': '1.0 Ang'})

        # Note that this  will break havoc with the band-k-points "pi/a"
        # option. The use of this option should be banned.

        # Note that the implicit coordinate convention of the Structure
        # class corresponds to the "Ang" convention in Siesta.
        # That is why the "atomic-coordinates-format" keyword is blocked
        # and reset.
        input_params.update({'atomic-coordinates-format': 'Ang'})

        # ============== Preparation of input data ===============
        #

        # ------------ CELL_PARAMETERS -----------
        cell_parameters_card = "%block lattice-vectors\n"
        for vector in structure.cell:
            cell_parameters_card += ("{0:18.10f} {1:18.10f} {2:18.10f}"
                                     "\n".format(*vector))
        cell_parameters_card += "%endblock lattice-vectors\n"

        # ------------- ATOMIC_SPECIES ------------
        # I create the subfolder that will contain the pseudopotentials
        tempfolder.get_subfolder(self._PSEUDO_SUBFOLDER, create=True)
        # I create the subfolder with the output data
        tempfolder.get_subfolder(self._OUTPUT_SUBFOLDER, create=True)

        atomic_species_card_list = []

        # Dictionary to get the atomic number of a given element
        datmn = dict([(v['symbol'], k) for k, v in elements.iteritems()])

        spind = {}
        spcount = 0
        for kind in structure.kinds:

            ps = pseudos[kind.name]

            # I add this pseudo file to the list of files to copy,
            # with the appropiate name
            local_copy_list.append((ps.get_file_abs_path(), os.path.join(
                self._PSEUDO_SUBFOLDER, kind.name + ".psf")))
            spcount += 1
            spind[kind.name] = spcount
            atomic_species_card_list.append("{0:5} {1:5} {2:5}\n".format(
                spind[kind.name], datmn[kind.symbol], kind.name.rjust(6)))

        atomic_species_card_list = (
            ["%block chemicalspecieslabel\n"] + list(atomic_species_card_list))
        atomic_species_card = "".join(atomic_species_card_list)
        atomic_species_card += "%endblock chemicalspecieslabel\n"
        # Free memory
        del atomic_species_card_list

        # ------------ ATOMIC_POSITIONS -----------
        atomic_positions_card_list = [
            "%block atomiccoordinatesandatomicspecies\n"
        ]
        countatm = 0
        for site in structure.sites:
            countatm += 1
            atomic_positions_card_list.append(
                "{0:18.10f} {1:18.10f} {2:18.10f} {3:4} {4:6} {5:6}\n".format(
                    site.position[0], site.position[1], site.position[2],
                    spind[site.kind_name], site.kind_name.rjust(6), countatm))
        atomic_positions_card = "".join(atomic_positions_card_list)
        del atomic_positions_card_list  # Free memory
        atomic_positions_card += "%endblock atomiccoordinatesandatomicspecies\n"

        # --------------- K-POINTS ----------------
        if kpoints is not None:
            #
            # Get a mesh for sampling
            # NOTE that there is not yet support for the 'kgrid-cutoff'
            # option in Siesta.
            #
            try:
                mesh, offset = kpoints.get_kpoints_mesh()
                has_mesh = True
            except AttributeError:
                raise InputValidationError("K-point sampling for scf "
                                           "must be given in mesh form")

            kpoints_card_list = ["%block kgrid_monkhorst_pack\n"]
            #
            # This will fail if has_mesh is False (for the case of a list),
            # since in that case 'offset' is undefined.
            #
            kpoints_card_list.append("{0:6} {1:6} {2:6} {3:18.10f}\n".format(
                mesh[0], 0, 0, offset[0]))
            kpoints_card_list.append("{0:6} {1:6} {2:6} {3:18.10f}\n".format(
                0, mesh[1], 0, offset[1]))
            kpoints_card_list.append("{0:6} {1:6} {2:6} {3:18.10f}\n".format(
                0, 0, mesh[2], offset[2]))

            kpoints_card = "".join(kpoints_card_list)
            kpoints_card += "%endblock kgrid_monkhorst_pack\n"
            del kpoints_card_list

        # --------------- K-POINTS-FOR-BANDS ----------------!
        #This part is computed only if flagbands=True
        #Two possibility are supported in Siesta: BandLines ad BandPoints
        #At the moment the user can't choose directly one of the two options
        #BandsLine is set automatically if bandskpoints has labels,
        #BandsPoints if bandskpoints has no labels
        #BandLinesScale =pi/a is not supported at the moment because currently
        #a=1 always. BandLinesScale ReciprocalLatticeVectors is always set
        if flagbands:
            bandskpoints_card_list = [
                "BandLinesScale ReciprocalLatticeVectors\n"
            ]
            if bandskpoints.labels == None:
                bandskpoints_card_list.append("%block BandPoints\n")
                for s in bandskpoints.get_kpoints():
                    bandskpoints_card_list.append(
                        "{0:8.3f} {1:8.3f} {2:8.3f} \n".format(
                            s[0], s[1], s[2]))
                fbkpoints_card = "".join(bandskpoints_card_list)
                fbkpoints_card += "%endblock BandPoints\n"
            else:
                bandskpoints_card_list.append("%block BandLines\n")
                savs = []
                listforbands = bandskpoints.get_kpoints()
                for s, m in bandskpoints.labels:
                    savs.append(s)
                rawindex = 0
                for s, m in bandskpoints.labels:
                    rawindex = rawindex + 1
                    x, y, z = listforbands[s]
                    if rawindex == 1:
                        bandskpoints_card_list.append(
                            "{0:3} {1:8.3f} {2:8.3f} {3:8.3f} {4:1} \n".format(
                                1, x, y, z, m))
                    else:
                        bandskpoints_card_list.append(
                            "{0:3} {1:8.3f} {2:8.3f} {3:8.3f} {4:1} \n".format(
                                s - savs[rawindex - 2], x, y, z, m))
                fbkpoints_card = "".join(bandskpoints_card_list)
                fbkpoints_card += "%endblock BandLines\n"
            del bandskpoints_card_list

        # ================ Namelists and cards ===================

        input_filename = tempfolder.get_abs_path(self._INPUT_FILE_NAME)

        with open(input_filename, 'w') as infile:
            # here print keys and values tp file

            for k, v in sorted(input_params.iteritems()):
                infile.write(get_input_data_text(k, v))
                # ,mapping=mapping_species))

            # Basis set info is processed just like the general
            # parameters section. Some discipline is needed to
            # put any basis-related parameters (including blocks)
            # in the basis dictionary in the input script.
            #
            if basis is not None:
                infile.write("#\n# -- Basis Set Info follows\n#\n")
                for k, v in input_basis.iteritems():
                    infile.write(get_input_data_text(k, v))

            # Write previously generated cards now
            infile.write("#\n# -- Structural Info follows\n#\n")
            infile.write(atomic_species_card)
            infile.write(cell_parameters_card)
            infile.write(atomic_positions_card)
            if kpoints is not None:
                infile.write("#\n# -- K-points Info follows\n#\n")
                infile.write(kpoints_card)
            if flagbands:
                infile.write("#\n# -- Bandlines/Bandpoints Info follows\n#\n")
                infile.write(fbkpoints_card)

            # Write max wall-clock time
            infile.write("#\n# -- Max wall-clock time block\n#\n")
            infile.write(
                "max.walltime {}".format(self.get_max_wallclock_seconds()))

        # ------------------------------------- END of fdf file creation

        # operations for restart

        # The presence of a 'parent_calc_folder' input node signals
        # that we want to get something from there, as indicated in the
        # self._restart_copy_from attribute.
        # In Siesta's case, for now, it is just the density-matrix file
        #
        # It will be copied to the current calculation's working folder.

        if parent_calc_folder is not None:
            remote_copy_list.append(
                (parent_calc_folder.get_computer().uuid, os.path.join(
                    parent_calc_folder.get_remote_path(),
                    self._restart_copy_from), self._restart_copy_to))

        calcinfo = CalcInfo()

        calcinfo.uuid = self.uuid
        #
        # Empty command line by default
        # Why use 'pop' ?
        cmdline_params = settings_dict.pop('CMDLINE', [])

        # Comment this paragraph better, if applicable to Siesta
        #
        #we commented calcinfo.stin_name and added it here in cmdline_params
        #in this way the mpirun ... pw.x ... < aiida.in
        #is replaced by mpirun ... pw.x ... -in aiida.in
        # in the scheduler, _get_run_line, if cmdline_params is empty, it
        # simply uses < calcinfo.stin_name

        if cmdline_params:
            calcinfo.cmdline_params = list(cmdline_params)
        calcinfo.local_copy_list = local_copy_list
        calcinfo.remote_copy_list = remote_copy_list

        calcinfo.stdin_name = self._INPUT_FILE_NAME
        calcinfo.stdout_name = self._OUTPUT_FILE_NAME
        calcinfo.xml_name = self._XML_FILE_NAME
        calcinfo.json_name = self._JSON_FILE_NAME
        calcinfo.messages_name = self._MESSAGES_FILE_NAME

        #
        # Code information object
        #
        codeinfo = CodeInfo()
        codeinfo.cmdline_params = list(cmdline_params)
        codeinfo.stdin_name = self._INPUT_FILE_NAME
        codeinfo.stdout_name = self._OUTPUT_FILE_NAME
        codeinfo.xml_name = self._XML_FILE_NAME
        codeinfo.json_name = self._JSON_FILE_NAME
        codeinfo.messages_name = self._MESSAGES_FILE_NAME
        codeinfo.code_uuid = code.uuid
        calcinfo.codes_info = [codeinfo]

        # Retrieve by default: the output file, the xml file, and the
        # messages file.
        # If flagbands=True we also add the bands file to the retrieve list!
        # This is extremely important because the parser parses the bands
        # only if aiida.bands is in the retrieve list!!

        calcinfo.retrieve_list = []
        calcinfo.retrieve_list.append(self._OUTPUT_FILE_NAME)
        calcinfo.retrieve_list.append(self._XML_FILE_NAME)
        calcinfo.retrieve_list.append(self._JSON_FILE_NAME)
        calcinfo.retrieve_list.append(self._MESSAGES_FILE_NAME)
        if flagbands:
            calcinfo.retrieve_list.append(self._BANDS_FILE_NAME)

        # Any other files specified in the settings dictionary
        settings_retrieve_list = settings_dict.pop('ADDITIONAL_RETRIEVE_LIST',
                                                   [])
        calcinfo.retrieve_list += settings_retrieve_list

        return calcinfo

    @classmethod
    def _get_linkname_pseudo_prefix(cls):
        """
        The prefix for the name of the link used for the pseudo for kind 'kind'

        :param kind: a string for the atomic kind for which we want
          to get the link name
        """
        return "pseudo_"

    @classmethod
    def _get_linkname_pseudo(cls, kind):
        """
        The name of the link used for the pseudo for kind 'kind'.
        It appends the pseudo name to the pseudo_prefix, as returned by the
        _get_linkname_pseudo_prefix() method.

        :note: if a list of strings is given, the elements are appended
          in the same order, separated by underscores

        :param kind: a string (or list of strings) for the atomic kind(s) for
            which we want to get the link name
        """
        # If it is a list of strings, and not a single string: join them
        # by underscore
        #
        # It might be better to use another character instead of '_'. As it
        # is now, it conflicts with species names of the form Symbol_extra.

        if isinstance(kind, (tuple, list)):
            suffix_string = "_".join(kind)
        elif isinstance(kind, basestring):
            suffix_string = kind
        else:
            raise TypeError("The parameter 'kind' of _get_linkname_pseudo can "
                            "only be a string or a list of strings")
        return "{}{}".format(cls._get_linkname_pseudo_prefix(), suffix_string)

    def use_pseudos_from_family(self, family_name):
        """
        Set the pseudo to use for all atomic kinds, picking pseudos from the
        family with name family_name.

        :note: The structure must already be set.

        :param family_name: the name of the group containing the pseudos
        """
        from collections import defaultdict

        try:
            ##  structure = inputdict.pop(self.get_linkname('structure'))
            structure = self.get_inputs_dict()[self.get_linkname('structure')]
        except AttributeError:
            raise ValueError(
                "Structure is not set yet! Therefore, the method "
                "use_pseudos_from_family cannot automatically set "
                "the pseudos")

        # A dict {kind_name: pseudo_object}
        kind_pseudo_dict = get_pseudos_from_structure(structure, family_name)

        # We have to group the species by pseudo, I use the pseudo PK
        # pseudo_dict will just map PK->pseudo_object
        pseudo_dict = {}
        # Will contain a list of all species of the pseudo with given PK
        pseudo_species = defaultdict(list)

        for kindname, pseudo in kind_pseudo_dict.iteritems():
            pseudo_dict[pseudo.pk] = pseudo
            pseudo_species[pseudo.pk].append(kindname)

        for pseudo_pk in pseudo_dict:
            pseudo = pseudo_dict[pseudo_pk]
            kinds = pseudo_species[pseudo_pk]
            # I set the pseudo for all species, sorting alphabetically
            self.use_pseudo(pseudo, sorted(kinds))

    def _set_parent_remotedata(self, remotedata):
        """
        Used to set a parent remotefolder in the restart of ph.
        """
        from aiida.common.exceptions import ValidationError

        if not isinstance(remotedata, RemoteData):
            raise ValueError('remotedata must be a RemoteData')

        # complain if another remotedata is already found
        input_remote = self.get_inputs(node_type=RemoteData)
        if input_remote:
            raise ValidationError(
                "Cannot set several parent calculation to a "
                "{} calculation".format(self.__class__.__name__))

        self.use_parent_folder(remotedata)

    def create_restart(self, use_output_structure=True, force_restart=True):
        """
        Simple Function to restart a calculation that was not completed
        (for example, due to max walltime reached, or lack of convergence)

        This version effectively requests that the density-matrix file be copied
        from the old calculation's output folder, and sets an fdf option to
        read it upon start. Other possibilites might be given by extra arguments
        in the future (for example, start from scratch)

        Returns a calculation c2, with all links prepared but not stored in DB.
        To submit it simply:
        c2.store_all()
        c2.submit()

        :param bool force_restart: restart also if parent is not in
           FINISHED state (e.g. FAILED, IMPORTED, etc.). Default=True.

        :param bool use_output_structure: if True, the output
           structure of the restarted calculation is used if
           available, rather than its input structure.
           Default=True.

        """
        from aiida.common.datastructures import calc_states

        # Check the calculation's state using ``from_attribute=True`` to
        # correctly handle IMPORTED calculations (so far this applies really only
        # to QE, but it is kept here for future use)

        # In the Siesta plugin, the parser marks non-converged calculations as FAILED
        # when the 'scf-must-converge' and/or 'geometry-must-converge' fdf flags
        # are set.
        # So in practice we will always need to use the force_restart=True
        #
        if self.get_state(from_attribute=True) != calc_states.FINISHED:
            if force_restart:
                pass
            else:
                raise InputValidationError(
                    "Calculation to be restarted must "
                    "be in the {} state. Otherwise, use the force_restart "
                    "flag".format(calc_states.FINISHED))

        # We start here the creation of the new calculation object, using
        # information from the current one
        # What exactly is involved in this 'copy'?
        #
        c2 = self.copy()

        calc_inp = self.get_inputs_dict()

        # The philosophy here is different from that of QE.

        # There is no 'restart' calculation mode in Siesta.
        # We can set the option to read and re-use the DM, if
        # the restart is due to lack of convergence.
        # There is no direct way to read a structure from the working folder
        # and restart a relaxation from it, so in practice we pick up
        # the latest structure from the output node list of the previous calculation.

        # As 'old_inp_dict' is a FDFDict object we can be sure that fdf options are effectively
        # canonicalized, so the following assignment will override any other values
        # for the re-use of DM flag, even if they are in mixed case, etc.
        # Note that options with aliases need to be handled with more care, by
        # setting all possible aliases.

        old_inp_dict = FDFDict(calc_inp['parameters'].get_dict())
        old_inp_dict['dm-use-save-dm'] = True
        c2.use_parameters(ParameterData(dict=old_inp_dict))

        remote_folders = self.get_outputs(node_type=RemoteData)
        if len(remote_folders) != 1:
            raise InputValidationError("More than one output RemoteData found "
                                       "in calculation {}".format(self.pk))
        remote_folder = remote_folders[0]
        c2._set_parent_remotedata(remote_folder)

        # Note that the items to copy from the parent folder are already specified
        # elsewhere in the plugin. We might re-define them here:
        #
        # c2._restart_copy_from = os.path.join(c2._OUTPUT_SUBFOLDER, '*.DM Rho.grid.nc')

        # Could we want to try with a new version of the code?
        c2.use_code(calc_inp['code'])

        # Pseudopotentials
        # This section could be done more cleanly with the following idiom
        # taken from a recent version of the QE plugin:
        #
        #   for linkname, input_node in calc_inp.iteritems():
        #         if isinstance(input_node, UpfData):
        #            c2.add_link_from(input_node, label=linkname)
        #
        # For Siesta, we need to use PsfData (or 'SiestaPseudoData' or similar
        # umbrella class, if we ever allow PsmlData as another kind of pseudo.
        #
        # But we need to make sure that the 'kinds' support is correctly handled
        # (it would be: for example, pseudo_C_Cred is the linkname that assigns
        # the pseudo to the C and Cred kinds).

        for link in calc_inp.keys():
            # Is it a pseudo node?
            if link.startswith(self._get_linkname_pseudo_prefix()):
                # Process the kinds associated to this pseudo
                kindstring = link[len(self._get_linkname_pseudo_prefix()):]
                kinds = kindstring.split('_')

                # Add the pseudo to the new calculation
                the_pseudo = calc_inp[link]
                c2.use_pseudo(the_pseudo, kind=kinds)

        # As explained above, by default we use the latest structure generated
        # by a possibly FAILED relaxation calculation, if available

        if use_output_structure:
            calc_out = self.get_outputs_dict()
            try:
                new_structure = calc_out['output_structure']
                c2.use_structure(new_structure)
            except KeyError:
                c2.use_structure(calc_inp['structure'])
        else:
            c2.use_structure(calc_inp[self.get_linkname('structure')])

        # These are optional...

        # But we could allow an optional argument 'new_basis' (in the form of
        # a basis (ParameterData) object, that would replace the old one. In
        # this case, DM re-use would not be possible.
        # Same for the k-points...
        #
        # In practice, this is probably better done in a workflow, and keep
        # this basic mechanism simple.

        try:
            old_basis = calc_inp['basis']
        except KeyError:
            old_basis = None
        if old_basis is not None:
            c2.use_basis(old_basis)

        try:
            old_kpoints = calc_inp['kpoints']
        except KeyError:
            old_kpoints = None
        if old_kpoints is not None:
            c2.use_kpoints(old_kpoints)

        # If the calculation needs to be restarted, it has probably not reached
        # the 'siesta_analysis' stage, so this link needs to be present.
        # ** Study the workflow implications

        try:
            old_bandskpoints = calc_inp['bandskpoints']
        except KeyError:
            old_bandskpoints = None
        if old_bandskpoints is not None:
            c2.use_bandskpoints(old_bandskpoints)

        try:
            old_settings_dict = calc_inp['settings'].get_dict()
        except KeyError:
            old_settings_dict = {}

        if old_settings_dict:  # if not empty dictionary
            settings = ParameterData(dict=old_settings_dict)
            c2.use_settings(settings)

        return c2


def get_input_data_text(key, val, mapping=None):
    """
    Given a key and a value, return a string (possibly multiline for arrays)
    with the text to be added to the input file.

    :param key: the flag name
    :param val: the flag value. If it is an array, a line for each element
            is produced, with variable indexing starting from 1.
            Each value is formatted using the conv_to_fortran function.
    :param mapping: Optional parameter, must be provided if val is a dictionary.
            It maps each key of the 'val' dictionary to the corresponding
            list index. For instance, if ``key='magn'``,
            ``val = {'Fe': 0.1, 'O': 0.2}`` and ``mapping = {'Fe': 2, 'O': 1}``,
            this function will return the two lines ``magn(1) = 0.2`` and
            ``magn(2) = 0.1``. This parameter is ignored if 'val'
            is not a dictionary.
    """
    from aiida.common.utils import conv_to_fortran
    # I check first the dictionary, because it would also match
    # hasattr(__iter__)
    if isinstance(val, dict):
        if mapping is None:
            raise ValueError("If 'val' is a dictionary, you must provide also "
                             "the 'mapping' parameter")

        list_of_strings = []
        for elemk, itemval in val.iteritems():
            try:
                idx = mapping[elemk]
            except KeyError:
                raise ValueError("Unable to find the key '{}' in the mapping "
                                 "dictionary".format(elemk))

            list_of_strings.append((idx, "  {0}({2}) = {1}\n".format(
                key, conv_to_fortran(itemval), idx)))

        # I first have to resort, then to remove the index from the first
        # column, finally to join the strings
        list_of_strings = zip(*sorted(list_of_strings))[1]
        return "".join(list_of_strings)
    elif hasattr(val, '__iter__'):
        # a list/array/tuple of values
        list_of_strings = [
            "{0}({2})  {1}\n".format(key, conv_to_fortran(itemval), idx + 1)
            for idx, itemval in enumerate(val)
        ]
        return "".join(list_of_strings)
    else:
        # single value
        if key[:6] == '%block':
            bname = key.split()[1]
            b1 = "{0}  {1}".format(key, my_conv_to_fortran(val))
            return b1 + "\n%endblock " + bname + "\n"
        else:
            return "{0}  {1}\n".format(key, my_conv_to_fortran(val))


def my_conv_to_fortran(val):
    """
    Special version to avoid surrounding strings with extra ' '. Otherwise the
    fdf tokenizer will not split values and units, for example.

    :param val: the value to be read and converted to a Fortran-friendly string.
    """
    # Note that bool should come before integer, because a boolean matches also
    # isinstance(...,int)
    if (isinstance(val, bool)):
        if val:
            val_str = '.true.'
        else:
            val_str = '.false.'
    elif (isinstance(val, (int, long))):
        val_str = "{:d}".format(val)
    elif (isinstance(val, float)):
        val_str = ("{:18.10e}".format(val)).replace('e', 'd')
    elif (isinstance(val, basestring)):
        val_str = "{!s}".format(val)
    else:
        raise ValueError("Invalid value passed, accepts only bools, ints, "
                         "floats and strings")

    return val_str


def _uppercase_dict(d, dict_name):
    from collections import Counter

    if isinstance(d, dict):
        new_dict = dict((str(k).upper(), v) for k, v in d.iteritems())
        if len(new_dict) != len(d):

            num_items = Counter(str(k).upper() for k in d.keys())
            double_keys = ",".join([k for k, v in num_items if v > 1])
            raise InputValidationError(
                "Inside the dictionary '{}' there are the following keys that "
                "are repeated more than once when compared case-insensitively: "
                "{}."
                "This is not allowed.".format(dict_name, double_keys))
        return new_dict
    else:
        raise TypeError(
            "_lowercase_dict accepts only dictionaries as argument")
