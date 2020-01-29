# -*- coding: utf-8 -*-
from __future__ import absolute_import

import os

import six
from aiida import orm
from aiida.common import CalcInfo, CodeInfo, InputValidationError
from aiida.common.constants import elements
from aiida.engine import CalcJob
from aiida.orm import Dict, StructureData, BandsData, ArrayData

from .tkdict import FDFDict
from aiida_siesta.data.psf import PsfData
from aiida_siesta.data.psml import PsmlData
# See the LICENSE.txt and AUTHORS.txt files.

###################################################################
## Few comments about aiida 1.0:                                 ##
## There is now a clear distinction between Nodes and Processes  ##
## A calculation is now a process and it is treated as a Process ##
## class similar to the WorkChains. Use of class variables &     ##
## the input spec is necessary                                   ##
###################################################################


class SiestaCalculation(CalcJob):
    """
    Siesta calculator class for AiiDA.
    """
    _siesta_plugin_version = '1.0.1'

    ###########################################################
    ## Important distinction between input.spec of the class ##
    ## (can be modified) and pure parameters, stored as      ##
    ## class variables only                                  ##
    ###########################################################

    # Parameters stored as class variables
    # 1) Keywords that cannot be set (need to canonoze this?)
    # 2) Filepaths of certain outputs
    _aiida_blocked_keywords = ['system-name', 'system-label']
    _aiida_blocked_keywords.append('number-of-species')
    _aiida_blocked_keywords.append('number-of-atoms')
    _aiida_blocked_keywords.append('latticeconstant')
    _aiida_blocked_keywords.append('lattice-constant')
    _aiida_blocked_keywords.append('atomic-coordinates-format')
    _aiida_blocked_keywords.append('atomiccoordinatesformat')
    _aiida_blocked_keywords.append('use-tree-timer')
    _aiida_blocked_keywords.append('xml-write')
    _aiida_blocked_keywords.append('dm-use-save-dm')
    _aiida_blocked_keywords.append('dmusesavedm')
    _PSEUDO_SUBFOLDER = './'
    _OUTPUT_SUBFOLDER = './'
    _PREFIX = 'aiida'
    _DEFAULT_XML_FILE = 'aiida.xml'
    _DEFAULT_JSON_FILE = 'time.json'
    _DEFAULT_MESSAGES_FILE = 'MESSAGES'
    _DEFAULT_BANDS_FILE = 'aiida.bands'


    # Default of the input.spec, it's just default, but user
    # could change the name
    _DEFAULT_INPUT_FILE = 'aiida.fdf'
    _DEFAULT_OUTPUT_FILE = 'aiida.out'

    # in restarts, it will copy from the parent the following
    # (fow now, just the density matrix file)
    _restart_copy_from = os.path.join(_OUTPUT_SUBFOLDER, '*.DM')

    # in restarts, it will copy the previous folder in the following one
    _restart_copy_to = _OUTPUT_SUBFOLDER

    @classmethod
    def define(cls, spec):
        super(SiestaCalculation, cls).define(spec)

        # Input nodes
        spec.input('code', valid_type=orm.Code, help='Input code')
        spec.input('structure',
                   valid_type=orm.StructureData,
                   help='Input structure')
        spec.input('kpoints',
                   valid_type=orm.KpointsData,
                   help='Input kpoints',
                   required=False)
        spec.input('bandskpoints',
                   valid_type=orm.KpointsData,
                   help='Input kpoints for bands',
                   required=False)
        spec.input('basis',
                   valid_type=orm.Dict,
                   help='Input basis',
                   required=False)
        spec.input('settings',
                   valid_type=orm.Dict,
                   help='Input settings',
                   required=False)
        spec.input('parameters',
                   valid_type=orm.Dict,
                   help='Input parameters')
        spec.input('parent_calc_folder',
                   valid_type=orm.RemoteData,
                   required=False,
                   help='Parent folder')
        spec.input_namespace(
                  'pseudos',
                  valid_type=(PsfData, PsmlData),
                  help='Input pseudo potentials',
                  dynamic=True)

        # Metadada.options host the inputs that are not stored
        # as a separate node, but attached to `CalcJobNode`
        # as attributes. They are optional, since a default is 
        # specified, but they might be changed by the user.
        spec.input('metadata.options.input_filename',
                   valid_type=six.string_types,
                   default=cls._DEFAULT_INPUT_FILE)
        spec.input('metadata.options.output_filename',
                   valid_type=six.string_types,
                   default=cls._DEFAULT_OUTPUT_FILE)
        spec.input('metadata.options.parser_name',
                   valid_type=six.string_types,
                   default='siesta.parser')

        # Output nodes
        spec.output('output_parameters',
                    valid_type=Dict,
                    required=True,
                    help='The calculation results')
        spec.output('output_structure',
                    valid_type=StructureData,
                    required=False,
                    help='Optional relaxed structure')
        # Note name change: bands_array --> bands
        spec.output('bands',
                    valid_type=BandsData,
                    required=False,
                    help='Optional band structure')
        # I don't know why the bands parameters are parsed as BandsData already 
        # contains the kpoints (Emanuele)
        # AG: Agreed, this will go soon.
        spec.output('bands_parameters',
                    valid_type=Dict,
                    required=False,
                    help='Optional parameters of bands')
        spec.output('forces_and_stress',
                    valid_type=ArrayData,
                    required=False,
                    help='Optional forces and stress')

        # Option that allows acces through node.res
        # should be existing output node and a Dict
        spec.default_output_node = 'output_parameters'

        # Error handeling
        spec.exit_code(
            100,
            'ERROR_NO_RETRIEVED_FOLDER',
            message='The retrieved folder data node could not be accessed.')
        spec.exit_code(
            120,
            'SCF_NOT_CONV',
            message='Calculation did not reach scf convergence!')
        spec.exit_code(
            130,
            'GEOM_NOT_CONV',
            message='Calculation did not reach geometry convergence!')


#to DO SOON: improve help for pseudo.

    def prepare_for_submission(self, tempfolder):
        """
        Create the input files from the input nodes passed to this instance of the `CalcJob`.

        :param tempfolder: an `aiida.common.folders.Folder` to temporarily write files on disk
        :return: `aiida.common.datastructures.CalcInfo` instance
        """

        #####################################################
        # BEGINNING OF INITIAL INPUT CHECK                  #
        # All input ports that are defined via spec.input   #
        # are checked by default, only need to asses their  #
        # presence in case they are optional                #
        #####################################################

        code = self.inputs.code
        structure = self.inputs.structure
        parameters = self.inputs.parameters

        if 'kpoints' in self.inputs:
            kpoints = self.inputs.kpoints
        else:
            kpoints = None

        if 'basis' in self.inputs:
            basis = self.inputs.basis
        else:
            basis = None

        if 'settings' in self.inputs:
            settings = self.inputs.settings.get_dict()
            settings_dict = _uppercase_dict(settings, dict_name='settings')
        else:
            settings_dict = {}

        if 'bandskpoints' in self.inputs:
            bandskpoints = self.inputs.bandskpoints
        else:
            bandskpoints = None

        if 'parent_calc_folder' in self.inputs:
            parent_calc_folder = self.inputs.parent_calc_folder
        else:
            parent_calc_folder = None

        pseudos = self.inputs.pseudos
        kinds = [kind.name for kind in structure.kinds]
        if set(kinds) != set(pseudos.keys()):
            raise ValueError(
                'Mismatch between the defined pseudos and the list of kinds of the structure.\n',
                'Pseudos: {} \n'.format(', '.join(list(pseudos.keys()))),
                'Kinds: {}'.format(', '.join(list(kinds))),
            )

        # List of the file to copy in the folder where the calculation
        # runs, for instance pseudo files
        local_copy_list = []

        # List of files for restart
        remote_copy_list = []

        ##############################
        # END OF INITIAL INPUT CHECK #
        ##############################

        # ============== Preprocess of input parameters ===============
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

        # ---------------- CELL_PARAMETERS ------------------------
        cell_parameters_card = "%block lattice-vectors\n"
        for vector in structure.cell:
            cell_parameters_card += ("{0:18.10f} {1:18.10f} {2:18.10f}"
                                     "\n".format(*vector))
        cell_parameters_card += "%endblock lattice-vectors\n"

        # --------------ATOMIC_SPECIES & PSEUDOS-------------------
        # I create the subfolder that will contain the pseudopotentials
        tempfolder.get_subfolder(self._PSEUDO_SUBFOLDER, create=True)
        # I create the subfolder with the output data
        tempfolder.get_subfolder(self._OUTPUT_SUBFOLDER, create=True)

        atomic_species_card_list = []

        # Dictionary to get the atomic number of a given element
        datmn = dict([(v['symbol'], k) for k, v in six.iteritems(elements)])

        spind = {}
        spcount = 0
        for kind in structure.kinds:

            spcount += 1  # species count
            spind[kind.name] = spcount
            atomic_species_card_list.append("{0:5} {1:5} {2:5}\n".format(
                spind[kind.name], datmn[kind.symbol], kind.name.rjust(6)))

            ps = pseudos[kind.name]

            # Add this pseudo file to the list of files to copy, with
            # the appropiate name. In the case of sub-species
            # (different kind.name but same kind.symbol, e.g.,
            # 'C_surf', sharing the same pseudo with 'C'), we will
            # copy the file ('C.psf') twice, once as 'C.psf', and once
            # as 'C_surf.psf'.  This is required by Siesta.

            # ... list of tuples with format ('node_uuid', 'filename', relativedestpath')
            # We probably should be pre-pending 'self._PSEUDO_SUBFOLDER' in the
            # last slot, for generality...
            if isinstance(ps, PsfData):
                local_copy_list.append((ps.uuid, ps.filename,
                                        kind.name + ".psf"))
            elif isinstance(ps, PsmlData):
                local_copy_list.append((ps.uuid, ps.filename,
                                        kind.name + ".psml"))
            else:
                pass
                
        atomic_species_card_list = (["%block chemicalspecieslabel\n"] +
                                    list(atomic_species_card_list))
        atomic_species_card = "".join(atomic_species_card_list)
        atomic_species_card += "%endblock chemicalspecieslabel\n"
        # Free memory
        del atomic_species_card_list

        # --------------------- ATOMIC_POSITIONS -----------------------
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

        # -------------------- K-POINTS ----------------------------
        # It is optional, if not specified, gamma point only is performed,
        # this is default of siesta
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

        # ----------------- K-POINTS-FOR-BANDS ----------------------
        #Two possibility are supported in Siesta: BandLines ad BandPoints
        #At the moment the user can't choose directly one of the two options
        #BandsLine is set automatically if bandskpoints has labels,
        #BandsPoints if bandskpoints has no labels
        #BandLinesScale =pi/a is not supported at the moment because currently
        #a=1 always. BandLinesScale ReciprocalLatticeVectors is always set
        if bandskpoints is not None:
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

        # ================ Operations for restart =======================

        # The presence of a 'parent_calc_folder' input node signals
        # that we want to get something from there, as indicated in the
        # self._restart_copy_from attribute.
        # In Siesta's case, for now, it is just the density-matrix file
        #
        # It will be copied to the current calculation's working folder.

        # NOTE: This mechanism is not flexible enough.
        # Maybe we should pass the information about which file(s) to
        # copy in the metadata 'options' dictionary
        if parent_calc_folder is not None:
            remote_copy_list.append(
                (parent_calc_folder.computer.uuid,
                 os.path.join(parent_calc_folder.get_remote_path(),
                              self._restart_copy_from), self._restart_copy_to))

            input_params.update({'dm-use-save-dm': "T"})

        # ====================== FDF file creation ========================

        # To have easy access to inputs metadata options
        metadataoption = self.inputs.metadata.options

        # input_filename = self.inputs.metadata.options.input_filename
        input_filename = tempfolder.get_abs_path(metadataoption.input_filename)

        with open(input_filename, 'w') as infile:
            # here print keys and values tp file

            # for k, v in sorted(six.iteritems(input_params)):
            for k, v in sorted(input_params.get_filtered_items()):
                infile.write("%s %s\n" % (k, v))

            # Basis set info is processed just like the general
            # parameters section. Some discipline is needed to
            # put any basis-related parameters (including blocks)
            # in the basis dictionary in the input script.
            #
            if basis is not None:
                infile.write("#\n# -- Basis Set Info follows\n#\n")
                for k, v in six.iteritems(basis.get_dict()):
                    infile.write("%s %s\n" % (k, v))

            # Write previously generated cards now
            infile.write("#\n# -- Structural Info follows\n#\n")
            infile.write(atomic_species_card)
            infile.write(cell_parameters_card)
            infile.write(atomic_positions_card)
            if kpoints is not None:
                infile.write("#\n# -- K-points Info follows\n#\n")
                infile.write(kpoints_card)
            if bandskpoints is not None:
                infile.write("#\n# -- Bandlines/Bandpoints Info follows\n#\n")
                infile.write(fbkpoints_card)

            # Write max wall-clock time
            infile.write("#\n# -- Max wall-clock time block\n#\n")
            infile.write("max.walltime {}".format(
                metadataoption.max_wallclock_seconds))

        # ====================== Code and Calc info ========================
        # Code information object and Calc information object are now
        # only used to set up the CMDLINE (the bash line that launches siesta)
        # and to set up the list of files to retrieve.

        cmdline_params = settings_dict.pop('CMDLINE', [])
        
        codeinfo = CodeInfo()
        codeinfo.cmdline_params = list(cmdline_params)
        codeinfo.stdin_name = metadataoption.input_filename
        codeinfo.stdout_name = metadataoption.output_filename
        codeinfo.code_uuid = code.uuid
        
        calcinfo = CalcInfo()
        calcinfo.uuid = str(self.uuid)
        if cmdline_params:
            calcinfo.cmdline_params = list(cmdline_params)
        calcinfo.local_copy_list = local_copy_list
        calcinfo.remote_copy_list = remote_copy_list
        calcinfo.stdin_name = metadataoption.input_filename
        calcinfo.stdout_name = metadataoption.output_filename
        calcinfo.codes_info = [codeinfo]
        # Retrieve by default: the output file, the xml file, the
        # messages file, and the json timing file.
        # If bandskpoints, also the bands file is added to the retrieve list.
        calcinfo.retrieve_list = []
        calcinfo.retrieve_list.append(metadataoption.output_filename)
        calcinfo.retrieve_list.append(self._DEFAULT_XML_FILE) 
        calcinfo.retrieve_list.append(self._DEFAULT_JSON_FILE) 
        calcinfo.retrieve_list.append(self._DEFAULT_MESSAGES_FILE) 
        if bandskpoints is not None:
            calcinfo.retrieve_list.append(self._DEFAULT_BANDS_FILE) 
        # Any other files specified in the settings dictionary
        settings_retrieve_list = settings_dict.pop('ADDITIONAL_RETRIEVE_LIST',
                                                   [])
        calcinfo.retrieve_list += settings_retrieve_list

        return calcinfo


def _uppercase_dict(d, dict_name):
    from collections import Counter

    if isinstance(d, dict):
        new_dict = dict((str(k).upper(), v) for k, v in six.iteritems(d))
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
