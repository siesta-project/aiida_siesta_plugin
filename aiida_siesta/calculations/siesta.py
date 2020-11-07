import os
from aiida import orm
from aiida.common import CalcInfo, CodeInfo, InputValidationError
from aiida.common.constants import elements
from aiida.engine import CalcJob
from aiida.orm import Dict, StructureData, BandsData, ArrayData
from aiida_siesta.calculations.tkdict import FDFDict
from aiida_siesta.data.psf import PsfData
from aiida_siesta.data.psml import PsmlData

# See the LICENSE.txt and AUTHORS.txt files.

###################################################################################
## Since aiida 1.0 There is now a clear distinction between Nodes and Processes. ##
## A calculation is now a process and it is treated as a Process class similar   ##
## to the WorkChains. Use of class variables & the input spec is necessary.      ##
###################################################################################


class SiestaCalculation(CalcJob):
    """
    Siesta calculator class for AiiDA.
    """
    #_version = '' Aiida gets the plugin version automatically from package __init__

    ###################################################################
    ## Important distinction between input.spec of the class (can be ##
    ## modified) and pure parameters, stored as class variables only ##
    ###################################################################

    # Parameters stored as class variables
    # 1) Keywords that cannot be set (already canonized by FDFDict)
    # 2) Filepaths of certain outputs
    _aiida_blocked_keywords = [FDFDict.translate_key('system-name'), FDFDict.translate_key('system-label')]
    _aiida_blocked_keywords.append(FDFDict.translate_key('number-of-species'))
    _aiida_blocked_keywords.append(FDFDict.translate_key('number-of-atoms'))
    _aiida_blocked_keywords.append(FDFDict.translate_key('lattice-constant'))
    _aiida_blocked_keywords.append(FDFDict.translate_key('atomic-coordinates-format'))
    _aiida_blocked_keywords.append(FDFDict.translate_key('use-tree-timer'))
    _aiida_blocked_keywords.append(FDFDict.translate_key('xml-write'))
    _aiida_blocked_keywords.append(FDFDict.translate_key('dm-use-save-dm'))
    _aiida_blocked_keywords.append(FDFDict.translate_key('geometry-must-converge'))
    _PSEUDO_SUBFOLDER = './'
    _OUTPUT_SUBFOLDER = './'
    _JSON_FILE = 'time.json'
    _MESSAGES_FILE = 'MESSAGES'

    # Default of the input.spec, it's just default, but user could change the name
    _DEFAULT_PREFIX = 'aiida'
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
        spec.input('structure', valid_type=orm.StructureData, help='Input structure')
        spec.input('kpoints', valid_type=orm.KpointsData, help='Input kpoints', required=False)
        spec.input('bandskpoints', valid_type=orm.KpointsData, help='Input kpoints for bands', required=False)
        spec.input('basis', valid_type=orm.Dict, help='Input basis', required=False)
        spec.input('settings', valid_type=orm.Dict, help='Input settings', required=False)
        spec.input('parameters', valid_type=orm.Dict, help='Input parameters')
        spec.input('parent_calc_folder', valid_type=orm.RemoteData, required=False, help='Parent folder')
        spec.input_namespace('pseudos', valid_type=(PsfData, PsmlData), help='Input pseudo potentials', dynamic=True)

        # Metadada.options host the inputs that are not stored as a separate node, but attached to `CalcJobNode`
        # as attributes. They are optional, since a default is specified, but they might be changed by the user.
        # The first one is siesta specific. The others are defined in the CalcJob, here we change the default.
        spec.input('metadata.options.prefix', valid_type=str, default=cls._DEFAULT_PREFIX)
        spec.inputs['metadata']['options']['input_filename'].default = cls._DEFAULT_INPUT_FILE
        spec.inputs['metadata']['options']['output_filename'].default = cls._DEFAULT_OUTPUT_FILE
        spec.inputs['metadata']['options']['parser_name'].default = 'siesta.parser'

        # Output nodes
        spec.output('output_parameters', valid_type=Dict, required=True, help='The calculation results')
        spec.output('output_structure', valid_type=StructureData, required=False, help='Optional relaxed structure')
        spec.output('bands', valid_type=BandsData, required=False, help='Optional band structure')
        #spec.output('bands_parameters', valid_type=Dict, required=False, help='Optional parameters of bands')
        spec.output('forces_and_stress', valid_type=ArrayData, required=False, help='Optional forces and stress')

        # Option that allows acces through node.res should be existing output node and a Dict
        spec.default_output_node = 'output_parameters'

        # Exit codes for specific errors. Useful for error handeling in workchains
        spec.exit_code(453, 'BANDS_PARSE_FAIL', message='Failure while parsing the bands file')
        spec.exit_code(452, 'BANDS_FILE_NOT_PRODUCED', message='Bands analysis was requested, but file is not present')
        spec.exit_code(450, 'SCF_NOT_CONV', message='Calculation did not reach scf convergence!')
        spec.exit_code(451, 'GEOM_NOT_CONV', message='Calculation did not reach geometry convergence!')
        spec.exit_code(350, 'UNEXPECTED_TERMINATION', message='Statement "Job completed" not detected, unknown error')
        spec.exit_code(449, 'SPLIT_NORM', message='Split_norm parameter too small')
        spec.exit_code(448, 'BASIS_POLARIZ', message='Problems in the polarization of a basis element')

    def prepare_for_submission(self, folder):  # noqa: MC0001  - is mccabe too complex funct -
        """
        Create the input files from the input nodes passed to this instance of the `CalcJob`.

        :param folder: an `aiida.common.folders.Folder` to temporarily write files on disk
        :return: `aiida.common.datastructures.CalcInfo` instance
        """

        ###########################################################################
        # BEGINNING OF INITIAL INPUT CHECK                                        #
        # All input ports that are defined via spec.input are checked by default, #
        # only need to asses their presence in case they are optional.            #
        ###########################################################################

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

        ##############################
        # END OF INITIAL INPUT CHECK #
        ##############################

        # ============== Initialization of some lists ===============
        # List of the file to copy in the folder where the calculation
        # runs, for instance pseudo files
        local_copy_list = []

        # List of files for restart
        remote_copy_list = []

        # ============== Preprocess of input parameters ===============

        input_params = FDFDict(parameters.get_dict())

        # Look for blocked keywords and add the proper values to the dictionary
        for key in input_params:
            if "pao" in key:
                raise InputValidationError(
                    "You can not put PAO options in the parameters input port "
                    "they belong to the basis input port "
                )
            if key in self._aiida_blocked_keywords:
                raise InputValidationError(
                    "You cannot specify explicitly the '{}' flag in the "
                    "input parameters".format(input_params.get_last_untranslated_key(key))
                )

        input_params.update({'system-name': self.inputs.metadata.options.prefix})
        input_params.update({'system-label': self.inputs.metadata.options.prefix})
        input_params.update({'use-tree-timer': 'T'})
        input_params.update({'xml-write': 'T'})
        input_params.update({'number-of-species': len(structure.kinds)})
        input_params.update({'number-of-atoms': len(structure.sites)})
        input_params.update({'geometry-must-converge': 'T'})
        input_params.update({'lattice-constant': '1.0 Ang'})
        input_params.update({'atomic-coordinates-format': 'Ang'})
        # NOTES:
        # 1) The lattice-constant parameter must be 1.0 Ang to impose the units and consider
        #   that the dimenstions of the lattice vectors are already correct with no need of alat.
        #   This breaks the band-k-points "pi/a" option. The use of this option is banned.
        # 2) The implicit coordinate convention of the StructureData class corresponds to the "Ang"
        #   convention in Siesta. That is why "atomic-coordinates-format" is blocked and reset.
        # 3) The Siesta code doesn't raise any warining if the geometry is not converged, unless
        #   the keyword geometry-must-converge is set. That's why it is always added.

        # ============== Preparation of input data ===============

        # ---------------- CELL_PARAMETERS ------------------------
        cell_parameters_card = "%block lattice-vectors\n"
        for vector in structure.cell:
            cell_parameters_card += ("{0:18.10f} {1:18.10f} {2:18.10f}" "\n".format(*vector))
        cell_parameters_card += "%endblock lattice-vectors\n"

        # --------------ATOMIC_SPECIES & PSEUDOS-------------------
        # Subfolder that will contain the pseudopotentials and output data
        folder.get_subfolder(self._PSEUDO_SUBFOLDER, create=True)
        folder.get_subfolder(self._OUTPUT_SUBFOLDER, create=True)
        atomic_species_card_list = []
        # Dictionary to get the atomic number of a given element
        #pylint: disable=consider-using-dict-comprehension
        datmn = dict([(v['symbol'], k) for k, v in elements.items()])
        spind = {}
        spcount = 0
        for kind in structure.kinds:
            spcount += 1  # species count
            spind[kind.name] = spcount
            atomic_species_card_list.append(
                "{0:5} {1:5} {2:5}\n".format(spind[kind.name], datmn[kind.symbol], kind.name.rjust(6))
            )
            psp = pseudos[kind.name]
            # Add this pseudo file to the list of files to copy, with the appropiate name.
            # In the case of sub-species (different kind.name but same kind.symbol, e.g.,
            # 'C_surf', sharing the same pseudo with 'C'), we copy the file ('C.psf')
            # twice, once as 'C.psf', and once as 'C_surf.psf'. This is required by Siesta.
            # It is passed in form of a list of tuples with format ('node_uuid', 'filename',
            # relativedestpath'). We probably should be pre-pending 'self._PSEUDO_SUBFOLDER'
            # in the last slot, for generality, even if is not necessary for siesta.
            if isinstance(psp, PsfData):
                local_copy_list.append((psp.uuid, psp.filename, kind.name + ".psf"))
            elif isinstance(psp, PsmlData):
                local_copy_list.append((psp.uuid, psp.filename, kind.name + ".psml"))
            else:
                pass
        atomic_species_card_list = (["%block chemicalspecieslabel\n"] + list(atomic_species_card_list))
        atomic_species_card = "".join(atomic_species_card_list)
        atomic_species_card += "%endblock chemicalspecieslabel\n"
        # Free memory
        del atomic_species_card_list

        # --------------------- ATOMIC_POSITIONS -----------------------
        atomic_positions_card_list = ["%block atomiccoordinatesandatomicspecies\n"]
        countatm = 0
        for site in structure.sites:
            countatm += 1
            atomic_positions_card_list.append(
                "{0:18.10f} {1:18.10f} {2:18.10f} {3:4} {4:6} {5:6}\n".format(
                    site.position[0], site.position[1], site.position[2], spind[site.kind_name],
                    site.kind_name.rjust(6), countatm
                )
            )
        atomic_positions_card = "".join(atomic_positions_card_list)
        del atomic_positions_card_list  # Free memory
        atomic_positions_card += "%endblock atomiccoordinatesandatomicspecies\n"

        # -------------------- K-POINTS ----------------------------
        # It is optional, if not specified, gamma point only is performed (default of siesta)
        if kpoints is not None:
            # There is not yet support for the 'kgrid-cutoff' option in Siesta. Only mesh accepted
            try:
                mesh, offset = kpoints.get_kpoints_mesh()
            except AttributeError:
                raise InputValidationError("K-point sampling for scf " "must be given in mesh form")
            kpoints_card_list = ["%block kgrid_monkhorst_pack\n"]
            # This would fail if kpoints is not a mash (for the case of a list),
            # since in that case 'offset' is undefined.
            kpoints_card_list.append("{0:6} {1:6} {2:6} {3:18.10f}\n".format(mesh[0], 0, 0, offset[0]))
            kpoints_card_list.append("{0:6} {1:6} {2:6} {3:18.10f}\n".format(0, mesh[1], 0, offset[1]))
            kpoints_card_list.append("{0:6} {1:6} {2:6} {3:18.10f}\n".format(0, 0, mesh[2], offset[2]))
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
            #first, we check that the user constracted the kpoints using the cell
            #of the input structure, and not a random cell. This helps parsing
            kpcell = bandskpoints.get_attribute("cell", None)
            if kpcell:
                if kpcell != structure.cell:
                    raise ValueError(
                        'The cell used for `bandskpoints` must be the same of the input structure.'
                        'Alternatively do not set any cell to the bandskpoints.'
                    )
            #second we rise a warning about consequences when the cell is relaxed
            var_cell_keys = [FDFDict.translate_key("md-variable-cell"), FDFDict.translate_key("md-constant-volume")]
            var_cell_keys.append(FDFDict.translate_key("md-relax-cell-only"))
            for key in input_params:
                if key in var_cell_keys:
                    logline = (
                        "Requested calculation of bands after a relaxation with variable cell! " +
                        "If the symmetry of the cell will change, the kpoints path for bands will be wrong. " +
                        "It is suggested to use the `BandGapWorkChain` instead."
                    )
                    if isinstance(input_params[key], str):
                        if FDFDict.translate_key(input_params[key]) in ["t", "true", "yes"]:
                            self.logger.warning(logline)
                            break
                    else:
                        if input_params[key] is True:
                            self.logger.warning(logline)
                            break
            #the band line scale
            bandskpoints_card_list = ["BandLinesScale ReciprocalLatticeVectors\n"]
            #set the BandPoints
            if bandskpoints.labels is None:
                bandskpoints_card_list.append("%block BandPoints\n")
                for kpo in bandskpoints.get_kpoints():
                    bandskpoints_card_list.append("{0:8.3f} {1:8.3f} {2:8.3f} \n".format(kpo[0], kpo[1], kpo[2]))
                fbkpoints_card = "".join(bandskpoints_card_list)
                fbkpoints_card += "%endblock BandPoints\n"
            #set the BandLines
            else:
                bandskpoints_card_list.append("%block BandLines\n")
                savindx = []
                listforbands = bandskpoints.get_kpoints()
                for indx, label in bandskpoints.labels:
                    savindx.append(indx)
                rawindex = 0
                for indx, label in bandskpoints.labels:
                    rawindex = rawindex + 1
                    x, y, z = listforbands[indx]
                    if rawindex == 1:
                        bandskpoints_card_list.append(
                            "{0:3} {1:8.3f} {2:8.3f} {3:8.3f} {4:1} \n".format(1, x, y, z, label)
                        )
                    else:
                        bandskpoints_card_list.append(
                            "{0:3} {1:8.3f} {2:8.3f} {3:8.3f} {4:1} \n".format(
                                indx - savindx[rawindex - 2], x, y, z, label
                            )
                        )
                fbkpoints_card = "".join(bandskpoints_card_list)
                fbkpoints_card += "%endblock BandLines\n"
            del bandskpoints_card_list

        # ================ Operations for restart =======================
        # The presence of a 'parent_calc_folder' input node signals that we want to
        # get something from there, as indicated in the self._restart_copy_from attribute.
        # In Siesta's case, for now, just the density-matrix file is copied
        # to the current calculation's working folder.
        # ISSUE: Is this mechanism flexible enough? An alternative would be to
        # pass the information about which file(s) to copy in the metadata.options dictionary
        if parent_calc_folder is not None:
            remote_copy_list.append((
                parent_calc_folder.computer.uuid,
                os.path.join(parent_calc_folder.get_remote_path(), self._restart_copy_from), self._restart_copy_to
            ))
            input_params.update({'dm-use-save-dm': "T"})

        # ====================== FDF file creation ========================

        # To have easy access to inputs metadata options
        metadataoption = self.inputs.metadata.options

        # input_filename = self.inputs.metadata.options.input_filename
        input_filename = folder.get_abs_path(metadataoption.input_filename)

        with open(input_filename, 'w') as infile:
            # here print keys and values tp file

            for k, v in sorted(input_params.get_filtered_items()):
                infile.write("%s %s\n" % (k, v))

            # Basis set info is processed just like the general
            # parameters section. Some discipline is needed to
            # put any basis-related parameters (including blocks)
            # in the basis dictionary in the input script.
            if basis is not None:
                infile.write("#\n# -- Basis Set Info follows\n#\n")
                for k, v in basis.get_dict().items():
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
            infile.write("max.walltime {}\n".format(metadataoption.max_wallclock_seconds))

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
        xml_file = str(metadataoption.prefix) + ".xml"
        bands_file = str(metadataoption.prefix) + ".bands"
        calcinfo.retrieve_list.append(metadataoption.output_filename)
        calcinfo.retrieve_list.append(xml_file)
        calcinfo.retrieve_list.append(self._JSON_FILE)
        calcinfo.retrieve_list.append(self._MESSAGES_FILE)
        if bandskpoints is not None:
            calcinfo.retrieve_list.append(bands_file)
        # Any other files specified in the settings dictionary
        settings_retrieve_list = settings_dict.pop('ADDITIONAL_RETRIEVE_LIST', [])
        calcinfo.retrieve_list += settings_retrieve_list

        return calcinfo

    @classmethod
    def inputs_generator(cls):  # pylint: disable=no-self-argument,no-self-use
        from aiida_siesta.utils.inputs_generators import SiestaCalculationInputsGenerator
        return SiestaCalculationInputsGenerator(cls)


def _uppercase_dict(indic, dict_name):
    from collections import Counter

    if not isinstance(indic, dict):
        raise TypeError("_uppercase_dict accepts only dictionaries as argument")

    new_dict = dict((str(k).upper(), v) for k, v in indic.items())

    if len(new_dict) != len(indic):
        num_items = Counter(str(k).upper() for k in indic.keys())
        double_keys = ",".join([k for k, v in num_items if v > 1])
        raise InputValidationError(
            "Inside the dictionary '{}' there are the following keys that "
            "are repeated more than once when compared case-insensitively: "
            "{}."
            "This is not allowed.".format(dict_name, double_keys)
        )

    return new_dict
