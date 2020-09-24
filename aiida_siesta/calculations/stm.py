import os
from aiida import orm
from aiida.common import CalcInfo, CodeInfo, InputValidationError
from aiida.engine import CalcJob
from aiida.orm import Dict, ArrayData

# See the LICENSE.txt and AUTHORS.txt files.


class STMCalculation(CalcJob):
    """
    Plugin for the "plstm" program in the Siesta distribution, which
    takes and .LDOS or .RHO file and generates a plot file to simulate
    an STM image.
    """

    _OUTPUT_SUBFOLDER = './'

    # Default input and output files
    _DEFAULT_INPUT_FILE = 'stm.in'
    _DEFAULT_OUTPUT_FILE = 'stm.out'

    # in restarts, it will copy from the parent the following
    _restart_copy_from = os.path.join(_OUTPUT_SUBFOLDER, '*.LDOS')

    # in restarts, it will copy the previous folder in the following one
    _restart_copy_to = _OUTPUT_SUBFOLDER

    @classmethod
    def define(cls, spec):
        super(STMCalculation, cls).define(spec)

        spec.input('code', valid_type=orm.Code, help='Input code')
        spec.input('settings', valid_type=orm.Dict, help='Input settings', required=False)
        spec.input(
            'spin_option',
            valid_type=orm.Str,
            default=lambda: orm.Str("q"),
            help='Spin option, follows plstm sintax: '
            '"q" no spin, "s" total spin, "x","y","z" the three '
            'spin components (only available in non-collinear case)'
        )
        spec.input('mode', valid_type=orm.Str, help='Allowed values are "constant-height" or "constant-current"')
        spec.input('value', valid_type=orm.Float, help='Value of height in Ang or value of current in e/bohr**3')
        spec.input('ldos_folder', valid_type=orm.RemoteData, required=True, help='Parent folder')

        # Metadata, defined in CalCJob, here we modify default
        spec.inputs['metadata']['options']['input_filename'].default = cls._DEFAULT_INPUT_FILE
        spec.inputs['metadata']['options']['output_filename'].default = cls._DEFAULT_OUTPUT_FILE
        spec.inputs['metadata']['options']['parser_name'].default = 'siesta.stm'

        #Outputs
        spec.output('stm_array', valid_type=ArrayData, required=True, help='The contour data for the STM image')
        spec.output(
            'output_parameters',
            valid_type=Dict,
            required=True,
            help='Other outputs, for the moment only parser version and name of .STM file'
        )

        # exit codes
        spec.exit_code(
            100, 'ERROR_NO_RETRIEVED_FOLDER', message='The retrieved folder data node could not be accessed.'
        )
        spec.exit_code(101, 'ERROR_OUTPUT_PLOT_MISSING', message='The retrieved folder does not contain a CH.STM file')
        spec.exit_code(102, 'ERROR_OUTPUT_PLOT_READ', message='The .STM file can not be read')
        spec.exit_code(102, 'ERROR_CREATION_STM_ARRAY', message='The array containing the STM data can not be produced')

    def prepare_for_submission(self, folder):  # noqa: MC0001  - is mccabe too complex funct -
        """
        Create the input files from the input nodes passed to this instance of the `CalcJob`.

        :param folder: an `aiida.common.folders.Folder` to temporarily write files on disk
        :return: `aiida.common.datastructures.CalcInfo` instance
        """

        code = self.inputs.code
        ldos_folder = self.inputs.ldos_folder

        allowedmodes = ["constant-height", "constant-current"]
        mode = self.inputs.mode
        if mode.value not in allowedmodes:
            raise ValueError("The allowed options for the port 'mode' are {}".format(allowedmodes))

        value = self.inputs.value
        spin_option = self.inputs.spin_option

        if 'settings' in self.inputs:
            settings = self.inputs.settings.get_dict()
            settings_dict = _uppercase_dict(settings, dict_name='settings')
        else:
            settings_dict = {}

        # List of the files to copy in the folder where the calculation
        # runs, for instance pseudo files
        local_copy_list = []

        # List of files for restart
        remote_copy_list = []

        # ============== Creation of input file =========================
        #Input file is only necessary for the old versions of plstm
        #For the new versions, all is done through command line (next section)

        # To have easy access to inputs metadata options
        metadataoption = self.inputs.metadata.options

        # input_filename access
        input_filename = folder.get_abs_path(metadataoption.input_filename)

        # Getting the prefix from ldos_folder
        if "prefix" in ldos_folder.creator.attributes:
            prefix = str(ldos_folder.creator.get_attribute("prefix"))
        else:
            self.report("No prefix detected from the remote folder, set 'aiida' as prefix")
            prefix = "aiida"
        ldosfile = prefix + ".LDOS"

        # Convert height to bohr...
        if mode.value == "constant-height":
            vvalue = value.value / 0.529177
        else:
            vvalue = value.value
        with open(input_filename, 'w') as infile:
            infile.write("{}\n".format(prefix))
            infile.write("ldos\n")
            infile.write("{}\n".format(mode.value))
            infile.write("{0:.5f}\n".format(vvalue))
            infile.write("unformatted\n")

        # ====================== Code and Calc info ========================
        # Code information object and Calc information object are now
        # only used to set up the CMDLINE (the bash line that launches siesta)
        # and to set up the list of files to retrieve.
        # The presence of a 'ldos_folder' is mandatory, to get the LDOS file
        # as indicated in the self._restart_copy_from attribute.
        # (this is not technically a restart, though)

        # It will be copied to the current calculation's working folder.
        remote_copy_list.append((
            ldos_folder.computer.uuid, os.path.join(ldos_folder.get_remote_path(),
                                                    self._restart_copy_from), self._restart_copy_to
        ))

        # Empty command line by default. Why use 'pop' ?
        cmdline_params = settings_dict.pop('CMDLINE', [])

        # Code information object. Sets the command line
        codeinfo = CodeInfo()
        if mode.value == "constant-height":
            cmdline_params = (list(cmdline_params) + ['-z', '{0:.5f}'.format(vvalue)])
        else:
            cmdline_params = (list(cmdline_params) + ['-i', '{0:.5f}'.format(vvalue)])
        if spin_option.value != "q":
            cmdline_params = (list(cmdline_params) + ['-s', str(spin_option.value), ldosfile])
        else:
            cmdline_params = (list(cmdline_params) + [ldosfile])
        codeinfo.cmdline_params = list(cmdline_params)
        codeinfo.stdin_name = metadataoption.input_filename
        codeinfo.stdout_name = metadataoption.output_filename
        codeinfo.code_uuid = code.uuid

        # Calc information object. Important for files to copy, retrieve, etc
        calcinfo = CalcInfo()
        calcinfo.uuid = str(self.uuid)
        calcinfo.local_copy_list = local_copy_list
        calcinfo.remote_copy_list = remote_copy_list
        calcinfo.codes_info = [codeinfo]
        #Next three are useless in my opinion, as they can be access
        #through calcinfo.codes_info and not used interally by AiiDA
        if cmdline_params:
            calcinfo.cmdline_params = list(cmdline_params)
        calcinfo.stdin_name = metadataoption.input_filename
        calcinfo.stdout_name = metadataoption.output_filename
        # Retrieve by default: the output file and the plot file,
        # Some logic to understand which is the plot file will be in
        # parser, here we put to retrieve every file ending in *CH.STM
        calcinfo.retrieve_list = []
        calcinfo.retrieve_list.append(metadataoption.output_filename)
        calcinfo.retrieve_list.append("*.STM")
        # Any other files specified in the settings dictionary
        settings_retrieve_list = settings_dict.pop('ADDITIONAL_RETRIEVE_LIST', [])
        calcinfo.retrieve_list += settings_retrieve_list

        return calcinfo


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
