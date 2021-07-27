import os
from aiida import orm
from aiida.common import CalcInfo, CodeInfo
from aiida.engine import CalcJob
from aiida.orm import Dict, ArrayData

# See the LICENSE.txt and AUTHORS.txt files.


def validate_mode(value, _):
    """
    Validate mode input port.
    """
    if value:
        allowedmodes = ["constant-height", "constant-current"]
        if value.value not in allowedmodes:
            return f"The allowed options for the port 'mode' are {allowedmodes}."


def validate_spin(value, _):
    """
    Validate spin_option input port.
    """
    if value:
        allowedspins = ["q", "s", "x", "y", "z"]
        if value.value not in allowedspins:
            return f"The allowed options for the port 'spin_option' are {allowedspins}."


class STMCalculation(CalcJob):
    """
    Plugin for the "plstm" program in the Siesta distribution, which takes the .LDOS file and
    generates a plot file to simulate an STM image.
    It supports both the old "plstm" versions (inputs in a files) and the new ones (inputs in the command
    line). Spin options are supported only in recent "plstm" versions, therefore ignored otherwise.
    """

    # Default input and output files
    _DEFAULT_INPUT_FILE = 'stm.in'
    _DEFAULT_OUTPUT_FILE = 'stm.out'

    # in restarts, it will copy from the parent the following
    _restart_copy_from = os.path.join('./', '*.LDOS')

    # in restarts, it will copy the previous folder in the following one
    _restart_copy_to = './'

    @classmethod
    def define(cls, spec):
        super(STMCalculation, cls).define(spec)

        spec.input('code', valid_type=orm.Code, help='Input code')
        spec.input('settings', valid_type=orm.Dict, help='Input settings', required=False)
        spec.input(
            'spin_option',
            valid_type=orm.Str,
            default=lambda: orm.Str("q"),
            help='Spin option follows plstm sintax: "q" no spin, "s" total spin, "x","y","z" the three spin components',
            validator=validate_spin
        )
        spec.input(
            'mode',
            valid_type=orm.Str,
            help='Allowed values are "constant-height" or "constant-current"',
            validator=validate_mode
        )
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
            help='For the moment only parser version and name of .STM file'
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
        value = self.inputs.value
        spin_option = self.inputs.spin_option
        mode = self.inputs.mode

        # As internal convention, the keys of the settings dict are uppercase
        if 'settings' in self.inputs:
            settings = self.inputs.settings.get_dict()
            settings_dict = {str(k).upper(): v for (k, v) in settings.items()}
        else:
            settings_dict = {}

        # List of files for restart
        remote_copy_list = []

        # ======================== Creation of input file =========================
        # Input file is only necessary for the old versions of plstm.
        # For the new versions, all is done through command line (next section).

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
            infile.write(f"{prefix}\n")
            infile.write("ldos\n")
            infile.write(f"{mode.value}\n")
            infile.write(f"{vvalue:.5f}\n")
            infile.write("unformatted\n")

        # ============================== Code and Calc info ===============================
        # Code information object is used to set up the the bash line that launches siesta
        # (CMDLINE and input output files).
        # Calc information object is to set up thee list of files to retrieve.
        # The presence of a 'ldos_folder' is mandatory, to get the LDOS file as indicated in
        # the self._restart_copy_from attribute. (this is not technically a restart, though)

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
        calcinfo.local_copy_list = []  #No local files to copy (no pseudo for instance)
        calcinfo.remote_copy_list = remote_copy_list
        calcinfo.codes_info = [codeinfo]
        # Retrieve by default: the output file and the plot file. Some logic to understand which
        # is the plot file will be in parser, here we put to retrieve every file ending in *.STM
        calcinfo.retrieve_list = []
        calcinfo.retrieve_list.append(metadataoption.output_filename)
        calcinfo.retrieve_list.append("*.STM")
        # Any other files specified in the settings dictionary
        settings_retrieve_list = settings_dict.pop('ADDITIONAL_RETRIEVE_LIST', [])
        calcinfo.retrieve_list += settings_retrieve_list

        return calcinfo
