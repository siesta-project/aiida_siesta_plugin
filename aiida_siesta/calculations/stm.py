# -*- coding: utf-8 -*-
from __future__ import absolute_import
import os
import six
from aiida import orm
from aiida.common import CalcInfo, CodeInfo, InputValidationError
from aiida.common.constants import elements
from aiida.engine import CalcJob, exceptions
from aiida.orm import Dict, RemoteData, ArrayData

# Module with fdf-aware dictionary
from .tkdict import FDFDict

# See the LICENSE.txt and AUTHORS.txt files.


class STMCalculation(CalcJob):
    """
    Plugin for the "plstm" program in the Siesta distribution, which
    takes and .LDOS or .RHO file and generates a plot file to simulate
    an STM image.
    """
    _stm_plugin_version = '1.0.1'

    # Keywords that cannot be set
    # We need to canonicalize this!
    _aiida_blocked_keywords = ['mode', 'system-label', 'extension']
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
        spec.input('settings',
                   valid_type=orm.Dict,
                   help='Input settings',
                   required=False)
        spec.input('parameters', valid_type=orm.Dict, help='Input parameters')
        spec.input('ldos_folder',
                   valid_type=orm.RemoteData,
                   required=True,
                   help='Parent folder')

        # Metadata, defined in CalCJob, here we modify default
        spec.inputs['metadata']['options']['input_filename'].default=cls._DEFAULT_INPUT_FILE
        spec.inputs['metadata']['options']['output_filename'].default=cls._DEFAULT_OUTPUT_FILE
        spec.inputs['metadata']['options']['parser_name'].default='siesta.stm'

        #Outputs
        spec.output('stm_array',
                    valid_type=ArrayData,
                    required=True,
                    help='The contour data for the image')
        spec.output('output_parameters',
                    valid_type=Dict,
                    required=True,
                    help='Other output')

        # exit codes
        spec.exit_code(
            100,
            'ERROR_NO_RETRIEVED_FOLDER',
            message='The retrieved folder data node could not be accessed.')
        spec.exit_code(
            101,
            'ERROR_OUTPUT_PLOT_MISSING',
            message='The retrieved folder does not contain a CH.STM file')
        spec.exit_code(
            102,
            'ERROR_OUTPUT_PLOT_READ',
            message='The CH.STM file can not be read')


    def prepare_for_submission(self, tempfolder):
        """
        Create the input files from the input nodes passed to this instance of the `CalcJob`.

        :param tempfolder: an `aiida.common.folders.Folder` to temporarily write files on disk
        :return: `aiida.common.datastructures.CalcInfo` instance
        """

        code = self.inputs.code
        parameters = self.inputs.parameters
        ldos_folder = self.inputs.ldos_folder
        
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


        #useless so far, input_params is not stored or used!
        #input_params.update({'system-label': self._PREFIX})
        input_params.update({'mode': 'constant-height'})
        input_params.update({'extension': 'ldos'})

        # Maybe check that the 'z' coordinate makes sense...


        # ============== Creation of input file ===============

        # To have easy access to inputs metadata options
        metadataoption = self.inputs.metadata.options

        # input_filename access
        input_filename = tempfolder.get_abs_path(metadataoption.input_filename)

        # Getting the prefix from ldos_folder
        try:
            prefix = str(ldos_folder.creator.get_attribute("prefix"))
        except:
            self.report("No prefix detected from the remote folder, set 'aiida' as prefix")
            prefix = "aiida"
        ldosfile = prefix + ".LDOS"

        # Convert height to bohr...
        zbohr = input_params['z'] / 0.529177
        #This can be improved a lot.
        with open(input_filename, 'w') as infile:
            infile.write("{}\n".format(prefix))
            infile.write("ldos\n")
            infile.write("constant-height\n")
            infile.write("{}\n".format(zbohr))
            infile.write("unformatted\n")

        # ====================== Code and Calc info ========================
        # Code information object and Calc information object are now
        # only used to set up the CMDLINE (the bash line that launches siesta)
        # and to set up the list of files to retrieve.

        # The presence of a 'ldos_folder' is mandatory, to get the LDOS file
        # as indicated in the self._restart_copy_from attribute.
        # (this is not technically a restart, though)

        # It will be copied to the current calculation's working folder.
        remote_copy_list.append((ldos_folder.computer.uuid, 
            os.path.join(ldos_folder.get_remote_path(), self._restart_copy_from), 
            self._restart_copy_to))

        # Empty command line by default
        # Why use 'pop' ?
        cmdline_params = settings_dict.pop('CMDLINE', [])
        #
        # Code information object
        #
        codeinfo = CodeInfo()
        codeinfo.cmdline_params = (list(cmdline_params) + ['-z', str(zbohr), ldosfile])
        codeinfo.stdin_name = metadataoption.input_filename
        codeinfo.stdout_name = metadataoption.output_filename
        codeinfo.code_uuid = code.uuid

        #
        # Calc information object
        #
        calcinfo = CalcInfo()
        calcinfo.uuid = str(self.uuid)
        if cmdline_params:
            calcinfo.cmdline_params = list(cmdline_params)
        calcinfo.local_copy_list = local_copy_list
        calcinfo.remote_copy_list = remote_copy_list
        calcinfo.stdin_name = metadataoption.input_filename
        calcinfo.stdout_name = metadataoption.output_filename
        calcinfo.codes_info = [codeinfo]
        #

        # Retrieve by default: the output file and the plot file
        calcinfo.retrieve_list = []
        calcinfo.retrieve_list.append(metadataoption.output_filename)
        # Some logic to understand which is the plot file will be in
        # parser, here we put to retrieve every file ending in *CH.STM
        calcinfo.retrieve_list.append("*.CH.STM")
        print(calcinfo.retrieve_list)

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
