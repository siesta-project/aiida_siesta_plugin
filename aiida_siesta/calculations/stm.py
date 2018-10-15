# -*- coding: utf-8 -*-
import os

# Module with fdf-aware dictionary
from tkdict import FDFDict

from aiida.orm.calculation.job import JobCalculation
from aiida.common.exceptions import InputValidationError
from aiida.common.datastructures import CalcInfo
from aiida.common.utils import classproperty
from aiida.common.datastructures import CodeInfo

from aiida.orm.data.parameter import ParameterData
from aiida.orm.data.remote import RemoteData 

__copyright__ = u"Copyright (c), 2015, ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE (Theory and Simulation of Materials (THEOS) and National Centre for Computational Design and Discovery of Novel Materials (NCCR MARVEL)), Switzerland and ROBERT BOSCH LLC, USA. All rights reserved."
__license__ = "MIT license, see LICENSE.txt file"
__version__ = "0.9.10"
__contributors__ = "Victor M. Garcia-Suarez, Alberto Garcia"

class STMCalculation(JobCalculation):
    """
    Plugin for the "plstm" program in the Siesta distribution, which
    takes and .LDOS or .RHO file and generates a plot file to simulate
    an STM image.
    """
    _stm_plugin_version = 'aiida-0.12.0--stm-0.9.10'
    
    def _init_internal_params(self):
        super(STMCalculation, self)._init_internal_params()

        # Default Siesta output parser provided by AiiDA
        self._default_parser = "siesta.stm"

        # Keywords that cannot be set
        # We need to canonicalize this!

        self._aiida_blocked_keywords = ['mode','system-label','extension']

        # Default input and output files                                        
        self._DEFAULT_INPUT_FILE = 'stm.in'
        self._DEFAULT_OUTPUT_FILE = 'stm.out'
	self._DEFAULT_PLOT_FILE = 'aiida.CH.STM'

        self._OUTPUT_SUBFOLDER = './'
        self._PREFIX = 'aiida'
        self._INPUT_FILE_NAME = 'stm.in'
        self._OUTPUT_FILE_NAME = 'stm.out'
	self._PLOT_FILE_NAME = 'aiida.CH.STM'

        # in restarts, it will copy from the parent the following
        self._restart_copy_from = os.path.join(self._OUTPUT_SUBFOLDER, '*.LDOS')

        # in restarts, it will copy the previous folder in the following one
        self._restart_copy_to = self._OUTPUT_SUBFOLDER


    @classproperty
    def _use_methods(cls):
        """
        Extend the parent _use_methods with further keys.
        """
        retdict = JobCalculation._use_methods

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

        return retdict

    def _prepare_for_submission(self,tempfolder,
                                    inputdict):        
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
        settings = inputdict.pop(self.get_linkname('settings'),None)
        if settings is None:
            settings_dict = {}
        else:
            if not isinstance(settings,  ParameterData):
                raise InputValidationError("settings, if specified, must be of "
                                           "type ParameterData")
            
            # Settings converted to UPPERCASE
            # Presumably to standardize the usage and avoid
            # ambiguities
            settings_dict = _uppercase_dict(settings.get_dict(),
                                            dict_name='settings')

        try:
            parameters = inputdict.pop(self.get_linkname('parameters'))
        except KeyError:
            raise InputValidationError("No parameters specified for this "
                "calculation")
        if not isinstance(parameters, ParameterData):
            raise InputValidationError("parameters is not of type "
                "ParameterData")

        try:
            parent_calc_folder = inputdict.pop(self.get_linkname('parent_folder'))
        except KeyError:
            raise InputValidationError("No parent_calc_folder specified for this "
                "calculation")
        if not isinstance(parent_calc_folder,  RemoteData):
            raise InputValidationError("parent_calc_folder, if specified,"
                                       "must be of type RemoteData")

        #
        # Important note: This program should NOT be run with MPI.
        # Scripts using this plugin should use:
        #
        # calc.set_withmpi(False)
        #
        # We do it right here, and hope that it will not be overriden
        #
        # self.set_withmpi(False)
        #
        try:
            code = inputdict.pop(self.get_linkname('code'))
        except KeyError:
            raise InputValidationError("No code specified for this calculation")
                                
        # Here, there should be no more parameters...
        if inputdict:
            raise InputValidationError("The following input data nodes are "
                "unrecognized: {}".format(inputdict.keys()))

        # END OF INITIAL INPUT CHECK #

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
                        "input parameters".format(input_params.get_last_key(key)))

        input_params.update({'system-label': self._PREFIX})
        input_params.update({'mode': 'constant-height'})
        input_params.update({'extension': 'ldos'})

        # Maybe check that the 'z' coordinate makes sense...
        
        input_filename = tempfolder.get_abs_path(self._INPUT_FILE_NAME)

        with open(input_filename,'w') as infile:
            infile.write("aiida\n")
            infile.write("ldos\n")
            infile.write("constant-height\n")
            # Convert height to bohr...
            infile.write("{}\n".format(input_params['z']/0.529177))
            infile.write("unformatted\n")

        # ------------------------------------- END of input file creation
        
        # The presence of a 'parent_calc_folder' is mandatory, to get the LDOS file
        # as indicated in the self._restart_copy_from attribute.
        # (this is not technically a restart, though)
        
        # It will be copied to the current calculation's working folder.
        
        if parent_calc_folder is not None:
            remote_copy_list.append(
                    (parent_calc_folder.get_computer().uuid,
                     os.path.join(parent_calc_folder.get_remote_path(),
                                  self._restart_copy_from),
                     self._restart_copy_to
                     ))

        calcinfo = CalcInfo()
        calcinfo.uuid = self.uuid
        #
        # Empty command line by default
        # Why use 'pop' ?
        cmdline_params = settings_dict.pop('CMDLINE', [])
        
        if cmdline_params: 
            calcinfo.cmdline_params = list(cmdline_params)
        calcinfo.local_copy_list = local_copy_list
        calcinfo.remote_copy_list = remote_copy_list
        
        calcinfo.stdin_name = self._INPUT_FILE_NAME
        calcinfo.stdout_name = self._OUTPUT_FILE_NAME
        #
        # Code information object
        #
        codeinfo = CodeInfo()
        codeinfo.cmdline_params = list(cmdline_params)
        codeinfo.stdin_name = self._INPUT_FILE_NAME
        codeinfo.stdout_name = self._OUTPUT_FILE_NAME
        codeinfo.code_uuid = code.uuid
        calcinfo.codes_info = [codeinfo]

        # Retrieve by default: the output file and the plot file
        
        calcinfo.retrieve_list = []         
        calcinfo.retrieve_list.append(self._OUTPUT_FILE_NAME)
        calcinfo.retrieve_list.append(self._PLOT_FILE_NAME)

        # Any other files specified in the settings dictionary
        settings_retrieve_list = settings_dict.pop('ADDITIONAL_RETRIEVE_LIST',[])
        calcinfo.retrieve_list += settings_retrieve_list

        return calcinfo

    def _set_parent_remotedata(self,remotedata):
        """
        Used to set a parent remotefolder that holds the LDOS file
        from a previous Siesta calculation
        """
        from aiida.common.exceptions import ValidationError
        
        if not isinstance(remotedata,RemoteData):
            raise ValueError('remotedata must be a RemoteData')
        
        # complain if another remotedata is already found
        input_remote = self.get_inputs(node_type=RemoteData)
        if input_remote:
            raise ValidationError("Cannot set several parent calculation to a "
                "{} calculation".format(self.__class__.__name__))

        self.use_parent_folder(remotedata)

def get_input_data_text(key,val, mapping=None):
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
            
            list_of_strings.append((idx,"  {0}({2}) = {1}\n".format(
                key, conv_to_fortran(itemval), idx)))
        
        # I first have to resort, then to remove the index from the first
        # column, finally to join the strings
        list_of_strings = zip(*sorted(list_of_strings))[1]
        return "".join(list_of_strings)                          
    elif hasattr(val,'__iter__'):
        # a list/array/tuple of values
        list_of_strings = [
            "{0}({2})  {1}\n".format(key, conv_to_fortran(itemval), idx+1)
            for idx, itemval in enumerate(val)]
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

    if isinstance(d,dict):
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
        raise TypeError("_lowercase_dict accepts only dictionaries as argument")

