# -*- coding: utf-8 -*-
from aiida.parsers.parser import Parser
from aiida_siesta.calculations.siesta import SiestaCalculation
from aiida_siesta.calculations.stm import STMCalculation
from aiida.orm.data.parameter import ParameterData

__copyright__ = u"Copyright (c), 2015, ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE (Theory and Simulation of Materials (THEOS) and National Centre for Computational Design and Discovery of Novel Materials (NCCR MARVEL)), Switzerland and ROBERT BOSCH LLC, USA. All rights reserved."
__license__ = "MIT license, see LICENSE.txt file"
__version__ = "0.7.0"
__contributors__ = "Andrius Merkys, Giovanni Pizzi, Victor Garcia-Suarez, Alberto Garcia, Emanuele Bosoni"

# -*- coding: utf-8 -*-

from aiida.parsers.exceptions import OutputParsingError

class STMOutputParsingError(OutputParsingError):
     pass
#---------------------------

class STMParser(Parser):
    """
    Parser for the output of the "plstm" program in the Siesta distribution.
    """
    def __init__(self,calc):
        """
        Initialize the instance of STMParser
        """
        # check for valid input
        self._check_calc_compatibility(calc)
        super(STMParser, self).__init__(calc)

    def _check_calc_compatibility(self,calc):
        if not isinstance(calc,STMCalculation):
            raise STMOutputParsingError("Input calc must be a STMCalculation")

    def _get_output_nodes(self, output_path, plot_path):
        """
        Extracts output nodes from the standard output and standard error
        files. (and plot file)
        """
        from aiida.orm.data.array.trajectory import TrajectoryData
        import re

        parser_version = 'aiida-0.9.0--stm-0.1'
        parser_info = {}
        parser_info['parser_info'] = 'AiiDA STM(Siesta) Parser V. {}'.format(parser_version)
        parser_info['parser_warnings'] = []

        result_dict = {}
        result_list = []

        #in_struc = self._calc.get_inputs_dict()['structure']
        #
        # Settings (optional)
        #
        try:
             in_settings = self._calc.get_inputs_dict()['settings']
        except KeyError:
             in_settings = None

        # Add parser info dictionary
        parsed_dict = dict(result_dict.items() + parser_info.items())

        output_data = ParameterData(dict=parsed_dict)
        
        link_name = self.get_linkname_outparams()
        result_list.append((link_name,output_data))

        # Save forces and stress in an ArrayData object
        stm_data = get_stm_data(plot_path)

        if stm_data is not None:
             from aiida.orm.data.array import ArrayData
             import numpy
             arraydata = ArrayData()
             arraydata.set_array('stm_data', numpy.array(stm_data))
             result_list.append((self.get_linkname_outarray(),arraydata))

        return successful, result_list

    def parse_with_retrieved(self,retrieved):
        """
        Receives in input a dictionary of retrieved nodes.
        Does all the logic here.
        """
        
        from aiida.common.exceptions import InvalidOperation
        import os

        output_path = None
        plot_path = None
        try:
            output_path, plot_path = self._fetch_output_files(retrieved)
        except InvalidOperation:
            raise
        except IOError as e:
            self.logger.error(e.message)
            return False, ()

        if output_path is None and plot_path is None:
            self.logger.error("No output files found")
            return False, ()

        successful, out_nodes = self._get_output_nodes(output_path, plot_path)
        
        return successful, out_nodes

    def _fetch_output_files(self, retrieved):
        """
        Checks the output folder for standard output and standard error
        files, returns their absolute paths on success.

        :param retrieved: A dictionary of retrieved nodes, as obtained from the
          parser.
        """
        from aiida.common.datastructures import calc_states
        from aiida.common.exceptions import InvalidOperation
        import os

        # check in order not to overwrite anything
#         state = self._calc.get_state()
#         if state != calc_states.PARSING:
#             raise InvalidOperation("Calculation not in {} state"
#                                    .format(calc_states.PARSING) )

        # Check that the retrieved folder is there
        try:
            out_folder = retrieved[self._calc._get_linkname_retrieved()]
        except KeyError:
            raise IOError("No retrieved folder found")

        list_of_files = out_folder.get_folder_list()

        output_path = None
	plot_path = None

        if self._calc._DEFAULT_OUTPUT_FILE in list_of_files:
            output_path = os.path.join( out_folder.get_abs_path('.'),
                                        self._calc._DEFAULT_OUTPUT_FILE )
        if self._calc._DEFAULT_PLOT_FILE in list_of_files:
            plot_path  = os.path.join( out_folder.get_abs_path('.'),
                                        self._calc._DEFAULT_PLOT_FILE )

        return output_path, plot_path

    def get_linkname_outarray(self):
        """                                                                     
        Returns the name of the link to the output_array                        
        In Siesta, Node exists to hold the final forces and stress,
        pending the implementation of trajectory data.
        """
        return 'stm_array'

   def get_stm_data(self,plot_path):
        """
        Parses the STM plot file to get a list of lists
        """
        import numpy as np
        from aiida.common.exceptions import InputValidationError
        from aiida.common.exceptions import ValidationError

        #f=open(plot_path)
        dummy=[ 3.14, 2.1718, 0.0, 1.0 ]
        return dummy
        
