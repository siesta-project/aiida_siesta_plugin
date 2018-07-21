# -*- coding: utf-8 -*-
import numpy as np
from aiida.orm.data.parameter import ParameterData
from aiida.parsers.parser import Parser
from aiida.parsers.exceptions import OutputParsingError
from aiida_siesta.calculations.vibra import VibraCalculation

__copyright__ = u"Copyright (c), 2015, ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE (Theory and Simulation of Materials (THEOS) and National Centre for Computational Design and Discovery of Novel Materials (NCCR MARVEL)), Switzerland and ROBERT BOSCH LLC, USA. All rights reserved."
__license__ = "MIT license, see LICENSE.txt file"
__version__ = "0.12.0"
__contributors__ = "Victor M. Garcia-Suarez"
# Based on the 0.9.0 version of the STM workflow developed by Alberto
# Garcia for the aiida_siesta plugin

class VibraOutputParsingError(OutputParsingError):
     pass

class VibraParser(Parser):
    """
    Parser for the output of a Vibra calculation.
    """
    def __init__(self,calc):
        """
        Initialize the instance of VibraParser
        """
        # check for valid input
        self._check_calc_compatibility(calc)
        super(VibraParser, self).__init__(calc)

    def _check_calc_compatibility(self,calc):
        if not isinstance(calc,VibraCalculation):
            raise VibraOutputParsingError("Input calc must be a VibraCalculation")

    def _get_output_nodes(self, output_path, bands_path):
        """
        Extracts output nodes from the standard output and standard error
        files. (And XML and JSON files)
        """
        from aiida.orm.data.array.trajectory import TrajectoryData
        import re

        result_list = []

        # Add errors
        successful = True
        if output_path is None:
            errors_list = ['WARNING: No aiida.out file...']
        else:
            successful, errors_list = self.get_errors_from_file(output_path)

        result_dict = {}
        result_dict["errors"] = errors_list

        # Add warnings
        warnings_list = self.get_warnings_from_file(output_path)
        result_dict["warnings"] = warnings_list

        # Add outuput data
        output_dict = self.get_output_from_file(output_path)
        result_dict.update(output_dict)

        # Add parser info dictionary
        parser_info = {}
        parser_version = 'aiida-0.11.0--plugin-0.11.5'
        parser_info['parser_info'] =\
            'AiiDA Vibra Parser V. {}'.format(parser_version)
        parser_info['parser_warnings'] = []
        parsed_dict = dict(result_dict.items() + parser_info.items())

        output_data = ParameterData(dict=parsed_dict)
        
        link_name = self.get_linkname_outparams()
        result_list.append((link_name,output_data))

        # Parse band-structure information if available
        if bands_path is not None:
             bands, coords = self.get_bands(bands_path)
             from aiida.orm.data.array.bands import BandsData
             arraybands = BandsData()
             arraybands.set_kpoints(self._calc.inp.bandskpoints.get_kpoints(cartesian=True))
             arraybands.set_bands(bands,units="eV")
             result_list.append((self.get_linkname_bandsarray(), arraybands))
             bandsparameters = ParameterData(dict={"kp_coordinates": coords})
             result_list.append((self.get_linkname_bandsparameters(), bandsparameters))

        return successful, result_list

    def parse_with_retrieved(self,retrieved):
        """
        Receives in input a dictionary of retrieved nodes.
        Does all the logic here.
        """
        
        from aiida.common.exceptions import InvalidOperation
        import os

        output_path = None
        bands_path = None
        try:
            output_path, bands_path=\
                self._fetch_output_files(retrieved)
        except InvalidOperation:
            raise
        except IOError as e:
            self.logger.error(e.message)
            return False, ()

        if output_path is None and bands_path is None:
            self.logger.error("No output files found")
            return False, ()

        successful, out_nodes = self._get_output_nodes(output_path, bands_path)
        
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

        # Check that the retrieved folder is there
        try:
            out_folder = retrieved[self._calc._get_linkname_retrieved()]
        except KeyError:
            raise IOError("No retrieved folder found")

        list_of_files = out_folder.get_folder_list()

        output_path = None
        bands_path = None

        if self._calc._DEFAULT_OUTPUT_FILE in list_of_files:
            output_path = os.path.join( out_folder.get_abs_path('.'),
                                        self._calc._DEFAULT_OUTPUT_FILE )
        if self._calc._DEFAULT_BANDS_FILE in list_of_files:
            bands_path = os.path.join( out_folder.get_abs_path('.'),
                                        self._calc._DEFAULT_BANDS_FILE )

        return output_path, bands_path

    def get_errors_from_file(self,output_path):
        """
        Generates a list of errors from the 'aiida.out' file.

        :param output_path: 

        Returns a boolean indicating success (True) or failure (False)
        and a list of strings.
        """
        f = open(output_path)
        lines = f.read().split('\n')   # There will be a final '' element

        import re
     
        # Search for 'Error' messages, log them, and return immediately
        lineerror = []
        there_are_fatals = False 
        for line in lines:
            if re.match('^.*Error.*$',line):
                self.logger.error(line)
                lineerror.append(line)
                there_are_fatals = True
               
        if there_are_fatals:
            lineeror.append(lines[-1])
            return False, lineerror

        # Make sure that the job did finish (and was not interrupted
        # externally)

        normal_end = False
        for line in lines:
            if re.match('^.*Zero point energy.*$',line):
                normal_end = True
               
        if normal_end == False:
            lines[-1] = 'FATAL: ABNORMAL_EXTERNAL_TERMINATION'
            self.logger.error("Calculation interrupted externally")
            return False, lines[-2:] # Return also last line of the file

        return True, lineerror

    def get_warnings_from_file(self,output_path):
        """
        Generates a list of warnings from the 'aiida.out' file.

        :param output_path: 

        Returns a list of strings.
        """
        f = open(output_path)
        lines = f.read().split('\n')   # There will be a final '' element

        import re

        # Find warnings
        linewarning = []
        for line in lines:
            if re.match('^.*Warning.*$',line):
                linewarning.append(line)
     
        return linewarning

    def get_output_from_file(self,output_path):
        """
        Generates a list of variables from the 'aiida.out' file.

        :param output_path: 

        Returns a list of strings.
        """
        f = open(output_path)
        lines = f.read().split('\n')   # There will be a final '' element

        import re

        # Find data
        output_dict = {}
        for line in lines:
            if re.match('^.*System Name.*$',line):
                output_dict['system_name'] = line.split()[-1]
            if re.match('^.*System Label.*$',line):
                output_dict['system_label'] = line.split()[-1]
            if re.match('^.*Number of unit cells in Supercell.*$',line):
                output_dict['number_of_unit_cells'] = int(line.split()[-1])
            if re.match('^.*Eigenvectors =.*$',line):
                output_dict['eigenvectors_calc'] = line.split()[-1]
            if re.match('^.*Zero point energy.*$',line):
                output_dict['zero_point_energy'] = float(line.split()[-2])
     
        return output_dict

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
        if self._calc.inp.bandskpoints.labels == None:
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
                    spinup[i, j] = (
                        float(tottx[i * (block_length + 3) + 6 + j + 3]))
                    if (nspins == 2):
                        # Probably wrong! - need to test!!!
                        spindown[i, j] = (float(
                            tottx[i * (nbands * 2 + 3) + 6 + j + 3 + nbands]))
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
                    spinup[i, j] = (
                        float(tottx[i * (block_length + 1) + 8 + j + 1]))
                    if (nspins == 2):
                        # Probably wrong! - need to test!!!
                        spindown[i, j] = (float(
                            tottx[i * (nbands * 2 + 1) + 8 + j + 1 + nbands]))
        if (nspins == 2):
            bands = (spinup, spindown)
        elif (nspins == 1):
            bands = spinup
        else:
            raise NotImplementedError(
                'nspins=4: non collinear bands not implemented yet')

        return (bands, coords)
