# -*- coding: utf-8 -*-
import numpy as np
from aiida.parsers.parser import Parser
from aiida_siesta.calculations.siesta import SiestaCalculation
from aiida.orm.data.parameter import ParameterData

# TODO Get modules metadata from setup script.
__copyright__ = u"Copyright (c), 2015, ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE (Theory and Simulation of Materials (THEOS) and National Centre for Computational Design and Discovery of Novel Materials (NCCR MARVEL)), Switzerland and ROBERT BOSCH LLC, USA. All rights reserved."
__license__ = "MIT license, see LICENSE.txt file"
__version__ = "0.9.10"
__contributors__ = "Andrius Merkys, Giovanni Pizzi, Victor Garcia-Suarez, Alberto Garcia, Emanuele Bosoni"


# These auxiliary functions should be put in another module...
#
# List of scalar values from CML to be transferred to AiiDA
#
standard_output_list = [ 'siesta:FreeE', 'siesta:E_KS',
                         'siesta:Ebs', 'siesta:E_Fermi',
                         'siesta:stot']  ## leave svec for later


def text_to_array(s, dtype):
    return np.array(s.replace("\n", "").split(), dtype=dtype)


def get_parsed_xml_doc(xml_path):

     from xml.dom import minidom

     try:
          xmldoc = minidom.parse(xml_path)
     except:
          xmldoc = None

     # We might want to add some extra consistency checks
     
     return xmldoc

def get_dict_from_xml_doc(xmldoc):

     # Scalar items
     
     scalar_dict = {}
     
     # Metadata items
     itemlist = xmldoc.getElementsByTagName('metadata')
     for s in itemlist:
         #
         # Maybe make sure that 'name' does not contain
         # forbidden characters
         #
         name = s.attributes['name'].value
         value = s.attributes['content'].value
         scalar_dict[name] = value

     # Scalar output items
     # From the last "SCF Finalization" module 
     # This means that we do not record non-converged values
     
     itemlist = xmldoc.getElementsByTagName('module')

     scf_final = None
     for m in itemlist:
       if 'title' in m.attributes.keys():
          # Get last scf finalization module
          if m.attributes['title'].value == "SCF Finalization":
               scf_final = m

     if scf_final is not None:

      # wrapped in <property> elements with a <scalar> child
      props = scf_final.getElementsByTagName('property')
     
      for s in props:
        if 'dictRef' in s.attributes.keys():
          name = s.attributes['dictRef'].value
          if name in standard_output_list:
             data = s.getElementsByTagName('scalar')[0]
             value = data.childNodes[0].nodeValue
             units = data.attributes['units'].value
             loc_colon = units.find(':')
             unit_name = units[loc_colon+1:]
             loc_colon = name.find(':')
             reduced_name = name[loc_colon+1:]

             # Put units in separate entries, as in QE
             # Use numbers (floats) instead of strings
             
             scalar_dict[reduced_name] = float(value)
             scalar_dict[reduced_name+"_units"] = unit_name

     scalar_dict['variable_geometry'] = is_variable_geometry(xmldoc)
     #
     # Sizes of orbital set (and non-zero interactions), and mesh 
     #
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
     for m in itemlist:
       # Check the type of the first "step" module, which is a "geometry" one
       if 'serial' in m.attributes.keys():
           if 'dictRef' in m.attributes.keys():
               if m.attributes['dictRef'].value == "Single-Point":
                   return False
               else:
                   return True

     # If we reach this point, something is very wrong
     return False
     
def get_sizes_info(xmldoc):
     """
     Gets the number of orbitals and non-zero interactions
     """
     no_u = None
     nnz = None
     mesh = None
     
     itemlist = xmldoc.getElementsByTagName('module')
     for m in itemlist:
       # Process the first "step" module, which is a "geometry" one
       if 'serial' in m.attributes.keys():
           # Get properties here
           props_list = m.getElementsByTagName('property')
           for p in props_list:
               if 'dictRef' in p.attributes.keys():
                   if p.attributes['dictRef'].value == "siesta:no_u":
                       scalar = p.getElementsByTagName('scalar')[0]
                       no_u = int(scalar.childNodes[0].data)
                   if p.attributes['dictRef'].value == "siesta:nnz":
                       scalar = p.getElementsByTagName('scalar')[0]
                       nnz = int(scalar.childNodes[0].data)
                   if p.attributes['dictRef'].value == "siesta:ntm":
                       array = p.getElementsByTagName('array')[0]
                       mesh = [int(s) for s in array.childNodes[0].data.split()]
               

           return no_u, nnz, mesh

def get_last_structure(xmldoc, input_structure):

    from aiida.orm import DataFactory

    itemlist = xmldoc.getElementsByTagName('module')
    #
    # Use the last "geometry" module, and not the
    # "Finalization" one.

    finalmodule = None
    for m in itemlist:
      # Get a "geometry" module by the criteria:
      if 'serial' in m.attributes.keys():
          if 'dictRef' in m.attributes.keys():
              if m.attributes['dictRef'].value != "SCF":
                  finalmodule = m

    # In case there is no appropriate data, fall back and
    # at least return the initial structure
    # (this should not be necessary, as the initial Geometry module
    # is opened very soon)
    if finalmodule is None:
         #self.logger.warning("Returning input structure in output_structure node")
         return input_structure
    
    atoms = finalmodule.getElementsByTagName('atom')
    cellvectors = finalmodule.getElementsByTagName('latticeVector')

    atomlist = []

    for a in atoms:
         kind = a.attributes['elementType'].value
         x = a.attributes['x3'].value
         y = a.attributes['y3'].value
         z = a.attributes['z3'].value
         atomlist.append([kind,[float(x),float(y),float(z)]])
    
    cell = []
    for l in cellvectors:
         data = l.childNodes[0].data.split()
         cell.append([float(s) for s in data])

    # Generally it is better to pass the input structure
    # and reset the data, since site 'names' are not handled by
    # the CML file (at least not in Siesta versions <= 4.0)
    #

    s = input_structure.copy()
    s.reset_cell(cell)
    new_pos = [atom[1] for atom in atomlist]
    s.reset_sites_positions(new_pos)
          
    return s

                        
def get_final_forces_and_stress(xmldoc):
 #
 # Extracts final forces and stress as lists of lists...
 #
 itemlist = xmldoc.getElementsByTagName('module')

 # Note: In modern versions of Siesta, forces and stresses
 # are written in the "SCF Finalization" modules at the end
 # of each geometry step.

 # Search for the last one of those modules
 
 scf_final = None
 for m in itemlist:
     if 'title' in m.attributes.keys():
          # Get last scf finalization module
          if m.attributes['title'].value == "SCF Finalization":
               scf_final = m

 forces = None
 stress = None

 if scf_final is not None:
      props = scf_final.getElementsByTagName('property')
      for p in props:
        if 'dictRef' in p.attributes.keys():

           if p.attributes['dictRef'].value=='siesta:forces':
                mat = p.getElementsByTagName('matrix')[0]
                # Get flat list and reshape as list of lists
                # using info on rows and columns in CML file
                rows = int(mat.attributes['rows'].value)
                cols = int(mat.attributes['columns'].value)
                f = mat.childNodes[0].data.split()
                f = [float(x) for x in f]
                forces = [ f[rows*i : rows*(i+1)] for i in range(cols)]

           if p.attributes['dictRef'].value=='siesta:stress':
                mat = p.getElementsByTagName('matrix')[0]
                # Get flat list and reshape as list of lists
                # using info on rows and columns in CML file
                rows = int(mat.attributes['rows'].value)
                cols = int(mat.attributes['columns'].value)
                s = mat.childNodes[0].data.split()
                s = [float(x) for x in s]
                stress = [ s[rows*i : rows*(i+1)] for i in range(cols)]

 return forces, stress


#----------------------------------------------------------------------

from aiida.parsers.exceptions import OutputParsingError

class SiestaOutputParsingError(OutputParsingError):
     pass
class SiestaCMLParsingError(OutputParsingError):
     pass

#---------------------------

class SiestaParser(Parser):
    """
    Parser for the output of Siesta.
    """
    def __init__(self,calc):
        """
        Initialize the instance of SiestaParser
        """
        # check for valid input
        self._check_calc_compatibility(calc)
        super(SiestaParser, self).__init__(calc)

    def _check_calc_compatibility(self,calc):
        if not isinstance(calc,SiestaCalculation):
            raise SiestaOutputParsingError("Input calc must be a SiestaCalculation")

    def _get_output_nodes(self, output_path, messages_path, xml_path, json_path, bands_path):
        """
        Extracts output nodes from the standard output and standard error
        files. (And XML and JSON files)
        """
        from aiida.orm.data.array.trajectory import TrajectoryData
        import re

        parser_version = 'aiida-0.12.0--plugin-0.9.10'
        parser_info = {}
        parser_info['parser_info'] = 'AiiDA Siesta Parser V. {}'.format(parser_version)
        parser_info['parser_warnings'] = []


        result_list = []

        if xml_path is None:
            self.logger.error("Could not find a CML file to parse")
            # NOTE aiida.xml is not there?
            raise SiestaOutputParsingError("Could not find a CML file to parse")

        # We get everything from the CML file

        xmldoc = get_parsed_xml_doc(xml_path)
        if xmldoc is None:
            self.logger.error("Malformed CML file: cannot parse")
            raise SiestaCMLParsingError("Malformed CML file: cannot parse")
        
        # These are examples of how we can access input items
        #
        # Structure (mandatory)
        #
        in_struc = self._calc.get_inputs_dict()['structure']
        #
        # Settings (optional)
        #
        try:
             in_settings = self._calc.get_inputs_dict()['settings']
        except KeyError:
             in_settings = None

        result_dict = get_dict_from_xml_doc(xmldoc)

        # Add timing information

        if json_path is None:
            self.logger.info("Could not find a time.json file to parse")
        else:
             from json_time import get_timing_info
             global_time, timing_decomp = get_timing_info(json_path)
             if global_time is None:
                  self.logger.info("Cannot fully parse the time.json file")
             else:
                  result_dict["global_time"] = global_time
                  result_dict["timing_decomposition"] = timing_decomp
        
        # Add warnings
        successful = True
        if messages_path is None:
             # Perhaps using an old version of Siesta
             warnings_list = ['WARNING: No MESSAGES file...']
        else:
             successful, warnings_list = self.get_warnings_from_file(messages_path)

        result_dict["warnings"] = warnings_list
        
        # Add parser info dictionary
        parsed_dict = dict(result_dict.items() + parser_info.items())

        output_data = ParameterData(dict=parsed_dict)
        
        link_name = self.get_linkname_outparams()
        result_list.append((link_name,output_data))

        # If the structure has changed, save it
        if is_variable_geometry(xmldoc):
             # Get the input structure to copy its site names,
             # as the CML file traditionally contained only the
             # atomic symbols.
             #
             struc = get_last_structure(xmldoc,in_struc)
             result_list.append((self.get_linkname_outstructure(),struc))

        # Save forces and stress in an ArrayData object
        forces, stress = get_final_forces_and_stress(xmldoc)

        if forces is not None and stress is not None:
             from aiida.orm.data.array import ArrayData
             arraydata = ArrayData()
             arraydata.set_array('forces', np.array(forces))
             arraydata.set_array('stress', np.array(stress))
             result_list.append((self.get_linkname_outarray(),arraydata))

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
        messages_path  = None
        xml_path  = None
        json_path  = None
        bands_path = None
        try:
            output_path, messages_path, xml_path, json_path, bands_path = \
                             self._fetch_output_files(retrieved)
        except InvalidOperation:
            raise
        except IOError as e:
            self.logger.error(e.message)
            return False, ()

        if output_path is None and messages_path is None and xml_path is None:
            self.logger.error("No output files found")
            return False, ()

        successful, out_nodes = self._get_output_nodes(output_path, messages_path,
                                                       xml_path, json_path, bands_path)
        
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
        messages_path  = None
        xml_path  = None
        json_path  = None
        bands_path = None

        if self._calc._DEFAULT_OUTPUT_FILE in list_of_files:
            output_path = os.path.join( out_folder.get_abs_path('.'),
                                        self._calc._DEFAULT_OUTPUT_FILE )
        if self._calc._DEFAULT_XML_FILE in list_of_files:
            xml_path = os.path.join( out_folder.get_abs_path('.'),
                                        self._calc._DEFAULT_XML_FILE )
        if self._calc._DEFAULT_JSON_FILE in list_of_files:
            json_path = os.path.join( out_folder.get_abs_path('.'),
                                        self._calc._DEFAULT_JSON_FILE )
        if self._calc._DEFAULT_MESSAGES_FILE in list_of_files:
            messages_path  = os.path.join( out_folder.get_abs_path('.'),
                                        self._calc._DEFAULT_MESSAGES_FILE )
        if self._calc._DEFAULT_BANDS_FILE in list_of_files:
            bands_path  = os.path.join( out_folder.get_abs_path('.'),
                                        self._calc._DEFAULT_BANDS_FILE )

        return output_path, messages_path, xml_path, json_path, bands_path

    def get_warnings_from_file(self,messages_path):
     """
     Generates a list of warnings from the 'MESSAGES' file, which
     contains a line per message, prefixed with 'INFO',
     'WARNING' or 'FATAL'.

     :param messages_path: 

     Returns a boolean indicating success (True) or failure (False)
     and a list of strings.
     """
     f=open(messages_path)
     lines=f.read().split('\n')   # There will be a final '' element

     import re
     
     # Search for 'FATAL:' messages, log them, and return immediately
     there_are_fatals = False 
     for line in lines:
          if re.match('^FATAL:.*$',line):
               self.logger.error(line)
               there_are_fatals = True
               
     if there_are_fatals:
          return False, lines[:-1]  # Remove last (empty) element

     # Make sure that the job did finish (and was not interrupted
     # externally)

     normal_end = False
     for line in lines:
          if re.match('^INFO: Job completed.*$',line):
               normal_end = True
               
     if normal_end == False:
          lines[-1] = 'FATAL: ABNORMAL_EXTERNAL_TERMINATION'
          self.logger.error("Calculation interrupted externally")
          return False, lines

     # (Insert any other "non-success" conditions before next section)
     # (e.g.: be very picky about (some) 'WARNING:' messages)
     
     # Return with success flag
     
     return True, lines[:-1]  # Remove last (empty) element

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


    def get_linkname_outstructure(self):
        """
        Returns the name of the link to the output_structure
        Node exists if positions or cell changed.
        """
        return 'output_structure'

    def get_linkname_outarray(self):
        """                                                                     
        Returns the name of the link to the output_array                        
        In Siesta, Node exists to hold the final forces and stress,
        pending the implementation of trajectory data.
        """
        return 'output_array'

    def get_linkname_bandsarray(self):
        """                                                                     
        Returns the name of the link to the bands_array                        
        In Siesta, Node exists to hold the bands,
        pending the implementation of trajectory data.
        """
        return 'bands_array'

    def get_linkname_bandsparameters(self):
        """
        Returns the name of the link to the bands_path.
        X-axis data for bands. Maybe should use ArrayData (db-integrity?).
        """
        return 'bands_parameters'
