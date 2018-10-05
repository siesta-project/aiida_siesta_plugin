#!/usr/bin/env runaiida
# -*- coding: utf-8 -*-

__copyright__ = u"Copyright (c), 2015, ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE (Theory and Simulation of Materials (THEOS) and National Centre for Computational Design and Discovery of Novel Materials (NCCR MARVEL)), Switzerland and ROBERT BOSCH LLC, USA. All rights reserved."
__license__ = "MIT license, see LICENSE.txt file"
__version__ = "0.7.0"
__contributors__ = "Andrea Cepellotti, Victor Garcia-Suarez, Alberto Garcia, Emanuele Bosoni"

import sys
import os
import pymatgen as mg
from aiida.tools import get_explicit_kpoints_path
from aiida.common.example_helpers import test_and_get_code
from aiida.common.exceptions import NotExistent


# This script will send a Siesta calculation on a structure taken from
# a cif file.
# The band structure is calculated and the kpoint path is automatically 
# generated using seekpath.

################################################################

PsfData = DataFactory('siesta.psf')
ParameterData = DataFactory('parameter')
KpointsData = DataFactory('array.kpoints')
StructureData = DataFactory('structure')

try:
    dontsend = sys.argv[1]
    if dontsend == "--dont-send":
        submit_test = True
    elif dontsend == "--send":
        submit_test = False
    else:
        raise IndexError
except IndexError:
    print >> sys.stderr, ("The first parameter can only be either "
                          "--send or --dont-send")
    sys.exit(1)

try:
    codename = sys.argv[2]
except IndexError:
    codename = 'siesta4.0.1@parsons'

code = test_and_get_code(codename, expected_code_type='siesta.siesta')
#
#  Set up calculation object first
#
calc = code.new_calc()
calc.label = "O_el_cell_spin"
calc.description = "Band calculation with Siesta. O elementary cell spin polarized"
calc.set_max_wallclock_seconds(30*60) # 30 min

calc.set_resources({"num_machines": 1})
#------------------
queue = None
# calc.set_queue_name(queue)
#------------------


print calc.description

#
#----Settings first  -----------------------------
#
settings_dict={'additional_retrieve_list': ['aiida.BONDS', 'aiida.EIG']}
settings = ParameterData(dict=settings_dict)
calc.use_settings(settings)
#---------------------------------------------------

#
# Structure -----------------------------------------
#
#Passing through SeeK-path first, to get the standardized cell.
#Necessary for the automatic choice of the bands path.
#
structure = mg.Structure.from_file("data/O2_ICSD_173933.cif", primitive=False)
s = StructureData(pymatgen_structure=structure)
elements = list(s.get_symbols_set())
seekpath_parameters = ParameterData(dict={'reference_distance': 0.02,'symprec': 0.0001})
result=get_explicit_kpoints_path(s, **seekpath_parameters.get_dict())
newstructure = result['primitive_structure']
calc.use_structure(newstructure)
#-------------------------------------------------------------

#
# Parameters ---------------------------------------------------
#
# Note the use of '.' in some entries. This will be fixed below.
# Note also that some entries have ':' as separator. This is not
# allowed in Siesta, and will be fixed by the plugin itself. The
# latter case is an unfortunate historical choice. It should not
# be used in modern scripts.
#
params_dict= {
'xc.functional': 'LDA',
'xc.authors': 'CA',
'spin:polarized': True,
'noncollinearspin': False,
'mesh-cutoff': '200.000 Ry',
'max-scfiterations': 1000,
'dm-numberpulay': 5,
'dm-mixingweight': 0.050,
'dm-tolerance': 1.e-4,
'dm-mixscf1': True,
'negl-nonoverlap-int': False,
'solution-method': 'diagon',
'electronic-temperature': '100.000 K',
'md-typeofrun': 'cg',
'md-numcgsteps': 2,
'md-maxcgdispl': '0.200 bohr',
'md-maxforcetol': '0.050 eV/Ang',
'writeforces': True,
'writecoorstep': True,
'write-mulliken-pop': 1
}
#
# Sanitize, as '.' is not kosher for the database handlers
#
params_dict = { k.replace('.','-') :v for k,v in params_dict.iteritems() }
#
parameters = ParameterData(dict=params_dict)
calc.use_parameters(parameters)
#
#----------------------------------------------------------
#
# Basis Set Info ------------------------------------------
# The basis dictionary follows the 'parameters' convention
#
basis_dict = {
'pao-basistype': 'split',
'pao-splitnorm': 0.150,
'pao-energyshift': '0.020 Ry',
'%block pao-basis-sizes' :"""
O    SZP  """,
}
#
basis_dict = { k.replace('.','-') :v for k,v in  basis_dict.iteritems() }
#
basis = ParameterData(dict=basis_dict)
calc.use_basis(basis)
#--------------------------------------------------------------


# Pseudopotentials ----------------------------------------------
#
# This exemplifies the handling of pseudos for different species
# Those sharing the same pseudo should be indicated.
# Families support is not yet available for this.
#
raw_pseudos = [ ("O.psf", 'O')]

for fname, kinds, in raw_pseudos:
    absname = os.path.realpath(os.path.join(os.path.dirname(__file__),
                                            "data",fname))
    pseudo, created = PsfData.get_or_create(absname,use_first=True)
    if created:
        print "Created the pseudo for {}".format(kinds)
    else:
        print "Using the pseudo for {} from DB: {}".format(kinds,pseudo.pk)
        
    # Attach pseudo node to the calculation
    calc.use_pseudo(pseudo,kind=kinds)
#-------------------------------------------------------------------

# K-points for scf cycle -------------------------------------------
kts = KpointsData()
kpoints_mesh = 4
kts.set_kpoints_mesh([kpoints_mesh,kpoints_mesh,kpoints_mesh])

# K-points for bands, uncomment your favourite  --------------------
# 
bandskpoints = KpointsData()

##Making use of SeeK-path for the automatic path
##The choice of the distance between kpoints is in the call seekpath_parameters
##All high symmetry points included, labels already included
bandskpoints=result['explicit_kpoints']


calc.use_kpoints(kts)
calc.use_bandskpoints(bandskpoints)

#from aiida.orm.data.remote import RemoteData
#calc.set_outdir(remotedata)

if submit_test:
    subfolder, script_filename = calc.submit_test()
    print "Test_submit for calculation (uuid='{}')".format(
        calc.uuid)
    print "Submit file in {}".format(os.path.join(
        os.path.relpath(subfolder.abspath),
        script_filename
        ))
else:
    calc.store_all()
    print "created calculation; calc=Calculation(uuid='{}') # ID={}".format(
        calc.uuid,calc.dbnode.pk)
    calc.submit()
    print "submitted calculation; calc=Calculation(uuid='{}') # ID={}".format(
        calc.uuid,calc.dbnode.pk)

