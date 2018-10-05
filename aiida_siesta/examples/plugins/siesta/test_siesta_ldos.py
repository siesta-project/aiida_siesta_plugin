#!/usr/bin/env runaiida
# -*- coding: utf-8 -*-

__copyright__ = u"Copyright (c), 2015, ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE (Theory and Simulation of Materials (THEOS) and National Centre for Computational Design and Discovery of Novel Materials (NCCR MARVEL)), Switzerland and ROBERT BOSCH LLC, USA. All rights reserved."
__license__ = "MIT license, see LICENSE.txt file"
__version__ = "0.7.0"
__contributors__ = "Andrea Cepellotti, Victor Garcia-Suarez, Alberto Garcia"

import sys
import os

from aiida.common.example_helpers import test_and_get_code
from aiida.common.exceptions import NotExistent

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
    codename = 'Siesta-4.0@rinaldo'

code = test_and_get_code(codename, expected_code_type='siesta.siesta')
#
#  Set up calculation object first
#
calc = code.new_calc()
calc.label = "Test Siesta. Benzene molecule"
calc.description = "Test calculation with the Siesta code. Benzene molecule"

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
alat = 15. # angstrom
cell = [[alat, 0., 0.,],
        [0., alat, 0.,],
        [0., 0., alat,],
       ]

# Benzene molecule
# Note an atom tagged (for convenience) with a different label
#
s = StructureData(cell=cell)
s.append_atom(position=(0.000,0.000,0.468),symbols=['H'])
s.append_atom(position=(0.000,0.000,1.620),symbols=['C'])
s.append_atom(position=(0.000,-2.233,1.754),symbols=['H'])
s.append_atom(position=(0.000,2.233,1.754),symbols=['H'])
s.append_atom(position=(0.000,-1.225,2.327),symbols='C',name="Cred")
s.append_atom(position=(0.000,1.225,2.327),symbols=['C'])
s.append_atom(position=(0.000,-1.225,3.737),symbols=['C'])
s.append_atom(position=(0.000,1.225,3.737),symbols=['C'])
s.append_atom(position=(0.000,-2.233,4.311),symbols=['H'])
s.append_atom(position=(0.000,2.233,4.311),symbols=['H'])
s.append_atom(position=(0.000,0.000,4.442),symbols=['C'])
s.append_atom(position=(0.000,0.000,5.604),symbols=['H'])

elements = list(s.get_symbols_set())
calc.use_structure(s)
#-------------------------------------------------------------

#
# Parameters ---------------------------------------------------
#
ldos_block_content="\n {e1} {e2} eV".format(e1=-5.0,e2=1.0)

params_dict= {
'xc-functional': 'LDA',
'xc-authors': 'CA',
'mesh-cutoff': '200.000 Ry',
'dm-numberpulay': 5,
'dm-mixingweight': 0.050,
'dm-tolerance': 1.e-4,
'electronic-temperature': '100.000 K',
'%block local-density-of-states': ldos_block_content,
}
parameters = ParameterData(dict=params_dict)
calc.use_parameters(parameters)
#
# Basis Set Info ------------------------------------------
# The basis dictionary follows the 'parameters' convention
#
basis_dict = {
'pao-basistype': 'split',
'pao-splitnorm': 0.150,
'pao-energyshift': '0.020 Ry',
'%block pao-basis-sizes' :"""
C    SZP
Cred SZ
H    SZP  """,
}
basis = ParameterData(dict=basis_dict)
calc.use_basis(basis)
#--------------------------------------------------------------


# Pseudopotentials ----------------------------------------------
#
# This exemplifies the handling of pseudos for different species
# Those sharing the same pseudo should be indicated.
# Families support is not yet available for this.
#
raw_pseudos = [ ("C.psf", ['C', 'Cred']),
                ("H.psf", 'H')]

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

## For remote codes, it is not necessary to manually set the computer,
## since it is set automatically by new_calc
#computer = code.get_remote_computer()
#calc = code.new_calc(computer=computer)

calc.set_max_wallclock_seconds(30*60) # 30 min

calc.set_resources({"num_machines": 1, "num_mpiprocs_per_machine": 2})
code_mpi_enabled =  False
try:
    code_mpi_enabled =  code.get_extra("mpi")
except AttributeError:
    pass
calc.set_withmpi(code_mpi_enabled)
#------------------
queue = None
# calc.set_queue_name(queue)
#------------------

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

