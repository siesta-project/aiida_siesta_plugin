#!/usr/bin/env runaiida
# -*- coding: utf-8 -*-

from __future__ import absolute_import
from __future__ import print_function
import six
__copyright__ = u"Copyright (c), 2015, ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE (Theory and Simulation of Materials (THEOS) and National Centre for Computational Design and Discovery of Novel Materials (NCCR MARVEL)), Switzerland and ROBERT BOSCH LLC, USA. All rights reserved."
__license__ = "MIT license, see LICENSE.txt file"
__version__ = "0.7.0"
__contributors__ = "Andrea Cepellotti, Victor Garcia-Suarez, Alberto Garcia, Emanuele Bosoni"

import sys
import os

import os.path as op
from aiida.engine import run
from aiida.orm import load_code
from aiida.common import NotExistent
from aiida_siesta.calculations.siesta import SiestaCalculation
from aiida.plugins import DataFactory

##########################################################
#                                                        #
#  Siesta calculation on benzene molecule, first to try  #
#                                                        #
##########################################################

PsfData = DataFactory('siesta.psf')
Dict = DataFactory('dict')
KpointsData = DataFactory('array.kpoints')
StructureData = DataFactory('structure')

##########################################################
# Unfortunately I'm not aware of any way to submit tests #
# that only create the folder, but don't store in the    #
# database. In other words, I don't think there is any-  #
# thing similar to the old submit test                   #
##########################################################

#try:
#    dontsend = sys.argv[1]
#    if dontsend == "--dont-send":
#        submit_test = True
#    elif dontsend == "--send":
#        submit_test = False
#    else:
#        raise IndexError
#except IndexError:
#    print(("The first parameter can only be either "
#                          "--send or --dont-send"), file=sys.stderr)
#    sys.exit(1)
#
try:
    codename = sys.argv[1]
except IndexError:
    codename = 'Siesta4.0.1@kelvin'

code=load_code(codename)



options = {
    "queue_name" : "debug",
    "max_wallclock_seconds" : 360,
    "resources" : {
        "num_machines": 1,
        "num_mpiprocs_per_machine": 1,
    }
}



# A Siesta executable compiled in serial mode might not work properly
# on a computer set up for MPI operation.
# This snippet can be used to check whether a code has been compiled
# with mpi support, and act accordingly
# For this to work, the user has to manually add the record in the
# database. In the verdi shell:
#
# code = load_node(code_PK)
# code.set_extra("mpi",True)
#---------------------------------------------------
#
#code_mpi_enabled =  False
#try:
#    code_mpi_enabled =  code.get_extra("mpi")
#except AttributeError:
#    pass
#calc.set_withmpi(code_mpi_enabled)
#------------------



#
#--------- Settings ---------------------------------
#
settings_dict={'additional_retrieve_list': ['aiida.BONDS', 'aiida.EIG']}
settings = Dict(dict=settings_dict)
#----------------------------------------------------

#
#--------- Structure --------------------------------
#
alat = 15. # angstrom
cell = [[alat, 0., 0.,],
        [0., alat, 0.,],
        [0., 0., alat,],
       ]

# Note an atom tagged (for convenience) with a different label

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
#-------------------------------------------------------------

#
# -------------Parameters -------------------------------------
#
# Note the use of '.' in some entries. This will be fixed below.
# Note also that some entries have ':' as separator. This is not
# allowed in Siesta, and will be fixed by the plugin itself. The
# latter case is an unfortunate historical choice. It should not
# be used in modern scripts.
#
params_dict= {
'xc-functional': 'LDA',
'xc-authors': 'CA',
'spin-polarized': True,
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
'write-mulliken-pop': 1,
}
#
parameters = Dict(dict=params_dict)
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
C    SZP
Cred SZ
H    SZP
%endblock pao-basis-sizes""",
}
#
# basis_dict = { k.replace('.','-') :v for k,v in  six.iteritems(basis_dict) }
#
basis = Dict(dict=basis_dict)
#--------------------------------------------------------------


# Pseudopotentials ----------------------------------------------
#
# This exemplifies the handling of pseudos for different species
# Those sharing the same pseudo should be indicated.
#
pseudos_list = []
raw_pseudos = [ ("C.psf", ['C', 'Cred']),
                ("H.psf", 'H')]

for fname, kinds, in raw_pseudos:
    absname = os.path.realpath(os.path.join(os.path.dirname(__file__),
                                            "data",fname))
    pseudo, created = PsfData.get_or_create(absname,use_first=True)
    if created:
        print("Created the pseudo for {}".format(kinds))
    else:
        print("Using the pseudo for {} from DB: {}".format(kinds,pseudo.pk))
    pseudos_list.append(pseudo)

#-------------------------------------------------------------------

####### Needed in new version??
#####calc.set_resources({"parallel_env": 'mpi', "tot_num_mpiprocs": 1})


inputs = {
    'structure': s,
    'parameters': parameters,
    'code': code,
    'basis': basis,
    'pseudos': {
        'C': pseudos_list[0],
        'Cred': pseudos_list[0],
        'H': pseudos_list[1],
    },
    'metadata': {
        'options': options,
        'label' : "Benzene molecule",
    }
}

run(SiestaCalculation, **inputs)

#if submit_test:
#    subfolder, script_filename = calc.submit_test()
#    print("Test_submit for calculation (uuid='{}')".format(
#        calc.uuid))
#    print("Submit file in {}".format(os.path.join(
#        os.path.relpath(subfolder.abspath),
#        script_filename
#        )))
#else:
#    calc.store_all()
#    print("created calculation; calc=Calculation(uuid='{}') # ID={}".format(
#        calc.uuid,calc.dbnode.pk))
#    calc.submit()
#    print("submitted calculation; calc=Calculation(uuid='{}') # ID={}".format(
#        calc.uuid,calc.dbnode.pk))
