#!/usr/bin/env runaiida
# -*- coding: utf-8 -*-

from __future__ import absolute_import
from __future__ import print_function

import sys

import pymatgen as mg

from aiida.engine import submit
from aiida.orm import load_code
from aiida_siesta.calculations.siesta import SiestaCalculation
from aiida_siesta.data.psf import get_pseudos_from_structure
from aiida.plugins import DataFactory
from aiida.tools import get_explicit_kpoints_path

# This script will send a Siesta calculation on a structure taken from
# a cif file.
# The band structure is calculated and the kpoint path is automatically
# generated using seekpath.

PsfData = DataFactory('siesta.psf')
Dict = DataFactory('dict')
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
    print(("The first parameter can only be either "
           "--send or --dont-send"),
          file=sys.stderr)
    sys.exit(1)

try:
    codename = sys.argv[2]
except IndexError:
    codename = 'siesta4.0.1@parsons'

#
#------------------Code and computer options ---------------------------
#
code = load_code(codename)

options = {
    "queue_name": "debug",
    "max_wallclock_seconds": 1700,
    "resources": {
        "num_machines": 1,
        "num_mpiprocs_per_machine": 1,
    }
}
#
settings_dict = {'additional_retrieve_list': ['aiida.BONDS', 'aiida.EIG']}
settings = Dict(dict=settings_dict)
#---------------------------------------------------------------------

#
# Structure -----------------------------------------
#
# Passing through SeeK-path first, to get the standardized cell.
# Necessary for the automatic choice of the bands path.
#
structure = mg.Structure.from_file("data/O2_ICSD_173933.cif", primitive=False)
s = StructureData(pymatgen_structure=structure)

seekpath_parameters = {'reference_distance': 0.02, 'symprec': 0.0001}
result = get_explicit_kpoints_path(s, **seekpath_parameters)
newstructure = result['primitive_structure']

#
# Parameters ---------------------------------------------------
#
params_dict = {
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
    'write-mulliken-pop': 1
}
#
parameters = Dict(dict=params_dict)

#----------------------------------------------------------
#
# Basis Set Info ------------------------------------------
# The basis dictionary follows the 'parameters' convention
#
basis_dict = {
    'pao-basistype': 'split',
    'pao-splitnorm': 0.150,
    'pao-energyshift': '0.020 Ry',
    '%block pao-basis-sizes': """
O    SZP  
%endblock pao-basis-sizes""",
}
#
basis = Dict(dict=basis_dict)
#--------------------------------------------------------------

#--------------------- Pseudopotentials ---------------------------------
#
# FIXME: The family name is hardwired
#
pseudos_dict = get_pseudos_from_structure(s, 'lda-ag')
print(pseudos_dict)
#-----------------------------------------------------------------------

# K-points for scf cycle -------------------------------------------
kts = KpointsData()
kpoints_mesh = 4
kts.set_kpoints_mesh([kpoints_mesh, kpoints_mesh, kpoints_mesh])

#
bandskpoints = KpointsData()

# Making use of SeeK-path for the automatic path
# The choice of the distance between kpoints is in the call seekpath_parameters
# All high symmetry points included, labels already included
bandskpoints = result['explicit_kpoints']

#
#--All the inputs of a Siesta calculations are listed in a dictionary--
#
inputs = {
    'structure': newstructure,
    'parameters': parameters,
    'code': code,
    'basis': basis,
    'kpoints': kts,
    'bandskpoints': bandskpoints,
    'pseudos': pseudos_dict,
    'metadata': {
        'options': options,
        'label': "O_el_cell_spin from CIF"
    }
}

if submit_test:
    inputs["metadata"]["dry_run"] = True
    inputs["metadata"]["store_provenance"] = False
    process = submit(SiestaCalculation, **inputs)
    print("Submited test for calculation (uuid='{}')".format(process.uuid))
    print("Check the folder submit_test for the result of the test")

else:
    process = submit(SiestaCalculation, **inputs)
    print("Submitted calculation; ID={}".format(process.pk))
    print("For information about this calculation type: verdi process show {}".
          format(process.pk))
    print("For a list of running processes type: verdi process list")
