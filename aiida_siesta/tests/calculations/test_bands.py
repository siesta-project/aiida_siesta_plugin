#!/usr/bin/env runaiida
# -*- coding: utf-8 -*-
from __future__ import absolute_import
from __future__ import print_function
import os.path as op
import sys
from aiida.orm import load_code
import numpy as np


def create_FCC_structure(alat):
    lvecs = np.array([1.0,
                      1.0,
                      1.0,]) * alat

    cell = np.array([[0.0, 0.5, 0.5,],
                     [0.5, 0.0, 0.5,],
                     [0.5, 0.5, 0.0,],])

    return np.multiply(cell, lvecs)


#def test_bands(siesta_develop):
"""Test workfunction runs and outputs results in serial mode."""
from aiida.orm import Code, Dict, StructureData,KpointsData
from aiida_siesta.calculations.siesta import SiestaCalculation
from aiida.plugins import DataFactory
from aiida.tools import get_explicit_kpoints_path

PsfData = DataFactory('siesta.psf')    

# codename = 'siesta@develop'
# code = Code.get_from_string(codename)
#code = siesta_develop['siesta_code']
code = load_code("Siesta4.0.1@kelvin")
options = {
 "max_wallclock_seconds" : 360,
 "resources" : {
     "num_machines": 1,
     "num_mpiprocs_per_machine": 1,
 }
}

# Define MgO structure
alat = 4.117  # MgO lattice constant, Angstroms
cell = create_FCC_structure(alat)  # Creating MgO FCC-cell
structure = StructureData(cell=cell)  # Creating structure from cell
# Placing basis atoms
structure.append_atom(
    position=(0.000 * alat, 0.000 * alat, 0.000 * alat), symbols=['Mg'])
structure.append_atom(
    position=(0.500 * alat, 0.500 * alat, 0.500 * alat), symbols=['O'])

seekpath_parameters = Dict(dict={'reference_distance': 0.02,'symprec': 0.0001})
result=get_explicit_kpoints_path(structure, **seekpath_parameters.get_dict())
newstructure = result['primitive_structure']

# Pseudopotentials
# from aiida_siesta.data.psf import get_pseudos_from_structure
# inputs.pseudo = get_pseudos_from_structure(structure, "test_psf_family")
raw_pseudos = [
    ("Mg.psf", 'Mg'),
    ("O.psf", 'O'),
]

pseudos_list = []
for fname, kind in raw_pseudos:
    absname = op.realpath(op.join(op.dirname(__file__), '..', 'pseudos', fname))
    pseudo, created = PsfData.get_or_create(absname, use_first=True)

    if created:
        print("\nCreated the pseudo for {}".format(kind))
    else:
        print("\nUsing the pseudo for {} from DB: {}".format(kind, pseudo.pk))
    # Attach pseudo node to the calculation
    pseudos_list.append(pseudo)

# K-points mesh
kpoints_mesh = KpointsData()
kpoints_mesh.set_kpoints_mesh([6, 6, 6], [0.5, 0.5, 0.5])

# Bands' k-points
bandskpoints = KpointsData()
bandskpoints=result['explicit_kpoints']

#bandskpoints.set_cell(structure.cell, structure.pbc)
#bandskpoints.set_kpoints_path([
#    ('K', 'G', 39),
#    ('G', 'X', 37),
#    ('X', 'W', 19),
#    ('W', 'L', 27),
#    ('L', 'G', 32),
#])

# Calculation parameters
parameters = Dict(
    dict={
        'xc-functional': 'LDA',
        'xc-authors': 'CA',
        'spin-polarized': False,
        'meshcutoff': '200 Ry',
        'max-scfiterations': 50,
#        'xml-write': True,
    })

# Create and run Siesta calculation process
from aiida.engine import run

inputs = {
 'structure': newstructure,
 'parameters': parameters,
 'code': code,
 # 'basis': basis,
 'kpoints': kpoints_mesh,
 'bandskpoints' : bandskpoints,
 'pseudos': {
     'Mg': pseudos_list[0],
     'O': pseudos_list[1],
 },
 'metadata': {
     "dry_run": True,
     "store_provenance": False,
     'options': options,
 }
}

#result = run(SiestaCalculation, **inputs)
run(SiestaCalculation, **inputs)
#assert result['bands_array'] is not None
