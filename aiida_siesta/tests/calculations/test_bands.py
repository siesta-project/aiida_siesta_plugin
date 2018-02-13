#!/usr/bin/env runaiida
# -*- coding: utf-8 -*-
import os
import sys

import numpy as np


def create_FCC_structure(alat):
    lvecs = np.array([1.0,
                      1.0,
                      1.0,]) * alat

    cell = np.array([[0.0, 0.5, 0.5,],
                     [0.5, 0.0, 0.5,],
                     [0.5, 0.5, 0.0,],])

    return np.multiply(cell, lvecs)


def test_bands(siesta_develop):
    """Test workfunction runs and outputs results in serial mode."""
    from aiida.orm import Code, DataFactory, CalculationFactory
    from aiida.orm.data.base import Float, Str

    SiestaCalculation = CalculationFactory('siesta.siesta')
    PsfData = DataFactory('siesta.psf')
    StructureData = DataFactory('structure')
    ParameterData = DataFactory('parameter')
    KpointsData = DataFactory('array.kpoints')

    codename = 'siesta@develop'
    code = Code.get_from_string(codename)

    inputs = SiestaCalculation.process().get_inputs_template()
    inputs.code = code
    inputs._options.resources = {
        "num_machines": 1,
        "num_mpiprocs_per_machine": 1,
    }
    inputs._options.max_wallclock_seconds = 30 * 60

    # Define MgO structure
    alat = 4.117  # MgO lattice constant, Angstroms
    cell = create_FCC_structure(alat)  # Creating MgO FCC-cell
    structure = StructureData(cell=cell)  # Creating structure from cell
    # Placing basis atoms
    structure.append_atom(
        position=(0.000 * alat, 0.000 * alat, 0.000 * alat), symbols=['Mg'])
    structure.append_atom(
        position=(0.500 * alat, 0.500 * alat, 0.500 * alat), symbols=['O'])
    inputs.structure = structure

    # Pseudopotentials
    # from aiida_siesta.data.psf import get_pseudos_from_structure
    # inputs.pseudo = get_pseudos_from_structure(structure, "test_psf_family")
    raw_pseudos = [
        ("Mg.psf", 'Mg'),
        ("O.psf", 'O'),
    ]
    pseudo_dict = {}
    for fname, kind in raw_pseudos:
        absname = os.path.realpath(
            os.path.join(os.path.dirname(__file__), '..', 'pseudos', fname))
        pseudo, created = PsfData.get_or_create(absname, use_first=True)

        if created:
            print "Created the pseudo for {}".format(kind)
        else:
            print "Using the pseudo for {} from DB: {}".format(kind, pseudo.pk)
        # Attach pseudo node to the calculation
        pseudo_dict[kind] = pseudo

    inputs.pseudo = pseudo_dict

    # K-points mesh
    kpoints_mesh = KpointsData()
    kpoints_mesh.set_kpoints_mesh([6, 6, 6], [0.5, 0.5, 0.5])
    inputs.kpoints = kpoints_mesh

    # Bands' k-points
    bandskpoints = KpointsData()
    bandskpoints.set_cell(structure.cell, structure.pbc)
    bandskpoints.set_kpoints_path([
        ('K', 'G', 39),
        ('G', 'X', 37),
        ('X', 'W', 19),
        ('W', 'L', 27),
        ('L', 'G', 32),
    ])
    inputs.bandskpoints = bandskpoints

    # Calculation parameters
    parameters = ParameterData(
        dict={
            'xc-functional': 'LDA',
            'xc-authors': 'CA',
            'spin-polarized': False,
            'meshcutoff': '200 Ry',
            'max-scfiterations': 50,
            'xml-write': True,
        })
    inputs.parameters = parameters

    # Create and run Siesta calculation process
    from aiida.work.run import run

    JobCalc = SiestaCalculation.process()
    result = run(JobCalc, **inputs)

    assert result['bands_array'] is not None
