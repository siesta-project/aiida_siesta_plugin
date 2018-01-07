#!/usr/bin/env runaiida
# -*- coding: utf-8 -*-
import os.path as op


def test_simple_submission(siesta_develop):
    """Test that single calculation is submitted."""
    from aiida.orm import DataFactory
    from aiida.common.example_helpers import test_and_get_code

    PsfData = DataFactory('siesta.psf')
    StructureData = DataFactory('structure')
    ParameterData = DataFactory('parameter')
    KpointsData = DataFactory('array.kpoints')

    code = test_and_get_code('siesta', expected_code_type='siesta.siesta')
    assert code is not None

    # Si diamond structure
    alat = 5.430  # angstrom
    cell = [
        [
            0.5 * alat,
            0.5 * alat,
            0.,
        ],
        [
            0.,
            0.5 * alat,
            0.5 * alat,
        ],
        [
            0.5 * alat,
            0.,
            0.5 * alat,
        ],
    ]

    # Si
    # This was originally given in the "ScaledCartesian" format
    #
    s = StructureData(cell=cell)
    s.append_atom(
        position=(0.000 * alat, 0.000 * alat, 0.000 * alat), symbols=['Si'])
    s.append_atom(
        position=(0.250 * alat, 0.250 * alat, 0.250 * alat), symbols=['Si'])

    parameters = ParameterData(dict={
        'xc-functional': 'LDA',
        'xc-authors': 'CA',
        'max-scfiterations': 50,
        'dm-numberpulay': 4,
        'dm-mixingweight': 0.3,
        'dm-tolerance': 1.e-3,
        'Solution-method': 'diagon',
        'electronic-temperature': '25 meV',
        'md-typeofrun': 'cg',
        'md-numcgsteps': 3,
        'md-maxcgdispl': '0.1 Ang',
        'md-maxforcetol': '0.04 eV/Ang',
        'xml-write': True
    })

    basis = ParameterData(dict={
        'pao-energy-shift': '300 meV',
        '%block pao-basis-sizes': 'Si DZP',
    })

    kpoints = KpointsData()
    kpoints.set_kpoints_mesh([4, 4, 4])

    calc = code.new_calc(
        max_wallclock_seconds=3600,
        resources={
            "num_machines": 1,
            "num_mpiprocs_per_machine": 1,
        })

    calc.label = "Si bulk"
    calc.description = "Test calculation with the Siesta code. Si bulk"

    # Use raw pseudos for this example:
    raw_pseudos = [("Si.psf", 'Si')]

    calc.use_structure(s)
    calc.use_code(code)
    calc.use_parameters(parameters)
    calc.use_basis(basis)
    calc.use_kpoints(kpoints)

    # Pseudo business
    #
    # TODO: understand how to work with pseudo families
    # calc.use_pseudos_from_family(pseudo_family)
    #
    for fname, kind in raw_pseudos:
        absname = op.realpath(op.join(op.dirname(__file__), "..", "pseudos", fname))
        pseudo, created = PsfData.get_or_create(absname, use_first=True)

        if created:
            print "\nCreated the pseudo for {}".format(kind)
        else:
            print "\nUsing the pseudo for {} from DB: {}".format(
                kind, pseudo.pk)

    # Attach pseudo node to the calculation
    calc.use_pseudo(pseudo, kind=kind)

    calc.store_all()
    print "created calculation with PK={}".format(calc.pk)
    subfolder, script_filename = calc.submit_test()
    print "Test_submit for calculation (uuid='{}')".format(calc.uuid)

    assert op.isfile(op.join(op.relpath(subfolder.abspath), script_filename))
    print "Submit file in {}".format(
        op.join(op.relpath(subfolder.abspath), script_filename))

    assert op.isfile(op.join(op.relpath(subfolder.abspath), 'aiida.fdf'))
    assert op.isfile(op.join(op.relpath(subfolder.abspath), 'Si.psf'))
