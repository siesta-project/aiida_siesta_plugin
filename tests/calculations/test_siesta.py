#!/usr/bin/env runaiida
# -*- coding: utf-8 -*-
from __future__ import absolute_import
from __future__ import print_function
import os.path as op
from aiida import orm
from aiida.common import datastructures
from aiida.tools import get_explicit_kpoints_path

#from aiida_quantumespresso.utils.resources import get_default_options


def test_base(aiida_profile, fixture_sandbox, generate_calc_job, 
    fixture_code, generate_structure, generate_kpoints_mesh, generate_basis,
    generate_param, generate_psf_data, file_regression):
    """Test that single calculation is submitted."""

    entry_point_name = 'siesta.siesta'

    psf = generate_psf_data('Si')

    inputs = {
        'code': fixture_code(entry_point_name),
        'structure': generate_structure(),
        'kpoints': generate_kpoints_mesh(2),
        'parameters': generate_param(),
        'basis': generate_basis(),
        'pseudos': {
            'Si': psf
        },
        'metadata': {
            'options': {
               'resources': {'num_machines': 1  },
               'max_wallclock_seconds': 1800,
               'withmpi': False,
               }
        }
    }

    calc_info = generate_calc_job(fixture_sandbox, entry_point_name, inputs)

    #cmdline_params = ['-in', 'aiida.fdf']
    local_copy_list = [(psf.uuid, psf.filename, 'Si.psf')]
    retrieve_list = ['MESSAGES', 'time.json', 'aiida.out', 'aiida.xml']

    # Check the attributes of the returned `CalcInfo`
    assert isinstance(calc_info, datastructures.CalcInfo)
    #assert sorted(calc_info.cmdline_params) == sorted(cmdline_params) 
    #we don't modify the command line
    #But maybe we should check the input and output filename
    assert sorted(calc_info.local_copy_list) == sorted(local_copy_list)
    assert sorted(calc_info.retrieve_list) == sorted(retrieve_list)

    with fixture_sandbox.open('aiida.fdf') as handle:
        input_written = handle.read()

    # Checks on the files written to the sandbox folder as raw input
    # Here it bothers me. Why Si.psf and _aiidasubmit.sh are not in sandbox?
    assert sorted(fixture_sandbox.get_content_list()) == sorted(['aiida.fdf'])
    file_regression.check(input_written, encoding='utf-8', extension='.fdf')


def test_bandslines(aiida_profile, fixture_sandbox, generate_calc_job,
    fixture_code, generate_structure, generate_kpoints_mesh, generate_basis,
    generate_param, generate_psf_data, file_regression):
    """Test that single calculation is submitted."""

    entry_point_name = 'siesta.siesta'

    psf = generate_psf_data('Si')

    s=generate_structure()
    seekpath_parameters = orm.Dict(dict={
    'reference_distance': 0.02,
    'symprec': 0.0001
    })
    result = get_explicit_kpoints_path(s, **seekpath_parameters.get_dict())
    structure = result['primitive_structure']

    bandskpoints = orm.KpointsData()
    bandskpoints = result['explicit_kpoints']
    
    
    inputs = {
        'bandskpoints': bandskpoints,
        'code': fixture_code(entry_point_name),
        'structure': structure,
        'kpoints': generate_kpoints_mesh(2),
        'parameters': generate_param(),
        'basis': generate_basis(),
        'pseudos': {
            'Si': psf
        },
        'metadata': {
            'options': {
               'resources': {'num_machines': 1  },
               'max_wallclock_seconds': 1800,
               'withmpi': False,
               }
        }
    }

    calc_info = generate_calc_job(fixture_sandbox, entry_point_name, inputs)

    retrieve_list = ['MESSAGES', 'time.json', 'aiida.out', 'aiida.xml','aiida.bands']

    assert sorted(calc_info.retrieve_list) == sorted(retrieve_list)

    with fixture_sandbox.open('aiida.fdf') as handle:
        input_written = handle.read()

    file_regression.check(input_written, encoding='utf-8', extension='.fdf')

def test_bandspoints(aiida_profile, fixture_sandbox, generate_calc_job,
    fixture_code, generate_structure, generate_kpoints_mesh, generate_basis,
    generate_param, generate_psf_data, file_regression):
    """Test that single calculation is submitted."""

    entry_point_name = 'siesta.siesta'

    psf = generate_psf_data('Si')

    structure=generate_structure()
    bandskpoints = orm.KpointsData()
    kpp = [(0.500,  0.250, 0.750), (0.500,  0.500, 0.500), (0., 0., 0.)]
    bandskpoints.set_cell(structure.cell, structure.pbc)
    bandskpoints.set_kpoints(kpp)


    inputs = {
        'bandskpoints': bandskpoints,
        'code': fixture_code(entry_point_name),
        'structure': structure,
        'kpoints': generate_kpoints_mesh(2),
        'parameters': generate_param(),
        'basis': generate_basis(),
        'pseudos': {
            'Si': psf
        },
        'metadata': {
            'options': {
               'resources': {'num_machines': 1  },
               'max_wallclock_seconds': 1800,
               'withmpi': False,
               }
        }
    }

    calc_info = generate_calc_job(fixture_sandbox, entry_point_name, inputs)

    retrieve_list = ['MESSAGES', 'time.json', 'aiida.out', 'aiida.xml','aiida.bands']

    assert sorted(calc_info.retrieve_list) == sorted(retrieve_list)

    with fixture_sandbox.open('aiida.fdf') as handle:
        input_written = handle.read()

    file_regression.check(input_written, encoding='utf-8', extension='.fdf')

