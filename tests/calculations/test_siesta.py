#!/usr/bin/env runaiida
# -*- coding: utf-8 -*-
import os.path as op
import pytest
from aiida import orm
from aiida.common import datastructures
from aiida.tools import get_explicit_kpoints_path

#from aiida_quantumespresso.utils.resources import get_default_options


def test_base(aiida_profile, fixture_sandbox, generate_calc_job, 
    fixture_code, generate_structure, generate_kpoints_mesh, generate_basis,
    generate_param, generate_psf_data, generate_psml_data, file_regression):
    """
    Test that single calculation is submitted with the right content of the 
    aiida.fdf file.
    """

    entry_point_name = 'siesta.siesta'

    psf = generate_psf_data('Si')
    psml = generate_psml_data('Si')

    inputs = {
        'code': fixture_code(entry_point_name),
        'structure': generate_structure(),
        'kpoints': generate_kpoints_mesh(2),
        'parameters': generate_param(),
        'basis': generate_basis(),
        'pseudos': {
            'Si': psf,
            'SiDiff': psml
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
    local_copy_list = [(psf.uuid, psf.filename, 'Si.psf'),(psml.uuid, psml.filename,'SiDiff.psml')]
    retrieve_list = ['aiida.EIG','aiida.KP','MESSAGES','time.json','aiida.out','aiida.xml','*.ion.xml']
    
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


# Still in progress ...
#def test_restart(aiida_profile, fixture_sandbox, fixture_localhost, generate_calc_job, 
#    fixture_code, generate_structure, generate_kpoints_mesh, generate_basis,
#    generate_param, generate_psf_data, generate_psml_data, file_regression):
#    """Test that single calculation is submitted."""
#
#    entry_point_name = 'siesta.siesta'
#
#    psf = generate_psf_data('Si')
#    psml = generate_psml_data('Si')
#
#    import os
#    basepath = os.path.dirname(os.path.abspath(__file__))
#    fake_remote_path=os.path.join(basepath, 'fixtures', 'restart') 
#    fake_remote_folder = orm.RemoteData(computer=fixture_localhost, remote_path=fake_remote_path)
#
#    #remote_folder.add_incoming(node, link_type=LinkType.CREATE, link_label='remote_folder')
#    #remote_folder.store()
#
#    inputs = {
#        'parent_calc_folder': fake_remote_folder,
#        'code': fixture_code(entry_point_name),
#        'structure': generate_structure(),
#        'kpoints': generate_kpoints_mesh(2),
#        'parameters': generate_param(),
#        'basis': generate_basis(),
#        'pseudos': {
#            'Si': psf,
#            'SiDiff': psml
#        },
#        'metadata': {
#            'options': {
#               'resources': {'num_machines': 1  },
#               'max_wallclock_seconds': 1800,
#               'withmpi': False,
#               }
#        }
#    }
#
#    calc_info = generate_calc_job(fixture_sandbox, entry_point_name, inputs)
#    remote_copy_list = ["as.DM"]
#    assert sorted(calc_info.remote_copy_list) == sorted(remote_copy_list)

# Should be extended to all the blocked_keyworld
def test_blocked_keyword(aiida_profile, fixture_sandbox, generate_calc_job, 
    fixture_code, generate_structure, generate_kpoints_mesh, generate_basis,
    generate_param, generate_psf_data, generate_psml_data, file_regression):
    """
    Test that the plugin is able to detect the forbidden keyword 'system-name'
    and "pao".
    """

    entry_point_name = 'siesta.siesta'

    psf = generate_psf_data('Si')
    psml = generate_psml_data('Si')

    parameters = generate_param()
    parameters.set_attribute('system-name',"whatever")
    parameters.set_attribute('pao-sp',"whatever")

    inputs = {
        'code': fixture_code(entry_point_name),
        'structure': generate_structure(),
        'kpoints': generate_kpoints_mesh(2),
        'parameters': parameters,
        'basis': generate_basis(),
        'pseudos': {
            'Si': psf,
            'SiDiff': psml
        },
        'metadata': {
            'options': {
               'resources': {'num_machines': 1  },
               'max_wallclock_seconds': 1800,
               'withmpi': False,
               }
        }
    }

    from aiida.common import InputValidationError
    import pytest
    with pytest.raises(InputValidationError):
        calc_info = generate_calc_job(fixture_sandbox, entry_point_name, inputs)


def test_cell_consistency(aiida_profile, fixture_sandbox, generate_calc_job,
    fixture_code, generate_structure, generate_kpoints_mesh, generate_basis,
    generate_param, generate_psf_data, generate_psml_data, file_regression):
    """
    Tests that is forbidden to define a cell in bandskpoints
    different from the structure in inputs.
    """

    entry_point_name = 'siesta.siesta'

    psf = generate_psf_data('Si')
    psml = generate_psml_data('Si')

    parameters = generate_param()

    s=generate_structure()
    seekpath_parameters = orm.Dict(dict={
    'reference_distance': 0.02,
    'symprec': 0.0001
    })
    result = get_explicit_kpoints_path(s, **seekpath_parameters.get_dict())
    structure = result['primitive_structure']

    bandskpoints = orm.KpointsData()
    bandskpoints = result['explicit_kpoints']
    bandskpoints.set_cell([[1., 1., 0.0], [0.0, 2., 2.], [2., 0.0, 2.]])

    inputs = {
        'bandskpoints' : bandskpoints,
        'code': fixture_code(entry_point_name),
        'structure': generate_structure(),
        'kpoints': generate_kpoints_mesh(2),
        'parameters': parameters,
        'basis': generate_basis(),
        'pseudos': {
            'Si': psf,
            'SiDiff': psml
        },
        'metadata': {
            'options': {
               'resources': {'num_machines': 1  },
               'max_wallclock_seconds': 1800,
               'withmpi': False,
               }
        }
    }

    import pytest
    with pytest.raises(ValueError):
        calc_info = generate_calc_job(fixture_sandbox, entry_point_name, inputs)


def test_bandslines(aiida_profile, fixture_sandbox, generate_calc_job,
    fixture_code, generate_structure, generate_kpoints_mesh, generate_basis,
    generate_param, generate_psf_data, generate_psml_data, file_regression):
    """
    Test that the fdf file contains the correct information when a band structure
    calculation is required by the user and he/she specifies a band path.
    """

    entry_point_name = 'siesta.siesta'

    psf = generate_psf_data('Si')
    psml = generate_psml_data('Si')

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
            'Si': psf,
            'SiDiff': psml
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

    retrieve_list = ['aiida.EIG','aiida.KP','MESSAGES','time.json','aiida.out','aiida.xml','aiida.bands','*.ion.xml']

    assert sorted(calc_info.retrieve_list) == sorted(retrieve_list)

    with fixture_sandbox.open('aiida.fdf') as handle:
        input_written = handle.read()

    file_regression.check(input_written, encoding='utf-8', extension='.fdf')

def test_bandspoints(aiida_profile, fixture_sandbox, generate_calc_job,
    fixture_code, generate_structure, generate_kpoints_mesh, generate_basis,
    generate_param, generate_psf_data, generate_psml_data, file_regression):
    """
    Test that the fdf file contains the correct information when a band structure
    calculation is required by the user and he/she specifies single k-points.
    """

    entry_point_name = 'siesta.siesta'

    psf = generate_psf_data('Si')
    psml = generate_psml_data('Si')

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
            'Si': psf,
            'SiDiff': psml
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

    retrieve_list = ['aiida.EIG','aiida.KP','MESSAGES','time.json','aiida.out','aiida.xml','aiida.bands','*.ion.xml']

    assert sorted(calc_info.retrieve_list) == sorted(retrieve_list)

    with fixture_sandbox.open('aiida.fdf') as handle:
        input_written = handle.read()

    file_regression.check(input_written, encoding='utf-8', extension='.fdf')


def test_floating_orbs(aiida_profile, fixture_sandbox, generate_calc_job, 
    fixture_code, generate_structure, generate_kpoints_mesh, generate_basis,
    generate_param, generate_psf_data, generate_psml_data, file_regression):
    """
    Test that single calculation is submitted with the right content of the 
    aiida.fdf file.
    """

    entry_point_name = 'siesta.siesta'

    psf = generate_psf_data('Si')
    psml = generate_psml_data('Si')

    inputs = {
        'code': fixture_code(entry_point_name),
        'structure': generate_structure(),
        'kpoints': generate_kpoints_mesh(2),
        'parameters': generate_param(),
        'pseudos': {
            'Si': psf,
            'SiDiff': psml
        },
        'metadata': {
            'options': {
               'resources': {'num_machines': 1  },
               'max_wallclock_seconds': 1800,
               'withmpi': False,
               }
        }
    }

    #fail because no list
    basis = generate_basis().get_dict()
    basis["floating_sites"] = "www"
    inputs['basis'] = orm.Dict(dict=basis)
    with pytest.raises(ValueError):
        calc_info = generate_calc_job(fixture_sandbox, entry_point_name, inputs)

    #fail because no dictionaries
    basis = generate_basis().get_dict()
    basis["floating_sites"] = ["ddd","www"]
    inputs['basis'] = orm.Dict(dict=basis)
    with pytest.raises(ValueError):
        calc_info = generate_calc_job(fixture_sandbox, entry_point_name, inputs)

    #fail because missing name
    basis = generate_basis().get_dict()
    basis["floating_sites"] = [{"symbols":'Si',"position": ( 0.125, 0.125, 0.125)}]
    inputs['basis'] = orm.Dict(dict=basis)
    with pytest.raises(ValueError):
        calc_info = generate_calc_job(fixture_sandbox, entry_point_name, inputs)

    #fail because missing symbols
    basis = generate_basis().get_dict()
    basis["floating_sites"] = [{"name":'SiDiff',"position": ( 0.125, 0.125, 0.125)}]
    inputs['basis'] = orm.Dict(dict=basis)
    with pytest.raises(ValueError):
        calc_info = generate_calc_job(fixture_sandbox, entry_point_name, inputs)

    #fail because missing position
    basis = generate_basis().get_dict()
    basis["floating_sites"] = [{"name":'SiDiff',"symbols":'Si'}]
    inputs['basis'] = orm.Dict(dict=basis)
    with pytest.raises(ValueError):
        calc_info = generate_calc_job(fixture_sandbox, entry_point_name, inputs)

    #fail because same name of a site in real structure
    basis = generate_basis().get_dict()
    basis["floating_sites"] = [{"name":'SiDiff',"symbols":'Si',"position": ( 0.125, 0.125, 0.125)}]
    inputs['basis'] = orm.Dict(dict=basis)
    with pytest.raises(ValueError):
        calc_info = generate_calc_job(fixture_sandbox, entry_point_name, inputs)
    
    #fail because no pseudo:
    basis = generate_basis().get_dict()
    basis["floating_sites"] = [{"name":'Si_bond',"symbols":'Si',"position": ( 0.125, 0.125, 0.125)}] 
    inputs['basis'] = orm.Dict(dict=basis)
    with pytest.raises(ValueError):
        calc_info = generate_calc_job(fixture_sandbox, entry_point_name, inputs)

    inputs['pseudos']['Si_bond'] = psf
    calc_info = generate_calc_job(fixture_sandbox, entry_point_name, inputs)

    with fixture_sandbox.open('aiida.fdf') as handle:
        input_written = handle.read()

    file_regression.check(input_written, encoding='utf-8', extension='.fdf')

