#!/usr/bin/env runaiida
# -*- coding: utf-8 -*-
import os.path as op

from aiida import orm
from aiida.common import datastructures
from aiida.tools import get_explicit_kpoints_path
import pytest

#from aiida_quantumespresso.utils.resources import get_default_options


def test_base(aiida_profile, fixture_sandbox, generate_calc_job,
    fixture_code, generate_structure, generate_kpoints_mesh, generate_basis,
    generate_param, generate_psf_data, generate_psml_data, file_regression):
    """
    Test that single calculation is submitted with the right content of the
    aiida.fdf file and with correct lists of options. Also checks that the
    `settings` inputs properly works.
    """

    entry_point_name = 'siesta.siesta'

    psf = generate_psf_data('Si')
    psml = generate_psml_data('Si')

    settings_dict = {'cmdline': ['-option1', '-option2'], 'ADDITIONAL_RETRIEVE_LIST': ["w"]}

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
        'settings': orm.Dict(dict=settings_dict),
        'metadata': {
            'options': {
               'resources': {'num_machines': 1  },
               'max_wallclock_seconds': 1800,
               'withmpi': False,
               }
        }
    }

    calc_info = generate_calc_job(fixture_sandbox, entry_point_name, inputs)

    cmdline_params = ['-option1', '-option2']
    local_copy_list = [(psf.uuid, psf.filename, 'Si.psf'),(psml.uuid, psml.filename,'SiDiff.psml')]
    retrieve_list = [
            'w','BASIS_HARRIS_ENTHALPY','BASIS_ENTHALPY', 'MESSAGES','time.json','aiida.out','aiida.xml','*.ion.xml'
            ]

    # Check the attributes of the returned `CalcInfo`
    assert isinstance(calc_info, datastructures.CalcInfo)
    assert sorted(calc_info.codes_info[0].cmdline_params) == sorted(cmdline_params)
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

def test_validators(aiida_profile, fixture_sandbox, generate_calc_job,
    fixture_code, generate_structure, generate_kpoints_mesh, generate_basis,
    generate_param, generate_psf_data, generate_psml_data, file_regression):
    """
    Test the validators of the input ports
    """

    entry_point_name = 'siesta.siesta'

    psf = generate_psf_data('Si')
    psml = generate_psml_data('Si')

    parameters = generate_param()
    kpoints = orm.KpointsData()

    inputs = {
        'code': fixture_code(entry_point_name),
        'pseudos': {
            'Si': psf,
        #    'SiDiff': psml
        },
        'metadata': {
            'options': {
               'resources': {'num_machines': 1  },
               'max_wallclock_seconds': 1800,
               'withmpi': False,
               }
        }
    }

    # Test the structur validator
    structure = generate_structure().clone()
    structure.append_atom(symbols=["Si","Al"],position=(1.0,2.0,0.0),weights=[0.5,0.5],name="w")
    inputs["structure"] = structure
    with pytest.raises(ValueError):
        calc_info = generate_calc_job(fixture_sandbox, entry_point_name, inputs)
    structure = generate_structure()
    inputs["structure"] = structure

    # Test the parameters validator
    parameters.set_attribute('system-name',"whatever")
    inputs["parameters"] = parameters.clone()
    with pytest.raises(ValueError):
        calc_info = generate_calc_job(fixture_sandbox, entry_point_name, inputs)
    parameters.delete_attribute('system-name')
    parameters.set_attribute('pao-sp',"whatever")
    inputs["parameters"] = parameters.clone()
    with pytest.raises(ValueError):
        calc_info = generate_calc_job(fixture_sandbox, entry_point_name, inputs)
    parameters.delete_attribute('pao-sp')
    inputs["parameters"] = parameters.clone()

    # Test the kpoints validator
    inputs["kpoints"] = kpoints
    with pytest.raises(ValueError):
        calc_info = generate_calc_job(fixture_sandbox, entry_point_name, inputs)
    inputs["kpoints"] = generate_kpoints_mesh(2)

    # Test missing pseudo
    with pytest.raises(ValueError):
        calc_info = generate_calc_job(fixture_sandbox, entry_point_name, inputs)

    # The basis validator is checked in the test of test_floating_orbs.
    # The bandskpoints validator is tested in test_bandslines
    # Test on ions is in the test_ions.
    # Validators on optical tested in test_optical


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

    #bandskpoints = orm.KpointsData()
    #bandskpoints = result['explicit_kpoints']


    inputs = {
        #'bandskpoints': bandskpoints,
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

    # Checks the bandkpoints validator
    bandskpoints = orm.KpointsData()
    inputs["bandskpoints"] = bandskpoints
    with pytest.raises(ValueError):
        generate_calc_job(fixture_sandbox, entry_point_name, inputs)

    inputs["bandskpoints"] = result['explicit_kpoints']

    calc_info = generate_calc_job(fixture_sandbox, entry_point_name, inputs)

    retrieve_list = [
            'BASIS_HARRIS_ENTHALPY','BASIS_ENTHALPY', 'MESSAGES','time.json','aiida.out','aiida.xml','aiida.bands','*.ion.xml'
            ]

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

    retrieve_list = [
            'BASIS_HARRIS_ENTHALPY', 'BASIS_ENTHALPY', 'MESSAGES','time.json','aiida.out','aiida.xml','aiida.bands','*.ion.xml'
            ]

    assert sorted(calc_info.retrieve_list) == sorted(retrieve_list)

    with fixture_sandbox.open('aiida.fdf') as handle:
        input_written = handle.read()

    file_regression.check(input_written, encoding='utf-8', extension='.fdf')


def test_floating_orbs(aiida_profile, fixture_sandbox, generate_calc_job,
    fixture_code, generate_structure, generate_kpoints_mesh, generate_basis,
    generate_param, generate_psf_data, generate_psml_data, file_regression):
    """
    Test all the parts related to the floating_sites: the validation of the
    basis inputs, the creation of the correct aiida.fdf and the copy of the
    correct pseudos.
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


def test_ions(aiida_profile, fixture_sandbox, generate_calc_job,
    fixture_code, generate_structure, generate_kpoints_mesh, generate_basis,
    generate_param, generate_ion_data, generate_psml_data, file_regression):
    """
    Test all the parts related to ion files:
    """

    entry_point_name = 'siesta.siesta'

    ion = generate_ion_data('Si')
    psml = generate_psml_data('Si')

    inputs = {
        'code': fixture_code(entry_point_name),
        'structure': generate_structure(),
        'kpoints': generate_kpoints_mesh(2),
        'parameters': generate_param(),
        'ions': {
            'Si': ion,
    #       'SiDiff': ion
        },
        'metadata': {
            'options': {
               'resources': {'num_machines': 1  },
               'max_wallclock_seconds': 1800,
               'withmpi': False,
               }
        }
    }

    # faile because missing ion for SiDiff
    with pytest.raises(ValueError):
        calc_info = generate_calc_job(fixture_sandbox, entry_point_name, inputs)

    inputs["ions"]['SiDiff'] = ion

    calc_info = generate_calc_job(fixture_sandbox, entry_point_name, inputs)

    with fixture_sandbox.open('aiida.fdf') as handle:
        input_written = handle.read()

    file_regression.check(input_written, encoding='utf-8', extension='.fdf')

    with fixture_sandbox.open('Si.ion') as ion_handle:
        ion_input_written = ion_handle.read()

    file_regression.check(ion_input_written, encoding='utf-8', extension='.ion')



def test_lua(aiida_profile, fixture_sandbox, generate_calc_job,
    fixture_code, generate_structure, generate_kpoints_mesh, generate_basis,
    generate_param, generate_psml_data, generate_lua_file, generate_lua_folder, file_regression):
    """
    Test that single calculation is submitted with the right content of the
    aiida.fdf file.
    """

    entry_point_name = 'siesta.siesta'

    psml = generate_psml_data('Si')

    lua_script = generate_lua_file()
    lua_folder = generate_lua_folder()
    lua_parameters = {
        'number_of_internal_images_in_path': 5,
        'neb_spring_constant': 0.45,
        'neb_image_file_prefix': "image-"
    }
    lua_retrieve_list = ['NEB.results' ]


    inputs = {
        'code': fixture_code(entry_point_name),
        'structure': generate_structure(),
        'kpoints': generate_kpoints_mesh(2),
        'parameters': generate_param(),
        'pseudos': {
            'Si': psml,
            'SiDiff': psml
        },
        'lua': {
            'script': lua_script,
            'input_files': lua_folder,
            'parameters': orm.Dict(dict=lua_parameters),
            'retrieve_list': orm.List(list=lua_retrieve_list)
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

    list_lua_fold = lua_folder.list_object_names()
    local_copy_list = [
            (psml.uuid, psml.filename, 'Si.psml'),
            (psml.uuid, psml.filename, 'SiDiff.psml'),
            (lua_script.uuid, lua_script.filename, lua_script.filename),
            (lua_folder.uuid, list_lua_fold[0], list_lua_fold[0]),
            (lua_folder.uuid, list_lua_fold[1], list_lua_fold[1]),
            (lua_folder.uuid, list_lua_fold[2], list_lua_fold[2]),
            (lua_folder.uuid, list_lua_fold[3], list_lua_fold[3]),
            ]

    retrieve_list = [
            'BASIS_HARRIS_ENTHALPY','BASIS_ENTHALPY', 'MESSAGES','time.json','aiida.out','aiida.xml','*.ion.xml','NEB.results'
            ]

    assert sorted(calc_info.local_copy_list) == sorted(local_copy_list)
    assert sorted(calc_info.retrieve_list) == sorted(retrieve_list)

    with fixture_sandbox.open('aiida.fdf') as handle:
        input_written = handle.read()

    file_regression.check(input_written, encoding='utf-8', extension='.fdf')

    with fixture_sandbox.open('config.lua') as conflua_handle:
        conflua_input_written = conflua_handle.read()

    file_regression.check(conflua_input_written, encoding='utf-8', extension='.lua')


def test_optical(aiida_profile, fixture_sandbox, generate_calc_job,
    fixture_code, generate_structure, generate_kpoints_mesh, generate_basis,
    generate_param, generate_psml_data, file_regression):
    """
    Test that single calculation is submitted with the right content of the
    aiida.fdf file.
    """

    entry_point_name = 'siesta.siesta'

    psml = generate_psml_data('Si')

    inputs = {
        'code': fixture_code(entry_point_name),
        'structure': generate_structure(),
        'kpoints': generate_kpoints_mesh(2),
        'parameters': generate_param(),
        'pseudos': {
            'Si': psml,
            'SiDiff': psml
        },
        'optical': orm.Dict(dict={}),
        'metadata': {
            'options': {
               'resources': {'num_machines': 1  },
               'max_wallclock_seconds': 1800,
               'withmpi': False,
               }
        }
    }

    #Fail because no optical-mesh block
    with pytest.raises(ValueError):
        calc_info = generate_calc_job(fixture_sandbox, entry_point_name, inputs)

    optical_parameters = {"%block optical-mesh" : "\n 1 1 1 \n%endblock optical-mesh"}
    inputs["optical"] = orm.Dict(dict=optical_parameters)
    calc_info = generate_calc_job(fixture_sandbox, entry_point_name, inputs)

    optical_parameters["optical-polarization-type"] = "polarized"
    inputs["optical"] = orm.Dict(dict=optical_parameters)
    #Fail because an optical vector is needed for polarized optical calculation
    with pytest.raises(ValueError):
        calc_info = generate_calc_job(fixture_sandbox, entry_point_name, inputs)

    optical_parameters["%block optical-vector"] = "\n 1.0 0.0 0.0 \n%endblock optical-vector"
    optical_parameters["optical-calculation"] = False
    inputs["optical"] = orm.Dict(dict=optical_parameters)
    calc_info = generate_calc_job(fixture_sandbox, entry_point_name, inputs)

    retrieve_list = ['BASIS_ENTHALPY', 'BASIS_HARRIS_ENTHALPY', 'MESSAGES','time.json','aiida.out','aiida.xml','*.ion.xml',"aiida.EPSIMG"]

    assert sorted(calc_info.retrieve_list) == sorted(retrieve_list)

    with fixture_sandbox.open('aiida.fdf') as handle:
        input_written = handle.read()

    file_regression.check(input_written, encoding='utf-8', extension='.fdf')
