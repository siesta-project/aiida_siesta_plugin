#!/usr/bin/env runaiida
import os.path as op
from aiida import orm
from aiida.common import datastructures
from aiida.tools import get_explicit_kpoints_path


def test_base(aiida_profile, fixture_sandbox, fixture_localhost, generate_calc_job, 
    fixture_code, file_regression, generate_remote_data):
    """Test that single calculation is submitted."""

    entry_point_name = 'siesta.stm'
    remote_ldos_folder = generate_remote_data(fixture_localhost, "/tmp/whatever", "siesta.siesta")

    inputs = {
        'code': fixture_code(entry_point_name),
        'ldos_folder' : remote_ldos_folder,
        'mode' : orm.Str("constant-current"),
        'value' : orm.Float(2),
        'spin_option' : orm.Str("s"),
        'metadata': {
            'options': {
               'resources': {'num_machines': 1  },
               'max_wallclock_seconds': 1800,
               'withmpi': False,
               }
        }
    }

    calc_info = generate_calc_job(fixture_sandbox, entry_point_name, inputs)
    subf='./'
    remote_copy_list = [(remote_ldos_folder.computer.uuid,
            op.join(remote_ldos_folder.get_remote_path(), subf, '*.LDOS'),
            subf)]
    cmdline_params = ['-i', '2.00000', '-s', 's', 'aiida.LDOS']
    retrieve_list = ['*.STM', 'stm.out']
    
    # Check the attributes of the returned `CalcInfo`
    assert isinstance(calc_info, datastructures.CalcInfo)
    #check command line
    assert calc_info.codes_info[0].cmdline_params == cmdline_params
    #check remote_copy_list and retrieve_list. The local_copy_list is eampty
    #assert sorted(calc_info.local_copy_list) == sorted(local_copy_list)
    assert sorted(calc_info.remote_copy_list) == sorted(remote_copy_list)
    assert sorted(calc_info.retrieve_list) == sorted(retrieve_list)

    with fixture_sandbox.open('stm.in') as handle:
        input_written = handle.read()
    # Checks on the files written to the sandbox folder as raw input
    # Here it bothers me. Why Si.psf and _aiidasubmit.sh are not in sandbox?
    assert sorted(fixture_sandbox.get_content_list()) == sorted(['stm.in'])
    file_regression.check(input_written, encoding='utf-8', extension='.fdf')


#def test_blocked_keyword(aiida_profile, fixture_sandbox, generate_calc_job, 
#    fixture_code, generate_structure, generate_kpoints_mesh, generate_basis,
#    generate_param, generate_psf_data, generate_psml_data, file_regression):
#    """Test that single calculation is submitted."""
#
#    entry_point_name = 'siesta.siesta'
#
#    psf = generate_psf_data('Si')
#    psml = generate_psml_data('Si')
#
#    parameters = generate_param()
#    parameters.set_attribute('system-name',"whatever")
#
#    inputs = {
#        'code': fixture_code(entry_point_name),
#        'structure': generate_structure(),
#        'kpoints': generate_kpoints_mesh(2),
#        'parameters': parameters,
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
#    from aiida.common import InputValidationError
#    import pytest
#    with pytest.raises(InputValidationError):
#        calc_info = generate_calc_job(fixture_sandbox, entry_point_name, inputs)
