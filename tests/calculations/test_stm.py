#!/usr/bin/env runaiida
# -*- coding: utf-8 -*-
import os.path as op

from aiida import orm
from aiida.common import datastructures
from aiida.tools import get_explicit_kpoints_path
import pytest


def test_base(aiida_profile, fixture_sandbox, fixture_localhost, generate_calc_job,
    fixture_code, file_regression, generate_remote_data):
    """
    Test that a single STM calculation is submitted.
    """

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
    assert sorted(calc_info.remote_copy_list) == sorted(remote_copy_list)
    assert sorted(calc_info.retrieve_list) == sorted(retrieve_list)

    with fixture_sandbox.open('stm.in') as handle:
        input_written = handle.read()
    # Checks on the files written to the sandbox folder as raw input
    # Here it bothers me. Why Si.psf and _aiidasubmit.sh are not in sandbox?
    assert sorted(fixture_sandbox.get_content_list()) == sorted(['stm.in'])
    file_regression.check(input_written, encoding='utf-8', extension='.fdf')



def test_validators(aiida_profile, fixture_sandbox, fixture_localhost, generate_calc_job,
    fixture_code, file_regression, generate_remote_data):
    """
    Test validators STM calculation is submitted.
    """

    entry_point_name = 'siesta.stm'
    remote_ldos_folder = generate_remote_data(fixture_localhost, "/tmp/whatever", "siesta.siesta")

    inputs = {
        'code': fixture_code(entry_point_name),
        'ldos_folder' : remote_ldos_folder,
        'value' : orm.Float(2),
        'metadata': {
            'options': {
               'resources': {'num_machines': 1  },
               'max_wallclock_seconds': 1800,
               'withmpi': False,
               }
        }
    }

    inputs["mode"] = orm.Str("wrong")
    with pytest.raises(ValueError):
        calc_info = generate_calc_job(fixture_sandbox, entry_point_name, inputs)
    inputs["mode"] = orm.Str("constant-height")

    inputs["spin_option"] = orm.Str("wrong")
    with pytest.raises(ValueError):
        calc_info = generate_calc_job(fixture_sandbox, entry_point_name, inputs)
