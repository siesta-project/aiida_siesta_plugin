#!/usr/bin/env runaiida
import pytest
from plumpy import ProcessState
from aiida.engine import ExitCode
from aiida import orm
from aiida.common import (LinkType, AttributeDict)
from aiida_siesta.workflows.stm import SiestaSTMWorkChain

@pytest.fixture
def generate_workchain_stm(generate_psml_data, fixture_code, generate_workchain, 
        generate_structure, generate_param, generate_basis, generate_kpoints_mesh,
        generate_calc_job_node, generate_parser):
    """Generate an instance of a `BandgapWorkChain`."""

    def _generate_workchain_stm():

        entry_point_code_siesta = 'siesta.siesta'
        entry_point_code = 'siesta.stm'
        entry_point_wc = 'siesta.stm'

        psml = generate_psml_data('Si')

        structure = generate_structure()

        inputs = {
            'code': fixture_code(entry_point_code_siesta),
            'stm_code' : fixture_code(entry_point_code),
            'stm_mode' : orm.Str("constant-height"),
            'stm_spin' : orm.Str("none"),
            'stm_value' : orm.Float(1),
            'emin' : orm.Float(-1),
            'emax' : orm.Float(1),
            'structure': structure,
            'kpoints': generate_kpoints_mesh(2),
            'parameters': generate_param(),
            'basis': generate_basis(),
            'pseudos': {
                'Si': psml,
                'SiDiff': psml
            },
            'options': orm.Dict(dict={
               'resources': {'num_machines': 1  },
               'max_wallclock_seconds': 1800,
               'withmpi': False,
               })
        }

        process = generate_workchain(entry_point_wc, inputs)

        return process

    return _generate_workchain_stm


def test_setup_end_error_stm(aiida_profile, generate_workchain_stm, 
        generate_wc_job_node, fixture_localhost):
    """Test `SiestaSTMWorkChain`."""
    
    process = generate_workchain_stm()
    process.checks()

    assert process.ctx.spinstm == "none"
    assert not process.ctx.ldosdefinedinparam

    ldos_basewc = generate_wc_job_node("siesta.base", fixture_localhost)
    ldos_basewc.set_process_state(ProcessState.FINISHED)
    #Just like this, the exit status in not zero, therefore not finished_ok
    #ldos_basewc.set_exit_status(ExitCode(0).status)
    process.ctx.siesta_ldos = ldos_basewc

    result = process.run_results()

    #assert not process.is_finished_ok
    assert result == SiestaSTMWorkChain.exit_codes.ERROR_STM_PLUGIN


def test_outputs(aiida_profile, generate_workchain_stm, generate_wc_job_node, 
         generate_calc_job_node, fixture_localhost):
    """Test `SiestaSTMWorkChain`."""

    process = generate_workchain_stm()
    process.checks()

    ldos_basewc = generate_wc_job_node("siesta.base", fixture_localhost)
    ldos_basewc.set_process_state(ProcessState.FINISHED)
    #Now we put the exit code finished_ok
    ldos_basewc.set_exit_status(ExitCode(0).status)
    process.ctx.siesta_ldos = ldos_basewc

    name = 'default'
    entry_point_calc_job = 'siesta.stm'
    inputs = AttributeDict({
        'spin_option' : orm.Str("q")
    })
    attributes=AttributeDict({'input_filename':'stm.in', 'output_filename':'stm.out'})
    stm_node = generate_calc_job_node(entry_point_calc_job, fixture_localhost, name, inputs, attributes)
    stm_array = orm.ArrayData()
    stm_array.add_incoming(stm_node, link_type=LinkType.CREATE, link_label='stm_array')
    stm_array.store()
    process.ctx.stm_calc = stm_node

    first_basewc = generate_wc_job_node("siesta.base", fixture_localhost)
    out_par = orm.Dict(dict={"variable_geometry":False})
    out_par.add_incoming(first_basewc, link_type=LinkType.CREATE, link_label='output_parameters')
    out_par.store()
    process.ctx.workchain_base = first_basewc

    result = process.run_results()

    assert result == ExitCode(0)
    assert isinstance(process.outputs["stm_array"], orm.ArrayData)

