#!/usr/bin/env runaiida
# -*- coding: utf-8 -*-
from aiida import orm
from aiida.common import AttributeDict, LinkType
from aiida.engine import ExitCode
from plumpy import ProcessState
import pytest

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

        params = generate_param().get_dict()
        params["%block local-density-of-states"] =  "\n -9.6  -1.6 eV \n %endblock local-density-of-states"

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
            'parameters': orm.Dict(dict=params),
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


def test_setup_mainrun_errorbasewc(aiida_profile, generate_workchain_stm,
        generate_wc_job_node, fixture_localhost):
    """Test `SiestaSTMWorkChain`."""

    process = generate_workchain_stm()
    process.checks()

    assert process.ctx.spinstm == "none"
    assert process.ctx.ldosdefinedinparam

    process.run_siesta_wc()

    #Here we could check that process.run_siesta_wc() returns ToContext(workchain_base=running)

    #Fake the base siesta wc in context
    first_basewc = generate_wc_job_node("siesta.base", fixture_localhost)
    first_basewc.set_process_state(ProcessState.FINISHED)
    #Just like this, the exit status in not zero, therefore not finished_ok
    #It will allow to check error in next step
    process.ctx.workchain_base = first_basewc

    result = process.run_siesta_with_ldos()

    assert result == SiestaSTMWorkChain.exit_codes.ERROR_BASE_WC

def test_runldos_errorldos(aiida_profile, generate_workchain_stm, generate_psml_data,
        generate_wc_job_node, fixture_localhost, fixture_code, generate_structure):

    process = generate_workchain_stm()
    process.checks()

    #Fake the base siesta wc in context, now with some inputs
    psml = generate_psml_data('Si')
    inputs = AttributeDict({
        'structure': generate_structure(),
        'code': fixture_code("siesta.siesta"),
        'parameters': orm.Dict(dict={"sm":"sm"}),
        'options': orm.Dict(dict={'resources': {'num_machines': 1  },'max_wallclock_seconds': 1800,'withmpi': False}),
        'pseudos': {'Si': psml,'SiDiff': psml},
        })
    first_basewc = generate_wc_job_node("siesta.base", fixture_localhost, inputs)
    first_basewc.set_process_state(ProcessState.FINISHED)
    first_basewc.set_exit_status(ExitCode(0).status)
    #Now is_finished_ok, next step will submit ldos calc
    #Also needed to fake outputs ports of the basewc: remote_folder and output_parameters
    #It is different respect to a CalcJob, nodes has to be stored before and link is "RETURN"
    remote_folder = orm.RemoteData(computer=fixture_localhost, remote_path='/tmp')
    remote_folder.store()
    remote_folder.add_incoming(first_basewc, link_type=LinkType.RETURN, link_label='remote_folder')
    out_par = orm.Dict(dict={"E_Fermi":-1})
    out_par.store()
    out_par.add_incoming(first_basewc, link_type=LinkType.RETURN, link_label='output_parameters')
    process.ctx.workchain_base = first_basewc

    process.run_siesta_with_ldos()

    #Here we might check ToContext(siesta_ldos=running)

    ldos_basewc = generate_wc_job_node("siesta.base", fixture_localhost)
    ldos_basewc.set_process_state(ProcessState.FINISHED)
    #We don't set exit status so that if appears as not is_finished_ok
    process.ctx.siesta_ldos = ldos_basewc

    result = process.run_stm()

    assert result == SiestaSTMWorkChain.exit_codes.ERROR_LDOS_WC

def test_runstm_failstm(aiida_profile, generate_workchain_stm, generate_wc_job_node,
         generate_calc_job_node, fixture_localhost):

    process = generate_workchain_stm()
    process.checks()

    ldos_basewc = generate_wc_job_node("siesta.base", fixture_localhost)
    ldos_basewc.set_process_state(ProcessState.FINISHED)
    ldos_basewc.set_exit_status(ExitCode(0).status)
    #Now is_finished_ok, but need to set outputs
    remote_folder = orm.RemoteData(computer=fixture_localhost, remote_path='/tmp')
    remote_folder.store()
    remote_folder.add_incoming(ldos_basewc, link_type=LinkType.RETURN, link_label='remote_folder')
    process.ctx.siesta_ldos = ldos_basewc

    process.run_stm()

    #Fake the stm calculation
    name = 'default'
    entry_point_calc_job = 'siesta.stm'
    inputs = AttributeDict({
        'spin_option' : orm.Str("q")
    })
    attributes=AttributeDict({'input_filename':'stm.in', 'output_filename':'stm.out'})
    stm_node = generate_calc_job_node(entry_point_calc_job, fixture_localhost, name, inputs, attributes)
    stm_node.set_process_state(ProcessState.FINISHED)
    process.ctx.stm_calc = stm_node

    result = process.run_results()

    assert result == SiestaSTMWorkChain.exit_codes.ERROR_STM_PLUGIN


def test_outputs(aiida_profile, generate_workchain_stm, generate_wc_job_node,
         generate_calc_job_node, fixture_localhost):
    """Test `SiestaSTMWorkChain`."""

    process = generate_workchain_stm()
    process.checks()

    name = 'default'
    entry_point_calc_job = 'siesta.stm'
    inputs = AttributeDict({
        'spin_option' : orm.Str("q")
    })
    attributes=AttributeDict({'input_filename':'stm.in', 'output_filename':'stm.out'})
    stm_node = generate_calc_job_node(entry_point_calc_job, fixture_localhost, name, inputs, attributes)
    stm_node.set_process_state(ProcessState.FINISHED)
    stm_node.set_exit_status(ExitCode(0).status)
    stm_array = orm.ArrayData()
    stm_array.add_incoming(stm_node, link_type=LinkType.CREATE, link_label='stm_array')
    stm_array.store()
    process.ctx.stm_calc = stm_node

    first_basewc = generate_wc_job_node("siesta.base", fixture_localhost)
    out_par = orm.Dict(dict={"variable_geometry":False})
    out_par.store()
    out_par.add_incoming(first_basewc, link_type=LinkType.RETURN, link_label='output_parameters')
    process.ctx.workchain_base = first_basewc

    result = process.run_results()

    assert result == ExitCode(0)
    assert isinstance(process.outputs["stm_array"], orm.ArrayData)
