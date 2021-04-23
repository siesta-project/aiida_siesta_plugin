#!/usr/bin/env runaiida
import pytest
from plumpy import ProcessState
from aiida.engine import ExitCode
from aiida import orm
from aiida.common import (LinkType, AttributeDict)
from aiida_siesta.workflows.iterate import set_up_parameters_dict 

@pytest.fixture
def generate_workchain_converge(generate_psml_data, fixture_code, generate_workchain,
        generate_structure, generate_param, generate_basis, generate_kpoints_mesh,
        generate_calc_job_node, generate_parser):
    """Generate an instance of a `BandgapWorkChain`."""

    def _generate_workchain_converge():

        entry_point_wc = 'siesta.converger'

        psml = generate_psml_data('Si')

        inputs = {
            'code': fixture_code("siesta.siesta"),
            'structure': generate_structure(),
            'parameters': generate_param(),
            'pseudos': {
                'Si': psml,
                'SiDiff': psml
            },
            'options': orm.Dict(dict={
               'resources': {'num_machines': 1  },
               'max_wallclock_seconds': 1800,
               'withmpi': False,
               }),
            'iterate_over' : {"pao":[1,2],"mesh":[2,3]},
            'batch_size' : orm.Int(2)
        }

        process = generate_workchain(entry_point_wc, inputs)

        return process

    return _generate_workchain_converge


def test_converged(aiida_profile, generate_workchain_converge):
    """Test `SiestaConvereger` convergence check"""
    
    process = generate_workchain_converge()
    process.initialize()

    assert "target_values" in process.ctx

    proceed = process._should_proceed()
    assert proceed == True

    process.ctx.target_values.append(10)
    process.ctx.target_values.append(11)
    process.ctx.last_step_processes = []
    proceed = process._should_proceed()
    assert proceed == True

    process.ctx.target_values.append(11.001)
    proceed = process._should_proceed()
    assert proceed == False


def test_analyze_process(aiida_profile, generate_workchain_converge, 
        fixture_localhost, generate_wc_job_node):
    """ 
    Test method `_analyze_process` of `SiestaConvereger`. Moreover it
    calls return_results but we no correct setup. Therefore it does
    only test the case when the simulation is not converged and as output
    only the Bool `converged` is returned.
    """
    
    process = generate_workchain_converge()
    process.initialize()

    basewc = generate_wc_job_node("siesta.base", fixture_localhost)
    basewc.set_process_state(ProcessState.FINISHED)
    basewc.set_exit_status(ExitCode(0).status)
    out_par = orm.Dict(dict={"E_KS":1111,"E_KS_units":"eV"})
    out_par.store()
    out_par.add_incoming(basewc, link_type=LinkType.RETURN, link_label='output_parameters')
    
    assert len(process.ctx.target_values) == 0

    process._analyze_process(basewc)

    assert len(process.ctx.target_values) == 1

    process.ctx.used_values = []

    process.return_results()

    assert "converged" in process.outputs 
    assert process.outputs["converged"].value == False

@pytest.fixture
def generate_workchain_seq_converger(generate_workchain, generate_structure, generate_psml_data):

    psml=generate_psml_data("Si")
    def _generate_workchain_seq_converge():
        entry_point_wc = 'siesta.sequential_converger'
        inputs = {
            'converger_inputs' : {'structure': generate_structure(), 'pseudos': {'Si': psml, 'SiDiff': psml}},
            'iterate_over' : [{"pao":[1,2]},{"mesh":[2,3]}]
        }
        process = generate_workchain(entry_point_wc, inputs)
        return process
    
    return _generate_workchain_seq_converge

def test_sequential(aiida_profile, generate_workchain_seq_converger, generate_wc_job_node,
        fixture_localhost):
    """
    We test here the SiestaSequentialConverger, just the main two methods
    that distinguish it from the BaseIterator, meaning the `initialize` and the 
    `_analyze_process`
    """

    from aiida.common.extendeddicts import AttributeDict

    process = generate_workchain_seq_converger()
    process.initialize()

    assert process.ctx.iteration_keys == ('iterate_over',)

    convergerwc = generate_wc_job_node("siesta.converger", fixture_localhost)
    convergerwc.set_process_state(ProcessState.FINISHED)
    convergerwc.set_exit_status(ExitCode(0).status)
    out_par = orm.Dict(dict={"test_par_1":1111,"test_par_2":"eV"})
    out_par.store()
    out_conv = orm.Bool(True)
    out_conv.store()
    val_conv = orm.Float(1)
    val_conv.store()
    out_par.add_incoming(convergerwc, link_type=LinkType.RETURN, link_label='converged_parameters')
    out_conv.add_incoming(convergerwc, link_type=LinkType.RETURN, link_label='converged')
    val_conv.add_incoming(convergerwc, link_type=LinkType.RETURN, link_label='converged_target_value')

    process.ctx.last_inputs = AttributeDict({})

    process._analyze_process(convergerwc)

    assert process.ctx.already_converged == {"test_par_1":1111,"test_par_2":"eV"}
    assert process.ctx.last_target_value == val_conv
    assert "parameters" in process.ctx.last_inputs
    assert "test_par_1" in process.ctx.last_inputs.parameters.attributes
    assert process.ctx.last_inputs.parameters.attributes["test_par_2"] == "eV"



def test_sequential_not_conv(aiida_profile, generate_workchain_seq_converger, generate_wc_job_node,
        fixture_localhost):
    """
    We test here the SiestaSequentialConverger, in the case a Converger fails
    to converge.
    """

    from aiida.common.extendeddicts import AttributeDict

    process = generate_workchain_seq_converger()
    process.initialize()

    assert process.ctx.iteration_keys == ('iterate_over',)

    inputs = {"iterate_over": orm.Dict(dict={"s":[2,2]})}
    convergerwc = generate_wc_job_node("siesta.converger", fixture_localhost, inputs)
    convergerwc.set_process_state(ProcessState.FINISHED)
    convergerwc.set_exit_status(ExitCode(0).status)
    out_conv = orm.Bool(False)
    out_conv.store()
    out_conv.add_incoming(convergerwc, link_type=LinkType.RETURN, link_label='converged')

    process.ctx.last_inputs = AttributeDict({"parameters":{"yo":"yo"}})

    process._analyze_process(convergerwc)

    assert process.ctx.already_converged == {}
    assert "parameters" in process.ctx.last_inputs
