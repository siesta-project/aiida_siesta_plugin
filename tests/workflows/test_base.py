#!/usr/bin/env runaiida
import pytest
from plumpy import ProcessState
from aiida import orm
from aiida.common import LinkType
from aiida.common import AttributeDict
from aiida.engine import ProcessHandlerReport
from aiida_siesta.calculations.siesta import SiestaCalculation
from aiida_siesta.workflows.base import SiestaBaseWorkChain

@pytest.fixture
def generate_workchain_base(generate_psml_data, fixture_code, fixture_localhost, generate_workchain, 
        generate_structure, generate_param, generate_basis, generate_kpoints_mesh,
        generate_calc_job_node, generate_parser):
    """Generate an instance of a `SiestaBaseWorkChain`."""

    def _generate_workchain_base(exit_code=None):

        entry_point_wc = 'siesta.base'
        entry_point_code = 'siesta.siesta'

        psml = generate_psml_data('Si')

        inputs = {
            'code': fixture_code(entry_point_code),
            'structure': generate_structure(),
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

        if exit_code is not None:
            #Through the generate_calc_job_node we generate a node of a SiestaCalculation
            #with already remote_folder and retrieved as outputs (triggered by "default",
            #with `inputs2` as inputs and `attributes` as attributes
            inputs2 = AttributeDict({'structure': generate_structure()})
            attributes=AttributeDict({'input_filename':'aiida.fdf', 'output_filename':'aiida.out', 'prefix':'aiida'})
            node = generate_calc_job_node("siesta.siesta", fixture_localhost, "default", inputs2, attributes)
            #Set process state of the node
            node.set_process_state(ProcessState.FINISHED)
            node.set_exit_status(exit_code.status)

            process.ctx.iteration = 1
            process.ctx.children = [node]

        return process

    return _generate_workchain_base

def test_setup(aiida_profile, generate_workchain_base):
    """Test `SiestaBaseWorkChain.setup`."""
    process = generate_workchain_base()
    process.setup()
    process.prepare_inputs()

    assert isinstance(process.ctx.inputs, dict)

    
def test_handle_error_geom_not_conv(aiida_profile, generate_workchain_base):
    """Test `SiestaBaseWorkChain.handle_error_geom_not_conv`."""
    process = generate_workchain_base(exit_code=SiestaCalculation.exit_codes.GEOM_NOT_CONV)
    process.setup()
    process.prepare_inputs()

    #Add another fake output to the SiestaCalculation node
    calculation = process.ctx.children[-1]
    out_par = orm.Dict(dict={"variable_geometry":False})
    out_par.add_incoming(calculation, link_type=LinkType.CREATE, link_label='output_parameters')
    out_par.store()

    result = process.handle_error_geom_not_conv(calculation)
    assert isinstance(result, ProcessHandlerReport)
    assert result.do_break
    #assert result.exit_code == SiestaBaseWorkChain.exit_codes.GEOM_NOT_CONV

    #result = process.inspect_process()
    #assert result == PwBaseWorkChain.exit_codes.ERROR_UNRECOVERABLE_FAILURE

def test_handle_error_scf_not_conv(aiida_profile, generate_workchain_base):
    """Test `SiestaBaseWorkChain.handle_error_scf_not_conv`."""
    process = generate_workchain_base(exit_code=SiestaCalculation.exit_codes.SCF_NOT_CONV)
    process.setup()
    process.prepare_inputs()

    #Add another fake output to the SiestaCalculation node
    calculation = process.ctx.children[-1]
    out_par = orm.Dict(dict={"variable_geometry":False})
    out_par.add_incoming(calculation, link_type=LinkType.CREATE, link_label='output_parameters')
    out_par.store()

    result = process.handle_error_scf_not_conv(calculation)
    assert isinstance(result, ProcessHandlerReport)
    assert result.do_break
    #assert result.exit_code == SiestaBaseWorkChain.exit_codes.GEOM_NOT_CONV

    #result = process.inspect_process()
    #assert result == PwBaseWorkChain.exit_codes.ERROR_UNRECOVERABLE_FAILURE


def test_handle_error_basis_pol(aiida_profile, generate_workchain_base):
    """Test `SiestaBaseWorkChain.handle_error_basis_pol`."""
    process = generate_workchain_base(exit_code=SiestaCalculation.exit_codes.BASIS_POLARIZ)
    process.setup()
    process.prepare_inputs()

    result = process.handle_error_basis_pol(process.ctx.children[-1])
    assert isinstance(result, ProcessHandlerReport)
    assert result.do_break
    assert result.exit_code == SiestaBaseWorkChain.exit_codes.ERROR_BASIS_POL


def test_handle_error_bands(aiida_profile, generate_workchain_base):
    """Test `SiestaBaseWorkChain.handle_error_bands`."""
    process = generate_workchain_base(exit_code=SiestaCalculation.exit_codes.BANDS_PARSE_FAIL)
    process.setup()
    process.prepare_inputs()

    result = process.handle_error_bands(process.ctx.children[-1])
    assert isinstance(result, ProcessHandlerReport)
    assert result.do_break
    assert result.exit_code == SiestaBaseWorkChain.exit_codes.ERROR_BANDS_PARSING

    #result = process.inspect_process()
    #assert result == PwBaseWorkChain.exit_codes.ERROR_UNRECOVERABLE_FAILURE


def test_handle_error_split_norm(aiida_profile, generate_workchain_base):
    """Test `SiestaBaseWorkChain.handle_error_split_norm`."""
    process = generate_workchain_base(exit_code=SiestaCalculation.exit_codes.SPLIT_NORM)
    process.setup()
    process.prepare_inputs()

    #Add log info at the SiestaCalculation
    calculation = process.ctx.children[-1]
    calculation.logger.error("Error in split_norm option. Minimum value is 1")

    result = process.handle_error_split_norm(calculation)
    assert isinstance(result, ProcessHandlerReport)
    assert result.do_break

