#!/usr/bin/env runaiida
import pytest
from plumpy import ProcessState
from aiida.engine import ExitCode
from aiida import orm
from aiida.common import (LinkType, AttributeDict)

@pytest.fixture
def generate_workchain_eos(generate_psml_data, fixture_code, fixture_localhost, generate_workchain, 
        generate_structure, generate_param, generate_basis, generate_kpoints_mesh,
        generate_calc_job_node, generate_parser):
    """Generate an instance of a `BandgapWorkChain`."""

    def _generate_workchain_eos():

        entry_point_wc = 'siesta.eos'
        entry_point_code = 'siesta.siesta'

        psml = generate_psml_data('Si')

        structure = generate_structure()

        inputs = {
            'code': fixture_code(entry_point_code),
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

    return _generate_workchain_eos


def test_setup_run(aiida_profile, generate_workchain_eos):
    """Test `EoSFixedCellShape`."""
    
    process = generate_workchain_eos()
    process.initio()

    assert isinstance(process.ctx.s0, orm.StructureData)

    process.run_base_wcs()


def test_results(aiida_profile, generate_workchain_eos,
        generate_wc_job_node, generate_structure, fixture_localhost):

    process = generate_workchain_eos()
    process.initio()

    values=[-13,-18,-21,-22,-21,-18,-13]

    for ind in range(len(process.ctx.scales)):
        inputs = AttributeDict({
            'structure': generate_structure(scale=process.ctx.scales[ind]),
            })
        basewc = generate_wc_job_node("siesta.base", fixture_localhost, inputs)
        basewc.set_process_state(ProcessState.FINISHED)
        basewc.set_exit_status(ExitCode(0).status)
        out_par = orm.Dict(dict={"E_KS":values[ind],"E_KS_units":"eV"})
        out_par.store()
        out_par.add_incoming(basewc, link_type=LinkType.RETURN, link_label='output_parameters')
        process.ctx[str(process.ctx.scales[ind])] = basewc

    result = process.return_results()

    assert result == ExitCode(0)
    assert isinstance(process.outputs["results_dict"], orm.Dict)
    assert (process.outputs["results_dict"]["fit_res"]['Vo(ang^3/atom)'] > 20)
