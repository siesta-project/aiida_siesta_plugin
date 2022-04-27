#!/usr/bin/env runaiida
import pytest
from plumpy import ProcessState
from aiida.engine import ExitCode
from aiida import orm
from aiida.common import (LinkType, AttributeDict)
from aiida_siesta.workflows.iterate import set_up_parameters_dict 

@pytest.fixture
def generate_workchain_two_steps_op(generate_psml_data, fixture_code, generate_workchain,
        generate_structure, generate_param, generate_basis, generate_kpoints_mesh,
        generate_calc_job_node, generate_parser):
    """Generate an instance of a `BandgapWorkChain`."""

    def _generate_workchain_two_steps_op():

        entry_point_wc = 'siesta.two_step_basis_opt'

        psml = generate_psml_data('Si')

        basis = orm.Dict(
            dict={
                '%block pao-basis': "\nSi   2\n n=3   0   2\n 4.99376   $sz2 \n n=3   1   2\n 6.2538    $pz2 \n%endblock pao-basis"
            })


        inputs = {
            'siesta_base': {
                'structure': generate_structure(),
                'parameters': generate_param(),
                'code': fixture_code("siesta.siesta"),
                'basis': basis,
                'pseudos': {'Si': psml, 'SiDiff': psml},
                'options': orm.Dict(dict={
                    'resources': {'num_machines': 1  },
                    'max_wallclock_seconds': 1800,
                    'withmpi': False,
                    }),
                },
            'simplex': {
                'output_name': orm.Str("basis_enthalpy"),
                'max_iters': orm.Int(4),
                'variables_dict': orm.Dict(dict={
                    "sz2":[2.0,4.8,3.0],
                    "pz2":[2.0,6.0,3.0]
                    }),
                },
            'macrostep':{
                'lambda_scaling_factor': orm.Float(0.5),
                'initial_lambda': orm.Float(0.4),
                'minimum_lambda': orm.Float(0.03)
                },
        }

        process = generate_workchain(entry_point_wc, inputs)

        return process

    return _generate_workchain_two_steps_op


def test_setup_next_step_smooth_run(aiida_profile, generate_workchain_two_steps_op):
    """Test `SiestaBaseWorkChain.setup`."""
    process = generate_workchain_two_steps_op()
    process.preprocess()

    assert process.ctx.current_lambda == 0.4
    assert process.ctx.vars_dict == {'pz2': [2.0, 6.0, 3.0], 'sz2': [2.0, 4.8, 3.0]}

    assert process.should_run() == True

    process.run_simplex()
