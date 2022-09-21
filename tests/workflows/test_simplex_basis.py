#!/usr/bin/env runaiida
# -*- coding: utf-8 -*-
from aiida import orm
from aiida.common import AttributeDict, LinkType
from aiida.engine import ExitCode
from plumpy import ProcessState
import pytest

from aiida_siesta.workflows.iterate import set_up_parameters_dict


@pytest.fixture
def generate_workchain_simplex_basis(generate_psml_data, fixture_code, generate_workchain,
        generate_structure, generate_param, generate_basis, generate_kpoints_mesh,
        generate_calc_job_node, generate_parser):
    """Generate an instance of a `BandgapWorkChain`."""

    def _generate_workchain_simplex_basis():

        entry_point_wc = 'siesta.simplex_basis'

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
        }

        process = generate_workchain(entry_point_wc, inputs)

        return process

    return _generate_workchain_simplex_basis


def test_setup_next_step_smooth_run(aiida_profile, generate_workchain_simplex_basis):
    """Test `SiestaBaseWorkChain.setup`."""
    process = generate_workchain_simplex_basis()
    process.preprocess()

    assert process.ctx.simplex[0] == [3.0, 3.0]

    process.run_optimizer()
