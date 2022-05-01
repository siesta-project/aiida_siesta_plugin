#!/usr/bin/env runaiida
# -*- coding: utf-8 -*-
from aiida import orm
from aiida.common import AttributeDict, LinkType
from aiida.engine import ExitCode
from plumpy import ProcessState
import pytest

from aiida_siesta.workflows.iterate import set_up_parameters_dict


@pytest.fixture
def generate_workchain_iterate(generate_psml_data, fixture_code, generate_workchain,
        generate_structure, generate_param, generate_basis, generate_kpoints_mesh,
        generate_calc_job_node, generate_parser):
    """Generate an instance of a `BandgapWorkChain`."""

    def _generate_workchain_iterate():

        entry_point_wc = 'siesta.iterator'

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

    return _generate_workchain_iterate


def test_setup_next_step_smooth_run(aiida_profile, generate_workchain_iterate):
    """Test `SiestaBaseWorkChain.setup`."""
    process = generate_workchain_iterate()
    process.initialize()

    assert process.ctx.iteration_keys == ("pao","mesh")
    a_val_list = process.ctx.iteration_vals[0]
    assert isinstance(a_val_list,list)
    assert isinstance(orm.load_node(a_val_list[0]), orm.Int)
    assert "mesh" in process.ctx._iteration_parsing
    assert 'input_key' in process.ctx._iteration_parsing["mesh"]
    assert process.ctx._iteration_parsing["mesh"]['input_key'] == "parameters"
    assert process.ctx._iteration_parsing["pao"]['input_key'] == "basis"
    #assert process.ctx._iteration_parsing["mesh"]['parse_func'] == set_up_parameters_dict
    assert isinstance(process.ctx.values_iterator, zip)
    assert isinstance(process.ctx.inputs, dict)

    assert len(process.ctx.used_values) == 0
    process.next_step()
    assert len(process.ctx.used_values) == 1

    process.run_batch()

    assert len(process.ctx.used_values) == 2
    assert "structure" in process.ctx.last_inputs
    assert "pao" in process.ctx.last_inputs["basis"].attributes
