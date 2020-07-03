#!/usr/bin/env runaiida
import pytest
from aiida import orm


@pytest.fixture
def generate_workchain_stm(generate_psml_data, fixture_code, fixture_localhost, generate_workchain, 
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


def test_setup(aiida_profile, generate_workchain_stm):
    """Test `SiestaSTMWorkChain`."""
    
    process = generate_workchain_stm()
    process.checks()

#    assert isinstance(process.ctx.s0, orm.StructureData)

