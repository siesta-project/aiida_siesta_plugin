#!/usr/bin/env runaiida
import pytest
from aiida import orm


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


def test_setup(aiida_profile, generate_workchain_eos):
    """Test `EoSFixedCellShape`."""
    
    process = generate_workchain_eos()
    process.initio()

    assert isinstance(process.ctx.s0, orm.StructureData)
