#!/usr/bin/env runaiida
import pytest
from plumpy import ProcessState
from aiida import orm
from aiida.common import (LinkType, AttributeDict)
from aiida.engine import ExitCode

@pytest.fixture
def generate_workchain_epsilon(generate_psml_data, fixture_code, fixture_localhost, generate_workchain, 
        generate_structure, generate_param, generate_basis, generate_kpoints_mesh,
        generate_calc_job_node, generate_parser):
    """Generate an instance of a `EpsilonWorkChain`."""

    def _generate_workchain_bandgap(optical=False):

        entry_point_wc = 'siesta.epsilon'
        entry_point_code = 'siesta.siesta'

        psml = generate_psml_data('Si')

        inputs = {
            'code': fixture_code(entry_point_code),
            'structure': generate_structure(),
            'parameters': generate_param(),
            'kpoints': generate_kpoints_mesh(2),
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

        if optical:
            optical = {
                "optical-calculation" : True,
                "optical-broaden": "0.5 eV",
                "optical-polarization-type": "polarized",
                "%block optical-vector": "\n 1.0 0.0 0.0 \n%endblock optical-vector",
                "%block optical-mesh": "\n 3 3 3 \n%endblock optical-mesh"
            }
            inputs["optical"] = orm.Dict(dict=optical)


        process = generate_workchain(entry_point_wc, inputs)

        return process

    return _generate_workchain_bandgap


def test_opt_inputs(aiida_profile, generate_workchain_epsilon):
    """Test `BangapWorkChain.preprocess`."""

    with pytest.raises(ValueError):
        process = generate_workchain_epsilon()

    process = generate_workchain_epsilon(optical=True)


#I do not know how to test the postprocess function.
