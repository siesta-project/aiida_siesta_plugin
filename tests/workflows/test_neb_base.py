#!/usr/bin/env runaiida
# -*- coding: utf-8 -*-
from aiida import orm
from aiida.common import AttributeDict, LinkType
from aiida.engine import ProcessHandlerReport
from aiida_pseudo.groups.family.pseudo import PseudoPotentialFamily
from plumpy import ProcessState
import pytest

from aiida_siesta.calculations.siesta import SiestaCalculation
from aiida_siesta.workflows.base import SiestaBaseWorkChain

#from aiida_siesta.groups.pseudos import PsmlFamily


@pytest.fixture
def generate_workchain_neb_base(generate_path_object, fixture_code, fixture_localhost, generate_workchain,
        generate_structure, generate_param, generate_psf_data, generate_kpoints_mesh, generate_lua_file,
        generate_calc_job_node, generate_parser):
    """
    Generate an instance of a `SiestaBaseNEBWorkChain`.
    The parameters:
    1) exit_code. This inputs triggers the generation of a SiestaCalculation through the
       `generate_calc_job_node` fixture with a particular exit_code and it is attached to
       the current SiestaBaseWorkChain as children.
    2) num_files. An input to determine the number of xyz files trasformed in structures of
       a TrajectoryData and passed in `starting_path` input. Useful to test the validator.
    3) kinds. Bool that triggesrs to attach the kinds as an attribute of the TrajectoryData
       in `starting_path` input. Again used to test validator.
    """

    def _generate_workchain_neb_base(exit_code=None, num_files=None, kinds=True):

        entry_point_wc = 'siesta.baseneb'
        entry_point_code = 'siesta.siesta'

        inputs = {
            'code': fixture_code(entry_point_code),
            'kpoints': generate_kpoints_mesh(2),
            'parameters': generate_param(),
            'pseudos': {
                'O': generate_psf_data('O'),
                'H': generate_psf_data('H'),
            },
            'neb_script': generate_lua_file(),
            'starting_path': generate_path_object(num_files, kinds),
            'options': orm.Dict(dict={
               'resources': {'num_machines': 1  },
               'max_wallclock_seconds': 1800,
               'withmpi': False,
               })
        }

        process = generate_workchain(entry_point_wc, inputs)

        if exit_code is not None:
            # Generate a SiestaCalculation node with exit code. The node comes
            # with `remote_folder` and `retrieved` as outputs (triggered by "neb",
            # which is a folder in `parsers/fixtures/siesta`)
            node = generate_calc_job_node("siesta.siesta", fixture_localhost, "neb")
            #Set process state of the node
            node.set_process_state(ProcessState.FINISHED)
            node.set_exit_status(exit_code.status)

            process.ctx.neb_wk = node

        return process

    return _generate_workchain_neb_base

def test_validators(aiida_profile, generate_workchain_neb_base):
    """
    Test SiestaBaseNEBWorkChain validator of the input port `starting_path`.
    """
    with pytest.raises(ValueError):
        process = generate_workchain_neb_base(num_files=1, kinds=True)
    with pytest.raises(ValueError):
        process = generate_workchain_neb_base(num_files=2, kinds=True)
    with pytest.raises(ValueError):
        process = generate_workchain_neb_base(num_files=None, kinds=False)

def test_create_reference_structure_and_run(aiida_profile, generate_workchain_neb_base):
    """
    Test `SiestaBaseNEBWorkChain.create_reference_structure` and
    `SiestaBaseNEBWorkChain.run_neb` functions.
    """
    process = generate_workchain_neb_base(num_files=None, kinds=True)
    process.create_reference_structure()

    assert isinstance(process.ctx.reference_structure, orm.StructureData)

    process.run_neb()

def test_run_results(aiida_profile, generate_workchain_neb_base):
    """
    Test `SiestaBaseNEBWorkChain.run_results` function.
    """
    from aiida.engine import ExitCode

    process = generate_workchain_neb_base(exit_code=ExitCode(0))
    process.create_reference_structure()
    process.run_results()

    assert ("neb_output_package" in process.outputs)
