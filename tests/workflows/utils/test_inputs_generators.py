#!/usr/bin/env runaiida
# -*- coding: utf-8 -*-
from aiida_siesta.groups.pseudos import PsmlFamily
from aiida.plugins import WorkflowFactory

def test_baseworkchain_inpgen(aiida_profile, fixture_code, generate_structure):
    """Test the validation of subclasses of `InputsGenerator`."""

    #Here I fake the pseudofamilies
    PsmlFamily.objects.get_or_create("nc-sr-04_pbe_standard_psml")
    PsmlFamily.objects.get_or_create("nc-sr-04_pbe_stringent_psml")

    from aiida_siesta.workflows.utils.inputs_generators import BaseWorkChainInputsGenerator

    inp_gen = BaseWorkChainInputsGenerator(WorkflowFactory("siesta.base"))
    structure = generate_structure()
    protocol = inp_gen.get_default_protocol_name()
    code = fixture_code("siesta.siesta")
    code.store()
    calc_engines = {"siesta": {'code': code.uuid, 'options': {"resources": {"num_mpiprocs_per_machine": 1}, "max_wallclock_seconds": 360}}} 

    build = inp_gen.get_filled_builder(structure, calc_engines, protocol)

    assert "parameters" in build
