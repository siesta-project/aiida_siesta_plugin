#!/usr/bin/env runaiida
# -*- coding: utf-8 -*-
import os
from aiida_siesta.workflows.utils import protocols #.protocols import ProtocolManager
from aiida_siesta.groups.pseudos import PsmlFamily
from importlib import reload 

def test_registries(aiida_profile):
    """
    Test that single calculation is submitted with the right content of the
    aiida.fdf file.
    """

    #Here I fake the pseudofamilies
    PsmlFamily.objects.get_or_create("nc-sr-04_pbe_standard_psml")
    PsmlFamily.objects.get_or_create("nc-sr-04_pbe_stringent_psml")

    basepath = os.path.dirname(os.path.abspath(__file__))
    filepath = os.path.join(basepath, 'fixtures/protocols/registries/custom_prot.yaml')

    os.environ["AIIDA_SIESTA_PROTOCOLS"] = filepath

    #I must reload protocols module because the "AIIDA_SIESTA_PROTOCOLS" is checked
    #in the class, not in its __init__ (when you instanciate). Is it a problem?
    reload(protocols)
    a=protocols.ProtocolManager()
    assert 'standard_my' in a.get_protocol_names()

    filepath = os.path.join(basepath, 'fixtures/protocols/registries/wrong_prot.yaml')
    os.environ["AIIDA_SIESTA_PROTOCOLS"] = filepath

    reload(protocols)

    import pytest
    with pytest.raises(RuntimeError):
        a=protocols.ProtocolManager()

    #assert 'standard_my' in a.get_protocol_names()


#def test_methods(aiida_profile):
    

