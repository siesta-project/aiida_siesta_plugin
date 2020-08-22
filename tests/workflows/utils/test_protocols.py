#!/usr/bin/env runaiida
# -*- coding: utf-8 -*-
import os
from aiida_siesta.workflows.utils.protocols import ProtocolManager
from aiida_siesta.groups.pseudos import PsmlFamily

def test_registries(aiida_profile):
    """
    Test that all the protocols are loaded, the one delivered in the package and
    the custom ones. It also loads a protocol with missing mandatory entries and it
    checks that an error is risen.
    """

    #Here I fake the pseudofamilies
    PsmlFamily.objects.get_or_create("nc-sr-04_pbe_standard-psf")
    PsmlFamily.objects.get_or_create("nc-sr-04_pbe_stringent-psf")

    basepath = os.path.dirname(os.path.abspath(__file__))
    filepath = os.path.join(basepath, 'fixtures/protocols/registries/custom_prot.yaml')

    os.environ["AIIDA_SIESTA_PROTOCOLS"] = filepath

    pmanager=ProtocolManager()
    assert 'standard_my' in pmanager.get_protocol_names()

    filepath = os.path.join(basepath, 'fixtures/protocols/registries/wrong_prot.yaml')
    os.environ["AIIDA_SIESTA_PROTOCOLS"] = filepath

    import pytest
    with pytest.raises(RuntimeError):
        pmanager=ProtocolManager()

    #A reset of the environment variable is needed at the end, otherwise it
    #messes up the other tests
    del os.environ['AIIDA_SIESTA_PROTOCOLS']


def test_methods(aiida_profile):
    """
    Test the 5 public methods of the class `ProtocolManager`
    """

    #Here I fake the pseudofamilies
    PsmlFamily.objects.get_or_create("nc-sr-04_pbe_standard-psf")
    PsmlFamily.objects.get_or_create("nc-sr-04_pbe_stringent-psf")

    pmanager=ProtocolManager()
        
    assert pmanager.is_valid_protocol("standard")
    assert not pmanager.is_valid_protocol("yoyo")

    reflist = ["standard","stringent"]
    plist = pmanager.get_protocol_names()
    assert plist == list(reflist)

    assert pmanager.get_default_protocol_name() == "standard"

    import pytest
    with pytest.raises(ValueError):
        pmanager.get_protocol_info("yoyo")
    assert pmanager.get_protocol_info("standard") == pmanager._protocols["standard"]["description"]

    import pytest
    with pytest.raises(ValueError):
        pmanager.get_protocol("yoyo")
    assert pmanager.get_protocol("standard") == pmanager._protocols["standard"]
