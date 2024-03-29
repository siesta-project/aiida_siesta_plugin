# -*- coding: utf-8 -*-
import pytest

from aiida_siesta.utils.pao_manager import PaoManager


def test_set_from_ion(generate_ion_data):

    pao_man = PaoManager()

    ion = generate_ion_data('Si')

    pao_man.set_from_ion(ion)

    assert pao_man.name == "Si"
    assert pao_man._gen_dict is not None
    assert pao_man._pol_dict == {3: {1: {1: 4.0531999999999995, 2: 3.1566}}}
    assert pao_man._conf_dict == {}

def test_validator_and_get_pao_block():

    pao_man = PaoManager()

    with pytest.raises(RuntimeError):
        pao_man.get_pao_block()

    pao_man.name = "Si"

    with pytest.raises(RuntimeError):
        pao_man.get_pao_block()

    pao_man._gen_dict = {3: {0: {1: 4.05}}}

    with pytest.raises(RuntimeError):
        pao_man.get_pao_block()

    pao_man._pol_dict = {}

    assert pao_man.get_pao_block() == "Si 1\n  n=3  0  1\n 7.65335"

    pao_man._gen_dict = {}

    with pytest.raises(RuntimeError):
        pao_man.get_pao_block()


def test_confinements_features(generate_ion_data):

    pao_man = PaoManager()

    ion = generate_ion_data('Si_with_conf')

    pao_man.set_from_ion(ion)

    assert pao_man.name == "Si"
    assert pao_man._gen_dict is not None
    assert pao_man._pol_dict == {3: {1: {1: 4.0531999999999995, 2: 3.1566}}}
    assert pao_man._conf_dict == {'Q': {3: {1: [3.0, 0.5, 0.01]}}, 'E': {3: {0: [2.0, 0.3]}}}

    assert pao_man.get_pao_block() == 'Si 2\n  n=3  0  2 E 2.0 0.3 \n 5.965078\t 4.419101\n  n=3  1  2 P 2 Q 3.0 0.5 0.01 \n 7.659398\t 5.13417'

def test_pao_size(generate_ion_data):

    pao_man = PaoManager()

    ion = generate_ion_data('Si')

    pao_man.set_from_ion(ion)

    assert pao_man.pao_size() == "DZDP"

def test_change_all_radius():

    pao_man = PaoManager()

    pao_man.name = "Si"
    pao_man._gen_dict = {3: {0: {1: 4.05}}}
    pao_man._pol_dict = {3: {0: {1: 4.05}}}

    pao_man.change_all_radius(2)

    assert pao_man._gen_dict == {3: {0: {1: 4.131}}}
    assert pao_man._pol_dict == {3: {0: {1: 4.131}}}


def test_reset_radius():

    pao_man = PaoManager()

    pao_man.name = "Si"
    pao_man._gen_dict = {3: {0: {1: 4.05}}}
    pao_man._pol_dict = {3: {0: {1: 4.05}}}

    with pytest.raises(ValueError):
        pao_man.reset_radius("Bohr",0.0,3,1,2)

    pao_man.reset_radius("Bohr",0.0,3,0,1)
    assert pao_man._gen_dict == {3: {0: {1: 0.0}}}
    assert pao_man._pol_dict == {3: {0: {1: 0.0}}}


def test_add_polarization():

    pao_man = PaoManager()

    pao_man.name = "Si"
    pao_man._gen_dict = {3: {0: {1: 4.05}}}
    pao_man._pol_dict = {3: {0: {1: 4.05}}}

    with pytest.raises(ValueError):
        pao_man.add_polarization(3,1)

    pao_man.add_polarization(3,0)

    assert pao_man._pol_dict == {3: {0: {1: 4.05, 2: 0.0}}}
    assert pao_man.pao_size() == "SZDP"


def test_remove_polarization():

    pao_man = PaoManager()

    pao_man.name = "Si"
    pao_man._gen_dict = {3: {0: {1: 4.05}}}
    pao_man._pol_dict = {3: {0: {1: 4.05, 2: 0.0}}}

    with pytest.raises(ValueError):
        pao_man.remove_polarization(3,1)

    pao_man.remove_polarization(3,0)
    assert pao_man._pol_dict == {3: {0: {1: 4.05}}}
    assert pao_man.pao_size() == "SZP"

    pao_man.remove_polarization(3,0)
    assert pao_man._pol_dict == {}
    assert pao_man.pao_size() == "SZ"


def test_add_orbital():

    pao_man = PaoManager()

    pao_man.name = "Si"
    pao_man._gen_dict = {3: {0: {1: 4.05}}}
    pao_man._pol_dict = {3: {0: {1: 4.05}}}

    with pytest.raises(ValueError):
        pao_man.add_orbital("Bohr",0.0,3,1,2)

    pao_man.add_orbital("Bohr",0.0,3,0,2)

    assert pao_man._gen_dict == {3: {0: {1: 4.05, 2: 0.0}}}
    assert pao_man.pao_size() == "DZP"


def test_remove_orbital():

    pao_man = PaoManager()

    pao_man.name = "Si"
    pao_man._gen_dict = {3: {0: {1: 4.05, 2: 0.0}}}
    pao_man._pol_dict = {3: {0: {1: 4.05}}}

    with pytest.raises(ValueError):
        pao_man.remove_orbital(3,1,1)
    with pytest.raises(ValueError):
        pao_man.remove_orbital(3,0,1)

    pao_man.remove_orbital(3,0,2)
    assert pao_man._gen_dict == {3: {0: {1: 4.05}}}
    assert pao_man._pol_dict == {3: {0: {1: 4.05}}}
    assert pao_man.pao_size() == "SZP"

    pao_man.remove_orbital(3,0,1)
    assert pao_man._gen_dict == {}
    assert pao_man._pol_dict == {}


def test_remove_polarization_occu(generate_ion_data):

    pao_man = PaoManager()

    ion = generate_ion_data('Si_with_conf')

    pao_man.set_from_ion(ion)

    assert pao_man._pol_occu == {3: {1: {1: 0.0, 2: 0.0}}}

    pao_man.remove_polarization(3,1)
    assert pao_man._pol_dict == {3: {1: {1: 4.0531999999999995}}}
    assert pao_man._pol_occu == {3: {1: {1: 0.0}}}

    pao_man.remove_polarization(3,1)
    assert pao_man._pol_dict == {}
    assert pao_man._pol_occu == {}


def test_remove_orbital_occu_and_conf(generate_ion_data):

    pao_man = PaoManager()

    ion = generate_ion_data('Si_with_conf')

    pao_man.set_from_ion(ion)

    assert pao_man._gen_occu == {3: {0: {1: 2.0, 2: 0.0}, 1: {1: 2.0, 2: 0.0}}}
    assert pao_man._pol_occu == {3: {1: {1: 0.0, 2: 0.0}}}
    assert pao_man._conf_dict == {'E': {3: {0: [2.0, 0.3]}}, 'Q': {3: {1: [3.0, 0.5, 0.01]}}}

    pao_man.remove_orbital(3,0,2)
    assert pao_man._gen_occu == {3: {0: {1: 2.0}, 1: {1: 2.0, 2: 0.0}}}
    assert pao_man._pol_occu == {3: {1: {1: 0.0, 2: 0.0}}}
    assert pao_man._conf_dict == {'E': {3: {0: [2.0, 0.3]}}, 'Q': {3: {1: [3.0, 0.5, 0.01]}}}

    pao_man.remove_orbital(3,0,1)
    assert pao_man._gen_occu == {3: {1: {1: 2.0, 2: 0.0}}}
    assert pao_man._pol_occu == {3: {1: {1: 0.0, 2: 0.0}}}
    assert pao_man._conf_dict == {'Q': {3: {1: [3.0, 0.5, 0.01]}}}

    pao_man.remove_orbital(3,1,2)
    assert pao_man._gen_occu == {3: {1: {1: 2.0}}}
    assert pao_man._pol_occu == {3: {1: {1: 0.0, 2: 0.0}}}
    assert pao_man._conf_dict == {'Q': {3: {1: [3.0, 0.5, 0.01]}}}

    pao_man.remove_orbital(3,1,1)
    assert pao_man._gen_occu == {}
    assert pao_man._pol_occu == {}
    assert pao_man._conf_dict == {}
