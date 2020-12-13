import pytest
from aiida_siesta.utils.pao_manager import PaoManager

def test_set_from_ion(generate_ion_data):

    pao_man = PaoManager()

    ion = generate_ion_data('Si')

    pao_man.set_from_ion(ion)

    assert pao_man.name == "Si"
    assert pao_man._gen_dict is not None
    assert pao_man._pol_dict == {3: {1: {1: 4.0531999999999995, 2: 3.1566}}}

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

    assert pao_man.get_pao_block() == "Si 1\n  n=3  0  1 \n 7.6533504667600045"

    pao_man._gen_dict = {}

    with pytest.raises(RuntimeError):
        pao_man.get_pao_block()


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


