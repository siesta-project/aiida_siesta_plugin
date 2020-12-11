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

def test_pao_size(generate_ion_data):
    
    pao_man = PaoManager()

    ion = generate_ion_data('Si')

    pao_man.set_from_ion(ion)

    pao_man.pao_size

def test_change_all_radius(generate_ion_data):

    pao_man = PaoManager()

    ion = generate_ion_data('Si')

    pao_man.set_from_ion(ion)

    pao_man.change_all_radius(2)

def test_add_polarization(generate_ion_data):

    pao_man = PaoManager()

    ion = generate_ion_data('Si')

    pao_man.set_from_ion(ion)

    pao_man.add_polarization(3,1)


