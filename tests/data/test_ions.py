import pytest
from aiida.common.exceptions import StoringNotAllowed

def test_ions(generate_ion_data):
    """
    Test the `set_file`, `_prepare_source` and `parse_ion` methods of IonData. They are called
    during instanciation of the class. Instanciation can be done through the passing of a stream of 
    a filepath. The parent class already allows that but, it is important to test the parsing 
    of the info on the file (element, name, atomic_number, mass) in both cases since if filepath 
    is passed, it is immediately converted to a stream.
    """

    ion = generate_ion_data('Si')
    assert ion.element == "Si"
    assert ion.name == "Si"
    assert ion.atomic_number == 14
    assert ion.mass == 28.09
    assert ion.is_stored == False
    ion.store()
    assert ion.is_stored == True

    ion = generate_ion_data('Si', stream=True)
    assert ion.element == "Si"
    assert ion.name == "Si"
    assert ion.atomic_number == 14
    assert ion.mass == 28.09
    assert ion.is_stored == False
    ion.store()
    assert ion.is_stored == True


def test_ions_validators(generate_ion_data):
    """
    Test the validators called before the storing of the IonData
    """

    ion = generate_ion_data('Si')

    ion.attributes["element"] = "new_name"
    assert ion.is_stored == False
    assert ion.element == "new_name"
    with pytest.raises(StoringNotAllowed):
        ion.store()

    ion = generate_ion_data('Si')

    ion.attributes["md5"] = "weee"
    assert ion.is_stored == False
    assert ion.md5 == "weee"
    with pytest.raises(StoringNotAllowed):
        ion.store()

    ion = generate_ion_data('Si')

    ion.attributes["atomic_number"] = 22
    assert ion.is_stored == False
    assert ion.atomic_number == 22
    with pytest.raises(StoringNotAllowed):
        ion.store()

    ion = generate_ion_data('Si')

    ion.attributes["name"] = "namee"
    assert ion.is_stored == False
    assert ion.name == "namee"
    with pytest.raises(StoringNotAllowed):
        ion.store()


def test_get_or_create(generate_ion_data):
    """
    The `get_or_create` method is tested. It is very important to test
    because it also needs to replicate the checks the instantiation does
    and it can not rely on the parent class
    """

    from aiida_siesta.data.ion import IonData
    import os,io
    filename = os.path.join('tests', 'fixtures', 'ions', 'SiDiff.ion.xml')
    filepath = os.path.abspath(filename)

    # Test that strings must be abs paths
    with pytest.raises(TypeError):
        IonData.get_or_create("ss")

    # Test that other apart from sting or stream is not allowed
    with pytest.raises(TypeError):
        IonData.get_or_create(1)

    # Test that strims of bytes only are accepted.
    with pytest.raises(TypeError):
        with io.open(filepath, 'r') as handle:
            IonData.get_or_create(handle)

    # Good stream
    with io.open(filepath, 'rb') as handle:
        test_i = IonData.get_or_create(handle)
    assert test_i.is_stored == False

    ion = IonData.get_or_create(filepath)
    assert ion.is_stored == False
    ion.store()
    assert ion.is_stored == True

    #Now it is stored! Therefore a new call will say that is stored already
    ion = IonData.get_or_create(filepath)
    assert ion.is_stored == True

    with io.open(filepath, 'rb') as handle:
        test_i = IonData.get_or_create(handle)
    assert test_i.is_stored == True


def test_get_orbitals(generate_ion_data):
    """
    Test the get_orbitals method
    """

    from aiida_siesta.data.atomic_orbitals import SislAtomicOrbital
    ion = generate_ion_data('Si')
    orbit_list = ion.get_orbitals()
    assert len(orbit_list) == 18
    assert isinstance(orbit_list[0],SislAtomicOrbital)


def test_analyze_basis_specs(generate_ion_data):
    """
    Test the method hidden _analyze_basis_specs
    """

    ion = generate_ion_data('Si')
    dict_conf = ion._analyze_basis_specs()
    assert dict_conf == {'3p': {'E': [0.0, 0.0], 'Q': [0.0, 0.0, 0.01]}, '3s': {'E': [0.0, 0.0], 'Q': [0.0, 0.0, 0.01]}}

    ion = generate_ion_data('Si_with_conf')
    dict_conf = ion._analyze_basis_specs()
    assert dict_conf == {'3p': {'E': [0.0, 0.0], 'Q': [3.0, 0.5, 0.01]}, '3s': {'E': [2.0, 0.3], 'Q': [0.0, 0.0, 0.01]}}


def test_get_info_charge_conf(generate_ion_data):
    """
    Test the method get_info_charge_confinement
    """

    ion = generate_ion_data('Si')
    dict_conf = ion.get_info_charge_confinement()
    assert dict_conf == {}

    ion = generate_ion_data('Si_with_conf')
    dict_conf = ion.get_info_charge_confinement()
    assert dict_conf == {'3p': [3.0, 0.5, 0.01]}


def test_get_info_soft_conf(generate_ion_data):
    """
    Test the method get_info_soft_confinement
    """

    ion = generate_ion_data('Si')
    dict_conf = ion.get_info_soft_confinement()
    assert dict_conf == {}

    ion = generate_ion_data('Si_with_conf')
    dict_conf = ion.get_info_soft_confinement()
    assert dict_conf == {'3s': [2.0, 0.3]}


# The `get_content_ascii_format` is tested is the calculations/test_siesta.py
# the rest of methods in utils/test_pao_manager.py
