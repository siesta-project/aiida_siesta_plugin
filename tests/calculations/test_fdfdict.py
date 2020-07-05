#!/usr/bin/env python
from aiida_siesta.calculations.tkdict import FDFDict


def test_fdfdict_wrong_argument():
    """
    Simple test of a FDFDict class instance with an unaccepted
    argument passed
    """

    import pytest
    with pytest.raises(RuntimeError):
        f = FDFDict("w")

def test_fdfdict_no_argument():
    """
    Simple test of a FDFDict class instance with no
    argument passed, also tests `get_last_untranslated_key`
    """

    f = FDFDict()

    # insertion and saving
    f["A-_.B"] = 1
    assert("a_b" in list(f._storage.keys()))
    assert f["a_:.b"] == 1
    assert f.get_last_untranslated_key("a::_b") == "A-_.B"

def test_fdfdict():
    """
    Simple test of a FDFDict class instance with a dictionary passed
    """

    inp_dict = {"w":3,"e":4,"w--":5}
    
    f = FDFDict(inp_dict)

    assert f["w"] == 5
    
    f["e-"] = 4
    assert f.get_last_untranslated_key("e") == "e-"

    assert "w" in f

    for i in f:
        assert i in ["w","e"]

def test_fdfdict_methods():
    """
    Simple test of a FDFDict class instance with a dictionary passed
    """

    inp_dict = {"w":3,"e":4,"w--":5}

    f = FDFDict(inp_dict)

    assert sorted(f.values()) == sorted([5,4])
    assert sorted(f.keys()) == sorted(["w","e"])
    assert sorted(f.untranslated_keys()) == sorted(["w--","e"])
    assert sorted(f.items()) == sorted([("w",5),("e",4)])
    assert sorted(f.untranslated_items()) == sorted([("w--",5),("e",4)])
    assert f.get_dict() == {"w":5,"e":4}
    assert f.get_untranslated_dict() == {"w--":5,"e":4}
