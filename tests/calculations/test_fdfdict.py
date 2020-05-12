#!/usr/bin/env python

def test_fdfdict():
    """
    Simple test of a FDFDict class instance.
    """
    from aiida_siesta.calculations.tkdict import FDFDict

    f = FDFDict()

    # insertion and saving
    f["A-_.B"] = 1
    assert("a_b" in list(f._storage.keys()))
    assert f["a_:.b"] == 1
    assert f.get_last_key("a::_b") == "A-_.B"
