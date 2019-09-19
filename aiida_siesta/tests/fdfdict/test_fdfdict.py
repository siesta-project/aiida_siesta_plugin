#!/usr/bin/env python
# -*- coding: utf-8 -*-


from __future__ import absolute_import
def test_fdfdict():
    from aiida_siesta.calculations.tkdict import FDFDict

    f = FDFDict()

    # insertion and saving
    f["A-_.B"] = 1
    assert("a_b" in list(f._storage.keys()))
    assert f["a_:.b"] == 1
    assert f.get_last_key("a::_b") == "A-_.B"
