#!/usr/bin/env python
# -*- coding: utf-8 -*-

from plum.util import load_class


def test_psf_command():
    pseudo_cmd = load_class('aiida_siesta.commands.psfdata._Psf')

    assert 'uploadfamily' in dir(pseudo_cmd)
    assert 'listfamilies' in dir(pseudo_cmd)
    assert '_import_psf' in dir(pseudo_cmd)
    assert 'exportfamily' in dir(pseudo_cmd)
