#!/usr/bin/env python
# -*- coding: utf-8 -*-

import pytest
from aiida import load_dbenv

load_dbenv()

@pytest.fixture(scope="session")
def siesta_develop():
    sd = {}
    from aiida.orm.code import Code
    sd["code"] = Code.get_from_string('siesta@develop')
    return sd
