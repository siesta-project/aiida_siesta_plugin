#!/usr/bin/env python
# -*- coding: utf-8 -*-

import pytest
import configparser
from aiida import load_dbenv

load_dbenv()

testconfig = configparser.ConfigParser()
testconfig.read("testconfig.ini")


@pytest.fixture(scope="session")
def siesta_develop():
    sd = {}
    from aiida.orm.code import Code
    sd["siesta_code"] = Code.get_from_string(testconfig['codes']['siesta'])
    return sd
