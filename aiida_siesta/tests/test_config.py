#!/usr/bin/env python
# -*- coding: utf-8 -*-


def test_code_config(configure):
    from aiida.orm.code import Code
    code = Code.get_from_string('tsiesta@localhost')
    assert code is not None
    assert code.description == 'test siesta code object'
    assert code.pk == 1
    assert code.get_remote_exec_path() == '/usr/bin/siesta'

    computer = code.get_remote_computer()
    assert computer is not None
    assert computer.pk == 1
