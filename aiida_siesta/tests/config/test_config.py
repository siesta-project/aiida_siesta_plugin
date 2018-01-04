#!/usr/bin/env python
# -*- coding: utf-8 -*-


def test_development_config(siesta_develop):
    code = siesta_develop["code"]
    assert code is not None
    assert code.description == 'SIESTA code built in the Docker image'
    assert code.pk == 1
    assert code.get_remote_exec_path() == '/usr/local/bin/siesta'

    computer = code.get_remote_computer()
    assert computer is not None
    assert computer.pk == 1
    assert computer.name == 'develop'
