#!/usr/bin/env python
# -*- coding: utf-8 -*-


def test_siesta_calculation_entry_point(siesta_develop):
    from aiida.orm import CalculationFactory
    siesta_calculation = CalculationFactory('siesta.siesta')
    assert siesta_calculation is not None


def test_siesta_parser_entry_point(siesta_develop):
    from aiida.parsers import ParserFactory
    siesta_parser = ParserFactory('siesta.parser')
    assert siesta_parser is not None


def test_siesta_psf_entry_point(siesta_develop):
    from aiida.orm import DataFactory
    siesta_psf = DataFactory('siesta.psf')
    assert siesta_psf is not None
