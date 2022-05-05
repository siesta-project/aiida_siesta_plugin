# -*- coding: utf-8 -*-
def test_siesta_calculation_entry_point():
    from aiida.plugins import CalculationFactory
    siesta_calculation = CalculationFactory('siesta.siesta')
    assert siesta_calculation is not None


def test_siesta_parser_entry_point():
    from aiida.plugins import ParserFactory
    siesta_parser = ParserFactory('siesta.parser')
    assert siesta_parser is not None
