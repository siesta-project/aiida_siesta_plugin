def test_siesta_calculation_entry_point():
    from aiida.plugins import CalculationFactory
    siesta_calculation = CalculationFactory('siesta.siesta')
    assert siesta_calculation is not None


def test_siesta_parser_entry_point():
    from aiida.plugins import ParserFactory
    siesta_parser = ParserFactory('siesta.parser')
    assert siesta_parser is not None


def test_siesta_psf_entry_point():
    from aiida.plugins import DataFactory
    siesta_psf = DataFactory('siesta.psf')
    assert siesta_psf is not None

def test_siesta_psml_entry_point():
    from aiida.plugins import DataFactory
    siesta_psml = DataFactory('siesta.psf')
    assert siesta_psml is not None

