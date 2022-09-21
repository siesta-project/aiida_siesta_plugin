# -*- coding: utf-8 -*-
"""Tests for the `SiestaParser`."""

from aiida import orm
from aiida.common import AttributeDict
import pytest


def test_stm_default(aiida_profile, fixture_localhost, generate_calc_job_node,
    generate_parser, generate_structure, data_regression):
    """Test a parser of a stm calculation.
    The output is created by running a dead simple SCF calculation for a silicon structure. This test should test the
    standard parsing of the stdout content and XML file stored in the standard results node.
    """
    name = 'default'
    entry_point_calc_job = 'siesta.stm'
    entry_point_parser = 'siesta.stm'

    inputs = AttributeDict({
        'spin_option' : orm.Str("q")
    })

    attributes=AttributeDict({'input_filename':'stm.in', 'output_filename':'stm.out'})

    node = generate_calc_job_node(entry_point_calc_job, fixture_localhost, name, inputs, attributes)
    parser = generate_parser(entry_point_parser)
    results, calcfunction = parser.parse_from_node(node, store_provenance=False)

    assert calcfunction.is_finished, calcfunction.exception
    assert calcfunction.is_finished_ok, calcfunction.exit_message
    assert not orm.Log.objects.get_logs_for(node)
    assert 'stm_array' in results
    assert 'output_parameters' in results

    data_regression.check({
        'stm_array': results['stm_array'].attributes,
        'output_parameters': results['output_parameters'].attributes,
    })
