# -*- coding: utf-8 -*-
# pylint: disable=invalid-name,redefined-outer-name
"""Tests for the `SiestaParser`."""
from __future__ import absolute_import

import pytest

from aiida import orm
from aiida.common import AttributeDict

def test_siesta_default(aiida_profile, fixture_localhost, generate_calc_job_node, 
    generate_parser, generate_structure, data_regression):
    """Test a parser of a siesta calculation.
    The output is created by running a dead simple SCF calculation for a silicon structure. This test should test the
    standard parsing of the stdout content and XML file stored in the standard results node.
    """
    name = 'default'
    entry_point_calc_job = 'siesta.siesta'
    entry_point_parser = 'siesta.parser'

    structure=generate_structure()

    inputs = AttributeDict({
        'structure': structure
    })

    attributes=AttributeDict({'input_filename':'aiida.fdf', 'output_filename':'aiida.out', 'prefix':'aiida'})

    node = generate_calc_job_node(entry_point_calc_job, fixture_localhost, name, inputs, attributes)
    parser = generate_parser(entry_point_parser)
    results, calcfunction = parser.parse_from_node(node, store_provenance=False)

    assert calcfunction.is_finished, calcfunction.exception
    assert calcfunction.is_finished_ok, calcfunction.exit_message
    assert not orm.Log.objects.get_logs_for(node)
    assert 'forces_and_stress' in results
    assert 'output_parameters' in results

    data_regression.check({
        'forces_and_stress': results['forces_and_stress'].attributes,
        'output_parameters': results['output_parameters'].get_dict()
    })

# As it is implemented now, there is no point to test also the case bandslines as
# I assert the attributes of bands, not the actual array!
def test_siesta_bandspoints(aiida_profile, fixture_localhost, generate_calc_job_node,
    generate_parser, generate_structure, data_regression):
    """Test parsing of bands in a siesta calculation when the bandspoints option is set in the submission file.
    """
    name = 'bandspoints'
    entry_point_calc_job = 'siesta.siesta'
    entry_point_parser = 'siesta.parser'

    structure=generate_structure()
    bandskpoints = orm.KpointsData()
    kpp = [(0.500,  0.250, 0.750), (0.500,  0.500, 0.500), (0., 0., 0.)]
    bandskpoints.set_cell(structure.cell, structure.pbc)
    bandskpoints.set_kpoints(kpp)

    inputs = AttributeDict({
        'structure': structure,
        'bandskpoints': bandskpoints
    })

    attributes=AttributeDict({'input_filename':'aiida.fdf', 'output_filename':'aiida.out', 'prefix':'aiida'})

    node = generate_calc_job_node(entry_point_calc_job, fixture_localhost, name, inputs, attributes)
    parser = generate_parser(entry_point_parser)
    results, calcfunction = parser.parse_from_node(node, store_provenance=False)

    assert calcfunction.is_finished, calcfunction.exception
    assert calcfunction.is_finished_ok, calcfunction.exit_message
    assert not orm.Log.objects.get_logs_for(node)
    assert 'forces_and_stress' in results
    assert 'output_parameters' in results
    assert 'bands' in results 

    data_regression.check({
        'forces_and_stress': results['forces_and_stress'].attributes,
        'output_parameters': results['output_parameters'].get_dict(),
        'bands': results['bands'].attributes
    })

