import pytest
from aiida import orm
from aiida.common import AttributeDict


def test_siesta_default(aiida_profile, fixture_localhost, generate_calc_job_node, 
    generate_parser, generate_structure, generate_basis, data_regression):
    """
    Test a parser of a siesta calculation.
    The output is created by running a dead simple SCF calculation for a silicon structure. 
    We test the standard parsing of the XML file stored in the standard results node.
    No other file (time.json, MESSAGES, BASIS_ENTHALPY) is present. Therefore no error check is done,
    but the appropriate warnings are issued.
    """

    name = 'default'
    entry_point_calc_job = 'siesta.siesta'
    entry_point_parser = 'siesta.parser'

    structure=generate_structure()
    basis=generate_basis().get_dict()
    basis["floating_sites"] = [{"name":'Si_bond',"symbols":'Si',"position": ( 0.125, 0.125, 0.125)}]

    inputs = AttributeDict({
        'structure': structure,
        'basis' : orm.Dict(dict=basis)
    })

    attributes=AttributeDict({'input_filename':'aiida.fdf', 'output_filename':'aiida.out', 'prefix':'aiida'})

    node = generate_calc_job_node(entry_point_calc_job, fixture_localhost, name, inputs, attributes)
    parser = generate_parser(entry_point_parser)
    results, calcfunction = parser.parse_from_node(node, store_provenance=False)

    assert calcfunction.is_finished
    assert calcfunction.exception is None
    assert calcfunction.is_finished_ok
    assert calcfunction.exit_message is None
    assert not orm.Log.objects.get_logs_for(node)
    assert 'forces_and_stress' in results
    assert 'output_parameters' in results
    assert 'output_structure' in results
    assert 'ion_files' in results
    assert 'Si' in results['ion_files']
    assert 'Si_bond' in results['ion_files']

    data_regression.check({
        'forces_and_stress': results['forces_and_stress'].attributes,
        'output_parameters': results['output_parameters'].get_dict()
    })



def test_siesta_no_ion(aiida_profile, fixture_localhost, generate_calc_job_node, 
    generate_parser, generate_basis, generate_structure, generate_ion_data, data_regression):
    """
    Test a parser of a siesta calculation, but the .ion.xml are not found
    """

    name = 'no_ion'
    entry_point_calc_job = 'siesta.siesta'
    entry_point_parser = 'siesta.parser'

    structure=generate_structure()
    basis=generate_basis().get_dict()
    basis["floating_sites"] = [{"name":'Si_bond',"symbols":'Si',"position": ( 0.125, 0.125, 0.125)}]

    inputs = AttributeDict({
        'structure': structure,
        'basis' : orm.Dict(dict=basis)
    })

    attributes=AttributeDict({'input_filename':'aiida.fdf', 'output_filename':'aiida.out', 'prefix':'aiida'})

    node = generate_calc_job_node(entry_point_calc_job, fixture_localhost, name, inputs, attributes)
    parser = generate_parser(entry_point_parser)
    results, calcfunction = parser.parse_from_node(node, store_provenance=False)

    assert calcfunction.is_finished
    assert calcfunction.exception is None
    assert calcfunction.is_finished_ok
    assert calcfunction.exit_message is None
    log = orm.Log.objects.get_logs_for(node)
    assert len(log) == 3
    sum_strings = log[0].message + log[1].message + log[2].message
    assert "no ion file retrieved for Si_bond" in sum_strings
    assert "no ion file retrieved for Si" in sum_strings

    #However the warning should not be present if ions are in input
    inputs["ions"] =  {"Si_bond":generate_ion_data("Si")}
    node = generate_calc_job_node(entry_point_calc_job, fixture_localhost, name, inputs, attributes)
    parser = generate_parser(entry_point_parser)
    results, calcfunction = parser.parse_from_node(node, store_provenance=False)
    log = orm.Log.objects.get_logs_for(node)
    assert len(log) == 0


# As it is implemented now, there is no point to test also the case bandslines as
# I assert the attributes of bands, not the actual array!
def test_siesta_bandspoints(aiida_profile, fixture_localhost, generate_calc_job_node,
    generate_parser, generate_structure, data_regression):
    """
    Test parsing of bands in a siesta calculation when the bandspoints option is set in the submission file.
    Also the time.json, MESSAGES and BASIS_ENTHALPY files are added, therefore their parsing is tested as well. The MESSAGES
    file is the standard containing only "INFO: Job completed". 
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

    assert calcfunction.is_finished
    assert calcfunction.exception is None
    assert calcfunction.is_finished_ok
    assert calcfunction.exit_message is None
    assert not orm.Log.objects.get_logs_for(node)
    assert 'forces_and_stress' in results
    assert 'output_parameters' in results
    assert 'output_structure' in results
    assert 'bands' in results 

    data_regression.check({
        'forces_and_stress': results['forces_and_stress'].attributes,
        'output_parameters': results['output_parameters'].get_dict(),
        'bands': results['bands'].attributes
    })


def test_siesta_optical(aiida_profile, fixture_localhost, generate_calc_job_node,
    generate_parser, generate_structure, data_regression):
    """
    Test parsing of eps2 in a siesta calculation when the bandspoints option is set in the submission file.
    """

    name = 'optical'
    entry_point_calc_job = 'siesta.siesta'
    entry_point_parser = 'siesta.parser'

    structure=generate_structure()

    optical = {
        "optical-calculation" : True,
        "optical-broaden": "0.5 eV",
        "optical-polarization-type": "polarized",
        "%block optical-vector": "\n 1.0 0.0 0.0 \n%endblock optical-vector",
        "%block optical-mesh": "\n 3 3 3 \n%endblock optical-mesh"
    }

    inputs = AttributeDict({
        'structure': structure,
        'optical': orm.Dict(dict=optical)
    })

    attributes=AttributeDict({'input_filename':'aiida.fdf', 'output_filename':'aiida.out', 'prefix':'aiida'})

    node = generate_calc_job_node(entry_point_calc_job, fixture_localhost, name, inputs, attributes)
    parser = generate_parser(entry_point_parser)
    results, calcfunction = parser.parse_from_node(node, store_provenance=False)

    assert calcfunction.is_finished
    assert calcfunction.exception is None
    #assert calcfunction.is_finished_ok
    assert calcfunction.exit_message is None
    assert not orm.Log.objects.get_logs_for(node)
    assert 'forces_and_stress' in results
    assert 'output_parameters' in results
    assert 'output_structure' in results
    assert 'optical_eps2' in results

    #data_regression.check({
    #    'forces_and_stress': results['forces_and_stress'].attributes,
    #    'output_parameters': results['output_parameters'].get_dict(),
    #    'bands': results['bands'].attributes
    #})


def test_siesta_empty_messages(aiida_profile, fixture_localhost, generate_calc_job_node, 
    generate_parser, generate_structure, data_regression):
    """
    An empty MESSAGES file is parsed. The parser goes through all the error checks but no
    recognizable error is detected. Therefore the calculation finishes with UNEXPECTED_TERMINATION
    error. However all the outputs are returned because we assume the .xml file was produced
    sucessfully. This example mimics the sitiuation when the calculation is interupted just before
    its completition (maybe externally).
    The file time.json is not present at this time.
    """

    name = 'empty_messages'
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

    assert calcfunction.is_finished
    assert calcfunction.exception is None
    assert not calcfunction.is_finished_ok
    assert calcfunction.exit_message == 'Statement "Job completed" not detected, unknown error'
    logs = orm.Log.objects.get_logs_for(node)[0]
    assert "Job completed" in logs.message
    assert 'forces_and_stress' in results
    assert 'output_parameters' in results
    assert 'output_structure' in results

    data_regression.check({
        'forces_and_stress': results['forces_and_stress'].attributes,
        'output_parameters': results['output_parameters'].get_dict()
    })


def test_siesta_no_scf_conv(aiida_profile, fixture_localhost, generate_calc_job_node,
    generate_parser, generate_structure, data_regression):
    """
    Test a parser in the situation when siesta stops with "SCF_NOT_CONV" situation. It is produced with
    modern versions of the code that return "FATAL" in this case (unless scf-must-converge F) is set by the
    user. An output_parameters node is always produced.
    The file time.json is present at this time.
    """

    name = 'no_scf_conv'
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

    assert calcfunction.is_finished
    assert calcfunction.exception is None
    assert not calcfunction.is_finished_ok
    assert calcfunction.exit_message == 'Calculation did not reach scf convergence!'
    logs = orm.Log.objects.get_logs_for(node)
    for log in logs:
        if "SCF_NOT_CONV" in log.message:
            mylog=log
    assert "SCF_NOT_CONV" in mylog.message
    assert 'output_parameters' in results

    data_regression.check({'output_parameters': results['output_parameters'].get_dict()})



def test_siesta_no_geom_conv(aiida_profile, fixture_localhost, generate_calc_job_node,
    generate_parser, generate_structure, data_regression):
    """
    Test a parser in the situation when siesta stops with "GEOM_NOT_CONV" situation. It is produced with
    modern versions of the code that return "FATAL" in this case (unless scf-must-converge F) is set by the
    user. An output_parameters node is always produced.
    The file time.json is present at this time.
    """

    name = 'no_geom_conv'
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

    assert calcfunction.is_finished
    assert calcfunction.exception is None
    assert not calcfunction.is_finished_ok
    assert calcfunction.exit_message == 'Calculation did not reach geometry convergence!'
    logs = orm.Log.objects.get_logs_for(node)
    for log in logs:
        if "GEOM_NOT_CONV" in log.message:
            mylog=log
    assert "GEOM_NOT_CONV" in mylog.message
    assert 'output_parameters' in results
    assert 'output_structure' in results

    data_regression.check({
        'output_structure': results['output_structure'].attributes,
        'output_parameters': results['output_parameters'].get_dict()
        })


def test_siesta_bands_error(aiida_profile, fixture_localhost, generate_calc_job_node,
    generate_parser, generate_structure, data_regression):
    """
    Test parsing of bands in a situation when the bands file is truncated.
    Also the time.json and MESSAGES file are added, therefore their parsing is tested as well. The MESSAGES
    file is the standard containing only "INFO: Job completed". 
    """

    name = 'bands_error'
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

    assert calcfunction.is_finished
    assert calcfunction.exception is None
    assert not calcfunction.is_finished_ok
    assert calcfunction.exit_message == 'Failure while parsing the bands file'
    assert 'output_parameters' in results
    assert 'output_structure' in results


def test_basis_entalpies(aiida_profile, fixture_localhost, generate_calc_job_node, 
    generate_parser, generate_structure, data_regression):
    """
    Test the parsing of the BASIS_ENTHALPY and BASIS_HARRIS_ENTHALPY files.
    """

    name = 'basis_enthalpies'

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

    assert calcfunction.is_finished
    assert calcfunction.exception is None
    assert not calcfunction.is_finished_ok

    assert calcfunction.exit_message == 'Statement "Job completed" not detected, unknown error'
    logs = orm.Log.objects.get_logs_for(node)[0]
    assert "Job completed" in logs.message
    assert 'forces_and_stress' in results
    assert 'output_parameters' in results
    assert 'output_structure' in results

    data_regression.check({
        'forces_and_stress': results['forces_and_stress'].attributes,
        'output_parameters': results['output_parameters'].get_dict()
    })


def test_siesta_optical_error(aiida_profile, fixture_localhost, generate_calc_job_node,
    generate_parser, generate_structure, data_regression):
    """
    Test parsing of eps2 in a siesta calculation when the bandspoints option is set in the submission file.
    """

    name = 'optical_error'
    entry_point_calc_job = 'siesta.siesta'
    entry_point_parser = 'siesta.parser'

    structure=generate_structure()

    optical = {
        "optical-calculation" : True,
        "optical-broaden": "0.5 eV",
        "optical-polarization-type": "polarized",
        "%block optical-vector": "\n 1.0 0.0 0.0 \n%endblock optical-vector",
        "%block optical-mesh": "\n 3 3 3 \n%endblock optical-mesh"
    }

    inputs = AttributeDict({
        'structure': structure,
        'optical': orm.Dict(dict=optical)
    })

    attributes=AttributeDict({'input_filename':'aiida.fdf', 'output_filename':'aiida.out', 'prefix':'aiida'})

    node = generate_calc_job_node(entry_point_calc_job, fixture_localhost, name, inputs, attributes)
    parser = generate_parser(entry_point_parser)
    results, calcfunction = parser.parse_from_node(node, store_provenance=False)

    assert calcfunction.is_finished
    assert calcfunction.exception is None
    assert not calcfunction.is_finished_ok
    assert calcfunction.exit_message == 'Optical calculation requested, but file is not present'
    assert 'forces_and_stress' in results
    assert 'output_parameters' in results
    assert 'output_structure' in results
    assert 'optical_eps2' not in results
