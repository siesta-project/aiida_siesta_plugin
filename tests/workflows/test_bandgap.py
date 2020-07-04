#!/usr/bin/env runaiida
import pytest
from aiida import orm


@pytest.fixture
def generate_workchain_bandgap(generate_psml_data, fixture_code, fixture_localhost, generate_workchain, 
        generate_structure, generate_param, generate_basis, generate_kpoints_mesh,
        generate_calc_job_node, generate_parser):
    """Generate an instance of a `BandgapWorkChain`."""

    def _generate_workchain_bandgap(bands=False):

        entry_point_wc = 'siesta.bandgap'
        entry_point_code = 'siesta.siesta'

        psml = generate_psml_data('Si')

        structure = generate_structure()

        inputs = {
            'code': fixture_code(entry_point_code),
            'structure': structure,
            'kpoints': generate_kpoints_mesh(2),
            'parameters': generate_param(),
            'basis': generate_basis(),
            'pseudos': {
                'Si': psml,
                'SiDiff': psml
            },
            'options': orm.Dict(dict={
               'resources': {'num_machines': 1  },
               'max_wallclock_seconds': 1800,
               'withmpi': False,
               })
        }

        if bands:
            bandskpoints = orm.KpointsData()
            kpp = [(0.500,  0.250, 0.750), (0.500,  0.500, 0.500), (0., 0., 0.)]
            bandskpoints.set_cell(structure.cell, structure.pbc)
            bandskpoints.set_kpoints(kpp)
            inputs["bandskpoints"] = bandskpoints


        process = generate_workchain(entry_point_wc, inputs)

        return process

    return _generate_workchain_bandgap

def test_setup_fail(aiida_profile, generate_workchain_bandgap):
    """Test `BangapWorkChain.setup`."""
    with pytest.raises(ValueError):
        process = generate_workchain_bandgap()
        process.preprocess()
        process.setup()
        process.prepare_inputs()

def test_setup(aiida_profile, generate_workchain_bandgap):
    """Test `BangapWorkChain.setup`."""
    
    process = generate_workchain_bandgap(bands=True)
    process.preprocess()
    process.setup()
    process.prepare_inputs()

    assert isinstance(process.ctx.inputs, dict)

def test_postprocess(aiida_profile, generate_workchain_bandgap):
#    """Test `BangapWorkChain.setup`."""
#
    process = generate_workchain_bandgap(bands=True)
    process.out('output_parameters', orm.Dict(dict={'E_Fermi' : -1}))
    bands = orm.BandsData()
    bkp = process.inputs.bandskpoints
    bands.set_kpoints(bkp.get_kpoints(cartesian=True))
    #bands.labels = bkp.labels
    import numpy as np
    bands.set_bands(np.array([[-2,1,3],[-2,1,3],[-2,1,3]]), units="eV")
    process.out('bands', bands)
    process.postprocess()

    res = process.outputs["band_gap_info"]
    assert isinstance(res,orm.Dict)
    assert res['is_insulator']
    assert res['band_gap'] == 3.0

    #assert isinstance(process.ctx.inputs, dict)

