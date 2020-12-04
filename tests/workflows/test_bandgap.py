#!/usr/bin/env runaiida
import pytest
from plumpy import ProcessState
from aiida import orm
from aiida.common import (LinkType, AttributeDict)
from aiida.engine import ExitCode

@pytest.fixture
def generate_workchain_bandgap(generate_psml_data, fixture_code, fixture_localhost, generate_workchain, 
        generate_structure, generate_param, generate_basis, generate_kpoints_mesh,
        generate_calc_job_node, generate_parser):
    """Generate an instance of a `BandgapWorkChain`."""

    def _generate_workchain_bandgap(bands=False,relax=False):

        entry_point_wc = 'siesta.bandgap'
        entry_point_code = 'siesta.siesta'

        psml = generate_psml_data('Si')

        inputs = {
            'code': fixture_code(entry_point_code),
            'structure': generate_structure(),
            'kpoints': generate_kpoints_mesh(2),
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

        if relax:
            inputs["parameters"] = generate_param()
        else:
            param = generate_param().get_dict()
            for item in param.copy():
                if item.startswith("md"):
                    param.pop(item)
            inputs["parameters"] = orm.Dict(dict=param)

        if bands:
            bandskpoints = orm.KpointsData()
            kpp = [(0.500,  0.250, 0.750), (0.500,  0.500, 0.500), (0., 0., 0.)]
            bandskpoints.set_kpoints(kpp)
            inputs["bandskpoints"] = bandskpoints


        process = generate_workchain(entry_point_wc, inputs)

        return process

    return _generate_workchain_bandgap

def test_preproc_and_add_kpb(aiida_profile, generate_workchain_bandgap):
    """Test `BangapWorkChain.preprocess`."""

    process = generate_workchain_bandgap(bands=True)
    process.preprocess()
    assert not process.ctx.need_to_generate_bandskp
    assert not process.ctx.need_fin_step

    process = generate_workchain_bandgap(bands=False,relax=False)
    process.preprocess()
    assert process.ctx.need_to_generate_bandskp
    assert not process.ctx.need_fin_step 

    res = process.run_siesta_wc()

    assert "bandskpoints" in res['workchain_base'].inputs


def test_preproc_and_relax(aiida_profile, fixture_localhost, fixture_code, generate_psml_data,
        generate_structure, generate_wc_job_node, generate_workchain_bandgap):
    """Test `BangapWorkChain.setup`."""

    process = generate_workchain_bandgap(bands=False,relax=True)
    process.preprocess()
    assert not process.ctx.need_to_generate_bandskp
    assert process.ctx.need_fin_step

    res = process.run_siesta_wc()

    assert not "bandskpoints" in res['workchain_base'].inputs

    psml = generate_psml_data("Si")

    inputs = AttributeDict({
        'structure': generate_structure(),
        'code': fixture_code("siesta.siesta"),
        'parameters': orm.Dict(dict={"md": 3, "ee":4}),
        'options': orm.Dict(dict={'resources': {'num_machines': 1  },'max_wallclock_seconds': 1800,'withmpi': False}),
        'pseudos': {'Si': psml,'SiDiff': psml},
    })
    fin_basewc = generate_wc_job_node("siesta.base", fixture_localhost, inputs)
    fin_basewc.set_process_state(ProcessState.FINISHED)
    fin_basewc.set_exit_status(ExitCode(0).status)
    out_struct = generate_structure()
    out_struct.store()
    out_struct.add_incoming(fin_basewc, link_type=LinkType.RETURN, link_label='output_structure')

    process.ctx.workchain_base = fin_basewc

    finwc = process.run_last()

    assert "bandskpoints" in finwc['final_run'].inputs
    assert finwc['final_run'].inputs.parameters.get_dict() == {"ee":4}


def test_final_run(aiida_profile, fixture_localhost, generate_workchain_bandgap, generate_wc_job_node):
    """Test `BangapWorkChain.setup`."""

    process = generate_workchain_bandgap(bands=True)

    process.ctx.need_fin_step = False

    fin_basewc = generate_wc_job_node("siesta.base", fixture_localhost)
    fin_basewc.set_process_state(ProcessState.FINISHED)
    fin_basewc.set_exit_status(ExitCode(0).status)
    
    out_par = orm.Dict(dict={'E_Fermi' : -1})
    out_par.store()
    out_par.add_incoming(fin_basewc, link_type=LinkType.RETURN, link_label='output_parameters')
    
    out_force_stress = orm.ArrayData()
    out_force_stress.store()
    out_force_stress.add_incoming(fin_basewc, link_type=LinkType.RETURN, link_label='forces_and_stress')
    
    remote_folder = orm.RemoteData(computer=fixture_localhost, remote_path='/tmp')
    remote_folder.store()
    remote_folder.add_incoming(fin_basewc, link_type=LinkType.RETURN, link_label='remote_folder')

    bands = orm.BandsData()
    bkp = process.inputs.bandskpoints
    bands.set_kpointsdata(bkp)
    import numpy as np
    bands.set_bands(np.array([[-2,1,3],[-2,1,3],[-2,1,3]]), units="eV")
    out_bands = bands
    out_bands.store()
    out_bands.add_incoming(fin_basewc, link_type=LinkType.RETURN, link_label='bands')
   
    process.ctx.workchain_base = fin_basewc
    process.run_results()

    res = process.outputs["band_gap_info"]
    assert isinstance(res,orm.Dict)
    assert res['is_insulator']
    assert res['band_gap'] == 3.0
