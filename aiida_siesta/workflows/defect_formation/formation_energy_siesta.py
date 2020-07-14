# -*- coding: utf-8 -*-
########################################################################################
# Copyright (c), The AiiDA-Defects authors. All rights reserved.                       #
#                                                                                      #
# AiiDA-Defects is hosted on GitHub at https://github.com/ConradJohnston/aiida-defects #
# For further information on the license, see the LICENSE.txt file                     #
########################################################################################
from __future__ import absolute_import

import numpy as np

from aiida import orm
from aiida.engine import WorkChain, calcfunction, ToContext, if_, submit
from aiida.plugins import WorkflowFactory
from aiida_siesta.workflows.base import SiestaBaseWorkChain

from .formation_energy_base import FormationEnergyWorkchainBase
from .utils import run_siesta_calculation
from .utils import get_raw_formation_energy, get_corrected_formation_energy, get_corrected_aligned_formation_energy

from aiida.common import AttributeDict
from aiida_siesta.calculations.tkdict import FDFDict

class FormationEnergyWorkchainSIESTA(FormationEnergyWorkchainBase):
    """
    Compute the formation energy for a given defect using SIESTA
    """
    @classmethod
    def define(cls, spec):
        super(FormationEnergyWorkchainSIESTA, cls).define(spec)

        # DFT and DFPT calculations with QuantumESPRESSO are handled with different codes, so here
        # we keep track of things with two separate namespaces. An additional code, and an additional
        # namespace, is used for postprocessing
        #spec.expose_inputs(SiestaBaseWorkChain, exclude=('metadata',))
        #spec.inputs._ports['pseudos'].dynamic = True  #Temporary fix to issue #135 plumpy
        # spec.inputs._ports["pseudos.defect"].dynamic = True
        # DFT inputs for Host (SIESTA)
        spec.input_namespace("siesta.dft.supercell_host",
            help="The siesta code to use for the calculations")
        spec.input_namespace("siesta.dft.supercell_defect_q0",
            help="The siesta code to use for the calculations")
        spec.input_namespace("siesta.dft.supercell_defect_q",
            help="The siesta code to use for the calculations")
        # HOST Inputs
        spec.input("siesta.dft.supercell_host.code",
            valid_type=orm.Code,
            help="The siesta code to use for the calculations")
        spec.input("siesta.dft.supercell_host.kpoints",
            valid_type=orm.KpointsData,
            help="The k-point grid to use for the calculations")
        spec.input("siesta.dft.supercell_host.basis",
            valid_type=orm.Dict,
            help="The siesta basis to use for the host calculations")
        spec.input("siesta.dft.supercell_host.parameters",
            valid_type=orm.Dict,
            help="Parameters for the SIESTA calcuations. Some will be set automatically")
        spec.input("siesta.dft.supercell_host.options",
            valid_type=orm.Dict,
            help="options for the SIESTA calcuations.")
        # Defect_q0 without charge
        spec.input("siesta.dft.supercell_defect_q0.code",
            valid_type=orm.Code,
            help="The siesta code to use for the calculations")
        spec.input("siesta.dft.supercell_defect_q0.kpoints",
            valid_type=orm.KpointsData,
            help="The k-point grid to use for the calculations")
        spec.input("siesta.dft.supercell_defect_q0.basis",
            valid_type=orm.Dict,
            help="The siesta basis to use for the host calculations")
        spec.input("siesta.dft.supercell_defect_q0.parameters",
            valid_type=orm.Dict,
            help="Parameters for the SIESTA calcuations. Some will be set automatically")
        spec.input("siesta.dft.supercell_defect_q0.options",
            valid_type=orm.Dict,
            help="options for the SIESTA calcuations.")
        # DFT inputs for Defect With Charge (SIESTA)
        spec.input("siesta.dft.supercell_defect_q.code",
            valid_type=orm.Code,
            help="The siesta code to use for the calculations")
        spec.input("siesta.dft.supercell_defect_q.kpoints",
            valid_type=orm.KpointsData,
            help="The k-point grid to use for the calculations")
        spec.input("siesta.dft.supercell_defect_q.basis",
            valid_type=orm.Dict,
            help="The siesta basis to use for the host calculations")
        spec.input("siesta.dft.supercell_defect_q.parameters",
            valid_type=orm.Dict,
            help="Parameters for the SIESTA calcuations. Some will be set automatically")
        spec.input("siesta.dft.supercell_defect_q.options",
            valid_type=orm.Dict,
            help="options for the SIESTA calcuations.")

        

        spec.outline(
                     cls.setup,
                     if_(cls.correction_required)(
                                                  if_(cls.is_gaussian_scheme)
                                                  (
                                                  cls.raise_not_implemented
                                                  #    cls.prep_dft_calcs_gaussian_correction,
                                                  #    cls.check_dft_calcs_gaussian_correction,
                                                  #    cls.get_dft_potentials_gaussian_correction,
                                                  #    cls.check_dft_potentials_gaussian_correction,
                                                  #    cls.run_gaussian_correction_workchain),
                                                  #if_(cls.is_point_scheme)(
                                                  #    cls.raise_not_implemented
                                                  #cls.check_correction_workchain
                                                  ))
                    #cls.prep_dft_calcs_no_correction
                    #cls.check_dft_calcs_no_correction
                    #cls.compute_formation_energy
                    )

    def prep_dft_calcs_no_correction(self):
 
        """
        Submit All DFT Calculations
        """
        self.report("Setting Up the No correction Formation Energy Workchain ")
        # For the Host
        siesta_inputs = self.inputs.siesta.dft.supercell_host.code.get_builder()
        siesta_inputs = AttributeDict(self.exposed_inputs(SiestaBaseWorkChain))
        siesta_inputs.structure = self.inputs.host_structure
        siesta_inputs.parameters = self.input.dft.supercell_host.parameters
        siesta_inputs.kpoints = self.inputs.siesta.dft.supercell_host.kpoints
        siesta_inputs.options = self.input.siesta.dft.supercell_host.options
        siesta_inputs.basis = self.input.siesta.dft.supercell_host.basis
        future = self.submit(SiestaBaseWorkChain,**siesta_inputs)
        #future = self.submit(**siesta_inputs)
        self.report(
            'Launching SIESTA for host structure (PK={})  (PK={})'
            .format(self.inputs.host_structure.pk, future.pk))
        self.to_context(**{'calc_host': future})        
        # For Defect structure; neutral charge state
        siesta_inputs = self.inputs.siesta.dft.supercell_defect_q0.code.get_builder()
        siesta_inputs = AttributeDict(self.exposed_inputs(SiestaBaseWorkChain))
        siesta_inputs.structure = self.inputs.defect_structure
        siesta_inputs.parameters = self.input.dft.supercell_defect_q0.parameters
        siesta_inputs.kpoints = self.inputs.siesta.dft.supercell_defect_q0.kpoints
        siesta_inputs.options = self.input.siesta.dft.supercell_defect_q0.options
        siesta_inputs.basis = self.input.siesta.dft.supercell_defect_q0.basis
        future = self.submit(SiestaBaseWorkChain,**siesta_inputs)
        #future = self.submit(**siesta_inputs)
        self.report(
            'Launching SIESTA for defect structure (PK={}) with charge {} (PK={})'
            .format(self.inputs.defect_structure.pk, "0.0", future.pk))
        self.to_context(**{'calc_defect_q0': future})
        # For Defect structure; target charge state
        siesta_inputs = self.inputs.siesta.dft.supercell_defect_q.code.get_builder()
        siesta_inputs = AttributeDict(self.exposed_inputs(SiestaBaseWorkChain))
        siesta_inputs.structure = self.inputs.defect_structure
        siesta_inputs.parameters = self.input.dft.supercell_defect_q.parameters
        siesta_inputs.kpoints = self.inputs.siesta.dft.supercell_defect_q.kpoints
        siesta_inputs.options = self.input.siesta.dft.supercell_defect_q.options
        siesta_inputs.basis = self.input.siesta.dft.supercell_defect_q.basis
        future = self.submit(SiestaBaseWorkChain,**siesta_inputs)
        #future = self.submit(**siesta_inputs)
        self.report(
            'Launching PWSCF for defect structure (PK={}) with charge {} (PK={})'
            .format(self.inputs.defect_structure.pk,
                    self.inputs.defect_charge.value, future.pk))
        self.to_context(**{'calc_defect_q': future})

    def prep_dft_calcs_gaussian_correction(self):
        """
        Get the required inputs for the Gaussian Countercharge correction workchain.
        This method runs the required calculations to generate the energies and potentials
        for the Gaussian scheme.
        """

        self.report("Setting up the Gaussian Countercharge correction workchain")

        siesta_inputs = self.inputs.siesta.dft.supercell.code.get_builder()
        siesta_inputs.pseudos = self.inputs.siesta.dft.supercell.pseudopotentials
        siesta_inputs.kpoints = self.inputs.siesta.dft.supercell.kpoints
#        siesta_inputs.metadata = self.inputs.qe.dft.supercell.scheduler_options.get_dict()

        parameters = self.inputs.siesta.dft.supercell.parameters.get_dict()

        # We set 'tot_charge' later so throw an error if the user tries to set it to avoid
        # any ambiguity or unseen modification of user input
        if 'tot_charge' in parameters['SYSTEM']:
            self.report('You cannot set the "tot_charge" PW.x parameter explicitly')
            return self.exit_codes.ERROR_PARAMETER_OVERRIDE

        # Host structure
        pw_inputs.structure = self.inputs.host_structure
        parameters['SYSTEM']['tot_charge'] = orm.Float(0.)
        pw_inputs.parameters = orm.Dict(dict=parameters)

        future = self.submit(pw_inputs)
        self.report(
            'Launching PWSCF for host structure (PK={}) with charge {} (PK={})'
            .format(self.inputs.host_structure.pk, "0.0", future.pk))
        self.to_context(**{'calc_host': future})

        # Defect structure; neutral charge state
        pw_inputs.structure = self.inputs.defect_structure
        parameters['SYSTEM']['tot_charge'] = orm.Float(0.)
        pw_inputs.parameters = orm.Dict(dict=parameters)

        future = self.submit(pw_inputs)
        self.report(
            'Launching PWSCF for defect structure (PK={}) with charge {} (PK={})'
            .format(self.inputs.defect_structure.pk, "0.0", future.pk))
        self.to_context(**{'calc_defect_q0': future})

        # Defect structure; target charge state
        pw_inputs.structure = self.inputs.defect_structure
        parameters['SYSTEM']['tot_charge'] = self.inputs.defect_charge
        pw_inputs.parameters = orm.Dict(dict=parameters)

        future = self.submit(pw_inputs)
        self.report(
            'Launching PWSCF for defect structure (PK={}) with charge {} (PK={})'
            .format(self.inputs.defect_structure.pk,
                    self.inputs.defect_charge.value, future.pk))
        self.to_context(**{'calc_defect_q': future})


    def check_dft_calcs_gaussian_correction(self):
        """
        Check if the required calculations for the Gaussian Countercharge correction workchain
        have finished correctly.
        """

        # Host
        host_calc = self.ctx['calc_host']
        if host_calc.is_finished_ok:
            self.ctx.host_energy = orm.Float(host_calc.outputs.output_parameters.get_dict()['energy']) # eV
            self.ctx.host_vbm = orm.Float(host_calc.outputs.output_band.get_array('bands')[0][-1]) # valence band maximum
        else:
            self.report(
                'PWSCF for the host structure has failed with status {}'.
                format(host_calc.exit_status))
            return self.exit_codes.ERROR_DFT_CALCULATION_FAILED

        # Defect (q=0)
        defect_q0_calc = self.ctx['calc_defect_q0']
        if not defect_q0_calc.is_finished_ok:
            self.report(
                'PWSCF for the defect structure (with charge 0) has failed with status {}'
                .format(defect_q0_calc.exit_status))
            return self.exit_codes.ERROR_DFT_CALCULATION_FAILED

        # Defect (q=q)
        defect_q_calc = self.ctx['calc_defect_q']
        if defect_q_calc.is_finished_ok:
            self.ctx.defect_energy = orm.Float(defect_q_calc.outputs.output_parameters.get_dict()['energy']) # eV
        else:
            self.report(
                'PWSCF for the defect structure (with charge {}) has failed with status {}'
                .format(self.inputs.defect_charge.value,
                        defect_q_calc.exit_status))
            return self.exit_codes.ERROR_DFT_CALCULATION_FAILED

#    def get_dft_potentials_gaussian_correction(self):
#        """
#        Obtain the electrostatic potentials from the SIESTA calculations.
#        """

#        # User inputs
#        pp_inputs = self.inputs.qe.pp.code.get_builder()
#        pp_inputs.metadata = self.inputs.qe.pp.scheduler_options.get_dict()

#        # Fixed settings
#        pp_inputs.plot_number = orm.Int(11)  # Elctrostatic potential
#        pp_inputs.plot_dimension = orm.Int(3)  # 3D

#        inputs.parent_folder = self.ctx['calc_host'].outputs.remote_folder
#        future = self.submit(pp_inputs)
#        self.report(
#            'Launching SIESTA VT retrieved for host structure (PK={}) with charge {} (PK={})'.
#            format(self.inputs.host_structure.pk, "0.0", future.pk))
#        self.to_context(**{'pp_host': future})

#        pp_inputs.parent_folder = self.ctx[
#            'calc_defect_q0'].outputs.remote_folder
#        future = self.submit(pp_inputs)
#        self.report(
#            'Launching PP.x for defect structure (PK={}) with charge {} (PK={})'
#            .format(self.inputs.defect_structure.pk, "0.0", future.pk))
#        self.to_context(**{'pp_defect_q0': future})
#
#        pp_inputs.parent_folder = self.ctx[
#            'calc_defect_q'].outputs.remote_folder
#        future = self.submit(pp_inputs)
#        self.report(
#            'Launching PP.x for defect structure (PK={}) with charge {} (PK={})'
#            .format(self.inputs.defect_structure.pk,
#                    self.inputs.defect_charge.value, future.pk))
#        self.to_context(**{'pp_defect_q': future})

    def check_dft_potentials_gaussian_correction(self):
        """
        Check if the required calculations for the Gaussian Countercharge correction workchain
        have finished correctly.
        """

        # Host
        host_pp = self.ctx['pp_host']
        if host_pp.is_finished_ok:
            data_array = host_pp.outputs.output_data.get_array('data')
            v_data = orm.ArrayData()
            v_data.set_array('data', data_array)
            self.ctx.v_host = v_data
        else:
            self.report(
                'Post processing for the host structure has failed with status {}'
                .format(host_pp.exit_status))
            return self.exit_codes.ERROR_PP_CALCULATION_FAILED

        # Defect (q=0)
        defect_q0_pp = self.ctx['pp_defect_q0']
        if defect_q0_pp.is_finished_ok:
            data_array = host_pp.outputs.output_data.get_array('data')
            v_data = orm.ArrayData()
            v_data.set_array('data', data_array)
            self.ctx.v_defect_q0 = v_data
        else:
            self.report(
                'Post processing for the defect structure (with charge 0) has failed with status {}'
                .format(defect_q0_pp.exit_status))
            return self.exit_codes.ERROR_PP_CALCULATION_FAILED

        # Defect (q=q)
        defect_q_pp = self.ctx['pp_defect_q']
        if defect_q_pp.is_finished_ok:
            data_array = host_pp.outputs.output_data.get_array('data')
            v_data = orm.ArrayData()
            v_data.set_array('data', data_array)
            self.ctx.v_defect_q = v_data
        else:
            self.report(
                'Post processing for the defect structure (with charge {}) has failed with status {}'
                .format(self.inputs.defect_charge.value,
                        defect_q_pp.exit_status))
            return self.exit_codes.ERROR_PP_CALCULATION_FAILED
