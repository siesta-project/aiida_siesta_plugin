# -*- coding: utf-8 -*-
########################################################################################
# Copyright (c), The AiiDA-SIESTA-LUA-NEB authors. All rights reserved.                #
#                                                                                      #
# AiiDA-SIESTA-LUA-Defects is hosted on GitHub at https:                               #
# For further information on the license, see the LICENSE.txt file                     #
########################################################################################
from __future__ import absolute_import

from aiida import orm
from aiida.engine import WorkChain, calcfunction, ToContext, if_, submit

#from .utils import (get_barrier_energy)


class BarrierEnergyWorkchainBase(WorkChain):
    """
    The base class to compute the Barrier energy for a given initial & final defect, containing the 
    generic, code-agnostic methods, error codes, etc.

    """
    @classmethod
    def define(cls, spec):
        super(BarrierEnergyWorkchainBase, cls).define(spec)
        # fmt: off
        # Structures 
        spec.input(
            "host_structure", 
            valid_type=orm.StructureData, 
            help="Pristine structure",
            required=True
        )
        spec.input(
            "initial_defect_structure", 
            valid_type=orm.StructureData, 
            help="Initial Defective structure",
            required=False

        )
        spec.input(
            "final_defect_structure",
            valid_type=orm.StructureData,
            help="Final Defective structure",
            required=False
        )

        # Defect details
        spec.input(
            "defect_charge", 
            valid_type=orm.Float, 
            help="Defect charge state",
            required=False
            )
        spec.input(
            "initial_defect_site",
            valid_type=orm.List,
            help="Initial Defect site position in crystal coordinates",
            required=True
        )
        spec.input(
            "final_defect_site",
            valid_type=orm.List,
            help="Final Defect site position in crystal coordinates",
            required=True
        )

        spec.input(
            "ghost",
             valid_type=orm.Str,
             default=orm.Str("yes"),
             help="Either using Ghost Atoms or not (Particular to SIESTA)",
            required=False,
                )
        # Methodology
        spec.input(
            "barrier_scheme",
            valid_type=orm.Str,
            help="The Scheme to Apply for Barrier",
            required=True
        )
        spec.input(
            "neb_scheme",
            valid_type=orm.Str,
            help="The NEB Scheme For Calculating Barrier Energy (Particular to SIESTA/LUA)",
            required=True

        )
        spec.input(
            "number_of_image",
            valid_type=orm.Int,
            default=orm.Int(6),
            help="The Number of Image to Apply for Image Generatation",
            required=True

        )
        spec.input(
            "spring_constant",
            valid_type=orm.Float,
            default=orm.Float(0.5),
            help="The spirng constant value (k) to Apply for NEB",
            required=True
        )
        spec.input(
            "image_interpolation_method",
            valid_type=orm.Str,
            default=orm.Str("idpp"),
            help="The Method for image generation Either (LI) Linear or (IDPP) Image dependetn pair potential to Apply for NEB",
            required=True
        )



        # Outputs
        spec.output(
                "barrier_energy", 
                valid_type=orm.Float, 
                 )

        # Error codes
        spec.exit_code( 401, "ERROR_INVALID_BARRIER_SCHEME",
            message="The requested correction scheme is not recognised",
        )
        spec.exit_code(402,"ERROR_INVALID_NEB_SCHEME",
            message="The requested neb scheme is not recognised",
        )
        spec.exit_code(403, "ERROR_BARRIER_WORKCHAIN_FAILED",
            message="The barrier scheme sub-workchain failed",
        )

        spec.exit_code(404, "ERROR_DFT_CALCULATION_FAILED", 
            message="DFT calculation failed",
        )
        spec.exit_code(405, "ERROR_PP_CALCULATION_FAILED",
            message="A post-processing calculation failed",
        )
        spec.exit_code(406,"ERROR_Pre_PROCESSING_CALCULATION_FAILED",
            message="A pre-processing calculation failed"
        )
        spec.exit_code(500, "ERROR_PARAMETER_OVERRIDE",
            message="Input parameter dictionary key cannot be set explicitly",
        )
        spec.exit_code(999, "ERROR_NOT_IMPLEMENTED",
            message="The requested method is not yet implemented",
        )
        # fmt: on

    def setup(self):
        """ 
        Setup the NEB workchain
        """

        # Check if barrier scheme is valid:
        self.report("Checking Barrier Scheme")
        barrier_schemes_available = ["vacancy-exchange","exchange","kick-out","ring","substitutional"]
        if self.inputs.barrier_scheme is not None:
            if self.inputs.barrier_scheme not in barrier_schemes_available:
                return self.exit_codes.ERROR_INVALID_BARRIER_SCHEME
        # Check if image scheme is valid:
        self.report("Checking NEB Scheme")
        neb_schemes_available = ["neb","dneb","vc-neb"]
        if self.inputs.neb_scheme is not None:
            if self.inputs.neb_scheme not in neb_schemes_available:
                return self.exit_codes.ERROR_INVALID_BARRIER_SCHEME


    def barrier_requested(self):
        """
        Check if barrier is requested
        """
        #if self.inputs.correction_scheme is not None:
        if self.inputs.barrier_scheme=="vacancy-exchange" or self.inputs.barrier_scheme=="exchange" or self.inputs.barrier_scheme=="kick-out" or self.inputs.barrier_scheme=="ring" or self.inputs.barrier_scheme=="substitutional": 
            self.report("There Barrier Scheme will be "+ str(self.inputs.barrier_scheme) )
            return True
        else:
            self.report("There will be no Corrections applied")
            return False

    #----------------------------------------------------------
    # Checking which type of Barrier want to calculate
    #----------------------------------------------------------
    def vacancy_exchange_barrier(self):
        """
        Checking exchange barrier
        """
        return self.inputs.barrier_scheme=="vacancy_exchange"

    def exchange_barrier(self):
        """
        Checking exchange barrier
        """
        return self.inputs.barrier_scheme=="exchange"

    def kick_out_barrier(self):
        """
        Checking kick out barrier
        """
        return self.inputs.barrier_scheme=="kick-out"

    def ring_barrier(self):
        """
        Checking exchange barrier
        """
        return self.inputs.barrier_scheme=="ring"

    def substitutional_barrier(self):
        """
        Checking exchange barrier
        """
        return self.inputs.barrier_scheme=="substitutional"


  
    def vacancy_exchange_image_genration_workchain(self):
        """
        Run the workchain for Unrelaxed vacancy exchange image Generation  
        """
        from .image_generation  import (ImageGenerationWorkChain)

        self.report("Computing unrelaxed images for vacancy exchange barrier scheme")

        inputs = {
            "number_of_image ": self.inputs.number_of_image,
            "image_interpolation_method": self.inputs.image_interpolation_method,
            "initial_defect_site":self.inputs.initial_defect_structure,
            "final_defect_site":self.inputs.final_defect_structure,
            "host_structure": self.inputs.host_structure,
        }

        workchain_future = self.submit(ImageGenerationWorkChain, **inputs)
        label = "image_generation_workchain"
        self.to_context(**{label: workchain_future})

    def relaxed_vacancy_exchange_image_generation_workchain(self):
        """
        Run the workchain for Relaxed vacancy exchange image Generation  
        """
        from .image_generation  import (ImageGenerationWorkChain)

        self.report("Computing relaxed image for vacancy exchange barrier scheme")

        inputs = {
            "number_of_image ": self.inputs.number_of_image,
            "image_generation_type": self.inputs.image_generation_type,
            "initial_defect_site":self.inputs.initial_defect_site,
            "final_defect_site":self.inputs.final_defect_site,
            
            "host_structure": self.inputs.host_structure,
        }

        workchain_future = self.submit(ImageGenerationWorkChain, **inputs)
        label = "image_generation_workchain"
        self.to_context(**{label: workchain_future})   

    
    #============================================
    # Compute Charged Formation Energy without Corrections
    #============================================    
#    def compute_barrier_energy(self):
#        """ 
#        Compute the Barrier Energy 
#        """
#        # Raw formation energy
#        self.ctx.barrier = get_barrier_energy(
#            self.ctx.barrier_energy
#        )
#        self.report(
#            "The computed barrier energy  is {} eV".format(
#                self.ctx.barrier.value
#            )
#        )

    def raise_not_implemented(self):
        """
        Raise a not-implemented error
        """
        return self.exit_codes.ERROR_NOT_IMPLEMENTED
