from aiida.orm import (Float, Str, ArrayData, Dict, StructureData)
from aiida.engine import submit
from aiida.engine import WorkChain, ToContext, calcfunction
from aiida_siesta.workflows.base import SiestaBaseWorkChain
from aiida_siesta.workflows.functions.relaxinputs import SiestaRelaxationInputsGenerator

@calcfunction
def get_energy(pardict):
    return Float(pardict["E_KS"])

@calcfunction
def get_for_stre(ForcAndStress):
    forces=ArrayData()
    forces.set_array(name="forces",array=ForcAndStress.get_array("forces"))
    stress=ArrayData()
    stress.set_array(name="stress",array=ForcAndStress.get_array("stress"))
    return {"forces" : forces,"stress" : stress}

class SiestaRelaxWorkChainProtocol(WorkChain):
    """
    Workchain to relax a structure through Siesta. It make use of 
    protocols, the inputs are generated through the method get_builder
    of the class  SiestaRelaxationInputsGenerator.
    """
    def __init__(self, *args, **kwargs):
        super(SiestaRelaxWorkChainProtocol, self).__init__(*args, **kwargs)

    @classmethod
    def define(cls, spec):
        super(SiestaRelaxWorkChainProtocol, cls).define(spec)
        spec.input('calc_engines', valid_type=Dict)
        spec.input('structure', valid_type=StructureData)
        spec.input('protocol', valid_type=Str, default=Str('standard'))
        spec.input('relaxation_type', valid_type=Str)
        spec.input('threshold_forces', required=False)
        spec.input('threshold_stress',required=False)
        spec.outline(
            cls.setup_protocol,
        #    cls.run_relax,
            cls.run_results,
        )
        # These the standard outouts agreed with other plugins
        spec.output('relaxed_structure',valid_type=StructureData)
        spec.output('forces', valid_type=ArrayData)
        spec.output('stress', valid_type=ArrayData)
        spec.output('total_energy', valid_type=Float, required=False)

    def setup_protocol(self):
        self.report("Setup protocol")

        if "threshold_forces" not in self.inputs:
            threshold_forces = None
        else:
            threshold_forces = self.inputs.threshold_forces

        if "threshold_stress" not in self.inputs:
            threshold_stress = None
        else:
            threshold_stress = self.inputs.threshold_stress

        instgen = SiestaRelaxationInputsGenerator()
        builder=instgen.get_builder(
                    structure=self.inputs.structure,
                    protocol=self.inputs.protocol.value,
                    relaxation_type=self.inputs.relaxation_type,
                    calc_engines=self.inputs.calc_engines,
                    threshold_forces=threshold_forces,
                    threshold_stress=threshold_stress
                    )
        
#        self.ctx.builder=builder

#        return

#    def run_relax(self):
        self.report("Run SiestaWC")
        future = self.submit(builder)
        return ToContext(calc=future)

    def run_results(self):
        self.report("Set outputs")
        self.out('relaxed_structure', self.ctx.calc.outputs.output_structure)
        self.out('total_energy', get_energy(self.ctx.calc.outputs.output_parameters))
        res_dict=get_for_stre(self.ctx.calc.outputs.forces_and_stress)
        self.out('forces', res_dict["forces"])
        self.out('stress', res_dict["stress"])

