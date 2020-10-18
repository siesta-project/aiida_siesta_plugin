from aiida.engine import run_get_node, WorkChain
from aiida.orm import StructureData, Int, Str   
from aiida_siesta.workflows.base import SiestaBaseWorkChain
from aiida_siesta.workflows.iterate import SiestaIterator

import ase
from ase.build import molecule, graphene

import inspect
import importlib
from pathlib import Path

GENERAL_INPUTS = {}

def check_class(stage, deliverable):

    if stage in [1,2,4]:
        assert issubclass(deliverable, WorkChain), ("You delivered something that is not a workchain. Please"
            " create your workflow inheriting from `WorkChain`: `class MyWorkflow(WorkChain):`"
        )
    elif stage == 3:
        assert issubclass(deliverable, SiestaIterator), ("Your iterator does not inherit from `SiestaIterator`." 
            " Please create your workflow inheriting from `SiestaIterator`: `class MyIterator(SiestaIterator):`"
            " so that you can benefit from its functionalities."
        )

def check_inputs(stage, deliverable):

    siestabase_inps = SiestaBaseWorkChain.spec().inputs
    deliverable_inps = deliverable.spec().inputs
    
    if stage == 1:
        missing_inps = set(siestabase_inps) - set(deliverable_inps)
        assert not missing_inps, (f"Some settings of SiestaBaseWorkChain are not exposed: {missing_inps}."
            "The user must be able to use these inputs!")
        assert "year" in deliverable_inps, ("There is no input called 'year', how will the user pass the year"
            " they want the energy for?")

        year_dtype = deliverable_inps["year"].valid_type
        if year_dtype:
            assert issubclass(Int, year_dtype), "The year input does not accept an integer?"
    elif stage == 2:
        check_inputs(1, deliverable)

        C_substitute = "C_substitute"
        assert C_substitute in deliverable_inps, (f"There is no input called '{C_substitute}'. How are the users"
            " going to specify the substituting atom if they want?")
        assert deliverable_inps[C_substitute].required == False, ("The '{C_substitute}' input should not be required."
            " You can change that by passing required=False to the input definition")

        default = deliverable_inps[C_substitute].default
        assert default, ("The '{C_substitute}' input does not have a default value.",
            " You can set one by using the 'default' keyword on input definition.")
        assert callable(default) and isinstance(default(), Str) and default() == "Si", (f"The default value for '{C_substitute}'"
            " should be 'Si'. However, aiida likes that the default is produced each time the workchain is called"
            ", so you can set 'default' to a function that returns Str('Si'), e.g.: `lambda: Str('Si')`")
    elif stage == 3:

        scfmixerweight = "scfmixerweight"
        assert scfmixerweight in deliverable_inps, (f"There is no input called '{scfmixerweight}'. How are the users"
            " going to specify what mixing weights to try?")
        assert scfmixerweight in deliverable._iteration_inputs, (f"'{scfmixerweight}' is not an iteration input, please"
            " use cls.iteration_input() instead of spec.input() inside define to inform that it is an iteration input"
            " so that all functionality can be added"
        )
        assert not deliverable_inps[scfmixerweight].required, f"It would be nice that {scfmixerweight} was not required."
        assert deliverable_inps[scfmixerweight].default is not None, (f"'{scfmixerweight}' should have some default values"
            " so that the user doesn't need to provide them. The default passed to cls.iteration_input should be a python list,"
            " don't worry about using aiida types here!"
        )

def check_outputs(stage, deliverable):
    deliverable_outs= deliverable.spec().outputs

    if stage in [1,2]:
        assert "corrected_E" in deliverable_outs, ("We want the corrected energy in the 'corrected_E' output"
            ", but we didn't find such an output in your workchain.")
    if stage == 3:
        assert "fastest_weight" in deliverable_outs, ("We need to have the fastest weight provided as the 'fastest_weight'"
            "output. Is it too much to ask?!!"
        )

def prepare_run(stage, deliverable, inputs):
    inputs["structure"] = StructureData(ase=graphene(vacuum=15))
    inputs["year"] = Int(2050)

    if stage == 2:
        inputs["C_substitute"] = Str("N")
    elif stage in [3, 4]:
        inputs["C_substitute"] = Str("H")
        inputs["max_iterations"] = Int(1)
        inputs["batch_size"] = Int(3)
    
def run(stage, deliverable, inputs):
    return run_get_node(deliverable, **inputs)

def check_run(stage, calc_outputs, calc_node, deliverable, inputs):

    if stage in [1,2]:

        siesta_calc = None
        for called in calc_node.called:
            if issubclass(called.process_class, SiestaBaseWorkChain):
                siesta_calc = called
        
        assert siesta_calc is not None, ("We didn't find a call to SiestaBaseWorkChain, do you call it"
            " inside your workflow?")
        
        uncorrected_E = siesta_calc.outputs.output_parameters["E_KS"]

        expected_E = uncorrected_E + (inputs["year"].value - 2020)
        output_E = calc_outputs["corrected_E"].value

        assert abs(output_E - expected_E) < 0.01, ("The energy was not properly corrected."
            f" We expected {uncorrected_E} + ({inputs['year'].value} - {2020}) = {expected_E} eV, but got {output_E} eV instead."
            " Did you use E_KS as the energy from siesta's output parameters?")
    if stage == 2:

        input_struct = inputs["structure"].get_ase()
        used_struct = siesta_calc.inputs["structure"].get_ase()

        assert list(input_struct.symbols) != list(used_struct.symbols), ("Carbons were not substituted and a simulation"
            f" with C was performed in year {inputs['year'].value}. How useless!")
        assert [inputs["C_substitute"].value if sym == "C" else sym for sym in input_struct.symbols] == list(used_struct.symbols), (
            f"All carbon atoms have not been properly substituted. The input structure contained {input_struct.symbols}"
            f" and it has been converted to {used_struct.symbols}"
        )
    if stage == 3:
        pass
    if stage == 4:
        pass


def deliver_stage(stage, deliverable, **kwargs):
    """Runs all tests on a deliverable to check if it works as expected"""
    inputs = {**GENERAL_INPUTS, **kwargs}

    check_class(stage, deliverable)

    check_inputs(stage, deliverable)

    check_outputs(stage, deliverable)

    prepare_run(stage, deliverable, inputs)

    calc_outputs, calc_node = run(stage, deliverable, inputs)

    check_run(stage, calc_outputs, calc_node, deliverable, inputs)
  
    print(f"ALL TESTS PASSED! Stage {stage} completed succesfully.")

def stage_solution(stage):

    solutions = importlib.import_module("_solutions")
    
    deliverable = getattr(solutions, f"Stage{stage}Deliverable")

    return Solution(deliverable)

class Solution:

    def __init__(self, deliverable):
        self.deliverable = deliverable
        self._code = inspect.getsource(deliverable)
    
    def __str__(self):
        return self._code
    
    __repr__ = __str__
    
    def _repr_markdown_(self):
        return f"```python\n{self._code}\n```"