import aiida
from aiida.engine import WorkChain, ToContext, calcfunction
from aiida.orm import Code, Str, Dict, Float, StructureData
from aiida_siesta.workflows.base import SiestaBaseWorkChain
from aiida_siesta.workflows.iterate import SiestaIterator

import numpy as np

class Stage1Deliverable(WorkChain):
    
    @classmethod
    def define(cls, spec):
        super().define(spec)

        # Define the outline of the workflow, i.e. the order in which methods are executed.
        spec.outline(
            cls.initialize,
            cls.run_calc,
            cls.process_calc
        )

        # We expose the inputs of the SiestaBaseWorkChain so that they can be passed to our workchain
        spec.expose_inputs(SiestaBaseWorkChain, exclude=["metadata"])
        spec.expose_outputs(SiestaBaseWorkChain)
        
        # In case we need to define more inputs, we call spec.input
        spec.input("year")
        
        # Similarly, if we want to define outputs, we do it like this
        spec.output("corrected_E")
    
    def initialize(self):
        """
        Here we will do some initialization/preprocessing if we need it.
        """
        
    def process_inputs(self, inputs):
        """
        This is a chance to process the inputs that will go into the calculation.
        
        It's called before submitting the calculation (see run_calc).
        """
        return inputs
        
    def run_calc(self):
        """
        Where we actually run the calculation
        """
        
        # Get the inputs that we exposed for SiestaBaseWorkChain
        inputs = self.exposed_inputs(SiestaBaseWorkChain)
        
        self.process_inputs(inputs)
        
        # Submit the calculation
        process_node = self.submit(SiestaBaseWorkChain, **inputs)

        # Wait for the process to finish before going to the next step
        # This will store the process into self.ctx
        return ToContext(calc=process_node)
    
    def process_calc(self):
        """
        Here we will process the outputs of the 
        """
        # Retrieve the process that we put to context at the end of run_calc
        process_node = self.ctx.calc
        
        # We return all the outputs from SiestaBaseWorkChain
        self.out_many(self.exposed_outputs(self.ctx.calc, SiestaBaseWorkChain))
        
        # Process the outputs to generate a new node (see process_output docs)
        # https://aiida-siesta-plugin.readthedocs.io/en/latest/plugins/siesta.html#outputs
        corrected_E = correct_energy(process_node.outputs.output_parameters, self.inputs.year)
        
        # When we are ready to output some value, we just need to call self.out 
        self.out("corrected_E", corrected_E)

@calcfunction
def correct_energy(output_parameters, year):
    """
    Processes the output and generates a new node.
    
    In aiida, when you generate a node (i.e. something that you want to use as
    an input or output), you have to do it inside a calcfunction.
    """ 
    return output_parameters["E_KS"] + (year - 2020)

class Stage2Deliverable(Stage1Deliverable):
    
    @classmethod
    def define(cls, spec):
        super().define(spec)

        spec.input("C_substitute", required=False, default=lambda: Str("Si"))
    
    def process_inputs(self, inputs):
        
        if self.inputs.year > 2040:
            inputs["structure"] = substitute_Cs(inputs["structure"], self.inputs.C_substitute)
            
        return inputs
        

@calcfunction
def substitute_Cs(input_struct, substitute):
    
    new_struct = input_struct.get_ase()
    
    new_struct.set_chemical_symbols([sym if sym != "C" else substitute.value for sym in new_struct.symbols])
    
    return StructureData(ase=new_struct)

class Stage3Deliverable(SiestaIterator):
    _process_class = Stage2Deliverable
    _iterate_over_input = False
    
    @classmethod
    def define(cls, spec):
        super().define(spec)
        
        cls.iteration_input("scfmixerweight", required=False, default=[0.25, 0.1, 0.05, 0.01, 0.005, 0.001])
        
        spec.output("fastest_weight")
    
    def initialize(self):
        super().initialize()
        
        self.ctx.times = []
    
    def _analyze_process(self, process_node):
        # Get the time it took for the calculation
        if "output_parameters" in process_node.outputs:
            time = process_node.outputs.output_parameters["global_time"]
        else:
            time = np.inf
        
        self.ctx.times.append(time)
        
    def return_results(self):
        i_max = np.argmin(self.ctx.times)
        
        fastest_weight = self.ctx.used_values[i_max][0]
        
        self.out("fastest_weight", fastest_weight)

class Stage4Deliverable(WorkChain):
    
    @classmethod
    def define(cls, spec):
        super().define(spec)
        
        # Define the outline
        spec.outline(
            cls.find_best_weight,
            cls.relax_structure,
            cls.return_results
        )
        
        # Expose the inputs of the iterator that finds the best weight
        spec.expose_inputs(Stage3Deliverable)
        
        # Expose the outputs of the final relaxation
        spec.expose_outputs(Stage2Deliverable)
        
    def find_best_weight(self):
        """
        Runs the iterator to find the fastest mixing weight.
        """
        # Submit the iterator
        weights_calc = self.submit(Stage3Deliverable, **self.exposed_inputs(Stage3Deliverable))
        
        # And wait for it to end
        return ToContext(weights_calc=weights_calc)
    
    def relax_structure(self):
        """
        Submits a relaxation of the structure using the fastest weight
        """
        # Grab the fastest weight from the outputs of the iterator
        fastest_weight = self.ctx.weights_calc.outputs["fastest_weight"]
        
        # In order to generate inputs for the relaxation, we grab them from the iterator.
        # However, we only need to grab those that are accepted by the relaxation.
        input_keys = set(Stage2Deliverable.spec().inputs).intersection(self.ctx.weights_calc.inputs)
        inputs = {key:self.ctx.weights_calc.inputs[key] for key in input_keys}
        
        # Now, we call a calcfunction to generate the new "parameters" input, given the
        # fastest weight. It will also set mdsteps so that a relaxation is run.
        inputs["parameters"] = setup_cg_parameters(inputs["parameters"], fastest_weight)
        
        # Submit the relaxation
        relaxation_calc = self.submit(Stage2Deliverable, **inputs)
        
        # And wait for it to end
        return ToContext(relaxation_calc=relaxation_calc)
    
    def return_results(self):
        """
        Returns the outputs of the workchain
        """
        # Return all the outputs of the relaxation.
        self.out_many(self.exposed_outputs(self.ctx.relaxation_calc, Stage2Deliverable))
        
@calcfunction
def setup_cg_parameters(old_params, fastest_weight):
    """
    Sets up the parameters input to run a relaxation using the fastest weight.
    """
    # Get the parameters dictionary so that we can modify it
    parameters = old_params.get_dict()
    
    # Set scfmixerweight to the optimal weight (provided)
    parameters["scfmixerweight"] = fastest_weight.value
    
    # Set the mdsteps to a number different than 0 so that it relaxes the structure
    parameters["mdsteps"] = 2000
    
    # Return the newly constructed node.
    return Dict(dict=parameters)