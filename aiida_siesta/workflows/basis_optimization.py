
class BandgapWorkChain(WorkChain):
    """
    Workchain to obtain the bands and bandgap of a structure through Siesta.
    If "bandskpoints" are set in inputs, it behaves like `SiestaBaseWorkChain`
    adding just the bandgap calculation at the end. If no bandskpoints
    was specified, the bands are computed anyway on a kpoints path automatically
    assigned using seekpath and the input (output) structure
    of the single-point (relaxation/md) calculation.
    """

    @classmethod
    def define(cls, spec):
        super().define(spec)
        spec.expose_inputs(SiestaBaseWorkChain, exclude=('metadata','basis'))
        spec.input('basis_sizes', valid_type=orm.List, default=lambda:orm.List(list=["DZ","DZP","TZ"]))
        spec.input('optimization_params', valid_type=orm.Dict, default=lambda: orm.Dict(dict={
                'max_iter': 40,
                'something_for_conv': 0.0001
            }),
            help='dictionary of parameters that are pased to `...`'
        )
        spec.output('optimal_basis_block', valid_type=orm.Dict)
        spec.outline(
            cls.preprocess,
            if_ (no_custom_basis):(
                cls.run_sizes,
                )
            cls.run_optimizer,
            cls.run_results,
        )
        spec.exit_code(200, 'ERROR_SIZES_RUN', message='The SiestaIterator running pao sizes failed')
        spec.exit_code(201, 'ERROR_FINAL_WC', message='The SiestaBaseWorkChain to obtain the bands failed')

    def preprocess(self):
        """
        In the preprocess, we make decisions on bandskpoints if they are not requested in input.
        In case of single point calculation, bandskpoints are added using seekpath.
        In case of relaxation, the relaxation is run, but an extra step at the end of
        the calculation will calculate the bands.
        """
        self.ctx.need_fin_step = False
        var_geom = False
        self.ctx.need_to_generate_bandskp = False

    def no_custom_basis(self):
        return "custom_basis" not in self.inputs

    def run_sizes(self):
        inputs = AttributeDict(self.exposed_inputs(SiestaBaseWorkChain))
        sizes_list = self.inputs.basis_sizes.get_list()
        run_sizes = self.submit(SiestaIterator,
                **inputs, 
                iterate_over={'pao-basis-size': sizes_list},
                batch_size=Int(len(sizes_list))
                )

        return ToContext(sizes_run=run_sizes)

    def set_optimization(self):
        if self.no_custom_basis:
            if not self.ctx.sizes_run.is_finished_ok:
                return self.exit_codes.ERROR_SIZES_RUN

            basis_enthalpy = 100
            good_calc = None
            for calc in self.ctx.sizes_run.called:
                if calc.outputs.output_parameters["basis_enthalpy"] <= basis_enthalpy:
                    basis_enthalpy = calc.outputs.output_parameters["basis_enthalpy"])
                    good_calc = calc
        
            #basis template, max, min
            template_basis = extract_pao_block(good_calc) #maybe some inputs. only first z, only second, spdf,

        

    def run_optimizer(self):
result_wf = run( 
    OptimizationWorkChain, 
    engine=NelderMead, 
    engine_kwargs=Dict(dict=dict(simplex=[[3.2,3.7,3.6,7.5],[5,3.7,3.6,7.5],[3.2,5.5,3.6,7.5],[3.2,3.7,4.8,7.5],[3.2,3.7,3.6,8.7]], input_key="list_of_values", result_key='ene',max_iter=2)), 
    evaluate_process=ForBasisOptWorkChain, 
    evaluate={"upper_boundaries":List(list=[6.5,6.5,6.5,8.7]),"lower_boundaries":List(list=[2,2,3.5,4]),"kpoints":kp,"structure":d.inputs.structure,"basis":basis,"code":load_code(61638),"options":Dict(dict=metadata),"parameters":d.inputs.parameters,"pseudos":{"Sc":d.inputs.pseudos__Sc}} 
)       
