from aiida.orm import (Str, Float, Code, Dict, ArrayData, StructureData)
from aiida.common import AttributeDict
from aiida.engine import WorkChain, calcfunction, ToContext
from aiida_siesta.workflows.base import SiestaBaseWorkChain
from aiida_siesta.calculations.stm import STMCalculation
from aiida_siesta.calculations.tkdict import FDFDict


#This is one of the situations where there is no obbligation to
#use @calcfunction (no error or warning from AiiDA), but it wouldn't
#be bad to use it. We are basically modifing one of the input nodes
#before submitting it to the SiestaBaseWorkChain. The old parameters
#are stored in the input port of SiestaSTMWorkChain, the new once
#in the input port of the SiestaBaseWorkChain, but without @calcfunction,
#the connection between them is lost.
def strip_ldosblock(param):

    param_copy = param.get_dict()
    translated_para = FDFDict(param_copy)
    translated_para.pop('%block localdensityofstates')

    return Dict(dict=translated_para)


#Here, instead, the use of @calcfunction is mandatory as we want to
#create something to return in output
@calcfunction
def create_non_coll_array(**arrays):
    arraydata = ArrayData()
    arraydata.set_array('grid_X', arrays["q"].get_array("grid_X"))
    arraydata.set_array('grid_Y', arrays["q"].get_array("grid_Y"))
    arraydata.set_array('STM_q', arrays["q"].get_array("STM"))
    for spinmod in ("x", "y", "z"):
        arraydata.set_array('STM_s{}'.format(spinmod), arrays[spinmod].get_array("STM"))

    return arraydata


class SiestaSTMWorkChain(WorkChain):
    """
    STM Workchain. This workchain runs a DFT calculation with siesta, calculates
    the local density of states in an energy window specified by the user (stored
    in a .LDOS file) and post-process it in order to produce simulated STM images.
    """

    def __init__(self, *args, **kwargs):
        super(SiestaSTMWorkChain, self).__init__(*args, **kwargs)

    @classmethod
    def define(cls, spec):
        super(SiestaSTMWorkChain, cls).define(spec)
        spec.expose_inputs(SiestaBaseWorkChain, exclude=('metadata',))
        #Temporary fix to issue #135 plumpy
        #spec.inputs._ports['pseudos'].dynamic = True  #pylint: disable=protected-access
        spec.input('emin', valid_type=Float, help='Lower boundary energy (in eV respect to Ef) for LDOS calculation')
        spec.input('emax', valid_type=Float, help='Higher boundary energy (in eV respect to Ef) for LDOS calculation')
        spec.input('stm_code', valid_type=Code, help='STM plstm code')
        spec.input('stm_options', valid_type=Dict, required=False, help='STM plstm code resources and options')
        spec.input('stm_mode', valid_type=Str, help='Allowed values are "constant-height" or "constant-current"')
        spec.input('stm_value', valid_type=Float, help='Value of height in Ang or value of current in e/bohr**3')
        spec.input('stm_spin', valid_type=Str, help='Allowed values are "none", "collinear" or "non-collinear"')
        spec.outline(
            cls.checks,
            cls.run_siesta_wc,
            cls.run_siesta_with_ldos,
            cls.run_stm,
            cls.run_results,
        )

        spec.output('stm_array', valid_type=ArrayData)
        spec.output('output_structure', valid_type=StructureData, required=False)
        #spec.output('output_parameters', valid_type=Dict)

        spec.exit_code(200, 'ERROR_BASE_WC', message='The main SiestaBaseWorkChain failed')
        spec.exit_code(201, 'ERROR_LDOS_WC', message='The SiestaBaseWorkChain to obtain the .LDOS file failed')
        spec.exit_code(202, 'ERROR_STM_PLUGIN', message='The STM post-process failed')

    def checks(self):  # noqa: MC0001  - is mccabe too complex funct -
        """
        Checks on inputs and definition of few variables useful in the next steps
        """

        stm_code = self.inputs.stm_code
        code = self.inputs.code
        mode = self.inputs.stm_mode.value
        spinstm = self.inputs.stm_spin.value
        param_dict = self.inputs.parameters.get_dict()
        translatedkey = FDFDict(param_dict)

        if code.computer.pk != stm_code.computer.pk:
            raise ValueError("The siesta code and the stm code must be on the same computer!")

        allowedmodes = ["constant-height", "constant-current"]
        if mode not in allowedmodes:
            raise ValueError("The allowed options for the port 'stm_mode' are {}".format(allowedmodes))

        #Some logic regarding the spin. Example: the user can request in input
        #of the workchain a stm_spin "non-collinear", but the underline
        #siesta calculation was performed with no spin, a warning must be issued
        #and the stm analysis is carried in non-spin mode. The opposite however
        #is allowed, stm_spin "none" can be generated also from the "non-collinear"
        #calculations. Sintax of both Siesta4.0 and Siesta>4.1 is supported.
        allowedspins = ["none", "collinear", "non-collinear"]
        if spinstm not in allowedspins:
            raise ValueError("The allowed options for the port 'stm_spin' are {}".format(allowedspins))

        spinsiesta = "none"
        newversintax = False
        oldversintax = False

        for k, v in sorted(translatedkey.get_filtered_items()):
            if k == "spinpolarized":
                oldversintax = True
                if v is True or v == "T" or v == "true" or v == ".true.":
                    spinsiesta = "coll"
        if oldversintax:
            for k, v in sorted(translatedkey.get_filtered_items()):
                if k == "noncollinearspin":
                    if v is True or v == "T" or v == "true" or v == ".true.":
                        spinsiesta = "noncoll"
        for k, v in sorted(translatedkey.get_filtered_items()):
            if k == "spin":
                newversintax = True
                translatevalue = FDFDict.translate_key(v)
                if translatevalue in ('noncollinear', 'spinorbit'):
                    spinsiesta = "noncoll"
                if translatevalue == "polarized":
                    spinsiesta = "coll"

        if newversintax and oldversintax:
            self.report(
                'WARNING: in the siesta input parameters, both keywork "spin" and '
                '"spinpolarized" have been detected. This might confuse the WorkChain and return '
                'unexpected outputs'
            )

        if spinsiesta == "none" and spinstm != "none":
            self.report(
                'WARNING: Requested STM with spin option, '
                'but the DFT run was performed with no spin. "stm_spin" is reset to "none"'
            )
            spinstm = "none"
        if spinsiesta == "coll" and spinstm == "non-collinear":
            self.report(
                'WARNING: Requested STM with spin non-collinear option, '
                'but the DFT run was performed with collinear spin. "stm_spin" is reset to "collinear"'
            )
            spinstm = "collinear"

        self.ctx.spinstm = spinstm

        #LDOS check, the inputs "emax" and "emin" define the energy range for the calculation of the
        #ldos. If a block "localdensityofstates" is found in the parameters of the siesta calculation,
        #a warining is issued and the block is stripped.
        self.ctx.ldosdefinedinparam = False
        for k, v in sorted(translatedkey.get_filtered_items()):
            if k == "%block localdensityofstates":
                self.ctx.ldosdefinedinparam = True
                self.report(
                    'WARNING: A local-density-of-state block was defined in the input parameters, however '
                    'this will be ignored as the input ports "emin" and "emax" define the energy of LDOS'
                )

    def run_siesta_wc(self):
        """
        Run the SiestaBaseWorkChain, might be a relaxation or a scf only.
        """

        self.report('Initial checks where succesfull')

        inputs = AttributeDict(self.exposed_inputs(SiestaBaseWorkChain))
        if self.ctx.ldosdefinedinparam:
            newpar = strip_ldosblock(inputs.parameters)
            inputs["parameters"] = newpar

        running = self.submit(SiestaBaseWorkChain, **inputs)
        self.report('Launched SiestaBaseWorkChain<{}> to perform the siesta calculation.'.format(running.pk))

        return ToContext(workchain_base=running)

    def run_siesta_with_ldos(self):
        """
        This step is only necessary because the old versions of siesta do not allow
        to specify the energy range in the local-density-of-states block refered to
        the Fermi energy. Therefore, knowing now the Ef from the previous step,
        we can more effectivly select an energy range for the LDOS file.
        """

        if not self.ctx.workchain_base.is_finished_ok:
            return self.exit_codes.ERROR_BASE_WC

        outwc = self.ctx.workchain_base.outputs

        if "output_structure" in outwc:
            self.report(
                'First Siesta calculation concluded succesfully. In case a restart of the WorkChain is needed, '
                'set node {} as parent_calc_folder and node {} as structure'.format(
                    outwc.remote_folder.pk, outwc.output_structure.pk
                )
            )
        else:
            self.report(
                'First Siesta calculation concluded succesfully. In case a restart of the '
                'WorkChain is needed, set node {} as parent_calc_folder'.format(outwc.remote_folder.pk)
            )

        efermi = outwc.output_parameters["E_Fermi"]
        okemax = self.inputs.emax.value + efermi
        okemin = self.inputs.emin.value + efermi
        restart = self.ctx.workchain_base.get_builder_restart()
        if "output_structure" in outwc:
            restart.structure = outwc.output_structure
        ldos_e = "\n{0:.5f} {1:.5f} eV \n%endblock local-density-of-states".format(okemin, okemax)
        param_dict = restart.parameters.get_dict()
        param_dict["%block local-density-of-states"] = ldos_e
        #pop the relax keys??
        restart.parameters = Dict(dict=param_dict)
        restart.parent_calc_folder = outwc.remote_folder
        settings_dict = {'additional_retrieve_list': ['aiida.BONDS', 'aiida.LDOS']}
        restart.settings = Dict(dict=settings_dict)

        running = self.submit(restart)
        self.report('Launched SiestaBaseWorkChain<{}> to obtain the .LDOS file.'.format(running.pk))

        return ToContext(siesta_ldos=running)

    def run_stm(self):
        """
        Run a STMCalculation with the relaxed_calculation parent folder
        """

        if not self.ctx.siesta_ldos.is_finished_ok:
            return self.exit_codes.ERROR_LDOS_WC

        base_ldos = self.ctx.siesta_ldos
        remote_folder = base_ldos.outputs.remote_folder

        self.report(
            'Finished siesta run to obtain .LDOS file. The remote folder hosting the file '
            'is in the node {}'.format(remote_folder.pk)
        )

        if 'stm_options' in self.inputs:
            optio = self.inputs.stm_options.get_dict()
        else:
            optio = self.inputs.options.get_dict()
            optio["resources"]["num_machines"] = 1
            optio["resources"]["num_mpiprocs_per_machine"] = 1
            optio['withmpi'] = False

        stm_inputs = {
            'value': self.inputs.stm_value,
            'mode': self.inputs.stm_mode,
            'code': self.inputs.stm_code,
            'ldos_folder': remote_folder,
            'metadata': {
                'options': optio
            }
        }

        if self.ctx.spinstm == "non-collinear":
            calcs = {}
            for spinmod in ("q", "x", "y", "z"):
                stm_inputs['spin_option'] = Str(spinmod)
                future = self.submit(STMCalculation, **stm_inputs)
                self.report('launching STMCalculation<{0}> in {1} spin mode'.format(future.pk, spinmod))
                calcs[spinmod] = future
            return ToContext(**calcs)
        if self.ctx.spinstm == "collinear":  #pylint: disable=no-else-return
            stm_inputs['spin_option'] = Str("s")
            running = self.submit(STMCalculation, **stm_inputs)
            self.report('launching STMCalculation<{}> in s spin mode'.format(running.pk))
            return ToContext(stm_calc=running)
        else:
            stm_inputs['spin_option'] = Str("q")
            running = self.submit(STMCalculation, **stm_inputs)
            self.report('launching STMCalculation<{}> in q spin mode'.format(running.pk))
            return ToContext(stm_calc=running)

    def run_results(self):
        """
        Attach the relevant output nodes
        """

        from aiida.engine import ExitCode

        if self.ctx.spinstm == "non-collinear":
            cumarray = {}
            for spinmod in ("q", "x", "y", "z"):
                stmnode = self.ctx[spinmod]
                if not stmnode.is_finished_ok:
                    return self.exit_codes.ERROR_STM_PLUGIN
                cumarray[spinmod] = stmnode.outputs.stm_array
            stm_array = create_non_coll_array(**cumarray)
            self.out('stm_array', stm_array)
        else:
            if not self.ctx.stm_calc.is_finished_ok:
                return self.exit_codes.ERROR_STM_PLUGIN
            stm_plot_calc = self.ctx.stm_calc
            stm_array = stm_plot_calc.outputs.stm_array
            self.out('stm_array', stm_array)

        if self.ctx.workchain_base.outputs.output_parameters.attributes["variable_geometry"]:
            output_structure = self.ctx.siesta_ldos.outputs.output_structure
            ##CAREFULL, if strip relax options ,this must be modified
            self.out('output_structure', output_structure)

        self.report('STM workchain succesfully completed')
        return ExitCode(0)

    @classmethod
    def inputs_generator(cls):  # pylint: disable=no-self-argument,no-self-use
        from aiida_siesta.utils.inputs_generators import StmWorkChainInputsGenerator
        return StmWorkChainInputsGenerator(cls)
