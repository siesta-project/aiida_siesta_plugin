from aiida.plugins import DataFactory
from aiida.engine import calcfunction
from aiida.orm import Float
from aiida_siesta.utils.tkdict import FDFDict
from aiida_siesta.workflows.base import SiestaBaseWorkChain

from ..utils.iterate_absclass import BaseIterator


@calcfunction
def get_info(outpar, struct):
    """
    Calcfunction creating a dictionary with selected
    inputs and results of a SiestaBaseWC, usefull for
    the creation of the Equation of State.
    :param struct: Aiida structure of a SiestaBaseWC
    :param outpar: The output_parameters of a SiestaBaseWC
    :return: A dictionary containing volume per atom and
             energy per atom.
    """

    evdict = {}
    evdict["vol"] = struct.get_cell_volume() / len(struct.sites)
    evdict["vol_units"] = 'ang^3/atom'
    evdict["en"] = outpar['E_KS'] / len(struct.sites)
    evdict["en_units"] = outpar['E_KS_units'] + '/atom'
    Dict = DataFactory('dict')
    resultdict = Dict(dict=evdict)

    return resultdict


def delta_project_BM_fit(volumes, energies):  #pylint: disable=invalid-name
    """
    The fitting procedure implemented in this function
    was copied from the Delta Project Code.
    https://github.com/molmod/DeltaCodesDFT/blob/master/eosfit.py
    It is introduced to fully uniform the delta test procedure
    with the one performed with other codes, moreover it
    has the upside to not use scypi.
    """

    import numpy as np

    #Does the fit always succeed?
    fitdata = np.polyfit(volumes**(-2. / 3.), energies, 3, full=True)
    ssr = fitdata[1]
    sst = np.sum((energies - np.average(energies))**2.)
    residuals0 = ssr / sst
    deriv0 = np.poly1d(fitdata[0])  #pylint: disable=invalid-name
    deriv1 = np.polyder(deriv0, 1)  #pylint: disable=invalid-name
    deriv2 = np.polyder(deriv1, 1)  #pylint: disable=invalid-name
    deriv3 = np.polyder(deriv2, 1)  #pylint: disable=invalid-name

    volume0 = 0
    x = 0
    for x in np.roots(deriv1):
        if x > 0 and deriv2(x) > 0:
            E0 = deriv0(x)  #pylint: disable=invalid-name
            volume0 = x**(-3. / 2.)
            break

    #Here something checking if the fit is good!
    #The choice of residuals0 > 0.01 it is not supported by a real scientific reason, just from experience.
    #Values ~ 0.1 are when fit random numbers, ~ 10^-5 appears for good fits. The check on the presence of
    #a minmum covers the situations when an almost linear dependence is fitted (very far from minimum)
    if volume0 == 0 or residuals0 > 0.01:
        return residuals0, volume0
    derivV2 = 4. / 9. * x**5. * deriv2(x)  #pylint: disable=invalid-name
    derivV3 = (-20. / 9. * x**(13. / 2.) * deriv2(x) - 8. / 27. * x**(15. / 2.) * deriv3(x))  #pylint: disable=invalid-name
    bulk_modulus0 = derivV2 / x**(3. / 2.)
    bulk_deriv0 = -1 - x**(-3. / 2.) * derivV3 / derivV2
    return E0, volume0, bulk_modulus0, bulk_deriv0


@calcfunction
def rescale(structure, scale):
    """
    Calcfunction to rescale a structure by a scaling factor.
    Uses ase.
    :param structure: An AiiDA structure to rescale
    :param scale: The scale factor
    :return: The rescaled structure
    """

    the_ase = structure.get_ase()
    new_ase = the_ase.copy()
    new_ase.set_cell(the_ase.get_cell() * pow(float(scale), 1 / 3), scale_atoms=True)
    new_structure = DataFactory('structure')(ase=new_ase)

    return new_structure


@calcfunction
def scale_to_vol(structure, vol):
    """
    Calcfunction to scale a structure to a target volume.
    Uses ase.
    :param stru: An aiida structure
    :param vol: The target volume per atom in angstroms
    :return: The new scaled AiiDA structure
    """

    in_structure = structure.get_ase()
    new = in_structure.copy()
    vol_ratio = vol * len(in_structure) / in_structure.get_volume()
    new.set_cell(in_structure.get_cell() * pow(vol_ratio, 1 / 3), scale_atoms=True)
    StructureData = DataFactory("structure")
    structure_new = StructureData(ase=new)

    return structure_new


def get_scaled(val, inputs):

    modified_struct = rescale(inputs.structure, val)

    return modified_struct


#def standard_BM_fit(volumes, energies):
#
#    from scipy.optimize import curve_fit
#    import numpy as np
#
#    def birch_murnaghan(V, E0, V0, B0, B01):
#        r = (V0 / V) ** (2. / 3.)
#        return E0 + 9. / 16. * B0 * V0 * (r - 1.) ** 2 * \
#                (2. + (B01 - 4.) * (r - 1.))
#
#    #Does the fit always succeed?
#    params, covariance = curve_fit(
#         birch_murnaghan,
#         xdata=volumes,
#         ydata=energies,
#         p0=(
#            energies.min(),  # E0
#            volumes.mean(),  # V0
#            0.1,  # B0
#            3.,  # B01
#            ),
#         sigma=None
#        )
#
#    #Implement here something checking if there is a reasonable minimum!
#    #...
#    #Implement here something checking if the covariance is small enough!
#    #...
#        return residuals0, volume0
#    else:
#        return params[0], params[1], params[2], params[3]


@calcfunction
def fit_and_final_dicts(**calcs):
    """
    Calcfunction that collects all the E vs V, performs
    the birch_murnaghan fit and creates a dictionary with all
    the relevant results. Uses scipy.optimize and numpy.
    :param clacs: Dictionaries result of get_info
    :return: A dictionary containing a list EvsV and
             the results of the murnagan fit.
    """

    import numpy as np

    eos = []
    volu = []
    ener = []
    for cal in calcs:
        arg = calcs[cal]
        volu.append(arg["vol"])
        ener.append(arg["en"])
        eos.append([arg["vol"], arg["en"], arg["vol_units"], arg["en_units"]])

    volumes = np.array(volu)
    energies = np.array(ener)
    try:
        #pylint: disable=invalid-name,unbalanced-tuple-unpacking
        E0, volume0, bulk_modulus0, bulk_deriv0 = delta_project_BM_fit(volumes, energies)
        #E0, volume0, bulk_modulus0, bulk_deriv0 = standard_BM_fit(volumes,energies)
        fit_res = {}
        fit_res["Eo(eV/atom)"] = E0
        fit_res["Vo(ang^3/atom)"] = volume0
        fit_res["Bo(eV/ang^3)"] = bulk_modulus0
        fit_res["Bo(GPa)"] = bulk_modulus0 * 160.21766208
        fit_res["B1"] = bulk_deriv0
    except:  # nopep8 #pylint: disable=bare-except
        fit_res = {}
        #residuals0, volume0 = delta_project_BM_fit(volumes, energies)
        #In the future we could use these info to improve help,
        #residuals0 is a np array

    Dict = DataFactory("dict")
    if fit_res:
        result_dict = Dict(dict={'eos_data': eos, "fit_res": fit_res})
    else:
        result_dict = Dict(dict={'eos_data': eos})

    return result_dict


class EqOfStateFixedCellShape(BaseIterator):
    """
    WorkChain to calculate the equation of state of a solid.
    The cell shape is fixed, only the volume is rescaled.
    In particular the volumes considered are 7 equidistant volumes
    around a starting volume. The starting volume is
    an optional input of the WorkChain (called volume_per_atom).
    If not specified, the input structure volume is used with no modifications.
    All the SiestaBaseWorkChain inputs are other inputs of the workchain.
    This WorkChain also tries to perform a Birch_Murnaghan fit
    on the calculatad E(V) data.
    """

    _process_class = SiestaBaseWorkChain
    # We remove the iterate_over port because we actually only want to expose
    # one kind of iteration: the structure scales.
    _iterate_over_port = False

    @classmethod
    def define(cls, spec):
        super().define(spec)

        spec.input(
            "volume_per_atom",
            valid_type=Float,
            required=False,
            help="Volume per atom around which to perform the EqOfState"
        )

        cls.iteration_input(
            "scales",
            default=(0.94, 0.96, 0.98, 1., 1.02, 1.04, 1.06),
            parse_func=get_scaled,
            input_key="structure",
            help="""
            Factors by which the structure should be scaled.
            """,
        )

        spec.output(
            'results_dict',
            valid_type=DataFactory("dict"),
            required=True,
            help="Containing the calculated E(V) data and, if the fit is sucessfull, "
            "the resulting fit parameters"
        )
        spec.output(
            'equilibrium_structure',
            valid_type=DataFactory("structure"),
            required=False,
            help="Equilibrium volume structure. Returned only if the fit is succesfull"
        )

    def initialize(self):
        super().initialize()

        self.ctx.collectwcinfo = []

        # We are going to overwrite the initial structure if volume_per_atom is provided
        if "volume_per_atom" in self.ctx.inputs:
            self.ctx.inputs.structure = scale_to_vol(self.ctx.inputs.structure, self.ctx.inputs.volume_per_atom)

        test_input_params = FDFDict(self.ctx.inputs.parameters.get_dict())
        for k, v in sorted(test_input_params.get_filtered_items()):
            if k in ('mdvariablecell', 'mdrelaxcellonly'):
                if v is True or v == "T" or v == "true" or v == ".true.":
                    self.report(
                        'WARNING: Relaxation with variable cell detected! '
                        'No action taken, but are you sure this is what you want?'
                    )

    def _analyze_process(self, process_node):

        if "output_structure" in process_node.outputs:
            out_struct = process_node.outputs.output_structure
        else:
            out_struct = process_node.inputs.structure

        info = get_info(process_node.outputs.output_parameters, out_struct)

        self.ctx.collectwcinfo.append(info)

    def return_results(self):

        from aiida.engine import ExitCode

        collectwcinfo = {
            f"s{scale.value}".replace(".", "_"): info
            for (scale,), info in zip(self.ctx.used_values, self.ctx.collectwcinfo)
        }

        res_dict = fit_and_final_dicts(**collectwcinfo)

        self.out('results_dict', res_dict)

        if "fit_res" in res_dict.attributes:
            self.report('Birch-Murnaghan fit was succesfull, creating the equilibrium structure output node')
            eq_structure = scale_to_vol(self.ctx.inputs.structure, Float(res_dict["fit_res"]["Vo(ang^3/atom)"]))
            self.out('equilibrium_structure', eq_structure)
        else:
            self.report("WARNING: Birch-Murnaghan fit failed, check your results_dict['eos_data']")

        return ExitCode(0)

    @classmethod
    def inputs_generator(cls):  # pylint: disable=no-self-argument,no-self-use
        from aiida_siesta.utils.protocols_system.input_generators import EosWorkChainInputGenerator
        return EosWorkChainInputGenerator(cls)
