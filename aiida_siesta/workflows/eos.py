# -*- coding: utf-8 -*-

from __future__ import absolute_import
from aiida_siesta.workflows.base import SiestaBaseWorkChain
from aiida.plugins import DataFactory
from aiida.common import AttributeDict
from aiida.engine import WorkChain, calcfunction, ToContext
from aiida.orm import Float

@calcfunction
def scale_to_vol(stru, vol):
    """
    Calcfunction to scale a structure to a target volume.
    Uses pymatgen.
    :param stru: An aiida structure
    :param vol: The target volume per atom in angstroms
    :return: The new scaled AiiDA structure
    """

    in_structure = stru.get_pymatgen()
    new = in_structure.copy()
    new.scale_lattice(float(vol)*in_structure.num_sites)
    StructureData = DataFactory("structure")
    structure  = StructureData(pymatgen_structure=new)

    return structure


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
    new_ase.set_cell(the_ase.get_cell() * float(scale), scale_atoms=True)
    new_structure = DataFactory('structure')(ase=new_ase)
    
    return new_structure


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

    ev = {}
    ev["vol"] = struct.get_cell_volume()/len(struct.sites)
    ev["vol_units"] = 'ang/atom'
    ev["en"] = outpar['E_KS']/len(struct.sites)
    ev["en_units"] = outpar['E_KS_units']+'/atom'
    Dict = DataFactory('dict')
    resultdict = Dict(dict=ev)
    
    return resultdict


def birch_murnaghan(V, E0, V0, B0, B01):
    r = (V0 / V) ** (2. / 3.)
    return E0 + 9. / 16. * B0 * V0 * (r - 1.) ** 2 * \
                (2. + (B01 - 4.) * (r - 1.))


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
    from scipy.optimize import curve_fit
    import numpy as np

    eos = []
    volu = []
    ener = []
    for cal in calcs:
        arg = calcs[cal]
        volu.append(arg["vol"])
        ener.append(arg["en"])
        eos.append([arg["vol"],arg["en"],arg["vol_units"],arg["en_units"]])

    try:
        fit_res = {}
        volumes = np.array(volu)
        energies = np.array(ener)
        params, covariance = curve_fit(
            birch_murnaghan, xdata=volumes, ydata=energies,
            p0=(
                energies.min(),  # E0
                volumes.mean(),  # V0
                0.1,  # B0
                3.,  # B01
            ),
            sigma=None
        )
        fit_res["Eo"] = params[0]
        fit_res["Vo"] = params[1]
        fit_res["Bo"] = params[2]
        fit_res["B1"] = params[3]
        perr = np.sqrt(np.diag(covariance))
        fit_res["EoErr"] = perr[0]
        fit_res["VoErr"] = perr[1]
        fit_res["BoErr"] = perr[2]
        fit_res["B1Err"] = perr[3]
    except:
        fit_res = None

    Dict = DataFactory("dict")
    if fit_res is None:
        result_dict = Dict(dict={'eos_data': eos})
    else:
        result_dict = Dict(dict={'eos_data': eos, "fit_res": fit_res})

    return result_dict


class IsotropicEosFast(WorkChain):
    """
    WorkChain to calculate the isotropic equation of state.
    This means that the cell shape is fixed, only the volume
    is rescaled. In particular the volumes considered are 7
    equidistant volumes around a starting volume.
    The starting volume is an optional input of the WorkChain.
    All the SiestaCalculation parameters are other inputs of the
    workchain.
    """

    @classmethod
    def define(cls, spec):
        super(IsotropicEosFast, cls).define(spec)
        spec.input("volume_per_atom",  valid_type=Float, required=False)
        spec.expose_inputs(SiestaBaseWorkChain)#, exclude=('kpoints',))
        spec.outline(
            cls.initio,
            cls.run_base_wcs,
            cls.return_results
        )
        spec.output('res_dict', valid_type=DataFactory("dict"), required=True)
        spec.output('equilibrium_structure', valid_type=DataFactory("structure"), required=False)

    def initio(self):
        self.ctx.scales = (0.94, 0.96, 0.98, 1., 1.02, 1.04, 1.06)
        self.ctx.expinputs = AttributeDict(self.exposed_inputs(SiestaBaseWorkChain))
        self.report("Starting IsotropicEosfast Workchain")
        if "pseudo_family" not in self.inputs:
            if not self.inputs.pseudos: 
                raise ValueError(
                'neither an explicit pseudos dictionary nor a pseudo_family was specified'
                )
        if "volume_per_atom" in self.inputs:
            self.ctx.s0 = scale_to_vol(self.ctx.expinputs.structure,self.inputs.volume_per_atom)
        else:
            self.ctx.s0 = self.ctx.expinputs.structure

    def run_base_wcs(self):
        calcs = {}
        for scale in self.ctx.scales:
                scaled = rescale(self.ctx.s0, Float(scale))
                inputs = self.ctx.expinputs
                inputs["structure"] = scaled
                future = self.submit(SiestaBaseWorkChain, **inputs)
                self.report('Launching SiestaBaseWorkChain<{}>'.format(future.pk))
                calcs[str(scale)] = future
        return ToContext(**calcs)  #Here it waits

    def return_results(self):
        self.report('All 7 calculations finished. Post process starts')
        collectwcinfo = {}
        for label in self.ctx.scales:
            wcnode = self.ctx[str(label)]
            # To performe the eos with variable cell is wrong. We should 
            # implement errors if there are siesta keywords that activate 
            # the relax cell in the parameter dict. For the moment the
            # use of output_structure should at least make the user realise
            # that something is wrong (we expect to have structure converged 
            # to similar volums if you do relax cell)
            if "output_structure" in wcnode.outputs: 
                info = get_info(wcnode.outputs.output_parameters,wcnode.outputs.output_structure)
            else:
                info = get_info(wcnode.outputs.output_parameters,wcnode.inputs.structure)
            collectwcinfo["s"+str(label).replace(".","_")] = info

        res_dict = fit_and_final_dicts(**collectwcinfo)

        self.out('res_dict', res_dict)

        if "fit_res" in res_dict:
            eq_structure = scale_to_vol(self.ctx.s0, res_dict["fit_res"]["Vo"])
            self.out('equilibrium_structure', eq_structure)


