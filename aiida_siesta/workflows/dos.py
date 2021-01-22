from aiida import orm
from aiida.engine import WorkChain, calcfunction, ToContext
from aiida.common import AttributeDict
from aiida_siesta.workflows.base import SiestaBaseWorkChain
from aiida_siesta.calculations.tkdict import FDFDict
from aiida_siesta.data.eig import EigData


@calcfunction
def get_dos(eig, dos_inp):
    """
    From an EigData instance, calles the `compute_dos` method
    and puts the results in an ArrayData.

    :param eig: the EigData instance
    :param dos_inp: a dictionary with the values to be unfolded as arguments of `compute_dos`
    :return array_data: ArratData containind the array with the energies (energies) and the one with
        dos (dos, a nspin x len(energies) array).
    """
    dos = eig.compute_dos(**dos_inp)
    print(dos)
    array_data = orm.ArrayData()
    array_data.set_array("energies", dos[0])
    array_data.set_array("dos", dos[1])
    return array_data


def drop_md_keys(param):
    """
    Drop the molecular dinamics keyword from the parameters dictionary.

    :param param: dict containing the parameters of the calculation
    :retuen param: dict where the md keywords have been removed
    """
    for item in param.copy().keys():
        trans_item = FDFDict.translate_key(item)
        if trans_item.startswith("md"):
            param.pop(item)
    return param


def set_dos_defaults(dos_inp, eig):
    """
    Method to set the defaults for arguments that have not been specified in input.
    Note that validation is already performed in input. Here we just fix the case when
    `e_max` has been set smaller than `e_min`.

    :param dos_inp: the dictionary passed as `dos_arguments` input (or an empty dict)
    :param eig: the EigData instance obtained from the main Siesta calculation.
    :return new_dict: a new dictionary containing the arguments for `compute_dos`
        where the defaults are set.
    """

    eigs = eig.get_eigs()

    d_ene = dos_inp.get("d_ene", None)
    if d_ene is None:
        d_ene = 0.005

    e_max = dos_inp.get("e_max", None)
    if e_max is None:
        e_max = eigs.max()

    e_min = dos_inp.get("e_min", None)
    if e_min is None:
        e_min = eigs.min()

    if e_min > e_max:
        e_min, e_max = e_max, e_min

    smearing = dos_inp.get("smearing", None)
    if smearing is None:
        smearing = 0.05

    new_dict = {"e_max": e_max, "e_min": e_min, "smearing": smearing, "d_ene": d_ene}

    return new_dict


def check_eigs(eig, dos_dict, kp_in):
    """
    Analyse the avarage distance between energies and takes action based on the
    current used kpoints mesh and the selected smearing for the future calculation of the dos.
    In particular, sets a threshold of "minimum required avarage separation" based on the
    smearing. calculates the separation in the eigenvalues in EigData and, if the
    separation is bigger that the threshold, creates a ne kpoints mesh to improve the
    separation. This new mesh will be used in a supplementary calculation.

    :param eig: the EigData instance hosting the eigenvalues.
    :param dos_dict: the arguments for `compute_dos of EigData, where, in particular,
        yhe smearing is stored.
    :param kp_in: the kpoints data (containing the mesh), used for the sista calculation
    :return need_cal: bool signaling whether or not a calculation with improved kp mesh is needed
    :return kpoints: KpointsData vontaining theimproved kpoints mesh.
    """

    from aiida.orm import KpointsData

    is_insulator, max_val, min_cond = get_gap(eig)

    thresh_ins = dos_dict["smearing"]
    thresh = dos_dict["smearing"] / 2

    eigs = eig.get_eigs()

    if is_insulator:
        count_cond = 0
        count_val = 0
    else:
        count = 0

    for spin in eigs:  #pylint: disable=too-many-nested-blocks
        for kpoint in spin:
            for ene in kpoint:
                if dos_dict["e_max"] > ene > dos_dict["e_min"]:
                    if is_insulator:
                        if ene < eig.e_fermi:
                            count_val = count_val + 1
                        else:
                            count_cond = count_cond + 1
                    else:
                        count = count + 1

    if is_insulator:
        mean_separation_cond = (dos_dict["e_max"] - min_cond) / count_cond
        mean_separation_val = -(dos_dict["e_min"] - max_val) / count_val
    else:
        mean_separation = (dos_dict["e_max"] - dos_dict["e_min"]) / count

    need_cal = False
    kpoints = KpointsData()
    if is_insulator:
        if mean_separation_cond > thresh_ins or mean_separation_val > thresh_ins:
            mean_separation = max(mean_separation_cond, mean_separation_val)
            need_cal = True
            new_mesh = [int(old * mean_separation / thresh_ins) for old in kp_in.get_kpoints_mesh()[0]]
            kpoints.set_kpoints_mesh(new_mesh)
    else:
        if mean_separation > thresh:
            need_cal = True
            new_mesh = [int(old * mean_separation / thresh) for old in kp_in.get_kpoints_mesh()[0]]
            kpoints.set_kpoints_mesh(new_mesh)

    return need_cal, kpoints


def get_gap(eig):
    """
    Return the band gap, from the eigenvalues (not so reliable?)

    :param eig: the EigData containing the eigenvalues
    :return is_insulator: bool signaling insulator vs metal
    :return val: valence level
    :return cond: conduction level
    """
    eigs = eig.get_eigs()
    cond = 200000
    val = -200000

    for spin in eigs:
        for kpoints in spin:
            for ene in kpoints:
                if val < ene <= eig.e_fermi:
                    val = ene
                if cond > ene >= eig.e_fermi:
                    cond = ene

    is_insulator = cond - val > 0.001

    return is_insulator, val, cond


class DosWorkChain(WorkChain):
    """
    Workchain to obtain the dos of a structure through Siesta.
    In situations where the kpoints mesh used for the calculation is considered too small
    (see `check_eigs` method), a further calculation with improved kpooints mesh is
    performed before returning the dos data.
    """

    @classmethod
    def define(cls, spec):
        super().define(spec)
        spec.expose_inputs(SiestaBaseWorkChain, exclude=('metadata',))
        spec.input_namespace(
            'dos_arguments',
            dynamic=True,
            required=False,
            help='The inputs that will be passed to `comput_dos` method.'
        )
        spec.input('dos_arguments.smearing', valid_type=(float, int), non_db=True, required=False)
        spec.input('dos_arguments.e_max', valid_type=(float, int), non_db=True, required=False)
        spec.input('dos_arguments.e_min', valid_type=(float, int), non_db=True, required=False)
        spec.input('dos_arguments.d_ene', valid_type=(float, int), non_db=True, required=False)
        spec.input('dos_arguments.distribution', valid_type=str, non_db=True, required=False)

        spec.output('dos_array', valid_type=orm.ArrayData)
        spec.output('eig', valid_type=EigData)

        spec.outline(
            #        cls.preprocess,
            cls.run_siesta_wc,
            cls.run_last,
            cls.run_results,
        )
        spec.exit_code(200, 'ERROR_MAIN_WC', message='The main SiestaBaseWorkChain failed')
        spec.exit_code(201, 'ERROR_FINAL_WC', message='The SiestaBaseWorkChain to obtain the bands failed')

    def run_siesta_wc(self):
        """
        Run the SiestaBaseWorkChain, might be a relaxation or a scf only.
        """

        self.report('Initial checks where succesfull')

        inputs = AttributeDict(self.exposed_inputs(SiestaBaseWorkChain))

        running = self.submit(SiestaBaseWorkChain, **inputs)
        self.report(f'Launched SiestaBaseWorkChain<{running.pk}> to perform the siesta calculation.')

        return ToContext(workchain_base=running)

    def run_last(self):
        """
        .
        """
        base_wc = self.ctx.workchain_base

        if not base_wc.is_finished_ok:
            return self.exit_codes.ERROR_MAIN_WC

        if "dos_arguments" in self.inputs:
            inp_dos_args = self.inputs.dos_arguments
        else:
            inp_dos_args = {}
        dos_dict = set_dos_defaults(inp_dos_args, base_wc.outputs.eig)

        self.ctx.dos_dict = dos_dict

        new_calc_needed, kpoints = check_eigs(base_wc.outputs.eig, dos_dict, base_wc.inputs.kpoints)

        self.ctx.need_fin_step = new_calc_needed

        if new_calc_needed:
            new_calc = self.ctx.workchain_base.get_builder_restart()
            if "output_structure" in self.ctx.workchain_base.outputs:
                new_calc.structure = self.ctx.workchain_base.outputs.output_structure
                new_param = drop_md_keys(new_calc.parameters.get_dict())
                new_calc.parameters = orm.Dict(dict=new_param)
            new_calc.parent_calc_folder = self.ctx.workchain_base.outputs.remote_folder
            new_calc.kpoints = kpoints
            running = self.submit(new_calc)
            self.report(f'Launched SiestaBaseWorkChain<{running.pk}> with denser kpoints mesh.')
            return ToContext(final_run=running)

    def run_results(self):
        if self.ctx.need_fin_step:
            if not self.ctx.final_run.is_finished_ok:
                return self.exit_codes.ERROR_FINAL_WC
            outps = self.ctx.final_run.outputs
        else:
            outps = self.ctx.workchain_base.outputs

        dos_array = get_dos(outps["eig"], orm.Dict(dict=self.ctx.dos_dict))

        self.out('eig', outps['eig'])
        self.out('dos_array', dos_array)


#    @classmethod
#    def inputs_generator(cls):  # pylint: disable=no-self-argument,no-self-use
#        from aiida_siesta.utils.inputs_generators import BaseWorkChainInputsGenerator
#        return BaseWorkChainInputsGenerator(cls)
