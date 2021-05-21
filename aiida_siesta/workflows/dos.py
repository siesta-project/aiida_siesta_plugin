from aiida import orm
from aiida.engine import WorkChain, calcfunction, ToContext, if_
from aiida.common import AttributeDict
from aiida_siesta.workflows.base import SiestaBaseWorkChain
from aiida_siesta.calculations.tkdict import FDFDict
from aiida_siesta.data.eig import EigData
from aiida_siesta.data.pdos import PdosData


@calcfunction
def get_dos(eig, dos_inp):
    """
    From an EigData instance, calls the `compute_dos` method
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
    Drop any molecular dynamics keywords from the parameters dictionary.

    :param param: dict containing the parameters of the calculation
    :return param: dict where the md keywords have been removed
    """
    for item in param.copy().keys():
        trans_item = FDFDict.translate_key(item)
        if trans_item.startswith("md"):
            param.pop(item)
    return param


def add_pdos_block(param, dos_dict):
    """
    Add the PDOS block.
    Also removes a possible previous set up pdos block
    and pdos kpoints. (maybe a warnng needed)
    """
    for item in param.copy().keys():
        trans_item = FDFDict.translate_key(item)
        if trans_item == "%block projecteddensityofstates":
            param.pop(item)
        if trans_item == "pdosgrid?":
            param.pop(item)

    n_ener = int(dos_dict['e_max'] - dos_dict['e_min'] / dos_dict['d_ene'])
    card = f"\n {float(dos_dict['e_min'])} {float(dos_dict['e_max'])} {dos_dict['smearing']} {n_ener} eV"
    card += "\n%endblock projected-density-of-states"
    param["%block projected-density-of-states"] = card

    return param


def set_dos_defaults(dos_inp, eig):
    """
    Method to set the defaults for arguments that have not been specified in input.
    Note that validation is already performed in input.

    :param dos_inp: the dictionary passed as `dos` inputs (or an empty dict)
    :param eig: the EigData instance obtained from the main Siesta calculation.
    :return new_dict: a new dictionary containing the arguments for computing the dos
        where the defaults are set, includes kpoints and `is_metal`. These two are set
        to None if not specified in input, a futher check is implemented in the main code.
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

    smearing = dos_inp.get("smearing", None)
    if smearing is None:
        smearing = 0.05

    also_pdos = dos_inp.get("also_pdos", None)
    if also_pdos is None:
        also_pdos = False

    kpoints_mesh = dos_inp.get("kpoints_mesh", None)

    is_metal = dos_inp.get("is_metal", None)

    new_dict = {
        "e_max": e_max,
        "e_min": e_min,
        "smearing": smearing,
        "d_ene": d_ene,
        "also_pdos": also_pdos,
        "kpoints_mesh": kpoints_mesh,
        "is_metal": is_metal
    }

    return new_dict


def get_val_cond(eig):
    """
    Return the max of the valence and minimum of conduction.
    :param eig: the EigData containing the eigenvalues
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

    return val, cond


def check_eigs(eig, dos_dict, kp_in):  #pylint: disable=too-many-branches
    """
    Analyses the average distance between energies and takes action based on the
    current used kpoints mesh and the selected smearing for the future calculation of the dos.
    In particular, sets a threshold of "minimum required average separation" based on the
    smearing. calculates the separation in the eigenvalues in EigData and, if the
    separation is bigger that the threshold, creates a denser kpoints mesh to reduce the
    separation. This new mesh will be used in a supplementary calculation.

    :param eig: the EigData instance hosting the eigenvalues.
    :param dos_dict: the arguments for `compute_dos of EigData, where, in particular,
        the smearing is stored.
    :param kp_in: the kpoints data (containing the mesh), used for the sista calculation
    :return need_cal: bool signaling whether or not a calculation with improved kp mesh is needed
    :return kpoints: KpointsData vontaining theimproved kpoints mesh.
    """

    from aiida.orm import KpointsData

    if not dos_dict["is_metal"]:
        max_val, min_cond = get_val_cond(eig)

    if dos_dict["is_metal"]:
        thresh = dos_dict["smearing"] / 5
    else:
        thresh = dos_dict["smearing"]

    #print(thresh)

    eigs = eig.get_eigs()

    if dos_dict["is_metal"]:
        count = 0
    else:
        count_cond = 0
        count_val = 0

    for spin in eigs:  #pylint: disable=too-many-nested-blocks
        for kpoint in spin:
            for ene in kpoint:
                if dos_dict["e_max"] > ene > dos_dict["e_min"]:
                    if dos_dict["is_metal"]:
                        count = count + 1
                    else:
                        if ene < eig.e_fermi:
                            count_val = count_val + 1
                        else:
                            count_cond = count_cond + 1

    if dos_dict["is_metal"]:
        mean_separation = (dos_dict["e_max"] - dos_dict["e_min"]) / count
        #print(mean_separation)
    else:
        mean_separation_cond = (dos_dict["e_max"] - min_cond) / count_cond
        mean_separation_val = -(dos_dict["e_min"] - max_val) / count_val
        #print(mean_separation_cond, mean_separation_val)

    need_new_calc = False
    kpoints = None
    if not dos_dict["is_metal"]:
        if mean_separation_cond > thresh or mean_separation_val > thresh:
            mean_separation = max(mean_separation_cond, mean_separation_val)
            need_new_calc = True
            new_mesh = [int(old * (mean_separation / thresh)**(1 / 3)) for old in kp_in.get_kpoints_mesh()[0]]
            kpoints = KpointsData()
            kpoints.set_kpoints_mesh(new_mesh)
    else:
        if mean_separation > thresh:
            need_new_calc = True
            new_mesh = [int(old * (mean_separation / thresh)**(1 / 3)) for old in kp_in.get_kpoints_mesh()[0]]
            kpoints = KpointsData()
            kpoints.set_kpoints_mesh(new_mesh)

    return need_new_calc, kpoints


def validate_dos_mesh(value, _):
    """Validate the `dos.kpoints_mesh` input."""
    if isinstance(value, orm.KpointsData):
        try:
            value.get_kpoints_mesh()
        except AttributeError:
            return '`dos.kpoints_mesh` must contain a kpoints mesh.'
    else:
        if not isinstance(value, str):
            return '`dos.kpoints_mesh` must be a KpointsData or str("automatic")'
        if not value == "automatic":
            return '`dos.kpoints_mesh` must be a KpointsData or str("automatic")'


def validate_d_ene(value, _):
    """Validate the `dos.d_ene` input."""
    if value:
        if not isinstance(value, (float, int)):
            return '`dos.d_ene` must the str or float'
        if value <= 0:
            return '`dos.d_ene` can not be smaller than zero'


def validate_smearing(value, _):
    """Validate the `dos.d_ene` input."""
    #print("ss",value)
    if value:
        if not isinstance(value, (float, int)):
            return '`dos.smearing` must the str or float'
        if value <= 0:
            return '`dos.smearing` can not be smaller than zero'


def validate_inputs(value, _):
    """Validate the entire input namespace."""
    if 'dos' in value:
        if 'e_max' in value['dos'] and 'e_min' in value['dos']:
            if value['dos']['e_max'] <= value['dos']['e_min']:
                return '`dos.e_max` can not be smaller than `dos.e_min`'


class DosPdosWorkChain(WorkChain):
    """
    Workchain to obtain the dos (and pdos if requested) of a structure through SIESTA.
    A separate (usually denser) kpoints mesh (`dos.kpoints_mesh`) can be specified in input for
    the calculation of the dos/pdos. After the end of the self consistent cycle (and possible
    relaxation), the siesta calculation is restarted with the new kpoints mesh (`dos.kpoints_mesh`)
    and the dos/pdos is calculated.
    The `dos.kpoints_mesh` also accepts the string "automatic" that triggers an internal check that
    guesses if the kpoints mesh used for the main calculation is sufficiently dense for the purpose
    of calculating the dos/pdos (see `check_eigs` method). If not, a new kpoints mesh is automatically
    set based on the requested smearing and current used kpoints.
    """

    @classmethod
    def define(cls, spec):
        super().define(spec)
        spec.expose_inputs(SiestaBaseWorkChain, exclude=('metadata',))
        spec.input_namespace('dos', dynamic=True, required=False, help='The inputs for the Dos.')
        spec.input('dos.smearing', validator=validate_smearing, non_db=True, required=False)
        spec.input('dos.e_max', valid_type=(float, int), non_db=True, required=False)
        spec.input('dos.e_min', valid_type=(float, int), non_db=True, required=False)
        spec.input('dos.d_ene', validator=validate_d_ene, non_db=True, required=False)
        spec.input('dos.also_pdos', valid_type=bool, non_db=True, required=False)
        spec.input('dos.is_metal', valid_type=bool, non_db=True, required=False)
        spec.input('dos.kpoints_mesh', validator=validate_dos_mesh, non_db=True, required=False)
        spec.inputs.validator = validate_inputs

        spec.output('dos_array', valid_type=orm.ArrayData)
        spec.output('eig', valid_type=EigData)
        spec.output('pdos', valid_type=PdosData, required=False)

        spec.outline(
            cls.run_main_siesta_wc,
            cls.analyse_dos_inputs,
            if_(cls.need_extra_run)(cls.run_last),
            cls.run_results,
        )
        spec.exit_code(200, 'ERROR_MAIN_WC', message='The main SiestaBaseWorkChain failed')
        spec.exit_code(201, 'ERROR_FINAL_WC', message='The SiestaBaseWorkChain to obtain the bands failed')

    def run_main_siesta_wc(self):
        """
        Run the SiestaBaseWorkChain, might be a relaxation or a scf only.
        """

        inputs = AttributeDict(self.exposed_inputs(SiestaBaseWorkChain))

        running = self.submit(SiestaBaseWorkChain, **inputs)
        self.report(f'Launched SiestaBaseWorkChain<{running.pk}> to perform the siesta calculation.')

        return ToContext(workchain_base=running)

    def analyse_dos_inputs(self):
        """
        Implement decision making based on the dos inputs and the eigenvalues of
        the calculation just run, including the decision to run an extra calculation.
        """
        base_wc = self.ctx.workchain_base
        bwc_outs = base_wc.outputs

        if not base_wc.is_finished_ok:
            return self.exit_codes.ERROR_MAIN_WC

        self.report('Concluded main WorkChain. Analysing Dos inputs.')

        if "dos" in self.inputs:
            inputs_dos = self.inputs.dos
        else:
            inputs_dos = {}

        #Sets defaults for `e_max`, `e_min`, `d_ene`, `smearing`, `also_pdos`
        #if they are not specified in input.
        #Sets `kpoints_mesh` to None if not passed in input. Same for `is_metal`.
        dos_dict = set_dos_defaults(inputs_dos, bwc_outs.eig)

        #Actions if `is_metal` is None.
        if dos_dict["is_metal"] is None:
            if "bands" in base_wc.outputs:
                from aiida.orm.nodes.data.array.bands import find_bandgap
                res = find_bandgap(fermi_energy=bwc_outs.eig.e_fermi, bandsdata=bwc_outs.bands)
                dos_dict["is_metal"] = not res[0]
            else:
                dos_dict["is_metal"] = True

        #Actions for `kpoints_mesh`. Three possibilities:
        #1) kpoints_mesh is None => no need for extra calc, kpoints_mesh stays None.
        #2) kpoints_mesh="automatic" => function called to determine if extra calc needed.
        #   If yes, kpoints_mesh becomes the KpointsData instance containing the mesh.
        #   Otherwise is set to None.
        #3) kpoints_mesh is already a KpointsData instance => only signal extra calc needed.
        new_calc_needed = False
        if dos_dict["kpoints_mesh"] is not None:
            if dos_dict["kpoints_mesh"] == "automatic":
                new_calc_needed, kpoints = check_eigs(bwc_outs.eig, dos_dict, base_wc.inputs.kpoints)
                dos_dict["kpoints_mesh"] = kpoints
            else:
                new_calc_needed = True

        #Even in case nothing was requested in terms of kpoints, an extra calculation must
        #be performed if PDOS is requested
        if dos_dict["also_pdos"]:
            new_calc_needed = True

        self.report(f'Concluded analyses. Need to run again? {new_calc_needed}')

        self.ctx.dos_dict = dos_dict
        self.ctx.need_fin_step = new_calc_needed

    def need_extra_run(self):
        """
        Return whether to run an extra calculation
        """
        return self.ctx.need_fin_step

    def run_last(self):
        """
        In case `need_extra_run` returns True, this function is called to run the
        extra calculation.
        """

        new_calc = self.ctx.workchain_base.get_builder_restart()

        if "output_structure" in self.ctx.workchain_base.outputs:
            new_calc.structure = self.ctx.workchain_base.outputs.output_structure
            new_param = drop_md_keys(new_calc.parameters.get_dict())
            new_calc.parameters = orm.Dict(dict=new_param)

        if self.ctx.dos_dict["also_pdos"]:
            new_param = add_pdos_block(new_calc.parameters.get_dict(), self.ctx.dos_dict)
            new_calc.parameters = orm.Dict(dict=new_param)

        if self.ctx.dos_dict["kpoints_mesh"] is not None:
            new_calc.kpoints = self.ctx.dos_dict["kpoints_mesh"]

        new_calc.parent_calc_folder = self.ctx.workchain_base.outputs.remote_folder
        running = self.submit(new_calc)
        self.report(f'Launched SiestaBaseWorkChain<{running.pk}> with denser kpoints mesh.')
        return ToContext(final_run=running)

    def run_results(self):
        """
        Return the results of the calculation.
        """
        if self.ctx.need_fin_step:
            if not self.ctx.final_run.is_finished_ok:
                return self.exit_codes.ERROR_FINAL_WC
            outps = self.ctx.final_run.outputs
        else:
            outps = self.ctx.workchain_base.outputs

        self.report('Creating outputs')

        #Drop the keys that are not accepted as arguments of `eig.compute_dos`
        self.ctx.dos_dict.pop("is_metal")
        self.ctx.dos_dict.pop("kpoints_mesh")
        self.ctx.dos_dict.pop('also_pdos')

        dos_array = get_dos(outps["eig"], orm.Dict(dict=self.ctx.dos_dict))

        self.out('eig', outps['eig'])
        self.out('dos_array', dos_array)
        if 'pdos' in outps:
            self.out('pdos', outps['pdos'])

    @classmethod
    def inputs_generator(cls):  # pylint: disable=no-self-argument,no-self-use
        from aiida_siesta.utils.protocols_system.input_generators import DosPdosWorkChainInputGenerator
        return DosPdosWorkChainInputGenerator(cls)
