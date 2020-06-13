from aiida import orm
from aiida.common.lang import classproperty
from aiida.engine import calcfunction
from aiida_siesta.workflows.base import SiestaBaseWorkChain


@calcfunction
def get_bandgap(e_fermi, band):
    """
    Takes a band object, and a Fermi energy, and extracts the band gap value and 'is_insulator' boolean
    :param band: (orm.BandsData): band-structure object
    :param e_energy: (orm.Float): value of the fermi energy.
    :return: An orm.Dict containing the keys:
        'band_gap': A float, or None in case of a metal. It is zero when the homo is
                          equal to the lumo (e.g. in semi-metals).
        'band_gap_units': A string, here 'eV'
        'is_insulator': A boolean
    """
    from aiida.orm.nodes.data.array.bands import find_bandgap

    is_insulator, bandgap = find_bandgap(fermi_energy=e_fermi.value, bandsdata=band)
    output = {}
    output['band_gap'] = bandgap
    output['band_gap_units'] = 'eV'
    output['is_insulator'] = is_insulator
    return orm.Dict(dict=output)


class BandgapWorkChain(SiestaBaseWorkChain):
    """
    Workchain to obtain the bandgap of a structure through Siesta.
    """

    @classmethod
    def define(cls, spec):
        super(BandgapWorkChain, cls).define(spec)
        spec.output('band_gap_info', valid_type=orm.Dict, required=False)

    def preprocess(self):
        """
        Only pre process is to check that the bandskpoints is defined in input.
        """
        if "bandskpoints" not in self.inputs:
            raise ValueError(
                'you are running the bandgap WorkChain without requesting the bands calculation, '
                'set the port bandskpoints in input'
            )

    def postprocess(self):
        """
        Only post process is the calculation of the band gap from the band, knowing the fermi energy.
        """
        self.report("Obtaining the band gap")
        out_par = self.outputs['output_parameters']
        e_fermi = out_par.get_dict()['E_Fermi']
        res_dict = get_bandgap(orm.Float(e_fermi), self.outputs["bands"])
        self.out('band_gap_info', res_dict)
