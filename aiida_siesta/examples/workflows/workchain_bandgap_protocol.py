"""
This example shows how to use the protocol technology inside a
workchain!
Go at the end of the script to understand how to use this workchain!
"""
import aiida.orm
from aiida.engine import submit, run
from aiida.engine import WorkChain, ToContext, calcfunction
from aiida_siesta.workflows.base import SiestaBaseWorkChain
from aiida_siesta.workflows.functions.bandsinputs import SiestaBandsInputsGenerator

@calcfunction
def get_bandgap(e_fermi, band):
    """
    Takes a band object, and a Fermi energy, 
    and extracts the band gap value and 'is_insulator' boolean
    :param band: (orm.BandsData): band-structure object
    :param e_energy: (orm.Float): value of the fermi energy.
    :return: An orm.Dict containing the keys:
             'band_gap':  A float, or None in case of a metal. It is zero when the homo is
                          equal to the lumo (e.g. in semi-metals).
             'band_gap_units': A string, here 'eV'
             'is_insulator': A boolean
    """
    from aiida.orm.nodes.data.array.bands import find_bandgap

    is_insulator, bandgap = find_bandgap(fermi_energy=e_fermi.value,
                                         bandsdata=band)
    output = {}
    output['band_gap'] = bandgap
    output['band_gap_units'] = 'eV'
    output['is_insulator'] = is_insulator
    return orm.Dict(dict=output)


class SiestaBandgapWorkChainProtocol(WorkChain):
    """
    Workchain to obtain the bandgap of a structure through Siesta. It make use of 
    protocols, the inputs are generated through the method get_builder
    of the class  SiestaBandsInputsGenerator. The outputs
    follows some standardization agreed at AiiDA Hackaton of Feb 2020.
    """
    def __init__(self, *args, **kwargs):
        super(SiestaBandgapWorkChainProtocol, self).__init__(*args, **kwargs)

    @classmethod
    def define(cls, spec):
        super(SiestaRelaxWorkChainProtocol, cls).define(spec)
        spec.input('calc_engines', valid_type=orm.Dict)
        spec.input('structure', valid_type=orm.StructureData)
        spec.input('protocol', valid_type=orm.Str, default=Str('standard'))
        spec.input('path_generator', valid_type=orm.Str, help='possible values are Str("seekpath") and Str("legacy")')
        spec.outline(
            cls.setup_and_run,
            cls.run_results,
        )
        # These the standard outouts agreed with other plugins
        spec.output('bands', valid_type=orm.BandsData)
        spec.output('band_gap_info', valid_type=orm.Dict)

    def setup_and_run(self):
        self.report("Setup protocol")

        instgen = SiestaBandsInputsGenerator()
        builder=instgen.get_builder(
                    structure=self.inputs.structure,
                    protocol=self.inputs.protocol.value,
                    path_generator=self.inputs.path_generator,
                    calc_engines=self.inputs.calc_engines,
                    )
       
        ################################################################
        #         HERE THE USER HAS THE FREEDOM TO CHANGE ANY          #
        #        PARAMETER HE/SHE WANTS BEFORE SUMBITTING THE WC       #
        ################################################################

        self.report("Run SiestaWC")
        future = self.submit(builder)
        return ToContext(calc=future)

    def run_results(self):
        self.report("Set outputs")
        self.out('bands', self.ctx.calc.outputs.bands)
        e_fermi = self.ctx.calc.outputs.output_parameters.get_dict()['E_Fermi']
        res_dict = get_bandgap(orm.Float(e_fermi),self.ctx.calc.outputs.bands)
        self.out('band_gap_info', res_dict)


#Here is the code to run the workchain! But be carefull, it can not
#be run from here!!!!!!!!
#Options:
#1) copy the commented code below in a verdi shell and it will run
#2) register this SiestaBandgapWorkChainProtocol as an entry point in
#   setup.json and then you can copy the code below in any file and
#   run it with runaiida, in that case you can replace "run" with "submit"


#from workchain_bandgap_protocol import SiestaBandgapWorkChainProtocol
#from aiida.orm import (Str, Dict, StructureData)
#from aiida.engine import submit, run
#
#calc_engines = {
#     'bands': {
#         'code': 'SiestaHere@localhost',
#         'options': {
#             'resources': {'num_machines': 1, "num_mpiprocs_per_machine": 1},
#             "max_wallclock_seconds": 360, #'queue_name': 'DevQ', 'withmpi': True, 'account': "tcphy113c"
#         }}}
#protocol="stringent"
#path_generator = "seekpath"
#alat = 5.430  # angstrom
#cell = [
#    [
#        0.5 * alat,
#        0.5 * alat,
#        0.,
#    ],
#    [
#        0.,
#        0.5 * alat,
#        0.5 * alat,
#    ],
#    [
#        0.5 * alat,
#        0.,
#        0.5 * alat,
#    ],
#]
#structure = StructureData(cell=cell)
#structure.append_atom(position=(0.000 * alat, 0.000 * alat, 0.000 * alat),
#                      symbols=['Si'])
#structure.append_atom(position=(0.250 * alat, 0.250 * alat, 0.250 * alat),
#                      symbols=['Si'])
#
#inputs={
#        "structure" : structure,
#        "calc_engines" : Dict(dict=calc_engines),
#        "protocol" : Str(protocol),
#        "path_generator" : Str(path_generator),
#    }
#
#run(SiestaBandgapWorkChainProtocol, **inputs)
