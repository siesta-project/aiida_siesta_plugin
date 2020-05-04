import sys
import os
from aiida.engine import submit
from aiida.orm import load_code
from aiida_siesta.calculations.siesta import SiestaCalculation
from aiida.plugins import DataFactory
from aiida_siesta.workflows.stm import SiestaSTMWorkChain

# This is an example for the submission of the STM WorkChain.
# It shows how to select all the inputs of the WorkChain
# and how to submit it.
# For a quick run, type:
# "runaiida example_stm.py SiestaCode@cmp StmCoode@cmp"
# but, as any other example here, it is made to make the user
# play with all the inputs of the WorkChain


PsfData = DataFactory('siesta.psf')
Dict = DataFactory('dict')
KpointsData = DataFactory('array.kpoints')
StructureData = DataFactory('structure')


try:
    codename = sys.argv[1]
except IndexError:
    codename = 'SiestaHere@localhost'
try:
    stm_codename = sys.argv[2]
except IndexError:
    stm_codename = 'STMhere@localhost'

#
#------------------Code and computer options ---------------------------
#
code = load_code(codename)
stm_code = load_code(stm_codename)

options = {
    "max_wallclock_seconds": 33360,
    'withmpi': True,
    #'account': "tcphy113c",
    #'queue_name': "DevQ",
    "resources": {
        "num_machines": 1,
        "num_mpiprocs_per_machine": 4,
    }
}
#
# Structure -----------------------------------------
#
from pymatgen import Lattice, Structure, Molecule

lattice = Lattice.from_parameters(a=5.77810, b=5.77810, c=5.77810*3.771897, alpha=90, beta=90, gamma=120)

coords = [[0.000078, -0.000047, -0.024346],
          [0.000163, 3.335886, 0.074130],
          [2.889124, 1.667948, -0.049784]]

truct = Structure(lattice, ["Cr", "Cr", "Cr"], coords, coords_are_cartesian=True)

s = StructureData(pymatgen_structure=truct)
#
# Parameters ---------------------------------------------------
#
params_dict = {
    'MeshCutoff': "1200. Ry",
    'spin':'spin-orbit',
    'md-maxcgdispl': '0.1 Ang',
    'md-maxforcetol': '0.04 eV/Ang',
    'write-forces': True,
    'xc-functional': 'GGA',
    'xc-authors': 'PW91',
    'dm-numberpulay': 8,
    'dm-mixingweight': 0.004,
    'dm-tolerance': 1.e-4,
    'electronic-temperature': '300.000 K',
    '%block DMInitSpin':
     """
     1 +3.73  90.0  60.0
     2 +3.73  90.0 -60.0
     3 +3.73  90.0 180.0\n%endblock DMInitSpin"""
}
parameters = Dict(dict=params_dict)
#
# Kpoints
#
kpoints = KpointsData()
kpoints.set_kpoints_mesh([15, 15, 1])
#
# Basis Set Info ------------------------------------------
# The basis dictionary follows the 'parameters' convention
#
basis_dict = {
    'pao-basistype':'split',
    'pao-splitnorm': 0.150,
    'pao-energyshift':'0.020 Ry',
    '%Block PAO-Basis':
    """
    Cr   3
    n=4   0   2
	0.0   0.0
    n=4   1   2
	0.0   0.0
    n=3   2   2
	0.0   0.0\n%endblock PAO.Basis"""
}
basis = Dict(dict=basis_dict)
#--------------------------------------------------------------

#--------------------- Pseudopotentials ---------------------------------
#
# This exemplifies the handling of pseudos for different species
# Those sharing the same pseudo should be indicated.
#
pseudos_dict = {}
raw_pseudos = [("Cr.psf", ['Cr'])]

for fname, kinds, in raw_pseudos:
    absname = os.path.realpath(
        os.path.join(os.path.dirname(__file__), "../plugins/siesta/data/sample-psf-family",
                     fname))
    pseudo, created = PsfData.get_or_create(absname, use_first=True)
    if created:
        print("Created the pseudo for {}".format(kinds))
    else:
        print("Using the pseudo for {} from DB: {}".format(kinds, pseudo.pk))
    for j in kinds:
        pseudos_dict[j]=pseudo


#-----------------------------------------------------------------------
#--All the inputs of a Siesta calculations are listed in a dictionary--
#
inputs = {
    'structure': s,
    'parameters': parameters,
    'code': code,
    'basis': basis,
    'pseudos': pseudos_dict,
    'kpoints': kpoints,
    'options': Dict(dict=options),
    'emin': Float(-1),
    'emax': Float(+0.1),
    'stm_code': stm_code,
    'stm_mode': Str("constant-height"),
    'stm_value': Float(3.0),
    'stm_spin': Str("non-collinear")
    }

process = submit(SiestaSTMWorkChain, **inputs)
print("Submitted workchain; ID={}".format(process.pk))
print(
    "For information about this workchain type: verdi process show {}".format(
        process.pk))
print("For a list of running processes type: verdi process list")

