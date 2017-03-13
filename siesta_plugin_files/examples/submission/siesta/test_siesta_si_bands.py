#!/usr/bin/env runaiida
# -*- coding: utf-8 -*-

__copyright__ = u"Copyright (c), 2015, ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE (Theory and Simulation of Materials (THEOS) and National Centre for Computational Design and Discovery of Novel Materials (NCCR MARVEL)), Switzerland and ROBERT BOSCH LLC, USA. All rights reserved."
__license__ = "MIT license, see LICENSE.txt file"
__version__ = "0.7.0"
__contributors__ = "Andrea Cepellotti, Victor Garcia-Suarez, Alberto Garcia, Emanuele Bosoni"

import sys
import os

from aiida.common.example_helpers import test_and_get_code
from aiida.common.exceptions import NotExistent

################################################################

PsfData = DataFactory('psf')
ParameterData = DataFactory('parameter')
KpointsData = DataFactory('array.kpoints')
StructureData = DataFactory('structure')

try:
    dontsend = sys.argv[1]
    if dontsend == "--dont-send":
        submit_test = True
    elif dontsend == "--send":
        submit_test = False
    else:
        raise IndexError
except IndexError:
    print >> sys.stderr, ("The first parameter can only be either "
                          "--send or --dont-send")
    sys.exit(1)

try:
    codename = sys.argv[2]
except IndexError:
    codename = 'Siesta-4.0@rinaldo'

# If True, load the pseudos from the family specified below
# Otherwise, use static files provided
auto_pseudos = True

queue = None
settings = None

code = test_and_get_code(codename, expected_code_type='siesta')

alat = 5.430 # angstrom
cell = [[0.5*alat, 0.5*alat, 0.,],
        [0., 0.5*alat, 0.5*alat,],
        [0.5*alat, 0., 0.5*alat,],
       ]

# Si
# This was originally given in the "ScaledCartesian" format
#
s = StructureData(cell=cell)
s.append_atom(position=(0.000*alat,0.000*alat,0.000*alat),symbols=['Si'])
s.append_atom(position=(0.250*alat,0.250*alat,0.250*alat),symbols=['Si'])


elements = list(s.get_symbols_set())

if auto_pseudos:
    valid_pseudo_groups = PsfData.get_psf_groups(filter_elements=elements)

    try:
        #pseudo_family = sys.argv[3]
        pseudo_family = 'lda-ag'
    except IndexError:
        print >> sys.stderr, "Error, auto_pseudos set to True. You therefore need to pass as second parameter"
        print >> sys.stderr, "the pseudo family name."
        print >> sys.stderr, "Valid PSF families are:"
        print >> sys.stderr, "\n".join("* {}".format(i.name) for i in valid_pseudo_groups)
        sys.exit(1)

    try:
        PsfData.get_psf_group(pseudo_family)
    except NotExistent:
        print >> sys.stderr, "auto_pseudos is set to True and pseudo_family='{}',".format(pseudo_family)
        print >> sys.stderr, "but no group with such a name found in the DB."
        print >> sys.stderr, "Valid PSF groups are:"
        print >> sys.stderr, ",".join(i.name for i in valid_pseudo_groups)
        sys.exit(1)

parameters = ParameterData(dict={
                'xc:functional': 'LDA',
                'xc:authors': 'CA',
                'spinpolarized': True,
                'meshcutoff': '40.000 Ry',
                'max-scfiterations': 50,
                'dm:numberpulay': 4,
                'dm:mixingweight': 0.3,
                'dm:tolerance': 1.e-3,
                'Solution-method': 'diagon',
                'electronic-temperature': '25 meV',
                'md:typeofrun': 'cg',
                'md:numcgsteps': 3,
                'md:maxcgdispl': '0.1 Ang',
                'md:maxforcetol': '0.04 eV/Ang',
                'writeforces': True,
                'writecoorstep': True,
                'xml:write': True,
                'dm:usesavedm': True
                })

basis = ParameterData(dict={
'pao-energy-shift': '300 meV',
'%block pao-basis-sizes': """
Si DZP                    """,
})

kpoints = KpointsData()

# method mesh
kpoints_mesh = 4
kpoints.set_kpoints_mesh([kpoints_mesh,kpoints_mesh,kpoints_mesh])


# (the object settings is optional)
settings_dict={'test_key': 'test_value'}
settings = ParameterData(dict=settings_dict)

## For remote codes, it is not necessary to manually set the computer,
## since it is set automatically by new_calc
#computer = code.get_remote_computer()
#calc = code.new_calc(computer=computer)

calc = code.new_calc()
calc.label = "Si_bulk"
calc.description = "Siesta test calculation. Si bulk + automatic bands"
calc.set_max_wallclock_seconds(30*60) # 30 min

#------------ clarify this
# Valid only for Slurm and PBS (using default values for the
# number_cpus_per_machine), change for SGE-like schedulers 
## Otherwise, to specify a given # of cpus per machine, uncomment the following:
# calc.set_resources({"num_machines": 1, "num_mpiprocs_per_machine": 8})
#calc.set_resources({"parallel_env": 'openmpi',"tot_num_mpiprocs": 1,"num_machines": 1,"num_cpus": 2})
#------------ clarify this
calc.set_resources({"num_machines": 1})

#calc.set_custom_scheduler_commands("#SBATCH --account=ch3")

if queue is not None:
    calc.set_queue_name(queue)

calc.use_structure(s)
calc.use_parameters(parameters)
calc.use_basis(basis)

if auto_pseudos:
    try:
        calc.use_pseudos_from_family(pseudo_family)
        print "Pseudos successfully loaded from family {}".format(pseudo_family)
    except NotExistent:
        print ("Pseudo or pseudo family not found. You may want to load the "
               "pseudo family, or set auto_pseudos to False.")
        raise
else:
    raw_pseudos = [("Si.psf", 'Si')]

    for fname, kinds, in raw_pseudos:
      absname = os.path.realpath(os.path.join(os.path.dirname(__file__),
                                            "data",fname))
      pseudo, created = PsfData.get_or_create(absname,use_first=True)
      if created:
        print "Created the pseudo for {}".format(kinds)
      else:
        print "Using the pseudo for {} from DB: {}".format(kinds,pseudo.pk)
        
      # Attach pseudo node to the calculation
      calc.use_pseudo(pseudo,kind=kinds)

calc.use_kpoints(kpoints)

# K-points for bands, uncomment your favourite  --------------------
# NOTE: bandskpoints.set_cell(s.cell, s.pbc) HAS TO BE SET ALWAYS ###
bandskpoints = KpointsData()

##..Set a path, label needed, 40 is number of kp between W-L and between L-G..##
#kpp = [('W',  (0.500,  0.250, 0.750), 'L', (0.500,  0.500, 0.500), 40),
#        ('L', (0.500,  0.500, 0.500), 'G', (0., 0., 0.), 40)]
#bandskpoints.set_cell(s.cell, s.pbc)
#bandskpoints.set_kpoints(kpp)


##..........................Only points, no labels............................##
#kpp = [(0.500,  0.250, 0.750), (0.500,  0.500, 0.500), (0., 0., 0.)]
#bandskpoints.set_cell(s.cell, s.pbc)
#bandskpoints.set_kpoints(kpp)

##..kp path automatically generated from structure (all high-simmetry point)..##
##.....labels automatically included, 0.05 is the distance between kpoints....##
# Not available in 0.5.0...
bandskpoints.set_cell(s.cell, s.pbc)
bandskpoints.set_kpoints_path(kpoint_distance = 0.05)

calc.use_bandskpoints(bandskpoints)

if settings is not None:
    calc.use_settings(settings)
#from aiida.orm.data.remote import RemoteData
#calc.set_outdir(remotedata)

if submit_test:
    subfolder, script_filename = calc.submit_test()
    print "Test_submit for calculation (uuid='{}')".format(
        calc.uuid)
    print "Submit file in {}".format(os.path.join(
        os.path.relpath(subfolder.abspath),
        script_filename
        ))
else:
    calc.store_all()
    print "created calculation; calc=Calculation(uuid='{}') # ID={}".format(
        calc.uuid,calc.dbnode.pk)
    calc.submit()
    print "submitted calculation; calc=Calculation(uuid='{}') # ID={}".format(
        calc.uuid,calc.dbnode.pk)

