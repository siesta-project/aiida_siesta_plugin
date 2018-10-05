#!/usr/bin/env runaiida
# -*- coding: utf-8 -*-

__copyright__ = u"Copyright (c), 2015, ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE (Theory and Simulation of Materials (THEOS) and National Centre for Computational Design and Discovery of Novel Materials (NCCR MARVEL)), Switzerland and ROBERT BOSCH LLC, USA. All rights reserved."
__license__ = "MIT license, see LICENSE.txt file"
__version__ = "0.7.0"
__contributors__ = "Andrea Cepellotti, Victor Garcia-Suarez, Alberto Garcia, Emanuele Bosoni"

import sys
import os
from aiida.tools import get_explicit_kpoints_path
from aiida.common.example_helpers import test_and_get_code
from aiida.common.exceptions import NotExistent

#In this example we will calculate the band structure of Si.
#Thanks to SeeK-path we can automatically generate the 
#high symmetry points path where to calculate the bands.
#Alternatively, a manual path or list of k-points can be set.

################################################################

PsfData = DataFactory('siesta.psf')
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
    codename = 'siesta4.0.1@parsons'

#
##----------Set calculation----------------------
##For remote codes, it is not necessary to manually set the computer,
##since it is set automatically by new_calc
#computer = code.get_remote_computer()
#calc = code.new_calc(computer=computer)
code = test_and_get_code(codename, expected_code_type='siesta.siesta')
calc = code.new_calc()
calc.label = "Si_bulk"
calc.description = "Siesta test calculation. Si bulk + automatic bands"
calc.set_max_wallclock_seconds(30*60) # 30 min
calc.set_resources({"num_machines": 1})

queue = None
if queue is not None:
    calc.set_queue_name(queue)

#
##--------------Settings------------------
##The object settings is optional.
##settings_dict={'test_key': 'test_value'}
##settings = ParameterData(dict=settings_dict)
settings = None
if settings is not None:
    calc.use_settings(settings)

#
##-------------------Structure-----------------------------------
##Manually set the structure, all the quantities must be in Ang.
##Then, we pass through SeeK-path, to get the standardized cell,
##necessary for the automatic choice of the bands path.

alat = 5.430 # angstrom
cell = [[0.5*alat, 0.5*alat, 0.,],
        [0., 0.5*alat, 0.5*alat,],
        [0.5*alat, 0., 0.5*alat,],
       ]

s = StructureData(cell=cell)
s.append_atom(position=(0.000*alat,0.000*alat,0.000*alat),symbols=['Si'])
s.append_atom(position=(0.250*alat,0.250*alat,0.250*alat),symbols=['Si'])
elements = list(s.get_symbols_set())
seekpath_parameters = ParameterData(dict={'reference_distance': 0.02,'symprec': 0.0001})
result=get_explicit_kpoints_path(s, **seekpath_parameters.get_dict())
newstructure = result['primitive_structure']
calc.use_structure(newstructure)


#
##---------------------Pseudos---------------------------
##If auto_pseudos = True, load the pseudos from the family specified
##below. Otherwise, use static files provided
auto_pseudos = False
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


#
##----------Calculation parameters----------------
parameters = ParameterData(dict={
                'xc-functional': 'LDA',
                'xc-authors': 'CA',
                'spin-polarized': True,
                'meshcutoff': '40.000 Ry',
                'max-scfiterations': 50,
                'dm-numberpulay': 4,
                'dm-mixingweight': 0.3,
                'dm-tolerance': 1.e-3,
                'Solution-method': 'diagon',
                'electronic-temperature': '25 meV',
                'md-typeofrun': 'cg',
                'md-numcgsteps': 3,
                'md-maxcgdispl': '0.1 Ang',
                'md-maxforcetol': '0.04 eV/Ang',
                'writeforces': True,
                'writecoorstep': True,
                'dm-usesavedm': True
                })
calc.use_parameters(parameters)


#
##-----------------Basis-set-------------------
basis = ParameterData(dict={
'pao-energy-shift': '300 meV',
'%block pao-basis-sizes': """
Si DZP                    """,
})
calc.use_basis(basis)


#
##------------KPoints for calculation---------
kpoints = KpointsData()
# method mesh
kpoints_mesh = 4
kpoints.set_kpoints_mesh([kpoints_mesh,kpoints_mesh,kpoints_mesh])
calc.use_kpoints(kpoints)

#
##-------------------K-points for bands --------------------
bandskpoints = KpointsData()
##Uncomment your favourite, three options:

##1)
##.....Making use of SeeK-path for the automatic path......
##The choice of the distance between kpoints is in the call seekpath_parameters
##All high symmetry points included, labels already included
#bandskpoints=result['explicit_kpoints']

##2)
##.....Only points, no labels.......
##Mandatory to set cell and pbc
#kpp = [(0.500,  0.250, 0.750), (0.500,  0.500, 0.500), (0., 0., 0.)]
#bandskpoints.set_cell(newstructure.cell, newstructure.pbc)
#bandskpoints.set_kpoints(kpp)

##3)
##.....Set a manual path.......
##Labels needed, 40 is number of kp between W-L and between L-G...
##With new version of aiida, the use of this functionality is a bit involved.
##First we have to call a function in tools.data.array.kpoints.legacy.
##The funcion is called get_explicit_kpoints_path, but is not the same
##function I used before.
##Second step is to assign the results of the funtion to a KpointsData()
##Mandatory to set cell and pbc
kpp = [('W',  (0.500,  0.250, 0.750), 'L', (0.500,  0.500, 0.500), 40),
        ('L', (0.500,  0.500, 0.500), 'G', (0., 0., 0.), 40)]
from aiida.tools.data.array.kpoints import legacy
fs=legacy.get_explicit_kpoints_path(kpp)
bandskpoints=KpointsData()
bandskpoints.set_cell(newstructure.cell, newstructure.pbc)
bandskpoints.set_kpoints(fs[3])
bandskpoints.labels=fs[4]

calc.use_bandskpoints(bandskpoints)

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

