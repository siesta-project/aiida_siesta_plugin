#!/usr/bin/env runaiida
# -*- coding: utf-8 -*-

#
# Siesta+Vibra workflow.
# It is advised to use this example as a template submission script.
#
from __future__ import absolute_import
from __future__ import print_function
import argparse
from aiida.common.exceptions import NotExistent
from aiida.orm import Int, Str, Float, Dict, StructureData, KpointsData, ArrayData
from aiida.engine.launch import run
from aiida_siesta.workflows.vibrawf import SiestaVibraWorkChain
import numpy as np
from six.moves import range


def parser_setup():
    """
    Setup the parser of command line arguments and return it. This is separated from the main
    execution body to allow tests to effectively mock the setup of the parser and the command line arguments
    """
    parser = argparse.ArgumentParser(
        description='Run the SiestaVibraWorkChain for a given input structure',
    )
    parser.add_argument(
        '-c', type=str, required=True, dest='codename',
        help='the name of the AiiDA code that references siesta.siesta plugin'
    )
    parser.add_argument(
        '-v', type=str, required=True, dest='vibra_codename',
        help='the name of the AiiDA code that references siesta.vibra plugin'
    )
    parser.add_argument(
        '-p', type=str, required=False, dest='protocol', default='standard',
        help='the protocol (default: %(default)s)'
    )
    parser.add_argument(
        '-s', type=int, required=False, dest='structure', default=0,
        help='the node id of the structure'
    )

    return parser


def execute(args):
    """
    The main execution of the script, which will run some preliminary checks on the command
    line arguments before passing them to the workchain and running it
    """
    try:
        code = Code.get_from_string(args.codename)
    except NotExistent as exception:
        print("Execution failed: could not retrieve the code '{}'".format(args.codename))
        print("Exception report: {}".format(exception))
        return

    try:
        vibra_code = Code.get_from_string(args.vibra_codename)
    except NotExistent as exception:
        print("Execution failed: could not retrieve the code '{}'".format(args.stm_codename))
        print("Exception report: {}".format(exception))
        return

    protocol = Str(args.protocol)

    # Structure. Bulk silicon

    # Supercell Size
    SuperCell_1=1
    SuperCell_2=1
    SuperCell_3=1

    scnumbers=np.array([SuperCell_1,SuperCell_2,SuperCell_3])
    scarray=ArrayData()
    scarray.set_array('sca',scnumbers)

    alat=5.43 # Angstrom. Not passed to the fdf file (only for internal use)
    cell=[[0., alat/2, alat/2,],
          [alat/2, 0., alat/2,],
          [alat/2, alat/2, 0.,]]
    pf=alat*0.125
    na=2
    x0=[[pf,pf,pf],
        [-pf,-pf,-pf]]

    s1=StructureData(cell=cell)
    for i in range(na):
        s1.append_atom(position=(x0[i][0],x0[i][1],x0[i][2]),symbols=['Si'])

    bandskpoints = KpointsData()
    # kpp = [(1,1.,1.,1.),
    #        (15,0.,0.5,0.5),
    #        (25,0.,0.,0.),
    #        (20,0.5,0.5,0.5),
    #        (20,0.,0.5,0.5),
    #        (15,0.25,0.5,0.75),
    #        (20,0.5,0.5,0.5)]
    # lpp = [[0,'\Gamma'],
    #        [1,'X'],
    #        [2,'\Gamma'],
    #        [3,'L'],
    #        [4,'X'],
    #        [5,'W'],
    #        [6,'L']]
    bandskpoints.set_cell(s1.cell, s1.pbc)
    # bandskpoints.set_kpoints(kpp,labels=lpp)
    bandskpoints.set_kpoints_path([
        ('G', 'X', 20),
        ('X', 'K', 10),
        ('K', 'G', 20),
        ('G', 'L', 20),
        ('L', 'X', 20),
        ('X', 'W', 5),
        ('W', 'L', 10),
    ])

    if args.structure > 0:
        structure = load_node(args.structure)
    else:
        structure = s1

    global_parameters =  Dict(dict={
        'atomicdispl': '0.0211672 Ang'
        # 'atomicdispl': '0.02  Ang'
    })

    siesta_parameters =  Dict(dict={
        # 'dm_convergence_threshold': 1.0e-5
    })

    vibra_parameters = Dict(dict={
        'eigenvectors': False
    })

    run(SiestaVibraWorkChain,
        code=code,
        vibra_code=vibra_code,
        scarray=scarray,
        structure=structure,
        protocol=protocol,
        bandskpoints=bandskpoints,
        global_parameters=global_parameters,
        siesta_parameters=siesta_parameters,
        vibra_parameters=vibra_parameters,
    )


def main():
    """
    Setup the parser to retrieve the command line arguments and pass them to the main execution function.
    """
    parser = parser_setup()
    args   = parser.parse_args()
    result = execute(args)


if __name__ == "__main__":
    main()
