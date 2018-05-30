#!/usr/bin/env runaiida
# -*- coding: utf-8 -*-

#
# An example of Workchain to perform geometry relaxation
# Note: The current input structure is non-optimal, in the
# sense that the structure is pulled from the database, while
# the parameters are set here. For example, the parameters are
# taken from the 'test_siesta_geom_fail.py' legacy test, which
# is for a water molecule.
#
import argparse
from aiida.common.exceptions import NotExistent
from aiida.orm.data.base import Int, Str, Float
from aiida.orm.data.parameter import ParameterData
from aiida.orm.data.structure import StructureData
from aiida.orm.data.array.kpoints import KpointsData
from aiida.work.run import run

from aiida_siesta.workflows.vibrawf import SiestaVibraWorkChain


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
        print "Execution failed: could not retrieve the code '{}'".format(args.codename)
        print "Exception report: {}".format(exception)
        return

    try:
        vibra_code = Code.get_from_string(args.vibra_codename)
    except NotExistent as exception:
        print "Execution failed: could not retrieve the code '{}'".format(args.stm_codename)
        print "Exception report: {}".format(exception)
        return

    protocol = Str(args.protocol)
    
    # Structure and supercell. Bulk silicon

    SuperCell_1=1
    SuperCell_2=1
    SuperCell_3=1

    try:
        lxmax=SuperCell_1
    except:
        lxmax=0
    try:
        lymax=SuperCell_2
    except:
        lymax=0
    try:
        lzmax=SuperCell_2
    except:
        lzmax=0

    alat=5.43 # Angstrom. Not passed to the fdf file (only for internal use)
    cell=[[0., alat/2, alat/2,],
          [alat/2, 0., alat/2,],
          [alat/2, alat/2, 0.,]]
    pf=alat*0.125
    na=2
    x0=[[pf,pf,pf],
        [-pf,-pf,-pf]]

    scell=[[0 for x in range(3)] for y in range(3)]
    for i in range(3):
        scell[i][0]=(2*lxmax+1)*cell[i][0]
        scell[i][1]=(2*lymax+1)*cell[i][1]
        scell[i][2]=(2*lzmax+1)*cell[i][2]

    s1=StructureData(cell=cell)
    for i in range(na):
        s1.append_atom(position=(x0[i][0],x0[i][1],x0[i][2]),symbols=['Si'])

    nna=na*(2*lxmax+1)*(2*lymax+1)*(2*lzmax+1)
    xa=[[0 for x in range(3)] for y in range(nna)]
    iatm=-1
    rr=[0 for x in range(3)]
    for i in range(-lxmax,lxmax+1):
        for j in range(-lymax,lymax+1):
            for k in range(-lzmax,lzmax+1):
                rr[0]=i*cell[0][0]+j*cell[0][1]+k*cell[0][2]
                rr[1]=i*cell[1][0]+j*cell[1][1]+k*cell[1][2]
                rr[2]=i*cell[2][0]+j*cell[2][1]+k*cell[2][2]
                for ia in range(na):
                    iatm=iatm+1
                    for ix in range(3):
                        xa[iatm][ix]=x0[ia][ix]+rr[ix]

    s2=StructureData(cell=scell)
    for i in range(nna):
        s2.append_atom(position=(xa[i][0],xa[i][1],xa[i][2]),symbols=['Si'])

    bandskpoints = KpointsData()
    kpp = [(1,1.,1.,1.),
           (15,0.,0.5,0.5),
           (25,0.,0.,0.),
           (20,0.5,0.5,0.5),
           (20,0.,0.5,0.5),
           (15,0.25,0.5,0.75),
           (20,0.5,0.5,0.5)]
    lpp = [[0,'\Gamma'],
           [1,'X'],
           [2,'\Gamma'],
           [3,'L'],
           [4,'X'],
           [5,'W'],
           [6,'L']]
    bandskpoints.set_cell(s1.cell, s1.pbc)
    bandskpoints.set_kpoints(kpp,labels=lpp)

    if args.structure > 0:
        structure = load_node(args.structure)
    else:
        structure = s1
    structure_sc = s2

    run(SiestaVibraWorkChain,
        code=code,
        vibra_code=vibra_code,
        structure=structure,
        structure_sc=structure_sc,
        protocol=protocol,
        bandskpoints=bandskpoints)


def main():
    """
    Setup the parser to retrieve the command line arguments and pass them to the main execution function.
    """
    parser = parser_setup()
    args   = parser.parse_args()
    result = execute(args)


if __name__ == "__main__":
    main()
