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

from aiida_siesta.workflows.vibra import SiestaVibraWorkChain


def parser_setup():
    """
    Setup the parser of command line arguments and return it. This is separated from the main
    execution body to allow tests to effectively mock the setup of the parser and the command line arguments
    """
    parser = argparse.ArgumentParser(
        description='Run the SiestaVibraWorkChain for a given input structure',
    )
    parser.add_argument(
        '-f', type=str, required=True, dest='fcbuild_codename',
        help='the name of the AiiDA code that references siesta.fcbuild plugin'
    )
    parser.add_argument(
        '-c', type=str, required=True, dest='codename',
        help='the name of the AiiDA code that references siesta.siesta plugin'
    )
    parser.add_argument(
        '-v', type=str, required=True, dest='vibrator_codename', 
        help='the name of the AiiDA code that references siesta.vibrator plugin'
    )
    parser.add_argument(
        '-l', nargs=3, type=int, default=[1, 1, 1], dest='scell', metavar='Q',
        help='define the q-points mesh. (default: %(default)s)'
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
        fcbuild_code = Code.get_from_string(args.fcbuild_codename)
    except NotExistent as exception:
        print "Execution failed: could not retrieve the code '{}'".format(args.fcbuild_codename)
        print "Exception report: {}".format(exception)
        return

    try:
        code = Code.get_from_string(args.codename)
    except NotExistent as exception:
        print "Execution failed: could not retrieve the code '{}'".format(args.codename)
        print "Exception report: {}".format(exception)
        return

    try:
        vibrator_code = Code.get_from_string(args.vibrator_codename)
    except NotExistent as exception:
        print "Execution failed: could not retrieve the code '{}'".format(args.stm_codename)
        print "Exception report: {}".format(exception)
        return

    protocol = Str(args.protocol)
    
    # Bulk silicon
    alat = 5.43 # Angstrom. Not passed to the fdf file (only for internal use)
    cell = [[0., alat/2, alat/2,],
            [alat/2, 0., alat/2,],
            [alat/2, alat/2, 0.,],]
    sicd = alat*0.125
    s = StructureData(cell=cell)
    s.append_atom(position=(sicd,sicd,sicd),symbols=['Si'])
    s.append_atom(position=(-sicd,-sicd,-sicd),symbols=['Si'])
    
    fcbparams = {'SuperCell_1': args.scell[0],
                 'SuperCell_2': args.scell[1],
                 'SuperCell_3': args.scell[2]}

    bandskpoints = KpointsData()
    kpp = [(1,2.,2.,2.),
           (15,2.,0.,0.),
           (25,0.,0.,0.),
           (20,1.,1.,1.),
           (20,2.,0.,0.),
           (15,2.,1.,0.),
           (20,1.,1.,1.)]
    lpp = [[0,'\Gamma'],
           [1,'X'],
           [2,'\Gamma'],
           [3,'L'],
           [4,'X'],
           [5,'W'],
           [6,'L']]
    bandskpoints.set_cell(s.cell, s.pbc)
    bandskpoints.set_kpoints(kpp,labels=lpp)

    if args.structure > 0:
        structure = load_node(args.structure)
    else:
        structure = s

    run(SiestaVibraWorkChain,
        fcbuild_code=fcbuild_code,
        code=code,
        vibrator_code=vibrator_code,
        fcbparams=ParameterData(dict=fcbparams),
        structure=structure,
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
