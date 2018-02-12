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
from aiida.orm.data.structure import StructureData
from aiida.work.run import run

from aiida_siesta.workflows.stm import SiestaSTMWorkChain


def parser_setup():
    """
    Setup the parser of command line arguments and return it. This is separated from the main
    execution body to allow tests to effectively mock the setup of the parser and the command line arguments
    """
    parser = argparse.ArgumentParser(
        description='Run the SiestaSTMWorkChain for a given input structure',
    )
    parser.add_argument(
        '-c', type=str, required=True, dest='codename',
        help='the name of the AiiDA code that references Siesta.siesta plugin'
    )
    parser.add_argument(
        '-t', type=str, required=True, dest='stm_codename', 
        help='the name of the AiiDA code that references Siesta.stm plugin'
    )
    parser.add_argument(
        '-p', type=str, required=False, dest='protocol', default='standard',
        help='the protocol (default: %(default)s)'
    )
    parser.add_argument(
        '-s', type=int, required=False, dest='structure', default=0,
        help='the node id of the structure'
    )
    parser.add_argument(
        '-z', type=float, required=False, dest='height', default=7.5,
        help='the height (in Ang) at which to compute the image'
    )
    parser.add_argument(
        '-e', type=float, required=False, dest='e1', default=-5.0,
        help='the lower limit of the energy window'
    )
    parser.add_argument(
        '-E', type=float, required=False, dest='e2', default=1.0,
        help='the upper limit of the energy window'
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
        stm_code = Code.get_from_string(args.stm_codename)
    except NotExistent as exception:
        print "Execution failed: could not retrieve the code '{}'".format(args.stm_codename)
        print "Exception report: {}".format(exception)
        return


    height = Float(args.height)
    e1 = Float(args.e1)
    e2 = Float(args.e2)

    protocol = Str(args.protocol)
    

            
    alat = 15. # angstrom
    cell = [[alat, 0., 0.,],
            [0., alat, 0.,],
            [0., 0., alat,],
    ]

    # Benzene molecule
    #
    s = StructureData(cell=cell)

    def perm(x,y,z):
        return (z,y+0.5*alat,0.5*alat)
    
    s.append_atom(position=perm(0.000,0.000,0.468),symbols=['H'])
    s.append_atom(position=perm(0.000,0.000,1.620),symbols=['C'])
    s.append_atom(position=perm(0.000,-2.233,1.754),symbols=['H'])
    s.append_atom(position=perm(0.000,2.233,1.754),symbols=['H'])
    s.append_atom(position=perm(0.000,-1.225,2.327),symbols=['C'])
    s.append_atom(position=perm(0.000,1.225,2.327),symbols=['C'])
    s.append_atom(position=perm(0.000,-1.225,3.737),symbols=['C'])
    s.append_atom(position=perm(0.000,1.225,3.737),symbols=['C'])
    s.append_atom(position=perm(0.000,-2.233,4.311),symbols=['H'])
    s.append_atom(position=perm(0.000,2.233,4.311),symbols=['H'])
    s.append_atom(position=perm(0.000,0.000,4.442),symbols=['C'])
    s.append_atom(position=perm(0.000,0.000,5.604),symbols=['H'])

    if args.structure > 0:
        structure = load_node(args.structure)
    else:
        structure = s

    run(
        SiestaSTMWorkChain,
        code=code,
        stm_code=stm_code,
        structure=structure,
        protocol=protocol,
        height=height,
        e1=e1,
        e2=e2
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
