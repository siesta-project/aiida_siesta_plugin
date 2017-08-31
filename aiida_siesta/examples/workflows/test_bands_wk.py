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
from aiida.orm.data.base import Int, Str
from aiida.orm.data.structure import StructureData
from aiida.work.run import run

from aiida_siesta.workflows.bands import SiestaBandsWorkChain


def parser_setup():
    """
    Setup the parser of command line arguments and return it. This is separated from the main
    execution body to allow tests to effectively mock the setup of the parser and the command line arguments
    """
    parser = argparse.ArgumentParser(
        description='Run the SiestaBandsWorkChain for a given input structure',
    )
    parser.add_argument(
        '-c', type=str, required=True, dest='codename',
        help='the name of the AiiDA code that references Siesta.siesta plugin'
    )
    parser.add_argument(
        '-p', type=str, required=False, dest='protocol', default='standard',
        help='the protocol (default: %(default)s)'
    )
    parser.add_argument(
        '-s', type=int, required=False, dest='structure',
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
        protocol = args.protocol
    except:
        print "Cannot seem to get protocol..."
        protocol = "standard"

    protocol = Str(protocol)
    
    try:
        structure = load_node(args.structure)
    except:
        #
        # Slightly distorted structure
        #
        alat = 5.430 # angstrom
        cell = [[0.5*alat, 0.5*alat, 0.,],
                [0., 0.5*alat, 0.5*alat,],
                [0.5*alat, 0., 0.5*alat,],
        ]

        # Si
        # This was originally given in the "ScaledCartesian" format
        #
        structure = StructureData(cell=cell)
        structure.append_atom(position=(0.000*alat,0.000*alat,0.000*alat),symbols=['Si'])
        structure.append_atom(position=(0.250*alat,0.245*alat,0.250*alat),symbols=['Si'])
        
        #print "Execution failed: failed to load the node for the given structure pk '{}'".format(args.structure)
        #print "Exception report: {}".format(exception)
        #return

    if not isinstance(structure, StructureData):
        print "The provided pk {} for the structure does not correspond to StructureData, aborting...".format(args.parent_calc)
        return

    run(
        SiestaBandsWorkChain,
        code=code,
        structure=structure,
        protocol=protocol
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
