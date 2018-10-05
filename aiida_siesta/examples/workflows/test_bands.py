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
from aiida.orm.data.parameter import ParameterData
from aiida.orm.data.structure import StructureData
from aiida.orm.data.array.kpoints import KpointsData
from aiida.work.run import run

from aiida_siesta.workflows.base import SiestaBaseWorkChain


def parser_setup():
    """
    Setup the parser of command line arguments and return it. This is separated from the main
    execution body to allow tests to effectively mock the setup of the parser and the command line arguments
    """
    parser = argparse.ArgumentParser(
        description='Run the SiestaBaseWorkChain for a given input structure',
    )
    parser.add_argument(
        '-m', type=int, default=5, dest='max_iterations',
        help='the maximum number of iterations to allow in the Workflow. (default: %(default)d)'
    )
    parser.add_argument(
        '-k', nargs=3, type=int, default=[1, 1, 1], dest='kpoints', metavar='Q',
        help='define the q-points mesh. (default: %(default)s)'
    )
    parser.add_argument(
        '-c', type=str, required=True, dest='codename',
        help='the name of the AiiDA code that references Siesta.siesta plugin'
    )
    parser.add_argument(
        '-p', type=str, required=True, dest='pseudo_family',
        help='the name of pseudo family to use'
    )
    parser.add_argument(
        '-s', type=int, required=False, dest='structure',
        help='the node id of the structure'
    )
    parser.add_argument(
        '-w', type=int, default=1800, dest='max_wallclock_seconds',
        help='the maximum wallclock time in seconds to set for the calculations. (default: %(default)d)'
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

    kpoints = KpointsData()
    kpoints.set_kpoints_mesh(args.kpoints)
    bandskpoints = KpointsData()

    bandskpoints.set_cell(structure.cell, structure.pbc)
    bandskpoints.set_kpoints_path(kpoint_distance = 0.05)


    parameters = {
                'xc-functional': 'LDA',
                'xc-authors': 'CA',
                'spinpolarized': False,
                'meshcutoff': '150.0 Ry',
                'max-scfiterations': 50,
                'dm-numberpulay': 4,
                'dm-mixingweight': 0.3,
                'dm-tolerance': 1.e-4,
                'Solution-method': 'diagon',
                'electronic-temperature': '25 meV',
                'md-typeofrun': 'cg',
                'md-numcgsteps': 10,
                'md-maxcgdispl': '0.1 Ang',
                'md-maxforcetol': '0.04 eV/Ang'
    }

    # default basis
    basis = {
        'pao-energy-shift': '100 meV',
        '%block pao-basis-sizes': """
        Si DZP                    """,
    }
    
    settings = {}
    options  = {
        'resources': {
            'num_machines': 1
        },
        'max_wallclock_seconds': args.max_wallclock_seconds,
    }

    run(
        SiestaBaseWorkChain,
        code=code,
        structure=structure,
        pseudo_family=Str(args.pseudo_family),
        kpoints=kpoints,
        bandskpoints=bandskpoints,
        parameters=ParameterData(dict=parameters),
        settings=ParameterData(dict=settings),
        options=ParameterData(dict=options),
        basis=ParameterData(dict=basis),
        max_iterations=Int(args.max_iterations),
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
