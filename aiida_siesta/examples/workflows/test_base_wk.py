#!/usr/bin/env runaiida
# -*- coding: utf-8 -*-

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
        '-m', type=int, default=20, dest='max_iterations',
        help='the maximum number of iterations to allow in the Workflow. (default: %(default)d)'
    )
    parser.add_argument(
        '-k', nargs=3, type=int, default=[4, 4, 4], dest='kpoints', metavar='Q',
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

    alat = 10.0 # angstrom
    cell = [[alat, 0., 0.,],
            [0., alat, 0.,],
            [0., 0., alat,],
    ]

    # Water molecule
    # One of the H atoms is sligthy moved

    s = StructureData(cell=cell)
    s.append_atom(position=(0.000,0.000,0.00),symbols=['O'])
    s.append_atom(position=(0.757,0.586,0.00),symbols=['H'])
    s.append_atom(position=(-0.780,0.600,0.00),symbols=['H'])


    structure = s


    kpoints = KpointsData()
    kpoints.set_kpoints_mesh(args.kpoints)

    parameters = {
        'meshcutoff': '80.000 Ry',
        'dm:numberpulay': 4,
        'dm:mixingweight': 0.2,
        'dm:tolerance': 1.e-3,
        'max-scfiterations': 30,
        'scf-must-converge': True,
        'geometry-must-converge': True,
        'electronic-temperature': '25 meV',
        'md-typeofrun': 'CG',
        'md-numcgsteps': 6,
        'md-maxcgdispl': '0.1 Ang',
        'md-maxforcetol': '0.03 eV/Ang'
    }
    basis = {
        'pao-energy-shift': '300 meV',
        'pao-basis-size': 'DZP'
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
