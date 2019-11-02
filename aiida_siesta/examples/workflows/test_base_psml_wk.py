#!/usr/bin/env runaiida
# -*- coding: utf-8 -*-

#
# Example to test the new dual psf-psml support.
# You can pass a "psf" or a "psml" family (or a mixed one, soon)
# To load the 'sample-psml-family' pseudos in ../plugins/siesta/data, use
# the 'verdi data psml uploadfamily' command.
#
# (The original test used a psml file with only relativistic projectors, hence
# the need for the 'spin:SO' setting) (see Si.psml in examples/plugin/siesta/data)

from __future__ import absolute_import
from __future__ import print_function
import six

import argparse
from aiida.common.exceptions import NotExistent
from aiida.orm import Int, Str
from aiida.plugins import DataFactory
from aiida.engine import run

from aiida_siesta.workflows.base import SiestaBaseWorkChain
from six.moves import range

Dict = DataFactory('dict')
KpointsData = DataFactory('array.kpoints')
StructureData = DataFactory('structure')


def parser_setup():
    """
    Setup the parser of command line arguments and return it. This is separated from the main
    execution body to allow tests to effectively mock the setup of the parser and the command line arguments
    """
    parser = argparse.ArgumentParser(
        description='Run the SiestaBaseWorkChain for a given input structure',
    )
    parser.add_argument(
        '-m',
        type=int,
        default=8,
        dest='max_iterations',
        help=
        'the maximum number of iterations to allow in the Workflow. (default: %(default)d)'
    )
    parser.add_argument(
        '-k',
        nargs=3,
        type=int,
        default=[4, 4, 4],
        dest='kpoints',
        metavar='Q',
        help='define the q-points mesh. (default: %(default)s)')
    parser.add_argument(
        '-c',
        type=str,
        required=True,
        dest='codename',
        help='the name of the AiiDA code that references Siesta.siesta plugin')
    parser.add_argument(
        '-p',
        type=str,
        required=True,
        dest='pseudo_family',
        help='the name of pseudo family to use')
    parser.add_argument(
        '-w',
        type=int,
        default=1800,
        dest='max_wallclock_seconds',
        help=
        'the maximum wallclock time in seconds to set for the calculations. (default: %(default)d)'
    )

    return parser


def execute(args):
    """
    The main execution of the script, which will run some preliminary checks on the command
    line arguments before passing them to the workchain and running it
    """
    # try:
    #     code = load_code(args.codename)
    # except NotExistent as exception:
    #     print("Execution failed: could not retrieve the code {}".format(args.codename))
    #     print("Exception report: {}".format(exception))
    #     return

    code = load_code(args.codename)

    # Structure. Bulk silicon

    alat = 5.43  # Angstrom. Not passed to the fdf file (only for internal use)
    cell = [[
        0.,
        alat / 2,
        alat / 2,
    ], [
        alat / 2,
        0.,
        alat / 2,
    ], [
        alat / 2,
        alat / 2,
        0.,
    ]]
    pf = alat * 0.125
    na = 2
    x0 = [[pf, pf, pf], [-pf, -pf, -pf]]

    structure = StructureData(cell=cell)
    for i in range(na):
        structure.append_atom(
            position=(x0[i][0], x0[i][1], x0[i][2]), symbols=['Si'])

    kpoints = KpointsData()
    kpoints.set_kpoints_mesh(args.kpoints)

    parameters = {
        'meshcutoff': '80.000 Ry',
        'dm-numberpulay': 4,
        'dm-mixingweight': 0.2,
        'dm-tolerance': 1.e-3,
        'spin': 'SO',
        'max-scfiterations': 30,
        'scf-must-converge': True,
        'geometry-must-converge': True,
        'electronic-temperature': '25 meV',
        'md-typeofrun': 'CG',
        'md-numcgsteps': 6,
        'md-maxcgdispl': '0.1 Ang',
        'md-maxforcetol': '0.03 eV/Ang'
    }
    basis = {'pao-energy-shift': '300 meV', 'pao-basis-size': 'DZP'}
    settings = {}
    options = {
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
        parameters=Dict(dict=parameters),
        settings=Dict(dict=settings),
        options=Dict(dict=options),
        basis=Dict(dict=basis),
        max_iterations=Int(args.max_iterations),
    )


def main():
    """
    Setup the parser to retrieve the command line arguments and pass them to the main execution function.
    """
    parser = parser_setup()
    args = parser.parse_args()
    result = execute(args)


if __name__ == "__main__":
    main()
