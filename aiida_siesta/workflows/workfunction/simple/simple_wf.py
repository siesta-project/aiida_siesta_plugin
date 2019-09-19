#!/usr/bin/env runaiida
# -*- coding: utf-8 -*-

from __future__ import absolute_import
from __future__ import print_function
import os, sys
import numpy as np
from aiida.work.workfunction import workfunction as wf
from aiida.common.example_helpers import test_and_get_code
from aiida.orm.data.base import Float, Str
from aiida.orm.calculation.job.siesta import SiestaCalculation
from aiida.work.process_registry import ProcessRegistry
# from aiida.work.run import run
from aiida.work.run import async
from six.moves import zip

PsfData = DataFactory('psf')
StructureData = DataFactory('structure')
ParameterData = DataFactory('parameter')
KpointsData = DataFactory('array.kpoints')

###############################
# Set Invocation Params
scale_facs = (0.96, 0.98, 1.0, 1.02, 1.04)
labels = ["c1", "c2", "c3", "c4", "c5"]
#codename = 'siesta_icn2work@icn2work'
codename = 'Siesta-4.0@rinaldo'
###############################

@wf
def create_structure():
  """ Si diamond structure """
  alat = 5.430 # angstrom
  cell = np.array([[0.5, 0.5, 0.,],
                   [0., 0.5, 0.5,],
                   [0.5, 0., 0.5,],]) * alat

  # Si
  # This was originally given in the "ScaledCartesian" format
  #
  # StructureData = DataFactory("structure")
  structure = StructureData(cell=cell)
  structure.append_atom(position=(0.000*alat,0.000*alat,0.000*alat),symbols=['Si'])
  structure.append_atom(position=(0.250*alat,0.250*alat,0.250*alat),symbols=['Si'])

  return structure

@wf
def rescale(structure, scale):
  """
  Workfunction to rescale a structure
  :param structure: An AiiDA structure to rescale
  :param scale: The scale factor (for the lattice constant)
  :return: The rescaled structure
  """
  the_ase = structure.get_ase()
  new_ase = the_ase.copy()
  new_ase.set_cell(the_ase.get_cell() * float(scale), scale_atoms=True)
  new_structure = DataFactory("structure")(ase=new_ase)

  return new_structure


@wf
def create_rescaled(scale):
  """
  Workfunction to create and immediately rescale
  a crystal structure of a given element.
  """
  s0 = create_structure()
  return rescale(s0, scale)


def geninputs(structure):
    # The inputs
    inputs = SiestaCalculation.process().get_inputs_template()

    # The structure
    inputs.structure = structure

    # inputs.code = Code.get_from_string(codename)
    inputs.code = Code.get_from_string(codename)
    # calc.label = "PW test"
    # calc.description = "My first AiiDA calculation of Silicon with Quantum ESPRESSO"
    inputs._options.resources = {"num_machines": 1, "num_mpiprocs_per_machine": 1,}
    inputs._options.max_wallclock_seconds = 30 * 60

    # Kpoints
    KpointsData = DataFactory("array.kpoints")
    kpoints = KpointsData()
    kpoints_mesh = 4
    kpoints.set_kpoints_mesh([kpoints_mesh, kpoints_mesh, kpoints_mesh])
    inputs.kpoints = kpoints

    # Calculation parameters
    parameters_dict = {
        'xc-functional': 'LDA',
        'xc-authors': 'CA',
        'max-scfiterations': 50,
        'dm-numberpulay': 4,
        'dm-mixingweight': 0.3,
        'dm-tolerance': 1.e-3,
        'Solution-method': 'diagon',
        'electronic-temperature': '25 meV',
        'md-typeofrun': 'cg',
        'md-numcgsteps': 3,
        'md-maxcgdispl': '0.1 Ang',
        'md-maxforcetol': '0.04 eV/Ang'
    }
    ParameterData = DataFactory("parameter")
    inputs.parameters = ParameterData(dict=parameters_dict)

    # Pseudopotentials
    raw_pseudos = [("Si.psf", 'Si')]
    pseudo_dict = {}
    for fname, kind in raw_pseudos:
      absname = os.path.realpath(os.path.join(os.path.dirname(__file__), fname))
      pseudo, created = PsfData.get_or_create(absname, use_first=True)

      if created:
        print("Created the pseudo for {}".format(kind))
      else:
        print("Using the pseudo for {} from DB: {}".format(kind, pseudo.pk))
      # Attach pseudo node to the calculation
      pseudo_dict[kind] = pseudo

    inputs.pseudo = pseudo_dict

    # Basis set
    inputs.basis = ParameterData(dict={
      'pao-energy-shift': '300 meV',
      '%block pao-basis-sizes': 'Si DZP',
    })

    return inputs


@wf
def run_wf():
  print("Workfunction node identifiers: {}".format(ProcessRegistry().current_calc_node))
  #Instantiate a JobCalc process and create basic structure
  JobCalc = SiestaCalculation.process()
  s0 = create_structure()
  calcs = {}
  for label, factor in zip(labels, scale_facs):
    s = rescale(s0,Float(factor))
    inputs = geninputs(s)
    print("Running a scf for Si with scale factor {}".format(factor))
    # result = run(JobCalc,**inputs)
    result = async(JobCalc,**inputs)


if __name__=="__main__":
  run_wf()
