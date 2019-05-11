#!/usr/bin/env runaiida
# -*- coding: utf-8 -*-
from __future__ import absolute_import
from __future__ import print_function

import os.path as op

from aiida.engine import run

#Loading calculation we want to restart (old_calc)
g=load_node(375)

#Set up the a new calculation with all
#the inputs of the old one (we use the builder)
restart=g.get_builder_restart()

#The inputs of old_calc attched to restart are
#already stored!!! Can't modify them straight away
#If you want to change something you make a clone
#and reassign it to the builder.
#Here we change dm-tollerance for example
newpar=restart.parameters.clone()
newpar.attributes["dm-tolerance"]=0.00001
restart.parameters=newpar


#We need to take care here of passing the
#output geometry of old_calc to the new calculation
if g.outputs.output_parameters.attributes["variable_geometry"]:
   restart.structure=g.outputs.output_structure

#The most important line. The presence of
#parent_calc_folder triggers the real restart
#meaning the copy of the .DM and the
#addition of use-saved-dm to the parameters
restart.parent_calc_folder=g.outputs.remote_folder


restart.metadata.dry_run=True
restart.metadata.store_provenance=False
run(restart)

