#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# This script might need to be run with a "framework enabled" python, since it
# uses matplotlib. See: https://matplotlib.org/faq/osx_framework.html
#
# If you are using a virtualenv with AiiDA, your best bet is to define a shell
# function:
#
# function frameworkpython {
#    if [[ ! -z "$VIRTUAL_ENV" ]]; then
#           PYTHONHOME=$VIRTUAL_ENV /usr/local/bin/python "$@"
#    else
#           /usr/local/bin/python "$@"
#    fi
# }
#
# and run the script as:
#
#   frameworkpython get_stm_image.py [ id of stm_array ]
#
#

# Example script to load an stm_array from a Siesta STM calculation and
# generate a contour plot

# Note that we are not using "aiidaenv" in the first line. We need 
# the following two lines before any other AiiDA loads:
#
from aiida import load_dbenv
load_dbenv()

from aiida.orm import load_node

import matplotlib
import numpy as np
import matplotlib.pyplot as plt
import sys

try:
    stm_id = int(sys.argv[1])
except:
    print >> sys.stderr, ("Must provide as parameter the stm_array ID")
    sys.exit(1)
    
arraydata = load_node(stm_id)

X = arraydata.get_array("X")
Y = arraydata.get_array("Y")
Z = arraydata.get_array("Z")

plt.figure()
cp = plt.contour(X, Y, Z)
plt.clabel(cp, inline=True, fontsize=10)
plt.title('Contour Plot')
plt.xlabel('x (bohr)')
plt.ylabel('y (bohr)')
plt.show()

        
        
                
