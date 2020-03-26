#!/usr/bin/env runaiida

# Example script to load an stm_array from a Siesta STM calculation and
# generate a contour plot.
# To run runaiida get_stm_image.py [ id of stm_array ]

# This script might need to be run with a "framework enabled" python, since it
# uses matplotlib. See: https://matplotlib.org/faq/osx_framework.html
# (Note: it seems to work with python 3)
#
# If you are using a virtualenv with AiiDA, you could define a shell
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



from aiida.orm import load_node
import matplotlib
import numpy as np
import matplotlib.pyplot as plt
import sys

try:
    stm_id = int(sys.argv[1])
except:
    print("Must provide as parameter the stm_array ID", sys.stderr)
    sys.exit(1)

arraydata = load_node(stm_id)

noncoll=False
X = arraydata.get_array("grid_X")
Y = arraydata.get_array("grid_Y")
try:
    Z = arraydata.get_array("STM")
except:
    Z = arraydata.get_array("STM_q")
    Sx = arraydata.get_array("STM_sx")
    Sy = arraydata.get_array("STM_sy")
    Sz = arraydata.get_array("STM_sz")
    noncoll=True


if noncoll:
    fig, ax = plt.subplots(2,2)
    cp = ax[0, 0].contour(X, Y, Z)
    ax[0, 0].clabel(cp, inline=True, fontsize=10)
    ax[0, 0].set_aspect('equal', 'box')
    ax[0, 0].set_title('Total charge (q)')
    ax[0, 0].set_xlabel('x (bohr)')
    ax[0, 0].set_ylabel('y (bohr)')
    cp = ax[0, 1].contour(X, Y, Sx)
    ax[0, 1].clabel(cp, inline=True, fontsize=10)
    ax[0, 1].set_aspect('equal', 'box')
    ax[0, 1].set_title('Spin direction x')
    ax[0, 1].set_xlabel('x (bohr)')
    ax[0, 1].set_ylabel('y (bohr)')
    cp = ax[1, 0].contour(X, Y, Sy)
    ax[1, 0].clabel(cp, inline=True, fontsize=10)
    ax[1, 0].set_aspect('equal', 'box')
    ax[1, 0].set_title('Spin direction y')
    ax[1, 0].set_xlabel('x (bohr)')
    ax[1, 0].set_ylabel('y (bohr)')
    cp = ax[1, 1].contour(X, Y, Sz)
    ax[1, 1].clabel(cp, inline=True, fontsize=10)
    ax[1, 1].set_aspect('equal', 'box')
    ax[1, 1].set_title('Spin direction z')
    ax[1, 1].set_xlabel('x (bohr)')
    ax[1, 1].set_ylabel('y (bohr)')
    fig.tight_layout()
    plt.show()
else:
    fig, ax = plt.subplots(1,1)
    cp = ax.contour(X, Y, Z)
    ax.clabel(cp, inline=True, fontsize=10)
    ax.set_aspect('equal', 'box')
    ax.set_title('Total charge (q)')
    ax.set_xlabel('x (bohr)')
    ax.set_ylabel('y (bohr)')
    fig.tight_layout()
    #plt.figure()
    #cp = plt.contour(X, Y, Z)
    #plt.clabel(cp, inline=True, fontsize=10)
    #plt.title('Contour Plot')
    #plt.xlabel('x (bohr)')
    #plt.ylabel('y (bohr)')
    plt.show()
