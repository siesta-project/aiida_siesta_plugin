
from aiida_siesta.workflows.neb_base import parse_neb
from aiida_siesta.utils.neb import plot_neb

def plot_neb_info(f,s):
    """
     Given a 'retrieved' folder f from a
     NEB Siesta calculation, and a reference
     structure 's', compute the associated
     NEB trajectory object and plot the
     barrier curve
    """
     
    t=parse_neb(f,s)
    plot_neb(t)
