from aiida_siesta.workflows.neb_base import parse_neb
from aiida_siesta.utils.neb import plot_neb


def plot_neb_info(folder, structure):
    """
     Given a 'retrieved' folder from a
     NEB Siesta calculation, and a reference
     structure, compute the associated
     NEB trajectory object and plot the
     barrier curve
    """

    traj = parse_neb(folder, structure)
    plot_neb(traj)
