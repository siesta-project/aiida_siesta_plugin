###### from aiida.engine import calcfunction
######  @calcfunction
#def parse_neb(retrieved, ref_structure):
#    """
#    Wrapper to preserve provenance.
#    :param: retrieved:  the retrieved folder from a NEB calculation
#                        (containing .xyz files and NEB data files)
#    :param: ref_structure: a reference structure
#    :return: a Trajectory object generated from the .xyz files, and
#             with extra arrays for NEB results.
#    """
#    import os
#    from aiida.orm import TrajectoryData
#    from aiida_siesta.utils.xyz_utils import get_structure_list_from_folder
#    from aiida_siesta.utils.neb import parse_neb_results
#
#    folder_path = retrieved._repository._get_base_folder().abspath
#    struct_list = get_structure_list_from_folder(folder_path, ref_structure)
#
#    traj = TrajectoryData(struct_list)
#
#    neb_results_file = 'NEB.results'
#    if neb_results_file in retrieved._repository.list_object_names():
#        neb_results_path = os.path.join(folder_path, neb_results_file)
#        annotated_traj = parse_neb_results(neb_results_path, traj)
#
#        _kinds_raw = [k.get_raw() for k in ref_structure.kinds]
#        annotated_traj.set_attribute('kinds', _kinds_raw)
#
#    return annotated_traj


def parse_neb_results(file, traj_in):
    """
    Parses NEB.results
    :param: file: NEB results
    :param: traj_in: TrajectoryData object with final MEP images

    :return: Extended trajectory object with NEB data arrays
             and estimation of barrier, and number of iterations.
    """
    import numpy as np

    n_images = traj_in.numsteps

    # digest the whole file

    data = np.loadtxt(file)

    number_of_neb_iterations = int(len(data) / n_images)

    # Get the data for the final iteration
    final = data[-n_images:]

    # Create a new object for hygiene
    traj = traj_in.clone()

    energies = final[:, 2]
    min_neb = max(energies)
    max_neb = min(energies)
    barrier = abs(max_neb - min_neb)

    traj.set_attribute('barrier', barrier)
    traj.set_attribute('neb_iterations', number_of_neb_iterations)
    traj.set_array('reaction_coordinates', final[:, 1])
    traj.set_array('energies', energies)
    traj.set_array('ediff', final[:, 3])
    traj.set_array('curvature', final[:, 4])
    traj.set_array('max_force', final[:, 5])

    return traj


def plot_neb(traj):
    """
    Plot the neb energies and the value of the computed barrier
    in a file NEB.png
    :param traj: a TrajectoryData poduced by neb_base workchain
    """

    import matplotlib.pyplot as plt
    import numpy as np
    from scipy.interpolate import interp1d

    #im = traj.get_array('steps')
    x = traj.get_array('reaction_coordinates')
    y = traj.get_array('ediff')
    #y2 = traj.get_array('energies')

    barrier = round(traj.get_attribute('barrier'), 3)

    xnew = np.linspace(0, x[len(x) - 1], num=1000, endpoint=True)
    fun1 = interp1d(x, y, kind='linear')
    fun2 = interp1d(x, y, kind='cubic')
    fun3 = interp1d(x, y, kind='quadratic')
    plt.plot(x, y, "o", xnew, fun1(xnew), "-", xnew, fun2(xnew), "--", xnew, fun3(xnew), 'r')
    plt.title("Barrier Energy = " + str(barrier) + " eV")
    plt.legend(['data', 'linear', 'cubic', 'quadratic'], loc='best')

    plt.savefig("NEB.png")
    plt.show()
