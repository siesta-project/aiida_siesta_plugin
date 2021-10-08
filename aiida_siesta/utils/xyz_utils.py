def get_positions_from_xyz_file(file):
    """
    Simple parser for xyz file and return positions in a list
    """

    positions = []
    with open(file, 'r') as fileh:
        lines = fileh.readlines()
        for line in lines[2:]:
            parts = line.split()
            # Support the case in which the species label is present
            if len(parts) == 4:
                start = 1
            else:
                start = 0
            pos = [float(i) for i in parts[start:]]
            positions.append(pos)

    return positions


def get_structure_list_from_folder(folder, ref_struct):
    """
    Return a list of StructureData with coordinates taken from .xyz
    files on a folder.
    :param folder: an absolute path to a folder
    :param ref_struct: a StructureData used as a reference, from it
                       the kinds and cells are taken.
    """

    import glob

    xyz_list = glob.glob("{}/*.xyz".format(folder))
    xyz_list.sort()

    # Compute number of expected (physical) sites and use it
    # below to discard ghost sites (which are always trailing the rest)
    nsites = len(ref_struct.sites)

    structure_list = []
    for file in xyz_list:
        positions = get_positions_from_xyz_file(file)
        struct = ref_struct.clone()
        struct.reset_sites_positions(positions[0:nsites])
        structure_list.append(struct)

    return structure_list


def write_xyz_file_from_structure(struct, filename, labels=True):
    """
    From a StructureData, returns an xyz file located in `filename`
    absolute path.
    """

    xyz_tuple = struct._prepare_xyz()

    if labels:
        # We add the labels
        open(filename, 'wb').write(xyz_tuple[0])
    else:
        # We need to remove the labels
        # First, turn bytes into string and split
        lines = xyz_tuple[0].decode().split("\n")
        with open(filename, "w") as fileo:
            #
            # Write first two lines
            #
            fileo.write("{}\n".format(lines[0]))
            fileo.write("{}\n".format(lines[1]))
            #
            # Write only positions
            #
            for line in lines[2:]:
                pos = line.split()[1:]
                fileo.write("{} {} {}\n".format(pos[0], pos[1], pos[2]))


#def write_xyz_files_from_trajectory(t, ref_struct, pattern, labels=True):
#    """
#   Given a trajectory and a file pattern, write
#   xyz files from trajectory
#
#   :param: t: TrajectoryData object
#   :param: ref_struct: Needed for kinds list
#   :pattern: string to form the filenames
#   :labels: whether to print the atomic labels or not
#   """
#    kinds = ref_struct.kinds
#
#    # loop over structures
#    for i in range(t.numsteps):
#        s = t.get_step_structure(i, custom_kinds=kinds)
#        # write a xyz file with a standard prefix in the folder
#        # Note that currently we do not want the labels in these files
#        filename = pattern + "{}.xyz".format(i)
#        write_xyz_file_from_structure(s, filename, labels=labels)
