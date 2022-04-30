# -*- coding: utf-8 -*-
"""
Functions managing xyz files.
"""


def get_positions_from_xyz_file(file):
    """
    Simple parser for xyz file and return positions in a list.
    """
    positions = []
    with open(file, 'r', encoding='utf8') as fileh:
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
    Return a list of StructureData with coordinates taken from .xyz files on a folder.

    :param folder: an absolute path to a folder
    :param ref_struct: a StructureData used as a reference, from it
                       the kinds and cells are taken.
    """
    import glob

    xyz_list = glob.glob(f"{folder}/*.xyz")
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
    From a StructureData, returns an xyz file located in `filename` absolute path.
    """
    xyz_tuple = struct._prepare_xyz()  #pylint: disable=protected-access

    if labels:
        # We add the labels
        with open(filename, 'wb', encoding='utf8') as fileo:
            fileo.write(xyz_tuple[0])
    else:
        # We need to remove the labels
        # First, turn bytes into string and split
        lines = xyz_tuple[0].decode().split("\n")
        with open(filename, "w", encoding='utf8') as fileo:
            #
            # Write first two lines
            #
            fileo.write(f"{lines[0]}\n")
            fileo.write(f"{lines[1]}\n")
            #
            # Write only positions
            #
            for line in lines[2:]:
                pos = line.split()[1:]
                fileo.write(f"{pos[0]} {pos[1]} {pos[2]}\n")
