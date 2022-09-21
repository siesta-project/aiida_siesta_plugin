# -*- coding: utf-8 -*-
"""
Collect function to process structures.
"""


def clone_aiida_structure(structure):
    """
    A cloned structure is not quite ready to store more atoms.

    Need to check if it is a bug and origin!
    """
    truct = structure.clone()
    truct._internal_kind_tags = {}  # pylint: disable=protected-access

    return truct


def add_ghost_sites_to_structure(original_structure, basis):
    """
    Add the ghost states to the structure.

    Returns a structure object, and a list of ghost species names,
    appending the relevant sites to the original structure, as
    directed by the floating_sites element of the basis dictionary.
    It currently accepts a full basis dictionary, but it would be
    more appropriate if it took only a list of ghost dictionaries.
    """
    structure = clone_aiida_structure(original_structure)
    floating_species_names = []
    #Add ghosts to the structure
    if basis is not None:
        basis_dict = basis.get_dict()
        floating = basis_dict.pop('floating_sites', None)
        if floating is not None:
            original_kind_names = [kind.name for kind in original_structure.kinds]
            for item in floating:
                if item["name"] in original_kind_names:
                    raise ValueError(
                        "It is not possibe to specify `floating_sites` "
                        "(ghosts states) with the same name of a structure kind."
                    )
                structure.append_atom(position=item["position"], symbols=item["symbols"], name=item["name"])
                floating_species_names.append(item["name"])

    return structure, floating_species_names


def aiida_struct_to_ase(aiida_struct):
    """
    Custom version to bypass the inappropriate implementation of the site.get_ase() routine.

    That routine does not set tags for sites whose names coincide with the atomic symbol...).
    Here we always set the  tag to the "species number", which is the "kind number" + 1.
    """
    import ase

    # Build a "species dictionary" mapping kind names to tags (starting at 1)
    _kinds = aiida_struct.kinds
    sp_index = {}
    for index, value in enumerate(_kinds):
        sp_index[value.name] = index + 1

    s_ase = ase.Atoms(cell=aiida_struct.cell, pbc=aiida_struct.pbc)

    for site in aiida_struct.sites:
        ase_atom = site.get_ase(kinds=_kinds)
        ase_atom.tag = sp_index[site.kind_name]
        s_ase.append(ase_atom)

    return s_ase


def ase_struct_to_aiida(s_ase, kinds):
    """
    Converts an ASE structure object to an equivalent AiiDA object, preserving the kind names.

    :param: s_ase: The ASE object
    :param: kinds: The kinds object of a reference AiiDA structure
    """

    from aiida.orm import StructureData

    aiida_struct = StructureData(cell=s_ase.cell)

    positions = s_ase.positions
    tags = s_ase.get_tags()

    for index, value in enumerate(positions):
        kind = kinds[tags[index] - 1]
        aiida_struct.append_atom(position=value, symbols=kind.symbol, name=kind.name)

    return aiida_struct
