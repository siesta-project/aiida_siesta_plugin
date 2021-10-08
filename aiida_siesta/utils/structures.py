def clone_aiida_structure(structure):
    """
    A cloned structure is not quite ready to store more atoms.
    This function fixes it
    """

    truct = structure.clone()
    truct._internal_kind_tags = {}

    return truct


def add_ghost_sites_to_structure(original_structure, basis):
    """
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
    This is a custom version to bypass the inappropriate implementation
    of the site.get_ase() routine (it does not set tags for sites whose
    names coincide with the atomic symbol...). Here we always set the
    tag to the "species number", which is the "kind number" + 1.
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
    Converts an ASE structure object to an equivalent AiiDA object,
    preserving the kind names.
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


#def exchange_sites_in_structure(s, i1, i2):
#    """
#    Given a structure s, return another structure with the coordinates of hte i1, i2 sites interchanged
#    :param: s  a StructureData object
#    :param: i1, i2     site indexes to be exchanged
#    """
#    from aiida.orm.nodes.data.structure import Site
#
#    sites = s.sites
#
#    s1 = sites[i1]
#    s2 = sites[i2]
#    p1 = s1.position
#    p2 = s2.position
#    n1 = Site(kind_name=s1.kind_name, position=p2)
#    n2 = Site(kind_name=s2.kind_name, position=p1)
#    sites[i1] = n1
#    sites[i2] = n2
#
#    t = clone_aiida_structure(s)
#    t.clear_sites()
#
#    for site in sites:
#        t.append_site(site)
#
#    return t
#
#
#def find_intermediate_structure(s, i1, intermediate_position, i2=None):
#    """
#    Given a structure and the indexes of two sites, and an
#    intermediate point, generate an intermediate structure in which
#    the i1 atom has moved to the intermediate point. If a second atom
#    index is input, the second atom is moved to an image point with an
#    opposite relative displacement. This is useful to 'prime' a path
#    to avoid head-on collissions, for example.
#
#    :param: s  a StructureData object
#    :param: i1, i2,  indexes of atoms in the structure s
#            If i2 is None (default), only the first atom is moved.
#    :param: intermediate_position, a list of floats
#
#    """
#
#    import numpy as np
#    from aiida.orm.nodes.data.structure import Site
#
#    i1_path_position = np.array(intermediate_position)
#
#    sites = s.sites
#
#    p1 = np.array(sites[i1].position)
#    s1 = sites[i1]
#    n1 = Site(kind_name=s1.kind_name, position=i1_path_position)
#    sites[i1] = n1
#
#    # The second atom's image position is obtained by
#    # reversing the sign of the above relative position for i1
#    if i2 is not None:
#        p2 = np.array(sites[i2].position)
#        # Relative position of the path point and the first atom
#        p_wrt_p1 = i1_path_position - p1
#        i2_path_position = p2 - p_wrt_p1
#        s2 = sites[i2]
#        n2 = Site(kind_name=s2.kind_name, position=i2_path_position)
#        sites[i2] = n2
#
#    intermediate_structure = clone_aiida_structure(s)
#    intermediate_structure.clear_sites()
#    for site in sites:
#        intermediate_structure.append_site(site)
#
#    return intermediate_structure
#
#
#def compute_mid_path_position(s, i1, i2, migration_direction):
#    """
#    The basic heuristic here is to avoid head-on collissions
#    by defining an "avoidance cylinder" around the line
#    joining the two atoms exchanged. The input "migration_direction"
#    serves to define a point on the surface of that cylinder, at
#    the mid-point, which is used as the mid-point of the starting path.
#
#    :param: s  a StructureData object
#    :param: i1, i2,  indexes of the atoms
#    :param: migration direction, in lattice coordinates
#    """
#
#    import numpy as np
#
#    AVOIDANCE_RADIUS = 1.00  # 1.0 angstrom
#
#    cell = np.array(s.cell)
#    cell_direction = np.array(migration_direction)
#
#    cart_direction = np.matmul(cell, cell_direction)
#
#    # Find positions of i1 and i2
#    sites = s.sites
#    p1 = np.array(sites[i1].position)
#    p2 = np.array(sites[i2].position)
#
#    #
#    # Find the unit vector parallel to the line i1-i2
#    #
#    d = p2 - p1
#    dmod = np.sqrt(d.dot(d))
#    u = d / dmod
#
#    # Sanity check: migration direction should not be near-parallel
#    # to the line joining the exchanged atoms...
#
#    cross_product = np.cross(d, cart_direction)
#    mod_cross_product = np.sqrt(cross_product.dot(cross_product))
#    mod_cd = np.sqrt(np.dot(cart_direction, cart_direction))
#
#    if np.abs(mod_cross_product / (mod_cd * dmod)) < 1.0e-2:
#        print("Migration direction near parallel to line of sight")
#        return None
#
#    # Find component of cart_direction perpendicular to the i1-i2 line,
#    # and unit vector
#    #
#    c_perp = cart_direction - u.dot(cart_direction) * u
#    c_perp_mod = np.sqrt(c_perp.dot(c_perp))
#    u_perp = c_perp / c_perp_mod
#
#    #
#    # The mid-point of the path is now determined by the vector sum
#    # of half of d and u_perp times the radius of the avoidance cylinder
#    #
#    path_mid_point = p1 + 0.5 * dmod * u + AVOIDANCE_RADIUS * u_perp
#
#    return path_mid_point.tolist()
#
#
#def find_mid_path_position(s, pos1, pos2, migration_direction):
#    """
#    The basic heuristic here is to avoid head-on collissions
#    by defining an "avoidance cylinder" around the line
#    joining the initial and final points . The input "migration_direction"
#    serves to define a point on the surface of that cylinder, at
#    the mid-point, which is used as the mid-point of the starting path.
#
#    :param: s  a StructureData object
#    :param: pos1, pos2,  initial and final positions
#    :param: migration direction, in lattice coordinates
#    """
#
#    import numpy as np
#
#    AVOIDANCE_RADIUS = 1.00  # 1.0 angstrom
#
#    cell = np.array(s.cell)
#    cell_direction = np.array(migration_direction)
#
#    cart_direction = np.matmul(cell, cell_direction)
#
#    # Find positions of i1 and i2
#    sites = s.sites
#    p1 = np.array(pos1)
#    p2 = np.array(pos2)
#
#    #
#    # Find the unit vector parallel to the line i1-i2
#    #
#    d = p2 - p1
#    dmod = np.sqrt(d.dot(d))
#    u = d / dmod
#
#    # Sanity check: migration direction should not be near-parallel
#    # to the line joining the two sites
#
#    cross_product = np.cross(d, cart_direction)
#    mod_cross_product = np.sqrt(cross_product.dot(cross_product))
#    mod_cd = np.sqrt(np.dot(cart_direction, cart_direction))
#
#    if np.abs(mod_cross_product / (mod_cd * dmod)) < 1.0e-2:
#        print("Migration direction near parallel to line of sight")
#        return None
#
#    # Find component of cart_direction perpendicular to the i1-i2 line,
#    # and unit vector
#    #
#    c_perp = cart_direction - u.dot(cart_direction) * u
#    c_perp_mod = np.sqrt(c_perp.dot(c_perp))
#    u_perp = c_perp / c_perp_mod
#
#    #
#    # The mid-point of the path is now determined by the vector sum
#    # of half of d and u_perp times the radius of the avoidance cylinder
#    #
#    path_mid_point = p1 + 0.5 * dmod * u + AVOIDANCE_RADIUS * u_perp
#
#    return path_mid_point.tolist()
