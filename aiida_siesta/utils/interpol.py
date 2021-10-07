#
# Module
#
from aiida_siesta.utils.structures import aiida_struct_to_ase
from aiida_siesta.utils.structures import ase_struct_to_aiida


#
def interpolate_two_structures(s1, s2, n_images):
    """
    Interpolate linearly the coordinates of two structures.
    Assume StructureData objects s1 and s2 are commensurate.
    :param: n_images is the number of internal points ("images")
            over which to interpolate

    :return: a list of structures
    """

    import numpy as np
    p1 = np.array([site.position for site in s1.sites])
    p2 = np.array([site.position for site in s2.sites])

    structure_list = [s1]
    delta = (p2 - p1) / (n_images + 1)
    for i in range(n_images):
        pi = p1 + (i + 1) * delta
        # Note that this is a simple clone.
        # No more atoms can be added to this structure
        # Be on guard for possible issues with 'ghost' atoms later
        si = s1.clone()
        #
        # If problems are encountered use
        # si._internal_kind_tags={}
        #
        si.reset_sites_positions(pi)
        structure_list.append(si)

    structure_list.append(s2)

    return structure_list


def interpolate_two_structures_ase(s1, s2, n_images, interp_method="idpp"):
    """
    Uses the NEB pre-optimizer to interpolate two structures
    :param: s1  initial AiiDA structure
    :param: s2  final AiiDA structure
    :param: n_images is the number of internal points ("images")
            over which to interpolate
    :param: interp_method: A string with the interpolation method
    """

    from ase.neb import NEB

    ase1 = aiida_struct_to_ase(s1)
    ase2 = aiida_struct_to_ase(s2)

    tags = ase2.get_tags()

    images = [ase1]
    for i in range(int(n_images)):
        images.append(ase1.copy())
        # 'copy' to get extra real slots

    images.append(ase2)

    neb = NEB(images)
    neb.interpolate(method=interp_method)

    slist = []
    # Structure list
    for i in range(n_images + 2):
        #
        # Careful with this. The NEB module does not
        # seem to maintain the tags
        #
        neb.images[i].set_tags(tags)
        s = ase_struct_to_aiida(neb.images[i], s1.kinds)
        slist.append(s)

    return slist
