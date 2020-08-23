from abc import ABC, abstractmethod

from aiida.plugins import DataFactory
from aiida.engine import calcfunction

class ParamPlugin(ABC):

    @classmethod
    @abstractmethod
    def get(cls, input_key):
        pass

@calcfunction
def rescale(structure, scale):
    """
    Calcfunction to rescale a structure by a scaling factor.
    Uses ase.
    :param structure: An AiiDA structure to rescale
    :param scale: The scale factor
    :return: The rescaled structure
    """

    the_ase = structure.get_ase()
    new_ase = the_ase.copy()
    new_ase.set_cell(the_ase.get_cell() * float(scale), scale_atoms=True)
    new_structure = DataFactory('structure')(ase=new_ase)

    return new_structure

@calcfunction
def scale_to_vol(structure, vol):
    """
    Calcfunction to scale a structure to a target volume.
    Uses pymatgen.
    :param stru: An aiida structure
    :param vol: The target volume per atom in angstroms
    :return: The new scaled AiiDA structure
    """

    in_structure = structure.get_pymatgen_structure()
    new = in_structure.copy()
    new.scale_lattice(float(vol) * in_structure.num_sites)
    StructureData = DataFactory("structure")
    structure = StructureData(pymatgen_structure=new)

    return structure

ops = {
    "scale": rescale,
    "scales": rescale,
    "volscale": scale_to_vol
}

def modify_structure(val, inputs, parameter, input_key="structure"):

    op = parameter.split("_")[-1]

    structure = getattr(inputs, input_key)

    modified_struct = ops[op](structure, val)

    return modified_struct

class StructureModificationsPlugin(ParamPlugin):

    _group_name = "Structure modifications"

    @classmethod
    def get(cls, input_key="structure"):

        group = {
            "group_name": cls._group_name,
            "input_key": input_key,
            "parse_func": modify_structure,
            "keys": {"scale": None, "scales": None}
        }

        group["keys"] = {f"{input_key}_{key}": val for key, val in group["keys"].items()}

        return group

    def get_operator(cls, key, default=None):
        return ops.get(key, default)

def get_plugin(name, **kwargs):
    for cls in ParamPlugin.__subclasses__():
        if cls._group_name == name:
            return cls.get(**kwargs)
