from aiida.plugins import DataFactory
from aiida.orm import KpointsData

from ..calculations.tkdict import FDFDict
from .base import SiestaBaseWorkChain
from ..utils.iterate_absclass import BaseIterator

# The following are helper functions to parse input values in the SiestaIterator. See
# the global dict SIESTA_ITERATION_PARAMS to know which parameters make use of them.


def set_up_parameters_dict(val, inputs, parameter, input_key, defaults=None):
    """
    Parsing function that sets an fdf parameter.

    This is used by the basis parameters (input_key='basis') and the rest of
    fdf parameters (input_key='parameters')

    Parameters
    -----------
    val: str, float, int, bool
        the value to set for the parameter.
    inputs: AttributeDict
        all the current inputs, so that we can extract the curren FDFdict.
    parameter: str
        the name of the fdf flag. Case and hyphen insensitive.
    input_key: str
        the input of the fdf dict that should be modified.
    defaults: dict, optional
        a dictionary of defaults. The only key that matters right now is "units",
        which is the units to set for the value if the value is a number. This may change.
    """

    val = getattr(val, "value", val)

    # Get the current FDFdict for the corresponding input
    parameters = getattr(inputs, input_key, DataFactory('dict')())
    parameters = FDFDict(parameters.get_dict())

    # Set the units for the value if needed
    if isinstance(val, (int, float)):
        if defaults is not None and "units" in defaults:
            val = f'{val} {defaults["units"]}'

    # Then set the value of the parameter in the FDF dict.
    parameters[parameter] = val

    # And then just translate it again to a dict to use it in the input
    return DataFactory('dict')(dict=parameters.get_dict())


def set_up_kpoint_grid(val, inputs, parameter, input_key='kpoints'):
    """
    Parsing function that sets a kpoint grid.

    This is used by the basis parameters (input_key='basis') and the rest of
    fdf parameters (input_key='parameters')

    Parameters
    -----------
    val: str, float, int, bool
        the value to set for the parameter.
    inputs: AttributeDict
        all the current inputs, so that we can extract the old kpoints.
    parameter: str, {'kpoints_density', 'kpoints_0', 'kpoints_1', 'kpoints_2'}
        used to understand how to parse the value.

        If 'kpoints_density': the value is interpreted as the maximum distance
        between two grid points along a reciprocal axis.

        Else, it is interpreted as the number of points for one of the components
        of the grid.
    input_key: str
        the input of the fdf dict that should be modified.
    """

    # If there is already a KpointsData() in inputs grab it to modify it
    old_kpoints = getattr(inputs, input_key, None)

    # Else define a new one
    if old_kpoints is None:
        old_kpoints = KpointsData()

    # Get the mesh and the offset
    try:
        mesh, offset = old_kpoints.get_kpoints_mesh()
    except (KeyError, AttributeError):
        mesh, offset = [1, 1, 1], [0, 0, 0]

    # Get also the cell
    if hasattr(old_kpoints, 'cell'):
        cell = old_kpoints.cell
        pbc = old_kpoints.pbc
    else:
        cell = inputs.structure.cell
        pbc = inputs.structure.pbc

    # And finally define the new KpointsData according to the required mode.

    # Change the density
    if parameter == 'kpoints_density':
        new_kpoints = KpointsData()
        new_kpoints.set_cell(cell)
        new_kpoints.pbc = pbc
        new_kpoints.set_kpoints_mesh_from_density(val.value, offset=offset)

    # Change an individual component
    else:

        component = int(parameter[-1])

        mesh[component] = val

        new_kpoints = KpointsData()
        new_kpoints.set_kpoints_mesh(mesh, offset)

    return new_kpoints


# This is the parameters' look up list for the siesta iterator, which enables iterating
# over extra parameters apart from the inputs of SiestaBaseWorkChain. This may be taken
# as an example to allow extra parameters in any input iterator that uses a different
# process.

SIESTA_ITERATION_PARAMS = ({
    "group_name":
    "Basis parameters",
    "help":
    """
    Supports iterating over parameters of the fdf that are related to the basis.
    """,
    "input_key":
    "basis",
    "parse_func":
    set_up_parameters_dict,
    "condition":
    lambda parameter: FDFDict.translate_key(parameter).startswith("pao"),
    "keys":
    FDFDict({
        "paobasissize": {
            'defaults': {
                'values_list': ['SZ', 'SZP', 'DZ', 'DZP', 'TZ', 'TZP']
            }
        },
        "paoenergyshift": {
            'defaults': {
                'units': 'Ry'
            }
        }
    })
}, {
    "group_name": "SCF Brillouin zone",
    "input_key": "kpoints",
    "help": """
    Group of parameters that support modifying the kpoint sampling for the SCF cycles.
    """,
    "parse_func": set_up_kpoint_grid,
    "keys": {
        'kpoints_density': None,
        'kpoints_0': None,
        'kpoints_1': None,
        'kpoints_2': None
    }
}, {
    "group_name":
    "FDF parameters",
    "input_key":
    "parameters",
    "help":
    """
    All the parameters that are used "raw" in the fdf file.

    WARNING: Even if the parameter does not make sense, if it has not been interpreted as belonging
    to any other group, it will be interpreted as an fdf parameter.
    """,
    "condition":
    lambda parameter: True,
    "parse_func":
    set_up_parameters_dict,
    "keys":
    FDFDict({"meshcutoff": {
        'defaults': {
            'units': 'Ry',
            'init_value': 100,
            'step': 100
        }
    }})
})


class SiestaIterator(BaseIterator):
    """
    Iterator for the SietaBaseWorkChain. The iterator is extended to iterate over any Siesta keyword.
    WARNING: if a keyword not recognized by Siesta is used in `iterate_over`, the iterator will not
    complain. It will just add the keyword to the parameters dict and run the calculation!
    """

    _process_class = SiestaBaseWorkChain
    _expose_inputs_kwargs = {'exclude': ('metadata',)}
    _params_lookup = SIESTA_ITERATION_PARAMS
