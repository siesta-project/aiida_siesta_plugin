from aiida.parsers import Parser
from aiida.orm import Dict
from aiida.common import OutputParsingError
from aiida.common import exceptions

# See the LICENSE.txt and AUTHORS.txt files.


class STMOutputParsingError(OutputParsingError):
    pass


#---------------------------


class STMParser(Parser):
    """
    Parser for the output of the "plstm" program in the Siesta distribution.
    """

    _version = "1.1.1"

    def parse(self, **kwargs):
        """
        Receives in input a dictionary of retrieved nodes.
        Does all the logic here.
        """
        from aiida.engine import ExitCode

        try:
            output_folder = self.retrieved
        except exceptions.NotExistent:
            return self.exit_codes.ERROR_NO_RETRIEVED_FOLDER

        filename_plot = None
        for element in output_folder.list_object_names():
            if ".STM" in element:
                filename_plot = element

        if filename_plot is None:
            return self.exit_codes.ERROR_OUTPUT_PLOT_MISSING

        old_version = True
        if "LDOS" in filename_plot:
            old_version = False
        if old_version and self.node.inputs.spin_option.value != "q":
            self.node.logger.warning(
                "Spin option different from 'q' was requested in "
                "input, but the plstm version used for the calculation do not implement spin support. "
                "The parsed STM data refers to the total charge option."
            )

        try:
            plot_contents = output_folder.get_object_content(filename_plot)
        except (IOError, OSError):
            return self.exit_codes.ERROR_OUTPUT_PLOT_READ

        # Save grid_X, grid_Y, and STM arrays in an ArrayData object
        try:
            stm_data = get_stm_data(plot_contents)
        except (IOError, OSError):
            return self.exit_codes.ERROR_CREATION_STM_ARRAY

        self.out('stm_array', stm_data)

        parser_info = {}
        parser_info['parser_info'] = 'AiiDA STM(Siesta) Parser V. {}'.format(self._version)
        parser_info['parser_warnings'] = []
        parser_info['output_data_filename'] = filename_plot

        # Possibly put here some parsed data from the stm.out file
        # (but it is not very interesting)
        result_dict = {}

        # Add parser info dictionary
        parsed_dict = dict(list(result_dict.items()) + list(parser_info.items()))

        output_data = Dict(dict=parsed_dict)
        self.out('output_parameters', output_data)

        return ExitCode(0)


def get_stm_data(plot_contents):
    """
    Parses the STM plot file to get an Array object with
    X, Y, and Z arrays in the 'meshgrid'
    setting, as in the example code:

    import numpy as np
    xlist = np.linspace(-3.0, 3.0, 3)
    ylist = np.linspace(-3.0, 3.0, 4)
    X, Y = np.meshgrid(xlist, ylist)
    Z = np.sqrt(X**2 + Y**2)

    X:
    [[-3.  0.  3.]
    [-3.  0.  3.]
    [-3.  0.  3.]
    [-3.  0.  3.]]

    Y:
    [[-3. -3. -3.]
    [-1. -1. -1.]
    [ 1.  1.  1.]
    [ 3.  3.  3.]]

    Z:
    [[ 4.24264069  3.          4.24264069]
    [ 3.16227766  1.          3.16227766]
    [ 3.16227766  1.          3.16227766]
    [ 4.24264069  3.          4.24264069]]

    These can then be used in matplotlib to get a contour plot.

    :param plot_contents: the contents of the *.STM file as a string
    :return: `aiida.orm.ArrayData` instance representing the STM contour.
    """

    import numpy as np
    from itertools import groupby
    from aiida.orm import ArrayData

    # aiida.CH.STM or aiida.CC.STM...
    data = plot_contents.split('\n')
    data = [i.split() for i in data]

    # The data in the file is organized in "lines" parallel to the Y axes
    # (that is, for constant X) separated by blank lines.
    # In the following we use the 'groupby' function to get at the individual
    # blocks one by one, and set the appropriate arrays.

    # I am not sure about the mechanics of groupby,
    # so repeat
    xx = []  #pylint: disable=invalid-name
    yy = []  #pylint: disable=invalid-name
    zz = []  #pylint: disable=invalid-name
    #
    # Function to separate the blocks
    h = lambda x: len(x) == 0  #pylint: disable=invalid-name
    #
    for k, v in groupby(data, h):
        if not k:
            xx.append([i[0] for i in v])
    for k, v in groupby(data, h):
        if not k:
            yy.append([i[1] for i in v])
    for k, v in groupby(data, h):
        if not k:
            zz.append([i[2] for i in v])

    # Now, transpose, since x runs fastest in our fortran code,
    # the opposite convention of the meshgrid paradigm.

    arrayx = np.array(xx, dtype=float).transpose()
    arrayy = np.array(yy, dtype=float).transpose()
    arrayz = np.array(zz, dtype=float).transpose()

    arraydata = ArrayData()
    arraydata.set_array('grid_X', np.array(arrayx))
    arraydata.set_array('grid_Y', np.array(arrayy))
    arraydata.set_array('STM', np.array(arrayz))

    return arraydata
