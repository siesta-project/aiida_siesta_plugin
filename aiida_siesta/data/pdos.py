"""
This module manages the .pdos files in the local repository.
"""

from aiida.orm.nodes import SinglefileData
from aiida.common.exceptions import StoringNotAllowed


def parse_pdos(fname):
    """
    Try to get few relevant information from the pdos. For the moment, only the
    element symbol, name, mass and atomic number.
    Raise a ParsingError exception if the file does not contain the element symbol
    or the atomic number. The presence in the file of mass and name is, instead, considered
    optional. If not present, None is returned.
    """
    from xml.etree.ElementTree import ElementTree
    from aiida.common.exceptions import ParsingError

    parsed_data = {}

    el_tr = ElementTree(None, fname)
    root = el_tr.getroot()

    if root.find('nspin') is None or root.find('norbitals') is None or root.find('fermi_energy') is None:
        raise ParsingError(f"Currupted PDOS file {fname}: nspin, norbitals or fermi_energy missing")

    parsed_data["n_spin"] = int(root.find('nspin').text)
    parsed_data["n_orbitals"] = int(root.find('norbitals').text)
    parsed_data["fermi_energy"] = float(root.find('fermi_energy').text)

    parsed_data["min_energy"] = float(root.find('energy_values').text.split()[0])
    parsed_data["max_energy"] = float(root.find('energy_values').text.split()[-1])
    parsed_data["bins"] = len(root.find('energy_values').text.split())

    return parsed_data


class PdosData(SinglefileData):
    """
    Handler for pdos files
    """

    def set_file(self, file_abs_path, filename=None):  #pylint: disable=arguments-differ
        """
        This is called in the __init__ of SingleFileData
        """
        # print("Called set_file","type of filename:",type(filename))
        parsed_data = parse_pdos(file_abs_path)
        #md5 = md5_file(file_abs_path)

        super().set_file(file_abs_path, filename)

        self.set_attribute('n_spin', parsed_data["n_spin"])
        self.set_attribute('n_orbitals', parsed_data["n_orbitals"])
        self.set_attribute('fermi_energy', parsed_data["fermi_energy"])
        self.set_attribute('min_energy', parsed_data["min_energy"])
        self.set_attribute('max_energy', parsed_data["max_energy"])
        self.set_attribute('bins', parsed_data["bins"])

    def store(self, *args, **kwargs):  # pylint: disable=arguments-differ
        """
        Store the node. It requires a previous check on the assigned attributes.
        In fact, the attributes of this particular class must just reflect the info
        in the ion file. However the design of `Data` class do not allow to
        make attributes immutable before storing and, therefore, a crazy user
        might think to change them before storing.
        Here we check that the attributes actually corresponds to the file info.
        """

        if self.is_stored:
            return self

        try:
            self.validate_atts(
                self.n_spin, self.n_orbitals, self.fermi_energy, self.min_energy, self.max_energy, self.bins
            )
        except ValueError as exception:
            raise StoringNotAllowed(exception) from exception

        return super().store(*args, **kwargs)

    def validate_atts(self, nspin, norbs, ef, min_e, max_e, bins):  #pylint: disable=invalid-name
        """
        Validate the given nspin, norbs, ef, min_e, max_e, bins are the one of the stored file.
        Unfortunately it requires to reparse the file.
        """
        with self.open(mode='r') as handle:
            parsed_data = parse_pdos(handle.name)
        if nspin != parsed_data["n_spin"] or norbs != parsed_data["n_orbitals"] or ef != parsed_data["fermi_energy"]:
            raise ValueError(
                'n_spin, n_orbitals or fermi_energy do not correspond to the one in the pdos file. '
                'The attributes of this class can not be modified manually.'
            )
        if min_e != parsed_data["min_energy"] or max_e != parsed_data["max_energy"] or bins != parsed_data["bins"]:
            raise ValueError(
                'min_energy, max_energy or bins do not correspond to the one in the pdos file. '
                'The attributes of this class can not be modified manually.'
            )

    @property
    def n_spin(self):
        return self.get_attribute('n_spin', None)

    @property
    def n_orbitals(self):
        return self.get_attribute('n_orbitals', None)

    @property
    def fermi_energy(self):
        return self.get_attribute('fermi_energy', None)

    @property
    def min_energy(self):
        return self.get_attribute('min_energy', None)

    @property
    def max_energy(self):
        return self.get_attribute('max_energy', None)

    @property
    def bins(self):
        return self.get_attribute('bins', None)

    def get_energy_array(self):

        import sisl
        import os

        tmp_file = open("tmp.PDOS", "w")
        tmp_file.write(self.get_content())
        tmp_file.close()
        sile = sisl.get_sile("tmp.PDOS")
        sisl_data = sile.read_data()
        os.remove("tmp.PDOS")

        return sisl_data[1]

    def get_pdoses_array(self):

        import sisl
        import os

        tmp_file = open("tmp.PDOS", "w")
        tmp_file.write(self.get_content())
        tmp_file.close()
        sile = sisl.get_sile("tmp.PDOS")
        sisl_data = sile.read_data()
        os.remove("tmp.PDOS")

        return sisl_data[2]
