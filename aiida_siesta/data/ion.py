"""
This module manages the .ion.xml files in the local repository.
"""

from aiida.common.files import md5_file, md5_from_filelike
from aiida.orm.nodes import SinglefileData
from aiida.common.exceptions import StoringNotAllowed
from aiida_siesta.utils.pao_manager import PaoManager


def xml_element_to_string(element, tail=True):
    string = "<" + element.tag + ">" + element.text + "</" + element.tag + ">"
    if tail:
        string = string + element.tail
    return string


def parse_ion(fname):
    """
    Try to get relevant information from the .ion. For the moment, only the
    element symbol, name, mass and atomic number.
    Raise a ParsingError exception if the file does not contain the element symbol
    or the atomic number. The presence in the file of mass and name is, instead, considered
    optional. If not present, None is returned.
    """
    from xml.etree.ElementTree import ElementTree
    from aiida.common.exceptions import ParsingError
    from aiida.orm.nodes.data.structure import _valid_symbols

    parsed_data = {}

    el_tr = ElementTree(None, fname)
    root = el_tr.getroot()

    if root.find('symbol') is None or root.find('z') is None or root.find('label') is None:
        raise ParsingError(f"Currupted ion file {fname}: element symbol or atomic number missing")

    parsed_data["element"] = str(root.find('symbol').text.strip())
    if parsed_data["element"] not in _valid_symbols:
        raise ParsingError(f"Unknown element symbol {parsed_data['element']} in file {fname}")

    parsed_data["name"] = str(root.find('label').text.strip())
    parsed_data["atomic_number"] = int(root.find('z').text)
    parsed_data["mass"] = float(root.find('mass').text)

    return parsed_data


class IonData(SinglefileData):
    """
    Handler for ion files
    """

    def set_file(self, file_abs_path, filename=None):  #pylint: disable=arguments-differ
        """
        This is called in the __init__ of SingleFileData
        """
        # print("Called set_file","type of filename:",type(filename))
        parsed_data = parse_ion(file_abs_path)
        md5 = md5_file(file_abs_path)

        super().set_file(file_abs_path, filename)

        self.set_attribute('element', parsed_data["element"])
        self.set_attribute('name', parsed_data["name"])
        self.set_attribute('atomic_number', parsed_data["atomic_number"])
        if parsed_data["mass"] is not None:
            self.set_attribute('mass', parsed_data["mass"])
        self.set_attribute('md5', md5)

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
            self.validate_md5(self.md5)
        except ValueError as exception:
            raise StoringNotAllowed(exception) from exception

        try:
            self.validate_others_atts(self.element, self.name, self.atomic_number)
        except ValueError as exception:
            raise StoringNotAllowed(exception) from exception

        return super().store(*args, **kwargs)

    def validate_others_atts(self, elem, name, atm_n):
        """
        Validate the given element, name, atomic_number are the one of the stored file.
        Unfortunately it requires to reparse the file.
        :param elem: the symbol of the element.
               name: the name assigned to the atom/site.
               atm_n: the atomic number of the atom/site.
        :raises ValueError: if the element symbol is invalid.
        """
        with self.open(mode='r') as handle:
            parsed_data = parse_ion(handle.name)
        if elem != parsed_data["element"] or name != parsed_data["name"] or atm_n != parsed_data["atomic_number"]:
            raise ValueError(
                'element, name or atomic_number do not correspond to the the one in the ion file. '
                'The attributes of this class can not be modified manually.'
            )

    def validate_md5(self, md5: str):
        """
        Validate that the md5 checksum matches that of the currently stored file.
        :param value: the md5 checksum.
        :raises ValueError: if the md5 does not match that of the currently stored file.
        """
        with self.open(mode='rb') as handle:
            md5_fil = md5_from_filelike(handle)
            if md5 != md5_fil:
                raise ValueError(
                    f'Th md5 does not match that of stored file: {md5} != {md5_fil}. '
                    'The attributes of this class can not be modified manually.'
                )

    @classmethod
    def get_or_create(cls, file_abs_path, filename=None, use_first=False, store_ion=False):
        """
        Pass the same parameter of the __init__; if a file with the same md5
        is found, that IonData is returned, otherwise a new IonFile instance is
        created. Optionally, the new instance can be directly stored.

        :param file_abs_path: an absolute path file on disk.
        :param filename: optional explicit filename to give to the file stored in the repository.
                         Ignored if a file with the same md5 has been found.
        :param use_first: if False (default), raise an exception if more than one ion with same
                          md5 is found. If it is True, instead, use the first available ion.
        :param bool store_ion: if False (default), the IonData object is not stored in the database.
                               If True, the `store()` method is called and the returned IonData is stored.
        :return (ion, created): where ion is the IonData object, and create is either True if the
                                object was created, or False if the object was retrieved from the DB.
        """
        import os

        if not os.path.abspath(file_abs_path):
            raise ValueError("filename must be an absolute path")

        md5 = md5_file(file_abs_path)

        ions = cls.from_md5(md5)
        if not ions:
            instance = cls(file=file_abs_path, filename=filename)
            if store_ion:
                instance.store()
            return (instance, True)

        if len(ions) > 1:
            if use_first:
                return (ions[0], False)

            all_in_string = ",".join([str(i.pk) for i in ions])
            raise ValueError(
                f"More than one copy of a ion file with the same MD5 has been found in the DB. pks={all_in_string}. "
                "Set argument `use_first` to True to return a random one betweeen them."
            )

        return (ions[0], False)

    @classmethod
    def from_md5(cls, md5):
        """
        Return a list of all ions files that match a given MD5 hash.

        Note that the hash has to be stored in a _md5 attribute, otherwise
        the pseudo will not be found.
        """
        from aiida.orm import QueryBuilder
        qb = QueryBuilder()
        qb.append(cls, filters={'attributes.md5': {'==': md5}})
        return [_ for [_] in qb.all()]

    @property
    def element(self):
        return self.get_attribute('element', None)

    @property
    def name(self):
        return self.get_attribute('name', None)

    @property
    def mass(self):
        return self.get_attribute('mass', None)

    @property
    def atomic_number(self):
        return self.get_attribute('atomic_number', None)

    @property
    def md5(self):
        return self.get_attribute('md5', None)

    def get_content_ascii_format(self):
        """
        from the content, write the ld format .ion file. Necessary since siesta reads it.
        """
        from xml.etree import ElementTree

        root = ElementTree.fromstring(self.get_content())
        string = ""

        preamble_el = root.find("preamble")

        string = string + "<" + preamble_el.tag + ">" + preamble_el.text
        string = string + xml_element_to_string(preamble_el[0])  #basis
        string = string + xml_element_to_string(preamble_el[1])  #pseudo_header
        string = string + "</" + preamble_el.tag + ">\n"
        string = string + root.find("symbol").text + "\n"
        string = string + root.find("label").text + "\n"
        string = string + root.find("z").text + "\n"
        string = string + root.find("valence").text + "\n"
        string = string + root.find("mass").text + "\n"
        string = string + root.find("self_energy").text + "\n"
        string = string + root.find("lmax_basis").text + root.find("norbs_nl").text + "\n"
        string = string + root.find("lmax_projs").text + root.find("nprojs_nl").text + "#\n"

        string = string + "# PAOs:__________________________\n"
        for orbital in root.find("paos"):
            string = string + orbital.attrib["l"] + orbital.attrib["n"] + orbital.attrib["z"] + orbital.attrib[
                "ispol"] + orbital.attrib["population"] + "\n"
            radfunc = orbital.find("radfunc")
            string = string + radfunc.find("npts").text + radfunc.find("delta").text + radfunc.find("cutoff").text
            string = string + radfunc.find("data").text

        string = string + "# KBs:__________________________\n"
        for projector in root.find("kbs"):
            string = string + projector.attrib["l"] + projector.attrib["n"] + projector.attrib["ref_energy"] + "\n"
            radfunc = projector.find("radfunc")
            string = string + radfunc.find("npts").text + radfunc.find("delta").text + radfunc.find("cutoff").text
            string = string + radfunc.find("data").text

        string = string + "# Vna:__________________________\n"
        radfunc = root.find("vna").find("radfunc")
        string = string + radfunc.find("npts").text + radfunc.find("delta").text + radfunc.find("cutoff").text
        string = string + radfunc.find("data").text

        string = string + "# Chlocal:__________________________\n"
        radfunc = root.find("chlocal").find("radfunc")
        string = string + radfunc.find("npts").text + radfunc.find("delta").text + radfunc.find("cutoff").text
        string = string + radfunc.find("data").text

        core_info = radfunc = root.find("core")
        if core_info is not None:
            string = string + "# Core:__________________________\n"
            radfunc = root.find("core").find("radfunc")
            string = string + radfunc.find("npts").text + radfunc.find("delta").text + radfunc.find("cutoff").text
            string = string + radfunc.find("data").text

        return string

    def get_orbitals(self):
        """
        Uses sisl to read the file and return the orbitals
        """

        import sisl
        import os
        from aiida_siesta.data.atomic_orbitals import SislAtomicOrbital

        tmp_file = open("tmp.ion.xml", "w")
        tmp_file.write(self.get_content())
        tmp_file.close()
        sile = sisl.get_sile("tmp.ion.xml")
        sisl_atom = sile.read_basis()
        listorb = []
        for i in sisl_atom.orbitals:
            atorb = SislAtomicOrbital(i.name(), (i.orb.__getstate__()["r"], i.orb.__getstate__()["f"]), q0=i.q0)
            listorb.append(atorb)
        os.remove("tmp.ion.xml")

        return listorb

    def get_pao_modifier(self):
        """
        Get the PaoModifiers, that makes available methods for the modification of
        a PAO block.
        """
        pao_manager = PaoManager()
        pao_manager.set_from_ion(self)

        return pao_manager

    def get_pao_block(self):
        """
        Returns the PAO block correspondent to the orbitals stored in this class.
        """

        pao_manager = PaoManager()
        pao_manager.set_from_ion(self)

        return pao_manager.get_pao_block()

    def pao_size(self):
        """
        Returns a string sumarazing the size of the basis in this class. Following the siesta convencion
        (SZ, DZ, DZP ...)
        Please note that the string is composed  checking the maximum Z registred by orbitals, for both the polarized
        and unpolarixed orbitals. This means that the algorithm is not able to actually detect if the orbitals
        here are generated by a "PAO.BasisSize" in siesta or it is a manual PAO block. Take this method carefully.
        """
        pao_manager = PaoManager()
        pao_manager.set_from_ion(self)

        return pao_manager.pao_size()
