"""
This module manages the .ion.xml files in the local repository.
"""

from aiida.common.files import md5_file
from aiida.orm.nodes import SinglefileData


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

    @classmethod
    def get_or_create(cls, filename, use_first=False, store_ion=True):
        """
        Pass the same parameter of the init; if a file with the same md5
        is found, that IonData is returned.

        :param filename: an absolute filename on disk
        :param use_first: if False (default), raise an exception if more than \
                one potential is found.\
                If it is True, instead, use the first available pseudopotential.
        :param bool store_ion: If false, the IonData objects are not stored in
                the database. default=True.
        :return (ion, created): where ion is the IonData object, and create is either\
            True if the object was created, or False if the object was retrieved\
            from the DB.
        """
        import os

        if not os.path.abspath(filename):
            raise ValueError("filename must be an absolute path")

        md5 = md5_file(filename)

        ions = cls.from_md5(md5)
        if not ions:
            instance = cls(file=filename)
            print(instance.element)
            if store_ion:
                instance.store()
            return (instance, True)

        if len(ions) > 1:
            if use_first:
                return (ions[0], False)

            raise ValueError(
                "More than one copy of a ion file "
                "with the same MD5 has been found in the "
                "DB. pks={}".format(",".join([str(i.pk) for i in ions]))
            )

        return (ions[0], False)

    def store(self, *args, **kwargs):  # pylint: disable=arguments-differ
        """
        Store the node, reparsing the file so that the md5 and the element
        are correctly reset. The reason is because there is nothing that
        prevents the user to modify attributes manually before storing.
        However the attributes should only reflect info in the file!
        Probably I will change the approach to a verification system.
        """
        from aiida.common.files import md5_from_filelike

        if self.is_stored:
            return self

        # Already done??
        with self.open(mode='r') as handle:
            parsed_data = parse_ion(handle.name)

        # Open in binary mode which is required for generating the md5 checksum
        with self.open(mode='rb') as handle:
            md5sum = md5_from_filelike(handle)

        self.set_attribute('element', parsed_data["element"])
        self.set_attribute('name', parsed_data["name"])
        if parsed_data["mass"] is not None:
            self.set_attribute('mass', parsed_data["mass"])
        self.set_attribute('atomic_number', parsed_data["atomic_number"])
        self.set_attribute('md5', md5sum)

        return super().store(*args, **kwargs)

    @classmethod
    def from_md5(cls, md5):
        """
        Return a list of all ions files that match a given MD5 hash.

        Note that the hash has to be stored in a _md5 attribute, otherwise
        the pseudo will not be found.
        """
        from aiida.orm.querybuilder import QueryBuilder
        qb = QueryBuilder()
        qb.append(cls, filters={'attributes.md5': {'==': md5}})
        return [_ for [_] in qb.all()]

    def set_file(self, filename):  # pylint: disable=arguments-differ
        """
        This is called in the __init__ of SingleFileData
        """
        # print("Called set_file","type of filename:",type(filename))
        parsed_data = parse_ion(filename)
        md5sum = md5_file(filename)

        super().set_file(filename)

        self.set_attribute('element', parsed_data["element"])
        self.set_attribute('name', parsed_data["name"])
        self.set_attribute('atomic_number', parsed_data["atomic_number"])
        if parsed_data["mass"] is not None:
            self.set_attribute('mass', parsed_data["mass"])
        self.set_attribute('md5', md5sum)

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
    def md5sum(self):
        return self.get_attribute('md5', None)

    def _validate(self):
        from aiida.common.exceptions import ValidationError
        from aiida.common.files import md5_from_filelike

        super()._validate()

        # Yet another parsing ???
        with self.open(mode='r') as handle:
            parsed_data = parse_ion(handle.name)

        # Open in binary mode which is required for generating the md5 checksum
        with self.open(mode='rb') as handle:
            md5 = md5_from_filelike(handle)

        # This is erroneous exception,
        # as it is in the `upf` module oin `aiida_core`
        try:
            element = parsed_data['element']
        except KeyError:
            raise ValidationError("No 'element' could be parsed in the PSML " "file {}".format(self.filename))

        try:
            attr_element = self.get_attribute('element')
        except AttributeError:
            raise ValidationError("attribute 'element' not set.")

        try:
            attr_md5 = self.get_attribute('md5')
        except AttributeError:
            raise ValidationError("attribute 'md5' not set.")

        if attr_element != element:
            raise ValidationError(
                "Attribute 'element' says '{}' but '{}' was "
                "parsed instead.".format(attr_element, element)
            )

        if attr_md5 != md5:
            raise ValidationError("Attribute 'md5' says '{}' but '{}' was " "parsed instead.".format(attr_md5, md5))

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

    def _get_pao_dict(self):
        """
        Internal function that goes through the orbitals and extracts two dictionaries.
        return: gen_dictl, `dict`
            containing information about the non polarized orbitals (the one with attr "P" False)
        return: gen_dictl, `dict`
            containing information about the polarized orbitals (the one with attr "P" True)
        The structure of both dict is
        {n1 : { l1 : { z1 : r1, z2 : r2, ...}, l2 : ...} ...}
        where n* l* z* are integers carring the value of the n,l shells and z (please note we do not
        keep track of quantum number m as in the PAO construction is not possible to specify different
        feature for m shells with same l) and r* are the maximum radius for the n,l,z orbital.
        Maximum radius for the siesta orbitals means the radius after which the radial part is zero.
        """
        gen_dict = {}
        pol_dict = {}
        for orbital in self.get_orbitals():
            i = orbital.attributes
            if not i["P"]:
                if i["n"] not in gen_dict:
                    gen_dict[i["n"]] = {i["l"]: {i["Z"]: i["R"]}}
                else:
                    if i["l"] not in gen_dict[i["n"]]:
                        gen_dict[i["n"]][i["l"]] = {i["Z"]: i["R"]}
                    else:
                        if i["Z"] not in gen_dict[i["n"]][i["l"]]:
                            gen_dict[i["n"]][i["l"]][i["Z"]] = i["R"]
            else:
                if i["n"] not in pol_dict:
                    pol_dict[i["n"]] = {i["l"] - 1: {i["Z"]: i["R"]}}
                else:
                    if i["l"] - 1 not in pol_dict[i["n"]]:
                        pol_dict[i["n"]][i["l"] - 1] = {i["Z"]: i["R"]}
                    else:
                        if i["Z"] not in pol_dict[i["n"]][i["l"] - 1]:
                            pol_dict[i["n"]][i["l"] - 1][i["Z"]] = i["R"]

        #The polarization of 2p is 3d (2d does not exist).
        #Same for polarization of 1s e 3d.
        for num in [2, 3, 4]:
            if num in pol_dict and num not in gen_dict:
                pol_dict[num - 1] = pol_dict[num]
                pol_dict.pop(num)

        return gen_dict, pol_dict

    def get_pao_modifier(self):
        """
        Get the PaoModifiers, that makes available methods for the modification of
        a PAO block.
        """

        dictl, pol = self._get_pao_dict()

        return PaoModifier(self.name, dictl, pol)

    def get_pao_block(self):
        """
        Returns the PAO block correspondent to the orbitals stored in this class.
        """

        dictl, pol = self._get_pao_dict()

        return PaoModifier(self.name, dictl, pol).get_pao_block()

    def pao_size(self):
        """
        Returns a string sumarazing the size of the basis in this class. Following the siesta convencion
        (SZ, DZ, DZP ...)
        Please note that the string is composed  checking the maximum Z registred by orbitals, for both the polarized
        and unpolarixed orbitals. This means that the algorithm is not able to actually detect if the orbitals
        here are generated by a "PAO.BasisSize" in siesta or it is a manual PAO block. Take this method carefully.
        """
        max_z = 0
        max_pz = 0

        for i in self.attributes["orbital_dicts"]:
            if not i["P"]:
                if i["Z"] > max_z:
                    max_z = i["Z"]
            else:
                if i["Z"] > max_pz:
                    max_pz = i["Z"]

        dictma = ["S", "D", "T"]
        dictpo = ["", "D", "T"]

        string = "{}Z{}P".format(dictma[max_z - 1], dictpo[max_pz - 1])

        self.logger.warning(
            "Carefull. This information is an approximation. It is extracted looking "
            "at the maximum Z for polarized and unpolarized orbitals. Not all the siesta basis "
            "choises can be condensed into a simple string. Check the `get_pao_block` for "
            "more detailed information"
        )

        return string


class PaoModifier:
    """
    class to help modifications to PAO basis block. Requires the gen_dict and pol_dict
    as defined in `_get_pao_dict` of the class IonData and also the name.
    """

    def __init__(self, name, dictl, pol=None):
        """
        Sets attributs with the required inputs
        """
        self.name = name
        self.gen_dict = dictl
        self.pol_dict = pol

    def change_all_radius(self, percentage):
        """
        increment decrement of all orbitas of a percentage
        """
        dictl = self.gen_dict.copy()
        pol = self.pol_dict.copy()
        for i in dictl:
            for j in dictl[i]:
                for l in dictl[i][j]:
                    dictl[i][j][l] = dictl[i][j][l] + percentage / 100 * dictl[i][j][l]
        for i in pol:
            for j in pol[i]:
                for l in pol[i][j]:
                    pol[i][j][l] = pol[i][j][l] + percentage / 100 * pol[i][j][l]

        self.gen_dict = dictl
        self.pol_dict = pol

    def add_polarization(self, n, l):  #pylint: disable=invalid-name
        """
        add polarization to the orbital with quantum numbers n, l
        """
        try:
            self.gen_dict[n][l]
        except KeyError:
            raise ValueError("no orbital with n = {0} and l = {1} is present in the basis".format(n, l))

        if l in self.pol_dict[n]:
            num_z = len(self.pol_dict[n][l])
            self.pol_dict[n][l][num_z + 1] = 0.000
        else:
            self.pol_dict[n][l] = {"1": self.gen_dict[n][l][1]}

    def get_pao_block(self):
        """
        From the info of the gen_dict, pol_dict, creates the PAO block.
        return a string card containing the block.
        """
        ang_to_bohr = 1.8897161646321
        number_of_l = 0
        dictl = self.gen_dict
        pol = self.pol_dict
        for i in dictl:
            number_of_l = number_of_l + len(dictl[i])

        atomic_paobasis_card = str(self.name) + " " + str(number_of_l) + "\n"
        for i in dictl:
            for j in dictl[i]:
                if i in pol:
                    if j in pol[i]:
                        atomic_paobasis_card += "  n={}  {}  {}  P {} \n".format(i, j, len(dictl[i][j]), len(pol[i][j]))
                        listi = [dictl[i][j][l] * ang_to_bohr for l in dictl[i][j]]
                        atomic_paobasis_card += '\t'.join([f' {val}' for val in listi]) + "\n"
                    else:
                        atomic_paobasis_card += "  n={}  {}  {} \n".format(i, j, len(dictl[i][j]))
                        #print("  n={}  {}  {}".format(i, j, len(dictl[i][j])))
                        listi = [dictl[i][j][l] * ang_to_bohr for l in dictl[i][j]]
                        atomic_paobasis_card += '\t'.join([f' {val}' for val in listi]) + "\n"
                else:
                    atomic_paobasis_card += "  n={}  {}  {} \n".format(i, j, len(dictl[i][j]))
                    listi = [dictl[i][j][l] * ang_to_bohr for l in dictl[i][j]]
                    atomic_paobasis_card += '\t'.join([f' {val}' for val in listi]) + "\n"

        return atomic_paobasis_card[:-1]
