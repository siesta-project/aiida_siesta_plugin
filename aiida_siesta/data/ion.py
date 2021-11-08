"""
This module manages the .ion.xml files in the local repository.
"""
import io
import pathlib
from aiida.common.files import md5_from_filelike
from aiida.orm.nodes import SinglefileData
from aiida.common.exceptions import StoringNotAllowed
from aiida_siesta.utils.pao_manager import PaoManager


def xml_element_to_string(element, tail=True):
    string = "<" + element.tag + ">" + element.text + "</" + element.tag + ">"
    if tail:
        string = string + element.tail
    return string


def parse_ion(stream):
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

    el_tr = ElementTree(None, stream)
    root = el_tr.getroot()

    if root.find('symbol') is None or root.find('z') is None or root.find('label') is None:
        raise ParsingError(f"Currupted ion file {stream.name}: element symbol or atomic number missing")

    parsed_data["element"] = str(root.find('symbol').text.strip())
    if parsed_data["element"] not in _valid_symbols:
        raise ParsingError(f"Unknown element symbol {parsed_data['element']} in file {stream.name}")

    parsed_data["name"] = str(root.find('label').text.strip())
    parsed_data["atomic_number"] = int(root.find('z').text)
    parsed_data["mass"] = float(root.find('mass').text)

    return parsed_data


class IonData(SinglefileData):
    """
    Handler for ion files
    """

    @staticmethod
    def _prepare_source(source):
        """
        if the ``source`` is a valid file on disk, its content is read and returned as a stream of bytes.
        """
        if isinstance(source, (str, pathlib.Path)):
            save_name = source
            with open(source, 'rb') as handle:
                source = io.BytesIO(handle.read())
                source.name = save_name

        return source

    def set_file(self, source, filename=None):  #pylint: disable=arguments-differ
        """
        This is called in the __init__ of SingleFileData. It supports both absolute path and file streams.
        It is convenient to convert possible absolute paths in the corresponding file streams so we then
        support for other methods called here only the file streams.
        Please note that this approach have problems if we create subclasses with `IonData` as a parent.
        This is because the call to super does not return anything.
        Therefore we can not have source = super().set_file(source, filename, **kwargs)
        """
        # Check we have a valid input and set the file and filename as attributes
        super().set_file(source, filename)
        # Transorm abs_paths in streams
        source = self._prepare_source(source)
        source.seek(0)
        # Set the md5 attribute
        self.set_attribute('md5', md5_from_filelike(source))
        source.seek(0)
        # Set other attributes extracted reading the source
        parsed_data = parse_ion(source)
        source.seek(0)
        self.set_attribute('element', parsed_data["element"])
        self.set_attribute('name', parsed_data["name"])
        self.set_attribute('atomic_number', parsed_data["atomic_number"])
        if parsed_data["mass"] is not None:
            self.set_attribute('mass', parsed_data["mass"])

    def store(self, **kwargs):  # pylint: disable=arguments-differ
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

        return super().store(**kwargs)

    def validate_others_atts(self, elem, name, atm_n):
        """
        Validate the given element, name, atomic_number are the one of the stored file.
        Unfortunately it requires to reparse the file.
        :param elem: the symbol of the element.
               name: the name assigned to the atom/site.
               atm_n: the atomic number of the atom/site.
        :raises ValueError: if the element symbol is invalid.
        """
        with self.open(mode='rb') as handle:
            parsed_data = parse_ion(handle)
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
    def get_or_create(cls, source, filename=None):
        """
        Pass the same parameter of the __init__; if a file with the same md5
        is found, that IonData is returned, otherwise a new IonFile instance is
        created.

        :param source: an absolute path file on disk or a filelike object.
        :param filename: optional explicit filename to give to the file stored in the repository.
                         Ignored if a file with the same md5 has been found.
        :return ion: the IonData object.
        """
        from aiida import orm

        if isinstance(source, (str, pathlib.Path)):
            if not pathlib.Path(source).is_file():
                raise TypeError(f'`source` should be a str or pathlib.Path of a filepath on disk, got: {source}')
            source = cls._prepare_source(source)
            source.seek(0)

        readable_bytes = (hasattr(source, 'read') and hasattr(source, 'mode') and 'b' in source.mode)
        bol = isinstance(source, io.BytesIO) or readable_bytes
        if not bol:
            raise TypeError(
                f'`source` should be a str or `pathlib.Path` of a filepath on disk or a stream of bytes, got: {source}'
            )

        query = orm.QueryBuilder()
        query.append(cls, subclassing=False, filters={'attributes.md5': md5_from_filelike(source)})

        existing = query.first()

        if existing:
            ion = existing[0]
        else:
            source.seek(0)
            ion = cls(source, filename)

        return ion

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

    def get_content_ascii_format(self):  #pylint: disable=too-many-statements
        """
        from the content, write the old format .ion file. Necessary since siesta only reads
        ion info in this format.
        """
        from xml.etree import ElementTree

        root = ElementTree.fromstring(self.get_content())

        #preliminary check on lj_projs, necessary due to a problem in siesta.
        #See "Add lj_projs and j support to ion xml files" commit to siesta in GitLab
        found_lj_proj = root.find("lj_projs")
        if found_lj_proj is not None:
            have_lj_proj = "T" in found_lj_proj.text or "t" in found_lj_proj.text
        else:
            ln_list = []
            for proj in root.find("kbs"):
                ln_list.append((int(proj.attrib["l"]), int(proj.attrib["n"])))
            if len(ln_list) == len(set(ln_list)):
                have_lj_proj = False
            else:
                have_lj_proj = True

        #start of the file content construction
        string = ""

        #Construct the preamble (basis spec and pseudo header)
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
        if have_lj_proj:
            string = string + root.find("lmax_projs").text + root.find("nprojs_nl").text + "T\n"
        else:
            string = string + root.find("lmax_projs").text + root.find("nprojs_nl").text + "#\n"

        #The Paos
        string = string + "# PAOs:__________________________\n"
        for orbital in root.find("paos"):
            string = string + orbital.attrib["l"] + orbital.attrib["n"] + orbital.attrib["z"] + orbital.attrib[
                "ispol"] + orbital.attrib["population"] + "\n"
            radfunc = orbital.find("radfunc")
            string = string + radfunc.find("npts").text + radfunc.find("delta").text + radfunc.find("cutoff").text
            string = string + radfunc.find("data").text

        #The KBs. Note that (in case of have_lj_proj) the j value is not read from the .ion.xml but calculated
        #on site. This is because the j value for each projector was added only in recent version of siesta.
        #The implementation assumes that j=l-1/2 is always the first listed, j=l+1/2 the second! Hope it is true!!
        string = string + "# KBs:__________________________\n"
        collect_ln = []
        for projector in root.find("kbs"):
            l_val = int(projector.attrib["l"])
            n_val = int(projector.attrib["n"])
            if have_lj_proj:
                if (l_val, n_val) in collect_ln:
                    j_val = "   " + str(l_val + 0.5) + "  "
                else:
                    j_val = "   " + str(abs(l_val - 0.5)) + "  "  #abs for the l=0 case
                string = string + " " + str(l_val) + j_val + str(n_val) + projector.attrib["ref_energy"] + "\n"
            else:
                string = string + " " + str(l_val) + "  " + str(n_val) + projector.attrib["ref_energy"] + "\n"
            collect_ln.append((l_val, n_val))
            radfunc = projector.find("radfunc")
            string = string + radfunc.find("npts").text + radfunc.find("delta").text + radfunc.find("cutoff").text
            string = string + radfunc.find("data").text

        #Other quantities
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

    def _analyze_basis_specs(self):
        """
        From the file content get the info on soft confinement and charge confinement.
        The detection is based on reading the `<basis_specs>` section of the .ion file,
        and it is pure string analysis. Therefore it is sensible to the change
        of the printed info in Siesta. The current implementation supports only Max Siesta versions.
        Monitor the basis_type.f90 file of Siesta to see if some changes
        are introduced.
        """
        #Extract first the `<basis_specs>` block from the content and put each word in a list
        content = self.get_content()
        str_start = content.index("<basis_specs>") + len("<basis_specs>")
        str_end = content.index("</basis_specs>", str_start)
        interest_content = content[str_start:str_end].split()

        #Loop over the list to get the needed info
        dict_q_and_e = {}
        for ind, val in enumerate(interest_content):
            if 'i=' in val:
                dict_key = interest_content[ind + 3].strip("(").strip(")")
                # Discart situation when no info on `vcte` are reported, meaning when
                # perturbative polarization orbital or empty shell. Note that this distintion
                # can be done only on MaX versions since the checked strings have been introduced
                # only recently. With previous version there was no way to distinguish these orbitals.
                if interest_content[ind + 4] != "(perturbative" and interest_content[ind + 4] != "(empty":
                    for ind2, val2 in enumerate(interest_content[ind:]):
                        saveind = ind2
                        if val2 == "vcte:":
                            break
                    dict_q_and_e[dict_key] = {
                        "Q": [
                            float(interest_content[ind + saveind + 5]),
                            float(interest_content[ind + saveind + 7]),
                            float(interest_content[ind + saveind + 9])
                        ],
                        "E": [float(interest_content[ind + saveind + 1]),
                              float(interest_content[ind + saveind + 3])]
                    }

        return dict_q_and_e

    def get_info_charge_confinement(self):
        """
        Return the orbitals with charge confinement
        """
        dict_q_and_e = self._analyze_basis_specs()
        collect_q = {}
        for basis_el_k, basis_el_v in dict_q_and_e.items():
            if basis_el_v["Q"][0] != 0.0:
                collect_q[basis_el_k] = basis_el_v["Q"]

        return collect_q

    def get_info_soft_confinement(self):
        """
        Return the orbitals with charge confinement
        """
        dict_q_and_e = self._analyze_basis_specs()
        collect_e = {}
        for basis_el_k, basis_el_v in dict_q_and_e.items():
            if basis_el_v["E"][0] != 0.0:
                collect_e[basis_el_k] = basis_el_v["E"]

        return collect_e
