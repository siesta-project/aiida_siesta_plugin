import copy
from aiida.common.exceptions import NotExistent
from aiida.plugins import DataFactory
from aiida.orm import OrbitalData, load_node
from aiida.tools.data.orbital import Orbital
from aiida.common.exceptions import ValidationError
from aiida.plugins.entry_point import get_entry_point_from_class
from sisl import AtomicOrbital, Atom


class SislAtomicOrbital(AtomicOrbital, Orbital):
    """
    This is a subclass of `aiida.tools.data.orbital.Orbital` (general class representing
    orbitals in aiida) that presents the fetures of sisl.AtomicOrbital.
    """

    def __init__(self, *args, **kwargs):
        """
        The super is sisl.AtomicOrbital, therefore the instantiation follows its rules, as explained here:
        https://github.com/zerothi/sisl/blob/e28cdfff68d444139d54a755596de2a0285b0fd6/sisl/orbital.py#L689
        In addition we call `_set_orbital_dict`, that stores the info of the orbital in an internal dictionary.
        This is required in order to use this class as an orbital of a OrbitalData (the entity storable in the
        database! In fact this class is not a Data!)
        """
        super().__init__(*args, **kwargs)

        self._orbital_dict = self._set_orbital_dict()

    def set_orbital_dict(self, init_dict):
        """
        Because we want the users to follow the sisl way of defining orbitals, we disallow the possibility
        of setting the orbital properties after instanciation.
        """
        raise NotImplementedError(
            "SislAtomicOrbital does not implement `set_orbital_dict`, the attributes defining the orbital must be "
            "passed directly as aguments during the class instanciation, following the sisl standards"
        )

    def _set_orbital_dict(self):
        """
        Function called in the __init__, applys some checks and returns the orbital fetures in dictionary
        """

        entry_point = get_entry_point_from_class(self.__class__.__module__, self.__class__.__name__)[1]
        if entry_point is None:
            raise ValidationError(
                'Unable to detect entry point for current class {}, maybe you did not register an entry point for it?'.
                format(self.__class__)
            )

        if self.orb._r_inp_array is None or self.orb._f_inp_array is None:
            raise ValidationError(
                'SislAtomicOrbital do not support the passing of the radial function explicitly. '
                'A tuple (r,f) must be passed'
            )

        validated_dict = {}
        validated_dict['_orbital_type'] = entry_point.name
        validated_dict["n"] = self.n
        validated_dict["l"] = self.l
        validated_dict["m"] = self.m
        validated_dict["Z"] = self.Z
        validated_dict["P"] = self.P
        validated_dict["R"] = self.R
        validated_dict["q0"] = self.q0
        validated_dict['spherical'] = (self.orb._r_inp_array, self.orb._f_inp_array)

        return validated_dict


class BasisAtomicElement(OrbitalData):
    """
    Class with the scope to host the orbitals (instances of SislAtomicOrbital) asociated to an atomic species
    and implement several methods to help the management of basis.
    It makes use of the logic of OrbitalData, that stores in a serialized object (storable) the info of each orbital
    and knnows the logic for reconstructing SislAtomicOrbital instances.
    """

    def __init__(self, *args, **kwargs):
        """
        The super is necessary to import all the features of a Data object. Moreover initializes
        to None the defining attributes of the class.
        """

        super().__init__(*args, **kwargs)

        self.set_attribute("name", None)
        self.set_attribute("atomic_number", None)
        self.set_attribute("mass", None)
        self.set_attribute("pseudo_pk", None)
        self._orbital_list = []

    @property
    def name(self):
        """
        Name should correspond to the name of the atom in the structure (not the symbol)
        """
        return copy.deepcopy(self.get_attribute('name'))

    @name.setter
    def name(self, value):
        self.set_name(value)

    def set_name(self, value):
        from aiida.common.exceptions import ModificationNotAllowed
        if self.is_stored:
            raise ModificationNotAllowed('The BasisAtomicElement object cannot be modified, it has already been stored')
        self.set_attribute('name', value)

    @property
    def atomic_number(self):
        """
        The atomic number identifies the element (symbol in structure)
        """
        return copy.deepcopy(self.get_attribute('atomic_number'))

    @atomic_number.setter
    def atomic_number(self, value):
        self.set_atomic_number(value)

    def set_atomic_number(self, value):
        from aiida.common.exceptions import ModificationNotAllowed
        if self.is_stored:
            raise ModificationNotAllowed('The BasisAtomicElement object cannot be modified, it has already been stored')
        self.set_attribute('atomic_number', value)

    @property
    def mass(self):
        """
        Mass allows to distinguish isotopes.
        """
        return copy.deepcopy(self.get_attribute('mass'))

    @mass.setter
    def mass(self, value):
        self.set_mass(value)

    def set_mass(self, value):
        from aiida.common.exceptions import ModificationNotAllowed
        if self.is_stored:
            raise ModificationNotAllowed('The BasisAtomicElement object cannot be modified, it has already been stored')
        self.set_attribute('mass', value)

    @property
    def pseudo_pk(self):
        """
        For the basis (and orbitals)  generation, siesta starts from the pseudopotential. We force here to store
        the pseudo pk that ia associated to the orbitals in this class.
        """
        return copy.deepcopy(self.get_attribute('pseudo_pk'))

    @pseudo_pk.setter
    def pseudo_pk(self, value):
        self.set_pseudo_pk(value)

    def set_pseudo_pk(self, value):
        from aiida.common.exceptions import ModificationNotAllowed
        if self.is_stored:
            raise ModificationNotAllowed('The BasisAtomicElement object cannot be modified, it has already been stored')
        if value is not None:
            try:
                pse_node = load_node(value)
            except NotExistent:
                raise ValueError("pseudo_pk must be a valid pk number.")
            if not isinstance(pse_node,
                              DataFactory("siesta.psf")) and not isinstance(pse_node, DataFactory("siesta.psml")):
                raise ValueError("pseudo_pk must associated to a pseudo data node.")
        self.set_attribute('pseudo_pk', value)

    def set_atom_info(self, name, atomic_number, mass, pseudo_pk=None):
        """
        Function to collectively define the defining attributes of the class.
        """
        self.set_name(name)
        self.set_mass(mass)
        self.set_atomic_number(atomic_number)
        self.set_pseudo_pk(pseudo_pk)

    def set_orbitals(self, orbitals):
        """
        After performing some checks, it makes use of the `set_orbitals` of OrbitalData to
        set the orbitals. This function sets the attributes defining each orbital in a dictionary.
        Then this dictionary is put in a list that is associated to the attribute 'orbital_dicts'
        of the present class.
        """

        attr = self.attributes
        if attr["name"] is None or attr["atomic_number"] is None or attr["mass"] is None:
            raise ValueError("need to set first the atom_info with `set_atom_info`")

        if attr["pseudo_pk"] is None:
            self.logger.warning("It is always better to define also a pseudo_pk")

        if not isinstance(orbitals, list):
            orbitals = [orbitals]
        for orb in orbitals:
            if not isinstance(orb, SislAtomicOrbital):
                raise ValueError("only SislAtomicOrbital are accepted as orbitals")

        super().set_orbitals(orbitals)

    def set_from_sisl_atoms(self, sislatom):
        """
        Sets the Orbital data from a sisl.Atom instance. This is the class returned
        from sisl.read_basis().
        """

        if not isinstance(sislatom, Atom):
            raise ValueError("only sisl.AtomicOrbital are accepted as orbitals")

        self.set_attribute("name", sislatom.symbol)
        self.set_attribute("mass", sislatom.mass)
        self.set_attribute("atomic_number", sislatom.Z)
        listorb = []
        for i in sislatom.orbitals:
            atorb = SislAtomicOrbital(i.name(), (i.orb._r_inp_array, i.orb._r_inp_array), q0=i.q0)
            listorb.append(atorb)
        self.set_orbitals(listorb)

    @property
    def number_of_orbitals(self):
        return len(self.attributes["orbital_dicts"])

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
        for i in self.attributes["orbital_dicts"]:
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

        return gen_dict, pol_dict

    def get_pao_modifier(self):
        """
        Get the PaoModifiers, that makes available methods for the modification of
        a PAO block.
        """

        dictl, pol = self._get_pao_dict()

        return PaoModifier(self.name, dictl, pol)

    def get_pao_basis(self):
        """
        Returns the PAO block correspondent to the orbitals stored in this class.
        """

        dictl, pol = self._get_pao_dict()

        return PaoModifier(self.name, dictl, pol).get_pao_basis()

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

        return string


class PaoModifier:
    """
    class to help modifications to PAO basis block. Requires the gen_dict and pol_dict
    as defined in `_get_pao_dict` of the class BasisAtomicElement and also the name.
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

    def get_pao_basis(self):
        number_of_l = 0
        dictl = self.gen_dict
        pol = self.pol_dict
        for i in dictl:
            number_of_l = number_of_l + len(dictl[i])

        print(self.name, number_of_l)
        for i in dictl:
            for j in dictl[i]:
                if i in pol:
                    if j in pol[i]:
                        print("  n={}  {}  {}  P {}".format(i, j, len(dictl[i][j]), len(pol[i][j])))
                        listi = [dictl[i][j][l] for l in dictl[i][j]]
                        print(listi)
                    else:
                        print("  n={}  {}  {}".format(i, j, len(dictl[i][j])))
                        listi = [dictl[i][j][l] for l in dictl[i][j]]
                        print(listi)
                else:
                    print("  n={}  {}  {}".format(i, j, len(dictl[i][j])))
                    listi = [dictl[i][j][l] for l in dictl[i][j]]
                    print(listi)

    #for index in range(self.number_of_orbitals):
    #    i = self._orbital_list[index]
    #    for jindex in range(index+1,self.number_of_orbitals):
    #        j = self._orbital_list[jindex]
    #        print(i.orb.R, j.orb.R)
    #        if i.n == j.n and i.orb.l == j.orb.l and i.Z == j.Z and i.P == j.P:
    #            if (i.orb.R != j.orb.R):
    #                print("error")

    #nlist = []
    #for i in self._orbital_list:
    #    if i.n not in nlist:
    #        nlist.append(i.n)
    #        llist["{}".format(i.n)] = [i.l]

    #atomic_paobasis_card = "{}{}".format(self.get_attribute("name"),nlshell)


#class BasisData(Data):
#    """
#    A collection of BasisAtomicElement
#    """
