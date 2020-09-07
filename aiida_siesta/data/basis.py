import copy
from aiida.orm import OrbitalData
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
        super().__init__(*args, **kwargs)

        self._orbital_dict = self._set_orbital_dict()

    def set_orbital_dict(self, init_dict):
        raise NotImplementedError(
            "SislAtomicOrbital does not implement `set_orbital_dict`, the attributes defining the orbital must be "
            "passed directly as aguments during the class instanciation, following the sisl standards"
        )

    def _set_orbital_dict(self):

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
        #validated_dict["_r_inp_array"] = self.orb._r_inp_array
        #validated_dict["_f_inp_array"] = self.orb._f_inp_array

        return validated_dict

    def get_orbital_dict(self):
        """
        returns the internal keys as a dictionary
        """
        return self._orbital_dict


class BasisAtomicElement(OrbitalData):
    """
    Class that contains the orbitas asociated to an atomic species.
    """

    def __init__(self, *args, **kwargs):

        super().__init__(*args, **kwargs)

        self.set_attribute("name", None)
        self.set_attribute("atomic_number", None)
        self.set_attribute("mass", None)
        self.set_attribute("pseudo_pk", None)
        self._orbital_list = []

    @property
    def name(self):
        """
        Returns the cell shape.

        :return: a 3x3 list of lists.
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
        Returns the cell shape.

        :return: a 3x3 list of lists.
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
        Returns the cell shape.

        :return: a 3x3 list of lists.
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

    def set_atom_info(self, name, atomic_number, mass):
        self.set_name(name)
        self.set_mass(mass)
        self.set_atomic_number(atomic_number)

    #def set_pseudo(self, pseudo):
    #    if not isinstance(pseudo,PsmlData) or isinstance(pseudo, PsfData):
    #        raise ValueError("pseudo needs to be either a PsfData or PsmlData")
    #    try:
    #        self.set_attribute("pseudo_pk", pseudo.pk)
    #    except AttributeError:
    #        raise ValueError("The pseudo must be stored")

    def set_orbitals(self, orbitals):

        attr = self.attributes
        if attr["name"] is None or attr["atomic_number"] is None or attr["mass"] is None:
            raise ValueError("need to set first the atom_info with `set_atom_info`")
    #    if attr["pseudo_pk"] is None:
    #        raise ValueError("need to set first the associated pseudo with `set_pseudo`")

    #if not isinstance(orbital, AtomicOrbital):
    #    raise ValueError("only sisl.AtomicOrbital are accepted as orbitals")

    #self._orbital_list.append(orbital)
        super().set_orbitals(orbitals)

    #def get_orbitals(self):
    #    return self._orbital_list

    @property
    def number_of_orbitals(self):
        return len(self.attributes["orbital_dicts"])

    def set_from_sisl_atoms(self, sislatom):
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

    def get_pao_dict(self):
        dictl = {}
        pol = {}
        for i in self.attributes["orbital_dicts"]:
            if not i["P"]:
                if i["n"] not in dictl:
                    dictl[i["n"]] = {i["l"]: {i["Z"]: i["R"]}}
                else:
                    if i["l"] not in dictl[i["n"]]:
                        dictl[i["n"]][i["l"]] = {i["Z"]: i["R"]}
                    else:
                        if i["Z"] not in dictl[i["n"]][i["l"]]:
                            dictl[i["n"]][i["l"]][i["Z"]] = i["R"]
            else:
                if i["n"] not in pol:
                    pol[i["n"]] = {i["l"] - 1: {i["Z"]: i["R"]}}
                else:
                    if i["l"] - 1 not in pol[i["n"]]:
                        pol[i["n"]][i["l"] - 1] = {i["Z"]: i["R"]}
                    else:
                        if i["Z"] not in pol[i["n"]][i["l"] - 1]:
                            pol[i["n"]][i["l"] - 1][i["Z"]] = i["R"]

        return dictl, pol

    def get_pao_modifier(self):

        dictl, pol = self.get_pao_dict()

        return PaoModifier(self.name, dictl, pol)


class PaoModifier:
    """
    class to help modifications to PAO
    """

    def __init__(self, name, dictl, pol=None):
        self.name = name
        self.dictl = dictl
        self.pol = pol

    def change_all_radius(self, percentage):
        """
        increment decrement of all orbitas of a percentage
        """

    #def add_polarization(self, n, l):
    #    """
    #    add polarization to the orbital with quantum numbers n, l
    #    """

    def get_pao_basis(self):
        number_of_l = 0
        dictl = self.dictl
        pol = self.pol
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
