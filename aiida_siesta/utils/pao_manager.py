"""
The PAO Manager
"""

from aiida.plugins import DataFactory

ANG_TO_BOHR = 1.8897161646321


class PaoManager:
    """
    Class to help modifications of PAO basis block. Also translates orbitals info
    contained in ion files into a PAO block. Core of the implementation are the internal
    dictionaries `_gen_dict` and `_pol_dict`, desribed in the __iniy__().
    For the moment can only treat one single site at the time.
    """

    def __init__(self):
        """
        Initialize to None the fundamental attributes of the class:
        1) self.name, the name of the site associated to this PAO.
        2) self._gen_dict, containing information about the non polarized orbitals
           (the one with attr "P" False) of the PAO block.
        3) self._pol_dict, containing information about the polarized orbitals (the one
           with attr "P" True)
        The structure of both dicts is
            {n1 : { l1 : { Z1 : r1, Z2 : r2, ...}, l2 : ...} ...}
        where n* l* Z* are integers carring the value of the n,l shells and Z (please note we do not
        keep track of quantum number m as in the PAO construction is not possible to specify different
        feature for m shells with same l) and r* are the maximum radius for the n,l,z orbital.
        Maximum radius for the siesta orbitals means the radius after which the radial part is zero.
        Also note that the radii are parsed in Å (from IonData) and treated internally as Å like
        any other aiida internals, however the Pao Block requires the radii in Bohr, therefore
        extra care is necessary for methods that return the block and that accept explicit radii.
        Also three other dict are initialized.
        4) self._gen_occu, containing the occupations of the orbitals self._gen_dict
        5) self._pol_occu, containing the occupations of the orbitals self._pol_dict
        6) self._conf_dict, containing info on soft confinements and charge confinements.
        These dicts are set to empty (not None) because, for the moment, they are not considered
        mandatory to specify, so they are not checked as a requirement to have set.
        They are however set by `set_from_ion` and can be used for some internal reasons.
        """
        self.name = None
        self._gen_dict = None
        self._pol_dict = None
        self._gen_occu = {}
        self._pol_occu = {}
        self._conf_dict = {}

    def _validate_attrs(self, raise_if_empty=False):
        """
        Checks that the attributes are set, used in any methods of the class, except setters.
        Note that "not-set attributes" carry value "None". Empty dictionaty are instead
        valid values for the attributes.
        The self._pol_dict = {} signals basis with no polarized orbital. The self._gen_dict = {}
        might appen during the construction, but it is not allowed when the results
        should be returned. Therefore the `raise_if_empty` argument signals to check
        also if self._gen_dict is empty and it is called in the getters methods.
        """
        if self.name is None or self._gen_dict is None or self._pol_dict is None:
            raise RuntimeError("You need to set first a PAO block")

        if raise_if_empty:
            if not self._gen_dict:
                raise RuntimeError("No orbitals set, nothing to return")

    @staticmethod
    def _validate_nlz(n, l, Z):  # noqa
        """
        Checks that the passed n, l, Z conform with what it expected.

        :param n: int, principal quantum number of the orbital to be polarized
        :param l: int, angular quantum number of the orbital to be polarized
        :param Z: int, the Z orbital to be changed.
        """
        for i in [n, l, Z]:
            if not isinstance(i, int):
                raise ValueError("n, l, Z must be integers")

        if n < 0:
            raise ValueError("n must be bigger or equal to zero")

        if l < 0:
            raise ValueError("l must be bigger or equal to zero")

        if Z <= 0:
            raise ValueError("Z must be bigger then zero")

    def set_from_ion(self, ion_data_instance):  # noqa
        """
        Sets the basic attributes of the class from an IonData.
        It goes through orbitals and extracts the two fundamental dictionaries of the class:
        `_gen_dict` and `_pol_dict`.

        :param ion_data_instance: the IonData instance from which to exctract the info.
        """

        if not isinstance(ion_data_instance, DataFactory("siesta.ion")):
            raise ValueError(f"{ion_data_instance} is not an IonData instance")

        self.name = ion_data_instance.name

        gen_dict = {}
        pol_dict = {}
        gen_occu = {}
        pol_occu = {}
        for orbital in ion_data_instance.get_orbitals():
            i = orbital.attributes
            if not i["P"]:
                if i["n"] not in gen_dict:
                    gen_dict[i["n"]] = {i["l"]: {i["Z"]: i["R"]}}
                    gen_occu[i["n"]] = {i["l"]: {i["Z"]: i["q0"]}}
                else:
                    if i["l"] not in gen_dict[i["n"]]:
                        gen_dict[i["n"]][i["l"]] = {i["Z"]: i["R"]}
                        gen_occu[i["n"]][i["l"]] = {i["Z"]: i["q0"]}
                    else:
                        if i["Z"] not in gen_dict[i["n"]][i["l"]]:
                            gen_dict[i["n"]][i["l"]][i["Z"]] = i["R"]
                            gen_occu[i["n"]][i["l"]][i["Z"]] = i["q0"]
                        else:
                            gen_occu[i["n"]][i["l"]][i["Z"]] += i["q0"]
            else:
                if i["n"] not in pol_dict:
                    pol_dict[i["n"]] = {i["l"] - 1: {i["Z"]: i["R"]}}
                    pol_occu[i["n"]] = {i["l"] - 1: {i["Z"]: i["q0"]}}
                else:
                    if i["l"] - 1 not in pol_dict[i["n"]]:
                        pol_dict[i["n"]][i["l"] - 1] = {i["Z"]: i["R"]}
                        pol_occu[i["n"]][i["l"] - 1] = {i["Z"]: i["q0"]}
                    else:
                        if i["Z"] not in pol_dict[i["n"]][i["l"] - 1]:
                            pol_dict[i["n"]][i["l"] - 1][i["Z"]] = i["R"]
                            pol_occu[i["n"]][i["l"] - 1][i["Z"]] = i["q0"]
                        else:
                            pol_occu[i["n"]][i["l"] - 1][i["Z"]] += i["q0"]

        #The polarization of 2p is 3d (2d does not exist).
        #Same for polarization of 1s e 3d.
        for num in [2, 3, 4]:
            if num in pol_dict and num not in gen_dict:
                pol_dict[num - 1] = pol_dict[num]
                pol_dict.pop(num)
                pol_occu[num - 1] = pol_occu[num]
                pol_dict.pop(num)

        map_l = {"s": 0, "p": 1, "d": 2, "f": 3}
        confinement_dict = {}
        if ion_data_instance.get_info_charge_confinement():
            confinement_dict["Q"] = {}
            for k, v in ion_data_instance.get_info_charge_confinement().items():
                confinement_dict["Q"][int(k[0])] = {map_l[k[1]]: v}
        if ion_data_instance.get_info_soft_confinement():
            confinement_dict["E"] = {}
            for k, v in ion_data_instance.get_info_soft_confinement().items():
                confinement_dict["E"][int(k[0])] = {map_l[k[1]]: v}

        self._conf_dict = confinement_dict  #It might be empty
        self._gen_dict = gen_dict
        self._pol_dict = pol_dict
        self._gen_occu = gen_occu
        self._pol_occu = pol_occu

    def change_all_radius(self, percentage):
        """
        Increment or decrement the radius of all orbitas of a percentage.

        :param percentage: positive (for inscreasing) or negative
            (for decrising) float representing the percentage
            of change of the radius.

        All radii are changed, also the polarized one.
        """
        self._validate_attrs()

        dictl = self._gen_dict.copy()
        for i in dictl:
            for j in dictl[i]:
                for l in dictl[i][j]:  # noqa
                    dictl[i][j][l] = dictl[i][j][l] + percentage / 100 * dictl[i][j][l]
        self._gen_dict = dictl

        #Note that for the sake of a new pao_block, the radius of polarized
        #orbitals does not matter. However we change it for consistency
        if self._pol_dict:
            pol = self._pol_dict.copy()
            for i in pol:
                for j in pol[i]:
                    for l in pol[i][j]:  # noqa
                        pol[i][j][l] = pol[i][j][l] + percentage / 100 * pol[i][j][l]
            self._pol_dict = pol

    def reset_radius(self, radius_units, new_radius, n, l, Z=1):  # noqa
        """
        Reset the radius of an orbital with n, l, Z quantum numbers.

        :param radius_units: either Bohr or Ang
        :param new_radius: new radius that will be set, in Bohr or Ang
        :param n: int, principal quantum number of the orbital to be polarized
        :param l: int, angular quantum number of the orbital to be polarized
        :param Z: int, the Z orbital to be changed.

        For consistency, if a Z=1 orbital is changed, also the corresponding
        polarized orbital is changed if present and the following Zs of the
        polarized orbital (if presents) are set to zero.
        Note: this last action is irrelevant for the the creation of Pao Blocks
        as the radii of polarized orbitals can not be set manually in the Pao block.
        """

        self._validate_attrs()
        self._validate_nlz(n, l, Z)

        if radius_units not in ["Bohr", "Ang"]:
            raise ValueError("`radius_units` only accepts 'Bohr' or 'Ang'")

        if radius_units == "Bohr":
            new_radius = new_radius / ANG_TO_BOHR

        if n not in self._gen_dict or l not in self._gen_dict[n] or Z not in self._gen_dict[n][l]:
            raise ValueError(f"No orbital defined with n={n}, l={l}, Z={Z}")

        self._gen_dict[n][l][Z] = new_radius

        if Z == 1:
            try:
                self._pol_dict[n][l][Z] = new_radius
                if len(self._pol_dict[n][l]) > 1:
                    for other_z in range(len(self._pol_dict[n][l]) - 1):
                        self._pol_dict[n][l][other_z + 2] = 0.000
            except KeyError:
                pass

    def add_polarization(self, n, l):  # noqa
        """
        Add polarization to the orbital with quantum numbers n, l.

        :param: n: principal quantum number of the orbital to be polarized
        :param: l: angular quantum number of the orbital to be polarized

        The polarized orbital is added in the internal `self._pol_dict` with
        the same radius of the corresponding orbital in the `self._gen_dict`.
        If more z are present, the radius is the one of the first zeta.
        Note that for the sake of a new pao_block, the radius of polarized
        orbitals does not matter. However we introduce it for consistency
        """

        self._validate_attrs()
        self._validate_nlz(n, l, 1)

        try:
            self._gen_dict[n][l]
        except KeyError:
            raise ValueError("no orbital with n = {0} and l = {1} is present in the basis".format(n, l))

        if n not in self._pol_dict:
            self._pol_dict[n] = {}

        if l in self._pol_dict[n]:
            num_z = len(self._pol_dict[n][l])
            self._pol_dict[n][l][num_z + 1] = 0.000
        else:
            self._pol_dict[n][l] = {1: self._gen_dict[n][l][1]}

    def remove_polarization(self, n, l):  # noqa
        """
        Add polarization to the orbital with quantum numbers n, l.

        :param: n: principal quantum number of the orbital to be polarized
        :param: l: angular quantum number of the orbital to be polarized
        """

        self._validate_attrs()
        self._validate_nlz(n, l, 1)

        try:
            self._pol_dict[n][l]
        except KeyError:
            raise ValueError(f"no polarized orbital with n={n} and l={l} is present in the basis")

        has_occu = self._pol_occu.get(n, {}).get(l)

        num_z = len(self._pol_dict[n][l])
        if num_z > 1:
            self._pol_dict[n][l].pop(num_z)
            if has_occu is not None:
                self._pol_occu[n][l].pop(num_z)
        else:
            self._pol_dict[n].pop(l)
            if has_occu is not None:
                self._pol_occu[n].pop(l)

        #remove eampty self._pol_dict[n]
        if not self._pol_dict[n]:
            self._pol_dict.pop(n)
            if has_occu is not None:
                self._pol_occu.pop(n)

    def add_orbital(self, radius_units, radius, n, l, Z=1):  # noqa
        """
        Add an orbital with n, l, Z quantum numbers.

        :param radius_units: either Bohr or Ang
        :param radius: new radius that will be set, in Bohr or Ang
        :param n: int, principal quantum number of the orbital to be polarized
        :param l: int, angular quantum number of the orbital to be polarized
        :param Z: int, the Z orbital to be changed.

        Note that the addition of Zs is sequential, meaning it is not possible
        to add Z=m if Z=m-1 is not present.
        """

        self._validate_attrs()
        self._validate_nlz(n, l, Z)

        if radius_units not in ["Bohr", "Ang"]:
            raise ValueError("`radius_units` only accepts 'Bohr' or 'Ang'")

        if radius_units == "Bohr":
            radius = radius / ANG_TO_BOHR

        if n in self._gen_dict and l in self._gen_dict[n] and Z in self._gen_dict[n][l]:
            raise ValueError(f"Orbital with n={n}, l={l}, Z={Z} is already present. Can not be added.")

        if Z != 1:
            try:
                self._gen_dict[n][l][Z - 1]
            except KeyError:
                raise ValueError(f"Orbital with n={n}, l={l}, Z={Z-1} not present. Can not add further Zs")

        try:
            self._gen_dict[n][l][Z] = radius
        except KeyError:
            try:
                self._gen_dict[n][l] = {Z: radius}
            except KeyError:
                self._gen_dict[n] = {l: {Z: radius}}

    def remove_orbital(self, n, l, Z):  # noqa
        """
        Add an orbital with n, l, Z quantum numbers.

        :param radius_units: either Bohr or Ang
        :param radius: new radius that will be set, in Bohr or Ang
        :param n: int, principal quantum number of the orbital to be polarized
        :param l: int, angular quantum number of the orbital to be polarized
        :param Z: int, the Z orbital to be changed.

        If no orbital with n, l is left, also the corresponding polarized
        orbital, if present, is deleted.
        """

        self._validate_attrs()
        self._validate_nlz(n, l, Z)

        if n not in self._gen_dict or l not in self._gen_dict[n] or Z not in self._gen_dict[n][l]:
            raise ValueError(f"Orbital with n={n}, l={l}, Z={Z} is not present. Can not be removed.")

        if Z + 1 in self._gen_dict[n][l]:
            raise ValueError(f"Orbital with n={n}, l={l}, Z={Z+1} found. This needs to be removed first")

        self._gen_dict[n][l].pop(Z)

        has_occu = self._gen_occu.get(n, {}).get(l, {}).get(Z)
        if has_occu is not None:
            self._gen_occu[n][l].pop(Z)

        #remove eampty dicts, including polarized and confinements
        if not self._gen_dict[n][l]:
            self._gen_dict[n].pop(l)
            if has_occu is not None:
                self._gen_occu[n].pop(l)
            try:
                self._pol_dict[n].pop(l)
            except KeyError:
                pass
            try:
                self._pol_occu[n].pop(l)
            except KeyError:
                pass
            have_q_conf = self._conf_dict.get("Q", {}).get(n, {}).get(l)
            if have_q_conf is not None:
                self._conf_dict["Q"][n].pop(l)
                if not self._conf_dict["Q"][n]:
                    self._conf_dict["Q"].pop(n)
                if not self._conf_dict["Q"]:
                    self._conf_dict.pop("Q")
            have_e_conf = self._conf_dict.get("E", {}).get(n, {}).get(l)
            if have_e_conf is not None:
                self._conf_dict["E"][n].pop(l)
                if not self._conf_dict["E"][n]:
                    self._conf_dict["E"].pop(n)
                if not self._conf_dict["E"]:
                    self._conf_dict.pop("E")
        if not self._gen_dict[n]:
            self._gen_dict.pop(n)
            if has_occu is not None:
                self._gen_occu.pop(n)
        try:
            if not self._pol_dict[n]:
                self._pol_dict.pop(n)
        except KeyError:
            pass
        try:
            if not self._pol_occu[n]:
                self._pol_occu.pop(n)
        except KeyError:
            pass

    def get_pao_block(self):  # noqa: MC0001  - is mccabe too complex funct -
        """
        From the info of the `_gen_dict`, `_pol_dict`, creates the PAO block.

        return: a string card containing the block.

        Conversion into Bohr is performed.
        Generally the radii values are floats. However we allow a hack in the
        BasisOptimizer to use stings for values. For this reason the conversion
        into Bohr is in an if statement.
        """
        import copy

        self._validate_attrs(raise_if_empty=True)

        number_of_l = 0
        dictl = copy.deepcopy(self._gen_dict)
        pol = copy.deepcopy(self._pol_dict)
        for i in dictl:
            number_of_l = number_of_l + len(dictl[i])

        #Conversion in bohr. Polarization is not necessary, but do for consistency
        for i in dictl:
            for j in dictl[i]:
                for l in dictl[i][j]:  # noqa
                    if isinstance(dictl[i][j][l], (float, int)):
                        dictl[i][j][l] = round(dictl[i][j][l] * ANG_TO_BOHR, 6)
        for i in pol:
            for j in pol[i]:
                for l in pol[i][j]:  # noqa
                    if isinstance(pol[i][j][l], (float, int)):
                        pol[i][j][l] = round(pol[i][j][l] * ANG_TO_BOHR, 6)

        atomic_paobasis_card = str(self.name) + " " + str(number_of_l) + "\n"
        for in_n in [0, 1, 2, 3, 4, 5, 6, 7, 8]:
            for i in dictl:
                if i == in_n:
                    for j in dictl[i]:
                        atomic_paobasis_card += f"  n={i}  {j}  {len(dictl[i][j])}"
                        ij_pol = pol.get(i, {}).get(j)  #This works only because pol[i] is always dict, same below
                        ij_q = self._conf_dict.get("Q", {}).get(i, {}).get(j)
                        ij_e = self._conf_dict.get("E", {}).get(i, {}).get(j)
                        if ij_pol:
                            atomic_paobasis_card += f" P {len(pol[i][j])}"
                        if ij_q:
                            values_q = "".join([f'{val} ' for val in self._conf_dict['Q'][i][j]])
                            atomic_paobasis_card += f" Q {values_q}"
                        if ij_e:
                            values_e = "".join([f'{val} ' for val in self._conf_dict['E'][i][j]])
                            atomic_paobasis_card += f" E {values_e}"
                        atomic_paobasis_card += "\n"
                        listi = [dictl[i][j][l] for l in dictl[i][j]]  # noqa
                        atomic_paobasis_card += '\t'.join([f' {val}' for val in listi]) + "\n"

        return atomic_paobasis_card[:-1]

    def pao_size(self):
        """
        Returns a string sumarazing the size of the basis in this class. Following the siesta convencion
        (SZ, DZ, DZP ...)
        Please note that the string is composed  checking the maximum Z registred by orbitals, for both the polarized
        and unpolarixed orbitals. This means that the algorithm is not able to actually detect if the orbitals
        here are generated by a "PAO.BasisSize" in siesta or it is a manual PAO block. Take this method carefully.
        """

        self._validate_attrs(raise_if_empty=True)

        max_z = 0
        for i in self._gen_dict:
            for j in self._gen_dict[i]:
                if len(self._gen_dict[i][j]) > max_z:
                    max_z = len(self._gen_dict[i][j])

        if self._pol_dict:
            max_pz = 0
            for i in self._pol_dict:
                for j in self._pol_dict[i]:
                    if len(self._pol_dict[i][j]) > max_pz:
                        max_pz = len(self._pol_dict[i][j])

        dictma = ["S", "D", "T"]
        dictpo = ["", "D", "T"]
        if self._pol_dict:
            string = "{}Z{}P".format(dictma[max_z - 1], dictpo[max_pz - 1])
        else:
            string = "{}Z".format(dictma[max_z - 1])

        import logging

        logging.warning(
            " Carefull. This information is an approximation. It is extracted looking "
            "at the maximum Z for polarized and unpolarized orbitals. Not all the siesta basis "
            "choises can be condensed into a simple string. Check the `get_pao_block` for "
            "more detailed information"
        )

        return string
