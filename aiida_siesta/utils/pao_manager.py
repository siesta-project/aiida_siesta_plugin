"""
The PAO Manager
"""

from aiida.plugins import DataFactory


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
            {n1 : { l1 : { z1 : r1, z2 : r2, ...}, l2 : ...} ...}
        where n* l* z* are integers carring the value of the n,l shells and z (please note we do not
        keep track of quantum number m as in the PAO construction is not possible to specify different
        feature for m shells with same l) and r* are the maximum radius for the n,l,z orbital.
        Maximum radius for the siesta orbitals means the radius after which the radial part is zero.
        """
        self.name = None
        self._gen_dict = None
        self._pol_dict = None

    def _validate_attrs(self):
        """
        Checks that the attributes are set, used in any function that
        is meant to return results.
        Plase note that we also require self._pol_dict to be not None, in case
        no polarized orbitals are present, self._pol_dict must be set to {}.
        """
        if self.name is None or self._gen_dict is None or self._pol_dict is None:
            raise RuntimeError("You need to set first a PAO block")

    def set_from_ion(self, ion_data_instance):
        """
        Sets the basic attributes of the class from an IonData.
        It goes through orbitals and extracts the two fundamental dictionaries of the class:
        `_gen_dict` and `_pol_dict`.
        """

        if not isinstance(ion_data_instance, DataFactory("siesta.ion")):
            raise ValueError(f"{ion_data_instance} is not an IonData instance")

        self.name = ion_data_instance.name

        gen_dict = {}
        pol_dict = {}
        for orbital in ion_data_instance.get_orbitals():
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

        self._gen_dict = gen_dict
        self._pol_dict = pol_dict

    def change_all_radius(self, percentage):
        """
        increment or decrement of all orbitas of a percentage.
        :param: percentage: positive (for inscreasing) or negative
                (for decrising) float representing the percentage
                of change of the radius
        Please nota that all radius are changed!
        """
        self._validate_attrs()

        dictl = self._gen_dict.copy()
        for i in dictl:
            for j in dictl[i]:
                for l in dictl[i][j]:
                    dictl[i][j][l] = dictl[i][j][l] + percentage / 100 * dictl[i][j][l]
        self._gen_dict = dictl

        if self._pol_dict:
            pol = self._pol_dict.copy()
            for i in pol:
                for j in pol[i]:
                    for l in pol[i][j]:
                        pol[i][j][l] = pol[i][j][l] + percentage / 100 * pol[i][j][l]
            self._pol_dict = pol

    def add_polarization(self, n, l):  #pylint: disable=invalid-name
        """
        Add polarization to the orbital with quantum numbers n, l.
        :param: n: principal quantum number of the orbital to be polarized
        :param: l: angular quantum number of the orbital to be polarized
        """

        self._validate_attrs()

        try:
            self._gen_dict[n][l]
        except KeyError:
            raise ValueError("no orbital with n = {0} and l = {1} is present in the basis".format(n, l))

        if l in self._pol_dict[n]:
            num_z = len(self._pol_dict[n][l])
            self._pol_dict[n][l][num_z + 1] = 0.000
        else:
            self._pol_dict[n][l] = {"1": self._gen_dict[n][l][1]}

    def get_pao_block(self):
        """
        From the info of the `_gen_dict`, `_pol_dict`, creates the PAO block.
        return: a string card containing the block.
        """

        self._validate_attrs()

        ang_to_bohr = 1.8897161646321
        number_of_l = 0
        dictl = self._gen_dict
        pol = self._pol_dict
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

    def pao_size(self):
        """
        Returns a string sumarazing the size of the basis in this class. Following the siesta convencion
        (SZ, DZ, DZP ...)
        Please note that the string is composed  checking the maximum Z registred by orbitals, for both the polarized
        and unpolarixed orbitals. This means that the algorithm is not able to actually detect if the orbitals
        here are generated by a "PAO.BasisSize" in siesta or it is a manual PAO block. Take this method carefully.
        """

        self._validate_attrs()

        max_z = 0
        for i in self._gen_dict:
            if len(self._gen_dict[i]) > max_z:
                max_z = len(self._gen_dict[i])

        if self._pol_dict:
            max_pz = 0
            for i in self._pol_dict:
                if len(self._pol_dict[i]) > max_pz:
                    max_pz = len(self._pol_dict[i])

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
