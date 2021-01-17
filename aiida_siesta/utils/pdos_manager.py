"""
The Pdos Manager
"""

from aiida.plugins import DataFactory

ANG_TO_BOHR = 1.8897161646321


class PdosManager:
    """
    Class to help plotting of Pdos. Help to select and sum contributions.
    Core of the implementation are is the internal
    attribute `_selected_orbitals` desribed in the __init__().
    """

    def __init__(self):
        """
        Initialize to None the fundamental attributes of the class:
        1) self._energies, the list of energies where the pdos was calculated.
        2) self._pdoses, the pdoses data
        3) self._original_orbitals, the list of all orbitals.
        The three attributes above are just taken form PdosData and never manipulated
        4) self._selected_orbitals, the list of orbitals currently selected by the manager.
           This is the key of the implementation, methods can add/remove orbitals from
           here, and then, only the orbitals in this attribute will be printed/plotted
           when requested.
           Each element of the list is a list contain only the fundamental quantities to allow the
           selection. in particular [orb, atom_tag, position, original_index]
           Where:
             orb is the sisl.AtomicOrbital object containing all the properties of the orbital (n,l,m,Z,P)
             atom_tag is the name of the atom the orbital belong to
             position is the cartesian coordinates of the atom
             original_index in the index in the orbital list of the PdosData
        """
        self._energies = None
        self._pdoses = None
        self._original_orbitals = None
        self._selected_orbitals = None

    def _validate_attrs(self):
        """
        Checks that the attributes are set, used in any methods of the class, except setters.
        """
        if self._energies is None or self._pdoses is None:
            raise RuntimeError("You need to set first a the pdos data with `set_from_pdos`")

        if self._selected_orbitals is None or self._original_orbitals is None:
            raise RuntimeError("You need to set first a the pdos data with `set_from_pdos`")

    def set_from_pdos(self, pdos_data_instance):
        """
        Sets the basic attributes of the class from an PdosData.
        :param pdos_data_instance: the PdosData instance from which to exctract the info.
        """

        if not isinstance(pdos_data_instance, DataFactory("siesta.pdos")):
            raise ValueError(f"{pdos_data_instance} is not an PdosData instance")

        self._energies = pdos_data_instance.get_energy_array()
        self._pdoses = pdos_data_instance.get_pdoses_array()
        self._original_orbitals = pdos_data_instance.get_orbitals_list()
        selected_orbitals = []
        for count, value in enumerate(pdos_data_instance.get_orbitals_list()):
            selected_orbitals.append(value + [count])
        self._selected_orbitals = selected_orbitals

    def reset_selection(self):
        """
        Reset `self._selected_orbitals` to its original value.
        """
        self._validate_attrs()

        selected_orbitals = []
        for count, value in enumerate(self._original_orbitals):
            selected_orbitals.append(value + [count])

        self._selected_orbitals = selected_orbitals

    def get_orbitals_list(self):
        """
        Return the `self._original_orbitals`
        """
        self._validate_attrs()

        return self._original_orbitals

    def get_selected_orbitals_list(self):
        """
        Return the `self._selected_orbitals`
        """
        self._validate_attrs()

        return self._selected_orbitals

    @staticmethod
    def _validate_after_selection(selected_orbitals, passed_vals, quantity_name):
        """
        The selection carried out in self.select() never intrinsecally fails, if the user
        wants to select an orbital that is not present, it won't be returned in the
        selected_orbitals dictionary.
        Here, we implement error if the call does not select any valid orbital,
        we implement warning if one of the request could not be satisfied.

        :param selected_orbitals: list, the list of selected orbitals
        :param passed_vals: list, the list passed by user, can be of any argument allowed by
             `self.select`. The argument name is in the next param.
        :param quantity_name: the quantity passed (species, n, l, orbitals, .. ) meaning one
             of the arguments of `self.select`.
        """

        import logging

        if not selected_orbitals:
            raise RuntimeError(
                "The selection process would resul in zero orbitals selected. Aborted. "
                "Use the method `self.get_selected_orbitals_list` to check the available orbitals"
            )

        if quantity_name == "orbital_indexes":
            if len(selected_orbitals) != len(passed_vals):
                logging.warning(
                    "One or more orbitals could not be selected. "
                    "Use the method `get_selected_orbitals_list` to check the selected orbitals. "
                    "The `orbital_index` is the last entry in each element of the list."
                )
        else:
            list_of = []
            for orb in selected_orbitals:
                if quantity_name == "species":
                    if orb[1] not in list_of:
                        list_of.append(orb[1])
                elif quantity_name == "n":
                    if orb[0].n not in list_of:
                        list_of.append(orb[0].n)
                elif quantity_name == "l":
                    if orb[0].l not in list_of:
                        list_of.append(orb[0].l)
            if len(list_of) != len(passed_vals):
                logging.warning(
                    f" One or more {quantity_name} could not be selected. "
                    "Use the method `self.get_selected_orbitals_list` to check the selected orbitals."
                )

    def select(self, species=None, n=None, l=None, orbital_indexes=None):  # noqa
        """
        Main mathod of the class. Accepts arguments to select the orbitals
        whose pdos we are interested in. Modifies `self._selected_orbitals`.
        Note that no validation is performed on the passed arguments. This is because
        the implementation selects the orbitals among the available ones (in
        self._selected_orbitals) with quantity (species, n, l, ..) equal to the passed
        values. Therefore passing a wrong value (for instance integers for `species`)
        just results in not selecting orbitals.

        :param species: the species of the atom the orbital belongs to
        :param n: principal quantum number
        :param l: angular quantum number
        :param orbital_indexes: index of the orbital
        """
        self._validate_attrs()

        if orbital_indexes is not None:

            #If orbital_indexes is set, no other argument must be passed
            if any([species, n, l]):
                raise ValueError("Either only the orbital_indexes or any other argument combination can be used")

            #Transform to list in case not iterable.
            if not isinstance(orbital_indexes, (tuple, list)):
                orbital_indexes = [orbital_indexes]

            selected_orbitals = []
            for sel_orb in self._selected_orbitals:
                if sel_orb[3] in orbital_indexes:
                    selected_orbitals.append(sel_orb)

            #Errors and warnings if some orbitals could not be selected
            self._validate_after_selection(selected_orbitals, orbital_indexes, "orbital_indexes")

            self._selected_orbitals = selected_orbitals

        if species is not None:
            if not isinstance(species, (tuple, list)):
                species = [species]

            selected_orbitals = []
            for sel_orb in self._selected_orbitals:
                if sel_orb[1] in species:
                    selected_orbitals.append(sel_orb)

            self._validate_after_selection(selected_orbitals, species, "species")

            self._selected_orbitals = selected_orbitals

        if n is not None:
            if not isinstance(n, (tuple, list)):
                n = [n]

            selected_orbitals = []
            for sel_orb in self._selected_orbitals:
                if sel_orb[0].n in n:
                    selected_orbitals.append(sel_orb)

            self._validate_after_selection(selected_orbitals, n, "n")

            self._selected_orbitals = selected_orbitals

        if l is not None:
            if not isinstance(l, (tuple, list)):
                l = [l]  # noqa

            selected_orbitals = []
            for sel_orb in self._selected_orbitals:
                if sel_orb[0].l in l:
                    selected_orbitals.append(sel_orb)

            self._validate_after_selection(selected_orbitals, l, "l")

            self._selected_orbitals = selected_orbitals

    def plot(self, sum_orbitals=False, custom_label=None, plot_show=False):
        """
        Plot the pdos of the selected orbitals. Also performs the sum of
        all the selected orbitals if requested. It does not show the plot
        unless `plot_show=True`. This is done to allow the selection of other
        pdoses to plot in the same graph.

        :param sum_orbitals: bool to activate the sum of pdoses
        :param custom_label: string to be passed as label for the curve. Allowed
           only if sum is active, because otherwise the label is automatically constracted
        :param: plot_show: bool to actevate `plt.show`
        """
        import matplotlib.pyplot as plt
        import numpy as np
        import logging

        n_spin = len(self._pdoses)

        if sum_orbitals:
            for spin in range(n_spin):
                sum_array = np.zeros(len(self._energies))
                for orbs in self._selected_orbitals:
                    sum_array = sum_array + self._pdoses[spin][orbs[3]]
                if custom_label is not None:
                    label = custom_label + " spin" + str(spin)
                else:
                    label = "spin" + str(spin)
                plt.plot(self._energies, sum_array, label=label)
                plt.legend()

        else:
            if custom_label is not None:
                logging.warning("`custom_label` implemented only for `sum_orbitals = True`.")
            for spin in range(n_spin):
                for orbs in self._selected_orbitals:
                    label = orbs[0].name() + " spin" + str(spin) + " " + orbs[1] + "@" + str(
                        np.around(orbs[2], decimals=4)
                    )
                    plt.plot(self._energies, self._pdoses[spin][orbs[3]], label=label)
                    plt.legend()

        if plot_show:
            plt.show()
        else:
            logging.warning(
                "The `plot` function does not call automatically `plt.show`. This allows "
                "to create other curves (with another selection process) to plot in the same graph. "
                "At the end, `import matplotlib.pyplot as plt` and call `plt.show()` to see the plot."
            )
