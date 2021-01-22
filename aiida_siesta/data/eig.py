from aiida.orm import KpointsData
K2EV = 0.025851984449230724


class EigData(KpointsData):
    """
    Class to handle eig data
    """

    @property
    def e_fermi(self):
        """
        The fermi energy, in eV
        """
        return self.get_attribute('e_fermi', None)

    @e_fermi.setter
    def e_fermi(self, e_fermi):
        if not isinstance(e_fermi, (float, int)):
            raise ValueError("The e_fermi must be an integer or a float")
        self.set_attribute('e_fermi', e_fermi)

    def set_eigs(self, eigs):
        """
        Set the eigenvalues, must be in eV. The `set_array` checks wether the passed
        quantity is a np.array, here we also ensure the array is of the right type.
        In fact, arrays of (for instance) strings would lead to crashes in other methos
        of this class.
        """
        self.set_array('eigs', eigs)

        if eigs.dtype.kind not in ["b", "i", "u", "f", "c", "m", "M"]:
            raise TypeError("the `eigs` array must be a numerical type")

        try:
            eigs[0][0][0]
        except IndexError:
            raise IndexError("the eigenvalues must be a np.array of dimension n_spin x n_kp x n_energies")

    def get_eigs(self, scale_to_ef=False):
        """
        Get the eigenvalues
        """
        import numpy as np
        try:
            eig_array = np.array(self.get_array('eigs'))
        except KeyError:
            raise AttributeError('No stored eigs has been found. Use method `set_eigs`.')

        if scale_to_ef:
            if self.e_fermi is None:
                raise ValueError("The fermi energy is not set. Assign a value to `e_fermi`")
            return eig_array - self.e_fermi

        return eig_array

    def _validate_args(self, d_ene, e_max, e_min, smearing, scale_to_ef):

        eigs = self.get_eigs(scale_to_ef=scale_to_ef)

        if not isinstance(d_ene, (float, int)):
            raise ValueError("The `d_ene` must be an integer or a float")

        if e_max is None:
            e_max = eigs.max()
        else:
            if not isinstance(e_max, (float, int)):
                raise ValueError("The `e_max` must be an integer or a float")

        if e_min is None:
            e_min = eigs.min()
        else:
            if not isinstance(e_min, (float, int)):
                raise ValueError("The `e_min` must be an integer or a float")

        if e_min >= e_max:
            raise ValueError("The `e_min` must be smaller than `e_max`")

        if smearing is None:
            smearing = K2EV * 300
        else:
            if not isinstance(smearing, (float, int)):
                raise ValueError("The `smearing` must be an integer or a float")

        return e_max, e_min, smearing

    def compute_dos(
        self, d_ene=0.005, e_max=None, e_min=None, distribution="gaussian", smearing=None, scale_to_ef=False
    ):
        """
        Compute the dos from eigenvalues
        """
        import numpy as np
        from sisl.physics import get_distribution

        eigs = self.get_eigs(scale_to_ef=scale_to_ef)

        n_spin = len(eigs)
        n_kp = len(eigs[0])

        e_max, e_min, smearing = self._validate_args(d_ene, e_max, e_min, smearing, scale_to_ef)
        print(e_max, e_min, smearing)

        energies_bins = np.arange(e_min - smearing * 4, e_max + smearing * 4, d_ene)

        try:
            weights = self.get_kpoints(also_weights=True)[1]
            if len(weights) != len(eigs[0]):
                raise ValueError("Number of k-points in `self.get_eigs` and number of k-point's weights differ")
        except AttributeError:
            import logging
            logging.warning(" The kpoint's weights are not set!!! An equal weight for every kpoint is assumed.")
            weights = np.full(n_kp, 1 / n_kp)

        distr = get_distribution(distribution, smearing=smearing)

        dos_list = []
        for spin in range(n_spin):
            dos = np.zeros(len(energies_bins))
            for kpoint in range(n_kp):
                for ene in eigs[spin, kpoint, :]:
                    dos += distr(energies_bins - ene) * weights[kpoint]
            dos_list.append(dos)

        return energies_bins, np.array(dos_list)

        #plt.figure()
        #ax = plt.gca()
        #ax.set_title('DOS kT={:.1f} K'.format(kT * units('eV', 'K')))
        #ax.set_xlabel('E - Ef [eV]')
        #ax.set_xlim(E.min(), E.max())
        #ax.set_ylabel('DOS [1/eV]')
        #if ns._eigs.shape[0] == 2:
        #    for i, ud in enumerate(['up', 'down']):
        #        myplot(ax, ud, E, ns._eigs[i, :, :], ns._weight)
        #    plt.legend()
        #else:
        #    myplot(ax, '', E, ns._eigs[0, :, :], ns._weight)
        #if out is None:
        #    plt.show()
        #else:
        #    plt.savefig(out)


#
#       return p, namespace
