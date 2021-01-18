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
        return self.get_attribute('ef', None)

    @e_fermi.setter
    def e_fermi(self, e_fermi):
        if not isinstance(e_fermi, (float, int)):
            raise ValueError("The ef must be an integer or a float")
        self.set_attribute('e_fermi', e_fermi)

    def set_eigs(self, eigs):
        """
        Set the eigenvalues, must be in eV
        """
        self.set_array('eigs', eigs)

    def get_eigs(self, scale_to_ef=False):
        """
        Get the eigenvalues
        """
        import numpy as np
        try:
            eig_array = np.array(self.get_array('eigs'))
        except KeyError:
            raise AttributeError('No stored eigs has been found')

        if scale_to_ef:
            if self.e_fermi is None:
                raise ValueError("The fermi energy is not set")
            return eig_array - self.e_fermi

        return eig_array

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
        if e_max is None:
            e_max = eigs.max()
        if e_min is None:
            e_min = eigs.min()
        if smearing is None:
            smearing = K2EV * 300

        energies_bins = np.arange(e_min - smearing * 4, e_max + smearing * 4, d_ene)

        try:
            weights = self.get_kpoints(also_weights=True)[1]
            #check wheigts lenght
        except AttributeError:
            #warining no kp set, weights
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
