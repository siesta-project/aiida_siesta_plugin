from aiida.plugins.entry_point import get_entry_point_from_class
from sisl import AtomicOrbital


class SislAtomicOrbital(AtomicOrbital):
    """
    This is a thin wrapper for sisl.AtomicOrbital. It just sets a dictionary "attributes"
    where the relevant info are summarized.
    Note: this class is not Data! It can not be stored! I is used in order to
    return in an easy way the orbitals from an IonData.
    """

    def __init__(self, *args, **kwargs):
        """
        The super is sisl.AtomicOrbital, therefore the instantiation follows its rules, as explained here:
        https://github.com/zerothi/sisl/blob/e28cdfff68d444139d54a755596de2a0285b0fd6/sisl/orbital.py#L689
        In addition we call `_set_orbital_dict`, that stores the info of the orbital in an internal dictionary.
        """
        super().__init__(*args, **kwargs)

        self.attributes = self._set_orbital_dict()

    def _set_orbital_dict(self):
        """
        Function called in the __init__, returns the orbital fetures in dictionary
        """

        orbital_dict = {}

        entry_point = get_entry_point_from_class(self.__class__.__module__, self.__class__.__name__)[1]
        if entry_point is not None:
            orbital_dict['orbital_type'] = entry_point.name

        orbital_dict["n"] = self.n
        orbital_dict["l"] = self.l
        orbital_dict["m"] = self.m
        orbital_dict["Z"] = self.Z
        orbital_dict["P"] = self.P
        orbital_dict["R"] = self.R
        orbital_dict["q0"] = self.q0
        #orbital_dict['spherical'] = (self.orb.__getstate__()["r"], self.orb.__getstate__()["f"])

        return orbital_dict
