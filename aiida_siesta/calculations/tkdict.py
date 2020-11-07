"""
Module with implementation of TKDict (translated-keys-dictionary) class.
It is actually a dictionary with 'translation insensitive' keys. For
example, in the FDFDict subclass:

    MD.TypeOfRun, md-type-of-run, mdTypeOfRun, mdtypeofrun

all represent the same key in the dictionary. The actual form of the
key returned by methods such as 'keys()' is the latest to be used in
a setting operation.

Vladimir Dikan and Alberto Garcia, 2017
Refined by Emanuele Bosoni in 2020

"""

from abc import abstractmethod
from collections.abc import MutableMapping


class TKDict(MutableMapping):
    """
    Abstract class that mocks a python dictionary but with some translation
    rules applied to the key of the dictionary.
    Its subclasses need to define the method `translate_key` that contains the translation rule.
    Once the translation rule is set. This abstract class
    provides all the methods for the management of the dictionary.
    Internally, it defines an attribute _storage that is a dictionary
    with the translated_keys as keys and the (value, original-key) tuples as values.
    """

    @classmethod
    @abstractmethod
    def translate_key(cls, key):
        """ Definition of a rule for key translation. """

    def __init__(self, inp_dict=None):  #*args, **kw):
        """
        Create translated-keys-dictionary from initial data. Initail data must be a dictionary.
        If several input keys translate to same string, only the last occurrence is saved.
        The class can also be instantiaced with no argument passed.
        """
        # _storage is internal dictionary stored as: {<translated_key>: (<value>, <initial_key>), }
        self._storage = {}

        if inp_dict is not None:
            if not isinstance(inp_dict, dict):
                message = 'invalid argument for `{}`: it only accepts a dictionary'.format(self.__class__.__name__)
                raise RuntimeError(message)

        if inp_dict is not None:
            for inp_key in inp_dict:
                self[inp_key] = inp_dict[inp_key]

    def __setitem__(self, key, value):
        """
        Store a (value, initial_key) tuple under translated key.
        """
        trans_key = self.translate_key(key)
        # check if we already have a translated key in _storage
        # if so, overwrite the value in tuple, but not the initial key
        self._storage.__setitem__(trans_key, (value, key))

    def __getitem__(self, key):
        """
        Translate the key, unpack value-tuple and return the value if exists.
        If the translated_key is not present in the dictionary, a KeyError
        will be reported, like in any normal python dictionary.
        """
        trans_key = self.translate_key(key)
        return self._storage[trans_key][0]

    def __delitem__(self, key):
        """ Translate the key, purge value-tuple """
        self._storage.__delitem__(self.translate_key(key))

    def __iter__(self):
        """We iter on the translated keys"""
        return iter(self._storage)

    def __len__(self):
        return len(self._storage)

    def __repr__(self):
        return self._storage.__repr__()

    def __str__(self):
        return self._storage.__str__()

    def values(self):
        """
        Return list of values.

        """
        return [self[k] for k in self]

    def untranslated_keys(self):
        """
        Return a list of last occurencies of the untranslated keys. For each translated_key,
        it is stored as second entry of the self._storage[translated_key] tuple.
        """
        return [self._storage[k][1] for k in self._storage]

    def keys(self):
        """
        Return list of translated keys. The self_storage keys are already translated.

        """
        return self._storage.keys()

    def untranslated_items(self):
        """
        Return a list of the items with the last occurencies of the untranslated keys
        as key.
        """
        return [(self._storage[key][1], self._storage[key][0]) for key in self]

    def items(self):
        """
        Return a list of the items with the translated key as key.

        """
        return [(key, self._storage[key][0]) for key in self]

    def get_dict(self):
        """
        Return a dictionary, where the key are the translated keys
        """
        return {key: val[0] for key, val in self._storage.items()}

    def get_untranslated_dict(self):
        """
        Return a dictionary, where the key are the translated keys
        """
        return {val[1]: val[0] for val in self._storage.values()}

    #Only for back compatibility!!! Does the same job of
    #items when iterated, but can't be called by itself
    def get_filtered_items(self):
        for k, v in self._storage.items():
            yield k, v[0]

    def get_last_untranslated_key(self, key):
        """
        Translate the key, unpack value-tuple and return
        the corresponding initial key if exists or None.
        """
        trans_key = self.translate_key(key)
        return self._storage[trans_key][1]
        #A KeyError here means the key is not defined


class FDFDict(TKDict):  # pylint: disable=too-many-ancestors
    """
    FDFDict class represents data from .fdf-file.

    This class follows a boilerplate raw python3-compatible translation rule.

    Behavior: drop dashes/dots/colons from key -> lowercase.
    (that keeps other special characters including underscores untouched)
    """

    @classmethod
    def translate_key(cls, key):
        to_remove = "-.:"

        if not isinstance(key, str):
            raise Exception("Key name error in FDFDict")

        # Unicode uses a single dictionary for translation
        table = {ord(char): None for char in to_remove}
        return key.translate(table).lower()
