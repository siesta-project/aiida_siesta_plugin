""" Module with implementation of TKDict (translated-keys-dictionary) class.
    It is actually a dictionary with 'translation insensitive' keys. For
    example, in the FDFDict subclass:

        MD.TypeOfRun, md-type-of-run, mdTypeOfRun, mdtypeofrun

    all represent the same key in the dictionary. The actual form of the
    key returned by methods such as 'keys()' is the latest to be used in
    a setting operation.

    Vladimir Dikan and Alberto Garcia, 2017
   
"""
# Not compatible with python3 due to string-specific translation.
# TODO: use regular expressions for key translation

from string import maketrans
from collections import MutableMapping


class TKDict(MutableMapping):
    """ Dictionary-like class that also contains character translation and deletion data.
    Stores (value, initial-key) tuples accessible by a translated key.
    """

    @classmethod
    def translate_key(self, key):
        """ Definition of a rule for key translation. """
        raise NotImplementedError

    def __init__(self, *args, **kw):
        """ Create translated-keys-dictionary from initial data.
        If several input keys translate to same string, only first occurrence is saved.
        """
        # _storage is internal dictionary stored as: {<translated_key>: (<value>, <initial_key>), }
        self._storage = {}
        inp_dict = dict(*args, **kw)

        for inp_key in inp_dict.keys():
            self[inp_key] = inp_dict[inp_key]

    def keys(self):
        """ Return list of last key occurences. """
        # _storage keys are translated
        return [self.get_last_key(k) for k in self._storage.keys()]

    def iterkeys(self):
        'D.iterkeys() -> an iterator over the keys of D *** UPDATE'
        for key in self:
            value, last_key = self._storage[key]
            yield last_key
                                    
    def __setitem__(self, key, value):
        """ Store a (value, initial_key) tuple under translated key. """
        trans_key = self.translate_key(key)
        # check if we already have a translated key in _storage
        # if so, overwrite the value in tuple, but not the initial key
        self._storage.__setitem__(trans_key, (value, key))

    def __getitem__(self, key):
        """ Translate the key, unpack value-tuple and return the value if exists or None. """
        trans_key = self.translate_key(key)
        try:
            value, last_key = self._storage[trans_key]
            #self._storage.__setitem__(trans_key, (value, key))
            return value
        except KeyError:
            return None

    def iteritems(self):
        'D.iteritems() -> an iterator over the (key, value) items of D'
        for key in self:
            value, last_key = self._storage[key]
            yield (last_key, value)

    def items(self):
        "D.items() -> list of D's (key, value) pairs, as 2-tuples"
        return [(self._storage[key][1],self._storage[key][0]) for key in self]
    
    def get_last_key(self, key):
        """ Translate the key, unpack value-tuple and return
        the corresponding initial key if exists or None.
        """
        trans_key = self.translate_key(key)
        try:
            value, last_key = self._storage[trans_key]
            return last_key
        except KeyError:
            return None

    def __delitem__(self, key):
        """ Translate the key, purge value-tuple """
        self._storage.__delitem__(self.translate_key(key))

    def __iter__(self):
        return iter(self._storage)

    def __len__(self):
        return len(self._storage)

    def __repr__(self):
        return self._storage.__repr__()

    def __str__(self):
        return self._storage.__str__()


class FDict(TKDict):
    """ FDict class represents data from .fdf-file. """
    @classmethod
    def translate_key(self, key):
        # This will not work for unicode strings. See below
        return key.strip('-_')\
                  .translate(maketrans('_.', '--'), '')\
                  .lower()

class FDFDict(TKDict):
    """ FDFDict class represents data from .fdf-file. """
    @classmethod
    def translate_key(self, key):
        # There are incompatible 'translate'
        # methods for str and unicode objects... 
        
        if isinstance(key, unicode):
            # Unicode uses a single dictionary for translation
            to_remove = "-."
            table = {ord(char): None for char in to_remove}
            return key.translate(table).lower()
        else:
            # str accepts a translation map and an optional list of chars
            # to remove
            return key.translate(None, '-.').lower()
        
