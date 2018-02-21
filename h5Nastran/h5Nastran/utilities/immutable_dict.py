"""
This is used mainly so the data contained in the dict isn't inadvertently mutated.

"""

from collections import Mapping

try:
    from types import MappingProxyType
except ImportError:
    # workaround for https://github.com/SteveDoyle2/pyNastran/issues/481
    # might just get rid of the ImmutableDict in the future, if you change the data in the dict, user beware
    def MappingProxyType(obj):
        return obj


class ImmutableDict(Mapping):
    """
    Copies a dict and proxies it via types.MappingProxyType to make it immutable.
    https://stackoverflow.com/questions/9997176/immutable-dictionary-only-use-as-a-key-for-another-dictionary/39673094#39673094
    """
    def __init__(self, somedict):
        dictcopy = dict(somedict)  # make a copy
        self._dict = MappingProxyType(dictcopy)  # lock it
        self._hash = None

    def __getitem__(self, key):
        return self._dict[key]

    def __len__(self):
        return len(self._dict)

    def __iter__(self):
        return iter(self._dict)

    def __hash__(self):
        if self._hash is None:
            self._hash = hash(frozenset(self._dict.items()))
        return self._hash

    def __eq__(self, other):
        try:
            return self._dict == other._dict
        except AttributeError:
            return False

    def __repr__(self):
        return str(self._dict)
