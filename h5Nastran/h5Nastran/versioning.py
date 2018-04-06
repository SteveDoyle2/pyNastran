from __future__ import print_function, absolute_import
from six import iteritems, itervalues, add_metaclass
from six.moves import range

from collections import defaultdict


class VersioningData(object):
    def __init__(self):
        self._data = defaultdict(lambda: defaultdict(lambda: defaultdict(dict)))
        # first key is nastran_type, first value is dict
        # second key is nastran_version, second value is dict
        # third key is h5n_version, third value is dict
        # fourth key is class name, fourth value is the class

    def register(self, nastran_type, nastran_version, h5n_version, kls):
        if isinstance(nastran_type, str):
            nastran_type = nastran_type.lower()

        tmp = self._data[nastran_type]

        if nastran_version not in tmp:
            assert nastran_version == (0, 0, 0), '%s: nastran type %r, version %r must be defined first!' % (
                kls.__name__, nastran_type, (0, 0, 0)
            )

        tmp = tmp[nastran_version]

        if h5n_version not in tmp:
            assert h5n_version == (0, 0, 0), '%s: nastran type %r, version %r, h5n version %r must be defined first!' % (
                kls.__name__, nastran_type, nastran_version, (0, 0, 0)
            )

        tmp = tmp[h5n_version]

        assert kls.__name__ not in tmp, '%s: nastran type %r, version %r, h5n version %r already defined!' % (
            kls.__name__, nastran_type, nastran_version, h5n_version
        )

        tmp[kls.__name__] = kls

    def get_class(self, nastran_type, nastran_version, h5n_version, classname):
        if isinstance(nastran_type, str):
            nastran_type = nastran_type.lower()

        try:
            tmp = self._data[nastran_type]
        except KeyError:
            tmp = self._data[None]

        keys = sorted(tmp.keys())
        last_key = None

        for key in keys:
            if key <= nastran_version:
                last_key = key
            else:
                break

        if last_key is None:
            raise Exception('Cannot find %s: nastran type %s, version %r!'  %(
                classname, nastran_type, nastran_version
            ))

        tmp = tmp[last_key]

        keys = sorted(tmp.keys())
        last_key = None

        for key in keys:
            if key <= h5n_version:
                last_key = key
            else:
                break

        if last_key is None:
            raise Exception('Cannot find %s: nastran type %s, version %r, h5n version %r!' % (
                classname, nastran_type, nastran_version, h5n_version
            ))

        tmp = tmp[last_key]

        return tmp[classname]


class VersioningMetaClass(type):
    data = None

    def __new__(cls, clsname, bases, attrs):
        newclass = super(VersioningMetaClass, cls).__new__(cls, clsname, bases, attrs)
        cls.data.register(newclass.nastran_type, newclass.nastran_version, newclass.h5n_version, newclass)
        return newclass

    @classmethod
    def get_class(cls, *args):
        return cls.data.get_class(*args)
