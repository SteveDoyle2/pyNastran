from __future__ import print_function, absolute_import
from six import iteritems, itervalues, add_metaclass
from six.moves import range

from .versioning import VersioningData, VersioningMetaClass



class H5NastranNodeMetaClass(VersioningMetaClass):
    data = VersioningData()


@add_metaclass(H5NastranNodeMetaClass)
class H5NastranNode(object):
    nastran_type = None
    nastran_version = (0, 0, 0)
    h5n_version = (0, 0, 0)

    def __new__(cls, h5n, *args):
        kls = cls.get_class(h5n.nastran_type, h5n.nastran_version, h5n.h5n_version)
        return object.__new__(kls)
    
    @classmethod
    def get_class(cls, nastran_type, nastran_version, h5n_version):
        return H5NastranNodeMetaClass.get_class(nastran_type, nastran_version, h5n_version, cls.__name__)
    
    def path(self):
        raise NotImplementedError

    def read(self):
        for key, item in iteritems(self.__dict__):
            if key.startswith('_'):
                continue
            try:
                item.read()
            except AttributeError:
                pass

    def update(self):
        raise NotImplementedError
