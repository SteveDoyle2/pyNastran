from __future__ import print_function, absolute_import
from six import iteritems, itervalues, add_metaclass
from six.moves import range


_registered_node = {}


class H5NastranNodeMetaClass(type):
    def __new__(cls, clsname, bases, attrs):
        newclass = super(H5NastranNodeMetaClass, cls).__new__(cls, clsname, bases, attrs)

        if newclass.__name__ not in _registered_input_tables:
            assert newclass.version == (0, 0, 0), '%s version %r must be defined first!' % (
            newclass.__name__, (0, 0, 0))
            tmp = _registered_node[newclass.__name__] = {}
        else:
            tmp = _registered_node[newclass.__name__]

        assert newclass.version not in tmp
        tmp[newclass.version] = newclass

        # return last version defined, this allows newer versions to sublass last defined version as long as versions
        # are defined in order, although better to explicitly use KLS.get_version((i, j, k))
        return newclass


class H5NastranNode(object):
    version = (0, 0, 0)
    
    @classmethod
    def get_version(cls, version):
        return _registered_node[cls.__name__][version]
    
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
