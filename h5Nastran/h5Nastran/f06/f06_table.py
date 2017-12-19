from __future__ import print_function, absolute_import
from six import iteritems, iterkeys, itervalues, add_metaclass
from six.moves import range


_tables = set()


def _register_class(kls):
    _tables.add(kls)


class RegisterClass(type):
    def __new__(cls, clsname, bases, attrs):
        newclass = super(RegisterClass, cls).__new__(cls, clsname, bases, attrs)
        _register_class(newclass)
        return newclass


@add_metaclass(RegisterClass)
class F06Table(object):

    _last_table = None

    @classmethod
    def find_table(cls, table_lines):
        if cls._last_table is not None:
            if cls._last_table.is_match(table_lines):
                return cls._last_table

        for table in _tables:
            if table.is_match(table_lines):
                cls._last_table = table
                return table

        return None

    @classmethod
    def is_match(cls, table_lines):
        return False

    def __init__(self):
        self.header = []
        self.data = []
        self.line_number = -1

    def to_punch(self):
        raise NotImplementedError

    def set_data(self, table_lines):
        raise NotImplementedError

    def _load_factor(self):
        raise NotImplementedError
