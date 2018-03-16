from collections import defaultdict, OrderedDict
from copy import deepcopy

import numpy as np
import tables

from ..msc import data_tables
from ..data_helper import DataHelper

_defaults = {
    '<f8': DataHelper.default_double,
    '<i8': DataHelper.default_int
}


def _get_dtype(table):
    try:
        return data_tables[table.same_as].dtype
    except KeyError:
        return table.dtype


def _convert_dtype(dtype):
    new_dtype = []

    for i in range(len(dtype)-1):
        _dtype = dtype[i]

        if not isinstance(_dtype[1], list):
            new_dtype.append(_dtype)
        else:
            _shape = _dtype[2]
            for j in range(len(_dtype[1])):
                _ = _dtype[1][j]
                assert not isinstance(_, list)  # not doing multiple levels
                _name = _[0]
                _type = _[1]
                new_dtype.append((_name, _type, _shape))

    new_dtype.append(dtype[-1])  # add DOMAIN_ID

    return new_dtype


class TableDef(object):

    @classmethod
    def create(cls, table_def, defaults=None, len_id=None, pos_id=None, subtables=None, rename=None):
        if isinstance(table_def, str):
            table_def = data_tables[table_def]
        if subtables is None:
            subtables = [TableDef.create(data_tables.get(_, _), rename=rename) for _ in table_def.subtables]
        return cls(table_def.name, table_def.path, _get_dtype(table_def), defaults, len_id, pos_id, subtables, rename)

    def __init__(self, table_id, group, dtype, defaults=None, len_id=None, pos_id=None, subtables=None, rename=None):
        self.implemented = True
        self.table_id = table_id
        self.group = group

        self.dtype = np.dtype(dtype)

        try:
            self.Format = tables.descr_from_dtype(self.dtype)[0]
        except NotImplementedError:
            dtype = _convert_dtype(dtype)
            self.dtype = np.dtype(dtype)
            self.Format = tables.descr_from_dtype(self.dtype)[0]

        self.subtables = subtables

        self.attrs = list(self.dtype.names)

        try:
            self.attrs.remove('DOMAIN_ID')
        except ValueError:
            pass

        def _zero():
            return 0

        self._pos = defaultdict(_zero)

        _dtypes = {_[0]: (_[1], _[2]) for _ in dtype}

        for subtable in self.subtables:
            try:
                self.attrs.remove(subtable.pos_id)
                self.attrs.remove(subtable.len_id)
            except ValueError:
                print(self.path(), self.attrs, subtable.pos_id, subtable.len_id)
                raise

        self.table = None

        self.h5f = None

        if defaults is not None:
            self.defaults = defaults
        else:
            self.defaults = {}
            for i in range(len(self.attrs)):
                name = self.attrs[i]
                _type = _dtypes[name][0]
                if _type.startswith('S'):
                    self.defaults[name] = DataHelper.default_str
                else:
                    self.defaults[name] = _defaults[_type]

        if len_id is None:
            len_id = '%s_LEN' % self.table_id.replace('_', '')

        if pos_id is None:
            pos_id = '%s_POS' % self.table_id.replace('_', '')

        if rename is None:
            rename = {}

        self.len_id = rename.get(len_id, len_id)
        self.pos_id = rename.get(pos_id, pos_id)

        cols = set(self.dtype.names)

        for subtable in self.subtables:
            assert subtable.len_id in cols
            assert subtable.pos_id in cols

    def not_implemented(self):
        self.implemented = False

    def set_h5f(self, h5f):
        self.h5f = h5f
        for subtable in self.subtables:
            subtable.set_h5f(h5f)

    def _make_table(self, expected_rows=100000):
        if self.implemented is False:
            return None
        try:
            return self.h5f.get_node(self.path())
        except tables.NoSuchNodeError:
            try:
                self.h5f.create_table(self.group, self.table_id, self.Format, self.table_id,
                                        expectedrows=expected_rows, createparents=True)
            except tables.FileModeError:
                return None
            return self.h5f.get_node(self.path())

    def get_table(self):
        if self.table is None:
            self.table = self._make_table()
        return self.table

    def path(self):
        return self.group + '/' + self.table_id
                
    def write_data(self, data):
        identity = data['IDENTITY']
        names = list(self.dtype.names)

        if names[-1] == 'DOMAIN_ID':
            names.pop()

        data_len = len(identity[names[0]])

        if data_len == 0:
            return

        _data = np.empty(data_len, dtype=self.dtype)
        for name in names:
            try:
                _data[name] = identity[name]
            except (ValueError, TypeError) as e:
                print(name, identity[name], self.dtype, data_len, names[0], identity[names[0]])
                raise e

        # FIXME: bdf input domains need to be corrected, for now, just assuming there are no super elements
        #        need an easy way to determine super element ids, does pyNastran have a method to do so?
    
        try:
            _data['DOMAIN_ID'] = 1
        except (ValueError, KeyError):
            pass

        table = self.get_table()
        table.append(_data)
    
        subtable_ids = data.get('_subtables', [])
    
        for j in range(len(self.subtables)):
            subtable = self.subtables[j]
            subtable.write_data(data[subtable_ids[j]])
            
        table.flush()

    def read(self):
        try:
            table = self.h5f.get_node(self.path())
        except tables.NoSuchNodeError:
            table = None

        if table is None:
            return np.empty(0, dtype=self.dtype)

        return table.read()


def _get_value(value, default):
    if isinstance(value, (list, tuple, np.ndarray)):
        return [_get_value(_, default) for _ in value]

    if value in ('', None):
        return default

    return value


class InputTable(object):
    card_id = ''
    table_def = None  # type: TableDef

    def to_bdf(bdf):
        raise NotImplementedError

    def from_bdf(self, cards):
        raise NotImplementedError

    def __init__(self, h5n, parent):
        self._h5n = h5n
        self._parent = parent
        self._table_def = deepcopy(self.table_def)
        self._table_def.set_h5f(self._h5n.h5f)
        self.data = OrderedDict()

        if self.card_id is not None:
            if self.card_id == '':
                self.__class__.card_id = self.__class__.__name__

            self._h5n.register_card_table(self)

    def write_data(self, cards):
        if self.__class__.from_bdf is not InputTable.from_bdf:
            self._table_def.write_data(self.from_bdf(cards))
        else:
            raise NotImplementedError

    def finalize(self):
        pass

    def read(self):
        if self.__class__.from_bdf is InputTable.from_bdf:
            # if not implemented, then bail now so the table isn't created in the h5 file
            return np.empty(0, dtype=self.table_def.dtype)

        self.data.clear()
        # main table should ALWAYS be identity, regardless of inconsistencies in MSC spec
        self.data['identity'] = self._table_def.read()
        self.data[self._table_def.table_id.lower()] = self.data['identity']  # set to actual table id for convenience

        def _get_subtable_data(subtable):
            assert subtable.table_id not in self.data, subtable.table_id
            self.data[subtable.table_id.lower()] = subtable.read()

            for _ in subtable.subtables:
                _get_subtable_data(_)

        for subtable in self._table_def.subtables:
            _get_subtable_data(subtable)

        return self.data['identity']

    def __getattr__(self, item):
        if len(self.data) == 0:
            self.read()
        try:
            return self.data[item]
        except KeyError:
            raise AttributeError('Attribute %s not found on CardTable.' % (item))
