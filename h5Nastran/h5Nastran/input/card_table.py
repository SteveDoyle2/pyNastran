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
            subtables = [TableDef.create(data_tables[_], rename=rename) for _ in table_def.subtables]
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

    def write_subdata(self, data):
        names = []

        is_dict = False

        if isinstance(data, TableData):
            data_len = len(data.data)
        elif isinstance(data, dict):
            names = list(self.dtype.names)
            data_len = len(data['IDENTITY'][names[0]])
            is_dict = True
        else:
            raise Exception

        if data_len == 0:
            return

        table = self.get_table()
        table_row = table.row

        if not is_dict:
            for i in range(len(data.data)):
                _write_data_to_table(self, table_row, data.data[i])

                for j in range(len(self.subtables)):
                    if data.subdata_len[i, j] == 0:
                        continue
                    subtable = self.subtables[j]
                    table_row[subtable.len_id] = data.subdata_len[i, j]
                    table_row[subtable.pos_id] = self._pos[subtable.pos_id]
                    self._pos[subtable.pos_id] += data.subdata_len[i, j]

                table_row.append()

            for j in range(len(self.subtables)):
                subdata = data.subdata[j]
                subtable = self.subtables[j]
                subtable.write_subdata(subdata)

        else:
            identity = data['IDENTITY']
            names = list(self.dtype.names)
            data_len = len(identity[names[0]])
            _data = np.empty(data_len, dtype=self.dtype)
            for name in names:
                _data[name] = identity[name]

            table.append(_data)

            subtable_ids = data.get('_subtables', [])

            for j in range(len(self.subtables)):
                subtable = self.subtables[j]
                subtable.write_subdata(data[subtable_ids[j]])

    def write_data(self, cards, domain, from_bdf):
        ids = sorted(cards.keys())

        table = self.get_table()
        table_row = table.row

        for card_id in ids:
            card = cards[card_id]

            # if isinstance(card, list):
            #     card_data = [from_bdf(card[i]).data[0] for i in range(len(card))]
            #     self.write_subdata(TableData(card_data))
            #     continue

            data = from_bdf(card)
            if isinstance(data, TableData):
                for i in range(len(data.data)):
                    _write_data_to_table(self, table_row, data.data[i])

                    for j in range(len(self.subtables)):
                        if data.subdata_len[i, j] == 0:
                            continue
                        subtable = self.subtables[j]
                        table_row[subtable.len_id] = data.subdata_len[i, j]
                        table_row[subtable.pos_id] = self._pos[subtable.pos_id]
                        self._pos[subtable.pos_id] += data.subdata_len[i, j]

                    table_row['DOMAIN_ID'] = domain
                    table_row.append()

                for j in range(len(self.subtables)):
                    subdata = data.subdata[j]
                    subtable = self.subtables[j]
                    subtable.write_subdata(subdata)
            elif isinstance(data, dict):
                identity = data['IDENTITY']
                names = list(self.dtype.names)[:-1]
                data_len = len(identity[names[0]])
                _data = np.empty(data_len, dtype=self.dtype)
                for name in names:
                    _data[name] = identity[name]

                _data['DOMAIN_ID'] = domain

                table.append(_data)

                subtable_ids = data.get('_subtables', [])

                for j in range(len(self.subtables)):
                    subtable = self.subtables[j]
                    subtable.write_subdata(data[subtable_ids[j]])
            else:
                raise Exception

        self.h5f.flush()

    def read(self):
        table = self.get_table()
        if table is None:
            return np.empty(0, dtype=self.dtype)
        return table.read()


def _get_value(value, default):
    if isinstance(value, (list, tuple, np.ndarray)):
        return [_get_value(_, default) for _ in value]

    if value in ('', None):
        return default

    return value


def _write_data_to_table(table_def, table_row, data, append=False):
    attrs = table_def.attrs
    defaults = table_def.defaults

    # if len(data) == 0:
    #     return

    for i in range(len(attrs)):
        attr = attrs[i]
        table_row[attr] = _get_value(data[i], defaults[attr])

    if append:
        table_row.append()


class TableData(object):
    def __init__(self, data=None, subdata_len=None, subdata=None):

        if data is None:
            data = []

        if subdata_len is None:
            subdata_len = np.zeros((1, 1))

        if subdata is None:
            subdata = []

        self.data = data
        self.subdata_len = subdata_len  # type: np.ndarray
        self.subdata = subdata  # type: list[TableData]

    def __repr__(self):
        result = [str(self.data)]
        for subdata in self.subdata:
            result.append(subdata.__repr__())
        return ';'.join(result)

    def validate(self):
        if len(self.subdata) > 0:
            shape = self.subdata_len.shape
            assert shape[0] == len(self.data), (shape[0], len(self.data))
            assert shape[1] == len(self.subdata), (shape[1], len(self.subdata))
        for subdata in self.subdata:
            subdata.validate()


class CardTable(object):
    card_id = ''
    table_def = None  # type: TableDef

    @staticmethod
    def from_bdf(card):
        raise NotImplementedError

    @staticmethod
    def to_bdf(data):
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

    def write_data(self, cards, domain):
        if self.from_bdf is CardTable.from_bdf:
            self._table_def.not_implemented()
            raise NotImplementedError
        self._table_def.write_data(cards, domain, self.from_bdf)

    def finalize(self):
        pass

    def read(self):
        if self.from_bdf is CardTable.from_bdf:
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
