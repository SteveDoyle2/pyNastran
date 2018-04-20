from __future__ import print_function, absolute_import

from collections import defaultdict
from copy import deepcopy

import numpy as np
import pandas as pd
import tables
from six import iteritems, add_metaclass
from typing import List

from h5Nastran.msc import data_tables
from h5Nastran.post_process.result_readers.punch import PunchTableData
from h5Nastran.versioning import VersioningMetaClass, VersioningData
from .result_data import ResultData

pd.options.mode.chained_assignment = None


########################################################################################################################


def _get_dtype(table):
    try:
        return data_tables[table.same_as].dtype
    except KeyError:
        return table.dtype


########################################################################################################################


class _Format(object):
    def __init__(self):
        raise NotImplementedError


class IndexFormat(tables.IsDescription):
    DOMAIN_ID = tables.Int64Col(pos=1)
    POSITION = tables.Int64Col(pos=2)
    LENGTH = tables.Int64Col(pos=3)


class PrivateIndexFormat(tables.IsDescription):
    LOCATION = tables.Int64Col(pos=1)
    LENGTH = tables.Int64Col(pos=2)
    OFFSET = tables.Int64Col(pos=3)


private_index_format_dtype = tables.dtype_from_descr(PrivateIndexFormat)


class PrivateIndexDataFormat(tables.IsDescription):
    ID = tables.Int64Col(pos=1)


private_index_data_format_dtype = tables.dtype_from_descr(PrivateIndexDataFormat)

########################################################################################################################


def _validator(data):
    return data

########################################################################################################################


def _convert_dtype(dtype, indices):
    new_dtype = []
    new_indices = DataGetter()

    for i in range(len(dtype)-1):
        _dtype = dtype[i]
        _indices = indices.indices[i]

        if not isinstance(_dtype[1], list):
            new_dtype.append(_dtype)
            new_indices.indices.append(_indices)
        else:
            for j in range(len(_dtype[1])):
                _ = _dtype[1][j]
                assert not isinstance(_, list)  # not doing multiple levels
                _name = _[0]
                _type = _[1]
                if len(_indices) == 1:
                    _shape = ()
                    new_indices.indices.append(_indices[0][j])
                else:
                    _shape = (len(_indices),)
                    new_indices.indices.append([_[j] for _ in _indices])

                new_dtype.append((_name, _type, _shape))

    new_dtype.append(dtype[-1])  # add DOMAIN_ID

    return new_dtype, new_indices

########################################################################################################################


class DefinedValue(object):
    def __init__(self, value):
        self.value = value


class DataGetter(object):
    def __init__(self, dtype=None, indices=None, indices_len=None):
        self.indices = []

        if dtype is not None:
            self.make_indices_from_dtype(dtype)

        if indices is not None:
            self.indices = deepcopy(indices)

        self._indices_len = indices_len

    def __len__(self):
        if self._indices_len is None:
            return len(self.indices)
        else:
            return self._indices_len

    def make_indices_from_dtype(self, dtype):
        del self.indices[:]

        def _make_indices(_dtype, i):
            indices = []

            for _ in _dtype:
                if i == 1:
                    i = 2  # skip 2nd field in punch file by default

                _d = _[1]
                shape = _[2]
                try:
                    size = shape[0]
                except IndexError:
                    size = 1

                if isinstance(_d, list):
                    _indices = _make_indices(_d, i)

                    if size > 1:
                        _list_indices = [_indices]
                        for j in range(size-1):
                            _list_indices.append((np.array(_list_indices[-1]) + len(_list_indices[-1])).tolist())
                        _indices = _list_indices
                    else:
                        _indices = [_indices]

                    def _count(_indices_):
                        _result = 0
                        for _i in _indices_:
                            if isinstance(_i, list):
                                _result += _count(_i)
                            else:
                                _result += 1

                        return _result

                    i += _count(_indices)
                    indices.append(_indices)
                else:
                    indices.append(i)
                    i += size

            return indices

        self.indices.extend(_make_indices(dtype, 0))
        self.indices.pop()

    def get_data(self, data, indices=None):
        result = []

        if indices is None:
            indices = self.indices

        for i in indices:
            try:
                # print(i)
                if isinstance(i, list):
                    result.append(self.get_data(data, i))
                elif isinstance(i, (int, slice)):
                    result.append(data[i])
                elif isinstance(i, DefinedValue):
                    result.append(i.value)
                else:
                    result.append(data[i])
            except IndexError:
                raise IndexError((i, data))

        return result


def _get_data(data, index):
    if isinstance(index, (int, slice)):
        return data[index]
    elif isinstance(index, (list, tuple)):
        return [_get_data(data, i) for i in index]
    elif isinstance(index, DefinedValue):
        return index.value
    elif isinstance(index, DataGetter):
        return index.get_data(data)
    else:
        raise TypeError('Unknown index type! %s' % str(type(index)))

########################################################################################################################


# TODO: merge TableDef into ResultTable
class TableDef(object):
    
    result_table_data = None  # type: ResultTableData

    @classmethod
    def create(cls, table_def, results_type, indices=None, validator=None, len_id=None, pos_id=None, subtables=None, rename=None,
               is_subtable=False):
        if isinstance(table_def, str):
            table_def = data_tables[table_def]
        try:
            dtype = data_tables[table_def.same_as].dtype
        except KeyError:
            dtype = table_def.dtype
        index_id = dtype[0][0]
        if indices is None:
            indices = DataGetter(dtype)
        for _ in dtype:
            if isinstance(_[1], list):
                # need to convert dtype since pytables doesn't support nested dtypes
                dtype, indices = _convert_dtype(dtype, indices)
                break
        if subtables is None:
            subtables = [TableDef.create(data_tables[_], '', rename=rename, is_subtable=True) for _ in table_def.subtables]
        return cls(table_def.name, table_def.path, results_type, index_id, dtype, indices, validator,
                   len_id, pos_id, subtables, rename, is_subtable=is_subtable)

    def __init__(self, table_id, group, results_type, index_id, dtype, indices, validator=None,
                 len_id=None, pos_id=None, subtables=None, rename=None, is_subtable=False):
        self.implemented = True
        self.table_id = table_id
        self.group = group

        self.dtype = np.dtype(dtype)

        if subtables is None:
            subtables = []

        assert len(subtables) == 0, 'Cannot have subtables defined for this result table type.'

        self.subtables = subtables

        self.attrs = list(self.dtype.names)

        try:
            self.attrs.remove('DOMAIN_ID')
        except ValueError:
            pass

        def _zero():
            return 0

        self._pos = defaultdict(_zero)

        for subtable in self.subtables:
            try:
                self.attrs.remove(subtable.pos_id)
                self.attrs.remove(subtable.len_id)
            except ValueError:
                print(self.path(), self.attrs, subtable.pos_id, subtable.len_id)
                raise

        self.table = None

        self.h5f = None
        self.h5n = None

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

        ################################################################################################################

        self.is_subtable = is_subtable

        self.index_id = index_id

        try:
            self.results_type = results_type[0]
        except IndexError:
            self.results_type = results_type

        if validator is None:
            validator = _validator

        self.validator = validator

        try:
            self.Format = tables.descr_from_dtype(self.dtype)[0]
        except NotImplementedError:
            dtype, indices = _convert_dtype(dtype, indices)
            self.dtype = np.dtype(dtype)
            self.Format = tables.descr_from_dtype(self.dtype)[0]

        self.indices = deepcopy(indices)

        assert isinstance(self.indices, DataGetter)

        assert len(self.indices) == len(self.dtype.names) - 1, str((len(self.indices), len(self.dtype.names)))

        self.domain_count = 0

        self.IndexFormat = IndexFormat
        self.PrivateIndexFormat = PrivateIndexFormat
        self.PrivateIndexDataFormat = PrivateIndexDataFormat

        self._index_table = None
        self._private_index_table = None

        self._index_data = []
        self._subcase_index = []
        self._index_offset = 0

        self._subcase_ids = set()
        
        self._index_options = {}
        
        self.result_table = None
    
    def add_index_option(self, option, indices):
        assert isinstance(indices, (DataGetter, type(None)))
        self._index_options[option] = indices

    def finalize(self):
        if self.is_subtable:
            return
        self._write_index()
        self._write_private_index()

    def get_table(self):
        if self.table is None:
            self.table = self._make_table()
        return self.table
    
    def set_result_table(self, result_table):
        self.result_table = result_table

    def not_implemented(self):
        self.implemented = False

    def read(self, indices):
        table = self.get_table()

        indices = np.array(indices, dtype='i8')

        indices_len = len(indices)

        data = np.empty(indices_len, dtype=table._v_dtype)

        if indices_len > 0:
            table._read_elements(indices, data)

        return data

    def search(self, data_ids, domains=(), **kwargs):
        private_index_table = self._get_private_index_table()
        
        if len(domains) == 0:
            domains = list(range(len(private_index_table)))

        indices = set()

        for domain in domains:
            try:
                index_dict, offset = private_index_table[domain-1]
            except IndexError:
                continue

            index_dict_get = index_dict.get

            for data_id in data_ids:
                _indices = index_dict_get(data_id, {})
                indices.update(set(index + offset for index in _indices))

        results = self.read(sorted(indices))
        
        return self.result_table_data.from_records(results)

    def path(self):
        return self.group + '/' + self.table_id

    def set_h5n(self, h5n):
        self.h5f = h5n.h5f
        self.h5n = h5n
        for subtable in self.subtables:
            subtable.set_h5n(h5n)

    def write_punch_data(self, data):
        assert isinstance(data, PunchTableData)

        options = data.header.options

        indices = self.indices

        for option in options:
            try:
                _indices = self._index_options[option]
                if _indices is not None:
                    indices = _indices
            except KeyError:
                msg = """
                Result table '%s' is not supported!
                A parameter in your bdf file directed Nastran to output the result table with an option that is
                currently not supported.  This option might affect the format of the results file.  H5Nastran
                needs to know the format.
                This requires the following to be added to the result table definition:
                table_def.add_index_option('%s', new_format)
                where new_format might simply be None (no change to default format).
                See RESULT/ELEMENTAL/STRESS/QUAD4 in the source code for an example.
                Please create a new issue on github.com/SteveDoyle2/pyNastran to request this format to be supported.
                """ % (data.header.results_type, option)
                raise H5NastranException(msg)

        data_ = data.data

        result = np.empty(len(data_), dtype=self.dtype)

        validator = self.validator

        names = list(self.dtype.names)

        _result = {name: result[name] for name in names}

        for i in range(len(data_)):
            try:
                _data = validator(indices.get_data(data_[i]))
            except IndexError:
                lens = [len(_) for _ in data_]
                try:
                    data_0 = data_[i-1]
                except IndexError:
                    data_0 = None

                data_1 = data_[i]

                try:
                    data_2 = data_[i+1]
                except IndexError:
                    data_2 = None

                raise IndexError((lens, data_0, data_1, data_2))

            for j in range(len(names)-1):
                _result[names[j]][i] = _data[j]

        result_data = ResultData()
        result_data.data = result
        result_data.options.update(options)

        result_data.set_result_type(data.header.results_type_basic)
        result_data.subcase_id = data.header.subcase_id
        result_data.punch_results()

        self.write_data(result_data)

    def write_op2_data(self, data):
        raise NotImplementedError

    def write_data(self, data):
        assert isinstance(data, ResultData)

        names = list(self.dtype.names)

        data_len = len(data.data[names[0]])

        if data_len == 0:
            return

        table = self.get_table()

        if isinstance(data.data, np.ndarray) and data.data.dtype == self.dtype:
            _data = data.data
        else:
            _data = np.empty(data_len, dtype=self.dtype)
            for name in names:
                _data[name][:] = data.data[name]

        subcase_id = data.subcase_id

        if subcase_id not in self._subcase_ids:
            self.domain_count += 1
            self._subcase_ids.add(subcase_id)

        _data['DOMAIN_ID'][:] = self.domain_count

        self.result_table.apply_options(data)

        table.append(_data)
        self._record_data_indices(_data)

        self.h5f.flush()

    def _get_index_table(self):
        if self.is_subtable:
            return None
        h5f = self.h5f
        if self._index_table is None:
            index_table = h5f.get_node('/INDEX%s' % self.path())
            index_table = index_table.read()

            self._index_table = [
                set(range(index_table['POSITION'][i], index_table['POSITION'][i] + index_table['LENGTH'][i]))
                for i in range(index_table.shape[0])
            ]

        return self._index_table

    def _get_private_index_table(self):
        if self.is_subtable:
            return None
        h5f = self.h5f
        if self._private_index_table is None:
            data = h5f.get_node(self._private_index_path + '/DATA')
            data = data.read()

            identity = h5f.get_node(self._private_index_path + '/IDENTITY')
            identity = identity.read()

            private_index_table = self._private_index_table = []

            _index_data = {}

            for i in range(identity.shape[0]):
                location, length, offset = identity[i]

                try:
                    private_index_table.append((_index_data[location][0], offset))
                except KeyError:
                    _data_dict = load_data_dict(data['ID'][location: location + length])
                    _index_data[location] = (_data_dict, offset)
                    private_index_table.append(_index_data[location])

        return self._private_index_table

    def _get_private_index_tables(self):
        if self.is_subtable:
            return None, None
        h5f = self.h5f
        try:
            identity = h5f.get_node(self._private_index_path + '/IDENTITY')
        except tables.NoSuchNodeError:
            identity = h5f.create_table(self._private_index_path, 'IDENTITY', self.PrivateIndexFormat, 'Private Index',
                                        expectedrows=len(self._subcase_index), createparents=True)
        try:
            data = h5f.get_node(self._private_index_path + '/DATA')
        except tables.NoSuchNodeError:
            data = h5f.create_table(self._private_index_path, 'DATA', self.PrivateIndexDataFormat, 'Private Index Data',
                                    expectedrows=sum([len(_) for _ in self._index_data]), createparents=True)

        return identity, data

    def _make_table(self, expected_rows=100000):
        if self.implemented is False:
            return None
        try:
            return self.h5f.get_node(self.path())
        except tables.NoSuchNodeError:
            try:
                self.h5f.create_table(self.group, self.table_id, self.Format, self.results_type,
                                        expectedrows=expected_rows, createparents=True)
            except tables.FileModeError:
                return None
            return self.h5f.get_node(self.path())

    @property
    def _private_index_path(self):
        return self.h5n.table_paths.private_index_path + self.path()

    def _record_data_indices(self, data):
        if self.is_subtable:
            return

        serialized_data = serialize_indices(data[self.index_id][:])

        index_data = serialized_data.astype(dtype=private_index_data_format_dtype)

        found_index = False

        index_data_offset = 0

        for i in range(len(self._index_data)):
            _index_data = self._index_data[i]
            if _index_data.shape == index_data.shape:
                if np.all(_index_data == index_data):
                    self._subcase_index.append((index_data_offset, index_data.shape[0], self._index_offset))
                    found_index = True
                    break
            index_data_offset += _index_data.shape[0]

        if not found_index:
            self._index_data.append(index_data)
            self._subcase_index.append((index_data_offset, index_data.shape[0], self._index_offset))

        self._index_offset += data.shape[0]

    def _write_index(self):
        if self.is_subtable:
            return

        table = self.get_table()

        # noinspection PyProtectedMember
        domain_id = table.cols._f_col('DOMAIN_ID')[:]

        unique, counts = np.unique(domain_id, return_counts=True)
        counts = dict(zip(unique, counts))

        domains = self.h5f.create_table('/INDEX' + self.group, self.table_id, self.IndexFormat, self.results_type,
                                   expectedrows=len(counts), createparents=True)

        row = domains.row

        pos = 0

        for i in range(len(counts)):
            d = i + 1
            try:
                count = counts[d]
            except KeyError:
                raise KeyError(str((d, counts)))

            row['DOMAIN_ID'] = d
            row['POSITION'] = pos
            row['LENGTH'] = count

            row.append()

            pos += count

        domains.flush()

        self._index_table = None

    def _write_private_index(self):
        if self.is_subtable:
            return

        identity, data = self._get_private_index_tables()

        identity.append(np.array(self._subcase_index, dtype=private_index_format_dtype))

        for index_data in self._index_data:
            data.append(index_data)

        self.h5f.flush()

        self._private_index_table = None

        del self._subcase_index[:]
        del self._index_data[:]

########################################################################################################################


def serialize_indices(data):
    return serialize_data_dict(get_data_dict(data))


def get_data_dict(data):
    from collections import OrderedDict
    data_dict = OrderedDict()

    for i in range(data.shape[0]):
        data_id = int(data[i])
        try:
            data_dict[data_id].append(i)
        except KeyError:
            data_dict[data_id] = [i]

    return data_dict


def serialize_data_dict(data_dict):
    data = []

    for key, _data in iteritems(data_dict):
        _data_ = [key, len(_data)]
        _data_.extend(_data)
        data.extend(_data_)

    return np.array(data)


def load_data_dict(serialize_data, offset=0):
    data_dict = {}

    last_i = serialize_data.shape[0] - 1

    i = 0
    while True:
        data_id = serialize_data[i]
        data_len = serialize_data[i + 1]

        _data = serialize_data[i + 2: i + 2 + data_len]

        data_dict[data_id] = _data + offset

        i += 2 + data_len

        if i >= last_i:
            break

    assert i == last_i + 1

    return data_dict


########################################################################################################################


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


########################################################################################################################


class ResultTableVersioningData(VersioningData):
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

        result_type = kls.result_type

        if isinstance(result_type, str):
            result_type = [result_type]

        for res_type in result_type:
            tmp[res_type] = kls


class ResultTableMetaClass(VersioningMetaClass):
    data = ResultTableVersioningData()


# class ResultTableMetaClass(type):
#     def __new__(cls, clsname, bases, attrs):
#         newclass = super(ResultTableMetaClass, cls).__new__(cls, clsname, bases, attrs)
#
#         result_type = newclass.result_type
#
#         if isinstance(result_type, str):
#             result_type = [result_type]
#
#         for res_type in result_type:
#             if res_type not in _registered_result_tables:
#                 assert newclass.version == (0, 0, 0), '%s version %r must be defined first!' % (
#                 newclass.__name__, (0, 0, 0))
#                 tmp = _registered_result_tables[res_type] = {}
#             else:
#                 tmp = _registered_result_tables[res_type]
#
#             assert newclass.version not in tmp
#             tmp[newclass.version] = newclass
#
#         # return last version defined, this allows newer versions to sublass last defined version as long as versions
#         # are defined in order, although better to explicitly use KLS.get_version((i, j, k))
#         return newclass


@add_metaclass(ResultTableMetaClass)
class ResultTable(object):
    result_type = ''
    table_def = None  # type: TableDef
    result_data_cols = []  # type: List[int]
    result_data_group_by = []  # type: List[str]
    nastran_type = None
    nastran_version = (0, 0, 0)
    h5n_version = (0, 0, 0)

    @classmethod
    def get_class(cls, nastran_type, nastran_version, h5n_version):
        # TODO: make result_type only a string
        # this will only work if cls.result_type is a string
        return InputTableMetaClass.get_class(nastran_type, nastran_version, h5n_version, cls.result_type)

    def __init__(self, h5n, parent):
        self._h5n = h5n
        self._parent = parent
        self._table_def = deepcopy(self.table_def)
        self._table_def.set_h5n(self._h5n)
        self._table_def.set_result_table(self)
        
        if len(self.result_data_cols) == 0:
            dtype = self.table_def.dtype
            names = dtype.names
            data_cols = []
            for name in names:
                if dtype[name] == 'float64':
                    data_cols.append(name)
            self.__class__.result_data_cols.extend(data_cols)

        if len(self.result_data_group_by) == 0:
            dtype = self.table_def.dtype
            names = dtype.names
            data_cols = []
            for name in names:
                if dtype[name] != 'float64':
                    data_cols.append(name)
            self.__class__.result_data_group_by.extend(data_cols)

        class _ResultTableData(ResultTableData):
            data_cols = pd.Index(self.result_data_cols)
            data_group_by = list(self.result_data_group_by)
        
        self._table_def.result_table_data = _ResultTableData

        if self.result_type is not None:
            self._h5n.register_result_table(self)

    def write_punch_data(self, data):
        # type: (PunchTableData) -> None
        self._table_def.write_punch_data(data)

    def write_op2_data(self, data):
        self._table_def.write_op2_data(data)

    def write_data(self, data):
        self._table_def.write_data(data)

    def finalize(self):
        self._table_def.finalize()

    def read(self, indices):
        return self._table_def.read(indices)

    def search(self, data_ids, domains=(), **kwargs):
        try:
            return self._table_def.search(data_ids, domains, **kwargs)
        except tables.exceptions.NoSuchNodeError:
            return self._table_def.result_table_data()

    @property
    def results_type(self):
        return self._table_def.results_type

    @results_type.setter
    def results_type(self, value):
        self._table_def.results_type = value

    def path(self):
        return self._parent.path() + [self.__class__.__name__]

    def read_h5_table(self):
        path = '/'.join(self.path())

        data = self._h5n.h5f.get_node(path).read()

        return data
    
    def apply_options(self, data):
        return data


########################################################################################################################


class ResultTableData(pd.DataFrame):
    _metadata = ['data_cols', 'data_group_by', 'the_dtype']

    data_cols = []  # List[str]
    data_group_by = []  # List[str]
    the_dtype = None
    
    @classmethod
    def from_records(cls, nparr):
        try:
            result = super(ResultTableData, cls).from_records(nparr)  # .tolist(), columns=nparr.dtype.names)
        except Exception:
            result = super(ResultTableData, cls).from_records(nparr.tolist(), columns=nparr.dtype.names)
        cls.the_dtype = nparr.dtype
        return result

    @property
    def _constructor(self):
        return self.__class__

    def to_numpy(self):
        data_size = len(self[list(self.keys())[0]])
        result = np.empty(data_size, dtype=self.the_dtype)

        for col in self.keys():
            res = result[col]
            data = self[col]

            for i in range(data_size):
                res[i] = data[i]

        return result

    def _data_group_by(self):
        if len(self.data_group_by) == 0:
            keys = list(self.keys())
            self.__class__.data_group_by = [keys[0], keys[-1]]

        return self.data_group_by

    def _data_cols(self):
        if len(self.data_cols) == 0:
            self.__class__.data_cols = self.select_dtypes(include=['float']).columns

        return self.data_cols

    def __truediv__(self, other):
        assert other != 0.
        return self.__mul__(1 / other)

    def __div__(self, other):
        return self.__truediv__(other)

    def __itruediv__(self, other):
        assert other != 0.
        return self.__imul__(1 / other)

    def __idiv__(self, other):
        return self.__itruediv__(other)

    def __mul__(self, other):
        assert isinstance(other, (float, int))
        result = self.copy(deep=True)

        for i in self._data_cols():
            result[i] *= other

        return result

    def __rmul__(self, other):
        return self.__mul__(other)

    def __imul__(self, other):
        assert isinstance(other, (float, int))

        for i in self._data_cols():
            self[i] *= other

        return self

    def __add__(self, other):
        assert np.all(other.dtypes == self.dtypes)
        result = self.copy(deep=True).append(other)
        data_group_by = self._data_group_by()
        result = result.groupby(data_group_by, as_index=False).sum()
        result.sort_values(list(reversed(data_group_by)), inplace=True)
        result = self.__class__(result[self.keys()])
        return result
    
    def __sub__(self, other):
        assert np.all(other.dtypes == self.dtypes)
        result = self.copy(deep=True).append(-1 * other)
        data_group_by = self._data_group_by()
        result = result.groupby(data_group_by, as_index=False).sum()
        result.sort_values(list(reversed(data_group_by)), inplace=True)
        result = self.__class__(result[self.keys()])
        return result

    def __rsub__(self, other):
        raise NotImplementedError

    def __iadd__(self, other):
        raise NotImplementedError
    
    def __isub__(self, other):
        raise NotImplementedError

    def __radd__(self, other):
        raise NotImplementedError

    def __rtruediv__(self, other):
        raise NotImplementedError

    def __rdiv__(self, other):
        return NotImplementedError


class H5NastranException(Exception):
    pass
