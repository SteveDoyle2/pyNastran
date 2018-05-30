"""
Default and unknown values, for bdf cards in particular.
If a value is blank or None coming from pyNastran, use the defaults.
If the h5Nastran developer doesn't know what the value should be, use the unknowns.
"""
import numpy as np
import tables


class Defaults(object):
    default_double = np.nan
    default_int = -10
    default_str = ''

    unknown_double = -10101010.
    unknown_int = -20
    unknown_str = '?'

    dtype = np.dtype([
        ('DEFAULT_DOUBLE', '<f8'),
        ('DEFAULT_INT', '<i8'),
        ('DEFAULT_STR', 'S8'),
        ('UNKNOWN_DOUBLE', '<f8'),
        ('UNKNOWN_INT', '<i8'),
        ('UNKNOWN_STR', 'S8')
    ])

    Format = tables.descr_from_dtype(dtype)[0]

    def get_value_int(self, val):
        if val is None:
            return self.default_int
        else:
            return val

    def get_value_double(self, val):
        if val is None:
            return self.default_double
        else:
            return val

    def get_value_str(self, val):
        if val is None:
            return self.default_str
        else:
            return val

    def get_list_double(self, vals):
        if vals[0] is None:
            return [self.default_double] * len(vals)
        else:
            return vals

    def get_list_int(self, vals):
        if vals[0] is None:
            return [self.default_int] * len(vals)
        else:
            return vals

    def get_list_str(self, vals):
        if vals[0] is None:
            return [self.default_str] * len(vals)
        else:
            return vals

    def to_value_int(self, val):
        if val == self.default_int:
            return None
        else:
            return val

    def to_value_double(self, val):
        if val == self.default_double:
            return None
        else:
            return val

    def to_value_str(self, val):
        if val == self.default_str:
            return None
        else:
            return val

    def to_list_double(self, vals):
        if vals[0] == self.default_double:
            return None
        else:
            return vals

    def to_list_int(self, vals):
        if vals[0] == self.default_int:
            return None
        else:
            return vals

    def to_list_str(self, vals):
        if vals[0] == self.default_str:
            return None
        else:
            return vals

    def save(self, h5n):       
        h5f = h5n.h5f
        
        h5f.create_table(h5n.table_paths.defaults_path, h5n.table_paths.defaults_table, self.Format,
                              'DATA DEFAULTS', expectedrows=1, createparents=True)

        table = h5f.get_node(h5n.table_paths.defaults)
        table.append(self._save())

    def load(self, h5n):
        self._load(h5n.h5f.get_node(h5n.table_paths.defaults).read())
        
    def _load(self, data):
        self.default_double = data['DEFAULT_DOUBLE'][0]
        self.default_int = data['DEFAULT_INT'][0]
        self.default_str = data['DEFAULT_STR'][0].decode()
        self.unknown_double = data['UNKNOWN_DOUBLE'][0]
        self.unknown_int = data['UNKNOWN_INT'][0]
        self.unknown_str = data['UNKNOWN_STR'][0].decode()
        
    def _save(self):
        data = np.empty(1, dtype=self.dtype)
        data['DEFAULT_DOUBLE'][0] = self.default_double
        data['DEFAULT_INT'][0] = self.default_int
        data['DEFAULT_STR'][0] = self.default_str
        data['UNKNOWN_DOUBLE'][0] = self.unknown_double
        data['UNKNOWN_INT'][0] = self.unknown_int
        data['UNKNOWN_STR'][0] = self.unknown_str
        return data
