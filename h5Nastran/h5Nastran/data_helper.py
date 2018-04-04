"""
Default and unknown values, for bdf cards in particular.
If a value is blank or None coming from pyNastran, use the defaults.
If the h5Nastran developer doesn't know what the value should be, use the unknowns.
"""
import numpy as np
import tables


class DataHelper(object):
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

    def save(self):
        data = np.empty(1, dtype=self.dtype)
        data['DEFAULT_DOUBLE'][0] = self.default_double
        data['DEFAULT_INT'][0] = self.default_int
        data['DEFAULT_STR'][0] = self.default_str
        data['UNKNOWN_DOUBLE'][0] = self.unknown_double
        data['UNKNOWN_INT'][0] = self.unknown_int
        data['UNKNOWN_STR'][0] = self.unknown_str
        return data

    def load(self, data):
        self.default_double = data['DEFAULT_DOUBLE'][0]
        self.default_int = data['DEFAULT_INT'][0]
        self.default_str = data['DEFAULT_STR'][0].decode()
        self.unknown_double = data['UNKNOWN_DOUBLE'][0]
        self.unknown_int = data['UNKNOWN_INT'][0]
        self.unknown_str = data['UNKNOWN_STR'][0].decode()
