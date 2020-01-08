# pylint: disable=R0902,R0904,R0914,C0111
"""
All table cards are defined in this file.  This includes:

* table_d
 * TABLED1 - Dynamic Table = f(Time, Frequency)
 * TABLED2
 * TABLED3
* table_m
 * TABLEM1 - Material table = f(Temperature)
 * TABLEM2
 * TABLEM3
 * TABLEM4
*tables
 * TABLEST - Material table = f(Stress)
 * TABLES1
 * TABLEHT - Material table = f(Temperature)
 * TABLEH1
*random_tables
 * TABRND1
 * TABRNDG

"""
from typing import Any
import numpy as np

from pyNastran.bdf.field_writer_8 import set_blank_if_default, print_card_8
from pyNastran.bdf.field_writer_16 import print_card_16
from pyNastran.bdf.field_writer_double import print_card_double

from pyNastran.bdf.cards.base_card import BaseCard
from pyNastran.bdf.bdf_interface.assign_type import (
    integer, integer_or_blank, double, string, string_or_blank,
    double_or_string, double_or_blank, integer_or_string)

def make_xy(table_id, table_type, xy):
    try:
        xy = np.array(xy, dtype='float64')
    except ValueError:
        msg = 'cannot parse %s table_id=%r\n' % (table_type, table_id)
        for xi, yi in xy:
            try:
                xi2 = float(xi)
            except ValueError:
                xi2 = '*' + xi
            try:
                yi2 = float(yi)
            except ValueError:
                yi2 = '*' + yi

            msg += '  %s %s\n' % (xi2, yi2)
        raise ValueError(msg)
    x = xy[:, 0]
    y = xy[:, 1]
    return x, y

class Table(BaseCard):
    def __init__(self):
        BaseCard.__init__(self)

    #def parse_fields(self, xy, nrepeated, is_data=False):
        #self.table = TableObj(xy, nrepeated, is_data)

    def write_card(self, size: int=8, is_double: bool=False) -> str:
        card = self.repr_fields()
        if size == 8:
            return self.comment + print_card_8(card)
        if is_double:
            return self.comment + print_card_double(card)
        return self.comment + print_card_16(card)

    #cxy = np.array(self.tc.table.table)
    #fc = cxy[:, 0]
    #yc = cxy[:, 1]
    #assert fc.shape == yc.shape, 'fc.shape=%s yc.shape=%s' % (str(fc.shape), str(yc.shape))
    #print('fc =', fc)
    #print('yc =', yc)
    #self.tc.interpolate(freq)
    #c = interp1d(fc, yc, freq)


class DTABLE(BaseCard):
    type = 'DTABLE'

    @classmethod
    def _init_from_empty(cls):
        default_values = {'CAT' : 1}
        return DTABLE(default_values, comment='')

    def _finalize_hdf5(self, encoding):
        """hdf5 helper function"""
        keys, values = self.default_values
        self.default_values = {key : value if not np.isnan(value) else None
                               for key, value in zip(keys, values)}

    def __init__(self, default_values, comment=''):
        """
        Creates a DTABLE card

        Parameters
        ----------
        default_values : dict
            key : str
                the parameter name
            value : float
                the value
        comment : str; default=''
            a comment for the card

        """
        BaseCard.__init__(self)
        if comment:
            self.comment = comment
        self.default_values = default_values
        #print('default_values = %s' % default_values)
        #for key, value in self.default_values.items():
            #print(key, type(key))
        assert len(self.default_values) > 0, self.default_values
        #print(self)

    @classmethod
    def add_card(cls, card, comment=''):
        """
        Adds a DTABLE card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card

        """
        nfields = len(card) - 1
        assert nfields % 2 == 0, nfields

        default_values = {}
        j = 1
        for i in range(1, nfields + 1, 2):
            label = string(card, i, 'label_%i' % j)
            value = double(card, i + 1, 'value_%i' % j)
            assert label not in default_values, 'label_%i=%r is not unique' % (j, label)
            default_values[label] = value
            j += 1
        assert j >= 2, j
        return DTABLE(default_values, comment=comment)

    def __getitem__(self, key):
        try:
            item = self.default_values[key]
        except KeyError:
            msg = 'expected_key=%r\n' % str(key)
            for keyi, value in self.default_values.items():
                msg += 'DTABLE; key=%r value=%r\n' % (keyi, value)
            raise KeyError(msg)
        return item

    def raw_fields(self):
        list_fields = ['DTABLE']
        #print('***default_values = %s' % self.default_values)
        assert len(self.default_values) > 0, self.default_values
        for label, value in sorted(self.default_values.items()):
            list_fields += [label, value]
        return list_fields

    #def repr_fields(self):
        #return self.raw_fields()

    def write_card(self, size: int=8, is_double: bool=False) -> str:
        card = self.repr_fields()
        if size == 8:
            return self.comment + print_card_8(card)
        if is_double:
            return self.comment + print_card_double(card)
        return self.comment + print_card_16(card)


class TABLED1(Table):
    """
    Dynamic Load Tabular Function, Form 1
    Defines a tabular function for use in generating frequency-dependent and
    time-dependent dynamic loads.

    +---------+------+-------+-------+--------+-----+-----+------+------+
    |    1    |   2  |   3   |   4   |    5   |  6  |  7  |  8   |   9  |
    +=========+======+=======+=======+========+=====+=====+======+======+
    | TABLED1 |  TID | XAXIS | YAXIS | EXTRAP |     |     |      |      |
    +---------+------+-------+-------+--------+-----+-----+------+------+
    |         |  x1  |  y1   |   x2  |   y2   | x3  | y3  | etc. | ENDT |
    +---------+------+-------+-------+--------+-----+-----+------+------+
    | TABLED1 |  32  |       |       |        |     |     |      |      |
    +---------+------+-------+-------+--------+-----+-----+------+------+
    |         | -3.0 |  6.9  |  2.0  |   5.6  | 3.0 | 5.6 | ENDT |      |
    +---------+------+-------+-------+--------+-----+-----+------+------+

    ..note:: EXTRAP is NX specific

    """
    type = 'TABLED1'

    @classmethod
    def _init_from_empty(cls):
        tid = 1
        x = [0., 1.]
        y = [0., 1.]
        return TABLED1(tid, x, y, xaxis='LINEAR', yaxis='LINEAR', extrap=0, comment='')

    def __init__(self, tid, x, y, xaxis='LINEAR', yaxis='LINEAR', extrap=0, comment=''):
        """
        Creates a TABLED1, which is a dynamic load card that is applied
        by the DAREA card

        Parameters
        ----------
        tid : int
            table id
        x : List[float]
            nvalues
        y : List[float]
            nvalues
        xaxis : str
            LINEAR, LOG
        yaxis : str
            LINEAR, LOG
        extrap : int; default=0
            Extrapolation method:
                0 : linear
                1 : constant
            .. note:: this is NX specific
        comment : str; default=''
            a comment for the card

        """
        Table.__init__(self)
        if comment:
            self.comment = comment
        self.tid = tid
        self.extrap = extrap
        self.x = np.asarray(x, dtype='float64')
        self.y = np.asarray(y, dtype='float64')
        self.xaxis = xaxis
        self.yaxis = yaxis
        assert self.xaxis in ['LINEAR', 'LOG'], 'xaxis=%r' % (self.xaxis)
        assert self.yaxis in ['LINEAR', 'LOG'], 'yaxis=%r' % (self.yaxis)

    @classmethod
    def add_card(cls, card, comment=''):
        """
        Adds a TABLED1 card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card

        """
        table_id = integer(card, 1, 'tid')
        xaxis = string_or_blank(card, 2, 'xaxis', 'LINEAR')
        yaxis = string_or_blank(card, 3, 'yaxis', 'LINEAR')
        extrap = integer_or_blank(card, 4, 'yaxis', 0)

        x, y = read_table(card, table_id, 'TABLED1')
        return TABLED1(table_id, x, y, xaxis=xaxis, yaxis=yaxis, extrap=extrap, comment=comment)

    @classmethod
    def add_op2_data(cls, data, comment=''):
        table_id = data[0]
        xaxis = _map_axis(data[1])
        yaxis = _map_axis(data[2])
        xy = data[3:]
        xy = np.array(xy, dtype='float64')
        xy = xy.reshape(xy.size // 2, 2)
        x = xy[:, 0]
        y = xy[:, 1]
        return TABLED1(table_id, x, y, xaxis=xaxis, yaxis=yaxis, comment=comment)

    def raw_fields(self):
        xy = []
        for xi, yi in zip(self.x, self.y):
            xy.extend([xi, yi])
        list_fields = ['TABLED1', self.tid, self.xaxis, self.yaxis, self.extrap,
                       None, None, None, None] + xy + ['ENDT']
        return list_fields

    def repr_fields(self):
        extrap = set_blank_if_default(self.extrap, 0)
        xy = []
        for xi, yi in zip(self.x, self.y):
            xy.extend([xi, yi])
        list_fields = ['TABLED1', self.tid, self.xaxis, self.yaxis, extrap,
                       None, None, None, None] + xy + ['ENDT']
        return list_fields

    def interpolate(self, x):
        if isinstance(x, float):
            x = [x]
        x = np.asarray(x)
        #nx = x.size
        #ny = self.y.size

        # xj follow xi
        i = np.searchsorted(self.x, x, side='left') - 1
        j = i + 1
        #k = np.where(j == ny)[0]

        # TODO: handle out of range errors
        xi = self.x[i]
        yi = self.y[i]
        try:
            xj = self.x[j]
            yj = self.y[j]
        except IndexError:
            #print('table.x = %s' % self.x)
            #print('table.y = %s' % self.y)
            #print('x = %s' % x)
            #print('yi = %s' % yi)
            return yi

        # TODO: could probably speed this up with log rules
        if self.xaxis == 'LINEAR' and self.yaxis == 'LINEAR':
            dx = xj - xi
            y = (xj - x) / dx * yi + (x - xi) / dx * yj
        elif self.xaxis == 'LOG' and self.yaxis == 'LINEAR':
            dx = np.log(xj / xi)
            y = np.log(xj / x) / dx * yi + np.log(x / xi) / dx * yj
        elif self.xaxis == 'LINEAR' and self.yaxis == 'LOG':
            dx = xj - xi
            lny = (xj - x) / dx * np.log(yi) + (x - xi) / dx * np.log(yj)
            y = np.exp(lny)
        elif self.xaxis == 'LOG' and self.yaxis == 'LOG':
            dx = np.log(xj / xi)
            lny = (xj - x) / dx * np.log(yi) + (x - xi) / dx * np.log(yj)
            y = np.exp(lny)
        else:
            raise NotImplementedError('xaxis=%r yaxis=%r' % (self.xaxis, self.yaxis))
        return y


class TABLED2(Table):
    """
    Dynamic Load Tabular Function, Form 2
    Defines a tabular function for use in generating frequency-dependent and
    time-dependent dynamic loads. Also contains parametric data for use with the
    table.
    """
    type = 'TABLED2'

    @classmethod
    def _init_from_empty(cls):
        tid = 1
        x1 = 1.
        x = [0., 1.]
        y = [0., 1.]
        return TABLED2(tid, x1, x, y, extrap=0, comment='')

    def __init__(self, tid, x1, x, y, extrap=0, comment=''):
        """
        Parameters
        ----------
        tid : int
            table id
        x1 : float
            y = yT(x - x1)
        x : List[float]
            the x values
        y : List[float]
            the y values
        extrap : int; default=0
            Extrapolation method:
                0 : linear
                1 : constant
            .. note:: this is NX specific
        comment : str; default=''
            a comment for the card

        """
        Table.__init__(self)
        if comment:
            self.comment = comment
        self.tid = tid
        self.x1 = x1
        self.extrap = extrap
        self.x = np.asarray(x, dtype='float64')
        self.y = np.asarray(y, dtype='float64')

    @classmethod
    def add_card(cls, card, comment=''):
        """
        Adds a TABLED2 card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card

        """
        table_id = integer(card, 1, 'tid')
        x1 = double(card, 2, 'x1')
        extrap = integer_or_blank(card, 3, 'extrap', default=0)
        x, y = read_table(card, table_id, 'TABLED2')
        return TABLED2(table_id, x1, x, y, extrap=extrap, comment=comment)

    @classmethod
    def add_op2_data(cls, data, comment=''):
        table_id = data[0]
        x1 = data[1]
        xy = data[2:]
        xy = np.array(xy, dtype='float64')
        xy = xy.reshape(xy.size // 2, 2)
        x = xy[:, 0]
        y = xy[:, 1]
        return TABLED2(table_id, x1, x, y, comment=comment)

    def interpolate(self, x):
        if isinstance(x, float):
            x = [x]
        x = np.asarray(x)
        #nx = x.size
        #ny = self.y.size

        # xj follow xi
        i = np.searchsorted(self.x, x, side='left') - 1
        j = i + 1
        #k = np.where(j == ny)[0]

        # TODO: handle out of range errors
        xi = self.x[i]
        yi = self.y[i]
        try:
            xj = self.x[j]
            yj = self.y[j]
        except IndexError:
            #print('table.x = %s' % self.x)
            #print('table.y = %s' % self.y)
            #print('x = %s' % x)
            #print('yi = %s' % yi)
            return yi

        dx = xj - xi
        y = (xj - x) / dx * yi + (x - xi) / dx * yj
        return y

    def raw_fields(self):
        xy = []
        for xi, yi in zip(self.x, self.y):
            xy.extend([xi, yi])
        list_fields = ['TABLED2', self.tid, self.x1, self.extrap, None, None,
                       None, None, None] + xy + ['ENDT']
        return list_fields

    def repr_fields(self):
        xy = []
        for xi, yi in zip(self.x, self.y):
            xy.extend([xi, yi])
        extrap = set_blank_if_default(self.extrap, 0)
        list_fields = ['TABLED2', self.tid, self.x1, extrap, None, None,
                       None, None, None] + xy + ['ENDT']
        return list_fields


class TABLED3(Table):
    """
    Dynamic Load Tabular Function, Form 3
    Defines a tabular function for use in generating frequency-dependent and
    time-dependent dynamic loads. Also contains parametric data for use with the
    table.
    """
    type = 'TABLED3'

    @classmethod
    def _init_from_empty(cls):
        tid = 1
        x1 = 1.
        x2 = 2.
        x = [0., 1.]
        y = [0., 1.]
        return TABLED3(tid, x1, x2, x, y, extrap=0, comment='')

    def __init__(self, tid, x1, x2, x, y, extrap=0, comment=''):
        """
        Parameters
        ----------
        tid : int
            table id
        x1 : float
            y = yT(x - x1)
        x2 : ???
            ???
        x : List[float]
            the x values
        y : List[float]
            the y values
        extrap : int; default=0
            Extrapolation method:
                0 : linear
                1 : constant
            .. note:: this is NX specific
        comment : str; default=''
            a comment for the card

        """
        Table.__init__(self)
        if comment:
            self.comment = comment
        self.tid = tid
        self.x1 = x1
        self.x2 = x2
        self.extrap = extrap
        self.x = np.asarray(x, dtype='float64')
        self.y = np.asarray(y, dtype='float64')
        assert self.x2 != 0.0

    @classmethod
    def add_card(cls, card, comment=''):
        """
        Adds a TABLED3 card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card

        """
        table_id = integer(card, 1, 'tid')
        x1 = double(card, 2, 'x1')
        x2 = double(card, 3, 'x2')
        extrap = integer_or_blank(card, 4, 'extrap', default=0)
        x, y = read_table(card, table_id, 'TABLED3')
        return TABLED3(table_id, x1, x2, x, y, extrap=extrap, comment=comment)

    @classmethod
    def add_op2_data(cls, data, comment=''):
        table_id = data[0]
        x1 = data[1]
        x2 = data[2]
        xy = data[3:]
        xy = np.array(xy, dtype='float64')
        xy = xy.reshape(xy.size // 2, 2)
        x = xy[:, 0]
        y = xy[:, 1]
        return TABLED3(table_id, x1, x2, x, y, comment=comment)

    def raw_fields(self):
        xy = []
        for xi, yi in zip(self.x, self.y):
            xy.extend([xi, yi])
        list_fields = ['TABLED3', self.tid, self.x1, self.x2, self.extrap,
                       None, None, None, None] + xy + ['ENDT']
        return list_fields

    def repr_fields(self):
        xy = []
        for xi, yi in zip(self.x, self.y):
            xy.extend([xi, yi])
        extrap = set_blank_if_default(self.extrap, 0)
        list_fields = ['TABLED3', self.tid, self.x1, self.x2, extrap,
                       None, None, None, None] + xy + ['ENDT']
        return list_fields


class TABLED4(Table):
    """
    Dynamic Load Tabular Function, Form 4
    Defines the coefficients of a power series for use in generating
    frequency-dependent and time-dependent dynamic loads. Also contains
    parametric data for use with the table.

    """
    type = 'TABLED4'

    @classmethod
    def _init_from_empty(cls):
        tid = 1
        x1 = 1.
        x2 = 1.
        x3 = 1.
        x4 = 1.
        a = [1., 2.]
        return TABLED4(tid, x1, x2, x3, x4, a, comment='')

    def __init__(self, tid, x1, x2, x3, x4, a, comment=''):
        Table.__init__(self)
        if comment:
            self.comment = comment
        self.tid = tid
        self.x1 = x1
        self.x2 = x2
        self.x3 = x3
        self.x4 = x4
        self.a = np.array(a)

        assert self.x2 != 0.0, 'x2=%s\n%s' % (self.x2, str(self))
        assert self.x3 <= self.x4, 'x3=%s x4=%s\n%s' % (self.x3, self.x4, str(self))

    @classmethod
    def add_card(cls, card, comment=''):
        """
        Adds a TABLED4 card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card

        """
        table_id = integer(card, 1, 'tid')
        x1 = double(card, 2, 'x1')
        x2 = double(card, 3, 'x2')
        x3 = double(card, 4, 'x3')
        x4 = double(card, 5, 'x4')

        nfields = len(card) - 1
        nterms = nfields - 9
        if nterms < 0:
            raise SyntaxError('%r card is too short' % cls.type)

        a = []
        j = 0
        for i in range(9, nfields):
            ai = double(card, i, 'a%i' % (j))
            a.append(ai)
            j += 1
        string(card, nfields, 'ENDT')
        return TABLED4(table_id, x1, x2, x3, x4, a, comment=comment)

    @classmethod
    def add_op2_data(cls, data, comment=''):
        table_id = data[0]
        x1 = data[1]
        x2 = data[2]
        x3 = data[3]
        x4 = data[4]
        a = data[5:]
        return TABLED4(table_id, x1, x2, x3, x4, a, comment=comment)

    def raw_fields(self):
        list_fields = ['TABLED4', self.tid, self.x1, self.x2, self.x3, self.x4,
                       None, None, None] + list(self.a) + ['ENDT']
        return list_fields

    def repr_fields(self):
        return self.raw_fields()

    def interpolate(self, x):
        """
        y = sum_{i=0}^N Ai * ((x-x1)/x2))^i
        """
        if isinstance(x, float):
            x = [x]
        x = np.asarray(x)
        nx = x.size

        na = self.a.size
        n = np.arange(0., na)

        x1 = np.ones(nx) * self.x1
        x2 = np.ones(nx) * self.x2

        i = np.where(x < self.x3)[0]
        x[i] = self.x3

        j = np.where(x > self.x4)[0]
        x[j] = self.x4

        #yi = np.zeros(x.shape, dtype=x.dtype)
        yi = self.a * ((x - x1) / x2) ** n
        return yi.sum()

class TABLED5(Table):
    """
    Dynamic Load Tabular Function, Form 5
    Defines a value as a function of two variables for use in generating
    frequency-dependent and time-dependent dynamic loads.

    """
    type = 'TABLED5'
    def __init__(self, tid, xs, table_ids, comment=''):
        Table.__init__(self)
        if comment:
            self.comment = comment
        self.tid = tid
        self.xs = xs
        self.table_ids = table_ids

    @classmethod
    def add_card(cls, card, comment=''):
        """
        Adds a TABLED5 card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card

        """
        table_id = integer(card, 1, 'tid')

        nfields = len(card) - 1
        nterms = nfields - 9
        if nterms < 0:
            raise SyntaxError('%r card is too short' % cls.type)

        nfields = len(card) - 1
        nterms = (nfields - 9) // 2
        if nterms < 0:
            raise SyntaxError('%r card is too short' % cls.type)
        xs = []
        table_ids = []
        for i in range(nterms):
            n = 9 + i * 2
            if card.field(n) == 'ENDT':
                break
            x = double_or_string(card, n, 'x' + str(i + 1))
            table_id = integer_or_string(card, n + 1, 'table_id' + str(i + 1))
            if x == 'SKIP' or table_id == 'SKIP':
                continue
            xs.append(x)
            table_ids.append(table_id)
        string(card, nfields, 'ENDT')
        return TABLED5(table_id, xs, table_ids, comment=comment)

    #@classmethod
    #def add_op2_data(cls, data, comment=''):
        #table_id = data[0]
        #x1 = data[1]
        #x2 = data[2]
        #x3 = data[3]
        #x4 = data[4]
        #a = data[5:]
        #return TABLED4(table_id, x1, x2, x3, x4, a, comment=comment)

    def raw_fields(self):
        x_table = []
        for xi, tablei in zip(self.xs, self.table_ids):
            x_table.extend([xi, tablei])
        list_fields = ['TABLED5', self.tid, None, None, None, None,
                       None, None, None] + x_table + ['ENDT']
        return list_fields

    def repr_fields(self):
        return self.raw_fields()

    #def interpolate(self, x):
        #"""
        #y = sum_{i=0}^N Ai * ((x-x1)/x2))^i
        #"""
        #if isinstance(x, float):
            #x = [x]
        #x = np.asarray(x)
        #nx = x.size

        #na = self.a.size
        #n = np.arange(0., na)

        #x1 = np.ones(nx) * self.x1
        #x2 = np.ones(nx) * self.x2

        #i = np.where(x < self.x3)[0]
        #x[i] = self.x3

        #j = np.where(x > self.x4)[0]
        #x[j] = self.x4

        #yi = self.a * ((x - x1) / x2) ** n
        #return yi.sum()


class TABDMP1(Table):
    type = 'TABDMP1'

    @classmethod
    def _init_from_empty(cls):
        tid = 1
        x = [0., 1.]
        y = [0., 1.]
        return TABDMP1(tid, x, y, Type='G', comment='')

    def __init__(self, tid: int, x: Any, y: Any, Type: str='G', comment: str='') -> None:
        Table.__init__(self)
        if comment:
            self.comment = comment
        self.tid = tid
        self.Type = Type
        self.x = np.asarray(x, dtype='float64')
        self.y = np.asarray(y, dtype='float64')
        assert self.Type in ['G', 'CRIT', 'Q'], 'Type=%r' % self.Type

    @classmethod
    def add_card(cls, card, comment=''):
        """
        Adds a TABDMP1 card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card

        """
        table_id = integer(card, 1, 'tid')
        Type = string_or_blank(card, 2, 'Type', 'G')
        x, y = read_table(card, table_id, 'TABDMP1')
        return TABDMP1(table_id, x, y, Type=Type, comment=comment)

    @classmethod
    def add_op2_data(cls, data, comment=''):
        table_id = data[0]
        unused_x1 = data[1]
        Type = data[2]
        xy = data[5:]
        xy = np.array(xy, dtype='float64')
        x = xy[:, 0]
        y = xy[:, 1]
        return TABDMP1(table_id, x, y, Type=Type, comment=comment)

    def raw_fields(self):
        xy = []
        for xi, yi in zip(self.x, self.y):
            xy.extend([xi, yi])
        list_fields = ['TABDMP1', self.tid, self.Type, None, None, None, None,
                       None, None] + xy + ['ENDT']
        return list_fields

    def repr_fields(self):
        return self.raw_fields()


class TABLEM1(Table):
    """
    MSC
    ===
    +---------+------+-------+-------+--------+-----+-----+------+------+
    |    1    |   2  |   3   |   4   |   5    |  6  |  7  |  8   |   9  |
    +=========+======+=======+=======+========+=====+=====+======+======+
    | TABLEM1 |  TID |       |       |        |     |     |      |      |
    +---------+------+-------+-------+--------+-----+-----+------+------+
    |         |  x1  |  y1   |   x2  |   y2   | x3  | y3  | etc. | ENDT |
    +---------+------+-------+-------+--------+-----+-----+------+------+
    | TABLEM1 |  32  |       |       |        |     |     |      |      |
    +---------+------+-------+-------+--------+-----+-----+------+------+
    |         | -3.0 |  6.9  |  2.0  |  5.6   | 3.0 | 5.6 | ENDT |      |
    +---------+------+-------+-------+--------+-----+-----+------+------+

    NX
    ==
    +---------+------+-------+-------+--------+-----+-----+------+------+
    |    1    |   2  |   3   |   4   |   5    |  6  |  7  |  8   |   9  |
    +=========+======+=======+=======+========+=====+=====+======+======+
    | TABLEM1 |  TID | XAXIS | YAXIS | EXTRAP |     |     |      |      |
    +---------+------+-------+-------+--------+-----+-----+------+------+
    |         |  x1  |  y1   |   x2  |   y2   | x3  | y3  | etc. | ENDT |
    +---------+------+-------+-------+--------+-----+-----+------+------+
    | TABLEM1 |  32  |       |       |        |     |     |      |      |
    +---------+------+-------+-------+--------+-----+-----+------+------+
    |         | -3.0 |  6.9  |  2.0  |  5.6   | 3.0 | 5.6 | ENDT |      |
    +---------+------+-------+-------+--------+-----+-----+------+------+

    """
    type = 'TABLEM1'
    @classmethod
    def _init_from_empty(cls):
        tid = 1
        x = [0., 1.]
        y = [0., 1.]
        return TABLEM1(tid, x, y, xaxis='LINEAR', yaxis='LINEAR', comment='')

    def __init__(self, tid, x, y, xaxis='LINEAR', yaxis='LINEAR', comment=''):
        Table.__init__(self)
        if comment:
            self.comment = comment
        self.tid = tid
        self.xaxis = xaxis # linear/log
        self.yaxis = yaxis # linear/log
        self.x = np.asarray(x, dtype='float64')
        self.y = np.asarray(y, dtype='float64')

    @classmethod
    def add_card(cls, card, comment=''):
        """
        Adds a TABLEM1 card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card

        """
        table_id = integer(card, 1, 'tid')
        xaxis = string_or_blank(card, 2, 'xaxis', 'LINEAR')
        yaxis = string_or_blank(card, 3, 'yaxis', 'LINEAR')
        x, y = read_table(card, table_id, 'TABLEM1')
        return TABLEM1(table_id, x, y, xaxis=xaxis, yaxis=yaxis, comment=comment)

    @classmethod
    def add_op2_data(cls, data, comment=''):
        table_id = data[0]
        xy = data[1:]
        xy = np.array(xy, dtype='float64')
        xy = xy.reshape(xy.size // 2, 2)
        x = xy[:, 0]
        y = xy[:, 1]
        return TABLEM1(table_id, x, y, comment=comment)

    def raw_fields(self):
        xy = []
        for xi, yi in zip(self.x, self.y):
            xy.extend([xi, yi])
        list_fields = ['TABLEM1', self.tid, self.xaxis, self.yaxis, None, None,
                       None, None, None] + xy + ['ENDT']
        return list_fields

    def repr_fields(self):
        xaxis = set_blank_if_default(self.xaxis, 'LINEAR')
        yaxis = set_blank_if_default(self.yaxis, 'LINEAR')
        xy = []
        for xi, yi in zip(self.x, self.y):
            xy.extend([xi, yi])
        list_fields = ['TABLEM1', self.tid, xaxis, yaxis, None, None,
                       None, None, None] + xy + ['ENDT']
        return list_fields


class TABLEM2(Table):
    """
    +---------+------+-------+--------+-----+-----+-----+------+------+
    |    1    |   2  |   3   |   4    |  5  |  6  |  7  |  8   |   9  |
    +=========+======+=======+========+=====+=====+=====+======+======+
    | TABLEM2 |  TID |   X1  | EXTRAP |     |     |     |      |      |
    +---------+------+-------+--------+-----+-----+-----+------+------+
    |         |  x1  |  y1   |   x2   | y2  | x3  | y3  | etc. | ENDT |
    +---------+------+-------+--------+-----+-----+-----+------+------+
    | TABLEM2 |  32  | -10.5 |        |     |     |     |      |      |
    +---------+------+-------+--------+-----+-----+-----+------+------+
    |         | -3.0 |  6.9  |  2.0   | 5.6 | 3.0 | 5.6 | ENDT |      |
    +---------+------+-------+--------+-----+-----+-----+------+------+

    ..note:: EXTRAP is NX specific

    """
    type = 'TABLEM2'

    @classmethod
    def _init_from_empty(cls):
        tid = 1
        x1 = 1.
        x = [0., 1.]
        y = [0., 1.]
        return TABLEM2(tid, x1, x, y, extrap=0, comment='')

    def __init__(self, tid, x1, x, y, extrap=0, comment=''):
        Table.__init__(self)
        if comment:
            self.comment = comment
        self.tid = tid
        self.x1 = x1
        self.extrap = extrap
        self.x = np.asarray(x, dtype='float64')
        self.y = np.asarray(y, dtype='float64')

    @classmethod
    def add_card(cls, card, comment=''):
        """
        Adds a TABLEM2 card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card
        """
        table_id = integer(card, 1, 'tid')

        # defined in MSC as an integer and used as a float...int > 0
        # defined in NX as a float; real
        # no default given in either, but from context, let's assume 0.0
        x1 = double_or_blank(card, 2, 'x1', 0.0)
        extrap = integer_or_blank(card, 3, 'EXTRAP', default=0)
        x, y = read_table(card, table_id, 'TABLEM2')
        return TABLEM2(table_id, x1, x, y, extrap=extrap, comment=comment)

    @classmethod
    def add_op2_data(cls, data, comment=''):
        table_id = data[0]
        x1 = data[1]
        xy = data[2:]
        xy = np.array(xy, dtype='float64')
        xy = xy.reshape(xy.size // 2, 2)
        x = xy[:, 0]
        y = xy[:, 1]
        return TABLEM2(table_id, x1, x, y, comment=comment)

    def raw_fields(self):
        xy = []
        for xi, yi in zip(self.x, self.y):
            xy.extend([xi, yi])
        list_fields = ['TABLEM2', self.tid, self.x1, self.extrap, None, None,
                       None, None, None] + xy + ['ENDT']
        return list_fields

    def repr_fields(self):
        extrap = set_blank_if_default(self.extrap, 0)
        xy = []
        for xi, yi in zip(self.x, self.y):
            xy.extend([xi, yi])
        list_fields = ['TABLEM2', self.tid, self.x1, extrap, None, None,
                       None, None, None] + xy + ['ENDT']
        return list_fields


class TABLEM3(Table):
    """
    +---------+------+-------+-------+--------+-----+-----+------+------+
    |    1    |   2  |   3   |   4   |    5   |  6  |  7  |  8   |   9  |
    +=========+======+=======+=======+========+=====+=====+======+======+
    | TABLEM3 |  TID |   X1  |   X2  | EXTRAP |     |     |      |      |
    +---------+------+-------+-------+--------+-----+-----+------+------+
    |         |  x1  |  y1   |   x2  |   y2   | x3  | y3  | etc. | ENDT |
    +---------+------+-------+-------+--------+-----+-----+------+------+
    | TABLEM3 |  32  | 126.9 | 30.0  |        |     |     |      |      |
    +---------+------+-------+-------+--------+-----+-----+------+------+
    |         | -3.0 |  6.9  |  2.0  |   5.6  | 3.0 | 5.6 | ENDT |      |
    +---------+------+-------+-------+--------+-----+-----+------+------+

    ..note:: EXTRAP is NX specific

    """
    type = 'TABLEM3'

    @classmethod
    def _init_from_empty(cls):
        tid = 1
        x1 = 1.
        x2 = 2.
        x = [0., 1.]
        y = [0., 1.]
        return TABLEM3(tid, x1, x2, x, y, extrap=0, comment='')

    def __init__(self, tid, x1, x2, x, y, extrap=0, comment=''):
        Table.__init__(self)
        if comment:
            self.comment = comment
        self.tid = tid
        self.x1 = x1
        self.x2 = x2
        self.extrap = extrap
        self.x = np.asarray(x, dtype='float64')
        self.y = np.asarray(y, dtype='float64')
        assert self.x2 != 0.0

    @classmethod
    def add_card(cls, card, comment=''):
        """
        Adds a TABLEM3 card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card

        """
        table_id = integer(card, 1, 'tid')
        x1 = double(card, 2, 'x1')
        x2 = double(card, 3, 'x2')
        extrap = integer_or_blank(card, 4, 'extrap', default=0)
        x, y = read_table(card, table_id, 'TABLEM3')
        return TABLEM3(table_id, x1, x2, x, y, extrap=extrap, comment=comment)

    @classmethod
    def add_op2_data(cls, data, comment=''):
        table_id = data[0]
        x1 = data[1]
        x2 = data[2]
        xy = data[3:]
        xy = np.array(xy, dtype='float64')
        xy = xy.reshape(xy.size // 2, 2)
        x = xy[:, 0]
        y = xy[:, 1]
        return TABLEM3(table_id, x1, x2, x, y, comment=comment)

    def raw_fields(self):
        xy = []
        for xi, yi in zip(self.x, self.y):
            xy.extend([xi, yi])
        list_fields = ['TABLEM3', self.tid, self.x1, self.x2, self.extrap,
                       None, None, None, None] + xy + ['ENDT']
        return list_fields

    def repr_fields(self):
        xy = []
        extrap = set_blank_if_default(self.extrap, 0)
        for xi, yi in zip(self.x, self.y):
            xy.extend([xi, yi])
        list_fields = ['TABLEM3', self.tid, self.x1, self.x2, extrap,
                       None, None, None, None] + xy + ['ENDT']
        return list_fields


class TABLEM4(Table):
    """
    +---------+------+---------+--------+-----+--------+------+------+
    |    1    |   2  |    3    |    4   |  5  |   6    |  7   |  8   |
    +=========+======+=========+========+=====+========+======+======+
    | TABLEM4 |  TID |   X1    |   X2   | X3  |   X4   |      |      |
    +---------+------+---------+--------+-----+--------+------+------+
    |         |  A1  |   A2    |   A3   | A4  |   A5   | etc. | ENDT |
    +---------+------+---------+--------+-----+--------+------+------+
    | TABLEM4 |  32  |   0.0   |   1.0  | 0.0 |  100.  |      |      |
    +---------+------+---------+--------+-----+--------+------+------+
    |         | 2.91 | -0.0329 | 6.51-5 | 0.0 | -3.4-7 | ENDT |      |
    +---------+------+---------+--------+-----+--------+------+------+
    """
    type = 'TABLEM4'

    @classmethod
    def _init_from_empty(cls):
        tid = 1
        x1 = 1.
        x2 = 1.
        x3 = 1.
        x4 = 2.
        a = [1., 2.]
        return TABLEM4(tid, x1, x2, x3, x4, a, comment='')

    def __init__(self, tid, x1, x2, x3, x4, a, comment=''):
        Table.__init__(self)
        if comment:
            self.comment = comment
        self.tid = tid
        self.x1 = x1
        self.x2 = x2
        self.x3 = x3
        self.x4 = x4
        self.a = np.asarray(a)
        assert self.x2 != 0.0, 'x2=%s\n%s' % (self.x2, str(self))
        assert self.x3 <= self.x4, 'x3=%s x4=%s\n%s' % (self.x3, self.x4, str(self))

    @classmethod
    def add_card(cls, card, comment=''):
        """
        Adds a TABLEM4 card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card

        """
        table_id = integer(card, 1, 'tid')
        x1 = double(card, 2, 'x1')
        x2 = double(card, 3, 'x2')
        x3 = double(card, 4, 'x3')
        x4 = double(card, 5, 'x4')

        nfields = len(card) - 1
        nterms = nfields - 9
        if nterms < 0:
            raise SyntaxError('%r card is too short' % cls.type)

        a = []
        j = 0
        for i in range(9, nfields):
            ai = double_or_blank(card, i, 'a%i' % (j), 0.0)
            a.append(ai)
            j += 1
        string(card, nfields, 'ENDT')
        return TABLEM4(table_id, x1, x2, x3, x4, a, comment=comment)

    @classmethod
    def add_op2_data(cls, data, comment=''):
        """
        Adds a TABLEM4 card from the OP2

        Parameters
        ----------
        data : List[varies]
            a list of fields defined in OP2 format
        comment : str; default=''
            a comment for the card

        """
        table_id = data[0]
        x1 = data[1]
        x2 = data[2]
        x3 = data[3]
        x4 = data[4]
        a = data[3:]
        return TABLEM4(table_id, x1, x2, x3, x4, a, comment=comment)

    def raw_fields(self):
        list_fields = ['TABLEM4', self.tid, self.x1, self.x2, self.x3, self.x4,
                       None, None, None] + list(self.a) + ['ENDT']
        return list_fields

    def repr_fields(self):
        return self.raw_fields()


class TABLES1(Table):
    """
    +---------+------+-------+-------+--------+-----+-------+------+------+
    |    1    |   2  |   3   |   4   |    5   |  6  |   7   |  8   |   9  |
    +=========+======+=======+=======+========+=====+=======+======+======+
    | TABLES1 |  TID | TYPE  |       |        |     |       |      |      |
    +---------+------+-------+-------+--------+-----+-------+------+------+
    |         |  x1  |  y1   |   x2  |   y2   | x3  |  y3   | etc. | ENDT |
    +---------+------+-------+-------+--------+-----+-------+------+------+
    | TABLES1 |  32  |       |       |        |     |       |      |      |
    +---------+------+-------+-------+--------+-----+-------+------+------+
    |         |  0.0 |  0.0  |  0.01 |  1000. | 0.2 | 1500. | ENDT |      |
    +---------+------+-------+-------+--------+-----+-------+------+------+
    """
    type = 'TABLES1'

    @classmethod
    def _init_from_empty(cls):
        tid = 1
        x = [0., 1.]
        y = [0., 1.]
        return TABLES1(tid, x, y, Type=1, comment='')

    def __init__(self, tid, x, y, Type=1, comment=''):
        """
        Adds a TABLES1 card, which defines a stress dependent material

        Parameters
        ----------
        tid : int
            Table ID
        Type : int; default=1
            Type of stress-strain curve (1 or 2)
            1 - Cauchy (true) stress vs. total true strain
            2 - Cauchy (true) stress vs. plastic true strain (MSC only)
            Type is MSC-specific and was added somewhere between
            2006 and 2016.
        x, y : List[float]
            table values
        comment : str; default=''
            a comment for the card

        """
        Table.__init__(self)
        if comment:
            self.comment = comment
        self.tid = tid
        self.Type = Type
        self.x = np.asarray(x, dtype='float64')
        self.y = np.asarray(y, dtype='float64')
        assert self.Type in [1, 2], 'TABLES1 Type=%s' % self.Type

    @classmethod
    def add_card(cls, card, comment=''):
        """
        Adds a TABLES1 card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card

        """
        table_id = integer(card, 1, 'tid')
        Type = integer_or_blank(card, 2, 'Type', 1)
        x, y = read_table(card, table_id, 'TABLES1')
        return TABLES1(table_id, x, y, Type=Type, comment=comment)

    @classmethod
    def add_op2_data(cls, data, comment=''):
        """
        Adds a TABLES1 card from the OP2

        Parameters
        ----------
        data : List[varies]
            a list of fields defined in OP2 format
        comment : str; default=''
            a comment for the card

        """
        table_id = data[0]
        xy = data[1:]
        xy = np.array(xy, dtype='float64')
        xy = xy.reshape(xy.size // 2, 2)
        x = xy[:, 0]
        y = xy[:, 1]
        return TABLES1(table_id, x, y, Type=1, comment=comment)

    def raw_fields(self):
        xy = []
        for xi, yi in zip(self.x, self.y):
            xy.extend([xi, yi])
        list_fields = ['TABLES1', self.tid, self.Type, None, None, None,
                       None, None, None] + xy + ['ENDT']
        return list_fields

    def repr_fields(self):
        xy = []
        for xi, yi in zip(self.x, self.y):
            xy.extend([xi, yi])

        # MSC 2005.2 doesn't support Type; 2016.1 does
        stress_strain_curve_type = set_blank_if_default(self.Type, 1)
        list_fields = ['TABLES1', self.tid, stress_strain_curve_type, None, None, None,
                       None, None, None] + xy + ['ENDT']
        return list_fields


class TABLEST(Table):
    """
    +---------+-------+-------+-------+--------+------+------+------+------+
    |    1    |   2   |   3   |   4   |    5   |  6  |   7   |  8   |   9  |
    +=========+=======+=======+=======+========+=====+=======+======+======+
    | TABLEST |  TID  |       |       |        |      |      |      |      |
    +---------+-------+-------+-------+--------+------+------+------+------+
    |         |   x1  |  y1   |   x2  |   y2   |  x3  |  y3  | etc. | ENDT |
    +---------+-------+-------+-------+--------+------+------+------+------+
    | TABLEST |   32  |       |       |        |      |      |      |      |
    +---------+-------+-------+-------+--------+------+------+------+------+
    |         | 150.0 |  10.0 | 175.0 |  20.   | ENDT |      |      |      |
    +---------+-------+-------+-------+--------+------+------+------+------+
    """
    type = 'TABLEST'

    def __init__(self, tid, x, y, comment=''):
        Table.__init__(self)
        if comment:
            self.comment = comment
        self.tid = tid
        self.x = np.asarray(x, dtype='float64')
        self.y = np.asarray(y, dtype='float64')

    @classmethod
    def add_card(cls, card, comment=''):
        """
        Adds a TABLEST card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card

        """
        table_id = integer(card, 1, 'tid')
        x, y = read_table(card, table_id, 'TABLEST')
        return TABLEST(table_id, x, y, comment=comment)

    @classmethod
    def add_op2_data(cls, data, comment=''):
        """
        Adds a TABLEST card from the OP2

        Parameters
        ----------
        data : List[varies]
            a list of fields defined in OP2 format
        comment : str; default=''
            a comment for the card

        """
        table_id = data[0]
        xy = data[1:]
        xy = np.array(xy, dtype='float64')
        xy = xy.reshape(xy.size // 2, 2)
        x = xy[:, 0]
        y = xy[:, 1]
        return TABLEST(table_id, x, y, comment=comment)

    def raw_fields(self):
        xy = []
        for xi, yi in zip(self.x, self.y):
            xy.extend([xi, yi])
        list_fields = ['TABLEST', self.tid, None, None, None, None,
                       None, None, None] + xy + ['ENDT']
        return list_fields

    def repr_fields(self):
        return self.raw_fields()


class TABLEH1(Table):
    """
    +---------+------+-------+-------+--------+-----+-------+------+------+
    |    1    |   2  |   3   |   4   |    5   |  6  |   7   |  8   |   9  |
    +=========+======+=======+=======+========+=====+=======+======+======+
    | TABLEH1 |  TID |       |       |        |     |       |      |      |
    +---------+------+-------+-------+--------+-----+-------+------+------+
    |         |  x1  |  y1   |   x2  |   y2   | x3  |  y3   | etc. | ENDT |
    +---------+------+-------+-------+--------+-----+-------+------+------+
    | TABLEH1 |  32  |       |       |        |     |       |      |      |
    +---------+------+-------+-------+--------+-----+-------+------+------+
    |         |  0.0 |  0.0  |  0.01 |  1000. | 0.2 | 1500. | ENDT |      |
    +---------+------+-------+-------+--------+-----+-------+------+------+
    """
    type = 'TABLEH1'

    @classmethod
    def _init_from_empty(cls):
        tid = 1
        x = [0., 1.]
        y = [0., 1.]
        return TABLEH1(tid, x, y, comment='')

    def __init__(self, tid, x, y, comment=''):
        """
        Adds a TABLEH1 card, which defines convection heat transfer coefficient.
        It's referenced by a TABLEHT.

        Parameters
        ----------
        tid : int
            Table ID
        x, y : List[float]
            table values
        comment : str; default=''
            a comment for the card

        """
        Table.__init__(self)
        if comment:
            self.comment = comment
        self.tid = tid
        self.x = np.asarray(x, dtype='float64')
        self.y = np.asarray(y, dtype='float64')

    @classmethod
    def add_card(cls, card, comment=''):
        """
        Adds a TABLEH1 card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card

        """
        table_id = integer(card, 1, 'tid')
        x, y = read_table(card, table_id, 'TABLEH1')
        return TABLEH1(table_id, x, y, comment=comment)

    @classmethod
    def add_op2_data(cls, data, comment=''):
        """
        Adds a TABLEH1 card from the OP2

        Parameters
        ----------
        data : List[varies]
            a list of fields defined in OP2 format
        comment : str; default=''
            a comment for the card

        """
        table_id = data[0]
        xy = data[1:]
        xy = np.array(xy, dtype='float64')
        xy = xy.reshape(xy.size // 2, 2)
        x = xy[:, 0]
        y = xy[:, 1]
        return TABLEH1(table_id, x, y, comment=comment)

    def raw_fields(self):
        xy = []
        for xi, yi in zip(self.x, self.y):
            xy.extend([xi, yi])
        list_fields = ['TABLEH1', self.tid, None, None, None, None,
                       None, None, None] + xy + ['ENDT']
        return list_fields

    def repr_fields(self):
        xy = []
        for xi, yi in zip(self.x, self.y):
            xy.extend([xi, yi])

        list_fields = ['TABLEH1', self.tid, None, None, None, None,
                       None, None, None] + xy + ['ENDT']
        return list_fields

class TABLEHT(Table):
    """
    +---------+-------+-------+-------+--------+------+------+------+------+
    |    1    |   2   |   3   |   4   |    5   |  6  |   7   |  8   |   9  |
    +=========+=======+=======+=======+========+=====+=======+======+======+
    | TABLEHT |  TID  |       |       |        |      |      |      |      |
    +---------+-------+-------+-------+--------+------+------+------+------+
    |         |   x1  |  tid1 |   x2  |  tid2  |  x3  | tid3 | etc. | ENDT |
    +---------+-------+-------+-------+--------+------+------+------+------+
    | TABLEHT |   32  |       |       |        |      |      |      |      |
    +---------+-------+-------+-------+--------+------+------+------+------+
    |         |  1.   |   10  |  5.   |   11   | ENDT |      |      |      |
    +---------+-------+-------+-------+--------+------+------+------+------+
    """
    type = 'TABLEHT'

    @classmethod
    def _init_from_empty(cls):
        tid = 1
        x = [0., 1.]
        y = [0., 1.]
        return TABLEHT(tid, x, y, comment='')

    def __init__(self, tid: int, x, y, comment=''):
        """
        Adds a TABLEHT card, which a function of two variables for
        convection heat transfer coefficient.

        Parameters
        ----------
        tid : int
            Table ID
        x, y : List[float]
            table values
        comment : str; default=''
            a comment for the card

        """
        Table.__init__(self)
        if comment:
            self.comment = comment
        self.tid = tid
        self.x = np.asarray(x, dtype='float64')
        self.y = np.asarray(y, dtype='int32')

    @classmethod
    def add_card(cls, card, comment=''):
        """
        Adds a TABLEHT card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card

        """
        table_id = integer(card, 1, 'tid')
        x, y = read_table_float_int(card, table_id, 'TABLEHT')
        return TABLEHT(table_id, x, y, comment=comment)

    @classmethod
    def add_op2_data(cls, data, comment=''):
        """
        Adds a TABLEHT card from the OP2

        Parameters
        ----------
        data : List[varies]
            a list of fields defined in OP2 format
        comment : str; default=''
            a comment for the card

        """
        table_id = data[0]
        xy = data[1:]
        xy = np.array(xy, dtype='float64')
        xy = xy.reshape(xy.size // 2, 2)
        x = xy[:, 0]
        y = xy[:, 1]
        return TABLEHT(table_id, x, y, comment=comment)

    def raw_fields(self):
        xy = []
        for xi, yi in zip(self.x, self.y):
            xy.extend([xi, yi])
        list_fields = ['TABLEHT', self.tid, None, None, None, None,
                       None, None, None] + xy + ['ENDT']
        return list_fields

    def repr_fields(self):
        return self.raw_fields()


#class RandomTable(Table):
    #type = 'TABLE??'

    #def __init__(self):
        #Table.__init__(self)


class TABRND1(Table):
    type = 'TABRND1'

    @classmethod
    def _init_from_empty(cls):
        tid = 1
        x = [0., 1.]
        y = [0., 1.]
        return TABRND1(tid, x, y, xaxis='LINEAR', yaxis='LINEAR', comment='')

    def __init__(self, tid, x, y, xaxis='LINEAR', yaxis='LINEAR', comment=''):
        Table.__init__(self)
        if comment:
            self.comment = comment

        self.tid = tid
        self.x = np.asarray(x, dtype='float64')
        self.y = np.asarray(y, dtype='float64')
        self.xaxis = xaxis
        self.yaxis = yaxis
        assert self.xaxis in ['LINEAR', 'LOG'], 'xaxis=%r' % (self.xaxis)
        assert self.yaxis in ['LINEAR', 'LOG'], 'yaxis=%r' % (self.yaxis)

    @classmethod
    def add_card(cls, card, comment=''):
        """
        Adds a TABRND1 card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card

        """
        table_id = integer(card, 1, 'tid')
        xaxis = string_or_blank(card, 2, 'xaxis', 'LINEAR')
        yaxis = string_or_blank(card, 3, 'yaxis', 'LINEAR')
        x, y = read_table(card, table_id, 'TABRND1')
        return TABRND1(table_id, x, y, xaxis=xaxis, yaxis=yaxis, comment=comment)

    @classmethod
    def add_op2_data(cls, data, comment=''):
        """
        Adds a TABRND1 card from the OP2

        Parameters
        ----------
        data : List[varies]
            a list of fields defined in OP2 format
        comment : str; default=''
            a comment for the card

        """
        table_id = data[0]
        xaxis = _map_axis(data[1])
        yaxis = _map_axis(data[2])
        xy = data[3:]
        xy = np.array(xy, dtype='float64')
        x = xy[:, 0]
        y = xy[:, 1]
        return TABRND1(table_id, x, y, xaxis=xaxis, yaxis=yaxis, comment=comment)

    #def parse_fields(self, xy, nrepeated, is_data=False):
        #self.table = TableObj(xy, nrepeated, is_data)

    def raw_fields(self):
        xy = []
        for xi, yi in zip(self.x, self.y):
            xy.extend([xi, yi])
        list_fields = ['TABRND1', self.tid, self.xaxis, self.yaxis, None, None,
                       None, None, None] + xy + ['ENDT']
        return list_fields

    def repr_fields(self):
        xy = []
        for xi, yi in zip(self.x, self.y):
            xy.extend([xi, yi])
        xaxis = set_blank_if_default(self.xaxis, 'LINEAR')
        yaxis = set_blank_if_default(self.yaxis, 'LINEAR')
        list_fields = ['TABRND1', self.tid, xaxis, yaxis, None, None,
                       None, None, None] + xy + ['ENDT']
        return list_fields


class TABRNDG(Table):
    r"""
    Gust Power Spectral Density

    Defines the power spectral density (PSD) of a gust for aeroelastic response
    analysis.

    """
    type = 'TABRNDG'

    @classmethod
    def _init_from_empty(cls):
        tid = 1
        Type = 1
        LU = 1.
        WG = 1.
        return TABRNDG(tid, Type, LU, WG, comment='')

    def __init__(self, tid, Type, LU, WG, comment=''):
        """
        Creates a TABRNDG card

        Parameters
        ----------
        tid : int
            table id
        Type : int
           PSD type
           1 : von Karman
           2 : Dryden
        LU : float
            Scale of turbulence divided by velocity (units of time)
        WG : float
            Root-mean-square gust velocity
        comment : str; default=''
            a comment for the card

        """
        Table.__init__(self)
        if comment:
            self.comment = comment
        #: Table identification number. (Integer >0)
        self.tid = tid
        #: PSD Type: 1. von Karman; 2. Dryden
        self.Type = Type
        #: Scale of turbulence divided by velocity (units of time; Real)
        self.LU = LU
        #: Root-mean-square gust velocity. (Real)
        self.WG = WG
        assert self.Type in [1, 2], ('Type must be 1 or 2.  '
                                     'Type=%s' % (self.Type))

    @classmethod
    def add_card(cls, card, comment=''):
        """
        Adds a TABRNDG card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card

        """
        table_id = integer(card, 1, 'tid')
        Type = integer(card, 2, 'Type')
        LU = double(card, 3, 'LU')
        WG = double(card, 4, 'WG')
        return TABRNDG(table_id, Type, LU, WG, comment=comment)

    def raw_fields(self):
        list_fields = ['TABRNDG', self.tid, self.Type, self.LU, self.WG]
        return list_fields

    def repr_fields(self):
        return self.raw_fields()


def _map_axis(axis):
    if axis == 0:
        axis_type = 'LINEAR'
    elif axis == 1:
        axis_type = 'LOG'
    else: # pragma: no cover
        raise ValueError('axis=%r' % axis)
    return axis_type

def read_table(card, table_id, table_type):
    """common method for reading tables that handles SKIP"""
    nfields = len(card) - 1
    nterms = (nfields - 9) // 2
    if nterms < 0:
        raise SyntaxError('%r card is too short' % table_type)

    xy = []
    for i in range(nterms):
        n = 9 + i * 2
        if card.field(n) == 'ENDT':
            break
        xi = double_or_string(card, n, 'x' + str(i + 1))
        yi = double_or_string(card, n + 1, 'y' + str(i + 1))
        if xi == 'SKIP' or yi == 'SKIP':
            continue
        xy.append([xi, yi])
    string(card, nfields, 'ENDT')
    x, y = make_xy(table_id, table_type, xy)
    return x, y

def read_table_float_int(card, table_id, table_type):
    """common method for reading tables that handles SKIP"""
    nfields = len(card) - 1
    nterms = (nfields - 9) // 2
    if nterms < 0:
        raise SyntaxError('%r card is too short' % table_type)

    xy = []
    for i in range(nterms):
        n = 9 + i * 2
        if card.field(n) == 'ENDT':
            break
        xi = double_or_string(card, n, 'x' + str(i + 1))
        yi = integer_or_string(card, n + 1, 'y' + str(i + 1))
        if xi == 'SKIP' or yi == 'SKIP':
            continue
        xy.append([xi, yi])
    string(card, nfields, 'ENDT')
    x, y = make_xy(table_id, table_type, xy)
    return x, y
