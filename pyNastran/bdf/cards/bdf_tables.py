# pylint: disable=R0902,R0904,R0914,C0111
"""
All table cards are defined in this file.  This includes:

* Table
 * TABLED1 - Dynamic Table = f(Time, Frequency)
 * TABLED2
 * TABLED3
 * TABLEM1 - Material table = f(Temperature)
 * TABLEM2
 * TABLEM3
 * TABLEM4
 * TABLES1 - Material table = f(Stress)
 * TABLEST
 * RandomTable
   * TABRND1
 * TABRNDG
 * TIC

All tables have a self.table parameter that is a TableObj
"""
from __future__ import (nested_scopes, generators, division, absolute_import,
                        print_function, unicode_literals)
from six import iteritems
from six.moves import range

import numpy as np

from pyNastran.bdf.field_writer_8 import set_blank_if_default, print_card_8
from pyNastran.bdf.field_writer_16 import print_card_16
from pyNastran.bdf.field_writer_double import print_card_double

from pyNastran.bdf.cards.base_card import BaseCard
#from pyNastran.utils.dev import list_print
from pyNastran.bdf.bdf_interface.assign_type import (integer,
    integer_or_blank, double, components, string, string_or_blank,
    double_or_string)

class Table(BaseCard):
    def __init__(self):
        pass

    def _map_axis(self, axis):
        if axis == 0:
            axis_type = 'LINEAR'
        else:
            raise ValueError('axis=%r' % axis)
        return axis_type

    #def parse_fields(self, xy, nrepeated, is_data=False):
        #self.table = TableObj(xy, nrepeated, is_data)

    def write_card(self, size=8, is_double=False):
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
    def __init__(self, default_values, comment=''):
        if comment:
            self._comment = comment
        self.default_values = default_values
        print('default_values = %s' % default_values)
        #for key, value in iteritems(self.default_values):
            #print(key, type(key))
        assert len(self.default_values) > 0, self.default_values
        print(self)

    @classmethod
    def add_card(cls, card, comment=''):
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
            for keyi, value in iteritems(self.default_values):
                msg += 'DTABLE; key=%r value=%r\n' % (keyi, value)
            raise KeyError(msg)
        return item

    def raw_fields(self):
        list_fields = ['DTABLE']
        print('***default_values = %s' % self.default_values)
        assert len(self.default_values) > 0, self.default_values
        for label, value in sorted(iteritems(self.default_values)):
            list_fields += [label, value]
        return list_fields

    #def repr_fields(self):
        #return self.raw_fields()

    def write_card(self, size=8, is_double=False):
        card = self.repr_fields()
        if size == 8:
            return self.comment + print_card_8(card)
        if is_double:
            return self.comment + print_card_double(card)
        return self.comment + print_card_16(card)


class TABLED1(Table):
    type = 'TABLED1'
    #def __init__(self, tid, xaxis, yaxis, xy, comment=''):
    def __init__(self, tid, x, y, xaxis='LINEAR', yaxis='LINEAR', comment=''):
        Table.__init__(self)
        if comment:
            self._comment = comment
        self.tid = tid
        self.x = x
        self.y = y
        self.xaxis = xaxis
        self.yaxis = yaxis
        assert self.xaxis in ['LINEAR', 'LOG'], 'xaxis=%r' % (self.xaxis)
        assert self.yaxis in ['LINEAR', 'LOG'], 'yaxis=%r' % (self.yaxis)

    @classmethod
    def add_card(cls, card, comment=''):
        tid = integer(card, 1, 'tid')
        xaxis = string_or_blank(card, 2, 'xaxis', 'LINEAR')
        yaxis = string_or_blank(card, 3, 'yaxis', 'LINEAR')

        nfields = len(card) - 1
        nterms = (nfields - 9) // 2
        if nterms < 0:
            raise SyntaxError('%r card is too short' % cls.type)
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
        xy = np.array(xy)
        x = xy[:, 0]
        y = xy[:, 1]
        return TABLED1(tid, x, y, xaxis=xaxis, yaxis=yaxis, comment=comment)

    @classmethod
    def add_op2_data(cls, data, comment=''):
        tid = data[0]
        xaxis = cls._map_axis(data[1])
        yaxis = cls._map_axis(data[2])
        xy = data[3:]
        xy = np.array(xy)
        xy.reshape(xy.size // 2, 2)
        x = xy[:, 0]
        y = xy[:, 1]
        return TABLED1(tid, x, y, xaxis=xaxis, yaxis=yaxis, comment=comment)

    def raw_fields(self):
        xy = []
        for xi, yi in zip(self.x, self.y):
            xy.extend([xi, yi])
        list_fields = ['TABLED1', self.tid, self.xaxis, self.yaxis, None,
                       None, None, None, None] + xy + ['ENDT']
        return list_fields

    def repr_fields(self):
        #xaxis = set_blank_if_default(self.xaxis, 'LINEAR')
        #yaxis = set_blank_if_default(self.yaxis, 'LINEAR')
        return self.raw_fields()

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
    type = 'TABLED2'
    def __init__(self, tid, x1, x, y, comment=''):
        Table.__init__(self)
        if comment:
            self._comment = comment
        self.tid = tid
        self.x1 = x1
        self.x = x
        self.y = y

    @classmethod
    def add_card(cls, card, comment=''):
        tid = integer(card, 1, 'tid')
        x1 = double(card, 2, 'x1')

        nfields = len(card) - 1
        nterms = (nfields - 9) // 2
        if nterms < 0:
            raise SyntaxError('%r card is too short' % cls.type)
        xy = []
        for i in range(nterms):
            n = 9 + i * 2
            if card.field(n) == 'ENDT':
                break
            x = double_or_string(card, n, 'x' + str(i + 1))
            y = double_or_string(card, n + 1, 'y' + str(i + 1))
            if x == 'SKIP' or y == 'SKIP':
                continue
            xy.append([x, y])
        string(card, nfields, 'ENDT')
        xy = np.array(xy)
        x = xy[:, 0]
        y = xy[:, 1]
        return TABLED2(tid, x1, x, y, comment=comment)

    @classmethod
    def add_op2_data(cls, data, comment=''):
        tid = data[0]
        x1 = data[1]
        xy = data[2:]
        xy = np.array(xy)
        xy.reshape(xy.size // 2, 2)
        x = xy[:, 0]
        y = xy[:, 1]
        return TABLED2(tid, x1, x, y, comment=comment)

    def raw_fields(self):
        xy = []
        for xi, yi in zip(self.x, self.y):
            xy.extend([xi, yi])
        list_fields = ['TABLED2', self.tid, self.x1, None, None, None,
                       None, None, None] + xy + ['ENDT']
        return list_fields

    def repr_fields(self):
        return self.raw_fields()


class TABLED3(Table):
    type = 'TABLED3'
    def __init__(self, tid, x1, x2, x, y, comment=''):
        Table.__init__(self)
        if comment:
            self._comment = comment
        self.tid = tid
        self.x1 = x1
        self.x2 = x2
        self.x = x
        self.y = y
        assert self.x2 != 0.0

    @classmethod
    def add_card(cls, card, comment=''):
        tid = integer(card, 1, 'tid')
        x1 = double(card, 2, 'x1')
        x2 = double(card, 3, 'x2')

        nfields = len(card) - 1
        nterms = (nfields - 9) // 2
        if nterms < 0:
            raise SyntaxError('%r card is too short' % cls.type)
        xy = []
        for i in range(nterms):
            n = 9 + i * 2
            if card.field(n) == 'ENDT':
                break
            x = double_or_string(card, n, 'x' + str(i + 1))
            y = double_or_string(card, n + 1, 'y' + str(i + 1))
            if x == 'SKIP' or y == 'SKIP':
                continue
            xy.append([x, y])
        string(card, nfields, 'ENDT')
        xy = np.array(xy)
        x = xy[:, 0]
        y = xy[:, 1]
        return TABLED3(tid, x1, x2, x, y, comment=comment)

    @classmethod
    def add_op2_data(cls, data, comment=''):
        tid = data[0]
        x1 = data[1]
        x2 = data[2]
        xy = data[3:]
        xy = np.array(xy)
        x = xy[:, 0]
        y = xy[:, 1]
        return TABLED3(tid, x1, x2, x, y, comment=comment)

    def raw_fields(self):
        xy = []
        for xi, yi in zip(self.x, self.y):
            xy.extend([xi, yi])
        list_fields = ['TABLED3', self.tid, self.x1, self.x2, None,
                       None, None, None, None] + xy + ['ENDT']
        return list_fields

    def repr_fields(self):
        return self.raw_fields()


class TABLED4(Table):
    type = 'TABLED4'
    def __init__(self, tid, x1, x2, x3, x4, a, comment=''):
        Table.__init__(self)
        if comment:
            self._comment = comment
        self.tid = tid
        self.x1 = x1
        self.x2 = x2
        self.x3 = x3
        self.x4 = x4
        self.a = np.array(a)

        assert self.x2 != 0.0
        assert self.x3 < self.x4

    @classmethod
    def add_card(cls, card, comment=''):
        tid = integer(card, 1, 'tid')
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
        return TABLED4(tid, x1, x2, x3, x4, a, comment=comment)

    @classmethod
    def add_op2_data(cls, data, comment=''):
        tid = data[0]
        x1 = data[1]
        x2 = data[2]
        x3 = data[3]
        x4 = data[4]
        a = data[5:]
        return TABLED4(tid, x1, x2, x3, x4, a, comment=comment)

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

        yi = self.a * ((x - x1) / x2) ** n
        return yi.sum()


class TABDMP1(Table):
    type = 'TABDMP1'
    def __init__(self, tid, Type, x, y, comment=''):
        Table.__init__(self)
        if comment:
            self._comment = comment
        self.tid = tid
        self.Type = Type
        self.x = x
        self.y = y
        assert self.Type in ['G', 'CRIT', 'Q'], 'Type=%r' % self.Type

    @classmethod
    def add_card(cls, card, comment=''):
        tid = integer(card, 1, 'tid')
        Type = string_or_blank(card, 2, 'Type', 'G')

        nfields = len(card) - 1
        nterms = (nfields - 9) // 2
        if nterms < 0:
            raise SyntaxError('%r card is too short' % cls.type)
        xy = []
        for i in range(nterms):
            n = 9 + i * 2
            if card.field(n) == 'ENDT':
                break
            x = double_or_string(card, n, 'x' + str(i + 1))
            y = double_or_string(card, n + 1, 'y' + str(i + 1))
            if x == 'SKIP' or y == 'SKIP':
                continue
            xy.append([x, y])
        string(card, nfields, 'ENDT')

        xy = np.array(xy)
        x = xy[:, 0]
        y = xy[:, 1]
        return TABDMP1(tid, Type, x, y, comment=comment)

    @classmethod
    def add_op2_data(cls, data, comment=''):
        tid = data[0]
        x1 = data[1]
        Type = data[2]
        xy = data[5:]
        xy = np.array(xy)
        x = xy[:, 0]
        y = xy[:, 1]
        return TABDMP1(tid, Type, x, y, comment=comment)

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
    type = 'TABLEM1'
    def __init__(self, tid, x, y, comment=''):
        Table.__init__(self)
        if comment:
            self._comment = comment
        self.tid = tid
        self.x = x
        self.y = y

    @classmethod
    def add_card(cls, card, comment=''):
        tid = integer(card, 1, 'tid')

        nfields = len(card) - 1
        nterms = (nfields - 9) // 2
        if nterms < 0:
            raise SyntaxError('%r card is too short' % cls.type)
        xy = []
        for i in range(nterms):
            n = 9 + i * 2
            if card.field(n) == 'ENDT':
                break
            x = double_or_string(card, n, 'x' + str(i + 1))
            y = double_or_string(card, n + 1, 'y' + str(i + 1))
            if x == 'SKIP' or y == 'SKIP':
                continue
            xy.append([x, y])
        string(card, nfields, 'ENDT')
        xy = np.array(xy)
        x = xy[:, 0]
        y = xy[:, 1]
        return TABLEM1(tid, x, y, comment=comment)

    @classmethod
    def add_op2_data(cls, data, comment=''):
        tid = data[0]
        xy = data[1:]
        xy = np.array(xy)
        xy.reshape(xy.size // 2, 2)
        x = xy[:, 0]
        y = xy[:, 1]
        return TABLEM1(tid, x, y, comment=comment)


    def raw_fields(self):
        xy = []
        for xi, yi in zip(self.x, self.y):
            xy.extend([xi, yi])
        list_fields = ['TABLEM1', self.tid, None, None, None, None,
                       None, None, None] + xy + ['ENDT']
        return list_fields


class TABLEM2(Table):
    type = 'TABLEM2'
    def __init__(self, tid, x1, x, y, comment=''):
        Table.__init__(self)
        if comment:
            self._comment = comment
        self.tid = tid
        self.x1 = x1
        self.x = x
        self.y = y

    @classmethod
    def add_card(cls, card, comment=''):
        tid = integer(card, 1, 'tid')
        x1 = double(card, 2, 'x1')

        nfields = len(card) - 1
        nterms = (nfields - 9) // 2
        if nterms < 0:
            raise SyntaxError('%r card is too short' % cls.type)
        xy = []
        for i in range(nterms):
            n = 9 + i * 2
            if card.field(n) == 'ENDT':
                break
            x = double_or_string(card, n, 'x' + str(i + 1))
            y = double_or_string(card, n + 1, 'y' + str(i + 1))
            if x == 'SKIP' or y == 'SKIP':
                continue
            xy.append([x, y])
        string(card, nfields, 'ENDT')
        xy = np.array(xy)
        x = xy[:, 0]
        y = xy[:, 1]
        return TABLEM2(tid, x1, x, y, comment=comment)

    @classmethod
    def add_op2_data(cls, data, comment=''):
        tid = data[0]
        x1 = data[1]
        xy = data[2:]
        xy = np.array(xy)
        xy.reshape(xy.size // 2, 2)
        x = xy[:, 0]
        y = xy[:, 1]
        return TABLEM2(tid, x1, x, y, comment=comment)

    def raw_fields(self):
        xy = []
        for xi, yi in zip(self.x, self.y):
            xy.extend([xi, yi])
        list_fields = ['TABLEM2', self.tid, self.x1, None, None, None,
                       None, None, None] + xy + ['ENDT']
        return list_fields

    def repr_fields(self):
        return self.raw_fields()


class TABLEM3(Table):
    type = 'TABLEM3'
    def __init__(self, tid, x1, x2, x, y, comment=''):
        Table.__init__(self)
        if comment:
            self._comment = comment
        self.tid = tid
        self.x1 = x1
        self.x2 = x2
        self.x = x
        self.y = y
        assert self.x2 != 0.0

    @classmethod
    def add_card(cls, card, comment=''):
        tid = integer(card, 1, 'tid')
        x1 = double(card, 2, 'x1')
        x2 = double(card, 3, 'x2')

        nfields = len(card) - 1
        nterms = (nfields - 9) // 2
        if nterms < 0:
            raise SyntaxError('%r card is too short' % cls.type)
        xy = []
        for i in range(nterms):
            n = 9 + i * 2
            if card.field(n) == 'ENDT':
                break
            x = double_or_string(card, n, 'x' + str(i + 1))
            y = double_or_string(card, n + 1, 'y' + str(i + 1))
            if x == 'SKIP' or y == 'SKIP':
                continue
            xy.append([x, y])
        string(card, nfields, 'ENDT')
        xy = np.array(xy)
        x = xy[:, 0]
        y = xy[:, 1]
        return TABLEM3(tid, x1, x2, x, y, comment=comment)

    @classmethod
    def add_op2_data(cls, data, comment=''):
        tid = data[0]
        x1 = data[1]
        x2 = data[2]
        xy = data[3:]
        xy = np.array(xy)
        xy.reshape(xy.size // 2, 2)
        x = xy[:, 0]
        y = xy[:, 1]
        return TABLEM3(tid, x1, x2, x, y, comment=comment)

    def raw_fields(self):
        xy = []
        for xi, yi in zip(self.x, self.y):
            xy.extend([xi, yi])
        list_fields = ['TABLEM3', self.tid, self.x1, self.x2, None,
                       None, None, None, None] + xy + ['ENDT']
        return list_fields

    def repr_fields(self):
        return self.raw_fields()


class TABLEM4(Table):
    type = 'TABLEM4'
    def __init__(self, tid, x1, x2, x3, x4, a, comment=''):
        Table.__init__(self)
        if comment:
            self._comment = comment
        self.tid = tid
        self.x1 = x1
        self.x2 = x2
        self.x3 = x3
        self.x4 = x4
        self.a = a
        assert self.x2 != 0.0
        assert self.x3 < self.x4

    @classmethod
    def add_card(cls, card, comment=''):
        tid = integer(card, 1, 'tid')
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
        return TABLEM4(tid, x1, x2, x3, x4, a, comment=comment)

    @classmethod
    def add_op2_data(cls, data, comment=''):
        tid = data[0]
        x1 = data[1]
        x2 = data[2]
        x3 = data[3]
        x4 = data[4]
        a = data[3:]
        return TABLEM4(tid, x1, x2, x3, x4, a, comment=comment)

    def raw_fields(self):
        list_fields = ['TABLEM4', self.tid, self.x1, self.x2, self.x3, self.x4,
                       None, None, None] + list(self.a) + ['ENDT']
        return list_fields

    def repr_fields(self):
        return self.raw_fields()


class TABLES1(Table):
    type = 'TABLES1'

    def __init__(self, tid, Type, x, y, comment=''):
        Table.__init__(self)
        if comment:
            self._comment = comment
        self.tid = tid
        self.Type = Type
        self.x = x
        self.y = y
        assert self.Type in [1, 2], 'TABLES1 Type=%s' % self.Type

    @classmethod
    def add_card(cls, card, comment=''):
        tid = integer(card, 1, 'tid')
        Type = integer_or_blank(card, 2, 'Type', 1)

        nfields = len(card) - 1
        nterms = (nfields - 9) // 2
        if nterms < 0:
            raise SyntaxError('%r card is too short' % cls.type)
        xy = []
        for i in range(nterms):
            n = 9 + i * 2
            if card.field(n) == 'ENDT':
                break
            x = double_or_string(card, n, 'x' + str(i + 1))
            y = double_or_string(card, n + 1, 'y' + str(i + 1))
            if x == 'SKIP' or y == 'SKIP':
                continue
            xy.append([x, y])
        string(card, nfields, 'ENDT')
        xy = np.array(xy)
        x = xy[:, 0]
        y = xy[:, 1]
        return TABLES1(tid, Type, x, y, comment=comment)

    @classmethod
    def add_op2_data(cls, data, comment=''):
        tid = data[0]
        xy = data[1:]
        xy = np.array(xy)
        xy.reshape(xy.size // 2, 2)
        x = xy[:, 0]
        y = xy[:, 1]
        return TABLES1(tid, Type, x, y, comment=comment)

    def raw_fields(self):
        xy = []
        for xi, yi in zip(self.x, self.y):
            xy.extend([xi, yi])
        list_fields = ['TABLES1', self.tid, self.Type, None, None, None,
                       None, None, None] + xy + ['ENDT']
        return list_fields

    def repr_fields(self):
        return self.raw_fields()


class TABLEST(Table):
    type = 'TABLEST'

    def __init__(self, tid, x, y, comment=''):
        Table.__init__(self)
        if comment:
            self._comment = comment
        self.tid = tid
        self.x = x
        self.y = y

    @classmethod
    def add_card(cls, card, comment=''):
        tid = integer(card, 1, 'tid')

        nfields = len(card) - 1
        nterms = (nfields - 9) // 2
        if nterms < 0:
            raise SyntaxError('%r card is too short' % cls.type)
        xy = []
        for i in range(nterms):
            n = 9 + i * 2
            if card.field(n) == 'ENDT':
                break
            x = double_or_string(card, n, 'x' + str(i + 1))
            y = double_or_string(card, n + 1, 'y' + str(i + 1))
            if x == 'SKIP' or y == 'SKIP':
                continue
            xy.append([x, y])
        string(card, nfields, 'ENDT')
        xy = np.array(xy)
        x = xy[:, 0]
        y = xy[:, 1]
        return TABLEST(tid, x, y, comment=comment)

    @classmethod
    def add_op2_data(cls, data, comment=''):
        tid = data[0]
        xy = data[1:]
        xy = np.array(xy)
        xy.reshape(xy.size // 2, 2)
        x = xy[:, 0]
        y = xy[:, 1]
        return TABLEST(tid, x, y, comment=comment)

    def raw_fields(self):
        xy = []
        for xi, yi in zip(self.x, self.y):
            xy.extend([xi, yi])
        list_fields = ['TABLEST', self.tid, None, None, None, None,
                       None, None, None] + xy + ['ENDT']
        return list_fields

    def repr_fields(self):
        return self.raw_fields()


class RandomTable(Table):
    type = 'TABLE??'

    def __init__(self):
        pass


class TABRND1(RandomTable):
    type = 'TABRND1'

    def __init__(self, tid, x, y, xaxis='LINEAR', yaxis='LINEAR', comment=''):
        RandomTable.__init__(self)
        if comment:
            self._comment = comment

        self.tid = tid
        self.x = x
        self.y = y
        self.xaxis = xaxis
        self.yaxis = yaxis
        assert self.xaxis in ['LINEAR', 'LOG'], 'xaxis=%r' % (self.xaxis)
        assert self.yaxis in ['LINEAR', 'LOG'], 'yaxis=%r' % (self.yaxis)

    @classmethod
    def add_card(cls, card, comment=''):
        tid = integer(card, 1, 'tid')
        xaxis = string_or_blank(card, 2, 'xaxis', 'LINEAR')
        yaxis = string_or_blank(card, 3, 'yaxis', 'LINEAR')

        nfields = len(card) - 1
        nterms = (nfields - 9) // 2
        if nterms < 0:
            raise SyntaxError('%r card is too short' % cls.type)
        xy = []
        for i in range(nterms):
            n = 9 + i * 2
            if card.field(n) == 'ENDT':
                break
            x = double_or_string(card, n, 'x' + str(i + 1))
            y = double_or_string(card, n + 1, 'y' + str(i + 1))
            if x == 'SKIP' or y == 'SKIP':
                continue
            xy.append([x, y])
        string(card, nfields, 'ENDT')
        xy = np.array(xy)
        x = xy[:, 0]
        y = xy[:, 1]
        return TABRND1(tid, x, y, xaxis=xaxis, yaxis=yaxis, comment=comment)

    @classmethod
    def add_op2_data(cls, data, comment=''):
        tid = data[0]
        xaxis = cls._map_axis(data[1])
        yaxis = cls._map_axis(data[2])
        xy = data[3:]
        xy = np.array(xy)
        x = xy[:, 0]
        y = xy[:, 1]
        return TABRND1(tid, x, y, xaxis=xaxis, yaxis=yaxis, comment=comment)

    #def parse_fields(self, xy, nrepeated, is_data=False):
        #self.table = TableObj(xy, nrepeated, is_data)

    def _map_axis(self, axis):
        if axis == 0:
            axis_type = 'LINEAR'
        else:
            raise ValueError('axis=%r' % (axis))
        return axis_type

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


class TABRNDG(RandomTable):
    r"""
    Gust Power Spectral Density

    Defines the power spectral density (PSD) of a gust for aeroelastic response
    analysis.
    """
    type = 'TABRNDG'

    def __init__(self, tid, Type, LU, WG, comment=''):
        RandomTable.__init__(self)
        if comment:
            self._comment = comment
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
        tid = integer(card, 1, 'tid')
        Type = integer(card, 2, 'Type')
        LU = double(card, 3, 'LU')
        WG = double(card, 4, 'WG')
        return TABRNDG(tid, Type, LU, WG, comment=comment)

    def raw_fields(self):
        list_fields = ['TABRNDG', self.tid, self.Type, self.LU, self.WG]
        return list_fields

    def repr_fields(self):
        return self.raw_fields()


class TIC(Table):
    """Transient Initial Condition"""
    type = 'TIC'

    def __init__(self, card=None, data=None, comment=''):
        """
        Defines values for the initial conditions of variables used in
        structural transient analysis. Both displacement and velocity values
        may be specified at independent degrees-of-freedom. This entry may not
        be used for heat transfer analysis.
        """
        Table.__init__(self, card, data)
        if comment:
            self._comment = comment
        if card:
            self.sid = integer(card, 1, 'sid')
            self.G = integer(card, 2, 'G')
            assert self.G > 0
            self.C = components(card, 3, 'C')
            self.U0 = double(card, 4, 'U0')
            self.V0 = double(card, 5, 'V0')
        else:
            self.sid = data[0]
            self.G = data[1]
            self.C = data[2]
            self.U0 = data[3]
            self.V0 = data[4]

    def raw_fields(self):
        list_fields = ['TIC', self.sid, self.G, self.C, self.U0, self.V0]
        return list_fields

    def repr_fields(self):
        return self.raw_fields()
