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
from six import string_types, iteritems
from six.moves import range

from pyNastran.bdf.field_writer_8 import set_blank_if_default, print_card_8
from pyNastran.bdf.field_writer_16 import print_card_16
from pyNastran.bdf.field_writer_double import print_card_double

from pyNastran.bdf.cards.base_card import BaseCard
from pyNastran.utils.dev import list_print
from pyNastran.bdf.bdf_interface.assign_type import (integer,
    integer_or_blank, double, components, string, string_or_blank,
    double_or_string)

class Table(BaseCard):
    def __init__(self, card, data):
        pass

    def map_axis(self, axis):
        if axis == 0:
            axisType = 'LINEAR'
        else:
            raise ValueError('axis=|%s|' % axis)
        return axisType

    def parse_fields(self, xy, nrepeated, is_data=False):
        self.table = TableObj(xy, nrepeated, is_data)

    def write_card(self, size=8, is_double=False):
        card = self.repr_fields()
        if size == 8:
            return self.comment + print_card_8(card)
        if is_double:
            return self.comment + print_card_double(card)
        return self.comment + print_card_16(card)


class TableObj(object):
    def __init__(self, xy, nrepeated, is_data=False):
        """
        Parameters
        ----------
        xy : List[float/str]
            the X/Y data with an ENDT appended
        nrepeated : int
            ???
        is_data : bool
            did this come from the OP2/BDF (True -> OP2)
        """
        self.table = []
        xy = self._cleanup_xy(xy, is_data)

        nxy = len(xy)

        if not is_data:
            if nxy % nrepeated != 0:
                self._crash_fields(xy, nrepeated, nxy)
            if nxy % nrepeated != 0:
                msg = 'invalid table length nrepeat=%s xy=%s' % (
                    nrepeated, list_print(xy))
                raise RuntimeError(msg)

        i = 0
        while i < nxy:
            pack = []
            for j in range(nrepeated):
                pack.append(xy[i + j])
            i += nrepeated
            self.table.append(pack)

    def _crash_fields(self, xy, nrepeated, nxy):
        """
        Creates the print message if there was an error

        :param xy:        the xy data as a table with alternating x, y entries
        :param nrepeated: ???
        :param nxy:       ???
        """
        try:
            msg = ''
            for i in range(nxy):
                for j in range(nrepeated):
                    try:
                        msg += '%-8g ' % xy[i * nrepeated + j]
                    except TypeError:
                        msg += '*%-8s ' % xy[i * nrepeated + j]
                    except IndexError:
                        msg += 'IndexError'
                msg += '\n'
        except:
            print(xy)
            print(msg)
            #assert nxy%nrepeated == 0, msg
            raise

    def _cleanup_xy(self, xy, is_data=False):
        """
        Removes the **ENDT** field.

        :param xy:     the xy data as a table with alternating x, y entries
        :param is_data: did this come from the OP2/BDF (True -> OP2)
        """
        xy2 = []  # remove extra ENDTs

        if 1:  # hardcoded b/c ENDT has been removed
            return xy

        foundENDT = False
        for value in xy:
            if isinstance(value, string_types) and 'ENDT' in value.upper():
                foundENDT = True
            else:
                xy2.append(value)

        if not is_data:
            assert foundENDT == True, xy
        return xy2

    def fields(self):
        list_fields = []
        for pack in self.table:
            list_fields += pack
        return list_fields


class DTABLE(BaseCard):
    type = 'DTABLE'
    def __init__(self, card=None, data=None, comment=''):
        if comment:
            self._comment = comment
        nfields = len(card) - 1
        assert nfields % 2 == 0, nfields

        self.default_values = {}
        j = 1
        for i in range(1, nfields + 1, 2):
            label = string(card, i, 'label_%i' % j)
            value = double(card, i + 1, 'value_%i' % j)
            assert label not in self.default_values, 'label_%i=%r is not unique' % (j, label)
            self.default_values[label] = value
            j += 1

    def __getitem__(self, key):
        return self.default_values[key]

    def raw_fields(self):
        list_fields = ['DTABLE']
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
    def __init__(self, card=None, data=None, comment=''):
        Table.__init__(self, card, data)
        if comment:
            self._comment = comment
        if card:
            self.tid = integer(card, 1, 'tid')
            self.xaxis = string_or_blank(card, 2, 'xaxis', 'LINEAR')
            self.yaxis = string_or_blank(card, 3, 'yaxis', 'LINEAR')

            nfields = len(card) - 1
            nterms = (nfields - 9) // 2
            if nterms < 0:
                raise SyntaxError('%r card is too short' % self.type)
            xy = []
            for i in range(nterms):
                n = 9 + i * 2
                if card.field(n) == 'ENDT':
                    break
                x = double_or_string(card, n, 'x' + str(i + 1))
                y = double_or_string(card, n + 1, 'y' + str(i + 1))
                if x == 'SKIP' or y == 'SKIP':
                    continue
                xy += [x, y]
            string(card, nfields, 'ENDT')
            is_data = False
        else:
            self.tid = data[0]
            self.xaxis = self.map_axis(data[1])
            self.yaxis = self.map_axis(data[2])
            xy = data[3:]
            is_data = True
        assert self.xaxis in ['LINEAR', 'LOG'], 'xaxis=%r' % (self.xaxis)
        assert self.yaxis in ['LINEAR', 'LOG'], 'yaxis=%r' % (self.yaxis)
        self.parse_fields(xy, nrepeated=2, is_data=is_data)

    def raw_fields(self):
        list_fields = ['TABLED1', self.tid, self.xaxis, self.yaxis, None,
                       None, None, None, None] + self.table.fields() + ['ENDT']
        return list_fields

    def repr_fields(self):
        #xaxis = set_blank_if_default(self.xaxis, 'LINEAR')
        #yaxis = set_blank_if_default(self.yaxis, 'LINEAR')
        return self.raw_fields()


class TABLED2(Table):
    type = 'TABLED2'
    def __init__(self, card=None, data=None, comment=''):
        Table.__init__(self, card, data)
        if comment:
            self._comment = comment
        if card:
            self.tid = integer(card, 1, 'tid')
            self.x1 = double(card, 2, 'x1')

            nfields = len(card) - 1
            nterms = (nfields - 9) // 2
            if nterms < 0:
                raise SyntaxError('%r card is too short' % self.type)
            xy = []
            for i in range(nterms):
                n = 9 + i * 2
                if card.field(n) == 'ENDT':
                    break
                x = double_or_string(card, n, 'x' + str(i + 1))
                y = double_or_string(card, n + 1, 'y' + str(i + 1))
                if x == 'SKIP' or y == 'SKIP':
                    continue
                xy += [x, y]
            string(card, nfields, 'ENDT')
            is_data = False
        else:
            self.tid = data[0]
            self.x1 = data[1]
            xy = data[2:]
            is_data = True
        self.parse_fields(xy, nrepeated=2, is_data=is_data)

    def raw_fields(self):
        list_fields = ['TABLED2', self.tid, self.x1, None, None, None,
                       None, None, None] + self.table.fields() + ['ENDT']
        return list_fields

    def repr_fields(self):
        return self.raw_fields()


class TABLED3(Table):
    type = 'TABLED3'
    def __init__(self, card=None, data=None, comment=''):
        Table.__init__(self, card, data)
        if comment:
            self._comment = comment
        if card:
            self.tid = integer(card, 1, 'tid')
            self.x1 = double(card, 2, 'x1')
            self.x2 = double(card, 3, 'x2')
            assert self.x2 != 0.0

            nfields = len(card) - 1
            nterms = (nfields - 9) // 2
            if nterms < 0:
                raise SyntaxError('%r card is too short' % self.type)
            xy = []
            for i in range(nterms):
                n = 9 + i * 2
                if card.field(n) == 'ENDT':
                    break
                x = double_or_string(card, n, 'x' + str(i + 1))
                y = double_or_string(card, n + 1, 'y' + str(i + 1))
                if x == 'SKIP' or y == 'SKIP':
                    continue
                xy += [x, y]
            string(card, nfields, 'ENDT')
            is_data = False
        else:
            self.tid = data[0]
            self.x1 = data[1]
            self.x2 = data[2]
            xy = data[3:]
            is_data = True
        self.parse_fields(xy, nrepeated=2, is_data=is_data)

    def raw_fields(self):
        list_fields = ['TABLED3', self.tid, self.x1, self.x2, None,
                       None, None, None, None] + self.table.fields() + ['ENDT']
        return list_fields

    def repr_fields(self):
        return self.raw_fields()


class TABLED4(Table):
    type = 'TABLED4'
    def __init__(self, card=None, data=None, comment=''):
        Table.__init__(self, card, data)
        if comment:
            self._comment = comment
        if card:
            self.tid = integer(card, 1, 'tid')
            self.x1 = double(card, 2, 'x1')
            self.x2 = double(card, 3, 'x2')
            self.x3 = double(card, 4, 'x3')
            self.x4 = double(card, 5, 'x4')

            assert self.x2 != 0.0
            assert self.x3 < self.x4

            nfields = len(card) - 1
            nterms = (nfields - 9) // 2
            if nterms < 0:
                raise SyntaxError('%r card is too short' % self.type)
            xy = []
            for i in range(nterms):
                n = 9 + i * 2
                if card.field(n) == 'ENDT':
                    break
                x = double_or_string(card, n, 'x' + str(i + 1))
                y = double_or_string(card, n + 1, 'y' + str(i + 1))
                if x == 'SKIP' or y == 'SKIP':
                    continue
                xy += [x, y]
            string(card, nfields, 'ENDT')
            is_data = False
        else:
            self.tid = data[0]
            self.x1 = data[1]
            self.x2 = data[2]
            self.x3 = data[3]
            self.x4 = data[4]
            xy = data[5:]
            is_data = True
        self.parse_fields(xy, nrepeated=2, is_data=is_data)

    def raw_fields(self):
        list_fields = ['TABLED4', self.tid, self.x1, self.x2, self.x3, self.x4,
                       None, None, None] + self.table.fields() + ['ENDT']
        return list_fields

    def repr_fields(self):
        return self.raw_fields()


class TABDMP1(Table):
    type = 'TABDMP1'
    def __init__(self, card=None, data=None, comment=''):
        Table.__init__(self, card, data)
        if comment:
            self._comment = comment
        if card:
            self.tid = integer(card, 1, 'tid')
            self.Type = string_or_blank(card, 2, 'Type', 'G')
            assert self.Type in ['G', 'CRIT', 'Q'], 'Type=%r' % self.Type

            nfields = len(card) - 1
            nterms = (nfields - 9) // 2
            if nterms < 0:
                raise SyntaxError('%r card is too short' % self.type)
            xy = []
            for i in range(nterms):
                n = 9 + i * 2
                if card.field(n) == 'ENDT':
                    break
                x = double_or_string(card, n, 'x' + str(i + 1))
                y = double_or_string(card, n + 1, 'y' + str(i + 1))
                if x == 'SKIP' or y == 'SKIP':
                    continue
                xy += [x, y]
            string(card, nfields, 'ENDT')
            is_data = False
        else:
            self.tid = data[0]
            self.x1 = data[1]
            self.Type = data[2]
            xy = data[5:]
            is_data = True
        self.parse_fields(xy, nrepeated=2, is_data=is_data)

    def raw_fields(self):
        list_fields = ['TABDMP1', self.tid, self.Type, None, None, None, None,
                       None, None] + self.table.fields() + ['ENDT']
        return list_fields

    def repr_fields(self):
        return self.raw_fields()


class TABLEM1(Table):
    type = 'TABLEM1'
    def __init__(self, card=None, data=None, comment=''):
        Table.__init__(self, card, data)
        if comment:
            self._comment = comment
        if card:
            self.tid = integer(card, 1, 'tid')

            nfields = len(card) - 1
            nterms = (nfields - 9) // 2
            if nterms < 0:
                raise SyntaxError('%r card is too short' % self.type)
            xy = []
            for i in range(nterms):
                n = 9 + i * 2
                if card.field(n) == 'ENDT':
                    break
                x = double_or_string(card, n, 'x' + str(i + 1))
                y = double_or_string(card, n + 1, 'y' + str(i + 1))
                if x == 'SKIP' or y == 'SKIP':
                    continue
                xy += [x, y]
            string(card, nfields, 'ENDT')
            is_data = False
        else:
            self.tid = data[0]
            xy = data[1:]
            is_data = True
        self.parse_fields(xy, nrepeated=2, is_data=is_data)

    def raw_fields(self):
        list_fields = ['TABLEM1', self.tid, None, None, None, None,
                       None, None, None] + self.table.fields() + ['ENDT']
        return list_fields


class TABLEM2(Table):
    type = 'TABLEM2'
    def __init__(self, card=None, data=None, comment=''):
        Table.__init__(self, card, data)
        if comment:
            self._comment = comment
        if card:
            self.tid = integer(card, 1, 'tid')
            self.x1 = double(card, 2, 'x1')

            nfields = len(card) - 1
            nterms = (nfields - 9) // 2
            if nterms < 0:
                raise SyntaxError('%r card is too short' % self.type)
            xy = []
            for i in range(nterms):
                n = 9 + i * 2
                if card.field(n) == 'ENDT':
                    break
                x = double_or_string(card, n, 'x' + str(i + 1))
                y = double_or_string(card, n + 1, 'y' + str(i + 1))
                if x == 'SKIP' or y == 'SKIP':
                    continue
                xy += [x, y]
            string(card, nfields, 'ENDT')
            is_data = False
        else:
            self.tid = data[0]
            self.x1 = data[1]
            xy = data[2:]
            is_data = True
        self.parse_fields(xy, nrepeated=2, is_data=is_data)

    def raw_fields(self):
        list_fields = ['TABLEM2', self.tid, self.x1, None, None, None,
                       None, None, None] + self.table.fields() + ['ENDT']
        return list_fields

    def repr_fields(self):
        return self.raw_fields()


class TABLEM3(Table):
    type = 'TABLEM3'
    def __init__(self, card=None, data=None, comment=''):
        Table.__init__(self, card, data)
        if comment:
            self._comment = comment
        if card:
            self.tid = integer(card, 1, 'tid')
            self.x1 = double(card, 2, 'x1')
            self.x2 = double(card, 3, 'x2')
            assert self.x2 != 0.0

            nfields = len(card) - 1
            nterms = (nfields - 9) // 2
            if nterms < 0:
                raise SyntaxError('%r card is too short' % self.type)
            xy = []
            for i in range(nterms):
                n = 9 + i * 2
                if card.field(n) == 'ENDT':
                    break
                x = double_or_string(card, n, 'x' + str(i + 1))
                y = double_or_string(card, n + 1, 'y' + str(i + 1))
                if x == 'SKIP' or y == 'SKIP':
                    continue
                xy += [x, y]
            string(card, nfields, 'ENDT')
            is_data = False
        else:
            self.tid = data[0]
            self.x1 = data[1]
            self.x2 = data[2]
            xy = data[3:]
            is_data = True
        self.parse_fields(xy, nrepeated=2, is_data=is_data)

    def raw_fields(self):
        list_fields = ['TABLEM3', self.tid, self.x1, self.x2, None,
                       None, None, None, None] + self.table.fields() + ['ENDT']
        return list_fields

    def repr_fields(self):
        return self.raw_fields()


class TABLEM4(Table):
    type = 'TABLEM4'
    def __init__(self, card=None, data=None, comment=''):
        Table.__init__(self, card, data)
        if comment:
            self._comment = comment
        if card:
            self.tid = integer(card, 1, 'tid')
            self.x1 = double(card, 2, 'x1')
            self.x2 = double(card, 3, 'x2')
            assert self.x2 != 0.0
            self.x3 = double(card, 4, 'x3')
            self.x4 = double(card, 5, 'x4')
            assert self.x3 < self.x4

            nfields = len(card) - 1
            nterms = (nfields - 9) // 2
            if nterms < 0:
                raise SyntaxError('%r card is too short' % self.type)
            xy = []
            for i in range(nterms):
                n = 9 + i * 2
                if card.field(n) == 'ENDT':
                    break
                x = double_or_string(card, n, 'x' + str(i + 1))
                y = double_or_string(card, n + 1, 'y' + str(i + 1))
                if x == 'SKIP' or y == 'SKIP':
                    continue
                xy += [x, y]
            string(card, nfields, 'ENDT')
            is_data = False
        else:
            self.tid = data[0]
            self.x1 = data[1]
            self.x2 = data[2]
            self.x3 = data[3]
            self.x4 = data[4]
            xy = data[3:]
            is_data = True
        self.parse_fields(xy, nrepeated=1, is_data=is_data)

    def raw_fields(self):
        list_fields = ['TABLEM4', self.tid, self.x1, self.x2, self.x3, self.x4,
                       None, None, None] + self.table.fields() + ['ENDT']
        return list_fields

    def repr_fields(self):
        return self.raw_fields()


class TABLES1(Table):
    type = 'TABLES1'

    def __init__(self, card=None, data=None, comment=''):
        Table.__init__(self, card, data)
        if comment:
            self._comment = comment
        if card:
            self.tid = integer(card, 1, 'tid')
            self.Type = integer_or_blank(card, 2, 'Type', 1)
            assert self.Type in [1, 2], 'TABLES1 Type=%s' % self.Type

            nfields = len(card) - 1
            nterms = (nfields - 9) // 2
            if nterms < 0:
                raise SyntaxError('%r card is too short' % self.type)
            xy = []
            for i in range(nterms):
                n = 9 + i * 2
                if card.field(n) == 'ENDT':
                    break
                x = double_or_string(card, n, 'x' + str(i + 1))
                y = double_or_string(card, n + 1, 'y' + str(i + 1))
                if x == 'SKIP' or y == 'SKIP':
                    continue
                xy += [x, y]
            string(card, nfields, 'ENDT')
            is_data = False
        else:
            self.tid = data[0]
            xy = data[1:]
            is_data = True
        self.parse_fields(xy, nrepeated=2, is_data=is_data)

    def raw_fields(self):
        list_fields = ['TABLES1', self.tid, self.Type, None, None, None,
                       None, None, None] + self.table.fields() + ['ENDT']
        return list_fields

    def repr_fields(self):
        return self.raw_fields()


class TABLEST(Table):
    type = 'TABLEST'

    def __init__(self, card=None, data=None, comment=''):
        Table.__init__(self, card, data)
        if comment:
            self._comment = comment
        if card:
            self.tid = integer(card, 1, 'tid')

            nfields = len(card) - 1
            nterms = (nfields - 9) // 2
            if nterms < 0:
                raise SyntaxError('%r card is too short' % self.type)
            xy = []
            for i in range(nterms):
                n = 9 + i * 2
                if card.field(n) == 'ENDT':
                    break
                x = double_or_string(card, n, 'x' + str(i + 1))
                y = double_or_string(card, n + 1, 'y' + str(i + 1))
                if x == 'SKIP' or y == 'SKIP':
                    continue
                xy += [x, y]
            string(card, nfields, 'ENDT')
            is_data = False
        else:
            self.tid = data[0]
            xy = data[1:]
            is_data = True
        self.parse_fields(xy, nrepeated=2, is_data=is_data)

    def raw_fields(self):
        list_fields = ['TABLEST', self.tid, None, None, None, None,
                       None, None, None] + self.table.fields() + ['ENDT']
        return list_fields

    def repr_fields(self):
        return self.raw_fields()


class RandomTable(Table):
    type = 'TABLE??'

    def __init__(self, card, data):
        pass


class TABRND1(RandomTable):
    type = 'TABRND1'

    def __init__(self, card=None, data=None, comment=''):
        RandomTable.__init__(self, card, data)
        if comment:
            self._comment = comment
        if card:
            self.tid = integer(card, 1, 'tid')
            self.xaxis = string_or_blank(card, 2, 'xaxis', 'LINEAR')
            self.yaxis = string_or_blank(card, 3, 'yaxis', 'LINEAR')

            nfields = len(card) - 1
            nterms = (nfields - 9) // 2
            if nterms < 0:
                raise SyntaxError('%r card is too short' % self.type)
            xy = []
            for i in range(nterms):
                n = 9 + i * 2
                if card.field(n) == 'ENDT':
                    break
                x = double_or_string(card, n, 'x' + str(i + 1))
                y = double_or_string(card, n + 1, 'y' + str(i + 1))
                if x == 'SKIP' or y == 'SKIP':
                    continue
                xy += [x, y]
            string(card, nfields, 'ENDT')
            is_data = False
        else:
            self.tid = data[0]
            self.xaxis = self.map_axis(data[1])
            self.yaxis = self.map_axis(data[2])
            xy = data[3:]
            is_data = True
        assert self.xaxis in ['LINEAR', 'LOG'], 'xaxis=%r' % (self.xaxis)
        assert self.yaxis in ['LINEAR', 'LOG'], 'yaxis=%r' % (self.yaxis)
        self.parse_fields(xy, nrepeated=2, is_data=is_data)

    def parse_fields(self, xy, nrepeated, is_data=False):
        self.table = TableObj(xy, nrepeated, is_data)

    def map_axis(self, axis):
        if axis == 0:
            axisType = 'LINEAR'
        else:
            raise ValueError('axis=|%s|' % (axis))
        return axisType

    def raw_fields(self):
        list_fields = ['TABRND1', self.tid, self.xaxis, self.yaxis, None, None,
                       None, None, None] + self.table.fields() + ['ENDT']
        return list_fields

    def repr_fields(self):
        xaxis = set_blank_if_default(self.xaxis, 'LINEAR')
        yaxis = set_blank_if_default(self.yaxis, 'LINEAR')
        list_fields = ['TABRND1', self.tid, xaxis, yaxis, None, None,
                       None, None, None] + self.table.fields() + ['ENDT']
        return list_fields


class TABRNDG(RandomTable):
    r"""
    Gust Power Spectral Density

    Defines the power spectral density (PSD) of a gust for aeroelastic response
    analysis.
    """
    type = 'TABRNDG'

    def __init__(self, card=None, data=None, comment=''):
        RandomTable.__init__(self, card, data)
        if comment:
            self._comment = comment
        if card:
            #: Table identification number. (Integer >0)
            self.tid = integer(card, 1, 'tid')
            #: PSD Type: 1. von Karman; 2. Dryden
            self.Type = integer(card, 2, 'Type')
            #: Scale of turbulence divided by velocity (units of time; Real)
            self.LU = double(card, 3, 'LU')
            #: Root-mean-square gust velocity. (Real)
            self.WG = double(card, 4, 'WG')
            assert self.Type in [1, 2], ('Type must be 1 or 2.  '
                                         'Type=%s' % (self.Type))
        else:
            raise NotImplementedError()

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
