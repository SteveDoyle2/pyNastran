"""
All bush properties are defined in this file.  This includes:

 *   PBUSH
 *   PBUSH1D
 *   PBUSH2D (not implemented)
 *   PBUSHT (not implemented)

All bush properties are BushingProperty and Property objects.
"""
from __future__ import (nested_scopes, generators, division, absolute_import,
                        print_function, unicode_literals)

from pyNastran.utils import integer_types
from pyNastran.bdf.cards.base_card import Property
from pyNastran.bdf.bdf_interface.assign_type import (
    integer, integer_or_blank, double, double_or_blank, string,
    string_or_blank, blank, fields)
from pyNastran.bdf.field_writer_8 import print_card_8
from pyNastran.bdf.field_writer_16 import print_card_16

class BushingProperty(Property):
    type = 'BushingProperty'

    def __init__(self):
        Property.__init__(self)

    def cross_reference(self, model):
        pass

    def uncross_reference(self):
        pass


class PBUSH(BushingProperty):
    type = 'PBUSH'
    _field_map = {
        1: 'pid',
    }

    def __init__(self, pid, k, b, ge, rcv, comment=''):
        BushingProperty.__init__(self)
        if comment:
            self._comment = comment

        #: Property ID
        self.pid = pid
        self.vars = []

        # K parameter
        self.Ki = k
        if k:
            self.vars.append('K')
        # B parameter
        self.Bi = b
        if b:
            self.vars.append('B')

        # GE parameter
        self.GEi = ge
        if ge:
            self.vars.append('GE')

        # RCV parameters
        sa, st, ea, et = rcv
        if sa is not None or st is not None or ea is not None or et is not None:
            self.vars.append('RCV')
            self.sa = rcv[0]
            self.st = rcv[1]
            self.ea = rcv[2]
            self.et = rcv[3]

    @classmethod
    def add_card(cls, card, comment=''):
        k_fields = []
        b_fields = []
        ge_fields = []
        rcv_fields = [None, None, None, None]

        pid = integer(card, 1, 'pid')

        nfields = card.nfields
        istart = 2
        while istart < nfields:
            pname = string(card, istart, 'pname')
            if   pname == 'K':
                k_fields = cls._read_k(card, istart)
            elif pname == 'B':
                b_fields = cls._read_b(card, istart)
            elif pname == 'GE':
                ge_fields = cls._read_ge(card, istart)
            elif pname == 'RCV':
                rcv_fields = cls._read_rcv(card, istart)
            else:
                break
            istart += 8
        return PBUSH(pid, k_fields, b_fields, ge_fields, rcv_fields,
                     comment=comment)

    #@classmethod
    #def add_op2_data(cls, data, comment=''):
        #pid = data[0]
        #b = data[1]
        #return PBUSH(pid, k_fields, b_fields, ge_fields, rcv_fields,
                     #comment=comment)

    def _verify(self, xref=False):
        pid = self.Pid()
        assert isinstance(pid, integer_types), 'pid=%r' % pid

    @classmethod
    def _read_k(cls, card, istart):
        # Flag indicating that the next 1 to 6 fields are stiffness values in
        # the element coordinate system.
        #self.k = string(card, istart, 'k')

        #: Nominal stiffness values in directions 1 through 6.
        #: See Remarks 2 and 3. (Real; Default = 0.0)
        Ki = fields(double_or_blank, card, 'Ki', istart + 1, istart + 7)
        return Ki

    @classmethod
    def _read_b(cls, card, istart):
        # Flag indicating that the next 1 to 6 fields are force-per-velocity
        # damping.
        #self.b = string(card, istart, 'b')

        #: Force per unit velocity (Real)
        #: Nominal damping coefficients in direction 1 through 6 in units of
        #: force per unit velocity. See Remarks 2, 3, and 9. (Real; Default=0.)
        Bi = fields(double_or_blank, card, 'Bi', istart + 1, istart + 7)
        return Bi

    @classmethod
    def _read_ge(cls, card, istart):
        # Flag indicating that the next fields, 1 through 6 are structural
        # damping constants. See Remark 7. (Character)
        #self.ge = string(card, istart, 'ge')

        #: Nominal structural damping constant in directions 1 through 6. See
        #: Remarks 2. and 3. (Real; Default = 0.0)
        GEi = fields(double_or_blank, card, 'GEi', istart + 1, istart + 7)
        return GEi

    @classmethod
    def _read_rcv(cls, card, istart):
        # Flag indicating that the next 1 to 4 fields are stress or strain
        # coefficients. (Character)
        #self.rcv = string(card, istart, 'rcv')

        sa = double_or_blank(card, istart + 1, 'sa', 1.)
        st = double_or_blank(card, istart + 2, 'st', 1.)
        ea = double_or_blank(card, istart + 3, 'ea', 1.)
        et = double_or_blank(card, istart + 4, 'et', 1.)
        return [sa, st, ea, et]

    def raw_fields(self):
        list_fields = ['PBUSH', self.pid]
        for var in self.vars:
            if var == 'K':
                list_fields += ['K'] + self.Ki
            elif var == 'B':
                list_fields += ['B'] + self.Bi
            elif var == 'GE':
                list_fields += ['GE'] + self.GEi
            elif var == 'RCV':
                list_fields += ['RCV', self.sa, self.st, self.ea, self.et]
            else:
                raise RuntimeError('not supported PBUSH field...')
            nspaces = 8 - (len(list_fields) - 1) % 8

            if nspaces == 8:
                list_fields += [None]
            elif nspaces < 8:
                list_fields += [None] * (nspaces + 1)
        return list_fields

    def repr_fields(self):
        return self.raw_fields()

    def write_card(self, size=8, is_double=False):
        card = self.repr_fields()
        if size == 8:
            return self.comment + print_card_8(card)
        return self.comment + print_card_16(card)


class PBUSH1D(BushingProperty):
    type = 'PBUSH1D'

    #def __init__(self, pid, k, c, m, sa, se, comment=''):
    def __init__(self):
        BushingProperty.__init__(self)
        #if comment:
            #self._comment = comment

        #: Property ID
        #self.pid = pid
        #self.k = k
        #self.c = c
        #self.m = m

        #self.sa = sa
        #self.se = se

        # SPRING parameters
        self.springType = None
        self.springIDT = None
        self.springIDC = None
        self.springIDTDU = None
        self.springIDCDU = None

        # DAMPER parameters
        self.damperType = None
        self.damperIDT = None
        self.damperIDC = None
        self.damperIDTDV = None
        self.damperIDCDV = None

        # GENER parameters
        self.generIDT = None
        self.generIDC = None
        self.generIDTDU = None
        self.generIDCDU = None
        self.generIDTDV = None
        self.generIDCDV = None

        # SHOCK parameters
        self.shockType = None
        self.shockCVT = None
        self.shockCVC = None
        self.shockExpVT = None
        self.shockExpVC = None
        self.shockIDTS = None

        self.shockIDETS = None
        self.shockIDECS = None
        self.shockIDETSD = None
        self.shockIDECSD = None

    #@classmethod
    def add_card(self, card, comment=''):
        if comment:
            self._comment = comment
        pid = integer(card, 1, 'pid')
        k = double_or_blank(card, 2, 'k', 0.0)
        c = double_or_blank(card, 3, 'c', 0.0)
        m = double_or_blank(card, 4, 'm', 0.0)

        sa = double_or_blank(card, 6, 'sa', 0.0)
        se = double_or_blank(card, 7, 'se', 0.0)

        nfields = card.nfields
        self.vars = []
        istart = 9
        while istart < nfields:
            pname = string(card, istart, 'pname')
            if pname == 'SHOCKA':
                istart = self._read_shock(card, istart)
            elif pname == 'SPRING':
                self._read_spring(card, istart)
            elif pname == 'DAMPER':
                self._read_damper(card, istart)
            elif pname == 'GENER':
                self._read_gener(card, istart)
            else:
                break
            istart += 8

        self.pid = pid
        self.k = k
        self.c = c
        self.m = m

        self.sa = sa
        self.se = se
        #return PBUSH1D(pid, k, c, m, sa, se, vars, comment=comment)

    #@classmethod
    #def add_op2_data(cls, data, comment=''):
        #pid = data[0]
        #b = data[1]
        #raise NotImplementedError('PBUSH1D data...')

    def _verify(self, xref=False):
        pid = self.Pid()
        assert isinstance(pid, int), 'pid=%r' % pid

    def _read_shock(self, card, istart):
        """
        F(u, v) = Cv * S(u) * sign(v) * |v|^ev
        """
        self.shockType = string_or_blank(card, istart + 1, 'shockType')
        self.shockCVT = double(card, istart + 2, 'shockCVT')
        self.shockCVC = double_or_blank(card, istart + 3, 'shockCVC')
        self.shockExpVT = double_or_blank(card, istart + 4, 'shockExpVT', 1.0)
        self.shockExpVC = double_or_blank(card, istart + 5,
                                          'shockExpVC', self.shockExpVT)

        if self.shockType == 'TABLE':
            pass
            # self.shockIDTS = integer(card, istart + 6, 'shockIDTS')
            # self.shockIDETS = blank(card, istart + 9, 'shockIDETS')
            # self.shockIDECS = blank(card, istart + 10, 'shockIDECS')
            # self.shockIDETSD = blank(card, istart + 11, 'shockIDETSD')
            # self.shockIDECSD = blank(card, istart + 12, 'shockIDECSD')
        elif self.shockType == 'EQUAT':
            self.shockIDTS = blank(card, istart + 6, 'shockIDTS')
            self.shockIDETS = integer(card, istart + 9, 'shockIDETS')
            self.shockIDECS = integer_or_blank(card, istart + 10,
                                               'shockIDECS', self.shockIDETS)
            self.shockIDETSD = integer(card, istart + 11, 'shockIDETSD')
            self.shockIDECSD = integer_or_blank(card, istart + 11,
                                                'shockIDECSD', self.shockIDETSD)

            #def DEquation(self):
                #if isinstance(self.dequation, int):
                    #return self.dequation
                #return self.dequation.equation_id
        else:
            raise RuntimeError('Invalid shockType=%r on card\n%s' %(self.shockType, card))
        istart += 8
        return istart

    def _read_spring(self, card, istart):
        """
        F(u) = Ft(u)
        """
        self.springType = string_or_blank(card, istart + 1, 'springType')
        self.springIDT = integer(card, istart + 2, 'springIDT')

        if self.springType == 'TABLE':
            self.springIDC = blank(card, istart + 3, 'springIDC')
            self.springIDTDU = blank(card, istart + 4, 'springIDTDU')
            self.springIDCDU = blank(card, istart + 5, 'springIDCDU')
        elif self.springType == 'EQUAT':
            self.springIDC = integer_or_blank(card, istart + 3,
                                              'springIDC', self.springIDT)
            self.springIDTDU = integer(card, istart + 4, 'springIDTDU')
            self.springIDCDU = integer_or_blank(card, istart + 5,
                                                'springIDCDU', self.springIDTDU)
        else:
            raise RuntimeError('Invalid springType=|%s| on card\n%s' %(self.springType, card))

        self.vars.append('SPRING')

    def _read_damper(self, card, istart):
        """
        F(v) = Ft(u)
        """
        self.damperType = string_or_blank(card, istart + 1, 'damperType')
        self.damperIDT = integer(card, istart + 2, 'damperIDT')
        if self.damperType == 'TABLE':
            self.damperIDC = blank(card, istart + 3, 'damperIDC')
            self.damperIDTDV = blank(card, istart + 4, 'damperIDTDV')
            self.damperIDCDV = blank(card, istart + 5, 'damperIDCDV')
        elif self.damperType == 'EQUAT':
            self.damperIDC = integer_or_blank(card, istart + 3, 'damperIDC')
            self.damperIDTDV = integer(card, istart + 4, 'damperIDTDV')
            self.damperIDCDV = integer_or_blank(card, istart + 5, 'damperIDCDV', self.damperIDTDV)
        else:
            raise RuntimeError('Invalid springType=|%s| on card\n%s' %(self.springType, card))

        self.vars.append('DAMPER')

    def _read_gener(self, card, istart):
        """
        F(u, v) = Ft(u, v)
        """
        self.generIDT = integer(card, istart + 2, 'generIDT')
        self.generIDC = integer_or_blank(card, istart + 3,
                                         'generIDC', self.generIDT)
        self.generIDTDU = integer(card, istart + 4, 'generIDTDU')
        self.generIDCDU = integer_or_blank(card, istart + 5,
                                           'generIDCDU', self.generIDTDU)
        self.generIDTDV = integer(card, istart + 6, 'generIDTDV')
        self.generIDCDV = integer_or_blank(card, istart + 7,
                                           'generIDCDV', self.generIDTDV)
        self.vars.append('GENER')

    def _shock_fields(self):
        list_fields = ['SHOCKA', self.shockType, self.shockCVT, self.shockCVC,
                       self.shockExpVT, self.shockExpVC, self.shockIDTS, None, None,
                       self.shockIDETS, self.shockIDECS, self.shockIDETSD,
                       self.shockIDECSD]
        return list_fields

    def _spring_fields(self):
        list_fields = ['SPRING', self.springType, self.springIDT,
                       self.springIDC, self.springIDTDU, self.springIDCDU]
        return list_fields

    def _damper_fields(self):
        list_fields = ['DAMPER', self.damperType, self.damperIDT,
                       self.damperIDC, self.damperIDTDV, self.damperIDCDV]
        return list_fields

    def _gener_fields(self):
        list_fields = ['GENER', None, self.generIDT, self.generIDC,
                       self.generIDTDU, self.generIDCDU, self.generIDTDV,
                       self.generIDCDV]
        return list_fields

    def raw_fields(self):
        list_fields = ['PBUSH1D', self.pid, self.k, self.c, self.m, None,
                       self.sa, self.se, None]
        for var in self.vars:
            if var == 'SHOCKA':
                list_fields += self._shock_fields()
            elif var == 'SPRING':
                list_fields += self._spring_fields()
            elif var == 'DAMPER':
                list_fields += self._damper_fields()
            elif var == 'GENER':
                list_fields += self._gener_fields()
            else:
                msg = 'var=%s not supported PBUSH1D field...' % var
                raise RuntimeError(msg)
            nspaces = 8 - (len(list_fields) - 1) % 8

            if nspaces < 8:
                list_fields += [None] * (nspaces)
        return list_fields

    def repr_fields(self):
        return self.raw_fields()

    def write_card(self, size=8, is_double=False):
        card = self.repr_fields()
        if size == 8:
            return self.comment + print_card_8(card)
        return self.comment + print_card_16(card)


class PBUSH2D(BushingProperty):
    type = 'PBUSH2D'

    def __init__(self, card=None, data=None, comment=''):
        BushingProperty.__init__(self, card, data)
        if comment:
            self._comment = comment
        if card:
            raise NotImplementedError()
        else:
            raise NotImplementedError()

    def write_card(self, size=8, is_double=False):
        """
        Writes the card with the specified width and precision

        Parameters
        ----------
        size : int (default=8)
            size of the field; {8, 16}
        is_double : bool (default=False)
            is this card double precision

        Returns
        -------
        msg : str
            the string representation of the card
        """
        card = self.repr_fields()
        if size == 8:
            return self.comment + print_card_8(card)
        return self.comment + print_card_16(card)


class PBUSHT(BushingProperty):
    type = 'PBUSHT'

    def __init__(self, pid, k_tables, b_tables,
                 ge_tables, kn_tables, comment=''):
        BushingProperty.__init__(self)
        if comment:
            self._comment = comment
        self.pid = pid
        self.k_tables = k_tables
        self.b_tables = b_tables
        self.ge_tables = ge_tables
        self.kn_tables = kn_tables

    @classmethod
    def add_card(cls, card, comment=''):
        k_tables = []
        b_tables = []
        ge_tables = []
        kn_tables = []

        pid = integer(card, 1, 'pid')
        nfields = len(card) - 1
        nrows = nfields // 8
        if nfields % 8 != 0:
            nrows += 1

        for irow in range(nrows):
            ifield = 1 + irow * 8
            param = string(card, ifield + 1, 'param_type')
            table = []
            for j in range(6):
                table_value = integer_or_blank(card, ifield + j + 2, param + '%i' % (j+1))
                table.append(table_value)
            if param == 'K':
                k_tables = table
            elif param == 'B':
                b_tables = table
            elif param == 'GE':
                ge_tables = table
            elif param == 'KN':
                kn_tables = table
            else:
                raise ValueError(param)
        return PBUSHT(pid, k_tables, b_tables, ge_tables, kn_tables,
                      comment=comment)

    def raw_fields(self):
        list_fields = ['PBUSHT', self.pid]
        if self.k_tables:
            list_fields += ['K'] + self.k_tables + [None]
        if self.b_tables:
            list_fields += ['B'] + self.b_tables + [None]
        if self.ge_tables:
            list_fields += ['GE'] + self.ge_tables + [None]
        if self.kn_tables:
            list_fields += ['KN'] + self.kn_tables
        return list_fields

    def write_card(self, size=8, is_double=False):
        card = self.repr_fields()
        if size == 8:
            return self.comment + print_card_8(card)
        return self.comment + print_card_16(card)
