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

    @classmethod
    def add_op2_data(cls, data, comment=''):
        (pid, k1, k2, k3, k4, k5, k6, b1, b2, b3, b4, b5, b6,
         g1, g2, g3, g4, g5, g6, sa, st, ea, et) = data
        k_fields = [k1, k2, k3, k4, k5, k6]
        b_fields = [b1, b2, b3, b4, b5, b6]
        ge_fields = [g1, g2, g3, g4, g5, g6]
        rcv_fields = [sa, st, ea, et]
        return PBUSH(pid, k_fields, b_fields, ge_fields, rcv_fields,
                     comment=comment)

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

    def __init__(self, pid, k, c, m, sa, se, optional_vars, comment=''):
        BushingProperty.__init__(self)
        if comment:
            self._comment = comment

        #: Property ID
        self.pid = pid
        self.k = k
        self.c = c
        self.m = m

        self.sa = sa
        self.se = se

        # SPRING parameters
        #self.spring_type = None
        #self.spring_idt = None
        #self.spring_idc = None
        #self.spring_idtdu = None
        #self.spring_idcdu = None

        # DAMPER parameters
        #self.damper_type = None
        #self.damper_idt = None
        #self.damper_idc = None
        #self.damper_idtdv = None
        #self.damper_idcdv = None

        # GENER parameters
        #self.gener_idt = None
        #self.gener_idc = None
        #self.gener_idtdu = None
        #self.gener_idcdu = None
        #self.gener_idtdv = None
        #self.gener_idcdv = None

        # SHOCK parameters
        #self.shockType = None
        #self.shockCVT = None
        #self.shockCVC = None
        #self.shockExpVT = None
        #self.shockExpVC = None
        #self.shockIDTS = None

        #self.shockIDETS = None
        #self.shockIDECS = None
        #self.shockIDETSD = None
        #self.shockIDECSD = None
        if 'SHOCKA' in optional_vars:
            (shock_type, shock_cvt, shock_cvc, shock_exp_vt, shock_exp_vc,
             shock_idts, shock_idets, shock_idecs, shock_idetsd, shock_idecsd
            ) = optional_vars['SHOCKA']
            self.shock_type = shock_type
            self.shock_cvt = shock_cvt
            self.shock_cvc = shock_cvc
            self.shock_exp_vt = shock_exp_vt
            self.shock_exp_vc = shock_exp_vc
            self.shock_idts = shock_idts
            self.shock_idets = shock_idets
            self.shock_idecs = shock_idecs
            self.shock_idetsd = shock_idetsd
            self.shock_idecsd = shock_idecsd

        if 'SPRING' in optional_vars:
            (spring_type, spring_idt, spring_idc, spring_idtdu,
             spring_idcdu) = optional_vars['SPRING']
            self.spring_type = spring_type
            self.spring_idt = spring_idt
            self.spring_idtc = spring_idc
            self.spring_idc = spring_idc
            self.spring_idtdu = spring_idtdu
            self.spring_idcdu = spring_idcdu

        if 'DAMPER' in optional_vars:
            (damper_type, damper_idt, damper_idc, damper_idtdv,
             damper_idcdv) = optional_vars['DAMPER']
            self.damper_type = damper_type
            self.damper_idt = damper_idt
            self.damper_idc = damper_idc
            self.damper_idtdv = damper_idtdv
            self.damper_idcdv = damper_idcdv

        if 'GENER' in optional_vars:
            (gener_idt, gener_idc,
             gener_idtdu, gener_idcdu,
             gener_idtdv, gener_idcdv
             ) = optional_vars['GENER']

            self.gener_idt = gener_idt
            self.gener_idc = gener_idc
            self.gener_idtdu = gener_idtdu
            self.gener_idcdu = gener_idcdu
            self.gener_idtdv = gener_idtdv
            self.gener_idcdv = gener_idcdv
        self.vars = list(optional_vars.keys())
        self.vars.sort()

    @classmethod
    def add_card(cls, card, comment=''):
        pid = integer(card, 1, 'pid')
        k = double_or_blank(card, 2, 'k', 0.0)
        c = double_or_blank(card, 3, 'c', 0.0)
        m = double_or_blank(card, 4, 'm', 0.0)

        sa = double_or_blank(card, 6, 'sa', 0.0)
        se = double_or_blank(card, 7, 'se', 0.0)

        nfields = card.nfields
        optional_vars = {}
        istart = 9
        while istart < nfields:
            pname = string(card, istart, 'pname')
            if pname == 'SHOCKA':
                istart, out = cls._read_shock(card, istart)
                optional_vars['SHOCKA'] = out

            elif pname == 'SPRING':
                out = cls._read_spring(card, istart)
                optional_vars['SPRING'] = out

            elif pname == 'DAMPER':
                out = cls._read_damper(card, istart)
                optional_vars['DAMPER'] = out

            elif pname == 'GENER':
                out = cls._read_gener(card, istart)
                optional_vars['GENER'] = out
            else:
                break
            istart += 8

        #self.pid = pid
        #self.k = k
        #self.c = c
        #self.m = m

        #self.sa = sa
        #self.se = se
        return PBUSH1D(pid, k, c, m, sa, se, optional_vars, comment=comment)

    #@classmethod
    #def add_op2_data(cls, data, comment=''):
        #pid = data[0]
        #b = data[1]
        #raise NotImplementedError('PBUSH1D data...')

    def _verify(self, xref=False):
        pid = self.Pid()
        assert isinstance(pid, int), 'pid=%r' % pid

    @staticmethod
    def _read_shock(card, istart):
        """
        F(u, v) = Cv * S(u) * sign(v) * |v|^ev
        """
        shock_type = string_or_blank(card, istart + 1, 'shockType')
        shock_cvt = double(card, istart + 2, 'shockCVT')
        shock_cvc = double_or_blank(card, istart + 3, 'shockCVC')
        shock_exp_vt = double_or_blank(card, istart + 4, 'shockExpVT', 1.0)
        shock_exp_vc = double_or_blank(card, istart + 5,
                                       'shockExpVC', shock_exp_vt)

        if shock_type == 'TABLE':
            shock_idts = None
            shock_idets = None
            shock_idecs = None
            shock_idetsd = None
            shock_idecsd = None
            # shock_idts = integer(card, istart + 6, 'shockIDTS')
            # shock_idets = blank(card, istart + 9, 'shockIDETS')
            # shock_idecs = blank(card, istart + 10, 'shockIDECS')
            # shock_idetsd = blank(card, istart + 11, 'shockIDETSD')
            # shock_idecsd = blank(card, istart + 12, 'shockIDECSD')
        elif shock_type == 'EQUAT':
            shock_idts = blank(card, istart + 6, 'shockIDTS')
            shock_idets = integer(card, istart + 9, 'shockIDETS')
            shock_idecs = integer_or_blank(card, istart + 10,
                                           'shockIDECS', shock_idets)
            shock_idetsd = integer(card, istart + 11, 'shockIDETSD')
            shock_idecsd = integer_or_blank(card, istart + 11,
                                            'shockIDECSD', shock_idetsd)
            #def DEquation(self):
                #if isinstance(self.dequation, int):
                    #return self.dequation
                #return self.dequation.equation_id
        else:
            msg = 'Invalid shockType=%r on card\n%s' % (shock_type, card)
            raise RuntimeError(msg)

        out = (
            shock_type, shock_cvt, shock_cvc, shock_exp_vt, shock_exp_vc,
            shock_idts, shock_idets, shock_idecs, shock_idetsd, shock_idecsd
        )
        istart += 8
        return istart, out

    @staticmethod
    def _read_spring(card, istart):
        """
        F(u) = Ft(u)
        """
        spring_type = string_or_blank(card, istart + 1, 'springType')
        spring_idt = integer(card, istart + 2, 'springIDT')

        if spring_type == 'TABLE':
            spring_idc = blank(card, istart + 3, 'springIDC')
            spring_idtdu = blank(card, istart + 4, 'springIDTDU')
            spring_idcdu = blank(card, istart + 5, 'springIDCDU')
        elif spring_type == 'EQUAT':
            spring_idc = integer_or_blank(card, istart + 3,
                                          'springIDC', spring_idt)
            spring_idtdu = integer(card, istart + 4, 'springIDTDU')
            spring_idcdu = integer_or_blank(card, istart + 5,
                                            'springIDCDU', spring_idtdu)
        else:
            msg = 'Invalid springType=%r on card\n%s' % (spring_type, card)
            raise RuntimeError(msg)

        return spring_type, spring_idt, spring_idc, spring_idtdu, spring_idcdu

    @staticmethod
    def _read_damper(card, istart):
        """
        F(v) = Ft(u)
        """
        damper_type = string_or_blank(card, istart + 1, 'damperType')
        damper_idt = integer(card, istart + 2, 'damperIDT')
        if damper_type == 'TABLE':
            damper_idc = blank(card, istart + 3, 'damperIDC')
            damper_idtdv = blank(card, istart + 4, 'damperIDTDV')
            damper_idcdv = blank(card, istart + 5, 'damperIDCDV')
        elif damper_type == 'EQUAT':
            damper_idc = integer_or_blank(card, istart + 3, 'damperIDC')
            damper_idtdv = integer(card, istart + 4, 'damperIDTDV')
            damper_idcdv = integer_or_blank(card, istart + 5, 'damperIDCDV', damper_idtdv)
        else:
            msg = 'Invalid damper_type=%r on card\n%s' % (damper_type, card)
            raise RuntimeError(msg)

        return damper_type, damper_idt, damper_idc, damper_idtdv, damper_idcdv

    @staticmethod
    def _read_gener(card, istart):
        """
        F(u, v) = Ft(u, v)
        """
        gener_idt = integer(card, istart + 2, 'generIDT')
        gener_idc = integer_or_blank(card, istart + 3,
                                    'generIDC', gener_idt)
        gener_idtdu = integer(card, istart + 4, 'generIDTDU')
        gener_idcdu = integer_or_blank(card, istart + 5,
                                      'generIDCDU', gener_idtdu)
        gener_idtdv = integer(card, istart + 6, 'generIDTDV')
        gener_idcdv = integer_or_blank(card, istart + 7,
                                       'generIDCDV', gener_idtdv)
        out = (
            gener_idt, gener_idc,
            gener_idtdu, gener_idcdu,
            gener_idtdv, gener_idcdv
        )
        return out

    def _shock_fields(self):
        list_fields = ['SHOCKA', self.shock_type, self.shock_cvt, self.shock_cvc,
                       self.shock_exp_vt, self.shock_exp_vc, self.shock_idts, None, None,
                       self.shock_idets, self.shock_idecs, self.shock_idetsd,
                       self.shock_idecsd]
        return list_fields

    def _spring_fields(self):
        list_fields = ['SPRING', self.spring_type, self.spring_idt,
                       self.spring_idc, self.spring_idtdu, self.spring_idcdu]
        return list_fields

    def _damper_fields(self):
        list_fields = ['DAMPER', self.damper_type, self.damper_idt,
                       self.damper_idc, self.damper_idtdv, self.damper_idcdv]
        return list_fields

    def _gener_fields(self):
        list_fields = ['GENER', None, self.gener_idt, self.gener_idc,
                       self.gener_idtdu, self.gener_idcdu, self.gener_idtdv,
                       self.gener_idcdv]
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
