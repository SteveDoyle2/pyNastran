# pylint: disable=C0103,R0902,R0904,R0914,C0111
"""
All damper properties are defined in this file.  This includes:

 *   PDAMP
 *   PDAMP5 (not implemented)
 *   PDAMPT
 *   PVISC

All damper properties are DamperProperty and Property objects.
"""
from __future__ import (nested_scopes, generators, division, absolute_import,
                        print_function, unicode_literals)
from six import integer_types

from pyNastran.bdf.field_writer_8 import set_blank_if_default, print_card_8
from pyNastran.bdf.field_writer_16 import print_card_16
from pyNastran.bdf.cards.baseCard import Property
from pyNastran.bdf.bdfInterface.assign_type import (integer, integer_or_blank,
    double, double_or_blank)


class DamperProperty(Property):
    def __init__(self, card, data):
        Property.__init__(self, card, data)

    def cross_reference(self, model):
        pass


class PDAMP(DamperProperty):
    type = 'PDAMP'
    _field_map = {
        1: 'pid', 2:'b',
    }

    def __init__(self, card=None, nPDAMP=0, data=None, comment=''):
        DamperProperty.__init__(self, card, data)
        if comment:
            self._comment = comment
        nOffset = nPDAMP * 2
        if card:
            # 3 PDAMP properties can be defined on 1 PDAMP card
            #: Property ID
            self.pid = integer(card, 1 + nOffset, 'pid')

            # these are split into 2 separate cards
            #: Force per unit velocity (Real)
            self.b = double(card, 2 + nOffset, 'b')
        else:
            self.pid = data[0]
            self.b = data[1]

    def B(self):
        return self.b

    def _verify(self, xref=True):
        pid = self.Pid()
        b = self.B()
        assert isinstance(pid, int), 'pid=%r\n%s' % (pid, str(self))
        assert isinstance(b, float), 'b=%r\n%s' % (b, str(self))

    def raw_fields(self):
        list_fields = ['PDAMP', self.pid, self.b]
        return list_fields

    def repr_fields(self):
        return self.raw_fields()

    def write_card(self, size=8, is_double=False):
        card = self.repr_fields()
        if size == 8:
            return self.comment() + print_card_8(card)
        return self.comment() + print_card_16(card)


class PDAMP5(DamperProperty):
    type = 'PDAMP5'
    _field_map = {
        1: 'pid', 2:'mid', 3:'b',
    }

    def __init__(self, card=None, data=None, comment=''):
        """
        Defines the damping multiplier and references the material properties
        for damping. CDAMP5 is intended for heat transfer analysis only.
        """
        DamperProperty.__init__(self, card, data)
        if comment:
            self._comment = comment
        if card:
            #: Property ID
            self.pid = integer(card, 1, 'pid')
            #: Material ID
            self.mid = integer(card, 2, 'mid')
            #: Damping multiplier. (Real > 0.0)
            #: B is the mass that multiplies the heat capacity CP on the MAT4
            #: or MAT5 entry.
            self.b = double(card, 3, 'b')
            assert len(card) == 4, 'len(PDAMP5 card) = %i' % len(card)
        else:
            self.pid = data[0]
            self.mid = data[1]
            self.b = data[2]

    def _verify(self, xref=True):
        pid = self.Pid()
        #b = self.B()
        assert isinstance(pid, integer_types), 'pid=%r\n%s' % (pid, str(self))
        #assert isinstance(b, float), 'b=%r\n%s' % (b, str(self))

    def cross_reference(self, model):
        self.mid = model.Material(self.mid)

    def Mid(self):
        if isinstance(self.mid, integer_types):
            return self.mid
        return self.mid.mid

    def repr_fields(self):
        return self.raw_fields()

    def raw_fields(self):
        list_fields = ['PDAMP5', self.pid, self.Mid(), self.b]
        return list_fields

    def write_card(self, size=8, is_double=False):
        card = self.repr_fields()
        if size == 8:
            return self.comment() + print_card_8(card)
        return self.comment() + print_card_16(card)


class PDAMPT(DamperProperty):
    type = 'PDAMPT'
    _field_map = {
        1: 'pid', 2:'tbid',
    }

    def __init__(self, card=None, data=None, comment=''):
        DamperProperty.__init__(self, card, data)
        if comment:
            self._comment = comment
        if card:
            #: Property ID
            self.pid = integer(card, 1, 'pid')
            #: Identification number of a TABLEDi entry that defines the
            #: damping force per-unit velocity versus frequency relationship
            self.tbid = integer_or_blank(card, 2, 'tbid', 0)
            assert len(card) <= 3, 'len(PDAMPT card) = %i' % len(card)
        else:
            self.pid = data[0]
            self.tbid = data[1]

    def _verify(self, xref=False):
        pid = self.Pid()
        assert isinstance(pid, int), 'pid=%r' % pid

    def cross_reference(self, model):
        self.tbid = model.Table(self.tbid)

    def Tbid(self):
        if self.tbid == 0:
            return None
        elif isinstance(self.tbid, integer_types):
            return self.tbid
        return self.tbid.tid

    def repr_fields(self):
        return self.raw_fields()

    def raw_fields(self):
        list_fields = ['PDAMPT', self.pid, self.Tbid()]
        return list_fields

    def write_card(self, size=8, is_double=False):
        card = self.repr_fields()
        if size == 8:
            return self.comment() + print_card_8(card)
        return self.comment() + print_card_16(card)


class PVISC(DamperProperty):
    type = 'PVISC'
    _field_map = {
        1: 'pid', 2:'ce', 3:'cr',
    }

    def __init__(self, card=None, nPVISC=0, data=None, comment=''):
        DamperProperty.__init__(self, card, data)
        if comment:
            self._comment = comment
        if card:
            self.pid = integer(card, 1 + 4 * nPVISC, 'pid')
            self.ce = double(card, 2 + 4 * nPVISC, 'ce')
            self.cr = double_or_blank(card, 3 + 4 * nPVISC, 'cr', 0.)
        else:
            self.pid = data[0]
            self.ce = data[1]
            self.cr = data[2]

    def cross_reference(self, model):
        pass

    def _verify(self, xref=False):
        pid = self.Pid()
        assert isinstance(pid, int), 'pid=%r' % pid

    def raw_fields(self):
        list_fields = ['PVISC', self.pid, self.ce, self.cr]
        return list_fields

    def repr_fields(self):
        cr = set_blank_if_default(self.cr, 0.)
        list_fields = ['PVISC', self.pid, self.ce, cr]
        return list_fields

    def write_card(self, size=8, is_double=False):
        card = self.repr_fields()
        if size == 8:
            return self.comment() + print_card_8(card)
        return self.comment() + print_card_16(card)
