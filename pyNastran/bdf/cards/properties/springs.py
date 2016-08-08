# pylint: disable=C0103,R0902,R0904,R0914
"""
All spring properties are defined in this file.  This includes:

 * PELAS
 * PELAST

All spring properties are SpringProperty and Property objects.
"""
from __future__ import (nested_scopes, generators, division, absolute_import,
                        print_function, unicode_literals)

from pyNastran.utils import integer_types
from pyNastran.bdf.field_writer_8 import set_blank_if_default
from pyNastran.bdf.cards.base_card import Property
from pyNastran.bdf.bdf_interface.assign_type import (
    integer, integer_or_blank, double, double_or_blank)
from pyNastran.bdf.field_writer_8 import print_card_8


class SpringProperty(Property):
    def __init__(self):
        Property.__init__(self)

    def write_card(self, size=8, is_double=False):
        card = self.repr_fields()
        return self.comment + print_card_8(card)


class PELAS(SpringProperty):
    """
    Specifies the stiffness, damping coefficient, and stress coefficient of a
    scalar elastic (spring) element (CELAS1 or CELAS3 entry).
    """
    type = 'PELAS'
    _field_map = {
        1: 'pid', 2:'k', 3:'ge', 4:'s',
    }

    def __init__(self, pid, k, ge, s, comment=''):
        if comment:
            self._comment = comment
        SpringProperty.__init__(self)
        # 2 PELAS properties can be defined on 1 PELAS card
        # these are split into 2 separate cards

        #: Property identification number. (Integer > 0)
        self.pid = pid
        #: Ki Elastic property value. (Real)
        self.k = k
        #: Damping coefficient, . See Remarks 5. and 6. (Real)
        #: To obtain the damping coefficient GE, multiply the
        #: critical damping ratio c/c0 by 2.0.
        self.ge = ge
        #: Stress coefficient. (Real)
        self.s = s

    @classmethod
    def add_card(cls, card, icard=0, comment=''):
        noffset = icard * 4
        pid = integer(card, 1 + noffset, 'pid')
        k = double(card, 2 + noffset, 'k')
        ge = double_or_blank(card, 3 + noffset, 'ge', 0.)
        s = double_or_blank(card, 4 + noffset, 's', 0.)
        return PELAS(pid, k, ge, s, comment=comment)

    @classmethod
    def add_op2_data(cls, data, comment=''):
        pid = data[0]
        k = data[1]
        ge = data[2]
        s = data[3]
        return PELAS(pid, k, ge, s, comment=comment)

    def cross_reference(self, model):
        #if self.sol in [108, 129]:
            #self.pid = self.pelasts[self.pid]
        pass

    def uncross_reference(self):
        pass

    def K(self):
        return self.k

    def _verify(self, xref=False):
        pid = self.Pid()
        k = self.K()
        ge = self.ge
        s = self.s
        assert isinstance(pid, int), 'pid=%r' % pid
        assert isinstance(k, float), 'k=%r' % k
        assert isinstance(ge, float), 'ge=%r' % ge
        assert isinstance(s, float), 'ge=%r' % s

    def raw_fields(self):
        list_fields = ['PELAS', self.pid, self.k, self.ge, self.s]
        return list_fields

    def repr_fields(self):
        ge = set_blank_if_default(self.ge, 0.)
        s = set_blank_if_default(self.s, 0.)
        list_fields = ['PELAS', self.pid, self.k, ge, s]
        return list_fields


class PELAST(SpringProperty):
    """
    Frequency Dependent Elastic Property
    Defines the frequency dependent properties for a PELAS Bulk Data entry.

    The PELAST entry is ignored in all solution sequences except frequency
    response (108) or nonlinear analyses (129).
    """
    type = 'PELAST'
    _field_map = {
        1: 'pid', 2:'tkid', 3:'tgeid', 4:'tknid',
    }

    def __init__(self, pid, tkid, tgeid, tknid, comment=''):
        SpringProperty.__init__(self)
        if comment:
            self._comment = comment

        #: Property identification number. (Integer > 0)
        self.pid = pid
        #: Identification number of a TABLEDi entry that defines the
        #: force per unit displacement vs. frequency relationship.
        #: (Integer > 0; Default = 0)
        self.tkid = tkid
        #: Identification number of a TABLEDi entry that defines the
        #: nondimensional structural damping coefficient vs. frequency
        #: relationship. (Integer > 0; Default = 0)
        self.tgeid = tgeid
        #: Identification number of a TABELDi entry that defines the nonlinear
        #: force vs. displacement relationship. (Integer > 0; Default = 0)
        self.tknid = tknid

    @classmethod
    def add_card(cls, card, comment=''):
        pid = integer(card, 1, 'pid')
        tkid = integer_or_blank(card, 2, 'tkid', 0)
        tgeid = integer_or_blank(card, 3, 'tgeid', 0)
        tknid = integer_or_blank(card, 4, 'tknid', 0)
        assert len(card) <= 5, 'len(PELAST card) = %i\ncard=%s' % (len(card), card)
        return PELAST(pid, tkid, tgeid, tknid, comment=comment)

    def cross_reference(self, model):
        """
        Cross links the card so referenced cards can be extracted directly

        Parameters
        ----------
        model : BDF()
            the BDF object
        """
        self.pid = model.Property(self.pid)
        self.pid_ref = self.pid
        if self.tkid > 0:
            self.tkid = model.Table(self.tkid)
            self.tkid_ref = self.tkid
        if self.tgeid > 0:
            self.tgeid = model.Table(self.tgeid)
            self.tgeid_ref = self.tgeid
        if self.tknid > 0:
            self.tknid = model.Table(self.tknid)
            self.tknid_ref = self.tknid

    def uncross_reference(self, model):
        self.pid = self.Pid()
        self.tkid = self.Tkid()
        self.tgeid = self.Tgeid()
        self.tknid = self.Tknid()
        del self.pid_ref, self.tkid_ref, self.tgeid_ref, self.tknid_ref

    def Pid(self):
        if isinstance(self.pid, integer_types):
            return self.pid
        return self.pid_ref.pid

    def Tkid(self):
        """
        Returns the table ID for force per unit displacement vs frequency
        (k=F/d vs freq)
        """
        if isinstance(self.tkid, integer_types):
            return self.tkid
        return self.tkid_ref.tid

    def Tknid(self):
        """
        Returns the table ID for nondimensional force vs. displacement
        """
        if isinstance(self.tknid, integer_types):
            return self.tknid
        return self.tknid_ref.tid

    def Tgeid(self):
        """
        Returns the table ID for nondimensional structural damping
        coefficient vs. frequency (c/c0 vs freq)
        """
        if isinstance(self.tgeid, integer_types):
            return self.tgeid
        return self.tgeid_ref.tid

    def raw_fields(self):
        return ['PELAST', self.pid, self.Tkid(), self.Tgeid(), self.Tknid()]

