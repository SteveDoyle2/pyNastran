# pylint: disable=C0103
"""
All spring properties are defined in this file.  This includes:

 * PELAS
 * PELAST

All spring properties are SpringProperty and Property objects.

"""
from __future__ import annotations
from typing import TYPE_CHECKING
from pyNastran.bdf.field_writer_8 import set_blank_if_default
from pyNastran.bdf.cards.base_card import Property
from pyNastran.bdf.bdf_interface.assign_type import (
    integer, integer_or_blank, double, double_or_blank)
from pyNastran.bdf.field_writer_8 import print_card_8
if TYPE_CHECKING:  # pragma: no cover
    from pyNastran.bdf.bdf import BDF


class SpringProperty(Property):
    def __init__(self):
        Property.__init__(self)

    def write_card(self, size: int=8, is_double: bool=False) -> str:
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
    pname_fid_map = {
        #2 : 'k', 'K' : 'k',
        #3 : 'ge', 'GE' : 'ge',
        #4 : 's', 'S' : 's',
        'K1' : 'k',
        'GE1' : 'ge',
    }

    @classmethod
    def _init_from_empty(cls):
        pid = 1
        k = 1.
        return PELAS(pid, k, ge=0., s=0., comment='')

    def __init__(self, pid, k, ge=0., s=0., comment=''):
        """
        Creates a PELAS card

        Parameters
        ----------
        pid : int
            property id
        k : float
            spring stiffness
        ge : int; default=0.0
            damping coefficient
        s : float; default=0.0
            stress coefficient
        comment : str; default=''
            a comment for the card
        """
        if comment:
            self.comment = comment
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

        self.pelast_ref = None

    @classmethod
    def add_card(cls, card, icard=0, comment=''):
        """
        Adds a PELAS card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        icard : int; default=0
            the index of the PELAS card that's being parsed
            must be 0 or 1
        comment : str; default=''
            a comment for the card
        """
        noffset = icard * 4
        pid = integer(card, 1 + noffset, 'pid')
        k = double(card, 2 + noffset, 'k')
        ge = double_or_blank(card, 3 + noffset, 'ge', 0.)
        s = double_or_blank(card, 4 + noffset, 's', 0.)
        return PELAS(pid, k, ge, s, comment=comment)

    @classmethod
    def add_op2_data(cls, data, comment=''):
        """
        Adds a PELAS card from the OP2

        Parameters
        ----------
        data : List[varies]
            a list of fields defined in OP2 format
        comment : str; default=''
            a comment for the card
        """
        pid = data[0]
        k = data[1]
        ge = data[2]
        s = data[3]
        return PELAS(pid, k, ge, s, comment=comment)

    def cross_reference(self, model: BDF) -> None:
        #if model.sol in [108, 129]:
        if self.pid in model.pelast:
            self.pelast_ref = model.pelast[self.pid]

    def uncross_reference(self) -> None:
        """Removes cross-reference links"""
        self.pelast_ref = None

    def K(self):
        return self.k

    def _verify(self, xref):
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
    pname_fid_map = {
    'TKID' : 'tknid',
    }

    @classmethod
    def _init_from_empty(cls):
        pid = 1
        return PELAST(pid, tkid=0, tgeid=0, tknid=0, comment='')

    def __init__(self, pid, tkid=0, tgeid=0, tknid=0, comment=''):
        """
        Creates a PELAST card

        Parameters
        ----------
        pid : int
            property id
        tkid : float
            TABLEDx that defines k vs. frequency
        tgeid : int; default=0
            TABLEDx that defines ge vs. frequency
        s : float; default=0.
            TABLEDx that defines force vs. displacement
        comment : str; default=''
            a comment for the card
        """
        SpringProperty.__init__(self)
        if comment:
            self.comment = comment

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
        self.pid_ref = None
        self.tkid_ref = None
        self.tgeid_ref = None
        self.tknid_ref = None

    @classmethod
    def add_card(cls, card, comment=''):
        """
        Adds a PELAST card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card
        """
        pid = integer(card, 1, 'pid')
        tkid = integer_or_blank(card, 2, 'tkid', 0)
        tgeid = integer_or_blank(card, 3, 'tgeid', 0)
        tknid = integer_or_blank(card, 4, 'tknid', 0)
        assert len(card) <= 5, 'len(PELAST card) = %i\ncard=%s' % (len(card), card)
        return PELAST(pid, tkid, tgeid, tknid, comment=comment)

    @classmethod
    def add_op2_data(cls, data, comment=''):
        """
        Adds a PELAST card from the OP2

        Parameters
        ----------
        data : List[varies]
            a list of fields defined in OP2 format
        comment : str; default=''
            a comment for the card
        """
        (pid, tkid, tgeid, tknid) = data
        return PELAST(pid, tkid, tgeid, tknid, comment=comment)

    def cross_reference(self, model: BDF) -> None:
        """
        Cross links the card so referenced cards can be extracted directly

        Parameters
        ----------
        model : BDF()
            the BDF object
        """
        self.pid_ref = model.Property(self.pid)
        if self.tkid > 0:
            self.tkid_ref = model.TableD(self.tkid)
        if self.tgeid > 0:
            self.tgeid_ref = model.TableD(self.tgeid)
        if self.tknid > 0:
            self.tknid_ref = model.TableD(self.tknid)

    def uncross_reference(self, model):
        self.pid = self.Pid()
        self.tkid = self.Tkid()
        self.tgeid = self.Tgeid()
        self.tknid = self.Tknid()
        self.pid_ref = None
        self.tkid_ref = None
        self.tgeid_ref = None
        self.tknid_ref = None

    def Pid(self):
        if self.pid_ref is None:
            return self.pid
        return self.pid_ref.pid

    def Tkid(self):
        """
        Returns the table ID for force per unit displacement vs frequency
        (k=F/d vs freq)
        """
        if self.tkid_ref is None:
            return self.tkid
        return self.tkid_ref.tid

    def Tknid(self):
        """
        Returns the table ID for nondimensional force vs. displacement
        """
        if self.tknid_ref is None:
            return self.tknid
        return self.tknid_ref.tid

    def Tgeid(self):
        """
        Returns the table ID for nondimensional structural damping
        coefficient vs. frequency (c/c0 vs freq)
        """
        if self.tgeid_ref is None:
            return self.tgeid
        return self.tgeid_ref.tid

    def raw_fields(self):
        return ['PELAST', self.pid, self.Tkid(), self.Tgeid(), self.Tknid()]
