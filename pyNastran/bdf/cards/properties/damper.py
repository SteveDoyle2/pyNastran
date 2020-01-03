# pylint: disable=C0103
"""
All damper properties are defined in this file.  This includes:

 *   PDAMP
 *   PDAMP5 (not implemented)
 *   PDAMPT
 *   PVISC

All damper properties are DamperProperty and Property objects.

"""
from __future__ import annotations
from typing import TYPE_CHECKING

from pyNastran.utils.numpy_utils import integer_types
from pyNastran.bdf.field_writer_8 import set_blank_if_default, print_card_8
from pyNastran.bdf.field_writer_16 import print_card_16
from pyNastran.bdf.cards.base_card import Property
from pyNastran.bdf.bdf_interface.assign_type import (
    integer, integer_or_blank, double, double_or_blank)
if TYPE_CHECKING:  # pragma: no cover
    from pyNastran.bdf.bdf import BDF


class DamperProperty(Property):
    def __init__(self):
        Property.__init__(self)

    def cross_reference(self, model: BDF) -> None:
        pass

    def uncross_reference(self) -> None:
        """Removes cross-reference links"""
        pass


class PDAMP(DamperProperty):
    """
    +-------+------+-----+------+----+------+----+------+----+
    |   1   |  2   |  3  |   4  | 5  |  6   |  7 |   8  |  9 |
    +=======+======+=====+======+====+======+====+======+====+
    | PDAMP | PID1 | B1  | PID2 | B2 | PID3 | B3 | PID4 | B4 |
    +-------+------+-----+------+----+------+----+------+----+
    | PDAMP |  1   | 2.0 |      |    |      |    |      |    |
    +-------+------+-----+------+----+------+----+------+----+
    """
    type = 'PDAMP'
    _field_map = {
        1: 'pid', 2:'b',
    }
    pname_fid_map = {
        # 1-based
        3 : 'b',
        'B1' : 'b',
    }

    @classmethod
    def _init_from_empty(cls):
        pid = 1
        b = 1.
        return PDAMP(pid, b, comment='')

    def __init__(self, pid, b, comment=''):
        DamperProperty.__init__(self)
        if comment:
            self.comment = comment
        # 3 PDAMP properties can be defined on 1 PDAMP card
        #: Property ID
        self.pid = pid
        # these are split into 2 separate cards
        #: Force per unit velocity (Real)
        self.b = b

    @classmethod
    def add_card(cls, card, icard=0, comment=''):
        """
        Adds a PDAMP card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        icard : int; default=0
            the index of the card that's being parsed
        comment : str; default=''
            a comment for the card
        """
        noffset = icard * 2
        pid = integer(card, 1 + noffset, 'pid')

        b = double(card, 2 + noffset, 'b')
        return PDAMP(pid, b, comment=comment)

    @classmethod
    def add_op2_data(cls, data, comment=''):
        """
        Adds a PDAMP card from the OP2

        Parameters
        ----------
        data : List[varies]
            a list of fields defined in OP2 format
        comment : str; default=''
            a comment for the card
        """
        pid = data[0]
        b = data[1]
        return PDAMP(pid, b, comment=comment)

    def B(self):
        return self.b

    def _verify(self, xref):
        pid = self.Pid()
        b = self.B()
        assert isinstance(pid, integer_types), 'pid=%r\n%s' % (pid, str(self))
        assert isinstance(b, float), 'b=%r\n%s' % (b, str(self))

    def raw_fields(self):
        list_fields = ['PDAMP', self.pid, self.b]
        return list_fields

    def repr_fields(self):
        return self.raw_fields()

    def write_card(self, size: int=8, is_double: bool=False) -> str:
        card = self.repr_fields()
        if size == 8:
            return self.comment + print_card_8(card)
        return self.comment + print_card_16(card)


class PDAMP5(DamperProperty):
    type = 'PDAMP5'
    _field_map = {
        1: 'pid', 2:'mid', 3:'b',
    }

    @classmethod
    def _init_from_empty(cls):
        pid = 1
        mid = 2
        b = 1.
        return PDAMP5(pid, mid, b, comment='')

    def __init__(self, pid, mid, b, comment=''):
        """
        Defines the damping multiplier and references the material properties
        for damping. CDAMP5 is intended for heat transfer analysis only.

        """
        DamperProperty.__init__(self)
        if comment:
            self.comment = comment
            #: Property ID
        self.pid = pid
        #: Material ID
        self.mid = mid
        #: Damping multiplier. (Real > 0.0)
        #: B is the mass that multiplies the heat capacity CP on the MAT4
        #: or MAT5 entry.
        self.b = b
        self.mid_ref = None

    @classmethod
    def add_card(cls, card, comment=''):
        """
        Adds a PDAMP5 card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card

        """
        pid = integer(card, 1, 'pid')
        mid = integer(card, 2, 'mid')
        b = double(card, 3, 'b')
        assert len(card) == 4, 'len(PDAMP5 card) = %i\ncard=%s' % (len(card), card)
        return PDAMP5(pid, mid, b, comment=comment)

    @classmethod
    def add_op2_data(cls, data, comment=''):
        """
        Adds a PDAMP5 card from the OP2

        Parameters
        ----------
        data : List[varies]
            a list of fields defined in OP2 format
        comment : str; default=''
            a comment for the card

        """
        pid = data[0]
        mid = data[1]
        b = data[2]
        return PDAMP5(pid, mid, b, comment=comment)

    def _verify(self, xref):
        pid = self.Pid()
        #b = self.B()
        assert isinstance(pid, integer_types), 'pid=%r\n%s' % (pid, str(self))
        #assert isinstance(b, float), 'b=%r\n%s' % (b, str(self))

    def cross_reference(self, model: BDF) -> None:
        """
        Cross links the card so referenced cards can be extracted directly

        Parameters
        ----------
        model : BDF()
            the BDF object

        """
        self.mid_ref = model.Material(self.mid)

    def uncross_reference(self) -> None:
        """Removes cross-reference links"""
        self.mid = self.Mid()
        self.mid_ref = None

    def Mid(self):
        if self.mid_ref is not None:
            return self.mid_ref.mid
        return self.mid

    def repr_fields(self):
        return self.raw_fields()

    def raw_fields(self):
        list_fields = ['PDAMP5', self.pid, self.Mid(), self.b]
        return list_fields

    def write_card(self, size: int=8, is_double: bool=False) -> str:
        card = self.repr_fields()
        if size == 8:
            return self.comment + print_card_8(card)
        return self.comment + print_card_16(card)


class PDAMPT(DamperProperty):
    type = 'PDAMPT'
    _field_map = {
        1: 'pid', 2:'tbid',
    }
    pname_fid_map = {
        # 1-based
        #3 : 'b',
        #'B1' : 'b',
    }

    @classmethod
    def _init_from_empty(cls):
        pid = 1
        tbid = 2
        return PDAMPT(pid, tbid, comment='')

    def __init__(self, pid, tbid, comment=''):
        DamperProperty.__init__(self)
        if comment:
            self.comment = comment
        #: Property ID
        self.pid = pid
        #: Identification number of a TABLEDi entry that defines the
        #: damping force per-unit velocity versus frequency relationship
        self.tbid = tbid
        self.tbid_ref = None

    @classmethod
    def add_card(cls, card, comment=''):
        """
        Adds a PDAMPT card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card

        """
        pid = integer(card, 1, 'pid')
        tbid = integer_or_blank(card, 2, 'tbid', 0)
        assert len(card) <= 3, 'len(PDAMPT card) = %i\ncard=%s' % (len(card), card)
        return PDAMPT(pid, tbid, comment=comment)

    @classmethod
    def add_op2_data(cls, data, comment=''):
        """
        Adds a PDAMPT card from the OP2

        Parameters
        ----------
        data : List[varies]
            a list of fields defined in OP2 format
        comment : str; default=''
            a comment for the card

        """
        pid = data[0]
        tbid = data[1]
        return PDAMPT(pid, tbid, comment=comment)

    def _verify(self, xref):
        pid = self.Pid()
        assert isinstance(pid, integer_types), 'pid=%r' % pid

    def cross_reference(self, model: BDF) -> None:
        """
        Cross links the card so referenced cards can be extracted directly

        Parameters
        ----------
        model : BDF()
            the BDF object

        """
        self.tbid_ref = model.TableD(self.tbid)

    def uncross_reference(self) -> None:
        """Removes cross-reference links"""
        self.tbid = self.Tbid()
        self.tbid_ref = None

    def Tbid(self):
        if self.tbid_ref is not None:
            return self.tbid_ref.tid
        elif self.tbid == 0:
            return None
        return self.tbid

    def repr_fields(self):
        return self.raw_fields()

    def raw_fields(self):
        list_fields = ['PDAMPT', self.pid, self.Tbid()]
        return list_fields

    def write_card(self, size: int=8, is_double: bool=False) -> str:
        card = self.repr_fields()
        if size == 8:
            return self.comment + print_card_8(card)
        return self.comment + print_card_16(card)


class PVISC(DamperProperty):
    """
    Viscous Damping Element Property
    Defines properties of a one-dimensional viscous damping element (CVISC entry).

    +-------+------+-----+------+------+-----+-----+
    |   1   |  2   |  3  |  4   |   5  |  6  |  7  |
    +=======+======+=====+======+======+=====+=====+
    | PVISC | PID1 | CE1 | CR1  | PID2 | CE2 | CR2 |
    +-------+------+-----+------+------+-----+-----+
    | PVISC |  3   | 6.2 | 3.94 |      |     |     |
    +-------+------+-----+------+------+-----+-----+
    """
    type = 'PVISC'
    _field_map = {
        1: 'pid', 2:'ce', 3:'cr',
    }
    pname_fid_map = {
        'CE1' : 'ce',
    }

    @classmethod
    def _init_from_empty(cls):
        pid = 1
        ce = 1.
        cr = 1.
        return PVISC(pid, ce, cr, comment='')

    def __init__(self, pid, ce, cr, comment=''):
        """
        Creates a PVISC card

        Parameters
        ----------
        pid : int
            property id for a CVISC
        ce : float
            Viscous damping values for extension in units of force per unit velocity
        cr : float
            Viscous damping values for rotation in units of moment per unit velocity.
        comment : str; default=''
            a comment for the card
        """
        DamperProperty.__init__(self)
        if comment:
            self.comment = comment
        self.pid = pid
        self.ce = ce
        self.cr = cr

    @classmethod
    def add_card(cls, card, icard=0, comment=''):
        """
        Adds a PMASS card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        icard : int; default=0
            the index of the card that's being parsed
        comment : str; default=''
            a comment for the card
        """
        pid = integer(card, 1 + 4 * icard, 'pid')
        ce = double(card, 2 + 4 * icard, 'ce')
        cr = double_or_blank(card, 3 + 4 * icard, 'cr', 0.)
        return PVISC(pid, ce, cr, comment=comment)

    @classmethod
    def add_op2_data(cls, data, comment=''):
        """
        Adds a PVISC card from the OP2

        Parameters
        ----------
        data : List[varies]
            a list of fields defined in OP2 format
        comment : str; default=''
            a comment for the card
        """
        pid = data[0]
        ce = data[1]
        cr = data[2]
        return PVISC(pid, ce, cr, comment=comment)

    def cross_reference(self, model: BDF) -> None:
        pass

    def uncross_reference(self) -> None:
        """Removes cross-reference links"""
        pass

    def _verify(self, xref):
        pid = self.Pid()
        assert isinstance(pid, integer_types), 'pid=%r' % pid

    def raw_fields(self):
        list_fields = ['PVISC', self.pid, self.ce, self.cr]
        return list_fields

    def repr_fields(self):
        cr = set_blank_if_default(self.cr, 0.)
        list_fields = ['PVISC', self.pid, self.ce, cr]
        return list_fields

    def write_card(self, size: int=8, is_double: bool=False) -> str:
        card = self.repr_fields()
        if size == 8:
            return self.comment + print_card_8(card)
        return self.comment + print_card_16(card)
