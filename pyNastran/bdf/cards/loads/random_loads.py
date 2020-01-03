# pylint: disable=R0902,R0904,R0914,W0231,R0201
"""
All static loads are defined in this file.  This includes:

 * LSEQ
 * DAREA
 * SLOAD
 * RFORCE
 * RANDPS

"""
from __future__ import annotations
from typing import TYPE_CHECKING

#from pyNastran.bdf.errors import CrossReferenceError
from pyNastran.bdf.cards.base_card import BaseCard
from pyNastran.bdf.bdf_interface.assign_type import (
    integer, integer_or_blank, double, double_or_blank,)
from pyNastran.bdf.field_writer_8 import print_card_8
from pyNastran.bdf.field_writer_16 import print_card_16
from pyNastran.bdf.field_writer_double import print_card_double
if TYPE_CHECKING:  # pragma: no cover
    from pyNastran.bdf.bdf import BDF


class RandomLoad(BaseCard):
    def __init__(self, card, data):
        pass


class RANDPS(RandomLoad):
    r"""
    Power Spectral Density Specification

    Defines load set power spectral density factors for use in random analysis
    having the frequency dependent form:

    .. math:: S_{jk}(F) = (X+iY)G(F)
    """
    type = 'RANDPS'

    @classmethod
    def _init_from_empty(cls):
        sid = 1
        j = 2
        k = 3
        return RANDPS(sid, j, k, x=0., y=0., tid=0, comment='')

    def __init__(self, sid, j, k, x=0., y=0., tid=0, comment=''):
        """
        Creates a RANDPS card

        Parameters
        ----------
        sid : int
            random analysis set id
            defined by RANDOM in the case control deck
        j : int
            Subcase id of the excited load set
        k : int
            Subcase id of the applied load set
            k > j
        x / y : float; default=0.0
            Components of the complex number
        tid : int; default=0
            TABRNDi id that defines G(F)
        comment : str; default=''
            a comment for the card

        """
        if comment:
            self.comment = comment
        #: Random analysis set identification number. (Integer > 0)
        #: Defined by RANDOM in the Case Control Deck.
        self.sid = sid

        #: Subcase identification number of the excited load set.
        #: (Integer > 0)
        self.j = j

        #: Subcase identification number of the applied load set.
        #: (Integer >= 0; K >= J)
        self.k = k

        #: Components of the complex number. (Real)
        self.x = x
        self.y = y

        #: Identification number of a TABRNDi entry that defines G(F).
        self.tid = tid
        assert self.sid > 0, 'sid=%s\n%s' % (self.sid, self)
        self.tid_ref = None

    def validate(self):
        assert self.k >= self.j, 'k=%s j=%s\n%s' % (self.k, self.j, self)

    @classmethod
    def add_card(cls, card, comment=''):
        """
        Adds a RANDPS card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card

        """
        sid = integer(card, 1, 'sid')
        j = integer(card, 2, 'j')
        k = integer(card, 3, 'k')
        x = double_or_blank(card, 4, 'x', 0.0)
        y = double_or_blank(card, 5, 'y', 0.0)
        tid = integer_or_blank(card, 6, 'tid', 0)
        assert len(card) <= 7, 'len(RANDPS card) = %i\ncard=%s' % (len(card), card)
        return RANDPS(sid, j, k, x, y, tid, comment=comment)

    def cross_reference(self, model: BDF) -> None:
        """
        Cross links the card so referenced cards can be extracted directly

        Parameters
        ----------
        model : BDF()
            the BDF object

        """
        if self.tid:
            msg = ', which is required by RANDPS sid=%s' % (self.sid)
            #self.tid = model.Table(self.tid, msg=msg)
            self.tid_ref = model.RandomTable(self.tid, msg=msg)

    def safe_cross_reference(self, model, xref_errors):
        return self.cross_reference(model)

    def uncross_reference(self) -> None:
        """Removes cross-reference links"""
        self.tid = self.Tid()
        self.tid_ref = None

    def get_loads(self):
        return [self]

    def Tid(self):
        if self.tid_ref is not None:
            return self.tid_ref.tid
        elif self.tid == 0:
            return None
        else:
            return self.tid

    def raw_fields(self):
        list_fields = ['RANDPS', self.sid, self.j, self.k, self.x, self.y,
                       self.Tid()]
        return list_fields

    def repr_fields(self):
        return self.raw_fields()

    def write_card(self, size: int=8, is_double: bool=False) -> str:
        card = self.repr_fields()
        if size == 8:
            return self.comment + print_card_8(card)
        if is_double:
            return self.comment + print_card_double(card)
        return self.comment + print_card_16(card)

class RANDT1(RandomLoad):
    type = 'RANDT1'

    @classmethod
    def _init_from_empty(cls):
        sid = 1
        n = 10
        t0 = 1
        tmax = 1.
        return RANDT1(sid, n, t0, tmax, comment='')

    def __init__(self, sid, n, t0, tmax, comment=''):
        # type: (int, int, int, float, str) -> None
        """
        Creates a RANDT1 card

        Parameters
        ----------
        sid : int
            random analysis set id
            defined by RANDOM in the case control deck
        n : int
            ???
        t0 : int
            ???
        tmax : float
            ???
        comment : str; default=''
            a comment for the card

        """
        if comment:
            self.comment = comment
        #: Random analysis set identification number. (Integer > 0)
        #: Defined by RANDOM in the Case Control Deck.
        self.sid = sid

        self.n = n
        self.t0 = t0
        self.tmax = tmax

        assert self.sid > 0, 'sid=%s\n%s' % (self.sid, self)

    #def validate(self):
        #assert self.k >= self.j, 'k=%s j=%s\n%s' % (self.k, self.j, self)

    @classmethod
    def add_card(cls, card, comment=''):
        """
        Adds a RANDT1 card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card

        """
        sid = integer(card, 1, 'sid')
        n = integer(card, 2, 'n')
        t0 = double(card, 3, 't0')
        tmax = double(card, 4, 'tmax')
        assert len(card) <= 5, 'len(RANDT1 card) = %i\ncard=%s' % (len(card), card)
        return RANDT1(sid, n, t0, tmax, comment=comment)

    def cross_reference(self, model: BDF) -> None:
        """
        Cross links the card so referenced cards can be extracted directly

        Parameters
        ----------
        model : BDF()
            the BDF object

        """
        pass

    def safe_cross_reference(self, model, xref_errors):
        pass

    def uncross_reference(self) -> None:
        """Removes cross-reference links"""
        pass

    def get_loads(self):
        return [self]

    def raw_fields(self):
        list_fields = ['RANDT1', self.sid, self.n, self.t0, self.tmax]
        return list_fields

    def repr_fields(self):
        return self.raw_fields()

    def write_card(self, size: int=8, is_double: bool=False) -> str:
        card = self.repr_fields()
        if size == 8:
            return self.comment + print_card_8(card)
        if is_double:
            return self.comment + print_card_double(card)
        return self.comment + print_card_16(card)
