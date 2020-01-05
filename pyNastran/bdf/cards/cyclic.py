"""
defines:
 * CYJOIN

All cards are BaseCard objects.

"""
from __future__ import annotations
from typing import TYPE_CHECKING
#import numpy as np

#from pyNastran.utils.numpy_utils import integer_types
#from pyNastran.bdf.field_writer_8 import set_blank_if_default
from pyNastran.bdf.cards.base_card import BaseCard, expand_thru
from pyNastran.bdf.bdf_interface.assign_type import (
    integer, # integer_or_blank, double, double_or_blank,
    #string_or_blank, blank,
    fields,
    #components_or_blank,
    #integer_string_or_blank, integer_or_double, #parse_components,
    #modal_components_or_blank,
    string, integer_or_string,
)
from pyNastran.bdf.field_writer_8 import print_card_8
from pyNastran.bdf.field_writer_16 import print_card_16
if TYPE_CHECKING:  # pragma: no cover
    from pyNastran.bdf.bdf import BDF


class CYAX(BaseCard):
    """
    Lists grid points that lie on the axis of symmetry in cyclic symmetry analysis.
    CYAX G1 G2 G3 G4 G5 G6 G7 G8
         G9 G10 -etc.-
    """
    type = 'CYAX'

    @classmethod
    def _init_from_empty(cls):
        nids = [1]
        return CYAX(nids, comment='')

    def __init__(self, nids, comment=''):
        if comment:
            self.comment = comment
        self.nids = nids

    @classmethod
    def add_card(cls, card, comment=''):
        """
        Adds a AXIF card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card

        """
        nids = fields(integer, card, 'n', i=1, j=len(card))
        return CYAX(nids, comment=comment)

    def raw_fields(self):
        list_fields = ['CYAX'] + list(self.nids)
        return list_fields

    def write_card(self, size: int=8, is_double: bool=False) -> str:
        card = self.repr_fields()
        msg = self.comment + print_card_8(card)
        return msg


class CYJOIN(BaseCard):
    """
    CYJOIN SIDE C G1 G2 G3 G4 G5 G6
        G7 G8 G9 -etc.-
    """
    type = 'CYJOIN'

    @classmethod
    def _init_from_empty(cls):
        side = 1
        coord = 'R'
        nids = [1, 2]
        return CYJOIN(side, coord, nids)

    def __init__(self, side, coord, nids, comment=''):
        BaseCard.__init__(self)
        if comment:
            self.comment = comment
        self.side = side
        self.coord = coord
        self.nids = expand_thru(nids)
        assert coord in ['T1', 'T2', 'T3', 'R', 'C', 'S'], coord

    def validate(self):
        pass
        #assert len(self.grids1) > 0, 'ngrids1=%s\n%s' % (len(self.grids1), str(self))

    def cross_reference(self, model: BDF) -> None:
        pass

    @classmethod
    def add_card(cls, card, comment=''):
        """
        Adds a CYJOIN card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card

        """
        side = integer(card, 1, 'side')
        coord = string(card, 2, 'coord')
        nids = fields(integer_or_string, card, 'ID', i=3, j=len(card))
        return CYJOIN(side, coord, nids, comment=comment)

    def raw_fields(self):
        #msg = self.comment
        self.nids.sort()
        nids = list(self.nids)
        list_fields = ['CYJOIN', self.side, self.coord] + nids
        return list_fields

    def write_card(self, size: int=8, is_double: bool=False) -> str:
        # double precision?
        card = self.repr_fields()
        if size == 8:
            return self.comment + print_card_8(card)
        return self.comment + print_card_16(card)
