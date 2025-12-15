# coding: utf-8
# pylint: disable=W0212,C0103
"""
All ZONA aero cards are defined in this file.  This includes:
 * SET1
 * SETADD

All cards are BaseCard objects.

"""
from __future__ import annotations
from typing import TYPE_CHECKING

# from pyNastran.bdf.cards.aero import zona
from pyNastran.bdf.field_writer_8 import print_card_8
from pyNastran.bdf.cards.base_card import BaseCard
from pyNastran.bdf.bdf_interface.assign_type import (
    integer,
)

if TYPE_CHECKING:  # pragma: no cover
    from pyNastran.bdf.bdf import BDF
    from pyNastran.bdf.bdf_interface.bdf_card import BDFCard


class SETADD(BaseCard):
    type = 'SETADD'
    # _field_map = {
    #     1: 'sid', 2: 'mach', 3: 'q', 8: 'aeqr',
    # }

    def __init__(self, setadd_id: int, ids: list[int], comment: str=''):
        BaseCard.__init__(self)
        if comment:
            self.comment = comment

        self.setadd_id = setadd_id
        self.ids = ids

    @classmethod
    def add_card(cls, card: BDFCard, comment: str=''):
        """
        Adds a GAINSET card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card

        """
        # GAINSET SID TF1 TF2 TF3 TF4
        # GAINSET 10  1   2   3
        setadd_id = integer(card, 1, 'setadd_id')

        # CJUNCT, MIMOSS, SISOTF
        ids = []
        for ifield in range(2, len(card)):
            idi = integer(card, ifield, 'OTFID, id')
            ids.append(idi)
        assert len(card) >= 3, f'len(SETADD card) = {len(card):d}\ncard={card}'
        return SETADD(setadd_id, ids, comment=comment)

    # def validate(self):
    #     assert self.true_g in ['TRUE', 'G'], 'true_g=%r' % self.true_g

    def cross_reference_set(self, model: BDF, method: str,
                            msg: str='') -> None:
        assert method == 'Node'
        ids_ref = []
        log = model.log
        zona = model.zona
        for idi in self.ids:
            if idi in model.sets:
                id_ref = model.sets[idi]
            else:
                sets = list(model.sets)
                sets.sort()
                msg = (
                    f'SETADD={self.setadd_id}{msg}: id={idi} is not [SETADD]\n'
                    f' - sets = {sets}\n'
                )
                log.warning(msg)
                id_ref = None
                # raise RuntimeError(msg)
            ids_ref.append(id_ref)
        self.ids_ref = ids_ref

    def safe_cross_reference(self, model: BDF, xref_errors):
        self.cross_reference(model)

    def uncross_reference(self) -> None:
        """Removes cross-reference links"""
        self.ids_ref = None

    def raw_fields(self):
        """
        Gets the fields in their unmodified form

        Returns
        -------
        fields : list[varies]
            the fields that define the card

        """
        list_fields = ['SETADD', self.setadd_id] + self.ids
        return list_fields

    def repr_fields(self):
        list_fields = self.raw_fields()
        return list_fields

    def write_card(self, size: int = 8, is_double: bool = False) -> str:
        card = self.repr_fields()
        return self.comment + print_card_8(card)
