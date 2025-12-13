from __future__ import annotations
from typing import TYPE_CHECKING

from pyNastran.bdf.field_writer_8 import print_card_8
from pyNastran.bdf.cards.base_card import BaseCard
if TYPE_CHECKING:  # pragma: no cover
    from pyNastran.bdf.bdf import BDF
    from pyNastran.bdf.bdf_interface.bdf_card import BDFCard


class DMIL(BaseCard):
    type = 'DMIL'

    # _field_map = {
    #     1: 'sid', 2: 'mach', 3: 'q', 8: 'aeqr',
    # }
    def __init__(self, name: str, col: int, row: int,
                 values: list[float],
                 comment: str=''):
        BaseCard.__init__(self)
        if comment:
            self.comment = comment

        self.name = name
        assert len(name) > 0, name
        self.col = col
        self.row = row
        self.values = values

    @classmethod
    def add_card(cls, card: BDFCard, comment: str=''):
        """
        Adds a TRIM card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card

        """
        fields = []
        field0 = ''
        for i, field in enumerate(card.fields()[1:]):
            if field is not None:
                field0 += field

            if i % 2 == 1:
                # print(f'adding field={field0}')
                fields.append(field0)
                field0 = ''
        if field0:
            fields.append(field0)

        name = fields[0]
        col = int(fields[1])
        row = int(fields[2])
        values = []
        for value_str in fields[3:]:
            value = float(value_str)
            values.append(value)
        return DMIL(name, col, row, values, comment=comment)

    # def validate(self):
    #     assert self.true_g in ['TRUE', 'G'], 'true_g=%r' % self.true_g

    def cross_reference(self, model: BDF) -> None:
        pass

    def safe_cross_reference(self, model: BDF):
        pass

    def uncross_reference(self) -> None:
        """Removes cross-reference links"""
        pass

    def raw_fields(self):
        """
        Gets the fields in their unmodified form

        Returns
        -------
        fields : list[varies]
            the fields that define the card

        """
        list_fields = [
            'DMIL', self.mldtime_id, self.tstart, self.tend,
            self.dt, self.out_dt,
            self.print_flag, self.method]
        return list_fields

    def repr_fields(self):
        list_fields = self.raw_fields()
        return list_fields

    def write_card(self, size: int=8, is_double: bool=False) -> str:
        card = self.repr_fields()
        return self.comment + print_card_8(card)
