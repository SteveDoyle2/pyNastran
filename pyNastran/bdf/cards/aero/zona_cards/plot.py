from __future__ import annotations
from typing import TYPE_CHECKING

from pyNastran.bdf.field_writer_8 import print_card_8
from pyNastran.bdf.cards.base_card import BaseCard
from pyNastran.bdf.bdf_interface.assign_type import (
    integer, integer_or_blank, double, string,
    string_or_blank, string_multifield,
)

if TYPE_CHECKING:  # pragma: no cover
    from pyNastran.bdf.bdf import BDF
    from pyNastran.bdf.bdf_interface.bdf_card import BDFCard


class PLTMODE(BaseCard):
    """
    """
    type = 'PLTMODE'

    def __init__(self, set_id: int, symmetry: str,
                 mode: int, max_disp: float,
                 output_format: str, filename: str,
                 comment: str=''):
        BaseCard.__init__(self)

        if comment:
            self.comment = comment
        self.set_id = set_id
        self.mode = mode
        self.symmetry = symmetry
        self.max_disp = max_disp
        self.output_format = output_format
        self.filename = filename

    @classmethod
    def add_card(cls, card: BDFCard, comment: str=''):
        # remove None
        fields = [field for field in card.card if field is not None]
        card.card = fields
        card.nfields = len(card.card)

        #['PLTMODE', '27', 'ASYM', '7', '0.4', 'tecplot', 'MODE07.p', 'lt']
        set_id = integer(card, 1, 'set_id')

        symmetry = string(card, 2, 'sym/asym')
        mode = integer(card, 3, 'mode')
        max_disp = double(card, 4, 'max_disp')
        assert symmetry in {'SYM', 'ASYM', 'ANTI'}, symmetry
        # if max_disp is None:
        #     ifield += 1
        #     max_disp = double(card, ifield, 'max_disp')
        #     ifield += 1
        output_format = string(card, 5, 'format')
        if output_format == 'TECP':
            output_format = 'TECPLOT'
        assert output_format in {'TECPLOT', 'FEMAP', 'NASTRAN', 'PATRAN'}, format

        filename = string_multifield(card, (6, 7), 'filename')
        assert len(card) <= 8, f'len(PLTMODE card) = {len(card):d}\ncard={card}'
        return PLTMODE(set_id, symmetry, mode, max_disp,
                       output_format, filename, comment=comment)

    def cross_reference(self, model: BDF) -> None:
        return

    def repr_fields(self):
        """
        Gets the fields in their simplified form

        Returns
        -------
        fields : list[varies]
          the fields that define the card

        """
        list_fields = ['PLTMODE', self.set_id, self.mode, self.max_disp,
                       self.output_format, self.filename]
        return list_fields

    def write_card(self, size: int=8, is_double: bool=False) -> str:
        card = self.repr_fields()
        return self.comment + print_card_8(card)


class PLTAERO(BaseCard):
    """
    PLTAERO SETID FEMGRID OFFSET FORM    FILENM   CELL VCT
    PLTAERO 100   YES     100000 TECPLOT AERO.PLT YES  YES
    """
    type = 'PLTAERO'

    def __init__(self, set_id: int, femgrid: str,
                 offset: int, out_format: str, filename: str,
                 cell: str='NO', vct: str='NO', comment: str=''):
        BaseCard.__init__(self)

        if comment:
            self.comment = comment

        if out_format == 'TECP':
            out_format = 'TECPLOT'

        self.set_id = set_id
        self.femgrid = femgrid
        self.offset = offset
        self.out_format = out_format
        self.filename = filename
        self.cell = cell
        self.vct = vct
        assert out_format in {'TECPLOT', 'PATRAN', 'IDEAS', 'FEMAP', 'NASTRAN', 'NASTL', 'ANSYS'}, out_format
        assert vct in {'YES', 'NO'}, vct
        assert cell in {'YES', 'NO'}, cell

    @classmethod
    def add_card(cls, card: BDFCard, comment: str=''):
        # remove None
        #print(card)
        #fields = [field for field in card.card if field is not None]
        #card.card = fields
        card.nfields = len(card.card)

        #['PLTMODE', '27', 'ASYM', '7', '0.4', 'tecplot', 'MODE07.p', 'lt']
        set_id = integer(card, 1, 'set_id')
        femgrid = string_or_blank(card, 2, 'id_flutter', default='')
        assert femgrid in {'YES', 'NO', ''}, f'femgrid={femgrid!r}'

        offset = integer_or_blank(card, 3, 'mode', default=-1)
        #ntime = integer(card, 4, 'ntime')
        #max_disp = double(card, 4, 'max_disp')

        out_format = string(card, 4, 'format')

        filename = string_multifield(card, (5, 6), 'filename')

        cell = string_or_blank(card, 7, 'cell', default='NO')
        vct = string_or_blank(card, 8, 'vct', default='NO')
        assert len(card) <= 9, f'len(PLTAERO card) = {len(card):d}\ncard={card}'
        return PLTAERO(set_id, femgrid, offset, out_format, filename,
                       cell=cell, vct=vct, comment=comment)

    def cross_reference(self, model: BDF) -> None:
        return

    def repr_fields(self):
        """
        Gets the fields in their simplified form

        Returns
        -------
        fields : list[varies]
          the fields that define the card

        """
        list_fields = ['PLTAERO', self.set_id, self.femgrid, self.offset,
                       self.out_format, self.filename, self.cell, self.vct]
        return list_fields

    def write_card(self, size: int=8, is_double: bool=False) -> str:
        # TODO: needs a better writer
        card = self.repr_fields()
        return self.comment + print_card_8(card)