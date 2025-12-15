from __future__ import annotations
from typing import TYPE_CHECKING

from pyNastran.bdf.field_writer_8 import print_card_8
from pyNastran.bdf.cards.base_card import BaseCard
from pyNastran.bdf.bdf_interface.assign_type import (
    integer, integer_or_blank, double, string,
    string_or_blank, string_multifield,
    blank, string_multifield_dollar_int_or_blank,
)
from .utils import split_filename_dollar

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

        if output_format == 'TECP':
            output_format = 'TECPLOT'
        self.set_id = set_id
        self.mode = mode
        self.symmetry = symmetry
        self.max_disp = max_disp
        self.output_format = output_format
        self.filename = filename
        assert output_format in {'TECPLOT', 'FEMAP', 'NASTRAN', 'PATRAN'}, output_format
        assert symmetry in {'SYM', 'ASYM', 'ANTI'}, symmetry

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
        # if max_disp is None:
        #     ifield += 1
        #     max_disp = double(card, ifield, 'max_disp')
        #     ifield += 1
        output_format = string(card, 5, 'format')

        filename = string_multifield(card, (6, 7), 'filename')
        assert len(card) <= 8, f'len(PLTMODE card) = {len(card):d}\ncard={card}'
        return PLTMODE(set_id, symmetry, mode, max_disp,
                       output_format, filename, comment=comment)

    def cross_reference(self, model: BDF) -> None:
        return

    def safe_cross_reference(self, model: BDF, xref_errors) -> None:
        self.cross_reference(model)

    def uncross_reference(self) -> None:
        pass

    def repr_fields(self):
        """
        Gets the fields in their simplified form

        Returns
        -------
        fields : list[varies]
          the fields that define the card

        """
        filenamea, filenameb = split_filename_dollar(self.filename)

        list_fields = [
            'PLTMODE', self.set_id, self.symmetry,
            self.mode, self.max_disp,
            self.output_format, filenamea, filenameb]
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

    def safe_cross_reference(self, model: BDF, xref_errors) -> None:
        self.cross_reference(model)

    def uncross_reference(self) -> None:
        pass

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

class PLTVG(BaseCard):
    type = 'PLTVG'

    def __init__(self, set_id: int, flutter_id: int,
                 filename: str | int,
                 comment: str=''):
        BaseCard.__init__(self)

        if comment:
            self.comment = comment

        self.set_id = set_id
        self.flutter_id = flutter_id
        self.filename = filename
        # assert out_format in {'TECPLOT', 'PATRAN', 'IDEAS', 'FEMAP', 'NASTRAN', 'NASTL', 'ANSYS'}, out_format
        # assert vct in {'YES', 'NO'}, vct
        # assert cell in {'YES', 'NO'}, cell

    @classmethod
    def add_card(cls, card: BDFCard, comment: str=''):
        # remove None
        # print(card)
        # fields = [field for field in card.card if field is not None]
        # card.card = fields
        card.nfields = len(card.card)

        # ['PLTVG', '10',  '1',   '4',  'Q', 'TABLE', 'WEKZTRAN', '_678.DAT']
        # ['PLTVG', '11',  '100', None, 'M', None,    'PLTVG.10', '0']
        # ['PLTVG', '11',  '100', None, 'V', None,    'VG7.PLT']
        # ['PLTVG', '100', '100', None, 'Q', None,    'PLTVG_FR', 'Q.DAT']
        set_id = integer(card, 1, 'set_id')
        flutter_id = integer(card, 2, 'flutter_id')
        field3 = integer_or_blank(card, 3, 'field3', default=0)
        field4 = string(card, 4, 'field4')
        assert field4 in {'M', 'Q', 'V'}, f'field4={field4!r}'
        field5 = string_or_blank(card, 5, 'field5', default='')
        # assert field4 == 0, field4
        filename = string_multifield_dollar_int_or_blank(
            card, (6, 7), 'filename', default='')
        # assert filename == 'PLTVG_FRQ.DAT', filename
        assert len(card) <= 8, f'len(PLTVG card) = {len(card):d}\ncard={card}'
        return PLTVG(set_id, flutter_id, filename,
                     comment=comment)

    def cross_reference(self, model: BDF) -> None:
        return

    def safe_cross_reference(self, model: BDF, xref_errors) -> None:
        self.cross_reference(model)

    def uncross_reference(self) -> None:
        pass

    def repr_fields(self):
        """
        Gets the fields in their simplified form

        Returns
        -------
        fields : list[varies]
          the fields that define the card

        """
        filenamea, filenameb = split_filename_dollar(self.filename)
        list_fields = ['PLTVG', self.set_id, self.flutter_id, None,
                       None, filenamea, filenameb]
        return list_fields

    def write_card(self, size: int=8, is_double: bool=False) -> str:
        # TODO: needs a better writer
        card = self.repr_fields()
        return self.comment + print_card_8(card)


class PLTCP(BaseCard):
    type = 'PLTCP'
    def __init__(self, set_id: int,
                 sym_flag: str,
                 out_format: str,
                 filename: str | int,
                 comment: str=''):
        BaseCard.__init__(self)

        if comment:
            self.comment = comment

        self.set_id = set_id
        self.sym_flag = sym_flag
        self.out_format = out_format
        self.filename = filename
        assert sym_flag in {'SYM', 'ANTI'}, sym_flag
        assert out_format in {'TECPLOT',}, out_format
        # assert vct in {'YES', 'NO'}, vct
        # assert cell in {'YES', 'NO'}, cell

    @classmethod
    def add_card(cls, card: BDFCard, comment: str=''):
        # remove None
        # print(card)
        # fields = [field for field in card.card if field is not None]
        # card.card = fields
        card.nfields = len(card.card)

        ['PLTCP', '30', 'SYM', '80', '5', '1', 'TECPLOT', 'CP7.PLT']
        set_id = integer(card, 1, 'set_id')
        sym_flag = string(card, 2, 'field2')
        field3 = integer(card, 3, 'field3')
        field4 = integer(card, 4, 'field4')
        field5 = integer(card, 5, 'field5')
        out_format = string(card, 6, 'field5')
        filename = string_multifield_dollar_int_or_blank(
            card, (7, 8), 'filename', default='')
        # assert filename == 'CP7.PLT', filename
        assert len(card) <= 9, f'len(PLTCP card) = {len(card):d}\ncard={card}'
        return PLTCP(set_id, sym_flag, out_format, filename,
                     comment=comment)

    def cross_reference(self, model: BDF) -> None:
        return

    def safe_cross_reference(self, model: BDF, xref_errors) -> None:
        self.cross_reference(model)

    def uncross_reference(self) -> None:
        pass

    def repr_fields(self):
        """
        Gets the fields in their simplified form

        Returns
        -------
        fields : list[varies]
          the fields that define the card

        """
        # ['PLTCP', '30', 'SYM', '80', '5', '1', 'TECPLOT', 'CP7.PLT']
        list_fields = ['PLTCP', self.set_id, self.sym_flag, None,
                       None, None, self.out_format, self.filename]
        return list_fields

    def write_card(self, size: int = 8, is_double: bool = False) -> str:
        # TODO: needs a better writer
        card = self.repr_fields()
        return self.comment + print_card_8(card)


class PLTMIST(BaseCard):
    type = 'PLTMIST'

    def __init__(self, set_id: int,
                 out_format: str, filename: str,
                 comment: str=''):
        BaseCard.__init__(self)

        if comment:
            self.comment = comment

        if out_format == 'TECP':
            out_format = 'TECPLOT'

        self.set_id = set_id
        self.out_format = out_format
        self.filename = filename
        assert out_format in {'TECPLOT', ''}, f'out_format={out_format!r}'

    @classmethod
    def add_card(cls, card: BDFCard, comment: str=''):
        # remove None
        #print(card)
        #fields = [field for field in card.card if field is not None]
        #card.card = fields
        card.nfields = len(card.card)

        ['PLTMIST', '150', '150', '1', '1',  None, None,      'ROGER11', '.DAT']
        ['PLTMIST', '31',  '30',  '1', '33', None, 'TECPLOT', 'QHG1.PL', 'T']
        set_id = integer(card, 1, 'set_id')
        field2 = integer(card, 2, 'field2')
        field3 = integer(card, 3, 'field3')
        field4 = integer(card, 4, 'field4')
        field5 = blank(card, 5, 'field5')
        # assert femgrid in {'YES', 'NO', ''}, f'femgrid={femgrid!r}'

        out_format = string_or_blank(card, 6, 'format', default='')

        filename = string_multifield(card, (7, 8), 'filename')
        # assert filename == 'ROGER11.DAT', f'filename={filename!r}'
        assert len(card) <= 9, f'len(PLTMIST card) = {len(card):d}\ncard={card}'
        return PLTMIST(set_id, out_format, filename,
                       comment=comment)

    def cross_reference(self, model: BDF) -> None:
        return

    def safe_cross_reference(self, model: BDF, xref_errors) -> None:
        self.cross_reference(model)

    def uncross_reference(self) -> None:
        pass

    def repr_fields(self):
        """
        Gets the fields in their simplified form

        Returns
        -------
        fields : list[varies]
          the fields that define the card

        """
        filenamea, filenameb = split_filename_dollar(self.filename)
        list_fields = ['PLTMIST', self.set_id, None, None, None, None,
                       self.out_format, filenamea, filenameb]
        return list_fields

    def write_card(self, size: int=8, is_double: bool=False) -> str:
        # TODO: needs a better writer
        card = self.repr_fields()
        return self.comment + print_card_8(card)


class PLTSURF(BaseCard):
    type = 'PLTSURF'
    def __init__(self, set_id: int, label: str,
                 out_format: str, filename: str,
                 scale_factor: float=1.0,
                 comment: str=''):
        BaseCard.__init__(self)

        if comment:
            self.comment = comment

        if out_format == 'TECP':
            out_format = 'TECPLOT'

        self.set_id = set_id
        self.label = label
        self.out_format = out_format
        self.scale_factor = scale_factor
        self.filename = filename
        assert out_format in {'TECPLOT', ''}, f'out_format={out_format!r}'

    @classmethod
    def add_card(cls, card: BDFCard, comment: str=''):
        # remove None
        #print(card)
        #fields = [field for field in card.card if field is not None]
        #card.card = fields
        card.nfields = len(card.card)

        ['PLTSURF', '201', 'AILERON', '10.', 'TECPLOT', 'AILERON.', 'PLT']
        set_id = integer(card, 1, 'set_id')
        label = string(card, 2, 'label')
        scale_factor = double(card, 3, 'scale_factor')
        out_format = string(card, 4, 'out_format')

        filename = string_multifield(card, (5, 6), 'filename')
        # assert filename == 'ROGER11.DAT', f'filename={filename!r}'
        assert len(card) <= 7, f'len(PLTSURF card) = {len(card):d}\ncard={card}'
        return PLTSURF(set_id, label, out_format, filename,
                       scale_factor=scale_factor, comment=comment)

    def cross_reference(self, model: BDF) -> None:
        return

    def safe_cross_reference(self, model: BDF, xref_errors) -> None:
        self.cross_reference(model)

    def uncross_reference(self) -> None:
        pass

    def repr_fields(self):
        """
        Gets the fields in their simplified form

        Returns
        -------
        fields : list[varies]
          the fields that define the card

        """
        filenamea, filenameb = split_filename_dollar(self.filename)
        list_fields = [
            'PLTSURF', self.set_id, self.label, self.scale_factor,
            self.out_format, filenamea, filenameb]
        return list_fields

    def write_card(self, size: int=8, is_double: bool=False) -> str:
        # TODO: needs a better writer
        card = self.repr_fields()
        return self.comment + print_card_8(card)


class PLTFLUT(BaseCard):
    type = 'PLTFLUT'
    def __init__(self, set_id: int,
                 out_format: str, filename: str,
                 scale_factor: float=1.0,
                 comment: str=''):
        BaseCard.__init__(self)

        if comment:
            self.comment = comment

        if out_format == 'TECP':
            out_format = 'TECPLOT'

        self.set_id = set_id
        self.out_format = out_format
        self.scale_factor = scale_factor
        self.filename = filename
        assert out_format in {'TECPLOT', ''}, f'out_format={out_format!r}'

    @classmethod
    def add_card(cls, card: BDFCard, comment: str=''):
        # remove None
        #print(card)
        #fields = [field for field in card.card if field is not None]
        #card.card = fields
        card.nfields = len(card.card)

        # ['PLTFLUT', '10', '10', '1', '25', '0.5', 'TECPLOT', 'OPENFLT', '.PLT']
        set_id = integer(card, 1, 'set_id')
        field2 = integer(card, 2, 'field2')
        field3 = integer(card, 3, 'field3')
        field4 = integer(card, 4, 'field4')
        scale_factor = double(card, 5, 'scale_factor')
        out_format = string(card, 6, 'out_format')

        filename = string_multifield(card, (7, 8), 'filename')
        # assert filename == 'ROGER11.DAT', f'filename={filename!r}'
        assert len(card) == 9, f'len(PLTFLUT card) = {len(card):d}\ncard={card}'
        return PLTFLUT(set_id, out_format, filename,
                       scale_factor=scale_factor, comment=comment)

    def cross_reference(self, model: BDF) -> None:
        return

    def safe_cross_reference(self, model: BDF, xref_errors) -> None:
        self.cross_reference(model)

    def uncross_reference(self) -> None:
        pass

    def repr_fields(self):
        """
        Gets the fields in their simplified form

        Returns
        -------
        fields : list[varies]
          the fields that define the card

        """
        filenamea, filenameb = split_filename_dollar(self.filename)
        list_fields = ['PLTFLUT', self.set_id, None, None, None,
                       self.scale_factor, self.out_format, filenamea, filenameb]
        return list_fields

    def write_card(self, size: int=8, is_double: bool=False) -> str:
        # TODO: needs a better writer
        card = self.repr_fields()
        return self.comment + print_card_8(card)
