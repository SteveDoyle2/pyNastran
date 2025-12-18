from __future__ import annotations
from typing import TYPE_CHECKING

from pyNastran.bdf.field_writer_8 import print_card_8
from pyNastran.bdf.cards.base_card import BaseCard
from pyNastran.bdf.bdf_interface.assign_type import (
    integer, integer_or_blank, double, string,
    string_or_blank, string_multifield,
    blank, string_multifield_dollar_int,
    string_multifield_dollar_int_or_blank, double_or_blank,
    integer_or_string,
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
        # fields = [field for field in card.card if field is not None]
        # card.card = fields
        # card.nfields = len(card.card)

        # ['PLTMODE', '27',  'ASYM', '7', '0.4', 'tecplot', 'MODE07.p', 'lt']
        # ['PLTMODE', '104', 'ANTI', '4', None, '2.0', 'TECPLOT', 'GAFA_MOD', 'E4.PLT']
        set_id = integer(card, 1, 'set_id')

        symmetry = string(card, 2, 'sym/asym')
        mode = integer(card, 3, 'mode')
        mode_type = '' if card.field(4) is None else card.field(4)
        max_disp = double(card, 5, 'max_disp')
        # if max_disp is None:
        #     ifield += 1
        #     max_disp = double(card, ifield, 'max_disp')
        #     ifield += 1
        output_format = string(card, 6, 'format')

        filename = string_multifield_dollar_int(card, (7, 8), 'filename')
        # assert filename == 'GAFA_MODE4.PLT', filename
        aero_filename = string_multifield_dollar_int_or_blank(
            card, (9, 10), 'aero_filename', default='AEROGEOM.PAT')
        assert len(card) <= 9, f'len(PLTMODE card) = {len(card):d}\ncard={card}'
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
                 xaxis: str, filename: str | int, nmode: int=0,
                 output_format: str='TABLE', rho_ref: float=1.0,
                 comment: str=''):
        BaseCard.__init__(self)

        if comment:
            self.comment = comment

        self.set_id = set_id
        self.flutter_id = flutter_id
        self.nmode = nmode
        self.xaxis = xaxis
        self.output_format = output_format
        self.filename = filename
        self.rho_ref = rho_ref
        # assert out_format in {'TECPLOT', 'PATRAN', 'IDEAS', 'FEMAP', 'NASTRAN', 'NASTL', 'ANSYS'}, out_format
        # assert vct in {'YES', 'NO'}, vct
        # assert cell in {'YES', 'NO'}, cell
        assert xaxis in {'M', 'R', 'Q', 'H', 'V/VR', 'V', 'EQUV'}, f'field4={xaxis!r}'
        assert output_format in {'TABLE', 'IDEAS', 'FEMAP', 'ESA'}, f'output_format={output_format!r}'

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
        nmode = integer_or_blank(card, 3, 'field3', default=0)
        xaxis = string(card, 4, 'field4')
        output_format = string_or_blank(card, 5, 'form', default='TABLE')
        filename = string_multifield_dollar_int_or_blank(
            card, (6, 7), 'filename', default='')
        rho_ref = double_or_blank(card, 8, 'rho_ref', default=1.0)
        # assert filename == 'PLTVG_FRQ.DAT', filename
        assert len(card) <= 8, f'len(PLTVG card) = {len(card):d}\ncard={card}'
        return PLTVG(set_id, flutter_id, xaxis, filename,
                     nmode=nmode, output_format=output_format,
                     rho_ref=rho_ref, comment=comment)

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
        list_fields = ['PLTVG', self.set_id, self.flutter_id, self.nmode,
                       self.xaxis, self.output_format, filenamea, filenameb]
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

        # PLTCP SETID SYM IDMK IK MODE FORM FILENM CONT
        #       AERONM
        # ['PLTCP', '30', 'SYM', '80', '5', '1', 'TECPLOT', 'CP7.PLT']
        set_id = integer(card, 1, 'set_id')
        sym_flag = string(card, 2, 'field2')
        mkaero_id = integer(card, 3, 'field3')
        ik = integer(card, 4, 'ik')
        mode = integer_or_string(card, 5, 'mode')
        out_format = string_or_blank(card, 6, 'out_format', default='TECPLOT')
        filename = string_multifield_dollar_int_or_blank(
            card, (7, 8), 'filename', default='')
        aero_filename = string_multifield_dollar_int_or_blank(
            card, (9, 10), 'aero_filename', default='')
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

    def __init__(self, set_id: int, ase_id: int,
                 irow: int, icol: int, klist: int,
                 out_format: str, filename: str,
                 comment: str=''):
        BaseCard.__init__(self)

        if comment:
            self.comment = comment

        if out_format == 'TECP':
            out_format = 'TECPLOT'

        self.set_id = set_id
        self.ase_id = ase_id
        self.irow = irow
        self.icol = icol
        self.klist = klist
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

        # ['PLTMIST', '150', '150', '1', '1',  None, None,      'ROGER11', '.DAT']
        # ['PLTMIST', '31',  '30',  '1', '33', None, 'TECPLOT', 'QHG1.PL', 'T']
        set_id = integer(card, 1, 'set_id')
        ase_id = integer(card, 2, 'ase_id')
        irow = integer(card, 3, 'irow')
        icol = integer(card, 4, 'icol')
        klist = integer_or_blank(card, 5, 'klist', default=0)
        out_format = string_or_blank(card, 6, 'format', default='')

        filename = string_multifield_dollar_int(
            card, (7, 8), 'filename')
        # assert filename == 'ROGER11.DAT', f'filename={filename!r}'
        assert len(card) <= 9, f'len(PLTMIST card) = {len(card):d}\ncard={card}'
        return PLTMIST(set_id, ase_id, irow, icol, klist, out_format, filename,
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
        list_fields = ['PLTMIST', self.set_id, self.ase_id,
                       self.irow, self.icol, self.klist,
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
        assert out_format in {'TECPLOT', 'PATRAN', 'IDEAS', 'FEMAP', 'ANSYS', 'NASTRAN', 'NASTL'}, f'out_format={out_format!r}'

    @classmethod
    def add_card(cls, card: BDFCard, comment: str=''):
        # remove None
        #print(card)
        #fields = [field for field in card.card if field is not None]
        #card.card = fields
        card.nfields = len(card.card)

        # ['PLTSURF', '201', 'AILERON', '10.', 'TECPLOT', 'AILERON.', 'PLT']
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

        # PLTFLUT SETID IDFLUT MODE NTIME MAXDISP FORM FILENM
        #         ---AERONM---
        # ['PLTFLUT', '10', '10', '1', '25', '0.5', 'TECPLOT', 'OPENFLT', '.PLT']
        set_id = integer(card, 1, 'set_id')
        flutter_id = integer(card, 2, 'flutter_id')
        mode = integer(card, 3, 'mode')
        ntime = integer_or_blank(card, 4, 'ntime', default=1)
        scale_factor = double_or_blank(card, 5, 'scale_factor', default=1.0)
        out_format = string_or_blank(card, 6, 'out_format', default='TECPLOT')

        filename = string_multifield(card, (7, 8), 'filename')
        aero_filename = string_multifield_dollar_int_or_blank(
            card, (9, 10), 'aero_filename', default='AEROGEOM.PAT')
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


class PLTBODE(BaseCard):
    type = 'PLTBODE'
    def __init__(self, set_id: int,
                 cmargin_id: int,
                 fmin: float,
                 fmax: float,
                 nf: int,
                 filename: str,
                 log_scale_flag: int=0,
                 draw_flag: int=0,
                 comment: str=''):
        BaseCard.__init__(self)

        if comment:
            self.comment = comment
        self.set_id = set_id
        self.cmargin_id = cmargin_id
        self.fmin = fmin
        self.fmax = fmax
        self.nf = nf
        self.log_scale_flag = log_scale_flag
        self.draw_flag = draw_flag
        self.filename = filename

    @classmethod
    def add_card(cls, card: BDFCard, comment: str=''):
        # remove None
        #print(card)
        #fields = [field for field in card.card if field is not None]
        #card.card = fields
        # card.nfields = len(card.card)

        # PLTBODE SETID IDCMAR FMIN FMAX NF LOGSCAL FILENM
        #         DRAW
        # PLTBODE 20    10     50.  200. 21 0       BODE.PLT
        # ['PLTBODE', '101', '1',
        #  '5.', '70.',
        #  '51', '0', 'STATE_BO', 'DE.PLT']
        set_id = integer(card, 1, 'set_id')
        cmargin_id = integer(card, 2, 'field2')
        fmin = double(card, 3, 'fmin')
        fmax = double(card, 3, 'fmax')
        nf = integer(card, 5, 'field7')
        log_scale_flag = integer_or_blank(card, 6, 'field8', default=0)

        filename = string_multifield(card, (7, 8), 'filename')
        draw_flag = integer_or_blank(card, 9, 'draw_flag', default=0)
        # assert filename == 'ROGER11.DAT', f'filename={filename!r}'
        assert len(card) <= 10, f'len(PLTBODE card) = {len(card):d}\ncard={card}'
        return PLTBODE(set_id, cmargin_id, fmin, fmax, nf, filename,
                       log_scale_flag=log_scale_flag, draw_flag=draw_flag, comment=comment)

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
        list_fields = ['PLTBODE', self.set_id, self.cmargin_id,
                       self.fmin, self.fmax, self.nf,
                       self.log_scale_flag, filenamea, filenameb, self.draw_flag]
        return list_fields

    def write_card(self, size: int=8, is_double: bool=False) -> str:
        # TODO: needs a better writer
        card = self.repr_fields()
        return self.comment + print_card_8(card)


class PLTTIME(BaseCard):
    type = 'PLTTIME'
    def __init__(self, set_id: int,
                 mloads_id: int, tstart: float, tend: float, ndt: int,
                 out_type: str, filename: str, aero_filename: str,
                 output_format: str='TECPLOT', scale_factor: float=1.0,
                 comment: str=''):
        BaseCard.__init__(self)

        if comment:
            self.comment = comment
        if out_type == 'ELAS':
            out_type = 'ELASTIC'
        self.set_id = set_id
        self.mloads_id = mloads_id
        self.tstart = tstart
        self.tend = tend
        self.ndt = ndt
        self.out_type = out_type
        self.output_format = output_format
        self.scale_factor = scale_factor
        self.filename = filename
        self.aero_filename = aero_filename
        assert out_type in {'FORCE', 'FORCESOF', 'FORCERFA', 'MANEUVER', 'ELASTIC', 'NORIGID', 'UCP'}, f'out_type={out_type!r}'
        assert output_format in {'TECPLOT', 'PATRAN', 'IDEAS', 'FEMAP', 'OUTPUT4', 'NASTRAN', 'NASTNL'}, f'output_format={output_format!r}'

    @classmethod
    def add_card(cls, card: BDFCard, comment: str=''):
        # remove None
        #print(card)
        #fields = [field for field in card.card if field is not None]
        #card.card = fields
        # card.nfields = len(card.card)

        # PLTTIME IDPLT  IDMLD TS    TE     NDT TYPE    FORM    SCALE
        #         ---FILENM--- ---AERONM---
        # PLTTIME 10     20    -2.0  1.0    10  ELASTIC TECPLOT 1.0
        #         TECPLOT.PLT
        set_id = integer(card, 1, 'set_id')
        mloads_id = integer(card, 2, 'mloads_id')
        tstart = double(card, 3, 'tstart')
        tend = double(card, 4, 'tend')
        ndt = integer(card, 5, 'ndt')
        out_type = string(card, 6, 'out_type')
        output_format = string_or_blank(card, 7, 'output_format', default='TECPLOT')
        scale_factor = double_or_blank(card, 8, 'draw_flag', default=1.0)

        filename = string_multifield_dollar_int(card, (9, 10), 'filename')
        aero_filename = string_multifield_dollar_int(card, (10, 11), 'aero_filename')
        # assert filename == 'ROGER11.DAT', f'filename={filename!r}'
        assert len(card) <= 11, f'len(PLTTIME card) = {len(card):d}\ncard={card}'
        return PLTTIME(set_id, mloads_id, tstart, tend, ndt, out_type, filename,
                       aero_filename, output_format=output_format, scale_factor=scale_factor, comment=comment)

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
        aerofilenamea, aerofilenameb = split_filename_dollar(self.aero_filename)
        list_fields = ['PLTTIME', self.set_id, self.mloads_id,
                       self.tstart, self.tend, self.ndt,
                       self.out_type, self.output_format, self.scale_factor,
                       filenamea, filenameb, aerofilenamea, aerofilenameb]
        return list_fields

    def write_card(self, size: int=8, is_double: bool=False) -> str:
        # TODO: needs a better writer
        card = self.repr_fields()
        return self.comment + print_card_8(card)


class PLTTRIM(BaseCard):
    type = 'PLTTRIM'
    def __init__(self, set_id: int,
                 trim_id: int,
                 out_type: str, filename: str, aero_filename: str,
                 flex: str='FLEX', output_format: str='TECPLOT',
                 scale_factor: float=1.0,
                 comment: str=''):
        BaseCard.__init__(self)

        if comment:
            self.comment = comment
        if out_type == 'ELAS':
            out_type = 'ELASTIC'
        self.set_id = set_id
        self.trim_id = trim_id
        self.flex = flex
        self.out_type = out_type
        self.output_format = output_format
        self.scale_factor = scale_factor
        self.filename = filename
        self.aero_filename = aero_filename
        assert flex in {'RIGID', 'FLEX'}, f'flex={flex!r}'
        assert out_type in {'FORCE', 'AERO', 'INERTIAL', 'CP', 'DEFORM', 'ELASTIC'}, f'out_type={out_type!r}'
        assert output_format in {'TECPLOT', 'PATRAN', 'IDEAS', 'FEMAP', 'OUTPUT4', 'NASTRAN', 'NASTNL'}, f'output_format={output_format!r}'

    @classmethod
    def add_card(cls, card: BDFCard, comment: str=''):
        # remove None
        #print(card)
        #fields = [field for field in card.card if field is not None]
        #card.card = fields
        # card.nfields = len(card.card)

        # PLTTRIM IDPLT IDTRIM FLEX TYPE   FORM    ---FILENM--- SCALE
        #         ---AERONM---
        # PLTTRIM 100   10     FLEX DEFORM TECPLOT PLTTRIM.DAT
        set_id = integer(card, 1, 'set_id')
        trim_id = integer(card, 2, 'trim_id')
        flex = string_or_blank(card, 3, 'flex', default='FLEX')
        out_type = string(card, 4, 'out_type')
        output_format = string_or_blank(card, 5, 'output_format', default='TECPLOT')
        filename = string_multifield_dollar_int(card, (6, 7), 'filename')
        scale_factor = double_or_blank(card, 8, 'draw_flag', default=1.0)
        aero_filename = string_multifield_dollar_int_or_blank(
            card, (9, 10), 'aero_filename')
        # assert filename == 'ROGER11.DAT', f'filename={filename!r}'
        assert len(card) <= 10, f'len(PLTTRIM card) = {len(card):d}\ncard={card}'
        return PLTTRIM(set_id, trim_id, out_type, filename,
                       aero_filename, flex=flex,
                       output_format=output_format, scale_factor=scale_factor,
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
        aerofilenamea, aerofilenameb = split_filename_dollar(self.aero_filename)
        list_fields = ['PLTTRIM', self.set_id, self.trim_id,
                       self.flex, self.out_type, self.output_format,
                       filenamea, filenameb, self.scale_factor,
                       aerofilenamea, aerofilenameb]
        return list_fields

    def write_card(self, size: int=8, is_double: bool=False) -> str:
        # TODO: needs a better writer
        card = self.repr_fields()
        return self.comment + print_card_8(card)
