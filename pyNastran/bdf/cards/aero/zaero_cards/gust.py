# coding: utf-8
# pylint: disable=W0212,C0103
"""
All ZONA aero cards are defined in this file.  This includes:
 * TRIM

All cards are BaseCard objects.

"""
from __future__ import annotations
from itertools import count
from typing import TYPE_CHECKING

from matplotlib import pyplot as plt

from pyNastran.utils.numpy_utils import integer_types
from pyNastran.bdf.field_writer_8 import print_card_8
from pyNastran.bdf.cards.base_card import BaseCard
from pyNastran.bdf.bdf_interface.assign_type import (
    integer, integer_or_blank, double_or_blank, string,
    string_or_blank, double, integer_or_string,
    integer_or_double,
    # integer_or_string, integer_string_or_blank,
    string_multifield, # parse_components as fcomponent
)
from pyNastran.bdf.cards.aero.zaero_cards.ase import (
    ASECONT, CJUNCT, )
from pyNastran.bdf.cards.aero.zaero_cards.cards import MLDCOMD

if TYPE_CHECKING:  # pragma: no cover
    from pyNastran.bdf.bdf import BDF
    from pyNastran.bdf.bdf_interface.bdf_card import BDFCard


class GLOADS(BaseCard):
    type = 'GLOADS'

    # _field_map = {
    #     1: 'sid', 2: 'mach', 3: 'q', 8: 'aeqr',
    # }

    def __init__(self, gloads_id: int, asecont_id: int,
                 flutter_id: int, minstat_id: int,
                 mldstat_id: int, mldcomd_id: int, mldtime_id: int,
                 mldprnt_id: int,
                 save_flag: str, form: str, filename: str,
                 save_freq: str, filename_freq: str, comment: str=''):
        BaseCard.__init__(self)
        if comment:
            self.comment = comment

        save_flag = save_flag.upper()
        if save_flag == 'ACQU':
            save_flag = 'ACQUIRE'

        self.gloads_id = gloads_id
        self.asecont_id = asecont_id
        self.flutter_id = flutter_id
        self.minstat_id = minstat_id

        self.mldstat_id = mldstat_id
        self.mldcomd_id = mldcomd_id
        self.mldtime_id = mldtime_id
        self.mldprnt_id = mldprnt_id

        self.form = form
        self.save_flag = save_flag
        self.save_freq = save_freq
        self.filename_freq = filename_freq
        self.filename = filename

        self.asecont_ref = None
        self.flutter_ref = None
        self.minstat_ref = None
        self.mldstat_ref = None
        self.mldcomd_ref = None
        self.mldtime_ref = None
        self.mldprnt_ref = None
        assert form in {'FORMAT', 'FORMAT23', 'UNFORM', ''}, f'form={form!r}'
        assert save_flag in {'SAVE', 'ACQUIRE', ''}, f'save_flag={save_flag!r}'
        assert save_freq in {'SAVE', 'ACQUIRE', ''}, f'save_freq={save_freq!r}'

    @classmethod
    def add_card(cls, card: BDFCard, comment: str = ''):
        """
        Adds a GLOADS card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card

        """
        filename = ''
        filename_freq = ''

        # GLOADS SID  CONID FLTID  RAAID   STATES  IDGUST TIME MLDPRNT CONT
        #        SAVE FORM  ----FILENM---- SAVEFRQ FILEFRQ
        gloads_id = integer(card, 1, 'gloads_id')
        asecont_id = integer_or_blank(card, 2, 'asecont_id', default=0)
        flutter_id = integer_or_blank(card, 3, 'flutter_id', default=0)
        minstat_id = integer_or_blank(card, 4, 'raaid/minstat_id', default=0)

        mldstat_id = integer_or_blank(card, 5, 'states/mldstat_id', default=0)
        mldcomd_id = integer(card, 6, 'command/dgust/mldcomd_id')
        mldtime_id = integer_or_blank(card, 7, 'time/mldtime_id', default=0)
        mldprnt_id = integer_or_blank(card, 8, 'mldprnt_id', default=0)

        save_flag = string_or_blank(card, 9, 'save_flag', default='')
        form = string_or_blank(card, 10, 'save_flag', default='')

        if len(card) > 12:
            filename = string_multifield(card, (12, 13), 'filename')
        save_freq = string_or_blank(card, 14, 'save_freq', default='')

        if len(card) > 14:
            filename_freq = string_multifield(card, (15, 16), 'filename')
        assert 7 <= len(card) < 14, f'len(GLOADS card) = {len(card):d}\ncard={card}'
        return GLOADS(gloads_id, asecont_id, flutter_id, minstat_id,
                      mldstat_id, mldcomd_id, mldtime_id,
                      mldprnt_id, save_flag, form, filename,
                      save_freq, filename_freq, comment=comment)

    # def validate(self):
    #     assert self.true_g in ['TRUE', 'G'], 'true_g=%r' % self.true_g

    def cross_reference(self, model: BDF) -> None:
        zaero = model.zaero
        if self.asecont_id > 0:
            self.asecont_ref = zaero.asecont[self.asecont_id]
        if self.flutter_id > 0:
            self.flutter_ref = model.flutters[self.flutter_id]
        if self.minstat_id:
            self.minstat_ref = zaero.minstat[abs(self.minstat_id)]

        if self.mldstat_id > 0:
            self.mldstat_ref = zaero.mldstat[self.mldstat_id]

        if self.mldcomd_id in zaero.dgust:
            self.mldcomd_ref = zaero.dgust[self.mldcomd_id]
        elif self.mldcomd_id in zaero.cgust:
            self.mldcomd_ref = zaero.cgust[self.mldcomd_id]
        else:
            dgust = list(zona.dgust)
            cgust = list(zona.cgust)
            msg = (
                f'{self.mldcomd_id} is not in [DGUST, CGUST]\n'
                f' - dgust = {dgust}\n'
                f' - cgust = {cgust}')
            raise RuntimeError(msg)

        if self.mldcomd_id in zaero.mldcomd:
            self.mldcomd_ref = zaero.mldcomd[self.mldcomd_id]
        if self.mldtime_id:
            self.mldtime_ref = zaero.mldtime[self.mldtime_id]
        self.mldprnt_ref = zaero.mldprnt[self.mldprnt_id]

    def safe_cross_reference(self, model: BDF, xref_errors):
        self.cross_reference(model)

    def uncross_reference(self) -> None:
        """Removes cross-reference links"""
        self.asecont_ref = None
        self.flutter_ref = None
        self.minstat_ref = None
        self.mldstat_ref = None
        self.mldcomd_ref = None
        self.mldtime_ref = None
        self.mldprnt_ref = None

    def plot(self, fig: plt.Figure):
        # if self.mldcomd_ref is None:
        #     return
        mldcomd = self.mldcomd_ref
        ntables = len(mldcomd.extinps_ref)
        assert ntables > 0, ntables
        axes = fig.subplots(nrows=ntables)
        if ntables == 1:
            axes = [axes]
        assert len(axes) == ntables, (axes, ntables)

        for i, extinp, table in zip(count(), mldcomd.extinps_ref, mldcomd.tables_ref):
            label = extinp.label
            ax = axes[i]
            ax.grid(True)
            ax.set_ylabel(label)
            ax.set_xlabel('Time (sec)')
            ax.plot(table.x, table.y, label=label)
            ax.legend()

    def raw_fields(self):
        """
        Gets the fields in their unmodified form

        Returns
        -------
        fields : list[varies]
            the fields that define the card

        """
        list_fields = [
            'GLOADS', self.gloads_id, self.asecont_id,
            self.flutter_id, self.minstat_id,
            self.mldstat_id, self.mldcomd_id, self.mldtime_id,
            self.mldprnt_id, self.save_flag, self.form,
            self.filename, self.save_freq, self.filename_freq]
        return list_fields

    def repr_fields(self):
        list_fields = self.raw_fields()
        return list_fields

    def write_card(self, size: int = 8, is_double: bool = False) -> str:
        card = self.repr_fields()
        return self.comment + print_card_8(card)


class DGUST(BaseCard):
    type = 'DGUST'

    # _field_map = {
    #     1: 'sid', 2: 'mach', 3: 'q', 8: 'aeqr',
    # }

    def __init__(self, dgust_id: int, gust_type: str,
                 length_gust: str | float, gust_velocity: float,
                 x0: float, fmax: float=0.0, df: float=0.01, nap: int=0,
                 comment: str=''):
        BaseCard.__init__(self)
        if comment:
            self.comment = comment

        self.dgust_id = dgust_id
        self.gust_type = gust_type
        self.length_gust = length_gust
        self.gust_velocity = gust_velocity

        self.x0 = x0
        self.fmax = fmax
        self.df = df
        self.nap = nap

    @classmethod
    def add_card(cls, card: BDFCard, comment: str = ''):
        """
        Adds a DGUST card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card

        """
        # DGUST SID TYPE  LG  WGV X0   FMAX DF    NAP
        # DGUST 10  OMCOS 1.0 0.1 -3.0 90.0 0.001
        dgust_id = integer(card, 1, 'dgust_id')
        gust_type = string(card, 2, 'gust_type')
        assert gust_type in {'STEP', 'SINE', 'OMCOS', 'RANDOM', 'TABLE', 'STEPP', 'SINEP', 'OMCOSP', 'RANDOMP', 'TABLEP'}, f'gust_type={gust_type!r}'
        length_gust = integer_or_double(card, 3, 'Lg, length_gust')
        gust_velocity = double(card, 4, 'wgv, gust_velocity')
        x0 = double(card, 5, 'x0')
        fmax = double_or_blank(card, 6, 'fmax', default=0)
        df = double_or_blank(card, 7, 'df', default=0.01)
        nap = integer_or_blank(card, 8, 'df', default=0)

        assert len(card) < 9, f'len(DGUST card) = {len(card):d}\ncard={card}'
        return DGUST(dgust_id, gust_type, length_gust, gust_velocity,
                     x0, fmax=fmax, df=df, nap=nap, comment=comment)

    # def validate(self):
    #     assert self.true_g in ['TRUE', 'G'], 'true_g=%r' % self.true_g

    def cross_reference(self, model: BDF) -> None:
        pass

    def safe_cross_reference(self, model: BDF, xref_errors):
        self.cross_reference(model)

    def uncross_reference(self) -> None:
        """Removes cross-reference links"""
        pass

    def plot(self, fig: plt.Figure):
        # if self.mldcomd_ref is None:
        #     return
        mldcomd = self.mldcomd_ref
        ntables = len(mldcomd.extinps_ref)
        assert ntables > 0, ntables
        axes = fig.subplots(nrows=ntables)
        if ntables == 1:
            axes = [axes]
        assert len(axes) == ntables, (axes, ntables)

        for i, extinp, table in zip(count(), mldcomd.extinps_ref, mldcomd.tables_ref):
            label = extinp.label
            ax = axes[i]
            ax.grid(True)
            ax.set_ylabel(label)
            ax.set_xlabel('Time (sec)')
            ax.plot(table.x, table.y, label=label)
            ax.legend()

    def raw_fields(self):
        """
        Gets the fields in their unmodified form

        Returns
        -------
        fields : list[varies]
            the fields that define the card

        """
        list_fields = [
            'DGUST', self.dgust_id, self.gust_type,
            self.length_gust, self.gust_velocity,
            self.x0, self.fmax, self.df, self.nap]
        return list_fields

    def repr_fields(self):
        list_fields = self.raw_fields()
        return list_fields

    def write_card(self, size: int = 8, is_double: bool = False) -> str:
        card = self.repr_fields()
        return self.comment + print_card_8(card)


class CGUST(BaseCard):
    type = 'CGUST'

    # _field_map = {
    #     1: 'sid', 2: 'mach', 3: 'q', 8: 'aeqr',
    # }

    def __init__(self, cgust_id: int, gust_type: str,
                 one_over_velocity: float,
                 x0: float, rms_velocity: float,
                 length_gust: float=2500.,
                 fmax: float=0.0, df: float=0.01, comment: str=''):
        BaseCard.__init__(self)
        if comment:
            self.comment = comment

        self.cgust_id = cgust_id
        self.gust_type = gust_type
        self.length_gust = length_gust
        self.one_over_velocity = one_over_velocity

        self.x0 = x0
        self.fmax = fmax
        self.df = df
        self.rms_velocity = rms_velocity

    @classmethod
    def add_card(cls, card: BDFCard, comment: str = ''):
        """
        Adds a DGUST card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card

        """
        # CGUST SID TYPE LG  VINV X0  FMAX DF    GURMS
        # CGUST 10  VK2  1.0 0.1 -3.0 90.0 0.001 100.0
        cgust_id = integer(card, 1, 'cgust_id')
        gust_type = string(card, 2, 'gust_type')
        length_gust = double_or_blank(card, 3, 'Lg, length_gust', default=2500.0)
        one_over_velocity = double(card, 4, 'vinv, 1/forward_velocity')
        x0 = double_or_blank(card, 5, 'x0')
        fmax = double_or_blank(card, 6, 'fmax', default=0)
        df = double_or_blank(card, 7, 'df', default=0.01)
        rms_velocity = double(card, 8, 'rms_velocity')
        assert gust_type in {'DRY', 'DRYP', 'VK2', 'VK2P', 'RANDOMP', 'TABLEP'}, f'gust_type={gust_type!r}'

        assert len(card) <= 9, f'len(CGUST card) = {len(card):d}\ncard={card}'
        return CGUST(cgust_id, gust_type, length_gust, one_over_velocity,
                     x0, rms_velocity, fmax=fmax, df=df, comment=comment)

    # def validate(self):
    #     assert self.true_g in ['TRUE', 'G'], 'true_g=%r' % self.true_g

    def cross_reference(self, model: BDF) -> None:
        pass

    def safe_cross_reference(self, model: BDF, xref_errors):
        self.cross_reference(model)

    def uncross_reference(self) -> None:
        """Removes cross-reference links"""
        pass

    def plot(self, fig: plt.Figure):
        # if self.mldcomd_ref is None:
        #     return
        mldcomd = self.mldcomd_ref
        ntables = len(mldcomd.extinps_ref)
        assert ntables > 0, ntables
        axes = fig.subplots(nrows=ntables)
        if ntables == 1:
            axes = [axes]
        assert len(axes) == ntables, (axes, ntables)

        for i, extinp, table in zip(count(), mldcomd.extinps_ref, mldcomd.tables_ref):
            label = extinp.label
            ax = axes[i]
            ax.grid(True)
            ax.set_ylabel(label)
            ax.set_xlabel('Time (sec)')
            ax.plot(table.x, table.y, label=label)
            ax.legend()

    def raw_fields(self):
        """
        Gets the fields in their unmodified form

        Returns
        -------
        fields : list[varies]
            the fields that define the card

        """
        list_fields = [
            'CGUST', self.cgust_id, self.gust_type,
            self.length_gust, self.one_over_velocity,
            self.x0, self.fmax, self.df, self.rms_velocity]
        return list_fields

    def repr_fields(self):
        list_fields = self.raw_fields()
        return list_fields

    def write_card(self, size: int = 8, is_double: bool = False) -> str:
        card = self.repr_fields()
        return self.comment + print_card_8(card)

MFTGUST = None
