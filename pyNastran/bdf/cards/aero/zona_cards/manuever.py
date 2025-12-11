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

from pyNastran.bdf.cards.aero.zona_cards.spline import cross_reference_set
from pyNastran.utils.numpy_utils import integer_types
from pyNastran.bdf.field_writer_8 import print_card_8
from pyNastran.bdf.cards.base_card import BaseCard
from pyNastran.bdf.bdf_interface.assign_type import (
    integer, integer_or_blank, double_or_blank, string,
    string_or_blank, double,
    integer_or_string, integer_string_or_blank,
    string_multifield, parse_components as fcomponent
)

if TYPE_CHECKING:  # pragma: no cover
    from pyNastran.bdf.bdf import BDF
    from pyNastran.bdf.bdf_interface.bdf_card import BDFCard
from pyNastran.bdf.cards.aero.zona_cards.ase import (
    ASECONT, CJUNCT, )
from pyNastran.bdf.cards.aero.zona_cards.cards import MLDCOMD



class MLOADS(BaseCard):
    """
    Defines the control system, aeroelastic system, airframe states,
    pilot input commands and time integration for transient maneuver load analysis.

    SID Unique set identification number. (Integer > 0) (See Remark 1)
    CONID Identification number of an ASECONT bulk data card specifying parameters of the
    control system. Required only for the closed-loop system. (Integer ≥ 0) (See Remark
    2)
    FLTID Identification number of a FLUTTER bulk data card specifying the flight condition
    and the associated structural and aerodynamic matrices. (Integer > 0) (See Remark 3)
    RAAID Identification number of a MINSTAT bulk data card specifying the parameters for
    rational function aerodynamic approximation. For RAAID > 0, The state space
    approach is used. For RAAID = 0 the frequency domain approach is used. (Integer ≥
    0)
    STATES Identification number of a MLDSTAT bulk data card specifying the parameters of the
    airframe states. (Integer ≥ 0)(See Remark 4)
    COMMAND Identification number of a MLDCOMD bulk data card specifying the parameters of
    the pilot’s input commands for the maneuver. (Integer > 0)
    TIME Identification number of a MLDTIME bulk data card specifying the parameters of the
    time integration for solving the transient response problem. (Integer > 0)
    MLDPRNT Identification number of a MLDPRNT bulk data card specifying the time history of
    parameters that are to be printed out. (Integer ≥ 0)
    FMAX Maximum frequency in cycles/sec. Used only for frequency domain approach or
    summation of forces by inverse Fourier transform to obtain the time domain
    aerodynamic forces (SOF = “YES” in MLDPRNT bulk data card). (Real ≥ 0.0 ) (See
    Remark 5).
    DF Frequency step size in cycles/sec. Note that based on the Nyquist criterion, the
    maximum time for the time domain response (specified by the entry TEND in the
    MLDTIME bulk data card) must be less than 1/DF. Otherwise, the time-domain
    response becomes periodic even the input excitation is a decay function. (Real > 0.0,
    default = 0.01Hz).
    """
    type = 'MLOADS'
    # _field_map = {
    #     1: 'sid', 2: 'mach', 3: 'q', 8: 'aeqr',
    # }

    def __init__(self, mloads_id: int, asecont_id: int,
                 flutter_id: int, minstat_id: int,
                 mldstat_id: int, mldcomd_id: int, mldtime_id: int,
                 mldprint_id: int,
                 fmax: float, save_freq: str, filename: str,
                 df: float=0.01, comment: str=''):
        BaseCard.__init__(self)
        if comment:
            self.comment = comment

        self.mloads_id = mloads_id
        assert mloads_id > 0, mloads_id
        self.asecont_id = asecont_id
        self.flutter_id = flutter_id
        self.minstat_id = minstat_id

        self.mldstat_id = mldstat_id
        self.mldcomd_id = mldcomd_id
        self.mldtime_id = mldtime_id
        self.mldprint_id = mldprint_id

        self.df = df
        self.fmax = fmax
        self.save_freq = save_freq
        self.filename = filename

        self.asecont_ref = None
        self.flutter_ref = None
        self.minstat_ref = None
        self.mldstat_ref = None
        self.mldcomd_ref = None
        self.mldtime_ref = None
        self.mldprint_ref = None

    @classmethod
    def add_card(cls, card: BDFCard, comment: str=''):
        """
        Adds a MLOADS card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card

        """
        # MLOADS SID  CONTID FLTID   RAAID STATES COMMAND TIME MLDPRNT CONT
        #        FMAX DF     SAVEFRQ --FILEFRQ---
        mloads_id = integer(card, 1, 'mloads_id')
        asecont_id = integer_or_blank(card, 2, 'asecont_id', default=0)
        flutter_id = integer_or_blank(card, 3, 'flutter_id', default=0)
        minstat_id = integer_or_blank(card, 4, 'raaid/minstat_id', default=0)

        mldstat_id = integer_or_blank(card, 5, 'states/mldstat_id', default=0)
        mldcomd_id = integer(card, 6, 'command/mldcomd_id')
        mldtime_id = integer(card, 7, 'time/mldtime_id')
        mldprint_id = integer_or_blank(card, 8, 'mldprint_id', default=0)

        fmax = double_or_blank(card, 9, 'fmax', default=None)
        df = double_or_blank(card, 10, 'df', default=0.01)
        save_freq = string_or_blank(card, 11, 'save_freq', default='')
        assert save_freq in {'SAVE', 'ACQUIRE', ''}, save_freq
        filename = ''
        if len(card) > 12:
            filename = string_multifield(card, (12, 13), 'filename')

        assert 7 <= len(card) < 14, f'len(MLOADS card) = {len(card):d}\ncard={card}'
        return MLOADS(mloads_id, asecont_id, flutter_id, minstat_id,
                      mldstat_id, mldcomd_id, mldtime_id,
                      mldprint_id, fmax, save_freq, filename, df=df, comment=comment)

    # def validate(self):
    #     assert self.true_g in ['TRUE', 'G'], 'true_g=%r' % self.true_g

    def cross_reference(self, model: BDF) -> None:
        zona = model.zona
        if self.asecont_id > 0:
            self.asecont_ref = zona.asecont[self.asecont_id]
        if self.flutter_id > 0:
            self.flutter_ref = model.flutters[self.flutter_id]
        # self.minstat_ref = zona.minstat[self.minstat_id]

        if self.mldstat_id > 0:
            self.mldstat_ref = zona.mldstat[self.mldstat_id]
        self.mldcomd_ref = zona.mldcomd[self.mldcomd_id]
        # self.mldtime_ref = zona.mldtime[self.mldtime_id]
        # self.mldprint_ref = zona.mldprnt[self.mldprint_id]

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
        self.mldprint_ref = None

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
            'MLOADS', self.mloads_id, self.asecont_id,
            self.flutter_id, self.minstat_id,
            self.mldstat_id, self.mldcomd_id, self.mldtime_id,
            self.mldprint_id, self.fmax, self.df, self.save_freq,
            self.filename]
        return list_fields

    def repr_fields(self):
        list_fields = self.raw_fields()
        return list_fields

    def write_card(self, size: int=8, is_double: bool=False) -> str:
        card = self.repr_fields()
        return self.comment + print_card_8(card)


class EXTINP(BaseCard):
    type = 'EXTINP'
    # _field_map = {
    #     1: 'sid', 2: 'mach', 3: 'q', 8: 'aeqr',
    # }

    def __init__(self, extinp_id: int, input_type: int,
                 itf_id: int, itf_component: int,
                 label: str, comment: str=''):
        BaseCard.__init__(self)
        if comment:
            self.comment = comment

        self.extinp_id = extinp_id
        self.input_type = input_type
        self.itf_id = itf_id
        self.itf_component = itf_component
        self.label = label
        self.itf_ref = None

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
        # EXTINP ID  TYPE ITFID CI LABEL
        # EXTINP 100      400   1  PILOT
        extinp_id = integer(card, 1, 'extinp_id')
        input_type = string_or_blank(card, 2, 'input_type', default='')
        itf_id = integer(card, 3, 'asecont_id')
        itf_component = integer(card, 4, 'itf_component')
        assert itf_component in {1, 2, 3, 4, 5, 6}, itf_component
        label = string(card, 5, 'label')
        assert len(card) == 6, f'len(EXTINP card) = {len(card):d}\ncard={card}'
        return EXTINP(extinp_id, input_type, itf_id, itf_component, label, comment=comment)

    # def validate(self):
    #     assert self.true_g in ['TRUE', 'G'], 'true_g=%r' % self.true_g

    def cross_reference(self, model: BDF, parent: ASECONT | MLDCOMD) -> None:
        self.itf_ref = get_external_obj(self, model, self.itf_id, parent)
        # self.itf_id = itf_id
        # self.itf_component = itf_component
        # self.asecont_ref = model.asecont[self.asecont_id]
        # self.flutter_ref = model.flutter[self.flutter_id]
        # self.minstat_ref = model.minstat[self.minstat_id]
        #
        # self.mldstat_ref = model.mldstat[self.mldstat_id]
        # self.mldcomd_ref = model.mldcomd[self.mldcomd_id]
        # self.mldtime_ref = model.mldtime[self.mldtime_id]
        # self.mldprint_ref = model.mldprnt[self.mldprint_id]

    def safe_cross_reference(self, model: BDF,
                             parent: ASECONT | MLDCOMD):
        self.cross_reference(model, parent)

    def uncross_reference(self) -> None:
        """Removes cross-reference links"""
        self.itf_ref = None

    def raw_fields(self):
        """
        Gets the fields in their unmodified form

        Returns
        -------
        fields : list[varies]
            the fields that define the card

        """
        list_fields = ['EXTINP', self.extinp_id, self.input_type,
                       self.itf_id, self.itf_component, self.label]
        return list_fields

    def repr_fields(self):
        list_fields = self.raw_fields()
        return list_fields

    def write_card(self, size: int=8, is_double: bool=False) -> str:
        card = self.repr_fields()
        return self.comment + print_card_8(card)


class EXTOUT(BaseCard):
    type = 'EXTOUT'
    # _field_map = {
    #     1: 'sid', 2: 'mach', 3: 'q', 8: 'aeqr',
    # }

    def __init__(self, extout_id: int, input_type: int,
                 itf_id: int, itf_component: int,
                 label: str, comment: str=''):
        BaseCard.__init__(self)
        if comment:
            self.comment = comment

        self.extout_id = extout_id
        self.input_type = input_type
        self.itf_id = itf_id
        self.itf_component = itf_component
        self.label = label
        self.itf_ref = None

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
        # EXTOUT ID  TYPE ITFID CI LABEL
        # EXTOUT 100      400   1  PILOT
        extout_id = integer(card, 1, 'extout_id')
        input_type = string_or_blank(card, 2, 'input_type', default='')
        itf_id = integer(card, 3, 'asecont_id')
        itf_component = integer(card, 4, 'itf_component')
        assert itf_component in {1, 2, 3, 4, 5, 6}, itf_component
        label = string(card, 5, 'label')
        assert len(card) == 6, f'len(EXTOUT card) = {len(card):d}\ncard={card}'
        return EXTOUT(extout_id, input_type, itf_id, itf_component, label, comment=comment)

    # def validate(self):
    #     assert self.true_g in ['TRUE', 'G'], 'true_g=%r' % self.true_g

    def cross_reference(self, model: BDF, parent: CJUNCT | MLDCOMD) -> None:
        self.itf_ref = get_external_obj(self, model, self.itf_id, parent)

    def safe_cross_reference(self, model: BDF,
                             parent: ASECONT | MLDCOMD):
        self.cross_reference(model)

    def uncross_reference(self) -> None:
        """Removes cross-reference links"""
        self.itf_ref = None

    def raw_fields(self):
        """
        Gets the fields in their unmodified form

        Returns
        -------
        fields : list[varies]
            the fields that define the card

        """
        list_fields = ['EXTOUT', self.extout_id, self.input_type,
                       self.itf_id, self.itf_component, self.label]
        return list_fields

    def repr_fields(self):
        list_fields = self.raw_fields()
        return list_fields

    def write_card(self, size: int=8, is_double: bool=False) -> str:
        card = self.repr_fields()
        return self.comment + print_card_8(card)


class ACTU(BaseCard):
    type = 'ACTU'
    # _field_map = {
    #     1: 'sid', 2: 'mach', 3: 'q', 8: 'aeqr',
    # }

    def __init__(self, actu_id: int,
                 a0: float, a1: float, a2: float, comment: str=''):
        BaseCard.__init__(self)
        if comment:
            self.comment = comment

        self.actu_id = actu_id
        self.a0 = a0
        self.a1 = a1
        self.a2 = a2

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
        # ACTU ID A0  A1  A2
        # ACTU 10 0.8 0.5 0.2
        actu_id = integer(card, 1, 'actu_id')
        a0 = double_or_blank(card, 2, 'a0', default=0.0)
        a1 = double_or_blank(card, 3, 'a1', default=0.0)
        a2 = double_or_blank(card, 4, 'a2', default=0.0)
        assert 2 <= len(card) <= 5, f'len(ACTU card) = {len(card):d}\ncard={card}'
        return ACTU(actu_id, a0, a1, a2, comment=comment)

    # def validate(self):
    #     assert self.true_g in ['TRUE', 'G'], 'true_g=%r' % self.true_g

    def cross_reference(self, model: BDF) -> None:
        pass

    def safe_cross_reference(self, model: BDF, xref_errors):
        self.cross_reference(model)

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
        list_fields = ['ACTU', self.actu_id, self.a0, self.a1, self.a2]
        return list_fields

    def repr_fields(self):
        list_fields = self.raw_fields()
        return list_fields

    def write_card(self, size: int=8, is_double: bool=False) -> str:
        card = self.repr_fields()
        return self.comment + print_card_8(card)


ALLOWED_TRIMFNC_LABELS = {
    'AERO': {
        'CDL',  # induced drag
        'CP',  # coefficient of pressure
        'CL', 'CY',
        'CR', 'CM', 'CN',  # CR is CL (moment)
        'NX', 'NY', 'NZ',
        'PDOT', 'QDOT', 'RDOT', 'LOADMOD', 'TRIMVAR'},
    'FEM': {
        'LOADMOD', 'LOADMOD1', 'GRIDDISP', 'FORCE'},
    'MODAL': {'AEFACT', 'DMI', 'PCHFILE'},
}

class TRIMFNC(BaseCard):
    type = 'TRIMFNC'
    # _field_map = {
    #     1: 'sid', 2: 'mach', 3: 'q', 8: 'aeqr',
    # }

    def __init__(self, trimfnc_id: int, fcn_type: str, label: str,
                 rhs_flag: str,
                 is_set: int | str,
                 ia_set: int | str,
                 remark: str, comment: str=''):
        BaseCard.__init__(self)
        if comment:
            self.comment = comment

        self.trimfnc_id = trimfnc_id
        self.fcn_type = fcn_type
        self.label = label
        self.rhs_flag = rhs_flag
        self.is_set = is_set
        self.ia_set = ia_set
        self.remark = remark
        self.is_set_ref = None
        self.ia_set_ref = None

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
        # TRIMFNC IDFNC TYPE  LABEL RHS ISSET   IASET   REMARK
        # TRIMFNC 10    MODAL DMI   LHS MATRIXR MATRIXL STRESS.AT.CBAR
        trimfnc_id = integer(card, 1, 'trimfnc_id')
        fcn_type = string(card, 2, 'fcn_type')
        label = string(card, 3, 'label')
        rhs_flag = string_or_blank(card, 4, 'rhs_flag', default='')
        is_set = integer_string_or_blank(card, 5, 'is_set', default='')
        ia_set = integer_string_or_blank(card, 6, 'ia_set', default='')

        remark = string_multifield(card, (7, 8), 'remark')

        allowed_trim_fnc_labels = ALLOWED_TRIMFNC_LABELS[fcn_type]
        assert fcn_type in ['FEM', 'AERO', 'MODAL'], fcn_type
        assert label in allowed_trim_fnc_labels, f'label={label!r}'
        assert len(card) <= 9, f'len(TRIMFNC card) = {len(card):d}\ncard={card}'
        return TRIMFNC(trimfnc_id, fcn_type, label, rhs_flag, is_set, ia_set, remark, comment=comment)

    # def validate(self):
    #     assert self.true_g in ['TRUE', 'G'], 'true_g=%r' % self.true_g

    def cross_reference(self, model: BDF) -> None:
        zona = model.zona
        if self.fcn_type == 'AERO' and self.label == 'TRIMVAR':
            self.is_set_ref = zona.trimvar[self.is_set]
            # if self.ia_set:
            #     self.ia_set_ref = zona.trimvar[self.ia_set]
        elif self.fcn_type == 'AERO' and self.label == 'LOADMOD':
            self.is_set_ref = zona.loadmod[self.is_set]
            # if self.ia_set:
            #     self.ia_set_ref = zona.trimvar[self.ia_set]
        elif self.fcn_type == 'FEM' and self.label in ['LOADMOD', 'LOADMOD1']:
            self.is_set_ref = zona.loadmod[self.is_set]
        elif self.fcn_type == 'MODAL' and self.label == 'DMI':
            self.is_set_ref = model.dmi[self.is_set]
        elif self.fcn_type == 'AERO' and self.label in {'CP', 'CDL',
                                                        'CY', 'CL',
                                                        'CR', 'CM', 'CN',
                                                        'NX', 'NY', 'NZ',
                                                        'PDOT', 'QDOT', 'RDOT'}:
            self.is_set_ref = None
        else:
            raise RuntimeError(f'fcn_type={self.fcn_type!r}, label={self.label!r} not implemented')

        # if isinstance(self.is_set, integer_types):
        #     self.is_set_ref = model.aefacts[self.is_set]
        # else:
        #     self.is_set_ref = model.dmi[self.is_set]

        # if isinstance(self.ia_set, integer_types):
        #     self.ia_set_ref = model.aefacts[self.ia_set]
        # else:
        #     self.ia_set_ref = model.dmi[self.is_set]

    def safe_cross_reference(self, model: BDF, xref_errors):
        self.cross_reference(model)

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
            'TRIMFNC', self.trimfnc_id, self.fcn_type, self.label,
            self.rhs_flag, self.is_set, self.ia_set, self.remark]
        return list_fields

    def repr_fields(self):
        list_fields = self.raw_fields()
        return list_fields

    def write_card(self, size: int=8, is_double: bool=False) -> str:
        msg = f'TRIMFNC {self.trimfnc_id:<8d}{self.fcn_type:<8s}{self.label:<8s}'\
              f'{self.rhs_flag:<8s}{str(self.is_set):<8s}{str(self.ia_set):<8s}{self.remark}\n'
        card = self.repr_fields()
        # return self.comment + print_card_8(card)
        return self.comment + msg


class LOADMOD(BaseCard):
    type = 'LOADMOD'
    # _field_map = {
    #     1: 'sid', 2: 'mach', 3: 'q', 8: 'aeqr',
    # }

    def __init__(self, loadmod_id: int, label, cid: int,
                 set_k: int, set_g: int,
                 comment: str=''):
        BaseCard.__init__(self)
        if comment:
            self.comment = comment

        self.loadmod_id = loadmod_id
        self.label = label
        self.cid = cid
        self.set_k = set_k
        self.set_g = set_g
        assert set_k + set_g > 0, (set_k, set_g)
        self.cid_ref = None
        self.set_k_ref = None
        self.set_g_ref = None

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
        # LOADMOD LID LABEL  CP SETK SETG
        # LOADMOD 10  XSHEAR 1  1
        loadmod_id = integer(card, 1, 'loadmod_id')
        label = string(card, 2, 'label')
        cid = integer_or_blank(card, 3, 'cp/cid', default=0)
        set_k = integer_or_blank(card, 4, 'set_k', default=0)
        set_g = integer_or_blank(card, 5, 'set_g', default=0)
        assert len(card) == 6, f'len(LOADMOD card) = {len(card):d}\ncard={card}'
        return LOADMOD(loadmod_id, label, cid, set_k, set_g, comment=comment)

    # def validate(self):
    #     assert self.true_g in ['TRUE', 'G'], 'true_g=%r' % self.true_g

    def cross_reference(self, model: BDF) -> None:
        which_msg = f', which is required by LOADMOD={self.loadmod_id}'
        self.cid_ref = model.Coord(self.cid, which_msg)
        zona = model.zona
        if self.set_k > 0:
            self.set_k_ref = zona.panlsts[self.set_k]
        if self.set_g > 0:
            msg = f'LOADMOD={self.loadmod_id}'
            self.set_g_ref = cross_reference_set(model, self.set_g, msg, which_msg)

    def safe_cross_reference(self, model: BDF, xref_errors):
        self.cross_reference(model)

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
            'LOADMOD', self.loadmod_id, self.label, self.cid,
            self.set_k, self.set_g]
        return list_fields

    def repr_fields(self):
        list_fields = self.raw_fields()
        return list_fields

    def write_card(self, size: int=8, is_double: bool=False) -> str:
        card = self.repr_fields()
        return self.comment + print_card_8(card)


class RBRED(BaseCard):
    type = 'RBRED'

    # _field_map = {
    #     1: 'sid', 2: 'mach', 3: 'q', 8: 'aeqr',
    # }

    def __init__(self, sid: int, id_ase: int,
                 component: str, node_id: int,
                 phugoid0: str, comment: str = ''):
        BaseCard.__init__(self)
        if comment:
            self.comment = comment

        self.sid = sid
        self.id_ase = id_ase
        self.component = component
        self.node_id = node_id
        self.phugoid0 = phugoid0
        self.id_ase_ref = None

    @classmethod
    def add_card(cls, card: BDFCard, comment: str = ''):
        """
        Adds a RBRED card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card

        """
        # RBRED SID IDASE C   GID PHUGO
        # RBRED 10  200   246 10
        sid = integer(card, 1, 'sid')
        id_ase = integer(card, 2, 'id_ase')
        component = fcomponent(card, 3, 'component')
        node_id = integer_or_blank(card, 4, 'node_id', default=0)
        phugoid0 = string_or_blank(card, 5, 'phugoid0', default='NO')
        assert 4 <= len(card) <= 6, f'len(RBRED) = {len(card):d}\ncard={card}'
        return RBRED(sid, id_ase, component, node_id, phugoid0, comment=comment)

    # def validate(self):
    #     assert self.true_g in ['TRUE', 'G'], 'true_g=%r' % self.true_g

    def cross_reference(self, model: BDF) -> None:
        # msg = f', which is required by RBRED={self.sid}'
        # ASE, MLOADS, ELOADS, GLOADS, DFS, or NLFLTR
        zona = model.zona
        if self.id_ase in zona.ase:
            self.id_ase_ref = zona.ase[self.id_ase]
        elif self.id_ase in zona.mloads:
            self.id_ase_ref = zona.mloads[self.id_ase]
        elif self.id_ase in zona.eloads:
            self.id_ase_ref = zona.eloads[self.id_ase]
        elif self.id_ase in zona.gloads:
            self.id_ase_ref = zona.gloads[self.id_ase]
        else:  # pragma: no cover
            ase = list(zona.ase)
            mloads = list(zona.mloads)
            eloads = list(zona.eloads)
            gloads = list(zona.gloads)
            dfs = list(zona.dfs)
            nlfltr = list(zona.nlfltr)
            msg = (
                f'{self.id_ase} is not in [ASE, MLOADS, ELOADS, GLOADS, DFS, NLFLTR]\n'
                f' - ase    = {ase}\n'
                f' - mloads = {mloads}\n'
                f' - eloads = {eloads}\n'
                f' - gloads = {gloads}\n'
                f' - dfs    = {dfs}\n'
                f' - nlfltr = {nlfltr}\n')
            raise RuntimeError(msg)
        # self.node_ref = model.Node(self.id_ase, msg)

    def safe_cross_reference(self, model: BDF, xref_errors):
        self.cross_reference(model)

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
            'RBRED', self.sid, self.id_ase,
            self.component, self.node_id, self.phugoid0]
        return list_fields

    def repr_fields(self):
        list_fields = self.raw_fields()
        return list_fields

    def write_card(self, size: int = 8, is_double: bool = False) -> str:
        card = self.repr_fields()
        return self.comment + print_card_8(card)


def get_external_obj(obj: EXTINP | EXTOUT, model: BDF,
                     itf_id: int, parent: MLDCOMD):
    """
    If there is no ASECONT bulk data card specified, and the EXTINP bulk data card is referred to by
    a MLDCOMD bulk data card, the ITFID entry must refer to an ACTU bulk data card (associated
    with aerodynamic control surface) or a SISOTF bulk data card (associated with GRIDFRC bulk
    data card). Also, CI = 1 is required.
    """
    zona = model.zona

    asecont = list(zona.asecont)
    actu = list(zona.actu)
    sisotf = list(zona.sisotf)
    cjunct = list(zona.cjunct)
    msg = (
        f'itf_id={itf_id} is not in the model\n'
        f' - asecont = {asecont}\n'
        f' - actu    = {actu}\n'
        f' - sisotf  = {sisotf}\n'
        f' - cjunct  = {cjunct}\n'
    )
    if itf_id in zona.asecont:
        itf_ref = zona.asecont[itf_id]
    elif obj.type == 'EXTINP':
        # If there is:
        #  - no ASECONT card specified
        #  - the EXTINP card is referred to by an MLDCOMD card
        #  -> ITFID entry must refer to an ACTU or SISOTF
        if parent.type == 'MLDCOMD' and len(zona.asecont) == 0:
            if itf_id in zona.actu:
                itf_ref = zona.actu[itf_id]
                assert obj.itf_component == 1, f'component={obj.itf_component}\n{str(obj)}\nparent:\n{str(parent)}'
            elif itf_id in zona.sisotf:
                itf_ref = zona.sisotf[itf_id]
            else:
                raise RuntimeError(msg)
        elif itf_id in zona.cjunct:
            itf_ref = zona.cjunct[itf_id]
        else:
            raise RuntimeError(msg)
    else:
        assert parent.type in ['ASECONT', 'MLDCOMD'], parent
        # if itf_id in zona.actu:
        #     itf_ref = zona.actu[itf_id]
        # elif itf_id in zona.sisotf:
        #     itf_ref = zona.sisotf[itf_id]
        if itf_id in zona.cjunct:
            itf_ref = zona.cjunct[itf_id]
        else:
            raise RuntimeError(msg)
    return itf_ref
