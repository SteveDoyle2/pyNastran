# coding: utf-8
# pylint: disable=W0212,C0103
"""
All ZONA aero cards are defined in this file.  This includes:
 * ATM
 * FIXMATM

All cards are BaseCard objects.

"""
from __future__ import annotations
from typing import TYPE_CHECKING

import numpy as np
from matplotlib import pyplot as plt

from pyNastran.utils.convert import convert_length
from pyNastran.utils.numpy_utils import integer_types, float_types
from pyNastran.bdf.field_writer_8 import print_card_8
from pyNastran.bdf.cards.base_card import BaseCard
from pyNastran.bdf.bdf_interface.assign_type import (
    integer, integer_or_blank, double_or_blank, string,
    blank, double, string_or_blank,
    # integer_or_string,
    # string_or_blank,
    # integer_or_double,
    # integer_or_string, integer_string_or_blank,
    # string_multifield,  # parse_components as fcomponent
)
from pyNastran.bdf.cards.aero.zaero_interface.get_card import get_atmos, get_mkaeroz
if TYPE_CHECKING:  # pragma: no cover
    from pyNastran.bdf.bdf import BDF
    from pyNastran.bdf.bdf_interface.bdf_card import BDFCard

class ATMOS(BaseCard):
    type = 'ATMOS'

    # _field_map = {
    #     1: 'sid', 2: 'mach', 3: 'q', 8: 'aeqr',
    # }

    def __init__(self, atmos_id: int,
                 mass_unit: str,
                 length_unit: str,
                 temperature_unit: str,
                 atmosphere_table: list[float], comment: str=''):
        BaseCard.__init__(self)
        if comment:
            self.comment = comment

        self.atmos_id = atmos_id
        length_unit = length_unit.upper()
        self.length_unit = length_unit
        self.temperature_unit = temperature_unit
        self.atmosphere_table = np.asarray(atmosphere_table)
        if self.atmosphere_table.ndim == 1:
            nalt = len(self.atmosphere_table) // 4
            self.atmosphere_table = self.atmosphere_table.reshape(nalt, 4)
        else:
            assert self.atmosphere_table.ndim == 2, self.atmosphere_table.shape
        MASS_MAP = {
            'SLINCH': 'SLIN',
            'LBF': 'LBF/',
            'N': 'N/',
        }
        mass_unit = MASS_MAP.get(mass_unit.upper(), mass_unit)
        self.mass_unit = mass_unit

        assert mass_unit in ['SLIN', 'SLUG', 'LBM', 'G', 'KG', 'LBF/', 'N/', 'NONE'], f'mass_unit={mass_unit}'
        assert length_unit in ['IN', 'FT', 'M', 'MM', 'CM', 'KM', 'NONE'], f'length_unit={length_unit}'
        assert temperature_unit in ['R', 'K', 'C', 'F'], temperature_unit
        assert len(atmosphere_table) > 0, atmosphere_table
        alt = self.alt
        sos = self.sos
        assert len(np.unique(alt)) == len(alt)
        assert np.array_equal(np.unique(alt), alt)
        # assert np.array_equal(np.unique(sos), sos)

    @property
    def alt(self) -> np.ndarray:
        return self.atmosphere_table[:, 0]

    @property
    def sos(self) -> np.ndarray:
        return self.atmosphere_table[:, 1]

    @property
    def density(self) -> np.ndarray:
        return self.atmosphere_table[:, 2]

    @property
    def temperature(self) -> np.ndarray:
        return self.atmosphere_table[:, 3]

    @classmethod
    def add_card(cls, card: BDFCard, comment: str=''):
        """
        Adds a ATMOS card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card

        """
        # ATMOS IDATM AMMUNIT AMLUNIT AMTUNIT
        #       ALT1  SOUND1  DEN1    TEMP1  ALT2 SOUND2 DEN2 TEMP2
        #       ALTi  SOUNDi  DENi    TEMPi -etc-
        # FIXMATM 100     10    12    slug ft 1.0 -1 +FIX1
        #         -10000. 0. 10000. 20000. 30000.
        atmos_id = integer(card, 1, 'atmos_id')
        mass_unit = string(card, 2, 'mass_unit')
        length_unit = string(card, 3, 'length_unit')
        temperature_unit = string(card, 4, 'temperature_unit')

        atmosphere_table = []
        j = 1
        for ifield in range(9, len(card), 4):
            alt = double(card, ifield, f'alt{j+1}')
            sos = double(card, ifield+1, f'sound{j+1}')
            rho = double(card, ifield+2, f'density{j+1}')
            temp = double(card, ifield+3, f'temperature{j+1}')
            atmosphere_table.append([alt, sos, rho, temp])
            j += 1
        assert len(card) > 8, f'len(ATMOS card) = {len(card):d}\ncard={card}'
        return ATMOS(atmos_id, mass_unit, length_unit,
                     temperature_unit, atmosphere_table, comment=comment)

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
        ntables = 1
        assert ntables > 0, ntables
        axes = fig.subplots(nrows=ntables)
        if ntables == 1:
            axes = [axes]
        assert len(axes) == ntables, (axes, ntables)

    def raw_fields(self):
        """
        Gets the fields in their unmodified form

        Returns
        -------
        fields : list[varies]
            the fields that define the card

        """
        assert self.temperature_unit in ['R', 'K', 'F', 'C'], self.temperature_unit
        list_fields = [
            'ATMOS', self.atmos_id, self.mass_unit, self.length_unit, self.temperature_unit,
            None, None, None, None,
        ] + self.atmosphere_table.ravel().tolist()
        return list_fields

    def repr_fields(self) -> list:
        list_fields = self.raw_fields()
        return list_fields

    def write_card(self, size: int = 8, is_double: bool = False) -> str:
        card = self.repr_fields()
        return self.comment + print_card_8(card)


class FIXHATM(BaseCard):
    type = 'FIXHATM'

    # _field_map = {
    #     1: 'sid', 2: 'mach', 3: 'q', 8: 'aeqr',
    # }

    def __init__(self, sid: int, alt: float, atm_id: int,
                 mass_unit: str, length_unit: str,
                 fluttf_id: int, print_flag: int, mkaeroz_ids: list[int],
                 vref: float=1.0, comment: str=''):
        BaseCard.__init__(self)
        if comment:
            self.comment = comment

        self.sid = sid
        self.mkaeroz_ids = mkaeroz_ids
        self.atm_id = atm_id
        self.mass_unit = mass_unit
        self.length_unit = length_unit
        self.vref = vref
        self.fluttf_id = fluttf_id
        self.print_flag = print_flag
        self.alt = alt
        self.atmos_ref = None
        self.mkaerozs_ref = None

    @classmethod
    def add_card(cls, card: BDFCard, comment: str = ''):
        """
        Adds a FIXMATM card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card

        """
        # FIXMATM SETID   IDMK  IDATM FTMUNIT FTLUNIT VREF FLUTTF PRINT CONT
        #         ALT1    ALT2  ALTi -etc-
        # FIXMATM 100     10    12    slug ft 1.0 -1 +FIX1
        #         -10000. 0. 10000. 20000. 30000.
        sid = integer(card, 1, 'cgust_id')
        alt = double(card, 2, 'alt')
        atm_id = integer(card, 3, 'atm_id')
        mass_unit = string(card, 4, 'mass_unit')
        length_unit = string(card, 5, 'length_unit')
        vref = double_or_blank(card, 6, 'vref', default=1.0)
        fluttf_id = integer_or_blank(card, 7, 'fluttf_id', default=0)
        print_flag = integer(card, 8, 'print_flag')

        j = 1
        mkaeroz_ids = []
        for ifield in range(9, len(card)):
            mkaeoz_id = integer(card, ifield, f'mkaeoz_id{j}')
            mkaeroz_ids.append(mkaeoz_id)
            j += 1
        assert len(card) > 8, f'len(FIXEMATM card) = {len(card):d}\ncard={card}'
        return FIXHATM(sid, alt, atm_id, mass_unit,
                       length_unit, fluttf_id, print_flag, mkaeroz_ids,
                       vref=vref, comment=comment)

    # def validate(self):
    #     assert self.true_g in ['TRUE', 'G'], 'true_g=%r' % self.true_g

    def cross_reference(self, model: BDF) -> None:
        if self.atm_id:
            self.atmos_ref = model.zaero.atmos[self.atm_id]
        if self.fluttf_id and 0:
            self.fluttf_ref = model.zaero.fluttf[self.fluttf_id]
        mkaerozs_ref = []
        for mkaeroz_id in self.mkaeroz_ids:
            mkaeroz_ref = model.zaero.mkaeroz[mkaeroz_id]
            mkaerozs_ref.append(mkaeroz_ref)
        self.mkaerozs_ref = mkaerozs_ref

    def safe_cross_reference(self, model: BDF, xref_errors):
        self.cross_reference(model)

    def uncross_reference(self) -> None:
        """Removes cross-reference links"""
        self.atmos_ref = None
        self.mkaerozs_ref = None

    def plot(self, fig: plt.Figure):
        ntables = 1
        assert ntables > 0, ntables
        axes = fig.subplots(nrows=ntables)
        if ntables == 1:
            axes = [axes]
        assert len(axes) == ntables, (axes, ntables)


    def raw_fields(self):
        """
        Gets the fields in their unmodified form

        Returns
        -------
        fields : list[varies]
            the fields that define the card

        """
        list_fields = [
            'FIXMATM', self.sid, self.alt, self.atm_id, self.mass_unit,
            self.length_unit, self.vref, self.fluttf_id, self.print_flag] + self.mkaeroz_ids
        return list_fields

    def repr_fields(self):
        list_fields = self.raw_fields()
        return list_fields

    def write_card(self, size: int = 8, is_double: bool = False) -> str:
        card = self.repr_fields()
        return self.comment + print_card_8(card)


class FIXMATM(BaseCard):
    type = 'FIXMATM'

    # _field_map = {
    #     1: 'sid', 2: 'mach', 3: 'q', 8: 'aeqr',
    # }

    def __init__(self, sid: int, mkaeroz_id: int,
                 mass_unit: str, length_unit: str,
                 alts: list[float],
                 atm_id: int = 0,
                 fluttf_id: int=0, print_flag: int=0,
                 vref: float = 1.0, comment: str=''):
        BaseCard.__init__(self)
        if comment:
            self.comment = comment

        self.sid = sid
        self.mkaeroz_id = mkaeroz_id
        self.atm_id = atm_id
        self.mass_unit = mass_unit
        self.length_unit = length_unit
        self.vref = vref
        self.fluttf_id = fluttf_id
        self.print_flag = print_flag
        self.alts = np.asarray(alts)
        self.atmos_ref = None
        self.fluttf_ref = None
        self.mkaeroz_ref = None
        assert isinstance(mkaeroz_id, integer_types), self.get_stats()
        assert len(alts) >= 0, alts
        if len(alts) > 1:
            assert alts[0] < alts[1], alts
        if atm_id == 0:
            # -100 kft to 260 kft
            alts2 = convert_length(self.alts, self.length_unit.lower(), 'ft')
            assert alts2.min() >= -100_000., (alts2.min(), alts2.max())  # ft
            assert alts2.max() <= 260_000., (alts2.min(), alts2.max())  # ft
        assert isinstance(alts[0], float_types), alts

    @classmethod
    def add_card(cls, card: BDFCard, comment: str=''):
        """
        Adds a FIXMATM card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card

        """
        # FIXMATM SETID   IDMK  IDATM FTMUNIT FTLUNIT VREF FLUTTF PRINT CONT
        #         ALT1    ALT2  ALTi -etc-
        # FIXMATM 100     10    12    slug ft 1.0 -1 +FIX1
        #         -10000. 0. 10000. 20000. 30000.
        sid = integer(card, 1, 'sid')
        mkaeroz_id = integer(card, 2, 'mkaeroz_id')
        atm_id = integer_or_blank(card, 3, 'atm_id', default=0)
        mass_unit = string(card, 4, 'mass_unit')
        length_unit = string(card, 5, 'length_unit')
        vref = double_or_blank(card, 6, 'vref', default=1.0)
        fluttf_id = integer_or_blank(card, 7, 'fluttf_id', default=0)
        print_flag = integer_or_blank(card, 8, 'print_flag', default=0)

        alts = []
        j = 1
        for ifield in range(9, len(card)):
            alt = double(card, ifield, f'alt{j}')
            alts.append(alt)
            j += 1
        assert len(card) > 8, f'len(FIXEMATM card) = {len(card):d}\ncard={card}'
        assert isinstance(mkaeroz_id, integer_types)
        return FIXMATM(sid, mkaeroz_id, mass_unit, length_unit, alts,
                       atm_id=atm_id,
                       print_flag=print_flag, fluttf_id=fluttf_id,
                       vref=vref, comment=comment)

    # def validate(self):
    #     assert self.true_g in ['TRUE', 'G'], 'true_g=%r' % self.true_g

    def cross_reference(self, model: BDF) -> None:
        cross_reference_atmos_mkaeroz_fluttf(self, model)
        # self.plot_atmosphere()

    def plot_atmosphere(self) -> None:
        atmos = self.atmos_ref
        from matplotlib import pyplot as plt
        fig = plt.figure(1)
        ax = fig.gca()
        ax.plot(atmos.alt, atmos.sos, 'k-')

        sos = np.interp(self.alts, atmos.alt, atmos.sos)
        ax.plot(self.alts, sos, 'o')
        ax.grid(True)
        ax.set_xlabel('Altitude (ft)')
        ax.set_ylabel('SO (ft/s)')
        plt.show()

    def safe_cross_reference(self, model: BDF, xref_errors):
        self.cross_reference(model)

    def uncross_reference(self) -> None:
        """Removes cross-reference links"""
        self.atmos_ref = None
        self.fluttf_ref = None
        self.mkaeroz_ref = None

    def plot(self, fig: plt.Figure):
        ntables = 1
        assert ntables > 0, ntables
        axes = fig.subplots(nrows=ntables)
        if ntables == 1:
            axes = [axes]
        assert len(axes) == ntables, (axes, ntables)

    def raw_fields(self):
        """
        Gets the fields in their unmodified form

        Returns
        -------
        fields : list[varies]
            the fields that define the card

        """
        list_fields = [
            'FIXMATM', self.sid, self.mkaeroz_id, self.atm_id, self.mass_unit,
            self.length_unit, self.vref, self.fluttf_id, self.print_flag] + self.alts.tolist()
        return list_fields

    def repr_fields(self):
        list_fields = self.raw_fields()
        return list_fields

    def write_card(self, size: int = 8, is_double: bool = False) -> str:
        card = self.repr_fields()
        return self.comment + print_card_8(card)


class FIXMACH(BaseCard):
    type = 'FIXMACH'

    # _field_map = {
    #     1: 'sid', 2: 'mach', 3: 'q', 8: 'aeqr',
    # }

    def __init__(self, sid: int, mkaeroz_id: int,
                 mass_unit: str, length_unit: str,
                 fluttf_id: int, print_flag: int,
                 velocity: list[float], rho: list[float],
                 vref: float=1.0, comment: str=''):
        BaseCard.__init__(self)
        if comment:
            self.comment = comment

        self.sid = sid
        self.mkaeroz_id = mkaeroz_id
        self.mass_unit = mass_unit
        self.length_unit = length_unit
        assert isinstance(mass_unit, str), mass_unit
        assert isinstance(length_unit, str), length_unit
        self.vref = vref
        self.fluttf_id = fluttf_id
        self.print_flag = print_flag
        self.velocity = np.asarray(velocity)
        self.rho = np.asarray(rho)
        self.atmos_ref = None
        self.mkaeroz_ref = None
        self.fluttf_ref = None
        assert isinstance(mkaeroz_id, integer_types), self.get_stats()

    @property
    def eas(self) -> np.ndarray:
        eas = 0.5 * self.rho * (self.velocity / self.vref) ** 2
        return eas

    @classmethod
    def add_card(cls, card: BDFCard, comment: str=''):
        """
        Adds a FIXMACH card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card

        """
        # FIXMACH SETID IDMK  ------  FTMUNIT FTLUNIT VREF FLUTTF PRINT
        #         V1    RHO1   V2     RHO2    etc.
        # FIXMACH 100   10            slug    ft      1.0         3
        #         900.  0.002         1000.  .00238 1200. .0030
        sid = integer(card, 1, 'sid')
        mkaeroz_id = integer(card, 2, 'mkaeroz_id')
        blank(card, 3, 'blank')
        mass_unit = string(card, 4, 'mass_unit')
        length_unit = string(card, 5, 'length_unit')
        vref = double_or_blank(card, 6, 'vref', default=1.0)
        fluttf_id = integer_or_blank(card, 7, 'fluttf_id', default=0)
        print_flag = integer(card, 8, 'print_flag')

        nfields_left = len(card) - 9
        assert nfields_left % 2 == 0, nfields_left
        assert nfields_left > 0, nfields_left
        rho = []
        velocity = []
        j = 1
        for ifield in range(9, len(card), 2):
            veli = double(card, ifield, f'velocity{j}')
            rhoi = double(card, ifield+1, f'rho{j}')
            rho.append(rhoi)
            velocity.append(veli)
            j += 1
        assert len(card) > 8, f'len(FIXEMATM card) = {len(card):d}\ncard={card}'
        assert isinstance(mkaeroz_id, integer_types)
        return FIXMACH(sid, mkaeroz_id, mass_unit, length_unit,
                       fluttf_id, print_flag, velocity, rho,
                       vref=vref, comment=comment)

    # def validate(self):
    #     assert self.true_g in ['TRUE', 'G'], 'true_g=%r' % self.true_g

    def cross_reference(self, model: BDF) -> None:
        assert isinstance(self.mkaeroz_id, integer_types), self.get_stats()
        if self.fluttf_id and 0:
            self.fluttf_ref = model.zaero.fluttf[self.fluttf_id]
        self.mkaeroz_ref = model.zaero.mkaeroz[self.mkaeroz_id]

    def safe_cross_reference(self, model: BDF, xref_errors):
        self.cross_reference(model)

    def uncross_reference(self) -> None:
        """Removes cross-reference links"""
        self.mkaeroz_ref = None
        self.fluttf_ref = None

    def plot(self, fig: plt.Figure):
        ntables = 1
        assert ntables > 0, ntables
        axes = fig.subplots(nrows=ntables)
        if ntables == 1:
            axes = [axes]
        assert len(axes) == ntables, (axes, ntables)

    def raw_fields(self):
        """
        Gets the fields in their unmodified form

        Returns
        -------
        fields : list[varies]
            the fields that define the card

        """
        list_fields = [
            'FIXMACH', self.sid, self.mkaeroz_id, '', self.mass_unit,
            self.length_unit, self.vref, self.fluttf_id, self.print_flag]
        for vel, rho in zip(self.velocity, self.rho):
            list_fields.extend([vel, rho])
        return list_fields

    def repr_fields(self):
        list_fields = self.raw_fields()
        return list_fields

    def write_card(self, size: int = 8, is_double: bool = False) -> str:
        card = self.repr_fields()
        return self.comment + print_card_8(card)


class FIXMDEN(BaseCard):
    type = 'FIXMDEN'
    # _field_map = {
    #     1: 'sid', 2: 'mach', 3: 'q', 8: 'aeqr',
    # }

    def __init__(self, sid: int, mkaeroz_id: int,
                 rho: float, mass_unit: str, length_unit: str,
                 fluttf_id: int, print_flag: int,
                 velocity: list[float],
                 vref: float=1.0, comment: str=''):
        BaseCard.__init__(self)
        if comment:
            self.comment = comment

        self.sid = sid
        self.mkaeroz_id = mkaeroz_id
        self.mass_unit = mass_unit
        self.length_unit = length_unit
        self.vref = vref
        self.fluttf_id = fluttf_id
        self.print_flag = print_flag
        self.velocity = np.asarray(velocity)
        self.rho = rho
        self.atmos_ref = None
        self.mkaeroz_ref = None
        assert isinstance(mkaeroz_id, integer_types), self.get_stats()

    @classmethod
    def add_card(cls, card: BDFCard, comment: str=''):
        """
        Adds a FIXMDEN card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card

        """
        # FIXMDEN SETID IDMK DEN    FTMUNIT FTLUNIT VREF FLUTTF PRINT CONT
        #         V1 V2
        # FIXMDEN 100   10   .00238 slug    ft      1.0         3
        #         900.  950. 1000.  1100.   1200.   2000.
        sid = integer(card, 1, 'sid')
        mkaeroz_id = integer(card, 2, 'mkaeroz_id')
        density = double(card, 3, 'density')
        mass_unit = string(card, 4, 'mass_unit')
        length_unit = string(card, 5, 'length_unit')
        vref = double_or_blank(card, 6, 'vref', default=1.0)
        fluttf_id = integer_or_blank(card, 7, 'fluttf_id', default=0)
        print_flag = integer_or_blank(card, 8, 'print_flag')

        velocity = []
        j = 1
        for ifield in range(9, len(card)):
            veli = double(card, ifield, f'velocity{j}')
            velocity.append(veli)
            j += 1
        assert len(card) > 8, f'len(FIXEMATM card) = {len(card):d}\ncard={card}'
        assert isinstance(mkaeroz_id, integer_types)
        return FIXMDEN(sid, mkaeroz_id, density, mass_unit, length_unit,
                       fluttf_id, print_flag, velocity,
                       vref=vref, comment=comment)

    # def validate(self):
    #     assert self.true_g in ['TRUE', 'G'], 'true_g=%r' % self.true_g

    def cross_reference(self, model: BDF) -> None:
        assert isinstance(self.mkaeroz_id, integer_types), self.get_stats()
        if self.fluttf_id and 0:
            self.fluttf_ref = model.zaero.fluttf[self.fluttf_id]
        self.mkaeroz_ref = model.zaero.mkaeroz[self.mkaeroz_id]

    def safe_cross_reference(self, model: BDF, xref_errors):
        self.cross_reference(model)

    def uncross_reference(self) -> None:
        """Removes cross-reference links"""
        self.mkaeroz_ref = None

    def plot(self, fig: plt.Figure):
        ntables = 1
        assert ntables > 0, ntables
        axes = fig.subplots(nrows=ntables)
        if ntables == 1:
            axes = [axes]
        assert len(axes) == ntables, (axes, ntables)

    def raw_fields(self):
        """
        Gets the fields in their unmodified form

        Returns
        -------
        fields : list[varies]
            the fields that define the card

        """
        list_fields = [
            'FIXMDEN', self.sid, self.mkaeroz_id, self.rho, self.mass_unit,
            self.length_unit, self.vref, self.fluttf_id, self.print_flag] + self.velocity.tolist()
        return list_fields

    def repr_fields(self):
        list_fields = self.raw_fields()
        return list_fields

    def write_card(self, size: int = 8, is_double: bool = False) -> str:
        card = self.repr_fields()
        return self.comment + print_card_8(card)


def cross_reference_atmos_mkaeroz_fluttf(card: FIXMATM, model):
    msg = f', which is required by {card.type}={card.sid}'
    assert isinstance(card.mkaeroz_id, integer_types), card.get_stats()
    if card.atm_id:
        card.atmos_ref = get_atmos(model.zaero, card.atm_id, msg=msg)
    if card.fluttf_id and 0:
        card.fluttf_ref = model.zaero.fluttf[card.fluttf_id]
    card.mkaeroz_ref = get_mkaeroz(model.zaero, card.mkaeroz_id, msg=msg)

