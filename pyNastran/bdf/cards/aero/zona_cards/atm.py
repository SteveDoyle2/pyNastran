# coding: utf-8
# pylint: disable=W0212,C0103
"""
All ZONA aero cards are defined in this file.  This includes:
 * ATM
 * FIXMATM

All cards are BaseCard objects.

"""
from __future__ import annotations
# from itertools import count
from typing import TYPE_CHECKING

from matplotlib import pyplot as plt

from pyNastran.utils.numpy_utils import integer_types
from pyNastran.bdf.field_writer_8 import print_card_8
from pyNastran.bdf.cards.base_card import BaseCard
from pyNastran.bdf.bdf_interface.assign_type import (
    integer, integer_or_blank, double_or_blank, string,
    blank, double,
    # integer_or_string,
    # string_or_blank,
    # integer_or_double,
    # integer_or_string, integer_string_or_blank,
    # string_multifield,  # parse_components as fcomponent
)

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
        self.mass_unit = mass_unit
        self.length_unit = length_unit
        self.temperature_unit = temperature_unit
        self.atmosphere_table = atmosphere_table

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
        temperature_unit = string(card, 3, 'temperature_unit')

        atmosphere_table = []
        j = 1
        for ifield in range(9, len(card), 4):
            alt = double(card, ifield, f'alt{j+1}')
            sos = double(card, ifield+1, f'sound{j+1}')
            rho = double(card, ifield+2, f'density{j+1}')
            temp = double(card, ifield+3, f'temperature{j+1}')
            atmosphere_table.extend([alt, sos, rho, temp])
            j += 1
        assert len(atmosphere_table) % 4 == 0
        assert len(atmosphere_table) // 4 > 0
        assert len(card) > 8, f'len(FIXEMATM card) = {len(card):d}\ncard={card}'
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
        list_fields = [
            'ATMOS', self.atmos_id, self.mass_unit, self.length_unit, self.temperature_unit,
            None, None, None, None,
        ] + self.atmosphere_table
        return list_fields

    def repr_fields(self):
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
                 mass_unit: str, length_unit: str, vref: float,
                 fluttf_id: int, print_flag: int, mkaeroz_ids: list[int],
                 comment: str=''):
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
        vref = double(card, 6, 'vref')
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
                       length_unit, vref, fluttf_id, print_flag, mkaeroz_ids, comment=comment)

    # def validate(self):
    #     assert self.true_g in ['TRUE', 'G'], 'true_g=%r' % self.true_g

    def cross_reference(self, model: BDF) -> None:
        if self.atm_id:
            self.atmos_ref = model.zona.atmos[self.atm_id]
        mkaerozs_ref = []
        for mkaeroz_id in self.mkaeroz_ids:
            mkaeroz_ref = model.zona.mkaeroz[mkaeroz_id]
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

    def __init__(self, sid: int, mkaeroz_id: int, atm_id: int,
                 mass_unit: str, length_unit: str, vref: float,
                 fluttf_id: int, print_flag: int, alts: list[float], comment: str=''):
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
        self.alts = alts
        self.atmos_ref = None
        self.mkaeroz_ref = None
        assert isinstance(mkaeroz_id, integer_types), self.get_stats()

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
        atm_id = integer(card, 3, 'atm_id')
        mass_unit = string(card, 4, 'mass_unit')
        length_unit = string(card, 5, 'length_unit')
        vref = double(card, 6, 'vref')
        fluttf_id = integer_or_blank(card, 7, 'fluttf_id', default=0)
        print_flag = integer(card, 8, 'print_flag')

        alts = []
        j = 1
        for ifield in range(9, len(card)):
            alt = double(card, ifield, f'alt{j}')
            alts.append(alt)
            j += 1
        assert len(card) > 8, f'len(FIXEMATM card) = {len(card):d}\ncard={card}'
        assert isinstance(mkaeroz_id, integer_types)
        return FIXMATM(sid, mkaeroz_id, atm_id, mass_unit,
                       length_unit, vref, fluttf_id, print_flag, alts, comment=comment)

    # def validate(self):
    #     assert self.true_g in ['TRUE', 'G'], 'true_g=%r' % self.true_g

    def cross_reference(self, model: BDF) -> None:
        if self.atm_id:
            self.atmos_ref = model.zona.atmos[self.atm_id]
        assert isinstance(self.mkaeroz_id, integer_types), self.get_stats()
        self.mkaeroz_ref = model.zona.mkaeroz[self.mkaeroz_id]

    def safe_cross_reference(self, model: BDF, xref_errors):
        self.cross_reference(model)

    def uncross_reference(self) -> None:
        """Removes cross-reference links"""
        self.atmos_ref = None
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
            self.length_unit, self.vref, self.fluttf_id, self.print_flag] + self.alts
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
        self.vref = vref
        self.fluttf_id = fluttf_id
        self.print_flag = print_flag
        self.velocity = velocity
        self.rho = rho
        assert isinstance(rho, list), rho
        self.atmos_ref = None
        self.mkaeroz_ref = None
        assert isinstance(mkaeroz_id, integer_types), self.get_stats()

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
        self.mkaeroz_ref = model.zona.mkaeroz[self.mkaeroz_id]

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
        self.velocity = velocity
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
        self.mkaeroz_ref = model.zona.mkaeroz[self.mkaeroz_id]

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
            self.length_unit, self.vref, self.fluttf_id, self.print_flag] + self.velocity
        return list_fields

    def repr_fields(self):
        list_fields = self.raw_fields()
        return list_fields

    def write_card(self, size: int = 8, is_double: bool = False) -> str:
        card = self.repr_fields()
        return self.comment + print_card_8(card)
