from __future__ import annotations
from itertools import zip_longest
from typing import TYPE_CHECKING
import numpy as np
#from pyNastran.bdf.field_writer_8 import print_card_8 # , print_float_8, print_field_8
#from pyNastran.bdf.field_writer_16 import print_card_16, print_scientific_16, print_field_16
#from pyNastran.bdf.field_writer_double import print_scientific_double
from pyNastran.bdf.bdf_interface.assign_type import (
    string, integer, double,
    integer_or_blank, double_or_blank, string_or_blank,
)
from pyNastran.bdf.bdf_interface.assign_type_force import force_double_or_blank
#from pyNastran.bdf.cards.materials import mat1_E_G_nu, get_G_default, set_blank_if_default

from pyNastran.dev.bdf_vectorized3.cards.base_card import Material, parse_material_check #get_print_card_8_16,
from pyNastran.dev.bdf_vectorized3.cards.write_utils import get_print_card, array_str, array_default_int, array_float, array_float_nan

if TYPE_CHECKING:  # pragma: no cover
    from pyNastran.bdf.bdf_interface.bdf_card import BDFCard
    #from pyNastran.dev.bdf_vectorized3.bdf import BDF
    from pyNastran.dev.bdf_vectorized3.types import TextIOLike


class MATT1(Material):
    """
    Specifies temperature-dependent material properties on MAT1 entry
    fields via TABLEMi entries.

    +-------+-------+-------+-------+-------+--------+------+------+-------+
    |   1   |   2   |   3   |   4   |   5   |    6   |  7   |  8   |   9   |
    +=======+=======+=======+=======+=======+========+======+======+=======+
    | MATT1 |  MID  |  T(E) |  T(G) | T(NU) | T(RHO) | T(A) |      | T(GE) |
    +-------+-------+-------+-------+-------+--------+------+------+-------+
    |       | T(ST) | T(SC) | T(SS) |       |        |      |      |       |
    +-------+-------+-------+-------+-------+--------+------+------+-------+

    """
    @Material.clear_check
    def clear(self) -> None:
        self.material_id = np.array([], dtype='int32')
        self.e_table = np.array([], dtype='int32')
        self.g_table = np.array([], dtype='int32')
        self.nu_table = np.array([], dtype='int32')
        self.rho_table = np.array([], dtype='int32')
        self.alpha_table = np.array([], dtype='int32')
        self.ge_table = np.array([], dtype='int32')
        self.st_table = np.array([], dtype='int32')
        self.sc_table = np.array([], dtype='int32')
        self.ss_table = np.array([], dtype='int32')

    def add(self, mid, e_table=None, g_table=None, nu_table=None, rho_table=None,
            a_table=None, ge_table=None, st_table=None, sc_table=None, ss_table=None,
            comment: str='') -> int:
        """Creates a MATT1 card"""
        self.cards.append((mid, e_table, g_table, nu_table, rho_table, a_table,
                           ge_table, st_table, sc_table, ss_table, comment))
        self.n += 1
        return self.n

    def add_card(self, card: BDFCard, comment: str='') -> int:
        """
        Adds a MATT1 card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card

        """
        mid = integer(card, 1, 'mid')
        e_table = integer_or_blank(card, 2, 'T(E)')
        g_table = integer_or_blank(card, 3, 'T(G)')
        nu_table = integer_or_blank(card, 4, 'T(nu)')
        rho_table = integer_or_blank(card, 5, 'T(rho)')
        a_table = integer_or_blank(card, 6, 'T(A)')
        ge_table = integer_or_blank(card, 8, 'T(ge)')
        st_table = integer_or_blank(card, 9, 'T(st)')
        sc_table = integer_or_blank(card, 10, 'T(sc)')
        ss_table = integer_or_blank(card, 11, 'T(ss)')

        assert len(card) <= 12, f'len(MATT1 card) = {len(card):d}\ncard={card}'
        #return MATT1(mid, e_table, g_table, nu_table, rho_table, a_table,
                     #ge_table, st_table, sc_table, ss_table, comment=comment)
        self.cards.append((mid, e_table, g_table, nu_table, rho_table, a_table,
                           ge_table, st_table, sc_table, ss_table, comment))
        self.n += 1
        return self.n

    @Material.parse_cards_check
    def parse_cards(self) -> None:
        ncards = len(self.cards)
        material_id = np.zeros(ncards, dtype='int32')
        e_table = np.zeros(ncards, dtype='int32')
        g_table = np.zeros(ncards, dtype='int32')
        nu_table = np.zeros(ncards, dtype='int32')
        rho_table = np.zeros(ncards, dtype='int32')
        alpha_table = np.zeros(ncards, dtype='int32')
        ge_table = np.zeros(ncards, dtype='int32')
        st_table = np.zeros(ncards, dtype='int32')
        sc_table = np.zeros(ncards, dtype='int32')
        ss_table = np.zeros(ncards, dtype='int32')

        for i, card in enumerate(self.cards):
            (mid, e_tablei, g_tablei, nu_tablei, rho_tablei, a_tablei,
             ge_tablei, st_tablei, sc_tablei, ss_tablei, comment) = card
            e_tablei = 0 if e_tablei is None else e_tablei
            g_tablei = 0 if g_tablei is None else g_tablei
            nu_tablei = 0 if nu_tablei is None else nu_tablei
            rho_tablei = 0 if rho_tablei is None else rho_tablei
            a_tablei = 0 if a_tablei is None else a_tablei
            ge_tablei = 0 if ge_tablei is None else ge_tablei
            st_tablei = 0 if st_tablei is None else st_tablei
            sc_tablei = 0 if sc_tablei is None else sc_tablei
            ss_tablei = 0 if ss_tablei is None else ss_tablei
            material_id[i] = mid
            e_table[i] = e_tablei
            g_table[i] = g_tablei
            nu_table[i] = nu_tablei
            rho_table[i] = rho_tablei
            alpha_table[i] = a_tablei
            ge_table[i] = ge_tablei
            st_table[i] = st_tablei
            sc_table[i] = sc_tablei
            ss_table[i] = ss_tablei
        self._save(material_id, e_table, g_table, nu_table,
                   rho_table, alpha_table,
                   ge_table,
                   st_table, sc_table, ss_table,)
        self.sort()
        self.cards = []

    def _save(self, material_id, e_table, g_table, nu_table,
              rho_table, alpha_table,
              ge_table,
              st_table, sc_table, ss_table,):
        if len(self.material_id):
            asdf
            #material_id = np.hstack([self.material_id, material_id])
            #E = np.hstack([self.E, E])
            #G = np.hstack([self.G, G])
            #nu = np.hstack([self.nu, nu])
            #rho = np.hstack([self.rho, rho])
            #alpha = np.hstack([self.alpha, alpha])
            #tref = np.hstack([self.tref, tref])
            #ge = np.hstack([self.ge, ge])
            #Ss = np.hstack([self.Ss, Ss])
            #St = np.hstack([self.St, St])
            #Sc = np.hstack([self.Sc, Sc])
        self.material_id = material_id
        self.e_table = e_table
        self.g_table = g_table
        self.nu_table = nu_table
        self.rho_table = rho_table
        self.alpha_table = alpha_table
        self.ge_table = ge_table
        self.st_table = st_table
        self.sc_table = sc_table
        self.ss_table = ss_table

    def set_used(self, used_dict: dict[str, list[np.ndarray]]) -> None:
        table0 = np.hstack([
            self.e_table, self.g_table, self.nu_table,
            self.rho_table, self.alpha_table,
            self.ge_table, self.st_table, self.ss_table, ])
        utable0 = np.unique(table0)
        table = np.setdiff1d(utable0, [0])
        used_dict['tablem_id'].append(table)

    def convert(self, stiffness_scale: float=1.0,
                density_scale: float=1.0,
                alpha_scale: float=1.0,
                temperature_scale: float=1.0,
                stress_scale: float=1.0, **kwargs) -> None:
        pass
        tables0 = [
            (self.e_table, stiffness_scale),
            (self.g_table, stiffness_scale),
            (self.rho_table, density_scale),
            # nu - nope
            (self.alpha_table, alpha_scale),
            (self.st_table, stress_scale),
            (self.sc_table, stress_scale),
            (self.ss_table, stress_scale),
        ]
        table_ids = {}
        for table0, scale in tables0:
            utable0 = np.unique(table0)
            table = np.setdiff1d(utable0, [0])
            #if table:
                #table_ids[scale].append(table)
            #model.tablem

    def __apply_slice__(self, mat: MATT1, i: np.ndarray) -> None:  # ignore[override]
        mat.n = len(i)
        mat.material_id = self.material_id[i]
        mat.e_table = self.e_table[i]
        mat.g_table = self.g_table[i]
        mat.nu_table = self.nu_table[i]
        mat.rho_table = self.rho_table[i]
        mat.alpha_table = self.alpha_table[i]
        mat.ge_table = self.ge_table[i]
        mat.st_table = self.st_table[i]
        mat.sc_table = self.sc_table[i]
        mat.ss_table = self.ss_table[i]

    def geom_check(self, missing: dict[str, np.ndarray]):
        pass
    @property
    def max_id(self):
        tables = np.hstack([self.material_id,
                            self.e_table, self.g_table, self.nu_table,
                            self.rho_table,
                            self.alpha_table,
                            self.st_table, self.sc_table, self.ss_table,
                            self.ge_table, ])
        return tables.max()

    @parse_material_check
    def write_file(self, bdf_file: TextIOLike,
                   size: int=8, is_double: bool=False,
                   write_card_header: bool=False) -> None:
        print_card = get_print_card(size, self.max_id)

        material_id = array_str(self.material_id, size=size)
        e_table = array_default_int(self.e_table, default=0, size=size)
        g_table = array_default_int(self.g_table, default=0, size=size)
        nu_table = array_default_int(self.nu_table, default=0, size=size)
        rho_table = array_default_int(self.rho_table, default=0, size=size)
        alpha_table = array_default_int(self.alpha_table, default=0, size=size)
        ge_table = array_default_int(self.ge_table, default=0, size=size)
        ss_table = array_default_int(self.ss_table, default=0, size=size)
        st_table = array_default_int(self.st_table, default=0, size=size)
        sc_table = array_default_int(self.sc_table, default=0, size=size)

        #tables = np.column_stack([self.e_table, self.g_table, self.nu_table, self.rho_table])
        for mid, e, g, nu, rho, alpha, ge, ss, st, sc, \
            in zip_longest(material_id, e_table, g_table, nu_table, rho_table,
                           alpha_table, ge_table,
                           ss_table, st_table, sc_table):
            list_fields = ['MATT1', mid, e, g, nu, rho, alpha, '', ge,
                           st, sc, ss]
            bdf_file.write(print_card(list_fields))
        return


class MATS1(Material):
    """
    Specifies stress-dependent material properties for use in applications
    involving nonlinear materials. This entry is used if a MAT1, MAT2 or MAT9
    entry is specified with the same MID in a nonlinear solution sequence
    (SOLs 106 and 129).

    #: Identification number of a MAT1, MAT2, or MAT9 entry.
    self.mid = mid

    #: Identification number of a TABLES1 or TABLEST entry. If H is
    #: given, then this field must be blank.
    self.tid = tid

    #: Type of material nonlinearity. ('NLELAST' for nonlinear elastic
    #: or 'PLASTIC' for elastoplastic.)
    self.Type = Type

    #: Work hardening slope (slope of stress versus plastic strain)
    #: in units of stress. For elastic-perfectly plastic cases,
    #: H=0.0.  For more than a single slope in the plastic range,
    #: the stress-strain data must be supplied on a TABLES1 entry
    #: referenced by TID, and this field must be blank
    self.h = h

    #: Hardening Rule, selected by one of the following values
    #: (Integer): (1) Isotropic (Default) (2) Kinematic
    #: (3) Combined isotropic and kinematic hardening
    self.hr = hr

    #: Yield function criterion, selected by one of the following
    #: values (1) Von Mises (2) Tresca (3) Mohr-Coulomb
    #: (4) Drucker-Prager
    self.yf = yf
    #: Initial yield point
    self.limit1 = limit1
    #: Internal friction angle, measured in degrees, for the
    #: Mohr-Coulomb and Drucker-Prager yield criteria
    self.limit2 = limit2
        if self.Type not in ['NLELAST', 'PLASTIC', 'PLSTRN']:
            raise ValueError('MATS1 Type must be [NLELAST, PLASTIC, PLSTRN]; Type=%r' % self.Type)

    """
    @Material.clear_check
    def clear(self) -> None:
        self.material_id = np.array([], dtype='int32')
        self.table_id = np.array([], dtype='int32')
        self.Type = np.array([], dtype='|U8')
        self.hardening_slope = np.array([], dtype='float64')
        self.hr = np.array([], dtype='int32')
        self.yf = np.array([], dtype='int32')
        self.limit1 = np.array([], dtype='float64')
        self.limit2 = np.array([], dtype='float64')
        self.stress_strain_measure = np.array([], dtype='|U8')

    def add(self, mid: int, tid: int, Type: str,
            h: float, hr: int, yf: int, limit1: float, limit2: float,
            stress_strain_measure: str='', comment: str='') -> int:
        """Creates a MATS1 card"""
        assert isinstance(Type, str), f'Type={Type!r}'
        self.cards.append((mid, tid, Type, h, hr, yf, limit1, limit2, stress_strain_measure, comment))
        self.n += 1
        return self.n

    def add_card(self, card: BDFCard, comment: str='') -> int:
        """
        Adds a MATS1 card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card

        """
        mid = integer(card, 1, 'mid')
        tid = integer_or_blank(card, 2, 'tables1_id')
        Type = string(card, 3, 'Type')

        if Type not in {'NLELAST', 'PLASTIC', 'PLSTRN'}:
            raise ValueError('MATS1 Type must be [NLELAST, PLASTIC, PLSTRN]; Type=%r' % Type)
        if Type == 'NLELAST':
            # should we even read these?
            hardening_slope = None
            hr = None
            yf = None
            limit1 = None
            limit2 = None
            #h = blank(card, 4, 'h')
            #hr = blank(card, 6, 'hr')
            #yf = blank(card, 5, 'yf')
            #limit1 = blank(card, 7, 'yf')
            #limit2 = blank(card, 8, 'yf')
        else:
            fdouble_or_blank = force_double_or_blank if self.model.is_lax_parser else double_or_blank
            hardening_slope = fdouble_or_blank(card, 4, 'H', default=0.0)
            yf = integer_or_blank(card, 5, 'yf', default=1)
            hr = integer_or_blank(card, 6, 'hr', default=1)
            limit1 = double(card, 7, 'limit1')

            if yf in [3, 4]:
                limit2 = double(card, 8, 'limit2')
            else:
                #limit2 = blank(card, 8, 'limit2')
                limit2 = None
        stress_strain_measure = string_or_blank(card, 10, 'stress/strain measure', default='')
        assert len(card) <= 9, f'len(MATS1 card) = {len(card):d}\ncard={card}'
        #return MATS1(mid, tid, Type, h, hr, yf, limit1, limit2, comment=comment)
        self.cards.append((mid, tid, Type, hardening_slope, hr, yf,
                           limit1, limit2, stress_strain_measure, comment))
        self.n += 1
        return self.n

    @Material.parse_cards_check
    def parse_cards(self) -> None:
        ncards = len(self.cards)
        material_id = np.zeros(ncards, dtype='int32')
        table_id = np.zeros(ncards, dtype='int32')
        Type = np.zeros(ncards, dtype='|U8')
        hardening_slope = np.zeros(ncards, dtype='float64')
        hr = np.zeros(ncards, dtype='int32')
        yf = np.zeros(ncards, dtype='int32')
        limit1 = np.zeros(ncards, dtype='float64')
        limit2 = np.zeros(ncards, dtype='float64')
        stress_strain_measure = np.zeros(ncards, dtype='|U8')

        for i, card in enumerate(self.cards):
            (mid, tid, typei, hardening_slopei, hri, yfi,
             limit1i, limit2i, stress_strain_measurei, comment) = card
            tid = 0 if tid is None else tid
            hri = 1 if hri is None else hri
            yfi = 1 if yfi is None else yfi
            material_id[i] = mid
            table_id[i] = tid
            Type[i] = typei
            hardening_slope[i] = hardening_slopei
            hr[i] = hri
            yf[i] = yfi
            limit1[i] = limit1i
            limit2[i] = limit2i
            stress_strain_measure[i] = stress_strain_measurei
        self._save(material_id, table_id, Type, hardening_slope, hr, yf,
                   limit1, limit2, stress_strain_measure)
        self.sort()
        self.cards = []

    def _save(self, material_id, table_id, Type, hardening_slope, hr, yf,
              limit1, limit2, stress_strain_measure):
        if len(self.material_id):
            asdf
            #material_id = np.hstack([self.material_id, material_id])
            #E = np.hstack([self.E, E])
        self.material_id = material_id
        self.table_id = table_id
        self.Type = Type
        self.hardening_slope = hardening_slope
        self.hr = hr
        self.yf = yf
        self.limit1 = limit1
        self.limit2 = limit2
        self.stress_strain_measure = stress_strain_measure

    def set_used(self, used_dict: dict[str, list[np.ndarray]]) -> None:
        pass
        #table0 = np.hstack([
            #self.e_table, self.g_table, self.nu_table,
            #self.rho_table, self.alpha_table,
            #self.ge_table, self.st_table, self.ss_table, ])
        #utable0 = np.unique(table0)
        #table = np.setdiff1d(utable0, [0])
        #used_dict['tablem_id'].append(table)

    #def convert(self, stiffness_scale: float=1.0,
                #density_scale: float=1.0,
                #alpha_scale: float=1.0,
                #temperature_scale: float=1.0,
                #stress_scale: float=1.0, **kwargs) -> None:
        #tables0 = [
            #(self.e_table, stiffness_scale),
            #(self.g_table, stiffness_scale),
            #(self.rho_table, density_scale),
            ## nu - nope
            #(self.alpha_table, alpha_scale),
            #(self.st_table, stress_scale),
            #(self.sc_table, stress_scale),
            #(self.ss_table, stress_scale),
        #]
        #table_ids = {}
        #for table0, scale in tables0:
            #utable0 = np.unique(table0)
            #table = np.setdiff1d(utable0, [0])
            #if table:
                #table_ids[scale].append(table)
            #model.tablem

    def __apply_slice__(self, mat: MATT1, i: np.ndarray) -> None:  # ignore[override]
        mat.n = len(i)
        mat.material_id = self.material_id[i]
        mat.table_id = self.table_id[i]
        mat.Type = self.Type[i]
        mat.hardening_slope = self.hardening_slope[i]
        mat.hr = self.hr[i]
        mat.yf = self.yf[i]
        mat.limit1 = self.limit1[i]
        mat.limit2 = self.limit2[i]
        mat.stress_strain_measure = self.stress_strain_measure[i]

    def geom_check(self, missing: dict[str, np.ndarray]):
        pass

    @property
    def max_id(self):
        return max(self.material_id.max(), self.table_id.max())

    @parse_material_check
    def write_file(self, bdf_file: TextIOLike,
                   size: int=8, is_double: bool=False,
                   write_card_header: bool=False) -> None:
        print_card = get_print_card(size, self.max_id)

        material_id = array_str(self.material_id, size=size)
        table_id = array_str(self.table_id, size=size)
        hardening_slope = array_float(self.hardening_slope, size=size)
        yf = array_str(self.yf, size=size)
        hr = array_str(self.hr, size=size)
        limit1 = array_float(self.limit1, size=size)
        limit2 = array_float_nan(self.limit2, size=size)

        #tables = np.column_stack([self.e_table, self.g_table, self.nu_table, self.rho_table])
        for mid, table_idi, typei, hi, yfi, hri, limit1i, limit2i, stress_strain_measure \
            in zip_longest(material_id, table_id, self.Type,
                           hardening_slope, yf, hr, limit1, limit2, self.stress_strain_measure):
            list_fields = ['MATS1', mid, table_idi, typei,
                           hi, yfi, hri, limit1i, limit2i, stress_strain_measure]
            bdf_file.write(print_card(list_fields))
        return

#MATS8 #  MSC card

class MATT2(Material):
    """
    Specifies temperature-dependent material properties on MAT2 entry
    fields via TABLEMi entries.

    +-------+-------+--------+--------+--------+--------+--------+--------+--------+
    |   1   |   2   |    3   |   4    |   5    |    6   |    7   |    8   |   9    |
    +=======+=======+========+========+========+========+========+========+========+
    | MATT2 |  MID  | T(G12) | T(G13) | T(G13) | T(G22) | T(G23) | T(G33) | T(RHO) |
    +-------+-------+--------+--------+--------+--------+--------+--------+--------+
    |       | T(A1) | T(A2)  | T(A3)  |        | T(GE)  | T(ST)  | T(SC)  |  T(SS) |
    +-------+-------+--------+--------+--------+--------+--------+--------+--------+

    """
    @Material.clear_check
    def clear(self) -> None:
        self.material_id = np.array([], dtype='int32')
        self.g11_table = np.array([], dtype='int32')
        self.g22_table = np.array([], dtype='int32')
        self.g33_table = np.array([], dtype='int32')
        self.g12_table = np.array([], dtype='int32')
        self.g13_table = np.array([], dtype='int32')
        self.g23_table = np.array([], dtype='int32')
        self.rho_table = np.array([], dtype='int32')
        self.a1_table = np.array([], dtype='int32')
        self.a2_table = np.array([], dtype='int32')
        self.a3_table = np.array([], dtype='int32')
        self.st_table = np.array([], dtype='int32')
        self.sc_table = np.array([], dtype='int32')
        self.ss_table = np.array([], dtype='int32')
        self.ge_table = np.array([], dtype='int32')

    def add(self, mid: int,
            g11_table=None, g22_table=None, g33_table=None,
            g12_table=None, g13_table=None, g23_table=None,
            rho_table=None,
            a1_table=None, a2_table=None, a3_table=None,
            st_table=None, sc_table=None, ss_table=None,
            ge_table=None, comment: str='') -> int:
        """Creates a MATT2 card"""
        self.cards.append((mid,
                           g11_table, g22_table, g33_table,
                           g12_table, g13_table, g23_table,
                           rho_table,
                           a1_table, a2_table, a3_table,
                           st_table, sc_table, ss_table,
                           ge_table, comment))
        self.n += 1
        return self.n

    def add_card(self, card: BDFCard, comment: str='') -> int:
        """
        Adds a MATT2 card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card

        """
        mid = integer(card, 1, 'mid')
        g11_table = integer_or_blank(card, 2, 'T(G11)', default=0)
        g12_table = integer_or_blank(card, 3, 'T(G12)', default=0)
        g13_table = integer_or_blank(card, 4, 'T(G13)', default=0)
        g22_table = integer_or_blank(card, 5, 'T(G22)', default=0)
        g23_table = integer_or_blank(card, 6, 'T(G23)', default=0)
        g33_table = integer_or_blank(card, 7, 'T(G33)', default=0)
        rho_table = integer_or_blank(card, 8, 'T(rho)', default=0)
        a1_table = integer_or_blank(card, 9, 'T(A1)', default=0)
        a2_table = integer_or_blank(card, 10, 'T(A2)', default=0)
        a3_table = integer_or_blank(card, 11, 'T(A3)', default=0)
        ge_table = integer_or_blank(card, 13, 'T(ge)', default=0)
        st_table = integer_or_blank(card, 14, 'T(st)', default=0)
        sc_table = integer_or_blank(card, 15, 'T(sc)', default=0)
        ss_table = integer_or_blank(card, 16, 'T(ss)', default=0)

        assert len(card) <= 17, f'len(MATT2 card) = {len(card):d}\ncard={card}'
        #return MATT2(mid, g11_table, g12_table, g13_table, g22_table, g23_table,
                     #g33_table, rho_table, a1_table,
                     #a2_table, a3_table, ge_table,
                     #st_table, sc_table, ss_table,
                     #comment=comment)
        self.cards.append((mid,
                           g11_table, g22_table, g33_table,
                           g12_table, g13_table, g23_table,
                           rho_table,
                           a1_table, a2_table, a3_table,
                           st_table, sc_table, ss_table,
                           ge_table, comment))
        self.n += 1
        return self.n

    @Material.parse_cards_check
    def parse_cards(self) -> None:
        ncards = len(self.cards)
        idtype = self.model.idtype
        material_id = np.zeros(ncards, dtype=idtype)

        g11_table = np.zeros(ncards, dtype=idtype)
        g22_table = np.zeros(ncards, dtype=idtype)
        g33_table = np.zeros(ncards, dtype=idtype)
        g12_table = np.zeros(ncards, dtype=idtype)
        g13_table = np.zeros(ncards, dtype=idtype)
        g23_table = np.zeros(ncards, dtype=idtype)
        rho_table = np.zeros(ncards, dtype=idtype)
        a1_table = np.zeros(ncards, dtype=idtype)
        a2_table = np.zeros(ncards, dtype=idtype)
        a3_table = np.zeros(ncards, dtype=idtype)
        st_table = np.zeros(ncards, dtype=idtype)
        sc_table = np.zeros(ncards, dtype=idtype)
        ss_table = np.zeros(ncards, dtype=idtype)
        ge_table = np.zeros(ncards, dtype=idtype)
        for i, card in enumerate(self.cards):
            (mid,
             g11_tablei, g22_tablei, g33_tablei,
             g12_tablei, g13_tablei, g23_tablei,
             rho_tablei,
             a1_tablei, a2_tablei, a3_tablei,
             st_tablei, sc_tablei, ss_tablei,
             ge_tablei, comment) = card

            g11_tablei = 0 if g11_tablei is None else g11_tablei
            g22_tablei = 0 if g22_tablei is None else g22_tablei
            g33_tablei = 0 if g33_tablei is None else g33_tablei

            g12_tablei = 0 if g12_tablei is None else g12_tablei
            g13_tablei = 0 if g13_tablei is None else g13_tablei
            g23_tablei = 0 if g23_tablei is None else g23_tablei
            rho_tablei = 0 if rho_tablei is None else rho_tablei

            a1_tablei = 0 if a1_tablei is None else a1_tablei
            a2_tablei = 0 if a2_tablei is None else a2_tablei
            a3_tablei = 0 if a3_tablei is None else a3_tablei

            st_tablei = 0 if st_tablei is None else st_tablei
            sc_tablei = 0 if sc_tablei is None else sc_tablei
            ss_tablei = 0 if ss_tablei is None else ss_tablei

            ge_tablei = 0 if ge_tablei is None else ge_tablei
            material_id[i] = mid
            g11_table[i] = g11_tablei
            g22_table[i] = g12_tablei
            g33_table[i] = g13_tablei

            g12_table[i] = g12_tablei
            g13_table[i] = g13_tablei
            g23_table[i] = g23_tablei

            rho_table[i] = rho_tablei
            a1_table[i] = a1_tablei
            a2_table[i] = a2_tablei
            a3_table[i] = a3_tablei

            st_table[i] = st_tablei
            sc_table[i] = sc_tablei
            ss_table[i] = ss_tablei
            ge_table[i] = ge_tablei
        self._save(material_id,
                   g11_table, g22_table, g33_table,
                   g12_table, g13_table, g23_table,
                   rho_table,
                   a1_table, a2_table, a3_table,
                   st_table, sc_table, ss_table,
                   ge_table)
        self.sort()
        self.cards = []

    def _save(self, material_id,
              g11_table, g22_table, g33_table,
              g12_table, g13_table, g23_table,
              rho_table,
              a1_table, a2_table, a3_table,
              st_table, sc_table, ss_table,
              ge_table):
        if len(self.material_id):
            asdf
            #material_id = np.hstack([self.material_id, material_id])
            #E = np.hstack([self.E, E])

        self.material_id = material_id
        self.g11_table = g11_table
        self.g22_table = g22_table
        self.g33_table = g33_table
        self.g12_table = g12_table
        self.g13_table = g13_table
        self.g23_table = g23_table
        self.rho_table = rho_table
        self.a1_table = a1_table
        self.a2_table = a2_table
        self.a3_table = a3_table
        self.st_table = st_table
        self.sc_table = sc_table
        self.ss_table = ss_table
        self.ge_table = ge_table

    def set_used(self, used_dict: dict[str, list[np.ndarray]]) -> None:
        table0 = np.hstack([
            self.g11_table, self.g22_table, self.g33_table,
            self.g12_table, self.g13_table, self.g23_table,
            self.rho_table,
            self.a1_table, self.a2_table, self.a3_table,
            self.ge_table, self.st_table, self.ss_table, ])
        utable0 = np.unique(table0)
        table = np.setdiff1d(utable0, [0])
        used_dict['tablem_id'].append(table)

    #def convert(self, stiffness_scale: float=1.0,
                #density_scale: float=1.0,
                #alpha_scale: float=1.0,
                #temperature_scale: float=1.0,
                #stress_scale: float=1.0, **kwargs) -> None:
        #tables0 = [
            #(self.e_table, stiffness_scale),
            #(self.g_table, stiffness_scale),
            #(self.rho_table, density_scale),
            ## nu - nope
            #(self.alpha_table, alpha_scale),
            #(self.st_table, stress_scale),
            #(self.sc_table, stress_scale),
            #(self.ss_table, stress_scale),
        #]
        #table_ids = {}
        #for table0, scale in tables0:
            #utable0 = np.unique(table0)
            #table = np.setdiff1d(utable0, [0])
            #if table:
                #table_ids[scale].append(table)
            #model.tablem

    def __apply_slice__(self, mat: MATT1, i: np.ndarray) -> None:  # ignore[override]
        mat.n = len(i)
        mat.material_id = self.material_id[i]
        mat.g11_table = self.g11_table[i]
        mat.g22_table = self.g22_table[i]
        mat.g33_table = self.g33_table[i]
        mat.g12_table = self.g12_table[i]
        mat.g13_table = self.g13_table[i]
        mat.g23_table = self.g23_table[i]
        mat.rho_table = self.rho_table[i]

        mat.a1_table = self.a1_table[i]
        mat.a2_table = self.a2_table[i]
        mat.a3_table = self.a3_table[i]

        mat.st_table = self.st_table[i]
        mat.sc_table = self.sc_table[i]
        mat.ss_table = self.ss_table[i]
        mat.ge_table = self.ge_table[i]

    def geom_check(self, missing: dict[str, np.ndarray]):
        pass

    @property
    def max_id(self):
        tables = np.hstack([self.material_id,
                            self.g11_table, self.g22_table, self.g33_table,
                            self.g12_table, self.g13_table, self.g23_table,
                            self.rho_table,
                            self.a1_table, self.a2_table, self.a3_table,
                            self.st_table, self.sc_table, self.ss_table,
                            self.ge_table, ])
        return tables.max()

    @parse_material_check
    def write_file(self, bdf_file: TextIOLike,
                   size: int=8, is_double: bool=False,
                   write_card_header: bool=False) -> None:
        print_card = get_print_card(size, self.max_id)

        material_id = array_str(self.material_id, size=size)
        g11_table = array_default_int(self.g11_table, default=0, size=size)
        g22_table = array_default_int(self.g22_table, default=0, size=size)
        g33_table = array_default_int(self.g33_table, default=0, size=size)
        g12_table = array_default_int(self.g12_table, default=0, size=size)
        g13_table = array_default_int(self.g13_table, default=0, size=size)
        g23_table = array_default_int(self.g23_table, default=0, size=size)
        rho_table = array_default_int(self.rho_table, default=0, size=size)
        a1_table = array_default_int(self.a1_table, default=0, size=size)
        a2_table = array_default_int(self.a2_table, default=0, size=size)
        a3_table = array_default_int(self.a3_table, default=0, size=size)
        st_table = array_default_int(self.st_table, default=0, size=size)
        sc_table = array_default_int(self.sc_table, default=0, size=size)
        ss_table = array_default_int(self.ss_table, default=0, size=size)
        ge_table = array_default_int(self.ge_table, default=0, size=size)

        for (mid,
             g11_tablei, g22_tablei, g33_tablei,
             g12_tablei, g13_tablei, g23_tablei,
             rho_tablei,
             a1_tablei, a2_tablei, a3_tablei,
             st_tablei, sc_tablei, ss_tablei, ge_tablei) in zip_longest(
                 material_id,
                 g11_table, g22_table, g33_table,
                 g12_table, g13_table, g23_table, rho_table,
                 a1_table, a2_table, a3_table,
                 st_table, sc_table, ss_table,
                 ge_table):

            list_fields = [
                'MATT2', mid, g11_tablei, g12_tablei,
                g13_tablei, g22_tablei, g23_tablei,
                g33_tablei, rho_tablei, a1_tablei,
                a2_tablei, a3_tablei, None, ge_tablei,
                st_tablei, sc_tablei, ss_tablei]
            bdf_file.write(print_card(list_fields))
        return


class MATT8(Material):
    """
    Specifies temperature-dependent material properties on MAT2 entry
    fields via TABLEMi entries.

    +-------+--------+--------+-------+---------+--------+--------+--------+--------+
    |   1   |   2    |   3    |   4   |    5    |   6    |   7    |    8   |   9    |
    +=======+========+========+=======+=========+========+========+========+========+
    | MATT8 |  MID   | T(E1)  | T(E2) | T(Nu12) | T(G12) | T(G1z) | T(G2z) | T(RHO) |
    +-------+--------+--------+-------+---------+--------+--------+--------+--------+
    |       |  T(A1) | T(A2)  |       |  T(Xt)  | T(Xc)  | T(Yt)  | T(Yc)  | T(S)   |
    +-------+--------+--------+-------+---------+--------+--------+--------+--------+
    |       |  T(GE) | T(F12) |       |         |        |        |        |        |
    +-------+--------+--------+-------+---------+--------+--------+--------+--------+

    """
    @Material.clear_check
    def clear(self) -> None:
        self.material_id = np.array([], dtype='int32')

    def add(self, mid: int, e1_table=None, e2_table=None, nu12_table=None,
            g12_table=None, g1z_table=None, g2z_table=None, rho_table=None,
            a1_table=None, a2_table=None,
            xt_table=None, xc_table=None, yt_table=None, yc_table=None,
            s_table=None, ge_table=None, f12_table=None, comment: str='') -> int:
        """Creates a MATT8 card"""
        self.cards.append((mid, e1_table, e2_table, nu12_table, g12_table,
                           g1z_table, g2z_table, rho_table,
                           a1_table, a2_table, xt_table,
                           xc_table, yt_table, yc_table,
                           s_table, ge_table, f12_table, comment))
        self.n += 1
        return self.n

    def add_card(self, card: BDFCard, comment: str='') -> int:
        """
        Adds a MATT8 card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card

        """
        mid = integer(card, 1, 'mid')
        e1_table = integer_or_blank(card, 2, 'T(E1)')
        e2_table = integer_or_blank(card, 3, 'T(E2)')
        nu12_table = integer_or_blank(card, 4, 'T(Nu12)')
        g12_table = integer_or_blank(card, 5, 'T(G12)')
        g1z_table = integer_or_blank(card, 6, 'T(G1z)')
        g2z_table = integer_or_blank(card, 7, 'T(G2z)')
        rho_table = integer_or_blank(card, 8, 'T(Rho)')
        a1_table = integer_or_blank(card, 9, 'T(A1)')
        a2_table = integer_or_blank(card, 10, 'T(A2)')

        xt_table = integer_or_blank(card, 12, 'T(Xt)')
        xc_table = integer_or_blank(card, 13, 'T(Xc)')
        yt_table = integer_or_blank(card, 14, 'T(Yt)')
        yc_table = integer_or_blank(card, 15, 'T(Yc)')
        s_table = integer_or_blank(card, 16, 'T(S)')
        ge_table = integer_or_blank(card, 17, 'T(GE)')
        f12_table = integer_or_blank(card, 18, 'T(F12)')

        assert len(card) <= 19, f'len(MATT8 card) = {len(card):d}\ncard={card}'
        #return MATT8(mid, e1_table, e2_table, nu12_table, g12_table,
                     #g1z_table, g2z_table, rho_table,
                     #a1_table, a2_table, xt_table,
                     #xc_table, yt_table, yc_table,
                     #s_table, ge_table, f12_table,
                     #comment=comment)
        self.cards.append((mid, e1_table, e2_table, nu12_table, g12_table,
                           g1z_table, g2z_table, rho_table,
                           a1_table, a2_table, xt_table,
                           xc_table, yt_table, yc_table,
                           s_table, ge_table, f12_table, comment))
        self.n += 1
        return self.n

    @Material.parse_cards_check
    def parse_cards(self) -> None:
        ncards = len(self.cards)
        idtype = self.model.idtype
        material_id = np.zeros(ncards, dtype=idtype)
        e1_table = np.zeros(ncards, dtype=idtype)
        e2_table = np.zeros(ncards, dtype=idtype)
        nu12_table = np.zeros(ncards, dtype=idtype)
        g12_table = np.zeros(ncards, dtype=idtype)
        g1z_table = np.zeros(ncards, dtype=idtype)
        g2z_table = np.zeros(ncards, dtype=idtype)
        rho_table = np.zeros(ncards, dtype=idtype)
        a1_table = np.zeros(ncards, dtype=idtype)
        a2_table = np.zeros(ncards, dtype=idtype)
        xt_table = np.zeros(ncards, dtype=idtype)
        xc_table = np.zeros(ncards, dtype=idtype)
        yt_table = np.zeros(ncards, dtype=idtype)
        yc_table = np.zeros(ncards, dtype=idtype)
        s_table = np.zeros(ncards, dtype=idtype)
        ge_table = np.zeros(ncards, dtype=idtype)
        f12_table = np.zeros(ncards, dtype=idtype)
        for i, card in enumerate(self.cards):
            (mid, e1_tablei, e2_tablei, nu12_tablei, g12_tablei,
             g1z_tablei, g2z_tablei, rho_tablei,
             a1_tablei, a2_tablei, xt_tablei,
             xc_tablei, yt_tablei, yc_tablei,
             s_tablei, ge_tablei, f12_tablei, comment) = card
            e1_tablei = 0 if e1_tablei  is None else e1_tablei
            e2_tablei = 0 if e2_tablei  is None else e2_tablei
            nu12_tablei = 0 if nu12_tablei is None else nu12_tablei
            g12_tablei = 0 if g12_tablei is None else g12_tablei
            g1z_tablei = 0 if g1z_tablei is None else g1z_tablei
            g2z_tablei = 0 if g2z_tablei is None else g2z_tablei
            rho_tablei = 0 if rho_tablei is None else rho_tablei
            a1_tablei = 0 if a1_tablei is None else a1_tablei
            a2_tablei = 0 if a2_tablei is None else a2_tablei
            xt_tablei = 0 if xt_tablei is None else xt_tablei
            xc_tablei = 0 if xc_tablei is None else xc_tablei
            yt_tablei = 0 if yt_tablei is None else yt_tablei
            yc_tablei = 0 if yc_tablei is None else yc_tablei
            s_tablei = 0 if s_tablei is None else s_tablei
            ge_tablei = 0 if ge_tablei is None else ge_tablei
            f12_tablei = 0 if f12_tablei is None else f12_tablei
            material_id[i] = mid
            e1_table[i] = e1_tablei
            e2_table[i] = e2_tablei
            nu12_table[i] = nu12_tablei
            g12_table[i] = g12_tablei
            g1z_table[i] = g1z_tablei
            g2z_table[i] = g2z_tablei
            rho_table[i] = rho_tablei
            a1_table[i] = a1_tablei
            a2_table[i] = a2_tablei
            xt_table[i] = xt_tablei
            xc_table[i] = yt_tablei
            yt_table[i] = yt_tablei
            yc_table[i] = yc_tablei
            s_table[i] = s_tablei
            ge_table[i] = ge_tablei
            f12_table[i] = f12_tablei
        self._save(material_id, e1_table, e2_table, nu12_table, g12_table,
                   g1z_table, g2z_table, rho_table,
                   a1_table, a2_table, xt_table,
                   xc_table, yt_table, yc_table,
                   s_table, ge_table, f12_table)
        self.sort()
        self.cards = []

    def _save(self, material_id, e1_table, e2_table, nu12_table, g12_table,
              g1z_table, g2z_table, rho_table,
              a1_table, a2_table, xt_table,
              xc_table, yt_table, yc_table,
              s_table, ge_table, f12_table):
        if len(self.material_id):
            asdf
            #material_id = np.hstack([self.material_id, material_id])
            #E = np.hstack([self.E, E])

        self.material_id = material_id
        self.e1_table = e1_table
        self.e2_table = e2_table
        self.nu12_table = nu12_table
        self.g12_table = g12_table
        self.g1z_table = g1z_table
        self.g2z_table = g2z_table
        self.rho_table = rho_table
        self.a1_table = a1_table
        self.a2_table = a2_table
        self.xt_table = xt_table
        self.xc_table = xc_table
        self.yt_table = yt_table
        self.yc_table = yc_table
        self.s_table = s_table
        self.ge_table = ge_table
        self.f12_table = f12_table

    def set_used(self, used_dict: dict[str, list[np.ndarray]]) -> None:
        pass
        #table0 = np.hstack([
            #self.e_table, self.g_table, self.nu_table,
            #self.rho_table, self.alpha_table,
            #self.ge_table, self.st_table, self.ss_table, ])
        #utable0 = np.unique(table0)
        #table = np.setdiff1d(utable0, [0])
        #used_dict['tablem_id'].append(table)

    #def convert(self, stiffness_scale: float=1.0,
                #density_scale: float=1.0,
                #alpha_scale: float=1.0,
                #temperature_scale: float=1.0,
                #stress_scale: float=1.0, **kwargs) -> None:
        #tables0 = [
            #(self.e_table, stiffness_scale),
            #(self.g_table, stiffness_scale),
            #(self.rho_table, density_scale),
            ## nu - nope
            #(self.alpha_table, alpha_scale),
            #(self.st_table, stress_scale),
            #(self.sc_table, stress_scale),
            #(self.ss_table, stress_scale),
        #]
        #table_ids = {}
        #for table0, scale in tables0:
            #utable0 = np.unique(table0)
            #table = np.setdiff1d(utable0, [0])
            #if table:
                #table_ids[scale].append(table)
            #model.tablem

    def __apply_slice__(self, mat: MATT1, i: np.ndarray) -> None:  # ignore[override]
        mat.n = len(i)
        mat.material_id = self.material_id[i]
        mat.e1_table = self.e1_table[i]
        mat.e2_table = self.e2_table[i]
        mat.nu12_table = self.nu12_table[i]
        mat.g12_table = self.g12_table[i]
        mat.g1z_table = self.g1z_table[i]
        mat.g2z_table = self.g2z_table[i]
        mat.rho_table = self.rho_table[i]
        mat.a1_table = self.a1_table[i]
        mat.a2_table = self.a2_table[i]
        mat.xt_table = self.xt_table[i]
        mat.xc_table = self.xc_table[i]
        mat.yt_table = self.yt_table[i]
        mat.yc_table = self.yc_table[i]
        mat.s_table = self.s_table[i]
        mat.ge_table = self.ge_table[i]
        mat.f12_table = self.f12_table[i]

    def geom_check(self, missing: dict[str, np.ndarray]):
        pass

    @property
    def max_id(self):
        tables = np.hstack([self.material_id, self.e1_table, self.e2_table, self.nu12_table,
                            self.g12_table, self.g1z_table, self.g2z_table, self.rho_table,
                            self.a1_table, self.a2_table,
                            self.xt_table, self.xc_table, self.yt_table, self.yc_table,
                            self.s_table, self.ge_table, self.f12_table, ],
        )
        return tables.max()

    @parse_material_check
    def write_file(self, bdf_file: TextIOLike,
                   size: int=8, is_double: bool=False,
                   write_card_header: bool=False) -> None:
        print_card = get_print_card(size, self.max_id)

        material_id = array_str(self.material_id, size=size)
        e1_table = array_default_int(self.e1_table, default=0, size=size)
        e2_table = array_default_int(self.e2_table, default=0, size=size)
        nu12_table = array_default_int(self.nu12_table, default=0, size=size)
        g12_table = array_default_int(self.g12_table, default=0, size=size)
        g1z_table = array_default_int(self.g1z_table, default=0, size=size)
        g2z_table = array_default_int(self.g2z_table, default=0, size=size)
        rho_table = array_default_int(self.rho_table, default=0, size=size)
        a1_table = array_default_int(self.a1_table, default=0, size=size)
        a2_table = array_default_int(self.a2_table, default=0, size=size)
        xt_table = array_default_int(self.xt_table, default=0, size=size)
        xc_table = array_default_int(self.xc_table, default=0, size=size)
        yt_table = array_default_int(self.yt_table, default=0, size=size)
        yc_table = array_default_int(self.yc_table, default=0, size=size)
        s_table = array_default_int(self.s_table, default=0, size=size)
        ge_table = array_default_int(self.ge_table, default=0, size=size)
        f12_table = array_default_int(self.f12_table, default=0, size=size)

        for (mid, e1_tablei, e2_tablei, nu12_tablei,
             g12_tablei, g1z_tablei, g2z_tablei, rho_tablei,
             a1_tablei, a2_tablei,
             xt_tablei, xc_tablei, yt_tablei, yc_tablei,
             s_tablei, ge_tablei, f12_tablei) in zip_longest(
                 material_id, e1_table, e2_table, nu12_table,
                 g12_table, g1z_table, g2z_table, rho_table,
                 a1_table, a2_table,
                 xt_table, xc_table, yt_table, yc_table,
                 s_table, ge_table, f12_table):
            list_fields = ['MATT8', mid, e1_tablei, e2_tablei, nu12_tablei,
                           g12_tablei, g1z_tablei, g2z_tablei, rho_tablei,
                           a1_tablei, a2_tablei, None,
                           xt_tablei, xc_tablei, yt_tablei, yc_tablei,
                           s_tablei, ge_tablei, f12_tablei]
            bdf_file.write(print_card(list_fields))
        return


class MATT9(Material):
    """
    Specifies temperature-dependent material properties on MAT2 entry
    fields via TABLEMi entries.

    +-------+--------+--------+-------+---------+--------+--------+--------+--------+
    |   1   |   2    |   3    |   4   |    5    |   6    |   7    |    8   |   9    |
    +=======+========+========+=======+=========+========+========+========+========+
    | MATT8 |  MID   | T(E1)  | T(E2) | T(Nu12) | T(G12) | T(G1z) | T(G2z) | T(RHO) |
    +-------+--------+--------+-------+---------+--------+--------+--------+--------+
    |       |  T(A1) | T(A2)  |       |  T(Xt)  | T(Xc)  | T(Yt)  | T(Yc)  | T(S)   |
    +-------+--------+--------+-------+---------+--------+--------+--------+--------+
    |       |  T(GE) | T(F12) |       |         |        |        |        |        |
    +-------+--------+--------+-------+---------+--------+--------+--------+--------+

    """
    @Material.clear_check
    def clear(self) -> None:
        self.material_id = np.array([], dtype='int32')

    def add(self, mid: int,
            g11_table=None, g12_table=None, g13_table=None, g14_table=None,
            g15_table=None, g16_table=None, g22_table=None, g23_table=None,
            g24_table=None, g25_table=None, g26_table=None, g33_table=None,
            g34_table=None, g35_table=None, g36_table=None, g44_table=None,
            g45_table=None, g46_table=None, g55_table=None, g56_table=None,
            g66_table=None, rho_table=None,
            a1_table=None, a2_table=None, a3_table=None,
            a4_table=None, a5_table=None, a6_table=None,
            ge_table=None, comment: str='') -> int:
        """Creates a MATT9 card"""
        self.cards.append((mid,
                           g11_table, g12_table, g13_table, g14_table, g15_table, g16_table,
                           g22_table, g23_table, g24_table, g25_table, g26_table,
                           g33_table, g34_table, g35_table, g36_table,
                           g44_table, g45_table, g46_table,
                           g55_table, g56_table, g66_table, rho_table,
                           a1_table, a2_table, a3_table,
                           a4_table, a5_table, a6_table, ge_table, comment))
        self.n += 1
        return self.n - 1

    def add_card(self, card: BDFCard, comment: str='') -> int:
        """
        Adds a MATT9 card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card

        """
        mid = integer(card, 1, 'mid')
        g11_table = integer_or_blank(card, 2, 'T(G11)', default=0)
        g12_table = integer_or_blank(card, 3, 'T(G12)', default=0)
        g13_table = integer_or_blank(card, 4, 'T(G13)', default=0)
        g14_table = integer_or_blank(card, 5, 'T(G14)', default=0)
        g15_table = integer_or_blank(card, 6, 'T(G15)', default=0)
        g16_table = integer_or_blank(card, 7, 'T(G16)', default=0)

        g22_table = integer_or_blank(card, 8, 'T(G22)', default=0)
        g23_table = integer_or_blank(card, 9, 'T(G23)', default=0)
        g24_table = integer_or_blank(card, 10, 'T(G24)', default=0)
        g25_table = integer_or_blank(card, 11, 'T(G25)', default=0)
        g26_table = integer_or_blank(card, 12, 'T(G26)', default=0)

        g33_table = integer_or_blank(card, 13, 'T(G33)', default=0)
        g34_table = integer_or_blank(card, 14, 'T(G34)', default=0)
        g35_table = integer_or_blank(card, 15, 'T(G35)', default=0)
        g36_table = integer_or_blank(card, 16, 'T(G36)', default=0)

        g44_table = integer_or_blank(card, 17, 'T(G44)', default=0)
        g45_table = integer_or_blank(card, 18, 'T(G45)', default=0)
        g46_table = integer_or_blank(card, 19, 'T(G46)', default=0)

        g55_table = integer_or_blank(card, 20, 'T(G55)', default=0)
        g56_table = integer_or_blank(card, 21, 'T(G56)', default=0)
        g66_table = integer_or_blank(card, 22, 'T(G66)', default=0)

        rho_table = integer_or_blank(card, 23, 'T(RHO)', default=0)
        a1_table = integer_or_blank(card, 24, 'T(A1)', default=0)
        a2_table = integer_or_blank(card, 25, 'T(A2)', default=0)
        a3_table = integer_or_blank(card, 26, 'T(A3)', default=0)
        a4_table = integer_or_blank(card, 27, 'T(A4)', default=0)
        a5_table = integer_or_blank(card, 28, 'T(A5)', default=0)
        a6_table = integer_or_blank(card, 29, 'T(A6)', default=0)

        ge_table = integer_or_blank(card, 31, 'T(GE)', default=0)

        assert len(card) <= 32, f'len(MATT9 card) = {len(card):d}\ncard={card}'
        #return MATT9(mid, g11_table, g12_table, g13_table, g14_table, g15_table, g16_table,
                     #g22_table, g23_table, g24_table, g25_table, g26_table,
                     #g33_table, g34_table, g35_table, g36_table,
                     #g44_table, g45_table, g46_table,
                     #g55_table, g56_table, g66_table,
                     #rho_table,
                     #a1_table, a2_table, a3_table, a4_table, a5_table, a6_table,
                     #ge_table, comment=comment)

        self.cards.append((mid, g11_table, g12_table, g13_table, g14_table, g15_table, g16_table,
                           g22_table, g23_table, g24_table, g25_table, g26_table,
                           g33_table, g34_table, g35_table, g36_table,
                           g44_table, g45_table, g46_table,
                           g55_table, g56_table, g66_table, rho_table,
                           a1_table, a2_table, a3_table,
                           a4_table, a5_table, a6_table, ge_table, comment))
        self.n += 1
        return self.n - 1

    @Material.parse_cards_check
    def parse_cards(self) -> None:
        ncards = len(self.cards)
        idtype = self.model.idtype
        fdtype = self.model.fdtype
        material_id = np.zeros(ncards, dtype=idtype)

        g11_table = np.zeros(ncards, dtype=fdtype)
        g12_table = np.zeros(ncards, dtype=fdtype)
        g13_table = np.zeros(ncards, dtype=fdtype)
        g14_table = np.zeros(ncards, dtype=fdtype)
        g15_table = np.zeros(ncards, dtype=fdtype)
        g16_table = np.zeros(ncards, dtype=fdtype)
        g22_table = np.zeros(ncards, dtype=fdtype)
        g23_table = np.zeros(ncards, dtype=fdtype)
        g24_table = np.zeros(ncards, dtype=fdtype)
        g25_table = np.zeros(ncards, dtype=fdtype)
        g26_table = np.zeros(ncards, dtype=fdtype)
        g33_table = np.zeros(ncards, dtype=fdtype)
        g34_table = np.zeros(ncards, dtype=fdtype)
        g35_table = np.zeros(ncards, dtype=fdtype)
        g36_table = np.zeros(ncards, dtype=fdtype)
        g44_table = np.zeros(ncards, dtype=fdtype)
        g45_table = np.zeros(ncards, dtype=fdtype)
        g46_table = np.zeros(ncards, dtype=fdtype)
        g55_table = np.zeros(ncards, dtype=fdtype)
        g56_table = np.zeros(ncards, dtype=fdtype)
        g66_table = np.zeros(ncards, dtype=fdtype)
        rho_table = np.zeros(ncards, dtype=fdtype)
        a1_table = np.zeros(ncards, dtype=fdtype)
        a2_table = np.zeros(ncards, dtype=fdtype)
        a3_table = np.zeros(ncards, dtype=fdtype)
        a4_table = np.zeros(ncards, dtype=fdtype)
        a5_table = np.zeros(ncards, dtype=fdtype)
        a6_table = np.zeros(ncards, dtype=fdtype)
        ge_table = np.zeros(ncards, dtype=fdtype)
        for i, card in enumerate(self.cards):
            (mid,
             g11_tablei, g12_tablei, g13_tablei, g14_tablei, g15_tablei, g16_tablei,
             g22_tablei, g23_tablei, g24_tablei, g25_tablei, g26_tablei,
             g33_tablei, g34_tablei, g35_tablei, g36_tablei,
             g44_tablei, g45_tablei, g46_tablei,
             g55_tablei, g56_tablei, g66_tablei,
             rho_tablei,
             a1_tablei, a2_tablei, a3_tablei,
             a4_tablei, a5_tablei, a6_tablei,
             ge_tablei, comment) = card

            g11_tablei = 0 if g11_tablei is None else g11_tablei
            g12_tablei = 0 if g12_tablei is None else g12_tablei
            g13_tablei = 0 if g13_tablei is None else g13_tablei
            g14_tablei = 0 if g14_tablei is None else g14_tablei
            g15_tablei = 0 if g15_tablei is None else g15_tablei
            g16_tablei = 0 if g16_tablei is None else g16_tablei

            g22_tablei = 0 if g22_tablei is None else g22_tablei
            g23_tablei = 0 if g23_tablei is None else g23_tablei
            g24_tablei = 0 if g24_tablei is None else g24_tablei
            g25_tablei = 0 if g25_tablei is None else g25_tablei
            g26_tablei = 0 if g26_tablei is None else g26_tablei

            g33_tablei = 0 if g33_tablei is None else g33_tablei
            g34_tablei = 0 if g34_tablei is None else g34_tablei
            g35_tablei = 0 if g35_tablei is None else g35_tablei
            g36_tablei = 0 if g36_tablei is None else g36_tablei

            g44_tablei = 0 if g44_tablei is None else g44_tablei
            g45_tablei = 0 if g45_tablei is None else g45_tablei
            g46_tablei = 0 if g46_tablei is None else g46_tablei
            g55_tablei = 0 if g55_tablei is None else g55_tablei
            g56_tablei = 0 if g56_tablei is None else g56_tablei
            g66_tablei = 0 if g66_tablei is None else g66_tablei
            rho_tablei = 0 if rho_tablei is None else rho_tablei

            a1_tablei = 0 if a1_tablei is None else a1_tablei
            a2_tablei = 0 if a2_tablei is None else a2_tablei
            a3_tablei = 0 if a3_tablei is None else a3_tablei
            a4_tablei = 0 if a4_tablei is None else a4_tablei
            a5_tablei = 0 if a5_tablei is None else a5_tablei
            a6_tablei = 0 if a6_tablei is None else a6_tablei

            ge_tablei = 0 if ge_tablei is None else ge_tablei
            material_id[i] = mid

            g11_table[i] = g11_tablei
            g12_table[i] = g12_tablei
            g13_table[i] = g13_tablei
            g14_table[i] = g14_tablei
            g15_table[i] = g15_tablei
            g16_table[i] = g16_tablei
            g22_table[i] = g22_tablei
            g23_table[i] = g23_tablei
            g24_table[i] = g24_tablei
            g25_table[i] = g25_tablei
            g26_table[i] = g26_tablei
            g33_table[i] = g33_tablei
            g34_table[i] = g34_tablei
            g35_table[i] = g35_tablei
            g36_table[i] = g36_tablei
            g44_table[i] = g44_tablei
            g45_table[i] = g45_tablei
            g46_table[i] = g46_tablei
            g55_table[i] = g55_tablei
            g56_table[i] = g56_tablei
            g66_table[i] = g66_tablei
            rho_table[i] = rho_tablei
            a1_table[i] = a1_tablei
            a2_table[i] = a2_tablei
            a3_table[i] = a3_tablei
            a4_table[i] = a4_tablei
            a5_table[i] = a5_tablei
            a6_table[i] = a6_tablei
            ge_table[i] = ge_tablei
        self._save(material_id,
                   g11_table, g12_table, g13_table, g14_table, g15_table, g16_table,
                   g22_table, g23_table, g24_table, g25_table, g26_table,
                   g33_table, g34_table, g35_table, g36_table,
                   g44_table, g45_table, g46_table,
                   g55_table, g56_table, g66_table, rho_table,
                   a1_table, a2_table, a3_table,
                   a4_table, a5_table, a6_table, ge_table)
        self.sort()
        self.cards = []

    def _save(self, material_id,
              g11_table, g12_table, g13_table, g14_table, g15_table, g16_table,
              g22_table, g23_table, g24_table, g25_table, g26_table,
              g33_table, g34_table, g35_table, g36_table,
              g44_table, g45_table, g46_table,
              g55_table, g56_table, g66_table, rho_table,
              a1_table, a2_table, a3_table,
              a4_table, a5_table, a6_table, ge_table):
        if len(self.material_id):
            asdf
            #material_id = np.hstack([self.material_id, material_id])
            #E = np.hstack([self.E, E])

        self.material_id = material_id
        self.g11_table = g11_table
        self.g12_table = g12_table
        self.g13_table = g13_table
        self.g14_table = g14_table
        self.g15_table = g15_table
        self.g16_table = g16_table
        self.g22_table = g22_table
        self.g23_table = g23_table
        self.g24_table = g24_table
        self.g25_table = g25_table
        self.g26_table = g26_table
        self.g33_table = g33_table
        self.g34_table = g34_table
        self.g35_table = g35_table
        self.g36_table = g36_table
        self.g44_table = g44_table
        self.g45_table = g45_table
        self.g46_table = g46_table
        self.g55_table = g55_table
        self.g56_table = g56_table
        self.g66_table = g66_table
        self.rho_table = rho_table

        self.a1_table = a1_table
        self.a2_table = a2_table
        self.a3_table = a3_table
        self.a4_table = a4_table
        self.a5_table = a5_table
        self.a6_table = a6_table
        self.ge_table = ge_table

    def set_used(self, used_dict: dict[str, list[np.ndarray]]) -> None:
        pass
        #table0 = np.hstack([
            #self.e_table, self.g_table, self.nu_table,
            #self.rho_table, self.alpha_table,
            #self.ge_table, self.st_table, self.ss_table, ])
        #utable0 = np.unique(table0)
        #table = np.setdiff1d(utable0, [0])
        #used_dict['tablem_id'].append(table)

    #def convert(self, stiffness_scale: float=1.0,
                #density_scale: float=1.0,
                #alpha_scale: float=1.0,
                #temperature_scale: float=1.0,
                #stress_scale: float=1.0, **kwargs) -> None:
        #tables0 = [
            #(self.e_table, stiffness_scale),
            #(self.g_table, stiffness_scale),
            #(self.rho_table, density_scale),
            ## nu - nope
            #(self.alpha_table, alpha_scale),
            #(self.st_table, stress_scale),
            #(self.sc_table, stress_scale),
            #(self.ss_table, stress_scale),
        #]
        #table_ids = {}
        #for table0, scale in tables0:
            #utable0 = np.unique(table0)
            #table = np.setdiff1d(utable0, [0])
            #if table:
                #table_ids[scale].append(table)
            #model.tablem

    def __apply_slice__(self, mat: MATT9, i: np.ndarray) -> None:  # ignore[override]
        mat.n = len(i)
        mat.material_id = self.material_id[i]

        mat.g11_table = self.g11_table[i]
        mat.g12_table = self.g12_table[i]
        mat.g13_table = self.g13_table[i]
        mat.g14_table = self.g14_table[i]
        mat.g15_table = self.g15_table[i]
        mat.g16_table = self.g16_table[i]
        mat.g22_table = self.g22_table[i]
        mat.g23_table = self.g23_table[i]
        mat.g24_table = self.g24_table[i]
        mat.g25_table = self.g25_table[i]
        mat.g26_table = self.g26_table[i]
        mat.g33_table = self.g33_table[i]
        mat.g34_table = self.g34_table[i]
        mat.g35_table = self.g35_table[i]
        mat.g36_table = self.g36_table[i]
        mat.g44_table = self.g44_table[i]
        mat.g45_table = self.g45_table[i]
        mat.g46_table = self.g46_table[i]
        mat.g55_table = self.g55_table[i]
        mat.g56_table = self.g56_table[i]
        mat.g66_table = self.g66_table[i]
        mat.rho_table = self.rho_table[i]

        mat.a1_table = self.a1_table[i]
        mat.a2_table = self.a2_table[i]
        mat.a3_table = self.a3_table[i]
        mat.a4_table = self.a4_table[i]
        mat.a5_table = self.a5_table[i]
        mat.a6_table = self.a6_table[i]

        mat.ge_table = self.ge_table[i]

    def geom_check(self, missing: dict[str, np.ndarray]):
        pass

    @property
    def max_id(self):
        tables = np.hstack([
            self.material_id,
            self.g11_table, self.g12_table, self.g13_table, self.g14_table, self.g15_table, self.g16_table,
            self.g22_table, self.g23_table, self.g24_table, self.g25_table, self.g26_table,
            self.g33_table, self.g34_table, self.g35_table, self.g36_table,
            self.g44_table, self.g45_table, self.g46_table,
            self.g55_table, self.g56_table, self.g66_table, self.rho_table,
            self.a1_table, self.a2_table, self.a3_table,
            self.a4_table, self.a5_table, self.a6_table, self.ge_table])
        return tables.max()

    @parse_material_check
    def write_file(self, bdf_file: TextIOLike,
                   size: int=8, is_double: bool=False,
                   write_card_header: bool=False) -> None:
        print_card = get_print_card(size, self.max_id)

        material_id = array_str(self.material_id, size=size)

        g11_table = array_default_int(self.g11_table, default=0, size=size)
        g12_table = array_default_int(self.g12_table, default=0, size=size)
        g13_table = array_default_int(self.g13_table, default=0, size=size)
        g14_table = array_default_int(self.g14_table, default=0, size=size)
        g15_table = array_default_int(self.g15_table, default=0, size=size)
        g16_table = array_default_int(self.g16_table, default=0, size=size)
        g22_table = array_default_int(self.g22_table, default=0, size=size)
        g23_table = array_default_int(self.g23_table, default=0, size=size)
        g24_table = array_default_int(self.g24_table, default=0, size=size)
        g25_table = array_default_int(self.g25_table, default=0, size=size)
        g26_table = array_default_int(self.g26_table, default=0, size=size)
        g33_table = array_default_int(self.g33_table, default=0, size=size)
        g34_table = array_default_int(self.g34_table, default=0, size=size)
        g35_table = array_default_int(self.g35_table, default=0, size=size)
        g36_table = array_default_int(self.g36_table, default=0, size=size)
        g44_table = array_default_int(self.g44_table, default=0, size=size)
        g45_table = array_default_int(self.g45_table, default=0, size=size)
        g46_table = array_default_int(self.g46_table, default=0, size=size)
        g55_table = array_default_int(self.g55_table, default=0, size=size)
        g56_table = array_default_int(self.g56_table, default=0, size=size)
        g66_table = array_default_int(self.g66_table, default=0, size=size)
        rho_table = array_default_int(self.rho_table, default=0, size=size)
        a1_table = array_default_int(self.a1_table, default=0, size=size)
        a2_table = array_default_int(self.a2_table, default=0, size=size)
        a3_table = array_default_int(self.a3_table, default=0, size=size)
        a4_table = array_default_int(self.a4_table, default=0, size=size)
        a5_table = array_default_int(self.a5_table, default=0, size=size)
        a6_table = array_default_int(self.a6_table, default=0, size=size)

        ge_table = array_default_int(self.ge_table, default=0, size=size)

        for (mid,
             g11_tablei, g12_tablei, g13_tablei, g14_tablei, g15_tablei, g16_tablei,
             g22_tablei, g23_tablei, g24_tablei, g25_tablei, g26_tablei,
             g33_tablei, g34_tablei, g35_tablei, g36_tablei,
             g44_tablei, g45_tablei, g46_tablei,
             g55_tablei, g56_tablei, g66_tablei, rho_tablei,
             a1_tablei, a2_tablei, a3_tablei, a4_tablei, a5_tablei, a6_tablei,
             ge_tablei,) in zip_longest(
                 material_id,
                 g11_table, g12_table, g13_table, g14_table, g15_table, g16_table,
                 g22_table, g23_table, g24_table, g25_table, g26_table,
                 g33_table, g34_table, g35_table, g36_table,
                 g44_table, g45_table, g46_table,
                 g55_table, g56_table, g66_table, rho_table,
                 a1_table, a2_table, a3_table, a4_table, a5_table, a6_table, ge_table):
            list_fields = [
                'MATT9', mid,
                g11_tablei, g12_tablei, g13_tablei, g14_tablei, g15_tablei, g16_tablei,
                g22_tablei, g23_tablei, g24_tablei, g25_tablei, g26_tablei,
                g33_tablei, g34_tablei, g35_tablei, g36_tablei,
                g44_tablei, g45_tablei, g46_tablei,
                g55_tablei, g56_tablei, g66_tablei,
                rho_tablei,
                a1_tablei, a2_tablei, a3_tablei, a4_tablei, a5_tablei, a6_tablei,
                ge_tablei,
            ]

            bdf_file.write(print_card(list_fields))
        return
