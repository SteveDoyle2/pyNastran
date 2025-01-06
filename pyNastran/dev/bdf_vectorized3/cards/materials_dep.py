from __future__ import annotations
from itertools import zip_longest
from typing import TYPE_CHECKING
import numpy as np

from pyNastran.bdf.bdf_interface.assign_type import (
    string, integer, double,
    integer_or_blank, double_or_blank, string_or_blank,
)
from pyNastran.bdf.bdf_interface.assign_type_force import force_double_or_blank
#from pyNastran.bdf.cards.materials import mat1_E_G_nu, get_G_default, set_blank_if_default

from pyNastran.dev.bdf_vectorized3.cards.base_card import Material, parse_check
from pyNastran.dev.bdf_vectorized3.cards.write_utils import (
    get_print_card_size, array_str, array_default_int, array_float_nan)

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
        return self.n - 1

    def add_card(self, card: BDFCard, ifile: int, comment: str='') -> int:
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
        return self.n - 1

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

        for icard, card in enumerate(self.cards):
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
            material_id[icard] = mid
            e_table[icard] = e_tablei
            g_table[icard] = g_tablei
            nu_table[icard] = nu_tablei
            rho_table[icard] = rho_tablei
            alpha_table[icard] = a_tablei
            ge_table[icard] = ge_tablei
            st_table[icard] = st_tablei
            sc_table[icard] = sc_tablei
            ss_table[icard] = ss_tablei
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
            raise RuntimeError(f'stacking of {self.type} is not supported')
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
    def max_id(self) -> int:
        tables = np.hstack([self.material_id,
                            self.e_table, self.g_table, self.nu_table,
                            self.rho_table,
                            self.alpha_table,
                            self.st_table, self.sc_table, self.ss_table,
                            self.ge_table, ])
        return tables.max()

    @parse_check
    def write_file(self, bdf_file: TextIOLike,
                   size: int=8, is_double: bool=False,
                   write_card_header: bool=False) -> None:
        print_card, size = get_print_card_size(size, self.max_id)

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
        return self.n - 1

    def add_card(self, card: BDFCard, ifile: int, comment: str='') -> int:
        """
        Adds a MATS1 card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card

        """
        map_type_dict = {
            'NLELAS': 'NLELAST',
        }
        mid = integer(card, 1, 'mid')
        tid = integer_or_blank(card, 2, 'tables1_id')
        nonlinear_type = string(card, 3, 'Type')

        nonlinear_type = map_type_dict.get(nonlinear_type, nonlinear_type)
        if nonlinear_type not in {'NLELAST', 'PLASTIC', 'PLSTRN'}:
            raise ValueError('MATS1 Type must be [NLELAST, PLASTIC, PLSTRN]; Type=%r' % nonlinear_type)

        if nonlinear_type == 'NLELAST':
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
            fdouble_or_blank = double_or_blank if self.model.is_strict_card_parser else force_double_or_blank
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
        self.cards.append((mid, tid, nonlinear_type, hardening_slope, hr, yf,
                           limit1, limit2, stress_strain_measure, comment))
        self.n += 1
        return self.n - 1

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

        for icard, card in enumerate(self.cards):
            (mid, tid, nonlinear_typei, hardening_slopei, hri, yfi,
             limit1i, limit2i, stress_strain_measurei, comment) = card
            tid = 0 if tid is None else tid
            hri = 1 if hri is None else hri
            yfi = 1 if yfi is None else yfi
            material_id[icard] = mid
            table_id[icard] = tid
            Type[icard] = nonlinear_typei
            hardening_slope[icard] = hardening_slopei
            hr[icard] = hri
            yf[icard] = yfi
            limit1[icard] = limit1i
            limit2[icard] = limit2i
            stress_strain_measure[icard] = stress_strain_measurei
        self._save(material_id, table_id, Type, hardening_slope, hr, yf,
                   limit1, limit2, stress_strain_measure)
        self.sort()
        self.cards = []

    def _save(self, material_id, table_id, Type, hardening_slope, hr, yf,
              limit1, limit2, stress_strain_measure):
        if len(self.material_id):
            raise RuntimeError(f'stacking of {self.type} is not supported')
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
    def max_id(self) -> int:
        return max(self.material_id.max(), self.table_id.max())

    @parse_check
    def write_file(self, bdf_file: TextIOLike,
                   size: int=8, is_double: bool=False,
                   write_card_header: bool=False) -> None:
        print_card, size = get_print_card_size(size, self.max_id)

        material_id = array_str(self.material_id, size=size)
        table_id = array_str(self.table_id, size=size)
        hardening_slope = array_float_nan(self.hardening_slope, size=size)
        yf = array_str(self.yf, size=size)
        hr = array_str(self.hr, size=size)
        limit1 = array_float_nan(self.limit1, size=size)
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
        return self.n - 1

    def add_card(self, card: BDFCard, ifile: int, comment: str='') -> int:
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
        return self.n - 1

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
        for icard, card in enumerate(self.cards):
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
            material_id[icard] = mid
            g11_table[icard] = g11_tablei
            g22_table[icard] = g12_tablei
            g33_table[icard] = g13_tablei

            g12_table[icard] = g12_tablei
            g13_table[icard] = g13_tablei
            g23_table[icard] = g23_tablei

            rho_table[icard] = rho_tablei
            a1_table[icard] = a1_tablei
            a2_table[icard] = a2_tablei
            a3_table[icard] = a3_tablei

            st_table[icard] = st_tablei
            sc_table[icard] = sc_tablei
            ss_table[icard] = ss_tablei
            ge_table[icard] = ge_tablei
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
            raise RuntimeError(f'stacking of {self.type} is not supported')
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
    def max_id(self) -> int:
        tables = np.hstack([self.material_id,
                            self.g11_table, self.g22_table, self.g33_table,
                            self.g12_table, self.g13_table, self.g23_table,
                            self.rho_table,
                            self.a1_table, self.a2_table, self.a3_table,
                            self.st_table, self.sc_table, self.ss_table,
                            self.ge_table, ])
        return tables.max()

    @parse_check
    def write_file(self, bdf_file: TextIOLike,
                   size: int=8, is_double: bool=False,
                   write_card_header: bool=False) -> None:
        print_card, size = get_print_card_size(size, self.max_id)

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


class MATT3(Material):
    @Material.clear_check
    def clear(self) -> None:
        self.material_id = np.array([], dtype='int32')
        #self.g11_table = np.array([], dtype='int32')
        #self.g22_table = np.array([], dtype='int32')
        #self.g33_table = np.array([], dtype='int32')
        #self.g12_table = np.array([], dtype='int32')
        #self.g13_table = np.array([], dtype='int32')
        #self.g23_table = np.array([], dtype='int32')
        #self.rho_table = np.array([], dtype='int32')
        #self.a1_table = np.array([], dtype='int32')
        #self.a2_table = np.array([], dtype='int32')
        #self.a3_table = np.array([], dtype='int32')
        #self.st_table = np.array([], dtype='int32')
        #self.sc_table = np.array([], dtype='int32')
        #self.ss_table = np.array([], dtype='int32')
        #self.ge_table = np.array([], dtype='int32')

    def add(self, mid: int,
            ex_table: int=0, eth_table: int=0, ez_table: int=0,
            nuth_table: int=0, nuxz_table: int=0, rho_table: int=0,
            gzx_table: int=0,
            ax_table: int=0, ath_table: int=0, az_table: int=0,
            ge_table: int=0, comment: str='') -> int:
        """Creates a MATT2 card"""
        self.cards.append((mid, ex_table, eth_table, ez_table,
                           nuth_table, nuxz_table, rho_table, gzx_table,
                           ax_table, ath_table, az_table, ge_table, comment))
        self.n += 1
        return self.n - 1

    def add_card(self, card: BDFCard, ifile: int, comment: str='') -> int:
        """
        Adds a MATT3 card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card

        """
        mid = integer(card, 1, 'mid')
        ex_table = integer_or_blank(card, 2, 'T(EX)', default=0)
        eth_table = integer_or_blank(card, 3, 'T(ETH)', default=0)
        ez_table = integer_or_blank(card, 5, 'T(EZ)', default=0)
        nuth_table = integer_or_blank(card, 6, 'T(NUTH)', default=0)
        nuxz_table = integer_or_blank(card, 7, 'T(NUXZ)', default=0)
        rho_table = integer_or_blank(card, 8, 'T(RHO)', default=0)

        gzx_table = integer_or_blank(card, 11, 'T(GZX)', default=0)
        ax_table = integer_or_blank(card, 12, 'T(AX)', default=0)
        ath_table = integer_or_blank(card, 13, 'T(ATH)', default=0)
        az_table = integer_or_blank(card, 14, 'T(AZ)', default=0)
        ge_table = integer_or_blank(card, 16, 'T(GE)', default=0)

        assert len(card) <= 16, f'len(MATT3 card) = {len(card):d}\ncard={card}'
        #return MATT3(mid, ex_table, eth_table, ez_table,
                     #nuth_table, nuxz_table, rho_table, gzx_table,
                     #ax_table, ath_table, az_table, ge_table, comment=comment)
        self.cards.append((mid, ex_table, eth_table, ez_table,
                           nuth_table, nuxz_table, rho_table, gzx_table,
                           ax_table, ath_table, az_table, ge_table, comment))
        self.n += 1
        return self.n - 1

    @Material.parse_cards_check
    def parse_cards(self) -> None:
        ncards = len(self.cards)
        idtype = self.model.idtype
        material_id = np.zeros(ncards, dtype=idtype)

        ex_table = np.zeros(ncards, dtype=idtype)
        eth_table = np.zeros(ncards, dtype=idtype)
        ez_table = np.zeros(ncards, dtype=idtype)
        nuth_table = np.zeros(ncards, dtype=idtype)
        nuxz_table = np.zeros(ncards, dtype=idtype)
        rho_table = np.zeros(ncards, dtype=idtype)
        gzx_table = np.zeros(ncards, dtype=idtype)
        ax_table = np.zeros(ncards, dtype=idtype)
        ath_table = np.zeros(ncards, dtype=idtype)
        az_table = np.zeros(ncards, dtype=idtype)
        ge_table = np.zeros(ncards, dtype=idtype)
        for icard, card in enumerate(self.cards):
            (mid,
             ex_tablei, eth_tablei, ez_tablei,
             nuth_tablei, nuxz_tablei, rho_tablei, gzx_tablei,
             ax_tablei, ath_tablei, az_tablei, ge_tablei, comment) = card

            ex_tablei = 0 if ex_tablei is None else ex_tablei
            eth_tablei = 0 if eth_tablei is None else eth_tablei
            ez_tablei = 0 if ez_tablei is None else ez_tablei

            nuth_tablei = 0 if nuth_tablei is None else nuth_tablei
            nuxz_tablei = 0 if nuxz_tablei is None else nuxz_tablei
            gzx_tablei = 0 if gzx_tablei is None else gzx_tablei
            rho_tablei = 0 if rho_tablei is None else rho_tablei

            ax_tablei = 0 if ax_tablei is None else ax_tablei
            ath_tablei = 0 if ath_tablei is None else ath_tablei
            az_tablei = 0 if az_tablei is None else az_tablei

            ge_tablei = 0 if ge_tablei is None else ge_tablei
            material_id[icard] = mid
            ex_table[icard] = ex_tablei
            eth_table[icard] = eth_tablei
            ez_table[icard] = ez_tablei

            nuth_table[icard] = nuth_tablei
            nuxz_table[icard] = nuxz_tablei
            gzx_table[icard] = gzx_tablei

            rho_table[icard] = rho_tablei
            ax_table[icard] = ax_tablei
            ath_table[icard] = ath_tablei
            az_table[icard] = az_tablei
            ge_table[icard] = ge_tablei
        self._save(material_id,
                   ex_table, eth_table, ez_table,
                   nuth_table, nuxz_table, rho_table, gzx_table,
                   ax_table, ath_table, az_table, ge_table)
        self.sort()
        self.cards = []

    def _save(self, material_id,
              ex_table, eth_table, ez_table,
              nuth_table, nuxz_table, rho_table, gzx_table,
              ax_table, ath_table, az_table, ge_table):
        if len(self.material_id):
            raise RuntimeError(f'stacking of {self.type} is not supported')
            #material_id = np.hstack([self.material_id, material_id])
            #E = np.hstack([self.E, E])

        self.material_id = material_id
        self.ex_table = ex_table
        self.eth_table = eth_table
        self.ez_table = ez_table
        self.nuth_table = nuth_table
        self.nuxz_table = nuxz_table
        self.gzx_table = gzx_table
        self.rho_table = rho_table
        self.ax_table = ax_table
        self.ath_table = ath_table
        self.az_table = az_table
        self.ge_table = ge_table

    def set_used(self, used_dict: dict[str, list[np.ndarray]]) -> None:
        table0 = np.hstack([
            self.ex_table, self.eth_table, self.ez_table,
            self.nuth_table, self.nuxz_table,
            self.gzx_table, self.rho_table,
            self.ax_table, self.ath_table, self.az_table,
            self.ge_table,])
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

    def __apply_slice__(self, mat: MATT3, i: np.ndarray) -> None:  # ignore[override]
        mat.n = len(i)
        mat.material_id = self.material_id[i]
        mat.ex_table = self.ex_table[i]
        mat.eth_table = self.eth_table[i]
        mat.ez_table = self.ez_table[i]
        mat.nuth_table = self.nuth_table[i]
        mat.nuxz_table = self.nuxz_table[i]
        mat.gzx_table = self.gzx_table[i]
        mat.rho_table = self.rho_table[i]
        mat.ax_table = self.ax_table[i]
        mat.ath_table = self.ath_table[i]
        mat.az_table = self.az_table[i]
        mat.ge_table = self.ge_table[i]

    def geom_check(self, missing: dict[str, np.ndarray]):
        pass

    @property
    def max_id(self) -> int:
        tables = np.hstack([self.material_id,
                            self.ex_table, self.eth_table, self.ez_table,
                            self.nuth_table, self.nuxz_table,
                            self.gzx_table, self.rho_table,
                            self.ax_table, self.ath_table, self.az_table,
                            self.ge_table,])
        return tables.max()

    @parse_check
    def write_file(self, bdf_file: TextIOLike,
                   size: int=8, is_double: bool=False,
                   write_card_header: bool=False) -> None:
        print_card, size = get_print_card_size(size, self.max_id)

        material_id = array_str(self.material_id, size=size)
        ex_table = array_default_int(self.ex_table, default=0, size=size)
        eth_table = array_default_int(self.eth_table, default=0, size=size)
        ez_table = array_default_int(self.ez_table, default=0, size=size)
        nuth_table = array_default_int(self.nuth_table, default=0, size=size)
        nuxz_table = array_default_int(self.nuxz_table, default=0, size=size)
        gzx_table = array_default_int(self.gzx_table, default=0, size=size)
        rho_table = array_default_int(self.rho_table, default=0, size=size)
        ax_table = array_default_int(self.ax_table, default=0, size=size)
        ath_table = array_default_int(self.ath_table, default=0, size=size)
        az_table = array_default_int(self.az_table, default=0, size=size)
        ge_table = array_default_int(self.ge_table, default=0, size=size)

        for (mid, ex_tablei, eth_tablei, ez_tablei,
             nuth_tablei, nuxz_tablei, rho_tablei, gzx_tablei,
             ax_tablei, ath_tablei, az_tablei, ge_tablei) in zip_longest(
                 material_id,
                 ex_table, eth_table, ez_table,
                 nuth_table, nuxz_table, rho_table, gzx_table,
                 ax_table, ath_table, az_table, ge_table):

            list_fields = [
                'MATT3', mid, ex_tablei, eth_tablei, ez_tablei,
                nuth_tablei, nuxz_tablei, rho_tablei, None, None,
                gzx_tablei, ax_tablei, ath_tablei, az_tablei,
                None, ge_tablei,]
            bdf_file.write(print_card(list_fields))
        return


class MATT4(Material):
    """
    Specifies temperature-dependent material properties on MAT2 entry
    fields via TABLEMi entries.

    +-------+-------+-------+-------+--------+-------+-------+---------+
    |   1   |   2   |   3   |   4   |   5    |   6   |   7   |    8    |
    +=======+=======+=======+=======+========+=======+=======+=========+
    | MATT4 |  MID  |  T(K) | T(CP) |        | T(H)  | T(mu) | T(HGEN) |
    +-------+-------+-------+-------+--------+-------+-------+---------+

    """
    @Material.clear_check
    def clear(self) -> None:
        self.material_id = np.array([], dtype='int32')
        self.k_table = np.array([], dtype='int32')
        self.cp_table = np.array([], dtype='int32')
        self.h_table = np.array([], dtype='int32')
        self.mu_table = np.array([], dtype='int32')
        self.hgen_table = np.array([], dtype='int32')

    def add(self, mid: int, k_table: int=0, cp_table: int=0, h_table: int=0,
            mu_table: int=0, hgen_table: int=0, comment: str='') -> int:
        """Creates a MATT4 card"""
        self.cards.append((mid, k_table, cp_table, h_table,
                           mu_table, hgen_table, comment))
        self.n += 1
        return self.n - 1

    def add_card(self, card: BDFCard, ifile: int, comment: str='') -> int:
        """
        Adds a MATT4 card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card

        """
        mid = integer(card, 1, 'mid')
        k_table = integer_or_blank(card, 2, 'T(K)', default=0)
        cp_table = integer_or_blank(card, 3, 'T(CP)', default=0)
        h_table = integer_or_blank(card, 5, 'T(H)', default=0)
        mu_table = integer_or_blank(card, 6, 'T(mu)', default=0)
        hgen_table = integer_or_blank(card, 7, 'T(HGEN)', default=0)

        assert len(card) <= 8, 'len(MATT4 card) = {len(card):d}\ncard={card}'
        #return MATT4(mid, k_table, cp_table, h_table, mu_table,
                     #hgen_table, comment=comment)
        self.cards.append((mid, k_table, cp_table, h_table,
                           mu_table, hgen_table, comment))
        self.n += 1
        return self.n - 1

    @Material.parse_cards_check
    def parse_cards(self) -> None:
        ncards = len(self.cards)
        idtype = self.model.idtype
        material_id = np.zeros(ncards, dtype=idtype)

        k_table = np.zeros(ncards, dtype=idtype)
        cp_table = np.zeros(ncards, dtype=idtype)
        h_table = np.zeros(ncards, dtype=idtype)
        mu_table = np.zeros(ncards, dtype=idtype)
        hgen_table = np.zeros(ncards, dtype=idtype)
        for icard, card in enumerate(self.cards):
            (mid, k_tablei, cp_tablei, h_tablei, mu_tablei, hgen_tablei, comment) = card

            k_tablei = 0 if k_tablei is None else k_tablei
            cp_tablei = 0 if cp_tablei is None else cp_tablei
            h_tablei = 0 if h_tablei is None else h_tablei
            mu_tablei = 0 if mu_tablei is None else mu_tablei
            hgen_tablei = 0 if hgen_tablei is None else hgen_tablei

            material_id[icard] = mid
            k_table[icard] = k_tablei
            cp_table[icard] = cp_tablei
            h_table[icard] = h_tablei

            mu_table[icard] = mu_tablei
            hgen_table[icard] = hgen_tablei
        self._save(material_id, k_table, cp_table, h_table, mu_table, hgen_table)
        self.sort()
        self.cards = []

    def _save(self, material_id, k_table, cp_table, h_table, mu_table, hgen_table):
        if len(self.material_id):
            raise RuntimeError(f'stacking of {self.type} is not supported')
            #material_id = np.hstack([self.material_id, material_id])
            #E = np.hstack([self.E, E])

        self.material_id = material_id
        self.k_table = k_table
        self.cp_table = cp_table
        self.h_table = h_table
        self.mu_table = mu_table
        self.hgen_table = hgen_table

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

    def __apply_slice__(self, mat: MATT3, i: np.ndarray) -> None:  # ignore[override]
        mat.n = len(i)
        mat.material_id = self.material_id[i]
        mat.k_table = self.k_table[i]
        mat.cp_table = self.cp_table[i]
        mat.h_table = self.h_table[i]
        mat.mu_table = self.mu_table[i]
        mat.hgen_table = self.hgen_table[i]

    def geom_check(self, missing: dict[str, np.ndarray]):
        pass

    @property
    def max_id(self) -> int:
        tables = np.hstack([self.material_id,
                            self.k_table, self.cp_table, self.h_table,
                            self.mu_table, self.hgen_table, ])
        return tables.max()

    @parse_check
    def write_file(self, bdf_file: TextIOLike,
                   size: int=8, is_double: bool=False,
                   write_card_header: bool=False) -> None:
        print_card, size = get_print_card_size(size, self.max_id)

        material_id = array_str(self.material_id, size=size)
        k_table = array_default_int(self.k_table, default=0, size=size)
        cp_table = array_default_int(self.cp_table, default=0, size=size)
        h_table = array_default_int(self.h_table, default=0, size=size)
        mu_table = array_default_int(self.mu_table, default=0, size=size)
        hgen_table = array_default_int(self.hgen_table, default=0, size=size)

        for (mid, k_tablei, cp_tablei, h_tablei, mu_tablei, hgen_tablei) in zip_longest(
            material_id, k_table, cp_table, h_table, mu_table, hgen_table):

            list_fields = [
                'MATT4', mid, k_tablei, cp_tablei,
                None, h_tablei, mu_tablei, hgen_tablei,
            ]
            bdf_file.write(print_card(list_fields))
        return


class MATT5(Material):
    """
    Specifies temperature-dependent material properties on MAT2 entry
    fields via TABLEMi entries.

    +-------+---------+---------+--------+--------+--------+--------+--------+-------+
    |   1   |    2    |    3    |   4    |   5    |   6    |   7    |    8   |   9   |
    +=======+=========+=========+========+========+========+========+========+=======+
    | MATT5 |   MID   | T(Kxx)  | T(Kxy) | T(Kxz) | T(Kyy) | T(Kyz) | T(Kzz) | T(CP) |
    +-------+---------+---------+--------+--------+--------+--------+--------+-------+
    |       |         | T(HGEN) |        |        |        |        |        |       |
    +-------+---------+---------+--------+--------+--------+--------+--------+-------+

    """
    @Material.clear_check
    def clear(self) -> None:
        self.material_id = np.array([], dtype='int32')
        self.kxx_table = np.array([], dtype='int32')
        self.kxy_table = np.array([], dtype='int32')
        self.kxz_table = np.array([], dtype='int32')
        self.kyy_table = np.array([], dtype='int32')
        self.kyz_table = np.array([], dtype='int32')
        self.kzz_table = np.array([], dtype='int32')
        self.cp_table = np.array([], dtype='int32')
        self.hgen_table = np.array([], dtype='int32')

    def add(self, mid, kxx_table: int=0, kxy_table: int=0, kxz_table: int=0,
            kyy_table: int=0, kyz_table: int=0, kzz_table: int=0,
            cp_table: int=0, hgen_table: int=0, comment: str='') -> int:
        """Creates a MATT5 card"""
        self.cards.append((mid, kxx_table, kxy_table, kxz_table, kyy_table,
                           kyz_table, kzz_table, cp_table, hgen_table,
                           comment))
        self.n += 1
        return self.n - 1

    def add_card(self, card: BDFCard, ifile: int, comment: str='') -> int:
        """
        Adds a MATT2 card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card

        """
        """
        Adds a MATT5 card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card

        """
        mid = integer(card, 1, 'mid')
        kxx_table = integer_or_blank(card, 2, 'T(Kxx)', default=0)
        kxy_table = integer_or_blank(card, 3, 'T(Kxy)', default=0)
        kxz_table = integer_or_blank(card, 5, 'T(Kxz)', default=0)
        kyy_table = integer_or_blank(card, 6, 'T(Kyy)', default=0)
        kyz_table = integer_or_blank(card, 7, 'T(Kyz)', default=0)
        kzz_table = integer_or_blank(card, 8, 'T(Kyz)', default=0)
        cp_table = integer_or_blank(card, 9, 'T(Kyz)', default=0)
        hgen_table = integer_or_blank(card, 11, 'T(HGEN)', default=0)

        assert len(card) <= 12, f'len(MATT5 card) = {len(card):d}\ncard={card}'
        #return MATT5(mid, kxx_table, kxy_table, kxz_table, kyy_table,
                     #kyz_table, kzz_table, cp_table, hgen_table,
                     #comment=comment)
        self.cards.append((mid, kxx_table, kxy_table, kxz_table, kyy_table,
                           kyz_table, kzz_table, cp_table, hgen_table,
                           comment))
        self.n += 1
        return self.n - 1

    @Material.parse_cards_check
    def parse_cards(self) -> None:
        ncards = len(self.cards)
        idtype = self.model.idtype
        material_id = np.zeros(ncards, dtype=idtype)
        kxx_table = np.zeros(ncards, dtype=idtype)
        kxy_table = np.zeros(ncards, dtype=idtype)
        kxz_table = np.zeros(ncards, dtype=idtype)
        kyy_table = np.zeros(ncards, dtype=idtype)
        kyz_table = np.zeros(ncards, dtype=idtype)
        kzz_table = np.zeros(ncards, dtype=idtype)
        cp_table = np.zeros(ncards, dtype=idtype)
        hgen_table = np.zeros(ncards, dtype=idtype)

        for icard, card in enumerate(self.cards):
            (mid, kxx_tablei, kxy_tablei, kxz_tablei, kyy_tablei,
             kyz_tablei, kzz_tablei, cp_tablei, hgen_tablei,
             comment) = card

            kxx_tablei = 0 if kxx_tablei is None else kxx_tablei
            kxy_tablei = 0 if kxy_tablei is None else kxy_tablei
            kxz_tablei = 0 if kxz_tablei is None else kxz_tablei

            kyy_tablei = 0 if kyy_tablei is None else kyy_tablei
            kyz_tablei = 0 if kyz_tablei is None else kyz_tablei
            kzz_tablei = 0 if kzz_tablei is None else kzz_tablei

            cp_tablei = 0 if cp_tablei is None else cp_tablei
            hgen_tablei = 0 if hgen_tablei is None else hgen_tablei
            material_id[icard] = mid
            kxx_table[icard] = kxx_tablei
            kxy_table[icard] = kxy_tablei
            kxz_table[icard] = kxz_tablei

            kyy_table[icard] = kyy_tablei
            kyz_table[icard] = kyz_tablei
            kzz_table[icard] = kzz_tablei

            cp_table[icard] = cp_tablei
            hgen_table[icard] = hgen_tablei
        self._save(material_id,
                   kxx_table, kxy_table, kxz_table, kyy_table,
                   kyz_table, kzz_table, cp_table, hgen_table)
        self.sort()
        self.cards = []

    def _save(self, material_id,
              kxx_table, kxy_table, kxz_table, kyy_table,
              kyz_table, kzz_table, cp_table, hgen_table):
        if len(self.material_id):
            raise RuntimeError(f'stacking of {self.type} is not supported')
            #material_id = np.hstack([self.material_id, material_id])
            #E = np.hstack([self.E, E])

        self.material_id = material_id
        self.kxx_table = kxx_table
        self.kxy_table = kxy_table
        self.kxz_table = kxz_table
        self.kyy_table = kyy_table
        self.kyz_table = kyz_table
        self.kzz_table = kzz_table
        self.cp_table = cp_table
        self.hgen_table = hgen_table

    def set_used(self, used_dict: dict[str, list[np.ndarray]]) -> None:
        table0 = np.hstack([
            self.kxx_table, self.kxy_table, self.kxz_table, self.kyy_table,
            kyz_table, self.kzz_table, self.cp_table, self.hgen_table, ])
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

    def __apply_slice__(self, mat: MATT5, i: np.ndarray) -> None:  # ignore[override]
        mat.n = len(i)
        mat.material_id = self.material_id[i]
        mat.kxx_table = self.kxx_table[i]
        mat.kxy_table = self.kxy_table[i]
        mat.kxz_table = self.kxz_table[i]
        mat.kyy_table = self.kyy_table[i]
        mat.kyz_table = self.kyz_table[i]
        mat.kzz_table = self.kzz_table[i]
        mat.cp_table = self.cp_table[i]
        mat.hgen_table = self.hgen_table[i]

    def geom_check(self, missing: dict[str, np.ndarray]):
        pass

    @property
    def max_id(self) -> int:
        tables = np.hstack([self.material_id,
                            self.kxx_table, self.kxy_table, self.kxz_table, self.kyy_table,
                            self.kyz_table, self.kzz_table, self.cp_table, self.hgen_table, ])
        return tables.max()

    @parse_check
    def write_file(self, bdf_file: TextIOLike,
                   size: int=8, is_double: bool=False,
                   write_card_header: bool=False) -> None:
        print_card, size = get_print_card_size(size, self.max_id)

        material_id = array_str(self.material_id, size=size)
        kxx_table = array_default_int(self.kxx_table, default=0, size=size)
        kxy_table = array_default_int(self.kxy_table, default=0, size=size)
        kxz_table = array_default_int(self.kxz_table, default=0, size=size)
        kyy_table = array_default_int(self.kyy_table, default=0, size=size)
        kyz_table = array_default_int(self.kyz_table, default=0, size=size)
        kzz_table = array_default_int(self.kzz_table, default=0, size=size)
        cp_table = array_default_int(self.cp_table, default=0, size=size)
        hgen_table = array_default_int(self.hgen_table, default=0, size=size)

        for (mid,
             kxx_tablei, kxy_tablei, kxz_tablei, kyy_tablei,
             kyz_tablei, kzz_tablei, cp_tablei, hgen_tablei) in zip_longest(
                 material_id,
                 kxx_table, kxy_table, kxz_table, kyy_table,
                 kyz_table, kzz_table, cp_table, hgen_table,):

            list_fields = [
                'MATT5', kxx_tablei, kxy_tablei, kxz_tablei,
                kyy_tablei, kyz_tablei, kzz_tablei,
                cp_tablei, None, hgen_tablei]
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
        return self.n - 1

    def add_card(self, card: BDFCard, ifile: int, comment: str='') -> int:
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
        return self.n - 1

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
        for icard, card in enumerate(self.cards):
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
            material_id[icard] = mid
            e1_table[icard] = e1_tablei
            e2_table[icard] = e2_tablei
            nu12_table[icard] = nu12_tablei
            g12_table[icard] = g12_tablei
            g1z_table[icard] = g1z_tablei
            g2z_table[icard] = g2z_tablei
            rho_table[icard] = rho_tablei
            a1_table[icard] = a1_tablei
            a2_table[icard] = a2_tablei
            xt_table[icard] = xt_tablei
            xc_table[icard] = yt_tablei
            yt_table[icard] = yt_tablei
            yc_table[icard] = yc_tablei
            s_table[icard] = s_tablei
            ge_table[icard] = ge_tablei
            f12_table[icard] = f12_tablei
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
            raise RuntimeError(f'stacking of {self.type} is not supported')
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

    #def convert(self, stiffness_scale: float=1.0, **kwargs) -> dict[int, np.ndarray]:
        #scale[('tablem', stiffness_scale)] = self.e1_table

    @property
    def max_id(self) -> int:
        tables = np.hstack([self.material_id, self.e1_table, self.e2_table, self.nu12_table,
                            self.g12_table, self.g1z_table, self.g2z_table, self.rho_table,
                            self.a1_table, self.a2_table,
                            self.xt_table, self.xc_table, self.yt_table, self.yc_table,
                            self.s_table, self.ge_table, self.f12_table, ],
        )
        return tables.max()

    @parse_check
    def write_file(self, bdf_file: TextIOLike,
                   size: int=8, is_double: bool=False,
                   write_card_header: bool=False) -> None:
        print_card, size = get_print_card_size(size, self.max_id)

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
    Solid Element Anisotropic Material Temperature Dependence
    Specifies temperature-dependent material properties on MAT9 entry
    fields via TABLEMi entries.

    +--------+---------+---------+---------+---------+---------+---------+---------+---------+
    |    1   |    2    |     3   |     4   |     5   |     6   |    7    |     8   |    9    |
    +========+=========+=========+=========+=========+=========+=========+=========+=========+
    | MATT9  |   MID   |  T(G11) |  T(G12) |  T(G13) |  T(G14) |  T(G15) |  T(G16) |  T(G22) |
    +--------+---------+---------+---------+---------+---------+---------+---------+---------+
    |        |  T(G23) |  T(G24) |  T(G25) |  T(G26) |  T(G33) |  T(G34) |  T(G35) |  T(G36) |
    +--------+---------+---------+---------+---------+---------+---------+---------+---------+
    |        |  T(G44) |  T(G45) |  T(G46) |  T(G55) |  T(G56) |  T(G66) |  T(RHO) |  T(A1)  |
    +--------+---------+---------+---------+---------+---------+---------+---------+---------+
    |        |  T(A2)  |  T(A3)  |  T(A4)  |  T(A5)  |  T(A6)  |         |   T(GE) |         |
    +--------+---------+---------+---------+---------+---------+---------+---------+---------+
    |        | T(GE11) | T(GE12) | T(GE13) | T(GE14) | T(GE15) | T(GE16) | T(GE22) | T(GE23) |
    +--------+---------+---------+---------+---------+---------+---------+---------+---------+
    |        | T(GE24) | T(GE25) | T(GE26) | T(GE33) | T(GE34) | T(GE35) | T(GE36) | T(GE44) |
    +--------+---------+---------+---------+---------+---------+---------+---------+---------+
    |        | T(GE45) | T(GE46) | T(GE55) | T(GE56) | T(GE66) |         |         |         |
    +--------+---------+---------+---------+---------+---------+---------+---------+---------+

    NX ends at T(GE)
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
            ge_table=None,
            ge11_table=None, ge12_table=None, ge13_table=None, ge14_table=None, ge15_table=None, ge16_table=None,
            ge22_table=None, ge23_table=None, ge24_table=None, ge25_table=None, ge26_table=None,
            ge33_table=None, ge34_table=None, ge35_table=None, ge36_table=None,
            ge44_table=None, ge45_table=None, ge46_table=None,
            ge55_table=None, ge56_table=None,
            ge66_table=None, comment: str='') -> int:
        """Creates a MATT9 card"""
        ges_tables = [
            ge11_table, ge12_table, ge13_table, ge14_table, ge15_table, ge16_table,
            ge22_table, ge23_table, ge24_table, ge25_table, ge26_table,
            ge33_table, ge34_table, ge35_table, ge36_table,
            ge44_table, ge45_table, ge46_table,
            ge55_table, ge56_table,
            ge66_table, ]
        self.cards.append((mid,
                           g11_table, g12_table, g13_table, g14_table, g15_table, g16_table,
                           g22_table, g23_table, g24_table, g25_table, g26_table,
                           g33_table, g34_table, g35_table, g36_table,
                           g44_table, g45_table, g46_table,
                           g55_table, g56_table, g66_table, rho_table,
                           a1_table, a2_table, a3_table,
                           a4_table, a5_table, a6_table, ge_table, ges_tables, comment))
        self.n += 1
        return self.n - 1

    def add_card(self, card: BDFCard, ifile: int, comment: str='') -> int:
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
        # end of NX

        #T(GE11) T(GE12) T(GE13) T(GE14) T(GE15) T(GE16) T(GE22) T(GE23)
        #T(GE24) T(GE25) T(GE26) T(GE33) T(GE34) T(GE35) T(GE36) T(GE44)
        #T(GE45) T(GE46) T(GE55) T(GE56) T(GE66)
        ge11_table = integer_or_blank(card, 33, 'T(GE11)', default=0)
        ge12_table = integer_or_blank(card, 34, 'T(GE12)', default=0)
        ge13_table = integer_or_blank(card, 35, 'T(GE13)', default=0)
        ge14_table = integer_or_blank(card, 36, 'T(GE14)', default=0)
        ge15_table = integer_or_blank(card, 37, 'T(GE15)', default=0)
        ge16_table = integer_or_blank(card, 38, 'T(GE16)', default=0)
        ge22_table = integer_or_blank(card, 39, 'T(GE22)', default=0)
        ge23_table = integer_or_blank(card, 40, 'T(GE23)', default=0)
        ge24_table = integer_or_blank(card, 41, 'T(GE24)', default=0)
        ge25_table = integer_or_blank(card, 42, 'T(GE25)', default=0)
        ge26_table = integer_or_blank(card, 43, 'T(GE26)', default=0)
        ge33_table = integer_or_blank(card, 44, 'T(GE33)', default=0)
        ge34_table = integer_or_blank(card, 45, 'T(GE34)', default=0)
        ge35_table = integer_or_blank(card, 46, 'T(GE35)', default=0)
        ge36_table = integer_or_blank(card, 47, 'T(GE36)', default=0)
        ge44_table = integer_or_blank(card, 48, 'T(GE44)', default=0)
        ge45_table = integer_or_blank(card, 49, 'T(GE45)', default=0)
        ge46_table = integer_or_blank(card, 50, 'T(GE46)', default=0)
        ge55_table = integer_or_blank(card, 51, 'T(GE55)', default=0)
        ge56_table = integer_or_blank(card, 52, 'T(GE56)', default=0)
        ge66_table = integer_or_blank(card, 53, 'T(GE66)', default=0)

        ges_tables = [
            ge11_table, ge12_table, ge13_table, ge14_table, ge15_table, ge16_table,
            ge22_table, ge23_table, ge24_table, ge25_table, ge26_table,
            ge33_table, ge34_table, ge35_table, ge36_table,
            ge44_table, ge45_table, ge46_table,
            ge55_table, ge56_table,
            ge66_table, ]

        #assert len(card) <= 32, f'len(MATT9 card) = {len(card):d}\ncard={card}'   # NX
        assert len(card) <= 54, f'len(MATT9 card) = {len(card):d}\ncard={card}'    # MSC
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
                           a4_table, a5_table, a6_table, ge_table, ges_tables, comment))
        self.n += 1
        return self.n - 1

    @Material.parse_cards_check
    def parse_cards(self) -> None:
        ncards = len(self.cards)
        idtype = self.model.idtype
        #fdtype = self.model.fdtype
        material_id = np.zeros(ncards, dtype=idtype)

        g11_table = np.zeros(ncards, dtype=idtype)
        g12_table = np.zeros(ncards, dtype=idtype)
        g13_table = np.zeros(ncards, dtype=idtype)
        g14_table = np.zeros(ncards, dtype=idtype)
        g15_table = np.zeros(ncards, dtype=idtype)
        g16_table = np.zeros(ncards, dtype=idtype)
        g22_table = np.zeros(ncards, dtype=idtype)
        g23_table = np.zeros(ncards, dtype=idtype)
        g24_table = np.zeros(ncards, dtype=idtype)
        g25_table = np.zeros(ncards, dtype=idtype)
        g26_table = np.zeros(ncards, dtype=idtype)
        g33_table = np.zeros(ncards, dtype=idtype)
        g34_table = np.zeros(ncards, dtype=idtype)
        g35_table = np.zeros(ncards, dtype=idtype)
        g36_table = np.zeros(ncards, dtype=idtype)
        g44_table = np.zeros(ncards, dtype=idtype)
        g45_table = np.zeros(ncards, dtype=idtype)
        g46_table = np.zeros(ncards, dtype=idtype)
        g55_table = np.zeros(ncards, dtype=idtype)
        g56_table = np.zeros(ncards, dtype=idtype)
        g66_table = np.zeros(ncards, dtype=idtype)
        rho_table = np.zeros(ncards, dtype=idtype)
        a1_table = np.zeros(ncards, dtype=idtype)
        a2_table = np.zeros(ncards, dtype=idtype)
        a3_table = np.zeros(ncards, dtype=idtype)
        a4_table = np.zeros(ncards, dtype=idtype)
        a5_table = np.zeros(ncards, dtype=idtype)
        a6_table = np.zeros(ncards, dtype=idtype)
        ge_table = np.zeros(ncards, dtype=idtype)
        ges_table = np.zeros((ncards, 21), dtype=idtype)
        for icard, card in enumerate(self.cards):
            (mid,
             g11_tablei, g12_tablei, g13_tablei, g14_tablei, g15_tablei, g16_tablei,
             g22_tablei, g23_tablei, g24_tablei, g25_tablei, g26_tablei,
             g33_tablei, g34_tablei, g35_tablei, g36_tablei,
             g44_tablei, g45_tablei, g46_tablei,
             g55_tablei, g56_tablei, g66_tablei,
             rho_tablei,
             a1_tablei, a2_tablei, a3_tablei,
             a4_tablei, a5_tablei, a6_tablei,
             ge_tablei, ges_tablei, comment) = card

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
            ges_tablei = [0 if ge_tablei is None else ge_tablei
                          for ge_tablei in ges_tablei]

            material_id[icard] = mid

            g11_table[icard] = g11_tablei
            g12_table[icard] = g12_tablei
            g13_table[icard] = g13_tablei
            g14_table[icard] = g14_tablei
            g15_table[icard] = g15_tablei
            g16_table[icard] = g16_tablei
            g22_table[icard] = g22_tablei
            g23_table[icard] = g23_tablei
            g24_table[icard] = g24_tablei
            g25_table[icard] = g25_tablei
            g26_table[icard] = g26_tablei
            g33_table[icard] = g33_tablei
            g34_table[icard] = g34_tablei
            g35_table[icard] = g35_tablei
            g36_table[icard] = g36_tablei
            g44_table[icard] = g44_tablei
            g45_table[icard] = g45_tablei
            g46_table[icard] = g46_tablei
            g55_table[icard] = g55_tablei
            g56_table[icard] = g56_tablei
            g66_table[icard] = g66_tablei
            rho_table[icard] = rho_tablei
            a1_table[icard] = a1_tablei
            a2_table[icard] = a2_tablei
            a3_table[icard] = a3_tablei
            a4_table[icard] = a4_tablei
            a5_table[icard] = a5_tablei
            a6_table[icard] = a6_tablei
            ge_table[icard] = ge_tablei
            ges_table[icard, :] = ges_tablei

        self._save(material_id,
                   g11_table, g12_table, g13_table, g14_table, g15_table, g16_table,
                   g22_table, g23_table, g24_table, g25_table, g26_table,
                   g33_table, g34_table, g35_table, g36_table,
                   g44_table, g45_table, g46_table,
                   g55_table, g56_table, g66_table, rho_table,
                   a1_table, a2_table, a3_table,
                   a4_table, a5_table, a6_table, ge_table, ges_table)
        self.sort()
        self.cards = []

    def _save(self, material_id,
              g11_table, g12_table, g13_table, g14_table, g15_table, g16_table,
              g22_table, g23_table, g24_table, g25_table, g26_table,
              g33_table, g34_table, g35_table, g36_table,
              g44_table, g45_table, g46_table,
              g55_table, g56_table, g66_table, rho_table,
              a1_table, a2_table, a3_table,
              a4_table, a5_table, a6_table, ge_table, ges_table):
        if len(self.material_id):
            raise RuntimeError(f'stacking of {self.type} is not supported')
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
        self.ges_table = ges_table

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
        mat.ges_table = self.ges_table[i, :]

    def geom_check(self, missing: dict[str, np.ndarray]):
        pass

    @property
    def max_id(self) -> int:
        tables = np.hstack([
            self.material_id,
            self.g11_table, self.g12_table, self.g13_table, self.g14_table, self.g15_table, self.g16_table,
            self.g22_table, self.g23_table, self.g24_table, self.g25_table, self.g26_table,
            self.g33_table, self.g34_table, self.g35_table, self.g36_table,
            self.g44_table, self.g45_table, self.g46_table,
            self.g55_table, self.g56_table, self.g66_table, self.rho_table,
            self.a1_table, self.a2_table, self.a3_table,
            self.a4_table, self.a5_table, self.a6_table,
            self.ge_table, self.ges_table.ravel(),
        ])
        return tables.max()

    @parse_check
    def write_file(self, bdf_file: TextIOLike,
                   size: int=8, is_double: bool=False,
                   write_card_header: bool=False) -> None:
        print_card, size = get_print_card_size(size, self.max_id)

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
        ges_table = array_default_int(self.ges_table, default=0, size=size).tolist()

        for (mid,
             g11_tablei, g12_tablei, g13_tablei, g14_tablei, g15_tablei, g16_tablei,
             g22_tablei, g23_tablei, g24_tablei, g25_tablei, g26_tablei,
             g33_tablei, g34_tablei, g35_tablei, g36_tablei,
             g44_tablei, g45_tablei, g46_tablei,
             g55_tablei, g56_tablei, g66_tablei, rho_tablei,
             a1_tablei, a2_tablei, a3_tablei, a4_tablei, a5_tablei, a6_tablei,
             ge_tablei, ges_tablei) in zip_longest(
                 material_id,
                 g11_table, g12_table, g13_table, g14_table, g15_table, g16_table,
                 g22_table, g23_table, g24_table, g25_table, g26_table,
                 g33_table, g34_table, g35_table, g36_table,
                 g44_table, g45_table, g46_table,
                 g55_table, g56_table, g66_table, rho_table,
                 a1_table, a2_table, a3_table, a4_table, a5_table, a6_table, ge_table,
                 ges_table):
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
            ] + ges_tablei

            bdf_file.write(print_card(list_fields))
        return
