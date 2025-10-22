"""
All material dependency cards are defined in this file.  This includes:

 * MATS1 (isotropic solid/shell)

 * MATT1 (isotropic solid/shell)
 * MATT2 (anisotropic)
 * MATT3 (linear orthotropic) - NA
 * MATT4 (thermal)
 * MATT5 (thermal)
 * MATT8 (orthotropic shell) - NA
 * MATT9 (anisotropic solid) - NA
 * MATT11 (orthotropic solid)

All cards are Material objects.

"""
#pylint: disable=E1103,C0103,C0111
from __future__ import annotations
from typing import Optional, TYPE_CHECKING

from pyNastran.bdf.cards.base_card import BaseCard
from pyNastran.bdf.bdf_interface.internal_get import material_id, table_id
from pyNastran.bdf.bdf_interface.assign_type import (
    integer, integer_or_blank, double, double_or_blank, string, string_or_blank)
from pyNastran.bdf.field_writer_8 import print_card_8
from pyNastran.bdf.field_writer_16 import print_card_16
if TYPE_CHECKING:  # pragma: no cover
    from pyNastran.bdf.bdf import BDF
    from pyNastran.bdf.bdf_interface.bdf_card import BDFCard


def _xref_table(self, model: BDF, key: str, msg: str) -> None:
    """
    Cross-references the TABLEM1s

    Parameters
    ----------
    key: str
        e_table
    """
    slot = getattr(self, key)  # etable
    if slot is not None and slot != 0:
        mid_ref = model.TableM(slot, msg + f' for {key}')
        setattr(self, key + '_ref', mid_ref)


class MaterialDependence(BaseCard):
    def __init__(self):
        self.mid = None

    def Mid(self) -> int:
        return material_id(self.mid_ref, self.mid)

    def _get_table(self, key: str) -> int:
        """internal method for accessing tables"""
        table = getattr(self, key)
        assert isinstance(table, int), key
        table_ref = getattr(self, key + '_ref')
        return table_id(table_ref, table)


class MaterialDependenceThermal(MaterialDependence):
    def __init__(self):
        MaterialDependence.__init__(self)

def _safe_xref_table(self, model: BDF, key: str,
                     xref_errors, msg: str):
    """
    Cross-references the TABLEM1s

    Parameters
    ----------
    key: str
        e_table
    """
    slot = getattr(self, key)
    if slot is not None and slot != 0 and slot in model.tables_m:
        mid_ref = model.safe_tablem(slot, self.mid, xref_errors, msg)
        setattr(self, key + '_ref', mid_ref)


class MATS1(MaterialDependence):
    """
    Specifies stress-dependent material properties for use in applications
    involving nonlinear materials. This entry is used if a MAT1, MAT2 or MAT9
    entry is specified with the same MID in a nonlinear solution sequence
    (SOLs 106 and 129).

    Format (NX Nastran):
    +--------+---------+-------+-------+------+-----+-----+--------+--------+
    |   1    |   2     |    3  |  4    |  5   |  6  |  7  |  8     |  9     |
    +========+=========+=======+=======+======+=====+=====+========+========+
    | MATS1  |  MID    | TID   | TYPE  |  H   | YF  | HR  | LIMIT1 | LIMIT2 |
    +--------+---------+-------+-------+------+-----+-----+--------+--------+
    |        | STRMEAS |       |       |      |     |     |        |        |
    +--------+---------+-------+-------+------+-----+-----+--------+--------+

    """
    type = 'MATS1'

    def __init__(self, mid: int, nl_type: Optional[str],
                 h: float, hr: float, yf: float,
                 limit1: Optional[float], limit2: Optional[float],
                 strmeas: Optional[str] = None,
                 tid: int=0, comment: str=''):
        MaterialDependence.__init__(self)
        if comment:
            self.comment = comment
        #: Identification number of a MAT1, MAT2, or MAT9 entry.
        self.mid = mid

        #: Identification number of a TABLES1 or TABLEST entry.
        # If H is given, then this field must be blank.
        self.tid = tid

        #: Type of material nonlinearity. ('NLELAST' for nonlinear elastic
        #: or 'PLASTIC' for elastoplastic.)
        if nl_type is None:
            pass
        else:
            assert isinstance(nl_type, str), nl_type
            nl_type = nl_type.upper()
            if nl_type == 'NLELAS':
                nl_type = 'NLELAST'
        self.nl_type = nl_type

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

        #: Stress/strain measure of the TABLES1 or TABLEST data referenced by the TID field.
        #: Valid for NX Nastran SOL 401 and SOL 402 only.
        self.strmeas = strmeas

        self.tid_ref = None
        self.mid_ref = None
        assert tid is not None

    @classmethod
    def _init_from_empty(cls):
        mid = 1
        tid = 1
        nl_type = None
        h = None
        hr = None
        yf = None
        limit1 = None
        limit2 = None
        strmeas = None
        return MATS1(mid, nl_type, h, hr, yf, limit1, limit2, strmeas, tid=tid, comment='')

    def validate(self) -> None:
        if self.nl_type not in ['NLELAST', 'PLASTIC', 'PLSTRN']:
            raise ValueError(f'MATS1 Type must be [NLELAST, PLASTIC, PLSTRN]; Type={self.nl_type!r}')

    @classmethod
    def add_card(cls, card: BDFCard, comment: str=''):
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
        tid = integer_or_blank(card, 2, 'tid', default=0)
        nl_type = string(card, 3, 'Type')

        if nl_type not in ['NLELAST', 'PLASTIC', 'PLSTRN']:
            raise ValueError('MATS1 Type must be [NLELAST, PLASTIC, PLSTRN]; Type=%r' % Type)
        if nl_type == 'NLELAST':
            # should we even read these?
            h = None
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
            h = double_or_blank(card, 4, 'H')
            yf = integer_or_blank(card, 5, 'yf', default=1)
            hr = integer_or_blank(card, 6, 'hr', default=1)
            limit1 = double_or_blank(card, 7, 'limit1')

            if yf in [3, 4]:
                limit2 = double(card, 8, 'limit2')
            else:
                #limit2 = blank(card, 8, 'limit2')
                limit2 = None

        if len(card) > 9:
            strmeas = string_or_blank(card, 9, 'strmeas')
        else:
            strmeas = None

        assert len(card) <= 10, f'len(MATS1 card) = {len(card):d}\ncard={card}'
        return MATS1(mid, nl_type, h, hr, yf, limit1, limit2, strmeas, tid=tid, comment=comment)

    @classmethod
    def add_op2_data(cls, data, comment: str=''):
        """
        Adds a MATS1 card from the OP2

        Parameters
        ----------
        data : list[varies]
            a list of fields defined in OP2 format
        comment : str; default=''
            a comment for the card

        """

        if len(data) < 9:
            (mid, tid, nl_type_int, h, yf, hr, limit1, limit2) = data
            strmeas = None
        else:
            (mid, tid, nl_type_int, h, yf, hr, limit1, limit2, strmeas_int) = data
            strmeas_map = {
                0: None,  # NULL
                1: 'UNDEF',
                2: 'ENG',
                3: 'TRUE',
                4: 'CAUCHY',
            }
            strmeas = strmeas_map[strmeas_int]

        if nl_type_int == 1:
            nl_type = 'NLELAST'
        elif nl_type_int == 2:
            nl_type = 'PLASTIC'
        elif nl_type_int == 3:
            nl_type = 'PLSTRN'
        else:  # pragma: no cover
            raise RuntimeError(f'Invalid Type:  mid={mid}; Type={nl_type_int}; must be 1=NLELAST, '
                               '2=PLASTIC, or 3=PLSTRN')

        return MATS1(mid, nl_type, h, hr, yf, limit1, limit2, strmeas, tid=tid, comment=comment)

    def Yf(self) -> str:
        d = {1: 'VonMises', 2: 'Tresca', 3: 'MohrCoulomb', 4: 'Drucker-Prager'}
        return d[self.yf]

    def Hf(self) -> str:
        d = {1: 'Isotropic', 2: 'Kinematic', 3: 'Combined'}
        return d[self.hr]

    def E(self, strain: float) -> None:
        """
        Gets E (Young's Modulus) for a given strain.

        Parameters
        ----------
        strain : float / None
            the strain (None -> linear E value)

        Returns
        -------
        E : float
            Young's Modulus

        """
        msg = "E (Young's Modulus) not implemented for MATS1"
        raise NotImplementedError(msg)
        #if self.tid:
            #E = self.tid_ref.Value(strain)
        #return E

    def cross_reference(self, model: BDF) -> None:
        """
        Cross links the card so referenced cards can be extracted directly

        Parameters
        ----------
        model : BDF()
            the BDF object

        """
        msg = ', which is required by MATS1 mid=%s' % self.mid
        self.mid_ref = model.Material(self.mid, msg=msg)
        if self.tid > 0:  # then self.h is used
            self.tid_ref = model.Table(self.tid, msg=msg) # TABLES1 or TABLEST

    def safe_cross_reference(self, model: BDF, xref_errors) -> None:
        """
        Cross links the card so referenced cards can be extracted directly

        Parameters
        ----------
        model : BDF()
            the BDF object

        """
        msg = ', which is required by MATS1 mid=%s' % self.mid
        self.mid_ref = model.safe_material(self.mid, self.mid, xref_errors, msg=msg)
        if self.tid > 0:  # then self.h is used
            self.tid_ref = model.Table(self.tid, msg=msg) # TABLES1 or TABLEST

    def uncross_reference(self) -> None:
        """Removes cross-reference links"""
        self.mid = self.Mid()
        if self.tid:
            self.tid = self.Tid()
        self.tid_ref = None
        self.mid_ref = None

    def Tid(self) -> int:
        return table_id(self.tid_ref, self.tid)

    def raw_fields(self) -> list:
        list_fields = ['MATS1', self.Mid(), self.Tid(), self.nl_type,
                       self.h, self.yf, self.hr, self.limit1, self.limit2, self.strmeas]
        return list_fields

    def repr_fields(self) -> list:
        return self.raw_fields()

    def write_card(self, size: int=8, is_double: bool=False) -> str:
        card = self.repr_fields()
        if size == 8:
            return self.comment + print_card_8(card)
        return self.comment + print_card_16(card)


class MATDMG(MaterialDependence):
    """
    Material Properties for Progressive Ply Failure

    Defines material properties and parameters for progressive ply failure in composite solid elements defined with the
    PCOMPS bulk entry. Used in combination with MAT11 entries that have the same MID. Valid for SOLs 401 and 402.

    +--------+-------+--------+-------+------+-----+-----+-----+-----+
    |   1    |   2   |    3   |  4    |  5   |  6  |  7  |  8  |  9  |
    +========+=======+========+=======+======+=====+=====+=====+=====+
    | MATDMG |  MID  | PPFMOD |       |      |     |     |     |     |
    +--------+-------+--------+-------+------+-----+-----+-----+-----+
    |        | COEF1 | COEF2  | -etc- |      |     |     |     |     |
    +--------+-------+--------+-------+------+-----+-----+-----+-----+

    Format for PPFMOD == "UD"
    +--------+--------+------+------+------+---------+---------+------+------+
    |   1    |   2    |   3  |  4   |  5   |    6    |    7    |  8   |  9   |
    +========+========+======+======+======+=========+=========+======+======+
    | MATDMG |  MID   |  UD  |      |      |         |         |      |      |
    +--------+--------+------+------+------+---------+---------+------+------+
    |        |  Y012  | YC12 | YS12 | YS22 | Y11LIMT | Y11LIMC | KSIT | KSIC |
    +--------+--------+------+------+------+---------+---------+------+------+
    |        |   B2   |  B3  |  A   | LITK |   BIGK  |   EXPN  | TAU | ADEL  |
    +--------+--------+------+------+------+---------+---------+------+------+
    |        | PLYUNI | TID  | HBAR | DMAX |   PE    |         |     |       |
    +--------+--------+------+------+------+---------+---------+------+------+

    Format for PPFMOD == "EUD"
    +--------+--------+------+------+-------+---------+---------+------+------+
    |   1    |   2    |   3  |  4   |  5    |    6    |    7    |  8   |  9   |
    +========+========+======+======+=======+=========+=========+======+======+
    | MATDMG |   MID  | EUD  |      |       |         |         |      |      |
    +--------+--------+------+------+-------+---------+---------+------+------+
    |        |  Y012  | YC12 |   K  | ALPHA | Y11LIMT | Y11LIMC | KSIT | KSIC |
    +--------+--------+------+------+-------+---------+---------+------+------+
    |        |   B2   |  B3  |   A  | LITK  |  BIGK   |  EXPN   | TAU  | ADEL |
    +--------+--------+------+------+-------+---------+---------+------+------+
    |        |  USER  | R01  | HBAR | DMAX  |   DS    |   GIC   | GIIC | GIIIC|
    +--------+--------+------+------+-------+---------+---------+------+------+

    """
    type = 'MATDMG'

    def __init__(self, mid: int, ppf_model: str,
                 y012: float, yc12: float, ys12: float,
                 ys22: float, y11limt: float, y11limc: float,
                 ksit=None, ksic=None,
                 b2=None, b3=None, a=None,
                 litk=None, bigk=None, expn=None,
                 tau=None, adel=None,
                 plyuni=None, tid=None, hbar=None, dmax=None, pe=None,
                 user=None, r01=None, ds=None,
                 gic=None, giic=None, giiic=None,
                 comment: str=''):
        MaterialDependence.__init__(self)
        if comment:
            self.comment = comment
        if ppf_model == 'EUD':
            assert gic is not None, f'Invalid GIC; {gic}; should be a float in case of PPFMOD == "EUD"'
            assert giic is not None, f'Invalid GIIC; {giic}; should be a float in case of PPFMOD == "EUD"'
            assert giiic is not None, f'Invalid GIIIC; {giiic}; should be a float in case of PPFMOD == "EUD"'

        self.mid = mid
        self.ppf_model = ppf_model

        self.y012 = y012
        self.yc12 = yc12
        self.ys12 = ys12
        self.ys22 = ys22
        self.y11limt = y11limt
        self.y11limc = y11limc
        self.ksit = ksit
        self.ksic = ksic
        self.b2 = b2
        self.b3 = b3
        self.a = a
        self.litk = litk
        self.bigk = bigk
        self.expn = expn
        self.tau = tau
        self.adel = adel
        self.plyuni = plyuni
        self.tid = tid
        self.hbar = hbar
        self.dmax = dmax
        self.pe = pe
        self.user = user
        self.r01 = r01
        self.ds = ds
        self.gic = gic
        self.giic = giic
        self.giiic = giiic

    @classmethod
    def add_card(cls, card: BDFCard, comment: str=''):
        """
        Adds a MATDMG card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card

        """
        mid = integer(card, 1, 'mid')
        ppf_model = string(card, 2, 'PPFMOD')

        y012 = double(card, 9, 'Y012')
        yc12 = double(card, 10, 'YC12')
        ys12 = double(card, 11, 'YS12')
        ys22 = double(card, 12, 'YS22')

        y11limt = double(card, 13, 'Y11LIMT')
        y11limc = double(card, 14, 'Y11LIMC')
        ksit = double_or_blank(card, 15, 'KSIT')
        ksic = double_or_blank(card, 16, 'KSIC')

        b2 = double(card, 17, 'B2')
        b3 = double(card, 18, 'B3')
        a = double(card, 19, 'A')

        litk = double(card, 20, 'LITK')
        bigk = double(card, 21, 'BIGK')
        expn = double(card, 22, 'EXPN')

        tau = double(card, 23, 'TAU')
        adel = double(card, 24, 'ADEL')

        if ppf_model == 'UD':
            plyuni = integer_or_blank(card, 25, 'PLYUNI')
            tid = integer_or_blank(card, 26, 'TID')
            hbar = double(card, 27, 'HBAR')
            dmax = double(card, 28, 'DMAX')
            pe = integer_or_blank(card, 29, 'PE')

            user = None
            r01 = None
            ds = None
            gic = None
            giic = None
            giiic = None
        elif ppf_model == 'EUD':
            user = integer(card, 25, 'USER')
            r01 = double(card, 26, 'R01')
            hbar = double(card, 27, 'HBAR')
            dmax = double(card, 28, 'DMAX')
            ds = double_or_blank(card, 29, 'DS')
            gic = double(card, 30, 'GIC')
            giic = double(card, 31, 'GIIC')
            giiic = double(card, 32, 'GIIIC')

            plyuni = None
            tid = None
            pe = None
        else:
            raise RuntimeError(f'Invalid PPFMOD: mid={mid}; PPFMOD={ppf_model}; must be UD or EUD.')

        return MATDMG(mid, ppf_model,
                      y012, yc12, ys12, ys22, y11limt, y11limc, ksit, ksic,
                      b2, b3, a, litk, bigk, expn, tau, adel,
                      plyuni, tid, hbar, dmax, pe,
                      user, r01, ds, gic, giic, giiic,
                      comment=comment)

    def raw_fields(self):

        if self.ppf_model == "UD":
            list_fields = ['MATDMG', self.mid, self.ppf_model, None, None, None, None, None, None,
                           self.y012, self.yc12, self.ys12, self.ys22, self.y11limt, self.y11limc, self.ksit, self.ksic,
                           self.b2, self.b3, self.a, self.litk, self.bigk, self.expn, self.tau, self.adel,
                           self.plyuni, self.tid, self.hbar, self.dmax, self.pe,
                           self.user, self.r01, self.ds, self.gic, self.giic, self.giiic]
        else:
            list_fields = ['MATDMG', self.mid, self.ppf_model, None, None, None, None, None, None,
                           self.y012, self.yc12, self.ys12, self.ys22, self.y11limt, self.y11limc, self.ksit, self.ksic,
                           self.b2, self.b3, self.a, self.litk, self.bigk, self.expn, self.tau, self.adel,
                           self.user, self.r01, self.hbar, self.dmax, self.ds,
                           self.gic, self.giic, self.giiic]

        return list_fields

    def repr_fields(self):
        """
        Gets the fields in their simplified form

        Returns
        -------
        fields : [varies, ...]
            the fields that define the card

        """
        return self.raw_fields()

    def write_card(self, size: int=8, is_double: bool=False) -> str:
        card = self.repr_fields()
        if size == 8:
            return self.comment + print_card_8(card)
        return self.comment + print_card_16(card)

class MATT1(MaterialDependenceThermal):
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
    type = 'MATT1'

    @classmethod
    def _init_from_empty(cls):
        mid = 1
        return MATT1(mid, e_table=None)

    def __init__(self, mid: int,
                 e_table: int=0, g_table: int=0, nu_table: int=0,
                 rho_table: int=0, a_table: int=0, ge_table: int=0,
                 st_table: int=0, sc_table: int=0, ss_table: int=0,
                 comment: str=''):
        MaterialDependenceThermal.__init__(self)
        if comment:
            self.comment = comment
        self.mid = mid
        assert e_table is not None, e_table
        # if e_table == 0:
        #     e_table = None
        # if g_table == 0:
        #     g_table = None
        # if nu_table == 0:
        #     nu_table = None
        # if rho_table == 0:
        #     rho_table = None
        # if a_table == 0:
        #     a_table = None
        # if ge_table == 0:
        #     ge_table = None
        # if st_table == 0:
        #     st_table = None
        # if sc_table == 0:
        #     sc_table = None
        # if ss_table == 0:
        #     ss_table = None
        self.e_table = e_table
        self.g_table = g_table
        self.nu_table = nu_table
        self.rho_table = rho_table
        self.a_table = a_table
        self.ge_table = ge_table
        self.st_table = st_table
        self.sc_table = sc_table
        self.ss_table = ss_table

        self.mid_ref = None
        self.ss_table_ref = None
        self.e_table_ref = None
        self.g_table_ref = None
        self.nu_table_ref = None
        self.rho_table_ref = None
        self.a_table_ref = None
        self.ge_table_ref = None
        self.st_table_ref = None
        self.sc_table_ref = None
        self.ss_table_ref = None

    @classmethod
    def add_card(cls, card: BDFCard, comment: str=''):
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
        e_table = integer_or_blank(card, 2, 'T(E)', default=0)
        g_table = integer_or_blank(card, 3, 'T(G)', default=0)
        nu_table = integer_or_blank(card, 4, 'T(nu)', default=0)
        rho_table = integer_or_blank(card, 5, 'T(rho)', default=0)
        a_table = integer_or_blank(card, 6, 'T(A)', default=0)
        ge_table = integer_or_blank(card, 8, 'T(ge)', default=0)
        st_table = integer_or_blank(card, 9, 'T(st)', default=0)
        sc_table = integer_or_blank(card, 10, 'T(sc)', default=0)
        ss_table = integer_or_blank(card, 11, 'T(ss)', default=0)

        assert len(card) <= 12, f'len(MATT1 card) = {len(card):d}\ncard={card}'
        return MATT1(mid, e_table, g_table, nu_table, rho_table, a_table,
                     ge_table, st_table, sc_table, ss_table, comment=comment)

    @classmethod
    def add_op2_data(cls, data: list[int], comment: str=''):
        """
        Adds a MATT1 card from the OP2

        Parameters
        ----------
        data : list[varies]
            a list of fields defined in OP2 format
        comment : str; default=''
            a comment for the card

        """
        (mid, E_table, G_table, nu_table, rho_table, A_table, dunno_a, ge_table,
         st_table, sc_table, ss_table, dunno_b) = data
        # if E_table == 0:
        #     E_table = None
        if E_table > 100000000:
            E_table = -(E_table - 100000000)

        # if G_table == 0:
        #     G_table = None
        if G_table > 100000000:
            G_table = -(G_table - 100000000)

        # if nu_table == 0:
        #     nu_table = None
        if nu_table > 100000000:
            nu_table = -(nu_table - 100000000)

        # if rho_table == 0:
        #     rho_table = None
        if rho_table > 100000000:
            rho_table = -(rho_table - 100000000)

        # if A_table == 0:
        #     A_table = None
        if A_table > 100000000:
            A_table = -(A_table - 100000000)

        mat = MATT1(mid, E_table, G_table, nu_table, rho_table, A_table,
                    ge_table, st_table, sc_table, ss_table, comment=comment)
        assert dunno_a == 0, '%s; dunno_a=%s\n%s' % (data, dunno_a, str(mat))
        assert dunno_b == 0, '%s; dunno_b=%s\n%s' % (data, dunno_b, str(mat))
        return mat

    def E(self, temperature: float) -> float:
        """
        Gets E (Young's Modulus) for a given temperature.

        Parameters
        ----------
        temperature : float; default=None
            the temperature (None -> linear E value)

        Returns
        -------
        E : float
            Young's Modulus

        """
        E = None
        if self.E_table:
            E = self.E_table.Value(temperature)
        return E

    def cross_reference(self, model: BDF) -> None:
        """
        Cross links the card so referenced cards can be extracted directly

        Parameters
        ----------
        model : BDF()
            the BDF object

        """
        msg = ', which is required by MATT1 mid=%s' % self.mid
        self.mid_ref = model.Material(self.mid, msg=msg)

        _xref_table(self, model, 'e_table', msg=msg)
        _xref_table(self, model, 'g_table', msg=msg)
        _xref_table(self, model, 'nu_table', msg=msg)
        _xref_table(self, model, 'rho_table', msg=msg)
        _xref_table(self, model, 'a_table', msg=msg)
        _xref_table(self, model, 'ge_table', msg=msg)
        _xref_table(self, model, 'st_table', msg=msg)
        _xref_table(self, model, 'sc_table', msg=msg)
        _xref_table(self, model, 'ss_table', msg=msg)

    def safe_cross_reference(self, model: BDF, xref_errors) -> None:
        """
        Cross links the card so referenced cards can be extracted directly

        Parameters
        ----------
        model : BDF()
            the BDF object

        """
        msg = ', which is required by MATT1 mid=%s' % self.mid
        self.mid_ref = model.safe_material(self.mid, self.mid, xref_errors, msg=msg)

        _safe_xref_table(self, model, 'e_table', xref_errors, msg=msg)
        _safe_xref_table(self, model, 'g_table', xref_errors, msg=msg)
        _safe_xref_table(self, model, 'nu_table', xref_errors, msg=msg)
        _safe_xref_table(self, model, 'rho_table', xref_errors, msg=msg)
        _safe_xref_table(self, model, 'a_table', xref_errors, msg=msg)
        _safe_xref_table(self, model, 'ge_table', xref_errors, msg=msg)
        _safe_xref_table(self, model, 'st_table', xref_errors, msg=msg)
        _safe_xref_table(self, model, 'sc_table', xref_errors, msg=msg)
        _safe_xref_table(self, model, 'ss_table', xref_errors, msg=msg)

    def uncross_reference(self) -> None:
        """Removes cross-reference links"""
        self.mid = self.Mid()
        ## TODO: remove refs

        self.e_table = self.E_table()
        self.g_table = self.G_table()
        self.nu_table = self.Nu_table()
        self.rho_table = self.Rho_table()
        self.a_table = self.A_table()
        self.ge_table = self.Ge_table()
        self.st_table = self.St_table()
        self.sc_table = self.Sc_table()
        self.ss_table = self.Ss_table()

        self.mid_ref = None
        self.mid_ref = None
        self.ss_table_ref = None
        self.e_table_ref = None
        self.g_table_ref = None
        self.nu_table_ref = None
        self.rho_table_ref = None
        self.a_table_ref = None
        self.ge_table_ref = None
        self.st_table_ref = None
        self.sc_table_ref = None
        self.ss_table_ref = None

    def E_table(self) -> int:
        return self._get_table('e_table')

    def G_table(self) -> int:
        return self._get_table('g_table')

    def Nu_table(self) -> int:
        return self._get_table('nu_table')

    def Rho_table(self) -> int:
        return self._get_table('rho_table')

    def A_table(self) -> int:
        return self._get_table('a_table')

    def Ge_table(self) -> int:
        return self._get_table('ge_table')

    def St_table(self) -> int:
        return self._get_table('st_table')

    def Sc_table(self) -> int:
        return self._get_table('sc_table')

    def Ss_table(self) -> int:
        return self._get_table('ss_table')

    def raw_fields(self) -> list:
        list_fields = [
            'MATT1', self.Mid(), self.E_table(), self.G_table(),
            self.Nu_table(), self.Rho_table(), self.A_table(), self.Ge_table(),
            self.St_table(), self.Sc_table(), self.Ss_table(),
        ]
        return list_fields

    def repr_fields(self):
        return self.raw_fields()

    def write_card(self, size: int=8, is_double: bool=False) -> str:
        card = self.repr_fields()
        if size == 8:
            return self.comment + print_card_8(card)
        return self.comment + print_card_16(card)

class MATT2(MaterialDependenceThermal):
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
    type = 'MATT2'

    @classmethod
    def _init_from_empty(cls):
        mid = 1
        return MATT2(mid, g11_table=None)

    def __init__(self, mid: int, g11_table: int=0, g12_table: int=0, g13_table: int=0,
                 g22_table: int=0, g23_table: int=0, g33_table: int=0, rho_table: int=0,
                 a1_table: int=0, a2_table: int=0, a3_table: int=0,
                 ge_table: int=0, st_table: int=0, sc_table: int=0, ss_table: int=0,
                 comment: str=''):
        MaterialDependenceThermal.__init__(self)
        if comment:
            self.comment = comment

        # if g11_table == 0:
        #     g11_table = None
        # if g12_table == 0:
        #     g12_table = None
        # if g13_table == 0:
        #     g13_table = None
        # if g22_table == 0:
        #     g22_table = None
        # if g23_table == 0:
        #     g23_table = None
        # if g33_table == 0:
        #     g33_table = None
        # if rho_table == 0:
        #     rho_table = None

        # if a1_table == 0:
        #     a1_table = None
        # if a2_table == 0:
        #     a2_table = None
        # if a3_table == 0:
        #     a3_table = None
        # if ge_table == 0:
        #     ge_table = None

        # if st_table == 0:
        #     st_table = None
        # if sc_table == 0:
        #     sc_table = None
        # if ss_table == 0:
        #     ss_table = None
        self.mid = mid
        self.g11_table = g11_table
        self.g12_table = g12_table
        self.g13_table = g13_table
        self.g22_table = g22_table
        self.g23_table = g23_table
        self.g33_table = g33_table
        self.rho_table = rho_table
        self.a1_table = a1_table
        self.a2_table = a2_table
        self.a3_table = a3_table
        self.ge_table = ge_table
        self.st_table = st_table
        self.sc_table = sc_table
        self.ss_table = ss_table

        self.mid_ref = None
        self.g11_table_ref = None
        self.g12_table_ref = None
        self.g13_table_ref = None
        self.g22_table_ref = None
        self.g23_table_ref = None
        self.g33_table_ref = None
        self.rho_table_ref = None
        self.a1_table_ref = None
        self.a2_table_ref = None
        self.a3_table_ref = None
        self.ge_table_ref = None
        self.st_table_ref = None
        self.sc_table_ref = None
        self.ss_table_ref = None

    def validate(self):
        assert self.g11_table >= 0, self.g11_table
        assert self.g12_table >= 0, self.g12_table
        assert self.g13_table >= 0, self.g13_table
        assert self.g22_table >= 0, self.g22_table
        assert self.g23_table >= 0, self.g23_table
        assert self.g33_table >= 0, self.g33_table
        assert self.rho_table >= 0, self.rho_table
        assert self.a1_table >= 0, self.a1_table
        assert self.a2_table >= 0, self.a2_table
        assert self.a3_table >= 0, self.a3_table
        assert self.ge_table >= 0, self.ge_table
        assert self.st_table >= 0, self.st_table
        assert self.sc_table >= 0, self.sc_table
        assert self.ss_table >= 0, self.ss_table

    @classmethod
    def add_card(cls, card: BDFCard, comment: str=''):
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
        return MATT2(mid, g11_table, g12_table, g13_table, g22_table, g23_table,
                     g33_table, rho_table, a1_table,
                     a2_table, a3_table, ge_table,
                     st_table, sc_table, ss_table,
                     comment=comment)

    def cross_reference(self, model: BDF) -> None:
        """
        Cross links the card so referenced cards can be extracted directly

        Parameters
        ----------
        model : BDF()
            the BDF object

        """
        msg = ', which is required by MATT2 mid=%s' % self.mid
        self.mid_ref = model.Material(self.mid, msg=msg)

        _xref_table(self, model, 'g11_table', msg=msg)
        _xref_table(self, model, 'g12_table', msg=msg)
        _xref_table(self, model, 'g13_table', msg=msg)
        _xref_table(self, model, 'g22_table', msg=msg)
        _xref_table(self, model, 'g23_table', msg=msg)
        _xref_table(self, model, 'g33_table', msg=msg)
        _xref_table(self, model, 'rho_table', msg=msg)
        _xref_table(self, model, 'a1_table', msg=msg)
        _xref_table(self, model, 'a2_table', msg=msg)
        _xref_table(self, model, 'a3_table', msg=msg)
        _xref_table(self, model, 'ge_table', msg=msg)
        _xref_table(self, model, 'st_table', msg=msg)
        _xref_table(self, model, 'sc_table', msg=msg)
        _xref_table(self, model, 'ss_table', msg=msg)

    def safe_cross_reference(self, model: BDF, xref_errors) -> None:
        """
        Cross links the card so referenced cards can be extracted directly

        Parameters
        ----------
        model : BDF()
            the BDF object

        """
        msg = ', which is required by MATT2 mid=%s' % self.mid
        self.mid_ref = model.safe_material(self.mid, self.mid, xref_errors, msg=msg)

        _safe_xref_table(self, model, 'g11_table', xref_errors, msg=msg)
        _safe_xref_table(self, model, 'g12_table', xref_errors, msg=msg)
        _safe_xref_table(self, model, 'g13_table', xref_errors, msg=msg)
        _safe_xref_table(self, model, 'g22_table', xref_errors, msg=msg)
        _safe_xref_table(self, model, 'g23_table', xref_errors, msg=msg)
        _safe_xref_table(self, model, 'g33_table', xref_errors, msg=msg)
        _safe_xref_table(self, model, 'rho_table', xref_errors, msg=msg)
        _safe_xref_table(self, model, 'a1_table', xref_errors, msg=msg)
        _safe_xref_table(self, model, 'a2_table', xref_errors, msg=msg)
        _safe_xref_table(self, model, 'a3_table', xref_errors, msg=msg)
        _safe_xref_table(self, model, 'ge_table', xref_errors, msg=msg)
        _safe_xref_table(self, model, 'st_table', xref_errors, msg=msg)
        _safe_xref_table(self, model, 'sc_table', xref_errors, msg=msg)
        _safe_xref_table(self, model, 'ss_table', xref_errors, msg=msg)

    def uncross_reference(self) -> None:
        """Removes cross-reference links"""
        self.mid = self.Mid()
        self.g11_table = self.G11_table()
        self.g12_table = self.G12_table()
        self.g13_table = self.G13_table()
        self.g22_table = self.G22_table()
        self.g23_table = self.G23_table()
        self.g33_table = self.G33_table()
        self.rho_table = self.Rho_table()
        self.a1_table = self.A1_table()
        self.a2_table = self.A2_table()
        self.a3_table = self.A3_table()
        self.ge_table = self.Ge_table()
        self.st_table = self.St_table()
        self.sc_table = self.Sc_table()
        self.ss_table = self.Ss_table()
        self.mid_ref = None

    def G11_table(self) -> int:
        return self._get_table('g11_table')

    def G12_table(self) -> int:
        return self._get_table('g12_table')

    def G13_table(self) -> int:
        return self._get_table('g13_table')

    def G22_table(self) -> int:
        return self._get_table('g22_table')

    def G23_table(self) -> int:
        return self._get_table('g23_table')

    def G33_table(self) -> int:
        return self._get_table('g33_table')

    def Rho_table(self) -> int:
        return self._get_table('rho_table')

    def A1_table(self) -> int:
        return self._get_table('a1_table')

    def A2_table(self) -> int:
        return self._get_table('a2_table')

    def A3_table(self) -> int:
        return self._get_table('a3_table')

    def Ge_table(self) -> int:
        return self._get_table('ge_table')

    def St_table(self) -> int:
        return self._get_table('st_table')

    def Sc_table(self) -> int:
        return self._get_table('sc_table')

    def Ss_table(self) -> int:
        return self._get_table('ss_table')

    def raw_fields(self) -> list:
        list_fields = [
            'MATT2', self.Mid(), self.G11_table(), self.G12_table(),
            self.G13_table(), self.G22_table(), self.G23_table(),
            self.G33_table(), self.Rho_table(), self.A1_table(),
            self.A2_table(), self.A3_table(), None, self.Ge_table(),
            self.St_table(), self.Sc_table(), self.Ss_table()
        ]
        return list_fields

    def repr_fields(self):
        return self.raw_fields()

    def write_card(self, size: int=8, is_double: bool=False) -> str:
        card = self.repr_fields()
        if size == 8:
            return self.comment + print_card_8(card)
        return self.comment + print_card_16(card)

#MATT3 - CTRIAX6 only
class MATT3(MaterialDependenceThermal):
    """
    Specifies temperature-dependent material properties on MAT3 entry fields via
    TABLEMi entries that are temperature dependent.

    +--------+-------+-------+--------+-------+----------+----------+---------+--------+
    |    1   |   2   |   3   |    4   |   5   |     6    |     7    |    8    |    9   |
    +========+=======+=======+========+=======+==========+==========+=========+========+
    | MATT3  |  MID  | T(EX) | T(ETH) | T(EZ) | T(NUXTH) | T(NUTHZ) | T(NUZX) | T(RHO) |
    +--------+-------+-------+--------+-------+----------+----------+---------+--------+
    |        |       |       | T(GZX) | T(AX) |  T(ATH)  |  T(AZ)   |         |  T(GE) |
    +--------+-------+-------+--------+-------+----------+----------+---------+--------+

    """
    type = 'MATT3'

    @classmethod
    def _init_from_empty(cls):
        mid = 1
        return MATT3(mid, ex_table=None)

    def __init__(self, mid: int, ex_table: int=0, eth_table: int=0, ez_table: int=0,
                 nuth_table: int=0, nuxz_table: int=0, rho_table: int=0,
                 gzx_table: int=0, ax_table: int=0, ath_table: int=0, az_table: int=0,
                 ge_table: int=0, comment: str=''):
        """
        Creates a MATT3 card
        """
        MaterialDependenceThermal.__init__(self)
        if comment:
            self.comment = comment

        # ex_table = None if ex_table == 0 else ex_table
        # eth_table = None if eth_table == 0 else eth_table
        # ez_table = None if ez_table == 0 else ez_table
        # nuth_table = None if nuth_table == 0 else nuth_table
        # nuxz_table = None if nuxz_table == 0 else nuxz_table
        # rho_table = None if rho_table == 0 else rho_table
        # gzx_table = None if gzx_table == 0 else gzx_table
        # ax_table = None if ax_table == 0 else ax_table
        # ath_table = None if ath_table == 0 else ath_table
        # az_table = None if az_table == 0 else az_table
        # ge_table = None if ge_table == 0 else ge_table

        self.mid = mid
        self.ex_table = ex_table
        self.eth_table = eth_table
        self.ez_table = ez_table
        self.nuth_table = nuth_table
        self.nuxz_table = nuxz_table
        self.rho_table = rho_table
        self.gzx_table = gzx_table
        self.ax_table = ax_table
        self.ath_table = ath_table
        self.az_table = az_table
        self.ge_table = ge_table

        self.ex_table_ref = None
        self.eth_table_ref = None
        self.ez_table_ref = None
        self.nuth_table_ref = None
        self.nuxz_table_ref = None
        self.rho_table_ref = None
        self.gzx_table_ref = None
        self.ax_table_ref = None
        self.ath_table_ref = None
        self.az_table_ref = None
        self.ge_table_ref = None
        self.mid_ref = None

    def validate(self) -> None:
        assert self.ex_table >= 0, self.ex_table
        assert self.eth_table >= 0, self.eth_table
        assert self.ez_table >= 0, self.ez_table
        assert self.nuth_table >= 0, self.nuth_table
        assert self.nuxz_table >= 0, self.nuxz_table
        assert self.rho_table >= 0, self.rho_table
        assert self.gzx_table >= 0, self.gzx_table
        assert self.ax_table >= 0, self.ax_table
        assert self.ath_table >= 0, self.ath_table
        assert self.az_table >= 0, self.az_table
        assert self.ge_table >= 0, self.ge_table

    def cross_reference(self, model: BDF) -> None:
        msg = ', which is required by MATT3 mid=%s' % self.mid
        self.mid_ref = model.Material(self.mid, msg=msg)

        #self._get_table('ex_table')
        _xref_table(self, model, 'ex_table', msg=msg)
        _xref_table(self, model, 'eth_table', msg=msg)
        _xref_table(self, model, 'ez_table', msg=msg)
        _xref_table(self, model, 'nuth_table', msg=msg)
        _xref_table(self, model, 'nuxz_table', msg=msg)
        _xref_table(self, model, 'rho_table', msg=msg)

        _xref_table(self, model, 'gzx_table', msg=msg)
        _xref_table(self, model, 'ax_table', msg=msg)
        _xref_table(self, model, 'ath_table', msg=msg)
        _xref_table(self, model, 'az_table', msg=msg)
        _xref_table(self, model, 'ge_table', msg=msg)

    def safe_cross_reference(self, model: BDF, xref_errors) -> None:
        msg = ', which is required by MATT3 mid=%s' % self.mid
        self.mid_ref = model.safe_material(self.mid, self.mid, xref_errors, msg=msg)

        #self._get_table('ex_table')
        _safe_xref_table(self, model, 'ex_table', xref_errors, msg=msg)
        _safe_xref_table(self, model, 'eth_table', xref_errors, msg=msg)
        _safe_xref_table(self, model, 'ez_table', xref_errors, msg=msg)
        _safe_xref_table(self, model, 'nuth_table', xref_errors, msg=msg)
        _safe_xref_table(self, model, 'nuxz_table', xref_errors, msg=msg)
        _safe_xref_table(self, model, 'rho_table', xref_errors, msg=msg)

        _safe_xref_table(self, model, 'gzx_table', xref_errors, msg=msg)
        _safe_xref_table(self, model, 'ax_table', xref_errors, msg=msg)
        _safe_xref_table(self, model, 'ath_table', xref_errors, msg=msg)
        _safe_xref_table(self, model, 'az_table', xref_errors, msg=msg)
        _safe_xref_table(self, model, 'ge_table', xref_errors, msg=msg)

    def uncross_reference(self) -> None:
        """Removes cross-reference links"""
        self.mid = self.Mid()
        self.mid_ref = None

        self.ex_table = self.Ex_table()
        self.eth_table = self.Eth_table()
        self.ez_table = self.Ez_table()
        self.nuth_table = self.Nuth_table()
        self.nuxz_table = self.Nuxz_table()
        self.rho_table = self.Rho_table()
        self.gzx_table = self.Gzx_table()
        self.ax_table = self.Ax_table()
        self.ath_table = self.Ath_table()
        self.az_table = self.Az_table()
        self.ge_table = self.Ge_table()

        self.ex_table_ref = None
        self.eth_table_ref = None
        self.ez_table_ref = None
        self.nuth_table_ref = None
        self.nuxz_table_ref = None
        self.rho_table_ref = None
        self.gzx_table_ref = None
        self.ax_table_ref = None
        self.ath_table_ref = None
        self.az_table_ref = None
        self.ge_table_ref = None

    @classmethod
    def add_card(cls, card: BDFCard, comment: str=''):
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
        return MATT3(mid, ex_table, eth_table, ez_table,
                     nuth_table, nuxz_table, rho_table, gzx_table,
                     ax_table, ath_table, az_table, ge_table, comment=comment)

    def Ex_table(self) -> int:
        return table_id(self.ex_table_ref, self.ex_table)

    def Eth_table(self) -> int:
        return table_id(self.eth_table_ref, self.eth_table)

    def Ez_table(self) -> int:
        return table_id(self.ez_table_ref, self.ez_table)

    def Nuth_table(self) -> int:
        return table_id(self.nuth_table_ref, self.nuth_table)

    def Nuxz_table(self) -> int:
        return table_id(self.nuxz_table_ref, self.nuxz_table)

    def Rho_table(self) -> int:
        return table_id(self.rho_table_ref, self.rho_table)

    def Gzx_table(self) -> int:
        return table_id(self.gzx_table_ref, self.gzx_table)

    def Ax_table(self) -> int:
        return table_id(self.ax_table_ref, self.ax_table)

    def Ath_table(self) -> int:
        return table_id(self.ath_table_ref, self.ath_table)

    def Az_table(self) -> int:
        return table_id(self.az_table_ref, self.az_table)

    def Ge_table(self) -> int:
        return table_id(self.ge_table_ref, self.ge_table)

    def raw_fields(self) -> list:
        list_fields = [
            'MATT3', self.Mid(), self.Ex_table(), self.Eth_table(), self.Ez_table(),
            self.Nuth_table(), self.Nuxz_table(), self.Rho_table(), None, None,
            self.Gzx_table(), self.Ax_table(), self.Ath_table(), self.Az_table(),
            None, self.Ge_table(),
        ]
        return list_fields

    def repr_fields(self):
        return self.raw_fields()

    def write_card(self, size: int=8, is_double: bool=False) -> str:
        card = self.repr_fields()
        if size == 8:
            return self.comment + print_card_8(card)
        return self.comment + print_card_16(card)


class MATT4(MaterialDependenceThermal):
    """
    Specifies temperature-dependent material properties on MAT2 entry
    fields via TABLEMi entries.

    +-------+-------+-------+-------+--------+-------+-------+---------+
    |   1   |   2   |   3   |   4   |   5    |   6   |   7   |    8    |
    +=======+=======+=======+=======+========+=======+=======+=========+
    | MATT4 |  MID  |  T(K) | T(CP) |        | T(H)  | T(mu) | T(HGEN) |
    +-------+-------+-------+-------+--------+-------+-------+---------+

    """
    type = 'MATT4'

    @classmethod
    def _init_from_empty(cls):
        mid = 1
        return MATT4(mid, k_table=None)

    def __init__(self, mid: int, k_table: int=0, cp_table: int=0, h_table=None,
                 mu_table: int=0, hgen_table: int=0, comment: str=''):
        MaterialDependenceThermal.__init__(self)
        if comment:
            self.comment = comment

        self.mid = mid
        self.k_table = k_table
        self.cp_table = cp_table
        self.h_table = h_table
        self.mu_table = mu_table
        self.hgen_table = hgen_table

        self.mid_ref = None
        self.k_table_ref = None
        self.cp_table_ref = None
        self.h_table_ref = None
        self.mu_table_ref = None
        self.hgen_table_ref = None

    @classmethod
    def add_card(cls, card: BDFCard, comment: str=''):
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
        return MATT4(mid, k_table, cp_table, h_table, mu_table,
                     hgen_table, comment=comment)

    @classmethod
    def add_op2_data(cls, data, comment=''):
        """
        Adds a MATT4 card from the OP2

        Parameters
        ----------
        data : list[varies]
            a list of fields defined in OP2 format
        comment : str; default=''
            a comment for the card

        """
        (mid, k_table, cp_table, null, h_table, mu_table, hgen_table) = data
        assert null == 0, data
        return MATT4(mid, k_table, cp_table, h_table, mu_table,
                     hgen_table, comment=comment)

    def cross_reference(self, model: BDF) -> None:
        """
        Cross links the card so referenced cards can be extracted directly

        Parameters
        ----------
        model : BDF()
            the BDF object

        """
        msg = ', which is required by MATT4 mid=%s' % self.mid
        self.mid_ref = model.Material(self.mid, msg=msg)

        _xref_table(self, model, 'k_table', msg=msg)
        _xref_table(self, model, 'cp_table', msg=msg)
        _xref_table(self, model, 'h_table', msg=msg)
        _xref_table(self, model, 'mu_table', msg=msg)
        _xref_table(self, model, 'hgen_table', msg=msg)

    def safe_cross_reference(self, model: BDF, xref_errors) -> None:
        """
        Cross links the card so referenced cards can be extracted directly

        Parameters
        ----------
        model : BDF()
            the BDF object

        """
        msg = ', which is required by MATT4 mid=%s' % self.mid
        self.mid_ref = model.safe_material(self.mid, self.mid, xref_errors, msg=msg)

        _safe_xref_table(self, model, 'k_table', xref_errors, msg=msg)
        _safe_xref_table(self, model, 'cp_table', xref_errors, msg=msg)
        _safe_xref_table(self, model, 'h_table', xref_errors, msg=msg)
        _safe_xref_table(self, model, 'mu_table', xref_errors, msg=msg)
        _safe_xref_table(self, model, 'hgen_table', xref_errors, msg=msg)

    def uncross_reference(self) -> None:
        """Removes cross-reference links"""
        self.mid = self.Mid()
        self.k_table = self.K_table()
        self.cp_table = self.Cp_table()
        self.h_table = self.H_table()
        self.mu_table = self.Mu_table()
        self.hgen_table = self.Hgen_table()

        self.mid_ref = None
        self.mid_ref = None
        self.k_table_ref = None
        self.cp_table_ref = None
        self.h_table_ref = None
        self.mu_table_ref = None
        self.hgen_table_ref = None

    def K_table(self) -> int:
        return self._get_table('k_table')

    def Cp_table(self) -> int:
        return self._get_table('cp_table')

    def H_table(self) -> int:
        return self._get_table('h_table')

    def Mu_table(self) -> int:
        return self._get_table('mu_table')

    def Hgen_table(self) -> int:
        return self._get_table('hgen_table')

    def raw_fields(self) -> list:
        list_fields = [
            'MATT4', self.Mid(), self.K_table(), self.Cp_table(),
            None,
            self.H_table(), self.Mu_table(), self.Hgen_table()
        ]
        return list_fields

    def repr_fields(self) -> list:
        return self.raw_fields()

    def write_card(self, size: int=8, is_double: bool=False) -> str:
        card = self.repr_fields()
        if size == 8:
            return self.comment + print_card_8(card)
        return self.comment + print_card_16(card)


class MATT5(MaterialDependenceThermal):
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
    type = 'MATT5'

    @classmethod
    def _init_from_empty(cls):
        mid = 1
        return MATT5(mid)

    def __init__(self, mid: int, kxx_table: int=0, kxy_table: int=0, kxz_table=None,
                 kyy_table: int=0, kyz_table: int=0, kzz_table: int=0,
                 cp_table: int=0, hgen_table: int=0, comment: str=''):
        MaterialDependenceThermal.__init__(self)
        if comment:
            self.comment = comment
        self.mid = mid
        self.kxx_table = kxx_table
        self.kxy_table = kxy_table
        self.kxz_table = kxz_table
        self.kyy_table = kyy_table
        self.kyz_table = kyz_table
        self.kzz_table = kzz_table
        self.cp_table = cp_table
        self.hgen_table = hgen_table

        self.mid_ref = None
        self.kxx_table_ref = None
        self.kxy_table_ref = None
        self.kxz_table_ref = None
        self.kyy_table_ref = None
        self.kyz_table_ref = None
        self.kzz_table_ref = None
        self.cp_table_ref = None
        self.hgen_table_ref = None

    @classmethod
    def add_card(cls, card: BDFCard, comment: str=''):
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
        return MATT5(mid, kxx_table, kxy_table, kxz_table, kyy_table,
                     kyz_table, kzz_table, cp_table, hgen_table,
                     comment=comment)

    @classmethod
    def add_op2_data(cls, data, comment: str=''):
        """
        Adds a MATT5 card from the OP2

        Parameters
        ----------
        data : list[varies]
            a list of fields defined in OP2 format
        comment : str; default=''
            a comment for the card

        """
        (mid, kxx_table, kxy_table, kxz_table, kyy_table, kyz_table, kzz_table,
         cp_table, null, hgen_table) = data

        if kxx_table == 0:
            kxx_table = None
        if kxy_table == 0:
            kxy_table = None
        if kxz_table == 0:
            kxz_table = None
        if kyy_table == 0:
            kyy_table = None
        if kyz_table == 0:
            kyz_table = None
        if kzz_table == 0:
            kzz_table = None
        if cp_table == 0:
            cp_table = None
        if hgen_table == 0:
            hgen_table = None
        assert null == 0, data
        return MATT5(mid, kxx_table, kxy_table, kxz_table, kyy_table,
                     kyz_table, kzz_table, cp_table, hgen_table,
                     comment=comment)

    def cross_reference(self, model: BDF) -> None:
        """
        Cross links the card so referenced cards can be extracted directly

        Parameters
        ----------
        model : BDF()
            the BDF object

        """
        msg = ', which is required by MATT5 mid=%s' % self.mid
        self.mid_ref = model.Material(self.mid, msg=msg)

        _xref_table(self, model, 'kxx_table', msg=msg)
        _xref_table(self, model, 'kxy_table', msg=msg)
        _xref_table(self, model, 'kxz_table', msg=msg)
        _xref_table(self, model, 'kyy_table', msg=msg)
        _xref_table(self, model, 'kyz_table', msg=msg)
        _xref_table(self, model, 'kzz_table', msg=msg)
        _xref_table(self, model, 'cp_table', msg=msg)
        _xref_table(self, model, 'hgen_table', msg=msg)

    def safe_cross_reference(self, model: BDF, xref_errors) -> None:
        """
        Cross links the card so referenced cards can be extracted directly

        Parameters
        ----------
        model : BDF()
            the BDF object

        """
        msg = ', which is required by MATT5 mid=%s' % self.mid
        self.mid_ref = model.safe_material(self.mid, self.mid, xref_errors, msg=msg)

        _safe_xref_table(self, model, 'kxx_table', xref_errors, msg=msg)
        _safe_xref_table(self, model, 'kxy_table', xref_errors, msg=msg)
        _safe_xref_table(self, model, 'kxz_table', xref_errors, msg=msg)
        _safe_xref_table(self, model, 'kyy_table', xref_errors, msg=msg)
        _safe_xref_table(self, model, 'kyz_table', xref_errors, msg=msg)
        _safe_xref_table(self, model, 'kzz_table', xref_errors, msg=msg)
        _safe_xref_table(self, model, 'cp_table', xref_errors, msg=msg)
        _safe_xref_table(self, model, 'hgen_table', xref_errors, msg=msg)

    def uncross_reference(self) -> None:
        """Removes cross-reference links"""
        self.mid = self.Mid()
        self.kxx_table = self.Kxx_table()
        self.kxy_table = self.Kxy_table()
        self.kxz_table = self.Kxz_table()
        self.kyy_table = self.Kyy_table()
        self.kyz_table = self.Kyz_table()
        self.kzz_table = self.Kzz_table()
        self.cp_table = self.Cp_table()
        self.hgen_table = self.Hgen_table()

        self.mid_ref = None
        self.kxx_table_ref = None
        self.kxy_table_ref = None
        self.kxz_table_ref = None
        self.kyy_table_ref = None
        self.kyz_table_ref = None
        self.kzz_table_ref = None
        self.cp_table_ref = None
        self.hgen_table_ref = None

    def Kxx_table(self) -> int:
        return self._get_table('kxx_table')

    def Kxy_table(self) -> int:
        return self._get_table('kxy_table')

    def Kxz_table(self) -> int:
        return self._get_table('kxz_table')

    def Kyy_table(self) -> int:
        return self._get_table('kyy_table')

    def Kyz_table(self) -> int:
        return self._get_table('kyz_table')

    def Kzz_table(self) -> int:
        return self._get_table('kzz_table')

    def Cp_table(self) -> int:
        return self._get_table('cp_table')

    def Hgen_table(self) -> int:
        return self._get_table('hgen_table')

    def raw_fields(self):
        list_fields = ['MATT5', self.Mid(),
                       self.Kxx_table(), self.Kxy_table(), self.Kxz_table(),
                       self.Kyy_table(), self.Kyz_table(), self.Kzz_table(),
                       self.Cp_table(), None, self.Hgen_table()]
        return list_fields

    def repr_fields(self):
        return self.raw_fields()

    def write_card(self, size: int=8, is_double: bool=False) -> str:
        card = self.repr_fields()
        if size == 8:
            return self.comment + print_card_8(card)
        return self.comment + print_card_16(card)


class MATT8(MaterialDependenceThermal):
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
    type = 'MATT8'

    @classmethod
    def _init_from_empty(cls):
        mid = 1
        return MATT8(mid, e1_table=None)

    def __init__(self, mid: int, e1_table: int=0, e2_table: int=0, nu12_table: int=0,
                 g12_table: int=0, g1z_table: int=0, g2z_table: int=0,
                 rho_table: int=0,
                 a1_table: int=0, a2_table: int=0,
                 xt_table: int=0, xc_table: int=0,
                 yt_table: int=0, yc_table: int=0,
                 s_table: int=0, ge_table: int=0,
                 f12_table: int=0, comment: str=''):
        MaterialDependenceThermal.__init__(self)
        if comment:
            self.comment = comment

        self.mid = mid
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

        self.mid_ref = None

        self.e1_table_ref = None
        self.e2_table_ref = None
        self.nu12_table_ref = None
        self.g12_table_ref = None
        self.g1z_table_ref = None
        self.g2z_table_ref = None
        self.rho_table_ref = None

        self.a1_table_ref = None
        self.a2_table_ref = None
        self.xt_table_ref = None
        self.xc_table_ref = None
        self.yt_table_ref = None
        self.yc_table_ref = None
        self.s_table_ref = None
        self.ge_table_ref = None
        self.f12_table_ref = None

    @classmethod
    def add_card(cls, card: BDFCard, comment: str=''):
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
        e1_table = integer_or_blank(card, 2, 'T(E1)', default=0)
        e2_table = integer_or_blank(card, 3, 'T(E2)', default=0)
        nu12_table = integer_or_blank(card, 4, 'T(Nu12)', default=0)
        g12_table = integer_or_blank(card, 5, 'T(G12)', default=0)
        g1z_table = integer_or_blank(card, 6, 'T(G1z)', default=0)
        g2z_table = integer_or_blank(card, 7, 'T(G2z)', default=0)
        rho_table = integer_or_blank(card, 8, 'T(Rho)', default=0)
        a1_table = integer_or_blank(card, 9, 'T(A1)', default=0)
        a2_table = integer_or_blank(card, 10, 'T(A2)', default=0)

        xt_table = integer_or_blank(card, 12, 'T(Xt)', default=0)
        xc_table = integer_or_blank(card, 13, 'T(Xc)', default=0)
        yt_table = integer_or_blank(card, 14, 'T(Yt)', default=0)
        yc_table = integer_or_blank(card, 15, 'T(Yc)', default=0)
        s_table = integer_or_blank(card, 16, 'T(S)', default=0)
        ge_table = integer_or_blank(card, 17, 'T(GE)', default=0)
        f12_table = integer_or_blank(card, 18, 'T(F12)', default=0)

        assert len(card) <= 19, f'len(MATT8 card) = {len(card):d}\ncard={card}'
        return MATT8(mid, e1_table, e2_table, nu12_table, g12_table,
                     g1z_table, g2z_table, rho_table,
                     a1_table, a2_table, xt_table,
                     xc_table, yt_table, yc_table,
                     s_table, ge_table, f12_table,
                     comment=comment)

    def cross_reference(self, model: BDF) -> None:
        """
        Cross links the card so referenced cards can be extracted directly

        Parameters
        ----------
        model : BDF()
            the BDF object

        """
        msg = ', which is required by MATT8 mid=%s' % self.mid
        self.mid_ref = model.Material(self.mid, msg=msg)

        _xref_table(self, model, 'e1_table', msg=msg)
        _xref_table(self, model, 'e2_table', msg=msg)
        _xref_table(self, model, 'nu12_table', msg=msg)
        _xref_table(self, model, 'g12_table', msg=msg)
        _xref_table(self, model, 'g1z_table', msg=msg)
        _xref_table(self, model, 'g2z_table', msg=msg)
        _xref_table(self, model, 'rho_table', msg=msg)
        _xref_table(self, model, 'a1_table', msg=msg)
        _xref_table(self, model, 'a2_table', msg=msg)
        _xref_table(self, model, 'xt_table', msg=msg)
        _xref_table(self, model, 'xc_table', msg=msg)
        _xref_table(self, model, 'yt_table', msg=msg)
        _xref_table(self, model, 'yc_table', msg=msg)
        _xref_table(self, model, 's_table', msg=msg)
        _xref_table(self, model, 'ge_table', msg=msg)
        _xref_table(self, model, 'f12_table', msg=msg)

    def safe_cross_reference(self, model: BDF, xref_errors) -> None:
        """
        Cross links the card so referenced cards can be extracted directly

        Parameters
        ----------
        model : BDF()
            the BDF object

        """
        msg = ', which is required by MATT8 mid=%s' % self.mid
        self.mid_ref = model.safe_material(self.mid, self.mid, xref_errors, msg=msg)

        _safe_xref_table(self, model, 'e1_table', xref_errors, msg=msg)
        _safe_xref_table(self, model, 'e2_table', xref_errors, msg=msg)
        _safe_xref_table(self, model, 'nu12_table', xref_errors, msg=msg)
        _safe_xref_table(self, model, 'g12_table', xref_errors, msg=msg)
        _safe_xref_table(self, model, 'g1z_table', xref_errors, msg=msg)
        _safe_xref_table(self, model, 'g2z_table', xref_errors, msg=msg)
        _safe_xref_table(self, model, 'rho_table', xref_errors, msg=msg)
        _safe_xref_table(self, model, 'a1_table', xref_errors, msg=msg)
        _safe_xref_table(self, model, 'a2_table', xref_errors, msg=msg)
        _safe_xref_table(self, model, 'xt_table', xref_errors, msg=msg)
        _safe_xref_table(self, model, 'xc_table', xref_errors, msg=msg)
        _safe_xref_table(self, model, 'yt_table', xref_errors, msg=msg)
        _safe_xref_table(self, model, 'yc_table', xref_errors, msg=msg)
        _safe_xref_table(self, model, 's_table', xref_errors, msg=msg)
        _safe_xref_table(self, model, 'ge_table', xref_errors, msg=msg)
        _safe_xref_table(self, model, 'f12_table', xref_errors, msg=msg)

    def uncross_reference(self) -> None:
        """Removes cross-reference links"""
        self.e1_table = self.E1_table()
        self.e2_table = self.E2_table()
        self.nu12_table = self.Nu12_table()
        self.g12_table = self.G12_table()
        self.g1z_table = self.G1z_table()
        self.g2z_table = self.G2z_table()
        self.rho_table = self.Rho_table()
        self.a1_table = self.A1_table()
        self.a2_table = self.A2_table()

        self.xt_table = self.Xt_table()
        self.xc_table = self.Xc_table()
        self.yt_table = self.Yt_table()
        self.yc_table = self.Yc_table()
        self.s_table = self.S_table()
        self.ge_table = self.Ge_table()
        self.f12_table = self.F12_table()

        self.e1_table_ref = None
        self.e2_table_ref = None
        self.nu12_table_ref = None
        self.g12_table_ref = None
        self.g1z_table_ref = None
        self.g2z_table_ref = None
        self.rho_table_ref = None

        self.a1_table_ref = None
        self.a2_table_ref = None
        self.xt_table_ref = None
        self.xc_table_ref = None
        self.yt_table_ref = None
        self.yc_table_ref = None
        self.s_table_ref = None
        self.ge_table_ref = None
        self.f12_table_ref = None

    def E1_table(self) -> int:
        return table_id(self.e1_table_ref, self.e1_table)

    def E2_table(self) -> int:
        return table_id(self.e2_table_ref, self.e2_table)

    def Nu12_table(self) -> int:
        return table_id(self.nu12_table_ref, self.nu12_table)

    def G12_table(self) -> int:
        return table_id(self.g12_table_ref, self.g12_table)

    def G1z_table(self) -> int:
        return table_id(self.g1z_table_ref, self.g1z_table)

    def G2z_table(self) -> int:
        return table_id(self.g2z_table_ref, self.g2z_table)

    def Rho_table(self) -> int:
        return table_id(self.rho_table_ref, self.rho_table)

    def A1_table(self) -> int:
        return table_id(self.a1_table_ref, self.a1_table)

    def A2_table(self) -> int:
        return table_id(self.a2_table_ref, self.a2_table)

    def S_table(self) -> int:
        return table_id(self.s_table_ref, self.s_table)

    def Ge_table(self) -> int:
        return table_id(self.ge_table_ref, self.ge_table)

    def F12_table(self) -> int:
        return table_id(self.f12_table_ref, self.f12_table)

    def Xt_table(self) -> int:
        return table_id(self.xt_table_ref, self.xt_table)

    def Xc_table(self) -> int:
        return table_id(self.xc_table_ref, self.xc_table)

    def Yt_table(self) -> int:
        return table_id(self.yt_table_ref, self.yt_table)

    def Yc_table(self) -> int:
        return table_id(self.yc_table_ref, self.yc_table)

    def raw_fields(self) -> list:
        list_fields = ['MATT8', self.mid, self.E1_table(), self.E2_table(), self.Nu12_table(),
                       self.G12_table(), self.G1z_table(), self.G2z_table(), self.Rho_table(),
                       self.A1_table(), self.A2_table(), None,
                       self.Xt_table(), self.Xc_table(), self.Yt_table(), self.Yc_table(),
                       self.S_table(), self.Ge_table(), self.F12_table()]
        return list_fields

    def write_card(self, size: int=8, is_double: bool=False) -> str:
        """
        +--------+--------+--------+--------+--------+--------+--------+--------+--------+
        |    1   |   2    |    3   |    4   |    5   |    6   |   7    |    8   |   9    |
        +========+========+========+========+========+========+========+========+========+
        | MATT9  |  MID   | T(G11) | T(G12) | T(G13) | T(G14) | T(G15) | T(G16) | T(G22) |
        +--------+--------+--------+--------+--------+--------+--------+--------+--------+
        |        | T(G23) | T(G24) | T(G25) | T(G26) | T(G33) | T(G34) | T(G35) | T(G36) |
        +--------+--------+--------+--------+--------+--------+--------+--------+--------+
        |        | T(G44) | T(G45) | T(G46) | T(G55) | T(G56) | T(G66) | T(RHO) | T(A1)  |
        +--------+--------+--------+--------+--------+--------+--------+--------+--------+
        |        | T(A2)  | T(A3)  | T(A4)  | T(A5)  | T(A6)  |        |  T(GE) |        |
        +--------+--------+--------+--------+--------+--------+--------+--------+--------+

        """
        list_fields = self.raw_fields()
        return self.comment + print_card_8(list_fields)


class MATT9(MaterialDependenceThermal):
    type = 'MATT9'

    @classmethod
    def _init_from_empty(cls):
        mid = 1
        return MATT9(mid, g11_table=None)

    def __init__(self, mid,
                 g11_table: int=0, g12_table: int=0, g13_table: int=0, g14_table: int=0,
                 g15_table: int=0, g16_table: int=0,
                 g22_table: int=0, g23_table: int=0, g24_table: int=0,
                 g25_table: int=0, g26_table: int=0,
                 g33_table: int=0, g34_table: int=0, g35_table: int=0, g36_table: int=0,
                 g44_table: int=0, g45_table: int=0, g46_table: int=0,
                 g55_table: int=0, g56_table: int=0,
                 g66_table: int=0,
                 rho_table: int=0,
                 a1_table: int=0, a2_table: int=0, a3_table: int=0,
                 a4_table: int=0, a5_table: int=0, a6_table: int=0,
                 ge_table: int=0,
                 comment: str=''):
        assert isinstance(ge_table, int), ge_table
        MaterialDependenceThermal.__init__(self)
        if comment:
            self.comment = comment

        self.mid = mid
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
        self.mid_ref = None

    @classmethod
    def add_card(cls, card: BDFCard, comment: str=''):
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
        return MATT9(mid, g11_table, g12_table, g13_table, g14_table, g15_table, g16_table,
                     g22_table, g23_table, g24_table, g25_table, g26_table,
                     g33_table, g34_table, g35_table, g36_table,
                     g44_table, g45_table, g46_table,
                     g55_table, g56_table, g66_table,
                     rho_table,
                     a1_table, a2_table, a3_table, a4_table, a5_table, a6_table,
                     ge_table, comment=comment)

    def cross_reference(self, model: BDF) -> None:
        """
        Cross links the card so referenced cards can be extracted directly

        Parameters
        ----------
        model : BDF()
            the BDF object

        """
        msg = ', which is required by MATT9 mid=%s' % self.mid
        self.mid_ref = model.Material(self.mid, msg=msg)

        #if self.e1_table is not None:
            #self.e1_table_ref = model.TableM(self.e1_table)
        #if self.e2_table is not None:
            #self.e2_table_ref = model.TableM(self.e2_table)

    def safe_cross_reference(self, model: BDF, xref_errors) -> None:
        self.cross_reference(model)

    def uncross_reference(self) -> None:
        """Removes cross-reference links"""
        pass
        #self.e1_table = self.E1_table()
        #self.e2_table = self.E2_table()
        #self.e1_table_ref = None
        #self.e2_table_ref = None

    #def E1_table(self):
        #if self.e1_table_ref is not None:
            #return self.e1_table_ref.tid
        #return self.e1_table

    def raw_fields(self):
        list_fields = [
            'MATT9', self.mid,
            self.g11_table, self.g12_table, self.g13_table, self.g14_table, self.g15_table, self.g16_table,
            self.g22_table, self.g23_table, self.g24_table, self.g25_table, self.g26_table,
            self.g33_table, self.g34_table, self.g35_table, self.g36_table,
            self.g44_table, self.g45_table, self.g46_table,
            self.g55_table, self.g56_table, self.g66_table,
            self.rho_table,
            self.a1_table, self.a2_table, self.a3_table,
            self.a4_table, self.a5_table, self.a6_table,
            self.ge_table,
        ]
        return list_fields

    def write_card(self, size: int=8, is_double: bool=False) -> str:
        list_fields = self.raw_fields()
        return self.comment + print_card_8(list_fields)


class MATT11(MaterialDependenceThermal):
    """
    +--------+-----+-----+-----+----+------+------+------+-----+
    |    1   |  2  |  3  |  4  |  5 |   6  |   7  |   8  |  9  |
    +========+=====+=====+=====+====+======+======+======+=====+
    | MATT11 | MID | E1  | E2  | E3 | NU12 | NU13 | NU23 | G12 |
    +--------+-----+-----+-----+----+------+------+------+-----+
    |        | G13 | G23 | RHO | A1 |  A2  |  A3  |      |  GE |
    +--------+-----+-----+-----+----+------+------+------+-----+
    """
    type = 'MATT11'

    @classmethod
    def _init_from_empty(cls):
        mid = 1
        return MATT11(mid, e1_table=None, e2_table=None, e3_table=None,
                      nu12_table=None, nu13_table=None, nu23_table=None,
                      g12_table=None, g13_table=None, g23_table=None,
                      rho_table=None,
                      a1_table=None, a2_table=None, a3_table=None,
                      ge_table=None, comment='')

    def __init__(self, mid: int,
                 e1_table: int=0, e2_table: int=0, e3_table: int=0,
                 nu12_table: int=0, nu13_table: int=0, nu23_table: int=0,
                 g12_table: int=0, g13_table: int=0, g23_table: int=0,
                 rho_table: int=0,
                 a1_table: int=0, a2_table: int=0, a3_table: int=0,
                 ge_table: int=0,
                 comment: str=''):
        MaterialDependenceThermal.__init__(self)
        if comment:
            self.comment = comment

        self.mid = mid
        self.e1_table = e1_table
        self.e2_table = e2_table
        self.e3_table = e3_table
        self.nu12_table = nu12_table
        self.nu13_table = nu13_table
        self.nu23_table = nu23_table
        self.g12_table = g12_table
        self.g13_table = g13_table
        self.g23_table = g23_table
        self.rho_table = rho_table
        self.a1_table = a1_table
        self.a2_table = a2_table
        self.a3_table = a3_table
        self.ge_table = ge_table
        self.mid_ref = None

        self.e1_table_ref = None
        self.e2_table_ref = None
        self.e3_table_ref = None
        self.nu12_table_ref = None
        self.nu13_table_ref = None
        self.nu23_table_ref = None
        self.g12_table_ref = None
        self.g13_table_ref = None
        self.g23_table_ref = None
        self.a1_table_ref = None
        self.a2_table_ref = None
        self.a3_table_ref = None

    @classmethod
    def add_card(cls, card: BDFCard, comment: str=''):
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
        e1_table = integer_or_blank(card, 2, 'T(E1)')
        e2_table = integer_or_blank(card, 3, 'T(E2)')
        e3_table = integer_or_blank(card, 4, 'T(E3)')
        nu12_table = integer_or_blank(card, 5, 'T(NU12)')
        nu13_table = integer_or_blank(card, 6, 'T(NU13)')
        nu23_table = integer_or_blank(card, 7, 'T(NU23)')
        g12_table = integer_or_blank(card, 8, 'T(G12)')

        g13_table = integer_or_blank(card, 9, 'T(G13)')
        g23_table = integer_or_blank(card, 10, 'T(G23)')
        rho_table = integer_or_blank(card, 11, 'T(RHO)')
        a1_table = integer_or_blank(card, 12, 'T(A1)')
        a2_table = integer_or_blank(card, 13, 'T(A3)')
        a3_table = integer_or_blank(card, 14, 'T(A3)')
        #15
        ge_table = integer_or_blank(card, 16, 'T(GE)')

        assert len(card) <= 17, f'len(MATT11 card) = {len(card):d}\ncard={card}'
        return MATT11(mid, e1_table=e1_table, e2_table=e2_table, e3_table=e3_table,
                      nu12_table=nu12_table, nu13_table=nu13_table, nu23_table=nu23_table,
                      g12_table=g12_table, g13_table=g13_table, g23_table=g23_table,
                      a1_table=a1_table, a2_table=a2_table, a3_table=a3_table,
                      rho_table=rho_table, ge_table=ge_table, comment=comment)

    def cross_reference(self, model: BDF) -> None:
        """
        Cross links the card so referenced cards can be extracted directly

        Parameters
        ----------
        model : BDF()
            the BDF object

        """
        msg = ', which is required by MATT9 mid=%s' % self.mid
        self.mid_ref = model.Material(self.mid, msg=msg)

        #if self.e1_table is not None:
            #self.e1_table_ref = model.TableM(self.e1_table)
        #if self.e2_table is not None:
            #self.e2_table_ref = model.TableM(self.e2_table)

    def safe_cross_reference(self, model: BDF, xref_errors) -> None:
        self.cross_reference(model)

    def uncross_reference(self) -> None:
        """Removes cross-reference links"""
        pass
        #self.e1_table = self.E1_table()
        #self.e2_table = self.E2_table()
        #self.e1_table_ref = None
        #self.e2_table_ref = None

    #def E1_table(self):
        #if self.e1_table_ref is not None:
            #return self.e1_table_ref.tid
        #return self.e1_table

    def raw_fields(self):
        list_fields = [
            'MATT11', self.mid,
            self.e1_table, self.e2_table, self.e3_table,
            self.nu12_table, self.nu13_table, self.nu23_table,
            self.g12_table, self.g13_table, self.g23_table,
            self.rho_table,
            self.a1_table, self.a2_table, self.a3_table,
            None, self.ge_table,
        ]
        return list_fields

    def write_card(self, size: int=8, is_double: bool=False) -> str:
        list_fields = self.raw_fields()
        return self.comment + print_card_8(list_fields)
