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

All cards are Material objects.

"""
#pylint: disable=E1103,C0103,C0111
from __future__ import annotations
from typing import TYPE_CHECKING

from pyNastran.bdf.cards.base_card import BaseCard
from pyNastran.bdf.bdf_interface.assign_type import (
    integer, integer_or_blank, double, double_or_blank, string)
from pyNastran.bdf.field_writer_8 import print_card_8
from pyNastran.bdf.field_writer_16 import print_card_16
if TYPE_CHECKING:  # pragma: no cover
    from pyNastran.bdf.bdf import BDF


class MaterialDependence(BaseCard):
    def __init__(self):
        self.mid = None

    def Mid(self):
        if self.mid_ref is None:
            return self.mid
        return self.mid_ref.mid  # TODO: is this something that should be supported?

    def _get_table(self, key):
        """internal method for accessing tables"""
        table = getattr(self, key)
        table_ref = getattr(self, key + '_ref')
        if table_ref is not None:
            return table_ref.tid
        return table

class MaterialDependenceThermal(MaterialDependence):
    def __init__(self):
        MaterialDependence.__init__(self)

    def _xref_table(self, model, key, msg):
        slot = getattr(self, key)
        if slot is not None:
            setattr(self, key + '_ref', model.TableM(slot, msg + f' for {key}'))

class MATS1(MaterialDependence):
    """
    Specifies stress-dependent material properties for use in applications
    involving nonlinear materials. This entry is used if a MAT1, MAT2 or MAT9
    entry is specified with the same MID in a nonlinear solution sequence
    (SOLs 106 and 129).

    """
    type = 'MATS1'

    def __init__(self, mid, tid, Type, h, hr, yf, limit1, limit2, comment=''):
        MaterialDependence.__init__(self)
        if comment:
            self.comment = comment
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
        self.tid_ref = None
        self.mid_ref = None

    @classmethod
    def _init_from_empty(cls):
        mid = 1
        tid = 1
        Type = None
        h = None
        hr = None
        yf = None
        limit1 = None
        limit2 = None
        return MATS1(mid, tid, Type, h, hr, yf, limit1, limit2, comment='')

    def validate(self):
        if self.Type not in ['NLELAST', 'PLASTIC', 'PLSTRN']:
            raise ValueError('MATS1 Type must be [NLELAST, PLASTIC, PLSTRN]; Type=%r' % self.Type)

    @classmethod
    def add_card(cls, card, comment=''):
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
        tid = integer_or_blank(card, 2, 'tid')
        Type = string(card, 3, 'Type')

        if Type not in ['NLELAST', 'PLASTIC', 'PLSTRN']:
            raise ValueError('MATS1 Type must be [NLELAST, PLASTIC, PLSTRN]; Type=%r' % Type)
        if Type == 'NLELAST':
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
            yf = integer_or_blank(card, 5, 'yf', 1)
            hr = integer_or_blank(card, 6, 'hr', 1)
            limit1 = double(card, 7, 'limit1')

            if yf in [3, 4]:
                limit2 = double(card, 8, 'limit2')
            else:
                #limit2 = blank(card, 8, 'limit2')
                limit2 = None
        assert len(card) <= 9, 'len(MATS1 card) = %i\ncard=%s' % (len(card), card)
        return MATS1(mid, tid, Type, h, hr, yf, limit1, limit2, comment=comment)

    @classmethod
    def add_op2_data(cls, data, comment=''):
        """
        Adds a MATS1 card from the OP2

        Parameters
        ----------
        data : List[varies]
            a list of fields defined in OP2 format
        comment : str; default=''
            a comment for the card

        """
        (mid, tid, Type, h, yf, hr, limit1, limit2) = data
        if Type == 1:
            Type = 'NLELAST'
        elif Type == 2:
            Type = 'PLASTIC'
        elif Type == 3:
            Type = 'PLSTRN'
        else:  # pragma: no cover
            raise RuntimeError(f'Invalid Type:  mid={mid}; Type={Type}; must be 1=NLELAST, '
                               '2=PLASTIC, or 3=PLSTRN')
        return MATS1(mid, tid, Type, h, hr, yf, limit1, limit2, comment=comment)

    def Yf(self):
        d = {1: 'VonMises', 2: 'Tresca', 3: 'MohrCoulomb', 4: 'Drucker-Prager'}
        return d[self.yf]

    def Hf(self):
        d = {1: 'Isotropic', 2: 'Kinematic', 3: 'Combined'}
        return d[self.hr]

    def E(self, strain):
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
        if self.tid:  # then self.h is used
            self.tid_ref = model.Table(self.tid, msg=msg) # TABLES1 or TABLEST

    def uncross_reference(self) -> None:
        """Removes cross-reference links"""
        self.mid = self.Mid()
        if self.tid:
            self.tid = self.Tid()
        self.tid_ref = None
        self.mid_ref = None

    def Tid(self):
        if self.tid_ref is None:
            return self.tid
        return self.tid_ref.tid

    def raw_fields(self):
        list_fields = ['MATS1', self.Mid(), self.Tid(), self.Type,
                       self.h, self.yf, self.hr, self.limit1, self.limit2]
        return list_fields

    def repr_fields(self):
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
        return MATT1(mid, e_table=None, g_table=None, nu_table=None, rho_table=None,
                     a_table=None, ge_table=None, st_table=None,
                     sc_table=None, ss_table=None, comment='')

    def __init__(self, mid, e_table=None, g_table=None, nu_table=None,
                 rho_table=None, a_table=None, ge_table=None, st_table=None,
                 sc_table=None, ss_table=None, comment=''):
        MaterialDependenceThermal.__init__(self)
        if comment:
            self.comment = comment
        self.mid = mid
        if e_table == 0:
            e_table = None
        if g_table == 0:
            g_table = None
        if nu_table == 0:
            nu_table = None
        if rho_table == 0:
            rho_table = None
        if a_table == 0:
            a_table = None
        if ge_table == 0:
            ge_table = None
        if st_table == 0:
            st_table = None
        if sc_table == 0:
            sc_table = None
        if ss_table == 0:
            ss_table = None

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
    def add_card(cls, card, comment=''):
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

        assert len(card) <= 12, 'len(MATT1 card) = %i\ncard=%s' % (len(card), card)
        return MATT1(mid, e_table, g_table, nu_table, rho_table, a_table,
                     ge_table, st_table, sc_table, ss_table, comment=comment)

    @classmethod
    def add_op2_data(cls, data, comment=''):
        """
        Adds a MATT1 card from the OP2

        Parameters
        ----------
        data : List[varies]
            a list of fields defined in OP2 format
        comment : str; default=''
            a comment for the card

        """
        (mid, E_table, G_table, nu_table, rho_table, A_table, dunno_a, ge_table,
         st_table, sc_table, ss_table, dunno_b) = data
        if E_table == 0:
            E_table = None
        elif E_table > 100000000:
            E_table = -(E_table - 100000000)

        if G_table == 0:
            G_table = None
        elif G_table > 100000000:
            G_table = -(G_table - 100000000)

        if nu_table == 0:
            nu_table = None
        elif nu_table > 100000000:
            nu_table = -(nu_table - 100000000)

        if rho_table == 0:
            rho_table = None
        elif rho_table > 100000000:
            rho_table = -(rho_table - 100000000)

        if E_table == 0:
            E_table = None
        if A_table > 100000000:
            A_table = -(A_table - 100000000)

        mat = MATT1(mid, E_table, G_table, nu_table, rho_table, A_table,
                    ge_table, st_table, sc_table, ss_table, comment=comment)
        assert dunno_a == 0, '%s; dunno_a=%s\n%s' % (data, dunno_a, str(mat))
        assert dunno_b == 0, '%s; dunno_b=%s\n%s' % (data, dunno_b, str(mat))
        return mat

    def E(self, temperature):
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

        self._xref_table(model, 'e_table', msg=msg)
        self._xref_table(model, 'g_table', msg=msg)
        self._xref_table(model, 'nu_table', msg=msg)
        self._xref_table(model, 'rho_table', msg=msg)
        self._xref_table(model, 'a_table', msg=msg)
        self._xref_table(model, 'ge_table', msg=msg)
        self._xref_table(model, 'st_table', msg=msg)
        self._xref_table(model, 'sc_table', msg=msg)
        self._xref_table(model, 'ss_table', msg=msg)

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

    def E_table(self):
        return self._get_table('e_table')

    def G_table(self):
        return self._get_table('g_table')

    def Nu_table(self):
        return self._get_table('nu_table')

    def Rho_table(self):
        return self._get_table('rho_table')

    def A_table(self):
        return self._get_table('a_table')

    def Ge_table(self):
        return self._get_table('ge_table')

    def St_table(self):
        return self._get_table('st_table')

    def Sc_table(self):
        return self._get_table('sc_table')

    def Ss_table(self):
        return self._get_table('ss_table')

    def raw_fields(self):
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
        return MATT2(mid, g11_table=None, g12_table=None, g13_table=None, g22_table=None,
                     g23_table=None, g33_table=None, rho_table=None,
                     a1_table=None, a2_table=None, a3_table=None, ge_table=None,
                     st_table=None, sc_table=None, ss_table=None, comment='')

    def __init__(self, mid, g11_table=None, g12_table=None, g13_table=None,
                 g22_table=None, g23_table=None, g33_table=None, rho_table=None,
                 a1_table=None, a2_table=None, a3_table=None,
                 ge_table=None, st_table=None, sc_table=None, ss_table=None, comment=''):
        MaterialDependenceThermal.__init__(self)
        if comment:
            self.comment = comment

        if g12_table == 0:
            g12_table = None
        if g13_table == 0:
            g13_table = None
        if g22_table == 0:
            g22_table = None
        if g23_table == 0:
            g23_table = None
        if g33_table == 0:
            g33_table = None
        if rho_table == 0:
            rho_table = None

        if a1_table == 0:
            a1_table = None
        if a2_table == 0:
            a2_table = None
        if a3_table == 0:
            a3_table = None
        if ge_table == 0:
            ge_table = None

        if st_table == 0:
            st_table = None
        if sc_table == 0:
            sc_table = None
        if ss_table == 0:
            ss_table = None
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

    @classmethod
    def add_card(cls, card, comment=''):
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
        g11_table = integer_or_blank(card, 2, 'T(G11)')
        g12_table = integer_or_blank(card, 3, 'T(G12)')
        g13_table = integer_or_blank(card, 4, 'T(G13)')
        g22_table = integer_or_blank(card, 5, 'T(G22)')
        g23_table = integer_or_blank(card, 6, 'T(G23)')
        g33_table = integer_or_blank(card, 7, 'T(G33)')
        rho_table = integer_or_blank(card, 8, 'T(rho)')
        a1_table = integer_or_blank(card, 9, 'T(A1)')
        a2_table = integer_or_blank(card, 10, 'T(A2)')
        a3_table = integer_or_blank(card, 11, 'T(A3)')
        ge_table = integer_or_blank(card, 13, 'T(ge)')
        st_table = integer_or_blank(card, 14, 'T(st)')
        sc_table = integer_or_blank(card, 15, 'T(sc)')
        ss_table = integer_or_blank(card, 16, 'T(ss)')

        assert len(card) <= 17, 'len(MATT2 card) = %i\ncard=%s' % (len(card), card)
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

        self._xref_table(model, 'g11_table', msg=msg)
        self._xref_table(model, 'g12_table', msg=msg)
        self._xref_table(model, 'g13_table', msg=msg)
        self._xref_table(model, 'g22_table', msg=msg)
        self._xref_table(model, 'g23_table', msg=msg)
        self._xref_table(model, 'g33_table', msg=msg)
        self._xref_table(model, 'rho_table', msg=msg)
        self._xref_table(model, 'a1_table', msg=msg)
        self._xref_table(model, 'a2_table', msg=msg)
        self._xref_table(model, 'a3_table', msg=msg)
        self._xref_table(model, 'ge_table', msg=msg)
        self._xref_table(model, 'st_table', msg=msg)
        self._xref_table(model, 'sc_table', msg=msg)
        self._xref_table(model, 'ss_table', msg=msg)

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

    def G11_table(self):
        return self._get_table('g11_table')

    def G12_table(self):
        return self._get_table('g12_table')

    def G13_table(self):
        return self._get_table('g13_table')

    def G22_table(self):
        return self._get_table('g22_table')

    def G23_table(self):
        return self._get_table('g23_table')

    def G33_table(self):
        return self._get_table('g33_table')

    def Rho_table(self):
        return self._get_table('rho_table')

    def A1_table(self):
        return self._get_table('a1_table')

    def A2_table(self):
        return self._get_table('a2_table')

    def A3_table(self):
        return self._get_table('a3_table')

    def Ge_table(self):
        return self._get_table('ge_table')

    def St_table(self):
        return self._get_table('st_table')

    def Sc_table(self):
        return self._get_table('sc_table')

    def Ss_table(self):
        return self._get_table('ss_table')

    def raw_fields(self):
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
        return MATT3(mid, ex_table=None, eth_table=None, ez_table=None, nuth_table=None,
                     nuxz_table=None, rho_table=None, gzx_table=None,
                     ax_table=None, ath_table=None, az_table=None, ge_table=None, comment='')

    def __init__(self, mid, ex_table=None, eth_table=None, ez_table=None,
                 nuth_table=None, nuxz_table=None, rho_table=None,
                 gzx_table=None, ax_table=None, ath_table=None, az_table=None,
                 ge_table=None, comment=''):
        """
        Creates a MATT3 card
        """
        MaterialDependenceThermal.__init__(self)
        if comment:
            self.comment = comment

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

    def cross_reference(self, model: BDF) -> None:
        msg = ', which is required by MATT3 mid=%s' % self.mid
        self.mid_ref = model.Material(self.mid, msg=msg)

        #self._get_table('ex_table')
        if self.ex_table is not None:
            self.ex_table_ref = model.TableM(self.ex_table)
        if self.eth_table is not None:
            self.eth_table_ref = model.TableM(self.eth_table)
        if self.ez_table is not None:
            self.ez_table_ref = model.TableM(self.ez_table)
        if self.nuth_table is not None:
            self.nuth_table_ref = model.TableM(self.nuth_table)
        if self.nuxz_table is not None:
            self.nuxz_table_ref = model.TableM(self.nuxz_table)
        if self.rho_table is not None:
            self.rho_table_ref = model.TableM(self.rho_table)

        if self.gzx_table is not None:
            self.gzx_table_ref = model.TableM(self.gzx_table)
        if self.ax_table is not None:
            self.ax_table_ref = model.TableM(self.ax_table)
        if self.ath_table is not None:
            self.ath_table_ref = model.TableM(self.ath_table)
        if self.az_table is not None:
            self.az_table_ref = model.TableM(self.az_table)
        if self.ge_table is not None:
            self.ge_table_ref = model.TableM(self.ge_table)

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
    def add_card(cls, card, comment=''):
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
        ex_table = integer_or_blank(card, 2, 'T(EX)')
        eth_table = integer_or_blank(card, 3, 'T(ETH)')
        ez_table = integer_or_blank(card, 5, 'T(EZ)')
        nuth_table = integer_or_blank(card, 6, 'T(NUTH)')
        nuxz_table = integer_or_blank(card, 7, 'T(NUXZ)')
        rho_table = integer_or_blank(card, 8, 'T(RHO)')

        gzx_table = integer_or_blank(card, 11, 'T(GZX)')
        ax_table = integer_or_blank(card, 12, 'T(AX)')
        ath_table = integer_or_blank(card, 13, 'T(ATH)')
        az_table = integer_or_blank(card, 14, 'T(AZ)')
        ge_table = integer_or_blank(card, 16, 'T(GE)')

        assert len(card) <= 16, 'len(MATT3 card) = %i\ncard=%s' % (len(card), card)
        return MATT3(mid, ex_table, eth_table, ez_table,
                     nuth_table, nuxz_table, rho_table, gzx_table,
                     ax_table, ath_table, az_table, ge_table, comment=comment)

    def Ex_table(self):
        if self.ex_table_ref is not None:
            return self.ex_table_ref.tid
        return self.ex_table

    def Eth_table(self):
        if self.eth_table_ref is not None:
            return self.eth_table_ref.tid
        return self.eth_table

    def Ez_table(self):
        if self.ez_table_ref is not None:
            return self.ez_table_ref.tid
        return self.eth_table

    def Nuth_table(self):
        if self.nuth_table_ref is not None:
            return self.nuth_table_ref.tid
        return self.nuth_table

    def Nuxz_table(self):
        if self.nuxz_table_ref is not None:
            return self.nuxz_table_ref.tid
        return self.nuxz_table

    def Rho_table(self):
        if self.rho_table_ref is not None:
            return self.rho_table_ref.tid
        return self.rho_table

    def Gzx_table(self):
        if self.gzx_table_ref is not None:
            return self.gzx_table_ref.tid
        return self.gzx_table

    def Ax_table(self):
        if self.ax_table_ref is not None:
            return self.ax_table_ref.tid
        return self.ax_table

    def Ath_table(self):
        if self.ath_table_ref is not None:
            return self.ath_table_ref.tid
        return self.ath_table

    def Az_table(self):
        if self.az_table_ref is not None:
            return self.az_table_ref.tid
        return self.az_table

    def Ge_table(self):
        if self.ge_table_ref is not None:
            return self.ge_table_ref.tid
        return self.ge_table

    def raw_fields(self):
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
        return MATT4(mid, k_table=None, cp_table=None, h_table=None,
                     mu_table=None, hgen_table=None, comment='')

    def __init__(self, mid, k_table=None, cp_table=None, h_table=None,
                 mu_table=None, hgen_table=None, comment=''):
        MaterialDependenceThermal.__init__(self)
        if comment:
            self.comment = comment
        if k_table == 0:
            k_table = None
        if cp_table == 0:
            cp_table = None
        if h_table == 0:
            h_table = None
        if mu_table == 0:
            mu_table = None
        if hgen_table == 0:
            hgen_table = None

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
    def add_card(cls, card, comment=''):
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
        k_table = integer_or_blank(card, 2, 'T(K)')
        cp_table = integer_or_blank(card, 3, 'T(CP)')
        h_table = integer_or_blank(card, 5, 'T(H)')
        mu_table = integer_or_blank(card, 6, 'T(mu)')
        hgen_table = integer_or_blank(card, 7, 'T(HGEN)')

        assert len(card) <= 8, 'len(MATT4 card) = %i\ncard=%s' % (len(card), card)
        return MATT4(mid, k_table, cp_table, h_table, mu_table,
                     hgen_table, comment=comment)

    @classmethod
    def add_op2_data(cls, data, comment=''):
        """
        Adds a MATT4 card from the OP2

        Parameters
        ----------
        data : List[varies]
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

        self._xref_table(model, 'k_table', msg=msg)
        self._xref_table(model, 'cp_table', msg=msg)
        self._xref_table(model, 'h_table', msg=msg)
        self._xref_table(model, 'mu_table', msg=msg)
        self._xref_table(model, 'hgen_table', msg=msg)

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

    def K_table(self):
        return self._get_table('k_table')

    def Cp_table(self):
        return self._get_table('cp_table')

    def H_table(self):
        return self._get_table('h_table')

    def Mu_table(self):
        return self._get_table('mu_table')

    def Hgen_table(self):
        return self._get_table('hgen_table')

    def raw_fields(self):
        list_fields = [
            'MATT4', self.Mid(), self.K_table(), self.Cp_table(),
            None,
            self.H_table(), self.Mu_table(), self.Hgen_table()
        ]
        return list_fields

    def repr_fields(self):
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
        return MATT5(mid, kxx_table=None, kxy_table=None, kxz_table=None, kyy_table=None,
                     kyz_table=None, kzz_table=None, cp_table=None, hgen_table=None, comment='')

    def __init__(self, mid, kxx_table=None, kxy_table=None, kxz_table=None,
                 kyy_table=None, kyz_table=None, kzz_table=None,
                 cp_table=None, hgen_table=None, comment=''):
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
    def add_card(cls, card, comment=''):
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
        kxx_table = integer_or_blank(card, 2, 'T(Kxx)')
        kxy_table = integer_or_blank(card, 3, 'T(Kxy)')
        kxz_table = integer_or_blank(card, 5, 'T(Kxz)')
        kyy_table = integer_or_blank(card, 6, 'T(Kyy)')
        kyz_table = integer_or_blank(card, 7, 'T(Kyz)')
        kzz_table = integer_or_blank(card, 8, 'T(Kyz)')
        cp_table = integer_or_blank(card, 9, 'T(Kyz)')
        hgen_table = integer_or_blank(card, 11, 'T(HGEN)')

        assert len(card) <= 12, 'len(MATT5 card) = %i\ncard=%s' % (len(card), card)
        return MATT5(mid, kxx_table, kxy_table, kxz_table, kyy_table,
                     kyz_table, kzz_table, cp_table, hgen_table,
                     comment=comment)

    @classmethod
    def add_op2_data(cls, data, comment=''):
        """
        Adds a MATT5 card from the OP2

        Parameters
        ----------
        data : List[varies]
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

        self._xref_table(model, 'kxx_table', msg=msg)
        self._xref_table(model, 'kxy_table', msg=msg)
        self._xref_table(model, 'kxz_table', msg=msg)
        self._xref_table(model, 'kyy_table', msg=msg)
        self._xref_table(model, 'kyz_table', msg=msg)
        self._xref_table(model, 'kzz_table', msg=msg)
        self._xref_table(model, 'cp_table', msg=msg)
        self._xref_table(model, 'hgen_table', msg=msg)

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

    def Kxx_table(self):
        return self._get_table('kxx_table')

    def Kxy_table(self):
        return self._get_table('kxy_table')

    def Kxz_table(self):
        return self._get_table('kxz_table')

    def Kyy_table(self):
        return self._get_table('kyy_table')

    def Kyz_table(self):
        return self._get_table('kyz_table')

    def Kzz_table(self):
        return self._get_table('kzz_table')

    def Cp_table(self):
        return self._get_table('cp_table')

    def Hgen_table(self):
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
        return MATT8(mid, e1_table=None, e2_table=None, nu12_table=None, g12_table=None,
                     g1z_table=None, g2z_table=None, rho_table=None,
                     a1_table=None, a2_table=None, xt_table=None, xc_table=None,
                     yt_table=None, yc_table=None, s_table=None, ge_table=None,
                     f12_table=None, comment='')

    def __init__(self, mid, e1_table=None, e2_table=None, nu12_table=None,
                 g12_table=None, g1z_table=None, g2z_table=None, rho_table=None,
                 a1_table=None, a2_table=None,
                 xt_table=None, xc_table=None, yt_table=None, yc_table=None,
                 s_table=None, ge_table=None, f12_table=None, comment=''):
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
    def add_card(cls, card, comment=''):
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

        assert len(card) <= 19, 'len(MATT8 card) = %i\ncard=%s' % (len(card), card)
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

        if self.e1_table is not None:
            self.e1_table_ref = model.TableM(self.e1_table)
        if self.e2_table is not None:
            self.e2_table_ref = model.TableM(self.e2_table)
        if self.nu12_table is not None:
            self.nu12_table_ref = model.TableM(self.nu12_table)
        if self.g12_table is not None:
            self.g12_table_ref = model.TableM(self.g12_table)
        if self.g1z_table is not None:
            self.g1z_table_ref = model.TableM(self.g1z_table)
        if self.g2z_table is not None:
            self.g2z_table_ref = model.TableM(self.g2z_table)
        if self.rho_table is not None:
            self.rho_table_ref = model.TableM(self.rho_table)

        if self.a1_table is not None:
            self.a1_table_ref = model.TableM(self.a1_table)
        if self.a2_table is not None:
            self.a2_table_ref = model.TableM(self.a2_table)
        if self.xt_table is not None:
            self.xt_table_ref = model.TableM(self.xt_table)
        if self.xc_table is not None:
            self.xc_table_ref = model.TableM(self.xc_table)
        if self.yt_table is not None:
            self.yt_table_ref = model.TableM(self.yt_table)
        if self.s_table is not None:
            self.s_table_ref = model.TableM(self.s_table)

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

    def E1_table(self):
        if self.e1_table_ref is not None:
            return self.e1_table_ref.tid
        return self.e1_table

    def E2_table(self):
        if self.e2_table_ref is not None:
            return self.e2_table_ref.tid
        return self.e2_table

    def Nu12_table(self):
        if self.nu12_table_ref is not None:
            return self.nu12_table_ref.tid
        return self.nu12_table

    def G12_table(self):
        if self.g12_table_ref is not None:
            return self.g12_table_ref.tid
        return self.g12_table

    def G1z_table(self):
        if self.g1z_table_ref is not None:
            return self.g1z_table_ref.tid
        return self.g1z_table

    def G2z_table(self):
        if self.g2z_table_ref is not None:
            return self.g2z_table_ref.tid
        return self.g2z_table

    def Rho_table(self):
        if self.rho_table_ref is not None:
            return self.rho_table_ref.tid
        return self.rho_table

    def A1_table(self):
        if self.a1_table_ref is not None:
            return self.a1_table_ref.tid
        return self.a1_table

    def A2_table(self):
        if self.a2_table_ref is not None:
            return self.a2_table_ref.tid
        return self.a2_table

    def S_table(self):
        if self.s_table_ref is not None:
            return self.s_table_ref.tid
        return self.s_table

    def Ge_table(self):
        if self.ge_table_ref is not None:
            return self.ge_table_ref.tid
        return self.ge_table

    def F12_table(self):
        if self.f12_table_ref is not None:
            return self.f12_table_ref.tid
        return self.f12_table

    def Xt_table(self):
        if self.xt_table_ref is not None:
            return self.xt_table_ref.tid
        return self.xt_table

    def Xc_table(self):
        if self.xc_table_ref is not None:
            return self.xc_table_ref.tid
        return self.xc_table

    def Yt_table(self):
        if self.yt_table_ref is not None:
            return self.yt_table_ref.tid
        return self.yt_table

    def Yc_table(self):
        if self.yc_table_ref is not None:
            return self.yc_table_ref.tid
        return self.yc_table

    def raw_fields(self):
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
        return MATT9(mid, g11_table=None, g12_table=None, g13_table=None, g14_table=None,
                     g15_table=None, g16_table=None, g22_table=None, g23_table=None,
                     g24_table=None, g25_table=None, g26_table=None, g33_table=None,
                     g34_table=None, g35_table=None, g36_table=None, g44_table=None,
                     g45_table=None, g46_table=None, g55_table=None, g56_table=None,
                     g66_table=None, rho_table=None,
                     a1_table=None, a2_table=None, a3_table=None,
                     a4_table=None, a5_table=None, a6_table=None, ge_table=None, comment='')

    def __init__(self, mid,
                 g11_table=None, g12_table=None, g13_table=None, g14_table=None,
                 g15_table=None, g16_table=None,
                 g22_table=None, g23_table=None, g24_table=None,
                 g25_table=None, g26_table=None,
                 g33_table=None, g34_table=None, g35_table=None, g36_table=None,
                 g44_table=None, g45_table=None, g46_table=None,
                 g55_table=None, g56_table=None,
                 g66_table=None,
                 rho_table=None,
                 a1_table=None, a2_table=None, a3_table=None,
                 a4_table=None, a5_table=None, a6_table=None,
                 ge_table=None,
                 comment=''):
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
    def add_card(cls, card, comment=''):
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
        g11_table = integer_or_blank(card, 2, 'T(G11)')
        g12_table = integer_or_blank(card, 3, 'T(G12)')
        g13_table = integer_or_blank(card, 4, 'T(G13)')
        g14_table = integer_or_blank(card, 5, 'T(G14)')
        g15_table = integer_or_blank(card, 6, 'T(G15)')
        g16_table = integer_or_blank(card, 7, 'T(G16)')

        g22_table = integer_or_blank(card, 8, 'T(G22)')
        g23_table = integer_or_blank(card, 9, 'T(G23)')
        g24_table = integer_or_blank(card, 10, 'T(G24)')
        g25_table = integer_or_blank(card, 11, 'T(G25)')
        g26_table = integer_or_blank(card, 12, 'T(G26)')

        g33_table = integer_or_blank(card, 13, 'T(G33)')
        g34_table = integer_or_blank(card, 14, 'T(G34)')
        g35_table = integer_or_blank(card, 15, 'T(G35)')
        g36_table = integer_or_blank(card, 16, 'T(G36)')

        g44_table = integer_or_blank(card, 17, 'T(G44)')
        g45_table = integer_or_blank(card, 18, 'T(G45)')
        g46_table = integer_or_blank(card, 19, 'T(G46)')

        g55_table = integer_or_blank(card, 20, 'T(G55)')
        g56_table = integer_or_blank(card, 21, 'T(G56)')
        g66_table = integer_or_blank(card, 22, 'T(G66)')

        rho_table = integer_or_blank(card, 23, 'T(RHO)')
        a1_table = integer_or_blank(card, 24, 'T(A1)')
        a2_table = integer_or_blank(card, 25, 'T(A2)')
        a3_table = integer_or_blank(card, 26, 'T(A3)')
        a4_table = integer_or_blank(card, 27, 'T(A4)')
        a5_table = integer_or_blank(card, 28, 'T(A5)')
        a6_table = integer_or_blank(card, 29, 'T(A6)')

        ge_table = integer_or_blank(card, 31, 'T(GE)')

        assert len(card) <= 32, 'len(MATT9 card) = %i\ncard=%s' % (len(card), card)
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
            self.g11_table,
            self.g12_table,
            self.g13_table,
            self.g14_table,
            self.g15_table,
            self.g16_table,

            self.g22_table,
            self.g23_table,
            self.g24_table,
            self.g25_table,
            self.g26_table,

            self.g33_table,
            self.g34_table,
            self.g35_table,
            self.g36_table,

            self.g44_table,
            self.g45_table,
            self.g46_table,

            self.g55_table,
            self.g56_table,
            self.g66_table,

            self.rho_table,
            self.a1_table,
            self.a2_table,
            self.a3_table,
            self.a4_table,
            self.a5_table,
            self.a6_table,
            self.ge_table,
        ]
        return list_fields

    def write_card(self, size: int=8, is_double: bool=False) -> str:
        list_fields = self.raw_fields()
        return self.comment + print_card_8(list_fields)
