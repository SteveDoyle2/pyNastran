#pylint: disable=E1103,C0103,C0111
from __future__ import (nested_scopes, generators, division, absolute_import,
                        print_function, unicode_literals)
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
from pyNastran.utils import integer_types
from pyNastran.bdf.cards.base_card import BaseCard
from pyNastran.bdf.cards.bdf_tables import Table
from pyNastran.bdf.bdf_interface.assign_type import (
    integer, integer_or_blank, double, double_or_blank, string)
from pyNastran.bdf.field_writer_8 import print_card_8
from pyNastran.bdf.field_writer_16 import print_card_16


class MaterialDependence(BaseCard):
    def __init__(self):
        self.mid = None

    def Mid(self):
        if isinstance(self.mid, integer_types):
            return self.mid
        return self.mid.mid  # TODO: is this something that should be supported?

    def _get_table(self, key):
        """internal method for accessing tables"""
        table = getattr(self, key)
        if table is None or isinstance(table, integer_types):
            return table
        return table.tid


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

    def validate(self):
        if self.Type not in ['NLELAST', 'PLASTIC']:
            raise ValueError('MATS1 Type must be [NLELAST, PLASTIC]; Type=%r' % self.Type)

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

        if Type not in ['NLELAST', 'PLASTIC']:
            raise ValueError('MATS1 Type must be [NLELAST, PLASTIC]; Type=%r' % Type)
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
        else:
            raise RuntimeError('Invalid Type:  Type=%s; must be 1=NLELAST '
                               'or 2=PLASTIC' % (Type))
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
            #E = self.tid.Value(strain)
        #return E

    def cross_reference(self, model):
        """
        Cross links the card so referenced cards can be extracted directly

        Parameters
        ----------
        model : BDF()
            the BDF object
        """
        msg = ', which is required by MATS1 mid=%s' % self.mid
        self.mid = model.Material(self.mid, msg=msg)
        self.mid_ref = self.mid
        if self.tid:  # then self.h is used
            self.tid = model.Table(self.tid, msg=msg) # TABLES1 or TABLEST
            self.tid_ref = self.tid

    def uncross_reference(self):
        self.mid = self.Mid()
        if self.tid:
            self.tid = self.Tid()
            del self.tid_ref
        del self.mid_ref

    def Tid(self):
        if isinstance(self.tid, Table):
            return self.tid_ref.tid
        return self.tid

    def raw_fields(self):
        list_fields = ['MATS1', self.Mid(), self.Tid(), self.Type,
                       self.h, self.yf, self.hr, self.limit1, self.limit2]
        return list_fields

    def repr_fields(self):
        return self.raw_fields()

    def write_card(self, size=8, is_double=False):
        card = self.repr_fields()
        if size == 8:
            return self.comment + print_card_8(card)
        return self.comment + print_card_16(card)


class MATT1(MaterialDependence):
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

    def __init__(self, mid, E_table=None, G_table=None, nu_table=None,
                 rho_table=None, A_table=None, ge_table=None, st_table=None,
                 sc_table=None, ss_table=None, comment=''):
        MaterialDependence.__init__(self)
        if comment:
            self.comment = comment
        self.mid = mid
        if E_table == 0:
            E_table = None
        if G_table == 0:
            G_table = None
        if nu_table == 0:
            nu_table = None
        if rho_table == 0:
            rho_table = None
        if A_table == 0:
            A_table = None
        if ge_table == 0:
            ge_table = None
        if st_table == 0:
            st_table = None
        if sc_table == 0:
            sc_table = None
        if ss_table == 0:
            ss_table = None
        self._E_table = E_table
        self._G_table = G_table
        self._nu_table = nu_table
        self._rho_table = rho_table
        self._A_table = A_table
        self._ge_table = ge_table
        self._st_table = st_table
        self._sc_table = sc_table
        self._ss_table = ss_table

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
        E_table = integer_or_blank(card, 2, 'T(E)')
        G_table = integer_or_blank(card, 3, 'T(G)')
        nu_table = integer_or_blank(card, 4, 'T(nu)')
        rho_table = integer_or_blank(card, 5, 'T(rho)')
        A_table = integer_or_blank(card, 6, 'T(A)')
        ge_table = integer_or_blank(card, 8, 'T(ge)')
        st_table = integer_or_blank(card, 9, 'T(st)')
        sc_table = integer_or_blank(card, 10, 'T(sc)')
        ss_table = integer_or_blank(card, 11, 'T(ss)')

        assert len(card) <= 12, 'len(MATT1 card) = %i\ncard=%s' % (len(card), card)
        return MATT1(mid, E_table, G_table, nu_table, rho_table, A_table,
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
        if self._E_table:
            E = self._E_table.Value(temperature)
        return E

    def cross_reference(self, model):
        """
        Cross links the card so referenced cards can be extracted directly

        Parameters
        ----------
        model : BDF()
            the BDF object
        """
        msg = ', which is required by MATT1 mid=%s' % self.mid
        self.mid = model.Material(self.mid, msg=msg)
        self.mid_ref = self.mid

        ## TODO: add refs
        self._xref_table(model, '_E_table', msg=msg)
        self._xref_table(model, '_G_table', msg=msg)
        self._xref_table(model, '_nu_table', msg=msg)
        self._xref_table(model, '_rho_table', msg=msg)
        self._xref_table(model, '_A_table', msg=msg)
        self._xref_table(model, '_ge_table', msg=msg)
        self._xref_table(model, '_st_table', msg=msg)
        self._xref_table(model, '_sc_table', msg=msg)
        self._xref_table(model, '_ss_table', msg=msg)

    def uncross_reference(self):
        self.mid = self.Mid()
        ## TODO: remove refs

        self._E_table = self.E_table()
        self._G_table = self.G_table()
        self._nu_table = self.nu_table()
        self._rho_table = self.rho_table()
        self._A_table = self.A_table()
        self._ge_table = self.ge_table()
        self._st_table = self.st_table()
        self._sc_table = self.sc_table()
        self._ss_table = self.ss_table()
        del self.mid_ref

    def _xref_table(self, model, key, msg):
        slot = getattr(self, key)
        if slot is not None:
            setattr(self, key, model.TableM(slot, msg))

    def E_table(self):
        return self._get_table('_E_table')

    def G_table(self):
        return self._get_table('_G_table')

    def nu_table(self):
        return self._get_table('_nu_table')

    def rho_table(self):
        return self._get_table('_rho_table')

    def A_table(self):
        return self._get_table('_A_table')

    def ge_table(self):
        return self._get_table('_ge_table')

    def st_table(self):
        return self._get_table('_st_table')

    def sc_table(self):
        return self._get_table('_sc_table')

    def ss_table(self):
        return self._get_table('_ss_table')

    def raw_fields(self):
        list_fields = [
            'MATT1', self.Mid(), self.E_table(), self.G_table(),
            self.nu_table(), self.rho_table(), self.A_table(), self.ge_table(),
            self.st_table(), self.sc_table(), self.ss_table(),
        ]
        return list_fields

    def repr_fields(self):
        return self.raw_fields()

    def write_card(self, size=8, is_double=False):
        card = self.repr_fields()
        if size == 8:
            return self.comment + print_card_8(card)
        return self.comment + print_card_16(card)

class MATT2(MaterialDependence):
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

    def __init__(self, mid, G11_table=None, G12_table=None, G13_table=None,
                 G22_table=None, G23_table=None, G33_table=None, rho_table=None,
                 A1_table=None, A2_table=None, A3_table=None,
                 ge_table=None, st_table=None, sc_table=None, ss_table=None, comment=''):
        MaterialDependence.__init__(self)
        if comment:
            self.comment = comment

        self.mid = mid
        self._G11_table = G11_table
        self._G12_table = G12_table
        self._G13_table = G13_table
        self._G22_table = G22_table
        self._G23_table = G23_table
        self._G33_table = G33_table
        self._rho_table = rho_table
        self._A1_table = A1_table
        self._A2_table = A2_table
        self._A3_table = A3_table
        self._ge_table = ge_table
        self._st_table = st_table
        self._sc_table = sc_table
        self._ss_table = ss_table

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
        G11_table = integer_or_blank(card, 2, 'T(G11)')
        G12_table = integer_or_blank(card, 3, 'T(G12)')
        G13_table = integer_or_blank(card, 4, 'T(G13)')
        G22_table = integer_or_blank(card, 5, 'T(G22)')
        G23_table = integer_or_blank(card, 6, 'T(G23)')
        G33_table = integer_or_blank(card, 7, 'T(G33)')
        rho_table = integer_or_blank(card, 8, 'T(rho)')
        A1_table = integer_or_blank(card, 9, 'T(A1)')
        A2_table = integer_or_blank(card, 10, 'T(A2)')
        A3_table = integer_or_blank(card, 11, 'T(A3)')
        ge_table = integer_or_blank(card, 13, 'T(ge)')
        st_table = integer_or_blank(card, 14, 'T(st)')
        sc_table = integer_or_blank(card, 15, 'T(sc)')
        ss_table = integer_or_blank(card, 16, 'T(ss)')

        assert len(card) <= 17, 'len(MATT2 card) = %i\ncard=%s' % (len(card), card)
        return MATT2(mid, G11_table, G12_table, G13_table, G22_table, G23_table,
                     G33_table, rho_table, A1_table,
                     A2_table, A3_table, ge_table,
                     st_table, sc_table, ss_table,
                     comment=comment)

    def cross_reference(self, model):
        """
        Cross links the card so referenced cards can be extracted directly

        Parameters
        ----------
        model : BDF()
            the BDF object
        """
        msg = ', which is required by MATT2 mid=%s' % self.mid
        self.mid = model.Material(self.mid, msg=msg)
        self.mid_ref = self.mid

        ## TODO: add refs
        self._xref_table(model, '_G11_table', msg=msg)
        self._xref_table(model, '_G12_table', msg=msg)
        self._xref_table(model, '_G13_table', msg=msg)
        self._xref_table(model, '_G22_table', msg=msg)
        self._xref_table(model, '_G23_table', msg=msg)
        self._xref_table(model, '_G33_table', msg=msg)
        self._xref_table(model, '_rho_table', msg=msg)
        self._xref_table(model, '_A1_table', msg=msg)
        self._xref_table(model, '_A2_table', msg=msg)
        self._xref_table(model, '_A3_table', msg=msg)
        self._xref_table(model, '_ge_table', msg=msg)
        self._xref_table(model, '_st_table', msg=msg)
        self._xref_table(model, '_sc_table', msg=msg)
        self._xref_table(model, '_ss_table', msg=msg)

    def uncross_reference(self):
        self.mid = self.Mid()
        self._G11_table = self.G11_table()
        self._G12_table = self.G12_table()
        self._G13_table = self.G13_table()
        self._G22_table = self.G22_table()
        self._G23_table = self.G23_table()
        self._G33_table = self.G33_table()
        self._rho_table = self.rho_table()
        self._A1_table = self.A1_table()
        self._A2_table = self.A2_table()
        self._A3_table = self.A3_table()
        self._ge_table = self.ge_table()
        self._st_table = self.st_table()
        self._sc_table = self.sc_table()
        self._ss_table = self.ss_table()
        del self.mid_ref

    def _xref_table(self, model, key, msg):
        slot = getattr(self, key)
        if slot is not None:
            setattr(self, key, model.TableM(slot, msg))

    def G11_table(self):
        return self._get_table('_G11_table')

    def G12_table(self):
        return self._get_table('_G12_table')

    def G13_table(self):
        return self._get_table('_G13_table')

    def G22_table(self):
        return self._get_table('_G22_table')

    def G23_table(self):
        return self._get_table('_G23_table')

    def G33_table(self):
        return self._get_table('_G33_table')

    def rho_table(self):
        return self._get_table('_rho_table')

    def A1_table(self):
        return self._get_table('_A1_table')

    def A2_table(self):
        return self._get_table('_A2_table')

    def A3_table(self):
        return self._get_table('_A3_table')

    def ge_table(self):
        return self._get_table('_ge_table')

    def st_table(self):
        return self._get_table('_st_table')

    def sc_table(self):
        return self._get_table('_sc_table')

    def ss_table(self):
        return self._get_table('_ss_table')

    def raw_fields(self):
        list_fields = [
            'MATT2', self.Mid(), self.G11_table(), self.G12_table(),
            self.G13_table(), self.G22_table(), self.G23_table(),
            self.G33_table(), self.rho_table(), self.A1_table(),
            self.A2_table(), self.A3_table(), None, self.ge_table(),
            self.st_table(), self.sc_table(), self.ss_table()
        ]
        return list_fields

    def repr_fields(self):
        return self.raw_fields()

    def write_card(self, size=8, is_double=False):
        card = self.repr_fields()
        if size == 8:
            return self.comment + print_card_8(card)
        return self.comment + print_card_16(card)

#MATT3 - CTRIAX6 only
class MATT3(MaterialDependence):
    """
    Specifies temperature-dependent material properties on MAT3 entry fields via
    TABLEMi entries that are temperature dependent.
    +--------+--------+-------+--------+-------+----------+----------+---------+--------+
    |    1   |    2   |   3   |    4   |   5   |     6    |     7    |    8    |    9   |
    +========+========+=======+========+=======+==========+==========+=========+========+
    | MATT3  |   MID  | T(EX) | T(ETH) | T(EZ) | T(NUXTH) | T(NUTHZ) | T(NUZX) | T(RHO) |
    +--------+--------+-------+--------+-------+----------+----------+---------+--------+
    |        |        |       | T(GZX) | T(AX) |  T(ATH)  |  T(AZ)   |         |  T(GE) |
    +--------+--------+-------+--------+-------+----------+----------+---------+--------+
    """
    type = 'MATT3'

    def __init__(self, mid, ex_table=None, eth_table=None, ez_table=None,
                 nuth_table=None, nuxz_table=None, rho_table=None,
                 gzx_table=None, ax_table=None, ath_table=None, az_table=None,
                 ge_table=None, comment=''):
        """
        Creates a MATT3 card
        """
        MaterialDependence.__init__(self)
        if comment:
            self.comment = comment

        self.mid = mid
        self._ex_table = ex_table
        self._eth_table = eth_table
        self._ez_table = ez_table
        self._nuth_table = nuth_table
        self._nuxz_table = nuxz_table
        self._rho_table = rho_table
        self._gzx_table = gzx_table
        self._ax_table = ax_table
        self._ath_table = ath_table
        self._az_table = az_table
        self._ge_table = ge_table

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

    def cross_reference(self, model):
        msg = ', which is required by MATT3 mid=%s' % self.mid
        self.mid = model.Material(self.mid, msg=msg)
        self.mid_ref = self.mid

        if self._ex_table is not None:
            self.ex_table_ref = model.TableM(self._ex_table)
        if self._eth_table is not None:
            self.eth_table_ref = model.TableM(self._eth_table)
        if self._ez_table is not None:
            self.ez_table_ref = model.TableM(self._ez_table)
        if self._nuth_table is not None:
            self.nuth_table_ref = model.TableM(self._nuth_table)
        if self._nuxz_table is not None:
            self.nuxz_table_ref = model.TableM(self._nuxz_table)
        if self._rho_table is not None:
            self.rho_table_ref = model.TableM(self._rho_table)

        if self._gzx_table is not None:
            self.gzx_table_ref = model.TableM(self._gzx_table)
        if self._ax_table is not None:
            self.ax_table_ref = model.TableM(self._ax_table)
        if self._ath_table is not None:
            self.ath_table_ref = model.TableM(self._ath_table)
        if self._az_table is not None:
            self.az_table_ref = model.TableM(self._az_table)
        if self._ge_table is not None:
            self.ge_table_ref = model.TableM(self._ge_table)



    def uncross_reference(self):
        self.mid = self.Mid()
        self.mid_ref = None

        self._ex_table = self.ex_table()
        self._eth_table = self.eth_table()
        self._ez_table = self.ez_table()
        self._nuth_table = self.nuth_table()
        self._nuxz_table = self.nuxz_table()
        self._rho_table = self.rho_table()
        self._gzx_table = self.gzx_table()
        self._ax_table = self.ax_table()
        self._ath_table = self.ath_table()
        self._az_table = self.az_table()
        self._ge_table = self.ge_table()

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

    def ex_table(self):
        if self.ex_table_ref is not None:
            return self.ex_table_ref.tid
        return self._ex_table

    def eth_table(self):
        if self.eth_table_ref is not None:
            return self.eth_table_ref.tid
        return self._eth_table

    def ez_table(self):
        if self.ez_table_ref is not None:
            return self.ez_table_ref.tid
        return self._eth_table

    def nuth_table(self):
        if self.nuth_table_ref is not None:
            return self.nuth_table_ref.tid
        return self._nuth_table

    def nuxz_table(self):
        if self.nuxz_table_ref is not None:
            return self.nuxz_table_ref.tid
        return self._nuxz_table

    def rho_table(self):
        if self.rho_table_ref is not None:
            return self.rho_table_ref.tid
        return self._rho_table

    def gzx_table(self):
        if self.gzx_table_ref is not None:
            return self.gzx_table_ref.tid
        return self._gzx_table

    def ax_table(self):
        if self.ax_table_ref is not None:
            return self.ax_table_ref.tid
        return self._ax_table

    def ath_table(self):
        if self.ath_table_ref is not None:
            return self.ath_table_ref.tid
        return self._ath_table

    def az_table(self):
        if self.az_table_ref is not None:
            return self.az_table_ref.tid
        return self._az_table

    def ge_table(self):
        if self.ge_table_ref is not None:
            return self.ge_table_ref.tid
        return self._ge_table

    def raw_fields(self):
        list_fields = [
            'MATT3', self.Mid(), self.ex_table(), self.eth_table(), self.ez_table(),
            self.nuth_table(), self.nuxz_table(), self.rho_table(), None, None,
            self.gzx_table(), self.ax_table(), self.ath_table(), self.az_table(),
            None, self.ge_table(),
        ]
        return list_fields

    def repr_fields(self):
        return self.raw_fields()

    def write_card(self, size=8, is_double=False):
        card = self.repr_fields()
        if size == 8:
            return self.comment + print_card_8(card)
        return self.comment + print_card_16(card)

class MATT4(MaterialDependence):
    """
    Specifies temperature-dependent material properties on MAT2 entry
    fields via TABLEMi entries.

    +-------+-------+-------+-------+--------+-------+-------+---------+--------+
    |   1   |   2   |   3   |   4   |   5    |   6   |   7   |    8    |   9    |
    +=======+=======+=======+=======+========+=======+=======+=========+========+
    | MATT4 |  MID  |  T(K) | T(CP) |        | T(H)  | T(mu) | T(HGEN) |        |
    +-------+-------+-------+-------+--------+-------+-------+---------+--------+
    """
    type = 'MATT4'

    def __init__(self, mid, k_table=None, cp_table=None, h_table=None,
                 mu_table=None, hgen_table=None, comment=''):
        MaterialDependence.__init__(self)
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
        self._k_table = k_table
        self._cp_table = cp_table
        self._h_table = h_table
        self._mu_table = mu_table
        self._hgen_table = hgen_table

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

    def cross_reference(self, model):
        """
        Cross links the card so referenced cards can be extracted directly

        Parameters
        ----------
        model : BDF()
            the BDF object
        """
        msg = ', which is required by MATT4 mid=%s' % self.mid
        self.mid = model.Material(self.mid, msg=msg)
        self.mid_ref = self.mid

        ## TODO: add refs
        self._xref_table(model, '_k_table', msg=msg)
        self._xref_table(model, '_cp_table', msg=msg)
        self._xref_table(model, '_h_table', msg=msg)
        self._xref_table(model, '_mu_table', msg=msg)
        self._xref_table(model, '_hgen_table', msg=msg)

        self.mid_ref = self.mid

    def uncross_reference(self):
        self.mid = self.Mid()
        self._k_table = self.K_table()
        self._cp_table = self.Cp_table()
        self._h_table = self.H_table()
        self._mu_table = self.mu_table()
        self._hgen_table = self.Hgen_table()
        del self.mid_ref

    def _xref_table(self, model, key, msg):
        slot = getattr(self, key)
        if slot is not None:
            setattr(self, key, model.TableM(slot, msg))

    def K_table(self):
        return self._get_table('_k_table')

    def Cp_table(self):
        return self._get_table('_cp_table')

    def H_table(self):
        return self._get_table('_h_table')

    def mu_table(self):
        return self._get_table('_mu_table')

    def Hgen_table(self):
        return self._get_table('_hgen_table')

    def raw_fields(self):
        list_fields = [
            'MATT4', self.Mid(), self.K_table(), self.Cp_table(),
            None,
            self.H_table(), self.mu_table(), self.Hgen_table()
        ]
        return list_fields

    def repr_fields(self):
        return self.raw_fields()

    def write_card(self, size=8, is_double=False):
        card = self.repr_fields()
        if size == 8:
            return self.comment + print_card_8(card)
        return self.comment + print_card_16(card)


class MATT5(MaterialDependence):
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

    def __init__(self, mid, kxx_table=None, kxy_table=None, kxz_table=None,
                 kyy_table=None, kyz_table=None, kzz_table=None,
                 cp_table=None, hgen_table=None, comment=''):
        MaterialDependence.__init__(self)
        if comment:
            self.comment = comment
        self.mid = mid
        self._kxx_table = kxx_table
        self._kxy_table = kxy_table
        self._kxz_table = kxz_table
        self._kyy_table = kyy_table
        self._kyz_table = kyz_table
        self._kzz_table = kzz_table
        self._cp_table = cp_table
        self._hgen_table = hgen_table

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

    def cross_reference(self, model):
        """
        Cross links the card so referenced cards can be extracted directly

        Parameters
        ----------
        model : BDF()
            the BDF object
        """
        msg = ', which is required by MATT5 mid=%s' % self.mid
        self.mid = model.Material(self.mid, msg=msg)
        self.mid_ref = self.mid

        ## TODO: add refs
        self._xref_table(model, '_kxx_table', msg=msg)
        self._xref_table(model, '_kxy_table', msg=msg)
        self._xref_table(model, '_kxz_table', msg=msg)
        self._xref_table(model, '_kyy_table', msg=msg)
        self._xref_table(model, '_kyz_table', msg=msg)
        self._xref_table(model, '_kzz_table', msg=msg)
        self._xref_table(model, '_cp_table', msg=msg)
        self._xref_table(model, '_hgen_table', msg=msg)

    def uncross_reference(self):
        self.mid = self.Mid()
        self._kxx_table = self.Kxx_table()
        self._kxy_table = self.Kxy_table()
        self._kxz_table = self.Kxz_table()
        self._kyy_table = self.Kyy_table()
        self._kyz_table = self.Kyz_table()
        self._kzz_table = self.Kzz_table()
        self._cp_table = self.Cp_table()
        self._hgen_table = self.Hgen_table()
        del self.mid_ref

    def _xref_table(self, model, key, msg):
        slot = getattr(self, key)
        if slot is not None:
            setattr(self, key, model.TableM(slot, msg))

    def Kxx_table(self):
        return self._get_table('_kxx_table')

    def Kxy_table(self):
        return self._get_table('_kxy_table')

    def Kxz_table(self):
        return self._get_table('_kxz_table')

    def Kyy_table(self):
        return self._get_table('_kyy_table')

    def Kyz_table(self):
        return self._get_table('_kyz_table')

    def Kzz_table(self):
        return self._get_table('_kzz_table')

    def Cp_table(self):
        return self._get_table('_cp_table')

    def Hgen_table(self):
        return self._get_table('_hgen_table')

    def raw_fields(self):
        list_fields = ['MATT5', self.Mid(),
                       self.Kxx_table(), self.Kxy_table(), self.Kxz_table(),
                       self.Kyy_table(), self.Kyz_table(), self.Kzz_table(),
                       self.Cp_table(), None, self.Hgen_table()]
        return list_fields

    def repr_fields(self):
        return self.raw_fields()

    def write_card(self, size=8, is_double=False):
        card = self.repr_fields()
        if size == 8:
            return self.comment + print_card_8(card)
        return self.comment + print_card_16(card)


class MATT8(MaterialDependence):
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

    def __init__(self, mid, E1_table=None, E2_table=None, Nu12_table=None,
                 G12_table=None, G1z_table=None, G2z_table=None, rho_table=None,
                 a1_table=None, a2_table=None,
                 xt_table=None, xc_table=None, yt_table=None, yc_table=None,
                 s_table=None, ge_table=None, f12_table=None, comment=''):
        MaterialDependence.__init__(self)
        if comment:
            self.comment = comment

        self.mid = mid
        self._E1_table = E1_table
        self._E2_table = E2_table
        self._Nu12_table = Nu12_table
        self._G12_table = G12_table
        self._G1z_table = G1z_table
        self._G2z_table = G2z_table
        self._rho_table = rho_table
        self._a1_table = a1_table
        self._a2_table = a2_table

        self._xt_table = xt_table
        self._xc_table = xc_table
        self._yt_table = yt_table
        self._yc_table = yc_table
        self._s_table = s_table
        self._ge_table = ge_table
        self._f12_table = f12_table

        self.E1_table_ref = None
        self.E2_table_ref = None
        self.Nu12_table_ref = None
        self.G12_table_ref = None
        self.G1z_table_ref = None
        self.G2z_table_ref = None
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
        E1_table = integer_or_blank(card, 2, 'T(E1)')
        E2_table = integer_or_blank(card, 3, 'T(E2)')
        Nu12_table = integer_or_blank(card, 3, 'T(Nu12)')
        G12_table = integer_or_blank(card, 5, 'T(G12)')
        G1z_table = integer_or_blank(card, 6, 'T(G1z)')
        G2z_table = integer_or_blank(card, 7, 'T(G2z)')
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
        return MATT8(mid, E1_table, E2_table, Nu12_table, G12_table,
                     G1z_table, G2z_table, rho_table,
                     a1_table, a2_table, xt_table,
                     xc_table, yt_table, yc_table,
                     s_table, ge_table, f12_table,
                     comment=comment)

    def cross_reference(self, model):
        """
        Cross links the card so referenced cards can be extracted directly

        Parameters
        ----------
        model : BDF()
            the BDF object
        """
        msg = ', which is required by MATT1 mid=%s' % self.mid
        self.mid_ref = model.Material(self.mid, msg=msg)

        if self._E1_table is not None:
            self.E1_table_ref = model.TableM(self._E1_table)
        if self._E2_table is not None:
            self.E2_table_ref = model.TableM(self._E2_table)
        if self._Nu12_table is not None:
            self.Nu12_table_ref = model.TableM(self._Nu12_table)
        if self._G12_table is not None:
            self.G12_table_ref = model.TableM(self._G12_table)
        if self._G1z_table is not None:
            self.G1z_table_ref = model.TableM(self._G1z_table)
        if self._G2z_table is not None:
            self.G2z_table_ref = model.TableM(self._G2z_table)
        if self._rho_table is not None:
            self.rho_table_ref = model.TableM(self._rho_table)

        if self._a1_table is not None:
            self.a1_table_ref = model.TableM(self._a1_table)
        if self._a2_table is not None:
            self.a2_table_ref = model.TableM(self._a2_table)
        if self._xt_table is not None:
            self.xt_table_ref = model.TableM(self._xt_table)
        if self._xc_table is not None:
            self.xc_table_ref = model.TableM(self._xc_table)
        if self._yt_table is not None:
            self.yt_table_ref = model.TableM(self._yt_table)
        if self._s_table is not None:
            self.s_table_ref = model.TableM(self._s_table)

    def uncross_reference(self):
        self._E1_table = self.E1_table()
        self._E2_table = self.E2_table()
        self._Nu12_table = self.Nu12_table()
        self._G12_table = self.G12_table()
        self._G1z_table = self.G1z_table()
        self._G2z_table = self.G2z_table()
        self._rho_table = self.rho_table()
        self._a1_table = self.a1_table()
        self._a2_table = self.a2_table()

        self._xt_table = self.xt_table()
        self._xc_table = self.xc_table()
        self._yt_table = self.yt_table()
        self._yc_table = self.yc_table()
        self._s_table = self.s_table()
        self._ge_table = self.ge_table()
        self._f12_table = self.f12_table()

        self.E1_table_ref = None
        self.E2_table_ref = None
        self.Nu12_table_ref = None
        self.G12_table_ref = None
        self.G1z_table_ref = None
        self.G2z_table_ref = None
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
        if self.E1_table_ref is not None:
            return self.E1_table_ref.tid
        return self._E1_table

    def E2_table(self):
        if self.E2_table_ref is not None:
            return self.E2_table_ref.tid
        return self._E2_table

    def Nu12_table(self):
        if self.Nu12_table_ref is not None:
            return self.Nu12_table_ref.tid
        return self._Nu12_table

    def G12_table(self):
        if self.G12_table_ref is not None:
            return self.G12_table_ref.tid
        return self._G12_table

    def G1z_table(self):
        if self.G1z_table_ref is not None:
            return self.G1z_table_ref.tid
        return self._G1z_table

    def G2z_table(self):
        if self.G2z_table_ref is not None:
            return self.G2z_table_ref.tid
        return self._G2z_table

    def G1z_table(self):
        if self.G1z_table_ref is not None:
            return self.G1z_table_ref.tid
        return self._G1z_table

    def rho_table(self):
        if self.rho_table_ref is not None:
            return self.rho_table_ref.tid
        return self._rho_table

    def a1_table(self):
        if self.a1_table_ref is not None:
            return self.a1_table_ref.tid
        return self._a1_table

    def a2_table(self):
        if self.a2_table_ref is not None:
            return self.a2_table_ref.tid
        return self._a2_table

    def s_table(self):
        if self.s_table_ref is not None:
            return self.s_table_ref.tid
        return self._s_table

    def ge_table(self):
        if self.ge_table_ref is not None:
            return self.ge_table_ref.tid
        return self._ge_table

    def f12_table(self):
        if self.f12_table_ref is not None:
            return self.f12_table_ref.tid
        return self._f12_table

    def xt_table(self):
        if self.xt_table_ref is not None:
            return self.xt_table_ref.tid
        return self._xt_table

    def xc_table(self):
        if self.xc_table_ref is not None:
            return self.xc_table_ref.tid
        return self._xc_table

    def yt_table(self):
        if self.yt_table_ref is not None:
            return self.yt_table_ref.tid
        return self._yt_table

    def yc_table(self):
        if self.yc_table_ref is not None:
            return self.yc_table_ref.tid
        return self._yc_table

    def raw_fields(self):
        list_fields = ['MATT8', self.mid, self.E1_table(), self.E2_table(), self.G12_table(),
                       self.G1z_table(), self.G2z_table(), self.rho_table(),
                       self.a1_table(), self.a2_table(), None,
                       self.xt_table(), self.xc_table(), self.yt_table(), self.yc_table(),
                       self.s_table(), self.ge_table(), self.f12_table()]
        return list_fields

    def write_card(self, size=8, is_double=False):
        """
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
        list_fields = self.raw_fields()
        return self.comment + print_card_8(list_fields)


#MATT9
