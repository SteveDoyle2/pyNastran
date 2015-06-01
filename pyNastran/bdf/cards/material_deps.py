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
from six import integer_types

from pyNastran.bdf.cards.baseCard import BaseCard
from pyNastran.bdf.cards.bdf_tables import Table
from pyNastran.bdf.bdfInterface.assign_type import (integer, integer_or_blank,
    double, double_or_blank, string, blank)
from pyNastran.bdf.field_writer_8 import print_card_8
from pyNastran.bdf.field_writer_16 import print_card_16


class MaterialDependence(BaseCard):
    def __init__(self, card, data):
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

    def __init__(self, card=None, data=None, comment=''):
        MaterialDependence.__init__(self, card, data)
        if comment:
            self._comment = comment

        if card:
            #: Identification number of a MAT1, MAT2, or MAT9 entry.
            self.mid = integer(card, 1, 'mid')
            #: Identification number of a TABLES1 or TABLEST entry. If H is
            #: given, then this field must be blank.
            self.tid = integer_or_blank(card, 2, 'tid')
            #: Type of material nonlinearity. ('NLELAST' for nonlinear elastic
            #: or 'PLASTIC' for elastoplastic.)
            self.Type = string(card, 3, 'Type')

            if self.Type not in ['NLELAST', 'PLASTIC']:
                raise ValueError('MATS1 Type must be [NLELAST, PLASTIC]; Type=%r' % self.Type)
            if self.Type == 'NLELAST':
                # should we even read these?
                self.h = None
                self.hr = None
                self.yf = None
                self.limit1 = None
                self.limit2 = None
                #self.h = blank(card, 4, 'h')
                #self.hr = blank(card, 6, 'hr')
                #self.yf = blank(card, 5, 'yf')
                #self.limit1 = blank(card, 7, 'yf')
                #self.limit2 = blank(card, 8, 'yf')
            else:
                #: Work hardening slope (slope of stress versus plastic strain)
                #: in units of stress. For elastic-perfectly plastic cases,
                #: H=0.0.  For more than a single slope in the plastic range,
                #: the stress-strain data must be supplied on a TABLES1 entry
                #: referenced by TID, and this field must be blank
                self.h = double_or_blank(card, 4, 'H')

                #: Yield function criterion, selected by one of the following
                #: values (1) Von Mises (2) Tresca (3) Mohr-Coulomb
                #: (4) Drucker-Prager
                self.yf = integer_or_blank(card, 5, 'yf', 1)

                #: Hardening Rule, selected by one of the following values
                #: (Integer): (1) Isotropic (Default) (2) Kinematic
                #: (3) Combined isotropic and kinematic hardening
                self.hr = integer_or_blank(card, 6, 'hr', 1)
                #: Initial yield point
                self.limit1 = double(card, 7, 'limit1')

                if self.yf == 3 or self.yf == 4:
                    #: Internal friction angle, measured in degrees, for the
                    #: Mohr-Coulomb and Drucker-Prager yield criteria
                    self.limit2 = double(card, 8, 'limit2')
                else:
                    #self.limit2 = blank(card, 8, 'limit2')
                    self.limit2 = None
            assert len(card) <= 9, 'len(MATS1 card) = %i' % len(card)
        else:
            (mid, tid, Type, h, yf, hr, limit1, limit2) = data
            self.mid = mid
            self.tid = tid
            if Type == 1:
                self.Type = 'NLELAST'
            elif Type == 2:
                self.Type = 'PLASTIC'
            else:
                raise RuntimeError('Invalid Type:  Type=%s; must be 1=NLELAST '
                                   'or 2=PLASTIC' % (Type))
            self.h = h
            self.yf = yf
            self.hr = hr
            self.limit1 = limit1
            self.limit2 = limit2

    def Yf(self):
        d = {1: 'VonMises', 2: 'Tresca', 3: 'MohrCoulomb', 4: 'Drucker-Prager'}
        return d[self.yf]

    def Hf(self):
        d = {1: 'Isotropic', 2: 'Kinematic', 3: 'Combined'}
        return d[self.hr]

    def E(self, strain):
        """
        Gets E (Young's Modulus) for a given strain.

        :param self:   the object pointer
        :param strain: the strain (None -> linear E value)
        :returns E:    Young's Modulus
        """
        msg = "E (Young's Modulus) not implemented for MATS1"
        raise NotImplementedError(msg)
        #if self.tid:
            #E = self.tid.Value(strain)
        #return E

    def cross_reference(self, model):
        msg = 'which is required by MATS1 mid=%s' % self.mid
        self.mid = model.Material(self.mid, msg=msg)
        if self.tid:  # then self.h is used
            self.tid = model.Table(self.tid, msg=msg)

    def Tid(self):
        if isinstance(self.tid, Table):
            return self.tid.tid
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
            return self.comment() + print_card_8(card)
        return self.comment() + print_card_16(card)


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

    def __init__(self, card=None, data=None, comment=''):
        MaterialDependence.__init__(self, card, data)
        if comment:
            self._comment = comment

        if card:
            self.mid = integer(card, 1, 'mid')
            self._E_table = integer_or_blank(card, 2, 'T(E)')
            self._G_table = integer_or_blank(card, 3, 'T(G)')
            self._nu_table = integer_or_blank(card, 4, 'T(nu)')
            self._rho_table = integer_or_blank(card, 5, 'T(rho)')
            self._A_table = integer_or_blank(card, 6, 'T(A)')
            self._ge_table = integer_or_blank(card, 8, 'T(ge)')
            self._st_table = integer_or_blank(card, 9, 'T(st)')
            self._sc_table = integer_or_blank(card, 10, 'T(sc)')
            self._ss_table = integer_or_blank(card, 11, 'T(ss)')

            assert len(card) <= 11, 'len(MATT1 card) = %i' % len(card)
        else:
            raise NotImplementedError()

    def E(self, temperature):
        """
        Gets E (Young's Modulus) for a given temperature.

        :param self:   the object pointer
        :param temperature: the temperature (None -> linear E value)
        :returns E:    Young's Modulus
        """
        E = None
        if self._E_table:
            E = self._E_table.Value(temperature)
        return E

    def cross_reference(self, model):
        msg = 'which is required by MATT1 mid=%s' % self.mid
        self.mid = model.Material(self.mid, msg=msg)
        self._xref_table(model, '_E_table', msg=msg)
        self._xref_table(model, '_G_table', msg=msg)
        self._xref_table(model, '_nu_table', msg=msg)
        self._xref_table(model, '_rho_table', msg=msg)
        self._xref_table(model, '_A_table', msg=msg)
        self._xref_table(model, '_ge_table', msg=msg)
        self._xref_table(model, '_st_table', msg=msg)
        self._xref_table(model, '_sc_table', msg=msg)
        self._xref_table(model, '_ss_table', msg=msg)

    def _xref_table(self, model, key, msg):
        slot = getattr(self, key)
        if slot is not None:
            setattr(self, key, model.Table(slot, msg))

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
            return self.comment() + print_card_8(card)
        return self.comment() + print_card_16(card)

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

    def __init__(self, card=None, data=None, comment=''):
        MaterialDependence.__init__(self, card, data)
        if comment:
            self._comment = comment

        if card:
            self.mid = integer(card, 1, 'mid')
            self._G11_table = integer_or_blank(card, 2, 'T(G11)')
            self._G12_table = integer_or_blank(card, 3, 'T(G12)')
            self._G13_table = integer_or_blank(card, 4, 'T(G13)')
            self._G22_table = integer_or_blank(card, 5, 'T(G22)')
            self._G23_table = integer_or_blank(card, 6, 'T(G23)')
            self._G33_table = integer_or_blank(card, 7, 'T(G33)')
            self._rho_table = integer_or_blank(card, 8, 'T(rho)')
            self._A1_table = integer_or_blank(card, 9, 'T(A1)')
            self._A2_table = integer_or_blank(card, 10, 'T(A2)')
            self._A3_table = integer_or_blank(card, 11, 'T(A3)')
            self._ge_table = integer_or_blank(card, 13, 'T(ge)')
            self._st_table = integer_or_blank(card, 14, 'T(st)')
            self._sc_table = integer_or_blank(card, 15, 'T(sc)')
            self._ss_table = integer_or_blank(card, 16, 'T(ss)')

            assert len(card) <= 17, 'len(MATT2 card) = %i' % len(card)
        else:
            raise NotImplementedError()

    def cross_reference(self, model):
        msg = 'which is required by MATT2 mid=%s' % self.mid
        self.mid = model.Material(self.mid, msg=msg)
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

    def _xref_table(self, model, key, msg):
        slot = getattr(self, key)
        if slot is not None:
            setattr(self, key, model.Table(slot, msg))

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
            return self.comment() + print_card_8(card)
        return self.comment() + print_card_16(card)

#MATT3 - CTRIAX6 only

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

    def __init__(self, card=None, data=None, comment=''):
        MaterialDependence.__init__(self, card, data)
        if comment:
            self._comment = comment

        if card:
            self.mid = integer(card, 1, 'mid')
            self._k_table = integer_or_blank(card, 2, 'T(K)')
            self._cp_table = integer_or_blank(card, 3, 'T(CP)')
            self._H_table = integer_or_blank(card, 5, 'T(H)')
            self._mu_table = integer_or_blank(card, 6, 'T(mu)')
            self._Hgen_table = integer_or_blank(card, 7, 'T(HGEN)')

            assert len(card) <= 8, 'len(MATT4 card) = %i' % len(card)
        else:
            raise NotImplementedError()

    def cross_reference(self, model):
        msg = 'which is required by MATT4 mid=%s' % self.mid
        self.mid = model.Material(self.mid, msg=msg)
        self._xref_table(model, '_k_table', msg=msg)
        self._xref_table(model, '_cp_table', msg=msg)
        self._xref_table(model, '_H_table', msg=msg)
        self._xref_table(model, '_mu_table', msg=msg)
        self._xref_table(model, '_Hgen_table', msg=msg)

    def _xref_table(self, model, key, msg):
        slot = getattr(self, key)
        if slot is not None:
            setattr(self, key, model.Table(slot, msg))

    def K_table(self):
        return self._get_table('_k_table')

    def Cp_table(self):
        return self._get_table('_cp_table')

    def H_table(self):
        return self._get_table('_H_table')

    def mu_table(self):
        return self._get_table('_mu_table')

    def Hgen_table(self):
        return self._get_table('_Hgen_table')

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
            return self.comment() + print_card_8(card)
        return self.comment() + print_card_16(card)

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

    def __init__(self, card=None, data=None, comment=''):
        MaterialDependence.__init__(self, card, data)
        if comment:
            self._comment = comment

        if card:
            self.mid = integer(card, 1, 'mid')
            self._kxx_table = integer_or_blank(card, 2, 'T(Kxx)')
            self._kxy_table = integer_or_blank(card, 3, 'T(Kxy)')
            self._kxz_table = integer_or_blank(card, 5, 'T(Kxz)')
            self._kyy_table = integer_or_blank(card, 6, 'T(Kyy)')
            self._kyz_table = integer_or_blank(card, 7, 'T(Kyz)')
            self._kzz_table = integer_or_blank(card, 8, 'T(Kyz)')
            self._cp_table = integer_or_blank(card, 9, 'T(Kyz)')
            self._hgen_table = integer_or_blank(card, 11, 'T(HGEN)')

            assert len(card) <= 12, 'len(MATT5 card) = %i' % len(card)
        else:
            raise NotImplementedError()

    def cross_reference(self, model):
        msg = 'which is required by MATT5 mid=%s' % self.mid
        self.mid = model.Material(self.mid, msg=msg)
        self._xref_table(model, '_kxx_table', msg=msg)
        self._xref_table(model, '_kxy_table', msg=msg)
        self._xref_table(model, '_kxz_table', msg=msg)
        self._xref_table(model, '_kyy_table', msg=msg)
        self._xref_table(model, '_kyz_table', msg=msg)
        self._xref_table(model, '_kzz_table', msg=msg)
        self._xref_table(model, '_cp_table', msg=msg)
        self._xref_table(model, '_hgen_table', msg=msg)

    def _xref_table(self, model, key, msg):
        slot = getattr(self, key)
        if slot is not None:
            setattr(self, key, model.Table(slot, msg))

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
            return self.comment() + print_card_8(card)
        return self.comment() + print_card_16(card)

#MATT8
#MATT9

class MATT8(MaterialDependence):
    """
    Specifies temperature-dependent material properties on MAT2 entry
    fields via TABLEMi entries.

    +-------+--------+--------+-------+---------+--------+--------+--------+--------+
    |   1   |   2    |   3    |   4   |    5    |   6    |   7    |    8   |   9    |
    +=======+========+========+=======+=========+========+========+========+========+
    | MATT8 |  MID   | T(E1)  | T(E2) | T(Nu12) | T(G12) | T(G1z) | T(G2z) | T(RHO) |
    +-------+--------+--------+-------+---------+--------+--------+--------+--------+
    |       |  T(A1) | T(A2)  |       |  T(Xt)  | T(Yc)  | T(Yt)  | T(Yc)  | T(S)   |
    +-------+--------+--------+-------+---------+--------+--------+--------+--------+
    |       |  T(GE) | T(F12) |       |         |        |        |        |        |
    +-------+--------+--------+-------+---------+--------+--------+--------+--------+
    """
    type = 'MATT8'

    def __init__(self, card=None, data=None, comment=''):
        MaterialDependence.__init__(self, card, data)
        if comment:
            self._comment = comment

        if card:
            self.mid = integer(card, 1, 'mid')
            self._E1_table = integer_or_blank(card, 2, 'T(E1)')
            self._E2_table = integer_or_blank(card, 3, 'T(E2)')
            self._Nu12_table = integer_or_blank(card, 3, 'T(Nu12)')
            self._G12_table = integer_or_blank(card, 5, 'T(G12)')
            self._G1z_table = integer_or_blank(card, 6, 'T(G1z)')
            self._G2z_table = integer_or_blank(card, 7, 'T(G2z)')
            self._rho_table = integer_or_blank(card, 8, 'T(Rho)')
            self._A1_table = integer_or_blank(card, 9, 'T(A1)')
            self._A2_table = integer_or_blank(card, 10, 'T(A2)')

            self._Xt_table = integer_or_blank(card, 12, 'T(Xt)')
            self._Xc_table = integer_or_blank(card, 13, 'T(Xc)')
            self._Yt_table = integer_or_blank(card, 14, 'T(Yt)')
            self._Yc_table = integer_or_blank(card, 15, 'T(Yc)')
            self._S_table = integer_or_blank(card, 16, 'T(S)')
            self._GE_table = integer_or_blank(card, 17, 'T(GE)')
            self._F12_table = integer_or_blank(card, 18, 'T(F12)')

            assert len(card) <= 19, 'len(MATT8 card) = %i' % len(card)
        else:
            raise NotImplementedError()

