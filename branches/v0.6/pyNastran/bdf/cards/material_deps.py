from __future__ import (nested_scopes, generators, division, absolute_import,
                        print_function, unicode_literals)
#from numpy import zeros, array

#from pyNastran.bdf.fieldWriter import set_blank_if_default
from pyNastran.bdf.cards.baseCard import BaseCard
from pyNastran.bdf.cards.tables import Table
from pyNastran.bdf.bdfInterface.assign_type import (integer, integer_or_blank)

#from pyNastran.bdf.bdfInterface.assign_type import (integer, integer_or_blank,
#    double, double_or_blank,
#    string, string_or_blank, blank)


class MaterialDependence(BaseCard):
    def __init__(self, card, data):
        pass


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

            if self.Type == 'NLELAST':
                self.h = blank(card, 4, 'h')
                self.hr = blank(card, 6, 'hr')
                self.yf = blank(card, 5, 'yf')
                self.limit1 = blank(card, 7, 'yf')
                self.limit2 = blank(card, 8, 'yf')
            else:
                #: Work hardening slope (slope of stress versus plastic strain) in
                #: units of stress. For elastic-perfectly plastic cases, H=0.0.
                #: For more than a single slope in the plastic range, the
                #: stress-strain data must be supplied on a TABLES1 entry
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

    def E(self, strain=None):
        """
        Gets E (Young's Modulus) for a given strain.

        :param self:   the object pointer
        :param strain: the strain (None -> linear E value)
        :returns E:    Young's Modulus
        """
        msg = "E (Young's Modulus) not implemented for MATS1"
        raise NotImplementedError(msg)
        if self.tid:
            E = self.tid.Value(strain)
        return E

    def cross_reference(self, model):
        self.mid = model.Material(self.mid)
        if self.tid:  # then self.h is used
            self.tid = model.Table(self.tid)

    def Mid(self):
        if isinstance(self.mid, int):
            return self.mid
        return self.mid.mid

    def Tid(self):
        if isinstance(self.tid, Table):
            return self.tid.tid
        return self.tid

    def rawFields(self):
        list_fields = ['MATS1', self.Mid(), self.Tid(), self.Type,
                  self.h, self.yf, self.hr, self.limit1, self.limit2]
        return list_fields

    def reprFields(self):
        return self.rawFields()


class MATT1(MaterialDependence):
    """
    Specifies temperature-dependent material properties on MAT1 entry
    fields via TABLEMi entries.
    """
    type = 'MATT1'

    def __init__(self, card=None, data=None, comment=''):
        MaterialDependence.__init__(self, card, data)
        if comment:
            self._comment = comment

        if card:
            self.mid = integer(card, 1, 'mid')
            self._E_table   = integer_or_blank(card, 2, 'T(E)')
            self._G_table   = integer_or_blank(card, 3, 'T(E)')
            self._nu_table  = integer_or_blank(card, 4, 'T(E)')
            self._rho_table = integer_or_blank(card, 5, 'T(rho)')
            self._A_table   = integer_or_blank(card, 6, 'T(A)')
            self._ge_table  = integer_or_blank(card, 8, 'T(ge)')
            self._st_table  = integer_or_blank(card, 9, 'T(st)')
            self._sc_table  = integer_or_blank(card, 10, 'T(sc)')
            self._ss_table  = integer_or_blank(card, 11, 'T(ss)')

            assert len(card) <= 11, 'len(MATT1 card) = %i' % len(card)
        else:
            raise NotImplementedError()

    def E(self, temperature):
        """
        Gets E (Young's Modulus) for a given strain.

        :param self:   the object pointer
        :param strain: the strain (None -> linear E value)
        :returns E:    Young's Modulus
        """
        aaa
        msg = "E (Young's Modulus) not implemented for MATS1"
        raise NotImplementedError(msg)
        if self.tid:
            E = self.tid.Value(strain)
        return E

    def cross_reference(self, model):
        self.mid = model.Material(self.mid)
        self._xref_table(model, '_E_table')
        self._xref_table(model, '_G_table')
        self._xref_table(model, '_nu_table')
        self._xref_table(model, '_rho_table')
        self._xref_table(model, '_A_table')
        self._xref_table(model, '_ge_table')
        self._xref_table(model, '_st_table')
        self._xref_table(model, '_sc_table')
        self._xref_table(model, '_ss_table')
    
    def _xref_table(self, model, key):
        slot = getattr(self, key)
        if slot is not None:
            setattr(self, key, model.Table(slot))

    def _get_table(self, key):
        slot = getattr(self, key)
        if slot is None or isinstance(slot, int):
            return slot
        return slot.tid

    def Mid(self):
        if isinstance(self.mid, int):
            return self.mid
        return self.mid.mid

    #def Tid(self):
        #if isinstance(self.tid, Table):
            #return self.tid.tid
        #return self.tid

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

    def rawFields(self):
        list_fields = ['MATT1', self.Mid(), self.E_table(), self.G_table(),
            self.nu_table(), self.rho_table(), self.A_table(), self.ge_table(),
            self.st_table(), self.sc_table(), self.ss_table(),
        ]
        return list_fields

    def reprFields(self):
        return self.rawFields()