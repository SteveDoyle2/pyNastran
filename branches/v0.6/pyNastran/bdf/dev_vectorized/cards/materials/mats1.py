from numpy import array

from pyNastran.bdf.fieldWriter import print_card
from pyNastran.bdf.bdfInterface.assign_type import (integer, integer_or_blank,
    double, double_or_blank, blank, integer_or_string, string)


class MATS1(object):
    """
    Specifies stress-dependent material properties for use in applications
    involving nonlinear materials. This entry is used if a MAT1, MAT2 or MAT9
    entry is specified with the same MID in a nonlinear solution sequence
    (SOLs 106 and 129).
    """
    type = 'MATS1'

    def __init__(self, model):
        self.model = model
        self.material_id = None
        self.table_id = None
        self.Type = None
        self.h = None
        self.hr = None
        self.yf = None
        self.limit1 = None
        self.limit2 = None
        self.n = 1

    def add(self, card=None, comment=''):
        if comment:
            self._comment = comment

        #: Identification number of a MAT1, MAT2, or MAT9 entry.
        self.material_id = integer(card, 1, 'mid')
        #: Identification number of a TABLES1 or TABLEST entry. If H is
        #: given, then this field must be blank.
        self.table_id = integer_or_blank(card, 2, 'tid')
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
            if self.h is None:
                self.hflag = False
            else:
                self.hflag = True

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

    def build(self):
        pass

    def Yf(self):
        d = {1: 'VonMises', 2: 'Tresca', 3: 'MohrCoulomb', 4: 'Drucker-Prager'}
        return d[self.yf]

    def Hf(self):
        d = {1: 'Isotropic', 2: 'Kinematic', 3: 'Combined'}
        return d[self.hr]

    def E(self, stress=None):
        """
        Gets E (Young's Modulus) for a given stress.

        :param self:   the object pointer
        :param stress: the stress (None -> linear E value)
        :returns E:    Young's Modulus
        """
        msg = "E (Young's Modulus) not implemented for MATS1"
        #raise NotImplementedError(msg)
        if self.h is None:
            if self.table_id:
                E = self.tables[table_id].Value(strain)
            else:
                if stress <= self.limit1:
                    return self.materials.mat1[mid].E()
                else:
                    return self.h
                #else:
                    yf = self.Yf()
                    if yf == 'VonMises':
                        E = None
                    elif yf == 'Tresca':
                        pass
                    elif yf == 'MohrCoulomb':
                        self.limit2
                    elif yf == 'Drucker-Prager':
                        self.limit2
                    else:
                        raise NotImplementedError()
                    
                E = self.h
        return E

    def rawFields(self):
        list_fields = ['MATS1', self.Mid(), self.Tid(), self.Type,
                  self.h, self.yf, self.hr, self.limit1, self.limit2]
        return list_fields

    def write_bdf(self, f, size=8, material_ids=None):
        card = ['MATS1', self.material_id, self.table_id, self.Type,
                self.h, self.yf, self.hr, self.limit1, self.limit2]
        f.write(print_card(card))
        
    def reprFields(self):
        return self.rawFields()
