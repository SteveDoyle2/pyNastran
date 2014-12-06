from numpy import array, zeros

from pyNastran.bdf.fieldWriter import print_card_8
from pyNastran.bdf.fieldWriter16 import print_card_16
from pyNastran.bdf.bdfInterface.assign_type import (integer, integer_or_blank,
    double, double_or_blank, blank, integer_or_string, string)
from pyNastran.bdf.dev_vectorized.cards.vectorized_card import VectorizedCard
from .mat1 import Material

class MATS1(Material):
    """
    Specifies stress-dependent material properties for use in applications
    involving nonlinear materials. This entry is used if a MAT1, MAT2 or MAT9
    entry is specified with the same MID in a nonlinear solution sequence
    (SOLs 106 and 129).
    """
    type = 'MATS1'

    def __init__(self, model):
        Material.__init__(self, model)
        self.n = 0
        self.material_id = None
        self.table_id = None
        self.Type = None
        self.h = None
        self.hr = None
        self.yf = None
        self.limit1 = None
        self.limit2 = None
        self.i = 0

    def slice_by_index(self, i):
        i = asarray(i)
        obj = MATS1(self.model)
        obj.n = len(i)
        obj.i = self.i
        #obj._cards = self._cards[i]
        #obj._comments = obj._comments[i]
        #obj.comments = obj.comments[i]
        obj.material_id = self.material_id[i]
        obj.table_id = self.table_id[i]
        obj.Type = self.Type[i]
        obj.h = self.h[i]
        obj.hflag = self.hflag[i]
        obj.yf = self.yf[i]
        obj.hr = self.hr[i]
        obj.limit1 = self.limit1[i]
        obj.limit2 = self.limit2[i]
        return obj

    def allocate(self, ncards):
        self.n = ncards
        float_fmt = self.model.float
        #: Identification number of a MAT1, MAT2, or MAT9 entry.
        self.material_id = zeros(ncards, dtype='int32')
        #: Identification number of a TABLES1 or TABLEST entry. If H is
        #: given, then this field must be blank.
        self.table_id = zeros(ncards, dtype='int32')
        #: Type of material nonlinearity. ('NLELAST' for nonlinear elastic
        #: or 'PLASTIC' for elastoplastic.)
        self.Type = zeros(ncards, dtype='|S8')

        #: Work hardening slope (slope of stress versus plastic strain) in
        #: units of stress. For elastic-perfectly plastic cases, H=0.0.
        #: For more than a single slope in the plastic range, the
        #: stress-strain data must be supplied on a TABLES1 entry
        #: referenced by TID, and this field must be blank
        self.h = zeros(ncards, dtype=float_fmt)
        self.hflag = zeros(ncards, dtype='bool')

        #: Yield function criterion, selected by one of the following
        #: values (1) Von Mises (2) Tresca (3) Mohr-Coulomb
        #: (4) Drucker-Prager
        self.yf = zeros(ncards, dtype='int32')

        #: Hardening Rule, selected by one of the following values
        #: (Integer): (1) Isotropic (Default) (2) Kinematic
        #: (3) Combined isotropic and kinematic hardening
        self.hr = zeros(ncards, dtype='int32')
        #: Initial yield point
        self.limit1 = zeros(ncards, dtype=float_fmt)
        self.limit2 = zeros(ncards, dtype=float_fmt)

    def add(self, card=None, comment=''):
        if comment:
            self._comment = comment
        i = self.i
        #: Identification number of a MAT1, MAT2, or MAT9 entry.
        self.material_id[i] = integer(card, 1, 'mid')
        #: Identification number of a TABLES1 or TABLEST entry. If H is
        #: given, then this field must be blank.
        self.table_id[i] = integer_or_blank(card, 2, 'tid')
        #: Type of material nonlinearity. ('NLELAST' for nonlinear elastic
        #: or 'PLASTIC' for elastoplastic.)
        self.Type[i] = string(card, 3, 'Type')

        if self.Type[i] == 'NLELAST':
            self.h[i] = blank(card, 4, 'h')
            self.hr[i] = blank(card, 6, 'hr')
            self.yf[i] = blank(card, 5, 'yf')
            self.limit1[i] = blank(card, 7, 'yf')
            self.limit2[i] = blank(card, 8, 'yf')
        else:
            #: Work hardening slope (slope of stress versus plastic strain) in
            #: units of stress. For elastic-perfectly plastic cases, H=0.0.
            #: For more than a single slope in the plastic range, the
            #: stress-strain data must be supplied on a TABLES1 entry
            #: referenced by TID, and this field must be blank
            h = double_or_blank(card, 4, 'H')
            self.h[i] = h
            if h is None:
                self.hflag[i] = False
            else:
                self.hflag[i] = True

            #: Yield function criterion, selected by one of the following
            #: values (1) Von Mises (2) Tresca (3) Mohr-Coulomb
            #: (4) Drucker-Prager
            self.yf[i] = integer_or_blank(card, 5, 'yf', 1)

            #: Hardening Rule, selected by one of the following values
            #: (Integer): (1) Isotropic (Default) (2) Kinematic
            #: (3) Combined isotropic and kinematic hardening
            self.hr[i] = integer_or_blank(card, 6, 'hr', 1)
            #: Initial yield point
            self.limit1[i] = double(card, 7, 'limit1')

            if self.yf[i] == 3 or self.yf[i] == 4:
                #: Internal friction angle, measured in degrees, for the
                #: Mohr-Coulomb and Drucker-Prager yield criteria
                self.limit2[i] = double(card, 8, 'limit2')
            else:
                #self.limit2[i] = blank(card, 8, 'limit2')
                #self.limit2[i] = None
                pass
        assert len(card) <= 9, 'len(MATS1 card) = %i' % len(card)
        self.i += 1

    def build(self):
        if self.n:
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

    def raw_fields(self):
        list_fields = ['MATS1', self.Mid(), self.Tid(), self.Type,
                  self.h, self.yf, self.hr, self.limit1, self.limit2]
        return list_fields

    def write_bdf(self, f, size=8, material_id=None):
        if size == 8:
            for mid, table_id, Type, h, hflag, yf, hr, limit1, limit2 in zip(
                self.material_id, self.table_id, self.Type,
                self.h, self.hflag, self.yf, self.hr, self.limit1, self.limit2):
                if not hflag:
                    h = None
                if limit2 is None:
                    card = ['MATS1', mid, table_id, Type, h, yf, hr, limit1]
                else:
                    card = ['MATS1', mid, table_id, Type, h, yf, hr, limit1, limit2]
                f.write(print_card_8(card))
        else:
            for mid, table_id, Type, h, hflag, yf, hr, limit1, limit2 in zip(
                self.material_id, self.table_id, self.Type,
                self.h, self.hflag, self.yf, self.hr, self.limit1, self.limit2):
                if not hflag:
                    h = None
                if limit2 is None:
                    card = ['MATS1', mid, table_id, Type, h, yf, hr, limit1]
                else:
                    card = ['MATS1', mid, table_id, Type, h, yf, hr, limit1, limit2]
                f.write(print_card_16(card))

    def reprFields(self):
        return self.raw_fields()
