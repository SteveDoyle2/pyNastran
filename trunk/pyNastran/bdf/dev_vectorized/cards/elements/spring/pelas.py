from six.moves import zip

from numpy import array, arange, dot, zeros, unique, searchsorted, asarray
from numpy.linalg import norm

from pyNastran.bdf.fieldWriter import print_card_8
from pyNastran.bdf.fieldWriter16 import print_card_16
from pyNastran.bdf.bdfInterface.assign_type import (integer, integer_or_blank,
    double, double_or_blank, integer_double_or_blank, blank)


class PELAS(object):
    """
    Specifies the stiffness, damping coefficient, and stress coefficient of a
    scalar elastic (spring) element (CELAS1 or CELAS3 entry).
    """
    type = 'PELAS'

    def __len__(self):
        return self.n

    def __init__(self, model):
        self.model = model
        #: Card count
        self.n = 0
        self._property_id = []
        self._K = []
        self._ge = []
        self._s = []

    def add(self, card, nPELAS=0, comment=''):
        self.n = 1
        #if comment:
            #self._comment = comment
        nOffset = nPELAS * 5
        # 2 PELAS properties can be defined on 1 PELAS card
        # these are split into 2 separate cards

        pid = integer_or_blank(card, 1 + nOffset, 'pid')
        if pid is not None:
            self._property_id.append(pid)
            self._K.append(double(card, 2 + nOffset, 'k'))
            self._ge.append(double_or_blank(card, 3 + nOffset, 'ge', 0.))
            self._s.append(double_or_blank(card, 4 + nOffset, 's', 0.))

    def allocate(self, ncards):
        float_fmt = self.model.float
        self.property_id = zeros(ncards, dtype='int32')
        self.K = zeros(ncards, dtype=float_fmt)
        self.ge = zeros(ncards, dtype=float_fmt)
        self.s = zeros(ncards, dtype=float_fmt)

    def build(self):
        if self.n:
            #: Property identification number. (Integer > 0)
            self.property_id = array(self._property_id)
            #: Ki Elastic property value. (Real)
            self.K = array(self._K)
            #: Damping coefficient, . See Remarks 5. and 6. (Real)
            #: To obtain the damping coefficient GE, multiply the
            #: critical damping ratio c/c0 by 2.0.
            self.ge = array(self._ge)
            #: Stress coefficient. (Real)
            self.s = array(self._s)
            self.n = len(self.K)

            i = self.property_id.argsort()
            self.property_id = self.property_id[i]
            self.K = self.K[i]
            self.ge = self.ge[i]
            self.s = self.s[i]

            self._property_id = []
            self._K = []
            self._ge = []
            self._s = []
        else:
            self.property_id = array([], dtype='int32')

    def write_bdf(self, f, size=8, property_id=None):
        if self.n:
            if property_id is None:
                i = arange(self.n)
            else:
                i = property_id

            for (pid, k, ge, s) in zip(self.property_id[i], self.K[i], self.ge[i], self.s[i]):
                card = ['PELAS', pid, k, ge, s]
                if size == 8:
                    f.write(print_card_8(card))
                else:
                    f.write(print_card_16(card))

    def __getitem__(self, property_ids):
        """
        Allows for slicing:
         - elements[1:10]
         - elements[4]
         - elements[1:10:2]
         - elements[[1,2,5]]
         - elements[array([1,2,5])]
        """
        i = searchsorted(self.property_id, property_ids)
        return self.slice_by_index(i)

    def slice_by_index(self, i):
        i = asarray(i)
        obj = PELAS(self.model)
        obj.n = len(i)
        #obj._cards = self._cards[i]
        #obj._comments = obj._comments[i]
        #obj.comments = obj.comments[i]
        obj.property_id = self.property_id[i]
        obj.K = self.K[i]
        obj.ge = self.ge[i]
        obj.s = self.s[i]
        return obj
