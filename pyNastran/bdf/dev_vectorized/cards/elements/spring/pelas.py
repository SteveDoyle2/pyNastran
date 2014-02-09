from itertools import izip

from numpy import array, arange, dot, zeros, unique, searchsorted
from numpy.linalg import norm

from pyNastran.bdf.fieldWriter import print_card
from pyNastran.bdf.bdfInterface.assign_type import (integer, integer_or_blank,
    double, double_or_blank, integer_double_or_blank, blank)

class PELAS(object):
    """
    Specifies the stiffness, damping coefficient, and stress coefficient of a
    scalar elastic (spring) element (CELAS1 or CELAS3 entry).
    """
    type = 'PELAS'

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

    def write_bdf(self, f, size=8, pids=None):
        if self.n:
            if pids is None:
                i = arange(self.n)
            else:
                i = pids

            for (pid, k, ge, s) in izip(self.property_id[i], self.K[i], self.ge[i], self.s[i]):
                card = ['PELAS', pid, k, ge, s]
                f.write(print_card(card, size=size))