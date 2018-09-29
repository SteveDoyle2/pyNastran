from numpy import array, zeros, searchsorted

from pyNastran.bdf.field_writer_8 import print_card_8
from pyNastran.bdf.field_writer_16 import print_card_16
from pyNastran.bdf.bdf_interface.assign_type import (integer_or_blank,
    double, double_or_blank)
from pyNastran.dev.bdf_vectorized.cards.elements.property import Property


class PELAS(Property):
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

    def add_card(self, card, nPELAS=0, comment=''):
        self.n = 1
        #if comment:
            # self.comment = comment
        noffset = nPELAS * 5
        # 2 PELAS properties can be defined on 1 PELAS card
        # these are split into 2 separate cards

        pid = integer_or_blank(card, 1 + noffset, 'pid')
        if pid is not None:
            self._property_id.append(pid)
            self._K.append(double(card, 2 + noffset, 'k'))
            self._ge.append(double_or_blank(card, 3 + noffset, 'ge', 0.))
            self._s.append(double_or_blank(card, 4 + noffset, 's', 0.))

    def allocate(self, card_count):
        ncards = card_count[self.type]
        if ncards:
            float_fmt = self.model.float_fmt
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

    def update(self, maps):
        """
        maps = {
            'node' : nid_map,
            'property' : pid_map,
        }
        """
        if self.n:
            nid_map = maps['node']
            pid_map = maps['property']
            for i, pid in enumerate(self.property_id):
                self.property_id[i] = pid_map[pid]

    def write_card_by_index(self, bdf_file, size=8, is_double=False, i=None):
        if self.n:
            for (pid, k, ge, s) in zip(self.property_id[i], self.K[i], self.ge[i], self.s[i]):
                card = ['PELAS', pid, k, ge, s]
                if size == 8:
                    bdf_file.write(print_card_8(card))
                else:
                    bdf_file.write(print_card_16(card))

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
        i = self._validate_slice(i)
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
