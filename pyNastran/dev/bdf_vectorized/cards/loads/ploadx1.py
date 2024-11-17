from numpy import arange, zeros, unique

from pyNastran.bdf.field_writer_8 import print_card_8
from pyNastran.bdf.field_writer_16 import print_card_16
from pyNastran.bdf.bdf_interface.assign_type import (integer,
    double, double_or_blank)


class PLOADX1:
    type = 'PLOADX1'
    def __init__(self, model):
        """
        Defines the PLOADX1 object.

        Parameters
        ----------
        model : BDF
           the BDF object
        """
        self.model = model
        self.n = 0
        self._cards = []
        self._comments = []

    def __getitem__(self, i):
        unique_lid = unique(self.load_id)
        if len(i):
            f = PLOADX1(self.model)
            f.load_id = self.load_id[i]
            f.element_id = self.element_id[i, :]
            f.p = self.p[i]
            f.node_ids = self.node_ids[i, :]
            f.theta = self.theta[i]
            f.n = len(i)
            return f
        raise RuntimeError('len(i) = 0')

    def __mul__(self, value):
        f = PLOADX1(self.model)
        f.load_id = self.load_id
        f.element_id = self.element_id
        f.p = self.p * value
        f.node_ids = self.node_ids
        f.theta = self.theta
        f.n = self.n
        return f

    def __rmul__(self, value):
        return self.__mul__(value)

    def add_card(self, card: BDFCard, comment: str=''):
        self._cards.append(card)
        self._comments.append(comment)

    def build(self):
        cards = self._cards
        ncards = len(cards)
        self.n = ncards
        if ncards:
            float_fmt = self.model.float_fmt
            #: Property ID
            self.load_id = zeros(ncards, 'int32')
            #: Element ID
            self.element_id = zeros(ncards, 'int32')

            # Surface traction at grid point GA. (Real)
            self.p = zeros((ncards, 2), float_fmt)

            #: Corner grid points. GA and GB are any two adjacent corner grid points of the
            #: element. (Integer > 0)
            self.node_ids = zeros((ncards, 2), 'int32')

            #: Angle between surface traction and inward normal to the line segment.
            #: (Real Default = 0.0)
            self.theta = zeros(ncards, float_fmt)

            for i, card in enumerate(cards):
                self.load_id[i] = integer(card, 1, 'load_id')
                self.element_id[i] = integer(card, 2, 'element_id')
                pa = double(card, 3, 'pa')
                pb = double_or_blank(card, 4, 'pb', pa)
                self.p[i, :] = [pa, pb]
                self.node_ids[i, :] = [integer(card, 5, 'ga'),
                                       integer(card, 6, 'gb')]
                self.theta[i] = double_or_blank(card, 7, 'theta', 0.)
                assert len(card) <= 8, 'len(PLOADX1 card) = %i\ncard=%s' % (len(card), card)

            i = self.load_id.argsort()
            self.load_id = self.load_id[i]
            self.element_id = self.element_id[i]
            self.node_ids = self.node_ids[i, :]
            self.p = self.p[i, :]
            self.theta = self.theta[i]
            self._cards = []
            self._comments = []

    def get_stats(self):
        msg = []
        if self.n:
            msg.append('  %-8s: %i' % ('PLOADX1', self.n))
        return msg

    def get_index(self, load_ids=None):
        #if load_ids:
        i = arange(self.n)
        #else:
        #    i = searchsorted(load_ids, self.load_id)
        return i

    def write_card(self, bdf_file, size=8, is_double=False, load_id=None):
        if self.n:
            i = self.get_index(load_ids)
            for (lid, eid, p, n, theta) in zip(self.load_id[i],
                    self.element_id[i], self.p[i], self.node_ids[i], self.theta[i]):
                card = ['PLOADX1', lid, eid, p[0], p[1], n[0], n[1], theta]
                if size == 8:
                    bdf_file.write(print_card_8(card))
                else:
                    bdf_file.write(print_card_16(card))
