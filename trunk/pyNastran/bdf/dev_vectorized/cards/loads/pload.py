from six.moves import zip
from numpy import arange, array, zeros, searchsorted, unique

from pyNastran.bdf.fieldWriter import set_blank_if_default
from pyNastran.bdf.fieldWriter import print_card
from pyNastran.bdf.bdfInterface.assign_type import (integer, integer_or_blank,
    double, double_or_blank, string, string_or_blank)


class PLOAD(object):
    type = 'PLOAD'

    def __init__(self, model):
        """
        Defines the PLOAD object.

        :param self: the PLOAD object
        :param model: the BDF object
        """
        self.model = model
        self.n = 0
        self._comments = []
        self.load_id = []
        self.node_ids = []
        self.pressure = []

    def __getitem__(self, i):
        unique_lid = unique(self.load_id)
        if len(i):
            f = PLOAD(self.model)
            f.load_id = self.load_id[i]
            f.element_id = self.element_id[i]
            f.pressure = self.pressure[i]
            f.node_ids = self.node_ids[i]
            f.n = len(i)
            return f
        raise RuntimeError('len(i) = 0')

    def __mul__(self, value):
        f = PLOAD(self.model)
        f.load_id = self.load_id
        f.element_id = self.element_id
        f.pressure = self.pressure
        f.node_ids = self.node_ids
        f.n = self.n
        return f

    def __rmul__(self, value):
        return self.__mul__(value)

    def add(self, card, comment=None):
        self.load_id.append(integer(card, 1, 'sid'))
        self.pressure.append(double(card, 2, 'p'))
        node_ids.append([integer(card, 3, 'n1'),
                      integer(card, 4, 'n2'),
                      integer(card, 5, 'n3'),
                      integer_or_blank(card, 6, 'n4', 0)])
        self.node_ids.append(node_ids)
        assert len(card) <= 7, 'len(PLOAD card) = %i' % len(card)
        self._comments.append(comment)

    def build(self):
        #if comment:
            #self._comment = comment
        #cards = self._cards
        #ncards = len(cards)

        float_fmt = self.model.float
        self.n = len(self.load_id)
        if self.n:
            self.load_id = zeros(self.load_id, 'int32')
            self.pressure = array(self.pressure, float_fmt)
            self.node_ids = array(self.node_ids, 'int32')

    def get_stats(self):
        msg = []
        if self.n:
            msg.append('  %-8s: %i' % ('PLOAD', self.n))
        return msg

    def write_bdf(self, f, size=8, load_ids=None):
        if self.n:
            if load_ids is None:
                i = arange(self.n)
            else:
                i = searchsorted(load_ids, self.load_id)

            n3 = ['' if n3i == 0 else n3i for n3i in self.node_ids[i, 3]]
            for (load_id, p, n, n3i) in zip(self.load_id[i], self.pressure[i],self.node_ids[i, :2], n3):
                card = ['PLOAD', load_id, p, n[0], n[1], n[2], n3i ]
                f.write(print_card(card))

