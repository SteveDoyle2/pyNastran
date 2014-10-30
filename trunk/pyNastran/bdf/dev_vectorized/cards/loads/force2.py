from six.moves import zip, StringIO

from numpy import zeros, searchsorted, unique

from pyNastran.bdf.fieldWriter import set_blank_if_default
from pyNastran.bdf.fieldWriter import print_card
from pyNastran.bdf.bdfInterface.assign_type import (integer, integer_or_blank,
    double, double_or_blank, string_or_blank)

from pyNastran.bdf.dev_vectorized.cards.vectorized_card import VectorizedCard

class FORCE2(VectorizedCard):
    type = 'FORCE2'
    def __init__(self, model):
        """
        Defines the FORCE2 object.

        :param self: the FORCE2 object
        :param model: the BDF object

        ..todo:: collapse loads
        """
        VectorizedCard.__init__(self, model)

    def allocate(self, ncards):
        float_fmt = self.model.float
        self.load_id = zeros(ncards, 'int32')
        self.node_id = zeros(ncards, 'int32')
        self.coord_id = zeros(ncards, 'int32')
        self.mag = zeros(ncards, float_fmt)
        self.xyz = zeros((ncards, 3), float_fmt)

    def __getitem__(self, i):
        unique_lid = unique(self.load_id)
        if len(i):
            f = FORCE2(self.model)
            f.load_id = self.load_id[i]
            f.node_id = self.node_id[i]
            f.coord_id = self.coord_id[i]
            f.mag = self.mag[i]
            f.xyz = self.xyz[i]
            f.n = len(i)
            return f
        raise RuntimeError('len(i) = 0')

    def __mul__(self, value):
        f = FORCE2(self.model)
        f.load_id = self.load_id
        f.node_id = self.node_id
        f.coord_id = self.coord_id
        f.mag = self.mag * value
        f.xyz = self.xyz * value
        f.n = self.n
        return f

    def __rmul__(self, value):
        return self.__mul__(value)

    def build(self):
        """
        :param self: the FORCE object
        :param cards: the list of FORCE cards
        """
        cards = self._cards
        ncards = len(cards)
        self.n = ncards
        if ncards:
            float_fmt = self.model.float
            #: Property ID
            self.load_id = zeros(ncards, 'int32')
            self.node_id = zeros(ncards, 'int32')
            self.coord_id = zeros(ncards, 'int32')
            self.mag = zeros(ncards, float_fmt)
            self.xyz = zeros((ncards, 3), float_fmt)

            for i, card in enumerate(cards):
                self.load_id[i] = integer(card, 1, 'sid')
                self.node_id[i] = integer(card, 2, 'node')
                self.coord_id[i] = integer_or_blank(card, 3, 'cid', 0)
                self.mag[i] = double(card, 4, 'mag')
                xyz = [double_or_blank(card, 5, 'X1', 0.0),
                       double_or_blank(card, 6, 'X2', 0.0),
                       double_or_blank(card, 7, 'X3', 0.0)]
                self.xyz[i] = xyz
                assert len(card) <= 8, 'len(FORCE card) = %i' % len(card)

            i = self.load_id.argsort()
            self.load_id = self.load_id[i]
            self.node_id = self.node_id[i]
            self.coord_id = self.coord_id[i]
            self.mag = self.mag[i]
            self.xyz = self.xyz[i]
            self._cards = []
            self._comments = []

    def write_bdf(self, f, size=8, lids=None):
        if self.n:
            for (lid, nid, cid, mag, xyz) in zip(
                 self.load_id, self.node_id, self.coord_id, self.mag, self.xyz):

                card = ['FORCE2', lid, nid, cid, mag, xyz[0], xyz[1], xyz[2] ]
                f.write(print_card(card))
