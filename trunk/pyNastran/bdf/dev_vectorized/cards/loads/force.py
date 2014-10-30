import StringIO
from six.moves import zip

from numpy import zeros, searchsorted, unique, where

from pyNastran.bdf.fieldWriter import set_blank_if_default
from pyNastran.bdf.fieldWriter import print_card
from pyNastran.bdf.bdfInterface.assign_type import (integer, integer_or_blank,
    double, double_or_blank, string_or_blank)


class FORCE(object):
    type = 'FORCE'
    def __init__(self, model):
        """
        Defines the FORCE object.

        :param self: the FORCE object
        :param model: the BDF object

        ..todo:: collapse loads
        """
        self.model = model
        self.n = 0
        self._cards = []
        self._comments = []

    def __contains__(self, key):
        """TODO: should check against unique values"""
        if key in self.load_id:
            return True
        return False
        #return dict.__contains__(self, self.__keytransform__(key))

    def get_load_ids(self):
        #print('load_id = %s' % self.load_id)
        return unique(self.load_id)

    def __getitem__(self, i):
        unique_lid = unique(self.load_id)
        if len(i):
            f = FORCE(self.model)
            f.load_id = self.load_id[i]
            f.node_id = self.node_id[i]
            f.coord_id = self.coord_id[i]
            f.mag = self.mag[i]
            f.xyz = self.xyz[i]
            f.n = len(i)
            return f
        raise RuntimeError('len(i) = 0')

    def __mul__(self, value):
        f = FORCE(self.model)
        f.load_id = self.load_id
        f.node_id = self.node_id
        f.coord_id = self.coord_id
        f.mag = self.mag * value
        f.xyz = self.xyz * value
        f.n = self.n
        return f

    def __rmul__(self, value):
        return self.__mul__(value)

    def add(self, card, comment):
        self._cards.append(card)
        self._comments.append(comment)

    def allocate(self, ncards):
        float_fmt = self.model.float
        self.load_id = zeros(ncards, 'int32')
        self.node_id = zeros(ncards, 'int32')
        self.coord_id = zeros(ncards, 'int32')
        self.mag = zeros(ncards, float_fmt)
        self.xyz = zeros((ncards, 3), float_fmt)

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

    def get_stats(self):
        msg = []
        if self.n:
            msg.append('  %-8s: %i' % ('FORCE', self.n))
        return msg

    def write_bdf(self, f, size=8, is_double=False, load_id=None):
        if self.n:
            if load_id is None:
                for (lid, nid, cid, mag, xyz) in zip(
                     self.load_id, self.node_id, self.coord_id, self.mag, self.xyz):

                    card = ['FORCE', lid, nid, cid, mag, xyz[0], xyz[1], xyz[2] ]
                    f.write(print_card(card))
            else:
                for lid in unique(load_id):
                    i = where(self.load_id == lid)[0]
                    for (lid, nid, cid, mag, xyz) in zip(
                         self.load_id[i], self.node_id[i], self.coord_id[i], self.mag[i], self.xyz[i]):

                        card = ['FORCE', lid, nid, cid, mag, xyz[0], xyz[1], xyz[2] ]
                        f.write(print_card(card))

    def __repr__(self):
        f = StringIO.StringIO()
        self.write_bdf(f)
        return f.getvalue().rstrip()