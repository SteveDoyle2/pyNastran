from itertools import izip

from numpy import zeros, searchsorted, unique, where, array

from pyNastran.bdf.fieldWriter import print_card
from pyNastran.bdf.bdfInterface.assign_type import (integer, integer_or_blank,
    double, double_or_blank, string_or_blank)


class GRAV(object):
    type = 'GRAV'
    def __init__(self, model):
        """
        Defines the GRAV object.

        :param self: the GRAV object
        :param model: the BDF object
        
        ..todo:: collapse loads
        """
        self.model = model
        self.n = 0
        self._cards = []
        self._comments = []

    def __getitem__(self, i):
        unique_lid = unique(self.load_id)
        if len(i):
            f = GRAV(self.model)
            f.load_id = self.load_id[i]
            f.coord_id = self.coord_id[i]
            f.scale = self.scale[i]
            f.N = self.N[i]
            f.mb = self.mb[i]
            f.n = len(i)
            return f
        raise RuntimeError('len(i) = 0')

    def __mul__(self, value):
        f = GRAV(self.model)
        f.load_id = self.load_id
        f.coord_id = self.coord_id
        f.scale = self.scale * value
        f.N = self.N
        f.mb = self.mb
        f.n = self.n
        return f

    def __rmul__(self, value):
        return self.__mul__(value)

    def add(self, card, comment):
        self._cards.append(card)
        self._comments.append(comment)

    def build(self):
        """
        :param self: the GRAV object
        :param cards: the list of GRAV cards
        """
        cards = self._cards
        ncards = len(cards)
        self.n = ncards
        if ncards:
            float_fmt = self.model.float
            #: Set identification number
            self.load_id = zeros(ncards, 'int32')
            #: Coordinate system identification number.
            self.coord_id = zeros(ncards, 'int32')
            #: scale factor
            self.scale = zeros(ncards, float_fmt)
            self.N = zeros((ncards, 3), float_fmt)
            self.mb = zeros(ncards, 'int32')

            for i, card in enumerate(cards):
                self.load_id[i] = integer(card, 1, 'sid')
                self.coord_id[i] = integer_or_blank(card, 2, 'cid', 0)
                self.scale[i] = double(card, 3, 'scale')

                #: Acceleration vector components measured in coordinate system CID
                self.N[i, :] = [double_or_blank(card, 4, 'N1', 0.0),
                                double_or_blank(card, 5, 'N2', 0.0),
                                double_or_blank(card, 6, 'N3', 0.0)]

                #: Indicates whether the CID coordinate system is defined in the
                #: main Bulk Data Section (MB = -1) or the partitioned superelement
                #: Bulk Data Section (MB = 0). Coordinate systems referenced in the
                #: main Bulk Data Section are considered stationary with respect to
                #: the assembly basic coordinate system. See Remark 10.
                #: (Integer; Default = 0)
                self.mb[i] = integer_or_blank(card, 7, 'mb', 0)
                assert len(card) <= 8, 'len(GRAV card) = %i' % len(card)
            
            i = self.load_id.argsort()
            self.load_id = self.load_id[i]
            self.node_id = self.node_id[i]
            self.coord_id = self.coord_id[i]
            self.scale = self.scale[i]
            self.N = self.N[i]
            self._cards = []
            self._comments = []
        
    def get_stats(self):
        msg = []
        if self.n:
            msg.append('  %-8s: %i' % ('GRAV', self.n))
        return msg

    def write_bdf(self, f, size=8, lids=None):
        if self.n:
            for (lid, cid, scale, N, mb) in izip(
                 self.load_id, self.coord_id, self.scale, self.N, self.mb):

                card = ['GRAV', lid, cid, scale, N[0], N[1], n[2], mb]
                f.write(print_card(card))