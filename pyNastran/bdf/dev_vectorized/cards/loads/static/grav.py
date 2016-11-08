from six.moves import zip

import numpy as np
from numpy import zeros, unique

from pyNastran.bdf.field_writer_8 import print_card_8
from pyNastran.bdf.field_writer_16 import print_card_16
from pyNastran.bdf.bdf_interface.assign_type import (integer, integer_or_blank,
    double, double_or_blank)
from pyNastran.bdf.dev_vectorized.cards.vectorized_card import VectorizedCard


class GRAV(VectorizedCard):
    """
    +------+-----+-----+------+-----+-----+------+-----+
    | GRAV | SID | CID | A    | N1  | N2  | N3   |  MB |
    +------+-----+-----+------+-----+-----+------+-----+
    | GRAV | 1   | 3   | 32.2 | 0.0 | 0.0 | -1.0 |     |
    +------+-----+-----+------+-----+-----+------+-----+
    """
    type = 'GRAV'
    def __init__(self, model):
        """
        Defines the GRAV object.

        Parameters
        ----------
        model : BDF
           the BDF object

        .. todo:: collapse loads
        """
        VectorizedCard.__init__(self, model)
        #self.model = model
        #self.n = 0
        #self._cards = []
        #self._comments = []

    def __getitem__(self, i):
        #unique_lid = unique(self.load_id)
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

    def allocate(self, card_count):
        ncards = card_count[self.type]
        if ncards:
            self.n = ncards
            float_fmt = self.model.float_fmt
            #: Set identification number
            self.load_id = zeros(ncards, 'int32')
            #: Coordinate system identification number.
            self.coord_id = zeros(ncards, 'int32')
            #: scale factor
            self.scale = zeros(ncards, float_fmt)
            self.N = zeros((ncards, 3), float_fmt)
            self.mb = zeros(ncards, 'int32')

    def add_card(self, card, comment=''):
        #self._cards.append(card)
        #self._comments.append(comment)
        i = self.i
        self.load_id[i] = integer(card, 1, 'sid')
        #self.node_id[i] = integer(card, 1, 'node_id')
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
        assert len(card) <= 8, 'len(GRAV card) = %i\ncard=%s' % (len(card), card)
        self.i += 1

    def build(self):
        """
        Parameters
        ----------
        :param cards: the list of GRAV cards
        """
        if self.n:
            i = self.load_id.argsort()
            self.load_id = self.load_id[i]
            #self.node_id = self.node_id[i]
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

    def write_card(self, bdf_file, size=8, is_double=False, load_id=None):
        if self.n:
            for (lid, cid, scale, N, mb) in zip(
                 self.load_id, self.coord_id, self.scale, self.N, self.mb):

                card = ['GRAV', lid, cid, scale, N[0], N[1], N[2], mb]
                if size == 8:
                    bdf_file.write(print_card_8(card))
                else:
                    bdf_file.write(print_card_16(card))

    def get_load_ids(self):
        return np.unique(self.load_id)

