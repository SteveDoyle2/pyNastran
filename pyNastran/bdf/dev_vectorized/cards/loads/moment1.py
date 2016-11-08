from __future__ import print_function
from six.moves import zip, StringIO

import numpy as np
from numpy import zeros, unique

#from pyNastran.bdf.field_writer_8 import set_blank_if_default
from pyNastran.bdf.field_writer_8 import print_card_8
from pyNastran.bdf.field_writer_16 import print_card_16
from pyNastran.bdf.bdf_interface.assign_type import (
    integer, integer_or_blank, double, double_or_blank)

from pyNastran.bdf.dev_vectorized.cards.vectorized_card import VectorizedCard

class MOMENT1(VectorizedCard):
    """
    Defines a static concentrated moment at a grid point by specifying a
    magnitude and two grid points that determine the direction.::

    +---------+-----+---+---+----+----+
    |    1    |  2  | 3 | 4 | 5  | 6  |
    +=========+=====+===+===+====+====+
    | MOMENT1 | SID | G | M | G1 | G2 |
    +---------+-----+---+---+----+----+
    """
    type = 'MOMENT1'
    def __init__(self, model):
        """
        Defines the MOMENT1 object.

        Parameters
        ----------
        model : BDF
           the BDF object
        """
        VectorizedCard.__init__(self, model)

    def get_load_ids(self):
        return unique(self.load_id)

    def __getitem__(self, i):
        unique_lid = unique(self.load_id)
        if len(i):
            obj = MOMENT1(self.model)
            obj.load_id = self.load_id[i]
            obj.node_id = self.node_id[i]
            obj.coord_id = self.coord_id[i]
            obj.mag = self.mag[i]
            obj.xyz = self.xyz[i]
            obj.n = len(i)
            return obj
        raise RuntimeError('len(i) = 0')

    def __mul__(self, value):
        obj = MOMENT1(self.model)
        obj.load_id = self.load_id
        obj.node_id = self.node_id
        obj.coord_id = self.coord_id
        obj.mag = self.mag * value
        obj.xyz = self.xyz * value
        obj.n = self.n
        return obj

    def __rmul__(self, value):
        return self.__mul__(value)

    def add_card(self, card, comment=''):
        i = self.i
        self.load_id[i] = integer(card, 1, 'sid')
        self.node_id[i] = integer(card, 2, 'node')
        self.mag[i] = double(card, 3, 'mag')
        n12 = [
            integer(card, 4, 'g1'),
            integer(card, 5, 'g2'),
        ]
        self.grids[i, :] = n12
        assert len(card) == 6, 'len(MOMENT1 card) = %i\ncard=%s' % (len(card), card)
        self.i += 1

    def allocate(self, card_count):
        ncards = card_count[self.type]
        if ncards:
            self.n = ncards
        float_fmt = self.model.float_fmt
        self.load_id = zeros(ncards, 'int32')
        self.node_id = zeros(ncards, 'int32')
        self.mag = zeros(ncards, float_fmt)
        self.grids = zeros((ncards, 2), 'int32')
        self.xyz = zeros((ncards, 3), float_fmt)

    def build(self):
        if self.n:
            i = self.load_id.argsort()
            self.load_id = self.load_id[i]
            self.node_id = self.node_id[i]
            self.mag = self.mag[i]
            self.grids = self.grids[i, :]
            self.xyz = self.xyz[i]
            #self._comments = []

    def cross_reference(self):
        g1 = self.grids[:, 0]
        g2 = self.grids[:, 1]
        nids = self.model.nids
        g1i = np.searchsorted(nids, g1)
        g2i = np.searchsorted(nids, g2)
        xyz1 = self.model.xyz_cid0[g1i, :]
        xyz2 = self.model.xyz_cid0[g2i, :]

        xyz = xyz2 - xyz1
        self.xyz = xyz / np.linalg.norm(xyz, axis=0)

    def get_stats(self):
        msg = []
        if self.n:
            msg.append('  %-8s: %i' % ('MOMENT', self.n))
        return msg

    def write_card(self, bdf_file, size=8, is_double=False, load_id=None):
        if self.n:
            for (lid, nid, cid, mag, xyz) in zip(
                 self.load_id, self.node_id, self.coord_id, self.mag, self.xyz):

                card = ['MOMENT1', lid, nid, cid, mag, xyz[0], xyz[1], xyz[2]]
                if size == 8:
                    bdf_file.write(print_card_8(card))
                else:
                    bdf_file.write(print_card_16(card))

    def __repr__(self):
        string_io = StringIO()
        self.write_card(string_io)
        return string_io.getvalue().rstrip()
