import numpy as np
from numpy import zeros, unique

#from pyNastran.bdf.field_writer_8 import set_blank_if_default
from pyNastran.bdf.field_writer_8 import print_card_8
from pyNastran.bdf.field_writer_16 import print_card_16
from pyNastran.bdf.bdf_interface.assign_type import (integer, integer_or_blank,
    double, double_or_blank)

from pyNastran.dev.bdf_vectorized.cards.loads.vectorized_load import VectorizedLoad

class MOMENT2(VectorizedLoad):
    """
    Defines a static concentrated moment at a grid point by specification
    of a magnitude and four grid points that determine the direction.::

    +---------+-----+---+---+----+----+----+----+
    | MOMENT2 | SID | G | M | G1 | G2 | G3 | G4 |
    +---------+-----+---+---+----+----+----+----+
    """
    type = 'MOMENT2'
    def __init__(self, model):
        """
        Defines the MOMENT2 object.

        Parameters
        ----------
        model : BDF
           the BDF object

        .. todo:: collapse loads
        """
        VectorizedLoad.__init__(self, model)

    def allocate(self, card_count):
        ncards = card_count[self.type]
        if ncards:
            self.n = ncards
            float_fmt = self.model.float_fmt
            self.load_id = zeros(ncards, 'int32')
            self.node_id = zeros(ncards, 'int32')
            self.coord_id = zeros(ncards, 'int32')
            self.mag = zeros(ncards, float_fmt)
            self.grids = zeros((ncards, 4), 'int32')
            self.xyz = zeros((ncards, 3), float_fmt)

    def cross_reference(self):
        g1 = self.grids[:, 0]
        g2 = self.grids[:, 1]
        g3 = self.grids[:, 2]
        g4 = self.grids[:, 3]
        nids = self.model.nids
        g1i = np.searchsorted(nids, g1)
        g2i = np.searchsorted(nids, g2)
        g3i = np.searchsorted(nids, g3)
        g4i = np.searchsorted(nids, g4)
        xyz1 = self.model.xyz_cid0[g1i, :]
        xyz2 = self.model.xyz_cid0[g2i, :]
        xyz3 = self.model.xyz_cid0[g3i, :]
        xyz4 = self.model.xyz_cid0[g4i, :]
        v21 = xyz2 - xyz1
        v43 = xyz4 - xyz3
        xyz = np.cross(v21, v43)
        self.xyz = xyz / np.linalg.norm(xyz, axis=0)
        #v12 = self.g2_ref.get_position() - self.g1_ref.get_position()
        #v34 = self.g4_ref.get_position() - self.g3_ref.get_position()
        #v12 = v12 / norm(v12)
        #v34 = v34 / norm(v34)
        #self.xyz = cross(v12, v34)

    def __getitem__(self, i):
        unique_lid = unique(self.load_id)
        if len(i):
            f = MOMENT2(self.model)
            f.load_id = self.load_id[i]
            f.node_id = self.node_id[i]
            f.mag = self.mag[i]
            self.grids = self.grids[i, :]
            f.xyz = self.xyz[i, :]
            f.n = len(i)
            return f
        raise RuntimeError('len(i) = 0')

    def __mul__(self, value):
        f = MOMENT2(self.model)
        f.load_id = self.load_id
        f.node_id = self.node_id
        self.grids = self.grids[i, :]
        f.mag = self.mag * value
        f.xyz = self.xyz * value
        f.n = self.n
        return f

    def __rmul__(self, value):
        return self.__mul__(value)

    def add_card(self, card, comment):
        i = self.i
        self.load_id[i] = integer(card, 1, 'sid')
        self.node_id[i] = integer(card, 2, 'node')
        self.mag[i] = double(card, 4, 'mag')
        n1234 = [
            integer(card, 4, 'g1'),
            integer(card, 5, 'g2'),
            integer(card, 6, 'g3'),
            integer(card, 7, 'g4'),
        ]
        self.grids = n1234
        #self.xyz[i] = xyz

        mag = double(card, 3, 'mag')
        assert len(card) <= 8, 'len(MOMENT2 card) = %i\ncard=%s' % (len(card), card)
        self.i += 1

    def build(self):
        if self.n:
            i = self.load_id.argsort()
            self.load_id = self.load_id[i]
            self.node_id = self.node_id[i]
            self.coord_id = self.coord_id[i]
            self.mag = self.mag[i]
            self.grids = self.grids[i, :]
            self.xyz = self.xyz[i, :]
            self._comments = []

    def write_card_by_index(self, bdf_file, size=8, is_double=False, i=None):
        for (lid, nid, cid, mag, n1234) in zip(
             self.load_id, self.node_id, self.mag, self.grids):

            card = ['MOMENT2', lid, nid, cid, mag] + list(n1234)
            if size == 8:
                bdf_file.write(print_card_8(card))
            else:
                bdf_file.write(print_card_16(card))
