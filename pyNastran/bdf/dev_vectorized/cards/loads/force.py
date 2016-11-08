from six.moves import zip, StringIO

from numpy import zeros, unique, where

from pyNastran.bdf.field_writer_8 import print_card_8
from pyNastran.bdf.field_writer_16 import print_card_16
from pyNastran.bdf.bdf_interface.assign_type import (
    integer, integer_or_blank, double, double_or_blank)
from pyNastran.bdf.dev_vectorized.cards.vectorized_card import VectorizedCard


class FORCE(VectorizedCard):
    type = 'FORCE'
    def __init__(self, model):
        """
        Defines the FORCE object.

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
            obj = FORCE(self.model)
            obj.load_id = self.load_id[i]
            obj.node_id = self.node_id[i]
            obj.coord_id = self.coord_id[i]
            obj.mag = self.mag[i]
            obj.xyz = self.xyz[i]
            obj.n = len(i)
            return obj
        raise RuntimeError('len(i) = 0')

    def __mul__(self, value):
        obj = FORCE(self.model)
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
        #self._comments.append(comment)
        i = self.i
        self.load_id[i] = integer(card, 1, 'sid')
        self.node_id[i] = integer(card, 2, 'node')
        self.coord_id[i] = integer_or_blank(card, 3, 'cid', 0)
        self.mag[i] = double(card, 4, 'mag')
        xyz = [double_or_blank(card, 5, 'X1', 0.0),
               double_or_blank(card, 6, 'X2', 0.0),
               double_or_blank(card, 7, 'X3', 0.0)]
        self.xyz[i, :] = xyz
        assert len(card) <= 8, 'len(FORCE card) = %i\ncard=%s' % (len(card), card)
        self.i += 1

    def allocate(self, card_count):
        ncards = card_count[self.type]
        if ncards:
            self.n = ncards
            float_fmt = self.model.float_fmt
            self.load_id = zeros(ncards, 'int32')
            self.node_id = zeros(ncards, 'int32')
            self.coord_id = zeros(ncards, 'int32')
            self.mag = zeros(ncards, float_fmt)
            self.xyz = zeros((ncards, 3), float_fmt)

    def build(self):
        """
        Parameters
        ----------
        cards : the list of FORCE cards
           cards
        """
        pass
        # if we argsort this, we screw up the order
        #if self.n:
            #i = self.load_id.argsort()
            #print('self.load_id =', self.load_id)
            #print('i =', i)
            #self.load_id = self.load_id[i]
            #self.node_id = self.node_id[i]
            #self.coord_id = self.coord_id[i]
            #self.mag = self.mag[i]
            #self.xyz = self.xyz[i, :]

    def get_stats(self):
        msg = []
        if self.n:
            msg.append('  %-8s: %i' % ('FORCE', self.n))
        return msg

    def write_card(self, bdf_file, size=8, is_double=False, load_id=None):
        if self.n:
            if load_id is None:
                for (lid, nid, cid, mag, xyz) in zip(
                        self.load_id, self.node_id, self.coord_id, self.mag, self.xyz):

                    card = ['FORCE', lid, nid, cid, mag, xyz[0], xyz[1], xyz[2]]
                    if size == 8:
                        bdf_file.write(print_card_8(card))
                    else:
                        bdf_file.write(print_card_16(card))
            else:
                for lid in unique(load_id):
                    i = where(self.load_id == lid)[0]
                    for (lid, nid, cid, mag, xyz) in zip(
                            self.load_id[i], self.node_id[i], self.coord_id[i],
                            self.mag[i], self.xyz[i]):

                        card = ['FORCE', lid, nid, cid, mag, xyz[0], xyz[1], xyz[2]]
                        if size == 8:
                            bdf_file.write(print_card_8(card))
                        else:
                            bdf_file.write(print_card_16(card))

    def __repr__(self):
        file_obj = StringIO()
        self.write_card(file_obj)
        return file_obj.getvalue().rstrip()
