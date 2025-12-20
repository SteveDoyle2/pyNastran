from numpy import zeros, unique

from pyNastran.bdf.field_writer_8 import set_blank_if_default
from pyNastran.bdf.field_writer_8 import print_card_8
from pyNastran.bdf.field_writer_16 import print_card_16
from pyNastran.bdf.bdf_interface.bdf_card import BDFCard
from pyNastran.bdf.bdf_interface.assign_type import (integer, integer_or_blank,
    double_or_blank)


class RFORCE:
    type = 'RFORCE'
    def __init__(self, model):
        """
        Defines the RFORCE object.

        Parameters
        ----------
        model : BDF
           the BDF object

        .. todo:: collapse loads
        """
        self.model = model
        self.n = 0
        self._cards = []
        self._comments = []

    def __getitem__(self, i):
        unique_lid = unique(self.load_id)
        if len(i):
            f = RFORCE(self.model)
            f.load_id = self.load_id[i]
            f.node_id = self.node_id[i]
            f.coord_id = self.coord_id[i]
            f.scale_vel = self.scale_vel[i]
            f.r = self.r[i]
            f.method = self.method[i]
            f.scale_acc = self.scale_acc[i]
            f.mb = self.mb[i]
            f.idrf = self.idrf[i]
            f.n = len(i)
            return f
        raise RuntimeError('len(i) = 0')

    def __mul__(self, value):
        f = RFORCE(self.model)
        f.load_id = self.load_id
        f.node_id = self.node_id
        f.coord_id = self.coord_id
        f.scale_vel = self.scale_vel * value
        f.r = self.r
        f.method = self.method
        f.scale_acc = self.scale_acc * value
        f.mb = self.mb
        f.idrf = self.idrf
        f.n = self.n
        return f

    def __rmul__(self, value):
        return self.__mul__(value)

    def allocate(self, ncards):
        float_fmt = self.model.float_fmt
        self.load_id = zeros(ncards, 'int32')
        self.node_id = zeros(ncards, 'int32')
        self.coord_id = zeros(ncards, 'int32')
        self.r = zeros((ncards, 3), float_fmt)
        self.scale_vel = zeros(ncards, float_fmt)
        self.method = zeros(ncards, 'int32')
        self.scale_acc = zeros(ncards, float_fmt)
        self.mb = zeros(ncards, 'int32')
        self.idrf = zeros(ncards, 'int32')

    def add_card(self, card: BDFCard, comment: str=''):
        self._cards.append(card)
        self._comments.append(comment)

    def build(self):
        """
        :param cards: the list of RFORCE cards
        """
        cards = self._cards
        ncards = len(cards)
        self.n = ncards
        if ncards:
            float_fmt = self.model.float_fmt
            #: Property ID
            self.load_id = zeros(ncards, 'int32')
            self.node_id = zeros(ncards, 'int32')
            self.coord_id = zeros(ncards, 'int32')
            self.r = zeros((ncards, 3), float_fmt)
            self.scale_vel = zeros(ncards, float_fmt)
            self.method = zeros(ncards, 'int32')
            self.scale_acc = zeros(ncards, float_fmt)
            self.mb = zeros(ncards, 'int32')
            self.idrf = zeros(ncards, 'int32')

            for i, card in enumerate(cards):
                self.load_id[i] = integer(card, 1, 'load_id')
                self.node_id[i] = integer(card, 2, 'node_id')
                self.coord_id[i] = integer_or_blank(card, 3, 'coord_id', 0)
                self.scale_vel[i] = double_or_blank(card, 4, 'scale', 1.)
                self.r[i] = [double_or_blank(card, 5, 'R1', 0.),
                             double_or_blank(card, 6, 'R2', 0.),
                             double_or_blank(card, 7, 'R3', 0.)]
                self.method[i] = integer_or_blank(card, 8, 'method', 1)
                self.scale_acc[i] = double_or_blank(card, 9, 'racc', 0.)
                self.mb[i] = integer_or_blank(card, 10, 'mb', 0)
                self.idrf[i] = integer_or_blank(card, 11, 'idrf', 0)
                assert len(card) <= 12, 'len(RFORCE card) = %i\ncard=%s' % (len(card), card)

            i = self.load_id.argsort()
            self.load_id = self.load_id[i]
            self.node_id = self.node_id[i]
            self.coord_id = self.coord_id[i]
            self.scale_vel = self.scale_vel[i]
            self.r = self.r[i]
            self.method = self.method[i]
            self.scale_acc = self.scale_acc[i]
            self.mb = self.mb[i]
            self.idrf = self.idrf[i]
            self._cards = []
            self._comments = []

    def get_stats(self):
        msg = []
        if self.n:
            msg.append('  %-8s: %i' % ('RFORCE', self.n))
        return msg

    def write_card(self, bdf_file, size=8, lids=None):
        if self.n:
            for (lid, nid, cid, scale_vel, r, method, scale_acc, mb, idrf) in zip(
                 self.load_id, self.node_id, self.coord_id, self.scale_vel,
                 self.r, self.scale_acc, self.mb, self.idrf):

                #method = set_blank_if_default(method, 1)
                scale_acc = set_blank_if_default(scale_acc, 0.)
                mb = set_blank_if_default(mb, 0)
                idrf = set_blank_if_default(idrf, 0)
                card = ['RFORCE', lid, nid, cid, scale_vel,
                        r[0], r[1], r[2], method, scale_acc, mb, idrf]
                if size == 8:
                    bdf_file.write(print_card_8(card))
                else:
                    bdf_file.write(print_card_16(card))
