from numpy import array, zeros, arange, concatenate, searchsorted, where, unique

from pyNastran.bdf.fieldWriter import print_card
from pyNastran.bdf.bdfInterface.assign_type import (integer, integer_or_blank,
    double_or_blank, integer_double_or_blank, blank)


class PBEAML(object):
    type = 'PBEAML'
    def __init__(self, model):
        """
        Defines the PCOMP object.

        :param self: the PBEAML object
        :param model: the BDF object
        :param cards: the list of PBEAML cards
        """
        self.model = model
        self.n = 0
        self._cards = []
        self._comments = []

    def add(self, card, comment):
        self._cards.append(card)
        self._comments.append(comment)

    def build(self):
        cards = self._cards
        ncards = len(cards)
        self.n = ncards

        if ncards:
            #: Property ID
            self.property_id = zeros(ncards, 'int32')
            self.material_id = zeros(ncards, 'int32')

            ncards = len(cards)
            for i, card in enumerate(cards):
                self.property_id[i] = integer(card, 1, 'property_id')
                self.material_id[i] = integer(card, 2, 'material_id')
            i = self.property_id.argsort()
            self.property_id = self.property_id[i]
            self.material_id = self.material_id[i]
            unique_pids = unique(self.property_id)

            if len(unique_pids) != len(self.property_id):
                raise RuntimeError('There are duplicate PCOMP IDs...')
            self._cards = []
            self._comments = []

    #=========================================================================
    def get_index(self, property_ids):
        if isinstance(property_ids, int):
            property_ids = array([property_ids])
        if property_ids is None:
            return arange(self.n)

        indexs = searchsorted(self.property_id, property_ids)
        assert len(indexs) == len(property_ids), 'indexs=%s pids=%s' % (indexs, property_ids)
        return indexs

    #=========================================================================
    def write_bdf(self, f, size=8, property_ids=None):
        if self.n:
            if property_ids is None:
                i = arange(self.n)
            else:
                i = searchsorted(self.property_id, property_ids)

            #cid = [cid if cid != 0 else '' for cid in self.coord_id]

            for (pid, mid) in zip(self.property_id[i], self.material_id[i]):
                card = ['PBEAML', pid, mid,]
                f.write(print_card(card))
