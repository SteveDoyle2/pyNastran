from numpy import array, zeros, arange, concatenate, searchsorted, where, unique, asarray

from pyNastran.bdf.dev_vectorized.cards.elements.property import Property

from pyNastran.bdf.fieldWriter import print_card
from pyNastran.bdf.bdfInterface.assign_type import (integer, integer_or_blank,
    double_or_blank, integer_double_or_blank, blank)
from pyNastran.bdf.dev_vectorized.utils import slice_to_iter


class PBARL(Property):
    type = 'PBARL'

    def __init__(self, model):
        """
        Defines the PCOMP object.

        :param self: the PBARL object
        :param model: the BDF object
        :param cards: the list of PBARL cards
        """
        Property.__init__(self, model)

    def allocate(self, ncards):
        float_fmt = self.model.float
        pass

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
        else:
            self.property_id = array([], dtype='int32')

    #=========================================================================
    def get_index(self, property_ids):
        if isinstance(property_ids, int):
            property_ids = array([property_ids])
        if property_ids is None:
            return arange(self.n)

        indexs = searchsorted(self.property_id, property_ids)
        assert len(indexs) == len(property_ids), 'indexs=%s pids=%s' % (indexs, property_ids)
        return indexs

    def __getitem__(self, property_ids):
        """
        Allows for slicing:
         - elements[1:10]
         - elements[4]
         - elements[1:10:2]
         - elements[[1,2,5]]
         - elements[array([1,2,5])]
        """
        i = searchsorted(self.property_id, property_ids)
        return self.slice_by_index(i)

    #def __getitem__(self, property_ids):
        #property_ids, int_flag = slice_to_iter(property_ids)
        #obj = PBARL(self.model)

        #properties = {}
        #for pid in sorted(property_ids):
            #properties[pid] = self.properties[pid]
        #obj.n = len(property_ids)
        #obj.properties = properties
        #obj.property_id = sorted(self.properties.keys())
        ##obj._comments = obj._comments[index]
        ##obj.comments = obj.comments[index]
        #return obj

    #=========================================================================
    def write_bdf(self, f, size=8, property_ids=None):
        if self.n:
            if property_ids is None:
                i = arange(self.n)
            else:
                i = searchsorted(self.property_id, property_ids)

            #cid = [cid if cid != 0 else '' for cid in self.coord_id]

            print('*pbarl write pids=%s' % self.property_id)
            for (pid, mid) in zip(self.property_id[i], self.material_id[i]):
                card = ['PBARL', pid, mid,]
                f.write(print_card(card))

    def slice_by_index(self, i):
        i = asarray(i)
        obj = PBARL(self.model)
        obj.n = len(i)
        #obj._cards = self._cards[i]
        #obj._comments = obj._comments[i]
        #obj.comments = obj.comments[i]
        obj.property_id = self.property_id[i]
        obj.material_id = self.material_id[i]
        #obj.area = self.area[i]
        #obj.I1 = self.I1[i]
        #obj.I2 = self.I2[i]
        #obj.J = self.J[i]
        #obj.nsm = self.nsm[i]
        return obj
