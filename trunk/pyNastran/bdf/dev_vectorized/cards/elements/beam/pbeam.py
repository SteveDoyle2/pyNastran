from numpy import array, zeros, arange, concatenate, searchsorted, where, unique, asarray

from pyNastran.bdf.fieldWriter import print_card_8
from pyNastran.bdf.fieldWriter16 import print_card_16
from pyNastran.bdf.dev_vectorized.cards.elements.property import Property

from pyNastran.bdf.cards.properties.beam import PBEAM as vPBEAM
from pyNastran.bdf.dev_vectorized.utils import slice_to_iter


class PBEAM(object):
    type = 'PBEAM'

    def __len__(self):
        return self.n

    def __iter__(self):
        pids = self.property_id
        for pid in pids:
            yield pid

    def values(self):
        pids = self.property_id
        for pid in pids:
            yield self.__getitem__(pid)

    def items(self):
        pids = self.property_id
        for pid in pids:
            yield pid, self.__getitem__(pid)

    def __getitem__(self, property_ids):
        property_ids, int_flag = slice_to_iter(property_ids)
        obj = PBEAM(self.model)

        properties = {}
        for pid in sorted(property_ids):
            properties[pid] = self.properties[pid]
        obj.n = len(property_ids)
        obj.properties = properties
        obj.property_id = sorted(self.properties.keys())
        #obj._comments = obj._comments[index]
        #obj.comments = obj.comments[index]
        return obj

    def slice_by_index(self, i):
        i = asarray(i)
        obj = PBEAM(self.model)
        asdf
        return obj

    def __init__(self, model):
        """
        Defines the PBEAM object.

        :param self: the PBEAM object
        :param model: the BDF object
        :param cards: the list of PBEAM cards
        """
        self.properties = {}
        self.model = model
        self.n = 0

    def allocate(self, ncards):
        pass

    def add(self, card, comment):
        prop = vPBEAM(card, comment=comment)
        self.properties[prop.pid] = prop

    def build(self):
        self.n = len(self.properties)
        self.property_id = array(sorted(self.properties.keys()), dtype='int32')

    #=========================================================================
    def write_bdf(self, f, size=8, property_ids=None):
        if size == 8:
            for pid, prop in sorted(self.properties.iteritems()):
                f.write(prop.write_bdf(size, print_card_8))
        else:
            for pid, prop in sorted(self.properties.iteritems()):
                f.write(prop.write_bdf(size, print_card_16))
