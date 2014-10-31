from six import iteritems
from six.moves import StringIO
from numpy import array, zeros, arange, searchsorted, where, unique, hstack, concatenate

from pyNastran.bdf.dev_vectorized.utils import slice_to_iter

from pyNastran.bdf.fieldWriter import print_card_8
from pyNastran.bdf.fieldWriter16 import print_card_16
from pyNastran.bdf.bdfInterface.assign_type import (integer, integer_or_blank,
    double_or_blank, integer_double_or_blank, blank)
from pyNastran.bdf.cards.properties.bush import PBUSH as vPBUSH

class PBUSH(object):
    type = 'PBUSH'
    def __len__(self):
        return self.n

    def __init__(self, model):
        """
        Defines the PCOMP object.

        :param self: the PBUSH object
        :param model: the BDF object
        :param cards: the list of PBUSH cards
        """
        self.model = model
        self.properties = {}
        #self._cards = []
        #self._comments = []

    def add(self, card, comment):
        prop = vPBUSH(card, comment=comment)
        self.properties[prop.pid] = prop
        #self._cards.append(card)
        #self._comments.append(comment)

    def allocate(self, ncards):
        pass

    def get_properties(self, property_id=None):
        props = []
        if property_id is None:
            property_id = self.property_id
        for pid in self.properties:
            prop = self.properties[pid]
            props.append(prop)
        return props

    def build(self):
        self.n = len(self.properties)
        self.property_id = array(sorted(self.properties.keys()), dtype='int32')
        return
        self.n = 0
        cards = self._cards
        ncards = len(cards)
        self.n = ncards
        #return
        if ncards:
            #: Property ID
            self.property_id = zeros(ncards, 'int32')

            ncards = len(cards)
            for i, card in enumerate(cards):
                self.property_id[i] = integer(card, 1, 'pid')
            i = self.property_id.argsort()
            self.property_id = self.property_id[i]
            unique_pids = unique(self.property_id)

            if len(unique_pids) != len(self.property_id):
                raise RuntimeError('There are duplicate PCOMP IDs...')
            self._cards = []
            self._comments = []

    def write_bdf(self, f, size=8, property_ids=None):
        if size == 8:
            for pid, pcomp in sorted(iteritems(self.properties)):
                f.write(pcomp.write_bdf(size, print_card_8))
        else:
            for pid, pcomp in sorted(iteritems(self.properties)):
                f.write(pcomp.write_bdf(size, print_card_16))

    #def slice_by_index(self, i):
        #i = asarray(i)
        #return self.__getitem__()

    def __getitem__(self, property_ids):
        property_ids, int_flag = slice_to_iter(property_ids)
        obj = PBUSH(self.model)

        properties = {}
        for pid in sorted(property_ids):
            properties[pid] = self.properties[pid]
        obj.n = len(property_ids)
        obj.properties = properties
        obj.property_id = sorted(self.properties.keys())
        #obj._comments = obj._comments[index]
        #obj.comments = obj.comments[index]
        return obj


    def __repr__(self):
        f = StringIO.StringIO()
        f.write('<PBUSH object> n=%s\n' % self.n)
        self.write_bdf(f)
        return f.getvalue()
