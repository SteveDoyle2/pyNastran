import numpy as np
from numpy import where

from pyNastran.dev.bdf_vectorized.cards.vectorized_card import VectorizedCard

class Property(VectorizedCard):
    def __init__(self, model):
        VectorizedCard.__init__(self, model)

    def shrink(self, refcheck=True):
        i = where(self.property_id == 0)[0]
        self.resize(i[0], refcheck=refcheck)

    def __iter__(self):
        for i in range(self.n):
            yield i

    def values(self):
        for i in range(self.n):
            yield self.__getitem__([i])

    def items(self):
        for i in range(self.n):
            yield i, self.__getitem__([i])

    def __getitem__(self, i):
        return self.slice_by_index(i)

    def slice_by_property_id(self, property_id):
        """
        Allows for slicing:
         - properties[1:10]
         - properties[4]
         - properties[1:10:2]
         - properties[[1,2,5]]
         - properties[array([1,2,5])]
        """
        i = self.get_property_index_by_property_id(property_id)
        return self.slice_by_index(i)

    def get_property_id_by_property_index(self, i=None):
        if i is None:
            property_id = self.property_id
        else:
            property_id = self.property_id[i]
        return property_id

    def get_property_index_by_property_id(self, property_id=None, msg=''):
        i = self._get_sorted_index(self.property_id, property_id, 'property_id', 'property_id in %s%s' % (self.type, msg), check=True)
        return i

    def write_card(self, bdf_file, size=8, is_double=False, property_id=None):
        if self.n:
            if property_id is None:
                i = np.arange(self.n)
            else:
                i = np.searchsorted(self.property_id, property_id)
            return self.write_card_by_index(bdf_file, size=size, is_double=is_double, i=i)

    #def get_property_index_by_property_id(self, property_id=None):
        #if property_id is None:
            #return arange(self.n)
        #return searchsorted(self.property_id, property_id)
