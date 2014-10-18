from numpy import searchsorted
from pyNastran.bdf.dev_vectorized.cards.vectorized_card import VectorizedCard

class Property(VectorizedCard):
    def __init__(self, model):
        VectorizedCard.__init__(self, model)

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
