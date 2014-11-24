from numpy import searchsorted
from pyNastran.bdf.dev_vectorized.cards.vectorized_card import VectorizedCard

class Property(VectorizedCard):
    def __init__(self, model):
        VectorizedCard.__init__(self, model)

    def shrink(self, refcheck=True):
        i = where(self.property_id==0)[0]
        self.resize(i[0], refcheck=refcheck)

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

    def __getitem__(self, property_id):
        """
        Allows for slicing:
         - elements[1:10]
         - elements[4]
         - elements[1:10:2]
         - elements[[1,2,5]]
         - elements[array([1,2,5])]
        """
        i = searchsorted(self.property_id, property_id)
        return self.slice_by_index(i)

    def get_property_id_by_property_index(self, i=None):
        if i is None:
            property_id = self.property_id
        else:
            property_id = self.property_id[i]
        return property_id

    def get_property_index_by_property_id(self, property_id=None, msg=''):
        i = self._get_sorted_index(self.property_id, property_id, self.n, 'property_id', 'property_id in %s%s' % (self.type, msg), check=True)
        return i
