from numpy import array, searchsorted
from pyNastran.dev.bdf_vectorized.cards.vectorized_card import VectorizedCard


class Material(VectorizedCard):
    def __init__(self, model):
        VectorizedCard.__init__(self, model)
        self.material_id = array([], dtype='int32')

    def __getitem__(self, i):
        return self.slice_by_index(i)

    def __iter__(self):
        for i in range(self.n):
            yield i

    def values(self):
        for i in range(self.n):
            yield self.__getitem__([i])

    def items(self):
        for i in range(self.n):
            yield i, self.__getitem__([i])

    def slice_by_material_id(self, material_id):
        """
        Allows for slicing:
         - materials[1:10]
         - materials[4]
         - materials[1:10:2]
         - materials[[1,2,5]]
         - materials[array([1,2,5])]
        """
        i = self.get_material_index_by_material_id(material_id)
        return self.slice_by_index(i)

    def get_material_index_by_material_id(self, material_id=None):
        if material_id is None:
            return slice(None)
        i = searchsorted(self.material_id, material_id)
        return i

