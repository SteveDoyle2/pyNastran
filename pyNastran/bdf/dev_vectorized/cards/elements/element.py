from pyNastran.bdf.dev_vectorized.cards.vectorized_card import VectorizedCard

class Element(VectorizedCard):
    def __init__(self, model):
        VectorizedCard.__init__(self, model)

    def __iter__(self):
        for i in self.n:
            yield i

    def values(self):
        for i in self.n:
            yield self.__getitem__([i])

    def items(self):
        for i in self.n:
            yield i, self.__getitem__([i])

    def __getitem__(self, i):
        return self.slice_by_index(i)

    def get_element_id_by_element_index(self, i=None):
        #i = self._get_sorted_index(self.element_id, element_id, self.n,
        #                           'element_id', 'element_id in %s' % self.type, check=True)
        return self.element_id[i]

