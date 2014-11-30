from pyNastran.bdf.dev_vectorized.cards.vectorized_card import VectorizedCard

class Element(VectorizedCard):
    def __init__(self, model):
        VectorizedCard.__init__(self, model)

    def __iter__(self):
        eids = self.element_id
        for eid in eids:
            yield eid

    def values(self):
        eids = self.element_id
        for eid in eids:
            yield self.__getitem__(eid)

    def items(self):
        eids = self.element_id
        for eid in eids:
            yield eid, self.__getitem__(eid)

    def get_element_id_by_element_index(self, i=None):
        #i = self._get_sorted_index(self.element_id, element_id, self.n,
        #                           'element_id', 'element_id in %s' % self.type, check=True)
        return self.element_id[i]

