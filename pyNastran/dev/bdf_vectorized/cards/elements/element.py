from pyNastran.utils.numpy_utils import integer_types
from pyNastran.dev.bdf_vectorized.cards.vectorized_card import VectorizedCard

class Element(VectorizedCard):
    """
    an Element is a subclass for the Elements that defines:
      - for
    """
    def __init__(self, model):
        VectorizedCard.__init__(self, model)

    def __iter__(self):
        """
        for ielem in range(model.elements):
           ...
        """
        for i in self.n:
            yield i

    def values(self):
        """
        for elem in model.elements.values():
           ...
        """
        for i in self.n:
            yield self.__getitem__([i])

    def items(self):
        """
        for ielem, elem in model.elements.items():
           ...
        """
        for i in self.n:
            yield i, self.__getitem__([i])

    def __getitem__(self, i):
        return self.slice_by_index(i)

    def slice_by_element_id(self, element_id):
        """
        Allows for slicing:
         - elements[1:10]
         - elements[4]
         - elements[1:10:2]
         - elements[[1,2,5]]
         - elements[array([1,2,5])]
        """
        i = self.get_element_index_by_element_id(element_id)
        return self.slice_by_index(i)

    def get_element_id_by_element_index(self, i=None):
        #i = self._get_sorted_index(self.element_id, element_id, self.n,
        #                           'element_id', 'element_id in %s' % self.type, check=True)
        return self.element_id[i]

    def get_element_index_by_element_id(self, element_id=None, msg=''):
        self.model.log.debug('Type=%s' % self.type)
        self.model.log.debug('element_id = %s' % element_id)
        if isinstance(element_id, integer_types):
            assert element_id > 0, element_id
        elif element_id is None:
            pass
        elif isinstance(element_id, list):
            assert min(element_id) > 0, element_id
        else:
            assert element_id.min() > 0, element_id
        i = self._get_sorted_index(self.element_id, element_id, 'element_id', 'element_id in %s%s' % (self.type, msg), check=True)
        return i
