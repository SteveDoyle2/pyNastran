import numpy as np
from pyNastran.utils.numpy_utils import integer_types
from pyNastran.dev.bdf_vectorized.cards.vectorized_card import VectorizedCard

class VectorizedLoad(VectorizedCard):
    def __init__(self, model):
        VectorizedCard.__init__(self, model)

    def get_load_index_by_load_id(self, load_id):
        if load_id is None:
            return np.arange(self.n)
        #msg = ''
        assert isinstance(load_id, integer_types), 'load_id=%s type=%s' % (load_id, type(load_id))
        return np.where(self.load_id == load_id)[0]
        #i = self._get_sorted_index(self.load_id, load_id, 'load_id', 'load_id in %s%s' % (self.type, msg), check=True)

    def write_card(self, bdf_file, size=8, is_double=False, load_id=None):
        if self.n:
            i = self.get_load_index_by_load_id(load_id)
            self.write_card_by_index(bdf_file, size=size, is_double=is_double, i=i)
