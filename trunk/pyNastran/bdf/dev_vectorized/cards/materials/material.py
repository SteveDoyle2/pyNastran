from numpy import array, where, searchsorted
from pyNastran.bdf.dev_vectorized.cards.vectorized_card import VectorizedCard


class Material(VectorizedCard):
    def __init__(self, model):
        VectorizedCard.__init__(self, model)
        self.material_id = array([], dtype='int32')

    def __getitem__(self, material_id):
        #self.model.log.debug('self.material_id = %s' % self.material_id)
        #self.model.log.debug('material_id = %s' % material_id)
        #material_id = slice_to_iter(material_id)
        i = where(self.material_id == material_id)[0]
        return self.slice_by_index(i)

    def __iter__(self):
        mids = self.material_id
        for mid in mids:
            yield mid

    def values(self):
        mids = self.material_id
        for mid in mids:
            yield self.__getitem__(mid)

    def items(self):
        mids = self.material_id
        #self.model.log.debug('mids = %s' % mids)
        for mid in mids:
            yield mid, self.__getitem__(mid)

    def iteritems(self):
        return self.items()

    def iterkeys(self):
        return self.keys()

    def __len__(self):
        return self.n

    def get_material_index_by_material_id(self, material_id=None):
        if material_id is None:
            return slice(None)
        i = searchsorted(self.material_id, material_id)
        return i


