from six.moves import zip, StringIO
from numpy import array, zeros, unique, asarray, searchsorted, arange
from pyNastran.bdf.fieldWriter import print_card
from pyNastran.bdf.bdfInterface.assign_type import (integer,
    double, double_or_blank)

from pyNastran.bdf.dev_vectorized.cards.elements.property import Property

class PSHEAR(Property):
    type = 'PSHEAR'

    def __init__(self, model):
        Property.__init__(self, model)

    def allocate(self, ncards):
        self.n = ncards
        float_fmt = self.model.float
        self.property_id = zeros(ncards, 'int32')
        #: Material ID
        self.material_id = zeros(ncards, 'int32')
        self.thickness = zeros(ncards, float_fmt)
        self.nsm = zeros(ncards, float_fmt)
        self.f1 = zeros(ncards, float_fmt)
        self.f2 = zeros(ncards, float_fmt)

    def add(self, card, comment):
        i = self.i
        self.property_id[i] = integer(card, 1, 'pid')
        self.material_id[i] = integer(card, 2, 'mid')
        self.thickness[i] = double(card, 3, 't')
        self.nsm[i] = double_or_blank(card, 4, 'nsm', 0.0)
        self.f1[i] = double_or_blank(card, 5, 'f1', 0.0)
        self.f2[i] = double_or_blank(card, 6, 'f2', 0.0)
        assert self.thickness[i] > 0.0
        #assert self.f1 >= 0.0
        #assert self.f2 >= 0.0
        assert len(card) <= 7, 'len(PSHEAR card) = %i' % len(card)
        self.i += 1

    def build(self):
        if self.n:
            i = self.property_id.argsort()
            self.property_id = self.property_id[i]
            self.material_id = self.material_id[i]
            self.thickness = self.thickness[i]
            self.nsm = self.nsm[i]
            self.f1 = self.f1[i]
            self.f2 = self.f2[i]

            unique_pids = unique(self.property_id)
            if len(unique_pids) != len(self.property_id):
                raise RuntimeError('There are duplicate PSHEAR IDs...')
            self._cards = []
            self._comments = []
        else:
            self.property_id = array([], dtype='int32')

    def get_mass_per_area(self, property_id=None):
        # A * (rho * t + nsm)
        if property_id is None:
            property_id = self.property_id
            n = len(property_id)
            i = arange(n)
        else:
            print('property_id =', property_id)
            n = len(property_ids)
            i = searchsorted(self.property_id, property_id)
        mpa = zeros(n, dtype='float64')
        rho = self.model.materials.get_density_by_material_id(self.material_id)
        mpa = rho * self.thickness[i] + self.nsm[i]
        return mpa

    def write_bdf(self, f, size=8, pids=None):
        if self.n:
            if pids is None:
                i = arange(self.n)
            else:
                assert len(unique(pids))==len(pids), unique(pids)
                i = searchsorted(self.property_id, pids)

                #self.property_id = zeros(ncards, 'int32')
                #self.material_id = zeros(ncards, 'int32')
                #self.thickness = zeros(ncards, float_fmt)
                #self.nsm = zeros(ncards, float_fmt)
                #self.f1 = zeros(ncards, float_fmt)
                #self.f2 = zeros(ncards, float_fmt)

            for (pid, mid, t, nsm, f1, f2) in zip(self.property_id[i], self.material_id[i], self.thickness[i], self.nsm[i], self.f1[i], self.f2[i]):
                card = ['PSHEAR', pid, mid, t, nsm, f1, f2]
                f.write(print_card(card, size=size))

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

    def slice_by_index(self, i):
        i = asarray(i)
        obj = PSHEAR(self.model)
        obj.n = len(i)
        #obj._cards = self._cards[i]
        #obj._comments = obj._comments[i]
        #obj.comments = obj.comments[i]
        obj.property_id = self.property_id[i]
        obj.material_id = self.material_id[i]
        obj.thickness = self.thickness[i]
        obj.nsm = self.nsm[i]
        obj.f1 = self.f1[i]
        obj.f2 = self.f2[i]
        return obj