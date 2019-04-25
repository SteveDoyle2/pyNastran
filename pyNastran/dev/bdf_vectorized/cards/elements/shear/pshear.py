from numpy import array, zeros, unique, searchsorted, arange
from pyNastran.bdf.field_writer_8 import print_card_8
from pyNastran.bdf.field_writer_16 import print_card_16
from pyNastran.bdf.bdf_interface.assign_type import (integer,
    double, double_or_blank)

from pyNastran.dev.bdf_vectorized.cards.elements.property import Property

class PSHEAR(Property):
    type = 'PSHEAR'

    def __init__(self, model):
        Property.__init__(self, model)

    def allocate(self, card_count):
        ncards = card_count[self.type]
        if ncards:
            self.n = ncards
            float_fmt = self.model.float_fmt
            self.property_id = zeros(ncards, 'int32')
            #: Material ID
            self.material_id = zeros(ncards, 'int32')
            self.thickness = zeros(ncards, float_fmt)
            self.nsm = zeros(ncards, float_fmt)
            self.f1 = zeros(ncards, float_fmt)
            self.f2 = zeros(ncards, float_fmt)

    def add_card(self, card, comment=''):
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
        assert len(card) <= 7, 'len(PSHEAR card) = %i\ncard=%s' % (len(card), card)
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

    def update(self, maps):
        """
        maps = {
            'property' : pid_map,
            'material' : mid_map,
        }
        """
        if self.n:
            pid_map = maps['property']
            mid_map = maps['material']
            for i, pid, mids in enumerate(zip(self.property_id, self.material_ids)):
                self.property_id[i] = pid_map[pid]

    def get_mass_per_area(self, property_id=None):
        # A * (rho * t + nsm)
        if property_id is None:
            property_id = self.property_id
            n = len(property_id)
            i = arange(n)
        else:
            self.model.log.debug('property_id = %r' % property_id)
            n = len(property_id)
            i = searchsorted(self.property_id, property_id)
        mpa = zeros(n, dtype='float64')
        rho = self.model.materials.get_density_by_material_id(self.material_id)
        mpa = rho * self.thickness[i] + self.nsm[i]
        return mpa

    def write_card(self, bdf_file, size=8, pids=None):
        if self.n:
            if pids is None:
                i = arange(self.n)
            else:
                assert len(unique(pids)) == len(pids), unique(pids)
                i = searchsorted(self.property_id, pids)

                #self.property_id = zeros(ncards, 'int32')
                #self.material_id = zeros(ncards, 'int32')
                #self.thickness = zeros(ncards, float_fmt)
                #self.nsm = zeros(ncards, float_fmt)
                #self.f1 = zeros(ncards, float_fmt)
                #self.f2 = zeros(ncards, float_fmt)

            for (pid, mid, t, nsm, f1, f2) in zip(self.property_id[i], self.material_id[i], self.thickness[i], self.nsm[i], self.f1[i], self.f2[i]):
                card = ['PSHEAR', pid, mid, t, nsm, f1, f2]
                if size == 8:
                    bdf_file.write(print_card_8(card))
                else:
                    bdf_file.write(print_card_16(card))

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
        i = self._validate_slice(i)
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
