from numpy import array, zeros, arange, searchsorted, unique

from pyNastran.utils.numpy_utils import integer_types

from pyNastran.bdf.field_writer_8 import print_card_8
from pyNastran.bdf.bdf_interface.bdf_card import BDFCard
from pyNastran.bdf.bdf_interface.assign_type import (integer,
    double_or_blank)

from pyNastran.dev.bdf_vectorized.cards.elements.property import Property


class PBAR(Property):
    type = 'PBAR'

    def __init__(self, model):
        """
        Defines the PBAR object.

        Parameters
        ----------
        model : BDF
           the BDF object
        """
        Property.__init__(self, model)

    #def add1(self, card, comment=''):
        #self._cards.append(card)
        #self._comments.append(comment)

    def allocate(self, card_count):
        ncards = card_count[self.type]
        if ncards:
            self.n = ncards
            float_fmt = self.model.float_fmt
            #: Property ID
            self.property_id = zeros(ncards, 'int32')
            #: material ID
            self.material_id = zeros(ncards, 'int32')
            #: Area
            self.area = zeros(ncards, float_fmt)
            #: I1
            self.I1 = zeros(ncards, float_fmt)
            #: I2
            self.I2 = zeros(ncards, float_fmt)

            #: Polar Moment of Inertia J -> use J()
            #: default=1/2(I1+I2) for SOL=600, otherwise 0.0
            self.J = zeros(ncards, float_fmt)
            self.nsm = zeros(ncards, float_fmt)

            self.C1 = zeros(ncards, float_fmt)
            self.C2 = zeros(ncards, float_fmt)
            self.D1 = zeros(ncards, float_fmt)
            self.D2 = zeros(ncards, float_fmt)
            self.E1 = zeros(ncards, float_fmt)
            self.E2 = zeros(ncards, float_fmt)
            self.F1 = zeros(ncards, float_fmt)
            self.F2 = zeros(ncards, float_fmt)

            #: default=infinite; assume 1e8
            self.K1 = zeros(ncards, float_fmt)

            #: default=infinite; assume 1e8
            self.K2 = zeros(ncards, float_fmt)

            #: I12 -> use I12()
            self.i12 = zeros(ncards, float_fmt)


    def add_card(self, card: BDFCard, comment: str=''):
        i = self.i

        self.property_id[i] = integer(card, 1, 'property_id')
        self.material_id[i] = integer(card, 2, 'material_id')
        self.area[i] = double_or_blank(card, 3, 'area', 0.0)
        self.I1[i] = double_or_blank(card, 4, 'I1', 0.0)
        self.I2[i] = double_or_blank(card, 5, 'I2', 0.0)

        #: .. todo:: support SOL 600 default
        Jdefault = 0.5 * (self.I1[i] + self.I2[i])
        self.J[i] = double_or_blank(card, 6, 'J', Jdefault)
        self.nsm[i] = double_or_blank(card, 7, 'non-structural_mass', 0.0)

        if 1:
            #c1 = double_or_blank(card, 9, 'C1', 0.0)
            #c2 = double_or_blank(card, 10, 'C2', 0.0)
            #d1 = double_or_blank(card, 11, 'D1', 0.0)
            #d2 = double_or_blank(card, 12, 'D2', 0.0)
            #e1 = double_or_blank(card, 13, 'E1', 0.0)
            #e2 = double_or_blank(card, 14, 'E2', 0.0)
            #f1 = double_or_blank(card, 15, 'F1', 0.0)
            #f2 = double_or_blank(card, 16, 'F2', 0.0)

            #i12 = double_or_blank(card, 19, 'I12', 0.0)

            #if A == 0.0:
                #k1 = blank(card, 17, 'K1')
                #k2 = blank(card, 18, 'K2')
            #elif i12 != 0.0:
                ## K1 / K2 are ignored
                #k1 = None
                #k2 = None
            #else:
                ##: default=infinite; assume 1e8
                #k1 = double_or_blank(card, 17, 'K1', 1e8)
                ##: default=infinite; assume 1e8
                #k2 = double_or_blank(card, 18, 'K2', 1e8)

            self.C1[i] = double_or_blank(card, 9, 'C1', 0.0)
            self.C2[i] = double_or_blank(card, 10, 'C2', 0.0)
            self.D1[i] = double_or_blank(card, 11, 'D1', 0.0)
            self.D2[i] = double_or_blank(card, 12, 'D2', 0.0)
            self.E1[i] = double_or_blank(card, 13, 'E1', 0.0)
            self.E2[i] = double_or_blank(card, 14, 'E2', 0.0)
            self.F1[i] = double_or_blank(card, 15, 'F1', 0.0)
            self.F2[i] = double_or_blank(card, 16, 'F2', 0.0)

            self.K1[i] = double_or_blank(card, 17, 'K1', 1e8)
            self.K2[i] = double_or_blank(card, 18, 'K2', 1e8)
            self.i12[i] = double_or_blank(card, 19, 'I12', 0.0)
            #if self.A == 0.0 and self.i12 == 0.0:
                #assert self.K1 is None, 'K1 must be blank if A=0.0 and I12=0.0; A=%r I12=%r K1=%r' % (self.A, self.i12, self.K1)
                #assert self.K2 is None, 'K2 must be blank if A=0.0 and I12=0.0; A=%r I12=%r K2=%r' % (self.A, self.i12, self.K2)
            #assert len(card) <= 20, 'len(PBAR card) = %i\ncard=%s' % (len(card), card)

        self.i += 1

    def build(self):
        if self.n:
            i = self.property_id.argsort()
            self.property_id = self.property_id[i]
            self.material_id = self.material_id[i]

            self.area = self.area[i]
            self.I1 = self.I1[i]
            self.I2 = self.I2[i]
            self.J = self.J[i]
            self.nsm = self.nsm[i]

            unique_pids = unique(self.property_id)
            if len(unique_pids) != len(self.property_id):
                raise RuntimeError('There are duplicate PCOMP IDs...')
        else:
            self.property_id = array([], dtype='int32')

    def update(self, maps):
        """
        maps = {
            'node' : nid_map,
            'property' : pid_map,
        }
        """
        if self.n:
            nid_map = maps['node']
            pid_map = maps['property']
            for i, pid in enumerate(self.property_id):
                self.property_id[i] = pid_map[pid]

    #=========================================================================
    def get_index(self, property_id):
        if isinstance(property_id, integer_types):
            property_ids = array([property_id])
        if property_ids is None:
            return arange(self.n)

        indexs = searchsorted(self.property_id, property_id)
        assert len(indexs) == len(property_id), 'indexs=%s pids=%s' % (indexs, property_id)
        return indexs

    #=========================================================================
    def write_card(self, bdf_file, size=8, property_id=None):
        if self.n:
            if property_id is None:
                i = arange(self.n)
            else:
                i = searchsorted(self.property_id, property_id)

            for (pid, mid, area, I1, I2, J, nsm, c1, c2, d1, d2, e1, e2, f1, f2, k1, k2,
                 i12) in zip(self.property_id[i], self.material_id[i],
                    self.area[i], self.I1[i], self.I2[i], self.J[i], self.nsm[i],
                    self.C1[i], self.C2[i], self.D1[i], self.D2[i],
                    self.E1[i], self.E2[i], self.F1[i], self.F2[i],
                    self.K1[i], self.K2[i], self.i12[i]):
                if pid in self._comments:
                    bdf_file.write(self._comments[pid])

                #card = ['PBAR', pid, mid, area, I1, I2, J]
                #c1 = set_blank_if_default(self.c1, 0.0)
                #c2 = set_blank_if_default(self.c2, 0.0)

                #d1 = set_blank_if_default(self.d1, 0.0)
                #d2 = set_blank_if_default(self.d2, 0.0)

                #e1 = set_blank_if_default(self.e1, 0.0)
                #e2 = set_blank_if_default(self.e2, 0.0)

                #f1 = set_blank_if_default(self.f1, 0.0)
                #f2 = set_blank_if_default(self.f2, 0.0)

                #k1 = set_blank_if_default(self.k1, 1e8)
                #k2 = set_blank_if_default(self.k2, 1e8)

                list_fields = ['PBAR', pid, mid, area, I1, I2, J, nsm,
                               None, c1, c2, d1, d2, e1, e2, f1, f2, k1, k2, i12]
                bdf_file.write(print_card_8(list_fields))

    def get_mass_per_length_by_property_id(self, property_id=None):
        if property_id is None:
            i = arange(self.n)
        else:
            i = searchsorted(self.property_id, property_id)
        A = self.area[i]
        mid = self.material_id[i]
        rho = self.model.materials.get_density_by_material_id(mid)
        nsm = self.nsm[i]
        return rho * A + nsm

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

    def slice_by_index(self, i):
        i = self._validate_slice(i)
        #print('pbar i = %s' % i)
        obj = PBAR(self.model)
        obj.n = len(i)
        #obj._cards = self._cards[i]
        #obj._comments = obj._comments[i]
        #obj.comments = obj.comments[i]
        #print(self)
        obj.property_id = self.property_id[i]
        obj.material_id = self.material_id[i]
        obj.area = self.area[i]
        obj.I1 = self.I1[i]
        obj.I2 = self.I2[i]
        obj.J = self.J[i]
        obj.nsm = self.nsm[i]
        return obj
