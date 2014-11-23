
from six.moves import zip

from numpy import array, zeros, unique, searchsorted, arange, asarray

from pyNastran.bdf.dev_vectorized.cards.elements.property import Property

from pyNastran.bdf.fieldWriter import set_blank_if_default
from pyNastran.bdf.fieldWriter import print_card_8
from pyNastran.bdf.bdfInterface.assign_type import (integer, integer_or_blank,
    double, double_or_blank, integer_double_or_blank, integer_string_or_blank,
    string_or_blank, blank)


class PROD(Property):
    type = 'PROD'

    def __init__(self, model):
        """
        Defines the PROD object.

        :param self: the PROD object
        :param model: the BDF object
        """
        Property.__init__(self, model)

    def allocate(self, ncards):
        self.n = ncards
        #print('%s ncards=%s' % (self.type, ncards))
        float_fmt = self.model.float
        #: Property ID
        self.property_id = zeros(ncards, 'int32')
        self.material_id = zeros(ncards, 'int32')
        self.A = zeros(ncards, float_fmt)
        self.J = zeros(ncards, float_fmt)
        self.c = zeros(ncards, float_fmt)
        self.nsm = zeros(ncards, float_fmt)

    def add(self, card, comment):
        i = self.i
        self.property_id[i] = integer(card, 1, 'property_id')
        self.material_id[i] = integer(card, 2, 'material_id')
        self.A[i] = double(card, 3, 'A')
        self.J[i] = double_or_blank(card, 4, 'J', 0.0)
        self.c[i] = double_or_blank(card, 5, 'c', 0.0)
        self.nsm[i] = double_or_blank(card, 6, 'nsm', 0.0)
        assert len(card) <= 7, 'len(PROD card) = %i' % len(card)
        self.i += 1

    def build(self):
        """
        :param self: the PROD object
        :param cards: the list of PROD cards
        """
        if self.n:
            i = self.property_id.argsort()
            self.property_id = self.property_id[i]
            self.material_id = self.material_id[i]
            self.A = self.A[i]
            self.J = self.J[i]
            self.c = self.c[i]
            self.nsm = self.nsm[i]

            unique_pids = unique(self.property_id)
            if len(unique_pids) != len(self.property_id):
                raise RuntimeError('There are duplicate PROD IDs...')
            self._cards = []
            self._comments = []
        else:
            self.property_id = array([], dtype='int32')


    def get_mass_per_length_by_property_id(self, property_id=None):
        # L * (A * rho + nsm)
        if property_id is None:
            i = arange(self.n)
        else:
            i = self.get_index(property_id)
        A = self.A[i]
        mid = self.material_id[i]
        #mat = self.model.materials.get_material(mid)
        rho = self.model.materials.get_density_by_material_id(mid)
        nsm = self.nsm[i]
        return A * rho + nsm

    def get_Area_by_property_id(self, property_id):
        i = self.get_index(property_id)
        A = self.A[i]
        return A

    def get_non_structural_mass_by_property_id(self, property_id):
        i = self.get_index(property_id)
        nsm = self.nsm[i]
        return nsm

    def get_E(self, property_id):
        i = self.get_index(property_id)
        material_id = self.material_id[i]
        E = self.model.materials.get_E_by_material_id(material_id)
        return E

    def get_G_by_property_id(self, property_id):
        i = self.get_index(property_id)
        material_id = self.material_id[i]
        G = self.model.materials.get_G_by_material_id(material_id)
        return G

    def get_J_by_property_id(self, property_id):
        i = self.get_index(property_id)
        J = self.J[i]
        return J

    def get_c_by_property_id(self, property_id):
        i = self.get_index(property_id)
        c = self.c[i]
        return c

    def get_index(self, property_id):
        if isinstance(property_id, int):
            property_id = array([property_id])
        if property_id is None:
            return arange(self.n)

        indexs = searchsorted(self.property_id, property_id)
        assert len(indexs) == len(property_id), 'indexs=%s pids=%s' % (indexs, property_id)
        return indexs

    def write_bdf(self, f, size=8, property_id=None):
        if self.n:
            if self.n:
                if property_id is None:
                    i = arange(self.n)
                else:
                    assert len(unique(property_id))==len(property_id), unique(property_id)
                    i = searchsorted(self.property_id, property_id)
            for (pid, mid, A, J, c, nsm) in zip(
                 self.property_id, self.material_id[i], self.A[i], self.J[i], self.c[i], self.nsm[i]):

                #self.mid = integer(card, 4, 'mid')
                #self.A = double(card, 5, 'A')
                #self.j = double_or_blank(card, 6, 'j', 0.0)
                #self.c = double_or_blank(card, 7, 'c', 0.0)
                #self.nsm = double_or_blank(card, 8, 'nsm', 0.0)

                card = ['PROD', pid, mid, A, J, c, nsm]
                f.write(print_card_8(card))

    def slice_by_index(self, i):
        i = asarray(i)
        obj = PROD(self.model)
        obj.n = len(i)
        #obj._cards = self._cards[i]
        #obj._comments = obj._comments[i]
        #obj.comments = obj.comments[i]
        obj.property_id = self.property_id[i]
        obj.material_id = self.material_id[i]
        obj.A = self.A[i]
        obj.J = self.J[i]
        obj.c = self.c[i]
        obj.nsm = self.nsm[i]
        return obj
