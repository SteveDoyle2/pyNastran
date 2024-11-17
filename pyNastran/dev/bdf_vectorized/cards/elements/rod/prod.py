from numpy import array, zeros, unique, searchsorted, arange

from pyNastran.dev.bdf_vectorized.cards.elements.property import Property

from pyNastran.bdf.field_writer_8 import print_card_8
from pyNastran.bdf.bdf_interface.assign_type import (integer,
    double, double_or_blank)


class PROD(Property):
    type = 'PROD'

    def __init__(self, model):
        """
        Defines the PROD object.

        Parameters
        ----------
        model : BDF
           the BDF object
        """

        Property.__init__(self, model)

    def allocate(self, card_count):
        ncards = card_count[self.type]
        if ncards:
            self.n = ncards
            self.model.log.debug('%s ncards=%s' % (self.type, ncards))
            float_fmt = self.model.float_fmt
            #: Property ID
            self.property_id = zeros(ncards, 'int32')
            self.material_id = zeros(ncards, 'int32')
            self.A = zeros(ncards, float_fmt)
            self.J = zeros(ncards, float_fmt)
            self.c = zeros(ncards, float_fmt)
            self.nsm = zeros(ncards, float_fmt)

    def add_card(self, card: BDFCard, comment: str=''):
        #self.model.log.debug('n=%s i=%s' % (self.n, self.i))
        i = self.i
        self.property_id[i] = integer(card, 1, 'property_id')
        self.material_id[i] = integer(card, 2, 'material_id')
        self.A[i] = double(card, 3, 'A')
        self.J[i] = double_or_blank(card, 4, 'J', 0.0)
        self.c[i] = double_or_blank(card, 5, 'c', 0.0)
        self.nsm[i] = double_or_blank(card, 6, 'nsm', 0.0)
        assert len(card) <= 7, 'len(PROD card) = %i\ncard=%s' % (len(card), card)
        self.i += 1

    def build(self):
        """
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

    #=========================================================================

    def get_mass_per_length_by_property_id(self, property_id=None):
        # L * (A * rho + nsm)
        if property_id is None:
            i = arange(self.n)
        else:
            i = self.get_property_index_by_property_id(property_id)
        A = self.A[i]
        mid = self.material_id[i]
        #mat = self.model.materials.get_material(mid)
        rho = self.model.materials.get_density_by_material_id(mid)
        nsm = self.nsm[i]
        return A * rho + nsm

    def get_area_by_property_id(self, property_id=None):
        i = self.get_property_index_by_property_id(property_id)
        return self.get_area_by_property_index(i)

    def get_area_by_property_index(self, i=None):
        A = self.A[i]
        return A

    def get_non_structural_mass_by_property_id(self, property_id=None):
        i = self.get_property_index_by_property_id(property_id)
        nsm = self.nsm[i]
        return nsm

    #def get_E_by_property_id(self, property_id=None):
        #i = self.get_property_index_by_property_id(property_id)
        #material_id = self.material_id[i]
        #E = self.model.materials.get_E_by_material_id(material_id)
        #return E

    def get_E_by_property_id(self, property_id=None):
        mid = self.get_material_id_by_property_id(property_id)
        E = self.model.materials.get_E_by_material_id(mid)
        return E

    #def get_G_by_property_id(self, property_id=None):
        #i = self.get_property_index_by_property_id(property_id)
        #material_id = self.material_id[i]
        #G = self.model.materials.get_G_by_material_id(material_id)
        #return G

    def get_G_by_property_id(self, property_id=None):
        mid = self.get_material_id_by_property_id(property_id)
        G = self.model.materials.get_G_by_material_id(mid)
        return G

    def get_J_by_property_id(self, property_id=None):
        i = self.get_property_index_by_property_id(property_id)
        return self.get_J_by_property_index(i)

    def get_J_by_property_index(self, i=None):
        J = self.J[i]
        return J

    def get_c_by_property_id(self, property_id=None):
        i = self.get_property_index_by_property_id(property_id)
        c = self.c[i]
        return c

    def get_material_id_by_property_id(self, property_id=None):
        i = self.get_property_index_by_property_id(property_id)
        mid = self.material_id[i]
        return mid

    #=========================================================================

    def get_density_by_property_id(self, property_id=None):
        mid = self.get_material_id_by_property_id(property_id)
        density = self.model.materials.get_density_by_material_id(mid)
        return density

    #def get_J_by_property_id(self, property_id=None):
        #mid = self.get_material_id_by_property_id(property_id)
        #J = self.model.materials.get_J_by_material_id(mid)
        #return J

    #def get_E_by_property_id(self, property_id=None):
        #mid = self.get_material_id_by_property_id(property_id)
        #E = self.model.materials.get_E_by_material_id(mid)
        #return E

    #=========================================================================

    def write_card(self, bdf_file, size=8, property_id=None):
        if self.n:
            if self.n:
                if property_id is None:
                    i = arange(self.n)
                else:
                    assert len(unique(property_id)) == len(property_id), unique(property_id)
                    i = searchsorted(self.property_id, property_id)
            for (pid, mid, A, J, c, nsm) in zip(
                 self.property_id, self.material_id[i], self.A[i], self.J[i], self.c[i], self.nsm[i]):

                #self.mid = integer(card, 4, 'mid')
                #self.A = double(card, 5, 'A')
                #self.j = double_or_blank(card, 6, 'j', 0.0)
                #self.c = double_or_blank(card, 7, 'c', 0.0)
                #self.nsm = double_or_blank(card, 8, 'nsm', 0.0)

                card = ['PROD', pid, mid, A, J, c, nsm]
                bdf_file.write(print_card_8(card))

    def slice_by_index(self, i):
        i = self._validate_slice(i)
        obj = PROD(self.model)
        n = len(i)
        obj.n = n
        obj.i = n
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
