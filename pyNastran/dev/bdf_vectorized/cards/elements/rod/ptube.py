from numpy import array, zeros, unique, searchsorted, arange, pi

from pyNastran.dev.bdf_vectorized.cards.elements.property import Property

#from pyNastran.bdf.field_writer_8 import set_blank_if_default
from pyNastran.bdf.field_writer_8 import print_card_8
from pyNastran.bdf.bdf_interface.assign_type import (integer,
    double, double_or_blank)
from pyNastran.bdf.bdf_interface.bdf_card import BDFCard


class PTUBE(Property):
    type = 'PTUBE'

    def __init__(self, model):
        """
        Defines the PTUBE object.

        Parameters
        ----------
        model : BDF
           the BDF object

        """
        Property.__init__(self, model)

    def allocate(self, ncards):
        self.n = ncards
        self.model.log.debug('%s ncards=%s' % (self.type, ncards))
        float_fmt = self.model.float_fmt
        #: Property ID
        self.property_id = zeros(ncards, 'int32')
        self.material_id = zeros(ncards, 'int32')

        self.OD = zeros((ncards, 2), float_fmt)
        self.t = zeros(ncards, float_fmt)
        self.nsm = zeros(ncards, float_fmt)

    def add_card(self, card: BDFCard, comment: str=''):
        self.model.log.debug('n=%s i=%s' % (self.n, self.i))
        i = self.i
        self.property_id[i] = integer(card, 1, 'property_id')
        self.material_id[i] = integer(card, 2, 'material_id')
        OD1 = double(card, 3, 'OD1')
        t = double_or_blank(card, 4, 't')
        if t is None:
            t = OD1 / 2.
        self.t[i] = t
        self.nsm[i] = double_or_blank(card, 5, 'nsm', 0.0)
        OD2 = double_or_blank(card, 6, 'OD2', OD1)
        self.OD[i, :] = [OD1, OD2]
        assert len(card) <= 7, 'len(PTUBE card) = %i\ncard=%s' % (len(card), card)
        self.i += 1

    def build(self):
        """
        :param cards: the list of PTUBE cards
        """
        if self.n:
            i = self.property_id.argsort()
            self.property_id = self.property_id[i]
            self.material_id = self.material_id[i]
            self.OD = self.OD[i, :]
            self.t = self.t[i]
            self.nsm = self.nsm[i]

            unique_pids = unique(self.property_id)
            if len(unique_pids) != len(self.property_id):
                raise RuntimeError('There are duplicate PTUBE IDs...')
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
            nid_map = maps['node']
            mid_map = maps['material']
            for i, (pid, mid) in enumerate(zip(self.property_id, self.material_id)):
                self.property_id[i] = pid_map[pid]
                self.material_id[i] = mid_map[mid]

    #=========================================================================

    def get_mass_per_length_by_property_id(self, property_id=None):
        # L * (A * rho + nsm)
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
        area = zeros(len(i), dtype='float64')
        for ni, ii in enumerate(i):
            A = (self._area1(ii) + self._area2(ii)) / 2.
            area[ni] = A
        return area

    def _area1(self, i):
        """Gets the Area of Section 1 of the CTUBE."""
        Dout = self.OD[i, 0]
        if self.t[i] == 0:
            return pi / 4. * Dout**2
        Din = Dout - 2 * self.t
        A1 = pi / 4. * (Dout * Dout - Din * Din)
        return A1

    def _area2(self, i):
        """Gets the Area of Section 2 of the CTUBE."""
        Dout = self.OD[i, 1]
        if self.t[i] == 0:
            return pi / 4. * Dout**2
        Din = Dout - 2 * self.t
        A2 = pi / 4. * (Dout * Dout - Din * Din)
        return A2

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
        J = []
        for ni, ii in enumerate(i):
            Ji = self._Ji(ii)
            J.append(Ji)
        return array(J, dtype='float64')

    def _Ji(self, i):
        Dout = self.OD[i, 0]
        if self.t[0] == 0.0:
            return pi / 8. * Dout**4
        Din = Dout - 2 * self.t[i]
        return pi / 8. * (Dout**4 - Din**2)

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

            for (pid, mid, (OD1, OD2), t, nsm) in zip(
                 self.property_id, self.material_id[i], self.OD[i, :], self.t[i], self.nsm[i]):

                #t = set_blank_if_default(t, OD1 / 2.)
                #nsm = set_blank_if_default(nsm, 0.0)
                #OD2 = set_blank_if_default(OD2, OD1)

                card = ['PTUBE', pid, mid, OD1, t, nsm, OD2]
                bdf_file.write(print_card_8(card))

    def slice_by_index(self, i):
        i = self._validate_slice(i)
        obj = PTUBE(self.model)
        n = len(i)
        obj.n = n
        obj.i = n
        #obj._cards = self._cards[i]
        #obj._comments = obj._comments[i]
        #obj.comments = obj.comments[i]
        obj.property_id = self.property_id[i]
        obj.material_id = self.material_id[i]
        obj.OD = self.OD[i, :]
        obj.t = self.t[i]
        obj.nsm = self.nsm[i]
        return obj
