import cStringIO
from numpy import array, zeros, arange, searchsorted, where, unique, hstack, concatenate

from pyNastran.bdf.dev_vectorized.utils import slice_to_iter

from pyNastran.bdf.fieldWriter import print_card_8
from pyNastran.bdf.fieldWriter16 import print_card_16
from pyNastran.bdf.bdfInterface.assign_type import (integer, integer_or_blank,
    double_or_blank, integer_double_or_blank, blank)
from pyNastran.bdf.cards.properties.shell import PCOMP as vPCOMP

class PCOMP(object):
    type = 'PCOMP'
    def __init__(self, model):
        """
        Defines the PCOMP object.

        :param self: the PCOMP object
        :param model: the BDF object
        :param cards: the list of PCOMP cards
        """
        self.model = model
        self.properties = {}
        #self._cards = []
        #self._comments = []

    def add(self, card, comment):
        prop = vPCOMP(card, comment=comment)
        self.properties[prop.pid] = prop
        #self._cards.append(card)
        #self._comments.append(comment)

    def get_thickness(self, property_id=None):
        props = self.get_properties(property_id)
        return array([prop.get_thickness() for prop in props])

    def get_nonstructural_mass(self, property_id=None):
        props = self.get_properties(property_id)
        d = array([prop.get_nonstructural_mass() for prop in props])
        return d

    def get_material_ids(self, property_id=None):
        props = self.get_properties(property_id)
        d = concatenate([prop.get_material_ids() for prop in props])
        #print('mids PCOMP = %s' % d)
        return d

    def get_mass_per_area(self, property_id=None):
        #massPerArea = self.nsm + self.Rho() * self.t

        props = self.get_properties(property_id)
        mids_all = unique(self.get_material_ids(property_id))
        mids_all.sort()
        densities = self.model.materials.get_density(mids_all)

        #print('nprops =', len(props))
        mpa = zeros(len(props), dtype='float64')
        for i, prop in enumerate(props):
            j = searchsorted(mids_all, prop.get_material_ids())
            densities2 = densities[j]
            #print('nsm = %s' % prop.get_nonstructural_mass())
            #print('rhos = %s' % densities2)
            mpai = prop.get_mass_per_area_rho(densities2)
            #print('mpai = %s' % mpai)
            mpa[i] = mpai
        return mpa

    def get_properties(self, property_id=None):
        props = []
        if property_id is None:
            property_id = self.property_id
        for pid in self.properties:
            prop = self.properties[pid]
            props.append(prop)
        return props

    def build(self):
        self.n = len(self.properties)
        self.property_id = sorted(self.properties.keys())
        return
        self.n = 0
        cards = self._cards
        ncards = len(cards)
        self.n = ncards
        #return
        if ncards:
            #: Property ID
            self.property_id = zeros(ncards, 'int32')

            ncards = len(cards)
            for i, card in enumerate(cards):
                self.property_id[i] = integer(card, 1, 'pid')
            i = self.property_id.argsort()
            self.property_id = self.property_id[i]
            unique_pids = unique(self.property_id)

            if len(unique_pids) != len(self.property_id):
                raise RuntimeError('There are duplicate PCOMP IDs...')
            self._cards = []
            self._comments = []

    def write_bdf(self, f, size=8, property_ids=None):
        if size == 8:
            for pid, pcomp in sorted(self.properties.iteritems()):
                f.write(pcomp.write_bdf(size, print_card_8))
        else:
            for pid, pcomp in sorted(self.properties.iteritems()):
                f.write(pcomp.write_bdf(size, print_card_16))


    #def slice_by_index(self, i):
        #i = asarray(i)
        #return self.__getitem__()


    def __getitem__(self, property_ids):
        property_ids, int_flag = slice_to_iter(property_ids)
        obj = PCOMP(self.model)

        properties = {}
        for pid in sorted(property_ids):
            properties[pid] = self.properties[pid]
        obj.n = len(property_ids)
        obj.properties = properties
        obj.property_id = sorted(self.properties.keys())
        #obj._comments = obj._comments[index]
        #obj.comments = obj.comments[index]
        return obj


    def __repr__(self):
        f = cStringIO.StringIO()
        f.write('<PCOMP object> n=%s\n' % self.n)
        self.write_bdf(f)
        return f.getvalue()
