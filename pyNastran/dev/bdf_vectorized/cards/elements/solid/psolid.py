from io import StringIO
from numpy import zeros, unique, where, searchsorted, asarray, array

from pyNastran.utils.numpy_utils import integer_types
from pyNastran.bdf.field_writer_8 import set_blank_if_default
from pyNastran.bdf.field_writer_8 import print_card_8
from pyNastran.bdf.bdf_interface.assign_type import (integer, integer_or_blank,
    integer_string_or_blank,
    string_or_blank)

from pyNastran.dev.bdf_vectorized.cards.elements.property import Property

class PSOLID(Property):
    type = 'PSOLID'
    def __init__(self, model):
        """
        Defines the PSOLID object.

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
            #float_fmt = self.model.float_fmt
            #: Property ID
            self.property_id = zeros(ncards, 'int32')
            #: Material ID
            self.material_id = zeros(ncards, 'int32')
            self.cordm = zeros(ncards, 'int32')
            self.integ = zeros(ncards, dtype='|U8')
            self.stress = zeros(ncards, dtype='|U8')
            self.isop = zeros(ncards, dtype='|U8')
            self.fctn = zeros(ncards, dtype='|U8')

    def add_card(self, card: BDFCard, comment: str=''):
        i = self.i
        self.property_id[i] = integer(card, 1, 'pid')
        self.material_id[i] = integer(card, 2, 'mid')
        self.cordm[i] = integer_or_blank(card, 3, 'cordm', 0)
        self.integ[i] = integer_string_or_blank(card, 4, 'integ', '')
        #validIntegration = ['THREE', 'TWO', 'FULL', 'BUBBLE',
        #                    2, 3, None, 'REDUCED']
        # ISOP
        # ------
        #    1.  FULL
        #    2.
        #    3.
        #    REDUCED

        # IN
        # ------
        #    1.
        #    2.      TWO
        #    3.      THREE
        #    BUBBLE - 2 for CTETRA, 3 for CHEXA/CPENTA

        # STRESS
        # ------
        #    1.  GAUSS (no midside nodes on CPENTA/CHEXA; ok on CTETRA)
        #    2.
        self.stress[i] = integer_string_or_blank(card, 5, 'stress', '')
        self.isop[i] = integer_string_or_blank(card, 6, 'isop', '')
        self.fctn[i] = string_or_blank(card, 7, 'fctn', 'SMECH')
        assert len(card) <= 8, 'len(PSOLID card) = %i\ncard=%s' % (len(card), card)
        self.i += 1

    def build(self):
        """
        :param cards: the list of PSOLID cards
        """
        #print("N[%s] = %s" % (self.type, self.n))
        if self.n:
            i = self.property_id.argsort()
            self.property_id = self.property_id[i]
            #print("PSOLID.property_id =", self.property_id)
            self.material_id = self.material_id[i]
            self.cordm = self.cordm[i]
            self.integ = self.integ[i]
            self.stress = self.stress[i]
            self.isop = self.isop[i]
            self.fctn = self.fctn[i]

            unique_pids = unique(self.property_id)
            if len(unique_pids) != len(self.property_id):
                raise RuntimeError('There are duplicate PSOLID IDs...')
        else:
            self.property_id = array([], dtype='int32')
            self.material_id = array([], dtype='int32')

    def update(self, maps):
        """
        maps = {
            'property' : pid_map,
            'material' : mid_map,
        }
        """
        if self.n:
            #nid_map = maps['node']
            pid_map = maps['property']
            mid_map = maps['material']
            for i, pid, mids in enumerate(zip(self.property_id, self.material_ids)):
                self.property_id[i] = pid_map[pid]

    def get_density_by_property_id(self, property_id=None):
        if property_id is None:
            property_id = self.property_id
        elif isinstance(property_id, integer_types):
            property_id = array([property_id], dtype='int32')
        property_id = asarray(property_id)
        n = len(property_id)
        rho = zeros(n, dtype='float64')
        upid = unique(property_id)
        for i, pid in enumerate(upid):
            # get the location of pid id in global props
            j = where(pid == self.property_id)[0]

            # get the location of pid in the local props
            k = where(pid == property_id)[0]

            if len(j) == 0:
                msg = 'pid=%s was not found in %s' % (upid, self.property_id)
                raise ValueError(msg)
            mid = self.material_id[j[0]]
            rhoi = self.model.materials.get_density_by_material_id([mid])
            #print('pid=%s rho[%s]=%s' % (pid, k, rhoi))
            rho[k] = rhoi

        if rho.min() == 0.0:
            msg = 'property_id = %s\n' % property_id
            msg += 'i = %s\n' % i
            msg += 'rho = %s\n' % rho
            msg += 'rhoj = %s' % rhoi
            raise ValueError(msg)

        assert rho.shape == (n, ), rho.shape
        self.model.log.debug('rho = %s' % rho)
        return rho

    def write_card(self, bdf_file, size=8, property_id=None):
        if self.n:
            #print("PSOLID.property_id =", self.property_id)
            for(pid, mid, cordm, integ, stress, isop, fctn) in zip(
                self.property_id, self.material_id, self.cordm,
                self.integ, self.stress, self.isop, self.fctn):
                if pid in self._comments:
                    bdf_file.write(self._comments[pid])

                cordm = set_blank_if_default(cordm, 0)
                fctn = set_blank_if_default(fctn, 'SMECH')
                card = ['PSOLID', pid, mid, cordm, integ,
                        stress, isop, fctn]
                bdf_file.write(print_card_8(card))

    def __getitem__(self, property_id):
        i = searchsorted(self.property_id, property_id)
        return self.slice_by_index(i)

    def slice_by_index(self, i):
        i = asarray(i)
        obj = PSOLID(self.model)
        obj.n = len(i)
        #obj._cards = self._cards[i]
        #obj._comments = obj._comments[i]
        #obj.comments = obj.comments[i]
        obj.property_id = self.property_id[i]
        obj.material_id = self.material_id[i]
        obj.cordm = self.cordm[i]
        obj.integ = self.integ[i]
        obj.stress = self.stress[i]
        obj.isop = self.isop[i]
        obj.fctn = self.fctn[i]
        return obj

    def __repr__(self):
        f = StringIO()
        f.write('<PSOLID object> n=%s\n' % self.n)
        self.write_card(f)
        #print f
        return f.getvalue()
