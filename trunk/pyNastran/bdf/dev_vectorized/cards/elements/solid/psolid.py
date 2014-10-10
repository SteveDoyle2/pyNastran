import cStringIO
from numpy import zeros, unique, where

from pyNastran.bdf.fieldWriter import set_blank_if_default
from pyNastran.bdf.fieldWriter import print_card_8
from pyNastran.bdf.bdfInterface.assign_type import (integer, integer_or_blank,
    double_or_blank, integer_double_or_blank, integer_string_or_blank,
    string_or_blank, blank)

class PSOLID(object):
    type = 'PSOLID'
    def __init__(self, model):
        """
        Defines the PSOLID object.

        :param self: the PSOLID object
        :param model: the BDF object
        """
        self.model = model
        self.n = 0
        self._cards = []
        self._comments = []

    def add(self, card, comment):
        self._cards.append(card)
        self._comments.append(comment)

    def build(self):
        """
        :param self: the PSOLID object
        :param cards: the list of PSOLID cards
        """

        cards = self._cards
        ncards = len(cards)
        self.n = ncards
        #print "N[%s] = %s" % (self.type, self.n)
        if ncards:
            float_fmt = self.model.float
            #: Property ID
            self.property_id = zeros(ncards, 'int32')
            #: Material ID
            self.material_id = zeros(ncards, 'int32')
            self.cordm = zeros(ncards, 'int32')
            self.integ = zeros(ncards, dtype='|S8')
            self.stress = zeros(ncards, dtype='|S8')
            self.isop = zeros(ncards, dtype='|S8')
            self.fctn = zeros(ncards, dtype='|S8')
            #print "isop", self.isop

            for i, card in enumerate(cards):
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
                assert len(card) <= 8, 'len(PSOLID card) = %i' % len(card)

            i = self.property_id.argsort()
            self.property_id = self.property_id[i]
            #print "PSOLID.property_id =", self.property_id
            self.material_id = self.material_id[i]
            self.cordm = self.cordm[i]
            self.integ = self.integ[i]
            self.stress = self.stress[i]
            self.isop = self.isop[i]
            self.fctn = self.fctn[i]

            unique_pids = unique(self.property_id)
            if len(unique_pids) != len(self.property_id):
                raise RuntimeError('There are duplicate PSOLID IDs...')
            self._cards = []
            self._comments = []
        else:
            self.property_id = array([], dtype='int32')
            aaa

    def get_density(self, property_ids=None):
        if property_ids is None:
            property_ids = self.property_id

        rho = []
        i = where(property_ids == self.property_id)[0]
        for mid in self.material_id[i]:
            rhoi = self.model.materials.mat1[mid].get_density()
            rho.append(rhoi)
        return rho

    def write_bdf(self, f, size=8, property_ids=None):
        if self.n:
            #print "PSOLID.property_id =", self.property_id
            for (pid, mid, cordm, integ, stress, isop, fctn) in zip(
                 self.property_id, self.material_id, self.cordm,
                 self.integ, self.stress, self.isop, self.fctn):

                cordm = set_blank_if_default(cordm, 0)
                fctn = set_blank_if_default(self.fctn, 'SMECH')
                card = ['PSOLID', pid, mid, cordm, integ,
                        stress, isop, fctn]
                #print card
                f.write(print_card_8(card))

    def __repr__(self):
        f = cStringIO.StringIO()
        f.write('<PSOLID object> n=%s\n' % self.n)
        self.write_bdf(f)
        #print f
        return f.getvalue()

