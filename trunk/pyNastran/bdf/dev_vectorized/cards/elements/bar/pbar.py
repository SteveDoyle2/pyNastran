from numpy import array, zeros, arange, concatenate, searchsorted, where, unique

from pyNastran.bdf.fieldWriter import print_card_8
from pyNastran.bdf.bdfInterface.assign_type import (integer, integer_or_blank,
    double_or_blank, integer_double_or_blank, blank)


class PBAR(object):
    type = 'PBAR'
    def __init__(self, model):
        """
        Defines the PCOMP object.

        :param self: the PCOMP object
        :param model: the BDF object
        :param cards: the list of PCOMP cards
        """
        self.model = model
        self.n = 0
        self._cards = []
        self._comments = []

    def add(self, card, comment):
        self._cards.append(card)
        self._comments.append(comment)

    def build(self):
        cards = self._cards
        ncards = len(cards)
        self.n = ncards

        if ncards:
            #: Property ID
            self.property_id = zeros(ncards, 'int32')
            self.material_id = zeros(ncards, 'int32')
            self.area = zeros(ncards, 'float64')
            self.I1 = zeros(ncards, 'float64')
            self.I2 = zeros(ncards, 'float64')
            self.J = zeros(ncards, 'float64')
            self.nsm = zeros(ncards, 'float64')

            for i, card in enumerate(cards):
                #: property ID
                self.property_id[i] = integer(card, 1, 'property_id')

                #: material ID
                self.material_id[i] = integer(card, 2, 'material_id')


                #: material ID
                self.area[i] = double_or_blank(card, 3, 'area', 0.0)

                #: I1
                self.I1[i] = double_or_blank(card, 4, 'I1', 0.0)

                #: I2
                self.I2[i] = double_or_blank(card, 5, 'I2', 0.0)

                #: Polar Moment of Inertia J -> use J()
                #: default=1/2(I1+I2) for SOL=600, otherwise 0.0
                #: .. todo:: support SOL 600 default

                Jdefault = 0.5 * (self.I1[i] + self.I2[i])
                self.J[i] = double_or_blank(card, 6, 'J', Jdefault)
                self.nsm[i] = double_or_blank(card, 7, 'non-structural_mass', 0.0)

                if 0:
                    self.C1 = double_or_blank(card, 9, 'C1', 0.0)
                    self.C2 = double_or_blank(card, 10, 'C2', 0.0)
                    self.D1 = double_or_blank(card, 11, 'D1', 0.0)
                    self.D2 = double_or_blank(card, 12, 'D2', 0.0)
                    self.E1 = double_or_blank(card, 13, 'E1', 0.0)
                    self.E2 = double_or_blank(card, 14, 'E2', 0.0)
                    self.F1 = double_or_blank(card, 15, 'F1', 0.0)
                    self.F2 = double_or_blank(card, 16, 'F2', 0.0)

                    #: default=infinite; assume 1e8
                    self.K1 = double_or_blank(card, 17, 'K1', 1e8)
                    #: default=infinite; assume 1e8
                    self.K2 = double_or_blank(card, 18, 'K2', 1e8)
                    #: I12 -> use I12()
                    self.i12 = double_or_blank(card, 19, 'I12', 0.0)
                    if self.A == 0.0 and self.i12 == 0.0:
                        assert self.K1 is None, 'K1 must be blank if A=0.0 and I12=0.0; A=%r I12=%r K1=%r' % (self.A, self.i12, self.K1)
                        assert self.K2 is None, 'K2 must be blank if A=0.0 and I12=0.0; A=%r I12=%r K2=%r' % (self.A, self.i12, self.K2)
                    assert len(card) <= 20, 'len(PBAR card) = %i' % len(card)

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
            self._cards = []
            self._comments = []
        else:
            self.property_id = array([], dtype='int32')

    #=========================================================================
    def get_index(self, property_ids):
        if isinstance(property_ids, int):
            property_ids = array([property_ids])
        if property_ids is None:
            return arange(self.n)

        indexs = searchsorted(self.property_id, property_ids)
        assert len(indexs) == len(property_ids), 'indexs=%s pids=%s' % (indexs, property_ids)
        return indexs

    #=========================================================================
    def write_bdf(self, f, size=8, property_ids=None):
        if self.n:
            if property_ids is None:
                i = arange(self.n)
            else:
                i = searchsorted(self.property_id, property_ids)

            for (pid, mid, area, I1, I2, J) in zip(self.property_id[i], self.material_id[i],
                    self.area[i], self.I1[i], self.I2[i], self.J[i]):
                card = ['PBAR', pid, mid, area, I1, I2, J]
                f.write(print_card_8(card))
