from numpy import array, arange, zeros, unique, searchsorted

from pyNastran.bdf.field_writer_8 import print_float_8
from pyNastran.bdf.field_writer_16 import print_float_16

from pyNastran.bdf.bdf_interface.assign_type import (integer, integer_or_blank,
    double_or_blank)

from pyNastran.dev.bdf_vectorized.cards.vectorized_card import VectorizedCard

class CONM2(VectorizedCard):
    type = 'CONM2'
    def __init__(self, model):
        """
        Defines the CONM2 object.

        Parameters
        ----------
        model : BDF
           the BDF object
        """
        VectorizedCard.__init__(self, model)

    def allocate(self, card_count):
        ncards = card_count[self.type]
        if ncards:
            self.n = ncards
            float_fmt = self.model.float_fmt

            #: Element ID
            self.element_id = zeros(ncards, 'int32')
            #: Property ID
            self.property_id = zeros(ncards, 'int32')
            self.node_id = zeros(ncards, 'int32')
            self.coord_id = zeros(ncards, 'int32')
            self.mass = zeros(ncards, float_fmt)
            self.x = zeros((ncards, 3), float_fmt)
            self.I = zeros((ncards, 6), float_fmt)

    def add_card(self, card, comment=''):
        eid = integer(card, 1, 'element_id')
        if comment:
            self.set_comment(eid, comment)
        i = self.i
        self.element_id[i] = eid
        self.node_id[i] = integer(card, 2, 'node_id')
        self.coord_id[i] = integer_or_blank(card, 3, 'coord_id', 0)
        self.mass[i] = double_or_blank(card, 4, 'mass', 0.)
        self.x[i, :] = [double_or_blank(card, 5, 'x1', 0.0),
                        double_or_blank(card, 6, 'x2', 0.0),
                        double_or_blank(card, 7, 'x3', 0.0)]

        self.I[i, :] = [double_or_blank(card, 9, 'I11', 0.0),
                        double_or_blank(card, 10, 'I21', 0.0),
                        double_or_blank(card, 11, 'I22', 0.0),
                        double_or_blank(card, 12, 'I31', 0.0),
                        double_or_blank(card, 13, 'I32', 0.0),
                        double_or_blank(card, 14, 'I33', 0.0)]
        assert len(card) <= 15, 'len(CONM2 card) = %i\ncard=%s' % (len(card), card)
        self.i += 1

    def build(self):
        if self.n:
            i = self.element_id.argsort()
            #print("iconm2 %s %s" % (i, type(i)))
            self.element_id = self.element_id[i]
            self.node_id = self.node_id[i]
            self.coord_id = self.coord_id[i]
            self.mass = self.mass[i]
            self.x = self.x[i, :]
            self.I = self.I[i, :]

            unique_eids = unique(self.element_id)
            if len(unique_eids) != len(self.element_id):
                raise RuntimeError('There are duplicate CONM2 IDs...')
            #self._cards = []
        else:
            self.element_id = array([], dtype='int32')
            self.property_id = array([], dtype='int32')

    def get_mass_by_element_id(self, element_id=None, total=False):
        """
        mass = rho * A * L + nsm
        """
        if element_id is None:
            element_id = arange(self.n)
        #grid_cid0 = self.model.grid.get_position_by_index()
        #p = grid_cid0[self.node_id]

        mass = self.mass
        if total:
            return mass.sum()
        else:
            return mass

    def get_stats(self):
        msg = []
        if self.n:
            msg.append('  %-8s: %i' % ('CONM2', self.n))
        return msg

    def write_card(self, bdf_file, size=8, is_double=False, element_id=None):
        assert self.n > 0, self.n
        if self.n:
            if element_id is None:
                i = arange(self.n)
            else:
                i = searchsorted(self.element_id, element_id)

            blank = ' ' * size
            #cid = [cid if cid != 0 else '' for cid in self.coord_id]


            #Mass = [x if x != 0.0 else '' for x in self.mass[i]]
            if size == 8:
                X0 = [print_float_8(x) if x != 0.0 else blank for x in self.x[i, 0]]
                X1 = [print_float_8(x) if x != 0.0 else blank for x in self.x[i, 1]]
                X2 = [print_float_8(x) if x != 0.0 else '' for x in self.x[i, 2]]

                I0 = [print_float_8(x) if x != 0.0 else blank for x in self.I[i, 0]]
                I1 = [print_float_8(x) if x != 0.0 else blank for x in self.I[i, 1]]
                I2 = [print_float_8(x) if x != 0.0 else blank for x in self.I[i, 2]]
                I3 = [print_float_8(x) if x != 0.0 else blank for x in self.I[i, 3]]
                I4 = [print_float_8(x) if x != 0.0 else blank for x in self.I[i, 4]]
                I5 = [print_float_8(x) if x != 0.0 else blank for x in self.I[i, 5]]
                for (eid, nid, cid, mass, x0, x1, x2, i0, i1, i2, i3, i4, i5) in zip(
                    self.element_id[i], self.node_id[i], self.coord_id[i], self.mass[i],
                    X0, X1, X2, I0, I1, I2, I3, I4, I5):
                    if eid in self._comments:
                        bdf_file.write(self._comments[eid])
                    carda = '%-8s%8i%8i%8i%8s%s%s%s\n' % ('CONM2', eid, nid, cid, print_float_8(mass), x0, x1, x2)
                    cardb = '%8s%8s%8s%8s%8s%8s%8s' % ('', i0, i1, i2, i3, i4, i5)
                    cardi = (carda + cardb).rstrip() + '\n'
                    bdf_file.write(cardi)
            else:
                assert size == 16, size
                X0 = [print_float_16(x) if x != 0.0 else blank for x in self.x[i, 0]]
                X1 = [print_float_16(x) if x != 0.0 else blank for x in self.x[i, 1]]
                X2 = [print_float_16(x) if x != 0.0 else '' for x in self.x[i, 2]]

                I0 = [print_float_16(x) if x != 0.0 else blank for x in self.I[i, 0]]
                I1 = [print_float_16(x) if x != 0.0 else blank for x in self.I[i, 1]]
                I2 = [print_float_16(x) if x != 0.0 else blank for x in self.I[i, 2]]
                I3 = [print_float_16(x) if x != 0.0 else blank for x in self.I[i, 3]]
                I4 = [print_float_16(x) if x != 0.0 else blank for x in self.I[i, 4]]
                I5 = [print_float_16(x) if x != 0.0 else blank for x in self.I[i, 5]]
                for (eid, nid, cid, mass, x0, x1, x2, i0, i1, i2, i3, i4, i5) in zip(
                    self.element_id[i], self.node_id[i], self.coord_id[i], self.mass[i],
                    X0, X1, X2, I0, I1, I2, I3, I4, I5):
                    if eid in self._comments:
                        bdf_file.write(self._comments[eid])
                    carda = '%-8s%16i%16i%16i%16s\n' % ('CONM2*', eid, nid, cid, print_float_8(mass), )
                    cardb = '%-8s%16s%16s%16s%16s\n' % ('*', x0, x1, x2, i0)
                    cardc = '%-8s%16s%16s%16s%16s\n' % ('*', i1, i2, i3, i4)
                    cardd = '%-8s%16s' % ('*', i5)
                    cardcd = (cardc + cardd).rstrip('* \n')
                    if cardcd:
                        cardi = carda + cardb + cardc + cardd
                    else:
                        cardi = carda + cardb
                    bdf_file.write(cardi.rstrip() + '\n')

    #def get_stiffness_matrix(self, model, node_ids, index0s, fnorm=1.0):
        #return(K, dofs, n_ijv)
