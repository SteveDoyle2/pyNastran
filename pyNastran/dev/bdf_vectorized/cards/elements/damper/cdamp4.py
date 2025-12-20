from numpy import arange, array, dot, zeros, unique, searchsorted, transpose
from numpy.linalg import norm  # type: ignore

from pyNastran.dev.bdf_vectorized.cards.elements.damper.damp_element import DamperElement

from pyNastran.bdf.field_writer_8 import print_card_8
from pyNastran.bdf.field_writer_16 import print_card_16
from pyNastran.bdf.bdf_interface.assign_type import (integer, integer_or_blank,
    double, double_or_blank)
from pyNastran.bdf.bdf_interface.bdf_card import BDFCard

class CDAMP4(DamperElement):
    type = 'CDAMP4'
    def __init__(self, model):
        """
        Defines the CDAMP4 object.

        Parameters
        ----------
        model : BDF
           the BDF object
        """
        DamperElement.__init__(self, model)

    def allocate(self, card_count):
        ncards = card_count[self.type]
        if ncards:
            self.n = ncards
            float_fmt = self.model.float_fmt
            #: Element ID
            self.element_id = zeros(ncards, 'int32')
            #: Property ID
            self.property_id = zeros(ncards, 'int32')
            # Node IDs
            self.node_ids = zeros((ncards, 2), 'int32')
            #: damping of the scalar damper
            self.B = zeros(ncards, float_fmt)
            #: component number
            self.components = zeros((ncards, 2), 'int32')
            #: damping coefficient
            self.ge = zeros(ncards, float_fmt)
            #: stress coefficient
            self.s = zeros(ncards, float_fmt)


    def add_card(self, card, comment=None):
        i = self.i
        self.element_id[i] = integer(card, 1, 'eid')
        self.B[i] = double(card, 2, 'k')
        self.node_ids[i, :] = [integer(card, 3, 'G1'),
                               integer(card, 5, 'G2')]
        self.components[i, :] = [integer_or_blank(card, 4, 'C1', 0),
                                 integer_or_blank(card, 6, 'C2', 0)]
        self.ge[i] = double_or_blank(card, 7, 'ge', 0.)
        self.s[i] = double_or_blank(card, 8, 's', 0.)
        assert len(card) <= 9, 'len(CDAMP4 card) = %i\ncard=%s' % (len(card), card) + str(card)
        self.i += 1

    def build(self):
        if self.n:
            i = self.element_id.argsort()
            self.element_id = self.element_id[i]
            self.B = self.B[i]
            self.node_ids = self.node_ids[i, :]
            self.components = self.components[i, :]
            self.ge = self.ge[i]
            self.s = self.s[i]

            unique_eids = unique(self.element_id)
            if len(unique_eids) != len(self.element_id):
                raise RuntimeError('There are duplicate CDAMP4 IDs...')
            self._cards = []
        else:
            self.element_id = array([], dtype='int32')
            self.property_id = array([], dtype='int32')

    def update(self, maps):
        """
        maps = {
            'node_id' : nid_map,
        }
        """
        if self.n:
            eid_map = maps['element']
            nid_map = maps['node']
            for i, eid, nids in enumerate(zip(self.element_id, self.node_ids)):
                self.element_id[i] = eid_map[eid]
                self.node_ids[i, 0] = nid_map[nids[0]]
                self.node_ids[i, 1] = nid_map[nids[1]]

    def write_card(self, bdf_file, size=8, eids=None):
        if self.n:
            if eids is None:
                i = arange(self.n)
            else:
                i = searchsorted(self.element_id, self.eid)

            N0 = self.node_ids[i, 0]
            N1 = self.node_ids[i, 1]
            C0 = self.components[i, 0]
            C1 = self.components[i, 1]
            for (eid, k, n0, n1, c0, c1, ge, s) in zip(self.element_id[i],
                    self.B[i], N0, N1, C0, C1, self.ge[i], self.s[i]):
                card = ['CDAMP4', eid, k, n0, c0, n1, c1, ge, s]
                if size == 8:
                    bdf_file.write(print_card_8(card))
                else:
                    bdf_file.write(print_card_16(card))

    def get_damping_matrix(self, i, model, positions, index0s, fnorm=1.0):
        """gets the damping matrix for CDAMP4"""
        bi = self.B[i]
        k = bi * array([[1, -1,],
                        [-1, 1]])

        n1, n2 = self.node_ids[i, :]

        n_ijv = [
            (n1, 1),
            (n2, 1),
        ]
        dofs = n_ijv
        return k, dofs, n_ijv

    def displacement_stress(self, model, positions, q, dofs,
                            ni, o1, e1, f1):

        n = self.n
        du_axial = zeros(n, 'float64')
        for i in range(self.n):
            n0, n1 = self.node_ids[i, :]

            n11 = dofs[(n0, 1)]
            n21 = dofs[(n1, 1)]

            q_axial = array([
                q[n11],
                q[n21],
            ])
            u_axial = q_axial
            du_axial[i] = u_axial[0] - u_axial[1]

        s = self.s
        bi = self.B

        e1[ni : ni+n] = du_axial * s
        f1[ni : ni+n] = bi * du_axial
        o1[ni : ni+n] = f1[ni: ni+n] * s
        #return (axial_strain, axial_stress, axial_force)

    def slice_by_index(self, i):
        obj = CDAMP4(self.model)
        obj.n = len(i)
        #obj._cards = self._cards[i]
        #obj._comments = obj._comments[i]
        #obj.comments = obj.comments[i]
        obj.element_id = self.element_id[i]
        obj.property_id = self.property_id[i]
        obj.node_ids = self.node_ids[i, :]
        return obj
