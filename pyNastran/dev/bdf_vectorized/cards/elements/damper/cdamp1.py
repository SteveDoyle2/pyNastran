from numpy import array, zeros, unique, searchsorted, where, arange

from pyNastran.dev.bdf_vectorized.cards.elements.damper.damp_element import DamperElement

from pyNastran.bdf.field_writer_8 import print_card_8
from pyNastran.bdf.field_writer_16 import print_card_16
from pyNastran.bdf.bdf_interface.assign_type import integer, integer_or_blank
from pyNastran.bdf.bdf_interface.bdf_card import BDFCard


class CDAMP1(DamperElement):
    type = 'CDAMP1'
    def __init__(self, model):
        """
        Defines the CDAMP1 object.

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
            #: Element ID
            self.element_id = zeros(ncards, 'int32')
            #: Property ID
            self.property_id = zeros(ncards, 'int32')
            #: Node IDs
            self.node_ids = zeros((ncards, 2), 'int32')
            #: component number
            self.components = zeros((ncards, 2), 'int32')

    def add_card(self, card: BDFCard, comment: str=''):
        i = self.i
        self.element_id[i] = integer(card, 1, 'eid')
        self.property_id[i] = integer_or_blank(card, 2, 'pid', self.element_id[i])
        self.node_ids[i, :] = [integer(card, 3, 'n1'),
                               integer(card, 5, 'n2')]
        self.components[i, :] = [integer_or_blank(card, 4, 'c1', 0),
                                 integer_or_blank(card, 6, 'c2', 0)]
        assert len(card) <= 7, 'len(CDAMP1 card) = %i\ncard=%s' % (len(card), card)
        self.i += 1

    def build(self):
        if self.n:
            i = self.element_id.argsort()
            self.element_id = self.element_id[i]
            self.property_id = self.property_id[i]
            self.node_ids = self.node_ids[i, :]
            self.components = self.components[i, :]

            unique_eids = unique(self.element_id)
            if len(unique_eids) != len(self.element_id):
                raise RuntimeError('There are duplicate CDAMP1 IDs...')
            self._cards = []
        else:
            self.element_id = array([], dtype='int32')
            self.property_id = array([], dtype='int32')

    def update(self, maps):
        """
        maps = {
            'node_id' : nid_map,
            'property' : pid_map,
        }
        """
        if self.n:
            eid_map = maps['element']
            pid_map = maps['property']
            nid_map = maps['node']
            for i, eid, pid, nids in enumerate(zip(self.element_id, self.property_id, self.node_ids)):
                self.element_id[i] = eid_map[eid]
                self.property_id[i] = pid_map[pid]
                self.node_ids[i, 0] = nid_map[nids[0]]
                self.node_ids[i, 1] = nid_map[nids[1]]

    def write_card(self, bdf_file, size=8, eids=None):
        if self.n:
            if eids is None:
                i = arange(self.n)
            else:
                i = searchsorted(self.element_id, self.eid)

            for (eid, pid, n, c) in zip(self.element_id[i], self.property_id[i], self.node_ids[i], self.components[i]):
                card = ['CDAMP1', eid, pid, n[0], n[1], c[0], c[1]]
                if size == 8:
                    bdf_file.write(print_card_8(card))
                else:
                    bdf_file.write(print_card_16(card))

    def get_damping_matrix(self, i, model, positions, index0s, fnorm=1.0):
        """gets the damping matrix for CDAMP1"""
        #print("----------------")
        ipid = where(self.model.pdamp.property_id==self.property_id[i])[0][0]
        prop = self.model.pdamp
        bi = prop.B[ipid]
        k = bi * array([[1, -1,],
                        [-1, 1]])

        #========================
        n1, n2 = self.node_ids[i, :]
        c1, c2 = self.components[i, :]
        #i0, i1 = index0s

        delta1 = 0 if c1 in [0, 1, 2, 3] else 3
        delta2 = 0 if c2 in [0, 1, 2, 3] else 3

        c1b = c1-1 if c1 > 0 else c1
        c2b = c2-1 if c2 > 0 else c2

        i1 = index0s[n1]
        i2 = index0s[n2]
        dofs = [
            i1 + c1b,
            i2 + c1b,
        ]

        n_ijv = [
            (n1, 1 + delta1),
            (n2, 1 + delta2),
        ]
        return k, dofs, n_ijv

    def displacement_stress(self, model, positions, q, dofs,
                            ni, o1, e1, f1):
        n = self.n

        du_axial = zeros(n, 'float64')
        for i in range(self.n):
            (n1, n2) = self.node_ids[i, :]

            n11 = dofs[(n1, 1)]
            n21 = dofs[(n2, 1)]

            q_axial = array([
                q[n11],
                q[n21],
            ])
            u_axial = q_axial
            du_axial[i] = u_axial[0] - u_axial[1]

        self.model.log.debug("len(pdamp) = %s" % self.model.pdamp.n)
        i = searchsorted(self.model.pdamp.property_id, self.property_id)
        k = self.model.pdamp.B[i]
        s = self.model.pdamp.s[i]
        self.model.log.debug("k=%s s=%s du_axial=%s" % (k, s, du_axial))

        e1[ni: ni+n] = du_axial * s
        f1[ni: ni+n] = k * du_axial
        o1[ni: ni+n] = f1[ni: ni+n] * s
        #return (axial_strain, axial_stress, axial_force)
