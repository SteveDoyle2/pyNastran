from numpy import dot, array, zeros, unique, searchsorted, transpose, where
from numpy.linalg import norm

from .conrod import _Lambda

from pyNastran.bdf.fieldWriter import print_card
from pyNastran.bdf.bdfInterface.assign_type import (integer, integer_or_blank,
    double_or_blank, integer_double_or_blank, blank)

class CELAS1(object):
    type = 'CELAS1'
    def __init__(self, model):
        """
        Defines the CELAS1 object.

        :param self: the CELAS1 object
        :param model: the BDF object
        """
        self.model = model
        self.n = 0
        self._cards = []
        self._comments = []

    def add(self, card, comment=None):
        self._cards.append(card)
        self._comments.append(comment)

    def build(self):
        """
        :param self: the CELAS1 object
        """
        cards = self._cards
        ncards = len(cards)
        self.n = ncards
        if ncards:
            #: Element ID
            self.element_id = zeros(ncards, 'int32')
            #: Property ID
            self.property_id = zeros(ncards, 'int32')
            #: Node IDs
            self.node_ids = zeros((ncards, 2), 'int32')
            #: component number
            self.components = zeros((ncards, 2), 'int32')

            for i, card in enumerate(cards):
                self.element_id[i] = integer(card, 1, 'eid')
                self.property_id[i] = integer_or_blank(card, 2, 'pid', self.element_id[i])
                self.node_ids[i] = [integer(card, 3, 'n1'),
                                    integer(card, 5, 'n2')]
                self.components[i] = [integer_or_blank(card, 4, 'c1', 0),
                                      integer_or_blank(card, 6, 'c2', 0)]
                assert len(card) <= 7, 'len(CELAS1 card) = %i' % len(card)

            i = self.element_id.argsort()
            self.element_id = self.element_id[i]
            self.property_id = self.property_id[i]
            self.node_ids = self.node_ids[i, :]
            self.components = self.components[i, :]

            unique_eids = unique(self.element_id)
            if len(unique_eids) != len(self.element_id):
                raise RuntimeError('There are duplicate CELAS1 IDs...')
            self._cards = []
            self._comments = []

    def get_stats(self):
        msg = []
        if self.n:
            msg.append('  %-8s: %i' % ('CELAS1', self.n))
        return msg

    def write_bdf(self, f, size=8, eids=None):
        if self.n:
            if eids is None:
                i = arange(self.n)
            else:
                i = searchsorted(self.element_id, self.eid)

            for (eid, pid, n, c) in zip(self.element_id[i], self.property_id[i], self.node_ids[i], self.components[i]):
                card = ['CELAS1', eid, pid, n[0], n[1], c[0], c[1] ]
                f.write(print_card(card))

    def get_stiffness(self, i, model, positions, index0s, fnorm=1.0):  # CELAS1/CONROD
        #print("----------------")
        ipid = where(self.model.pelas.property_id==self.property_id[i])[0][0]
        prop = self.model.pelas
        ki = prop.K[ipid]

        #========================
        n0, n1 = self.node_ids[i, :]

        i0, i1 = index0s
        p0 = positions[n0]
        p1 = positions[n1]

        v1 = p0 - p1
        L = norm(v1)
        if L == 0.0:
            msg = 'invalid CELAS1 length=0.0\n%s' % (self.__repr__())
            raise ZeroDivisionError(msg)

        try:
            Lambda = _Lambda(v1, debug=True)
        except ZeroDivisionError:
            raise ZeroDivisionError("CELAS2 xyz[%i]=%s; xyz[%i]=%s" % (n0, p0, n1, p1))

        #========================
        k = ki * array([[1, -1,],
                        [-1, 1]])
        K = dot(dot(transpose(Lambda), k), Lambda)

        c1, c2 = self.components[i, :]
        delta1 = 0 if c1 in [1, 2, 3] else 3
        delta2 = 0 if c2 in [1, 2, 3] else 3

        nIJV = [
            (n0, 1 + delta1), (n0, 2 + delta1), (n0, 3 + delta1),
            (n1, 1 + delta2), (n1, 2 + delta2), (n1, 3 + delta2),
        ]
        dofs = nIJV
        return (K, dofs, nIJV)

    def displacement_stress(self, model, positions, q, dofs):
        n = self.n
        o1 = zeros(n, 'float64')
        e1 = zeros(n, 'float64')
        f1 = zeros(n, 'float64')


        du_axial = zeros(n, 'float64')
        for i in xrange(self.n):
            (n0, n1) = self.node_ids[i, :]
            if n0 == n1:
                raise RuntimeError('CELAS2 eid=%s n1=%s n2=%s' % (self.element_id[i], n0, n1))
            p0 = positions[n0]
            p1 = positions[n1]

            v1 = p0 - p1
            L = norm(v1)
            try:
                Lambda = _Lambda(v1, debug=True)
            except ZeroDivisionError:
                raise ZeroDivisionError("CELAS2 xyz[%i]=%s; xyz[%i]=%s" % (n0, p0, n1, p1))

            n01 = dofs[(n0, 1)]
            n11 = dofs[(n1, 1)]

            n02 = dofs[(n0, 2)]
            n12 = dofs[(n1, 2)]

            n03 = dofs[(n0, 3)]
            n13 = dofs[(n1, 3)]

            q_axial = array([
                q[n01], q[n02], q[n03],
                q[n11], q[n12], q[n13]
            ])
            u_axial = dot(Lambda, q_axial)
            du_axial[i] = u_axial[0] - u_axial[1]

        print "len(pelas) =", self.model.pelas.n
        i = searchsorted(self.model.pelas.property_id, self.property_id)
        k = self.model.pelas.K[i]
        s = self.model.pelas.s[i]
        print "k=%s s=%s du_axial=%s" % (k, s, du_axial)

        axial_strain = du_axial * s
        axial_force = k * du_axial
        axial_stress = axial_force * s

        return (axial_strain, axial_stress, axial_force)