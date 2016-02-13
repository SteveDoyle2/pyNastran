from six.moves import zip, range

from numpy import arange, array, dot, zeros, unique, searchsorted, transpose
from numpy.linalg import norm

from pyNastran.bdf.dev_vectorized.cards.elements.rod.conrod import _Lambda
from pyNastran.bdf.dev_vectorized.cards.elements.spring.spring_element import SpringElement

from pyNastran.bdf.field_writer_8 import print_card_8
from pyNastran.bdf.field_writer_16 import print_card_16
from pyNastran.bdf.bdfInterface.assign_type import (integer, integer_or_blank,
    double, double_or_blank)

class CELAS4(SpringElement):
    type = 'CELAS4'
    def __init__(self, model):
        """
        Defines the CELAS4 object.

        :param model: the BDF object
        """
        SpringElement.__init__(self, model)

    def allocate(self, ncards):
        self.n = ncards
        float_fmt = self.model.float
        #: Element ID
        self.element_id = zeros(ncards, 'int32')
        #: Property ID
        self.property_id = zeros(ncards, 'int32')
        # Node IDs
        self.node_ids = zeros((ncards, 2), 'int32')
        #: stiffness of the scalar spring
        self.K = zeros(ncards, float_fmt)
        #: component number
        self.components = zeros((ncards, 2), 'int32')
        #: damping coefficient
        self.ge = zeros(ncards, float_fmt)
        #: stress coefficient
        self.s = zeros(ncards, float_fmt)


    def add(self, card, comment=None):
        i = self.i
        self.element_id[i] = integer(card, 1, 'eid')
        self.K[i] = double(card, 2, 'k')
        self.node_ids[i, :] = [integer(card, 3, 'G1'),
                               integer(card, 5, 'G2')]
        self.components[i, :] = [integer_or_blank(card, 4, 'C1', 0),
                                 integer_or_blank(card, 6, 'C2', 0)]
        self.ge[i] = double_or_blank(card, 7, 'ge', 0.)
        self.s[i] = double_or_blank(card, 8, 's', 0.)
        assert len(card) <= 9, 'len(CELAS4 card) = %i' % len(card) + str(card)
        self.i += 1

    def build(self):
        if self.n:
            i = self.element_id.argsort()
            self.element_id = self.element_id[i]
            self.K = self.K[i]
            self.node_ids = self.node_ids[i, :]
            self.components = self.components[i, :]
            self.ge = self.ge[i]
            self.s = self.s[i]

            unique_eids = unique(self.element_id)
            if len(unique_eids) != len(self.element_id):
                raise RuntimeError('There are duplicate CELAS4 IDs...')
            self._cards = []
            self._comments = []
        else:
            self.element_id = array([], dtype='int32')
            self.property_id = array([], dtype='int32')

    def write_card(self, f, size=8, eids=None):
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
                    self.K[i], N0, N1, C0, C1, self.ge[i], self.s[i]):
                card = ['CELAS4', eid, k, n0, c0, n1, c1, ge, s]
                if size == 8:
                    f.write(print_card_8(card))
                else:
                    f.write(print_card_16(card))

    def get_stiffness(self, i, model, positions, index0s, fnorm=1.0):  # CELAS4
        ki = self.K[i]

        k = ki * array([[1, -1,],
                        [-1, 1]])

        n0, n1 = self.node_ids[i, :]

        p0 = positions[n0]
        p1 = positions[n1]

        v1 = p0 - p1
        L = norm(v1)
        if L == 0.0:
            msg = 'invalid CELAS4 length=0.0\n%s' % (self.__repr__())
            raise ZeroDivisionError(msg)

        try:
            Lambda = _Lambda(v1, debug=True)
        except ZeroDivisionError:
            raise ZeroDivisionError("CELAS4 xyz[%i]=%s; xyz[%i]=%s" % (n0, p0, n1, p1))

        #========================
        K = dot(dot(transpose(Lambda), k), Lambda)

        c0, c1 = self.components[i, :]
        n0, n1 = self.node_ids[i, :]
        delta0 = 0 if c0 in [1, 2, 3] else 3
        delta1 = 0 if c1 in [1, 2, 3] else 3

        nIJV = [
            (n0, 1 + delta0), (n0, 2 + delta0), (n0, 3 + delta0),
            (n1, 1 + delta1), (n1, 2 + delta1), (n1, 3 + delta1),
        ]
        dofs = nIJV
        return (K, dofs, nIJV)

    def displacement_stress(self, model, positions, q, dofs,
                            ni, e1, f1, o1):

        n = self.n
        du_axial = zeros(n, 'float64')
        for i in range(self.n):
            n0, n1 = self.node_ids[i, :]
            if n0 == n1:
                raise RuntimeError('CELAS4 eid=%s n1=%s n2=%s' % (self.element_id[i], n0, n1))
            p0 = positions[n0]
            p1 = positions[n1]

            v1 = p0 - p1
            L = norm(v1)
            try:
                Lambda = _Lambda(v1, debug=True)
            except ZeroDivisionError:
                raise ZeroDivisionError("CELAS4 xyz[%i]=%s; xyz[%i]=%s" % (n0, p0, n1, p1))

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

        s = self.s
        ki = self.K

        e1[ni : ni+n] = du_axial * s
        f1[ni : ni+n] = ki * du_axial
        o1[ni : ni+n] = f1[ni: ni+n] * s

        #return (axial_strain, axial_stress, axial_force)

    def slice_by_index(self, i):
        obj = CELAS4(self.model)
        obj.n = len(i)
        #obj._cards = self._cards[i]
        #obj._comments = obj._comments[i]
        #obj.comments = obj.comments[i]
        obj.element_id = self.element_id[i]
        obj.property_id = self.property_id[i]
        obj.node_ids = self.node_ids[i, :]
        return obj
