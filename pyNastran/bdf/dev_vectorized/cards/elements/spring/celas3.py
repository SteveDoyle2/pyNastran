from itertools import izip

from numpy import arange, array, dot, zeros, unique, searchsorted, transpose
from numpy.linalg import norm

from ..rod.conrod import _Lambda

from pyNastran.bdf.fieldWriter import print_card
from pyNastran.bdf.bdfInterface.assign_type import (integer, integer_or_blank,
    double, double_or_blank, integer_double_or_blank, blank)

class CELAS3(object):
    type = 'CELAS3'
    def __init__(self, model):
        """
        Defines the CELAS3 object.

        :param self: the CELAS3 object
        :param model: the BDF object
        """
        self.model = model
        self.n = 0
        self._cards = []
        self._comments = []

        self.element_id = []
        #: Property ID
        self.property_id = []

        # Node IDs
        self.node_ids = []

        #: stiffness of the scalar spring
        self.K = []

        #: component number
        self.components = []

        #: damping coefficient
        self.ge = []

        #: stress coefficient
        self.s = []


    def add(self, card, comment=None):
        self._cards.append(card)
        self._comments.append(comment)

        self.element_id.append(integer(card, 1, 'eid'))
        self.K.append(double(card, 2, 'k'))
        self.node_ids.append([integer(card, 3, 'G1'),
                              integer_or_blank(card, 5, 'G2')])
        self.components.append([integer_or_blank(card, 4, 'C1', 0),
                              integer_or_blank(card, 6, 'C2', 0)])
        self.ge.append(double_or_blank(card, 7, 'ge', 0.))
        self.s.append(double_or_blank(card, 8, 's', 0.))
        assert len(card) <= 9, 'len(CELAS3 card) = %i' % len(card) + str(card)

    def build(self):
        """
        :param self: the CELAS3 object
        """
        cards = self._cards
        ncards = len(cards)
        self.n = ncards
        
        if ncards:
            float_fmt = self.model.float
            self.element_id = array(self.element_id, 'int32')
            self.node_ids = array(self.node_ids)
            self.K = array(self.K, float_fmt)
            self.components = array(self.components, 'int32')
            self.ge = array(self.ge, float_fmt)
            self.s = array(self.s, float_fmt)
            return

            float_fmt = self.model.float

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

            for i, card in enumerate(cards):
                self.element_id[i] = integer(card, 1, 'eid')
                self.K[i] = double(card, 2, 'k')
                self.node_ids[i] = [integer(card, 3, 'G1'),
                                    integer(card, 5, 'G2')]
                self.components[i] = [integer_or_blank(card, 4, 'C1', 0),
                                      integer_or_blank(card, 6, 'C2', 0)]
                self.ge[i] = double_or_blank(card, 7, 'ge', 0.)
                self.s[i] = double_or_blank(card, 8, 's', 0.)
                assert len(card) <= 9, 'len(CELAS3 card) = %i' % len(card) + str(card)

            i = self.element_id.argsort()
            self.element_id = self.element_id[i]
            self.K = self.K[i]
            self.node_ids = self.node_ids[i, :]
            self.components = self.components[i, :]
            self.ge = self.ge[i]
            self.s = self.s[i]

            unique_eids = unique(self.element_id)
            if len(unique_eids) != len(self.element_id):
                raise RuntimeError('There are duplicate CELAS3 IDs...')
            self._cards = []
            self._comments = []

    def get_stats(self):
        msg = []
        if self.n:
            msg.append('  %-8s: %i' % ('CELAS3', self.n))
        return msg

    def write_bdf(self, f, size=8, eids=None):
        if self.n:
            if eids is None:
                i = arange(self.n)
            else:
                i = searchsorted(self.element_id, self.eid)

            N0 = self.node_ids[i, 0]
            N1 = self.node_ids[i, 1]
            C0 = self.components[i, 0]
            C1 = self.components[i, 1]
            for (eid, k, n0, n1, c0, c1, ge, s) in izip(self.element_id[i],
                    self.K[i], N0, N1, C0, C1, self.ge[i], self.s[i]):
                card = ['CELAS3', eid, k, n0, c0, n1, c1, ge, s ]
                f.write(print_card(card))

    def get_stiffness(self, i, model, positions, index0s, fnorm=1.0):  # CELAS3
        ki = self.K[i]
        
        k = ki * array([[1, -1,],
                        [-1, 1]])

        n0, n1 = self.node_ids[i, :]

        p0 = positions[n0]
        p1 = positions[n1]
        
        v1 = p0 - p1
        L = norm(v1)
        if L == 0.0:
            msg = 'invalid CELAS3 length=0.0\n%s' % (self.__repr__())
            raise ZeroDivisionError(msg)

        try:
            Lambda = _Lambda(v1, debug=True)
        except ZeroDivisionError:
            raise ZeroDivisionError("CELAS3 xyz[%i]=%s; xyz[%i]=%s" % (n0, p0, n1, p1))

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
        for i in xrange(self.n):
            n0, n1 = self.node_ids[i, :]
            if n0 == n1:
                raise RuntimeError('CELAS3 eid=%s n1=%s n2=%s' % (self.element_id[i], n0, n1))
            p0 = positions[n0]
            p1 = positions[n1]

            v1 = p0 - p1
            L = norm(v1)
            try:
                Lambda = _Lambda(v1, debug=True)
            except ZeroDivisionError:
                raise ZeroDivisionError("CELAS3 xyz[%i]=%s; xyz[%i]=%s" % (n0, p0, n1, p1))

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
