# pylint: disable=R0904,R0902,E1101,E1103,C0111,C0302,C0103,W0101
from __future__ import (nested_scopes, generators, division, absolute_import,
                        print_function, unicode_literals)
from six.moves import range

from numpy import matrix, zeros, array, transpose, dot, ones
from numpy import eye, allclose, cross
from numpy.linalg import norm

from pyNastran.utils.dev import list_print
from pyNastran.bdf.fieldWriter import set_blank_if_default
from pyNastran.bdf.cards.baseCard import Element #, Mid
from pyNastran.bdf.bdfInterface.assign_type import (integer, integer_or_blank,
    integer_double_or_blank, double, double_or_blank,
    string_or_blank, integer_double_string_or_blank, integer_or_double)
from pyNastran.bdf.fieldWriter import print_card_8
from pyNastran.bdf.fieldWriter16 import print_card_16


def _Lambda(model, n1, n2, debug=True):
    """
    ::
      3d  [l,m,n,0,0,0]  2x6
          [0,0,0,l,m,n]
    """
    #R = self.Rmatrix(model,is3D)

    p1 = model.Node(n1).Position()
    p2 = model.Node(n2).Position()
    v1 = p2 - p1
    #if debug:
        #print("v1=%s" % (v1))
    v1 = v1 / norm(v1)
    (l, m, n) = v1
    #l = 1
    #m = 2
    #n = 3
    Lambda = matrix(zeros((2, 6), 'd'))
    Lambda[0, 0] = Lambda[1, 3] = l
    Lambda[0, 1] = Lambda[1, 4] = m
    Lambda[0, 2] = Lambda[1, 5] = n

    #print("R = \n",R)
    #debug = True
    #if debug:
        #print("Lambda = \n" + str(Lambda))
        #sys.exit('asdf')
    return Lambda


class RodElement(Element):  # CROD, CONROD, CTUBE

    def __init__(self, card, data):
        Element.__init__(self, card, data)

    def cross_reference(self, model):
        msg = ' which is required by %s eid=%s' % (self.type, self.eid)
        self.nodes = model.Nodes(self.nodes, msg=msg)
        self.pid = model.Property(self.pid, msg=msg)

    def nodeIDs(self):
        return self._nodeIDs(allowEmptyNodes=False)

    def Rho(self):
        r"""returns the material density  \f$ \rho \f$"""
        return self.pid.mid.rho

    def Length(self):
        r"""
        Gets the length of the element.

        .. math:: L = \sqrt{  (n_{x2}-n_{x1})^2+(n_{y2}-n_{y1})^2+(n_{z2}-n_{z1})^2  }
        :param self: the CROD/CONROD/CTUBE element
        """
        L = norm(self.nodes[1].Position() - self.nodes[0].Position())
        return L

    def Mass(self):
        r"""
        get the mass of the element.

        .. math:: m = \left( \rho A + nsm \right) L
        """
        L = self.Length()
        mass = (self.Rho() * self.Area() + self.Nsm()) * L
        return mass

    def Rmatrix(self, model, is3D):
        r"""
        where   :math:`[R]_{ij}` is the tranformation matrix

        .. math::
          [R]_{ij} = \left[
          \begin{array}{ccc}
              g_x \dot e_x & g_x \dot e_y &  g_x \dot e_z    \\
              g_y \dot e_x & g_y \dot e_y &  g_y \dot e_z    \\
              g_z \dot e_x & g_z \dot e_y &  g_z \dot e_z
          \end{array} \right]
        """
        (n1, n2) = self.nodeIDs()
        v1 = model.Node(n2).xyz - model.Node(n1).xyz
        v1 /= norm(v1)

        v1x = array([v1[0], 0., 0.])
        v1y = array([0., v1[1], 0.])
        v1z = array([0., 0., v1[2]])

        g1x = array([1., 0., 0.])
        g1y = array([0., 1., 0.])
        g1z = array([0., 0., 1.])

        if is3D:
            return matrix([  # global rod
                           [dot(v1x, g1x), dot(v1y, g1x), dot(v1z, g1x)],
                           [dot(v1x, g1y), dot(v1y, g1y), dot(v1z, g1y)],
                           [dot(v1x, g1z), dot(v1y, g1z), dot(v1z, g1z)],
                          ])  # rod

        return matrix([  # there can be no z component
                        [dot(v1x, g1x), dot(v1y, g1x)],
                        [dot(v1x, g1y), dot(v1y, g1y)],
                      ])  # rod

    def Lambda(self, model, debug=True):
        (n1, n2) = self.nodeIDs()
        return _Lambda(model, n1, n2, debug=debug)

    #def Stiffness1(self,model):
        #nIJV = [(nodes[0],1),(nodes[0],2),(nodes[1],1),]

    def Stiffness(self, model, node_ids, index0s, fnorm=1.0):  # CROD/CONROD
        #print("----------------")
        A = self.Area()
        E = self.E()
        G = self.G()
        J = self.J()

        #========================
        #(n1, n2) = self.nodeIDs()
        n0, n1 = self.nodeIDs()

        i0, i1 = index0s
        node0 = self.nodes[0]
        node1 = self.nodes[1]

        p0 = model.Node(n0).xyz
        p1 = model.Node(n1).xyz
        L = norm(p0 - p1)
        if L == 0.0:
            msg = 'invalid CROD length=0.0\n%s' % (self.__repr__())
            raise ZeroDivisionError(msg)
        #========================
        print("A=%g E=%g G=%g J=%g L=%g" % (A, E, G, J, L))
        k_axial = A * E / L
        k_torsion = G * J / L
        #k_axial = 1.0
        #k_torsion = 2.0

        k = matrix([[1., -1.], [-1., 1.]])  # 1D rod

        Lambda = self.Lambda(model)
        K = dot(dot(transpose(Lambda), k), Lambda)
        Ki, Kj = K.shape

        # for testing
        #K = ones((Ki, Ki), 'float64')

        K2 = zeros((Ki*2, Kj*2), 'float64')
        if k_axial == 0.0 and k_torsion == 0.0:
            dofs = []
            nIJV = []
            K2 = []
        elif k_torsion == 0.0: # axial; 2D or 3D
            K2 = K * k_axial
            dofs = array([
                i0, i0+1, i0+2,
                i1, i1+1, i1+2,
            ], 'int32')
            nIJV = [
                # axial
                (n0, 1), (n0, 2), (n0, 3),
                (n1, 1), (n1, 2), (n1, 3),
            ]
        elif k_axial == 0.0: # torsion; assume 3D
            K2 = K * k_torsion
            dofs = array([
                i0+3, i0+4, i0+5,
                i1+3, i1+4, i1+5,
            ], 'int32')
            nIJV = [
                # torsion
                (n0, 4), (n0, 5), (n0, 6),
                (n1, 4), (n1, 5), (n1, 6),
            ]

        else:  # axial + torsion; assume 3D
            # u1fx, u1fy, u1fz, u2fx, u2fy, u2fz
            K2[:Ki, :Ki] = K * k_axial

            # u1mx, u1my, u1mz, u2mx, u2my, u2mz
            K2[Ki:, Ki:] = K * k_torsion

            dofs = array([
                i0, i0+1, i0+2,
                i1, i1+1, i1+2,

                i0+3, i0+4, i0+5,
                i1+3, i1+4, i1+5,
            ], 'int32')
            nIJV = [
                # axial
                (n0, 1), (n0, 2), (n0, 3),
                (n1, 1), (n1, 2), (n1, 3),

                # torsion
                (n0, 4), (n0, 5), (n0, 6),
                (n1, 4), (n1, 5), (n1, 6),
            ]

        #Fg = dot(dot(transpose(Lambda), grav), Lambda)
        #print("K=\n", K / fnorm)
        #print("K2=\n", K2 / fnorm)

        #========================

        #print(K / fnorm)
        #print("K[%s] = \n%s\n" % (self.eid, list_print(K/fnorm)))

        #print('dofs =', dofs)
        #print('K =\n', list_print(K / fnorm))

        return(K2, dofs, nIJV)

    def Fg(self, model, grav, fnorm):
        n0, n1 = self.nodeIDs()
        mass = self.Mass()
        mg = -grav * mass

        #gx = array([1., 0., 0.])
        #gy = array([0., 1., 0.])
        #gz = array([0., 0., 1.])
        Fg = [mg[0], mg[1], mg[2],
              mg[0], mg[1], mg[2]]
        # bad
        nGrav = [(n0, 1), (n0, 2), (n0, 3),
                 (n1, 1), (n1, 2), (n1, 3)]

        #print("Fg = ",Fg)
        #print("mass = ", mass)
        #print("Fg =\n", Fg)
        #print("Fg[%s] = %s\n" % (self.eid, Fg))
        return (Fg, nGrav)

    def displacement_stress(self, model, q, dofs): # CROD/CONROD
        (n1, n2) = self.nodeIDs()
        Lambda = self.Lambda(model, debug=False)

        #print("**dofs =", dofs)
        n11 = dofs[(n1, 1)]
        n21 = dofs[(n2, 1)]

        n12 = dofs[(n1, 2)]
        n22 = dofs[(n2, 2)]

        n13 = dofs[(n1, 3)]
        n23 = dofs[(n2, 3)]

        # moments
        n14 = dofs[(n1, 4)]
        n24 = dofs[(n2, 4)]

        n15 = dofs[(n1, 5)]
        n25 = dofs[(n2, 5)]

        n16 = dofs[(n1, 6)]
        n26 = dofs[(n2, 6)]

        q_axial = array([
            q[n11], q[n12], q[n13],
            q[n21], q[n22], q[n23]
        ])
        q_torsion = array([
            q[n14], q[n15], q[n16],
            q[n24], q[n25], q[n26]
        ])
        #print("type=%s n1=%s n2=%s" % (self.type, n1, n2))
        #print("n11=%s n12=%s n21=%s n22=%s" %(n11,n12,n21,n22))

        #print("q2[%s] = %s" % (self.eid, q2))
        #print("Lambda = \n"+str(Lambda))

        #print "Lsize = ",Lambda.shape
        #print "qsize = ",q.shape
        u_axial = dot(array(Lambda), q_axial)
        du_axial = -u_axial[0] + u_axial[1]
        u_torsion = dot(array(Lambda), q_torsion)
        du_torsion = -u_torsion[0] + u_torsion[1]

        L = self.Length()
        E = self.E()
        A = self.Area()

        C = self.C()
        J = self.J()
        G = self.G()

        axial_strain = du_axial / L
        torsional_strain = du_torsion * C / L

        axial_stress = E * axial_strain
        torsional_stress = G * torsional_strain

        axial_force = axial_stress * A
        torsional_moment = du_torsion * G * J / L
        #print("axial_strain = %s [psi]" % axial_strain)
        #print("axial_stress = %s [psi]" % axial_stress)
        #print("axial_force  = %s [lb]\n" % axial_force)
        return (axial_strain, torsional_strain,
                axial_stress, torsional_stress,
                axial_force, torsional_moment)


class CROD(RodElement):
    type = 'CROD'
    _field_map = {
        1: 'eid', 2:'pid',
    }

    def _update_field_helper(self, n, value):
        if n == 3:
            self.nodes[0] = value
        elif n == 4:
            self.nodes[1] = value
        else:
            raise KeyError('Field %r=%r is an invalid %s entry.' % (n, value, self.type))

    def __init__(self, card=None, data=None, comment=''):
        RodElement.__init__(self, card, data)
        if comment:
            self._comment = comment
        if card:
            self.eid = integer(card, 1, 'eid')
            self.pid = integer_or_blank(card, 2, 'pid', self.eid)
            nids = [integer(card, 3, 'n1'),
                    integer(card, 4, 'n2')]
            assert len(card) == 5, 'len(CROD card) = %i' % len(card)
        else:
            self.eid = data[0]
            self.pid = data[1]
            nids = data[2:4]
        self.prepareNodeIDs(nids)
        assert len(self.nodes) == 2

    def _verify(self, xref=False):
        eid = self.Eid()
        pid = self.Pid()
        assert isinstance(pid, int), 'pid=%r' % pid
        if xref:  # True
            mid = self.Mid()
            L = self.Length()
            A = self.Area()
            nsm = self.Nsm()
            mpa = self.MassPerLength()
            mass = self.Mass()
            assert isinstance(mid, int), 'mid=%r' % mid
            assert isinstance(L, float), 'L=%r' % L
            assert isinstance(A, float), 'A=%r' % A
            assert isinstance(nsm, float), 'nsm=%r' % nsm
            assert isinstance(mpa, float), 'mass_per_length=%r' % mpa
            assert isinstance(mass, float), 'mass=%r' % mass

        c = self.Centroid()
        for i in range(3):
            assert isinstance(c[i], float), 'centroid[%i]=%r' % (i, c[i])

    def Centroid(self):
        return (self.nodes[0].Position() + self.nodes[1].Position()) / 2.

    def Eid(self):
        return self.eid

    def Mid(self):
        if isinstance(self.pid, int):
            raise RuntimeError('Element eid=%i has not been cross referenced.\n%s' % (self.eid, str(self)))
        return self.pid.Mid()

    def Area(self):
        if isinstance(self.pid, int):
            raise RuntimeError('Element eid=%i has not been cross referenced.\n%s' % (self.eid, str(self)))
        return self.pid.A

    def Nsm(self):
        if isinstance(self.pid, int):
            raise RuntimeError('Element eid=%i has not been cross referenced.\n%s' % (self.eid, str(self)))
        return self.pid.nsm

    def E(self):
        if isinstance(self.pid, int):
            raise RuntimeError('Element eid=%i has not been cross referenced.\n%s' % (self.eid, str(self)))
        return self.pid.mid.E()

    def G(self):
        if isinstance(self.pid, int):
            raise RuntimeError('Element eid=%i has not been cross referenced.\n%s' % (self.eid, str(self)))
        return self.pid.mid.G()

    def J(self):
        if isinstance(self.pid, int):
            raise RuntimeError('Element eid=%i has not been cross referenced.\n%s' % (self.eid, str(self)))
        return self.pid.J()

    def C(self):
        if isinstance(self.pid, int):
            raise RuntimeError('Element eid=%i has not been cross referenced.\n%s' % (self.eid, str(self)))
        return self.pid.c

    def MassPerLength(self):
        if isinstance(self.pid, int):
            raise RuntimeError('Element eid=%i has not been cross referenced.\n%s' % (self.eid, str(self)))
        massPerLength = self.pid.mid.rho * self.pid.A + self.pid.nsm
        return massPerLength

    def raw_fields(self):
        list_fields = ['CROD', self.eid, self.Pid()] + self.nodeIDs()
        return list_fields

    def repr_fields(self):
        return self.raw_fields()

    def write_bdf(self, size, card_writer):
        card = self.raw_fields()
        return self.comment() + print_card_8(card)


class CTUBE(RodElement):
    type = 'CTUBE'
    _field_map = {
        1: 'eid', 2:'pid',
    }

    def _update_field_helper(self, n, value):
        if n == 3:
            self.nodes[0] = value
        elif n == 4:
            self.nodes[1] = value
        else:
            raise KeyError('Field %r=%r is an invalid %s entry.' % (n, value, self.type))

    def __init__(self, card=None, data=None, comment=''):
        RodElement.__init__(self, card, data)
        if comment:
            self._comment = comment
        if card:
            self.eid = integer(card, 1, 'eid')
            self.pid = integer_or_blank(card, 2, 'pid', self.eid)
            nids = [integer(card, 3, 'n1'),
                    integer(card, 4, 'n2')]
            assert len(card) == 5, 'len(CTUBE card) = %i' % len(card)
        else:
            self.eid = data[0]
            self.pid = data[1]
            nids = data[2:4]
        self.prepareNodeIDs(nids)
        assert len(self.nodes) == 2

    def _verify(self, xref=False):
        pid = self.Pid()
        A = self.Area()
        assert isinstance(pid, int), 'pid=%r' % pid
        assert isinstance(A, float), 'A=%r' % A
        if xref:
            L = self.Length()
            nsm = self.Nsm()
            assert isinstance(L, float), 'L=%r' % L
            assert isinstance(nsm, float), 'nsm=%r' % nsm
            if self.pid.mid.type == 'MAT1':
                mpa = self.pid.MassPerLength()
                mass = self.Mass()
                assert isinstance(mpa, float), 'mass_per_length=%r' % mpa
                assert isinstance(mass, float), 'mass=%r' % mass
            elif self.pid.mid.type == 'MAT4':
                pass
            else:
                raise NotImplementedError('_verify does not support self.pid.mid.type=%s' % self.pid.mid.type)

        c = self.Centroid()
        for i in range(3):
            assert isinstance(c[i], float), 'centroid[%i]=%r' % (i, c[i])

    def Eid(self):
        return self.eid

    def Mid(self):
        return self.pid.Mid()

    def Mass(self):
        return self.pid.MassPerLength() * self.Length()

    def Nsm(self):
        if isinstance(self.pid, int):
            raise RuntimeError('Element eid=%i has not been cross referenced.\n%s' % (self.eid, str(self)))
        return self.pid.Nsm()

    def Area(self):
        if isinstance(self.pid, int):
            raise RuntimeError('Element eid=%i has not been cross referenced.\n%s' % (self.eid, str(self)))
        return self.pid.Area()

    def E(self):
        if isinstance(self.pid, int):
            raise RuntimeError('Element eid=%i has not been cross referenced.\n%s' % (self.eid, str(self)))
        return self.pid.mid.E()

    def G(self):
        if isinstance(self.pid, int):
            raise RuntimeError('Element eid=%i has not been cross referenced.\n%s' % (self.eid, str(self)))
        return self.pid.mid.G()

    def J(self):
        if isinstance(self.pid, int):
            raise RuntimeError('Element eid=%i has not been cross referenced.\n%s' % (self.eid, str(self)))
        return self.pid.J()

    def Centroid(self):
        return (self.nodes[0].Position() + self.nodes[1].Position()) / 2.

    def raw_fields(self):
        list_fields = ['CTUBE', self.eid, self.Pid()] + self.nodeIDs()
        return list_fields

    def write_bdf(self, size, card_writer):
        card = self.repr_fields()
        return self.comment() + print_card_8(card)


class CONROD(RodElement):
    type = 'CONROD'
    _field_map = {
        1: 'eid', 4:'mid', 5:'A', 6:'j', 7:'c', 8:'nsm',
    }

    def _update_field_helper(self, n, value):
        if n == 2:
            self.nodes[0] = value
        elif n == 3:
            self.nodes[1] = value
        else:
            raise KeyError('Field %r=%r is an invalid %s entry.' % (n, value, self.type))

    def __init__(self, card=None, data=None, comment=''):
        RodElement.__init__(self, card, data)
        if comment:
            self._comment = comment
        if card:
            self.eid = integer(card, 1, 'eid')
            nids = [integer(card, 2, 'n1'),
                    integer(card, 3, 'n2')]
            self.mid = integer(card, 4, 'mid')
            self.A = double(card, 5, 'A')
            self.j = double_or_blank(card, 6, 'j', 0.0)
            self.c = double_or_blank(card, 7, 'c', 0.0)
            self.nsm = double_or_blank(card, 8, 'nsm', 0.0)
        else:
            self.eid = data[0]
            nids = data[1:3]
            self.mid = data[3]
            self.A = data[4]
            self.j = data[5]
            self.c = data[6]
            self.nsm = data[7]
        self.prepareNodeIDs(nids)
        assert len(self.nodes) == 2
        #print self.nodes

    def cross_reference(self, model):
        msg = ' which is required by %s eid=%s' % (self.type, self.eid)
        self.nodes = model.Nodes(self.nodes, msg=msg)
        self.mid = model.Material(self.mid, msg=msg)

    def _verify(self, xref=False):
        pid = self.Pid()
        assert pid is None, 'pid=%r' % pid

        if xref:  # True
            mid = self.Mid()
            L = self.Length()
            A = self.Area()
            nsm = self.Nsm()
            mpa = self.MassPerLength()
            mass = self.Mass()
            assert isinstance(mid, int), 'mid=%r' % mid
            assert isinstance(L, float), 'L=%r' % L
            assert isinstance(A, float), 'A=%r' % A
            assert isinstance(nsm, float), 'nsm=%r' % nsm
            assert isinstance(mpa, float), 'mass_per_length=%r' % mpa
            assert isinstance(mass, float), 'mass=%r' % mass

            c = self.Centroid()
            for i in range(3):
                assert isinstance(c[i], float), 'centroid[%i]=%r' % (i, c[i])

    def Centroid(self):
        return (self.nodes[0].Position() + self.nodes[1].Position()) / 2.

    def Mid(self):
        if isinstance(self.mid, int):
            return self.mid
        #elif self.mid is None:
            #print ("No material defined for element ", self.eid)
            #return None
        else:
            return self.mid.mid

    def Eid(self):
        return self.eid

    def Pid(self):
        return None

    def MassPerLength(self):
        if isinstance(self.mid, int):
            raise RuntimeError('Element eid=%i has not been cross referenced.\n%s' % (self.eid, str(self)))
        massPerLength = self.mid.rho * self.A + self.nsm
        return massPerLength

    def C(self):
        """torsional constant"""
        return self.c

    def Area(self):
        return self.A

    def J(self):
        r"""returns the Polar Moment of Inertia, :math:`J`"""
        return self.j

    def Nsm(self):
        """Placeholder method for the non-structural mass"""
        return self.nsm

    def E(self):
        r"""returns the Young's Modulus, :math:`E`$"""
        if isinstance(self.mid, int):
            raise RuntimeError('Element eid=%i has not been cross referenced.\n%s' % (self.eid, str(self)))
        return self.mid.E()

    def G(self):
        r"""returns the Shear Modulus, :math:`G`"""
        if isinstance(self.mid, int):
            raise RuntimeError('Element eid=%i has not been cross referenced.\n%s' % (self.eid, str(self)))
        return self.mid.G()

    def Rho(self):
        r"""returns the material density, :math:`\rho`"""
        if isinstance(self.mid, int):
            raise RuntimeError('Element eid=%i has not been cross referenced.\n%s' % (self.eid, str(self)))
        return self.mid.rho

    def writeCodeAster(self):
        msg = ''
        msg += "    POUTRE=_F(GROUP_MA='CONROD_%s',\n" % self.eid
        msg += "              SECTION='CERCLE',  # circular section\n"
        if self.Thickness():
            msg += "              CARA=('R','EP'),   # radius, thickness\n"
            msg += "              VALE=(%g,%g),\n" % (
                self.Radius(), self.Thickness())
        else:
            msg += "              CARA=('R')   # radius\n"
            msg += "              VALE=(%g),\n" % self.Radius()
        return msg

    def raw_fields(self):
        list_fields = ['CONROD', self.eid] + self.nodeIDs() + [
                  self.Mid(), self.A, self.j, self.c, self.nsm]
        return list_fields

    def repr_fields(self):
        j = set_blank_if_default(self.j, 0.0)
        c = set_blank_if_default(self.c, 0.0)
        nsm = set_blank_if_default(self.nsm, 0.0)
        list_fields = ['CONROD', self.eid] + self.nodeIDs() + [self.Mid(),
                  self.A, j, c, nsm]
        return list_fields

    def write_bdf(self, size, card_writer):
        card = self.repr_fields()
        if size == 8:
            return self.comment() + print_card_8(card)
        return self.comment() + print_card_16(card)
        #return self.comment() + card_writer(card)
