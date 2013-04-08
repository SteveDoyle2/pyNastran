# pylint: disable=R0904,R0902,E1101,E1103,C0111,C0302,C0103,W0101
from __future__ import (nested_scopes, generators, division, absolute_import,
                        print_function, unicode_literals)
#import sys

from numpy import matrix, zeros, ones, array, transpose, dot
from numpy.linalg import norm

from pyNastran.bdf.fieldWriter import set_blank_if_default
from pyNastran.bdf.cards.baseCard import Element, Mid
from pyNastran.bdf.assign_type import (integer, integer_or_blank,
    integer_double_or_blank, double, double_or_blank,
    string_or_blank, integer_double_string_or_blank)

class RodElement(Element):  # CROD, CONROD, CTUBE
    def __init__(self, card, data):
        Element.__init__(self, card, data)

    def cross_reference(self, model):
        msg = ' which is required by %s eid=%s' % (self.type, self.eid)
        self.nodes = model.Nodes(self.nodes, msg=msg)
        self.pid = model.Property(self.pid, msg=msg)

    def Rho(self):
        r"""returns the material density  \f$ \rho \f$"""
        return self.pid.mid.rho

    def Length(self):
        r"""
        Returns the length of the element
        \f[ \large \sqrt{  (n_{x2}-n_{x1})^2+(n_{y2}-n_{y1})^2+(n_{z2}-n_{z1})^2  } \f]
        @param self the CROD/CONROD/CTUBE element
        """
        L = norm(self.nodes[1].Position() - self.nodes[0].Position())
        return L

    def Mass(self):
        r"""
        returns the mass of the element
        \f[ \large  mass = \left( \rho A + nsm \right) L  \f]
        """
        L = self.Length()
        mass = (self.Rho() * self.Area() + self.Nsm()) * L
        return mass

    def Rmatrix(self, model, is3D):
        r"""
        where   \f$ [R]_{ij} \f$ is the tranformation matrix
        \f[ \large  [R]_{ij} \left[
          \begin{array}{ccc}
              g_x \dot e_x & g_x \dot e_y &  g_x \dot e_z    \\
              g_y \dot e_x & g_y \dot e_y &  g_y \dot e_z    \\
              g_z \dot e_x & g_z \dot e_y &  g_z \dot e_z
          \end{array} \right]
        \f]
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
        """
        2d  [l,m,0,0]
            [0,0,l,m]

        3d  [l,m,n,0,0,0]
            [0,0,0,l,m,n]
        """
        is3D = False
        #R = self.Rmatrix(model,is3D)

        (n1, n2) = self.nodeIDs()
        p1 = model.Node(n1).Position()
        p2 = model.Node(n2).Position()
        v1 = p2 - p1
        if debug:
            print("v1=%s" % (v1))
        v1 = v1 / norm(v1)
        (l, m, n) = v1
        if is3D:
            Lambda = matrix(zeros((2, 6), 'd'))  # 3D
        else:
            Lambda = matrix(zeros((2, 4), 'd'))

        #print("R = \n",R)
        Lambda[0, 0] = Lambda[1, 2] = l
        Lambda[0, 1] = Lambda[1, 3] = m

        if is3D:
            Lambda[0, 2] = Lambda[1, 5] = n  # 3D
        if debug:
            print("Lambda = \n" + str(Lambda))
        return Lambda

    #def Stiffness1(self,model):
        #nIJV = [(nodes[0],1),(nodes[0],2),(nodes[1],1),]

    def Stiffness(self, model, grav, is3D):  # CROD/CONROD
        print("----------------")
        Lambda = self.Lambda(model)
        #print("Lambda = \n"+str(Lambda))

        k = self.Stiffness1D(model)  # /250000.
        #print R
        #print(k)
        #print("Lambda.shape = ",Lambda.shape)
        K = dot(dot(transpose(Lambda), k), Lambda)
        #Fg = dot(dot(transpose(Lambda),grav),Lambda)

        #print('size(grav) =',grav.shape)
        mass = self.Mass()
        mg = -grav * mass
        if is3D:
            Fg = [mg[0], mg[1], mg[2], mg[0], mg[1], mg[2]]
        else:
            Fg = [mg[0], mg[1], mg[0], mg[1]]
        #print("Fg = ",Fg)
        print("mass = ", mass)
        print("Fg   = ", Fg)
        #print(K)
        print("K[%s] = \n%s\n" % (self.eid, K))
        #print("Fg[%s] = %s\n" %(self.eid,Fg))

        nodes = self.nodeIDs()
        nIJV = [(nodes[0], 1), (nodes[0], 2),
                (nodes[1], 1), (nodes[1], 2), ]
        nGrav = [(nodes[0], 1), (nodes[0], 2), (nodes[0], 3),
                 (nodes[1], 1), (nodes[1], 2), (nodes[1], 3)]

        return(K, nIJV, Fg, nGrav)

    def displacementStressF06(self, model, q, dofs):
        pass

    def displacementStress(self, model, q, dofs):
        (n1, n2) = self.nodeIDs()
        Lambda = self.Lambda(model, debug=False)
        n11 = dofs[(n1, 1)]
        n12 = dofs[(n1, 2)]
        n21 = dofs[(n2, 1)]
        n22 = dofs[(n2, 2)]
        print("type=%s n1=%s n2=%s" % (self.type, n1, n2))
        #print("n11=%s n12=%s n21=%s n22=%s" %(n11,n12,n21,n22))

        q2 = array([q[n11], q[n12], q[n21], q[n22]])
        print("q[%s] = %s" % (self.eid, q2))
        #print("Lambda = \n"+str(Lambda))

        #print "Lsize = ",Lambda.shape
        #print "qsize = ",q.shape
        u = dot(array(Lambda), q2)
        #L = self.Length()
        EL = self.E() / self.Length()

        #stressX = -EL*u[0]+EL*u[1]
        stressX = EL * (-u[0] + u[1])
        print("stressX = %s [psi]\n" % (stressX))
        return stressX

    def Stiffness1D(self, model):  # CROD/CONROD
        """
        @todo remove this method after making sure it still works
        """
        #L = norm(r)
        (n1, n2) = self.nodeIDs()
        #node1 = model.Node(n1)
        #node2 = model.Node(n2)

        p1 = model.Node(n1).xyz
        p2 = model.Node(n2).xyz
        #print "node1 = ",node1
        #print "node2 = ",node2
        L = norm(p1 - p2)

        if L == 0.0:
            msg = 'invalid CROD length=0.0\n%s' % (self.__repr__())
            raise ZeroDivisionError(msg)

        A = self.Area()
        #mat = self.mid
        E = self.E()
        print("A = ", A)
        print("E = ", E)
        print("L = ", L)
        #ki = 1.
        ki = A * E / L
        #ki = 250000.
        #knorm = 250000.
        K = ki * matrix([[1., -1.], [-1., 1.]])  # rod

        print("A=%g E=%g L=%g  AE/L=%g" % (A, E, L, A * E / L))
        #print "K = \n",K
        return K


class LineElement(Element):  # CBAR, CBEAM, CBEAM3, CBEND
    def __init__(self, card, data):
        Element.__init__(self, card, data)

    def C(self):
        """torsional constant"""
        return self.pid.C()

    def Area(self):
        """returns the area of the element face"""
        raise NotImplementedError('implement self.Area() for %s' % self.type)

    def E(self):
        r"""returns the Young's Modulus  \f$ E \f$"""
        return self.pid.mid.E()

    def G(self):
        r"""returns the Shear Modulus   \f$ G \f$"""
        return self.pid.mid.G()

    def J(self):
        r"""returns the Polar Moment of Inertia.   \f$ J \f$"""
        return self.pid.J()

    def I11(self):
        r"""returns the Moment of Inertia.   \f$ I_{11} \f$"""
        return self.pid.I11()

    def I22(self):
        r"""returns the Moment of Inertia.   \f$ I_{22} \f$"""
        return self.pid.I22()

    def I12(self):
        r"""returns the Moment of Inertia.   \f$ I_{12} \f$"""
        return self.pid.I12()

    def Nu(self):
        r"""returns Poisson's Ratio  \f$ \nu \f$"""
        return self.pid.mid.nu

    def Rho(self):
        r"""returns the material density  \f$ \rho \f$"""
        #print str(self.pid),type(self.pid)
        #raise NotImplementedError('implement self.Rho() for %s' % self.type)
        return self.pid.mid.rho

    def Nsm(self):
        """Placeholder method for the non-structural mass"""
        raise NotImplementedError('implement self.Area() for %s' % self.type)

    def MassPerLength(self):
        """Returns the mass per unit length"""
        return self.pid.MassPerLength()

    def Mass(self):
        r"""
        returns the mass of the element

        \f[ \large  mass = \left( \rho A + nsm \right) L  \f]
        """
        L = self.Length()
        try:
            mass = (self.Rho() * self.Area() + self.Nsm()) * L
        except TypeError:
            msg = 'TypeError on eid=%s pid=%s:\n' % (self.eid, self.Pid())
            msg += 'rho = %s\narea = %s\nnsm = %s\nL = %s' % (self.Rho(),
                                                              self.Area(),
                                                              self.Nsm(), L)
            raise TypeError(msg)

        return mass

    def cross_reference(self, model):
        msg = ' which is required by %s eid=%s' % (self.type, self.eid)
        self.nodes = model.Nodes(self.nodes, msg=msg)
        self.pid = model.Property(self.pid, msg=msg)

    def Length(self):
        r"""
        Returns the length of the element
        \f[ \large \sqrt{  (n_{x2}-n_{x1})^2+(n_{y2}-n_{y1})^2+(n_{z2}-n_{z1})^2  } \f]
        @param self the object pointer
        """
        L = norm(self.nodes[1].Position() - self.nodes[0].Position())
        return L

    # def k_Axial(self):
    #     r"""
    #     Returns the axial stiffness matrix.
    #
    #     \f[ \large   k_{Axial} = \frac{AE}{2L}
    #       \left[
    #       \begin{array}{cc}
    #           1 & -1 \\
    #          -1 &  1
    #       \end{array} \right]
    #     \f]
    #     """
    #     raise NotImplementedError()
    #     L = self.Length()
    #     E = self.E()
    #     A = self.Area()
    #     kMag = A * E / (2 * L)
    #     K = ones(1, 1)
    #     K[0, 1] = K[1, 0] = -1
    #     return kMag * K

    # def k_Torsion(self):  # not done
    #     r"""
    #     Returns the torsional stiffness matrix.
    #
    #     \f[ \large   k_{Axial} = \frac{L}{GJ}
    #       \left[
    #       \begin{array}{cc}
    #           1 & -1 \\
    #          -1 &  1
    #       \end{array} \right]
    #     \f]
    #     @warning formula not verified
    #     """
    #     raise NotImplementedError()
    #     L = self.Length()
    #     G = self.G()
    #     J = self.J()
    #     #A = self.Area()
    #     #kMag = A*E/(2*L)
    #     kMag = L / (G * J)
    #     K = ones(1, 1)
    #     K[0, 1] = K[1, 0] = -1
    #     return kMag * K

    # def k_Bending(self):
    #     r"""
    #     Returns the bending stiffness matrix.
    #
    #     \f[ \large  k_{Bending} = \frac{EI}{L^3}
    #       \left[
    #       \begin{array}{cccc}
    #          12 &  6L   & -12 &  6L    \\
    #          6L &  4L^2 & -6L &  2L^2  \\
    #         -12 & -6L   &  12 & -6L    \\
    #          6L &  2L^2 & -6L &  4L^2
    #       \end{array} \right]
    #     \f]
    #     """
    #     raise NotImplementedError()
    #     L = self.Length()
    #     E = self.E()
    #     I = self.I()
    #     LL = L * L
    #     LLL = L * LL
    #     sL = 6 * L
    #     tLL = 2 * LL
    #     fLL = 4 * LL
    #     kMag = E * I / LLL
    #
    #     #K = Matrix(zeros(4,4))
    #     K = matrix([[12., sL, -12, sL],
    #                 [sL, fLL, -sL, tLL],
    #                 [-12, -sL, 12., -sL],
    #                 [sL, tLL, -sL, fLL]])
    #     #M[1,0] = sL
    #     #M[2,0] = -12.
    #     #M[3,0] = sL
    #
    #     #M[2,4] =  -sL
    #     #M[1,1] = M[3,3] = fLL
    #     return kMag * K


class CROD(RodElement):
    type = 'CROD'
    def __init__(self, card=None, data=None, comment=''):
        RodElement.__init__(self, card, data)
        if comment:
            self._comment = comment
        if card:
            self.eid = integer(card, 1, 'eid')
            self.pid = integer_or_blank(card, 2, 'pid', self.eid)
            nids = [integer(card, 3, 'n1'),
                    integer(card, 4, 'n2')]
            assert len(card) <= 4, 'len(CROD card) = %i' % len(card)
        else:
            self.eid = data[0]
            self.pid = data[1]
            nids = data[2:4]
        self.prepareNodeIDs(nids)
        assert len(self.nodes) == 2

    def Centroid(self):
        return (self.nodes[0].Position() + self.nodes[1].Position()) / 2.

    def Mid(self):
        return self.pid.Mid()

    def Area(self):
        return self.pid.A

    def Nsm(self):
        return self.pid.nsm

    def MassPerLength(self):
        massPerLength = self.pid.mid.rho * self.pid.A + self.pid.nsm
        return massPerLength

    def rawFields(self):
        list_fields = ['CROD', self.eid, self.Pid()] + self.nodeIDs()
        return list_fields

    def reprFields(self):
        return self.rawFields()


class CTUBE(RodElement):
    type = 'CTUBE'

    def __init__(self, card=None, data=None, comment=''):
        RodElement.__init__(self, card, data)
        if comment:
            self._comment = comment
        if card:
            self.eid = integer(card, 1, 'eid')
            self.pid = integer_or_blank(card, 2, 'pid', self.eid)
            nids = [integer(card, 3, 'n1'),
                    integer(card, 4, 'n2')]
            assert len(card) <= 4, 'len(CTUBE card) = %i' % len(card)
        else:
            self.eid = data[0]
            self.pid = data[1]
            nids = data[2:4]
        self.prepareNodeIDs(nids)
        assert len(self.nodes) == 2

    def Area(self):
        return self.pid.Area()

    def Centroid(self):
        """@todo improve the formuala for CTUBE centroid"""
        return (self.nodes[0].Position() + self.nodes[1].Position()) / 2.

    def rawFields(self):
        list_fields = ['CTUBE', self.eid, self.Pid()] + self.nodeIDs()
        return list_fields


class CONROD(RodElement):
    type = 'CONROD'

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

    def Centroid(self):
        return (self.nodes[0].Position() + self.nodes[1].Position()) / 2.

    def Mid(self):
        return Mid(self)

    def Pid(self):
        return None

    def MassPerLength(self):
        massPerLength = self.mid.rho * self.A + self.nsm
        return massPerLength

    def C(self):
        """torsional constant"""
        return self.c

    def Area(self):
        return self.A

    def J(self):
        r"""returns the Polar Moment of Inertia.   \f$ J \f$"""
        return self.j

    def Nsm(self):
        """Placeholder method for the non-structural mass"""
        return self.nsm

    def E(self):
        r"""returns the Young's Modulus  \f$ E \f$"""
        return self.mid.E()

    def G(self):
        r"""returns the Shear Modulus   \f$ G \f$"""
        return self.mid.G()

    def Nu(self):
        r"""returns Poisson's Ratio  \f$ \nu \f$"""
        return self.mid.nu

    def Rho(self):
        r"""returns the material density  \f$ \rho \f$"""
        return self.mid.rho

    def writeCodeAster(self):
        msg = ''
        msg += "    POUTRE=_F(GROUP_MA='CONROD_%s',\n" % (self.eid)
        msg += "              SECTION='CERCLE',  # circular section\n"
        if self.Thickness():
            msg += "              CARA=('R','EP'),   # radius, thickness\n"
            msg += "              VALE=(%g,%g),\n" % (
                self.Radius(), self.Thickness())
        else:
            msg += "              CARA=('R')   # radius\n"
            msg += "              VALE=(%g),\n" % (self.Radius())
        return msg

    def rawFields(self):
        list_fields = ['CONROD', self.eid] + self.nodeIDs() + [
                  self.Mid(), self.A, self.j, self.c, self.nsm]
        return list_fields

    def reprFields(self):
        j = set_blank_if_default(self.j, 0.0)
        c = set_blank_if_default(self.c, 0.0)
        nsm = set_blank_if_default(self.nsm, 0.0)
        list_fields = ['CONROD', self.eid] + self.nodeIDs() + [self.Mid(),
                  self.A, j, c, nsm]
        return list_fields


class CBAR(LineElement):
    """
    CBAR EID PID GA GB X1 X2 X3 OFFT
    PA PB W1A W2A W3A W1B W2B W3B
    or
    CBAR EID PID GA GB G0 OFFT
    PA PB W1A W2A W3A W1B W2B W3B

    CBAR       22062       4   21798   21799   0.0+0   0.0+0     -1.
                               0.0+0   0.0+0     -9.   0.0+0   0.0+0     -9.
    """
    type = 'CBAR'
    asterType = 'CBAR'
    def __init__(self, card=None, data=None, comment=''):
        LineElement.__init__(self, card, data)
        if comment:
            self._comment = comment
        if card:
            self.eid = integer(card, 1, 'eid')
            self.pid = integer_or_blank(card, 2, 'pid', self.eid)
            self.ga = integer(card, 3, 'ga')
            self.gb = integer(card, 4, 'gb')
            self.initX_G0(card)

            self.offt = string_or_blank(card, 8, 'offt', 'GGG')
            #print 'self.offt = |%s|' %(self.offt)

            self.pa = integer_or_blank(card, 9, 'pa', 0)
            self.pb = integer_or_blank(card, 10, 'pb', 0)

            self.w1a = double_or_blank(card, 11, 'w1a', 0.0)
            self.w2a = double_or_blank(card, 12, 'w2a', 0.0)
            self.w3a = double_or_blank(card, 13, 'w3a', 0.0)

            self.w1b = double_or_blank(card, 14, 'w1b', 0.0)
            self.w2b = double_or_blank(card, 15, 'w2b', 0.0)
            self.w3b = double_or_blank(card, 16, 'w3b', 0.0)
            assert len(card) <= 16, 'len(CBAR card) = %i' % len(card)
        else:  ## @todo verify
            #data = [[eid,pid,ga,gb,pa,pb,w1a,w2a,w3a,w1b,w2b,w3b],[f,g0]]
            #data = [[eid,pid,ga,gb,pa,pb,w1a,w2a,w3a,w1b,w2b,w3b],[f,x1,x2,x3]]

            main = data[0]
            flag = data[1][0]
            if flag in [0, 1]:
                self.g0 = None
                self.x1 = data[1][1]
                self.x2 = data[1][2]
                self.x3 = data[1][3]
            else:
                self.g0 = data[1][1]
                self.x1 = None
                self.x2 = None
                self.x3 = None

            self.eid = main[0]
            self.pid = main[1]
            self.ga = main[2]
            self.gb = main[3]
            #self.offt = str(data[4]) # GGG
            self.offt = 'GGG'  ## @todo offt can be an integer; translate to char
            self.pa = main[4]
            self.pb = main[5]

            self.w1a = main[6]
            self.w2a = main[7]
            self.w3a = main[8]

            self.w1b = main[9]
            self.w2b = main[10]
            self.w3b = main[11]
        #print("offt = %s" %(self.offt))
        if not isinstance(self.offt, unicode) or isinstance(self.offt, str):
            raise SyntaxError('invalid offt expected a string of length 3 '
                              'offt=|%s|' % self.offt)

        msg = 'invalid offt parameter of %s...offt=%s' % (self.type, self.offt)
        # B,G,O
        assert self.offt[0] in ['G', 'B'], msg
        assert self.offt[1] in ['G', 'O', 'E'], msg
        assert self.offt[2] in ['G', 'O', 'E'], msg

    def Mid(self):
        return self.pid.Mid()

    def Area(self):
        A = self.pid.Area()
        assert isinstance(A, float)
        return A

    def Length(self):
        L = norm(self.gb.Position() - self.ga.Position())
        assert isinstance(L, float)
        return L

    def Nsm(self):
        nsm = self.pid.Nsm()
        assert isinstance(nsm, float)
        return nsm

    def I1(self):
        return self.pid.I1()

    def I2(self):
        return self.pid.I2()

    def Centroid(self):
        return (self.ga.Position() + self.gb.Position()) / 2.

    def initX_G0(self, card):
        field5 = integer_double_or_blank(card, 5, 'g0_x1')
        if isinstance(field5, int):
            self.g0 = field5
            self.x1 = None
            self.x2 = None
            self.x3 = None
        elif isinstance(field5, float):
            self.g0 = None
            self.x1 = field5
            self.x2 = double_or_blank(card, 6, 'x2', 0.0)
            self.x3 = double_or_blank(card, 7, 'x3', 0.0)
        else:
            #msg = 'field5 on %s is the wrong type...id=%s field5=%s '
            #      'type=%s' %(self.type,self.eid,field5,type(field5))
            #raise InvalidFieldError(msg)
            self.g0 = None
            self.x1 = 0.
            self.x2 = 0.
            self.x3 = 0.
        #if self.eid==14100238:
            #print "g0=%s x1=%s x2=%s x3=%s" %(self.g0, self.x1, self.x2,
            #                                  self.x3)

    def cross_reference(self, mesh):
        """
        set g0-ga to x1,x2,x3
        """
        #if self.g0:
        #    v = nodes[self.g0].Position()-nodes[self.ga].Position()
        #    self.x1 = v[0]
        #    self.x2 = v[1]
        #    self.x3 = v[2]
        msg = ' which is required by %s eid=%s' % (self.type, self.eid)
        self.ga = mesh.Node(self.ga, msg=msg)
        self.gb = mesh.Node(self.gb, msg=msg)
        self.pid = mesh.Property(self.pid, msg=msg)

    #def updateNodes(self,nodes):
    #    """@todo maybe improve"""
    #    self.cross_reference(self,nodes)

    def Ga(self):
        if isinstance(self.ga, int):
            return self.ga
        else:
            return self.ga.nid

    def Gb(self):
        if isinstance(self.gb, int):
            return self.gb
        else:
            return self.gb.nid

    def getX_G0_defaults(self):
        if self.g0:
            return (self.g0, None, None)
        else:
            #x1 = set_blank_if_default(self.x1, 0.0)
            #x2 = set_blank_if_default(self.x2, 0.0)
            #x3 = set_blank_if_default(self.x3, 0.0)
            return (self.x1, self.x2, self.x3)

    def nodeIDs(self):
        return [self.Ga(), self.Gb()]

    def Stiffness(self, model, grav, is3D):  # CBAR
        print("----------------")
        Lambda = self.Lambda(model)
        #print("Lambda = \n"+str(Lambda))

        k = self.Stiffness1D(model)  # /250000.
        #print R
        #print(k)
        print("Lambda.shape = ", Lambda.shape)
        K = dot(dot(transpose(Lambda), k), Lambda)
        #Fg = dot(dot(transpose(Lambda),grav),Lambda)

        #print('size(grav) =',grav.shape)
        mass = self.Mass()
        mg = -grav * mass
        if is3D:
            Fg = [mg[0], mg[1], mg[2], mg[0], mg[1], mg[2]]
        else:
            Fg = [mg[0], mg[1], mg[0], mg[1]]
        #print("Fg = ",Fg)
        print("mass = ", mass)
        print("Fg   = ", Fg)
        #print(K)
        print("K[%s] = \n%s\n" % (self.eid, K))
        #print("Fg[%s] = %s\n" %(self.eid,Fg))

        n0, n1 = self.nodeIDs()
                 # u1          v1          theta1         u2,v2,w2
                 
        nIJV = [
                 (n0, 1), (n0, 2), (n0, 5),  (n1, 1), (n1, 2), (n1, 5),  # X1
                 (n0, 1), (n0, 2), (n0, 5),  (n1, 1), (n1, 2), (n1, 5),  # Y1
                 (n0, 1), (n0, 2), (n0, 5),  (n1, 1), (n1, 2), (n1, 5),  # M1

                 (n0, 1), (n0, 2), (n0, 5),  (n1, 1), (n1, 2), (n1, 5),  # X2
                 (n0, 1), (n0, 2), (n0, 5),  (n1, 1), (n1, 2), (n1, 5),  # Y2
                 (n0, 1), (n0, 2), (n0, 5),  (n1, 1), (n1, 2), (n1, 5),  # M2
                ]
        #print("nIJV = ",nIJV)
        nGrav = [(n0, 1), (n0, 2), (n0, 3),
                 (n1, 1), (n1, 2), (n1, 3)]

        return(K, nIJV, Fg, nGrav)

    def Stiffness1D(self, model):  # CBAR
        #L = norm(r)
        (n1, n2) = self.nodeIDs()
        #node1 = model.Node(n1)
        #node2 = model.Node(n2)

        p1 = model.Node(n1).xyz
        p2 = model.Node(n2).xyz
        #print "node1 = ",node1
        #print "node2 = ",node2
        L = norm(p1 - p2)

        if L == 0.0:
            msg = 'invalid CBAR length=0.0\n%s' % (self.__repr__())
            raise ZeroDivisionError(msg)

        A = self.Area()
        #mat = self.mid
        E = self.E()
        I1 = self.I1()
        I2 = self.I2()
        print("A  = ", A)
        print("E  = ", E)
        print("L  = ", L)
        print("I1 = ", I1)
        print("I2 = ", I2)
        #ki = 1.
        #ki = A*E/L
        #ki = 250000.
        #knorm = 250000.
        #K = ki*matrix([[1.,-1.],[-1.,1.]]) # rod
        #P = 10.416666667
        P = 1.
        k1 = A * E / L
        k2 = P * L ** 3 / (E * I1)
        print("A=%g E=%g L=%g I1=%g I2=%G AE/L=%g L^3/E*Iz=%g" % (
            A, E, L, I1, I2, k1, k2))
        #print "K = \n",K
        #return K

        L = self.Length()
        E = self.E()
        #I = self.I1()

        LL = L * L
        #LLL = L*LL
        #kMag = E*I/LLL
        sL = 6 * L
        tLL = 2 * LL
        fLL = 4 * LL

        #K = zeros(4,4))
        K = array([[0., 12., sL, 0., -12, sL],
                   [0., sL, fLL, 0., -sL, tLL],
                   [0., 0, 0, 0., 0, 0],
                   [0., -12, -sL, 0., 12., -sL],
                   [0., sL, tLL, 0., -sL, fLL],
                   [0., 0, 0, 0., 0, 0],
                   ])
        print("k =\n" + str(K))
        return K

    def Lambda(self, model, debug=True):  # CBAR from CROD/CONROD
        """
        2d  [l,m,0,0,  l,m,0,0]
            [0,0,l,m,  0,0,l,m] L*k = 2x4*W4x4

        3d  [l,m,n,0,0,0,  l,m,n,0,0,0]
            [0,0,0,l,m,n,  0,0,0,l,m,n]
        """
        is3D = False
        #R = self.Rmatrix(model,is3D)

        (n1, n2) = self.nodeIDs()
        p1 = model.Node(n1).Position()
        p2 = model.Node(n2).Position()
        v1 = p2 - p1
        if debug:
            print("v1=%s" % (v1))
        v1 = v1 / norm(v1)
        #v2 = v2/norm(v2)
        #v3 = v3/norm(v3)
        (l, m, n) = v1
        #(o,p,q) = v2
        #(r,s,t) = v3
        print("l=%s m=%s n=%s" % (l, m, n))
        if is3D:
            Lambda = matrix([[l, m, n, 0, 0, 0, ],
                             [m, l, n, 0, 0, 0, ],
                             [0, 0, 1, 0, 0, 0, ],
                             [0, 0, 0, l, m, n, ],
                             [0, 0, 0, m, l, n, ],
                             [0, 0, 0, 0, 0, 1, ],
                             ])
        else:
            Lambda = matrix([[l, m, n, 0, 0, 0, ],
                             [m, l, n, 0, 0, 0, ],
                             [0, 0, 1, 0, 0, 0, ],
                             [0, 0, 0, l, m, n, ],
                             [0, 0, 0, m, l, n, ],
                             [0, 0, 0, 0, 0, 1, ],
                             ])
        if debug:
            print("Lambda* = \n" + str(Lambda))
        return Lambda

    def rawFields(self):
        """@todo not perfectly accurate"""
        (x1, x2, x3) = self.getX_G0_defaults()
        offt = set_blank_if_default(self.offt, 'GGG')
        list_fields = ['CBAR', self.eid, self.Pid(), self.Ga(), self.Gb(), x1, x2,
                  x3, offt, self.pa, self.pb, self.w1a, self.w2a, self.w3a,
                  self.w1b, self.w2b, self.w3b]
        return list_fields

    def reprFields(self):
        pa = set_blank_if_default(self.pa, 0)
        pb = set_blank_if_default(self.pb, 0)

        w1a = set_blank_if_default(self.w1a, 0.0)
        w2a = set_blank_if_default(self.w2a, 0.0)
        w3a = set_blank_if_default(self.w3a, 0.0)
        w1b = set_blank_if_default(self.w1b, 0.0)
        w2b = set_blank_if_default(self.w2b, 0.0)
        w3b = set_blank_if_default(self.w3b, 0.0)
        (x1, x2, x3) = self.getX_G0_defaults()
        offt = set_blank_if_default(self.offt, 'GGG')
        list_fields = ['CBAR', self.eid, self.Pid(), self.Ga(), self.Gb(), x1, x2,
                  x3, offt, pa, pb, w1a, w2a, w3a, w1b, w2b, w3b]
        return list_fields


class CBEAM3(CBAR):
    """
    Defines a three-node beam element
    """
    type = 'CBEAM3'

    def __init__(self, card=None, data=None, comment=''):
        LineElement.__init__(self, card, data)
        if comment:
            self._comment = comment
        if card:
            self.eid = integer(card, 1, 'eid')
            self.pid = integer(card, 2, 'pid', self.eid)
            self.ga = integer(card, 3, 'ga')
            self.gb = integer(card, 4, 'gb')
            self.gc = integer(card, 5, 'gc')

            self.initX_G0(card)

            self.w1a = double_or_blank(card, 9, 'w1a', 0.0)
            self.w2a = double_or_blank(card, 10, 'w2a', 0.0)
            self.w3a = double_or_blank(card, 11, 'w3a', 0.0)

            self.w1b = double_or_blank(card, 12, 'w1b', 0.0)
            self.w2b = double_or_blank(card, 13, 'w2b', 0.0)
            self.w3b = double_or_blank(card, 14, 'w3b', 0.0)

            self.w1c = double_or_blank(card, 15, 'w1c', 0.0)
            self.w2c = double_or_blank(card, 16, 'w2c', 0.0)
            self.w3c = double_or_blank(card, 17, 'w3c', 0.0)

            self.twa = double_or_blank(card, 18, 0., 'twa')
            self.twb = double_or_blank(card, 19, 0., 'twb')
            self.twc = double_or_blank(card, 20, 0., 'twc')

            self.sa = integer_or_blank(card, 21, 'sa')
            self.sb = integer_or_blank(card, 22, 'sb')
            self.sc = integer_or_blank(card, 23, 'sc')
            assert len(card) <= 23, 'len(CBEAM3 card) = %i' % len(card)
        else:
            raise NotImplementedError(data)

    def cross_reference(self, model):
        msg = ' which is required by %s eid=%s' % (self.type, self.eid)
        self.ga = model.Node(self.ga, msg=msg)
        self.gb = model.Node(self.gb, msg=msg)
        self.gc = model.Node(self.gc, msg=msg)
        self.pid = model.Property(self.pid, msg=msg)

    def Length(self):
        """
        L = gb-ga
        @todo improve formula
        """
        L = norm(self.gb.Position() - self.ga.Position())
        assert isinstance(L, float)
        return L

    def rawFields(self):
        (x1, x2, x3) = self.getX_G0_defaults()
        (ga, gb, gc) = self.nodeIDs()
        list_fields = ['CBEAM3', self.eid, self.Pid(), ga, gb, gc, x1, x2, x3,
                  self.w1a, self.w2a, self.w3a, self.w1b, self.w2b, self.w3b,
                  self.w1c, self.w2c, self.w3c, self.twa, self.twb, self.twc,
                  self.sa, self.sb, self.sc]
        return list_fields

    def reprFields(self):
        w1a = set_blank_if_default(self.w1a, 0.0)
        w2a = set_blank_if_default(self.w2a, 0.0)
        w3a = set_blank_if_default(self.w3a, 0.0)
        w1b = set_blank_if_default(self.w1b, 0.0)
        w2b = set_blank_if_default(self.w2b, 0.0)
        w3b = set_blank_if_default(self.w3b, 0.0)
        w1c = set_blank_if_default(self.w1c, 0.0)
        w2c = set_blank_if_default(self.w2c, 0.0)
        w3c = set_blank_if_default(self.w3c, 0.0)

        twa = set_blank_if_default(self.twa, 0.0)
        twb = set_blank_if_default(self.twb, 0.0)
        twc = set_blank_if_default(self.twc, 0.0)

        (x1, x2, x3) = self.getX_G0_defaults()
        (ga, gb, gc) = self.nodeIDs()
        list_fields = ['CBEAM3', self.eid, self.Pid(), ga, gb, x1, x2, x3,
                  w1a, w2a, w3a, w1b, w2b, w3b, w1c, w2c, w3c,
                  twa, twb, twc, self.sa, self.sb, self.sc]
        return list_fields


class CBEAM(CBAR):
    """
    CBEAM EID PID GA GB X1 X2 X3 OFFT/BIT
    PA PB W1A W2A W3A W1B W2B W3B
    SA SB
    or
    CBEAM EID PID GA GB G0 - - OFFT/BIT
    PA PB W1A W2A W3A W1B W2B W3B
    SA SB
    """
    type = 'CBEAM'

    def __init__(self, card=None, data=None, comment=''):
        LineElement.__init__(self, card, data)
        if comment:
            self._comment = comment
        if card:
            self.eid = integer(card, 1, 'eid')
            self.pid = integer_or_blank(card, 2, 'pid', self.eid)
            self.ga = integer(card, 3, 'ga')
            self.gb = integer(card, 4, 'gb')

            self.initX_G0(card)
            self.initOfftBit(card)
            self.pa = integer_or_blank(card, 9, 'pa')
            self.pb = integer_or_blank(card, 10, 'pb')

            self.w1a = double_or_blank(card, 11, 'w1a', 0.0)
            self.w2a = double_or_blank(card, 12, 'w2a', 0.0)
            self.w3a = double_or_blank(card, 13, 'w3a', 0.0)

            self.w1b = double_or_blank(card, 14, 'w1b', 0.0)
            self.w2b = double_or_blank(card, 15, 'w2b', 0.0)
            self.w3b = double_or_blank(card, 16, 'w3b', 0.0)

            self.sa = integer_or_blank(card, 17, 'sa')
            self.sb = integer_or_blank(card, 18, 'sb')
            assert len(card) <= 18, 'len(CBEAM card) = %i' % len(card)
        else:  ## @todo verify
            #data = [[eid,pid,ga,gb,sa,sb, pa,pb,w1a,w2a,w3a,w1b,w2b,w3b],
            #        [f,g0]]
            #data = [[eid,pid,ga,gb,sa,sb, pa,pb,w1a,w2a,w3a,w1b,w2b,w3b],
            #        [f,x1,x2,x3]]

            main = data[0]

            flag = data[1][0]
            if flag in [0, 1]:
                self.g0 = None
                self.x1 = data[1][1]
                self.x2 = data[1][2]
                self.x3 = data[1][3]
            else:
                self.g0 = data[1][1]
                self.x1 = None
                self.x2 = None
                self.x3 = None

            self.eid = main[0]
            self.pid = main[1]
            self.ga = main[2]
            self.gb = main[3]
            self.sa = main[4]
            self.sb = main[5]

            self.isOfft = True  ## @todo is this correct???
            #self.offt = str(data[6]) # GGG
            self.offt = 'GGG'  ## @todo is this correct???

            self.pa = main[6]
            self.pb = main[7]

            self.w1a = main[8]
            self.w2a = main[9]
            self.w3a = main[10]

            self.w1b = main[11]
            self.w2b = main[12]
            self.w3b = main[13]

    def initOfftBit(self, card):
        field8 = integer_double_string_or_blank(card, 8, 'field8')
        if isinstance(field8, float):
            self.isOfft = False
            self.offt = None
            self.bit = field8
        elif field8 is None:
            self.isOfft = True
            self.offt = 'GGG'  # default
            self.bit = None
        elif isinstance(field8, unicode):
            self.isOfft = True
            self.bit = None
            self.offt = field8
            #print("self.offt = ", self.offt)
            msg = 'invalid offt parameter of CBEAM...offt=%s' % self.offt
            assert self.offt[0] in ['G', 'B', 'O', 'E'], msg
            assert self.offt[1] in ['G', 'B', 'O', 'E'], msg
            assert self.offt[2] in ['G', 'B', 'O', 'E'], msg
        else:
            msg = ('field8 on %s card is not a string(offt) or bit '
                   '(float)...field8=%s\n' % (self.type, field8))
            raise RuntimeError("Card Instantiation: %s" % msg)

    def Mid(self):
        return self.pid.Mid()

    def Area(self):
        return self.pid.Area()

    def Nsm(self):
        #print "CBEAM pid = ",str(self.pid)
        return self.pid.Nsm()

    def getOfft_Bit_defaults(self):
        if self.isOfft:
            field8 = self.offt
        else:
            field8 = set_blank_if_default(self.bit, 0.0)
        return field8

    def cross_reference(self, model):
        msg = ' which is required by %s eid=%s' % (self.type, self.eid)
        self.ga = model.Node(self.ga, msg=msg)
        self.gb = model.Node(self.gb, msg=msg)
        self.pid = model.Property(self.pid, msg=msg)

    def Stiffness(self, model, r, A, E, I):  # CBEAM
        """
        from makeTruss???
        """
        Ke = zeros((6, 6), 'd')
        L = r
        AE = A * E
        EI = E * I

        if 1:
            Ke[0, 0] = AE / L
            Ke[3, 0] = -AE / L
            Ke[0, 3] = -AE / L
            Ke[3, 3] = AE / L

            Ke[1, 1] = 12 * EI / L ** 3
            Ke[1, 2] = 6 * EI / L ** 2
            Ke[2, 1] = Ke[1, 2]  # 6*EI/L**2
            Ke[2, 2] = 4 * EI / L

            Ke[1, 4] = -Ke[1, 1]  # -12*EI/L**3
            Ke[1, 5] = Ke[1, 2]  # 6*EI/L**2
            Ke[2, 4] = -Ke[1, 2]  # -6*EI/L**2
            Ke[2, 5] = 2 * EI / L

            Ke[4, 1] = -Ke[1, 4]  # -12*EI/L**3
            Ke[4, 2] = Ke[2, 4]  # -6*EI/L**2
            Ke[5, 1] = Ke[1, 2]  # 6*EI/L**2
            Ke[5, 2] = Ke[2, 5]  # 2*EI/L

            Ke[4, 4] = Ke[1, 1]  # 12*EI/L**3
            Ke[4, 5] = Ke[2, 4]  # -6*EI/L**2
            Ke[5, 4] = Ke[2, 4]  # -6*EI/L**2
            Ke[5, 5] = Ke[2, 2]  # 4*EI/L
        else:
            Ke[0, 0] = AE
            Ke[3, 0] = -AE
            Ke[0, 3] = -AE
            Ke[3, 3] = AE

            Ke[1, 1] = 12 * EI / L ** 2
            Ke[1, 2] = 6 * EI / L
            Ke[2, 1] = Ke[1, 2]  # 6*EI/L**2
            Ke[2, 2] = 4 * EI

            Ke[1, 4] = -Ke[1, 1]  # -12*EI/L**3
            Ke[1, 5] = Ke[1, 2]  # 6*EI/L**2
            Ke[2, 4] = -Ke[1, 2]  # -6*EI/L**2
            Ke[2, 5] = 2 * EI

            Ke[4, 1] = -Ke[1, 4]  # -12*EI/L**3
            Ke[4, 2] = Ke[2, 4]  # -6*EI/L**2
            Ke[5, 1] = Ke[1, 2]  # 6*EI/L**2
            Ke[5, 2] = Ke[2, 5]  # 2*EI/L

            Ke[4, 4] = Ke[1, 1]  # 12*EI/L**3
            Ke[4, 5] = Ke[2, 4]  # -6*EI/L**2
            Ke[5, 4] = Ke[2, 4]  # -6*EI/L**2
            Ke[5, 5] = Ke[2, 2]  # 4*EI/L

            Ke = Ke / L
        return Ke

    def rawFields(self):
        (x1, x2, x3) = self.getX_G0_defaults()
        offt = self.getOfft_Bit_defaults()
        ga, gb = self.nodeIDs()
        list_fields = ['CBEAM', self.eid, self.Pid(), ga, gb, x1, x2, x3, offt,
                  self.pa, self.pb, self.w1a, self.w2a, self.w3a,
                  self.w1b, self.w2b, self.w3b, self.sa, self.sb]
        return list_fields

    def reprFields(self):
        w1a = set_blank_if_default(self.w1a, 0.0)
        w2a = set_blank_if_default(self.w2a, 0.0)
        w3a = set_blank_if_default(self.w3a, 0.0)
        w1b = set_blank_if_default(self.w1b, 0.0)
        w2b = set_blank_if_default(self.w2b, 0.0)
        w3b = set_blank_if_default(self.w3b, 0.0)
        (x1, x2, x3) = self.getX_G0_defaults()
        offt = self.getOfft_Bit_defaults()
        ga, gb = self.nodeIDs()
        list_fields = ['CBEAM', self.eid, self.Pid(), ga, gb, x1, x2, x3, offt,
                  self.pa, self.pb, w1a, w2a, w3a,
                  w1b, w2b, w3b, self.sa, self.sb]
        return list_fields


class CBEND(LineElement):
    type = 'CBEND'

    def __init__(self, card=None, data=None, comment=''):
        LineElement.__init__(self, card, data)
        if comment:
            self._comment = comment
        if card:
            self.eid = integer(card, 1, 'eid')
            self.pid = integer_or_blank(card, 2, 'pid', self.eid)
            self.ga = integer(card, 3, 'ga')
            self.gb = integer(card, 4, 'gb')
            x1Go = integer_double_or_blank(card, 5, 'x1_g0')
            if isinstance(x1Go, int):
                self.g0 = x1Go
                self.x1 = None
                self.x2 = None
                self.x3 = None
            elif isinstance(x1Go, float):
                self.g0 = None
                self.x1 = x1Go
                self.x2 = double(card, 6, 'x2')
                self.x3 = double(card, 7, 'x3')
            else:
                raise ValueError('invalid x1Go=|%s| on CBEND' % x1Go)
            self.geom = integer(card, 8, 'geom')
            assert len(card) <= 8, 'len(CBEND card) = %i' % len(card)
            assert self.geom in [1, 2, 3, 4], 'geom is invalid geom=|%s|' % self.geom
        else:
            raise NotImplementedError(data)

    def Area(self):
        return self.pid.Area()

    def rawFields(self):
        (x1, x2, x3) = self.getX_G0_defaults()
        list_fields = ['CBEND', self.eid, self.Pid(), self.Ga(), self.Gb(),
                  x1, x2, x3, self.geom]
        return list_fields

    def reprFields(self):
        return self.rawFields()