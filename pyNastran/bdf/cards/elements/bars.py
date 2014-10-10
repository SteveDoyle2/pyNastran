# pylint: disable=R0904,R0902,E1101,E1103,C0111,C0302,C0103,W0101
from __future__ import (nested_scopes, generators, division, absolute_import,
                        print_function, unicode_literals)

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


class LineElement(Element):  # CBAR, CBEAM, CBEAM3, CBEND
    def __init__(self, card, data):
        Element.__init__(self, card, data)

    def Lambda(self, model, debug=True):
        (n1, n2) = self.nodeIDs()
        return _Lambda(model, n1, n2, debug=debug)

    def C(self):
        """torsional constant"""
        if isinstance(self.pid, int):
            raise RuntimeError('Element eid=%i has not been cross referenced.\n%s' % (self.eid, str(self)))
        return self.pid.C()

    def Area(self):
        """returns the area of the element face"""
        raise NotImplementedError('implement self.Area() for %s' % self.type)

    def E(self):
        """returns the Young's Modulus, :math:`E`"""
        if isinstance(self.pid, int):
            raise RuntimeError('Element eid=%i has not been cross referenced.\n%s' % (self.eid, str(self)))
        return self.pid.mid.E()

    def G(self):
        """returns the Shear Modulus, :math:`G`"""
        if isinstance(self.pid, int):
            raise RuntimeError('Element eid=%i has not been cross referenced.\n%s' % (self.eid, str(self)))
        return self.pid.mid.G()

    def J(self):
        """returns the Polar Moment of Inertia, :math:`J`"""
        if isinstance(self.pid, int):
            raise RuntimeError('Element eid=%i has not been cross referenced.\n%s' % (self.eid, str(self)))
        return self.pid.J()

    def I11(self):
        """returns the Moment of Inertia, :math:`I_{11}`"""
        if isinstance(self.pid, int):
            raise RuntimeError('Element eid=%i has not been cross referenced.\n%s' % (self.eid, str(self)))
        return self.pid.I11()

    def I22(self):
        """returns the Moment of Inertia, :math:`I_{22}`"""
        if isinstance(self.pid, int):
            raise RuntimeError('Element eid=%i has not been cross referenced.\n%s' % (self.eid, str(self)))
        return self.pid.I22()

    def I12(self):
        """returns the Moment of Inertia, :math:`I_{12}`"""
        if isinstance(self.pid, int):
            raise RuntimeError('Element eid=%i has not been cross referenced.\n%s' % (self.eid, str(self)))
        return self.pid.I12()

    def Nu(self):
        """Get Poisson's Ratio, :math:`\nu`"""
        if isinstance(self.pid, int):
            raise RuntimeError('Element eid=%i has not been cross referenced.\n%s' % (self.eid, str(self)))
        return self.pid.mid.nu

    def Rho(self):
        """Get the material density, :math:`\rho`"""
        #print str(self.pid),type(self.pid)
        #raise NotImplementedError('implement self.Rho() for %s' % self.type)
        if isinstance(self.pid, int):
            raise RuntimeError('Element eid=%i has not been cross referenced.\n%s' % (self.eid, str(self)))
        return self.pid.mid.rho

    def Nsm(self):
        """Placeholder method for the non-structural mass, :math:`nsm`"""
        raise NotImplementedError('implement self.Area() for %s' % self.type)

    def MassPerLength(self):
        """
        Get the mass per unit length, :math:`\frac{m}{L}`
        """
        if isinstance(self.pid, int):
            raise RuntimeError('Element eid=%i has not been cross referenced.\n%s' % (self.eid, str(self)))
        return self.pid.MassPerLength()

    def Mass(self):
        r"""
        Get the mass of the element.

        .. math:: m = \left( \rho A + nsm \right) L
        """
        L = self.Length()
        mass = L * self.MassPerLength()
        #try:
            #mass = (self.Rho() * self.Area() + self.Nsm()) * L
        #except TypeError:
            #msg = 'TypeError on eid=%s pid=%s:\n' % (self.eid, self.Pid())
            #msg += 'rho = %s\narea = %s\nnsm = %s\nL = %s' % (self.Rho(),
            #                                                  self.Area(),
            #                                                  self.Nsm(), L)
            #raise TypeError(msg)

        return mass

    def cross_reference(self, model):
        msg = ' which is required by %s eid=%s' % (self.type, self.eid)
        self.nodes = model.Nodes(self.nodes, msg=msg)
        self.pid = model.Property(self.pid, msg=msg)
        #self.g0 = model.nodes[self.g0]

    def Length(self):
        r"""
        Gets the length, :math:`L`, of the element.

        .. math:: L = \sqrt{  (n_{x2}-n_{x1})^2+(n_{y2}-n_{y1})^2+(n_{z2}-n_{z1})^2  }

        :param self: the object pointer
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
    #     .. warning:: formula not verified
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


class CBAR(LineElement):
    """
    +------+-----+-----+-----+-----+-----+-----+-----+------+
    | CBAR | EID | PID | GA  | GB  | X1  | X2  | X3  | OFFT |
    +------+-----+-----+-----+-----+-----+-----+-----+------+
    | PA   | PB  | W1A | W2A | W3A | W1B | W2B | W3B |
    +------+-----+-----+-----+-----+-----+-----+-----+

    or

    +------+-----+-----+-----+-----+-----+------+
    | CBAR | EID | PID | GA  | GB  | G0  | OFFT |
    +------+-----+-----+-----+-----+-----+------+-----+
    | PA   | PB  | W1A | W2A | W3A | W1B | W2B  | W3B |
    +------+-----+-----+-----+-----+-----+------+-----+

    ::

      CBAR       22062       4   21798   21799   0.0+0   0.0+0     -1.
                                 0.0+0   0.0+0     -9.   0.0+0   0.0+0     -9.
    """
    type = 'CBAR'
    asterType = 'CBAR'
    _field_map = {
        1: 'eid', 2:'pid', 3:'ga', 4:'gb',
        8:'offt', 9:'pa', 10:'pb',
    }

    def _update_field_helper(self, n, value):
        if n == 11:
            self.wa[0] = value
        elif n == 12:
            self.wa[1] = value
        elif n == 13:
            self.wa[2] = value
        elif n == 14:
            self.wb[0] = value
        elif n == 15:
            self.wb[1] = value
        elif n == 16:
            self.wb[2] = value
        else:
            if self.g0 is not None:
                if n == 5:
                    self.g0 = value
                else:
                    raise KeyError('Field %r=%r is an invalid %s entry.' % (n, value, self.type))
            else:
                if n == 5:
                    self.x[0] = value
                elif n == 6:
                    self.x[1] = value
                elif n == 7:
                    self.x[2] = value
                else:
                    raise KeyError('Field %r=%r is an invalid %s entry.' % (n, value, self.type))

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

            self.wa = array([double_or_blank(card, 11, 'w1a', 0.0),
                             double_or_blank(card, 12, 'w2a', 0.0),
                             double_or_blank(card, 13, 'w3a', 0.0)], dtype='float64')

            self.wb = array([double_or_blank(card, 14, 'w1b', 0.0),
                             double_or_blank(card, 15, 'w2b', 0.0),
                             double_or_blank(card, 16, 'w3b', 0.0)], dtype='float64')
            assert len(card) <= 17, 'len(CBAR card) = %i' % len(card)
        else:  #: .. todo:: verify
            #data = [[eid,pid,ga,gb,pa,pb,w1a,w2a,w3a,w1b,w2b,w3b],[f,g0]]
            #data = [[eid,pid,ga,gb,pa,pb,w1a,w2a,w3a,w1b,w2b,w3b],[f,x1,x2,x3]]

            main = data[0]
            flag = data[1][0]
            if flag in [0, 1]:
                self.g0 = None
                self.x = array([data[1][1],
                                data[1][2],
                                data[1][3]], dtype='float64')
            else:
                self.g0 = data[1][1]
                self.x = None

            self.eid = main[0]
            self.pid = main[1]
            self.ga = main[2]
            self.gb = main[3]
            #self.offt = str(data[4]) # GGG
            self.offt = 'GGG'  #: .. todo:: offt can be an integer; translate to char
            self.pa = main[4]
            self.pb = main[5]

            self.wa = array([main[6], main[7], main[8]], dtype='float64')
            self.wb = array([main[9], main[10], main[11]], dtype='float64')
        #print("offt = %s" %(self.offt))

        if not isinstance(self.offt, basestring):
            raise SyntaxError('invalid offt expected a string of length 3 '
                              'offt=|%r|; Type=%s' % (self.offt, type(self.offt)))

        if self.g0 in [self.ga, self.gb]:
            msg = 'G0=%s cannot be GA=%s or GB=%s' % (self.g0, self.ga, self.gb)
            raise RuntimeError(msg)

        msg = 'invalid offt parameter of %s...offt=%s' % (self.type, self.offt)
        # B,G,O
        assert self.offt[0] in ['G', 'B'], msg
        assert self.offt[1] in ['G', 'O', 'E'], msg
        assert self.offt[2] in ['G', 'O', 'E'], msg

    def _verify(self, xref=False):
        pid = self.Pid()
        if xref:  # True
            mid = self.Mid()
            A = self.Area()
            nsm = self.Nsm()
            mpl = self.MassPerLength()
            L = self.Length()
            mass = self.Mass()
        assert isinstance(mid, int), 'mid=%r' % mid
        assert isinstance(A, float), 'A=%r' % A
        assert isinstance(L, float), 'L=%r' % L
        assert isinstance(nsm, float), 'nsm=%r' % nsm
        assert isinstance(mpl, float), 'mass_per_length=%r' % mpl
        assert isinstance(mass, float), 'nass=%r' % mass

    def Mid(self):
        if isinstance(self.pid, int):
            raise RuntimeError('Element eid=%i has not been cross referenced.\n%s' % (self.eid, str(self)))
        return self.pid.Mid()

    def Area(self):
        if isinstance(self.pid, int):
            raise RuntimeError('Element eid=%i has not been cross referenced.\n%s' % (self.eid, str(self)))
        A = self.pid.Area()
        assert isinstance(A, float)
        return A

    def J(self):
        if isinstance(self.pid, int):
            raise RuntimeError('Element eid=%i has not been cross referenced.\n%s' % (self.eid, str(self)))
        j = self.pid.J()
        assert isinstance(j, float), 'J=%r for CBAR eid=%s pid=%s pidType=%s' % (j, self.eid, self.pid.pid, self.pid.type)
        return j

    def Length(self):
        L = norm(self.gb.Position() - self.ga.Position())
        assert isinstance(L, float)
        return L

    def Nsm(self):
        if isinstance(self.pid, int):
            raise RuntimeError('Element eid=%i has not been cross referenced.\n%s' % (self.eid, str(self)))
        nsm = self.pid.Nsm()
        assert isinstance(nsm, float)
        return nsm

    def I1(self):
        if isinstance(self.pid, int):
            raise RuntimeError('Element eid=%i has not been cross referenced.\n%s' % (self.eid, str(self)))
        return self.pid.I1()

    def I2(self):
        return self.pid.I2()

    def Centroid(self):
        return (self.ga.Position() + self.gb.Position()) / 2.

    def initX_G0(self, card):
        field5 = integer_double_or_blank(card, 5, 'g0_x1', 0.0)
        if isinstance(field5, int):
            self.g0 = field5
            self.x = None
        elif isinstance(field5, float):
            self.g0 = None
            self.x = array([field5,
                            double_or_blank(card, 6, 'x2', 0.0),
                            double_or_blank(card, 7, 'x3', 0.0)], dtype='float64')
            if norm(self.x) == 0.0:
                msg = 'G0 vector defining plane 1 is not defined.\n'
                msg += 'G0 = %s\n' % self.g0
                msg += 'X  = %s\n' % self.x
                raise RuntimeError(msg)
        else:
            msg = ('field5 on %s (G0/X1) is the wrong type...id=%s field5=%s '
                   'type=%s' % (self.type, self.eid, field5, type(field5)))
            raise RuntimeError(msg)

    def cross_reference(self, model):
        #if self.g0:
        #    self.x = nodes[self.g0].Position() - nodes[self.ga].Position()
        msg = ' which is required by %s eid=%s' % (self.type, self.eid)
        self.ga = model.Node(self.ga, msg=msg)
        self.gb = model.Node(self.gb, msg=msg)
        self.pid = model.Property(self.pid, msg=msg)

    #def updateNodes(self,nodes):
    #    """
    #    .. todo:: maybe improve
    #    """
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
            #x1 = set_blank_if_default(self.x[0], 0.0)
            #x2 = set_blank_if_default(self.x[1], 0.0)
            #x3 = set_blank_if_default(self.x[2], 0.0)
            return list(self.x)

    def get_orientation_vector(self, xyz):
        if self.g0:
            v = xyz[self.g0] - xyz[self.Ga()]
        else:
            v = self.x
        assert self.offt == 'GGG', self.offt
        return v

    def nodeIDs(self):
        return [self.Ga(), self.Gb()]

    def Stiffness(self, model, grav, is3D=False):  # CBAR
        #print("----------------")
        Lambda = self.Lambda(model, is3D=is3D)
        Lambda_beam = self.Lambda_beam(model)
        #print("Lambda = \n"+str(Lambda))

        #k_beam = self.Stiffness1D(model)
        #k_unit = array([[1.0, -1.0],
        #                [-1.0, 1.0]], dtype='float64')
        #print R
        #print(k)
        #print("Lambda.shape = ", Lambda.shape)
        #K_unit = dot(dot(transpose(Lambda), k_unit), Lambda)

        A = self.Area()
        E = self.E()
        G = self.G()
        J = self.J()
        L = self.L()

        I1 = self.I1()
        I2 = self.I2()

        #k_axial = A * E / L
        #k_torsion = G * J / L
        #k_bend = A * E / L**3

        #K_axial = k_axial * K_unit
        #K_torsion = k_torsion * K_unit
        #K_beam = dot(dot(transpose(Lambda_beam), k_beam), Lambda)

        Lambda_beam2 = Lambda_beam
        K_beam = dot(Lambda_beam2, k_beam)
        K_beam = dot(dot(transpose(Lambda_beam2), k_beam), Lambda_beam2)
        K = K_beam
        #Fg = dot(dot(transpose(Lambda),grav),Lambda)

        #print('size(grav) =',grav.shape)
        mass = self.Mass()
        mg = -grav * mass
        if is3D:
            Fg = [mg[0], mg[1], mg[2], mg[0], mg[1], mg[2]]
        else:
            Fg = [mg[0], mg[1], mg[0], mg[1]]
        #print("Fg = ",Fg)
        #print("mass = ", mass)
        #print("Fg = \n", Fg)
        #print(K)
        #print("K[%s] = \n%s\n" % (self.eid, K/10*6))
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
        #print("A  = ", A)
        #print("E  = ", E)
        #print("L  = ", L)
        #print("I1 = ", I1)
        #print("I2 = ", I2)
        #ki = 1.
        #ki = A*E/L
        #K = ki*matrix([[1.,-1.],[-1.,1.]]) # rod
        #P = 10.416666667
        P = 1.
        k1 = A * E / L
        k2 = P * L ** 3 / (E * I1)
        #print("A=%g E=%g L=%g I1=%g I2=%G AE/L=%g L^3/E*Iz=%g" % (
        #    A, E, L, I1, I2, k1, k2))
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
        #print("k =\n" + str(K))
        return K, None

    def Lambda(self, model, debug=True):  # CBAR
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
            print("v1=%s" % v1)
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
        """
        .. todo:: not perfectly accurate b/c ???
        """
        (x1, x2, x3) = self.getX_G0_defaults()
        offt = set_blank_if_default(self.offt, 'GGG')
        list_fields = ['CBAR', self.eid, self.Pid(), self.Ga(), self.Gb(), x1, x2,
                  x3, offt, self.pa, self.pb] + list(self.wa) + list(self.wb)
        return list_fields

    def reprFields(self):
        pa = set_blank_if_default(self.pa, 0)
        pb = set_blank_if_default(self.pb, 0)

        w1a = set_blank_if_default(self.wa[0], 0.0)
        w2a = set_blank_if_default(self.wa[1], 0.0)
        w3a = set_blank_if_default(self.wa[2], 0.0)

        w1b = set_blank_if_default(self.wb[0], 0.0)
        w2b = set_blank_if_default(self.wb[1], 0.0)
        w3b = set_blank_if_default(self.wb[2], 0.0)
        (x1, x2, x3) = self.getX_G0_defaults()
        offt = set_blank_if_default(self.offt, 'GGG')
        list_fields = ['CBAR', self.eid, self.Pid(), self.Ga(), self.Gb(), x1, x2,
                  x3, offt, pa, pb, w1a, w2a, w3a, w1b, w2b, w3b]
        return list_fields

    def write_bdf(self, size, card_writer):
        card = self.reprFields()
        if size == 8:
            return self.comment() + print_card_8(card)
        return self.comment() + print_card_16(card)
        #return self.comment() + card_writer(card)


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

            self.wa = array([double_or_blank(card, 9, 'w1a', 0.0),
                             double_or_blank(card, 10, 'w2a', 0.0),
                             double_or_blank(card, 11, 'w3a', 0.0)], 'float64')

            self.wb = array([double_or_blank(card, 12, 'w1b', 0.0),
                             double_or_blank(card, 13, 'w2b', 0.0),
                             double_or_blank(card, 14, 'w3b', 0.0)], 'float64')

            self.wc = array([double_or_blank(card, 15, 'w1c', 0.0),
                             double_or_blank(card, 16, 'w2c', 0.0),
                             double_or_blank(card, 17, 'w3c', 0.0)], 'float64')

            self.twa = double_or_blank(card, 18, 0., 'twa')
            self.twb = double_or_blank(card, 19, 0., 'twb')
            self.twc = double_or_blank(card, 20, 0., 'twc')

            self.sa = integer_or_blank(card, 21, 'sa')
            self.sb = integer_or_blank(card, 22, 'sb')
            self.sc = integer_or_blank(card, 23, 'sc')
            assert len(card) <= 24, 'len(CBEAM3 card) = %i' % len(card)
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
        .. math:: L = g_b - g_a
        """
        L = norm(self.gb.Position() - self.ga.Position())
        return L

    def rawFields(self):
        (x1, x2, x3) = self.getX_G0_defaults()
        (ga, gb, gc) = self.nodeIDs()
        list_fields = ['CBEAM3', self.eid, self.Pid(), ga, gb, gc, x1, x2, x3] + \
                  list(self.wa) + list(self.wb) + list(self.wc) + [self.twa,
                  self.twb, self.twc, self.sa, self.sb, self.sc]
        return list_fields

    def reprFields(self):
        w1a = set_blank_if_default(self.wa[0], 0.0)
        w2a = set_blank_if_default(self.wa[1], 0.0)
        w3a = set_blank_if_default(self.wa[2], 0.0)
        w1b = set_blank_if_default(self.wb[0], 0.0)
        w2b = set_blank_if_default(self.wb[1], 0.0)
        w3b = set_blank_if_default(self.wb[2], 0.0)
        w1c = set_blank_if_default(self.wc[0], 0.0)
        w2c = set_blank_if_default(self.wc[1], 0.0)
        w3c = set_blank_if_default(self.wc[2], 0.0)

        twa = set_blank_if_default(self.twa, 0.0)
        twb = set_blank_if_default(self.twb, 0.0)
        twc = set_blank_if_default(self.twc, 0.0)

        (x1, x2, x3) = self.getX_G0_defaults()
        (ga, gb, gc) = self.nodeIDs()
        list_fields = ['CBEAM3', self.eid, self.Pid(), ga, gb, x1, x2, x3,
                  w1a, w2a, w3a, w1b, w2b, w3b, w1c, w2c, w3c,
                  twa, twb, twc, self.sa, self.sb, self.sc]
        return list_fields

    def write_bdf(self, size, card_writer):
        card = self.reprFields()
        return self.comment() + card_writer(card)


class CBEAM(CBAR):
    """
    ::

      CBEAM EID PID GA GB X1 X2 X3 OFFT/BIT
      PA PB W1A W2A W3A W1B W2B W3B
      SA SB

    or::

      CBEAM EID PID GA GB G0 - - OFFT/BIT
      PA PB W1A W2A W3A W1B W2B W3B
      SA SB
    """
    type = 'CBEAM'
    _field_map = {
        1: 'eid', 2:'pid', 3:'ga', 4:'gb', #5:'x_g0', 6:'g1', 7:'g2',
        #8:'offt',
        9:'pa', 10:'pb',
        17:'sa', 18:'sb',
    }

    def _update_field_helper(self, n, value):
        if n == 11:
            self.wa[0] = value
        elif n == 12:
            self.wa[1] = value
        elif n == 13:
            self.wa[2] = value

        elif n == 14:
            self.wb[0] = value
        elif n == 15:
            self.wb[1] = value
        elif n == 16:
            self.wb[2] = value
        else:
            if self.g0 is not None:
                if n == 5:
                    self.g0 = value
                else:  # offt
                    raise KeyError('Field %r=%r is an invalid %s entry or is unsupported.' % (n, value, self.type))
            else:
                if n == 5:
                    self.x[0] = value
                elif n == 6:
                    self.x[1] = value
                elif n == 7:
                    self.x[2] = value
                else:
                    raise KeyError('Field %r=%r is an invalid %s entry or is unsupported.' % (n, value, self.type))

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

            self.wa = array([double_or_blank(card, 11, 'w1a', 0.0),
                             double_or_blank(card, 12, 'w2a', 0.0),
                             double_or_blank(card, 13, 'w3a', 0.0)], 'float64')

            self.wb = array([double_or_blank(card, 14, 'w1b', 0.0),
                             double_or_blank(card, 15, 'w2b', 0.0),
                             double_or_blank(card, 16, 'w3b', 0.0)], 'float64')

            self.sa = integer_or_blank(card, 17, 'sa')
            self.sb = integer_or_blank(card, 18, 'sb')
            assert len(card) <= 19, 'len(CBEAM card) = %i' % len(card)
        else:  #: .. todo:: verify
            #data = [[eid,pid,ga,gb,sa,sb, pa,pb,w1a,w2a,w3a,w1b,w2b,w3b],
            #        [f,g0]]
            #data = [[eid,pid,ga,gb,sa,sb, pa,pb,w1a,w2a,w3a,w1b,w2b,w3b],
            #        [f,x1,x2,x3]]

            main = data[0]

            flag = data[1][0]
            if flag in [0, 1]:
                self.g0 = None
                self.x = array([data[1][1],
                                data[1][2],
                                data[1][3]], dtype='float64')
            else:
                self.g0 = data[1][1]
                self.x = None

            self.eid = main[0]
            self.pid = main[1]
            self.ga = main[2]
            self.gb = main[3]
            self.sa = main[4]
            self.sb = main[5]

            self.isOfft = True  #: .. todo:: is this correct???
            #self.offt = str(data[6]) # GGG
            self.offt = 'GGG'  #: .. todo:: is this correct???

            self.pa = main[6]
            self.pb = main[7]

            self.w1a = main[8]
            self.w2a = main[9]
            self.w3a = main[10]

            self.w1b = main[11]
            self.w2b = main[12]
            self.w3b = main[13]

        if self.g0 in [self.ga, self.gb]:
            msg = 'G0=%s cannot be GA=%s or GB=%s' % (self.g0, self.ga, self.gb)
            raise RuntimeError(msg)

    def Nodes(self):
        return [self.ga, self.gb]

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
        elif isinstance(field8, basestring):
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
        if isinstance(self.pid, int):
            raise RuntimeError('Element eid=%i has not been cross referenced.\n%s' % (self.eid, str(self)))
        return self.pid.Mid()

    def Area(self):
        if isinstance(self.pid, int):
            raise RuntimeError('Element eid=%i has not been cross referenced.\n%s' % (self.eid, str(self)))
        return self.pid.Area()

    def Nsm(self):
        if isinstance(self.pid, int):
            raise RuntimeError('Element eid=%i has not been cross referenced.\n%s' % (self.eid, str(self)))
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
        if self.g0:
            self.g0_vector = model.nodes[self.g0].Position() - self.ga.Position()
        else:
            self.g0_vector = self.x

    def Lambda(self, model, is3D=False, debug=True):  # CBAR from CROD/CONROD
        """
        2d  [l,m,0,0,  l,m,0,0]
            [0,0,l,m,  0,0,l,m] L*k = 2x4*W4x4

        3d  [l,m,n,0,0,0,  l,m,n,0,0,0]
            [0,0,0,l,m,n,  0,0,0,l,m,n]
        """
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

    def displacement_stress(self, model, q, dofs, is3D=False):  # CBEAM
        (n1, n2) = self.nodeIDs()
        Lambda = self.Lambda_beam2(model, n1, n2, debug=False)

        #print("**dofs =", dofs)
        n11 = dofs[(n1, 1)]
        n12 = dofs[(n1, 2)]

        n21 = dofs[(n2, 1)]
        n22 = dofs[(n2, 2)]

        n13 = dofs[(n1, 3)]
        n23 = dofs[(n2, 3)]

        n14 = dofs[(n1, 4)]
        n24 = dofs[(n2, 4)]

        n15 = dofs[(n1, 5)]
        n25 = dofs[(n2, 5)]

        n16 = dofs[(n1, 6)]
        n26 = dofs[(n2, 6)]

        q2 = array([
            q[n11], q[n12], q[n13],
            q[n21], q[n22], q[n23],

            q[n14], q[n15], q[n16],
            q[n24], q[n25], q[n26],
        ])

        #print("type=%s n1=%s n2=%s" % (self.type, n1, n2))
        #print("n11=%s n12=%s n21=%s n22=%s" %(n11,n12,n21,n22))

        #print("q2[%s] = %s" % (self.eid, q2))
        #print("Lambda = \n"+str(Lambda))

        #print "Lsize = ",Lambda.shape
        #print "qsize = ",q.shape
        L = self.Length()
        A = self.Area()
        E = self.E()
        Iz = self.I1()
        Iy = self.I2()
        I1_I2_I12 = self.pid.I1_I2_I12()

        G = self.G()
        J = self.J()
        C = 1.0

        #K = self._stiffness(L, A, E, Iz, Iy, G, J)
        #Kf = dot(transpose(Lambda), dot(K, Lambda))
        u = dot(Lambda, q2)

        du_axial   = -u[0] + u[6]
        du_lateral = -u[1] + u[7]
        du_normal  = -u[2] + u[8]
        du_torsion = -u[3] + u[9]
        du_bend_y  = -u[4] + u[10]
        du_bend_z  = -u[5] + u[11]
        #==============================
        # strain
        axial_strain = du_axial / L
        torsional_strain = du_torsion * C / L

        # stress
        axial_stress = E * axial_strain
        torsional_stress = G * torsional_strain

        # forces
        axial_force = axial_stress * A
        torsional_moment = du_torsion * G * J / L


        #print("axial_strain = %s [psi]" % axial_strain)
        #print("axial_stress = %s [psi]" % axial_stress)
        #print("axial_force  = %s [lb]\n" % axial_force)
        return axial_strain, axial_stress, axial_force

    def Rmatrix(self, model, n1, n2, debug=True):
        p1 = model.Node(n1).Position()
        p2 = model.Node(n2).Position()
        v1 = p2 - p1
        v1 = v1 / norm(v1)

        #v2 = self.g0_vector
        #v2 = p3 - p1

        v3 = cross(v1, self.g0_vector)
        v3 /= norm(v3)

        v2 = cross(v1, v3)

        R = array([ v1, v2, v3 ])
        if debug:
            print("v1 =", v1)
            print("v2 =", v2)
            print("v3 =", v3)
            print('R =\n', R)
            #print('R.shape =', R.shape)
        return R, v1, v2, v3
        (l1, m1, n1) = v1
        (l2, m2, n2) = v2
        (l3, m3, n3) = v3

    def Lambda_beam(self, model, n1, n2, debug=True):
        """
        ::

          2d  [l,m,0,0]
              [0,0,l,m]

        ::

          3d  [l,m,n,0,0,0]  2x6
              [0,0,0,l,m,n]
        """
        R = self.Rmatrix(model, n1, n2)

        p1 = model.Node(n1).Position()
        p2 = model.Node(n2).Position()
        v1 = p2 - p1
        if debug:
            print("v1=%s" % (v1))
        v1 = v1 / norm(v1)
        (l, m, n) = v1
        #l = 1
        #m = 2
        #n = 3
        Lambda = matrix(zeros((4, 12), 'd'))  # 3D
        Lambda[0, 0] = Lambda[1, 3] = l
        Lambda[0, 1] = Lambda[1, 4] = m
        Lambda[0, 2] = Lambda[1, 5] = n

        Lambda[2, 6] = Lambda[3, 9] = l
        Lambda[2, 7] = Lambda[3, 10] = m
        Lambda[2, 8] = Lambda[3, 11] = n

        #print("R = \n",R)
        #print("is3D = \n",is3D)
        #debug = True
        if debug:
            print("Lambda = \n" + str(Lambda))
            #sys.exit('asdf')
        return Lambda

    def Lambda_beam2(self, model, n1, n2, debug=True):
        """
        ::

          2d  [l,m,0,0]
              [0,0,l,m]

        ::

          3d  [l,m,n,0,0,0]  2x6
              [0,0,0,l,m,n]
        """
        R, v1, v2, v3 = self.Rmatrix(model, n1, n2)

        Lambda = array(zeros((12, 12), 'd'))  # 3D
        Lambda[0,0:3] = v1
        Lambda[1,0:3] = v2
        Lambda[2,0:3] = v3

        Lambda[3,3:6] = v1
        Lambda[4,3:6] = v2
        Lambda[5,3:6] = v3

        Lambda[6,6:9] = v1
        Lambda[7,6:9] = v2
        Lambda[8,6:9] = v3

        Lambda[9, 9:] = v1
        Lambda[10,9:] = v2
        Lambda[11,9:] = v3

        if debug:
            print("Lambda = \n" + str(Lambda))
            #sys.exit('asdf')
        return Lambda

    def Fg(self, model, grav, fnorm):
        (n1, n2) = self.nodeIDs()
        m = self.Mass()
        mg = m * grav
        Fg = [mg[0], mg[1], mg[2], mg[0], mg[1], mg[2]]

        nGrav = [(n1, 1), (n1, 2), (n2, 3),
                 (n2, 1), (n2, 2), (n1, 3)]
        return(Fg, nGrav)

    def Stiffness(self, model, node_ids, index0s, fnorm=1.0):  # CBEAM
        """
        from makeTruss???
        http://www.engr.sjsu.edu/ragarwal/ME273/pdf/Chapter%204%20-%20Beam%20Element.pdf
        """
        L = self.Length()
        A = self.Area()
        E = self.E()
        I11 = self.I11()
        I22 = self.I22()
        G = self.G()
        J = self.J()
        print("A=%s E=%s L=%s G=%s J=%s" % (A, E, L, G, J))
        #G = J = A = E = 1.0
        #A = J = G = 0.0

        (n1, n2) = self.nodeIDs()

        #R
        #Lambda = matrix(zeros((4, 12), 'd'))  # 3D
        #Lambda[0, 0] = Lambda[1, 3] = l
        #Lambda[0, 1] = Lambda[1, 4] = m
        #Lambda[0, 2] = Lambda[1, 5] = n

        #Lambda[2, 6] = Lambda[3, 9] = l
        #Lambda[2, 7] = Lambda[3, 10] = m
        #Lambda[2, 8] = Lambda[3, 11] = n

        Ke = self._stiffness(L, A, E, I11, I22, G, J)

        #Lambda = self.Lambda(model)
        Lambda_beam2 = self.Lambda_beam2(model, n1, n2)
        #print("Lambda = \n"+str(Lambda))

        #k_rod = matrix([[1., -1.], [-1., 1.]])  # 1D rod

        #k, ki = self.Stiffness1D(model)  # CBEAM
        #print R
        #print(k)
        #print("Lambda.shape = ", Lambda.shape)
        #K = dot(dot(transpose(Lambda), k), Lambda)
        K = dot(dot(transpose(Lambda_beam2), Ke), Lambda_beam2)
        #
        #= dot(Lambda_beam2, Ke)

        #k_axial = A * E / L
        #k_torsion = G * J / L
        #k_bendiug = k_axial / L**2

        i0, i1 = index0s
        dofs = array([ i0, i0+1, i0+2, i0+3, i0+4, i0+5,
                       i1, i1+1, i1+2, i1+3, i1+4, i1+5], 'int32')
        nIJV = [
            (n1, 1), (n1, 2), (n1, 3), (n1, 4), (n1, 5), (n1, 6),
            (n2, 1), (n2, 2), (n2, 3), (n2, 4), (n2, 5), (n2, 6),
        ]
        return(K, dofs, nIJV)

    def _stiffness(self, L, A, E, Iz, Iy, G, J):  # CBEAM
        AE = A * E
        #EI = E * I

        AE_L = AE / L
        #EI = E * I
        L2 = L * L
        L3 = L2 * L
        #EIz_L3 = EIz / L3
        K = zeros((12, 12), 'd')
        # top-left diagonal going right
        K[0,0] =  A * E / L
        K[1,1] = 12 * E * Iz / L3
        K[2,2] = 12 * E * Iy / L3
        K[3,3] = G * J / L
        K[4,4] = 4 * E * Iy / L
        K[5,5] = 4 * E * Iz / L
        K[6,6] = K[0,0] # A * E / L
        K[7,7] = K[1,1] # 12 * E * Iz / L3
        K[8,8] = K[2,2] # 12 * E * Iy / L3
        K[9,9] = K[3,3] # G * J / L
        K[10,10]= K[4,4]# 4 * E * Iy / L
        K[11,11]= K[5,5]# 4 * E * Iz / L

        # top-center diagonal going right
        #K[0,6] = -A * E / L
        K[1,7]  = -12 * E * Iz / L3
        K[2,8]  = -12 * E * Iy / L3
        K[3,9]  = -G * J / L
        K[4,10] = 2 * E * Iy / L
        K[5,11] = 2 * E * Iz / L

        # top-center diagonal going left
        K[0,6] = -K[0,0] # -A * E / L
        K[1,5] =  6 * E * Iz / L2
        K[2,4] = -6 * E * Iy / L2
        #K[3,3]
        K[4,2] = -6 * E * Iy / L2
        K[5,1] =  6 * E * Iz / L2
        K[6,0] = K[0,6] # -A * E / L

        # top-right diagnoal going left
        K[1,11] =  6 * E * Iz / L2
        K[2,10] = -6 * E * Iy / L2
        #K[3,9] = - G * J / L
        K[4,8] =  6 * E * Iy / L2
        K[5,7] = -6 * E * Iz / L2
        #K[6,6] = A * E / L
        K[7,5] = -6 * E * Iz / L2
        K[8,4] =  6 * E * Iy / L2
        #K[9,3] = -G * J / L
        K[10,2] = -6 *E * Iy / L2
        K[11,1] = 6 * E * Iz / L2

        # center-left diagonal going right
        #K[6,0] = -A * E / L
        K[7,1] = -12 * E * Iz / L3
        K[8,2] = -12 * E * Iy / L3
        K[9,3] = -G * J / L
        K[10,4] = 2 * E * Iy / L
        K[11,5] = 2 * E * Iz / L

        # center-right diagnoal going left
        K[7,11] = -6 * E * Iz / L2
        K[8,10] =  6 * E * Iy / L2
        #K[9, 9]= G * J / L
        K[10,8] =  6 * E * Iy / L2
        K[11,7] = -6 * E * Iz / L2
        return K

    def rawFields(self):
        (x1, x2, x3) = self.getX_G0_defaults()
        offt = self.getOfft_Bit_defaults()
        ga, gb = self.nodeIDs()
        list_fields = ['CBEAM', self.eid, self.Pid(), ga, gb, x1, x2, x3, offt,
                  self.pa, self.pb] + list(self.wa) + list(self.wb) + [self.sa, self.sb]
        return list_fields

    def reprFields(self):
        w1a = set_blank_if_default(self.wa[0], 0.0)
        w2a = set_blank_if_default(self.wa[1], 0.0)
        w3a = set_blank_if_default(self.wa[2], 0.0)
        w1b = set_blank_if_default(self.wb[0], 0.0)
        w2b = set_blank_if_default(self.wb[1], 0.0)
        w3b = set_blank_if_default(self.wb[2], 0.0)
        (x1, x2, x3) = self.getX_G0_defaults()
        offt = self.getOfft_Bit_defaults()
        ga, gb = self.nodeIDs()
        list_fields = ['CBEAM', self.eid, self.Pid(), ga, gb, x1, x2, x3, offt,
                  self.pa, self.pb, w1a, w2a, w3a,
                  w1b, w2b, w3b, self.sa, self.sb]
        return list_fields

    def write_bdf(self, size, card_writer):
        card = self.reprFields()
        if size == 8:
            return self.comment() + print_card_8(card)
        return self.comment() + print_card_16(card)
        #return self.comment() + card_writer(card)


class CBEND(LineElement):
    type = 'CBEND'
    _field_map = {
        1: 'eid', 2:'pid', 3:'ga', 4:'gb', 8:'geom',
    }

    def _update_field_helper(self, n, value):
        if self.g0 is not None:
            if n == 5:
                self.g0 = value
            else:
                raise KeyError('Field %r=%r is an invalid %s entry.' % (n, value, self.type))
        else:
            if n == 5:
                self.x[0] = value
            elif n == 6:
                self.x[1] = value
            elif n == 7:
                self.x[2] = value
            else:
                raise KeyError('Field %r=%r is an invalid %s entry.' % (n, value, self.type))

    def __init__(self, card=None, data=None, comment=''):
        LineElement.__init__(self, card, data)
        if comment:
            self._comment = comment
        if card:
            self.eid = integer(card, 1, 'eid')
            self.pid = integer_or_blank(card, 2, 'pid', self.eid)
            self.ga = integer(card, 3, 'ga')
            self.gb = integer(card, 4, 'gb')
            x1Go = integer_double_or_blank(card, 5, 'x1_g0', 0.0)
            if isinstance(x1Go, int):
                self.g0 = x1Go
                self.x = None
            elif isinstance(x1Go, float):
                self.g0 = None
                self.x = array([double_or_blank(card, 5, 'x1', 0.0),
                                double_or_blank(card, 6, 'x2', 0.0),
                                double_or_blank(card, 7, 'x3', 0.0)], dtype='float64')
                if norm(self.x) == 0.0:
                    msg = 'G0 vector defining plane 1 is not defined.\n'
                    msg += 'G0 = %s\n' % self.g0
                    msg += 'X  = %s\n' % self.x
                    raise RuntimeError(msg)
            else:
                raise ValueError('invalid x1Go=|%s| on CBEND' % x1Go)
            self.geom = integer(card, 8, 'geom')

            assert len(card) == 9, 'len(CBEND card) = %i' % len(card)
            assert self.geom in [1, 2, 3, 4], 'geom is invalid geom=|%s|' % self.geom
        else:
            raise NotImplementedError(data)

        if self.g0 in [self.ga, self.gb]:
            msg = 'G0=%s cannot be GA=%s or GB=%s' % (self.g0, self.ga, self.gb)
            raise RuntimeError(msg)

    def Area(self):
        if isinstance(self.pid, int):
            raise RuntimeError('Element eid=%i has not been cross referenced.\n%s' % (self.eid, str(self)))
        return self.pid.Area()

    def rawFields(self):
        (x1, x2, x3) = self.getX_G0_defaults()
        list_fields = ['CBEND', self.eid, self.Pid(), self.Ga(), self.Gb(),
                  x1, x2, x3, self.geom]
        return list_fields

    def reprFields(self):
        return self.rawFields()

    def write_bdf(self, size, card_writer):
        card = self.reprFields()
        return self.comment() + card_writer(card)
