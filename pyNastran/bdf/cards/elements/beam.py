# pylint: disable=R0904,R0902,E1101,E1103,C0111,C0302,C0103,W0101
from six import string_types
from numpy import matrix, zeros, array, transpose, dot, ones
from numpy import eye, allclose, cross
from numpy.linalg import norm


from .bars import CBAR, LineElement
from pyNastran.bdf.bdfInterface.assign_type import (integer, integer_or_blank,
    integer_double_or_blank, double, double_or_blank,
    string_or_blank, integer_double_string_or_blank, integer_or_double)
from pyNastran.bdf.fieldWriter import set_blank_if_default
from pyNastran.bdf.fieldWriter import print_card_8
from pyNastran.bdf.fieldWriter16 import print_card_16

class CBEAM(CBAR):
    """
    +-------+-----+-----+-----+-----+-----+-----+-----+----------+
    | CBEAM | EID | PID | GA  | GB  | X1  | X2  | X3  | OFFT/BIT |
    +-------+-----+-----+-----+-----+-----+-----+-----+----------+
    |       | PA  | PB  | W1A | W2A | W3A | W1B | W2B | W3B      |
    +-------+-----+-----+-----+-----+-----+-----+-----+----------+
    |       | SA  | SB  |
    +-------+-----+-----+

    or

    +-------+-----+-----+-----+-----+-----+-----+-----+----------+
    | CBEAM | EID | PID | GA  | GB  | G0  |     |     | OFFT/BIT |
    +-------+-----+-----+-----+-----+-----+-----+-----+----------+
    |       | PA  | PB  | W1A | W2A | W3A | W1B | W2B | W3B      |
    +-------+-----+-----+-----+-----+-----+-----+-----+----------+
    |       | SA  | SB  |
    +-------+-----+-----+

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

            self.sa = integer_or_blank(card, 17, 'sa', 0)
            self.sb = integer_or_blank(card, 18, 'sb', 0)
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
        elif isinstance(field8, string_types):
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

        #print("Lsize = ",Lambda.shape)
        #print("qsize = ",q.shape)
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
        #print(R)
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

        sa = set_blank_if_default(self.sa, 0)
        sb = set_blank_if_default(self.sb, 0)
        (x1, x2, x3) = self.getX_G0_defaults()
        offt = self.getOfft_Bit_defaults()
        ga, gb = self.nodeIDs()
        list_fields = ['CBEAM', self.eid, self.Pid(), ga, gb, x1, x2, x3, offt,
                  self.pa, self.pb, w1a, w2a, w3a,
                  w1b, w2b, w3b, sa, sb]
        return list_fields

    def write_bdf(self, size, card_writer):
        card = self.reprFields()
        if size == 8:
            return self.comment() + print_card_8(card)
        return self.comment() + print_card_16(card)
        #return self.comment() + card_writer(card)
