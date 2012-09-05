# pylint: disable=C0103,R0902,R0904,R0914
from __future__ import (nested_scopes, generators, division, absolute_import,
                        print_function, unicode_literals)
from numpy import dot, cross, array, matrix, zeros
from numpy.linalg import solve

from pyNastran.bdf.cards.elements.elements import Element
from pyNastran.general.mathematics import Area, gauss


def Volume4(n1, n2, n3, n4):
    r"""
    V = (a-d) * ((b-d) x (c-d))/6   where x is cross and * is dot
    \f[ \large V = \frac{(a-d) \cdot \left( (b-d) \times (c-d) \right) }{6} \f]
    """
    V = dot((n1 - n4), cross(n2 - n4, n3 - n4)) / 6.
    return V


class SolidElement(Element):
    def __init__(self, card, data):
        Element.__init__(self, card, data)

    def cross_reference(self, model):
        self.nodes = model.Nodes(self.nodes)
        self.pid = model.Property(self.pid)

    def Mass(self):
        return self.Rho() * self.Volume()

    def Mid(self):
        return self.pid.Mid()

    def Rho(self):
        try:
            return self.pid.mid.rho
        except AttributeError:
            print("self.pid = %s" % (self.pid))
            print("self.pid.mid = %s" % (str(self.pid.mid)))
            #print("self.pid = %s" %(self.pid))
            #print("self.pid = %s" %(self.pid))
            raise

    def isSameCard(self, elem, debug=False):
        if self.type != elem.type:
            return False
        fields1 = [self.eid, self.Pid()] + self.nodes
        fields2 = [elem.eid, elem.Pid()] + elem.nodes
        if debug:
            print("fields1=%s fields2=%s" % (fields1, fields2))
        return self.isSameFields(fields1, fields2)

    def rawFields(self):
        fields = [self.type, self.eid, self.Pid()] + self.nodeIDs()
        return fields


class CHEXA8(SolidElement):
    """
    @code
    CHEXA EID PID G1 G2 G3 G4 G5 G6
    G7 G8
    @endcode
    """
    type = 'CHEXA'
    asterType = 'HEXA8'
    calculixType = 'C3D8'

    def __init__(self, card=None, data=None):
        SolidElement.__init__(self, card, data)
        if card:
            #print "hexa = ",card
            self.eid = card.field(1)
            self.pid = card.field(2)
            nids = card.fields(3, 11)
        else:
            self.eid = data[0]
            self.pid = data[1]
            nids = data[2:]
            assert len(data) == 10, 'len(data)=%s data=%s' % (len(data), data)
        #print "nids = ",nids
        self.prepareNodeIDs(nids)
        assert len(self.nodes) == 8

    def areaCentroid(self, n1, n2, n3, n4):
        a = n1 - n2
        b = n2 - n4
        area1 = Area(a, b)
        c1 = (n1 + n2 + n4) / 3.

        a = n2 - n4
        b = n2 - n3
        area2 = Area(a, b)
        c2 = (n2 + n3 + n4) / 3.

        area = area1 + area2
        centroid = (c1 * area1 + c2 * area2) / area
        return(area, centroid)

    def Centroid(self):
        (n1, n2, n3, n4, n5, n6, n7, n8) = self.nodePositions()
        A1, c1 = self.areaCentroid(n1, n2, n3, n4)
        A2, c2 = self.areaCentroid(n5, n6, n7, n8)
        c = (c1 + c2) / 2.
        return c

    def Volume(self):
        (n1, n2, n3, n4, n5, n6, n7, n8) = self.nodePositions()
        (A1, c1) = self.areaCentroid(n1, n2, n3, n4)
        (A2, c2) = self.areaCentroid(n5, n6, n7, n8)
        V = (A1 + A2) / 2. * (c1 - c2)
        return abs(V)


class CHEXA20(CHEXA8):
    """
    @code
    CHEXA EID PID G1 G2 G3 G4 G5 G6
    G7 G8 G9 G10 G11 G12 G13 G14
    G15 G16 G17 G18 G19 G20
    @endcode
    """
    type = 'CHEXA'
    asterType = 'HEXA20'
    calculixType = 'C3D20'

    def __init__(self, card=None, data=None):
        SolidElement.__init__(self, card, data)

        if card:
            self.eid = card.field(1)
            self.pid = card.field(2)
            nids = card.fields(3, 23)
        else:
            self.eid = data[0]
            self.pid = data[1]
            nids = data[2:]
        self.prepareNodeIDs(nids, allowEmptyNodes=True)
        msg = 'len(nids)=%s nids=%s' % (len(nids), nids)
        assert len(self.nodes) <= 20, msg

    def Centroid(self):
        (n1, n2, n3, n4, n5,
         n6, n7, n8, n9, n10,
         n11, n12, n13, n14, n15,
         n16, n17, n18, n19, n20) = self.nodePositions()
        (A1, c1) = self.areaCentroid(n1, n2, n3, n4)
        (A2, c2) = self.areaCentroid(n5, n6, n7, n8)
        c = (c1 + c2) / 2.
        return c

    def Volume(self):
        (n1, n2, n3, n4, n5,
         n6, n7, n8, n9, n10,
         n11, n12, n13, n14, n15,
         n16, n17, n18, n19, n20) = self.nodePositions()
        (A1, c1) = self.areaCentroid(n1, n2, n3, n4)
        (A2, c2) = self.areaCentroid(n5, n6, n7, n8)
        V = (A1 + A2) / 2. * (c1 - c2)
        return abs(V)


class CPENTA6(SolidElement):
    """
    @code
    CPENTA EID PID G1 G2 G3 G4 G5 G6
      *----------*
     / \        / \
    / A \      / c \
    *---*-----*-----*
    V = (A1+A2)/2  * (c1-c2)
    C = (c1-c2)/2
    @endcode
    """
    type = 'CPENTA'
    asterType = 'PENTA6'
    calculixType = 'C3D6'

    def __init__(self, card=None, data=None):
        SolidElement.__init__(self, card, data)

        if card:
            self.eid = card.field(1)
            self.pid = card.field(2)
            nids = card.fields(3, 9)
        else:
            self.eid = data[0]
            self.pid = data[1]
            nids = data[2:]
            assert len(data) == 8, 'len(data)=%s data=%s' % (len(data), data)
        self.prepareNodeIDs(nids)
        assert len(self.nodes) == 6

    def getFaceNodesAndArea(self, nid, nidOpposite):
        nids = self.nodeIDs()[:6]
        indx1 = nids.index(nid)
        indx2 = nids.index(nidOpposite)
        ##  offset so it's easier to map the nodes with the QRG
        pack = [indx1 + 1, indx2 + 1]
        pack.sort()
        mapper = {  # reverse points away from the element
                    [1, 2]: [1, 2, 3],  # close
                    [2, 3]: [1, 2, 3],
                    [1, 3]: [1, 2, 3],

                    [4, 5]: [4, 5, 6],  # far-reverse
                    [5, 6]: [4, 5, 6],
                    [4, 6]: [4, 5, 6],

                    [1, 5]: [1, 2, 5, 4],  # bottom
                    [2, 4]: [1, 2, 5, 4],

                    [1, 6]: [1, 3, 6, 4],  # left-reverse
                    [3, 4]: [1, 3, 6, 4],

                    [2, 6]: [2, 5, 6, 3],  # right
                    [3, 5]: [2, 3, 6, 5],
        }
        pack2 = mapper[pack]
        if len(pack2) == 3:
            (n1, n2, n3) = pack2
            faceNodeIDs = [n1, n2, n3]
            n1i = nids.index(n1 - 1)
            n2i = nids.index(n2 - 1)
            n3i = nids.index(n3 - 1)
            p1 = self.nodes[n1i].Position()
            p2 = self.nodes[n2i].Position()
            p3 = self.nodes[n3i].Position()
            a = p1 - p2
            b = p1 - p3
            A = Area(a, b)
        else:
            (n1, n2, n3, n4) = pack2
            n1i = nids.index(n1 - 1)
            n2i = nids.index(n2 - 1)
            n3i = nids.index(n3 - 1)
            n4i = nids.index(n4 - 1)
            faceNodeIDs = [n1, n2, n3, n4]
            p1 = self.nodes[n1i].Position()
            p2 = self.nodes[n2i].Position()
            p3 = self.nodes[n3i].Position()
            p4 = self.nodes[n4i].Position()
            a = p1 - p3
            b = p2 - p4
            A = Area(a, b)
        return [faceNodeIDs, A]

    def Centroid(self):
        (n1, n2, n3, n4, n5, n6) = self.nodePositions()
        c1 = (n1 + n2 + n3) / 3.
        c2 = (n4 + n5 + n6) / 3.
        c = (c1 + c2) / 2.
        return c

    def Volume(self):
        (n1, n2, n3, n4, n5, n6) = self.nodePositions()
        A1 = Area(n3 - n1, n2 - n1)
        A2 = Area(n6 - n4, n5 - n4)
        c1 = (n1 + n2 + n3) / 3.
        c2 = (n4 + n5 + n6) / 3.

        V = (A1 + A2) / 2. * (c1 - c2)
        return abs(V)


class CPENTA15(CPENTA6):
    """
    @code
    CPENTA EID PID G1 G2 G3 G4 G5 G6
    G7 G8 G9 G10 G11 G12 G13 G14
    G15
    @endcode
    """
    type = 'CPENTA'
    asterType = 'PENTA15'
    calculixType = 'C3D15'

    def __init__(self, card=None, data=None):
        SolidElement.__init__(self, card, data)

        if card:
            self.eid = card.field(1)
            self.pid = card.field(2)
            nids = card.fields(3, 18)
        else:
            self.eid = data[0]
            self.pid = data[1]
            nids = data[2:]
            assert len(data) == 17, 'len(data)=%s data=%s' % (len(data), data)
        self.prepareNodeIDs(nids, allowEmptyNodes=True)
        assert len(self.nodes) <= 15

    def Centroid(self):
        (n1, n2, n3, n4, n5,
         n6, n7, n8, n9, n10,
         n11, n12, n13, n14, n15) = self.nodePositions()
        c1 = (n1 + n2 + n3) / 3.
        c2 = (n4 + n5 + n6) / 3.
        c = (c1 - c2) / 2.
        return c

    def Volume(self):
        (n1, n2, n3, n4, n5,
         n6, n7, n8, n9, n10,
         n11, n12, n13, n14, n15) = self.nodePositions()
        A1 = Area(n3 - n1, n2 - n1)
        A2 = Area(n6 - n4, n5 - n4)
        c1 = (n1 + n2 + n3) / 3.
        c2 = (n4 + n5 + n6) / 3.

        V = (A1 + A2) / 2. * (c1 - c2)
        return abs(V)


class CTETRA4(SolidElement):
    """
    @code
    CTETRA EID PID G1 G2 G3 G4
    @endcode
    """
    type = 'CTETRA'
    asterType = 'TETRA4'
    calculixType = 'C3D4'

    def __init__(self, card=None, data=None):
        SolidElement.__init__(self, card, data)
        if card:
            self.eid = card.field(1)
            self.pid = card.field(2)
            nids = card.fields(3, 7)
        else:
            self.eid = data[0]
            self.pid = data[1]
            nids = data[2:]
            assert len(data) == 6, 'len(data)=%s data=%s' % (len(data), data)
        self.prepareNodeIDs(nids)
        assert len(self.nodes) == 4

    def Volume(self):
        (n1, n2, n3, n4) = self.nodePositions()
        return Volume4(n1, n2, n3, n4)

    def Centroid(self):
        (n1, n2, n3, n4) = self.nodePositions()
        return (n1 + n2 + n3 + n4) / 4.

    def getFaceNodes(self, nid, nidOpposite):
        nids = self.nodeIDs()[:4]
        indx = nids.index(nidOpposite)
        nids.pop(indx)
        return nids

    def Stiffness(self):
        ng = 1
        (P, W) = gauss(1)

        K = zeros((6, 6))
        for i in xrange(ng):
            for j in xrange(ng):
                for k in xrange(ng):
                    p = [P[i], P[j], P[k]]
                    w = W[i] * W[j] * W[k]
                    pz = self.zeta(p)
                    J = self.Jacobian(p)
                    #N = self.N()
                    #B = self.B()
                    K += w * J * self.BtEB(pz)
                    #K += w*J*B.T*E*B
                ###
            ###
        ###
        return K

    def BtEB(self, pz):
        #E = self.Dsolid()

        #o = E*e
        #e = [exx,eyy,ezz,2exy,2eyz,2ezx]
        #C = array([[Bxi*E11+Byi*E14+Bzi*E16, Byi*E12+Bxi*E14+Bzi*E15, Bzi*E13+Byi*E15+Bxi*E16]
        #           [Bxi*E12+Byi*E24+Bzi*E26, Byi*E22+Bxi*E24+Bzi*E25, Bzi*E23+Byi*E25+Bxi*E26]
        #           [Bxi*E13+Byi*E34+Bzi*E36, Byi*E23+Bxi*E34+Bzi*E35, Bzi*E33+Byi*E35+Bxi*E36]
        #           [Bxi*E14+Byi*E44+Bzi*E46, Byi*E24+Bxi*E44+Bzi*E45, Bzi*E34+Byi*E45+Bxi*E46]
        #           [Bxi*E15+Byi*E45+Bzi*E56, Byi*E25+Bxi*E45+Bzi*E55, Bzi*E35+Byi*E55+Bxi*E56]
        #           [Bxi*E16+Byi*E46+Bzi*E16, Byi*E26+Bxi*E46+Bzi*E56, Bzi*E36+Byi*E56+Bxi*E66]])

        #Qij = array([[Bxj*C11+ Byj*C41+Bzj+Bzj*C61, Bxj*C12+Byj*C42+Bzj*C62, Bxj*C13+Byj*C43+Bzj*C63],
        #             [Byj*C21+ Bxj*C41+Bzj+Bzj*C51, Byj*C22+Bxj*C42+Bzj*C52, Byj*C23+Bxj*C43+Bzj*C53],
        #             [Bzj*C31+ Byj*C51+Bzj+Bxj*C61, Bzj*C32+Byj*C52+Bxj*C62, Bzj*C33+Byj*C53+Bxj*C63]])
        Qij = None
        return Qij

    def zeta(self, p):
        p2 = array([1, p[0], p[1], p[2]])
        J = self.Jacobian()
        zeta = solve(J, p2)
        return zeta

    def Jacobian(self):
        r"""
        \f[ \large   [J] =
          \left[
          \begin{array}{ccc}
              1   & 1   & 1   \\
              x_1 & y_1 & z_1 \\
              x_2 & y_2 & z_2 \\
              x_3 & y_3 & z_3 \\
              x_4 & y_4 & z_4 \\
          \end{array} \right]
        \f]
         @todo
            this has got to be wrong
         @warning
            this has got to be wrong
        """
        m = matrix((6, 6), 'd')
        (n1, n2, n3, n4) = self.nodePositions()
        m[0][0] = m[0][1] = m[0][2] = m[0][2] = 1.
        m[1][0] = n1[0]
        m[2][0] = n1[1]
        m[3][0] = n1[2]
        m[1][1] = n2[0]
        m[2][1] = n2[1]
        m[3][1] = n2[2]
        m[1][2] = n3[0]
        m[2][2] = n3[1]
        m[3][2] = n3[2]
        m[1][3] = n4[0]
        m[2][3] = n4[1]
        m[3][3] = n4[2]
        return m


class CTETRA10(CTETRA4):
    """
    @code
    CTETRA EID PID G1 G2 G3 G4 G5 G6
    G7 G8 G9 G10
    CTETRA   1       1       239     229     516     99      335     103
             265     334     101     102
    @endcode
    """
    type = 'CTETRA'
    asterType = 'TETRA10'
    calculixType = 'C3D10'

    def __init__(self, card=None, data=None):
        SolidElement.__init__(self, card, data)
        if card:
            self.eid = card.field(1)
            self.pid = card.field(2)
            nids = card.fields(3, 13)
        else:
            self.eid = data[0]
            self.pid = data[1]
            nids = data[2:]
            assert len(data) == 12, 'len(data)=%s data=%s' % (len(data), data)
        self.prepareNodeIDs(nids, allowEmptyNodes=True)
        assert len(self.nodes) <= 10

    def N_10(self, g1, g2, g3, g4):
        N1 = g1 * (2 * g1 - 1)
        N2 = g2 * (2 * g2 - 1)
        N3 = g3 * (2 * g3 - 1)
        N4 = g4 * (2 * g4 - 1)
        N5 = 4 * g1 * g2
        N6 = 4 * g2 * g3
        N7 = 4 * g3 * g1
        N8 = 4 * g1 * g4
        N9 = 4 * g2 * g4
        N10 = 4 * g3 * g4
        return (N1, N2, N3, N4, N5, N6, N7, N8, N9, N10)

    def Volume(self):
        (n1, n2, n3, n4, n5, n6, n7, n8, n9, n10) = self.nodePositions()
        return Volume4(n1, n2, n3, n4)

    def Centroid(self):
        (n1, n2, n3, n4, n5, n6, n7, n8, n9, n10) = self.nodePositions()
        return (n1 + n2 + n3 + n4) / 4.

    def getFaceNodes(self, nid, nidOpposite):
        nids = self.nodeIDs()[:4]
        indx = nids.index(nidOpposite)
        nids.pop(indx)
        return nids
