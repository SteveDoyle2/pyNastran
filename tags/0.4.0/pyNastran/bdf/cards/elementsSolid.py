## GNU Lesser General Public License
## 
## Program pyNastran - a python interface to NASTRAN files
## Copyright (C) 2011-2012  Steven Doyle, Al Danial
## 
## Authors and copyright holders of pyNastran
## Steven Doyle <mesheb82@gmail.com>
## Al Danial    <al.danial@gmail.com>
## 
## This file is part of pyNastran.
## 
## pyNastran is free software: you can redistribute it and/or modify
## it under the terms of the GNU Lesser General Public License as published by
## the Free Software Foundation, either version 3 of the License, or
## (at your option) any later version.
## 
## pyNastran is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
## 
## You should have received a copy of the GNU Lesser General Public License
## along with pyNastran.  If not, see <http://www.gnu.org/licenses/>.
## 
from numpy import dot,cross,matrix,sqrt
from numpy.linalg import solve

# pyNastran
from elements import Element
from pyNastran.general.generalMath import Area

def Volume4(n1,n2,n3,n4):
    """
    V = (a-d) * ((b-d) x (c-d))/6   where x is cross and * is dot
    \f[ \large V = {(a-d) \dot \left( (b-d) \times (c-d) \right) }{6} \f]
    """
    V = dot((n1-n4),cross(n2-n4,n3-n4))/6.
    return V

class SolidElement(Element):
    def __init__(self,card,data):
        Element.__init__(self,card,data)

    def crossReference(self,mesh):
        self.nodes = mesh.Nodes(self.nodes)
        self.pid   = mesh.Property(self.pid)

    def Mass(self):
        return self.Rho()*self.Volume()
    
    def Mid(self):
        return self.pid.Mid()

    def Rho(self):
        return self.pid.mid.rho

    def isSameCard(self,elem):
        if self.type!=elem.type:  return False
        fields1 = [self.eid,self.Pid()]+self.nodes
        fields2 = [elem.eid,elem.Pid()]+elem.nodes
        return self.isSameFields(fields1,fields2)

    def rawFields(self):
        fields = [self.type,self.eid,self.Pid()]+self.nodeIDs()
        return fields

class CHEXA8(SolidElement):
    """
    CHEXA EID PID G1 G2 G3 G4 G5 G6
    G7 G8
    """
    type = 'CHEXA'
    def __init__(self,card=None,data=None):
        SolidElement.__init__(self,card,data)
        if card:
            #print "hexa = ",card
            self.eid = card.field(1)
            self.pid = card.field(2)
            nids = card.fields(3,11)
        else:
            self.eid = data[0]
            self.pid = data[1]
            nids     = data[2:]
            assert len(data)==10,'len(data)=%s data=%s' %(len(data),data)
        #print "nids = ",nids
        self.prepareNodeIDs(nids)
        assert len(self.nodes)==8

    def areaCentroid(self,n1,n2,n3,n4):
        a = n1-n2
        b = n2-n4
        area1 = Area(a,b)
        c1 = (n1+n2+n4)/3.

        a = n2-n4
        b = n2-n3
        area2 = Area(a,b)
        c2 = (n2+n3+n4)/3.
        
        area = area1+area2
        centroid = (c1*area1+c2*area2)/area
        return(area,centroid)

    def Centroid(self):
        (n1,n2,n3,n4,n5,n6,n7,n8) = self.nodePositions()
        A1,c1 = self.areaCentroid(n1,n2,n3,n4)
        A2,c2 = self.areaCentroid(n5,n6,n7,n8)
        c = (c1+c2)/2.
        return c

    def Volume(self):
        (n1,n2,n3,n4,n5,n6,n7,n8) = self.nodePositions()
        A1,c1 = self.areaCentroid(n1,n2,n3,n4)
        A2,c2 = self.areaCentroid(n5,n6,n7,n8)
        V = (A1+A2)/2.*(c1-c2)
        return abs(V)

class CHEXA20(CHEXA8):
    """
    CHEXA EID PID G1 G2 G3 G4 G5 G6
    G7 G8 G9 G10 G11 G12 G13 G14
    G15 G16 G17 G18 G19 G20
    """
    type = 'CHEXA'
    def __init__(self,card=None,data=None):
        SolidElement.__init__(self,card,data)

        if card:
            self.eid = card.field(1)
            self.pid = card.field(2)
            nids = card.fields(3,23)
        else:
            self.eid = data[0]
            self.pid = data[1]
            nids     = data[2:]
        self.prepareNodeIDs(nids,allowEmptyNodes=True)
        msg = 'len(nids)=%s nids=%s' %(len(nids),nids)
        assert len(self.nodes)<=20,msg

    def Centroid(self):
        (n1,n2,n3,n4,n5,n6,n7,n8,n9,n10,n11,n12,n13,n14,n15,n16,n17,n18,n19,n20) = self.nodePositions()
        A1,c1 = self.areaCentroid(n1,n2,n3,n4)
        A2,c2 = self.areaCentroid(n5,n6,n7,n8)
        c = (c1+c2)/2.
        return c

    def Volume(self):
        (n1,n2,n3,n4,n5,n6,n7,n8,n9,n10,n11,n12,n13,n14,n15,n16,n17,n18,n19,n20) = self.nodePositions()
        A1,c1 = self.areaCentroid(n1,n2,n3,n4)
        A2,c2 = self.areaCentroid(n5,n6,n7,n8)
        V = (A1+A2)/2.*(c1-c2)
        return abs(V)


class CPENTA6(SolidElement):
    """
    CPENTA EID PID G1 G2 G3 G4 G5 G6
      *----------*
     / \        / \
    / A \      / c \
    *---*-----*-----*
    V = (A1+A2)/2  * (c1-c2)
    C = (c1-c2)/2
    """
    type = 'CPENTA'
    def __init__(self,card=None,data=None):
        SolidElement.__init__(self,card,data)

        if card:
            self.eid = card.field(1)
            self.pid = card.field(2)
            nids = card.fields(3,9)
        else:
            self.eid = data[0]
            self.pid = data[1]
            nids     = data[2:]
            assert len(data)==8,'len(data)=%s data=%s' %(len(data),data)
        self.prepareNodeIDs(nids)
        assert len(self.nodes)==6

    def Centroid(self):
        n = self.nodePositions()
        (n1,n2,n3,n4,n5,n6) = n
        c1 = (n1+n2+n3)/3.
        c2 = (n4+n5+n6)/3.
        c = (c1+c2)/2.
        return c

    def Volume(self):
        n = self.nodePositions()
        (n1,n2,n3,n4,n5,n6) = n
        A1 = Area(n3-n1,n2-n1)
        A2 = Area(n6-n4,n5-n4)
        c1 = (n1+n2+n3)/3.
        c2 = (n4+n5+n6)/3.
       
        V = (A1+A2)/2. * (c1-c2)
        return abs(V)

class CPENTA15(CPENTA6):
    """
    CPENTA EID PID G1 G2 G3 G4 G5 G6
    G7 G8 G9 G10 G11 G12 G13 G14
    G15
    """
    type = 'CPENTA'
    def __init__(self,card=None,data=None):
        SolidElement.__init__(self,card,data)

        if card:
            self.eid = card.field(1)
            self.pid = card.field(2)
            nids = card.fields(3,18)
        else:
            self.eid = data[0]
            self.pid = data[1]
            nids = data[2:]
            assert len(data)==17,'len(data)=%s data=%s' %(len(data),data)
        self.prepareNodeIDs(nids,allowEmptyNodes=True)
        assert len(self.nodes)<=15

    def Centroid(self):
        n = self.nodePositions()
        (n1,n2,n3,n4,n5,n6,n7,n8,n9,n10,n11,n12,n13,n14,n15) = n
        c1 = (n1+n2+n3)/3.
        c2 = (n4+n5+n6)/3.
        c = (c1-c2)/2.
        return c

    def Volume(self):
        n = self.nodePositions()
        (n1,n2,n3,n4,n5,n6,n7,n8,n9,n10,n11,n12,n13,n14,n15) = n
        A1 = Area(n3-n1,n2-n1)
        A2 = Area(n6-n4,n5-n4)
        c1 = (n1+n2+n3)/3.
        c2 = (n4+n5+n6)/3.
       
        V = (A1+A2)/2. * (c1-c2)
        return abs(V)

class CTETRA4(SolidElement):
    """
    CTETRA EID PID G1 G2 G3 G4
    """
    type = 'CTETRA'
    def __init__(self,card=None,data=None):
        SolidElement.__init__(self,card,data)
        if card:
            self.eid = card.field(1)
            self.pid = card.field(2)
            nids = card.fields(3,7)
        else:
            self.eid = data[0]
            self.pid = data[1]
            nids = data[2:]
            assert len(data)==6,'len(data)=%s data=%s' %(len(data),data)
        self.prepareNodeIDs(nids)
        assert len(self.nodes)==4

    def Volume(self):
        (n1,n2,n3,n4) = self.nodePositions()
        return Volume4(n1,n2,n3,n4)

    def Centroid(self):
        (n1,n2,n3,n4) = self.nodePositions()
        return (n1+n2+n3+n4)/4.

    def Stiffness(self):
        ng = 1
        (P,W) = gauss(1)

        K = zeros((6,6))
        for i in range(ng):
            for j in range(ng):
                for k in range(ng):
                    p = [ P[i],P[j],P[k] ]
                    w = W[i]*W[j]*W[k]
                    pz = zeta(p)
                    J = self.Jacobian(p)
                    #N = self.N()
                    #B = self.B()
                    K += w*J*BtEB(pz)
                    #K += w*J*B.T*E*B
                ###
            ###
        ###
        return K

    def BtEB(self,pz):
        E = self.Dsolid()

        #o = E*e
        #e = [exx,eyy,ezz,2exy,2eyz,2ezx]
        C = array([[Bxi*E11+Byi*E14+Bzi*E16, Byi*E12+Bxi*E14+Bzi*E15, Bzi*E13+Byi*E15+Bxi*E16]
                   [Bxi*E12+Byi*E24+Bzi*E26, Byi*E22+Bxi*E24+Bzi*E25, Bzi*E23+Byi*E25+Bxi*E26]
                   [Bxi*E13+Byi*E34+Bzi*E36, Byi*E23+Bxi*E34+Bzi*E35, Bzi*E33+Byi*E35+Bxi*E36]
                   [Bxi*E14+Byi*E44+Bzi*E46, Byi*E24+Bxi*E44+Bzi*E45, Bzi*E34+Byi*E45+Bxi*E46]
                   [Bxi*E15+Byi*E45+Bzi*E56, Byi*E25+Bxi*E45+Bzi*E55, Bzi*E35+Byi*E55+Bxi*E56]
                   [Bxi*E16+Byi*E46+Bzi*E16, Byi*E26+Bxi*E46+Bzi*E56, Bzi*E36+Byi*E56+Bxi*E66]])

        Qij = array([[Bxj*C11+ Byj*C41+Bzj+Bzj*C61, Bxj*C12+Byj*C42+Bzj*C62, Bxj*C13+Byj*C43+Bzj*C63],
                     [Byj*C21+ Bxj*C41+Bzj+Bzj*C51, Byj*C22+Bxj*C42+Bzj*C52, Byj*C23+Bxj*C43+Bzj*C53],
                     [Bzj*C31+ Byj*C51+Bzj+Bxj*C61, Bzj*C32+Byj*C52+Bxj*C62, Bzj*C33+Byj*C53+Bxj*C63]])
        return Qij
            
    def zeta(self,p):
        p2 = array([1,p[0],p[1],p[2]])
        zeta = solve(J,p2)
        return zeta

    def Jacobian(self):
        """
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
        m = matrix((6,6),'d')
        m[0][0] = m[0][1] = m[0][2] = m[0][2] = 1.
        m[1][0]=n1[0]; m[2][0]=n1[1]; m[3][0]=n1[2];
        m[1][1]=n2[0]; m[2][1]=n2[1]; m[3][1]=n2[2];
        m[1][2]=n3[0]; m[2][2]=n3[1]; m[3][2]=n3[2];
        m[1][3]=n4[0]; m[2][3]=n4[1]; m[3][3]=n4[2];
        return m

class CTETRA10(CTETRA4):
    """
    CTETRA EID PID G1 G2 G3 G4 G5 G6
    G7 G8 G9 G10
    CTETRA   1       1       239     229     516     99      335     103
             265     334     101     102
    """
    type = 'CTETRA'
    def __init__(self,card=None,data=None):
        SolidElement.__init__(self,card,data)
        if card:
            self.eid = card.field(1)
            self.pid = card.field(2)
            nids = card.fields(3,13)
        else:
            self.eid = data[0]
            self.pid = data[1]
            nids = data[2:]
            assert len(data)==12,'len(data)=%s data=%s' %(len(data),data)
        self.prepareNodeIDs(nids,allowEmptyNodes=True)
        assert len(self.nodes)<=10

    def N_10(self,g1,g2,g3,g4):
        N1 = g1*(2*g1-1)
        N2 = g2*(2*g2-1)
        N3 = g3*(2*g3-1)
        N4 = g4*(2*g4-1)
        N5 = 4*g1*g2
        N6 = 4*g2*g3
        N7 = 4*g3*g1
        N8 = 4*g1*g4
        N9 = 4*g2*g4
        N10 = 4*g3*g4

    def Volume(self):
        (n1,n2,n3,n4,n5,n6,n7,n8,n9,n10) = self.nodePositions()
        return Volume4(n1,n2,n3,n4)

    def Centroid(self):
        (n1,n2,n3,n4,n5,n6,n7,n8,n9,n10) = self.nodePositions()
        return (n1+n2+n3+n4)/4.
