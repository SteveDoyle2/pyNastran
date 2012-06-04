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
import sys

#from numpy import zeros,dot
#from numpy.linalg import inv

from pyNastran.bdf.cards.baseCard import Element
from pyNastran.general.generalMath import Area,Triangle_AreaCentroidNormal,Normal
from pyNastran.bdf.errors import *

class ShellElement(Element):
    def __init__(self,card,data):
        Element.__init__(self,card,data)

    def Area(self):
        raise NotImplementedError('Area undefined for %s' %(self.type))

    def Thickness(self):
        #raise Exception('area undefined for %s' %(self.type))
        return self.pid.Thickness()

    def mid(self):
        return self.pid.mid()

    def Mid(self):
        return self.pid.Mid()

    def Rho(self):
        return self.pid.mid().rho

    def Nsm(self):
        return self.pid.Nsm()

    def MassPerArea(self):
        return self.pid.MassPerArea()

    def Mass(self):
        """
        \f[ \large  mass = \frac{mass}{area} area  \f]
        """
        return self.pid.MassPerArea()*self.Area()

    def crossReference(self,mesh):
        self.nodes = mesh.Nodes(self.nodes)
        self.pid   = mesh.Property(self.pid)

class TriShell(ShellElement):
    def __init__(self,card,data):
        ShellElement.__init__(self,card,data)

    def Thickness(self):
        return self.pid.Thickness()

    def AreaCentroidNormal(self):
        """
        returns area,centroid, normal as it's more efficient to do them together
        """
        (n0,n1,n2) = self.nodePositions()
        return Triangle_AreaCentroidNormal([n0,n1,n2])

    def Area(self):
        """
        returns the normal vector
        \f[ \large A = \frac{1}{2} (n_0-n_1) \cross (n_0-n_2)  \f]
        """
        (n0,n1,n2) = self.nodePositions()
        a = n0-n1
        b = n0-n2
        area = Area(a,b)
        return area

    def Normal(self):
        """
        returns the normal vector
        \f[ \large a = (n_0-n_1) \cross (n_0-n_2)  \f]
        \f[ \large n = \frac{n}{norm(N)}           \f]
        """
        (n0,n1,n2) = self.nodePositions()
        a = n0-n1
        b = n0-n2
        return Normal(a,b)

    def Centroid(self,debug=False):
        """
        returns the centroid
        \f[ \large CG = \frac{1}{3} (n_0+n_1+n_2)  \f]
        """
        (n0,n1,n2) = self.nodePositions()
        centroid = self.CentroidTriangle(n0,n1,n2,debug)
        return centroid

    def MassMatrix(self):
        """
        6x6 mass matrix triangle
        http://www.colorado.edu/engineering/cas/courses.d/IFEM.d/IFEM.Ch32.d/IFEM.Ch32.pdf
        """
        mass = self.Mass() # rho*A*t
        if isLumped:  # lumped mass matrix
            Mass = mass/3*eye(6);
        else: # consistent mass
            M = eye(6)*2.
            M[2,0] = M[4,0] = M[0,2] = M[0,4] = 1.
            M[3,1] = M[5,1] = M[1,3] = M[1,5] = 1.
            M[4,2] = M[2,4] = 1.
            M[5,3] = M[5,3] = 1.
            Mass = mass/12*M
        return Mass

class CTRIA3(TriShell):
    type = 'CTRIA3'
    def __init__(self,card=None,data=None):
        TriShell.__init__(self,card,data)
        if card:
            ## element ID number
            self.eid = int(card.field(1))
            self.pid = card.field(2)

            nids = card.fields(3,6)

            self.thetaMcid = card.field(6,0.0)
            self.zOffset   = card.field(7,0.0)

            self.TFlag = card.field(10,0)
            self.T1 = card.field(11,1.0)
            self.T2 = card.field(12,1.0)
            self.T3 = card.field(13,1.0)
        else:
            self.eid = data[0]
            self.pid = data[1]
            nids = data[2:5]

            self.thetaMcid = data[5]
            self.zOffset   = data[6]
            self.TFlag     = data[7]
            self.T1 = data[8]
            self.T2 = data[9]
            self.T3 = data[10]
            if self.T1==-1.0: self.T1=1.0
            if self.T2==-1.0: self.T2=1.0
            if self.T3==-1.0: self.T3=1.0
        ###
        self.prepareNodeIDs(nids)
        assert len(self.nodes)==3

    def Interp(self,un):
        """
        Interpolation based on the area coordinates
        """
        (n0,n1,n2) = self.nodePositions()
        nc = (n0+n1+n2)/3.

        a = n0-nc
        b = n1-nc
        c = n2-nc

        tA1 = det(cross(b,c))   # 2*A1
        tA2 = det(cross(c,a))   # 2*A2
        tA3 = det(cross(a,b))   # 2*A3
        otA = 1./(tA1+tA2+tA3)  # 1/2A

        S = array([tA1,tA2,tA3])*otA  # Ai/A
        u = S*un
        return un

    def Jacob(self):
        (n0,n1,n2) = self.nodePositions()
        (nx0,ny0,nz0) = n0
        (nx1,ny1,nz1) = n1
        (nx2,ny2,nz2) = n2
        #J = matrix([n0,n1-n0,n2-n0])
        
        J = matrix([ [nx0,nx1-nx0,nx2-nx0],
                     [ny0,ny1-ny0,ny2-ny0],
                     [nz0,nz1-nz0,nz2-nz0], ])
        detJ = J.det()
        return J

    def getReprDefaults(self):
        zOffset   = self.setBlankIfDefault(self.zOffset,0.0)
        TFlag     = self.setBlankIfDefault(self.TFlag,0)
        thetaMcid = self.setBlankIfDefault(self.thetaMcid,0.0)

        T1 = self.setBlankIfDefault(self.T1,1.0)
        T2 = self.setBlankIfDefault(self.T2,1.0)
        T3 = self.setBlankIfDefault(self.T3,1.0)
        return (thetaMcid,zOffset,TFlag,T1,T2,T3)

    def rawFields(self):
        fields = ['CTRIA3',self.eid,self.Pid()]+self.nodeIDs()+[self.thetaMcid,self.zOffset,None]+[
            None,self.TFlag,self.T1,self.T2,self.T3]
        return fields

    def reprFields(self):
        (thetaMcid,zOffset,TFlag,T1,T2,T3) = self.getReprDefaults()
        fields = [self.type,self.eid,self.Pid()]+self.nodeIDs()+[thetaMcid,zOffset,None]+[
        None,TFlag,T1,T2,T3]
        return fields

class CTRIA6(TriShell):
    type = 'CTRIA6'
    def __init__(self,card=None,data=None):
        TriShell.__init__(self,card,data)
        if card:
            ## element ID number
            self.eid = int(card.field(1))
            self.pid = card.field(2)

            nids = card.fields(3,9)
            self.thetaMcid = card.field(9 ,0.0)
            self.zOffset   = card.field(10,0.0)

            self.T1 = card.field(11,1.0)
            self.T2 = card.field(12,1.0)
            self.T3 = card.field(13,1.0)
            self.TFlag = card.field(14,0)
        else:
            self.eid = data[0]
            self.pid = data[1]
            nids     = data[2:8]
            self.thetaMcid = data[8]
            self.zOffset   = data[8]
            self.T1    = data[9]
            self.T2    = data[10]
            self.T3    = data[11]
            self.TFlag = data[12]
        self.prepareNodeIDs(nids,allowEmptyNodes=True)
        assert len(nids)==6,'error on CTRIA6'

        #print self.thetaMcid
        #print card

        #print "self.xi = ",self.xi
        #raise

    def Thickness(self):
        return self.pid.Thickness()

    def getReprDefaults(self):
        zOffset   = self.setBlankIfDefault(self.zOffset,0.0)
        TFlag     = self.setBlankIfDefault(self.TFlag,0)
        thetaMcid = self.setBlankIfDefault(self.thetaMcid,0.0)

        T1 = self.setBlankIfDefault(self.T1,1.0)
        T2 = self.setBlankIfDefault(self.T2,1.0)
        T3 = self.setBlankIfDefault(self.T3,1.0)
        return (thetaMcid,zOffset,TFlag,T1,T2,T3)

    def rawFields(self):
        fields = ['CTRIA6',self.eid,self.Pid()]+self.nodeIDs()+[self.thetaMcid,self.zOffset,None]+[
        None,self.TFlag,self.T1,self.T2,self.T3]
        return fields

    def reprFields(self):
        (thetaMcid,zOffset,TFlag,T1,T2,T3) = self.getReprDefaults()
        fields = ['CTRIA6',self.eid,self.Pid()]+self.nodeIDs()+[thetaMcid,zOffset,None]+[
        None,TFlag,T1,T2,T3]
        return fields

class CTRIAR(TriShell):
    type = 'CTRIAR'
    def __init__(self,card=None,data=None):
        TriShell.__init__(self,card,data)
        ## element ID number
        self.eid = int(card.field(1))
        self.pid = card.field(2)

        nids = card.fields(3,6)
        self.prepareNodeIDs(nids)
        assert len(self.nodes)==3

        self.thetaMcid = card.field(6,0.0)
        self.zOffset   = card.field(7,0.0)

        self.TFlag = card.field(10,0)
        self.T1 = card.field(11,1.0)
        self.T2 = card.field(12,1.0)
        self.T3 = card.field(13,1.0)

    def Thickness(self):
        return self.pid.Thickness()

    def getReprDefaults(self):
        zOffset   = self.setBlankIfDefault(self.zOffset,0.0)
        TFlag     = self.setBlankIfDefault(self.TFlag,0)
        thetaMcid = self.setBlankIfDefault(self.thetaMcid,0.0)

        T1 = self.setBlankIfDefault(self.T1,1.0)
        T2 = self.setBlankIfDefault(self.T2,1.0)
        T3 = self.setBlankIfDefault(self.T3,1.0)
        return (thetaMcid,zOffset,TFlag,T1,T2,T3)

    def rawFields(self):
        fields = [self.eid,self.Pid()]+self.nodeIDs()+[self.thetaMcid,self.zOffset,self.TFlag,self.T1,self.T2,self.T3]
        return fields

    def reprFields(self):
        (thetaMcid,zOffset,TFlag,T1,T2,T3) = self.getReprDefaults()
        fields = ['CTRIAR',self.eid,self.Pid()]+self.nodeIDs()+[thetaMcid,zOffset,None,
                  None,TFlag,T1,T2,T3]
        return fields

class CTRIAX(TriShell):
    type = 'CTRIAX'
    def __init__(self,card=None,data=None):
        TriShell.__init__(self,card,data)
        ## element ID number
        self.eid = int(card.field(1))

        nids = card.fields(3,9)
        self.prepareNodeIDs(nids,allowEmptyNodes=True)
        assert len(nids)==6,'error on CTRIAX'

    def rawFields(self):
        fields = ['CTRIAX',self.eid,self.Pid()]+self.nodeIDs()
        return fields

    def reprFields(self):
        return self.rawFields()

class CTRIAX6(TriShell):
    """
    Nodes defined in a non-standard way
         5
        / \
       6   4
     /       \
    1----2----3
    """
    type = 'CTRIAX6'
    def __init__(self,card=None,data=None):
        TriShell.__init__(self,card,data)
        ## element ID number
        self.eid = int(card.field(1))
        self.mid = int(card.field(2))

        nids = card.fields(3,9)
        self.prepareNodeIDs(nids,allowEmptyNodes=True)
        assert len(nids)==6,'error on CTRIAX6'

        self.th = self.setDefaultIfBlank(card.fields(10),0.0)

    def Area(self):
        """
        returns the normal vector
        \f[ \large A = \frac{1}{2} (n_0-n_1) \cross (n_0-n_2)  \f]
        """
        (n1,n2,n3,n4,n5,n6) = self.nodePositions()
        a = n1-n3
        b = n1-n5
        area = Area(a,b)
        return area

    def Thickness(self):
        raise InvalidRequestError('CTRIAX6 does not have a thickness')

    def crossReference(self,model):
        self.nodes = model.Nodes(self.nodes)
        self.mid   = model.Material(self.mid)
    
    def Mid(self):
        if isinstance(self.mid,int):
            return self.mid
        return self.mid.mid

    def rawFields(self):
        fields = ['CTRIAX6',self.eid,self.Mid(),self.Pid()]+self.nodeIDs()+[
                  self.th]
        return fields

    def reprFields(self):
        th = self.th
        fields = ['CTRIAX6',self.eid,self.Mid()]+self.nodeIDs()+[
                  th]
        return fields

class QuadShell(ShellElement):
    def __init__(self,card=None,data=None):
        ShellElement.__init__(self,card,data)

    def Thickness(self):
        return self.pid.Thickness()

    def Normal(self):
        (n1,n2,n3,n4) = self.nodePositions()
        a = n1-n3
        b = n2-n4
        return Normal(a,b)

    def AreaCentroidNormal(self):
        (area,centroid) = self.AreaCentroid()
        normal = self.Normal()
        return (area,centroid,normal)

    def AreaCentroid(self,debug=False):
        """
        1-----2
        |    /|
        | A1/ |
        |  /  |
        |/ A2 |
        4-----3
        
        centroid
           c = sum(ci*Ai)/sum(A)
           where:
             c=centroid
             A=area
        """
        #if debug:
        #    print "nodes = ",self.nodes
        (n1,n2,n3,n4) = self.nodePositions()
        a = n1-n2
        b = n2-n4
        area1 = Area(a,b)
        c1 = self.CentroidTriangle(n1,n2,n4)

        a = n2-n4
        b = n2-n3
        area2 = Area(a,b)
        c2 = self.CentroidTriangle(n2,n3,n4)
        
        area = area1+area2
        centroid = (c1*area1+c2*area2)/area
        if debug:
            print "c1=%s \n c2=%s \n a1=%s a2=%s" %(c1,c2,area1,area2)
            print "type(centroid=%s centroid=%s \n" %(type(centroid),centroid)
        return(area,centroid)
    ###

    def Centroid(self,debug=False):
        nodes = self.nodePositions()
        (area,centroid) = self.AreaCentroid(debug)
        return centroid

    def Area(self):
        """
        \f[ A = \frac{1}{2} \lvert (n_1-n_3) \times (n_2-n_4) \rvert \f]
        where a and b are the quad's cross node point vectors
        """
        (n1,n2,n3,n4) = self.nodePositions()
        a = n1-n3
        b = n2-n4
        area = Area(a,b)
        return area
    ###

    def MassMatrix(self):
        """
        6x6 mass matrix triangle
        http://www.colorado.edu/engineering/cas/courses.d/IFEM.d/IFEM.Ch32.d/IFEM.Ch32.pdf
        """
        mass = self.Mass() # rho*A*t
        if isLumped:  # lumped mass matrix
            Mass = mass/3*eye(6);
        else: # consistent mass
            if G==1:
                M = eye(8)*1. # 1x1 Gauss Rule
                M[2,0] = M[3,1] = M[4,2] = M[5,3] = M[6,4] = M[7,5] = 1.
                M[4,0] = M[5,1] = M[6,2] = M[7,3] = 1.
                M[6,0] = M[7,1] = 1.

                M[0,2] = M[1,3] = M[2,4] = M[3,5] = M[4,6] = M[5,7] = 1.
                M[0,4] = M[1,5] = M[2,6] = M[3,7] = 1.
                M[0,6] = M[1,7] = 1.
                Mass = mass/32*M
            if G==2:
                M = eye(8)*4. # 2x2 Gauss Rule
                M[2,0] = M[3,1] = M[4,2] = M[5,3] = M[6,4] = M[7,5] = 2.
                M[4,0] = M[5,1] = M[6,2] = M[7,3] = 1.
                M[6,0] = M[7,1] = 2.

                M[0,2] = M[1,3] = M[2,4] = M[3,5] = M[4,6] = M[5,7] = 2.
                M[0,4] = M[1,5] = M[2,6] = M[3,7] = 1.
                M[0,6] = M[1,7] = 2.
                Mass = mass/72*M
        return Mass

    def writeAsCTRIA3(self,newID):
        """
        triangle - 012
        triangle - 023
        """
        zOffset = self.setBlankIfDefault(self.zOffset,0.0)
        nodes1 = [self.nodes[0],self.nodes[1],self.nodes[2]]
        nodes2 = [self.nodes[0],self.nodes[2],self.nodes[3]]
        fields1 = ['CTRIA3',self.eid,self.mid1]+nodes1+[self.thetaMcid,zOffset]
        fields2 = ['CTRIA3',newID,   self.mid1]+nodes2+[self.thetaMcid,zOffset]
        return self.printCard(fields1)+printCard(fields2)

    def getReprDefaults(self,debug=False):
        zOffset   = self.setBlankIfDefault(self.zOffset,0.0)
        TFlag     = self.setBlankIfDefault(self.TFlag,0)
        thetaMcid = self.setBlankIfDefault(self.thetaMcid,0.0)

        T1 = self.setBlankIfDefault(self.T1,1.0)
        T2 = self.setBlankIfDefault(self.T2,1.0)
        T3 = self.setBlankIfDefault(self.T3,1.0)
        T4 = self.setBlankIfDefault(self.T4,1.0)
        if debug:
        #if 1:
            print "eid       = ",self.eid
            print "nodes     = ",self.nodes

            print "self.zOffset   = ",self.zOffset
            print "self.TFlag     = ",self.TFlag
            print "self.thetaMcid = ",self.thetaMcid

            print "zOffset   = ",zOffset
            print "TFlag     = ",TFlag
            print "thetaMcid = ",thetaMcid

            print "T1 = ",T1
            print "T2 = ",T2
            print "T3 = ",T3
            print "T4 = ",T4,'\n'
        return (thetaMcid,zOffset,TFlag,T1,T2,T3,T4)

class CSHEAR(QuadShell):
    type = 'CSHEAR'
    def __init__(self,card=None,data=None):
        QuadShell.__init__(self,card,data)
        if card:
            self.eid = card.field(1)
            self.pid = card.field(2)
            nids = card.fields(3,7)
        else:
            self.eid = data[0]
            self.pid = data[1]
            nids = data[2:]
        ###
        self.prepareNodeIDs(nids)
        assert len(self.nodes)==4

    def Normal(self):
        (n1,n2,n3,n4) = self.nodePositions()
        a = n1-n3
        b = n2-n4
        return Normal(a,b)

    def AreaCentroidNormal(self):
        (area,centroid) = self.AreaCentroid()
        normal = self.Normal()
        return (area,centroid,normal)

    def AreaCentroid(self,debug=False):
        """
        1-----2
        |    /|
        | A1/ |
        |  /  |
        |/ A2 |
        4-----3
        
        centroid
           c = sum(ci*Ai)/sum(A)
           where:
             c=centroid
             A=area
        """
        #if debug:
        #    print "nodes = ",self.nodes
        (n1,n2,n3,n4) = self.nodePositions()
        a = n1-n2
        b = n2-n4
        area1 = Area(a,b)
        c1 = self.CentroidTriangle(n1,n2,n4)

        a = n2-n4
        b = n2-n3
        area2 = Area(a,b)
        c2 = self.CentroidTriangle(n2,n3,n4)
        
        area = area1+area2
        centroid = (c1*area1+c2*area2)/area
        if debug:
            print "c1=%s \n c2=%s \n a1=%s a2=%s" %(c1,c2,area1,area2)
            print "type(centroid=%s centroid=%s \n" %(type(centroid),centroid)
        return(area,centroid)
    ###

    def Centroid(self,debug=False):
        (area,centroid) = self.AreaCentroid(debug)
        return centroid

    def Area(self):
        """
        \f[ A = \frac{1}{2} \lvert (n_1-n_3) \times (n_2-n_4) \rvert \f]
        where a and b are the quad's cross node point vectors
        """
        (n1,n2,n3,n4) = self.nodePositions()
        a = n1-n3
        b = n2-n4
        area = Area(a,b)
        return area
    ###

    def rawFields(self):
        fields = ['CSHEAR',self.eid,self.Pid()]+self.nodeIDs()
        return fields
    
    def reprFields(self):
        return self.rawFields()


class CQUAD4(QuadShell):
    type = 'CQUAD4'
    def __init__(self,card=None,data=None):
        QuadShell.__init__(self,card,data)
        if card:
            ## element ID number
            self.eid = int(card.field(1))
            self.pid = card.field(2)

            nids = card.fields(3,7)

            self.thetaMcid = card.field(7,0.0)
            self.zOffset   = card.field(8,0.0)

            self.TFlag = card.field(10,0)
            self.T1 = card.field(11,1.0)
            self.T2 = card.field(12,1.0)
            self.T3 = card.field(13,1.0)
            self.T4 = card.field(14,1.0)
        else:
            self.eid = data[0]
            self.pid = data[1]
            nids = data[2:6]

            self.thetaMcid = data[6]
            self.zOffset   = data[7]
            self.TFlag     = data[8]
            self.T1 = data[9]
            self.T2 = data[10]
            self.T3 = data[11]
            self.T4 = data[12]
            if self.T1==-1.0: self.T1=1.0
            if self.T2==-1.0: self.T2=1.0
            if self.T3==-1.0: self.T3=1.0
            if self.T4==-1.0: self.T4=1.0
        ###            
        self.prepareNodeIDs(nids)
        assert len(self.nodes)==4,'CQUAD4'

        #print "self.xi = ",self.xi
        #print "nodes = ",self.nodes
        #for nid in nids:
        #    self.nodes.append(int(nid))

        #print 'self.T1 = ',self.T1
        #sys.exit()
        #if self.id==20020:
            #print "thetaMcid = ",card.field(7)
            #print "actual = ",self.thetaMcid
            #print str(self)
            #sys.exit()
        
    def rawFields(self):
        fields = [self.type,self.eid,self.Pid()]+self.nodeIDs()+[self.thetaMcid,self.zOffset,self.TFlag,self.T1,self.T2,self.T3,self.T4]
        return fields

    def reprFields(self):
        debug = False
        #if self.id==20020:
        #    debug = True
        (thetaMcid,zOffset,TFlag,T1,T2,T3,T4) = self.getReprDefaults(debug=debug)

        fields = ['CQUAD4',self.eid,self.Pid()]+self.nodeIDs()+[thetaMcid,zOffset,
                  None,TFlag,T1,T2,T3,T4]

        #if self.id==20020:
        #    print "thetaMcid = ",thetaMcid
        #    print "actual = ",self.thetaMcid
        #    #print str(self)
        #    print fields
        #    sys.exit()
        return fields

class CQUADR(QuadShell):
    type = 'CQUADR'
    def __init__(self,card=None,data=None):
        QuadShell.__init__(self,card,data)
        if card:
            ## element ID number
            self.eid = int(card.field(1))
            self.pid = card.field(2)

            nids = card.fields(3,7)

            self.thetaMcid = card.field(7,0.0)
            self.zOffset   = card.field(8,0.0)

            self.TFlag = card.field(10,0)
            self.T1 = card.field(11,1.0)
            self.T2 = card.field(12,1.0)
            self.T3 = card.field(13,1.0)
            self.T4 = card.field(14,1.0)
        else:
            self.eid = data[0]
            self.pid = data[1]
            nids = data[2:6]

            self.thetaMcid = data[6]
            self.zOffset   = data[7]
            self.TFlag     = data[8]
            self.T1 = data[9]
            self.T2 = data[10]
            self.T3 = data[11]
            self.T4 = data[12]
            if self.T1==-1.0: self.T1=1.0
            if self.T2==-1.0: self.T2=1.0
            if self.T3==-1.0: self.T3=1.0
            if self.T4==-1.0: self.T4=1.0
        ###            
        self.prepareNodeIDs(nids)
        assert len(self.nodes)==4,'CQUADR'

    def Thickness(self):
        return self.pid.Thickness()

    def rawFields(self):
        fields = ['CQUADR',self.eid,self.Pid()]+self.nodeIDs()+[self.thetaMcid,self.zOffset,
                  None,self.TFlag,self.T1,self.T2,self.T3,self.T4,]
        return fields

    def reprFields(self):
        return self.rawFields()

class CQUAD(QuadShell):
    type = 'CQUAD'
    def __init__(self,card=None,data=None):
        QuadShell.__init__(self,card,data)
        ## element ID number
        self.eid = int(card.field(1))
        self.pid = card.field(2)

        nids = card.fields(3,12)
        self.prepareNodeIDs(nids)
        assert len(self.nodes)==9

        self.thetaMcid = card.field(7,0.0)
        self.zOffset   = card.field(8,0.0)

        self.TFlag = card.field(10,0)
        self.T1 = card.field(11,1.0)
        self.T2 = card.field(12,1.0)
        self.T3 = card.field(13,1.0)
        self.T4 = card.field(14,1.0)

    def Thickness(self):
        return self.pid.Thickness()

    def rawFields(self):
        fields = ['CQUAD',self.eid,self.Pid()]+self.nodeIDs()
        return fields

    def reprFields(self):
        (thetaMcid,zOffset,TFlag,T1,T2,T3,T4) = self.getReprDefaults()
        fields = ['CQUAD',self.eid,self.Pid()]+self.nodeIDs()
        return fields

class CQUAD8(QuadShell):
    type = 'CQUAD8'
    def __init__(self,card=None,data=None):
        QuadShell.__init__(self,card,data)
        if card:
            ## element ID number
            self.eid = int(card.field(1))
            self.pid = card.field(2)
            nids     = card.fields(3,11)
            self.T1  = card.field(11,1.0)
            self.T2  = card.field(12,1.0)
            self.T3  = card.field(13,1.0)
            self.T4  = card.field(14,1.0)
            self.thetaMcid = card.field(15,0.0)
            self.zOffset   = card.field(16,0.0)
            self.TFlag     = card.field(17,0)
        else:
            #print "CQUAD8 = ",data
            #(6401, 
            #6400, 
            #6401, 6402, 6405, 6403, 0, 0, 6404, 0, 
            #-1.0, -1.0, -1.0, -1.0, 
            #0.0, 0)
            self.eid = data[0]
            self.pid = data[1]
            nids     = data[2:10]
            self.T1  = data[10]
            self.T2  = data[11]
            self.T3  = data[12]
            self.T4  = data[13]
            self.thetaMcid = data[14]
            self.zOffset   = data[14]
            self.TFlag     = data[15]
        ###
        self.prepareNodeIDs(nids,allowEmptyNodes=True)
        assert len(self.nodes)==8

    def Thickness(self):
        return self.pid.Thickness()

    def rawFields(self):
        fields = ['CQUAD8',self.eid,self.Pid()]+self.nodeIDs()+[
                  self.T1,self.T2,self.T3,self.T4,self.thetaMcid,self.zOffset,
                  self.TFlag]
        return fields

    def reprFields(self):
        (thetaMcid,zOffset,TFlag,T1,T2,T3,T4) = self.getReprDefaults()
        fields = ['CQUAD8',self.eid,self.Pid()]+self.nodeIDs()+[
                  T1,T2,T3,T4,thetaMcid,zOffset,
                  TFlag]
        return fields

class CQUADX(QuadShell):
    type = 'CQUADX'
    def __init__(self,card=None,data=None):
        QuadShell.__init__(self,card,data)
        ## element ID number
        self.eid = int(card.field(1))
        self.pid = card.field(2)

        nids = card.fields(3,12)
        self.prepareNodeIDs(nids,allowEmptyNodes=True)
        assert len(self.nodes)==9

    def Thickness(self):
        return self.pid.Thickness()

    def rawFields(self):
        fields = ['CQUADX',self.eid,self.Pid()]+self.nodeIDs()
        return fields
    
    def reprFields(self):
        return self.rawFields()
