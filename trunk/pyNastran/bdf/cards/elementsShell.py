import sys

from baseCard import Element
from pyNastran.general.generalMath import Area

class ShellElement(Element):
    def __init__(self,card,data):
        Element.__init__(self,card,data)

    def Area(self):
        raise Exception('area undefined for %s' %(self.type))

    def Mass(self):
        """
        \f[ \large  mass = \frac{mass}{area} area  \f]
        """
        return self.pid.massPerArea()*self.Area()

    def crossReference(self,mesh):
        self.nodes = mesh.Nodes(self.nodes)
        self.pid   = mesh.Property(self.pid)

class CTRIA3(ShellElement):
    type = 'CTRIA3'
    def __init__(self,card=None,data=None):
        ShellElement.__init__(self,card,data)
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

    def getReprDefaults(self):
        zOffset   = self.setBlankIfDefault(self.zOffset,0.0)
        TFlag     = self.setBlankIfDefault(self.TFlag,0)
        thetaMcid = self.setBlankIfDefault(self.thetaMcid,0.0)

        T1 = self.setBlankIfDefault(self.T1,1.0)
        T2 = self.setBlankIfDefault(self.T2,1.0)
        T3 = self.setBlankIfDefault(self.T3,1.0)
        return (thetaMcid,zOffset,TFlag,T1,T2,T3)

    def __repr__(self):
        (thetaMcid,zOffset,TFlag,T1,T2,T3) = self.getReprDefaults()

        fields = [self.type,self.eid,self.Pid()]+self.nodeIDs()+[thetaMcid,zOffset,None]+[
        None,TFlag,T1,T2,T3]
        return self.printCard(fields)

    def AreaCentroidNormal(self):
        """
        returns area,centroid, normal as it's more efficient to do them together
        """
        (n0,n1,n2) = self.nodePositions()
        return Triangle_AreaCentroidNormal(nodes)

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

class CTRIA6(CTRIA3):
    type = 'CTRIA6'
    def __init__(self,card=None,data=None):
        ShellElement.__init__(self,card,data)
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
        self.prepareNodeIDs(nids)
        assert len(nids)==6,'error on CTRIA6'

        #print self.thetaMcid
        #print card

        #print "self.xi = ",self.xi
        #raise

    def __repr__(self):
        (thetaMcid,zOffset,TFlag,T1,T2,T3) = self.getReprDefaults()
        fields = [self.type,self.eid,self.Pid()]+self.nodeIDs()+[thetaMcid,zOffset,None]+[
        None,TFlag,T1,T2,T3]
        return self.printCard(fields)

class CTRIAR(CTRIA3):
    type = 'CTRIAR'
    def __init__(self,card=None,data=None):
        ShellElement.__init__(self,card,data)
        ## element ID number
        self.eid = int(card.field(1))
        self.pid = card.field(2)

        nids = card.fields(3,6)
        self.prepareNodeIDs(nids)
        assert len(self.nodes)==4

        self.thetaMcid = card.field(6,0.0)
        self.zOffset   = card.field(7,0.0)

        self.TFlag = card.field(10,0)
        self.T1 = card.field(11,1.0)
        self.T2 = card.field(12,1.0)
        self.T3 = card.field(13,1.0)

    def __repr__(self):
        (thetaMcid,zOffset,TFlag,T1,T2,T3) = self.getReprDefaults()
        fields = [self.type,self.eid,self.Pid()]+self.nodeIDs()+[thetaMcid,zOffset,None,
                  None,TFlag,T1,T2,T3]
        return self.printCard(fields)

class CTRIAX(CTRIA3):
    type = 'CTRIAX'
    def __init__(self,card=None,data=None):
        ShellElement.__init__(self,card,data)
        ## element ID number
        self.eid = int(card.field(1))

        nids = card.fields(3,9)
        self.prepareNodeIDs(nids)
        assert len(nids)==6,'error on CTRIAX'

    def __repr__(self):
        fields = [self.type,self.eid,self.Pid()]+self.nodeIDs()
        return self.printCard(fields)

class CTRIAX6(CTRIA3):
    type = 'CTRIAX6'
    def __init__(self,card=None,data=None):
        ShellElement.__init__(self,card,data)
        ## element ID number
        self.eid = int(card.field(1))

        nids = card.fields(3,9)
        self.prepareNodeIDs(nids)
        assert len(nids)==6,'error on CTRIAX6'

        self.th = self.setDefaultIfBlank(card.fields(10),0.0)

    def __repr__(self):
        th = self.th
        fields = [self.type,self.eid,self.Pid()]+self.nodeIDs()+[
                  th]
        return self.printCard(fields)

class CQUAD4(ShellElement):
    type = 'CQUAD4'
    def __init__(self,card=None,data=None):
        ShellElement.__init__(self,card,data)
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
            print "T4 = ",T4
        return (thetaMcid,zOffset,TFlag,T1,T2,T3,T4)

    def __repr__(self):
        debug = False
        #if self.id==20020:
        #    debug = True
        (thetaMcid,zOffset,TFlag,T1,T2,T3,T4) = self.getReprDefaults(debug=debug)


        fields = [self.type,self.eid,self.Pid()]+self.nodeIDs()+[thetaMcid,zOffset,
                  None,TFlag,T1,T2,T3,T4]

        #if self.id==20020:
        #    print "thetaMcid = ",thetaMcid
        #    print "actual = ",self.thetaMcid
        #    #print str(self)
        #    print fields
        #    sys.exit()

        return self.printCard(fields)

    def writeAsCTRIA3(self,newID):
        """
        triangle - 012
        triangle - 023
        """
        zOffset = self.setBlankIfDefault(self.zOffset,0.0)
        nodes1 = [self.nodes[0],self.nodes[1],self.nodes[2]]
        nodes2 = [self.nodes[0],self.nodes[2],self.nodes[3]]
        fields1 = ['CTRIA3',self.eid,self.mid]+nodes1+[self.thetaMcid,zOffset]
        fields2 = ['CTRIA3',newID,   self.mid]+nodes2+[self.thetaMcid,zOffset]
        return self.printCard(fields1)+printCard(fields2)
    
    def Normal(self,nodes):
        (n1,n2,n3,n4) = self.nodePositions()
        a = n1-n3
        b = n2-n4
        return Normal(a,b)

    def AreaCentroidNormal(self):
        (n1,n2,n3,n4) = self.nodePositions()
        (area,centroid) = self.AreaCentroid(nodes)
        normal = self.Normal(nodes)
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

    def Centroid(self,nodes,debug=False):
        (area,centroid) = self.AreaCentroid(nodes,debug)
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

class CQUADR(CQUAD4):
    type = 'CQUADR'
    def __init__(self,card=None,data=None):
        CQUAD4.__init__(self,card,data)

    def __repr__(self):
        fields = [self.type,self.eid,self.Pid()]+self.nodeIDs()+[
                  self.T1,self.T2,self.T3,self.T4,self.thetaMcid,self.zOffset,
                  self.TFlag]
        return self.printCard(fields)

class CQUAD(CQUAD4):
    type = 'CQUAD'
    def __init__(self,card=None,data=None):
        ShellElement.__init__(self,card,data)
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

    def __repr__(self):
        (thetaMcid,zOffset,TFlag,T1,T2,T3,T4) = self.getReprDefaults()
        fields = [self.type,self.eid,self.Pid()]+self.nodeIDs()+[thetaMcid,zOffset,
                  None,TFlag,T1,T2,T3,T4]
        return self.printCard(fields)

class CQUAD8(CQUAD4):
    type = 'CQUAD8'
    def __init__(self,card=None,data=None):
        ShellElement.__init__(self,card,data)
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
            print "CQUAD8 = ",data
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
        self.prepareNodeIDs(nids)
        assert len(self.nodes)==8

    def __repr__(self):
        (thetaMcid,zOffset,TFlag,T1,T2,T3,T4) = self.getReprDefaults()
        fields = [self.type,self.eid,self.Pid()]+self.nodeIDs()+[
                  T1,T2,T3,T4,thetaMcid,zOffset,
                  TFlag]
        return self.printCard(fields)

class CQUADX(CQUAD4):
    type = 'CQUADX'
    def __init__(self,card=None,data=None):
        ShellElement.__init__(self,card,data)
        ## element ID number
        self.eid = int(card.field(1))
        self.pid = card.field(2)

        nids = card.fields(3,12)
        self.prepareNodeIDs(nids)
        assert len(self.nodes)==9

    def __repr__(self):
        fields = [self.type,self.eid,self.Pid()]+self.nodeIDs()
        return self.printCard(fields)
