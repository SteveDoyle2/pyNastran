from numpy import zeros,array

from pyNastran.bdf.cards.baseCard import Element


class PointElement(Element):
    def __init__(self,card=None,data=None):
        Element.__init__(self,card,data)
        
    def Mass(self):
        return self.mass

class CMASS1(PointElement):
    """
    Defines a scalar mass element.
    CMASS2 EID M   G1 C1 G2 C2
    CMASS1 EID PID G1 C1 G2 C2
    """
    type = 'CMASS2'
    def __init__(self,card=None,data=None):
        PointElement.__init__(self,card,data)
        if card:
            self.eid = card.field(1)
            self.pid = card.field(2,self.eid)
            self.g1 = card.field(3)
            self.c1 = card.field(4)
            self.g2 = card.field(5)
            self.c2 = card.field(6)
        else:
            self.eid = data[0]
            self.pid = data[1]
            self.g1  = data[2]
            self.c1  = data[3]
            self.g2  = data[4]
            self.c2  = data[5]
        ###

    def Mass(self):
        return self.pid.mass

    def crossReference(self,mesh):
        #self.g1 = mesh.Node(self.g1)
        #self.g2 = mesh.Node(self.g2)
        self.pid = mesh.Property(self.pid)

    def rawFields(self):
        fields = ['CMASS1',self.eid,self.Pid(),self.g1,self.c1,self.g2,self.c2]
        return fields

class CMASS2(PointElement):
    """
    Defines a scalar mass element without reference to a property entry.
    CMASS2 EID M G1 C1 G2 C2
    """
    type = 'CMASS2'
    def __init__(self,card=None,data=None):
        PointElement.__init__(self,card,data)
        if card:
            self.eid  = card.field(1)
            self.mass = card.field(2,0.)
            self.g1   = card.field(3)
            self.c1   = card.field(4)
            self.g2   = card.field(5)
            self.c2   = card.field(6)
        else:
            self.eid  = data[0]
            self.mass = data[1]
            self.g1   = data[2]
            self.c1   = data[3]
            self.g2   = data[4]
            self.c2   = data[5]
        ###

    def Mass(self):
        return self.mass

    def crossReference(self,mesh):
        #self.g1 = mesh.Node(self.g1)
        #self.g2 = mesh.Node(self.g2)
        pass

    def rawFields(self):
        fields = ['CMASS2',self.eid,self.mass,self.g1,self.c1,self.g2,self.c2]
        return fields

    def reprFields(self):
        mass = self.setBlankIfDefault(self.mass,0.)
        fields = ['CMASS2',self.eid,mass,self.g1,self.c1,self.g2,self.c2]
        #print "cmass2 fields = ",fields
        return fields

class CMASS3(PointElement):
    """
    Defines a scalar mass element that is connected only to scalar points.
    CMASS3 EID PID S1 S2
    """
    type = 'CMASS3'
    def __init__(self,card=None,data=None):
        PointElement.__init__(self,card,data)

        if card:
            self.eid  = card.field(1)
            self.pid  = card.field(2,self.eid)
            self.s1   = card.field(3)
            self.s2   = card.field(4)
        else:
            self.eid = data[0]
            self.pid = data[1]
            self.s1  = data[2]
            self.s2  = data[3]
        ###

    def Mass(self):
        return self.pid.mass

    def isSameCard(self,elem):
        if self.type!=elem.type:  return False
        fields1 = [self.eid,self.Pid(),self.s1,self.s2]
        fields2 = [elem.eid,elem.Pid(),elem.s1,elem.s2]
        return self.isSameFields(fields1,fields2)

    def crossReference(self,mesh):
        """
        links up the propertiy ID
        """
        #self.s1 = mesh.Node(self.s1)
        #self.s2 = mesh.Node(self.s2)
        self.pid = mesh.Property(self.pid)

    def rawFields(self):
        fields = ['CMASS3',self.eid,self.Pid(),self.s1,self.s2]
        return fields

class CMASS4(PointElement):
    """
    Defines a scalar mass element that is connected only to scalar points, without reference to a property entry
    CMASS4 EID M S1 S2
    """
    type = 'CMASS4'
    def __init__(self,card=None,data=None):
        PointElement.__init__(self,card,data)
        
        if card:
            self.eid  = card.field(1)
            self.mass = card.field(2,0.)
            self.s1 = card.field(3)
            self.s2 = card.field(4)
        else:
            self.eid  = data[0]
            self.mass = data[1]
            self.s1   = data[2]
            self.s2   = data[3]
        ###

    def Mass(self):
        return self.mass

    def isSameCard(self,elem):
        if self.type!=elem.type:  return False
        fields1 = [self.eid,self.Pid(),self.s1,self.s2]
        fields2 = [elem.eid,elem.Pid(),elem.s1,elem.s2]
        return self.isSameFields(fields1,fields2)

    def crossReference(self,mesh):
        """
        """
        #self.s1 = mesh.Node(self.s1)
        #self.s2 = mesh.Node(self.s2)
        pass

    def rawFields(self):
        fields = ['CMASS4',self.eid,self.mass,self.s1,self.s2]
        return fields

   
class CONM1(PointElement):
    type = 'CONM1'
    def __init__(self,card=None,data=None):
        """
        Concentrated Mass Element Connection, General Form
        Defines a 6 x 6 symmetric mass matrix at a geometric grid point
        
        CONM1 EID G CID M11 M21 M22 M31 M32
              M33 M41 M42 M43 M44 M51 M52 M53
              M54 M55 M61 M62 M63 M64 M65 M66

        [M] = [M11 M21 M31 M41 M51 M61]
              [    M22 M32 M42 M52 M62]
              [        M33 M43 M53 M63]
              [            M44 M54 M64]
              [    Sym         M55 M65]
              [                    M66]
        """
        PointElement.__init__(self,card,data)
        if card:
            #self.nids  = [ card[1] ]
            #del self.nids
            #self.pid = None
            self.eid = card.field(1)
            self.nid = card.field(2)
            self.cid = card.field(3,0)
            
            m = zeros((6,6))
            m[0,0] = card.field(4,0.)   # M11
            m[1,0] = card.field(5,0.)   # M21
            m[1,1] = card.field(6,0.)   # M22
            m[2,0] = card.field(7,0.)   # M31
            m[2,1] = card.field(8,0.)   # M32
            m[2,2] = card.field(9,0.)   # M33
            m[3,0] = card.field(10,0.)  # M41
            m[3,1] = card.field(11,0.)  # M42
            m[3,2] = card.field(12,0.)  # M43
            m[3,3] = card.field(13,0.)  # M44
            m[4,0] = card.field(14,0.)  # M51
            m[4,1] = card.field(15,0.)  # M52
            m[4,2] = card.field(16,0.)  # M53
            m[4,3] = card.field(17,0.)  # M54
            m[4,4] = card.field(18,0.)  # M55
            m[5,0] = card.field(19,0.)  # M61
            m[5,1] = card.field(20,0.)  # M62
            m[5,2] = card.field(21,0.)  # M63
            m[5,3] = card.field(22,0.)  # M64
            m[5,4] = card.field(23,0.)  # M65
            m[5,5] = card.field(24,0.)  # M66
        ###
        else:
            raise NotImplementedError('COMN1 data')
        ###
        self.massMatrix = m

    def Nid(self):
        if isinstance(self.nid,int):
            return self.nid
        return self.nid.nid

    def Cid(self):
        if isinstance(self.cid,int):
            return self.cid
        return self.cid.cid

    def crossReference(self,model):
        self.nid = model.Node(self.nid)
        self.cid = model.Coord(self.cid)

    def MassMatrix(self):
        return self.massMatrix

    def rawFields(self):
        cid  = self.setBlankIfDefault(self.Cid(),0)
        nid = self.Nid()
        m = self.massMatrix
        fields = ['CONM1',eid,nid,cid,m[0,0],m[1,0],m[1,1],m[2,0],m[2,1],
                 m[2,2],m[3,0],m[3,1],m[3,2],m[3,3],m[4,0],m[4,1],m[4,2],
                 m[4,3],m[4,4],m[5,0],m[5,1],m[5,2],m[5,3],m[5,4],m[5,5] ]
        return fields

    def reprFields(self):
        fields2 = fields[0:4]
        for field in fields[4:]:
            val = self.setBlankIfDefault(field,0.)
            fields2.append(val)
        return fields2

class CONM2(PointElement): ## @todo not done
    """
    @todo support cid != 0
    """
    type = 'CONM2'
    # 'CONM2    501274  11064          132.274'
    def __init__(self,card=None,data=None):
        PointElement.__init__(self,card,data)
        if card:
            #self.nids  = [ card[1] ]
            #del self.nids
            #self.pid = None
            self.eid  = card.field(1)
            self.nid  = card.field(2)
            self.cid  = card.field(3,0)
            self.mass = card.field(4,0.)
            self.X    = array(card.fields(5,8,[0.,0.,0.]))
            self.I    = card.fields(9,15,[0.]*6)
        else:
            self.eid  = data[0]
            self.nid  = data[1]
            self.cid  = data[2]
            self.mass = data[3]
            self.X    = data[4:7]
            self.I    = data[7:]
        ###

    def Mass(self):
        if self.cid==0:
            raise NotImplementedError('cid=0 is not supported for CONM2...')
        return self.mass

    def Centroid(self):
        if self.cid==0:
            raise NotImplementedError('cid=0 is not supported for CONM2...')
        if self.cid==-1:
            return self.X
        else:
            dx,matrix = self.cid.transformToGlobal(self.X)
            #print "dx = ",dx
            #print "X  = ",self.nid.Position()
            X2  = self.nid.Position()+dx
        return X2

    def crossReference(self,model):
        """
        @warning only supports cid=0
        """
        if self.cid!=-1:
            self.cid = model.Coord(self.cid)
        self.nid = model.Node(self.nid)

    def Nid(self):
        if isinstance(self.nid,int):
            return self.nid
        return self.nid.nid

    def Cid(self):
        if isinstance(self.cid,int):
            return self.cid
        return self.cid.cid

    def writeCodeAster(self):
        msg += "    DISCRET=_F(\n"
        msg += "             'CARA='M_T_D_N'\n"
        msg += "              NOEUD=N%s\n" %(self.Nid())
        msg += "              VALE=%g),\n" %(self.mass)
        return msg

    def rawFields(self):
        #print "X=%r" %(self.X)
        #print "I=%r" %(self.I)
        fields = ['CONM2',self.eid,self.Nid(),self.Cid(),self.mass]+list(self.X)+[None]+list(self.I)
        return fields

    def reprFields(self):
        I = []
        for i in self.I:
            if i==0.:
                I.append(None)
            else:
                I.append(i)
            ###
        ###
        X = []
        for x in self.X:
            if x==0.:
                X.append(None)
            else:
                X.append(x)
            ###
        ###

        cid  = self.setBlankIfDefault(self.Cid(),0)
        mass = self.setBlankIfDefault(self.mass,0.)
        fields = ['CONM2',self.eid,self.Nid(),cid,mass]+X+[None]+I
        return fields

