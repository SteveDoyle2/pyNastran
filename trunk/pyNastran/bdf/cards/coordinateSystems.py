#import sys
from numpy import array,cross,dot
from numpy.linalg import norm

# my code
from pyNastran.bdf.BDF_Card import BDF_Card
from baseCard import BaseCard

class Coord(BaseCard):
    def __init__(self,card,data):
        self.isCrossReferenced = False
        self.isResolved = False

        #self.type = card[0]
        #print "cid = ",self.cid

    def setup(self):
        assert len(self.eo)==3,self.eo
        assert len(self.ez)==3,self.ez
        assert len(self.ex)==3,self.ex
        
        #print "eo = ",self.eo
        #print "ez = ",self.ez
        #print "ex = ",self.ex

        ## the normalized version of self.ez
        self.ez0 = self.normalize(self.ez-self.eo)
        ## the normalized version of self.ex
        self.ex0 = self.normalize(self.ex-self.eo)
        ## the normalized version of self.ey
        self.ey0 = cross(self.ez0,self.ex0)

    def normalize(self,v):
        #print "v = ",v
        return v/norm(v)

    def __repr__(self):
        fields = [self.type,self.cid]
        return self.printCard(fields)

class Cord2x(Coord):
    def __init__(self,card,data):
        Coord.__init__(self,card,data)
    
    def Rid(self):
        """Returns the reference coordinate system self.rid"""
        if isinstance(self.rid,int):
            return self.rid
        return self.rid.cid

    def resolveCid(self):
        """
        Turns the coordinate system from being a coordinate system of
        type 1 depending on a type 2 to a type 1 depending on a type 1.

        More generally, takes a coordinate system that depends on multiple
        levels of coordinate systems and resolves them in order to resolve
        it's coordinate frame.  So, if a coordinate system is of type 2, this
        will effectively set rid to 0 with a type 2.
        
        This should handle any number of coordinate systems or coordinate
        system types assuming there is no circular references.
        """
        #print str(self)
        #pritn self.rid
        if self.cid==0 or isinstance(self.rid,int) or self.rid.isResolved:  # rid=0
            return
        elif self.rid.isResolved==False: # rid
            assert self.rid.isCrossReferenced==False,'there is a circular reference between Coord %s and Coord %s' %(self.cid,self.Rid())
            self.rid.resolveCid()
        
        ## rid coordinate system is now resolved, time to resolve the cid coordinate system
        ## rid may be in a different coordinate system than 
        self.eo = self.transformToGlobal(self.eo)
        self.isResolved = True
        
        ## the axes are normalized, so we just assume they're points and resolve them
        self.eo  = self.rid.transformToGlobal(self.eo )
        self.ex0 = self.rid.transformToGlobal(self.ex0)
        self.ey0 = self.rid.transformToGlobal(self.ey0)
        self.ez0 = self.rid.transformToGlobal(self.ez0)
        
        self._transformGlobalToSelf(self)


class Cord1x(Coord):
    def __init__(self,card,data):
        Coord.__init__(self,card,data)

    def resolveCid(self):
        ## the origin
        self.eo = self.g1.Position()
        ## a point on the z-axis
        self.ez = self.g2.Position()
        ## a point on the xz-plane
        self.ex = self.g3.Position()

        self.setup()

class CORD3G(Coord):
    """
    Defines a general coordinate system using three rotational angles as functions of
    coordinate values in the reference coordinate system. The CORD3G entry is used with
    the MAT9 entry to orient material principal axes for 3-D composite analysis

    CORD3G CID METHOD FORM THETAID1 THETAID2 THETAID3 CIDREF
    CORD3G 100 E313   EQN  110      111      112      0
    """
    type = 'CORD3G'
    def __init__(self,card=['CORD3G',0,0,0,0,0,0,0]):
        """
        Intilizes the CORD3G
        @param self   the object pointer
        @param card   a list version of the fields
        """
        if isinstance(card,list):
            assert len(card)==8
            card = BDF_Card(card)
        Coord.__init__(self,card)

        self.cid = card.field(1)
        method   = card.field(2)
        self.methodES  = method[0]
        self.methodInt = int(method[1:])
        assert self.ESmethod in ['E','S']
        assert 0 < self.methodInt < 1000
        
        self.form   = card.field(3,'EQN')
        self.theta1 = card.field(4)
        self.theta2 = card.field(5)
        self.theta3 = card.field(6)
        self.cidRef = card.field(7)

        assert self.form in ['EQN','TABLE'] # EQN for DEQATN, TABLE for TABLE3D
    
    def crossReference(self):
        pass

    def transformToGlobal(self,p,debug=False):
        """
        @warning not done, just setting up how you'd do this
        @note per http://en.wikipedia.org/wiki/Euler_angles
         "This means for example that a convention named (YXZ) is the result
          of performing first an intrinsic Z rotation, followed by X and 
          Y rotations, in the moving axes (Note: the order of multiplication 
          of matrices is the opposite of the order in which they're
          applied to a vector)."
        """
        for rotation,theta in zip(self.rotations,self.thetas):
            ct = cos(radians(theta))
            st = sin(radians(theta))
            if   rotation==1:  p = dot(RotationX(ct,st),p)
            elif rotation==2:  p = dot(RotationX(ct,st),p)
            elif rotation==3:  p = dot(RotationX(ct,st),p)
        ###
        return p

    def RotationX(self,ct,st):
        matrix = array([[ 1.,  0., 0.],
                        [ ct,  0.,-st],
                        [-st,  0., ct]])

    def RotationY(self,ct,st):
        matrix = array([[ ct,  0., st],
                        [ 0.,  1., 0.],
                        [-st,  0., ct]])

    def RotationZ(self,ct,st):
        matrix = array([[ ct, st,  0.],
                        [-st, ct,  0.],
                        [  0., 0., 1.]])

    def __repr__(self):
        fields = ['CORD1R',self.g1,self.g2,self.g3]
        self.printCard(fields)

class CORD1R(Cord1x):
    type = 'CORD1R'
    """
    CORD1R CIDA G1A G2A G3A CIDB G1B G2B G3B
    """
    def __init__(self,nCoord=0,card=['CORD1R',0,0,0,0]):
        """
        Intilizes the CORD1R
        @param self   the object pointer
        @param nCoord the coordinate location on the line (there are possibly 2 coordinates on 1 card)
        @param card   a list version of the fields (1 CORD1R only)
        """
        if isinstance(card,list):
            assert len(card)==5
            card = BDF_Card(card)
        Cord1x.__init__(self,card)
        
        assert nCoord==0 or nChord==1
        nCoord *= 4  # 0 if the 1st coord, 4 if the 2nd

        ## reference coordinate system ID
        self.rid = card.field(2,0)
        if self.rid==0:
            self.isResolved = True
        
        ## the coordinate ID
        self.cid = card.field(1+nCoord)
        ## a Node at the origin
        self.g1  = card.field(2+nCoord)
        ## a Node on the z-axis
        self.g2  = card.field(3+nCoord)
        ## a Node on the xz-plane
        self.g3  = card.field(4+nCoord)
        assert self.g1 != self.g2
        assert self.g1 != self.g3
        assert self.g2 != self.g3

    def __repr__(self):
        fields = ['CORD1R',self.g1,self.g2,self.g3]
        self.printCard(fields)

    def crossReference(self,model):
        """
        Links self.rid to a coordinate system.
        @param self  the object pointer
        @param model the BDF object
        @warning
            Doesn't set rid to the coordinate system if it's in the global.
            This isn't a problem, it's meant to speed up the code in order
            to resolve extra coordinate systems.
        """
        self.isCrossReferenced = True
        self.g1 = model.Node(self.g1)
        self.g2 = model.Node(self.g2)
        self.g3 = model.Node(self.g3)

    def _transformGlobalToSelf(self):
        """
        this function takes a point in the global XYZ coordinate system with type=rectangular 
        and puts it in the rectangular coordinate system.
        """
        pass

    def resolveCid():
        self.setup()

class CORD2R(Cord2x):  # working for simple cases...
    type = 'CORD2R'
    def __init__(self,card=None,data=[0,0,  0.,0.,0.,  0.,0.,1., 1.,0.,0.]):
        #if isinstance(card,list):
        #    card = BDF_Card(card)
        Cord2x.__init__(self,card,data)
        #self.isResolved = False
        
        if card:
            ## coordinate system ID
            self.cid  = card.field(1)
            ## reference coordinate system ID
            self.rid = card.field(2,0)

            ## origin in a point relative to the rid coordinate system
            self.eo = array( card.fields(3,6 ,[0.,0.,0.]) )
            ## z-axis in a point relative to the rid coordinate system
            self.ez = array( card.fields(6,9 ,[0.,0.,0.]) )
            ## a point on the xz-plane relative to the rid coordinate system
            self.ex = array( card.fields(9,12,[0.,0.,0.]) )
        else:
            self.cid = data[0]
            self.rid = data[1]
            self.eo  = array(data[2:5])
            self.ez  = array(data[5:8])
            self.ex  = array(data[8:11])
        ###
        assert len(self.eo)==3
        assert len(self.ez)==3
        assert len(self.ex)==3
        if self.rid==0:
            self.isResolved = True
        self.setup()
        
    def crossReference(self,model):
        """
        Links self.rid to a coordinate system.
        @param self  the object pointer
        @param model the BDF object
        @warning
            Doesn't set rid to the coordinate system if it's in the global.
            This isn't a problem, it's meant to speed up the code in order
            to resolve extra coordinate systems.
        """
        self.isCrossReferenced = True
        if self.rid != 0:
            self.rid = model.Coord(self.rid)
        ###

    def _transformGlobalToSelf(self):
        """
        this function takes a point in the global XYZ coordinate system with type=rectangular 
        and puts it in the rectangular coordinate system.
        """
        pass

    def transformToGlobal(self,p,debug=False):
        """
        Transforms a point from the local coordinate system to the reference coordinate
        frames "global" coordinate system.  
        @note that this doesn't transform it from a coordinate system 
        \f[ \large [p_{global}]_{1x3} = [p_{local} -p_{origin}]_{1x3}[\Beta_{ij}]_{3x3}    \f]
        
        where   \f$ [\Beta]_{ij} \f$ is the tranformation matrix
        \f[ \large  [\Beta]_{ij} \left[ 
          \begin{array}{ccc}
              g_x \dot e_x & g_x \dot e_y &  g_x \dot e_z    \\
              g_y \dot e_x & g_y \dot e_y &  g_y \dot e_z    \\
              g_z \dot e_x & g_z \dot e_y &  g_z \dot e_z
          \end{array} \right]
        \f] 
        
        \f$ g \f$ is the global directional vector (e.g. \f$ g_x = [1,0,0]\f$ 
        \f$ e \f$ is the ith direction in the local coordinate system
        """
        if not self.isResolved:
            self.resolveCid()
        if self.cid==0:
            return p

        #p2 = p-self.eo
        
        # Bij = Bip*j
        ex = self.ex0
        ey = self.ey0
        ez = self.ez0
        if isinstance(self.rid,int):
            gx = array([1.,0.,0.])
            gy = array([0.,1.,0.])
            gz = array([0.,0.,1.])
        else:
            gx = self.rid.ex0
            gy = self.rid.ey0
            gz = self.rid.ez0
        ###
        
        print "gx = ",gx
        print "gy = ",gy
        print "gz = ",gz
        matrix = array([[dot(gx,ex),dot(gx,ey),dot(gx,ez)],
                        [dot(gy,ex),dot(gy,ey),dot(gy,ez)],
                        [dot(gz,ex),dot(gz,ey),dot(gz,ez)]])
        print "p = ",p
        print "matrix = \n",matrix
        p2 = dot(p+self.eo,matrix)
        p3 = p2#+self.eo
        print "eo = ",self.eo
        print "p2 = ",p2
        print '------------------------'
        print "p3 = ",p3
        
        #print str(self)
        if isinstance(self.rid,int):
            return p2
        else:
            return self.rid.transformToGlobal(p3)
        ###

    def __repr__(self):
        #eo = self.eo
        #ez = self.ez
        #ex = self.ex
        fields = ['CORD2R',self.cid,self.Rid()] +list(self.eo)+list(self.ez)+list(self.ex)

        #print "z1=%s z2=%s z3=%s" %(self.z1,self.z2,self.z3)
        #print "fields = ",fields
        return self.printCard(fields)


class CORD2C(Cord2x):  # not done...
    type = 'CORD2C'
    def __init__(self,card=None,data=[0,0,  0.,0.,0.,  0.,0.,1., 1.,0.,0.]):
        #if isinstance(card,list):
        #    card = BDF_Card(card)

        Cord2x.__init__(self,card,data)
        if card:
            self.cid  = card.field(1)
            self.rid = card.field(2,0)

            #print card
            self.eo = array( card.fields(3,6 ,[0.,0.,0.]) )
            self.ez = array( card.fields(6,9 ,[0.,0.,0.]) )
            self.ex = array( card.fields(9,12,[0.,0.,0.]) )
        else:
            self.cid = data[0]
            self.rid = data[0]
            self.eo  = array(data[2:5])
            self.ez  = array(data[5:8])
            self.ex  = array(data[8:11])
        ###
        
        assert len(self.eo)==3,self.eo
        assert len(self.ez)==3,self.ez
        assert len(self.ex)==3,self.ex
        
        #print "eo = ",self.eo
        #print "ez = ",self.ez
        #print "ex = ",self.ex

        self.ez0 = self.normalize(self.ez-self.eo)
        self.ex0 = self.normalize(self.ex-self.eo)
        self.ey0 = cross(self.ez0,self.ex0)

        #print card
        #print str(self)

        
    def _transformGlobalToSelf(self):
        """
        this function takes a point in the global XYZ coordinate system with type=rectangular 
        and puts it in a cylindrical coordinate system.
        """
        pass

    def xyz_To_RThetaZ(ex,ey,ez):
        """
        @code
        y       R
        |     / 
        |   /
        | / theta
        *------------x
        @endcode
        
        \f[ \large x = R \cos(\theta) \f]
        \f[ \large y = R \sin(\theta) \f]
        """
        pass

    def RThetaZ_To_xyz(er,et,ez):
        pass

    def point_RThetaZ_To_xyz(p,Rct,Rst):
        p2 = copy.deepcopy(p)
        p2[0] = p[0]*Rct+p[1]*Rst
        p2[1] = p[0]*Rst+p[1]*Rct
        return p2

    def transformToGlobal(self,p):
        (R,theta,z) = p
        raise Exception('cylindrical coordinate system...point R=%s theta=%s z=%s' %(R,theta,z))
        #assert R!=0.0
        Rct = R*cos(theta)
        Rst = R*sin(theta)
        originXYZ = point_RThetaZ_To_xyz(p,Rct,Rst)
        x2 = originXYZ[0] + (self.ex-self.eo[0])*Rct + (self.ey-self.eo[1])*Rst
        y2 = originXYZ[1] + (self.ex-self.eo[0])*Rst + (self.ey-self.eo[1])*Rct
        z2 = originXYZ[2]                                                       + (self.ez-self.eo[2])*z
        #return array([x2,y2,z2])
        
        p = array(x2,y2,z2)
        p2 = p-self.eo
        
        #Bij = Bip*j
        ex = self.ex0
        ey = self.ey0
        ez = self.ez0
        gx = array([1.,0.,0.])
        gy = array([0.,1.,0.])
        gz = array([0.,0.,1.])
        
        matrix = array([[dot(gx,ex),dot(gx,ey),dot(gx,ez)],
                        [dot(gy,ex),dot(gy,ey),dot(gy,ez)],
                        [dot(gz,ex),dot(gz,ey),dot(gz,ez)]])
        print "p = ",p
        print "matrix = ",matrix
        p2 = dot(p,matrix)
        p3 = p2+self.eo
        print "p2 = ",p2
        
        #print str(self)
        return p

    def __repr__(self):
        eo = self.eo
        ez = self.ez
        ex = self.ex
        fields = ['CORD2C',self.cid,self.Rid()] +list(eo)+list(ez) +list(ex)
        print "fields = ",fields
        raise Exception('not coded...')

        #print "z1=%s z2=%s z3=%s" %(self.z1,self.z2,self.z3)
        return self.printCard(fields)
