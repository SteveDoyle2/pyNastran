#import sys
from numpy import array,cross,dot
from numpy.linalg import norm

# my code
from pyNastran.bdf.errors import *
#from pyNastran.bdf.BDF_Card import BDF_Card
from baseCard import BaseCard

class Coord(BaseCard):
    def __init__(self,card,data):
        self.isCrossReferenced = False
        self.isResolved = False

    def setup(self):
        """
        \f[ e_{13}  = e_3 - e_1                \f]
        \f[ e_{12}  = e_2 - e_1                \f]
        \f[ k       = \frac{e_{12}}{|e_{12}|}  \f]
        \f[ j_{dir} = k \cross e_{13}          \f]
        \f[ j = \frac{j_{dir}}{|j_{dir}|}      \f]
        \f[ i = j \cross k                     \f]
        """
        try:
            assert len(self.e1)==3,self.e1
            assert len(self.e2)==3,self.e2
            assert len(self.e3)==3,self.e3
        except TypeError:
            msg = ''
            msg += "\ntype = %s\n" %(self.type)
            msg += "e1 = %s\n" %(self.e1)
            msg += "e2 = %s\n" %(self.e2)
            msg += "e3 = %s\n" %(self.e3)
            raise TypeError(msg)
        
        #print "e1 = ",self.e1
        #print "e2 = ",self.e2
        #print "e3 = ",self.e3
        
        #try:
        ## e_{13}
        e13 = self.e3-self.e1
        ## e_{12}
        e12 = self.e2-self.e1
        ## k = (G3 cross G1) normalized
        self.k = self.normalize(e12)
        #print "k = %s" %(self.k)
        #print "e13 = %s" %(e13)
        ## j = (k cross e13) normalized
        self.j = self.normalize(cross(self.k,e13))
        ## i = j cross k
        self.i = cross(self.j,self.k)
        #except TypeError:
        #    msg  = 'There is a problem handling these lines:\n'
        #    msg += '    self.k = self.normalize(self.e3-self.e1)\n'
        #    msg += '    self.ex0 = self.normalize(self.e2-self.e1)\n'
        #    msg += 'e1=%s Type=%s\n' %(self.e1,type(self.e1))
        #    msg += 'e2=%s Type=%s\n' %(self.e2,type(self.e2))
        #    msg += 'e3=%s Type=%s\n' %(self.e3,type(self.e3))
        #    #print msg
        #    raise CoordTypeError(msg)

    def normalize(self,v):
        #print "v = ",v
        normV = norm(v)
        assert normV>0.,'v=%s norm(v)=%s' %(v,normV)
        #print "normV = ",normV
        return v/normV

    def reprFields(self):
        return self.rawFields()

    #def resolveCid(self):
        #pass

class RectangularCoord(object):
    def _transformGlobalToSelf(self):
        """
        this function takes a point in the global XYZ coordinate system with
        type=rectangular and puts it in the rectangular coordinate system.
        """
        pass

class CylindricalCoord(object):
    """
    \f[ r        = \sqrt(x^2+y^2)      \f]
    \f[ \theta   = tan^-1(\frac{y}{x}) \f]
    \f[ z = z\f]

    \f[ x = r cos(\theta) \f]
    \f[ y = r sin(\theta) \f]
    \f[ z = z             \f]
    http://en.wikipedia.org/wiki/Cylindrical_coordinate_system
    @note \f$ \phi \f$ and \f$ \theta \f$ are flipped per wikipedia to be consistent with nastran's documentation
    @see refman.pdf
    """
    def _transformGlobalToSelf(self):
        """
        this function takes a point in the global XYZ coordinate system with
        type=rectangular and puts it in the cylindrical coordinate system.
        """
        raise Exception('not done...')

class SphericalCoord(object):
    """
    \f[ r = \rho = \sqrt(x^2+y^2+z^2)  \f]
    \f[ \theta   = tan^-1(\frac{y}{x}) \f]
    \f[ \phi     = cos^-1(\frac{z}{r}) \f]

    \f[ x = r cos(\theta)sin(\phi) \f]
    \f[ y = r sin(\theta)sin(\phi) \f]
    \f[ z = r cos(\phi)            \f]
    http://en.wikipedia.org/wiki/Spherical_coordinate_system
    @note \f$ \phi \f$ and \f$ \theta \f$ are flipped per wikipedia to be consistent with nastran's documentation
    @see refman.pdf
    """
    def _transformGlobalToSelf(self):
        """
        this function takes a point in the global XYZ coordinate system with
        type=rectangular and puts it in the spherical coordinate system.
        """
        raise Exception('not done...')
    
class Cord2x(Coord):
    def __init__(self,card,data):
        self.isResolved = False
        Coord.__init__(self,card,data)

        if card:
            ## coordinate system ID
            self.cid  = card.field(1)
            ## reference coordinate system ID
            self.rid = card.field(2,0)

            ## origin in a point relative to the rid coordinate system
            self.e1 = array( card.fields(3,6 ,[0.,0.,0.]) )
            ## z-axis in a point relative to the rid coordinate system
            self.e2 = array( card.fields(6,9 ,[0.,0.,0.]) )
            ## a point on the xz-plane relative to the rid coordinate system
            self.e3 = array( card.fields(9,12,[0.,0.,0.]) )
        else:
            self.cid = data[0]
            self.rid = data[1]
            self.e1  = array(data[2:5])
            self.e2  = array(data[5:8])
            self.e3  = array(data[8:11])
            assert len(data)==11,'data = %s' %(data)
        ###
        assert len(self.e1)==3
        assert len(self.e2)==3
        assert len(self.e3)==3
        if self.rid==0:
            self.isResolved = True
        self.setup()
    
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
        ## rid may be in a different coordinate system than ???
        self.e1 = self.transformToGlobal(self.e1)
        self.isResolved = True
        
        ## the axes are normalized, so we just assume they're points and resolve them
        self.e1 = self.rid.transformToGlobal(self.e1)
        self.i  = self.rid.transformToGlobal(self.i)
        self.j  = self.rid.transformToGlobal(self.j)
        self.k  = self.rid.transformToGlobal(self.k)
        self._transformGlobalToSelf()

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

    def transformToGlobal(self,p,debug=False):
        """
        Transforms a point from the local coordinate system to the reference coordinate
        frames "global" coordinate system.  
        @note that this doesn't transform it from a coordinate system 
        \f[ \large [p_{global}]_{1x3} = [p_{local} -p_{origin}]_{1x3}[\Beta_{ij}]_{3x3}    \f]
        
        where   \f$ [\Beta]_{ij} \f$ is the tranformation matrix
        \f[ \large  [\Beta]_{ij} \left[ 
          \begin{array}{ccc}
              g_x \dot i  &  g_x \dot j  &  g_x \dot k    \\
              g_y \dot i  &  g_y \dot j  &  g_y \dot k    \\
              g_z \dot i  &  g_z \dot j  &  g_z \dot k
          \end{array} \right]
        \f] 
        
        \f$ g   \f$ is the global directional vector (e.g. \f$ g_x = [1,0,0]\f$ 
        \f$ ijk \f$ is the ith direction in the local coordinate system
        """
        if not self.isResolved:
            self.resolveCid()
        if self.cid==0:
            return p

        #p2 = p-self.eo
        
        # Bij = Bip*j
        i = self.i
        j = self.j
        k = self.k
        if isinstance(self.rid,int):
            gx = array([1.,0.,0.])
            gy = array([0.,1.,0.])
            gz = array([0.,0.,1.])
        else:
            gx = self.rid.i
            gy = self.rid.j
            gz = self.rid.k
        ###
        
        print "gx = ",gx
        print "gy = ",gy
        print "gz = ",gz
        matrix = array([[dot(gx,i),dot(gx,j),dot(gx,k)],
                        [dot(gy,i),dot(gy,j),dot(gy,k)],
                        [dot(gz,i),dot(gz,j),dot(gz,k)]])
        print "p = ",p
        print "matrix = \n",matrix
        p2 = dot(p+self.e1,matrix)
        p3 = p2#+self.eo
        print "e1 = ",self.e1
        print "p2 = ",p2
        print '------------------------'
        print "p3 = ",p3
        
        #print str(self)
        if isinstance(self.rid,int):
            return p2
        else:
            return self.rid.transformToGlobal(p3)
        ###

    def Rid(self):
        """Returns the reference coordinate system self.rid"""
        if isinstance(self.rid,int):
            return self.rid
        return self.rid.cid


class Cord1x(Coord):
    def __init__(self,card,nCoord,data):
        Coord.__init__(self,card,data)

        self.isResolved = False
        if nCoord is not None:
            assert nCoord==0 or nCoord==1,'nCoord=|%s|' %(nCoord)
            nCoord *= 4  # 0 if the 1st coord, 4 if the 2nd

            ## the coordinate ID
            self.cid = card.field(1+nCoord)
            ## a Node at the origin
            self.g1  = card.field(2+nCoord)
            ## a Node on the z-axis
            self.g2  = card.field(3+nCoord)
            ## a Node on the xz-plane
            self.g3  = card.field(4+nCoord)
        else:
            self.cid = data[0]
            self.g1  = data[1]
            self.g2  = data[2]
            self.g3  = data[3]
            assert len(data)==4,'data = %s' %(data)
        ###
        assert self.g1 != self.g2
        assert self.g1 != self.g3
        assert self.g2 != self.g3

        self.e1=None; self.e2=None; self.e3=None
        self.i=None;  self.j=None;  self.k=None

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

        ## the origin
        self.e1 = self.g1.Position()
        ## a point on the z-axis
        self.e2 = self.g2.Position()
        ## a point on the xz-plane
        self.e3 = self.g3.Position()

    #def crossReference():
    #    ## the origin
    #    self.e1 = self.g1.Position()
    #    ## a point on the z-axis
    #    self.e2 = self.g2.Position()
    #    ## a point on the xz-plane
    #    self.e3 = self.g3.Position()
    #    self.setup()

    def resolveCid(self):
        self.setup()

    def G1(self):
        if isinstance(self.g1,int):
            return self.g1
        return self.g1.nid

    def G2(self):
        if isinstance(self.g2,int):
            return self.g2
        return self.g2.nid

    def G3(self):
        if isinstance(self.g3,int):
            return self.g3
        return self.g3.nid

    def GridIDs(self):
        """
        returns [g1,g2,g3]
        """
        grids = [self.G1(),self.G2(),self.G3()]
        return grids

class CORD3G(Coord):
    """
    Defines a general coordinate system using three rotational angles as functions of
    coordinate values in the reference coordinate system. The CORD3G entry is used with
    the MAT9 entry to orient material principal axes for 3-D composite analysis

    CORD3G CID METHOD FORM THETAID1 THETAID2 THETAID3 CIDREF
    CORD3G 100 E313   EQN  110      111      112      0
    """
    type = 'CORD3G'
    def __init__(self,card=[0,0,0,0,0,0,0],data=None):
        """
        Intilizes the CORD3G
        @param self   the object pointer
        @param card   a list version of the fields
        """
        if isinstance(card,list):
            assert len(card)==8
            card = BDF_Card(card)
        Coord.__init__(self,card,data)

        self.cid = card.field(1)
        method   = card.field(2)
        self.methodES  = method[0]
        self.methodInt = int(method[1:])
        assert self.ESmethod in ['E','S']
        assert 0 < self.methodInt < 1000

        self.form   = card.field(3,'EQN')
        self.thetas = card.field(4,7)
        assert len(thetas)==3,'thetas=%s' %(self.thetas)
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

    def rawFields(self):
        method = self.methodES+str(self.methodInt)
        fields = ['CORD3G',self.cid,method,self.form]+self.thetas+[self.cidRef]
        return fields

class CORD1R(Cord1x,RectangularCoord):
    type = 'CORD1R'
    """
    CORD1R CIDA G1A G2A G3A CIDB G1B G2B G3B
    @todo not done
    """

    def __init__(self,card=None,nCoord=0,data=None):
        """
        Intilizes the CORD1R
        @param self   the object pointer
        @param nCoord the coordinate location on the line (there are possibly 2 coordinates on 1 card)
        @param card   a list version of the fields (1 CORD1R only)
        """
        Cord1x.__init__(self,card,nCoord,data)

    def rawFields(self):
        fields = ['CORD1R',self.cid]+self.GridIDs()
        return fields

class CORD1C(Cord1x,CylindricalCoord):
    type = 'CORD1C'
    """
    CORD1C CIDA G1A G2A G3A CIDB G1B G2B G3B
    @todo not done
    """
    def __init__(self,card=None,nCoord=0,data=None):
        """
        Intilizes the CORD1R
        @param self   the object pointer
        @param nCoord the coordinate location on the line (there are possibly 2 coordinates on 1 card)
        @param card   a list version of the fields (1 CORD1R only)
        """
        Cord1x.__init__(self,card,nCoord,data)

    def rawFields(self):
        fields = ['CORD1C',self.cid]+self.GridIDs()
        return fields

class CORD1S(Cord1x,SphericalCoord):
    type = 'CORD1S'
    """
    CORD1S CIDA G1A G2A G3A CIDB G1B G2B G3B
    @todo not done
    """
    def __init__(self,card=None,nCoord=0,data=None):
        """
        Intilizes the CORD1S
        @param self   the object pointer
        @param nCoord the coordinate location on the line (there are possibly 2 coordinates on 1 card)
        @param data   a list version of the fields (1 CORD1S only)
        """
        Cord1x.__init__(self,card,nCoord,data)

    def rawFields(self):
        fields = ['CORD1S',self.cid]+self.GridIDs()
        return fields


class CORD2R(Cord2x):  # working for simple cases...
    type = 'CORD2R'
    def __init__(self,card=None,data=[0,0,  0.,0.,0.,  0.,0.,1., 1.,0.,0.]):
        Cord2x.__init__(self,card,data)
        
    def rawFields(self):
        rid = self.setBlankIfDefault(self.Rid(),0)
        fields = ['CORD2R',self.cid,rid] +list(self.e1)+list(self.e2)+list(self.e3)
        return fields


class CORD2S(Cord2x,SphericalCoord):  # not done...
    type = 'CORD2S'
    def __init__(self,card=None,data=[0,0,  0.,0.,0.,  0.,0.,1., 1.,0.,0.]):
        Cord2x.__init__(self,card,data)

    def rawFields(self):
        rid = self.setBlankIfDefault(self.Rid(),0)
        fields = ['CORD2R',self.cid,rid] +list(self.e1)+list(self.e2)+list(self.e3)
        return fields

class CORD2C(Cord2x,CylindricalCoord):  # not done...
    type = 'CORD2C'

    def __init__(self,card=None,data=[0,0,  0.,0.,0.,  0.,0.,1., 1.,0.,0.]):
        Cord2x.__init__(self,card,data)
        if card:
            self.cid = card.field(1)
            self.rid = card.field(2,0)

            #print card
            self.e1 = array( card.fields(3,6 ,[0.,0.,0.]) )
            self.e2 = array( card.fields(6,9 ,[0.,0.,0.]) )
            self.e3 = array( card.fields(9,12,[0.,0.,0.]) )
        else:
            self.cid = data[0]
            self.rid = data[1]
            self.e1  = array(data[2:5])
            self.e2  = array(data[5:8])
            self.e3  = array(data[8:11])
            assert len(data)==11,'data = %s' %(data)
        ###
        
        assert len(self.e1)==3,self.e1
        assert len(self.e2)==3,self.e2
        assert len(self.e3)==3,self.e3
        
        #print "e1 = ",self.e1
        #print "e2 = ",self.e2
        #print "e3 = ",self.e3
        self.setup()
        #print card
        #print str(self)

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
        i = self.i
        j = self.j
        k = self.k
        x2 = originXYZ[0] + (i-self.e1[0])*Rct + (j-self.e1[1])*Rst
        y2 = originXYZ[1] + (i-self.e1[0])*Rst + (j-self.e1[1])*Rct
        z2 = originXYZ[2]                                           + (k-self.e1[2])*z
        #return array([x2,y2,z2])
        
        p = array(x2,y2,z2)
        p2 = p-self.e1
        
        #Bij = Bip*j
        gx = array([1.,0.,0.])
        gy = array([0.,1.,0.])
        gz = array([0.,0.,1.])
        
        matrix = array([[dot(gx,i),dot(gx,j),dot(gx,k)],
                        [dot(gy,i),dot(gy,j),dot(gy,k)],
                        [dot(gz,i),dot(gz,j),dot(gz,k)]])
        print "p = ",p
        print "matrix = ",matrix
        p2 = dot(p,matrix)
        p3 = p2+self.e1
        print "p2 = ",p2
        
        #print str(self)
        return p

    def rawFields(self):
        rid = self.setBlankIfDefault(self.Rid(),0)
        fields = ['CORD2C',self.cid,rid] +list(self.e1)+list(self.e2)+list(self.e3)
        return fields
