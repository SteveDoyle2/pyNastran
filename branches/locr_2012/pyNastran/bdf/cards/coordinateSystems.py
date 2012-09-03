# pylint: disable=C0103,R0902,R0904,R0914
from __future__ import (nested_scopes, generators, division, absolute_import,
                        print_function, unicode_literals)
from math import sqrt, degrees, radians, atan2, acos, sin, cos
from itertools import izip

from numpy import array, cross, dot, transpose, zeros
from numpy.linalg import norm

from pyNastran.bdf.fieldWriter import set_blank_if_default
from pyNastran.bdf.cards.baseCard import BaseCard, BDFCard
from pyNastran.general.general import ListPrint


class Coord(BaseCard):
    type = 'COORD'

    def __init__(self, card, data):
        """
        Defines a general CORDxx object
        @param self the object pointer
        @param card a BDFCard object
        @param data a list analogous to the card
        """
        ## has the coordinate system been linked yet
        self.isCrossReferenced = False
        ## have all the transformation matricies been determined
        self.isResolved = False
        self.cid = None
        self.e1 = None
        self.e2 = None
        self.e3 = None

    def Cid(self):
        """returns the coordinate ID"""
        return self.cid

    def setup(self, debug=False):
        r"""
        \f[ e_{13}  = e_3 - e_1                \f]
        \f[ e_{12}  = e_2 - e_1                \f]
        \f[ k       = \frac{e_{12}}{|e_{12}|}  \f]
        \f[ j_{dir} = k \times e_{13}          \f]
        \f[ j = \frac{j_{dir}}{|j_{dir}|}      \f]
        \f[ i = j \times k                     \f]
        """
        try:
            assert len(self.e1) == 3, self.e1
            assert len(self.e2) == 3, self.e2
            assert len(self.e3) == 3, self.e3
            ## e_{13}
            e13 = self.e3 - self.e1
            ## e_{12}
            e12 = self.e2 - self.e1
            #print "e13 = %s" %(e13)
            #print "e12 = %s" %(e12)
        except TypeError:
            msg = ''
            msg += "\ntype = %s\n" % (self.type)
            msg += "\ncid  = %s\n" % (self.Cid())
            msg += "e1 = %s\n" % (self.e1)
            msg += "e2 = %s\n" % (self.e2)
            msg += "e3 = %s\n" % (self.e3)
            raise TypeError(msg)

        #print self
        #print "e1 = ",self.e1
        #print "e2 = ",self.e2
        #print "e3 = ",self.e3

        try:
            ## k = (G3 cross G1) normalized
            self.k = self.normalize(e12)
            ## j = (k cross e13) normalized
            self.j = self.normalize(cross(self.k, e13))
        except RuntimeError:
            print("---InvalidUnitVectorError---")
            print("Cp  = %s" % (self.Cid()))
            print("e1  = %s" % (self.e1))
            print("e2  = %s" % (self.e2))
            print("e3  = %s" % (self.e3))
            print("e13 = %s" % (e13))
            print("e12 = %s" % (e12))
            print("k = norm(e12)")
            print("k   = %s\n" % (self.k))
            print("j = norm(cross(k,e13))")
            raise
        try:
            ## i = j cross k
            self.i = cross(self.j, self.k)
        except RuntimeError:
            print("---InvalidUnitVectorError---")
            print("Cp  = %s" % (self.Cid()))
            print("e1  = %s" % (self.e1))
            print("e2  = %s" % (self.e2))
            print("e3  = %s" % (self.e3))
            print("e13 = %s" % (e13))
            print("e12 = %s" % (e12))
            print("k = norm(e12)")
            print("k   = %s\n" % (self.k))
            print("j = norm(cross(k,e13))")
            print("j   = %s" % (self.j))
            raise

        if debug:
            print("Cp = %s" % (self.Cid()))
            print("e1 = %s" % (self.e1))
            print("e2 = %s" % (self.e2))
            print("e3 = %s" % (self.e3))
            print('-----')
            print("e13 = %s" % (e13))
            print("e12 = %s" % (e12))
            print('-----')
            print("i   = %s" % (self.i))
            print("j   = %s" % (self.j))
            print("k   = %s\n" % (self.k))
            print('-----')

        #except TypeError:
        #    msg  = 'There is a problem handling these lines:\n'
        #    msg += '    self.k = self.normalize(self.e3-self.e1)\n'
        #    msg += '    self.ex0 = self.normalize(self.e2-self.e1)\n'
        #    msg += 'e1=%s Type=%s\n' %(self.e1,type(self.e1))
        #    msg += 'e2=%s Type=%s\n' %(self.e2,type(self.e2))
        #    msg += 'e3=%s Type=%s\n' %(self.e3,type(self.e3))
        #    #print msg
        #    raise CoordTypeError(msg)
        #print "k = %s" %(self.k)
        #print "e13 = %s" %(e13)

    def transformToLocal(self, p, matrix, debug=False):
        r"""
        Transforms the global point p to the local coordinate system
        @param self
          the object pointer
        @param p
          the point to transform
        @param matrix
          the transformation matrix to apply - created by transformToGlobal
        @param debug
          developer debug
        @note
          uses the matrix as there is no linking from a global coordinate
          system to the local
        @note
          the matrix that comes in is the local to global, so we need to invert
          the matrix. Luckily the inverse of a tranformation matrix
          \f$ [\phi] \f$ is the transpose of the matrix.
          \f[ p_{Global} = (p_{Local}-e_1 )[\phi]+e_1 \f]
          \f[ [phi]^{-1} = [phi]^T \f]
          (pc-e1) =(pG-e1)mT
          (pc-e1)*m = pG-e1
          (pc-e1)*m+e1 = pG

        @note
          be very careful of when you apply e1.  It gets removed whenever
          rotations are applied.  These equations need some TLC, but the
          methods are ok.
        """
        #pGlobal = self.transformToGlobal(p, debug=False)
        pCoord = dot(p - self.e1, transpose(matrix))
        pLocal = self.XYZtoCoord(pCoord)
        if debug:
            print("p      = %s" % (p))
            print("p-e1   = %s" % (p - self.e1))
            print("pLocal = %s\n" % (pLocal))
            print("pCoord = %s" % (pCoord))
        return pLocal
        #return pGlobal

    def normalize(self, v):
        """
        Normalizes v into a unit vector
        @param self
          the object pointer
        @param v
          the vector to normalize
        @retval
          nNorm v has been normalized
        """
        normV = norm(v)
        if not normV > 0.:
            raise RuntimeError('v=%s norm(v)=%s' % (v, normV))
        return v / normV

    def T(self):
        r"""
        Returns the 6x6 transformation
        \f[ \large  [\lambda] = [B_{ij}] \f]
        \f[
        [T] =
        \left[
          \begin{array}{cc}
          \lambda  & 0 \\
          0  & \lambda \\
          \end{array}
        \right]
        \f]
        """
        (a, matrix) = self.transformToGlobal(self.e1)
        t = zeros((6, 6))  # transformation matrix
        t[0:2, 0:2] = matrix
        t[3:5, 3:5] = matrix
        return t

    def reprFields(self):
        return self.rawFields()

    #def resolveCid(self):
        #pass


class RectangularCoord(object):
    def coordToXYZ(self, p):
        #print("p = %s" %(p))
        #print("e1 = %s" %(self.e1))
        return p + self.e1

    def XYZtoCoord(self, p):
        return p


class CylindricalCoord(object):
    r"""
    \f[ r        = \sqrt(x^2+y^2)      \f]
    \f[ \theta   = tan^-1(\frac{y}{x}) \f]
    \f[ z        = z                   \f]

    \f[ x = r cos(\theta) \f]
    \f[ y = r sin(\theta) \f]
    \f[ z = z             \f]
    \f[ p = [x,y,z] + e_1 \f]
    http://en.wikipedia.org/wiki/Cylindrical_coordinate_system
    @note
      \f$ \phi \f$ and \f$ \theta \f$ are flipped per wikipedia to be
      consistent with nastran's documentation
    @see refman.pdf
    """
    def coordToXYZ(self, p):
        r"""
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
        R = p[0]
        theta = radians(p[1])
        x = R * cos(theta)
        y = R * sin(theta)
        return array([x, y, p[2]]) + self.e1

    def XYZtoCoord(self, p):
        (x, y, z) = p
        theta = degrees(atan2(y, x))
        R = sqrt(x * x + y * y)
        return array([R, theta, z])


class SphericalCoord(object):
    r"""
    \f[ r = \rho = \sqrt(x^2+y^2+z^2)  \f]
    \f[ \theta   = tan^-1(\frac{y}{x}) \f]
    \f[ \phi     = cos^-1(\frac{z}{r}) \f]

    \f[ x = r cos(\theta)sin(\phi) \f]
    \f[ y = r sin(\theta)sin(\phi) \f]
    \f[ z = r cos(\phi)            \f]
    \f[ p = [x,y,z] + e_1          \f]
    http://en.wikipedia.org/wiki/Spherical_coordinate_system
    @note
      \f$ \phi \f$ and \f$ \theta \f$ are flipped per wikipedia to be
      consistent with nastran's documentation
    @see refman.pdf
    """
    def coordToXYZ(self, p):
        R = p[0]
        theta = radians(p[1])
        phi = radians(p[2])
        x = R * cos(theta) * sin(phi)
        y = R * sin(theta) * sin(phi)
        z = R * cos(phi)
        return array([x, y, z]) + self.e1

    def XYZtoCoord(self, p):
        (x, y, z) = p
        R = sqrt(x * x + y * y + z * z)
        theta = degrees(atan2(y, x))
        if R > 0:
            phi = degrees(acos(z / R))
        else:
            phi = 0.
        return array([R, theta, phi])


class Cord2x(Coord):
    def __init__(self, card, data):
        """
        defines the CORD2x class
        @param self the object pointer
        @param card a BDFCard object
        @param data a list analogous to the card
        """
        self.isResolved = False
        Coord.__init__(self, card, data)

        if card:
            ## coordinate system ID
            self.cid = card.field(1)
            ## reference coordinate system ID
            self.rid = card.field(2, 0)

            ## origin in a point relative to the rid coordinate system
            self.e1 = array(card.fields(3, 6, [0., 0., 0.]))
            ## z-axis in a point relative to the rid coordinate system
            self.e2 = array(card.fields(6, 9, [0., 0., 0.]))
            ## a point on the xz-plane relative to the rid coordinate system
            self.e3 = array(card.fields(9, 12, [0., 0., 0.]))
        else:
            self.cid = data[0]
            self.rid = data[1]
            self.e1 = array(data[2:5])
            self.e2 = array(data[5:8])
            self.e3 = array(data[8:11])
            assert len(data) == 11, 'data = %s' % (data)

        assert len(self.e1) == 3
        assert len(self.e2) == 3
        assert len(self.e3) == 3
        if self.rid == 0:
            self.isResolved = True
        self.setup()

    def resolveCid(self):
        """
        Turns the coordinate system from being a coordinate system of
        type 1 depending on a type 2 to a type 1 depending on nothing.

        More generally, takes a coordinate system that depends on multiple
        levels of coordinate systems and resolves them in order to resolve
        it's coordinate frame.  So, if a coordinate system is of type 2, this
        will effectively set rid to 0 with a type 2.

        This should handle any number of coordinate systems or coordinate
        system types assuming there is no circular references.
        """
        #print str(self)
        #print self.rid
        #print "cid=%s rid=%s"%(self.cid, self.Rid())
        if self.cid == 0 or isinstance(self.rid, int) or self.rid.isResolved:
            return  # rid=0 so already resolved
        elif self.rid.isResolved is False:  # rid
            msg  = ('there is a circular reference between Coord %s and '
                    'Coord %s' %(self.cid,self.Rid()))
            #assert self.rid.isCrossReferenced==False,msg)
            #print "  resolving cid=%s rid=%s" %(self.cid,self.Rid())
            self.rid.resolveCid()

        ## rid coordinate system is now resolved, time to resolve the cid
        ## coordinate system. rid may be in a different coordinate system
        ## than cid
        self.isResolved = True
        self.e1, matrix = self.transformToGlobal(self.e1)

        ## the axes are normalized, so assume they're points and
        ## resolve them in the XYZ system, but dont subtract e1 off
        ## (hence the False)
        self.e1, matrix = self.rid.transformToGlobal(self.e1)  # origin
        i, matrix = self.rid.transformToGlobal(self.i, False)
        j, matrix = self.rid.transformToGlobal(self.j, False)
        k, matrix = self.rid.transformToGlobal(self.k, False)

        ## the axes are global, so now we put them in the cid
        self.i = i
        self.j = j
        self.k = k

    def cross_reference(self, model):
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

    def transformToGlobal(self, p, resolveAltCoord=True, debug=False):
        r"""
        Transforms a point from the local coordinate system to the reference
        coordinate frames "global" coordinate system.

        \f[ \large [p_{global}]_{1\times 3} =
            [p_{local} -p_{origin}]_{1\times 3}[\beta_{ij}]_{3\times 3}  \f]

        where   \f$ [\beta]_{ij} \f$ is the transformation matrix
        \f[ \large  [\beta]_{ij} \left[
          \begin{array}{ccc}
              g_x \cdot i  &  g_x \cdot j  &  g_x \cdot k    \\
              g_y \cdot i  &  g_y \cdot j  &  g_y \cdot k    \\
              g_z \cdot i  &  g_z \cdot j  &  g_z \cdot k
          \end{array} \right]
        \f]

        * \f$ g  \f$ is the global directional vector (e.g. \f$ g_x = [1,0,0]\f$)
        * \f$ ijk \f$ is the ith direction in the local coordinate system
        
        @warning
          make sure you cross-reference before calling this
        @warning
          you probably shouldnt call this, call the Node methods Position
          and PositionWRT
        """
        if debug:
            print("p = %s" % (p))
            print("p-e1 = %s" % (p - self.e1))

        if not self.isResolved:
            self.resolveCid()
        if self.cid == 0:
            return p, array([[1., 0., 0.],
                             [0., 1., 0.],
                             [0., 0., 1.]])

        # the ijk axes arent resolved as R-theta-z, only points
        if resolveAltCoord:
            #print("p* = %s" %(p))
            p = self.coordToXYZ(p)
        #p2 = p-self.eo

        # Bij = Bip*j
        i = self.i
        j = self.j
        k = self.k
        if isinstance(self.rid, int):  # rid=0
            gx = array([1., 0., 0.])
            gy = array([0., 1., 0.])
            gz = array([0., 0., 1.])
        else:
            gx = self.rid.i
            gy = self.rid.j
            gz = self.rid.k
        ###

        matrix = array([[dot(gx, i), dot(gy, i), dot(gz, i)],
                        [dot(gx, j), dot(gy, j), dot(gz, j)],
                        [dot(gx, k), dot(gy, k), dot(gz, k)]])
        p2 = dot(p - self.e1, matrix)
        p3 = p2 + self.e1

        if debug:
            print("Cp = ", self.Cid())
            print("gx = %s" % (gx))
            print("gy = %s" % (gy))
            print("gz = %s" % (gz))
            print("p = %s" % (ListPrint(p)))
            print("matrix = \n", matrix)
            print("e1 = %s" % (ListPrint(self.e1)))
            print("p2 = %s" % (ListPrint(p2)))
            print('------------------------')
            print("p3 = %s\n" % (ListPrint(p3)))

        #print str(self)
        if isinstance(self.rid, int):
            return (p3, matrix)
        else:
              ## @todo do i need to multiply rid.transform(p3)[1]*matrix
            return (self.rid.transformToGlobal(p3)[0], matrix)
        ###

    def Rid(self):
        """Returns the reference coordinate system self.rid"""
        if isinstance(self.rid, int):
            return self.rid
        return self.rid.cid


class Cord1x(Coord):
    rid = 0  # used only for transform to global

    def __init__(self, card, nCoord, data):
        Coord.__init__(self, card, data)

        self.isResolved = False
        if nCoord is not None:
            assert nCoord == 0 or nCoord == 1, 'nCoord=|%s|' % (nCoord)
            nCoord *= 4  # 0 if the 1st coord, 4 if the 2nd

            ## the coordinate ID
            self.cid = card.field(1 + nCoord)
            ## a Node at the origin
            self.g1 = card.field(2 + nCoord)
            ## a Node on the z-axis
            self.g2 = card.field(3 + nCoord)
            ## a Node on the xz-plane
            self.g3 = card.field(4 + nCoord)
        else:
            self.cid = data[0]
            self.g1 = data[1]
            self.g2 = data[2]
            self.g3 = data[3]
            assert len(data) == 4, 'data = %s' % (data)
        ###
        assert self.g1 != self.g2
        assert self.g1 != self.g3
        assert self.g2 != self.g3

        self.e1 = None
        self.e2 = None
        self.e3 = None
        self.i = None
        self.j = None
        self.k = None

    def cross_reference(self, model):
        """
        Links self.rid to a coordinate system.
        @param self  the object pointer
        @param model the BDF object
        """
        self.isCrossReferenced = True
        ## grid point 1
        self.g1 = model.Node(self.g1)
        ## grid point 2
        self.g2 = model.Node(self.g2)
        ## grid point 3
        self.g3 = model.Node(self.g3)

    def resolveCid(self):
        """
        finds the position of the nodes used define the coordinate system
        and sets the ijk vectors
        """
        ## the origin
        self.e1 = self.g1.Position()
        ## a point on the z-axis
        self.e2 = self.g2.Position()
        ## a point on the xz-plane
        self.e3 = self.g3.Position()
        self.setup()

    def G1(self):
        if isinstance(self.g1, int):
            return self.g1
        return self.g1.nid

    def G2(self):
        if isinstance(self.g2, int):
            return self.g2
        return self.g2.nid

    def G3(self):
        if isinstance(self.g3, int):
            return self.g3
        return self.g3.nid

    def NodeIDs(self):
        """
        returns [g1,g2,g3]
        """
        grids = [self.G1(), self.G2(), self.G3()]
        return grids


class CORD3G(Coord):  # not done
    """
    Defines a general coordinate system using three rotational angles as
    functions of coordinate values in the reference coordinate system.
    The CORD3G entry is used with the MAT9 entry to orient material principal
    axes for 3-D composite analysis

    @code
    CORD3G CID METHOD FORM THETAID1 THETAID2 THETAID3 CIDREF
    CORD3G 100 E313   EQN  110      111      112      0
    @endcode
    """
    
    type = 'CORD3G'

    def __init__(self, card=[0, 0, 0, 0, 0, 0, 0], data=None):
        """
        Intilizes the CORD3G
        @param self   the object pointer
        @param card   a list version of the fields
        """
        if isinstance(card, list):
            assert len(card) == 8
            card = BDFCard(card)
        Coord.__init__(self, card, data)

        self.cid = card.field(1)
        method = card.field(2)
        self.methodES = method[0]
        self.methodInt = int(method[1:])
        assert self.methodES in ['E', 'S']
        assert 0 < self.methodInt < 1000

        self.form = card.field(3, 'EQN')
        self.thetas = card.field(4, 7)
        assert len(self.thetas) == 3, 'thetas=%s' % (self.thetas)
        self.cidRef = card.field(7)

        # EQN for DEQATN, TABLE for TABLE3D
        assert self.form in ['EQN', 'TABLE']

    def cross_reference(self, model):
        self.cidRef = model.Coord(self.cidRef)

    def CidRef(self):
        if isinstance(self.cidRef, int):
            return self.cidRef
        return self.cidRef.cid

    def transformToGlobal(self, p, debug=False):
        """
        @warning not done, just setting up how you'd do this
        @note per http://en.wikipedia.org/wiki/Euler_angles
         "This means for example that a convention named (YXZ) is the result
         of performing first an intrinsic Z rotation, followed by X and
         Y rotations, in the moving axes (Note: the order of multiplication
         of matrices is the opposite of the order in which they're
         applied to a vector)."
        """
        for (rotation, theta) in izip(self.rotations, self.thetas):
            ct = cos(radians(theta))
            st = sin(radians(theta))
            if   rotation == 1:
                p = dot(self.RotationX(ct, st), p)
            elif rotation == 2:
                p = dot(self.RotationY(ct, st), p)
            elif rotation == 3:
                p = dot(self.RotationZ(ct, st), p)
        return p

    def RotationX(self, ct, st):
        matrix = array([[1., 0., 0.],
                        [ct, 0., -st],
                        [-st, 0., ct]])
        return matrix

    def RotationY(self, ct, st):
        matrix = array([[ct, 0., st],
                        [0., 1., 0.],
                        [-st, 0., ct]])
        return matrix

    def RotationZ(self, ct, st):
        matrix = array([[ct, st, 0.],
                        [-st, ct, 0.],
                        [0., 0., 1.]])
        return matrix

    def rawFields(self):
        method = self.methodES + str(self.methodInt)
        fields = (['CORD3G', self.cid, method, self.form] + self.thetas +
                  [self.CidRef()])
        return fields


class CORD1R(Cord1x, RectangularCoord):
    """
    CORD1R CIDA G1A G2A G3A CIDB G1B G2B G3B
    """
    type = 'CORD1R'

    def __init__(self, card=None, nCoord=0, data=None):
        """
        Intilizes the CORD1R
        @param self
          the object pointer
        @param nCoord
          the coordinate location on the line (there are possibly 2 coordinates
          on 1 card)
        @param card
          a list version of the fields (1 CORD1R only)
        """
        Cord1x.__init__(self, card, nCoord, data)

    def rawFields(self):
        fields = ['CORD1R', self.cid] + self.NodeIDs()
        return fields


class CORD1C(Cord1x, CylindricalCoord):
    """
    CORD1C CIDA G1A G2A G3A CIDB G1B G2B G3B
    """
    type = 'CORD1C'

    def __init__(self, card=None, nCoord=0, data=None):
        """
        Intilizes the CORD1R
        @param self
          the object pointer
        @param card
          a BDFCard object
        @param nCoord
          the coordinate location on the line (there are possibly 2 coordinates
          on 1 card)
        @param data
          a list version of the fields (1 CORD1R only)

        """
        Cord1x.__init__(self, card, nCoord, data)

    def rawFields(self):
        fields = ['CORD1C', self.cid] + self.NodeIDs()
        return fields


class CORD1S(Cord1x, SphericalCoord):
    type = 'CORD1S'
    """
    CORD1S CIDA G1A G2A G3A CIDB G1B G2B G3B
    """
    def __init__(self, card=None, nCoord=0, data=None):
        """
        Intilizes the CORD1S
        @param self
          the object pointer
        @param card
          a BDFCard object
        @param nCoord
          the coordinate location on the line (there are possibly 2 coordinates
          on 1 card)
        @param data
          a list version of the fields (1 CORD1S only)
        """
        Cord1x.__init__(self, card, nCoord, data)

    def rawFields(self):
        fields = ['CORD1S', self.cid] + self.NodeIDs()
        return fields


class CORD2R(Cord2x, RectangularCoord):
    type = 'CORD2R'

    def __init__(self, card=None,
                 data=[0, 0, 0., 0., 0., 0., 0., 1., 1., 0., 0.]):
        """
        Intilizes the CORD2R
        @param self
          the object pointer
        @param card
          a BDFCard object
        @param data
          a list version of the fields (1 CORD2R only)
        """
        Cord2x.__init__(self, card, data)

    def rawFields(self):
        rid = set_blank_if_default(self.Rid(), 0)
        fields = ['CORD2R', self.cid, rid] + list(self.e1) + list(
            self.e2) + list(self.e3)
        return fields


class CORD2S(Cord2x, SphericalCoord):
    type = 'CORD2S'

    def __init__(self, card=None, data=None):
        """
        Intilizes the CORD2R
        @param self
          the object pointer
        @param card
          a BDFCard object
        @param data
          a list version of the fields (1 CORD2S only)
        """
        Cord2x.__init__(self, card, data)

    def rawFields(self):
        rid = set_blank_if_default(self.Rid(), 0)
        fields = (['CORD2S', self.cid, rid] + list(self.e1) + list(self.e2) +
                  list(self.e3))
        return fields


class CORD2C(Cord2x, CylindricalCoord):
    type = 'CORD2C'

    def __init__(self, card=None, data=None):
        """
        Intilizes the CORD2C
        @param self
          the object pointer
        @param card
          a BDFCard object
        @param data
          a list version of the fields (1 CORD2C only)
        """
        Cord2x.__init__(self, card, data)

    def rawFields(self):
        rid = set_blank_if_default(self.Rid(), 0)
        fields = (['CORD2C', self.cid, rid] + list(self.e1) + list(self.e2) +
                  list(self.e3))
        return fields
