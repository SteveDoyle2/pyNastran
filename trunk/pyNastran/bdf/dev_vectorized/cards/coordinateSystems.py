# pylint: disable=C0103,R0902,R0904,R0914
"""
All coordinate cards are defined in this file.  This includes:

 * CORD1R
 * CORD1C
 * CORD1S
 * CORD2R
 * CORD2C
 * CORD2S
"""
from __future__ import (nested_scopes, generators, division, absolute_import,
                        print_function, unicode_literals)
from math import sqrt, degrees, radians, atan2, acos, sin, cos
from itertools import izip

from numpy import array, cross, dot, transpose, zeros, vstack
from numpy.linalg import norm

from pyNastran.bdf.fieldWriter import set_blank_if_default
from pyNastran.bdf.cards.baseCard import BaseCard
from pyNastran.utils.dev import list_print
from pyNastran.bdf.bdfInterface.assign_type import (integer, integer_or_blank,
    double_or_blank, string_or_blank)
from pyNastran.bdf.fieldWriter import print_card_8


def normalize(v):
    r"""
    Normalizes v into a unit vector.

    :param self: the coordinate system object
    :param v:    the vector to normalize

    :returns:  normalized v

    .. math:: v_{norm} = \frac{v}{\lvert v \lvert}
    """
    normV = norm(v)
    if not normV > 0.:
        raise RuntimeError('v=%s norm(v)=%s' % (v, normV))
    return v / normV


class Coord(BaseCard):
    type = 'COORD'

    def __repr__(self):
        aaa
        return write_bdf(size=8, is_double=False)

    def __str__(self):
        return self.__repr__()

    def __init__(self, card, data):
        """
        Defines a general CORDxx object

        :param self: the coordinate system object
        :param card: a BDFCard object
        :param data: a list analogous to the card
        """
        #: have all the transformation matricies been determined
        self.isResolved = False
        self.cid = None
        self.e1 = None
        self.e2 = None
        self.e3 = None
        self.i = None
        self.j = None
        self.k = None
        self.origin = None

    def Cid(self):
        """Gets the coordinate ID"""
        return self.cid

    def setup(self, debug=False):
        """
        .. math:: e_{13} = e_3 - e_1

        .. math:: e_{12} = e_2 - e_1

        .. math:: k = \frac{e_{12}}{|e_{12}|}

        .. math:: j_{dir} = k \times e_{13}

        .. math:: j = \frac{j_{dir}}{|j_{dir}|}

        .. math:: i = j \times k
        """
        try:
            assert len(self.e1) == 3, self.e1
            assert len(self.e2) == 3, self.e2
            assert len(self.e3) == 3, self.e3
        except AssertionError:
            msg = 'Invalid Vector Length\n'
            msg += "type = %s\n" % (self.type)
            msg += "cid  = %s\n" % (self.Cid())
            msg += "e1 = %s\n" % str(self.e1)
            msg += "e2 = %s\n" % str(self.e2)
            msg += "e3 = %s\n" % str(self.e3)
            raise RuntimeError(msg)

        if self.Cid() == 0:
            self.origin = array([0., 0., 0.], dtype='float64')
            self.i = array([1., 0., 0.], dtype='float64')
            self.j = array([0., 1., 0.], dtype='float64')
            self.k = array([0., 0., 1.], dtype='float64')
            return

        if not self.isResolved:
            self.rid.setup()

        if self.Rid() == 0:
            self.origin = self.e1
            e1 = self.e1
            e2 = self.e2
            e3 = self.e3
        else:
            self.origin, rid_matrix = self.rid.transformToGlobal(self.e1)
            e1 = self.origin
            e2, rid_matrix = self.rid.transformToGlobal(self.e2)
            e3, rid_matrix = self.rid.transformToGlobal(self.e3)

        try:
            # e_{13}
            e13 = e3 - e1
            # e_{12}
            e12 = e2 - e1
            #print("e13 = %s" % e13)
            #print("e12 = %s" % e12)
        except TypeError:
            msg = ''
            msg += "\ntype = %s\n" % (self.type)
            msg += "\ncid  = %s\n" % (self.Cid())
            msg += "e1 = %s\n" % str(e1)
            msg += "e2 = %s\n" % str(e2)
            msg += "e3 = %s\n" % str(e3)
            raise TypeError(msg)

        try:
            #: k = (G3 cross G1) normalized
            self.k = normalize(e12)
        except RuntimeError:
            print("---InvalidUnitVectorError---")
            print("Cp  = %s" % (self.Cid()))
            print("e1  = %s" % (self.e1))
            print("e2  = %s" % (self.e2))
            print("e3  = %s" % (self.e3))
            print("e1* = %s" % (e1))
            print("e2* = %s" % (e2))
            print("e3* = %s" % (e3))
            print("e13 = %s" % (e13))
            print("e12 = %s" % (e12))
            print("k = normalize(e12)")
            raise

        try:
            # j = (k cross e13) normalized
            self.j = normalize(cross(self.k, e13))
        except RuntimeError:
            print("---InvalidUnitVectorError---")
            print("Cp  = %s" % (self.Cid()))
            print("e1  = %s" % (self.e1))
            print("e2  = %s" % (self.e2))
            print("e3  = %s" % (self.e3))
            print("e1* = %s" % (e1))
            print("e2* = %s" % (e2))
            print("e3* = %s" % (e3))
            print("e13 = %s" % (e13))
            print("e12 = %s" % (e12))
            print("k = norm(e12)")
            print("k   = %s\n" % (self.k))
            print("j* = cross(k,e13)")
            print("j*  = %s" % (cross(self.k, e13)))
            print("j = norm(cross(k,e13))\n")
            raise

        try:
            #: i = j cross k
            self.i = cross(self.j, self.k)
        except RuntimeError:
            print("---InvalidUnitVectorError---")
            print("Cp  = %s" % (self.Cid()))
            print("e1  = %s" % (self.e1))
            print("e2  = %s" % (self.e2))
            print("e3  = %s" % (self.e3))
            print("e13 = %s" % (e13))
            print("e12 = %s" % (e12))
            print("k = normalize(e12)")
            print("k   = %s\n" % (self.k))
            print("j = norm(cross(k,e13))")
            print("j   = %s" % (self.j))
            raise

        if debug:
            print("\nCid = %s" % (self.Cid()))
            print("Rid = %s" % (self.Rid()))
            print("e1 = %s" % (self.e1))
            print("e2 = %s" % (self.e2))
            print("e3 = %s" % (self.e3))
            print("e1* = [%g, %g, %g]" % tuple(e1))
            print("e2* = [%g, %g, %g]" % tuple(e2))
            print("e3* = [%g, %g, %g]" % tuple(e3))
            print('-----')
            print("e13 = %s" % (e13))
            print("e12 = %s" % (e12))
            print('-----')
            print("i   = %s len=%s"   % (str(self.i), norm(self.i)))
            print("j   = %s len=%s"   % (str(self.j), norm(self.j)))
            print("k   = %s len=%s\n" % (str(self.k), norm(self.k)))
            print('-----')

    def transformToGlobal(self, p, debug=False):
        r"""
        Transforms a point from the local coordinate system to the reference
        coordinate frames "global" coordinate system.

        .. math::
            [p_{global}]_{1\times 3} =
            [p_{local}]_{1\times 3}[\beta_{ij}]_{3\times 3} + [p_{origin}]

        where :math:`[\beta]_{ij}` is the transformation matrix

        .. math::
          [\beta]_{ij} = \left[
          \begin{array}{ccc}
              g_x \cdot i  &  g_x \cdot j  &  g_x \cdot k    \\
              g_y \cdot i  &  g_y \cdot j  &  g_y \cdot k    \\
              g_z \cdot i  &  g_z \cdot j  &  g_z \cdot k
          \end{array} \right]

        * :math:`g` is the global directional vector (e.g. :math:`g_x = [1,0,0]`)
        * :math:`ijk` is the math:`i^{th}` direction in the local coordinate system

        :param self:    the coordinate system object
        :param p:       the point to be transformed in the local frame.  Type=1x3 NUMPY.NDARRAY
        :param debug:   developer debug (default=False)
        :returns p2:  the point in the global frame.  Type=1x3 NUMPY.NDARRAY
        :returns beta:  the rotation matrix.  Type=6x6 NUMPY.NDARRAY

        .. warning:: make sure you cross-reference before calling this
        .. warning:: you probably shouldnt call this, call the Node methods
                     Position and PositionWRT
        """
        if debug:
            print("p = %s" % p)

        if self.cid == 0:
            return p, array([[1., 0., 0.],
                             [0., 1., 0.],
                             [0., 0., 1.]], dtype='float64')
        if not self.isResolved:
            self.rid.setup()

        # the ijk axes arent resolved as R-theta-z, only points
        p2 = self.coordToXYZ(p)

        if self.i is None:
            raise RuntimeError("Local unit vectors haven't been set.\nType=%r cid=%s rid=%s" % (self.type, self.cid, self.rid))
        # Bij = Bip*j
        #i = self.i
        #j = self.j
        #k = self.k
        #gx = array([1., 0., 0.], dtype='float64')
        #gy = array([0., 1., 0.], dtype='float64')
        #gz = array([0., 0., 1.], dtype='float64')

        #matrix = array([[dot(gx, i), dot(gy, i), dot(gz, i)],
        #                [dot(gx, j), dot(gy, j), dot(gz, j)],
        #                [dot(gx, k), dot(gy, k), dot(gz, k)]], dtype='float64')
        matrix = vstack([self.i, self.j, self.k])

        # rotate point p2 from the local frame to the global frame
        # shift point p by origin is in global coordinates
        #p3 = p2[0]*self.i + p2[1]*self.j + p2[2]*self.k + self.origin
        p3 = dot(p2, matrix) + self.origin

        if debug:
            print("Cp = ", self.Cid())
            print("p = %s" % (list_print(p)))
            print("matrix = \n", matrix)
            print("e1 = %s" % (list_print(self.e1)))
            print("p2 = %s" % (list_print(p2)))
            print('------------------------')
            print("p3 = %s" % (list_print(p3)))

        return (p3, matrix)

    def transformToLocal(self, p, beta, debug=False):
        r"""
        Transforms the global point p to the local coordinate system

        :param self:   the coordinate system object
        :param p:      the point to transform
        :param beta:   the transformation matrix to apply - created by
                       transformToGlobal
        :param debug:  developer debug

        .. note::  uses the matrix as there is no linking from a global
                   coordinate system to the local

        .. note::  the matrix that comes in is the local to global, so we need
                   to invert the matrix. The inverse of the tranformation
                   matrix :math:`[\beta]` is the transpose of the matrix.

        .. math:: p_{global} = (p_{coord})[\beta] + p_{origin}

        .. math:: [\beta]^{-1} = [\beta]^T

        .. math:: p_{coord} = (p_{global} -p_{origin}) [\beta]^T

        .. math:: p_{local} = transform(p_{coord})

        Where transform(x) depends on the rectangular, cylindrical, or
        spherical coordinate system
        """
        #betaT = hstack([self.i,self.j,self.k])  # verify
        #pGlobal = self.transformToGlobal(p, debug=False)
        if self.origin is None:
            raise RuntimeError('Origin=%s; Cid=%s Rid=%s' % (self.origin, self.cid, self.Rid()))
        pCoord = dot(p - self.origin, transpose(beta))
        pLocal = self.XYZtoCoord(pCoord)
        if debug:
            print("p        = %s" % p)
            print("p-origin = %s" % (p - self.origin))
            print("pCoord = %s" % pCoord)
            print("pLocal = %s\n" % pLocal)
        return pLocal

    def beta(self):
        r"""
        Gets the 3 x 3 transformation

        .. math:: [\lambda] = [B_{ij}]
        """
        if self.cid == 0:
            return p, array([[1., 0., 0.],
                             [0., 1., 0.],
                             [0., 0., 1.]], dtype='float64')
        matrix = vstack([self.i, self.j, self.k])
        return matrix

    def beta_n(self, n):
        r"""
        Gets the 3n x 3n transformation

        .. math:: [\lambda] = [B_{ij}]
        """
        assert n < 10, 'n=%r' % n
        matrix = self.beta()
        t = zeros((3*n, 3*n), dtype='float64')  # transformation matrix
        for i in xrange(n):
            t[i*3:i*3+2, i*3:i*3+2] = matrix[0:2, 0:2]
        return t

    def T(self):
        r"""
        Gets the 6 x 6 transformation

        .. math:: [\lambda] = [B_{ij}]

        .. math::
          [T] =
          \left[
            \begin{array}{cc}
            \lambda  & 0 \\
            0  & \lambda \\
            \end{array}
          \right]
        """
        return self.beta_n(2)

    def reprFields(self):
        return self.rawFields()


class RectangularCoord(object):
    def coordToXYZ(self, p):
        """
        :param self:  the coordinate system object
        :returns xyz: the point in the local coordinate system
        """
        return p

    def XYZtoCoord(self, p):
        """
        :param self:  the coordinate system object
        :returns xyz: the delta xyz point in the local coordinate system
        """
        return p


class CylindricalCoord(object):
    r"""
    .. math:: r      = \sqrt(x^2+y^2)
    .. math:: \theta = tan^{-1}\left(\frac{y}{x}\right)
    .. math:: z      = z

    .. math:: x = r cos(\theta)
    .. math:: y = r sin(\theta)
    .. math:: z = z
    .. math:: p = [x,y,z] + e_1

    http://en.wikipedia.org/wiki/Cylindrical_coordinate_system

    .. note:: :math`\phi` and :math:`\theta` are flipped per wikipedia
              to be consistent with nastran's documentation

    .. _msc:  http://simcompanion.mscsoftware.com/resources/sites/MSC/content/meta/DOCUMENTATION/9000/DOC9188/~secure/refman.pdf?token=WDkwz5Q6v7LTw9Vb5p+nwkbZMJAxZ4rU6BoR7AHZFxi2Tl1QdrbVvWj00qmcC4+S3fnbL4WUa5ovbpBwGDBt+zFPzsGyYC13zvGPg0j/5SrMF6bnWrQoTGyJb8ho1ROYsm2OqdSA9jVceaFHQVc+tJq4b49VogM4dZBxyi/QrHgdUgPFos8BAL9mgju5WGk8yYcFtRzQIxU=
    .. seealso:: `MSC Reference Manual (pdf) <`http://simcompanion.mscsoftware.com/resources/sites/MSC/content/meta/DOCUMENTATION/9000/DOC9188/~secure/refman.pdf?token=WDkwz5Q6v7LTw9Vb5p+nwkbZMJAxZ4rU6BoR7AHZFxi2Tl1QdrbVvWj00qmcC4+S3fnbL4WUa5ovbpBwGDBt+zFPzsGyYC13zvGPg0j/5SrMF6bnWrQoTGyJb8ho1ROYsm2OqdSA9jVceaFHQVc+tJq4b49VogM4dZBxyi/QrHgdUgPFos8BAL9mgju5WGk8yYcFtRzQIxU=>`_.
    """
    def coordToXYZ(self, p):
        r"""
        ::

          y       R
          |     /
          |   /
          | / theta
          *------------x

        .. math:: x = R \cos(\theta)
        .. math:: y = R \sin(\theta)

        :param self:  the coordinate system object
        :returns xyz: the point in the local coordinate system
        """
        R = p[0]
        theta = radians(p[1])
        x = R * cos(theta)
        y = R * sin(theta)
        return array([x, y, p[2]], dtype='float64')

    def XYZtoCoord(self, p):
        """
        :param self:  the coordinate system object
        :returns xyz: the delta xyz point in the local coordinate system
        """
        (x, y, z) = p
        theta = degrees(atan2(y, x))
        R = sqrt(x * x + y * y)
        return array([R, theta, z], dtype='float64')


class SphericalCoord(object):
    r"""
    .. math:: r = \rho = \sqrt(x^2+y^2+z^2)

    .. math:: \theta   = \cos^{-1}\left(\frac{z}{r}\right)

    .. math:: \phi     = \tan^{-1}\left(\frac{y}{x}\right)

    .. math:: x = r \sin(\theta)\cos(\phi)

    .. math:: y = r \sin(\theta)\sin(\phi)

    .. math:: z = r \cos(\theta)

    .. math:: p = [x,y,z]

    .. seealso:: http://en.wikipedia.org/wiki/Spherical_coordinate_system

    .. seealso:: `MSC Reference Manual (pdf) <`http://simcompanion.mscsoftware.com/resources/sites/MSC/content/meta/DOCUMENTATION/9000/DOC9188/~secure/refman.pdf?token=WDkwz5Q6v7LTw9Vb5p+nwkbZMJAxZ4rU6BoR7AHZFxi2Tl1QdrbVvWj00qmcC4+S3fnbL4WUa5ovbpBwGDBt+zFPzsGyYC13zvGPg0j/5SrMF6bnWrQoTGyJb8ho1ROYsm2OqdSA9jVceaFHQVc+tJq4b49VogM4dZBxyi/QrHgdUgPFos8BAL9mgju5WGk8yYcFtRzQIxU=>`_.
    """
    def XYZtoCoord(self, p):
        """
        :param self:  the coordinate system object
        :returns xyz: the loca XYZ point in the R, \theta, \phi coordinate system
        """
        (x, y, z) = p
        R = sqrt(x * x + y * y + z * z)
        phi = degrees(atan2(y, x))
        if R > 0:
            theta = degrees(acos(z / R))
        else:
            theta = 0.
        return array([R, theta, phi], dtype='float64')

    def coordToXYZ(self, p):
        """
        :param self:  the coordinate system object
        :returns xyz: the R, \theta, \phi point in the local XYZ coordinate system
        """
        R = p[0]
        theta = radians(p[1])
        phi = radians(p[2])
        x = R * sin(theta) * cos(phi)
        y = R * sin(theta) * sin(phi)
        z = R * cos(theta)
        return array([x, y, z], dtype='float64')


class Cord2x(Coord):
    def __init__(self, card, data):
        """
        Defines the CORD2x class

        :param self: the coordinate system object
        :param card: a BDFCard object
        :param data: a list analogous to the card
        """
        self.isResolved = False
        Coord.__init__(self, card, data)

        if card:
            #: coordinate system ID
            self.cid = integer(card, 1, 'cid')
            #: reference coordinate system ID
            self.rid = integer_or_blank(card, 2, 'rid', 0)

            #: origin in a point relative to the rid coordinate system
            self.e1 = array([double_or_blank(card, 3, 'e1x', 0.0),
                             double_or_blank(card, 4, 'e1y', 0.0),
                             double_or_blank(card, 5, 'e1z', 0.0)],
                            dtype='float64')
            #: z-axis in a point relative to the rid coordinate system
            self.e2 = array([double_or_blank(card, 6, 'e2x', 0.0),
                             double_or_blank(card, 7, 'e2y', 0.0),
                             double_or_blank(card, 8, 'e2z', 0.0)],
                            dtype='float64')
            #: a point on the xz-plane relative to the rid coordinate system
            self.e3 = array([double_or_blank(card, 9,  'e3x', 0.0),
                             double_or_blank(card, 10, 'e3y', 0.0),
                             double_or_blank(card, 11, 'e3z', 0.0)],
                            dtype='float64')
        else:
            self.cid = data[0]
            self.rid = data[1]
            self.e1 = array(data[2:5], dtype='float64')
            self.e2 = array(data[5:8], dtype='float64')
            self.e3 = array(data[8:11], dtype='float64')
            assert len(data) == 11, 'data = %s' % (data)

        assert len(self.e1) == 3
        assert len(self.e2) == 3
        assert len(self.e3) == 3

        #: the global axes
        self.i = None
        self.j = None
        self.k = None

        if self.rid == 0:
            self.isResolved = True
            self.setup()

    def _verify(self, xref):
        """
        Verifies all methods for this object work

        :param self: the CORD2x object pointer
        :param xref: has this model been cross referenced
        :type xref:  bool
        """
        cid = self.Cid()
        rid = self.Rid()
        assert isinstance(cid, int), 'cid=%r' % cid
        assert isinstance(rid, int), 'rid=%r' % rid

    def write_bdf(self, size=8, is_double=False):
        card = self.reprFields()
        if size == 8:
            return self.comment() + print_card_8(card)
        return self.comment() + print_card_16(card)

    def cross_reference(self, model):
        """
        Links self.rid to a coordinate system.

        :param self:  the coordinate system object
        :param model: the BDF object
        ..warning:: Doesn't set rid to the coordinate system if it's in the
                    global.  This isn't a problem, it's meant to speed up the
                    code in order to resolve extra coordinate systems.
        """
        if self.rid != 0:
            msg = ' which is required by %s cid=%s' % (self.type, self.cid)
            self.rid = model.Coord(self.rid, msg=msg)

    def Rid(self):
        """Gets the reference coordinate system self.rid"""
        if isinstance(self.rid, int):
            return self.rid
        return self.rid.cid


class Cord1x(Coord):
    rid = 0  # used only for transform to global

    def Rid(self):
        """Gets the reference coordinate system self.rid"""
        return self.rid

    def __init__(self, card, nCoord, data):
        Coord.__init__(self, card, data)
        self.isResolved = False
        if nCoord is not None:
            assert nCoord == 0 or nCoord == 1, 'nCoord=|%s|' % (nCoord)
            nCoord *= 4  # 0 if the 1st coord, 4 if the 2nd

            #: the coordinate ID
            self.cid = integer(card, 1 + nCoord, 'cid')
            #: a Node at the origin
            self.g1 = integer(card, 2 + nCoord, 'g1')
            #: a Node on the z-axis
            self.g2 = integer(card, 3 + nCoord, 'g2')
            #: a Node on the xz-plane
            self.g3 = integer(card, 4 + nCoord, 'g3')
        else:
            self.cid = data[0]
            self.g1 = data[1]
            self.g2 = data[2]
            self.g3 = data[3]
            assert len(data) == 4, 'data = %s' % (data)

        assert self.g1 != self.g2
        assert self.g1 != self.g3
        assert self.g2 != self.g3

        self.e1 = None
        self.e2 = None
        self.e3 = None
        self.i = None
        self.j = None
        self.k = None

    def to_CORD2x(self, model):
        """
        Converts a coordinate system from a CORD1x to a CORD2x

        :param self:  the coordinate system object
        :param model: a BDF model
        """
        rid1 = self.g1.cid
        rid2 = self.g2.cid
        rid3 = self.g2.cid
        type1 = self.g1.type
        #type2 = self.g2.type
        #type3 = self.g3.type

        p1 = self.g1.xyz
        if rid2 == rid1 and rid3 == rid1:
            # if all the same coordinate system, make them coord 2
            # the types are automatically the same then
            p2 = self.g2.xyz
            p3 = self.g3.xyz
        else:
            # otherwise make them the same rid as g1
            p2 = self.g2.PositionWRT(model, rid1)
            p3 = self.g3.PositionWRT(model, rid1)

        data = [type1, self.cid, rid1, list(p1) + list(p2) + list(p3)]

        if type1 == 'CORD2R':
            coord = CORD2R(card=None, data=data, comment=self.comment())
        elif type1 == 'CORD2C':
            coord = CORD2R(card=None, data=data, comment=self.comment())
        elif type1 == 'CORD2C':
            coord = CORD2R(card=None, data=data, comment=self.comment())
        else:
            raise RuntimeError('coordinate type of \n%s is %s' % (str(self), type1))
        model.coords[self.cid] = coord

    def _verify(self, xref):
        """
        Verifies all methods for this object work

        :param self: the CORD1x object pointer
        :param xref: has this model been cross referenced
        :type xref:  bool
        """
        cid = self.Cid()
        assert isinstance(cid, int), 'cid=%r' % cid

    def cross_reference(self, model):
        """
        Links self.rid to a coordinate system.

        :param self:  the coordinate system object
        :param model: the BDF object
        """
        msg = ' which is required by %s cid=%s' % (self.type, self.cid)
        #: grid point 1
        self.g1 = model.Node(self.g1, msg=msg)
        #: grid point 2
        self.g2 = model.Node(self.g2, msg=msg)
        #: grid point 3
        self.g3 = model.Node(self.g3, msg=msg)

    def setup(self, debug=False):
        """
        Finds the position of the nodes used define the coordinate system
        and sets the ijk vectors

        :param self: the coordinate system object
        """
        if self.isResolved:
            return

        #: the origin
        self.e1 = self.g1.Position()

        #: a point on the z-axis
        self.e2 = self.g2.Position()

        #: a point on the xz-plane
        self.e3 = self.g3.Position()

        # rid is resolved b/c e1, e2, & e3 are in global coordinates
        self.isResolved = True
        self.Coord.setup(debug=debug)

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
        Gets the integers for the node [g1,g2,g3]

        :param self: the coordinate system object
        """
        grids = [self.G1(), self.G2(), self.G3()]
        return grids

    def write_bdf(self, size=8, is_double=False):
        card = self.reprFields()
        return self.comment() + print_card_8(card)


class CORD3G(Coord):  # not done
    """
    Defines a general coordinate system using three rotational angles as
    functions of coordinate values in the reference coordinate system.
    The CORD3G entry is used with the MAT9 entry to orient material principal
    axes for 3-D composite analysis.::

      CORD3G CID METHOD FORM THETAID1 THETAID2 THETAID3 CIDREF
      CORD3G 100 E313   EQN  110      111      112      0
    """

    type = 'CORD3G'

    def __init__(self, card=None, data=None, comment=''):
        """
        Intilizes the CORD3G

        :param self: the CORD3G coordinate system object
        :param card: a list version of the fields
        """
        if comment:
            self._comment = comment

        Coord.__init__(self, card, data)

        self.cid = integer(card, 1, 'cid')
        method = string_or_blank(card, 2, 'E313')
        self.methodES = method[0]
        self.methodInt = int(method[1:])
        assert self.methodES in ['E', 'S'] # Euler / Space-Fixed
        assert 0 < self.methodInt < 1000

        self.form = string_or_blank(card, 3, 'form', 'EQN')
        self.thetas = [integer(card, 4, 'theta1'),
                       integer(card, 5, 'theta2'),
                       integer(card, 6, 'theta3')]
        assert len(self.thetas) == 3, 'thetas=%s' % (self.thetas)
        self.rid = integer_or_blank(card, 7, 'cidRef')
        assert len(card) <= 8, 'len(CORD3G card) = %i' % len(card)

        # EQN for DEQATN, TABLE for TABLE3D
        assert self.form in ['EQN', 'TABLE']

    def cross_reference(self, model):
        msg = ' which is required by %s cid=%s' % (self.type, self.cid)
        self.rid = model.Coord(self.rid, msg=msg)

    def Rid(self):
        if isinstance(self.rid, int):
            return self.rid
        return self.rid.cid

    def coord3g_transformToGlobal(self, p, debug=False):
        """
        :param self:  the coordinate system object
        :param p:     the point to transform.  TYPE=NUMPY.NDARRAY.
        :param debug: should debug messages be printed

        .. warning:: not done, just setting up how you'd do this
        .. note::    per http://en.wikipedia.org/wiki/Euler_angles
         "This means for example that a convention named (YXZ) is the result
         of performing first an intrinsic Z rotation, followed by X and
         Y rotations, in the moving axes (Note: the order of multiplication
         of matrices is the opposite of the order in which they're
         applied to a vector)."
        """
        if self.methodES == 'E':
            rotations = [int(i) for i in str(self.methodInt)]
            for (rotation, theta) in izip(rotations, self.thetas):
                ct = cos(radians(theta))
                st = sin(radians(theta))
                if   rotation == 1:
                    p = dot(self.RotationX(ct, st), p)
                elif rotation == 2:
                    p = dot(self.RotationY(ct, st), p)
                elif rotation == 3:
                    p = dot(self.RotationZ(ct, st), p)
                else:
                    raise RuntimeError('rotation=%s rotations=%s' % (rotation, rotations))
        elif self.methodES == 'S':
            raise RuntimeError('Space-Fixed rotation hasnt been implemented')
        else:
            raise RuntimeError('Invalid method; Use Euler or Space-Fixed.  MethodES=%r' % self.methodES)
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
        list_fields = (['CORD3G', self.cid, method, self.form] + self.thetas +
                  [self.Rid()])
        return list_fields


class CORD1R(Cord1x, RectangularCoord):
    """
    ::

    +-------+------+-----+-----+------+------+-----+------+-----+
    |   1   |   2  |  3  |  4  |   5  |  6   |  7  |  8   |  9  |
    +=======+======+=====+=====+======+======+=====+======+=====+
    |CORD1R | CIDA | G1A | G2A | CIDB | G1B  | G2B | G3B  |     |
    +-------+------+-----+-----+------+------+-----+------+-----+
    """
    type = 'CORD1R'
    Type = 'R'

    def __init__(self, card=None, nCoord=0, data=None, comment=''):
        """
        Intilizes the CORD1R

        :param self:   the CORD1R coordinate system object
        :param nCoord: the coordinate location on the line
                       (there are possibly 2 coordinates on 1 card)
        :param card:   a list version of the fields (1 CORD1R only)
        """
        Cord1x.__init__(self, card, nCoord, data)
        if comment:
            self._comment = comment

    def rawFields(self):
        list_fields = ['CORD1R', self.cid] + self.NodeIDs()
        return list_fields


class CORD1C(Cord1x, CylindricalCoord):
    """
    ::

    +-------+------+-----+-----+------+------+-----+------+-----+
    |   1   |   2  |  3  |  4  |   5  |  6   |  7  |  8   |  9  |
    +=======+======+=====+=====+======+======+=====+======+=====+
    |CORD1C | CIDA | G1A | G2A | CIDB | G1B  | G2B | G3B  |     |
    +-------+------+-----+-----+------+------+-----+------+-----+
    """
    type = 'CORD1C'
    Type = 'C'

    def __init__(self, card=None, nCoord=0, data=None, comment=''):
        """
        Intilizes the CORD1R

        :param self:   the CORD1C coordinate system object
        :param card:   a BDFCard object
        :param nCoord: the coordinate location on the line
                       (there are possibly 2 coordinates on 1 card)
        :param data:   a list version of the fields (1 CORD1R only)
        """
        Cord1x.__init__(self, card, nCoord, data)
        if comment:
            self._comment = comment

    def rawFields(self):
        list_fields = ['CORD1C', self.cid] + self.NodeIDs()
        return list_fields


class CORD1S(Cord1x, SphericalCoord):
    """
    ::

    +-------+------+-----+-----+------+------+-----+------+-----+
    |   1   |   2  |  3  |  4  |   5  |  6   |  7  |  8   |  9  |
    +=======+======+=====+=====+======+======+=====+======+=====+
    |CORD1S | CIDA | G1A | G2A | CIDB | G1B  | G2B | G3B  |     |
    +-------+------+-----+-----+------+------+-----+------+-----+
    """
    type = 'CORD1S'
    Type = 'S'

    def __init__(self, card=None, nCoord=0, data=None, comment=''):
        """
        Intilizes the CORD1S

        :param self:   the CORD1S coordinate system object
        :param card:   a BDFCard object
        :param nCoord: the coordinate location on the line
                       (there are possibly 2 coordinates on 1 card)
        :param data:   a list version of the fields (1 CORD1S only)
        """
        Cord1x.__init__(self, card, nCoord, data)
        if comment:
            self._comment = comment

    def rawFields(self):
        list_fields = ['CORD1S', self.cid] + self.NodeIDs()
        return list_fields


class CORD2R(Cord2x, RectangularCoord):
    type = 'CORD2R'
    Type = 'R'

    def __init__(self, card=None,
                 data=[0, 0, 0., 0., 0., 0., 0., 1., 1., 0., 0.], comment=''):
        """
        Intilizes the CORD2R

        +--------+-----+-----+-----+----+-----+----+----+-----+
        |    1   |   2 |  3  |  4  |  5 |  6  |  7 |  8 |  9  |
        +========+=====+=====+=====+====+=====+====+====+=====+
        | CORD2R | CID | RID | A1  | A2 | A3  | B1 | B2 |     |
        +--------+-----+-----+-----+----+-----+----+----+-----+
        |        | B3  | C1  | C2  | C3 |
        +--------+-----+-----+-----+----+

        :param self: the CORD2R coordinate system object
        :param card: a BDFCard object
        :param data: a list version of the fields (1 CORD2R only)
        """
        Cord2x.__init__(self, card, data)
        if comment:
            self._comment = comment

    def _verify(self, xref):
        """
        Verifies all methods for this object work

        :param self: the CORD2R object pointer
        :param xref: has this model been cross referenced
        :type xref:  bool
        """
        cid = self.Cid()
        rid = self.Rid()
        assert isinstance(cid, int), 'cid=%r' % cid
        assert isinstance(rid, int), 'rid=%r' % rid

    def rawFields(self):
        rid = set_blank_if_default(self.Rid(), 0)
        list_fields = ['CORD2R', self.cid, rid] + list(self.e1) + list(
            self.e2) + list(self.e3)
        return list_fields


class CORD2C(Cord2x, CylindricalCoord):
    type = 'CORD2C'
    Type = 'C'

    def __init__(self, card=None, data=None, comment=''):
        """
        Intilizes the CORD2C

        +--------+-----+-----+-----+----+-----+----+----+-----+
        |    1   |   2 |  3  |  4  |  5 |  6  |  7 |  8 |  9  |
        +========+=====+=====+=====+====+=====+====+====+=====+
        | CORD2C | CID | RID | A1  | A2 | A3  | B1 | B2 |     |
        +--------+-----+-----+-----+----+-----+----+----+-----+
        |        | B3  | C1  | C2  | C3 |
        +--------+-----+-----+-----+----+

        :param self: the CORD2C coordinate system object
        :param card: a BDFCard object
        :param data: a list version of the fields (1 CORD2C only)
        """
        Cord2x.__init__(self, card, data)
        if comment:
            self._comment = comment

    def rawFields(self):
        rid = set_blank_if_default(self.Rid(), 0)
        list_fields = (['CORD2C', self.cid, rid] + list(self.e1) +
                       list(self.e2) + list(self.e3))
        return list_fields


class CORD2S(Cord2x, SphericalCoord):
    type = 'CORD2S'
    Type = 'S'

    def __init__(self, card=None, data=None, comment=''):
        """
        Intilizes the CORD2R

        +--------+-----+-----+-----+----+-----+----+----+-----+
        |    1   |   2 |  3  |  4  |  5 |  6  |  7 |  8 |  9  |
        +========+=====+=====+=====+====+=====+====+====+=====+
        | CORD2S | CID | RID | A1  | A2 | A3  | B1 | B2 |     |
        +--------+-----+-----+-----+----+-----+----+----+-----+
        |        | B3  | C1  | C2  | C3 |
        +--------+-----+-----+-----+----+

        :param self: the CORD2S coordinate system object
        :param card: a BDFCard object
        :param data: a list version of the fields (1 CORD2S only)
        """
        Cord2x.__init__(self, card, data)
        if comment:
            self._comment = comment

    def rawFields(self):
        rid = set_blank_if_default(self.Rid(), 0)
        list_fields = (['CORD2S', self.cid, rid] + list(self.e1) +
                       list(self.e2) + list(self.e3))
        return list_fields
