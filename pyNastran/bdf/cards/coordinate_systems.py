# pylint: disable=R0902,R0904,R0914
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
from six.moves import zip, range
from math import sqrt, degrees, radians, atan2, acos, sin, cos

import numpy as np
from numpy.linalg import norm

from pyNastran.utils import integer_types
from pyNastran.bdf.field_writer_8 import set_blank_if_default
from pyNastran.bdf.cards.base_card import BaseCard
from pyNastran.bdf.bdf_interface.assign_type import (integer, integer_or_blank,
    double_or_blank, string_or_blank, string)
from pyNastran.bdf.field_writer_8 import print_card_8
from pyNastran.bdf.field_writer_16 import print_card_16
from pyNastran.bdf.field_writer_double import print_card_double


def normalize(v):
    r"""
    Normalizes v into a unit vector.

    :param v:    the vector to normalize

    :returns:  normalized v

    .. math:: v_{norm} = \frac{v}{\lvert v \lvert}
    """
    norm_v = norm(v)
    if not norm_v > 0.:
        raise RuntimeError('v=%s norm(v)=%s' % (v, norm_v))
    return v / norm_v


class Coord(BaseCard):
    type = 'COORD'

    def __init__(self):
        """
        Defines a general CORDxx object

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


    def XYZtoCoord(self, p):
        self.deprecated('XYZtoCoord', 'xyz_to_coord', '0.8')
        return self.xyz_to_coord(p)

    def coordToXYZ(self, p):
        self.deprecated('coordToXYZ', 'coord_to_xyz', '0.8')
        return self.coord_to_xyz(p)

    def Cid(self):
        """Gets the coordinate ID"""
        return self.cid

    def setup_global_cord2x(self):
        """
        Sets up a global CORD2R, CORD2S, CORD2C
        """
        #if self.Cid() == 0:
            #self.origin = np.array([0., 0., 0.], dtype='float64')
            #return
        #self.i = np.array([1., 0., 0.], dtype='float64')
        #self.j = np.array([0., 1., 0.], dtype='float64')
        #self.k = np.array([0., 0., 1.], dtype='float64')

        self.origin = self.e1
        e1 = self.e1
        e2 = self.e2
        e3 = self.e3

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
            print("k   = normalize(e12)")
            raise

        try:
            # j = (k cross e13) normalized
            self.j = normalize(np.cross(self.k, e13))
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
            print("k   = norm(e12)")
            print("k   = %s\n" % (self.k))
            print("j*  = cross(k, e13)")
            print("j*  = %s" % (np.cross(self.k, e13)))
            print("j   = norm(cross(k, e13))\n")
            raise

        try:
            #: i = j cross k
            self.i = np.cross(self.j, self.k)
        except RuntimeError:
            print("---InvalidUnitVectorError---")
            print("Cp  = %s" % (self.Cid()))
            print("Rid = %s" % (self.Rid()))
            print("e1  = %s" % (self.e1))
            print("e2  = %s" % (self.e2))
            print("e3  = %s" % (self.e3))
            print("e13 = %s" % (e13))
            print("e12 = %s" % (e12))
            print("k   = normalize(e12)")
            print("k   = %s\n" % (self.k))
            print("j   = norm(cross(k,e13))")
            print("j   = %s" % (self.j))
            raise

        #print('done setting up cid=%s rid=%s' % (self.cid, self.Rid()))

    def setup(self):
        r"""
        .. math::
          e_{13} = e_3 - e_1

        .. math::
          e_{12} = e_2 - e_1

        .. math::
          k = \frac{e_{12}}{\lvert e_{12} \rvert}

        .. math::
          j_{dir} = k \times e_{13}

        .. math::
          j = \frac{j_{dir}}{\lvert j_{dir} \rvert}

        .. math::
          i = j \times k
        """
        #if self.isResolved:
            #return
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
            self.origin = np.array([0., 0., 0.], dtype='float64')
            self.i = np.array([1., 0., 0.], dtype='float64')
            self.j = np.array([0., 1., 0.], dtype='float64')
            self.k = np.array([0., 0., 1.], dtype='float64')
            return

        #print('setting up cid=%s rid=%s' % (self.cid, self.Rid()))
        if not self.isResolved and self.type in ['CORD2R', 'CORD2C', 'CORD2S']:
            self.rid.setup()

        if self.Rid() == 0:
            self.origin = self.e1
            e1 = self.e1
            e2 = self.e2
            e3 = self.e3
        else:
            self.origin = self.rid_ref.transform_node_to_global(self.e1)
            e1 = self.origin
            e2 = self.rid_ref.transform_node_to_global(self.e2)
            e3 = self.rid_ref.transform_node_to_global(self.e3)

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
            print("k   = normalize(e12)")
            raise

        try:
            # j = (k cross e13) normalized
            self.j = normalize(np.cross(self.k, e13))
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
            print("k   = norm(e12)")
            print("k   = %s\n" % (self.k))
            print("j*  = cross(k, e13)")
            print("j*  = %s" % (np.cross(self.k, e13)))
            print("j   = norm(cross(k, e13))\n")
            raise

        try:
            #: i = j cross k
            self.i = np.cross(self.j, self.k)
        except RuntimeError:
            print("---InvalidUnitVectorError---")
            print("Cp  = %s" % (self.Cid()))
            print("Rid = %s" % (self.Rid()))
            print("e1  = %s" % (self.e1))
            print("e2  = %s" % (self.e2))
            print("e3  = %s" % (self.e3))
            print("e13 = %s" % (e13))
            print("e12 = %s" % (e12))
            print("k   = normalize(e12)")
            print("k   = %s\n" % (self.k))
            print("j   = norm(cross(k,e13))")
            print("j   = %s" % (self.j))
            raise

        if 0:
            print("\nCid = %s" % (self.Cid()))
            print("Rid = %s" % (self.Rid()))
            print("e1  = %s" % (self.e1))
            print("e2  = %s" % (self.e2))
            print("e3  = %s" % (self.e3))
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
        #print('done setting up cid=%s rid=%s' % (self.cid, self.Rid()))

    #def transform_force_to_global(self, F, M):
        #raise NotImplementedError('transform_force_to_global')
        #Fg = self.transform_vector_to_global(self, F)
        #Mg = self.transform_vector_to_global(self, M)

        #r = self.origin #  maybe a minus sign?
        #Mdelta = cross(r, Fg)
        #return Fg, Mg + Mdelta

    def transform_vector_to_global(self, p):
        """
        Transforms a generalized vector from the local frame to the
        global frame.  A generalized vector is unchanged when you shift
        it's point of application.  So:
          - Generalized Vectors (Force, Moment about the origin)
          - Not Generalized Vectors (node xyz, displacement, Moment)

        Parameters
        ----------
        p : (1,3) ndarray
            the vector in the local frame

        Returns
        -------
        p3 : (1,3) ndarray
            the vector in the global frame

        .. note:: Shifting the load application point of a force creates
                  a moment, but the force will be the same.
        """
        if self.cid == 0:
            return p

        if not self.isResolved:
            if isinstance(self.rid, int) and self.rid != 0:
                raise RuntimeError("BDF has not been cross referenced.")
            if self.type in ['CORD2R', 'CORD2C', 'CORD2S']:
                self.rid_ref.setup()
            else:
                self.setup()

        # the ijk axes arent resolved as R-theta-z, only points
        p2 = self.coord_to_xyz(p)

        if self.i is None:
            msg = "Local unit vectors haven't been set.\nType=%r cid=%s rid=%s" % (
                self.type, self.cid, self.rid)
            raise RuntimeError(msg)
        matrix = np.vstack([self.i, self.j, self.k])

        # rotate point p2 from the local frame to the global frame
        p3 = np.dot(p2, matrix)
        return p3

    def resolve(self):
        if not self.isResolved:
            if isinstance(self.rid, int) and self.rid != 0:
                raise RuntimeError("BDF has not been cross referenced.")
            if self.type in ['CORD2R', 'CORD2C', 'CORD2S']:
                self.rid_ref.setup()
            else:
                self.setup()

    def transform_vector_to_global_array(self, p):
        """
        Transforms a generalized vector from the local frame to the
        global frame.
        """
        if self.cid == 0:
            return p
        self.resolve()

        # the ijk axes arent resolved as R-theta-z, only points
        p2 = self.coord_to_xyz_array(p)

        if self.i is None:
            msg = "Local unit vectors haven't been set.\nType=%r cid=%s rid=%s" % (
                self.type, self.cid, self.rid)
            raise RuntimeError(msg)
        matrix = np.vstack([self.i, self.j, self.k])

        # rotate point p2 from the local frame to the global frame
        p3 = np.dot(p2, matrix)
        return p3

    def transform_node_to_global(self, xyz):
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

        Parameters
        ----------
        xyz : (1,3) ndarray
            the point in the local frame to be transformed

        Returns
        -------
        xyz_global : (1,3) ndarray
            the point in the global frame

        .. warning:: make sure you cross-reference before calling this
        .. warning:: you probably shouldnt call this, call the Node methods
                     Position and PositionWRT
        """
        if self.cid == 0:
            return xyz
        return self.transform_vector_to_global(xyz) + self.origin

    def _transform_node_to_global_array(self, xyz):
        if self.cid == 0:
            return xyz
        return self.transform_vector_to_global_array(xyz) + self.origin

    def _transform_node_to_local(self, xyz, beta):
        """
        Parameters
        ----------
        xyz : (1,3) ndarray
            the point in the global frame

        Returns
        -------
        xyz_local : (1,3) ndarray
            the point in the local frame
        """
        if self.origin is None:
            raise RuntimeError('Origin=%s; Cid=%s Rid=%s' % (self.origin, self.cid, self.Rid()))
        xyz_coord = np.dot(xyz - self.origin, np.transpose(beta))
        xyz_local = self.xyz_to_coord(xyz_coord)
        return xyz_local

    def _transform_node_to_local_array(self, xyz, beta):
        """
        Parameters
        ----------
        xyz : (n,3) ndarray
            the points in the global frame

        Returns
        -------
        xyz_local : (1,3) ndarray
            the point in the local frame
        """
        if self.origin is None:
            raise RuntimeError('Origin=%s; Cid=%s Rid=%s' % (self.origin, self.cid, self.Rid()))
        xyz_coord = np.dot(xyz - self.origin, np.transpose(beta))
        xyz_local = self.xyz_to_coord(xyz_coord)
        return xyz_local

    def transform_node_to_local(self, xyz):
        r"""
        Transforms the global point p to the local coordinate system

        Parameterse
        -----------
        xyz : (1,3) ndarray
            the point to transform

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
        beta = self.beta()
        return self._transform_node_to_local(xyz, beta)

    def transform_node_to_local_array(self, xyz):
        """
        Transforms the global point p to the local coordinate system
        """
        beta = self.beta()
        return self._transform_node_to_local_array(xyz, beta)

    def transform_vector_to_local(self, xyz):
        """
        see transform_node_to_local, but set the origin to <0, 0, 0>
        """
        beta = self.beta()
        xyz_coord = np.dot(xyz, np.transpose(beta))
        xyz_local = self.xyz_to_coord(xyz_coord)
        return xyz_local

    @property
    def global_to_local(self):
        r"""
        Gets the 3 x 3 global to local transform
        """
        return self.beta().T

    @property
    def local_to_global(self):
        """
        Gets the 3 x 3 local to global transform
        """
        return self.beta()

    def beta(self):
        r"""
        Gets the 3 x 3 transformation

        .. math:: [\lambda] = [B_{ij}]
        """
        if self.cid == 0:
            return np.array([[1., 0., 0.],
                             [0., 1., 0.],
                             [0., 0., 1.]], dtype='float64')
        matrix = np.vstack([self.i, self.j, self.k])
        return matrix

    def beta_n(self, n):
        r"""
        Gets the 3n x 3n transformation

        .. math:: [\lambda] = [B_{ij}]
        """
        assert n < 10, 'n=%r' % n
        matrix = self.beta()
        t = np.zeros((3*n, 3*n), dtype='float64')  # transformation matrix
        for i in range(n):
            t[i*3:i*3+2, i*3:i*3+2] = matrix[0:2, 0:2]
        return t

    def repr_fields(self):
        return self.raw_fields()

    def move_origin(self, xyz):
        """
        Move the coordinate system to a new origin while maintaining
        the orientation

        :param xyz: the new origin point to move the coordinate to in
                    the global coordinate system
        """
        if self.i == None:
            self.setup()

        xyz = _fix_xyz_shape(xyz)
        if self.type in ['CORD2R', 'CORD2C', 'CORD2S']:
            e1 = self.rid_ref.transform_node_to_global(self.e1)
            e2 = self.rid_ref.transform_node_to_global(self.e2)
            e3 = self.rid_ref.transform_node_to_global(self.e3)
            e12 = e2 - e1
            e13 = e3 - e1
            self.e1 = self.rid_ref.transform_node_to_local(xyz)
            self.e2 = self.rid_ref.transform_node_to_local(xyz + e12)
            self.e3 = self.rid_ref.transform_node_to_local(xyz + e13)
        else:
            raise RuntimeError('Cannot move %s; cid=%s' % (self.type, self.cid))
        self.origin = xyz


def _fix_xyz_shape(xyz, name='xyz'):
    """
    Checks the shape of a grid point location and fixes it if possible
    """
    xyz = np.asarray(xyz)
    if not isinstance(xyz, np.ndarray):
        msg = '%s must be type ndarray; type=%s' % (name, type(xyz))
        raise TypeError(msg)
    if xyz.shape == (1, 3):
        xyz.reshape(3, 1)

    if xyz.shape != (3,):
        msg = '%s must be 3 by 1 dim\n%s is %s' % (name, name, xyz.shape)
        raise RuntimeError(msg)
    return xyz


def define_spherical_cutting_plane(model, origin, rid, cids, thetas, phis):
    r"""
    Creates a series of coordinate systems defined as constant origin,
    with a series of theta and phi angles, which are defined about the
    aerodynamic axis <1, 0, 0>.  This is intended to be with a
    supersonic mach plane for calculating wave drag where:
        .. math::  \theta = \mu = \frac{1}{\sqrt(Mach^2 - 1)}
        .. math::  \phi = [-\pi, \pi]

    :param model: a BDF object
    :param origin: a (3,) ndarray defining the location of the origin in
                   the global coordinate frame
    :param rid:  the new spherical reference coordinate system id
    :param cids:  list of new coordinate system ids
    :param thetas:  list of thetas (in radians)
    :param phis:  list of phis (in radians)

    .. note:: creates 1 CORD2S and ncid CORD2R coordinate systems

    .. todo:: hasn't been tested...
    """
    if len(cids) != len(thetas):
        msg = 'len(cids)=%s len(thetas)=%s; must be equal' % (len(cids, len(thetas)))
        raise RuntimeError(msg)
    if len(cids) != len(phis):
        msg = 'len(cids)=%s len(phis)=%s; must be equal' % (len(cids, len(phis)))
        raise RuntimeError(msg)

    # check for dupliate coords
    assert rid not in model.coords
    for cid in cids:
        assert cid not in model.coords

    # create the spherical coordinate system
    origin = _fix_xyz_shape(origin, 'origin')
    e2 = origin + np.array([0., 0., 1.])
    e3 = origin + np.array([1., 0., 0.])
    card = ['CORD2S', rid, 0] + list(origin) + list(e2) + list(e3)
    model.add_card(card, card[0], is_list=True)

    # create the mach planes
    for cid, theta, phi in zip(cids, thetas, phis):
        e2 = [1.0, theta, 0.]
        e3 = [1.0, theta, phi]
        card = ['CORD2R', cid, rid] + [0., 0., 0.] + e2 + e3
        model.add_card(card, card[0], is_list=True)


def define_coord_e123(model, Type, cid, origin, rid=0,
                      xaxis=None, yaxis=None, zaxis=None,
                      xyplane=None, yzplane=None, xzplane=None):
    """
    Create a coordinate system based on a defined axis and point on the
    plane.  This is the generalized version of the CORDx card.

    Parameters
    ----------
    model : BDF()
        a BDF object
    Type : str
        'CORD2R', 'CORD2C', 'CORD2S'
    cid : int
        the new coordinate system id
    origin : (3,) ndarray
         defines the location of the origin in the global coordinate frame
    rid : int; default=0
        the new reference coordinate system id
    xaxis : (3,) ndarray
        defines the x axis (default=None)
    yaxis : (3,) ndarray
        defines the y axis (default=None)
    zaxis : (3,) ndarray
        defines the z axis (default=None)

    .. note:: one axis (xaxis, yaxis, zaxis) and one plane
              (xyplane, yzplane, xz plane) must be defined; the others
              must be None
    .. note:: the axes and planes are defined in the rid coordinate system

    TODO: hasn't been tested...
    """
    assert Type in ['CORD2R', 'CORD2C', 'CORD2S'], Type
    origin = _fix_xyz_shape(origin, 'origin')
    rcoord = model.Coord(rid)

    # check for overdefined axes
    if xaxis is not None:
        assert yaxis is None and zaxis is None, 'yaxis=%s zaxis=%s' % (yaxis, zaxis)
        xaxis = _fix_xyz_shape(xaxis, 'xaxis')
        xaxis = rcoord.transform_node_to_global(xaxis)

    elif yaxis is not None:
        assert zaxis is None, 'zaxis=%s' % (zaxis)
        yaxis = _fix_xyz_shape(yaxis, 'yaxis')
        yaxis = rcoord.transform_node_to_global(yaxis)
    else:
        zaxis = _fix_xyz_shape(zaxis, 'zaxis')
        zaxis = rcoord.transform_node_to_global(zaxis)

    # check for invalid planes
    if xyplane is not None:
        assert yzplane is None and xzplane is None, 'yzplane=%s xzplane=%s' % (yzplane, xzplane)
        assert xaxis is not None or yaxis is not None, 'xaxis=%s yaxis=%s' % (xaxis, yaxis)
        xyplane = _fix_xyz_shape(xyplane, 'xyplane')
        xyplane = rcoord.transform_node_to_global(xyplane)
    elif yzplane is not None:
        assert xzplane is None, 'xzplane=%s' % (xzplane)
        assert yaxis is not None or zaxis is not None, 'yaxis=%s zaxis=%s' % (yaxis, zaxis)
        yzplane = _fix_xyz_shape(yzplane, 'yzplane')
        yzplane = rcoord.transform_node_to_global(yzplane)
    else:
        assert xaxis is not None or zaxis is not None, 'xaxis=%s zaxis=%s' % (xaxis, zaxis)
        xzplane = _fix_xyz_shape(xzplane, 'xzplane')
        xzplane = rcoord.transform_node_to_global(xzplane)

    if xyplane is not None:
        if xaxis is not None:
            i = xaxis / norm(xaxis)
            khat = np.cross(i, xyplane)  # xyplane is "defining" yaxis
            k = khat / norm(khat)
            j = np.cross(k, i)
        elif yaxis is not None:
            j = yaxis / norm(yaxis)
            khat = np.cross(xyplane, j)  # xyplane is "defining" xaxis
            k = khat / norm(khat)
            i = np.cross(j, k)

    elif yzplane is not None:
        if yaxis is not None:
            j = yaxis / norm(yaxis)
            ihat = np.cross(j, yzplane)  # yzplane is "defining" zaxis
            i = ihat / norm(ihat)
            k = np.cross(i, j)
        elif zaxis is not None:
            k = zaxis / norm(zaxis)
            ihat = np.cross(yzplane, zaxis)  # yzplane is "defining" yaxis
            i = ihat / norm(ihat)
            j = np.cross(k, i)

    elif xzplane is not None:
        if xaxis is not None:
            i = xaxis / norm(xaxis)
            jhat = np.cross(xzplane, i)  # xzplane is "defining" zaxis
            j = jhat / norm(jhat)
            k = np.cross(i, j)
        elif zaxis is not None:
            # standard
            k = zaxis / norm(zaxis)
            jhat = np.cross(k, xzplane) # xzplane is "defining" xaxis
            j = jhat / norm(jhat)
            i = np.cross(j, k)
    define_coord_ijk(model, Type, cid, origin, rid, i, j, k)


def define_coord_ijk(model, Type, cid, origin, rid=0, i=None, j=None, k=None):
    """
    Create a coordinate system based on 2 or 3 perpendicular unit vectors

    Parameters
    ----------
    model : BDF()
        a BDF object
    Type : str
        'CORD2R', 'CORD2C', 'CORD2S'
    cid : int
        the new coordinate system id
    origin : (3,) ndarray
         defines the location of the origin in the global coordinate frame
    rid : int; default=0
        the new reference coordinate system id
    i : (3,) ndarray
        defines the i unit vector
    j : (3,) ndarray
        defines the j unit vector
    k : (3,) ndarray
        defines the k unit vector

    TODO: hasn't been tested...
    """
    assert Type in ['CORD2R', 'CORD2C', 'CORD2S'], Type
    origin = _fix_xyz_shape(origin, 'origin')

    # create cross vectors
    if i is None:
        if j is not None and k is not None:
            i = np.cross(k, j)
        else:
            raise RuntimeError('i, j and k are None')
    else:
        # i is defined
        if j is not None and k is not None:
            # all 3 vectors are defined
            pass
        elif j is None:
            j = np.cross(k, i)
        elif k is None:
            k = np.cross(i, j)
        else:
            raise RuntimeError('j or k are None; j=%s k=%s' % (j, k))

    # define e1, e2, e3
    rcoord = model.Coord(rid, 'which is required to create cid=%s' % cid)
    e1 = rcoord.transform_node_to_local(origin)
    e2 = rcoord.transform_node_to_local(origin + k) # point on z axis
    e3 = rcoord.transform_node_to_local(origin + i) # point on x-z plane / point on x axis
    card = [Type, cid, rid] + list(e1) + list(e2) + list(e3)
    model.add_card(Type, card, is_list=True)

    coord = model.Coord(cid)
    if model.xref:
        coord.cross_reference(model)


class RectangularCoord(object):

    @classmethod
    def coord_to_xyz(cls, p):
        """
        Returns
        -------
        xyz : (3,) ndarray
            the point in the local coordinate system
        """
        return p

    @classmethod
    def xyz_to_coord(cls, p):
        """
        Returns
        -------
        xyz : (3,) ndarray
            the delta xyz point in the local coordinate system
        """
        return p

    @classmethod
    def coord_to_xyz_array(cls, p):
        """
        Returns
        -------
        xyz : (n, 3) ndarray
            the point in the local coordinate system
        """
        return p

    @classmethod
    def xyz_to_coord_array(cls, p):
        """
        Returns
        -------
        xyz : (n, 3) ndarray
            the delta xyz point in the local coordinate system
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

    .. _msc:  http://simcompanion.mscsoftware.com/resources/sites/MSC/content/meta/DOCUMENTATION/9000/DOC9188/~secure/refman.pdf?token=WDkwz5Q6v7LTw9Vb5p+nwkbZMJAxZ4rU6BoR7AHZFxi2Tl1QdrbVvWj00qmcC4+S3fnbL4WUa5ovbpBwGDBt+zFPzsGyYC13zvGPg0j/5SrMF6bnWrQoTGyJb8ho1ROYsm2OqdSA9jVceaFHQVc+tJq4b49VogM4dZBxyi/QrHgdUgPFos8BAL9mgju5WGk8yYcFtRzQIxU=
    .. seealso:: `MSC Reference Manual (pdf) <`http://simcompanion.mscsoftware.com/resources/sites/MSC/content/meta/DOCUMENTATION/9000/DOC9188/~secure/refman.pdf?token=WDkwz5Q6v7LTw9Vb5p+nwkbZMJAxZ4rU6BoR7AHZFxi2Tl1QdrbVvWj00qmcC4+S3fnbL4WUa5ovbpBwGDBt+zFPzsGyYC13zvGPg0j/5SrMF6bnWrQoTGyJb8ho1ROYsm2OqdSA9jVceaFHQVc+tJq4b49VogM4dZBxyi/QrHgdUgPFos8BAL9mgju5WGk8yYcFtRzQIxU=>`_.
    """

    @classmethod
    def coord_to_xyz(cls, p):
        r"""
        ::

          y       R
          |     /
          |   /
          | / theta
          *------------x

        .. math:: x = R \cos(\theta)
        .. math:: y = R \sin(\theta)

        Returns
        -------
        xyz : (3,) float ndarray
            the point in the local coordinate system
        """
        R = p[0]
        theta = radians(p[1])
        x = R * cos(theta)
        y = R * sin(theta)
        return np.array([x, y, p[2]], dtype='float64')

    @classmethod
    def coord_to_xyz_array(cls, p):
        r"""
        ::

          y       R
          |     /
          |   /
          | / theta
          *------------x

        .. math:: x = R \cos(\theta)
        .. math:: y = R \sin(\theta)

        Returns
        -------
        xyz : (3,) float ndarray
            the point in the local coordinate system
        """
        R = p[:, 0]
        theta = np.radians(p[:, 1])
        x = R * np.cos(theta)
        y = R * np.sin(theta)
        return np.array([x, y, p[:, 2]], dtype='float64').T

    @classmethod
    def xyz_to_coord(cls, p):
        """
        Returns
        -------
        xyz : (3,) float ndarray
            the delta xyz point in the local coordinate system
        """
        (x, y, z) = p
        theta = degrees(atan2(y, x))
        R = sqrt(x * x + y * y)
        return np.array([R, theta, z], dtype='float64')

    @classmethod
    def xyz_to_coord_array(cls, p):
        r"""
        ::

          y       R
          |     /
          |   /
          | / theta
          *------------x

        .. math:: x = R \cos(\theta)
        .. math:: y = R \sin(\theta)

        Returns
        -------
        xyz : (3,) float ndarray
            the point in the local coordinate system
        """
        x = p[:, 0]
        y = p[:, 1]
        theta = np.degrees(np.arctan(y, x))
        R = np.sqrt(x * x + y * y)
        return np.array([R, theta, p[:, 2]], dtype='float64').T


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
    @classmethod
    def xyz_to_coord(cls, p):
        r"""
        :returns xyz: the loca XYZ point in the R, \theta, \phi coordinate system
        """
        (x, y, z) = p
        R = sqrt(x * x + y * y + z * z)
        phi = degrees(atan2(y, x))
        if R > 0:
            theta = degrees(acos(z / R))
        else:
            theta = 0.
        return np.array([R, theta, phi], dtype='float64')

    @classmethod
    def xyz_to_coord_array(cls, p):
        r"""
        :returns xyz: the loca XYZ point in the R, \theta, \phi coordinate system
        """
        x = p[:, 0]
        y = p[:, 1]
        z = p[:, 2]
        R = np.sqrt(x * x + y * y + z * z)
        phi = np.degrees(np.arctan2(y, x))
        theta = np.degrees(np.arccos(z / R))

        i = np.where(R == 0.0)
        if len(i):
            theta[i] = 0.0
        return np.array([R, theta, phi], dtype='float64').T

    @classmethod
    def coord_to_xyz(cls, p):
        """
        Returns
        -------
        xyz : (3,) float ndarray
            the R, \theta, \phi in the local coordinate system
        """
        R = p[0]
        theta = radians(p[1])
        phi = radians(p[2])
        x = R * sin(theta) * cos(phi)
        y = R * sin(theta) * sin(phi)
        z = R * cos(theta)
        return np.array([x, y, z], dtype='float64')

    @classmethod
    def coord_to_xyz_array(cls, p):
        """
        Returns
        -------
        xyz : (3,) float ndarray
            the R, \theta, \phi in the local coordinate system
        """
        R = p[:, 0]
        theta = np.radians(p[:, 1])
        phi = np.radians(p[:, 2])
        x = R * np.sin(theta) * np.cos(phi)
        y = R * np.sin(theta) * np.sin(phi)
        z = R * np.cos(theta)
        return np.array([x, y, z], dtype='float64').T

class Cord2x(Coord):

    def __init__(self, cid, rid=0, origin=None, zaxis=None, xzplane=None, comment=''):
        """
        This method emulates the CORD2x card.

        Parameters
        ----------
        cid : int
            coord id
        rid : int; default=0
            reference coord id
        origin : ndarray/None; default=None -> [0., 0., 0.]
            the origin
        zaxis : ndarray/None; default=None -> [0., 0., 1.]
            a point on the z-axis
        xzplane : ndarray/None; default=None -> [1., 0., 0.]
            a point on the xz-plane

        .. note :: no type checking
        """
        Coord.__init__(self)
        if comment:
            self._comment = comment
        self.cid = cid
        self.rid = rid
        if origin is None:
            self.e1 = np.array([0., 0., 0.], dtype='float64')
        else:
            self.e1 = np.asarray(origin)

        if zaxis is None:
            self.e2 = np.array([0., 0., 1.], dtype='float64')
        else:
            self.e2 = np.asarray(zaxis)

        if xzplane is None:
            self.e3 = np.array([1., 0., 0.], dtype='float64')
        else:
            self.e3 = np.asarray(xzplane)
        self._finish_setup()

    @classmethod
    def _add(cls, cid, rid=0, origin=None, zaxis=None, xzplane=None, comment=''):
        cid = cid
        rid = rid

        if origin is None:
            e1 = np.array([0., 0., 0.], dtype='float64')
        else:
            e1 = np.asarray(origin)

        if zaxis is None:
            e2 = np.array([0., 0., 1.], dtype='float64')
        else:
            e2 = np.asarray(zaxis)

        if xzplane is None:
            e3 = np.array([1., 0., 0.], dtype='float64')
        else:
            e3 = np.asarray(xzplane)
        return cls(cid, rid, e1, e2, e3, comment=comment)
        #self._finish_setup()

    @classmethod
    def add_axes(cls, cid, rid=0, origin=None,
                 xaxis=None, yaxis=None, zaxis=None,
                 xyplane=None, yzplane=None, xzplane=None,
                 comment=''):
        """
        Create a coordinate system based on a defined axis and point on the
        plane.  This is the generalized version of the CORD2x card.

        Parameters
        ----------
        cid : int
            the new coordinate system id
        rid : int; default=0
            the new reference coordinate system id
        origin : (3,) ndarray
             defines the location of the origin in the global coordinate frame
        xaxis : (3,) ndarray
            defines the x axis (default=None)
        yaxis : (3,) ndarray
            defines the y axis (default=None)
        zaxis : (3,) ndarray
            defines the z axis (default=None)

        .. note:: one axis (xaxis, yaxis, zaxis) and one plane
                  (xyplane, yzplane, xz plane) must be defined; the others
                  must be None
        .. note:: the axes and planes are defined in the rid coordinate system
        """
        assert cls.type in ['CORD2R', 'CORD2C', 'CORD2S'], cls.type
        if origin is None:
            origin = np.array([0., 0., 0.], dtype='float64')
        else:
            origin = _fix_xyz_shape(origin, 'origin')

        # check for overdefined axes
        if xaxis is not None:
            assert yaxis is None and zaxis is None, 'yaxis=%s zaxis=%s' % (yaxis, zaxis)
            xaxis = _fix_xyz_shape(xaxis, 'xaxis')
            xaxis = cls.coord_to_xyz(xaxis)
        elif yaxis is not None:
            assert zaxis is None, 'zaxis=%s' % (zaxis)
            yaxis = _fix_xyz_shape(yaxis, 'yaxis')
            yaxis = cls.coord_to_xyz(yaxis)
        else:
            zaxis = _fix_xyz_shape(zaxis, 'zaxis')
            zaxis = cls.coord_to_xyz(zaxis)

        # check for invalid planes
        if xyplane is not None:
            assert yzplane is None and xzplane is None, 'yzplane=%s xzplane=%s' % (yzplane, xzplane)
            assert xaxis is not None or yaxis is not None, 'xaxis=%s yaxis=%s' % (xaxis, yaxis)
            xyplane = _fix_xyz_shape(xyplane, 'xyplane')
            xyplane = cls.coord_to_xyz(xyplane)
        elif yzplane is not None:
            assert xzplane is None, 'xzplane=%s' % (xzplane)
            assert yaxis is not None or zaxis is not None, 'yaxis=%s zaxis=%s' % (yaxis, zaxis)
            yzplane = _fix_xyz_shape(yzplane, 'yzplane')
            yzplane = cls.coord_to_xyz(yzplane)
        else:
            assert xaxis is not None or zaxis is not None, 'xaxis=%s zaxis=%s' % (xaxis, zaxis)
            xzplane = _fix_xyz_shape(xzplane, 'xzplane')
            xzplane = cls.coord_to_xyz(xzplane)

        if xyplane is not None:
            if xaxis is not None:
                i = xaxis / norm(xaxis)
                khat = np.cross(i, xyplane)  # xyplane is "defining" yaxis
                k = khat / norm(khat)
                j = np.cross(k, i)
            elif yaxis is not None:
                j = yaxis / norm(yaxis)
                khat = np.cross(xyplane, j)  # xyplane is "defining" xaxis
                k = khat / norm(khat)
                i = np.cross(j, k)

        elif yzplane is not None:
            if yaxis is not None:
                j = yaxis / norm(yaxis)
                ihat = np.cross(j, yzplane)  # yzplane is "defining" zaxis
                i = ihat / norm(ihat)
                k = np.cross(i, j)
            elif zaxis is not None:
                k = zaxis / norm(zaxis)
                ihat = np.cross(yzplane, zaxis)  # yzplane is "defining" yaxis
                i = ihat / norm(ihat)
                j = np.cross(k, i)

        elif xzplane is not None:
            if xaxis is not None:
                i = xaxis / norm(xaxis)
                jhat = np.cross(xzplane, i)  # xzplane is "defining" zaxis
                j = jhat / norm(jhat)
                k = np.cross(i, j)
            elif zaxis is not None:
                # standard
                k = zaxis / norm(zaxis)
                jhat = np.cross(k, xzplane) # xzplane is "defining" xaxis
                j = jhat / norm(jhat)
                i = np.cross(j, k)
        return cls.add_ijk(cid, rid, origin, i, j, k, comment=comment)

    @classmethod
    def add_ijk(cls, cid, rid=0, origin=None, i=None, j=None, k=None, comment=''):
        """
        Create a coordinate system based on 2 or 3 perpendicular unit vectors

        Parameters
        ----------
        cid : int
            the new coordinate system id
        origin : (3,) ndarray
             defines the location of the origin in the global coordinate frame
        rid : int; default=0
            the new reference coordinate system id
        i : (3,) ndarray
            defines the i unit vector
        j : (3,) ndarray
            defines the j unit vector
        k : (3,) ndarray
            defines the k unit vector
        """
        Type = cls.type
        assert Type in ['CORD2R', 'CORD2C', 'CORD2S'], Type
        if origin is None:
            origin = np.array([0., 0., 0.], dtype='float64')
        else:
            origin = _fix_xyz_shape(origin, 'origin')

        # create cross vectors
        if i is None:
            if j is not None and k is not None:
                i = np.cross(k, j)
            else:
                raise RuntimeError('i, j and k are None')
        else:
            # i is defined
            if j is not None and k is not None:
                # all 3 vectors are defined
                pass
            elif j is None:
                j = np.cross(k, i)
            elif k is None:
                k = np.cross(i, j)
            else:
                raise RuntimeError('j or k are None; j=%s k=%s' % (j, k))

        # origin
        e1 = origin
        # point on z axis
        e2 = origin + k

        # point on x-z plane / point on x axis
        e3 = origin + i
        return cls(cid, rid, e1, e2, e3, comment=comment)

    @classmethod
    def add_op2_data(cls, data, comment=''):
        cid = data[0]
        rid = data[1]
        e1 = np.array(data[2:5], dtype='float64')
        e2 = np.array(data[5:8], dtype='float64')
        e3 = np.array(data[8:11], dtype='float64')
        assert len(data) == 11, 'data = %s' % (data)
        return cls(cid, rid, e1, e2, e3, comment=comment)

    @classmethod
    def add_card(cls, card, comment=''):
        """
        Defines the CORD2x class
        """
        #: coordinate system ID
        cid = integer(card, 1, 'cid')
        #: reference coordinate system ID
        rid = integer_or_blank(card, 2, 'rid', 0)

        #: origin in a point relative to the rid coordinate system
        origin = np.array([double_or_blank(card, 3, 'e1x', 0.0),
                           double_or_blank(card, 4, 'e1y', 0.0),
                           double_or_blank(card, 5, 'e1z', 0.0)],
                          dtype='float64')
        #: z-axis in a point relative to the rid coordinate system
        zaxis = np.array([double_or_blank(card, 6, 'e2x', 0.0),
                          double_or_blank(card, 7, 'e2y', 0.0),
                          double_or_blank(card, 8, 'e2z', 0.0)],
                         dtype='float64')
        #: a point on the xz-plane relative to the rid coordinate system
        xzplane = np.array([double_or_blank(card, 9, 'e3x', 0.0),
                            double_or_blank(card, 10, 'e3y', 0.0),
                            double_or_blank(card, 11, 'e3z', 0.0)],
                           dtype='float64')
        return cls(cid, rid, origin, zaxis, xzplane, comment=comment)
        #self._finish_setup()

    def _finish_setup(self):
        assert len(self.e1) == 3, self.e1
        assert len(self.e2) == 3, self.e2
        assert len(self.e3) == 3, self.e3

        #: the global axes
        self.i = None
        self.j = None
        self.k = None

        if self.rid == 0:
            self.isResolved = True
            self.setup()
            #self.setup_global_cord2x()

    def _verify(self, xref):
        """
        Verifies all methods for this object work

        Parameters
        ----------
        xref : bool
            has this model been cross referenced
        """
        cid = self.Cid()
        rid = self.Rid()
        assert isinstance(cid, int), 'cid=%r' % cid
        assert isinstance(rid, int), 'rid=%r' % rid

    def write_card(self, size=8, is_double=False):
        card = self.repr_fields()
        if size == 8:
            return self.comment + print_card_8(card)
        elif is_double:
            return self.comment + print_card_double(card)
        return self.comment + print_card_16(card)

    def cross_reference(self, model):
        """
        Cross links the card so referenced cards can be extracted directly

        Parameters
        ----------
        model : BDF()
            the BDF object

        .. warning:: Doesn't set rid to the coordinate system if it's in the
                    global.  This isn't a problem.  It's meant to speed up the
                    code in order to resolve extra coordinate systems.
        """
        if self.Rid() != 0:
            msg = ' which is required by %s cid=%s' % (self.type, self.cid)
            self.rid = model.Coord(self.rid, msg=msg)
            self.rid_ref = self.rid

    def uncross_reference(self):
        if self.rid == 0:
            return
        self.rid = self.Rid()
        del self.rid_ref

    def Rid(self):
        """Gets the reference coordinate system self.rid"""
        if isinstance(self.rid, integer_types):
            return self.rid
        return self.rid_ref.cid


class Cord1x(Coord):
    rid = 0  # used only for transform to global

    def Rid(self):
        """Gets the reference coordinate system self.rid"""
        return self.rid

    def __init__(self, cid, g1, g2, g3, comment=''):
        Coord.__init__(self)
        if comment:
            self._comment = comment

        #: the coordinate ID
        self.cid = cid
        #: a Node at the origin
        self.g1 = g1
        #: a Node on the z-axis
        self.g2 = g2
        #: a Node on the xz-plane
        self.g3 = g3

    def validate(self):
        assert self.g1 != self.g2, str(self)
        assert self.g1 != self.g3, str(self)
        assert self.g2 != self.g3, str(self)

    @classmethod
    def add_card(cls, card, icard=0, comment=''):
        #self.isResolved = False
        assert icard in (0, 1), 'icard=%r' % (icard)
        ncoord = 4 * icard  # 0 if the 1st coord, 4 if the 2nd

        cid = integer(card, 1 + ncoord, 'cid')
        g1 = integer(card, 2 + ncoord, 'g1')
        g2 = integer(card, 3 + ncoord, 'g2')
        g3 = integer(card, 4 + ncoord, 'g3')

        #self.e1 = None
        #self.e2 = None
        #self.e3 = None
        #self.i = None
        #self.j = None
        #self.k = None
        return cls(cid, g1, g2, g3, comment=comment)

    @classmethod
    def add_op2_data(cls, data, comment=''):
        cid = data[0]
        g1 = data[1]
        g2 = data[2]
        g3 = data[3]
        assert len(data) == 4, 'data = %s' % (data)
        print('cid =', cid)
        return cls(cid, g1, g2, g3, comment=comment)

    def to_CORD2x(self, model, rid=0):
        """
        Converts a coordinate system from a CORD1x to a CORD2x

        Parameters
        ----------
        model : BDF()
            a BDF model
        rid : int; default=0
            The relative coordinate system
        """
        rid1 = self.g1.Cid()
        rid2 = self.g2.Cid()
        rid3 = self.g2.Cid()

        # assume the points are in rid
        p1 = self.g1_ref.xyz
        p2 = self.g2_ref.xyz
        p3 = self.g3_ref.xyz

        # move the nodes in necessary into rid system
        if rid != rid1:
            p1 = self.g1_ref.get_position_wrt(model, rid)
        if rid != rid2:
            p2 = self.g2_ref.get_position_wrt(model, rid)
        if rid != rid3:
            p3 = self.g3_ref.get_position_wrt(model, rid)

        type1 = self.type.replace('1', '2')
        data = [type1, self.cid, rid1, list(p1) + list(p2) + list(p3)]

        if self.type == 'CORD1R':
            coord = CORD2R.add_op2_data(data=data, comment=self.comment)
        elif self.type == 'CORD1C':
            coord = CORD2C.add_op2_data(data=data, comment=self.comment)
        elif self.type == 'CORD1S':
            coord = CORD2S.add_op2_data(data=data, comment=self.comment)
        else:
            raise RuntimeError('coordinate type of \n%s is %s' % (str(self), type1))
        model.coords[self.cid] = coord

    def _verify(self, xref):
        """
        Verifies all methods for this object work

        Parameters
        ----------
        xref : bool
            has this model been cross referenced
        """
        cid = self.Cid()
        assert isinstance(cid, int), 'cid=%r' % cid

    def cross_reference(self, model):
        """
        Cross links the card so referenced cards can be extracted directly

        Parameters
        ----------
        model : BDF()
            the BDF object
        """
        msg = ' which is required by %s cid=%s' % (self.type, self.cid)
        #: grid point 1
        self.g1 = model.Node(self.g1, msg=msg)
        self.g1_ref = self.g1
        #: grid point 2
        self.g2 = model.Node(self.g2, msg=msg)
        self.g2_ref = self.g2
        #: grid point 3
        self.g3 = model.Node(self.g3, msg=msg)
        self.g3_ref = self.g3

    def uncross_reference(self):
        self.g1 = self.G1()
        self.g2 = self.G2()
        self.g3 = self.G3()
        del self.g1_ref, self.g2_ref, self.g3_ref

    def setup(self):
        """
        Finds the position of the nodes used define the coordinate system
        and sets the ijk vectors
        """
        if self.isResolved:
            return

        self.g1_ref.cp_ref.setup()
        self.g2_ref.cp_ref.setup()
        self.g3_ref.cp_ref.setup()

        #: the origin in the local frame
        self.e1 = self.g1_ref.get_position()

        #: a point on the z-axis
        self.e2 = self.g2_ref.get_position()

        #: a point on the xz-plane
        self.e3 = self.g3_ref.get_position()

        # rid is resolved b/c e1, e2, & e3 are in global coordinates
        self.isResolved = False

        #print('setting up cid=%s' % self.cid)
        # call the Coord class' setup method
        super(Cord1x, self).setup()

    def G1(self):
        if isinstance(self.g1, integer_types):
            return self.g1
        return self.g1_ref.nid

    def G2(self):
        if isinstance(self.g2, integer_types):
            return self.g2
        return self.g2_ref.nid

    def G3(self):
        if isinstance(self.g3, integer_types):
            return self.g3
        return self.g3_ref.nid

    def nodeIDs(self):
        """
        Gets the integers for the node [g1,g2,g3]
        """
        self.deprecated('self.nodeIDs()', 'self.node_ids', '0.8')
        return self.node_ids

    @property
    def node_ids(self):
        grids = [self.G1(), self.G2(), self.G3()]
        return grids

    def write_card(self, size=8, is_double=False):
        card = self.repr_fields()
        return self.comment + print_card_8(card)


class GMCORD(BaseCard):
    type = 'GMCORD'

    def __init__(self, cid, entity, gm_ids, comment=''):
        if comment:
            self._comment = comment
        self.cid = cid
        self.entity = entity
        self.gm_ids = gm_ids

    @classmethod
    def add_card(cls, card, comment=''):
        cid = integer(card, 1, 'cid')
        entity = string(card, 2, 'entity')
        gm_ids = [
            integer(card, 3, 'GM_ID1'),
            integer_or_blank(card, 4, 'GM_ID2'),
        ]
        return GMCORD(cid, entity, gm_ids, comment=comment)

    def cross_reference(self, model):
        pass

    def setup(self):
        pass

    def raw_fields(self):
        list_fields = ['GMCORD', self.cid, self.entity] + self.gm_ids
        return list_fields

    def write_card(self, size=8, is_double=False):
        card = self.repr_fields()
        return self.comment + print_card_8(card)

class CORD3G(Coord):  # not done
    """
    Defines a general coordinate system using three rotational angles as
    functions of coordinate values in the reference coordinate system.
    The CORD3G entry is used with the MAT9 entry to orient material principal
    axes for 3-D composite analysis.

    +--------+-----+--------+------+----------+----------+----------+--------+
    | CORD3G | CID | METHOD | FORM | THETAID1 | THETAID2 | THETAID3 | CIDREF |
    +--------+-----+--------+------+----------+----------+----------+--------+
    | CORD3G | 100 |  E313  | EQN  |    110   |    111   |    112   |   0    |
    +--------+-----+--------+------+----------+----------+----------+--------+
    """
    type = 'CORD3G'

    def __init__(self, cid, method_es, method_int, form, thetas, rid, comment=''):
        """
        Intilizes the CORD3G

        :param card: a list version of the fields
        """
        Coord.__init__(self, card, comment)
        self.cid = cid
        self.method_es = method_es
        self.method_int = method_int
        self.form = form
        self.thetas = thetas
        self.rid = rid

        assert 0 < self.method_int < 1000
        assert len(self.thetas) == 3, 'thetas=%s' % (self.thetas)

        # EQN for DEQATN, TABLE for TABLE3D
        assert self.form in ['EQN', 'TABLE']

    @classmethod
    def add_op2_data(cls, data, comment=''):
        raise NotImplementedError(data)

    @classmethod
    def add_card(cls, card, comment=''):
        self.cid = integer(card, 1, 'cid')
        method = string_or_blank(card, 2, 'E313')
        method_es = method[0]
        method_int = int(method[1:])
        assert self.methodES in ['E', 'S'] # Euler / Space-Fixed

        form = string_or_blank(card, 3, 'form', 'EQN')
        thetas = [integer(card, 4, 'theta1'),
                  integer(card, 5, 'theta2'),
                  integer(card, 6, 'theta3')]
        rid = integer_or_blank(card, 7, 'cidRef')
        assert len(card) <= 8, 'len(CORD3G card) = %i\ncard=%s' % (len(card), card)

        return CORD3G(cid, method_es, method_int, form, thetas, rid,
                      comment=comment)

    def cross_reference(self, model):
        """
        Cross links the card so referenced cards can be extracted directly

        Parameters
        ----------
        model : BDF()
            the BDF object
        """
        msg = ' which is required by %s cid=%s' % (self.type, self.cid)
        self.rid = model.Coord(self.rid, msg=msg)
        self.rid_ref = self.rid

    def uncross_reference(self):
        self.rid = self.Rid()
        del self.rid_ref

    def Rid(self):
        if isinstance(self.rid, integer_types):
            return self.rid
        return self.rid_ref.cid

    def coord3g_transform_to_global(self, p):
        """
        :param p:     the point to transform.  TYPE=NUMPY.NDARRAY.

        .. warning:: not done, just setting up how you'd do this
        .. note::    per http://en.wikipedia.org/wiki/Euler_angles
         "This means for example that a convention named (YXZ) is the result
         of performing first an intrinsic Z rotation, followed by X and
         Y rotations, in the moving axes (Note: the order of multiplication
         of matrices is the opposite of the order in which they're
         applied to a vector)."
        """
        if self.method_es == 'E':
            rotations = [int(i) for i in str(self.method_int)]
            for (rotation, theta) in zip(rotations, self.thetas):
                ct = cos(radians(theta))
                st = sin(radians(theta))
                if rotation == 1:
                    p = np.dot(self.rotation_x(ct, st), p)
                elif rotation == 2:
                    p = np.dot(self.rotation_y(ct, st), p)
                elif rotation == 3:
                    p = np.dot(self.rotation_z(ct, st), p)
                else:
                    raise RuntimeError('rotation=%s rotations=%s' % (rotation, rotations))
        elif self.method_es == 'S':
            raise RuntimeError('Space-Fixed rotation hasnt been implemented')
        else:
            msg = 'Invalid method; Use Euler or Space-Fixed.  method_es=%r' % self.method_es
            raise RuntimeError(msg)
        return p

    def rotation_x(self, ct, st):
        matrix = np.array([[1., 0., 0.],
                           [ct, 0., -st],
                           [-st, 0., ct]])
        return matrix

    def rotation_y(self, ct, st):
        matrix = np.array([[ct, 0., st],
                           [0., 1., 0.],
                           [-st, 0., ct]])
        return matrix

    def rotation_z(self, ct, st):
        matrix = np.array([[ct, st, 0.],
                           [-st, ct, 0.],
                           [0., 0., 1.]])
        return matrix

    def raw_fields(self):
        method = self.method_es + str(self.method_int)
        list_fields = (['CORD3G', self.cid, method, self.form] + self.thetas +
                       [self.Rid()])
        return list_fields


class CORD1R(Cord1x, RectangularCoord):
    type = 'CORD1R'
    Type = 'R'
    int_type = 0

    def __init__(self, cid, g1, g2, g3, comment=''):
        """
        Intilizes the CORD1R

        +-------+------+-----+-----+------+------+-----+------+-----+
        |   1   |   2  |  3  |  4  |   5  |  6   |  7  |  8   |  9  |
        +=======+======+=====+=====+======+======+=====+======+=====+
        |CORD1R | CIDA | G1A | G2A | CIDB | G1B  | G2B | G3B  |     |
        +-------+------+-----+-----+------+------+-----+------+-----+

        Parameters
        ----------
        self :CORD1R()
            the CORD1R coordinate system object
        icard : int
            the coordinate location on the line
            (there are possibly 2 coordinates on 1 card)
        card : List
            a list version of the fields (1 CORD1R only)
        """
        Cord1x.__init__(self, cid, g1, g2, g3, comment=comment)

    def raw_fields(self):
        list_fields = ['CORD1R', self.cid] + self.node_ids
        return list_fields


class CORD1C(Cord1x, CylindricalCoord):
    type = 'CORD1C'
    Type = 'C'

    def __init__(self, cid, g1, g2, g3, comment=''):
        """
        Intilizes the CORD1R

        +-------+------+-----+-----+------+------+-----+------+-----+
        |   1   |   2  |  3  |  4  |   5  |  6   |  7  |  8   |  9  |
        +=======+======+=====+=====+======+======+=====+======+=====+
        |CORD1C | CIDA | G1A | G2A | CIDB | G1B  | G2B | G3B  |     |
        +-------+------+-----+-----+------+------+-----+------+-----+

        :param icard:  the coordinate location on the line
                       (there are possibly 2 coordinates on 1 card)
        :param data:   a list version of the fields (1 CORD1R only)
        """
        Cord1x.__init__(self, cid, g1, g2, g3, comment=comment)

    #@classmethod
    #def add_op2_data(cls, data, comment):
        #self.cid = data[0]
        #self.g1 = data[1]
        #self.g2 = data[2]
        #self.g3 = data[3]
        #assert len(data) == 4, 'data = %s' % (data)
        #bbbb

    def raw_fields(self):
        list_fields = ['CORD1C', self.cid] + self.node_ids
        return list_fields


class CORD1S(Cord1x, SphericalCoord):
    type = 'CORD1S'
    Type = 'S'

    def __init__(self, cid, g1, g2, g3, comment=''):
        """
        Intilizes the CORD1S

        +-------+------+-----+-----+------+------+-----+------+-----+
        |   1   |   2  |  3  |  4  |   5  |  6   |  7  |  8   |  9  |
        +=======+======+=====+=====+======+======+=====+======+=====+
        |CORD1S | CIDA | G1A | G2A | CIDB | G1B  | G2B | G3B  |     |
        +-------+------+-----+-----+------+------+-----+------+-----+

        :param card:   a BDFCard object
        :param icard: the coordinate location on the line
                       (there are possibly 2 coordinates on 1 card)
        :param data:   a list version of the fields (1 CORD1S only)
        """
        Cord1x.__init__(self, cid, g1, g2, g3, comment=comment)

    def raw_fields(self):
        list_fields = ['CORD1S', self.cid] + self.node_ids
        return list_fields


class CORD2R(Cord2x, RectangularCoord):
    type = 'CORD2R'
    Type = 'R'

    def __init__(self, cid, rid=0, origin=None, zaxis=None, xzplane=None, comment=''):
        """
        Intilizes the CORD2R

        +--------+-----+-----+-----+----+-----+----+----+-----+
        |    1   |   2 |  3  |  4  |  5 |  6  |  7 |  8 |  9  |
        +========+=====+=====+=====+====+=====+====+====+=====+
        | CORD2R | CID | RID | A1  | A2 | A3  | B1 | B2 |     |
        +--------+-----+-----+-----+----+-----+----+----+-----+
        |        | B3  | C1  | C2  | C3 |     |    |    |     |
        +--------+-----+-----+-----+----+-----+----+----+-----+

        :param card: a BDFCard object
        :param data: a list version of the fields (1 CORD2R only)
                     default=None -> [0, 0, 0., 0., 0., 0., 0., 1., 1., 0., 0.]

        .. note :: no type checking
        """
        Cord2x.__init__(self, cid, rid, origin, zaxis, xzplane, comment=comment)

    def _verify(self, xref):
        """
        Verifies all methods for this object work

        :param xref: has this model been cross referenced
        :type xref:  bool
        """
        cid = self.Cid()
        rid = self.Rid()
        assert isinstance(cid, int), 'cid=%r' % cid
        assert isinstance(rid, int), 'rid=%r' % rid

    def raw_fields(self):
        rid = set_blank_if_default(self.Rid(), 0)
        list_fields = ['CORD2R', self.cid, rid] + list(self.e1) + list(
            self.e2) + list(self.e3)
        return list_fields


class CORD2C(Cord2x, CylindricalCoord):
    type = 'CORD2C'
    Type = 'C'

    def __init__(self, cid, rid=0, origin=None, zaxis=None, xzplane=None, comment=''):
        """
        Intilizes the CORD2C

        +--------+-----+-----+-----+----+-----+----+----+-----+
        |    1   |   2 |  3  |  4  |  5 |  6  |  7 |  8 |  9  |
        +========+=====+=====+=====+====+=====+====+====+=====+
        | CORD2C | CID | RID | A1  | A2 | A3  | B1 | B2 |     |
        +--------+-----+-----+-----+----+-----+----+----+-----+
        |        | B3  | C1  | C2  | C3 |     |    |    |     |
        +--------+-----+-----+-----+----+-----+----+----+-----+

        :param card: a BDFCard object
        :param data: a list version of the fields (1 CORD2C only)
        """
        Cord2x.__init__(self, cid, rid, origin, zaxis, xzplane, comment=comment)

    def raw_fields(self):
        rid = set_blank_if_default(self.Rid(), 0)
        list_fields = (['CORD2C', self.cid, rid] + list(self.e1) +
                       list(self.e2) + list(self.e3))
        return list_fields


class CORD2S(Cord2x, SphericalCoord):
    type = 'CORD2S'
    Type = 'S'

    def __init__(self, cid, rid=0, origin=None, zaxis=None, xzplane=None, comment=''):
        """
        Intilizes the CORD2S

        +--------+-----+-----+-----+----+-----+----+----+-----+
        |    1   |   2 |  3  |  4  |  5 |  6  |  7 |  8 |  9  |
        +========+=====+=====+=====+====+=====+====+====+=====+
        | CORD2S | CID | RID | A1  | A2 | A3  | B1 | B2 |     |
        +--------+-----+-----+-----+----+-----+----+----+-----+
        |        | B3  | C1  | C2  | C3 |     |    |    |     |
        +--------+-----+-----+-----+----+-----+----+----+-----+

        :param card: a BDFCard object
        :param data: a list version of the fields (1 CORD2S only)
        """
        Cord2x.__init__(self, cid, rid, origin, zaxis, xzplane, comment=comment)

    def raw_fields(self):
        rid = set_blank_if_default(self.Rid(), 0)
        list_fields = (['CORD2S', self.cid, rid] + list(self.e1) +
                       list(self.e2) + list(self.e3))
        return list_fields
