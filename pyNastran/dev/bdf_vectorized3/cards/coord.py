from __future__ import annotations
from itertools import count
from typing import Any, TYPE_CHECKING
import numpy as np
#from pyNastran.bdf.field_writer_8 import print_float_8 # print_card_8
#from pyNastran.bdf.field_writer_16 import print_card_16 #, print_scientific_16, print_field_16
#from pyNastran.bdf.field_writer_double import print_scientific_double
from pyNastran.utils.numpy_utils import integer_types
from pyNastran.bdf.bdf_interface.assign_type import (
    integer, integer_or_blank, double_or_blank,
)
from pyNastran.femutils.coord_transforms import (
    #xyz_to_rtz_array,
    #xyz_to_rtp_array, # xyz to xxx transforms
    rtz_to_xyz_array as cylindrical_to_rectangular,
    rtp_to_xyz_array as spherical_to_rectangular, # xxx to xyz transforms
    #rtz_to_rtp_array, rtp_to_rtz_array, # rtp/rtz and rtz/rtp transforms
)

from pyNastran.bdf.bdf_interface.bdf_card import BDFCard
from pyNastran.bdf.cards.coordinate_systems import _fix_xyz_shape
from pyNastran.dev.bdf_vectorized3.cards.base_card import VectorizedBaseCard
from pyNastran.dev.bdf_vectorized3.cards.write_utils import (
    get_print_card_size, array_float, array_str, array_default_int)
if TYPE_CHECKING:  # pragma: no cover
    from pyNastran.dev.bdf_vectorized3.bdf import BDF
    from pyNastran.dev.bdf_vectorized3.types import TextIOLike
    from .grid import GRID

class COORD(VectorizedBaseCard):
    """
    TODO: cleanup API

    methods:
        transform_node_to_*_xyz (local/global)
        transform_node_to_*_coord_id (local)
        transform_local_xyz_to_global
        transform_global_xyz_to_local_coord_id
    """
    _id_name = 'coord_id'
    def __init__(self, model: BDF):
        super().__init__(model)
        self.cards1: list = []
        self.cards2: list = []
        self.coord_type = np.array(['R'], dtype='|U1')
        self.icoord = np.array([2], dtype='int8')
        self.is_resolved = np.array([True], dtype='bool')
        self.nodes = np.full((1, 3), -1, dtype='int32')

        self.e1 = np.array([[0., 0., 0.]], dtype='float64')
        self.e2 = np.array([[0., 0., 1.]], dtype='float64')
        self.e3 = np.array([[1., 0., 0.]], dtype='float64')
        self.origin = np.array([[0., 0., 0.]], dtype='float64')
        self.i = np.array([[1., 0., 0.]], dtype='float64')
        self.j = np.array([[0., 1., 0.]], dtype='float64')
        self.k = np.array([[0., 0., 1.]], dtype='float64')
        ones = np.eye(3, dtype='float64')
        self.T = np.array([ones], dtype='float64')
        self.xyz_to_global_transform = {0: ones}
        assert ones.shape == (3, 3), ones
        assert self.T.shape == (1, 3, 3), self.T.shape

        self.coord_id = np.array([0])
        self.ref_coord_id = np.array([0])
        self.origin = np.zeros((1, 3), dtype='float64')
        self.z_axis = np.zeros((1, 3), dtype='float64')
        self.xz_plane = np.zeros((1, 3), dtype='float64')
        #self.origin[0, :] = [0., 0., 0.]
        #self.z_axis[0, :] = [0., 0., 1.]
        #self.xz_plane[0, :] = [1., 0., 0.]
        self.n = 1
        self.setup()

    @classmethod
    def slice_card_by_coord_id(cls, coord_id: np.ndarray) -> COORD:
        i = cls.index(coord_id)
        coord = cls.slice_card_by_index(i)
        return coord

    def set_used(self, used_dict: dict[str, list[np.ndarray]]) -> None:
        is_cord1 = (self.icoord == 1)
        is_cord2 = ~is_cord1
        if is_cord1.sum():
            used_dict['node_id'].append(self.nodes[is_cord1, :].ravel())
        if is_cord2.sum():
            used_dict['coord_id'].append(self.ref_coord_id)

    def remove_unused(self, used_dict: dict[str, np.ndarray]) -> int:
        coord_id = used_dict['coord_id']
        ncards_removed = len(self.coord_id) - len(coord_id)
        if ncards_removed:
            self.slice_card_by_id(coord_id, assume_sorted=True, sort_ids=False)
        return ncards_removed

    def convert(self, xyz_scale: float=1.0, **kwargs) -> None:
        is_xyz = (self.icoord == 2) & (self.coord_type == 'R')
        is_rtz = (self.icoord == 2) & (self.coord_type == 'C')
        is_rtp = (self.icoord == 2) & (self.coord_type == 'S')
        for res in [self.origin, self.z_axis, self.xz_plane, self.e1, self.e2, self.e3]:
            res[is_xyz, :] *= xyz_scale # xyz
            res[is_rtz, 0] *= xyz_scale # R
            res[is_rtz, 2] *= xyz_scale # z
            res[is_rtp, 0] *= xyz_scale # rho

    def __apply_slice__(self, coord: COORD, i: np.ndarray) -> None:
        # general
        coord.n = len(i)
        coord.icoord = self.icoord[i]
        coord.coord_type = self.coord_type[i]
        coord.coord_id = self.coord_id[i]

        # CORD1
        coord.nodes = self.nodes[i, :]

        # CORD2
        coord.ref_coord_id = self.ref_coord_id[i]
        coord.e1 = self.e1[i, :]
        coord.e2 = self.e2[i, :]
        coord.e3 = self.e3[i, :]

        # general
        coord.origin = self.origin[i, :]
        coord.z_axis = self.z_axis[i, :]
        coord.xz_plane = self.xz_plane[i, :]
        coord.i = self.i[i, :]
        coord.j = self.j[i, :]
        coord.k = self.k[i, :]
        coord.T = self.T[i, :, :]
        coord.is_resolved = self.is_resolved[i]
        coord.xyz_to_global_transform = {cid: self.xyz_to_global_transform[cid]
                                         for cid in coord.coord_id}
        #assert np.all(np.isfinite(coord.e1)), coord.e1
        #assert np.all(np.isfinite(coord.origin)), coord.origin
        #assert np.all(np.isfinite(coord.T)), coord.T

    def _set_global_coord(self):
        """being lazy with setting coordinate systems with NANs means we have
        to fix the global coordinate system..."""
        icoords = np.where(self.coord_id == 0)[0]
        for icoord in icoords:
            self.e1[icoord, :] = [0., 0., 0.]
            self.e2[icoord, :] = [0., 0., 1.]
            self.e3[icoord, :] = [1., 0., 0.]

            self.origin[icoord, :] = [0., 0., 0.]
            self.z_axis[icoord, :] = [0., 0., 1.]
            self.xz_plane[icoord, :] = [1., 0., 0.]

            self.i[icoord, :] = [1., 0., 0.]
            self.j[icoord, :] = [0., 1., 0.]
            self.k[icoord, :] = [0., 0., 1.]
            self.T[icoord, :] = np.eye(3)
            self.is_resolved[icoord] = True
            self.xyz_to_global_transform[0] = self.T[icoord, :]

    def add_axes(self, cid: int, coord_type: str, rid: int=0, origin=None,
                 xaxis=None, yaxis=None, zaxis=None,
                 xyplane=None, yzplane=None, xzplane=None,
                 comment: str=''):
        """
        Create a coordinate system based on a defined axis and point on the
        plane.  This is the generalized version of the CORD2x card.

        Parameters
        ----------
        cid : int
            the new coordinate system id
        coord_type: str
            R, C, S
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

        Notes
        -----
        One axis (xaxis, yaxis, zaxis) and one plane
        (xyplane, yzplane, xz plane) must be defined; the others
        must be None

        The axes and planes are defined in the rid coordinate system

        """
        assert coord_type in ['R', 'C', 'S'], coord_type
        if origin is None:
            origin = np.array([0., 0., 0.], dtype='float64')
        else:
            origin = _fix_xyz_shape(origin, 'origin')

        # check for overdefined axes
        if xaxis is not None:
            assert yaxis is None and zaxis is None, 'yaxis=%s zaxis=%s' % (yaxis, zaxis)
            xaxis = _fix_xyz_shape(xaxis, 'xaxis')
            xaxis = coord_to_xyz(xaxis, coord_type)
        elif yaxis is not None:
            assert zaxis is None, 'zaxis=%s' % (zaxis)
            yaxis = _fix_xyz_shape(yaxis, 'yaxis')
            yaxis = coord_to_xyz(yaxis, coord_type)
        else:
            zaxis = _fix_xyz_shape(zaxis, 'zaxis')
            zaxis = coord_to_xyz(zaxis, coord_type)

        # check for invalid planes
        if xyplane is not None:
            assert yzplane is None and xzplane is None, 'yzplane=%s xzplane=%s' % (yzplane, xzplane)
            assert xaxis is not None or yaxis is not None, 'xaxis=%s yaxis=%s' % (xaxis, yaxis)
            xyplane = _fix_xyz_shape(xyplane, 'xyplane')
            xyplane = coord_to_xyz(xyplane, coord_type)
        elif yzplane is not None:
            assert xzplane is None, 'xzplane=%s' % (xzplane)
            assert yaxis is not None or zaxis is not None, 'yaxis=%s zaxis=%s' % (yaxis, zaxis)
            yzplane = _fix_xyz_shape(yzplane, 'yzplane')
            yzplane = coord_to_xyz(yzplane, coord_type)
        else:
            assert xaxis is not None or zaxis is not None, 'xaxis=%s zaxis=%s' % (xaxis, zaxis)
            xzplane = _fix_xyz_shape(xzplane, 'xzplane')
            xzplane = coord_to_xyz(xzplane, coord_type)

        if xyplane is not None:
            if xaxis is not None:
                i = xaxis / np.linalg.norm(xaxis)
                khat = np.cross(i, xyplane)  # xyplane is "defining" yaxis
                k = khat / np.linalg.norm(khat)
                j = np.cross(k, i)
            elif yaxis is not None:
                j = yaxis / np.linalg.norm(yaxis)
                khat = np.cross(xyplane, j)  # xyplane is "defining" xaxis
                k = khat / norm(khat)
                i = np.cross(j, k)

        elif yzplane is not None:
            if yaxis is not None:
                j = yaxis / np.linalg.norm(yaxis)
                ihat = np.cross(j, yzplane)  # yzplane is "defining" zaxis
                i = ihat / np.linalg.norm(ihat)
                k = np.cross(i, j)
            elif zaxis is not None:
                k = zaxis / np.linalg.norm(zaxis)
                ihat = np.cross(yzplane, zaxis)  # yzplane is "defining" yaxis
                i = ihat / np.linalg.norm(ihat)
                j = np.cross(k, i)

        elif xzplane is not None:
            if xaxis is not None:
                i = xaxis / np.linalg.norm(xaxis)
                jhat = np.cross(xzplane, i)  # xzplane is "defining" zaxis
                j = jhat / np.linalg.norm(jhat)
                k = np.cross(i, j)
            elif zaxis is not None:
                # standard
                k = zaxis / np.linalg.norm(zaxis)
                jhat = np.cross(k, xzplane) # xzplane is "defining" xaxis
                j = jhat / np.linalg.norm(jhat)
                i = np.cross(j, k)
        self.add_ijk(cid, coord_type, origin, i, j, k, rid=rid, comment=comment)

    def add_ijk(self, cid: int, coord_type: str, origin=None, i=None, j=None, k=None,
                rid: int=0, comment: str=''):
        """
        Create a coordinate system based on 2 or 3 perpendicular unit vectors

        Parameters
        ----------
        cid : int
            the new coordinate system id
        origin : (3,) float ndarray
            defines the location of the origin in the global coordinate frame
        rid : int; default=0
            the new reference coordinate system id
        i : (3,) float ndarray
            defines the i unit vector
        j : (3,) float ndarray
            defines the j unit vector
        k : (3,) float ndarray
            defines the k unit vector

        """
        assert coord_type in ['R', 'C', 'S'], coord_type
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
                raise RuntimeError(f'j or k are None; j={j} k={k}')

        if np.abs(k).max() == 0.0 or np.abs(i).max() == 0.0:
            msg = (
                'coordinate vectors arent perpendicular\n'
                '  origin = %s\n'
                '  i = %s\n'
                '  j = %s\n'
                '  k = %s\n' % (origin, i, j, k))
            raise RuntimeError(msg)
        # origin
        e1 = origin
        # point on z axis
        e2 = origin + k

        # point on x-z plane / point on x axis
        e3 = origin + i
        cardi = (2, coord_type, cid), (rid, e1, e2, e3), comment
        self.cards2.append(cardi)
        self.n += 1
        #return cls(cid, e1, e2, e3, rid=rid, comment=comment)

    def sort(self) -> None:
        if np.array_equal(self.coord_id, np.unique(self.coord_id)):
            return
        isort = np.argsort(self.coord_id)
        self.__apply_slice__(self, isort)

    def slice_card_by_coord_id(self, coord_id: np.ndarray) -> COORD:
        assert len(self.coord_id) > 0, self.coord_id
        i = self.index(coord_id)
        coord = COORD(self.model)
        self.__apply_slice__(coord, i)
        return coord

    def add_cord1r(self, cid: int, g1: int, g2: int, g3: int, comment: str='') -> int:
        cardi = (1, 'R', cid), (g1, g2, g3), comment
        self.cards1.append(cardi)
        self.n += 1
        return self.n

    def add_cord1c(self, cid: int, g1: int, g2: int, g3: int, comment: str='') -> int:
        cardi = (1, 'C', cid), (g1, g2, g3), comment
        self.cards1.append(cardi)
        self.n += 1
        return self.n

    def add_cord1s(self, cid: int, g1: int, g2: int, g3: int, comment: str='') -> int:
        cardi = (1, 'S', cid), (g1, g2, g3), comment
        self.cards1.append(cardi)
        self.n += 1
        return self.n

    def add_cord1r_bdf(self, card: BDFCard, icard: int, comment: str='') -> int:
        ncoord = icard * 5
        cid = integer(card, 1 + ncoord, 'cid')
        grid_origin = integer(card, 2 + ncoord, 'g1')
        grid_zaxis = integer(card, 3 + ncoord, 'g2')
        grid_xzplane = integer(card, 4 + ncoord, 'g3')
        self.nodes[icard, :] = [grid_origin, grid_zaxis, grid_xzplane]
        cardi = (1, 'R', cid), (grid_origin, grid_zaxis, grid_xzplane), comment
        self.cards1.append(cardi)
        self.n += 1
        return self.n

    def add_cord1c_bdf(self, card: BDFCard, icard: int, comment: str='') -> int:
        ncoord = icard * 5
        cid = integer(card, 1 + ncoord, 'cid')
        grid_origin = integer(card, 2 + ncoord, 'g1')
        grid_zaxis = integer(card, 3 + ncoord, 'g2')
        grid_xzplane = integer(card, 4 + ncoord, 'g3')
        self.nodes[icard, :] = [grid_origin, grid_zaxis, grid_xzplane]
        cardi = (1, 'C', cid), (grid_origin, grid_zaxis, grid_xzplane), comment
        self.cards1.append(cardi)
        self.n += 1
        return self.n

    def add_cord1s_bdf(self, card: BDFCard, icard: int, comment: str='') -> int:
        ncoord = icard * 5
        cid = integer(card, 1 + ncoord, 'cid')
        grid_origin = integer(card, 2 + ncoord, 'g1')
        grid_zaxis = integer(card, 3 + ncoord, 'g2')
        grid_xzplane = integer(card, 4 + ncoord, 'g3')
        self.nodes[icard, :] = [grid_origin, grid_zaxis, grid_xzplane]
        cardi = (1, 'S', cid), (grid_origin, grid_zaxis, grid_xzplane), comment
        self.cards1.append(cardi)
        self.n += 1
        return self.n

    def add_cord2r(self, cid: int,
                   origin: np.ndarray | list[float],
                   zaxis: np.ndarray | list[float],
                   xzplane: np.ndarray | list[float],
                   rid: int=0, setup: bool=True, comment: str='') -> int:
        """
        Creates the CORD2R card, which defines a rectangular coordinate
        system using 3 vectors.

        Parameters
        ----------
        cid : int
            coordinate system id
        rid : int; default=0
            the referenced coordinate system that defines the system the
            vectors
        origin : list[float, float, float]
            the origin of the coordinate system
        zaxis : list[float, float, float]
            the z-axis of the coordinate system
        xzplane : list[float, float, float]
            a point on the xz plane
        comment : str; default=''
            a comment for the card

        """
        origin, zaxis, xzplane = _default_cord2_axes(origin, zaxis, xzplane)
        cardi = (2, 'R', cid), (rid, origin, zaxis, xzplane), comment
        self.cards2.append(cardi)
        self.n += 1
        return self.n

    def add_cord2c(self, cid: int,
                   origin: np.ndarray | list[float],
                   zaxis: np.ndarray | list[float],
                   xzplane: np.ndarray | list[float],
                   rid: int=0, setup: bool=True, comment: str='') -> int:
        origin, zaxis, xzplane = _default_cord2_axes(origin, zaxis, xzplane)
        cardi = (2, 'C', cid), (rid, origin, zaxis, xzplane), comment
        self.cards2.append(cardi)
        self.n += 1
        return self.n

    def add_cord2s(self, cid: int,
                   origin: np.ndarray | list[float],
                   zaxis: np.ndarray | list[float],
                   xzplane: np.ndarray | list[float],
                   rid: int=0, setup: bool=True, comment: str='') -> int:
        origin, zaxis, xzplane = _default_cord2_axes(origin, zaxis, xzplane)
        cardi = (2, 'S', cid), (rid, origin, zaxis, xzplane), comment
        self.cards2.append(cardi)
        self.n += 1
        return self.n

    #def _add_cord2(self, coordtype: str, card: BDFCard, comment: str='') -> int:
        #assert isinstance(coordtype, str), coordtype
        #assert isinstance(card, BDFCard), card
        #self.cards.append((2, coordtype, 0, card, comment))
        #i = len(self.cards) - 1
        #self.n += 1
        #return i

    def add_cord2r_bdf(self, card: BDFCard, comment: str='') -> int:
        cid, rid, origin, zaxis, xzplane = _parse_cord2x(card, 'CORD2R')
        cardi = (2, 'R', cid), (rid, origin, zaxis, xzplane), comment
        self.cards2.append(cardi)
        self.n += 1
        return self.n

    def add_cord2c_bdf(self, card: BDFCard, comment: str='') -> int:
        cid, rid, origin, zaxis, xzplane = _parse_cord2x(card, 'CORD2C')
        cardi = (2, 'C', cid), (rid, origin, zaxis, xzplane), comment
        self.cards2.append(cardi)
        self.n += 1
        return self.n

    def add_cord2s_bdf(self, card: BDFCard, comment: str='') -> int:
        cid, rid, origin, zaxis, xzplane = _parse_cord2x(card, 'CORD2S')
        cardi = (2, 'S', cid), (rid, origin, zaxis, xzplane), comment
        self.cards2.append(cardi)
        self.n += 1
        return self.n

    #def add_cord2r(self, card: BDFCard, comment: str='') -> int:
        #return self._add_cord2('R', card, comment)
    #def add_cord2c(self, card: BDFCard, comment: str='') -> int:
        #return self._add_cord2('C', card, comment)
    #def add_cord2s(self, card: BDFCard, comment: str='') -> int:
        #return self._add_cord2('S', card, comment)

    def parse_cards(self) -> None:
        if self.n == 0:
            return
        ncoords1 = len(self.cards1)
        ncoords2 = len(self.cards2)
        ncoords = ncoords1 + ncoords2
        if ncoords == 0:
            return
        assert ncoords > 0, ncoords

        dn = self.n - ncoords
        i0 = len(self.icoord)
        if dn:
            self.coord_type = np.hstack([self.coord_type, np.zeros(ncoords, dtype='|U1')])
            self.icoord = np.hstack([self.icoord, np.zeros(ncoords, dtype='int8')])

            #self.coord_id1 = np.hstack([self.coord_id, np.zeros(ncoords1, dtype='int32')])
            self.nodes = np.vstack([self.nodes, np.full((ncoords, 3), -1, dtype='int32')])

            self.coord_id = np.hstack([self.coord_id, np.full(ncoords, -1, dtype='int32')])
            self.ref_coord_id = np.hstack([self.ref_coord_id, np.full(ncoords, -1, dtype='int32')])

            self.e1 = np.vstack([self.e1, np.full((ncoords, 3), np.nan, dtype='float64')])
            self.e2 = np.vstack([self.e2, np.full((ncoords, 3), np.nan, dtype='float64')])
            self.e3 = np.vstack([self.e3, np.full((ncoords, 3), np.nan, dtype='float64')])

            self.origin   = np.vstack([self.origin,   np.full((ncoords, 3), np.nan, dtype='float64')])
            self.z_axis   = np.vstack([self.z_axis,   np.full((ncoords, 3), np.nan, dtype='float64')])
            self.xz_plane = np.vstack([self.xz_plane, np.full((ncoords, 3), np.nan, dtype='float64')])

            #TT = np.array([ones], dtype='float64')

            self.is_resolved = np.hstack([self.is_resolved, np.full(ncoords, False, dtype='bool')])
            self.i = np.vstack([self.i, np.full((ncoords, 3), np.nan, dtype='float64')])
            self.j = np.vstack([self.j, np.full((ncoords, 3), np.nan, dtype='float64')])
            self.k = np.vstack([self.k, np.full((ncoords, 3), np.nan, dtype='float64')])
            self.T = np.vstack([self.T, np.full((ncoords, 3, 3), np.nan, dtype='float64')])
            #self.is_resolved[0] = True
            ncoords_actual = len(self.coord_id)
            assert self.T.shape == (ncoords_actual, 3, 3), self.T.shape
            del ncoords_actual
        else:
            self.coord_type = np.zeros(ncoords, dtype='|U1')
            self.icoord = np.zeros(ncoords, dtype='int8')

            self.nodes = np.full((ncoords, 3), -1, dtype='int32')
            self.coord_id = np.full(ncoords, -1, dtype='int32')
            self.ref_coord_id = np.full(ncoords, -1, dtype='int32')

            self.e1 = np.full((ncoords, 3), np.nan, dtype='float64')
            self.e2 = np.full((ncoords, 3), np.nan, dtype='float64')
            self.e3 = np.full((ncoords, 3), np.nan, dtype='float64')

            self.origin = np.full((ncoords, 3), np.nan, dtype='float64')
            #self.z_axis = np.full((ncoords, 3), np.nan, dtype='float64')
            #self.xz_plane = np.full((ncoords, 3), np.nan, dtype='float64')
            self.is_resolved = np.full(ncoords, False, dtype='bool')
            self.i = np.full((ncoords, 3), np.nan, dtype='float64')
            self.j = np.full((ncoords, 3), np.nan, dtype='float64')
            self.k = np.full((ncoords, 3), np.nan, dtype='float64')
        assert self.nodes.shape[0] == len(self.is_resolved)
        assert self.nodes.shape[0] == len(self.i)
        assert self.nodes.shape[0] == len(self.e1)
        assert self.nodes.shape[0] == self.n

        icard = i0
        for card in self.cards1:
            (icoord, coord_type, cid), cardi, comment = card
            #print((icoord, coord_type, cid))
            (grid_origin, grid_zaxis, grid_xzplane) = cardi
            assert icoord == 1, icoord

            self.coord_id[icard] = cid
            self.icoord[icard] = icoord
            self.coord_type[icard] = coord_type

            self.nodes[icard, :] = [grid_origin, grid_zaxis, grid_xzplane]
            self.xyz_to_global_transform[cid] = np.full((3, 3), np.nan, dtype='float64')
            icard += 1

        icard = i0 + len(self.cards1)
        for card in self.cards2:
            (icoord, coord_type, cid), cardi, comment = card
            #print((icoord, coord_type, cid))

            assert icoord == 2, icoord
            rid, origin, zaxis, xzplane = cardi

            self.coord_id[icard] = cid
            self.icoord[icard] = icoord
            self.coord_type[icard] = coord_type

            #self.icoord[icard] = icoord
            self.ref_coord_id[icard] = rid
            #self.origin[icard, :] = origin
            #self.z_axis[icard, :] = zaxis
            #self.xz_plane[icard, :] = xzplane
            self.e1[icard, :] = origin
            self.e2[icard, :] = zaxis
            self.e3[icard, :] = xzplane

            self.xyz_to_global_transform[cid] = np.full((3, 3), np.nan, dtype='float64')
            icard += 1

        assert self.coord_id.min() >= 0, self.coord_id
        assert self.icoord.min() >= 1 and self.icoord.max() <= 3, (self.icoord.min(), self.icoord.max())
        #assert len(np.unique(self.coord_id)) == len(self.coord_id), self.coord_id


        #assert self.nodes.shape[0] == len(self.is_resolved)
        #assert self.nodes.shape[0] == len(self.i)
        #assert self.nodes.shape[0] == len(self.e1)
        #assert self.nodes.shape[0] == self.n
        self.sort()
        ucoords = np.unique(self.coord_id)
        if len(self.coord_id) != len(ucoords):
            #print('slicing for unique coords')
            i = self.index(ucoords)
            #self.slice_card_by_coord_id(ucoords)
            self.__apply_slice__(self, i)

        np.searchsorted(self.coord_id, self.coord_id)

        #print('self.coord_id', self.coord_id)
        assert self.nodes.shape[0] == len(self.is_resolved)
        assert self.nodes.shape[0] == len(self.i)
        assert self.nodes.shape[0] == len(self.e1)
        assert self.nodes.shape[0] == self.n

        self.cards1 = []
        self.cards2 = []
        self.setup()

    #def sort(self):
        #isort = np.argsort(self.coord_id)
        #self.coord_id = self.coord_id[isort]
        #self.ref_coord_id = self.ref_coord_id[isort]
        #self.is_resolved = self.is_resolved[isort]
        #self.origin = self.origin[isort, :]
        #self.z_axis = self.z_axis[isort, :]
        #self.xz_plane = self.xz_plane[isort, :]

    def setup(self):
        model = self.model
        log = model.log
        assert len(np.unique(self.coord_id)) == len(self.coord_id)
        self.xyz_to_global_transform = {
            0: np.eye(3, dtype='float64')
        }

        grid = model.grid
        nresolved = self.is_resolved.sum()
        ncoords = len(self.is_resolved)

        resolved = set()
        if self.coord_id[0] == 0:
            # cid=0 definition
            self.i[0, :] = [1., 0., 0.]
            self.j[0, :] = [0., 1., 0.]
            self.k[0, :] = [0., 0., 1.]
            resolved.add(0)
        self.resolve0(resolved, nresolved, ncoords)

        ## unresolved coords
        unresolved_cids = set(self.coord_id.tolist())
        if 0 in unresolved_cids:
            unresolved_cids.remove(0)

        rid_to_i_icoord_coordtype = {cid: (i, icoord, coord_type) for i, cid, icoord, coord_type in
                                     zip(count(), self.coord_id, self.icoord, self.coord_type)}

        grid = self.model.grid
        if grid.n:
            grid.sort()
        assert len(np.unique(self.coord_id)) == len(self.coord_id)

        while nresolved < ncoords:
            i1 = np.where(self.icoord == 1)[0]
            i2 = np.where(self.icoord == 2)[0]
            resolvable_rids = {rid for rid in self.ref_coord_id[i2] if rid in resolved}
            coord2_cids_to_resolve = [cid for cid, rid in zip(self.coord_id[i2], self.ref_coord_id[i2])
                                      if rid in resolvable_rids and cid not in resolved]

            coord1_cids_to_resolve = self._find_cord1s_to_resolve(
                grid, i1, resolved)

            assert len(np.unique(self.coord_id)) == len(self.coord_id)
            coords_to_resolve = coord1_cids_to_resolve + coord2_cids_to_resolve
            if len(coords_to_resolve) == 0:
                raise RuntimeError(f'cannot resolve any coordinate systems...unresolved_cids={unresolved_cids}')

            log.debug(f'to_resolve: coord1={np.array(coord1_cids_to_resolve)}; coord2={np.array(coord2_cids_to_resolve)}')
            #nresolved0 = nresolved
            nresolved = self._resolve_cord1(coord1_cids_to_resolve, nresolved, resolved, grid,
                                            unresolved_cids)
            #if nresolved > nresolved0:
                #log.debug(f'resolved CORD1x: {coord1_cids_to_resolve}')

            #nresolved0 = nresolved
            nresolved = self._resolve_cord2(
                coord2_cids_to_resolve, resolved, nresolved, unresolved_cids,
                rid_to_i_icoord_coordtype)

            #if nresolved > nresolved0:
                #log.debug(f'resolved CORD2x: {coord2_cids_to_resolve}')
            if unresolved_cids:
                log.debug(f'unresolved_cids = {unresolved_cids}')

        if 0 in resolved:
            # just limiting log messages
            resolved.remove(0)
        if resolved:
            log.info(f'resolved = {np.array(list(resolved))}')
        if unresolved_cids:
            raise RuntimeError(f'unresolved_cids = {unresolved_cids}')
        assert np.all(np.isfinite(self.origin)), self.origin
        assert np.all(np.isfinite(self.e1)), self.e1
        assert np.all(np.isfinite(self.T)), self.T
        #self.model.log.warning(f'finished setting up coords; resolved={resolved}')
        return

    def resolve0(self, resolved: list[int], nresolved: int, ncoords: int) -> int:
        if nresolved != ncoords and 0:
            ## find the coords that reference cid=0
            ## make sure to drop cid=0
            izero = np.where(self.ref_coord_id[1:] == 0)[0]

            ## resolve the coordinate systems that reference cid=0
            nrid_zero = len(izero)
            if nrid_zero:
                # correct for dropping cid=0
                izero += 1

                k = self.z_axis[izero, :] - self.origin[izero, :]
                knorm = np.linalg.norm(k, axis=1)
                #assert len(knorm) == nrid_zero
                k /= knorm[:, np.newaxis]

                ik = self.xz_plane[izero, :] - self.origin[izero, :]
                j = np.cross(k, ik, axis=1)
                jnorm = np.linalg.norm(j, axis=1)
                j /= jnorm[:, np.newaxis]
                #assert len(jnorm) == nrid_zero

                #assert j.shape[0] == nrid_zero
                i = np.cross(k, j, axis=1)
                self.i[izero, :] = i
                self.j[izero, :] = j
                self.k[izero, :] = k
                self.is_resolved[izero] = True
                nresolved += nrid_zero
                resolved += self.coord_id[izero].tolist()
        return nresolved

    def _find_cord1s_to_resolve(self, grid: GRID, i1: np.ndarray, resolved: Set[int]) -> list[int]:
        coord1_cids_to_resolve = []
        for cid, nodes in zip(self.coord_id[i1], self.nodes[i1, :]):
            inid = grid.index(nodes)
            cps_for_nodes = grid.cp[inid]
            is_resolved = [cpi in resolved for cpi in cps_for_nodes]
            if all(is_resolved):
                coord1_cids_to_resolve.append(cid)
        return coord1_cids_to_resolve

    def _resolve_cord1(self, coord1_cids_to_resolve: list[int],
                       nresolved: int, resolved: Set[int], grid: GRID,
                       unresolved_cids: Set[int]) -> int:
        """resolve CORD1R, CORD1S, CORD1C"""
        if not coord1_cids_to_resolve:
            return nresolved
        coord1_cids_to_resolve.sort()
        inids = np.searchsorted(self.coord_id, coord1_cids_to_resolve)
        for cid, i in zip(coord1_cids_to_resolve, inids):
            cid = self.coord_id[i]
            #rid = self.ref_coord_id[i]
            coord_type = self.coord_type[i]
            icoord = self.icoord[i]
            node_id = self.nodes[i, :]
            assert coord_type in {'R', 'C', 'S'}, coord_type
            assert icoord == 1, icoord

            assert node_id.min() > 0, node_id
            gridi = grid.slice_card_by_node_id(node_id, sort_ids=True)
            grid_xyz_cid0 = gridi.xyz_cid0()
            grid_node_id = gridi.node_id
            inid = np.searchsorted(grid_node_id, node_id)
            #assert np.array_equal(grid_node_id, node_id)
            grid_xyz_cid0 = grid_xyz_cid0[inid, :]
            del inid, grid_node_id

            #these all should be in the global frame
            origin = grid_xyz_cid0[0, :]
            zaxis = grid_xyz_cid0[1, :]
            xzplane = grid_xyz_cid0[2, :]

            self.e1[i, :] = origin
            self.e2[i, :] = zaxis
            self.e3[i, :] = xzplane

            self.origin[i, :] = origin
            self.z_axis[i, :] = zaxis
            self.xz_plane[i, :] = xzplane

            msg = f'CORD1{coord_type} cid={cid}; nodes={node_id}\n'
            iaxis, jaxis, kaxis = axes_to_coord_vectors(origin, zaxis, xzplane, msg=msg)
            self.i[i, :] = iaxis
            self.j[i, :] = jaxis
            self.k[i, :] = kaxis
            beta = np.vstack([iaxis, jaxis, kaxis])
            #self.model.log.info(f'saving T{i} for cid={cid}\n  beta={beta}')
            self.T[i, :, :] = beta
            self.xyz_to_global_transform[cid] = beta
            self.is_resolved[i] = True
            resolved.add(cid)
            unresolved_cids.remove(cid)
            nresolved += 1
        return nresolved

    def _resolve_cord2(self, coord2_cids_to_resolve: list[int],
                       resolved: Set[int], nresolved: int,
                       unresolved_cids: Set[int],
                       rid_to_i_icoord_coordtype: dict[int, Any]) -> int:
        """resolve CORD2R, CORD2S, CORD2C"""
        #log = self.model.log
        icid = np.searchsorted(self.coord_id, coord2_cids_to_resolve)
        for i in icid:
            cid = self.coord_id[i]
            rid = self.ref_coord_id[i]
            coord_type = self.coord_type[i]
            icoord = self.icoord[i]
            assert coord_type in {'R', 'C', 'S'}, coord_type
            assert icoord == 2, icoord
            irid, ricoord, rcoord_type = rid_to_i_icoord_coordtype[rid]

            assert rcoord_type in {'R', 'C', 'S'}, rcoord_type
            assert ricoord == 2, ricoord

            # these will be used by the transformation
            #ri = self.i[irid, :]
            #rj = self.j[irid, :]
            #rk = self.k[irid, :]
            #beta = np.vstack([ri, rj, rk])
            #transform_local_xyz_to_global(beta, rid)
            #betar = self.xyz_to_global_transform[rid]
            #originr = self.origin[irid]

            origin0 = self.e1[i, :]
            zaxis0 = self.e2[i, :]
            xzplane0 = self.e3[i, :]
            msg = (
                f'\norigin0: {origin0}\n'
                f'zaxis0: {zaxis0}\n'
                f'xzplane0: {xzplane0}'
            )
            if 1:
                origin = self._transform_local_xyz_to_global(origin0, rid, irid)[0, :]
                zaxis = self._transform_local_xyz_to_global(zaxis0, rid, irid)[0, :]
                xzplane = self._transform_local_xyz_to_global(xzplane0, rid, irid)[0, :]
            else:
                #if coord_type == 'R':
                    #xyz_cidr = xyz_cid
                #elif coord_type == 'C':
                    #xyz_cidr = cylindrical_to_rectangular(xyz_cid)
                #elif coord_type == 'S':
                    #xyz_cidr = spherical_to_rectangular(xyz_cid)
                #else:
                    #raise NotImplementedError(coord_type)

                #xyz_cid0 = xyz_cid.copy()
                #for i, xyzi in enumerate(xyz_cidr):
                    #xyz_cid0[i, :] = xyzi @ betar + origin

                #msg = f'CORD2{coord_type} cid={cid}; rid={rid}'
                #origin = transform_local_xyz_to_global(self.e1)
                #e1i = origin
                #e2i = self.transform_local_xyz_to_global(self.e2)
                #e3i = self.transform_local_xyz_to_global(self.e3)

                if rcoord_type == 'R':
                    origin = origin0
                    zaxis = zaxis0
                    xzplane = xzplane0
                elif rcoord_type == 'C':
                    origin = transform_cylindrical_to_rectangular(origin0)
                    zaxis = transform_cylindrical_to_rectangular(zaxis0)
                    xzplane = transform_cylindrical_to_rectangular(xzplane0)
                    msg += (
                        f'\norigin0: {origin0}\n'
                        f'zaxis0: {zaxis0}\n'
                        f'xzplane0: {xzplane0}'
                    )
                elif rcoord_type == 'S':
                    origin = transform_spherical_to_rectangular(origin0)
                    zaxis = transform_spherical_to_rectangular(zaxis0)
                    xzplane = transform_spherical_to_rectangular(xzplane0)
                    msg += (
                        f'\norigin0: {origin0}\n'
                        f'zaxis0: {zaxis0}\n'
                        f'xzplane0: {xzplane0}'
                    )
                else:
                    raise RuntimeError(rcoord_type)
            #msg = (
                #f'cid={cid} rid={rid} rcoord_type={rcoord_type}\n'
                #f'  origin0: {origin0}\n'
                #f'  zaxis0: {zaxis0}\n'
                #f'  xzplane0: {xzplane0}\n\n'

                #f'  origin: {origin}\n'
                #f'  zaxis: {zaxis}\n'
                #f'  xzplane: {xzplane}'
            #)
            #self.model.log.info(msg)

            self.origin[i, :] = origin
            self.z_axis[i, :] = zaxis
            self.xz_plane[i, :] = xzplane
            iaxis, jaxis, kaxis = axes_to_coord_vectors(origin, zaxis, xzplane, msg=msg)
            self.i[i, :] = iaxis
            self.j[i, :] = jaxis
            self.k[i, :] = kaxis
            beta = np.vstack([iaxis, jaxis, kaxis])
            #self.model.log.info(f'saving T{i} for cid={cid}\n  beta={beta}')
            self.T[i, :, :] = beta

            self.xyz_to_global_transform[cid] = beta
            self.is_resolved[i] = True
            nresolved += 1

            #log.debug(f'resolved cid={cid}')
            resolved.add(cid)
            unresolved_cids.remove(cid)
        return nresolved

    #@parse_load_check
    def write_file(self, bdf_file: TextIOLike,
                   size: int=8, is_double: bool=False,
                   write_card_header: bool=False) -> None:
        if len(self.coord_id) == 0:
            return ''

        max_int = max(self.coord_id.max(), self.nodes.max(), self.ref_coord_id.max())
        print_card, size = get_print_card_size(size, max_int)

        class_name = self.type
        nodes = array_str(self.nodes, size=size)
        ref_coord_ids = array_default_int(self.ref_coord_id, default=0, size=size)
        #array_str
        e1s = array_float(self.e1, size=size, is_double=is_double)
        e2s = array_float(self.e2, size=size, is_double=is_double)
        e3s = array_float(self.e3, size=size, is_double=is_double)
        if size == 8:
            for icoord, coord_type, cid, rid, nodesi, origin, zaxis, xzplane in zip(
                    self.icoord, self.coord_type,
                    self.coord_id, ref_coord_ids, nodes,
                    e1s, e2s, e3s):
                if cid == 0:
                    continue
                if icoord == 1:
                    class_name = 'CORD1%s' % coord_type
                    msg = '%-8s%8d%8s%8s%8s\n' % (
                        class_name, cid, nodesi[0], nodesi[1], nodesi[2], )
                else:
                    class_name = 'CORD2%s' % coord_type
                    msg = '%-8s%8d%8s%8s%8s%8s%8s%8s%s\n        %8s%8s%8s\n' % (
                        class_name, cid, rid,
                        origin[0], origin[1], origin[2],
                        zaxis[0], zaxis[1], zaxis[2],
                        xzplane[0], xzplane[1], xzplane[2],
                    )
                    assert rid != -1
                bdf_file.write(msg)
        else:
            for icoord, coord_type, cid, rid, nodesi, origin, zaxis, xzplane in zip(
                    self.icoord, self.coord_type,
                    self.coord_id, ref_coord_ids, nodes,
                    e1s, e2s, e3s):
                if cid == 0:
                    continue
                if icoord == 1:
                    class_name = 'CORD1%s' % coord_type
                    fields = [class_name, cid, nodesi[0], nodesi[1], nodesi[2]]
                else:
                    class_name = 'CORD2%s' % coord_type
                    fields = [class_name, cid, rid,
                              origin[0], origin[1], origin[2],
                              zaxis[0], zaxis[1], zaxis[2],
                              xzplane[0], xzplane[1], xzplane[2],
                              ]
                    assert rid != -1
                bdf_file.write(print_card(fields))
        return

    def transform_force_local_to_global(self, force: np.ndarray, local_coord_id: int=0) -> np.ndarray:
        force2 = np.full(force.shape, np.nan, dtype=force.dtype)
        T = self.xyz_to_global_transform[local_coord_id]

        icid = self.index(local_coord_id).squeeze()
        coord_type = self.coord_type[icid]
        assert coord_type == 'R', coord_type
        #self.model.log.warning(f'force.shape={force.shape}')
        #self.model.log.warning(f'T.shape={T.shape}')
        for i, forcei in enumerate(force):
            force2[i, :] = forcei @ T # reversed?
        return force2

    def check_missing_ids(self, coord_id: np.ndarray):
        missing_coords = np.setdiff1d(coord_id, self.coord_id)
        if len(missing_coords):
            raise RuntimeError(f'coords={missing_coords} not found in {self.coord_id}')

    def transform_local_xyz_to_global_coords(self, xyz: np.ndarray, coord_id: np.ndarray) -> np.ndarray:
        """takes a consistent set of xyz and cd values and transforms them"""
        self.check_missing_ids(coord_id)
        xyz2 = np.zeros(xyz.shape, dtype=xyz.dtype)
        for i, xyzi, cid in zip(count(), xyz, coord_id):
            xyz2[i, :] = self.transform_local_xyz_to_global(xyzi, cid)
        return xyz2

    def transform_xyz_to_global_assuming_rectangular(self, xyz: np.ndarray) -> np.ndarray:
        assert len(xyz) == len(self.coord_id)
        xform = self.T
        xyz2 = np.einsum('ni,nij->nj', xyz, xform)
        assert xyz.shape == xyz2.shape, (xyz.shape, xyz2.shape)
        return xyz2

    def transform_node_to_local_xyz(self, node_id: np.ndarray, cid: int) -> np.ndarray:
        xyz_cid0 = self.transform_node_to_global_xyz(node_id)
        #beta = self.xyz_to_global_transform[cid]
        xyz_cid = self.model.coord.transform_global_xyz_to_local_coord_id(xyz_cid0, cid)
        return xyz_cid

    def transform_node_to_global_xyz(self, node_id: np.ndarray) -> np.ndarray:
        grid = self.model.grid
        grid2 = grid.slice_card_by_node_id(node_id, sort_ids=False)
        xyz_cid0 = grid2.xyz_cid0()
        return xyz_cid0

    def transform_offset_xyz_to_global_xyz(self, xyz: np.ndarray, coord_id: np.ndarray) -> np.ndarray:
        self.check_missing_ids(coord_id)

        icoord = np.searchsorted(self.coord_id, coord_id)
        origin = self.origin[icoord, :]

        T = self.T[icoord, :, :]
        xyz2 = np.einsum('ni,nij->nj', xyz, T)
        assert xyz2.shape == xyz.shape

        xyz2 += origin
        assert xyz2.shape == xyz.shape
        #xyz2 = xyz @ T + origin
        return xyz2

    def transform_node_to_local_coord_id(self, node_id: np.ndarray, local_coord_id: int) -> np.ndarray:
        # transform_node_to_local_xyz
        # transform_node_to_local_coord_id
        ## TODO: this is redundant
        grid2 = self.model.grid.slice_card_by_node_id(node_id, sort_ids=False)
        xyz_cid0 = grid2.xyz_cid0()
        xyz_cid = self.transform_global_xyz_to_local_coord_id(xyz_cid0, local_coord_id)
        return xyz_cid

    def transform_local_xyz_to_global(self,
                                      xyz_cid: np.ndarray,
                                      local_coord_id: int) -> np.ndarray:
        """transforms a series of xyz values from cid=local to cid=0"""
        xyz_cid = np.atleast_2d(xyz_cid)
        assert isinstance(local_coord_id, integer_types), local_coord_id
        if local_coord_id == 0:
            return xyz_cid.copy()

        if local_coord_id not in self.coord_id:
            raise KeyError(f'local_coord_id={local_coord_id} cannot be found')
        icid = self.index(local_coord_id).squeeze()
        if not self.coord_id[icid] == local_coord_id:
            raise KeyError(f'icid={icid} local_coord_id={local_coord_id} actual={self.coord_id[icid]}\n'
                           f'coords={self.coord_id}')
        xyz_cid0 = self._transform_local_xyz_to_global(xyz_cid, local_coord_id, icid)
        return xyz_cid0

    def _transform_local_xyz_to_global(self,
                                       xyz_cid: np.ndarray,
                                       rid: int,
                                       ilocal_coord_id: int) -> np.ndarray:
        """transforms a series of xyz values in local_coord_id"""
        xyz_cid = np.atleast_2d(xyz_cid)
        origin = self.origin[ilocal_coord_id, :]
        beta = self.T[ilocal_coord_id, :, :]
        assert np.all(np.isfinite(origin)), origin
        assert np.all(np.isfinite(beta)), beta
        coord_type = self.coord_type[ilocal_coord_id]
        if coord_type == 'R':
            xyz_cidr = xyz_cid
        elif coord_type == 'C':
            xyz_cidr = cylindrical_to_rectangular(xyz_cid)
        elif coord_type == 'S':
            xyz_cidr = spherical_to_rectangular(xyz_cid)
        else:
            raise NotImplementedError(coord_type)

        xyz_cid0 = xyz_cid.copy()

        ## TODO: vectorize this
        for i, xyzi in enumerate(xyz_cidr):
            #print(xyzi.shape, beta.shape)
            xyzi @ beta
            xyz_cid0[i, :] = xyzi @ beta
            #xyz_cid0 = xyz_cidr @ beta + origin

        xyz_cid0 += origin
        #nnodes = xyz_cidr.shape[0]

        xyz_cid0_new = xyz_cidr @ beta + origin
        for xyz, xyz_new in zip(xyz_cid0, xyz_cid0_new):
            assert np.allclose(xyz, xyz_new)
        #if nnodes > 1:
            ## why is this shape wrong?
            #xyz_cid0_new = np.einsum('ij,ijk->ik', xyz_cidr, beta[np.newaxis, :, :])
            #asdf
        return xyz_cid0

    def transform_global_xyz_to_local_coord_id(self, xyz_cid0: np.ndarray, local_coord_id: int) -> np.ndarray:
        xyz_cid0 = np.atleast_2d(xyz_cid0)

        xyz_cid = np.full(xyz_cid0.shape, np.nan, dtype=xyz_cid0.dtype)
        Ti = self.xyz_to_global_transform[local_coord_id]
        assert Ti.shape == (3, 3), Ti.shape
        T = Ti.T

        icid = self.index(local_coord_id)
        assert self.coord_id[icid] == local_coord_id
        origin = self.origin[icid, :]

        xyz_local = xyz_cid0 - origin # [np.newaxis, :]
        assert xyz_local.shape == xyz_cid.shape
        for i, xyzr in enumerate(xyz_local):
            xyzi_cid = xyzr @ T # reversed?
            xyz_cid[i, :] = xyzi_cid

        coord_type = self.coord_type[icid]
        if coord_type == 'R':
            pass
            #xyz_cidr = xyz[icp, :]
        #elif coord_type == 'C':
            #xyz_cidr = cylindrical_to_rectangular(xyz[icp, :])
        #elif coord_type == 'S':
            #xyz_cidr = spherical_to_rectangular(xyz[icp, :])
        else:
            raise NotImplementedError(coord_type)
        #xyz_cid0i = xyz_cidr @ beta + origin
        #xyz_cid0i - origin = xyz_cidr @ beta


        assert coord_type == 'R', coord_type
        return xyz_cid

def axes_to_coord_vectors(origin, zaxis, xzplane, msg: str) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    kaxis = zaxis - origin
    ik = xzplane - origin

    knorm = np.linalg.norm(kaxis)
    if knorm == 0.:
        msg2 = (
            f'Cannot normalize k axis; {msg}\n'
            f'origin={origin}\n'
            f'zaxis={zaxis}\n'
            f'kaxis={kaxis}\n'
        )
        raise RuntimeError(msg2)
    kaxis /= knorm

     #j = k x i
    jaxis = np.cross(kaxis, ik)
    jnorm = np.linalg.norm(jaxis)
    if jnorm == 0.:
        msg2 = (
            f'Cannot normalize j axis; {msg}\n'
            f'origin={origin}\n'
            f'zaxis={zaxis}\n'
            f'kaxis={kaxis}\n'
            f'jaxis={jaxis}\n'
        )
        raise RuntimeError(msg2)
    assert jnorm > 0., jnorm
    jaxis /= jnorm

    iaxis = np.cross(jaxis, kaxis)
    return iaxis, jaxis, kaxis

def transform_cylindrical_to_rectangular(rtz: np.ndarray) -> np.ndarray:
    R, t, z = rtz
    theta = np.radians(t)
    x = R * np.cos(theta)
    y = R * np.sin(theta)
    xyz = np.array([x, y, z], dtype=rtz.dtype)
    return xyz

def transform_spherical_to_rectangular(rtp: np.ndarray) -> np.ndarray:
    radius = rtp[0]
    theta = np.radians(rtp[1])
    phi = np.radians(rtp[2])
    x = radius * np.sin(theta) * np.cos(phi)
    y = radius * np.sin(theta) * np.sin(phi)
    z = radius * np.cos(theta)
    xyz = np.array([x, y, z], dtype=rtp.dtype)
    return xyz

#def _parse_cord1x(card: BDFCard) -> tuple[int, int, int, int]:
    #cid = integer(card, 1, 'cid')
    #return cid, grid_origin, grid_zaxis, grid_xzplane

def _parse_cord2x(card: BDFCard,
                  card_type: str) -> tuple[int, int,
                                           np.ndarray, np.ndarray, np.ndarray]:
    cid = integer(card, 1, 'cid')

    #: reference coordinate system ID
    rid = integer_or_blank(card, 2, 'rid', default=0)

    #: origin in a point relative to the rid coordinate system
    origin = np.array([double_or_blank(card, 3, 'e1x', default=0.0),
                       double_or_blank(card, 4, 'e1y', default=0.0),
                       double_or_blank(card, 5, 'e1z', default=0.0)],
                      dtype='float64')
    #: z-axis in a point relative to the rid coordinate system
    zaxis = np.array([double_or_blank(card, 6, 'e2x', default=0.0),
                      double_or_blank(card, 7, 'e2y', default=0.0),
                      double_or_blank(card, 8, 'e2z', default=0.0)],
                     dtype='float64')
    #: a point on the xz-plane relative to the rid coordinate system
    xzplane = np.array([double_or_blank(card, 9, 'e3x', default=0.0),
                        double_or_blank(card, 10, 'e3y', default=0.0),
                        double_or_blank(card, 11, 'e3z', default=0.0)],
                       dtype='float64')
    assert len(card) <= 12, f'len({card_type} card) = {len(card):d}\ncard={card}'
    return cid, rid, origin, zaxis, xzplane

def _default_cord2_axes(origin: np.ndarray | list[float],
                        zaxis: np.ndarray | list[float],
                        xzplane: np.ndarray | list[float]) -> tuple[np.ndarray | list[float],
                                                                    np.ndarray | list[float],
                                                                    np.ndarray | list[float]]:
    if origin is None:
        origin = np.array([0., 0., 0.], dtype='float64')
    if zaxis is None:
        zaxis = np.array([0., 0., 1.], dtype='float64')
    if xzplane is None:
        xzplane = np.array([1., 0., 0.], dtype='float64')
    return origin, zaxis, xzplane

def coord_to_xyz(xyz: np.ndarray, coord_type: str) -> np.ndarray:
    if coord_type == 'R':
        pass
    elif coord_type == 'C':
        xyz = transform_cylindrical_to_rectangular(xyz)
    elif coord_type == 'S':
        xyz = transform_spherical_to_rectangular(xyz)
    else:
        raise RuntimeError(coord_type)
    return xyz
