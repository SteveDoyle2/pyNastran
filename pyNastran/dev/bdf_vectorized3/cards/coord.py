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
from pyNastran.bdf.bdf_interface.assign_type_force import force_double_or_blank
from pyNastran.femutils.coord_transforms import (
    #xyz_to_rtz_array,
    #xyz_to_rtp_array, # xyz to xxx transforms
    rtz_to_xyz_array as cylindrical_to_rectangular,
    rtp_to_xyz_array as spherical_to_rectangular, # xxx to xyz transforms
    #rtz_to_rtp_array, rtp_to_rtz_array, # rtp/rtz and rtz/rtp transforms
)

from pyNastran.bdf.bdf_interface.bdf_card import BDFCard
from pyNastran.bdf.cards.coordinate_systems import (
    #_fix_xyz_shape,
    setup_add_axes,
    setup_add_ijk,
    transform_cylindrical_to_rectangular,
    transform_spherical_to_rectangular,
)
# ------------------------------------------------------------------------------
from pyNastran.dev.bdf_vectorized3.cards.base_card import parse_check
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
    def clear(self) -> None:
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

    #@classmethod
    #def slice_card_by_coord_id(cls, coord_id: np.ndarray) -> COORD:
        #i = cls.index(coord_id)
        #coord = cls.slice_card_by_index(i)
        #return coord

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
                 ifile: int=0, comment: str=''):
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
        origin, i, j, k = setup_add_axes(
            cid, coord_type, rid=rid, origin=origin,
            xaxis=xaxis, yaxis=yaxis, zaxis=zaxis,
            xyplane=xyplane, yzplane=yzplane, xzplane=xzplane,)
        n = self.add_ijk(cid, coord_type, origin, i, j, k, rid=rid, comment=comment)
        return n

    def add_ijk(self, cid: int, coord_type: str,
                origin=None, i=None, j=None, k=None,
                rid: int=0, ifile: int=0, comment: str=''):
        """
        Create a coordinate system based on 2 or 3 perpendicular unit vectors

        Parameters
        ----------
        cid : int
            the new coordinate system id
        coord_type : str
            R, C, S
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
        origin, e1, e2, e3 = setup_add_ijk(coord_type, origin, i, j, k)
        cardi = (2, coord_type, cid), (rid, e1, e2, e3), ifile, comment
        self.cards2.append(cardi)
        self.n += 1
        return self.n - 1
        #return cls(cid, e1, e2, e3, rid=rid, comment=comment)

    def sort(self) -> None:
        if np.array_equal(self.coord_id, np.unique(self.coord_id)):
            return
        isort = np.argsort(self.coord_id)
        self.__apply_slice__(self, isort)

    #def slice_card_by_coord_id(self, coord_id: np.ndarray) -> COORD:
        #assert len(self.coord_id) > 0, self.coord_id
        #i = self.index(coord_id)
        #coord = COORD(self.model)
        #self.__apply_slice__(coord, i)
        #return coord

    def add_cord1r(self, cid: int, g1: int, g2: int, g3: int,
                   ifile: int=0, comment: str='') -> int:
        cardi = (1, 'R', cid), (g1, g2, g3), ifile, comment
        self.cards1.append(cardi)
        self.n += 1
        return self.n - 1

    def add_cord1c(self, cid: int, g1: int, g2: int, g3: int,
                   ifile: int=0, comment: str='') -> int:
        cardi = (1, 'C', cid), (g1, g2, g3), ifile, comment
        self.cards1.append(cardi)
        self.n += 1
        return self.n - 1

    def add_cord1s(self, cid: int, g1: int, g2: int, g3: int,
                   ifile: int=0, comment: str='') -> int:
        cardi = (1, 'S', cid), (g1, g2, g3), ifile, comment
        self.cards1.append(cardi)
        self.n += 1
        return self.n - 1

    def add_cord1r_bdf(self, card: BDFCard, icard: int, ifile: int, comment: str='') -> int:
        """
        Creates the CORD1R card, which defines a rectangular coordinate
        system using 3 GRIDs.
        """
        ncoord = icard * 5
        cid = integer(card, 1 + ncoord, 'cid')
        grid_origin = integer(card, 2 + ncoord, 'g1')
        grid_zaxis = integer(card, 3 + ncoord, 'g2')
        grid_xzplane = integer(card, 4 + ncoord, 'g3')
        cardi = (1, 'R', cid), (grid_origin, grid_zaxis, grid_xzplane), ifile, comment
        self.cards1.append(cardi)
        self.n += 1
        return self.n - 1

    def add_cord1c_bdf(self, card: BDFCard, icard: int,
                       ifile: int, comment: str='') -> int:
        """
        Creates the CORD1C card, which defines a rectangular coordinate
        system using 3 GRIDs.
        """
        ncoord = icard * 5
        cid = integer(card, 1 + ncoord, 'cid')
        grid_origin = integer(card, 2 + ncoord, 'g1')
        grid_zaxis = integer(card, 3 + ncoord, 'g2')
        grid_xzplane = integer(card, 4 + ncoord, 'g3')
        cardi = (1, 'C', cid), (grid_origin, grid_zaxis, grid_xzplane), ifile, comment
        self.cards1.append(cardi)
        self.n += 1
        return self.n - 1

    def add_cord1s_bdf(self, card: BDFCard, icard: int,
                       ifile: int, comment: str='') -> int:
        """
        Creates the CORD1S card, which defines a rectangular coordinate
        system using 3 GRIDs.
        """
        ncoord = icard * 5
        cid = integer(card, 1 + ncoord, 'cid')
        grid_origin = integer(card, 2 + ncoord, 'g1')
        grid_zaxis = integer(card, 3 + ncoord, 'g2')
        grid_xzplane = integer(card, 4 + ncoord, 'g3')
        cardi = (1, 'S', cid), (grid_origin, grid_zaxis, grid_xzplane), ifile, comment
        self.cards1.append(cardi)
        self.n += 1
        return self.n - 1

    def add_cord2r(self, cid: int,
                   origin: np.ndarray | list[float],
                   zaxis: np.ndarray | list[float],
                   xzplane: np.ndarray | list[float],
                   rid: int=0, setup: bool=True,
                   ifile: int=0, comment: str='') -> int:
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
        cardi = (2, 'R', cid), (rid, origin, zaxis, xzplane), ifile, comment
        self.cards2.append(cardi)
        self.n += 1
        return self.n - 1

    def add_cord2c(self, cid: int,
                   origin: np.ndarray | list[float],
                   zaxis: np.ndarray | list[float],
                   xzplane: np.ndarray | list[float],
                   rid: int=0, setup: bool=True,
                   ifile: int=0, comment: str='') -> int:
        """
        Creates the CORD2C card, which defines a rectangular coordinate
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
        cardi = (2, 'C', cid), (rid, origin, zaxis, xzplane), ifile, comment
        self.cards2.append(cardi)
        self.n += 1
        return self.n - 1

    def add_cord2s(self, cid: int,
                   origin: np.ndarray | list[float],
                   zaxis: np.ndarray | list[float],
                   xzplane: np.ndarray | list[float],
                   rid: int=0, setup: bool=True,
                   ifile: int=0, comment: str='') -> int:
        """
        Creates the CORD2S card, which defines a rectangular coordinate
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
        cardi = (2, 'S', cid), (rid, origin, zaxis, xzplane), ifile, comment
        self.cards2.append(cardi)
        self.n += 1
        return self.n - 1

    def add_cord2r_bdf(self, card: BDFCard, ifile: int, comment: str='') -> int:
        cid, rid, origin, zaxis, xzplane = _parse_cord2x(self.model, card, 'CORD2R')
        cardi = (2, 'R', cid), (rid, origin, zaxis, xzplane), ifile, comment
        self.cards2.append(cardi)
        self.n += 1
        return self.n - 1

    def add_cord2c_bdf(self, card: BDFCard, ifile: int, comment: str='') -> int:
        cid, rid, origin, zaxis, xzplane = _parse_cord2x(self.model, card, 'CORD2C')
        cardi = (2, 'C', cid), (rid, origin, zaxis, xzplane), ifile, comment
        self.cards2.append(cardi)
        self.n += 1
        return self.n - 1

    def add_cord2s_bdf(self, card: BDFCard, ifile: int, comment: str='') -> int:
        cid, rid, origin, zaxis, xzplane = _parse_cord2x(self.model, card, 'CORD2S')
        cardi = (2, 'S', cid), (rid, origin, zaxis, xzplane), ifile, comment
        self.cards2.append(cardi)
        self.n += 1
        return self.n - 1

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
        idtype = self.model.idtype
        if dn:
            self.coord_type = np.hstack([self.coord_type, np.zeros(ncoords, dtype='|U1')])
            self.icoord = np.hstack([self.icoord, np.zeros(ncoords, dtype='int8')])

            #self.coord_id1 = np.hstack([self.coord_id, np.zeros(ncoords1, dtype='int32')])
            self.nodes = np.vstack([self.nodes, np.full((ncoords, 3), -1, dtype=idtype)])

            self.coord_id = np.hstack([self.coord_id, np.full(ncoords, -1, dtype='int32')])
            self.ref_coord_id = np.hstack([self.ref_coord_id, np.full(ncoords, -1, dtype='int32')])

            null = np.full((ncoords, 3), np.nan, dtype='float64')
            self.e1 = np.vstack([self.e1, null])
            self.e2 = np.vstack([self.e2, null])
            self.e3 = np.vstack([self.e3, null])

            self.origin   = np.vstack([self.origin, null])
            self.z_axis   = np.vstack([self.z_axis, null])
            self.xz_plane = np.vstack([self.xz_plane, null])

            #TT = np.array([ones], dtype='float64')
            self.is_resolved = np.hstack([self.is_resolved, np.full(ncoords, False, dtype='bool')])
            self.i = np.vstack([self.i, null])
            self.j = np.vstack([self.j, null])
            self.k = np.vstack([self.k, null])
            self.T = np.vstack([self.T, np.full((ncoords, 3, 3), np.nan, dtype='float64')])
            #self.is_resolved[0] = True
            ncoords_actual = len(self.coord_id)
            assert self.T.shape == (ncoords_actual, 3, 3), self.T.shape
            del ncoords_actual
        else:
            self.coord_type = np.zeros(ncoords, dtype='|U1')
            self.icoord = np.zeros(ncoords, dtype='int8')

            self.nodes = np.full((ncoords, 3), -1, dtype=idtype)
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
            (icoord, coord_type, cid), cardi, ifilei, comment = card
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
            (icoord, coord_type, cid), cardi, ifile, comment = card
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
            #self.slice_card_by_id(ucoords)
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

        rid_to_i_icoord_coordtype = {cid: (i, icoord, str(coord_type)) for i, cid, icoord, coord_type in
                                     zip(count(), self.coord_id, self.icoord, self.coord_type)}

        grid = self.model.grid
        if grid.n:
            grid.sort()
        assert len(np.unique(self.coord_id)) == len(self.coord_id)

        debug = False
        if debug:  # pragma: no cover
            print(f'nresolved = {nresolved}')
        while nresolved < ncoords:
            i1 = np.where(self.icoord == 1)[0]
            i2 = np.where(self.icoord == 2)[0]
            coord2_resolvable_rids = {rid for rid in self.ref_coord_id[i2]
                                      if rid in resolved}
            coord2_cids_to_resolve = [cid for cid, rid in zip(self.coord_id[i2], self.ref_coord_id[i2])
                                      if rid in coord2_resolvable_rids and cid not in resolved]

            coord1_cids_to_resolve = self._find_cord1s_to_resolve(
                grid, i1, resolved)

            assert len(np.unique(self.coord_id)) == len(self.coord_id)
            coords_to_resolve = coord1_cids_to_resolve + coord2_cids_to_resolve
            if len(coords_to_resolve) == 0:
                raise RuntimeError(f'cannot resolve any coordinate systems...unresolved_cids={unresolved_cids}')

            if debug:  # pragma: no cover
                print('-----------------------')
                n_to_resolve = len(coord1_cids_to_resolve) + len(coord2_cids_to_resolve)
                log.debug(f'n={n_to_resolve} resolve this cycle: coord1={np.array(coord1_cids_to_resolve)}; '
                          f'coord2={np.array(coord2_cids_to_resolve)}')
            nresolved1 = nresolved
            nresolved, cord1s_resolved = self._resolve_cord1(
                coord1_cids_to_resolve, nresolved, resolved, grid,
                unresolved_cids)
            if debug:  # pragma: no cover
                if nresolved > nresolved1:
                    log.debug(f'n={len(coord1_cids_to_resolve)}; resolved CORD1x={cord1s_resolved} -> {coord1_cids_to_resolve}')
                    log.debug(f'n={len(unresolved_cids)}; unresolved_cids = {unresolved_cids}')

            nresolved2 = nresolved
            nresolved, cord2s_resolved = self._resolve_cord2(
                coord2_cids_to_resolve, resolved, nresolved, unresolved_cids,
                rid_to_i_icoord_coordtype)

            if debug:  # pragma: no cover
                if nresolved > nresolved2:
                    log.debug(f'n={len(coord2_cids_to_resolve)}; resolved CORD2x={cord2s_resolved} -> {coord2_cids_to_resolve}')
                if unresolved_cids:
                    log.debug(f'n={len(unresolved_cids)}; unresolved_cids = {unresolved_cids}')

        if 0 in resolved:
            # just limiting log messages
            resolved.remove(0)
        if debug and resolved:  # pragma: no cover
            log.info(f'n={len(resolved)}; resolved={np.array(list(resolved))}')
        if unresolved_cids:
            #print(f'unresolved_cids = {unresolved_cids}\n{self.write()}')
            raise RuntimeError(f'unresolved_cids = {unresolved_cids}\n{self.write()}')
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

    def _find_cord1s_to_resolve(self, grid: GRID, i1: np.ndarray,
                                resolved: set[int]) -> list[int]:
        coord1_cids_to_resolve = []
        for cid, nodes in zip(self.coord_id[i1], self.nodes[i1, :]):
            if cid in resolved:
                continue
            inid = grid.index(nodes)
            cps_for_nodes = grid.cp[inid]
            is_resolved = [cpi in resolved for cpi in cps_for_nodes]
            if all(is_resolved):
                coord1_cids_to_resolve.append(cid)
        return coord1_cids_to_resolve

    def _resolve_cord1(self, coord1_cids_to_resolve: list[int],
                       nresolved: int, resolved: set[int], grid: GRID,
                       unresolved_cids: set[int]) -> tuple[int, set[int]]:
        """resolve CORD1R, CORD1S, CORD1C"""
        resolved1 = set()
        if not coord1_cids_to_resolve:
            return nresolved, resolved1
        coord1_cids_to_resolve.sort()
        inids = np.searchsorted(self.coord_id, coord1_cids_to_resolve)
        for cidi, i in zip(coord1_cids_to_resolve, inids):
            cid = self.coord_id[i]
            assert cid == cidi
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
            resolved1.add(cid)
            if cid in unresolved_cids:
                unresolved_cids.remove(cid)
            nresolved += 1
        return nresolved, resolved1

    def _resolve_cord2(self, coord2_cids_to_resolve: list[int],
                       resolved: set[int], nresolved: int,
                       unresolved_cids: set[int],
                       rid_to_i_icoord_coordtype: dict[int, Any]) -> tuple[int, set[int]]:
        """resolve CORD2R, CORD2S, CORD2C"""
        #log = self.model.log
        icid = np.searchsorted(self.coord_id, coord2_cids_to_resolve)
        resolved2 = set()
        if len(icid) == 0:
            return nresolved, resolved2

        for i in icid:
            cid = self.coord_id[i]
            if cid in resolved:
                continue
            rid = self.ref_coord_id[i]
            coord_type = self.coord_type[i]
            icoord = self.icoord[i]
            assert coord_type in {'R', 'C', 'S'}, coord_type
            assert icoord == 2, icoord

            # get info about the reference's reference coord:
            #  - CORD2R -> CORD2R
            #  - CORD2R -> CORD1R
            irid, ricoord, rcoord_type = rid_to_i_icoord_coordtype[rid]

            assert rcoord_type in {'R', 'C', 'S'}, rcoord_type
            assert ricoord in {1, 2}, ricoord

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
            resolved2.add(cid)
            unresolved_cids.remove(cid)
        return nresolved, resolved2

    #@property
    #def cards(self) -> list[Any]:
        #return self.cards1 + self.cards2

    @property
    def max_id(self) -> int:
        return max(self.coord_id.max(), self.nodes.max(), self.ref_coord_id.max())

    @parse_check
    def write_file(self, bdf_file: TextIOLike,
                   size: int=8, is_double: bool=False,
                   write_card_header: bool=False) -> None:
        print_card, size = get_print_card_size(size, self.max_id)

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
                    assert rid != '-1'
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
                    assert rid != '-1'
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

    def equivalence_nodes(self, nid_old_to_new: dict[int, int]) -> None:
        """helper for bdf_equivalence_nodes"""
        nodes = self.nodes.ravel()
        for i, nid1 in enumerate(nodes):
            nid2 = nid_old_to_new.get(nid1, nid1)
            nodes[i] = nid2

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

#def _parse_cord1x(card: BDFCard) -> tuple[int, int, int, int]:
    #cid = integer(card, 1, 'cid')
    #return cid, grid_origin, grid_zaxis, grid_xzplane

def _parse_cord2x(model: BDF,
                  card: BDFCard,
                  card_type: str) -> tuple[int, int,
                                           np.ndarray, np.ndarray, np.ndarray]:
    fdouble_or_blank = double_or_blank if model.is_strict_card_parser else force_double_or_blank

    cid = integer(card, 1, 'cid')

    #: reference coordinate system ID
    rid = integer_or_blank(card, 2, 'rid', default=0)

    #: origin in a point relative to the rid coordinate system
    origin = np.array([fdouble_or_blank(card, 3, 'e1x', default=0.0),
                       fdouble_or_blank(card, 4, 'e1y', default=0.0),
                       fdouble_or_blank(card, 5, 'e1z', default=0.0)],
                      dtype='float64')
    #: z-axis in a point relative to the rid coordinate system
    zaxis = np.array([fdouble_or_blank(card, 6, 'e2x', default=0.0),
                      fdouble_or_blank(card, 7, 'e2y', default=0.0),
                      fdouble_or_blank(card, 8, 'e2z', default=0.0)],
                     dtype='float64')
    #: a point on the xz-plane relative to the rid coordinate system
    xzplane = np.array([fdouble_or_blank(card, 9, 'e3x', default=0.0),
                        fdouble_or_blank(card, 10, 'e3y', default=0.0),
                        fdouble_or_blank(card, 11, 'e3z', default=0.0)],
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
