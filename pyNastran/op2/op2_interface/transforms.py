"""
Defines:
 - transform_displacement_to_global(subcase, result, icd_transform, coords, xyz_cid0,
                                    log, debug=False)
 - transform_gpforce_to_globali(subcase, result,
                                 nids_all, nids_transform,
                                 i_transform, coords, xyz_cid0, log)

"""
from __future__ import annotations
import sys
from typing import Optional, TYPE_CHECKING
import numpy as np

from pyNastran.femutils.coord_transforms import (
    cylindrical_rotation_matrix, spherical_rotation_matrix)
from pyNastran.femutils.matrix3d import (
    dot_n33_33,
    #dot_n33_n33,
    #dot_33_n33,
    dot_n33_n3)
if TYPE_CHECKING:  # pragma: no cover
    from cpylog import SimpleLogger
    from pyNastran.bdf.subcase import Subcase
    from pyNastran.bdf.cards.coordinate_systems import Coord
    from pyNastran.op2.result_objects.table_object import RealTableArray
    from pyNastran.op2.tables.ogf_gridPointForces.ogf_objects import RealGridPointForces

def transform_displacement_global_to_local(subcase: Subcase,
                                           result: RealTableArray,
                                           cid_out: int,
                                           icd_transform: dict[int, np.ndarray],
                                           coords: dict[int, Coord],
                                           xyz_cid0: np.ndarray,
                                           log: SimpleLogger, debug: bool=False) -> None:
    """
    Performs an inplace operation to transform the DISPLACMENT, VELOCITY,
    ACCELERATION result into the global (cid=0) frame

    Assuming you've called ``transform_displacement_to_global`` first,
    set:

    nnodes = xyz_cid0.shape[0]
    icd_transform = {
        0: np.arange(nnodes),
    }

    """
    data = result.data
    nnodesi = data.shape[1]

    coord_out = coords.coord[cid_out]
    assert coord_out.type in {'CORD2R', 'CORD1R'}, str(coord_out)
    coord_out_transform = coord_out.beta().T

    # assert 0 in icd_transform, icd_transform
    for cid, inode in icd_transform.items():
        if cid == -1 or len(inode) == 0:
            continue
        # if cid == 0 and transform_to_global:
        #     continue
        coord = coords[cid]
        coord_type = coord.type
        cid_transform = coord_out_transform

        # a global coordinate system has 1.0 along the main diagonal
        is_global_cid = False
        if np.array_equal([1., 1., 1.], np.diagonal(cid_transform)):
            is_global_cid = True

        if not is_global_cid and debug:
            log.debug('coord\n%s' % coord)
            log.debug(cid_transform)
            log.debug('inode = [%s]' % ', '.join([str(val).rstrip('L') for val in inode.tolist()]))
            log.debug('data.shape = %s' % str(data.shape))
            log.debug('len(inode) = %s' % len(inode))
            assert np.array_equal(inode, np.unique(inode))

        _transform_coord_nids(
            inode, data, nnodesi,
            cid_transform, coord, coord_type, is_global_cid,
            log,
            xyz_cid0=xyz_cid0)

def transform_displacement_to_global(subcase: Subcase,
                                     result: RealTableArray,
                                     icd_transform: dict[int, np.ndarray],
                                     coords: dict[int, Coord],
                                     xyz_cid0: np.ndarray,
                                     log: SimpleLogger, debug: bool=False) -> None:
    """
    Performs an inplace operation to transform the DISPLACMENT, VELOCITY,
    ACCELERATION result into the global (cid=0) frame

    """
    data = result.data
    nnodesi = data.shape[1]

    # if 0 not in icd_transform:
    #     icd_transform[0] = np.zeros(0, dtype=np.float64)

    for cid, inode in icd_transform.items():
        if cid in {-1, 0} or len(inode) == 0:
            continue

        coord = coords[cid]
        coord_type = coord.type
        cid_transform = coord.beta()

        # a global coordinate system has 1.0 along the main diagonal
        is_global_cid = False
        if np.array_equal([1., 1., 1.], np.diagonal(cid_transform)):
            is_global_cid = True

        if not is_global_cid and debug:
            log.debug('coord\n%s' % coord)
            log.debug(cid_transform)
            log.debug('inode = [%s]' % ', '.join([str(val).rstrip('L') for val in inode.tolist()]))
            log.debug('data.shape = %s' % str(data.shape))
            log.debug('len(inode) = %s' % len(inode))
            assert np.array_equal(inode, np.unique(inode))

        _transform_coord_nids(
            inode, data, nnodesi,
            cid_transform, coord, coord_type, is_global_cid,
            log,
            xyz_cid0=xyz_cid0)


def _transform_coord_nids(inode: np.ndarray,
                          data: np.ndarray,
                          nnodesi: int,
                          cid_transform: np.ndarray,
                          coord: Coord,
                          coord_type: str, is_global_cid: bool,
                          log: SimpleLogger,
                          xyz_cid0: Optional[np.ndarray]) -> None:
        if coord_type in ['CORD2R', 'CORD1R']:
            if is_global_cid:
                #self.log.debug('is_global_cid')
                return
            #self.log.debug('rectangular')


            # isat_tran.op2
            #  - nspoint = 4
            #  - ngrid = 5379
            #
            #  - ntotal = 5383
            #  data.shape = (101, 8, 6)
            try:
                translation = data[:, inode, :3]
                rotation = data[:, inode, 3:]
            except IndexError:
                log.warning('shape of inode is incorrect')
                #translation = np.full(inode.shape, np.nan, dtype=data.dtype)
                #rotation = np.full(inode.shape, np.nan, dtype=data.dtype)
                print(data.shape, inode.shape, nnodesi)
                print(inode[:nnodesi].shape)
                return
                #translation = data[:, inode[:nnodesi], :3]
                #rotation = data[:, inode[:nnodesi], 3:]
            data[:, inode, :3] = translation @ cid_transform
            data[:, inode, 3:] = rotation @ cid_transform

        elif coord_type in ['CORD2C', 'CORD1C']:
            #self.log.debug('cylindrical')
            if xyz_cid0 is None:
                msg = 'xyz_cid0 is required for cylindrical coordinate transforms'
                raise RuntimeError(msg)
            _transform_cylindrical_displacement(inode, data, coord, xyz_cid0, cid_transform,
                                                is_global_cid)

        elif coord_type in ['CORD2S', 'CORD1S']:
            #print('spherical')
            if xyz_cid0 is None:
                msg = ('xyz_cid is required for spherical '
                       'coordinate transforms')
                raise RuntimeError(msg)
            _transform_spherical_displacement(inode, data, coord, xyz_cid0, cid_transform,
                                              is_global_cid)
        else:
            raise RuntimeError(coord)


def _transform_cylindrical_displacement(inode: np.ndarray,
                                        data: np.ndarray,
                                        coord: Coord,
                                        xyz_cid0: np.ndarray,
                                        cid_transform: np.ndarray,
                                        is_global_cid: bool) -> None:
    """helper method for transform_displacement_to_global

    Transforms displacement/velocity/acceleration from cylindrical CD
    to global. The theta for each node is computed via
    coord.transform_node_to_local_array (global → coord-local).

    Full transform (column-vector): v_global = beta^T @ R_theta @ v_cyl

    """
    xyzi = xyz_cid0[inode, :]
    rtz_cid = coord.transform_node_to_local_array(xyzi)

    thetad = rtz_cid[:, 1]
    thetar = np.radians(thetad)

    xforms = cylindrical_rotation_matrix(thetar, dtype='float64')

    # Full transform (column-vector): v_global = beta^T @ R_theta @ v_cyl
    # combined[i] = beta^T @ R_theta[i]
    if is_global_cid:
        xforms_combined = xforms
    else:
        beta_T = cid_transform.T
        xforms_combined = np.einsum('ij,njk->nik', beta_T, xforms)

    for itime in range(data.shape[0]):
        translation = data[itime, inode, :3]
        rotation = data[itime, inode, 3:]
        data[itime, inode, :3] = dot_n33_n3(xforms_combined, translation)
        data[itime, inode, 3:] = dot_n33_n3(xforms_combined, rotation)

def _transform_spherical_displacement(inode: np.ndarray,
                                      data: np.ndarray,
                                      coord: Coord,
                                      xyz_cid0: np.ndarray,
                                      cid_transform: np.ndarray,
                                      is_global_cid: bool) -> None:
    """helper method for transform_displacement_to_global

    Transforms displacement/velocity/acceleration from spherical CD
    to global. The rotation matrix depends on each node's (theta, phi)
    in the coord's local frame.

    Full transform (column-vector): v_global = beta^T @ R_sph @ v_sph

    """
    xyzi = xyz_cid0[inode, :]
    rtp_cid = coord.transform_node_to_local_array(xyzi)

    thetar = np.radians(rtp_cid[:, 1])
    phir = np.radians(rtp_cid[:, 2])

    xforms = spherical_rotation_matrix(thetar, phir, dtype='float64')

    # Full transform (column-vector): v_global = beta^T @ R_sph @ v_sph
    if is_global_cid:
        xforms_combined = xforms
    else:
        beta_T = cid_transform.T
        xforms_combined = np.einsum('ij,njk->nik', beta_T, xforms)

    for itime in range(data.shape[0]):
        translation = data[itime, inode, :3]
        rotation = data[itime, inode, 3:]
        data[itime, inode, :3] = dot_n33_n3(xforms_combined, translation)
        data[itime, inode, 3:] = dot_n33_n3(xforms_combined, rotation)

def transform_gpforce_to_globali(subcase: Subcase,
                                 result: RealGridPointForces,
                                 nids_all: np.ndarray,
                                 nids_transform: np.ndarray,
                                 i_transform: dict[int, np.ndarray],
                                 coords: dict[int, Coord],
                                 xyz_cid0: np.ndarray,
                                 log: SimpleLogger):
    """
    Performs an inplace operation to transform the GPFORCE result
    into the global (cid=0) frame

    """
    log.debug('result.name = %s' % result.class_name)
    data = result.data

    # inode_xyz :
    #    the indices of the nodes in the model grid point list
    for cid, inode_xyz in i_transform.items():
        log.debug('cid = %s' % cid)
        if cid in [-1, 0]:
            continue
        coord = coords[cid]
        coord_type = coord.type
        cid_transform = coord.beta()

        # a global coordinate system has 1s along the main diagonal
        is_global_cid = False
        if np.array_equal([1., 1., 1.], np.diagonal(cid_transform)):
            is_global_cid = True
        nids = np.array(nids_transform[cid])

        #from pyNastran.op2.tables.ogf_gridPointForces.ogf_objects import RealGridPointForcesArray
        #result = RealGridPointForcesArray()
        if result.is_unique: # TODO: doesn't support preload
            #self.node_element = np.zeros((self.ntimes, self.ntotal, 2), dtype='int32')
            nids_all_gp = result.node_element[0, :, 0]

            log.debug('nids_all_gp = %s' % list(nids_all_gp))
            log.debug('nids = %s' % list(nids))

            # the indices of the grid points that we're transforming
            inode_gp = np.where(np.isin(nids_all_gp, nids))[0]
            log.debug('inode_gp = %s' % list(inode_gp))
            nids_gp = nids_all_gp[inode_gp]
            log.debug('nids_gp = %s' % list(nids_gp))
            if not np.array_equal(np.unique(nids_gp), np.unique(nids)):
                msg = 'nids_gp=%s nids=%s' % (nids_gp, nids)
                raise RuntimeError(msg)

            log.debug('---------')
            # the transformation index to go from xyz to the grid point forces
            inode_gp_xyz = np.searchsorted(nids_all, nids_all_gp)
            log.debug('nids_all = %s' % list(nids_all))
            log.debug('nids_all_gp = %s' % list(nids_all_gp))
            log.debug('inode_xyz = %s' % list(inode_xyz))
            log.debug('inode_gp_xyz = %s' % list(inode_gp_xyz))
            sys.stdout.flush()
            assert len(inode_gp_xyz) == len(nids_all_gp), len(nids_all_gp)
        else:
            raise NotImplementedError(result)

        if coord_type in ['CORD2R', 'CORD1R']:
            if is_global_cid:
                continue
            #print('coord\n', coord)
            #print(cid_transform)
            #print('inode = %s' % inode)
            #print('rectangular')
            translation = data[:, inode_gp, :3]
            rotation = data[:, inode_gp, 3:]
            data[:, inode_gp, :3] = translation @ cid_transform
            data[:, inode_gp, 3:] = rotation @ cid_transform
        elif coord_type in ['CORD2C', 'CORD1C']:
            #print('cylindrical')
            if xyz_cid0 is None:
                msg = ('xyz_cid is required for cylindrical '
                       'coordinate transforms')
                raise RuntimeError(msg)
            _transform_cylindrical_gpforce(inode_gp_xyz, inode_gp, data,
                                           cid_transform, coord,
                                           xyz_cid0, log)

        elif coord_type in ['CORD2S', 'CORD1S']:
            log.debug('spherical')
            if xyz_cid0 is None:
                msg = ('xyz_cid is required for spherical '
                       'coordinate transforms')
                raise RuntimeError(msg)

            _transform_spherical_gpforce(inode_gp_xyz, inode_gp, data,
                                         cid_transform, coord,
                                         xyz_cid0, log)
        else:  # pragma: no cover
            raise RuntimeError(coord)
    return

def _transform_cylindrical_gpforce(inode_xyz: np.ndarray,
                                   inode_gp: np.ndarray,
                                   data: np.ndarray,
                                   cid_transform: np.ndarray,
                                   coord: Coord,
                                   xyz_cid0: np.ndarray,
                                   log: SimpleLogger) -> None:
    """Transform grid point forces from cylindrical CD to global.

    In a cylindrical coordinate system, the radial/tangential directions
    are position-dependent (vary with theta). Each GPF entry must be rotated
    by the theta of its grid point before applying the coord-to-global transform.

    Full transform (column-vector): v_global = beta^T @ R_theta @ v_cyl

    """
    xyz_indices = inode_xyz[inode_gp]
    xyzi = xyz_cid0[xyz_indices, :]
    rtz_cid = coord.transform_node_to_local_array(xyzi)
    thetad = rtz_cid[:, 1]
    thetar = np.radians(thetad)

    # R_theta: (n, 3, 3) — rotates cylindrical (r,t,z) to coord rectangular
    xforms = cylindrical_rotation_matrix(thetar, dtype='float64')

    # Full transform (column-vector): v_global = beta^T @ R_theta @ v_cyl
    # combined[i] = beta^T @ R_theta[i]
    is_global_cid = np.array_equal([1., 1., 1.], np.diagonal(cid_transform))
    if is_global_cid:
        xforms_combined = xforms
    else:
        beta_T = cid_transform.T
        xforms_combined = np.einsum('ij,njk->nik', beta_T, xforms)

    for itime in range(data.shape[0]):
        translation = data[itime, inode_gp, :3]
        rotation = data[itime, inode_gp, 3:]
        data[itime, inode_gp, :3] = dot_n33_n3(xforms_combined, translation)
        data[itime, inode_gp, 3:] = dot_n33_n3(xforms_combined, rotation)

def _transform_spherical_gpforce(inode_xyz: np.ndarray,
                                 inode_gp: np.ndarray,
                                 data: np.ndarray,
                                 cid_transform: np.ndarray,
                                 coord: Coord,
                                 xyz_cid0: np.ndarray,
                                 log: SimpleLogger) -> None:
    """Transform grid point forces from spherical CD to global.

    Full transform (column-vector): v_global = beta^T @ R_sph @ v_sph

    """
    xyz_indices = inode_xyz[inode_gp]
    xyzi = xyz_cid0[xyz_indices, :]
    rtp_cid = coord.transform_node_to_local_array(xyzi)

    thetar = np.radians(rtp_cid[:, 1])
    phir = np.radians(rtp_cid[:, 2])

    xforms = spherical_rotation_matrix(thetar, phir, dtype='float64')

    is_global_cid = np.array_equal([1., 1., 1.], np.diagonal(cid_transform))
    if is_global_cid:
        xforms_combined = xforms
    else:
        beta_T = cid_transform.T
        xforms_combined = np.einsum('ij,njk->nik', beta_T, xforms)

    for itime in range(data.shape[0]):
        translation = data[itime, inode_gp, :3]
        rotation = data[itime, inode_gp, 3:]
        data[itime, inode_gp, :3] = dot_n33_n3(xforms_combined, translation)
        data[itime, inode_gp, 3:] = dot_n33_n3(xforms_combined, rotation)
