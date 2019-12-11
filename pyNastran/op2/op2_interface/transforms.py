"""
Defines:
 - transform_displacement_to_global(subcase, result, icd_transform, coords, xyz_cid0,
                                    log, debug=False)
 - transform_gpforce_to_globali(subcase, result,
                                 nids_all, nids_transform,
                                 i_transform, coords, xyz_cid0, log)

"""
import sys
import numpy as np

from pyNastran.femutils.coord_transforms import cylindrical_rotation_matrix
from pyNastran.femutils.matrix3d import (
    dot_n33_33,
    #dot_n33_n33,
    #dot_33_n33,
    dot_n33_n3)


def transform_displacement_to_global(subcase, result, icd_transform, coords, xyz_cid0,
                                     log, debug=False):
    """
    Performs an inplace operation to transform the DISPLACMENT, VELOCITY,
    ACCELERATION result into the global (cid=0) frame

    """
    #print('result.name = ', result.class_name)
    data = result.data
    nnodesi = data.shape[1]
    for cid, inode in icd_transform.items():
        if cid in [-1, 0]:
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

        if coord_type in ['CORD2R', 'CORD1R']:
            if is_global_cid:
                #self.log.debug('is_global_cid')
                continue
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
                translation = np.full(inode.shape, np.nan, dtype=data.dtype)
                rotation = np.full(inode.shape, np.nan, dtype=data.dtype)
                print(data.shape, inode.shape, nnodesi)
                print(inode[:nnodesi].shape)
                continue
                translation = data[:, inode[:nnodesi], :3]
                rotation = data[:, inode[:nnodesi], 3:]
            data[:, inode, :3] = translation.dot(cid_transform)
            data[:, inode, 3:] = rotation.dot(cid_transform)

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


def _transform_cylindrical_displacement(inode, data, coord, xyz_cid0, cid_transform, is_global_cid):
    """helper method for transform_displacement_to_global"""
    xyzi = xyz_cid0[inode, :]
    rtz_cid = coord.xyz_to_coord_array(xyzi)

    #print('rtz_cid:')
    #print(rtz_cid)
    thetad = rtz_cid[:, 1]
    thetar = np.radians(thetad)
    #print('thetad[cid=%s] = %s' % (coord.cid, thetad))
    #print('thetar = ', thetar)
    #self.log.debug('thetad = %s' % list(thetad))

    #np.set_printoptions(precision=None, threshold=None, edgeitems=None, linewidth=None, suppress=True,
                        #nanstr=None, infstr=None, formatter=None, sign=None, floatmode=None)
    xforms = cylindrical_rotation_matrix(thetar, dtype='float64')
    #print('xforms.shape = %s' % str(xforms.shape))

    for itime in range(data.shape[0]):
        translation = data[itime, inode, :3]
        rotation = data[itime, inode, 3:]
        #theta_max1 = translation[:, 1].max()
        #theta_max2 = rotation[:, 1].max()
        #print('theta_max = ', max(theta_max1, theta_max2))

        if is_global_cid:
            translation2 = dot_n33_n3(xforms, translation)
            rotation2 = dot_n33_n3(xforms, rotation)
        else:
            #print('cid_transform[cid=%s]:\n%s' % (coord.cid, cid_transform))
            #print('xforms:\n%s' % (xforms))
            #xforms2 = dot_33_n33(xforms, cid_transform)  # wrong shape
            #xforms2a = dot_33_n33(cid_transform, xforms)
            #xforms2b = dot_n33_33(xforms, cid_transform)
            #xforms2b = xforms
            # triple products
            #xforms2c = dot_33_n33(cid_transform.T, dot_n33_33(xforms, cid_transform))
            #xforms2d = dot_33_n33(cid_transform, dot_n33_33(xforms, cid_transform.T))

            #print('xform A:\n%s'  % xforms2a)
            #print('xform B:\n%s'  % xforms2b)
            #translation2 = dot_n33_n3(xforms2a, translation)
            #xforms2 = xforms2b
            if 1:
                xforms2b = dot_n33_33(xforms, cid_transform)
                translation2 = dot_n33_n3(xforms2b, translation)
                rotation2 = dot_n33_n3(xforms2b, rotation)
            elif 0:  ## pragma: no cover
                xforms2b = dot_n33_33(xforms, cid_transform)
                translation2a = dot_n33_n3(xforms2b, translation)
                rotation2a = dot_n33_n3(xforms2b, rotation)
                print('translation2a.shape =', translation2a.shape)
                assert translation.shape == translation2a.shape

                translation2 = translation2a.dot(cid_transform)
                rotation2 = rotation2a.dot(cid_transform)

            elif 0:  ## pragma: no cover
                # bad shape
                translation2a = dot_n33_n3(xforms, translation)
                rotation2a = dot_n33_n3(xforms, rotation)
                print('translation2a.shape =', translation2a.shape)
                assert translation.shape == translation2a.shape

                translation2 = cid_transform.dot(translation2a.dot(cid_transform))
                rotation2 = cid_transform.dot(rotation2a.dot(cid_transform))
            else:  ## pragma: no cover
                raise NotImplementedError('pick an option...')
        #print('translation.shape', translation.shape)
        #print('translation2.shape', translation2.shape)
        data[itime, inode, :3] = translation2
        data[itime, inode, 3:] = rotation2

def _transform_spherical_displacement(inode, data, coord, xyz_cid0, cid_transform, is_global_cid):
    """helper method for transform_displacement_to_global"""
    xyzi = xyz_cid0[inode, :]
    #_transform_spherical_displacement_func
    #rtp_cid = coord.xyz_to_coord_array(xyzi)
    #theta = rtp_cid[:, 1]
    #phi = rtp_cid[:, 2]
    for itime in range(data.shape[0]):
        translation = data[itime, inode, :3]
        rotation = data[itime, inode, 3:]
        #if 0:
            #translation[:, 1] += theta
            #translation[:, 2] += phi
            #rotation[:, 1] += theta
            #rotation[:, 2] += phi
        translation = coord.coord_to_xyz_array(translation)
        rotation = coord.coord_to_xyz_array(rotation)
        if is_global_cid:
            data[itime, inode, :3] = translation
            data[itime, inode, 3:] = rotation
            return
        data[itime, inode, :3] = translation.dot(cid_transform)
        data[itime, inode, 3:] = rotation.dot(cid_transform)

def transform_gpforce_to_globali(subcase, result,
                                 nids_all, nids_transform,
                                 i_transform, coords, xyz_cid0, log):
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
            #self.node_element = zeros((self.ntimes, self.ntotal, 2), dtype='int32')
            nids_all_gp = result.node_element[0, :, 0]

            log.debug('nids_all_gp = %s' % list(nids_all_gp))
            log.debug('nids = %s' % list(nids))

            # the indices of the grid points that we're transforming
            inode_gp = np.where(np.in1d(nids_all_gp, nids))[0]
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
            data[:, inode_gp, :3] = translation.dot(cid_transform)
            data[:, inode_gp, 3:] = rotation.dot(cid_transform)
        elif coord_type in ['CORD2C', 'CORD1C']:
            #print('cylindrical')
            if xyz_cid0 is None:
                msg = ('xyz_cid is required for cylindrical '
                       'coordinate transforms')
                raise RuntimeError(msg)
            #_transform_cylindrical_displacement(inode, data, coord, xyz_cid0, cid_transform)
            _transform_cylindrical_gpforce(inode_xyz, inode_gp, data, cid_transform, coord,
                                           xyz_cid0, log)

        elif coord_type in ['CORD2S', 'CORD1S']:
            log.debug('spherical')
            if xyz_cid0 is None:
                msg = ('xyz_cid is required for spherical '
                       'coordinate transforms')
                raise RuntimeError(msg)

            _transform_spherical_gpforce(inode_xyz, inode_gp, data, cid_transform, coord,
                                         xyz_cid0, log)
        else:
            raise RuntimeError(coord)
    return

def _transform_cylindrical_gpforce(unused_inode_xyz, inode_gp, data, cid_transform, coord,
                                   xyz_cid0, log):
    """helper method for transform_gpforce_to_globali"""
    #xyzi = xyz_cid0[inode_xyz, :]s
    #rtz_cid = coord.xyz_to_coord_array(xyzi)
    #theta_xyz = rtz_cid[inode_xyz, 1]
    #theta = rtz_cid[inode_gp_xyz, 1]
    log.debug('coord\n%s' % coord)
    log.debug(cid_transform)
    #print('theta_xyz = %s' % list(theta_xyz))
    #print('theta     = %s' % list(theta))

    for itime in range(data.shape[0]):
        #inode = np.where(np.in1d(nids_all, 2))[0]
        #print('start', data[itime, inode_gp, :3])
        translation = data[itime, inode_gp, :3]
        rotation = data[itime, inode_gp, 3:]
        if 0:
            # horrible
            translation = coord.coord_to_xyz_array(data[itime, inode_gp, :3])
            rotation = coord.coord_to_xyz_array(data[itime, inode_gp, 3:])
            data[itime, inode_gp, :3] = translation.dot(cid_transform)
            data[itime, inode_gp, 3:] = rotation.dot(cid_transform)
        #elif 0:
            ## expectedly horrible
            #translation[:, 1] += theta
            #rotation[:, 1] += theta
            #translation = coord.coord_to_xyz_array(data[itime, inode_gp, :3])
            #rotation = coord.coord_to_xyz_array(data[itime, inode_gp, 3:])
            #data[itime, inode_gp, :3] = translation.dot(cid_transform)
            #data[itime, inode_gp, 3:] = rotation.dot(cid_transform)
        #elif 0:
            ## actually not bad...
            #translation = coord.xyz_to_coord_array(translation)
            #rotation = coord.xyz_to_coord_array(rotation)
            #translation[:, 1] += theta
            #rotation[:, 1] += theta
            #translation = coord.coord_to_xyz_array(translation)
            #rotation = coord.coord_to_xyz_array(rotation)
            #data[itime, inode_gp, :3] = translation.dot(cid_transform)
            #data[itime, inode_gp, 3:] = rotation.dot(cid_transform)
        #elif 0:
            ## not that bad, worse than previous
            #translation = coord.xyz_to_coord_array(translation)
            #rotation = coord.xyz_to_coord_array(rotation)
            #translation[:, 1] += theta
            #rotation[:, 1] += theta
            #translation = coord.coord_to_xyz_array(translation)
            #rotation = coord.coord_to_xyz_array(rotation)
            #data[itime, inode_gp, :3] = translation.dot(cid_transform.T)
            #data[itime, inode_gp, 3:] = rotation.dot(cid_transform.T)
        elif 1:
            # doesn't work...actually pretty close
            data[itime, inode_gp, :3] = translation.dot(cid_transform)
            data[itime, inode_gp, 3:] = rotation.dot(cid_transform)
        #elif 0:
            ## very, very close
            #data[itime, inode_gp, :3] = translation.dot(cid_transform.T)
            #data[itime, inode_gp, 3:] = rotation.dot(cid_transform.T)
        #elif 0:
            ## is this just the same as one of the previous?
            #data[itime, inode_gp, :3] = cid_transform.T.dot(translation.T).T
            #data[itime, inode_gp, 3:] = cid_transform.T.dot(rotation.T).T
        else:
            raise RuntimeError('no option selected...')

        #if is_global_cid:
            #data[itime, inode, :3] = translation
            #data[itime, inode, 3:] = rotation
            #continue
        #data[itime, inode_gp, :3] = translation.dot(cid_transform)
        #data[itime, inode_gp, 3:] = rotation.dot(cid_transform)
        #print('end', data[itime, inode_gp, :3])

def _transform_spherical_gpforce(inode_xyz, unused_inode_gp, data, cid_transform, coord,
                                 unused_xyz_cid0, unused_log):
    """helper method for transform_gpforce_to_globali"""
    #xyzi = xyz_cid0[inode_xyz, :]
    #rtp_cid = coord.xyz_to_coord_array(xyzi)
    #theta = rtp_cid[:, 1]
    #phi = rtp_cid[:, 2]
    for itime in range(data.shape[0]):
        translation = data[itime, inode_xyz, :3]
        rotation = data[itime, inode_xyz, 3:]
        #if 0:
            #translation[:, 1] += theta
            #translation[:, 2] += phi
            #rotation[:, 1] += theta
            #rotation[:, 2] += phi
        translation = coord.coord_to_xyz_array(data[itime, inode_xyz, :3])
        rotation = coord.coord_to_xyz_array(data[itime, inode_xyz, 3:])
        #if is_global_cid:
            #data[itime, inode_xyz, :3] = translation
            #data[itime, inode_xyz, 3:] = rotation
            #continue
        data[itime, inode_xyz, :3] = translation.dot(cid_transform)
        data[itime, inode_xyz, 3:] = rotation.dot(cid_transform)
