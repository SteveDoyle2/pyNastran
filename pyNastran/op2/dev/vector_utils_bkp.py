import numpy as np
from numpy import unique


def transform_force_moment(force_in_local, moment_in_local,
                           coord_out, coords,
                           nid_cd, i_transform,
                           xyz_cid0, summation_point_cid0=None,
                           consider_rxf=True,
                           debug=False, log=None):
    """
    Transforms force/moment from global to local and returns all the forces.

    Parameters
    ----------
    force : (N, 3) ndarray
        forces in the local frame
    moment : (N, 3) ndarray
        moments in the local frame
    coord_out : CORD()
        the desired local frame
    xyz_cid0 : (n, 3) ndarray
        the nodes in the global frame
    summation_point_cid0 : (3, ) ndarray
        the summation point in the global frame
    consider_rxf : bool; default=True
        considers the r x F term
    debug : bool; default=False
        debugging flag
    log : log; default=None
        a logger object that gets used when debug=True

    Returns
    -------
    force_out : (n, 3) float ndarray
        the ith float components in the coord_out coordinate frame
    moment_out : (n, 3) float ndarray
        the ith moment components about the summation point in the coord_out coordinate frame

    .. todo:: doesn't seem to handle cylindrical/spherical systems

    https://flexiblelearning.auckland.ac.nz/sportsci303/13_3/forceplate_manual.pdf

    xyz0 = T_1_to_0 @ xyz1
    xyz1 = T_1_to_0.T @ xyz0

    xyz2 = T_2_to_0.T @ xyz0
    xyz2 = T_2_to_0.T @ T_1_to_0 @ xyz1


    xyz_g = T_a2g @ xyz_a
    xyz_g = T_b2g @ xyz_b
    T_b2g @ xyz_b = T_a2g @ xyz_a
    xyz_b = T_b2g.T @ T_a2g @ xyz_a = T_g2b @ T_a2g @ xyz_a
    """
    debug = True
    assert logger is not None
    #print('nid_cd*\n', nid_cd)
    dtype = force_in_local.dtype
    #dtype = 'float64'
    force_out = np.zeros(force_in_local.shape, dtype=dtype)
    moment_out = np.zeros(force_in_local.shape, dtype=dtype)

    #assert coord_out.Type == 'R', 'Only rectangular coordinate systems are supported'

    nids = nid_cd[:, 0]
    cds = nid_cd[:, 1]
    ucds = unique(cds)

    coord_out_cid = coord_out.cid
    coord_out_T = coord_out.beta()

    if debug:
        log.debug(coord_out)
        if consider_rxf:
            for ii in range(xyz_cid0.shape[0]):
                log.debug('***i=%s xyz=%s nid=%s cd=%s' % (
                    ii, xyz_cid0[ii, :], nid_cd[ii, 0], nid_cd[ii, 1]))
        log.debug('------------')

    if consider_rxf and summation_point_cid0 is None:
        summation_point_cid0 = np.array([0., 0., 0.])

    for cd in ucds:
        i = np.where(cds == cd)[0]
        nidsi = nids[i]
        analysis_coord = coords[cd]
        if debug:
            log.debug('i = %s' % i)
            log.debug('nids = %s' % nidsi)
        #summation_point_cid0 = analysis_coord.origin
        #print('summation_point =', summation_point_cid0)
        #print(analysis_coord)
        cd_T = analysis_coord.beta()
        #print(cd_T)

        force_in_locali = force_in_local[i, :]
        moment_in_locali = moment_in_local[i, :]
        force_in_locali.astype('float64')
        moment_in_locali.astype('float64')

        if debug:
            log.debug('force_input =%s' % force_in_locali)

        #force_in_locali = analysis_coord.coord_to_xyz_array(force_in_locali)
        #moment_in_locali = analysis_coord.coord_to_xyz_array(moment_in_locali)

        #print('coord_out_T\n', coord_out_T)

        # is this:
        #  1.  cd_T.T and coord_out_T
        #  2.  cd_T and coord_out_T.T
        #  3.  reverse it all...
        if 0:
            force_in_globali = cd_T @ force_in_locali.T
            moment_in_globali = cd_T @ moment_in_locali.T
            force_outi = (coord_out_T.T @ force_in_globali).T
            moment_outi = (coord_out_T.T @ moment_in_globali).T
        elif 0:
            force_in_globali = cd_T.T @ force_in_locali.T
            moment_in_globali = cd_T.T @ moment_in_locali.T
            force_outi = (coord_out_T @ force_in_globali).T
            moment_outi = (coord_out_T @ moment_in_globali).T
        else:
            # old
            force_in_globali = (force_in_locali @ cd_T).T
            moment_in_globali = (moment_in_locali @ cd_T).T
            #force_outi = (coord_out_T @ force_in_globali).T
            #moment_outi = (coord_out_T @ moment_in_globali).T
            force_outi = (coord_out_T @ force_in_globali).T
            moment_outi = (coord_out_T @ moment_in_globali).T

        if debug:
            log.debug('force_in_locali =%s' % force_in_locali.T)
            log.debug('force_in_globali = %s' % force_in_globali.T)
            log.debug('force_outi = %s' % force_outi)
            #log.debug('moment_outi = %s' % moment_outi)

        force_out[i, :] = force_outi
        moment_out[i, :] = moment_outi

        # rxF from local_in to global to local_out
        if consider_rxf:
            delta = xyz_cid0[i, :] - summation_point_cid0[np.newaxis, :]
            #delta_minus_origin = delta - analysis_coord.origin
            #delta_in_cid = dot(delta_minus_origin, cd_T.T)
            if debug:
                log.debug('delta = %s' % delta)
            rxf = np.cross(delta, force_in_globali.T)

            if debug:
                log.debug('rxf = %s' % rxf)
            rxf_in_cid = (coord_out_T @ rxf.T).T
            if debug:
                log.debug('rxf_in_cid = %s' % rxf_in_cid)

            moment_out[i, :] += rxf_in_cid
        #print('------------')
    #force_out = coord_out.xyz_to_coord_array(force_out)
    #moment_out = coord_out.xyz_to_coord_array(moment_out)
    return -force_out, -moment_out
