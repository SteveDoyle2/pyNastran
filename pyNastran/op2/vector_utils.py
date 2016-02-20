from __future__ import print_function
from struct import calcsize

import numpy as np
from numpy import array, where, zeros, asarray, dot, arccos, sqrt, pi
from numpy import cos, unique, cross, abs as npabs


def iformat(format_old, precision=2):
    """
    Converts binary data types to size vector arrays.

    Parameters
    ----------
    format_old : str
        the int/float data types in single precision format
    precision : int; default=2
        the precision to convert to
        1 : single precision (no conversion)
        2 : double precision

    Returns
    -------
    format_new : str
        the int/float data types in single/double precision format

    Example
    -------
    >>> iformat('8i6f10s', precision=1)
    '8i6f10s'
    >>> iformat('8i6f10s', precision=2)
    '8l6d10q'
    """
    if precision == 2:  # double
        format_new = format_old.replace('i', 'l').replace('f', 'd')
    elif precision == 1: # single
        format_new = format_old.replace('l', 'i').replace('d', 'f')
    else:
        raise NotImplementedError(precision)
    ndata = calcsize(format_new)
    return format_new, ndata


def abs_max_min_global(values):
    """
    This is useful for figuring out absolute max or min principal stresses
    across single/multiple elements and finding a global max/min value.

    Parameters
    ----------
    values: ndarray/listtuple
        an ND-array of values;
        common NDARRAY/list/tuple shapes:
            1. [nprincipal_stresses]
            2. [nelements, nprincipal_stresses]

    Returns
    -------
    abs_max_mins: int/float
        an array of the max or min principal stress
        don't input mixed types

    nvalues >= 1
      >>> element1 = [0.0, -1.0, 2.0]  # 2.0
      >>> element2 = [0.0, -3.0, 2.0]  # -3.0
      >>> values = abs_max_min_global([element1, element2])
      >>> values
      -3.0

      >>> element1 = [0.0, -1.0, 2.0]  # 2.0
      >>> values = abs_max_min_global([element1])
      >>> values
      2.0

    .. note:: [3.0,  2.0, -3.0] will return 3.0, and
             [-3.0, 2.0,  3.0] will return 3.0
    """
    # support lists/tuples
    values = asarray(values)

    # find the [max,
    #           min]
    # we organize it as [max, min], which is why the note applies
    # we could make both of the edge cases return -3.0, but if you're using
    # this function it shouldn't matter
    values2 = array([values.max(),
                     values.min()])

    # we figure out the absolute max/min
    abs_vals = np.abs(values2)
    abs_val = abs_vals.max()

    # find the location of the absolute max value
    # 1.  we take the first value (the where[0]) to chop the return value
    #     since there is no else conditional
    # 2.  we take the first value (the where[0][0]) to only get the max
    #     value if 2+ values are returned
    j = where(abs_val == abs_vals)[0][0]

    # get the raw value from the absoluted value, so:
    # value = npabs(raw_value)
    return values2[j]


def abs_max_min_vector(values):
    """
    This is useful for figuring out principal stresses across multiple
    elements.

    Parameters
    ----------
    values: ndarray/listtuple
        an array of values, where the rows are interated over
        and the columns are going to be compressed

        common NDARRAY/list/tuple shapes:
            1. [nprincipal_stresses]
            2. [nelements, nprincipal_stresses]

    Returns
    -------
    abs_max_mins: NDARRAY shape=[nelements] with dtype=values.dtype
        an array of the max or min principal stress
        don't input mixed types

    ::
       >>> element1 = [0.0,  1.0, 2.0]  # 2.0
       >>> element2 = [0.0, -1.0, 2.0]  # 2.0
       >>> element3 = [0.0, -3.0, 2.0]  # -3.0
       >>> values = [element1 element2, element3]
       >>> values0 = abs_max_min_vectorized(values)
       >>> values0
       [2.0, 2.0, -3.0]

    .. note:: [3.0,  2.0, -3.0] will return 3.0, and
              [-3.0, 2.0,  3.0] will return 3.0
    """
    # support lists/tuples
    values = asarray(values)

    # find the [maxs,
    #           mins]
    # we organize it as [maxs, mins], which is why the note applies
    # we could make both of the edge cases return -3.0, but if you're using
    # this function it shouldn't matter
    maxs_mins = array([values.max(axis=1),
                       values.min(axis=1)])

    # we figure out the absolute max/min for each row
    abs_vals = npabs(maxs_mins)
    absolute_maxs = abs_vals.max(axis=0)

    outs = zeros(absolute_maxs.shape[0], dtype=values.dtype)
    for i, absolute_max in enumerate(absolute_maxs):
        # find the location of the absolute max value
        # 1.  we take the first value (the where[0]) to chop the return value
        #     since there is no else conditional
        # 2.  we take the first value (the where[0][0]) to only get the max
        #     value if 2+ values are returned
        j = where(absolute_max == abs_vals[:, i])[0][0]

        # get the raw value from the absoluted value, so:
        # value = npabs(raw_value)
        outs[i] = maxs_mins[j, i]
    return outs


def abs_max_min(values, global_abs_max=True):
    if global_abs_max:
        return abs_max_min_global(values)
    return abs_max_min_vector(values)


# def principal_2d(o11, o22, o12):
    # oxx = 5
    # return oxx, oy


def principal_3d(o11, o22, o33, o12, o23, o13):
    """http://www.continuummechanics.org/cm/principalstrain.html"""
    # e = a
    i1 = o11 + o22 + o33
    i2 = o11*o22 + o22*o33 + o11*o33 - o12**2 - o13**2 - o23**2
    i3 = o11*o22*o33 - o11*o23**2 - o22*o13**2 + 2*o12*o13*o23
    Q = 3 * i2 - i1**2
    R = (2*i1**3 - 9*i1*i2 + 27*i3) / 54.
    theta = arccos(R / sqrt(-Q**3))

    q2 = 2 * sqrt(-Q)
    i13 = 1./3. * i1
    p1 = q2 * cos(theta/3) + i13
    p2 = q2 * cos(theta/3 + 2*pi/3.) + i13
    p3 = q2 * cos(theta/3 + 4*pi/3.) + i13
    max_min_mid = array([p1, p2, p3])
    pmax = max_min_mid.max(axis=0)
    pmin = max_min_mid.min(axis=0)
    return pmax, pmin


def test_abs_max_min_global():
    #print(iformat('4si3f', 2))
    print(abs_max_min_global([0.0, 2.0, 1.0]))
    print(abs_max_min_global([0.0, 2.0, -1.0]))
    print(abs_max_min_global([0.0, 2.0, -3.0]))
    print(abs_max_min_global(array([0.0, 2.0, -3.0])))
    print(abs_max_min_global([1.0]))

    # gets the global max/min value
    print(abs_max_min_global([
        [0.0, 2.0, -3.0],
        [0.0, 2.0, -4.0],
    ]))
    print(abs_max_min_global(array([
        [0.0, 2.0, -3.0],
        [0.0, 2.0, -4.0],
    ])))


def test_abs_max_min_vector():
    print(abs_max_min_vector(array([
        [0.0, 2.0, 1.0],
        [0.0, 2.0, -1.0],
        [0.0, 2.0, -3.0],
    ])))

    print(abs_max_min_vector([
        [0.0, 2.0, 1.0],
        [0.0, 2.0, -1.0],
        [0.0, 2.0, -3.0],
        [0.0, 2.0, 4.0],
    ]))
    print(abs_max_min_vector(array([
        [0.0, 2.0, 1.0],
        [0.0, 2.0, -1.0],
        [0.0, 2.0, -3.0],
        [0.0, 2.0, 4.0],
    ])))

    print(abs_max_min_vector(array([
        [3.0, 2.0, -3.0],
        [-3.0, 2.0, 3.0],
    ])))

    # not an array
    #print(abs_max_min([
        #[0.0, 2.0, 1.0],
        #[0.0, 2.0, -1.0],
        #[0.0, 2.0, -3.0],
        #[0.0, 2.0, 4.0],
    #]))


def transform_force_from_local_to_global(force_in_local, gpforce_nids,
                                         nid_cd, i_transform, beta_transforms):
    """
    Transforms force (not moment) from global to local.

    Supports cylindrical/spherical coordinate systems.

    Parameters
    ----------
    force_in_local : (N, 3) ndarray
        forces in the mixed CD (output) coordinate frame
    gpforce_nids : (N, ) int ndarray
        node ids (no unique assumption)
    nid_cd : (Nn, 2) int ndarray
        [node_id, cp] or [node_id, cd] ndarray
    i_transform : dict[int] = (Nn, ) int ndarray
        key : cp / cd
        value : indices of node ids corresponding to nid_cd
    beta_transforms : dict[int] = (3, 3) float ndarray
        key : cp / cd
        value : transformation matrix

    .. warning:: the function signature will change...
    .. todo:: doesn't support origin shifts
    .. todo:: doesn't support moment
    .. todo:: doesn't support coordinate transforms (e.g. rectangular to cylindrical)
    .. todo:: sum of moments about a point must have an rxF term to get the
               same value as Patran.
    """
    force_in_global = zeros(force_in_local.shape, dtype=force_in_local.dtype)
    nids = nid_cd[:, 0]
    cds = nid_cd[:, 1]
    ucds = unique(cds)
    for cd in ucds:
        if cd == 0:
            i = where(cds == cd)[0]
            nidsi = nids[i]
            #force_in_global[:, :3] = force_in_local[i, :3]
            #force_in_global[:, 3:] = force_in_local[i, 3:]
        else:
            raise NotImplementedError(ucds)
    return force_in_local


def transform_force_from_global_to_local(force_in_global, coord_in, coord_out):
    """
    Transforms force (not moment) from global to local.

    Supports cylindrical/spherical coordinate systems.

    Parameters
    ----------
    force_in_global : (N, 3) ndarray
        forces in the global frame
    coord_in : CORD()
        the global coordinate system object
    coord_out : CORD()
        the desired local frame

    .. warning:: the function signature will change...
    .. todo:: doesn't support origin shifts
    .. todo:: doesn't support moment
    .. todo:: sum of moments about a point must have an rxF term to get the
               same value as Patran.
    """
    F = force_in_global
    M = zeros(3, dtype=force_in_global.dtype)
    cp = coord_in
    coord_to = coord_out  # this is really cd
    r = cp.origin - coord_to.origin

    #cp_cid = cp.cid
    if cp.cid != 0:
        Fxyz_local_1 = cp.coord_to_xyz(F)
        Mxyz_local_1 = cp.coord_to_xyz(M)
        Fxyz_global = dot(Fxyz_local_1, cp.beta())
        Fxyz_local_2 = dot(dot(Fxyz_local_1, cp.beta()), coord_to.beta().T)
    else:
        Fxyz_local_1 = F
        Mxyz_local_1 = M
        Fxyz_global = Fxyz_local_1
        Fxyz_local_2 = dot(Fxyz_local_1, coord_to.beta().T)


    # find the moment about the new origin due to the force
    if npabs(r).max() > 0.:
        Mxyz_global = cross(r, Fxyz_global)
        dMxyz_local_2 = cross(r, Fxyz_local_2)
        Mxyz_local_2 = Mxyz_local_1 + dMxyz_local_2
    else:
        Mxyz_local_2 = r


    # rotate the delta moment into the local frame
    M_local = coord_to.xyz_to_coord(Mxyz_local_2)

    #return Fxyz_local_2, Mxyz_local_2
    return Fxyz_local_2


def transform_force_moment_from_global_to_local(force_in_global, moment_in_global,
                                                coord_in, coord_out, coords,
                                                xyz_cid0, offset):
    """
    Transforms force/moment from global to local.

    Supports cylindrical/spherical coordinate systems.

    Parameters
    ----------
    force_in_global : (N, 3) ndarray
        forces in the global frame
    moment_in_global : (N, 3) ndarray
        moments in the global frame
    coord_in : CORD()
        the global coordinate system object
    coord_out : CORD()
        the desired local frame
    xyz_cid0 : (n, 3) ndarray
        the nodes in the global frame
    offset : (3, ) ndarray
        the offset vector

    .. warning:: the function signature will change...
    .. todo:: doesn't support origin shifts
    .. todo:: sum of moments about a point must have an rxF term to get the
               same value as Patran.
    """
    delta = xyz_cid0 - offset

    F = force_in_global
    M = zeros(3, dtype=force_in_global.dtype)
    cp = coord_in
    coord_to = coord_out  # this is really cd
    r = cp.origin - coord_to.origin

    #cp_cid = cp.cid
    if cp.cid != 0:
        assert cp.Type == 'R', 'Only rectangular coordinate systems are supported'
        Fxyz_local_1 = cp.coord_to_xyz(F)
        Mxyz_local_1 = cp.coord_to_xyz(M)
        Fxyz_global = dot(Fxyz_local_1, cp.beta())
        Fxyz_local_2 = dot(dot(Fxyz_local_1, cp.beta()), coord_to.beta().T)
    else:
        Fxyz_local_1 = F
        Mxyz_local_1 = M
        Fxyz_global = Fxyz_local_1
        Fxyz_local_2 = dot(Fxyz_local_1, coord_to.beta().T)


    # find the moment about the new origin due to the force
    if npabs(r).max() > 0.:
        Mxyz_global = cross(r, Fxyz_global)
        dMxyz_local_2 = cross(r, Fxyz_local_2)
        Mxyz_local_2 = Mxyz_local_1 + dMxyz_local_2
    else:
        Mxyz_local_2 = r

    # rotate the delta moment into the local frame
    M_local = coord_to.xyz_to_coord(Mxyz_local_2)

    #delta = xyz_cid0 - offset
    #moment_offset = np.cross(delta, Fxyz_global) # is this right?


    #assert coord_to.Type == 'R', 'Only rectangular coordinate systems are supported'

    #return Fxyz_local_2, Mxyz_local_2
    return Fxyz_local_2, M_local # + moment_offset


def transform_force_moment_from_local_to_global(force_in_local, moment_in_local,
                                                coord_in, coord_out, coords,
                                                nid_cd, i_transform, beta_transforms): # , xyz_cid0, offset
    """
    Transforms force/moment from global to local.

    Supports cylindrical/spherical coordinate systems.

    Parameters
    ----------
    force_in_local : (N, 3) ndarray
        forces in the local frame
    moment_in_local : (N, 3) ndarray
        moments in the local frame
    coord_in : CORD()
        the global coordinate system object
    coord_out : CORD()
        the desired local frame
    #xyz_cid0 : (n, 3) ndarray
    #    the nodes in the global frame
    #offset : (3, ) ndarray
    #    the offset vector

    .. warning:: the function signature will change...
    .. todo:: sum of moments about a point must have an rxF term to get the
               same value as Patran.
    """
    force_in_global = zeros(force_in_local.shape, dtype=force_in_local.dtype)
    moment_in_global = zeros(moment_in_local.shape, dtype=moment_in_local.dtype)

    #assert coord_out.Type == 'R', 'Only rectangular coordinate systems are supported'

    nids = nid_cd[:, 0]
    cds = nid_cd[:, 1]
    ucds = unique(cds)
    for cd in ucds:
        if cd == 0:
            i = where(cds == cd)[0]
            nidsi = nids[i]
            #force_in_global[:, :3] = force_in_local[i, :3]
            #force_in_global[:, 3:] = force_in_local[i, 3:]
        else:
            raise NotImplementedError(ucds)
    return force_in_global, moment_in_global

if 0:
    delta = xyz_cid0 - offset
    moment_offset = np.cross(delta, Fxyz_global)

    F = force_in_global
    M = zeros(3, dtype=force_in_global.dtype)
    cp = coord_in
    coord_to = coord_out  # this is really cd
    r = cp.origin - coord_to.origin

    #cp_cid = cp.cid
    if cp.cid != 0:
        # local coord
        assert cp.Type == 'R', 'Only rectangular coordinate systems are supported'
        Fxyz_local_1 = cp.coord_to_xyz(F)
        Mxyz_local_1 = cp.coord_to_xyz(M)
        Fxyz_global = dot(Fxyz_local_1, cp.beta())
        Fxyz_local_2 = dot(dot(Fxyz_local_1, cp.beta()), coord_to.beta().T)
    else:
        # global coord
        Fxyz_local_1 = F
        Mxyz_local_1 = M
        Fxyz_global = Fxyz_local_1
        Fxyz_local_2 = dot(Fxyz_local_1, coord_to.beta().T)


    # find the moment about the new origin due to the force
    if npabs(r).max() > 0.:
        Mxyz_global = cross(r, Fxyz_global)
        dMxyz_local_2 = cross(r, Fxyz_local_2)
        Mxyz_local_2 = Mxyz_local_1 + dMxyz_local_2
    else:
        Mxyz_local_2 = r


    # rotate the delta moment into the local frame
    M_local = coord_to.xyz_to_coord(Mxyz_local_2)

    return Fxyz_local_2, Mxyz_local_2


def transform_force_moment(force_in_local, moment_in_local,
                           coord_out, coords,
                           nid_cd, i_transform, beta_transforms,
                           xyz_cid0, summation_point_cid0=None):
    """
    Transforms force/moment from global to local and returns all the forces.

    Supports cylindrical/spherical coordinate systems.

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

    .. warning:: the function signature will change...
    .. todo:: sum of moments about a point must have an rxF term to get the
               same value as Patran.

    Fglobal = Flocal @ T
    Flocal = T.T @ Fglobal
    Flocal2 = T2.T @ (Flocal1 @ T1)

    """
    force_out = zeros(force_in_local.shape, dtype=force_in_local.dtype)
    moment_out = zeros(force_in_local.shape, dtype=force_in_local.dtype)

    #assert coord_out.Type == 'R', 'Only rectangular coordinate systems are supported'

    nids = nid_cd[:, 0]
    cds = nid_cd[:, 1]
    ucds = unique(cds)

    coord_out_cid = coord_out.cid
    coord_out_T = coord_out.beta()
    #print(coord_out_T)

    #if summation_point_cid0 is None:
    summation_point_cid0 = np.array([0., 0., 0.])

    for cd in ucds:
        i = where(cds == cd)[0]
        #nidsi = nids[i]
        #summation_point_cid0 = np.array([10., 10., 10.])
        analysis_coord = coords[cd]
        #summation_point_cid0 = analysis_coord.origin
        #print('summation_point =', summation_point_cid0)
        print(analysis_coord)
        cd_T = analysis_coord.beta()
        #print(cd_T)

        # rxF from local_in to global to local_out
        force_in_locali = force_in_local[i, :]
        print('force_in_locali =\n', force_in_locali)
        moment_in_locali = moment_in_local[i, :]

        #if 0:
            #if analysis_coord.Type == 'R':
                #pass
            #elif analysis_coord.Type == 'C':
                ##R = force_in_locali[:, 0]
                ##theta = force_in_locali[:, 1]
                #force_in_locali = analysis_coord.coord_to_xyz_array(force_in_locali)
                #print(force_in_locali.shape)
            #else:
                #raise RuntimeError(analysis_coord)
        #force_in_locali = analysis_coord.coord_to_xyz_array(force_in_locali)

        force_in_globali = dot(force_in_locali, cd_T.T)
        moment_in_globali = dot(moment_in_locali, cd_T.T)
        force_outi = dot(coord_out_T, force_in_globali.T).T
        moment_outi = dot(coord_out_T, moment_in_globali.T).T

        # matrix multiplcation is NOT communitive
        #transformi = dot(coord_out_T, cd_T.T)
        #force_outi = dot(transformi, force_in_locali.T).T
        print('force_in_globali =\n', force_in_globali)
        print('force_outi =\n', force_outi)
        print('moment_outi =\n', moment_outi)

        force_out[i, :] = force_outi
        moment_out[i, :] = moment_outi
        if summation_point_cid0 is not None:
            print('xyz =', xyz_cid0[i, :])
            delta = xyz_cid0[i, :] - summation_point_cid0[np.newaxis, :]
            #delta_minus_origin = delta - analysis_coord.origin
            #delta_in_cid = dot(delta_minus_origin, cd_T.T)
            rxf = cross(delta, force_in_globali)
            rxf_in_cid = dot(coord_out_T, rxf.T).T
            print('rxf =\n', rxf)
            print('rxf_sum =', rxf.sum(axis=0), np.linalg.norm(rxf.sum(axis=0)))
            print('rxf_in_cid =', rxf_in_cid.sum(axis=0), np.linalg.norm(rxf_in_cid.sum(axis=0)))
            moment_out[i, :] += rxf_in_cid

    print('end')
    return -force_out, -moment_out

def transform_force_moment_sum(force_in_local, moment_in_local,
                               coord_out, coords,
                               nid_cd, i_transform, beta_transforms,
                               xyz_cid0, summation_point_cid0=None):
    """
    Transforms force/moment from global to local and returns a sum of forces/moments.

    Supports cylindrical/spherical coordinate systems.

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

    .. warning:: the function signature will change...
    .. todo:: sum of moments about a point must have an rxF term to get the
               same value as Patran.

    Fglobal = Flocal @ T
    Flocal = T.T @ Fglobal
    Flocal2 = T2.T @ (Flocal1 @ T1)

    """
    force_out = zeros(force_in_local.shape, dtype=force_in_local.dtype)
    moment_out = zeros(force_in_local.shape, dtype=force_in_local.dtype)
    force_sum = np.zeros(3)
    moment_sum = np.zeros(3)
    #assert coord_out.Type == 'R', 'Only rectangular coordinate systems are supported'

    nids = nid_cd[:, 0]
    cds = nid_cd[:, 1]
    ucds = unique(cds)

    coord_out_cid = coord_out.cid
    coord_out_T = coord_out.beta().T

    if summation_point_cid0 is None:
        summation_point_cid0 = coord_out.origin
        #summation_point_cid0 = np.array([0., 0., 0.])

    #summation_point_cid0 = np.array([10., 10., 10.])
    for cd in ucds:
        i = where(cds == cd)[0]
        #nidsi = nids[i]
        analysis_coord = coords[cd]
        cd_T = analysis_coord.beta()

        # rxF from local_in to global to local_out
        force_in_locali = force_in_local[i, :]
        moment_in_locali = moment_in_local[i, :]
        force_in_globali = dot(force_in_locali, cd_T)

        # matrix multiplcation is NOT communitive
        forcei = dot(coord_out_T, dot(force_in_globali, cd_T).T).T
        momenti = dot(coord_out_T, dot(moment_in_locali, cd_T).T).T
        force_out[i, :] = forcei
        moment_out[i, :] = momenti
        force_sum += forcei.sum(axis=0)
        moment_sum += momenti.sum(axis=0)

        if summation_point_cid0 is not None:
            delta = xyz_cid0[i, :] - summation_point_cid0[np.newaxis, :]
            rxf = cross(delta, force_in_globali)
            moment_out[i, :] += rxf
            moment_sum += rxf.sum(axis=0)

    return force_out, moment_out, force_sum, moment_sum


if __name__ == '__main__':  # pragma: no cover
    #print(iformat('4si3f'))
    test_abs_max_min_global()
    test_abs_max_min_vector()
