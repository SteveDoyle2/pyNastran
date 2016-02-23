from __future__ import print_function
from struct import calcsize
from itertools import count

import numpy as np
from numpy import array, where, zeros, asarray, dot, arccos, sqrt, pi
from numpy import cos, unique, cross, abs as npabs

def sortedsum1d(ids, values, axis=None):
    """
    Sums the values in a sorted 1d/2d array.

    Parameters
    ----------
    ids : (n, ) int ndarray
        the integer values in a sorted order that indicate what is to
        be summed
    values : (n, m) float ndarray
        the values to be summed
        Presumably m is 2 because this function is intended to be used
        for summing forces and moments
    axis : None or int or tuple of ints, optional
        Axis or axes along which a sum is performed. The default
        (axis=None) is perform a sum over all the dimensions of the
        input array. axis may be negative, in which case it counts
        from the last to the first axis.
        If this is a tuple of ints, a sum is performed on multiple
        axes, instead of a single axis or all the axes as before.

        Not actually supported...

    Returns
    -------
    out : (nunique, m) float ndarray
        the summed values

    1D Example
    ----------
    For an set  of arrays, there are 5 values in sorted order.

    ..code-block :: python

       ids = [1, 1, 2, 2, 3, 3, 3, 4, 5, 5, 5]
       values = 1.0 * ids

    We want to sum the values such that:
    ..code-block :: python
        out = [2.0, 4.0, 9.0, 4.0, 15.0]

    2D Example
    ----------
    ..code-block :: python

       ids = [1, 1, 2, 2, 3, 3, 3, 4, 5, 5, 5]
       values = [
           [1.0, 1.0, 2.0, 2.0, 3.0, 3.0],
           [1.0, 1.0, 2.0, 2.0, 3.0, 3.0],
       ]


    values = 1.0 * ids
    For 2D
    Todo
    ----
    This could probably be more efficient
    Doesn't support axis
    """
    uids = np.unique(ids)
    i1 = np.searchsorted(ids, uids, side='left') # left is the default
    i2 = np.searchsorted(ids, uids, side='right')
    print(i1, i2)
    out = zeros(values.shape, dtype=values.dtype)

    for i, i1i, i2i in zip(count(), i1, i2):
        print(i, i1i, i2i)
        out[i, :] = values[i1i:i2i, :].sum(axis=axis)
    return out

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


def transform_force(force_in_local,
                    coord_out, coords,
                    nid_cd, i_transform, beta_transforms):
    """
    Transforms force/moment from global to local and returns all the forces.

    Supports cylindrical/spherical coordinate systems.

    Parameters
    ----------
    force : (N, 3) ndarray
        forces in the local frame
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

    nids = nid_cd[:, 0]
    cds = nid_cd[:, 1]
    ucds = unique(cds)

    coord_out_cid = coord_out.cid
    coord_out_T = coord_out.beta()

    for cd in ucds:
        i = where(cds == cd)[0]
        nidsi = nids[i]
        analysis_coord = coords[cd]
        cd_T = analysis_coord.beta()

        # rxF from local_in to global to local_out
        force_in_locali = force_in_local[i, :]

        force_in_globali = dot(force_in_locali, cd_T)
        force_outi = dot(coord_out_T, force_in_globali.T).T
        force_out[i, :] = force_outi
    return -force_out


def transform_force_moment(force_in_local, moment_in_local,
                           coord_out, coords,
                           nid_cd, i_transform, beta_transforms,
                           xyz_cid0, summation_point_cid0=None,
                           consider_rxf=True,
                           debug=False, logger=None):
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
    logger : logger; default=None
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
    """
    debug = True
    assert logger is not None
    #print('nid_cd*\n', nid_cd)
    dtype = force_in_local.dtype
    dtype = 'float64'
    force_out = zeros(force_in_local.shape, dtype=dtype)
    moment_out = zeros(force_in_local.shape, dtype=dtype)

    #assert coord_out.Type == 'R', 'Only rectangular coordinate systems are supported'

    nids = nid_cd[:, 0]
    cds = nid_cd[:, 1]
    ucds = unique(cds)

    coord_out_cid = coord_out.cid
    coord_out_T = coord_out.beta()

    if debug:
        logger.debug(coord_out)
        if consider_rxf:
            for ii in range(xyz_cid0.shape[0]):
                logger.debug('***i=%s xyz=%s nid=%s cd=%s' % (
                    ii, xyz_cid0[ii, :], nid_cd[ii, 0], nid_cd[ii, 1]))
        logger.debug('------------')

    if consider_rxf and summation_point_cid0 is None:
        summation_point_cid0 = np.array([0., 0., 0.])

    for cd in ucds:
        i = where(cds == cd)[0]
        nidsi = nids[i]
        analysis_coord = coords[cd]
        if debug:
            logger.debug('i = %s' % i)
            logger.debug('nids = %s' % nidsi)
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
            logger.debug('force_input =%s' % force_in_locali)

        if 1:
            if analysis_coord.Type == 'R':
                pass
            elif analysis_coord.Type == 'C':
                #R = force_in_locali[:, 0]
                #theta = force_in_locali[:, 1]
                #print(analysis_coord.origin)
                force_in_locali = analysis_coord.coord_to_xyz_array(force_in_locali)
                moment_in_locali = analysis_coord.coord_to_xyz_array(moment_in_locali)
                logger.debug('cylindrical')
            else:
                raise RuntimeError(analysis_coord)

        #print('coord_out_T\n', coord_out_T)

        # is this:
        #  1.  cd_T.T and coord_out_T
        #  2.  cd_T and coord_out_T.T
        #  3.  reverse it all...
        if 0:
            force_in_globali = dot(cd_T, force_in_locali.T)
            moment_in_globali = dot(cd_T, moment_in_locali.T)
            force_outi = dot(coord_out_T.T, force_in_globali).T
            moment_outi = dot(coord_out_T.T, moment_in_globali).T
        elif 0:
            force_in_globali = dot(cd_T.T, force_in_locali.T)
            moment_in_globali = dot(cd_T.T, moment_in_locali.T)
            force_outi = dot(coord_out_T, force_in_globali).T
            moment_outi = dot(coord_out_T, moment_in_globali).T
        else:
            # old
            force_in_globali = dot(force_in_locali, cd_T).T
            moment_in_globali = dot(moment_in_locali, cd_T).T
            #force_outi = dot(coord_out_T, force_in_globali).T
            #moment_outi = dot(coord_out_T, moment_in_globali).T
            force_outi = dot(coord_out_T, force_in_globali).T
            moment_outi = dot(coord_out_T, moment_in_globali).T

        if debug:
            logger.debug('force_in_locali =%s' % force_in_locali.T)
            logger.debug('force_in_globali = %s' % force_in_globali.T)
            logger.debug('force_outi = %s' % force_outi)
            #logger.debug('moment_outi = %s' % moment_outi)

        force_out[i, :] = force_outi
        moment_out[i, :] = moment_outi

        # rxF from local_in to global to local_out
        if consider_rxf:
            delta = xyz_cid0[i, :] - summation_point_cid0[np.newaxis, :]
            #delta_minus_origin = delta - analysis_coord.origin
            #delta_in_cid = dot(delta_minus_origin, cd_T.T)
            if debug:
                logger.debug('delta = %s' % delta)
            rxf = cross(delta, force_in_globali.T)

            if debug:
                logger.debug('rxf = %s' % rxf)
            rxf_in_cid = dot(coord_out_T, rxf.T).T
            if debug:
                logger.debug('rxf_in_cid = %s' % rxf_in_cid)

            moment_out[i, :] += rxf_in_cid
        #print('------------')
    force_out = coord_out.xyz_to_coord_array(force_out)
    moment_out = coord_out.xyz_to_coord_array(moment_out)
    return -force_out, -moment_out


def transform_force_moment_sum(force_in_local, moment_in_local,
                               coord_out, coords,
                               nid_cd, i_transform, beta_transforms,
                               xyz_cid0, summation_point_cid0=None,
                               consider_rxf=True,
                               debug=False, logger=None):
    """
    Transforms force/moment from global to local and returns a sum of forces/moments.

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
    logger : logger; default=None
        a logger object that gets used when debug=True

    Returns
    -------
    force_out : (n, 3) float ndarray
        the ith float components in the coord_out coordinate frame
    moment_out : (n, 3) float ndarray
        the ith moment components about the summation point in the coord_out coordinate frame
    force_out_sum : (3, ) float ndarray
        the sum of forces in the coord_out coordinate frame
    moment_out_sum : (3, ) float ndarray
        the sum of moments about the summation point in the coord_out coordinate frame

    .. todo:: doesn't seem to handle cylindrical/spherical systems
    """
    assert logger is not None
    out = transform_force_moment(
        force_in_local, moment_in_local,
        coord_out, coords,  nid_cd,
        i_transform, beta_transforms, xyz_cid0,
        summation_point_cid0=summation_point_cid0, consider_rxf=consider_rxf,
        debug=debug, logger=logger)
    force_out, moment_out = out
    return force_out, moment_out, force_out.sum(axis=0), moment_out.sum(axis=0)


if __name__ == '__main__':  # pragma: no cover
    #print(iformat('4si3f'))
    test_abs_max_min_global()
    test_abs_max_min_vector()
