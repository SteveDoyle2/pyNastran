from __future__ import print_function
from struct import calcsize
from itertools import count

import numpy as np
from numpy import array, where, zeros, asarray, dot, arccos, sqrt, pi
from numpy import cos, unique, cross, abs as npabs

def filter1d(a, b=None, zero_tol=0.001):
    """

    Filters a 1d numpy array of values near 0.

    Parameters
    ----------
    a : (n, ) float ndarray
        a vector to compare
    b : (n, ) float ndarray; default=None
        another vector to compare
        If b is defined, both a and b must be near 0 for a value to be removed.
    zero_tol : float; default=0.001
        the zero tolerance value

    Returns
    -------
    k : (m, ) int ndarray
        the indices of the removed values

    a = [1., 2., 0.1]
    >>> i = filter(a, zero_tol=0.5)
    a[i]
    >>> [0, 1]

    a = [1., 2., 0.1]
    b = [1., -0.1, 0.1]
    >>> i = filter(a, b, zero_tol=0.5)
    [0, 1]
    """
    a = np.asarray(a)
    i = np.where(np.abs(a) > zero_tol)[0]
    if b is None:
        return i
    b = np.asarray(b)
    assert a.shape == b.shape, 'a.shape=%s b.shape=%s' % (str(a.shape), str(b.shape))
    assert a.size == b.size, 'a.size=%s b.size=%s' % (str(a.size), str(b.size))

    j = np.where(np.abs(b) > zero_tol)[0]
    k = np.unique(np.hstack([i, j]))
    return k

def where_searchsorted(a, v, side='left', x=None, y=None):
    """ where(a, b, [x, y])

        Return elements, either from `x` or `y`, depending on `condition`.

        If only `condition` is given, return ``condition.nonzero()``.

        Parameters
        ----------
        condition : array_like, bool
            When True, yield `x`, otherwise yield `y`.
        x, y : array_like, optional
            Values from which to choose. `x` and `y` need to have the same
            shape as `condition`.

        Returns
        -------
        out : ndarray or tuple of ndarrays
            If both `x` and `y` are specified, the output array contains
            elements of `x` where `condition` is True, and elements from
            `y` elsewhere.

            If only `condition` is given, return the tuple
            ``condition.nonzero()``, the indices where `condition` is True.

        See Also
        --------
        nonzero, choose

        Notes
        -----
        If `x` and `y` are given and input arrays are 1-D, `where` is
        equivalent to::

            [xv if c else yv for (c,xv,yv) in zip(condition,x,y)]

        Examples
        --------
        >>> np.where([[True, False], [True, True]],
        ...          [[1, 2], [3, 4]],
        ...          [[9, 8], [7, 6]])
        array([[1, 8],
               [3, 4]])

        >>> np.where([[0, 1], [1, 0]])
        (array([0, 1]), array([1, 0]))

        >>> x = np.arange(9.).reshape(3, 3)
        >>> np.where( x > 5 )
        (array([2, 2, 2]), array([0, 1, 2]))
        >>> x[np.where( x > 3.0 )]               # Note: result is 1D.
        array([ 4.,  5.,  6.,  7.,  8.])
        >>> np.where(x < 5, x, -1)               # Note: broadcasting.
        array([[ 0.,  1.,  2.],
               [ 3.,  4., -1.],
               [-1., -1., -1.]])

        Find the indices of elements of `x` that are in `goodvalues`.

        >>> goodvalues = [3, 4, 7]
        >>> ix = np.in1d(x.ravel(), goodvalues).reshape(x.shape)
        >>> ix
        array([[False, False, False],
               [ True,  True, False],
               [False,  True, False]], dtype=bool)
        >>> np.where(ix)
        (array([1, 1, 2]), array([0, 1, 1])) """
    # TODO: take advantage of searchsorted
    assert x is None, x
    assert y is None, y
    assert side == 'left', side
    i = np.where(in1d(a, v), x=x, y=y)
    return i

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
    uids = unique(ids)
    i1 = np.searchsorted(ids, uids, side='left') # left is the default
    i2 = np.searchsorted(ids, uids, side='right')
    print(i1, i2)
    out = np.zeros(values.shape, dtype=values.dtype)

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
    values = np.asarray(values)

    # find the [max,
    #           min]
    # we organize it as [max, min], which is why the note applies
    # we could make both of the edge cases return -3.0, but if you're using
    # this function it shouldn't matter
    values2 = np.array([values.max(), values.min()])

    # we figure out the absolute max/min
    abs_vals = np.abs(values2)
    abs_val = abs_vals.max()

    # find the location of the absolute max value
    # 1.  we take the first value (the where[0]) to chop the return value
    #     since there is no else conditional
    # 2.  we take the first value (the where[0][0]) to only get the max
    #     value if 2+ values are returned
    j = np.where(abs_val == abs_vals)[0][0]

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
    values = np.asarray(values)

    # find the [maxs,
    #           mins]
    # we organize it as [maxs, mins], which is why the note applies
    # we could make both of the edge cases return -3.0, but if you're using
    # this function it shouldn't matter
    maxs_mins = np.array([values.max(axis=1), values.min(axis=1)])

    # we figure out the absolute max/min for each row
    abs_vals = np.abs(maxs_mins)
    absolute_maxs = abs_vals.max(axis=0)

    outs = np.zeros(absolute_maxs.shape[0], dtype=values.dtype)
    for i, absolute_max in enumerate(absolute_maxs):
        # find the location of the absolute max value
        # 1.  we take the first value (the where[0]) to chop the return value
        #     since there is no else conditional
        # 2.  we take the first value (the where[0][0]) to only get the max
        #     value if 2+ values are returned
        j = np.where(absolute_max == abs_vals[:, i])[0][0]

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
    max_min_mid = np.array([p1, p2, p3])
    pmax = max_min_mid.max(axis=0)
    pmin = max_min_mid.min(axis=0)
    return pmax, pmin


def test_abs_max_min_global():
    #print(iformat('4si3f', 2))
    print(abs_max_min_global([0.0, 2.0, 1.0]))
    print(abs_max_min_global([0.0, 2.0, -1.0]))
    print(abs_max_min_global([0.0, 2.0, -3.0]))
    print(abs_max_min_global(np.array([0.0, 2.0, -3.0])))
    print(abs_max_min_global([1.0]))

    # gets the global max/min value
    print(abs_max_min_global([
        [0.0, 2.0, -3.0],
        [0.0, 2.0, -4.0],
    ]))
    print(abs_max_min_global(np.array([
        [0.0, 2.0, -3.0],
        [0.0, 2.0, -4.0],
    ])))


def test_abs_max_min_vector():
    print(abs_max_min_vector(np.array([
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
    print(abs_max_min_vector(np.array([
        [0.0, 2.0, 1.0],
        [0.0, 2.0, -1.0],
        [0.0, 2.0, -3.0],
        [0.0, 2.0, 4.0],
    ])))

    print(abs_max_min_vector(np.array([
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
    force_out = np.zeros(force_in_local.shape, dtype=force_in_local.dtype)

    nids = nid_cd[:, 0]
    cds = nid_cd[:, 1]
    ucds = unique(cds)

    coord_out_cid = coord_out.cid
    coord_out_T = coord_out.beta()

    for cd in ucds:
        i = np.where(cds == cd)[0]
        nidsi = nids[i]
        analysis_coord = coords[cd]
        cd_T = analysis_coord.beta()

        # rxF from local_in to global to local_out
        force_in_locali = force_in_local[i, :]

        force_in_globali = np.dot(force_in_locali, cd_T)
        force_outi = np.dot(coord_out_T, force_in_globali.T).T
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

    Method
    ------
    xyz_g = T_a2g @ xyz_a
    xyz_g = T_b2g @ xyz_b
    T_b2g @ xyz_b = T_a2g @ xyz_a
    xyz_b = T_b2g.T @ T_a2g @ xyz_a = T_g2b @ T_a2g @ xyz_a
    """
    #debug = True
    assert logger is not None
    dtype = force_in_local.dtype
    #dtype = 'float64'
    force_out = zeros(force_in_local.shape, dtype=dtype)
    moment_out = zeros(force_in_local.shape, dtype=dtype)

    nids = nid_cd[:, 0]
    cds = nid_cd[:, 1]
    ucds = unique(cds)

    #coord_out_cid = coord_out.cid
    beta_out = coord_out.beta().T
    print('beta_out =\n%s' % beta_out)

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
        beta_cd = analysis_coord.beta()
        print('beta_cd =\n%s' % beta_cd)

        print('i =', i)
        print('force_in_local =', force_in_local)
        #force_in_locali.astype('float64')
        #moment_in_locali.astype('float64')

        # rotate loads from an arbitrary coordinate system to local xyz
        if 0:
            force_in_locali = -force_in_local[i, :]
            moment_in_locali = -moment_in_local[i, :]
        else:
            force_in_locali = -force_in_local[i, :]
            moment_in_locali = -moment_in_local[i, :]
            print(analysis_coord)
            print(force_in_locali)
            force_in_locali = analysis_coord.coord_to_xyz_array(force_in_locali)
            moment_in_locali = analysis_coord.coord_to_xyz_array(moment_in_locali)


            if debug:
                logger.debug('i = %s' % i)
                logger.debug('nids = %s' % nidsi)
                logger.debug('force_input =%s' % force_in_locali)

        # rotate the forces/moments into a coordinate system coincident
        # with the local frame and with the same primary directions
        # as the global frame
        force_in_globali = dot(force_in_locali, beta_cd)
        print('force_in_globali = %s' % force_in_globali)
        moment_in_globali = dot(moment_in_locali, beta_cd)

        # rotate the forces and moments into a coordinate system coincident
        # with the output frame and with the same primary directions
        # as the output frame
        force_outi = dot(force_in_globali, beta_out)
        moment_outi = dot(moment_in_globali, beta_out)

        if debug:
            logger.debug('force_in_locali =%s' % force_in_locali.T)
            logger.debug('force_in_globali = %s' % force_in_globali.T)
            logger.debug('force_outi = %s' % force_outi)
            logger.debug('moment_outi = %s' % moment_outi)

        # these are in the local XYZ coordinates
        # we'll do the final transform later
        force_out[i, :] = force_outi
        moment_out[i, :] = moment_outi

        # Now we need to consider the "r x F" term.
        # We calculate it in the global xyz frame about the summation
        # point.  Then we transform it to the output XYZ frame
        if consider_rxf:
            delta = xyz_cid0[i, :] - summation_point_cid0[np.newaxis, :]
            rxf = cross(delta, force_in_globali)

            rxf_in_cid = dot(rxf, beta_out)
            if debug:
                logger.debug('delta = %s' % delta)
                logger.debug('rxf = %s' % rxf.T)
                logger.debug('rxf_in_cid = %s' % rxf_in_cid)

            moment_out[i, :] += rxf_in_cid

    # rotate loads from the local XYZ to an arbitrary coordinate system
    # flip the sign of the output to be consistent with Patran
    #force_out2 = -coord_out.xyz_to_coord_array(force_out)
    #moment_out2 = -coord_out.xyz_to_coord_array(moment_out)
    #return force_out2, moment_out2
    return force_out, moment_out


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
        coord_out, coords, nid_cd,
        i_transform, beta_transforms, xyz_cid0,
        summation_point_cid0=summation_point_cid0, consider_rxf=consider_rxf,
        debug=debug, logger=logger)
    force_out, moment_out = out
    if debug:
        logger.debug('force_sum = %s' % force_out.sum(axis=0))
        logger.debug('moment_sum = %s' % moment_out.sum(axis=0))
    return force_out, moment_out, force_out.sum(axis=0), moment_out.sum(axis=0)


if __name__ == '__main__':  # pragma: no cover
    #print(iformat('4si3f'))
    test_abs_max_min_global()
    test_abs_max_min_vector()
