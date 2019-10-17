"""
defines some methods for working with arrays:
 - filter1d(a, b=None, zero_tol=0.001)
 - where_searchsorted(a, v, side='left', x=None, y=None)
 - sortedsum1d(ids, values, axis=None)
 - iformat(format_old, precision=2)
 - abs_max_min_global(values)
 - abs_max_min_vector(values)
 - abs_max_min(values, global_abs_max=True)
 - principal_3d(o11, o22, o33, o12, o23, o13)
 - transform_force(force_in_local,
                   coord_out, coords,
                   nid_cd, i_transform)
 - transform_force_moment(force_in_local, moment_in_local,
                          coord_out, coords,
                          nid_cd, i_transform,
                          xyz_cid0, summation_point_cid0=None,
                          consider_rxf=True,
                          debug=False, log=None)
 - transform_force_moment_sum(force_in_local, moment_in_local,
                              coord_out, coords,
                              nid_cd, i_transform,
                              xyz_cid0, summation_point_cid0=None,
                              consider_rxf=True,
                              debug=False, log=None)

"""
from __future__ import annotations
from struct import calcsize
from itertools import count
from typing import Optional, Dict, TYPE_CHECKING
import numpy as np
from numpy import arccos, sqrt, pi, in1d, cos, unique, cross, ndarray
if TYPE_CHECKING:  # pragma: no cover
    from pyNastran.bdf.bdf import CORDx # , CORD1R, CORD1C, CORD1S, CORD2R, CORD2C, CORD2S


def filter1d(a: ndarray, b: Optional[ndarray]=None, zero_tol: float=0.001):
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
    """Implements a np.where that assumes a sorted array set."""
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

    Examples
    --------
    **1D Example**
    For an set  of arrays, there are 5 values in sorted order.

    ..code-block :: python

       ids = [1, 1, 2, 2, 3, 3, 3, 4, 5, 5, 5]
       values = 1.0 * ids

    We want to sum the values such that:
    ..code-block :: python
        out = [2.0, 4.0, 9.0, 4.0, 15.0]

    **2D Example**
    ..code-block :: python

       ids = [1, 1, 2, 2, 3, 3, 3, 4, 5, 5, 5]
       values = [
           [1.0, 1.0, 2.0, 2.0, 3.0, 3.0],
           [1.0, 1.0, 2.0, 2.0, 3.0, 3.0],
       ]


    values = 1.0 * ids
    For 2D

    .. todo::  This could probably be more efficient
    .. todo::  Doesn't support axis

    """
    uids = unique(ids)
    i1 = np.searchsorted(ids, uids, side='left') # left is the default
    i2 = np.searchsorted(ids, uids, side='right')
    out = np.zeros(values.shape, dtype=values.dtype)

    for i, i1i, i2i in zip(count(), i1, i2):
        out[i, :] = values[i1i:i2i, :].sum(axis=axis)
    return out

def iformat(format_old: str, precision: int=2) -> str:
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

    Examples
    --------
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
    """
    Gets the maximum value of x and -x.
    This is used for getting the max/min principal stress.
    """
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


def transform_force(force_in_local,
                    coord_out: CORDx, coords: Dict[int, CORDx],
                    nid_cd: int, unused_icd_transform):
    """
    Transforms force/moment from global to local and returns all the forces.

    Supports cylindrical/spherical coordinate systems.

    Parameters
    ----------
    force : (N, 3) ndarray
        forces in the local frame
    coord_out : CORD()
        the desired local frame
    coords : dict[int] = CORDx
        all the coordinate systems
        key : int
        value : CORDx
    nid_cd : (M, 2) int ndarray
        the (BDF.point_ids, cd) array
    icd_transform : dict[cd] = (Mi, ) int ndarray
        the mapping for nid_cd

    .. warning:: the function signature will change...
    .. todo:: sum of moments about a point must have an rxF term to get the
               same value as Patran.

    Fglobal = Flocal @ T
    Flocal = T.T @ Fglobal
    Flocal2 = T2.T @ (Flocal1 @ T1)

    """
    force_out = np.zeros(force_in_local.shape, dtype=force_in_local.dtype)

    #nids = nid_cd[:, 0]
    cds = nid_cd[:, 1]
    ucds = unique(cds)

    coord_out_cid = coord_out.cid
    coord_out_T = coord_out.beta()

    for cd in ucds:
        i = np.where(cds == cd)[0]
        #nidsi = nids[i]
        analysis_coord = coords[cd]
        cd_T = analysis_coord.beta()

        # rxF from local_in to global to local_out
        force_in_locali = force_in_local[i, :]

        force_in_globali = force_in_locali @ cd_T
        force_outi = (coord_out_T @ force_in_globali.T).T
        force_out[i, :] = force_outi
    return -force_out


def transform_force_moment(force_in_local, moment_in_local,
                           coord_out: CORDx, coords: Dict[int, CORDx],
                           nid_cd: int, icd_transform: Dict[int, ndarray],
                           xyz_cid0: ndarray,
                           summation_point_cid0: Optional[ndarray]=None,
                           consider_rxf: bool=True,
                           debug: bool=False, log=None):
    """
    Transforms force/moment from global to local and returns all the forces.

    Parameters
    ----------
    force_in_local : (N, 3) ndarray
        forces in the local frame
    moment_in_local : (N, 3) ndarray
        moments in the local frame
    coord_out : CORDx()
        the desired local frame
    coords : dict[int] = CORDx
        all the coordinate systems
        key : int
        value : CORDx
    nid_cd : (M, 2) int ndarray
        the (BDF.point_ids, cd) array
    icd_transform : dict[cd] = (Mi, ) int ndarray
        the mapping for nid_cd
    xyz_cid0 : (n, 3) ndarray
        the nodes in the global frame
    summation_point_cid0 : (3, ) ndarray
        the summation point in the global frame???
    consider_rxf : bool; default=True
        considers the r x F term
    debug : bool; default=False
        debugging flag
    log : log; default=None
        a log object that gets used when debug=True

    Returns
    -------
    force_out : (n, 3) float ndarray
        the ith float components in the coord_out coordinate frame
    moment_out : (n, 3) float ndarray
        the ith moment components about the summation point in the
        coord_out coordinate frame

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
    #print('consider_rxf =', consider_rxf)
    #debug = True
    assert log is not None
    assert nid_cd.shape[0] == force_in_local.shape[0]
    dtype = force_in_local.dtype
    #dtype = 'float64'

    force_in_local_sum = force_in_local.sum(axis=0)
    force_out = np.zeros(force_in_local.shape, dtype=dtype)
    moment_out = np.zeros(force_in_local.shape, dtype=dtype)

    nids = nid_cd[:, 0]
    cds = nid_cd[:, 1] #* 0
    ucds = unique(cds)

    #coord_out_cid = coord_out.cid
    beta_out = coord_out.beta().T

    if debug:
        log.debug('beta_out =\n%s' % beta_out)
        log.debug(coord_out)
        if consider_rxf:
            for ii in range(xyz_cid0.shape[0]):
                log.debug('***i=%s xyz=%s nid=%s cd=%s' % (
                    ii, xyz_cid0[ii, :], nid_cd[ii, 0], nid_cd[ii, 1]))
        log.debug('------------')
        log.debug('ucds = %s' % ucds)

    if consider_rxf and summation_point_cid0 is None:
        summation_point_cid0 = np.array([0., 0., 0.])

    #eye = np.eye(3, dtype=beta_cd.dtype)
    for cd in ucds:
        #log.debug('cd = %s' % cd)
        i = np.where(cds == cd)[0]
        nidsi = nids[i]
        analysis_coord = coords[cd]
        beta_cd = analysis_coord.beta()

        force_in_locali = -force_in_local[i, :]
        moment_in_locali = -moment_in_local[i, :]

        #if 0 and np.array_equal(beta_cd, eye):
        force_in_globali = force_in_locali
        moment_in_globali = moment_in_locali
        if debug:
            #log.debug('analysis_coord =\n%s' % analysis_coord)
            log.debug('beta_cd =\n%s' % beta_cd)
            log.debug('i = %s' % i)
            log.debug('force_in_local = %s' % force_in_local)
            log.debug('force_in_local.sum() = %s' % force_in_local_sum)
        #force_in_locali.astype('float64')
        #moment_in_locali.astype('float64')

        # rotate loads from an arbitrary coordinate system to local xyz
        if 0:
            log.debug(analysis_coord)
            log.debug(force_in_locali)
            force_in_locali = analysis_coord.coord_to_xyz_array(force_in_locali)
            moment_in_locali = analysis_coord.coord_to_xyz_array(moment_in_locali)

            if debug:
                log.debug('i = %s' % i)
                log.debug('nids = %s' % nidsi)
                log.debug('force_input = %s' % force_in_locali)

        # rotate the forces/moments into a coordinate system coincident
        # with the local frame and with the same primary directions
        # as the global frame
        force_in_globali = force_in_locali.dot(beta_cd)
        moment_in_globali = moment_in_locali.dot(beta_cd)

        # rotate the forces and moments into a coordinate system coincident
        # with the output frame and with the same primary directions
        # as the output frame

        #if 0 and np.array_equal(beta_out, eye):
            #force_outi = force_in_globali
            #moment_outi = moment_in_globali
        force_outi = force_in_globali.dot(beta_out)
        moment_outi = moment_in_globali.dot(beta_out)

        if debug:
            #if show_local:
            log.debug('force_in_locali = \n%s' % force_in_locali.T)
            log.debug('force_in_globali = \n%s' % force_in_globali.T)
            log.debug('force_in_locali.sum()  = %s' % force_in_locali.sum(axis=0))
            log.debug('force_in_globali.sum() = %s' % force_in_globali.sum(axis=0))
            log.debug('force_outi = %s' % force_outi)
            #log.debug('moment_outi = %s' % moment_outi)

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

            rxf_in_cid = rxf.dot(beta_out)
            if debug:
                log.debug('delta_moment = %s' % delta)
                #log.debug('rxf = %s' % rxf.T)
                #log.debug('rxf_in_cid = %s' % rxf_in_cid)

            moment_out[i, :] += rxf_in_cid

    # rotate loads from the local XYZ to an arbitrary coordinate system
    # flip the sign of the output to be consistent with Patran
    #force_out2 = -coord_out.xyz_to_coord_array(force_out)
    #moment_out2 = -coord_out.xyz_to_coord_array(moment_out)
    #return force_out2, moment_out2
    return force_out, moment_out


def transform_force_moment_sum(force_in_local, moment_in_local,
                               coord_out, coords,
                               nid_cd, icd_transform,
                               xyz_cid0, summation_point_cid0=None,
                               consider_rxf=True,
                               debug=False, log=None):
    """
    Transforms force/moment from global to local and returns a sum of forces/moments.

    Parameters
    ----------
    force_in_local : (N, 3) ndarray
        forces in the local frame
    moment_in_local : (N, 3) ndarray
        moments in the local frame
    coord_out : CORDx()
        the desired local frame
    coords : dict[int] = CORDx
        all the coordinate systems
        key : int
        value : CORDx
    nid_cd : (M, 2) int ndarray
        the (BDF.point_ids, cd) array
    icd_transform : dict[cd] = (Mi, ) int ndarray
        the mapping for nid_cd
    xyz_cid0 : (nnodes + nspoints + nepoints, 3) ndarray
        the grid locations in coordinate system 0
    summation_point_cid0 : (3, ) ndarray
        the summation point in the global frame
    consider_rxf : bool; default=True
        considers the r x F term
    debug : bool; default=False
        debugging flag
    log : log; default=None
        a log object that gets used when debug=True

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
    assert log is not None
    out = transform_force_moment(
        force_in_local, moment_in_local,
        coord_out, coords, nid_cd,
        icd_transform, xyz_cid0,
        summation_point_cid0=summation_point_cid0, consider_rxf=consider_rxf,
        debug=debug, log=log)
    force_out, moment_out = out
    if debug:
        log.debug('force_sum = %s' % force_out.sum(axis=0))
        if consider_rxf:
            log.debug('moment_sum = %s' % moment_out.sum(axis=0))
    return force_out, moment_out, force_out.sum(axis=0), moment_out.sum(axis=0)
