from struct import calcsize
from numpy import array, abs, where, zeros, asarray

def iformat(Format, precision=2):
    if precision == 2:  # double
        f = Format.replace('i', 'l').replace('f', 'd')
    elif precision == 1: # single
        f = Format.replace('l', 'i').replace('d', 'f')
    else:
        raise NotImplementedError(precision)
    ndata = calcsize(f)
    return f, ndata

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
    abs_vals = abs(values2)
    abs_val = abs_vals.max()

    # find the location of the absolute max value
    # 1.  we take the first value (the where[0]) to chop the return value
    #     since there is no else conditional
    # 2.  we take the first value (the where[0][0]) to only get the max
    #     value if 2+ values are returned
    j = where(abs_val == abs_vals)[0][0]

    # get the raw value from the absoluted value, so:
    # value = abs(raw_value)
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
    abs_vals = abs(maxs_mins)
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
        # value = abs(raw_value)
        outs[i] = maxs_mins[j, i]
    return outs


def abs_max_min(values, global_abs_max=True):
    if global_abs_max:
        return abs_max_min_global(values)
    return abs_max_min_vector(values)


def principal_2d(o11, o22, o12):
    oxx = 5
    return oxx, oy
def principal_3d(o11, o22, o33, o12, o23, o13):
    """http://www.continuummechanics.org/cm/principalstrain.html"""
    e = a
    i1 = o11 + o22 + o33
    i2 = o11*o22 + o22*o33 + o11*o33 - o12**2 - o13**2 - o23**2
    i3 = o11*o22*o33 - o111*o23**2 - o22*o13**2 + 2*o12*o13*o23
    Q = 3 * i2 - i1**2
    R = (2*i1**3 - 9*i1*i2 + 27*i3) / 54.
    theta = arccos(R / sqrt(-Q**3))

    q2 = 2 * sqrt(-Q)
    i13 = 1./3. * i1
    p1 = q2 * cos(theta/3) + i13
    p2 = q2 * cos(theta/3 + 2*pi/3.) + i13
    p3 = q2 * cos(theta/3 + 4*pi/3.) + i13
    max_min_mid = array([p1, p2, p3])
    pmax = max_mid_mid.max(axis=0)
    pmax = max_mid_mid.min(axis=0)
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
    force_in_global = zeros(force_in_local.shape, dtype='float32')
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
    F = force_in_global
    M = zeros(3, dtype='float32')
    cp = coord_in
    coord_to = coord_out  # this is really cd
    r = cp.origin - coord_to.origin

    #cp_cid = cp.cid
    if cp.cid != 0:
        Fxyz_local_1 = cp.coordToXYZ(F)
        Mxyz_local_1 = cp.coordToXYZ(M)
        Fxyz_global = dot(Fxyz_local_1, cp.beta())
        Fxyz_local_2 = dot(dot(Fxyz_local_1, cp.beta()), coord_to.beta().T)
    else:
        Fxyz_local_1 = F
        Mxyz_local_1 = M
        Fxyz_global = Fxyz_local_1
        Fxyz_local_2 = dot(Fxyz_local_1, coord_to.beta().T)


    # find the moment about the new origin due to the force
    if abs(r).max() > 0.:
        Mxyz_global = cross(r, Fxyz_global)
        dMxyz_local_2 = cross(r, Fxyz_local_2)
        Mxyz_local_2 = Mxyz_local_1 + dMxyz_local_2
    else:
        Mxyz_local_2 = r


    # rotate the delta moment into the local frame
    M_local = coord_to.XYZtoCoord(Mxyz_local_2)

    #return Fxyz_local_2, Mxyz_local_2
    return Fxyz_local_2

if __name__ == '__main__':  # pragma: no cover
    #print(iformat('4si3f'))
    test_abs_max_min_global()
    test_abs_max_min_vector()
