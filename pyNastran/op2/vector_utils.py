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

    :param values: an ND-array of values
    :type values:  common NDARRAY/list/tuple shapes:
                      1. [nprincipal_stresses]
                      2. [nelements, nprincipal_stresses]

    :returns abs_max_mins: an array of the max or min principal stress
    :type abs_max_mins:    int/float depending on input type
                           (don't input mixed types...)

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

    :param values: an array of values, where the rows are interated over
                   and the columns are going to be compressed
    :type values:  NDARRAY shape=[nelements, nprincipal_stresses]

    :returns abs_max_mins: an array of the max or min principal stress
    :type abs_max_mins:    NDARRAY shape=[nelements] with dtype=values.dtype

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

if __name__ == '__main__':  # pragma: no cover
    #print(iformat('4si3f'))
    test_abs_max_min_global()
    test_abs_max_min_vector()
