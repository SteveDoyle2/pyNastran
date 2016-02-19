from six.moves import range
from numpy import unique, array, ndarray

def slice_to_iter(ids):
    """
    Simplifies possible slice-like inputs

    Parameters
    ----------
    ids : list
        a list of values that will be iterated over and may be a slice

    Returns
    -------
    ids2 : list/ndarray/range
        an iterable as list, ndarray, range
    int_flag : bool
        when integers are passed in, we can return a single object at the end,
        but need a flag to remember

    A slice is something like:
        elements[1:10:2]

    and represented as:
        slice(start, stop, step)

    The __getitem__ method is called with the following as input:
      - slice
      - integer
      - list
      - numpy ndarray

    All these can be used to iterate over an object
    """
    int_flag = False
    if isinstance(ids, int):
        ids2 = array([ids])
        int_flag = True
    elif isinstance(ids, slice):
        if ids.step is None:
            ids2 = range(ids.start, ids.stop)
        else:
            ids2 = range(ids.start, ids.stop, ids.step)
    elif ids is None:
        raise TypeError('cannot turn None into an iterable')
    elif isinstance(ids, (list, ndarray)):
        ids2 = ids
    else:
        raise KeyError(ids)
    return ids2, int_flag

#def unique_rows(data):
    #"""
    #finds the unique rows of a numpy array
    #"""
    #uniq = unique(data.view(data.dtype.descr * data.shape[1]))
    #return uniq.view(data.dtype).reshape(-1, data.shape[1])

def unique2d(a):
    """
    Gets the unique pairs in a 2D vector where the pairs are defined:
    (column 0, column 1).

    Parameters
    ----------
    a : (n,2) ndarray
        the input data

    Returns
    -------
    u : (m,2)
        the unique values in a

    .. note:: this is intended to be used to find unique rows of
              element-id/property-id or property-id/material-id pairs
    .. note:: it works by finding the unique complex numbers and doesn't
              extend well to a 3 column pair
    """
    x, y = a.T
    b = x + y*1.0j
    idx = unique(b, return_index=True)[1]
    return a[idx]
