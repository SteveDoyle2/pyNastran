from numpy import array, ndarray
#from pyNastran.femutils.utils import unique2d

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
