import numpy as np


def build_groups(objs, name, is_element=False):
    groups = {}
    Types = []
    for obj in objs:
        if obj is None:
            continue
        if hasattr(obj, '_get_types'):
            #print(obj.__class__.__name__)
            Types2 = obj._get_types(nlimit=False)
            Types += Types2
            #Types += [Type.type for Type in Types2]
        else:
            #Types += [obj.type]
            Types += [obj]
    for Type in Types:
        #if Type is None:
            #continue
        group_data = getattr(Type, name)
        if not isinstance(group_data, np.ndarray):
            msg = 'Type %s does not return an ndarray when %s is requested' % (Type.type, name)
            raise RuntimeError(msg)

        msg = 'class %s has a type of %r' % (Type.__class__.__name__, Type.type)
        assert Type.__class__.__name__ == Type.type, msg
        #if is_element:
            #assert hasattr(Type, 'op2_id'), 'class %s has no attribute op2_id' % Type.type

        if len(group_data):
            groups[Type.type] = group_data
    #print("groups = %s" % groups)
    return groups

def asarray(ids, dtype='float32', order=None):
    """
    np.asarray method that actually works...

    Parameters
    ----------
    ids : array_like
        Input data, in any form that can be converted to an array.  This
        includes lists, lists of tuples, tuples, tuples of tuples, tuples
        of lists and ndarrays.
    dtype : data-type, optional (default=float32)
        By default, the data-type is inferred from the input data.
    order : {'C', 'F'}, optional
        Whether to use row-major (C-style) or
        column-major (Fortran-style) memory representation.
        Defaults to 'C'.

    Returns
    -------
    out : ndarray
        Array interpretation of `ids`.  No copy is performed if the input
        is already an ndarray.  If `ids` is a subclass of ndarray, a base
        class ndarray is returned.
    """
    if isinstance(ids, np.ndarray):
        return ids
    elif isinstance(ids, (list, tuple)):
        return np.asarray(ids, dtype=dtype)
    elif isinstance(ids, set):
        return np.asarray(list(ids), dtype=dtype)
    else:
        raise NotImplementedError(type(ids))
