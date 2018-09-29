import numpy as np


def build_groups(objs, name, is_element=False):
    """
    objs : List[CTETRA4, CPENTA6, CHEXA8, ...]
        the objects that are being considered
    name : str
        element_id, etc.
    is_element : bool; default=False
        is this an element
    """
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

    for class_obj in Types:
        #if element_obj is None:
            #continue
        if class_obj.n == 0:
            continue
        group_data = getattr(class_obj, name)
        if not isinstance(group_data, np.ndarray):
            msg = 'Type %s does not return an ndarray when %s is requested' % (class_obj.type, name)
            raise RuntimeError(msg)

        msg = 'class %s has a type of %r' % (class_obj.__class__.__name__, class_obj.type)
        assert class_obj.__class__.__name__ == class_obj.type, msg
        #if is_element:
            #assert hasattr(class_obj, 'op2_id'), 'class %s has no attribute op2_id' % class_obj.type

        if len(group_data):
            groups[class_obj.type] = group_data
    #print("groups = %s" % groups)
    return NdarrayDict(groups)

class NdarrayDict(dict):
    """
    Dictionary subclass to print it in a nicer format.  It puts a return
    line after each key, because the values are ndarrays.
    """
    def __repr__(self):
        """writes the dictionary"""
        msg = 'NdarrayDict(\n'
        for key in sorted(self.keys()):
            value = self[key]
            msg += '%r : %r,\n' % (key, value)
        msg += ')'
        return msg

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
