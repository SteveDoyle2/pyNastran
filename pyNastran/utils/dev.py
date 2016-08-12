from __future__ import print_function
import os
from six import string_types, iteritems

from numpy import array, ndarray

from pyNastran.utils import object_attributes


def get_files_of_type(dirname, extension='.txt', max_size=100., limit_file='no_dig.txt'):
    """
    Gets the list of all the files with a given extension in the specified directory

    Parameters
    ----------
    dirname : str
        the directory name
    extension : str; default='.txt'
        list of filetypes to get
    max_size : float; default=100.0
        size in MB for max file size
    limit_file : str; default=no_dig.txt
        the presence of this file indicates no folder digging
        should be done on this folder

    Returns
    -------
    files : List[str]
        list of all the files with a given extension in the specified directory
    """
    if not os.path.exists(dirname):
        return []

    filenames2 = []
    filenames = os.listdir(dirname)
    allow_digging = True
    if limit_file in filenames:
        allow_digging = False
    for f in filenames:
        filename = os.path.join(dirname, f)
        if os.path.isdir(filename):
            if allow_digging:
                filenames2 += get_files_of_type(filename, extension, max_size)
                #assert len(filenames2) > 0, dirnamei
            else:
                print('no digging in filename=%s; dirname=%s' % (filename, dirname))
        elif (os.path.isfile(filename) and
              os.path.splitext(f)[1].endswith(extension) and
              os.path.getsize(filename) / 1048576. <= max_size):
            filenames2.append(filename)
    return filenames2
    #return [os.path.join(dirname, f) for f in os.listdir(dirname)
    #        if os.path.splitext(f)[1].endswith(extension)
    #         and os.path.getsize(os.path.join(dirname, f)) / (1048576.) <= max_size]


def list_print(lst, float_fmt='%-4.2f'):
    """
    Prints a list, numpy array, or numpy matrix in an abbreviated format.
    Supported element types: None, string, numbers. Useful for debugging.

    Parameters
    ----------
    lst : list / numpy array / numpy matrix

    Returns
    -------
    msg : str
        the clean string representation of the object
    """
    def _print(val):
        if val is None or isinstance(val, string_types):
            return str(val)
        if isinstance(val, float):
            return float_fmt % val
        try:
            return '%g' % val
        except TypeError:
            print("parameter = %r" % val)
            raise

    try: # TODO: remove try block and fix bug in OP2 code or add a warning message
        if len(lst) == 0:
            return '[]'

        if isinstance(lst, ndarray) and lst.ndim == 2:
            row, col = lst.shape
            return ("["+",\n ".join(["["+",".join([float_fmt % lst[i, j]
                    for j in range(col)])+"]" for i in range(row)])+"]")
        return "[" + ", ".join([_print(a) for a in lst]) + "]"
    except: # not a list
        return _print(lst)


def write_class(name, obj, nspaces=0, nbase=0):
    object_type = obj.__class__.__name__
    obj_attrs = object_attributes(obj, 'both')
    if not obj_attrs:
        return "%s()" % object_type

    spaces = ' ' * nspaces
    nspaces2 = nspaces + 4
    #spaces2 = nspaces2 * ' '
    msg = "%s(\n" % object_type
    for attr in obj_attrs[:-1]:
        value = getattr(obj, attr)
        #msg += '?'
        msg += write_object_attributes(attr, value, nspaces2, nbase, is_class=True)
    attr = obj_attrs[-1]
    value = getattr(obj, attr)
    msg += write_object_attributes(attr, value, nspaces2, nbase, is_class=True)
    msg += '%s)' % spaces

    #print("dir(obj) =", dir(obj))
    #print("obj_attrs =", obj_attrs)
    return msg


def write_object_attributes(attr, obj, nspaces, nbase=0, is_class=False):
    msg = ''
    if isinstance(obj, (int, float)) or obj is None:
        msg += '%s' % (str(obj))
    elif isinstance(obj, string_types):
        msg += "'%s'" % obj
    elif isinstance(obj, unicode):
        msg += "u'%s'" % obj
    elif isinstance(obj, list):
        msg += write_list(obj, nspaces, nbase, is_class)
    elif isinstance(obj, tuple):
        msg += write_tuple(obj, nspaces, nbase, is_class)
    elif isinstance(obj, dict):
        msg += write_dict(obj, nspaces, nbase, is_class)
    else:
        object_type = type(obj)
        #raise RuntimeError('object_type=%s is not supported; value=%s' % (objec_type, obj))
    return msg


def write_dict(obj, nspaces, nbase, is_class):
    spaces = (nbase+nspaces) * ' '
    nspaces2 = nspaces + 4
    if len(obj) == 0:
        return '{}'

    msg = '{\n'
    for key, value in sorted(iteritems(obj)):
        #msg += '#'
        msg += write_object_attributes(key, value, nspaces2, nbase, is_class=False)
    msg += '%s}' % spaces
    return msg


def write_list(obj, nspaces, nbase, is_class):
    if len(obj) == 0:
        return '[]'
    return ' ' * (nspaces + nbase) + "???"  # don't choke on long numpy arrays

    #spaces = ' ' * (nspaces + nbase)
    #msg = '[\n%s    ' % spaces
    #for value in obj[:-1]:
        #msg += write_value(value, nspaces+4, nbase, is_class) + ', '
    #msg += write_value(obj[-1], nspaces+4, nbase, is_class) + '\n%s]' % spaces
    #return msg


def write_tuple(obj, nspaces, nbase, is_class):
    msg = '('
    for value in obj:
        msg += write_value(value, nspaces, nbase, is_class) + ', '
    msg += ')'
    return msg


def write_array(a, nspaces=0):
    return ' '*nspaces + '[???]'
    shape = a.shape
    dtype = a.dtype
    if len(shape) == 1:
        msg = 'array(['
        #print("a = ",a)
        for ai in a[:-1]:
            #print("ai = ",ai)
            if isinstance(ai, (int, float)):
                msg += '%s, ' % ai
            elif isinstance(ai, string_types):
                msg += "'%s'," % ai
            else:
                object_type = type(ai)
                msg = 'object_type=%s is not supported; value=%s' % (object_type, ai)
                raise RuntimeError(msg)
                #return "'array(.not supported type.)'"
            msg += '%s, ' % ai
        if len(a) > 0:
            if isinstance(a[-1], (int, float)):
                msg += "%s], dtype='%s')" % (a[-1], dtype)
            elif isinstance(a[-1], string_types):
                msg += "'%s'], dtype='%s')" % (a[-1], dtype)
            else:
                object_type = type(ai)
                msg = 'object_type=%s is not supported; value=%s' % (object_type, ai)
                raise RuntimeError(msg)
                #return "'array(.not supported type.)'"
        else:
            msg += '], dtype=%s)' % dtype
    elif len(shape) == 2:
        spaces = ' '*nspaces
        msg = 'array(['
        for i, ai in enumerate(a):
            if i > 0:
                msg += '%s[' % spaces
            for bi in ai[:-1]:
                msg += '%s, ' % bi
            msg += '%s' % ai[-1]

            if i+1 == len(a):
                msg += '], dtype=%s)' % dtype
            else:
                msg += '],\n'
    elif len(shape) == 3:
        return "'array(.not supported shape.)'"
    return msg

def main():
    from numpy import zeros
    class C(object):
        def __init__(self):
            pass
    class B(object):
        def __init__(self, x=None, e=None):
            self.x = 4
            self.e = C()

    class A(object):
        def __init__(self, a=None, b=None, c=None, d=None):
            self.a = a
            self.b = b
            self.c = c
            self.d = {
                'a' : 4,
                'b' : [1, 2, 3],
                'c' : {1 : 2},
                'd' : B(),
                (1, 2) : 4,
            }

    z = zeros(2, dtype='float64')
    #print(z)
    #print(z.dtype)
    dict_a = {
        'strString' : 'a string',
        'strFloat' : 1.0,
        'strInt': 2,
        'strTuple': (1, 2),
        'strNone' : None,
        'strClass' : A('a', 'b', 'c'),
        'strList' : [1, 2, 3],
        'nullList' : [],
        'nullArray' : array([]),
        'stringArray' : array(['s']),
        'stringArray2' : array(['a', 'b']),
        'nullDict' : {},
        u'unicodStr' : u'',
        'ListOfLists' : [[[], [[]], 2, {'a':3}]],
        1 : 1,
        None : 4,
        1.0 : 5,
        (1, 2) : 6,
        'strArray' : array([4, 5, 6]),
        'strArray2' : zeros((2, 2)),
        'strArray3' : zeros((2, 2, 2)),
    }
    dict_b = {
        'string2' : 'a string',
        'float2' : 1.0,
        'int2': 2,
        'dictA' : dict_a,
    }


    dict_c = {
        'dictA' : {
            None : 4,
            1 : 5,
            'strClass' : A(a='a', b='b', c='c'),
            'strFloat' : 1.0,
            'strInt' : 2,
            'strNone' : None,
            'strString' : 'a string',
            'strTuple' : (1, 2),
            (1, 2) : 6,
        },
        'float2' : 1.0,
        'int2' : 2,
        'string2' : 'a string',
    }
    #assert sorted(dictB.items())==sorted(dictC.items())
    #print(write_object_attributes('dictA', dictA, is_class=False))
    nspaces = 0
    msg = write_object_attributes('dictB', dict_b, nspaces, nbase=0)
    print(msg)
    with open('junk.py', 'wb') as file_obj:
        file_obj.write(msg)

    import junk
    #dictB2 = eval(msg)

if __name__ == '__main__':  # pragma: no cover
    main()
