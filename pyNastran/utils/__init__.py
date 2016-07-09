# -*- coding: utf-8 -*-
from __future__ import print_function
from types import MethodType
import os
import io
from codecs import open as codec_open
from itertools import count

from six import PY2, string_types, iteritems, StringIO

import numpy as np
#from numpy import ndarray, array, loadtxt, int32, int64


if PY2:
    integer_types = (int, long, np.int32, np.int64)
    integer_float_types = (int, long, np.int32, np.int64, float)
else:
    integer_types = (int, np.int32, np.int64)
    integer_float_types = (int, np.int32, np.int64, float)

def is_file_obj(filename):
    """does this object behave like a file object?"""
    #if not (hasattr(out_filename, 'read') and hasattr(out_filename, 'write')) or isinstance(out_filename, file) or isinstance(out_filename, StringIO):
    return ((hasattr(filename, 'read') and hasattr(filename, 'write'))
            or isinstance(filename, file)
            or isinstance(filename, StringIO))

def b(string):
    """reimplementation of six.b(...) to work in Python 2"""
    return string.encode('latin-1')

def merge_dicts(dict_list, strict=True):
    """merges two or more dictionaries"""
    assert isinstance(dict_list, list), type(dict_list)
    dict_out = {}
    for adict in dict_list:
        assert isinstance(adict, dict), adict
        for key, value in iteritems(adict):
            if key not in dict_out:
                dict_out[key] = value
            elif strict:
                raise RuntimeError('key=%r exists in multiple dictionaries' % key)
            else:
                print('key=%r is dropped?' % key)
    return dict_out

def loadtxt_nice(filename, delimiter=None, skiprows=0, comment='#', dtype=np.float64,
                 converters=None, usecols=None, unpack=False,
                 ndmin=0,
                 ):
    """
    Reimplmenentation of numpy's loadtxt that doesn't complain about
    training commas (or other delimiter) that vary from  one line to
    the other.  It also provides better error messages.

    Parameters
    ----------
    filename : varies
        str : the filename to load
        file : the file object to load
        cStringIO/StringIO : a file-like object
    delimiter : str; default=None (any whitespace)
        the field splitter (e.g. comma or tab)
    skiprows : int; default=1
        the number of rows to skip
    comment : str, default='#'
        the comment line
    dtype : numpy.dtype; default=None (float)
        allows for alternate casting
        int32, float32, ...
        dtype = {
            names : ('A', 'B', 'C'),
            formats : ('int32', 'float32', 'float64'),
        }
    usecols : sequence; default=None
        Which columns to read, with 0 being the first.  For example,
        ``usecols = (1,4,5)`` will extract the 2nd, 5th and 6th columns.
        The default, None, results in all columns being read.
    unpack : bool, optional
        If True, the returned array is transposed, so that arguments may be
        unpacked using ``x, y, z = loadtxt(...)``.  When used with a structured
        data-type, arrays are returned for each field.  Default is False.

    converters : dict; default=None
        not supported
        crashes if not None
        A dictionary mapping column number to a function that will convert
        that column to a float.  E.g., if column 0 is a date string:
        ``converters = {0: datestr2num}``.  Converters can also be used to
        provide a default value for missing data (but see also `genfromtxt`):
        ``converters = {3: lambda s: float(s.strip() or 0)}``.  Default: None.

    ndmin : int, optional
        crashes if not 0
        The returned array will have at least `ndmin` dimensions.
        Otherwise mono-dimensional axes will be squeezed.
        Legal values: 0 (default), 1 or 2.

    Returns
    -------
    data : (nrows, ncols) ndarray
        the data object
    """
    if converters is not None:
        raise NotImplementedError('converters=%r must be None' % converters)
    #if ndmin is not [0, 2]: ## TODO: remove 2
        #raise NotImplementedError('ndmin=%r must be 0' % ndmin)

    if delimiter is None:
        ending_characters = '\n\r \t'
    else:
        ending_characters = '\n\r \t' + delimiter

    data = []
    if isinstance(filename, StringIO):
        lines = filename.getvalue().split('\n')[skiprows:]
        filename = None
    elif is_file_obj(filename):
        lines = filename.readlines()[skiprows:]
        filename = filename.name
    else:
        with codec_open(_filename(filename), 'r') as file_obj:
            if skiprows:
                lines = file_obj.readlines()[skiprows:]
            else:
                lines = file_obj.readlines()

    if usecols:
        for usecol in usecols:
            assert isinstance(usecol, int), 'usecol=%s usecols=%s' % (usecol, usecols)
        assert len(np.unique(usecols)), 'usecols=%s must be unique' % (usecols)
        for line in lines:
            if line.startswith(comment):
                continue
            sline = line.strip(delimiter).split(delimiter)
            data.append([sline[i] for i in usecols])
    else:
        for line in lines:
            if line.startswith(comment):
                continue
            sline = line.strip(delimiter).split(delimiter)
            data.append(sline)
    del lines

    #print(data)
    allowed_dtypes = ['float32', 'float64', 'float128', np.float64, 'int32', 'int64', 'int128']
    if dtype in allowed_dtypes:
        assert dtype in allowed_dtypes, 'dtype=%r allowed_dtypes=[%s]' % (dtype, ', '.join(allowed_dtypes))
        X = np.array(data, dtype=dtype)
    elif isinstance(dtype, dict):
        a = np.array(data, dtype=object)

        X = {}
        names = dtype['names']
        nnames = len(names)
        assert len(set(names)) == nnames, 'non-unique headers in %s' % str(names)
        for icol, name, dtypei in zip(count(), dtype['names'], dtype['formats']):
            assert dtypei in allowed_dtypes, 'dtype=%r allowed_dtypes=[%s]' % (dtypei, ', '.join(allowed_dtypes))
            try:
                X[name] = np.asarray(a[:, icol], dtype=dtypei)
            except IndexError:
                # the number of columns in A is not consistent
                ncols = [len(datai) for datai in data]
                ucols = np.unique(ncols)
                msg = 'The number of columns is not consistent; expected=%s; actual=%s' % (nnames, ucols)
                raise IndexError(msg)
            except ValueError:
                print(a)
                # we only allow floats
                msg = ''
                if dtypei in ['float32', 'float64', 'float128', np.float64]:
                    for irow, val in zip(count(), a[:, icol]):
                        try:
                            float(val)
                        except:
                            msg += 'for name=%r, row=%s -> val=%r (expected float)\n' % (name, irow, val)
                            is_failed = True
                elif dtypei in ['int32', 'int64', 'int128']:
                    for irow, val in zip(count(), a[:, icol]):
                        try:
                            int(val)
                        except:
                            msg += 'for name=%r, row=%s -> val=%r (expected int)\n' % (name, irow, val)
                            is_failed = True
                else:
                    raise NotImplementedError(dtype)
                if is_failed:
                    raise RuntimeError(msg)

        #print('A =', A)
    else:
        raise NotImplementedError('dtype_else=%s' % dtype)
        #return np.array(data, dtype=dtype)

    #if usecols is not None:
        #raise NotImplementedError('usecols=%s must be None' % str(usecols))
    #if unpack is not False:
        #raise NotImplementedError('unpack=%r must be False' % unpack)

    if not isinstance(dtype, dict):
        # Tweak the size and shape of the arrays - remove extraneous dimensions
        if X.ndim > ndmin:
            X = np.squeeze(X)
        # and ensure we have the minimum number of dimensions asked for
        # - has to be in this order for the odd case ndmin=1, X.squeeze().ndim=0
        if X.ndim < ndmin:
            if ndmin == 1:
                X = np.atleast_1d(X)
            elif ndmin == 2:
                X = np.atleast_2d(X).T


    if unpack:
        #print(X)
        if isinstance(dtype, dict) > 1:
            if ndmin == 0:
                # For structured arrays, return an array for each field.
                return (np.squeeze(X[name]) for name in dtype.names)
            else:
                raise RuntimeError('I think this can never happen...type(dtype)=dict; ndmin=%s' % ndmin)
        else:
            #print('X = ', X)
            #raise NotImplementedError('unpack=%s dtypes=%s' % (unpack, dtype))
            #if ndmin == 0: # and A.shape[0] == 1
                #out = (np.squeeze(X[:, i]) for i in range(X.shape[1]))
                #return out
            #else:
                #return (X[:, i] for i in range(X.shape[1]))
            #return (X[:, i] for i in range(X.shape[1]))
            return X.T
    else:
        return X
    #return np.array(data)


def is_binary_file(filename):
    """
    Return true if the given filename is binary.

    Parameters
    ----------
    filename : str
        the filename to test

    Returns
    -------
    binary_flag : bool
        True if filename is a binary file (contains null byte)
        and False otherwise.

    :raises:  IOError if the file cannot be opened.

    Based on the idea (.. seealso:: http://bytes.com/topic/python/answers/21222-determine-file-type-binary-text)
    that file is binary if it contains null.

    .. warning:: this may not work for unicode."""
    assert isinstance(filename, string_types), '%r is not a valid filename' % filename
    assert os.path.exists(filename), '%r does not exist\n%s' % (filename, print_bad_path(filename))
    with io.open(filename, mode='rb') as fil:
        for chunk in iter(lambda: fil.read(1024), bytes()):
            if b'\0' in chunk:  # found null byte
                return True
    return False


def print_bad_path(path):
    """
    Prints information about the existence (access possibility) of the parts
    of the given path. Useful for debugging when the path to a given file
    is wrong.

    Parameters
    ----------
    path : str
        path to check

    Returns
    -------
    msg : str
        string with informations whether access to parts of the path
        is possible
    """
    raw_path = path
    if len(path) > 255:
        path = os.path.abspath(_filename(path))
        npath = os.path.dirname(path)
        res = [path]
        while path != npath:
            path, npath = npath, os.path.dirname(npath)
            res.append(path)
        msg = {True: 'passed', False: 'failed'}
        return '\n'.join(['%s: %s' % (msg[os.path.exists(i)], i[4:]) for i in res])
    else:
        path = os.path.abspath(path)
        npath = os.path.dirname(path)
        res = [path]
        while path != npath:
            path, npath = npath, os.path.dirname(npath)
            res.append(path)
        msg = {True: 'passed', False: 'failed'}
        return '\n'.join(['%s: %s' % (msg[os.path.exists(i)], i) for i in res])

def _filename(filename):
    """
    Prepends some magic data to a filename in order to have long filenames.

    .. warning:: This might be Windows specific.
    """
    if len(filename) > 255:
        return '\\\\?\\' + filename
    return filename

def __object_attr(obj, mode, keys_to_skip, attr_type):
    """list object attributes of a given type"""
    #print('keys_to_skip=%s' % keys_to_skip)
    keys_to_skip = [] if keys_to_skip is None else keys_to_skip
    test = {
        'public':  lambda k: (not k.startswith('_') and k not in keys_to_skip),
        'private': lambda k: (k.startswith('_') and not k.startswith('__') and k not in keys_to_skip),
        'both': lambda k: (not k.startswith('__') and k not in keys_to_skip),
        'all':  lambda k: (k not in keys_to_skip),
    }

    if not mode in test:
        print('Wrong mode! Accepted modes: public, private, both, all.')
        return None
    check = test[mode]

    out = []
    for k in dir(obj):
        if k in keys_to_skip:
            continue
        if check(k) and attr_type(getattr(obj, k)):
            out.append(k)
    out.sort()
    return out
    #return sorted([k for k in dir(obj) if (check(k) and
    #                                           attr_type(getattr(obj, k)))])


def object_methods(obj, mode='public', keys_to_skip=None):
    """
    List the names of methods of a class as strings. Returns public methods
    as default.

    Parameters
    ----------
    obj : instance
        the object for checking
    mode : str
        defines what kind of methods will be listed
        * "public" - names that do not begin with underscore
        * "private" - names that begin with single underscore
        * "both" - private and public
        * "all" - all methods that are defined for the object
    keys_to_skip : List[str]; default=None -> []
        names to not consider to avoid deprecation warnings

    Returns
    -------
    method : List[str]
        sorted list of the names of methods of a given type
        or None if the mode is wrong
    """
    return __object_attr(obj, mode, keys_to_skip, lambda x: isinstance(x, MethodType))


def object_attributes(obj, mode='public', keys_to_skip=None):
    """
    List the names of attributes of a class as strings. Returns public
    attributes as default.

    Parameters
    ----------
    obj : instance
        the object for checking
    mode : str
        defines what kind of attributes will be listed
        * 'public' - names that do not begin with underscore
        * 'private' - names that begin with single underscore
        * 'both' - private and public
        * 'all' - all attributes that are defined for the object

    Returns
    -------
    attribute_names : List[str]
        sorted list of the names of attributes of a given type or None
        if the mode is wrong
    """
    return __object_attr(obj, mode, keys_to_skip, lambda x: not isinstance(x, MethodType))


#def write_object_attributes(name, obj, nspaces=0, nbase=0, is_class=True, debug=False):
    #"""
    #Writes a series of nested objects
    #"""
    #spaces = (nbase + nspaces) * ' '
    #msg = spaces
    #xml = spaces
    #if is_class:
        #equals = '='
    #else:
        #equals = ':'

    #if debug:
        #print('attr=%s equals=%r' % (name, equals))
    ## name
    #if isinstance(obj, dict):
        #if nspaces == 0:
            #msg += '%s %s ' % (name, equals)
        #else:
            #if isinstance(name, tuple):
                #msg += '%s %s ' % (str(name), equals)
            #else:
                #msg += "'%s' %s " % (name, equals)
    #elif isinstance(name, string_types):
        #if is_class:
            #key = '%s' % name
        #else:
            #key = "'%s'" % name
    ## elif isinstance(name, unicode):
        ## if is_class:
            ## key = u'%s' % name
        ## else:
            ## key = "u'%s'" % name
    #elif isinstance(name, (int, float, tuple)) or name is None:
        #key = '%s' % str(name)
    #else:
        #raise RuntimeError('key=%s is not a string.  Type=%s' % (name, type(name)))

    #if debug:
        #print('name=%s type=%s' % (name, type(obj)))

    ## write the object
    #if isinstance(obj, (int, float)) or obj is None:
        #xml += '<name=%s value=%s type=%s>' % (name, obj, type(obj))
        #msg += '%s %s %s,\n' % (key, equals, write_value(obj, nspaces, nbase, is_class))
    #elif is_string(obj):
        #msg += '%s %s %s,\n' % (key, equals, write_value(obj, nspaces, nbase, is_class))

    #elif isinstance(obj, dict):
        #msg += write_dict(obj, nspaces, nbase, is_class) + ',\n'
    #elif isinstance(obj, (tuple, list)):
        #msg += '%s %s %s,\n' % (key, equals, write_value(obj, nspaces, nbase, is_class))

    #elif isinstance(obj, np.ndarray):
        #starter = '%s%s %s' % (nspaces, key, equals)
        #msg += '%s %s %s,\n' % (key, equals, write_array(obj, nspaces + 6 + len(starter)))
    #else:  # generic class
        #objectType = obj.__class__.__name__
        ##raise RuntimeError('objectType=%s is not supported' % objectType)
        #msg += '%s %s ' % (key, equals)
        #msg += write_class(name, obj, nspaces, nbase) + ',\n'  # comma for class
    #if nspaces == 0:
        #msg = msg[:-2]
    #if debug:
        #print('%r' % msg)
    #return msg
