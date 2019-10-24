"""
Numpy io functions are defined in this file.
This includes:
 - loadtxt_nice(filename, delimiter=None, skiprows=0, comments='#', dtype=np.float64,
                converters=None, usecols=None, unpack=False,
                ndmin=0,)

"""
import sys
from io import StringIO
from itertools import count

import numpy as np
from numpy.lib._iotools import _is_string_like
from numpy.compat import asstr, asunicode#, is_pathlib_path, asunicode

#try:
    #from pathlib import Path
#except ImportError:
    #Path = None

#def is_pathlib_path(obj):
    #"""
    #Check whether obj is a pathlib.Path object.
    #"""
    #return Path is not None and isinstance(obj, Path)


from pyNastran.utils import is_file_obj, _filename

__all__ = ['loadtxt_nice']

def loadtxt_nice(filename, delimiter=None, skiprows=0, comments='#', dtype=np.float64,
                 converters=None, usecols=None, unpack=False,
                 ndmin=0,):
    """
    Reimplementation of numpy's loadtxt that doesn't complain about
    training commas (or other delimiter) that vary from  one line to
    the other.  It also provides better error messages when failing to
    load files.

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
    comments : str, default='#'
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
    complex_dtypes = [
        'complex64', 'complex128', np.complex128, np.complex64,
    ]
    if dtype in complex_dtypes:
        return np.loadtxt(
            filename, dtype=dtype, comments=comments, delimiter=delimiter,
            converters=converters,
            skiprows=skiprows, usecols=usecols,
            unpack=unpack, ndmin=ndmin,
            encoding='bytes')

    if converters is not None:
        raise NotImplementedError('converters=%r must be None' % converters)

    #if ndmin is not [0, 2]: ## TODO: remove 2
        #raise NotImplementedError('ndmin=%r must be 0' % ndmin)

    #if delimiter is None:
        #ending_characters = '\n\r \t'
    #else:
        #ending_characters = '\n\r \t' + delimiter

    data = []
    if isinstance(filename, StringIO):
        lines = filename.getvalue().split('\n')[skiprows:]
        filename = None
    elif isinstance(filename, str):
        with open(_filename(filename), 'r') as file_obj:
            if skiprows:
                lines = file_obj.readlines()[skiprows:]
            else:
                lines = file_obj.readlines()
    elif is_file_obj(filename):
        lines = filename.readlines()[skiprows:]
        filename = filename.name
    else:  # pragma: no cover
        raise TypeError('filename=%s is not a file-like object; type=%s' % (filename, type(filename)))

    if usecols:
        for usecol in usecols:
            assert isinstance(usecol, int), 'usecol=%s usecols=%s' % (usecol, usecols)
        assert len(np.unique(usecols)), 'usecols=%s must be unique' % (usecols)
        for line in lines:
            if line.startswith(comments):
                continue
            sline = line.strip(delimiter).split(delimiter)
            data.append([sline[i] for i in usecols])
    else:
        for line in lines:
            if line.startswith(comments):
                continue
            sline = line.strip(delimiter).split(delimiter)
            data.append(sline)
    del lines

    allowed_float_dtypes = [
        np.float64, np.float32, 'float32', 'float64', 'f4',
        np.int32, np.int64, 'int32', 'int64', 'i4',
        #'float128', np.float128,
        #'int128', np.int128,
    ]
    #if dtype not in allowed_float_dtypes:  # pragma: no cover
        #'dtype=%r allowed_float_dtypes=[%s]' % (
            #dtype, ', '.join(allowed_float_dtypes))
        #raise RuntimeError(msg)

    if dtype in allowed_float_dtypes:
        X = np.array(data, dtype=dtype)
    #elif dtype in complex_dtypes:
        #return np.loadtxt(fname, dtype=dtype, comments=comments, delimiter=delimiter,
                          #converters=converters,
                          #skiprows=skiprows, usecols=usecols,
                          #unpack=unpack, ndmin=ndmin,
                          #encoding='bytes')
        #if dtype not in allowed_complex_dtypes:  # pragma: no cover
            #'dtype=%r allowed_complex_dtypes=[%s]' % (
                #dtype, ', '.join(allowed_complex_dtypes))
            #raise RuntimeError(msg)
        #data2 = (d.strip('()'))
    elif isinstance(dtype, dict):
        X = _loadtxt_as_dict(data, dtype, allowed_float_dtypes)
        #print('A =', A)
    elif isinstance(dtype, list): # tuple
        assert len(data) == len(dtype)
        X = []
        for name_dtypei, datai in zip(dtype, data):
            unused_name, dtypei = name_dtypei
            xi = np.array(datai, dtype=dtypei)
            X.append(xi)
        unpack = False
        dtype = {name_dtypei[0] : name_dtypei[1] for name_dtypei in dtype}
        #X = np.vstack(X)#.astype(dtype)
        #X = X.astype(dtype)
        #print('X =', X)
        #dtype = np.dtype(dtype)
        #X = np.array(data, dtype=dtype)
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
        # - has o be in this order for the odd case ndmin=1, X.squeeze().ndim=0
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
                msg = 'I think this can never happen...type(dtype)=dict; ndmin=%s' % ndmin
                raise RuntimeError(msg)
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
    return X
    #return np.array(data)

def _loadtxt_as_dict(data, dtype, allowed_dtypes):
    """helper method for ``loadtxt_nice``"""
    a = np.array(data, dtype=object)
    X = {}
    names = dtype['names']
    nnames = len(names)
    assert len(set(names)) == nnames, 'non-unique headers in %s' % str(names)
    for icol, name, dtypei in zip(count(), dtype['names'], dtype['formats']):
        if isinstance(dtypei, str) and 's' in dtypei.lower():
            pass
        elif dtypei not in allowed_dtypes:
            allowed_dtypes_str = ', '.join([str(allowed_dtype)
                                            for allowed_dtype in allowed_dtypes])
            raise RuntimeError('dtype=%r allowed_dtypes=[%s]' % (
                dtypei, allowed_dtypes_str))

        try:
            X[name] = np.asarray(a[:, icol], dtype=dtypei)
        except IndexError:
            # the number of columns in A is not consistent
            ncols = [len(datai) for datai in data]
            ucols = np.unique(ncols)
            msg = 'The number of columns is not consistent; expected=%s; actual=%s' % (
                nnames, ucols)
            raise IndexError(msg)
        except ValueError:
            # we only allow floats
            msg = ''
            if dtypei in ['float32', 'float64', 'float128', np.float64, 'f4']:
                for irow, val in zip(count(), a[:, icol]):
                    try:
                        float(val)
                    except ValueError:
                        msg += 'for name=%r, row=%s -> val=%r (expected float)\n' % (
                            name, irow, val)
                        is_failed = True
            elif dtypei in ['int32', 'int64', 'int128', 'i4']:
                for irow, val in zip(count(), a[:, icol]):
                    try:
                        int(val)
                    except ValueError:
                        msg += 'for name=%r, row=%s -> val=%r (expected int)\n' % (
                            name, irow, val)
                        is_failed = True
            else:
                raise NotImplementedError('dtypei=%s dtype=%s' % (dtypei, dtype))
            if is_failed:
                raise RuntimeError(msg)
    return X
