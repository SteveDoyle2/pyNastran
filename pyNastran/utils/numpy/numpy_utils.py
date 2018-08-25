"""
Various numpy dependent mathematical functions are defined in this file.
This includes:
 - is_array = isfinite(array)
 - is_array = isfinite_and_greater_than(array, value)
 - is_array = isfinite_and_nonzero(array)
 - c = cross2d(a, b)
 - row_col_pairs = unique2d(array)
 - row_col_pairs, optional_index, optional_inverse = unique_rows(
       return_index=False, return_inverse=False):
 - augmented_identity(A)
 - loadtxt_nice(filename, delimiter=None, skiprows=0, comment='#', dtype=np.float64,
                converters=None, usecols=None, unpack=False,
                ndmin=0,)
 - savetxt_nice(fname, X, fmt='%.18e', delimiter=' ', newline='\n', header='',
                footer='', comments='# ')
"""
from __future__ import print_function
import sys
from codecs import open as codec_open
from itertools import count

from six import StringIO
import numpy as np
#from numpy._iotools import _is_string_like
from numpy.lib._iotools import _is_string_like
from numpy.compat import asstr, asbytes

from pyNastran.utils import is_file_obj, _filename

def isfinite(_array):
    """are any of the values finite?"""
    return np.any(np.isfinite(_array))

def isfinite_and_greater_than(_array, value):
    """are any of the values finite and greater than some value?"""
    return isfinite(_array) and abs(np.nanmax(_array) > value)

def isfinite_and_nonzero(_array):
    """are any of the values finite and a value non-zero?"""
    return isfinite_and_greater_than(_array, 0.)

#ver = np.lib.NumpyVersion(np.__version__)
#if ver < '1.13.0':

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
    idx = np.unique(b, return_index=True)[1]
    return a[idx]

#def unique_rows(data):
    #"""
    #finds the unique rows of a numpy array
    #"""
    #uniq = unique(data.view(data.dtype.descr * data.shape[1]))
    #return uniq.view(data.dtype).reshape(-1, data.shape[1])

def unique_rows(A, return_index=False, return_inverse=False):
    """
    Similar to MATLAB's unique(A, 'rows'), this returns B, I, J
    where B is the unique rows of A and I and J satisfy
    A = B[J,:] and B = A[I,:]

    Returns
    -------
    I : ndarray?
        the index array;
        returns if return_index=True
    J : ndarray?
        the inverse array;
        returns if return_inverse=True

    Example (not tested)
    --------------------
    >>> B       = unique_rows(A, return_index=False, return_inverse=False)
    >>> B, I    = unique_rows(A, return_index=True,  return_inverse=False)
    >>> B, J    = unique_rows(A, return_index=False, return_inverse=True)
    >>> B, I, J = unique_rows(A, return_index=True,  return_inverse=True)

    per https://github.com/numpy/numpy/issues/2871
    """
    A = np.require(A, requirements='C')
    assert A.ndim == 2, "array must be 2-dim'l"

    B = np.unique(A.view([('', A.dtype)] * A.shape[1]),
                  return_index=return_index,
                  return_inverse=return_inverse)

    if return_index or return_inverse:
        return (B[0].view(A.dtype).reshape((-1, A.shape[1]), order='C'),) \
            + B[1:]
    else:
        return B.view(A.dtype).reshape((-1, A.shape[1]), order='C')

def cross2d(a, b):
    """
    Interface to np.cross for 2d matrices

    [cx0, cy0, cz0] = [ax0, ay0, az0]   [bx0, by0, bz0]
    |cx1, cy1, cz1| = |ax1, ay1, az1| x |bx1, by1, bz1|
    |cx2, cy2, cz2| = |ax2, ay2, az2|   |bx2, by2, bz2|
    [cx3, cy3, cz3] = [ax3, ay3, az3]   [bx3, by3, bz3]
    """
    # axisa=-1, axisb=-1, axisc=-1,
    return np.cross(a, b, axis=1)

def augmented_identity(nx, ny):
    """
    Creates an Identity Matrix augmented with zeros.
    The location of the extra zeros depends on nx/ny.

    .. code-block:: python

      [ 1, 0, 0, 0 ]
      [ 0, 1, 0, 0 ]
      [ 0, 0, 1, 0 ]
    """
    eye = np.eye(max(nx, ny), dtype='float64')
    return eye[:nx, :ny]

def perpendicular_vector(v):
    """
    Finds an arbitrary perpendicular vector to *v*.

    Parameters
    -----------
    v : (3, ) float ndarray
        Vector

    Returns
    --------
    u : (3, ) float ndarray
        Perpendicular vector

    """
    # for two vectors (x, y, z) and (a, b, c) to be perpendicular,
    # the following equation has to be fulfilled
    #     xyz dot abc = 0
    #     0 = ax + by + cz = 0

    # x = y = z = 0 is not an acceptable solution
    if v[0] == v[1] == v[2] == 0.:
        raise ValueError('zero-vector')

    # If one dimension is zero, this can be solved by setting that to
    # non-zero and the others to zero. Example: (4, 2, 0) lies in the
    # x-y-Plane, so (0, 0, 1) is orthogonal to the plane.
    if v[0] == 0:
        return np.array([1., 0., 0.])
    if v[1] == 0:
        return np.array([0., 1., 0.])
    if v[2] == 0:
        return np.array([0., 0., 1.])

    # arbitrarily set a = b = 1
    # then the equation simplifies to
    #     c = -(x + y)/z
    return np.array([1., 1., -1.0 * (v[0] + v[1]) / v[2]])

def norm2d(v):
    mag = np.linalg.norm(v, axis=1)
    assert v.shape[0] == len(mag)
    return mag

def normalize_vector2d(v):
    mag = norm2d(v)
    nmag = len(mag)
    i = v / mag[:, np.newaxis]
    return i, nmag

def axes_stack(i, j, k, nmag):
    i.shape = (nmag, 1, 3)
    j.shape = (nmag, 1, 3)
    k.shape = (nmag, 1, 3)
    coords = np.hstack([i, j, k])
    return coords

def perpendicular_vector2d(v_array):
    """
    Finds an array of arbitrary perpendicular vector to *v_array*.

    Parameters
    -----------
    v_array : (3, ) or (n, 3) float ndarray
        Vector

    Returns
    --------
    v_array : (n, 3) float ndarray
        Perpendicular vector

    """
    v = np.atleast_2d(v_array)
    # for two vectors (x, y, z) and (a, b, c) to be perpendicular,
    # the following equation has to be fulfilled
    #     xyz dot abc = 0
    #     0 = ax + by + cz = 0

    # x = y = z = 0 is not an acceptable solution
    v1 = v[:, 0]
    v2 = v[:, 1]
    v3 = v[:, 2]

    v1_equals_0 = v1 == 0.
    v2_equals_0 = v2 == 0.
    v3_equals_0 = v3 == 0.
    ibad = np.where(v1_equals_0 & v2_equals_0 & v3_equals_0)[0]
    if ibad:
        raise ValueError('zero-vector')

    vout = np.full(v.shape, np.nan, dtype=None, order='C')
    # If one dimension is zero, this can be solved by setting that to
    # non-zero and the others to zero. Example: (4, 2, 0) lies in the
    # x-y-Plane, so (0, 0, 1) is orthogonal to the plane.
    iv1 = np.where(v1_equals_0)[0]
    iv2 = np.where(v2_equals_0)[0]
    iv3 = np.where(v3_equals_0)[0]
    if len(iv1):
        vout[iv1] = [1., 0., 0.]
    if len(iv2):
        vout[iv2] = [0., 1., 0.]
    if len(iv3):
        vout[iv3] = [0., 0., 1.]

    all_true = v1_equals_0 | v2_equals_0 | v3_equals_0
    #print('iv1 =', iv1)
    #print('iv2 =', iv2)
    #print('iv3 =', iv3)
    #print('all_true =', all_true)
    if np.all(all_true):
        return vout
    is_3d = np.where(~all_true)[0]
    #print('----------')
    #print('is_3d =', is_3d)
    #print(vout)

    # arbitrarily set a = b = 1
    # then the equation simplifies to
    #     c = -(x + y)/z
    vout[is_3d, :2] = 1.
    vout[is_3d, 2] = -1. * (v[is_3d, 0] + v[is_3d, 1]) / v[is_3d, 2]
    return vout

def dot3d(A, B, debug=True):
    """
    Multiplies two N x 3 x 3 matrices together

    Parameters
    ----------
    A, B : (n, 3, 3) float ndarray
        the set of matrices to multiply

    Returns
    -------
    C : (n, 3, 3) float ndarray
        the set of 3 x 3 matrix multiplies
    """
    C = np.matmul(A, B)
    if debug:
        dtype = A.dtype
        print('------------------------')
        D = np.zeros(A.shape, dtype=dtype)
        for i, Ai, Bi in zip(count(), A, B):
            D[i, :, :] = Ai.dot(Bi)
            #print(D[i, :, :])
            #print('------------------------')
    assert np.all(np.allclose(C, D))
    return C

def transpose3d(T):
    """
    Returns the transpose in 3d

    Parameters
    ----------
    A : (n, a, b)
        the set of matrices to transpose

    Returns
    -------
    transpose : (n, b, a) float ndarray
        the transposed matrix
    """
    return np.transpose(T, axes=(1, 2))

def triple(A, T):
    """
    Calculates the matrix triple product  for a series of::

        triple[n, :, :] = T.T @ A @ T

    Parameters
    ----------
    A, T : (n, 3, 3)
        the set of matrices to multiply

    Returns
    -------
    triple : (n, 3, 3) float ndarray
        the set of 3 x 3 matrix triples

    """
    return np.matmul(transpose3d(T), np.matmul(A, T))

def triple_transpose(A, T):
    """
    Calculates the matrix triple product for a series of::

        triple[n, :, :] = T @ A @ T.T

    Parameters
    ----------
    A, T : (n, 3, 3)
        the set of matrices to multiply

    Returns
    -------
    triple : (n, 3, 3) float ndarray
        the set of 3 x 3 matrix triples

    """
    return np.matmul(T, np.matmul(A, transpose3d(T)))

def loadtxt_nice(filename, delimiter=None, skiprows=0, comment='#', dtype=np.float64,
                 converters=None, usecols=None, unpack=False,
                 ndmin=0,):
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
    complex_dtypes = [
        'complex128', np.complex128,
    ]
    if dtype in complex_dtypes:
        return np.loadtxt(
            filename, dtype=dtype, comments=comment, delimiter=delimiter,
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
        with codec_open(_filename(filename), 'r') as file_obj:
            if skiprows:
                lines = file_obj.readlines()[skiprows:]
            else:
                lines = file_obj.readlines()
    elif is_file_obj(filename):
        lines = filename.readlines()[skiprows:]
        filename = filename.name
    else:  # pragma: no cover
        raise TypeError('filename=%s is not a file-like object' % filename)

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
    allowed_float_dtypes = [
        'int32', 'int64', 'int128',
        'float32', 'float64', 'float128', np.float64,
    ]
    #if dtype not in allowed_float_dtypes:  # pragma: no cover
        #'dtype=%r allowed_float_dtypes=[%s]' % (
            #dtype, ', '.join(allowed_float_dtypes))
        #raise RuntimeError(msg)

    if dtype in allowed_float_dtypes:
        X = np.array(data, dtype=dtype)
    #elif dtype in complex_dtypes:
        #return np.loadtxt(fname, dtype=dtype, comments=comment, delimiter=delimiter,
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
    else:
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
        if dtypei not in allowed_dtypes:
            raise RuntimeError('dtype=%r allowed_dtypes=[%s]' % (
                dtypei, ', '.join(allowed_dtypes)))
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
            if dtypei in ['float32', 'float64', 'float128', np.float64]:
                for irow, val in zip(count(), a[:, icol]):
                    try:
                        float(val)
                    except ValueError:
                        msg += 'for name=%r, row=%s -> val=%r (expected float)\n' % (
                            name, irow, val)
                        is_failed = True
            elif dtypei in ['int32', 'int64', 'int128']:
                for irow, val in zip(count(), a[:, icol]):
                    try:
                        int(val)
                    except ValueError:
                        msg += 'for name=%r, row=%s -> val=%r (expected int)\n' % (
                            name, irow, val)
                        is_failed = True
            else:
                raise NotImplementedError(dtype)
            if is_failed:
                raise RuntimeError(msg)
    return X

def savetxt_nice(fname, X, fmt='%.18e', delimiter=' ', newline='\n', header='',
                 footer='', comments='# '):
    """
    Reimplmenentation of numpy's savetxt that doesn't complain about
    bytes when saving to unicode files in Python 3.

    Save an array to a text file.

    Parameters
    ----------
    fname : filename or file handle
        If the filename ends in ``.gz``, the file is automatically saved in
        compressed gzip format.  `loadtxt` understands gzipped files
        transparently.
    X : array_like
        Data to be saved to a text file.
    fmt : str or sequence of strs, optional
        A single format (%10.5f), a sequence of formats, or a
        multi-format string, e.g. 'Iteration %d -- %10.5f', in which
        case `delimiter` is ignored. For complex `X`, the legal options
        for `fmt` are:
            a) a single specifier, `fmt='%.4e'`, resulting in numbers formatted
                like `' (%s+%sj)' % (fmt, fmt)`
            b) a full string specifying every real and imaginary part, e.g.
                `' %.4e %+.4j %.4e %+.4j %.4e %+.4j'` for 3 columns
            c) a list of specifiers, one per column - in this case, the real
                and imaginary part must have separate specifiers,
                e.g. `['%.3e + %.3ej', '(%.15e%+.15ej)']` for 2 columns
    delimiter : str, optional
        String or character separating columns.
    newline : str, optional
        String or character separating lines.

        .. versionadded:: 1.5.0
    header : str, optional
        String that will be written at the beginning of the file.

        .. versionadded:: 1.7.0
    footer : str, optional
        String that will be written at the end of the file.

        .. versionadded:: 1.7.0
    comments : str, optional
        String that will be prepended to the ``header`` and ``footer`` strings,
        to mark them as comments. Default: '# ',  as expected by e.g.
        ``numpy.loadtxt``.

        .. versionadded:: 1.7.0


    See Also
    --------
    save : Save an array to a binary file in NumPy ``.npy`` format
    savez : Save several arrays into an uncompressed ``.npz`` archive
    savez_compressed : Save several arrays into a compressed ``.npz`` archive

    Notes
    -----
    Further explanation of the `fmt` parameter
    (``%[flag]width[.precision]specifier``):

    flags:
        ``-`` : left justify

        ``+`` : Forces to precede result with + or -.

        ``0`` : Left pad the number with zeros instead of space (see width).

    width:
        Minimum number of characters to be printed. The value is not truncated
        if it has more characters.

    precision:
        - For integer specifiers (eg. ``d,i,o,x``), the minimum number of
          digits.
        - For ``e, E`` and ``f`` specifiers, the number of digits to print
          after the decimal point.
        - For ``g`` and ``G``, the maximum number of significant digits.
        - For ``s``, the maximum number of characters.

    specifiers:
        ``c`` : character

        ``d`` or ``i`` : signed decimal integer

        ``e`` or ``E`` : scientific notation with ``e`` or ``E``.

        ``f`` : decimal floating point

        ``g,G`` : use the shorter of ``e,E`` or ``f``

        ``o`` : signed octal

        ``s`` : string of characters

        ``u`` : unsigned decimal integer

        ``x,X`` : unsigned hexadecimal integer

    This explanation of ``fmt`` is not complete, for an exhaustive
    specification see [1]_.

    References
    ----------
    .. [1] `Format Specification Mini-Language
           <http://docs.python.org/library/string.html#
           format-specification-mini-language>`_, Python Documentation.

    Examples
    --------
    >>> x = y = z = np.arange(0.0,5.0,1.0)
    >>> np.savetxt('test.out', x, delimiter=',')   # X is an array
    >>> np.savetxt('test.out', (x,y,z))   # x,y,z equal sized 1D arrays
    >>> np.savetxt('test.out', x, fmt='%1.4e')   # use exponential notation

    """
    # Py3 conversions first
    if isinstance(fmt, bytes):
        fmt = asstr(fmt)
    delimiter = asstr(delimiter)

    own_fh = False
    if _is_string_like(fname):
        own_fh = True
        if fname.endswith('.gz'):
            import gzip
            fh = gzip.open(fname, 'wb')
        else:
            if sys.version_info[0] >= 3:
                fh = open(fname, 'wb')
            else:
                fh = open(fname, 'w')
    elif hasattr(fname, 'write'):
        fh = fname
    else:
        raise ValueError('fname must be a string or file handle')

    try:
        X = np.asarray(X)

        # Handle 1-dimensional arrays
        if X.ndim == 1:
            # Common case -- 1d array of numbers
            if X.dtype.names is None:
                X = np.atleast_2d(X).T
                ncol = 1

            # Complex dtype -- each field indicates a separate column
            else:
                ncol = len(X.dtype.descr)
        else:
            ncol = X.shape[1]

        iscomplex_X = np.iscomplexobj(X)
        # `fmt` can be a string with multiple insertion points or a
        # list of formats.  E.g. '%10.5f\t%10d' or ('%10.5f', '$10d')
        if isinstance(fmt, (list, tuple)):
            if len(fmt) != ncol:
                raise AttributeError('fmt has wrong shape.  %s' % str(fmt))
            txt_format = asstr(delimiter).join(map(asstr, fmt))
        elif isinstance(fmt, str):
            n_fmt_chars = fmt.count('%')
            error = ValueError('fmt has wrong number of %% formats:  %s' % fmt)
            if n_fmt_chars == 1:
                if iscomplex_X:
                    fmt = [' (%s+%sj)' % (fmt, fmt), ] * ncol
                else:
                    fmt = [fmt, ] * ncol
                txt_format = delimiter.join(fmt)
            elif iscomplex_X and n_fmt_chars != (2 * ncol):
                raise error
            elif ((not iscomplex_X) and n_fmt_chars != ncol):
                raise error
            else:
                txt_format = fmt
        else:
            raise ValueError('invalid fmt: %r' % (fmt,))

        if len(header) > 0:
            header = header.replace('\n', '\n' + comments)
            fh.write(asbytes(comments + header + newline))
        if iscomplex_X:
            for row in X:
                row2 = []
                for number in row:
                    row2.append(number.real)
                    row2.append(number.imag)
                fh.write(asbytes(txt_format % tuple(row2) + newline))
        else:
            for row in X:
                try:
                    #print('txt_format = %r' % txt_format, type(txt_format))
                    fh.write(asbytes(txt_format % tuple(row) + newline))
                except TypeError:
                    raise TypeError("Mismatch between array dtype ('%s') and "
                                    "format specifier ('%s')"
                                    % (str(X.dtype), txt_format))
        if len(footer) > 0:
            footer = footer.replace('\n', '\n' + comments)
            fh.write(asbytes(comments + footer + newline))
    finally:
        if own_fh:
            fh.close()
