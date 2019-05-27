"""
Numpy io functions are defined in this file.
This includes:
 - loadtxt_nice(filename, delimiter=None, skiprows=0, comments='#', dtype=np.float64,
                converters=None, usecols=None, unpack=False,
                ndmin=0,)
 - savetxt_nice(fname, X, fmt='%.18e', delimiter=' ', newline='\n', header='',
                footer='', comments='# ')

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

__all__ = ['loadtxt_nice', 'savetxt_nice']

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

def savetxt_nice(fname, X, fmt='%.18e', delimiter=' ', newline='\n', header='',
                 footer='', comments='# ', encoding=None):
    """
    Save an array to a text file.  This is 95% a backport of numpy 1.15.1's
    savetxt.  It does not support:
      - DataSource's URL
      - pathlib fnames
      - gz files

    Parameters
    ----------
    fname : filename or file handle
        If the filename ends in ``.gz``, the file is automatically saved in
        compressed gzip format.  `loadtxt` understands gzipped files
        transparently.
    X : 1D or 2D array_like
        Data to be saved to a text file.
    fmt : str or sequence of strs, optional
        A single format (%10.5f), a sequence of formats, or a
        multi-format string, e.g. 'Iteration %d -- %10.5f', in which
        case `delimiter` is ignored. For complex `X`, the legal options
        for `fmt` are:

        * a single specifier, `fmt='%.4e'`, resulting in numbers formatted
          like `' (%s+%sj)' % (fmt, fmt)`
        * a full string specifying every real and imaginary part, e.g.
          `' %.4e %+.4ej %.4e %+.4ej %.4e %+.4ej'` for 3 columns
        * a list of specifiers, one per column - in this case, the real
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
    encoding : {None, str}, optional
        Encoding used to encode the outputfile. Does not apply to output
        streams. If the encoding is something other than 'bytes' or 'latin1'
        you will not be able to load the file in NumPy versions < 1.14. Default
        is 'latin1'.

        .. versionadded:: 1.14.0


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

    class WriteWrap:
        """Convert to unicode in py2 or to bytes on bytestream inputs."""
        def __init__(self, fh, encoding):
            self.fh = fh
            self.encoding = encoding
            self.do_write = self.first_write

        def close(self):
            self.fh.close()

        def write(self, v):
            self.do_write(v)

        def write_bytes(self, v):
            if isinstance(v, bytes):
                self.fh.write(v)
            else:
                self.fh.write(v.encode(self.encoding))

        def write_normal(self, v):
            self.fh.write(asunicode(v))

        def first_write(self, v):
            try:
                self.write_normal(v)
                self.write = self.write_normal
            except TypeError:
                # input is probably a bytestream
                self.write_bytes(v)
                self.write = self.write_bytes

    own_fh = False
    #if is_pathlib_path(fname):
        #fname = str(fname)

    if _is_string_like(fname):
        # datasource doesn't support creating a new file ...
        open(fname, 'wt').close()
        fh = open(fname, 'wt', encoding=encoding)
        own_fh = True
        # need to convert str to unicode for text io output
        if sys.version_info[0] == 2:
            fh = WriteWrap(fh, encoding or 'latin1')
    elif hasattr(fname, 'write'):
        # wrap to handle byte output streams
        fh = WriteWrap(fname, encoding or 'latin1')

    #if _is_string_like(fname):
        #own_fh = True
        #if fname.endswith('.gz'):
            #import gzip
            #fh = gzip.open(fname, 'wb')
        #else:
            #if sys.version_info[0] >= 3:
                #fh = open(fname, 'wb')
            #else:
                #fh = open(fname, 'w')
    #elif hasattr(fname, 'write'):
        #fh = fname
    else:
        raise ValueError('fname must be a string or file handle')

    try:
        X = np.asarray(X)

        # Handle 1-dimensional arrays
        if X.ndim == 0 or X.ndim > 2:
            raise ValueError(
                "Expected 1D or 2D array, got %dD array instead" % X.ndim)
        elif X.ndim == 1:
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
            fh.write(comments + header + newline)
        if iscomplex_X:
            for row in X:
                row2 = []
                for number in row:
                    row2.append(number.real)
                    row2.append(number.imag)
                s = txt_format % tuple(row2) + newline
                fh.write(s.replace('+-', '-'))
        else:
            for row in X:
                try:
                    v = txt_format % tuple(row) + newline
                except TypeError:
                    raise TypeError("Mismatch between array dtype ('%s') and "
                                    "format specifier ('%s')"
                                    % (str(X.dtype), txt_format))
                fh.write(v)

        if len(footer) > 0:
            footer = footer.replace('\n', '\n' + comments)
            fh.write(comments + footer + newline)
    finally:
        if own_fh:
            fh.close()
