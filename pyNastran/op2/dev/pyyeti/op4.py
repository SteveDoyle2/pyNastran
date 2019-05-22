# -*- coding: utf-8 -*-
"""
Python tools for reading/writing Nastran .op4 files.  Can read and
write all formats (as far as I know) with the restriction that the
output files created by this class are always double precision. The
binary files can be in big or little endian format.

Currently, all matrices are read into dense (non-sparse format)
matrices even if they were written a sparse format.

@author: Tim Widrick
"""
import struct
import warnings
import sys
import itertools as it
import numpy as np


class OP4:
    """
    Class for reading/writing Nastran output4 (.op4) files.

    See demo below and refer to the help on these functions for more
    information: :func:`write`, :func:`dctload`, :func:`listload`,
    :func:`dir`.

    Examples
    --------
    Instantiate the class and create matrices for demo:

    >>> import op4
    >>> o4 = op4.OP4()
    >>> import numpy as np
    >>> r = np.random.randn(3, 5)
    >>> c = 1j*np.random.randn(3, 5) + r

    Write binary op4 file, with 'r' first:

    >>> o4.write('testbin.op4', ['r', 'c'], [r, c])

    Write ascii op4 file without caring about order:

    >>> o4.write('testascii.op4', dict(r=r, c=c), binary=False)

    To read an op4 file into a dictionary (indexed by the name in
    lower case):

    >>> dct = o4.dctload('testbin.op4')

    Or, to read into a list:

    >>> names, mats, forms, mtypes = o4.listload('testascii.op4')

    Check some results:

    >>> print(np.all(r == dct['r'][0]))
    True

    >>> if names[0] == 'c':
    ...     print(np.all(c == mats[0]))
    ... else:
    ...     print(np.all(c == mats[1]))
    True

    To print a 'directory' of an op4 file:

    >>> d = o4.dir('testbin.op4')
    r       ,      3 x 5     , form=2, mtype=2
    c       ,      3 x 5     , form=2, mtype=4

    Clean up:

    >>> import os
    >>> os.remove('testbin.op4')
    >>> os.remove('testascii.op4')

    """

    def __init__(self):
        self._fileh = None
        string = '%.1E' % 1.2
        self._expdigits = len(string) - (string.find('E') + 2)
        self._rows4bigmat = 65536
        # Tunable value ... if number of values exceeds this, read
        # with numpy.fromfile instead of struct.unpack.
        self._rowsCutoff = 3000

    def __del__(self):
        if self._fileh:
            self._fileh.close()
            self._fileh = None

    def _op4close(self):
        if self._fileh:
            self._fileh.close()
            self._fileh = None

    def _op4open_read(self, filename):
        """
        Open binary or ascii op4 file for reading.

        Sets these class variables:

        _fileh : file handle
            Value returned by open().  File is opened in 'r' mode if
            ascii, 'rb' mode if binary.
        _ascii : bool
            True if file is ascii.
        _dformat : bool
            True if an ascii file uses 'D' instead of 'E' (eg, 1.4D3
            instead of 1.4E3)
        _bit64 : True or False
            True if 'key' integers are 64-bit in binary op4 files.
        _endian : string
            Will be '' if no byte-swapping is required; otherwise,
            either '>' or '<' for big-endian and little-endian,
            respectively.  Only used for binary files.
        _Str_i4 : struct.Struct object
            Precompiled for reading 4 byte integers
        _Str_i : struct.Struct object
            Precompiled for reading 4 or 8 byte integers
        _bytes_i : integer
            Either 4 or 8, to go with Str_i.
        _str_sr : string
            Either self._endian + b'%df' or self._endian + b'%dd',
            depending on self._bit64; for reading single precision
            reals.
        _bytes_sr : integer
            Number of bytes in single real.
        _str_dr : string
            self._endian + b'%dd', for reading double precision reals.
        _wordsperdouble : integer
            Either 1 or 2; 2 if self._bit64 is False.
        """
        self._fileh = open(filename, 'rb')
        header_bytes = self._fileh.read(16)
        self._endian = b''
        self._uendian = ''
        self._dformat = False

        # Assuming binary, check for a zero byte in the 'type' field;
        # will have one at front or back if binary:
        if header_bytes[12] == 0 or header_bytes[15] == 0:
            self._ascii = False
            if sys.byteorder == 'little':
                if header_bytes[12] == 0:
                    self._endian = b'>'
                    self._uendian = '>'
            else:
                if header_bytes[12] != 0:
                    self._endian = b'<'
                    self._uendian = '<'
            self._Str_i4 = struct.Struct(self._endian + b'i')
            reclen = self._Str_i4.unpack(header_bytes[:4])[0]
            if reclen == 48:
                self._bit64 = True
                self._Str_i = struct.Struct(self._endian + b'q')
                self._bytes_i = 8
                self._Str_ii = struct.Struct(self._endian + b'qq')
                self._bytes_ii = 16
                self._Str_iii = struct.Struct(self._endian + b'3q')
                self._bytes_iii = 24
                self._Str_iiii = struct.Struct(self._endian + b'4q')
                self._bytes_iiii = 32
                self._str_sr = self._endian + b'%dd'
                self._str_sr_fromfile = np.dtype(self._uendian + 'f8')
                self._bytes_sr = 8
                self._wordsperdouble = 1
            else:
                self._bit64 = False
                self._Str_i = self._Str_i4
                self._bytes_i = 4
                self._Str_ii = struct.Struct(self._endian + b'ii')
                self._bytes_ii = 8
                self._Str_iii = struct.Struct(self._endian + b'3i')
                self._bytes_iii = 12
                self._Str_iiii = struct.Struct(self._endian + b'4i')
                self._bytes_iiii = 16
                self._str_sr = self._endian + b'%df'
                self._str_sr_fromfile = np.dtype(self._uendian + 'f4')
                self._bytes_sr = 4
                self._wordsperdouble = 2
            self._str_dr = self._endian + b'%dd'
            self._str_dr_fromfile = np.dtype(self._uendian + 'f8')
            self._fileh.seek(0)
        else:
            self._ascii = True
            self._fileh.readline()
            self._fileh.readline()
            line = self._fileh.readline().decode()
            # sparse formats have integer header line:
            if line.find('.') == -1:
                line = self._fileh.readline().decode()
            if line.find('D') > -1:
                self._dformat = True
            self._fileh.close()
            self._fileh = open(filename, 'r')

    def _skipop4_ascii(self, perline, rows, cols, mtype, numlen):
        """
        Skip an op4 matrix - ascii.

        Parameters
        ----------
        perline : integer
            Number of elements per line in the file.
        rows : integer
            Number of rows in matrix.
        cols : integer
            Number of columns in matrix.
        mtype : integer
            Nastran matrix type.
        numlen : integer
            String length for each number.

        Returns
        -------
        None

        On entry, file is positioned after the title line, but before
        the first column is printed.  On exit, the file is positioned
        so the next readline will get the next title line.
        """
        # read until next matrix:
        if rows < 0 or rows >= self._rows4bigmat:
            bigmat = True
        else:
            bigmat = False
        if mtype & 1:
            wper = 1
        else:
            wper = 2
        line = self._fileh.readline()
        c = int(line[:8]) - 1
        r = int(line[8:16])
        if r > 0:
            while c < cols:
                elems = int(line[16:24])
                nlines = (elems + perline - 1) // perline
                for _ in it.repeat(None, nlines):
                    self._fileh.readline()
                line = self._fileh.readline()
                c = int(line[:8]) - 1
        elif bigmat:
            while c < cols:
                elems = int(line[16:24])
                while elems > 0:
                    line = self._fileh.readline()
                    L = int(line[:8])-1    # L
                    elems -= L + 2
                    L //= wper
                    # read column as a long string
                    nlines = (L + perline - 1) // perline
                    for _ in it.repeat(None, nlines):
                        self._fileh.readline()
                line = self._fileh.readline()
                c = int(line[:8]) - 1
        else:
            while c < cols:
                elems = int(line[16:24])
                while elems > 0:
                    line = self._fileh.readline()
                    IS = int(line[:8])
                    L = (IS >> 16) - 1        # L
                    elems -= L + 1
                    L //= wper
                    # read column as a long string
                    nlines = (L + perline - 1) // perline
                    for _ in it.repeat(None, nlines):
                        self._fileh.readline()
                line = self._fileh.readline()
                c = int(line[:8]) - 1
        self._fileh.readline()

    def _check_name(self, name):
        """
        Check name read from op4 file and put '_' on front if needed.

        Returns new name (usually the same as the input name).
        """
        if not (name[0].isalpha() or name[0] == '_'):
            oldname, name = name, '_'+name
            warnings.warn('Output4 file has matrix name: {0}.  '
                          'Changing to {1}.'.format(oldname, name),
                          RuntimeWarning)
        return name

    def _loadop4_ascii(self, patternlist=None, listonly=False):
        """
        Reads next matching matrix or returns information on the next
        matrix in the ascii op4 file.

        Parameters
        ----------
        patternlist : list
            List of string patterns; each matrix name is matched
            against this list:  if it matches any of the patterns, it
            is read in.
        listonly : bool
            True if only reading name.

        Returns
        -------
        tuple: (name, matrix, form, mtype)
            name : string
                Lower-case name of matrix.
            matrix : 2d ndarray
                The matrix.
            form : integer
                Nastran form of matrix.
            mtype : integer
                Nastran matrix type.

        .. note::  All outputs will be None if reached EOF.
        .. note::  The `matrix` output will be [rows, cols] of the
                   matrix if the matrix is skipped.
        """
        while 1:
            line = self._fileh.readline()
            line = line.rstrip()
            if line == '':
                return None, None, None, None
            cols = int(line[:8])
            rows = int(line[8:16])
            form = int(line[16:24])
            mtype = int(line[24:32])
            length = len(line)
            if length > 40:
                name = line[32:40].strip().lower()
            else:
                name = line[32:].lower()
            name = self._check_name(name)
            perline = 5
            numlen = 16
            if length > 44:
                # 1P,3E24.16  <-- starts at position 40
                numformat = line[43:]
                p = numformat.find('E')
                if p < 0:
                    p = numformat.find('D')
                if p > 0:
                    perline = int(numformat[:p])
                    numlen = int(numformat[p+1:].split('.')[0])
            if patternlist and name not in patternlist:
                skip = 1
            else:
                skip = 0
            if listonly or skip:
                self._skipop4_ascii(perline, rows, cols,
                                    mtype, numlen)
                if listonly:
                    return name, (abs(rows), cols), form, mtype
            else:
                break

        if rows < 0 or rows >= self._rows4bigmat:
            rows = abs(rows)
            bigmat = True
        else:
            bigmat = False

        if mtype > 2:
            # complex matrix ... just read as if it's real rows *= 2
            # (must also use fortran ordering for this to work)
            multiplier = 2
        else:
            multiplier = 1

        # real matrix
        X = np.zeros((rows*multiplier, cols), dtype=float, order='F')
        if mtype & 1:
            wper = 1
        else:
            wper = 2

        line = self._fileh.readline()
        linelen = perline * numlen
        c = int(line[:8]) - 1
        r = int(line[8:16])
        if r > 0:
            while c < cols:
                elems = int(line[16:24])
                r = (r-1)*multiplier
                # read column as a long string
                nlines = (elems - 1) // perline
                blocklist = [ln[:linelen]
                             for ln in it.islice(self._fileh,
                                                 nlines)]
                s = ''.join(blocklist) + self._fileh.readline()
                if self._dformat:
                    s = s.replace('D', 'E')
                a = 0
                for i in range(elems):
                    b = a + numlen
                    X[r+i, c] = s[a:b]
                    a = b
                line = self._fileh.readline()
                c = int(line[:8]) - 1
                r = int(line[8:16])
        elif bigmat:
            while c < cols:
                elems = int(line[16:24])
                while elems > 0:
                    line = self._fileh.readline()
                    L = int(line[:8])-1    # L
                    r = int(line[8:16])-1  # irow-1
                    elems -= L + 2
                    r *= multiplier
                    L //= wper
                    nlines = (L - 1) // perline
                    blocklist = [ln[:linelen]
                                 for ln in it.islice(self._fileh,
                                                     nlines)]
                    s = ''.join(blocklist) + self._fileh.readline()
                    if self._dformat:
                        s = s.replace('D', 'E')
                    a = 0
                    for i in range(L):
                        b = a + numlen
                        X[r+i, c] = s[a:b]
                        a = b
                line = self._fileh.readline()
                c = int(line[:8]) - 1
                r = int(line[8:16])
        else:
            while c < cols:
                elems = int(line[16:24])
                while elems > 0:
                    line = self._fileh.readline()
                    IS = int(line)  # [:8])
                    L = (IS >> 16) - 1        # L
                    r = IS - ((L+1) << 16)-1  # irow-1
                    elems -= L + 1
                    r *= multiplier
                    L //= wper
                    nlines = (L - 1) // perline
                    blocklist = [ln[:linelen]
                                 for ln in it.islice(self._fileh,
                                                     nlines)]
                    s = ''.join(blocklist) + self._fileh.readline()
                    if self._dformat:
                        s = s.replace('D', 'E')
                    a = 0
                    for i in range(L):
                        b = a + numlen
                        X[r+i, c] = s[a:b]
                        a = b
                line = self._fileh.readline()
                c = int(line[:8]) - 1
                r = int(line[8:16])
        self._fileh.readline()
        if mtype > 2:
            X.dtype = complex
        return name, X, form, mtype

    def _skipop4_binary(self, cols):
        """
        Skip a binary op4 matrix.

        Parameters
        ----------
        cols : integer
            Number of columns in matrix.
        """
        # Scan matrix by column
        icol = 1
        bi = self._bytes_i
        delta = 4 - bi
        while icol <= cols:
            # Read record length at start of record:
            reclen = self._Str_i4.unpack(self._fileh.read(4))[0]
            # Read column header
            icol = self._Str_i.unpack(self._fileh.read(bi))[0]
            self._fileh.seek(reclen + delta, 1)

    def _loadop4_binary(self, patternlist=None, listonly=False):
        """
        Reads next matching matrix or returns information on the next
        matrix in the binary op4 file.

        Parameters
        ----------
        patternlist : list
            List of string patterns; each matrix name is matched
            against this list:  if it matches any of the patterns, it
            is read in.
        listonly : bool
            True if only reading name.

        Returns
        -------
        tuple: (name, matrix, form, mtype)
            name : string
                Lower-case name of matrix.
            matrix : 2d ndarray
                The matrix.
            form : integer
                Nastran form of matrix.
            mtype : integer
                Nastran matrix type.

        .. note::  All outputs will be None if reached EOF.
        .. note::  The `matrix` output will be [rows, cols] of the matrix
                   if the matrix is skipped.
        """
        fp = self._fileh
        while 1:
            if len(fp.read(4)) == 0:
                return None, None, None, None
            cols, rows, form, mtype =\
                self._Str_iiii.unpack(fp.read(self._bytes_iiii))
            # Read ascii name of matrix:
            if self._bit64:
                name = fp.read(16).decode()
            else:
                name = fp.read(8).decode()
            name = self._check_name(name.lower().strip())
            fp.read(4)
            if patternlist and name not in patternlist:
                skip = 1
            else:
                skip = 0
            if listonly or skip:
                self._skipop4_binary(cols)
                if listonly:
                    return name, (abs(rows), cols), form, mtype
            else:
                break

        if rows < 0 or rows >= self._rows4bigmat:
            rows = abs(rows)
            bigmat = True
        else:
            bigmat = False
        if mtype > 2:
            # complex matrix ... just read as if it's real rows *= 2
            # (must also use fortran ordering for this to work)
            multiplier = 2
        else:
            multiplier = 1

        # real matrix
        X = np.zeros((rows*multiplier, cols), dtype=float, order='F')
        if mtype & 1:
            numform = self._str_sr
            numform2 = self._str_sr_fromfile
            bytesreal = self._bytes_sr
            wper = 1
        else:
            numform = self._str_dr
            numform2 = self._str_dr_fromfile
            bytesreal = 8
            wper = self._wordsperdouble

        reclen = self._Str_i4.unpack(fp.read(4))[0]
        c, r, nwords =\
            self._Str_iii.unpack(fp.read(self._bytes_iii))
        c -= 1
        if r > 0:  # non sparse format
            while c < cols:
                r = (r-1)*multiplier
                nwords //= wper
                if nwords < self._rowsCutoff:
                    X[r:r+nwords, c] =\
                        struct.unpack(numform % nwords,
                                      fp.read(bytesreal*nwords))
                else:
                    X[r:r+nwords, c] = np.fromfile(fp, numform2,
                                                   nwords)
                fp.read(4)
                reclen = self._Str_i4.unpack(fp.read(4))[0]
                c, r, nwords =\
                    self._Str_iii.unpack(fp.read(self._bytes_iii))
                c -= 1
        elif bigmat:
            while c < cols:
                # bigmat sparse format
                # Read column data, one string of numbers at a time
                # (strings of zeros are skipped)
                while nwords > 0:
                    L, r = self._Str_ii.unpack(fp.read(self._bytes_ii))
                    nwords -= L + 1
                    L = (L-1) // wper
                    r = (r-1) * multiplier
                    if L < self._rowsCutoff:
                        X[r:r+L, c] =\
                            struct.unpack(numform % L,
                                          fp.read(bytesreal*L))
                    else:
                        X[r:r+L, c] = np.fromfile(fp, numform2, L)
                fp.read(4)
                reclen = self._Str_i4.unpack(fp.read(4))[0]
                c, r, nwords =\
                    self._Str_iii.unpack(fp.read(self._bytes_iii))
                c -= 1
        else:
            while c < cols:
                # non-bigmat sparse format
                # Read column data, one string of numbers at a time
                # (strings of zeros are skipped)
                while nwords > 0:
                    IS = self._Str_i.unpack(fp.read(self._bytes_i))[0]
                    L = (IS >> 16) - 1             # L
                    r = IS - ((L+1) << 16) - 1     # irow-1
                    nwords -= L + 1                # words left
                    L //= wper
                    r *= multiplier
                    if L < self._rowsCutoff:
                        X[r:r+L, c] =\
                            struct.unpack(numform % L,
                                          fp.read(bytesreal*L))
                    else:
                        X[r:r+L, c] = np.fromfile(fp, numform2, L)
                fp.read(4)
                reclen = self._Str_i4.unpack(fp.read(4))[0]
                c, r, nwords =\
                    self._Str_iii.unpack(fp.read(self._bytes_iii))
                c -= 1
        # read final bytes of record and record marker
        fp.read(reclen-3*self._bytes_i+4)
        if mtype > 2:
            X.dtype = complex
        return name, X, form, mtype

    def _sparse_col_stats(self, v):
        """
        Returns locations of non-zero values and length of each
        series.

        Parameters
        ----------
        v : ndarray
            1d ndarray (the column of the matrix).

        Returns
        -------
        ind : ndarray
            m x 2 ndarray.  m is number of non-zero sequences in v.
            First column contains the indices to the start of each
            sequence and the second column contains the length of
            the sequence.

        For example, if v is:
        ::
          v = [ 0.,  0.,  0.,  7.,  5.,  0.,  6.,  0.,  2.,  3.]

        Then, ind will be:
        ::
          ind = [[3 2]
                 [6 1]
                 [8 2]]
        """
        pv = np.nonzero(v)[0]
        dpv = np.diff(pv)
        starts = np.nonzero(dpv != 1)[0] + 1
        nrows = len(starts)+1
        ind = np.zeros((nrows, 2), int)
        ind[0, 0] = pv[0]
        if nrows > 1:
            ind[1:, 0] = pv[starts]
            ind[0, 1] = starts[0] - 1
            ind[1:-1, 1] = np.diff(starts) - 1
        ind[-1, 1] = len(dpv) - len(starts) - sum(ind[:, 1])
        ind[:, 1] += 1
        return ind

    def _write_ascii_header(self, op4_file, name, matrix, digits, bigmat=False):
        """
        Utility routine that writes the header for ascii matrices.

        Parameters
        ----------
        op4_file : file handle
            Output of open() using binary mode.
        name : string
            Name of matrix.
        matrix : matrix
            Matrix to write.
        digits : integer
            Number of significant digits after the decimal to include
            in the ascii output.
        bigmat : bool
            If true, matrix is to be written in 'bigmat' format.

        Returns
        -------
        tuple: (cols, multiplier, perline, numlen, numform)
            cols : integer
                Number of columns in matrix.
            multiplier : integer
                2 for complex, 1 for real.
            perline : integer
                Number of values written per row.
            numlen : integer
                Number of characters per value.
            numform : string
                Format string for numbers, eg: '%16.9E'.
        """
        numlen = digits + 5 + self._expdigits  # -1.digitsE-009
        perline = 80 // numlen
        rows, cols = matrix.shape
        if rows == cols:
            if np.allclose(matrix.T, matrix):
                form = 6
            else:
                form = 1
        else:
            form = 2
        if np.iscomplexobj(matrix):
            mtype = 4
            multiplier = 2
        else:
            mtype = 2
            multiplier = 1
        if bigmat:
            if rows < self._rows4bigmat:
                rows = -rows
        op4_file.write('{0:8}{1:8}{2:8}{3:8}{4:8s}1P,{5}E{6}.{7}\n'.format(
            cols, rows, form, mtype, name.upper(),
            perline, numlen, digits))
        numform = '%{0}.{1}E'.format(numlen, digits)
        return cols, multiplier, perline, numlen, numform

    def _write_ascii(self, op4_file, name, matrix, digits):
        """
        Write a matrix to a file in ascii, non-sparse format.

        Parameters
        ----------
        op4_file : file handle
            Output of open() using text mode.
        name : string
            Name of matrix.
        matrix : matrix
            Matrix to write.
        digits : integer
            Number of significant digits after the decimal to include
            in the ascii output.
        """
        (cols, multiplier, perline, unused_numlen, numform) = self._write_ascii_header(
             op4_file, name, matrix, digits, bigmat=False)
        for c in range(cols):
            v = matrix[:, c]
            if np.any(v):
                pv = np.nonzero(v)[0]
                s = pv[0, 0]
                e = pv[0, -1]
                elems = (e - s + 1) * multiplier
                op4_file.write('{0:8}{1:8}{2:8}\n'.format(c+1, s+1, elems))
                v = np.asarray(v[s:e+1]).flatten()
                v.dtype = float
                neven = ((elems - 1) // perline) * perline
                for i in range(0, neven, perline):
                    for j in range(perline):
                        op4_file.write(numform % v[i+j])
                    op4_file.write('\n')
                for i in range(neven, elems):
                    op4_file.write(numform % v[i])
                op4_file.write('\n')
        op4_file.write('{0:8}{1:8}{2:8}\n'.format(cols+1, 1, 1))
        op4_file.write(numform % 2**.5)
        op4_file.write('\n')

    def _write_ascii_sparse_nonbigmat(self, op4_file, name, matrix, digits):
        """
        Write a matrix to a file in ascii, non-bigmat sparse format.

        Parameters
        ----------
        op4_file : file handle
            Output of open() using binary mode.
        name : string
            Name of matrix.
        matrix : matrix
            Matrix to write.
        digits : integer
            Number of significant digits after the decimal to include
            in the ascii output.

        .. note:: if rows > 65535, bigmat is turned on and the
                  :func:`write_ascii_sparse_bigmat` function is
                  called ...that's a Nastran rule.
        """
        rows, cols = matrix.shape
        if rows >= self._rows4bigmat:
            self._write_ascii_sparse_bigmat(op4_file, name, matrix, digits)
            return
        (cols, multiplier, perline,  numform) = self._write_ascii_header(
             op4_file, name, matrix, digits, bigmat=False)
        for c in range(cols):
            v = matrix[:, c]
            if np.any(v):
                v = np.asarray(v).flatten()
                ind = self._sparse_col_stats(v)
                nwords = ind.shape[0] + 2*sum(ind[:, 1])*multiplier
                op4_file.write('{0:8}{1:8}{2:8}\n'.format(c+1, 0, nwords))
                for row in ind:
                    r = row[0]
                    L = row[1]*2*multiplier
                    IS = (r+1) + ((L+1) << 16)
                    op4_file.write('{0:12}\n'.format(IS))
                    string = v[r:r+row[1]]
                    string.dtype = float
                    elems = L // 2
                    neven = ((elems - 1) // perline) * perline
                    for i in range(0, neven, perline):
                        for j in range(perline):
                            op4_file.write(numform % string[i+j])
                        op4_file.write('\n')
                    for i in range(neven, elems):
                        op4_file.write(numform % string[i])
                    op4_file.write('\n')
        op4_file.write('{0:8}{1:8}{2:8}\n'.format(cols+1, 1, 1))
        op4_file.write(numform % 2**.5)
        op4_file.write('\n')

    def _write_ascii_sparse_bigmat(self, op4_file, name, matrix, digits):
        """
        Write a matrix to a file in ascii, bigmat sparse format.

        Parameters
        ----------
        op4_file : file handle
            Output of open() using binary mode.
        name : string
            Name of matrix.
        matrix : matrix
            Matrix to write.
        digits : integer
            Number of significant digits after the decimal to include
            in the ascii output.
        """
        (cols, multiplier, perline, unused_numlen, numform) = self._write_ascii_header(
             op4_file, name, matrix, digits, bigmat=True)
        for c in range(cols):
            v = matrix[:, c]
            if np.any(v):
                v = np.asarray(v).flatten()
                ind = self._sparse_col_stats(v)
                nwords = 2*ind.shape[0] + 2*sum(ind[:, 1])*multiplier
                op4_file.write('{0:8}{1:8}{2:8}\n'.format(c+1, 0, nwords))
                for row in ind:
                    r = row[0]
                    L = row[1]*2*multiplier
                    op4_file.write('{0:8}{1:8}\n'.format(L+1, r+1))
                    string = v[r:r+row[1]]
                    string.dtype = float
                    elems = L // 2
                    neven = ((elems - 1) // perline) * perline
                    for i in range(0, neven, perline):
                        for j in range(perline):
                            op4_file.write(numform % string[i+j])
                        op4_file.write('\n')
                    for i in range(neven, elems):
                        op4_file.write(numform % string[i])
                    op4_file.write('\n')
        op4_file.write('{0:8}{1:8}{2:8}\n'.format(cols+1, 1, 1))
        op4_file.write(numform % 2**.5)
        op4_file.write('\n')

    def _write_binary_header(self, op4_file, name, matrix,
                             endian, bigmat=False):
        """
        Utility routine that writes the header for binary matrices.

        Parameters
        ----------
        op4_file : file handle
            Output of open() using binary mode.
        name : string
            Name of matrix.
        matrix : matrix
            Matrix to write.
        endian : string
            Endian setting for binary output:  '' for native, '>' for
            big-endian and '<' for little-endian.
        bigmat : bool
            If true, matrix is to be written in 'bigmat' format.

        Returns
        -------
        tuple: (cols, multiplier)
            cols : integer
                Number of columns in matrix.
            multiplier : integer
                2 for complex, 1 for real.
        """
        rows, cols = matrix.shape
        if rows == cols:
            if np.allclose(matrix.T, matrix):
                form = 6
            else:
                form = 1
        else:
            form = 2
        if np.iscomplexobj(matrix):
            mtype = 4
            multiplier = 2
        else:
            mtype = 2
            multiplier = 1

        # write 1st record (24 bytes: 4 4-byte ints, 1 8-byte string)
        name = ('{0:<8}'.format(name.upper())).encode()
        if bigmat:
            if rows < self._rows4bigmat:
                rows = -rows
        op4_file.write(struct.pack(endian+'5i8si', 24, cols, rows,
                                   form, mtype, name, 24))
        return cols, multiplier

    def _write_binary(self, op4_file, name, matrix, endian):
        """
        Write a matrix to a file in double precision binary format.

        Parameters
        ----------
        op4_file : file handle
            Output of open() using binary mode.
        name : string
            Name of matrix.
        matrix : matrix
            Matrix to write.
        endian : string
            Endian setting for binary output:  '' for native, '>' for
            big-endian and '<' for little-endian.
        """
        cols, multiplier = self._write_binary_header(
            op4_file, name, matrix, endian)
        col_header = struct.Struct(endian+'4i')
        col_trailer = struct.Struct(endian+'i')
        for c in range(cols):
            v = matrix[:, c]
            if np.any(v):
                pv = np.nonzero(v)[0]
                s = pv[0, 0]
                e = pv[0, -1]
                elems = (e - s + 1) * multiplier
                reclen = 3*4 + elems*8
                op4_file.write(col_header.pack(reclen, c+1, s+1, 2*elems))
                v = np.asarray(v[s:e+1]).flatten()
                v.dtype = float
                op4_file.write(struct.pack(endian+('%dd' % elems), *v))
                op4_file.write(col_trailer.pack(reclen))
        reclen = 3*4 + 8
        op4_file.write(col_header.pack(reclen, cols+1, 1, 2))
        op4_file.write(struct.pack(endian+'d', 2**.5))
        op4_file.write(col_trailer.pack(reclen))

    def _write_binary_sparse_nonbigmat(self, op4_file, name, matrix, endian):
        """
        Write a matrix to a file in double precision binary, non-bigmat
        sparse format.

        Parameters
        ----------
        op4_file : file handle
            Output of open() using binary mode.
        name : string
            Name of matrix.
        matrix : matrix
            Matrix to write.
        endian : string
            Endian setting for binary output:  '' for native, '>' for
            big-endian and '<' for little-endian.

        .. note::  if rows > 65535, bigmat is turned on and the
                   :func:`write_binary_sparse_bigmat` function is
                   called ...that's a Nastran rule.
        """
        rows, cols = matrix.shape
        if rows >= self._rows4bigmat:
            self._write_binary_sparse_bigmat(op4_file, name, matrix, endian)
            return
        cols, multiplier = self._write_binary_header(
            op4_file, name, matrix, endian)
        col_header = struct.Struct(endian+'4i')
        col_trailer = struct.Struct(endian+'i')
        for c in range(cols):
            v = matrix[:, c]
            if np.any(v):
                v = np.asarray(v).flatten()
                ind = self._sparse_col_stats(v)
                nwords = ind.shape[0] + 2*sum(ind[:, 1])*multiplier
                reclen = (3 + nwords)*4
                op4_file.write(col_header.pack(reclen, c+1, 0, nwords))
                for row in ind:
                    r = row[0]
                    L = row[1]*2*multiplier
                    IS = (r+1) + ((L+1) << 16)
                    op4_file.write(col_trailer.pack(IS))
                    string = v[r:r+row[1]]
                    string.dtype = float
                    op4_file.write(struct.pack(endian+('%dd' % len(string)), *string))
                op4_file.write(col_trailer.pack(reclen))
        reclen = 3*4 + 8
        op4_file.write(col_header.pack(reclen, cols+1, 1, 2))
        op4_file.write(struct.pack(endian+'d', 2**.5))
        op4_file.write(col_trailer.pack(reclen))

    def _write_binary_sparse_bigmat(self, op4_file, name, matrix, endian):
        """
        Write a matrix to a file in double precision binary, bigmat
        sparse format.

        Parameters
        ----------
        op4_file : file handle
            Output of open() using binary mode.
        name : string
            Name of matrix.
        matrix : matrix
            Matrix to write.
        endian : string
            Endian setting for binary output:  '' for native, '>' for
            big-endian and '<' for little-endian.
        """
        cols, multiplier = self._write_binary_header(
            op4_file, name, matrix, endian, True)
        colHeader = struct.Struct(endian+'4i')
        colTrailer = struct.Struct(endian+'i')
        LrStruct = struct.Struct(endian+'ii')
        for c in range(cols):
            v = matrix[:, c]
            if np.any(v):
                v = np.asarray(v).flatten()
                ind = self._sparse_col_stats(v)
                nwords = 2*ind.shape[0] + 2*sum(ind[:, 1])*multiplier
                reclen = (3 + nwords)*4
                op4_file.write(colHeader.pack(reclen, c+1, 0, nwords))
                for row in ind:
                    r = row[0]
                    L = row[1]*2*multiplier
                    op4_file.write(LrStruct.pack(L+1, r+1))
                    string = v[r:r+row[1]]
                    string.dtype = float
                    op4_file.write(struct.pack(endian+('%dd' % len(string)), *string))
                op4_file.write(colTrailer.pack(reclen))
        reclen = 3*4 + 8
        op4_file.write(colHeader.pack(reclen, cols+1, 1, 2))
        op4_file.write(struct.pack(endian+'d', 2**.5))
        op4_file.write(colTrailer.pack(reclen))

    def dctload(self, filename, namelist=None):
        """
        Read all matching matrices from op4 file into dictionary.

        Parameters
        ----------
        filename : string
            Name of op4 file to read.
        namelist : list, string, or None
            List of variable names to read in, or string with name of
            the single variable to read in, or None.  If None, all
            matrices are read in.

        Returns
        -------
        dct : dictionary
            Keys are the lower-case matrix names and the values are a
            tuple of:  (matrix, form, mtype).

        See also :func:`listload`, :func:`write`, :func:`dir`.
        """
        if isinstance(namelist, str):
            namelist = [namelist]
        self._op4open_read(filename)
        dct = {}
        try:
            if self._ascii:
                loadfunc = self._loadop4_ascii
            else:
                loadfunc = self._loadop4_binary
            while 1:
                name, X, form, mtype =\
                    loadfunc(patternlist=namelist)
                if not name:
                    break
                dct[name] = X, form, mtype
        finally:
            self._op4close()
        return dct

    def listload(self, filename, namelist=None):
        """
        Read all matching matrices from op4 file into a list; useful
        if op4 file has duplicate names.

        Parameters
        ----------
        filename : string
            Name of op4 file to read.
        namelist : list, string, or None
            List of variable names to read in, or string with name of
            the single variable to read in, or None.  If None, all
            matrices are read in.

        Returns
        -------
        tuple: (names, matrices, forms, mtypes)
            names : list
                Lower-case list of matrix names in order as read.
            matrices : list
                List of matrices in order as read.
            forms : list
                List of integers specifying the Nastran form of each
                matrix.
            mtypes : list
                List of integers specifying the Nastran type of each
                matrix.

        See also :func:`dctload`, :func:`write`, :func:`dir`.
        """
        if isinstance(namelist, str):
            namelist = [namelist]
        self._op4open_read(filename)
        names = []
        matrices = []
        forms = []
        mtypes = []
        try:
            if self._ascii:
                loadfunc = self._loadop4_ascii
            else:
                loadfunc = self._loadop4_binary
            while 1:
                name, X, form, mtype =\
                    loadfunc(patternlist=namelist)
                if not name:
                    break
                names.append(name)
                matrices.append(X)
                forms.append(form)
                mtypes.append(mtype)
        finally:
            self._op4close()
        return names, matrices, forms, mtypes

    def dir(self, filename, verbose=True):
        """
        Directory of all matrices in op4 file.

        Parameters
        ----------
        filename : string
            Name of op4 file to read.
        verbose : bool
            If true, directory will be printed to screen.

        Returns
        -------
        tuple: (names, sizes, forms, mtypes)
            names : list
                Lower-case list of matrix names in order as read.
            sizes : list
                List of sizes [(r1, c1), (r2, c2), ...], for each
                matrix.
            forms : list
                List of integers specifying the Nastran form of each
                matrix.
            mtypes : list
                List of integers specifying the Nastran type of each
                matrix.

        See also :func:`dctload`, :func:`listload`, :func:`write`.
        """
        self._op4open_read(filename)
        names = []
        sizes = []
        forms = []
        mtypes = []
        try:
            if self._ascii:
                loadfunc = self._loadop4_ascii
            else:
                loadfunc = self._loadop4_binary
            while 1:
                name, X, form, mtype =\
                    loadfunc(listonly=True)
                if not name:
                    break
                names.append(name)
                sizes.append(X)
                forms.append(form)
                mtypes.append(mtype)
            if verbose:
                for n, s, f, m in zip(names, sizes, forms, mtypes):
                    print('{0:8}, {1:6} x {2:<6}, form={3}, mtype={4}'
                          .format(n, s[0], s[1], f, m))
        finally:
            self._op4close()
        return names, sizes, forms, mtypes

    def write(self, filename, names, matrices=None,
              binary=True, digits=16, endian='',
              sparse=''):
        """
        Write op4 file.

        Parameters
        ----------
        filename : string
            Name of file.
        names : string or list or dictionary
            Matrix name or list of matrix names or dictionary indexed
            by the names.
        matrices : array or list
            2d ndarray or list of 2d ndarrays.  Ignored if `names` is
            a dictionary.
        binary : bool
            If true, a double precision binary file is written;
            otherwise an ascii file is created.
        digits : integer
            Number of significant digits after the decimal to include
            in the ascii output.  Ignored for binary files.
        endian : string
            Endian setting for binary output:  '' for native, '>' for
            big-endian and '<' for little-endian.
        sparse : string
            Empty or 'bigmat' or 'nonbigmat'.  If set to 'bigmat' or
            'nonbigmat', that sparse format is selected.  Note that if
            the number of rows is > 65535, then both the 'bigmat' and
            'nonbigmat' options become 'bigmat'.

        Returns
        -------
        None.

        .. note::  To write multiple matrices that have the same name
                   or to write the matrices in a specific order, `names`
                   must be a list, not a dictionary.  If a dictionary,
                   the matrices are written in alphabetical order.

        See also :func:`dctload`, :func:`listload`, :func:`dir`.

        Examples
        --------
        To write m, k, b, in that order to an ascii file::

            import numpy as np
            import op4
            o4 = op4.OP4()
            m = np.array([1, 2])
            k = np.array([3, 5])
            b = np.array([4, 6])
            names = ['m', 'k', 'b']
            values = [eval(v) for v in names]
            o4.write('mkb.op4', names, values, False, 9)

        Or, if you don't care about the order, you could create the
        dictionary input:
        ::
            o4.write('mkb.op4', dict(m=m, k=k, b=b),
                     binary=False, digits=9)
        """
        if isinstance(names, dict):
            k = sorted(names.keys())
            matrices = [names[j] for j in k]
            names = k
        else:
            if not isinstance(names, list):
                names = [names]
            if not isinstance(matrices, list):
                matrices = [matrices]

        # ensure double precision 2d arrays:
        def ensure_2d_dp(m):
            """Ensures 2d double precision array"""
            m = np.asmatrix(m)
            if np.iscomplexobj(m):
                if m.dtype != np.complex128:
                    return m.astype(np.complex128)
            elif m.dtype != np.float64:
                return m.astype(np.float64)
            return m
        matrices = [ensure_2d_dp(m) for m in matrices]

        if binary:
            if sparse == '':
                wrtfunc = self._write_binary
            elif sparse == 'bigmat':
                wrtfunc = self._write_binary_sparse_bigmat
            elif sparse == 'nonbigmat':
                wrtfunc = self._write_binary_sparse_nonbigmat
            else:
                raise ValueError('invalid sparse option')
            with open(filename, 'wb') as op4_file:
                for name, matrix in zip(names, matrices):
                    wrtfunc(op4_file, name, matrix, endian)
        else:
            if sparse == '':
                wrtfunc = self._write_ascii
            elif sparse == 'bigmat':
                wrtfunc = self._write_ascii_sparse_bigmat
            elif sparse == 'nonbigmat':
                wrtfunc = self._write_ascii_sparse_nonbigmat
            else:
                raise ValueError('invalid sparse option')
            with open(filename, 'w') as op4_file:
                for name, matrix in zip(names, matrices):
                    wrtfunc(op4_file, name, matrix, digits)
