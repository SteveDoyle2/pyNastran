# -*- coding: utf-8 -*-
"""
Some Python tools for reading select data from Nastran .op2 files.
Converted from the Yeti version.

Can read files in big or little endian format.

@author: Tim Widrick
"""
import sys
import struct
import itertools as it
import warnings

import numpy as np

import pyNastran.op2.dev.pyyeti.n2y as n2y

#  Notes on the op2 format.
#
#  DATA BLOCK:
#      All data blocks (including header) start with header 3 elements:
#      [reclen, key, endrec]
#        - reclen = 1 32-bit integer that specifies number of bytes in
#          key (either 4 or 8)
#        - key = 4 or 8 byte integer specifying number of words in next
#          record
#        - endrec = reclen
#
#      DATA SET, can be multiple records:
#          Next is [reclen, data, endrec]
#            - reclen = 1 32-bit integer that specifies number of bytes
#              in data
#            - data = reclen bytes long, variable format; may be part of
#              a data set or the complete set
#            - endrec = reclen
#
#          Next is info about whether we're done with current data set:
#          [reclen, key, endrec]
#            - reclen = 1 32-bit integer that specifies number of bytes
#              in key (either 4 or 8)
#            - key = 4 or 8 byte integer specifying number of words in
#              next record; if 0, done with data set
#            - endrec = reclen
#
#          If not done, we have [reclen, data, endrec] for part 2 (and
#          so on) for the record.
#
#      Once data set is complete, we have: [reclen, key, endrec]
#        - reclen = 1 32-bit integer that specifies number of bytes in
#          key (either 4 or 8)
#        - key = 4 or 8 byte integer specifying number of words in next
#          record (I think ... not useful?)
#        - endrec = reclen
#
#      Then: [reclen, rec_type, endrec]
#        - reclen = 1 32-bit integer that specifies number of bytes in
#          rec_type (either 4 or 8)
#        - rec_type = 0 if table (4 or 8 bytes)
#        - endrec = reclen
#
#      Then, info on whether we're done with data block:
#      [reclen, key, endrec]
#        - reclen = 1 32-bit integer that specifies number of bytes in
#          key (either 4 or 8)
#        - key = 4 or 8 byte integer specifying number of words in next
#          record; if 0, done with data block
#        - endrec = reclen
#
#      If not done, we have [reclen, data, endrec] for record 2 and so
#      on, until data block is read in.


def expand_dof(ids, pvgrids):
    """
    Expands vector of ids to [id, dof].

    Parameters
    ----------
    ids : 1d array-like
        Vector of node ids
    pvgrids : 1d array-like
        True/False vector same length as `ids`. The True entries
        indicate which elements in `ids` are grids; those will get all
        6 DOF while all other ids will just get 0 for the DOF.

    Returns
    -------
    dof : 2d ndarray
        2 column matrix: [id, dof]

    Examples
    --------
    >>> import numpy as np
    >>> import op2
    >>> ids = [1, 2, 3, 4]
    >>> pvgrids = [True, False, False, True]
    >>> expand_dof(ids, pvgrids)
    array([[1, 1],
           [1, 2],
           [1, 3],
           [1, 4],
           [1, 5],
           [1, 6],
           [2, 0],
           [3, 0],
           [4, 1],
           [4, 2],
           [4, 3],
           [4, 4],
           [4, 5],
           [4, 6]])

    """
    ids, pvgrids = np.atleast_1d(ids, pvgrids)
    n = len(ids)
    dof = np.zeros((n, 6), int)
    dof[pvgrids] = np.arange(1, 7)
    V = np.zeros((n, 6), bool)
    V[:, 0] = True
    V[pvgrids, 1:] = True
    expids = np.reshape(ids, (-1, 1)) * V
    V = V.flatten()
    expids = expids.flatten()
    dof = dof.flatten()
    return np.vstack((expids[V], dof[V])).T


class OP2:
    """Class for reading Nastran op2 files and nas2cam data files."""

    def __init__(self, filename=None):
        self._fileh = None
        self._CodeFuncs = None
        if isinstance(filename, str):
            self._op2_open(filename)

    def __del__(self):
        if self._fileh:
            self._fileh.close()
            self._fileh = None

    def __enter__(self):
        return self

    def __exit__(self, type, value, traceback):
        if self._fileh:
            self._fileh.close()
            self._fileh = None
        return False

    @property
    def CodeFuncs(self):
        """See :func:`_check_code`."""
        if self._CodeFuncs is None:
            def func1(item_code):
                if item_code // 1000 in [2, 3, 6]:
                    return 2
                return 1

            def func2(item_code):
                return item_code % 100

            def func3(item_code):
                return item_code % 1000

            def func4(item_code):
                return item_code // 10

            def func5(item_code):
                return item_code % 10

            def func6(item_code):
                warnings.warn('Function code 6 method not verified',
                              RuntimeWarning)
                if item_code & 8:
                    return 0
                return 1

            def func7(item_code):
                v = item_code // 1000
                if v in [0, 2]:
                    return 0
                if v in [1, 3]:
                    return 1
                return 2

            def funcbig(func_code, item_code):
                return item_code & (func_code & 65535)

            self._CodeFuncs = {
                1: func1, 2: func2, 3: func3, 4: func4,
                5: func5, 6: func6, 7: func7,
                'big': funcbig,
            }
        return self._CodeFuncs

    def _op2_open(self, filename):
        """
        Open op2 file in correct endian mode.

        Sets these class variables:

        _fileh : file handle
            Value returned by open().
        _swap : bool
            True if bytes must be swapped to correct endianness.
        _bit64 : True or False
            True if 'key' integers are 64-bit.
        _endian : string
            Will be '=' if `swap` is False; otherwise, either '>' or '<'
            for big-endian and little-endian, respectively.
        _intstr : string
            Either `endian` + 'i4' or `endian` + 'i8'.
        _ibytes : integer
            Either 4 or 8 (corresponds to `intstr`)
        _int32str : string
           `endian` + 'i4'.
        _label : string
            The op2 header label or, if none, None.
        _date : vector
            Three element date vector, or None.
        _nastheader : string
            Nastran header for file, or None.
        _postheaderpos : integer
            File position after header.
        dbnames : dictionary
            See :func:`directory` for description.  Contains data block
            names, bytes in file, file positions, and for matrices, the
            matrix size.
        dblist : list
            See :func:`directory` for description.  Contains same info
            as dbnames, but in a list of ordered and formatted strings.
        _Str4 : struct.Struct object
            Precompiled for reading 4 byte integers (corresponds to
            `int32str`).
        _Str : struct.Struct object
            Precompiled for reading 4 or 8 byte integers (corresponds
            to `intstr`).

        File is positioned after the header label (at `postheaderpos`).
        """
        self._fileh = open(filename, 'rb')
        self.dbnames = []
        self.dblist = []
        reclen = struct.unpack('i', self._fileh.read(4))[0]
        self._fileh.seek(0)

        reclen = np.array(reclen, dtype=np.int32)
        if not np.any(reclen == [4, 8]):
            self._swap = True
            reclen = reclen.byteswap()
            if not np.any(reclen == [4, 8]):
                self._fileh.close()
                self._fileh = None
                raise RuntimeError('Could not decipher file.  First'
                                   '4-byte integer should be 4 or 8.')
            if sys.byteorder == 'little':
                self._endian = b'>'
            else:
                self._endian = b'<'
        else:
            self._swap = False
            self._endian = b'='

        self._Str4 = struct.Struct(self._endian + b'i')
        if reclen == 4:
            self._bit64 = False
            self._intstr = self._endian + b'i4'
            self._intstru = self._endian + b'%di'
            self._ibytes = 4
            self._Str = self._Str4
        else:
            self._bit64 = True
            self._intstr = self._endian + b'i8'
            self._intstru = self._endian + b'%dq'
            self._ibytes = 8
            self._Str = struct.Struct(self._endian + b'q')
        #print('bit64 = ', self._bit64)

        self._rowsCutoff = 3000
        self._int32str = self._endian + b'i4'
        self._int32stru = self._endian + b'%di'
        self._read_op2_header()
        self._postheaderpos = self._fileh.tell()
        self.directory(verbose=False)

    def _get_key(self):
        """Reads [reclen, key, endrec] triplet and returns key."""
        self._fileh.read(4)
        key = self._Str.unpack(self._fileh.read(self._ibytes))[0]
        self._fileh.read(4)
        return key

    def _skip_key(self, n):
        """Skips `n` key triplets ([reclen, key, endrec])."""
        self._fileh.read(n*(8+self._ibytes))

    def _read_op2_header(self):
        """
        Returns Nastran output2 header label (or 'no header').
        """
        key = self._get_key()
        if key != 3:
            self._fileh.seek(0)
            self._date = self._nastheader = self._label = None
            return

        self._fileh.read(4)  # reclen
        frm = self._intstru % key
        kbytes = self._ibytes * key
        self._date = struct.unpack(frm, self._fileh.read(kbytes))
        # self._date = np.fromfile(self._fileh, self._intstr, key)
        self._fileh.read(4)  # endrec
        self._get_key()

        reclen = self._Str4.unpack(self._fileh.read(4))[0]
        self._nastheader = self._fileh.read(reclen).decode()
        self._fileh.read(4)  # endrec
        self._get_key()

        reclen = self._Str4.unpack(self._fileh.read(4))[0]
        self._label = self._fileh.read(reclen).decode().strip().replace(' ', '')
        self._fileh.read(4)  # endrec
        self._skip_key(2)

    def _valid_name(self, bstr):
        """
        Returns a valid variable name from the byte string `bstr`.
        """
        return ''.join(chr(c) for c in bstr if (
            47 < c < 58 or 64 < c < 91 or c == 95 or 96 < c < 123))

    def _read_op2_end_of_table(self):
        """Read Nastran output2 end-of-table marker.

        Returns
        -------
        tuple: (eot, key)
            eot : integer
                1 if end-of-file has been reached and 0 otherwise.
            key : integer
                0 of eot is 1; next key value otherwise.
        """
        bstr = self._fileh.read(4)  # reclen
        if len(bstr) == 4:
            key = self._Str.unpack(self._fileh.read(self._ibytes))[0]
            self._fileh.read(4)  # endrec
        else:
            key = 0
        if key == 0:
            return 1, 0
        return 0, key

    def _read_op2_name_trailer(self):
        """Read Nastran output2 datablock name and trailer.

        Returns
        -------
        tuple: (name, trailer, type)
            name : string
                Name of upcoming data block (upper case).
            trailer : tuple
                Data block trailer.
            type : 0 or 1
                0 means table, 1 means matrix.  I think.

        All outputs will be None for end-of-file.
        """
        eot, key = self._read_op2_end_of_table()
        if key == 0:
            #print('return None, None, None')
            return None, None, None

        reclen = self._Str4.unpack(self._fileh.read(4))[0]
        db_binary_name = self._fileh.read(reclen)
        db_name = self._valid_name(db_binary_name)
        self._fileh.read(4)  # endrec
        self._get_key()
        key = self._get_key()

        self._fileh.read(4)  # reclen
        frm = self._intstru % key
        nbytes = self._ibytes * key

        # prevents a giant read
        assert nbytes > 0, nbytes
        trailer = struct.unpack(frm, self._fileh.read(nbytes))
        # trailer = np.fromfile(self._fileh, self._intstr, key)
        self._fileh.read(4)  # endrec
        self._skip_key(4)

        reclen = self._Str4.unpack(self._fileh.read(4))[0]
        unused_db_name2 = self._valid_name(self._fileh.read(reclen))
        self._fileh.read(4)  # endrec

        self._skip_key(2)
        rec_type = self._get_key()
        return db_name, trailer, rec_type

    def read_op2_matrix(self, name, trailer):
        """
        Read and return Nastran op2 matrix at current file position.

        It is assumed that the name has already been read in via
        :func:`_read_op2_name_trailer`.

        The size of the matrix is read from trailer:
             nrows = trailer[2]
             ncols = trailer[1]
        """
        dtype = 1
        nrows = trailer[2]
        ncols = trailer[1]
        print('    %s (%s, %s)' % (name, nrows, ncols))
        matrix = np.zeros((nrows, ncols), order='F')
        if self._bit64:
            intsize = 8
        else:
            intsize = 4
        col = 0
        frm = self._endian + b'%dd'
        #print('frm =', frm)
        while dtype > 0:  # read in matrix columns
            # key is number of elements in next record (row # followed
            # by key-1 real numbers)
            key = self._get_key()
            # read column
            while key > 0:
                reclen = self._Str4.unpack(self._fileh.read(4))[0]
                r = self._Str.unpack(self._fileh.read(self._ibytes))[0]-1
                n = (reclen - intsize) // 8
                if n < self._rowsCutoff:
                    matrix[r:r+n, col] = struct.unpack(
                        frm % n, self._fileh.read(n*8))
                else:
                    matrix[r:r+n, col] = np.fromfile(
                        self._fileh, np.float64, n)
                self._fileh.read(4)  # endrec
                key = self._get_key()
            col += 1
            self._get_key()
            dtype = self._get_key()
        self._read_op2_end_of_table()
        if self._swap:
            matrix = matrix.byteswap()

        if name in ['EFMFSMS', 'EFMASSS', 'RBMASSS']:
            print(matrix)
        return matrix

    def skip_op2_matrix(self, trailer):
        """
        Skip Nastran op2 matrix at current position.

        It is assumed that the name has already been read in via
        :func:`_read_op2_name_trailer`.

        The size of the matrix is read from trailer:
             rows = trailer[2]
             cols = trailer[1]
        """
        dtype = 1
        while dtype > 0:  # read in matrix columns
            # key is number of elements in next record (row # followed
            # by key-1 real numbers)
            key = self._get_key()
            # skip column
            while key > 0:
                reclen = self._Str4.unpack(self._fileh.read(4))[0]
                self._fileh.seek(reclen, 1)
                self._fileh.read(4)  # endrec
                key = self._get_key()
            self._get_key()
            dtype = self._get_key()
        self._read_op2_end_of_table()

    def skip_op2_table(self):
        """Skip over Nastran output2 table."""
        eot, key = self._read_op2_end_of_table()
        if key == 0:
            return
        while key > 0:
            while key > 0:
                reclen = self._Str4.unpack(self._fileh.read(4))[0]
                self._fileh.seek(8+reclen, 1)
                key = self._Str.unpack(self._fileh.read(self._ibytes))[0]
                self._fileh.read(4)  # endrec
            self._skip_key(2)
            eot, key = self._read_op2_end_of_table()

    def read_op2_matrices(self):
        """Read all matrices from Nastran output2 file.

        Returns dictionary containing all matrices in the op2 file:
        {'NAME1': matrix1, 'NAME2': matrix2, ...}

        The keys are the names as stored (upper case).
        """
        self._fileh.seek(self._postheaderpos)
        mats = {}
        while 1:
            name, trailer, rectype = self._read_op2_name_trailer()
            if name is None:
                break
            if rectype > 0:
                print("Reading matrix {}...".format(name))
                mats[name] = self.read_op2_matrix(trailer)
            else:
                self.skip_op2_table()
        return mats

    def print_data_block_directory(self):
        """
        Prints op2 data block directory.  See also :func:`directory`.
        """
        if len(self.dblist) == 0:
            self.directory(verbose=False)
        for s in self.dblist:
            print(s)

    def directory(self, verbose=True, redo=False): # TODO: _read_op2_name_trailer
        """
        Return list of data block names in op2 file.

        Parameters
        ----------
        verbose : bool (or any true/false variable)
            If True, print names, sizes, and file offsets to screen.
        redo : bool
            If True, scan through file and redefine self.dbnames even
            if it is already set.

        Returns tuple: (dbnames, dblist)
        --------------------------------
        dbnames : Dictionary
            Dictionary indexed by data block name.  Each value is a
            list, one element per occurrence of the data block in the
            op2 file.  Each element is another list that has 3
            elements: [fpos, bytes, size]:
            ::
               fpos : 2-element list; file position start and stop
                      (stop value is start of next data block)
               bytes: number of bytes data block consumes in file
               size : 2-element list; for matrices, [rows, cols],
                      for tables [0, 0]
        dblist : list
            List of strings for printing.  Contains the info above
            in formatted and sorted (in file position order) strings.

        As an example of using dbnames, to get a list of all sizes of
        matrices named 'KAA':
        ::
            o2 = op2.OP2('mds.op2')
            s = [item[2] for item in o2.dbnames['KAA']]

        For another example, to read in first matrix named 'KAA':
        ::
            o2 = op2.OP2('mds.op2')
            fpos = o2.dbnames['KAA'][0][0][0]
            o2._fileh.seek(fpos)
            name, trailer, rectype = o2._read_op2_name_trailer()
            kaa = o2.read_op2_matrix(trailer)

        This routine also sets self.dbnames = dbnames.
        """
        if len(self.dbnames) > 0 and not redo:
            return self.dbnames
        dbnames = {}
        dblist = []
        self._fileh.seek(self._postheaderpos)
        pos = self._postheaderpos
        while 1:
            name, trailer, dbtype = self._read_op2_name_trailer()
            if name is None:
                break
            if dbtype > 0:
                self.skip_op2_matrix(trailer)
                size = [trailer[2], trailer[1]]
                s = 'Matrix {0:8}'.format(name)
            else:
                self.skip_op2_table()
                size = [0, 0]
                s = 'Table  {0:8}'.format(name)
            cur = self._fileh.tell()
            s += (', bytes = {0:10} [{1:10} to {2:10}]'.
                  format(cur-pos-1, pos, cur))
            if size != [0, 0]:
                s += (', {0:6} x {1:<}'.
                      format(size[0], size[1]))
            if name not in dbnames:
                dbnames[name] = []
            dbnames[name].append([[pos, cur], cur-pos-1, size])
            dblist.append(s)
            pos = cur
        self.dbnames = dbnames
        self.dblist = dblist
        if verbose:
            self.print_data_block_directory()
        return dbnames, dblist

    def read_op2_dynamics(self):
        """
        Reads the TLOAD data from a DYNAMICS datablock.

        Returns matrix of TLOADS.  Rows = 5 or 6, Cols = number of
        TLOADs.  TLOAD ids are in first row; other data in matrix may
        not be useful.
        """
        key = self._get_key()
        if self._ibytes == 4:
            header_str = struct.Struct(self._endian + b'iii')
            hbytes = 12
        else:
            header_str = struct.Struct(self._endian + b'qqq')
            hbytes = 24

        eot = 0
        print('self._intstr = %r' % self._intstr)
        data = np.zeros(0, dtype=self._intstr)
        while not eot:
            while key > 0:
                self._fileh.read(4)  # reclen
                header = header_str.unpack(self._fileh.read(hbytes))
                if header == (7107, 71, 138):
                    if key < self._rowsCutoff:
                        bytes = (key-3)*self._ibytes
                        ndata = struct.unpack(self._intstru % (key-3),
                                              self._fileh.read(bytes))
                    else:
                        ndata = np.fromfile(self._fileh,
                                            self._intstr, key-3)
                    data = np.hstack((data, ndata))
                else:
                    self._fileh.seek((key-3)*self._ibytes, 1)
                self._fileh.read(4)  # endrec
                key = self._get_key()
            self._skip_key(2)
            eot, key = self._read_op2_end_of_table()

        if np.any(data):
            L = len(data)
            mult5 = L == 5*(L // 5)
            mult6 = L == 6*(L // 6)
            err1 = ('Could not determine if TLOADs are 5 or 6 rows!  '
                    'Both work.  Routine needs updating.')
            err2 = ('Could not determine if TLOADs are 5 or 6 rows!  '
                    'Neither work.  Routine needs updating.')
            if mult5:
                mindelta5 = np.min(np.diff(data[0::5]))
            if mult6:
                mindelta6 = np.min(np.diff(data[0::6]))
            if mult5:
                if mult6:
                    # L is multiple of both 5 and 6:
                    if mindelta5 > 0:
                        if mindelta6 > 0:
                            raise ValueError(err1)
                        rows = 5
                    else:
                        if mindelta6 > 0:
                            rows = 6
                        else:
                            raise ValueError(err2)
                else:
                    if mindelta5 > 0:
                        rows = 5
                    else:
                        raise ValueError(err2)
            elif mult6:
                if mindelta6 > 0:
                    rows = 6
                else:
                    raise ValueError(err2)
            else:
                raise ValueError(err2)
            data = np.reshape(data, (rows, -1), order='F')
        return data

    def read_op2_tload(self):
        """
        Returns the TLOAD data from an op2 file.

        This routine scans the op2 file for the DYNAMICS datablock and
        then calls :func:`read_op2_dynamics` to read the data.
        """
        if len(self.dbnames) == 0:
            self.directory(verbose=False)
        fpos = self.dbnames['DYNAMICS'][0][0][0]
        self._fileh.seek(fpos)
        name, trailer, unused_dbtype = self._read_op2_name_trailer()
        return self.read_op2_dynamics()

    def read_op2_record(self, form=None, N=0):
        """
        Read Nastran output2 data record.

        Parameters
        ----------
        form : string or None
            String specifying format, or None to read in signed integers.
            One of::
               'int' (same as None)
               'uint'
               'single'
               'double'
               'bytes' -- raw bytes from file
        N : integer
            Number of elements in final data record; use 0 if unknown.

        Returns numpy 1-d vector or, if form=='bytes', a bytes string.

        This routine will read in a 'super' record if the data spans
        more than one logical record.
        """
        key = self._get_key()
        f = self._fileh
        if not form or form == 'int':
            frm = self._intstr
            frmu = self._intstru
            bytes_per = self._ibytes
        elif form == 'uint':
            frm = self._intstr.replace('i', 'u')
            frmu = self._intstru.replace('i', 'I')
            bytes_per = self._ibytes
        elif form == 'double':
            frm = self._endian + b'f8'
            frmu = self._endian + b'%dd'
            bytes_per = 8
        elif form == 'single':
            frm = self._endian + b'f4'
            frmu = self._endian + b'%df'
            bytes_per = 4
        elif form == 'bytes':
            data = b''
            while key > 0:
                reclen = self._Str4.unpack(f.read(4))[0]
                data += f.read(reclen)
                f.read(4)  # endrec
                key = self._get_key()
            self._skip_key(2)
            return data
        else:
            raise ValueError("form must be one of:  None, 'int', "
                             "'uint', 'double', 'single' or 'bytes'")
        if N:
            #print('frm=%r' % frm)
            data = np.zeros(N, dtype=frm)
            i = 0
            while key > 0:
                reclen = self._Str4.unpack(f.read(4))[0]
                # f.read(4)  # reclen
                n = reclen // bytes_per
                if n < self._rowsCutoff:
                    b = n * bytes_per
                    #print('frmu=%r' % frmu)
                    data[i:i+n] = struct.unpack(frmu % n, f.read(b))
                else:
                    data[i:i+n] = np.fromfile(f, frm, n)
                i += n
                f.read(4)  # endrec
                key = self._get_key()
        else:
            data = np.zeros(0, dtype=frm)
            while key > 0:
                reclen = self._Str4.unpack(f.read(4))[0]
                # f.read(4)  # reclen
                n = reclen // bytes_per
                if n < self._rowsCutoff:
                    b = n * bytes_per
                    cur = struct.unpack(frmu % n, f.read(b))
                else:
                    cur = np.fromfile(f, frm, n)
                data = np.hstack((data, cur))
                f.read(4)  # endrec
                key = self._get_key()
        self._skip_key(2)
        return data

    def skip_op2_record(self):
        """
        Skip over Nastran output2 data record (or super-record).
        """
        key = self._get_key()
        while key > 0:
            reclen = self._Str4.unpack(self._fileh.read(4))[0]
            self._fileh.seek(reclen+4, 1)
            key = self._get_key()
        self._skip_key(2)

    def read_op2_table_headers(self, name):
        """
        Read op2 table headers and echo them to the screen.

        Parameters
        ----------
        name : string
            Name of data block that headers are being read for.

        File must be positioned after name and trailer block.  For
        example, to read the table headers of the last GEOM1S data
        block::

            o2 = op2.OP2('modes.op2')
            fpos = o2.dbnames['GEOM1S'][-1][0][0]
            o2._fileh.seek(fpos)
            name, trailer, dbtype = o2._read_op2_name_trailer()
            o2.read_op2_table_headers('GEOM1S')

        """
        key = self._get_key()
        print("{0} Headers:".format(name))
        Frm = struct.Struct(self._intstru % 3)
        eot = 0
        while not eot:
            while key > 0:
                reclen = self._Str4.unpack(self._fileh.read(4))[0]
                head = Frm.unpack(self._fileh.read(3*self._ibytes))
                print(np.hstack((head, reclen)))
                self._fileh.seek((key-3)*self._ibytes, 1)
                self._fileh.read(4)
                key = self._get_key()
            self._skip_key(2)
            eot, key = self._read_op2_end_of_table()

    def _check_code(self, item_code, funcs, vals, name):
        """
        Checks that the code (ACODE or TCODE probably) value is
        acceptable.

        Parameters
        ----------
        item_code : integer
            The ACODE or TCODE (or similar) value for the record.
        funcs : list of integers
            These are the function code values to check for `code`
        vals : list of lists of integers
            These are the acceptable values for the `code` functions;
            ignored if `acode` is None.
        name : string
            Name for message; eg: 'TCODE'

        Returns
        -------
        True if all values are acceptable, False otherwise.

        Notes
        -----
        The function codes in `funcs` are:
            ======  ==========================================
            Code    Operation
            ======  ==========================================
               1    if (item_code//1000 = 2,3,6) then return 2
                    else return 1
               2    mod(item_code,100)
               3    mod(item_code,1000)
               4    item_code//10
               5    mod(item_code,10)
               6    if iand(item_code,8)!=0??? then set to 0,
                    else set to 1
               7    if item_code//1000
                        = 0 or 2, then set to 0
                        = 1 or 3, then set to 1
                        > 3, then set to 2.
            >65535  iand(item_code,iand(func_code,65535))
            ======  ==========================================
        where `iand` is the bit-wise AND operation. For example, ACODE,4
        means that the ACODE value should be integer divided it by 10.
        So, if ACODE is 22, ACODE,4 is 2.
        """
        if len(funcs) != len(vals):
            raise ValueError('len(funcs) != len(vals)!')
        for func, val in zip(funcs, vals):
            if 1 <= func <= 7:
                if self.CodeFuncs[func](item_code) not in val:
                    warnings.warn('{0} value {1} not acceptable; func={2}; allowed={3}'.
                                  format(name, item_code, func, val),
                                  RuntimeWarning)
                    return False
            elif func > 65535:
                if self.CodeFuncs['big'](func, item_code) not in val:
                    warnings.warn('{0} value {1} not acceptable'.
                                  format(name, item_code),
                                  RuntimeWarning)
                    return False
            else:
                raise ValueError('Unknown function code: {0}'.
                                 format(func))
        return True

    def _read_op2_ougv1(self, name):
        """
        Read op2 OUGV1 mode shape data block.

        Parameters
        ----------
        name : string
            Name of OUGV1 data block.

        Returns
        -------
        ougv1 : dict
            Dictionary with::
               'ougv1' : the OUGV1 matrix
               'lambda' : the eigenvalues; len(lambda) = size(ougv1,2)
               'dof' : 2-column matrix of:  [id, dof];
                       size(dof,1) = size(ougv1,1)

        Notes
        -----
        Can currently only read a real eigenvalue table (ACODE,4 = 2,
        TCODE,1 = 1, TCODE,2 = 7, and TCODE,7 in [0, 2]).
        """
        float2_str = struct.Struct(self._endian + b'ff')
        iif6_int = np.dtype(self._endian+'i4')
        iif6_bytes = 32
        if self._ibytes == 4:
            i4_Str = struct.Struct(self._endian + b'iiii')
            i4_bytes = 16
        else:
            i4_Str = struct.Struct(self._endian + b'qqqq')
            i4_bytes = 32
        pos = self._fileh.tell()
        key = self._get_key()
        lam = np.zeros(1, float)
        ougv1 = None
        J = 0
        eot = 0
        while not eot:
            if J == 1:
                # compute number of modes by files bytes:
                startpos = pos + 8 + self._ibytes
                bytes_per_mode = self._fileh.tell() - startpos
                dbdir = self.dbnames[name]
                for i, dbdiri in enumerate(dbdir):
                    if dbdiri[0][0] < startpos < dbdiri[0][1]:
                        endpos = dbdiri[0][1]
                        break
                nmodes = (endpos - startpos) // bytes_per_mode
                print('Number of modes in OUGV1 is {0:d}'.format(nmodes))
                keep = lam
                lam = np.zeros(nmodes, float)
                lam[0] = keep
                keep = ougv1
                ougv1 = np.zeros((keep.shape[0], nmodes), float,
                                 order='F')
                ougv1[:, 0] = keep[:, 0]
            # IDENT record:
            reclen = self._Str4.unpack(self._fileh.read(4))[0]
            header = i4_Str.unpack(self._fileh.read(i4_bytes))
            # header = (ACODE, TCODE, ...)
            achk = self._check_code(header[0], [4], [[2]], 'ACODE')

            # item_code, funcs, vals, name
            tchk = self._check_code(header[1], [1, 2, 7],
                                    [[1], [7], [0, 2]], 'TCODE')
            if not (achk and tchk):
                self._fileh.seek(pos)
                self.skip_op2_table()
                return
            self._fileh.read(self._ibytes)  # mode bytes
            lam[J] = float2_str.unpack(self._fileh.read(8))[0]
            # ttl bytes = reclen + 4 + 3*(4+ibytes+4)
            #           = reclen + 28 - 3*ibytes
            # read bytes = 4*ibytes + ibytes + 8 = 8 + 5*ibytes
            # self._fileh.seek(reclen-2*self._ibytes+20, 1)  # ... or:
            self._fileh.read(reclen-2*self._ibytes+20)

            # DATA record:
            if ougv1 is None:
                #print('masking')
                # - process DOF information on first column only
                # - there are 8 elements per node:
                #   id*10, type, x, y, z, rx, ry, rz
                data = self.read_op2_record('bytes')  # 1st column
                n = len(data) // iif6_bytes
                print('iif6_int =', iif6_int)  # int32
                data = np.fromstring(data, iif6_int)
                data1 = (data.reshape(n, 8))[:, :2]
                pvgrids = data1[:, 1] == 1
                dof = expand_dof(data1[:, 0] // 10, pvgrids)
                # form partition vector for modeshape data:
                V = np.zeros((n, 8), bool)
                V[:, 2] = True          # all nodes have 'x'
                V[pvgrids, 3:] = True   # only grids have all 6
                 #print('V =\n', V)
                V = V.flatten()
                # initialize ougv1 with first mode shape:
                data.dtype = np.float32  # reinterpret as floats
                ougv1 = data[V].reshape(-1, 1)
            else:
                data = self.read_op2_record('single', V.shape[0])
                ougv1[:, J] = data[V]
            J += 1
             #print('Finished reading mode {0:3d}, Frequency ={1:6.2f}'.format(
            #    J, np.sqrt(lam[J-1])/(2*np.pi)))
            eot, key = self._read_op2_end_of_table()
        return {'ougv1': ougv1, 'lambda': lam, 'dof': dof}

    def _read_op2_emap(self, nas, nse, trailer):
        """
        Read Nastran output2 EMAP data block.

        Parameters
        ----------
        nas : dict
            Dictionary; has at least {'dnids': {}}.
        nse : integer
            Number of superelements.
        trailer : 1-d array
            The trailer for the EMAP data block.

        Fills in the dnids member of nas.

        See :func:`read_nas2cam_op2`.
        """
        words4bits = trailer[4]
        data1 = self.read_op2_record()
        # [se bitpos proc_order dnse bitpos_dnse prim_se se_type]
        data1 = np.reshape(data1[:7*nse], (-1, 7))

        # read 2nd record:
        key = self._get_key()
        data2 = np.zeros(0, dtype='u4')
        frm = self._uendian + 'u4'
        frmu = self._endian + b'%dI'
        if self._ibytes == 8:
            mult = 2
        else:
            mult = 1
        while key > 0:
            self._fileh.read(4)  # reclen
            if mult*key < self._rowsCutoff:
                cur = struct.unpack(frmu % (mult*key),
                                    self._fileh.read(4*mult*key))
            else:
                cur = np.fromfile(self._fileh, frm, mult*key)
            data2 = np.hstack((data2, cur))
            self._fileh.read(4)  # endrec
            key = self._get_key()
        if self._ibytes == 8:
            data2 = np.reshape(data2, (4, -1))
            data2 = data2[[0, 3], :].flatten()
        self._skip_key(2)

        # [ grid_id [bitmap] ]
        data2 = np.reshape(data2, (-1, words4bits))
        # 1 in front need to skip over grid_id (vars are col indices)
        word4bit_up = 1 + data1[:, 1] // 32
        word4bit_dn = 1 + data1[:, 4] // 32
        bitpos_up = 31 - data1[:, 1] % 32
        bitpos_dn = 31 - data1[:, 4] % 32
        for j in range(nse-1):
            se = data1[j, 0]
            bitdn = 1 << bitpos_dn[j]
            bitup = 1 << bitpos_up[j]
            connected = np.logical_and(data2[:, word4bit_dn[j]] & bitdn,
                                       data2[:, word4bit_up[j]] & bitup)
            grids = data2[connected, 0]
            nas['dnids'][se] = grids
        for j in range(nse):  # = 1 to nse:
            self.skip_op2_record()
        self._get_key()

    def _read_op2_bgpdt(self):
        """
        Read record 1 of the Nastran output2 BGPDT data block.

        Returns vector of the BGPDT data or [] if no data found.
        Vector is 9*ngrids in length.  For each grid:
        ::
          [ coord_id
            internal_id
            external_id
            dof_type;
            permanent_set_constraint
            boundary_grid_id
            x
            y
            z ]

        The x, y, z values are the grid location in basic.

        See :func:`rdn2cop2`.
        """
        if self._ibytes == 4:
            Str = struct.Struct(self._endian + b'iiiiiiddd')
            sbytes = 24 + 24
            wpg = 12  # words per grid
            wpd = 2   # words per double
        else:
            Str = struct.Struct(self._endian + b'qqqqqqddd')
            sbytes = 48 + 24
            wpg = 9   # words per grid
            wpd = 1   # words per double
        rfrm = self._endian + b'%dd'
        key = self._get_key()

        datarec = []
        ileft = 0     # remaining left over
        dleft = 0     # remaining doubles left over

        a = np.arange(6)
        b = np.arange(6, 9)
        v = np.arange(9)
        A = grids = 0
        while key > 0:
            self._fileh.read(4)  # reclen
            if ileft > 0:
                i = A + a[6-ileft:] + grids*9
                # datarec[i] = np.fromfile(self._fileh, self._intstr, ileft)
                bytes = self._ibytes * ileft
                datarec[i] = struct.unpack(self._intstru % ileft,
                                           self._fileh.read(bytes))
            if dleft > 0:
                i = A + b[3-dleft:] + grids*9
                # datarec[i] = np.fromfile(self._fileh, rfrm, dleft)
                datarec[i] = struct.unpack(rfrm % dleft,
                                           self._fileh.read(8*dleft))

            key = key - ileft - dleft*wpd
            # number of complete grids remaining in this record:
            grids = key // wpg
            A = len(datarec)
            # round up for memory allocation (for possible partial):
            n = (key + wpg - 1) // wpg
            datarec = np.hstack((datarec, np.zeros(n*9)))

            Av = A + v
            for i in range(grids):
                datarec[Av + i*9] = Str.unpack(self._fileh.read(sbytes))

            # read in remainder of record if any
            ileft = 0
            dleft = 0
            if key > grids*wpg:
                # number of words left (1 word/int, 2 words/double)
                n = key - grids*wpg
                if n >= 6:
                    i = A+a+grids*9
                    # datarec[i] = np.fromfile(self._fileh, self._intstr, 6)
                    bytes = self._ibytes * 6
                    datarec[i] = struct.unpack(self._intstru % 6,
                                               self._fileh.read(bytes))
                    # divide by wpd to get number of doubles
                    n = (n - 6) // wpd
                    dleft = 3-n
                    if n >= 1:
                        i = A + b[:n] + grids*9
                        # datarec[i] = np.fromfile(self._fileh, rfrm, n)
                        datarec[i] = struct.unpack(rfrm % n,
                                                   self._fileh.read(8*n))
                else:
                    i = A + a[:n] + grids*9
                    # datarec[i] = np.fromfile(self._fileh, self._intstr, n)
                    bytes = self._ibytes * n
                    datarec[i] = struct.unpack(self._intstru % n,
                                               self._fileh.read(bytes))
                    ileft = 6-n
                    dleft = 3
            self._fileh.read(4)  # endrec
            key = self._get_key()
        self._skip_key(2)
        return datarec

    def _read_op2_bgpdt68(self):
        """
        Read record 1 of the Nastran output2 BGPDT68 data block.

        Returns vector of the BGPDT data or [] if no data found.
        Vector is 4*ngrids in length.  For each grid:
        ::
          [ coord_id
            x
            y
            z ]

        The x, y, z values are the grid location in basic.
        """
        struc = struct.Struct(self._endian + b'ifff')
        sbytes = 16
        wpg = 4  # words per grid
        wpd = 1  # words per single
        rfrm = self._endian + b'%df'
        key = self._get_key()

        datarec = []
        ileft = 0     # remaining left over
        dleft = 0     # remaining doubles left over

        a = np.arange(1)
        b = np.arange(1, 4)
        v = np.arange(4)
        A = grids = 0
        while key > 0:
            self._fileh.read(4)  # reclen
            if ileft > 0:
                i = A + grids*4
                # datarec[i] = np.fromfile(self._fileh, self._int32str, ileft)
                bytes = 4 * ileft
                datarec[i] = struct.unpack(self._int32stru % ileft,
                                           self._fileh.read(bytes))
            if dleft > 0:
                i = A + b[3-dleft:] + grids*4
                # datarec[i] = np.fromfile(self._fileh, rfrm, dleft)
                datarec[i] = struct.unpack(rfrm % dleft,
                                           self._fileh.read(4*dleft))
            key = key - ileft - dleft*wpd
            # number of complete grids remaining in this record:
            grids = key // wpg
            A = len(datarec)
            # round up for memory allocation (for possible partial):
            n = (key + wpg - 1) // wpg
            datarec = np.hstack((datarec, np.zeros(n*4)))

            Av = A + v
            for i in range(grids):
                datarec[Av + i*4] = struc.unpack(self._fileh.read(sbytes))

            # read in remainder of record if any
            ileft = 0
            dleft = 0
            if key > grids*wpg:
                # number of words left (1 word/int, 2 words/double)
                n = key - grids*wpg
                if n >= 1:
                    i = A + a + grids*4
                    datarec[i] = self._Str4.unpack(self._fileh.read(4))[0]
                    # divide by wpd to get number of doubles
                    n = (n - 1) // wpd
                    dleft = 3-n
                    if n >= 1:
                        i = A + b[:n] + grids*4
                        # datarec[i] = np.fromfile(self._fileh, rfrm, n)
                        datarec[i] = struct.unpack(rfrm % n,
                                                   self._fileh.read(4*n))
            self._fileh.read(4)  # endrec
            key = self._get_key()
        self._skip_key(2)
        return datarec

    def _read_op2_cstm(self):
        """
        Read Nastran output2 CSTM data block.

        Returns 14-column matrix 2-d array of the CSTM data:
        ::
          [
           [ id1 type xo yo zo T(1,1:3) T(2,1:3) T(3,1:3) ]
           [ id2 type xo yo zo T(1,1:3) T(2,1:3) T(3,1:3) ]
           ...
          ]

        T is transformation from local to basic for the coordinate
        system.

        See :func:`read_nas2cam_op2`.
        """
        cstm_rec1 = self.read_op2_record()
        cstm_rec2 = self.read_op2_record('double')
        self._read_op2_end_of_table()

        # assemble coordinate system table
        length = len(cstm_rec1)
        cstm = np.zeros((length/4, 14))
        cstm[:, 0] = cstm_rec1[::4]
        cstm[:, 1] = cstm_rec1[1::4]
        # start index into rec2 for xo, yo, zo, T (12 values) is in
        # last (4th) position in rec1 for each coordinate system:
        pv = range(12)
        for i, j in enumerate(cstm_rec1[3::4]):
            cstm[i, 2:] = cstm_rec2[j+pv-1]  # -1 for 0 offset
        return cstm

    def _read_op2_cstm68(self):
        """
        Read record 1 of Nastran output2 CSTM68 data block.

        Returns vector of the CSTM data or [] if no data found.  Vector
        is 14 * number of coordinate systems in length.  For each
        coordinate system:
        ::
          [ id type xo yo zo T(1,1:3) T(2,1:3) T(3,1:3) ]

        T is transformation from local to basic for the coordinate
        system.
        """
        Str = struct.Struct(self._endian + b'ii' + b'f'*12)
        sbytes = 4 *14
        wpg = 14   # words per grid
        wpd = 1    # words per single
        key = self._get_key()
        rfrm = self._endian + b'%df'
        datarec = []
        ileft = 0    # integers to read that are left over
        dleft = 0    # singles left

        a = np.arange(2)
        b = np.arange(2, 14)
        v = np.arange(14)
        A = grids = 0
        while key > 0:
            self._fileh.read(4)  # reclen
            if ileft > 0:
                i = A + a[2-ileft:] + grids*14
                # datarec[i] = np.fromfile(self._fileh,
                #                          self._int32str, ileft)
                bytes = 4 * ileft
                datarec[i] = struct.unpack(self._int32stru % ileft,
                                           self._fileh.read(bytes))
            if dleft > 0:
                i = A + b[12-dleft:] + grids*14
                # datarec[i] = np.fromfile(self._fileh, rfrm, dleft)
                datarec[i] = struct.unpack(rfrm % dleft,
                                           self._fileh.read(4*dleft))

            key = key - ileft - dleft*wpd
            # number of complete grids remaining in this record
            grids = key // wpg
            A = len(datarec)
            # round up for memory allocation (for possible partial):
            n = (key + wpg - 1) // wpg
            datarec = np.hstack((datarec, np.zeros(n*14)))

            Av = A + v
            for i in range(grids):
                datarec[Av + i*14] = Str.unpack(self._fileh.read(sbytes))

            # read in remainder of record if any
            ileft = 0
            dleft = 0
            if key > grids*wpg:
                # number of words left (1 word/int, 2 words/single)
                n = key - grids*wpg
                if n >= 2:
                    i = A + a + grids*14
                    # datarec[i] = np.fromfile(self._fileh,
                    #                          self._int32str, 2)
                    datarec[i] = struct.unpack(self._int32stru % 2,
                                               self._fileh.read(8))
                    # divide by wpd to get number of singles
                    n = (n - 2) // wpd
                    dleft = 12-n
                    if n >= 1:
                        i = A + b[:n] + grids*14
                        # datarec[i] = np.fromfile(self._fileh, rfrm, n)
                        datarec[i] = struct.unpack(rfrm % n,
                                                   self._fileh.read(4*n))
                else:
                    # n must be 1 here
                    i = A + grids*14
                    datarec[i] = self._Str4.unpack(self._fileh.read(4))[0]
                    ileft = 2-n
                    dleft = 12
            self._fileh.read(4)  # endrec
            key = self._get_key()
        self._skip_key(2)
        self._read_op2_end_of_table()
        return datarec

    def _read_op2_geom1_cord2(self):
        if self._ibytes == 4:
            header_Str = struct.Struct(self._endian + b'iii')
            cord2_Str = struct.Struct(self._endian + b'4i9f')
            sebulk_Str = struct.Struct(self._endian + b'4if3i')
            hbytes = 12
            cbytes = 4*13
            bbytes = 4*8
        else:
            header_Str = struct.Struct(self._endian + b'qqq')
            cord2_Str = struct.Struct(self._endian + b'4q9d')
            sebulk_Str = struct.Struct(self._endian + b'4qd3q')
            hbytes = 24
            cbytes = 8*13
            bbytes = 8*8

        CORD2R = (2101, 21, 8)
        CORD2C = (2001, 20, 9)
        CORD2S = (2201, 22, 10)
        SEBULK = (1427, 14, 465)
        SECONCT = (427, 4, 453)

        cord2 = np.zeros((0, 13))
        sebulk = np.zeros((1, 8))
        selist = np.array([[0, 0]], int)
        key = self._get_key()
        eot = 0
        # data = np.zeros(0, dtype=self._intstr)
        while not eot:
            while key > 0:
                self._fileh.read(4)  # reclen
                # reclen = self._Str4.unpack(self._fileh.read(4))[0]
                head = header_Str.unpack(self._fileh.read(hbytes))
                if head in [CORD2R, CORD2C, CORD2S]:
                    n = (key-3) // 13
                    data = np.zeros((n, 13))
                    for i in range(n):
                        data[i] = cord2_Str.unpack(self._fileh.read(cbytes))
                    cord2 = np.vstack((cord2, data))
                elif head == SEBULK:
                    n = (key-3) // 8
                    sebulk = np.zeros((n, 8))
                    for i in range(n):
                        sebulk[i] = sebulk_Str.unpack(self._fileh.read(bbytes))
                elif head == SECONCT:
                    n = key - 3
                    if n < self._rowsCutoff:
                        nbytes = n * self._ibytes
                        seconct = np.zeros(n, int)
                        seconct[:] = struct.unpack(self._intstru % n,
                                                   self._fileh.read(nbytes))
                    else:
                        seconct = np.fromfile(self._fileh, self._intstr, n)
                    pv = np.nonzero(seconct == -1)[0][1:-2:2] + 1
                    pv = np.hstack((0, pv))
                    u = np.unique(seconct[pv], return_index=True)[1]
                    pv = pv[u]
                    selist = np.vstack((seconct[pv], seconct[pv+1])).T
                    selist = np.vstack((selist, [0, 0]))
                else:
                    self._fileh.seek((key-3)*self._ibytes, 1)
                self._fileh.read(4)  # endrec
                key = self._get_key()
            self._skip_key(2)
            eot, key = self._read_op2_end_of_table()
        cord2 = np.delete(cord2, 2, axis=1)
        return n2y.build_coords(cord2), sebulk, selist

    def _read_op2_selist(self):
        """
        Read SLIST data block and return `selist` for :func:`read_nas2cam_op2`.

        See :func:`read_nas2cam_op2`.
        """
        slist = self.read_op2_record()
        slist[1::7] = 0
        self.skip_op2_record()
        self._read_op2_end_of_table()
        return np.vstack((slist[::7], slist[4::7])).T

    def _read_op2_uset(self):
        """
        Read the USET data block.

        Returns 1-d USET array.  The 2nd bit is cleared for the S-set.

        See :func:`rdn2cop2`.
        """
        uset = self.read_op2_record('uint')
        # clear the 2nd bit for all S-set:
        sset = 0 != (uset & n2y.mkusetmask("s"))
        if any(sset):
            uset[sset] = uset[sset] & ~2
        self._read_op2_end_of_table()
        return uset

    def _read_op2_eqexin(self):
        """
        Read the EQEXIN data block.

        Returns (EQEXIN1, EQEXIN) tuple.

        See :func:`read_nas2cam_op2`.
        """
        eqexin1 = self.read_op2_record()
        eqexin = self.read_op2_record()
        self._read_op2_end_of_table()
        return eqexin1, eqexin

    def _proc_bgpdt(self, eqexin1, eqexin, ver68=False, bgpdtin=None):
        """
        Reads and processes the BGPDT data block for :func:`read_nas2cam_op2`
        and :func:`read_post_op2`.

        Returns (bgpdt, dof, doftype, nid, upids)

        See :func:`read_nas2cam_op2`, :func:`read_post_op2`.
        """
        if ver68:
            bgpdt_rec1 = bgpdtin
        else:
            bgpdt_rec1 = self._read_op2_bgpdt()
            self.read_op2_record()
            self.skip_op2_table()

        # assemble coordinates table
        # bgpdt: [x, y, z, cid]
        if ver68:
            bgpdt_rec1 = bgpdt_rec1.reshape((-1, 4))
            bgpdt = bgpdt_rec1[:, [1, 2, 3, 0]]
        else:
            bgpdt_rec1 = bgpdt_rec1.reshape((-1, 9))
            bgpdt = bgpdt_rec1[:, [6, 7, 8, 0]]

        # assemble dof table:
        dof = eqexin[1::2] // 10
        doftype = eqexin[1::2] - 10*dof
        nid = eqexin[::2]

        # eqexin is in external sort, so sort it
        i = eqexin1[1::2].argsort()
        dof = dof[i]
        doftype = doftype[i]
        nid = nid[i]
        if ver68:
            upids = None
        else:
            upids = bgpdt_rec1[:, 5].astype(int)
        return bgpdt, dof, doftype, nid, upids

    def _build_Uset(self, se, dof, doftype, nid, uset, bgpdt,
                    cstm=None, cstm2=None):
        """
        Builds the 6-column uset table for :func:`rdn2cop2` and
        :func:`rdpostop2`.

        Returns: (uset, cstm, cstm2).

        See :func:`read_nas2cam_op2`.
        """
        # Fill in all dof use -1 as default and set dof as
        # appropriate ... make it big enough for grids (6 cols).
        # The -1s will be partitioned out below.
        rd = len(dof)
        rb = np.size(bgpdt, 0)
        if rd != rb:
            raise RuntimeError(
                'RDOP2USET:  BGPDTS incompatible with '
                'EQEXINS for superelement {}.\n'
                '  Guess:  residual run clobbered EQEXINS\n'
                '    Fix:  add the "fxphase0" alter to your '
                'residual run'.format(se))
        coordinfo = np.zeros((rd, 18))
        coordinfo[:, :4] = bgpdt
        if cstm is None:
            n = len(cstm2)
            cstm = np.zeros((n, 14))
            for i, key in enumerate(cstm2):
                cstm[i, :2] = cstm2[key][0, :2]
                cstm[i, 2:] = (cstm2[key].flatten())[3:]
        cref = cstm[:, 0].astype(int)
        c_all = bgpdt[:, 3].astype(int)
        i = np.argsort(cref)
        pv = i[np.searchsorted(cref, c_all, sorter=i)]
        coordinfo[:, 4] = cstm[pv, 1]
        coordinfo[:, 6:] = cstm[pv, 2:]

        grids = doftype == 1
        ngrids = np.sum(grids)
        nongrids = rd - ngrids
        doflist = np.zeros((rd, 6)) - 1
        if ngrids > 0:
            doflist[grids, :] = np.arange(1, 7)
        if nongrids > 0:
            doflist[grids == False, 0] = 0
        doflist = doflist.flatten()
        idlist = np.dot(nid.reshape(-1, 1), np.ones((1, 6))).flatten()
        coordinfo = coordinfo.reshape((rd*6, 3))

        # partition out -1s:
        pv = doflist != -1
        doflist = doflist[pv]
        idlist = idlist[pv]
        coordinfo = coordinfo[pv, :]
        if uset is None:
            warnings.warn('uset information not found.  Putting all '
                          'DOF in b-set.', RuntimeWarning)
            #import n2y
            b = n2y.mkusetmask('b')
            uset = np.zeros(len(doflist), int) + b
        uset = np.hstack((np.vstack((idlist, doflist, uset)).T,
                          coordinfo))
        if cstm2 is None:
            cstm2 = {}
            for row in cstm:
                m = np.zeros((5, 3))
                m[0, :2] = row[:2]
                m[1:, :] = row[2:].reshape((4, 3))
                cstm2[int(row[0])] = m
        return uset, cstm, cstm2

    def _read_op2_maps(self):
        """
        Reads and returns the MAPS information for :func:`read_nas2cam_op2`.
        """
        if self._ibytes == 4:
            id_Str = struct.Struct(self._endian + b'id')
            id_bytes = 12
        else:
            id_Str = struct.Struct(self._endian + b'qd')
            id_bytes = 16
        key = 1
        maps = np.zeros((0, 2))
        while key:
            key = self._get_key()  # 2 (1 integer, 1 double)
            self._fileh.read(4)  # reclen 12 or 16 bytes
            curmap = id_Str.unpack(self._fileh.read(id_bytes))
            maps = np.vstack((maps, curmap))
            self._fileh.read(4)  # endrec
            self._skip_key(2)  # 1st key is mystery negative
            key = self._get_key()  # 1 if cont, 0 if done
        self._get_key()
        maps[:, 0] -= 1
        return maps

    def _read_op2_drm(self):
        """
        Read Nastran output2 DRM data block (table).

        Returns tuple:  (drm, iddof)
        ----------------------------
        drm : ndarray
            The drm matrix.
        iddof : ndarray
            2-column matrix of [id, dof].

        This routine is beta -- check output carefully.
        """
        def get_str(iprev, elemtype, ir_str, ir_bytes):
            """
            ir_str : ???
            """
            if np.any(elemtype == np.array([4, 5])):
                ints_rec2 = 1
            else:
                ints_rec2 = 2
            if ints_rec2 != iprev:
                if self._bit64:
                    ir_str = struct.Struct(self._endian + b'q'*ints_rec2)
                    ir_bytes = 8*ints_rec2
                else:
                    ir_str = struct.Struct(self._endian + b'i'*ints_rec2)
                    ir_bytes = 4*ints_rec2
            return ir_str, ir_bytes, ints_rec2

        if self._bit64:
            rfrm = self._endian + b'f8'
            rfrmu = self._endian + b'%dd'
            rsize = 8
        else:
            rfrm = self._endian + b'f4'
            rfrmu = self._endian + b'%df'
            rsize = 4
        u1 = self.read_op2_record()
        elemtype = u1[1]
        elemid = u1[2]
        ir_str, ir_bytes, ints_rec2 = get_str(0, elemtype, None, None)
        nwords = u1[9]
        key = self._get_key()
        block = 7*4+3*self._ibytes

        # determine records/column by scanning first column:
        rpc = 0
        fp = self._fileh
        pos = fp.tell()
        id1 = -1
        drmrow = 0
        blocksize = 500   # number of rows or cols to grow iddof and drm
        drmrows = blocksize
        iddof = np.zeros((drmrows, 2), int)
        KN = key, nwords
        while key >= nwords:
            L = nwords - ints_rec2
            fp.read(4)   # reclen
            dataint = ir_str.unpack(fp.read(ir_bytes))
            id_cur = dataint[0] // 10
            if id1 == -1:
                id1 = id_cur
            elif id1 == id_cur:
                break
            rpc += 1
            if drmrow+L >= drmrows:
                iddof = np.vstack((iddof,
                                   np.zeros((blocksize, 2), int)))
                drmrows += blocksize
            iddof[drmrow:drmrow+L, 0] = id_cur
            iddof[drmrow:drmrow+L, 1] = elemid
            fp.seek(self._ibytes*L, 1)
            drmrow += L

            # read rest of record:
            for unused_i in range(1, key // nwords):
                dataint = ir_str.unpack(fp.read(ir_bytes))
                id_cur = dataint[0] // 10
                if drmrow+L >= drmrows:
                    iddof = np.vstack((iddof,
                                       np.zeros((blocksize, 2), int)))
                    drmrows += blocksize
                iddof[drmrow:drmrow+L, 0] = id_cur
                iddof[drmrow:drmrow+L, 1] = elemid
                fp.seek(self._ibytes*L, 1)
                drmrow += L
            fp.seek(block, 1)
            key = self._get_key()
            if key > 0:
                fp.read(4)   # reclen
                if key < self._rowsCutoff:
                    u1 = struct.unpack(self._intstru % key,
                                       fp.read(key*self._ibytes))
                else:
                    u1 = np.fromfile(fp, self._intstr, key)
                if u1[1] != elemtype:
                    raise RuntimeError('u1[1] != elemtype')
                # above check precludes next two lines:
                # elemtype = u1[1]
                # ir_str, ir_bytes, ints_rec2 = get_str(ints_rec2,
                #                                       elemtype,
                #                                       ir_str, ir_bytes)
                # if u1[2] != elemid:
                #     raise RuntimeError('u1[2] != elemid ... should it?')
                elemid = u1[2]
                if u1[9] != nwords:
                    raise RuntimeError('u1[9] != nwords ... should it?')
                # nwords = u1[9]
                fp.seek(block, 1)
                key = self._get_key()

        drmrows = drmrow
        iddof = iddof[:drmrows]
        drmcols = blocksize
        fp.seek(pos)
        B = np.zeros((drmrows, drmcols), order='F')
        drm = B.copy()
        drmcol = 0
        key, nwords = KN
        while key >= nwords:
            drmrow = 0
            if drmcol == drmcols:
                drm = np.asfortranarray(np.hstack((drm, B)))
                drmcols += blocksize
            for unused_ in it.repeat(None, rpc):
                L = nwords - ints_rec2
                fp.read(4)   # reclen
                for unused_i in range(key // nwords):
                    # dataint = ir_Str.unpack(fp.read(ir_bytes))
                    fp.read(ir_bytes)
                    if L < self._rowsCutoff:
                        drm[drmrow:drmrow+L, drmcol] = struct.unpack(
                            rfrmu % L, fp.read(rsize*L))
                    else:
                        drm[drmrow:drmrow+L, drmcol] = np.fromfile(fp, rfrm, L)
                    drmrow += L
                fp.seek(block, 1)
                key = self._get_key()
                if key > 0:
                    fp.read(4)   # reclen
                    if key < self._rowsCutoff:
                        u1 = struct.unpack(self._intstru % key,
                                           fp.read(key*self._ibytes))
                    else:
                        u1 = np.fromfile(fp, self._intstr, key)
                else:
                    break
                fp.seek(block, 1)
                key = self._get_key()
            drmcol += 1
        return drm[:, :drmcol], iddof

    def read_drm2op2(self, verbose=False):
        """
        Read op2 file output by DRM2 DMAP.

        Parameters
        ----------
        verbose : bool
            If true, echo names of tables and matrices to screen.

        Returns
        -------
        drmkeys : dictionary
            - 'dr' : data recovery items in order requested (from
            XYCDBDRS)
            - 'drs' : sorted version of 'dr' (from XYCDBDRS)
            - 'tougv1', 'tougs1', etc : directories corresponding to
            the data recovery matrices (which are written to op4).
            All of these start with 'to' (lower case).

        File is created with a header and then these data blocks are
        written:
        ::
          OUTPUT2  XYCDBDRS//0/OP2UNIT $
          OUTPUT2  TOUGV1,TOUGS1,TOUGD1//0/OP2UNIT $
          OUTPUT2  TOQGS1,TOQGD1,TOEFS1,TOEFD1//0/OP2UNIT $
          OUTPUT2  TOESS1,TOESD1//0/OP2UNIT $
        """
        self._fileh.seek(self._postheaderpos)
        drmkeys = {}
        self.verbose = verbose
        while 1:
            name, trailer, rectype = self._read_op2_name_trailer()
            if name is None:
                break
            if rectype > 0:
                if verbose:
                    print("Skipping matrix %r..." % name)
                self.skip_op2_matrix(trailer)
                # matrix = self.rdop2matrix(trailer)
            elif len(name) > 2 and name.find('TO') == 0:
                if verbose:
                    print("Reading %r..." % name)
                # self.skipop2table()
                # skip record 1
                self.read_op2_record()
                # record 2 contains directory
                # - starting at 10: type, id, number, row, 0
                info = self.read_op2_record()[10:]
                drmkeys[name.lower()] = (info.reshape(-1, 5).T)[:-1]
                self._read_op2_end_of_table()
            elif len(name) > 4 and name[:4] == 'XYCD':
                if verbose:
                    print("Reading %r..." % name)
                # record 1 contains order of request info
                drmkeys['dr'] = self.read_op2_record()
                # record 2 contains sorted list
                drmkeys['drs'] = self.read_op2_record().reshape(-1, 6).T
                self._read_op2_end_of_table()
            else:
                if verbose:
                    print("Skipping table %r..." % name)
                self.skip_op2_table()
        return drmkeys

    def read_nas2cam_op2(self):
        """
        Read Nastran output2 file written by DMAP NAS2CAM; usually
        called by :func:`rdnas2cam`.

        Returns dictionary with the following members:

        'selist' : array
            2-columns matrix:  [ seid, dnseid ] where, for each row,
            dnseid is the downstream superelement for seid. (dnseid = 0
            if seid = 0).
        'uset' : dictionary
            Indexed by the SE number.  Each member is a 6-column matrix
            described below.
        'cstm' : dictionary
            Indexed by the SE number.  Each member is a 14-column matrix
            containing the coordinate system transformation matrix for
            each coordinate system.  See description below.
        'cstm2' : dictionary
            Indexed by the SE number.  Each member is another dictionary
            indexed by the coordinate system id number.  This has the
            same information as 'cstm', but in a different format.  See
            description below.
        'maps' : dictionary
            Indexed by the SE number.  Each member is a mapping table
            for mapping the A-set order from upstream to downstream;
            see below.
        'dnids' : dictionary
            Indexed by the SE number.  Each member is a vector of ids of
            the A-set ids of grids and spoints for SE in the downstream
            superelement.  When using the CSUPER entry, these will be
            the ids on that entry.  (Does not have each DOF, just ids.)
        'upids' : dictionary
            Indexed by the SE number.  Each member is a vector of ids of
            the A-set grids and spoints for upstream se's that connect
            to SE.  These ids are internally generated and should match
            with 'dnids'.  This allows, for example, the routine
            :func:`n2y.upasetpv` to work.  (Does not have each DOF, just
            ids.)

        The module n2y has many routines that use the data created by
        this routine.

        'uset' description
        ------------------
        Each USET variable is a 6-column matrix where the rows
        correspond to the DOF in Nastran internal sort, and the columns
        are:
        ::
            USET = [ ID DOF Set_Membership Coord_Info ]

        where, for grids, Coord_Info is a 6 row by 3 column matrix:
        ::
            Coord_Info = [[x   y    z]    # location of node in basic
                          [id  type 0]    # coord. id and type
                          [xo  yo  zo]    # origin of coord. system
                          [    T     ]]   # 3x3 transformation to basic
                                          #  for coordinate system

            Coord_Info = [ 0 0 0 ]        # for SPOINTs

        'cstm' description
        ------------------
        Each CSTM contains all the coordinate system information for
        the superelement.  Some or all of this info is in the USET
        table, but if a coordinate system is not used as an output
        system of any grid, it will not show up in USET.  That is why
        CSTM is here.  CSTM has 14 columns:
        ::
            CSTM = [ id type xo yo zo T(1,:) T(2,:) T(3,:) ]

        Note that each CSTM always starts with the two ids 0 and -1.
        The 0 is the basic coordinate system and the -1 is a dummy for
        SPOINTs.  Note the T is transformation between coordinate
        systems as defined (not necessarily the same as the
        transformation for a particular grid ... which, for
        cylindrical and spherical, depends on grid location).  This is
        the same T as in the USET table.

        For example, to convert coordinates from global to basic:
        ::
          Rectangular (type = 1):
             [x; y; z] = T*[xg; yg; zg] + [xo; yo; zo]

          Cylindrical (type = 2):
             % c = cos(theta); s = sin(theta)
             [x; y; z] = T*[R c; R s; zg] + [xo; yo; zo]

          Spherical (type = 3):
             % s1 = sin(theta); s2 = sin(phi)
             [x; y; z] = T*[r s1 c2; r s1 s2; r c1] + [xo; yo; zo]


        'cstm2' description
        ------------------
        Each CSTM2 is a dictionary with the same 5x3 that the
        'Coord_Info' listed above has (doesn't include the first row
        which is the node location).  The dictionary is indexed by the
        coordinate id.

        'maps' description
        ------------------
        MAPS will be [] for superelements whose A-set dof did not get
        rearranged going downstream (on the CSUPER entry.)  For other
        superelements, MAPS will contain two columns: [order, scale].
        The first column reorders upstream A-set to be in the order
        that they appear in the downstream: Down = Up(MAPS(:,1)).  The
        second column is typically 1.0; if not, these routines will
        print an error message and stop.  Together with DNIDS, a
        partition vector can be formed for the A-set of an upstream
        superelement (see :func:`n2y.upasetpv`).

        The op2 file that this routine reads is written by the Nastran
        DMAP NAS2CAM.  The data in the file are expected to be in this
        order:
        ::
            SLIST & EMAP or  SUPERID
            For each superelement:
              USET
              EQEXINS
              CSTMS    (if required)
              BGPDTS
              MAPS     (if required)

        .. note:: The 2nd bit for the DOF column of all USET tables is
        cleared for all S-set.  See :func:`n2y.mkusetmask` for more
        information.

        See rdnas2cam, n2y.

        Example:
        ::
          import op2
          import n2y
          # list superelement 100 DOF that are in the B set:
          o2 = op2.OP2('nas2cam.op2')
          nas = op2.rdn2cop2()
          bset = n2y.mksetpv(nas['uset'][100], 'p', 'b')
          print('bset of se100 = ', nas['uset'][100][bset, :2])
        """
        # setup basic coordinate system info and a dummy for spoints:
        bc = np.array([[+0, 1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1],
                       [-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]])
        nas = {'uset': {},
               'cstm': {},
               'cstm2': {},
               'maps': {},
               'dnids': {},
               'upids': {}}
        self._fileh.seek(self._postheaderpos)
        # read datablock (slist) header record:
        name, trailer, dbtype = self._read_op2_name_trailer()
        if dbtype > 0:
            selist = np.hstack((self.read_op2_matrix(trailer), [[0]]))
            selist = selist.astype(int)
            name, trailer, dbtype = self._read_op2_name_trailer()
        else:
            selist = self._read_op2_selist()
            nse = np.size(selist, 0)
            name, trailer, dbtype = self._read_op2_name_trailer()
            if name == "EMAP":
                self._read_op2_emap(nas, nse, trailer)
                name, trailer, dbtype = self._read_op2_name_trailer()

        # read uset and eqexins tables and do some processing:
        for se in selist[:, 0]:
            if not name:
                break
            uset = self._read_op2_uset()
            name, trailer, dbtype = self._read_op2_name_trailer()
            eqexin1, eqexin = self._read_op2_eqexin()
            name, trailer, dbtype = self._read_op2_name_trailer()
            if name == "CSTMS":
                cstm = np.vstack((bc, self._read_op2_cstm()))
                name, trailer, dbtype = self._read_op2_name_trailer()
            else:
                cstm = bc
            bgpdt, dof, doftype, nid, upids = self._proc_bgpdt(eqexin1, eqexin)
            nas['upids'][se] = upids
            Uset, cstm, cstm2 = self._build_Uset(se, dof, doftype, nid,
                                                 uset, bgpdt, cstm, None)
            nas['uset'][se] = Uset
            nas['cstm'][se] = cstm
            nas['cstm2'][se] = cstm2
            name, trailer, dbtype = self._read_op2_name_trailer()
            if name == "MAPS":
                nas['maps'][se] = self._read_op2_maps()
                name, trailer, dbtype = self._read_op2_name_trailer()
            else:
                nas['maps'][se] = []
        nas['selist'] = selist
        return nas


def read_nas2cam(op2file='nas2cam', op4file=None):
    """
    Read op2/op4 data written by the DMAP NAS2CAM.

    Parameters
    ----------
    op2file : string
        Either the basename of the .op2 and .op4 files, or the full
        name of the .op2 file
    op4file : string or None
        The name of the .op4 file or, if None, builds name from the
        `op2file` input.

    Returns dictionary with these members:
    -------------------------------------
        All members created by :func:`rdn2cop2` (see that routine's
        help).

        'nrb' : integer
            The number of rigid-body modes for residual.

        All the following members are dictionaries indexed by SE number
        (note that this routine will read all matrices for each SE, not
        just those listed here):

        'ulvs' : dictionary
            The ULVS matrices (row partitions of residual modes to the
            A-set DOF of the SE).
        'lambda' : dictionary
            The eigenvalues for each SE.
        'gm' : dictionary
            N-set to M-set transformation matrix GM:  M = GM N.
        'got' : dictionary
            constraint modes
        'goq' : dictionary
            normal modes
        'rfmodes' : dictionary
            index partition vector for res-flex modes
        'maa' : dictionary
            A-set mass
        'baa' : dictionary
            A-set damping
        'kaa' : dictionary
            A-set stiffness
        'pha' : dictionary
            A-set modes
        'mdd' : dictionary
            D-set mass
        'bdd' : dictionary
            D-set damping
        'kdd' : dictionary
            D-set stiffness
        'pdt' : dictionary
            D-set loads
        'mgg' : dictionary
            G-set mass
        'kgg' : dictionary
            G-set stiffness
        'phg' : dictionary
            G-set mode shape matrix
        'rbg' : dictionary
            G-set rigid-body modes; see also drg output and rbgeom_uset
        'drg' : dictionary
            G-set transpose of rigid-body modes; see also 'rbg' and
            :func:n2y.`rbgeom_uset`.  `drg` = `rbg.T` if both are
            present.
        'pg' : dictionary
            G-set loads

        And any other "extra" matrices that were written to the op4
        file.  Some common extras are:

        'fgravh' : array
            gravity on generalized dof for se 0
        'fgravg' : array
            gravity on G-set physical dof for se 0

    See :func:`rdn2cop2` for a description of what is expected of the
    `op2file`.  The `op4file` is expected to contain certain marker
    matrices.  Scalar SE_START starts each superelement and can be
    followed by any matrices for that superelement.  The end of the
    superelement input is marked by a matrix named LOOP_END.

    See also the Nastran DMAP NAS2CAM.
    """
    if not op4file:
        op4file = op2file + '.op4'
        op2file = op2file + '.op2'

    # read op2 file:
    with OP2(op2file) as o2:
        nas = o2.rdn2cop2()

    # read op4 file:
    from pyNastran.op2.dev.op4 import OP4
    o4 = OP4()
    #op4names, op4vars, *_ = o4.listload(op4file)
    op4names, op4vars = o4.listload(op4file)[:1]

    # loop over superelements:
    j = 0
    for se in nas['selist'][:, 0]:
        if op4names[j] != "se_start":
            raise RuntimeError("matrices are not in understandable"
                               " order.  Expected 'se_start', got "
                               "'{0}'".format(op4names[j]))
        # read all matrices for this se
        j += 1
        while 1:
            name = op4names[j]
            if name in ("loop_end", "se_start"):
                # go on to next se or to residual
                break
            if name not in nas:
                nas[name] = {}
            if se == 0 and name == "lambda":
                # count number of rigid body modes
                nrb = sum(op4vars[j] < .005)[0]
                nas['nrb'] = nrb
                nas['lambda'][0] = abs(op4vars[j].flatten())
            elif name == 'lambda':
                nas[name][se] = op4vars[j].flatten()
            elif name == 'rfmodes':
                nas[name][se] = np.nonzero(op4vars[j])[0]
            else:
                nas[name][se] = op4vars[j]
            j += 1
        if name == "loop_end":
            j += 1
            break
    while j < len(op4vars):
        nas[op4names[j]] = op4vars[j]
        j += 1
    return nas


def get_dof_descs():
    """
    Returns dictionary of descriptions for Nastran data recovery items.
    Normally called by :func:`procdrm12`.

    Returns
    -------
    desc : dictionary
        Has keys: 'acce', 'spcf', 'stress', 'force'

    desc['acce'] : numpy string array
        ['T1', 'T2', 'T3',  'R1', 'R2', 'R3']
    desc['spcf'] : numpy string array
        ['Fx', 'Fy', 'Fz',  'Mx', 'My', 'Mz']
    desc['stress'] : dict
        Dictionary with element numbers as keys to numpy string arrays.
    desc['stress'] : dict
        Dictionary with element numbers as keys to numpy string arrays.

    The stress and force returns are dictionaries indexed by the
    element id.  For example, for the CBAR (which is element 34):

    desc['stress'][34] = ['CBAR Bending Stress 1 - End A',
                          'CBAR Bending Stress 2 - End A',
                          ...]
    desc['force'][34] = ['CBAR Bending Moment 1 - End A',
                         'CBAR Bending Moment 2 - End A',
                         ...]
    """
    #   Acceleration, Velocity, Displacement Recovery Items:
    accedesc = ["T1", "T2", "T3", "R1", "R2", "R3"]
    spcfdesc = ["Fx", "Fy", "Fz", "Mx", "My", "Mz"]
    stress = {}
    force = {}

    #  CBAR Recovery Items (element 34):                        Item code
    stress[34] = ["CBAR Bending Stress 1 - End A",         # 2
                  "CBAR Bending Stress 2 - End A",         # 3
                  "CBAR Bending Stress 3 - End A",         # 4
                  "CBAR Bending Stress 4 - End A",         # 5
                  "CBAR Axial Stress",                     # 6
                  "CBAR Max. Bend. Stress -End A",         # 7
                  "CBAR Min. Bend. Stress -End A",         # 8
                  "CBAR M.S. Tension",                     # 9
                  "CBAR Bending Stress 1 - End B",         # 10
                  "CBAR Bending Stress 2 - End B",         # 11
                  "CBAR Bending Stress 3 - End B",         # 12
                  "CBAR Bending Stress 4 - End B",         # 13
                  "CBAR Max. Bend. Stress -End B",         # 14
                  "CBAR Min. Bend. Stress -End B",         # 15
                  "CBAR M.S. Compression"]                 # 16

    force[34] = ["CBAR Bending Moment 1 - End A",         # 2
                 "CBAR Bending Moment 2 - End A",         # 3
                 "CBAR Bending Moment 1 - End B",         # 4
                 "CBAR Bending Moment 2 - End B",         # 5
                 "CBAR Shear 1",                          # 6
                 "CBAR Shear 2",                          # 7
                 "CBAR Axial Force",                      # 8
                 "CBAR Torque"]                           # 9

    #   CBEAM Recovery Items (element 2):                        Item code
    stress2_main = ["CBEAM External grid pt. ID",          # 2
                    "CBEAM Station dist./length",          # 3
                    "CBEAM Long. Stress at Pt. C",         # 4
                    "CBEAM Long. Stress at Pt. D",         # 5
                    "CBEAM Long. Stress at Pt. E",         # 6
                    "CBEAM Long. Stress at Pt. F",         # 7
                    "CBEAM Maximum stress",                # 8
                    "CBEAM Minimum stress",                # 9
                    "CBEAM M.S. Tension",                  # 10
                    "CBEAM M.S. Compression"]              # 11

    # expand and append station id for all 11 stations:
    stress2 = [i+' End-A' for i in stress2_main]
    for K in range(2, 11):
        id_string = ' K={0:2}'.format(K)
        stress2 += [i+id_string for i in stress2_main]
    stress2 += [i+' End-B' for i in stress2_main]
    stress[2] = stress2

    force2_main = ["CBEAM External grid pt. ID",             # 2
                   "CBEAM Station dist./length",             # 3
                   "CBEAM Bending moment plane 1",           # 4
                   "CBEAM Bending moment plane 2",           # 5
                   "CBEAM Web shear plane 1",                # 6
                   "CBEAM Web shear plane 2",                # 7
                   "CBEAM Axial force",                      # 8
                   "CBEAM Total torque",                     # 9
                   "CBEAM Warping torque"]                   # 10

    # expand and append station id for all 11 stations:
    force2 = [i+' End-A' for i in force2_main]
    for K in range(2, 11):
        id_string = ' K={0:2}'.format(K)
        force2 += [i+id_string for i in force2_main]
    force2 += [i+' End-B' for i in force2_main]
    force[2] = force2

    #   CBUSH Recovery Items (element 102):                        Item code
    stress[102] = ["CBUSH Translation-x",         # 2
                   "CBUSH Translation-y",         # 3
                   "CBUSH Translation-z",         # 4
                   "CBUSH Rotation-x",            # 5
                   "CBUSH Rotation-y",            # 6
                   "CBUSH Rotation-z"]            # 7

    force[102] = ["CBUSH Force-x",          # 2
                  "CBUSH Force-y",          # 3
                  "CBUSH Force-z",          # 4
                  "CBUSH Moment-x",         # 5
                  "CBUSH Moment-y",         # 6
                  "CBUSH Moment-z"]         # 7

    #   CROD Recovery Items (element 10=CONROD, 1=CROD):
    stress1 = ["Axial Stress",              # 2
               "M.S. Axial Stress",         # 3
               "Torsional Stress",          # 4
               "M.S. Torsional Stress"]     # 5
    force1 = ["Axial Force",        # 2
              "Torque"]             # 3
    stress[1] = ['CROD '+ i + '  ' for i in stress1]
    force[1] = ['CROD '+ i + '  ' for i in force1]
    stress[10] = ['CONROD ' + i for i in stress1]
    force[10] = ['CONROD ' + i for i in force1]

    #   CELAS1, 2, 3 Recovery Items (elements 11, 12, 13):
    stress[11] = 'CELAS1 Stress'
    stress[12] = 'CELAS2 Stress'
    stress[13] = 'CELAS3 Stress'
    force[11] = 'CELAS1 Force'
    force[12] = 'CELAS2 Force'
    force[13] = 'CELAS3 Force'

    #   CQUAD4 Recovery Items (element 33):
    stress[33] = ["CQUAD4 Fiber distance Z1",           # 2
                  "CQUAD4 Z1 Normal x",                 # 3
                  "CQUAD4 Z1 Normal y",                 # 4
                  "CQUAD4 Z1 Shear xy",                 # 5
                  "CQUAD4 Z1 Shear angle",              # 6
                  "CQUAD4 Z1 Major principal",          # 7
                  "CQUAD4 Z1 Minor principal",          # 8
                  "CQUAD4 Z1 von Mises or max shear",   # 9
                  "CQUAD4 Fiber distance Z2",           # 10
                  "CQUAD4 Z2 Normal x",                 # 11
                  "CQUAD4 Z2 Normal y",                 # 12
                  "CQUAD4 Z2 Shear xy",                 # 13
                  "CQUAD4 Z2 Shear angle",              # 14
                  "CQUAD4 Z2 Major principal",          # 15
                  "CQUAD4 Z2 Minor principal",          # 16
                  "CQUAD4 Z2 von Mises or max shear"]   # 17

    force[33] = ["CQUAD4 Membrane force x",         # 2
                 "CQUAD4 Membrane force y",         # 3
                 "CQUAD4 Membrane force xy",        # 4
                 "CQUAD4 Bending moment x",         # 5
                 "CQUAD4 Bending moment y",         # 6
                 "CQUAD4 Bending moment xy",        # 7
                 "CQUAD4 Shear x",                  # 8
                 "CQUAD4 Shear y"]                  # 9

    #   CQUADR Recovery Items (element 82, and CQUAD8-64):
    stress[82] = ["CQUADR EID                         ",      # 1
                  "CQUADR CEN/                        ",      # 2
                  "CQUADR 4                           ",      # 3
                  "CQUADR Fiber distance Z1           ",      # 4
                  "CQUADR Z1 Normal x                 ",      # 5
                  "CQUADR Z1 Normal y                 ",      # 6
                  "CQUADR Z1 Shear xy                 ",      # 7
                  "CQUADR Z1 Shear angle              ",      # 8
                  "CQUADR Z1 Major principal          ",      # 9
                  "CQUADR Z1 Minor principal          ",      # 10
                  "CQUADR Z1 von Mises or max shear   ",      # 11
                  "CQUADR Fiber distance Z2           ",      # 12
                  "CQUADR Z2 Normal x                 ",      # 13
                  "CQUADR Z2 Normal y                 ",      # 14
                  "CQUADR Z2 Shear xy                 ",      # 15
                  "CQUADR Z2 Shear angle              ",      # 16
                  "CQUADR Z2 Major principal          ",      # 17
                  "CQUADR Z2 Minor principal          ",      # 18
                  "CQUADR Z2 von Mises or max shear   ",      # 19

                  "CQUADR Grid 1                      ",      # 20
                  "CQUADR Fiber distance Z1         c1",      # 21
                  "CQUADR Z1 Normal x               c1",      # 22
                  "CQUADR Z1 Normal y               c1",      # 23
                  "CQUADR Z1 Shear xy               c1",      # 24
                  "CQUADR Z1 Shear angle            c1",      # 25
                  "CQUADR Z1 Major principal        c1",      # 26
                  "CQUADR Z1 Minor principal        c1",      # 27
                  "CQUADR Z1 von Mises or max shear c1",      # 28
                  "CQUADR Fiber distance Z2         c1",      # 29
                  "CQUADR Z2 Normal x               c1",      # 30
                  "CQUADR Z2 Normal y               c1",      # 31
                  "CQUADR Z2 Shear xy               c1",      # 32
                  "CQUADR Z2 Shear angle            c1",      # 33
                  "CQUADR Z2 Major principal        c1",      # 34
                  "CQUADR Z2 Minor principal        c1",      # 35
                  "CQUADR Z2 von Mises or max shear c1",      # 36

                  "CQUADR Grid 2                      ",      # 37
                  "CQUADR Fiber distance Z1         c2",      # 38
                  "CQUADR Z1 Normal x               c2",      # 39
                  "CQUADR Z1 Normal y               c2",      # 40
                  "CQUADR Z1 Shear xy               c2",      # 41
                  "CQUADR Z1 Shear angle            c2",      # 42
                  "CQUADR Z1 Major principal        c2",      # 43
                  "CQUADR Z1 Minor principal        c2",      # 44
                  "CQUADR Z1 von Mises or max shear c2",      # 45
                  "CQUADR Fiber distance Z2         c2",      # 46
                  "CQUADR Z2 Normal x               c2",      # 47
                  "CQUADR Z2 Normal y               c2",      # 48
                  "CQUADR Z2 Shear xy               c2",      # 49
                  "CQUADR Z2 Shear angle            c2",      # 50
                  "CQUADR Z2 Major principal        c2",      # 51
                  "CQUADR Z2 Minor principal        c2",      # 52
                  "CQUADR Z2 von Mises or max shear c2",      # 53

                  "CQUADR Grid 3                      ",      # 54
                  "CQUADR Fiber distance Z1         c3",      # 55
                  "CQUADR Z1 Normal x               c3",      # 56
                  "CQUADR Z1 Normal y               c3",      # 57
                  "CQUADR Z1 Shear xy               c3",      # 58
                  "CQUADR Z1 Shear angle            c3",      # 59
                  "CQUADR Z1 Major principal        c3",      # 60
                  "CQUADR Z1 Minor principal        c3",      # 61
                  "CQUADR Z1 von Mises or max shear c3",      # 62
                  "CQUADR Fiber distance Z2         c3",      # 63
                  "CQUADR Z2 Normal x               c3",      # 64
                  "CQUADR Z2 Normal y               c3",      # 65
                  "CQUADR Z2 Shear xy               c3",      # 66
                  "CQUADR Z2 Shear angle            c3",      # 67
                  "CQUADR Z2 Major principal        c3",      # 68
                  "CQUADR Z2 Minor principal        c3",      # 69
                  "CQUADR Z2 von Mises or max shear c3",      # 70

                  "CQUADR Grid 4                      ",      # 71
                  "CQUADR Fiber distance Z1         c4",      # 72
                  "CQUADR Z1 Normal x               c4",      # 73
                  "CQUADR Z1 Normal y               c4",      # 74
                  "CQUADR Z1 Shear xy               c4",      # 75
                  "CQUADR Z1 Shear angle            c4",      # 76
                  "CQUADR Z1 Major principal        c4",      # 77
                  "CQUADR Z1 Minor principal        c4",      # 78
                  "CQUADR Z1 von Mises or max shear c4",      # 79
                  "CQUADR Fiber distance Z2         c4",      # 80
                  "CQUADR Z2 Normal x               c4",      # 81
                  "CQUADR Z2 Normal y               c4",      # 82
                  "CQUADR Z2 Shear xy               c4",      # 83
                  "CQUADR Z2 Shear angle            c4",      # 84
                  "CQUADR Z2 Major principal        c4",      # 85
                  "CQUADR Z2 Minor principal        c4",      # 86
                  "CQUADR Z2 von Mises or max shear c4"]      # 87

    force[82] = ["CQUADR Membrane force x            ",      # 4
                 "CQUADR Membrane force y            ",      # 5
                 "CQUADR Membrane force xy           ",      # 6
                 "CQUADR Bending moment x            ",      # 7
                 "CQUADR Bending moment y            ",      # 8
                 "CQUADR Bending moment xy           ",      # 9
                 "CQUADR Shear x                     ",      # 10
                 "CQUADR Shear y                     ",      # 11

                 "CQUADR   (non-documented item)     ",      # 12

                 "CQUADR Membrane force x          c1",      # 13
                 "CQUADR Membrane force y          c1",      # 14
                 "CQUADR Membrane force xy         c1",      # 15
                 "CQUADR Bending moment x          c1",      # 16
                 "CQUADR Bending moment y          c1",      # 17
                 "CQUADR Bending moment xy         c1",      # 18
                 "CQUADR Shear x                   c1",      # 19
                 "CQUADR Shear y                   c1",      # 20

                 "CQUADR   (non-documented item)     ",      # 21

                 "CQUADR Membrane force x          c2",      # 22
                 "CQUADR Membrane force y          c2",      # 23
                 "CQUADR Membrane force xy         c2",      # 24
                 "CQUADR Bending moment x          c2",      # 25
                 "CQUADR Bending moment y          c2",      # 26
                 "CQUADR Bending moment xy         c2",      # 27
                 "CQUADR Shear x                   c2",      # 28
                 "CQUADR Shear y                   c2",      # 29

                 "CQUADR   (non-documented item)     ",      # 30

                 "CQUADR Membrane force x          c3",      # 31
                 "CQUADR Membrane force y          c3",      # 32
                 "CQUADR Membrane force xy         c3",      # 33
                 "CQUADR Bending moment x          c3",      # 34
                 "CQUADR Bending moment y          c3",      # 35
                 "CQUADR Bending moment xy         c3",      # 36
                 "CQUADR Shear x                   c3",      # 37
                 "CQUADR Shear y                   c3",      # 38

                 "CQUADR   (non-documented item)     ",      # 39

                 "CQUADR Membrane force x          c4",      # 40
                 "CQUADR Membrane force y          c4",      # 41
                 "CQUADR Membrane force xy         c4",      # 42
                 "CQUADR Bending moment x          c4",      # 43
                 "CQUADR Bending moment y          c4",      # 44
                 "CQUADR Bending moment xy         c4",      # 45
                 "CQUADR Shear x                   c4",      # 46
                 "CQUADR Shear y                   c4"]      # 47
    stress[64] = [i.replace('CQUADR', 'CQ8-64') for i in stress[82]]
    force[64] = [i.replace('CQUADR', 'CQ8-64') for i in force[82]]

    #   CTRIAR Recovery Items (element 70, and CTRIA6-75):
    stress[70] = ["CTRIAR Z1 Normal x                 ",       # 5
                  "CTRIAR Z1 Normal y                 ",       # 6
                  "CTRIAR Z1 Shear xy                 ",       # 7
                  "CTRIAR Z1 Q shear angle            ",       # 8
                  "CTRIAR Z1 Major principal          ",       # 9
                  "CTRIAR Z1 Minor principal          ",       # 10
                  "CTRIAR Z1 von Mises or max shear   ",       # 11
                  "CTRIAR   (non-documented item)     ",       # 12
                  "CTRIAR Z2 Normal x                 ",       # 13
                  "CTRIAR Z2 Normal y                 ",       # 14
                  "CTRIAR Z2 Shear xy                 ",       # 15
                  "CTRIAR Z2 Q shear angle            ",       # 16
                  "CTRIAR Z2 Major principal          ",       # 17
                  "CTRIAR Z2 Minor principal          ",       # 18
                  "CTRIAR Z2 von Mises or max shear   ",       # 19

                  "CTRIAR   (non-documented item)     ",       # 20
                  "CTRIAR   (non-documented item)     ",       # 21

                  "CTRIAR Z1 Normal x               c1",       # 22
                  "CTRIAR Z1 Normal y               c1",       # 23
                  "CTRIAR Z1 Shear xy               c1",       # 24
                  "CTRIAR Z1 Q shear angle          c1",       # 25
                  "CTRIAR Z1 Major principal        c1",       # 26
                  "CTRIAR Z1 Minor principal        c1",       # 27
                  "CTRIAR Z1 von Mises or max shear c1",       # 28
                  "CTRIAR   (non-documented item)   c1",       # 29
                  "CTRIAR Z2 Normal x               c1",       # 30
                  "CTRIAR Z2 Normal y               c1",       # 31
                  "CTRIAR Z2 Shear xy               c1",       # 32
                  "CTRIAR Z2 Q shear angle          c1",       # 33
                  "CTRIAR Z2 Major principal        c1",       # 34
                  "CTRIAR Z2 Minor principal        c1",       # 35
                  "CTRIAR Z2 von Mises or max shear c1",       # 36

                  "CTRIAR   (non-documented item)     ",       # 37
                  "CTRIAR   (non-documented item)     ",       # 38

                  "CTRIAR Z1 Normal x               c2",       # 39
                  "CTRIAR Z1 Normal y               c2",       # 40
                  "CTRIAR Z1 Shear xy               c2",       # 41
                  "CTRIAR Z1 Q shear angle          c2",       # 42
                  "CTRIAR Z1 Major principal        c2",       # 43
                  "CTRIAR Z1 Minor principal        c2",       # 44
                  "CTRIAR Z1 von Mises or max shear c2",       # 45
                  "CTRIAR   (non-documented item)   c2",       # 46
                  "CTRIAR Z2 Normal x               c2",       # 47
                  "CTRIAR Z2 Normal y               c2",       # 48
                  "CTRIAR Z2 Shear xy               c2",       # 49
                  "CTRIAR Z2 Q shear angle          c2",       # 50
                  "CTRIAR Z2 Major principal        c2",       # 51
                  "CTRIAR Z2 Minor principal        c2",       # 52
                  "CTRIAR Z2 von Mises or max shear c2",       # 53

                  "CTRIAR   (non-documented item)     ",       # 54
                  "CTRIAR   (non-documented item)     ",       # 55

                  "CTRIAR Z1 Normal x               c3",       # 56
                  "CTRIAR Z1 Normal y               c3",       # 57
                  "CTRIAR Z1 Shear xy               c3",       # 58
                  "CTRIAR Z1 Q shear angle          c3",       # 59
                  "CTRIAR Z1 Major principal        c3",       # 60
                  "CTRIAR Z1 Minor principal        c3",       # 61
                  "CTRIAR Z1 von Mises or max shear c3",       # 62
                  "CTRIAR   (non-documented item)   c3",       # 63
                  "CTRIAR Z2 Normal x               c3",       # 64
                  "CTRIAR Z2 Normal y               c3",       # 65
                  "CTRIAR Z2 Shear xy               c3",       # 66
                  "CTRIAR Z2 Q shear angle          c3",       # 67
                  "CTRIAR Z2 Major principal        c3",       # 68
                  "CTRIAR Z2 Minor principal        c3",       # 69
                  "CTRIAR Z2 von Mises or max shear c3"]       # 70

    force[70] = ["CTRIAR Membrane force x            ",      # 4
                 "CTRIAR Membrane force y            ",      # 5
                 "CTRIAR Membrane force xy           ",      # 6
                 "CTRIAR Bending moment x            ",      # 7
                 "CTRIAR Bending moment y            ",      # 8
                 "CTRIAR Bending moment xy           ",      # 9
                 "CTRIAR Shear x                     ",      # 10
                 "CTRIAR Shear y                     ",      # 11

                 "CTRIAR   (non-documented item)     ",      # 12

                 "CTRIAR Membrane force x          c1",      # 13
                 "CTRIAR Membrane force y          c1",      # 14
                 "CTRIAR Membrane force xy         c1",      # 15
                 "CTRIAR Bending moment x          c1",      # 16
                 "CTRIAR Bending moment y          c1",      # 17
                 "CTRIAR Bending moment xy         c1",      # 18
                 "CTRIAR Shear x                   c1",      # 19
                 "CTRIAR Shear y                   c1",      # 20

                 "CTRIAR   (non-documented item)     ",      # 21

                 "CTRIAR Membrane force x          c2",      # 22
                 "CTRIAR Membrane force y          c2",      # 23
                 "CTRIAR Membrane force xy         c2",      # 24
                 "CTRIAR Bending moment x          c2",      # 25
                 "CTRIAR Bending moment y          c2",      # 26
                 "CTRIAR Bending moment xy         c2",      # 27
                 "CTRIAR Shear x                   c2",      # 28
                 "CTRIAR Shear y                   c2",      # 29

                 "CTRIAR   (non-documented item)     ",      # 30

                 "CTRIAR Membrane force x          c3",      # 31
                 "CTRIAR Membrane force y          c3",      # 32
                 "CTRIAR Membrane force xy         c3",      # 33
                 "CTRIAR Bending moment x          c3",      # 34
                 "CTRIAR Bending moment y          c3",      # 35
                 "CTRIAR Bending moment xy         c3",      # 36
                 "CTRIAR Shear x                   c3",      # 37
                 "CTRIAR Shear y                   c3"]      # 38

    stress[75] = [i.replace('CTRIAR', 'CT6-75') for i in stress[70]]
    force[75] = [i.replace('CTRIAR', 'CT6-75') for i in force[70]]
    for i in stress:
        stress[i] = np.array(stress[i])
        force[i] = np.array(force[i])
    return {'acce': np.array(accedesc),
            'spcf': np.array(spcfdesc),
            'stress': stress,
            'force': force}


def _get_tinr(iddof, idj):
    """
    Called by get_drm.

    Parameters
    ----------
    iddof : 2d array
        Each col has [type, id, number of rows, start row]
    idj : integer
        Id to return info for.

    Returns tuple of (type, start row)

    .. note:: start row return value starts at 0, not at 1.
    """
    i = np.nonzero(iddof[1] == idj)[0]
    tinr = iddof[:, i]
    return tinr[0, 0], tinr[3, 0]-1


def get_drm(drminfo, otm, drms, drmkeys, dr, desc):
    """
    Called by :func:`procdrm12` to add displacement-dependent data
    recovery items to the otm input.

    Parameters
    ----------
    drminfo : tuple
        DRM Information; (output drm name, 3 or 5 character Nastran
        name, description index).

        - if the second input is 3 chars, say '---', this routine
          uses the following members of `drms` and `drmkeys`::
            'm---d1', 'm---s1' and 't---d1' if available (mode-acce), or
            'm---x1', 't---x1' if not (mode-disp)
        - if the second input is 5 chars, say '-----', this routine
          uses 'm-----' and 't-----'
        - the description index is used to get info from `desc`.
    otm : input/output dictionary
        Filled in with 'DTM' (or 'DTMA', 'DTMD') and 'DTM_id_dof',
        'DTM_desc'.
    drms : dictionary
        Contains all drms from op4 file.
    drmkeys : dictionary
        Contains the keys (directories) to the drms.
    dr : array
        Matrix 3 x number of data recovery items: [type; id; dof].
        Type is 1 for displacements.
    desc : dictionary
        Output of :func:`get_dof_descs`.

    Examples usages::

        get_drm(('DTM', 'oug', 'acce'), otm, drms, drmkeys, dr, desc)
        get_drm(('ATM', 'ougv1', 'acce'), ...)
        get_drm(('LTM', 'oef', 'force'), ...)
        get_drm(('SPCF', 'oqg', 'spcf'), ...)
        get_drm(('STM', 'oes', 'stress'), ...)

    """
    drc = dr.shape[1]
    ID = dr[1, :]
    DOF = dr[2, :]
    nm, nasnm, desci = drminfo
    otm[nm+'_id_dof'] = np.vstack((ID, DOF)).T

    # arg offset is for translating between Nastran argument to
    # matrix index; eg 'x' recovery for a grid is arg 3, so offset
    # is 3
    if nasnm.find('oug') > -1 or nasnm.find('oqg') > -1:
        offset = 3
        otm[nm+'_id_dof'][:, 1] -= 2
    else:
        offset = 2

    if not isinstance(desc[desci], dict):
        otm[nm+'_desc'] = desc[desci][DOF-offset]
        getdesc = False
    else:
        getdesc = True
        _desc = nm+'_desc'
        otm[_desc] = [''] * drc
        _dct = desc[desci]
        _name = desci.capitalize()

    if len(nasnm) == 3 and 'm'+nasnm+'d1' in drms:
        d1 = drms['m'+nasnm+'d1'][0]
        s1 = drms['m'+nasnm+'s1'][0]
        iddof = drmkeys['t'+nasnm+'d1']
        acce = nm+'A'
        disp = nm+'D'
        otm[acce] = np.zeros((drc, d1.shape[1]))
        otm[disp] = np.zeros((drc, s1.shape[1]))
        lastid = -1
        for j in range(drc):  # loop over requests
            # find rows corresponding to requested grid
            if ID[j] != lastid:
                eltype, srow = _get_tinr(iddof, ID[j])
                lastid = ID[j]
            otm[acce][j] = d1[srow+DOF[j]-offset]
            otm[disp][j] = s1[srow+DOF[j]-offset]
            if getdesc:
                if eltype in _dct:
                    otm[_desc][j] = _dct[eltype][DOF[j]-offset]
                else:
                    otm[_desc][j] = ('EL-{0}, El. Type {1:3}, '
                                     'Code {2:3}  ').format(_name,
                                                            eltype,
                                                            DOF[j])
    else:
        if len(nasnm) == 3:
            matname = 'm'+nasnm+'x1'
            tabname = 't'+nasnm+'x1'
        else:
            matname = 'm'+nasnm
            tabname = 't'+nasnm
        x1 = drms[matname][0]
        iddof = drmkeys[tabname]
        otm[nm] = np.zeros((drc, x1.shape[1]))
        lastid = -1
        for j in range(drc):  # loop over requests
            # find rows corresponding to requested grid
            if ID[j] != lastid:
                eltype, srow = _get_tinr(iddof, ID[j])
                lastid = ID[j]
            otm[nm][j] = x1[srow+DOF[j]-offset]
            if getdesc:
                if eltype in _dct:
                    otm[_desc][j] = _dct[eltype][DOF[j]-offset]
                else:
                    otm[_desc][j] = ('EL-{0}, El. Type {1:3}, '
                                     'Code {2:3}  ').format(_name,
                                                            eltype,
                                                            DOF[j])


def proccess_drm1_drm2(op2file, op4file=None, dosort=True):
    """
    Process op2/op4 file2 output from DRM1/DRM2 DMAPs to form data
    recovery matrices.

    Parameters
    ----------
    op2file : string
        Either the basename of the .op2 and .op4 files, or the full
        name of the .op2 file
    op4file : string or None
        The name of the .op4 file or, if None, builds name from the
        `op2file` input.
    dosort : bool
        If True, sort data recovery rows in ascending order by ID/DOF.
        Otherwise, return in order requested in Nastran run.

    Returns
    -------
    otm : dictionary
        Has data recovery matrices (DRMs), id/dof info, and generic
        descriptions.  The potential DRM keys are:
        ::
            'ATM'  : acceleration DRM

          For mode-displacement:
            'DTM'  : displacement DRM
            'LTM'  : element force (loads) DRM
            'SPCF' : SPC forces DRM
            'STM'  : element stress DRM

          For mode-acceleration:
            'DTMD' : displacement-dependent part of displacement DRM
            'DTMA' : acceleration-dependent part of displacement DRM
            'LTMD' : displacement-dependent part of element force DRM
            'LTMA' : acceleration-dependent part of element force DRM
            'SPCFD': displacement-dependent part of SPCF forces DRM
            'SPCFA': acceleration-dependent part of SPCF forces DRM
            'STMD' : displacement-dependent part of element stress DRM
            'STMA' : displacement-dependent part of element stress DRM

        The id/dof matrices are each 2 columns of [id, dof] with number
        of rows equal to the number of rows in corresponding DRM.  The
        keys are the applicable strings from:
        ::
            'ATM_id_dof'
            'DTM_id_dof'
            'LTM_id_dof' - dof is actually the Nastran item code
            'SPCF_id_dof'
            'STM_id_dof' - dof is actually the Nastran item code

        The descriptions are arrays of strings with generic descriptions
        for each data recovery item.  Length is equal to number of rows
        in corresponding DRM. See :func:`get_dof_descs` for more
        information.  The keys are the applicable strings from:
        ::
            'ATM_desc'
            'DTM_desc'
            'LTM_desc',
            'SPCF_desc'
            'STM_desc'.

    Currently, only displacements, accelerations, SPC forces, element
    forces and element stresses (for some elements) are implemented.

    Example usage::

        import op2
        otm = op2.proccess_drm1_drm2('drm2')

    """
    if not op4file:
        op4file = op2file + '.op4'
        op2file = op2file + '.op2'

    # read op4 file:
    from pyNastran.op2.dev.op4 import OP4
    o4 = OP4()
    drms = o4.dctload(op4file)

    with OP2(op2file) as o2:
        drm_keys = o2.rddrm2op2()
    N = drm_keys['drs'].shape[1]

    # drs format:
    # 6 elements per recovery item:
    #    1  -  Subcase number (0 for all)
    #    2  -  Vector request type
    #    3  -  Point or Element ID
    #    4  -  Component
    #    5  -  XY output type
    #    6  -  Destination code

    # Vector request type:
    Vreq = ["Displacement",     # 1
            "Velocity",         # 2
            "Acceleration",     # 3
            "SPC Force",        # 4
            "Load",             # 5
            "Stress",           # 6
            "Element Force",    # 7
            "SDisplacement",    # 8
            "SVelocity",        # 9
            "SAcceleration",    # 10
            "Nonlinear Force",  # 11
            "Total"]            # 12

    #   XY output type:
    #      1 = Response
    #      2 = PSDF
    #      3 = AUTO
    #
    #   Destination code:
    #      0 = XYpeak only   (from DRMEXT)
    #      1 = Print
    #      2 = Plot
    #      3 = Print, Plot
    #      4 = Punch
    #      5 = Print, Punch
    #      6 = Plot, Punch
    #      7 = Print, Plot, Punch

    if not dosort:
        # reshape dr:
        dr = drm_keys['dr']
        r = np.nonzero(dr == dr[0])[0]
        r = np.hstack((r, len(dr)))
        n = len(r) - 1
        # dr(r) = ? -- starts every XYPEAK card
        # dr(r+1:3) = 0, 0, 0  ?
        # dr(r+4) = 1  ?
        # dr(r+5) = request type
        # dr(r+6:8) = 0, 0, #(?)
        # dr(r+9) = id 1
        # dr(r+10) = dof 1
        # dr(r+11) = 0
        #     ... r + 9, 10, 11 can repeat until three -1's are reached
        #        These 3 values repeat when there is a comma: 1(T1),1(T2)
        # dr(X:X+2) = -1, -1, -1
        # 8-X+2 repeat until all dof for an XYPEAK are listed
        #        This section repeats when there is a slash: 1(T1)/1(T2)
        DR = np.zeros((3, N), dtype=int)  # [type; id; dof]
        R = 0  # index into DR columns
        for j in range(n):  # loop over XYPEAK cards
            curtype = dr[r[j] + 5]
            J = r[j] + 9  # index to first id
            while J < r[j+1]:
                while dr[J] != -1:
                    DR[:, R] = curtype, dr[J], dr[J+1]
                    R += 1
                    J += 3
                J += 4  # jump over [-1,-1,-1,#]
    else:
        DR = drm_keys['drs'][1:4]  # use sorted version

    desc = get_dof_descs()
    drm_info = {
        1: ('DTM', 'oug', 'acce'),
        3: ('ATM', 'ougv1', 'acce'),
        4: ('SPCF', 'oqg', 'spcf'),
        6: ('STM', 'oes', 'stress'),
        7: ('LTM', 'oef', 'force'),
    }
    otm = {}
    types = np.array([1, 3, 4, 6, 7])
    for drtype in range(1, 13):
        pv = np.nonzero(DR[0] == drtype)[0]
        if pv.size > 0:
            if np.any(drtype == types):
                print('Processing "{0}" requests...'.format(Vreq[drtype-1]))
                get_drm(drm_info[drtype], otm, drms,
                        drm_keys, DR[:, pv], desc)
            else:
                print('Skipping %r requests.  Needs to be added '
                      'to proccess_drm1_drm2().' % Vreq[drtype-1])
    return otm


def read_post_op2(op2_filename, verbose=False, getougv1=False):
    """
    Reads PARAM,POST,-1 op2 file and returns dictionary of data.

    Parameters
    ----------
    op2_filename : string
        Name of op2 file.
    verbose : bool
        If true, echo names of tables and matrices to screen
    getougv1 : bool
        If true, read the OUGV1 matrices, if any.

    Returns dictionary with following members
    -----------------------------------------
    'uset' : array
        6-column matrix as described in class OP2, member function
        :func:`readd_nas2cam_op2`.
    'cstm' : array
        14-column matrix containing the coordinate system
        transformation matrix for each coordinate system.  See
        description in class OP2, member function
        :func:`readd_nas2cam_op2`.
    'cstm2' : dictionary
        Dictionary indexed by the coordinate system id number.  This
        has the same information as 'cstm', but in a different format.
        See description in class OP2, member function
        :func:`readd_nas2cam_op2`.
    'mats' : dictionary
        Dictionary of matrices read from op2 file and indexed by the
        name.  The 'tload' entry is a typical entry.  If `getougv1` is
        true, `mats` will contain a list of all 'OUGV1' and 'BOPHIG'
        matrices.
    """
    # read op2 file:
    with OP2(op2_filename) as o2:
        mats = {}
        selist = uset = cstm2 = None
        se = 0
        if getougv1:
            mats['ougv1'] = []
        o2._fileh.seek(o2._postheaderpos)

        eqexin1 = None
        dof = None
        Uset = None
        cstm = None
        while 1:
            name, trailer, dbtype = o2._read_op2_name_trailer()
             #print('name = %r' % name)
             #print('trailer = %s' % str(trailer))
             #print('dbtype = %r' % dbtype)
            if name is None:
                break
            if name == '':
                raise RuntimeError('name=%r' % name)
            if dbtype > 0:
                if verbose:
                    print("Reading matrix {0}...".format(name))
                if name not in mats:
                    mats[name] = []
                mats[name] += [o2.read_op2_matrix(name, trailer)]
            else:
                if name.find('BGPDT') == 0:
                    if verbose:
                        print("Reading table {0}...".format(name))
                    bgpdt_rec1 = o2._read_op2_bgpdt68()
                    o2.skip_op2_table()
                    continue

                # if name.find('CSTM') == 0:
                #     if verbose:
                #         print("Reading table {}...".format(name))
                #     cstm = o2._rdop2cstm68().reshape((-1, 14))
                #     cstm = np.vstack((bc, cstm))
                #     continue

                elif name.find('GEOM1') == 0:
                    if verbose:
                        print("Reading table {0}...".format(name))
                    cords, unused_sebulk, selist = o2._read_op2_geom1_cord2()
                    if 0 not in cords:
                        cords[0] = np.array([[0., 1., 0.],
                                             [0., 0., 0.],
                                             [1., 0., 0.],
                                             [0., 1., 0.],
                                             [0., 0., 1.]])
                    if -1 not in cords:
                        cords[-1] = np.zeros((5, 3))  # dummy for spoints
                        cords[-1][0, 0] = -1
                    cstm2 = cords
                    continue

                elif name.find('DYNAMIC') == 0:
                    if verbose:
                        print("Reading DYNAMIC table {0}...".format(name))
                    mats['tload'] = o2.read_op2_dynamics()
                    continue

                elif name.find('EQEXIN') == 0:
                    if verbose:
                        print("Reading EQEXIN table {0}...".format(name))
                    eqexin1, eqexin = o2._read_op2_eqexin()
                    continue

                elif name.find('USET') == 0:
                    if verbose:
                        print("Reading USET table {0}...".format(name))
                    uset = o2._read_op2_uset()
                    continue

                elif getougv1 and (name.find('OUGV1') == 0 or
                                   name.find('BOPHIG') == 0):
                    if verbose:
                        print("Reading OUG table {0}...".format(name))
                    mats['ougv1'] += [o2._read_op2_ougv1(name)]
                    continue

                # if name.find('OEF1X') == 0:
                #    if verbose:
                #        print("Reading table {}...\n".format(name))
                #     mats['oef1x'] = o2._rdop2drm()
                #     continue

                elif verbose:
                    print("Skipping table %r..." % name)
                o2.skip_op2_table()

        if eqexin1 is not None:
            (bgpdt, dof, doftype, nid, upids) = o2._proc_bgpdt(
                eqexin1, eqexin, True, bgpdt_rec1)
        if dof is not None:
            Uset, cstm, cstm2 = o2._build_Uset(
                se, dof, doftype, nid, uset, bgpdt, None, cstm2)
    return {'uset': Uset,
            'cstm': cstm,
            'cstm2': cstm2,
            'mats': mats,
            'selist': selist}
