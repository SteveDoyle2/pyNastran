"""Main OP4 class"""
import sys
import os
from struct import pack, unpack, Struct

import numpy as np
from numpy import array, zeros, float32, float64, complex64, complex128, ndarray
from scipy.sparse import coo_matrix  # type: ignore
from cpylog import get_logger2

from pyNastran.utils import is_binary_file as file_is_binary
from pyNastran.utils.mathematics import print_matrix #, print_annotated_matrix


def read_op4(op4_filename=None, matrix_names=None, precision='default',
             debug=False, log=None):
    """
    Reads a NASTRAN OUTPUT4 file, and stores the
    matrices as the output arguments.  The number of
    matrices read is defined by the list matrix_names.  By default, all
    matrices will be read.  The resulting output is a dictionary of
    matrices that are accessed by their name.

    .. code-block:: python

       >>> from pyNastran.op4.op4 import OP4
       >>> op4 = OP4()

       # get all the matrices
       >>> matrices = op4.read_op4(op4_filename)
       >>> (formA, A) = matrices['A']
       >>> (formB, B) = matrices['B']
       >>> (formC, C) = matrices['C']

       # or to reduce memory usage
       >>> matrices = op4.read_op4(op4_filename, matrix_names=['A', 'B'])
       >>> (formA, A) = matrices['A']
       >>> (formB, B) = matrices['B']

       # or because you only want A
       >>> matrices = op4.read_op4(op4_filename, matrix_names='A')
       >>> (formA, A) = matrices['A']

       # get all the matrices, but select the file using a file dialog
       >>> matrices = op4.read_op4()
       >>>

    Parameters
    ----------
    op4_filename : str / None
        an OP4 filename.  Type=STRING.
    matrix_names : List[str], str / None
        matrix name(s) (None -> all)
    precision : str; {'default', 'single', 'double'}
        specifies if the matrices are in single or double precsion
        which means the format will be whatever the file is in

    Returns
    -------
    matricies : dict[str] = (int, matrix)
        dictionary of matrices where the key is the name and the value is [form, matrix]

        +------+----------------+
        | Form |   Definition   |
        +======+================+
        |  1   | Square         |
        +------+----------------+
        |  2   | Rectangular    |
        +------+----------------+
        |  3   | Diagonal       |
        +------+----------------+
        |  6   | Symmetric      |
        +------+----------------+
        |  8   | Id entity      |
        +------+----------------+
        |  9   | Pseudoidentity |
        +------+----------------+

        +--------+-------------------------+
        |  Type  | Object                  |
        +========+=========================+
        | Dense  | NUMPY.NDARRAY           |
        +--------+-------------------------+
        | Sparse | SCIPY.SPARSE.COO_MATRIX |
        +--------+-------------------------+

    .. note:: based off the MATLAB code SAVEOP4 developed by ATA-E and
              later UCSD.
    .. note:: it's strongly recommended that you convert sparse matrices to
              another format before doing math on them.  This is standard
              with sparse matrices.

    """
    op4 = OP4(log=log, debug=debug)
    return op4.read_op4(op4_filename, matrix_names, precision)


class OP4:
    """
    todo:: add endian checking
    todo:: test on big matrices
    todo:: finish write_op4

    """
    def __init__(self, log=None, debug=False):
        self.n = 0
        self._endian = ''
        self.debug = debug
        #assert debug == True, debug
        self.log = get_logger2(log, debug)
        self._new = False
        self.large = None

    def read_op4(self, op4_filename=None, matrix_names=None, precision='default'):
        """See ``read_op4``"""
        if precision not in ('default', 'single', 'double'):
            msg = "precision=%r and must be 'single', 'double', or 'default'" % precision
            raise ValueError(msg)

        if op4_filename is None:
            from pyNastran.utils.gui_io import load_file_dialog
            wildcard_wx = "Nastran OP4 (*.op4)|*.op4|" \
                "All files (*.*)|*.*"
            wildcard_qt = "Nastran OP4 (*.op4);;All files (*)"
            title = 'Please select a OP4 to load'
            op4_filename = load_file_dialog(title, wildcard_wx, wildcard_qt)[0]
            assert op4_filename is not None, op4_filename

        if not os.path.exists(op4_filename):
            raise IOError('cannot find op4_filename=%r' % op4_filename)

        if isinstance(matrix_names, str):
            matrix_names = [matrix_names]
        #assert isinstance(matrix_names, list), 'type(matrix_names)=%s' % type(matrix_names)

        if file_is_binary(op4_filename):
            return self.read_op4_binary(op4_filename, matrix_names, precision)
        return self.read_op4_ascii(op4_filename, matrix_names, precision)

#--------------------------------------------------------------------------
    def read_op4_ascii(self, op4_filename, matrix_names=None, precision='default'):
        """matrix_names must be a list or None, but basically the same"""
        with open(op4_filename, 'r') as op4:
            matrices = {}
            name = 'dummyName'
            while name is not None:
                (name, form, matrix) = self._read_matrix_ascii(op4, matrix_names, precision)
                if name is not None:
                    if matrix_names is None or name in matrix_names:
                        _save_matrix(matrices, name, form, matrix)
        return matrices

    def _read_matrix_ascii(self, op4, matrix_names=None, precision='default'):
        """Reads an ASCII matrix"""
        iline = 0
        line = op4.readline().rstrip()
        iline += 1
        if line == '':
            op4.close()
            return None, None, None
        ncols, nrows, form, matrix_type = line[0:32].split()
        nrows = int(nrows)

        is_big_mat, nrows = get_big_mat_nrows(nrows)
        if self.debug:
            self.log.info('is_big_matrix = %s' % is_big_mat)

        ncols = int(ncols)
        form = int(form)
        matrix_type = int(matrix_type)
        dtype = get_dtype(matrix_type, precision)

        name = line[32:40].strip()

        if self.debug:
            self.log.info('name=%s shape=(%s,%s) form=%s Type=%s' % (
                name, nrows, ncols, form, matrix_type))
        assert ncols > 0, 'ncols=%s' % ncols
        size = line[40:].strip()
        line_size = size.split(',')[1].split('E')[1].split('.')[0]  # 3E23.16 to 23
        line_size = int(line_size)

        line = op4.readline().rstrip()
        iline += 1
        (_icol, irow, _nwords) = line.split()

        is_sparse = False
        if irow == '0':
            is_sparse = True

        if matrix_type in [1, 2]:  # real
            A, iline = self._read_real_ascii(op4, iline, nrows, ncols, line_size, line,
                                             dtype, is_sparse, is_big_mat)
        elif matrix_type in [3, 4]:  # complex
            A, iline = self._read_complex_ascii(op4, iline, nrows, ncols, line_size, line,
                                                dtype, is_sparse, is_big_mat)
        else:
            raise RuntimeError('invalid matrix type.  matrix_type=%s' % matrix_type)

        if not(matrix_names is None or name in matrix_names):  # kill the matrix
            A = None

        if self.debug:
            self.log.info("form=%s name=%s A=\n%s" % (form, name, str(A)))
        return (name, form, A)

    def _read_real_sparse_ascii(self, op4, iline, nrows, ncols, line_size, line, dtype, is_big_mat):
        """Reads a sparse real ASCII matrix"""
        self.log.debug('_read_real_sparse_ascii')
        rows = []
        cols = []
        entries = []
        nloops = 0
        was_broken = False
        while 1:
            if nloops > 0 and not was_broken:
                line = op4.readline().rstrip()
                iline += 1
            was_broken = False

            icol, irow, nwords = line.split()
            icol = int(icol)

            if icol > ncols:
                break

            irow = int(irow)
            nwords = int(nwords)

            # This loop condition is overly complicated, but the first time
            # it will always execute.
            # Later if there is a sparse continuation line marker of
            # 1 (very large) integer, there will be no scientific notation value.
            # There also may be another sparse marker with 2 values.  These are not large.
            # The scientific check prevents you from getting stuck in an infinite
            # loop b/c no lines are read if there was one float value.
            # The check for 1 (or 2) integers is to prevent the check for 3 integers
            # which starts a new column.  We only want to continue a column.
            run_loop = True
            sline = line.strip().split()
            # next sparse entry
            while (len(sline) == 1 or len(sline) == 2) and 'E' not in line or run_loop:
                if is_big_mat:
                    irow, iline = self._get_irow_big_ascii(op4, iline, line, sline, irow)
                else:
                    irow, iline = self._get_irow_small_ascii(op4, iline, line, sline, irow)

                run_loop = False
                #iword = 0
                #is_done_reading_row = False
                while nwords:
                    n = 0
                    line = op4.readline().rstrip()
                    iline += 1
                    nwords_in_line = line.count('E')
                    if nwords_in_line == 0:
                        was_broken = True
                        break

                    for unused_i in range(nwords_in_line):
                        word = line[n:n + line_size]
                        rows.append(irow)
                        cols.append(icol)
                        entries.append(word)
                        if self.debug:
                            self.log.debug('  irow=%s icol=%s word=%.4g' % (
                                irow - 1, icol - 1, float(word)))
                        n += line_size
                        irow += 1
                    #iword += nwords_in_line
                    nwords -= nwords_in_line
                sline = line.strip().split()
                nloops += 1

        op4.readline()
        iline += 1

        #if rows == []:  # NULL matrix
            #raise NotImplementedError()

        rows = array(rows, dtype='int32') - 1
        cols = array(cols, dtype='int32') - 1
        A = coo_matrix((entries, (rows, cols)), shape=(nrows, ncols), dtype=dtype)
        #print("type = %s %s" % (type(A),type(A.toarray())))
        #A = A.toarray()
        return A, iline

    def _read_real_sparse_ascii_new(self, op4, iline, nrows, ncols, line_size, line,
                                    dtype, is_big_mat):
        """Reads a sparse real ASCII matrix"""
        self.log.debug('_read_real_sparse_ascii')
        rows = []
        cols = []
        entries = []
        nloops = 0
        was_broken = False
        while 1:
            if nloops > 0 and not was_broken:
                line = op4.readline().rstrip()
                iline += 1
            was_broken = False

            icol, irow, nwords = line.split()
            icol = int(icol)

            if icol > ncols:
                break

            irow = int(irow)
            nwords = int(nwords)

            # This loop condition is overly complicated, but the first time
            # it will always execute.
            # Later if there is a sparse continuation line marker of
            # 1 (very large) integer, there will be no scientific notation value.
            # There also may be another sparse marker with 2 values.  These are not large.
            # The scientific check prevents you from getting stuck in an infinite
            # loop b/c no lines are read if there was one float value.
            # The check for 1 (or 2) integers is to prevent the check for 3 integers
            # which starts a new column.  We only want to continue a column.
            run_loop = True
            sline = line.strip().split()
            # next sparse entry
            jrow1 = len(rows)
            while (len(sline) == 1 or len(sline) == 2) and 'E' not in line or run_loop:
                if is_big_mat:
                    irow, iline = self._get_irow_big_ascii(op4, iline, line, sline, irow)
                else:
                    irow, iline = self._get_irow_small_ascii(op4, iline, line, sline, irow)

                run_loop = False
                #is_done_reading_row = False
                while nwords:
                    line = op4.readline().rstrip()
                    iline += 1
                    nwords_in_line = line.count('E')
                    if nwords_in_line == 0:
                        was_broken = True
                        break

                    irows = list(range(irow, irow + nwords_in_line))
                    n = 0
                    for unused_i in range(nwords_in_line):
                        word = line[n:n + line_size]
                        entries.append(word)
                        n += line_size
                    rows.extend(irows)
                    #icols = [icol] * nwords_in_line
                    #cols.extend(icols)
                    irow += nwords_in_line
                    #assert len(rows) == len(cols), 'rows=%s\ncols=%s' % (rows, cols)
                    #assert len(rows) == len(entries)
                    nwords -= nwords_in_line
                sline = line.strip().split()
                nloops += 1
            jrow2 = len(rows)

            icols = [icol] * (jrow2 - jrow1)
            cols.append(icols)

        op4.readline()
        iline += 1

        #if rows == []:  # NULL matrix
            #raise NotImplementedError()

        cols = np.hstack(cols)
        rows = array(rows, dtype='int32') - 1
        cols = array(cols, dtype='int32') - 1
        A = coo_matrix((entries, (rows, cols)), shape=(nrows, ncols), dtype=dtype)
        #print("type = %s %s" % (type(A),type(A.toarray())))
        #A = A.toarray()
        return A, iline

    def _read_real_dense_ascii(self, op4, iline, nrows, ncols, line_size, line, dtype, is_big_mat):
        """Reads a real dense ASCII matrix"""
        self.log.debug('_read_real_dense_ascii')
        A = zeros((nrows, ncols), dtype=dtype)  # Initialize a real matrix
        nloops = 0
        was_broken = False
        while 1:
            if nloops > 0 and not was_broken:
                line = op4.readline().rstrip()
                iline += 1
            was_broken = False

            (icol, irow, nwords) = line.split()
            icol = int(icol)

            if icol > ncols:
                break

            irow = int(irow)
            nwords = int(nwords)

            # This loop condition is overly complicated, but the first time
            # it will always execute.
            # Later if there is a sparse continuation line marker of
            # 1 (very large) integer, there will be no scientific notation value.
            # There also may be another sparse marker with 2 values.  These are not large.
            # The scientific check prevents you from getting stuck in an infinite
            # loop b/c no lines are read if there was one float value.
            # The check for 1 (or 2) integers is to prevent the check for 3 integers
            # which starts a new column.  We only want to continue a column.
            run_loop = True
            sline = line.strip().split()
            # next dense entry
            while (len(sline) == 1 or len(sline) == 2) and 'E' not in line or run_loop:
                run_loop = False
                #i = 0
                #iword = 0
                #is_done_reading_row = False
                while nwords:
                    n = 0
                    line = op4.readline().rstrip()
                    iline += 1
                    nwords_in_line = line.count('E')
                    if nwords_in_line == 0:
                        was_broken = True
                        break

                    for unused_i in range(nwords_in_line):
                        word = line[n:n + line_size]
                        A[irow - 1, icol - 1] = word
                        n += line_size
                        irow += 1
                    #iword += nwords_in_line
                    nwords -= nwords_in_line
                sline = line.strip().split()
                nloops += 1
        op4.readline()
        iline += 1
        return A, iline

    def _read_real_ascii(self, op4, iline, nrows, ncols, line_size, line, dtype,
                         is_sparse, is_big_mat):
        """Reads a real ASCII matrix"""
        if is_sparse:
            if self._new:
                A, iline = self._read_real_sparse_ascii_new(op4, iline, nrows, ncols,
                                                            line_size, line, dtype, is_big_mat)
            else:
                A, iline = self._read_real_sparse_ascii(op4, iline, nrows, ncols,
                                                        line_size, line, dtype, is_big_mat)
        else:
            A, iline = self._read_real_dense_ascii(op4, iline, nrows, ncols,
                                                   line_size, line, dtype, is_big_mat)
        return A, iline

    def _read_complex_sparse_ascii(self, op4, iline, nrows, ncols, line_size,
                                   line, dtype, is_big_mat):
        """Reads a sparse complex ASCII matrix"""
        rows = []
        cols = []
        entries = []
        nloops = 0
        was_broken = False
        while 1:
            if nloops > 0 and not was_broken:
                line = op4.readline().rstrip()
                iline += 1
            was_broken = False

            (icol, irow, nwords) = line.split()
            icol = int(icol)

            if icol > ncols:
                break

            irow = int(irow)
            nwords = int(nwords)

            run_loop = True
            sline = line.strip().split()
            # next sparse entry
            while (len(sline) == 1 or len(sline) == 2) and 'E' not in line or run_loop:
                if is_big_mat:
                    irow, iline = self._get_irow_big_ascii(op4, iline, line, sline, irow)
                else:
                    irow, iline = self._get_irow_small_ascii(op4, iline, line, sline, irow)
                run_loop = False

                #i = 0
                is_real = True
                #is_done_reading_row = False
                while nwords:
                    n = 0
                    line = op4.readline().rstrip()
                    iline += 1
                    nwords_in_line = line.count('E')
                    if nwords_in_line == 0:
                        was_broken = True
                        break

                    for unused_i in range(nwords_in_line):
                        value = float(line[n:n + line_size])

                        if is_real:
                            real_value = value
                            is_real = False
                        else:
                            rows.append(irow)
                            cols.append(icol)
                            entries.append(complex(real_value, value))
                            irow += 1
                            is_real = True
                        n += line_size
                    nwords -= nwords_in_line
                sline = line.strip().split()
                nloops += 1

        rows = array(rows, dtype='int32') - 1
        cols = array(cols, dtype='int32') - 1
        A = coo_matrix((entries, (rows, cols)), shape=(nrows, ncols), dtype=dtype)
        op4.readline()
        iline += 1
        return A, iline

    def _read_complex_ascii(self, op4, iline, nrows, ncols, line_size, line,
                            dtype, is_sparse, is_big_mat):
        """Reads a complex ASCII matrix"""
        if is_sparse:
            A, iline = self._read_complex_sparse_ascii(op4, iline, nrows, ncols,
                                                       line_size, line, dtype, is_big_mat)
        else:
            A, iline = self._read_complex_dense_ascii(op4, iline, nrows, ncols,
                                                      line_size, line, dtype, is_big_mat)
        return A, iline

    def _read_complex_dense_ascii(self, op4, iline, nrows, ncols, line_size, line, dtype,
                                  is_big_mat):
        """Reads a dense complex ASCII matrix"""
        A = zeros((nrows, ncols), dtype=dtype)  # Initialize a complex matrix

        nloops = 0
        was_broken = False
        while 1:
            if nloops > 0 and not was_broken:
                line = op4.readline().rstrip()
                iline += 1
            was_broken = False

            (icol, irow, nwords) = line.split()
            icol = int(icol)

            if icol > ncols:
                break

            irow = int(irow)
            nwords = int(nwords)

            run_loop = True
            sline = line.strip().split()

            # next sparse entry
            while (len(sline) == 1 or len(sline) == 2) and 'E' not in line or run_loop:
                run_loop = False

                iword = 0
                #unused_is_done_reading_row = False
                while nwords:
                    n = 0
                    line = op4.readline().rstrip()
                    iline += 1
                    nwords_in_line = line.count('E')
                    if nwords_in_line == 0:
                        was_broken = True
                        break

                    for unused_i in range(nwords_in_line):
                        value = float(line[n:n + line_size])

                        if iword % 2 == 0:
                            #A[irow - 1, icol - 1].real = value
                            real_value = value
                        else:
                            A[irow - 1, icol - 1] = complex(real_value, value)
                            irow += 1
                        iword += 1
                        n += line_size
                    nwords -= nwords_in_line
                sline = line.strip().split()
                iline += 1
                nloops += 1

        op4.readline()
        iline += 1
        return A, iline

    def _get_irow_small_ascii(self, op4, iline, line, sline, irow):
        sline = line.strip().split()
        if len(sline) == 1:
            IS = int(line)
        else:
            line = op4.readline().strip()
            try:
                IS = int(line)
            except ValueError:
                msg = 'Line %i: Failed getting IROW from %r' % (iline, line)
                raise ValueError(msg)
            iline += 1
        L = IS // 65536 - 1
        irow = IS - 65536 * (L + 1)
        if self.debug:
            self.log.info('small_mat-next row')
            self.log.info('  IS=%s L=%s irow=%s' % (IS, L, irow))
        return irow, iline

    def _get_irow_small_binary(self, op4, data):
        """
        Returns
        -------
        irow : int
           the row id
        L : int
            the row length

        """
        if len(data) == 0:
            data = op4.read(4)
            self.n += 4

        IS, = unpack(self._endian + 'i', data)
        L = IS // 65536 - 1
        irow = IS - 65536 * (L + 1)
        if self.debug:
            self.log.info('small_mat-next row')
            self.log.info("  IS=%s L=%s irow=%s" % (IS, L, irow))
            assert IS > 0, IS
            assert L > 0, L
        return irow, L

    def _get_irow_big_ascii(self, op4, iline, line, sline, irow):
        sline = line.strip().split()
        if len(sline) == 2:
            pass
        else:
            sline = op4.readline().strip().split()
            iline += 1
        assert len(sline) == 2, 'sline=%s len(sline)=%s' % (sline, len(sline))
        (idummy, irow) = sline
        irow = int(irow)
        if self.debug:
            self.log.debug("idummy=%s irow=%s" % (idummy, irow))
        return irow, iline

    def _get_irow_big_binary(self, op4, data):
        """
        Returns
        -------
        irow : int
           the row id
        L : int
            ???

        """
        if len(data) == 0:
            data = op4.read(8)
            self.n += 8
        idummy, irow = unpack(self._endian + '2i', data)
        if self.debug:
            self.log.debug("idummy=%s irow=%s" % (idummy, irow))
            assert irow < 100, irow
        return (irow, idummy - 1)

#--------------------------------------------------------------------------
    def read_op4_binary(self, op4_filename, matrix_names=None, precision='default'):
        """matrix_names must be a list or None, but basically the same"""
        with open(op4_filename, mode='rb') as op4:
            self.n = 0
            self._endian = self._determine_endian(op4)

            matrices = {}
            name = 'dummyName'
            while name is not None:
                #print "i =", i
                # checks for the end of the file
                assert self.n == op4.tell(), 'n=%s tell=%s' % (self.n, op4.tell())
                n = self.n
                data1 = op4.read(1)
                op4.seek(n)
                if len(data1) == 0:
                    break
                #self.show(f, 60)

                (name, form, matrix) = self._read_matrix_binary(op4, precision, matrix_names)
                #print(print_matrix(matrix))
                if name is not None:
                    if matrix_names is None or name in matrix_names:  # save the matrix
                        name = name.decode('ascii')
                        _save_matrix(matrices, name, form, matrix)

                #print("not op4.closed = ",not op4.closed,form,name)
                # if not op4.closed or form is not None:
                #     data = op4.read(4)
                #     self.n += 4
                #     if len(data) == 0:
                #         break
                #     (record_length,) = unpack(self._endian + 'i', data)
                ##     print("record_length = %s" % record_length)
                #     if record_length == 24:
                #         self.n -= 4
                #         op4.seek(self.n)
                #     else:
                #         data = op4.read(4)
                #         if len(data) == 0:
                #             break
                #         (record_length2,) = unpack(self._endian + 'i', data)
                #         assert record_length2 == 24
                #         op4.seek(self.n)
                #
        return matrices

    def read_start_marker(self, op4):
        if self.debug:
            self.log.info('--------------------------------------')
        #self.show(op4, 60)
        data = op4.read(4)
        self.n += 4
        record_length, = unpack(self._endian + 'i', data)
        #print('record_length =', record_length)

        record_length = 16
        data = op4.read(record_length)
        self.n += record_length
        assert self.n == op4.tell(), 'n=%s tell=%s' % (self.n, op4.tell())

        if record_length == 16:
            a, icol, irow, nwords = unpack(self._endian + '4i', data)
            if self.debug:
                self.log.info("a=%s icol=%s irow=%s nwords=%s" % (a, icol, irow, nwords))
        else:
            raise NotImplementedError('record_length=%s' % record_length)
        return (a, icol, irow, nwords)

    def _read_matrix_binary(self, op4, precision, matrix_names):
        """Reads a binary matrix"""
        #self.show(f, 60)
        if self.debug:
            self.log.info("*************************")
        data = op4.read(4)
        self.n += 4
        (record_length,) = unpack(self._endian + 'i', data)
        assert self.n == op4.tell(), 'n=%s tell=%s' % (self.n, op4.tell())
        if self.debug:
            self.log.info("record_length = %s" % record_length)

        if record_length == 24:
            data = op4.read(record_length)
            self.n += record_length

            (ncols, nrows, form, Type, name) = unpack(
                self._endian + '4i8s', data)
            if self.debug:
                self.log.info("nrows=%s ncols=%s form=%s Type=%s name=%r" % (
                    nrows, ncols, form, Type, name))
        elif record_length == 48:
            data = op4.read(record_length)
            self.n += record_length

            #self._show_data(data, types='ifdlqILQs', endian=None)
            (ncols, nrows, form, Type, name) = unpack(
                self._endian + '4Q16s', data)
            if self.debug:
                self.log.info("nrows=%s ncols=%s form=%s Type=%s name=%r" % (
                    nrows, ncols, form, Type, name))
        else:
            #msg = record_length #+ self.print_block(data)
            msg = 'record_length=%s filename=%r' % (record_length, op4.name)
            raise NotImplementedError(msg)

        name = name.strip()
        if self.debug:
            if Type == 1:
                self.log.info("Type = Real, Single Precision")
            elif Type == 2:
                self.log.info("Type = Real, Double Precision")
            elif Type == 3:
                self.log.info("Type = Complex, Single Precision")
            elif Type == 4:
                self.log.info("Type = Complex, Double Precision")

        is_big_mat, nrows = get_big_mat_nrows(nrows)

        if self.debug:
            self.log.info('is_big_matrix = %s' % is_big_mat)

        # jump forward to get irow (needed for check on is_sparse),
        # then jump back
        nsave = self.n
        irow = self.read_start_marker(op4)[2]
        op4.seek(nsave)
        self.n = nsave

        #(nwords_per_value, nbytes_per_value, data_format, dtype) = self._get_matrix_info(Type)
        data_format, dtype = self._get_matrix_info(Type)[2:]

        is_sparse = False
        if irow == 0:
            is_sparse = True

        assert self.n == op4.tell(), 'n=%s tell=%s' % (self.n, op4.tell())
        if Type in [1, 2]:  # real
            A = self._read_real_binary(op4, nrows, ncols, Type, is_sparse, is_big_mat)
        elif Type in [3, 4]:  # complex
            A = self._read_complex_binary(op4, nrows, ncols, Type, is_sparse, is_big_mat)
        else:
            self.log.error('is_sparse=%s data_format=%s dtype=%s' % (is_sparse, data_format, dtype))
            raise TypeError("Type=%s" % Type)

        #try:
            #print_matrix(A.toarray())
        #except:
            #pass

        if data_format in ['d', 'dd']:
            op4.read(8)
            self.n += 8
        elif data_format in ['f', 'ff']:
            op4.read(4)
            self.n += 4
        else:
            raise NotImplementedError(data_format)
        #f.read(record_length); self.n+=record_length
        #self.show(f, 10)
        #f.read(4); self.n+=4

        assert self.n == op4.tell(), 'n=%s op4.tell=%s' % (self.n, op4.tell())
        return (name, form, A)

    def _get_matrix_info(self, matrix_type, debug=True):
        if matrix_type == 1:
            nwords_per_value = 1  # number words per value
            nbytes_per_value = 4
            data_format = 'f'
        elif matrix_type == 2:
            nwords_per_value = 2
            nbytes_per_value = 8
            data_format = 'd'
        elif matrix_type == 3:
            nwords_per_value = 2
            nbytes_per_value = 8
            data_format = 'f'
        elif matrix_type == 4:
            nwords_per_value = 4
            nbytes_per_value = 16
            data_format = 'd'
        else:
            raise RuntimeError("matrix_type=%s" % matrix_type)
        dtype = get_dtype(matrix_type)
        if self.debug and debug:
            self.log.info('matrix_type = %s' % matrix_type)
            self.log.info('  nwords_per_value = %s' % nwords_per_value)
            self.log.info('  nbytes_per_value = %s' % nbytes_per_value)
            self.log.info('  dtype = %s ' % dtype)
        return (nwords_per_value, nbytes_per_value, data_format, dtype)

    def _read_real_dense_binary(self, op4, nrows, ncols, matrix_type, is_big_mat):
        if self.debug:
            self.log.info('_read_real_dense_binary')
        out = self._get_matrix_info(matrix_type, debug=False)
        (nwords_per_value, _nbytes_per_value, data_format, dtype) = out
        A = zeros((nrows, ncols), dtype=dtype)

        icol = -1  # dummy value so the loop starts
        while icol < ncols + 1:  # if isDense
            assert self.n == op4.tell(), 'n=%s tell=%s' % (self.n, op4.tell())
            (icol, irow, nwords) = self.get_markers_dense(op4)
            L = nwords
            if icol == ncols + 1:
                break
            if L == -1:
                break

            record_length = 4 * nwords
            data = op4.read(record_length)
            self.n += record_length
            nvalues = L // nwords_per_value
            str_values = self._endian + '%i%s' % (nvalues, data_format)
            A[irow-1:irow-1+nvalues, icol-1] = unpack(str_values, data)
            if self.debug:
                self.log.info('A[%s:%s, %s] = %s' % (
                    irow - 1,
                    irow - 1 + nvalues,
                    icol-1,
                    A[irow-1:irow-1+nvalues, icol-1]))
        #assert self.n == op4.tell(), 'n=%s tell=%s' % (self.n, op4.tell())
        op4.read(4)
        self.n += 4
        return A


    def _read_real_binary(self, op4, nrows, ncols, matrix_type, is_sparse, is_big_mat):
        if is_sparse:
            A = self._read_real_sparse_binary(op4, nrows, ncols, matrix_type, is_big_mat)
        else:
            A = self._read_real_dense_binary(op4, nrows, ncols, matrix_type, is_big_mat)
        return A

    def _read_real_sparse_binary(self, op4, nrows, ncols, matrix_type, is_big_mat):
        if self.debug:
            self.log.info('_read_real_sparse_binary')
        out = self._get_matrix_info(matrix_type, debug=False)
        (nwords_per_value, nbytes_per_value, data_format, dtype) = out
        rows = []
        cols = []
        entries = []

        data = ''
        icol = -1  # dummy value so the loop starts
        while icol < ncols + 1:  # if isDense
            #if record_length==0:
            assert self.n == op4.tell(), 'n=%s tell=%s' % (self.n, op4.tell())
            (icol, irow, nwords) = self.get_markers_sparse(op4, is_big_mat)
            L = nwords

            if icol == ncols + 1:
                if self.debug:
                    self.log.info('breaking on icol=%s ncol+1=%s' % (icol, ncols + 1))
                break

            if is_big_mat:
                irow, L = self._get_irow_big_binary(op4, data[:8])
                data = data[8:]
            else:
                irow, L = self._get_irow_small_binary(op4, data[:4])
                data = data[4:]

            if L == -1:
                if self.debug:
                    self.log.info('breaking on L=-1')
                break

            if self.debug:
                self.log.info("  next icol")
                self.log.info("    n=%s icol=%s irow=%s nwords=%s" % (
                    self.n, icol, irow, nwords))
                self._show(op4, 100, types='qd')
                self.log.info('**************************************************')

            #if nwords == 0 and is_big_mat:
                #self.n -= 4
                #op4.seek(self.n)
                #break

            record_length = 4 * nwords
            data = op4.read(record_length)
            self.n += record_length
            if self.debug:
                self.log.info("  data_format=%s record_length=%s n_next=%s" % (
                    data_format, record_length, self.n))
            #if icol == ncols + 1:
                #break

            i = 0
            while len(data) > 0:
                if i > 0:
                    if is_big_mat:
                        (irow, L) = self._get_irow_big_binary(op4, data[0:8])
                        data = data[8:]
                    else:
                        (irow, L) = self._get_irow_small_binary(op4, data[0:4])
                        data = data[4:]
                    assert irow > 0
                nvalues = L // nwords_per_value
                str_values = self._endian + '%i%s' % (nvalues, data_format)

                if self.debug:
                    self.log.info('irow=%s L=%s nwords_per_value=%s nvalues=%s '
                                  'nbytes_per_value=%s' % (irow, L, nwords_per_value,
                                                           nvalues, nbytes_per_value))
                    self.log.info('str_values = %r' % str_values)

                value_list = unpack(str_values, data[0:nvalues * nbytes_per_value])
                assert self.n == op4.tell(), 'n=%s tell=%s' % (self.n, op4.tell())

                #irow -= 1
                #icol -= 1
                if self.debug:
                    self.log.info('rows = %s' % list(i+irow-1 for i in range(nvalues)))
                    self.log.info('cols = %s ' % ([icol-1] * nvalues))
                    self.log.info('value_list = %s' % str(value_list))

                rows.extend([i+irow-1 for i in range(nvalues)])
                irow += nvalues
                cols.extend([icol-1] * nvalues)
                entries.extend(value_list)

                record_length -= nvalues * nbytes_per_value
                data = data[nvalues * nbytes_per_value:]
                if self.debug:
                    self.log.info("  record_length=%s nbytes_per_value=%s len(data)=%s" % (
                        record_length, nbytes_per_value, len(data)))
                    ##print(A)
                    #print("********")  # ,data
                    #print(self.print_block(data))
                i += 1
            #print "-------------------------------"

        #if rows == []:  # NULL matrix
            #raise NotImplementedError()

        A = coo_matrix((entries, (rows, cols)), shape=(nrows, ncols),
                       dtype=dtype)
        op4.read(4)
        self.n += 4

        return A

    def _show(self, op4, n, types='ifs', endian=None):
        """Shows binary data"""
        assert self.n == op4.tell()
        nints = n // 4
        data = op4.read(4 * nints)
        strings, ints, floats = self._show_data(data, types=types, endian=endian)
        op4.seek(self.n)
        return strings, ints, floats

    def _show_data(self, data, types='ifs', endian=None):
        """
        Shows a data block as various types

        Parameters
        ----------
        data : bytes
            the binary string bytes
        types : str; default='ifs'
            i - int
            f - float
            s - string
            d - double (float; 8 bytes)

            l - long (int; 4 bytes)
            q - long long (int; int; 8 bytes)
            I - unsigned int (int; 4 bytes)
            L - unsigned long (int; 4 bytes)
            Q - unsigned long long (int; 8 bytes)
        endian : str; default=None -> auto determined somewhere else in the code
            the big/little endian {>, <}

        .. warning:: 's' is apparently not Python 3 friendly

        """
        return self._write_data(sys.stdout, data, types=types, endian=endian)

    def _write_data(self, f, data, types='ifs', endian=None):
        """
        Useful function for seeing what's going on locally when debugging.

        Parameters
        ----------
        data : bytes
            the binary string bytes
        types : str; default='ifs'
            i - int
            f - float
            s - string
            d - double (float; 8 bytes)

            l - long (int; 4 bytes)
            q - long long (int; int; 8 bytes)
            I - unsigned int (int; 4 bytes)
            L - unsigned long (int; 4 bytes)
            Q - unsigned long long (int; 8 bytes)
        endian : str; default=None -> auto determined somewhere else in the code
            the big/little endian {>, <}
        types = 'ifdlqILQ'

        """
        n = len(data)
        nints = n // 4
        ndoubles = n // 8
        strings = None
        ints = None
        floats = None
        longs = None

        if endian is None:
            endian = self._endian
            assert endian is not None, endian

        if 's' in types:
            strings = unpack(b'%s%is' % (endian, n), data)
            f.write("strings(s) = %s\n" % str(strings))
        if 'i' in types:
            ints = unpack(b'%s%ii' % (endian, nints), data)
            f.write("ints(i)    = %s\n" % str(ints))
        if 'f' in types:
            floats = unpack(b'%s%if' % (endian, nints), data)
            f.write("floats(f)  = %s\n" % str(floats))
        if 'd' in types:
            doubles = unpack(b'%s%id' % (endian, ndoubles), data[:ndoubles*8])
            f.write("doubles(d)  = %s\n" % str(doubles))

        if 'l' in types:
            longs = unpack(b'%s%il' % (endian, nints), data)
            f.write("long(l)  = %s\n" % str(longs))
        if 'I' in types:
            ints2 = unpack(b'%s%iI' % (endian, nints), data)
            f.write("unsigned int(I) = %s\n" % str(ints2))
        if 'L' in types:
            longs2 = unpack(b'%s%iL' % (endian, nints), data)
            f.write("unsigned long(L) = %s\n" % str(longs2))
        if 'q' in types:
            longs = unpack(b'%s%iq' % (endian, ndoubles), data[:ndoubles*8])
            f.write("long long(q) = %s\n" % str(longs))
        if 'Q' in types:
            longs = unpack(b'%s%iq' % (endian, ndoubles), data[:ndoubles*8])
            f.write("unsigned long long(Q) = %s\n" % str(longs))
        return strings, ints, floats

    def _show_ndata(self, f, n, types='ifs'):
        return self._write_ndata(sys.stdout, f, n, types=types)

    def _write_ndata(self, fout, f, n, types='ifs'):
        """Useful function for seeing what's going on locally when debugging."""
        nold = self.n
        data = f.read(n)
        self.n = nold
        f.seek(self.n)
        return self._write_data(fout, data, types=types)

    def _read_complex_dense_binary(self, op4, nrows, ncols, matrix_type, is_big_mat):
        """reads a dense complex binary matrix"""
        if self.debug:
            self.log.info('_read_complex_dense_binary')
        out = self._get_matrix_info(matrix_type, debug=False)
        (nwords_per_value, nbytes_per_value, data_format, dtype) = out

        A = zeros((nrows, ncols), dtype=dtype)
        record_length = 0
        icol = -1  # dummy value so the loop starts
        while icol < ncols + 1:  # if isDense
            #if record_length==0:
            assert self.n == op4.tell(), 'n=%s tell=%s' % (self.n, op4.tell())
            (icol, irow, nwords) = self.get_markers_dense(op4)
            if self.debug:
                self.log.info("N=%s icol=%s irow=%s nwords=%s" % (
                    self.n, icol, irow, nwords))
                self.log.info("-----------")

            L = nwords
            if icol == ncols + 1:
                if self.debug:
                    self.log.info('breaking...icol=%s ncols+1=%s' % (
                        icol, ncols + 1))
                break

            if L == -1:
                if self.debug:
                    self.log.info('breaking...nwords (L) = %s' % nwords)
                break

            if self.debug:
                self.log.info("  n=%s icol=%s irow=%s nwords=%s" % (
                    self.n, icol, irow, nwords))

            #if nwords == 0 and is_big_mat:
                #self.n -=4
                #f.seek(self.n)
                #break

            record_length = 4 * nwords
            data = op4.read(record_length)
            self.n += record_length
            if self.debug:
                self.log.info("data_format=%s record_length=%s n_next=%s" % (
                    data_format, record_length, self.n))
            if icol == ncols + 1:
                continue

            nvalues = nwords // nwords_per_value
            #nread = nwords // 4
            while record_length >= nbytes_per_value:
                if self.debug:
                    self.log.info("inner while...")
                    self.log.info("nwords  = %s" % nwords)
                    self.log.info("nvalues = %s" % nvalues)
                    self.log.info("nwords_per_value = %s" % nwords_per_value)

                #if nvalues == 0:
                    #assert icol == ncols + 1
                    #break

                # we have more 2x values for complex numbers
                str_values = self._endian + '%i%s' % (nvalues * 2, data_format)
                if self.debug:
                    self.log.info("str_values = %s" % str_values)
                    self.log.info("nvalues*nbytes_per_value=%s len(data)=%s" % (
                        nvalues * nbytes_per_value, len(data)))
                value_list = unpack(str_values, data[0:nvalues * nbytes_per_value])
                assert self.n == op4.tell(), 'n=%s tell=%s' % (self.n, op4.tell())
                #self.show(op4, 4)
                #print self.print_block(data)
                if self.debug:
                    self.log.info("value_list = %s" % str(value_list))

                #irow -= 1
                #icol -= 1
                irow -= 1
                icol -= 1
                for i, value in enumerate(value_list):
                    if i % 2 == 0:
                        real_value = value
                    else:
                        ai = complex(real_value, value)
                        if self.debug:
                            self.log.info("A[%s,%s] = %s" % (irow, icol, ai))
                        A[irow, icol] = ai
                        irow += 1

                record_length -= nvalues * nbytes_per_value
                data = data[nvalues * nbytes_per_value:]
                if self.debug:
                    self.log.info("record_length=%s nbytes_per_value=%s" % (
                        record_length, nbytes_per_value))
                    self.log.info(print_matrix(A))
                    self.log.info("******** %s" % data)

        op4.read(4)
        self.n += 4
        return A

    def _read_complex_binary(self, op4, nrows, ncols, matrix_type, is_sparse, is_big_mat):
        """Reads a complex binary matrix"""
        if is_sparse:
            A = self._read_complex_sparse_binary(op4, nrows, ncols, matrix_type, is_big_mat)
        else:
            A = self._read_complex_dense_binary(op4, nrows, ncols, matrix_type, is_big_mat)
        return A

    def _read_complex_sparse_binary(self, op4, nrows, ncols, matrix_type, is_big_mat):
        """Reads a sparse complex binary matrix"""
        if self.debug:
            self.log.info('_read_complex_sparse_binary')
        out = self._get_matrix_info(matrix_type, debug=False)
        (nwords_per_value, nbytes_per_value, data_format, dtype) = out
        rows = []
        cols = []
        entries = []
        record_length = 0
        data = ''
        icol = -1  # dummy value so the loop starts
        while icol < ncols + 1:  # if isDense
            #if record_length == 0:
            assert self.n == op4.tell(), 'n=%s tell=%s' % (self.n, op4.tell())
            (icol, irow, nwords) = self.get_markers_sparse(op4, is_big_mat)
            if self.debug:
                self.log.info("n=%s icol=%s irow=%s nwords=%s" % (
                    self.n, icol, irow, nwords))
                self.log.info("-----------")

            L = nwords
            if icol == ncols + 1:
                break

            if is_big_mat:
                (irow, L) = self._get_irow_big_binary(op4, data[:8])
                data = data[8:]
            else:
                (irow, L) = self._get_irow_small_binary(op4, data[:4])
                data = data[4:]

            if L == -1:
                if self.debug:
                    self.log.info('breaking on L=-1')
                break

            if self.debug:
                self.log.info("n=%s icol=%s irow=%s nwords=%s" % (
                    self.n, icol, irow, nwords))
                self._show(op4, 100, types='qf')
                self.log.info('\n\n')

            #if nwords == 0 and is_big_mat:
                #self.n -= 4
                #op4.seek(self.n)
                #break

            record_length = 4 * nwords
            data = op4.read(record_length)
            self.n += record_length
            if self.debug:
                self.log.info("data_format=%s record_length=%s n_next=%s" % (
                    data_format, record_length, self.n))
            if icol == ncols + 1:
                continue

            nvalues = nwords // nwords_per_value
            while record_length >= nbytes_per_value:
                if self.debug:
                    self.log.info("inner while...")
                    self.log.info("nwords  = %s" % nwords)
                    self.log.info("nvalues = %s" % nvalues)
                    self.log.info("nwords_per_value = %s" % nwords_per_value)

                #if nvalues == 0:
                    #assert icol == ncols + 1
                    #break

                # we have 2x values for complex
                str_values = self._endian + '%i%s' % (nvalues * 2, data_format)
                if self.debug:
                    self.log.info("  str_values = %s" % str_values)
                    self.log.info("  nvalues*nbytes_per_value=%s len(data)=%s" % (
                        nvalues * nbytes_per_value, len(data)))
                value_list = unpack(str_values, data[0:nvalues * nbytes_per_value])
                assert self.n == op4.tell(), 'n=%s tell=%s' % (self.n, op4.tell())
                #self.show(op4, 4)
                #print(self.print_block(data))
                if self.debug:
                    self.log.info("  value_list = %s" % str(value_list))

                #irow -= 1
                #icol -= 1
                irow -= 1
                icol -= 1

                cols += [icol] * nvalues
                rows += [i + irow for i in range(nvalues)]
                for i, value in enumerate(value_list):
                    if i % 2 == 0:
                        real_value = value
                    else:
                        if self.debug:
                            self.log.info("  A[%s,%s] = %s" % (
                                irow, icol, complex(real_value, value)))
                        #A[irow, icol] = complex(real_value, value)
                        entries.append(complex(real_value, value))
                        irow += 1

                record_length -= nvalues * nbytes_per_value
                data = data[nvalues * nbytes_per_value:]
                #print("record_length=%s nbytes_per_value=%s" % (record_length, nbytes_per_value))
                #print(print_matrix(A))
                #print("********", data)

        A = coo_matrix((entries, (rows, cols)), shape=(nrows, ncols), dtype=dtype)
        op4.read(4)
        self.n += 4
        return A

    def get_markers_sparse(self, op4, is_big_mat):
        if is_big_mat:
            (unused_a, icol, irow, nwords) = self.read_start_marker(op4)
            #irow = self._get_irow_big(op4)
            nwords -= 2
            #if nwords > 1:
               #nwords -= 2
            #else:
               #print("nwords0 = %s" % nwords)
               #nwords = 0
        else:
            unused_a, icol, irow, nwords = self.read_start_marker(op4)
            #if irow != 0:
                #assert nwords == 1, 'nwords=%s' % nwords

            #irow = self._get_irow_small(f)
            nwords -= 1
        return icol, irow, nwords

    def get_markers_dense(self, op4):
        a, icol, irow, nwords = self.read_start_marker(op4)
        if self.debug:
            self.log.info("n=%s a=%s icol=%s irow=%s nwords=%s"% (
                self.n, a, icol, irow, nwords))
        return icol, irow, nwords

    def write_op4(self, op4_filename, matrices, name_order=None,
                  precision='default', is_binary=True):
        """
        Writes the OP4

        Parameters
        ----------
        op4_filename : str/file
            The filename to write
            String -> opens a file (closed at the end)
            file   -> no file is opened and it's not closed
        matrices : Dict[str] = (form, np.ndarray)
            the matrices to write

        name_order: str / List[str]; default=None -> sorted based on name
            List of the names of the matrices that should be
            written or string
        is_binary : bool; default=True
            Should a binary file be written
        precision : str; default='default'
            Overwrite the default precision ('single', 'double', 'default')
            Applies to all matrices

        Examples
        --------
        # simple
        >>> write_op4(op4_filename, name_order=['A', 'B', 'C'],
                      precision='default', is_binary=True)

        # another method
        >>> matrices = {
            'A' : (formA, matrixA),
            'B' : (formB, matrixB),
        }

        .. todo::  This method is not even close to being done

        """
        if precision not in ('single', 'double', 'default'):
            msg = "precision=%r and must be 'single', 'double', or 'default'" % precision
            raise ValueError(msg)
        if is_binary not in (True, False):
            raise ValueError('is_binary=%r and must be True or False' % is_binary)
        #if nR == nC: op4_form = 1   # square
        #else:        op4_form = 2   # rectangular

        if isinstance(op4_filename, str):
            with open(op4_filename, 'w') as op4:
                self._write_op4_file(op4, name_order, is_binary, precision, matrices)
        else:
            op4 = op4_filename
            self._write_op4_file(op4, name_order, is_binary, precision, matrices)

    def _write_op4_file(self, op4, name_order, is_binary, precision, matrices):
        """Helper method for OP4 writing"""
        if name_order is None:
            name_order = sorted(matrices.keys())
        elif isinstance(name_order, str):
            name_order = [name_order]
        elif isinstance(name_order, bytes):
            name_order = [name_order]

        is_big_mat = False  ## .. todo:: hardcoded
        for name in name_order:
            try:
                (form, matrix) = matrices[name]
            except KeyError:
                raise KeyError('key=%r is an invalid matrix; keys=%s' % (
                    str(name), matrices.keys()))
            if not form in (1, 2, 3, 6, 8, 9):
                raise ValueError('form=%r and must be in [1, 2, 3, 6, 8, 9]' % form)

            if isinstance(matrix, coo_matrix):
                #write_DMIG(f, name, matrix, form, precision='default')
                if is_binary:
                    raise NotImplementedError('sparse binary op4 writing not implemented')
                else:
                    _write_sparse_matrix_ascii(
                        op4, name, matrix, form=form,
                        precision=precision, is_big_mat=is_big_mat)
            elif isinstance(matrix, ndarray):
                if is_binary:
                    self._write_dense_matrix_binary(
                        op4, name, matrix, form=form, precision=precision)
                else:
                    self._write_dense_matrix_ascii(
                        op4, name, matrix, form=form, precision=precision)
            else:
                msg = ('Matrix type=%r is not supported.  '
                       'types=[coo_matrix, ndarray]' % type(matrix))
                raise NotImplementedError(msg)


    def __backup(self, name: str, matrix: np.ndarray, form: int=2, precision: str='default'):
        """
        Put this documentation somewhere else...

        Parameters
        ----------
        name : str
            the name of the matrix
        matrix : ndarray
            a two-dimensional NUMPY.NDARRAY
        form : int (default=2)
            Form is defined as one of the following:
        precision : str; default=True
            {'default', 'single', 'double'}

        +======+================+
        | Form |   Definition   |
        +======+================+
        |  1   | Square         |
        +------+----------------+
        |  2   | Rectangular    |
        +------+----------------+
        |  3   | Diagonal       |
        +------+----------------+
        |  6   | Symmetric      |
        +------+----------------+
        |  8   | Id entity      |
        +------+----------------+
        |  9   | Pseudoidentity |
        +------+----------------+

        Not Supported by all OP4s (this is not a restriction of the OP4
        reader/writer)

        +======+================================+
        | Form |         Definition             |
        +======+================================+
        |  4   | Lower triangular factor        |
        +------+--------------------------------+
        |  5   | Upper triangular factor        |
        +------+--------------------------------+
        |  10  | Cholesky factor                |
        +------+--------------------------------+
        |  11  | Trapezoidal factor             |
        +------+--------------------------------+
        |  13  | Sparse lower triangular factor |
        +------+--------------------------------+
        |  15  | Sparse upper triangular factor |
        +------+--------------------------------+

        .. note:: form defaults to 2, but 1 can be easily determined.
                  Any others must be specified.

        """
        assert isinstance(name, str), name
        assert isinstance(form, int), form

    def _write_dense_matrix_binary(self, op4, name, matrix, form=2,
                                   precision='default', encoding='utf-8'):
        """
        24 bytes is the record length

        +-------------+--------------+
        | Word Number | Variable     |
        +-------------+--------------+
        |      1      |  ncols       |
        +-------------+--------------+
        |      2      |  nrows       |
        +-------------+--------------+
        |      3      |  form        |
        +-------------+--------------+
        |      4      |  matrix_type |
        +-------------+--------------+
        |    5, 6     |  name        |
        +-------------+--------------+
        6 words * 4 bytes/word = 24 bytes

        .. todo:: support precision

        """
        A = matrix
        matrix_type, nwords_per_value = _get_type_nwv(A[0, 0], precision)

        (nrows, ncols) = A.shape
        #if nrows==ncols and form==2:
            #form = 1
        name2 = '%-8s' % name
        name_bytes = name2.encode('ascii')
        assert len(name2) == 8, 'name=%r is too long; 8 characters max' % name
        s = Struct(self._endian + '5i8s')
        msg = s.pack(24, ncols, nrows, form, matrix_type, name_bytes)
        op4.write(msg)

        for icol in range(ncols):
            (istart, iend) = _get_start_end_row(A[:, icol], nrows)

            # write the column
            if istart is not None and iend is not None:
                iend += 1
                msg = pack(self._endian + '4i', 24, icol +
                           1, istart + 1, (iend - istart) * nwords_per_value)

                if matrix_type in [1, 2]: # real
                    if matrix_type == 1: # real, single
                        fmt = '%if' % (iend - istart)
                    else:         # real, double
                        fmt = '%id' % (iend - istart)
                    op4.write(pack(self._endian + fmt, *A[istart:iend+1, icol]))

                else:  # complex
                    if matrix_type == 3: # complex, single
                        fmt = '2f'
                    else:         # complex, double
                        fmt = '2d'
                    for irow in range(istart, iend):
                        msg += pack(self._endian + fmt, A[irow, icol].real,
                                    A[irow, icol].imag)
            op4.write(msg)
        if matrix_type in [1, 3]: # single precision
            # .. todo:: is this right???
            msg = pack(self._endian + '4if', 24, ncols + 1, 1, 1, 1.0)
        else: # double precision
            msg = pack(self._endian + '4id', 24, ncols + 1, 1, 1, 1.0)
        op4.write(msg)

    def _write_dense_matrix_ascii(self, op4, name, A, form=2, precision='default'):
        """Writes a dense ASCII matrix"""
        if self.debug:
            self.log.info('_write_dense_matrix_ascii')
        matrix_type, nwords_per_value = _get_type_nwv(A[0, 0], precision)

        (nrows, ncols) = A.shape
        msg = u'%8i%8i%8i%8i%-8s1P,3E23.16\n' % (ncols, nrows, form, matrix_type, name)
        op4.write(msg)

        if matrix_type in [1, 2]: # real
            for icol in range(ncols):
                value_str = ''
                (istart, iend) = _get_start_end_row(A[:, icol], nrows)

                # write the column
                if istart is not None and iend is not None:  # not a null column
                    iend += 1
                    msg = '%8i%8i%8i\n' % (icol + 1, istart + 1,
                                           (iend - istart) * nwords_per_value)
                    op4.write(msg)
                    for i, irow in enumerate(range(istart, iend)):
                        value_str += '%23.16E' % A[irow, icol]
                        if (i + 1) % 3 == 0:
                            op4.write(value_str + '\n')
                            value_str = ''
                if value_str:
                    op4.write(value_str + '\n')
        else: # complex
            for icol in range(ncols):
                value_str = ''
                (istart, iend) = _get_start_end_row(A[:, icol], nrows)

                # write the column
                if istart is not None and iend is not None:  # not a null column
                    iend += 1
                    msg = '%8i%8i%8i\n' % (icol + 1, istart + 1,
                                           (iend - istart) * nwords_per_value)
                    op4.write(msg)
                    i = 0
                    for irow in range(istart, iend):
                        value_str += '%23.16E' % A[irow, icol].real
                        if (i + 1) % 3 == 0:
                            op4.write(value_str + '\n')
                            value_str = ''

                        value_str += '%23.16E' % A[irow, icol].imag
                        if (i + 2) % 3 == 0:
                            op4.write(value_str + '\n')
                            value_str = ''
                        i += 2
                if value_str:
                    op4.write(value_str + '\n')

        # end of the matrix?
        msg = '%8i%8i%8i\n' % (ncols + 1, 1, 1)
        msg += ' 1.0000000000000000E+00\n'
        op4.write(msg)


    def _determine_endian(self, op4):
        """Get the endian"""
        data = op4.read(8)
        (record_length_big, unused_dum_a) = unpack('>ii', data)
        (record_length_little, unused_dum_b) = unpack('<ii', data)

        # 64-bit
        record_length_big2, = unpack('>Q', data)
        record_length_little2, = unpack('<Q', data)
        if record_length_big == 24:
            endian = '>'
            self.large = False
        elif record_length_little == 24:
            endian = '<'
            self.large = False
        elif record_length_big == 48:
            endian = '<'
            self.large = True
        elif record_length_little == 48:
            endian = '<'
            self.large = True
        else:
            self.log.info('rc2 big=%s, little=%s' % (record_length_big2, record_length_little2))

            msg = 'a 24 could not be found as the first word...endian error\n'
            msg += "record_length_big=%s record_length_little=%s" % (
                record_length_big, record_length_little)
            op4.seek(0) # types = 'ifdlqILQ'
            self._show_ndata(op4, 40, types='ifdlqILQs')
            raise RuntimeError(msg)
        op4.seek(0)
        return endian

def _save_matrix(matrices, name, form, matrix):
    """save the matrix"""
    if name in matrices:
        # there are duplicate matrices with the same name (e.g., the QHH)
        form0, matrix0 = matrices[name]
        if isinstance(form0, int):
            assert isinstance(form0, int), form0
            form2 = [form0, form]
            matrix2 = [matrix0, matrix]
            matrices[name] = (form2, matrix2)
        elif isinstance(form0, list):
            form2 = form0
            matrix2 = matrix0
            form2.append(form)
            matrix2.append(matrix)
        else:  # pragma: no cover
            raise TypeError('form0=%r' % form)
    else:
        # typical case
        matrices[name] = (form, matrix)

def _get_start_end_row(A, nrows):
    """Find the starting and ending points of the matrix"""
    istart = None
    for irow in range(nrows):
        if abs(A[irow]) > 0.0:
            istart = irow
            break
    iend = None
    for irow in reversed(range(nrows)):
        if abs(A[irow]) > 0.0:
            iend = irow
            break
    return (istart, iend)

def _write_sparse_matrix_ascii(op4, name, A, form: int=2, is_big_mat: bool=False,
                               precision: str='default'):
    """
    .. todo:: Does this work for complex matrices?
    """
    msg = ''
    if isinstance(name, bytes):
        name = name.decode('ascii')
    assert isinstance(name, str), 'name=%s' % name
    #A = A.tolil() # list-of-lists sparse matrix
    #print dir(A)
    matrix_type, nwords_per_value = _get_type_nwv(A.data[0], precision)
    if matrix_type in [3, 4]:
        complex_factor = 2
    else: # 1, 2
        complex_factor = 1
    (nrows, ncols) = A.shape

    #if nrows == ncols and form == 2:
        #form = 1
    #print("matrix_type=%s" % matrix_type)
    if is_big_mat:
        msg += '%8i%8i%8i%8i%-8s1P,3E23.16\n' % (ncols, -nrows, form, matrix_type, name)
    else:
        msg += '%8i%8i%8i%8i%-8s1P,3E23.16\n' % (ncols, nrows, form, matrix_type, name)

    #print("A.row = ", A.row)
    #print("A.col = ", A.col)

    cols = {}
    for j in A.col:
        cols[j] = []
    for i, jcol in enumerate(A.col):
        cols[jcol].append(i)
    #print("cols = ", cols)

    op4.write(msg)
    msg = ''
    for j, col in cols.items():
        #print("***********")
        #print("j=%s col=%s" % (j, col))
        #col.sort()

        #print('A =', A)
        irows = [A.row[jj] for jj in col]
        #print "irows = ",irows
        dpacks = compress_column(irows)
        #print("dpacks = %s" % (dpacks))

        npacks = len(dpacks)
        nrows = len(irows)
        if is_big_mat:
            #L = complex_factor * (2 * len(irows)) + 1
            L = 2 * npacks * nwords_per_value + nrows
            msg = '%8i%8i%8i\n' % (j + 1, 0, L)
        else:
            L = complex_factor * (2 * len(irows))
            msg = '%8i%8i%8i\n' % (j + 1, 0, L + 1)
        op4.write(msg)

        for (unused_ipack, dpack) in enumerate(dpacks):
            msg = ''
            #print("pack = ",pack)

            irow = A.row[col[dpack[0]]]
            if is_big_mat:
                #L = complex_factor * (2 * len(pack)) + 1
                #L = (nPacks+1) + nRows * complex_factor
                L = (len(dpack) + 1) * nwords_per_value
                #if iPack==0:
                    #L+=1

                #L = complex_factor * (2 + npacks) + 1
                #L = len(pack) + complex_factor * 2
                #msg = '%8i%8i%8i\n' % (j+1, 0, L+1)
                msg += '%8i%8i\n' % (L, irow + 1)
            else:
                #L = complex_factor * (2 * len(pack))
                #msg = '%8i%8i%8i\n' % (j+1, 0, L+1)

                IS = irow + 65536 * (L + 1) + 1
                msg += '%8i\n' % IS

            i = 0
            value_str = ''
            #print("ipack=%s rowPack=%s" % (ipack, [A.row[p] for p in dpack]))
            for p in dpack:
                irow = col[p]
                val = A.data[irow]
                irow = A.row[irow]

                if matrix_type in [1, 2]:
                    value_str += '%23.16E' % val
                    if (i + 1) % 3 == 0:
                        msg += value_str + '\n'
                        value_str = ''
                else:
                    value_str += '%23.16E' % val.real
                    if (i + 1) % 3 == 0:
                        msg += value_str + '\n'
                        value_str = ''
                    i += 1
                    value_str += '%23.16E' % val.imag
                    if (i + 1) % 3 == 0:
                        msg += value_str + '\n'
                        value_str = ''
                i += 1
            if value_str:
                msg += value_str + '\n'
            op4.write(msg)
    op4.write('%8i%8i%8i\n' % (ncols + 1, 1, 1))
    op4.write(' 1.0000000000000000E+00\n')

def get_big_mat_nrows(nrows: int):
    """
    Parameters
    ----------
    nrows : int
        the number of rows in the matrix

    Returns
    -------
    nrows : int
        the number of rows in the matrix
    BIGMAT : Input-logical-default=FALSE. BIGMAT is applicable only
        when IUNIT < 0. BIGMAT=FALSE selects the format that uses a
        string header as described under Remark 1. But, if the
        matrix has more than 65535 rows, then BIGMAT will
        automatically be set to TRUE regardless of the value
        specified.

    """
    if nrows < 0:  # if less than 0, big
        is_big_mat = True
        nrows = abs(nrows)
    elif nrows > 0:
        is_big_mat = False
        if nrows > 65535:
            is_big_mat = True
            nrows = abs(nrows)
    else:
        raise RuntimeError('unknown BIGMAT.  nrows=%s' % nrows)
    return is_big_mat, nrows


def get_dtype(matrix_type: int, precision: str='default'):
    """Reset the type if 'default' not selected"""
    if precision == 'single':
        if matrix_type in [1, 2]:
            dtype = 'float32'
        else:
            dtype = 'complex64'
    elif precision == 'double':
        if matrix_type in [1, 2]:
            dtype = 'float64'
        else:
            dtype = 'complex128'
    else:  # default
        if matrix_type == 1:
            dtype = 'float32'
        elif matrix_type == 2:
            dtype = 'float64'
        elif matrix_type == 3:
            dtype = 'complex64'
        else:
            dtype = 'complex128'
    return dtype


def _get_type_nwv(A, precision: str='default'):
    """
    Determines the Type and number of words per value
    an entry in the matrix takes up.

    Parameters
    ----------
    A : matrix
        a matrix or entry in a matrix (to save memory)
    precision : str
        data precision ='default', 'single', 'double'

    Returns
    -------
    matrix_type : int
        the dtype of the matrix as an integer
    nwords_per_value : int
        Number of words per value

    A word is 4 bytes

    +-------------+------------+-----------+--------+--------+
    | Matrix Type |   dtype    | precision | nwords | nbytes |
    +=============+============+===========+========+========+
    |      1      | float32    |   single  |    1   |    4   |
    +-------------+------------+-----------+--------+--------+
    |      2      | complex64  |   single  |    2   |    8   |
    +-------------+------------+-----------+--------+--------+
    |      3      | float64    |   double  |    2   |    8   |
    +-------------+------------+-----------+--------+--------+
    |      4      | complex128 |   double  |    4   |   16   |
    +-------------+------------+-----------+--------+--------+

    """
    # real
    if isinstance(A.dtype.type(), float32):
        nwords_per_value = 1
        if precision != 'double':
            matrix_type = 1
        else:
            matrix_type = 2
    elif isinstance(A.dtype.type(), float64):
        nwords_per_value = 1
        if precision != 'single':
            matrix_type = 2
        else:
            matrix_type = 1

    # complex
    elif isinstance(A.dtype.type(), complex64):
        nwords_per_value = 2
        if precision != 'double':
            matrix_type = 3
        else:
            matrix_type = 4
    elif isinstance(A.dtype.type(), complex128):
        nwords_per_value = 2
        if precision != 'single':
            matrix_type = 4
        else:
            matrix_type = 3
    else:
        msg = ('invalid matrix_type, only float32, float64, '
               'complex64, complex128; dtype=%r' % A.dtype)
        raise TypeError(msg)
    return matrix_type, nwords_per_value


def compress_column(col):
    """takes a dense matrix column and puts it into OP4 format"""
    packs = []

    n = 0
    i = 0
    packi = []
    while i < len(col):
        #print("i=%s n=%s col[i]=%s" % (i, n, col[i]))
        if col[i] == n + 1:
            #print("i=n=%s" % i)
            packi.append(i)
            n += 1
        else:
            if packi:
                packs.append(packi)
                #print("pack = ", pack)
            packi = [i]
            n = col[i]
        #print("pack = ", pack)
        i += 1

    if packi:
        packs.append(packi)
    #print("packs = ", packs)
    return packs
