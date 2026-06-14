"""Main OP4 class"""

import sys
import os
from struct import pack, unpack, Struct
from typing import TextIO, BinaryIO, Optional, Any, cast

import numpy as np
from numpy import float32, float64, complex64, complex128
from scipy.sparse import coo_matrix  # type: ignore
from cpylog import SimpleLogger, __version__ as CPYLOG_VERSION

from pyNastran.utils import is_binary_file as file_is_binary, PathLike, PurePath
from pyNastran.utils.mathematics import print_matrix #, print_annotated_matrix
from pyNastran.op2.result_objects.matrix import Matrix
if CPYLOG_VERSION > '1.6.0':
    from cpylog import get_logger
else:  # pragma: no cover
    from cpylog import get_logger2 as get_logger


class EmptyMatrixError(RuntimeError):
    """Raised when an OP4 matrix header has nrows=0 (empty/uninitialized)."""


class OP4:
    """
    todo:: add endian checking
    todo:: test on big matrices
    todo:: finish write_op4

    """
    def __init__(self, log=None, debug: bool=False):
        self.n = 0
        self._endian = ''
        self.debug = debug
        #assert debug == True, debug
        self.log = get_logger(log, debug)
        self._new = False
        self.large = False

    def read_op4(self, op4_filename: Optional[PathLike]=None,
                 matrix_names: Optional[list[str]]=None,
                 precision: str='default') -> dict[str, Matrix]:
        """See ``read_op4``"""
        if precision not in {'default', 'single', 'double'}:
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
            matrices = self.read_op4_binary(
                op4_filename, matrix_names, precision)
        else:
            matrices = self.read_op4_ascii(
                op4_filename, matrix_names, precision)
        return matrices

#--------------------------------------------------------------------------
    def read_op4_ascii(self, op4_filename: PathLike,
                       matrix_names: Optional[list[str]]=None,
                       precision: str='default') -> dict[str, Matrix]:
        """matrix_names must be a list or None, but basically the same"""
        matrices: dict[str, Matrix] = {}
        name = 'dummyName'
        with open(op4_filename, 'r') as op4:
            while name is not None:
                name, amat = self._read_matrix_ascii(op4, matrix_names, precision)
                if name is None:
                    assert amat is None
                    break

                assert isinstance(amat, Matrix), amat
                if is_saved_matrix(name, matrix_names):
                    assert amat is not None, amat
                    _save_matrix(matrices, name, amat)
        return matrices

    def _read_matrix_ascii(self, op4: TextIO,
                           matrix_names: Optional[list[str]]=None,
                           precision: str='default') -> tuple[Optional[str], Optional[Matrix]]:
        """Reads an ASCII matrix"""
        iline = 0
        line = op4.readline().rstrip()
        iline += 1
        if line == '':
            op4.close()
            return None, None
        ncols_str, nrows_str, form_str, matrix_type_str = line[0:32].split()
        nrows = int(nrows_str)

        try:
            is_big_mat, nrows = get_big_mat_nrows(nrows)
        except EmptyMatrixError:
            return None, None
        if self.debug:
            self.log.info('is_big_matrix = %s' % is_big_mat)

        ncols = int(ncols_str)
        form = int(form_str)
        matrix_type = int(matrix_type_str)
        dtype = get_dtype(matrix_type, precision)

        name = line[32:40].strip()

        if self.debug:
            self.log.info('name=%s shape=(%s,%s) form=%s Type=%s' % (
                name, nrows, ncols, form, matrix_type))
        assert ncols > 0, 'ncols=%s' % ncols
        size = line[40:].strip()
        line_size_str = size.split(',')[1].split('E')[1].split('.')[0]  # 3E23.16 to 23
        line_size = int(line_size_str)

        line = op4.readline().rstrip()
        iline += 1
        (_icol, irow, _nwords) = line.split()

        is_sparse = False
        if irow == '0':
            is_sparse = True

        if matrix_type in {1, 2}:  # real
            data_mat, iline = self._read_real_ascii(op4, iline, nrows, ncols, line_size, line,
                                                dtype, is_sparse, is_big_mat)
        elif matrix_type in {3, 4}:  # complex
            if is_sparse:
                data_mat, iline = self._read_complex_sparse_ascii(op4, iline, nrows, ncols,
                                                              line_size, line, dtype, is_big_mat)
            else:
                data_mat, iline = self._read_complex_dense_ascii(op4, iline, nrows, ncols,
                                                             line_size, line, dtype, is_big_mat)
        else:
            raise RuntimeError('invalid matrix type.  matrix_type=%d' % matrix_type)

        if self.debug:
            self.log.info("form=%s name=%s data_mat=\n%s" % (form, name, str(data_mat)))
        amat = Matrix(name, form, data=data_mat)
        return name, amat

    def _read_real_sparse_ascii(self, op4: TextIO, iline: int, nrows: int, ncols: int,
                                line_size: int, line: str, dtype: str,
                                is_big_mat: bool) -> tuple[coo_matrix, int]:
        """Reads a sparse real ASCII matrix"""
        self.log.debug('_read_real_sparse_ascii')
        rows = []
        cols = []
        entries: list[str] = []
        nloops = 0
        was_broken = False
        while 1:
            if nloops > 0 and not was_broken:
                line = op4.readline().rstrip()
                iline += 1
            was_broken = False

            icol_str, irow_str, nwords_str = line.split()
            icol = int(icol_str)

            if icol > ncols:
                break

            irow = int(irow_str)
            nwords = int(nwords_str)

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

        rows_array = np.array(rows, dtype='int32') - 1
        cols_array = np.array(cols, dtype='int32') - 1
        data_mat = coo_matrix((entries, (rows_array, cols_array)), shape=(nrows, ncols), dtype=dtype)
        #print("type = %s %s" % (type(data_mat),type(data_mat.toarray())))
        #data_mat = data_mat.toarray()
        return data_mat, iline

    def _read_real_sparse_ascii_new(self, op4: TextIO,
                                    iline: int, nrows: int, ncols: int,
                                    line_size: int, line: str,
                                    dtype: str, is_big_mat: bool) -> tuple[coo_matrix, int]:
        """Reads a sparse real ASCII matrix"""
        self.log.debug('_read_real_sparse_ascii')
        rows: list[int] = []
        cols: list[list[int]] = []
        entries: list[str] = []
        nloops = 0
        was_broken = False
        while 1:
            if nloops > 0 and not was_broken:
                line = op4.readline().rstrip()
                iline += 1
            was_broken = False

            icol_str, irow_str, nwords_str = line.split()
            icol = int(icol_str)

            if icol > ncols:
                break

            irow = int(irow_str)
            nwords = int(nwords_str)

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

                    irows: list[int] = list(range(irow, irow + nwords_in_line))
                    n = 0
                    for unused_i in range(nwords_in_line):
                        word: str = line[n:n + line_size]
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

        cols_array = np.hstack(cols)
        rows_array = np.array(rows, dtype='int32') - 1
        cols_array = np.array(cols_array, dtype='int32') - 1
        data_mat = coo_matrix((entries, (rows_array, cols_array)), shape=(nrows, ncols), dtype=dtype)
        #print("type = %s %s" % (type(data_mat), type(data_mat.toarray())))
        #data_mat = data_mat.toarray()
        return data_mat, iline

    def _read_real_dense_ascii(self, op4: TextIO, iline: int, nrows: int, ncols: int,
                               line_size: int, line: str, dtype: str,
                               is_big_mat: bool) -> tuple[np.ndarray, int]:
        """Reads a real dense ASCII matrix"""
        self.log.debug('_read_real_dense_ascii')
        data_mat = np.zeros((nrows, ncols), dtype=dtype)  # Initialize a real matrix
        nloops = 0
        was_broken = False
        while 1:
            if nloops > 0 and not was_broken:
                line = op4.readline().rstrip()
                iline += 1
            was_broken = False

            (icol_str, irow_str, nwords_str) = line.split()
            icol = int(icol_str)

            if icol > ncols:
                break

            irow = int(irow_str)
            nwords = int(nwords_str)

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
                        data_mat[irow - 1, icol - 1] = word
                        n += line_size
                        irow += 1
                    #iword += nwords_in_line
                    nwords -= nwords_in_line
                sline = line.strip().split()
                nloops += 1
        op4.readline()
        iline += 1
        return data_mat, iline

    def _read_real_ascii(self, op4: TextIO, iline: int, nrows: int, ncols: int,
                         line_size: int, line: str, dtype: str,
                         is_sparse: bool, is_big_mat: bool) -> tuple[np.ndarray, int]:
        """Reads a real ASCII matrix"""
        if is_sparse:
            if self._new:
                data_mat, iline = self._read_real_sparse_ascii_new(op4, iline, nrows, ncols,
                                                               line_size, line, dtype, is_big_mat)
            else:
                data_mat, iline = self._read_real_sparse_ascii(op4, iline, nrows, ncols,
                                                           line_size, line, dtype, is_big_mat)
        else:
            data_mat, iline = self._read_real_dense_ascii(op4, iline, nrows, ncols,
                                                      line_size, line, dtype, is_big_mat)
        return data_mat, iline

    def _read_complex_sparse_ascii(self, op4: TextIO, iline: int, nrows: int, ncols: int,
                                   line_size: int, line: str,
                                   dtype: str, is_big_mat: bool) -> tuple[coo_matrix, int]:
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

            (icol_str, irow_str, nwords_str) = line.split()
            icol = int(icol_str)

            if icol > ncols:
                break

            irow = int(irow_str)
            nwords = int(nwords_str)

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

        rows_array = np.array(rows, dtype='int32') - 1
        cols_array = np.array(cols, dtype='int32') - 1
        data_mat = coo_matrix((entries, (rows_array, cols_array)), shape=(nrows, ncols), dtype=dtype)
        op4.readline()
        iline += 1
        return data_mat, iline

    def _read_complex_dense_ascii(self, op4: TextIO, iline: int, nrows: int, ncols: int,
                                  line_size: int, line: str,
                                  dtype: str, is_big_mat: bool) -> tuple[np.ndarray, int]:
        """Reads a dense complex ASCII matrix"""
        data_mat = np.zeros((nrows, ncols), dtype=dtype)  # Initialize a complex matrix

        nloops = 0
        was_broken = False
        while 1:
            if nloops > 0 and not was_broken:
                line = op4.readline().rstrip()
                iline += 1
            was_broken = False

            (icol_str, irow_str, nwords_str) = line.split()
            icol = int(icol_str)

            if icol > ncols:
                break

            irow = int(irow_str)
            nwords = int(nwords_str)

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
                            #data_mat[irow - 1, icol - 1].real = value
                            real_value = value
                        else:
                            data_mat[irow - 1, icol - 1] = complex(real_value, value)
                            irow += 1
                        iword += 1
                        n += line_size
                    nwords -= nwords_in_line
                sline = line.strip().split()
                iline += 1
                nloops += 1

        op4.readline()
        iline += 1
        return data_mat, iline

    def _get_irow_small_ascii(self, op4: TextIO, iline: int, line: str, sline: list[str],
                              irow: int) -> tuple[int, int]:
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

    def _get_irow_small_binary(self, op4: BinaryIO, data_bytes: bytes) -> tuple[int, int]:
        """
        Returns
        -------
        irow : int
           the row id
        L : int
            the row length

        """
        if len(data_bytes) == 0:
            data_bytes = op4.read(4)
            self.n += 4

        IS, = unpack(self._endian + 'i', data_bytes)
        L = IS // 65536 - 1
        irow = IS - 65536 * (L + 1)
        if self.debug:
            self.log.info('small_mat-next row')
            self.log.info("  IS=%s L=%s irow=%s" % (IS, L, irow))
            assert IS > 0, IS
            assert L > 0, L
        return irow, L

    def _get_irow_big_ascii(self, op4: TextIO, iline: int,
                            line: str,
                            sline: list[str], irow: int) -> tuple[int, int]:
        sline = line.strip().split()
        if len(sline) == 2:
            pass
        else:
            sline = op4.readline().strip().split()
            iline += 1
        assert len(sline) == 2, 'sline=%s len(sline)=%s' % (sline, len(sline))
        (idummy, irow_str) = sline
        irow = int(irow_str)
        if self.debug:
            self.log.debug("idummy=%s irow=%s" % (idummy, irow))
        return irow, iline

    def _get_irow_big_binary(self, op4: BinaryIO, data_bytes: bytes) -> tuple[int, int]:
        """
        Returns
        -------
        irow : int
           the row id
        L : int
            ???

        """
        if len(data_bytes) == 0:
            data_bytes = op4.read(8)
            self.n += 8
        idummy, irow = unpack(self._endian + '2i', data_bytes)
        if self.debug:
            self.log.debug("idummy=%s irow=%s" % (idummy, irow))
            assert irow < 100, irow
        return irow, idummy - 1

#--------------------------------------------------------------------------
    def read_op4_binary(self, op4_filename: PathLike,
                        matrix_names: Optional[list[str]]=None,
                        precision: str='default',
                       use_matrix_class=False):
        """matrix_names must be a list or None, but basically the same"""
        self.n = 0
        matrices: dict[str, Matrix] = {}
        name = 'dummyName'

        with open(op4_filename, mode='rb') as op4:
            self._endian = self._determine_endian(op4)
            while name is not None:
                # checks for the end of the file
                assert self.n == op4.tell(), 'n=%s tell=%s' % (self.n, op4.tell())
                n = self.n
                data1 = op4.read(1)
                op4.seek(n)
                if len(data1) == 0:
                    break
                #self.show(f, 60)

                (name, amat) = self._read_matrix_binary(op4, precision, matrix_names)
                if name is None:
                    break
                #print(print_matrix(amat.matrix))
                if is_saved_matrix(name, matrix_names):
                    _save_matrix(matrices, name, amat)

                #print("not op4.closed = ",not op4.closed,form,name)
                # if not op4.closed or form is not None:
                #     data_bytes = op4.read(4)
                #     self.n += 4
                #     if len(data_bytes) == 0:
                #         break
                #     (record_length,) = unpack(self._endian + 'i', data_bytes)
                ##     print("record_length = %s" % record_length)
                #     if record_length == 24:
                #         self.n -= 4
                #         op4.seek(self.n)
                #     else:
                #         data_bytes = op4.read(4)
                #         if len(data_bytes) == 0:
                #             break
                #         (record_length2,) = unpack(self._endian + 'i', data_bytes)
                #         assert record_length2 == 24
                #         op4.seek(self.n)
                #
        return matrices

    def read_start_marker(self, op4: BinaryIO) -> tuple[int, int, int, int]:
        if self.debug:
            self.log.info('--------------------------------------')
        #self.show(op4, 60)
        data_bytes: bytes = op4.read(4)
        self.n += 4
        record_length, = unpack(self._endian + 'i', data_bytes)
        #print('record_length =', record_length)

        record_length = 16
        data_bytes = op4.read(record_length)
        self.n += record_length
        assert self.n == op4.tell(), 'n=%s tell=%s' % (self.n, op4.tell())

        if record_length == 16:
            a, icol, irow, nwords = unpack(self._endian + '4i', data_bytes)
            if self.debug:
                self.log.info("a=%s icol=%s irow=%s nwords=%s" % (a, icol, irow, nwords))
        else:
            raise NotImplementedError('record_length=%s' % record_length)
        return a, icol, irow, nwords

    def _read_matrix_binary(self, op4: BinaryIO, precision: str,
                            matrix_names: list[str]) -> tuple[str, Matrix]:
        """Reads a binary matrix"""
        #self.show(f, 60)
        log = self.log
        if self.debug:
            log.info("*************************")
        data_bytes = op4.read(4)
        self.n += 4
        (record_length,) = unpack(self._endian + 'i', data_bytes)
        assert self.n == op4.tell(), 'n=%s tell=%s' % (self.n, op4.tell())
        if self.debug:
            log.info("record_length = %s" % record_length)

        if record_length == 24:
            fmt = self._endian + '4i8s'
        elif record_length == 48:
            fmt = self._endian + '4Q16s'
        else:
            #msg = record_length #+ self.print_block(data_bytes)
            msg = 'record_length=%s filename=%r' % (record_length, op4.name)
            raise NotImplementedError(msg)

        data_bytes = op4.read(record_length)
        self.n += record_length
        (ncols, nrows, form, matrix_type, name) = unpack(fmt, data_bytes)
        if self.debug:
            log.info("nrows=%s ncols=%s form=%s matrix_type=%s name=%r" % (
                nrows, ncols, form, matrix_type, name))

        name = name.strip()
        name = name.decode('ascii')
        if self.debug:
            if matrix_type == 1:
                log.info("matrix_type = Real, Single Precision")
            elif matrix_type == 2:
                log.info("matrix_type = Real, Double Precision")
            elif matrix_type == 3:
                log.info("matrix_type = Complex, Single Precision")
            elif matrix_type == 4:
                log.info("matrix_type = Complex, Double Precision")

        try:
            is_big_mat, nrows = get_big_mat_nrows(nrows)
        except EmptyMatrixError:
            return None, None

        if self.debug:
            log.info('is_big_matrix = %s' % is_big_mat)

        # jump forward to get irow (needed for check on is_sparse),
        # then jump back
        nsave = self.n
        irow = self.read_start_marker(op4)[2]
        op4.seek(nsave)
        self.n = nsave

        #(nwords_per_value, nbytes_per_value, data_format, dtype) = self._get_matrix_info(matrix_type)
        data_format, dtype = _get_matrix_info(matrix_type, self.log, debug=self.debug)[2:]

        is_sparse = False
        if irow == 0:
            is_sparse = True

        assert self.n == op4.tell(), 'n=%s tell=%s' % (self.n, op4.tell())
        if matrix_type in {1, 2}:  # real
            if is_sparse:
                data_mat = self._read_real_sparse_binary(op4, nrows, ncols, matrix_type, is_big_mat)
            else:
                data_mat = self._read_real_dense_binary(op4, nrows, ncols, matrix_type, is_big_mat)
        elif matrix_type in {3, 4}:  # complex
            if is_sparse:
                data_mat = self._read_complex_sparse_binary(op4, nrows, ncols, matrix_type, is_big_mat)
            else:
                data_mat = self._read_complex_dense_binary(op4, nrows, ncols, matrix_type, is_big_mat)
        else:
            log.error('is_sparse=%s data_format=%s dtype=%s' % (is_sparse, data_format, dtype))
            raise TypeError(f'matrix_type={matrix_type}')

        #try:
            #print_matrix(A.toarray())
        #except Exception:
            #pass

        if data_format in {'d', 'dd'}:
            op4.read(8)
            self.n += 8
        elif data_format in {'f', 'ff'}:
            op4.read(4)
            self.n += 4
        else:
            raise NotImplementedError(data_format)
        #f.read(record_length); self.n+=record_length
        #self.show(f, 10)
        #f.read(4); self.n+=4

        assert self.n == op4.tell(), 'n=%s op4.tell=%s' % (self.n, op4.tell())
        amat = Matrix(name, form, data=data_mat)
        return name, amat

    def _read_real_dense_binary(self, op4: BinaryIO, nrows: int, ncols: int,
                                matrix_type: int, is_big_mat: bool) -> np.ndarray:
        if self.debug:
            self.log.info('_read_real_dense_binary')
        out = _get_matrix_info(matrix_type, self.log, debug=False)
        (nwords_per_value, _nbytes_per_value, data_format, dtype) = out
        data_mat = np.zeros((nrows, ncols), dtype=dtype)

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
            data_bytes = op4.read(record_length)
            self.n += record_length
            nvalues = L // nwords_per_value
            str_values = self._endian + '%i%s' % (nvalues, data_format)
            data_mat[irow-1:irow-1+nvalues, icol-1] = unpack(str_values, data_bytes)
            if self.debug:
                self.log.info('A[%s:%s, %s] = %s' % (
                    irow - 1,
                    irow - 1 + nvalues,
                    icol-1,
                    data_mat[irow-1:irow-1+nvalues, icol-1]))
        #assert self.n == op4.tell(), 'n=%s tell=%s' % (self.n, op4.tell())
        op4.read(4)
        self.n += 4
        return data_mat

    def _read_real_sparse_binary(self, op4: BinaryIO,
                                 nrows: int, ncols: int, matrix_type: int,
                                 is_big_mat: bool) -> coo_matrix:
        if self.debug:
            self.log.info('_read_real_sparse_binary')
        #self._show(op4, 200, types='ifsdq', endian=None)
        out = _get_matrix_info(matrix_type, self.log, debug=False)
        (nwords_per_value, nbytes_per_value, data_format, dtype) = out
        log = self.log
        rows = []
        cols = []
        entries: list[float] = []

        data_bytes = b''
        icol = -1  # dummy value so the loop starts
        while icol < ncols + 1:  # if isDense
            #if record_length==0:
            assert self.n == op4.tell(), 'n=%s tell=%s' % (self.n, op4.tell())
            (icol, irow, nwords) = self.get_markers_sparse(op4, is_big_mat)
            L = nwords

            if icol == ncols + 1:
                if self.debug:
                    log.info('breaking on icol=%s ncol+1=%s' % (icol, ncols + 1))
                break

            if is_big_mat:
                irow, L = self._get_irow_big_binary(op4, data_bytes[:8])
                data_bytes = data_bytes[8:]
            else:
                irow, L = self._get_irow_small_binary(op4, data_bytes[:4])
                data_bytes = data_bytes[4:]

            if L == -1:
                if self.debug:
                    log.info('breaking on L=-1')
                break

            if self.debug:
                log.info("  next icol")
                log.info("    n=%s icol=%s irow=%s nwords=%s" % (
                    self.n, icol, irow, nwords))
                self._show(op4, 100, types='qd')
                log.info('**************************************************')

            #if nwords == 0 and is_big_mat:
                #self.n -= 4
                #op4.seek(self.n)
                #break

            record_length = 4 * nwords
            data_bytes = op4.read(record_length)
            self.n += record_length
            if self.debug:
                log.info("  data_format=%s record_length=%s n_next=%s" % (
                    data_format, record_length, self.n))
            #if icol == ncols + 1:
                #break

            i = 0
            while len(data_bytes) > 0:
                if i > 0:
                    if is_big_mat:
                        (irow, L) = self._get_irow_big_binary(op4, data_bytes[0:8])
                        data_bytes = data_bytes[8:]
                    else:
                        (irow, L) = self._get_irow_small_binary(op4, data_bytes[0:4])
                        data_bytes = data_bytes[4:]
                    assert irow > 0
                nvalues = L // nwords_per_value
                str_values = self._endian + '%i%s' % (nvalues, data_format)

                if self.debug:
                    log.info('irow=%s L=%s nwords_per_value=%s nvalues=%s '
                             'nbytes_per_value=%s' % (irow, L, nwords_per_value,
                                                      nvalues, nbytes_per_value))
                    log.info('str_values = %r' % str_values)

                value_list = unpack(str_values, data_bytes[0:nvalues * nbytes_per_value])
                assert self.n == op4.tell(), 'n=%s tell=%s' % (self.n, op4.tell())

                #irow -= 1
                #icol -= 1
                if self.debug:
                    log.info('rows = %s' % list(i+irow-1 for i in range(nvalues)))
                    log.info('cols = %s ' % ([icol-1] * nvalues))
                    log.info('value_list = %s' % str(value_list))

                rows.extend([i+irow-1 for i in range(nvalues)])
                irow += nvalues
                cols.extend([icol-1] * nvalues)
                entries.extend(value_list)

                record_length -= nvalues * nbytes_per_value
                data_bytes = data_bytes[nvalues * nbytes_per_value:]
                if self.debug:
                    log.info("  record_length=%s nbytes_per_value=%s len(data)=%s" % (
                        record_length, nbytes_per_value, len(data_bytes)))
                    ##print(A)
                    #print("********")  # ,data_bytes
                    #print(self.print_block(data_bytes))
                i += 1
            #print "-------------------------------"

        #if rows == []:  # NULL matrix
            #raise NotImplementedError()

        A = coo_matrix((entries, (rows, cols)), shape=(nrows, ncols),
                       dtype=dtype)
        op4.read(4)
        self.n += 4
        return A

    def _show(self, op4: BinaryIO, n, types: str='ifs', endian: Optional[str]=None):
        """Shows binary data"""
        assert self.n == op4.tell()
        nints = n // 4
        data_bytes = op4.read(4 * nints)
        strings, ints, floats = self._show_data(data_bytes, types=types, endian=endian)
        op4.seek(self.n)
        return strings, ints, floats

    def _show_data(self, data_bytes: bytes, types: str='ifs', endian: Optional[str]=None):
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
        if endian is None:
            endian = self._endian
        return _write_data(sys.stdout, data_bytes, types=types, endian=endian)

    def _show_ndata(self, f: BinaryIO, n: int, types: str='ifs') -> None:
        #endian = self._endian
        return self._write_ndata(sys.stdout, f, n, types=types)

    def _write_ndata(self, fout: TextIO, f: BinaryIO, n: int, types: str='ifs') -> None:
        """Useful function for seeing what's going on locally when debugging."""
        endian = self._endian
        assert endian is not None, endian
        nold = self.n
        data_bytes = f.read(n)
        self.n = nold
        f.seek(self.n)
        return _write_data(fout, data_bytes, endian=endian, types=types)

    def _read_complex_dense_binary(self, op4: BinaryIO, nrows: int, ncols: int,
                                   matrix_type: int, is_big_mat: bool) -> coo_matrix:
        """reads a dense complex binary matrix"""
        if self.debug:
            self.log.info('_read_complex_dense_binary')
        out = _get_matrix_info(matrix_type, self.log, debug=False)
        (nwords_per_value, nbytes_per_value, data_format, dtype) = out

        A = np.zeros((nrows, ncols), dtype=dtype)
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
            data_bytes = op4.read(record_length)
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
                        nvalues * nbytes_per_value, len(data_bytes)))
                value_list = unpack(str_values, data_bytes[0:nvalues * nbytes_per_value])
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
                data_bytes = data_bytes[nvalues * nbytes_per_value:]
                if self.debug:
                    self.log.info("record_length=%s nbytes_per_value=%s" % (
                        record_length, nbytes_per_value))
                    self.log.info(print_matrix(A))
                    self.log.info("******** %r" % data_bytes)

        op4.read(4)
        self.n += 4
        return A

    #def _read_complex_binary(self, op4: BinaryIO, nrows: int, ncols: int,
                             #matrix_type: int, is_sparse: bool, is_big_mat: bool) -> coo_matrix:
        #"""Reads a complex binary matrix"""
        #if is_sparse:
            #A = self._read_complex_sparse_binary(op4, nrows, ncols, matrix_type, is_big_mat)
        #else:
            #A = self._read_complex_dense_binary(op4, nrows, ncols, matrix_type, is_big_mat)
        #return A

    def _read_complex_sparse_binary(self, op4: BinaryIO, nrows: int, ncols: int,
                                    matrix_type: int, is_big_mat: bool) -> coo_matrix:
        """Reads a sparse complex binary matrix"""
        if self.debug:
            self.log.info('_read_complex_sparse_binary')
        out = _get_matrix_info(matrix_type, self.log, debug=False)
        (nwords_per_value, nbytes_per_value, data_format, dtype) = out
        rows = []
        cols = []
        entries = []
        record_length = 0
        data_bytes = b''
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
                (irow, L) = self._get_irow_big_binary(op4, data_bytes[:8])
                data_bytes = data_bytes[8:]
            else:
                (irow, L) = self._get_irow_small_binary(op4, data_bytes[:4])
                data_bytes = data_bytes[4:]

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
            data_bytes = op4.read(record_length)
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
                    self.log.info("  nvalues*nbytes_per_value=%s len(data_bytes)=%s" % (
                        nvalues * nbytes_per_value, len(data_bytes)))
                value_list = unpack(str_values, data_bytes[0:nvalues * nbytes_per_value])
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
                data_bytes = data_bytes[nvalues * nbytes_per_value:]
                #print("record_length=%s nbytes_per_value=%s" % (record_length, nbytes_per_value))
                #print(print_matrix(A))
                #print("********", data)

        data_mat = coo_matrix((entries, (rows, cols)), shape=(nrows, ncols), dtype=dtype)
        op4.read(4)
        self.n += 4
        return data_mat

    def get_markers_sparse(self, op4: BinaryIO, is_big_mat: bool) -> tuple[int, int, int]:
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

    def get_markers_dense(self, op4: BinaryIO) -> tuple[int, int, int]:
        a, icol, irow, nwords = self.read_start_marker(op4)
        if self.debug:
            self.log.info("n=%s a=%s icol=%s irow=%s nwords=%s"% (
                self.n, a, icol, irow, nwords))
        return icol, irow, nwords

    def write_op4(self, op4_filename: Optional[PathLike],
                  matrices: dict[str, Matrix],
                  name_order=None,
                  precision: str='default',
                  is_binary: bool=True) -> None:
        """
        Writes the OP4

        Parameters
        ----------
        op4_filename : str/file
            The filename to write
            String -> opens a file (closed at the end)
            file   -> no file is opened and it's not closed
        matrices : dict[str] = (form, np.ndarray)
            the matrices to write

        name_order: str / list[str]; default=None -> sorted based on name
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
        >>> write_op4(op4_filename, matrices, name_order=['A', 'B', 'C'],
                      precision='default', is_binary=True)

        # another method
        >>> matrices = {
            'A' : Matrix('A', formA, data=matrixA),
            'B' : (formB, matrixB),
            'C' : (formC, matrixC),
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

        if op4_filename is None:
            from pyNastran.utils.gui_io import save_file_dialog
            wildcard_wx = "Nastran OP4 (*.op4)|*.op4|" \
                "All files (*.*)|*.*"
            wildcard_qt = "Nastran OP4 (*.op4);;All files (*)"
            title = 'Please select a OP4 to save'
            op4_filename = save_file_dialog(title, wildcard_wx, wildcard_qt)
            assert op4_filename is not None, op4_filename

        name_order2 = _prepare_name_order(matrices, name_order)
        if isinstance(op4_filename, (str, PurePath)):
            if is_binary:
                with open(op4_filename, 'wb') as op4b:
                    self._write_op4_file_binary(op4b, name_order2, precision, matrices)
            else:
                with open(op4_filename, 'w') as op4:
                    self._write_op4_file_ascii(op4, name_order2, precision, matrices)
        else:
            op4 = op4_filename
            if is_binary:
                op4 = cast(BinaryIO, op4)
                self._write_op4_file_binary(op4, name_order2, precision, matrices)
            else:
                op4 = cast(TextIO, op4)
                self._write_op4_file_ascii(op4, name_order2, precision, matrices)

    def _write_op4_file_ascii(self, op4: TextIO, name_order: list[str],
                              precision: str,
                              matrices: dict[str, Matrix]) -> None:
        """Helper method for OP4 writing"""
        for name in name_order:
            form, matrix = _write_form_matrix_helper(matrices, name)

            if isinstance(matrix, coo_matrix):
                nrows = matrix.shape[0]
                is_big_mat = (nrows > 65535)
                _write_sparse_matrix_ascii(
                    op4, name, matrix, form=form,
                    precision=precision, is_big_mat=is_big_mat)
            elif isinstance(matrix, np.ndarray):
                _write_dense_matrix_ascii(
                    self.log, op4, name, matrix, form=form, precision=precision, debug=self.debug)
            else:
                msg = ('Matrix type=%r is not supported.  '
                       'types=[coo_matrix, ndarray]' % type(matrix))
                raise NotImplementedError(msg)

    def _write_op4_file_binary(self, op4: BinaryIO,
                               name_order: list[str],
                               precision: str,
                               matrices: dict[str, Matrix]) -> None:
        """Helper method for OP4 writing"""

        #is_big_mat = False  ## .. todo:: hardcoded
        for name in name_order:
            form, matrix = _write_form_matrix_helper(matrices, name)

            if isinstance(matrix, coo_matrix):
                _write_sparse_matrix_binary(
                    op4, name, matrix, form=form, precision=precision, endian=self._endian)
            elif isinstance(matrix, np.ndarray):
                _write_dense_matrix_binary(
                    op4, name, matrix, form=form, precision=precision, endian=self._endian)
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

    def _determine_endian(self, op4: BinaryIO) -> str:
        """Get the endian"""
        data_bytes = op4.read(8)
        (record_length_big, unused_dum_a) = unpack('>ii', data_bytes)
        (record_length_little, unused_dum_b) = unpack('<ii', data_bytes)

        # 64-bit
        record_length_big2, = unpack('>Q', data_bytes)
        record_length_little2, = unpack('<Q', data_bytes)
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
            self._show_ndata(op4, 80, types='ifdlqILQs')
            raise RuntimeError(msg)
        op4.seek(0)
        return endian

def _save_matrix(matrices, name: str, amat: Matrix) -> None:
    """save the matrix"""
    assert isinstance(name, str), name
    assert isinstance(amat, Matrix), amat

    if name in matrices:
        # there are duplicate matrices with the same name (e.g., the QHH)
        amat0: Matrix = matrices[name]
        form0 = amat0.form
        data0 = amat0.data

        if isinstance(form0, int):
            assert isinstance(form0, int), form0
            form2 = [form0, amat.form]
            data2 = [data0, amat.data]
            amat0.form = form2
            amat0.data = data2
        elif isinstance(form0, list):
            form0.append(amat.form)
            data0.append(amat.data)
        else:  # pragma: no cover
            raise TypeError('form0=%r' % form0)
    else:
        # typical case
        matrices[name] = amat

def _get_start_end_row(A: np.ndarray, nrows: int) -> tuple[Optional[int], Optional[int]]:
    """Find the starting and ending points of the matrix"""
    nz = np.flatnonzero(A)
    if len(nz) == 0:
        return None, None
    return int(nz[0]), int(nz[-1])


def _fmt_values_e23(values: np.ndarray) -> str:
    """Format a 1-D float array as lines of 3 values in '%23.16E' format."""
    n = len(values)
    if n == 0:
        return ''
    # np.array2string / savetxt are slow; manual vectorization with format_float_scientific
    # is also slow. The fastest approach is a format string with numpy's vectorized str conversion.
    parts = []
    for i in range(0, n, 3):
        chunk = values[i:i+3]
        parts.append(''.join('%23.16E' % v for v in chunk))
    return '\n'.join(parts) + '\n'


def _write_dense_matrix_ascii(log: SimpleLogger,
                              op4: TextIO, name: str, A: np.ndarray,
                              form: int=2,
                              precision: str='default',
                              debug: bool=False) -> None:
    """Writes a dense ASCII matrix using vectorized column start/end detection and buffered output."""
    if debug:
        log.info('_write_dense_matrix_ascii')
    matrix_type, nwords_per_value = _get_type_nwv(A[0, 0], precision)

    (nrows, ncols) = A.shape
    # big_mat: write negative nrows to signal nrows > 65535
    nrows_header = -nrows if nrows > 65535 else nrows
    msg = '%8i%8i%8i%8i%-8s1P,3E23.16\n' % (ncols, nrows_header, form, matrix_type, name)
    op4.write(msg)

    is_complex = matrix_type in (3, 4)

    # Ensure column-major for efficient column slicing
    if not A.flags['F_CONTIGUOUS']:
        A = np.asfortranarray(A)

    # Vectorized column start/end computation
    nonzero_mask = (A != 0)
    has_data = nonzero_mask.any(axis=0)
    col_starts = np.argmax(nonzero_mask, axis=0)
    col_ends = nrows - 1 - np.argmax(nonzero_mask[::-1, :], axis=0)

    # Buffer output in chunks to reduce write() calls
    buf_parts: list[str] = []
    buf_size = 0
    flush_threshold = 256 * 1024  # 256KB buffer before flushing

    for icol in range(ncols):
        if not has_data[icol]:
            continue

        istart = int(col_starts[icol])
        iend = int(col_ends[icol]) + 1
        nvalues = iend - istart

        # Column header
        buf_parts.append('%8i%8i%8i\n' % (icol + 1, istart + 1,
                                           nvalues * nwords_per_value))

        # Format values
        if is_complex:
            segment = A[istart:iend, icol]
            floats = np.empty(nvalues * 2, dtype=np.float64)
            floats[0::2] = segment.real
            floats[1::2] = segment.imag
        else:
            floats = A[istart:iend, icol].astype(np.float64, copy=False)

        formatted = _fmt_values_e23(floats)
        buf_parts.append(formatted)
        buf_size += len(formatted) + 25

        if buf_size > flush_threshold:
            op4.write(''.join(buf_parts))
            buf_parts.clear()
            buf_size = 0

    # Terminator
    buf_parts.append('%8i%8i%8i\n' % (ncols + 1, 1, 1))
    buf_parts.append(' 1.0000000000000000E+00\n')
    op4.write(''.join(buf_parts))


def _write_dense_matrix_binary(op4: BinaryIO, name: str, matrix: np.ndarray, form: int=2,
                               precision: str='default', encoding: str='utf-8', endian: str='<') -> None:
    """Writes a dense matrix in binary OP4 format with proper Fortran record framing.

    Binary format per record: [4:reclen][reclen bytes of data][4:reclen]
    Header record (reclen=24): ncols(i4), nrows(i4), form(i4), type(i4), name(8s)
    Header record (reclen=48): ncols(u8), nrows(u8), form(u8), type(u8), name(16s)
    Column records (reclen=12+nwords*4): icol(i4), irow(i4), nwords(i4), values
    Terminator: column record with icol=ncols+1, irow=1, nwords=1, value=1.0

    big_mat is signaled by writing nrows as negative when nrows > 65535.
    48-byte headers are used when name > 8 chars or dimensions exceed int32.
    """
    A = matrix
    matrix_type, nwords_per_value = _get_type_nwv(A[0, 0], precision)
    (nrows, ncols) = A.shape

    dt_i4 = np.dtype(f'{endian}i4')
    i4_pack = Struct(endian + 'i')

    # Determine header size: 48-byte if name > 8 chars or dims exceed int32
    use_large_header = (len(name) > 8 or ncols > 2147483647 or nrows > 2147483647)

    if use_large_header:
        # 48-byte header uses uint64: reader detects big_mat via nrows > 65535 (no negation needed)
        name2 = '%-16s' % name
        assert len(name2) == 16, 'name=%r is too long; 16 characters max' % name
        name_bytes = name2.encode('ascii')
        dt_u8 = np.dtype(f'{endian}u8')
        hdr_reclen = i4_pack.pack(48)
        hdr_data = np.array([ncols, nrows, form, matrix_type], dtype=dt_u8).tobytes() + name_bytes
        op4.write(hdr_reclen + hdr_data + hdr_reclen)
    else:
        # 24-byte header uses int32: negate nrows to signal big_mat when nrows > 65535
        nrows_header = -nrows if nrows > 65535 else nrows
        name2 = '%-8s' % name
        assert len(name2) == 8, 'name=%r is too long; 8 characters max' % name
        name_bytes = name2.encode('ascii')
        hdr_reclen = i4_pack.pack(24)
        hdr_data = np.array([ncols, nrows_header, form, matrix_type], dtype=dt_i4).tobytes() + name_bytes
        op4.write(hdr_reclen + hdr_data + hdr_reclen)

    # OP4 nwords is in 4-byte word units: {type1: 1, type2: 2, type3: 2, type4: 4}
    if matrix_type == 1:
        dt_val = np.dtype(f'{endian}f4')
        nwords_per_value = 1
    elif matrix_type == 2:
        dt_val = np.dtype(f'{endian}f8')
        nwords_per_value = 2
    elif matrix_type == 3:
        dt_val = np.dtype(f'{endian}f4')
        nwords_per_value = 2
    else:  # matrix_type == 4
        dt_val = np.dtype(f'{endian}f8')
        nwords_per_value = 4
    is_complex = matrix_type in (3, 4)

    # Ensure column-major (Fortran) order for efficient column slicing
    if not A.flags['F_CONTIGUOUS']:
        A = np.asfortranarray(A)

    # Vectorized column start/end computation: find first/last nonzero row per column
    nonzero_mask = (A != 0)
    has_data = nonzero_mask.any(axis=0)  # bool array of length ncols
    # argmax on bool gives first True; flip for last True
    col_starts = np.argmax(nonzero_mask, axis=0)  # first nonzero row per col
    col_ends = nrows - 1 - np.argmax(nonzero_mask[::-1, :], axis=0)  # last nonzero row per col

    # Build output buffer
    chunks: list[bytes] = []
    buf_size = 0
    flush_threshold = 4 * 1024 * 1024  # 4MB

    for icol in range(ncols):
        if not has_data[icol]:
            continue

        istart = int(col_starts[icol])
        iend = int(col_ends[icol]) + 1
        nvalues = iend - istart
        nwords = nvalues * nwords_per_value

        if is_complex:
            segment = A[istart:iend, icol]
            interleaved = np.empty(nvalues * 2, dtype=dt_val)
            interleaved[0::2] = segment.real
            interleaved[1::2] = segment.imag
            val_bytes = interleaved.tobytes()
        else:
            val_bytes = A[istart:iend, icol].astype(dt_val, copy=False).tobytes()

        col_reclen = 12 + nwords * 4
        reclen_bytes = i4_pack.pack(col_reclen)
        col_hdr = np.array([icol + 1, istart + 1, nwords], dtype=dt_i4).tobytes()
        chunks.append(reclen_bytes + col_hdr + val_bytes + reclen_bytes)
        buf_size += len(chunks[-1])

        if buf_size > flush_threshold:
            op4.write(b''.join(chunks))
            chunks.clear()
            buf_size = 0

    # Terminator column: icol=ncols+1, irow=1, nwords depends on precision
    if matrix_type in (1, 3):
        term_val = np.array([1.0], dtype=np.dtype(f'{endian}f4')).tobytes()
        term_nwords = 1
    else:
        term_val = np.array([1.0], dtype=np.dtype(f'{endian}f8')).tobytes()
        term_nwords = 2
    term_reclen = 12 + term_nwords * 4
    reclen_bytes = i4_pack.pack(term_reclen)
    term_hdr = np.array([ncols + 1, 1, term_nwords], dtype=dt_i4).tobytes()
    chunks.append(reclen_bytes + term_hdr + term_val + reclen_bytes)
    op4.write(b''.join(chunks))


def _write_sparse_matrix_binary(op4: BinaryIO, name: str, A: coo_matrix, form: int=2,
                                precision: str='default', endian: str='<') -> None:
    """Writes a sparse matrix in binary OP4 format.

    Sparse column record layout:
      [4:col_reclen][icol(i4), irow=0(i4), nwords(i4), segment_data...][4:col_reclen]

    big_mat segments: [L+1(i4), irow(i4), values...]
      where L = nvalues * nwords_per_value
    small_mat segments: [IS(i4), values...]
      where IS = (L+1)*65536 + irow, L = nvalues * nwords_per_value
    """
    (nrows, ncols) = A.shape
    if A.nnz > 0:
        matrix_type, _ = _get_type_nwv(A.data[0], precision)
    else:
        if A.dtype in (float32, np.dtype('float32')):
            matrix_type = 1
        elif A.dtype in (complex64, np.dtype('complex64')):
            matrix_type = 3
        elif A.dtype in (complex128, np.dtype('complex128')):
            matrix_type = 4
        else:
            matrix_type = 2
    is_big_mat = (nrows > 65535)

    dt_i4 = np.dtype(f'{endian}i4')
    i4_pack = Struct(endian + 'i')

    # Header
    use_large_header = (len(name) > 8 or ncols > 2147483647 or nrows > 2147483647)
    if use_large_header:
        name2 = '%-16s' % name
        assert len(name2) == 16, 'name=%r is too long; 16 characters max' % name
        name_bytes = name2.encode('ascii')
        dt_u8 = np.dtype(f'{endian}u8')
        hdr_reclen = i4_pack.pack(48)
        hdr_data = np.array([ncols, nrows, form, matrix_type], dtype=dt_u8).tobytes() + name_bytes
    else:
        name2 = '%-8s' % name
        assert len(name2) == 8, 'name=%r is too long; 8 characters max' % name
        name_bytes = name2.encode('ascii')
        nrows_header = -nrows if is_big_mat else nrows
        hdr_reclen = i4_pack.pack(24)
        hdr_data = np.array([ncols, nrows_header, form, matrix_type], dtype=dt_i4).tobytes() + name_bytes
    op4.write(hdr_reclen + hdr_data + hdr_reclen)

    if A.nnz == 0:
        # Empty matrix: just write terminator
        if matrix_type in (1, 3):
            term_val = np.array([1.0], dtype=np.dtype(f'{endian}f4')).tobytes()
            term_nwords = 1
        else:
            term_val = np.array([1.0], dtype=np.dtype(f'{endian}f8')).tobytes()
            term_nwords = 2
        term_reclen = 12 + term_nwords * 4
        reclen_bytes = i4_pack.pack(term_reclen)
        term_hdr = np.array([ncols + 1, 1, term_nwords], dtype=dt_i4).tobytes()
        op4.write(reclen_bytes + term_hdr + term_val + reclen_bytes)
        return

    # Value dtype and nwords_per_value (in 4-byte word units)
    if matrix_type == 1:
        dt_val = np.dtype(f'{endian}f4')
        nwords_per_value = 1
        val_bytes_per_entry = 4
    elif matrix_type == 2:
        dt_val = np.dtype(f'{endian}f8')
        nwords_per_value = 2
        val_bytes_per_entry = 8
    elif matrix_type == 3:
        dt_val = np.dtype(f'{endian}f4')
        nwords_per_value = 2
        val_bytes_per_entry = 8
    else:  # matrix_type == 4
        dt_val = np.dtype(f'{endian}f8')
        nwords_per_value = 4
        val_bytes_per_entry = 16
    is_complex = matrix_type in (3, 4)

    # Convert to CSC for column-ordered access
    csc = A.tocsc()
    if not csc.has_sorted_indices:
        csc.sort_indices()
    indices = csc.indices
    indptr = csc.indptr

    # Pre-convert all values to target dtype bytes
    if is_complex:
        raw_data = csc.data
        interleaved = np.empty(len(raw_data) * 2, dtype=dt_val)
        interleaved[0::2] = raw_data.real
        interleaved[1::2] = raw_data.imag
        data_bytes_all = interleaved.tobytes()
    else:
        data_bytes_all = csc.data.astype(dt_val, copy=False).tobytes()

    # Vectorized segment detection across entire matrix
    # A segment starts at: beginning of each non-empty column, or within-column row discontinuity
    is_seg_start = np.zeros(len(indices), dtype=bool)
    nonempty_cols = np.where(np.diff(indptr) > 0)[0]
    is_seg_start[indptr[nonempty_cols]] = True

    if len(indices) > 1:
        all_diff = np.diff(indices)
        col_boundary = np.zeros(len(indices), dtype=bool)
        # Filter boundary indices to avoid out-of-bounds when many trailing columns are empty
        boundary_idx = indptr[1:-1]
        boundary_idx = boundary_idx[boundary_idx < len(indices)]
        col_boundary[boundary_idx] = True
        is_within_col = ~col_boundary[1:]
        intra_breaks = is_within_col & (all_diff != 1)
        is_seg_start[np.where(intra_breaks)[0] + 1] = True

    seg_start_pos = np.where(is_seg_start)[0]
    nseg = len(seg_start_pos)

    # Segment end positions
    seg_end_pos = np.empty(nseg, dtype=np.int64)
    seg_end_pos[:-1] = seg_start_pos[1:]
    seg_end_pos[-1] = len(indices)
    # Clip to column boundaries
    seg_cols = np.searchsorted(indptr[1:], seg_start_pos, side='right')
    seg_end_pos = np.minimum(seg_end_pos, indptr[seg_cols + 1])

    # Segment properties
    seg_nvalues = (seg_end_pos - seg_start_pos).astype(np.int32)
    seg_irow = indices[seg_start_pos]  # 0-based
    L_vals = seg_nvalues * nwords_per_value

    # Segment headers as bytes
    if is_big_mat:
        # big_mat: [L+1(i4), irow+1(i4)] per segment
        seg_hdrs = np.empty((nseg, 2), dtype=dt_i4)
        seg_hdrs[:, 0] = L_vals + 1
        seg_hdrs[:, 1] = seg_irow + 1
        seg_hdr_bytes_all = seg_hdrs.tobytes()
        seg_hdr_size = 8
        seg_words = (2 + L_vals).astype(np.int32)
    else:
        # small_mat: [IS(i4)] per segment
        IS_all = ((L_vals + 1) * 65536 + (seg_irow + 1)).astype(dt_i4)
        seg_hdr_bytes_all = IS_all.tobytes()
        seg_hdr_size = 4
        seg_words = (1 + L_vals).astype(np.int32)

    # Group segments by column and compute per-column total_nwords
    unique_cols_arr, col_seg_start_idx = np.unique(seg_cols, return_index=True)
    col_seg_end_idx = np.concatenate([col_seg_start_idx[1:], [nseg]])
    col_seg_counts = col_seg_end_idx - col_seg_start_idx
    col_total_nwords = np.add.reduceat(seg_words, col_seg_start_idx)

    # Pre-compute per-column record headers
    col_reclens = (12 + col_total_nwords * 4).astype(dt_i4)
    col_hdrs = np.zeros((len(unique_cols_arr), 3), dtype=dt_i4)
    col_hdrs[:, 0] = unique_cols_arr + 1  # 1-based icol
    # col_hdrs[:, 1] = 0  (irow=0 for sparse)
    col_hdrs[:, 2] = col_total_nwords

    # Serialize: pre-convert to byte arrays for slicing
    col_reclen_bytes = col_reclens.tobytes()
    col_hdr_bytes = col_hdrs.tobytes()

    # Build output with tight loop
    out = bytearray()
    seg_idx = 0
    for i in range(len(unique_cols_arr)):
        r_off = i * 4
        h_off = i * 12
        reclen_b = col_reclen_bytes[r_off:r_off + 4]
        out.extend(reclen_b)
        out.extend(col_hdr_bytes[h_off:h_off + 12])

        n_segs = int(col_seg_counts[i])
        for _ in range(n_segs):
            # Segment header
            sh_off = seg_idx * seg_hdr_size
            out.extend(seg_hdr_bytes_all[sh_off:sh_off + seg_hdr_size])
            # Values
            val_start = int(seg_start_pos[seg_idx]) * val_bytes_per_entry
            val_end = int(seg_end_pos[seg_idx]) * val_bytes_per_entry
            out.extend(data_bytes_all[val_start:val_end])
            seg_idx += 1

        out.extend(reclen_b)

    # Terminator
    if matrix_type in (1, 3):
        term_val = np.array([1.0], dtype=np.dtype(f'{endian}f4')).tobytes()
        term_nwords = 1
    else:
        term_val = np.array([1.0], dtype=np.dtype(f'{endian}f8')).tobytes()
        term_nwords = 2
    term_reclen = 12 + term_nwords * 4
    reclen_bytes = i4_pack.pack(term_reclen)
    term_hdr = np.array([ncols + 1, 1, term_nwords], dtype=dt_i4).tobytes()
    out.extend(reclen_bytes + term_hdr + term_val + reclen_bytes)

    op4.write(bytes(out))


def _write_sparse_matrix_ascii(op4: TextIO, name: str, A: coo_matrix,
                               form: int=2, is_big_mat: bool=False,
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

    cols: dict[int, list[int]] = {}
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

def get_big_mat_nrows(nrows: int) -> tuple[bool, int]:
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
        # nrows=0 indicates an empty/uninitialized matrix (e.g., from ASSIGN
        # with no OUTPUT4 data written, or a null matrix like BHH with no
        # structural damping)
        raise EmptyMatrixError('nrows=0')
    return is_big_mat, nrows


def _get_matrix_info(matrix_type: int,
                     log: SimpleLogger,
                     debug: bool=True) -> tuple[int, int, str, str]:
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
    if debug:
        log.info('matrix_type = %s' % matrix_type)
        log.info('  nwords_per_value = %s' % nwords_per_value)
        log.info('  nbytes_per_value = %s' % nbytes_per_value)
        log.info('  dtype = %s ' % dtype)
    return nwords_per_value, nbytes_per_value, data_format, dtype

def get_dtype(matrix_type: int, precision: str='default') -> str:
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


def _get_type_nwv(A: np.ndarray, precision: str='default') -> tuple[int, int]:
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

def _prepare_name_order(matrices: dict[str, Matrix], name_order) -> list[str]:
    """Helper method for OP4 writing"""
    if name_order is None:
        name_order = sorted(matrices.keys())
    elif isinstance(name_order, str):
        name_order = [name_order]
    elif isinstance(name_order, bytes):
        name_order = [name_order]
    return name_order

def _write_form_matrix_helper(matrices: dict[str, Matrix],
                              name: str) -> tuple[int, np.ndarray | coo_matrix]:
    """Helper method for OP4 writing"""
    try:
        mat_form = matrices[name]
    except KeyError:
        raise KeyError(f'key={name!r} is an invalid matrix; keys={matrices.keys()}')

    if isinstance(mat_form, Matrix):
        form = mat_form.form
        matrix = mat_form.data
    else:
        (form, matrix) = mat_form

    if not form in {1, 2, 3, 6, 8, 9}:
        raise ValueError(f'form={form!r} and must be in [1, 2, 3, 6, 8, 9]')
    return form, matrix

def _write_data(f: TextIO, data: bytes, endian: str, types: str='ifs'):
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

    assert endian is not None, endian

    data4 = data[:nints * 4]
    if 's' in types:
        strings = unpack('%s%is' % (endian, n), data[:n])
        f.write("strings(s) = %s\n" % str(strings))
    if 'i' in types:
        ints = unpack('%s%ii' % (endian, nints), data4)
        f.write("ints(i)    = %s\n" % str(ints))
    if 'f' in types:
        floats = unpack('%s%if' % (endian, nints), data4)
        f.write("floats(f)  = %s\n" % str(floats))
    if 'd' in types:
        doubles = unpack('%s%id' % (endian, ndoubles), data[:ndoubles*8])
        f.write("doubles(d)  = %s\n" % str(doubles))

    if 'l' in types:
        longs = unpack('%s%il' % (endian, nints), data4)
        f.write("long(l)  = %s\n" % str(longs))
    if 'I' in types:
        ints2 = unpack('%s%iI' % (endian, nints), data4)
        f.write("unsigned int(I) = %s\n" % str(ints2))
    if 'L' in types:
        longs2 = unpack('%s%iL' % (endian, nints), data4)
        f.write("unsigned long(L) = %s\n" % str(longs2))
    if 'q' in types:
        longs = unpack('%s%iq' % (endian, ndoubles), data[:ndoubles*8])
        f.write("long long(q) = %s\n" % str(longs))
    if 'Q' in types:
        longs = unpack('%s%iq' % (endian, ndoubles), data[:ndoubles*8])
        f.write("unsigned long long(Q) = %s\n" % str(longs))
    return strings, ints, floats

def is_saved_matrix(name: str, matrix_names: Optional[list[str]]) -> bool:
    #name = name.decode('ascii')
    if name is not None:
        if matrix_names is None or name in matrix_names:
            return True
    return False

def get_dtype_dtval(matrix_type: int,
                    precision: str) -> tuple[Any, Any]:
    matrix_type = MATRIX_TYPE_MAP[(matrix_type, precision)]
    if matrix_type in (1, 3):
        dtype = float32 if not is_complex else complex64
        dt_val = float32
    else:
        dtype = float64 if not is_complex else complex128
        dt_val = float64
    return dtype, dt_val
            
def _read_op4_ascii(op4_filename: PathLike,
                    precision: str='default',
                    matrix_names: Optional[list[str]]=None) -> dict[str, Matrix]:
    """
    Fast ASCII OP4 reader. Minimizes per-line Python
    overhead by deferring string-to-float conversion
    to numpy (C-level batch).

    Parameters
    ----------
    precision : str; {'default', 'single', 'double'}
        specifies if the matrices are in single or double precsion
        which means the format will be whatever the file is in
    """
    assert precision in {'default', 'single', 'double'}, f'precision={precision!r}'
    with open(op4_filename, 'r') as f:
        lines = f.readlines()

    matrices: dict[str, Matrix] = {}

    i = 0
    nlines = len(lines)
    while i < nlines:
        line = lines[i]
        i += 1
        if len(line) < 32:
            if not line.strip():
                break
            continue

        # Matrix header: ncols, nrows, form, type, name, format_spec
        # Non-header lines (terminator data, trailing values) are skipped here.
        # TODO: replace try/except flow control with explicit header validation
        parts = line[0:32].split()
        if len(parts) != 4:
            continue
        try:
            ncols = int(parts[0])
            nrows_raw = int(parts[1])
            form = int(parts[2])
            matrix_type = int(parts[3])
        except ValueError:
            log.error(str(parts))
            continue
        name = line[32:40].strip()

        try:
            is_big_mat, nrows = get_big_mat_nrows(nrows_raw)
        except EmptyMatrixError:
            log.error(f'nrows_raw={nrows_raw}')
            break

        # Parse format spec: e.g. "1P,3E23.16" or "1P,5E16.9"
        size_str = line[40:].strip()
        line_size = int(size_str.split(',')[1].split('E')[1].split('.')[0])

        is_complex = matrix_type in (3, 4)
        dtype, dt_val = get_dtype_dtval(matrix_type, precision)

        # First column header
        line = lines[i]
        i += 1
        parts = line.split()
        irow_first = int(parts[1])
        is_sparse = (irow_first == 0)

        # Skip this matrix if not requested
        if matrix_names is not None and name not in matrix_names:
            icol = int(parts[0])
            while icol <= ncols and i < nlines:
                line = lines[i]
                i += 1
                if 'E' in line:
                    continue
                sline = line.split()
                if len(sline) == 3:
                    icol = int(sline[0])
            continue

        rows_list: list[int] = []
        cols_list: list[int] = []
        vals_strs: list[str] = []

        icol = int(parts[0])
        irow = int(parts[1])
        nwords = int(parts[2])

        while icol <= ncols:
            if is_sparse:
                seg_irow = 0
                seg_strs: list[str] = []
                while i < nlines:
                    line = lines[i]
                    if 'E' in line:
                        # Data line — fast path (most common)
                        nw = line.count('E')
                        seg_strs.extend(line[j*line_size:(j+1)*line_size] for j in range(nw))
                        i += 1
                    else:
                        sline = line.split()
                        nsline = len(sline)
                        if nsline == 1:
                            # Flush previous segment
                            if seg_strs and seg_irow > 0:
                                if is_complex:
                                    nv = len(seg_strs) // 2
                                else:
                                    nv = len(seg_strs)
                                rows_list.extend(range(seg_irow - 1, seg_irow - 1 + nv))
                                cols_list.extend([icol - 1] * nv)
                                vals_strs.extend(seg_strs)
                                seg_strs = []
                            IS = int(sline[0])
                            L = IS // 65536 - 1
                            seg_irow = IS - 65536 * (L + 1)
                            i += 1
                        elif nsline == 2:
                            if seg_strs and seg_irow > 0:
                                if is_complex:
                                    nv = len(seg_strs) // 2
                                else:
                                    nv = len(seg_strs)
                                rows_list.extend(range(seg_irow - 1, seg_irow - 1 + nv))
                                cols_list.extend([icol - 1] * nv)
                                vals_strs.extend(seg_strs)
                                seg_strs = []
                            L = int(sline[0]) - 1
                            seg_irow = int(sline[1])
                            i += 1
                        elif nsline >= 3:
                            break
                        else:
                            i += 1

                # Flush final segment
                if seg_strs and seg_irow > 0:
                    if is_complex:
                        nv = len(seg_strs) // 2
                    else:
                        nv = len(seg_strs)
                    rows_list.extend(range(seg_irow - 1, seg_irow - 1 + nv))
                    cols_list.extend([icol - 1] * nv)
                    vals_strs.extend(seg_strs)

            else:
                # Dense column — tight loop on data lines
                col_start = len(vals_strs)
                while nwords > 0 and i < nlines:
                    line = lines[i]
                    if 'E' not in line:
                        break
                    nw = line.count('E')
                    vals_strs.extend(line[j*line_size:(j+1)*line_size] for j in range(nw))
                    nwords -= nw
                    i += 1

                nfloats = len(vals_strs) - col_start
                if nfloats > 0:
                    if is_complex:
                        nv = nfloats // 2
                    else:
                        nv = nfloats
                    rows_list.extend(range(irow - 1, irow - 1 + nv))
                    cols_list.extend([icol - 1] * nv)

            # Read next column header (skip data/marker lines)
            while i < nlines:
                line = lines[i]
                i += 1
                # Data lines contain 'E' and '.'
                if 'E' in line and '.' in line:
                    continue
                sline = line.split()
                if len(sline) >= 3 and '.' not in line:
                    icol = int(sline[0])
                    irow = int(sline[1])
                    nwords = int(sline[2])
                    break
                # 1 or 2 integer lines (segment markers) - skip
            else:
                break

        # Skip terminator column data lines
        while i < nlines:
            line = lines[i].rstrip()
            if not line:
                i += 1
                break
            # Matrix header has 4 integers in first 32 chars
            try:
                hdr_fields = line[0:32].split()
                if len(hdr_fields) == 4:
                    int(hdr_fields[0]); int(hdr_fields[1])
                    int(hdr_fields[2]); int(hdr_fields[3])
                    break
            except (ValueError, IndexError):
                pass
            i += 1

        if rows_list:
            all_rows = np.array(rows_list, dtype=np.int32)
            all_cols = np.array(cols_list, dtype=np.int32)
            if is_complex:
                raw = np.array(vals_strs, dtype=dt_val)
                all_vals = (raw[0::2] + 1j * raw[1::2]).astype(dtype)
            else:
                all_vals = np.array(vals_strs, dtype=dtype)
            data_mat = coo_matrix((all_vals, (all_rows, all_cols)),
                                  shape=(nrows, ncols), dtype=dtype)
        else:
            data_mat = coo_matrix((nrows, ncols), dtype=dtype)

        amat = Matrix(name, form, data=data_mat)
        _save_matrix(matrices, name, amat)

    return matrices


MATRIX_TYPE_MAP = {
    (1, 'default'): 1,
    (2, 'default'): 2,
    (3, 'default'): 3,
    (4, 'default'): 4,

    (1, 'single'): 1,
    (2, 'single'): 1,
    (3, 'single'): 3,
    (4, 'single'): 3,

    (1, 'double'): 2,
    (2, 'double'): 2,
    (3, 'double'): 4,
    (4, 'double'): 4,
}

def read_op4_fast(op4_filename: PathLike | None,
                  precision: str='default',
                  matrix_names: Optional[list[str]]=None) -> dict[str, Matrix]:
    """
    Fast OP4 reader using bulk buffer reads and numpy parsing.
    Supports both binary and ASCII formats.

    Performance vs read_op4 (typical speedups):
      - Dense binary full columns (200MB): ~7x faster (returns ndarray directly)
      - Dense binary partial columns (symmetric): ~3-4x faster
      - Sparse binary: ~2-3x faster

    Key optimizations:
      - Reads entire file into memory buffer once
      - Scans Fortran record positions with int.from_bytes on memoryview
        (avoids per-record np.frombuffer overhead)
      - For dense full-column matrices: reads directly into F-contiguous ndarray
        (skips COO intermediate and concatenation)
      - Vectorized row/col array construction via np.repeat + cumsum trick
        (eliminates per-column np.arange/np.full calls)
      - Memoryview byte-copy for value extraction into pre-allocated buffer
        (single np.frombuffer at end vs per-column)

    Parameters
    ----------
    op4_filename : PathLike
        Path to the OP4 file (binary or ASCII).
    matrix_names : list[str] or None
        If specified, only read these matrices (others are skipped).
    precision : str; {'default', 'single', 'double'}
        specifies if the matrices are in single or double precsion
        which means the format will be whatever the file is in

    Returns
    -------
    dict[str, Matrix]
        Dictionary mapping matrix names to Matrix objects.
        Dense full-column matrices: Matrix.data is np.ndarray (F-order).
        Partial-column or sparse matrices: Matrix.data is scipy.sparse.coo_matrix.
    """
    if precision not in {'default', 'single', 'double'}:
        msg = f"precision={precision!r} and must be 'single', 'double', or 'default'"
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
        raise IOError('cannot find op4_filename=%r' % str(op4_filename))

    if not file_is_binary(op4_filename):
        return _read_op4_ascii(
            op4_filename, precision=precision,
            matrix_names=matrix_names)

    with open(op4_filename, 'rb') as f:
        buf = f.read()
    if len(buf) == 0:
        return {}

    # Determine endianness from first record length marker
    rl_le = int.from_bytes(buf[0:4], 'little', signed=True)
    if rl_le in (24, 48):
        endian = '<'
    else:
        rl_be = int.from_bytes(buf[0:4], 'big', signed=True)
        if rl_be in (24, 48):
            endian = '>'
        else:
            raise RuntimeError(
                f'unrecognized binary OP4 header record length '
                f'(le={rl_le}, be={rl_be}) in {op4_filename}')

    dt_i4 = np.dtype(f'{endian}i4')
    dt_f4 = np.dtype(f'{endian}f4')
    dt_f8 = np.dtype(f'{endian}f8')

    matrices: dict[str, Matrix] = {}
    # Header: [4:reclen][reclen bytes: ncols,nrows,form,type,name][4:reclen]
    pos = 0
    while pos < len(buf):
        if pos + 4 > len(buf):
            break
        hdr_reclen = np.frombuffer(buf, dt_i4, count=1, offset=pos)[0]
        if hdr_reclen not in (24, 48):
            break
        pos += 4  # past leading Fortran marker

        if hdr_reclen == 24:
            hdr = np.frombuffer(buf, dt_i4, count=4, offset=pos)
            ncols, nrows, form, matrix_type = int(hdr[0]), int(hdr[1]), int(hdr[2]), int(hdr[3])
            name = buf[pos+16:pos+24].decode('ascii').strip()
            pos += 24 + 4  # header data + trailing Fortran marker
        else:  # hdr_reclen == 48
            dt_u8 = np.dtype(f'{endian}u8')
            hdr = np.frombuffer(buf, dt_u8, count=4, offset=pos)
            ncols, nrows, form, matrix_type = int(hdr[0]), int(hdr[1]), int(hdr[2]), int(hdr[3])
            name = buf[pos+32:pos+48].decode('ascii').strip()
            pos += 48 + 4  # header data + trailing Fortran marker

        # big_mat determination
        if nrows < 0:
            is_big_mat = True
            nrows = -nrows
        elif nrows > 65535:
            is_big_mat = True
        else:
            is_big_mat = False

        if matrix_names is not None and name not in matrix_names:
            # Skip this matrix by scanning past its column records
            while pos < len(buf):
                col_reclen = int(np.frombuffer(buf, dt_i4, count=1, offset=pos)[0])
                pos += 4
                rec_hdr = np.frombuffer(buf, dt_i4, count=3, offset=pos)
                icol = int(rec_hdr[0])
                pos += col_reclen + 4  # record data + trailing marker
                if icol == ncols + 1:
                    break
            continue

        if matrix_type == 1:
            nbytes_per_value = 4
            dt_val = dt_f4
            dtype = float32
            nwords_per_value = 1
            is_complex = False
        elif matrix_type == 2:
            nbytes_per_value = 8
            dt_val = dt_f8
            dtype = float64
            nwords_per_value = 2
            is_complex = False
        elif matrix_type == 3:
            nbytes_per_value = 8
            dt_val = dt_f4
            dtype = complex64
            nwords_per_value = 2
            is_complex = True
        elif matrix_type == 4:
            nbytes_per_value = 16
            dt_val = dt_f8
            dtype = complex128
            nwords_per_value = 4
            is_complex = True
        else:
            raise RuntimeError(
                f'unsupported matrix_type={matrix_type} for matrix '
                f'{name!r} in {op4_filename}')

        # Column format (Fortran unformatted records):
        # Each column is one Fortran record:
        #   [4: col_reclen] [col_reclen bytes: icol, irow, nwords, data...] [4: col_reclen]
        # Followed by terminator column with icol = ncols+1.
        # The 'a' field in the original code = col_reclen.
        # icol, irow, nwords are 4-byte ints at the start of the record data.
        # Data follows: nwords*4 bytes. col_reclen = 12 + nwords*4.

        # Peek at irow of first column to determine sparse vs dense
        # pos is at the leading Fortran marker of first column record
        # irow is at pos + 4 (leading marker) + 4 (icol) = pos + 8
        peek_irow = int(np.frombuffer(buf, dt_i4, count=1, offset=pos + 8)[0])
        is_sparse = (peek_irow == 0)

        if is_sparse:
            # Optimized sparse reader: scan segments collecting metadata,
            # then batch-construct row/col/val arrays
            mv = memoryview(buf)
            from_bytes = int.from_bytes
            byteorder = 'little' if endian == '<' else 'big'

            # Collect segment metadata: (icol_0based, irow_0based, nvalues, val_offset)
            seg_meta = []  # list of tuples
            total_values = 0

            while True:
                if pos + 4 > len(buf):
                    break
                col_reclen = from_bytes(mv[pos:pos+4], byteorder, signed=True)
                pos += 4
                icol = from_bytes(mv[pos:pos+4], byteorder, signed=True)
                if icol == ncols + 1:
                    pos += col_reclen + 4
                    break
                pos += 12  # skip icol, irow=0, nwords

                data_bytes = col_reclen - 12
                col_end = pos + data_bytes
                icol_0 = icol - 1

                if is_big_mat:
                    while pos < col_end:
                        if pos + 8 > col_end:
                            break
                        L = from_bytes(mv[pos:pos+4], byteorder, signed=True) - 1
                        irow = from_bytes(mv[pos+4:pos+8], byteorder, signed=True)
                        pos += 8
                        if L <= 0 or irow <= 0:
                            break
                        nvalues = L // nwords_per_value
                        seg_meta.append((icol_0, irow - 1, nvalues, pos))
                        total_values += nvalues
                        pos += nvalues * nbytes_per_value
                else:
                    while pos < col_end:
                        if pos + 4 > col_end:
                            break
                        IS = from_bytes(mv[pos:pos+4], byteorder, signed=True)
                        pos += 4
                        L = IS // 65536 - 1
                        irow = IS - 65536 * (L + 1)
                        if L <= 0 or irow <= 0:
                            break
                        nvalues = L // nwords_per_value
                        seg_meta.append((icol_0, irow - 1, nvalues, pos))
                        total_values += nvalues
                        pos += nvalues * nbytes_per_value

                pos = col_end + 4  # trailing Fortran marker

            if seg_meta:
                n_segs = len(seg_meta)
                # Build col array: repeat each icol by nvalues
                seg_cols = np.empty(n_segs, dtype=np.int32)
                seg_irows = np.empty(n_segs, dtype=np.int32)
                seg_nvals = np.empty(n_segs, dtype=np.int32)

                for si in range(n_segs):
                    seg_cols[si] = seg_meta[si][0]
                    seg_irows[si] = seg_meta[si][1]
                    seg_nvals[si] = seg_meta[si][2]

                all_cols = np.repeat(seg_cols, seg_nvals)

                # Build row array: for each segment, arange(irow, irow+nvalues)
                # Use cumsum trick: fill with 1s, then at segment boundaries reset
                all_rows = np.ones(total_values, dtype=np.int32)
                seg_starts = np.empty(n_segs, dtype=np.int64)
                seg_starts[0] = 0
                if n_segs > 1:
                    np.cumsum(seg_nvals[:-1], out=seg_starts[1:])
                all_rows[seg_starts] = seg_irows
                # Subtract previous end value at each boundary (except first)
                if n_segs > 1:
                    prev_ends = seg_irows[:-1] + seg_nvals[:-1] - 1
                    all_rows[seg_starts[1:]] -= prev_ends
                np.cumsum(all_rows, out=all_rows)

                # Build values array via memoryview byte copy
                if is_complex:
                    val_buf = bytearray(total_values * nbytes_per_value)
                    bi = 0
                    for si in range(n_segs):
                        nv = seg_meta[si][2]
                        off = seg_meta[si][3]
                        nb = nv * nbytes_per_value
                        val_buf[bi:bi+nb] = mv[off:off+nb]
                        bi += nb
                    raw = np.frombuffer(val_buf, dtype=dt_val)
                    all_vals = (raw[0::2] + 1j * raw[1::2]).astype(dtype)
                else:
                    val_buf = bytearray(total_values * nbytes_per_value)
                    bi = 0
                    for si in range(n_segs):
                        nv = seg_meta[si][2]
                        off = seg_meta[si][3]
                        nb = nv * nbytes_per_value
                        val_buf[bi:bi+nb] = mv[off:off+nb]
                        bi += nb
                    all_vals = np.frombuffer(val_buf, dtype=dt_val)

                data_mat = coo_matrix((all_vals, (all_rows, all_cols)),
                                      shape=(nrows, ncols), dtype=dtype)
            else:
                data_mat = coo_matrix((nrows, ncols), dtype=dtype)
        else:
            # Dense column format: optimized two-pass approach
            # Pass 1: scan record positions using int.from_bytes on memoryview
            #   (avoids per-column np.frombuffer overhead for 12-byte headers)
            # Pass 2: extract values based on column fullness:
            #   - All columns full (irow=1, nwords=nrows*nwpv): read into F-order ndarray
            #   - Partial columns (symmetric, banded): vectorized COO construction
            #     using np.repeat for cols and cumsum trick for rows
            mv = memoryview(buf)
            rec_positions = []  # (data_start, icol, irow, nwords) for each column record
            scan_pos = pos
            buf_len = len(buf)
            from_bytes = int.from_bytes
            byteorder = 'little' if endian == '<' else 'big'
            while scan_pos + 4 <= buf_len:
                col_reclen = from_bytes(mv[scan_pos:scan_pos+4], byteorder, signed=True)
                data_start = scan_pos + 4
                icol = from_bytes(mv[data_start:data_start+4], byteorder, signed=True)
                if icol == ncols + 1:
                    pos = data_start + col_reclen + 4
                    break
                irow = from_bytes(mv[data_start+4:data_start+8], byteorder, signed=True)
                nwords = from_bytes(mv[data_start+8:data_start+12], byteorder, signed=True)
                rec_positions.append((data_start + 12, icol, irow, nwords))
                scan_pos = data_start + col_reclen + 4  # past data + trailing marker
            else:
                pos = scan_pos

            # Pass 2: extract values — check if all columns are full (common case)
            if rec_positions:
                n_recs = len(rec_positions)
                all_full = True
                expected_nwords = nrows * nwords_per_value
                for _, icol_r, irow_r, nwords_r in rec_positions:
                    if irow_r != 1 or nwords_r != expected_nwords:
                        all_full = False
                        break

                if all_full and not is_complex:
                    # Fast path: all columns are full — read directly into ndarray
                    # TODO: could eliminate per-column loop if records are contiguous
                    #   (reshape buf slice directly), but Fortran framing bytes prevent this
                    result = np.empty((nrows, n_recs), dtype=dtype, order='F')
                    for col_idx in range(n_recs):
                        d_start = rec_positions[col_idx][0]
                        result[:, col_idx] = np.frombuffer(
                            buf, dt_val, count=nrows, offset=d_start)
                    data_mat = result
                elif all_full and is_complex:
                    # Fast path for complex: all columns full
                    result = np.empty((nrows, n_recs), dtype=dtype, order='F')
                    for col_idx in range(n_recs):
                        d_start = rec_positions[col_idx][0]
                        raw = np.frombuffer(buf, dt_val, count=nrows * 2, offset=d_start)
                        result[:, col_idx] = raw[0::2] + 1j * raw[1::2]
                    data_mat = result
                else:
                    # General path: partial columns — vectorized row/col construction
                    rec_icols = np.empty(n_recs, dtype=np.int32)
                    rec_irows = np.empty(n_recs, dtype=np.int32)
                    rec_nvals = np.empty(n_recs, dtype=np.int32)
                    for ri in range(n_recs):
                        rec_icols[ri] = rec_positions[ri][1] - 1
                        rec_irows[ri] = rec_positions[ri][2] - 1
                        rec_nvals[ri] = rec_positions[ri][3] // nwords_per_value

                    total_values = int(rec_nvals.sum())
                    all_cols = np.repeat(rec_icols, rec_nvals)

                    # Build row indices with cumsum trick
                    all_rows = np.ones(total_values, dtype=np.int32)
                    seg_starts = np.empty(n_recs, dtype=np.int64)
                    seg_starts[0] = 0
                    if n_recs > 1:
                        np.cumsum(rec_nvals[:-1], out=seg_starts[1:])
                    all_rows[seg_starts] = rec_irows
                    if n_recs > 1:
                        prev_ends = rec_irows[:-1] + rec_nvals[:-1] - 1
                        all_rows[seg_starts[1:]] -= prev_ends
                    np.cumsum(all_rows, out=all_rows)

                    # Extract values via memoryview byte copy
                    if is_complex:
                        val_buf = bytearray(total_values * nbytes_per_value)
                        bi = 0
                        for ri in range(n_recs):
                            nv = int(rec_nvals[ri])
                            d_start = rec_positions[ri][0]
                            nb = nv * nbytes_per_value
                            val_buf[bi:bi+nb] = mv[d_start:d_start+nb]
                            bi += nb
                        raw = np.frombuffer(val_buf, dtype=dt_val)
                        all_vals = (raw[0::2] + 1j * raw[1::2]).astype(dtype)
                    else:
                        val_buf = bytearray(total_values * nbytes_per_value)
                        bi = 0
                        for ri in range(n_recs):
                            nv = int(rec_nvals[ri])
                            d_start = rec_positions[ri][0]
                            nb = nv * nbytes_per_value
                            val_buf[bi:bi+nb] = mv[d_start:d_start+nb]
                            bi += nb
                        all_vals = np.frombuffer(val_buf, dtype=dt_val)

                    data_mat = coo_matrix((all_vals, (all_rows, all_cols)),
                                          shape=(nrows, ncols), dtype=dtype)
            else:
                data_mat = coo_matrix((nrows, ncols), dtype=dtype)

        matrix_type2 = MATRIX_TYPE_MAP[(matrix_type, precision)]
        if matrix_type != matrix_type2:
            if matrix_type2 == 1:
                dtypei = float32
            elif matrix_type2 == 2:
                dtypei = float64
            elif matrix_type2 == 3:
                dtypei = complex64
            elif matrix_type2 == 4:
               dtypei = complex128
            else:  # pragma: no cover
                raise RuntimeError(matrix_type)
            data_mat.astype(dtypei)

        amat = Matrix(name, form, data=data_mat)
        _save_matrix(matrices, name, amat)
    return matrices


def read_op4(op4_filename: Optional[PathLike]=None,
             matrix_names: Optional[list[str]]=None,
             precision: str='default',
             debug: bool=False, log=None) -> dict[str, Matrix]:
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
       >>> MatrixA = matrices['A']
       >>> MatrixB = matrices['B']
       >>> MatrixC = matrices['C']

       # or to reduce memory usage
       >>> matrices = op4.read_op4(op4_filename, matrix_names=['A', 'B'])
       >>> MatrixA = matrices['A']
       >>> MatrixB = matrices['B']

       # or because you only want A
       >>> matrices = op4.read_op4(op4_filename, matrix_names='A')
       >>> MatrixA = matrices['A']

       # get all the matrices, but select the file using a file dialog
       >>> matrices = op4.read_op4()
       >>>

    Parameters
    ----------
    op4_filename : str / None
        an OP4 filename.  Type=STRING.
    matrix_names : list[str], str / None
        matrix name(s) (None -> all)
    precision : str; {'default', 'single', 'double'}
        specifies if the matrices are in single or double precsion
        which means the format will be whatever the file is in

    Returns
    -------
    matricies : dict[str] = (int, Matrix)
        dictionary of matrices where the key is the name and the value is a matrix.
        To get the form: matrix.form
        To get the data: matrix.data

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
    matrices = op4.read_op4(
        op4_filename, matrix_names, precision)
    return matrices

def write_op4(op4_filename: Optional[PathLike],
              matrices: dict[str, Matrix],
              name_order=None,
              precision: str='default', is_binary: bool=True,
              log: Optional[SimpleLogger]=None, debug: bool=False) -> None:
    op4 = OP4(log=log, debug=debug)
    op4.write_op4(op4_filename, matrices,
                  name_order=name_order, precision=precision,
                  is_binary=is_binary)
