from __future__ import print_function
from six import string_types, iteritems, PY2, PY3
from six.moves import range
import os
import io
from struct import pack, unpack, Struct

from numpy import array, zeros, float32, float64, complex64, complex128, ndarray
from scipy.sparse import coo_matrix

from pyNastran.utils import is_binary_file as file_is_binary
from pyNastran.utils.mathematics import print_matrix, print_annotated_matrix


class OP4(object):
    """
     todo:: add endian checking
     todo:: test on big matrices
     todo:: finish write_op4
    """
    def __init__(self, log=None):
        self.n = 0
        self.endian = ''

    def read_op4(self, op4_filename=None, matrix_names=None, precision='default'):
        """
        Reads a NASTRAN OUTPUT4 file, regular or sparse, and stores the
        matrices as the output arguments of the function.  The number of
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
          >>> matrices = op4.read_op4(op4_filename, matrix_names=['A','B'])
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

              ==== ===============
              Form Definition
              ==== ===============
              1.   Square
              2.   Rectangular
              3.   Diagonal
              6.   Symmetric
              8.   Id entity
              9.   Pseudoidentity
              ==== ===============

            Matrix:
              Dense Type:  NUMPY.NDARRAY
              Sparse Type: SCIPY.SPARSE.COO_MATRIX

        .. note:: based off the MATLAB code SAVEOP4 developed by ATA-E and
                  later UCSD.
        .. note:: it's strongly recommended that you convert sparse matrices to
                  another format before doing math on them.  This is standard
                  with sparse matrices.

        .. warning:: sparse binary is buggy right now
        """
        if not precision in ('default', 'single', 'double'):
            raise ValueError("precision=%r and must be 'single', 'double', or 'default'" % precision)

        if op4_filename is None:
            from pyNastran.utils.gui_io import load_file_dialog
            wildcard_wx = "Nastran OP4 (*.op4)|*.op4|" \
                "All files (*.*)|*.*"
            wildcard_qt = "Nastran OP4 (*.op4);;All files (*)"
            title = 'Please select a OP4 to load'
            op4_filename, wildcard_level = load_file_dialog(title, wildcard_wx, wildcard_qt)
            assert op4_filename is not None, op4_filename

        if not os.path.exists(op4_filename):
            raise IOError('cannot find op4_filename=%r' % op4_filename)

        if isinstance(matrix_names, string_types):
            matrix_names = [matrix_names]
        #assert isinstance(matrix_names, list), 'type(matrix_names)=%s' % type(matrix_names)

        if file_is_binary(op4_filename):
            return self.read_op4_binary(op4_filename, matrix_names, precision)
        else:
            return self.read_op4_ascii(op4_filename, matrix_names, precision)

#--------------------------------------------------------------------------
    def read_op4_ascii(self, op4_filename, matrix_names=None, precision='default'):
        """matrix_names must be a list or None, but basically the same"""
        with open(op4_filename, 'r') as f:
            matrices = {}
            name = 'dummyName'
            while name is not None:
                (name, form, matrix) = self._read_matrix_ascii(f, matrix_names, precision)
                if name is not None:
                    if matrix_names is None or name in matrix_names:  # save the matrix
                        matrices[name] = (form, matrix)
        return matrices

    def _read_matrix_ascii(self, f, matrix_names=None, precision='default'):
        """reads a matrix"""
        line = f.readline().rstrip()
        if line == '':
            f.close()
            return None, None, None
        ncols, nrows, form, Type = line[0:32].split()
        nrows = int(nrows)

        if nrows < 0:  # if less than 0, big
            is_big_mat = True
        elif nrows > 0:
            is_big_mat = False
        else:
            raise RuntimeError('unknown BIGMAT.  nRows=%s' % nrows)

        nrows = abs(nrows)
        ncols = int(ncols)
        form = int(form)
        Type = int(Type)
        dtype = get_dtype(Type, precision)

        name = line[32:40].strip()
        size = line[40:].strip()
        line_size = size.split(',')[1].split('E')[1].split('.')[0]  # 3E23.16 to 23
        line_size = int(line_size)

        line = f.readline().rstrip()
        (icol, irow, nwords) = line.split()

        is_sparse = False
        if irow == '0':
            is_sparse = True

        if Type in [1, 2]:  # real
            A = self._read_real_ascii(f, nrows, ncols, line_size, line, dtype, is_sparse, is_big_mat)
        elif Type in [3, 4]:  # complex
            A = self._read_complex_ascii(f, nrows, ncols, line_size, line, dtype, is_sparse, is_big_mat)
        else:
            raise RuntimeError('invalid matrix type.  Type=%s' % Type)

        if not(matrix_names is None or name in matrix_names):  # kill the matrix
            A = None
        #print("form=%s name=%s A=\n%s" % (form, name, str(A)))
        return (name, form, A)

    def _read_real_sparse_ascii(self, f, nrows, ncols, line_size, line, dtype, is_big_mat):
        rows = []
        cols = []
        entries = []
        nloops = 0
        was_broken = False
        while 1:
            if nloops > 0 and not was_broken:
                line = f.readline().rstrip()
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
            while (len(sline) == 1 or len(sline) == 2) and 'E' not in line or run_loop:  # next sparse entry
                irow = self._get_irow_sparse_ascii(f, line, sline, nwords, irow,
                                                   is_big_mat)

                run_loop = False
                i = 0
                iword = 0
                is_done_reading_row = False
                while nwords:
                    n = 0
                    line = f.readline().rstrip()
                    nwords_in_line = line.count('E')
                    if nwords_in_line == 0:
                        was_broken = True
                        break

                    for i in range(nwords_in_line):
                        word = line[n:n + line_size]
                        rows.append(irow)
                        cols.append(icol)
                        entries.append(word)
                        n += line_size
                        irow += 1
                    iword += nwords_in_line
                    nwords -= nwords_in_line
                sline = line.strip().split()
                nloops += 1

        f.readline()

        #if rows == []:  # NULL matrix
            #raise NotImplementedError()

        rows = array(rows, dtype='int32') - 1
        cols = array(cols, dtype='int32') - 1
        A = coo_matrix((entries, (rows, cols)), shape=(nrows, ncols), dtype=dtype)
        #print("type = %s %s" % (type(A),type(A.todense())))
        #A = A.todense()
        return A

    def _read_real_dense_ascii(self, f, nrows, ncols, line_size, line, dtype, is_big_mat):
        A = zeros((nrows, ncols), dtype=dtype)  # Initialize a real matrix
        nloops = 0
        was_broken = False
        while 1:
            if nloops > 0 and not was_broken:
                line = f.readline().rstrip()
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
            while (len(sline) == 1 or len(sline) == 2) and 'E' not in line or run_loop:  # next dense entry
                run_loop = False
                i = 0
                iword = 0
                is_done_reading_row = False
                while nwords:
                    n = 0
                    line = f.readline().rstrip()
                    nwords_in_line = line.count('E')
                    if nwords_in_line == 0:
                        was_broken = True
                        break

                    for i in range(nwords_in_line):
                        word = line[n:n + line_size]
                        A[irow - 1, icol - 1] = word
                        n += line_size
                        irow += 1
                    iword += nwords_in_line
                    nwords -= nwords_in_line
                sline = line.strip().split()
                nloops += 1
        f.readline()
        return A

    def _read_real_ascii(self, f, nrows, ncols, line_size, line, dtype, is_sparse, is_big_mat):
        if is_sparse:
            A = self._read_real_sparse_ascii(f, nrows, ncols, line_size, line, dtype, is_big_mat)
        else:
            A = self._read_real_dense_ascii(f, nrows, ncols, line_size, line, dtype, is_big_mat)
        return A

    def _read_complex_sparse_ascii(self, f, nrows, ncols, line_size, line, dtype, is_big_mat):
        rows = []
        cols = []
        entries = []
        nloops = 0
        was_broken = False
        while 1:
            if nloops > 0 and not was_broken:
                line = f.readline().rstrip()
            was_broken = False

            (icol, irow, nwords) = line.split()
            icol = int(icol)

            if icol > ncols:
                break

            irow = int(irow)
            nwords = int(nwords)

            run_loop = True
            sline = line.strip().split()
            while (len(sline) == 1 or len(sline) == 2) and 'E' not in line or run_loop:  # next sparse entry
                irow = self._get_irow_sparse_ascii(f, line, sline, nwords, irow,
                                                   is_big_mat)
                run_loop = False

                i = 0
                is_real = True
                is_done_reading_row = False
                while nwords:
                    n = 0
                    line = f.readline().rstrip()
                    nwords_in_line = line.count('E')
                    if nwords_in_line == 0:
                        was_broken = True
                        break

                    for i in range(nwords_in_line):
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
        f.readline()
        return A

    def _read_complex_ascii(self, f, nrows, ncols, line_size, line, dtype, is_sparse, is_big_mat):
        if is_sparse:
            A = self._read_complex_sparse_ascii(f, nrows, ncols, line_size, line, dtype, is_big_mat)
        else:
            A = self._read_complex_dense_ascii(f, nrows, ncols, line_size, line, dtype, is_big_mat)
        return A

    def _read_complex_dense_ascii(self, f, nrows, ncols, line_size, line, dtype, is_big_mat):
        A = zeros((nrows, ncols), dtype=dtype)  # Initialize a complex matrix

        nloops = 0
        was_broken = False
        while 1:
            if nloops > 0 and not was_broken:
                line = f.readline().rstrip()
            was_broken = False

            (icol, irow, nwords) = line.split()
            icol = int(icol)

            if icol > ncols:
                break

            irow = int(irow)
            nwords = int(nwords)

            run_loop = True
            sline = line.strip().split()
            while (len(sline) == 1 or len(sline) == 2) and 'E' not in line or run_loop:  # next sparse entry
                run_loop = False

                i = 0
                iword = 0
                is_done_reading_row = False
                while nwords:
                    n = 0
                    line = f.readline().rstrip()
                    nwords_in_line = line.count('E')
                    if nwords_in_line == 0:
                        was_broken = True
                        break

                    for i in range(nwords_in_line):
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
                nloops += 1

        f.readline()
        return A

    def _get_irow_sparse_ascii(self, f, line, sline, nwords, irow, is_big_mat):
        #nwords = (nwords-1)//2  # TODO this cant be right...
        sline = line.strip().split()
        if is_big_mat:
            if len(sline) == 2:
                pass
            else:
                sline = f.readline().strip().split()
            assert len(sline) == 2, 'sline=%s len(sline)=%s' % (sline, len(sline))
            (idummy, irow) = sline
            irow = int(irow)
        else:
            if len(sline) == 1:
                IS = int(line)
            else:
                IS = int(f.readline().strip())
            L = IS // 65536 - 1
            irow = IS - 65536 * (L + 1)

        return irow

#--------------------------------------------------------------------------
    def read_op4_binary(self, op4_filename, matrix_names=None, precision='default'):
        """matrix_names must be a list or None, but basically the same"""
        f = io.open(op4_filename, mode='rb')
        self.n = 0

        # get the endian
        data = f.read(4)
        (recordLengthBig,) = unpack('>i', data)
        (recordLengthLittle,) = unpack('<i', data)

        if recordLengthBig == 24:
            self.endian = '>'
        elif recordLengthLittle == 24:
            self.endian = '<'
        else:
            msg = 'a 24 could not be found as the first word...endian error\n'
            msg += "RL_Big=%s RL_Little=%s" % (
                recordLengthBig, recordLengthLittle)
            raise RuntimeError(msg)
        f.seek(0)

        i = 0
        matrices = {}
        name = 'dummyName'
        while name is not None:
            #print "i =", i
            # checks for the end of the file
            assert self.n == f.tell(), 'n=%s tell=%s' % (self.n, f.tell())
            n = self.n
            data1 = f.read(1)
            f.seek(n)
            if len(data1) == 0:
                break
            #self.show(f, 60)

            (name, form, matrix) = self._read_matrix_binary(f, precision, matrix_names)
            #print print_matrix(matrix)
            if name is not None:
                if matrix_names is None or name in matrix_names:  # save the matrix
                    matrices[name] = (form, matrix)

            #print "not f.closed = ",not f.closed,form,name
            # if not f.closed or form is not None:
            #     data = f.read(4)
            #     self.n+=4
            #     if len(data) == 0:
            #         break
            #     (record_length,) = unpack(self.endian+'i', data)
            #     print("RL = %s" % record_length)
            #     if record_length == 24:
            #         self.n-=4
            #         f.seek(self.n)
            #     else:
            #         data = f.read(4)
            #         if len(data) == 0:
            #             break
            #         (recordLength2,) = unpack(self.endian+'i', data)
            #         assert recordLength2 == 24
            #         f.seek(self.n)
            #
            i += 1
        return matrices

    def read_start_marker(self, f):
        #print '--------------------------------------'
        #self.show(f, 60)
        data = f.read(4)
        self.n += 4
        (record_length,) = unpack(self.endian + 'i', data)

        record_length = 16
        data = f.read(record_length)
        self.n += record_length
        assert self.n == f.tell(), 'n=%s tell=%s' % (self.n, f.tell())

        if record_length == 16:  # b,icol,irow,nwords,
            (a, icol, irow, nwords) = unpack(self.endian + '4i', data)
            #print("a=%s icol=%s irow=%s nwords=%s" % (a, icol, irow, nwords))
        else:
            raise NotImplementedError('record_length=%s' % record_length)
        return (a, icol, irow, nwords)

    def _get_irow_small(self, f, data):
        if len(data) == 0:
            data = f.read(4)
            self.n += 4

        IS, = unpack(self.endian + 'i', data)
        L = IS // 65536 - 1
        irow = IS - 65536 * (L + 1)
        #print("IS=%s L=%s irow=%s" % (IS, L, irow))
        #assert IS>0
        #assert L>0
        return irow, L

    def _get_irow_big(self, f, data):
        if len(data) == 0:
            data = f.read(8)
            self.n += 8
        (idummy, irow) = unpack(self.endian + '2i', data)
        #print("idummy=%s irow=%s" % (idummy, irow))
        #assert irow<100
        return (irow, idummy - 1)

    def _read_matrix_binary(self, f, precision, matrix_names):
        """reads a matrix"""
        #self.show(f, 60)
        #print("*************************")
        data = f.read(4)
        self.n += 4
        (record_length,) = unpack(self.endian + 'i', data)
        assert self.n == f.tell(), 'n=%s tell=%s' % (self.n, f.tell())
        #print("RL = %s" % record_length)

        if record_length == 24:
            data = f.read(record_length)
            self.n += record_length

            (ncols, nrows, form, Type, name) = unpack(
                self.endian + '4i8s', data)
            #print("nrows=%s ncols=%s form=%s Type=%s name=%r" % (nrows, ncols, form, Type, name))
        else:
            msg = record_length #+ self.print_block(data)
            raise NotImplementedError('record_length=%s filename=%r' % (record_length, f.name))

        name = name.strip()
        if 0:
            if Type == 1:
                print("Type = Real, Single Precision")
            elif Type == 2:
                print("Type = Real, Double Precision")
            elif Type == 3:
                print("Type = Complex, Single Precision")
            elif Type == 4:
                print("Type = Complex, Double Precision")

        if nrows < 0:  # if less than 0, big
            is_big_mat = True
            nrows = abs(nrows)
        elif nrows > 0:
            is_big_mat = False
        else:
            raise RuntimeError('unknown BIGMAT.  nRows=%s' % nrows)

        # jump forward to get if is_sparse, then jump back
        nsave = self.n
        (_a, _icol, _irow, _nwords) = self.read_start_marker(f)
        f.seek(nsave)
        self.n = nsave

        (NWV, NBW, d, dtype) = self._get_matrix_info(Type)

        is_sparse = False
        if _irow == 0:
            is_sparse = True

        assert self.n == f.tell(), 'n=%s tell=%s' % (self.n, f.tell())
        if Type in [1, 2]:  # real
            A = self._read_real_binary(f, nrows, ncols, Type, is_sparse, is_big_mat)
        elif Type in [3, 4]:  # complex
            A = self._read_complex_binary(f, nrows, ncols, Type, is_sparse, is_big_mat)
        else:
            raise TypeError("Type=%s" % Type)

        #try:
            #print_matrix(A.todense())
        #except:
            #pass

        if d in ['d', 'dd']:
            f.read(8)
            self.n += 8
        elif d in ['f', 'ff']:
            f.read(4)
            self.n += 4
        else:
            raise NotImplementedError(d)
        #f.read(record_length); self.n+=record_length
        #self.show(f, 10)
        #f.read(4); self.n+=4

        assert self.n == f.tell(), 'n=%s f.tell=%s' % (self.n, f.tell())
        return (name, form, A)

    def _get_matrix_info(self, Type):
        if Type == 1:
            NWV = 1  # number words per value
            NBW = 4
            d = 'f'
        elif Type == 2:
            NWV = 2
            NBW = 8
            d = 'd'
        elif Type == 3:
            NWV = 2
            NBW = 8
            d = 'ff'
        elif Type == 4:
            NWV = 4
            NBW = 16
            d = 'dd'
        else:
            raise RuntimeError("Type=%s" % Type)
        dtype = get_dtype(Type)
        return (NWV, NBW, d, dtype)

    def _read_real_dense_binary(self, f, nrows, ncols, Type, is_sparse, is_big_mat):
        (NWV, NBW, d, dtype) = self._get_matrix_info(Type)
        A = zeros((nrows, ncols), dtype=dtype)

        icol = -1  # dummy value so the loop starts
        while icol < ncols + 1:  # if isDense
            assert self.n == f.tell(), 'n=%s tell=%s' % (self.n, f.tell())
            (icol, irow, nwords) = self.get_markers_dense(f)
            L = nwords
            if icol == ncols + 1:
                break
            if L == -1:
                break

            record_length = 4 * nwords
            data = f.read(record_length)
            self.n += record_length
            nvalues = L // NWV
            str_values = self.endian + '%i%s' % (nvalues * len(d), d[0])
            A[irow-1:irow-1+nvalues, icol-1] = unpack(str_values, data)
        #assert self.n == f.tell(), 'n=%s tell=%s' % (self.n, f.tell())
        f.read(4)
        self.n += 4
        return A


    def _read_real_binary(self, f, nrows, ncols, Type, is_sparse, is_big_mat):
        if is_sparse:
            A = self._read_real_sparse_binary(f, nrows, ncols, Type, is_sparse, is_big_mat)
        else:
            A = self._read_real_dense_binary(f, nrows, ncols, Type, is_sparse, is_big_mat)
        return A

    def _read_real_sparse_binary(self, f, nrows, ncols, Type, is_sparse, is_big_mat):
        (NWV, NBW, d, dtype) = self._get_matrix_info(Type)
        rows = []
        cols = []
        entries = []

        data = ''
        icol = -1  # dummy value so the loop starts
        while icol < ncols + 1:  # if isDense
            #if record_length==0:
            assert self.n == f.tell(), 'n=%s tell=%s' % (self.n, f.tell())
            (icol, irow, nwords) = self.get_markers_sparse(f, is_big_mat)
            L = nwords

            if icol == ncols + 1:
                break

            if is_big_mat:
                (irow, L) = self._get_irow_big(f, data[:8])
                data = data[8:]
            else:
                (irow, L) = self._get_irow_small(f, data[:4])
                data = data[4:]

            if L == -1:
                break

            if 0:
                print("next icol")
                print("N=%s icol=%s irow=%s nwords=%s" % (
                    self.n, icol, irow, nwords))

            #if nwords==0 and is_big_mat:
                #self.n-=4; f.seek(self.n)
                #break

            record_length = 4 * nwords
            data = f.read(record_length)
            self.n += record_length
            #print("dataFormat=%s RL=%s NNext=%s" % (d, record_length, self.n))
            #if icol==ncols+1:
                #break

            nvalues = L // NWV

            str_values = self.endian + '%i%s' % (nvalues * len(d), d[0])

            i = 0
            while len(data) > 0:
                if i == 0:
                    pass
                else:
                    if is_big_mat:
                        (irow, L) = self._get_irow_big(f, data[0:8])
                        data = data[8:]
                    else:
                        (irow, L) = self._get_irow_small(f, data[0:4])
                        data = data[4:]
                    assert irow > 0
                    nvalues = L // NWV

                value_list = unpack(str_values, data[0:nvalues * NBW])
                assert self.n == f.tell(), 'n=%s tell=%s' % (self.n, f.tell())

                #irow -= 1
                #icol -= 1
                #print "is_sparse = ",is_sparse

                rows.extend([i+irow-1 for i in range(nvalues)])
                irow += nvalues
                cols.extend([icol-1] * nvalues)
                entries.extend(value_list)

                record_length -= nvalues * NBW
                data = data[nvalues * NBW:]
                #if 0:
                    #print("record_length=%s NBW=%s len(data)=%s" %
                          #(record_length, NBW, len(data)))
                    ##print(A)
                    #print("********")  # ,data
                    #print(self.print_block(data))
                i += 1
            #print "-------------------------------"

        #if rows == []:  # NULL matrix
            #raise NotImplementedError()

        A = coo_matrix((entries, (rows, cols)), shape=(nrows, ncols),
                       dtype=dtype)
        f.read(4)
        self.n += 4

        return A

    def show(self, f, n):
        #assert self.n == f.tell()
        nold = f.tell()
        nints = n // 4
        data = f.read(n)
        strings, ints, floats = self._show_data(data)
        f.seek(nold)
        return strings, ints, floats

    def _show_data(self, data):
        n = len(data)
        nints = n // 4
        strings = unpack(self.endian + '%is' % n, data)
        ints = unpack(self.endian + '%ii' % nints, data)
        floats = unpack(self.endian + '%if' % nints, data)
        #print("strings =", strings)
        print("ints    =", ints)
        print("floats  =", floats)
        return strings, ints, floats

    def _read_complex_dense_binary(self, f, nrows, ncols, Type, is_big_mat):
        (NWV, NBW, d, dtype) = self._get_matrix_info(Type)
        A = zeros((nrows, ncols), dtype=dtype)
        record_length = 0
        icol = -1  # dummy value so the loop starts
        while icol < ncols + 1:  # if isDense
            #if record_length==0:
            assert self.n == f.tell(), 'n=%s tell=%s' % (self.n, f.tell())
            (icol, irow, nwords) = self.get_markers_dense(f)
            if 0:
                print("N=%s icol=%s irow=%s nwords=%s" % (self.n, icol, irow,
                                                          nwords))
                print("-----------")

            L = nwords

            if icol == ncols + 1:
                break

            if L == -1:
                break

            if 0:
                print("next icol")
                print("N=%s icol=%s irow=%s nwords=%s" % (
                    self.n, icol, irow, nwords))

            #if nwords==0 and is_big_mat:
                #self.n-=4; f.seek(self.n)
                #break

            record_length = 4 * nwords
            data = f.read(record_length)
            self.n += record_length
            #print("dataFormat=%s RL=%s NNext=%s" % (d, record_length, self.n))
            if icol == ncols + 1:
                continue

            nvalues = nwords // NWV
            #nRead = nwords//4
            while record_length >= NBW:
                if 0:
                    print("inner while...")
                    print("nwords  = %s" % nwords)
                    print("nvalues = %s" % nvalues)
                    print("NWV     = %s" % NWV)

                #if nvalues == 0:
                    #assert icol == ncols+1
                    #break
                str_values = self.endian + '%i%s' % (nvalues * len(d), d[0])
                if 0:
                    print("str_values = %s" % str_values)
                    print("nvalues*NBW=%s len(data)=%s" %
                          (nvalues * NBW, len(data)))
                value_list = unpack(str_values, data[0:nvalues * NBW])
                assert self.n == f.tell(), 'n=%s tell=%s' % (self.n, f.tell())
                #self.show(f, 4)
                #print self.print_block(data)
                if 0:
                    print("value_list = %s" % (value_list))

                #irow-=1
                #icol-=1
                #print("is_sparse = ", is_sparse)
                irow -= 1
                icol -= 1

                for i, value in enumerate(value_list):
                    if i % 2 == 0:
                        real_value = value
                    else:
                        #print("A[%s,%s] = %s" % (irow, icol, complex(real_value, value)))
                        A[irow, icol] = complex(real_value, value)
                        irow += 1

                record_length -= nvalues * NBW
                data = data[nvalues * NBW:]
                #print("record_length=%s NBW=%s" % (record_length, NBW))
                #print(print_matrix(A))
                #print("********", data)

        f.read(4)
        self.n += 4
        return A

    def _read_complex_binary(self, f, nrows, ncols, Type, is_sparse, is_big_mat):
        if is_sparse:
            A = self._read_complex_sparse_binary(f, nrows, ncols, Type, is_big_mat)
        else:
            A = self._read_complex_dense_binary(f, nrows, ncols, Type, is_big_mat)
        return A

    def _read_complex_sparse_binary(self, f, nrows, ncols, Type, is_big_mat):
        (NWV, NBW, d, dtype) = self._get_matrix_info(Type)
        rows = []
        cols = []
        entries = []
        record_length = 0
        data = ''
        icol = -1  # dummy value so the loop starts
        while icol < ncols + 1:  # if isDense
            #if record_length==0:
            assert self.n == f.tell(), 'n=%s tell=%s' % (self.n, f.tell())
            (icol, irow, nwords) = self.get_markers_sparse(f, is_big_mat)
            if 0:
                print("N=%s icol=%s irow=%s nwords=%s" % (self.n, icol, irow,
                                                          nwords))
                print("-----------")

            L = nwords

            if icol == ncols + 1:
                break

            if is_big_mat:
                (irow, L) = self._get_irow_big(f, data[:8])
                data = data[8:]
            else:
                (irow, L) = self._get_irow_small(f, data[:4])
                data = data[4:]

            if L == -1:
                break

            if 0:
                print("next icol")
                print("N=%s icol=%s irow=%s nwords=%s" % (
                    self.n, icol, irow, nwords))

            #if nwords==0 and is_big_mat:
                #self.n-=4; f.seek(self.n)
                #break

            record_length = 4 * nwords
            data = f.read(record_length)
            self.n += record_length
            #print("dataFormat=%s RL=%s NNext=%s" % (d, record_length, self.n))
            if icol == ncols + 1:
                continue

            nvalues = nwords // NWV
            #nRead = nwords//4
            while record_length >= NBW:
                if 0:
                    print("inner while...")
                    print("nwords  = %s" % nwords)
                    print("nvalues = %s" % nvalues)
                    print("NWV     = %s" % NWV)

                #if nvalues == 0:
                    #assert icol==ncols+1
                    #break
                str_values = self.endian + '%i%s' % (nvalues * len(d), d[0])
                if 0:
                    print("str_values = %s" % str_values)
                    print("nvalues*NBW=%s len(data)=%s" %
                          (nvalues * NBW, len(data)))
                value_list = unpack(str_values, data[0:nvalues * NBW])
                assert self.n == f.tell(), 'n=%s tell=%s' % (self.n, f.tell())
                #self.show(f, 4)
                #print(self.print_block(data))
                if 0:
                    print("value_list = %s" % (value_list))

                #irow-=1
                #icol-=1
                irow -= 1
                icol -= 1

                cols += [icol] * nvalues
                rows += [i + irow for i in range(nvalues)]
                for i, value in enumerate(value_list):
                    if i % 2 == 0:
                        real_value = value
                    else:
                        #print("A[%s,%s] = %s" % (irow, icol, complex(real_value, value)))
                        #A[irow, icol] = complex(real_value, value)
                        entries.append(complex(real_value, value))
                        irow += 1

                record_length -= nvalues * NBW
                data = data[nvalues * NBW:]
                #print("record_length=%s NBW=%s" % (record_length, NBW))
                #print(print_matrix(A))
                #print("********", data)

        A = coo_matrix((entries, (rows, cols)), shape=(nrows, ncols), dtype=dtype)
        f.read(4)
        self.n += 4
        return A

    def get_markers_sparse(self, f, is_big_mat):
        if is_big_mat:
            (a, icol, irow, nwords) = self.read_start_marker(f)
            #irow = self._get_irow_big(f)
            nwords -= 2
            #if nwords>1:
            #    nwords -= 2
            #else:
            #    print("nwords0 = %s" % nwords)
            #    nwords = 0
        else:
            (a, icol, irow, nwords) = self.read_start_marker(f)
            if irow != 0:
                assert nwords == 1, 'nwords=%s' % nwords

            #(irow) = self._get_irow_small(f)
            nwords -= 1
        return (icol, irow, nwords)

    def get_markers_dense(self, f):
        (a, icol, irow, nwords) = self.read_start_marker(f)
        #print("N=%s a=%s icol=%s irow=%s nwords=%s"% (self.n, a, icol, irow, nwords))
        return (icol, irow, nwords)

#--------------------------------------------------------------------------
    def _get_type_nwv(self, A, precision='default'):
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
            1 = real, single precision
            2 = real, double precision
            3 = complex, single precision
            4 = complex, double precision
        NWV : int
            Number of Words per Value

        .. note:: a word is 4 bytes
                  nwords(float32)=1;    single precison
                  nwords(complex64)=2;  single precison
                  nwords(float64)=2;    double precision
                  nwords(complex128)=4; double precision
        """
        # real
        if isinstance(A.dtype.type(), float32):
            NWV = 1
            if precision != 'double':
                Type = 1
            else:
                Type = 2
        elif isinstance(A.dtype.type(), float64):
            NWV = 1
            if precision != 'single':
                Type = 2
            else:
                Type = 1

        # complex
        elif isinstance(A.dtype.type(), complex64):
            NWV = 2
            if precision != 'double':
                Type = 3
            else:
                Type = 4
        elif isinstance(A.dtype.type(), complex128):
            NWV = 2
            if precision != 'single':
                Type = 4
            else:
                Type = 3
        else:
            raise TypeError('invalid Type, only float32, float64, complex64, complex128; dtype=%r' % A.dtype)
        return (Type, NWV)

    def write_op4(self, op4_filename, matrices, name_order=None, precision='default', is_binary=True):
        """
        Writes the OP4

        :param op4_filename: The filename to write
        :type op4_filename:  String -> opens a file (closed at the end)
                             file   -> no file is opened and it's not closed

        Method 1:
        ---------
        :param name_order:  List of the names of the matrices that should be
                            written or string (default=None -> sorted based on matrices)
        :param is_binary: Should a binary file be written (True, False)
        :param precision: Overwrite the default precision ('single', 'double', 'default')
                          Applies to all matrices

        Method 2:
        ---------
        matrices = {
            'A' : (formA, matrixA),
            'B' : (formB, matrixB),
        }

        .. todo::  This method is not even close to being done
        """
        if not precision in ('single', 'double', 'default'):
            raise ValueError("precision=%r and must be 'single', 'double', or 'default'" % precision)
        if not is_binary in (True, False):
            raise ValueError('is_binary=%r and must be True or False' % is_binary)
        #if nR == nC: op4_form = 1   # square
        #else:        op4_form = 2   # rectangular

        if isinstance(op4_filename, string_types):
            if PY2 or is_binary:
                f = open(op4_filename, 'wb')
            else:
                f = open(op4_filename, 'w')
            with open(op4_filename, 'w') as f:
                self._write_op4_f(f, name_order, is_binary, precision, matrices)
        else:
            f = op4_filename
            self._write_op4_f(f, name_order, is_binary, precision, matrices)

    def _write_op4_f(self, f, name_order, is_binary, precision, matrices):
        if name_order is None:
            name_order = sorted(matrices.keys())
        elif isinstance(name_order, string_types):
            name_order = [name_order]

        is_big_mat = False  ## .. todo:: hardcoded
        for name in name_order:
            (form, matrix) = matrices[name]
            if not form in (1, 2, 3, 6, 8, 9):
                raise ValueError('form=%r and must be in [1, 2, 3, 6, 8, 9]' % form)

            if isinstance(matrix, coo_matrix):
                #write_DMIG(f, name, matrix, form, precision='default')
                if is_binary:
                    raise NotImplementedError('sparse binary op4 writing not implemented')
                else:
                    self._write_sparse_matrix_ascii(f, name, matrix, form=form,
                                                    precision=precision, is_big_mat=is_big_mat)
            elif isinstance(matrix, ndarray):
                if is_binary:
                    self._write_dense_matrix_binary(f, name, matrix, form=form, precision=precision)
                else:
                    self._write_dense_matrix_ascii(f, name, matrix, form=form, precision=precision)
            else:
                raise NotImplementedError('Matrix type=%r is not supported.  types=[coo_matrix, ndarray]' % type(matrix))


    def __backup(self, name, matrix, form=2, precision='default'):
        """
        Put this documentation somewhere else...

        :param name: the name of the matrix
        :param matrix: a two-dimensional NUMPY.NDARRAY
        :param form: Form is defined as one of the following:

        ==== ===============
        Form Definition
        ==== ===============
        1.   Square
        2.   Rectangular
        3.   Diagonal
        6.   Symmetric
        8.   Id entity
        9.   Pseudoidentity
        ==== ===============

        Not Supported by all OP4s (this is not a restriction of the OP4 reader/writer):

        ==== ===============================
        Form Definition
        ==== ===============================
        4.   Lower triangular factor
        5.   Upper triangular factor
        10.  Cholesky factor
        11.  Trapezoidal factor
        13.  Sparse lower triangular factor
        15.  Sparse upper triangular factor
        ==== ===============================

        .. note:: form defaults to 2, but 1 can be easily determined.  Any others must be specified.
        """
        assert isinstance(name, string_types), name
        assert isinstance(form, int), form

    def _write_sparse_matrix_ascii(self, f, name, A, form=2, is_big_mat=False, precision='default'):
        """
        .. todo:: Does this work for complex matrices?
        """
        msg = ''
        assert isinstance(name, string_types), 'name=%s' % name
        #A = A.tolil() # list-of-lists sparse matrix
        #print dir(A)
        (Type, NWV) = self._get_type_nwv(A.data[0], precision)
        if Type in [3, 4]:
            complex_factor = 2
        else: # 1, 2
            complex_factor = 1
        (nrows, ncols) = A.shape

        #if nrows == ncols and form == 2:
            #form = 1
        #print("Type=%s" % Type)
        if is_big_mat:
            msg += '%8i%8i%8i%8i%-8s1P,3E23.16\n' % (ncols, -nrows, form, Type, name)
        else:
            msg += '%8i%8i%8i%8i%-8s1P,3E23.16\n' % (ncols, nrows, form, Type, name)

        #print "A.row = ",A.row
        #print "A.col = ",A.col

        cols = {}
        for j in A.col:
            cols[j] = []
        for i, jcol in enumerate(A.col):
            cols[jcol].append(i)
        #print("cols = ", cols)

        f.write(msg)
        msg = ''
        for j, col in iteritems(cols):
            #print("***********")
            #print("j=%s col=%s" % (j, col))
            #col.sort()

            #print('A =', A)
            irows = [A.row[jj] for jj in col]
            #print "irows = ",irows
            (dpacks) = compress_column(irows)
            #print("dpacks = %s" % (dpacks))

            npacks = len(dpacks)
            nrows = len(irows)
            if is_big_mat:
                #L = complex_factor * (2 * len(irows)) + 1
                L = 2 * npacks * NWV + nrows
                msg = '%8i%8i%8i\n' % (j + 1, 0, L)
            else:
                L = complex_factor * (2 * len(irows))
                msg = '%8i%8i%8i\n' % (j + 1, 0, L + 1)
            f.write(msg)

            for (ipack, dpack) in enumerate(dpacks):
                msg = ''
                #print "pack = ",pack

                irow = A.row[col[dpack[0]]]
                if is_big_mat:
                    #L = complexFactor * (2 * len(pack)) + 1
                    #L = (nPacks+1) + nRows*complexFactor
                    L = (len(dpack) + 1) * NWV
                    #if iPack==0:
                        #L+=1

                    #L = complexFactor * (2 + nPacks) + 1
                    #L = len(pack) + complexFactor * 2
                    #msg = '%8i%8i%8i\n' % (j+1, 0, L+1)
                    msg += '%8i%8i\n' % (L, irow + 1)
                else:
                    #L = complexFactor * (2 * len(pack))
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

                    if Type in [1, 2]:
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
                f.write(msg)
        f.write('%8i%8i%8i\n' % (ncols + 1, 1, 1))
        f.write(' 1.0000000000000000E+00\n')

    def _write_dense_matrix_binary(self, f, name, matrix, form=2,
                                   precision='default', encoding='utf-8'):
        """
        24 bytes is the record length
          - (1) ncols
          - (2) nrows
          - (3) form
          - (4) Type
          - (5,6) name2
          6 words * 4 bytes/word = 24 bytes

        .. todo:: support precision
        """
        if PY3:
            raise NotImplementedError('_write_dense_matrix_binary is not supported')
        A = matrix
        (Type, NWV) = self._get_type_nwv(A[0, 0], precision)

        (nrows, ncols) = A.shape
        #if nrows==ncols and form==2:
            #form = 1
        name2 = '%-8s' % name
        assert len(name2) == 8, 'name=%r is too long; 8 characters max' % name
        s = Struct(self.endian + '5i8s')
        msg = s.pack(24, ncols, nrows, form, Type, name2)
        f.write(msg)

        for icol in range(ncols):
            (istart, iend) = self._get_start_end_row(A[:, icol], nrows)

            # write the column
            if istart is not None and iend is not None:
                iend += 1
                msg = pack(self.endian + '4i', 24, icol +
                           1, istart + 1, (iend - istart) * NWV)

                if Type in [1, 2]: # real
                    if Type == 1: # real, single
                        fmt = '%if' % (iend - istart)
                    else:         # real, double
                        fmt = '%id' % (iend - istart)
                    f.write(pack(fmt, *A[istart:iend+1, icol]))

                else:  # complex
                    if Type == 1: # complex, single
                        fmt = '2f'
                    else:         # complex, double
                        fmt = '2d'
                    for irow in range(istart, iend):
                        msg += pack(fmt, A[irow, icol].real,
                                    A[irow, icol].imag)
            f.write(msg)
        if Type in [1, 3]: # single precision
            msg = pack(self.endian + '4if', 24, ncols + 1, 1, 1, 1.0)  # .. todo:: is this right???
        else: # double precision
            msg = pack(self.endian + '4id', 24, ncols + 1, 1, 1, 1.0)
        f.write(msg)

    def _get_start_end_row(self, A, nrows):
        """find the starting and ending points of the matrix"""
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

    def _write_dense_matrix_ascii(self, f, name, A, form=2, precision='default'):
        (Type, NWV) = self._get_type_nwv(A[0, 0], precision)

        (nrows, ncols) = A.shape
        msg = u'%8i%8i%8i%8i%-8s1P,3E23.16\n' % (ncols, nrows, form, Type, name)
        f.write(msg)

        if Type in [1, 2]: # real
            for icol in range(ncols):
                value_str = ''
                (istart, iend) = self._get_start_end_row(A[:, icol], nrows)

                # write the column
                if istart is not None and iend is not None:  # not a null column
                    iend += 1
                    msg = '%8i%8i%8i\n' % (icol + 1, istart + 1,
                                           (iend - istart) * NWV)
                    f.write(msg)
                    for i, irow in enumerate(range(istart, iend)):
                        value_str += '%23.16E' % A[irow, icol]
                        if (i + 1) % 3 == 0:
                            f.write(value_str + '\n')
                            value_str = ''
                if value_str:
                    f.write(value_str + '\n')
        else: # complex
            for icol in range(ncols):
                value_str = ''
                (istart, iend) = self._get_start_end_row(A[:, icol], nrows)

                # write the column
                if istart is not None and iend is not None:  # not a null column
                    iend += 1
                    msg = '%8i%8i%8i\n' % (icol + 1, istart + 1,
                                           (iend - istart) * NWV)
                    f.write(msg)
                    i = 0
                    for irow in range(istart, iend):
                        value_str += '%23.16E' % A[irow, icol].real
                        if (i + 1) % 3 == 0:
                            f.write(value_str + '\n')
                            value_str = ''

                        value_str += '%23.16E' % A[irow, icol].imag
                        if (i + 2) % 3 == 0:
                            f.write(value_str + '\n')
                            value_str = ''
                        i += 2
                if value_str:
                    f.write(value_str + '\n')

        # end of the matrix?
        msg = '%8i%8i%8i\n' % (ncols + 1, 1, 1)
        msg += ' 1.0000000000000000E+00\n'
        f.write(msg)


def get_dtype(Type, precision='default'):
    """reset the type if 'default' not selected"""
    if precision == 'single':
        if Type in [1, 2]:
            dtype = 'float32'
        else:
            dtype = 'complex64'
    elif precision == 'double':
        if Type in [1, 2]:
            dtype = 'float64'
        else:
            dtype = 'complex128'
    else:  # default
        if Type == 1:
            dtype = 'float32'
        elif Type == 2:
            dtype = 'float64'
        elif Type == 3:
            dtype = 'complex64'
        else:
            dtype = 'complex128'
    return dtype


def matrices():
    strings = array([[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                     [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                     [1, 0, 3, 0, 5, 0, 7, 0, 9, 0, 11, 0, 13, 0, 15, 0, 17, 0, 19, 0],
                     [1, 0, 3, 0, 5, 0, 7, 0, 9, 0, 11, 0, 13, 0, 15, 0, 17, 0, 19, 0],
                     [1, 0, 3, 0, 5, 0, 7, 0, 9, 0, 11, 0, 13, 0, 15, 0, 17, 0, 19, 0],
                     [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                     [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                     [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                     [1, 0, 3, 0, 5, 0, 7, 0, 9, 0, 11, 0, 13, 0, 15, 0, 17, 0, 19, 0],
                     [1, 0, 3, 0, 5, 0, 7, 0, 9, 0, 11, 0, 13, 0, 15, 0, 17, 0, 19, 0],
                     [1, 0, 3, 0, 5, 0, 7, 0, 9, 0, 11, 0, 13, 0, 15, 0, 17, 0, 19, 0],
                     [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                     [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                     [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                     [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20],
                     [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20],
                     [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20],
                     [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                     [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                     [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                     [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20],
                     [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20],
                     [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20],
                     [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                     [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                     [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                     [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20],
                     [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20],
                     [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20],
                     [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]], dtype='float32') # f?
    return strings


def compress_column(col):
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

def main():
    from pyNastran.op4.utils import write_DMIG

    #compress_column([14, 15, 16, 20, 21, 22, 26, 27, 28])
    filenames = [
        'test/mat_t_dn.op4',
        'test/mat_t_s1.op4',
        'test/mat_t_s2.op4',
        'test/mat_b_dn.op4',
        'test/mat_b_s1.op4',
        'test/mat_b_s2.op4',
        #'test/b_sample.op4',
        #'binary.op4',
    ]

    #matrix_names = 'EYE10' # identity
    #matrix_names = 'LOW'
    #matrix_names = 'RND1RS' # real,single
    #matrix_names = 'RND1RD' # real,double
    #matrix_names = 'RND1CS' # complex,single
    #matrix_names = 'RND1CD' # complex,double
    #matrix_names = 'STRINGS'
    #matrix_names = 'EYE5CD' # complex identity
    matrix_names = None
    #strings = matrices()

    is_big_mat = True
    if PY2:
        f = open('ascii.op4', 'wb')
    else:
        f = open('ascii.op4', 'w')
    for fname in filenames:
        op4 = OP4()
        op4.endian = '>'
        #if 't' in fname:
        #else:
            #f = open('binary.op4','wb')

        matrices = op4.read_op4(fname, matrix_names=matrix_names,
                                precision='default')
        print("keys = %s" % matrices.keys())
        print("fname=%s" % fname)
        for name, (form, matrix) in sorted(iteritems(matrices)):
            op4.write_op4(f, matrices, name_order=name)
    f.close()
    print("-----------------------------")
    print("done")
    print("-----------------------------")

if __name__ == '__main__':  # pragma: no cover
    main()
