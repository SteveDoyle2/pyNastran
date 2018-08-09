"""
Defines various tables that don't fit in other sections:
  - MinorTables
    - _get_matrix_row_fmt_nterms_nfloats(self, nvalues, tout)
    - _read_matrix(self, table_name)
    - _read_matpool_matrix(self)
    - _read_matrix_mat(self)

  - grids_comp_array_to_index(grids1, comps1, grids2, comps2,
                              make_matrix_symmetric)
"""

from __future__ import print_function, unicode_literals
from struct import unpack
from itertools import count
from six import b
import numpy as np
import scipy  # type: ignore
from pyNastran.op2.tables.matrix import Matrix
from pyNastran.op2.tables.design_response import (
    WeightResponse, StressResponse, StrainResponse, ForceResponse,
    FlutterResponse, Convergence)
from pyNastran.op2.op2_interface.op2_common import OP2Common
from pyNastran.op2.errors import FortranMarkerError


class MinorTables(OP2Common):
    """reads various tables that don't fit into a larger category"""
    def __init__(self):
        OP2Common.__init__(self)

        #: should a MATPOOL "symmetric" matrix be stored as symmetric
        #: it takes double the RAM, but is easier to use
        self.apply_symmetry = True

    def _get_matrix_row_fmt_nterms_nfloats(self, nvalues, tout):
        """
        +------+---------------------------+
        | Type | Meaning                   |
        +------+---------------------------+
        |  1   | Real, single precision    |
        |  2   | Real, double precision    |
        |  3   | Complex, single precision |
        |  4   | Complex, double precision |
        +------+---------------------------+
        """
        if tout == 1:
            nfloats = nvalues
            nterms = nvalues
            fmt = b(self._uendian + 'i %if' % nfloats)
        elif tout == 2:
            nfloats = nvalues // 2
            nterms = nvalues // 2
            fmt = b(self._uendian + 'i %id' % nfloats)
        elif tout == 3:
            nfloats = nvalues
            nterms = nvalues // 2
            fmt = b(self._uendian + 'i %if' % nfloats)
        elif tout == 4:
            nfloats = nvalues // 2
            nterms = nvalues // 4
            fmt = b(self._uendian + 'i %id' % nfloats)
        else:
            raise RuntimeError('tout = %s' % tout)
        return fmt, nfloats, nterms


    def _read_matrix(self, table_name):
        """
        general method for reading matrices and MATPOOL matrices

        .. todo:: Doesn't support checking matrices vs. MATPOOLs
        .. todo:: MATPOOLs are disabled because they're not parsed properly
        """
        i = self.f.tell()
        # if we skip on read_mode=1, we don't get debugging
        # if we just use read_mode=2, some tests fail
        #
        if self.read_mode == 2 and not self.debug_file:
            try:
                self.op2_reader._skip_matrix_mat()  # doesn't work for matpools
            except MemoryError:
                raise
            except(RuntimeError, AssertionError, ValueError):
                self.op2_reader._goto(i)
                self.op2_reader._skip_table(table_name)
            return

        try:
            self._read_matrix_mat()
        except MemoryError:
            raise
        except(RuntimeError, AssertionError, ValueError):
            # read matpool matrix
            self.op2_reader._goto(i)
            try:
                self._read_matrix_matpool()
            except(RuntimeError, AssertionError, ValueError):
                #raise
                self.op2_reader._goto(i)
                self.op2_reader._skip_table(self.table_name)

    def _read_matrix_matpool(self):
        """
        Reads a MATPOOL matrix

        MATPOOL matrices are always sparse

        +------+-----------------+
        | Form | Meaning         |
        +======+=================+
        |  1   | Square          |
        |  2   | Rectangular     |
        |  6   | Symmetric       |
        |  9   | Pseudo identity |
        +------+-----------------+
        """
        table_name = self.op2_reader._read_table_name(rewind=False, stop_on_failure=True)
        utable_name = table_name.decode('utf-8')
        self.op2_reader.read_markers([-1])
        data = self._read_record()

        self.op2_reader.read_markers([-2, 1, 0])
        data = self._read_record()

        self.op2_reader.read_markers([-3, 1, 0])
        data = self._read_record()

        #nvalues = len(data) // 4
        assert len(data) % 4 == 0, len(data) / 4.

        header = unpack(self._endian + b'3i 8s 7i', data[:48]) # 48=4*12
        assert header[:3] == (114, 1, 120), 'header[:3]=%s header=%s' % (header[:3], header)

        # ncols_gset is needed for form=9
        #  list of header values:
        #    4:5   matrix name
        #    6     placeholder
        #    7     matrix shape (1=square, 2 or 9 = rectangular, 6=symmetric)
        #    8     input type flag (1=single, 2=double, 3=complex single,
        #                           4=complex double)
        #    9     output type flag (0=precision set by system cell,
        #                            1=single, 2=double, 3=complex single,
        #                            4=complex double)
        #   10     complex flag (0=real/imaginary, >0=magnitude/phase)
        #   11     placeholder
        #   12     number of columns in the G set
        #          (only necessary for matrix shape 9)
        matrix_name, junk1, matrix_shape, tin, tout, is_phase, junk2, ncols_gset = header[3:]
        matrix_name = matrix_name.strip()

        #self.log.debug('matrix_name=%s, junk1=%s, matrix_shape=%s, tin=%s, tout=%s, '
                       #'is_phase=%s, junk2=%s, ncols_gset=%s' % (
                           #matrix_name, junk1, matrix_shape, tin, tout,
                           #is_phase, junk2, ncols_gset))

        is_complex = False
        if tin > 2 or tout > 2:
            is_complex = True
            assert is_phase == 0, 'is_phase=%s' % is_phase
            imags = []

        if tout == 1:
            dtype = 'float32'
            fdtype = self.fdtype
        elif tout == 2:
            dtype = 'float64'
            fdtype = self.double_dtype
        elif tout == 3:
            dtype = 'complex64'
            fdtype = self.fdtype
        elif tout == 4:
            dtype = 'complex128'
            fdtype = self.double_dtype
        else:
            dtype = '???'
            msg = ('matrix_name=%s, junk1=%s, matrix_shape=%s, tin=%s, tout=%s, '
                   'is_phase=%s, junk2=%s, ncols_gset=%s' % (
                       matrix_name, junk1, matrix_shape, tin, tout,
                       is_phase, junk2, ncols_gset))
            self.log.warning(msg)
            raise RuntimeError(msg)

        is_symmetric = matrix_shape == 6
        #is_phase_flag = is_phase > 0

        if tout in [1, 3]:
            # works for float32, complex64
            ints = np.frombuffer(data[48:], dtype=self.idtype).copy()
            floats = np.frombuffer(data[48:], dtype=self.fdtype).copy()
            temp_ints = ints
        else:
            # works for float64, complex128
            temp_ints = np.frombuffer(data[48:], dtype=self.idtype).copy()

        # find the first index with ()-1,-1)
        iminus1 = np.where(temp_ints[:-1] == -1)[0]
        double_minus1 = (iminus1[:-1] + 1 == iminus1[1:])[:-1]

        # the field after our stop
        # we'll handle the off by 1 later with arange
        istop = iminus1[:-2][double_minus1]

        # 2 fields after is the start position
        # add on a 0 to the beginning to account for the starting position
        # istart defines icol
        istart = np.hstack([0, istop[:-1] + 2])

        col_nids_short = temp_ints[istart]
        col_dofs_short = temp_ints[istart+1]
        #nj2 = len(istart)  ## TODO: why is this wrong???

        row_nids = []
        row_dofs = []
        col_nids = []
        col_dofs = []
        reals = []
        imags = []
        for col_nidi, col_dofi, istarti, istopi in zip(
            col_nids_short, col_dofs_short, istart + 2, istop):

            ## TODO: preallocate arrays
            imag = None
            # The float32/complex64 blocks are really simple
            # because we can just use the data is from the temp_ints block.
            # We calculate istart/istop and directly access the float data.
            #
            # The float64/complex128 blocks are more complicated.
            # It's easier to just use istart/istop to calculate datai
            # case that, and then slice it.
            #
            # In all cases, we use istart/istop to calculate row/col nid/dof
            # because they are always int32.  Only the real/imag types change.
            if dtype == 'float32':
                irow = np.arange(istarti, istopi-1, step=3, dtype='int32')
                real = floats[irow + 2]
            elif dtype == 'complex64':
                irow = np.arange(istarti, istopi-1, step=4, dtype='int32')
                real = floats[irow + 2]
                imag = floats[irow + 3]

            elif dtype == 'float64':
                datai = data[48+(istarti*4) : 48+(istopi*4)]
                irow = np.arange(istarti, istopi-1, step=4, dtype='int32')
                assert len(datai) % 8 == 0, len(datai) / 8
                real = np.frombuffer(datai, dtype=fdtype)[1::2].copy()

            elif dtype == 'complex128':
                datai = data[48+(istarti*4) : 48+(istopi*4)]

                # iword
                # -----
                #   0    1    3     5   <---- iword
                #   1    1    2     2   <---- nwords
                # (nid, dof, real, imag)
                irow = np.arange(istarti, istopi-1, step=6, dtype='int32')
                assert len(datai) % 8 == 0, len(datai) / 8
                floats = np.frombuffer(datai, dtype=fdtype).copy()

                # ndoubles
                # --------
                #  <---0--->   1     2    <----- iword
                #      1       1     1    <----- nwords
                # (nid, dof, real, imag)
                real = floats[1::3]
                imag = floats[2::3]
            else:
                msg = '%s is not supported' % dtype
                self.log.error(msg)
                raise RuntimeError(msg)

            if len(irow) != len(real):
                msg = 'nrow=%s nreal=%s nimag=%s' % (len(irow), len(real), len(imag))
                raise RuntimeError(msg)

            # the row index; [1, 2, ..., 43]
            row_nid = temp_ints[irow]

            # the dof; [0, 0, ..., 0.]
            row_dof = temp_ints[irow + 1]
            urow_dof = np.unique(row_dof)
            for udofi in urow_dof:
                if udofi not in [0, 1, 2, 3, 4, 5, 6]:
                    msg = 'udofi=%s is invalid; must be in [0, 1, 2, 3, 4, 5, 6]; dofs=%s' % (
                        udofi, np.asarray(urow_dof, dtype='int32').tolist())
                    raise ValueError(msg)

            ni = len(irow)
            col_nid = np.ones(ni, dtype='int32') * col_nidi
            col_dof = np.ones(ni, dtype='int32') * col_dofi

            row_nids.append(row_nid)
            row_dofs.append(row_dof)
            col_nids.append(col_nid)
            col_dofs.append(col_dof)
            reals.append(real)
            imags.append(imag)

        row_nids_array = np.hstack(row_nids)
        row_dofs_array = np.hstack(row_dofs)

        col_nids_array = np.hstack(col_nids)
        col_dofs_array = np.hstack(col_dofs)
        real_array = np.hstack(reals)
        if is_complex:
            complex_array = np.hstack(imags)
            assert len(real_array) == len(complex_array)
            real_imag_array = real_array + 1.j * complex_array
        else:
            real_imag_array = real_array

        self._cast_matrix_matpool(utable_name, real_imag_array,
                                  col_nids_array, col_dofs_array,
                                  row_nids_array, row_dofs_array,
                                  matrix_shape, dtype, is_symmetric)

    def _cast_matrix_matpool(self, table_name, real_imag_array,
                             col_nids_array, col_dofs_array,
                             row_nids_array, row_dofs_array,
                             matrix_shape, dtype, is_symmetric):
        """helper method for _read_matpool_matrix"""

        make_matrix_symmetric = self.apply_symmetry and matrix_shape == 'symmetric'

        # TODO: this is way slower than it should be
        #       because we didn't preallocate the data and the
        #       grids_comp_array_to_index function needs work
        grids1 = col_nids_array
        comps1 = col_dofs_array
        grids2 = row_nids_array
        comps2 = row_dofs_array
        assert len(grids1) == len(comps1), 'ngrids1=%s ncomps1=%s' % (len(grids1), len(comps1))
        assert len(grids1) == len(grids2), 'ngrids1=%s ngrids2=%s' % (len(grids1), len(grids2))
        assert len(comps1) == len(comps2), 'ncomps1=%s ncomps2=%s' % (len(comps1), len(comps2))

        j1, j2, nj1, nj2, nj = grids_comp_array_to_index(
            grids1, comps1, grids2, comps2, make_matrix_symmetric)
        assert len(j1) == len(j2), 'nj1=%s nj2=%s' % (len(j1), len(j2))
        assert len(grids1) == len(real_imag_array), 'ngrids1=%s nreals=%s' % (len(j1), len(real_imag_array))

        # not 100% on these, they might be flipped
        #ncols = len(np.unique(j1))
        #mrows = len(np.unique(j2))

        if is_symmetric:
            mrows = nj
            ncols = nj
            #print('  j1 =', j1)
            #print('  j2 =', j2)
        else:
            ncols = nj1
            mrows = nj2

        try:
            matrix = scipy.sparse.coo_matrix(
                (real_imag_array, (j2, j1)),
                shape=(mrows, ncols), dtype=dtype)
        except ValueError:
            msg = 'Passed all the checks; cannot build MATPOOL sparse matrix...\n'
            spaces = '                                          '
            msg += '%sname=%s dtype=%s nrows=%s ncols=%s nj1=%s nj2=%s nj=%s' % (
                spaces, table_name, dtype, mrows, ncols, nj1, nj2, nj)
            self.log.error(msg)
            raise


        # enforce symmetry if necessary
        if make_matrix_symmetric:
            # get the upper and lower triangular matrices
            upper_tri = scipy.sparse.triu(matrix)
            lower_tri = scipy.sparse.tril(matrix)

            # extracts a [1, 2, 3, ..., n] off the diagonal of the matrix
            # and make it a diagonal matrix
            diagi = scipy.sparse.diags(scipy.sparse.diagional(upper_tri))

            # Check to see which triangle is populated.
            # If they both are, make sure they're equal
            # or average them and throw a warning
            lnnz = (lower_tri - diagi).nnz
            unnz = (upper_tri - diagi).nnz
            assert isinstance(lnnz, int), type(lnnz)
            assert isinstance(unnz, int), type(unnz)

            # both upper and lower triangle are populated
            if lnnz > 0 and unnz > 0:
                upper_tri_t = upper_tri.T
                if lower_tri == upper_tri_t:
                    matrix = upper_tri + upper_tri_t - diagi
                else:
                    self.log.warning(
                        'Matrix %r marked as symmetric does not contain '
                        'symmetric data.  Data will be symmetrized by averaging.' % table_name)
                    matrix = (matrix + matrix.T) / 2.
            elif lnnz > 0:
                #  lower triangle is populated
                matrix = lower_tri + lower_tri.T - diagi
            elif unnz > 0:
                #  upper triangle is populated
                matrix = upper_tri + upper_tri_t - diagi
            else:
                # matrix is diagonal (or null)
                matrix = diagi
            data = matrix

            # matrix is symmetric, but is not stored as symmetric
            matrix_shape = 'rectangular'

        m = Matrix(table_name, is_matpool=True, form=matrix_shape)
        m.data = matrix
        m.col_nid = col_nids_array
        m.col_dof = col_dofs_array
        m.row_nid = row_nids_array
        m.row_dof = row_dofs_array
        m.form = matrix_shape
        self.matrices[table_name] = m
        self.log.debug(m)

        self.op2_reader.read_markers([-4, 1, 0])
        data = self._read_record()

        if len(data) == 12:
            self.op2_reader.read_markers([-5, 1, 0, 0])
            return
        raise RuntimeError('failed on _read_matpool_matrix')

    def _read_matrix_mat(self):
        """
        Reads a matrix in "standard" form.  The forms are::
            standard:
                Return a matrix that looks similar to a matrix found
                in the OP4.  Created by:
                ``ASSIGN output2='model.op2', UNIT=12,UNFORMATTED,DELETE``
                ``OUTPUT2 KGG//0/12``
            matpool:
                Return a matrix that looks similar to a DMIG matrix
                (e.g., it contains the node id and DOF).  Created by:
                ``ASSIGN output2='model.op2', UNIT=12,UNFORMATTED,DELETE``
                ``TODO: add the magic keyword...``
                ``OUTPUT2 KGG//0/12``

        Matrix Trailer:
        +------+---------------------------------------------------+
        | Word | Contents                                          |
        +======+===================================================+
        |  1   | Number of columns in matrix                       |
        |  2   | Number of rows in matrix                          |
        |  3   | Form of the matrix                                |
        |  4   | Type of matrix                                    |
        |  5   | Largest number of nonzero words among all columns |
        |  6   | Density of the matrix multiplied by 10000         |
        |  7   | Size in blocks                                    |
        |  8   | Maximum string length over all strings            |
        |  9   | Number of strings                                 |
        |  10  | Average bandwidth                                 |
        |  11  | Maximum bandwidth                                 |
        |  12  | Number of null columns                            |
        +------+---------------------------------------------------+

        +------+--------------------------------+
        | Form | Meaning                        |
        +======+================================+
        |  1   | Square                         |
        |  2   | Rectangular                    |
        |  3   | Diagonal                       |
        |  4   | Lower triangular factor        |
        |  5   | Upper triangular factor        |
        |  6   | Symmetric                      |
        |  8   | Identity                       |
        |  9   | Pseudo identity                |
        |  10  | Cholesky factor                |
        |  11  | Trapezoidal factor             |
        |  13  | Sparse lower triangular factor |
        |  15  | Sparse upper triangular factor |
        +------+--------------------------------+

        +------+---------------------------+
        | Type | Meaning                   |
        +======+===========================+
        |  1   | Real, single precision    |
        |  2   | Real, double precision    |
        |  3   | Complex, single precision |
        |  4   | Complex, double precision |
        +------+---------------------------+
        """
        allowed_forms = [1, 2, 3, 4, 5, 6, 8, 9, 10, 11, 13, 15]
        #self.log.debug('----------------------------------------------------------------')
        table_name = self.op2_reader._read_table_name(rewind=False, stop_on_failure=True)
        self.op2_reader.read_markers([-1])
        data = self._read_record()

        # old-bad
        #matrix_num, form, mrows, ncols, tout, nvalues, g = unpack(self._endian + b'7i', data)

        #           good   good   good  good  ???    ???
        matrix_num, ncols, mrows, form, tout, nvalues, g = unpack(self._endian + b'7i', data)
        #print('g =', g)

        m = Matrix(table_name, form=form)
        self.matrices[table_name.decode('utf-8')] = m

        # matrix_num is a counter (101, 102, 103, ...)
        # 101 will be the first matrix 'A' (matrix_num=101),
        # then we'll read a new matrix 'B' (matrix_num=102),
        # etc.
        #
        # the matrix is Mrows x Ncols
        #
        # it has nvalues in it
        #
        # tout is the precision of the matrix
        # 0 - set precision by cell
        # 1 - real, single precision (float32)
        # 2 - real, double precision (float64)
        # 3 - complex, single precision (complex64)
        # 4 - complex, double precision (complex128)

        # form (bad name)
        # 1 - column matrix
        # 2 - factor matrix
        # 3 - factor matrix
        if tout == 1:
            dtype = 'float32'
        elif tout == 2:
            dtype = 'float64'
        elif tout == 3:
            dtype = 'complex64'
        elif tout == 4:
            dtype = 'complex128'
        else:
            dtype = '???'
            msg = ('unexpected tout for %s: matrix_num=%s form=%s '
                   'mrows=%s ncols=%s tout=%s nvalues=%s g=%s'  % (
                       table_name, matrix_num, form, mrows, ncols, tout, nvalues, g))
            self.log.warning(msg)
            raise RuntimeError(msg)

        #self.log.error('name=%r matrix_num=%s form=%s mrows=%s '
        #               'ncols=%s tout=%s nvalues=%s g=%s' % (
        #                   table_name, matrix_num, form, mrows, ncols, tout, nvalues, g))
        if form == 1:
            if ncols != mrows:
                self.log.warning('unexpected size for %s; form=%s mrows=%s ncols=%s' % (
                    table_name, form, mrows, ncols))
        elif form not in allowed_forms:
            self.log.error('name=%r matrix_num=%s form=%s mrows=%s '
                           'ncols=%s tout=%s nvalues=%s g=%s' % (
                               table_name, matrix_num, form, mrows, ncols, tout, nvalues, g))
            raise RuntimeError('form=%s; allowed=%s' % (form, allowed_forms))
        #self.log.debug('name=%r matrix_num=%s form=%s mrows=%s ncols=%s tout=%s nvalues=%s g=%s' % (
            #table_name, matrix_num, form, mrows, ncols, tout, nvalues, g))

        self.op2_reader.read_markers([-2, 1, 0])
        data = self._read_record()

        if len(data) == 16:
            name, ai, bi = unpack(self._endian + b'8s 2i', data)
            assert ai == 170, ai
            assert bi == 170, bi
        else:
            self.log.warning('unexpected matrix length=%s' % len(data))
            self.log.warning(self.show_data(data, types='if'))

        itable = -3
        j = None

        niter = 0
        niter_max = 100000000

        GCi = []
        GCj = []
        reals = []
        jj = 1
        while niter < niter_max:
            #nvalues = self.op2_reader.get_marker1(rewind=True)
            self.op2_reader.read_markers([itable, 1])
            one = self.op2_reader.get_marker1(rewind=False)

            if one:  # if keep going
                nvalues = self.op2_reader.get_marker1(rewind=True)

                while nvalues >= 0:
                    nvalues = self.op2_reader.get_marker1(rewind=False)
                    fmt, nfloats, nterms = self._get_matrix_row_fmt_nterms_nfloats(nvalues, tout)
                    GCjj = [jj] * nterms
                    GCj += GCjj

                    #-----------
                    data = self.read_block()
                    #self.show_data(data)
                    #print('***itable=%s nvalues=%s fmt=%r' % (itable, nvalues, fmt))
                    out = unpack(fmt, data)
                    ii = out[0]
                    values = out[1:]

                    GCii = list(range(ii, ii + nterms))
                    GCi += GCii
                    reals += values
                    nvalues = self.op2_reader.get_marker1(rewind=True)
                    if self.debug_file:
                        self.binary_debug.write('  GCi = %s\n' % GCii)
                        self.binary_debug.write('  GCj = %s\n' % GCjj)
                        self.binary_debug.write('  reals/imags = %s\n' % str(values))
                assert len(GCi) == len(GCj), 'nGCi=%s nGCj=%s' % (len(GCi), len(GCj))
                if tout in [1, 2]:
                    assert len(GCi) == len(reals), 'tout=%s nGCi=%s nreals=%s' % (tout, len(GCi), len(reals))
                else:
                    assert len(GCi)*2 == len(reals), 'tout=%s nGCi=%s nreals=%s' % (tout, len(GCi)*2, len(reals))
                jj += 1
            else:
                nvalues = self.op2_reader.get_marker1(rewind=False)
                assert nvalues == 0, nvalues

                matrix = self._cast_matrix_mat(GCi, GCj, mrows, ncols, reals, tout, dtype)
                m.data = matrix
                if matrix is not None:
                    self.matrices[table_name.decode('utf-8')] = m
                #nvalues = self.op2_reader.get_marker1(rewind=True)
                return
            itable -= 1
            niter += 1
        raise RuntimeError('MaxIteration: this should never happen; n=%s' % niter_max)

    def _cast_matrix_mat(self, GCi, GCj, mrows, ncols, reals, tout, dtype):
        """helper method for _read_matrix_mat"""
        #assert max(GCi) <= mrows, 'GCi=%s GCj=%s mrows=%s' % (GCi, GCj, mrows)
        #assert max(GCj) <= ncols, 'GCi=%s GCj=%s ncols=%s' % (GCi, GCj, ncols)

        # we subtract 1 to the indicides to account for Fortran
        GCi = np.array(GCi, dtype='int32') - 1
        GCj = np.array(GCj, dtype='int32') - 1
        try:
            if dtype == '???':
                matrix = None
                self.log.warning('what is the dtype?')
            elif tout in [1, 2]:
                # real
                real_array = np.array(reals, dtype=dtype)
                matrix = scipy.sparse.coo_matrix(
                    (real_array, (GCi, GCj)),
                    shape=(mrows, ncols), dtype=dtype)
                #self.log.info('created %s (real)' % self.table_name)
            elif tout in [3, 4]:
                # complex
                real_array = np.array(reals, dtype=dtype)
                nvalues_matrix = real_array.shape[0] // 2
                real_complex = real_array.reshape((nvalues_matrix, 2))
                real_imag = real_complex[:, 0] + real_complex[:, 1]*1j
                if self.binary_debug:
                    #self.binary_debug.write('reals = %s' % real_complex[:, 0])
                    #self.binary_debug.write('imags = %s' % real_complex[:, 1])
                    self.binary_debug.write('real_imag = %s' % real_imag)
                matrix = scipy.sparse.coo_matrix(
                    (real_imag, (GCi, GCj)),
                    shape=(mrows, ncols), dtype=dtype)
                #msg = 'created %s (complex)' % self.table_name
                #self.log.debug(msg)
                #raise RuntimeError(msg)
            else:
                raise RuntimeError('this should never happen')
        except ValueError:
            self.log.warning('shape=(%s, %s)' % (mrows, ncols))
            self.log.warning('cant make a coo/sparse matrix...trying dense')

            if dtype == '???':
                matrix = None
                self.log.warning('what is the dtype?')
            else:
                real_array = np.array(reals, dtype=dtype)
                self.log.debug('shape=%s mrows=%s ncols=%s' % (
                    str(real_array.shape), mrows, ncols))
                if len(reals) == mrows * ncols:
                    real_array = real_array.reshape(mrows, ncols)
                    self.log.info('created %s' % self.table_name)
                else:
                    self.log.warning('cant reshape because invalid sizes : created %s' %
                                     self.table_name)

                matrix = real_array
        return matrix

def grids_comp_array_to_index(grids1, comps1, grids2, comps2,
                              make_matrix_symmetric):
    """maps the dofs"""
    #from pyNastran.utils.mathematics import unique2d
    ai = np.vstack([grids1, comps1]).T
    bi = np.vstack([grids2, comps2]).T
    #print('grids2 =', grids2)
    #print('comps2 =', comps2)
    #c = np.vstack([a, b])
    #assert c.shape[1] == 2, c.shape
    #urows = unique2d(c)
    #urows = c

    nid_comp_to_dof_index = {}
    j = 0
    a_keys = set()
    for nid_dof in ai:
        #nid_dof = (int(nid), int(dof))
        nid_dof = tuple(nid_dof)
        if nid_dof not in a_keys:
            a_keys.add(nid_dof)
            if nid_dof not in nid_comp_to_dof_index:
                nid_comp_to_dof_index[nid_dof] = j
                j += 1
    nja = len(a_keys)
    del a_keys

    b_keys = set()
    for nid_dof in bi:
        nid_dof = tuple(nid_dof)
        if nid_dof not in b_keys:
            b_keys.add(nid_dof)
        if nid_dof not in nid_comp_to_dof_index:
            nid_comp_to_dof_index[nid_dof] = j
            j += 1
    njb = len(b_keys)
    del b_keys


    nj = len(nid_comp_to_dof_index)
    if make_matrix_symmetric:
        ja = np.zeros(nj, dtype='int32')
        for i, nid_dof in zip(count(), ai):
            j[i] = nid_comp_to_dof_index[tuple(nid_dof)]
        for i, nid_dof in zip(count(), bi):
            j[i] = nid_comp_to_dof_index[tuple(nid_dof)]
        return j, j, nj, nj, nj
    else:
        ja = np.zeros(grids1.shape, dtype='int32')
        for i, nid_dof in zip(count(), ai.tolist()):
            ja[i] = nid_comp_to_dof_index[tuple(nid_dof)]

        jb = np.zeros(grids2.shape, dtype='int32')
        for i, nid_dof in zip(count(), bi.tolist()):
            jb[i] = nid_comp_to_dof_index[tuple(nid_dof)]

        return ja, jb, nja, njb, nj
