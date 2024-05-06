from __future__ import annotations
from struct import unpack
from typing import TYPE_CHECKING

import numpy as np
import scipy  # type: ignore

#from pyNastran.utils.numpy_utils import integer_types

from pyNastran.op2.op2_interface.read_matrix_matpool import read_matrix_matpool
from pyNastran.op2.result_objects.matrix import Matrix
#from pyNastran.op2.result_objects.matrix_dict import MatrixDict

from pyNastran.op2.op2_interface.utils import (
    mapfmt, # reshape_bytes_block,
    #reshape_bytes_block_size,
)


if TYPE_CHECKING:
    from cpylog import SimpleLogger
    from pyNastran.op2.op2_interface.op2_reader import OP2Reader
    from pyNastran.op2.op2 import OP2


DENSE_MATRICES = [
    b'KELM', b'MELM', b'BELM',
    b'KELMP', b'MELMP',

    b'EFMASSS', b'EFMFACS', b'EFMFSMS',
    b'MEFMASS', b'MEFWTS', b'MPFACS', b'RBMASSS',
]
def read_matrix_mat(op2_reader: OP2Reader) -> None:
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
    op2: OP2 = op2_reader.op2
    log: SimpleLogger = op2_reader.log
    endian = op2_reader._endian
    size = op2_reader.size

    allowed_forms = [1, 2, 3, 4, 5, 6, 8, 9, 10, 11, 13, 15]
    #op2_reader.log.debug('----------------------------------------------------------------')
    table_name = op2_reader._read_table_name(rewind=False, stop_on_failure=True)
    op2_reader.read_markers([-1])
    data = op2_reader._read_record()

    # old-bad
    #matrix_num, form, mrows, ncols, tout, nvalues, g = unpack(endian + b'7i', data)

    fmt1 = mapfmt(endian + b'7i', size)
    #           good   good   good  good  ???    ???
    matrix_num, ncols, mrows, form, tout, nvalues, g = unpack(fmt1, data)
    #print('g =', g)

    utable_name = table_name.decode('utf-8')

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
        log.warning(msg)
        raise RuntimeError(msg)

    m = Matrix(utable_name, form=form)
    op2.matrices[utable_name] = m

    #op2_reader.log.error('name=%r matrix_num=%s form=%s mrows=%s '
    #               'ncols=%s tout=%s nvalues=%s g=%s' % (
    #                   table_name, matrix_num, form, mrows, ncols, tout, nvalues, g))
    if form == 1:
        if ncols != mrows:
            log.warning('unexpected size for %s; form=%s mrows=%s ncols=%s' % (
                table_name, form, mrows, ncols))
    elif form not in allowed_forms:
        log.error('name=%r matrix_num=%s form=%s mrows=%s '
                  'ncols=%s tout=%s nvalues=%s g=%s' % (
                      table_name, matrix_num, form, mrows, ncols,
                      tout, nvalues, g))
        raise RuntimeError('form=%s; allowed=%s' % (form, allowed_forms))

    if size == 4:
        log.debug('name=%r matrix_num=%s form=%s mrows=%s ncols=%s tout=%s '
                  'nvalues=%s g=%s' % (
                      table_name, matrix_num, form, mrows, ncols, tout, nvalues, g))
    else:
        #if tout == 1:
            #tout = 2
        log.info('name=%r matrix_num=%s form=%s mrows=%s ncols=%s tout=%s '
                 'nvalues=%s g=%s' % (
                     table_name, matrix_num, form, mrows, ncols, tout, nvalues, g))

    op2_reader.read_3_markers([-2, 1, 0])
    data = op2_reader._read_record()
    if size == 4:
        if len(data) == 16:
            unused_name, ai, bi = unpack(endian + b'8s 2i', data)
            assert ai == 170, ai
            assert bi == 170, bi
        else:
            log.warning('unexpected matrix length=%s' % len(data))
            #op2_reader.log.warning(op2_reader.show_data(data, types='if'))
    elif size == 8:
        if len(data) == 32:
            unused_name, ai, bi = unpack(endian + b'16s 2q', data)
            # name isn't mapped
            assert ai == 170, ai
            assert bi == 170, bi
        else:
            log.warning('unexpected matrix length=%s' % len(data))
            #op2_reader.log.warning(op2_reader.show_data(data, types='ifsqd', endian=None))
    else:
        raise RuntimeError(size)

    itable = -3
    unused_j = None

    niter = 0
    niter_max = 100000000

    GCi = []
    GCj = []
    reals = []
    jj = 1
    while niter < niter_max:
        #nvalues = op2_reader.get_marker1(rewind=True)
        op2_reader.read_markers([itable, 1])
        one = op2_reader.get_marker1(rewind=False)

        if one:  # if keep going
            nvalues = op2_reader.get_marker1(rewind=True)

            while nvalues >= 0:
                nvalues = op2_reader.get_marker1(rewind=False)
                fmt, unused_nfloats, nterms = _get_matrix_row_fmt_nterms_nfloats(
                    nvalues, tout, endian)
                GCjj = [jj] * nterms
                GCj += GCjj

                #-----------
                data = op2_reader._read_block()
                if size == 8:
                    #op2_reader.log.warning('skipping matrix')
                    #op2_reader.show_data(data, types='ifqd')
                    fmt = mapfmt(fmt, size)
                    #op2_reader.log.warning(fmt)
                    #print('***itable=%s nvalues=%s fmt=%r' % (itable, nvalues, fmt))
                    #continue
                out = unpack(fmt, data)
                #print(out)
                ii = out[0]
                values = out[1:]

                #list(range(2, 10))
                #[2, 3, 4, 5, 6, 7, 8, 9]
                GCii = list(range(ii, ii + nterms))
                GCi += GCii
                reals += values
                nvalues = op2_reader.get_marker1(rewind=True)
                if op2_reader.debug_file:
                    op2_reader.binary_debug.write('  GCi = %s\n' % GCii)
                    op2_reader.binary_debug.write('  GCj = %s\n' % GCjj)
                    op2_reader.binary_debug.write('  reals/imags = %s\n' % str(values))
            assert len(GCi) == len(GCj), 'nGCi=%s nGCj=%s' % (len(GCi), len(GCj))
            if size == 4:
                if tout in [1, 2]:
                    assert len(GCi) == len(reals), 'tout=%s nGCi=%s nreals=%s' % (tout, len(GCi), len(reals))
                else:
                    assert len(GCi)*2 == len(reals), 'tout=%s nGCi=%s nreals=%s' % (tout, len(GCi)*2, len(reals))
            jj += 1
        else:
            nvalues = op2_reader.get_marker1(rewind=False)
            assert nvalues == 0, nvalues

            matrix = _cast_matrix_mat(GCi, GCj, mrows, ncols, reals, tout, dtype, log)
            if table_name in DENSE_MATRICES:
                matrix = matrix.toarray()
            m.data = matrix
            if matrix is not None:
                op2.matrices[table_name.decode('utf-8')] = m
            #nvalues = op2_reader.get_marker1(rewind=True)
            return
        itable -= 1
        niter += 1
    raise RuntimeError('MaxIteration: this should never happen; n=%s' % niter_max)

def _skip_matrix_mat(op2_reader: OP2Reader) -> None:
    """
    Reads a matrix in "standard" form.

    Notes
    -----
    see read_matrix_mat

    """
    unused_table_name = op2_reader._read_table_name(rewind=False, stop_on_failure=True)
    op2_reader.read_markers([-1])
    unused_data = op2_reader._skip_record()

    op2_reader.read_3_markers([-2, 1, 0])
    unused_data = op2_reader._skip_record()

    itable = -3
    niter = 0
    niter_max = 100000000

    #jj = 1
    while niter < niter_max:
        #nvalues = op2_reader.get_marker1(rewind=True)
        #print('nvalues4a =', nvalues)
        op2_reader.read_markers([itable, 1])
        one = op2_reader.get_marker1(rewind=False)

        if one:  # if keep going
            nvalues = op2_reader.get_marker1(rewind=True)
            while nvalues >= 0:
                nvalues = op2_reader.get_marker1(rewind=False)
                unused_data = op2_reader._skip_block()
                nvalues = op2_reader.get_marker1(rewind=True)
            #jj += 1
        else:
            nvalues = op2_reader.get_marker1(rewind=False)
            assert nvalues == 0, nvalues
            return
        itable -= 1
        niter += 1
    raise RuntimeError('this should never happen; n=%s' % niter_max)

def read_matrix(op2_reader, table_name: bytes) -> None:
    """
    General method for reading matrices and MATPOOL matrices

    Note
    ----
    Matrices are read on read_mode = 1

    .. todo:: Doesn't support checking matrices vs. MATPOOLs

    """
    # it'd be nice to read in read_mode=2 just for faster dev time
    # for adding new tables when matrices exist, but the code fails
    #read_mode_to_read_matrix = 1

    #op2: OP2 = op2_reader.op2
    #i = op2.f.tell()
    # if we skip on read_mode=1, we don't get debugging
    # if we just use read_mode=2, some tests fail
    #
    mat_type = _check_matrix_type(op2_reader)
    if mat_type == 'matrix':
        read_matrix_mat(op2_reader)
    else:
        read_matrix_matpool(op2_reader)

    return
    #from traceback import format_exc
    #if op2.read_mode != read_mode_to_read_matrix and not op2_reader.debug_file:
        #try:
            #op2_reader._skip_matrix_mat()  # doesn't work for matpools
        #except MemoryError:  # pragma: no cover
            #raise
        #except(RuntimeError, AssertionError, ValueError):
            #raise
            #op2_reader._goto(i)
            #op2_reader._skip_table(table_name)
        #return

    #try:
    #    op2_reader._read_matrix_mat()
    #    return
    #except MemoryError:  # pragma: no cover
    #    raise
    #except(RuntimeError, AssertionError, ValueError):
    #    pass # op2_reader.log.error(str(format_exc()))

    # read matpool matrix
    #op2_reader._goto(i)
    #try:
    #    read_matrix_matpool(op2_reader)
    #    return
    #except(RuntimeError, AssertionError, ValueError):
    #    op2_reader.log.error(str(format_exc()))

    # I give up
    #op2_reader._goto(i)
    #op2_reader._skip_table(op2.table_name)

def _check_matrix_type(op2_reader: OP2Reader) -> str:
    op2: OP2 = op2_reader.op2
    i = op2.f.tell()
    table_name = op2_reader._read_table_name(rewind=False, stop_on_failure=True)
    unused_utable_name = table_name.decode('utf-8')
    #print(utable_name)
    op2_reader.read_markers([-1])

    # (104, 32768, 0, 0, 0, 0, 0)
    data = op2_reader._read_record()
    ints = np.frombuffer(data, dtype=op2.idtype8)
    op2_reader._goto(i)

    zeros = ints[2:]
    if np.abs(zeros).sum() == 0:
        return 'matpool'
    return 'matrix'

def _cast_matrix_mat(GCi: np.ndarray, GCj: np.ndarray,
                     mrows: int, ncols: int,
                     reals: np.ndarray,
                     tout: int,
                     dtype: str, log: SimpleLogger) -> np.ndarray:
    """helper method for _read_matrix_mat"""
    #assert max(GCi) <= mrows, 'GCi=%s GCj=%s mrows=%s' % (GCi, GCj, mrows)
    #assert max(GCj) <= ncols, 'GCi=%s GCj=%s ncols=%s' % (GCi, GCj, ncols)

    # we subtract 1 to the indicides to account for Fortran
    GCi = np.array(GCi, dtype='int32') - 1
    GCj = np.array(GCj, dtype='int32') - 1
    try:
        if dtype == '???':
            matrix = None
            log.warning('what is the dtype?')
        elif tout in {1, 2}:
            # real
            real_array = np.array(reals, dtype=dtype)
            matrix = scipy.sparse.coo_matrix(
                (real_array, (GCi, GCj)),
                shape=(mrows, ncols), dtype=dtype)
            #log.info(f'created {op2_reader.table_name} (real)')
        elif tout in {3, 4}:
            # complex
            real_array = np.array(reals, dtype=dtype)
            nvalues_matrix = real_array.shape[0] // 2
            real_complex = real_array.reshape((nvalues_matrix, 2))
            real_imag = real_complex[:, 0] + real_complex[:, 1]*1j
            #if op2_reader.binary_debug:
                #op2_reader.binary_debug.write('reals = %s' % real_complex[:, 0])
                #op2_reader.binary_debug.write('imags = %s' % real_complex[:, 1])
                #op2_reader.binary_debug.write('real_imag = %s' % real_imag)
            matrix = scipy.sparse.coo_matrix(
                (real_imag, (GCi, GCj)),
                shape=(mrows, ncols), dtype=dtype)
            #msg = 'created %s (complex)' % op2_reader.table_name
            #log.debug(msg)
            #raise RuntimeError(msg)
        else:
            raise RuntimeError('this should never happen')
    except ValueError:
        log.warning('shape=(%s, %s)' % (mrows, ncols))
        log.warning('cant make a coo/sparse matrix...trying dense')

        if dtype == '???':
            matrix = None
            log.warning('what is the dtype?')
        else:
            real_array = np.array(reals, dtype=dtype)
            log.debug('shape=%s mrows=%s ncols=%s' % (
                str(real_array.shape), mrows, ncols))
            if len(reals) == mrows * ncols:
                real_array = real_array.reshape(mrows, ncols)
                log.info(f'created {op2.table_name}')
            else:
                log.warning(f'cant reshape because invalid sizes : created {op2.table_name}')

            matrix = real_array
    return matrix


def _get_matrix_row_fmt_nterms_nfloats(nvalues: int,
                                       tout: int,
                                       endian: bytes) -> tuple[bytes, int, int]:
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
        fmt = endian + b'i %if' % nfloats
    elif tout == 2:
        nfloats = nvalues // 2
        nterms = nvalues // 2
        fmt = endian + b'i %id' % nfloats
    elif tout == 3:
        nfloats = nvalues
        nterms = nvalues // 2
        fmt = endian + b'i %if' % nfloats
    elif tout == 4:
        nfloats = nvalues // 2
        nterms = nvalues // 4
        fmt = endian + b'i %id' % nfloats
    else:
        raise RuntimeError(f'tout = {tout}')
    return fmt, nfloats, nterms
