from __future__ import annotations
from struct import Struct
from itertools import count
from typing import Any, TYPE_CHECKING

import numpy as np
import scipy  # type: ignore

from cpylog import SimpleLogger
from pyNastran.nptyping_interface import NDArrayNint
from pyNastran.op2.result_objects.matrix import Matrix
from pyNastran.op2.op2_interface.utils import reshape_bytes_block

if TYPE_CHECKING:  # pragma: no cover
    from pyNastran.op2.op2 import OP2

def read_matpool_dmig_4(op2: OP2, data: bytes,
                        unused_utable_name: str, debug: bool=False):
    """
    ncols_gset is needed for form=9
    list of header values:
      4:5   matrix name
      6     placeholder
      7     matrix shape (1=square, 2 or 9=rectangular, 6=symmetric)
      8     input type flag (1=single, 2=double, 3=complex single,
                             4=complex double)
      9     output type flag (0=precision set by system cell,
                              1=single, 2=double, 3=complex single,
                              4=complex double)
     10     complex flag (0=real/imaginary, >0=magnitude/phase)
     11     placeholder
     12     number of columns in the G set
            (only necessary for matrix shape 9)
    """
    #aaa
    log = op2.log
    log.debug('================================================')
    #matrix_name=b'B0', junk1=0, matrix_shape=6, tin=2, tout=2, is_phase=0, junk2=0, ncols_gset=0

    #ints    = (3, 1, 4, 1, 1717986918, 1072064102, -1, -1, -1, -1,
    #           540029250, 538976288, 0, 6, 1, 1, 0, 0, 0, 10402, 1, 10402, 3, -1717986918, 1071225241, -1, -1, -1, -1,
    #           540029762, 538976288, 0, 6, 1, 1, 0, 0, 0, 30301, 1, 30301, 3, 0, 1071644672, -1, -1, -1, -1)
    #floats  = (3, 1, 4, 1, 2.7200830220753216e+23, 1.7999999523162842, nan, nan, nan, nan,
    #           1.4924077914682395e-19, 1.3563156426940112e-19, 0.0, 6, 1, 1, 0.0, 0.0, 0.0, 10402, 1, 10402, 3, -1.5881868392106856e-23, 1.6999999284744263, nan, nan, nan, nan,
    #           1.4924739659172437e-19, 1.3563156426940112e-19, 0.0, 6, 1, 1, 0.0, 0.0, 0.0, 30301, 1, 30301, 3, 0.0, 1.75, nan, nan, nan, nan)
    #doubles (float64) = (2.1219957924e-314, 2.121995793e-314, 0.7, nan, nan, 6.013470018394104e-154, 1.2731974746e-313, 2.1219957915e-314, 0.0, 2.2073000217621e-310, 2.20730002176213e-310, -2.3534379293677296e-185, nan, nan, 1.208269179683613e-153, 2.66289668e-315, 2.121995794e-314, 5e-324, 0.0, 2.1220107616e-314, 6.3660023436e-314, 0.5, nan, nan)

    #size = op2.size
    factor = op2.factor
    endian = op2._endian
    header_fmt = endian + b'8s 7i'
    #idtype = 'int32'

    #nheader = 48 * factor
    #header = unpack(header_fmt, data[:nheader]) # 48=4*12
    #assert header[:3] == (114, 1, 120), 'header[:3]=%s header=%s' % (header[:3], header)

    # 1 NAME(2) CHAR4
    # 3 UNDEF none
    # 4 MATFORM  I Matrix Form
    # 5 MATRIX T I
    # MATRIX T=1
    # 6 MATTYPE I Matrix Type, repeat of previous word
    # 7 UNDEF(2 ) none
    # 9 MATCOLS I Matrix Columns
    # 10 GJ I
    # 11 CJ I
    # 12 GI I
    # 13 CI I
    # 14 VRS RS
    #matrix_name, junk1, matrix_shape, tin, tout, is_phase, junk2, ncols_gset = header[3:]
    #if factor == 2:
        #matrix_name = reshape_bytes_block(matrix_name)
    #matrix_name = matrix_name.strip()
    #self.log.debug('matrix_name=%s, junk1=%s, matrix_shape=%s, tin=%s, tout=%s, '
                   #'is_phase=%s, junk2=%s, ncols_gset=%s' % (
                       #matrix_name, junk1, matrix_shape, tin, tout,
                       #is_phase, junk2, ncols_gset))
    #self.show_data(data[48*factor:], types='ifsd')

    #is_complex = False
    #if tin > 2 or tout > 2:
        #is_complex = True
        #assert is_phase == 0, 'is_phase=%s' % is_phase
        #imags = []

    #if tout == 1:
        #dtype = 'float32'
        #fdtype = op2.fdtype
    #elif tout == 2:
        #dtype = 'float64'
        #fdtype = op2.double_dtype
    #elif tout == 3:
        #dtype = 'complex64'
        #fdtype = op2.fdtype
    #elif tout == 4:
        #dtype = 'complex128'
        #fdtype = op2.double_dtype
    #else:
        #dtype = '???'
        #msg = ('matrix_name=%s, junk1=%s, matrix_shape=%s, tin=%s, tout=%s, '
               #'is_phase=%s, junk2=%s, ncols_gset=%s' % (
                   #matrix_name, junk1, matrix_shape, tin, tout,
                   #is_phase, junk2, ncols_gset))
        #self.log.warning(msg)
        #raise RuntimeError(msg)

    #is_symmetric = matrix_shape == 6
    #is_phase_flag = is_phase > 0

    n = 12 * factor
    datai = data[n:]
    ints = np.frombuffer(datai, dtype=op2.idtype8).copy()
    #print(ints.tolist())

    # find the first index with ()-1,-1)
    iminus1 = np.where(ints == -1)[0]
    #istart, istop, kstops = _find_dmig_start_stop(datai, iminus1, header_fmt, self.size, self.log)
    size = 4
    istarts, istops, outs, kstarts, kstops = find_all_dmigs_start_stop(
        datai, header_fmt, size, iminus1,
        debug=debug)
    del size

    ioffset = 9
    # istop : the end of the matrix(s)
    # istart : the start of the matrix(s)

    # kstart : the start of the columns
    # kstop  : the end of the columns
    assert len(kstarts) > 0, kstarts
    log.info(f'  outs = {outs}')
    log.info(f'  istarts = {istarts}')
    log.info(f'  istops  = {istops}')
    log.info(f'  kstarts = {kstarts}')
    log.info(f'  kstops  = {kstops}')
    nmatrices = len(outs)
    for i, istart, istop, out, kstart, kstop in zip(count(), istarts, istops, outs, kstarts, kstops):
        # istop : the end of the matrix(s)
        # istart : the start of the matrix(s)
        #
        # kstart : the start of the columns
        # kstop  : the end of the columns
        log.info(f'-------------------------------')
        log.debug(f'n={n} istart={istart} istopi={istop} kstart={kstart} kstop={kstop}')
        kstart0 = kstart[0]
        kstop0 = kstop[0]

        #datai = data[n+istart*size:n+istop*size]
        #op2.show_data(datai, types='ifsd')

        datak = data[n+kstart0*4:n+kstop0*4]
        #op2.show_data(datak, types='ifqd')

        #datak = data[n+kstart0*size+4:n+kstop0*size+4]
        #op2.show_data(datak, types='ifqd')

        matrix_name, junk1, matrix_shape, tin, tout, is_phase, junk2, ncols_gset = out # [3:]
        matrix_name = matrix_name.strip()
        matrix_name_str = matrix_name.decode('latin1')
        #self.show_data(data[48*factor:], types='ifsd')

        is_complex = False
        if tin > 2 or tout > 2:
            #is_complex = True
            assert is_phase == 0, f'is_phase={is_phase}'
            #imags = []

        #if size == 8 and tin == 1:
            #tin = 2
        #if size == 8 and tout == 1:
            #tout = 2
        out_dtype, out_fdtype = get_dtype_fdtype_from_tout(op2, tout)
        in_dtype, in_fdtype = get_dtype_fdtype_from_tout(op2, tin)
        dtype, fdtype = in_dtype, in_fdtype

        is_symmetric = matrix_shape == 6
        is_phase_flag = is_phase > 0
        log.info(f'matrix_name={matrix_name} junk1={junk1} matrix_shape={matrix_shape} '
                 f'tin={tin} ({in_dtype} {in_fdtype}) tout={tout} ({out_dtype} {out_fdtype}) \n    '
                 f'is_phase={is_phase_flag} junk2={junk2} ncols_gset={ncols_gset}')

        #if self.size == 4:
            #if tout == 1:  # float32
                #nvalues = 3
            #elif tout == 2:  # float64
                #nvalues = 4
            #elif tout == 3:  # complex64
                #nvalues = 4
            #elif tout == 4:  # complex128
                #nvalues = 5
            #else:
                #raise NotImplementedError(tout)

        #----------------------------------------
        #print(kstart, kstop)
        #kstart = np.array(kstart
        kstop = np.array(kstop)
        kstart = np.hstack([ioffset, kstop[:-1] + 2])
        #self.show_data(data[kstart*4:kstop*4])
        #kfirst = ioffset
        #klast = kstop[-1]
        log.info(f'{matrix_name_str}: kstart = {kstart}; n={len(kstart)}')
        log.info(f'{matrix_name_str}: kstop  = {kstop}')
        #----------------------------------------
        log.info(f'{matrix_name_str}: casting floats')
        #nheader = ioffset * 4
        #nend = istop * self.size
        #floats = get_floats_4(datai[nheader:nend], fdtype, op2, tout)

        assert tout in [1, 2, 3, 4]
        #print(floats.tolist())
        #print(ints[kfirst:klast].tolist())
        #print(floats)
        #----------------------------------------
        log.info(f'getting col (nid,dof); kstart+1={kstart+1}')

        #C:\MSC.Software\simcenter_nastran_2019.2\tpl_post2\elsum15.op2
        # ints    = ('MATK    ', 0,   6, 2, 2, 0,   0,   0,
        #               10000, 3, 10000, 3, 0,   10.9, 10001, 3, 0  , -10.9, -1, -1,
        #               10001, 3, 10001, 3, 0, 10.9, -1, -1, -1, -1,
        #            'MATLOAD ', 0, 9, 1, 1, 0, 0, 0)
        # floats  = ('MATK    ', 0.0, 6, 2, 2, 0.0, 0.0, 0.0,
        #               10000, 3, 10000, 3, 0.0, 10.9, 10001, 3, 0.0, -10.9, nan, nan,
        #               10001, 3, 10001, 3, 0.0, 10.9, nan, nan, nan, nan,
        #             'MATLOAD ', 0.0, 9, 1, 1, 0.0, 0.0, 0.0)

        col_nids_short = ints[kstart]
        assert col_nids_short.min() > 0, col_nids_short.tolist()

        col_dofs_short = ints[kstart+1]
        ucol_dofs = np.unique(col_dofs_short)
        log.debug(f'  col_nids_short = {col_nids_short}; {ints.dtype}')
        log.debug(f'  col_dofs = {ucol_dofs}')
        for udofi in ucol_dofs:
            if udofi not in [0, 1, 2, 3, 4, 5, 6]:
                datak = data[n+kstart0*4:n+kstop0*4]
                op2.show_data(data[n:n+kstart0*4], types='ifs')
                op2.show_data(datak, types='if')
                #else:
                    #op2.show_data(data[n:n+kstart0*4], types='ifsqd')
                    #op2.show_data(datak, types='ifqd')
                print(ints.tolist())
                msg = 'udofi=%s is invalid; must be in [0, 1, 2, 3, 4, 5, 6]; dofs=%s' % (
                    udofi, np.asarray(ucol_dofs, dtype='int32').tolist())
                raise ValueError(msg)

        #nj2 = len(istart)  ## TODO: why is this wrong???
        # -------------------------------------------------
        log.info(f'extracting rows')
        row_nids = []
        row_dofs = []
        col_nids = []
        col_dofs = []
        reals = []
        imags = []
        log.debug(f'  dtype={dtype} fdtype={fdtype}')
        log.debug(f'  istart = {istart}')
        log.debug(f'  istop = {istop}')
        log.debug(f'  kstart = {kstart}')
        log.debug(f'  kstop = {kstop}')

        for col_nidi, col_dofi, istarti, istopi in zip(
                col_nids_short, col_dofs_short, kstart, kstop):

            nstart2 = istarti * 4
            nend_float2 = istopi * 4
            #nend_int2 = istopi * 4

            assert nend_float2 > nstart2
            #log.debug(f'  nstart2={nstart2} nend_float2={nend_float2}')

            #ints2 = np.frombuffer(datai[istarti*4:(istopi-2)*4], dtype=op2.idtype8).copy()
            #print(f'ints test = {ints2}')
            #log.debug(f'  nid={col_nidi} dof={col_dofi} istarti={istarti} istopi={istopi}')
            ## TODO: preallocate arrays
            imag = None
            # The float32/complex64 blocks are really simple
            # because we can just use the data is from the ints block.
            # We calculate istart/istop and directly access the float data.
            #
            # The float64/complex128 blocks are more complicated.
            # It's easier to just use istart/istop to calculate datai
            # case that, and then slice it.
            #
            # In all cases, we use istart/istop to calculate row/col nid/dof
            # because they are always int32.  Only the real/imag types change.
            #print(datai[nstart2:nend2])
            #self.show_data(datai[nstart2:nend2], types='ifd')
            #print('fdtype =', fdtype)
            #nints = len(ints2)
            #ints2 = ints2.reshape(nints//nvalues, nvalues)

            # irow brings us to the (grid, dof)
            #print(size, dtype)
            #  the +2 on nstart2 is for the grid,comp offset

            # nstart2: The start of the (g,c) data
            #          (including the matrix and the code)
            #          thus we add 2 to get to the (g,c)
            #nstart2 += 9 * 4
            assert nend_float2 > nstart2

            #print(nstart2, nend_int2, nend_float2)
            dataii = datai[nstart2:nend_float2]
            ints2 = np.frombuffer(dataii, dtype='int32').copy()
            floats = get_floats_4(dataii, fdtype, op2, tout)
            #op2.show_data(datai, types='ifs', endian=None, force=False)

            nints = len(ints2)

            #ints2 = ints2[:nints//3*3]
            #print('ints2 =', ints2, nints)
            #print('floats =', floats, nints)
            if dtype == 'float32':
                # [int, int, float]  -> get float
                real = floats[4::3]
                #print('real =', real)

                # [int, x, x
                #  int, x, x] -> get int
                irow = np.arange(2, nints, step=3, dtype='int32')
            elif dtype == 'float64':
                # [long, double]  -> get double
                real = floats[1::2]

                # [int, x, x, x
                #  int, x, x, x] -> get int
                irow = np.arange(2, nints, step=4, dtype='int32')
            else:
                raise RuntimeError((size, dtype))

            #log.debug(f'ints2 = {ints2}')
            #log.debug(f'real = {real}')
            #log.debug(f'irow = {irow}')
            if len(irow) == 0:
                print(ints2)
                msg = f'irow={irow} nints={nints}'
                raise RuntimeError(msg)
            assert len(irow) > 0, irow
            assert len(real) > 0, real
            assert len(irow) == len(real), len(irow)

            # the row index; [1, 2, ..., 43]
            row_nid = ints2[irow]

            # the dof; [0, 0, ..., 0.]
            row_dof = ints2[irow + 1]
            #log.debug(f'row nid,dof = =({row_nid}, {row_dof})')
            #log.debug(f'row_nid = {row_nid}')
            #log.debug(f'row_dof = {row_dof}')
            urow_dof = np.unique(row_dof)
            #log.debug(f'urow_dof = {urow_dof}')
            for udofi in urow_dof:
                if udofi not in [0, 1, 2, 3, 4, 5, 6]:
                    msg = 'udofi=%s is invalid; must be in [0, 1, 2, 3, 4, 5, 6]; dofs=%s' % (
                        udofi, np.asarray(urow_dof, dtype='int32').tolist())
                    raise ValueError(msg)
            ni = len(irow)
            #log.debug(f'real = {real}')
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
            #print(real_array)
        ioffset = kstop[-1] + 4

        if is_complex:
            complex_array = np.hstack(imags)
            assert len(real_array) == len(complex_array)
            real_imag_array = real_array + 1.j * complex_array
        else:
            real_imag_array = real_array

        #print('matrix_name_str =', matrix_name_str)
        #print('col_nids_array =', col_nids_array)
        #print('col_dofs_array =', col_dofs_array)
        #print('row_nids_array =', row_nids_array)
        #print('row_dofs_array =', row_dofs_array)
        #print('real_imag_array =', real_imag_array)
        m = _cast_matrix_matpool(
            matrix_name_str, real_imag_array,
            col_nids_array, col_dofs_array,
            row_nids_array, row_dofs_array,
            matrix_shape, dtype, is_symmetric, log,
            apply_symmetry=op2.apply_symmetry)
        str(m)
        op2.matrices[matrix_name_str] = m

        if nmatrices > 1:
            log.warning('breaking...')
        #break

def read_matpool_dmig(op2: OP2, data: bytes,
                      unused_utable_name: str, debug: bool=False):  # pragma: no cover
    """
    ncols_gset is needed for form=9
    list of header values:
      4:5   matrix name
      6     placeholder
      7     matrix shape (1=square, 2 or 9=rectangular, 6=symmetric)
      8     input type flag (1=single, 2=double, 3=complex single,
                             4=complex double)
      9     output type flag (0=precision set by system cell,
                              1=single, 2=double, 3=complex single,
                              4=complex double)
     10     complex flag (0=real/imaginary, >0=magnitude/phase)
     11     placeholder
     12     number of columns in the G set
            (only necessary for matrix shape 9)
    """
    log = op2.log
    log.debug('================================================')
    #matrix_name=b'B0', junk1=0, matrix_shape=6, tin=2, tout=2, is_phase=0, junk2=0, ncols_gset=0

    #ints    = (3, 1, 4, 1, 1717986918, 1072064102, -1, -1, -1, -1,
    #           540029250, 538976288, 0, 6, 1, 1, 0, 0, 0, 10402, 1, 10402, 3, -1717986918, 1071225241, -1, -1, -1, -1,
    #           540029762, 538976288, 0, 6, 1, 1, 0, 0, 0, 30301, 1, 30301, 3, 0, 1071644672, -1, -1, -1, -1)
    #floats  = (3, 1, 4, 1, 2.7200830220753216e+23, 1.7999999523162842, nan, nan, nan, nan,
    #           1.4924077914682395e-19, 1.3563156426940112e-19, 0.0, 6, 1, 1, 0.0, 0.0, 0.0, 10402, 1, 10402, 3, -1.5881868392106856e-23, 1.6999999284744263, nan, nan, nan, nan,
    #           1.4924739659172437e-19, 1.3563156426940112e-19, 0.0, 6, 1, 1, 0.0, 0.0, 0.0, 30301, 1, 30301, 3, 0.0, 1.75, nan, nan, nan, nan)
    #doubles (float64) = (2.1219957924e-314, 2.121995793e-314, 0.7, nan, nan, 6.013470018394104e-154, 1.2731974746e-313, 2.1219957915e-314, 0.0, 2.2073000217621e-310, 2.20730002176213e-310, -2.3534379293677296e-185, nan, nan, 1.208269179683613e-153, 2.66289668e-315, 2.121995794e-314, 5e-324, 0.0, 2.1220107616e-314, 6.3660023436e-314, 0.5, nan, nan)

    size = op2.size
    factor = op2.factor
    endian = op2._endian
    if size == 4:
        header_fmt = endian + b'8s 7i'
        #idtype = 'int32'
    else:
        #idtype = 'int64'
        header_fmt = endian + b'16s 7q'

    #nheader = 48 * factor
    #header = unpack(header_fmt, data[:nheader]) # 48=4*12
    #assert header[:3] == (114, 1, 120), 'header[:3]=%s header=%s' % (header[:3], header)

    # 1 NAME(2) CHAR4
    # 3 UNDEF none
    # 4 MATFORM  I Matrix Form
    # 5 MATRIX T I
    # MATRIX T=1
    # 6 MATTYPE I Matrix Type, repeat of previous word
    # 7 UNDEF(2 ) none
    # 9 MATCOLS I Matrix Columns
    # 10 GJ I
    # 11 CJ I
    # 12 GI I
    # 13 CI I
    # 14 VRS RS
    #matrix_name, junk1, matrix_shape, tin, tout, is_phase, junk2, ncols_gset = header[3:]
    #if factor == 2:
        #matrix_name = reshape_bytes_block(matrix_name)
    #matrix_name = matrix_name.strip()
    #self.log.debug('matrix_name=%s, junk1=%s, matrix_shape=%s, tin=%s, tout=%s, '
                   #'is_phase=%s, junk2=%s, ncols_gset=%s' % (
                       #matrix_name, junk1, matrix_shape, tin, tout,
                       #is_phase, junk2, ncols_gset))
    #self.show_data(data[48*factor:], types='ifsd')

    #is_complex = False
    #if tin > 2 or tout > 2:
        #is_complex = True
        #assert is_phase == 0, 'is_phase=%s' % is_phase
        #imags = []

    #if tout == 1:
        #dtype = 'float32'
        #fdtype = op2.fdtype
    #elif tout == 2:
        #dtype = 'float64'
        #fdtype = op2.double_dtype
    #elif tout == 3:
        #dtype = 'complex64'
        #fdtype = op2.fdtype
    #elif tout == 4:
        #dtype = 'complex128'
        #fdtype = op2.double_dtype
    #else:
        #dtype = '???'
        #msg = ('matrix_name=%s, junk1=%s, matrix_shape=%s, tin=%s, tout=%s, '
               #'is_phase=%s, junk2=%s, ncols_gset=%s' % (
                   #matrix_name, junk1, matrix_shape, tin, tout,
                   #is_phase, junk2, ncols_gset))
        #self.log.warning(msg)
        #raise RuntimeError(msg)

    #is_symmetric = matrix_shape == 6
    #is_phase_flag = is_phase > 0

    n = 12 * factor
    datai = data[n:]
    ints = np.frombuffer(datai, dtype=op2.idtype8).copy()
    #print(ints.tolist())

    # find the first index with ()-1,-1)
    iminus1 = np.where(ints == -1)[0]
    #istart, istop, kstops = _find_dmig_start_stop(datai, iminus1, header_fmt, self.size, self.log)
    istarts, istops, outs, kstarts, kstops = find_all_dmigs_start_stop(
        datai, header_fmt, size, iminus1,
        debug=debug)

    ioffset = 9
    # istop : the end of the matrix(s)
    # istart : the start of the matrix(s)

    # kstart : the start of the columns
    # kstop  : the end of the columns
    assert len(kstarts) > 0, kstarts
    log.info(f'  outs = {outs}')
    log.info(f'  istarts = {istarts}')
    log.info(f'  istops  = {istops}')
    log.info(f'  kstarts = {kstarts}')
    log.info(f'  kstops  = {kstops}')
    nmatrices = len(outs)
    for i, istart, istop, out, kstart, kstop in zip(count(), istarts, istops, outs, kstarts, kstops):
        # istop : the end of the matrix(s)
        # istart : the start of the matrix(s)
        #
        # kstart : the start of the columns
        # kstop  : the end of the columns
        log.info(f'-------------------------------')
        log.debug(f'n={n} istart={istart} istopi={istop} kstart={kstart} kstop={kstop}')
        kstart0 = kstart[0]
        kstop0 = kstop[0]

        #datai = data[n+istart*size:n+istop*size]
        #op2.show_data(datai, types='ifsd')

        datak = data[n+kstart0*size:n+kstop0*size]
        #op2.show_data(datak, types='ifqd')

        #datak = data[n+kstart0*size+4:n+kstop0*size+4]
        #op2.show_data(datak, types='ifqd')

        matrix_name, junk1, matrix_shape, tin, tout, is_phase, junk2, ncols_gset = out # [3:]
        if size == 8:
            matrix_name = reshape_bytes_block(matrix_name)
        matrix_name = matrix_name.strip()
        matrix_name_str = matrix_name.decode('latin1')
        #self.show_data(data[48*factor:], types='ifsd')

        is_complex = False
        if tin > 2 or tout > 2:
            #is_complex = True
            assert is_phase == 0, f'is_phase={is_phase}'
            #imags = []

        #if size == 8 and tin == 1:
            #tin = 2
        #if size == 8 and tout == 1:
            #tout = 2
        out_dtype, out_fdtype = get_dtype_fdtype_from_tout(op2, tout)
        in_dtype, in_fdtype = get_dtype_fdtype_from_tout(op2, tin)
        dtype, fdtype = in_dtype, in_fdtype

        is_symmetric = matrix_shape == 6
        is_phase_flag = is_phase > 0
        log.info(f'matrix_name={matrix_name} junk1={junk1} matrix_shape={matrix_shape} '
                 f'tin={tin} ({in_dtype} {in_fdtype}) tout={tout} ({out_dtype} {out_fdtype}) \n    '
                 f'is_phase={is_phase_flag} junk2={junk2} ncols_gset={ncols_gset}')

        #if self.size == 4:
            #if tout == 1:  # float32
                #nvalues = 3
            #elif tout == 2:  # float64
                #nvalues = 4
            #elif tout == 3:  # complex64
                #nvalues = 4
            #elif tout == 4:  # complex128
                #nvalues = 5
            #else:
                #raise NotImplementedError(tout)

        #----------------------------------------
        #print(kstart, kstop)
        #kstart = np.array(kstart
        kstop = np.array(kstop)
        kstart = np.hstack([ioffset, kstop[:-1] + 2])
        #self.show_data(data[kstart*4:kstop*4])
        #kfirst = ioffset
        #klast = kstop[-1]
        log.info(f'kstart = {kstart}')
        log.info(f'kstop  = {kstop}')
        #----------------------------------------
        log.info('casting floats')
        nheader = ioffset * size
        #nend = istop * self.size
        #floats = get_floats(datai[nheader:nend], fdtype, op2, tout, self.size)

        assert tout in [1, 2, 3, 4]
        #print(floats.tolist())
        #print(ints[kfirst:klast].tolist())
        #print(floats)
        #----------------------------------------
        log.info('getting col (nid,dof)')

        col_nids_short = ints[kstart]
        assert col_nids_short.min() > 0, col_nids_short.tolist()

        col_dofs_short = ints[kstart+1]
        ucol_dofs = np.unique(col_dofs_short)
        log.debug(f'  col_nids_short = {col_nids_short}; {ints.dtype}')
        log.debug(f'  col_dofs = {ucol_dofs}')
        for udofi in ucol_dofs:
            if udofi not in [0, 1, 2, 3, 4, 5, 6]:
                datak = data[n+kstart0*size:n+kstop0*size]
                if size == 4:
                    op2.show_data(data[n:n+kstart0*size], types='ifs')
                    op2.show_data(datak, types='if')
                else:
                    op2.show_data(data[n:n+kstart0*size], types='ifsqd')
                    op2.show_data(datak, types='ifqd')
                print(ints.tolist())
                msg = 'udofi=%s is invalid; must be in [0, 1, 2, 3, 4, 5, 6]; dofs=%s' % (
                    udofi, np.asarray(ucol_dofs, dtype='int32').tolist())
                raise ValueError(msg)

        #nj2 = len(istart)  ## TODO: why is this wrong???
        # -------------------------------------------------
        log.info(f'extracting rows')
        row_nids: list[np.ndarray] = []
        row_dofs: list[np.ndarray] = []
        col_nids: list[np.ndarray] = []
        col_dofs: list[np.ndarray] = []
        reals: list[np.ndarray] = []
        imags: list[np.ndarray] = []
        log.debug(f'  dtype={dtype} fdtype={fdtype}')
        log.debug(f'  istart = {istart}')
        log.debug(f'  istop = {istop}')
        log.debug(f'  kstart = {kstart}')
        log.debug(f'  kstop = {kstop}')

        for col_nidi, col_dofi, istarti, istopi in zip(
                col_nids_short, col_dofs_short, kstart, kstop):

            if size == 4:
                nstart2 = istarti * size
                nend_float2 = istopi * size
                #nend_int2 = istopi * size
            else:
                nstart2 = istarti * size
                nend_float2 = istopi * size
                #nend_int2 = istopi * size
                data
                #self.show_data(data[n:], types='ifsdq')

            assert nend_float2 > nstart2
            log.debug(f'  nstart2={nstart2} nend_float2={nend_float2}')

            #ints2 = np.frombuffer(datai[istarti*4:(istopi-2)*4], dtype=op2.idtype8).copy()
            #print(f'ints test = {ints2}')
            log.debug(f'  nid={col_nidi} dof={col_dofi} istarti={istarti} istopi={istopi}')
            ## TODO: preallocate arrays
            imag = None
            # The float32/complex64 blocks are really simple
            # because we can just use the data is from the ints block.
            # We calculate istart/istop and directly access the float data.
            #
            # The float64/complex128 blocks are more complicated.
            # It's easier to just use istart/istop to calculate datai
            # case that, and then slice it.
            #
            # In all cases, we use istart/istop to calculate row/col nid/dof
            # because they are always int32.  Only the real/imag types change.
            #print(datai[nstart2:nend2])
            #self.show_data(datai[nstart2:nend2], types='ifd')
            #print('fdtype =', fdtype)
            #nints = len(ints2)
            #ints2 = ints2.reshape(nints//nvalues, nvalues)

            # irow brings us to the (grid, dof)
            #print(size, dtype)
            #  the +2 on nstart2 is for the grid,comp offset

            # nstart2: The start of the (g,c) data
            #          (including the matrix and the code)
            #          thus we add 2 to get to the (g,c)
            #nstart2 += 9 * 4
            assert nend_float2 > nstart2
            if size == 4:
                #print(nstart2, nend_int2, nend_float2)
                dataii = datai[nstart2:nend_float2]
                ints2 = np.frombuffer(dataii, dtype='int32').copy()
                floats = get_floats_4(dataii, fdtype, op2, tout)
                #op2.show_data(datai, types='ifs', endian=None, force=False)

                nints = len(ints2)

                #ints2 = ints2[:nints//3*3]
                #print('ints2 =', ints2, nints)
                #print('floats =', floats, nints)
                if dtype == 'float32':
                    # [int, int, float]  -> get float
                    real = floats[4::3]
                    #print('real =', real)

                    # [int, x, x
                    #  int, x, x] -> get int
                    irow = np.arange(2, nints, step=3, dtype='int32')
                elif dtype == 'float64':
                    # [long, double]  -> get double
                    real = floats[1::2]

                    # [int, x, x, x
                    #  int, x, x, x] -> get int
                    irow = np.arange(2, nints, step=4, dtype='int32')
                else:
                    raise RuntimeError((size, dtype))
            else:
                assert nend_float2 - nstart2 > 1, nend_float2 - nstart2
                dataii = datai[nstart2:nend_float2]
                op2.show_data(datai, types='ifsqd', endian=None, force=False)
                ints2 = np.frombuffer(dataii, dtype='int64').copy()
                floats = np.frombuffer(dataii, dtype='float32').copy()
                #doubles = np.frombuffer(dataii, dtype='float64').copy()
                #print('floats', floats)
                #print('doubles', doubles)
                nints = len(ints2)
                if dtype == 'float32':
                    #0   1   2     3   4    5
                    # g1, c1,  (gi, ci, ri), ...
                    # (g, c,  # 0, 1
                    #  g, c, r, # 2, 3, 4
                    #  )
                    irow = np.arange(2, nints, step=3, dtype='int64')
                    #real = floats[irow + 4]  # nope...
                    #real = floats[irow + 2]  # ???
                    real = floats[3::3]
                    #print(real_old, real)
                    print(real, real.dtype)
                    #sss
                #elif dtype == 'float64':
                    #print(ints)
                    #irow = np.arange(2, nints, step=4, dtype='int64')
                    #print(floats, floats.dtype)
                    #real = floats[2::3]
                    #print(real)
                else:
                    raise RuntimeError((size, dtype))
            log.debug(f'ints2 = {ints2}')
            log.debug(f'real = {real}')
            log.debug(f'irow = {irow}')
            if len(irow) == 0:
                print(ints2)
                msg = f'irow={irow} nints={nints}'
                raise RuntimeError(msg)
            assert len(irow) > 0, irow
            assert len(real) > 0, real
            assert len(irow) == len(real), len(irow)

            if 0:  # pragma: no cover
                if dtype == 'float32':
                    irow = np.arange(istarti, istopi-1, step=3, dtype='int32')
                    real = floats[irow + 2]
                elif dtype == 'complex64':
                    irow = np.arange(istarti, istopi-1, step=4, dtype='int32')
                    real = floats[irow + 2]
                    imag = floats[irow + 3]

                elif dtype == 'float64':
                    datai = data[nheader+(istarti*4) : nheader+(istopi*4)]
                    #self.show_data(datai)
                    irow = np.arange(istarti, istopi-1, step=4, dtype='int32')
                    log.debug(f'irow float64: irow={irow}')
                    assert len(datai) % 8 == 0, len(datai) / 8
                    real = np.frombuffer(datai, dtype=fdtype)[1::2].copy()

                elif dtype == 'complex128':
                    datai = data[nheader+(istarti*4) : nheader+(istopi*4)]

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
                    msg = f'{dtype} is not supported'
                    log.error(msg)
                    raise RuntimeError(msg)

                if len(irow) != len(real):
                    msg = 'nrow=%s nreal=%s' % (len(irow), len(real))
                    raise RuntimeError(msg)
                #elif len(irow) == 0:
                    #msg = 'nrow=%s nreal=%s' % (len(irow), len(real))
                    #raise RuntimeError(msg)

            # the row index; [1, 2, ..., 43]
            row_nid = ints2[irow]

            # the dof; [0, 0, ..., 0.]
            row_dof = ints2[irow + 1]
            log.debug(f'row nid,dof = =({row_nid}, {row_dof})')
            log.debug(f'row_nid = {row_nid}')
            log.debug(f'row_dof = {row_dof}')
            urow_dof = np.unique(row_dof)
            log.debug(f'urow_dof = {urow_dof}')
            for udofi in urow_dof:
                if udofi not in [0, 1, 2, 3, 4, 5, 6]:
                    msg = 'udofi=%s is invalid; must be in [0, 1, 2, 3, 4, 5, 6]; dofs=%s' % (
                        udofi, np.asarray(urow_dof, dtype='int32').tolist())
                    raise ValueError(msg)
            ni = len(irow)
            log.debug(f'real = {real}')
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
            print(real_array)
        ioffset = kstop[-1] + 4

        if is_complex:
            complex_array = np.hstack(imags)
            assert len(real_array) == len(complex_array)
            real_imag_array = real_array + 1.j * complex_array
        else:
            real_imag_array = real_array

        #print('matrix_name_str =', matrix_name_str)
        #print('col_nids_array =', col_nids_array)
        #print('col_dofs_array =', col_dofs_array)
        #print('row_nids_array =', row_nids_array)
        #print('row_dofs_array =', row_dofs_array)
        #print('real_imag_array =', real_imag_array)
        m = _cast_matrix_matpool(
            matrix_name_str, real_imag_array,
            col_nids_array, col_dofs_array,
            row_nids_array, row_dofs_array,
            matrix_shape, dtype, is_symmetric, log,
            apply_symmetry=op2.apply_symmetry)
        str(m)
        op2.matrices[matrix_name_str] = m

        if nmatrices > 1:
            log.warning('breaking...')
        #break
    return

def read_matpool_dmig_8(op2: OP2, data: bytes,
                        unused_utable_name: str, debug: bool=False):
    """
    ncols_gset is needed for form=9
    list of header values:
      4:5   matrix name
      6     placeholder
      7     matrix shape (1=square, 2 or 9=rectangular, 6=symmetric)
      8     input type flag (1=single, 2=double, 3=complex single,
                             4=complex double)
      9     output type flag (0=precision set by system cell,
                              1=single, 2=double, 3=complex single,
                              4=complex double)
     10     complex flag (0=real/imaginary, >0=magnitude/phase)
     11     placeholder
     12     number of columns in the G set
            (only necessary for matrix shape 9)
    """
    #aaa
    log = op2.log
    log.debug('================================================')
    #matrix_name=b'B0', junk1=0, matrix_shape=6, tin=2, tout=2, is_phase=0, junk2=0, ncols_gset=0

    #ints    = (3, 1, 4, 1, 1717986918, 1072064102, -1, -1, -1, -1,
    #           540029250, 538976288, 0, 6, 1, 1, 0, 0, 0, 10402, 1, 10402, 3, -1717986918, 1071225241, -1, -1, -1, -1,
    #           540029762, 538976288, 0, 6, 1, 1, 0, 0, 0, 30301, 1, 30301, 3, 0, 1071644672, -1, -1, -1, -1)
    #floats  = (3, 1, 4, 1, 2.7200830220753216e+23, 1.7999999523162842, nan, nan, nan, nan,
    #           1.4924077914682395e-19, 1.3563156426940112e-19, 0.0, 6, 1, 1, 0.0, 0.0, 0.0, 10402, 1, 10402, 3, -1.5881868392106856e-23, 1.6999999284744263, nan, nan, nan, nan,
    #           1.4924739659172437e-19, 1.3563156426940112e-19, 0.0, 6, 1, 1, 0.0, 0.0, 0.0, 30301, 1, 30301, 3, 0.0, 1.75, nan, nan, nan, nan)
    #doubles (float64) = (2.1219957924e-314, 2.121995793e-314, 0.7, nan, nan, 6.013470018394104e-154, 1.2731974746e-313, 2.1219957915e-314, 0.0, 2.2073000217621e-310, 2.20730002176213e-310, -2.3534379293677296e-185, nan, nan, 1.208269179683613e-153, 2.66289668e-315, 2.121995794e-314, 5e-324, 0.0, 2.1220107616e-314, 6.3660023436e-314, 0.5, nan, nan)

    factor = op2.factor
    endian = op2._endian
    #idtype = 'int64'
    header_fmt = endian + b'16s 7q'

    #nheader = 48 * factor
    #header = unpack(header_fmt, data[:nheader]) # 48=4*12
    #assert header[:3] == (114, 1, 120), 'header[:3]=%s header=%s' % (header[:3], header)

    # 1 NAME(2) CHAR4
    # 3 UNDEF none
    # 4 MATFORM  I Matrix Form
    # 5 MATRIX T I
    # MATRIX T=1
    # 6 MATTYPE I Matrix Type, repeat of previous word
    # 7 UNDEF(2 ) none
    # 9 MATCOLS I Matrix Columns
    # 10 GJ I
    # 11 CJ I
    # 12 GI I
    # 13 CI I
    # 14 VRS RS
    #matrix_name, junk1, matrix_shape, tin, tout, is_phase, junk2, ncols_gset = header[3:]
    #if factor == 2:
        #matrix_name = reshape_bytes_block(matrix_name)
    #matrix_name = matrix_name.strip()
    #self.log.debug('matrix_name=%s, junk1=%s, matrix_shape=%s, tin=%s, tout=%s, '
                   #'is_phase=%s, junk2=%s, ncols_gset=%s' % (
                       #matrix_name, junk1, matrix_shape, tin, tout,
                       #is_phase, junk2, ncols_gset))
    #self.show_data(data[48*factor:], types='ifsd')

    #is_complex = False
    #if tin > 2 or tout > 2:
        #is_complex = True
        #assert is_phase == 0, 'is_phase=%s' % is_phase
        #imags = []

    #if tout == 1:
        #dtype = 'float32'
        #fdtype = op2.fdtype
    #elif tout == 2:
        #dtype = 'float64'
        #fdtype = op2.double_dtype
    #elif tout == 3:
        #dtype = 'complex64'
        #fdtype = op2.fdtype
    #elif tout == 4:
        #dtype = 'complex128'
        #fdtype = op2.double_dtype
    #else:
        #dtype = '???'
        #msg = ('matrix_name=%s, junk1=%s, matrix_shape=%s, tin=%s, tout=%s, '
               #'is_phase=%s, junk2=%s, ncols_gset=%s' % (
                   #matrix_name, junk1, matrix_shape, tin, tout,
                   #is_phase, junk2, ncols_gset))
        #self.log.warning(msg)
        #raise RuntimeError(msg)

    #is_symmetric = matrix_shape == 6
    #is_phase_flag = is_phase > 0

    n = 12 * factor
    datai = data[n:]
    ints = np.frombuffer(datai, dtype=op2.idtype8).copy()
    #print(ints.tolist())

    # find the first index with ()-1,-1)
    iminus1 = np.where(ints == -1)[0]
    #istart, istop, kstops = _find_dmig_start_stop(datai, iminus1, header_fmt, self.size, self.log)
    size = 8
    istarts, istops, outs, kstarts, kstops = find_all_dmigs_start_stop(
        datai, header_fmt, size, iminus1,
        debug=debug)

    ioffset = 9
    # istop : the end of the matrix(s)
    # istart : the start of the matrix(s)

    # kstart : the start of the columns
    # kstop  : the end of the columns
    assert len(kstarts) > 0, kstarts
    log.info(f'  outs = {outs}')
    log.info(f'  istarts = {istarts}')
    log.info(f'  istops  = {istops}')
    log.info(f'  kstarts = {kstarts}')
    log.info(f'  kstops  = {kstops}')
    nmatrices = len(outs)
    for i, istart, istop, out, kstart, kstop in zip(count(), istarts, istops, outs, kstarts, kstops):
        # istop : the end of the matrix(s)
        # istart : the start of the matrix(s)
        #
        # kstart : the start of the columns
        # kstop  : the end of the columns
        log.info(f'-------------------------------')
        log.debug(f'n={n} istart={istart} istopi={istop} kstart={kstart} kstop={kstop}')
        kstart0 = kstart[0]
        kstop0 = kstop[0]

        #datai = data[n+istart*size:n+istop*size]
        #op2.show_data(datai, types='ifsd')

        datak = data[n+kstart0*size:n+kstop0*size]
        #op2.show_data(datak, types='ifqd')

        #datak = data[n+kstart0*size+4:n+kstop0*size+4]
        #op2.show_data(datak, types='ifqd')

        matrix_name, junk1, matrix_shape, tin, tout, is_phase, junk2, ncols_gset = out # [3:]
        matrix_name = reshape_bytes_block(matrix_name).strip()
        matrix_name_str = matrix_name.decode('latin1')
        #self.show_data(data[48*factor:], types='ifsd')

        is_complex = False
        if tin > 2 or tout > 2:
            #is_complex = True
            assert is_phase == 0, f'is_phase={is_phase}'
            #imags = []

        #if size == 8 and tin == 1:
            #tin = 2
        #if size == 8 and tout == 1:
            #tout = 2
        out_dtype, out_fdtype = get_dtype_fdtype_from_tout(op2, tout)
        in_dtype, in_fdtype = get_dtype_fdtype_from_tout(op2, tin)
        dtype, fdtype = in_dtype, in_fdtype

        is_symmetric = matrix_shape == 6
        is_phase_flag = is_phase > 0
        log.info(f'matrix_name={matrix_name} junk1={junk1} matrix_shape={matrix_shape} '
                 f'tin={tin} ({in_dtype} {in_fdtype}) tout={tout} ({out_dtype} {out_fdtype}) \n    '
                 f'is_phase={is_phase_flag} junk2={junk2} ncols_gset={ncols_gset}')

        #if self.size == 4:
            #if tout == 1:  # float32
                #nvalues = 3
            #elif tout == 2:  # float64
                #nvalues = 4
            #elif tout == 3:  # complex64
                #nvalues = 4
            #elif tout == 4:  # complex128
                #nvalues = 5
            #else:
                #raise NotImplementedError(tout)

        #----------------------------------------
        #print(kstart, kstop)
        #kstart = np.array(kstart
        kstop = np.array(kstop)
        kstart = np.hstack([ioffset, kstop[:-1] + 2])
        #self.show_data(data[kstart*4:kstop*4])
        #kfirst = ioffset
        #klast = kstop[-1]
        log.info(f'kstart = {kstart}')
        log.info(f'kstop  = {kstop}')
        #----------------------------------------
        log.info('casting floats')
        #nheader = ioffset * size
        #nend = istop * self.size
        #floats = get_floats(datai[nheader:nend], fdtype, op2, tout, self.size)

        assert tout in [1, 2, 3, 4]
        #print(floats.tolist())
        #print(ints[kfirst:klast].tolist())
        #print(floats)
        #----------------------------------------
        log.info('getting col (nid,dof)')

        col_nids_short = ints[kstart]
        assert col_nids_short.min() > 0, col_nids_short.tolist()

        col_dofs_short = ints[kstart+1]
        ucol_dofs = np.unique(col_dofs_short)
        log.debug(f'  col_nids_short = {col_nids_short}; {ints.dtype}')
        log.debug(f'  col_dofs = {ucol_dofs}')
        for udofi in ucol_dofs:
            if udofi not in [0, 1, 2, 3, 4, 5, 6]:
                datak = data[n+kstart0*size:n+kstop0*size]
                if size == 4:
                    op2.show_data(data[n:n+kstart0*size], types='ifs')
                    op2.show_data(datak, types='if')
                else:
                    op2.show_data(data[n:n+kstart0*size], types='ifsqd')
                    op2.show_data(datak, types='ifqd')
                print(ints.tolist())
                msg = 'udofi=%s is invalid; must be in [0, 1, 2, 3, 4, 5, 6]; dofs=%s' % (
                    udofi, np.asarray(ucol_dofs, dtype='int32').tolist())
                raise ValueError(msg)

        #nj2 = len(istart)  ## TODO: why is this wrong???
        # -------------------------------------------------
        log.info(f'extracting rows')
        row_nids = []
        row_dofs = []
        col_nids = []
        col_dofs = []
        reals = []
        imags = []
        log.debug(f'  dtype={dtype} fdtype={fdtype}')
        log.debug(f'  istart = {istart}')
        log.debug(f'  istop = {istop}')
        log.debug(f'  kstart = {kstart}')
        log.debug(f'  kstop = {kstop}')

        for col_nidi, col_dofi, istarti, istopi in zip(
                col_nids_short, col_dofs_short, kstart, kstop):

            nstart2 = istarti * size
            nend_float2 = istopi * size
            #nend_int2 = istopi * size
            data
            #self.show_data(data[n:], types='ifsdq')

            assert nend_float2 > nstart2
            log.debug(f'  nstart2={nstart2} nend_float2={nend_float2}')

            #ints2 = np.frombuffer(datai[istarti*4:(istopi-2)*4], dtype=op2.idtype8).copy()
            #print(f'ints test = {ints2}')
            log.debug(f'  nid={col_nidi} dof={col_dofi} istarti={istarti} istopi={istopi}')
            ## TODO: preallocate arrays
            imag = None
            # The float32/complex64 blocks are really simple
            # because we can just use the data is from the ints block.
            # We calculate istart/istop and directly access the float data.
            #
            # The float64/complex128 blocks are more complicated.
            # It's easier to just use istart/istop to calculate datai
            # case that, and then slice it.
            #
            # In all cases, we use istart/istop to calculate row/col nid/dof
            # because they are always int32.  Only the real/imag types change.
            #print(datai[nstart2:nend2])
            #self.show_data(datai[nstart2:nend2], types='ifd')
            #print('fdtype =', fdtype)
            #nints = len(ints2)
            #ints2 = ints2.reshape(nints//nvalues, nvalues)

            # irow brings us to the (grid, dof)
            #print(size, dtype)
            #  the +2 on nstart2 is for the grid,comp offset

            # nstart2: The start of the (g,c) data
            #          (including the matrix and the code)
            #          thus we add 2 to get to the (g,c)
            #nstart2 += 9 * 4
            assert nend_float2 > nstart2
            assert nend_float2 - nstart2 > 1, nend_float2 - nstart2
            dataii = datai[nstart2:nend_float2]

            log.debug('showing data')
            op2.show_data(datai[:16], types='s', endian=None, force=False)
            op2.show_data(datai[16:], types='ifqd', endian=None, force=False)

            ints2 = np.frombuffer(dataii, dtype='int64').copy()
            floats = np.frombuffer(dataii, dtype='float64').copy()
            #doubles = np.frombuffer(dataii, dtype='float64').copy()
            #print('floats', floats)
            #print('doubles', doubles)
            nints = len(ints2)
            log.debug(f'  dtype={dtype}')
            if dtype == 'float32':  # it's a float64 for 64 bit :(
                #0   1   2     3   4    5
                # g1, c1,  (gi, ci, ri), ...
                # (g, c,  # 0, 1
                #  g, c, r, # 2, 3, 4
                #  )
                irow = np.arange(2, nints, step=3, dtype='int64')
                #real = floats[irow + 4]  # nope...
                #real = floats[irow + 2]  # ???
                real = floats[3::3]
                #print(real_old, real)
                print(real, real.dtype)
                #sss
            #elif dtype == 'float64':
                #print(ints)
                #irow = np.arange(2, nints, step=4, dtype='int64')
                #print(floats, floats.dtype)
                #real = floats[2::3]
                #print(real)
            else:
                raise RuntimeError((size, dtype))
            log.debug(f'ints2 = {ints2}')
            log.debug(f'real = {real}')
            log.debug(f'irow = {irow}')
            if len(irow) == 0:
                print(ints2)
                msg = f'irow={irow} nints={nints}'
                raise RuntimeError(msg)
            assert len(irow) > 0, irow
            assert len(real) > 0, real
            assert len(irow) == len(real), len(irow)

            # the row index; [1, 2, ..., 43]
            row_nid = ints2[irow]

            # the dof; [0, 0, ..., 0.]
            row_dof = ints2[irow + 1]
            log.debug(f'row nid,dof = =({row_nid}, {row_dof})')
            log.debug(f'row_nid = {row_nid}')
            log.debug(f'row_dof = {row_dof}')
            urow_dof = np.unique(row_dof)
            log.debug(f'urow_dof = {urow_dof}')
            for udofi in urow_dof:
                if udofi not in [0, 1, 2, 3, 4, 5, 6]:
                    msg = 'udofi=%s is invalid; must be in [0, 1, 2, 3, 4, 5, 6]; dofs=%s' % (
                        udofi, np.asarray(urow_dof, dtype='int32').tolist())
                    raise ValueError(msg)
            ni = len(irow)
            log.debug(f'real = {real}')
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
            print(real_array)
        ioffset = kstop[-1] + 4

        if is_complex:
            complex_array = np.hstack(imags)
            assert len(real_array) == len(complex_array)
            real_imag_array = real_array + 1.j * complex_array
        else:
            real_imag_array = real_array

        #print('matrix_name_str =', matrix_name_str)
        #print('col_nids_array =', col_nids_array)
        #print('col_dofs_array =', col_dofs_array)
        #print('row_nids_array =', row_nids_array)
        #print('row_dofs_array =', row_dofs_array)
        #print('real_imag_array =', real_imag_array)
        m = _cast_matrix_matpool(
            matrix_name_str, real_imag_array,
            col_nids_array, col_dofs_array,
            row_nids_array, row_dofs_array,
            matrix_shape, dtype, is_symmetric, log,
            apply_symmetry=op2.apply_symmetry)
        str(m)
        op2.matrices[matrix_name_str] = m

        if nmatrices > 1:
            log.warning('breaking...')
        #break
    return

def _cast_matrix_matpool(table_name: str,
                         real_imag_array,
                         col_nids_array, col_dofs_array,
                         row_nids_array, row_dofs_array,
                         matrix_shape: tuple[int, int],
                         dtype: str, is_symmetric: bool,
                         log: SimpleLogger,
                         apply_symmetry: bool=False) -> Matrix:
    """helper method for _read_matpool_matrix"""
    #op2 = self.op2
    make_matrix_symmetric = apply_symmetry and matrix_shape == 'symmetric'

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
        grids1, comps1, grids2, comps2, make_matrix_symmetric, idtype='int32')
    assert len(j1) == len(j2), 'nj1=%s nj2=%s' % (len(j1), len(j2))
    assert len(grids1) == len(real_imag_array), 'ngrids1=%s nreals=%s' % (len(j1), len(real_imag_array))

    # not 100% on these, they might be flipped
    #ncols = len(np.unique(j1))
    #mrows = len(np.unique(j2))

    if is_symmetric and make_matrix_symmetric:
        mrows = nj
        ncols = nj
        #print('  j1 =', j1)
        #print('  j2 =', j2)
    else:
        ncols = nj1
        mrows = nj2

    # C:\MSC.Software\simcenter_nastran_2019.2\tpl_post2\tr1071x.op2
    #print(real_imag_array)
    #print(j1)
    #print(j2)
    #print(mrows, ncols)
    try:
        matrix = scipy.sparse.coo_matrix(
            (real_imag_array, (j2, j1)),
            shape=(mrows, ncols), dtype=dtype)
    except ValueError:
        print('gc1', grids1, comps1)
        print('gc2', grids2, comps2)
        print(f'j1={j1}')
        print(f'j2={j2}')
        print(f'shape=({mrows}, {ncols})')
        msg = 'Passed all the checks; cannot build MATPOOL sparse matrix...\n'
        spaces = '                                          '
        msg += '%sname=%s dtype=%s nrows=%s ncols=%s nj1=%s nj2=%s nj=%s' % (
            spaces, table_name, dtype, mrows, ncols, nj1, nj2, nj)
        log.error(msg)
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
                log.warning(
                    f'Matrix {table_name!r} marked as symmetric does not contain '
                    'symmetric data.  Data will be symmetrized by averaging.')
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
        #data = matrix

        # matrix is symmetric, but is not stored as symmetric
        matrix_shape = 'rectangular'

    m = Matrix(table_name, form=matrix_shape)
    m.set_matpool_data(matrix,
                       col_nids_array, col_dofs_array,
                       row_nids_array, row_dofs_array)
    #m.form = matrix_shape
    #print(m)
    log.debug(str(m))
    return m


def find_all_dmigs_start_stop(data: bytes, header_fmt: bytes, size: int,
                              iminus1: NDArrayNint, debug: bool=True):
    """
    Parameters
    ----------
    data : bytes

    header_fmt : bytes
        the struct format
    size : int
        4/8
    iminus1: (n,) int ndarray
        flags
    debug: bool; default=True
        turn debugging on

    Returns
    -------
    istop : the end of the matrix(s)
    istart : the start of the matrix(s)

    kstart : the start of the columns
    kstop  : the end of the columns

    [14, 15,
    24, 25,
    37, 38,
    53, 54,
    72, 73,
    94, 95,
    119, 120,
    147, 148,
    178, 179, 212, 213, 249, 250, 289, 290,
    332, 333, 378, 379, 427, 428, 479, 480,
    534, 535, 592, 593, 653, 654, 717, 718,
    784, 785, 854, 855, 927, 928,
    1003, 1004,
    1005, 1006]

    >>> [15, 16, 17, 18,
         34, 35, 36, 37,
         53, 54, 55, 56]

    """
    #debug = True
    #print(iminus1.tolist())
    iends = []
    for i, ivalue in enumerate(iminus1[:-3]):
        if iminus1[i+1] == ivalue + 1 and iminus1[i+2] == ivalue + 2 and iminus1[i+3] == ivalue + 3:
            iends.append(ivalue)
            #iend = i
            #break
    #iends.append(iminus1[-1])
    iminus1 = iminus1[:i+4]
    #if ivalue

    istop = np.array(iends)
    istart = np.hstack([0, istop[:-1] + 4])
    if debug:
        print(f'iends = {iends}')
        print(f'istop = {istop}')
        print(f'istart = {istart}')

    structi = Struct(header_fmt)
    nheader = 9 * size
    outs = []
    kstarts = []
    kstops = []
    ig = 9
    for istarti, istopi in zip(istart, istop):
        ig = istarti + 2
        if debug:  # pragma: no cover
            print(ig)
        datai = data[istarti*size : istarti*size + nheader]
        out = structi.unpack(datai)
        outs.append(out)

        matrix_name, junk1, matrix_shape, tin, tout, is_phase, junk2, ncols_gset = out

        if tout == 1:
            nvalues = 3
        elif tout == 2:
            nvalues = 4
        else:  # pragma: no cover
            raise RuntimeError(tout)
        kstart, kstop = _get_dmig_kstop(ig, nvalues, istopi, iminus1, debug=debug)
        kstarts.append(kstart)
        kstops.append(kstop)
    #kstarts = np.array(kstarts)
    #kstops = np.array(kstops)

    return istart, istop, outs, kstarts, kstops

def _get_dmig_kstop(ig: int, nvalues: int, istop: int, iminus1,
                    debug: bool=True) -> tuple[np.ndarray, np.ndarray]:
    assert isinstance(debug, bool), debug
    kstart = []
    kstop = []
    #ig = 2
    #print('*ig = ', ig)
    kstart = [9 - 2 + ig]
    for i, ivalue in enumerate(iminus1[:-1]):
        if ivalue < ig:
            continue
        if ivalue > istop + 3:
            break

        leftover = (ivalue - ig) % nvalues
        if debug:  # pragma: no cover
            if leftover:
                print('  ', i, ivalue, ig, nvalues, 'leftover')
            else:
                print('  ', i, ivalue, ig, nvalues)
        kstop.append(ivalue)
        ig = ivalue + 4
        #if debug:
            #print('  ig = ', ig)
        kstart.append(ig)
    kstart.pop()
    assert len(kstart) > 0, kstart
    assert len(kstop) > 0, kstop
    kstart_array = np.array(kstart, dtype='int32')
    kstop_array = np.array(kstop, dtype='int32')
    return kstart_array, kstop_array

def get_dtype_fdtype_from_tout(op2: OP2, tout: int) -> tuple[str, str]:
    if tout == 1:
        dtype = 'float32'
        fdtype = op2.fdtype
    elif tout == 2:
        dtype = 'float64'
        fdtype = op2.double_dtype
    elif tout == 3:
        dtype = 'complex64'
        fdtype = op2.fdtype
    elif tout == 4:
        dtype = 'complex128'
        fdtype = op2.double_dtype
    else:  # pragma: no cover
        dtype = '???'
        fdtype = '???'
        #msg = ('matrix_name=%s, junk1=%s, matrix_shape=%s, tin=%s, tout=%s, '
               #'is_phase=%s, junk2=%s, ncols_gset=%s dtype=%s' % (
                   #matrix_name, junk1, matrix_shape, tin, tout,
                   #is_phase, junk2, ncols_gset, dtype))
        #op2.log.warning(msg)
        raise RuntimeError(f'matrix tout={tout}; expected 1=float32, 2=float64, 3=complex64, 4=complex128')
    return dtype, fdtype

def get_ints_4(data: bytes, idtype: str, op2: OP2) -> np.ndarray:
    ints = np.frombuffer(data, dtype=idtype).copy()
    return ints

def get_ints_8(data: bytes, idtype: str, op2: OP2) -> np.ndarray:
    ints = np.frombuffer(data, dtype=op2.idtype8).copy()
    return ints

def get_floats_4(data: bytes, fdtype: str, op2: OP2, tout: int) -> np.ndarray:
    if tout in {1, 3}:
        # works for float32, complex64
        floats = np.frombuffer(data, dtype=fdtype).copy()
    else:
        # works for float64, complex128
        floats = np.frombuffer(data, dtype=fdtype).copy()
    return floats

def get_floats_8(data: bytes, fdtype: str, op2: OP2, tout: int, size: int) -> np.ndarray:
    if tout == 1:
        # C:\MSC.Software\simcenter_nastran_2019.2\tpl_post2\mnfexam32_0.op2
        floats = np.frombuffer(data, dtype=op2.fdtype8).copy()
    return floats

def grids_comp_array_to_indexi(nid_comp_to_dof_index: Any,
                               a_keys: Any, ai, j: int) -> int:
    for nid_dof in ai:
        #nid_dof = (int(nid), int(dof))
        nid_dof = tuple(nid_dof)
        if nid_dof not in a_keys:
            a_keys.add(nid_dof)
            if nid_dof not in nid_comp_to_dof_index:
                nid_comp_to_dof_index[nid_dof] = j
                j += 1
    nj = len(a_keys)
    del a_keys
    assert nj == len(nid_comp_to_dof_index)
    return j

def grids_comp_array_to_index(grids1, comps1, grids2, comps2,
                              make_matrix_symmetric: bool,
                              idtype: str='int32') -> tuple[Any, Any, int, int, int]:
    """
    Maps the dofs

    Returns
    -------
    ja : (nja,) np.ndarray
        The row indicies
    jb : (njb,) np.ndarray
        The column indicies
    nja : int
        The number of rows of the matrix
    njb : int
        The number of columns of the matrix
    nj_sym : int
        Number of
    """
    #from pyNastran.femutils.utils import unique2d
    ai = np.vstack([grids1, comps1]).T
    bi = np.vstack([grids2, comps2]).T
    #print('grids2 =', grids2)
    #print('comps2 =', comps2)
    #c = np.vstack([a, b])
    #assert c.shape[1] == 2, c.shape
    #urows = unique2d(c)
    #urows = c

    #j = 0
    #b_keys = set()
    #for nid_dof in bi:
        #nid_dof = tuple(nid_dof)
        #if nid_dof not in b_keys:
            #b_keys.add(nid_dof)
        #if nid_dof not in nid_comp_to_dof_indexb:
            #nid_comp_to_dof_indexb[nid_dof] = j
            #j += 1
    #njb = len(b_keys)
    #del b_keys

    a_keys: set[str] = set()
    if make_matrix_symmetric:
        nid_comp_to_dof_indexa = {}
        nid_comp_to_dof_indexb = {}
        j = grids_comp_array_to_indexi(nid_comp_to_dof_indexa, a_keys, ai, j=0)
        nj_sym = grids_comp_array_to_indexi(nid_comp_to_dof_indexb, a_keys, bi, j=j)

        j_index = np.zeros(nj_sym, dtype=idtype)
        for i, nid_dof in zip(count(), ai):
            j_index[i] = nid_comp_to_dof_indexa[tuple(nid_dof)]
        for i, nid_dof in zip(count(), bi):
            j_index[i] = nid_comp_to_dof_indexb[tuple(nid_dof)]
        return j_index, j_index, nj_sym, nj_sym, nj_sym
    else:
        nid_comp_to_dof_indexa = {}
        nid_comp_to_dof_indexb = {}
        b_keys: set[str] = set()
        nja = grids_comp_array_to_indexi(nid_comp_to_dof_indexa, a_keys, ai, j=0)
        njb = grids_comp_array_to_indexi(nid_comp_to_dof_indexb, b_keys, bi, j=0)
        nja = len(nid_comp_to_dof_indexa)
        njb = len(nid_comp_to_dof_indexb)
        nj_sym = -1

        ja_index = np.zeros(grids1.shape, dtype=idtype)
        for i, nid_dof in zip(count(), ai.tolist()):
            ja_index[i] = nid_comp_to_dof_indexa[tuple(nid_dof)]

        jb_index = np.zeros(grids2.shape, dtype='int32')
        for i, nid_dof in zip(count(), bi.tolist()):
            jb_index[i] = nid_comp_to_dof_indexb[tuple(nid_dof)]

        return ja_index, jb_index, nja, njb, nj_sym

