from __future__ import annotations
from struct import Struct
from typing import Any, TYPE_CHECKING
import numpy as np

import pyNastran
from pyNastran.op2.op2_interface.op2_reader import mapfmt
from pyNastran.op2.op2_helper import polar_to_real_imag
from pyNastran.op2.op2_interface.utils import (
    # mapfmt, mapfmt8,
    apply_mag_phase,
)

from pyNastran.op2.tables.utils import get_eid_dt_from_eid_device
from pyNastran.op2.op2_interface.op2_codes import TABLES_BYTES

from pyNastran.op2.tables.oes_stressStrain.real.oes_beams import (
    RealBeamStressArray, RealBeamStrainArray,
    # RealNonlinearBeamStressArray
)
from pyNastran.op2.tables.oes_stressStrain.real.oes_beams import (
    RealBeamStressArray, RealBeamStrainArray,
    # RealNonlinearBeamStressArray
)
from pyNastran.op2.tables.oes_stressStrain.complex.oes_beams import ComplexBeamStressArray, ComplexBeamStrainArray
from pyNastran.op2.tables.oes_stressStrain.random.oes_beams import RandomBeamStressArray, RandomBeamStrainArray

if TYPE_CHECKING:  # pragma: no cover
    from pyNastran.op2.op2 import OP2
IS_DEV = 'dev' in pyNastran.__version__


def oes_cbeam_2(op2: OP2, data: bytes, ndata: int,
                dt: int | float, is_magnitude_phase: bool,
                result_type: int,
                prefix: str, postfix: str) -> tuple[int, int, int]:
    """
    reads stress/strain for element type:
     - 2 : CBEAM

    """
    # op2 = self.op2
    n = 0
    ## TODO: fix method to follow correct pattern...regarding???

    stress_strain = 'stress' if op2.is_stress else 'strain'
    result_name = f'{prefix}cbeam_{stress_strain}{postfix}'
    table_name_bytes = op2.table_name
    assert isinstance(table_name_bytes, bytes), table_name_bytes
    assert table_name_bytes in TABLES_BYTES, table_name_bytes

    if op2._results.is_not_saved(result_name):
        return ndata, None, None
    op2._results._found_result(result_name)

    slot = op2.get_result(result_name)
    factor = op2.factor
    if result_type == 0 and op2.num_wide == 111:  # real
        # TODO: vectorize
        ntotal = 444 * op2.factor  # 44 + 10*40  (11 nodes)

        obj_vector_real = RealBeamStressArray if op2.is_stress else RealBeamStrainArray

        nelements = ndata // ntotal
        nlayers = nelements * 11
        auto_return, is_vectorized = op2._create_oes_object4(
            nlayers, result_name, slot, obj_vector_real)
        if auto_return:
            op2._data_factor = 11
            # assert ntotal == op2.num_wide * 4
            return nelements * ntotal, None, None
        obj = op2.obj

        ntotal = op2.num_wide * 4 * op2.factor
        nelements = ndata // ntotal
        if op2.use_vector and is_vectorized and 0:
            raise NotImplementedError('CBEAM-2-real not vectorized')
        else:
            if op2.use_vector and is_vectorized and op2.sort_method == 1:
                obj._times[obj.itime] = dt

                n = nelements * ntotal
                itotal = obj.itotal
                itotal2 = itotal + nelements * 11

                ints = np.frombuffer(data, dtype=op2.idtype8).reshape(nelements, 111)
                floats = np.frombuffer(data, dtype=op2.fdtype8).reshape(nelements, 111)
                # print(ints[:2, :].tolist())
                # CBEAM    6       2       6       8       0.      1.      0.
                # CBEAM    7       2       8       9       0.      1.      0.
                # CBEAM    8       2       9       10      0.      1.      0.
                # CBEAM    9       2       10      11      0.      1.      0.
                # CBEAM    10      2       11      12      0.      1.      0.
                # [[61,
                #      6, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                #      8, 1065353216, 0, 0, 0, 0, 0, 0, 1, 1],
                # [71,
                #      8, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                #      9, 1065353216, 0, 0, 0, 0, 0, 0, 1, 1]]
                eids = ints[:, 0] // 10
                ints2 = ints[:, 1:].reshape(nelements * 11, 10)
                # print('floats[:, 1:].shape', floats[:, 1:].shape)  # (5,110)
                floats2 = floats[:, 1:].reshape(nelements * 11, 10)

                xxb = floats2[:, 1]
                # ints2 = ints[:, :2]
                # print(ints2[0, :])
                nids = ints2[:, 0]
                # print("eids =", eids)
                # print("nids =", nids.tolist())
                # print("xxb =", xxb)

                eids2 = np.array([eids] * 11, dtype=op2.idtype8).T.ravel()
                assert len(eids2) == len(nids)
                obj.element_node[itotal:itotal2, 0] = eids2
                obj.element_node[itotal:itotal2, 1] = nids
                obj.xxb[itotal:itotal2] = xxb
                obj.data[obj.itime, itotal:itotal2, :] = floats2[:, 2:]
                # self.data[self.itime, self.itotal, :] = [sxc, sxd, sxe, sxf,
                # smax, smin, mst, msc]
            else:
                if op2.use_vector:
                    op2.log.debug('vectorize CBEAM real SORT%s' % op2.sort_method)
                n = oes_cbeam_real_111(op2, data,
                                       obj,
                                       nelements, dt)

    elif result_type == 1 and op2.num_wide == 111:  # imag and random?
        # definitely complex results for MSC Nastran 2016.1

        ntotal = 444 * factor  # 44 + 10*40  (11 nodes)
        nelements = ndata // ntotal

        if op2.is_stress:
            obj_vector_complex = ComplexBeamStressArray
        else:
            obj_vector_complex = ComplexBeamStrainArray

        nlayers = nelements * 11
        auto_return, is_vectorized = op2._create_oes_object4(
            nlayers, result_name, slot, obj_vector_complex)
        if auto_return:
            op2._data_factor = 11
            return nelements * ntotal, None, None

        obj = op2.obj

        nnodes = 10  # 11-1
        # ntotal = op2.num_wide * 4
        if op2.use_vector and is_vectorized and op2.sort_method == 1:
            n = nelements * ntotal
            itotal = obj.itotal
            itotal2 = itotal + nelements * 11

            # chop off eid
            floats = np.frombuffer(data, dtype=op2.fdtype8).reshape(nelements, 111)[:, 1:]
            floats2 = floats.reshape(nelements * 11, 10).copy()

            obj._times[obj.itime] = dt
            if obj.itime == 0:
                ints = np.frombuffer(data, dtype=op2.idtype8).reshape(nelements, 111)
                eids = ints[:, 0] // 10
                eids2 = np.array([eids] * 11, dtype=op2.idtype8).T.ravel()

                ints2 = ints[:, 1:].reshape(nelements * 11, 10)

                nids = ints2[:, 0]
                assert eids.min() > 0, eids.min()
                assert nids.min() >= 0, nids.min()
                obj.element_node[itotal:itotal2, 0] = eids2
                obj.element_node[itotal:itotal2, 1] = nids

            #  0    1   2  3  4  5  6   7   8   9
            # grid, sd, c, d, e, f, c2, d2, e2, f2
            real_imag = apply_mag_phase(floats2, is_magnitude_phase, [2, 3, 4, 5], [6, 7, 8, 9])
            obj.data[obj.itime, itotal:itotal2, :] = real_imag
            obj.sd[itotal:itotal2] = floats2[:, 1]

            obj.itotal = itotal2
            # obj.ielement = ielement2
        else:
            if is_vectorized and op2.use_vector:  # pragma: no cover
                op2.log.debug('vectorize CBEAM imag SORT%s' % op2.sort_method)
            n = oes_cbeam_complex_111(op2, data, obj,
                                      nelements, nnodes,
                                      is_magnitude_phase)

    elif result_type == 2 and op2.num_wide == 67:  # random
        # TODO: vectorize
        ntotal = 268 * factor  # 1 + 11*6  (11 nodes)

        obj_vector_random = RandomBeamStressArray if op2.is_stress else RandomBeamStrainArray

        nelements = ndata // ntotal
        nlayers = nelements * 11
        auto_return, is_vectorized = op2._create_oes_object4(
            nlayers, result_name, slot, obj_vector_random)
        if auto_return:
            op2._data_factor = 11
            assert ntotal == op2.num_wide * 4
            return nelements * ntotal, None, None
        obj = op2.obj

        nnodes = 10  # 11-1
        ntotal = op2.num_wide * 4
        nelements = ndata // ntotal
        if op2.use_vector and is_vectorized and 0:
            raise NotImplementedError('CBEAM-2-random not vectorized')
        else:
            if is_vectorized and op2.use_vector:  # pragma: no cover
                op2.log.debug('vectorize CBEAM random SORT%s' % op2.sort_method)
            n = oes_cbeam_random_67(op2, data, obj,
                                    nelements, nnodes, dt)

    elif result_type == 1 and op2.num_wide in [67] and table_name_bytes in [b'OESXNO1']:
        # C:\MSC.Software\simcenter_nastran_2019.2\tpl_post2\tr1081x.op2
        msg = 'skipping freq CBEAM; numwide=67'
        n = op2._not_implemented_or_skip(data, ndata, msg)
        nelements = None
        ntotal = None
    elif result_type == 2 and op2.num_wide in [67] and table_name_bytes in [b'OESXNO1']:
        # C:\MSC.Software\simcenter_nastran_2019.2\tpl_post2\tr1081x.op2
        msg = 'skipping random CBEAM; numwide=67'
        n = op2._not_implemented_or_skip(data, ndata, msg)
        nelements = None
        ntotal = None
    elif result_type == 1 and op2.num_wide in [67] and prefix == 'psd.':
        # [b'OESPSD1'] and op2.table_code == 605:
        # test_nx_beam_psd
        #   device_code = 0   ???
        #   analysis_code = 5 Frequency
        #   table_code = 605
        #   OESPSD1 - OSTPSD2C - Composite Strain Power Spectral Density
        #   format_code = 2 Real / Imaginary
        #   result_type = 1 Complex
        #   sort_method = 1
        #   sort_code = 4
        #   sort_bits = (0, 0, 1)
        #   data_format = 0 Real
        #   sort_type = 0 Sort1
        #   is_random = 1 Random Responses
        # random_code = 0
        ntotal = (67 * 4 * factor)
        nelements = ndata // ntotal
        assert ntotal % ntotal == 0
        assert nelements > 0, nelements

        cls = RandomBeamStressArray if op2.is_stress else RandomBeamStrainArray
        nlayers = nelements * 11
        auto_return, is_vectorized = op2._create_oes_object4(
            nlayers, result_name, slot, cls)

        if auto_return:
            op2._data_factor = 11
            assert ntotal == op2.num_wide * 4
            return nelements * ntotal, None, None
        obj = op2.obj
        nnodes = 10  # 11-1
        n = oes_cbeam_random_67(op2, data, obj,
                                nelements, nnodes, dt)
    else:  # pragma: no cover
        # raise RuntimeError(table_name_bytes)
        if IS_DEV:
            raise RuntimeError(f'prefix={prefix!r}\n{op2.code_information()}')
        msg = 'skipping freq CBEAM; numwide=67'
        n = op2._not_implemented_or_skip(data, ndata, msg)
        nelements = None
        ntotal = None
    return n, nelements, ntotal


def oes_cbeam_real_111(op2: OP2, data: bytes,
                       obj: RealBeamStressArray | RealBeamStrainArray,
                       nelements: int, dt: Any) -> int:
    n = 0
    nnodes = 10  # 11-1
    n1 = 44 * op2.factor
    n2 = 40 * op2.factor
    fmt1 = mapfmt(op2._endian + op2._analysis_code_fmt + b'i9f', op2.size)
    fmt2 = mapfmt(op2._endian + b'i9f', op2.size)
    s1 = Struct(fmt1)
    s2 = Struct(fmt2)
    add_sort_x = getattr(obj, 'add_sort' + str(op2.sort_method))
    add_new_eid_sort_x = getattr(obj, 'add_new_eid_sort' + str(op2.sort_method))
    for unused_i in range(nelements):
        edata = data[n:n+n1]
        n += n1

        out = s1.unpack(edata)
        eid_device = out[0]
        eid, dt = get_eid_dt_from_eid_device(
            eid_device, op2.nonlinear_factor, op2.sort_method)
        if op2.is_debug_file:
            op2.binary_debug.write('CBEAM-2 - eid=%i out=%s\n' % (eid, str(out)))

        #(grid, sd, sxc, sxd, sxe, sxf, smax, smin, mst, msc) = out
        add_new_eid_sort_x(dt, eid, *out[1:])

        for unused_inode in range(nnodes):
            edata = data[n:n+n2]
            n += n2
            out = s2.unpack(edata)
            # (grid, sd, sxc, sxd, sxe, sxf, smax, smin, mst, msc) = out
            add_sort_x(dt, eid, *out)
    return n


def oes_cbeam_complex_111(op2: OP2, data: bytes,
                          obj: ComplexBeamStressArray | ComplexBeamStrainArray,
                          nelements: int, nnodes: int,
                          is_magnitude_phase: bool) -> int:
    n = 0
    #itotal = obj.itotal
    n1 = 44 * op2.factor
    n2 = 40 * op2.factor

    s1 = Struct(op2._endian + mapfmt(op2._analysis_code_fmt + b'i9f', op2.size))
    s2 = Struct(mapfmt(op2._endian + b'i9f', op2.size))

    add_sort_x = getattr(obj, 'add_sort' + str(op2.sort_method))
    for unused_i in range(nelements):
        edata = data[n:n+n1]
        n += n1
        out1 = s1.unpack(edata)
        eid_device = out1[0]
        eid, dt = get_eid_dt_from_eid_device(
            eid_device, op2.nonlinear_factor, op2.sort_method)
        if op2.is_debug_file:
            op2.binary_debug.write('CBEAM-2 - eid=%i out1=%s\n' % (eid, str(out1)))

        (grid, sd,
         excr, exdr, exer, exfr,
         exci, exdi, exei, exfi) = out1[1:]

        if is_magnitude_phase:
            exc = polar_to_real_imag(excr, exci)
            exd = polar_to_real_imag(exdr, exdi)
            exe = polar_to_real_imag(exer, exei)
            exf = polar_to_real_imag(exfr, exfi)
        else:
            exc = complex(excr, exci)
            exd = complex(exdr, exdi)
            exe = complex(exer, exei)
            exf = complex(exfr, exfi)

        add_sort_x(dt, eid, grid, sd,
                   exc, exd, exe, exf)

        for unused_inode in range(nnodes):
            edata = data[n:n+n2]
            n += n2
            out2 = s2.unpack(edata)
            (grid, sd,
             excr, exdr, exer, exfr,
             exci, exdi, exei, exfi) = out2

            if is_magnitude_phase:
                exc = polar_to_real_imag(excr, exci)
                exd = polar_to_real_imag(exdr, exdi)
                exe = polar_to_real_imag(exer, exei)
                exf = polar_to_real_imag(exfr, exfi)
            else:
                exc = complex(excr, exci)
                exd = complex(exdr, exdi)
                exe = complex(exer, exei)
                exf = complex(exfr, exfi)

            add_sort_x(dt, eid, grid, sd,
                       exc, exd, exe, exf)
            if op2.is_debug_file:
                op2.binary_debug.write('CBEAM-2 - eid=%i out2=%s\n' % (eid, str(out2)))
    return n


def oes_cbeam_random_67(op2: OP2, data: bytes,
                        obj: RandomBeamStressArray | RandomBeamStrainArray,
                        nelements: int, nnodes: int, dt: Any) -> int:
    n = 0
    n1 = 28 * op2.factor
    n2 = 24 * op2.factor # 6*4
    s1 = Struct(op2._endian + op2._analysis_code_fmt + b'i5f')
    s2 = Struct(op2._endian + b'i5f')
    add_sort_x = getattr(obj, 'add_sort' + str(op2.sort_method))
    add_eid_sort_x = getattr(obj, 'add_eid_sort' + str(op2.sort_method))
    for unused_i in range(nelements):
        edata = data[n:n+n1]
        n += n1

        out = s1.unpack(edata)
        eid_device = out[0]

        eid, dt = get_eid_dt_from_eid_device(
            eid_device, op2.nonlinear_factor, op2.sort_method)
        if op2.is_debug_file:
            op2.binary_debug.write('CBEAM-2 - eid=%i out=%s\n' % (eid, str(out)))

        #(grid, sd, sxc, sxd, sxe, sxf, smax, smin, mst, msc) = out
        #(grid, sd, sxc, sxd, sxe, sxf) = out

        # grid, sd, sxc, sxd, sxe, sxf = out[1:]
        add_eid_sort_x(dt, eid, *out[1:])

        for unused_inode in range(nnodes):
            edata = data[n:n+n2]
            n += n2
            out = s2.unpack(edata)
            # (grid, sd, sxc, sxd, sxe, sxf, smax, smin, mst, msc) = out
            # grid, sd, sxc, sxd, sxe, sxf = out
            add_sort_x(dt, eid, *out)
    return n
