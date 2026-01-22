from __future__ import annotations
from struct import Struct
from typing import TYPE_CHECKING
import numpy as np

from pyNastran.op2.op2_interface.utils import mapfmt, real_imag_from_list
from pyNastran.op2.op2_helper import polar_to_real_imag

from pyNastran.op2.tables.utils import get_is_slot_saved, get_eid_dt_from_eid_device

from pyNastran.op2.tables.oes_stressStrain.real.oes_bend import RealBendStressArray, RealBendStrainArray
from pyNastran.op2.tables.oes_stressStrain.complex.oes_bend import ComplexBendStressArray, ComplexBendStrainArray
from pyNastran.op2.tables.oes_stressStrain.random.oes_bend import RandomBendStressArray, RandomBendStrainArray

if TYPE_CHECKING:  # pragma: no cover
    from pyNastran.op2.op2 import OP2


def oes_cbend_69(op2: OP2, data, ndata, dt,
                 is_magnitude_phase,
                 result_type, prefix, postfix):
    """
    reads stress/strain for element type:
     - 69 : CBEND

    """
    stress_strain = 'stress' if op2.is_stress else 'strain'
    result_name = f'{prefix}cbend_{stress_strain}{postfix}'
    is_saved, slot = get_is_slot_saved(op2, result_name)
    if not is_saved:
        return ndata, None, None

    if op2.is_stress:
        obj_vector_real = RealBendStressArray
        obj_vector_complex = ComplexBendStressArray
        obj_vector_random = RandomBendStressArray
    else:
        obj_vector_real = RealBendStrainArray
        obj_vector_complex = ComplexBendStrainArray
        obj_vector_random = RandomBendStrainArray

    # print(op2.code_information())
    factor = op2.factor
    if result_type == 0 and op2.num_wide == 21:  # real
        # TCODE,7 =0 Real
        # 2 GRID I External Grid Point identification number
        # 3 CA RS Circumferential Angle
        # 4 EC RS Long. strain at Point C
        # 5 ED RS Long. strain at Point D
        # 6 EE RS Long. strain at Point E
        # 7 EF RS Long. strain at Point F
        # 8 EMAX RS Maximum strain
        # 9 EMIN RS Minimum strain
        # 10 MST RS Margin of Safety in Tension
        # 11 MSC RS Margin of Safety in Compression
        # Words 2 through 11 repeat 002 times
        # n = 0
        ntotal = 84 * factor  # 4*21
        nelements = ndata // ntotal
        assert ndata % ntotal == 0, 'ndata=%s ntotal=%s nelements=%s error=%s' % (
        ndata, ntotal, nelements, ndata % ntotal)

        # nlayers = nelements * 2
        if op2.sort_method == 2:
            op2.log.warning('real cbend stress/strain for SORT2 is not supported')
            # print(op2.code_information())
            return nelements * ntotal, None, None

        auto_return, is_vectorized = op2._create_oes_object4(
            nelements, result_name, slot, obj_vector_real)
        if auto_return:
            return nelements * ntotal, None, None

        obj = op2.obj
        assert obj is not None
        if op2.use_vector and is_vectorized and op2.sort_method == 1 and 0:
            n = nelements * ntotal
            itotal = obj.ielement
            ielement2 = obj.itotal + nelements
            itotal2 = ielement2

            floats = np.frombuffer(data, dtype=op2.fdtype).reshape(nelements, 4)
            itime = obj.itime
            obj._times[itime] = dt
            if itime == 0:
                ints = np.frombuffer(data, dtype=op2.idtype).reshape(nelements, 4)
                eids = ints[:, 0] // 10
                assert eids.min() > 0, eids.min()
                obj.element[itotal:itotal2] = eids

            # [max_strain, avg_strain, margin]
            obj.data[itime, itotal:itotal2, :] = floats[:, 1:].copy()
            obj.itotal = itotal2
            obj.ielement = ielement2
        else:
            n = oes_cbend_real_21(op2, data, obj,
                                  nelements, dt)

        # msg = ''
        # if op2.read_mode == 2:
        # msg = op2.code_information()
        # n = op2._not_implemented_or_skip(data, ndata, msg)
        # return n, None, None
    elif result_type == 1 and op2.num_wide == 21:  # complex
        n = 0
        ntotal = 84 * factor  # 4*21
        nelements = ndata // ntotal
        assert ndata % ntotal == 0, 'ndata=%s ntotal=%s nelements=%s error=%s' % (
        ndata, ntotal, nelements, ndata % ntotal)
        # TCODE,7 =1 Real / Imaginary
        # 2 GRID  I External Grid Point identification number
        # 3   CA RS Circumferential Angle
        # 4  SCR RS Long. Stress at Point C
        # 5  SDR RS Long. Stress at Point D
        # 6  SER RS Long. Stress at Point E
        # 7  SFR RS Long. Stress at Point F
        # 8  SCI RS Long. Stress at Point C
        # 9  SDI RS Long. Stress at Point D
        # 10 SEI RS Long. Stress at Point E
        # 11 SFI RS Long. Stress at Point F
        # Words 2 through 11 repeat 002 times

        auto_return, is_vectorized = op2._create_oes_object4(
            nelements, result_name, slot, obj_vector_complex)
        if auto_return:
            assert ntotal == op2.num_wide * 4
            return nelements * ntotal, None, None

        obj = op2.obj
        assert obj is not None
        if op2.use_vector and is_vectorized and op2.sort_method == 1 and 0:
            n = nelements * 4 * op2.num_wide
            itotal = obj.ielement
            ielement2 = obj.itotal + nelements
            itotal2 = ielement2

            floats = np.frombuffer(data, dtype=op2.fdtype).reshape(nelements, 4)
            itime = obj.itime
            obj._times[itime] = dt

            if itime == 0:
                ints = np.frombuffer(data, dtype=op2.idtype).reshape(nelements, 4)
                eids = ints[:, 0] // 10
                assert eids.min() > 0, eids.min()
                obj.element[itotal:itotal2] = eids

            # [max_strain, avg_strain, margin]
            obj.data[itime, itotal:itotal2, :] = floats[:, 1:].copy()
            obj.itotal = itotal2
            obj.ielement = ielement2
        else:
            n = oes_cbend_complex_21(op2, data, obj, nelements, ntotal,
                                     is_magnitude_phase)

    elif result_type == 2 and op2.num_wide == 13:
        ntotal = 52 * factor  # 4*13
        nelements = ndata // ntotal
        # TCODE,7 =2 Real
        # 2 GRID I External Grid Point identification number
        # 3 CA RS Circumferential Angle
        # 4 SC RS Long. Stress at Point C
        # 5 SD RS Long. Stress at Point D
        # 6 SE RS Long. Stress at Point E
        # 7 SF RS Long. Stress at Point F
        # Words 2 through 7 repeat 002 times
        # if op2.table_name != "OESPSD2":
        # msg = ''
        # if op2.read_mode == 2:
        # msg = op2.code_information()
        # n = op2._not_implemented_or_skip(data, ndata, msg)
        # return n, None, None

        auto_return, is_vectorized = op2._create_oes_object4(
            nelements, result_name, slot, obj_vector_random)
        if auto_return:
            assert ntotal == op2.num_wide * 4
            return nelements * ntotal, None, None

        obj = op2.obj
        if op2.use_vector and is_vectorized and op2.sort_method == 1 and 0:
            n = nelements * 4 * op2.num_wide
            itotal = obj.ielement
            ielement2 = obj.itotal + nelements
            itotal2 = ielement2

            floats = np.frombuffer(data, dtype=op2.fdtype).reshape(nelements, 4)
            itime = obj.itime
            obj._times[itime] = dt
            if itime == 0:
                ints = np.frombuffer(data, dtype=op2.idtype).reshape(nelements, 4)
                eids = ints[:, 0] // 10
                assert eids.min() > 0, eids.min()
                obj.element[itotal:itotal2] = eids

            # [max_strain, avg_strain, margin]
            obj.data[itime, itotal:itotal2, :] = floats[:, 1:].copy()
            obj.itotal = itotal2
            obj.ielement = ielement2
        else:
            n = oes_cbend_random_13(
                op2, data, obj,
                nelements, dt)
    else:  # pragma: no cover
        raise RuntimeError(op2.code_information())
    return n, nelements, ntotal

def oes_cbend_real_21(op2: OP2, data: bytes,
                      obj: RealBendStressArray | RealBendStrainArray,
                      nelements: int, dt) -> int:
    n = 0
    factor = op2.factor
    size = op2.size
    ntotal1 = 4 * factor
    ntotal2 = 40 * factor
    struct1 = Struct(op2._endian + mapfmt(op2._analysis_code_fmt, size))
    struct2 = Struct(op2._endian + mapfmt(b'i9f', size))
    add_sort_x = getattr(obj, 'add_sort' + str(op2.sort_method))

    #print('ntimes =', nelements)
    for unused_i in range(nelements):
        edata = data[n:n + ntotal1]
        eid_device, = struct1.unpack(edata)
        eid, dt = get_eid_dt_from_eid_device(
            eid_device, op2.nonlinear_factor, op2.sort_method)
        n += ntotal1

        for unused_j in range(2):
            edata = data[n:n + ntotal2]
            out = struct2.unpack(edata)
            if op2.is_debug_file:
                op2.binary_debug.write('BEND-69 - eid=%s %s\n' % (eid, str(out)))
            #print('BEND-69 - eid=%s %s\n' % (eid, str(out)))

            (grid, angle, sc, sd, se, sf, omax, omin, mst, msc) = out
            add_sort_x(dt, eid, grid, angle, sc, sd, se, sf, omax, omin, mst, msc)
            n += ntotal2
    return n

def oes_cbend_complex_21(op2: OP2, data: bytes,
                         obj: ComplexBendStressArray | ComplexBendStrainArray,
                         nelements: int, ntotal: int, is_magnitude_phase: bool) -> int:
    size = op2.size
    factor = op2.factor
    ntotal1 = 4 * factor
    ntotal2 = 40 * factor
    struct1 = Struct(op2._endian + mapfmt(op2._analysis_code_fmt, size))
    struct2 = Struct(op2._endian + mapfmt(b'i9f', size))
    add_sort_x = getattr(obj, 'add_sort' + str(op2.sort_method))

    n = 0
    for unused_i in range(nelements):
        edata = data[n:n + ntotal1]
        eid_device, = struct1.unpack(edata)
        eid, dt = get_eid_dt_from_eid_device(
            eid_device, op2.nonlinear_factor, op2.sort_method)

        n += ntotal1
        for unused_j in range(2):
            edata = data[n:n + ntotal2]
            out = struct2.unpack(edata)
            if op2.is_debug_file:
                op2.binary_debug.write('BEND-69 - eid=%s %s\n' % (eid, str(out)))
            #print('BEND-69 - eid=%s %s\n' % (eid, str(out)))

            (grid, angle,
             scr, sdr, ser, sfr,
             sci, sdi, sei, sfi) = out

            sc, sd, se, sf = real_imag_from_list([
                scr, sdr, ser, sfr,
                sci, sdi, sei, sfi], is_magnitude_phase)
            add_sort_x(dt, eid, grid, angle, sc, sd, se, sf)
            n += ntotal2
    return n


def oes_cbend_random_13(op2: OP2, data: bytes,
                        obj: RandomBendStressArray | RandomBendStrainArray,
                        nelements: int, dt) -> int:
    n = 0
    factor = op2.factor
    size = op2.size
    ntotal1 = 4 * factor
    ntotal2 = 24 * factor
    struct1 = Struct(op2._endian + mapfmt(op2._analysis_code_fmt, size))
    struct2 = Struct(op2._endian + mapfmt(b'i5f', size))

    for unused_ielem in range(nelements):
        edata = data[n:n + ntotal1]
        # self.show_data(edata)
        eid_device, = struct1.unpack(edata)
        eid, dt = get_eid_dt_from_eid_device(
            eid_device, op2.nonlinear_factor, op2.sort_method)

        n += ntotal1
        for unused_inode in range(2):
            edata = data[n:n + ntotal2]
            out = struct2.unpack(edata)
            if op2.is_debug_file:
                op2.binary_debug.write('BEND-69 - eid=%s dt=%s %s\n' % (eid, dt, str(out)))
            # print('BEND-69 - eid=%s dt=%s %s\n' % (eid, dt, str(out)))

            (grid, angle, sc, sd, se, sf) = out
            obj.add_sort1(dt, eid, grid, angle, sc, sd, se, sf)
            n += ntotal2
    return n
