from __future__ import annotations
from struct import Struct
from typing import TYPE_CHECKING
import numpy as np

from pyNastran.op2.op2_interface.op2_reader import mapfmt
from pyNastran.op2.op2_interface.utils import apply_mag_phase
from pyNastran.op2.tables.utils import get_is_slot_saved, get_eid_dt_from_eid_device
from pyNastran.op2.op2_helper import polar_to_real_imag

from pyNastran.op2.tables.oef_forces.oef_force_objects import (
    RealBendForceArray,
)
from pyNastran.op2.tables.oef_forces.oef_complex_force_objects import (
    ComplexCBendForceArray,
)
if TYPE_CHECKING:  # pragma: no cover
    from pyNastran.op2.op2 import OP2


def oef_cbend(op2: OP2, data, ndata, dt, is_magnitude_phase,
              result_type, prefix, postfix):
    """69-CBEND"""
    result_name = prefix + 'cbend_force' + postfix
    is_saved, slot = get_is_slot_saved(op2, result_name)
    if not is_saved:
        return ndata, None, None

    factor = op2.factor
    if op2.format_code == 1 and op2.num_wide == 15:  # real
        ntotal = 60 * factor  # 15*4
        nelements = ndata // ntotal
        auto_return, is_vectorized = op2._create_oes_object4(
            nelements, result_name, slot, RealBendForceArray)
        if auto_return:
            return nelements * op2.num_wide * 4, None, None

        obj = op2.obj
        if op2.use_vector and is_vectorized and op2.sort_method == 1:
            n = nelements * 4 * op2.num_wide
            itotal = obj.ielement
            ielement2 = obj.itotal + nelements
            itotal2 = ielement2

            floats = np.frombuffer(data, dtype=op2.fdtype8).reshape(nelements, 15).copy()
            obj._times[obj.itime] = dt
            if obj.itime == 0:
                ints = np.frombuffer(data, dtype=op2.idtype8).reshape(nelements, 15).copy()
                eids = ints[:, 0] // 10
                nids_a = ints[:, 1]
                nids_b = ints[:, 8]
                assert eids.min() > 0, eids.min()
                obj.element_node[itotal:itotal2, 0] = eids
                obj.element_node[itotal:itotal2, 1] = nids_a
                obj.element_node[itotal:itotal2, 2] = nids_b

            # [nid_a, bm1_a, bm2_a, ts1_a, ts2_a, af_a, trq_a,
            #  nid_b, bm1_b, bm2_b, ts1_b, ts2_b, af_b, trq_b]
            assert floats[:, 2:8].shape[1] == 6, floats[:, 2:8].shape
            assert floats[:, 9:].shape[1] == 6, floats[:, 9:].shape
            obj.data[obj.itime, itotal:itotal2, :6] = floats[:, 2:8]
            obj.data[obj.itime, itotal:itotal2, 6:] = floats[:, 9:]
            obj.itotal = itotal2
            obj.ielement = ielement2
        else:
            n = oef_cbend_real_15(op2, data, obj, nelements, ntotal)
    elif op2.format_code in [2, 3] and op2.num_wide == 27:  # imag
        # TODO: vectorize
        ntotal = 108 * factor  # 27*4
        nelements = ndata // ntotal

        auto_return, is_vectorized = op2._create_oes_object4(
            nelements, result_name, slot, ComplexCBendForceArray)
        if auto_return:
            return nelements * op2.num_wide * 4, None, None

        obj = op2.obj
        if op2.use_vector and is_vectorized and op2.sort_method == 1:
            n = nelements * 4 * op2.num_wide
            itotal = obj.ielement
            ielement2 = obj.itotal + nelements
            itotal2 = ielement2

            ireal = [2, 3, 4, 5, 6, 7, 15, 16, 17, 18, 19, 20]
            iimag = [8, 9, 10, 11, 12, 13, 21, 22, 23, 24, 25, 26]
            # 0    1
            # eid, nidA,
            # 2      3      4       5     6     7
            # 8      9      10      11    12    13
            # bm1Ar, bm2Ar, ts1Ar, ts2Ar, afAr, trqAr,
            # bm1Ai, bm2Ai, ts1Ai, ts2Ai, afAi, trqAi,
            # 14
            # nidB
            # 15     16     17     18     19    20
            # 21     22     23     24     25    26
            # bm1Br, bm2Br, ts1Br, ts2Br, afBr, trqBr,
            # bm1Bi, bm2Bi, ts1Bi, ts2Bi, afBi, trqBi
            floats = np.frombuffer(data, dtype=op2.fdtype8).reshape(nelements, 27).copy()
            obj._times[obj.itime] = dt
            if obj.itime == 0:
                ints = np.frombuffer(data, dtype=op2.idtype8).reshape(nelements, 27).copy()
                eids = ints[:, 0] // 10
                nids_a = ints[:, 1]
                nids_b = ints[:, 14]
                assert nids_a.min() > 0, nids_a
                assert nids_b.min() > 0, nids_b
                assert eids.min() > 0, eids.min()
                #print(nids_b)
                obj.element_node[itotal:itotal2, 0] = eids
                obj.element_node[itotal:itotal2, 1] = nids_a
                obj.element_node[itotal:itotal2, 2] = nids_b

            real_imag = apply_mag_phase(floats, is_magnitude_phase, ireal, iimag)
            obj.data[obj.itime, itotal:itotal2, :] = real_imag
            #assert floats[:, 1:6].shape[1] == 12, floats[:, 1:6].shape
            #assert floats[:, 7:].shape[1] == 6, floats[:, 7:].shape
            #obj.data[obj.itime, itotal:itotal2, :6] = floats[:, 1:6]
            #obj.data[obj.itime, itotal:itotal2, 6:] = floats[:, 7:]
            obj.itotal = itotal2
            obj.ielement = ielement2
        else:
            n = oef_cbend_imag_27(op2, data, obj,
                                  nelements, ntotal, is_magnitude_phase)
    else:  # pragma: no cover
        raise RuntimeError(op2.code_information())
        #msg = op2.code_information()
        #return op2._not_implemented_or_skip(data, ndata, msg), None, None
    return n, nelements, ntotal


def oef_cbend_real_15(op2: OP2, data: bytes,
                      obj: RealBendForceArray,
                      nelements: int, ntotal: int) -> int:
    n = 0
    fmt = mapfmt(op2._analysis_code_fmt + b' i6fi6f', op2.size)
    structi = Struct(op2._endian + fmt)
    add_sort_x = getattr(obj, 'add_sort' + str(op2.sort_method))
    for unused_i in range(nelements):
        edata = data[n:n + ntotal]

        out = structi.unpack(edata)
        if op2.is_debug_file:
            op2.binary_debug.write('OEF_BEND-69 - %s\n' % (str(out)))
        (eid_device,
         nid_a, bm1_a, bm2_a, ts1_a, ts2_a, af_a, trq_a,
         nid_b, bm1_b, bm2_b, ts1_b, ts2_b, af_b, trq_b) = out
        eid, dt = get_eid_dt_from_eid_device(
            eid_device, op2.nonlinear_factor, op2.sort_method)

        add_sort_x(
            dt, eid,
            nid_a, bm1_a, bm2_a, ts1_a, ts2_a, af_a, trq_a,
            nid_b, bm1_b, bm2_b, ts1_b, ts2_b, af_b, trq_b)
        n += ntotal
    return n

def oef_cbend_imag_27(op2: OP2, data: bytes,
                      obj: ComplexCBendForceArray,
                      nelements: int, ntotal: int,
                      is_magnitude_phase: bool) -> int:
    n = 0
    fmt = mapfmt(op2._analysis_code_fmt + b' i12f i12f', op2.size)
    structi = Struct(op2._endian + fmt)
    add_sort_x = getattr(obj, 'add_sort' + str(op2.sort_method))
    for unused_i in range(nelements):
        edata = data[n:n + 108]
        n += ntotal
        out = structi.unpack(edata)
        if op2.is_debug_file:
            op2.binary_debug.write('OEF_BEND-69 - %s\n' % (str(out)))
        (eid_device, nid_a,
         bm1_ar, bm2_ar, ts1_ar, ts2_ar, af_ar, trq_ar,
         bm1_ai, bm2_ai, ts1_ai, ts2_ai, af_ai, trq_ai,
         nid_b,
         bm1_br, bm2_br, ts1_br, ts2_br, af_br, trq_br,
         bm1_bi, bm2_bi, ts1_bi, ts2_bi, af_bi, trq_bi) = out
        eid, dt = get_eid_dt_from_eid_device(
            eid_device, op2.nonlinear_factor, op2.sort_method)

        if is_magnitude_phase:
            bm1_a = polar_to_real_imag(bm1_ar, bm1_ai)
            bm1_b = polar_to_real_imag(bm1_br, bm1_bi)
            bm2_a = polar_to_real_imag(bm2_ar, bm2_ai)
            bm2_b = polar_to_real_imag(bm2_br, bm2_bi)
            ts1_a = polar_to_real_imag(ts1_ar, ts1_ai)
            ts1_b = polar_to_real_imag(ts1_br, ts1_bi)
            ts2_a = polar_to_real_imag(ts2_ar, ts2_ai)
            ts2_b = polar_to_real_imag(ts2_br, ts2_bi)
            af_a = polar_to_real_imag(af_ar, af_ai)
            af_b = polar_to_real_imag(af_br, af_bi)
            trq_a = polar_to_real_imag(trq_ar, trq_ai)
            trq_b = polar_to_real_imag(trq_br, trq_bi)
        else:
            bm1_a = complex(bm1_ar, bm1_ai)
            bm1_b = complex(bm1_br, bm1_bi)
            bm2_a = complex(bm2_ar, bm2_ai)
            bm2_b = complex(bm2_br, bm2_bi)
            ts1_a = complex(ts1_ar, ts1_ai)
            ts1_b = complex(ts1_br, ts1_bi)
            ts2_a = complex(ts2_ar, ts2_ai)
            ts2_b = complex(ts2_br, ts2_bi)
            af_a = complex(af_ar, af_ai)
            af_b = complex(af_br, af_bi)
            trq_a = complex(trq_ar, trq_ai)
            trq_b = complex(trq_br, trq_bi)

        add_sort_x(dt, eid,
                   nid_a, bm1_a, bm2_a, ts1_a, ts2_a, af_a, trq_a,
                   nid_b, bm1_b, bm2_b, ts1_b, ts2_b, af_b, trq_b)
    return n
