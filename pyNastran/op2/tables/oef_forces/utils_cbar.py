from __future__ import annotations
from struct import Struct
from typing import Any, TYPE_CHECKING
import numpy as np

from pyNastran.op2.op2_interface.utils import mapfmt, real_imag_from_list
from pyNastran.op2.tables.utils import get_is_slot_saved, get_eid_dt_from_eid_device
from pyNastran.op2.op2_helper import polar_to_real_imag

from pyNastran.op2.tables.oef_forces.oef_force_objects import (
    RealCBarForceArray, RealCBar100ForceArray,
    RealCFastForceArrayNX, RealCWeldForceArray,
    RealCWeldForceArrayMSC,
)
from pyNastran.op2.tables.oef_forces.oef_complex_force_objects import (
    ComplexCBarForceArray,
    ComplexCWeldForceArray,
    ComplexCWeldForceArrayMSC,
)

if TYPE_CHECKING:  # pragma: no cover
    from pyNastran.op2.op2 import OP2


def oef_cbar_34(op2: OP2, data: bytes, ndata: int, dt: Any,
                is_magnitude_phase: bool,
                result_type: str, prefix: str, postfix: str) -> tuple[int, int, int]:
    """
    34-CBAR
    117-CWELDC
    119-CFAST-NX  (126-CFAST-MSC is like the CBUSH)
    118-WELDP-MSC

    """
    factor = op2.factor
    n = 0
    if op2.element_type == 34:
        result_name = prefix + 'cbar_force' + postfix
        obj_real = RealCBarForceArray
        obj_complex = ComplexCBarForceArray
    elif op2.element_type in [117, 200]:
        result_name = prefix + 'cweld_force' + postfix
        obj_real = RealCWeldForceArray
        obj_complex = ComplexCWeldForceArray
        assert op2.num_wide in [9, 17], op2.code_information()
    elif op2.element_type == 118:  # WELDP
        result_name = prefix + 'cweld_force' + postfix
        obj_real = RealCWeldForceArrayMSC
        obj_complex = ComplexCWeldForceArrayMSC
    elif op2.element_type == 119:
        result_name = prefix + 'cfast_force' + postfix
        obj_real = RealCFastForceArrayNX
        assert op2.num_wide == 9, op2.code_information()
    else:
        raise NotImplementedError(op2.element_type)

    is_saved, slot = get_is_slot_saved(op2, result_name)
    if not is_saved:
        return ndata, None, None
    # print(result_type in [0, 2], op2.num_wide == 9)
    # print(result_type == 1, op2.num_wide == 17)
    if result_type in [0, 2] and op2.num_wide == 9:  # real/random
        # real - format_code == 1
        # random - format_code == 3
        # result_name, is_random = self._apply_oef_ato_crm_psd_rms_no(result_name)
        slot = op2.get_result(result_name)

        ntotal = 36 * factor  # 9*4
        nelements = ndata // ntotal
        auto_return, is_vectorized = op2._create_oes_object4(
            nelements, result_name, slot, obj_real)
        if auto_return:
            return nelements * ntotal, None, None

        obj = op2.obj
        # return nelements * op2.num_wide * 4
        if op2.use_vector and is_vectorized and op2.sort_method == 1:
            n = nelements * ntotal
            itotal = obj.ielement
            ielement2 = obj.itotal + nelements
            itotal2 = ielement2

            floats = np.frombuffer(data, dtype=op2.fdtype8).reshape(nelements, 9)
            obj._times[obj.itime] = dt
            if obj.itime == 0:
                ints = np.frombuffer(data, dtype=op2.idtype8).reshape(nelements, 9)
                eids = ints[:, 0] // 10
                assert eids.min() > 0, eids.min()
                obj.element[itotal:itotal2] = eids

            # [bm1a, bm2a, bm1b, bm2b, ts1, ts2, af, trq]
            obj.data[obj.itime, itotal:itotal2, :] = floats[:, 1:].copy()
            obj.itotal = itotal2
            obj.ielement = ielement2
        # elif op2.use_vector and is_vectorized and op2.sort_method == 1:
        else:
            n = oef_cbar_real_9(op2, data, obj, nelements, ntotal)
    elif result_type == 1 and op2.num_wide == 17:  # imag
        # TODO: vectorize
        ntotal = 68 * factor  # 17*4
        nelements = ndata // ntotal
        assert ndata % ntotal == 0

        auto_return, is_vectorized = op2._create_oes_object4(
            nelements, result_name, slot, obj_complex)
        if auto_return:
            return nelements * ntotal, None, None

        obj = op2.obj
        n = oef_cbar_imag_17(op2, data, obj, nelements, ntotal, is_magnitude_phase)
    else:
        raise RuntimeError(op2.code_information())
        # print(op2.table_name)
        # msg = op2.code_information()
        # print(msg)
        # print(result_type)
        # raise NotImplementedError(op2.code_information())
        # return op2._not_implemented_or_skip(data, ndata, msg), None, None
    # print self.barForces
    return n, nelements, ntotal


def oef_cbar_100(op2: OP2, data, ndata, dt, unused_is_magnitude_phase,
                 result_type, prefix, postfix):
    """100-BARS"""
    result_name = prefix + 'cbar_force' + postfix  # _10nodes
    is_saved, slot = get_is_slot_saved(op2, result_name)
    if not is_saved:
        return ndata, None, None

    factor = op2.factor
    if op2.format_code == 1 and op2.num_wide == 8:  # real
        ntotal = 32 * factor  # 8*4
        nelements = ndata // ntotal
        auto_return, is_vectorized = op2._create_oes_object4(
            nelements, result_name, slot, RealCBar100ForceArray)
        if auto_return:
            return nelements * ntotal, None, None

        obj = op2.obj
        if op2.use_vector and is_vectorized and op2.sort_method == 1:
            n = nelements * ntotal
            itotal = obj.ielement
            ielement2 = obj.itotal + nelements
            itotal2 = ielement2

            floats = np.frombuffer(data, dtype=op2.fdtype8).reshape(nelements, 8)
            obj._times[obj.itime] = dt
            if obj.itime == 0:
                ints = np.frombuffer(data, dtype=op2.idtype8).reshape(nelements, 8)
                eids = ints[:, 0] // 10
                assert eids.min() > 0, eids.min()
                obj.element[itotal:itotal2] = eids

            # [axial, torsion, SMa, SMt]
            obj.data[obj.itime, itotal:itotal2, :] = floats[:, 1:].copy()
            obj.itotal = itotal2
            obj.ielement = ielement2
        else:
            n = oef_cbar_100_real_8(op2, data, obj,
                                    nelements, ntotal)

    # elif op2.format_code in [2, 3] and op2.num_wide == 14:  # imag
    else:  # pragma: no cover
        raise RuntimeError(op2.code_information())
        # msg = op2.code_information()
        # print(msg)
        # return op2._not_implemented_or_skip(data, ndata, msg)
    return n, nelements, ntotal


def oef_cbar_real_9(op2: OP2, data: bytes,
                    obj: RealCBarForceArray,
                    nelements: int, ntotal: int) -> int:
    n = 0
    fmt = mapfmt(op2._endian + op2._analysis_code_fmt + b'8f', op2.size)
    s = Struct(fmt)
    add_sort_x = getattr(obj, 'add_sort' + str(op2.sort_method))
    for unused_i in range(nelements):
        edata = data[n:n + ntotal]

        out = s.unpack(edata)
        if op2.is_debug_file:
            op2.binary_debug.write('OEF_CBar - %s\n' % (str(out)))
        (eid_device, bm1a, bm2a, bm1b, bm2b, ts1, ts2, af, trq) = out
        eid, dt = get_eid_dt_from_eid_device(
            eid_device, op2.nonlinear_factor, op2.sort_method)
        #data_in = [eid, bm1a, bm2a, bm1b, bm2b, ts1, ts2, af, trq]
        add_sort_x(dt, eid, bm1a, bm2a, bm1b, bm2b, ts1, ts2, af, trq)
        n += ntotal
    return n

def oef_cbar_imag_17(op2: OP2, data: bytes,
                     obj: ComplexCBarForceArray,
                     nelements: int, ntotal: int,
                     is_magnitude_phase: bool) -> int:
    n = 0
    fmt = mapfmt(op2._endian + op2._analysis_code_fmt + b'16f', op2.size)
    s = Struct(fmt)
    add_sort_x = getattr(obj, 'add_sort' + str(op2.sort_method))
    for unused_i in range(nelements):
        edata = data[n:n + ntotal]

        out = s.unpack(edata)
        (eid_device,
         bm1ar, bm2ar, bm1br, bm2br, ts1r, ts2r, afr, trqr,
         bm1ai, bm2ai, bm1bi, bm2bi, ts1i, ts2i, afi, trqi) = out
        if op2.is_debug_file:
            op2.binary_debug.write('OEF_CBar - %s\n' % (str(out)))
        eid, dt = get_eid_dt_from_eid_device(
            eid_device, op2.nonlinear_factor, op2.sort_method)

        bm1a, bm2a, bm1b, bm2b, ts1, ts2, af, trq = real_imag_from_list([
            bm1ar, bm2ar, bm1br, bm2br, ts1r, ts2r, afr, trqr,
            bm1ai, bm2ai, bm1bi, bm2bi, ts1i, ts2i, afi, trqi], is_magnitude_phase)

        #data_in = [bm1a, bm2a, bm1b, bm2b, ts1, ts2, af, trq]
        #print("eid_device=%s eid=%s dt=%s %s" % (
            #eid_device, eid, dt, self.get_element_type(op2.element_type)), data_in)
        add_sort_x(dt, eid, bm1a, bm2a, bm1b, bm2b, ts1, ts2, af, trq)
        n += ntotal
    return n

def oef_cbar_100_real_8(self, data: bytes,
                        obj: RealCBar100ForceArray,
                        nelements: int, ntotal: int) -> int:
    op2 = self
    n = 0
    s = Struct(op2._endian + op2._analysis_code_fmt + b'7f')
    add_sort_x = getattr(obj, 'add_sort' + str(op2.sort_method))
    for unused_i in range(nelements):
        edata = data[n:n+ntotal]

        out = s.unpack(edata)
        if op2.is_debug_file:
            op2.binary_debug.write('OEF_CBar100 - %s\n' % (str(out)))
        (eid_device, sd, bm1, bm2, ts1, ts2, af, trq) = out
        eid, dt = get_eid_dt_from_eid_device(
            eid_device, op2.nonlinear_factor, op2.sort_method)
        add_sort_x(dt, eid, sd, bm1, bm2, ts1, ts2, af, trq)
        n += ntotal
    return n
