from __future__ import annotations
from struct import Struct
from typing import Any, TYPE_CHECKING
import numpy as np

from pyNastran.op2.op2_interface.op2_reader import mapfmt
from pyNastran.op2.op2_helper import polar_to_real_imag

from pyNastran.op2.tables.utils import get_is_slot_saved, get_eid_dt_from_eid_device
from pyNastran.op2.op2_interface.utils import (
    apply_mag_phase,
)
from pyNastran.op2.tables.oes_stressStrain.utils import obj_set_element
from pyNastran.op2.tables.oes_stressStrain.real.oes_shear import RealShearStrainArray, RealShearStressArray

from pyNastran.op2.tables.oes_stressStrain.complex.oes_shear import ComplexShearStressArray, ComplexShearStrainArray
from pyNastran.op2.tables.oes_stressStrain.random.oes_shear import RandomShearStressArray, RandomShearStrainArray

if TYPE_CHECKING:  # pragma: no cover
    from pyNastran.op2.op2 import OP2



def oes_cshear_4(op2: OP2, data, ndata, dt, is_magnitude_phase,
                 result_type, prefix, postfix):
    """
    reads stress/strain for element type:
     - 4 : CSHEAR
    """
    # n = 0
    factor = op2.factor
    # 4-CSHEAR
    obj_vector_real = RealShearStressArray if op2.is_stress else RealShearStrainArray
    obj_vector_complex = ComplexShearStressArray if op2.is_stress else ComplexShearStressArray
    obj_vector_random = RandomShearStressArray if op2.is_stress else RandomShearStrainArray
    stress_strain = 'stress' if op2.is_stress else 'strain'
    result_name = f'{prefix}cshear_{stress_strain}{postfix}'

    is_saved, slot = get_is_slot_saved(op2, result_name)
    if not is_saved:
        return ndata, None, None

    if result_type == 0 and op2.num_wide == 4:  # real
        ntotal = 16 * factor  # 4*4
        nelements = ndata // ntotal
        auto_return, is_vectorized = op2._create_oes_object4(
            nelements, result_name, slot, obj_vector_real)
        if auto_return:
            return nelements * ntotal, None, None

        obj: RealShearStressArray | RealShearStrainArray = op2.obj
        assert obj is not None
        if op2.use_vector and is_vectorized and op2.sort_method == 1:
            n = nelements * ntotal
            itotal = obj.ielement
            ielement2 = obj.itotal + nelements
            itotal2 = ielement2

            floats = np.frombuffer(data, dtype=op2.fdtype8).reshape(nelements, 4)
            itime = obj.itime
            obj._times[itime] = dt
            obj_set_element(op2, obj, itotal, itotal2, data, nelements)

            # [max_strain, avg_strain, margin]
            obj.data[itime, itotal:itotal2, :] = floats[:, 1:].copy()
            obj.itotal = itotal2
            obj.ielement = ielement2
        elif op2.use_vector and is_vectorized and op2.sort_method == 2:
            op2.log.debug('vectorize CSHEAR real SORT%s' % op2.sort_method)
            n = oes_cshear_real_4(op2, data, obj,
                                  ntotal, nelements, dt)
        else:
            if is_vectorized and op2.use_vector:  # pragma: no cover
                op2.log.debug('vectorize CSHEAR real SORT%s' % op2.sort_method)
            n = oes_cshear_real_4(op2, data, obj,
                                  ntotal, nelements, dt)

    elif result_type == 1 and op2.num_wide == 5:  # imag
        ntotal = 20 * factor  # 4*5
        nelements = ndata // ntotal
        auto_return, is_vectorized = op2._create_oes_object4(
            nelements, result_name, slot, obj_vector_complex)
        if auto_return:
            return nelements * ntotal, None, None

        obj = op2.obj
        if op2.use_vector and is_vectorized and op2.sort_method == 1:
            n = nelements * ntotal
            itotal = obj.ielement
            ielement2 = obj.itotal + nelements
            itotal2 = ielement2

            floats = np.frombuffer(data, dtype=op2.fdtype8).reshape(nelements, 5).copy()
            obj._times[obj.itime] = dt
            obj_set_element(op2, obj, itotal, itotal2, data, nelements)

            # (eid_device, etmaxr, etmaxi, etavgr, etavgi)
            real_imag = apply_mag_phase(floats, is_magnitude_phase, [1, 3], [2, 4])
            obj.data[obj.itime, itotal:itotal2, :] = real_imag
            obj.itotal = itotal2
            obj.ielement = ielement2
        else:
            if is_vectorized and op2.use_vector:  # pragma: no cover
                op2.log.debug('vectorize CSHEAR imag SORT%s' % op2.sort_method)
            n = oes_cshear_complex_5(op2, data, obj,
                                     nelements, ntotal, is_magnitude_phase)
    elif result_type == 2 and op2.num_wide == 3:  # random
        ntotal = 12 * factor # 3*4
        nelements = ndata // ntotal
        auto_return, is_vectorized = op2._create_oes_object4(
            nelements, result_name, slot, obj_vector_random)
        if auto_return:
            assert ntotal == op2.num_wide * 4 * factor
            return nelements * ntotal, None, None

        obj = op2.obj
        assert obj is not None
        if op2.use_vector and is_vectorized and op2.sort_method == 1:
            n = nelements * 4 * op2.num_wide * factor
            itotal = obj.ielement
            ielement2 = obj.itotal + nelements
            itotal2 = ielement2

            floats = np.frombuffer(data, dtype=op2.fdtype).reshape(nelements, 3)
            itime = obj.itime
            obj._times[itime] = dt
            obj_set_element(op2, obj, itotal, itotal2, data, nelements)

            # [max_strain, avg_strain, margin]
            obj.data[itime, itotal:itotal2, :] = floats[:, 1:].copy()
            obj.itotal = itotal2
            obj.ielement = ielement2
        else:
            if is_vectorized and op2.use_vector:  # pragma: no cover
                op2.log.debug('vectorize CSHEAR random SORT%s' % op2.sort_method)
            n = oes_cshear_random_3(op2, data, obj,
                                    nelements, ntotal)
    else:  # pragma: no cover
        raise RuntimeError(op2.code_information())
    return n, nelements, ntotal


def oes_cshear_real_4(op2: OP2, data: bytes,
                      obj: RealShearStressArray | RealShearStrainArray,
                      ntotal: int, nelements: int, dt: Any) -> int:
    n = 0
    fmt = mapfmt(op2._analysis_code_fmt + b'3f', op2.size)
    struct1 = Struct(op2._endian + fmt)
    add_sort_x = getattr(obj, 'add_sort' + str(op2.sort_method))
    for unused_i in range(nelements):
        edata = data[n:n + ntotal]
        out = struct1.unpack(edata)  # num_wide=5
        if op2.is_debug_file:
            op2.binary_debug.write('CSHEAR-4 - %s\n' % str(out))

        (eid_device, max_strain, avg_strain, margin) = out
        #print(eid_device, max_strain, avg_strain, margin)
        eid, dt = get_eid_dt_from_eid_device(
            eid_device, op2.nonlinear_factor, op2.sort_method)
        add_sort_x(dt, eid, max_strain, avg_strain, margin)
        n += ntotal
    return n

def oes_cshear_complex_5(op2: OP2,
                         data: bytes,
                         obj: ComplexShearStressArray | ComplexShearStrainArray,
                         nelements: int, ntotal: int,
                         is_magnitude_phase: bool) -> int:
    n = 0
    fmt = mapfmt(op2._analysis_code_fmt + b'4f', op2.size)
    struct1 = Struct(op2._endian + fmt)
    add_sort_x = getattr(obj, 'add_sort' + str(op2.sort_method))
    for unused_i in range(nelements):
        edata = data[n:n + ntotal]
        out = struct1.unpack(edata)  # num_wide=5
        if op2.is_debug_file:
            op2.binary_debug.write('CSHEAR-4 - %s\n' % str(out))
        (eid_device, etmaxr, etmaxi, etavgr, etavgi) = out
        eid, dt = get_eid_dt_from_eid_device(
            eid_device, op2.nonlinear_factor, op2.sort_method)

        if is_magnitude_phase:
            etmax = polar_to_real_imag(etmaxr, etmaxi)
            etavg = polar_to_real_imag(etavgr, etavgi)
        else:
            etmax = complex(etmaxr, etmaxi)
            etavg = complex(etavgr, etavgi)
        add_sort_x(dt, eid, etmax, etavg)
        n += ntotal
    return n


def oes_cshear_random_3(op2: OP2, data: bytes,
                        obj: RandomShearStressArray | RandomShearStrainArray,
                        nelements: int, ntotal: int) -> int:
    n = 0
    add_sort_x = getattr(obj, 'add_sort' + str(op2.sort_method))
    fmt = mapfmt(op2._analysis_code_fmt + b'2f', op2.size)
    struct1 = Struct(op2._endian + fmt)
    for unused_i in range(nelements):
        edata = data[n:n + ntotal]
        out = struct1.unpack(edata)  # num_wide=5
        if op2.is_debug_file:
            op2.binary_debug.write('CSHEAR-4 - %s\n' % str(out))

        (eid_device, max_strain, avg_strain) = out
        eid, dt = get_eid_dt_from_eid_device(
            eid_device, op2.nonlinear_factor, op2.sort_method)
        #print(eid, dt)
        add_sort_x(dt, eid, max_strain, avg_strain)
        n += ntotal
    return n
