from __future__ import annotations
from struct import Struct
from typing import TYPE_CHECKING

import numpy as np
from pyNastran.op2.op2_interface.utils import (
    apply_mag_phase,
)
from pyNastran.op2.op2_interface.op2_reader import mapfmt
from pyNastran.op2.op2_helper import polar_to_real_imag
from pyNastran.op2.tables.utils import get_eid_dt_from_eid_device
from pyNastran.op2.tables.oes_stressStrain.utils import obj_set_element

from pyNastran.op2.tables.oes_stressStrain.real.oes_bush import RealBushStressArray, RealBushStrainArray
from pyNastran.op2.tables.oes_stressStrain.complex.oes_bush import (ComplexCBushStressArray, ComplexCBushStrainArray)
if TYPE_CHECKING:  # pragma: no cover
    from pyNastran.op2.op2 import OP2


def oes_cbush(op2: OP2,
              data, ndata, dt, is_magnitude_phase,
              result_type, prefix, postfix):
    """
    reads stress/strain for element type:
     - 102 : CBUSH

    """
    # n = 0
    factor = op2.factor
    if op2.is_stress:
        result_name = prefix + 'cbush_stress' + postfix
    else:
        result_name = prefix + 'cbush_strain' + postfix

    if op2._results.is_not_saved(result_name):
        return ndata, None, None
    op2._results._found_result(result_name)

    slot = op2.get_result(result_name)
    table_name = op2.table_name_str
    if result_type in [0, 2] and op2.num_wide == 7:  # real, random
        obj_vector_real = RealBushStressArray if op2.is_stress else RealBushStrainArray

        assert op2.num_wide == 7, "num_wide=%s not 7" % op2.num_wide
        ntotal = 28 * factor  # 4*7

        nelements = ndata // ntotal
        auto_return, is_vectorized = op2._create_oes_object4(
            nelements, result_name, slot, obj_vector_real)
        if auto_return:
            return nelements * ntotal, None, None
        obj = op2.obj

        if op2.use_vector and is_vectorized and op2.sort_method == 1:
            n = nelements * ntotal

            istart = obj.ielement
            iend = istart + nelements
            obj._times[obj.itime] = dt

            obj_set_element(op2, obj, istart, iend, data, nelements)

            floats = np.frombuffer(data, dtype=op2.fdtype8).reshape(nelements, 7)
            # [tx, ty, tz, rx, ry, rz]
            obj.data[obj.itime, istart:iend, :] = floats[:, 1:].copy()
        else:
            n = oes_cbush_real_7(op2, data, obj,
                                 nelements, ntotal, dt)
    elif result_type == 1 and op2.num_wide == 13:  # imag
        obj_complex = ComplexCBushStressArray if op2.is_stress else ComplexCBushStrainArray

        ntotal = 52 * factor  # 4*13
        nelements = ndata // ntotal
        auto_return, is_vectorized = op2._create_oes_object4(
            nelements, result_name, slot, obj_complex)
        if auto_return:
            return nelements * ntotal, None, None

        obj = op2.obj
        if op2.use_vector and is_vectorized and op2.sort_method == 1:
            n = nelements * ntotal
            itotal = obj.ielement
            ielement2 = obj.itotal + nelements
            itotal2 = ielement2

            floats = np.frombuffer(data, dtype=op2.fdtype8).reshape(nelements, 13).copy()
            obj._times[obj.itime] = dt
            obj_set_element(op2, obj, itotal, itotal2, data, nelements)

            isave1 = [1, 2, 3, 4, 5, 6]
            isave2 = [7, 8, 9, 10, 11, 12]
            real_imag = apply_mag_phase(floats, is_magnitude_phase, isave1, isave2)
            obj.data[obj.itime, itotal:itotal2, :] = real_imag

            obj.itotal = itotal2
            obj.ielement = ielement2
        else:
            n = oes_cbush_complex_13(op2, data, obj,
                                     nelements, ntotal,
                                     is_magnitude_phase)
    elif result_type == 1 and op2.num_wide == 7 and table_name == 'OESPSD1':
        result_type = 2
        n, nelements, ntotal = oes_cbush(op2, data, ndata, dt, is_magnitude_phase,
                                         result_type, prefix, postfix)
    else:  # pragma: no cover
        raise RuntimeError(op2.code_information())
        # msg = op2.code_information()
        # raise NotImplementedError(msg)
        # return op2._not_implemented_or_skip(data, ndata, msg)
    return n, nelements, ntotal

def oes_cbush_real_7(op2: OP2, data: bytes,
                     obj: RealBushStressArray | RealBushStrainArray,
                     nelements: int, ntotal: int, dt, debug: bool=False) -> int:
    n = 0
    struct1 = Struct(op2._endian + mapfmt(op2._analysis_code_fmt + b'6f', op2.size))
    add_sort_x = getattr(obj, 'add_sort' + str(op2.sort_method))
    nonlinear_factor = op2.nonlinear_factor
    #print(add_sort_x)
    #print('obj.is_sort1 =', obj.is_sort1, obj.table_name)
    for unused_i in range(nelements):
        edata = data[n:n + ntotal]
        n += ntotal

        out = struct1.unpack(edata)  # num_wide=7
        if op2.is_debug_file:
            op2.binary_debug.write('CBUSH-102 - %s\n' % str(out))

        (eid_device, tx, ty, tz, rx, ry, rz) = out
        eid, dt = get_eid_dt_from_eid_device(
            eid_device, nonlinear_factor, op2.sort_method)
        if debug:  # pragma: no cover
            print(f'CBUSH: eid_device={eid_device} eid={eid} dt={nonlinear_factor} nf={nonlinear_factor} -> {obj.data.shape}')

        add_sort_x(dt, eid, tx, ty, tz, rx, ry, rz)
    return n


def oes_cbush_complex_13(op2: OP2,
                         data: bytes,
                         obj: ComplexCBushStressArray | ComplexCBushStrainArray,
                         nelements: int, ntotal: int,
                         is_magnitude_phase: bool) -> int:
    n = 0
    struct1 = Struct(op2._endian + mapfmt(op2._analysis_code_fmt + b'12f', op2.size))
    add_sort_x = getattr(obj, 'add_sort' + str(op2.sort_method))
    for unused_i in range(nelements):
        edata = data[n:n + ntotal]
        out = struct1.unpack(edata)  # num_wide=7
        if op2.is_debug_file:
            op2.binary_debug.write('CBUSH-102 - %s\n' % str(out))

        (eid_device,
         txr, tyr, tzr, rxr, ryr, rzr,
         txi, tyi, tzi, rxi, ryi, rzi) = out
        eid, dt = get_eid_dt_from_eid_device(
            eid_device, op2.nonlinear_factor, op2.sort_method)

        if is_magnitude_phase:
            tx = polar_to_real_imag(txr, txi)
            ty = polar_to_real_imag(tyr, tyi)
            tz = polar_to_real_imag(tzr, tzi)
            rx = polar_to_real_imag(rxr, rxi)
            ry = polar_to_real_imag(ryr, ryi)
            rz = polar_to_real_imag(rzr, rzi)
        else:
            tx = complex(txr, txi)
            ty = complex(tyr, tyi)
            tz = complex(tzr, tzi)
            rx = complex(rxr, rxi)
            ry = complex(ryr, ryi)
            rz = complex(rzr, rzi)
        add_sort_x(dt, eid, tx, ty, tz, rx, ry, rz)
        n += ntotal
    return ntotal
