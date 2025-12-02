from __future__ import annotations
from struct import Struct
from typing import TYPE_CHECKING

import numpy as np
from pyNastran.op2.op2_interface.utils import (
    apply_mag_phase,
)
from pyNastran.op2.op2_interface.op2_reader import mapfmt
# from pyNastran.op2.op2_interface.op2_reader import mapfmt
from pyNastran.op2.op2_helper import polar_to_real_imag
from pyNastran.op2.tables.utils import get_is_slot_saved, get_eid_dt_from_eid_device
from pyNastran.op2.tables.oes_stressStrain.utils import obj_set_element
from pyNastran.op2.tables.oes_stressStrain.real.oes_bush1d import RealBush1DStressArray
from pyNastran.op2.tables.oes_stressStrain.complex.oes_bush1d import ComplexCBush1DStressArray

if TYPE_CHECKING:  # pragma: no cover
    from pyNastran.op2.op2 import OP2




def oes_cbush1d(op2: OP2, data, ndata, dt, is_magnitude_phase,
                result_type, prefix, postfix):
    """
    reads stress/strain for element type:
     - 40 : CBUSH1D
    """
    factor = op2.factor
    if op2.is_stress:
        result_name = f'{prefix}cbush1d_stress_strain{postfix}'
    else:
        result_name = f'{prefix}cbush1d_stress_strain{postfix}'

    is_saved, slot = get_is_slot_saved(op2, result_name)
    if not is_saved:
        return ndata, None, None

    if result_type == 0 and op2.num_wide == 8:  # real
        if op2.is_stress:
            obj_vector_real = RealBush1DStressArray
        else:
            # op2.create_transient_object(self.cbush1d_stress_strain, Bush1DStrain)  # undefined
            raise NotImplementedError('cbush1d_stress_strain; numwide=8')

        ntotal = 32 * factor  # 4*8
        nelements = ndata // ntotal
        auto_return, is_vectorized = op2._create_oes_object4(
            nelements, result_name, slot, obj_vector_real)
        if auto_return:
            return nelements * ntotal, None, None

        obj = op2.obj
        if op2.use_vector and is_vectorized and op2.sort_method == 1:
            n = nelements * ntotal

            itotal = obj.itotal
            itotal2 = itotal + nelements
            itime = obj.itime
            obj._times[itime] = dt

            if 1:  # obj.itime == 0:
                ints = np.frombuffer(data, dtype=op2.idtype).reshape(nelements, 8).copy()
                eids = ints[:, 0] // 10
                fail = ints[:, 7]
                obj.element[itotal:itotal2] = eids
                obj.is_failed[itime, itotal:itotal2, 0] = fail

            floats = np.frombuffer(data, dtype=op2.fdtype).reshape(nelements, 8)
            # [xxx, fe, ue, ve, ao, ae, ep, xxx]
            obj.data[itime, itotal:itotal2, :] = floats[:, 1:7].copy()

            obj.ielement = itotal2
            obj.itotal = itotal2
        else:
            n = oes_cbush1d_real_8(op2, data, obj, nelements, ntotal)
    elif result_type == 1 and op2.num_wide == 9:  # imag
        # TODO: vectorize object
        ntotal = 36 * factor  # 4*9
        nelements = ndata // ntotal

        if op2.is_stress:
            auto_return, is_vectorized = op2._create_oes_object4(
                nelements, result_name, slot, ComplexCBush1DStressArray)
        else:
            raise NotImplementedError('self.cbush1d_stress_strain; complex strain')

        if auto_return:
            assert ntotal == op2.num_wide * 4
            return nelements * ntotal, None, None

        obj = op2.obj
        if op2.use_vector and is_vectorized and op2.sort_method == 1:
            n = nelements * op2.num_wide * 4

            itotal = obj.itotal
            itotal2 = itotal + nelements
            itime = obj.itime
            obj._times[itime] = dt
            obj_set_element(op2, obj, itotal, itotal2, data, nelements)

            floats = np.frombuffer(data, dtype=op2.fdtype).reshape(nelements, 9).copy()
            # [fer, uer, aor, aer,
            # fei, uei, aoi, aei]
            isave1 = [1, 3, 5, 7]
            isave2 = [2, 4, 6, 8]
            real_imag = apply_mag_phase(floats, is_magnitude_phase, isave1, isave2)
            obj.data[obj.itime, itotal:itotal2, :] = real_imag

            obj.ielement = itotal2
            obj.itotal = itotal2
        else:
            n = oes_cbush1d_complex_9(
                op2, data, obj, nelements, ntotal,
                is_magnitude_phase)
    else:  # pragma: no cover
        raise RuntimeError(op2.code_information())
        # msg = op2.code_information()
        # raise NotImplementedError(msg)
        # return op2._not_implemented_or_skip(data, ndata, msg)
    return n, nelements, ntotal


def oes_cbush1d_real_8(op2: OP2, data: bytes,
                       obj: RealBush1DStressArray,
                       nelements: int, ntotal: int) -> int:
    n = 0
    fmt = mapfmt(op2._analysis_code_fmt + b'6fi', op2.size)
    struct1 = Struct(op2._endian + fmt)
    for unused_i in range(nelements):
        edata = data[n:n + ntotal]
        out = struct1.unpack(edata)  # num_wide=25
        if op2.is_debug_file:
            op2.binary_debug.write('CBUSH1D-40 - %s\n' % (str(out)))
        (eid_device, fe, ue, ve, ao, ae, ep, fail) = out
        eid, dt = get_eid_dt_from_eid_device(
            eid_device, op2.nonlinear_factor, op2.sort_method)

        # axial_force, axial_displacement, axial_velocity, axial_stress,
        # axial_strain, plastic_strain, is_failed
        obj.add_sort1(dt, eid, fe, ue, ve, ao, ae, ep, fail)
        n += ntotal
    return n

def oes_cbush1d_complex_9(op2: OP2, data: bytes,
                          obj: RealBush1DStressArray,
                          nelements: int, ntotal: int,
                          is_magnitude_phase: bool) -> int:
    n = 0
    fmt = mapfmt(op2._analysis_code_fmt + b'8f', op2.size)
    struct1 = Struct(op2._endian + fmt)
    for unused_i in range(nelements):
        edata = data[n:n + ntotal]

        out = struct1.unpack(edata)  # num_wide=25
        (eid_device,
         fer, uer, aor, aer,
         fei, uei, aoi, aei) = out
        eid, dt = get_eid_dt_from_eid_device(
            eid_device, op2.nonlinear_factor, op2.sort_method)

        if is_magnitude_phase:
            fe = polar_to_real_imag(fer, fei)
            ue = polar_to_real_imag(uer, uei)
            ao = polar_to_real_imag(aor, aoi)
            ae = polar_to_real_imag(aer, aei)
        else:
            fe = complex(fer, fei)
            ue = complex(uer, uei)
            ao = complex(aor, aoi)
            ae = complex(aer, aei)
        obj.add_new_eid(op2.element_type, dt, eid, fe, ue, ao, ae)
    return n
