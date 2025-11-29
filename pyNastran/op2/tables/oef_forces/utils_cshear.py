from __future__ import annotations
from struct import Struct
from typing import Any, TYPE_CHECKING
import numpy as np

from pyNastran.op2.op2_interface.op2_reader import mapfmt
from pyNastran.op2.op2_interface.utils import apply_mag_phase
from pyNastran.op2.tables.utils import get_eid_dt_from_eid_device
from pyNastran.op2.op2_helper import polar_to_real_imag

from pyNastran.op2.tables.oef_forces.oef_force_objects import (
    RealCShearForceArray,
)
from pyNastran.op2.tables.oef_forces.oef_complex_force_objects import (
    ComplexCShearForceArray,
)

if TYPE_CHECKING:  # pragma: no cover
    from pyNastran.op2.op2 import OP2


def oef_cshear(op2: OP2, data, ndata, dt, is_magnitude_phase,
               result_type, prefix, postfix):
    """4-CSHEAR"""
    result_name = prefix + 'cshear_force' + postfix
    if op2._results.is_not_saved(result_name):
        return ndata, None, None
    op2._results._found_result(result_name)
    slot = op2.get_result(result_name)

    n = 0
    factor = op2.factor
    if op2.format_code == 1 and op2.num_wide == 17:  # real
        ntotal = 68  # 17*4
        nelements = ndata // ntotal

        obj_real = RealCShearForceArray
        auto_return, is_vectorized = op2._create_oes_object4(
            nelements, result_name, slot, obj_real)
        if auto_return:
            return nelements * op2.num_wide * 4, None, None

        obj = op2.obj
        if op2.use_vector and is_vectorized and op2.sort_method == 1:
            n = nelements * 4 * op2.num_wide
            itotal = obj.ielement
            ielement2 = obj.itotal + nelements
            itotal2 = ielement2

            floats = np.frombuffer(data, dtype=op2.fdtype8).reshape(nelements, 17)
            obj._times[obj.itime] = dt
            if obj.itime == 0:
                ints = np.frombuffer(data, dtype=op2.idtype8).reshape(nelements, 17)
                eids = ints[:, 0] // 10
                assert eids.min() > 0, eids.min()
                obj.element[itotal:itotal2] = eids

            # [f41, f21, f12, f32, f23, f43, f34, f14, kf1,
            #  s12, kf2, s23, kf3, s34, kf4, s41]
            obj.data[obj.itime, itotal:itotal2, :] = floats[:, 1:].copy()
            obj.itotal = itotal2
            obj.ielement = ielement2
        else:
            n = oef_cshear_real_17(op2, data, obj,
                                   nelements, ntotal, dt)

    elif op2.format_code in [2, 3] and op2.num_wide == 33:  # imag
        ntotal = 132 * factor  # 33*4
        nelements = ndata // ntotal

        obj_complex = ComplexCShearForceArray
        auto_return, is_vectorized = op2._create_oes_object4(
            nelements, result_name, slot, obj_complex)
        if auto_return:
            return nelements * ntotal, None, None

        obj = op2.obj
        # if op2.is_debug_file:
        # op2.binary_debug.write('  [cap, element1, element2, ..., cap]\n')
        # op2.binary_debug.write('  cap = %i  # assume 1 cap when there could have been multiple\n' % ndata)
        # op2.binary_debug.write('  #elementi = [eid_device, axial, torque]\n')
        # op2.binary_debug.write('  nelements=%i; nnodes=1 # centroid\n' % nelements)

        if op2.use_vector and is_vectorized and op2.sort_method == 1:
            n = nelements * 4 * op2.num_wide
            itotal = obj.ielement
            ielement2 = obj.itotal + nelements
            itotal2 = ielement2

            floats = np.frombuffer(data, dtype=op2.fdtype8).reshape(nelements, 33).copy()
            obj._times[obj.itime] = dt
            if obj.itime == 0:
                ints = np.frombuffer(data, dtype=op2.idtype8).reshape(nelements, 33)
                eids = ints[:, 0] // 10
                assert eids.min() > 0, eids.min()
                obj.element[itotal:itotal2] = eids

            # [f41r, f21r, f12r, f32r, f23r, f43r, f34r, f14r
            # kf1r, s12r, kf2r, s23r, kf3r, s34r, kf4r, s41r
            # f41i, f21i, f12i, f32i, f23i, f43i, f34i, f14i
            # kf1i, s12i, kf2i, s23i, kf3i, s34i, kf4i, s41i]
            isave1 = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16]
            isave2 = [17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32]
            real_imag = apply_mag_phase(floats, is_magnitude_phase, isave1, isave2)
            obj.data[obj.itime, itotal:itotal2, :] = real_imag

            # obj.data[obj.itime, itotal:itotal2, :] = floats[:, 1:]
            obj.itotal = itotal2
            obj.ielement = ielement2
        else:
            n = oef_cshear_imag_33(op2, data, obj,
                                   nelements, ntotal,
                                   is_magnitude_phase)
    else:  # pragma: no cover
        raise RuntimeError(op2.code_information())
        # msg = op2.code_information()
        # return op2._not_implemented_or_skip(data, ndata, msg), None, None
    return n, nelements, ntotal



def oef_cshear_real_17(op2: OP2, data: bytes,
                       obj: RealCShearForceArray,
                       nelements: int, ntotal: int, dt: Any) -> int:
    n = 0
    s = Struct(op2._endian + mapfmt(op2._analysis_code_fmt + b'16f', op2.size))
    add_sort_x = getattr(obj, 'add_sort' + str(op2.sort_method))
    for unused_i in range(nelements):
        edata = data[n:n+ntotal]

        out = s.unpack(edata)
        if op2.is_debug_file:
            op2.binary_debug.write('OEF_Shear - %s\n' % (str(out)))
        (eid_device,
         f41, f21, f12, f32, f23, f43, f34, f14, kf1,
         s12, kf2, s23, kf3, s34, kf4, s41) = out
        eid, dt = get_eid_dt_from_eid_device(
            eid_device, op2.nonlinear_factor, op2.sort_method)

        #data_in = [eid,
                   #f41, f21, f12, f32, f23, f43, f34,
                   #f14, kf1, s12, kf2, s23, kf3, s34, kf4, s41]
        #print "%s" % (self.get_element_type(op2.element_type)), data_in
        add_sort_x(dt, eid,
                   f41, f21, f12, f32, f23, f43, f34,
                   f14, kf1, s12, kf2, s23, kf3, s34, kf4, s41)
        n += ntotal
    return n

def oef_cshear_imag_33(op2: OP2, data: bytes,
                       obj: ComplexCShearForceArray,
                       nelements: int, ntotal: int,
                       is_magnitude_phase: bool) -> int:
    #ntotal1 = 132 * self.factor
    n = 0
    s = Struct(mapfmt(op2._endian + op2._analysis_code_fmt + b'32f', op2.size))
    add_sort_x = getattr(obj, 'add_sort' + str(op2.sort_method))
    for unused_i in range(nelements):
        edata = data[n:n+ntotal]
        n += ntotal
        out = s.unpack(edata)
        if op2.is_debug_file:
            op2.binary_debug.write('OEF_Shear - %s\n' % (str(out)))
        (eid_device,
         f41r, f21r, f12r, f32r, f23r, f43r, f34r, f14r, # 8
         kf1r, s12r, kf2r, s23r, kf3r, s34r, kf4r, s41r, # 16
         f41i, f21i, f12i, f32i, f23i, f43i, f34i, f14i,
         kf1i, s12i, kf2i, s23i, kf3i, s34i, kf4i, s41i) = out
        if is_magnitude_phase:
            f41 = polar_to_real_imag(f41r, f41i)
            kf1 = polar_to_real_imag(kf1r, kf1i)
            f21 = polar_to_real_imag(f21r, f21i)
            kf2 = polar_to_real_imag(kf2r, kf2i)
            f12 = polar_to_real_imag(f12r, f12i)
            kf3 = polar_to_real_imag(kf3r, kf3i)
            f23 = polar_to_real_imag(f23r, f23i)
            kf4 = polar_to_real_imag(kf4r, kf4i)
            f32 = polar_to_real_imag(f32r, f32i)
            s12 = polar_to_real_imag(s12r, s12i)
            f43 = polar_to_real_imag(f43r, f43i)
            s23 = polar_to_real_imag(s23r, s23i)
            f34 = polar_to_real_imag(f34r, f34i)
            s34 = polar_to_real_imag(s34r, s34i)
            f14 = polar_to_real_imag(f14r, f14i)
            s41 = polar_to_real_imag(s41r, s41i)
        else:
            f41 = complex(f41r, f41i)
            kf1 = complex(kf1r, kf1i)
            f21 = complex(f21r, f21i)
            kf2 = complex(kf2r, kf2i)
            f12 = complex(f12r, f12i)
            kf3 = complex(kf3r, kf3i)
            f23 = complex(f23r, f23i)
            kf4 = complex(kf4r, kf4i)
            f32 = complex(f32r, f32i)
            s12 = complex(s12r, s12i)
            f43 = complex(f43r, f43i)
            s23 = complex(s23r, s23i)
            f34 = complex(f34r, f34i)
            s34 = complex(s34r, s34i)
            f14 = complex(f14r, f14i)
            s41 = complex(s41r, s41i)

        eid, dt = get_eid_dt_from_eid_device(
            eid_device, op2.nonlinear_factor, op2.sort_method)
        add_sort_x(dt, eid,
                   f41, f21, f12, f32, f23, f43, f34, f14,
                   kf1, s12, kf2, s23, kf3, s34, kf4, s41)
    return n
