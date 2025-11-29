from __future__ import annotations
from struct import Struct
from typing import TYPE_CHECKING
import numpy as np

from pyNastran.op2.op2_interface.op2_reader import mapfmt
from pyNastran.op2.tables.utils import get_is_slot_saved, get_eid_dt_from_eid_device
from pyNastran.op2.op2_helper import polar_to_real_imag
from pyNastran.op2.op2_interface.utils import apply_mag_phase

from pyNastran.op2.tables.oef_forces.oef_force_objects import (
    RealRodForceArray,
)
from pyNastran.op2.tables.oef_forces.oef_complex_force_objects import (
    ComplexRodForceArray,
)
if TYPE_CHECKING:  # pragma: no cover
    from pyNastran.op2.op2 import OP2


def oef_crod(op2: OP2,
             data, ndata, dt, is_magnitude_phase,
             result_type, prefix, postfix):
    """
    1-CROD
    3-CTUBE
    10-CONROD

    """
    factor = op2.factor
    size = op2.size
    n = 0
    obj_real = RealRodForceArray
    obj_complex = ComplexRodForceArray
    if op2.element_type == 1:  # CROD
        result_name = prefix + 'crod_force' + postfix
    elif op2.element_type == 3:  # CTUBE
        result_name = prefix + 'ctube_force' + postfix
    elif op2.element_type == 10:  # CONROD
        result_name = prefix + 'conrod_force' + postfix
    else:
        raise NotImplementedError(op2.code_information())
        # msg = 'sort1 Type=%s num=%s' % (op2.element_name, op2.element_type)
        # return op2._not_implemented_or_skip(data, ndata, msg)

    if op2._results.is_not_saved(result_name):
        return ndata, None, None
    op2._results._found_result(result_name)

    slot = op2.get_result(result_name)
    if op2.format_code == 1 and op2.num_wide == 3:  # real
        ntotal = 3 * 4 * factor
        nelements = ndata // ntotal

        auto_return, is_vectorized = op2._create_oes_object4(
            nelements, result_name, slot, obj_real)
        if auto_return:
            return nelements * op2.num_wide * 4, None, None
        obj = op2.obj
        if op2.is_debug_file:
            op2.binary_debug.write('  [cap, element1, element2, ..., cap]\n')
            op2.binary_debug.write('  cap = %i  # assume 1 cap when there could have been multiple\n' % ndata)
            op2.binary_debug.write('  #elementi = [eid_device, axial, torque]\n')
            op2.binary_debug.write('  nelements=%i; nnodes=1 # centroid\n' % nelements)

        if op2.use_vector and is_vectorized and op2.sort_method == 1:
            n = nelements * ntotal
            itotal = obj.ielement
            ielement2 = obj.itotal + nelements
            itotal2 = ielement2

            floats = np.frombuffer(data, dtype=op2.fdtype8).reshape(nelements, 3)
            obj._times[obj.itime] = dt
            if obj.itime == 0:
                ints = np.frombuffer(data, dtype=op2.idtype8).reshape(nelements, 3)
                eids = ints[:, 0] // 10
                assert eids.min() > 0, eids.min()
                obj.element[itotal:itotal2] = eids

            # [axial, torsion]
            obj.data[obj.itime, itotal:itotal2, :] = floats[:, 1:].copy()
            obj.itotal = itotal2
            obj.ielement = ielement2
        else:
            n = oef_crod_real_3(op2, data, obj,
                                nelements, ntotal)

    elif op2.format_code in [2, 3] and op2.num_wide == 5:  # imag
        ntotal = 20 * factor  # 5*4
        nelements = ndata // ntotal

        auto_return, is_vectorized = op2._create_oes_object4(
            nelements, result_name, slot, obj_complex)
        if auto_return:
            return nelements * ntotal, None, None

        obj = op2.obj
        if op2.is_debug_file:
            op2.binary_debug.write('  [cap, element1, element2, ..., cap]\n')
            op2.binary_debug.write('  cap = %i  # assume 1 cap when there could have been multiple\n' % ndata)
            op2.binary_debug.write('  #elementi = [eid_device, axial, torque]\n')
            op2.binary_debug.write('  nelements=%i; nnodes=1 # centroid\n' % nelements)

        if op2.use_vector and is_vectorized and op2.sort_method == 1:
            n = nelements * ntotal
            itotal = obj.ielement
            ielement2 = obj.itotal + nelements
            itotal2 = ielement2

            floats = np.frombuffer(data, dtype=op2.fdtype8).reshape(nelements, 5).copy()
            obj._times[obj.itime] = dt
            if obj.itime == 0:
                ints = np.frombuffer(data, dtype=op2.idtype8).reshape(nelements, 5)
                eids = ints[:, 0] // 10
                assert eids.min() > 0, eids.min()
                obj.element[itotal:itotal2] = eids

            # [axial_force, torque]
            # (eid_device, axial_real, torque_real, axial_imag, torque_imag) = out
            real_imag = apply_mag_phase(floats, is_magnitude_phase, [1, 2], [3, 4])
            obj.data[obj.itime, itotal:itotal2, :] = real_imag
            obj.itotal = itotal2
            obj.ielement = ielement2
        else:
            n = oef_crod_imag_5(op2, data, obj,
                                nelements, ntotal,
                                is_magnitude_phase)

    else:  # pragma: no cover
        raise RuntimeError(op2.code_information())
        # msg = op2.code_information()
        # print(msg)
        # return op2._not_implemented_or_skip(data, ndata, msg), None, None
    return n, nelements, ntotal

def oef_crod_real_3(op2: OP2, data: bytes,
                    obj: RealRodForceArray,
                    nelements: int, ntotal: int) -> int:
    n = 0
    fmt = mapfmt(op2._endian + op2._analysis_code_fmt + b'ff', op2.size)  # 3
    s = Struct(fmt)
    add_sort_x = getattr(obj, 'add_sort' + str(op2.sort_method))
    for unused_i in range(nelements):
        edata = data[n:n+ntotal]
        out = s.unpack(edata)
        (eid_device, axial, torque) = out
        eid, dt = get_eid_dt_from_eid_device(
            eid_device, op2.nonlinear_factor, op2.sort_method)
        if op2.is_debug_file:
            op2.binary_debug.write('OEF_Rod - %s\n' % (str(out)))
        add_sort_x(dt, eid, axial, torque)
        n += ntotal
    return n

def oef_crod_imag_5(op2: OP2, data: bytes,
                    obj: ComplexRodForceArray,
                    nelements: int, ntotal: int,
                    is_magnitude_phase: bool) -> int:
    n = 0
    fmt = mapfmt(op2._endian + op2._analysis_code_fmt + b'4f', op2.size)
    s = Struct(fmt)
    add_sort_x = getattr(obj, 'add_sort' + str(op2.sort_method))
    for unused_i in range(nelements):
        edata = data[n:n+ntotal]

        out = s.unpack(edata)
        if op2.is_debug_file:
            op2.binary_debug.write('OEF_CRod - %s\n' % (str(out)))
        (eid_device, axial_real, torque_real, axial_imag, torque_imag) = out
        eid, dt = get_eid_dt_from_eid_device(
            eid_device, op2.nonlinear_factor, op2.sort_method)
        if is_magnitude_phase:
            axial = polar_to_real_imag(axial_real, axial_imag)
            torque = polar_to_real_imag(torque_real, torque_imag)
        else:
            axial = complex(axial_real, axial_imag)
            torque = complex(torque_real, torque_imag)

        add_sort_x(dt, eid, axial, torque)
        n += ntotal
    return n
