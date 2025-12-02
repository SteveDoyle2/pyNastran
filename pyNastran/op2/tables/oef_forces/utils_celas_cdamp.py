from __future__ import annotations
from struct import Struct
from typing import Any, TYPE_CHECKING
import numpy as np

from pyNastran.op2.op2_interface.op2_reader import mapfmt
from pyNastran.op2.op2_interface.utils import apply_mag_phase
from pyNastran.op2.tables.utils import get_is_slot_saved, get_eid_dt_from_eid_device
from pyNastran.op2.op2_helper import polar_to_real_imag

from pyNastran.op2.tables.oef_forces.oef_force_objects import (
    RealSpringForceArray,
    RealDamperForceArray,
)
from pyNastran.op2.tables.oef_forces.oef_complex_force_objects import (
    ComplexSpringForceArray,
    ComplexDamperForceArray,
)

if TYPE_CHECKING:  # pragma: no cover
    from pyNastran.op2.op2 import OP2


def oef_celas_cdamp(op2: OP2, data, ndata, dt, is_magnitude_phase,
                     result_type, prefix, postfix):
    """
    11-CELAS1
    12-CELAS2
    13-CELAS3
    14-CELAS4

    20-CDAMP1
    21-CDAMP2
    22-CDAMP3
    23-CDAMP4

    """
    element_type = op2.element_type
    if element_type == 11:
        result_name = prefix + 'celas1_force' + postfix
        obj_real = RealSpringForceArray
        obj_complex = ComplexSpringForceArray
    elif element_type == 12:
        result_name = prefix + 'celas2_force' + postfix
        obj_real = RealSpringForceArray
        obj_complex = ComplexSpringForceArray
    elif element_type == 13:
        result_name = prefix + 'celas3_force' + postfix
        obj_real = RealSpringForceArray
        obj_complex = ComplexSpringForceArray
    elif element_type == 14:
        result_name = prefix + 'celas4_force' + postfix
        obj_real = RealSpringForceArray
        obj_complex = ComplexSpringForceArray

    elif element_type == 20:
        result_name = prefix + 'cdamp1_force' + postfix
        obj_real = RealDamperForceArray
        obj_complex = ComplexDamperForceArray
    elif element_type == 21:
        result_name = prefix + 'cdamp2_force' + postfix
        obj_real = RealDamperForceArray
        obj_complex = ComplexDamperForceArray
    elif element_type == 22:
        result_name = prefix + 'cdamp3_force' + postfix
        obj_real = RealDamperForceArray
        obj_complex = ComplexDamperForceArray
    elif element_type == 23:
        result_name = prefix + 'cdamp4_force' + postfix
        obj_real = RealDamperForceArray
        obj_complex = ComplexDamperForceArray
    else:
        raise NotImplementedError(op2.code_information())

    is_saved, slot = get_is_slot_saved(op2, result_name)
    if not is_saved:
        return ndata, None, None

    factor = op2.factor
    if op2.format_code == 1 and op2.num_wide == 2:  # real
        ntotal = 8 * factor  # 2 * 4
        nelements = ndata // ntotal
        auto_return, is_vectorized = op2._create_oes_object4(
            nelements, result_name, slot, obj_real)
        if auto_return:
            return nelements * op2.num_wide * 4, None, None

        obj = op2.obj
        if op2.use_vector and is_vectorized and op2.sort_method == 1:
            n = nelements * ntotal
            itotal = obj.ielement
            ielement2 = obj.itotal + nelements
            itotal2 = ielement2

            floats = np.frombuffer(data, dtype=op2.fdtype8).reshape(nelements, 2)
            obj._times[obj.itime] = dt
            if obj.itime == 0:
                ints = np.frombuffer(data, dtype=op2.idtype8).reshape(nelements, 2)
                eids = ints[:, 0] // 10
                assert eids.min() > 0, eids.min()
                obj.element[itotal:itotal2] = eids

            # (eid_device, force)
            obj.data[obj.itime, itotal:itotal2, 0] = floats[:, 1].copy()
            obj.itotal = itotal2
            obj.ielement = ielement2
        else:
            n = oef_celas_cdamp_real_2(op2, data, obj,
                                       nelements, ntotal, dt)

    elif op2.format_code in [2, 3] and op2.num_wide == 3:  # imag
        ntotal = 12 * factor  # 3*4
        nelements = ndata // ntotal

        auto_return, is_vectorized = op2._create_oes_object4(
            nelements, result_name, slot, obj_complex)
        if auto_return:
            return nelements * ntotal, None, None

        obj = op2.obj
        if op2.is_debug_file:
            op2.binary_debug.write('  [cap, element1, element2, ..., cap]\n')
            op2.binary_debug.write('  cap = %i  # assume 1 cap when there could have been multiple\n' % ndata)
            op2.binary_debug.write('  #elementi = [eid_device, force]\n')
            op2.binary_debug.write('  nelements=%i; nnodes=1 # centroid\n' % nelements)

        if op2.use_vector and is_vectorized and op2.sort_method == 1:
            n = nelements * ntotal
            itotal = obj.ielement
            ielement2 = obj.itotal + nelements
            itotal2 = ielement2

            floats = np.frombuffer(data, dtype=op2.fdtype8).reshape(nelements, 3).copy()
            obj._times[obj.itime] = dt
            if obj.itime == 0:
                ints = np.frombuffer(data, dtype=op2.idtype8).reshape(nelements, 3)
                eids = ints[:, 0] // 10
                assert eids.min() > 0, eids.min()
                obj.element[itotal:itotal2] = eids

            # [spring_force]
            real_imag = apply_mag_phase(floats, is_magnitude_phase, 1, 2)
            obj.data[obj.itime, itotal:itotal2, 0] = real_imag
            obj.itotal = itotal2
            obj.ielement = ielement2
        else:
            n = oef_celas_cdamp_imag_3(op2, data, obj,
                                       nelements, ntotal,
                                       is_magnitude_phase)
    else:
        raise RuntimeError(op2.code_information())
        # msg = 'OEF: element_name=%s element_type=%s' % (op2.element_name, op2.element_type)
        # msg = op2.code_information()
        # print(msg)
        # return op2._not_implemented_or_skip(data, ndata, msg), None, None
    return n, nelements, ntotal


def oef_celas_cdamp_real_2(op2: OP2, data: bytes,
                           obj: RealSpringForceArray | RealDamperForceArray,
                           nelements: int, ntotal: int, dt: Any) -> int:
    n = 0
    fmt = mapfmt(op2._endian + op2._analysis_code_fmt + b'f', op2.size)
    s = Struct(fmt)  # 2
    add_sort_x = getattr(obj, 'add_sort' + str(op2.sort_method))
    for unused_i in range(nelements):
        edata = data[n:n + ntotal]
        out = s.unpack(edata)
        if op2.is_debug_file:
            op2.binary_debug.write('OEF_SpringDamper - %s\n' % str(out))
        (eid_device, force) = out
        eid, dt = get_eid_dt_from_eid_device(
            eid_device, op2.nonlinear_factor, op2.sort_method)
        add_sort_x(dt, eid, force)
        n += ntotal
    return n


def oef_celas_cdamp_imag_3(op2: OP2, data: bytes,
                           obj: ComplexSpringForceArray | ComplexDamperForceArray,
                           nelements: int, ntotal: int,
                           is_magnitude_phase: bool) -> int:
    n = 0
    fmt = mapfmt(op2._endian + op2._analysis_code_fmt + b'2f', op2.size)
    structi = Struct(fmt)
    add_sort_x = getattr(obj, 'add_sort' + str(op2.sort_method))
    for unused_i in range(nelements):
        edata = data[n:n + ntotal]
        out = structi.unpack(edata)
        if op2.is_debug_file:
            op2.binary_debug.write('OEF_SpringDamper - %s\n' % str(out))
        (eid_device, force_real, force_imag) = out
        eid, dt = get_eid_dt_from_eid_device(
            eid_device, op2.nonlinear_factor, op2.sort_method)
        if is_magnitude_phase:
            force = polar_to_real_imag(force_real, force_imag)
        else:
            force = complex(force_real, force_imag)
        add_sort_x(dt, eid, force)
        n += ntotal
    return n
