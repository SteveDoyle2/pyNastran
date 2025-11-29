from __future__ import annotations
from struct import Struct
from typing import Any, TYPE_CHECKING
import numpy as np

from pyNastran.op2.op2_interface.op2_reader import mapfmt
from pyNastran.op2.op2_helper import polar_to_real_imag

from pyNastran.op2.tables.oes_stressStrain.utils import obj_set_element
from pyNastran.op2.tables.utils import get_is_slot_saved, get_eid_dt_from_eid_device

from pyNastran.op2.tables.oes_stressStrain.real.oes_springs import (
    RealSpringStressArray, RealSpringStrainArray,
    # RealNonlinearSpringStressArray
)
from pyNastran.op2.tables.oes_stressStrain.complex.oes_springs import ComplexSpringStressArray, \
    ComplexSpringStrainArray


if TYPE_CHECKING:  # pragma: no cover
    from pyNastran.op2.op2 import OP2


def oes_celas(op2: OP2, data, ndata, dt, is_magnitude_phase,
              result_type, prefix, postfix):
    """
    reads stress/strain for element type:
     - 11 : CELAS1
     - 12 : CELAS2
     - 13 : CELAS3
     - 14 : CELAS4
    """
    # n = 0
    factor = op2.factor
    if op2.element_type == 11:
        etype = 'celas1'
    elif op2.element_type == 12:
        etype = 'celas2'
    elif op2.element_type == 13:
        etype = 'celas3'
    elif op2.element_type == 14:
        etype = 'celas4'
    else:  # pragma: no cover
        raise RuntimeError(op2.element_type)

    if op2.is_stress:
        if prefix == '' and postfix == '':
            prefix = 'stress.'
        obj_real = RealSpringStressArray
        obj_complex = ComplexSpringStressArray
        result_name = f'{prefix}{etype}_stress{postfix}'
    else:
        if prefix == '' and postfix == '':
            prefix = 'strain.'
        obj_real = RealSpringStrainArray
        obj_complex = ComplexSpringStrainArray
        result_name = f'{prefix}{etype}_strain{postfix}'

    is_saved, slot = get_is_slot_saved(op2, result_name)
    if not is_saved:
        return ndata, None, None
    log = op2.log

    if op2.format_code == 1 and op2.num_wide == 2:  # real
        ntotal = 8 * factor  # 2 * 4
        nelements = ndata // ntotal
        auto_return, is_vectorized = op2._create_oes_object4(
            nelements, result_name, slot, obj_real)
        if auto_return:
            # assert ntotal == op2.num_wide * 4
            return nelements * ntotal, None, None

        obj = op2.obj
        if op2.use_vector and is_vectorized and op2.sort_method == 1:
            n = nelements * 4 * op2.num_wide
            itotal = obj.ielement
            ielement2 = obj.itotal + nelements
            itotal2 = ielement2

            floats = np.frombuffer(data, dtype=op2.fdtype8).reshape(nelements, 2)
            obj._times[obj.itime] = dt
            obj_set_element(op2, obj, itotal, itotal2, data, nelements)

            # (eid_device, stress)
            obj.data[obj.itime, itotal:itotal2, 0] = floats[:, 1].copy()
            obj.itotal = itotal2
            obj.ielement = ielement2
        else:
            if is_vectorized and op2.use_vector:  # pragma: no cover
                log.debug('vectorize CELASx real SORT%s' % op2.sort_method)
            n = oes_celas_real_2(op2, data, obj, nelements, ntotal, dt)

    elif op2.format_code in [2, 3] and op2.num_wide == 3:  # imag
        ntotal = 12 * factor
        nelements = ndata // ntotal

        nelements = ndata // ntotal
        auto_return, is_vectorized = op2._create_oes_object4(
            nelements, result_name, slot, obj_complex)
        if auto_return:
            return nelements * ntotal, None, None

        obj = op2.obj
        assert obj is not None, op2.code_information()
        if op2.use_vector and is_vectorized and op2.sort_method == 1:
            n = nelements * ntotal
            itotal = obj.ielement
            ielement2 = obj.itotal + nelements
            itotal2 = ielement2

            floats = np.frombuffer(data, dtype=op2.fdtype8).reshape(nelements, 3).copy()
            obj._times[obj.itime] = dt
            obj_set_element(op2, obj, itotal, itotal2, data, nelements)

            if is_magnitude_phase:
                mag = floats[:, 1]
                phase = floats[:, 2]
                rtheta = np.radians(phase)
                real_imag = mag * (np.cos(rtheta) + 1.j * np.sin(rtheta))
            else:
                real = floats[:, 1]
                imag = floats[:, 2]
                real_imag = real + 1.j * imag
            obj.data[obj.itime, itotal:itotal2, 0] = real_imag

            obj.itotal = itotal2
            obj.ielement = ielement2
        else:
            if is_vectorized and op2.use_vector:  # pragma: no cover
                log.debug('vectorize CELASx imag SORT%s' % op2.sort_method)

            n = oes_celas_complex_3(op2, data, obj,
                                    nelements, ntotal,
                                    dt, is_magnitude_phase)
    elif op2.format_code == 1 and op2.num_wide == 3:  # random
        raise RuntimeError(op2.code_information())
        # msg = op2.code_information()
        # return op2._not_implemented_or_skip(data, ndata, msg)
    else:  # pragma: no cover
        raise RuntimeError(op2.code_information())
    return n, nelements, ntotal


def oes_celas_real_2(op2: OP2, data: bytes,
                     obj: RealSpringStressArray | RealSpringStrainArray,
                     nelements: int, ntotal: int, dt: Any) -> int:
    n = 0
    fmt1 = mapfmt(op2._endian + op2._analysis_code_fmt + b'f', op2.size)
    struct1 = Struct(fmt1)
    add_sort_x = getattr(obj, 'add_sort' + str(op2.sort_method))
    for i in range(nelements):
        edata = data[n:n+ntotal]
        out = struct1.unpack(edata)
        (eid_device, ox) = out
        eid, dt = get_eid_dt_from_eid_device(
            eid_device, op2.nonlinear_factor, op2.sort_method)
        if eid <= 0: # pragma: no cover
            msg = 'table_name=%s sort_method=%s eid_device=%s nonlinear_factor=%s' % (
                op2.table_name_str, op2.sort_method,
                eid_device, op2.nonlinear_factor)
            raise RuntimeError(msg)
        if op2.is_debug_file:
            op2.binary_debug.write('  eid=%i result%i=[%i, %f]\n' % (
                eid, i, eid_device, ox))
        add_sort_x(dt, eid, ox)
        n += ntotal
    return n


def oes_celas_complex_3(op2: OP2, data: bytes, obj: ComplexSpringStressArray | ComplexSpringStrainArray,
                        nelements: int, ntotal: int,
                        dt, is_magnitude_phase: bool):
    n = 0
    struct1 = Struct(op2._endian + mapfmt(op2._analysis_code_fmt + b'2f', op2.size))
    add_sort_x = getattr(obj, 'add_sort' + str(op2.sort_method))
    for i in range(nelements):
        edata = data[n:n + ntotal]
        out = struct1.unpack(edata)
        (eid_device, axial_real, axial_imag) = out
        eid, dt = get_eid_dt_from_eid_device(
            eid_device, op2.nonlinear_factor, op2.sort_method)

        if is_magnitude_phase:
            axial = polar_to_real_imag(axial_real, axial_imag)
        else:
            axial = complex(axial_real, axial_imag)

        if op2.is_debug_file:
            op2.binary_debug.write('  eid=%i result%i=[%i, %f, %f]\n' % (
                eid, i, eid_device, axial_real, axial_imag))
        add_sort_x(dt, eid, axial)
        n += ntotal
    return n
