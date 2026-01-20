from __future__ import annotations
from struct import Struct
from typing import TYPE_CHECKING
import numpy as np

from pyNastran.op2.op2_interface.utils import mapfmt, real_imag_from_list, apply_mag_phase
from pyNastran.op2.op2_helper import polar_to_real_imag

from pyNastran.op2.tables.utils import get_is_slot_saved, get_eid_dt_from_eid_device
from pyNastran.op2.tables.oes_stressStrain.utils import obj_set_element
from pyNastran.op2.tables.oes_stressStrain.real.oes_rods import RealRodStressArray, RealRodStrainArray
from pyNastran.op2.tables.oes_stressStrain.complex.oes_rods import ComplexRodStressArray, ComplexRodStrainArray
from pyNastran.op2.tables.oes_stressStrain.random.oes_rods import RandomRodStressArray, RandomRodStrainArray


if TYPE_CHECKING:  # pragma: no cover
    from pyNastran.op2.op2 import OP2


def oes_crod(op2: OP2, data, ndata, dt, is_magnitude_phase,
             result_type, prefix, postfix):
    """
    reads stress/strain for element type:
     - 1 : CROD
     - 3 : CTUBE
     - 10 : CONROD

    """
    n = 0
    if op2.element_type == 1:  # CROD
        etype = 'crod'
    elif op2.element_type == 3:  # CTUBE
        etype = 'ctube'
    elif op2.element_type == 10:  # CONROD
        etype = 'conrod'
    else:  # pragma: no cover
        raise ValueError(op2.code_information())

    if op2.is_stress:
        obj_vector_real = RealRodStressArray
        obj_vector_complex = ComplexRodStressArray
        obj_vector_random = RandomRodStressArray
        result_name = f'{prefix}{etype}_stress{postfix}'
    else:
        obj_vector_real = RealRodStrainArray
        obj_vector_complex = ComplexRodStrainArray
        obj_vector_random = RandomRodStrainArray
        result_name = f'{prefix}{etype}_strain{postfix}'

    is_saved, slot = get_is_slot_saved(op2, result_name)
    if not is_saved:
        return ndata, None, None

    # result_name, unused_is_random = self._apply_oes_ato_crm_psd_rms_no(result_name)
    factor = op2.factor
    if result_type == 0 and op2.num_wide == 5:  # real
        ntotal = 5 * 4 * factor
        nelements = ndata // ntotal

        auto_return, is_vectorized = op2._create_oes_object4(
            nelements, result_name, slot, obj_vector_real)
        if auto_return:
            # assert ntotal == op2.num_wide * 4
            return nelements * ntotal, None, None

        obj = op2.obj
        if op2.use_vector and is_vectorized and op2.sort_method == 1:
            n = nelements * 4 * op2.num_wide
            itotal = obj.ielement
            ielement2 = obj.itotal + nelements
            itotal2 = ielement2

            floats = np.frombuffer(data, dtype=op2.fdtype8).reshape(nelements, 5)
            obj._times[obj.itime] = dt
            obj_set_element(op2, obj, itotal, itotal2, data, nelements)

            # [axial, torsion, SMa, SMt]
            obj.data[obj.itime, itotal:itotal2, :] = floats[:, 1:].copy()
            obj.itotal = itotal2
            obj.ielement = ielement2
        else:
            if is_vectorized and op2.use_vector:  # pragma: no cover
                op2.log.debug('vectorize CROD real SORT%s' % op2.sort_method)
            if op2.is_debug_file:
                op2.binary_debug.write('  [cap, element1, element2, ..., cap]\n')
                op2.binary_debug.write('  cap = %i  # assume 1 cap when there could have been multiple\n' % ndata)
                op2.binary_debug.write('  #elementi = [eid_device, axial, axial_margin, torsion, torsion_margin]\n')
                op2.binary_debug.write('  nelements=%i; nnodes=1 # centroid\n' % nelements)
            n = oes_crod_real_5(op2, data, obj, nelements, ntotal, dt)

    elif result_type == 1 and op2.num_wide == 5:  # imag
        ntotal = 20 * factor
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

            floats = np.frombuffer(data, dtype=op2.fdtype8).reshape(nelements, 5)
            obj._times[obj.itime] = dt
            obj_set_element(op2, obj, itotal, itotal2, data, nelements)

            real_imag = apply_mag_phase(floats, is_magnitude_phase, [1, 3], [2, 4])
            obj.data[obj.itime, itotal:itotal2, :] = real_imag

            obj.itotal = itotal2
            obj.ielement = ielement2
        else:
            if is_vectorized and op2.use_vector:  # pragma: no cover
                op2.log.debug('vectorize CROD imag SORT%s' % op2.sort_method)
            n = oes_crod_complex_5(op2, data, obj, nelements, ntotal, dt, is_magnitude_phase)

    # elif op2.format_code in [2, 3] and op2.num_wide == 8:  # is this imag ???
    # ntotal = 32
    # s = self.self.struct_i
    # nelements = ndata // ntotal
    # for i in range(nelements):
    # edata = data[n:n + 4]
    # eid_device, = s.unpack(edata)
    # assert eid > 0, eid
    # n += ntotal
    elif result_type == 2 and op2.num_wide == 3:  # random
        ntotal = 3 * 4 * factor
        nelements = ndata // ntotal

        auto_return, is_vectorized = op2._create_oes_object4(
            nelements, result_name, slot, obj_vector_random)
        if auto_return:
            return nelements * ntotal, None, None

        obj = op2.obj
        if op2.use_vector and is_vectorized and op2.sort_method == 1:
            if is_vectorized and op2.use_vector:  # pragma: no cover
                op2.log.debug('vectorize CROD random SORT%s' % op2.sort_method)
            n = nelements * ntotal
            itotal = obj.ielement
            ielement2 = obj.itotal + nelements
            itotal2 = ielement2

            floats = np.frombuffer(data, dtype=op2.fdtype).reshape(nelements, 3)
            obj._times[obj.itime] = dt
            obj_set_element(op2, obj, itotal, itotal2, data, nelements)

            # [axial, torsion, SMa, SMt]
            obj.data[obj.itime, itotal:itotal2, :] = floats[:, 1:].copy()
            obj.itotal = itotal2
            obj.ielement = ielement2
        else:
            n = oes_crod_random_3(op2, data, ndata, obj, nelements, ntotal)

    else:  # pragma: no cover
        raise RuntimeError(op2.code_information())
    assert op2.num_wide * 4 * factor == ntotal, f'numwide*4={op2.num_wide * 4} ntotal={ntotal} element_name={op2.element_name!r}\n{op2.code_information()}'
    return n, nelements, ntotal


def oes_crod_real_5(op2: OP2, data: bytes,
                    obj: RealRodStressArray | RealRodStrainArray,
                    nelements: int, ntotal: int, dt) -> int:
    n = 0
    struct1 = Struct(op2._endian + mapfmt(op2._analysis_code_fmt + b'4f', op2.size))
    add_sort_x = getattr(obj, 'add_sort' + str(op2.sort_method))
    for unused_i in range(nelements):
        edata = data[n:n+ntotal]
        out = struct1.unpack(edata)
        (eid_device, axial, axial_margin, torsion, torsion_margin) = out
        eid, dt = get_eid_dt_from_eid_device(
            eid_device, op2.nonlinear_factor, op2.sort_method)

        if op2.is_debug_file:
            op2.binary_debug.write('  eid=%i; C=[%s]\n' % (
                eid, ', '.join(['%r' % di for di in out])))
        add_sort_x(dt, eid, axial, axial_margin, torsion, torsion_margin)
        n += ntotal

    return n


def oes_crod_complex_5(op2: OP2, data: bytes,
                       obj: ComplexRodStressArray | ComplexRodStrainArray,
                       nelements: int, ntotal: int, dt,
                       is_magnitude_phase: bool) -> int:
    n = 0
    fmt = mapfmt(op2._endian + op2._analysis_code_fmt + b'4f', op2.size)
    struct1 = Struct(fmt)
    add_sort_x = getattr(obj, 'add_sort' + str(op2.sort_method))
    for unused_i in range(nelements):
        edata = data[n:n + ntotal]
        out = struct1.unpack(edata)
        (eid_device, axial_real, axial_imag, torsion_real, torsion_imag) = out
        eid, dt = get_eid_dt_from_eid_device(
            eid_device, op2.nonlinear_factor, op2.sort_method)

        axial, torsion = real_imag_from_list([
            axial_real, torsion_real,
            axial_imag, torsion_imag], is_magnitude_phase)
        add_sort_x(dt, eid, axial, torsion)
        n += ntotal
    return n


def oes_crod_random_3(op2: OP2, data: bytes, ndata: int,
             obj: RandomRodStressArray | RandomRodStrainArray,
             nelements: int, ntotal: int) -> int:
    n = 0
    if op2.is_debug_file:
        op2.binary_debug.write('  [cap, element1, element2, ..., cap]\n')
        op2.binary_debug.write('  cap = %i  # assume 1 cap when there could have been multiple\n' % ndata)
        op2.binary_debug.write('  #elementi = [eid_device, axial, torsion]\n')
        op2.binary_debug.write('  nelements=%i; nnodes=1 # centroid\n' % nelements)

    struct1 = Struct(op2._endian + op2._analysis_code_fmt + b'2f')
    add_sort_x = getattr(obj, 'add_sort' + str(op2.sort_method))
    for unused_i in range(nelements):
        edata = data[n:n+ntotal]
        out = struct1.unpack(edata)
        (eid_device, axial, torsion) = out
        eid, dt = get_eid_dt_from_eid_device(
            eid_device, op2.nonlinear_factor, op2.sort_method)
        if op2.is_debug_file:
            op2.binary_debug.write('  eid=%i; C=[%s]\n' % (
                eid, ', '.join(['%r' % di for di in out])))
        add_sort_x(dt, eid, axial, torsion)
        n += ntotal
    return n
