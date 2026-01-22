from __future__ import annotations
from struct import Struct
from typing import Any, TYPE_CHECKING
import numpy as np

from pyNastran.op2.op2_interface.utils import mapfmt, real_imag_from_list
from pyNastran.op2.op2_interface.utils import apply_mag_phase
from pyNastran.op2.op2_helper import polar_to_real_imag

from pyNastran.op2.tables.utils import get_is_slot_saved, get_eid_dt_from_eid_device
from pyNastran.op2.tables.oes_stressStrain.utils import obj_set_element

from pyNastran.op2.tables.oes_stressStrain.real.oes_bars import RealBarStressArray, RealBarStrainArray
from pyNastran.op2.tables.oes_stressStrain.complex.oes_bars import ComplexBarStressArray, ComplexBarStrainArray
from pyNastran.op2.tables.oes_stressStrain.random.oes_bars import RandomBarStressArray, RandomBarStrainArray

from pyNastran.op2.tables.oes_stressStrain.real.oes_bars100 import RealBar10NodesStressArray, RealBar10NodesStrainArray


if TYPE_CHECKING:  # pragma: no cover
    from pyNastran.op2.op2 import OP2


def oes_cbar_34(op2: OP2, data: bytes, ndata: int, dt: Any,
                is_magnitude_phase: bool,
                result_type: str, prefix: str, postfix: str) -> tuple[int, int, int]:
    """
    reads stress/strain for element type:
     - 34 : CBAR

    """
    # if isinstance(op2.nonlinear_factor, float):
    # op2.sort_bits[0] = 1 # sort2
    # op2.sort_method = 2

    n = 0
    stress_strain = 'stress' if op2.is_stress else 'strain'
    result_name = f'{prefix}cbar_{stress_strain}{postfix}'

    is_saved, slot = get_is_slot_saved(op2, result_name)
    if not is_saved:
        return ndata, None, None

    factor = op2.factor
    if result_type == 0 and op2.num_wide == 16:  # real
        obj_vector_real = RealBarStressArray if op2.is_stress else RealBarStrainArray

        ntotal = 64 * factor  # 16*4
        nelements = ndata // ntotal
        # print('CBAR nelements =', nelements)

        auto_return, is_vectorized = op2._create_oes_object4(
            nelements, result_name, slot, obj_vector_real)
        if auto_return:
            return ndata, None, None

        if op2.is_debug_file:
            op2.binary_debug.write('  [cap, element1, element2, ..., cap]\n')
            # op2.binary_debug.write('  cap = %i  # assume 1 cap when there could have been multiple\n' % ndata)
            op2.binary_debug.write('  #elementi = [eid_device, s1a, s2a, s3a, s4a, axial, smaxa, smina, MSt,\n')
            op2.binary_debug.write('                           s1b, s2b, s3b, s4b, smaxb, sminb,        MSc]\n')
            op2.binary_debug.write('  nelements=%i; nnodes=1 # centroid\n' % nelements)

        obj = op2.obj
        if op2.use_vector and is_vectorized and op2.sort_method == 1:
            # self.itime = 0
            # self.ielement = 0
            # self.itotal = 0
            # self.ntimes = 0
            # self.nelements = 0
            n = nelements * op2.num_wide * 4

            ielement = obj.ielement
            ielement2 = ielement + nelements
            obj._times[obj.itime] = dt
            obj_set_element(op2, obj, ielement, ielement2, data, nelements)

            floats = np.frombuffer(data, dtype=op2.fdtype8).reshape(nelements, 16)

            # [s1a, s2a, s3a, s4a, axial, smaxa, smina, margin_tension,
            # s1b, s2b, s3b, s4b,        smaxb, sminb, margin_compression]
            obj.data[obj.itime, ielement:ielement2, :] = floats[:, 1:].copy()
            obj.itotal = ielement2
            obj.ielement = ielement2
        else:
            if is_vectorized and op2.use_vector:  # pragma: no cover
                op2.log.debug('vectorize CBAR real SORT%s' % op2.sort_method)
            n = oes_cbar_real_16(op2, data, obj, nelements, ntotal, dt)
    elif result_type == 1 and op2.num_wide == 19:  # imag
        obj_vector_complex = ComplexBarStressArray if op2.is_stress else ComplexBarStrainArray

        ntotal = 76 * factor
        nelements = ndata // ntotal
        auto_return, is_vectorized = op2._create_oes_object4(
            nelements, result_name, slot, obj_vector_complex)
        if auto_return:
            return ndata, None, None

        if op2.is_debug_file:
            op2.binary_debug.write('  [cap, element1, element2, ..., cap]\n')
            # op2.binary_debug.write('  cap = %i  # assume 1 cap when there could have been multiple\n' % ndata)
            op2.binary_debug.write('  #elementi = [eid_device, s1a, s2a, s3a, s4a, axial,\n')
            op2.binary_debug.write('                           s1b, s2b, s3b, s4b]\n')
            op2.binary_debug.write('  nelements=%i; nnodes=1 # centroid\n' % nelements)

        obj = op2.obj
        if op2.use_vector and is_vectorized and op2.sort_method == 1:
            n = nelements * ntotal
            itotal = obj.itotal
            itotal2 = itotal + nelements
            ielement2 = itotal2

            floats = np.frombuffer(data, dtype=op2.fdtype8).reshape(nelements, 19).copy()
            obj._times[obj.itime] = dt
            obj_set_element(op2, obj, itotal, itotal2, data, nelements)

            isave1 = [1, 2, 3, 4, 5, 11, 12, 13, 14]
            isave2 = [6, 7, 8, 9, 10, 15, 16, 17, 18]
            real_imag = apply_mag_phase(floats, is_magnitude_phase, isave1, isave2)
            obj.data[obj.itime, itotal:itotal2, :] = real_imag

            obj.itotal = itotal2
            obj.ielement = ielement2
        else:
            if is_vectorized and op2.use_vector:  # pragma: no cover
                op2.log.debug('vectorize CBAR imag SORT%s' % op2.sort_method)
            n = oes_cbar_complex_19(op2, data, obj, nelements, ntotal, is_magnitude_phase)
    elif result_type == 2 and op2.num_wide == 19:  # random strain?
        raise RuntimeError(op2.code_information())

    elif result_type in [1, 2] and op2.num_wide == 10:  # random
        # random stress/strain per example
        #
        # DMAP says random stress has num_wide=10 and
        # random strain has numwide=19, but it's wrong...maybe???
        #
        # format_code = 1 - NO/RMS (SORT1 regardless of whether this is a SORT2 table or not)
        # format_code = 2 - ATO/PSD/CRM (actually SORT2)
        #
        element_id = op2.nonlinear_factor
        obj_vector_random = RandomBarStressArray if op2.is_stress else RandomBarStrainArray
        op2.data_code['nonlinear_factor'] = element_id

        ntotal = 10 * op2.size
        nelements = ndata // ntotal
        # print(f'CBAR* nelements={nelements}')
        auto_return, is_vectorized = op2._create_oes_object4(
            nelements, result_name, slot, obj_vector_random)
        if auto_return:
            return ndata, None, None

        if op2.is_debug_file:
            op2.binary_debug.write('  [cap, element1, element2, ..., cap]\n')
            # op2.binary_debug.write('  cap = %i  # assume 1 cap when there could have been multiple\n' % ndata)
            op2.binary_debug.write('  #elementi = [eid_device, s1a, s2a, s3a, s4a, axial,\n')
            op2.binary_debug.write('                           s1b, s2b, s3b, s4b]\n')
            op2.binary_debug.write('  nelements=%i; nnodes=1 # centroid\n' % nelements)

        obj = op2.obj
        if op2.use_vector and is_vectorized and 0:  # pragma: no cover
            # self.itime = 0
            # self.ielement = 0
            # self.itotal = 0
            # self.ntimes = 0
            # self.nelements = 0
            n = nelements * ntotal

            ielement = obj.ielement
            ielement2 = ielement + nelements
            obj._times[obj.itime] = dt
            obj_set_element(op2, obj, ielement, ielement2, data, nelements)

            floats = np.frombuffer(data, dtype=op2.fdtype).reshape(nelements, 10)

            # [s1a, s2a, s3a, s4a, axial,
            # s1b, s2b, s3b, s4b]
            obj.data[obj.itime, ielement:ielement2, :] = floats[:, 1:].copy()
            obj.itotal = ielement2
            obj.ielement = ielement2
        else:
            if is_vectorized and op2.use_vector and obj.itime == 0:  # pragma: no cover
                op2.log.debug('vectorize CBAR random SORT%s' % op2.sort_method)
            n = oes_cbar_random_10(op2, data, obj,
                                   nelements, ntotal)
    else:  # pragma: no cover
        raise RuntimeError(op2.code_information())
    return n, nelements, ntotal


def oes_cbar_real_16(op2: OP2, data: bytes,
                     obj: RealBarStressArray | RealBarStrainArray,
                     nelements: int, ntotal: int, dt: Any) -> int:
    n = 0
    fmt = mapfmt(op2._endian + op2._analysis_code_fmt + b'15f', op2.size)
    struct1 = Struct(fmt)
    add_sort_x = getattr(obj, 'add_new_eid_sort' + str(op2.sort_method))
    for i in range(nelements):
        edata = data[n:n+ntotal]
        out = struct1.unpack(edata)
        (eid_device,
         s1a, s2a, s3a, s4a, axial, smaxa, smina, margin_tension,
         s1b, s2b, s3b, s4b, smaxb, sminb, margin_compression) = out
        eid, dt = get_eid_dt_from_eid_device(
            eid_device, op2.nonlinear_factor, op2.sort_method)
        if op2.is_debug_file:
            op2.binary_debug.write('  eid=%i; C%i=[%s]\n' % (
                eid, i, ', '.join(['%r' % di for di in out])))
        n += ntotal
        add_sort_x(
            dt, eid,
            s1a, s2a, s3a, s4a, axial, smaxa, smina, margin_tension,
            s1b, s2b, s3b, s4b, smaxb, sminb, margin_compression)

    if op2.is_sort2:
        #print(add_sort_x)
        #print(''.join(obj.get_stats()))
        #print(f'{self.table_name} sort_method={op2.sort_method}', obj._times)
        assert len(np.unique(obj._times)) == len(obj._times), obj._times.tolist()
    return n


def oes_cbar_complex_19(op2: OP2,
                        data: bytes,
                        obj: ComplexBarStressArray | ComplexBarStrainArray,
                        nelements: int, ntotal: int,
                        is_magnitude_phase: bool) -> int:
    n = 0
    struct1 = Struct(mapfmt(op2._endian + op2._analysis_code_fmt + b'18f', op2.size))
    for i in range(nelements):
        edata = data[n:n+ntotal]
        n += ntotal
        out = struct1.unpack(edata)
        (eid_device,
         s1ar, s2ar, s3ar, s4ar, axialr,
         s1ai, s2ai, s3ai, s4ai, axiali,
         s1br, s2br, s3br, s4br,
         s1bi, s2bi, s3bi, s4bi) = out

        eid, dt = get_eid_dt_from_eid_device(
            eid_device, op2.nonlinear_factor, op2.sort_method)
        if op2.is_debug_file:
            op2.binary_debug.write('  eid=%i; C%i=[%s]\n' % (
                eid, i, ', '.join(['%r' % di for di in out])))
        (s1a, s2a, s3a, s4a, axial,
         s1b, s2b, s3b, s4b) = real_imag_from_list([
            s1ar, s2ar, s3ar, s4ar, axialr,
            s1br, s2br, s3br, s4br,
            s1ai, s2ai, s3ai, s4ai, axiali,
            s1bi, s2bi, s3bi, s4bi], is_magnitude_phase)
        obj.add_new_eid_sort1(dt, eid,
                              s1a, s2a, s3a, s4a, axial,
                              s1b, s2b, s3b, s4b)
    return n


def oes_cbar_random_10(op2: OP2, data: bytes,
                       obj: RandomBarStressArray | RandomBarStrainArray,
                       nelements: int, ntotal: int) -> int:
    n = 0
    #print(op2.code_information())
    #print('op2._analysis_code_fmt =', op2._analysis_code_fmt)
    struct1 = Struct(op2._endian + op2._analysis_code_fmt + b'9f')
    add_sort_x = getattr(obj, 'add_sort' + str(op2.sort_method))

    #self.log.info('op2.nonlinear_factor = %s' % op2.nonlinear_factor)
    #assert op2.sort_method == 2, op2.code_information()
    #if sort_method == 2:
        #obj.node_id = 42
    nonlinear_factor = op2.nonlinear_factor
    #print(f'CBAR: nelements={nelements}')
    for i in range(nelements):
        edata = data[n:n+ntotal]
        n += ntotal

        out = struct1.unpack(edata)
        (eid_device,
         s1a, s2a, s3a, s4a, axial,
         s1b, s2b, s3b, s4b) = out
        eid, dt = get_eid_dt_from_eid_device(
            eid_device, nonlinear_factor, op2.sort_method)
        #print(f'  eid_device={eid_device} eid={eid} dt={nonlinear_factor} nf={nonlinear_factor} -> {obj.data.shape}')
        #continue
        #print('  eid=%i; C%i=[%s]\n' % (eid, i, ', '.join(['%r' % di for di in out])))
        if op2.table_name_str == 'OESXRMS1':
            #assert sort_method == 2
            assert op2.sort_method == 1, op2.code_information()

        if op2.is_debug_file:
            op2.binary_debug.write('  eid=%i; C%i=[%s]\n' % (
                eid, i, ', '.join(['%r' % di for di in out])))

        assert eid > 0, "dt=%s eid=%s" % (dt, eid)
        add_sort_x(dt, eid,
                   s1a, s2a, s3a, s4a, axial,
                   s1b, s2b, s3b, s4b)
    return n


def oes_cbar_100(op2: OP2, data: bytes, ndata: int, dt, is_magnitude_phase: bool,
                 result_type: str, prefix: str, postfix: str):
    """
    reads stress/strain for element type:
     - 100 : BARS
    """
    stress_strain = 'stress' if op2.is_stress else 'strain'
    result_name = f'{prefix}cbar_{stress_strain}_10nodes{postfix}'

    is_saved, slot = get_is_slot_saved(op2, result_name)
    if not is_saved:
        return ndata, None, None

    factor = op2.factor
    if result_type == 0 and op2.num_wide == 10:  # real
        if op2.is_stress:
            obj_vector_real = RealBar10NodesStressArray
        else:
            obj_vector_real = RealBar10NodesStrainArray

        ntotal = 10 * 4 * factor
        nelements = ndata // ntotal

        auto_return, is_vectorized = op2._create_oes_object4(
            nelements, result_name, slot, obj_vector_real)
        if auto_return:
            return ndata, None, None

        if op2.is_debug_file:
            op2.binary_debug.write('  [cap, element1, element2, ..., cap]\n')
            # op2.binary_debug.write('  cap = %i  # assume 1 cap when there could have been multiple\n' % ndata)
            op2.binary_debug.write('  #elementi = [eid_device, sd, sxc, sxd, sxe, sxf, axial, smax, smin, MS]\n')
            op2.binary_debug.write('  nelements=%i; nnodes=1 # centroid\n' % nelements)
        obj = op2.obj

        if op2.use_vector and is_vectorized and op2.sort_method == 1:
            # self.itime = 0
            # self.ielement = 0
            # self.itotal = 0
            # self.ntimes = 0
            # self.nelements = 0
            n = nelements * ntotal

            istart = obj.itotal
            iend = istart + nelements
            obj._times[obj.itime] = dt

            obj_set_element(op2, obj, istart, iend, data, nelements)

            floats = np.frombuffer(data, dtype=op2.fdtype8).reshape(nelements, 10)
            # [sd, sxc, sxd, sxe, sxf, axial, smax, smin, MS]
            obj.data[obj.itime, istart:iend, :] = floats[:, 1:].copy()
        else:
            n = oes_cbar100_real_10(op2, data, obj, nelements, ntotal, dt)

    elif result_type == 1 and op2.num_wide == 16:  # complex
        msg = op2.code_information()
        return op2._not_implemented_or_skip(data, ndata, msg), None, None
    else:  # pragma: no cover
        raise RuntimeError(op2.code_information())
    return n, nelements, ntotal


def oes_cbar100_real_10(op2: OP2, data: bytes, obj, nelements: int, ntotal: int, dt: Any) -> int:
    n = 0
    struct1 = Struct(op2._endian + mapfmt(op2._analysis_code_fmt + b'9f', op2.size))
    for i in range(nelements):
        edata = data[n:n+ntotal]
        out = struct1.unpack(edata)
        (eid_device, sd, sxc, sxd, sxe, sxf, axial, smax, smin, MS) = out
        eid, dt = get_eid_dt_from_eid_device(
            eid_device, op2.nonlinear_factor, op2.sort_method)
        if op2.is_debug_file:
            op2.binary_debug.write('  eid=%i; C%i=[%s]\n' % (
                eid, i, ', '.join(['%r' % di for di in out])))
        n += ntotal
        obj.add_new_eid_sort1(op2.element_name, dt, eid,
                              sd, sxc, sxd, sxe, sxf, axial, smax, smin, MS)
    return n
