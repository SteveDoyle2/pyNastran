from __future__ import annotations
from struct import Struct
from typing import Any, TYPE_CHECKING

from pyNastran.op2.op2_interface.op2_reader import mapfmt
# from pyNastran.op2.op2_helper import polar_to_real_imag

from pyNastran.op2.tables.utils import get_eid_dt_from_eid_device
if TYPE_CHECKING:  # pragma: no cover
    from pyNastran.op2.tables.oes_stressStrain.real.oes_composite_plates import RealCompositePlateStressArray, \
        RealCompositePlateStrainArray
    from pyNastran.op2.tables.oes_stressStrain.complex.oes_composite_plates import (
        ComplexLayeredCompositeStrainArray, ComplexLayeredCompositeStressArray,
        ComplexLayeredCompositesArray)
    from pyNastran.op2.op2 import OP2


def oes_comp_shell_real_11(op2: OP2, data: bytes, ndata: int,
                           obj: RealCompositePlateStressArray | RealCompositePlateStrainArray,
                           ntotal: int, nelements: int, etype: str, dt: Any) -> int:
    n = 0
    eid_old = 0
    struct1 = Struct(op2._endian + mapfmt(op2._analysis_code_fmt + b'i9f', op2.size)) # 11
    sort_method = op2.sort_method
    add_eid_sort_x = getattr(obj, 'add_eid_sort' + str(op2.sort_method))
    add_sort_x = getattr(obj, 'add_sort' + str(op2.sort_method))

    for unused_i in range(nelements):
        edata = data[n:n+ntotal]  # 4*11
        out = struct1.unpack(edata)
        (eid_device, layer, o1, o2, t12, t1z, t2z, angle, major, minor, ovm) = out
        eid, dt = get_eid_dt_from_eid_device(
            eid_device, op2.nonlinear_factor, sort_method)

        if op2.is_debug_file:
            op2.binary_debug.write('  eid=%i; layer=%i; C=[%s]\n' % (eid, layer, ', '.join(['%r' % di for di in out])))

        if eid != eid_old:
            # originally initialized to None, the buffer doesnt reset it, so it is the old value
            add_eid_sort_x(etype, dt, eid, layer, o1, o2,
                           t12, t1z, t2z, angle, major, minor, ovm)
        else:
            add_sort_x(dt, eid, layer, o1, o2,
                       t12, t1z, t2z, angle, major, minor, ovm)
        eid_old = eid
        n += ntotal
    return n


def oes_shell_composite_complex_11(op2: OP2,
                                   data: bytes,
                                   obj: ComplexLayeredCompositeStressArray | ComplexLayeredCompositeStrainArray,
                                   ntotal: int, nelements: int, sort_method: int,
                                   dt: Any, is_magnitude_phase: bool) -> int:
    """OESCP, OESTRCP"""
    n = 0
    struct1 = Struct(op2._endian + op2._analysis_code_fmt + b'i9f')
    add_sort_x = getattr(obj, 'add_sort' + str(op2.sort_method))
    for unused_i in range(nelements):
        edata = data[n:n+ntotal]
        out = struct1.unpack(edata)

        (eid_device, ply_id, oxx, oyy, txy, txz, tyz, angle, omax, omin, max_shear) = out
        eid, dt = get_eid_dt_from_eid_device(
            eid_device, op2.nonlinear_factor, sort_method)
        #print(eid, out)

        #if op2.is_debug_file:
            #op2.binary_debug.write('%s-%s - (%s) + %s\n' % (op2.element_name, op2.element_type, eid_device, str(out)))
        #add_sort_x(dt, eid, theory, lamid, fp, fm, fb, fmax, fflag)
        add_sort_x(dt, eid, ply_id, oxx, oyy, txy, txz, tyz, angle, omax, omin, max_shear)
        n += ntotal
    return n


def oes_shell_composite_complex_13(op2: OP2,
                                   data: bytes,
                                   obj: ComplexLayeredCompositesArray,
                                   ntotal: int, nelements: int, sort_method: int,
                                   dt: Any, is_magnitude_phase: bool) -> int:
    """OESVM1C, OSTRVM1C"""
    # OESCP - STRAINS IN LAYERED COMPOSITE ELEMENTS (QUAD4)
    n = 0
    struct1 = Struct(op2._endian + op2._analysis_code_fmt + b'i9f ff')
    add_sort_x = getattr(obj, 'add_sort' + str(op2.sort_method))
    for unused_i in range(nelements):
        edata = data[n:n+ntotal]
        out = struct1.unpack(edata)

        (eid_device, ply_id,
         o1a, o2a, t12a, o1za, o2za,
         o1b, o2b, t12b, o1zb, e2zb, ovm,) = out
        eid, dt = get_eid_dt_from_eid_device(
            eid_device, op2.nonlinear_factor, sort_method)
        #print(eid, out)

        #print('%s-%s - (%s) + %s\n' % (op2.element_name, op2.element_type, eid_device, str(out)))
        #if op2.is_debug_file:
            #op2.binary_debug.write('%s-%s - (%s) + %s\n' % (op2.element_name, op2.element_type, eid_device, str(out)))
        add_sort_x(dt, eid, ply_id,
                   o1a, o2a, t12a, o1za, o2za,
                   o1b, o2b, t12b, o1zb, e2zb, ovm)
        n += ntotal
    return n
