from __future__ import annotations
from struct import Struct
from typing import Any, TYPE_CHECKING

from pyNastran.op2.op2_interface.op2_reader import mapfmt
from pyNastran.op2.op2_helper import polar_to_real_imag

from pyNastran.op2.tables.utils import get_eid_dt_from_eid_device

if TYPE_CHECKING:  # pragma: no cover
    from pyNastran.op2.tables.oes_stressStrain.real.oes_springs import (
        RealSpringStressArray, RealSpringStrainArray,
        # RealNonlinearSpringStressArray
    )
    from pyNastran.op2.tables.oes_stressStrain.complex.oes_springs import ComplexSpringStressArray, \
        ComplexSpringStrainArray
    from pyNastran.op2.op2 import OP2


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
