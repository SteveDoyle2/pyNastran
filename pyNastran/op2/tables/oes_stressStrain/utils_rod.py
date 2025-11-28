from __future__ import annotations
from struct import Struct
from typing import TYPE_CHECKING

from pyNastran.op2.op2_interface.op2_reader import mapfmt
from pyNastran.op2.op2_helper import polar_to_real_imag

from pyNastran.op2.tables.utils import get_eid_dt_from_eid_device
if TYPE_CHECKING:  # pragma: no cover
    from pyNastran.op2.tables.oes_stressStrain.real.oes_rods import RealRodStressArray, RealRodStrainArray
    from pyNastran.op2.tables.oes_stressStrain.complex.oes_rods import ComplexRodStressArray, ComplexRodStrainArray
    from pyNastran.op2.tables.oes_stressStrain.random.oes_rods import RandomRodStressArray, RandomRodStrainArray
    from pyNastran.op2.op2 import OP2


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

        if is_magnitude_phase:
            axial = polar_to_real_imag(axial_real, axial_imag)
            torsion = polar_to_real_imag(torsion_real, torsion_imag)
        else:
            axial = complex(axial_real, axial_imag)
            torsion = complex(torsion_real, torsion_imag)

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
