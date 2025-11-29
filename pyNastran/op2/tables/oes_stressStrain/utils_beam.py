from __future__ import annotations
from struct import Struct
from typing import Any, TYPE_CHECKING

from pyNastran.op2.op2_interface.op2_reader import mapfmt
from pyNastran.op2.op2_helper import polar_to_real_imag

from pyNastran.op2.tables.utils import get_eid_dt_from_eid_device

if TYPE_CHECKING:  # pragma: no cover
    from pyNastran.op2.tables.oes_stressStrain.real.oes_beams import (
        RealBeamStressArray, RealBeamStrainArray,
        # RealNonlinearBeamStressArray
    )
    from pyNastran.op2.tables.oes_stressStrain.real.oes_beams import (
        RealBeamStressArray, RealBeamStrainArray,
        # RealNonlinearBeamStressArray
    )
    from pyNastran.op2.tables.oes_stressStrain.complex.oes_beams import ComplexBeamStressArray, ComplexBeamStrainArray
    from pyNastran.op2.tables.oes_stressStrain.random.oes_beams import RandomBeamStressArray, RandomBeamStrainArray
    from pyNastran.op2.op2 import OP2


def oes_cbeam_real_111(op2: OP2, data: bytes,
                       obj: RealBeamStressArray | RealBeamStrainArray,
                       nelements: int, dt: Any) -> int:
    n = 0
    nnodes = 10  # 11-1
    n1 = 44 * op2.factor
    n2 = 40 * op2.factor
    fmt1 = mapfmt(op2._endian + op2._analysis_code_fmt + b'i9f', op2.size)
    fmt2 = mapfmt(op2._endian + b'i9f', op2.size)
    s1 = Struct(fmt1)
    s2 = Struct(fmt2)
    add_sort_x = getattr(obj, 'add_sort' + str(op2.sort_method))
    add_new_eid_sort_x = getattr(obj, 'add_new_eid_sort' + str(op2.sort_method))
    for unused_i in range(nelements):
        edata = data[n:n+n1]
        n += n1

        out = s1.unpack(edata)
        eid_device = out[0]
        eid, dt = get_eid_dt_from_eid_device(
            eid_device, op2.nonlinear_factor, op2.sort_method)
        if op2.is_debug_file:
            op2.binary_debug.write('CBEAM-2 - eid=%i out=%s\n' % (eid, str(out)))

        #(grid, sd, sxc, sxd, sxe, sxf, smax, smin, mst, msc) = out
        add_new_eid_sort_x(dt, eid, *out[1:])

        for unused_inode in range(nnodes):
            edata = data[n:n+n2]
            n += n2
            out = s2.unpack(edata)
            # (grid, sd, sxc, sxd, sxe, sxf, smax, smin, mst, msc) = out
            add_sort_x(dt, eid, *out)
    return n


def oes_cbeam_complex_111(op2: OP2, data: bytes,
                          obj: ComplexBeamStressArray | ComplexBeamStrainArray,
                          nelements: int, nnodes: int,
                          is_magnitude_phase: bool) -> int:
    n = 0
    #itotal = obj.itotal
    n1 = 44 * op2.factor
    n2 = 40 * op2.factor

    s1 = Struct(op2._endian + mapfmt(op2._analysis_code_fmt + b'i9f', op2.size))
    s2 = Struct(mapfmt(op2._endian + b'i9f', op2.size))

    add_sort_x = getattr(obj, 'add_sort' + str(op2.sort_method))
    for unused_i in range(nelements):
        edata = data[n:n+n1]
        n += n1
        out1 = s1.unpack(edata)
        eid_device = out1[0]
        eid, dt = get_eid_dt_from_eid_device(
            eid_device, op2.nonlinear_factor, op2.sort_method)
        if op2.is_debug_file:
            op2.binary_debug.write('CBEAM-2 - eid=%i out1=%s\n' % (eid, str(out1)))

        (grid, sd,
         excr, exdr, exer, exfr,
         exci, exdi, exei, exfi) = out1[1:]

        if is_magnitude_phase:
            exc = polar_to_real_imag(excr, exci)
            exd = polar_to_real_imag(exdr, exdi)
            exe = polar_to_real_imag(exer, exei)
            exf = polar_to_real_imag(exfr, exfi)
        else:
            exc = complex(excr, exci)
            exd = complex(exdr, exdi)
            exe = complex(exer, exei)
            exf = complex(exfr, exfi)

        add_sort_x(dt, eid, grid, sd,
                   exc, exd, exe, exf)

        for unused_inode in range(nnodes):
            edata = data[n:n+n2]
            n += n2
            out2 = s2.unpack(edata)
            (grid, sd,
             excr, exdr, exer, exfr,
             exci, exdi, exei, exfi) = out2

            if is_magnitude_phase:
                exc = polar_to_real_imag(excr, exci)
                exd = polar_to_real_imag(exdr, exdi)
                exe = polar_to_real_imag(exer, exei)
                exf = polar_to_real_imag(exfr, exfi)
            else:
                exc = complex(excr, exci)
                exd = complex(exdr, exdi)
                exe = complex(exer, exei)
                exf = complex(exfr, exfi)

            add_sort_x(dt, eid, grid, sd,
                       exc, exd, exe, exf)
            if op2.is_debug_file:
                op2.binary_debug.write('CBEAM-2 - eid=%i out2=%s\n' % (eid, str(out2)))
    return n


def oes_cbeam_random_67(op2: OP2, data: bytes,
                        obj: RandomBeamStressArray | RandomBeamStrainArray,
                        nelements: int, nnodes: int, dt: Any) -> int:
    n = 0
    n1 = 28 * op2.factor
    n2 = 24 * op2.factor # 6*4
    s1 = Struct(op2._endian + op2._analysis_code_fmt + b'i5f')
    s2 = Struct(op2._endian + b'i5f')
    add_sort_x = getattr(obj, 'add_sort' + str(op2.sort_method))
    add_eid_sort_x = getattr(obj, 'add_eid_sort' + str(op2.sort_method))
    for unused_i in range(nelements):
        edata = data[n:n+n1]
        n += n1

        out = s1.unpack(edata)
        eid_device = out[0]

        eid, dt = get_eid_dt_from_eid_device(
            eid_device, op2.nonlinear_factor, op2.sort_method)
        if op2.is_debug_file:
            op2.binary_debug.write('CBEAM-2 - eid=%i out=%s\n' % (eid, str(out)))

        #(grid, sd, sxc, sxd, sxe, sxf, smax, smin, mst, msc) = out
        #(grid, sd, sxc, sxd, sxe, sxf) = out

        # grid, sd, sxc, sxd, sxe, sxf = out[1:]
        add_eid_sort_x(dt, eid, *out[1:])

        for unused_inode in range(nnodes):
            edata = data[n:n+n2]
            n += n2
            out = s2.unpack(edata)
            # (grid, sd, sxc, sxd, sxe, sxf, smax, smin, mst, msc) = out
            # grid, sd, sxc, sxd, sxe, sxf = out
            add_sort_x(dt, eid, *out)
    return n
