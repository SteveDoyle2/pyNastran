from __future__ import annotations
from struct import Struct
from typing import Union, Any, TYPE_CHECKING
import numpy as np

from pyNastran.op2.op2_interface.op2_reader import mapfmt
from pyNastran.op2.op2_helper import polar_to_real_imag

from pyNastran.op2.tables.utils import get_eid_dt_from_eid_device
from pyNastran.op2.tables.oes_stressStrain.real.oes_bars import RealBarStressArray, RealBarStrainArray
#from pyNastran.op2.tables.oes_stressStrain.real.oes_bars100 import RealBar10NodesStressArray, RealBar10NodesStrainArray

from pyNastran.op2.tables.oes_stressStrain.real.oes_beams import (RealBeamStressArray, RealBeamStrainArray,
                                                                  #RealNonlinearBeamStressArray
                                                                  )
from pyNastran.op2.tables.oes_stressStrain.real.oes_bush import RealBushStressArray, RealBushStrainArray
#from pyNastran.op2.tables.oes_stressStrain.real.oes_bush1d import RealBush1DStressArray
from pyNastran.op2.tables.oes_stressStrain.complex.oes_fast import ComplexFastStressArray, ComplexFastStrainArray
from pyNastran.op2.tables.oes_stressStrain.complex.oes_weld import ComplexWeldStressArray, ComplexWeldStrainArray
from pyNastran.op2.tables.oes_stressStrain.real.oes_composite_plates import RealCompositePlateStressArray, RealCompositePlateStrainArray
#RealCompositePlateStressStrengthRatioArray, RealCompositePlateStrainStrengthRatioArray = None, None
#RealCompositePlateStrainStrengthRatioArray = None
from pyNastran.op2.tables.oes_stressStrain.real.oes_composite_plates_strength_ratio import RealCompositePlateStressStrengthRatioArray  # , RealCompositePlateStrainStrengthRatioArray
#from pyNastran.op2.tables.oes_stressStrain.real.oes_gap import NonlinearGapStressArray
from pyNastran.op2.tables.oes_stressStrain.real.oes_plates import RealPlateStressArray, RealPlateStrainArray
#from pyNastran.op2.tables.oes_stressStrain.real.oes_plate_strain import RealCPLSTRNPlateStressArray, RealCPLSTRNPlateStrainArray
from pyNastran.op2.tables.oes_stressStrain.real.oes_rods import RealRodStressArray, RealRodStrainArray
from pyNastran.op2.tables.oes_stressStrain.real.oes_shear import RealShearStrainArray, RealShearStressArray
from pyNastran.op2.tables.oes_stressStrain.real.oes_solids import RealSolidStrainArray, RealSolidStressArray
from pyNastran.op2.tables.oes_stressStrain.real.oes_solids_nx import RealSolidStressArrayNx, RealSolidStrainArrayNx
from pyNastran.op2.tables.oes_stressStrain.real.oes_springs import (RealSpringStressArray, RealSpringStrainArray,
                                                                    #RealNonlinearSpringStressArray
                                                                    )
from pyNastran.op2.tables.oes_stressStrain.real.oes_triax import RealTriaxStressArray, RealTriaxStrainArray
from pyNastran.op2.tables.oes_stressStrain.real.oes_bend import RealBendStressArray, RealBendStrainArray


from pyNastran.op2.tables.oes_stressStrain.complex.oes_bars import ComplexBarStressArray, ComplexBarStrainArray
from pyNastran.op2.tables.oes_stressStrain.complex.oes_beams import ComplexBeamStressArray, ComplexBeamStrainArray
from pyNastran.op2.tables.oes_stressStrain.complex.oes_bush import (ComplexCBushStressArray, ComplexCBushStrainArray)
#from pyNastran.op2.tables.oes_stressStrain.complex.oes_bush1d import ComplexCBush1DStressArray
from pyNastran.op2.tables.oes_stressStrain.complex.oes_plates import (
    ComplexPlateStressArray, ComplexPlateStrainArray)
from pyNastran.op2.tables.oes_stressStrain.complex.oes_composite_plates import (
    ComplexLayeredCompositeStrainArray, ComplexLayeredCompositeStressArray,
    ComplexLayeredCompositesArray)

from pyNastran.op2.tables.oes_stressStrain.complex.oes_plates_vm import (
    ComplexPlateVMStressArray, ComplexPlateVMStrainArray)
from pyNastran.op2.tables.oes_stressStrain.complex.oes_triax import ComplexTriaxStressArray
from pyNastran.op2.tables.oes_stressStrain.complex.oes_rods import ComplexRodStressArray, ComplexRodStrainArray
from pyNastran.op2.tables.oes_stressStrain.complex.oes_shear import ComplexShearStressArray, ComplexShearStrainArray
from pyNastran.op2.tables.oes_stressStrain.complex.oes_solids import ComplexSolidStressArray, ComplexSolidStrainArray
from pyNastran.op2.tables.oes_stressStrain.complex.oes_springs import ComplexSpringStressArray, ComplexSpringStrainArray
#from pyNastran.op2.tables.oes_stressStrain.complex.oes_bend import ComplexBendStressArray, ComplexBendStrainArray

from pyNastran.op2.tables.oes_stressStrain.random.oes_rods import RandomRodStressArray, RandomRodStrainArray
from pyNastran.op2.tables.oes_stressStrain.random.oes_bars import RandomBarStressArray, RandomBarStrainArray
from pyNastran.op2.tables.oes_stressStrain.random.oes_beams import RandomBeamStressArray, RandomBeamStrainArray
#from pyNastran.op2.tables.oes_stressStrain.random.oes_bend import RandomBendStressArray, RandomBendStrainArray
from pyNastran.op2.tables.oes_stressStrain.random.oes_plates import RandomPlateStressArray, RandomPlateStrainArray
from pyNastran.op2.tables.oes_stressStrain.random.oes_plates_vm import RandomPlateVMStressArray, RandomPlateVMStrainArray
from pyNastran.op2.tables.oes_stressStrain.random.oes_solids import RandomSolidStressArray, RandomSolidStrainArray
from pyNastran.op2.tables.oes_stressStrain.random.oes_shear import RandomShearStressArray, RandomShearStrainArray
#from pyNastran.op2.tables.oes_stressStrain.random.oes_composite_plates import RandomCompositePlateStressArray, RandomCompositePlateStrainArray

#from pyNastran.op2.tables.oes_stressStrain.oes_nonlinear_rod import RealNonlinearRodArray
#from pyNastran.op2.tables.oes_stressStrain.oes_nonlinear_bush import RealNonlinearBushArray
#from pyNastran.op2.tables.oes_stressStrain.oes_hyperelastic import (
    #HyperelasticQuadArray)
#from pyNastran.op2.tables.oes_stressStrain.oes_nonlinear import RealNonlinearPlateArray, RealNonlinearSolidArray
if TYPE_CHECKING:  # pragma: no cover
    from pyNastran.op2.op2 import OP2


def oes_celas_real_2(op2: OP2, data: bytes,
                     obj: Union[RealSpringStressArray, RealSpringStrainArray],
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


def oes_celas_complex_3(op2: OP2, data: bytes, obj: Union[ComplexSpringStressArray, ComplexSpringStrainArray],
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


def oes_cquad4_33_complex_vm_17(op2: OP2, data: bytes,
                                obj: Union[ComplexPlateVMStressArray, ComplexPlateVMStrainArray],
                                nelements: int, ntotal: int,
                                is_magnitude_phase: bool) -> int:
    """
    OESVM1/2 - Table of element stresses or strains with von Mises
    OSTRVM1/2  for frequency response results.

    """
    #self.log.info(f'nelements = {nelements}')
    n = 0
    struct1 = Struct(mapfmt(op2._endian + op2._analysis_code_fmt + b'16f', op2.size))
    cen = 0 # CEN/4
    add_sort_x = getattr(obj, 'add_sort' + str(op2.sort_method))

    for unused_i in range(nelements):
        edata = data[n:n+ntotal]
        out = struct1.unpack(edata)

        (eid_device,
         fd1, sx1r, sx1i, sy1r, sy1i, txy1r, txy1i, von_mises1,
         fd2, sx2r, sx2i, sy2r, sy2i, txy2r, txy2i, von_mises2) = out

        eid, dt = get_eid_dt_from_eid_device(
            eid_device, op2.nonlinear_factor, op2.sort_method)
        if op2.is_debug_file:
            op2.binary_debug.write('  eid=%i C=[%s]\n' % (
                eid, ', '.join(['%r' % di for di in out])))

        if is_magnitude_phase:
            sx1 = polar_to_real_imag(sx1r, sx1i)
            sx2 = polar_to_real_imag(sx2r, sx2i)
            sy1 = polar_to_real_imag(sy1r, sy1i)
            sy2 = polar_to_real_imag(sy2r, sy2i)
            txy1 = polar_to_real_imag(txy1r, txy1i)
            txy2 = polar_to_real_imag(txy2r, txy2i)
        else:
            sx1 = complex(sx1r, sx1i)
            sx2 = complex(sx2r, sx2i)
            sy1 = complex(sy1r, sy1i)
            sy2 = complex(sy2r, sy2i)
            txy1 = complex(txy1r, txy1i)
            txy2 = complex(txy2r, txy2i)
        #print(dt, eid, cen, sx1, sy1, txy1, max_shear1)
        #print(dt, eid, cen, sx2, sy2, txy2, max_shear2)
        add_sort_x(dt, eid, cen,
                   fd1, sx1, sy1, txy1, von_mises1,
                   fd2, sx2, sy2, txy2, von_mises2)
        n += ntotal
    return n


def oes_cbeam_real_111(op2: OP2, data: bytes,
                       obj: Union[RealBeamStressArray, RealBeamStrainArray],
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


def oes_crod_real_5(op2: OP2, data: bytes, obj: Union[RealRodStressArray, RealRodStrainArray],
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


def oes_crod_complex_5(op2: OP2, data: bytes, obj: Union[ComplexRodStressArray, ComplexRodStrainArray],
                       nelements: int, ntotal: int, dt, is_magnitude_phase: bool) -> int:
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
             obj: Union[RandomRodStressArray, RandomRodStrainArray],
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


def oes_cbeam_complex_111(op2: OP2, data: bytes,
                          obj: Union[ComplexBeamStressArray, ComplexBeamStrainArray],
                          nelements: int, nnodes: int,
                          is_magnitude_phase: bool) -> int:
    n = 0
    #itotal = obj.itotal
    n1 = 44 * op2.factor
    n2 = 40 * op2.factor

    s1 = Struct(op2._endian + mapfmt(op2._analysis_code_fmt + b'i9f', op2.size))
    if op2.size == 4:
        s2 = Struct(op2._endian + b'i9f')
    else:
        s2 = Struct(op2._endian + b'q9d')

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
                        obj: Union[RandomBeamStressArray, RandomBeamStrainArray],
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
        add_eid_sort_x(dt, eid, *out[1:])

        for unused_inode in range(nnodes):
            edata = data[n:n+n2]
            n += n2
            out = s2.unpack(edata)
            # (grid, sd, sxc, sxd, sxe, sxf, smax, smin, mst, msc) = out
            add_sort_x(dt, eid, *out)
    return n


def oes_cquad4_33_complex_15(op2: OP2,
                             data: bytes,
                             obj: Union[ComplexPlateStressArray, ComplexPlateStrainArray],
                             nelements: int, ntotal: int,
                             is_magnitude_phase: bool) -> int:
    """
    NX/MSC 2018.2
    2 FD1    RS Z1 = Fibre distance
    3 EX1R   RS Normal in x at Z1
    4 EX1I   RS Normal in x at Z1
    5 EY1R   RS Normal in y at Z1
    6 EY1I   RS Normal in y at Z1
    7 EXY1R  RS Shear in xy at Z1
    8 EXY1I  RS Shear in xy at Z1
    9 FD2    RS Z2 = Fibre distance
    10 EX2R  RS Normal in x at Z2
    11 EX2I  RS Normal in x at Z2
    12 EY2R  RS Normal in y at Z2
    13 EY2I  RS Normal in y at Z2
    14 EXY2R RS Shear in xy at Z2
    15 EXY2I RS Shear in xy at Z2
    """
    n = 0
    s1 = Struct(op2._endian + op2._analysis_code_fmt + b'14f')
    cen = 0  # 'CEN/4'
    add_sort_x = getattr(obj, 'add_sort' + str(op2.sort_method))

    for unused_i in range(nelements):
        edata = data[n:n+ntotal]  # 4*15=60
        n += ntotal
        out = s1.unpack(edata)  # 15
        (eid_device,
         fd1, sx1r, sx1i, sy1r, sy1i, txy1r, txy1i,
         fd2, sx2r, sx2i, sy2r, sy2i, txy2r, txy2i) = out

        eid, dt = get_eid_dt_from_eid_device(
            eid_device, op2.nonlinear_factor, op2.sort_method)
        if op2.is_debug_file:
            op2.binary_debug.write('  eid=%i C=%s\n' % (eid, str(out)))

        if is_magnitude_phase:
            sx1 = polar_to_real_imag(sx1r, sx1i)
            sx2 = polar_to_real_imag(sx2r, sx2i)
            sy1 = polar_to_real_imag(sy1r, sy1i)
            sy2 = polar_to_real_imag(sy2r, sy2i)
            txy1 = polar_to_real_imag(txy1r, txy1i)
            txy2 = polar_to_real_imag(txy2r, txy2i)
        else:
            sx1 = complex(sx1r, sx1i)
            sx2 = complex(sx2r, sx2i)
            sy1 = complex(sy1r, sy1i)
            sy2 = complex(sy2r, sy2i)
            txy1 = complex(txy1r, txy1i)
            txy2 = complex(txy2r, txy2i)

        add_sort_x(dt, eid, cen,
                   fd1, sx1, sy1, txy1,
                   fd2, sx2, sy2, txy2)
    return n


def oes_cquad4_random_vm_57(op2: OP2,
                            data: bytes,
                            obj: Union[ComplexPlateVMStressArray, ComplexPlateVMStrainArray],
                            nelements: int, ntotal: int, nnodes: int,
                            dt: Any) -> int:
    """
    name    ntotal numwide
    CQUAD8  228    57
    CTRIA6  184
    """

    #ntotal = 228 # 57*4
    n = 0
    cen = 0
    if op2.size == 4:
        struct1 = Struct(op2._endian + b'i4si 10f')
        struct2 = Struct(op2._endian + b'i 10f')
    else:
        print('ntotal =', ntotal, 228*2)
        struct1 = Struct(op2._endian + b'i4si 10f')
        struct2 = Struct(op2._endian + b'i 10f')
        raise NotImplementedError('64-bit')

    nnodes_corners = nnodes - 1
    #print('random-57: nelements =', nelements)
    factor = op2.factor

    ntotal1 = 52 * factor
    ntotal2 = 44 * factor
    add_sort_x = getattr(obj, 'add_sort' + str(op2.sort_method))
    for unused_i in range(nelements):
        ni = 0
        edata = data[n:n+ntotal]
        #self.show_data(edata)
        #out = struct1.unpack(edata)  # len=17*4
        (eid_device, grid1, four,
         fd1, ox1, oy1, txy1, vm1,
         fd2, ox2, oy2, txy2, vm2) = struct1.unpack(edata[ni:ni+ntotal1])  # 13*4=52
        eid, dt = get_eid_dt_from_eid_device(
            eid_device, op2.nonlinear_factor, op2.sort_method)
        #print(eid, cen)
        add_sort_x(dt, eid, cen,
                   fd1, ox2, oy2, txy2, vm1,
                   fd2, ox2, oy2, txy2, vm2)
        ni += ntotal1
        #print(eid, grid1,
              #fd1, ox1, oy1, txy1, vm1,
              #fd2, ox2, oy2, txy2, vm2)

        for unused_jnode in range(nnodes_corners): # 4*44=176
            (grid,
             fd1, ox1, oy1, txy1, vm1,
             fd2, ox2, oy2, txy2, vm2) = struct2.unpack(edata[ni:ni+ntotal2])  # 11*4=44
            #print(eid, grid)
            ni += ntotal2
            #print(grid,
                  #fd1, ox1, oy1, txy1, vm1,
                  #fd2, ox2, oy2, txy2, vm2)
            add_sort_x(dt, eid, grid,
                       fd1, ox2, oy2, txy2, vm1,
                       fd2, ox2, oy2, txy2, vm2)
        assert ni == ntotal, f'ni={ni}'
        n += ntotal
        #print()
    #assert n == ndata
    #n = self._not_implemented_or_skip(data, ndata, msg)
    #nelements = None
    #ntotal = None
    return n


def oes_cquad4_144_complex_77(op2: OP2,
                              data: bytes,
                              obj: Union[ComplexPlateStressArray, ComplexPlateStrainArray],
                              nelements: int, nnodes: int,
                              dt: Any,
                              is_magnitude_phase: bool) -> int:
    #print('nelements =', nelements)
    n = 0
    grid_center = 0
    if op2.size == 4:
        s1 = op2.struct_2i  # 2
    else:
        s1 = op2.struct_2q  # 2
    s2 = Struct(op2._endian + mapfmt(b'i14f', op2.size)) # 15

    ntotal1 = 8 * op2.factor
    ntotal2 = 60 * op2.factor  # 4*15

    add_sort_x = getattr(obj, 'add_sort' + str(op2.sort_method))
    for unused_i in range(nelements):
        (eid_device, _) = s1.unpack(data[n:n+ntotal1])
        n += ntotal1

        edata = data[n:n+ntotal2]
        n += ntotal2
        out = s2.unpack(edata)  # len=15*4

        eid = op2._check_id(eid_device, 'element', out)
        if op2.is_debug_file:
            op2.binary_debug.write('%s\n' % (str(out)))
        (grid,
         fd1, sx1r, sx1i, sy1r, sy1i, txy1r, txy1i,
         fd2, sx2r, sx2i, sy2r, sy2i, txy2r, txy2i) = out
        #grid_center = 'CEN/%i' % grid   # this is correct, but fails

        if is_magnitude_phase:
            sx1 = polar_to_real_imag(sx1r, sx1i)
            sy1 = polar_to_real_imag(sy1r, sy1i)
            sx2 = polar_to_real_imag(sx2r, sx2i)
            sy2 = polar_to_real_imag(sy2r, sy2i)
            txy1 = polar_to_real_imag(txy1r, txy1i)
            txy2 = polar_to_real_imag(txy2r, txy2i)
        else:
            sx1 = complex(sx1r, sx1i)
            sy1 = complex(sy1r, sy1i)
            sx2 = complex(sx2r, sx2i)
            sy2 = complex(sy2r, sy2i)
            txy1 = complex(txy1r, txy1i)
            txy2 = complex(txy2r, txy2i)

        add_sort_x(dt, eid, grid_center,
                   fd1, sx1, sy1, txy1,
                   fd2, sx2, sy2, txy2)

        for unused_node_id in range(nnodes):  # nodes pts
            edata = data[n:n + ntotal2]
            n += ntotal2
            out = s2.unpack(edata)
            if op2.is_debug_file:
                op2.binary_debug.write('%s\n' % (str(out)))
            (grid,
             fd1, sx1r, sx1i, sy1r, sy1i, txy1r, txy1i,
             fd2, sx2r, sx2i, sy2r, sy2i, txy2r, txy2i) = out

            if is_magnitude_phase:
                sx1 = polar_to_real_imag(sx1r, sx1i)
                sx2 = polar_to_real_imag(sx2r, sx2i)
                sy1 = polar_to_real_imag(sy1r, sy1i)
                sy2 = polar_to_real_imag(sy2r, sy2i)
                txy1 = polar_to_real_imag(txy1r, txy1i)
                txy2 = polar_to_real_imag(txy2r, txy2i)
            else:
                sx1 = complex(sx1r, sx1i)
                sx2 = complex(sx2r, sx2i)
                sy1 = complex(sy1r, sy1i)
                sy2 = complex(sy2r, sy2i)
                txy1 = complex(txy1r, txy1i)
                txy2 = complex(txy2r, txy2i)

            add_sort_x(dt, eid, grid,
                       fd1, sx1, sy1, txy1,
                       fd2, sx2, sy2, txy2)
    return n


def oes_cquad4_complex_vm_87(op2: OP2, data: bytes,
                             obj: Union[ComplexPlateVMStressArray, ComplexPlateVMStrainArray],
                             nelements: int, nnodes: int,
                             is_magnitude_phase: bool) -> int:
    """
    etype name      numwide nodes
    33  CQUAD4-33 -> 17      (not in this function)
    82  CQUADR    -> 87      (1+4 nodes)
    144 CQUAD144  -> 87      (1+4 nodes)
    ??? CTRIA6    -> 70?     (1+3 nodes)

    2 TERM    CHAR4
    3 GRID    I
    4 FD1     RS Fiber distance at z1
    5 EX1R    RS Normal in x at z1
    6 EX1I    RS Normal in x at z1
    7 EY1R    RS Normal in y at z1
    8 EY1I    RS Normal in y at z1
    9 ETXY1R  RS Shear in xy at z1
    10 ETXY1I RS Shear in xy at z1
    11 VM1    RS von Mises at z1
    12 FD2    RS Fiber distance at z2
    13 EX2R   RS Normal in x at z2
    14 EX2I   RS Normal in x at z2
    15 EY2R   RS Normal in y at z2
    16 EY2I   RS Normal in y at z2
    17 ETXY2R RS Shear in xy at z2
    18 ETXY2I RS Shear in xy at z2
    19 VM2    RS von Mises at z2

    """
    #strings = (b'20.0 CEN/ 4 \xcd\xccL\xbd\x91@\x1aD\xebX\xac\xc1b\x80FB\x96\xc9\xdd\xbf\xe7\xbe\n,\xaah\xc5\xa9b\x87\x14D\xcd\xccL=\x91@\x1a\xc4\xebX\xacAb\x80F\xc2\x96\xc9\xdd?\xe7\xbe\n\xac\xaah\xc5)b\x87\x14D\x01\x00\x00\x00\xcd\xccL\xbd\xc7\xd5\x9fD\xc9\x952\xc2b\x80FC\x96\xc9\xdd\xc0\xe7\xbe\n,\xaah\xc5\xa9\xca\x13\x95D\xcd\xccL=\xc7\xd5\x9f\xc4\xc9\x952Bb\x80F\xc3\x96\xc9\xdd@\xe7\xbe\n\xac\xaah\xc5)\xca\x13\x95D\x02\x00\x00\x00\xcd\xccL\xbd\xbf\xa62\xc2\xd4\x9b\xc7?b\x80\xc6\xc2\x96\xc9]@\xe7\xbe\n,\xaah\xc5\xa9nL\xacB\xcd\xccL=\xbf\xa62B\xd4\x9b\xc7\xbfb\x80\xc6B\x96\xc9]\xc0\xe7\xbe\n\xac\xaah\xc5)nL\xacB\x03\x00\x00\x00\xcd\xccL\xbd\xbf\xa62\xc2\xd4\x9b\xc7?b\x80\xc6\xc2\x96\xc9]@\xe7\xbe\n,\xaah\xc5\xa9nL\xacB\xcd\xccL=\xbf\xa62B\xd4\x9b\xc7\xbfb\x80\xc6B\x96\xc9]\xc0\xe7\xbe\n\xac\xaah\xc5)nL\xacB\x04\x00\x00\x00\xcd\xccL\xbd\xc7\xd5\x9fD\xc9\x952\xc2b\x80FC\x96\xc9\xdd\xc0\xe7\xbe\n,\xaah\xc5\xa9\xca\x13\x95D\xcd\xccL=\xc7\xd5\x9f\xc4\xc9\x952Bb\x80F\xc3\x96\xc9\xdd@\xe7\xbe\n\xac\xaah\xc5)\xca\x13\x95D',)
    # ints    = (20.0, CEN/, 4,
    #           -0.05, 617.0, -21.5, 49.6, -1.7, 738901735,              -1446680406,            1142196066,
    #           0.05, -617.0,  21.5, -49.6, 1.7, -1408581913, 700803242, 1142196066, 1, -1119040307, 1151325639, -1036872247, 1128693858, -1059206762, 738901735, -1446680406, 1150620618, 1028443341, -996158009, 1110611401, -1018789790, 1088276886, -1408581913, 700803242, 1150620618, 2, -1119040307, -1036867905, 1070046164, -1027178398, 1079888278, 738901735, -1446680406, 1118588014, 1028443341, 1110615743, -1077437484, 1120305250, -1067595370, -1408581913, 700803242, 1118588014, 3, -1119040307, -1036867905, 1070046164, -1027178398, 1079888278, 738901735, -1446680406, 1118588014, 1028443341, 1110615743, -1077437484, 1120305250, -1067595370, -1408581913, 700803242, 1118588014, 4, -1119040307, 1151325639, -1036872247, 1128693858, -1059206762, 738901735, -1446680406, 1150620618, 1028443341, -996158009, 1110611401, -1018789790, 1088276886, -1408581913, 700803242, 1150620618)
    #floats  = (20.0, CEN/, 4,
    #           -0.05,  617.0,   -21.5,   49.6,  -1.7,  1.9716951595721843e-12, -8.766713754677219e-14,  594.1,
    #            0.05, -617.0,    21.5,  -49.6,   1.7, -1.9716951595721843e-12,  8.766713754677219e-14,  594.1,
    #        1.4012984e-45,
    #           -0.05,  1278.68, -44.6,  198.5,  -6.9,  1.9716951595721843e-12, -8.766713754677219e-14, 1192.6,
    #            0.05, -1278.68,  44.6, -198.5,   6.9, -1.9716951595721843e-12,  8.766713754677219e-14, 1192.6,
    #        2.802596928649634e-45,
    #           -0.05, -44.6,      1.5,  -99.2,   3.4,  1.9716951595721843e-12, -8.766713754677219e-14, 86.149,
    #            0.05,  44.6,     -1.5,   99.2,  -3.4, -1.9716951595721843e-12,  8.766713754677219e-14, 86.149,
    #        4.203895392974451e-45,
    #           -0.05, -44.6,      1.5,  -99.2,   3.4,  1.9716951595721843e-12, -8.766713754677219e-14, 86.149,
    #            0.05,  44.6,     -1.5,   99.2,  -3.4, -1.9716951595721843e-12,  8.766713754677219e-14, 86.149,
    #        5.605193857299268e-45,
    #           -0.05, 1278.6,    -44.6, 198.5,  -6.9,  1.9716951595721843e-12, -8.766713754677219e-14, 1192.6,
    #           0.05, -1278.6,     44.6,-198.5,   6.9, -1.9716951595721843e-12,  8.766713754677219e-14, 1192.6)
    # -------------------------
    #ints    = (0,    'CEN/',
    #           5,     -0.025, 1131762990, -1056779073, 1137493831, -1051039598, 1133740074, -1055043205, 1142680183,
    #                   0.025, -1016740689, 1089802185, -1009242998, 1097240578, -1013972178, 1092196624, 1142689762,
    #           80000, -0.025, -1039337422, 1066818106, -1022560206, 1083595322, 1118381010, -1070307877, 1128400256,
    #                   0.025, 1111837962, -1076728084, 1128615178, -1059950868, -1029871408, 1076355832, 1130378350,
    #           80002, -0.025, 1144778127, -1043828805, 1161555343, -1027051589, 1134339017, -1054404396, 1160567412,
    #                   0.025, -1002241635, 1104149606, -985464419, 1120926822, -1014451576, 1091685317, 1160927605,
    #           80202, -0.025, -1013123432, 1093101862, -996473336, 1109743497, 1143660503, -1045020817, 1154190671,
    #                   0.025, 1133381240, -1055425923, 1150079339, -1038733090, -1002862704, 1103487198, 1154274046,
    #           80200, -0.025, 1141944942, -1046850567, 1095277553, -1093274063, 1120227570, -1068338410, 1142258188,
    #                   0.025, -1007929795, 1098641187, 1116060214, -1072783145, -1033402643, 1073147904, 1140969037,
    #           2.5,  'CEN/',
    #           5,     -0.025, 1131832225, -1056704978, 1137555116, -1050908425, 1133787322, -1054942075, 1142729036,
    #                   0.025, -1016671483, 1089950313, -1009181671, 1097371840, -1013924936, 1092297741, 1142738537,
    #           80000, -0.025, -1039289241, 1066921230, -1022512025, 1083698446, 1118448079, -1070164312, 1128472712,
    #                   0.025, 1111886248, -1076624736, 1128663464, -1059847520, -1029804359, 1076499355, 1130448349,
    #           80002, -0.025, 1144846119, -1043683269, 1161623335, -1026906053, 1134386659, -1054302422, 1160629636,
    #                   0.025, -1002173624, 1104295182, -985396408, 1121072398, -1014403979, 1091787195, 1160989873,
    #           80202, -0.025, -1013031023, 1093299683, -996403997, 1109891922, 1143714425, -1044905405, 1154266576,
    #                   0.025, 1133473626, -1055228152, 1150148641, -1038584742, -1002808755, 1103602670, 1154348580,
    #           80200, -0.025, 1141995401, -1046742558, 1095362185, -1093092785, 1120294524, -1068195089, 1142310151,
    #                   0.025, -1007828983, 1098856980, 1116049937, -1072805159, -1033268970, 1073434045, 1141020308,
    #           10.0, 'CEN/',
    #           5,     -0.025, 1132706983, -1055459734, 1138546383, -1048703317, 1134552107, -1053240134, 1143519784,
    #                   0.025, -1015552254, 1091479459, -1008189778, 1099242966, -1013160238, 1093999495, 1143528113,
    #           80000, -0.025, -1038510246, 1068653763, -1021733030, 1085430979, 1119535604, -1067741982, 1129646715,
    #                   0.025, 1112666820, -1074888834, 1129444036, -1058111618, -1028717126, 1078921059, 1131585244,
    #           80002, -0.025,  1145947464, -1041231474, 1162724680, -1024454258, 1135157762, -1052586493, 1161637710,
    #                   0.025, -1001071999, 1106747574, -984294783, 1123524790, -1013633545, 1093501695, 1161998610,
    #           80202, -0.025, -1011529318,  1096648297, -995279005, 1112398467, 1144586718, -1042964811, 1155498288,
    #                   0.025, 1134974980, -1051880286, 1151273087, -1036079365, -1001936046, 1105544152, 1155559868,
    #           80200, -0.025,  1142813256, -1044921274, 1096762940, -1090229644, 1121380799, -1065774838, 1143152359,
    #                   0.025, -1006413905, 1100701913, 1115879378, -1073189541, -1031449359, 1076004431, 1141852280)
    #floats  = (0.0,  'CEN/',
    #           5,       -0.025,  245.33273315429688, -8.176939964294434,   409.5568542480469, -13.650529861450195,  295.00128173828125,  -9.832392692565918,
    #           623.8668, 0.025, -229.76829528808594,  7.658176898956299, -432.34796142578125,  14.410158157348633, -288.02484130859375,   9.599868774414062, 624.2481689453125,
    #           1.12-40, -0.025, -35.24237823486328,   1.174628496170044, -140.96951293945312,   4.698513984680176,   84.56996154785156, -2.8187167644500732,
    #           194.4375, 0.025,  49.325233459472656, -1.644010066986084,  197.30093383789062,  -6.576040267944336,   -78.7047119140625,   2.623228073120117, 224.20480346679688,
    #           1.12106, -0.025,  751.7118530273438, -25.054555892944336,   3006.847412109375, -100.21822357177734,   313.2795715332031, -10.441608428955078,
    #           2765.65,  0.025, -780.0252075195312, 25.998241424560547, -3120.100830078125, 103.99296569824219, -273.394775390625, 9.112248420715332, 2853.591064453125,
    #           1.1238,  -0.025, -313.926513671875, 10.463171005249023, -1240.1884765625, 41.33548355102539, 683.4974975585938, -22.780973434448242,
    #           1628.4,   0.025, 284.050537109375, -9.46740436553955, 1126.5443115234375, -37.54772186279297, -742.1181640625, 24.734798431396484, 1638.593505859375,
    #           1.12,    -0.025, 578.7879638671875, -19.29100227355957, 12.538071632385254, -0.4178939163684845, 98.65809631347656, -3.2882742881774902,
    #           597.90,   0.025, -472.4237365722656, 15.745882987976074, 66.86369323730469, -2.228566884994507, -57.88176345825195, 1.92919921875, 519.2234497070312,
    #           2.5,  'CEN/',
    #           5,          -0.025, 246.38917541503906, -8.247602462768555, 411.4271240234375, -13.775626182556152, 296.44317626953125, -9.928837776184082,
    #           626.671875,  0.025, -230.8242950439453, 7.728809833526611, -434.2195129394531, 14.53533935546875, -289.466552734375, 9.696301460266113, 627.2251586914062,
    #           1.1217e-40, -0.025, -35.42617416381836, 1.1869218349456787, -141.70469665527344, 4.747687339782715, 85.08165740966797, -2.852945327758789,
    #           195.153125,  0.025, 49.509429931640625, -1.6563301086425781, 198.0377197265625, -6.6253204345703125, -79.21625518798828, 2.6574466228485107, 225.2729034423828,
    #           1.1212e-40, -0.025, 755.8617553710938, -25.332143783569336, 3023.447021484375, -101.32857513427734, 314.7334899902344, -10.538858413696289,
    #           2780.65625,  0.025, -784.17626953125, 26.27590560913086, -3136.705078125, 105.10362243652344, -274.8473205566406, 9.209406852722168, 2868.793212890625,
    #           1.1358e-40, -0.025, -316.7466125488281, 10.651827812194824, -1248.6527099609375, 41.90167999267578, 686.7886352539062, -23.0011043548584,
    #           1637.40625,  0.025, 286.86993408203125, -9.656013488769531, 1135.0040283203125, -38.113624572753906, -745.4109497070312, 24.95504379272461, 1647.69189453125,
    #           1.1383e-40, -0.025, 581.8677368164062, -19.497013092041016, 12.618782997131348, -0.42329642176628113, 99.16891479492188, -3.3224446773529053,
    #           601.471875,  0.025, -475.5002746582031, 15.951679229736328, 66.78528594970703, -2.223318338394165, -58.391685485839844, 1.96330988407135, 522.352783203125,
    #           10.0, 'CEN/',
    #           5,          -0.025, 263.4738464355469, -9.435159683227539, 441.6781921386719, -15.878581047058105, 319.7825622558594, -11.551935195922852,
    #           674.140625,  0.025, -247.90237426757812, 8.915926933288574, -464.48968505859375, 16.639568328857422, -312.80328369140625, 11.319220542907715, 675.4170532226562,
    #           1.1337e-40, -0.025, -38.397804260253906, 1.3934558629989624, -153.9121704101562, 5.57382345199585, 93.37881469726562, -3.430473804473877,
    #           213.310547,  0.025, 52.48707580566406, -1.8632657527923584, 209.94830322265625, -7.453063011169434, -87.51118469238281, 3.234825849533081, 242.62054443359375,
    #           1.1202e-40, -0.025, 823.08251953125, -30.00857162475586, 3292.330078125, -120.03428649902344, 338.26568603515625, -12.17529582977295,
    #           3026296875,  0.025, -851.4141235351562, 30.953472137451172, -3405.656494140625, 123.81388854980469, -298.3591003417969, 10.844481468200684, 3115.06689453125,
    #           1.1298e-40, -0.025, -362.57501220703125, 13.845314979553223, -1385.9808349609375, 51.4633903503418, 740.0291748046875, -26.70249366760254,
    #           1788109375,  0.025, 332.6876220703125, -12.848787307739258, 1272.2655029296875, -47.67087936401367, -798.6768798828125, 28.658126831054688, 1795.55419921875,
    #           1.1293e-40, -0.025, 631.78564453125, -22.970836639404297, 13.954647064208984, -0.517249345779419, 107.45653533935547, -3.899477481842041,
    #           652.410938,  0.025, -525.3700561523438, 19.42228889465332, 65.48402404785156, -2.1316745281219482, -66.66590118408203, 2.539447546005249, 573.13232421875)
    n = 0
    ntotal1 = 8 * op2.factor  # 2*4
    ntotal2 = 68 * op2.factor # 17*4
    add_sort_x = getattr(obj, 'add_sort' + str(op2.sort_method))
    #print('ndata', len(data))
    #print('ntotal', ntotal1+nnodes*ntotal2)
    #print('nelements =', nelements)
    for i in range(nelements):
        #print(f'{i}/{nelements}')
        edata1 = data[n:n+ntotal1]
        #self.show_data(edata1)
        n += ntotal1

        fmt1 = mapfmt(op2._endian + op2._analysis_code_fmt + b'4s', op2.size)
        #fmt2 = mapfmt(op2._endian + b'i 5f2if 5f2if', op2.size)
        fmt2 = mapfmt(op2._endian + b'i 16f', op2.size)
        s1 = Struct(fmt1)
        s2 = Struct(fmt2)

        eid_device, cen = s1.unpack(edata1)
        #self.show_data(data, types='ifs')
        assert cen == b'CEN/', (eid_device, cen)

        eid, dt = get_eid_dt_from_eid_device(
            eid_device, op2.nonlinear_factor, op2.sort_method)

        for inode in range(nnodes):
            edata2 = data[n:n+ntotal2]

            #DISTANCE             NORMAL-X                     NORMAL-Y                    SHEAR-XY             VON MISES
            out = s2.unpack(edata2)
            (nid,
             fd1, oxx1r, oxx1i, oyy1r, oyy1i, txy1r, txy1i, ovm1,
             fd2, oxx2r, oxx2i, oyy2r, oyy2i, txy2r, txy2i, ovm2) = out
            fd1, oxx1r, oxx1i, oyy1r, oyy1i, txy1r, txy1i, ovm1,
            if inode == 0:
                nid = 0

            if is_magnitude_phase:
                oxx1 = polar_to_real_imag(oxx1r, oxx1i)
                oxx2 = polar_to_real_imag(oxx2r, oxx2i)
                oyy1 = polar_to_real_imag(oyy1r, oyy1i)
                oyy2 = polar_to_real_imag(oyy2r, oyy2i)
                txy1 = polar_to_real_imag(txy1r, txy1i)
                txy2 = polar_to_real_imag(txy2r, txy2i)
            else:
                oxx1 = complex(oxx1r, oxx1i)
                oxx2 = complex(oxx2r, oxx2i)
                oyy1 = complex(oyy1r, oyy1i)
                oyy2 = complex(oyy2r, oyy2i)
                txy1 = complex(txy1r, txy1i)
                txy2 = complex(txy2r, txy2i)

            #print((eid, dt), nid, f'{fd1:g}, {fd2:g}')
            #print('real', f'{oxx1r:g}, {oxx1i:g}, {oyy1r:g}, {oyy1i:g}, {txy1r:g}, {txy1i:g}, {ovm1:g}')
            #print('imag', oxx2i, oxx2i, oyy2r, oyy2i, txy2r, txy2i, ovm2)
            add_sort_x(dt, eid, nid,
                       fd1, oxx1, oyy1, txy1, ovm1,
                       fd2, oxx2, oyy2, txy2, ovm2)
            #print(out)
            #assert nid in (1, 2, 3, 4, 17), ((eid, dt), nid, fd1, fd2)
            n += ntotal2
        #print('-----------')
    return n


def oes_cquad4_33_random_9(op2: OP2, data: bytes,
                           obj: Union[RandomPlateStressArray, RandomPlateStrainArray],
                           nelements: int, ntotal: int) -> int:
    n = 0
    add_sort_x = getattr(obj, 'add_sort' + str(op2.sort_method))
    #print('cquad33_9 - SORT1')
    struct1 = Struct(op2._endian + op2._analysis_code_fmt + b'8f')
    #cen = 0 # CEN/4
    for unused_i in range(nelements):
        edata = data[n:n+ntotal]
        out = struct1.unpack(edata)

        (eid_device,
         fd1, sx1, sy1, txy1,
         fd2, sx2, sy2, txy2,) = out
        eid, dt = get_eid_dt_from_eid_device(
            eid_device, op2.nonlinear_factor, op2.sort_method)
        #self.log.info(f'SORT1 no VM: eid={eid} dt={dt}')

        if op2.is_debug_file:
            op2.binary_debug.write('  eid=%i C=[%s]\n' % (
                eid, ', '.join(['%r' % di for di in out])))
        add_sort_x(dt, eid, 0,
                   fd1, sx1, sy1, txy1,
                   fd2, sx2, sy2, txy2)
        n += ntotal
    return n


def oes_ctria3_real_17(op2: OP2, data: bytes,
                       obj: Union[RealPlateStressArray, RealPlateStrainArray],
                       ntotal: int, nelements: int, dt: Any) -> int:
    n = 0
    cen = 0  # 'CEN/3'
    #assert op2.sort_method == 1, op2.code_information()

    add_new_eid_sort_x = getattr(obj, 'add_new_eid_sort' + str(op2.sort_method))
    #print(add_new_eid_sort_x)
    struct1 = Struct(op2._endian + mapfmt(op2._analysis_code_fmt + b'16f', op2.size))
    #print('nelements =', nelements, nelements*2, obj.ntotal)
    for unused_i in range(nelements):
        edata = data[n:n + ntotal]
        out = struct1.unpack(edata)
        (eid_device,
         fd1, sx1, sy1, txy1, angle1, major1, minor1, vm1,
         fd2, sx2, sy2, txy2, angle2, major2, minor2, vm2,) = out
        eid, dt = get_eid_dt_from_eid_device(
            eid_device, op2.nonlinear_factor, op2.sort_method)
        if op2.is_debug_file:
            op2.binary_debug.write('  OES CTRIA3-74 - eid=%i; C=[%s]\n' % (
                eid, ', '.join(['%r' % di for di in out])))
        #print('  OES CTRIA3-74 - eid=%i; C=[%s]\n' % (
            #eid, ', '.join(['%r' % di for di in out])))

        add_new_eid_sort_x(dt, eid, cen,
                           fd1, sx1, sy1, txy1, angle1, major1, minor1, vm1,
                           fd2, sx2, sy2, txy2, angle2, major2, minor2, vm2)
        n += ntotal
    return n


def oes_quad4_33_real_17(op2: OP2, data: bytes,
                         obj: Union[RealPlateStressArray, RealPlateStrainArray],
                         ntotal: int, nelements: int, dt: Any) -> int:
    n = 0
    cen = 0 # CEN/4
    #assert op2.sort_method == 1, op2.code_information()
    add_new_eid_sort_x = getattr(obj, 'add_new_eid_sort' + str(op2.sort_method))
    fmt = mapfmt(op2._endian + op2._analysis_code_fmt + b'16f', op2.size)
    struct1 = Struct(fmt)
    #if op2.sort_method == 1:
        #print('**nelements =', nelements, nelements*2, obj.ntotal)
    #else:
        #print('**ntimes =', nelements)
    for unused_i in range(nelements):
        edata = data[n:n+ntotal]
        out = struct1.unpack(edata)

        (eid_device,
         fd1, sx1, sy1, txy1, angle1, major1, minor1, max_shear1,
         fd2, sx2, sy2, txy2, angle2, major2, minor2, max_shear2) = out

        eid, dt = get_eid_dt_from_eid_device(
            eid_device, op2.nonlinear_factor, op2.sort_method)
        if op2.is_debug_file:
            op2.binary_debug.write('  eid=%i C=[%s]\n' % (
                eid, ', '.join(['%r' % di for di in out])))

        add_new_eid_sort_x(
            dt, eid, cen,
            fd1, sx1, sy1, txy1, angle1, major1, minor1, max_shear1,
            fd2, sx2, sy2, txy2, angle2, major2, minor2, max_shear2)
        n += ntotal
    return n

def oes_ctria3_complex_vm_17(op2: OP2,
                             data: bytes,
                             obj: Union[ComplexPlateVMStressArray, ComplexPlateVMStrainArray],
                             nelements: int, ntotal: int,
                             dt,
                             is_magnitude_phase: bool) -> int:
    """
    #                   C O M P L E X   S T R E S S E S   I N   T R I A N G U L A R   E L E M E N T S   ( T R I A 3 )
    #                                                          (REAL/IMAGINARY)
    #
    #  ELEMENT       FIBER                                     - STRESSES IN ELEMENT  COORDINATE SYSTEM -
    #    ID.        DISTANCE              NORMAL-X                       NORMAL-Y                      SHEAR-XY               VON MISES
    #0       1  -4.359080E+00  -1.391918E+00 /  2.474756E-03  -1.423926E+00 /  2.530494E-03   2.655153E-02 / -5.158625E-05   1.408948E+00
    #            4.359080E+00   1.391918E+00 / -2.474756E-03   1.423926E+00 / -2.530494E-03  -2.655153E-02 /  5.158625E-05   1.408948E+00

    """
    n = 0
    cen = 0
    struct1 = Struct(mapfmt(op2._endian + op2._analysis_code_fmt + b'16f', op2.size))
    add_sort_x = getattr(obj, 'add_sort' + str(op2.sort_method))

    for unused_i in range(nelements):
        edata = data[n:n + ntotal]
        out = struct1.unpack(edata)
        (eid_device, fd1, oxx1r, oxx1i, oyy1r, oyy1i, txy1r, txy1i, ovm1,
                     fd2, oxx2r, oxx2i, oyy2r, oyy2i, txy2r, txy2i, ovm2, ) = out

        eid, dt = get_eid_dt_from_eid_device(
            eid_device, op2.nonlinear_factor, op2.sort_method)

        if is_magnitude_phase:
            oxx1 = polar_to_real_imag(oxx1r, oxx1i)
            oxx2 = polar_to_real_imag(oxx2r, oxx2i)
            oyy1 = polar_to_real_imag(oyy1r, oyy1i)
            oyy2 = polar_to_real_imag(oyy2r, oyy2i)
            txy1 = polar_to_real_imag(txy1r, txy1i)
            txy2 = polar_to_real_imag(txy2r, txy2i)
        else:
            oxx1 = complex(oxx1r, oxx1i)
            oxx2 = complex(oxx2r, oxx2i)
            oyy1 = complex(oyy1r, oyy1i)
            oyy2 = complex(oyy2r, oyy2i)
            txy1 = complex(txy1r, txy1i)
            txy2 = complex(txy2r, txy2i)
        #print(dt, eid, cen, sx1, sy1, txy1, max_shear1)
        #print(dt, eid, cen, sx2, sy2, txy2, max_shear2)
        add_sort_x(dt, eid, cen,
                   fd1, oxx1, oyy1, txy1, ovm1,
                   fd2, oxx2, oyy2, txy2, ovm2)
        n += ntotal

    #msg = op2.code_information()
    #self.log.warning(f'  skipping {op2.table_name_str} {op2.element_name}')
    #msg = '%s-%s' % (op2.table_name_str, op2.element_name)
    return n


def oes_ctria3_random_9(op2: OP2, data: bytes,
                        obj: Union[RandomPlateStressArray, RandomPlateStrainArray],
                        nelements: int, ntotal: int) -> int:
    n = 0
    struct1 = Struct(op2._endian + op2._analysis_code_fmt + b'8f')
    cen = 0 # CEN/4
    add_sort_x = getattr(obj, 'add_sort' + str(op2.sort_method))
    #print('ctria3_9 - SORT1')
    for unused_i in range(nelements):
        edata = data[n:n+ntotal]
        out = struct1.unpack(edata)

        (eid_device,
         fd1, sx1, sy1, txy1,
         fd2, sx2, sy2, txy2,) = out

        eid, dt = get_eid_dt_from_eid_device(
            eid_device, op2.nonlinear_factor, op2.sort_method)

        if op2.is_debug_file:
            op2.binary_debug.write('  eid=%i C=[%s]\n' % (
                eid, ', '.join(['%r' % di for di in out])))

        add_sort_x(dt, eid, cen,
                   fd1, sx1, sy1, txy1,
                   fd2, sx2, sy2, txy2)
        n += ntotal
    return n


def oes_cquad4_33_random_vm_11(op2: OP2, data: bytes,
                               obj: Union[RandomPlateVMStressArray, RandomPlateVMStrainArray],
                               nelements: int, ntotal: int) -> int:
    struct1 = Struct(op2._endian + op2._analysis_code_fmt + b'10f')
    cen = 0 # CEN/4
    n = 0
    add_sort_x = getattr(obj, 'add_sort' + str(op2.sort_method))
    for unused_i in range(nelements):
        edata = data[n:n+ntotal]
        out = struct1.unpack(edata)

        (eid_device,
         fd1, sx1, sy1, txy1, ovm1,
         fd2, sx2, sy2, txy2, ovm2,) = out

        eid, dt = get_eid_dt_from_eid_device(
            eid_device, op2.nonlinear_factor, op2.sort_method)
        #self.log.info(f'SORT1 VM: eid={eid} dt={dt}')

        if op2.is_debug_file:
            op2.binary_debug.write('  eid=%i C=[%s]\n' % (
                eid, ', '.join(['%r' % di for di in out])))

        add_sort_x(dt, eid, cen,
                   fd1, sx1, sy1, txy1, ovm1,
                   fd2, sx2, sy2, txy2, ovm2)
        n += ntotal
    return n


def oes_ctria3_random_vm_11(op2: OP2, data: bytes,
                            obj: Union[RandomPlateStressArray, RandomPlateStrainArray],
                            nelements: int, ntotal: int) -> int:
    n = 0
    struct1 = Struct(op2._endian + op2._analysis_code_fmt + b'10f')

    add_sort_x = getattr(obj, 'add_sort' + str(op2.sort_method))
    cen = 0 # CEN/4
    for unused_i in range(nelements):
        edata = data[n:n+ntotal]
        out = struct1.unpack(edata)

        (eid_device,
         fd1, sx1, sy1, txy1, ovm1,
         fd2, sx2, sy2, txy2, ovm2,) = out
        #print('CTRIA3', op2.nonlinear_factor, out)

        eid, dt = get_eid_dt_from_eid_device(
            eid_device, op2.nonlinear_factor, op2.sort_method)
        #self.log.info(f'SORT1 VM: eid={eid} dt={dt}')

        if op2.is_debug_file:
            op2.binary_debug.write('  eid=%i C=[%s]\n' % (
                eid, ', '.join(['%r' % di for di in out])))

        add_sort_x(dt, eid, cen,
                  fd1, sx1, sy1, txy1, ovm1,
                  fd2, sx2, sy2, txy2, ovm2)
        n += ntotal
    return n


def oes_cquad4_144_real(op2: OP2, data: bytes, ndata: int,
                        obj: RealPlateStrainArray,
                        nelements: int, nnodes: int, dt: Any) -> int:
    n = 0
    if op2.size == 4:
        center_format = op2._endian + op2._analysis_code_fmt + b'4si16f'
        node_format = op2._endian + b'i16f'
    else:
        center_format = op2._endian + mapfmt(op2._analysis_code_fmt, op2.size) + b'8sq16d'
        node_format = op2._endian + b'q16d'
    cs = Struct(center_format)
    ns = Struct(node_format)

    add_new_eid_sort_x = getattr(obj, 'add_new_eid_sort' + str(op2.sort_method))
    add_sort_x = getattr(obj, 'add_sort' + str(op2.sort_method))
    if op2.is_debug_file:
        op2.binary_debug.write(
            '  [cap, element1, element2, ..., cap]\n'
            '  cap = %i  # assume 1 cap when there could have been multiple\n'
            '  #elementi = [centeri, node1i, node2i, node3i, node4i]\n'
            '  #centeri = [eid_device, j, grid, fd1, sx1, sy1, txy1, angle1, major1, minor1, vm1,\n'
            '  #                                fd2, sx2, sy2, txy2, angle2, major2, minor2, vm2,)]\n'
            '  #nodeji = [grid, fd1, sx1, sy1, txy1, angle1, major1, minor1, vm1,\n'
            '  #                fd2, sx2, sy2, txy2, angle2, major2, minor2, vm2,)]\n'
            '  nelements=%i; nnodes=%i # +1 centroid\n' % (ndata, nelements, nnodes))

    grid_center = 0
    n76 = 76 * op2.factor
    n68 = 68 * op2.factor
    for unused_i in range(nelements):
        edata = data[n:n+n76]

        out = cs.unpack(edata)  # len=17*4
        # j is the number of nodes, so CQUAD4 -> 4, but we don't need to save it...
        (eid_device, unused_j,
         grid,
         fd1, sx1, sy1, txy1, angle1, major1, minor1, vm1,
         fd2, sx2, sy2, txy2, angle2, major2, minor2, vm2,) = out
        eid, dt = get_eid_dt_from_eid_device(
            eid_device, op2.nonlinear_factor, op2.sort_method)
        #print(f'tri eid={eid} dt={dt} nid={grid_center} unused_j={unused_j} grid={grid}')
        #print(out[:3])
        if op2.is_debug_file:
            op2.binary_debug.write('  eid=%i; C=[%s]\n' % (eid, ', '.join(['%r' % di for di in out])))
        add_new_eid_sort_x(dt, eid, grid_center,
                           fd1, sx1, sy1, txy1, angle1, major1, minor1, vm1,
                           fd2, sx2, sy2, txy2, angle2, major2, minor2, vm2)
        n += n76
        for inode in range(nnodes):
            out = ns.unpack(data[n:n + n68])
            (grid,
             fd1, sx1, sy1, txy1, angle1, major1, minor1, vm1,
             fd2, sx2, sy2, txy2, angle2, major2, minor2, vm2,) = out

            if op2.is_debug_file:
                d = tuple([grid,
                           fd1, sx1, sy1, txy1, angle1, major1, minor1, vm1,
                           fd2, sx2, sy2, txy2, angle2, major2, minor2, vm2])
                op2.binary_debug.write('  node%i = [%s]\n' % (inode+1, ', '.join(['%r' % di for di in d])))
            assert isinstance(grid, int), grid
            assert grid > 0, grid

            #print(f' {inode}: eid={eid} dt={dt} nid={grid}')
            add_sort_x(dt, eid, grid,
                       fd1, sx1, sy1, txy1, angle1, major1, minor1, vm1,
                       fd2, sx2, sy2, txy2, angle2, major2, minor2, vm2)
            n += n68
    return n


def oes_cquad4_144_random(op2: OP2, data: bytes,
                          obj: Union[RandomPlateStressArray, RandomPlateStrainArray],
                          nelements: int, nnodes: int, ndata: int) -> int:
    n = 0
    center_format = op2._endian + op2._analysis_code_fmt + b'4s i8f'
    node_format = op2._endian + b'i8f'
    cs = Struct(center_format)
    ns = Struct(node_format)

    if op2.is_debug_file:
        op2.binary_debug.write(
            '  [cap, element1, element2, ..., cap]\n'
            '  cap = %i  # assume 1 cap when there could have been multiple\n'
            '  #elementi = [centeri, node1i, node2i, node3i, node4i]\n'
            '  #centeri = [eid_device, j, grid, fd1, sx1, sy1, txy1,\n'
            '  #                                fd2, sx2, sy2, txy2,)]\n'
            '  #nodeji = [grid, fd1, sx1, sy1, txy1,\n'
            '  #                fd2, sx2, sy2, txy2,)]\n'
            '  nelements=%i; nnodes=%i # +1 centroid\n' % (ndata, nelements, nnodes))

    grid_center = 0
    add_sort_x = getattr(obj, 'add_sort' + str(op2.sort_method))

    ntotal1 = 44 * op2.factor
    ntotal2 = 36 * op2.factor
    for unused_i in range(nelements):
        edata = data[n:n+ntotal1]
        #self.show_data(edata)
        out = cs.unpack(edata)  # len=17*4
        #print(out)
        # j is the number of nodes, so CQUAD4 -> 4, but we don't need to save it...
        eid_device = out[0]
        (eid_device, unused_j,
         grid,
         fd1, sx1, sy1, txy1,
         fd2, sx2, sy2, txy2,) = out
        #print(out)

        eid, dt = get_eid_dt_from_eid_device(
            eid_device, op2.nonlinear_factor, op2.sort_method)

        if op2.is_debug_file:
            op2.binary_debug.write('  eid=%i; C=[%s]\n' % (eid, ', '.join(['%r' % di for di in out])))

        add_sort_x(dt, eid, grid_center,
                   fd1, sx1, sy1, txy1, fd2, sx2, sy2, txy2)
        n += ntotal1
        for inode in range(nnodes):
            out = ns.unpack(data[n:n + ntotal2])
            #print(out)
            (grid,
             fd1, sx1, sy1, txy1,
             fd2, sx2, sy2, txy2,) = out
            #print(out)
            if op2.is_debug_file:
                d = tuple([grid,
                           fd1, sx1, sy1, txy1,
                           fd2, sx2, sy2, txy2])
                op2.binary_debug.write('  node%i = [%s]\n' % (inode+1, ', '.join(['%r' % di for di in d])))
            assert isinstance(grid, int), grid
            assert grid > 0, grid

            # leaving off grid
            add_sort_x(dt, eid, grid,
                       fd1, sx1, sy1, txy1,
                       fd2, sx2, sy2, txy2)
            n += ntotal2
    return n


def oes_cbar_real_16(op2: OP2, data: bytes,
                     obj: Union[RealBarStressArray, RealBarStrainArray],
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


def oes_weldp_msc_real_8(op2: OP2, data: bytes,
                         obj: Union[int, float],
                         nelements: int, ntotal: int, dt: Any) -> int:
    #'    ELEMENT          AXIAL         MAX  STRESS      MIN  STRESS      MAX  STRESS      MIN  STRESS        MAXIMUM          BEARING '
    #'      ID             STRESS           END-A            END-A            END-B            END-B        SHEAR  STRESS       STRESS'
    #'        179      -3.153108E+00     8.089753E+02    -8.152815E+02     7.946552E+02    -8.009614E+02     2.852777E+01     1.179798E+01'

    n = 0
    fmt = mapfmt(op2._endian + op2._analysis_code_fmt + b'7f', op2.size)
    struct1 = Struct(fmt)
    add_sort_x = getattr(obj, 'add_sort' + str(op2.sort_method))
    for i in range(nelements):
        edata = data[n:n+ntotal]
        out = struct1.unpack(edata)
        (eid_device,
         axial, maxa, mina, maxb, minb, max_shear, bearing) = out
        eid, dt = get_eid_dt_from_eid_device(
            eid_device, op2.nonlinear_factor, op2.sort_method)
        if op2.is_debug_file:
            op2.binary_debug.write('  eid=%i; C%i=[%s]\n' % (
                eid, i, ', '.join(['%r' % di for di in out])))
        n += ntotal
        add_sort_x(
            dt, eid,
            axial, maxa, mina, maxb, minb, max_shear, bearing)

    if op2.is_sort2:
        #print(add_sort_x)
        #print(''.join(obj.get_stats()))
        #print(f'{self.table_name} sort_method={op2.sort_method}', obj._times)
        assert len(np.unique(obj._times)) == len(obj._times), obj._times.tolist()
    return n


def oes_weldp_msc_complex_15(op2: OP2,
                             data: bytes,
                             obj: Union[ComplexWeldStressArray, ComplexWeldStrainArray],
                             nelements: int, ntotal: int,
                             is_magnitude_phase: bool, dt: Any) -> int:
    n = 0
    struct1 = Struct(op2._endian + mapfmt(op2._analysis_code_fmt + b'14f', op2.size))
    add_sort_x = getattr(obj, 'add_sort' + str(op2.sort_method))
    for unused_i in range(nelements):
        edata = data[n:n + ntotal]
        out = struct1.unpack(edata)  # num_wide=13
        if op2.is_debug_file:
            op2.binary_debug.write('WELD-118 - %s\n' % str(out))
        (eid_device,
         raxial, rmax_a, rmin_a, rmax_b, rmin_b, rmax_shear, rbearing,
         iaxial, imax_a, imin_a, imax_b, imin_b, imax_shear, ibearing) = out
        eid, dt = get_eid_dt_from_eid_device(
            eid_device, op2.nonlinear_factor, op2.sort_method)
        assert eid > 0

        if is_magnitude_phase:
            axial = polar_to_real_imag(raxial, iaxial)
            max_a = polar_to_real_imag(rmax_a, imax_a)
            min_a = polar_to_real_imag(rmin_a, imin_a)

            max_b = polar_to_real_imag(rmax_b, imax_b)
            min_b = polar_to_real_imag(rmin_b, imin_b)
            max_shear = polar_to_real_imag(rmax_shear, imax_shear)
            bearing = polar_to_real_imag(rbearing, ibearing)
        else:
            axial = complex(raxial, iaxial)
            max_a = complex(rmax_a, imax_a)
            min_a = complex(rmin_a, imin_a)

            max_b = complex(rmax_b, imax_b)
            min_b = complex(rmin_b, imin_b)
            max_shear = complex(rmax_shear, imax_shear)
            bearing = complex(rbearing, ibearing)
        #print(out)
        add_sort_x(dt, eid, axial, max_a, min_a, max_b, min_b, max_shear, bearing)
        n += ntotal
    return n


def oes_fastp_msc_real_7(op2: OP2, data: bytes,
                         obj: Union[int, float],
                         nelements: int, ntotal: int, dt: Any) -> int:
    n = 0
    fmt = mapfmt(op2._endian + op2._analysis_code_fmt + b'6f', op2.size)
    struct1 = Struct(fmt)
    add_sort_x = getattr(obj, 'add_sort' + str(op2.sort_method))
    for i in range(nelements):
        edata = data[n:n+ntotal]
        out = struct1.unpack(edata)
        (eid_device,
         force_x, force_y, force_z, moment_x, moment_y, moment_z) = out
        eid, dt = get_eid_dt_from_eid_device(
            eid_device, op2.nonlinear_factor, op2.sort_method)
        if op2.is_debug_file:
            op2.binary_debug.write('  eid=%i; C%i=[%s]\n' % (
                eid, i, ', '.join(['%r' % di for di in out])))
        n += ntotal
        add_sort_x(
            dt, eid,
            force_x, force_y, force_z, moment_x, moment_y, moment_z)

    if op2.is_sort2:
        #print(add_sort_x)
        #print(''.join(obj.get_stats()))
        #print(f'{self.table_name} sort_method={op2.sort_method}', obj._times)
        assert len(np.unique(obj._times)) == len(obj._times), obj._times.tolist()
    return n


def oes_fastp_msc_complex_13(op2: OP2,
                             data: bytes,
                             obj: Union[ComplexFastStressArray, ComplexFastStrainArray],
                             nelements: int, ntotal: int,
                             is_magnitude_phase: bool, dt: Any) -> int:
    n = 0
    struct1 = Struct(op2._endian + mapfmt(op2._analysis_code_fmt + b'12f', op2.size))
    add_sort_x = getattr(obj, 'add_sort' + str(op2.sort_method))
    for unused_i in range(nelements):
        edata = data[n:n + ntotal]
        out = struct1.unpack(edata)  # num_wide=13
        if op2.is_debug_file:
            op2.binary_debug.write('FASTP-118 - %s\n' % str(out))
        (eid_device,
         rforce_x, rforce_y, rforce_z, rmoment_x, rmoment_y, rmoment_z,
         iforce_x, iforce_y, iforce_z, imoment_x, imoment_y, imoment_z) = out
        eid, dt = get_eid_dt_from_eid_device(
            eid_device, op2.nonlinear_factor, op2.sort_method)
        assert eid > 0

        if is_magnitude_phase:
            force_x = polar_to_real_imag(rforce_x, iforce_x)
            force_y = polar_to_real_imag(rforce_y, iforce_y)
            force_z = polar_to_real_imag(rforce_z, iforce_z)

            moment_x = polar_to_real_imag(rmoment_x, imoment_x)
            moment_y = polar_to_real_imag(rmoment_y, imoment_y)
            moment_z = polar_to_real_imag(rmoment_z, imoment_z)
        else:
            force_x = complex(rforce_x, iforce_x)
            force_y = complex(rforce_y, iforce_y)
            force_z = complex(rforce_z, iforce_z)

            moment_x = complex(rmoment_x, imoment_x)
            moment_y = complex(rmoment_y, imoment_y)
            moment_z = complex(rmoment_z, imoment_z)
        add_sort_x(dt, eid, force_x, force_y, force_z, moment_x, moment_y, moment_z)
        n += ntotal
    return n


def oes_cshear_real_4(op2: OP2, data: bytes,
                      obj: Union[RealShearStressArray, RealShearStrainArray],
                      ntotal: int, nelements: int, dt: Any) -> int:
    n = 0
    struct1 = Struct(op2._endian + mapfmt(op2._analysis_code_fmt + b'3f', op2.size))
    add_sort_x = getattr(obj, 'add_sort' + str(op2.sort_method))
    for unused_i in range(nelements):
        edata = data[n:n + ntotal]
        out = struct1.unpack(edata)  # num_wide=5
        if op2.is_debug_file:
            op2.binary_debug.write('CSHEAR-4 - %s\n' % str(out))

        (eid_device, max_strain, avg_strain, margin) = out
        #print(eid_device, max_strain, avg_strain, margin)
        eid, dt = get_eid_dt_from_eid_device(
            eid_device, op2.nonlinear_factor, op2.sort_method)
        add_sort_x(dt, eid, max_strain, avg_strain, margin)
        n += ntotal
    return n

def oes_cshear_complex_5(op2: OP2,
                         data: bytes,
                         obj: Union[ComplexShearStressArray, ComplexShearStrainArray],
                         nelements: int, ntotal: int,
                         is_magnitude_phase: bool) -> int:
    n = 0
    struct1 = Struct(op2._endian + mapfmt(op2._analysis_code_fmt + b'4f', op2.size))
    add_sort_x = getattr(obj, 'add_sort' + str(op2.sort_method))
    for unused_i in range(nelements):
        edata = data[n:n + ntotal]
        out = struct1.unpack(edata)  # num_wide=5
        if op2.is_debug_file:
            op2.binary_debug.write('CSHEAR-4 - %s\n' % str(out))
        (eid_device, etmaxr, etmaxi, etavgr, etavgi) = out
        eid, dt = get_eid_dt_from_eid_device(
            eid_device, op2.nonlinear_factor, op2.sort_method)

        if is_magnitude_phase:
            etmax = polar_to_real_imag(etmaxr, etmaxi)
            etavg = polar_to_real_imag(etavgr, etavgi)
        else:
            etmax = complex(etmaxr, etmaxi)
            etavg = complex(etavgr, etavgi)
        add_sort_x(dt, eid, etmax, etavg)
        n += ntotal
    return n


def oes_cbar_complex_19(op2: OP2,
                        data: bytes,
                        obj: Union[ComplexBarStressArray, ComplexBarStrainArray],
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
        if is_magnitude_phase:
            s1a = polar_to_real_imag(s1ar, s1ai)
            s1b = polar_to_real_imag(s1br, s1bi)
            s2a = polar_to_real_imag(s2ar, s2ai)
            s2b = polar_to_real_imag(s2br, s2bi)
            s3a = polar_to_real_imag(s3ar, s3ai)
            s3b = polar_to_real_imag(s3br, s3bi)
            s4a = polar_to_real_imag(s4ar, s4ai)
            s4b = polar_to_real_imag(s4br, s4bi)
            axial = polar_to_real_imag(axialr, axiali)
        else:
            s1a = complex(s1ar, s1ai)
            s1b = complex(s1br, s1bi)
            s2a = complex(s2ar, s2ai)
            s2b = complex(s2br, s2bi)
            s3a = complex(s3ar, s3ai)
            s3b = complex(s3br, s3bi)
            s4a = complex(s4ar, s4ai)
            s4b = complex(s4br, s4bi)
            axial = complex(axialr, axiali)

        obj.add_new_eid_sort1(dt, eid,
                              s1a, s2a, s3a, s4a, axial,
                              s1b, s2b, s3b, s4b)
    return n


def oes_cbar_random_10(op2: OP2, data: bytes,
                       obj: Union[RandomBarStressArray, RandomBarStrainArray],
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


def oes_cbush_real_7(op2: OP2, data: bytes,
                     obj: Union[RealBushStressArray, RealBushStrainArray],
                     nelements: int, ntotal: int, dt, debug: bool=False) -> int:
    n = 0
    struct1 = Struct(op2._endian + mapfmt(op2._analysis_code_fmt + b'6f', op2.size))
    add_sort_x = getattr(obj, 'add_sort' + str(op2.sort_method))
    nonlinear_factor = op2.nonlinear_factor
    #print(add_sort_x)
    #print('obj.is_sort1 =', obj.is_sort1, obj.table_name)
    for unused_i in range(nelements):
        edata = data[n:n + ntotal]
        n += ntotal

        out = struct1.unpack(edata)  # num_wide=7
        if op2.is_debug_file:
            op2.binary_debug.write('CBUSH-102 - %s\n' % str(out))

        (eid_device, tx, ty, tz, rx, ry, rz) = out
        eid, dt = get_eid_dt_from_eid_device(
            eid_device, nonlinear_factor, op2.sort_method)
        if debug:  # pragma: no cover
            print(f'CBUSH: eid_device={eid_device} eid={eid} dt={nonlinear_factor} nf={nonlinear_factor} -> {obj.data.shape}')

        add_sort_x(dt, eid, tx, ty, tz, rx, ry, rz)
    return n


def oes_cbush_complex_13(op2: OP2,
                         data: bytes,
                         obj: Union[ComplexCBushStressArray, ComplexCBushStrainArray],
                         nelements: int, ntotal: int,
                         is_magnitude_phase: bool) -> int:
    n = 0
    struct1 = Struct(op2._endian + mapfmt(op2._analysis_code_fmt + b'12f', op2.size))
    add_sort_x = getattr(obj, 'add_sort' + str(op2.sort_method))
    for unused_i in range(nelements):
        edata = data[n:n + ntotal]
        out = struct1.unpack(edata)  # num_wide=7
        if op2.is_debug_file:
            op2.binary_debug.write('CBUSH-102 - %s\n' % str(out))

        (eid_device,
         txr, tyr, tzr, rxr, ryr, rzr,
         txi, tyi, tzi, rxi, ryi, rzi) = out
        eid, dt = get_eid_dt_from_eid_device(
            eid_device, op2.nonlinear_factor, op2.sort_method)

        if is_magnitude_phase:
            tx = polar_to_real_imag(txr, txi)
            ty = polar_to_real_imag(tyr, tyi)
            tz = polar_to_real_imag(tzr, tzi)
            rx = polar_to_real_imag(rxr, rxi)
            ry = polar_to_real_imag(ryr, ryi)
            rz = polar_to_real_imag(rzr, rzi)
        else:
            tx = complex(txr, txi)
            ty = complex(tyr, tyi)
            tz = complex(tzr, tzi)
            rx = complex(rxr, rxi)
            ry = complex(ryr, ryi)
            rz = complex(rzr, rzi)
        add_sort_x(dt, eid, tx, ty, tz, rx, ry, rz)
        n += ntotal
    return ntotal

def oes_ctriax6_real_33(op2: OP2, data: bytes,
                        obj: Union[RealTriaxStressArray, RealTriaxStrainArray],
                        nelements: int, ntotal: int, dt: Any) -> int:
    n = 0
    ntotal1 = 36 * op2.factor
    ntotal2 = 32 * op2.factor
    s1 = Struct(op2._endian + b'2i7f')  # 36
    s2 = Struct(op2._endian + b'i7f')
    add_sort_x = getattr(obj, 'add_sort' + str(op2.sort_method))
    for unused_i in range(nelements):
        out = s1.unpack(data[n:n + ntotal1])
        (eid_device, loc, rs, azs, As, ss, maxp, tmax, octs) = out
        if op2.is_debug_file:
            op2.binary_debug.write('CTRIAX6-53A - %s\n' % (str(out)))
        eid, dt = get_eid_dt_from_eid_device(
            eid_device, op2.nonlinear_factor, op2.sort_method)

        add_sort_x(dt, eid, loc, rs, azs, As, ss, maxp, tmax, octs)
        n += ntotal1
        for unused_j in range(3):
            out = s2.unpack(data[n:n + ntotal2])
            (loc, rs, azs, As, ss, maxp, tmax, octs) = out
            if op2.is_debug_file:
                op2.binary_debug.write('CTRIAX6-53B - %s\n' % (str(out)))
            add_sort_x(dt, eid, loc, rs, azs, As, ss, maxp, tmax, octs)
            n += ntotal2
    return n


def oes_ctriax_complex_37(op2: OP2,
                          data: bytes,
                          obj: ComplexTriaxStressArray,
                          nelements: int,
                          is_magnitude_phase: bool) -> int:
    n = 0
    ntotal1 = 40 * op2.factor
    ntotal2 = 36 * op2.factor
    s1 = Struct(op2._endian + mapfmt(op2._analysis_code_fmt + b'i8f', op2.size)) # 10*4 = 40
    s2 = Struct(op2._endian + mapfmt(b'i8f', op2.size))  #  9*4 = 36

    add_sort_x = getattr(obj, 'add_sort' + str(op2.sort_method))
    for unused_i in range(nelements):
        out = s1.unpack(data[n:n + ntotal1])
        (eid_device, loc, rsr, rsi, azsr, azsi, Asr, Asi, ssr, ssi) = out
        eid, dt = get_eid_dt_from_eid_device(
            eid_device, op2.nonlinear_factor, op2.sort_method)
        if op2.is_debug_file:
            op2.binary_debug.write('CTRIAX6-53 eid=%i\n    %s\n' % (eid, str(out)))
        #print('CTRIAX6-53 eid=%i\n    %s\n' % (eid, str(out)))

        if is_magnitude_phase:
            rs = polar_to_real_imag(rsr, rsi)
            azs = polar_to_real_imag(azsr, azsi)
            As = polar_to_real_imag(Asr, Asi)
            ss = polar_to_real_imag(ssr, ssi)
        else:
            rs = complex(rsr, rsi)
            azs = complex(azsr, azsi)
            As = complex(Asr, Asi)
            ss = complex(ssr, ssi)
        obj.add_element_sort1(dt, eid)
        add_sort_x(dt, eid, loc, rs, azs, As, ss)

        n += ntotal1
        for unused_j in range(3):
            out = s2.unpack(data[n:n + ntotal2])
            (loc, rsr, rsi, azsr, azsi, Asr, Asi, ssr, ssi) = out
            if op2.is_debug_file:
                op2.binary_debug.write('    %s\n' % (str(out)))
            #print("eid=%s loc=%s rs=%s azs=%s as=%s ss=%s" % (
                #eid, loc, rs, azs, As, ss))

            if is_magnitude_phase:
                rs = polar_to_real_imag(rsr, rsi)
                azs = polar_to_real_imag(azsr, azsi)
                As = polar_to_real_imag(Asr, Asi)
                ss = polar_to_real_imag(ssr, ssi)
            else:
                rs = complex(rsr, rsi)
                azs = complex(azsr, azsi)
                As = complex(Asr, Asi)
                ss = complex(ssr, ssi)
            add_sort_x(dt, eid, loc, rs, azs, As, ss)
            n += ntotal2  # 4*8
    return n


def oes_csolid_real(op2: OP2, data: bytes,
                    obj: Union[RealSolidStressArray, RealSolidStrainArray],
                    nelements: int, dt: Any,
                    element_name: str, nnodes_expected: int,
                    preline1: str, preline2: str) -> int:
    """
    reads stress/strain for element type:
     - 39 : CTETRA
     - 67 : CHEXA
     - 68 : CPENTA

    """
    add_eid_sort_x = getattr(obj, 'add_eid_sort' + str(op2.sort_method))
    add_node_sort_x = getattr(obj, 'add_node_sort' + str(op2.sort_method))
    n = 0
    if op2.size == 4:
        fmt1 = op2._endian + op2._analysis_code_fmt + b'i4si'
        fmt2 = op2._endian + b'i20f'
    else:
        fmt1 = op2._endian + mapfmt(op2._analysis_code_fmt, op2.size) + b'q8sq'
        fmt2 = op2._endian + b'q20d'
    struct1 = Struct(fmt1)
    struct2 = Struct(fmt2)
    if op2.is_debug_file:
        msg = '%s-%s nelements=%s nnodes=%s; C=[sxx, sxy, s1, a1, a2, a3, pressure, svm,\n' % (
            op2.element_name, op2.element_type, nelements, nnodes_expected)
        msg += '                                 syy, syz, s2, b1, b2, b3,\n'
        msg += '                                 szz, sxz, s3, c1, c2, c3]\n'
        op2.binary_debug.write(msg)

    n16 = 16 * op2.factor
    n84 = 84 * op2.factor
    #if op2.is_sort1:
        #print('nelements =', nelements)
    #else:
        #print('ntimes =', nelements)
    for unused_i in range(nelements):
        edata = data[n:n+n16]
        out = struct1.unpack(edata)
        (eid_device, cid, unused_abcd, nnodes) = out
        eid, dt = get_eid_dt_from_eid_device(
            eid_device, op2.nonlinear_factor, op2.sort_method)
        #print(f'eid_device={eid_device} eid={eid} dt={dt}')
        if op2.is_debug_file:
            op2.binary_debug.write('%s - eid=%i; %s\n' % (preline1, eid, str(out)))

        assert nnodes < 21, 'print_block(data[n:n+16])'  #self.print_block(data[n:n+16])

        n += n16
        for inode in range(nnodes_expected):  # nodes pts, +1 for centroid (???)
            out = struct2.unpack(data[n:n + n84]) # 4*21 = 84
            if op2.is_debug_file:
                op2.binary_debug.write('%s - %s\n' % (preline2, str(out)))
            (grid_device,
             sxx, sxy, s1, a1, a2, a3, pressure, svm,
             syy, syz, s2, b1, b2, b3,
             szz, sxz, s3, c1, c2, c3) = out

            if op2.is_debug_file:
                op2.binary_debug.write('  eid=%s inode=%i; C=[%s]\n' % (
                    eid, grid_device, ', '.join(['%r' % di for di in out])))

            #if grid_device == 0:
                #grid = 'CENTER'
            #else:
                ##grid = (grid_device - device_code) // 10
                #grid = grid_device

            grid = grid_device
            a_cos = [a1, a2, a3]
            b_cos = [b1, b2, b3]
            c_cos = [c1, c2, c3]

            if inode == 0:
                #  this is correct, but fails
                #element_name = op2.element_name + str(nnodes)
                add_eid_sort_x(element_name, cid, dt, eid, grid,
                               sxx, syy, szz, sxy, syz, sxz, s1, s2, s3,
                               a_cos, b_cos, c_cos, pressure, svm)
            else:
                add_node_sort_x(dt, eid, inode, grid,
                                sxx, syy, szz, sxy, syz, sxz, s1, s2, s3,
                                a_cos, b_cos, c_cos, pressure, svm)
            n += n84
    return n


def oes_csolid_complex(op2: OP2, data: bytes,
                       obj: Union[ComplexSolidStressArray, ComplexSolidStrainArray],
                       nelements: int, # nnodes: int,
                       element_name: str, nnodes_expected: int,
                       is_magnitude_phase: bool) -> int:
    n = 0
    if op2.size == 4:
        s1 = Struct(op2._endian + b'2i4si')
    else:
        s1 = Struct(op2._endian + b'2q8sq')
    s2 = Struct(op2._endian + mapfmt(b'i12f', op2.size))
    ntotal1 = 16 * op2.factor
    ntotal2 = 52 * op2.factor
    for unused_i in range(nelements):
        edata = data[n:n+ntotal1]
        n += ntotal1
        out = s1.unpack(edata)
        (eid_device, cid, ctype, nodef) = out
        eid, dt = get_eid_dt_from_eid_device(
            eid_device, op2.nonlinear_factor, op2.sort_method)
        if op2.is_debug_file:
            op2.binary_debug.write('  eid=%i C=[%s]\n' % (
                eid, ', '.join(['%r' % di for di in out])))

        #element_name = op2.element_name + str(nodef)  # this is correct, but has problems...
        obj.add_eid_sort1(op2.element_type, element_name, dt, eid, cid, ctype, nodef)
        for inode in range(nnodes_expected):
            edata = data[n:n+ntotal2]
            n += ntotal2
            out = s2.unpack(edata)
            (grid,
             exr, eyr, ezr, etxyr, etyzr, etzxr,
             exi, eyi, ezi, etxyi, etyzi, etzxi) = out
            #if grid == 0:
                #grid = 'CENTER'

            if is_magnitude_phase:
                ex = polar_to_real_imag(exr, exi)
                ey = polar_to_real_imag(eyr, eyi)
                ez = polar_to_real_imag(ezr, ezi)
                etxy = polar_to_real_imag(etxyr, etxyi)
                etyz = polar_to_real_imag(etyzr, etyzi)
                etzx = polar_to_real_imag(etzxr, etzxi)
            else:
                ex = complex(exr, exi)
                ey = complex(eyr, eyi)
                ez = complex(ezr, ezi)
                etxy = complex(etxyr, etxyi)
                etyz = complex(etyzr, etyzi)
                etzx = complex(etzxr, etzxi)

            if op2.is_debug_file:
                op2.binary_debug.write('       node%s=[%s]\n' % (
                    grid, ', '.join(['%r' % di for di in out])))
            obj.add_node_sort1(dt, eid, grid, inode,
                               ex, ey, ez, etxy, etyz, etzx)
    return n


def oes_csolid_random(op2: OP2, data: bytes,
                      obj: Union[RandomSolidStressArray, RandomSolidStrainArray],
                      nelements: int,
                      element_name: str, nnodes_expected: int,
                      preline1: str, preline2: str) -> int:
    n = 0
    struct1 = Struct(op2._endian + op2._analysis_code_fmt + b'i4si')
    struct2 = Struct(op2._endian + b'i6f')
    if op2.is_debug_file and 0:
        msg = '%s-%s nelements=%s nnodes=%s; C=[sxx, sxy, s1, a1, a2, a3, pressure, svm,\n' % (
            op2.element_name, op2.element_type, nelements, nnodes_expected)
        msg += '                                 syy, syz, s2, b1, b2, b3,\n'
        msg += '                                 szz, sxz, s3, c1, c2, c3]\n'
        op2.binary_debug.write(msg)

    ntotal1 = 16 * op2.factor
    ntotal2 = 28 * op2.factor
    for unused_i in range(nelements):
        edata = data[n:n+ntotal1]
        out = struct1.unpack(edata)
        (eid_device, cid, unused_abcd, grid) = out

        eid, dt = get_eid_dt_from_eid_device(
            eid_device, op2.nonlinear_factor, op2.sort_method)
        assert eid > 0, eid

        if op2.is_debug_file and 0:
            op2.binary_debug.write('%s - eid=%i; %s\n' % (preline1, eid, str(out)))

        #assert nnodes < 21, 'print_block(data[n:n+16])'  #self.print_block(data[n:n+16])

        n += ntotal1
        for inode in range(nnodes_expected):  # nodes pts, +1 for centroid (???)
            #self.show_data(data[n:n+48])
            out = struct2.unpack(data[n:n + ntotal2]) # 4*7 = 28
            if op2.is_debug_file:
                op2.binary_debug.write('%s - %s\n' % (preline2, str(out)))
            (grid_device, sxx, syy, szz, txy, tyz, txz) = out

            if op2.is_debug_file:
                op2.binary_debug.write('  eid=%s inode=%i; C=[%s]\n' % (
                    eid, grid_device, ', '.join(['%r' % di for di in out])))

            #if grid_device == 0:
                #grid = 'CENTER'
            #else:
                ##grid = (grid_device - device_code) // 10
                #grid = grid_device

            grid = grid_device
            if inode == 0:
                #  this is correct, but fails
                #element_name = op2.element_name + str(nnodes)
                obj.add_eid_sort1(element_name, cid, dt, eid, grid,
                                  sxx, syy, szz, txy, tyz, txz)
            else:
                obj.add_node_sort1(dt, eid, inode, grid,
                                   sxx, syy, szz, txy, tyz, txz)
            n += ntotal2
    return n


def oes_comp_shell_real_11(op2: OP2, data: bytes, ndata: int,
                           obj: Union[RealCompositePlateStressArray, RealCompositePlateStrainArray],
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


def oesrt_comp_shell_real_9(op2: OP2, data: bytes, ndata: int,
                            obj: RealCompositePlateStressStrengthRatioArray,
                            ntotal: int, nelements: int, etype: str, dt: Any) -> int:
    """
    # (1, b'TSAI-WU ', 1, 9.118781089782715, 1.3563156426940112e-19, 310.244140625, nan, 1.3563156426940112e-19)
    # (-1, b'        ', 2, 11.119629859924316, 1.3563156426940112e-19, 202.90777587890625, nan, 1.3563156426940112e-19)
    # (-1, b'        ', 3, 9.412755012512207, 1.3563156426940112e-19, 124.20950317382812, nan, 1.3563156426940112e-19)
    # (-1, b'        ', 4, 9.56696605682373, 1.3563156426940112e-19, 89.90276336669922, nan, 1.3563156426940112e-19)
    # (-1, b'        ', 5, 11.674513816833496, 1.3563156426940112e-19, 78.44961547851562, nan, 1.3563156426940112e-19)
    # (-1, b'        ', 6, 9.891059875488281, 1.3563156426940112e-19, 63.61910629272461, nan, 1.3563156426940112e-19)
    # (-1, b'        ', 7, 8.669649959069946e+22, 1.3563156426940112e-19, 63.619102478027344, nan, 1.3563156426940112e-19)
    # (-1, b'        ', 8, 1.0372562759530309e+23, 1.3563156426940112e-19, 63.61910629272461, nan, 1.3563156426940112e-19)
    # (-1, b'        ', 9, 26.741111755371094, 1.3563156426940112e-19, 78.44961547851562, nan, 1.3563156426940112e-19)
    # (-1, b'        ', 10, 48.35708236694336, 1.3563156426940112e-19, 89.90276336669922, nan, 1.3563156426940112e-19)
    # (-1, b'        ', 11, 25.865455627441406, 1.3563156426940112e-19, 124.20950317382812, nan, 1.3563156426940112e-19)
    # (-1, b'        ', 12, 25.44878578186035, 1.3563156426940112e-19, 202.90777587890625, nan, 1.3563156426940112e-19)
    # (-1, b'        ', 13, 46.05017852783203, 1.3563156426940112e-19, 310.244140625, nan, 1.3563156426940112e-19)
    # (-1, b'        ', 14, 24.654462814331055, 1.3563156426940112e-19, nan, 9.118781089782715, 1.3563156426940112e-19)
    # (2, b'TSAI-WU ', 1, 12.055706024169922, 1.3563156426940112e-19, 310.244140625, nan, 1.3563156426940112e-19)
    """
    n = 0
    eid_old = 0
                                                                # ft ply  sr,fi,sr
    structs = Struct(op2._endian + op2._analysis_code_fmt + b' 8s  i    f  4s  f  f 4s')
    structf = Struct(op2._endian + op2._analysis_code_fmt + b' 8s  i    f  f   f  f 4s')
    sort_method = op2.sort_method
    assert sort_method == 1, 'oesrt_comp_shell_real_9'
    add_eid_sort_x = getattr(obj, 'add_eid_sort' + str(op2.sort_method))
    add_sort_x = getattr(obj, 'add_sort' + str(op2.sort_method))

    for unused_i in range(nelements):
        edata = data[n:n+ntotal]  # 4*9
        outs = structs.unpack(edata)

        # strings = (b'\x01\x00\x00\x00HFAIL   \x01\x00\x00\x00\x00\x00\xc8B\x00\x00\x00\x00\xff\xff\xff\xff\xff\xff\xff\xff    ',)
        #sr = 100.0; fi = 0.; sr = nan;
        #           eid ft         ply sr    fi     sr min  flag
        #ints    = (1, 'HFAIL   ', 1, 100.0, 0.0,   -1, -1, '    ')
        #floats  = (1, 'HFAIL   ', 1, 100.0, 0.0, nan, nan, '    ')
        (eid, failure_theory_bytes, ply_id, strength_ratio_ply, failure_index_bonding_bytes, strength_ratio_bonding,
         min_sr_bonding_fi_bonding, flag_bytes) = outs

        out = outs
        if failure_index_bonding_bytes != b'    ':
            outf = structf.unpack(edata)
            (eid, failure_theory_bytes, ply_id, strength_ratio_ply, failure_index_bonding, strength_ratio_bonding,
             min_sr_bonding_fi_bonding, flag_bytes) = outf
            out = outf
        else:
            failure_index_bonding = np.nan
        failure_theory = failure_theory_bytes.rstrip().decode('latin1')
        flag = flag_bytes.rstrip().decode('latin1')
        assert len(flag) <= 3, 'flag=%r' % flag

        #print(f'eid={eid}; ft={failure_theory!r}; sr_ply={strength_ratio_ply:g}; fi_bonding={failure_index_bonding!r} '
              #f'sr_bonding={strength_ratio_bonding}; min(sr_bonding, fi_bonding)={min_sr_bonding_fi_bonding:g} flag={flag!r}')
        #self.show_data(edata, types='ifs')

        #if op2.is_debug_file:
            #op2.binary_debug.write('  eid=%i; layer=%i; C=[%s]\n' % (eid, layer, ', '.join(['%r' % di for di in out])))
        if op2.is_debug_file:
            op2.binary_debug.write('%s-%s - (%s) + %s\n' % (op2.element_name, op2.element_type, eid, str(out)))
        #print(f'{op2.element_name}-{op2.element_type} - eid={eid}; ply_id={ply_id}; out={out}\n')

        #if ply_id == 1:  # not sufficient; ply_id may start at 101

        if eid == -1:
            add_sort_x(dt, eid_old, failure_theory, ply_id, strength_ratio_ply,
                       failure_index_bonding, strength_ratio_bonding, min_sr_bonding_fi_bonding, flag)
        #if eid != -1:
            # originally initialized to None, the buffer doesnt reset it, so it is the old value
            #add_eid_sort_x(etype, dt, eid, layer, o1, o2,
                           #t12, t1z, t2z, angle, major, minor, ovm)
        else:
            add_eid_sort_x(etype, dt, eid, failure_theory, ply_id, strength_ratio_ply,
                           failure_index_bonding, strength_ratio_bonding, min_sr_bonding_fi_bonding, flag)
            eid_old = eid
        n += ntotal
    return n


def oes_cshear_random_3(op2: OP2, data: bytes,
                        obj: Union[RandomShearStressArray, RandomShearStrainArray],
                        nelements: int, ntotal: int) -> int:
    n = 0
    add_sort_x = getattr(obj, 'add_sort' + str(op2.sort_method))
    struct1 = Struct(op2._endian + op2._analysis_code_fmt + b'2f')
    for unused_i in range(nelements):
        edata = data[n:n + ntotal]
        out = struct1.unpack(edata)  # num_wide=5
        if op2.is_debug_file:
            op2.binary_debug.write('CSHEAR-4 - %s\n' % str(out))

        (eid_device, max_strain, avg_strain) = out
        eid, dt = get_eid_dt_from_eid_device(
            eid_device, op2.nonlinear_factor, op2.sort_method)
        #print(eid, dt)
        add_sort_x(dt, eid, max_strain, avg_strain)
        n += ntotal
    return n


def oes_cbend_real_21(op2: OP2, data: bytes,
                      obj: Union[RealBendStressArray, RealBendStrainArray],
                      nelements: int, ntotal: int, dt) -> int:
    n = 0
    ntotali = 40 * op2.factor
    struct1 = Struct(op2._endian + op2._analysis_code_fmt)
    struct2 = Struct(op2._endian + b'i9f')
    add_sort_x = getattr(obj, 'add_sort' + str(op2.sort_method))

    #print('ntimes =', nelements)
    for unused_i in range(nelements):
        edata = data[n:n + 4]
        eid_device, = struct1.unpack(edata)
        eid, dt = get_eid_dt_from_eid_device(
            eid_device, op2.nonlinear_factor, op2.sort_method)
        n += 4

        for unused_j in range(2):
            edata = data[n:n + ntotali]
            out = struct2.unpack(edata)
            if op2.is_debug_file:
                op2.binary_debug.write('BEND-69 - eid=%s %s\n' % (eid, str(out)))
            #print('BEND-69 - eid=%s %s\n' % (eid, str(out)))

            (grid, angle, sc, sd, se, sf, omax, omin, mst, msc) = out

            add_sort_x(dt, eid, grid, angle, sc, sd, se, sf, omax, omin, mst, msc)
            n += ntotali
    return n


def oesrt_cquad4_95(op2: OP2, data: bytes, ndata: int) -> int:
    """unsupported element"""
    assert op2.num_wide == 9, f'num_wide={op2.num_wide} not 9'
    ntotal = 36 * op2.factor # 4*9
    #oesrt_cquad4_95

    n = 0
    struct1 = Struct(op2._endian + op2._analysis_code_fmt + b'8si3fi4s')
    nelements = ndata // ntotal
    #obj = op2.obj
    for unused_i in range(nelements):
        edata = data[n:n + ntotal]
        out = struct1.unpack(edata)  # num_wide=9
        if op2.is_debug_file:
            op2.binary_debug.write('CQUAD4-95 - %s\n' % str(out))
        #eid, failure, ply, failureIndexPly, failureIndexBonding, failureIndexMax, flag
        # 3,TSAIWU,1,8.5640,0.0,None

        (eid, failure, ply, strength_ratio_ply, failure_index_bonding, strength_ratio_bonding, flag, flag2) = out
        str((eid, failure, ply, strength_ratio_ply, failure_index_bonding, strength_ratio_bonding, flag, flag2))
        #strength_ratio_ply
        #print("eid=%s failure=%r ply=%s failureIndexPly=%s  failure_index_bonding=%s strength_ratio_bonding=%s flag=%s flag2=%s" % (
        #    eid, failure.strip(), ply, failureIndexPly, failure_index_bonding, strength_ratio_bonding, flag, flag2))
        #print("eid=%s strength_ratio_ply=%g failure_index_bonding=%s strength_ratio_bonding=%s" % (
            #eid, strength_ratio_ply, failure_index_bonding, strength_ratio_bonding))
        #obj.add_new_eid(element_name, dt, eid, force, stress)
        n += ntotal
    return n


def oes_csolid_composite_real(op2: OP2, data: bytes,
                              obj,
                              nelements: int, nedges: int,
                              element_name: str, preline1: str, preline2: str,
                              dt: Any) -> int:
    """
    1 PLY I Lamina number
    2 FLOC CHAR4 Fiber location (BOT, MID, TOP)
    3 GRID I Edge grid ID (center=0)

    4 EX1 RS Normal strain in the 1-direction
    5 EY1 RS Normal strain in the 2-direction
    6 EZ1 RS Normal strain in the 3-direction
    7 ET1 RS Shear strain in the 12-plane
    8 EL2 RS Shear strain in the 23-plane
    9 EL1 RS Shear strain in the 13-plane
    10 ETMAX1 RS von Mises strain
    For each fiber location requested (PLSLOC), words 3 through 10 repeat 4 times.
    """
    n = 0
    struct1 = Struct(op2._endian + op2._analysis_code_fmt + b'i4s')
    struct2 = Struct(op2._endian + b'i7f')
    if op2.is_debug_file:
        msg = '%s-%s nelements=%s nnodes=%s; C=[sxx, syy, szz, txy, tyz, txz, ovm,\n' % (
            op2.element_name, op2.element_type, nelements, nedges)
        op2.binary_debug.write(msg)

    ntotal1 = 12 * op2.factor
    ntotal2 = 32 * op2.factor
    for unused_i in range(nelements):
        edata = data[n:n+ntotal1]
        out = struct1.unpack(edata)
        (eid_device, ply, fiber_location) = out
        #print(out)
        eid, dt = get_eid_dt_from_eid_device(
            eid_device, op2.nonlinear_factor, op2.sort_method)

        if op2.is_debug_file:
            op2.binary_debug.write('%s - eid=%i; %s\n' % (preline1, eid, str(out)))

        n += ntotal1
        for inode in range(nedges):  # nodes pts, no centroid
            out = struct2.unpack(data[n:n + ntotal2])  # 4*8 = 32
            if op2.is_debug_file:
                op2.binary_debug.write('%s - %s\n' % (preline2, str(out)))
            (grid_device, sxx, syy, szz, txy, tyz, txz, ovm) = out

            if op2.is_debug_file:
                op2.binary_debug.write('  eid=%s inode=%i; C=[%s]\n' % (
                    eid, grid_device, ', '.join(['%r' % di for di in out])))

            #if grid_device == 0:
                #grid = 'CENTER'
            #else:
                ##grid = (grid_device - device_code) // 10
                #grid = grid_device

            grid = grid_device
            #a_cos = [a1, a2, a3]
            #b_cos = [b1, b2, b3]
            #c_cos = [c1, c2, c3]
            if 0:
                if inode == 0:
                    #  this is correct, but fails
                    #element_name = op2.element_name + str(nnodes)
                    obj.add_eid_sort1(element_name, cid, dt, eid, grid,
                                      sxx, syy, szz, txy, tyz, txz, ovm)
                else:
                    obj.add_node_sort1(dt, eid, inode, grid,
                                       sxx, syy, szz, txy, tyz, txz, ovm)
            n += ntotal2
    return n


def oes_shell_composite_complex_11(op2: OP2,
                                   data: bytes,
                                   obj: Union[ComplexLayeredCompositeStressArray, ComplexLayeredCompositeStrainArray],
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


def _oes_csolid2_real(op2: OP2, data: bytes,
                      n: int,
                      obj: Union[RealSolidStressArrayNx, RealSolidStrainArrayNx],
                      nnodes_expected: int,
                      nelements: int,
                      element_name: str,
                      stress_strain='stress') -> tuple[int, int, int]:

    obj = op2.obj
    preline1 = '%s-%s' % (op2.element_name, op2.element_type)
    preline2 = ' ' * len(preline1)
    #if is_vectorized and self.use_vector:  # pragma: no cover
        #self.log.debug('vectorize CSolid real SORT%s' % op2.sort_method)

    # 2 CID I Coordinate System
    # 3 CTYPE CHAR4 Grid or Gauss
    #
    # 4 GRID I Corner grid ID
    # 5 EX RS Strain in X
    # 6 EY RS Strain in Y
    # 7 EZ RS Strain in Z
    # 8 EXY RS Strain in XY
    # 9 EYZ RS Strain in YZ
    # 10 EZX RS Strain in ZX
    # 11 EVM RS Von Mises strain
    # Words 4 through 11 repeat nnodes times.
    if op2.size == 4:
        struct1 = Struct(op2._endian + op2._analysis_code_fmt + b'i4s')
        struct2 = Struct(op2._endian + b'i7f')
    else:
        struct1 = Struct(op2._endian + mapfmt(op2._analysis_code_fmt, op2.size) + b'q8s')
        struct2 = Struct(op2._endian + b'q7d')
    ntotal1 = 12 * op2.factor
    ntotal2 = 32 * op2.factor

    if op2.is_debug_file:
        msg = '%s-%s nelements=%s nnodes=%s; C=[sxx, syy, szz, txy, tyz, txz, ovm,\n' % (
            op2.element_name, op2.element_type, nelements, nnodes_expected)
        op2.binary_debug.write(msg)

    assert op2.is_sort1
    obj._times[obj.itime] = op2.nonlinear_factor
    for unused_i in range(nelements):
        edata = data[n:n+ntotal1]
        out = struct1.unpack(edata)
        (eid_device, cid, unused_abcd) = out
        eid, dt = get_eid_dt_from_eid_device(
            eid_device, op2.nonlinear_factor, op2.sort_method)
        #print(eid, dt, cid, unused_abcd)
        if op2.is_debug_file:
            op2.binary_debug.write('%s - eid=%i; %s\n' % (preline1, eid, str(out)))

        #assert nnodes < 21, 'print_block(data[n:n+16])'  #self.print_block(data[n:n+16])
        n += ntotal1
        for inode in range(nnodes_expected):  # nodes pts, no centroid
            out = struct2.unpack(data[n:n + ntotal2]) # 4*8 = 32
            if op2.is_debug_file:
                op2.binary_debug.write('%s - %s\n' % (preline2, str(out)))
            (grid_device, sxx, syy, szz, txy, tyz, txz, ovm) = out

            if op2.is_debug_file:
                op2.binary_debug.write('  eid=%s inode=%i; C=[%s]\n' % (
                    eid, grid_device, ', '.join(['%r' % di for di in out])))

            #grid = (grid_device - device_code) // 10
            grid = grid_device
            #print(eid, inode, grid)

            obj.add_node_sort1(dt, eid, inode, grid,
                               sxx, syy, szz, txy, tyz, txz, ovm)
            n += ntotal2
    return n


def oes_csolid_linear_hyperelastic_cosine_real(op2: OP2, data: bytes,
                                               nelements: int, nnodes_expected: int,
                                               preline1: str, preline2: str) -> int:
    n = 0
    # ELTYPE =140 Hyperelastic 8-noded hexahedron element linear format
    # (HEXAFD)
    # 2 TYPE CHAR4 Gaus
    #
    # 3 ID I
    # 4 SX RS
    # 5 SXY RS
    # 6 PA RS
    # 7 AX RS
    # 8 AY RS
    # 9 AZ RS
    # 10 PRESSURE RS
    # 11 SY RS
    # 12 SYZ RS
    # 13 PB RS
    # 14 BX RS
    # 15 BY RS
    # 16 BZ RS
    # 17 SZ RS
    # 18 SZX RS
    # 19 PC RS
    # 20 CX RS
    # 21 CY RS
    # 22 CZ RS
    # Words 3 through 22 repeat 008 times
    struct1 = Struct(op2._endian + op2._analysis_code_fmt + b'4s')
    struct2 = Struct(op2._endian + b'i19f')
    if op2.is_debug_file:
        msg = (
            f'{op2.element_name}-{op2.element_type} nelements={nelements} '
            f'nnodes={nnodes_expected}; '
            'C=[oxx, oxy, pa, ax, ay, az, pressure, '
            'oyy, oyz, pb, bx, by, bz, '
            'ozz, oxz, pc, cx, cy, cz]\n')
        op2.binary_debug.write(msg)

    for unused_i in range(nelements):
        edata = data[n:n+8]
        out = struct1.unpack(edata)
        (eid_device, grid_gauss, ) = out
        #print(out)
        eid, dt = get_eid_dt_from_eid_device(
            eid_device, op2.nonlinear_factor, op2.sort_method)

        if op2.is_debug_file:
            op2.binary_debug.write('%s - eid=%i; %s\n' % (preline1, eid, str(out)))
        #assert nnodes < 21, 'print_block(data[n:n+16])'  #self.print_block(data[n:n+16])

        n += 8
        for unused_inode in range(nnodes_expected):  # nodes pts, no centroid
            out = struct2.unpack(data[n:n + 80])  # 4*20 = 80
            if op2.is_debug_file:
                op2.binary_debug.write('%s - %s\n' % (preline2, str(out)))
            # nid, oxx, oxy, pa, ax, ay, az, pressure,
            #      oyy, oyz, pb, bx, by, bz,
            #      ozz, oxz, pc, cx, cy, cz
            (grid_device,
             oxx, oxy, pa, ax, ay, az, pressure,
             oyy, oyz, pb, bx, by, bz,
             ozz, oxz, pc, cx, cy, cz) = out
            #print(out)

            if op2.is_debug_file:
                op2.binary_debug.write('  eid=%s inode=%i; C=[%s]\n' % (
                    eid, grid_device, ', '.join(['%r' % di for di in out])))

            #if 0:
                #if inode == 0:
                    #  this is correct, but fails
                    #element_name = op2.element_name + str(nnodes)
                    #obj.add_eid_sort1(element_name, cid, dt, eid, grid,
                                      #sxx, syy, szz, txy, tyz, txz, ovm)
                #else:
                    #obj.add_node_sort1(dt, eid, inode, grid,
                                       #sxx, syy, szz, txy, tyz, txz, ovm)
            n += 80
    return n


def oes_csolid_linear_hyperelastic_real(op2: OP2, data: bytes, obj,
                                        nelements: int, nnodes_expected: int,
                                        preline1: str, preline2: str):
    n = 0
    #ELTYPE =163 Linear form for hyperelastic 20 node HEXAFD
    #2 TYPE CHAR4 Gaus
    #
    #3 ID I
    #4 SX RS
    #5 SXY RS
    #6 PA RS
    #7 AX RS
    #8 AY RS
    #9 AZ RS
    #10 PRESSURE RS
    #11 SY RS
    #12 SYZ RS
    #13 PB RS
    #14 BX RS
    #15 BY RS
    #16 BZ RS
    #17 SZ RS
    #18 SZX RS
    #19 PC RS
    #20 CX RS
    #21 CY RS
    #22 CZ RS
    #Words 3 through 22 repeat 027 times
    struct1 = Struct(op2._endian + op2._analysis_code_fmt + b'4s')
    struct2 = Struct(op2._endian + b'i19f')
    #if op2.is_debug_file:
        #msg = (
            #f'{op2.element_name}-{op2.element_type} nelements={nelements} '
            #f'nnodes={nnodes_expected}; '
            #'C=[sxx, syy, szz, txy, tyz, txz, pressure, '
            #'evol, exx, eyy, ezz, exy, eyz, exz]\n')
        #op2.binary_debug.write(msg)

    for unused_i in range(nelements):
        edata = data[n:n+8]
        out = struct1.unpack(edata)
        (eid_device, unused_abcd, ) = out
        eid, dt = get_eid_dt_from_eid_device(
            eid_device, op2.nonlinear_factor, op2.sort_method)

        if op2.is_debug_file:
            op2.binary_debug.write('%s - eid=%i; %s\n' % (preline1, eid, str(out)))
        #assert nnodes < 21, 'print_block(data[n:n+16])'  #self.print_block(data[n:n+16])

        n += 8
        for unused_inode in range(nnodes_expected):  # nodes pts, no centroid
            out = struct2.unpack(data[n:n + 80]) # 4*20 = 80
            if op2.is_debug_file:
                op2.binary_debug.write('%s - %s\n' % (preline2, str(out)))
            #(grid_device, sxx, syy, szz, txy, tyz, txz, pressure,
             #evol, exx, eyy, ezz, exy, eyz, exz) = out
            #print(out)

            if op2.is_debug_file:
                op2.binary_debug.write('  eid=%s inode=%i; C=[%s]\n' % (
                    eid, grid_device, ', '.join(['%r' % di for di in out])))

            #grid = grid_device
            #if 0:
                #if inode == 0:
                    ##  this is correct, but fails
                    ##element_name = op2.element_name + str(nnodes)
                    #obj.add_eid_sort1(element_name, cid, dt, eid, grid,
                                      #sxx, syy, szz, txy, tyz, txz, ovm)
                #else:
                    #obj.add_node_sort1(dt, eid, inode, grid,
                                       #sxx, syy, szz, txy, tyz, txz, ovm)
            n += 80
    return n
