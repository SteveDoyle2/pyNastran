from __future__ import annotations
from struct import Struct
from typing import Any, TYPE_CHECKING
import numpy as np

from pyNastran.op2.op2_interface.op2_reader import mapfmt
from pyNastran.op2.op2_helper import polar_to_real_imag

from pyNastran.op2.tables.utils import get_eid_dt_from_eid_device
from pyNastran.op2.tables.oes_stressStrain.real.oes_bars import RealBarStressArray, RealBarStrainArray
#from pyNastran.op2.tables.oes_stressStrain.real.oes_bars100 import RealBar10NodesStressArray, RealBar10NodesStrainArray

from pyNastran.op2.tables.oes_stressStrain.real.oes_bush import RealBushStressArray, RealBushStrainArray
from pyNastran.op2.tables.oes_stressStrain.real.oes_bush1d import RealBush1DStressArray
#from pyNastran.op2.tables.oes_stressStrain.real.oes_bush1d import RealBush1DStressArray
from pyNastran.op2.tables.oes_stressStrain.complex.oes_fast import ComplexFastStressArray, ComplexFastStrainArray
from pyNastran.op2.tables.oes_stressStrain.complex.oes_weld import ComplexWeldStressArray, ComplexWeldStrainArray
from pyNastran.op2.tables.oes_stressStrain.real.oes_composite_plates import RealCompositePlateStressArray, RealCompositePlateStrainArray
#RealCompositePlateStressStrengthRatioArray, RealCompositePlateStrainStrengthRatioArray = None, None
#RealCompositePlateStrainStrengthRatioArray = None
from pyNastran.op2.tables.oes_stressStrain.real.oes_composite_plates_strength_ratio import RealCompositePlateStressStrengthRatioArray  # , RealCompositePlateStrainStrengthRatioArray
#from pyNastran.op2.tables.oes_stressStrain.real.oes_gap import NonlinearGapStressArray
#from pyNastran.op2.tables.oes_stressStrain.real.oes_plate_strain import RealCPLSTRNPlateStressArray, RealCPLSTRNPlateStrainArray
from pyNastran.op2.tables.oes_stressStrain.real.oes_shear import RealShearStrainArray, RealShearStressArray
from pyNastran.op2.tables.oes_stressStrain.real.oes_solids import RealSolidStrainArray, RealSolidStressArray
from pyNastran.op2.tables.oes_stressStrain.real.oes_solids_nx import RealSolidStressArrayNx, RealSolidStrainArrayNx
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
from pyNastran.op2.tables.oes_stressStrain.complex.oes_shear import ComplexShearStressArray, ComplexShearStrainArray
from pyNastran.op2.tables.oes_stressStrain.complex.oes_solids import ComplexSolidStressArray, ComplexSolidStrainArray
from pyNastran.op2.tables.oes_stressStrain.complex.oes_springs import ComplexSpringStressArray, ComplexSpringStrainArray
#from pyNastran.op2.tables.oes_stressStrain.complex.oes_bend import ComplexBendStressArray, ComplexBendStrainArray

from pyNastran.op2.tables.oes_stressStrain.random.oes_bars import RandomBarStressArray, RandomBarStrainArray
#from pyNastran.op2.tables.oes_stressStrain.random.oes_bend import RandomBendStressArray, RandomBendStrainArray
from pyNastran.op2.tables.oes_stressStrain.random.oes_plates import RandomPlateStressArray, RandomPlateStrainArray
from pyNastran.op2.tables.oes_stressStrain.random.oes_plates_vm import RandomPlateVMStressArray, RandomPlateVMStrainArray
from pyNastran.op2.tables.oes_stressStrain.random.oes_solids import RandomSolidStressArray, RandomSolidStrainArray
from pyNastran.op2.tables.oes_stressStrain.random.oes_shear import RandomShearStressArray, RandomShearStrainArray
#from pyNastran.op2.tables.oes_stressStrain.random.oes_composite_plates import RandomCompositePlateStressArray, RandomCompositePlateStrainArray

#from pyNastran.op2.tables.oes_stressStrain.oes_nonlinear_rod import RealNonlinearRodArray
#from pyNastran.op2.tables.oes_stressStrain.oes_nonlinear_bush import RealNonlinearBushArray
# from pyNastran.op2.tables.oes_stressStrain.oes_hyperelastic import (
#     HyperelasticQuadArray)
#from pyNastran.op2.tables.oes_stressStrain.oes_nonlinear import RealNonlinearPlateArray, RealNonlinearSolidArray
if TYPE_CHECKING:  # pragma: no cover
    from pyNastran.op2.op2 import OP2


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


def oes_weldp_msc_real_8(op2: OP2, data: bytes,
                         obj: int | float,
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
                             obj: ComplexWeldStressArray | ComplexWeldStrainArray,
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
                         obj: int | float,
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
                             obj: ComplexFastStressArray | ComplexFastStrainArray,
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
                      obj: RealShearStressArray | RealShearStrainArray,
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
                         obj: ComplexShearStressArray | ComplexShearStrainArray,
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


def oes_cbush_real_7(op2: OP2, data: bytes,
                     obj: RealBushStressArray | RealBushStrainArray,
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
                         obj: ComplexCBushStressArray | ComplexCBushStrainArray,
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
                        obj: RealTriaxStressArray | RealTriaxStrainArray,
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
                    obj: RealSolidStressArray | RealSolidStrainArray,
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
                       obj: ComplexSolidStressArray | ComplexSolidStrainArray,
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
                      obj: RandomSolidStressArray | RandomSolidStrainArray,
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
                        obj: RandomShearStressArray | RandomShearStrainArray,
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
                      obj: RealBendStressArray | RealBendStrainArray,
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


def _oes_csolid2_real(op2: OP2, data: bytes,
                      n: int,
                      obj: RealSolidStressArrayNx | RealSolidStrainArrayNx,
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
                grid_device = out[0]
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


def oes_cbush1d_real_8(op2: OP2, data: bytes,
                       obj: RealBush1DStressArray,
                       nelements: int, ntotal: int,
                       is_magnitude_phase: bool) -> int:
    struct1 = Struct(op2._endian + op2._analysis_code_fmt + b'6fi')
    for unused_i in range(nelements):
        edata = data[n:n + ntotal]
        out = struct1.unpack(edata)  # num_wide=25
        if op2.is_debug_file:
            op2.binary_debug.write('CBUSH1D-40 - %s\n' % (str(out)))
        (eid_device, fe, ue, ve, ao, ae, ep, fail) = out
        eid, dt = get_eid_dt_from_eid_device(
            eid_device, op2.nonlinear_factor, op2.sort_method)

        # axial_force, axial_displacement, axial_velocity, axial_stress,
        # axial_strain, plastic_strain, is_failed
        obj.add_sort1(dt, eid, fe, ue, ve, ao, ae, ep, fail)
        n += ntotal
    return n
