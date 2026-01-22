from __future__ import annotations
from struct import Struct
from typing import Any, TYPE_CHECKING
import numpy as np

from pyNastran.op2.op2_interface.utils import mapfmt, real_imag_from_list
from pyNastran.op2.op2_helper import polar_to_real_imag

from pyNastran.op2.tables.utils import get_is_slot_saved, get_eid_dt_from_eid_device

from pyNastran.op2.tables.oes_stressStrain.complex.oes_fast import ComplexFastStressArray, ComplexFastStrainArray
from pyNastran.op2.tables.oes_stressStrain.complex.oes_weld import ComplexWeldStressArray, ComplexWeldStrainArray

#from pyNastran.op2.tables.oes_stressStrain.real.oes_solids import RealSolidStrainArray, RealSolidStressArray
from pyNastran.op2.tables.oes_stressStrain.real.oes_solids_nx import RealSolidStressArrayNx, RealSolidStrainArrayNx

from pyNastran.op2.tables.oes_stressStrain.complex.oes_solids import ComplexSolidStressArray, ComplexSolidStrainArray
from pyNastran.op2.tables.oes_stressStrain.random.oes_solids import RandomSolidStressArray, RandomSolidStrainArray

#from pyNastran.op2.tables.oes_stressStrain.oes_nonlinear_rod import RealNonlinearRodArray
#from pyNastran.op2.tables.oes_stressStrain.oes_nonlinear_bush import RealNonlinearBushArray
# from pyNastran.op2.tables.oes_stressStrain.oes_hyperelastic import (
#     HyperelasticQuadArray)
#from pyNastran.op2.tables.oes_stressStrain.oes_nonlinear import RealNonlinearPlateArray, RealNonlinearSolidArray
if TYPE_CHECKING:  # pragma: no cover
    from pyNastran.op2.op2 import OP2


def obj_set_element(op2: OP2, obj,
                    ielement: int, ielement2: int,
                    data: bytes, nelements: int):
    if obj.itime == 0:
        ints = np.frombuffer(data, dtype=op2.idtype8).reshape(nelements, op2.num_wide).copy()
        eids = ints[:, 0] // 10
        assert eids.min() > 0, eids.min()
        obj.element[ielement:ielement2] = eids

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

        axial, max_a, min_a, max_b, min_b, max_shear, bearing = real_imag_from_list([
            raxial, rmax_a, rmin_a, rmax_b, rmin_b, rmax_shear, rbearing,
            iaxial, imax_a, imin_a, imax_b, imin_b, imax_shear, ibearing], is_magnitude_phase)

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

        force_x, force_y, force_z, moment_x, moment_y, moment_z = real_imag_from_list([
            rforce_x, rforce_y, rforce_z, rmoment_x, rmoment_y, rmoment_z,
            iforce_x, iforce_y, iforce_z, imoment_x, imoment_y, imoment_z], is_magnitude_phase)
        add_sort_x(dt, eid, force_x, force_y, force_z, moment_x, moment_y, moment_z)
        n += ntotal
    return n


def oesrt_cquad4_95(op2: OP2, data: bytes, ndata: int) -> int:
    """unsupported element"""
    assert op2.num_wide == 9, f'num_wide={op2.num_wide} not 9'
    ntotal = 36 * op2.factor # 4*9
    #oesrt_cquad4_95

    n = 0
    if op2.size == 4:
        struct1 = Struct(op2._endian + op2._analysis_code_fmt + b'8si3fi4s')
    else:
        struct1 = Struct(op2._endian + mapfmt(op2._analysis_code_fmt, 8) + b'16sq3dq8s')
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
    if op2.size == 4:
        struct1 = Struct(op2._endian + op2._analysis_code_fmt + b'i4s')
        struct2 = Struct(op2._endian + b'i7f')
    else:
        struct1 = Struct(op2._endian + mapfmt(op2._analysis_code_fmt, 8) + b'q8s')
        struct2 = Struct(op2._endian + b'q7d')
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
    if op2.size == 4:
        struct1 = Struct(op2._endian + op2._analysis_code_fmt + b'4s')
        struct2 = Struct(op2._endian + b'i19f')
    else:
        struct1 = Struct(op2._endian + mapfmt(op2._analysis_code_fmt, 8) + b'8s')
        struct2 = Struct(op2._endian + b'q19d')
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
    if op2.size == 4:
        struct1 = Struct(op2._endian + op2._analysis_code_fmt + b'4s')
        struct2 = Struct(op2._endian + b'i19f')
    else:
        struct1 = Struct(op2._endian + mapfmt(op2._analysis_code_fmt, 8) + b'8s')
        struct2 = Struct(op2._endian + b'q19d')
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
