from __future__ import annotations
from struct import Struct
from typing import Union, TYPE_CHECKING

from pyNastran.op2.tables.oes_stressStrain.real.oes_plate_strain_nx import (
    RealCPLSTRNPlateStressNXArray, RealCPLSTRNPlateStrainNXArray)
from pyNastran.op2.tables.utils import get_eid_dt_from_eid_device
if TYPE_CHECKING:
    from pyNastran.op2.op2 import OP2

def oes_cplstn3_real_6(op2: OP2, data: bytes,
                       obj: Union[RealCPLSTRNPlateStressNXArray, RealCPLSTRNPlateStrainNXArray],
                       nelements: int, ntotal: int, dt) -> int:
    n = 0
    nid = 0
    name = 'CPLSTN3'
    struct_cen = Struct(op2._endian + op2._analysis_code_fmt + b'5f')
    for unused_i in range(nelements):
        edata = data[n:n + ntotal]
        out = struct_cen.unpack(edata)
        n += ntotal
        (eid_device, oxx, oyy, ozz, txy, von_mises) = out
        eid, dt = get_eid_dt_from_eid_device(
            eid_device, op2.nonlinear_factor, op2.sort_method)
        obj.add_new_eid_sort1(dt, eid, nid, oxx, oyy, ozz, txy, von_mises)
        if op2.is_debug_file:
            op2.binary_debug.write('%s - %s\n' % (name, str(out)))
    return n

def oes_cplstn4_real_32(op2: OP2, data: bytes,
                        obj: Union[RealCPLSTRNPlateStressNXArray, RealCPLSTRNPlateStrainNXArray],
                        nelements: int, ntotal: int, dt) -> int:
    n = 0
    name = 'CPLSTN4'
    nnodes_cen = 5
    struct_cen = Struct(op2._endian + op2._analysis_code_fmt + b'4s i5f')
    struct_node = Struct(op2._endian + b'i5f')
    ntotal1 = (2 + 6) * op2.size
    ntotal2 = 6 * op2.size

    for unused_i in range(nelements):
        edata = data[n:n + ntotal1]
        out = struct_cen.unpack(edata)
        n += ntotal1

        (eid_device, cen, nid, oxx, oyy, ozz, txy, von_mises) = out
        assert cen == b'CEN ', cen
        eid, dt = get_eid_dt_from_eid_device(
            eid_device, op2.nonlinear_factor, op2.sort_method)
        if op2.is_debug_file:
            op2.binary_debug.write('%s - %s\n' % (name, str(out)))
        obj.add_new_eid_sort1(dt, eid, nid, oxx, oyy, ozz, txy, von_mises)
        for i in range(nnodes_cen - 1):
            edata = data[n:n + ntotal2]
            nid, oxx, oyy, ozz, txy, von_mises = struct_node.unpack(edata)
            obj.add_sort1(dt, eid, nid, oxx, oyy, ozz, txy, von_mises)
            n += ntotal2
    return n

def oes_cplstn6_real_26(op2: OP2, data: bytes,
                        obj: Union[RealCPLSTRNPlateStressNXArray, RealCPLSTRNPlateStrainNXArray],
                        nelements: int, ntotal: int, dt) -> int:
    n = 0
    nnodes_cen = 4
    name = 'CPLSTN6'
    struct_cen = Struct(op2._endian + op2._analysis_code_fmt + b'4s i5f')
    struct_node = Struct(op2._endian + b'i5f')
    ntotal1 = (2 + 6) * op2.size
    ntotal2 = 6 * op2.size
    for unused_i in range(nelements):
        edata = data[n:n + ntotal1]
        out = struct_cen.unpack(edata)
        n += ntotal1

        (eid_device, cen, nid, oxx, oyy, ozz, txy, von_mises) = out
        assert cen == b'CEN ', cen
        eid, dt = get_eid_dt_from_eid_device(
            eid_device, op2.nonlinear_factor, op2.sort_method)
        if op2.is_debug_file:
            op2.binary_debug.write('%s - %s\n' % (name, str(out)))
        obj.add_new_eid_sort1(dt, eid, nid, oxx, oyy, ozz, txy, von_mises)
        for i in range(nnodes_cen - 1):
            edata = data[n:n + ntotal2]
            out = struct_node.unpack(edata)
            #print('  ', out)
            nid, oxx, oyy, ozz, txy, von_mises = out
            if op2.is_debug_file:
                op2.binary_debug.write(f'  {out}\n')
            obj.add_sort1(dt, eid, nid, oxx, oyy, ozz, txy, von_mises)
            n += ntotal2
    return n
