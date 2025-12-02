from __future__ import annotations
from struct import Struct
from typing import TYPE_CHECKING
import numpy as np

from pyNastran.op2.op2_interface.op2_reader import mapfmt
from pyNastran.op2.tables.oes_stressStrain.real.oes_plate_strain_nx import (
    RealCPLSTRNPlateStressNXArray, RealCPLSTRNPlateStrainNXArray)
from pyNastran.op2.tables.utils import get_is_slot_saved, get_eid_dt_from_eid_device
if TYPE_CHECKING:
    from pyNastran.op2.op2 import OP2


def oes_cplstn_nx(op2: OP2, data, ndata: int, dt, is_magnitude_phase: bool,
                  result_type: str, prefix: str, postfix: str):
    """
    reads stress/strain for element type:
    - 271 CPLSTN3
    - 272 CPLSTN4
    - 273 CPLSTN6
    - 274 CPLSTN8

    """
    factor = op2.factor
    n = 0
    if op2.is_stress:
        stress_strain = 'stress'
        obj_vector_real = RealCPLSTRNPlateStressNXArray
    else:
        stress_strain = 'strain'
        obj_vector_real = RealCPLSTRNPlateStrainNXArray

    if op2.element_type == 271:  # CPLSTN3
        nnodes = 3
        nnodes_cen = 1
    elif op2.element_type == 272:  # CPLSTN4
        nnodes = 4
        nnodes_cen = 5
    elif op2.element_type == 273:  # CPLSTN6
        nnodes = 6
        nnodes_cen = 4
    elif op2.element_type == 274:  # CPLSTN8
        nnodes = 8
        nnodes_cen = 5
    else:  # pragma: no cover
        # raise RuntimeError(op2.element_type)
        raise RuntimeError(op2.code_information())
    result_name = f'{stress_strain}.cplstn{nnodes}_{stress_strain}'
    name = f'CPLSTN{nnodes}'

    is_saved, slot = get_is_slot_saved(op2, result_name)
    if not is_saved:
        return ndata, None, None

    num_wide = op2.num_wide
    ntotal = num_wide * op2.size
    num_wide_real = 32  # CPLSTN4, CPLSTN3
    num_wide_real = 26  # CPLSTN6
    nelements = ndata // ntotal
    assert ndata % ntotal == 0
    # print(f'result_name = {result_name}')

    if result_type == 0 and num_wide == 6 and nnodes == 3:  # real; CPLSTN3
        nlayers = nelements
        auto_return, is_vectorized = op2._create_oes_object4(
            nlayers, result_name, slot, obj_vector_real)
        if auto_return:
            assert ntotal == op2.num_wide * 4
            return nelements * ntotal, None, None

        obj: RealCPLSTRNPlateStressNXArray = op2.obj
        if op2.use_vector and is_vectorized and op2.sort_method == 1:
            n = nelements * op2.num_wide * 4
            istart = obj.itotal
            iend = istart + nlayers
            obj._times[obj.itime] = dt

            if obj.itime == 0:
                ints1 = np.frombuffer(data, dtype=op2.idtype).reshape(nelements, 6).copy()
                eids = ints1[:, 0] // 10
                obj.element[istart:iend] = eids
                obj.element_node[istart:iend, 0] = eids
                obj.element_node[istart:iend, 1] = 0
            floats = np.frombuffer(data, dtype=op2.fdtype).reshape(nelements, 6)[:, 1:]
            # [oxx, oyy, ozz, txy, von_mises]
            obj.data[obj.itime, istart:iend, :] = floats.copy()
        else:
            n = oes_cplstn3_real_6(
                op2, data, obj,
                nelements, ntotal, dt)

    elif result_type == 0 and num_wide == 32 and nnodes in {4, 8}:  # real; CPLSTN4
        nlayers = nelements * nnodes_cen
        auto_return, is_vectorized = op2._create_oes_object4(
            nlayers, result_name, slot, obj_vector_real)
        if auto_return:
            op2._data_factor = nnodes_cen
            assert ntotal == op2.num_wide * 4
            return nelements * ntotal, None, None

        obj: RealCPLSTRNPlateStressNXArray = op2.obj
        # if op2.is_debug_file:
        # op2.binary_debug.write('  [cap, element1, element2, ..., cap]\n')
        # op2.binary_debug.write('  cap = %i  # assume 1 cap when there could have been multiple\n' % ndata)
        # op2.binary_debug.write('  element1 = [eid_device, layer, o1, o2, t12, t1z, t2z, angle, major, minor, ovm)]\n')
        # op2.binary_debug.write('  nelements=%i; nnodes=1 # centroid\n' % nelements)

        if op2.use_vector and is_vectorized and op2.sort_method == 1:
            n = nelements * op2.num_wide * 4
            iestart = obj.ielement
            ieend = iestart + nelements
            istart = obj.itotal
            iend = istart + nlayers
            obj._times[obj.itime] = dt

            if obj.itime == 0:
                ints1 = np.frombuffer(data, dtype=op2.idtype).reshape(nelements, 32).copy()
                ints2 = ints1[:, 2:].reshape(nlayers, 6)
                eids = ints1[:, 0] // 10
                nids = ints2[:, 0]
                obj.element[iestart:ieend] = eids
                eids2 = np.vstack([eids] * nnodes_cen).T.ravel()
                obj.element_node[istart:iend, 0] = eids2
                obj.element_node[istart:iend, 1] = nids
            floats = np.frombuffer(data, dtype=op2.fdtype).reshape(nelements, 32)[:, 2:]
            floats2 = floats.reshape(nlayers, 6)
            # [oxx, oyy, ozz, txy, von_mises]
            obj.data[obj.itime, istart:iend, :] = floats2[:, 1:].copy()
        else:
            n = oes_cplstn4_real_32(
                op2, data, obj,
                nelements, ntotal, dt)

    elif result_type == 0 and num_wide == 26 and nnodes == 6:  # real; CPLSTN6
        nlayers = nelements * nnodes_cen
        auto_return, is_vectorized = op2._create_oes_object4(
            nlayers, result_name, slot, obj_vector_real)
        if auto_return:
            op2._data_factor = nnodes_cen
            assert ntotal == op2.num_wide * 4
            return nelements * ntotal, None, None

        obj: RealCPLSTRNPlateStressNXArray = op2.obj
        if op2.use_vector and is_vectorized and op2.sort_method == 1:
            n = nelements * op2.num_wide * 4
            iestart = obj.ielement
            ieend = iestart + nelements
            istart = obj.itotal
            iend = istart + nlayers
            obj._times[obj.itime] = dt

            if obj.itime == 0:
                ints1 = np.frombuffer(data, dtype=op2.idtype).reshape(nelements, 26).copy()
                ints2 = ints1[:, 2:].reshape(nlayers, 6)
                eids = ints1[:, 0] // 10
                nids = ints2[:, 0]
                obj.element[iestart:ieend] = eids
                eids2 = np.vstack([eids] * nnodes_cen).T.ravel()
                obj.element_node[istart:iend, 0] = eids2
                obj.element_node[istart:iend, 1] = nids
            floats = np.frombuffer(data, dtype=op2.fdtype).reshape(nelements, 26)[:, 2:]
            floats2 = floats.reshape(nlayers, 6)
            # [oxx, oyy, ozz, txy, von_mises]
            obj.data[obj.itime, istart:iend, :] = floats2[:, 1:].copy()
        else:
            n = oes_cplstn6_real_26(
                op2, data, obj,
                nelements, ntotal, dt)

    else:  # pragma: no cover
        raise RuntimeError(op2.code_information())
    return n, nelements, ntotal


def oes_cplstn3_real_6(op2: OP2, data: bytes,
                       obj: RealCPLSTRNPlateStressNXArray | RealCPLSTRNPlateStrainNXArray,
                       nelements: int, ntotal: int, dt) -> int:
    n = 0
    nid = 0
    name = 'CPLSTN3'
    fmt = mapfmt(op2._analysis_code_fmt + b'5f', op2.size)
    struct_cen = Struct(op2._endian + fmt)
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
                        obj: RealCPLSTRNPlateStressNXArray | RealCPLSTRNPlateStrainNXArray,
                        nelements: int, ntotal: int, dt) -> int:
    n = 0
    name = 'CPLSTN4'
    nnodes_cen = 5
    struct_cen = Struct(op2._endian + mapfmt(op2._analysis_code_fmt + b'4s i5f', op2.size))
    struct_node = Struct(op2._endian + mapfmt(b'i5f', op2.size))
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
                        obj: RealCPLSTRNPlateStressNXArray | RealCPLSTRNPlateStrainNXArray,
                        nelements: int, ntotal: int, dt) -> int:
    n = 0
    nnodes_cen = 4
    name = 'CPLSTN6'
    fmt = op2._analysis_code_fmt + b'4s i5f'
    struct_cen = Struct(op2._endian + fmt)
    struct_node = Struct(op2._endian + mapfmt(b'i5f', op2.fmt))
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
