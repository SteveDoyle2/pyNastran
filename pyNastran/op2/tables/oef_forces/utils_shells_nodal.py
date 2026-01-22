from __future__ import annotations
from struct import Struct
from typing import TYPE_CHECKING
import numpy as np

from pyNastran.op2.op2_interface.utils import mapfmt, real_imag_from_list, apply_mag_phase
from pyNastran.op2.tables.utils import get_is_slot_saved, get_eid_dt_from_eid_device
from pyNastran.op2.op2_helper import polar_to_real_imag

from pyNastran.op2.tables.oef_forces.oef_force_objects import (
    RealPlateBilinearForceArray,
)
from pyNastran.op2.tables.oef_forces.oef_complex_force_objects import (
    ComplexPlate2ForceArray,
)

if TYPE_CHECKING:  # pragma: no cover
    from pyNastran.op2.op2 import OP2

def oef_shells_nodal(op2: OP2, data, ndata, dt, is_magnitude_phase,
                     result_type, prefix, postfix):
    """
    64-CQUAD8
    70-CTRIAR
    75-CTRIA6
    82-CQUADR
    144-CQUAD4-bilinear

    """
    n = 0
    element_type = op2.element_type
    if element_type == 64:
        result_name = prefix + 'cquad8_force' + postfix
    elif element_type == 70:
        result_name = prefix + 'ctriar_force' + postfix
    elif element_type == 75:
        result_name = prefix + 'ctria6_force' + postfix
    elif element_type == 82:
        result_name = prefix + 'cquadr_force' + postfix
    elif element_type == 144:
        result_name = prefix + 'cquad4_force' + postfix
    else:
        raise NotImplementedError(op2.code_information())
        #msg = op2.code_information()
        #return op2._not_implemented_or_skip(data, ndata, msg)

    if element_type in [70, 75]:  # CTRIAR,CTRIA6
        nnodes = 3
    elif element_type in [64, 82, 144]:  # CQUAD8,CQUADR,CQUAD4-bilinear
        nnodes = 4
    else:
        raise NotImplementedError(op2.code_information())
        #msg = 'name=%r type=%r' % (op2.element_name, op2.element_type)
        #return op2._not_implemented_or_skip(data, ndata, msg), None, None

    is_saved, slot = get_is_slot_saved(op2, result_name)
    if not is_saved:
        return ndata, None, None

    nnodes_all = nnodes + 1
    numwide_real = 2 + nnodes_all * 9 # centroidal node is the + 1
    numwide_imag = 2 + nnodes_all * 17

    factor = op2.factor
    if op2.format_code == 1 and op2.num_wide == numwide_real:  # real
        obj_real = RealPlateBilinearForceArray

        ntotal = (8 + nnodes_all * 36) * factor # centroidal node is the + 1
        assert ntotal == op2.num_wide * 4 * factor, 'ntotal=%s numwide=%s' % (ntotal, op2.num_wide * 4)

        nelements = ndata // ntotal
        auto_return, is_vectorized = op2._create_oes_object4(
            nelements, result_name, slot, obj_real)
        if auto_return:
            op2._data_factor = nnodes_all
            return nelements * ntotal, None, None

        obj = op2.obj
        if op2.use_vector and is_vectorized and op2.sort_method == 1:
            nlayers = nelements * nnodes_all
            n = nelements * op2.num_wide * 4

            istart = obj.itotal
            iend = istart + nlayers
            obj._times[obj.itime] = dt

            if obj.itime == 0:
                ints = np.frombuffer(data, dtype=op2.idtype8).reshape(nelements, numwide_real).copy()
                # Nastran makes this a 4 for CQUAD4s instead
                # of 0 like the bilinear stress element...
                ints[:, 2] = 0

                nids = ints[:, 2:].reshape(nlayers, 9)[:, 0]
                eids = ints[:, 0] // 10
                eids2 = np.vstack([eids] * nnodes_all).T.ravel()
                obj.element_node[istart:iend, 0] = eids2
                obj.element_node[istart:iend, 1] = nids

            floats = np.frombuffer(data, dtype=op2.fdtype8).reshape(nelements, numwide_real)
            results = floats[:, 2:].reshape(nlayers, 9)[:, 1:].copy()
            #[mx, my, mxy, bmx, bmy, bmxy, tx, ty]
            obj.data[obj.itime, istart:iend, :] = results
        else:
            n = oef_cquad4_144_real_9(op2, data, obj,
                                      nelements, nnodes)

    elif op2.format_code in [2, 3] and op2.num_wide == numwide_imag: # complex
        ntotal = numwide_imag * 4 * factor
        nelements = ndata // ntotal

        auto_return, is_vectorized = op2._create_oes_object4(
            nelements, result_name, slot, ComplexPlate2ForceArray)
        if auto_return:
            op2._data_factor = nnodes_all
            return nelements * ntotal, None, None

        obj = op2.obj
        if op2.use_vector and is_vectorized and op2.sort_method == 1:
            n = nelements * ntotal
            itotal = obj.itotal
            ielement = obj.ielement
            ielement2 = obj.ielement + nelements
            itotal2 = obj.itotal + nelements * nnodes_all

            floats = np.frombuffer(data, dtype=op2.fdtype8).reshape(nelements, numwide_imag)
            obj._times[obj.itime] = dt
            if obj.itime == 0:
                ints = np.frombuffer(data, dtype=op2.idtype8).reshape(nelements, numwide_imag).copy()
                ints[:, 2] = 0
                ints2 = ints[:, 2:].reshape(nelements * nnodes_all, 17)

                eids = ints[:, 0] // 10
                nids = ints2[:, 0]
                assert eids.min() > 0, eids.min()
                eids2 = np.vstack([eids] * nnodes_all).T.ravel()
                obj.element[ielement:ielement2] = eids
                obj.element_node[itotal:itotal2, 0] = eids2
                obj.element_node[itotal:itotal2, 1] = nids

            #[mx, my, mxy, bmx, bmy, bmxy, tx, ty]
            floats2 = floats[:, 2:].reshape(nelements * nnodes_all, 17).copy()
            isave1 = [1, 2, 3, 4, 5, 6, 7, 8]
            isave2 = [9, 10, 11, 12, 13, 14, 15, 16]
            real_imag = apply_mag_phase(floats2, is_magnitude_phase, isave1, isave2)
            obj.data[obj.itime, itotal:itotal2, :] = real_imag
            obj.itotal = itotal2
            obj.ielement = ielement2
        else:
            n = oef_cquad4_imag_17(op2, data, ndata, obj,
                                   nelements, nnodes,
                                   is_magnitude_phase)

    else:  # pragma: no cover
        raise RuntimeError(op2.code_information())
        #msg = op2.code_information()
        #print(msg)
        #return op2._not_implemented_or_skip(data, ndata, msg), None, None
    return n, nelements, ntotal


def oef_cquad4_144_real_9(op2: OP2, data: bytes,
                          obj: RealPlateBilinearForceArray,
                          nelements: int, nnodes: int) -> int:
    n = 0
    n44 = 44 * op2.factor
    n36 = 36 * op2.factor
    if op2.size == 4:
        fmt1 = op2._endian + op2._analysis_code_fmt + b'4si8f'  # 8+36
        fmt2 = op2._endian + b'i8f' # 36
    else:
        fmt1 = op2._endian + mapfmt(op2._analysis_code_fmt, 8) + b'8sq8d'
        fmt2 = op2._endian + b'q8d'
    s1 = Struct(fmt1)
    s2 = Struct(fmt2)

    add_sort_x = getattr(obj, 'add_sort' + str(op2.sort_method))
    for unused_i in range(nelements):
        edata = data[n:n + n44]

        out = s1.unpack(edata)
        if op2.is_debug_file:
            op2.binary_debug.write('OEF_Plate2-%s - %s\n' % (op2.element_type, str(out)))
        (eid_device, term, _nid, mx, my, mxy, bmx, bmy, bmxy, tx, ty) = out
        #term= 'CEN\'
        #_nid = 4
        # -> CEN/4
        nid = 0
        inode = 0
        eid, dt = get_eid_dt_from_eid_device(
            eid_device, op2.nonlinear_factor, op2.sort_method)
        add_sort_x(dt, eid, term,
                   inode, nid, mx, my, mxy, bmx, bmy, bmxy, tx, ty)
        n += n44
        for jnode in range(nnodes):
            edata = data[n : n + n36]
            out = s2.unpack(edata)
            if op2.is_debug_file:
                op2.binary_debug.write('    %s\n' % (str(out)))
            (nid, mx, my, mxy, bmx, bmy, bmxy, tx, ty) = out
            assert nid > 0, 'nid=%s' % nid
            add_sort_x(dt, eid, term,
                       jnode+1, nid, mx, my, mxy, bmx, bmy, bmxy, tx, ty)
            n += n36
    return n

def oef_cquad4_imag_17(op2: OP2, data: bytes, ndata: int,
                       obj: ComplexPlate2ForceArray,
                       nelements: int, nnodes: int,
                       is_magnitude_phase: bool) -> int:
    n = 0
    factor = op2.factor
    size = op2.size
    if size == 4:
        s1 = Struct(op2._endian + op2._analysis_code_fmt + b'4si16f')  # 2+17=19 * 4 = 76
        s2 = Struct(op2._endian + b'i16f')  # 17 * 4 = 68
    else:
        s1 = Struct(op2._endian + mapfmt(op2._analysis_code_fmt, 8) + b'8sq16d')  # 2+17=19 * 4 = 768
        s2 = Struct(op2._endian + b'q16d')  # 17 * 4 = 68
    ntotal = (8 + (nnodes + 1) * 68) * factor
    ntotal1 = 76 * factor
    ntotal2 = 68 * factor

    nelements = ndata // ntotal
    obj = op2.obj
    add_new_element_sort_x = getattr(obj, 'add_new_element_sort' + str(op2.sort_method))
    add_sort_x = getattr(obj, 'add_sort' + str(op2.sort_method))
    for ielem in range(nelements):
        edata = data[n:n + ntotal1]
        #op2.show_data(edata)
        n += ntotal1

        out = s1.unpack(edata)
        if op2.is_debug_file:
            op2.binary_debug.write('OEF_Plate2-%s - %s\n' % (op2.element_type, str(out)))
        (eid_device, term, nid,
         mxr, myr, mxyr, bmxr, bmyr, bmxyr, txr, tyr,
         mxi, myi, mxyi, bmxi, bmyi, bmxyi, txi, tyi) = out
        #print('term =', term)
        #print('nid =', nid)
        #term = 'CEN\'

        eid, dt = get_eid_dt_from_eid_device(
            eid_device, op2.nonlinear_factor, op2.sort_method)
        mx, my, mxy, bmx, bmy, bmxy, tx, ty = real_imag_from_list([
            mxr, myr, mxyr, bmxr, bmyr, bmxyr, txr, tyr,
            mxi, myi, mxyi, bmxi, bmyi, bmxyi, txi, tyi], is_magnitude_phase)
        # this nid is just from CEN/4 <---- 4
        nid = 0
        add_new_element_sort_x(dt, eid, term, nid, mx, my, mxy,
                               bmx, bmy, bmxy, tx, ty)

        for inid in range(nnodes):  # .. todo:: fix crash...
            edata = data[n:n+ntotal2]
            n += ntotal2
            out = s2.unpack(edata)
            (nid,
             mxr, myr, mxyr, bmxr, bmyr, bmxyr, txr, tyr,
             mxi, myi, mxyi, bmxi, bmyi, bmxyi, txi, tyi) = out
            mx, my, mxy, bmx, bmy, bmxy, tx, ty = real_imag_from_list([
                mxr, myr, mxyr, bmxr, bmyr, bmxyr, txr, tyr,
                mxi, myi, mxyi, bmxi, bmyi, bmxyi, txi, tyi], is_magnitude_phase)
            if op2.is_debug_file:
                op2.binary_debug.write('OEF_Plate2 - eid=%i nid=%s out=%s\n' % (
                    eid, nid, str(out)))
            add_sort_x(dt, eid, inid, nid, mx, my, mxy, bmx, bmy, bmxy, tx, ty)
    #aaa
    return n
