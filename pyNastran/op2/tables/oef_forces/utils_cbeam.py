from __future__ import annotations
from struct import Struct
from typing import Any, TYPE_CHECKING
import numpy as np

from pyNastran.op2.op2_interface.op2_reader import mapfmt
from pyNastran.op2.tables.utils import get_is_slot_saved, get_eid_dt_from_eid_device
from pyNastran.op2.op2_helper import polar_to_real_imag
from pyNastran.op2.op2_interface.utils import apply_mag_phase

from pyNastran.op2.tables.oef_forces.oef_force_objects import (
    RealCBeamForceArray,
)
from pyNastran.op2.tables.oef_forces.oef_complex_force_objects import (
    ComplexCBeamForceArray,
)
if TYPE_CHECKING:  # pragma: no cover
    from pyNastran.op2.op2 import OP2


def oef_cbeam(op2: OP2, data: bytes, ndata: int, dt, is_magnitude_phase: bool,
              result_type: str, prefix: str, postfix: str) -> int:
    """2-CBEAM"""
    n = 0
    result_name = prefix + 'cbeam_force' + postfix

    if op2._results.is_not_saved(result_name):
        return ndata, None, None
    op2._results._found_result(result_name)

    # if op2.format_code == 1 and op2.num_wide == 9:  # real centroid ???
    # raise RuntimeError('is this used?')
    # auto_return, is_vectorized = self._create_oes_object4(
    # nelements, result_name, slot, RealCBeamForceArray)
    # if auto_return:
    # return nelements * op2.num_wide * 4

    # obj = op2.obj
    ##is_vectorized = False
    # if op2.use_vector and is_vectorized and op2.sort_method == 1:
    # n = nelements * 4 * op2.num_wide
    # itotal = obj.itotal
    # itotal2 = obj.itotal + nelements
    # ielement2 = obj.ielement + nelements

    # floats = frombuffer(data, dtype=op2.fdtype8).reshape(nelements, 9)[:, 1:]
    # obj._times[obj.itime] = dt
    # if obj.itime == 0:
    # ints = frombuffer(data, dtype=op2.idtype8).reshape(nelements, 9)
    # eids = ints[:, 0] // 10
    # assert eids.min() > 0, eids.min()
    # assert 0 not in eids, eids

    # obj.element[itotal:itotal2] = eids
    # obj.element_node[itotal:itotal2, 0] = eids
    ##obj.element_node[itotal:itotal2, 1] = nids

    ##[sd, bm1, bm2, ts1, ts2, af, ttrq, wtrq]
    # obj.data[obj.itime, itotal:itotal2, :] = floats.copy()
    # obj.itotal = itotal2
    # obj.ielement = ielement2
    # else:
    # s = Struct(op2._endian + b'i8f')  # 36
    # ntotal = 36
    # nelements = ndata // ntotal
    # obj = op2.obj
    # for i in range(nelements):
    # edata = data[n:n+36]
    # out = s.unpack(edata)
    # if op2.is_debug_file:
    # op2.binary_debug.write('OEF_Beam - %s\n' % (str(out)))
    # (eid_device, sd, bm1, bm2, ts1, ts2, af, ttrq, wtrq) = out
    # eid, dt = get_eid_dt_from_eid_device(
    # eid_device, op2.nonlinear_factor, op2.sort_method)
    # n += 36

    factor = op2.factor
    if result_type in [0, 2] and op2.num_wide == 100:  # real/random
        if op2.sort_method == 2:
            msg = op2.code_information()
            if op2.read_mode == 2:
                return op2._not_implemented_or_skip(data, ndata, msg), None, None
            return ndata, None, None
        # real - format_code == 1
        # random - format_code == 2
        # result_name, is_random = self._apply_oef_ato_crm_psd_rms_no(result_name)
        slot = op2.get_result(result_name)
        ntotal = 400 * factor  # 1+(10-1)*11=100 ->100*4 = 400
        nelements = ndata // ntotal
        auto_return, is_vectorized = op2._create_oes_object4(
            nelements, result_name, slot, RealCBeamForceArray)
        if auto_return:
            op2._data_factor = 11
            return nelements * op2.num_wide * 4, None, None
        obj = op2.obj

        if op2.use_vector and is_vectorized and op2.sort_method == 1:
            n = nelements * ntotal
            itotal = obj.itotal
            itotal2 = obj.itotal + nelements * 11
            # ielement = obj.ielement
            ielement2 = obj.ielement + nelements

            floats = np.frombuffer(data, dtype=op2.fdtype8).reshape(nelements, 100)[:, 1:]
            obj._times[obj.itime] = dt
            if obj.itime == 0:
                ints = np.frombuffer(data, dtype=op2.idtype8).reshape(nelements, 100)
                eids = ints[:, 0] // 10
                assert eids.min() > 0, eids.min()
                assert 0 not in eids, eids
                eids2 = np.repeat(eids, 11)

                ints2 = ints[:, 1:].reshape(nelements * 11, 9)
                nids = ints2[:, 0]

                obj.element[itotal:itotal2] = eids2
                obj.element_node[itotal:itotal2, 0] = eids2
                obj.element_node[itotal:itotal2, 1] = nids

            # [nid, sd, bm1, bm2, ts1, ts2, af, ttrq, wtrq]
            floats2 = floats.reshape(nelements * 11, 9)[:, 1:].copy()
            # sd = floats2[:, 0]
            # obj.data[obj.itime, itotal:itotal2, :] = sd
            obj.data[obj.itime, itotal:itotal2, :] = floats2

            obj.itotal = itotal2
            obj.ielement = ielement2
        else:
            n = oef_cbeam_real_100(op2, data, obj,
                                   nelements, ntotal, dt)

    elif result_type == 1 and op2.num_wide == 177:  # imag
        slot = op2.get_result(result_name)
        ntotal = 708 * factor  # 177*4
        nelements = ndata // ntotal

        auto_return, is_vectorized = op2._create_oes_object4(
            nelements, result_name, slot, ComplexCBeamForceArray)
        if auto_return:
            op2._data_factor = 11
            return nelements * op2.num_wide * 4, None, None

        obj = op2.obj
        if op2.use_vector and is_vectorized and op2.sort_method == 1:
            n = nelements * ntotal
            itotal = obj.itotal
            itotal2 = obj.itotal + nelements * 11
            # ielement = obj.ielement
            ielement2 = obj.ielement + nelements

            floats = np.frombuffer(data, dtype=op2.fdtype8).reshape(nelements, 177)[:, 1:]
            obj._times[obj.itime] = dt
            if obj.itime == 0:
                ints = np.frombuffer(data, dtype=op2.idtype8).reshape(nelements, 177)
                eids = ints[:, 0] // 10
                assert eids.min() > 0, eids.min()
                assert 0 not in eids, eids
                eids2 = np.repeat(eids, 11)

                ints2 = ints[:, 1:].reshape(nelements * 11, 16)
                nids = ints2[:, 0].copy()

                obj.element[itotal:itotal2] = eids2
                obj.element_node[itotal:itotal2, 0] = eids2
                obj.element_node[itotal:itotal2, 1] = nids

            # [nid, sd, bm1r, bm2r, ts1r, ts2r, afr, ttrqr, wtrqr,
            #          bm1i, bm2i, ts1i, ts2i, afi, ttrqi, wtrqi]
            floats2 = floats.reshape(nelements * 11, 16)[:, 1:].copy()
            isave1 = slice(1, 8)
            isave2 = slice(8, None)
            real_imag = apply_mag_phase(floats2, is_magnitude_phase, isave1, isave2)

            sd = floats2[:, 0]
            obj.data[obj.itime, itotal:itotal2, 0] = sd
            obj.data[obj.itime, itotal:itotal2, 1:] = real_imag

            obj.itotal = itotal2
            obj.ielement = ielement2
        else:
            if op2.is_debug_file:
                op2.binary_debug.write('  [cap, element1, element2, ..., cap]\n')
                op2.binary_debug.write('  cap = %i  # assume 1 cap when there could have been multiple\n' % ndata)
                # op2.binary_debug.write('  #elementi = [eid_device, force]\n')
                # op2.binary_debug.write('  nelements=%i; nnodes=1 # centroid\n' % nelements)

            ntotal = 708 * factor  # (16*11+1)*4 = 177*4
            nelements = ndata // ntotal
            n = oef_cbeam_imag_177(op2, data, obj, nelements, ntotal, is_magnitude_phase)
    else:
        raise RuntimeError(op2.code_information())
        # msg = op2.code_information()
        # raise RuntimeError(msg)
        # print(msg)
        # return op2._not_implemented_or_skip(data, ndata, msg), None, None
    return n, nelements, ntotal


def oef_cbeam_real_100(op2: OP2, data: bytes, obj: RealCBeamForceArray,
                       nelements: int, ntotal: int, dt: Any) -> int:
    n = 0
    ntotal1 = 4 * op2.factor
    ntotal2 = 36 * op2.factor
    if op2.size == 4:
        s1 = op2.struct_i
        s2 = Struct(op2._endian + b'i8f')  # 36
    else:
        s1 = op2.struct_q
        s2 = Struct(op2._endian + b'q8d')  # 36

    sort_method = op2.sort_method
    add_sort_x = getattr(obj, 'add_sort' + str(sort_method))
    for unused_i in range(nelements):
        edata = data[n:n+ntotal1]
        eid_device, = s1.unpack(edata)
        eid, dt = get_eid_dt_from_eid_device(
            eid_device, op2.nonlinear_factor, sort_method)
        n += ntotal1

        for istation in range(11):
            edata = data[n:n+ntotal2]
            out = s2.unpack(edata)
            if op2.is_debug_file:
                op2.binary_debug.write('OEF_Beam - %s\n' % (str(out)))
            (nid, sd, bm1, bm2, ts1, ts2, af, ttrq, wtrq) = out

            if istation == 0 or sd > 0:
                add_sort_x(dt, eid, nid, sd, bm1, bm2, ts1, ts2,
                           af, ttrq, wtrq)
            n += ntotal2
    return n

def oef_cbeam_imag_177(op2: OP2, data: bytes,
                       obj: ComplexCBeamForceArray,
                       nelements: int, ntotal: int,
                       is_magnitude_phase: bool) -> int:
    n = 0
    #s1 = self.struct_i
    ntotal1 = 4 * op2.factor
    ntotal2 = 64 * op2.factor
    s1 = Struct(mapfmt(op2._endian + b'i', op2.size)) # self.struct_i
    s2 = Struct(mapfmt(op2._endian + b'i15f', op2.size))
    add_sort_x = getattr(obj, 'add_sort' + str(op2.sort_method))
    for unused_i in range(nelements):
        edata = data[n:n + ntotal1]
        eid_device, = s1.unpack(edata)
        eid, dt = get_eid_dt_from_eid_device(
            eid_device, op2.nonlinear_factor, op2.sort_method)

        n += ntotal1
        for unused_istation in range(11):
            edata = data[n:n+ntotal2]
            n += ntotal2
            out = s2.unpack(edata)
            if op2.is_debug_file:
                op2.binary_debug.write('OEF_Beam - %s\n' % (str(out)))
            (nid, sd,
             bm1r, bm2r, ts1r, ts2r, afr, ttrqr, wtrqr,
             bm1i, bm2i, ts1i, ts2i, afi, ttrqi, wtrqi) = out

            if is_magnitude_phase:
                bm1 = polar_to_real_imag(bm1r, bm1i)
                bm2 = polar_to_real_imag(bm2r, bm2i)
                ts1 = polar_to_real_imag(ts1r, ts1i)
                ts2 = polar_to_real_imag(ts2r, ts2i)
                af = polar_to_real_imag(afr, afi)
                ttrq = polar_to_real_imag(ttrqr, ttrqi)
                wtrq = polar_to_real_imag(wtrqr, wtrqi)
            else:
                bm1 = complex(bm1r, bm1i)
                bm2 = complex(bm2r, bm2i)
                ts1 = complex(ts1r, ts1i)
                ts2 = complex(ts2r, ts2i)
                af = complex(afr, afi)
                ttrq = complex(ttrqr, ttrqi)
                wtrq = complex(wtrqr, wtrqi)

            #if i == 0:
                #obj.add_new_element_sort1(
                    #dt, eid, nid, sd, bm1, bm2, ts1, ts2,
                    #af, ttrq, wtrq)
            #elif sd > 0.:
            add_sort_x(
                dt, eid, nid, sd, bm1, bm2, ts1, ts2,
                af, ttrq, wtrq)
            #else:
                ## don't add this field
                #pass
                #raise RuntimeError('CBEAM error; i=%s sd=%s' % (i, sd))
    return n
