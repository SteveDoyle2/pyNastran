from __future__ import annotations
from struct import Struct
from typing import TYPE_CHECKING
import numpy as np

from pyNastran.op2.op2_interface.op2_reader import mapfmt
from pyNastran.op2.tables.utils import get_is_slot_saved, get_eid_dt_from_eid_device
# from pyNastran.op2.op2_helper import polar_to_real_imag
# from pyNastran.op2.op2_interface.utils import apply_mag_phase

from pyNastran.op2.tables.oef_forces.oef_force_objects import (
    RealCGapForceArray,
)
if TYPE_CHECKING:  # pragma: no cover
    from pyNastran.op2.op2 import OP2


def oef_cgap(op2: OP2, data, ndata, dt, unused_is_magnitude_phase,
             result_type, prefix, postfix):
    """38-GAP"""
    result_name = prefix + 'cgap_force' + postfix
    if op2._results.is_not_saved(result_name):
        return ndata, None, None
    op2._results._found_result(result_name)
    slot = op2.get_result(result_name)

    factor = op2.factor
    if op2.format_code == 1 and op2.num_wide == 9:  # real
        ntotal = 36 * factor  # 9*4
        nelements = ndata // ntotal
        obj = op2.obj

        auto_return, is_vectorized = op2._create_oes_object4(
            nelements, result_name, slot, RealCGapForceArray)
        if auto_return:
            return nelements * ntotal, None, None

        obj = op2.obj
        if op2.use_vector and is_vectorized and op2.sort_method == 1:
            n = nelements * ntotal
            itotal = obj.ielement
            ielement2 = obj.itotal + nelements
            itotal2 = ielement2

            floats = np.frombuffer(data, dtype=op2.fdtype8).reshape(nelements, 9)
            obj._times[obj.itime] = dt
            if obj.itime == 0:
                ints = np.frombuffer(data, dtype=op2.idtype8).reshape(nelements, 9)
                eids = ints[:, 0] // 10
                assert eids.min() > 0, eids.min()
                obj.element[itotal:itotal2] = eids

            # [fx, sfy, sfz, u, v, w, sv, sw]
            obj.data[obj.itime, itotal:itotal2, :] = floats[:, 1:].copy()
            obj.itotal = itotal2
            obj.ielement = ielement2
        else:
            n = oef_cgap_real_9(op2, data, obj,
                                nelements, ntotal)

    else:  # pragma: no cover
        raise RuntimeError(op2.code_information())
        # msg = op2.code_information()
        # print(msg)
        # return op2._not_implemented_or_skip(data, ndata, msg), None, None
    return n, nelements, ntotal


def oef_cgap_real_9(op2, data: bytes,
                    obj: RealCGapForceArray,
                    nelements: int, ntotal: int) -> int:
    n = 0
    fmt = mapfmt(op2._analysis_code_fmt + b'8f', op2.size)
    s = Struct(op2._endian + fmt)
    add_sort_x = getattr(obj, 'add_sort' + str(op2.sort_method))
    for unused_i in range(nelements):
        edata = data[n:n+ntotal]

        out = s.unpack(edata)
        if op2.is_debug_file:
            op2.binary_debug.write('OEF_CGAP-38 - %s\n' % (str(out)))
        (eid_device, fx, sfy, sfz, u, v, w, sv, sw) = out
        eid, dt = get_eid_dt_from_eid_device(
            eid_device, op2.nonlinear_factor, op2.sort_method)
        #data_in = [eid, fx, sfy, sfz, u, v, w, sv, sw]
        #print "%s" %(self.get_element_type(op2.element_type)), data_in
        #eid = obj.add_new_eid_sort1(out)
        add_sort_x(dt, eid, fx, sfy, sfz, u, v, w, sv, sw)
        n += ntotal
    return n
