from __future__ import annotations
from struct import Struct
from typing import TYPE_CHECKING
import numpy as np

from pyNastran.op2.op2_interface.op2_reader import mapfmt
from pyNastran.op2.op2_interface.utils import apply_mag_phase
from pyNastran.op2.tables.utils import get_is_slot_saved, get_eid_dt_from_eid_device
from pyNastran.op2.op2_helper import polar_to_real_imag

from pyNastran.op2.tables.oef_forces.oef_force_objects import (
    RealConeAxForceArray,
)

if TYPE_CHECKING:  # pragma: no cover
    from pyNastran.op2.op2 import OP2




def oef_cconeax(op2: OP2, data, ndata, dt, unused_is_magnitude_phase,
                result_type, prefix, postfix):
    """35-CONEAX"""
    result_name = prefix + 'cconeax_force' + postfix
    is_saved, slot = get_is_slot_saved(op2, result_name)
    if not is_saved:
        return ndata, None, None

    n = 0
    factor = op2.factor
    if op2.format_code == 1 and op2.num_wide == 7:  # real
        ntotal = 28 * factor # 7*4
        nelements = ndata // ntotal

        auto_return, is_vectorized = op2._create_oes_object4(
            nelements, result_name, slot, RealConeAxForceArray)
        if auto_return:
            return nelements * op2.num_wide * 4, None, None

        obj = op2.obj
        if op2.use_vector and is_vectorized and op2.sort_method == 1:
            n = nelements * 4 * op2.num_wide
            itotal = obj.ielement
            ielement2 = obj.itotal + nelements
            itotal2 = ielement2

            floats = np.frombuffer(data, dtype=op2.fdtype8).reshape(nelements, 7)
            obj._times[obj.itime] = dt
            if obj.itime == 0:
                ints = np.frombuffer(data, dtype=op2.idtype8).reshape(nelements, 7)
                eids = ints[:, 0] // 10
                assert eids.min() > 0, eids.min()
                obj.element[itotal:itotal2] = eids

            # [hopa, bmu, bmv, tm, su, sv]
            obj.data[obj.itime, itotal:itotal2, :] = floats[:, 1:].copy()
            obj.itotal = itotal2
            obj.ielement = ielement2
        else:
            n = oef_cconeax_real_7(op2, data, obj,
                                   nelements, ntotal, dt)
    else:  # pragma: no cover
        raise RuntimeError(op2.code_information())
        #msg = op2.code_information()
        #return op2._not_implemented_or_skip(data, ndata, msg), None, None
    return n, nelements, ntotal


def oef_cconeax_real_7(self, data: bytes,
                       obj: RealConeAxForceArray,
                       nelements: int, ntotal: int, dt: Any) -> int:
    op2 = self
    n = 0
    s = Struct(op2._endian + mapfmt(op2._analysis_code_fmt + b'6f', op2.size))
    add_sort_x = getattr(obj, 'add_sort' + str(op2.sort_method))
    for unused_i in range(nelements):
        edata = data[n:n+ntotal]
        out = s.unpack(edata)
        if op2.is_debug_file:
            op2.binary_debug.write('OEF_CONEAX-35 - %s\n' % (str(out)))
        (eid_device, hopa, bmu, bmv, tm, su, sv) = out
        eid, dt = get_eid_dt_from_eid_device(
            eid_device, op2.nonlinear_factor, op2.sort_method)
        add_sort_x(dt, eid, hopa, bmu, bmv, tm, su, sv)
        n += ntotal
    return n
