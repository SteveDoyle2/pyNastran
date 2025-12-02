from __future__ import annotations
from struct import Struct
from typing import TYPE_CHECKING
import numpy as np

from pyNastran.op2.op2_interface.op2_reader import mapfmt
from pyNastran.op2.op2_interface.utils import apply_mag_phase
from pyNastran.op2.tables.utils import get_is_slot_saved, get_eid_dt_from_eid_device
from pyNastran.op2.op2_helper import polar_to_real_imag

from pyNastran.op2.tables.oef_forces.oef_force_objects import (
    RealPlateForceArray,
)
from pyNastran.op2.tables.oef_forces.oef_complex_force_objects import (
    ComplexPlateForceArray,
)

if TYPE_CHECKING:  # pragma: no cover
    from pyNastran.op2.op2 import OP2

def oef_shells_centroidal(op2: OP2, data, ndata, dt, is_magnitude_phase,
                          result_type, prefix, postfix):
    """
    33-CQUAD4
    74-CTRIA3
    227-CTRIAR
    228-CQUADR

    """
    assert op2.element_name != 'RBAR', op2.code_information()
    n = 0
    if op2.element_type == 33:
        result_name = prefix + 'cquad4_force' + postfix
    elif op2.element_type in [74, 83]:
        result_name = prefix + 'ctria3_force' + postfix
    elif op2.element_type == 227:
        result_name = prefix + 'ctriar_force' + postfix
    elif op2.element_type == 228:
        result_name = prefix + 'cquadr_force' + postfix
    else:
        #msg = 'sort1 Type=%s num=%s' % (op2.element_name, op2.element_type)
        #return op2._not_implemented_or_skip(data, ndata, msg)
        raise NotImplementedError(op2.code_information())
    #result_name, is_random = self._apply_oef_ato_crm_psd_rms_no(result_name)

    is_saved, slot = get_is_slot_saved(op2, result_name)
    if not is_saved:
        return ndata, None, None

    assert op2._data_factor == 1, op2._data_factor

    factor = op2.factor
    if op2.format_code in [1, 2] and op2.num_wide == 9:
        # real - format_code == 1
        # random - format_code == 2
        ntotal = 36 * factor # 9*4
        nelements = ndata // ntotal
        auto_return, is_vectorized = op2._create_oes_object4(
            nelements, result_name, slot, RealPlateForceArray)
        if auto_return:
            return nelements * ntotal, None, None

        obj = op2.obj
        if op2.use_vector and is_vectorized and op2.sort_method == 1:
            n = nelements * ntotal
            ielement = obj.ielement
            ielement2 = ielement + nelements

            floats = np.frombuffer(data, dtype=op2.fdtype8).reshape(nelements, 9)
            obj._times[obj.itime] = dt
            if obj.itime == 0:
                ints = np.frombuffer(data, dtype=op2.idtype8).reshape(nelements, 9)
                eids = ints[:, 0] // 10
                assert eids.min() > 0, eids.min()
                obj.element[ielement:ielement2] = eids

            #[mx, my, mxy, bmx, bmy, bmxy, tx, ty]
            obj.data[obj.itime, ielement:ielement2, :] = floats[:, 1:].copy()
            obj.itotal = ielement2
            obj.ielement = ielement2
        else:
            n = oef_cquad4_33_real_9(op2, data, obj,
                                     nelements, ntotal)

    elif op2.format_code in [2, 3] and op2.num_wide == 17:  # imag
        ntotal = 68 * factor
        nelements = ndata // ntotal

        auto_return, is_vectorized = op2._create_oes_object4(
            nelements, result_name, slot, ComplexPlateForceArray)
        if auto_return:
            return nelements * ntotal, None, None

        obj = op2.obj
        if op2.use_vector and is_vectorized and op2.sort_method == 1:
            n = nelements * ntotal
            ielement = obj.ielement
            ielement2 = ielement + nelements
            itotal = obj.itotal
            itotal2 = itotal + nelements

            floats = np.frombuffer(data, dtype=op2.fdtype8).reshape(nelements, 17).copy()
            obj._times[obj.itime] = dt
            if obj.itime == 0:
                ints = np.frombuffer(data, dtype=op2.idtype8).reshape(nelements, 17)
                eids = ints[:, 0] // 10
                assert eids.min() > 0, eids.min()
                obj.element[itotal:itotal2] = eids

            #[mx, my, mxy, bmx, bmy, bmxy, tx, ty]
            isave1 = [1, 2, 3, 4, 5, 6, 7, 8]
            isave2 = [9, 10, 11, 12, 13, 14, 15, 16]
            real_imag = apply_mag_phase(floats, is_magnitude_phase, isave1, isave2)
            obj.data[obj.itime, itotal:itotal2, :] = real_imag
            obj.itotal = itotal2
            obj.ielement = ielement2
        else:
            n = oef_cquad4_33_imag_17(op2, data, ndata, obj,
                                      nelements, ntotal,
                                      is_magnitude_phase)

    else:  # pragma: no cover
        raise RuntimeError(op2.code_information())
        #msg = op2.code_information()
        #print(msg)
        #return op2._not_implemented_or_skip(data, ndata, msg), None, None
    return n, nelements, ntotal


def oef_cquad4_33_real_9(op2: OP2, data: bytes,
                         obj: RealPlateForceArray,
                         nelements: int, ntotal: int) -> int:
    n = 0
    s = Struct(mapfmt(op2._endian + op2._analysis_code_fmt + b'8f', op2.size))
    add_sort_x = getattr(obj, 'add_sort' + str(op2.sort_method))
    for unused_i in range(nelements):
        edata = data[n:n+ntotal]
        out = s.unpack(edata)
        if op2.is_debug_file:
            op2.binary_debug.write('real_OEF_Plate-%s - %s\n' % (op2.element_type, str(out)))
        (eid_device, mx, my, mxy, bmx, bmy, bmxy, tx, ty) = out
        eid, dt = get_eid_dt_from_eid_device(
            eid_device, op2.nonlinear_factor, op2.sort_method)
        add_sort_x(dt, eid, mx, my, mxy, bmx, bmy, bmxy, tx, ty)
        n += ntotal
    return n

def oef_cquad4_33_imag_17(op2: OP2, data: bytes, ndata: int,
                          obj: ComplexPlateForceArray,
                          nelements: int, ntotal: int,
                          is_magnitude_phase: bool) -> int:
    n = 0
    s = Struct(mapfmt(op2._endian + op2._analysis_code_fmt + b'16f', op2.size))
    add_sort_x = getattr(obj, 'add_sort' + str(op2.sort_method))
    for unused_i in range(nelements):
        edata = data[n:n+ntotal]
        out = s.unpack(edata)
        (eid_device,
         mxr, myr, mxyr, bmxr, bmyr, bmxyr, txr, tyr,
         mxi, myi, mxyi, bmxi, bmyi, bmxyi, txi, tyi) = out
        eid, dt = get_eid_dt_from_eid_device(
            eid_device, op2.nonlinear_factor, op2.sort_method)
        if op2.is_debug_file:
            op2.binary_debug.write('complex_OEF_Plate-%s - %s\n' % (op2.element_type, str(out)))

        if is_magnitude_phase:
            mx = polar_to_real_imag(mxr, mxi)
            my = polar_to_real_imag(myr, myi)
            mxy = polar_to_real_imag(mxyr, mxyi)
            bmx = polar_to_real_imag(bmxr, bmxi)
            bmy = polar_to_real_imag(bmyr, bmyi)
            bmxy = polar_to_real_imag(bmxyr, bmxyi)
            tx = polar_to_real_imag(txr, txi)
            ty = polar_to_real_imag(tyr, tyi)
        else:
            mx = complex(mxr, mxi)
            my = complex(myr, myi)
            mxy = complex(mxyr, mxyi)
            bmx = complex(bmxr, bmxi)
            bmy = complex(bmyr, bmyi)
            bmxy = complex(bmxyr, bmxyi)
            tx = complex(txr, txi)
            ty = complex(tyr, tyi)
        add_sort_x(dt, eid, mx, my, mxy, bmx, bmy, bmxy, tx, ty)
        n += ntotal
    return n
