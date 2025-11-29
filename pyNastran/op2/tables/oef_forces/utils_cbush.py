from __future__ import annotations
from struct import Struct
from typing import Any, TYPE_CHECKING
import numpy as np

from pyNastran.op2.op2_interface.op2_reader import mapfmt
from pyNastran.op2.tables.utils import get_eid_dt_from_eid_device
from pyNastran.op2.op2_helper import polar_to_real_imag
# from pyNastran.op2.op2_interface.utils import apply_mag_phase

from pyNastran.op2.tables.oef_forces.oef_force_objects import (
    RealCBushForceArray,
    RealCFastForceArrayMSC,
    RealCBearForceArray,
)
from pyNastran.op2.tables.oef_forces.oef_complex_force_objects import (
    ComplexCBushForceArray,
    ComplexCFastForceArrayMSC,
    ComplexCBearForceArray,
)
if TYPE_CHECKING:  # pragma: no cover
    from pyNastran.op2.op2 import OP2


def oef_cbush(op2: OP2, data, ndata, dt, is_magnitude_phase,
              result_type, prefix, postfix):
    """
    102-CBUSH
    126-CFAST-MSC
    280-CBEAR

    """
    num_wide = op2.num_wide
    factor = op2.factor
    if op2.element_type == 102:
        result_name = prefix + 'cbush_force' + postfix
        real_obj = RealCBushForceArray
        complex_obj = ComplexCBushForceArray
        assert num_wide in [7, 13], op2.code_information()
    elif op2.element_type == 126:
        result_name = prefix + 'cfast_force' + postfix
        real_obj = RealCFastForceArrayMSC
        complex_obj = ComplexCFastForceArrayMSC
        assert num_wide in [7, 13], op2.code_information()
    elif op2.element_type == 280:
        result_name = prefix + 'cbear_force' + postfix
        assert num_wide in [7, 13], op2.code_information()
        real_obj = RealCBearForceArray
        complex_obj = ComplexCBearForceArray
    else:
        raise NotImplementedError(op2.code_information())

    # print(op2.code_information())
    if op2._results.is_not_saved(result_name):
        return ndata, None, None
    # result_name, is_random = self._apply_oef_ato_crm_psd_rms_no(result_name)
    op2._results._found_result(result_name)
    slot = op2.get_result(result_name)

    n = 0
    # op2.log.warning('dt=%s num_wide=%s result_type=%s', dt, num_wide, result_type)
    if result_type in {0, 2} and num_wide == 7:  # real/random
        numwide_real = 7
        # real - format_code == 1
        # random - format_code == 3
        ntotal = 28 * factor  # 7*4
        nelements = ndata // ntotal

        auto_return, is_vectorized = op2._create_oes_object4(
            nelements, result_name, slot, real_obj)
        if auto_return:
            return nelements * ntotal, None, None

        obj = op2.obj
        if op2.use_vector and is_vectorized and op2.sort_method == 1:
            # self.itime = 0
            # self.ielement = 0
            # self.itotal = 0
            # self.ntimes = 0
            # self.nelements = 0
            n = nelements * ntotal

            istart = obj.itotal
            iend = istart + nelements
            obj._times[obj.itime] = dt

            if obj.itime == 0:
                ints = np.frombuffer(data, dtype=op2.idtype8).reshape(nelements, numwide_real).copy()
                eids = ints[:, 0] // 10
                obj.element[istart:iend] = eids
            results = np.frombuffer(data, dtype=op2.fdtype8).reshape(nelements, numwide_real)

            # [fx, fy, fz, mx, my, mz]
            obj.data[obj.itime, istart:iend, :] = results[:, 1:].copy()
        else:
            n = oef_cbush_real_7(op2, data, obj,
                                 nelements, ntotal, dt)
    elif result_type == 1 and op2.num_wide == 13:  # imag
        # TCODE,7 =1 Real/imaginary or magnitude/phase
        # 2 FXR RS Force x - real/mag. part
        # 3 FYR RS Force y - real/mag. part
        # 4 FZR RS Force z - real/mag. part
        # 5 MXR RS Moment x - real/mag. part
        # 6 MYR RS Moment y - real/mag. part
        # 7 MZR RS Moment z - real/mag. part
        # 8 FXI RS Force x - imag./phase part
        # 9 FYI RS Force y - imag./phase part
        # 10 FZI RS Force z - imag./phase part
        # 11 MXI RS Moment x - imag./phase part
        # 12 MYI RS Moment y - imag./phase part
        # 13 MZI RS Moment z - imag./phase part

        # TODO: vectorize
        ntotal = 52 * factor  # 13*4
        nelements = ndata // ntotal
        # result_name = prefix + 'cbush_force' + postfix
        auto_return, is_vectorized = op2._create_oes_object4(
            nelements, result_name, slot, complex_obj)
        if auto_return:
            return nelements * ntotal, None, None

        obj = op2.obj
        n = oef_cbush_imag_13(op2, data, obj,
                              nelements, ntotal,
                              is_magnitude_phase)
    # elif op2.format_code == 2 and op2.num_wide == 7:
    # op2.log.warning(op2.code_information())
    else:
        raise NotImplementedError(op2.code_information())
    # else:  # pragma: no cover
    # msg = op2.code_information()
    # print(msg)
    # return op2._not_implemented_or_skip(data, ndata, msg), None, None
    return n, nelements, ntotal


def oef_cbush_real_7(self, data: bytes,
                     obj: RealCBushForceArray | RealCFastForceArrayMSC | RealCBearForceArray,
                     nelements: int, ntotal: int, dt: Any) -> int:
    op2 = self
    n = 0
    struct1 = Struct(op2._endian + mapfmt(op2._analysis_code_fmt + b'6f', op2.size))
    add_sort_x = getattr(obj, 'add_sort' + str(op2.sort_method))
    for unused_i in range(nelements):
        edata = data[n:n+ntotal]
        out = struct1.unpack(edata)
        if op2.is_debug_file:
            op2.binary_debug.write('OEF_CBUSH-102 - %s\n' % (str(out)))
        (eid_device, fx, fy, fz, mx, my, mz) = out
        eid, dt = get_eid_dt_from_eid_device(
            eid_device, op2.nonlinear_factor, op2.sort_method)
        #op2.log.debug('eid=%s dt=%s', eid, dt)
        add_sort_x(dt, eid, fx, fy, fz, mx, my, mz)
        n += ntotal
    return n

def oef_cbush_imag_13(self, data: bytes,
                      obj: ComplexCBushForceArray | ComplexCFastForceArrayMSC,
                      nelements: int, ntotal: int,
                      is_magnitude_phase: bool) -> int:
    """
    102-CBUSH
    126-CFAST-MSC
    280-CBEAR

    """
    op2 = self
    n = 0
    s = Struct(op2._endian + mapfmt(op2._analysis_code_fmt + b'12f', op2.size))
    add_sort_x = getattr(obj, 'add_sort' + str(op2.sort_method))
    for unused_i in range(nelements):
        edata = data[n:n + ntotal]
        out = s.unpack(edata)
        if op2.is_debug_file:
            op2.binary_debug.write('OEF_CBUSH-102 - %s\n' % (str(out)))
        (eid_device,
         fxr, fyr, fzr, mxr, myr, mzr,
         fxi, fyi, fzi, mxi, myi, mzi) = out
        eid, dt = get_eid_dt_from_eid_device(
            eid_device, op2.nonlinear_factor, op2.sort_method)

        if is_magnitude_phase:
            fx = polar_to_real_imag(fxr, fxi)
            mx = polar_to_real_imag(mxr, mxi)
            fy = polar_to_real_imag(fyr, fyi)
            my = polar_to_real_imag(myr, myi)
            fz = polar_to_real_imag(fzr, fzi)
            mz = polar_to_real_imag(mzr, mzi)
        else:
            fx = complex(fxr, fxi)
            mx = complex(mxr, mxi)
            fy = complex(fyr, fyi)
            my = complex(myr, myi)
            fz = complex(fzr, fzi)
            mz = complex(mzr, mzi)

        add_sort_x(dt, eid, fx, fy, fz, mx, my, mz)
        n += ntotal
    return n
