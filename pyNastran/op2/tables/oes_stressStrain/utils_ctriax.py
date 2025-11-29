from __future__ import annotations
from struct import Struct
from typing import Any, TYPE_CHECKING
import numpy as np

from pyNastran.op2.op2_interface.op2_reader import mapfmt
from pyNastran.op2.op2_helper import polar_to_real_imag
from pyNastran.op2.tables.utils import get_eid_dt_from_eid_device
from pyNastran.op2.op2_interface.utils import (
    apply_mag_phase,
)
from pyNastran.op2.tables.oes_stressStrain.real.oes_triax import RealTriaxStressArray, RealTriaxStrainArray
from pyNastran.op2.tables.oes_stressStrain.complex.oes_triax import ComplexTriaxStressArray, ComplexTriaxStrainArray

if TYPE_CHECKING:  # pragma: no cover
    from pyNastran.op2.op2 import OP2


def oes_ctriax6(op2: OP2, data, ndata, dt, is_magnitude_phase,
                result_type, prefix, postfix):
    """
    reads stress/strain for element type:
     - 53 : CTRIAX6
    """
    factor = op2.factor
    # n = 0
    if op2.is_stress:
        result_name = f'{prefix}ctriax_stress{postfix}'
    else:
        result_name = f'{prefix}ctriax_strain{postfix}'

    if op2._results.is_not_saved(result_name):
        return ndata, None, None
    op2._results._found_result(result_name)

    slot = op2.get_result(result_name)
    if result_type == 0 and op2.num_wide == 33:  # real
        obj_vector_real = RealTriaxStressArray if op2.is_stress else RealTriaxStrainArray

        ntotal = 132 * factor  # (1+8*4)*4 = 33*4 = 132
        nelements = ndata // ntotal
        auto_return, is_vectorized = op2._create_oes_object4(
            nelements, result_name, slot, obj_vector_real)
        if auto_return:
            op2._data_factor = 4
            assert ntotal == op2.num_wide * 4
            return nelements * ntotal, None, None

        obj = op2.obj
        nnodes_all = 4
        if op2.use_vector and is_vectorized and op2.sort_method == 1:
            n = nelements * op2.num_wide * 4

            itotal = obj.itotal
            itotal2 = itotal + nelements * nnodes_all
            ielement = obj.ielement
            ielement2 = ielement + nelements

            floats = np.frombuffer(data, dtype=op2.fdtype).reshape(nelements, 33)
            floats1 = floats[:, 1:].reshape(nelements * nnodes_all, 8).copy()

            obj._times[obj.itime] = dt
            if obj.itime == 0:
                ints = np.frombuffer(data, dtype=op2.idtype).reshape(nelements, 33).copy()
                ints1 = ints[:, 1:].reshape(nelements * nnodes_all, 8)
                eids = ints[:, 0] // 10
                ints[:, 0] = 0
                nids = ints1[:, 0]
                eids2 = np.vstack([eids] * nnodes_all).T.ravel()
                assert eids.min() > 0, eids.min()
                obj.element_node[itotal:itotal2, 0] = eids2
                obj.element_node[itotal:itotal2, 1] = nids

            # [loc, rs, azs, As, ss, maxp, tmax, octs]
            obj.data[obj.itime, itotal:itotal2, :] = floats1[:, 1:]
            obj.ielement = ielement2
            obj.itotal = itotal2
        else:
            n = oes_ctriax6_real_33(op2, data, obj,
                                    nelements, ntotal, dt)

    elif result_type == 1 and op2.num_wide == 37:  # imag
        # TODO: vectorize object
        if op2.is_stress:
            # print('op2.element_type', op2.element_type)
            # print('op2.element_name', op2.element_name)
            obj_vector_complex = ComplexTriaxStressArray
        else:
            obj_vector_complex = ComplexTriaxStrainArray

        num_wide = 1 + 4 * 9
        ntotal = num_wide * 4 * factor
        assert num_wide == op2.num_wide, num_wide
        nelements = ndata // ntotal  # (1+8*4)*4 = 33*4 = 132
        leftover = ndata % ntotal
        assert leftover == 0, 'ntotal=%s nelements=%s leftover=%s' % (ntotal, nelements, leftover)

        auto_return, is_vectorized = op2._create_oes_object4(
            nelements, result_name, slot, obj_vector_complex)

        if data is None:
            return ndata, None, None
        auto_return = False
        is_vectorized = False
        if auto_return:
            op2._data_factor = 4
            return nelements * ntotal, None, None

        obj = op2.obj
        nnodes_all = 4
        if op2.use_vector and is_vectorized and 0:
            n = nelements * ntotal
            itotal = obj.itotal
            itotal2 = itotal + nelements * nnodes_all
            ielement = obj.ielement
            ielement2 = ielement + nelements

            numwide_imag = 37
            floats = np.frombuffer(data, dtype=op2.fdtype).reshape(nelements, numwide_imag)
            floats1 = floats[:, 1:].reshape(nelements * nnodes_all, 9).copy()

            obj._times[obj.itime] = dt
            if obj.itime == 0:
                ints = np.frombuffer(data, dtype=op2.idtype).reshape(nelements, numwide_imag)
                ints1 = ints[:, 1:].reshape(nelements * nnodes_all, 9).copy()
                eids = ints[:, 0] // 10
                ints[:, 0] = 0
                nids = ints1[:, 0]
                eids2 = np.vstack([eids] * nnodes_all).T.ravel()
                assert eids.min() > 0, eids.min()
                obj.element_node[itotal:itotal2, 0] = eids2
                obj.element_node[itotal:itotal2, 1] = nids

            # [loc, rsr, rsi, azsr, azsi, Asr, Asi, ssr, ssi]
            isave1 = [1, 3, 5, 7]
            isave2 = [2, 4, 6, 9]
            real_imag = apply_mag_phase(floats1, is_magnitude_phase, isave1, isave2)

            obj.data[obj.itime, itotal:itotal2, :] = real_imag
            obj.itotal = itotal2
            obj.ielement = ielement2
        else:
            n = oes_ctriax_complex_37(op2, data, obj,
                                      nelements,
                                      is_magnitude_phase)
    else:  # pragma: no cover
        raise RuntimeError(op2.code_information())
        # msg = op2.code_information()
        # raise NotImplementedError(msg)
        # return op2._not_implemented_or_skip(data, ndata, msg)
    return n, nelements, ntotal



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
