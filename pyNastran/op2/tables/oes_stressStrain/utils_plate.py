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

from pyNastran.op2.tables.oes_stressStrain.real.oes_plates import RealPlateStressArray, RealPlateStrainArray
from pyNastran.op2.tables.oes_stressStrain.complex.oes_plates import (
    ComplexPlateStressArray, ComplexPlateStrainArray)

from pyNastran.op2.tables.oes_stressStrain.complex.oes_plates_vm import (
    ComplexPlateVMStressArray, ComplexPlateVMStrainArray)

from pyNastran.op2.tables.oes_stressStrain.random.oes_plates import RandomPlateStressArray, RandomPlateStrainArray
from pyNastran.op2.tables.oes_stressStrain.random.oes_plates_vm import RandomPlateVMStressArray, \
    RandomPlateVMStrainArray

if TYPE_CHECKING:  # pragma: no cover
    from pyNastran.op2.op2 import OP2


def oes_cquad4_33(op2: OP2, data, ndata: int, dt, is_magnitude_phase: bool,
                  result_type: int, prefix: str, postfix: str) -> tuple[int, Any, Any]:
    """
    reads stress/strain for element type:
     - 33 : CQUAD4-centroidal
     - 64 : CQUAD4-centroidal (nastran95)
     - 228 : CQUADR-centroidal

    """
    # n = 0
    factor = op2.factor
    size = op2.size
    # print('_oes_cquad4_33')
    if op2.element_type == 33:
        etype = 'cquad4'
    elif op2.element_type == 228:
        etype = 'cquadr'
        assert op2.num_wide in [17, 15], op2.code_information()
    else:  # pragma: no cover
        raise NotImplementedError(op2.code_information())

    if op2.is_stress:
        obj_vector_real = RealPlateStressArray
        obj_vector_complex = ComplexPlateStressArray
        result_name = f'{prefix}{etype}_stress{postfix}'
    else:
        obj_vector_real = RealPlateStrainArray
        obj_vector_complex = ComplexPlateStrainArray
        result_name = f'{prefix}{etype}_strain{postfix}'

    if op2._results.is_not_saved(result_name):
        return ndata, None, None
    op2._results._found_result(result_name)
    slot = op2.get_result(result_name)

    numwide_real = 17
    sort_method = op2.sort_method
    if result_type == 0 and op2.num_wide == 17:  # real
        ntotal = 68 * factor  # 4*17
        nelements = ndata // ntotal
        nlayers = nelements * 2  # 2 layers per node
        # op2.log.info(f'CQUAD4-33: len(data)={ndata} numwide={op2.num_wide} nelements={nelements} nlayers={nlayers}')

        auto_return, is_vectorized = op2._create_oes_object4(
            nlayers, result_name, slot, obj_vector_real)
        if auto_return:
            op2._data_factor = 2  # number of "layers" for an element
            return nelements * ntotal, None, None

        obj = op2.obj
        if op2.use_vector and is_vectorized and sort_method == 1:
            nnodes_expected = 2
            n = nelements * ntotal
            ielement = obj.ielement
            ielement2 = ielement + nelements
            itotal = obj.itotal
            itotal2 = itotal + nelements * nnodes_expected
            obj._times[obj.itime] = dt
            if obj.itime == 0:
                ints = np.frombuffer(data, dtype=op2.idtype8)
                ints1 = ints.reshape(nelements, numwide_real)
                eids = ints1[:, 0] // 10
                eids = np.vstack([eids, eids]).T.ravel()
                assert eids.min() > 0, eids.min()
                obj.element_node[itotal:itotal2, 0] = eids

            floats = np.frombuffer(data, dtype=op2.fdtype8).reshape(nelements, numwide_real)[:, 1:]

            # fd, sx, sy, txy, angle, major, minor, max_shear
            floats1 = floats.reshape(nelements * nnodes_expected, 8)
            obj.data[obj.itime, itotal:itotal2, :] = floats1.copy()
            obj.itotal = itotal2
            obj.ielement = ielement2
        else:
            if is_vectorized and op2.use_vector:  # pragma: no cover
                op2.log.debug(f'vectorize centroidal quad: {op2.element_name}-{op2.element_type} real '
                              f'SORT{sort_method}')
            n = oes_quad4_33_real_17(op2, data, obj, ntotal, nelements, dt)
        if op2.is_sort1:
            assert obj.element_node[:, 0].min() > 0, obj.element_node[:, 0].shape

    elif result_type == 1 and op2.num_wide == 15:  # imag
        # op2.to_nx(f' because CQUAD4-33 (numwide=15) was found')
        # nnodes = 0  # centroid + 4 corner points
        ntotal = op2.num_wide * size
        # op2.log.info(f'CQUAD4-33: len(data)={ndata} numwide={op2.num_wide} nelements={nelements} nlayers={nlayers}')

        nelements = ndata // ntotal
        nlayers = nelements * 2
        auto_return, is_vectorized = op2._create_oes_object4(
            nlayers, result_name, slot, obj_vector_complex)
        if auto_return:
            op2._data_factor = 2
            return nelements * ntotal, None, None

        obj = op2.obj
        if op2.use_vector and is_vectorized and sort_method == 1:
            n = nelements * ntotal
            nnodes_all = 1
            itotal = obj.itotal
            itotal2 = itotal + 2 * nelements * nnodes_all
            ielement = obj.ielement
            ielement2 = ielement + nelements

            floats = np.frombuffer(data, dtype=op2.fdtype).reshape(nelements, 15 * nnodes_all)
            floats1 = floats[:, 1:].reshape(nelements * nnodes_all * 2, 7).copy()
            obj._times[obj.itime] = dt
            if obj.itime == 0:
                ints = np.frombuffer(data, dtype=op2.idtype).reshape(nelements, 15 * nnodes_all).copy()
                eids = ints[:, 0] // 10
                ints[:, 0] = 0
                ints1 = ints.reshape(nelements * nnodes_all, 15)
                nids = ints[:, 0]
                assert eids.min() > 0, eids.min()
                eids2 = np.vstack([eids, eids]).T.ravel()
                nids2 = np.vstack([nids, nids]).T.ravel()
                obj.element_node[itotal:itotal2, 0] = eids2
                obj.element_node[itotal:itotal2, 1] = nids2

            # [fd, sxr, sxi, syr, syi, txyr, txyi]
            isave1 = [1, 3, 5]
            isave2 = [2, 4, 6]
            real_imag = apply_mag_phase(floats1, is_magnitude_phase, isave1, isave2)

            obj.fiber_distance[itotal:itotal2] = floats1[:, 0]
            obj.data[obj.itime, itotal:itotal2, :] = real_imag
            obj.itotal = itotal2
            obj.ielement = ielement2
        else:
            if is_vectorized and op2.use_vector:  # pragma: no cover
                op2.log.debug(f'vectorize CQUAD4-33 imag SORT{sort_method}')

            n = oes_cquad4_33_complex_15(
                op2, data, obj,
                nelements, ntotal,
                is_magnitude_phase)

    elif result_type in [1, 2] and op2.num_wide == 9:  # random msc
        # _oes_cquad4 is the same as _oes_ctria3
        element_id = op2.nonlinear_factor
        obj_vector_random = RandomPlateStressArray if op2.is_stress else RandomPlateStrainArray
        op2.data_code['nonlinear_factor'] = element_id

        if op2._results.is_not_saved(result_name):
            op2._data_factor = 2
            return ndata, None, None
        op2._results._found_result(result_name)
        slot = op2.get_result(result_name)

        ntotal = 36 * factor  # 4*9
        nelements = ndata // ntotal
        nlayers = nelements * 2
        nnodes_expected = 1
        # if op2.table_name_str.startswith('OSTRRMS'):
        # print(f'{op2.table_name_str} {result_name}: {nelements} {ntotal}; sort_method={sort_method}')

        auto_return, is_vectorized = op2._create_oes_object4(
            nlayers, result_name, slot, obj_vector_random)
        if auto_return:
            op2._data_factor = 2
            return nelements * ntotal, None, None

        obj = op2.obj
        if op2.use_vector and is_vectorized:
            n = nelements * ntotal
            itotal = obj.itotal
            itotal2 = itotal + nelements * 2

            if sort_method == 1:
                ielement = obj.ielement
                ielement2 = ielement + nelements

                obj._times[obj.itime] = dt
                if obj.itime == 0:
                    ints = np.frombuffer(data, dtype=op2.idtype)
                    ints1 = ints.reshape(nelements, 9)
                    eids = ints1[:, 0] // 10
                    assert eids.min() > 0, eids.min()
                    # print(eids)
                    eids2 = np.vstack([eids, eids]).T.ravel()
                    # print(eids2)

                    # TODO: what about layer 1/2?
                    # print(op2.code_information())
                    # print(f'eids.shape={eids.shape} obj.element_node.shape={obj.element_node.shape}')
                    obj.element_node[itotal:itotal2, 0] = eids2

                floats = np.frombuffer(data, dtype=op2.fdtype).reshape(nelements, 9)[:, 1:]
                # fd1, sx1, sy1, txy1, fd2, fx2, fy2, txy2
                floats2 = floats.reshape(nelements * nnodes_expected, 8)

                # [eid_device, fd1, sx1, sy1, txy1,
                #             fd2, sx2, sy2, txy2,]
                nf2 = floats2.shape[0]
                floats3 = floats2.reshape(nf2 * 2, 4)

                # print(obj)
                # print(op2.code_information())
                obj.fiber_distance[itotal:itotal2] = floats3[:, 0].copy()
                obj.data[obj.itime, itotal:itotal2, :] = floats3[:, 1:].copy()
                obj.itotal = itotal2
                obj.ielement = ielement2

            elif sort_method == 2 and op2._analysis_code_fmt == b'f':
                ielement = obj.itime
                ie_upper = 2 * ielement
                ie_lower = 2 * ielement + 1

                obj.element_node[ie_upper, 0] = dt
                obj.element_node[ie_lower, 0] = dt
                # obj._times[obj.itime] = dt

                floats = np.frombuffer(data, dtype=op2.fdtype).reshape(nelements, 9)  # [:, 1:]
                # itime is actually ielement
                # we grab the element id from the ints for all times
                if op2._analysis_code_fmt == b'i' and obj.itime == 0:
                    # print('analysis_code ', op2.analysis_code, op2._analysis_code_fmt)
                    ints = np.frombuffer(data, dtype=op2.idtype).reshape(nelements, 9)
                    eids = ints[:, 0] // 10
                    # nids = np.zeros(len(eids), dtype='int32')
                    # print(eids)
                    # eids = np.vstack([eids, nids]).T.ravel()
                    # print(eids.shape)
                    # print(eids)
                    # print(obj.element)
                    # assert eids.min() > 0, eids.min()
                    # obj.element[itotal:itotal2, 0] = eids
                    obj._times[itotal:itotal2] = eids
                    aaa
                elif op2._analysis_code_fmt == b'f' and obj.itime == 0:
                    # print(floats[:, 0])
                    # print(floats[:, 0].shape, obj._times.shape)
                    obj._times[itotal:itotal2] = floats[:, 0]

                floats1 = floats[:, 1:]
                # print(floats1)
                # print(floats1.shape)
                # fd, sx, sy, txy,
                floats2 = floats1.reshape(nelements * nnodes_expected, 8)
                nf2 = floats2.shape[0]
                # reshape it into 2 layers
                floats3 = floats2.reshape(nf2 * 2, 4)
                # we only need to grab the first two fiber/curvature values
                # as they're duplicated many times for the same element
                obj.fiber_distance[2 * obj.itime:2 * obj.itime + 2] = floats3[:2, 0].copy()
                # we apply the data across 2 rows because we have 2 layers
                obj.data[:, ie_upper, :] = floats3[::2, 1:].copy()
                obj.data[:, ie_lower, :] = floats3[1::2, 1:].copy()
            else:
                raise NotImplementedError(op2.code_information())
            obj.itotal = itotal2
            # obj.ielement = ielement2
        else:
            if is_vectorized and op2.use_vector:  # pragma: no cover
                op2.log.debug(f'vectorize CQUAD4-33 random numwide=9 SORT{sort_method}')
            n = oes_cquad4_33_random_9(op2, data, obj, nelements, ntotal)

    elif result_type in [1, 2] and op2.num_wide == 11:  # random
        # 2 FD1 RS Z1 = Fibre Distance
        # 3 SX1 RS Normal in x at Z1
        # 4 SY1 RS Normal in y at Z1
        # 5 TXY1 RS Shear in xy at Z1
        # 6 RMSVM1 RS RMS von Mises at Z1

        # 7 FD2 RS Z2 = Fibre Distance
        # 8 SX2 RS Normal in x at Z2
        # 9 SY2 RS Normal in y at Z2
        # 10 TXY2 RS Shear in xy at Z2
        # 11 RMSVM2 RS RMS von Mises at Z2

        element_id = op2.nonlinear_factor
        obj_vector_random = RandomPlateVMStressArray if op2.is_stress else RandomPlateVMStrainArray
        op2.data_code['nonlinear_factor'] = element_id

        if op2._results.is_not_saved(result_name):
            return ndata, None, None
        op2._results._found_result(result_name)
        slot = op2.get_result(result_name)

        ntotal = 44 * factor  # 4*11
        nelements = ndata // ntotal
        nlayers = nelements * 2
        nnodes_expected = 1

        auto_return, is_vectorized = op2._create_oes_object4(
            nlayers, result_name, slot, obj_vector_random)
        if auto_return:
            op2._data_factor = 2
            return nelements * ntotal, None, None

        obj = op2.obj
        if op2.use_vector and is_vectorized and 0:  # pragma: no cover
            n = nelements * 4 * op2.num_wide
            ielement = obj.ielement
            ielement2 = ielement + nelements
            itotal = obj.itotal
            itotal2 = itotal + nelements * nnodes_expected
            obj._times[obj.itime] = dt
            if obj.itime == 0:
                ints = np.frombuffer(data, dtype=op2.idtype)
                ints1 = ints.reshape(nelements, 9)
                eids = ints1[:, 0] // 10
                print(eids)
                eids = np.vstack([eids, eids]).T.ravel()
                print(eids.shape)
                print(eids)
                print(obj.element)
                assert eids.min() > 0, eids.min()
                obj.element[itotal:itotal2, 0] = eids

            floats = np.frombuffer(data, dtype=op2.fdtype).reshape(nelements, 11)[:, 1:]
            print(floats.shape)
            # fd, sx, sy, txy,
            floats1 = floats.reshape(nelements * nnodes_expected, 10)
            obj.data[obj.itime, itotal:itotal2, :] = floats1.copy()
            obj.itotal = itotal2
            obj.ielement = ielement2
        else:
            n = oes_cquad4_33_random_vm_11(op2, data, obj, nelements, ntotal)

    elif result_type == 1 and op2.num_wide == 17 and op2.table_name in [b'OESVM1', b'OESVM2', b'OSTRVM1',
                                                                        b'OSTRVM2']:  # freq
        # Table of element stresses for frequency response analysis that includes
        # von Mises stress output in SORT1 format.
        element_id = op2.nonlinear_factor
        if op2.is_stress:  # TODO: add new complex type
            obj_vector_complex = ComplexPlateVMStressArray
        else:
            obj_vector_complex = ComplexPlateVMStrainArray
        op2.data_code['nonlinear_factor'] = element_id

        if op2._results.is_not_saved(result_name):
            return ndata, None, None
        op2._results._found_result(result_name)
        slot = op2.get_result(result_name)

        ntotal = 68 * factor  # 4*17
        nelements = ndata // ntotal
        nlayers = nelements * 2
        nnodes_expected = 1
        # op2.log.info(f'CQUAD4-33: len(data)={ndata} numwide={op2.num_wide} nelements={nelements} nlayers={nlayers}')

        auto_return, is_vectorized = op2._create_oes_object4(
            nlayers, result_name, slot, obj_vector_complex)
        if auto_return:
            op2._data_factor = 2
            return nelements * ntotal, None, None

        # self.show_data(data)
        #   ELEMENT      FIBER                                     - STRESSES IN ELEMENT  COORDINATE SYSTEM -
        #      ID.       DISTANCE              NORMAL-X                       NORMAL-Y                      SHEAR-XY               VON MISES
        # 0     101  -5.000000E-01  -8.152692E-01 /  0.0           -1.321875E+00 /  0.0           -3.158517E+00 /  0.0            5.591334E+00
        #            5.000000E-01   1.728573E+00 /  0.0           -7.103837E+00 /  0.0            2.856040E+00 /  0.0            9.497519E+00
        # floats  = (1011,
        # -0.5, -0.8152692317962646, 0.0, -1.321874737739563, 0.0, -3.1585168838500977, 0.0, 5.591334342956543,
        # 0.5,   1.7285730838775635, 0.0, -7.103837490081787, 0.0,  2.8560397624969482, 0.0, 9.497518539428711)
        obj = op2.obj
        if is_vectorized and op2.use_vector and op2.sort_method == 1 and 0:  # pragma: no cover
            raise NotImplementedError(op2.table_name_str)
        else:
            if is_vectorized and op2.use_vector:  # pragma: no cover
                op2.log.debug('vectorize CQUAD4-33 complex '
                              f'{op2.table_name_str} SORT{op2.sort_method}')
            n = oes_cquad4_33_complex_vm_17(op2, data, obj, nelements, ntotal,
                                            is_magnitude_phase)

    else:  # pragma: no cover
        raise RuntimeError(op2.code_information())
    assert op2.obj.element_name == op2.element_name, op2.obj
    assert n > 0
    return n, nelements, ntotal


def oes_cquad4_144(op2: OP2, data: bytes, ndata: int, dt, is_magnitude_phase: bool,
                   result_type: int, prefix: str, postfix: str) -> tuple[int, Any, Any]:
    """
    reads stress/strain for element type:
     - 64 : CQUAD8
     - 70 : CTRIAR
     - 75 : CTRIA6
     - 82 : CQUADR
     - 144 : CQUAD4-bilinear

    """
    n = 0
    size = op2.size
    factor = op2.factor

    etype_map = {
        # element_type : (element_base, nnodes_expected, element_name)
        64: ('cquad8', 4, 'CQUAD8'),
        70: ('ctriar', 3, 'CTRIAR'),
        75: ('ctria6', 3, 'CTRIA6'),
        82: ('cquadr', 4, 'CQUADR'),
        144: ('cquad4', 4, 'CQUAD4-bilinear'),
    }
    if op2.is_stress:
        stress_strain = 'stress'
        obj_vector_real = RealPlateStressArray
        obj_vector_complex = ComplexPlateStressArray
        obj_vector_random = RandomPlateStressArray
    else:
        stress_strain = 'strain'
        obj_vector_real = RealPlateStrainArray
        obj_vector_complex = ComplexPlateStrainArray
        obj_vector_random = RandomPlateStrainArray

    # centroid not included in nnodes
    element_base, nnodes, element_name = etype_map[op2.element_type]
    # if prefix == '' and postfix == '':
    # prefix = stress_strain + '.'

    # stress.cquad4_stress
    result_name = prefix + f'{element_base}_{stress_strain}' + postfix

    if op2._results.is_not_saved(result_name):
        return ndata, None, None
    op2._results._found_result(result_name)
    log = op2.log

    nnodes_all = nnodes + 1  # adding the centroid

    slot = op2.get_result(result_name)
    numwide_real = 2 + 17 * nnodes_all
    numwide_imag = 2 + 15 * nnodes_all
    numwide_random = 2 + 9 * nnodes_all

    # numwide_imag2 = 2 + 16 * nnodes_all
    # print('%s real=%s imag=%s imag2=%s random=%s' % (
    # op2.element_name, numwide_real, numwide_imag, numwide_imag2, numwide_random
    # ))
    # etype = op2.element_name
    # grid_center = 'CEN/%i' % nnodes

    # OESVM1/2 (complex)
    #   87 - CQUAD8

    # OSTRNO1
    #   47 - CQUAD4

    # CQUAD8 real=87 imag=77 imag2=82 random=47
    # CQUAD4 ???=
    sort_method = op2.sort_method
    element_name_type = f'{op2.element_name}-{op2.element_type}'
    # print(op2.code_information())

    # if result_type == 0 and op2.num_wide == numwide_real and op2._nastran_format == 'optistruct' and op2.element_type == 75:  # real
    # could be a quad or tri
    # x = 1

    if result_type == 0 and op2.num_wide == numwide_real:  # real
        ntotal = 4 * (2 + 17 * nnodes_all) * factor
        nelements = ndata // ntotal
        assert ndata % ntotal == 0
        nlayers = 2 * nelements * nnodes_all  # 2 layers per node

        auto_return, is_vectorized = op2._create_oes_object4(
            nlayers, result_name, slot, obj_vector_real)
        if auto_return:
            op2._data_factor = 2 * nnodes_all  # number of "layers" for an element
            return nelements * ntotal, None, None

        obj = op2.obj
        # print('dt=%s, itime=%s' % (obj.itime, dt))
        if op2.use_vector and is_vectorized and sort_method == 1:
            # self.itime = 0
            # self.ielement = 0
            # self.itotal = 0
            # self.ntimes = 0
            # self.nelements = 0
            n = nelements * ntotal

            istart = obj.itotal
            iend = istart + nlayers
            obj._times[obj.itime] = dt

            if obj.itime == 0:
                ints = np.frombuffer(data, dtype=op2.idtype8).reshape(nelements, numwide_real)
                ints1 = ints[:, 2:].reshape(nlayers // 2, 17)[:, 0].reshape(nelements, nnodes_all).copy()
                ints1[:, 0] = 0.
                nids = ints1.ravel()

                eids = ints[:, 0] // 10
                eids2 = np.array([eids] * (nnodes_all * 2), dtype=op2.idtype8).T.ravel()
                nids2 = np.vstack([nids, nids]).T.ravel()
                obj.element_node[istart:iend, 0] = eids2
                obj.element_node[istart:iend, 1] = nids2
                # assert obj.element_node[:iend, 0].min() > 0, eids2
                if obj.nonlinear_factor is not None:
                    float_mask = np.arange(nelements * numwide_real, dtype=np.int32).reshape(nelements, numwide_real)
                    float_mask1 = float_mask[:, 2:].reshape(nlayers // 2, 17)[:, 1:].reshape(nlayers, 8)
                    obj.float_mask = float_mask1

            if obj.nonlinear_factor is not None:
                results = np.frombuffer(data, dtype=op2.fdtype8)[obj.float_mask].copy()
            else:
                floats = np.frombuffer(data, dtype=op2.fdtype8).reshape(nelements, numwide_real)
                floats1 = floats[:, 2:].reshape(nlayers // 2, 17)
                results = floats1[:, 1:].reshape(nlayers, 8).copy()

            # [fiber_dist, oxx, oyy, txy, angle, majorP, minorP, ovm]
            obj.data[obj.itime, istart:iend, :] = results
            assert obj.element_node[:, 0].min() > 0, obj.element_node[:, 0]
        else:
            if is_vectorized and op2.use_vector:  # pragma: no cover
                log.debug(f'vectorize nodal shell: {element_name_type}... real SORT{sort_method}')

            n = oes_cquad4_144_real(op2, data, ndata, obj,
                                    nelements, nnodes, dt)
        if op2.is_sort1:
            assert obj.element_node[:, 0].min() > 0, obj.element_node[:, 0]

    elif result_type == 1 and op2.num_wide == numwide_imag:  # complex
        ntotal = numwide_imag * 4 * factor
        # assert op2.num_wide * 4 == ntotal, 'numwide*4=%s ntotal=%s' % (op2.num_wide*4, ntotal)
        nelements = ndata // ntotal
        nlayers = nelements * 2 * nnodes_all
        # print(element_name_type)
        # print('ndata', ndata)
        # print('ntotal', ntotal)
        # print('nelements', nelements)
        # print('nlayers', nlayers)

        auto_return, is_vectorized = op2._create_oes_object4(
            nlayers, result_name, slot, obj_vector_complex)
        if auto_return:
            op2._data_factor = 2 * nnodes_all
            # if op2.num_wide == 77:
            # print('ntotal =', ndata, ntotal, op2.num_wide)
            # print('nelements * ntotal =', nelements * ntotal)
            return nelements * ntotal, None, None
        obj = op2.obj

        if op2.use_vector and is_vectorized and sort_method == 1:
            n = nelements * ntotal
            itotal = obj.itotal
            itotal2 = itotal + nelements * (nnodes_all * 2)
            ielement = obj.ielement
            ielement2 = ielement + nelements

            floats = np.frombuffer(data, dtype=op2.fdtype8).reshape(nelements, numwide_imag)
            floats1 = floats[:, 2:].reshape(nelements * nnodes_all, 15)
            floats2 = floats1[:, 1:].reshape(nelements * nnodes_all * 2, 7).copy()
            obj._times[obj.itime] = dt
            if obj.itime == 0:
                ints = np.frombuffer(data, dtype=op2.idtype8).reshape(nelements, numwide_imag).copy()
                ints[:, 2] = 0  # set center node to 0
                ints1 = ints[:, 2:].reshape(nelements * nnodes_all, 15)
                eids = ints[:, 0] // 10
                nids = ints1[:, 0]
                eids2 = np.vstack([eids] * (nnodes_all * 2)).T.ravel()
                nids2 = np.vstack([nids, nids]).T.ravel()
                assert eids.min() > 0, eids.min()
                obj.element_node[itotal:itotal2, 0] = eids2
                obj.element_node[itotal:itotal2, 1] = nids2

            # [fd, sxr, sxi, syr, syi, txyr, txyi]
            isave1 = [1, 3, 5]
            isave2 = [2, 4, 6]
            real_imag = apply_mag_phase(floats2, is_magnitude_phase, isave1, isave2)

            obj.fiber_distance[itotal:itotal2] = floats2[:, 0]
            obj.data[obj.itime, itotal:itotal2, :] = real_imag
            obj.itotal = itotal2
            obj.ielement = ielement2
        else:
            if is_vectorized and op2.use_vector:  # pragma: no cover
                log.debug(f'vectorize CQUAD4-144/{element_name_type}... imag SORT{sort_method}')
            # nnodes_cquad4 = 5
            # nelements = 3
            # nlayers = nelements * nodes_cquad4 * 2 = 3*5*2 = 30
            # ntotal = nlayers
            n = oes_cquad4_144_complex_77(op2, data, obj,
                                          nelements, nnodes,
                                          dt, is_magnitude_phase)
    # elif op2.format_code == 1 and op2.num_wide == numwide_random: # random
    # msg = op2.code_information()
    # msg += '  numwide=%s numwide_real=%s numwide_imag=%s numwide_random=%s' % (
    # op2.num_wide, numwide_real, numwide_imag, numwide_random)
    # return op2._not_implemented_or_skip(data, ndata, msg), None, None
    elif result_type == 2 and op2.num_wide == numwide_random:  # random
        # 47 - CQUAD8-64
        # 38 - CTRIAR-70
        ntotal = op2.num_wide * 4 * factor
        nelements = ndata // ntotal
        assert ndata % ntotal == 0
        nlayers = 2 * nelements * nnodes_all  # 2 layers per node
        # op2.log.info(f'random quad-144 ntotal={ntotal} ndata={ndata} ntotal={ntotal} numwide={op2.num_wide} -> nelements={nelements}')

        # if op2.read_mode == 1:
        # msg = ''
        # return op2._not_implemented_or_skip(data, ndata, msg), None, None

        auto_return, is_vectorized = op2._create_oes_object4(
            nlayers, result_name, slot, obj_vector_random)
        if auto_return:
            op2._data_factor = 2 * nnodes_all  # number of "layers" for an element
            return nelements * ntotal, None, None

        obj = op2.obj
        # print('dt=%s, itime=%s' % (obj.itime, dt))
        if op2.use_vector and is_vectorized and 0:
            # self.itime = 0
            # self.ielement = 0
            # self.itotal = 0
            # self.ntimes = 0
            # self.nelements = 0
            n = nelements * ntotal

            istart = obj.itotal
            iend = istart + nlayers
            obj._times[obj.itime] = dt

            if obj.itime == 0:
                ints = np.frombuffer(data, dtype=op2.idtype).reshape(nelements, numwide_real)
                ints1 = ints[:, 2:].reshape(nlayers // 2, 17)[:, 0].reshape(nelements, nnodes_all).copy()
                ints1[:, 0] = 0.
                nids = ints1.ravel()

                eids = ints[:, 0] // 10
                eids2 = np.array([eids] * (nnodes_all * 2), dtype='int32').T.ravel()
                nids2 = np.vstack([nids, nids]).T.ravel()
                obj.element_node[istart:iend, 0] = eids2
                obj.element_node[istart:iend, 1] = nids2
                # assert obj.element_node[:iend, 0].min() > 0, eids2
                if obj.nonlinear_factor is not None:
                    float_mask = np.arange(nelements * numwide_real, dtype=np.int32).reshape(nelements, numwide_real)
                    float_mask1 = float_mask[:, 2:].reshape(nlayers // 2, 17)[:, 1:].reshape(nlayers, 8)
                    obj.float_mask = float_mask1

            if obj.nonlinear_factor is not None:
                results = frombuffer(data, dtype=op2.fdtype)[obj.float_mask].copy()
            else:
                floats = frombuffer(data, dtype=op2.fdtype).reshape(nelements, numwide_real)
                floats1 = floats[:, 2:].reshape(nlayers // 2, 17)
                results = floats1[:, 1:].reshape(nlayers, 8).copy()

            # [fiber_dist, oxx, oyy, txy, angle, majorP, minorP, ovm]
            obj.data[obj.itime, istart:iend, :] = results
        else:
            if is_vectorized and op2.use_vector:  # pragma: no cover
                log.debug(f'vectorize CQUAD4-144/{element_name_type}... random SORT{sort_method}')
            # numwide_random = 2 + 9 * nnodes_all
            n = oes_cquad4_144_random(op2, data, obj, nelements, nnodes, ndata)

        # if op2.read_mode == 1:
        # msg = ''
        # return op2._not_implemented_or_skip(data, ndata, msg), None, None
        ##self.show_data(data[:44])
        # ntotal = 44
        # struct1 = Struct(op2._endian + b'i4si 8f')
        # for i in ra
        # for i in range(20):
        # edata = data[n:n+ntotal]
        # out = struct1.unpack(edata)
        # self.show_data(edata)
        # print(out)
        # n += ntotal

        ## 47 - CQUAD8-64
        ##msg = op2.code_information()
        # msg = '%s-CQUAD4-numwide=%s format_code=%s;\n numwide_real=%s numwide_imag=%s numwide_random=%s' % (
        # op2.table_name_str, op2.num_wide, op2.format_code,
        # numwide_real, numwide_imag, numwide_random)
        # return op2._not_implemented_or_skip(data, ndata, msg), None, None
    elif op2.table_name in [b'OESVM1', b'OESVM2', b'OSTRVM1', b'OSTRVM2']:
        # 82  CQUADR -> 87
        # CQUAD8 sort_method=2 ntotal=348 nelements=3
        # CTRIA6 sort_method=2 ntotal=348 nelements=2

        msg = op2.code_information()
        if result_type == 1:  # complex
            # ndata = 3828
            # ???
            #
            # ndata=1044
            assert op2.num_wide in [70, 87], op2.code_information()
            ntotal = op2.num_wide * self.size  # 87*4
            nelements = ndata // ntotal
            nlayers = nelements * 2
            assert ndata % ntotal == 0

            obj_vector_complex = ComplexPlateVMStressArray if op2.is_stress else ComplexPlateVMStrainArray

            auto_return, is_vectorized = op2._create_oes_object4(
                nlayers, result_name, slot, obj_vector_complex)

            if auto_return:
                op2._data_factor = 2 * nnodes_all  # number of "layers" for an element
                # op2._data_factor = 2
                return nelements * ntotal, None, None
            # if op2.read_mode == 1:
            # return op2._not_implemented_or_skip(data, ndata, msg), None, None

            # self.show_data(data)
            # print(ndata, ntotal)
            obj = op2.obj
            n = oes_cquad4_complex_vm_87(op2, data, obj, nelements, nnodes_all,
                                         is_magnitude_phase)

        # if result_type == 1 and op2.num_wide in [70, 87]:
        # 70 - CTRIA6-75
        # 87 - CQUAD4-144
        # pass
        else:
            # msg = (f'skipping {op2.table_name_str}-{op2.element_name}: numwide={op2.num_wide} '
            # f'result_type={op2.result_type} (complex);\n numwide_real={numwide_real} '
            # f'numwide_imag={numwide_imag} numwide_random={numwide_random}')
            # print(msg)
            raise NotImplementedError(msg)
        # return op2._not_implemented_or_skip(data, ndata, msg), None, None

    # elif op2.format_code in [2, 3] and op2.num_wide == 70:
    ## 87 - CQUAD4-144
    ##msg = op2.code_information()
    # msg = '%s-CTRIA6-numwide=%s numwide_real=%s numwide_imag=%s numwide_random=%s' % (
    # op2.table_name_str, op2.num_wide, numwide_real, numwide_imag, numwide_random)
    # return op2._not_implemented_or_skip(data, ndata, msg), None, None

    elif result_type == 2 and op2.table_name in [b'OESXRMS1', b'OESXNO1']:
        # CTRIA6-75  numwide=46  C:\MSC.Software\simcenter_nastran_2019.2\tpl_post2\tr1081x.op2
        # CQUAD8-64  numwide=64  C:\MSC.Software\simcenter_nastran_2019.2\tpl_post2\tr1081x.op2

        # corners + centroid
        if op2.element_type == 75:  # CTRIA6  (numwide=46)
            # C:\MSC.Software\simcenter_nastran_2019.2\tpl_post2\tr1081x.op2
            nnodes = 4
        elif op2.element_type in [64, 144]:
            # CQUAD8     (numwide=57) C:\MSC.Software\simcenter_nastran_2019.2\tpl_post2\tr1081x.op2
            # CQUAD4-144 (numwide=57) C:\MSC.Software\simcenter_nastran_2019.2\tpl_post2\plate_111o.op2
            nnodes = 5
        else:
            raise RuntimeError(op2.code_information())

        assert 2 + 11 * nnodes == op2.num_wide, op2.code_information()
        obj_vector_random = RandomPlateVMStressArray if op2.is_stress else RandomPlateVMStrainArray

        ntotal = op2.num_wide * size
        nelements = ndata // ntotal
        nlayers = nelements * nnodes_all * 2
        assert ndata % ntotal == 0
        # print(f'selement_name={op2.element_name} sort_method={op2.sort_method} ntotal={ntotal} num_wide={op2.num_wide} nelements={nelements} ndata={ndata} nlayers={nlayers}')
        # print('nnodes_all =', nnodes_all)
        # 57 = 2 + 11 * 5
        # ints    = (64011, 'CEN/',
        #           5, -1127428915, 1168048724, 1174717287, 1170166698, 1179261635, 0.025, 1167390785, 1175199013, 1169871798, 1179273875,
        #           80000, -1127428915, 1144450062, 1161227278, 1155248064, 1165738662, 0.025, 1149212546, 1165989762, 1154256422, 1167150968,
        #           80002, -1127428915, 1182003448, 1198780664, 1170939895, 1197489549, 0.025, 1182601858, 1199379074, 1169253917, 1197953979,
        #           80202, -1127428915, 1171034189, 1187612615, 1180541053, 1191448609, 0.025, 1169771716, 1186411799, 1181780048, 1191501390,
        #           80200, -1127428915, 1178323054, 1132616020, 1157628875, 1178722187, 0.025, 1176075105, 1152142782, 1150736003, 1177059240,
        #           64021, 'CEN/',
        #           5,
        # -1127428915, 1180051065, 1181683390, 1180538582, 1189027141,
        #  0.025, 1176941913, 1182077843, 1179498554, 1188072152,
        #           80200, -1127428915, 1179821552, 1170149790, 1168886008, 1181456870, 0.025, 1178168847, 1171555762, 1167658184, 1179691422,
        #           80202, -1127428915, 1180244047, 1197120002, 1184242951, 1198270294, 0.025, 1180604356, 1197449345, 1183498771, 1198177490,
        #           80402, -1127428915, 1192586977, 1184669070, 1184219971, 1194803525, 0.025, 1193976276, 1184750091, 1182770117, 1194901195,
        #           80400, -1127428915, 1199592178, 1184219800, 1179170385, 1198821390, 0.025, 1198216710, 1184079548, 1179988640, 1197735836,
        #           64111, 'CEN/',
        #           8,     -1127428915, 1176450101, 1187579168, 1179109877, 1190548891, 0.025, 1176466945, 1187779876, 1179033188, 1190633996,
        #           90000, -1127428915, 1163199692, 1169764570, 1159841955, 1171614247, 0.025, 1163848412, 1172543886, 1158998260, 1173065416,
        #           90002, -1127428915, 1183839433, 1201603141, 1152572706, 1200652206, 0.025, 1184231034, 1201902140, 1154931831, 1200916849,
        #           90202, -1127428915, 1156439069, 1187515149, 1193045713, 1200691352, 0.025, 1130954657, 1188796022, 1193252160, 1200929208,
        #           90200, -1127428915, 1188552094, 1155874499, 1172127542, 1189476821, 0.025, 1187333567, 1120409103, 1169694799, 1188368246,
        #           64121, 'CEN/',
        #           8,     -1127428915, 1188544334, 1178584960, 1188275834, 1196282134, 0.025, 1186800673, 1178575754, 1187794850, 1195571664,
        #           90200, -1127428915, 1189840387, 1170594118, 1173995469, 1190183343, 0.025, 1188877086, 1165971753, 1173556714, 1189607807,
        #           90202, -1127428915, 1140886249, 1189485297, 1193242194, 1200973710, 0.025, 1161620989, 1191030157, 1193534822, 1201299068,
        #           90402, -1127428915, 1176972259, 1174260793, 1195749260, 1202470249, 0.025, 1185547735, 1170629093, 1194910619, 1201965268,
        #           90400, -1127428915, 1202771906, 1185377334, 1175183109, 1201858966, 0.025, 1202359853, 1184055805, 1173072162, 1201506363)
        # floats  = (8.969851599989587e-41, 'CEN/',
        #           5,                      -0.025, 5088.291015625, 8496.8505859375, 6122.4580078125, 12934.6904296875,
        #                                    0.025, 4767.03173828125, 8967.2861328125, 5978.4638671875, 12946.6435546875,
        #           1.1210387714598537e-40, -0.025, 731.6883544921875, 2926.75341796875, 1757.4921875, 4028.16552734375,
        #                                    0.025, 1022.3673095703125, 4089.46923828125, 1636.442138671875, 4649.93359375,
        #           1.1210667974291402e-40, -0.025, 15612.2421875, 62448.96875, 6499.99560546875, 57405.55078125,
        #                                    0.025, 16196.626953125, 64786.5078125, 5676.76416015625, 59219.73046875,
        #           1.1238693943577898e-40, -0.025, 6546.03759765625, 25795.888671875, 14184.1220703125, 33808.12890625,
        #                                    0.025, 5929.595703125, 23450.544921875, 15394.078125, 34014.3046875,
        #           1.1238413683885033e-40, -0.025, 12018.107421875, 260.6978759765625, 2048.237060546875, 12407.8857421875,
        #                                    0.025, 9822.8447265625, 1378.429443359375, 1206.7034912109375, 10783.9140625,
        #           8.971252898453911e-41, 'CEN/',
        #           5,                      -0.025, 13705.6181640625, 15299.685546875, 14181.708984375, 28558.634765625,
        #                                    0.025, 10669.3369140625, 15684.8935546875, 13166.056640625, 26693.421875,
        #           1.1238413683885033e-40, -0.025, 13481.484375, 6114.2021484375, 5497.12109375, 15078.474609375,
        #                                    0.025, 11867.5146484375, 6800.7119140625, 4897.59765625, 13354.404296875,
        #           1.1238693943577898e-40, -0.025, 13894.0771484375, 55962.0078125, 19214.513671875, 60455.3359375,
        #                                    0.025, 14245.94140625, 57248.50390625, 17761.037109375, 60092.8203125,
        #           1.1266719912864394e-40, -0.025, 38254.87890625, 20046.77734375, 19169.630859375, 46913.26953125,
        #                                    0.025, 43681.828125, 20205.021484375, 16360.9423828125, 47294.79296875,
        #           1.126643965317153e-40,  -0.025, 65701.890625, 19169.296875, 12845.5791015625, 62608.0546875,
        #                                    0.025, 60246.0234375, 18895.3671875, 13644.65625, 58367.609375,
        #           8.983864584632835e-41, 'CEN/',
        #           8,                      -0.025, 10189.0517578125, 25730.5625, 12786.4892578125, 31530.802734375,
        #                                    0.025, 10205.5009765625, 26122.5703125, 12711.59765625, 31697.0234375,
        #           1.2611686178923354e-40, -0.025, 3408.2998046875, 5926.1064453125, 2588.539794921875, 6829.26904296875,
        #                                    0.025, 3566.6787109375, 7283.1943359375, 2382.5595703125, 7537.84765625,
        #           1.2611966438616219e-40, -0.025, 18426.392578125, 81412.5390625, 1430.910400390625, 73983.359375,
        #                                    0.025, 19191.23828125, 83748.46875, 1718.8895263671875, 76050.8828125,
        #           1.2639992407902715e-40, -0.025, 1902.8785400390625, 25605.525390625, 40046.81640625, 74289.1875,
        #                                    0.025, 232.99855041503906, 28107.23046875, 40853.25, 76147.4375,
        #           1.263971214820985e-40,  -0.025, 27630.80859375, 1833.9613037109375, 7079.9013671875, 29436.916015625,
        #                                    0.025, 25250.873046875, 100.04308319091797, 5892.03857421875, 27271.73046875,
        #           8.98526588309716e-41, 'CEN/',
        #           8,                     -0.025, 27615.65234375, 12273.875, 27091.23828125, 52689.0859375,
        #                                   0.025, 24210.064453125, 12264.884765625, 26151.81640625, 49913.8125,
        #           1.263971214820985e-40, -0.025, 30147.005859375, 6331.1591796875, 7991.97509765625, 30816.841796875,
        #                                   0.025, 28265.55859375, 4085.072509765625, 7777.7392578125, 29692.748046875,
        #           1.2639992407902715e-40,-0.025, 514.1704711914062, 29453.470703125, 40814.3203125, 76495.109375,
        #                                   0.025, 3022.874267578125, 32470.775390625, 41957.3984375, 79036.96875,
        #           1.2668018377189211e-40,-0.025, 10698.9716796875, 8121.52783203125, 50607.546875, 88186.8203125,
        #                                   0.025, 21762.919921875, 6348.23681640625, 47331.60546875, 84241.65625,
        #           1.2667738117496346e-40,-0.025, 90543.515625, 21430.10546875, 8951.7548828125, 83411.171875,
        #                                   0.025, 87324.3515625, 18848.994140625, 7541.1416015625, 80656.4609375)
        # if op2.read_mode == 2:
        # self.show_data(data)
        # ddd
        auto_return, is_vectorized = op2._create_oes_object4(
            nlayers, result_name, slot, obj_vector_random)
        if auto_return:
            op2._data_factor = 2 * nnodes_all  # number of "layers" for an element
            return nelements * ntotal, None, None

        obj = op2.obj
        n = oes_cquad4_random_vm_57(op2, data, op2.obj, nelements, ntotal, nnodes,
                                    dt)
    elif result_type == 1 and op2.table_name in [b'OESPSD1'] and op2.element_type == 82 and op2.num_wide == 47:  # QUADR
        msg = op2.code_information()
        return op2._not_implemented_or_skip(data, ndata, msg), None, None
    else:  # pragma: no cover
        raise RuntimeError(op2.code_information())
    return n, nelements, ntotal


def oes_ctria3_74(op2: OP2, data, ndata: int, dt, is_magnitude_phase: bool,
                  result_type: int, prefix: str, postfix: str) -> tuple[int, Any, Any]:
    """
    reads stress/strain for element type:
     - 74 : CTRIA3-centroidal
     - 83 : CTRIA3-centroidal  (NASA 95)
     - 227: TRIAR-centroidal

    """
    # print('_oes_ctria3')
    # n = 0
    factor = op2.factor

    etype_map = {
        # element_type : (element_base, element_name)
        74: ('ctria3', 'CTRIA3'),
        83: ('ctria3', 'CTRIA3'),  # NASA-95
        227: ('ctriar', 'CTRIAR'),
    }
    if op2.is_stress:
        stress_strain = 'stress'
        obj_vector_real = RealPlateStressArray
        obj_vector_complex = ComplexPlateStressArray
    else:
        stress_strain = 'strain'
        obj_vector_real = RealPlateStrainArray
        obj_vector_complex = ComplexPlateStrainArray

    # if prefix == '' and postfix == '':
    # prefix = stress_strain + '.'

    element_base, element_name = etype_map[op2.element_type]
    # stress.ctria3_stress
    result_name = prefix + f'{element_base}_{stress_strain}' + postfix

    if op2._results.is_not_saved(result_name):
        return ndata, None, None
    op2._results._found_result(result_name)

    slot = op2.get_result(result_name)
    # print(op2.element_name, result_name, op2.format_code, op2.num_wide)
    # table_names = [
    #    b'OES1', b'OES1X', b'OES1X1', b'OSTR1X',
    #    b'OES2', b'OSTR2',
    #    b'OESVM1', b'OESVM2', b'OSTRVM1', b'OSTRVM2',
    #    b'OESPSD1', b'OESRMS1',
    #    b'OESPSD2',  b'OESATO2',  b'OESCRM2',  b'OESNO1',  b'OESXRMS1', b'OESXNO1',
    #    b'OSTRPSD2', b'OSTRATO2', b'OSTRCRM2', b'OSTRNO1', b'OSTRRMS1',
    # ]
    # assert op2.table_name in table_names, op2.table_name
    sort_method = op2.sort_method
    element_name_type = f'{op2.element_name}-{op2.element_type}'

    if op2.format_code in [1, 3] and op2.num_wide == 17:  # real
        ntotal = 68 * factor  # 4*17
        nelements = ndata // ntotal
        nlayers = nelements * 2  # 2 layers per node
        # if self.code[0] == 100:
        # print(f'\nnelements={nelements} nlayers={nlayers} {self.code in slot}')
        auto_return, is_vectorized = op2._create_oes_object4(
            nlayers, result_name, slot, obj_vector_real)
        if auto_return:
            op2._data_factor = 2
            return nelements * ntotal, None, None

        if op2.is_debug_file:
            op2.binary_debug.write('  [cap, element1, element2, ..., cap]\n')
            op2.binary_debug.write('  cap = %i  # assume 1 cap when there could have been multiple\n' % ndata)
            op2.binary_debug.write('  #elementi = [eid_device, fd1, sx1, sy1, txy1, angle1, major1, minor1, vm1,\n')
            op2.binary_debug.write('  #                        fd2, sx2, sy2, txy2, angle2, major2, minor2, vm2,]\n')
            op2.binary_debug.write('  nelements=%i; nnodes=1 # centroid\n' % nelements)

        obj = op2.obj
        if op2.use_vector and is_vectorized and sort_method == 1:
            nfields = 17 * nelements
            nbytes = nfields * 4
            itotal = obj.itotal
            iend = obj.itotal + nlayers

            itime = obj.itime
            if itime == 0:
                ints = frombuffer(data, dtype=op2.idtype8).reshape(nelements, 17)
                eids = ints[:, 0] // 10
                # ilayers = ints[:, 1]
                # ints2 = ints[:, 1:].reshape(nlayers, 8)
                assert eids.min() > 0, eids
                obj._times[obj.itime] = dt
                obj.element_node[itotal:iend:2, 0] = eids
                obj.element_node[itotal + 1:iend + 1:2, 0] = eids
                # obj.element_node[itotal:iend, 1] = 0

            floats = frombuffer(data, dtype=op2.fdtype8).reshape(nelements, 17)
            floats1 = floats[:, 1:].reshape(nlayers, 8).copy()
            obj.data[obj.itime, itotal:iend, :] = floats1
            obj._times[obj.itime] = dt
            obj.itotal += nlayers
            n = nbytes
        else:
            if is_vectorized and op2.use_vector:  # pragma: no cover
                op2.log.debug(f'vectorize centroidal tri: {element_name_type} real SORT{sort_method}')
            n = oes_ctria3_real_17(op2, data, obj,
                                   ntotal, nelements, dt)
        if op2.is_sort1:
            assert obj.element_node[:, 0].min() > 0, obj.element_node[:, 0]

    elif op2.format_code in [2, 3] and op2.num_wide == 15:  # imag
        ntotal = 60 * factor  # 4*15
        nelements = ndata // ntotal
        auto_return, is_vectorized = op2._create_oes_object4(
            nelements, result_name, slot, obj_vector_complex)
        if auto_return:
            op2._data_factor = 2
            # assert ntotal == op2.num_wide * 4
            return nelements * ntotal, None, None
        obj = op2.obj

        if op2.use_vector and is_vectorized and sort_method == 1:
            n = nelements * ntotal
            itotal = obj.itotal
            itotal2 = itotal + nelements * 2
            ielement = obj.ielement
            ielement2 = ielement + nelements

            floats = frombuffer(data, dtype=op2.fdtype8).reshape(nelements, 15)
            floats1 = floats[:, 1:].reshape(nelements * 2, 7).copy()
            obj._times[obj.itime] = dt
            if obj.itime == 0:
                ints = frombuffer(data, dtype=op2.idtype8).reshape(nelements, 15).copy()
                eids = ints[:, 0] // 10
                ints[:, 0] = 0
                unused_ints1 = ints.reshape(nelements, 15)
                nids = ints[:, 0]

                assert eids.min() > 0, eids.min()
                eids2 = np.vstack([eids, eids]).T.ravel()
                nids2 = np.vstack([nids, nids]).T.ravel()
                obj.element_node[itotal:itotal2, 0] = eids2
                obj.element_node[itotal:itotal2, 1] = nids2

            # [fd, sxr, sxi, syr, syi, txyr, txyi]
            isave1 = [1, 3, 5]
            isave2 = [2, 4, 6]
            real_imag = apply_mag_phase(floats1, is_magnitude_phase, isave1, isave2)

            obj.fiber_distance[itotal:itotal2] = floats1[:, 0]
            obj.data[obj.itime, itotal:itotal2, :] = real_imag
            obj.itotal = itotal2
            obj.ielement = ielement2
        else:
            if is_vectorized and op2.use_vector:  # pragma: no cover
                op2.log.debug(f'vectorize CTRIA3 imag SORT{sort_method}')

            n = oes_ctria3_complex_15(
                op2, data, obj,
                nelements, ntotal, is_magnitude_phase,
                element_name_type)
    # elif op2.format_code == 1 and op2.num_wide == 9: # random?
    # msg = op2.code_information()
    # return op2._not_implemented_or_skip(data, ndata, msg), None, None
    elif op2.format_code in [2, 3] and op2.num_wide == 17 and op2.table_name in [b'OESVM1', b'OESVM2', b'OSTRVM1',
                                                                                 b'OSTRVM2']:
        # freq:
        # # random; CTRIA3
        assert op2.table_name in [b'OESVM1', b'OESVM2', b'OSTRVM1', b'OSTRVM2'], op2.code_information()

        element_id = op2.nonlinear_factor
        if op2.is_stress:  # TODO: add new complex type
            obj_vector_complex = ComplexPlateVMStressArray
        else:
            obj_vector_complex = ComplexPlateVMStrainArray
        op2.data_code['nonlinear_factor'] = element_id

        if op2._results.is_not_saved(result_name):
            return ndata, None, None
        op2._results._found_result(result_name)
        slot = op2.get_result(result_name)

        ntotal = 68 * factor  # 17*4
        nelements = ndata // ntotal
        nlayers = nelements * 2
        nnodes_expected = 1

        auto_return, is_vectorized = op2._create_oes_object4(
            nlayers, result_name, slot, obj_vector_complex)
        if auto_return:
            op2._data_factor = 2
            return nelements * ntotal, None, None

        obj = op2.obj
        # C:\MSC.Software\simcenter_nastran_2019.2\tpl_post2\shlthk14.op2
        nelements = ndata // ntotal
        assert ndata % ntotal == 0

        #                   C O M P L E X   S T R E S S E S   I N   T R I A N G U L A R   E L E M E N T S   ( T R I A 3 )
        #                                                          (REAL/IMAGINARY)
        #
        #  ELEMENT       FIBER                                     - STRESSES IN ELEMENT  COORDINATE SYSTEM -
        #    ID.        DISTANCE              NORMAL-X                       NORMAL-Y                      SHEAR-XY               VON MISES
        # 0       1  -4.359080E+00  -1.391918E+00 /  2.474756E-03  -1.423926E+00 /  2.530494E-03   2.655153E-02 / -5.158625E-05   1.408948E+00
        #            4.359080E+00   1.391918E+00 / -2.474756E-03   1.423926E+00 / -2.530494E-03  -2.655153E-02 /  5.158625E-05   1.408948E+00
        n = oes_ctria3_complex_vm_17(op2, data, obj, nelements, ntotal, dt,
                                     is_magnitude_phase)
        assert n is not None, n

    elif op2.format_code in [1, 2, 3] and op2.num_wide == 11:  # random; CTRIA3
        # 2 FD1 RS Z1 = Fibre Distance
        # 3 SX1 RS Normal in x at Z1
        # 4 SY1 RS Normal in y at Z1
        # 5 TXY1 RS Shear in xy at Z1
        # 6 RMSVM1 RS RMS von Mises at Z1

        # 7 FD2 RS Z2 = Fibre Distance
        # 8 SX2 RS Normal in x at Z2
        # 9 SY2 RS Normal in y at Z2
        # 10 TXY2 RS Shear in xy at Z2
        # 11 RMSVM2 RS RMS von Mises at Z2

        element_id = op2.nonlinear_factor
        obj_vector_random = RandomPlateVMStressArray if op2.is_stress else RandomPlateVMStrainArray
        op2.data_code['nonlinear_factor'] = element_id

        if op2._results.is_not_saved(result_name):
            return ndata, None, None
        op2._results._found_result(result_name)
        slot = op2.get_result(result_name)

        ntotal = 44 * factor  # 4*11
        nelements = ndata // ntotal
        nlayers = nelements * 2
        nnodes_expected = 1

        auto_return, is_vectorized = op2._create_oes_object4(
            nlayers, result_name, slot, obj_vector_random)
        if auto_return:
            op2._data_factor = 2
            return nelements * ntotal, None, None

        obj = op2.obj
        if op2.use_vector and is_vectorized and 0:  # pragma: no cover
            n = nelements * 4 * op2.num_wide
            ielement = obj.ielement
            ielement2 = ielement + nelements
            itotal = obj.itotal
            itotal2 = itotal + nelements * nnodes_expected
            obj._times[obj.itime] = dt
            if obj.itime == 0:
                ints = frombuffer(data, dtype=op2.idtype)
                ints1 = ints.reshape(nelements, 9)
                eids = ints1[:, 0] // 10
                print(eids)
                eids = np.vstack([eids, eids]).T.ravel()
                print(eids.shape)
                print(eids)
                print(obj.element)
                assert eids.min() > 0, eids.min()
                obj.element[itotal:itotal2, 0] = eids

            floats = frombuffer(data, dtype=op2.fdtype).reshape(nelements, 11)[:, 1:]
            print(floats.shape)
            # fd, sx, sy, txy,
            floats1 = floats.reshape(nelements * nnodes_expected, 10)
            obj.data[obj.itime, itotal:itotal2, :] = floats1.copy()
            obj.itotal = itotal2
            obj.ielement = ielement2
        else:
            if is_vectorized and op2.use_vector and obj.itime == 0:  # pragma: no cover
                op2.log.debug(f'vectorize {element_name_type} random numwide=11 SORT{sort_method}')
            n = oes_ctria3_random_vm_11(op2, data, obj, nelements, ntotal)

    elif op2.format_code in [1, 2] and op2.num_wide == 9:  # random MSC stress/strain; CTRIA3
        # _oes_cquad4 is the same as _oes_ctria3
        element_id = op2.nonlinear_factor
        if op2.is_stress:
            obj_vector_random = RandomPlateStressArray
            # result_name = prefix + 'ctria3_stress' + postfix
        else:
            obj_vector_random = RandomPlateStrainArray
            # result_name = prefix + 'ctria3_strain' + postfix

        op2.data_code['nonlinear_factor'] = element_id

        if op2._results.is_not_saved(result_name):
            return ndata, None, None
        op2._results._found_result(result_name)
        slot = op2.get_result(result_name)

        ntotal = 36 * factor  # 4*9
        nelements = ndata // ntotal
        nlayers = nelements * 2
        nnodes_expected = 1

        auto_return, is_vectorized = op2._create_oes_object4(
            nlayers, result_name, slot, obj_vector_random)
        if auto_return:
            op2._data_factor = 2
            return nelements * ntotal, None, None

        obj = op2.obj
        if op2.use_vector and is_vectorized and 0:  # pragma: no cover
            n = nelements * 4 * op2.num_wide
            ielement = obj.ielement
            ielement2 = ielement + nelements
            itotal = obj.itotal
            itotal2 = itotal + nelements * nnodes_expected
            obj._times[obj.itime] = dt
            if obj.itime == 0:
                ints = frombuffer(data, dtype=op2.idtype)
                ints1 = ints.reshape(nelements, 9)
                eids = ints1[:, 0] // 10
                # print(eids)
                eids = np.vstack([eids, eids]).T.ravel()
                # print(eids.shape)
                # print(eids)
                # print(obj.element)
                assert eids.min() > 0, eids.min()
                obj.element[itotal:itotal2, 0] = eids

            floats = frombuffer(data, dtype=op2.fdtype).reshape(nelements, 9)[:, 1:]
            # print(floats.shape)
            # fd, sx, sy, txy,
            floats1 = floats.reshape(nelements * nnodes_expected, 8)
            obj.data[obj.itime, itotal:itotal2, :] = floats1.copy()
            obj.itotal = itotal2
            obj.ielement = ielement2
        else:
            if is_vectorized and op2.use_vector:  # pragma: no cover
                op2.log.debug(f'vectorize {element_name_type} random2 SORT{sort_method}')
            n = oes_ctria3_random_9(op2, data, obj, nelements, ntotal)

    else:  # pragma: no cover
        raise RuntimeError(op2.code_information())
    assert n is not None, op2.code_information()
    return n, nelements, ntotal


def oes_cquad4_33_complex_vm_17(op2: OP2, data: bytes,
                                obj: ComplexPlateVMStressArray | ComplexPlateVMStrainArray,
                                nelements: int, ntotal: int,
                                is_magnitude_phase: bool) -> int:
    """
    OESVM1/2 - Table of element stresses or strains with von Mises
    OSTRVM1/2  for frequency response results.

    """
    #self.log.info(f'nelements = {nelements}')
    n = 0
    struct1 = Struct(mapfmt(op2._endian + op2._analysis_code_fmt + b'16f', op2.size))
    cen = 0 # CEN/4
    add_sort_x = getattr(obj, 'add_sort' + str(op2.sort_method))

    for unused_i in range(nelements):
        edata = data[n:n+ntotal]
        out = struct1.unpack(edata)

        (eid_device,
         fd1, sx1r, sx1i, sy1r, sy1i, txy1r, txy1i, von_mises1,
         fd2, sx2r, sx2i, sy2r, sy2i, txy2r, txy2i, von_mises2) = out

        eid, dt = get_eid_dt_from_eid_device(
            eid_device, op2.nonlinear_factor, op2.sort_method)
        if op2.is_debug_file:
            op2.binary_debug.write('  eid=%i C=[%s]\n' % (
                eid, ', '.join(['%r' % di for di in out])))

        if is_magnitude_phase:
            sx1 = polar_to_real_imag(sx1r, sx1i)
            sx2 = polar_to_real_imag(sx2r, sx2i)
            sy1 = polar_to_real_imag(sy1r, sy1i)
            sy2 = polar_to_real_imag(sy2r, sy2i)
            txy1 = polar_to_real_imag(txy1r, txy1i)
            txy2 = polar_to_real_imag(txy2r, txy2i)
        else:
            sx1 = complex(sx1r, sx1i)
            sx2 = complex(sx2r, sx2i)
            sy1 = complex(sy1r, sy1i)
            sy2 = complex(sy2r, sy2i)
            txy1 = complex(txy1r, txy1i)
            txy2 = complex(txy2r, txy2i)
        #print(dt, eid, cen, sx1, sy1, txy1, max_shear1)
        #print(dt, eid, cen, sx2, sy2, txy2, max_shear2)
        add_sort_x(dt, eid, cen,
                   fd1, sx1, sy1, txy1, von_mises1,
                   fd2, sx2, sy2, txy2, von_mises2)
        n += ntotal
    return n


def oes_cquad4_33_complex_15(op2: OP2,
                             data: bytes,
                             obj: ComplexPlateStressArray | ComplexPlateStrainArray,
                             nelements: int, ntotal: int,
                             is_magnitude_phase: bool) -> int:
    """
    NX/MSC 2018.2
    2 FD1    RS Z1 = Fibre distance
    3 EX1R   RS Normal in x at Z1
    4 EX1I   RS Normal in x at Z1
    5 EY1R   RS Normal in y at Z1
    6 EY1I   RS Normal in y at Z1
    7 EXY1R  RS Shear in xy at Z1
    8 EXY1I  RS Shear in xy at Z1
    9 FD2    RS Z2 = Fibre distance
    10 EX2R  RS Normal in x at Z2
    11 EX2I  RS Normal in x at Z2
    12 EY2R  RS Normal in y at Z2
    13 EY2I  RS Normal in y at Z2
    14 EXY2R RS Shear in xy at Z2
    15 EXY2I RS Shear in xy at Z2
    """
    n = 0
    s1 = Struct(op2._endian + op2._analysis_code_fmt + b'14f')
    cen = 0  # 'CEN/4'
    add_sort_x = getattr(obj, 'add_sort' + str(op2.sort_method))

    for unused_i in range(nelements):
        edata = data[n:n+ntotal]  # 4*15=60
        n += ntotal
        out = s1.unpack(edata)  # 15
        (eid_device,
         fd1, sx1r, sx1i, sy1r, sy1i, txy1r, txy1i,
         fd2, sx2r, sx2i, sy2r, sy2i, txy2r, txy2i) = out

        eid, dt = get_eid_dt_from_eid_device(
            eid_device, op2.nonlinear_factor, op2.sort_method)
        if op2.is_debug_file:
            op2.binary_debug.write('  eid=%i C=%s\n' % (eid, str(out)))

        if is_magnitude_phase:
            sx1 = polar_to_real_imag(sx1r, sx1i)
            sx2 = polar_to_real_imag(sx2r, sx2i)
            sy1 = polar_to_real_imag(sy1r, sy1i)
            sy2 = polar_to_real_imag(sy2r, sy2i)
            txy1 = polar_to_real_imag(txy1r, txy1i)
            txy2 = polar_to_real_imag(txy2r, txy2i)
        else:
            sx1 = complex(sx1r, sx1i)
            sx2 = complex(sx2r, sx2i)
            sy1 = complex(sy1r, sy1i)
            sy2 = complex(sy2r, sy2i)
            txy1 = complex(txy1r, txy1i)
            txy2 = complex(txy2r, txy2i)

        add_sort_x(dt, eid, cen,
                   fd1, sx1, sy1, txy1,
                   fd2, sx2, sy2, txy2)
    return n


def oes_cquad4_random_vm_57(op2: OP2,
                            data: bytes,
                            obj: ComplexPlateVMStressArray | ComplexPlateVMStrainArray,
                            nelements: int, ntotal: int, nnodes: int,
                            dt: Any) -> int:
    """
    name    ntotal numwide
    CQUAD8  228    57
    CTRIA6  184
    """

    #ntotal = 228 # 57*4
    n = 0
    cen = 0
    if op2.size == 4:
        struct1 = Struct(op2._endian + b'i4si 10f')
        struct2 = Struct(op2._endian + b'i 10f')
    else:
        print('ntotal =', ntotal, 228*2)
        struct1 = Struct(op2._endian + b'i4si 10f')
        struct2 = Struct(op2._endian + b'i 10f')
        raise NotImplementedError('64-bit')

    nnodes_corners = nnodes - 1
    #print('random-57: nelements =', nelements)
    factor = op2.factor

    ntotal1 = 52 * factor
    ntotal2 = 44 * factor
    add_sort_x = getattr(obj, 'add_sort' + str(op2.sort_method))
    for unused_i in range(nelements):
        ni = 0
        edata = data[n:n+ntotal]
        #self.show_data(edata)
        #out = struct1.unpack(edata)  # len=17*4
        (eid_device, grid1, four,
         fd1, ox1, oy1, txy1, vm1,
         fd2, ox2, oy2, txy2, vm2) = struct1.unpack(edata[ni:ni+ntotal1])  # 13*4=52
        eid, dt = get_eid_dt_from_eid_device(
            eid_device, op2.nonlinear_factor, op2.sort_method)
        #print(eid, cen)
        add_sort_x(dt, eid, cen,
                   fd1, ox2, oy2, txy2, vm1,
                   fd2, ox2, oy2, txy2, vm2)
        ni += ntotal1
        #print(eid, grid1,
              #fd1, ox1, oy1, txy1, vm1,
              #fd2, ox2, oy2, txy2, vm2)

        for unused_jnode in range(nnodes_corners): # 4*44=176
            (grid,
             fd1, ox1, oy1, txy1, vm1,
             fd2, ox2, oy2, txy2, vm2) = struct2.unpack(edata[ni:ni+ntotal2])  # 11*4=44
            #print(eid, grid)
            ni += ntotal2
            #print(grid,
                  #fd1, ox1, oy1, txy1, vm1,
                  #fd2, ox2, oy2, txy2, vm2)
            add_sort_x(dt, eid, grid,
                       fd1, ox2, oy2, txy2, vm1,
                       fd2, ox2, oy2, txy2, vm2)
        assert ni == ntotal, f'ni={ni}'
        n += ntotal
        #print()
    #assert n == ndata
    #n = self._not_implemented_or_skip(data, ndata, msg)
    #nelements = None
    #ntotal = None
    return n


def oes_cquad4_144_complex_77(op2: OP2,
                              data: bytes,
                              obj: ComplexPlateStressArray | ComplexPlateStrainArray,
                              nelements: int, nnodes: int,
                              dt: Any,
                              is_magnitude_phase: bool) -> int:
    #print('nelements =', nelements)
    n = 0
    grid_center = 0
    if op2.size == 4:
        s1 = op2.struct_2i  # 2
    else:
        s1 = op2.struct_2q  # 2
    s2 = Struct(op2._endian + mapfmt(b'i14f', op2.size)) # 15

    ntotal1 = 8 * op2.factor
    ntotal2 = 60 * op2.factor  # 4*15

    add_sort_x = getattr(obj, 'add_sort' + str(op2.sort_method))
    for unused_i in range(nelements):
        (eid_device, _) = s1.unpack(data[n:n+ntotal1])
        n += ntotal1

        edata = data[n:n+ntotal2]
        n += ntotal2
        out = s2.unpack(edata)  # len=15*4

        eid = op2._check_id(eid_device, 'element', out)
        if op2.is_debug_file:
            op2.binary_debug.write('%s\n' % (str(out)))
        (grid,
         fd1, sx1r, sx1i, sy1r, sy1i, txy1r, txy1i,
         fd2, sx2r, sx2i, sy2r, sy2i, txy2r, txy2i) = out
        #grid_center = 'CEN/%i' % grid   # this is correct, but fails

        if is_magnitude_phase:
            sx1 = polar_to_real_imag(sx1r, sx1i)
            sy1 = polar_to_real_imag(sy1r, sy1i)
            sx2 = polar_to_real_imag(sx2r, sx2i)
            sy2 = polar_to_real_imag(sy2r, sy2i)
            txy1 = polar_to_real_imag(txy1r, txy1i)
            txy2 = polar_to_real_imag(txy2r, txy2i)
        else:
            sx1 = complex(sx1r, sx1i)
            sy1 = complex(sy1r, sy1i)
            sx2 = complex(sx2r, sx2i)
            sy2 = complex(sy2r, sy2i)
            txy1 = complex(txy1r, txy1i)
            txy2 = complex(txy2r, txy2i)

        add_sort_x(dt, eid, grid_center,
                   fd1, sx1, sy1, txy1,
                   fd2, sx2, sy2, txy2)

        for unused_node_id in range(nnodes):  # nodes pts
            edata = data[n:n + ntotal2]
            n += ntotal2
            out = s2.unpack(edata)
            if op2.is_debug_file:
                op2.binary_debug.write('%s\n' % (str(out)))
            (grid,
             fd1, sx1r, sx1i, sy1r, sy1i, txy1r, txy1i,
             fd2, sx2r, sx2i, sy2r, sy2i, txy2r, txy2i) = out

            if is_magnitude_phase:
                sx1 = polar_to_real_imag(sx1r, sx1i)
                sx2 = polar_to_real_imag(sx2r, sx2i)
                sy1 = polar_to_real_imag(sy1r, sy1i)
                sy2 = polar_to_real_imag(sy2r, sy2i)
                txy1 = polar_to_real_imag(txy1r, txy1i)
                txy2 = polar_to_real_imag(txy2r, txy2i)
            else:
                sx1 = complex(sx1r, sx1i)
                sx2 = complex(sx2r, sx2i)
                sy1 = complex(sy1r, sy1i)
                sy2 = complex(sy2r, sy2i)
                txy1 = complex(txy1r, txy1i)
                txy2 = complex(txy2r, txy2i)

            add_sort_x(dt, eid, grid,
                       fd1, sx1, sy1, txy1,
                       fd2, sx2, sy2, txy2)
    return n


def oes_cquad4_complex_vm_87(op2: OP2, data: bytes,
                             obj: ComplexPlateVMStressArray | ComplexPlateVMStrainArray,
                             nelements: int, nnodes: int,
                             is_magnitude_phase: bool) -> int:
    """
    etype name      numwide nodes
    33  CQUAD4-33 -> 17      (not in this function)
    82  CQUADR    -> 87      (1+4 nodes)
    144 CQUAD144  -> 87      (1+4 nodes)
    ??? CTRIA6    -> 70?     (1+3 nodes)

    2 TERM    CHAR4
    3 GRID    I
    4 FD1     RS Fiber distance at z1
    5 EX1R    RS Normal in x at z1
    6 EX1I    RS Normal in x at z1
    7 EY1R    RS Normal in y at z1
    8 EY1I    RS Normal in y at z1
    9 ETXY1R  RS Shear in xy at z1
    10 ETXY1I RS Shear in xy at z1
    11 VM1    RS von Mises at z1
    12 FD2    RS Fiber distance at z2
    13 EX2R   RS Normal in x at z2
    14 EX2I   RS Normal in x at z2
    15 EY2R   RS Normal in y at z2
    16 EY2I   RS Normal in y at z2
    17 ETXY2R RS Shear in xy at z2
    18 ETXY2I RS Shear in xy at z2
    19 VM2    RS von Mises at z2

    """
    #strings = (b'20.0 CEN/ 4 \xcd\xccL\xbd\x91@\x1aD\xebX\xac\xc1b\x80FB\x96\xc9\xdd\xbf\xe7\xbe\n,\xaah\xc5\xa9b\x87\x14D\xcd\xccL=\x91@\x1a\xc4\xebX\xacAb\x80F\xc2\x96\xc9\xdd?\xe7\xbe\n\xac\xaah\xc5)b\x87\x14D\x01\x00\x00\x00\xcd\xccL\xbd\xc7\xd5\x9fD\xc9\x952\xc2b\x80FC\x96\xc9\xdd\xc0\xe7\xbe\n,\xaah\xc5\xa9\xca\x13\x95D\xcd\xccL=\xc7\xd5\x9f\xc4\xc9\x952Bb\x80F\xc3\x96\xc9\xdd@\xe7\xbe\n\xac\xaah\xc5)\xca\x13\x95D\x02\x00\x00\x00\xcd\xccL\xbd\xbf\xa62\xc2\xd4\x9b\xc7?b\x80\xc6\xc2\x96\xc9]@\xe7\xbe\n,\xaah\xc5\xa9nL\xacB\xcd\xccL=\xbf\xa62B\xd4\x9b\xc7\xbfb\x80\xc6B\x96\xc9]\xc0\xe7\xbe\n\xac\xaah\xc5)nL\xacB\x03\x00\x00\x00\xcd\xccL\xbd\xbf\xa62\xc2\xd4\x9b\xc7?b\x80\xc6\xc2\x96\xc9]@\xe7\xbe\n,\xaah\xc5\xa9nL\xacB\xcd\xccL=\xbf\xa62B\xd4\x9b\xc7\xbfb\x80\xc6B\x96\xc9]\xc0\xe7\xbe\n\xac\xaah\xc5)nL\xacB\x04\x00\x00\x00\xcd\xccL\xbd\xc7\xd5\x9fD\xc9\x952\xc2b\x80FC\x96\xc9\xdd\xc0\xe7\xbe\n,\xaah\xc5\xa9\xca\x13\x95D\xcd\xccL=\xc7\xd5\x9f\xc4\xc9\x952Bb\x80F\xc3\x96\xc9\xdd@\xe7\xbe\n\xac\xaah\xc5)\xca\x13\x95D',)
    # ints    = (20.0, CEN/, 4,
    #           -0.05, 617.0, -21.5, 49.6, -1.7, 738901735,              -1446680406,            1142196066,
    #           0.05, -617.0,  21.5, -49.6, 1.7, -1408581913, 700803242, 1142196066, 1, -1119040307, 1151325639, -1036872247, 1128693858, -1059206762, 738901735, -1446680406, 1150620618, 1028443341, -996158009, 1110611401, -1018789790, 1088276886, -1408581913, 700803242, 1150620618, 2, -1119040307, -1036867905, 1070046164, -1027178398, 1079888278, 738901735, -1446680406, 1118588014, 1028443341, 1110615743, -1077437484, 1120305250, -1067595370, -1408581913, 700803242, 1118588014, 3, -1119040307, -1036867905, 1070046164, -1027178398, 1079888278, 738901735, -1446680406, 1118588014, 1028443341, 1110615743, -1077437484, 1120305250, -1067595370, -1408581913, 700803242, 1118588014, 4, -1119040307, 1151325639, -1036872247, 1128693858, -1059206762, 738901735, -1446680406, 1150620618, 1028443341, -996158009, 1110611401, -1018789790, 1088276886, -1408581913, 700803242, 1150620618)
    #floats  = (20.0, CEN/, 4,
    #           -0.05,  617.0,   -21.5,   49.6,  -1.7,  1.9716951595721843e-12, -8.766713754677219e-14,  594.1,
    #            0.05, -617.0,    21.5,  -49.6,   1.7, -1.9716951595721843e-12,  8.766713754677219e-14,  594.1,
    #        1.4012984e-45,
    #           -0.05,  1278.68, -44.6,  198.5,  -6.9,  1.9716951595721843e-12, -8.766713754677219e-14, 1192.6,
    #            0.05, -1278.68,  44.6, -198.5,   6.9, -1.9716951595721843e-12,  8.766713754677219e-14, 1192.6,
    #        2.802596928649634e-45,
    #           -0.05, -44.6,      1.5,  -99.2,   3.4,  1.9716951595721843e-12, -8.766713754677219e-14, 86.149,
    #            0.05,  44.6,     -1.5,   99.2,  -3.4, -1.9716951595721843e-12,  8.766713754677219e-14, 86.149,
    #        4.203895392974451e-45,
    #           -0.05, -44.6,      1.5,  -99.2,   3.4,  1.9716951595721843e-12, -8.766713754677219e-14, 86.149,
    #            0.05,  44.6,     -1.5,   99.2,  -3.4, -1.9716951595721843e-12,  8.766713754677219e-14, 86.149,
    #        5.605193857299268e-45,
    #           -0.05, 1278.6,    -44.6, 198.5,  -6.9,  1.9716951595721843e-12, -8.766713754677219e-14, 1192.6,
    #           0.05, -1278.6,     44.6,-198.5,   6.9, -1.9716951595721843e-12,  8.766713754677219e-14, 1192.6)
    # -------------------------
    #ints    = (0,    'CEN/',
    #           5,     -0.025, 1131762990, -1056779073, 1137493831, -1051039598, 1133740074, -1055043205, 1142680183,
    #                   0.025, -1016740689, 1089802185, -1009242998, 1097240578, -1013972178, 1092196624, 1142689762,
    #           80000, -0.025, -1039337422, 1066818106, -1022560206, 1083595322, 1118381010, -1070307877, 1128400256,
    #                   0.025, 1111837962, -1076728084, 1128615178, -1059950868, -1029871408, 1076355832, 1130378350,
    #           80002, -0.025, 1144778127, -1043828805, 1161555343, -1027051589, 1134339017, -1054404396, 1160567412,
    #                   0.025, -1002241635, 1104149606, -985464419, 1120926822, -1014451576, 1091685317, 1160927605,
    #           80202, -0.025, -1013123432, 1093101862, -996473336, 1109743497, 1143660503, -1045020817, 1154190671,
    #                   0.025, 1133381240, -1055425923, 1150079339, -1038733090, -1002862704, 1103487198, 1154274046,
    #           80200, -0.025, 1141944942, -1046850567, 1095277553, -1093274063, 1120227570, -1068338410, 1142258188,
    #                   0.025, -1007929795, 1098641187, 1116060214, -1072783145, -1033402643, 1073147904, 1140969037,
    #           2.5,  'CEN/',
    #           5,     -0.025, 1131832225, -1056704978, 1137555116, -1050908425, 1133787322, -1054942075, 1142729036,
    #                   0.025, -1016671483, 1089950313, -1009181671, 1097371840, -1013924936, 1092297741, 1142738537,
    #           80000, -0.025, -1039289241, 1066921230, -1022512025, 1083698446, 1118448079, -1070164312, 1128472712,
    #                   0.025, 1111886248, -1076624736, 1128663464, -1059847520, -1029804359, 1076499355, 1130448349,
    #           80002, -0.025, 1144846119, -1043683269, 1161623335, -1026906053, 1134386659, -1054302422, 1160629636,
    #                   0.025, -1002173624, 1104295182, -985396408, 1121072398, -1014403979, 1091787195, 1160989873,
    #           80202, -0.025, -1013031023, 1093299683, -996403997, 1109891922, 1143714425, -1044905405, 1154266576,
    #                   0.025, 1133473626, -1055228152, 1150148641, -1038584742, -1002808755, 1103602670, 1154348580,
    #           80200, -0.025, 1141995401, -1046742558, 1095362185, -1093092785, 1120294524, -1068195089, 1142310151,
    #                   0.025, -1007828983, 1098856980, 1116049937, -1072805159, -1033268970, 1073434045, 1141020308,
    #           10.0, 'CEN/',
    #           5,     -0.025, 1132706983, -1055459734, 1138546383, -1048703317, 1134552107, -1053240134, 1143519784,
    #                   0.025, -1015552254, 1091479459, -1008189778, 1099242966, -1013160238, 1093999495, 1143528113,
    #           80000, -0.025, -1038510246, 1068653763, -1021733030, 1085430979, 1119535604, -1067741982, 1129646715,
    #                   0.025, 1112666820, -1074888834, 1129444036, -1058111618, -1028717126, 1078921059, 1131585244,
    #           80002, -0.025,  1145947464, -1041231474, 1162724680, -1024454258, 1135157762, -1052586493, 1161637710,
    #                   0.025, -1001071999, 1106747574, -984294783, 1123524790, -1013633545, 1093501695, 1161998610,
    #           80202, -0.025, -1011529318,  1096648297, -995279005, 1112398467, 1144586718, -1042964811, 1155498288,
    #                   0.025, 1134974980, -1051880286, 1151273087, -1036079365, -1001936046, 1105544152, 1155559868,
    #           80200, -0.025,  1142813256, -1044921274, 1096762940, -1090229644, 1121380799, -1065774838, 1143152359,
    #                   0.025, -1006413905, 1100701913, 1115879378, -1073189541, -1031449359, 1076004431, 1141852280)
    #floats  = (0.0,  'CEN/',
    #           5,       -0.025,  245.33273315429688, -8.176939964294434,   409.5568542480469, -13.650529861450195,  295.00128173828125,  -9.832392692565918,
    #           623.8668, 0.025, -229.76829528808594,  7.658176898956299, -432.34796142578125,  14.410158157348633, -288.02484130859375,   9.599868774414062, 624.2481689453125,
    #           1.12-40, -0.025, -35.24237823486328,   1.174628496170044, -140.96951293945312,   4.698513984680176,   84.56996154785156, -2.8187167644500732,
    #           194.4375, 0.025,  49.325233459472656, -1.644010066986084,  197.30093383789062,  -6.576040267944336,   -78.7047119140625,   2.623228073120117, 224.20480346679688,
    #           1.12106, -0.025,  751.7118530273438, -25.054555892944336,   3006.847412109375, -100.21822357177734,   313.2795715332031, -10.441608428955078,
    #           2765.65,  0.025, -780.0252075195312, 25.998241424560547, -3120.100830078125, 103.99296569824219, -273.394775390625, 9.112248420715332, 2853.591064453125,
    #           1.1238,  -0.025, -313.926513671875, 10.463171005249023, -1240.1884765625, 41.33548355102539, 683.4974975585938, -22.780973434448242,
    #           1628.4,   0.025, 284.050537109375, -9.46740436553955, 1126.5443115234375, -37.54772186279297, -742.1181640625, 24.734798431396484, 1638.593505859375,
    #           1.12,    -0.025, 578.7879638671875, -19.29100227355957, 12.538071632385254, -0.4178939163684845, 98.65809631347656, -3.2882742881774902,
    #           597.90,   0.025, -472.4237365722656, 15.745882987976074, 66.86369323730469, -2.228566884994507, -57.88176345825195, 1.92919921875, 519.2234497070312,
    #           2.5,  'CEN/',
    #           5,          -0.025, 246.38917541503906, -8.247602462768555, 411.4271240234375, -13.775626182556152, 296.44317626953125, -9.928837776184082,
    #           626.671875,  0.025, -230.8242950439453, 7.728809833526611, -434.2195129394531, 14.53533935546875, -289.466552734375, 9.696301460266113, 627.2251586914062,
    #           1.1217e-40, -0.025, -35.42617416381836, 1.1869218349456787, -141.70469665527344, 4.747687339782715, 85.08165740966797, -2.852945327758789,
    #           195.153125,  0.025, 49.509429931640625, -1.6563301086425781, 198.0377197265625, -6.6253204345703125, -79.21625518798828, 2.6574466228485107, 225.2729034423828,
    #           1.1212e-40, -0.025, 755.8617553710938, -25.332143783569336, 3023.447021484375, -101.32857513427734, 314.7334899902344, -10.538858413696289,
    #           2780.65625,  0.025, -784.17626953125, 26.27590560913086, -3136.705078125, 105.10362243652344, -274.8473205566406, 9.209406852722168, 2868.793212890625,
    #           1.1358e-40, -0.025, -316.7466125488281, 10.651827812194824, -1248.6527099609375, 41.90167999267578, 686.7886352539062, -23.0011043548584,
    #           1637.40625,  0.025, 286.86993408203125, -9.656013488769531, 1135.0040283203125, -38.113624572753906, -745.4109497070312, 24.95504379272461, 1647.69189453125,
    #           1.1383e-40, -0.025, 581.8677368164062, -19.497013092041016, 12.618782997131348, -0.42329642176628113, 99.16891479492188, -3.3224446773529053,
    #           601.471875,  0.025, -475.5002746582031, 15.951679229736328, 66.78528594970703, -2.223318338394165, -58.391685485839844, 1.96330988407135, 522.352783203125,
    #           10.0, 'CEN/',
    #           5,          -0.025, 263.4738464355469, -9.435159683227539, 441.6781921386719, -15.878581047058105, 319.7825622558594, -11.551935195922852,
    #           674.140625,  0.025, -247.90237426757812, 8.915926933288574, -464.48968505859375, 16.639568328857422, -312.80328369140625, 11.319220542907715, 675.4170532226562,
    #           1.1337e-40, -0.025, -38.397804260253906, 1.3934558629989624, -153.9121704101562, 5.57382345199585, 93.37881469726562, -3.430473804473877,
    #           213.310547,  0.025, 52.48707580566406, -1.8632657527923584, 209.94830322265625, -7.453063011169434, -87.51118469238281, 3.234825849533081, 242.62054443359375,
    #           1.1202e-40, -0.025, 823.08251953125, -30.00857162475586, 3292.330078125, -120.03428649902344, 338.26568603515625, -12.17529582977295,
    #           3026296875,  0.025, -851.4141235351562, 30.953472137451172, -3405.656494140625, 123.81388854980469, -298.3591003417969, 10.844481468200684, 3115.06689453125,
    #           1.1298e-40, -0.025, -362.57501220703125, 13.845314979553223, -1385.9808349609375, 51.4633903503418, 740.0291748046875, -26.70249366760254,
    #           1788109375,  0.025, 332.6876220703125, -12.848787307739258, 1272.2655029296875, -47.67087936401367, -798.6768798828125, 28.658126831054688, 1795.55419921875,
    #           1.1293e-40, -0.025, 631.78564453125, -22.970836639404297, 13.954647064208984, -0.517249345779419, 107.45653533935547, -3.899477481842041,
    #           652.410938,  0.025, -525.3700561523438, 19.42228889465332, 65.48402404785156, -2.1316745281219482, -66.66590118408203, 2.539447546005249, 573.13232421875)
    n = 0
    ntotal1 = 8 * op2.factor  # 2*4
    ntotal2 = 68 * op2.factor # 17*4
    add_sort_x = getattr(obj, 'add_sort' + str(op2.sort_method))
    #print('ndata', len(data))
    #print('ntotal', ntotal1+nnodes*ntotal2)
    #print('nelements =', nelements)
    for i in range(nelements):
        #print(f'{i}/{nelements}')
        edata1 = data[n:n+ntotal1]
        #self.show_data(edata1)
        n += ntotal1

        fmt1 = mapfmt(op2._endian + op2._analysis_code_fmt + b'4s', op2.size)
        #fmt2 = mapfmt(op2._endian + b'i 5f2if 5f2if', op2.size)
        fmt2 = mapfmt(op2._endian + b'i 16f', op2.size)
        s1 = Struct(fmt1)
        s2 = Struct(fmt2)

        eid_device, cen = s1.unpack(edata1)
        #self.show_data(data, types='ifs')
        assert cen == b'CEN/', (eid_device, cen)

        eid, dt = get_eid_dt_from_eid_device(
            eid_device, op2.nonlinear_factor, op2.sort_method)

        for inode in range(nnodes):
            edata2 = data[n:n+ntotal2]

            #DISTANCE             NORMAL-X                     NORMAL-Y                    SHEAR-XY             VON MISES
            out = s2.unpack(edata2)
            (nid,
             fd1, oxx1r, oxx1i, oyy1r, oyy1i, txy1r, txy1i, ovm1,
             fd2, oxx2r, oxx2i, oyy2r, oyy2i, txy2r, txy2i, ovm2) = out
            fd1, oxx1r, oxx1i, oyy1r, oyy1i, txy1r, txy1i, ovm1,
            if inode == 0:
                nid = 0

            if is_magnitude_phase:
                oxx1 = polar_to_real_imag(oxx1r, oxx1i)
                oxx2 = polar_to_real_imag(oxx2r, oxx2i)
                oyy1 = polar_to_real_imag(oyy1r, oyy1i)
                oyy2 = polar_to_real_imag(oyy2r, oyy2i)
                txy1 = polar_to_real_imag(txy1r, txy1i)
                txy2 = polar_to_real_imag(txy2r, txy2i)
            else:
                oxx1 = complex(oxx1r, oxx1i)
                oxx2 = complex(oxx2r, oxx2i)
                oyy1 = complex(oyy1r, oyy1i)
                oyy2 = complex(oyy2r, oyy2i)
                txy1 = complex(txy1r, txy1i)
                txy2 = complex(txy2r, txy2i)

            #print((eid, dt), nid, f'{fd1:g}, {fd2:g}')
            #print('real', f'{oxx1r:g}, {oxx1i:g}, {oyy1r:g}, {oyy1i:g}, {txy1r:g}, {txy1i:g}, {ovm1:g}')
            #print('imag', oxx2i, oxx2i, oyy2r, oyy2i, txy2r, txy2i, ovm2)
            add_sort_x(dt, eid, nid,
                       fd1, oxx1, oyy1, txy1, ovm1,
                       fd2, oxx2, oyy2, txy2, ovm2)
            #print(out)
            #assert nid in (1, 2, 3, 4, 17), ((eid, dt), nid, fd1, fd2)
            n += ntotal2
        #print('-----------')
    return n


def oes_cquad4_33_random_9(op2: OP2, data: bytes,
                           obj: RandomPlateStressArray | RandomPlateStrainArray,
                           nelements: int, ntotal: int) -> int:
    n = 0
    add_sort_x = getattr(obj, 'add_sort' + str(op2.sort_method))
    #print('cquad33_9 - SORT1')
    struct1 = Struct(op2._endian + op2._analysis_code_fmt + b'8f')
    #cen = 0 # CEN/4
    for unused_i in range(nelements):
        edata = data[n:n+ntotal]
        out = struct1.unpack(edata)

        (eid_device,
         fd1, sx1, sy1, txy1,
         fd2, sx2, sy2, txy2,) = out
        eid, dt = get_eid_dt_from_eid_device(
            eid_device, op2.nonlinear_factor, op2.sort_method)
        #self.log.info(f'SORT1 no VM: eid={eid} dt={dt}')

        if op2.is_debug_file:
            op2.binary_debug.write('  eid=%i C=[%s]\n' % (
                eid, ', '.join(['%r' % di for di in out])))
        add_sort_x(dt, eid, 0,
                   fd1, sx1, sy1, txy1,
                   fd2, sx2, sy2, txy2)
        n += ntotal
    return n


def oes_ctria3_real_17(op2: OP2, data: bytes,
                       obj: RealPlateStressArray | RealPlateStrainArray,
                       ntotal: int, nelements: int, dt: Any) -> int:
    n = 0
    cen = 0  # 'CEN/3'
    #assert op2.sort_method == 1, op2.code_information()

    add_new_eid_sort_x = getattr(obj, 'add_new_eid_sort' + str(op2.sort_method))
    #print(add_new_eid_sort_x)
    struct1 = Struct(op2._endian + mapfmt(op2._analysis_code_fmt + b'16f', op2.size))
    #print('nelements =', nelements, nelements*2, obj.ntotal)
    for unused_i in range(nelements):
        edata = data[n:n + ntotal]
        out = struct1.unpack(edata)
        (eid_device,
         fd1, sx1, sy1, txy1, angle1, major1, minor1, vm1,
         fd2, sx2, sy2, txy2, angle2, major2, minor2, vm2,) = out
        eid, dt = get_eid_dt_from_eid_device(
            eid_device, op2.nonlinear_factor, op2.sort_method)
        if op2.is_debug_file:
            op2.binary_debug.write('  OES CTRIA3-74 - eid=%i; C=[%s]\n' % (
                eid, ', '.join(['%r' % di for di in out])))
        #print('  OES CTRIA3-74 - eid=%i; C=[%s]\n' % (
            #eid, ', '.join(['%r' % di for di in out])))

        add_new_eid_sort_x(dt, eid, cen,
                           fd1, sx1, sy1, txy1, angle1, major1, minor1, vm1,
                           fd2, sx2, sy2, txy2, angle2, major2, minor2, vm2)
        n += ntotal
    return n


def oes_quad4_33_real_17(op2: OP2, data: bytes,
                         obj: RealPlateStressArray | RealPlateStrainArray,
                         ntotal: int, nelements: int, dt: Any) -> int:
    n = 0
    cen = 0 # CEN/4
    #assert op2.sort_method == 1, op2.code_information()
    add_new_eid_sort_x = getattr(obj, 'add_new_eid_sort' + str(op2.sort_method))
    fmt = mapfmt(op2._endian + op2._analysis_code_fmt + b'16f', op2.size)
    struct1 = Struct(fmt)
    #if op2.sort_method == 1:
        #print('**nelements =', nelements, nelements*2, obj.ntotal)
    #else:
        #print('**ntimes =', nelements)
    for unused_i in range(nelements):
        edata = data[n:n+ntotal]
        out = struct1.unpack(edata)

        (eid_device,
         fd1, sx1, sy1, txy1, angle1, major1, minor1, max_shear1,
         fd2, sx2, sy2, txy2, angle2, major2, minor2, max_shear2) = out

        eid, dt = get_eid_dt_from_eid_device(
            eid_device, op2.nonlinear_factor, op2.sort_method)
        if op2.is_debug_file:
            op2.binary_debug.write('  eid=%i C=[%s]\n' % (
                eid, ', '.join(['%r' % di for di in out])))

        add_new_eid_sort_x(
            dt, eid, cen,
            fd1, sx1, sy1, txy1, angle1, major1, minor1, max_shear1,
            fd2, sx2, sy2, txy2, angle2, major2, minor2, max_shear2)
        n += ntotal
    return n

def oes_ctria3_complex_vm_17(op2: OP2,
                             data: bytes,
                             obj: ComplexPlateVMStressArray | ComplexPlateVMStrainArray,
                             nelements: int, ntotal: int,
                             dt,
                             is_magnitude_phase: bool) -> int:
    """
    #                   C O M P L E X   S T R E S S E S   I N   T R I A N G U L A R   E L E M E N T S   ( T R I A 3 )
    #                                                          (REAL/IMAGINARY)
    #
    #  ELEMENT       FIBER                                     - STRESSES IN ELEMENT  COORDINATE SYSTEM -
    #    ID.        DISTANCE              NORMAL-X                       NORMAL-Y                      SHEAR-XY               VON MISES
    #0       1  -4.359080E+00  -1.391918E+00 /  2.474756E-03  -1.423926E+00 /  2.530494E-03   2.655153E-02 / -5.158625E-05   1.408948E+00
    #            4.359080E+00   1.391918E+00 / -2.474756E-03   1.423926E+00 / -2.530494E-03  -2.655153E-02 /  5.158625E-05   1.408948E+00

    """
    n = 0
    cen = 0
    struct1 = Struct(mapfmt(op2._endian + op2._analysis_code_fmt + b'16f', op2.size))
    add_sort_x = getattr(obj, 'add_sort' + str(op2.sort_method))

    for unused_i in range(nelements):
        edata = data[n:n + ntotal]
        out = struct1.unpack(edata)
        (eid_device, fd1, oxx1r, oxx1i, oyy1r, oyy1i, txy1r, txy1i, ovm1,
                     fd2, oxx2r, oxx2i, oyy2r, oyy2i, txy2r, txy2i, ovm2, ) = out

        eid, dt = get_eid_dt_from_eid_device(
            eid_device, op2.nonlinear_factor, op2.sort_method)

        if is_magnitude_phase:
            oxx1 = polar_to_real_imag(oxx1r, oxx1i)
            oxx2 = polar_to_real_imag(oxx2r, oxx2i)
            oyy1 = polar_to_real_imag(oyy1r, oyy1i)
            oyy2 = polar_to_real_imag(oyy2r, oyy2i)
            txy1 = polar_to_real_imag(txy1r, txy1i)
            txy2 = polar_to_real_imag(txy2r, txy2i)
        else:
            oxx1 = complex(oxx1r, oxx1i)
            oxx2 = complex(oxx2r, oxx2i)
            oyy1 = complex(oyy1r, oyy1i)
            oyy2 = complex(oyy2r, oyy2i)
            txy1 = complex(txy1r, txy1i)
            txy2 = complex(txy2r, txy2i)
        #print(dt, eid, cen, sx1, sy1, txy1, max_shear1)
        #print(dt, eid, cen, sx2, sy2, txy2, max_shear2)
        add_sort_x(dt, eid, cen,
                   fd1, oxx1, oyy1, txy1, ovm1,
                   fd2, oxx2, oyy2, txy2, ovm2)
        n += ntotal

    #msg = op2.code_information()
    #self.log.warning(f'  skipping {op2.table_name_str} {op2.element_name}')
    #msg = '%s-%s' % (op2.table_name_str, op2.element_name)
    return n


def oes_ctria3_random_9(op2: OP2, data: bytes,
                        obj: RandomPlateStressArray | RandomPlateStrainArray,
                        nelements: int, ntotal: int) -> int:
    n = 0
    struct1 = Struct(op2._endian + op2._analysis_code_fmt + b'8f')
    cen = 0 # CEN/4
    add_sort_x = getattr(obj, 'add_sort' + str(op2.sort_method))
    #print('ctria3_9 - SORT1')
    for unused_i in range(nelements):
        edata = data[n:n+ntotal]
        out = struct1.unpack(edata)

        (eid_device,
         fd1, sx1, sy1, txy1,
         fd2, sx2, sy2, txy2,) = out

        eid, dt = get_eid_dt_from_eid_device(
            eid_device, op2.nonlinear_factor, op2.sort_method)

        if op2.is_debug_file:
            op2.binary_debug.write('  eid=%i C=[%s]\n' % (
                eid, ', '.join(['%r' % di for di in out])))

        add_sort_x(dt, eid, cen,
                   fd1, sx1, sy1, txy1,
                   fd2, sx2, sy2, txy2)
        n += ntotal
    return n


def oes_cquad4_33_random_vm_11(op2: OP2, data: bytes,
                               obj: RandomPlateVMStressArray | RandomPlateVMStrainArray,
                               nelements: int, ntotal: int) -> int:
    struct1 = Struct(op2._endian + op2._analysis_code_fmt + b'10f')
    cen = 0 # CEN/4
    n = 0
    add_sort_x = getattr(obj, 'add_sort' + str(op2.sort_method))
    for unused_i in range(nelements):
        edata = data[n:n+ntotal]
        out = struct1.unpack(edata)

        (eid_device,
         fd1, sx1, sy1, txy1, ovm1,
         fd2, sx2, sy2, txy2, ovm2,) = out

        eid, dt = get_eid_dt_from_eid_device(
            eid_device, op2.nonlinear_factor, op2.sort_method)
        #self.log.info(f'SORT1 VM: eid={eid} dt={dt}')

        if op2.is_debug_file:
            op2.binary_debug.write('  eid=%i C=[%s]\n' % (
                eid, ', '.join(['%r' % di for di in out])))

        add_sort_x(dt, eid, cen,
                   fd1, sx1, sy1, txy1, ovm1,
                   fd2, sx2, sy2, txy2, ovm2)
        n += ntotal
    return n


def oes_ctria3_random_vm_11(op2: OP2, data: bytes,
                            obj: RandomPlateStressArray | RandomPlateStrainArray,
                            nelements: int, ntotal: int) -> int:
    n = 0
    struct1 = Struct(op2._endian + op2._analysis_code_fmt + b'10f')

    add_sort_x = getattr(obj, 'add_sort' + str(op2.sort_method))
    cen = 0 # CEN/4
    for unused_i in range(nelements):
        edata = data[n:n+ntotal]
        out = struct1.unpack(edata)

        (eid_device,
         fd1, sx1, sy1, txy1, ovm1,
         fd2, sx2, sy2, txy2, ovm2,) = out
        #print('CTRIA3', op2.nonlinear_factor, out)

        eid, dt = get_eid_dt_from_eid_device(
            eid_device, op2.nonlinear_factor, op2.sort_method)
        #self.log.info(f'SORT1 VM: eid={eid} dt={dt}')

        if op2.is_debug_file:
            op2.binary_debug.write('  eid=%i C=[%s]\n' % (
                eid, ', '.join(['%r' % di for di in out])))

        add_sort_x(dt, eid, cen,
                  fd1, sx1, sy1, txy1, ovm1,
                  fd2, sx2, sy2, txy2, ovm2)
        n += ntotal
    return n


def oes_cquad4_144_real(op2: OP2, data: bytes, ndata: int,
                        obj: RealPlateStrainArray,
                        nelements: int, nnodes: int, dt: Any) -> int:
    n = 0
    if op2.size == 4:
        center_format = op2._endian + op2._analysis_code_fmt + b'4si16f'
        node_format = op2._endian + b'i16f'
    else:
        center_format = op2._endian + mapfmt(op2._analysis_code_fmt, op2.size) + b'8sq16d'
        node_format = op2._endian + b'q16d'
    cs = Struct(center_format)
    ns = Struct(node_format)

    add_new_eid_sort_x = getattr(obj, 'add_new_eid_sort' + str(op2.sort_method))
    add_sort_x = getattr(obj, 'add_sort' + str(op2.sort_method))
    if op2.is_debug_file:
        op2.binary_debug.write(
            '  [cap, element1, element2, ..., cap]\n'
            '  cap = %i  # assume 1 cap when there could have been multiple\n'
            '  #elementi = [centeri, node1i, node2i, node3i, node4i]\n'
            '  #centeri = [eid_device, j, grid, fd1, sx1, sy1, txy1, angle1, major1, minor1, vm1,\n'
            '  #                                fd2, sx2, sy2, txy2, angle2, major2, minor2, vm2,)]\n'
            '  #nodeji = [grid, fd1, sx1, sy1, txy1, angle1, major1, minor1, vm1,\n'
            '  #                fd2, sx2, sy2, txy2, angle2, major2, minor2, vm2,)]\n'
            '  nelements=%i; nnodes=%i # +1 centroid\n' % (ndata, nelements, nnodes))

    grid_center = 0
    n76 = 76 * op2.factor
    n68 = 68 * op2.factor
    for unused_i in range(nelements):
        edata = data[n:n+n76]

        out = cs.unpack(edata)  # len=17*4
        # j is the number of nodes, so CQUAD4 -> 4, but we don't need to save it...
        (eid_device, unused_j,
         grid,
         fd1, sx1, sy1, txy1, angle1, major1, minor1, vm1,
         fd2, sx2, sy2, txy2, angle2, major2, minor2, vm2,) = out
        eid, dt = get_eid_dt_from_eid_device(
            eid_device, op2.nonlinear_factor, op2.sort_method)
        #print(f'tri eid={eid} dt={dt} nid={grid_center} unused_j={unused_j} grid={grid}')
        #print(out[:3])
        if op2.is_debug_file:
            op2.binary_debug.write('  eid=%i; C=[%s]\n' % (eid, ', '.join(['%r' % di for di in out])))
        add_new_eid_sort_x(dt, eid, grid_center,
                           fd1, sx1, sy1, txy1, angle1, major1, minor1, vm1,
                           fd2, sx2, sy2, txy2, angle2, major2, minor2, vm2)
        n += n76
        for inode in range(nnodes):
            out = ns.unpack(data[n:n + n68])
            (grid,
             fd1, sx1, sy1, txy1, angle1, major1, minor1, vm1,
             fd2, sx2, sy2, txy2, angle2, major2, minor2, vm2,) = out

            if op2.is_debug_file:
                d = tuple([grid,
                           fd1, sx1, sy1, txy1, angle1, major1, minor1, vm1,
                           fd2, sx2, sy2, txy2, angle2, major2, minor2, vm2])
                op2.binary_debug.write('  node%i = [%s]\n' % (inode+1, ', '.join(['%r' % di for di in d])))
            assert isinstance(grid, int), grid
            assert grid > 0, grid

            #print(f' {inode}: eid={eid} dt={dt} nid={grid}')
            add_sort_x(dt, eid, grid,
                       fd1, sx1, sy1, txy1, angle1, major1, minor1, vm1,
                       fd2, sx2, sy2, txy2, angle2, major2, minor2, vm2)
            n += n68
    return n


def oes_cquad4_144_random(op2: OP2, data: bytes,
                          obj: RandomPlateStressArray | RandomPlateStrainArray,
                          nelements: int, nnodes: int, ndata: int) -> int:
    n = 0
    center_format = op2._endian + op2._analysis_code_fmt + b'4s i8f'
    node_format = op2._endian + b'i8f'
    cs = Struct(center_format)
    ns = Struct(node_format)

    if op2.is_debug_file:
        op2.binary_debug.write(
            '  [cap, element1, element2, ..., cap]\n'
            '  cap = %i  # assume 1 cap when there could have been multiple\n'
            '  #elementi = [centeri, node1i, node2i, node3i, node4i]\n'
            '  #centeri = [eid_device, j, grid, fd1, sx1, sy1, txy1,\n'
            '  #                                fd2, sx2, sy2, txy2,)]\n'
            '  #nodeji = [grid, fd1, sx1, sy1, txy1,\n'
            '  #                fd2, sx2, sy2, txy2,)]\n'
            '  nelements=%i; nnodes=%i # +1 centroid\n' % (ndata, nelements, nnodes))

    grid_center = 0
    add_sort_x = getattr(obj, 'add_sort' + str(op2.sort_method))

    ntotal1 = 44 * op2.factor
    ntotal2 = 36 * op2.factor
    for unused_i in range(nelements):
        edata = data[n:n+ntotal1]
        #self.show_data(edata)
        out = cs.unpack(edata)  # len=17*4
        #print(out)
        # j is the number of nodes, so CQUAD4 -> 4, but we don't need to save it...
        eid_device = out[0]
        (eid_device, unused_j,
         grid,
         fd1, sx1, sy1, txy1,
         fd2, sx2, sy2, txy2,) = out
        #print(out)

        eid, dt = get_eid_dt_from_eid_device(
            eid_device, op2.nonlinear_factor, op2.sort_method)

        if op2.is_debug_file:
            op2.binary_debug.write('  eid=%i; C=[%s]\n' % (eid, ', '.join(['%r' % di for di in out])))

        add_sort_x(dt, eid, grid_center,
                   fd1, sx1, sy1, txy1, fd2, sx2, sy2, txy2)
        n += ntotal1
        for inode in range(nnodes):
            out = ns.unpack(data[n:n + ntotal2])
            #print(out)
            (grid,
             fd1, sx1, sy1, txy1,
             fd2, sx2, sy2, txy2,) = out
            #print(out)
            if op2.is_debug_file:
                d = tuple([grid,
                           fd1, sx1, sy1, txy1,
                           fd2, sx2, sy2, txy2])
                op2.binary_debug.write('  node%i = [%s]\n' % (inode+1, ', '.join(['%r' % di for di in d])))
            assert isinstance(grid, int), grid
            assert grid > 0, grid

            # leaving off grid
            add_sort_x(dt, eid, grid,
                       fd1, sx1, sy1, txy1,
                       fd2, sx2, sy2, txy2)
            n += ntotal2
    return n


def oes_ctria3_complex_15(op2: OP2, data: bytes,
                          obj: ComplexPlateVMStressArray | ComplexPlateVMStrainArray,
                          nelements: int, ntotal: int,
                          is_magnitude_phase: bool,
                          element_name_type: str) -> int:
    struct1 = Struct(op2._endian + mapfmt(op2._analysis_code_fmt + b'14f', op2.size))
    cen = 0  # CEN/3
    n = 0
    for unused_i in range(nelements):
        edata = data[n:n + ntotal]
        out = struct1.unpack(edata)
        (eid_device,
         fd1, sx1r, sx1i, sy1r, sy1i, txy1r, txy1i,
         fd2, sx2r, sx2i, sy2r, sy2i, txy2r, txy2i,) = out
        eid, dt = get_eid_dt_from_eid_device(
            eid_device, op2.nonlinear_factor, op2.sort_method)

        if op2.is_debug_file:
            op2.binary_debug.write('  OESC %s - eid=%i; C=[%s]\n' % (
                element_name_type, eid,
                ', '.join(['%r' % di for di in out])))

        if is_magnitude_phase:
            sx1 = polar_to_real_imag(sx1r, sx1i)
            sy1 = polar_to_real_imag(sy1r, sy1i)
            sx2 = polar_to_real_imag(sx2r, sx2i)
            sy2 = polar_to_real_imag(sy2r, sy2i)
            txy1 = polar_to_real_imag(txy1r, txy1i)
            txy2 = polar_to_real_imag(txy2r, txy2i)
        else:
            sx1 = complex(sx1r, sx1i)
            sy1 = complex(sy1r, sy1i)
            sx2 = complex(sx2r, sx2i)
            sy2 = complex(sy2r, sy2i)
            txy1 = complex(txy1r, txy1i)
            txy2 = complex(txy2r, txy2i)
        obj.add_sort1(dt, eid, cen,
                      fd1, sx1, sy1, txy1,
                      fd2, sx2, sy2, txy2)
        n += ntotal
    return n
