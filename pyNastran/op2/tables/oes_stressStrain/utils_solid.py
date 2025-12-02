from __future__ import annotations
from struct import Struct
from typing import Any, TYPE_CHECKING

import numpy as np

from pyNastran.op2.op2_interface.op2_reader import mapfmt
from pyNastran.op2.op2_interface.utils import (
    apply_mag_phase,
)
from pyNastran.op2.tables.utils import get_is_slot_saved, get_eid_dt_from_eid_device
from pyNastran.op2.op2_helper import polar_to_real_imag

from pyNastran.op2.tables.oes_stressStrain.real.oes_solids import RealSolidStrainArray, RealSolidStressArray
from pyNastran.op2.tables.oes_stressStrain.complex.oes_solids import ComplexSolidStressArray, ComplexSolidStrainArray
from pyNastran.op2.tables.oes_stressStrain.complex.oes_solids_vm import ComplexSolidStressVMArray, ComplexSolidStrainVMArray
from pyNastran.op2.tables.oes_stressStrain.random.oes_solids import RandomSolidStressArray, RandomSolidStrainArray
from pyNastran.op2.tables.oes_stressStrain.random.oes_solids_vm import RandomSolidVMStressArray, RandomSolidVMStrainArray


if TYPE_CHECKING:  # pragma: no cover
    from pyNastran.op2.op2 import OP2
    from pyNastran.op2.tables.oes_stressStrain.oes import OES

def oes_csolid(oes: OES,
               op2: OP2, data, ndata: int, dt, is_magnitude_phase: bool,
               result_type: str, prefix: str, postfix: str):
    """
    reads stress/strain for element type:
     - 39 : CTETRA
     - 67 : CHEXA
     - 68 : CPENTA
     - 255: CPYRAM

    """
    factor = op2.factor
    size = op2.size
    n = 0
    etype_map = {
        # element_type : (element_base, nnodes_expected, element_name)
        39: ('ctetra', 5, 'CTETRA4'),
        67: ('chexa', 9, 'CHEXA8'),
        68: ('cpenta', 7, 'CPENTA6'),
        255: ('cpyram', 6, 'CPYRAM5'),
    }
    element_base, nnodes_expected, element_name = etype_map[op2.element_type]

    if op2.is_stress:
        stress_strain = 'stress'
        obj_vector_real = RealSolidStressArray
        obj_vector_complex = ComplexSolidStressArray
        obj_vector_complex_vm = ComplexSolidStressVMArray
        obj_vector_random = RandomSolidStressArray
        obj_vector_random_vm = RandomSolidVMStressArray
    else:
        stress_strain = 'strain'
        obj_vector_real = RealSolidStrainArray
        obj_vector_complex = ComplexSolidStrainArray
        obj_vector_complex_vm = ComplexSolidStrainVMArray
        obj_vector_random = RandomSolidStrainArray
        obj_vector_random_vm = RandomSolidVMStrainArray

    if prefix == '' and postfix == '':
        prefix = stress_strain + '.'

    # stress.chexa_stress
    result_name = prefix + f'{element_base}_{stress_strain}' + postfix

    is_saved, slot = get_is_slot_saved(op2, result_name)
    if not is_saved:
        return ndata, None, None

    numwide_real = 4 + 21 * nnodes_expected
    numwide_imag = 4 + (17 - 4) * nnodes_expected
    numwide_imag2 = 4 + 14 * nnodes_expected
    numwide_random = 4 + (11 - 4) * nnodes_expected
    numwide_random2_vm = 4 + 14 * nnodes_expected
    numwide_random3 = 4 + 8 * nnodes_expected
    preline1 = '%s-%s' % (op2.element_name, op2.element_type)
    preline2 = ' ' * len(preline1)

    # print('numwide real=%s imag=%s random=%s' % (numwide_real, numwide_imag, numwide_random2))
    op2._data_factor = nnodes_expected
    log = op2.log
    # print(op2.format_code, op2.table_name, op2.num_wide, op2.element_name)
    # print(op2.code_information())
    if op2.format_code == 1 and op2.num_wide == numwide_real:  # real
        op2.log.debug(f'numwide_real={numwide_real}')
        ntotal = (16 + 84 * nnodes_expected) * factor
        nelements = ndata // ntotal
        auto_return, is_vectorized = op2._create_oes_object4(
            nelements, result_name, slot, obj_vector_real)
        if auto_return:
            return nelements * ntotal, None, None

        obj = op2.obj
        # sort_method = op2.sort_method2()
        # sort_method = self.sort_method2
        sort_method = op2.sort_method
        if op2.use_vector and is_vectorized and op2.sort_method == 1:
            n = nelements * ntotal
            itotal = obj.ielement
            itotali = obj.itotal + nelements
            itotal2 = obj.itotal + nelements * nnodes_expected
            obj._times[obj.itime] = dt
            if obj.itime == 0:
                # (eid_device, cid, abcd, nnodes)
                ints = np.frombuffer(data, dtype=op2.idtype8).copy()
                try:
                    ints1 = ints.reshape(nelements, numwide_real)
                except ValueError:
                    msg = 'ints.shape=%s; size=%s ' % (str(ints.shape), ints.size)
                    msg += 'nelements=%s numwide_real=%s nelements*numwide=%s' % (
                        nelements, numwide_real, nelements * numwide_real)
                    raise ValueError(msg)
                eids = ints1[:, 0] // 10
                cids = ints1[:, 1]
                # nids = ints1[:, 4]
                assert eids.min() > 0, eids.min()
                obj.element_node[itotal:itotal2, 0] = np.repeat(eids, nnodes_expected)
                ints2 = ints1[:, 4:].reshape(nelements * nnodes_expected, 21)
                grid_device = ints2[:, 0]  # .reshape(nelements, nnodes_expected)

                # print('%s-grid_device=%s' % (op2.element_name, grid_device))
                unused_grid_device2 = np.repeat(grid_device, nnodes_expected)
                try:
                    obj.element_node[itotal:itotal2, 1] = grid_device
                except ValueError:
                    msg = '%s; nnodes=%s\n' % (op2.element_name, nnodes_expected)
                    msg += 'itotal=%s itotal2=%s\n' % (itotal, itotal2)
                    msg += 'grid_device.shape=%s; size=%s\n' % (str(grid_device.shape), grid_device.size)
                    # msg += 'nids=%s' % nids
                    raise ValueError(msg)
                # log.debug(f'cids = {np.unique(cids)}')
                obj.element_cid[itotal:itotali, 0] = eids
                obj.element_cid[itotal:itotali, 1] = cids

            floats = np.frombuffer(data, dtype=op2.fdtype8).reshape(nelements, numwide_real)[:, 4:]
            # 1     9    15   2    10   16  3   11  17   8
            # [oxx, oyy, ozz, txy, tyz, txz, o1, o2, o3, ovm]
            # isave = [1, 9, 15, 2, 10, 16, 3, 11, 17, 8]
            # (grid_device,
            # sxx, sxy, s1, a1, a2, a3, pressure, svm,
            # syy, syz, s2, b1, b2, b3,
            # szz, sxz, s3, c1, c2, c3)
            floats1 = floats.reshape(nelements * nnodes_expected, 21)  # [:, 1:] # drop grid_device

            # o1/o2/o3 is not max/mid/min.  They are not consistently ordered, so we force it.
            max_mid_min = np.vstack([
                floats1[:, 3],
                floats1[:, 11],
                floats1[:, 17],
            ]).T
            max_mid_min.sort(axis=1)
            assert max_mid_min.shape == (nelements * nnodes_expected, 3), max_mid_min.shape
            obj.data[obj.itime, itotal:itotal2, 6:9] = max_mid_min[:, [2, 1, 0]]

            # obj.data[obj.itime, itotal:itotal2, :] = floats1[:, isave]
            obj.data[obj.itime, itotal:itotal2, :6] = floats1[:, [1, 9, 15, 2, 10, 16]]
            obj.data[obj.itime, itotal:itotal2, 9] = floats1[:, 8]
            obj.itotal = itotal2
            obj.ielement = itotali
        else:
            if is_vectorized and op2.use_vector:  # pragma: no cover
                log.debug(f'vectorize CSolid real SORT{sort_method} from {op2.table_name}')
            n = oes_csolid_real(op2, data, obj,
                                nelements, dt,
                                element_name, nnodes_expected,
                                preline1, preline2)

    elif op2.format_code in [2, 3] and op2.num_wide == numwide_imag:  # complex
        op2.log.debug(f'numwide_imag={numwide_imag}')
        ntotal = numwide_imag * 4 * factor
        nelements = ndata // ntotal
        oes.ntotal += nelements * nnodes_expected
        auto_return, is_vectorized = op2._create_oes_object4(
            nelements, result_name, slot, obj_vector_complex)
        if auto_return:
            return nelements * ntotal, None, None

        obj = op2.obj

        if op2.use_vector and is_vectorized and op2.sort_method == 1:
            n = nelements * ntotal
            ielement = obj.ielement
            ielement2 = ielement + nelements
            itotal = obj.itotal
            itotal2 = itotal + nelements * nnodes_expected

            floats = np.frombuffer(data, dtype=op2.fdtype8).reshape(nelements, numwide_imag)
            floats1 = floats[:, 4:].reshape(nelements * nnodes_expected, 13).copy()
            obj._times[obj.itime] = dt
            if obj.itime == 0:
                ints = np.frombuffer(data, dtype=op2.idtype8).reshape(nelements, numwide_imag)
                ints1 = ints[:, 4:].reshape(nelements * nnodes_expected, 13)
                eids = ints[:, 0] // 10
                cids = ints[:, 1]
                nids = ints1[:, 0]
                # TODO: ctype, nodef not considered
                assert eids.min() > 0, eids.min()
                assert nids.min() >= 0, nids.min()
                eids2 = np.vstack([eids] * nnodes_expected).T.ravel()
                # nids2 = np.vstack([nids] * nnodes_expected).T.ravel()
                # print(nids2)
                obj.element_node[itotal:itotal2, 0] = eids2
                obj.element_node[itotal:itotal2, 1] = nids

                obj.element_cid[ielement:ielement2, 0] = eids
                obj.element_cid[ielement:ielement2, 1] = cids

            # 0 is nid
            isave1 = [1, 2, 3, 4, 5, 6]
            isave2 = [7, 8, 9, 10, 11, 12]
            real_imag = apply_mag_phase(floats1, is_magnitude_phase, isave1, isave2)
            obj.data[obj.itime, itotal:itotal2, :] = real_imag

            obj.itotal = itotal2
            obj.ielement = ielement2
        else:
            if is_vectorized and op2.use_vector:  # pragma: no cover
                log.debug(f'vectorize CSolid imag SORT{op2.sort_method}')
            n = oes_csolid_complex(op2, data, obj,
                                   nelements,  # nnodes,
                                   element_name, nnodes_expected,
                                   is_magnitude_phase)

    elif op2.format_code == 1 and op2.num_wide == numwide_random:  # random
        op2.log.debug(f'numwide_random={numwide_random}')
        if not op2.is_sort1:
            log.debug(f'support CSolid random SORT{op2.sort_method}')
            return ndata, None, None

        ntotal = numwide_random * 4
        nelements = ndata // ntotal
        assert ndata % ntotal == 0, ndata
        auto_return, is_vectorized = op2._create_oes_object4(
            nelements, result_name, slot, obj_vector_random)
        if auto_return:
            return nelements * ntotal, None, None

        obj = op2.obj
        if op2.use_vector and is_vectorized and 0:  # pragma: no cover
            n = nelements * 4 * op2.num_wide
            itotal = obj.ielement
            itotali = obj.itotal + nelements
            itotal2 = obj.itotal + nelements * nnodes_expected
            obj._times[obj.itime] = dt
            if obj.itime == 0:
                # (eid_device, cid, abcd, nnodes)
                ints = np.frombuffer(data, dtype=op2.idtype).copy()
                try:
                    ints1 = ints.reshape(nelements, numwide_real)
                except ValueError:
                    msg = 'ints.shape=%s; size=%s ' % (str(ints.shape), ints.size)
                    msg += 'nelements=%s numwide_real=%s nelements*numwide=%s' % (
                        nelements, numwide_real, nelements * numwide_real)
                    raise ValueError(msg)
                eids = ints1[:, 0] // 10
                cids = ints1[:, 1]
                # nids = ints1[:, 4]
                assert eids.min() > 0, eids.min()
                obj.element_node[itotal:itotal2, 0] = np.repeat(eids, nnodes_expected)
                ints2 = ints1[:, 4:].reshape(nelements * nnodes_expected, 21)
                grid_device = ints2[:, 0]  # .reshape(nelements, nnodes_expected)

                # print('%s-grid_device=%s' % (op2.element_name, grid_device))
                unused_grid_device2 = np.repeat(grid_device, nnodes_expected)
                try:
                    obj.element_node[itotal:itotal2, 1] = grid_device
                except ValueError:
                    msg = '%s; nnodes=%s\n' % (op2.element_name, nnodes_expected)
                    msg += 'itotal=%s itotal2=%s\n' % (itotal, itotal2)
                    msg += 'grid_device.shape=%s; size=%s\n' % (str(grid_device.shape), grid_device.size)
                    # msg += 'nids=%s' % nids
                    raise ValueError(msg)
                obj.element_cid[itotal:itotali, 0] = eids
                obj.element_cid[itotal:itotali, 1] = cids

            floats = np.frombuffer(data, dtype=op2.fdtype).reshape(nelements, numwide_real)[:, 4:]
            # 1     9    15   2    10   16  3   11  17   8
            # [oxx, oyy, ozz, txy, tyz, txz, o1, o2, o3, ovm]
            # isave = [1, 9, 15, 2, 10, 16, 3, 11, 17, 8]
            # (grid_device,
            # sxx, sxy, s1, a1, a2, a3, pressure, svm,
            # syy, syz, s2, b1, b2, b3,
            # szz, sxz, s3, c1, c2, c3)
            floats1 = floats.reshape(nelements * nnodes_expected, 21)  # [:, 1:] # drop grid_device

            # o1/o2/o3 is not max/mid/min.  They are not consistently ordered, so we force it.
            max_mid_min = np.vstack([
                floats1[:, 3],
                floats1[:, 11],
                floats1[:, 17],
            ]).T
            max_mid_min.sort(axis=1)
            assert max_mid_min.shape == (nelements * nnodes_expected, 3), max_mid_min.shape
            obj.data[obj.itime, itotal:itotal2, 6:9] = max_mid_min[:, [2, 1, 0]]

            # obj.data[obj.itime, itotal:itotal2, :] = floats1[:, isave]
            obj.data[obj.itime, itotal:itotal2, :6] = floats1[:, [1, 9, 15, 2, 10, 16]]
            obj.data[obj.itime, itotal:itotal2, 9] = floats1[:, 8]
            obj.itotal = itotal2
            obj.ielement = itotali
        else:
            if is_vectorized and op2.use_vector and obj.itime == 0:  # pragma: no cover
                log.debug(f'vectorize CSolid random SORT{op2.sort_method}')
            n = oes_csolid_random(op2, data, obj, nelements,
                                  element_name, nnodes_expected,
                                  preline1, preline2)

    elif op2.format_code in [2, 3] and op2.num_wide == numwide_random2_vm and \
            op2.table_name in {b'OESVM1', b'OSTRVM1', b'OESVM2', b'OSTRVM2'}:
        # op2.log.debug(f'numwide_random2={numwide_random2}')
        # raise RuntimeError(op2.code_information())
        ## a = 18
        ## b = 14
        ## a + b * nnodes = numwide_random3
        ## a + b * 4 = 74  # CTETRA
        ## a + b * 6 = 102 # CPENTA
        ## a + b * 8 = 130 # CHEXA-67
        # msg = 'OES-CHEXA-random-numwide=%s numwide_real=%s numwide_imag=%s numwide_random=%s' % (
        # op2.num_wide, numwide_real, numwide_imag, numwide_random)
        # return op2._not_implemented_or_skip(data, ndata, msg)

        # print('numwide real=%s imag=%s random=%s' % (numwide_real, numwide_imag, numwide_random))
        unused_num_wide_random = 4 + nnodes_expected * (17 - 4)

        # print('random2=%s' % num_wide_random)
        # print(op2.code_information())

        # print('numwide=%s numwide_random=%s attempt2=%s subcase=%s' % (
        # op2.num_wide, numwide_random, num_wide_random, op2.isubcase))
        ##msg = op2.code_information()
        # ntotal = 130
        ntotal = op2.num_wide * size
        nelements = ndata // ntotal

        auto_return, is_vectorized = op2._create_oes_object4(
            nelements, result_name, slot, obj_vector_complex_vm)
        if auto_return:
            return nelements * ntotal, None, None
        ## cid, coord_type, nactive_pnts,
        ##      nid, oxx, oyy, ozz, txy, tyz, txz
        ntotal1 = 4 * size  # 3*4
        ntotal2 = 14 * size
        ntotali = ntotal1 + nnodes_expected * ntotal2
        assert ntotali == ntotal, (ntotali, ntotal)

        # 2 CID I Stress Coordinate System
        # 3 CTYPE CHAR4 Coordinate System Type (BCD)
        # 4 NODEF I Number of Active Points
        #
        # 5 GRID I External grid identification number (0=center)
        # 6 SXR RS Normal in x
        # 7 SYR RS Normal in y
        # 8 SZR RS Normal in z
        # 9 TXYR RS Shear in xy
        # 10 TYZR RS Shear in yz
        # 11 TZXR RS Shear in zx
        # 12 SXI RS Normal in x
        # 13 SYI RS Normal in y
        # 14 SZI RS Normal in z
        # 15 TXYI RS Shear in xy
        # 16 TYZI RS Shear in yz
        # 17 TZXI RS Shear in zx
        # 18 VM RS von Mises
        # Words 5 through 18 repeat 005 times
        obj = op2.obj
        if op2.use_vector and is_vectorized and op2.sort_method == 1:
            n = nelements * ntotal
            itotal = obj.ielement
            itotali = obj.itotal + nelements
            itotal2 = obj.itotal + nelements * nnodes_expected
            obj._times[obj.itime] = dt
            if obj.itime == 0:
                # (eid_device, cid, abcd, nnodes)
                ints = np.frombuffer(data, dtype=op2.idtype8).copy()
                try:
                    ints1 = ints.reshape(nelements, numwide_random2_vm)
                except ValueError:
                    msg = 'ints.shape=%s; size=%s ' % (str(ints.shape), ints.size)
                    msg += 'nelements=%s numwide_real=%s nelements*numwide=%s' % (
                        nelements, numwide_real, nelements * numwide_real)
                    raise ValueError(msg)

                eids = ints1[:, 0] // 10
                cids = ints1[:, 1]
                # print(f'eids = {eids}')
                assert eids.min() > 0, eids.min()
                obj.element_node[itotal:itotal2, 0] = np.repeat(eids, nnodes_expected)
                ints2 = ints1[:, 4:].reshape(nelements * nnodes_expected, 14)
                grid_device = ints2[:, 0]  # .reshape(nelements, nnodes_expected)
                # print(ints1)
                # print(ints2)
                # print(nelements, nnodes_expected, nelements*nnodes_expected)
                # print(f'nids = {grid_device}')
                # print('ints2.shape =', ints2.shape)

                # print('%s-grid_device=%s' % (op2.element_name, grid_device))
                unused_grid_device2 = np.repeat(grid_device, nnodes_expected)
                try:
                    obj.element_node[itotal:itotal2, 1] = grid_device
                except ValueError:
                    msg = '%s; nnodes=%s\n' % (op2.element_name, nnodes_expected)
                    msg += 'itotal=%s itotal2=%s\n' % (itotal, itotal2)
                    msg += 'grid_device.shape=%s; size=%s\n' % (str(grid_device.shape), grid_device.size)
                    # msg += 'nids=%s' % nids
                    raise ValueError(msg)
                obj.element_cid[itotal:itotali, 0] = eids
                obj.element_cid[itotal:itotali, 1] = cids

            floats = np.frombuffer(data, dtype=op2.fdtype8).reshape(nelements, numwide_random2_vm)[:, 4:]
            # 1    2    3    4    5    6    7 - verify...
            # [oxx, oyy, ozz, txy, tyz, txz, ovm]
            # sxr, syr, szr, txyr, tyzr, txzr, \
            #     sxi, syi, szi, txyi, tyzi, txzi, ovm = out
            isave1 = [1, 2, 3, 4, 5, 6]
            isave2 = [7, 8, 9, 10, 11, 12]
            floats1 = floats.reshape(nelements * nnodes_expected, 14)
            real_imag = apply_mag_phase(floats1, is_magnitude_phase, isave1, isave2)
            # floats2 = floats1[:, 1:] # drop grid_device

            # o1/o2/o3 is not max/mid/min.  They are not consistently ordered, so we force it.
            obj.data[obj.itime, itotal:itotal2, -1] = floats1[:, -1]
            obj.data[obj.itime, itotal:itotal2, :-1] = real_imag
            obj.itotal = itotal2
            obj.ielement = itotali
        else:
            struct1 = Struct(op2._endian + b'2i 4s i')
            struct2 = Struct(op2._endian + b'i 13f')
            element_num = op2.element_type
            element_type = op2.element_name
            for i in range(nelements):
                edata = data[n:n + ntotal1]
                out = struct1.unpack(edata)
                (eid_device, cid, ctype, nodef) = out
                # print(eid_device, cid, ctype, nodef)
                eid, dt = get_eid_dt_from_eid_device(
                    eid_device, op2.nonlinear_factor, op2.sort_method)
                # obj.add_sort1(eid, dt, cid)
                obj.add_eid_sort1(element_num, element_type, dt, eid, cid, ctype, nodef)

                if op2.is_debug_file:
                    op2.binary_debug.write('%s - eid=%i; %s\n' % (preline1, eid, str(out)))
                n += ntotal1
                for inode in range(nnodes_expected):  # nodes pts, +1 for centroid (???)
                    grid, *out = struct2.unpack(data[n:n + ntotal2])
                    # print(' ', grid, out)
                    # if op2.is_debug_file:
                    # op2.binary_debug.write('%s - %s\n' % (preline2, str(out)))
                    # (grid_device, sxx, syy, sz, txy, tyz, txz) = out
                    #     ELEMENT-ID    GRID-ID    NORMAL-X       NORMAL-Y       NORMAL-Z       SHEAR-XY       SHEAR-YZ       SHEAR-ZX      VON MISES
                    # 0                     1048   2.823925E+03  -8.389856E+03  -1.373890E+03  -3.886024E+03  -9.440013E+02  -8.069132E+02
                    #                             -1.142355E+02   3.635427E+02   6.020642E+01   1.693545E+02   4.121148E+01   3.410083E+01   1.210361E+04

                    sxr, syr, szr, txyr, tyzr, txzr, \
                        sxi, syi, szi, txyi, tyzi, txzi, ovm = out
                    if is_magnitude_phase:
                        sx = polar_to_real_imag(sxr, sxi)
                        sy = polar_to_real_imag(syr, syi)
                        sz = polar_to_real_imag(szr, szi)
                        txy = polar_to_real_imag(txyr, txyi)
                        tyz = polar_to_real_imag(tyzr, tyzi)
                        txz = polar_to_real_imag(txzr, txzi)
                    else:
                        sx = complex(sxr, sxi)
                        sy = complex(syr, syi)
                        sz = complex(szr, szi)
                        txy = complex(txyr, txyi)
                        tyz = complex(tyzr, tyzi)
                        txz = complex(txzr, txzi)
                    # del grid, sx, sy, sz, txy, tyz, txz
                    # if eid == 1048:
                    # print(f'grid={grid}', sx, sy, sz, txy, tyz, txz)
                    # print(grid, out)
                    obj.add_node_sort1(dt, eid, grid, inode, sx, sy, sz, txy, tyz, txz, ovm)
                    n += ntotal2
        # etype_num = f'{op2.element_name}-{op2.element_type:d};'
        # msg = f'OES-{etype_num:<11s} complex2; table_name={op2.table_name}; numwide={op2.num_wide}'
        # return op2._not_implemented_or_skip(data, ndata, msg), None, None
    elif op2.format_code in [1, 2] and op2.num_wide == 67 and op2.element_name == 'CHEXA':  # CHEXA
        op2.log.debug(f'numwide complex HEXA?')
        # CTETRA:
        # CPYRAM:
        # CPENTA:
        # CHEXA:  67
        # 44 = 5 * 8 + 4  (TETRA)
        # 52 = 6 * 8 + 4  (PYRAM)
        # 60 = 7 * 8 + 4  (PENTA)
        # 76 = 9 * 8 + 4  (HEXA)
        msg = 'skipping random CHEXA; numwide=67'
        # print(op2.code_information())
        n = op2._not_implemented_or_skip(data, ndata, msg)
        nelements = None
        ntotal = None
    elif op2.format_code in [1, 2, 3] and op2.num_wide in [60] and op2.table_name in [b'OESXRMS1',
                                                                                      b'OESXNO1']:  # CPENTA
        op2.log.debug(f'numwide complex PENTA?')
        # bad
        # if op2.read_mode == 2:
        # ints    = (68011, 0, 805.28, 6,
        # 0,   0, 0, 0, 0, 0, 0, 0,
        # 120000, 0, 0, 0, 0, 0, 0, 0,
        # 120001, 0, 0, 0, 0, 0, 0, 0,
        # 120010, 0, 0, 0, 0, 0, 0, 0,
        # 120100, 0, 0, 0, 0, 0, 0, 0,
        # 120101, 0, 0, 0, 0, 0, 0, 0,
        # 120110, 0, 0, 0, 0, 0, 0, 0,

        # 68111, 0, 1145655879, 15,
        # 0,      1080284864, 1080296481, 1080284990, 1080279426, 1080285072, 1080285570, 1080285750,
        # 130000, 1080285656, 1080285656, 1080285656, 1080285162, 1080287537, 1080285308, 1080285794,
        # 130002, 1080285656, 1080285656, 1080285656, 1080284551, 1080287537, 1080285308, 1080285794,
        # 130020, 1080285656, 1080285656, 1080285656, 1080289401, 1080287537, 1080285308, 1080285794,
        # 130200, 1080285656, 1080285656, 1080285656, 1080285269, 1080287537, 1080285308, 1080285794,
        # 130202, 1080288409, 1080287759, 1080323139, 1080285308, 1080285512, 1080285308, 1080285874,
        # 130220, 1080285333, 1080285373, 1080285450, 1080287537, 1080287537, 1080285625, 1080285771)
        # floats  = (68011, 0.0, 805.28, 6,
        # 0.0,                   0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
        # 1.681558157189705e-40, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
        # 1.681572170174437e-40, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
        # 1.681698287036213e-40, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
        # 1.682959455654153e-40, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
        # 1.682973468638786e-40, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
        # 1.683099585500578e-40, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,

        # 9.544383970362762e-41, 0.0, 805.28, 2.1019476964872256e-44,
        # 0.0,                    3.559982299804687, 3.562752008438110, 3.5600123405456, 3.558685779571533, 3.560031890869140, 3.560150623321533, 3.56019353866577,
        # 1.8216880036222622e-40, 3.560171127319336, 3.560171127319336, 3.5601711273193, 3.560053348541259, 3.560619592666626, 3.560088157653808, 3.56020402908325,
        # 1.8217160295915487e-40, 3.560171127319336, 3.560171127319336, 3.5601711273193, 3.559907674789428, 3.560619592666626, 3.560088157653808, 3.56020402908325,
        # 1.8219682633151272e-40, 3.560171127319336, 3.560171127319336, 3.5601711273193, 3.561064004898071, 3.560619592666626, 3.560088157653808, 3.56020402908325,
        # 1.8244906005509118e-40, 3.560171127319336, 3.560171127319336, 3.5601711273193, 3.560078859329223, 3.560619592666626, 3.560088157653808, 3.56020402908325,
        # 1.8245186265201983e-40, 3.560827493667602, 3.560672521591186, 3.5691077709198, 3.560088157653808, 3.560136795043945, 3.560088157653808, 3.56022310256958,
        # 1.8247708602437768e-40, 3.560094118118286, 3.56010365486145,  3.5601220130920, 3.560619592666626, 3.560619592666626, 3.560163736343384, 3.56019854545593)
        # self.show_data(data, types='ifs')

        # C:\MSC.Software\simcenter_nastran_2019.2\tpl_post2\tr1081x.op2
        msg = 'skipping random CPENTA; numwide=60'
        n = op2._not_implemented_or_skip(data, ndata, msg)
        nelements = None
        ntotal = None
    elif op2.format_code in [1, 2, 3] and op2.num_wide == numwide_random3 and op2.table_name in [b'OESXRMS1',
                                                                                                 b'OESXNO1']:
        op2.log.debug(f'numwide_random3={numwide_random3}')
        # CTETRA: 44
        # CPYRAM: 52
        # CHEXA:  76
        # '      FREQUENCY =  0.000000E+00'
        # '                      S T R E S S E S   I N   H E X A H E D R O N   S O L I D   E L E M E N T S   ( H E X A )'
        # '                                    ( ROOT MEAN SQUARE; RMSSF SCALE FACTOR =  1.00E+00 )'
        # '0                   CORNER      --------------------------CENTER AND CORNER POINT STRESSES---------------------------'
        # '      ELEMENT-ID   GRID-ID      NORMAL-X    NORMAL-Y    NORMAL-Z      SHEAR-XY    SHEAR-YZ    SHEAR-ZX   VON MISES'
        # '0            6701        0GRID CS  8 GP'
        # '0                   CENTER     6.413E+01   2.752E-01   6.413E+01     6.496E+00   6.496E+00   3.842E+01   9.359E+01'
        # '0                    30000     5.367E+00   1.069E+01   5.367E+00     1.134E-15   1.134E-15   4.473E-15   5.321E+00'
        # C:\MSC.Software\simcenter_nastran_2019.2\tpl_post2\tr1081x.op2

        # 2 CID       I Stress Coordinate System
        # 3 CTYPE CHAR4 Coordinate System Type (BCD)
        # 4 NODEF     I Number of Active Points
        # 5 GRID      I External grid identification number (0=center)
        msg = f'skipping random {element_name}; numwide={numwide_random3}'

        ntotal = op2.num_wide * 4
        nelements = ndata // ntotal
        assert nelements > 0, nelements
        assert ndata % ntotal == 0, ndata

        auto_return, is_vectorized = op2._create_oes_object4(
            nelements, result_name, slot, obj_vector_random_vm)
        if auto_return:
            return nelements * ntotal, None, None
        if op2.read_mode == 1:
            return ndata, None, None

        obj = op2.obj
        if op2.use_vector and is_vectorized:  # pragma: no cover
            is_finished = True
            n = nelements * 4 * op2.num_wide
            itotal = obj.ielement
            itotali = obj.itotal + nelements
            itotal2 = obj.itotal + nelements * nnodes_expected
            obj._times[obj.itime] = dt
            if obj.itime == 0:
                # (eid_device, cid, abcd, nnodes)
                ints = np.frombuffer(data, dtype=op2.idtype).copy()
                try:
                    ints1 = ints.reshape(nelements, numwide_random3)
                except ValueError:
                    msg = 'ints.shape=%s; size=%s ' % (str(ints.shape), ints.size)
                    msg += 'nelements=%s numwide_random3=%s nelements*numwide=%s' % (
                        nelements, numwide_random3, nelements * numwide_random3)
                    raise ValueError(msg)
                eids = ints1[:, 0] // 10
                cids = ints1[:, 1]
                # nids = ints1[:, 4]
                assert eids.min() > 0, eids.min()
                obj.element_node[itotal:itotal2, 0] = np.repeat(eids, nnodes_expected)
                ints2 = ints1[:, 4:].reshape(nelements * nnodes_expected, 8)
                grid_device = ints2[:, 0]  # .reshape(nelements, nnodes_expected)

                # print('%s-grid_device=%s' % (op2.element_name, grid_device))
                unused_grid_device2 = np.repeat(grid_device, nnodes_expected)
                print(f'grid_device = {grid_device.tolist()}')
                print(f'eids = {eids.tolist()}')
                print(f'cids = {cids.tolist()}')
                try:
                    obj.element_node[itotal:itotal2, 1] = grid_device
                except ValueError:
                    msg = '%s; nnodes=%s\n' % (op2.element_name, nnodes_expected)
                    msg += 'itotal=%s itotal2=%s\n' % (itotal, itotal2)
                    msg += 'grid_device.shape=%s; size=%s\n' % (str(grid_device.shape), grid_device.size)
                    # msg += 'nids=%s' % nids
                    raise ValueError(msg)
                obj.element_cid[itotal:itotali, 0] = eids
                obj.element_cid[itotal:itotali, 1] = cids

            floats = np.frombuffer(data, dtype=op2.fdtype).reshape(nelements, numwide_random3)[:, 4:]

            # (grid_device, oxx, oyy, ozz, txy, tyz, txz, von_mises)
            floats1 = floats.reshape(nelements * nnodes_expected, 8)[:, 1:]  # drop grid_device
            assert floats1.shape == (nelements * nnodes_expected, 7), floats1.shape

            obj.data[obj.itime, itotal:itotal2, :] = floats1
            obj.itotal = itotal2
            obj.ielement = itotali
        else:
            if is_vectorized and op2.use_vector and obj.itime == 0:  # pragma: no cover
                log.debug(f'vectorize CSolid random SORT{op2.sort_method}')
            n = oes_csolid_random3(op2, data, obj, nelements,
                                   element_name, nnodes_expected,
                                   ntotal)
    # elif op2.format_code in [2, 3] and op2.num_wide == 76:  # imag
    # C:\MSC.Software\simcenter_nastran_2019.2\tpl_post2\tr1081x.op2
    # analysis_code = 5   Frequency
    # table_code    = 905 OESXNO1-OESXNO1C - Cumulative Root Mean Square output
    # format_code   = 3   Magnitude/Phase
    # result_type   = 3   Magnitude/Phase
    # sort_method   = 1
    # random_code   = 0
    # element_type  = 67  CHEXA
    # msg = op2.code_information()
    # return op2._not_implemented_or_skip(data, ndata, msg), None, None
    else:  # pragma: no cover
        numwide_random2 = 0
        raise RuntimeError(op2.code_information() +
                           f'\nnumwide real={numwide_real} imag={numwide_imag} '
                           f'random2={numwide_random2} random3={numwide_random3}')
    return n, nelements, ntotal


def oes_csolid_real(op2: OP2, data: bytes,
                    obj: RealSolidStressArray | RealSolidStrainArray,
                    nelements: int, dt: Any,
                    element_name: str, nnodes_expected: int,
                    preline1: str, preline2: str) -> int:
    """
    reads stress/strain for element type:
     - 39 : CTETRA
     - 67 : CHEXA
     - 68 : CPENTA

    """
    add_eid_sort_x = getattr(obj, 'add_eid_sort' + str(op2.sort_method))
    add_node_sort_x = getattr(obj, 'add_node_sort' + str(op2.sort_method))
    n = 0
    if op2.size == 4:
        fmt1 = op2._endian + op2._analysis_code_fmt + b'i4si'
        fmt2 = op2._endian + b'i20f'
    else:
        fmt1 = op2._endian + mapfmt(op2._analysis_code_fmt, op2.size) + b'q8sq'
        fmt2 = op2._endian + b'q20d'
    struct1 = Struct(fmt1)
    struct2 = Struct(fmt2)
    if op2.is_debug_file:
        msg = '%s-%s nelements=%s nnodes=%s; C=[sxx, sxy, s1, a1, a2, a3, pressure, svm,\n' % (
            op2.element_name, op2.element_type, nelements, nnodes_expected)
        msg += '                                 syy, syz, s2, b1, b2, b3,\n'
        msg += '                                 szz, sxz, s3, c1, c2, c3]\n'
        op2.binary_debug.write(msg)

    n16 = 16 * op2.factor
    n84 = 84 * op2.factor
    #if op2.is_sort1:
        #print('nelements =', nelements)
    #else:
        #print('ntimes =', nelements)
    for unused_i in range(nelements):
        edata = data[n:n+n16]
        out = struct1.unpack(edata)
        (eid_device, cid, unused_abcd, nnodes) = out
        eid, dt = get_eid_dt_from_eid_device(
            eid_device, op2.nonlinear_factor, op2.sort_method)
        #print(f'eid_device={eid_device} eid={eid} dt={dt}')
        if op2.is_debug_file:
            op2.binary_debug.write('%s - eid=%i; %s\n' % (preline1, eid, str(out)))

        assert nnodes < 21, 'print_block(data[n:n+16])'  #self.print_block(data[n:n+16])

        n += n16
        for inode in range(nnodes_expected):  # nodes pts, +1 for centroid (???)
            out = struct2.unpack(data[n:n + n84]) # 4*21 = 84
            if op2.is_debug_file:
                op2.binary_debug.write('%s - %s\n' % (preline2, str(out)))
            (grid_device,
             sxx, sxy, s1, a1, a2, a3, pressure, svm,
             syy, syz, s2, b1, b2, b3,
             szz, sxz, s3, c1, c2, c3) = out

            if op2.is_debug_file:
                op2.binary_debug.write('  eid=%s inode=%i; C=[%s]\n' % (
                    eid, grid_device, ', '.join(['%r' % di for di in out])))

            #if grid_device == 0:
                #grid = 'CENTER'
            #else:
                ##grid = (grid_device - device_code) // 10
                #grid = grid_device

            grid = grid_device
            a_cos = [a1, a2, a3]
            b_cos = [b1, b2, b3]
            c_cos = [c1, c2, c3]

            if inode == 0:
                #  this is correct, but fails
                #element_name = op2.element_name + str(nnodes)
                add_eid_sort_x(element_name, cid, dt, eid, grid,
                               sxx, syy, szz, sxy, syz, sxz, s1, s2, s3,
                               a_cos, b_cos, c_cos, pressure, svm)
            else:
                add_node_sort_x(dt, eid, inode, grid,
                                sxx, syy, szz, sxy, syz, sxz, s1, s2, s3,
                                a_cos, b_cos, c_cos, pressure, svm)
            n += n84
    return n


def oes_csolid_complex(op2: OP2, data: bytes,
                       obj: ComplexSolidStressArray | ComplexSolidStrainArray,
                       nelements: int, # nnodes: int,
                       element_name: str, nnodes_expected: int,
                       is_magnitude_phase: bool) -> int:
    n = 0
    if op2.size == 4:
        s1 = Struct(op2._endian + b'2i4si')
    else:
        s1 = Struct(op2._endian + b'2q8sq')
    s2 = Struct(op2._endian + mapfmt(b'i12f', op2.size))
    ntotal1 = 16 * op2.factor
    ntotal2 = 52 * op2.factor
    for unused_i in range(nelements):
        edata = data[n:n+ntotal1]
        n += ntotal1
        out = s1.unpack(edata)
        (eid_device, cid, ctype, nodef) = out
        eid, dt = get_eid_dt_from_eid_device(
            eid_device, op2.nonlinear_factor, op2.sort_method)
        if op2.is_debug_file:
            op2.binary_debug.write('  eid=%i C=[%s]\n' % (
                eid, ', '.join(['%r' % di for di in out])))

        #element_name = op2.element_name + str(nodef)  # this is correct, but has problems...
        obj.add_eid_sort1(op2.element_type, element_name, dt, eid, cid, ctype, nodef)
        for inode in range(nnodes_expected):
            edata = data[n:n+ntotal2]
            n += ntotal2
            out = s2.unpack(edata)
            (grid,
             exr, eyr, ezr, etxyr, etyzr, etzxr,
             exi, eyi, ezi, etxyi, etyzi, etzxi) = out
            #if grid == 0:
                #grid = 'CENTER'

            if is_magnitude_phase:
                ex = polar_to_real_imag(exr, exi)
                ey = polar_to_real_imag(eyr, eyi)
                ez = polar_to_real_imag(ezr, ezi)
                etxy = polar_to_real_imag(etxyr, etxyi)
                etyz = polar_to_real_imag(etyzr, etyzi)
                etzx = polar_to_real_imag(etzxr, etzxi)
            else:
                ex = complex(exr, exi)
                ey = complex(eyr, eyi)
                ez = complex(ezr, ezi)
                etxy = complex(etxyr, etxyi)
                etyz = complex(etyzr, etyzi)
                etzx = complex(etzxr, etzxi)

            if op2.is_debug_file:
                op2.binary_debug.write('       node%s=[%s]\n' % (
                    grid, ', '.join(['%r' % di for di in out])))
            obj.add_node_sort1(dt, eid, grid, inode,
                               ex, ey, ez, etxy, etyz, etzx)
    return n


def oes_csolid_random(op2: OP2, data: bytes,
                      obj: RandomSolidStressArray | RandomSolidStrainArray,
                      nelements: int,
                      element_name: str, nnodes_expected: int,
                      preline1: str, preline2: str) -> int:
    n = 0
    struct1 = Struct(op2._endian + op2._analysis_code_fmt + b'i4si')
    struct2 = Struct(op2._endian + b'i6f')
    if op2.is_debug_file and 0:
        msg = '%s-%s nelements=%s nnodes=%s; C=[sxx, sxy, s1, a1, a2, a3, pressure, svm,\n' % (
            op2.element_name, op2.element_type, nelements, nnodes_expected)
        msg += '                                 syy, syz, s2, b1, b2, b3,\n'
        msg += '                                 szz, sxz, s3, c1, c2, c3]\n'
        op2.binary_debug.write(msg)

    ntotal1 = 16 * op2.factor
    ntotal2 = 28 * op2.factor
    for unused_i in range(nelements):
        edata = data[n:n+ntotal1]
        out = struct1.unpack(edata)
        (eid_device, cid, unused_abcd, grid) = out

        eid, dt = get_eid_dt_from_eid_device(
            eid_device, op2.nonlinear_factor, op2.sort_method)
        assert eid > 0, eid

        if op2.is_debug_file and 0:
            op2.binary_debug.write('%s - eid=%i; %s\n' % (preline1, eid, str(out)))

        #assert nnodes < 21, 'print_block(data[n:n+16])'  #self.print_block(data[n:n+16])

        n += ntotal1
        for inode in range(nnodes_expected):  # nodes pts, +1 for centroid (???)
            #self.show_data(data[n:n+48])
            out = struct2.unpack(data[n:n + ntotal2]) # 4*7 = 28
            if op2.is_debug_file:
                op2.binary_debug.write('%s - %s\n' % (preline2, str(out)))
            (grid_device, sxx, syy, szz, txy, tyz, txz) = out

            if op2.is_debug_file:
                op2.binary_debug.write('  eid=%s inode=%i; C=[%s]\n' % (
                    eid, grid_device, ', '.join(['%r' % di for di in out])))

            #if grid_device == 0:
                #grid = 'CENTER'
            #else:
                ##grid = (grid_device - device_code) // 10
                #grid = grid_device

            grid = grid_device
            if inode == 0:
                #  this is correct, but fails
                #element_name = op2.element_name + str(nnodes)
                obj.add_eid_sort1(element_name, cid, dt, eid, grid,
                                  sxx, syy, szz, txy, tyz, txz)
            else:
                obj.add_node_sort1(dt, eid, inode, grid,
                                   sxx, syy, szz, txy, tyz, txz)
            n += ntotal2
    return n



def oes_csolid_random3(op2: OP2, data: bytes, obj: RandomSolidVMStressArray, nelements: int,
                       element_name: str, nnodes_expected: int,
                       ntotal: int) -> int:
    # print(msg)
    n = 0
    assert op2.size == 4, op2.size
    fmt = b'ii4si' + b'i7f' * nnodes_expected
    structi = Struct(fmt)

    for ielement in range(nelements):
        datai = data[n:n + ntotal]
        # op2.show_data(datai, types='if')
        out = structi.unpack(datai)
        eid_device, cid, grid, nnodes, *other = out
        eid, dt = get_eid_dt_from_eid_device(
            eid_device, op2.nonlinear_factor, op2.sort_method)
        #print(eid, cid, grid, nnodes)
        # [oxx, oyy, ozz, txy, tyz, txz, von_mises]
        # 6 SX   RS Normal in x
        # 7 SY   RS Normal in y
        # 8 SZ   RS Normal in z
        # 9 TXY  RS Shear in xy
        # 10 TYZ RS Shear in yz
        # 11 TZX RS Shear in zx
        # 12 RMSVM RS RMS von Mises
        for inode in range(nnodes_expected):
            #print('   ', other[8 * inode], other[8 * inode + 1:8 * (inode + 1)])
            node_id = other[8 * inode]
            oxx, oyy, ozz, txy, tyz, txz, vm = other[8 * inode + 1:8 * (inode + 1)]
            obj.add_node_sort1(dt, eid, inode, node_id, oxx, oyy, ozz, txy, tyz, txz, vm)
        n += ntotal
    return n
