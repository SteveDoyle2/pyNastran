"""
Defines the Real/Complex Forces created by:
    GPFORCE = ALL

"""
from __future__ import annotations
from struct import Struct
from typing import TYPE_CHECKING
import numpy as np

from pyNastran.op2.op2_interface.op2_reader import mapfmt
from pyNastran.op2.op2_helper import polar_to_real_imag
from pyNastran.op2.tables.utils import get_is_slot_saved, get_eid_dt_from_eid_device
from pyNastran.op2.tables.ogf_gridPointForces.ogf_objects import (
    RealGridPointForcesArray, ComplexGridPointForcesArray)

if TYPE_CHECKING:  # pragma: no cover
    from pyNastran.op2.op2 import OP2

class OGPF:
    def __init__(self, op2: OP2):
        self.op2 = op2

    def _read_ogpf1_3(self, data: bytes, ndata: int):
        self.op2._op2_readers.reader_opg._read_opg1_3(data, ndata)  # TODO: this is wrong...

    def _read_ogpf2_3(self, data: bytes, ndata: int):
        self.op2._op2_readers.reader_opg._read_opg2_3(data, ndata)  # TODO: this is wrong...

    def _read_ogpf1_4(self, data: bytes, ndata: int) -> int:
        op2 = self.op2
        prefix = ''
        if op2.table_code == 19:  # grid point force balance
            if op2.table_name == b'OGPFB1':
                pass
            elif op2.table_name == b'RAGEATC':
                prefix = 'RAGEATC.'
            elif op2.table_name == b'RAGCONS':
                prefix = 'RAGCONS.'
            else:
                msg = f'table_name={op2.table_name} table_code={op2.table_code}'
                raise RuntimeError(msg)
            n = self.read_grid_point_forces1(data, ndata, prefix=prefix)
        else:
            raise NotImplementedError(op2.table_code)
        return n

    def _read_ogpf2_4(self, data: bytes, ndata: int) -> int:
        op2 = self.op2
        prefix = ''
        if op2.table_code == 19:  # grid point force balance
            assert op2.table_name == b'OGPFB2', op2.code_information()
            n = self.read_grid_point_forces2(data, ndata, prefix=prefix)
        else:
            raise NotImplementedError(op2.table_code)
        return n

    def read_grid_point_forces1(self, data: bytes, ndata: int, prefix: str='') -> int:
        """table_code = 19"""
        op2 = self.op2
        op2._setup_op2_subcase('GPFORCE')
        dt = op2.nonlinear_factor
        n = 0
        is_magnitude_phase = op2.is_magnitude_phase()

        factor = op2.factor
        if op2.thermal == 0:
            result_name = prefix + 'grid_point_forces'
            is_saved, slot = get_is_slot_saved(op2, result_name)
            if not is_saved:
                return ndata

            if op2.num_wide == 10:
                ntotal = 40 * factor # 4*10
                nnodes = ndata // ntotal
                obj_vector_real = RealGridPointForcesArray
                auto_return, is_vectorized = op2._op2_readers.reader_oes._create_ntotal_object(
                    nnodes, result_name, slot, obj_vector_real)
                if auto_return:
                    return nnodes * ntotal

                obj = op2.obj
                if op2.is_debug_file:
                    op2.binary_debug.write('  GPFORCE\n')
                    op2.binary_debug.write('  [cap, gpforce1, gpforce2, ..., cap]\n')
                    op2.binary_debug.write('  cap = %i  # assume 1 cap when there could have been multiple\n' % len(data))
                    op2.binary_debug.write('  gpforce1 = [nid_device, eid, elem_name, f1, f2, f3, m1, m2, m3]\n')
                    op2.binary_debug.write('  nnodes=%i\n' % nnodes)

                if op2.use_vector and is_vectorized:
                    # self.itime = 0
                    # self.ielement = 0
                    # self.itotal = 0
                    #self.ntimes = 0
                    #self.nelements = 0
                    n = nnodes * ntotal

                    istart = obj.itotal
                    iend = istart + nnodes
                    obj._times[obj.itime] = dt

                    itime = obj.itime
                    if itime == 0 or obj.is_unique:
                        ints = np.frombuffer(data, dtype=op2.idtype8).reshape(nnodes, 10).copy()

                        nids = ints[:, 0] // 10
                        eids = ints[:, 1]
                        if op2.size == 4:
                            strings = np.frombuffer(data, dtype=op2._uendian + 'S8').reshape(nnodes, 5).copy()
                            strings_save = strings[:, 1]
                        else:
                            strings = np.frombuffer(data, dtype=op2._uendian + 'S4').reshape(nnodes, 20).copy()
                            strings_save = np.core.defchararray.add(strings[:, 4], strings[:, 6])

                        if obj.is_unique:
                            obj.node_element[itime, istart:iend, 0] = nids
                            obj.node_element[itime, istart:iend, 1] = eids
                            obj.element_names[itime, istart:iend] = strings_save
                        else:
                            obj.node_element[istart:iend, 0] = nids
                            obj.node_element[istart:iend, 1] = eids
                            obj.element_names[istart:iend] = strings_save


                    floats = np.frombuffer(data, dtype=op2.fdtype8).reshape(nnodes, 10)
                    #[f1, f2, f3, m1, m2, m3]
                    obj.data[itime, istart:iend, :] = floats[:, 4:].copy()
                    #obj._times[obj.itime] = dt
                    #obj.itotal = itotal2
                    if op2.is_debug_file:
                        if itime != 0:
                            ints = np.frombuffer(data, dtype=op2.idtype).reshape(nnodes, 10)
                            strings = np.frombuffer(data, dtype=op2._uendian + 'S8').reshape(nnodes, 5)
                        for i in range(iend - istart):
                            op2.binary_debug.write('  nid=%s - (%s, %s, %s, %s, %s, %s, %s, %s, %s)\n' % (
                                ints[i, 0] // 10,
                                ints[i, 0], ints[i, 1], strings[i, 1],
                                floats[i, 4], floats[i, 5], floats[i, 6],
                                floats[i, 7], floats[i, 8], floats[i, 9], ))
                else:
                    if op2.size == 4:
                        fmt = op2._endian + b'ii8s6f'
                    else:
                        fmt = op2._endian + b'qq16s6d'
                    s = Struct(fmt)
                    for i in range(nnodes):
                        edata = data[n:n+ntotal]
                        out = s.unpack(edata)
                        (nid_device, eid, elem_name, f1, f2, f3, m1, m2, m3) = out
                        nid = nid_device // 10
                        elem_name = elem_name.strip()
                        if op2.is_debug_file:
                            op2.binary_debug.write('  nid=%s - %s\n' % (nid, str(out)))
                        op2.obj.add_sort1(dt, nid, eid, elem_name, f1, f2, f3, m1, m2, m3)
                        n += ntotal
            elif op2.num_wide == 16:
                # complex
                ntotal = 64
                nnodes = ndata // ntotal
                assert op2.size == 4, op2.size
                obj_vector_real = ComplexGridPointForcesArray
                auto_return, is_vectorized = op2._op2_readers.reader_oes._create_ntotal_object(
                    nnodes, result_name, slot, obj_vector_real)
                if auto_return:
                    return nnodes * op2.num_wide * 4

                obj = op2.obj
                is_vectorized = False
                if op2.use_vector and is_vectorized:
                    # self.itime = 0
                    # self.ielement = 0
                    # self.itotal = 0
                    #self.ntimes = 0
                    #self.nelements = 0
                    n = nnodes * op2.num_wide * 4

                    istart = obj.itotal
                    iend = istart + nnodes
                    obj._times[obj.itime] = dt

                    if obj.itime == 0:
                        ints = np.frombuffer(data, dtype=op2.idtype).reshape(nnodes, 16)
                        nids = ints[:, 0] // 10
                        eids = ints[:, 1]
                        obj.node_element[istart:iend, 0] = nids
                        obj.node_element[istart:iend, 1] = eids
                        strings = np.frombuffer(data, dtype=op2._uendian + 'S8').reshape(nnodes, 8)
                        obj.element_names[istart:iend] = strings[:, 1].copy()

                    floats = np.frombuffer(data, dtype=op2.fdtype).reshape(nnodes, 16)
                    #[f1, f2, f3, m1, m2, m3]
                    obj.data[obj.itime, istart:iend, :] = floats[:, 4:].copy()
                else:
                    s = Struct(op2._endian + b'ii8s12f')

                    #if op2.is_debug_file:
                        #op2.binary_debug.write('  GPFORCE\n')
                        #op2.binary_debug.write('  [cap, gpforce1, gpforce2, ..., cap]\n')
                        #op2.binary_debug.write('  cap = %i  # assume 1 cap when there could have been multiple\n' % len(data))
                        #op2.binary_debug.write('  gpforce1 = [nid_device, eid, elem_name, f1, f2, f3, m1, m2, m3]\n')
                        #op2.binary_debug.write('  nnodes=%i\n' % nnodes)

                    for i in range(nnodes):
                        edata = data[n:n+ntotal]
                        out = s.unpack(edata)
                        (nid_device, eid, elem_name,
                         f1r, f2r, f3r, m1r, m2r, m3r,
                         f1i, f2i, f3i, m1i, m2i, m3i) = out
                        nid = nid_device // 10
                        elem_name = elem_name.strip()
                        if op2.is_debug_file:
                            op2.binary_debug.write('  nid=%s - %s\n' % (nid, str(out)))

                        if is_magnitude_phase:
                            f1 = polar_to_real_imag(f1r, f1i)
                            f2 = polar_to_real_imag(f2r, f2i)
                            f3 = polar_to_real_imag(f3r, f3i)
                            m1 = polar_to_real_imag(m1r, m1i)
                            m2 = polar_to_real_imag(m2r, m2i)
                            m3 = polar_to_real_imag(m3r, m3i)
                        else:
                            f1 = complex(f1r, f1i)
                            f2 = complex(f2r, f2i)
                            f3 = complex(f3r, f3i)
                            m1 = complex(m1r, m1i)
                            m2 = complex(m2r, m2i)
                            m3 = complex(m3r, m3i)

                        op2.obj.add_sort1(dt, nid, eid, elem_name, f1, f2, f3, m1, m2, m3)
                        n += ntotal
            else:
                raise NotImplementedError(op2.code_information())
        else:
            raise NotImplementedError(op2.code_information())
            #msg = op2.code_information()
            #return self._not_implemented_or_skip(data, ndata, msg)
        return n

    def read_grid_point_forces2(self, data: bytes, ndata: int, prefix: str='') -> int:
        """table_code = 19"""
        op2 = self.op2
        #op2._setup_op2_subcase('GPFORCE')
        dt = op2.nonlinear_factor
        n = 0
        is_magnitude_phase = op2.is_magnitude_phase()

        assert op2.thermal == 0, op2.code_information()
        #if op2.thermal == 0:

        result_name = prefix + 'grid_point_forces'
        #print(op2.code_information())
        is_saved, slot = get_is_slot_saved(result_name)
        if not is_saved:
            return ndata

        if op2.num_wide == 8:
            ntotal = 32 * op2.factor # 4*8
            nnodes = ndata // ntotal
            #obj_vector_real = RealGridPointForcesArray
            #auto_return, is_vectorized = op2._op2_readers.reader_oes._create_ntotal_object(
                #nnodes, result_name, slot, obj_vector_real)
            auto_return = data is None
            if auto_return:
                return nnodes * ntotal
            if op2.size == 4:
                fmt = op2._endian + op2._analysis_code_fmt + b'4s 6f'
            else:
                fmt = op2._endian + mapfmt(op2._analysis_code_fmt, 8) + b'8s 6d'

            eid = op2.nonlinear_factor
            s = Struct(fmt)
            for i in range(nnodes):
                edata = data[n:n+ntotal]
                out = s.unpack(edata)
                #print(out)
                (dt, word, f1, f2, f3, m1, m2, m3) = out
                #print(dt, eid, word)
                assert word in {b'G   '}, out
                #eid, dt = get_eid_dt_from_eid_device(
                    #eid_device, op2.nonlinear_factor, op2.sort_method)
                #nid = nid_device // 10
                #elem_name = elem_name.strip()
                #if op2.is_debug_file:
                    #op2.binary_debug.write('  nid=%s - %s\n' % (nid, str(out)))
                #op2.obj.add_sort1(dt, nid, eid, elem_name, f1, f2, f3, m1, m2, m3)
                n += ntotal
        else:
            raise RuntimeError(op2.code_information())

        return n
