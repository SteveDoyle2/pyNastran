"""
Defines the Real/Complex Forces created by:
    GPFORCE = ALL

"""
from struct import Struct
import numpy as np

from pyNastran.op2.op2_helper import polar_to_real_imag
from pyNastran.op2.op2_interface.op2_common import OP2Common
from pyNastran.op2.tables.ogf_gridPointForces.ogf_objects import (
    RealGridPointForcesArray, ComplexGridPointForcesArray)


class OGPF(OP2Common):
    def __init__(self):
        OP2Common.__init__(self)

    def _read_ogpf1_3(self, data, ndata):
        self._read_opg1_3(data, ndata)  # TODO: this is wrong...

    def _read_ogpf1_4(self, data, ndata):
        prefix = ''
        if self.table_code == 19:  # grid point force balance
            if self.table_name == b'OGPFB1':
                pass
            elif self.table_name == b'RAGEATC':
                prefix = 'RAGEATC.'
            elif self.table_name == b'RAGCONS':
                prefix = 'RAGCONS.'
            else:
                msg = 'table_name=%s table_code=%s' % (self.table_name, self.table_code)
                raise RuntimeError(msg)
            n = self._read_grid_point_forces(data, ndata, prefix=prefix)
        else:
            raise NotImplementedError(self.table_code)
        return n

    def _read_grid_point_forces(self, data, ndata, prefix=''):
        """
        table_code = 19
        """
        self._setup_op2_subcase('GPFORCE')
        dt = self.nonlinear_factor
        n = 0
        is_magnitude_phase = self.is_magnitude_phase()

        if self.thermal == 0:
            result_name = prefix + 'grid_point_forces'
            if self._results.is_not_saved(result_name):
                return ndata
            self._results._found_result(result_name)
            slot = self.get_result(result_name)

            if self.num_wide == 10:
                ntotal = 40 * self.factor # 4*10
                nnodes = ndata // ntotal
                obj_vector_real = RealGridPointForcesArray
                auto_return, is_vectorized = self._create_ntotal_object(
                    nnodes, result_name, slot, obj_vector_real)
                if auto_return:
                    return nnodes * ntotal

                obj = self.obj
                if self.is_debug_file:
                    self.binary_debug.write('  GPFORCE\n')
                    self.binary_debug.write('  [cap, gpforce1, gpforce2, ..., cap]\n')
                    self.binary_debug.write('  cap = %i  # assume 1 cap when there could have been multiple\n' % len(data))
                    self.binary_debug.write('  gpforce1 = [nid_device, eid, elem_name, f1, f2, f3, m1, m2, m3]\n')
                    self.binary_debug.write('  nnodes=%i\n' % nnodes)

                if self.use_vector and is_vectorized:
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
                        ints = np.frombuffer(data, dtype=self.idtype8).reshape(nnodes, 10).copy()

                        nids = ints[:, 0] // 10
                        eids = ints[:, 1]
                        if self.size == 4:
                            strings = np.frombuffer(data, dtype=self._uendian + 'S8').reshape(nnodes, 5).copy()
                            strings_save = strings[:, 1]
                        else:
                            strings = np.frombuffer(data, dtype=self._uendian + 'S4').reshape(nnodes, 20).copy()
                            strings_save = np.core.defchararray.add(strings[:, 4], strings[:, 6])

                        if obj.is_unique:
                            obj.node_element[itime, istart:iend, 0] = nids
                            obj.node_element[itime, istart:iend, 1] = eids
                            obj.element_names[itime, istart:iend] = strings_save
                        else:
                            obj.node_element[istart:iend, 0] = nids
                            obj.node_element[istart:iend, 1] = eids
                            obj.element_names[istart:iend] = strings_save


                    floats = np.frombuffer(data, dtype=self.fdtype8).reshape(nnodes, 10)
                    #[f1, f2, f3, m1, m2, m3]
                    obj.data[itime, istart:iend, :] = floats[:, 4:].copy()
                    #obj._times[obj.itime] = dt
                    #obj.itotal = itotal2
                    if self.is_debug_file:
                        if itime != 0:
                            ints = np.frombuffer(data, dtype=self.idtype).reshape(nnodes, 10)
                            strings = np.frombuffer(data, dtype=self._uendian + 'S8').reshape(nnodes, 5)
                        for i in range(iend - istart):
                            self.binary_debug.write('  nid=%s - (%s, %s, %s, %s, %s, %s, %s, %s, %s)\n' % (
                                ints[i, 0] // 10,
                                ints[i, 0], ints[i, 1], strings[i, 1],
                                floats[i, 4], floats[i, 5], floats[i, 6],
                                floats[i, 7], floats[i, 8], floats[i, 9], ))
                else:
                    if self.size == 4:
                        fmt = self._endian + b'ii8s6f'
                    else:
                        fmt = self._endian + b'qq16s6d'
                    s = Struct(fmt)
                    for i in range(nnodes):
                        edata = data[n:n+ntotal]
                        out = s.unpack(edata)
                        (nid_device, eid, elem_name, f1, f2, f3, m1, m2, m3) = out
                        nid = nid_device // 10
                        elem_name = elem_name.strip()
                        if self.is_debug_file:
                            self.binary_debug.write('  nid=%s - %s\n' % (nid, str(out)))
                        self.obj.add_sort1(dt, nid, eid, elem_name, f1, f2, f3, m1, m2, m3)
                        n += ntotal
            elif self.num_wide == 16:
                # complex
                ntotal = 64
                nnodes = ndata // ntotal
                assert self.size == 4, self.size
                obj_vector_real = ComplexGridPointForcesArray
                auto_return, is_vectorized = self._create_ntotal_object(
                    nnodes, result_name, slot, obj_vector_real)
                if auto_return:
                    return nnodes * self.num_wide * 4

                obj = self.obj
                is_vectorized = False
                if self.use_vector and is_vectorized:
                    # self.itime = 0
                    # self.ielement = 0
                    # self.itotal = 0
                    #self.ntimes = 0
                    #self.nelements = 0
                    n = nnodes * self.num_wide * 4

                    istart = obj.itotal
                    iend = istart + nnodes
                    obj._times[obj.itime] = dt

                    if obj.itime == 0:
                        ints = np.frombuffer(data, dtype=self.idtype).reshape(nnodes, 16)
                        nids = ints[:, 0] // 10
                        eids = ints[:, 1]
                        obj.node_element[istart:iend, 0] = nids
                        obj.node_element[istart:iend, 1] = eids
                        strings = np.frombuffer(data, dtype=self._uendian + 'S8').reshape(nnodes, 8)
                        obj.element_names[istart:iend] = strings[:, 1].copy()

                    floats = np.frombuffer(data, dtype=self.fdtype).reshape(nnodes, 16)
                    #[f1, f2, f3, m1, m2, m3]
                    obj.data[obj.itime, istart:iend, :] = floats[:, 4:].copy()
                else:
                    s = Struct(self._endian + b'ii8s12f')

                    #if self.is_debug_file:
                        #self.binary_debug.write('  GPFORCE\n')
                        #self.binary_debug.write('  [cap, gpforce1, gpforce2, ..., cap]\n')
                        #self.binary_debug.write('  cap = %i  # assume 1 cap when there could have been multiple\n' % len(data))
                        #self.binary_debug.write('  gpforce1 = [nid_device, eid, elem_name, f1, f2, f3, m1, m2, m3]\n')
                        #self.binary_debug.write('  nnodes=%i\n' % nnodes)

                    for i in range(nnodes):
                        edata = data[n:n+ntotal]
                        out = s.unpack(edata)
                        (nid_device, eid, elem_name,
                         f1r, f2r, f3r, m1r, m2r, m3r,
                         f1i, f2i, f3i, m1i, m2i, m3i) = out
                        nid = nid_device // 10
                        elem_name = elem_name.strip()
                        if self.is_debug_file:
                            self.binary_debug.write('  nid=%s - %s\n' % (nid, str(out)))

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

                        self.obj.add_sort1(dt, nid, eid, elem_name, f1, f2, f3, m1, m2, m3)
                        n += ntotal
            else:
                raise NotImplementedError(self.code_information())
        else:
            raise NotImplementedError(self.code_information())
            #msg = self.code_information()
            #return self._not_implemented_or_skip(data, ndata, msg)
        return n
