from six import b
from six.moves import range
from struct import Struct
from numpy import fromstring

from pyNastran.op2.op2_common import OP2Common
from pyNastran.op2.tables.ogf_gridPointForces.ogf_Objects import (
    RealGridPointForcesArray, ComplexGridPointForces)


class OGPF(OP2Common):
    def __init__(self):
        OP2Common.__init__(self)

    def _read_ogpf1_3(self, data, ndata):
        self._read_opg1_3(data, ndata)  # TODO: this is wrong...

    def _read_ogpf1_4(self, data, ndata):
        if self.table_code == 19:  # grid point force balance
            assert self.table_name in [b'OGPFB1'], 'table_name=%s table_code=%s' % (self.table_name, self.table_code)
            n = self._read_grid_point_forces(data, ndata)
        else:
            raise NotImplementedError(self.table_code)
        return n

    def _read_grid_point_forces(self, data, ndata):
        """
        table_code = 19
        """
        dt = self.nonlinear_factor
        n = 0
        if self.thermal == 0:
            result_name = 'grid_point_forces'
            if self._results.is_not_saved(result_name):
                return len(data)
            self._results._found_result(result_name)
            slot = getattr(self, result_name)

            if self.num_wide == 10:
                ntotal = 40
                nnodes = len(data) // ntotal
                obj_vector_real = RealGridPointForcesArray
                auto_return, is_vectorized = self._create_ntotal_object(
                    nnodes, result_name, slot, obj_vector_real)
                if auto_return:
                    return nnodes * self.num_wide * 4

                obj = self.obj
                if self.use_vector and is_vectorized: #  and self.element_type in [144]
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
                        ints = fromstring(data, dtype=self.idtype).reshape(nnodes, 10)
                        nids = ints[:, 0] // 10
                        eids = ints[:, 1]
                        obj.node_element[istart:iend, 0] = nids
                        obj.node_element[istart:iend, 1] = eids
                        strings = fromstring(data, dtype=self._endian + 'S8').reshape(nnodes, 5)#[:, 2:3]
                        #a = strings[:, 0]
                        #aa = strings[:, 1]
                        obj.element_names[istart:iend] = strings[:, 1]

                    floats = fromstring(data, dtype=self.fdtype).reshape(nnodes, 10)
                    #[o1, o2, t12, t1z, t2z, angle, major, minor, ovm]
                    obj.data[obj.itime, istart:iend, :] = floats[:, 4:]
                else:
                    s = Struct(b(self._endian + 'ii8s6f'))

                    if self.is_debug_file:
                        self.binary_debug.write('  GPFORCE\n')
                        self.binary_debug.write('  [cap, gpforce1, gpforce2, ..., cap]\n')
                        self.binary_debug.write('  cap = %i  # assume 1 cap when there could have been multiple\n' % len(data))
                        self.binary_debug.write('  gpforce1 = [ekey, eid, elemName, f1, f2, f3, m1, m2, m3]\n')
                        self.binary_debug.write('  nnodes=%i\n' % nnodes)

                    for i in range(nnodes):
                        eData = data[n:n+ntotal]
                        out = s.unpack(eData)
                        (ekey, eid, elemName, f1, f2, f3, m1, m2, m3) = out
                        ekey = ekey // 10
                        elemName = elemName.strip()
                        #data = (eid, elemName, f1, f2, f3, m1, m2, m3)
                        if self.is_debug_file:
                            self.binary_debug.write('  nid=%s - %s\n' % (ekey, str(out)))

                        self.obj.add(dt, ekey, eid, elemName, f1, f2, f3, m1, m2, m3)
                        #print "eid/dt/freq=%s eid=%-6s eName=%-8s f1=%g f2=%g f3=%g m1=%g m2=%g m3=%g" %(ekey,eid,elemName,f1,f2,f3,m1,m2,m3)
                        n += ntotal
            else:
                raise NotImplementedError(self.code_information())
                #msg = self.code_information()
                #return self._not_implemented_or_skip(data, msg)

            #complex_obj = complexGridPointForcesObject

            #self._read_table(data, storage_obj, real_obj, complex_obj, 'node')
        #elif self.thermal == 1:
            #result_name = 'thermal_load_vectors'
            #storage_obj = self.thermal_load_vectors
            #real_obj = ThermalLoadVectorObject
            #complex_obj = None
            #self._read_table(data, storage_obj, real_obj, complex_obj, 'node')
        else:
            raise NotImplementedError(self.code_information())
            #msg = self.code_information()
            #return self._not_implemented_or_skip(data, msg)
        return n
