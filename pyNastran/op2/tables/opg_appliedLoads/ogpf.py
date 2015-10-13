from six.moves import range
from struct import Struct

from pyNastran.op2.op2_common import OP2Common
from pyNastran.op2.tables.ogf_gridPointForces.ogf_Objects import (
    RealGridPointForces, ComplexGridPointForces)


class OGPF(OP2Common):
    def __init__(self):
        OP2Common.__init__(self)

    def _read_ogpf1_3(self, data):
        self._read_opg1_3(data)  # TODO: this is wrong...

    def _read_ogpf1_4(self, data):
        if self.read_mode == 1:
            return len(data)
        if self.table_code == 19:  # grid point force balance
            assert self.table_name in [b'OGPFB1'], 'table_name=%s table_code=%s' % (self.table_name, self.table_code)
            n = self._read_grid_point_forces(data)
        else:
            raise NotImplementedError(self.table_code)
        return n

    def _read_grid_point_forces(self, data):
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
            if self.num_wide == 10:
                self.create_transient_object(self.grid_point_forces, RealGridPointForces)
                s = Struct(self._endian + b'ii8s6f')
                ntotal = 40
                nnodes = len(data) // ntotal

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
                    ekey = (ekey - self.device_code) // 10
                    elemName = elemName.strip()
                    #data = (eid, elemName, f1, f2, f3, m1, m2, m3)
                    if self.is_debug_file:
                        self.binary_debug.write('  nid=%s - %s\n' % (ekey, str(out)))

                    self.obj.add(dt, ekey, eid, elemName, f1, f2, f3, m1, m2, m3)
                    #print "eid/dt/freq=%s eid=%-6s eName=%-8s f1=%g f2=%g f3=%g m1=%g m2=%g m3=%g" %(ekey,eid,elemName,f1,f2,f3,m1,m2,m3)
                    n += ntotal
            else:
                msg = self.code_information()
                return self._not_implemented_or_skip(data, msg)

            #complex_obj = complexGridPointForcesObject

            #self._read_table(data, storage_obj, real_obj, complex_obj, 'node')
        #elif self.thermal == 1:
            #result_name = 'thermal_load_vectors'
            #storage_obj = self.thermal_load_vectors
            #real_obj = ThermalLoadVectorObject
            #complex_obj = None
            #self._read_table(data, storage_obj, real_obj, complex_obj, 'node')
        else:
            msg = self.code_information()
            return self._not_implemented_or_skip(data, msg)
        return n
