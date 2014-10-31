from six.moves import range
from struct import Struct

from pyNastran.op2.tables.ogf_gridPointForces.ogf_Objects import (
    RealGridPointForces, ComplexGridPointForces)


class OGPF(object):
    def __init__(self):
        pass

    def _read_ogpf1_3(self, data):
        self._read_opg1_3(data)  # TODO: this is wrong...

    def _read_ogpf1_4(self, data):
        if self.read_mode == 1:
            return len(data)
        if self.table_code == 19:  # grid point force balance
            assert self.table_name in ['OGPFB1'], 'table_name=%s table_code=%s' % (self.table_name, self.table_code)
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
            result_name = 'gridPointForces'
            if result_name not in self._saved_results:
                return len(data)
            self._found_results.add(result_name)
            if self.num_wide == 10:
                self.create_transient_object(self.gridPointForces, RealGridPointForces)
                s = Struct(b'ii8s6f')
                ntotal = 40
                nnodes = len(data) // ntotal
                for i in range(nnodes):
                    eData = data[n:n+ntotal]
                    out = s.unpack(eData)
                    (ekey, eid, elemName, f1, f2, f3, m1, m2, m3) = out
                    ekey = (ekey - self.device_code) // 10
                    elemName = elemName.strip()
                    #data = (eid, elemName, f1, f2, f3, m1, m2, m3)
                    self.obj.add(dt, ekey, eid, elemName, f1, f2, f3, m1, m2, m3)
                    #print "eid/dt/freq=%s eid=%-6s eName=%-8s f1=%g f2=%g f3=%g m1=%g m2=%g m3=%g" %(ekey,eid,elemName,f1,f2,f3,m1,m2,m3)
                    n += ntotal
            else:
                raise NotImplementedError('num_wide = %s' % (self.num_wide))

            #complex_obj = complexGridPointForcesObject

            #self._read_table(data, storage_obj, real_obj, complex_obj, 'node')
        #elif self.thermal == 1:
            #result_name = 'thermalLoadVectors'
            #storage_obj = self.thermalLoadVectors
            #real_obj = ThermalLoadVectorObject
            #complex_obj = None
            #self._read_table(data, storage_obj, real_obj, complex_obj, 'node')
        else:
            raise NotImplementedError(self.thermal)
        return n
