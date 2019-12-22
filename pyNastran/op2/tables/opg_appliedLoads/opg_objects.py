from typing import List
import numpy as np
from pyNastran.utils.numpy_utils import integer_types
from pyNastran.op2.result_objects.op2_objects import ScalarObject
from pyNastran.f06.f06_formatting import write_floats_13e, write_imag_floats_13e


class AppliedLoadsVectorArray(ScalarObject):
    def __init__(self, data_code, is_sort1, isubcase, dt):
        ScalarObject.__init__(self, data_code, isubcase)
        #self.dt = dt
        self.ntimes = 0
        self.itotal = 0
        self.eids = None
        self.sources = None
        self.data = None # forces/moments

    def data_type(self):
        raise NotImplementedError()

    def _reset_indices(self):
        self.itotal = 0

    def build(self):
        """sizes the vectorized attributes of the AppliedLoadsVectorArray"""
        self.eids = np.zeros(self.itotal, dtype='int32')
        self.sources = np.zeros(self.itotal, dtype='|S8')
        #[f1, f2, f3, m1, m2, m3]
        self.data = np.zeros((self.ntimes, self.itotal, 6), dtype=self.data_type())

    def get_stats(self, short=False) -> List[str]:
        if not self.is_built:
            return [
                '<%s>\n' % self.__class__.__name__,
                '  ntimes: %i\n' % self.ntimes,
                '  ntotal: %i\n' % self.ntotal,
            ]
        #ngrids = len(self.gridTypes)
        msg = []

        ntimes = len(self._times)
        #len(self.node_gridtype)
        #nnodes, two = self.node_gridtype.shape
        nelements = 0
        ntimes = self.data.shape[0]
        if self.nonlinear_factor not in (None, np.nan):  # transient
            msg.append('  type=%s ntimes=%s nelements=%s\n'
                       % (self.__class__.__name__, ntimes, nelements))
        else:
            msg.append('  type=%s nelements=%s\n'
                       % (self.__class__.__name__, nelements))
        msg.append('  data: [f1, f2, f3, m1, m2, m3] shape=%s dtype=%s\n'
                   % ([int(i) for i in self.data.shape], self.data.dtype))
        msg.append('  sources, eids\n')
        msg += self.get_data_code()
        return msg

    def add_sort1(self, node_id, eid, source, v1, v2, v3, v4, v5, v6):
        """unvectorized method for adding SORT1 transient data"""
        #raise NotImplementedError('AppliedLoadsVector')
        assert isinstance(eid, integer_types) and eid > 0, 'dt=%s eid=%s' % (dt, eid)
        msg = "node_id=%s v1=%s v2=%s v3=%s" % (node_id, v1, v2, v3)
        assert 0 < node_id < 1000000000, msg
        #assert nodeID not in self.forces

        #[f1, f2, f3, m1, m2, m3]
        self.data[self.itime, self.itotal, :] = [v1, v2, v3, v4, v5, v6]


class RealAppliedLoadsVectorArray(AppliedLoadsVectorArray):
    def __init__(self, data_code, is_sort1, isubcase, dt):
        AppliedLoadsVectorArray.__init__(self, data_code, isubcase, dt)

    def data_type(self):
        raise 'float32'

    def write_f06(self, f06_file, header=None, page_stamp='PAGE %s',
                  page_num=1, is_mag_phase=False, is_sort1=True):
        if header is None:
            header = []
        words = ['                      APPLIED LOADS VECTOR\n',
                 '\n',
                 '      EID SOURCE FX FY FZ MX MY MZ\n']
        #ntimes, ntotal = self.data.shape[:2]

        eids = self.eids
        for itime, dt in enumerate(self._times):
            if self.nonlinear_factor not in (None, np.nan):
                if isinstance(dt, float):
                    header[1] = ' %s = %10.4E\n' % (self.data_code['name'], dt)
                else:
                    header[1] = ' %s = %10i\n' % (self.data_code['name'], dt)
            f06_file.write(''.join(header + words))

            f1 = self.data[itime, :, 0]
            f2 = self.data[itime, :, 1]
            f3 = self.data[itime, :, 2]
            m1 = self.data[itime, :, 3]
            m2 = self.data[itime, :, 4]
            m3 = self.data[itime, :, 5]
            source = ''
            for eid, f1i, f2i, f3i, m1i, m2i, m3i in zip(eids, f1, f2, f3, m1, m2, m3):
                vals = [f1i, f2i, f3i, m1i, m2i, m3i]
                vals2 = write_floats_13e(vals)
                (dx, dy, dz, rx, ry, rz) = vals2
                f06_file.write('%14i %6s     %-13s  %-13s  %-13s  %-13s  %-13s  %s\n' % (
                    # TODO: fix this...
                    eid, source, dx, dy, dz, rx, ry, rz))
            f06_file.write(page_stamp % page_num)
            page_num += 1
        return page_num-1


class ComplexAppliedLoadsVectorArray(AppliedLoadsVectorArray):
    def __init__(self, data_code, is_sort1, isubcase, dt):
        AppliedLoadsVectorArray.__init__(self, data_code, isubcase, dt)

    def data_type(self):
        raise 'float32'

    def write_f06(self, f06_file, header=None, page_stamp='PAGE %s',
                  page_num=1, is_mag_phase=False, is_sort1=True):
        if header is None:
            header = []
        words = ['                      APPLIED LOADS VECTOR\n',
                 '\n',
                 '      EID SOURCE FX FY FZ MX MY MZ\n']
        #ntimes, ntotal, size = self.data.shape

        eids = self.eids
        for itime, dt in enumerate(self._times):
            if self.nonlinear_factor not in (None, np.nan):
                if isinstance(dt, float):
                    header[1] = ' %s = %10.4E\n' % (self.data_code['name'], dt)
                else:
                    header[1] = ' %s = %10i\n' % (self.data_code['name'], dt)
            f06_file.write(''.join(header + words))

            f1 = self.data[itime, :, 0]
            f2 = self.data[itime, :, 1]
            f3 = self.data[itime, :, 2]
            m1 = self.data[itime, :, 3]
            m2 = self.data[itime, :, 4]
            m3 = self.data[itime, :, 5]
            source = ''
            #node_id = ''
            for eid, f1i, f2i, f3i, m1i, m2i, m3i in zip(eids, f1, f2, f3, m1, m2, m3):
                vals = [f1i, f2i, f3i, m1i, m2i, m3i]
                vals2 = write_imag_floats_13e(vals, is_mag_phase)
                (dxr, dxi, dyr, dyi, dzr, dzi,
                 rxr, rxi, ryr, ryi, rzr, rzi) = vals2  # TODO :verify
                f06_file.write('%14i %6s     %-13s  %-13s  %-13s  %-13s  %-13s  %s\n'
                               '%14s %6s     %-13s  %-13s  %-13s  %-13s  %-13s  %s\n' % (
                                   eid, source, dxr, dyr, dzr, rxr, ryr, rzr,
                                   '', '', dxi, dyi, dzi, rxi, ryi, rzi))
            f06_file.write(page_stamp % page_num)
            page_num += 1
        return page_num-1
