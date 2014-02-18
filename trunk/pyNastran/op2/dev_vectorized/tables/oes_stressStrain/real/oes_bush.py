from __future__ import (nested_scopes, generators, division, absolute_import,
                        print_function, unicode_literals)

from numpy import zeros
import pandas as pd
from .oes_objects import StressObject, StrainObject
from pyNastran.f06.f06_formatting import writeFloats13E

class RealBushResults(object):

    def __init__(self):
        self.shape = {}
        self.eType = {}
        self.data = None
        self._ncount = 0
        self._ielement_start = None
        self._ielement_end = None

    def is_real(self):
        return True

    def is_imaginary(self):
        return False

    def get_transients(self):
        k = self.translations.keys()
        k.sort()
        return k

    def _increase_size(self, dt, nelements):
        if dt in self.shape:  # default dictionary
            self.shape[dt] += nelements
        else:
            self.shape[dt] = nelements
        #print("shape =", self.shape)

    def _get_shape(self):
        ndt = len(self.shape)
        dts = self.shape.keys()
        shape0 = dts[0]
        nelements = self.shape[shape0]
        #print("ndt=%s nnodes=%s dts=%s" % (str(ndt), str(nnodes), str(dts)))
        return ndt, nelements, dts

    def _increment(self, nelements):
        self._ielement_start += nelements
        self._ielement_end += nelements
        return self._ielement_start, self._ielement_end

    def _preallocate(self, dt, nelements):
        ndt, nelements_size, dts = self._get_shape()
        #print("ndt=%s nelements_size=%s dts=%s" % (ndt, nelements_size, str(dts)))

        if self._ielement_start is not None:
            return (self._ielement_start, self._ielement_start + nelements)
        print('----definition----')
        n = ndt * nelements_size
        if self._ncount != 0:
            asfd
        self._ncount += 1
        self._ielement_start = 0
        self._ielement_end = nelements

        data = {}
        columns = []
        if dts[0] is not None:
            name = self.data_code['name']
            if isinstance(dt, int):
                data[name] = pd.Series(zeros((n), dtype='int32'))
            else:
                data[name] = pd.Series(zeros((n), dtype='float32'))
            columns.append(name)

        #element_data['element_id'] = pd.Series(zeros((nelements_size), dtype='int32'))
        #element_data['element_type'] = pd.Series(zeros(nelements_size, dtype='str'))

        data['element_id'] = pd.Series(zeros((n), dtype='int32'))
        #columns.append('element_type')

        #data['grid_type'] = pd.Series(zeros(ndt), dtype='int32'))
        #data['grid_type_str'] = pd.Series(zeros(nnodes), dtype='str'))
        #print('n =', n)

        data['element_type'] = pd.Series(zeros(nelements_size, dtype='str'))
        data['T1'] = pd.Series(zeros((n), dtype='float32'))
        data['T2'] = pd.Series(zeros((n), dtype='float32'))
        data['T3'] = pd.Series(zeros((n), dtype='float32'))

        data['R1'] = pd.Series(zeros((n), dtype='float32'))
        data['R2'] = pd.Series(zeros((n), dtype='float32'))
        data['R3'] = pd.Series(zeros((n), dtype='float32'))

        columns += ['element_id', 'element_type', 'T1', 'T2', 'T3', 'R1', 'R2', 'R3',]

        self.data = pd.DataFrame(data, columns=columns)
        return (self._ielement_start, self._ielement_end)

    def _finalize(self, dt):
        ndt, nelements, dts = self._get_shape()

        if dt is not None and dt != max(dts):
            return
        #print("----finalize----")

        #grid_type_str = []
        #for grid_type in self.grid_type:
            #grid_type_str.append('C' if grid_type==0 else grid_type)
        #self.grid_type_str = pd.Series(grid_type_str, dtype='str')

        if dts[0] is not None:
            name = self.data_code['name'] # dt
            self.data = self.data.set_index([name, 'element_id'])
        else:
            self.data = self.data.set_index(['element_id'])
        #print("final\n", self.data.to_string())
        #print('---BushStressObject---')
        #print(self.data.to_string())

    def get_stats(self):
        msg = self._get_data_code()
        ndt, nelements, dts = self._get_shape()

        if self.nonlinear_factor is not None:  # transient
            name = self.data_code['name']
            dt_string = '%s, ' % name
            msg.append('  real type=%s n%ss=%s nelements=%s\n'
                       % (self.__class__.__name__, name, ndt, nelements))
        else:
            dt_string = ''
            msg.append('  real type=%s nelements=%s\n' % (self.__class__.__name__,
                                                          nelements))
        msg.append('  data        : index :  %selement_id\n' % dt_string)
        msg.append('              : result:  element_type, T1, T2, T3, R1, R2, R3\n')
        return msg

    def _write_f06_helper(self, words, pageStamp, f, pageNum, is_mag_phase):
        ndt, nelements, dts = self._get_shape()

        eids = self.data['element_id']
        T1 = self.data['T1']
        T2 = self.data['T2']
        T3 = self.data['T3']

        R1 = self.data['R1']
        R2 = self.data['R2']
        R3 = self.data['R3']

        msg = words
        for i in xrange(nelements):
            #eType = self.eType[eid]
            eid = eids[i]
            vals = [T1[i], T2[i], T3[i],
                    R1[i], R2[i], R3[i], ]
            (vals2, isAllZeros) = writeFloats13E(vals)
            [t1, t2, t3, r1, r2, r3] = vals2
            msg.append('0%27i     %13s %13s %13s %13s %13s %-s\n' % (eid, t1, t2, t3, r1, r2, r3.rstrip()))
        msg.append(pageStamp + str(pageNum) + '\n')
        f.write(''.join(msg))
        return pageNum

    def _write_f06_transient_helper(self, words, header, pageStamp, f, pageNum, is_mag_phase):
        ndt, nelements, dts = self._get_shape()
        ntotal = ndt * nelements

        eids = self.data['element_id']
        T1 = self.data['T1']
        T2 = self.data['T2']
        T3 = self.data['T3']

        R1 = self.data['R1']
        R2 = self.data['R2']
        R3 = self.data['R3']

        msg = []
        for n in xrange(ntotal):
            header[1] = ' %s = %10.4E\n' % (self.data_code['name'], dt)
            msg += header + words
            for ii in xrange(nelements):
                i = ii + n * nelements
                #eType = self.eType[eid]
                eid = eids[i]
                vals = [T1[i], T2[i], T3[i],
                        R1[i], R2[i], R3[i], ]
                (vals2, isAllZeros) = writeFloats13E(vals)
                [t1, t2, t3, r1, r2, r3] = vals2
                msg.append('0%8i   %13s  %13s  %13s  %13s  %13s  %-s\n' % (eid, t1, t2, t3, r1, r2, r3.rstrip()))
            msg.append(pageStamp + str(pageNum) + '\n')
            f.write(''.join(msg))
            msg = ['']
            pageNum += 1
        return pageNum - 1


class BushStressObject(RealBushResults, StressObject):

    def __init__(self, data_code, is_sort1, isubcase, dt, read_mode):
        RealBushResults.__init__(self)
        StressObject.__init__(self, data_code, isubcase, read_mode)

    def add_f06_data(self, data, transient):
        if transient is None:
            for line in data:
                (eType, eid, tx, ty, tz, rx, ry, rz) = line
                self.eType[eid] = 'CBUSH'
                self.translations[eid] = [tx, ty, tz]
                self.rotations[eid] = [rx, ry, rz]
            return

        (dtName, dt) = transient
        self.data_code['name'] = dtName

        if dt not in self.translations:
            self.update_dt(self.data_code, dt)

        for line in data:
            (eType, eid, tx, ty, tz, rx, ry, rz) = line
            self.eType[eid] = 'CBUSH'
            self.translations[dt][eid] = [tx, ty, tz]
            self.rotations[dt][eid] = [rx, ry, rz]

    def write_f06(self, header, pageStamp, f, pageNum=1, is_mag_phase=False):
        if self.nonlinear_factor is not None:
            return self._write_f06_transient(header, pageStamp, f, pageNum, is_mag_phase)

        words = header + [
            '                                  S T R E S S E S   I N   B U S H   E L E M E N T S        ( C B U S H )\n\n',
            '                  ELEMENT-ID        STRESS-TX     STRESS-TY     STRESS-TZ    STRESS-RX     STRESS-RY     STRESS-RZ \n',
        ]
        return self._write_f06_helper(words, pageStamp, f, pageNum, is_mag_phase)

    def _write_f06_transient(self, header, pageStamp, f, pageNum=1, is_mag_phase=False):
        words = [
            '                                  S T R E S S E S   I N   B U S H   E L E M E N T S        ( C B U S H )\n\n',
            '                  ELEMENT-ID        STRESS-TX     STRESS-TY     STRESS-TZ    STRESS-RX     STRESS-RY     STRESS-RZ \n',
        ]
        return self._write_f06_transient_helper(words, header, pageStamp, f, pageNum, is_mag_phase)


class BushStrainObject(RealBushResults, StrainObject):

    def __init__(self, data_code, is_sort1, isubcase, dt, read_mode):
        RealBushResults.__init__(self)
        StrainObject.__init__(self, data_code, isubcase, read_mode)

    def add_f06_data(self, data, transient):
        if transient is None:
            for line in data:
                (eType, eid, tx, ty, tz, rx, ry, rz) = line
                self.eType[eid] = 'CBUSH'
                self.translations[eid] = [tx, ty, tz]
                self.rotations[eid] = [rx, ry, rz]
            return

        (dtName, dt) = transient
        self.data_code['name'] = dtName

        if dt not in self.translations:
            self.update_dt(self.data_code, dt)

        for line in data:
            (eType, eid, tx, ty, tz, rx, ry, rz) = line
            self.eType[eid] = 'CBUSH'
            self.translations[dt][eid] = [tx, ty, tz]
            self.rotations[dt][eid] = [rx, ry, rz]

    def write_f06(self, header, pageStamp, f, pageNum=1, is_mag_phase=False):
        if self.nonlinear_factor is not None:
            return self._write_f06_transient(header, pageStamp, pageNum, f, is_mag_phase)

        words = header + [
            '                                  S T R A I N S   I N   B U S H   E L E M E N T S        ( C B U S H )\n\n',
            '                  ELEMENT-ID        STRESS-TX     STRESS-TY     STRESS-TZ    STRESS-RX     STRESS-RY     STRESS-RZ \n',
        ]
        return self._write_f06_helper(words, pageStamp, f, pageNum, is_mag_phase)

    def _write_f06_transient(self, header, pageStamp, f, pageNum=1, is_mag_phase=False):
        words = [
            '                                  S T R A I N S   I N   B U S H   E L E M E N T S        ( C B U S H )\n\n',
            '                  ELEMENT-ID        STRESS-TX     STRESS-TY     STRESS-TZ    STRESS-RX     STRESS-RY     STRESS-RZ \n',
        ]
        return self._write_f06_helper(words, header, pageStamp, f, pageNum, is_mag_phase)