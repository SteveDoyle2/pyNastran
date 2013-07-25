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

        if dt != max(dts):
            return
        print("----finalize----")

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
        print('---BushStressObject---')
        #print(self.data.to_string())

    def get_stats(self):
        msg = self.get_data_code()
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

    def write_f06(self, header, pageStamp, pageNum=1, f=None, is_mag_phase=False):
        if self.nonlinear_factor is not None:
            return self._write_f06_transient(header, pageStamp, pageNum, f, is_mag_phase)

        msg = header + [
            '                                  S T R E S S E S   I N   B U S H   E L E M E N T S        ( C B U S H )\n\n',
            '                  ELEMENT-ID        STRESS-TX     STRESS-TY     STRESS-TZ    STRESS-RX     STRESS-RY     STRESS-RZ \n',
        ]

        for eid, (tx, ty, tz) in sorted(self.translations.iteritems()):
            eType = self.eType[eid]
            (rx, ry, rz) = self.rotations[eid]

            vals = [tx, ty, tz, rx, ry, rz]
            (vals2, isAllZeros) = writeFloats13E(vals)
            [tx, ty, tz, rx, ry, rz] = vals2
            msg.append('0%8i   %13s  %13s  %13s  %13s  %13s  %-s\n' %
                       (eid, tx, ty, tz, rx, ry, rz.rstrip()))

        msg.append(pageStamp + str(pageNum) + '\n')
        return (''.join(msg), pageNum)

    def _write_f06_transient(self, header, pageStamp, pageNum=1, f=None, is_mag_phase=False):
        words = [
            '                                  S T R E S S E S   I N   B U S H   E L E M E N T S        ( C B U S H )\n\n',
            '                  ELEMENT-ID        STRESS-TX     STRESS-TY     STRESS-TZ    STRESS-RX     STRESS-RY     STRESS-RZ \n',
        ]
        msg = []
        for dt, translations in sorted(self.translations.iteritems()):
            header[1] = ' %s = %10.4E\n' % (self.data_code['name'], dt)
            msg += header + words
            for eid, (tx, ty, tz) in sorted(translations.iteritems()):
                eType = self.eType[eid]
                (rx, ry, rz) = self.rotations[dt][eid]

                vals = [tx, ty, tz, rx, ry, rz]
                (vals2, isAllZeros) = writeFloats13E(vals)
                [tx, ty, tz, rx, ry, rz] = vals2
                msg.append('0%8i   %13s  %13s  %13s  %13s  %13s  %-s\n' %
                           (eid, tx, ty, tz, rx, ry, rz.rstrip()))

            msg.append(pageStamp + str(pageNum) + '\n')
            pageNum += 1
        return (''.join(msg), pageNum - 1)

    def __repr__(self):
        raise NotImplementedError('CBUSH')
        if self.nonlinear_factor is not None:
            return self.__reprTransient__()

        msg = '---BUSH STRESS---\n'
        msg += '%-6s %6s ' % ('EID', 'eType')
        headers = ['s1', 's2', 's3', 's4', 'Axial']
        for header in headers:
            msg += '%8s ' % header
        msg += '\n'

        for eid, S1s in sorted(self.s1.iteritems()):
            eType = self.eType[eid]
            axial = self.axial[eid]

            s1 = self.s1[eid]
            s2 = self.s2[eid]
            s3 = self.s3[eid]
            s4 = self.s4[eid]
            msg += '%-6i %6s ' % (eid, eType)
            vals = [s1[0], s2[0], s3[0], s4[0], axial]
            for val in vals:
                if abs(val) < 1e-6:
                    msg += '%8s ' % '0'
                else:
                    msg += '%8i ' % val
            msg += '\n'

            msg += '%s ' % (' ' * 13)
            vals = [s1[1], s2[1], s3[1], s4[1]]
            for val in vals:
                if isinstance(val, str):
                    msg += '%8s ' % val
                elif abs(val) < 1e-6:
                    msg += '%8s ' % '0'
                else:
                    msg += '%8i ' % val
            msg += '\n'

            #msg += "eid=%-4s eType=%s s1=%-4i s2=%-4i s3=%-4i s4=%-4i axial=-%5i\n" %(eid,eType,s1[0],s2[0],s3[0],s4[0],axial)
            #msg += "%s                s1=%-4i s2=%-4i s3=%-4i s4=%-4i %s\n"         %(' '*4,    s1[1],s2[1],s3[1],s4[1])
        return msg

    def __reprTransient__(self):
        raise NotImplementedError('CBUSH')
        msg = '---BUSH STRESS---\n'
        msg += '%-6s %6s ' % ('EID', 'eType')
        headers = ['s1', 's2', 's3', 's4', 'Axial', 'sMax', 'sMin']
        for header in headers:
            msg += '%8s ' % header
        msg += '\n'

        for dt, S1ss in sorted(self.s1.iteritems()):
            msg += '%s = %g\n' % (self.data_code['name'], dt)
            for eid, S1s in sorted(S1ss.iteritems()):
                eType = self.eType[eid]
                axial = self.axial[dt][eid]

                s1 = self.s1[dt][eid]
                s2 = self.s2[dt][eid]
                s3 = self.s3[dt][eid]
                s4 = self.s4[dt][eid]
                msg += '%-6i %6s ' % (eid, eType)
                vals = [s1[0], s2[0], s3[0], s4[0], axial]
                for val in vals:
                    msg += '%8s %8s' % (val.real, val.imag)
                msg += '\n'

                msg += '%s ' % (' ' * 13)
                vals = [s1[1], s2[1], s3[1], s4[1]]
                for val in vals:
                    if isinstance(val, str):
                        msg += '%8s ' % val
                    elif abs(val) < 1e-6:
                        msg += '%8s %8s' % (val.real, val.imag)
                msg += '\n'

                #msg += "eid=%-4s eType=%s s1=%-4i s2=%-4i s3=%-4i s4=%-4i axial=-%5i\n" %(eid,eType,s1[0],s2[0],s3[0],s4[0],axial)
                #msg += "%s                s1=%-4i s2=%-4i s3=%-4i s4=%-4i %s\n"         %(' '*4,    s1[1],s2[1],s3[1],s4[1])
        return msg


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

    def write_f06(self, header, pageStamp, pageNum=1, f=None, is_mag_phase=False):
        raise NotImplementedError('CBUSH')
        return 'BushStress write_f06 not implemented...', pageNum
        if self.nonlinear_factor is not None:
            return self._write_f06_transient(header, pageStamp, pageNum, f, is_mag_phase)

        msg = header + [
            '                                  S T R A I N S    I N   B A R   E L E M E N T S          ( C B A R )\n',
            '  ELEMENT        SA1            SA2            SA3            SA4           AXIAL          SA-MAX         SA-MIN     M.S.-T\n',
            '    ID.          SB1            SB2            SB3            SB4           STRAIN         SB-MAX         SB-MIN     M.S.-C\n',
        ]
        for eid, E1s in sorted(self.e1.iteritems()):
            eType = self.eType[eid]
            axial = self.axial[eid]

            e1 = self.e1[eid]
            e2 = self.e2[eid]
            e3 = self.e3[eid]
            e4 = self.e4[eid]
            vals = [e1[0], e2[0], e3[0], e4[0], axial,
                    e1[1], e2[1], e3[1], e4[1]]
            (vals2, isAllZeros) = writeFloats13E(vals)
            [e10, e20, e30, e40, axial,
             e11, e21, e31, e41] = vals2

            msg.append('0%8i   %13s  %13s  %13s  %13s  %13s  %13s  %13s %-s\n' % (eid, e10, e20, e30, e40, axial.rstrip()))
            msg.append(' %8s   %13s  %13s  %13s  %13s  %13s  %13s  %13s %-s\n' % ('', e11, e21, e31, e41.rstrip()))
        msg.append(pageStamp + str(pageNum) + '\n')
        return (''.join(msg), pageNum)

    def _write_f06_transient(self, header, pageStamp, pageNum=1, f=None, is_mag_phase=False):
        raise NotImplementedError('CBUSH')
        words = [
            '                                  S T R A I N S    I N   B A R   E L E M E N T S           ( C B A R )\n',
            '  ELEMENT        SA1            SA2            SA3            SA4           AXIAL          SA-MAX         SA-MIN     M.S.-T\n',
            '    ID.          SB1            SB2            SB3            SB4           STRAIN         SB-MAX         SB-MIN     M.S.-C\n',
        ]
        msg = []
        for dt, E1s in sorted(self.e1.iteritems()):
            header[1] = ' %s = %10.4E\n' % (self.data_code['name'], dt)
            msg += header + words
            for eid, e1s in sorted(E1s.iteritems()):
                eType = self.eType[eid]
                axial = self.axial[eid]

                e1 = self.e1[eid]
                e2 = self.e2[eid]
                e3 = self.e3[eid]
                e4 = self.e4[eid]
                vals = [e1[0], e2[0], e3[0], e4[0], axial,
                        e1[1], e2[1], e3[1], e4[1]]
                (vals2, isAllZeros) = writeFloats13E(vals)
                [e10, e20, e30, e40,
                 e11, e21, e31, e41] = vals2

                msg.append('0%8i   %13s  %13s  %13s  %13s  %13s  %13s  %13s %-s\n' % (eid, e10, e20, e30, e40, axial.rstrip()))
                msg.append(' %8s   %13s  %13s  %13s  %13s  %13s  %13s  %13s %-s\n' % ('', e11, e21, e31, e41.rstrip()))
            msg.append(pageStamp + str(pageNum) + '\n')
            pageNum += 1
        return (''.join(msg), pageNum - 1)

    def __repr__(self):
        raise NotImplementedError('CBUSH')
        if self.nonlinear_factor is not None:
            return self.__reprTransient__()

        msg = '---BUSH STRAIN---\n'
        msg += '%-8s %6s ' % ('EID', 'eType')
        headers = ['e1', 'e2', 'e3', 'e4', 'Axial', 'eMax', 'eMin']
        for header in headers:
            msg += '%10s ' % header
        msg += '\n'

        for eid, E1s in sorted(self.e1.iteritems()):
            eType = self.eType[eid]
            axial = self.axial[eid]
            e1 = self.e1[eid]
            e2 = self.e2[eid]
            e3 = self.e3[eid]
            e4 = self.e4[eid]
            msg += '%-8i %6s ' % (eid, eType)
            vals = [e1[0], e2[0], e3[0], e4[0], axial]
            for val in vals:
                if abs(val) < 1e-6:
                    msg += '%10s ' % '0'
                else:
                    msg += '%10.3g ' % val
            msg += '\n'

            msg += '%s ' % (' ' * 17)
            vals = [e1[1], e2[1], e3[1], e4[1]]
            for val in vals:
                if isinstance(val, str):
                    msg += '%10s ' % val
                elif abs(val) < 1e-6:
                    msg += '%10s ' % '0'
                else:
                    msg += '%10.3g ' % val
            msg += '\n'

            #msg += "eid=%-4s eType=%s s1=%-4i s2=%-4i s3=%-4i s4=%-4i axial=-%5i\n" %(eid,eType,s1[0],s2[0],s3[0],s4[0],axial)
            #msg += "%s                s1=%-4i s2=%-4i s3=%-4i s4=%-4i %s\n"         %(' '*4,    s1[1],s2[1],s3[1],s4[1])

        return msg

    def __reprTransient__(self):
        raise NotImplementedError('CBUSH')
        msg = '---BUSH STRAIN---\n'
        msg += '%-8s %6s ' % ('EID', 'eType')
        headers = ['e1', 'e2', 'e3', 'e4', 'Axial', 'eMax', 'eMin']
        for header in headers:
            msg += '%10s ' % header
        msg += '\n'

        for dt, E1s in sorted(self.e1.iteritems()):
            msg += "%s = %g\n" % (self.data_code['name'], dt)
            for eid, e1s in sorted(Els.iteritems()):
                eType = self.eType[eid]
                axial = self.axial[dt][eid]
                e1 = self.e1[dt][eid]
                e2 = self.e2[dt][eid]
                e3 = self.e3[dt][eid]
                e4 = self.e4[dt][eid]
                msg += '%-8i %6s ' % (eid, eType)
                vals = [e1[0], e2[0], e3[0], e4[0], axial]
                for val in vals:
                    if abs(val) < 1e-6:
                        msg += '%10s ' % '0'
                    else:
                        msg += '%10.3g ' % val
                msg += '\n'

                msg += '%s ' % (' ' * 17)
                vals = [e1[1], e2[1], e3[1], e4[1]]
                for val in vals:
                    if isinstance(val, str):
                        msg += '%10s ' % val
                    elif abs(val) < 1e-6:
                        msg += '%10s ' % '0'
                    else:
                        msg += '%10.3g ' % val
                msg += '\n'

                #msg += "eid=%-4s eType=%s s1=%-4i s2=%-4i s3=%-4i s4=%-4i axial=-%5i\n" %(eid,eType,s1[0],s2[0],s3[0],s4[0],axial)
                #msg += "%s                s1=%-4i s2=%-4i s3=%-4i s4=%-4i %s\n"         %(' '*4,    s1[1],s2[1],s3[1],s4[1])

        return msg
