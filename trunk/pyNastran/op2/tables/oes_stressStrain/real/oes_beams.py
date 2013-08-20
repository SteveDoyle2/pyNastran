from __future__ import (nested_scopes, generators, division, absolute_import,
                        print_function, unicode_literals)
import pandas as pd
from numpy import zeros

from .oes_objects import StressObject, StrainObject
from pyNastran.f06.f06_formatting import writeFloats13E

class BeamResultsObject(object):
    def __init__(self):
        self.eType = 'CBEAM'
        self.shape = {}
        self._ncount = 0
        self._inode_start = None
        self._inode_end = None
        self._ielement_start = None
        self._ielement_end = None
        self.data = None
        self.element_data = None

    def getLengthTotal(self):
        return 444  # 44+10*40   (11 nodes)

    def getLength1(self):
        return (44, 'ifffffffff')

    def getLength2(self):
        return (40, 'ifffffffff')

    def _increase_size(self, dt, nelements, nnodes):
        if dt in self.shape:  # default dictionary
            self.shape[dt][0] += nelements
            self.shape[dt][1] += nnodes
        else:
            self.shape[dt] = [nelements, nnodes]
        #print("shape =", self.shape)

    def _get_shape(self):
        ndt = len(self.shape)
        dts = self.shape.keys()
        shape0 = dts[0]
        nelements = self.shape[shape0][0]
        nnodes = self.shape[shape0][1]
        #print("ndt=%s nnodes=%s dts=%s" % (str(ndt), str(nnodes), str(dts)))
        return ndt, nelements, nnodes, dts

    def _increment(self, nnodes, nelements):
        self._inode_start += nnodes
        self._inode_end += nnodes
        self._ielement_start += nelements
        self._ielement_end += nelements
        return self._inode_start, self._inode_end, self._ielement_start, self._ielement_end

    def _preallocate(self, dt, nnodes, nelements):
        ndt, nelements_size, nnodes_size, dts = self._get_shape()
        #print("ndt=%s nelements_size=%s nnodes_size=%s dts=%s" % (ndt, nelements_size, nnodes_size, str(dts)))

        if self._inode_start is not None:
            return (self._inode_start, self._inode_start + nnodes,
                    self._ielement_start, self._ielement_start + nelements)
        #print('----definition----')
        n = ndt * nnodes_size
        if self._ncount != 0:
            asfd
        self._ncount += 1
        self._inode_start = 0
        self._inode_end = nnodes

        self._ielement_start = 0
        self._ielement_end = nelements

        data = {}
        element_data = {}
        columns = []
        if dts[0] is not None:
            name = self.data_code['name']
            if isinstance(dt, int):
                data[name] = pd.Series(zeros((n), dtype='int32'))
            else:
                data[name] = pd.Series(zeros((n), dtype='float32'))
            columns.append(name)

        element_data['element_id'] = pd.Series(zeros((nelements_size), dtype='int32'))
        element_data['element_type'] = pd.Series(zeros(nelements_size, dtype='str'))

        data['element_id'] = pd.Series(zeros((n), dtype='int32'))
        data['grid'] = pd.Series(zeros((n), dtype='float32'))
        data['xxb'] = pd.Series(zeros((n), dtype='float32'))

        #columns.append('element_type')

        #data['grid_type'] = pd.Series(zeros(ndt), dtype='int32'))
        #data['grid_type_str'] = pd.Series(zeros(nnodes), dtype='str'))
        #print('n =', n)

        headers = self._get_headers()
        (ksxc, ksxd, ksxe, ksxf, ksmax, ksmin, mst, msc) = headers

        data[ksxc] = pd.Series(zeros((n), dtype='float32'))
        data[ksxd] = pd.Series(zeros((n), dtype='float32'))
        data[ksxe] = pd.Series(zeros((n), dtype='float32'))
        data[ksxf] = pd.Series(zeros((n), dtype='float32'))

        data[ksmax] = pd.Series(zeros((n), dtype='float32'))
        data[ksmin] = pd.Series(zeros((n), dtype='float32'))
        data['MS_tension'] = pd.Series(zeros((n), dtype='float32'))
        data['MS_compression'] = pd.Series(zeros((n), dtype='float32'))
        # element_type

        columns += ['element_id', 'grid', 'xxb'] + headers

        self.data = pd.DataFrame(data, columns=columns)
        self.element_data = pd.DataFrame(element_data, columns=['element_id', 'element_type'])
        return (self._inode_start, self._inode_end, self._ielement_start, self._ielement_end)

    def _finalize(self, dt):
        ndt, nelements, nnodes, dts = self._get_shape()

        if dt is not None and dt != max(dts):
            return
        #print("----finalize----")

        #grid_type_str = []
        #for grid_type in self.grid_type:
            #grid_type_str.append('C' if grid_type==0 else grid_type)
        #self.grid_type_str = pd.Series(grid_type_str, dtype='str')

        if dts[0] is not None:
            name = self.data_code['name']
            self.data = self.data.set_index([name, 'element_id', 'grid'])
        else:
            self.data = self.data.set_index(['element_id', 'grid'])
        #print("final\n", self.data.to_string())
        del self._inode_start
        del self._inode_end
        #print('---BeamStrainObject---')
        #print(self.data.to_string())

    def get_stats(self):
        msg = self.get_data_code()
        if None not in self.shape:  # transient
            name = self.data_code['name']
            dt_string = name + ', '
            #print("self.data.index[name] = ", self.data.index)
            ntimes = len(getattr(self, name+'s'))
            #ntimes = len(self.data.index[name])
            #s0 = self.smax.keys()[0]
            #nelements = len(self.smax[s0])
            nelements = len(self.data['sxc']) // ntimes // 2
            msg.append('  type=%s n%s=%s nelements=%s\n'
                       % (self.__class__.__name__, name, ntimes, nelements))
        else:
            dt_string = ''
            nelements = len(self.data['smax'])
            msg.append('  real type=%s nelements=%s\n' % (self.__class__.__name__,
                                                     nelements))
        #headers = self._get_headers()
        #(ksxc, ksxd, ksxe, ksxf, ksmax, ksmin, mst, msc) = headers
        msg.append('  data        : index :  %selement_id\n' % dt_string)
        msg.append('              : result:  xxb, %s, %s, %s, %s, %s, %s, %s, %s\n' %(ksxc, ksxd, ksxe, ksxf,
                                                                                      ksxmax, ksxmin, mst, msc))
        #msg.append('  eType, xxb, grids, smax, smin, MS_tension, '
        #           'MS_compression, sxc, sxd, sxe, sxf\n')
        return msg

    #def get_stats(self):
        #nelements = len(self.eType)

        #msg = self.get_data_code()
        #if self.dt is not None:  # transient
            #ntimes = len(self.smax)
            #s0 = self.smax.keys()[0]
            #nelements = len(self.smax[s0])
            #msg.append('  type=%s ntimes=%s nelements=%s\n'
                       #% (self.__class__.__name__, ntimes, nelements))
        #else:
            #nelements = len(self.smax)
            #msg.append('  type=%s nelements=%s\n' % (self.__class__.__name__,
                                                     #nelements))
        #msg.append('  eType, xxb, grids, smax, smin, MS_tension, '
                   #'MS_compression, sxc, sxd, sxe, sxf\n')
        #return msg

    def __repr__(self):
        return self.get_stats()


class BeamStressObject(BeamResultsObject, StressObject):
    """
    ::

      [1,0,0]
                   S T R E S S E S   I N   B E A M   E L E M E N T S        ( C B E A M )
                        STAT DIST/
       ELEMENT-ID  GRID   LENGTH    SXC           SXD           SXE           SXF           S-MAX         S-MIN         M.S.-T   M.S.-C
              1       1   0.000   -3.125000E+04 -3.125000E+04 -3.125000E+04 -3.125000E+04 -3.125000E+04 -3.125000E+04
                      2   1.000   -3.125000E+04 -3.125000E+04 -3.125000E+04 -3.125000E+04 -3.125000E+04 -3.125000E+04
    """
    def __init__(self, data_code, is_sort1, isubcase, dt, read_mode):
        BeamResultsObject.__init__(self)
        StressObject.__init__(self, data_code, isubcase, read_mode)

    def write_f06(self, header, pageStamp, f, pageNum=1, is_mag_phase=False):
        if self.nonlinear_factor is not None:
            return self._write_f06_transient(header, pageStamp, f, pageNum)

        msg = header + ['                                  S T R E S S E S   I N   B E A M   E L E M E N T S        ( C B E A M )\n',
                        '                    STAT DIST/\n',
                        '   ELEMENT-ID  GRID   LENGTH    SXC           SXD           SXE           SXF           S-MAX         S-MIN         M.S.-T   M.S.-C\n']

        for eid in sorted(self.smax):
            msg.append('0  %8i\n' % (eid))
            #print self.xxb[eid]
            for i, nid in enumerate(self.grids[eid]):
                #print i,nid
                xxb = self.xxb[eid][i]
                sxc = self.sxc[eid][i]
                sxd = self.sxd[eid][i]
                sxe = self.sxe[eid][i]
                sxf = self.sxf[eid][i]
                sMax = self.smax[eid][i]
                sMin = self.smin[eid][i]
                SMt = self.MS_tension[eid][i]
                SMc = self.MS_compression[eid][i]
                (vals2, isAllZeros) = writeFloats13E([sxc, sxd, sxe, sxf, sMax, sMin, SMt, SMc])
                (sxc, sxd, sxe, sxf, sMax, sMin, SMt, SMc) = vals2
                msg.append('%19s   %4.3f   %12s %12s %12s %12s %12s %12s %12s %s\n' % (nid, xxb, sxc, sxd, sxe, sxf, sMax, sMin, SMt, SMc.strip()))

        msg.append(pageStamp + str(pageNum) + '\n')
        f.write(''.join(msg))
        return pageNum

    def __write_f06_transient(self, header, pageStamp, f, pageNum=1, is_mag_phase=False):
        words = ['                                  S T R E S S E S   I N   B E A M   E L E M E N T S        ( C B E A M )\n',
                 '                    STAT DIST/\n',
                 '   ELEMENT-ID  GRID   LENGTH    SXC           SXD           SXE           SXF           S-MAX         S-MIN         M.S.-T   M.S.-C\n']
        msg = []
        for dt, SMaxs in sorted(self.smax.iteritems()):
            header[1] = ' %s = %10.4E\n' % (self.data_code['name'], dt)
            msg += header + words
            for eid, Smax in sorted(SMaxs.iteritems()):
                msg.append('0  %8i\n' % (eid))
                for i, nid in enumerate(self.grids[eid]):
                    xxb = self.xxb[eid][i]
                    #sd  = self.sd[eid][i]
                    sxc = self.sxc[dt][eid][i]
                    sxd = self.sxd[dt][eid][i]
                    sxe = self.sxe[dt][eid][i]
                    sxf = self.sxf[dt][eid][i]
                    sMax = self.smax[dt][eid][i]
                    sMin = self.smin[dt][eid][i]
                    SMt = self.MS_tension[dt][eid][i]
                    SMc = self.MS_compression[dt][eid][i]
                    (vals2, isAllZeros) = writeFloats13E([sxc, sxd,
                                                          sxe, sxf, sMax, sMin, SMt, SMc])
                    (sxc, sxd, sxe, sxf, sMax, sMin, SMt, SMc) = vals2
                    msg.append('%19s   %4.3f   %12s %12s %12s %12s %12s %12s %12s %s\n' % (nid, xxb, sxc, sxd, sxe, sxf, sMax, sMin, SMt, SMc.strip()))

            msg.append(pageStamp + str(pageNum) + '\n')
            f.write(''.join(msg))
            msg = ['']
            pageNum += 1
        return pageNum - 1


class BeamStrainObject(BeamResultsObject, StrainObject):
    def __init__(self, data_code, is_sort1, isubcase, dt, read_mode):
        BeamResultsObject.__init__(self)
        StrainObject.__init__(self, data_code, isubcase, read_mode)

    def _get_headers(self):
        return ['sxc', 'sxd', 'sxe', 'sxf', 'smax', 'smin', 'MS_tension', 'MS_compression',]

    def __preallocate(self, dt, nnodes, nelements):
        ndt, nelements_size, nnodes_size, dts = self._get_shape()
        #print("ndt=%s nelements_size=%s nnodes_size=%s dts=%s" % (ndt, nelements_size, nnodes_size, str(dts)))

        if self._inode_start is not None:
            return (self._inode_start, self._inode_start + nnodes,
                    self._ielement_start, self._ielement_start + nelements)
        #print('----definition----')
        n = ndt * nnodes_size
        if self._ncount != 0:
            asfd
        self._ncount += 1
        self._inode_start = 0
        self._inode_end = nnodes

        self._ielement_start = 0
        self._ielement_end = nelements

        data = {}
        element_data = {}
        columns = []
        if dts[0] is not None:
            name = self.data_code['name']
            if isinstance(dt, int):
                data[name] = pd.Series(zeros((n), dtype='int32'))
            else:
                data[name] = pd.Series(zeros((n), dtype='float32'))
            columns.append(name)

        element_data['element_id'] = pd.Series(zeros((nelements_size), dtype='int32'))
        element_data['element_type'] = pd.Series(zeros(nelements_size, dtype='str'))

        data['element_id'] = pd.Series(zeros((n), dtype='int32'))
        data['grid'] = pd.Series(zeros((n), dtype='float32'))
        data['xxb'] = pd.Series(zeros((n), dtype='float32'))

        #columns.append('element_type')

        #data['grid_type'] = pd.Series(zeros(ndt), dtype='int32'))
        #data['grid_type_str'] = pd.Series(zeros(nnodes), dtype='str'))
        #print('n =', n)

        data['exc'] = pd.Series(zeros((n), dtype='float32'))
        data['exd'] = pd.Series(zeros((n), dtype='float32'))
        data['exe'] = pd.Series(zeros((n), dtype='float32'))
        data['exf'] = pd.Series(zeros((n), dtype='float32'))

        data['emax'] = pd.Series(zeros((n), dtype='float32'))
        data['emin'] = pd.Series(zeros((n), dtype='float32'))
        data['MS_tension'] = pd.Series(zeros((n), dtype='float32'))
        data['MS_compression'] = pd.Series(zeros((n), dtype='float32'))
        # element_type

        headers = self._get_headers()
        columns += ['element_id', 'grid', 'xxb'] + headers

        self.data = pd.DataFrame(data, columns=columns)
        self.element_data = pd.DataFrame(element_data, columns=['element_id', 'element_type'])
        return (self._inode_start, self._inode_end, self._ielement_start, self._ielement_end)

    def _get_headers(self):
        return ['exc', 'exd', 'exe', 'exf', 'emax', 'emin', 'MS_tension', 'MS_compression',]

    def write_f06(self, header, pageStamp, f, pageNum=1, is_mag_phase=False):
        if self.nonlinear_factor is not None:
            return self._write_f06_transient(header, pageStamp, f, pageNum)

        msg = header + ['                                  S T R A I N S   I N   B E A M   E L E M E N T S        ( C B E A M )\n',
                        '                    STAT DIST/\n',
                        '   ELEMENT-ID  GRID   LENGTH    SXC           SXD           SXE           SXF           S-MAX         S-MIN         M.S.-T   M.S.-C\n']

        for eid in sorted(self.smax):
            msg.append('0  %8i\n' % (eid))
            #print self.xxb[eid]
            for i, nid in enumerate(self.grids[eid]):
                xxb = self.xxb[eid][i]
                sxc = self.sxc[eid][i]
                sxd = self.sxd[eid][i]
                sxe = self.sxe[eid][i]
                sxf = self.sxf[eid][i]
                sMax = self.smax[eid][i]
                sMin = self.smin[eid][i]
                SMt = self.MS_tension[eid][i]
                SMc = self.MS_compression[eid][i]
                (vals2, isAllZeros) = writeFloats13E([
                    sxc, sxd, sxe, sxf, sMax, sMin, SMt, SMc])
                (sxc, sxd, sxe, sxf, sMax, sMin, SMt, SMc) = vals2
                msg.append('%19s   %4.3f   %12s %12s %12s %12s %12s %12s %12s %s\n' % (nid, xxb, sxc, sxd, sxe, sxf, sMax, sMin, SMt, SMc.strip()))

        msg.append(pageStamp + str(pageNum) + '\n')
        f.write(''.join(msg))
        return pageNum

    def _write_f06_transient(self, header, pageStamp, f, pageNum=1, is_mag_phase=False):
        words = ['                                  S T R A I N S   I N   B E A M   E L E M E N T S        ( C B E A M )\n',
                 '                    STAT DIST/\n',
                 '   ELEMENT-ID  GRID   LENGTH    SXC           SXD           SXE           SXF           S-MAX         S-MIN         M.S.-T   M.S.-C\n']
        msg = []
        for dt, SMaxs in sorted(self.smax.iteritems()):
            header[1] = ' %s = %10.4E\n' % (self.data_code['name'], dt)
            msg += header + words
            for eid, Smax in sorted(SMaxs.iteritems()):
                msg.append('0  %8i\n' % (eid))
                for i, nid in enumerate(self.grids[eid]):
                    xxb = self.xxb[eid][i]
                    sxc = self.sxc[dt][eid][i]
                    sxd = self.sxd[dt][eid][i]
                    sxe = self.sxe[dt][eid][i]
                    sxf = self.sxf[dt][eid][i]
                    sMax = self.smax[dt][eid][i]
                    sMin = self.smin[dt][eid][i]
                    SMt = self.MS_tension[dt][eid][i]
                    SMc = self.MS_compression[dt][eid][i]
                    (vals2, isAllZeros) = writeFloats13E([sxc, sxd,
                                                          sxe, sxf, sMax, sMin, SMt, SMc])
                    (sxc, sxd, sxe, sxf, sMax, sMin, SMt, SMc) = vals2
                    msg.append('%19s   %4.3f   %12s %12s %12s %12s %12s %12s %12s %s\n' % (nid, xxb, sxc, sxd, sxe, sxf, sMax, sMin, SMt, SMc.strip()))

            msg.append(pageStamp + str(pageNum) + '\n')
            f.write(''.join(msg))
            msg = ['']
            pageNum += 1
        return pageNum - 1