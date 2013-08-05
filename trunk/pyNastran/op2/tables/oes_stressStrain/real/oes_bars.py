from __future__ import (nested_scopes, generators, division, absolute_import,
                        print_function, unicode_literals)
import pandas as pd
from numpy import zeros

from pyNastran.utils import is_string
from .oes_objects import StressObject, StrainObject
from pyNastran.f06.f06_formatting import writeFloats13E, writeImagFloats13E


class RealBarResults(object):
    def __init__(self):
        self.shape = {}
        self._ncount = 0
        self._inode_start = None
        self._inode_end = None
        self._ielement_start = None
        self._ielement_end = None
        self.data = None
        self.element_data = None

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

        n = ndt * nnodes_size
        if self._ncount != 0:
            asfd
        self._ncount += 1
        self._inode_start = 0
        self._inode_end = nnodes

        self._ielement_start = 0
        self._ielement_end = nelements

        data = {}
        #element_data = {}
        columns = []
        if dts[0] is not None:
            name = self.data_code['name']
            if isinstance(dt, int):
                data[name] = pd.Series(zeros((n), dtype='int32'))
            else:
                data[name] = pd.Series(zeros((n), dtype='float32'))
            columns.append(name)

        #element_data['element_type'] = pd.Series(zeros(nelements_size, dtype='str'))
        #element_data['element_id'] = pd.Series(zeros((nelements_size), dtype='int32'))

        data['element_id'] = pd.Series(zeros((n), dtype='int32'))
        #data['element_type'] = pd.Series(zeros((n), dtype='str'))
        #data['node_id'] = pd.Series(zeros((n), dtype='int32'))

        #columns.append('element_type')
        #columns.append('element_id')
        #columns.append('node_id')

        #data['grid_type'] = pd.Series(zeros(ndt), dtype='int32'))
        #data['grid_type_str'] = pd.Series(zeros(nnodes), dtype='str'))
        #print('n =', n)

        headers = self._get_headers()
        (s1a, s2a, s3a, s4a, s1b, s2b, s3b, s4b,
              smaxa, smina, smaxb, sminb) = headers

        data[s1a] = pd.Series(zeros((n), dtype='float32'))
        data[s2a] = pd.Series(zeros((n), dtype='float32'))
        data[s3a] = pd.Series(zeros((n), dtype='float32'))
        data[s4a] = pd.Series(zeros((n), dtype='float32'))
        data[smaxa] = pd.Series(zeros((n), dtype='float32'))
        data[smina] = pd.Series(zeros((n), dtype='float32'))

        data[s1b] = pd.Series(zeros((n), dtype='float32'))
        data[s2b] = pd.Series(zeros((n), dtype='float32'))
        data[s3b] = pd.Series(zeros((n), dtype='float32'))
        data[s4b] = pd.Series(zeros((n), dtype='float32'))
        data[smaxb] = pd.Series(zeros((n), dtype='float32'))
        data[sminb] = pd.Series(zeros((n), dtype='float32'))

        data['axial'] = pd.Series(zeros((n), dtype='float32'))
        data['MS_tension'] = pd.Series(zeros((n), dtype='float32'))
        data['MS_compression'] = pd.Series(zeros((n), dtype='float32'))
        # element_type

        columns += ['element_id', s1a, s2a, s3a, s4a, s1b, s2b, s3b, s4b,
                    'axial', smaxa, smina, smaxb, sminb,
                    'MS_tension', 'MS_compression']

        self.data = pd.DataFrame(data, columns=columns)
        #self.element_data = pd.DataFrame(element_data, columns=['element_id', 'element_type', 'cid'])
        return (self._inode_start, self._inode_end, self._ielement_start, self._ielement_end)

    def _finalize(self, dt):
        ndt, nelements, nnodes, dts = self._get_shape()

        if dt != max(dts):
            return
        print("----finalize----")

        #grid_type_str = []
        #for grid_type in self.grid_type:
            #grid_type_str.append('C' if grid_type==0 else grid_type)
        #self.grid_type_str = pd.Series(grid_type_str, dtype='str')

        if dts[0] is not None:
            name = self.data_code['name']
            self.data = self.data.set_index([name, 'element_id'])
        else:
            self.data = self.data.set_index('element_id')
        #print "final\n", self.data
        del self._inode_start
        del self._inode_end
        #print("---BarStressObject---")
        #print(self.data.to_string())

    def get_stats(self):
        nelements = len(self.data['axial'])
        msg = self.get_data_code()
        if self.nonlinear_factor is not None:  # transient
            ntimes = len(self.shape)
            msg.append('  real type=%s ntimes=%s nelements=%s\n'
                       % (self.__class__.__name__, ntimes, nelements // ntimes))
        else:
            msg.append('  real type=%s nelements=%s\n' % (self.__class__.__name__,
                                                          nelements))
        headers = self._get_headers()
        (s1a, s2a, s3a, s4a,
         s1b, s2b, s3b, s4b,
         smaxa, smina, smaxb, sminb) = headers

        msg.append('  element data: index :  element_id\n')
        msg.append('              : result:  element_type\n')
        msg.append('  data        : index :  element_id\n')
        msg.append('              : result:  %s, %s, %s, %s, %s, %s,\n' % (s1a, s2a, s3a, s4a, smaxa, smina) )
        msg.append('              : result:  %s, %s, %s, %s, %s, %s,\n' % (s1b, s2b, s3b, s4b, smaxb, sminb) )
        msg.append('                         axial, MS_tension, MS_compression\n')
        return msg

    def __repr__(self):
        return self.get_stats()


class BarStressObject(StressObject, RealBarResults):
    """
    ::
      # s_code=0
                                 S T R E S S E S   I N   B A R   E L E M E N T S          ( C B A R )
      ELEMENT        SA1            SA2            SA3            SA4           AXIAL          SA-MAX         SA-MIN     M.S.-T
        ID.          SB1            SB2            SB3            SB4           STRESS         SB-MAX         SB-MIN     M.S.-C
    """
    def __init__(self, data_code, is_sort1, isubcase, dt, read_mode):
        StressObject.__init__(self, data_code, isubcase, read_mode)
        RealBarResults.__init__(self)

    def add_f06_data(self, data, transient):
        if transient is None:
            for line in data:
                (eType, eid, s1A, s2A, s3A, s4A, axialA, smaxA, sminA, MSt,
                 s1B, s2B, s3B, s4B, smaxB, sminB, MSc) = line
                self.eType[eid] = 'CBAR'
                self.s1[eid] = [s1A, s1B]
                self.s2[eid] = [s2A, s2B]
                self.s3[eid] = [s3A, s3B]
                self.s4[eid] = [s4A, s4B]

                self.axial[eid] = axialA
                self.smax[eid] = [smaxA, smaxB]
                self.smin[eid] = [sminA, sminB]
                #self.MS_tension[eid]     = MSt
                #self.MS_compression[eid] = MSc
            return

        (dtName, dt) = transient
        self.data_code['name'] = dtName
        #print "dt = ",dt
        #print "dtName = ",dtName
        if dt not in self.s1:
            self.update_dt(self.data_code, dt)
            self.isTransient = True

        for line in data:
            (eType, eid, s1A, s2A, s3A, s4A, axialA, smaxA, sminA, MSt,
             s1B, s2B, s3B, s4B, smaxB, sminB, MSc) = line
            self.eType[eid] = 'CBAR'
            self.s1[dt][eid] = [s1A, s1B]
            self.s2[dt][eid] = [s2A, s2B]
            self.s3[dt][eid] = [s3A, s3B]
            self.s4[dt][eid] = [s4A, s4B]

            self.axial[dt][eid] = axialA
            self.smax[dt][eid] = [smaxA, smaxB]
            self.smin[dt][eid] = [sminA, sminB]
            #self.MS_tension[dt][eid]     = MSt
            #self.MS_compression[dt][eid] = MSc

    def getLength(self):
        return (68, 'iffffffffffffffff')

    def _get_headers(self):
        return ('s1a', 's2a', 's3a', 's4a', 's1b', 's2b', 's3b', 's4b',
                'smaxa', 'smina', 'smaxb', 'sminb')

    def write_f06(self, header, pageStamp, f, pageNum=1, is_mag_phase=False):
        if self.nonlinear_factor is not None:
            return self._write_f06_transient(header, pageStamp, f, pageNum)

        msg = header + [
                '                                 S T R E S S E S   I N   B A R   E L E M E N T S          ( C B A R )\n',
                '  ELEMENT        SA1            SA2            SA3            SA4           AXIAL          SA-MAX         SA-MIN     M.S.-T\n',
                '    ID.          SB1            SB2            SB3            SB4           STRESS         SB-MAX         SB-MIN     M.S.-C\n',
              ]
        for eid, S1s in sorted(self.s1.iteritems()):
            #eType = self.eType[eid]
            axial = self.axial[eid]
            #MSt = self.MSt[eid]
            #MSc = self.MSc[eid]
            MSt = ''
            MSc = ''

            s1 = self.s1[eid]
            s2 = self.s2[eid]
            s3 = self.s3[eid]
            s4 = self.s4[eid]
            smax = self.smax[eid]
            smin = self.smin[eid]
            vals = [s1[0], s2[0], s3[0], s4[0], axial, smax[0], smin[0],
                    s1[1], s2[1], s3[1], s4[1], smax[1], smin[1]]
            (vals2, isAllZeros) = writeFloats13E(vals)
            [s1a, s2a, s3a, s4a, axial, smaxa, smina,
             s1b, s2b, s3b, s4b, smaxb, sminb] = vals2
            msg.append('0%8i   %13s  %13s  %13s  %13s  %13s  %13s  %13s %-s\n' % (eid, s1a, s2a, s3a, s4a, axial, smaxa, smina, MSt.rstrip()))
            msg.append(' %8s   %13s  %13s  %13s  %13s  %13s  %13s  %13s %-s\n' % ('', s1b, s2b, s3b, s4b, '', smaxb, sminb, MSc.rstrip()))

        msg.append(pageStamp + str(pageNum) + '\n')
        return (''.join(msg), pageNum)

    def _write_f06_transient(self, header, pageStamp, f, pageNum=1, is_mag_phase=False):
        words = [
                '                                 S T R E S S E S   I N   B A R   E L E M E N T S          ( C B A R )\n',
                '  ELEMENT        SA1            SA2            SA3            SA4           AXIAL          SA-MAX         SA-MIN     M.S.-T\n',
                '    ID.          SB1            SB2            SB3            SB4           STRESS         SB-MAX         SB-MIN     M.S.-C\n',
              ]
        msg = []
        for dt, S1s in sorted(self.s1.iteritems()):
            header[1] = ' %s = %10.4E\n' % (self.data_code['name'], dt)
            msg += header + words
            for eid, S1 in sorted(S1s.iteritems()):
                #eType = self.eType[eid]
                axial = self.axial[dt][eid]
                #MSt = self.MSt[eid]
                #MSc = self.MSc[eid]
                MSt = ''
                MSc = ''

                s1 = self.s1[dt][eid]
                s2 = self.s2[dt][eid]
                s3 = self.s3[dt][eid]
                s4 = self.s4[dt][eid]
                smax = self.smax[dt][eid]
                smin = self.smin[dt][eid]
                vals = [s1[0], s2[0], s3[0], s4[0], axial, smax[0], smin[0],
                        s1[1], s2[1], s3[1], s4[1], smax[1], smin[1]]
                (vals2, isAllZeros) = writeFloats13E(vals)
                [s1a, s2a, s3a, s4a, axial, smaxa, smina,
                 s1b, s2b, s3b, s4b, smaxb, sminb] = vals2
                msg.append('0%8i   %13s  %13s  %13s  %13s  %13s  %13s  %13s %-s\n' % (eid, s1a, s2a, s3a, s4a, axial, smaxa, smina, MSt.rstrip()))
                msg.append(' %8s   %13s  %13s  %13s  %13s  %13s  %13s  %13s %-s\n' % ('', s1b, s2b, s3b, s4b, '', smaxb, sminb, MSc.rstrip()))

            msg.append(pageStamp + str(pageNum) + '\n')
            pageNum += 1
            f.write(''.join(msg))
            msg = ['']
        return pageNum - 1


class BarStrainObject(StrainObject, RealBarResults):
    """
    ::

      # s_code=10
                                       S T R A I N S   I N   B A R   E L E M E N T S          ( C B A R )
      ELEMENT        SA1            SA2            SA3            SA4           AXIAL          SA-MAX         SA-MIN     M.S.-T
        ID.          SB1            SB2            SB3            SB4           STRAIN         SB-MAX         SB-MIN     M.S.-C
    """
    def __init__(self, data_code, is_sort1, isubcase, dt, read_mode):
        StrainObject.__init__(self, data_code, isubcase, read_mode)
        RealBarResults.__init__(self)

    def _get_headers(self):
        return ('e1a', 'e2a', 'e3a', 'e4a', 'e1b', 'e2b', 'e3b', 'e4b',
                'emaxa', 'emina', 'emaxb', 'eminb')

    def add_f06_data(self, data, transient):
        if transient is None:
            for line in data:
                (eType, eid, e1A, e2A, e3A, e4A, axialA, emaxA, eminA, MSt,
                           e1B, e2B, e3B, e4B, emaxB, eminB, MSc) = line
                self.eType[eid] = 'CBAR'
                self.e1[eid] = [e1A, e1B]
                self.e2[eid] = [e2A, e2B]
                self.e3[eid] = [e3A, e3B]
                self.e4[eid] = [e4A, e4B]

                self.axial[eid] = axialA
                self.emax[eid] = [emaxA, emaxB]
                self.emin[eid] = [eminA, eminB]
                #self.MS_tension[eid]     = MSt
                #self.MS_compression[eid] = MSc
            return

        (dtName, dt) = transient
        self.data_code['name'] = dtName
        if dt not in self.s1:
            self.update_dt(self.data_code, dt)
            self.isTransient = True

        for line in data:
            (eType, eid, e1A, e2A, e3A, e4A, axialA, emaxA, eminA, MSt,
                       e1B, e2B, e3B, e4B, emaxB, eminB, MSc) = line
            self.eType[eid] = 'CBAR'
            self.e1[dt][eid] = [e1A, e1B]
            self.e2[dt][eid] = [e2A, e2B]
            self.e3[dt][eid] = [e3A, e3B]
            self.e4[dt][eid] = [e4A, e4B]

            self.axial[dt][eid] = axialA
            self.emax[dt][eid] = [emaxA, emaxB]
            self.emin[dt][eid] = [eminA, eminB]
            #self.MS_tension[dt][eid]     = MSt
            #self.MS_compression[dt][eid] = MSc

    def get_transients(self):
        k = self.e1.keys()
        k.sort()
        return k

    def write_f06(self, header, pageStamp, f, pageNum=1, is_mag_phase=False):
        if self.isTransient:
            return self._write_f06_transient(header, pageStamp, f, pageNum)

        msg = header + [
                '                                  S T R A I N S    I N   B A R   E L E M E N T S          ( C B A R )\n',
                '  ELEMENT        SA1            SA2            SA3            SA4           AXIAL          SA-MAX         SA-MIN     M.S.-T\n',
                '    ID.          SB1            SB2            SB3            SB4           STRAIN         SB-MAX         SB-MIN     M.S.-C\n',
              ]
        for eid, E1s in sorted(self.e1.iteritems()):
            #eType = self.eType[eid]
            axial = self.axial[eid]
            #MSt = self.MSt[eid]
            #MSc = self.MSc[eid]
            MSt = ''
            MSc = ''

            e1 = self.e1[eid]
            e2 = self.e2[eid]
            e3 = self.e3[eid]
            e4 = self.e4[eid]
            emax = self.emax[eid]
            emin = self.emin[eid]
            vals = [e1[0], e2[0], e3[0], e4[0], axial, emax[0], emin[0],
                    e1[1], e2[1], e3[1], e4[1], emax[1], emin[1]]
            (vals2, isAllZeros) = writeFloats13E(vals)
            [e10, e20, e30, e40, axial, emax0, emin0,
             e11, e21, e31, e41, emax1, emin1] = vals2

            msg.append('0%8i   %13s  %13s  %13s  %13s  %13s  %13s  %13s %-s\n' % (eid, e10, e20, e30, e40, axial, emax0, emin0, MSt.rstrip()))
            msg.append(' %8s   %13s  %13s  %13s  %13s  %13s  %13s  %13s %-s\n' % ('', e11, e21, e31, e41, '', emax1, emin1, MSc.rstrip()))

        msg.append(pageStamp + str(pageNum) + '\n')
        f.write(''.join(msg))
        return pageNum

    def _write_f06_transient(self, header, pageStamp, f, pageNum=1, is_mag_phase=False):
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
                #eType = self.eType[eid]
                axial = self.axial[eid]
                #MSt = self.MSt[eid]
                #MSc = self.MSc[eid]
                MSt = ''
                MSc = ''

                e1 = self.e1[eid]
                e2 = self.e2[eid]
                e3 = self.e3[eid]
                e4 = self.e4[eid]
                emax = self.emax[eid]
                emin = self.emin[eid]
                vals = [e1[0], e2[0], e3[0], e4[0], axial, emax[0], emin[0],
                        e1[1], e2[1], e3[1], e4[1], emax[1], emin[1]]
                (vals2, isAllZeros) = writeFloats13E(vals)
                [e10, e20, e30, e40, axial, emax0, emin0,
                 e11, e21, e31, e41, emax1, emin1] = vals2

                msg.append('0%8i   %13s  %13s  %13s  %13s  %13s  %13s  %13s %-s\n' % (eid, e10, e20, e30, e40, axial, emax0, emin0, MSt.rstrip()))
                msg.append(' %8s   %13s  %13s  %13s  %13s  %13s  %13s  %13s %-s\n' % ('', e11, e21, e31, e41, '', emax1, emin1, MSc.rstrip()))
            msg.append(pageStamp + str(pageNum) + '\n')
            f.write(''.join(msg))
            msg = ['']
            pageNum += 1
        return pageNum - 1