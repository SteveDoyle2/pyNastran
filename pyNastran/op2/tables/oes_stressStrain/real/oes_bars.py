from __future__ import (nested_scopes, generators, division, absolute_import,
                        print_function, unicode_literals)
from six import iteritems
from itertools import count
from numpy import zeros, searchsorted, ravel

from pyNastran.op2.tables.oes_stressStrain.real.oes_objects import StressObject, StrainObject, OES_Object
from pyNastran.f06.f06_formatting import writeFloats13E, _eigenvalue_header


class RealBarArray(OES_Object):
    def __init__(self, data_code, is_sort1, isubcase, dt):
        OES_Object.__init__(self, data_code, isubcase, apply_data_code=False)
        self.eType = {}
        #self.code = [self.format_code, self.sort_code, self.s_code]

        #self.ntimes = 0  # or frequency/mode
        #self.ntotal = 0
        self.ielement = 0
        self.nelements = 0  # result specific

        if is_sort1:
            if dt is not None:
                #self.add = self.add_sort1
                self.add_new_eid = self.add_new_eid_sort1
                #self.addNewNode = self.addNewNodeSort1
        else:
            raise NotImplementedError('SORT2')
            #assert dt is not None
            #self.add = self.addSort2
            #self.add_new_eid = self.add_new_eid_sort2
            #self.addNewNode = self.addNewNodeSort2

    def is_real(self):
        return True

    def is_complex(self):
        return False

    def _reset_indices(self):
        self.itotal = 0
        self.ielement = 0

    def _get_msgs(self):
        raise NotImplementedError('%s needs to implement _get_msgs' % self.__class__.__name__)

    def get_headers(self):
        raise NotImplementedError('%s needs to implement get_headers' % self.__class__.__name__)

    def build(self):
        if self.is_built:
            return
        #print("self.ielement =", self.ielement)
        #print('ntimes=%s nelements=%s ntotal=%s' % (self.ntimes, self.nelements, self.ntotal))

        assert self.ntimes > 0, 'ntimes=%s' % self.ntimes
        assert self.nelements > 0, 'nelements=%s' % self.nelements
        assert self.ntotal > 0, 'ntotal=%s' % self.ntotal

        if self.element_type == 34:
            nnodes_per_element = 1
        else:
            raise NotImplementedError(self.element_type)

        self.itime = 0
        self.ielement = 0
        self.itotal = 0
        #self.ntimes = 0
        #self.nelements = 0
        self.is_built = True

        #print("***name=%s type=%s nnodes_per_element=%s ntimes=%s nelements=%s ntotal=%s" % (
            #self.element_name, self.element_type, nnodes_per_element, self.ntimes, self.nelements, self.ntotal))
        dtype = 'float32'
        if isinstance(self.nonlinear_factor, int):
            dtype = 'int32'
        self._times = zeros(self.ntimes, dtype=dtype)
        self.element = zeros(self.ntotal, dtype='int32')

        #[s1a, s2a, s3a, s4a, axial, smaxa, smina, MS_tension,
        # s1b, s2b, s3b, s4b,        sminb, sminb, MS_compression]
        self.data = zeros((self.ntimes, self.ntotal, 15), dtype='float32')

    def add_new_eid(self, eType, dt, eid, s1a, s2a, s3a, s4a, axial, smaxa, smina, MSt,
                    s1b, s2b, s3b, s4b, smaxb, sminb, MSc):
        self.add_new_eid_sort1(eType, dt, eid,
                               s1a, s2a, s3a, s4a, axial, smaxa, smina, MSt,
                               s1b, s2b, s3b, s4b, smaxb, sminb, MSc)

    def add_new_eid_sort1(self, eType, dt, eid,
                          s1a, s2a, s3a, s4a, axial, smaxa, smina, MSt,
                          s1b, s2b, s3b, s4b, smaxb, sminb, MSc):

        assert isinstance(eid, int)
        self._times[self.itime] = dt
        self.element[self.itotal] = eid
        self.data[self.itime, self.itotal, :] = [s1a, s2a, s3a, s4a, axial, smaxa, smina, MSt,
                                                 s1b, s2b, s3b, s4b,        smaxb, sminb, MSc]
        self.itotal += 1
        self.ielement += 1

    #def add_sort1(self, dt, eid, nodeID, fd, oxx, oyy, txy, angle, majorP, minorP, ovm):
        #assert eid is not None
        #msg = "i=%s dt=%s eid=%s nodeID=%s fd=%g oxx=%g oyy=%g \ntxy=%g angle=%g major=%g minor=%g ovmShear=%g" % (
            #self.itotal, dt, eid, nodeID, fd, oxx, oyy, txy, angle, majorP, minorP, ovm)
        ##print(msg)
        #if isinstance(nodeID, string_types):
            #nodeID = 0
        ##assert isinstance(nodeID, int), nodeID
        #self.element_node[self.itotal, :] = [eid, nodeID]
        #self.data[self.itime, self.itotal, :] = [fd, oxx, oyy, txy, angle, majorP, minorP, ovm]
        #self.itotal += 1

    def get_stats(self):
        if not self.is_built:
            return ['<%s>\n' % self.__class__.__name__,
                    '  ntimes: %i\n' % self.ntimes,
                    '  ntotal: %i\n' % self.ntotal,
                    ]

        nelements = self.ntotal
        ntimes = self.ntimes
        ntotal = self.ntotal
        nelements = self.ntotal

        msg = []
        if self.nonlinear_factor is not None:  # transient
            msg.append('  type=%s ntimes=%i nelements=%i\n'
                       % (self.__class__.__name__, ntimes, nelements))
            ntimes_word = 'ntimes'
        else:
            msg.append('  type=%s nelements=%i\n'
                       % (self.__class__.__name__, nelements))
            ntimes_word = 1
        headers = self.get_headers()

        n = len(headers)
        assert n == self.data.shape[2], 'nheaders=%s shape=%s' % (n, str(self.data.shape))
        msg.append('  data: [%s, ntotal, %i] where %i=[%s]\n' % (ntimes_word, n, n, str(', '.join(headers))))
        msg.append('  data.shape = %s\n' % str(self.data.shape).replace('L', ''))
        msg.append('  element type: %s\n  ' % self.element_name)
        msg += self.get_data_code()
        return msg

    def get_element_index(self, eids):
        # elements are always sorted; nodes are not
        itot = searchsorted(eids, self.element)  #[0]
        return itot

    def eid_to_element_node_index(self, eids):
        ind = ravel([searchsorted(self.element == eid) for eid in eids])
        #ind = searchsorted(eids, self.element)
        #ind = ind.reshape(ind.size)
        #ind.sort()
        return ind

    def write_f06(self, header, page_stamp, page_num=1, f=None, is_mag_phase=False, is_sort1=True):
        msg = self._get_msgs()
        (ntimes, ntotal) = self.data.shape[:2]
        eids = self.element
        #print('CBAR ntimes=%s ntotal=%s' % (ntimes, ntotal))
        for itime in range(ntimes):
            dt = self._times[itime]
            header = _eigenvalue_header(self, header, itime, ntimes, dt)
            f.write(''.join(header + msg))

            s1a = self.data[itime, :, 0]
            s2a = self.data[itime, :, 1]
            s3a = self.data[itime, :, 2]
            s4a = self.data[itime, :, 3]

            axial = self.data[itime, :, 4]
            smaxa = self.data[itime, :, 5]
            smina = self.data[itime, :, 6]
            MSt = self.data[itime, :, 7]

            s1b = self.data[itime, :, 8]
            s2b = self.data[itime, :, 9]
            s3b = self.data[itime, :, 10]
            s4b = self.data[itime, :, 11]

            smaxb = self.data[itime, :, 12]
            sminb = self.data[itime, :, 13]
            MSc = self.data[itime, :, 14]

            # loop over all the elements
            for (i, eid, s1ai, s2ai, s3ai, s4ai, axiali, smaxai, sminai, MSti,
                         s1bi, s2bi, s3bi, s4bi,         smaxbi, sminbi, MSci) in zip(
                count(), eids, s1a, s2a, s3a, s4a, axial, smaxa, smina, MSt,
                               s1b, s2b, s3b, s4b,        smaxb, sminb, MSc):

                vals = [s1ai, s2ai, s3ai, s4ai, axiali, smaxai, sminai, MSti,
                        s1bi, s2bi, s3bi, s4bi,         smaxbi, sminbi, MSci]
                (vals2, is_all_zeros) = writeFloats13E(vals)
                [s1ai, s2ai, s3ai, s4ai, axiali, smaxai, sminai, MSti,
                 s1bi, s2bi, s3bi, s4bi,         smaxbi, sminbi, MSci] = vals2
                f.write('0%8i   %-13s  %-13s  %-13s  %-13s  %-13s  %-13s  %-13s %s\n'
                        ' %8s   %-13s  %-13s  %-13s  %-13s  %-13s  %-13s  %-13s %s\n'
                        % (eid, s1ai, s2ai, s3ai, s4ai, axiali, smaxai, sminai, MSti,
                            '', s1bi, s2bi, s3bi, s4bi, '',     smaxbi, sminbi, MSci))

            f.write(page_stamp % page_num)
            page_num += 1

        if self.nonlinear_factor is None:
            page_num -= 1
        return page_num


class RealBarStressArray(RealBarArray, StressObject):
    def __init__(self, data_code, is_sort1, isubcase, dt):
        RealBarArray.__init__(self, data_code, is_sort1, isubcase, dt)
        StressObject.__init__(self, data_code, isubcase)

    def get_headers(self):
        headers = ['s1a', 's2a', 's3a', 's4a', 'axial', 'smaxa', 'smina', 'MS_tension',
                   's1b', 's2b', 's3b', 's4b',          'smaxb', 'sminb', 'MS_compression']
        return headers

    def _get_msgs(self):
        if self.element_type == 34:
            pass
        else:
            raise NotImplementedError(self.element_type)

        msg = [
            '                                 S T R E S S E S   I N   B A R   E L E M E N T S          ( C B A R )\n',
            '  ELEMENT        SA1            SA2            SA3            SA4           AXIAL          SA-MAX         SA-MIN     M.S.-T\n',
            '    ID.          SB1            SB2            SB3            SB4           STRESS         SB-MAX         SB-MIN     M.S.-C\n',
        ]
        return msg


class RealBarStrainArray(RealBarArray, StrainObject):
    def __init__(self, data_code, is_sort1, isubcase, dt):
        RealBarArray.__init__(self, data_code, is_sort1, isubcase, dt)
        StrainObject.__init__(self, data_code, isubcase)

    def get_headers(self):
        headers = ['e1a', 'e2a', 'e3a', 'e4a', 'axial', 'emaxa', 'emina', 'MS_tension',
                   'e1b', 'e2b', 'e3b', 'e4b',          'emaxb', 'eminb', 'MS_compression']
        return headers

    def _get_msgs(self):
        if self.element_type == 34:
            pass
        else:
            raise NotImplementedError(self.element_type)

        msg = [
            '                                  S T R A I N S    I N   B A R   E L E M E N T S          ( C B A R )\n',
            '  ELEMENT        SA1            SA2            SA3            SA4           AXIAL          SA-MAX         SA-MIN     M.S.-T\n',
            '    ID.          SB1            SB2            SB3            SB4           STRAIN         SB-MAX         SB-MIN     M.S.-C\n',
        ]
        return msg


class RealBarStress(StressObject):
    """
    ::

      # s_code=0
                                 S T R E S S E S   I N   B A R   E L E M E N T S          ( C B A R )
      ELEMENT        SA1            SA2            SA3            SA4           AXIAL          SA-MAX         SA-MIN     M.S.-T
        ID.          SB1            SB2            SB3            SB4           STRESS         SB-MAX         SB-MIN     M.S.-C
    """
    def __init__(self, data_code, is_sort1, isubcase, dt):
        StressObject.__init__(self, data_code, isubcase)
        self.eType = {}

        self.code = [self.format_code, self.sort_code, self.s_code]
        self.s1 = {}
        self.s2 = {}
        self.s3 = {}
        self.s4 = {}
        self.axial = {}
        self.smax = {}
        self.smin = {}
        self.MS_tension = {}
        self.MS_compression = {}

        #if self.element_type==100:
            #self.getLength = self.getLength100_format1_sort0
            #self.add_new_eid = self.addNewEid100

        self.dt = dt
        #print "BAR dt=%s" %(dt)
        if is_sort1:
            if dt is not None:
                #self.add = self.add_sort1
                self.add_new_eid = self.add_new_eid_sort1
        else:
            assert dt is not None
            #self.add = self.addSort2
            #self.add_new_eid = self.add_new_eid_sort2

    def get_stats(self):
        nelements = len(self.eType)

        msg = self.get_data_code()
        if self.dt is not None:  # transient
            ntimes = len(self.s1)
            msg.append('  type=%s ntimes=%s nelements=%s\n'
                       % (self.__class__.__name__, ntimes, nelements))
        else:
            msg.append('  type=%s nelements=%s\n' % (self.__class__.__name__,
                                                     nelements))
        msg.append('  eType, s1, s2, s3, s4, axial, smax, smin, '
                   'MS_tension, MS_compression\n')
        return msg

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

    def delete_transient(self, dt):
        del self.s1[dt]
        del self.s2[dt]
        del self.s3[dt]
        del self.s4[dt]
        del self.axial[dt]
        del self.smax[dt]
        del self.smin[dt]

    def get_transients(self):
        k = self.s1.keys()
        k.sort()
        return k

    def add_new_transient(self, dt):
        """
        initializes the transient variables
        """
        #print "****add new transient****"
        self.dt = dt
        self.s1[dt] = {}
        self.s2[dt] = {}
        self.s3[dt] = {}
        self.s4[dt] = {}
        self.axial[dt] = {}
        self.smax[dt] = {}
        self.smin[dt] = {}
        #self.MS_tension[dt]     = {}
        #self.MS_compression[dt] = {}

    def addNewEid100(self, dt, out):
        #print "out = ",out
        #return
        (eid, s1, s2, s3, s4, axial, smax, smin, MSt, MSc) = out
        #print "Bar Stress add..."
        self.eType[eid] = 'CBAR'  # eType

        if eid in self.s1:
            self.s1[eid].append(s1)
            self.s2[eid].append(s2)
            self.s3[eid].append(s3)
            self.s4[eid].append(s4)
            self.axial[eid].append(axial)
            self.smax[eid].append(smax)
            self.smin[eid].append(smin)
            #self.MS_tension[eid].append(MSt)
            #self.MS_compression[eid].append(MSc)
        else:
            self.s1[eid] = [s1]
            self.s2[eid] = [s2]
            self.s3[eid] = [s3]
            self.s4[eid] = [s4]
            self.axial[eid] = axial
            self.smax[eid] = [smax]
            self.smin[eid] = [smin]
            #self.MS_tension[eid]     = MSt
            #self.MS_compression[eid] = MSc

        #msg = "eid=%s nodeID=%s fd=%g oxx=%g oyy=%g \ntxy=%g angle=%g major=%g minor=%g vm=%g" %(eid,nodeID,fd,oxx,oyy,txy,angle,majorP,minorP,ovm)
        #print msg
        #if nodeID==0: raise Exception(msg)

    def add_new_eid(self, eType, dt, eid, s1a, s2a, s3a, s4a, axial, smaxa, smina, MSt,
        s1b, s2b, s3b, s4b, smaxb, sminb, MSc):
        #print "Bar Stress add..."
        self.eType[eid] = eType

        self.s1[eid] = [s1a, s1b]
        self.s2[eid] = [s2a, s2b]
        self.s3[eid] = [s3a, s3b]
        self.s4[eid] = [s4a, s4b]
        self.axial[eid] = axial
        self.smax[eid] = [smaxa, smaxb]
        self.smin[eid] = [smina, sminb]
        #self.MS_tension[eid]     = MSt
        #self.MS_compression[eid] = MSc

        #msg = "eid=%s nodeID=%s fd=%g oxx=%g oyy=%g \ntxy=%g angle=%g major=%g minor=%g vm=%g" %(eid,nodeID,fd,oxx,oyy,txy,angle,majorP,minorP,ovm)
        #print msg
        #if nodeID==0: raise Exception(msg)

    def add_new_eid_sort1(self, eType, dt, eid, s1a, s2a, s3a, s4a, axial, smaxa, smina, MSt,
                                             s1b, s2b, s3b, s4b, smaxb, sminb, MSc):
        #msg = "dt=%s eid=%s s1a=%s" % (dt, eid, s1a)
        #print msg
        if dt not in self.s1:
            self.add_new_transient(dt)
        self.eType[eid] = eType
        #print self.s1
        self.s1[dt][eid] = [s1a, s1b]
        self.s2[dt][eid] = [s2a, s2b]
        self.s3[dt][eid] = [s3a, s3b]
        self.s4[dt][eid] = [s4a, s4b]
        self.axial[dt][eid] = axial
        self.smax[dt][eid] = [smaxa, smaxb]
        self.smin[dt][eid] = [smina, sminb]
        #self.MS_tension[dt][eid]     = MSt
        #self.MS_compression[dt][eid] = MSc

        #if nodeID==0: raise Exception(msg)

    def write_f06(self, header, page_stamp, page_num=1, f=None, is_mag_phase=False, is_sort1=True):
        if self.nonlinear_factor is not None:
            return self._write_f06_transient(header, page_stamp, page_num, f, is_sort1=is_sort1)

        msg = header + [
                '                                 S T R E S S E S   I N   B A R   E L E M E N T S          ( C B A R )\n',
                '  ELEMENT        SA1            SA2            SA3            SA4           AXIAL          SA-MAX         SA-MIN     M.S.-T\n',
                '    ID.          SB1            SB2            SB3            SB4           STRESS         SB-MAX         SB-MIN     M.S.-C\n',
              ]

        f.write(''.join(msg))
        for eid, S1s in sorted(iteritems(self.s1)):
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
            (vals2, is_all_zeros) = writeFloats13E(vals)
            [s1a, s2a, s3a, s4a, axial, smaxa, smina,
             s1b, s2b, s3b, s4b, smaxb, sminb] = vals2
            f.write('0%8i   %-13s  %-13s  %-13s  %-13s  %-13s  %-13s  %-13s %s\n'
                    ' %8s   %-13s  %-13s  %-13s  %-13s  %-13s  %-13s  %-13s %s\n'
                    % (eid, s1a, s2a, s3a, s4a, axial, smaxa, smina, MSt,
                       '', s1b, s2b, s3b, s4b, '', smaxb, sminb, MSc))
        f.write(page_stamp % page_num)
        return page_num

    def _write_f06_transient(self, header, page_stamp, page_num=1, f=None, is_mag_phase=False, is_sort1=True):
        words = [
                '                                 S T R E S S E S   I N   B A R   E L E M E N T S          ( C B A R )\n',
                '  ELEMENT        SA1            SA2            SA3            SA4           AXIAL          SA-MAX         SA-MIN     M.S.-T\n',
                '    ID.          SB1            SB2            SB3            SB4           STRESS         SB-MAX         SB-MIN     M.S.-C\n',
              ]
        i = 0
        ntimes = len(self.s1)
        for dt, S1s in sorted(iteritems(self.s1)):
            header = _eigenvalue_header(self, header, i, ntimes, dt)
            f.write(''.join(header + words))
            for eid, S1 in sorted(iteritems(S1s)):
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
                (vals2, is_all_zeros) = writeFloats13E(vals)
                [s1a, s2a, s3a, s4a, axial, smaxa, smina,
                 s1b, s2b, s3b, s4b, smaxb, sminb] = vals2
                f.write('0%8i   %-13s  %-13s  %-13s  %-13s  %-13s  %-13s  %-13s %s\n'
                        ' %8s   %-13s  %-13s  %-13s  %-13s  %-13s  %-13s  %-13s %s\n'
                        % (eid, s1a, s2a, s3a, s4a, axial, smaxa, smina, MSt,
                            '', s1b, s2b, s3b, s4b, '', smaxb, sminb, MSc))
            f.write(page_stamp % page_num)
            msg = ['']
            page_num += 1
            i += 1
        return page_num - 1


class RealBarStrain(StrainObject):
    """
    ::

      # s_code=10
                                       S T R A I N S   I N   B A R   E L E M E N T S          ( C B A R )
      ELEMENT        SA1            SA2            SA3            SA4           AXIAL          SA-MAX         SA-MIN     M.S.-T
        ID.          SB1            SB2            SB3            SB4           STRAIN         SB-MAX         SB-MIN     M.S.-C
    """
    def __init__(self, data_code, is_sort1, isubcase, dt):
        StrainObject.__init__(self, data_code, isubcase)
        self.eType = {}

        self.code = [self.format_code, self.sort_code, self.s_code]
        self.e1 = {}
        self.e2 = {}
        self.e3 = {}
        self.e4 = {}
        self.axial = {}
        self.emax = {}
        self.emin = {}
        self.MS_tension = {}
        self.MS_compression = {}

        if is_sort1:
            if dt is not None:
                self.add_new_eid = self.add_new_eid_sort1
        else:
            assert dt is not None
            self.add_new_eid = self.add_new_eid_sort2

    def get_stats(self):
        nelements = len(self.eType)

        msg = self.get_data_code()
        if self.dt is not None:  # transient
            ntimes = len(self.e1)
            msg.append('  type=%s ntimes=%s nelements=%s\n'
                       % (self.__class__.__name__, ntimes, nelements))
        else:
            msg.append('  type=%s nelements=%s\n' % (self.__class__.__name__,
                                                     nelements))
        msg.append('  eType, e1, e2, e3, e4, axial, emax, emin, '
                   'MS_tension, MS_compression\n')
        return msg

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
        if dt not in self.e1:
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

    def delete_transient(self, dt):
        del self.e1[dt]
        del self.e2[dt]
        del self.e3[dt]
        del self.e4[dt]
        del self.axial[dt]
        del self.emax[dt]
        del self.emin[dt]

    def get_transients(self):
        k = self.e1.keys()
        k.sort()
        return k

    def add_new_transient(self, dt):
        """
        initializes the transient variables
        """
        self.e1[dt] = {}
        self.e2[dt] = {}
        self.e3[dt] = {}
        self.e4[dt] = {}
        self.axial[dt] = {}
        self.emax[dt] = {}
        self.emin[dt] = {}
        #self.MS_tension[dt]     = {}
        #self.MS_compression[dt] = {}

    def add_new_eid(self, eType, dt, eid,
                    e1a, e2a, e3a, e4a, axial, emaxa, emina, MSt,
                    e1b, e2b, e3b, e4b, emaxb, eminb, MSc):
        #print "Bar Stress add..."
        self.eType[eid] = eType
        self.e1[eid] = [e1a, e1b]
        self.e2[eid] = [e2a, e2b]
        self.e3[eid] = [e3a, e3b]
        self.e4[eid] = [e4a, e4b]
        self.axial[eid] = axial
        self.emax[eid] = [emaxa, emaxb]
        self.emin[eid] = [emina, eminb]
        #self.MS_tension[eid]     = MSt
        #self.MS_compression[eid] = MSc

        #msg = "eid=%s nodeID=%s fd=%g oxx=%g oyy=%g \ntxy=%g angle=%g major=%g minor=%g vm=%g" %(eid,nodeID,fd,oxx,oyy,txy,angle,majorP,minorP,ovm)
        #print msg
        #if nodeID==0: raise Exception(msg)

    def add_new_eid_sort1(self, eType, dt, eid,
                          e1a, e2a, e3a, e4a, axial, emaxa, emina, MSt,
                          e1b, e2b, e3b, e4b, emaxb, eminb, MSc):
        #print "Bar Stress add..."

        self.eType[eid] = eType
        if dt not in self.e1:
            self.add_new_transient(dt)

        self.e1[dt][eid] = [e1a, e1b]
        self.e2[dt][eid] = [e2a, e2b]
        self.e3[dt][eid] = [e3a, e3b]
        self.e4[dt][eid] = [e4a, e4b]
        self.axial[dt][eid] = axial
        self.emax[dt][eid] = [emaxa, emaxb]
        self.emin[dt][eid] = [emina, eminb]
        #self.MS_tension[dt][eid]     = MSt
        #self.MS_compression[dt][eid] = MSc

        #msg = "eid=%s nodeID=%s fd=%g oxx=%g oyy=%g \ntxy=%g angle=%g major=%g minor=%g vm=%g" %(eid,nodeID,fd,oxx,oyy,txy,angle,majorP,minorP,ovm)
        #print msg
        #if nodeID==0: raise Exception(msg)

    def write_f06(self, header, page_stamp, page_num=1, f=None, is_mag_phase=False, is_sort1=True):
        if self.nonlinear_factor is not None:
            return self._write_f06_transient(header, page_stamp, page_num, f, is_sort1=is_sort1)

        msg = header + [
                '                                  S T R A I N S    I N   B A R   E L E M E N T S          ( C B A R )\n',
                '  ELEMENT        SA1            SA2            SA3            SA4           AXIAL          SA-MAX         SA-MIN     M.S.-T\n',
                '    ID.          SB1            SB2            SB3            SB4           STRAIN         SB-MAX         SB-MIN     M.S.-C\n',
              ]
        for eid, E1s in sorted(iteritems(self.e1)):
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
            (vals2, is_all_zeros) = writeFloats13E(vals)
            [e10, e20, e30, e40, axial, emax0, emin0,
             e11, e21, e31, e41, emax1, emin1] = vals2

            msg.append('0%8i   %-13s  %-13s  %-13s  %-13s  %-13s  %-13s  %-13s %s\n' % (eid, e10, e20, e30, e40, axial, emax0, emin0, MSt))
            msg.append(' %8s   %-13s  %-13s  %-13s  %-13s  %-13s  %-13s  %-13s %s\n' % ('', e11, e21, e31, e41, '', emax1, emin1, MSc))

        msg.append(page_stamp % page_num)
        f.write(''.join(msg))
        return page_num

    def _write_f06_transient(self, header, page_stamp, page_num=1, f=None, is_mag_phase=False, is_sort1=True):
        words = [
                '                                  S T R A I N S    I N   B A R   E L E M E N T S           ( C B A R )\n',
                '  ELEMENT        SA1            SA2            SA3            SA4           AXIAL          SA-MAX         SA-MIN     M.S.-T\n',
                '    ID.          SB1            SB2            SB3            SB4           STRAIN         SB-MAX         SB-MIN     M.S.-C\n',
              ]
        msg = []
        i = 0
        ntimes = len(self.e1)
        for dt, E1s in sorted(iteritems(self.e1)):
            header = _eigenvalue_header(self, header, i, ntimes, dt)
            msg += header + words
            for eid, e1s in sorted(iteritems(E1s)):
                #eType = self.eType[eid]
                axial = self.axial[dt][eid]
                #MSt = self.MSt[eid]
                #MSc = self.MSc[eid]
                MSt = ''
                MSc = ''

                e1 = self.e1[dt][eid]
                e2 = self.e2[dt][eid]
                e3 = self.e3[dt][eid]
                e4 = self.e4[dt][eid]
                emax = self.emax[dt][eid]
                emin = self.emin[dt][eid]
                vals = [e1[0], e2[0], e3[0], e4[0], axial, emax[0], emin[0],
                        e1[1], e2[1], e3[1], e4[1], emax[1], emin[1]]
                (vals2, is_all_zeros) = writeFloats13E(vals)
                [e10, e20, e30, e40, axial, emax0, emin0,
                 e11, e21, e31, e41, emax1, emin1] = vals2

                msg.append('0%8i   %-13s  %-13s  %-13s  %-13s  %-13s  %-13s  %-13s %s\n' % (eid, e10, e20, e30, e40, axial, emax0, emin0, MSt))
                msg.append(' %8s   %-13s  %-13s  %-13s  %-13s  %-13s  %-13s  %-13s %s\n' % ('', e11, e21, e31, e41, '', emax1, emin1, MSc))

            msg.append(page_stamp % page_num)
            f.write(''.join(msg))
            msg = ['']
            page_num += 1
            i += 1
        return page_num - 1
