from __future__ import (nested_scopes, generators, division, absolute_import,
                        print_function, unicode_literals)
from six import iteritems
from pyNastran.op2.tables.oes_stressStrain.real.oes_objects import StressObject, StrainObject, OES_Object
from pyNastran.f06.f06_formatting import writeFloats13E, writeImagFloats13E
from numpy import concatenate, zeros


class ComplexBarArray(OES_Object):
    def __init__(self, data_code, is_sort1, isubcase, dt):
        OES_Object.__init__(self, data_code, isubcase, apply_data_code=False)
        self.eType = {}
        self.result_flag = 0
        #self.code = [self.format_code, self.sort_code, self.s_code]

        #self.ntimes = 0  # or frequency/mode
        #self.ntotal = 0
        self.itime = 0
        self.nelements = 0  # result specific
        #self.cid = {}  # gridGauss

        if is_sort1:
            #sort1
            pass
        else:
            raise NotImplementedError('SORT2')

    def get_headers(self):
        raise RuntimeError()

    def _reset_indices(self):
        self.itotal = 0
        self.ielement = 0

    #def get_nnodes(self):
        #if self.element_type == 39: # CTETRA
            #nnodes = 5
        #elif self.element_type == 68: # CPENTA
            #nnodes = 7
        #elif self.element_type == 67: # CHEXA
            #nnodes = 9
        #else:
            #raise NotImplementedError(self.element_name)
        #return nnodes

    def build(self):
        #print('ntimes=%s nelements=%s ntotal=%s subtitle=%s' % (self.ntimes, self.nelements, self.ntotal, self.subtitle))
        if self.is_built:
            return
        nnodes = 1

        #self.names = []
        #self.nelements //= nnodes
        self.nelements /= self.ntimes
        #self.ntotal //= self.ntimes
        self.itime = 0
        self.ielement = 0
        self.itotal = 0
        self.is_built = True
        #print('ntotal=%s ntimes=%s nelements=%s' % (self.ntotal, self.ntimes, self.nelements))
        #print("ntimes=%s nelements=%s ntotal=%s" % (self.ntimes, self.nelements, self.ntotal))
        self._times = zeros(self.ntimes, 'float32')
        #self.element = array(self.nelements, dtype='|S8')

        #self.ntotal = self.nelements * nnodes
        self.element = zeros(self.ntotal, 'int32')

        # the number is messed up because of the offset for the element's properties
        if not self.nelements * nnodes == self.ntotal:
            msg = 'ntimes=%s nelements=%s nnodes=%s ne*nn=%s ntotal=%s' % (self.ntimes,
                                                                           self.nelements, nnodes,
                                                                           self.nelements * nnodes,
                                                                           self.ntotal)
            raise RuntimeError(msg)

        #[s1a, s2a, s3a, s4a, axial, s2a, s2b, s2c, s2d]
        self.data = zeros((self.ntimes, self.ntotal, 9), 'complex64')

    def add_new_eid_sort1(self, eType, dt, eid, e1a, e2a, e3a, e4a, axial,
                          e1b, e2b, e3b, e4b,):
        #self.e1[dt][eid] = [e1a, e1b]
        #self.e2[dt][eid] = [e2a, e2b]
        #self.e3[dt][eid] = [e3a, e3b]
        #self.e4[dt][eid] = [e4a, e4b]
        #self.axial[dt][eid] = axial

        #[sa1, sa2, sa3, sa4, axial, sb1, sb2, sb3, sb4]
        self.data[self.itime, self.itotal, :] = [e1a, e2a, e3a, e4a, axial,
                                                 e1b, e2b, e3b, e4b,]
        self.element[self.itotal] = eid
        self.itotal += 1

    def get_stats(self):
        if not self.is_built:
            return [
                '<%s>\n' % self.__class__.__name__,
                '  ntimes: %i\n' % self.ntimes,
                '  ntotal: %i\n' % self.ntotal,
            ]

        nelements = self.nelements
        ntimes = self.ntimes
        #ntotal = self.ntotal
        nnodes = self.element_node.shape[0]
        msg = []

        if self.nonlinear_factor is not None:  # transient
            msg.append('  type=%s ntimes=%i nelements=%i nnodes=%i\n'
                       % (self.__class__.__name__, ntimes, nelements, nnodes))
        else:
            msg.append('  type=%s nelements=%i nnodes=%i\n' % (self.__class__.__name__, nelements, nnodes))
        msg.append('  eType, cid\n')
        msg.append('  data: [ntimes, nnodes, 6] where 6=[%s]\n' % str(', '.join(self.get_headers())))
        msg.append('  data.shape = %s\n' % str(self.data.shape).replace('L', ''))

        msg.append(', '.join(set(self.eType.values())))
        msg += self.get_data_code()
        return msg

    def write_f06(self, header, page_stamp, page_num=1, f=None, is_mag_phase=False, is_sort1=True):
        #msg_temp, nnodes = get_f06_header(self, is_mag_phase, is_sort1)

        # write the f06
        ntimes = self.data.shape[0]

        # TODO: these are not
        if is_mag_phase:
            mag_phase = '                                                          (MAGNITUDE/PHASE)'
        else:
            mag_phase = '                                                          (REAL/IMAGINARY)'

        # force
        #msg = header + [
            #'                             C O M P L E X   F O R C E S   I N   B A R   E L E M E N T S   ( C B A R )',
            #mag_phase,
            #' ',
            #'                     BEND-MOMENT-END-A            BEND-MOMENT-END-B                  SHEAR',
            #'   FREQUENCY       PLANE 1       PLANE 2        PLANE 1       PLANE 2        PLANE 1       PLANE 2        FORCE          TORQUE',
        #]

        if is_sort1:
            line1 = '            ELEMENT                    LOCATION       LOCATION       LOCATION       LOCATION             AVERAGE\n',
            line2 = '              ID.                          1              2              3              4             AXIAL STRESS\n',
            #line1 = '    ELEMENT          BEND-MOMENT-END-A            BEND-MOMENT-END-B                  SHEAR'
            #line2 = '      ID.          PLANE 1       PLANE 2        PLANE 1       PLANE 2        PLANE 1       PLANE 2        FORCE          TORQUE'
        else:
            #line1 = '                     BEND-MOMENT-END-A            BEND-MOMENT-END-B                  SHEAR'
            line1 = '                                       LOCATION       LOCATION       LOCATION       LOCATION             AVERAGE\n',
            if name == 'freq':
                line2 = '           FREQUENCY                       1              2              3              4             AXIAL STRESS\n',
                #line2 = '   FREQUENCY       PLANE 1       PLANE 2        PLANE 1       PLANE 2        PLANE 1       PLANE 2        FORCE          TORQUE'
            else:
                raise RuntimeError(name)

        if self.is_stress():
            stress_strain = '                             C O M P L E X   S T R E S S E S   I N   B A R   E L E M E N T S   ( C B A R )'
        else:
            stress_strain = '                             C O M P L E X   S T R A I N S   I N   B A R   E L E M E N T S   ( C B A R )'

        msg_temp = header + [
            stress_strain,
            mag_phase,
            ' ',
            line1,
            line2,
        ]

        for itime in range(ntimes):
            dt = self._times[itime]

            dt_line = ' %14s = %12.5E\n' % (self.data_code['name'], dt)
            header[1] = dt_line
            msg = header + msg_temp
            f.write('\n'.join(msg))

            # TODO: can I get this without a reshape?
            sa1 = self.data[itime, :, 0]
            sa2 = self.data[itime, :, 1]
            sa3 = self.data[itime, :, 2]
            sa4 = self.data[itime, :, 3]
            axial = self.data[itime, :, 4]
            sb1 = self.data[itime, :, 5]
            sb2 = self.data[itime, :, 6]
            sb3 = self.data[itime, :, 7]
            sb4 = self.data[itime, :, 8]
            #[sa1, sa2, sa3, sa4, axial, sb1, sb2, sb3, sb4]

            eids = self.element_node[:, 0]
            for eid, s1ai, s2ai, s3ai, s4ai, axiali, s2ai, s2bi, s2ci, s2di in zip(eids, sa1, sa2, sa3, sa4, axial, sb1, sb2, sb3, sb4):
                vals = (s1ai, s2ai, s3ai, s4ai, axiali,
                        s2ai, s2bi, s2ci, s2di)
                (vals2, is_all_zeros) = writeImagFloats13E(vals, is_mag_phase)
                (s1ar, s2ar, s3ar, s4ar, axialr,
                 s1br, s2br, s3br, s4br,
                 s1ai, s2ai, s3ai, s4ai, axiali,
                 s1bi, s2bi, s3bi, s4bi) = vals2

                msg.append('0%8i   %-13s  %-13s  %-13s  %-13s  %s\n' % (eid, s1ar, s2ar, s3ar, s4ar, axialr))
                msg.append(' %8s   %-13s  %-13s  %-13s  %-13s  %s\n' % ('', s1ai, s2ai, s3ai, s4ai, axiali))

                msg.append(' %8s   %-13s  %-13s  %-13s  %s\n' % ('', s1br, s2br, s3br, s4br))
                msg.append(' %8s   %-13s  %-13s  %-13s  %s\n' % ('', s1bi, s2bi, s3bi, s4bi))
            f.write(page_stamp % page_num)
            page_num += 1
        return page_num - 1


class ComplexBarStressArray(ComplexBarArray, StressObject):
    def __init__(self, data_code, is_sort1, isubcase, dt):
        ComplexBarArray.__init__(self, data_code, is_sort1, isubcase, dt)
        StressObject.__init__(self, data_code, isubcase)

    def get_headers(self):
        headers = ['s1a', 's1b', 's1c', 's1d', 's1e', 'axial',
                   's2a', 's2b', 's2c', 's2d', 's2e', ]
        return headers

class ComplexBarStrainArray(ComplexBarArray, StrainObject):
    def __init__(self, data_code, is_sort1, isubcase, dt):
        ComplexBarArray.__init__(self, data_code, is_sort1, isubcase, dt)
        StrainObject.__init__(self, data_code, isubcase)

    def get_headers(self):
        headers = ['e1a', 'e1b', 'e1c', 'e1d', 'e1e', 'axial',
                   'e2a', '2b', 'e2c', 'e2d', 'e2e', ]
        return headers


class ComplexBarStress(StressObject):
    """
    # s_code=0
                           C O M P L E X   S T R E S S E S   I N   B A R   E L E M E N T S   ( C B A R )
                                                         (MAGNITUDE/PHASE)

            ELEMENT                    LOCATION       LOCATION       LOCATION       LOCATION             AVERAGE
              ID.                          1              2              3              4             AXIAL STRESS

                  1     ENDA          9.331276E+04   9.331276E+04   9.331276E+04   9.331276E+04        0.0
                                      180.0000         0.0            0.0          180.0000              0.0
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

        #if self.element_type==100:
            #self.getLength = self.getLength100_format1_sort0
            #self.add_new_eid = self.addNewEid100

        self.dt = dt

    def get_stats(self):
        nelements = len(self.eType)

        msg = self.get_data_code()
        if self.dt is not None:  # transient
            ntimes = len(self.axial)
            msg.append('  type=%s ntimes=%s nelements=%s\n'
                       % (self.__class__.__name__, ntimes, nelements))
        else:
            msg.append('  type=%s nelements=%s\n' % (self.__class__.__name__,
                                                     nelements))
        msg.append('  eType, s1, s2, s3, s4, axial\n')
        return msg

    def add_f06_data(self, data, transient):
        if transient is None:
            raise RuntimeError(transient)
            #for line in data:
                #(eType, eid, s1A, s2A, s3A, s4A, axialA,
                 #s1B, s2B, s3B, s4B,) = line
                #self.eType[eid] = 'CBAR'
                #self.s1[eid] = [s1A, s1B]
                #self.s2[eid] = [s2A, s2B]
                #self.s3[eid] = [s3A, s3B]
                #self.s4[eid] = [s4A, s4B]
                #self.axial[eid] = axialA
            #return

        (dtName, dt) = transient
        self.data_code['name'] = dtName
        #print "dt = ",dt
        #print "dtName = ",dtName
        if dt not in self.s1:
            self.update_dt(self.data_code, dt)

        for line in data:
            (eType, eid, s1A, s2A, s3A, s4A, axialA,
             s1B, s2B, s3B, s4B) = line
            self.eType[eid] = 'CBAR'
            self.s1[dt][eid] = [s1A, s1B]
            self.s2[dt][eid] = [s2A, s2B]
            self.s3[dt][eid] = [s3A, s3B]
            self.s4[dt][eid] = [s4A, s4B]
            self.axial[dt][eid] = axialA

    def delete_transient(self, dt):
        del self.s1[dt]
        del self.s2[dt]
        del self.s3[dt]
        del self.s4[dt]
        del self.axial[dt]

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

    def add_new_eid_sort1(self, eType, dt, eid, s1a, s2a, s3a, s4a, axial,
                       s1b, s2b, s3b, s4b):
        msg = "dt=%s eid=%s s1a=%s" % (dt, eid, s1a)
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

        #if nodeID==0: raise Exception(msg)

    def write_f06(self, header, page_stamp, page_num=1, f=None, is_mag_phase=False, is_sort1=True):
        if self.nonlinear_factor is not None:
            return self._write_f06_transient(header, page_stamp, page_num, f, is_mag_phase, is_sort1=is_sort1)

        msg = header + [
            '                                 S T R E S S E S   I N   B A R   E L E M E N T S          ( C B A R )\n',
            '  ELEMENT        SA1            SA2            SA3            SA4           AXIAL          SA-MAX         SA-MIN     M.S.-T\n',
            '    ID.          SB1            SB2            SB3            SB4           STRESS         SB-MAX         SB-MIN     M.S.-C\n',
        ]

        for eid, S1s in sorted(iteritems(self.s1)):
            eType = self.eType[eid]
            axial = self.axial[eid]
            s1 = self.s1[eid]
            s2 = self.s2[eid]
            s3 = self.s3[eid]
            s4 = self.s4[eid]

            vals = [s1[0], s2[0], s3[0], s4[0], axial,
                    s1[1], s2[1], s3[1], s4[1], ]
            (vals2, is_all_zeros) = writeImagFloats13E(vals, is_mag_phase)
            [s1ar, s2ar, s3ar, s4ar, axialr,
             s1br, s2br, s3br, s4br,
             s1ai, s2ai, s3ai, s4ai, axiali,
             s1bi, s2bi, s3bi, s4bi, ] = vals2
            msg.append('0%8i   %-13s  %-13s  %-13s  %-13s  %s\n' % (eid, s1ar, s2ar, s3ar, s4ar, axialr))
            msg.append(' %8s   %-13s  %-13s  %-13s  %-13s  %s\n' % ('', s1ai, s2ai, s3ai, s4ai, axiali))

            msg.append(' %8s   %-13s  %-13s  %-13s  %s\n' % ('', s1br, s2br, s3br, s4br))
            msg.append(' %8s   %-13s  %-13s  %-13s  %s\n' % ('', s1bi, s2bi, s3bi, s4bi))

        msg.append(page_stamp % page_num)
        f.write(''.join(msg))
        return page_num

    def _write_f06_transient(self, header, page_stamp, page_num=1, f=None, is_mag_phase=False, is_sort1=True):
        words = [
            '                                 S T R E S S E S   I N   B A R   E L E M E N T S          ( C B A R )\n',
            '  ELEMENT        SA1            SA2            SA3            SA4           AXIAL          SA-MAX         SA-MIN     M.S.-T\n',
            '    ID.          SB1            SB2            SB3            SB4           STRESS         SB-MAX         SB-MIN     M.S.-C\n',
        ]
        msg = []
        for dt, S1s in sorted(iteritems(self.s1)):
            header[1] = ' %s = %10.4E\n' % (self.data_code['name'], dt)
            msg += header + words
            for eid, S1 in sorted(iteritems(S1s)):
                eType = self.eType[eid]
                axial = self.axial[dt][eid]
                s1 = self.s1[dt][eid]
                s2 = self.s2[dt][eid]
                s3 = self.s3[dt][eid]
                s4 = self.s4[dt][eid]
                vals = [s1[0], s2[0], s3[0], s4[0], axial,
                        s1[1], s2[1], s3[1], s4[1], ]
                (vals2, is_all_zeros) = writeImagFloats13E(vals, is_mag_phase)
                [s1ar, s2ar, s3ar, s4ar, axialr,
                 s1br, s2br, s3br, s4br,
                 s1ai, s2ai, s3ai, s4ai, axiali,
                 s1bi, s2bi, s3bi, s4bi, ] = vals2
                msg.append('0%8i   %-13s  %-13s  %-13s  %-13s  %s\n' % (eid, s1ar, s2ar, s3ar, s4ar, axialr))
                msg.append(' %8s   %-13s  %-13s  %-13s  %-13s  %s\n' % ('', s1ai, s2ai, s3ai, s4ai, axiali))

                msg.append(' %8s   %-13s  %-13s  %-13s  %s\n' % ('', s1br, s2br, s3br, s4br))
                msg.append(' %8s   %-13s  %-13s  %-13s  %s\n' % ('', s1bi, s2bi, s3bi, s4bi))

            msg.append(page_stamp % page_num)
            f.write(''.join(msg))
            msg = ['']
            page_num += 1
        return page_num - 1


class ComplexBarStrain(StrainObject):
    """
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

    def get_stats(self):
        nelements = len(self.eType)

        msg = self.get_data_code()
        if self.dt is not None:  # transient
            ntimes = len(self.axial)
            msg.append('  type=%s ntimes=%s nelements=%s\n'
                       % (self.__class__.__name__, ntimes, nelements))
        else:
            msg.append('  type=%s nelements=%s\n' % (self.__class__.__name__,
                                                     nelements))
        msg.append('  eType, e1, e2, e3, e4, axial\n')
        return msg

    def add_f06_data(self, data, transient):
        if transient is None:
            raise RuntimeError(transient)

        (dtName, dt) = transient
        self.data_code['name'] = dtName
        if dt not in self.e1:
            self.update_dt(self.data_code, dt)

        for line in data:
            (eType, eid, e1A, e2A, e3A, e4A, axialA,
             e1B, e2B, e3B, e4B,) = line
            self.eType[eid] = 'CBAR'
            self.e1[dt][eid] = [e1A, e1B]
            self.e2[dt][eid] = [e2A, e2B]
            self.e3[dt][eid] = [e3A, e3B]
            self.e4[dt][eid] = [e4A, e4B]
            self.axial[dt][eid] = axialA

    def delete_transient(self, dt):
        del self.e1[dt]
        del self.e2[dt]
        del self.e3[dt]
        del self.e4[dt]
        del self.axial[dt]

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

    #def add_new_eid(self, eType, dt, eid, e1a, e2a, e3a, e4a, axial,
                    #e1b, e2b, e3b, e4b,):
        ##print "Bar Stress add..."
        #self.eType[eid] = eType
        #self.e1[eid] = [e1a, e1b]
        #self.e2[eid] = [e2a, e2b]
        #self.e3[eid] = [e3a, e3b]
        #self.e4[eid] = [e4a, e4b]
        #self.axial[eid] = axial

        #msg = "eid=%s nodeID=%s fd=%g oxx=%g oyy=%g \ntxy=%g angle=%g major=%g minor=%g vm=%g" %(eid,nodeID,fd,oxx,oyy,txy,angle,majorP,minorP,ovm)
        #print msg
        #if nodeID==0: raise Exception(msg)

    def add_new_eid_sort1(self, eType, dt, eid, e1a, e2a, e3a, e4a, axial,
                          e1b, e2b, e3b, e4b,):
        #print "Bar Stress add..."

        self.eType[eid] = eType
        if dt not in self.e1:
            self.add_new_transient(dt)

        self.e1[dt][eid] = [e1a, e1b]
        self.e2[dt][eid] = [e2a, e2b]
        self.e3[dt][eid] = [e3a, e3b]
        self.e4[dt][eid] = [e4a, e4b]
        self.axial[dt][eid] = axial

        #print msg
        #if nodeID==0: raise Exception(msg)

    def write_f06(self, header, page_stamp, page_num=1, f=None, is_mag_phase=False, is_sort1=True):
        f.write('%s write_f06 not implemented...\n' % self.__class__.__name__)
        return page_num

        if self.nonlinear_factor is not None:
            return self._write_f06_transient(header, page_stamp, page_num, f, is_mag_phase, is_sort1=is_sort1)

        msg = header + [
            '                                  S T R A I N S    I N   B A R   E L E M E N T S          ( C B A R )\n',
            '  ELEMENT        SA1            SA2            SA3            SA4           AXIAL          SA-MAX         SA-MIN     M.S.-T\n',
            '    ID.          SB1            SB2            SB3            SB4           STRAIN         SB-MAX         SB-MIN     M.S.-C\n',
        ]
        for eid, E1s in sorted(iteritems(self.e1)):
            eType = self.eType[eid]
            axial = self.axial[eid]

            e1 = self.e1[eid]
            e2 = self.e2[eid]
            e3 = self.e3[eid]
            e4 = self.e4[eid]
            vals = [e1[0], e2[0], e3[0], e4[0], axial,
                    e1[1], e2[1], e3[1], e4[1]]
            #(vals2, is_all_zeros) = writeFloats13E(vals)
            #[e10, e20, e30, e40, axial,
             #e11, e21, e31, e41] = vals2

            #vals = (s1ai, s2ai, s3ai, s4ai, axiali,
                    #s2ai, s2bi, s2ci, s2di)
            (vals2, is_all_zeros) = writeImagFloats13E(vals, is_mag_phase)
            (s1ar, s2ar, s3ar, s4ar, axialr,
             s1br, s2br, s3br, s4br,
             s1ai, s2ai, s3ai, s4ai, axiali,
             s1bi, s2bi, s3bi, s4bi) = vals2

            msg.append('0%8i   %-13s  %-13s  %-13s  %-13s  %s\n' % (eid, s1ar, s2ar, s3ar, s4ar, axialr))
            msg.append(' %8s   %-13s  %-13s  %-13s  %-13s  %s\n' % ('', s1ai, s2ai, s3ai, s4ai, axiali))

            msg.append(' %8s   %-13s  %-13s  %-13s  %s\n' % ('', s1br, s2br, s3br, s4br))
            msg.append(' %8s   %-13s  %-13s  %-13s  %s\n' % ('', s1bi, s2bi, s3bi, s4bi))

            #msg.append('0%8i   %-13s  %-13s  %-13s  %-13s  %-13s  %-13s  %-13s %s\n' % (eid, e10, e20, e30, e40, axial))
            #msg.append(' %8s   %-13s  %-13s  %-13s  %-13s  %-13s  %-13s  %-13s %s\n' % ('', e11, e21, e31, e41))
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
        for dt, E1s in sorted(iteritems(self.e1)):
            header[1] = ' %s = %10.4E\n' % (self.data_code['name'], dt)
            msg += header + words
            for eid, e1s in sorted(iteritems(E1s)):
                eType = self.eType[eid]
                axial = self.axial[eid]

                e1 = self.e1[eid]
                e2 = self.e2[eid]
                e3 = self.e3[eid]
                e4 = self.e4[eid]
                vals = [e1[0], e2[0], e3[0], e4[0], axial,
                        e1[1], e2[1], e3[1], e4[1]]
                (vals2, is_all_zeros) = writeImagFloats13E(vals, is_mag_phase)
                (s1ar, s2ar, s3ar, s4ar, axialr,
                 s1br, s2br, s3br, s4br,
                 s1ai, s2ai, s3ai, s4ai, axiali,
                 s1bi, s2bi, s3bi, s4bi) = vals2

                msg.append('0%8i   %-13s  %-13s  %-13s  %-13s  %s\n' % (eid, s1ar, s2ar, s3ar, s4ar, axialr))
                msg.append(' %8s   %-13s  %-13s  %-13s  %-13s  %s\n' % ('', s1ai, s2ai, s3ai, s4ai, axiali))

                msg.append(' %8s   %-13s  %-13s  %-13s  %s\n' % ('', s1br, s2br, s3br, s4br))
                msg.append(' %8s   %-13s  %-13s  %-13s  %s\n' % ('', s1bi, s2bi, s3bi, s4bi))

            msg.append(page_stamp % page_num)
            f.write(''.join(msg))
            msg = ['']
            page_num += 1
        return page_num - 1
