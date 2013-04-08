from __future__ import (nested_scopes, generators, division, absolute_import,
                        print_function, unicode_literals)

from ..real.oes_objects import StressObject


class ComplexBush1DStressObject(StressObject):
    """
    # s_code=0
                           C O M P L E X   S T R E S S E S   I N   B A R   E L E M E N T S   ( C B A R )
                                                         (MAGNITUDE/PHASE)

            ELEMENT                    LOCATION       LOCATION       LOCATION       LOCATION             AVERAGE
              ID.                          1              2              3              4             AXIAL STRESS

                  1     ENDA          9.331276E+04   9.331276E+04   9.331276E+04   9.331276E+04        0.0
                                      180.0000         0.0            0.0          180.0000              0.0
    """
    def __init__(self, data_code, is_sort1, isubcase, dt=None):
        StressObject.__init__(self, data_code, isubcase)
        self.eType = {}

        self.code = [self.format_code, self.sort_code, self.s_code]

        self.element_force = {}
        self.axial_displacement = {}
        self.axial_stress = {}
        self.axial_strain = {}

        self.dt = dt
        if is_sort1:
            if dt is not None:
                #self.add = self.add_sort1
                self.add_new_eid = self.add_new_eid_sort1
        else:
            assert dt is not None
            #self.add = self.addSort2
            self.add_new_eid = self.add_new_eid_sort2

    def get_stats(self):
        nelements = len(self.eType)

        msg = self.get_data_code()
        if self.dt is not None:  # transient
            ntimes = len(self.element_force)
            msg.append('  type=%s ntimes=%s nelements=%s\n'
                       % (self.__class__.__name__, ntimes, nelements))
        else:
            msg.append('  imaginary type=%s nelements=%s\n' % (self.__class__.__name__,
                                                     nelements))
        msg.append('  eType, element_force, axial_displacement, axial_velocity, axial_stress, axial_strain\n')
        return msg

    def add_f06_data(self, data, transient):
        raise NotImplementedError('CBUSH1D')
        if transient is None:
            for line in data:
                (eType, eid, fe, ue, ao, ae) = line
                self.eType[eid] = 'CBUSH1D'
                self.element_force[eid] = fe
                self.axial_displacement[eid] = ue
                self.axial_stress[eid] = ao
                self.axial_strain[eid] = ae
            return

        (dtName, dt) = transient
        self.data_code['name'] = dtName

        if dt not in self.element_force:
            self.update_dt(self.data_code, dt)

        for line in data:
            (eType, eid, fe, ue, ao, ae) = line
            self.eType[eid] = 'CBUSH1D'
            self.element_force[dt][eid] = fe
            self.axial_displacement[dt][eid] = ue
            self.axial_stress[dt][eid] = ao
            self.axial_strain[dt][eid] = ae

    def delete_transient(self, dt):
        del self.element_force[dt]
        del self.axial_displacement[dt]
        del self.axial_stress[dt]
        del self.axial_strain[dt]

    def get_transients(self):
        k = self.element_force.keys()
        k.sort()
        return k

    def add_new_transient(self, dt):
        """
        initializes the transient variables
        """
        #print "****add new transient****"
        self.dt = dt
        self.element_force[dt] = {}
        self.axial_displacement[dt] = {}
        self.axial_stress[dt] = {}
        self.axial_strain[dt] = {}

    def add_new_eid(self, eType, dt, eid, fe, ue, ao, ae):
        #print "Bush1d Stress add..."
        self.eType[eid] = eType

        self.element_force[eid] = fe
        self.axial_displacement[eid] = ue
        self.axial_stress[eid] = ao
        self.axial_strain[eid] = ae

    def add_new_eid_sort1(self, eType, dt, eid, fe, ue, ao, ae):
        #msg = "dt=%s eid=%s fe=%s" % (dt, eid, fe)
        #print msg
        if dt not in self.element_force:
            self.add_new_transient(dt)
        self.eType[eid] = eType
        self.element_force[dt][eid] = fe
        self.axial_displacement[dt][eid] = ue
        self.axial_stress[dt][eid] = ao
        self.axial_strain[dt][eid] = ae

    def write_f06(self, header, pageStamp, pageNum=1, f=None, isMagPhase=False):
        raise NotImplementedError('CBUSH1D')
        if self.nonlinear_factor is not None:
            return self._write_f06_transient(header, pageStamp, pageNum, f, isMagPhase)

        msg = header + [
            '                                 S T R E S S E S   I N   B A R   E L E M E N T S          ( C B A R )\n',
            '  ELEMENT        SA1            SA2            SA3            SA4           AXIAL          SA-MAX         SA-MIN     M.S.-T\n',
            '    ID.          SB1            SB2            SB3            SB4           STRESS         SB-MAX         SB-MIN     M.S.-C\n',
        ]

        for eid, S1s in sorted(self.s1.iteritems()):
            eType = self.eType[eid]
            axial = self.axial[eid]
            s1 = self.s1[eid]
            s2 = self.s2[eid]
            s3 = self.s3[eid]
            s4 = self.s4[eid]

            vals = [s1[0], s2[0], s3[0], s4[0], axial,
                    s1[1], s2[1], s3[1], s4[1], ]
            (vals2, isAllZeros) = writeImagFloats13E(vals, isMagPhase)
            [s1ar, s2ar, s3ar, s4ar, axialr,
             s1br, s2br, s3br, s4br,
             s1ai, s2ai, s3ai, s4ai, axiali,
             s1bi, s2bi, s3bi, s4bi, ] = vals2
            msg.append('0%8i   %13s  %13s  %13s  %13s  %-s\n' %
                       (eid, s1ar, s2ar, s3ar, s4ar, axialr.rstrip()))
            msg.append(' %8s   %13s  %13s  %13s  %13s  %-s\n' %
                       ('', s1ai, s2ai, s3ai, s4ai, axiali.rstrip()))

            msg.append(' %8s   %13s  %13s  %13s  %-s\n' % (
                '', s1br, s2br, s3br, s4br.rstrip()))
            msg.append(' %8s   %13s  %13s  %13s  %-s\n' % (
                '', s1bi, s2bi, s3bi, s4bi.rstrip()))

        msg.append(pageStamp + str(pageNum) + '\n')
        return (''.join(msg), pageNum)

    def _write_f06_transient(self, header, pageStamp, pageNum=1, f=None, isMagPhase=False):
        raise NotImplementedError('CBUSH1D')
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
                eType = self.eType[eid]
                axial = self.axial[dt][eid]
                s1 = self.s1[dt][eid]
                s2 = self.s2[dt][eid]
                s3 = self.s3[dt][eid]
                s4 = self.s4[dt][eid]
                vals = [s1[0], s2[0], s3[0], s4[0], axial,
                        s1[1], s2[1], s3[1], s4[1], ]
                (vals2, isAllZeros) = writeImagFloats13E(vals, isMagPhase)
                [s1ar, s2ar, s3ar, s4ar, axialr,
                 s1br, s2br, s3br, s4br,
                 s1ai, s2ai, s3ai, s4ai, axiali,
                 s1bi, s2bi, s3bi, s4bi, ] = vals2
                msg.append('0%8i   %13s  %13s  %13s  %13s  %-s\n' % (eid,
                                                                     s1ar, s2ar, s3ar, s4ar, axialr.rstrip()))
                msg.append(' %8s   %13s  %13s  %13s  %13s  %-s\n' % ('',
                                                                     s1ai, s2ai, s3ai, s4ai, axiali.rstrip()))

                msg.append(' %8s   %13s  %13s  %13s  %-s\n' %
                           ('', s1br, s2br, s3br, s4br.rstrip()))
                msg.append(' %8s   %13s  %13s  %13s  %-s\n' %
                           ('', s1bi, s2bi, s3bi, s4bi.rstrip()))

            msg.append(pageStamp + str(pageNum) + '\n')
            pageNum += 1
        return (''.join(msg), pageNum - 1)

    def __repr__(self):
        raise NotImplementedError('CBUSH1D')
        if self.nonlinear_factor is not None:
            return self.__reprTransient__()

        msg = '---BAR STRESS---\n'
        msg += '%-6s %6s ' % ('EID', 'eType')
        headers = ['s1', 's2', 's3', 's4', 'Axial']
        for header in headers:
            msg += '%8s ' % (header)
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
                    msg += '%8s ' % ('0')
                else:
                    msg += '%8i ' % (val)
            msg += '\n'

            msg += '%s ' % (' ' * 13)
            vals = [s1[1], s2[1], s3[1], s4[1]]
            for val in vals:
                if isinstance(val, str):
                    msg += '%8s ' % (val)
                elif abs(val) < 1e-6:
                    msg += '%8s ' % ('0')
                else:
                    msg += '%8i ' % (val)
            msg += '\n'

            #msg += "eid=%-4s eType=%s s1=%-4i s2=%-4i s3=%-4i s4=%-4i axial=-%5i\n" %(eid,eType,s1[0],s2[0],s3[0],s4[0],axial)
            #msg += "%s                s1=%-4i s2=%-4i s3=%-4i s4=%-4i %s\n"         %(' '*4,    s1[1],s2[1],s3[1],s4[1])
        return msg

    def __reprTransient__(self):
        raise NotImplementedError('CBUSH1D')
        msg = '---BAR STRESS---\n'
        msg += '%-6s %6s ' % ('EID', 'eType')
        headers = ['s1', 's2', 's3', 's4', 'Axial', 'sMax', 'sMin']
        for header in headers:
            msg += '%8s ' % (header)
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
                        msg += '%8s ' % (val)
                    elif abs(val) < 1e-6:
                        msg += '%8s %8s' % (val.real, val.imag)
                msg += '\n'

                #msg += "eid=%-4s eType=%s s1=%-4i s2=%-4i s3=%-4i s4=%-4i axial=-%5i\n" %(eid,eType,s1[0],s2[0],s3[0],s4[0],axial)
                #msg += "%s                s1=%-4i s2=%-4i s3=%-4i s4=%-4i %s\n"         %(' '*4,    s1[1],s2[1],s3[1],s4[1])
        return msg


