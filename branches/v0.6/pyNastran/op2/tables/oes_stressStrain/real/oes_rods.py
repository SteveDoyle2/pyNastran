## GNU Lesser General Public License
##
## Program pyNastran - a python interface to NASTRAN files
## Copyright (C) 2011-2012  Steven Doyle, Al Danial
##
## Authors and copyright holders of pyNastran
## Steven Doyle <mesheb82@gmail.com>
## Al Danial    <al.danial@gmail.com>
##
## This file is part of pyNastran.
##
## pyNastran is free software: you can redistribute it and/or modify
## it under the terms of the GNU Lesser General Public License as published by
## the Free Software Foundation, either version 3 of the License, or
## (at your option) any later version.
##
## pyNastran is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
##
## You should have received a copy of the GNU Lesser General Public License
## along with pyNastran.  If not, see <http://www.gnu.org/licenses/>.
##
from __future__ import (nested_scopes, generators, division, absolute_import,
                        print_function, unicode_literals)

from .oes_objects import StressObject, StrainObject
from pyNastran.f06.f06_formatting import writeFloats13E


class RodDamperObject(StressObject):
    def __init__(self, data_code, is_sort1, isubcase, dt=None):
        StressObject.__init__(self, data_code, isubcase)
        self.eType = 'CBUSH'

        self.code = [self.format_code, self.sort_code, self.s_code]
        self.axial = {}
        self.torsion = {}

    def get_stats(self):
        msg = self.get_data_code()
        eTypes = list(set(self.eType.values()))
        if self.nonlinear_factor is not None:  # transient
            ntimes = len(self.axial)
            a0 = self.stress.keys()[0]
            nelements = len(self.axial[a0])
            msg.append('  type=%s ntimes=%s nelements=%s\n'
                       % (self.__class__.__name__, ntimes, nelements))
        else:
            nelements = len(self.axial)
            msg.append('  type=%s nelements=%s\n' % (self.__class__.__name__,
                                                     nelements))
        msg.append('  eType, axial, torsion\n')
        msg.append('  eTypes = %s\n' %(', '.join(eTypes)))
        return msg


class RodStressObject(StressObject):
    """
    ::

      # format_code=1 stressCode=0
                                       S T R E S S E S   I N   R O D   E L E M E N T S      ( C R O D )
         ELEMENT       AXIAL       SAFETY      TORSIONAL     SAFETY       ELEMENT       AXIAL       SAFETY      TORSIONAL     SAFETY
           ID.        STRESS       MARGIN        STRESS      MARGIN         ID.        STRESS       MARGIN        STRESS      MARGIN
               1    5.000000E+03              0.0                               2    0.0                       0.0

      # format_code=1 stressCode=0
                                       S T R E S S E S   I N   R O D   E L E M E N T S      ( C R O D )
        ELEMENT       AXIAL       SAFETY      TORSIONAL     SAFETY       ELEMENT       AXIAL       SAFETY      TORSIONAL     SAFETY
          ID.        STRESS       MARGIN        STRESS      MARGIN         ID.        STRESS       MARGIN        STRESS      MARGIN
    """
    def __init__(self, data_code, is_sort1, isubcase, dt=None):
        StressObject.__init__(self, data_code, isubcase)
        self.eType = 'CROD'

        self.code = [self.format_code, self.sort_code, self.s_code]
        self.axial = {}
        self.torsion = {}

        self.MS_axial = {}
        self.MS_torsion = {}
        self.isImaginary = False

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
        msg = self.get_data_code()
        if self.dt is not None:  # transient
            ntimes = len(self.axial)
            a0 = self.axial.keys()[0]
            nelements = len(self.axial[a0])
            msg.append('  type=%s ntimes=%s nelements=%s\n'
                       % (self.__class__.__name__, ntimes, nelements))
        else:
            nelements = len(self.axial)
            msg.append('  type=%s nelements=%s\n' % (self.__class__.__name__,
                                                     nelements))
        msg.append('  eType, axial, torsion, MS_axial, MS_torsion\n')
        return msg

    def getLength(self):
        return (20, 'ffff')

    def add_f06_data(self, data, transient):
        if transient is None:
            for line in data:
                (eid, axial, MSa, torsion, MSt) = line
                if MSa is None:
                    MSa = 0.
                if MSt is None:
                    MSt = 0.
                self.axial[eid] = axial
                self.MS_axial[eid] = MSa
                self.torsion[eid] = torsion
                self.MS_torsion[eid] = MSt
            return

        (dtName, dt) = transient
        self.dt = dt
        self.data_code['name'] = dtName
        if dt not in self.axial:
            self.update_dt(self.data_code, dt)

        for line in data:
            (eid, axial, MSa, torsion, MSt) = line
            if MSa is None:
                MSa = 0.
            if MSt is None:
                MSt = 0.
            self.axial[dt][eid] = axial
            self.MS_axial[dt][eid] = MSa
            self.torsion[dt][eid] = torsion
            self.MS_torsion[dt][eid] = MSt

    def getLength(self):
        return (20, 'ffff')

    def delete_transient(self, dt):
        del self.axial[dt]
        del self.torsion[dt]
        del self.MS_axial[dt]
        del self.MS_torsion[dt]

    def get_transients(self):
        k = self.axial.keys()
        k.sort()
        return k

    def add_new_transient(self, dt):
        """
        initializes the transient variables
        """
        self.dt = dt
        self.axial[dt] = {}
        self.MS_axial[dt] = {}
        self.torsion[dt] = {}
        self.MS_torsion[dt] = {}

    def add_new_eid(self, dt, eid, out):
        #print "Rod Stress add..."
        (axial, SMa, torsion, SMt) = out
        assert isinstance(eid, int)
        self.axial[eid] = axial
        self.MS_axial[eid] = SMa
        self.torsion[eid] = torsion
        self.MS_torsion[eid] = SMt

    def add_new_eid_sort1(self, dt, eid, out):
        (axial, SMa, torsion, SMt) = out

        if dt not in self.axial:
            self.add_new_transient(dt)
        self.axial[dt][eid] = axial
        self.MS_axial[dt][eid] = SMa
        self.torsion[dt][eid] = torsion
        self.MS_torsion[dt][eid] = SMt

    def add_new_eid_sort2(self, eid, dt, out):
        (axial, SMa, torsion, SMt) = out

        if dt not in self.axial:
            self.add_new_transient(dt)
        self.axial[dt][eid] = axial
        self.MS_axial[dt][eid] = SMa
        self.torsion[dt][eid] = torsion
        self.MS_torsion[dt][eid] = SMt

    def __reprTransient__(self):
        msg = '---ROD STRESSES---\n'
        msg += '%-6s %6s ' % ('EID', 'eType')
        headers = ['axial', 'torsion', 'MS_axial', 'MS_torsion']
        for header in headers:
            msg += '%10s ' % header
        msg += '\n'

        for dt, axial in sorted(self.axial.iteritems()):
            msg += '%s = %g\n' % (self.data_code['name'], dt)
            for eid in sorted(axial):
                axial = self.axial[dt][eid]
                torsion = self.torsion[dt][eid]
                SMa = self.MS_axial[dt][eid]
                SMt = self.MS_torsion[dt][eid]
                msg += '%-6i %6s ' % (eid, self.eType)
                vals = [axial, torsion, SMa, SMt]
                for val in vals:
                    if abs(val) < 1e-6:
                        msg += '%10s ' % '0'
                    else:
                        msg += '%10i ' % val
                msg += '\n'
                #msg += "eid=%-4s eType=%s axial=%-4i torsion=%-4i\n" %(eid,self.eType,axial,torsion)
        return msg

    def write_f06(self, header, pageStamp, pageNum=1, f=None, is_mag_phase=False):
        if self.nonlinear_factor is not None:
            return self._write_f06_transient(header, pageStamp, pageNum, f)

        msg = header + ['                                     S T R E S S E S   I N   R O D   E L E M E N T S      ( C R O D )\n',
                        '       ELEMENT       AXIAL       SAFETY      TORSIONAL     SAFETY       ELEMENT       AXIAL       SAFETY      TORSIONAL     SAFETY\n',
                        '         ID.        STRESS       MARGIN        STRESS      MARGIN         ID.        STRESS       MARGIN        STRESS      MARGIN\n']
        out = []
        for eid in sorted(self.axial):
            axial = self.axial[eid]
            MSa = self.MS_axial[eid]
            torsion = self.torsion[eid]
            MSt = self.MS_torsion[eid]
            (vals2, isAllZeros) = writeFloats13E([axial, torsion])
            (axial, torsion) = vals2
            out.append([eid, axial, MSa, torsion, MSt])

        nOut = len(out)
        nWrite = nOut
        if nOut % 2 == 1:
            nWrite = nOut - 1
        for i in xrange(0, nWrite, 2):
            #print i,out[i:]
            outLine = '      %8i   %13s  %10.4E %13s  %10.4E   %8i   %13s  %10.4E %13s  %10.4E\n' % (tuple(out[i] + out[i + 1]))
            msg.append(outLine)

        if nOut % 2 == 1:
            outLine = '      %8i   %13s  %10.4E %13s  %10.4E\n' % (
                tuple(out[-1]))
            msg.append(outLine)
        msg.append(pageStamp + str(pageNum) + '\n')
        if f is not None:
            f.write(''.join(msg))
            msg = ['']
        return(''.join(msg), pageNum)

    def _write_f06_transient(self, header, pageStamp, pageNum=1, f=None, is_mag_phase=False):
        words = ['                                     S T R E S S E S   I N   R O D   E L E M E N T S      ( C R O D )\n',
                 '       ELEMENT       AXIAL       SAFETY      TORSIONAL     SAFETY       ELEMENT       AXIAL       SAFETY      TORSIONAL     SAFETY\n',
                 '         ID.        STRESS       MARGIN        STRESS      MARGIN         ID.        STRESS       MARGIN        STRESS      MARGIN\n']
        msg = []
        for dt, axials in sorted(self.axial.iteritems()):
            dtLine = '%14s = %12.5E\n' % (self.data_code['name'], dt)
            header[2] = dtLine
            msg += header + words
            out = []
            for eid in sorted(axials):
                axial = self.axial[dt][eid]
                MSa = self.MS_axial[dt][eid]
                torsion = self.torsion[dt][eid]
                MSt = self.MS_torsion[dt][eid]

                (vals2, isAllZeros) = writeFloats13E([axial, torsion])
                (axial, torsion) = vals2
                out.append([eid, axial, MSa, torsion, MSt])

            nOut = len(out)
            nWrite = nOut
            if nOut % 2 == 1:
                nWrite = nOut - 1
            for i in xrange(0, nWrite, 2):
                outLine = '      %8i   %13s  %10.4E %13s  %10.4E   %8i   %13s  %10.4E %13s  %10.4E\n' % (tuple(out[i] + out[i + 1]))
                msg.append(outLine)

            if nOut % 2 == 1:
                outLine = '      %8i   %13s  %10.4E %13s  %10.4E\n' % (
                    tuple(out[-1]))
                msg.append(outLine)
            msg.append(pageStamp + str(pageNum) + '\n')
            if f is not None:
                f.write(''.join(msg))
                msg = ['']
            pageNum += 1
        return(''.join(msg), pageNum - 1)

    def __repr__(self):
        if self.nonlinear_factor is not None:
            return self.__reprTransient__()

        #print 'axial = ',self.axial
        msg = '---ROD STRESSES---\n'
        msg += '%-6s %6s ' % ('EID', 'eType')
        headers = ['axial', 'torsion', 'MS_axial', 'MS_torsion']
        for header in headers:
            msg += '%10s ' % header
        msg += '\n'
        #print "self.code = ",self.code
        for eid in sorted(self.axial):
            #print self.__dict__.keys()
            axial = self.axial[eid]
            torsion = self.torsion[eid]
            SMa = self.MS_axial[eid]
            SMt = self.MS_torsion[eid]
            msg += '%-6i %6s ' % (eid, self.eType)
            vals = [axial, torsion, SMa, SMt]
            for val in vals:
                if abs(val) < 1e-6:
                    msg += '%10s ' % '0'
                else:
                    msg += '%10i ' % val
            msg += '\n'
            #msg += "eid=%-4s eType=%s axial=%-4i torsion=%-4i\n" %(eid,self.eType,axial,torsion)
        return msg


class RodStrainObject(StrainObject):
    """
    ::

      # s_code=1
                                       S T R A I N S   I N   R O D   E L E M E N T S      ( C R O D )
      ELEMENT       AXIAL       SAFETY      TORSIONAL     SAFETY
        ID.        STRAIN       MARGIN        STRAIN      MARGIN

      # s_code=10
                                         S T R A I N S   I N   R O D   E L E M E N T S      ( C O N R O D )
      ELEMENT       AXIAL       SAFETY      TORSIONAL     SAFETY       ELEMENT       AXIAL       SAFETY      TORSIONAL     SAFETY
        ID.        STRAIN       MARGIN        STRAIN      MARGIN         ID.        STRAIN       MARGIN        STRAIN      MARGIN
         1001    1.000000E+00   1.0E+00    1.250000E+00   3.0E+00         1007    1.000000E+00   1.0E+00    1.250000E+00   3.0E+00
    """
    def __init__(self, data_code, is_sort1, isubcase, dt=None):
        StrainObject.__init__(self, data_code, isubcase)
        self.eType = 'CROD'  # {} # 'CROD/CONROD/CTUBE'

        self.code = [self.format_code, self.sort_code, self.s_code]

        self.axial = {}
        self.torsion = {}

        self.MS_axial = {}
        self.MS_torsion = {}
        self.isImaginary = False

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
        msg = self.get_data_code()
        if self.dt is not None:  # transient
            ntimes = len(self.axial)
            a0 = self.axial.keys()[0]
            nelements = len(self.axial[a0])
            msg.append('  type=%s ntimes=%s nelements=%s\n'
                       % (self.__class__.__name__, ntimes, nelements))
        else:
            nelements = len(self.axial)
            msg.append('  type=%s nelements=%s\n' % (self.__class__.__name__,
                                                     nelements))
        msg.append('  eType, axial, torsion, MS_axial, MS_torsion\n')
        return msg

    def add_f06_data(self, data, transient):
        if transient is None:
            for line in data:
                (eid, axial, MSa, torsion, MSt) = line
                if MSa is None:
                    MSa = 0.
                if MSt is None:
                    MSt = 0.
                self.axial[eid] = axial
                self.MS_axial[eid] = MSa
                self.torsion[eid] = torsion
                self.MS_torsion[eid] = MSt
            return

        (dtName, dt) = transient
        self.dt = dt
        self.data_code['name'] = dtName
        if dt not in self.axial:
            self.update_dt(self.data_code, dt)

        for line in data:
            (eid, axial, MSa, torsion, MSt) = line
            if MSa is None:
                MSa = 0.
            if MSt is None:
                MSt = 0.
            self.axial[dt][eid] = axial
            self.MS_axial[dt][eid] = MSa
            self.torsion[dt][eid] = torsion
            self.MS_torsion[dt][eid] = MSt

    def getLength(self):
        return (20, 'ffff')

    def delete_transient(self, dt):
        del self.axial[dt]
        del self.torsion[dt]
        del self.MS_axial[dt]
        del self.MS_torsion[dt]

    def get_transients(self):
        k = self.axial.keys()
        k.sort()
        return k

    def add_new_transient(self, dt):
        """
        initializes the transient variables
        """
        self.dt = dt
        self.axial[self.dt] = {}
        self.MS_axial[self.dt] = {}
        self.torsion[self.dt] = {}
        self.MS_torsion[self.dt] = {}

    def add_new_eid(self, dt, eid, out):
        (axial, SMa, torsion, SMt) = out
        assert eid >= 0
        #self.eType = self.eType
        self.axial[eid] = axial
        self.MS_axial[eid] = SMa
        self.torsion[eid] = torsion
        self.MS_torsion[eid] = SMt

    def add_new_eid_sort1(self, dt, eid, out):
        (axial, SMa, torsion, SMt) = out
        assert eid >= 0
        #self.eType[eid] = self.element_type
        if dt not in self.axial:
            self.add_new_transient(dt)
        self.axial[dt][eid] = axial
        self.MS_axial[dt][eid] = SMa
        self.torsion[dt][eid] = torsion
        self.MS_torsion[dt][eid] = SMt

    def add_new_eid_sort2(self, eid, dt, out):
        (axial, SMa, torsion, SMt) = out
        assert eid >= 0
        #self.eType[eid] = self.element_type
        if dt not in self.axial:
            self.add_new_transient(dt)
        self.axial[dt][eid] = axial
        self.MS_axial[dt][eid] = SMa
        self.torsion[dt][eid] = torsion
        self.MS_torsion[dt][eid] = SMt

    def __reprTransient__(self):
        msg = '---ROD STRAINS---\n'
        msg += '%-6s %6s ' % ('EID', 'eType')
        headers = ['axial', 'torsion', 'MS_axial', 'MS_torsion']
        for header in headers:
            msg += '%10s ' % header
        msg += '\n'

        for dt, axial in sorted(self.axial.iteritems()):
            msg += '%s = %g\n' % (self.data_code['name'], dt)
            for eid in sorted(axial):
                axial = self.axial[dt][eid]
                torsion = self.torsion[dt][eid]
                SMa = self.MS_axial[dt][eid]
                SMt = self.MS_torsion[dt][eid]
                msg += '%-6i %6s ' % (eid, self.eType)
                vals = [axial, torsion, SMa, SMt]
                for val in vals:
                    if abs(val) < 1e-6:
                        msg += '%10s ' % '0'
                    else:
                        msg += '%10g ' % val
                msg += '\n'
                #msg += "eid=%-4s eType=%s axial=%-4i torsion=%-4i\n" %(eid,self.eType,axial,torsion)
        return msg

    def write_f06(self, header, pageStamp, pageNum=1, f=None, is_mag_phase=False):
        if self.dt is not None:
            print("dt =", self.dt)
            return self._write_f06_transient(header, pageStamp, pageNum, f)

        msg = header + ['                                       S T R A I N S   I N   R O D   E L E M E N T S      ( C R O D )\n',
                        '       ELEMENT       AXIAL       SAFETY      TORSIONAL     SAFETY       ELEMENT       AXIAL       SAFETY      TORSIONAL     SAFETY\n',
                        '         ID.        STRAIN       MARGIN        STRAIN      MARGIN         ID.        STRAIN       MARGIN        STRAIN      MARGIN\n']
        out = []
        for eid in sorted(self.axial):
            axial = self.axial[eid]
            MSa = self.MS_axial[eid]
            torsion = self.torsion[eid]
            MSt = self.MS_torsion[eid]
            (vals2, isAllZeros) = writeFloats13E([axial, torsion])
            (axial, torsion) = vals2
            out.append([eid, axial, MSa, torsion, MSt])

        nOut = len(out)
        nWrite = nOut
        if nOut % 2 == 1:
            nWrite = nOut - 1
        for i in xrange(0, nWrite, 2):
            outLine = '      %8i   %13s  %10.4E %13s  %10.4E   %8i   %13s  %10.4E %13s  %10.4E\n' % (tuple(out[i] + out[i + 1]))
            msg.append(outLine)

        if nOut % 2 == 1:
            outLine = '      %8i   %13s  %10.4E %13s  %10.4E\n' % (
                tuple(out[-1]))
            msg.append(outLine)
        msg.append(pageStamp + str(pageNum) + '\n')
        if f is not None:
            f.write(''.join(msg))
            msg = ['']
        return(''.join(msg), pageNum)

    def _write_f06_transient(self, header, pageStamp, pageNum=1, f=None, is_mag_phase=False):
        words = ['                                       S T R A I N S   I N   R O D   E L E M E N T S      ( C R O D )\n',
                 '       ELEMENT       AXIAL       SAFETY      TORSIONAL     SAFETY       ELEMENT       AXIAL       SAFETY      TORSIONAL     SAFETY\n',
                 '         ID.        STRAIN       MARGIN        STRAIN      MARGIN         ID.        STRAIN       MARGIN        STRAIN      MARGIN\n']
        msg = []
        for dt, axials in sorted(self.axial.iteritems()):
            dtLine = '%14s = %12.5E\n' % (self.data_code['name'], dt)
            header[2] = dtLine
            msg += header + words
            out = []
            for eid in sorted(axials):
                axial = self.axial[dt][eid]
                MSa = self.MS_axial[dt][eid]
                torsion = self.torsion[dt][eid]
                MSt = self.MS_torsion[dt][eid]

                out.append([eid, axial, MSa, torsion, MSt])

            nOut = len(out)
            nWrite = nOut
            if nOut % 2 == 1:
                nWrite = nOut - 1
            for i in xrange(0, nWrite, 2):
                outLine = '      %8i   %13.6E  %10.4E %13.6E  %10.4E   %8i   %13.6E  %10.4E %13.6E  %10.4E\n' % (tuple(out[i] + out[i + 1]))
                msg.append(outLine)

            if nOut % 2 == 1:
                outLine = '      %8i   %13.6E  %10.4E %13.6E  %10.4E\n' % (
                    tuple(out[-1]))
                msg.append(outLine)
            msg.append(pageStamp + str(pageNum) + '\n')
            if f is not None:
                f.write(''.join(msg))
                msg = ['']
            pageNum += 1
        return(''.join(msg), pageNum - 1)

    def __repr__(self):
        if self.dt is not None:
            return self.__reprTransient__()

        msg = '---ROD STRAINS---\n'
        msg += '%-6s %6s ' % ('EID', 'eType')
        headers = ['axial', 'torsion', 'MS_tension', 'MS_compression']
        for header in headers:
            msg += '%8s ' % header
        msg += '\n'

        for eid in sorted(self.axial):
            axial = self.axial[eid]
            torsion = self.torsion[eid]
            SMa = self.MS_axial[eid]
            SMt = self.MS_torsion[eid]
            msg += '%-6i %6s ' % (eid, self.eType)
            vals = [axial, torsion, SMa, SMt]
            for val in vals:
                if abs(val) < 1e-7:
                    msg += '%8s ' % '0'
                else:
                    msg += '%8.3g ' % val
            msg += '\n'
        return msg