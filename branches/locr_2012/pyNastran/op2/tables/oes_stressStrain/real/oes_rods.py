from __future__ import (nested_scopes, generators, division, absolute_import,
                        print_function, unicode_literals)
import sys

from .oes_objects import stressObject, strainObject


class RodDamperObject(stressObject):
    def __init__(self, dataCode, isSort1, iSubcase, dt=None):
        stressObject.__init__(self, dataCode, iSubcase)
        self.eType = 'CBUSH'

        self.code = [self.formatCode, self.sortCode, self.sCode]
        self.axial = {}
        self.torsion = {}


class RodStressObject(stressObject):
    """
    @code
    # formatCode=1 stressCode=0
                                     S T R E S S E S   I N   R O D   E L E M E N T S      ( C R O D )
       ELEMENT       AXIAL       SAFETY      TORSIONAL     SAFETY       ELEMENT       AXIAL       SAFETY      TORSIONAL     SAFETY
         ID.        STRESS       MARGIN        STRESS      MARGIN         ID.        STRESS       MARGIN        STRESS      MARGIN
             1    5.000000E+03              0.0                               2    0.0                       0.0

    # formatCode=1 stressCode=0
                                     S T R E S S E S   I N   R O D   E L E M E N T S      ( C R O D )
      ELEMENT       AXIAL       SAFETY      TORSIONAL     SAFETY       ELEMENT       AXIAL       SAFETY      TORSIONAL     SAFETY
        ID.        STRESS       MARGIN        STRESS      MARGIN         ID.        STRESS       MARGIN        STRESS      MARGIN
    @endcode
    """
    def __init__(self, dataCode, isSort1, iSubcase, dt=None):
        stressObject.__init__(self, dataCode, iSubcase)
        self.eType = 'CROD'

        self.code = [self.formatCode, self.sortCode, self.sCode]
        self.axial = {}
        self.torsion = {}

        self.MS_axial = {}
        self.MS_torsion = {}
        self.isImaginary = False

        self.dt = dt
        if isSort1:
            if dt is not None:
                #self.add = self.addSort1
                self.addNewEid = self.addNewEidSort1
            ###
        else:
            assert dt is not None
            #self.add = self.addSort2
            self.addNewEid = self.addNewEidSort2
        ###

    def getLength(self):
        return (20, 'ffff')

    def addF06Data(self, data, transient):
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
        self.dataCode['name'] = dtName
        if dt not in self.s1:
            self.updateDt(self.dataCode, dt)

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

    def deleteTransient(self, dt):
        del self.axial[dt]
        del self.torsion[dt]
        del self.MS_axial[dt]
        del self.MS_torsion[dt]

    def getTransients(self):
        k = self.axial.keys()
        k.sort()
        return k

    def addNewTransient(self, dt):
        """
        initializes the transient variables
        """
        self.dt = dt
        self.axial[dt] = {}
        self.MS_axial[dt] = {}
        self.torsion[dt] = {}
        self.MS_torsion[dt] = {}

    def addNewEid(self, dt, eid, out):
        #print "Rod Stress add..."
        (axial, SMa, torsion, SMt) = out
        assert isinstance(eid, int)
        self.axial[eid] = axial
        self.MS_axial[eid] = SMa
        self.torsion[eid] = torsion
        self.MS_torsion[eid] = SMt

    def addNewEidSort1(self, dt, eid, out):
        (axial, SMa, torsion, SMt) = out

        if dt not in self.axial:
            self.addNewTransient(dt)
        self.axial[dt][eid] = axial
        self.MS_axial[dt][eid] = SMa
        self.torsion[dt][eid] = torsion
        self.MS_torsion[dt][eid] = SMt

    def addNewEidSort2(self, eid, dt, out):
        (axial, SMa, torsion, SMt) = out

        if dt not in self.axial:
            self.addNewTransient(dt)
        self.axial[dt][eid] = axial
        self.MS_axial[dt][eid] = SMa
        self.torsion[dt][eid] = torsion
        self.MS_torsion[dt][eid] = SMt

    def __reprTransient__(self):
        msg = '---ROD STRESSES---\n'
        msg += '%-6s %6s ' % ('EID', 'eType')
        headers = ['axial', 'torsion', 'MS_axial', 'MS_torsion']
        for header in headers:
            msg += '%10s ' % (header)
        msg += '\n'

        for dt, axial in sorted(self.axial.iteritems()):
            msg += '%s = %g\n' % (self.dataCode['name'], dt)
            for eid in sorted(axial):
                axial = self.axial[dt][eid]
                torsion = self.torsion[dt][eid]
                SMa = self.MS_axial[dt][eid]
                SMt = self.MS_torsion[dt][eid]
                msg += '%-6i %6s ' % (eid, self.eType)
                vals = [axial, torsion, SMa, SMt]
                for val in vals:
                    if abs(val) < 1e-6:
                        msg += '%10s ' % ('0')
                    else:
                        msg += '%10i ' % (val)
                    ###
                msg += '\n'
                #msg += "eid=%-4s eType=%s axial=%-4i torsion=%-4i\n" %(eid,self.eType,axial,torsion)
            ###
        return msg

    def writeF06(self, header, pageStamp, pageNum=1, f=None, isMagPhase=False):
        if self.nonlinearFactor is not None:
            return self.writeF06Transient(header, pageStamp, pageNum, f)

        msg = header + ['                                     S T R E S S E S   I N   R O D   E L E M E N T S      ( C R O D )\n',
                        '       ELEMENT       AXIAL       SAFETY      TORSIONAL     SAFETY       ELEMENT       AXIAL       SAFETY      TORSIONAL     SAFETY\n',
                        '         ID.        STRESS       MARGIN        STRESS      MARGIN         ID.        STRESS       MARGIN        STRESS      MARGIN\n']
        out = []
        for eid in sorted(self.axial):
            axial = self.axial[eid]
            MSa = self.MS_axial[eid]
            torsion = self.torsion[eid]
            MSt = self.MS_torsion[eid]
            (vals2, isAllZeros) = self.writeFloats13E([axial, torsion])
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

    def writeF06Transient(self, header, pageStamp, pageNum=1, f=None, isMagPhase=False):
        words = ['                                     S T R E S S E S   I N   R O D   E L E M E N T S      ( C R O D )\n',
                 '       ELEMENT       AXIAL       SAFETY      TORSIONAL     SAFETY       ELEMENT       AXIAL       SAFETY      TORSIONAL     SAFETY\n',
                 '         ID.        STRESS       MARGIN        STRESS      MARGIN         ID.        STRESS       MARGIN        STRESS      MARGIN\n']
        msg = []
        for dt, axials in sorted(self.axial.iteritems()):
            dtLine = '%14s = %12.5E\n' % (self.dataCode['name'], dt)
            header[2] = dtLine
            msg += header + words
            out = []
            for eid in sorted(axials):
                axial = self.axial[dt][eid]
                MSa = self.MS_axial[dt][eid]
                torsion = self.torsion[dt][eid]
                MSt = self.MS_torsion[dt][eid]

                (vals2, isAllZeros) = self.writeFloats13E([axial, torsion])
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
        if self.nonlinearFactor is not None:
            return self.__reprTransient__()

        #print 'axial = ',self.axial
        msg = '---ROD STRESSES---\n'
        msg += '%-6s %6s ' % ('EID', 'eType')
        headers = ['axial', 'torsion', 'MS_axial', 'MS_torsion']
        for header in headers:
            msg += '%10s ' % (header)
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
                    msg += '%10s ' % ('0')
                else:
                    msg += '%10i ' % (val)
                ###
            msg += '\n'
            #msg += "eid=%-4s eType=%s axial=%-4i torsion=%-4i\n" %(eid,self.eType,axial,torsion)
        return msg


class RodStrainObject(strainObject):
    """
    @code
    # sCode=1
                                     S T R A I N S   I N   R O D   E L E M E N T S      ( C R O D )
    ELEMENT       AXIAL       SAFETY      TORSIONAL     SAFETY
      ID.        STRAIN       MARGIN        STRAIN      MARGIN

    # sCode=10
                                       S T R A I N S   I N   R O D   E L E M E N T S      ( C O N R O D )
    ELEMENT       AXIAL       SAFETY      TORSIONAL     SAFETY       ELEMENT       AXIAL       SAFETY      TORSIONAL     SAFETY
      ID.        STRAIN       MARGIN        STRAIN      MARGIN         ID.        STRAIN       MARGIN        STRAIN      MARGIN
       1001    1.000000E+00   1.0E+00    1.250000E+00   3.0E+00         1007    1.000000E+00   1.0E+00    1.250000E+00   3.0E+00
    @endcode
    """
    def __init__(self, dataCode, isSort1, iSubcase, dt=None):
        strainObject.__init__(self, dataCode, iSubcase)
        self.eType = 'CROD'  # {} # 'CROD/CONROD/CTUBE'

        self.code = [self.formatCode, self.sortCode, self.sCode]

        self.axial = {}
        self.torsion = {}

        self.MS_axial = {}
        self.MS_torsion = {}
        self.isImaginary = False

        self.dt = dt
        if isSort1:
            if dt is not None:
                #self.add = self.addSort1
                self.addNewEid = self.addNewEidSort1
        else:
            assert dt is not None
            #self.add = self.addSort2
            self.addNewEid = self.addNewEidSort2

    def addF06Data(self, data, transient):
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
        self.dataCode['name'] = dtName
        if dt not in self.s1:
            self.updateDt(self.dataCode, dt)

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

    def deleteTransient(self, dt):
        del self.axial[dt]
        del self.torsion[dt]
        del self.MS_axial[dt]
        del self.MS_torsion[dt]

    def getTransients(self):
        k = self.axial.keys()
        k.sort()
        return k

    def addNewTransient(self, dt):
        """
        initializes the transient variables
        """
        self.dt = dt
        self.axial[self.dt] = {}
        self.MS_axial[self.dt] = {}
        self.torsion[self.dt] = {}
        self.MS_torsion[self.dt] = {}

    def addNewEid(self, dt, eid, out):
        (axial, SMa, torsion, SMt) = out
        assert eid >= 0
        #self.eType = self.eType
        self.axial[eid] = axial
        self.MS_axial[eid] = SMa
        self.torsion[eid] = torsion
        self.MS_torsion[eid] = SMt

    def addNewEidSort1(self, dt, eid, out):
        (axial, SMa, torsion, SMt) = out
        assert eid >= 0
        #self.eType[eid] = self.elementType
        if dt not in self.axial:
            self.addNewTransient(dt)
        self.axial[dt][eid] = axial
        self.MS_axial[dt][eid] = SMa
        self.torsion[dt][eid] = torsion
        self.MS_torsion[dt][eid] = SMt

    def addNewEidSort2(self, eid, dt, out):
        (axial, SMa, torsion, SMt) = out
        assert eid >= 0
        #self.eType[eid] = self.elementType
        if dt not in self.axial:
            self.addNewTransient(dt)
        self.axial[dt][eid] = axial
        self.MS_axial[dt][eid] = SMa
        self.torsion[dt][eid] = torsion
        self.MS_torsion[dt][eid] = SMt

    def __reprTransient__(self):
        msg = '---ROD STRAINS---\n'
        msg += '%-6s %6s ' % ('EID', 'eType')
        headers = ['axial', 'torsion', 'MS_axial', 'MS_torsion']
        for header in headers:
            msg += '%10s ' % (header)
        msg += '\n'

        for dt, axial in sorted(self.axial.iteritems()):
            msg += '%s = %g\n' % (self.dataCode['name'], dt)
            for eid in sorted(axial):
                axial = self.axial[dt][eid]
                torsion = self.torsion[dt][eid]
                SMa = self.MS_axial[dt][eid]
                SMt = self.MS_torsion[dt][eid]
                msg += '%-6i %6s ' % (eid, self.eType)
                vals = [axial, torsion, SMa, SMt]
                for val in vals:
                    if abs(val) < 1e-6:
                        msg += '%10s ' % ('0')
                    else:
                        msg += '%10g ' % (val)
                msg += '\n'
                #msg += "eid=%-4s eType=%s axial=%-4i torsion=%-4i\n" %(eid,self.eType,axial,torsion)
        return msg

    def writeF06(self, header, pageStamp, pageNum=1, f=None, isMagPhase=False):
        if self.dt is not None:
            return self.writeF06Transient(header, pageStamp, pageNum, f)

        msg = header + ['                                       S T R A I N S   I N   R O D   E L E M E N T S      ( C R O D )\n',
                        '       ELEMENT       AXIAL       SAFETY      TORSIONAL     SAFETY       ELEMENT       AXIAL       SAFETY      TORSIONAL     SAFETY\n',
                        '         ID.        STRAIN       MARGIN        STRAIN      MARGIN         ID.        STRAIN       MARGIN        STRAIN      MARGIN\n']
        out = []
        for eid in sorted(self.axial):
            axial = self.axial[eid]
            MSa = self.MS_axial[eid]
            torsion = self.torsion[eid]
            MSt = self.MS_torsion[eid]
            (vals2, isAllZeros) = self.writeFloats13E([axial, torsion])
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

    def writeF06Transient(self, header, pageStamp, pageNum=1, f=None, isMagPhase=False):
        words = ['                                       S T R A I N S   I N   R O D   E L E M E N T S      ( C R O D )\n',
                 '       ELEMENT       AXIAL       SAFETY      TORSIONAL     SAFETY       ELEMENT       AXIAL       SAFETY      TORSIONAL     SAFETY\n',
                 '         ID.        STRAIN       MARGIN        STRAIN      MARGIN         ID.        STRAIN       MARGIN        STRAIN      MARGIN\n']
        msg = []
        for dt, axials in sorted(self.axial.iteritems()):
            dtLine = '%14s = %12.5E\n' % (self.dataCode['name'], dt)
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
            msg += '%8s ' % (header)
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
                    msg += '%8s ' % ('0')
                else:
                    msg += '%8.3g ' % (val)
            msg += '\n'
        return msg
