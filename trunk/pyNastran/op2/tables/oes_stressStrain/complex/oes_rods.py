from __future__ import (nested_scopes, generators, division, absolute_import,
                        print_function, unicode_literals)
#import sys

from ..real.oes_objects import stressObject, strainObject


class ComplexRodDamperObject(stressObject):
    def __init__(self, dataCode, isSort1, iSubcase, dt=None):
        stressObject.__init__(self, dataCode, iSubcase)
        self.eType = 'CBUSH'

        self.code = [self.formatCode, self.sortCode, self.sCode]
        self.axial = {}
        self.torsion = {}


class ComplexRodStressObject(stressObject):
    """
    """
    def __init__(self, dataCode, isSort1, iSubcase, dt=None):
        stressObject.__init__(self, dataCode, iSubcase)
        self.eType = 'CROD'

        self.code = [self.formatCode, self.sortCode, self.sCode]
        self.axial = {}
        self.torsion = {}

        self.dt = dt
        if isSort1:
            if dt is not None:
                #self.add = self.addSort1
                self.addNewEid = self.addNewEidSort1
        else:
            assert dt is not None
            #self.add = self.addSort2
            self.addNewEid = self.addNewEidSort2

    def getLength(self):
        return (20, 'ffff')

    def addF06Data(self, data, transient):
        if transient is None:
            for line in data:
                (eid, axial, torsion) = line
                self.axial[eid] = axial
                self.torsion[eid] = torsion
            return

        (dtName, dt) = transient
        self.dt = dt
        self.dataCode['name'] = dtName
        if dt not in self.s1:
            self.updateDt(self.dataCode, dt)

        for line in data:
            (eid, axial, torsion) = line
            self.axial[dt][eid] = axial
            self.torsion[dt][eid] = torsion

    def deleteTransient(self, dt):
        del self.axial[dt]
        del self.torsion[dt]

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
        self.torsion[dt] = {}

    def addNewEid(self, dt, eid, out):
        #print "Rod Stress add..."
        (axial, torsion) = out
        assert isinstance(eid, int)
        self.axial[eid] = axial
        self.torsion[eid] = torsion

    def addNewEidSort1(self, dt, eid, out):
        (axial, torsion) = out

        if dt not in self.axial:
            self.addNewTransient(dt)
        self.axial[dt][eid] = axial
        self.torsion[dt][eid] = torsion

    def addNewEidSort2(self, eid, dt, out):
        (axial, torsion) = out

        if dt not in self.axial:
            self.addNewTransient(dt)
        self.axial[dt][eid] = axial
        self.torsion[dt][eid] = torsion

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
                msg += '%-6i %6s ' % (eid, self.eType)
                vals = [axial, torsion]
                for val in vals:
                    if abs(val) < 1e-6:
                        msg += '%10s ' % ('0')
                    else:
                        msg += '%10i ' % (val)
                msg += '\n'
                #msg += "eid=%-4s eType=%s axial=%-4i torsion=%-4i\n" %(eid,self.eType,axial,torsion)
        return msg

    def writeF06(self, header, pageStamp, pageNum=1, f=None, isMagPhase=False):
        return 'ComplexRodStress writeF06 not implemented...', pageNum
        raise NotImplementedError()
        if self.nonlinearFactor is not None:
            return self.writeF06Transient(header, pageStamp, pageNum, f, isMagPhase)

        msg = header + ['                                     S T R E S S E S   I N   R O D   E L E M E N T S      ( C R O D )\n',
                        '       ELEMENT       AXIAL       SAFETY      TORSIONAL     SAFETY       ELEMENT       AXIAL       SAFETY      TORSIONAL     SAFETY\n',
                        '         ID.        STRESS       MARGIN        STRESS      MARGIN         ID.        STRESS       MARGIN        STRESS      MARGIN\n']
        out = []
        for eid in sorted(self.axial):
            axial = self.axial[eid]
            torsion = self.torsion[eid]
            (vals2, isAllZeros) = self.writeFloatsImag13E([axial, torsion])
            (axial, torsion) = vals2
            out.append([eid, axial, torsion])

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
        return(''.join(msg), pageNum)

    def writeF06Transient(self, header, pageStamp, pageNum=1, f=None, isMagPhase=False):
        return 'ComplexRodStress writeF06Transient not implemented...', pageNum
        raise NotImplementedError()
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
                torsion = self.torsion[dt][eid]

                (vals2, isAllZeros) = self.writeFloatsImag13E([axial, torsion])
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
            msg += '%-6i %6s ' % (eid, self.eType)
            vals = [axial, torsion]
            for val in vals:
                if abs(val) < 1e-6:
                    msg += '%10s ' % ('0')
                else:
                    msg += '%10i ' % (val)
            msg += '\n'
            #msg += "eid=%-4s eType=%s axial=%-4i torsion=%-4i\n" %(eid,self.eType,axial,torsion)
        return msg


class ComplexRodStrainObject(strainObject):
    """
    """
    def __init__(self, dataCode, isSort1, iSubcase, dt=None):
        strainObject.__init__(self, dataCode, iSubcase)
        self.eType = 'CROD'  # {} # 'CROD/CONROD/CTUBE'

        self.code = [self.formatCode, self.sortCode, self.sCode]

        self.axial = {}
        self.torsion = {}

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
                (eid, axial, torsion) = line
                self.axial[eid] = axial
                self.torsion[eid] = torsion
            return

        (dtName, dt) = transient
        self.dt = dt
        self.dataCode['name'] = dtName
        if dt not in self.s1:
            self.updateDt(self.dataCode, dt)

        for line in data:
            (eid, axial, torsion) = line
            self.axial[dt][eid] = axial
            self.torsion[dt][eid] = torsion

    def deleteTransient(self, dt):
        del self.axial[dt]
        del self.torsion[dt]

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
        self.torsion[self.dt] = {}

    def addNewEid(self, dt, eid, out):
        (axial, torsion) = out
        assert eid >= 0
        #self.eType = self.eType
        self.axial[eid] = axial
        self.torsion[eid] = torsion

    def addNewEidSort1(self, dt, eid, out):
        (axial, torsion) = out
        assert eid >= 0
        #self.eType[eid] = self.elementType
        if dt not in self.axial:
            self.addNewTransient(dt)
        self.axial[dt][eid] = axial
        self.torsion[dt][eid] = torsion

    def addNewEidSort2(self, eid, dt, out):
        (axial, torsion) = out
        assert eid >= 0
        #self.eType[eid] = self.elementType
        if dt not in self.axial:
            self.addNewTransient(dt)
        self.axial[dt][eid] = axial
        self.torsion[dt][eid] = torsion

    def __reprTransient__(self):
        msg = '---ROD STRAINS---\n'
        msg += '%-6s %6s ' % ('EID', 'eType')
        headers = ['axial', 'torsion']
        for header in headers:
            msg += '%10s ' % (header)
        msg += '\n'

        for dt, axial in sorted(self.axial.iteritems()):
            msg += '%s = %g\n' % (self.dataCode['name'], dt)
            for eid in sorted(axial):
                axial = self.axial[dt][eid]
                torsion = self.torsion[dt][eid]
                msg += '%-6i %6s ' % (eid, self.eType)
                vals = [axial, torsion]
                for val in vals:
                    if abs(val) < 1e-6:
                        msg += '%10s ' % ('0')
                    else:
                        msg += '%10g ' % (val)
                msg += '\n'
                #msg += "eid=%-4s eType=%s axial=%-4i torsion=%-4i\n" %(eid,self.eType,axial,torsion)
        return msg

    def writeF06(self, header, pageStamp, pageNum=1, f=None, isMagPhase=False):
        return 'ComplexRodStrain writeF06 not implemented...', pageNum
        raise NotImplementedError()
        if self.dt is not None:
            return self.writeF06Transient(header, pageStamp, pageNum, f, isMagPhase)

        msg = header + ['                                       S T R A I N S   I N   R O D   E L E M E N T S      ( C R O D )\n',
                        '       ELEMENT       AXIAL       SAFETY      TORSIONAL     SAFETY       ELEMENT       AXIAL       SAFETY      TORSIONAL     SAFETY\n',
                        '         ID.        STRAIN       MARGIN        STRAIN      MARGIN         ID.        STRAIN       MARGIN        STRAIN      MARGIN\n']
        out = []
        for eid in sorted(self.axial):
            axial = self.axial[eid]
            torsion = self.torsion[eid]
            (vals2, isAllZeros) = self.writeFloatsImag13E([axial, torsion])
            (axial, torsion) = vals2
            out.append([eid, axial, torsion])

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
        return(''.join(msg), pageNum)

    def writeF06Transient(self, header, pageStamp, pageNum=1, f=None, isMagPhase=False):
        return 'ComplexRodStress writeF06Transient not implemented...', pageNum
        raise NotImplementedError()
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
                torsion = self.torsion[dt][eid]

                out.append([eid, axial, torsion])

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
            pageNum += 1
        return(''.join(msg), pageNum - 1)

    def __repr__(self):
        if self.dt is not None:
            return self.__reprTransient__()

        msg = '---ROD STRAINS---\n'
        msg += '%-6s %6s ' % ('EID', 'eType')
        headers = ['axial', 'torsion']
        for header in headers:
            msg += '%8s ' % (header)
        msg += '\n'

        for eid in sorted(self.axial):
            axial = self.axial[eid]
            torsion = self.torsion[eid]
            msg += '%-6i %6s ' % (eid, self.eType)
            vals = [axial, torsion]
            for val in vals:
                if abs(val) < 1e-7:
                    msg += '%8s ' % ('0')
                else:
                    msg += '%8.3g ' % (val)
            msg += '\n'
        return msg
