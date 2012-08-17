from __future__ import (nested_scopes, generators, division, absolute_import,
                        print_function, unicode_literals)
import sys

from .oes_objects import stressObject, strainObject


class ShearStressObject(stressObject):
    """
    # formatCode=1 sortCode=0 stressCode=0
                                   S T R E S S E S   I N   S H E A R   P A N E L S      ( C S H E A R )
    ELEMENT            MAX            AVG        SAFETY         ELEMENT            MAX            AVG        SAFETY
      ID.             SHEAR          SHEAR       MARGIN           ID.             SHEAR          SHEAR       MARGIN
        328        1.721350E+03   1.570314E+03   7.2E+01
    """
    def __init__(self, dataCode, isSort1, iSubcase, dt=None):
        stressObject.__init__(self, dataCode, iSubcase)
        self.eType = 'CSHEAR'

        self.code = [self.formatCode, self.sortCode, self.sCode]
        self.maxShear = {}
        self.avgShear = {}
        self.margin = {}

        self.getLength = self.getLength
        self.isImaginary = False
        #if dt is not None:
        #    self.addNewTransient = self.addNewTransient
        #    self.addNewEid       = self.addNewEidTransient
        #else:
        #    self.addNewEid = self.addNewEid
        ###

        self.dt = dt
        if isSort1:
            if dt is not None:
                self.add = self.addSort1
                self.addNewEid = self.addNewEidSort1
        else:
            assert dt is not None
            self.add = self.addSort2
            self.addNewEid = self.addNewEidSort2

    def getLength(self):
        return (16, 'fff')

    def deleteTransient(self, dt):
        del self.maxShear[dt]
        del self.avgShear[dt]
        del self.margin[dt]

    def getTransients(self):
        k = self.maxShear.keys()
        k.sort()
        return k

    def addNewTransient(self, dt):
        """
        initializes the transient variables
        """
        self.dt = dt
        self.maxShear[dt] = {}
        self.avgShear[dt] = {}
        self.margin[dt] = {}

    def addNewEid(self, dt, eid, out):
        #print "Rod Stress add..."
        (maxShear, avgShear, margin) = out
        assert isinstance(eid, int)
        self.maxShear = {}
        self.avgShear = {}
        self.margin = {}
        self.maxShear[eid] = maxShear
        self.avgShear[eid] = avgShear
        self.margin[eid] = margin

    def addNewEidSort1(self, dt, eid, out):
        (maxShear, avgShear, margin) = out
        if dt not in self.maxShear:
            self.addNewTransient(dt)
        assert isinstance(eid, int)
        assert eid >= 0
        self.maxShear[dt][eid] = maxShear
        self.avgShear[dt][eid] = avgShear
        self.margin[dt][eid] = margin

    def addNewEidSort1(self, eid, dt, out):
        (maxShear, avgShear, margin) = out
        if dt not in self.maxShear:
            self.addNewTransient(dt)
        assert isinstance(eid, int)
        assert eid >= 0
        self.maxShear[dt][eid] = maxShear
        self.avgShear[dt][eid] = avgShear
        self.margin[dt][eid] = margin

    def __reprTransient__(self):
        msg = '---TRANSIENT CSHEAR STRESSES---\n'
        msg += '%-6s %6s ' % ('EID', 'eType')
        headers = ['maxShear', 'avgShear', 'Margin']
        for header in headers:
            msg += '%10s ' % (header)
        msg += '\n'

        for dt, maxShears in sorted(self.maxShear.iteritems()):
            msg += '%s = %g\n' % (self.dataCode['name'], dt)
            for eid in sorted(maxShears):
                maxShear = self.maxShear[dt][eid]
                avgShear = self.avgShear[dt][eid]
                margin = self.margin[dt][eid]
                msg += '%-6i %6s ' % (eid, self.eType)
                vals = [maxShear, avgShear, margin]
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

    def __repr__(self):
        if self.dt is not None:
            return self.__reprTransient__()

        msg = '---CSHEAR STRESSES---\n'
        msg += '%-6s %6s ' % ('EID', 'eType')
        headers = ['maxShear', 'avgShear', 'margin']
        for header in headers:
            msg += '%10s ' % (header)
        msg += '\n'
        #print "self.code = ",self.code
        for eid in sorted(self.maxShear):
            #print self.__dict__.keys()
            maxShear = self.maxShear[eid]
            avgShear = self.avgShear[eid]
            margin = self.margin[eid]
            msg += '%-6i %6s ' % (eid, self.eType)
            vals = [maxShear, avgShear, margin]
            for val in vals:
                if abs(val) < 1e-6:
                    msg += '%10s ' % ('0')
                else:
                    msg += '%10i ' % (val)
                ###
            msg += '\n'
            #msg += "eid=%-4s eType=%s axial=%-4i torsion=%-4i\n" %(eid,self.eType,axial,torsion)
        return msg


class ShearStrainObject(strainObject):
    """
    """
    def __init__(self, dataCode, isSort1, iSubcase, dt=None):
        strainObject.__init__(self, dataCode, iSubcase)
        self.eType = 'CSHEAR'
        raise Exception('not supported...CSHEAR strain')
        self.code = [self.formatCode, self.sortCode, self.sCode]
        self.maxShear = {}
        self.avgShear = {}
        self.margin = {}

        self.dt = dt
        if isSort1:
            if dt is not None:
                self.add = self.addSort1
                self.addNewEid = self.addNewEidSort1
        else:
            assert dt is not None
            self.add = self.addSort2
            self.addNewEid = self.addNewEidSort2

    def getLength(self):
        return (16, 'fff')

    def deleteTransient(self, dt):
        del self.maxShear[dt]
        del self.avgShear[dt]
        del self.margin[dt]

    def getTransients(self):
        k = self.maxShear.keys()
        k.sort()
        return k

    def addNewTransient(self, dt):
        """
        initializes the transient variables
        @note make sure you set self.dt first
        """
        self.dt = dt
        self.maxShear[dt] = {}
        self.avgShear[dt] = {}
        self.margin[dt] = {}

    def addNewEid(self, dt, eid, out):
        (axial, SMa, torsion, SMt) = out
        #print "Rod Strain add..."
        assert eid >= 0
        #self.eType = self.eType
        self.maxShearl[eid] = axial
        self.avgShear[eid] = SMa
        self.margin[eid] = torsion

    def addNewEidSort1(self, dt, eid, out):
        (maxShear, avgShear, margin) = out
        if dt not in self.maxShear:
            self.addNewTransient(dt)
        assert eid >= 0

        #self.eType[eid] = self.elementType
        self.maxShear[dt][eid] = maxShear
        self.avgShear[dt][eid] = avgShear
        self.margin[dt][eid] = margin

    def addNewEidSort2(self, eid, dt, out):
        (maxShear, avgShear, margin) = out
        if dt not in self.maxShear:
            self.addNewTransient(dt)
        assert eid >= 0

        #self.eType[eid] = self.elementType
        self.maxShear[dt][eid] = maxShear
        self.avgShear[dt][eid] = avgShear
        self.margin[dt][eid] = margin

    def __reprTransient__(self):
        msg = '---TRANSIENT CSHEAR STRAINS---\n'
        msg += '%-6s %6s ' % ('EID', 'eType')
        headers = ['maxShear', 'avgShear', 'Margin']
        for header in headers:
            msg += '%10s ' % (header)
        msg += '\n'

        for dt, maxShears in sorted(self.maxShear.iteritems()):
            msg += '%s = %g\n' % (self.dataCode['name'], dt)
            for eid in sorted(maxShears):
                maxShear = self.maxShear[dt][eid]
                avgShear = self.avgShear[dt][eid]
                margin = self.margin[dt][eid]
                msg += '%-6i %6s ' % (eid, self.eType)
                vals = [maxShear, avgShear, margin]
                for val in vals:
                    if abs(val) < 1e-6:
                        msg += '%10s ' % ('0')
                    else:
                        msg += '%10g ' % (val)
                msg += '\n'
                #msg += "eid=%-4s eType=%s axial=%-4i torsion=%-4i\n" %(eid,self.eType,axial,torsion)
        return msg

    def __repr__(self):
        if self.dt is not None:
            return self.__reprTransient__()

        msg = '---CSHEAR STRAINS---\n'
        msg += '%-6s %6s ' % ('EID', 'eType')
        headers = ['maxShear', 'avgShear', 'margin']
        for header in headers:
            msg += '%10s ' % (header)
        msg += '\n'

        #print "self.code = ",self.code
        for eid in sorted(self.maxShear):
            #print self.__dict__.keys()
            maxShear = self.maxShear[eid]
            avgShear = self.avgShear[eid]
            margin = self.margin[eid]
            msg += '%-6i %6s ' % (eid, self.eType)
            vals = [maxShear, avgShear, margin]

            for val in vals:
                if abs(val) < 1e-7:
                    msg += '%8s ' % ('0')
                else:
                    msg += '%8.3g ' % (val)
            msg += '\n'
        return msg
