from __future__ import (nested_scopes, generators, division, absolute_import,
                        print_function, unicode_literals)
import sys

from .oes_objects import stressObject, strainObject


class CelasStressObject(stressObject):
    """
    @code
                              S T R E S S E S   I N   S C A L A R   S P R I N G S        ( C E L A S 2 )
        TIME         STRESS              TIME         STRESS              TIME         STRESS              TIME         STRESS
    0.0            0.0               5.000000E-02   0.0               1.000000E-01   0.0               1.500000E-01   0.0
    2.000000E-01   0.0               2.500000E-01   0.0               3.000000E-01   0.0               3.500000E-01   0.0
    @endcode
    """
    def __init__(self, dataCode, isSort1, iSubcase, dt=None):
        stressObject.__init__(self, dataCode, iSubcase)
        self.eType = {}
        self.elementName = self.dataCode['elementName']

        self.code = [self.formatCode, self.sortCode, self.sCode]
        self.stress = {}

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
        return (8, 'f')

    def deleteTransient(self, dt):
        del self.stress[dt]

    def getTransients(self):
        k = self.stress.keys()
        k.sort()
        return k

    def addNewTransient(self, dt):
        """initializes the transient variables"""
        self.elementName = self.dataCode['elementName']
        self.dt = dt
        self.stress[dt] = {}

    def addNewEid(self, dt, eid, out):
        (stress,) = out
        self.eType[eid] = self.elementName
        self.stress[eid] = stress

    def addNewEidSort1(self, dt, eid, out):
        if dt not in self.stress:
            self.addNewTransient(dt)
        (stress,) = out
        self.eType[eid] = self.elementName
        self.stress[dt][eid] = stress

    def addNewEidSort2(self, eid, dt, out):
        if dt not in self.stress:
            self.addNewTransient(dt)
        (stress,) = out
        self.eType[eid] = self.elementName
        self.stress[dt][eid] = stress

    def __reprTransient__(self):
        msg = '---CELASx STRESSES---\n'
        msg += '%-6s %6s ' % ('EID', 'eType')
        headers = ['stress']
        for header in headers:
            msg += '%10s ' % (header)
        msg += '\n'

        for dt, stress in sorted(self.stress.iteritems()):
            msg += '%s = %g\n' % (self.dataCode['name'], dt)
            for eid, istress in sorted(stress.iteritems()):
                msg += '%-6g %6s ' % (eid, self.eType[eid])
                if abs(istress) < 1e-6:
                    msg += '%10s ' % ('0')
                else:
                    msg += '%10g ' % (istress)
                msg += '\n'
        return msg

    def __repr__(self):
        #print "spring dt=%s" %(self.dt)
        if self.dt is not None:
            return self.__reprTransient__()

        msg = '---CELASx STRESSES---\n'
        msg += '%-8s %6s ' % ('EID', 'eType')
        headers = ['stress']
        for header in headers:
            msg += '%10s ' % (header)
        msg += '\n'
        #print "self.code = ",self.code
        for eid, istress in sorted(self.stress.iteritems()):
            #print "eid=",eid
            #print "eType",self.eType
            msg += '%-8i %6s ' % (eid, self.eType[eid])
            if abs(istress) < 1e-6:
                msg += '%10s ' % ('0')
            else:
                msg += '%10i ' % (istress)
            msg += '\n'
            #msg += "eid=%-4s eType=%s axial=%-4i torsion=%-4i\n" %(eid,self.eType,axial,torsion)
        return msg


class CelasStrainObject(strainObject):
    def __init__(self, dataCode, isSort1, iSubcase, dt=None):
        strainObject.__init__(self, dataCode, iSubcase)
        self.eType = {}
        self.elementName = self.dataCode['elementName']

        self.code = [self.formatCode, self.sortCode, self.sCode]

        self.isTransient = False
        self.strain = {}

        self.dt = dt
        if isSort1:
            if dt is not None:
                #self.add = self.addSort1
                self.addNewEid = self.addNewEidSort1
        else:
            assert dt is not None
            #self.add = self.addSort2
            self.addNewEid = self.addNewEidSort2
        ###

    def getLength(self):
        return (8, 'f')

    def deleteTransient(self, dt):
        del self.strain[dt]

    def getTransients(self):
        k = self.strain.keys()
        k.sort()
        return k

    def addNewTransient(self, dt):
        """
        initializes the transient variables
        """
        self.strain[dt] = {}

    def addNewEid(self, dt, eid, out):
        (strain,) = out
        assert eid >= 0
        #self.eType = self.eType
        self.eType[eid] = self.elementName
        self.strain[eid] = strain

    def addNewEidSort1(self, dt, eid, out):
        #print out
        (strain,) = out
        assert eid >= 0

        self.eType[eid] = self.elementType
        self.strain[dt][eid] = strain

    def __repr__(self):
        if self.dt is not None:
            return self.__reprTransient__()

        msg = '---CELASx STRAINS---\n'
        msg += '%-8s %6s ' % ('EID', 'eType')
        headers = ['strain']
        for header in headers:
            msg += '%8s ' % (header)
        msg += '\n'

        for eid, strain in sorted(self.strain.iteritems()):
            #strain = self.strain[eid]
            msg += '%-8i %6s ' % (eid, self.eType[eid])

            if abs(strain) < 1e-7:
                msg += '%8s ' % ('0')
            else:
                msg += '%8.3g ' % (strain)
            msg += '\n'
        return msg
