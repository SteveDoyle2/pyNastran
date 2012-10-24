from __future__ import (nested_scopes, generators, division, absolute_import,
                        print_function, unicode_literals)
import sys

from .oes_objects import stressObject


class NonlinearGapStressObject(stressObject):
    """
    """
    def __init__(self, dataCode, isSort1, iSubcase, dt=None):
        stressObject.__init__(self, dataCode, iSubcase)
        self.eType = {}
        self.elementName = self.dataCode['elementName']

        self.code = [self.formatCode, self.sortCode, self.sCode]
        self.compX = {}
        self.shearY = {}
        self.shearZ = {}
        self.axialU = {}
        self.shearV = {}
        self.shearW = {}
        self.slipV = {}
        self.slipW = {}

        self.dt = dt
        if isSort1:
            if dt is not None:
                #self.add = self.addSort1
                self.addNewEid = self.addNewEidSort1
        else:
            assert dt is not None
            #self.add = self.addSort2
            self.addNewEid = self.addNewEidSort2

    def get_stats(self):
        nelements = len(self.eType)
        eTypes = list(set(self.eType.values()))
        msg = self.get_data_code()
        if self.dt is not None:  # transient
            ntimes = len(self.compX)
            msg.append('  type=%s ntimes=%s nelements=%s\n'
                       % (self.__class__.__name__, ntimes, nelements))
        else:
            msg.append('  type=%s nelements=%s\n' % (self.__class__.__name__,
                                                     nelements))
        msg.append('  eType, compX, shearY, shearZ, axialU, shearV, shearW, slipV, slipW\n')
        msg.append('  eTypes = %s\n' %(', '.join(eTypes)))
        return msg

    def getLength(self):
        return (8, 'f')

    def deleteTransient(self, dt):
        del self.compX[dt]

    def getTransients(self):
        k = self.compX.keys()
        k.sort()
        return k

    def addNewTransient(self, dt):
        """initializes the transient variables"""
        self.elementName = self.dataCode['elementName']
        self.dt = dt
        self.compX[dt] = {}
        self.shearY[dt] = {}
        self.shearZ[dt] = {}
        self.axialU[dt] = {}
        self.shearV[dt] = {}
        self.shearW[dt] = {}
        self.slipV[dt] = {}
        self.slipW[dt] = {}


    def addNewEid(self, dt, eid, cpx, shy, shz, au, shv, shw, slv, slp, form1, form2):
        (stress,) = out
        self.eType[eid] = self.elementName
        self.compX[eid] = cpx
        self.shearY[eid] = shy
        self.shearZ[eid] = shz
        self.axialU[eid] = au
        self.shearV[eid] = shv
        self.shearW[eid] = shw
        self.slipV[eid] = slv
        self.slipW[eid] = slp

    def addNewEidSort1(self, dt, eid, cpx, shy, shz, au, shv, shw, slv, slp, form1, form2):
        if dt not in self.compX:
            self.addNewTransient(dt)
        self.eType[eid] = self.elementName
        self.compX[dt][eid] = cpx
        self.shearY[dt][eid] = shy
        self.shearZ[dt][eid] = shz
        self.axialU[dt][eid] = au
        self.shearV[dt][eid] = shv
        self.shearW[dt][eid] = shw
        self.slipV[dt][eid] = slv
        self.slipW[dt][eid] = slp

    def addNewEidSort2(self, eid, dt, cpx, shy, shz, au, shv, shw, slv, slp, form1, form2):
        if dt not in self.compX:
            self.addNewTransient(dt)
        (stress,) = out
        self.eType[eid] = self.elementName
        self.compX[dt][eid] = cpx
        self.shearY[dt][eid] = shy
        self.shearZ[dt][eid] = shz
        self.axialU[dt][eid] = au
        self.shearV[dt][eid] = shv
        self.shearW[dt][eid] = shw
        self.slipV[dt][eid] = slv
        self.slipW[dt][eid] = slp

    def __reprTransient__(self):
        raise NotImplementedError('GAPNL')
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
        raise NotImplementedError('GAPNL')
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
