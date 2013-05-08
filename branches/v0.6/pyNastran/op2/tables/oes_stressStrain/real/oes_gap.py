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

from .oes_objects import StressObject


class NonlinearGapStressObject(StressObject):

    def __init__(self, data_code, is_sort1, isubcase, dt=None):
        StressObject.__init__(self, data_code, isubcase)
        self.eType = {}
        self.element_name = self.data_code['element_name']

        self.code = [self.format_code, self.sort_code, self.s_code]
        self.compX = {}
        self.shearY = {}
        self.shearZ = {}
        self.axialU = {}
        self.shearV = {}
        self.shearW = {}
        self.slipV = {}
        self.slipW = {}

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

    def delete_transient(self, dt):
        del self.compX[dt]

    def get_transients(self):
        k = self.compX.keys()
        k.sort()
        return k

    def add_new_transient(self, dt):
        """initializes the transient variables"""
        self.element_name = self.data_code['element_name']
        self.dt = dt
        self.compX[dt] = {}
        self.shearY[dt] = {}
        self.shearZ[dt] = {}
        self.axialU[dt] = {}
        self.shearV[dt] = {}
        self.shearW[dt] = {}
        self.slipV[dt] = {}
        self.slipW[dt] = {}


    def add_new_eid(self, dt, eid, cpx, shy, shz, au, shv, shw, slv, slp, form1, form2):
        (stress,) = out
        self.eType[eid] = self.element_name
        self.compX[eid] = cpx
        self.shearY[eid] = shy
        self.shearZ[eid] = shz
        self.axialU[eid] = au
        self.shearV[eid] = shv
        self.shearW[eid] = shw
        self.slipV[eid] = slv
        self.slipW[eid] = slp

    def add_new_eid_sort1(self, dt, eid, cpx, shy, shz, au, shv, shw, slv, slp, form1, form2):
        if dt not in self.compX:
            self.add_new_transient(dt)
        self.eType[eid] = self.element_name
        self.compX[dt][eid] = cpx
        self.shearY[dt][eid] = shy
        self.shearZ[dt][eid] = shz
        self.axialU[dt][eid] = au
        self.shearV[dt][eid] = shv
        self.shearW[dt][eid] = shw
        self.slipV[dt][eid] = slv
        self.slipW[dt][eid] = slp

    def add_new_eid_sort2(self, eid, dt, cpx, shy, shz, au, shv, shw, slv, slp, form1, form2):
        if dt not in self.compX:
            self.add_new_transient(dt)
        (stress,) = out
        self.eType[eid] = self.element_name
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
            msg += '%10s ' % header
        msg += '\n'

        for dt, stress in sorted(self.stress.iteritems()):
            msg += '%s = %g\n' % (self.data_code['name'], dt)
            for eid, istress in sorted(stress.iteritems()):
                msg += '%-6g %6s ' % (eid, self.eType[eid])
                if abs(istress) < 1e-6:
                    msg += '%10s ' % '0'
                else:
                    msg += '%10g ' % istress
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
            msg += '%10s ' % header
        msg += '\n'
        #print "self.code = ",self.code
        for eid, istress in sorted(self.stress.iteritems()):
            #print "eid=",eid
            #print "eType",self.eType
            msg += '%-8i %6s ' % (eid, self.eType[eid])
            if abs(istress) < 1e-6:
                msg += '%10s ' % '0'
            else:
                msg += '%10i ' % istress
            msg += '\n'
            #msg += "eid=%-4s eType=%s axial=%-4i torsion=%-4i\n" %(eid,self.eType,axial,torsion)
        return msg
