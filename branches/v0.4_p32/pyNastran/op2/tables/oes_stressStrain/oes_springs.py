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
import sys
from pyNastran.op2.tables.oes_stressStrain.oes_objects import stressObject,strainObject #,array
from pyNastran.op2.op2Errors import *

class celasStressObject(stressObject):
    """
                              S T R E S S E S   I N   S C A L A R   S P R I N G S        ( C E L A S 2 )
        TIME         STRESS              TIME         STRESS              TIME         STRESS              TIME         STRESS
    0.0            0.0               5.000000E-02   0.0               1.000000E-01   0.0               1.500000E-01   0.0
    2.000000E-01   0.0               2.500000E-01   0.0               3.000000E-01   0.0               3.500000E-01   0.0
    """
    def __init__(self,dataCode,iSubcase,dt=None):
        stressObject.__init__(self,dataCode,iSubcase)
        self.eType = {}
        self.elementName = self.dataCode['elementName']

        self.code = [self.formatCode,self.sortCode,self.sCode]
        self.stress = {}
        if self.code in [[1,0,0]]:
            self.getLength = self.getLength_format1_sort0
            #self.isImaginary = False
            if dt is not None:
                self.addNewTransient = self.addNewTransient_format1_sort0
                self.addNewEid       = self.addNewEidTransient_format1_sort0
                self.isTransient = True
            else:
                self.addNewEid = self.addNewEid_format1_sort0
            ###
        elif 0: #self.code==[2,1,0]:
            self.getLength       = self.getLength_format1_sort0
            self.addNewTransient = self.addNewTransient_format2_sort1
            self.addNewEid       = self.addNewEidTransient_format2_sort1
            self.isImaginary = True
        else:
            raise InvalidCodeError('rodStress - get the format/sort/stressCode=%s' %(self.code))
        ###
        if dt is not None:
            #self.isTransient = True
            self.dt = self.nonlinearFactor
            self.addNewTransient()
        ###

    def getLength_format1_sort0(self):
        return (8,'if')

    def deleteTransient(self,dt):
        del self.stress[dt]

    def getTransients(self):
        k = self.stress.keys()
        k.sort()
        return k

    def addNewTransient_format1_sort0(self):
        """initializes the transient variables"""
        self.elementName = self.dataCode['elementName']
        if self.dt not in self.stress:
            self.stress[self.dt] = {}

    def addNewTransient_format2_sort1(self):
        """initializes the transient variables"""
        raise Exception('not implemented')
        #print self.dataCode
        self.elementName = self.dataCode['elementName']
        if self.dt not in self.stress:
            self.axial[self.dt]     = {}
            self.torsion[self.dt]   = {}

    def addNewEid_format1_sort0(self,out):
        #print "Rod Stress add..."
        (eid,stress) = out
        eid = (eid-self.deviceCode) // 10
        #assert isinstance(eid,int)
        self.eType[eid]  = self.elementName
        self.stress[eid] = stress

    def addNewEid_format2_sort1(self,out):
        raise Exception('not implemented')
        (eid,axialReal,axialImag,torsionReal,torsionImag) = out
        eid = (eid-self.deviceCode) // 10
        assert eid >= 0
        self.eType[eid]      = self.elementName
        self.axial[eid]      = [axialReal,axialImag]
        self.torsion[eid]    = [torsionReal,torsionImag]

    def addNewEidTransient_format1_sort0(self,out):
        (eid,stress) = out
        eid = (eid-self.deviceCode) // 10
        self.eType[eid] = self.elementName
        dt = self.dt
        #assert isinstance(eid,int)
        #assert eid >= 0
        self.eType[eid]           = self.elementName
        self.stress[dt][eid]      = stress

    def addNewEidTransient_format2_sort1(self,out):
        raise Exception('not implemented')
        (eid,axialReal,axialImag,torsionReal,torsionImag) = out
        eid = (eid-self.deviceCode) // 10
        dt = self.dt
        assert eid >= 0
        self.eType[eid]          = self.elementName
        self.axial[dt][eid]      = [axialReal,axialImag]
        self.torsion[dt][eid]    = [torsionReal,torsionImag]

    def __reprTransient_format1_sort0__(self):
        msg = '---CELASx STRESSES---\n'
        msg += '%-6s %6s ' %('EID','eType')
        headers = ['stress']
        for header in headers:
            msg += '%10s ' %(header)
        msg += '\n'

        for dt,stress in sorted(self.stress.items()):
            msg += '%s = %g\n' %(self.dataCode['name'],dt)
            for eid,istress in sorted(stress.items()):
                msg += '%-6g %6s ' %(eid,self.eType[eid])
                if abs(istress)<1e-6:
                    msg += '%10s ' %('0')
                else:
                    msg += '%10g ' %(istress)
                ###
                msg += '\n'
            ###
        return msg

    def __reprTransient_format2_sort1__(self):
        raise Exception('not implemented')
        msg = '---COMPLEX CELASx STRESSES---\n'
        msg += '%-10s %10s ' %('EID','eType')
        headers = ['axialReal','axialImag','torsionReal','torsionImag']
        for header in headers:
            msg += '%10s ' %(header)
        msg += '\n'

        for dt,axial in sorted(self.axial.items()):
            msg += '%s = %g\n' %(self.dataCode['name'],dt)
            for eid in sorted(axial):
                axial   = self.axial[dt][eid]
                torsion = self.torsion[dt][eid]
                msg += '%-6g %6s ' %(eid,self.eType)
                vals = axial + torsion # concatination
                for val in vals:
                    if abs(val)<1e-6:
                        msg += '%10s ' %('0')
                    else:
                        msg += '%10i ' %(val)
                    ###
                msg += '\n'
                #msg += "eid=%-4s eType=%s axial=%-4i torsion=%-4i\n" %(eid,self.eType,axial,torsion)
            ###
        return msg

    def __repr__(self):
        if   self.isTransient and self.code==[1,0,0]:
            return self.__reprTransient_format1_sort0__()
        elif self.code==[2,1,0]:
            return self.__reprTransient_format2_sort1__()

        msg = '---CELASx STRESSES---\n'
        msg += '%-8s %6s ' %('EID','eType')
        headers = ['stress']
        for header in headers:
            msg += '%10s ' %(header)
        msg += '\n'
        #print "self.code = ",self.code
        for eid,istress in sorted(self.stress.items()):
            #print "eid=",eid
            #print "eType",self.eType
            msg += '%-8i %6s ' %(eid,self.eType[eid])
            if abs(istress)<1e-6:
                msg += '%10s ' %('0')
            else:
                msg += '%10i ' %(istress)
            ###
            msg += '\n'
            #msg += "eid=%-4s eType=%s axial=%-4i torsion=%-4i\n" %(eid,self.eType,axial,torsion)
        return msg

class celasStrainObject(strainObject):
    """
    """
    def __init__(self,dataCode,iSubcase,dt=None):
        strainObject.__init__(self,dataCode,iSubcase)
        self.eType = {}
        self.elementName = self.dataCode['elementName']

        self.code = [self.formatCode,self.sortCode,self.sCode]
        
        self.isTransient = False
        if self.code == [1,0,10]:
            self.getLength = self.getLength_format1_sort0
            self.strain = {}
            #self.isImaginary = False
            if dt is not None:
                self.addNewTransient = self.addNewTransient_format1_sort0
                self.addNewEid       = self.addNewEidTransient_format1_sort0
            ###
            else:
                self.addNewEid = self.addNewEid_format1_sort0
            ###

        #elif self.code in [[2,1,10]]:
        #    self.getLength = self.getLength_format1_sort0
        #    self.axial      = {}
        #    self.torsion    = {}
        #    self.addNewTransient = self.addNewTransient_format2_sort1
        #    self.addNewEid       = self.addNewEidTransient_format2_sort1
        #    #self.isImaginary = True
        else:
            raise InvalidCodeError('celasStrain - get the format/sort/stressCode=%s' %(self.code))
        ###
        if dt is not None:
            self.dt = self.nonlinearFactor
            self.isTransient = True
            self.addNewTransient()
        ###

    def getLength_format1_sort0(self):
        return (8,'if')

    def deleteTransient(self,dt):
        del self.strain[dt]

    def getTransients(self):
        k = self.strain.keys()
        k.sort()
        return k

    def addNewTransient_format1_sort0(self):
        """
        initializes the transient variables
        @note make sure you set self.dt first
        """
        if self.dt not in self.strain:
            self.strain[self.dt] = {}

    def addNewTransient_format2_sort1(self):
        """
        initializes the transient variables
        """
        raise Exception('not implemented')
        #print self.dataCode
        if self.dt not in self.strain:
            self.axial[self.dt]     = {}
            self.torsion[self.dt]   = {}

    def addNewEid_format1_sort0(self,out):
        (eid,strain) = out
        eid = (eid-self.deviceCode) // 10
        #print "Rod Strain add..."
        assert eid >= 0
        #self.eType = self.eType
        self.eType[eid]  = self.elementName
        self.strain[eid] = strain

    def addNewEid_format2_sort1(self,out):
        raise Exception('not implemented')
        assert eid >= 0
        (eid,axialReal,axialImag,torsionReal,torsionImag) = out
        eid = (eid-self.deviceCode) // 10
        self.eType[eid] = self.elementType
        self.axial[eid]      = [axialReal,axialImag]
        self.torsion[eid]    = [torsionReal,torsionImag]

    def addNewEidTransient_format1_sort0(self,out):
        #print "Rod Strain add..."
        #print out
        (eid,strain) = out
        eid = (eid-self.deviceCode) // 10
        assert eid >= 0
        dt = self.dt

        self.eType[eid] = self.elementType
        self.strain[dt][eid] = strain

    def addNewEidTransient_format2_sort1(self,out):
        raise Exception('not implemented')
        (eid,axialReal,axialImag,torsionReal,torsionImag) = out
        eid = (eid-self.deviceCode) // 10
        assert eid >= 0
        dt = self.dt
        self.axial[dt][eid]      = [axialReal,axialImag]
        self.torsion[dt][eid]    = [torsionReal,torsionImag]

    def __reprTransient_format2_sort1__(self):
        raise Exception('not implemented')
        msg = '---COMPLEX CELASx STRAINS---\n'
        msg += '%-8s %10s ' %('EID','eType')
        headers = ['axialReal','axialImag','torsionReal','torsionImag']
        for header in headers:
            msg += '%10s ' %(header)
        msg += '\n'

        for dt,axial in sorted(self.axial.items()):
            msg += '%s = %g\n' %(self.dataCode['name'],dt)
            for eid in sorted(axial):
                axial   = self.axial[dt][eid]
                torsion = self.torsion[dt][eid]
                msg += '%-8i %6s ' %(eid,self.eType[eid])
                vals = axial + torsion # concatination
                for val in vals:
                    if abs(val)<1e-6:
                        msg += '%10s ' %('0')
                    else:
                        msg += '%10g ' %(val)
                    ###
                msg += '\n'
                #msg += "eid=%-4s eType=%s axial=%-4i torsion=%-4i\n" %(eid,self.eType,axial,torsion)
            ###
        return msg

    def __repr__(self):
        if   self.isTransient and self.code==[1,0,10]:
            return self.__reprTransient_format1_sort0__()
        elif self.code==[2,1,10]:
            return self.__reprTransient_format2_sort1__()

        msg = '---CELASx STRAINS---\n'
        msg += '%-8s %6s ' %('EID','eType')
        headers = ['strain']
        for header in headers:
            msg += '%8s ' %(header)
        msg += '\n'

        for eid,strain in sorted(self.strain.items()):
            #strain = self.strain[eid]
            msg += '%-8i %6s ' %(eid,self.eType[eid])

            if abs(strain)<1e-7:
                msg += '%8s ' %('0')
            else:
                msg += '%8.3g ' %(strain)
            ###
            msg += '\n'
        return msg
