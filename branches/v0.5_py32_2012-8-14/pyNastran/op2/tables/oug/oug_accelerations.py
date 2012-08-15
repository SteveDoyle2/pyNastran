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
#import copy

# pyNastran
from pyNastran.op2.resultObjects.tableObject import TableObject,ComplexTableObject

class AccelerationObject(TableObject): # approachCode=11, thermal=0
    def __init__(self,dataCode,isSort1,iSubcase,dt=None):
        TableObject.__init__(self,dataCode,isSort1,iSubcase,dt)

    def writeMatlab(self,iSubcase,f=None,isMagPhase=False):
        name = 'accelerations'
        if self.nonlinearFactor is None:
            return self._writeMatlab(name,iSubcase,f)
        else:
            return self._writeMatlabTransient(name,iSubcase,f)

    def writeF06(self,header,pageStamp,pageNum=1,f=None,isMagPhase=False):
        if self.nonlinearFactor is not None:
            return self.writeF06Transient(header,pageStamp,pageNum,f)
        words = ['                                             A C C E L E R A T I O N   V E C T O R\n',
               ' \n',
               '      POINT ID.   TYPE          T1             T2             T3             R1             R2             R3\n']
        words += self.getTableMarker()
        return self._writeF06Block(words,header,pageStamp,pageNum,f)

    def writeF06Transient(self,header,pageStamp,pageNum=1,f=None,isMagPhase=False):
        words = ['                                             A C C E L E R A T I O N   V E C T O R\n',
                 ' \n',
                 '      POINT ID.   TYPE          T1             T2             T3             R1             R2             R3\n']
        words += self.getTableMarker()
        return self._writeF06TransientBlock(words,header,pageStamp,pageNum,f)

    def __repr__(self):
        #return ''
        if self.nonlinearFactor is not None:
            return self.__reprTransient__()

        msg = '---ACCELERATIONS---\n'
        msg += self.writeHeader()

        for nodeID,translation in sorted(self.translations.items()):
            rotation = self.rotations[nodeID]
            gridType = self.gridTypes[nodeID]

            (dx,dy,dz) = translation
            (rx,ry,rz) = rotation

            msg += '%-10i %-8s ' %(nodeID,gridType)
            vals = [dx,dy,dz,rx,ry,rz]
            for val in vals:
                if abs(val)<1e-6:
                    msg += '%10s ' %(0)
                else:
                    msg += '%10.3e ' %(val)
                ###
            msg += '\n'
        return msg

    def __reprTransient__(self):
        msg = '---TRANSIENT ACCELERATIONS---\n'
        msg += self.writeHeader()
        
        for dt,translations in sorted(self.translations.items()):
            msg += '%s = %g\n' %(self.dataCode['name'],dt)
            for nodeID,translation in sorted(translations.items()):
                rotation = self.rotations[dt][nodeID]
                gridType = self.gridTypes[nodeID]
                (dx,dy,dz) = translation
                (rx,ry,rz) = rotation

                msg += '%-10i %8s ' %(nodeID,gridType)
                vals = [dx,dy,dz,rx,ry,rz]
                for val in vals:
                    if abs(val)<1e-6:
                        msg += '%10s ' %(0)
                    else:
                        msg += '%10.3e ' %(val)
                    ###
                msg += '\n'
            ###
        return msg

class ComplexAccelerationObject(ComplexTableObject): # tableCode=11, approachCode=???
    def __init__(self,dataCode,isSort1,iSubcase,dt=None):
        ComplexTableObject.__init__(self,dataCode,isSort1,iSubcase,dt)

    def writeMatlab(self,iSubcase,f=None,isMagPhase=False):
        name = 'accelerations'
        if self.nonlinearFactor is None:
            return self._writeMatlab(name,iSubcase,f,isMagPhase)
        else:
            return self._writeMatlabTransient(name,iSubcase,f,isMagPhase)

    def writeF06(self,header,pageStamp,pageNum=1,f=None,isMagPhase=False):
        if self.nonlinearFactor is not None:
            return self.writeF06Transient(header,pageStamp,pageNum,f,isMagPhase)

        words = ['                                       C O M P L E X   A C C E L E R A T I O N   V E C T O R\n']
        return self._writeF06Block(words,header,pageStamp,pageNum,f,isMagPhase)

    def writeF06Transient(self,header,pageStamp,pageNum=1,f=None,isMagPhase=False):
        words = ['                                       C O M P L E X   A C C E L E R A T I O N   V E C T O R\n']
        return self._writeF06TransientBlock(words,header,pageStamp,pageNum,f,isMagPhase)

    def __repr__(self):
        return self.writeF06(['','',''],'PAGE ',1)[0]
        #if self.nonlinearFactor is not None:
            #return self.__reprTransient__()

        msg = '---COMPLEX ACCELERATIONS---\n'
        #if self.nonlinearFactor is not None:
        #    msg += '%s = %g\n' %(self.dataCode['name'],self.dt)
        headers = ['DxReal','DxImag','DyReal','DyImag','DzReal','DyImag','RxReal','RxImag','RyReal','RyImag','RzReal','RzImag']
        msg += '%-10s ' %('nodeID')
        for header in headers:
            msg += '%10s ' %(header)
        msg += '\n'

        for freq,translations in sorted(self.translations.items()):
            msg += '%s = %g\n' %(self.dataCode['name'],freq)

            for nodeID,translation in sorted(translations.items()):
                rotation = self.rotations[freq][nodeID]

                msg += '%-10i ' %(nodeID)
                vals = translation+rotation
                for val in vals:
                    if abs(val)<1e-6:
                        msg += '%10s ' %(0)
                    else:
                        msg += '%10.3e ' %(val)
                    ###
                msg += '\n'
            ###
        return msg
