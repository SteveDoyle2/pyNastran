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
import copy

# pyNastran
from pyNastran.op2.resultObjects.tableObject import TableObject,complexTableObject

class displacementObject(TableObject): # approachCode=1, sortCode=0, thermal=0
    def __init__(self,dataCode,iSubcase,dt=None):
        TableObject.__init__(self,dataCode,iSubcase,dt)

    def writeF06(self,header,pageStamp,pageNum=1):
        if self.dt is not None:
            return self.writeF06Transient(header,pageStamp,pageNum)
        msg = header+['                                             D I S P L A C E M E N T   V E C T O R\n',
               ' \n',
               '      POINT ID.   TYPE          T1             T2             T3             R1             R2             R3\n']
        for nodeID,translation in sorted(self.translations.items()):
            rotation = self.rotations[nodeID]
            gridType = self.gridTypes[nodeID]

            (dx,dy,dz) = translation
            (rx,ry,rz) = rotation
            vals = [dx,dy,dz,rx,ry,rz]
            (vals2,isAllZeros) = self.writeF06Floats13E(vals)
            [dx,dy,dz,rx,ry,rz] = vals2
            msg.append('%14i %6s     %13s  %13s  %13s  %13s  %13s  %-s\n' %(nodeID,gridType,dx,dy,dz,rx,ry,rz.rstrip()))
        ###
        msg.append(pageStamp+str(pageNum)+'\n')
        return (''.join(msg),pageNum)

    def writeF06Transient(self,header,pageStamp,pageNum=1):
        words = ['                                             D I S P L A C E M E N T   V E C T O R\n',
                 ' \n',
                 '      POINT ID.   TYPE          T1             T2             T3             R1             R2             R3\n']
        msg = []
        for dt,translations in sorted(self.translations.items()):
            header[1] = ' %s = %10.4E\n' %(self.dataCode['name'],dt)
            msg += header+words
            for nodeID,translation in sorted(translations.items()):
                rotation = self.rotations[dt][nodeID]
                gridType = self.gridTypes[nodeID]

                (dx,dy,dz) = translation
                (rx,ry,rz) = rotation
                vals = [dx,dy,dz,rx,ry,rz]
                (vals2,isAllZeros) = self.writeF06Floats13E(vals)
                [dx,dy,dz,rx,ry,rz] = vals2
                msg.append('%14i %6s     %13s  %13s  %13s  %13s  %13s  %-s\n' %(nodeID,gridType,dx,dy,dz,rx,ry,rz.rstrip()))
            ###
            msg.append(pageStamp+str(pageNum)+'\n')
        return (''.join(msg),pageNum)

    def __repr__(self):
        if self.dt is not None:
            return self.__reprTransient__()

        msg = '---DISPLACEMENTS---\n'
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
        msg = '---TRANSIENT DISPLACEMENTS---\n'
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

class complexDisplacementObject(complexTableObject): # approachCode=1, sortCode=0, thermal=0
    def __init__(self,dataCode,iSubcase,dt=None):
        complexTableObject.__init__(self,dataCode,iSubcase,dt)

    def writeF06(self,header,pageStamp,pageNum=1):
        words = ['                                       C O M P L E X   D I S P L A C E M E N T   V E C T O R\n',
                 '                                                          (REAL/IMAGINARY)\n',
                 ' \n',
                 '      POINT ID.   TYPE          T1             T2             T3             R1             R2             R3\n']
        msg = []
        for dt,translations in sorted(self.translations.items()):
            header[2] = ' %s = %10.4E\n' %(self.dataCode['name'],dt)
            msg += header+words
            for nodeID,translation in sorted(translations.items()):
                rotation = self.rotations[dt][nodeID]
                gridType = self.gridTypes[nodeID]

                (dxr,dxi,dyr,dyi,dzr,dzi) = translation
                (rxr,rxi,ryr,ryi,rzr,rzi) = rotation
                vals = [dxr,dyr,dzr,rxr,ryr,rzr,dxi,dyi,dzi,rxi,ryi,rzi]
                (vals2,isAllZeros) = self.writeF06Floats13E(vals)
                [dxr,dyr,dzr,rxr,ryr,rzr,dxi,dyi,dzi,rxi,ryi,rzi] = vals2
                msg.append('0 %12i %6s     %13s  %13s  %13s  %13s  %13s  %-s\n' %(nodeID,gridType,dxr,dyr,dzr,rxr,ryr,rzr.rstrip()))
                msg.append('  %12s %6s     %13s  %13s  %13s  %13s  %13s  %-s\n' %('','',          dxi,dyi,dzi,rxi,ryi,rzi.rstrip()))
            ###
            msg.append(pageStamp+str(pageNum)+'\n')
        return (''.join(msg),pageNum-1)

    def __repr__(self):
        return self.writeF06(['','',''],'PAGE ',1)[0]

        msg = '---COMPLEX DISPLACEMENTS---\n'
        #if self.dt is not None:
        #    msg += '%s = %g\n' %(self.dataCode['name'],self.dt)
        headers = ['DxReal','DxImag','DyReal','DyImag','DzReal','DyImag','RxReal','RxImag','RyReal','RyImag','RzReal','RzImag']
        msg += '%-10s ' %('nodeID')
        for header in headers:
            msg += '%10s ' %(header)
        msg += '\n'

        for freq,translations in sorted(self.translations.items()):
            msg += '%s = %g\n' %(self.dataCode['name'],dt)

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
