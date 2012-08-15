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
#from numpy import array
from pyNastran.op2.resultObjects.tableObject import TableObject,ComplexTableObject

class SPCForcesObject(TableObject):
    def __init__(self,dataCode,isSort1,iSubcase,dt=None):
        TableObject.__init__(self,dataCode,isSort1,iSubcase,dt)

    def writeMatlab(self,iSubcase,f=None,isMagPhase=False):
        name = 'spcForces'
        if self.nonlinearFactor is None:
            return self._writeMatlab(name,iSubcase,f)
        else:
            return self._writeMatlabTransient(name,iSubcase,f)

    def writeF06(self,header,pageStamp,pageNum=1,f=None,isMagPhase=False):
        if self.nonlinearFactor is not None:
            return self.writeF06Transient(header,pageStamp,pageNum,f)
        msg = header+['                               F O R C E S   O F   S I N G L E - P O I N T   C O N S T R A I N T\n',
               ' \n',
               '      POINT ID.   TYPE          T1             T2             T3             R1             R2             R3\n']
        for nodeID,translation in sorted(self.translations.items()):
            rotation = self.rotations[nodeID]
            gridType = self.gridTypes[nodeID]

            (dx,dy,dz) = translation
            (rx,ry,rz) = rotation
            vals = [dx,dy,dz,rx,ry,rz]
            (vals2,isAllZeros) = self.writeFloats13E(vals)
            if not isAllZeros:
                [dx,dy,dz,rx,ry,rz] = vals2
                msg.append('%14i %6s     %13s  %13s  %13s  %13s  %13s  %-s\n' %(nodeID,gridType,dx,dy,dz,rx,ry,rz.rstrip()))
        ###
        msg.append(pageStamp+str(pageNum)+'\n')
        if f is not None:
            f.write(''.join(msg))
            msg = ['']
        return (''.join(msg),pageNum)

    def writeF06Transient(self,header,pageStamp,pageNum=1,f=None,isMagPhase=False):
        words = ['                               F O R C E S   O F   S I N G L E - P O I N T   C O N S T R A I N T\n',
                 ' \n',
                 '      POINT ID.   TYPE          T1             T2             T3             R1             R2             R3\n']
        return self._writeF06TransientBlock(words,header,pageStamp,pageNum,f)

    def __reprTransient__(self):
        msg = '---SPC FORCES---\n'
        if self.nonlinearFactor is not None:
            msg += 'dt = %g\n' %(self.dt)

        headers = ['T1','T2','T3','R1','R2','R3']
        msg += '%-8s ' %('GRID')
        for header in headers:
            msg += '%10s ' %(header)
        msg += '\n'

        for dt,translations in sorted(self.translations.items()):
            msg += 'dt = %s' %(dt)
            for nodeID,translation in sorted(translations.items()):
                rotation = self.rotations[dt][nodeID]
                (Fx,Fy,Fz) = translation
                (Mx,My,Mz) = rotation

                msg += '%-8i ' %(nodeID)
                vals = [Fx,Fy,Fz,Mx,My,Mx]
                for val in vals:
                    if abs(val)<1e-6:
                        msg += '%10s ' %(0)
                    else:
                        msg += '%10.2f ' %(val)
                    ###
                msg += '\n'
            ###
        return msg

    def __repr__(self):
        return self.writeF06(['',''],'PAGE ',1)[0]

        if self.nonlinearFactor is not None:
            return self.__reprTransient__()

        msg = '---SPC FORCES---\n'
        if self.dt is not None:
            msg += 'dt = %g\n' %(self.dt)

        headers = ['T1','T2','T3','R1','R2','R3']
        msg += '%-8s ' %('GRID')
        for header in headers:
            msg += '%10s ' %(header)
        msg += '\n'

        for nodeID,translation in sorted(self.translations.items()):
            rotation = self.rotations[nodeID]
            (Fx,Fy,Fz) = translation
            (Mx,My,Mz) = rotation

            msg += '%-8i ' %(nodeID)
            vals = [Fx,Fy,Fz,Mx,My,Mx]
            for val in vals:
                if abs(val)<1e-6:
                    msg += '%10s ' %(0)
                else:
                    msg += '%10.2f ' %(val)
                ###
            msg += '\n'
        return msg

class ComplexSPCForcesObject(ComplexTableObject):
    def __init__(self,dataCode,isSort1,iSubcase,dt=None):
        ComplexTableObject.__init__(self,dataCode,isSort1,iSubcase,dt)

    def writeMatlab(self,iSubcase,f=None,isMagPhase=False):
        name = 'spcForces'
        if self.nonlinearFactor is None:
            return self._writeMatlab(name,iSubcase,f)
        else:
            return self._writeMatlabTransient(name,iSubcase,f)

    def writeF06(self,header,pageStamp,pageNum=1,f=None,isMagPhase=False):
        if self.nonlinearFactor is not None:
            return self.writeF06Transient(header,pageStamp,pageNum,f,isMagPhase)
        msg = header+['                               F O R C E S   O F   S I N G L E - P O I N T   C O N S T R A I N T\n',
               ' \n',
               '      POINT ID.   TYPE          T1             T2             T3             R1             R2             R3\n']
        raise RuntimeError('is this valid...')
        for nodeID,translation in sorted(self.translations.items()):
            rotation = self.rotations[nodeID]
            gridType = self.gridTypes[nodeID]

            (dx,dy,dz) = translation
            #dxr=dx.real; dyr=dy.real; dzr=dz.real; 
            #dxi=dx.imag; dyi=dy.imag; dzi=dz.imag

            (rx,ry,rz) = rotation
            #rxr=rx.real; ryr=ry.real; rzr=rz.real
            #rxi=rx.imag; ryi=ry.imag; rzi=rz.imag

            #vals = [dxr,dyr,dzr,rxr,ryr,rzr,dxi,dyi,dzi,rxi,ryi,rzi]
            vals = list(translation)+list(rotation)
            (vals2,isAllZeros) = self.writeFloats13E(vals)
            if not isAllZeros:
                [dx,dy,dz,rx,ry,rz] = vals2
                msg.append('%14i %6s     %13s  %13s  %13s  %13s  %13s  %-s\n' %(nodeID,gridType,dx,dy,dz,rx,ry,rz.rstrip()))
        ###
        msg.append(pageStamp+str(pageNum)+'\n')
        if f is not None:
            f.write(''.join(msg))
            msg = ['']
        return (''.join(msg),pageNum)

    def writeF06Transient(self,header,pageStamp,pageNum=1,f=None,isMagPhase=False):
        words = ['                         C O M P L E X   F O R C E S   O F   S I N G L E   P O I N T   C O N S T R A I N T\n']
        return self._writeF06TransientBlock(words,header,pageStamp,pageNum,f,isMagPhase)

    def __reprTransient__(self):
        msg = '---COMPLEX SPC FORCES---\n'
        if self.nonlinearFactor is not None:
            msg += 'dt = %g\n' %(self.dt)

        raise RuntimeError('is this valid...')
        headers = ['T1','T2','T3','R1','R2','R3']
        msg += '%-8s ' %('GRID')
        for header in headers:
            msg += '%10s ' %(header)
        msg += '\n'

        for dt,translations in sorted(self.translations.items()):
            msg += 'dt = %s' %(dt)
            for nodeID,translation in sorted(translations.items()):
                rotation = self.rotations[dt][nodeID]
                msg += '%-8i ' %(nodeID)
                vals = translation+rotation
                for val in vals:
                    if abs(val)<1e-6:
                        msg += '%10s ' %(0)
                    else:
                        msg += '%10.2f ' %(val)
                    ###
                msg += '\n'
            ###
        return msg

    def __repr__(self):
        return self.writeF06(['',''],'PAGE ',1)[0]
        if self.nonlinearFactor is not None:
            return self.__reprTransient__()

        msg = '---COMPLEX SPC FORCES---\n'
        if self.nonlinearFactor is not None:
            msg += 'dt = %g\n' %(self.dt)

        raise RuntimeError('is this valid...')
        headers = ['T1','T2','T3','R1','R2','R3']
        msg += '%-8s ' %('GRID')
        for header in headers:
            msg += '%10s ' %(header)
        msg += '\n'

        for nodeID,translation in sorted(self.translations.items()):
            rotation = self.rotations[nodeID]
            (Fx,Fy,Fz) = translation
            (Mx,My,Mz) = rotation

            msg += '%-8i ' %(nodeID)
            vals = [Fx,Fy,Fz,Mx,My,Mx]
            for val in vals:
                if abs(val)<1e-6:
                    msg += '%10s ' %(0)
                else:
                    msg += '%10.2f ' %(val)
                ###
            msg += '\n'
        return msg

