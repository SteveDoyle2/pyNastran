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

class MPCForcesObject(TableObject):
    def __init__(self,dataCode,isSort1,iSubcase,dt=None):
        TableObject.__init__(self,dataCode,isSort1,iSubcase,dt)

    def writeMatlab(self,iSubcase,f=None,isMagPhase=False):
        name = 'mpcForces'
        if self.nonlinearFactor is None:
            return self._writeMatlab(name,iSubcase,f)
        else:
            return self._writeMatlabTransient(name,iSubcase,f)

    def writeF06(self,header,pageStamp,pageNum=1,f=None,isMagPhase=False):
        if self.nonlinearFactor is not None:
            return self.writeF06Transient(header,pageStamp,pageNum)
        msg = header+['                               F O R C E S   O F   M U L T I - P O I N T   C O N S T R A I N T\n',
               ' \n',
               '      POINT ID.   TYPE          T1             T2             T3             R1             R2             R3\n']
        for nodeID,translation in sorted(self.translations.iteritems()):
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
        ###
        msg.append(pageStamp+str(pageNum)+'\n')
        if f is not None:
            f.write(''.join(msg))
            msg = ['']
        return (''.join(msg),pageNum)

    def writeF06Transient(self,header,pageStamp,pageNum=1,f=None,isMagPhase=False):
        words = ['                               F O R C E S   O F   M U L T I - P O I N T   C O N S T R A I N T\n',
                 ' \n',
                 '      POINT ID.   TYPE          T1             T2             T3             R1             R2             R3\n']
        msg = []
        for dt,translations in sorted(self.translations.iteritems()):
            header[1] = ' %s = %10.4E\n' %(self.dataCode['name'],dt)
            msg += header+words
            for nodeID,translation in sorted(translations.iteritems()):
                rotation = self.rotations[dt][nodeID]
                gridType = self.gridTypes[nodeID]

                (dx,dy,dz) = translation
                (rx,ry,rz) = rotation
                vals = [dx,dy,dz,rx,ry,rz]
                (vals2,isAllZeros) = self.writeFloats13E(vals)
                if not isAllZeros:
                    [dx,dy,dz,rx,ry,rz] = vals2
                    msg.append('%14i %6s     %13s  %13s  %13s  %13s  %13s  %-s\n' %(nodeID,gridType,dx,dy,dz,rx,ry,rz.rstrip()))
                ###
            ###
            msg.append(pageStamp+str(pageNum)+'\n')
            if f is not None:
                f.write(''.join(msg))
                msg = ['']
            pageNum+=1
        return (''.join(msg),pageNum-1)

    def __reprTransient__(self):
        msg = '---MPC FORCES---\n'
        if self.nonlinearFactor is not None:
            msg += 'dt = %g\n' %(self.dt)

        headers = ['T1','T2','T3','R1','R2','R3']
        msg += '%-8s ' %('GRID')
        for header in headers:
            msg += '%10s ' %(header)
        msg += '\n'

        for dt,translations in sorted(self.translations.iteritems()):
            msg += 'dt = %s' %(dt)
            for nodeID,translation in sorted(translations.iteritems()):
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
        if self.nonlinearFactor is not None:
            return self.__reprTransient__()

        msg = '---MPC FORCES---\n'
        if self.dt is not None:
            msg += 'dt = %g\n' %(self.dt)

        headers = ['T1','T2','T3','R1','R2','R3']
        msg += '%-8s ' %('GRID')
        for header in headers:
            msg += '%10s ' %(header)
        msg += '\n'

        for nodeID,translation in sorted(self.translations.iteritems()):
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

class ComplexMPCForcesObject(ComplexTableObject):
    def __init__(self,dataCode,isSort1,iSubcase,dt=None):
        ComplexTableObject.__init__(self,dataCode,isSort1,iSubcase,dt)

    def writeMatlab(self,iSubcase,f=None,isMagPhase=False):
        name = 'mpcForces'
        if self.nonlinearFactor is None:
            return self._writeMatlab(name,iSubcase,f,isMagPhase)
        else:
            return self._writeMatlabTransient(name,iSubcase,f,isMagPhase)

    def writeF06(self,header,pageStamp,pageNum=1,f=None,isMagPhase=False):
        if self.nonlinearFactor is not None:
            return self.writeF06Transient(header,pageStamp,pageNum,f,isMagPhase)
        msg = header+['                               F O R C E S   O F   M U L T I - P O I N T   C O N S T R A I N T\n',
               ' \n',
               '      POINT ID.   TYPE          T1             T2             T3             R1             R2             R3\n']
        raise RuntimeError('is this valid...')
        for nodeID,translation in sorted(self.translations.iteritems()):
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
        words = ['                         C O M P L E X   F O R C E S   O F   M U L T I   P O I N T   C O N S T R A I N T\n',
                 '                                                          (REAL/IMAGINARY)\n',
                 ' \n',
                 '      POINT ID.   TYPE          T1             T2             T3             R1             R2             R3\n']
        msg = []
        for dt,translations in sorted(self.translations.iteritems()):
            header[1] = ' %s = %10.4E\n' %(self.dataCode['name'],dt)
            msg+= header+words
            for nodeID,translation in sorted(translations.iteritems()):
                rotation = self.rotations[dt][nodeID]
                gridType = self.gridTypes[nodeID]

                (dx,dy,dz) = translation
                (rx,ry,rz) = rotation

                vals = [dx,dy,dz,rx,ry,rz]
                (vals2,isAllZeros) = self.writeImagFloats13E(vals,isMagPhase)
                if not isAllZeros:
                    [v1r,v2r,v3r,v4r,v5r,v6r,  v1i,v2i,v3i,v4i,v5i,v6i] = vals2
                    msg.append('0%13i %6s     %13s  %13s  %13s  %13s  %13s  %-s\n' %(nodeID,gridType,v1r,v2r,v3r,v4r,v5r,v6r.rstrip()))
                    msg.append(' %13i %6s     %13s  %13s  %13s  %13s  %13s  %-s\n' %(nodeID,gridType,v1i,v2i,v3i,v4i,v5i,v6i.rstrip()))
                ###
            ###
            msg.append(pageStamp+str(pageNum)+'\n')
            if f is not None:
                f.write(''.join(msg))
                msg = ['']
            pageNum+=1
        return (''.join(msg),pageNum-1)

    def __reprTransient__(self):
        msg = '---COMPLEX MPC FORCES---\n'
        if self.nonlinearFactor is not None:
            msg += 'dt = %g\n' %(self.dt)

        raise RuntimeError('is this valid...')
        headers = ['T1','T2','T3','R1','R2','R3']
        msg += '%-8s ' %('GRID')
        for header in headers:
            msg += '%10s ' %(header)
        msg += '\n'

        for dt,translations in sorted(self.translations.iteritems()):
            msg += 'dt = %s' %(dt)
            for nodeID,translation in sorted(translations.iteritems()):
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
        if self.nonlinearFactor is not None:
            return self.__reprTransient__()

        msg = '---COMPLEX MPC FORCES---\n'
        if self.dt is not None:
            msg += 'dt = %g\n' %(self.dt)

        raise RuntimeError('is this valid...')
        headers = ['T1','T2','T3','R1','R2','R3']
        msg += '%-8s ' %('GRID')
        for header in headers:
            msg += '%10s ' %(header)
        msg += '\n'

        for nodeID,translation in sorted(self.translations.iteritems()):
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

