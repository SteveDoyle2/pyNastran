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
from struct import pack
from pyNastran.op2.resultObjects.op2_Objects import scalarObject,array

class gridPointForcesObject(scalarObject):
    def __init__(self,dataCode,iSubcase,dt=None):
        scalarObject.__init__(self,dataCode,iSubcase)
        self.forces = {}
        self.moments = {}
        self.elemName = {}
        self.eids = {}
        if dt is not None:
            self.dt = dt
            self.isTransient = True
            self.addNewTransient()
            self.add       = self.addTransient
            #self.addNewEid = self.addNewEidTransient
        ###

    def addNewTransient(self,eKey):
        """initializes the transient variables"""
        if self.dt not in self.forces:
            self.forces[self.dt] = {}
            self.moments[self.dt] = {}
            self.elemName[eKey] = []
            self.eids[eKey] = []

    def add(self,eKey,eid,elemName,f1,f2,f3,m1,m2,m3):
        if eKey not in self.forces:
            self.eids[eKey] = []
            self.forces[eKey] = []
            self.moments[eKey] = []
            self.elemName[eKey] = []
        self.forces[eKey].append( [f1,f2,f3]) # Fx,Fy,Fz
        self.moments[eKey].append([m1,m2,m3]) # Mx,My,Mz
        self.elemName[eKey].append(elemName)
        self.eids[eKey].append(eid)

    def addTransient(self,eKey,eid,elemName,f1,f2,f3,m1,m2,m3):
        if eKey not in self.forces[self.dt]:
            self.eids[self.dt][eKey] = []
            self.forces[self.dt][eKey] = []
            self.moments[self.dt][eKey] = []
            self.elemName[self.dt][eKey] = []
        self.forces[self.dt][eKey].append( [f1,f2,f3]) # Mx,My,Mz
        self.moments[self.dt][eKey].append([m1,m2,m3]) # Fx,Fy,Fz
        self.elemName[self.dt][eKey].append(elemName)
        self.eids[self.dt][eKey].append(eid)

    def updateDt(self,dataCode,freq):
        self.dataCode = dataCode
        self.applyDataCode()
        if freq is not None:
            self.log.debug("updating %s...%s=%s  iSubcase=%s" %(self.dataCode['name'],self.dataCode['name'],freq,self.iSubcase))
            self.dt = dt
            self.addNewTransient()
        ###

    def deleteTransient(self,dt):
        del self.forces[dt]
        del self.moments[dt]
        del self.elemName[dt]

    def getTransients(self):
        k = self.forces.keys()
        k.sort()
        return k
    
    #def cleanupObj(self):
        #k = self.elemName.keys()
        #self.elemName = self.elemName[k[0]]
        #self.eids = self.eids[k[0]]

    def writeF06(self,header,pageStamp,pageNum=1):
        if self.isTransient:
            return self.writeF06Transient(header,pageStamp,pageNum)

        msg = header+['                                          G R I D   P O I N T   F O R C E   B A L A N C E\n',
               ' \n',
               '   POINT-ID    ELEMENT-ID     SOURCE             T1             T2             T3             R1             R2             R3\n',]
              #'0     13683          3736    TRIAX6         4.996584E+00   0.0            1.203093E+02   0.0            0.0            0.0'
              #'      13683          3737    TRIAX6        -4.996584E+00   0.0           -1.203093E+02   0.0            0.0            0.0'
              #'      13683                  *TOTALS*       6.366463E-12   0.0           -1.364242E-12   0.0            0.0            0.0'
        for eKey,force in sorted(self.forces.items()):
            zero = '0'
            for iLoad,f in enumerate(force):
                (f1,f2,f3) = f
                (m1,m2,m3) = self.moments[eKey][iLoad]
                (elemName) = self.elemName[eKey][iLoad]
                eid = self.eids[eKey][iLoad]
                vals = [f1,f2,f3,m1,m2,m3]
                (vals2,isAllZeros) = self.writeF06Floats13E(vals)
                [f1,f2,f3,m1,m2,m3] = vals2
                if eid==0:
                    eid=''

                msg.append('%s  %8s    %10s    %8s      %s  %s  %s  %s  %s  %-s\n' %(zero,eKey,eid,elemName,f1,f2,f3,m1,m2,m3))
                zero=' '
            ###
        ###
        msg.append(pageStamp+str(pageNum)+'\n')
        return (''.join(msg),pageNum)
    
    def writeF06Transient(self,header,pageStamp,pageNum=1):
        msg = header+['                                          G R I D   P O I N T   F O R C E   B A L A N C E\n',
               ' \n',
               '   POINT-ID    ELEMENT-ID     SOURCE             T1             T2             T3             R1             R2             R3\n',]
              #'0     13683          3736    TRIAX6         4.996584E+00   0.0            1.203093E+02   0.0            0.0            0.0'
              #'      13683          3737    TRIAX6        -4.996584E+00   0.0           -1.203093E+02   0.0            0.0            0.0'
              #'      13683                  *TOTALS*       6.366463E-12   0.0           -1.364242E-12   0.0            0.0            0.0'
        for dt,Forces in sorted(self.forces.items()):
            for eKey,force in sorted(Forces.items()):
                zero = '0'
                for iLoad,f in enumerate(force):
                    (f1,f2,f3) = f
                    (m1,m2,m3) = self.moments[dt][eKey][iLoad]
                    (elemName) = self.elemName[dt][eKey][iLoad]
                    eid = self.eids[dt][eKey][iLoad]

                    vals = [f1,f2,f3,m1,m2,m3]
                    (vals2,isAllZeros) = self.writeF06Floats13E(vals)
                    [f1,f2,f3,m1,m2,m3] = vals2
                    if eid==0:
                        eid=''

                    msg.append('%s  %8s    %10s    %8s      %s  %s  %s  %s  %s  %-s\n' %(zero,eKey,eid,elemName,f1,f2,f3,m1,m2,m3))
                    zero=' '
                ###
            ###
            msg.append(pageStamp+str(pageNum)+'\n')
            pageNum+=1
        return (''.join(msg),pageNum-1)
    
    def __repr__(self):
        return self.writeF06([],'PAGE ',1)[0]
        #return '---gridPointForceObject---'

class complexGridPointForcesObject(scalarObject):
    def __init__(self,dataCode,iSubcase,freq=None):
        scalarObject.__init__(self,dataCode,iSubcase)
class complexDisplacementObject(scalarObject): # approachCode=1, sortCode=0, thermal=0
    def __init__(self,dataCode,iSubcase,freq=None):
        scalarObject.__init__(self,dataCode,iSubcase)
        self.freq = freq
        #print "complexDisplacementObject - self.freq=|%s|" %(self.freq)
        self.gridType     = {}
        self.translations = {}
        self.rotations    = {}
        self.addNewTransient()

    def updateDt(self,dataCode,freq):
        self.dataCode = dataCode
        self.applyDataCode()
        if freq is not None:
            self.log.debug("updating %s...%s=%s  iSubcase=%s" %(self.dataCode['name'],self.dataCode['name'],freq,self.iSubcase))
            self.freq = freq
            self.addNewTransient()
        ###

    def deleteTransient(self,dt):
        del self.translations[dt]
        del self.rotations[dt]

    def getTransients(self):
        k = self.translations.keys()
        k.sort()
        return k

    def addNewTransient(self):
        """initializes the transient variables"""
        if self.dt not in self.translations:
            self.translations[self.freq] = {}
            self.rotations[self.freq]    = {}

    def add(self,nodeID,gridType,v1r,v1i,v2r,v2i,v3r,v3i,v4r,v4i,v5r,v5i,v6r,v6i):
        msg = "nodeID=%s v1r=%s v2r=%s v3r=%s" %(nodeID,v1r,v2r,v3r)
        #print msg
        #msg = ''
        assert 0<nodeID<1000000000, msg
        #assert nodeID not in self.translations,'complexDisplacementObject - static failure'

        self.translations[self.freq][nodeID] = [[v1r,v1i],[v2r,v2i],[v3r,v3i]] # dx,dy,dz
        self.rotations[self.freq][nodeID]    = [[v4r,v4i],[v5r,v5i],[v6r,v6i]] # rx,ry,rz
    ###

    def __repr__(self):
        msg = '---COMPLEX DISPLACEMENTS---\n'
        #if self.dt is not None:
        #    msg += '%s = %g\n' %(self.dataCode['name'],self.dt)
        headers = ['DxReal','DxImag','DyReal','DyImag','DzReal','DyImag','RxReal','RxImag','RyReal','RyImag','RzReal','RzImag']
        msg += '%-10s ' %('nodeID')
        for header in headers:
            msg += '%10s ' %(header)
        msg += '\n'

        for freq,translations in sorted(self.translations.items()):
            msg += 'freq = %g\n' %(freq)
            #print "freq = ",freq
            #print translations
            for nodeID,translation in sorted(translations.items()):
                rotation = self.rotations[freq][nodeID]
                (dx,dy,dz) = translation
                (rx,ry,rz) = rotation

                msg += '%-10i ' %(nodeID)
                vals = dx+dy+dz+rx+ry+rz
                for val in vals:
                    if abs(val)<1e-6:
                        msg += '%10s ' %(0)
                    else:
                        msg += '%10.3e ' %(val)
                    ###
                msg += '\n'
            ###
        return msg

#---------------------------------------------------------------------------------
#class staticFluxObj(scalarObject): # approachCode=1, tableCode=3 - whatever the static version of this is...

class fluxObject(scalarObject): # approachCode=1, tableCode=3, thermal=1
    def __init__(self,dataCode,iSubcase,dt=None):
        scalarObject.__init__(self,dataCode,iSubcase)

        self.dt = dt
        self.fluxes = {}
        if dt is not None:
            self.fluxes = {}
            self.isTransient = True
            raise Exception('transient is supported yet...')

    def deleteTransient(self,dt):
        del self.fluxes[dt]

    def getTransients(self):
        k = self.fluxes.keys()
        k.sort()
        return k

    def add(self,nodeID,gridType,v1,v2,v3,v4=None,v5=None,v6=None):
        assert 0<nodeID<1000000000, 'nodeID=%s' %(nodeID)
        assert nodeID not in self.fluxes
        self.fluxes[nodeID] = array([v1,v2,v3])

    def writeOp2(self,block3,deviceCode=1):
        """
        creates the binary data for writing the table
        @warning hasnt been tested...
        """
        msg = block3
        for nodeID,flux in sorted(self.fluxes.items()):
            grid = nodeID*10+deviceCode
            msg += pack('iffffff',grid,flux[0],flux[1],flux[2],0,0,0)
        ###
        return msg

    def __repr__(self):
        if self.isTransient:
            return self.__reprTransient__()

        msg = '---HEAT FLUX---\n'
        msg += '%-10s %-8s %-8s %-8s\n' %('NodeID','xFlux','yFlux','zFlux')
        for nodeID,flux in sorted(self.fluxes.items()):
            msg += '%10i ' %(nodeID)

            for val in flux:
                if abs(val)<1e-6:
                    msg += '%10s' %(0)
                else:
                    msg += '%10.3e ' %(val)
                ###
            msg += '\n'
        return msg

