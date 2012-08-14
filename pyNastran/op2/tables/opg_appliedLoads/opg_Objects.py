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
from numpy import array
from pyNastran.op2.resultObjects.op2_Objects import scalarObject

class AppliedLoadsObject(scalarObject): # approachCode=1, sortCode=0
    def __init__(self,dataCode,isSort1,iSubcase,dt=None):
        scalarObject.__init__(self,dataCode,iSubcase)
        self.dt = dt
        
        ## this could get very bad very quick, but it could be great!
        ## basically it's a way to handle transients without making
        ## a whole new class
        self.eids    = {}
        self.sources = {}
        self.forces  = {}
        self.moments = {}
        if self.dt is not None:
            assert dt>=0.
            raise NotImplementedError('transient appliedLoads not implemented...')
            self.eids    = {dt: []}
            self.sources = {dt: []}
            self.forces  = {dt: []}
            self.moments = {dt: []}
            self.add = self.addTransient
        ###

    def updateDt(self):
        """
        this method is called if the object
        already exits and a new time step is found
        """
        assert dt>=0.
        self.dt = dt
        self.eids[dt]    = []
        self.sources[dt] = []
        self.forces[dt]  = []
        self.moments[dt] = []

    def addNode(self,nodeID,eid,source,v1,v2,v3,v4,v5,v6):
        msg = "nodeID=%s eid=%s source=|%s| v1=%i v2=%i v3=%i v4=%i v5=%i v6=%i" %(nodeID,eid,source,v1,v2,v3,v4,v5,v6)
        assert 0<nodeID<1000000000, msg
        assert nodeID not in self.forces,msg

        self.eids[nodeID]    = [eid]
        self.sources[nodeID] = [source]
        self.forces[ nodeID] = [ array([v1,v2,v3]) ] # Fx,Fy,Fz
        self.moments[nodeID] = [ array([v4,v5,v6]) ] # Mx,My,Mz
    ###

    def add(self,nodeID,eid,source,v1,v2,v3,v4,v5,v6):
        msg = "nodeID=%s eid=%s source=|%s| v1=%i v2=%i v3=%i v4=%i v5=%i v6=%i" %(nodeID,eid,source,v1,v2,v3,v4,v5,v6)
        assert 0<nodeID<1000000000,msg
        if nodeID not in self.forces:
            self.addNode(nodeID,eid,source,v1,v2,v3,v4,v5,v6)
            return None
        #assert nodeID not in self.forces,msg

        self.eids[nodeID].append(eid)
        self.sources[nodeID].append(source)
        self.forces[nodeID].append( array([v1,v2,v3])) # Fx,Fy,Fz
        self.moments[nodeID].append(array([v4,v5,v6])) # Mx,My,Mz
    ###

    def addTransient(self,nodeID,eid,source,v1,v2,v3,v4,v5,v6):
        raise Exception('no implemented')
        msg = "nodeID=%s v1=%s v2=%s v3=%s" %(nodeID,v1,v2,v3)
        assert 0<nodeID<1000000000, msg
        assert nodeID not in self.forces
        
        self.forces[ self.dt][nodeID] = array([v1,v2,v3]) # Fx,Fy,Fz
        self.moments[self.dt][nodeID] = array([v4,v5,v6]) # Mx,My,Mz
    ###

    def __repr__(self):
        msg = '---APPLIED LOADS---\n'
        if self.dt is not None:
            msg += 'dt = %g\n' %(self.dt)
        headers = ['Fx','Fy','Fz','Mx','My','Mz']
        msg += '%-8s %-8s %8s ' %('nodeID','eID','source')
        for header in headers:
            msg += '%11s ' %(header)
        msg += '\n'

        msg += '-'*100+'\n'
        for nodeID,forces in sorted(self.forces.iteritems()):
            for i in xrange(len(forces)):
                force  = forces[i]
                moment = self.moments[nodeID][i]
                source = self.sources[nodeID][i]
                eid    = self.eids[nodeID][i]
                if '*TOTALS*' in source:
                    eid = '='
                elif eid==0:
                    eid = '*'

                (Fx,Fy,Fz) = force
                (Mx,My,Mz) = moment

                msg += '%-8i %-8s %8s ' %(nodeID,eid,source)
                vals = [Fx,Fy,Fz,Mx,My,Mz]
                for val in vals:
                    if abs(val)<1e-6:
                        msg += '%11s ' %(0)
                    else:
                        msg += '%11.3f ' %(val)
                    ###
                msg += '\n'
                if '*TOTALS*' in source:
                    msg += '-'*100+'\n'
                ###
            ###
        return msg

