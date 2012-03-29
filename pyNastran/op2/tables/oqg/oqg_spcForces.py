from numpy import array
from pyNastran.op2.resultObjects.op2_Objects import scalarObject

class spcForcesObject(scalarObject):
    def __init__(self,dataCode,iSubcase,dt=None):
        scalarObject.__init__(self,dataCode,iSubcase)
        self.dt = dt

        self.forces  = {}
        self.moments = {}
        if self.dt is not None:
            assert dt>=0.
            self.addNewTransient()
            self.add = self.addTransient
            #self.__repr__ = self.__reprTransient__  # why cant i do this...
        ###

    def updateDt(self,dataCode,dt):
        self.dataCode = dataCode
        self.applyDataCode()
        #assert dt>=0.
        #print "updating dt...dt=%s" %(dt)
        if dt is not None:
            self.dt = dt
            self.addNewTransient()
        ###

    def addNewTransient(self):
        self.forces[self.dt]  = {}
        self.moments[self.dt] = {}

    def updateDt(self,dataCode,dt):
        self.dataCode = dataCode
        self.applyDataCode()
        #assert dt>=0.
        #print "updating dt...dt=%s" %(dt)
        if dt is not None:
            self.dt = dt
            self.addNewTransient()
        ###

    #def addBinary(self,deviceCode,data):
        #print "*******add********"
        #(nodeID,v1,v2,v3,v4,v5,v6) = unpack('iffffff',data)

    def add(self,nodeID,gridType,v1,v2,v3,v4,v5,v6):
        msg = 'nodeID=%s' %(nodeID)
        assert 0<nodeID<1000000000,msg
        assert nodeID not in self.forces
        self.forces[ nodeID] = array([v1,v2,v3]) # Fx,Fy,Fz
        self.moments[nodeID] = array([v4,v5,v6]) # Mx,My,Mz

    def addTransient(self,nodeID,gridType,v1,v2,v3,v4,v5,v6):
        msg = 'nodeID=%s' %(nodeID)
        assert 0<nodeID<1000000000,msg
        assert nodeID not in self.forces[self.dt]
        self.forces[ self.dt][nodeID] = array([v1,v2,v3]) # Fx,Fy,Fz
        self.moments[self.dt][nodeID] = array([v4,v5,v6]) # Mx,My,Mz

    def __reprTransient__(self):
        msg = '---SPC FORCES---\n'
        if self.dt is not None:
            msg += 'dt = %g\n' %(self.dt)

        headers = ['Fx','Fy','Fz','Mx','My','Mz']
        msg += '%-8s ' %('GRID')
        for header in headers:
            msg += '%10s ' %(header)
        msg += '\n'

        for dt,forces in sorted(self.forces.items()):
            msg += 'dt = %s' %(dt)
            for nodeID,force in sorted(forces.items()):
                moment = self.moments[dt][nodeID]
                (Fx,Fy,Fz) = force
                (Mx,My,Mz) = moment

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
        if self.dt is not None:
            return self.__reprTransient__()

        msg = '---SPC FORCES---\n'
        if self.dt is not None:
            msg += 'dt = %g\n' %(self.dt)

        headers = ['Fx','Fy','Fz','Mx','My','Mz']
        msg += '%-8s ' %('GRID')
        for header in headers:
            msg += '%10s ' %(header)
        msg += '\n'

        for nodeID,force in sorted(self.forces.items()):
            moment = self.moments[nodeID]
            (Fx,Fy,Fz) = force
            (Mx,My,Mz) = moment

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

class mpcForcesObject(scalarObject):
    def __init__(self,dataCode,iSubcase,dt=None):
        scalarObject.__init__(self,dataCode,iSubcase)
        self.dt = dt

        self.forces  = {}
        self.moments = {}
        if self.dt is not None:
            assert dt>=0.
            self.addNewTransient()
            self.add = self.addTransient
            #self.__repr__ = self.__reprTransient__  # why cant i do this...
        ###

    def updateDt(self,dataCode,dt):
        self.dataCode = dataCode
        self.applyDataCode()
        #assert dt>=0.
        #print "updating dt...dt=%s" %(dt)
        if dt is not None:
            self.dt = dt
            self.addNewTransient()
        ###

    def addNewTransient(self):
        self.forces[self.dt]  = {}
        self.moments[self.dt] = {}

    def updateDt(self,dataCode,dt):
        self.dataCode = dataCode
        self.applyDataCode()
        #assert dt>=0.
        #print "updating dt...dt=%s" %(dt)
        if dt is not None:
            self.dt = dt
            self.addNewTransient()
        ###

    #def addBinary(self,deviceCode,data):
        #print "*******add********"
        #(nodeID,v1,v2,v3,v4,v5,v6) = unpack('iffffff',data)

    def add(self,nodeID,gridType,v1,v2,v3,v4,v5,v6):
        msg = 'nodeID=%s' %(nodeID)
        assert 0<nodeID<1000000000,msg
        assert nodeID not in self.forces
        self.forces[ nodeID] = array([v1,v2,v3]) # Fx,Fy,Fz
        self.moments[nodeID] = array([v4,v5,v6]) # Mx,My,Mz

    def addTransient(self,nodeID,gridType,v1,v2,v3,v4,v5,v6):
        msg = 'nodeID=%s' %(nodeID)
        assert 0<nodeID<1000000000,msg
        assert nodeID not in self.forces[self.dt]
        self.forces[ self.dt][nodeID] = array([v1,v2,v3]) # Fx,Fy,Fz
        self.moments[self.dt][nodeID] = array([v4,v5,v6]) # Mx,My,Mz

    def __reprTransient__(self):
        msg = '---MPC FORCES---\n'
        if self.dt is not None:
            msg += 'dt = %g\n' %(self.dt)

        headers = ['Fx','Fy','Fz','Mx','My','Mz']
        msg += '%-8s ' %('GRID')
        for header in headers:
            msg += '%10s ' %(header)
        msg += '\n'

        for dt,forces in sorted(self.forces.items()):
            msg += 'dt = %s' %(dt)
            for nodeID,force in sorted(forces.items()):
                moment = self.moments[dt][nodeID]
                (Fx,Fy,Fz) = force
                (Mx,My,Mz) = moment

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
        if self.dt is not None:
            return self.__reprTransient__()

        msg = '---MPC FORCES---\n'
        if self.dt is not None:
            msg += 'dt = %g\n' %(self.dt)

        headers = ['Fx','Fy','Fz','Mx','My','Mz']
        msg += '%-8s ' %('GRID')
        for header in headers:
            msg += '%10s ' %(header)
        msg += '\n'

        for nodeID,force in sorted(self.forces.items()):
            moment = self.moments[nodeID]
            (Fx,Fy,Fz) = force
            (Mx,My,Mz) = moment

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

