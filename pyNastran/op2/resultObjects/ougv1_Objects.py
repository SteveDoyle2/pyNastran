from struct import pack
from pyNastran.op2.resultObjects.op2_Objects import scalarObject,array

class displacementObject(scalarObject): # approachCode=1, tableCode=1
    def __init__(self,iSubcase,dt=None):
        scalarObject.__init__(self,iSubcase)
        self.dt = dt
        
        ## this could get very bad very quick, but it could be great!
        ## basically it's a way to handle transients without making
        ## a whole new class
        if self.dt is None:
            self.displacements = {}
            self.rotations     = {}
        else:
            assert dt>=0.
            self.displacements = {dt: {}}
            self.rotations     = {dt: {}}
            self.add = self.addTransient
            self.addBinary = self.addBinaryTransient
            self.__repr__ = self.__reprTransient__  # why cant i do this...
            #self.writeOp2 = self.writeOp2Transient
        ###

    def updateDt(self,dt=None):
        """
        this method is called if the object
        already exits and a new time step is found
        """
        assert dt>=0.
        self.dt = dt
        self.displacements[dt] = {}

    def addBinary(self,deviceCode,data):
        (nodeID,v1,v2,v3,v4,v5,v6) = unpack('iffffff',data)
        msg = "nodeID=%s v1=%s v2=%s v3=%s" %(nodeID,v1,v2,v3)
        assert 0<nodeID<1000000000, msg
        assert nodeID not in self.displacements

        self.displacements[nodeID] = array([v1,v2,v3]) # dx,dy,dz
        self.rotations[nodeID]     = array([v4,v5,v6]) # rx,ry,rz
    ###

    def addBinaryTransient(self,deviceCode,data):
        raise Exception('not implemented...')
        (nodeID,v1,v2,v3,v4,v5,v6) = unpack('iffffff',data)
        msg = "nodeID=%s v1=%s v2=%s v3=%s" %(nodeID,v1,v2,v3)
        assert 0<nodeID<1000000000, msg
        assert nodeID not in self.displacements

        self.displacements[nodeID] = array([v1,v2,v3]) # dx,dy,dz
        self.rotations[nodeID]     = array([v4,v5,v6]) # rx,ry,rz
    ###

    def add(self,nodeID,v1,v2,v3,v4,v5,v6):
        msg = "nodeID=%s v1=%s v2=%s v3=%s" %(nodeID,v1,v2,v3)
        assert 0<nodeID<1000000000, msg
        assert nodeID not in self.displacements

        self.displacements[nodeID] = array([v1,v2,v3]) # dx,dy,dz
        self.rotations[nodeID]     = array([v4,v5,v6]) # rx,ry,rz
    ###

    def addTransient(self,nodeID,v1,v2,v3,v4,v5,v6):
        msg = "nodeID=%s v1=%s v2=%s v3=%s" %(nodeID,v1,v2,v3)
        assert 0<nodeID<1000000000, msg
        assert nodeID not in self.displacements
        
        self.displacements[self.dt][nodeID] = array([v1,v2,v3]) # dx,dy,dz
        self.rotations[self.dt][nodeID]     = array([v4,v5,v6]) # rx,ry,rz
    ###

    def __reprTransient__(self):
        print "Transient..."
        raise Exception('this could be cool...')
        return self.__repr__()

    def writeOp2(self,block3,deviceCode=1):
        """
        creates the binary data for writing the table
        @warning hasnt been tested...
        """
        msg = block3
        for nodeID,displacement in sorted(self.displacements.items()):
            rotation = self.rotations[nodeID]
            (dx,dy,dz) = displacement
            (rx,ry,rz) = rotation
            
            grid = nodeID*10+deviceCode
            msg += pack('iffffff',grid,dx,dy,dz,rx,ry,rz)
        ###
        return msg

    #def writeOp2Transient(self,block3,deviceCode=1):
    #    """
    #    creates the binary data for writing the table
    #    @warning hasnt been tested...
    #    @warning dt slot needs to be fixed...
    #    """
    #    msg = ''
    #    for dt,displacements in sorted(self.displacements.items()):
    #        XXX = 50 ## this isnt correct... @todo update dt
    #        msg += block3[0:XXX] + pack('i',dt) + block3[XXX+4:]
    #        #msg += 'dt = %g\n' %(dt)
    #
    #        for nodeID,displacement in sorted(displacements.items()):
    #            rotation = self.rotations[nodeID]
    #            (dx,dy,dz) = displacement
    #            (rx,ry,rz) = rotation
    #
    #            grid = nodeID*10+deviceCode
    #            msg += pack('iffffff',grid,dx,dy,dz,rx,ry,rz)
    #        ###
    #    ###
    #    return msg

    def __repr__(self):
        msg = '---DISPLACEMENTS---\n'
        if self.dt is not None:
            msg += 'dt = %g\n' %(self.dt)
        headers = ['Dx','Dy','Dz','Rx','Ry','Rz']
        msg += '%-8s ' %('nodeID')
        for header in headers:
            msg += '%10s ' %(header)
        msg += '\n'

        for nodeID,displacement in sorted(self.displacements.items()):
            rotation = self.rotations[nodeID]
            (dx,dy,dz) = displacement
            (rx,ry,rz) = rotation

            msg += '%-8i ' %(nodeID)
            vals = [dx,dy,dz,rx,ry,rz]
            for val in vals:
                if abs(val)<1e-6:
                    msg += '%10s ' %(0)
                else:
                    msg += '%10.3e ' %(val)
                ###
            msg += '\n'
        return msg

class nonlinearDisplacementObject(scalarObject): # approachCode=10, tableCode=1
    def __init__(self,iSubcase,loadStep):
        scalarObject.__init__(self,iSubcase)
        self.loadStep = loadStep
        
        assert loadStep>=0.
        self.displacements = {loadStep: {}}
        raise Exception('not implemented')

#---------------------------------------------------------------------------------

class temperatureObject(scalarObject): # approachCode=1, tableCode=1
    def __init__(self,iSubcase,dt=None):
        scalarObject.__init__(self,iSubcase)
        self.dt = dt
        self.temperatures = {}
        
        print "dt = ",self.dt
        if self.dt is None:
            self.temperatures = {}
        else:
            self.temperatures = {self.dt:{}}
            self.add = self.addTransient
            self.addBinary = self.addBinaryTransient
            #self.__repr__ = self.__reprTransient__  # why cant i do this...            
        ###

    def updateDt(self,dt):
        """
        this method is called if the object
        already exits and a new time step is found
        """
        assert dt>=0.
        self.dt = dt
        self.temperatures[dt] = {}

    def addBinary(self,deviceCode,data):
        (nodeID,v1) = unpack('if',data[0:8])
        assert 0<nodeID<1000000000, 'nodeID=%s' %(nodeID)
        assert nodeID not in self.temperatures
        self.temperatures[nodeID] = v1

    def addBinaryTransient(self,deviceCode,data):
        raise Exception('not implemented...')
        (nodeID,v1) = unpack('if',data[0:8])
        assert 0<nodeID<1000000000, 'nodeID=%s' %(nodeID)
        assert nodeID not in self.temperatures
        self.temperatures[nodeID] = v1

    def add(self,nodeID,v1,v2=None,v3=None,v4=None,v5=None,v6=None):
        assert 0<nodeID<1000000000, 'nodeID=%s' %(nodeID)
        assert nodeID not in self.temperatures
        self.temperatures[nodeID] = v1

    def addTransient(self,nodeID,v1,v2=None,v3=None,v4=None,v5=None,v6=None):
        assert 0<nodeID<1000000000, 'nodeID=%s' %(nodeID)
        assert nodeID not in self.temperatures
        self.temperatures[self.dt][nodeID] = v1

    def __reprTransient__(self):
        msg = '---TEMPERATURE---\n'
        msg += '%-8s %-8s\n' %('GRID','TEMPERATURE')

        for dt,temperatures in sorted(self.temperatures.items()):
            msg += 'dt = %g\n' %(dt)
            for nodeID,T in sorted(temperatures.items()):
                msg += '%9i ' %(nodeID)

                if abs(T)<1e-6:
                    msg += '%10s\n' %(0)
                else:
                    msg += '%10g\n' %(T)
                ###
            ###
        return msg

    def writeOp2(self,block3,deviceCode=1):
        """
        creates the binary data for writing the table
        @warning hasnt been tested...
        """
        msg = block3
        for nodeID,T in sorted(self.temperatures.items()):
            grid = nodeID*10+deviceCode
            msg += pack('iffffff',grid,T,0,0,0,0,0)
        ###
        return msg

    def writeOp2Transient(self,block3,deviceCode=1):
        """
        creates the binary data for writing the table
        @warning hasnt been tested...
        @warning dt slot needs to be fixed...
        """
        msg = ''
        for dt,temperatures in sorted(self.temperatures.items()):
            XXX = 50 ## this isnt correct... @todo update dt
            msg += block3[0:XXX] + pack('i',dt) + block3[XXX+4:]
            #msg += 'dt = %g\n' %(dt)
    
            for nodeID,T in sorted(temperatures.items()):
                grid = nodeID*10+deviceCode
                msg += pack('iffffff',grid,T,0,0,0,0,0)
            ###
        ###
        return msg

    def __repr__(self):
        if self.dt is not None:
            return self.__reprTransient__()

        msg = '---TEMPERATURE---\n'
        msg += '%-8s %-8s\n' %('GRID','TEMPERATURE')
        for nodeID,T in sorted(self.temperatures.items()):
            msg += '%9i ' %(nodeID)

            if abs(T)<1e-6:
                msg += '%10s\n' %(0)
            else:
                msg += '%10g\n' %(T)
            ###
        return msg

class fluxObject(scalarObject): # approachCode=1, tableCode=3
    def __init__(self,iSubcase,dt=None):
        scalarObject.__init__(self,iSubcase)

        self.dt = dt
        if dt is not None:
            self.fluxes = {}
        else:
            raise Exception('transient not supported for flux yet...')
        

    def add(self,nodeID,v1,v2,v3,v4=None,v5=None,v6=None):
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
        msg = '---HEAT FLUX---\n'
        msg += '%-8s %-8s %-8s %-8s\n' %('GRID','xFlux','yFlux','zFlux')
        for nodeID,flux in sorted(self.fluxes.items()):
            msg += '%9i ' %(nodeID)

            for val in flux:
                if abs(val)<1e-6:
                    msg += '%10s' %(0)
                else:
                    msg += '%10.3e ' %(val)
                ###
            msg += '\n'
        return msg

class nonlinearTemperatureObject(scalarObject): # approachCode=10, tableCode=1
    def __init__(self,iSubcase,loadStep):
        scalarObject.__init__(self,iSubcase)
        self.loadStep = loadStep
        
        assert loadStep>=0.
        self.temperatures = {loadStep: {}}

    def updateDt(self,loadStep):
        """
        this method is called if the object
        already exits and a new time step is found
        """
        assert loadStep>=0.
        self.loadStep = loadStep
        self.temperatures[loadStep] = {}

    def add(self,nodeID,v1,v2,v3,v4,v5,v6): # addTransient
        #msg = "nodeID=%s v1=%s v2=%s v3=%s v4=%s v5=%s v6=%s" %(nodeID,v1,v2,v3,v4,v5,v6)
        msg = "nodeID=%s v1=%s" %(nodeID,v1)
        #print msg
        assert 0<nodeID<1000000000, msg
        assert nodeID not in self.temperatures[self.loadStep]
        
        self.temperatures[self.loadStep][nodeID] = v1 # T
    ###

    def __repr__(self):
        msg = '---TEMPERATURE VECTOR---\n'
        if self.loadStep is not None:
            msg += 'loadStep = %g\n' %(self.loadStep)
        headers = ['Temperature']
        msg += '%9s ' %('GRID')
        for header in headers:
            msg += '%10s ' %(header)
        msg += '\n'

        for dt,temps in sorted(self.temperatures.items()):
            for nodeID,T in sorted(temps.items()):
                msg += '%9i ' %(nodeID)
                if abs(T)<1e-6:
                    msg += '%10s ' %(0)
                else:
                    msg += '%10.3f ' %(T)
                ###
                msg += '\n'
            ###
        return msg
