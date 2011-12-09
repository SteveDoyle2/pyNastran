from struct import pack
from pyNastran.op2.resultObjects.op2_Objects import scalarObject,array

class displacementObject(scalarObject): # approachCode=1, sortCode=0, thermal=0
    def __init__(self,dataCode,iSubcase,dt=None):
        scalarObject.__init__(self,dataCode,iSubcase)
        self.dt = dt
        print "displacementObject - self.dt=|%s|" %(self.dt)
        ## this could get very bad very quick, but it could be great!
        ## basically it's a way to handle transients without making
        ## a whole new class
        self.displacements = {}
        self.rotations     = {}
        if dt is not None:
            #assert dt>=0.
            self.addNewTransient()
            self.add = self.addTransient
            #self.addBinary = self.addBinaryTransient
            #self.__repr__ = self.__reprTransient__  # why cant i do this...
            #self.writeOp2 = self.writeOp2Transient
        ###

    #def setLoadID(self,loadID):
    #    self.loadID = loadID

    def addNewTransient(self):
        """
        initializes the transient variables
        @note make sure you set self.dt first
        """
        self.displacements[self.dt] = {}
        self.rotations[self.dt]     = {}

    def addBinary(self,deviceCode,data):
        (nodeID,v1,v2,v3,v4,v5,v6) = unpack('iffffff',data)
        msg = "nodeID=%s v1=%s v2=%s v3=%s" %(nodeID,v1,v2,v3)
        assert 0<nodeID<1000000000, msg
        assert nodeID not in self.displacements

        self.displacements[nodeID] = array([v1,v2,v3]) # dx,dy,dz
        self.rotations[nodeID]     = array([v4,v5,v6]) # rx,ry,rz
    ###

    def add(self,nodeID,gridType,v1,v2,v3,v4,v5,v6):
        msg = "nodeID=%s v1=%s v2=%s v3=%s" %(nodeID,v1,v2,v3)
        assert 0<nodeID<1000000000, msg
        #assert nodeID not in self.displacements,'displacementObject - static failure'

        self.displacements[nodeID] = array([v1,v2,v3]) # dx,dy,dz
        self.rotations[nodeID]     = array([v4,v5,v6]) # rx,ry,rz
    ###

    def addTransient(self,nodeID,gridType,v1,v2,v3,v4,v5,v6):
        msg  = "nodeID=%s v1=%s v2=%s v3=%s\n" %(nodeID,v1,v2,v3)
        msg += "          v4=%s v5=%s v6=%s"   %(       v4,v5,v6)
        assert 0<nodeID<1000000000, msg
        assert nodeID not in self.displacements[self.dt],'displacementObject - transient failure'
        
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

    def __reprTransient__(self):
        msg = '---DISPLACEMENTS---\n'
        #if self.dt is not None:
        #    msg += 'dt = %g\n' %(self.dt)
        headers = ['Dx','Dy','Dz','Rx','Ry','Rz']
        msg += '%-8s ' %('nodeID')
        for header in headers:
            msg += '%10s ' %(header)
        msg += '\n'

        for dt,displacements in sorted(self.displacements.items()):
            msg += 'dt = %g\n' %(dt)
            for nodeID,displacement in sorted(displacements.items()):
                rotation = self.rotations[dt][nodeID]
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
            ###
        return msg

    def __repr__(self):
        if self.dt is not None:
            return self.__reprTransient__()

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

class temperatureObject(scalarObject): # approachCode=1, sortCode=0, thermal=1
    def __init__(self,dataCode,iSubcase,dt=None):
        scalarObject.__init__(self,dataCode,iSubcase)
        self.dt = dt
        self.temperatures = {}
        
        #print "dt = ",self.dt
        self.temperatures = {}
        if dt is not None:
            #assert dt>=0.
            self.addNewTransient()
            self.isTransient = True
            self.temperatures = {self.dt:{}}
            self.add = self.addTransient
            self.addBinary = self.addBinaryTransient
            #self.__repr__ = self.__reprTransient__  # why cant i do this...            
        ###

    def addNewTransient(self):
        """
        initializes the transient variables
        @note make sure you set self.dt first
        """
        self.temperatures[self.dt] = {}

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

    def add(self,nodeID,gridType,v1,v2=None,v3=None,v4=None,v5=None,v6=None):
        assert 0<nodeID<1000000000, 'nodeID=%s' %(nodeID)
        assert nodeID not in self.temperatures
        self.temperatures[nodeID] = v1

    def addTransient(self,nodeID,gridType,v1,v2=None,v3=None,v4=None,v5=None,v6=None):
        assert 0<nodeID<1000000000, 'nodeID=%s' %(nodeID)
        assert nodeID not in self.temperatures[self.dt]
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
        if self.isTransient:
            return self.__reprTransient__()

        msg = '---TEMPERATURE---\n'
        msg += '%-8s %-8s\n' %('GRID','TEMPERATURE')
        #print "self.dataCode=",self.dataCode
        #print "self.temperatures=",self.temperatures
        for nodeID,T in sorted(self.temperatures.items()):
            msg += '%9i ' %(nodeID)

            if abs(T)<1e-6:
                msg += '%10s\n' %(0)
            else:
                msg += '%10g\n' %(T)
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

#---------------------------------------------------------------------------------
class eigenVectorObject(scalarObject): # approachCode=2, sortCode=0, thermal=0
    """
    EIGENVALUE =  6.158494E+07
        CYCLES =  1.248985E+03         R E A L   E I G E N V E C T O R   N O .          1

    POINT ID.   TYPE          T1             T2             T3             R1             R2             R3
           1      G      2.547245E-17  -6.388945E-16   2.292728E+00  -1.076928E-15   2.579163E-17   0.0
        2002      G     -6.382321E-17  -1.556607E-15   3.242408E+00  -6.530917E-16   1.747180E-17   0.0
        2003      G      0.0            0.0            0.0            0.0            0.0            0.0
    """
    def __init__(self,dataCode,iSubcase,eigReal):
        scalarObject.__init__(self,dataCode,iSubcase)
        #self.eigReal = int(eigReal)
        self.eigReal = eigReal
        self.updateDt = self.updateEigReal
        #print "eigReal = %s" %(eigReal)
        
        #assert eigReal>=0.
        self.gridTypes = {}
        self.displacements = {self.eigReal: {}}
        self.rotations     = {self.eigReal: {}}

    def updateEigReal(self,eigReal):
        """
        this method is called if the object
        already exits and a new time step is found
        """
        #assert eigReal>=0.
        #self.eigReal = int(eigReal)
        self.eigReal = eigReal
        print "eigReal = %s" %(str(eigReal))
        self.displacements[self.eigReal] = {}
        self.rotations[self.eigReal] = {}

        self.appendDataMember('modes','mode')
        self.appendDataMember('eigrs','eigr')
        self.appendDataMember('eigis','eigi')

    def add(self,nodeID,gridType,v1,v2,v3,v4,v5,v6):
        msg = "nodeID=%s v1=%s v2=%s v3=%s" %(nodeID,v1,v2,v3)
        assert 0<nodeID<1000000000, msg
        #assert nodeID not in self.displacements

        #if gridType==0:
        #    Type = 'S'
        if gridType==1:
            Type = 'G'
        elif gridType==2:
            Type = 'S'
        else:
            raise Exception('invalid grid type...gridType=%s' %(gridType))

        self.gridTypes[nodeID] = Type
        #print 'self.eigReal = %s' %(self.eigReal),type(self.eigReal)
        #print "d = ",self.displacements
        self.displacements[self.eigReal][nodeID] = array([v1,v2,v3]) # dx,dy,dz
        self.rotations[self.eigReal][nodeID]     = array([v4,v5,v6]) # rx,ry,rz
    ###

    def eigenvalues(self):
        return sorted(self.displacements.keys())

    def __repr__(self):
        msg = '---EIGENVECTORS---\n'
        msg += '-eigenvalues-\n'
        for i,eigenvalue in enumerate(self.eigenvalues()):
            msg += '%-2s %f\n' %(i,eigenvalue)
        msg += '\n'

        #if self.eigReal is not None:
        #    msg += 'eigReal = %g\n' %(self.eigReal)
        headers = ['Tx','Ty','Tz','Rx','Ry','Rz']
        headerLine = '%-8s %8s ' %('nodeID','GridType',)
        for header in headers:
            headerLine += '%10s ' %(header)
        headerLine += '\n'

        for freq,eigVals in sorted(self.displacements.items()):
            msg += 'eigenvalueReal = %f\n' %(freq)
            msg += headerLine
            for nodeID,displacement in sorted(eigVals.items()):
                rotation = self.rotations[freq][nodeID]
                Type = self.gridTypes[nodeID]
                (dx,dy,dz) = displacement
                (rx,ry,rz) = rotation

                msg += '%-8i %8s ' %(nodeID,Type)
                vals = [dx,dy,dz,rx,ry,rz]
                for val in vals:
                    if abs(val)<1e-6:
                        msg += '%10s ' %(0)
                    else:
                        msg += '%10.3g ' %(val)
                    ###
                msg += '\n'
            ###
            msg += '\n'
        return msg

class eigenVector14Object(scalarObject): # approachCode=2, sortCode=0, thermal=0
    """
    @todo add table data
    """
    def __init__(self,dataCode,iSubcase,eigReal):
        scalarObject.__init__(self,dataCode,iSubcase)
        #self.eigReal = int(eigReal)
        self.eigReal = eigReal
        self.updateDt = self.updateEigReal
        #print "eigReal = %s" %(eigReal)
        
        #assert eigReal>=0.
        self.gridTypes = {}
        self.displacements = {self.eigReal: {}}
        self.rotations     = {self.eigReal: {}}

    def updateEigReal(self,eigReal):
        """
        this method is called if the object
        already exits and a new time step is found
        """
        #assert eigReal>=0.
        #self.eigReal = int(eigReal)
        self.eigReal = eigReal
        print "eigReal = %s" %(str(eigReal))
        self.displacements[self.eigReal] = {}
        self.rotations[self.eigReal] = {}

        self.appendDataMember('modes','mode')
        self.appendDataMember('eigrs','eigr')
        self.appendDataMember('eigis','eigi')

    def add(self,nodeID,gridType,v1,v2,v3,v4,v5,v6):
        msg = "nodeID=%s v1=%s v2=%s v3=%s v4=%s v5=%s v6=%s" %(nodeID,v1,v2,v3,v4,v5,v6)
        #assert 0<nodeID<1000000000, msg
        #assert nodeID not in self.displacements

        #if gridType==0:
        #    Type = 'S??'
        if gridType==1:
            Type = 'G'
        elif gridType==2:
            Type = 'S'
        else:
            raise Exception('invalid grid type...gridType=%s' %(gridType))

        self.gridTypes[nodeID] = Type
        #print 'self.eigReal = %s' %(self.eigReal),type(self.eigReal)
        #print "d = ",self.displacements
        self.displacements[self.eigReal][nodeID] = array([v1,v2,v3]) # dx,dy,dz
        self.rotations[self.eigReal][nodeID]     = array([v4,v5,v6]) # rx,ry,rz
    ###

    def eigenvalues(self):
        return sorted(self.displacements.keys())

    def __repr__(self):
        msg = '---EIGENVECTORS---\n'
        msg += '-eigenvalues-\n'
        for i,eigenvalue in enumerate(self.eigenvalues()):
            msg += '%-2s %f\n' %(i,eigenvalue)
        msg += '\n'

        #if self.eigReal is not None:
        #    msg += 'eigReal = %g\n' %(self.eigReal)
        headers = ['Tx','Ty','Tz','Rx','Ry','Rz']
        headerLine = '%-8s %8s ' %('nodeID','GridType',)
        for header in headers:
            headerLine += '%10s ' %(header)
        headerLine += '\n'

        for freq,eigVals in sorted(self.displacements.items()):
            msg += 'eigenvalueReal = %f\n' %(freq)
            msg += headerLine
            for nodeID,displacement in sorted(eigVals.items()):
                rotation = self.rotations[freq][nodeID]
                Type = self.gridTypes[nodeID]
                (dx,dy,dz) = displacement
                (rx,ry,rz) = rotation

                msg += '%-8i %8s ' %(nodeID,Type)
                vals = [dx,dy,dz,rx,ry,rz]
                for val in vals:
                    if abs(val)<1e-6:
                        msg += '%10s ' %(0)
                    else:
                        msg += '%10.3g ' %(val)
                    ###
                msg += '\n'
            ###
            msg += '\n'
        return msg

#---------------------------------------------------------------------------------

class nonlinearDisplacementObject(scalarObject): # approachCode=10, sortCode=0, thermal=0
    def __init__(self,dataCode,iSubcase,loadStep):
        scalarObject.__init__(self,dataCode,iSubcase)
        self.loadStep = loadStep
        
        assert loadStep>=0.
        self.displacements = {loadStep: {}}
        raise Exception('not implemented')

class nonlinearTemperatureObject(scalarObject): # approachCode=10, sortCode=0, thermal=1
    def __init__(self,dataCode,iSubcase,loadStep):
        scalarObject.__init__(self,dataCode,iSubcase)
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

    def add(self,nodeID,gridType,v1,v2,v3,v4,v5,v6): # addTransient
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
