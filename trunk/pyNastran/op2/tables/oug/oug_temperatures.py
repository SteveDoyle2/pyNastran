import sys
from struct import pack
from pyNastran.op2.resultObjects.op2_Objects import scalarObject,array

class temperatureObject(scalarObject): # approachCode=1, sortCode=0, thermal=1
    def __init__(self,dataCode,iSubcase,dt=None):
        scalarObject.__init__(self,dataCode,iSubcase)
        self.dt = dt
        self.gridTypes    = {}
        self.temperatures = {}
        
        #print "dt = ",self.dt
        if dt is not None:
            #assert dt>=0.
            self.addNewTransient()
            self.isTransient = True
            self.addNewTransient()
            #self.temperatures = {self.dt:{}}
            self.add = self.addTransient
            #self.addBinary = self.addBinaryTransient
            #self.__repr__ = self.__reprTransient__  # why cant i do this...            
        ###

    def updateDt(self,dataCode,dt):
        self.dataCode = dataCode
        self.applyDataCode()
        if dt is not None:
            self.log.debug("updating %s...%s=%s  iSubcase=%s" %(self.name,self.name,dt,self.iSubcase))
            self.dt = dt
            self.addNewTransient()
        ###

    def addNewTransient(self):
        """initializes the transient variables"""
        if self.dt not in self.temperatures:
            self.temperatures[self.dt] = {}

    def add(self,nodeID,gridType,v1,v2=None,v3=None,v4=None,v5=None,v6=None):
        assert 0<nodeID<1000000000, 'nodeID=%s' %(nodeID)
        #assert nodeID not in self.temperatures

        gridType = self.recastGridType(gridType)
        self.gridTypes[nodeID] = gridType
        self.temperatures[nodeID] = v1

    def addTransient(self,nodeID,gridType,v1,v2=None,v3=None,v4=None,v5=None,v6=None):
        assert 0<nodeID<1000000000, 'nodeID=%s' %(nodeID)
        #assert nodeID not in self.temperatures[self.dt]

        gridType = self.recastGridType(gridType)
        self.gridTypes[nodeID] = gridType
        self.temperatures[self.dt][nodeID] = v1

    def __reprTransient__(self):
        msg = '---TRANSIENT TEMPERATURE---\n'
        msg += '%-10s %8s %-8s\n' %('NodeID','GridType','Temperature')

        for dt,temperatures in sorted(self.temperatures.items()):
            msg += '%s = %g\n' %(self.dataCode['name'],dt)
            for nodeID,T in sorted(temperatures.items()):
                gridType = self.gridTypes[nodeID]
                msg += '%10s %8s ' %(nodeID,gridType)

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
            #msg += '%s = %g\n' %(self.dataCode['name'],dt)
    
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
        msg += '%-10s %8s %-8s\n' %('NodeID','GridType','Temperature')
        #print "self.dataCode=",self.dataCode
        #print "self.temperatures=",self.temperatures
        for nodeID,T in sorted(self.temperatures.items()):
            gridType = self.gridTypes[nodeID]
            msg += '%10s %8s ' %(nodeID,gridType)

            if abs(T)<1e-6:
                msg += '%10s\n' %(0)
            else:
                msg += '%10g\n' %(T)
            ###
        return msg

