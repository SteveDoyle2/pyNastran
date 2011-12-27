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
            self.addNewTransient()
            #assert dt>=0.
            self.isTransient = True
            self.add = self.addTransient
            self.addF = self.addTransientF
            #self.addBinary = self.addBinaryTransient
            #self.__repr__ = self.__reprTransient__  # why cant i do this...            
        ###
        self.parseLength()

    def parseLength(self):
        self.mainHeaders = []
        self.strFormat = ''
        if self.analysisCode==6:
            self.mainHeaders.append('Freq')
            self.strFormat += 'fi'
            self.add = self.addF
        elif self.analysisCode in[5,10]:
            self.mainHeaders.append('Time')
            self.strFormat += 'fi'
            self.add = self.addF
        elif self.analysisCode in [1,2,3,4,7,8,9,11,12]:
            self.mainHeaders.append('NodeID')
            self.strFormat += 'ii'
        else:
            raise Exception('invalid analysisCode=%s' %(self.analysisCode))
        self.mainHeaders.append('GridType')

        if self.isImaginary():  # elif self.dataFormat==1:
            self.strFormat += 'ffffffffffff'
            self.headers = ['Temperature']
            raise Exception('verify...add imaginary...')
        else:
            self.strFormat += 'ffffff'         # if self.dataFormat in [0,2]:
            self.headers = ['Temperature']
        
        self.mainHeaders = tuple(self.mainHeaders)

    def getLength(self):
        return (4*len(self.strFormat),self.strFormat)

    def updateDt(self,dataCode,dt):
        self.dataCode = dataCode
        self.applyDataCode()
        if dt is not None:
            self.log.debug("updating %s...%s=%s  iSubcase=%s" %(self.dataCode['name'],self.dataCode['name'],dt,self.iSubcase))
            self.dt = dt
            self.addNewTransient()
        ###

    def addNewTransient(self):
        """initializes the transient variables"""
        if self.dt not in self.temperatures:
            self.temperatures[self.dt] = {}

    def add(self,out):
        (nodeID,gridType,v1,v2,v3,v4,v5,v6) = out  # v2-v6 are 0
        nodeID = (nodeID-self.deviceCode) // 10
        assert 0<nodeID<1000000000, 'nodeID=%s' %(nodeID)
        #assert nodeID not in self.temperatures

        gridType = self.recastGridType(gridType)
        self.gridTypes[nodeID] = gridType
        self.temperatures[nodeID] = v1
    ###

    def addF(self,out):
        (freq,gridType,v1,v2,v3,v4,v5,v6) = out
        msg = "dt=%g %s=%s gridType=%s v1=%s v2=%s v3=%s" %(self.dt,self.mainHeaders[0],time,gridType,v1,v2,v3)
        #print msg
        #assert 0<nodeID<1000000000, msg
        #assert nodeID not in self.translations,'temperatrueObject - static failure'
        
        gridType = self.recastGridType(gridType)
        self.gridTypes[freq]    = gridType
        self.temperatures[nodeID] = v1
    ###

    def addTransient(self,out):
        (nodeID,gridType,v1,v2,v3,v4,v5,v6) = out  # v2-v6 are 0
        nodeID = (nodeID-self.deviceCode) // 10
        assert 0<nodeID<1000000000, 'nodeID=%s' %(nodeID)
        #assert nodeID not in self.temperatures[self.dt]

        gridType = self.recastGridType(gridType)
        self.gridTypes[nodeID] = gridType
        self.temperatures[self.dt][nodeID] = v1

    def addTransientF(self,out):
        (freq,gridType,v1,v2,v3,v4,v5,v6) = out  # v2-v6 are 0
        msg = "dt=%g %s=%s gridType=%s v1=%s v2=%s v3=%s" %(self.dt,self.mainHeaders[0],freq,gridType,v1,v2,v3)
        #print msg
        gridType = self.recastGridType(gridType)
        self.gridTypes[freq]             = gridType
        self.temperatures[self.dt][freq] = v1
    ###

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

    def writeHeader(self):
        (mainHeaders,headers) = self.getHeaders()
        msg = '%-10s %-8s ' %(mainHeaders)
        for header in headers:
            msg += '%10s ' %(header)
        msg += '\n'
        return msg

    def getHeaders(self):
        return (self.mainHeaders,self.headers)

    def __reprTransient__(self):
        msg = '---TRANSIENT TEMPERATURE---\n'
        msg += self.writeHeader()

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

    def __repr__(self):
        if self.isTransient:
            return self.__reprTransient__()

        msg = '---TEMPERATURE---\n'
        msg += self.writeHeader()
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

