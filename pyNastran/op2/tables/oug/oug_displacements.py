import sys
from struct import pack
from pyNastran.op2.resultObjects.op2_Objects import scalarObject,array

class displacementObject(scalarObject): # approachCode=1, sortCode=0, thermal=0
    def __init__(self,dataCode,iSubcase,dt=None):
        scalarObject.__init__(self,dataCode,iSubcase)
        self.dt = dt
        #print "displacementObject - self.dt=|%s|" %(self.dt)
        ## this could get very bad very quick, but it could be great!
        ## basically it's a way to handle transients without making
        ## a whole new class
        self.gridTypes    = {}
        self.translations = {}
        self.rotations    = {}
        
        if dt is not None:
            self.addNewTransient()
            self.add  = self.addTransient
            self.addF = self.addTransientF
            #self.addBinary = self.addBinaryTransient
            #self.__repr__ = self.__reprTransient__  # why cant i do this...
            #self.writeOp2 = self.writeOp2Transient
        ###
        self.parseLength()
    
    def addF06Data(self,data,transient):
        if transient is None:
            for line in data:
                (nodeID,gridType,t1,t2,t3,r1,r2,r3) = line
                self.gridTypes[nodeID]    = gridType
                self.translations[nodeID] = array([t1,t2,t3])
                self.rotations[nodeID]    = array([r1,r2,r3])
            ###
            return

        (dtName,dt) = transient
        self.dataCode['name'] = dtName
        if dt not in self.translations:
            self.updateDt(self.dataCode,dt)

        for line in data:
            (nodeID,gridType,t1,t2,t3,r1,r2,r3) = line
            self.gridTypes[nodeID]    = gridType
            self.translations[dt][nodeID] = array([t1,t2,t3])
            self.rotations[dt][nodeID]    = array([r1,r2,r3])
        ###

    def parseLength(self):
        self.mainHeaders = []
        self.strFormat = ''
        #print self.dataCode
        if self.analysisCode==5:
            self.mainHeaders.append('Freq')
            self.strFormat += 'fi'
            self.add = self.addF
            #raise Exception('???A')
        #elif self.analysisCode in[6]: # 10
            #self.mainHeaders.append('Time')
            #self.strFormat += 'fi'
            #self.add = self.addF
            #print self.dataCode
            #raise Exception('???B')
        elif self.analysisCode in [1,2,3,4,6,7,8,9,10,11,12]:
            self.mainHeaders.append('NodeID')
            self.strFormat += 'ii'
        else:
            raise Exception('invalid analysisCode=%s' %(self.analysisCode))
        self.mainHeaders.append('GridType')

        if self.isImaginary():  # elif self.dataFormat==1:
            self.strFormat += 'ffffffffffff'
            self.headers = ['Tx','Ty','Tz','Rx','Ry','Rz']
            raise Exception('verify...add imaginary...')
        else:
            self.strFormat += 'ffffff'         # if self.dataFormat in [0,2]:
            self.headers = ['Tx','Ty','Tz','Rx','Ry','Rz']
        
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

    def deleteTransient(self,dt):
        del self.translations[dt]
        del self.rotations[dt]

    def getTransients(self):
        k = self.translations.keys()
        k.sort()
        return k

    def addNewTransient(self):
        """
        initializes the transient variables
        """
        if self.dt not in self.translations:
            self.translations[self.dt] = {}
            self.rotations[self.dt]    = {}

    def addBinary(self,deviceCode,data):
        (nodeID,v1,v2,v3,v4,v5,v6) = unpack('iffffff',data)
        msg = "nodeID=%s v1=%s v2=%s v3=%s" %(nodeID,v1,v2,v3)
        assert 0<nodeID<1000000000, msg
        assert nodeID not in self.translations

        self.translations[nodeID] = array([v1,v2,v3]) # dx,dy,dz
        self.rotations[nodeID]    = array([v4,v5,v6]) # rx,ry,rz
    ###

    def add(self,out):
        (nodeID,gridType,v1,v2,v3,v4,v5,v6) = out
        msg = "nodeID=%s gridType=%s v1=%s v2=%s v3=%s" %(nodeID,gridType,v1,v2,v3)
        #print msg
        nodeID = (nodeID-self.deviceCode) // 10
        assert 0<nodeID<1000000000, msg
        #assert nodeID not in self.displacements,'displacementObject - static failure'
        
        gridType = self.recastGridType(gridType)
        self.gridTypes[nodeID] = gridType
        self.translations[nodeID] = array([v1,v2,v3]) # dx,dy,dz
        self.rotations[nodeID]    = array([v4,v5,v6]) # rx,ry,rz
    ###

    def addF(self,out):
        (freq,gridType,v1,v2,v3,v4,v5,v6) = out
        msg = "dt=%g %s=%s gridType=%s v1=%s v2=%s v3=%s" %(self.dt,self.mainHeaders[0],time,gridType,v1,v2,v3)
        #print msg
        #assert 0<nodeID<1000000000, msg
        #assert nodeID not in self.translations,'displacementObject - static failure'
        
        gridType = self.recastGridType(gridType)
        self.gridTypes[freq]    = gridType
        self.translations[freq] = array([v1,v2,v3]) # dx,dy,dz
        self.rotations[freq]    = array([v4,v5,v6]) # rx,ry,rz
    ###

    def addTransient(self,out):
        (nodeID,gridType,v1,v2,v3,v4,v5,v6) = out
        nodeID = (nodeID-self.deviceCode) // 10
        msg  = "nodeID=%s v1=%s v2=%s v3=%s\n" %(nodeID,v1,v2,v3)
        msg += "          v4=%s v5=%s v6=%s"   %(       v4,v5,v6)
        #print msg
        assert 0<nodeID<1000000000, msg
        #assert nodeID not in self.displacements[self.dt],'displacementObject - transient failure'

        gridType = self.recastGridType(gridType)
        self.gridTypes[nodeID]             = gridType
        self.translations[self.dt][nodeID] = array([v1,v2,v3]) # dx,dy,dz
        self.rotations[self.dt][nodeID]    = array([v4,v5,v6]) # rx,ry,rz
    ###

    def addTransientF(self,out):
        (freq,gridType,v1,v2,v3,v4,v5,v6) = out
        msg = "dt=%g %s=%s gridType=%s v1=%s v2=%s v3=%s" %(self.dt,self.mainHeaders[0],freq,gridType,v1,v2,v3)
        #print msg
        gridType = self.recastGridType(gridType)
        self.gridTypes[freq]             = gridType
        self.translations[self.dt][freq] = array([v1,v2,v3]) # dx,dy,dz
        self.rotations[self.dt][freq]    = array([v4,v5,v6]) # rx,ry,rz
    ###

    def writeOp2(self,block3,deviceCode=1):
        """
        creates the binary data for writing the table
        @warning hasnt been tested...
        """
        msg = block3
        for nodeID,translation in sorted(self.translations.items()):
            rotation = self.rotations[nodeID]
            (dx,dy,dz) = translation
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
    #        #msg += '%s = %g\n' %(self.dataCode['name'],dt)
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

    def writeHeader(self):
        (mainHeaders,headers) = self.getHeaders()
        msg = '%-10s %8s ' %(mainHeaders)
        for header in headers:
            msg += '%10s ' %(header)
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

    def GetAsSort1(self):
        return (self.translations,self.rotations)

    def GetAsSort2(self):
        """returns translations and rotations in sort2 format"""
        translations2 = {}
        rotations2 = {}
        if self.dt is not None:
            return self.__reprTransient__()

            for dt,translations in sorted(self.translations.items()):
                nodeIDs = translations.keys()
                for nodeID in nodeIDs:
                    translations2[nodeID] = {}
                    rotations2[nodeID] = {}

                for nodeID,translation in sorted(translations.items()):
                    rotation = self.rotations[dt][nodeID]
                    translations2[nodeID][dt] = translation
                    rotations2[nodeID][dt]    = rotation
                ###
        else:
            for nodeID,translation in sorted(self.translations.items()):
                rotation = self.rotations[nodeID]
                translations2[nodeID] = translation
                rotations2[nodeID]    = rotation
            ###
        ###
        return translations2,rotations2
        
    def getHeaders(self):
        return (self.mainHeaders,self.headers)

    def writeF06(self,header,pageStamp,pageNum=1):
        msg = ['                                             D I S P L A C E M E N T   V E C T O R',
               ' ',
               '      POINT ID.   TYPE          T1             T2             T3             R1             R2             R3']
        for nodeID,translation in sorted(self.translations.items()):
            rotation = self.rotations[nodeID]
            gridType = self.gridTypes[nodeID]

            (dx,dy,dz) = translation
            (rx,ry,rz) = rotation
            msg.append('%14i %6s     %13E  %13E  %13E  %13E  %13E  %13E' %(nodeID,gridType,dx,dy,dz,rx,ry,rz))
        ###
        msg.append(pageStamp+str(pageNum))
        msg.append('\n')
        return ('\n'.join(msg),pageNum)
            
        
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
