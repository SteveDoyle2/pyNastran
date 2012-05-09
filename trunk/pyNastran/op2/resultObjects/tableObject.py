import sys
from struct import pack
from numpy import array
from op2_Objects import scalarObject

class TableObject(scalarObject):  # displacement style table
    def __init__(self,dataCode,isSort1,iSubcase,dt):
        scalarObject.__init__(self,dataCode,iSubcase)
        self.dt = dt
        self.gridTypes    = {}
        self.translations = {}
        self.rotations    = {}

        if isSort1:
            if dt is not None:
                self.add = self.addSort1
            ###
        else:
            assert dt is not None
            self.add = self.addSort2
        ###

        #if dt is not None:
            #self.addNewTransient()
            #self.add  = self.addTransient
            #self.addF = self.addTransientF
            #self.addBinary = self.addBinaryTransient
            #self.__repr__ = self.__reprTransient__  # why cant i do this...
            #self.writeOp2 = self.writeOp2Transient
        ###
        #self.parseLength()

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

    def addNewTransient(self,dt):
        """initializes the transient variables"""
        self.translations[dt] = {}
        self.rotations[dt]    = {}

    #def addBinary(self,deviceCode,data):
        #(nodeID,v1,v2,v3,v4,v5,v6) = unpack('iffffff',data)
        #msg = "nodeID=%s v1=%s v2=%s v3=%s" %(nodeID,v1,v2,v3)
        #assert 0<nodeID<1000000000, msg
        #assert nodeID not in self.translations

        #self.translations[nodeID] = array([v1,v2,v3]) # dx,dy,dz
        #self.rotations[nodeID]    = array([v4,v5,v6]) # rx,ry,rz
    ###

    def add(self,dt,out):
        (nodeID,gridType,v1,v2,v3,v4,v5,v6) = out
        msg = "nodeID=%s gridType=%s v1=%s v2=%s v3=%s" %(nodeID,gridType,v1,v2,v3)
        #print msg
        assert 0<nodeID<1000000000, msg
        #assert nodeID not in self.displacements,'displacementObject - static failure'
        
        gridType = self.recastGridType(gridType)
        self.gridTypes[nodeID] = gridType
        self.translations[nodeID] = array([v1,v2,v3]) # dx,dy,dz
        self.rotations[nodeID]    = array([v4,v5,v6]) # rx,ry,rz
    ###

    def addSort1(self,dt,out):
        (nodeID,gridType,v1,v2,v3,v4,v5,v6) = out
        if dt not in self.translations:
            self.addNewTransient(dt)
        msg  = "nodeID=%s v1=%s v2=%s v3=%s\n" %(nodeID,v1,v2,v3)
        msg += "          v4=%s v5=%s v6=%s"   %(       v4,v5,v6)
        #print msg
        #assert 0<nodeID<1000000000, msg
        #assert nodeID not in self.translations[self.dt],'displacementObject - transient failure'

        gridType = self.recastGridType(gridType)
        self.gridTypes[nodeID]             = gridType
        self.translations[dt][nodeID] = array([v1,v2,v3]) # dx,dy,dz
        self.rotations[dt][nodeID]    = array([v4,v5,v6]) # rx,ry,rz
    ###

    def addSort2(self,nodeID,out):
        (dt,gridType,v1,v2,v3,v4,v5,v6) = out
        if dt not in self.translations:
            self.addNewTransient(dt)
        msg  = "nodeID=%s v1=%s v2=%s v3=%s\n" %(nodeID,v1,v2,v3)
        msg += "          v4=%s v5=%s v6=%s"   %(       v4,v5,v6)
        #print msg
        #assert 0<nodeID<1000000000, msg
        #assert nodeID not in self.translations[self.dt],'displacementObject - transient failure'

        gridType = self.recastGridType(gridType)
        self.gridTypes[nodeID]             = gridType
        self.translations[dt][nodeID] = array([v1,v2,v3]) # dx,dy,dz
        self.rotations[dt][nodeID]    = array([v4,v5,v6]) # rx,ry,rz
    ###

    #def writeOp2(self,block3,deviceCode=1):
        #"""
        #creates the binary data for writing the table
        #@warning hasnt been tested...
        #"""
        #msg = block3
        #for nodeID,translation in sorted(self.translations.iteritems()):
            #rotation = self.rotations[nodeID]
            #(dx,dy,dz) = translation
            #(rx,ry,rz) = rotation
            #grid = nodeID*10+deviceCode
            #msg += pack('iffffff',grid,dx,dy,dz,rx,ry,rz)
        ###
        #return msg

    #def writeOp2Transient(self,block3,deviceCode=1):
    #    """
    #    creates the binary data for writing the table
    #    @warning hasnt been tested...
    #    @warning dt slot needs to be fixed...
    #    """
    #    msg = ''
    #    for dt,displacements in sorted(self.displacements.iteritems()):
    #        XXX = 50 ## this isnt correct... @todo update dt
    #        msg += block3[0:XXX] + pack('i',dt) + block3[XXX+4:]
    #        #msg += '%s = %g\n' %(self.dataCode['name'],dt)
    #
    #        for nodeID,displacement in sorted(displacements.iteritems()):
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
        #(mainHeaders,headers) = self.getHeaders()
        mainHeaders = ('nodeID','gridType')
        headers = ('T1','T2','T3','R1','R2','R3')
        msg = '%-10s %8s ' %(mainHeaders)
        for header in headers:
            msg += '%10s ' %(header)
        msg += '\n'
        return msg

    def GetAsSort1(self):
        return (self.translations,self.rotations)

    def GetAsSort2(self):
        """returns translations and rotations in sort2 format"""
        translations2 = {}
        rotations2 = {}
        if self.dt is not None:
            return self.__reprTransient__()

            for dt,translations in sorted(self.translations.iteritems()):
                nodeIDs = translations.keys()
                for nodeID in nodeIDs:
                    translations2[nodeID] = {}
                    rotations2[nodeID] = {}

                for nodeID,translation in sorted(translations.iteritems()):
                    rotation = self.rotations[dt][nodeID]
                    translations2[nodeID][dt] = translation
                    rotations2[nodeID][dt]    = rotation
                ###
        else:
            for nodeID,translation in sorted(self.translations.iteritems()):
                rotation = self.rotations[nodeID]
                translations2[nodeID] = translation
                rotations2[nodeID]    = rotation
            ###
        ###
        return translations2,rotations2
        
    def getHeaders(self):
        return (self.mainHeaders,self.headers)

class complexTableObject(scalarObject):
    def __init__(self,dataCode,isSort1,iSubcase,dt):
        scalarObject.__init__(self,dataCode,iSubcase)
        self.dt = dt
        self.gridTypes    = {}
        self.translations = {}
        self.rotations    = {}

        if isSort1:
            if dt is not None:
                self.add = self.addSort1
            ###
        else:
            assert dt is not None
            self.add = self.addSort2
        ###

        #if dt is not None:
            #self.addNewTransient()
            #self.add  = self.addTransient
            #self.addF = self.addTransientF
            #self.addBinary = self.addBinaryTransient
            #self.__repr__ = self.__reprTransient__  # why cant i do this...
            #self.writeOp2 = self.writeOp2Transient
        ###
        
    def addF06Data(self,data,transient):
        if transient is None:
            for line in data:
                (nodeID,gridType,v1r,v2r,v3r,v4r,v5r,v6r, v1i,v2i,v3i,v4i,v5i,v6i) = line
                self.gridTypes[nodeID]    = gridType
                self.translations[self.dt][nodeID] = [complex(v1r,v1i),complex(v2r,v2i),complex(v3r,v3i)] # dx,dy,dz
                self.rotations[self.dt][nodeID]    = [complex(v4r,v4i),complex(v5r,v5i),complex(v6r,v6i)] # rx,ry,rz
            ###
            return

        (dtName,dt) = transient
        self.dataCode['name'] = dtName
        if dt not in self.translations:
            self.updateDt(self.dataCode,dt)

        for line in data:
            (nodeID,gridType,v1r,v2r,v3r,v4r,v5r,v6r, v1i,v2i,v3i,v4i,v5i,v6i) = line
            self.gridTypes[nodeID]    = gridType
            self.translations[self.dt][nodeID] = [complex(v1r,v1i),complex(v2r,v2i),complex(v3r,v3i)] # dx,dy,dz
            self.rotations[self.dt][nodeID]    = [complex(v4r,v4i),complex(v5r,v5i),complex(v6r,v6i)] # rx,ry,rz
        ###

    def updateDt(self,dataCode,dt):
        self.dataCode = dataCode
        self.applyDataCode()
        if dt is not None:
            self.log.debug("updating %s...%s=%s  iSubcase=%s" %(self.dataCode['name'],self.dataCode['name'],dt,self.iSubcase))
            self.addNewTransient()
        ###

    def deleteTransient(self,dt):
        del self.translations[dt]
        del self.rotations[dt]

    def getTransients(self):
        k = self.translations.keys()
        k.sort()
        return k

    def addNewTransient(self,dt):
        """initializes the transient variables"""
        self.translations[dt] = {}
        self.rotations[dt]    = {}

    def addSort2(self,eid,data):
        [dt,gridType,v1r,v1i,v2r,v2i,v3r,v3i,v4r,v4i,v5r,v5i,v6r,v6i] = data

        if dt not in self.translations:
            self.addNewTransient(dt)

        assert isinstance(eid,int),eid
        assert 0<nodeID<1000000000, msg
        if gridType==1:
            Type = 'G'
        elif gridType==2:
            Type = 'S'
        elif gridType==7:
            Type = 'L'
        else:
            raise Exception('invalid grid type...gridType=%s' %(gridType))

        self.gridTypes[nodeID] = Type
        self.translations[dt][nodeID] = [complex(v1r,v1i),complex(v2r,v2i),complex(v3r,v3i)] # dx,dy,dz
        self.rotations[dt][nodeID]    = [complex(v4r,v4i),complex(v5r,v5i),complex(v6r,v6i)] # rx,ry,rz
        
    def addSort1(self,dt,out):
        (nodeID,gridType,v1r,v1i,v2r,v2i,v3r,v3i,v4r,v4i,v5r,v5i,v6r,v6i) = out
        msg = "dt=%s nodeID=%s v1r=%s v2r=%s v3r=%s" %(dt,nodeID,v1r,v2r,v3r)
        #print msg
        if dt not in self.translations:
            self.addNewTransient(dt)
        #print msg
        #msg = ''
        #assert isinstance(nodeID,int),nodeID
        assert 0<nodeID<1000000000, msg
        #assert nodeID not in self.translations,'complexDisplacementObject - static failure'

        #if gridType==0:
        #    Type = 'S'
        if gridType==1:
            Type = 'G'
        elif gridType==2:
            Type = 'S'
        elif gridType==7:
            Type = 'L'
        else:
            raise Exception('invalid grid type...gridType=%s' %(gridType))

        self.gridTypes[nodeID] = Type
        self.translations[dt][nodeID] = [complex(v1r,v1i),complex(v2r,v2i),complex(v3r,v3i)] # dx,dy,dz
        self.rotations[dt][nodeID]    = [complex(v4r,v4i),complex(v5r,v5i),complex(v6r,v6i)] # rx,ry,rz

    def add(self,dt,out):
        (nodeID,gridType,v1r,v1i,v2r,v2i,v3r,v3i,v4r,v4i,v5r,v5i,v6r,v6i) = out
        msg = "dt=%s nodeID=%s v1r=%s v2r=%s v3r=%s" %(dt,nodeID,v1r,v2r,v3r)
        #assert isinstance(nodeID,int),nodeID
        assert 0<nodeID<1000000000, msg
        #assert nodeID not in self.translations,'complexDisplacementObject - static failure'

        if gridType==1:
            Type = 'G'
        elif gridType==2:
            Type = 'S'
        elif gridType==7:
            Type = 'L'
        else:
            raise Exception('invalid grid type...gridType=%s' %(gridType))

        self.gridTypes[nodeID] = Type
        self.translations[nodeID] = [complex(v1r,v1i),complex(v2r,v2i),complex(v3r,v3i)] # dx,dy,dz
        self.rotations[nodeID]    = [complex(v4r,v4i),complex(v5r,v5i),complex(v6r,v6i)] # rx,ry,rz
    ###
