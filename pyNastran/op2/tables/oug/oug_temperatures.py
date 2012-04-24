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

class temperatureObject(scalarObject): # approachCode=1, sortCode=0, thermal=1
    def __init__(self,dataCode,iSubcase,dt=None):
        scalarObject.__init__(self,dataCode,iSubcase)
        self.dt = dt
        self.gridTypes    = {}
        self.temperatures = {}
        
        if dt is not None:
            #raise Exception('potential test of temperature bug...')
            self.addNewTransient()
            #assert dt>=0.
            self.isTransient = True
            self.add = self.addTransient
            self.addF = self.addTransientF
            #self.addBinary = self.addBinaryTransient
            #self.__repr__ = self.__reprTransient__  # why cant i do this...
        ###
        self.parseLength()

    def addF06Data(self,data,transient):
        if transient is None:
            for line in data:
                (gridID,gridType) = line[0:2]
                temps = line[2:]
                for (i,temp) in enumerate(temps):
                    nodeID = gridID+i
                    self.gridTypes[nodeID] = gridType
                    self.temperatures[nodeID] = temp
                ###
            ###
            return

        (dtName,dt) = transient
        self.dataCode['name'] = dtName
        if dt not in self.temperatures:
            self.updateDt(self.dataCode,dt)
            self.isTransient = True

        for line in data:
            (gridID,gridType) = line[0:2]
            temps = line[2:]
            for (i,temp) in enumerate(temps):
                nodeID = gridID+i
                self.gridTypes[nodeID] = gridType
                self.temperatures[dt][nodeID] = temp
            ###
        ###

    def parseLength(self):
        self.mainHeaders = []
        self.strFormat = ''
        print "analysisCode = %s" %(self.analysisCode)
        if 0: #self.analysisCode==6:  # disabled
            self.mainHeaders.append('Freq')
            self.strFormat += 'fi'
            self.add = self.addF
        elif self.analysisCode in[5]:
            self.mainHeaders.append('Time')
            self.strFormat += 'fi'
            self.add = self.addF
        elif self.analysisCode in [1,2,3,4,6,7,8,9,10,11,12]:
            self.mainHeaders.append('NodeID')
            self.strFormat += 'ii'
        else:
            raise NotImplementedError('invalid analysisCode=%s' %(self.analysisCode))
        self.mainHeaders.append('GridType')

        if self.isImaginary():  # elif self.dataFormat==1:
            self.strFormat += 'ffffffffffff'
            self.headers = ['Temperature']
            raise NotImplementedError('verify...add imaginary...')
        else:
            self.strFormat += 'ffffff'         # if self.dataFormat in [0,2]
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

    def deleteTransient(self,dt):
        del self.temperatures[dt]

    def getTransients(self):
        k = self.temperatures.keys()
        k.sort()
        return k
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

    def writeF06(self,header,pageStamp,pageNum=1):
        #header = ['     DEFAULT                                                                                                                        ',
        #           '0                                                                                                            SUBCASE 1              ',
        #           '     LOAD STEP =  1.00000E+00']

        msgOrig = ['                                              T E M P E R A T U R E   V E C T O R\n',
               ' \n',
               '      POINT ID.   TYPE      ID   VALUE     ID+1 VALUE     ID+2 VALUE     ID+3 VALUE     ID+4 VALUE     ID+5 VALUE\n']
        msg = []
        if self.isTransient:
            for dt,temperatures in sorted(self.temperatures.items()):
                dtLine = '%14s = %12.5E\n'%(self.dataCode['name'],dt)
                header[2] = dtLine
                msg += header+msgOrig
                msg += self.printTempLines(temperatures)
                msg.append(pageStamp+str(pageNum)+'\n')
                pageNum += 1
            ###
            return(''.join(msg),pageNum) # transient

        msg += self.printTempLines(self.temperatures)
        msg.append(pageStamp+str(pageNum)+'\n')
        return(''.join(msg),pageNum)  # static
    
    def printTempLines(self,temperatures):
        msg = []
        pack = []
        oldNodeID = -1
        oldGridType = None
        for nodeID,T in sorted(temperatures.items()):
            gridType = self.gridTypes[nodeID]

            if oldNodeID+1==nodeID and gridType==oldGridType:
                oldNodeID = nodeID
                pack.append(T)
            else:
                if oldNodeID>0:
                    msg += self.printPack(pack)
                oldGridType = gridType
                oldNodeID = nodeID
                pack = [nodeID,gridType,T]
            ###
        ###
        if pack:
            msg += self.printPack(pack)
        ###
        return msg

    def printPack(self,pack):
        msg = []
        nID   = pack[0]
        gType = pack[1]
        while len(pack)>8:
            nID = pack[0]
            packOut = pack[:8]
            pack = [nID+6,gType]+pack[8:]
            msg.append('      %8i   %4s      %10.6E   %10.6E   %10.6E   %10.6E   %10.6E   %10.6E\n' %(tuple(packOut)))
        ###
        if pack:
            fmt = '      %8i   %4s   '+'   %10.6E'*(len(pack)-2)+'\n'
            out = fmt %(tuple(pack))
            msg.append(out)
        ###
        return msg

    def __repr__(self):
        if self.isTransient:
            return self.__reprTransient__()

        msg = '---TEMPERATURE---\n'
        msg += self.writeHeader()
        #print "self.dataCode=",self.dataCode
        for nodeID,T in sorted(self.temperatures.items()):
            gridType = self.gridTypes[nodeID]
            msg += '%10s %8s ' %(nodeID,gridType)
            #print "nodeID=%s T=%s" %(nodeID,T)
            if abs(T)<1e-6:
                msg += '%10s\n' %(0)
            else:
                msg += '%10g\n' %(T)
            ###
        return msg

