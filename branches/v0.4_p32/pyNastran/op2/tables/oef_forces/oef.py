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
from struct import unpack
from pyNastran.op2.op2Errors import *

# pyNastran
from pyNastran.op2.tables.oug.oug_displacements import displacementObject
#from pyNastran.op2.resultObjects.ougv1_Objects import (
#    temperatureObject,displacementObject,
#    nonlinearTemperatureObject,
#    fluxObject,nonlinearFluxObject)

#from pyNastran.op2.tables.oug.ougv1_Objects import (
#    displacementObject,temperatureObject)
#from pyNastran.op2.tables.oug.oug_eigenvectors import (
#    eigenVectorObject)
from .oef_Objects import (nonlinearFluxObject)

class OEF(object):
    """Table of element forces"""
    def readTable_OEF(self):
        table3     = self.readTable_OEF_3
        table4Data = self.readOEF1_Data
        self.readResultsTable(table3,table4Data)
        self.deleteAttributes_OEF()
        

    def deleteAttributes_OEF(self):
        params = ['elementType','dLoadID','loadID','obj','markerStart','oCode',
                  'eigr','eigi','eign','mode','freq','time','thermal',]
        self.deleteAttributes(params)
        #print self.obj
        #sys.exit('stopping in oef.py')

    def readTable_OEF_3(self,iTable): # iTable=-3
        bufferWords = self.getMarker()
        #print "2-bufferWords = ",bufferWords,bufferWords*4,'\n'

        data = self.getData(4)
        bufferSize, = unpack('i',data)
        data = self.getData(4*50)
        #self.printBlock(data)
        
        aCode = self.getBlockIntEntry(data,1)
        
        self.parseApproachCode(data)
        self.addDataParameter(data,'elementType', 'i',3,False)   ## element type
        self.addDataParameter(data,'dLoadID',     'i',8,False)   ## dynamic load set ID/random code
        self.addDataParameter(data,'formatCode',  'i',9,False)   ## format code
        self.addDataParameter(data,'numWide',     'i',10,False)  ## number of words per entry in record; @note is this needed for this table ???
        self.addDataParameter(data,'oCode',       'i',11,False)  ## undefined in DMAP...
        self.addDataParameter(data,'thermal',     'i',23,False)  ## thermal flag; 1 for heat ransfer, 0 otherwise

        #print "dLoadID(8)=%s formatCode(9)=%s numwde(10)=%s oCode(11)=%s thermal(23)=%s" %(self.dLoadID,self.formatCode,self.numWide,self.oCode,self.thermal)
        #print "thermal(23)=%s elementType(3)=%s" %(self.thermal,self.elementType)


        ## assuming tCode=1
        if self.analysisCode==1:   # statics
            self.addDataParameter(data,'loadID','i',5,False)   ## load set ID number
        elif self.analysisCode==2: # normal modes/buckling (real eigenvalues)
            self.addDataParameter(data,'mode','i',5)   ## mode number
            self.addDataParameter(data,'eign','f',6,False)   ## eigenvalue
        elif self.analysisCode==3: # differential stiffness 0
            self.addDataParameter(data,'loadID','i',5)   ## load set ID number
        elif self.analysisCode==4: # differential stiffness 1
            self.addDataParameter(data,'loadID','i',5)   ## load set ID number
        elif self.analysisCode==5:   # frequency
            self.addDataParameter(data,'freq','f',5)   ## frequency

        elif self.analysisCode==6: # transient
            self.addDataParameter(data,'time','f',5)   ## time step
            #print "time(5)=%s" %(self.time)
        elif self.analysisCode==7: # pre-buckling
            self.addDataParameter(data,'loadID','i',5)   ## load set ID number
            #print "loadID(5)=%s" %(self.loadID)
        elif self.analysisCode==8: # post-buckling
            self.addDataParameter(data,'loadID','i',5)       ## load set ID number
            self.addDataParameter(data,'eigr','f',6,False)   ## real eigenvalue
            #print "loadID(5)=%s  eigr(6)=%s" %(self.loadID,self.eigr)
        elif self.analysisCode==9: # complex eigenvalues
            self.addDataParameter(data,'mode','i',5)         ## mode number
            self.addDataParameter(data,'eigr','f',6,False)   ## real eigenvalue
            self.addDataParameter(data,'eigi','f',7,False)   ## imaginary eigenvalue
            #print "mode(5)=%s  eigr(6)=%s  eigi(7)=%s" %(self.mode,self.eigr,self.eigi)
        elif self.analysisCode==10: # nonlinear statics
            self.addDataParameter(data,'loadStep','f',5)   ## load step
            #print "loadStep(5) = %s" %(self.loadStep)
        elif self.analysisCode==11: # geometric nonlinear statics
            self.addDataParameter(data,'loadID','i',5)   ## load set ID number
            #print "loadID(5)=%s" %(self.loadID)
        else:
            raise InvalidAnalysisCodeError('invalid analysisCode...analysisCode=%s' %(str(self.analysisCode)+'\n'+self.codeInformation()))

        # tCode=2
        #if self.analysisCode==2: # sort2
        #    self.loadID = self.getValues(data,'i',5) ## load set ID number
        
        #print "*iSubcase=%s"%(self.iSubcase)
        #print "analysisCode=%s tableCode=%s thermal=%s" %(self.analysisCode,self.tableCode,self.thermal)
        #print self.codeInformation()

        #self.printBlock(data)
        self.readTitle()

    def readOEF1_Data(self):
        tfsCode = [self.tableCode,self.formatCode,self.sortCode]
        self.skipOES_Element() # skipping entire table
        self.skippedCardsFile.write('skipping atfsCode=%s in %s-OEF\n' %(self.atfsCode,self.tableName))
        return

        #print "tfsCode = %s" %(tfsCode)
        # element forces & moments / flux
        if 1:
            if tfsCode==[4,1,0]:
                self.readOEF1_Data_format1_sort0()
            elif tfsCode==[4,1,1]:
                self.readOEF1_Data_format1_sort1()
            elif tfsCode==[4,2,1]:
                self.readOEF1_Data_format2_sort1()
            elif tfsCode==[4,3,1]:
                self.readOEF1_Data_format3_sort1()

            # composite failure indicies
            elif tfsCode==[25,1,0]:
                self.readOEF1_Data_format1_sort0()
            else:
                raise Exception('bad tableCode/formatCode/sortCode=%s on OEF table' %(tfsCode))
            ###
        else:
            self.skipOES_Element()
        ###

    def readOEF1_Data_format3_sort1(self):
        self.skipOES_Element()

    def readOEF1_Data_format1_sort0(self):
        assert self.formatCode==1
        assert self.sortCode==0

        self.log.debug("self.analysisCode=%s tableCode(1)=%s sortCode=%s thermal(23)=%g" %(self.analysisCode,self.tableCode,self.sortCode,self.thermal))

        if self.thermal==0:
            if self.analysisCode==1: # displacement
                #print "isForces"
                self.createTransientObject(self.forces,displacementObject)
            elif self.analysisCode==2 and self.sortCode==1: # buckling forces
                #print "isBucklingForces"
                self.createTransientObject(self.forces,displacementObject)
            elif self.analysisCode==5: # frequency forces
                #print "isFrequencyForces"
                self.createTransientObject(self.modalSPCForces,eigenVectorObject)
            elif self.analysisCode==6: # transient displacement
                #print "isTransientForces"
                self.createTransientObject(self.forces,displacementObject)
            elif self.analysisCode==9: # complex eigenvalue forces
                #print "isComplexEigenvalues"
                self.createTransientObject(self.modalSPCForces,eigenVectorObject)
            elif self.analysisCode==10: # nonlinear static displacement
                #print "isNonlinearStaticForces"
                self.createTransientObject(self.forces,displacementObject)
            else:
                pass
                raise NotImplementedError('not supported OEF static solution...')
            ###

        elif self.thermal==1:
            #if self.analysisCode==1: # temperature
                #print "isTemperature"
                #raise Exception('verify...')
                #self.temperatures[self.iSubcase] = temperatureObject(self.iSubcase)
            #elif self.analysisCode==1 and self.sortCode==1: # heat fluxes
                #print "isFluxes"
                #raise Exception('verify...')
                #self.createTransientObject(self.fluxes,fluxObject)
            if self.analysisCode==5: # frequency forces
                #print "isFrequencyForces"
                self.createTransientObject(self.modalSPCForces,eigenVectorObject)
            elif self.analysisCode==6: # transient temperature
                #print "isTransientTemperature"
                self.createTransientObject(self.temperatureForces,temperatureObject)
            elif self.analysisCode==10: # nonlinear static displacement
                #print "isNonlinearStaticTemperatures"
                self.createTransientObject(self.fluxes,nonlinearFluxObject)
            else:
                msg = 'OEF_thermal format1_sort0 elementType=%-3s -> %-6s is not supported - fname=%s\n' %(self.elementType,self.ElementType(self.elementType),self.op2FileName)
                self.skippedCardsFile.write(msg)
                raise NotImplementedError('not supported OEF thermal solution...')
            ###
        else:
            msg = 'invalid thermal flag...not 0 or 1...flag=%s\n' %(self.thermal)
            sys.stderr.write(msg)
            raise Exception(msg)
        ###
        if self.obj:
            #self.skipOES_Element()
            self.readScalars8(debug=True)
            #self.readForces(data,self.obj)
            #self.skipOES_Element()
        else:
            self.skipOES_Element()
        ###

    def readOEF1_Data_format1_sort1(self):
        assert self.formatCode==1
        assert self.sortCode==1

        #print "self.analysisCode=%s tableCode(1)=%s thermal(23)=%g" %(self.analysisCode,self.tableCode,self.thermal)
        if self.thermal==0:
            if self.analysisCode==5: # frequency forces
                #print "isFrequencyForces"
                self.createTransientObject(self.modalSPCForces,eigenVectorObject)
                #self.readForces(self.obj)
            elif self.analysisCode==9: # complex eigenvalue forces
                #print "isComplexEigenvalues"
                self.createTransientObject(self.modalSPCForces,eigenVectorObject)
            else:
                raise NotImplementedError('not supported OEF static solution...')
            ###
        else:
            raise NotImplementedError('invalid thermal flag...not 0 or 1...flag=%s' %(self.thermal))
        ###
        self.skipOES_Element()
        #self.readForces(data,self.obj)
        

    def readOEF1_Data_format2_sort1(self):
        assert self.formatCode==2
        assert self.sortCode==1

        #print "self.analysisCode=%s tableCode(1)=%s thermal(23)=%g" %(self.analysisCode,self.tableCode,self.thermal)
        if self.thermal==0:
            if self.analysisCode==5: # frequency forces
                #print "isFrequencyForces"
                self.createTransientObject(self.modalSPCForces,eigenVectorObject)
            elif self.analysisCode==9: # complex eigenvalue forces
                #print "isComplexEigenvalues"
                self.createTransientObject(self.modalSPCForces,eigenVectorObject)
            else:
                raise NotImplementedError('not supported OEF static solution...')
            ###

        else:
            raise NotImplementedError('invalid thermal flag...not 0 or 1...flag=%s' %(self.thermal))
        ###
        if self.obj:
            self.skipOES_Element()
        else:
            self.skipOES_Element()
        #self.readForces(data,self.obj)

        
    def setupOEF(self,elementType,elementName,numWide,deviceCode):
        """length of result"""
        #print "elementType=%s numWide=%s type=%s" %(elementType,numWide,elementName)
        
        fCode = self.formatCode
        eType = self.elementType
        if self.thermal==0:
            if eType==1: # Rod
                if fCode in [0,2]:  recordLength = 12   # 3 fields
                else:               recordLength = 20   # 4 fields
            elif eType==2: # Beam
                if fCode in [0,2]:  recordLength = 40   # 10 fields
                else:               recordLength = 76   # 17 fields
            elif eType==10: # Conrod
                if fCode in [0,2]:  recordLength = 12   # 3 fields
                else:               recordLength = 16   # 4 fields
            elif eType==33: # CQUAD4
                if fCode in [0,2]:  recordLength = 36   # 9  fields
                else:               recordLength = 76   # 17 fields
            elif eType==34: # BAR
                if fCode in [0,2]:  recordLength = 36   # 9  fields
                else:               recordLength = 76   # 17 fields
            elif eType==74: # TRIA3
                if fCode in [0,2]:  recordLength = 36   # 9  fields
                else:               recordLength = 76   # 17 fields
            elif eType==100: # Bars
                if fCode in [0,2]:  recordLength = 24   # 8  fields
                else:               recordLength = 56   # 14 fields
            else:
                recordLength = None
                #raise Exception('need to define the word size for static elementType=%s=%s' %(self.elementType,elementName))
            ###
        else:
            recordLength = None
            #raise Exception('need to define the word size for thermal elementType=%s=%s' %(self.elementType,elementName))

        isSkipped = True
        func = self.skipMe
        ###
        #if numWide in [9,10]:
        #    func = self.readOEF_2D_3D
        #    isSkipped = False
        #    recordLength = 40
        #elif numWide==8:
        #    func = self.readOEF_CHBDY
        #    isSkipped = True
        #    recordLength = 40
        #else:
        #    raise Exception('need to define the word size for elementType=%s=%s' %(self.elementType,elementName))
        ###
        self.deviceCode = deviceCode
        return(recordLength,func,isSkipped) # 8*4

    def skipMe(self,data):
        pass

    def readOEF_2D_3D(self,data):
        #print "read_2D_3D"
        gridDevice, = unpack('i',data[0:4])
        grid = (gridDevice-self.deviceCode)/10
        eType = ''.join(unpack('cccccccc',data[4:12]))
        (xGrad,yGrad,zGrad,xFlux,yFlux,zFlux) = unpack('ffffff',data[12:36])
        #print "grid=%g dx=%i dy=%i dz=%i rx=%i ry=%i rz=%i" %(grid,xGrad,yGrad,zGrad,xFlux,yFlux,zFlux)
        
        return(grid,eType,xGrad,yGrad,zGrad,xFlux,yFlux,zFlux)

    def readOEF_CHBDY(self,data):
        #print "read_CHBDYx"
        gridDevice, = unpack('i',data[0:4])
        grid = (gridDevice-self.deviceCode)/10
        eType = ''.join(unpack('cccccccc',data[4:12]))
        (fApplied,freeConv,forceConv,fRad,fTotal) = unpack('fffff',data[12:32])
        return(grid,eType,fApplied,freeConv,forceConv,fRad,fTotal)

    def readForces(self,scalarObject):
        #print "readForces..."
        #print type(scalarObject)
        data = self.data
        #self.printBlock(data[0:self.numWide*4])
        
        (reqLen,func,isSkipped) = self.setupOEF(self.elementType,self.ElementType(self.elementType),self.numWide,self.deviceCode)
        if reqLen==None:
            reqLen = 1
            data = ''

        while len(data)>=reqLen:
            eData = data[:reqLen]
            #print "len(data) = ",len(data)
            #self.printBlock(data[:self.numWide*4])
            #print "eType = ",eType
            #print "len(data[8:40]"

            if not isSkipped:
                out = func(eData)
                scalarObject.add(*out)
            #print "gridDevice = ",gridDevice
            #print "deviceCode = ",deviceCode
            #grid = (gridDevice-self.deviceCode)/10
            #print "grid=%g dx=%g dy=%g dz=%g rx=%g ry=%g rz=%g" %(grid,xGrad,yGrad,zGrad,xFlux,yFlux,zFlux)
            #print type(scalarObject)
            #scalarObject.add(grid,xGrad,yGrad,zGrad,xFlux,yFlux,zFlux)
            data = data[reqLen:]
        ###
        #print self.obj
        self.handleResultsBuffer(self.readForces,scalarObject,debug=False)


    def passFunc(self,data):
        return

    #def readOEF_2D_3D(self,data):
    #    gridDevice, = unpack('i',data[0:4])
    #    grid = (gridDevice-self.deviceCode)/10
    #    eType = ''.join(unpack('cccccccc',data[4:12]))
    #    (xGrad,yGrad,zGrad,xFlux,yFlux,zFlux) = unpack('ffffff',data[12:36])
    #    return(grid,eType,xGrad,yGrad,zGrad,xFlux,yFlux,zFlux)

    def getOEF_nWords(self):
        if self.thermal==0:
            if self.elementType==12: # CELAS2
                if self.tableCode in [0,2]:  nWords = 2
                else:                        nWords = 3
            ###
            else:
                raise Exception('need to define the word size for elementType=%s=%s' %(self.elementType,self.ElementType(self.elementType)))
            ###
        else: # thermal
            nWords = self.numWide
        ###
        return nWords

    def readForcesNonlinear(self,scalarObject):
        #print "readForcesNonlinear..."
        data = self.data

        #print 'thermal skipping elementType=%s=%s' %(self.elementType,self.ElementType(self.elementType))
        sys.stderr.write('thermal skipping elementType=%s=%s\n' %(self.elementType,self.ElementType(self.elementType)))
        
        nWords = self.getOEF_nWords()
        reqLen = 4*nWords

        while len(data)>=reqLen:
            #print "len(data) = ",len(data),reqLen
            #self.printBlock(data[32:])
            gridDevice, = unpack('i',data[0:4])
            #eType = ''.join(unpack('cccccccc',data[4:12]))
            #print "eType = ",eType
            #print "len(data[8:40]"
            #print "elementType=%s" %(self.elementType)

            if self.elementType==12:
                if self.tableCode in [0,2]:
                    force = unpack('f',data[4:8])
                else:
                    (forceReal,forceImag) = unpack('ff',data[4:12])
                ###
            ###
            else:
                pass
                #raise Exception('elementType=%s' %(self.elementType))
            ###

            #if self.numWide in [9,10]:
            #    (xGrad,yGrad,zGrad,xFlux,yFlux,zFlux) = unpack('ffffff',data[12:36])
            #elif self.numWide==2: ## @todo CHBDY - how do i add this to the case...
            #    (fApplied,freeConv,forceConv,fRad,fTotal) = unpack('fffff',data[12:32])
            #    sys.stderr.write('skipping CHBDY\n')
            #    data = data[self.numWide*4:]
            #    continue
            #elif self.numWide==8: ## @todo CHBDY - how do i add this to the case...
            #    (fApplied,freeConv,forceConv,fRad,fTotal) = unpack('fffff',data[12:32])
            #    sys.stderr.write('skipping CHBDY\n')
            #    data = data[self.numWide*4:]
            #    continue
            #else:
            #    raise Exception('only CBEAM/CBAR/CTUBE/2D/3D elements supported...so no special thermal elements...numwde=%s' %(self.numWide))
            
            #print "gridDevice = ",gridDevice
            #print "deviceCode = ",deviceCode
            grid = (gridDevice-self.deviceCode)/10
            #print "grid=%g dx=%i dy=%i dz=%i rx=%i ry=%i rz=%i" %(grid,xGrad,yGrad,zGrad,xFlux,yFlux,zFlux)
            #print type(scalarObject)
            #scalarObject.add(grid,eType,xGrad,yGrad,zGrad,xFlux,yFlux,zFlux)
            data = data[self.numWide*4:]
        ###
        #print self.obj
        self.handleResultsBuffer(self.readForcesNonlinear,scalarObject,debug=False)
        
