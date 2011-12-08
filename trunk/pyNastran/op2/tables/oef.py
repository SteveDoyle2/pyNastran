import sys
from struct import unpack

# pyNastran
#from pyNastran.op2.resultObjects.ougv1_Objects import (
#    temperatureObject,displacementObject,
#    nonlinearTemperatureObject,
#    fluxObject,nonlinearFluxObject)

from pyNastran.op2.resultObjects.ougv1_Objects import (
    displacementObject,temperatureObject,
    eigenVectorObject)
from pyNastran.op2.resultObjects.oef_Objects import (
    nonlinearFluxObject)

class OEF(object):
    """Table of element forces"""
    def readTable_OEF1(self):
        table3     = self.readTable_OEF_3
        table4Data = self.readTable_OEF_4_Data
        self.readResultsTable(table3,table4Data)
        self.deleteAttributes_OEF()

    def deleteAttributes_OEF(self):
        params = ['elementType','dLoadID','loadID','obj','markerStart','oCode',
                  'eigr','eigi','eign','mode','freq','time','thermal',]
        self.deleteAttributes(params)

    def readTable_OEF_3(self,iTable): # iTable=-3
        bufferWords = self.getMarker()
        #print "2-bufferWords = ",bufferWords,bufferWords*4,'\n'

        data = self.getData(4)
        bufferSize, = unpack('i',data)
        data = self.getData(4*50)
        #self.printBlock(data)
        
        aCode = self.getBlockIntEntry(data,1)
        print "aCode = ",aCode
        
        self.parseApproachCode(data)
        self.elementType = self.getValues(data,'i',3)  ## element type
        self.dLoadID     = self.getValues(data,'i',8)  ## dynamic load set ID/random code
        self.formatCode  = self.getValues(data,'i',9)  ## format code
        self.numWide     = self.getValues(data,'i',10) ## number of words per entry in record; @note is this needed for this table ???
        self.oCode       = self.getValues(data,'i',11) ## undefined in DMAP...
        self.thermal     = self.getValues(data,'i',23) ## thermal flag; 1 for heat ransfer, 0 otherwise
        print "dLoadID(8)=%s formatCode(9)=%s numwde(10)=%s oCode(11)=%s thermal(23)=%s" %(self.dLoadID,self.formatCode,self.numWide,self.oCode,self.thermal)
        
        ## assuming tCode=1
        if self.approachCode==1:   # statics
            self.loadID = self.getValues(data,'i',5) ## load set ID number
        elif self.approachCode==2: # normal modes/buckling (real eigenvalues)
            self.mode      = self.getValues(data,'i',5) ## mode number
            self.eign      = self.getValues(data,'f',6) ## eigenvalue
        elif self.approachCode==3: # differential stiffness 0
            self.loadID = self.getValues(data,'i',5) ## load set ID number
        elif self.approachCode==4: # differential stiffness 1
            self.loadID = self.getValues(data,'i',5) ## load set ID number
        elif self.approachCode==5:   # frequency
            self.freq = self.getValues(data,'f',5) ## frequency

        elif self.approachCode==6: # transient
            self.time = self.getValues(data,'f',5) ## time step
            print "time(5)=%s" %(self.time)
        elif self.approachCode==7: # pre-buckling
            self.loadID = self.getValues(data,'i',5) ## load set ID number
            print "loadID(5)=%s" %(self.loadID)
        elif self.approachCode==8: # post-buckling
            self.loadID = self.getValues(data,'i',5) ## load set ID number
            self.eigr   = self.getValues(data,'f',6) ## real eigenvalue
            print "loadID(5)=%s  eigr(6)=%s" %(self.loadID,self.eigr)
        elif self.approachCode==9: # complex eigenvalues
            self.mode   = self.getValues(data,'i',5) ## mode
            self.eigr   = self.getValues(data,'f',6) ## real eigenvalue
            self.eigi   = self.getValues(data,'f',7) ## imaginary eigenvalue
            print "mode(5)=%s  eigr(6)=%s  eigi(7)=%s" %(self.mode,self.eigr,self.eigi)
        elif self.approachCode==10: # nonlinear statics
            self.loadStep = self.getValues(data,'f',5) ## load step
            print "loadStep(5) = %s" %(self.loadStep)
        elif self.approachCode==11: # geometric nonlinear statics
            self.loadID = self.getValues(data,'i',5) ## load set ID number
            print "loadID(5)=%s" %(self.loadID)
        else:
            raise RuntimeError('invalid approach code...approachCode=%s' %(self.approachCode))

        # tCode=2
        #if self.analysisCode==2: # sort2
        #    self.loadID = self.getValues(data,'i',5) ## load set ID number
        
        print "*iSubcase=%s"%(self.iSubcase)
        print "approachCode=%s tableCode=%s thermal=%s" %(self.approachCode,self.tableCode,self.thermal)
        print self.codeInformation()

        #self.printBlock(data)
        self.readTitle()

        #return (analysisCode,tableCode,thermal)

        #if self.j==3:
        #    #print str(self.obj)
        #    sys.exit('checkA...j=%s dt=6E-2 dx=%s dtActual=%f' %(self.j,'1.377e+01',self.dt))
        ###

    def readTable_OEF_4_Data(self,iTable): # iTable=-4
        isTable4Done = False
        isBlockDone  = False

        bufferWords = self.getMarker('OEF')
        #print len(bufferWords)
        self.data = self.readBlock()
        #self.printBlock(data)

        if bufferWords==146:  # table -4 is done, restarting table -3
            isTable4Done = True
            return isTable4Done,isBlockDone
        elif bufferWords==0:
            #print "bufferWords 0 - done with Table4"
            isTable4Done = True
            isBlockDone = True
            return isTable4Done,isBlockDone

        isBlockDone = not(bufferWords)
        self.readOEF1_Data()
        #print self.obj
        #del self.obj
        #print self.printSection(120)

        #print "-------finished OEF----------"
        return (isTable4Done,isBlockDone)


    def readOEF1_Data(self):
        tfsCode = [self.tableCode,self.formatCode,self.sortCode]
        if self.thermal==2:
            self.skipOES_Element(None)

        else: # thermal=0,1
            if   tfsCode==[4,1,0]:
                self.readOEF1_Data_format1_sort0()
            elif tfsCode==[4,1,1]:
                self.readOEF1_Data_format1_sort1()
            elif tfsCode==[4,2,1]:
                self.readOEF1_Data_format2_sort1()
            else:
                raise Exception('bad tableCode/formatCode/sortCode=%s on OEF table' %(tfsCode))
            ###
        ###
    def readOEF1_Data_format1_sort0(self):
        assert self.formatCode==1
        assert self.sortCode==0

        print "self.approachCode=%s tableCode(1)=%s sortCode=%s thermal(23)=%g" %(self.approachCode,self.tableCode,self.sortCode,self.thermal)
        if self.thermal==0:
            if self.approachCode==1: # displacement
                print "isForces"
                self.obj = displacementObject(self.iSubcase)
                self.displacementForces[self.iSubcase] = self.obj
                self.readForces(self.obj)

            elif self.approachCode==2 and self.sortCode==1: # buckling forces
                print "isBucklingForces"
                #self.obj = spcForcesObject(self.iSubcase)
                self.createTransientObject(self.bucklingForces,displacementObject,self.eigr)
                self.bucklingForces[self.iSubcase] = self.obj
            elif self.approachCode==5: # frequency forces
                print "isFrequencyForces"
                self.createTransientObject(self.freqForces,eigenVectorObject,self.freq)
                self.freqForces[self.iSubcase] = self.obj
            elif self.approachCode==6: # transient displacement
                print "isTransientForces"
                self.createTransientObject(self.displacementForces,displacementObject,self.time)
                self.displacementForces[self.iSubcase] = self.obj
                self.readForces(self.obj)

            elif self.approachCode==9: # complex eigenvalue forces
                print "isComplexEigenvalues"
                #self.obj = eigenVectorObject(self.iSubcase,self.eigr)
                self.createTransientObject(self.complexEigenvalueForces,eigenVectorObject,(self.mode,self.eigr,self.eigi))
                self.complexEigenvalueForces[self.iSubcase] = self.obj
                #print "****self", type(self.obj)
            elif self.approachCode==10: # nonlinear static displacement
                print "isNonlinearStaticForces"
                self.createTransientObject(self.nonlinearForces,displacementObject,self.loadStep)
                self.nonlinearForces[self.iSubcase] = self.obj
                ##self.readForcesNonlinear(self.obj)
                self.skipOES_Element(None)
            else:
                self.skipOES_Element(None)
                #raise Exception('not supported OEF static solution...')
            ###

        elif self.thermal==1:
            #if self.approachCode==1: # temperature
            #    print "isTemperature"
            #    raise Exception('verify...')
            #    self.temperatures[self.iSubcase] = temperatureObject(self.iSubcase)
            #elif self.approachCode==1 and self.sortCode==1: # heat fluxes
            #    print "isFluxes"
            #    raise Exception('verify...')
            #    self.obj = fluxObject(self.iSubcase,dt=self.dt)
            #    self.fluxes[self.iSubcase] = self.obj
            if self.approachCode==5: # frequency forces
                print "isFrequencyForces"
                self.createTransientObject(self.freqForces,eigenVectorObject,self.freq)
                self.freqForces[self.iSubcase] = self.obj
            elif self.approachCode==6: # transient temperature
                print "isTransientTemperature"
                #raise Exception('verify...')
                self.createTransientObject(self.temperatureForces,temperatureObject,self.time)
                self.temperatureForces[self.iSubcase] = self.obj  ## @todo modify the name of this...
                self.readForces(self.obj)
            elif self.approachCode==10: # nonlinear static displacement
                print "isNonlinearStaticTemperatures"
                self.createTransientObject(self.nonlinearFluxes,nonlinearFluxObject,self.loadStep)
                self.nonlinearFluxes[self.iSubcase] = self.obj
                self.readForcesNonlinear(self.obj)
            else:
                raise Exception('not supported OEF thermal solution...')
            ###
        else:
            msg = 'invalid thermal flag...not 0 or 1...flag=%s' %(self.thermal)
            self.obj = None
            sys.stderr.write(msg+'\n')
            #raise Exception(msg)
        ###
        #self.readForces(data,self.obj)
        #return

    def readOEF1_Data_format1_sort1(self):
        assert self.formatCode==1
        assert self.sortCode==1

        print "self.approachCode=%s tableCode(1)=%s thermal(23)=%g" %(self.approachCode,self.tableCode,self.thermal)
        if self.thermal==0:
            if self.approachCode==5: # frequency forces
                print "isFrequencyForces"
                self.createTransientObject(self.freqForces,eigenVectorObject,self.freq)
                self.freqForces[self.iSubcase] = self.obj
            elif self.approachCode==9: # complex eigenvalue forces
                print "isComplexEigenvalues"
                #self.obj = eigenVectorObject(self.iSubcase,self.eigr)
                self.createTransientObject(self.complexEigenvalueForces,eigenVectorObject,(self.mode,self.eigr,self.eigi))
                self.complexEigenvalueForces[self.iSubcase] = self.obj
                #print "****self", type(self.obj)
            else:
                raise Exception('not supported OEF static solution...')
            ###
        else:
            raise Exception('invalid thermal flag...not 0 or 1...flag=%s' %(self.thermal))
        ###
        #self.readForces(data,self.obj)
        #return
        

    def readOEF1_Data_format2_sort1(self):
        assert self.formatCode==2
        assert self.sortCode==1

        print "self.approachCode=%s tableCode(1)=%s thermal(23)=%g" %(self.approachCode,self.tableCode,self.thermal)
        if self.thermal==0:
            if self.approachCode==5: # frequency forces
                print "isFrequencyForces"
                self.createTransientObject(self.freqForces,eigenVectorObject,self.freq)
                self.freqForces[self.iSubcase] = self.obj
            elif self.approachCode==9: # complex eigenvalue forces
                print "isComplexEigenvalues"
                #self.obj = eigenVectorObject(self.iSubcase,self.eigr)
                self.createTransientObject(self.complexEigenvalueForces,eigenVectorObject,(self.mode,self.eigr,self.eigi))
                self.complexEigenvalueForces[self.iSubcase] = self.obj
                #print "****self", type(self.obj)
            else:
                raise Exception('not supported OEF static solution...')
            ###

        else:
            raise Exception('invalid thermal flag...not 0 or 1...flag=%s' %(self.thermal))
        ###
        #self.readForces(data,self.obj)
        #return

        
    def setupOEF(self,elementType,elementName,numWide,deviceCode):
        """length of result"""
        print "elementType=%s numWide=%s type=%s" %(elementType,numWide,elementName)
        
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
        return (recordLength,func,isSkipped) # 8*4

    def skipMe(self,data):
        pass

    def readOEF_2D_3D(self,data):
        print "read_2D_3D"
        gridDevice, = unpack('i',data[0:4])
        grid = (gridDevice-self.deviceCode)/10
        eType = ''.join(unpack('cccccccc',data[4:12]))
        (xGrad,yGrad,zGrad,xFlux,yFlux,zFlux) = unpack('ffffff',data[12:36])
        print "grid=%g dx=%i dy=%i dz=%i rx=%i ry=%i rz=%i" %(grid,xGrad,yGrad,zGrad,xFlux,yFlux,zFlux)
        
        return (grid,eType,xGrad,yGrad,zGrad,xFlux,yFlux,zFlux)

    def readOEF_CHBDY(self,data):
        print "read_CHBDYx"
        gridDevice, = unpack('i',data[0:4])
        grid = (gridDevice-self.deviceCode)/10
        eType = ''.join(unpack('cccccccc',data[4:12]))
        (fApplied,freeConv,forceConv,fRad,fTotal) = unpack('fffff',data[12:32])
        return (grid,eType,fApplied,freeConv,forceConv,fRad,fTotal)

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
    #    return (grid,eType,xGrad,yGrad,zGrad,xFlux,yFlux,zFlux)

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
        print "readForcesNonlinear..."
        data = self.data

        print 'thermal skipping elementType=%s=%s' %(self.elementType,self.ElementType(self.elementType))
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
        
