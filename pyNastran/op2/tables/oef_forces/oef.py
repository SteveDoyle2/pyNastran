from __future__ import division
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
from oef_Objects import (nonlinearFluxObject)

class Thermal_VU(object):
    def __init__(self,isSort1,dt):
        #self.eType = {}
        self.parent = {}
        self.coord = {}
        self.icord = {}
        self.theta = {}
        
        self.grad  = {}
        self.flux  = {}

        ## @todo if dt=None, handle SORT1 case
        if isSort1:
            self.add = self.addSort1
        else:
            self.add = self.addSort2
        ###

    def addNewTransient(self,dt):
        self.grad[dt] = {}
        self.flux[dt] = {}

    def addSort1(self,nNodes,dt,data):
        [eid,parent,coord,icord,theta,gradFluxes] = data
        if dt not in self.grad:
            self.addNewTransient(dt)
        self.parent[eid] = parent
        self.coord[eid] = coord
        self.icord[eid] = icord
        self.theta[eid] = theta
        #self.eType[eid]    = eType
        
        #self.grad[dt][eid] = [None,None,None]
        #self.flux[dt][eid] = [None,None,None]
        self.grad[dt][eid] = {}
        self.flux[dt][eid] = {}
        for gradFlux in gradFluxes:
            [nid,xGrad,yGrad,zGrad,xFlux,yFlux,zFlux] = gradFlux
            self.grad[dt][eid][nid] = [xGrad,yGrad,zGrad]
            self.flux[dt][eid][nid] = [xFlux,yFlux,zFlux]

    def addSort2(self,nNodes,eid,data):
        [dt,parent,coord,icord,theta,gradFluxes] = data
        if dt not in self.fApplied:
            self.addNewTransient(dt)
        self.parent[eid] = parent
        self.coord[eid] = coord
        self.icord[eid] = icord
        self.theta[eid] = theta
        #self.eType[eid]    = eType

        #self.grad[dt][eid] = [None,None,None]
        #self.flux[dt][eid] = [None,None,None]
        self.grad[dt][eid] = {}
        self.flux[dt][eid] = {}
        for gradFlux in gradFluxes:
            [nid,xGrad,yGrad,zGrad,xFlux,yFlux,zFlux] = gradFlux
            self.grad[dt][eid][nid] = [xGrad,yGrad,zGrad]
            self.flux[dt][eid][nid] = [xFlux,yFlux,zFlux]

    def __repr__(self):
        return str(self.grad)
        
class Thermal_1D(object):
    def __init__(self,isSort1,dt):
        self.eType = {}
        self.grad  = {}
        self.flux  = {}

        ## @todo if dt=None, handle SORT1 case
        if isSort1:
            self.add = self.addSort1
        else:
            self.add = self.addSort2
        ###

    def addNewTransient(self,dt):
        self.grad[dt] = {}
        self.flux[dt] = {}

    def addSort1(self,dt,data):
        [eid,eType,xGrad,yGrad,zGrad,xFlux,yFlux,zFlux] = data
        if dt not in self.grad:
            self.addNewTransient(dt)
        self.eType[eid]    = eType
        self.grad[dt][eid] = [xGrad,yGrad,zGrad]
        self.flux[dt][eid] = [xFlux,yFlux,zFlux]

    def addSort2(self,eid,data):
        [dt,eType,xGrad,yGrad,zGrad,xFlux,yFlux,zFlux] = data
        if dt not in self.fApplied:
            self.addNewTransient(dt)
        self.eType[eid]    = eType
        self.grad[dt][eid] = [xGrad,yGrad,zGrad]
        self.flux[dt][eid] = [xFlux,yFlux,zFlux]

    def __repr__(self):
        return str(self.grad)
        
class Thermal_2D_3D(object):
    def __init__(self,isSort1,dt):
        self.eType = {}
        self.grad  = {}
        self.flux  = {}

        ## @todo if dt=None, handle SORT1 case
        if isSort1:
            self.add = self.addSort1
        else:
            self.add = self.addSort2
        ###

    def addNewTransient(self,dt):
        self.grad[dt] = {}
        self.flux[dt] = {}

    def addSort1(self,dt,data):
        [eid,eType,xGrad,yGrad,zGrad,xFlux,yFlux,zFlux] = data
        if dt not in self.grad:
            self.addNewTransient(dt)
        self.eType[eid]    = eType
        self.grad[dt][eid] = [xGrad,yGrad,zGrad]
        self.flux[dt][eid] = [xFlux,yFlux,zFlux]

    def addSort2(self,eid,data):
        [dt,eType,xGrad,yGrad,zGrad,xFlux,yFlux,zFlux] = data
        if dt not in self.fApplied:
            self.addNewTransient(dt)
        self.eType[eid]    = eType
        self.grad[dt][eid] = [xGrad,yGrad,zGrad]
        self.flux[dt][eid] = [xFlux,yFlux,zFlux]

    def __repr__(self):
        return str(self.grad)
        
class CHBDYx(object): # [107,108,109] # CHBDYE, CHBDYG, CHBDYP
    def __init__(self,isSort1,dt):
        self.eType     = {}
        self.fApplied  = {}
        self.freeConv  = {}
        self.forceConv = {}
        self.fRad      = {}
        self.fTotal    = {}
        
        ## @todo if dt=None, handle SORT1 case
        if isSort1:
            self.add = self.addSort1
        else:
            self.add = self.addSort2
        ###

    def addNewTransient(self,dt):
        self.fApplied[dt]  = {}
        self.freeConv[dt]  = {}
        self.forceConv[dt] = {}
        self.fRad[dt]      = {}
        self.fTotal[dt]    = {}
        
    def addSort1(self,dt,data):
        [eid,eType,fApplied,freeConv,forceConv,fRad,fTotal] = data
        if dt not in self.fApplied:
            self.addNewTransient(dt)
        self.eType[eid]     = eType
        self.fApplied[dt][eid]  = fApplied
        self.freeConv[dt][eid]  = freeConv
        self.forceConv[dt][eid] = forceConv
        self.fRad[dt][eid]      = fRad
        self.fTotal[dt][eid]    = fTotal

    def addSort2(self,eid,data):
        [dt,eType,fApplied,freeConv,forceConv,fRad,fTotal] = data
        if dt not in self.fApplied:
            self.addNewTransient(dt)
        self.eType[eid]     = eType
        self.fApplied[dt][eid]  = fApplied
        self.freeConv[dt][eid]  = freeConv
        self.forceConv[dt][eid] = forceConv
        self.fRad[dt][eid]      = fRad
        self.fTotal[dt][eid]    = fTotal

    def __repr__(self):
        return str(self.fApplied)

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
        #print '-'*80
        self.readTitle()

    def updateSort1(self):
        extract = self.extractSort1
        return 'i',extract

    def updateSort2(self):
        extract = self.extractSort2
        return 'f'

    def extractSort1(self,eidDevice,dt):
        #eidDevice, = unpack('i',data)
        #print "eidDevice=%s dt=%s" %(eidDevice,dt)
        return (eidDevice-self.deviceCode)//10

    def extractSort2(self,timeFreq,eid):
        #print "timeFreq=%s eid=%s" %(timeFreq,eid)
        #gridDevice, = unpack('i',data)
        return timeFreq
    
    def OEF_CHBDYx(self): # [107,108,109]  CHBDYE, CHBDYG, CHBDYP
        if self.makeOp2Debug:
            self.op2Debug.write('---OEF_CHBDYx---\n')
        deviceCode = self.deviceCode
        
        dt = self.nonlinearFactor
        isSort1 = self.isSort1()
        if isSort1:
            #print "SORT1 - %s" %(self.ElementType(self.elementType))
            format1 = 'iccccccccfffff' # SORT1
            extract = self.extractSort1
            #dt = self.nonlinearFactor
        else:
            #print "SORT2 - %s" %(self.ElementType(self.elementType))
            format1 = 'fccccccccfffff' # SORT2
            extract = self.extractSort2
            #eid = self.nonlinearFactor

        #self.obj = CHBDYx(isSort1,dt)
        self.createThermalTransientObject(self.thermalLoad_CHBDY,CHBDYx,isSort1)

        while len(self.data)>=32: # 8*4
            eData     = self.data[0:32]
            self.data = self.data[32: ]
            #print "len(data) = ",len(eData)

            out = unpack(format1, eData)
            (eid,a,b,c,d,e,f,g,h,fApplied,freeConv,forceConv,fRad,fTotal) = out
            eid2  = extract(eid,dt)
            eType = a+b+c+d+e+f+g+h
            #print "eType=%s" %(eType)
            
            dataIn = [eid2,eType,fApplied,freeConv,forceConv,fRad,fTotal]
            #print dataIn
            #eid = self.obj.addNewEid(out)
            self.obj.add(dt,dataIn)
            #print "len(data) = ",len(self.data)
        ###
        self.handleResultsBuffer(self.OEF_CHBDYx)
        if self.makeOp2Debug:
            print "done with OEF_CHBDYx"
        
        #print self.thermalLoad_CHBDY
        #sys.exit('end of CHBDYx')
        #raise Exception('add OEF_CHBDYx...')

    #def makeResult
    def createThermalTransientObject(self,resultName,objClass,isSort1):
        #print resultName
        if self.iSubcase in resultName:
            self.obj = resultName[self.iSubcase]
            print "returning iSubcase result=%s" %(self.iSubcase)
        else:
            self.obj = objClass(isSort1,self.nonlinearFactor)
            resultName[self.iSubcase] = self.obj
            print "creating iSubcase result=%s" %(self.iSubcase)
        ###
        #return self.obj
    
    def OEF_VU_Element(self): # 189-VUQUAD 190-VUTRIA,191-VUBEAM
        dt = self.nonlinearFactor
        isSort1 = self.isSort1()
        print "numWide = ",self.numWide
        n = 36
        if self.elementType in [189]:
            nNodes = 4
        elif self.elementType in [190]:
            nNodes = 3
        elif self.elementType in [191]:
            nNodes = 2
        else:
            raise NotImplementedError(self.codeInformation())

        if isSort1:
            print "SORT1 - %s" %(self.ElementType(self.elementType))
            format1 = 'iiiccccii' # SORT1
            extract = self.extractSort1
            #dt = self.nonlinearFactor
        else:
            print "SORT2 - %s" %(self.ElementType(self.elementType))
            format1 = 'fiiccccii' # SORT2
            extract = self.extractSort2
            #eid = self.nonlinearFactor
        ###
        formatAll = 'iffffff'#*nNodes
        self.createThermalTransientObject(self.thermalLoad_2D_3D,Thermal_VU,isSort1)

        n = 24+28*nNodes
        while len(self.data)>=n: # 10*4
            eData     = self.data[0:6*4]
            self.data = self.data[6*4: ]
            print "len(data) = ",len(eData)

            out = unpack(format1, eData)
            (eid,parent,coord,icordA,icordB,icordC,icordD,theta,null) = out
            icord = icordA+icordB+icordC+icordD

            eid2  = extract(eid,dt)
            dataIn = [eid2,parent,coord,icord,theta]

            gradFluxes = []
            for i in range(nNodes):
                eData     = self.data[0:7*4]
                self.data = self.data[7*4: ]
                print "i=%s len(data)=%s" %(i,len(eData))
                out = unpack(formatAll, eData)
                gradFluxes.append(out)
            dataIn.append(gradFluxes)
            #eType = a+b+c+d+e+f+g+h
            #print "eType=%s" %(eType)
            
            #dataIn = [eid2,eType,xGrad,yGrad,zGrad,xFlux,yFlux,zFlux]
            print dataIn
            #eid = self.obj.addNewEid(out)            
            self.obj.add(nNodes,dt,dataIn)
            #print "len(data) = ",len(self.data)
        ###
        self.handleResultsBuffer(self.OEF_1D)
        if self.makeOp2Debug:
            print "done with OEF_1D"
        
        print self.thermalLoad_2D_3D


    def OEF_1D(self): # 1-ROD, 2-BEAM, 3-TUBE, 10-CONROD, 34-BAR, 69-BEND
        dt = self.nonlinearFactor
        isSort1 = self.isSort1()
        print "numWide = ",self.numWide
        n = 36
        if isSort1:
            print "SORT1 - %s" %(self.ElementType(self.elementType))
            format1 = 'iccccccccffffff' # SORT1
            extract = self.extractSort1
            #dt = self.nonlinearFactor
        else:
            print "SORT2 - %s" %(self.ElementType(self.elementType))
            format1 = 'fccccccccffffff' # SORT2
            extract = self.extractSort2
            #eid = self.nonlinearFactor
        ###
        self.createThermalTransientObject(self.thermalLoad_2D_3D,Thermal_1D,isSort1)
 
        n = 36
        while len(self.data)>=n: # 10*4
            eData     = self.data[0:n]
            self.data = self.data[n: ]
            #print "len(data) = ",len(eData)

            out = unpack(format1, eData)
            (eid,a,b,c,d,e,f,g,h,xGrad,yGrad,zGrad,xFlux,yFlux,zFlux) = out
            eid2  = extract(eid,dt)
            eType = a+b+c+d+e+f+g+h
            #print "eType=%s" %(eType)
            
            dataIn = [eid2,eType,xGrad,yGrad,zGrad,xFlux,yFlux,zFlux]
            print dataIn
            #eid = self.obj.addNewEid(out)            
            self.obj.add(dt,dataIn)
            #print "len(data) = ",len(self.data)
        ###
        self.handleResultsBuffer(self.OEF_1D)
        if self.makeOp2Debug:
            print "done with OEF_1D"
        
        print self.thermalLoad_2D_3D

    def OEF_2D_3D(self): # 33-QUAD4, 39-TETRA, 53-TRIAX6,64-QUAD8, 67-HEXA, 68-PENTA, 74-TRIA3, 75-TRIA6
        """numWide==10"""
        dt = self.nonlinearFactor
        isSort1 = self.isSort1()
        print "numWide = ",self.numWide
        if self.elementType in [39,67,68]: # HEXA,PENTA
            n = 40
            if isSort1:
                print "SORT1 - %s" %(self.ElementType(self.elementType))
                format1 = 'iccccccccffffffi' # SORT1
                extract = self.extractSort1
                #dt = self.nonlinearFactor
            else:
                print "SORT2 - %s" %(self.ElementType(self.elementType))
                format1 = 'fccccccccffffffi' # SORT2
                extract = self.extractSort2
                #eid = self.nonlinearFactor
        elif self.elementType in [33,53,64,74,75]: # no zed on this element for some reason...
            n = 36
            if isSort1:
                print "SORT1 - %s" %(self.ElementType(self.elementType))
                format1 = 'iccccccccffffff' # SORT1
                extract = self.extractSort1
                #dt = self.nonlinearFactor
            else:
                print "SORT2 - %s" %(self.ElementType(self.elementType))
                format1 = 'fccccccccffffff' # SORT2
                extract = self.extractSort2
                #eid = self.nonlinearFactor
            ###
        else:
            raise NotImplementedError(self.codeInformation())
        ###
        self.createThermalTransientObject(self.thermalLoad_2D_3D,Thermal_2D_3D,isSort1)

        #n = 36
        while len(self.data)>=n: # 10*4
            eData     = self.data[0:n]
            self.data = self.data[n: ]

            out = unpack(format1, eData)
            (eid,a,b,c,d,e,f,g,h,xGrad,yGrad,zGrad,xFlux,yFlux,zFlux) = out[0:15]
            eid2  = extract(eid,dt)
            eType = a+b+c+d+e+f+g+h
            #print "eType=%s" %(eType)
            
            dataIn = [eid2,eType,xGrad,yGrad,zGrad,xFlux,yFlux,zFlux]
            print dataIn
            #eid = self.obj.addNewEid(out)            
            self.obj.add(dt,dataIn)
            #print "len(data) = ",len(self.data)
        ###
        self.handleResultsBuffer(self.OEF_2D_3D)
        if self.makeOp2Debug:
            print "done with OEF_2D_3D"
        
        print self.thermalLoad_2D_3D

    def readOEF_Thermal(self):
        #print "self.elementType = ",self.elementType
        if self.elementType in [107,108,109]: # CHBDYE, CHBDYG, CHBDYP
            assert self.numWide==8,self.codeInformation()
            self.OEF_CHBDYx()
        elif self.elementType in [33,39,67,68]: # QUAD4,TETRA,HEXA,PENTA
            assert self.numWide in [9,10],self.codeInformation()
            self.OEF_2D_3D()
        elif self.elementType in [53,64,74,75]: # TRIAX6,QUAD8,TRIA3,TRIA6
            assert self.numWide==9,self.codeInformation()
            self.OEF_2D_3D()
        elif self.elementType in [1,2,3,10,34,69]: # ROD,BEAM,TUBE,CONROD,BAR,BEND
            assert self.numWide==9,self.codeInformation()
            self.OEF_1D()
        elif self.elementType in [189,190,191]: # VUQUAD,VUTRIA,VUBEAM
            #assert self.numWide==27,self.codeInformation()
            self.OEF_VU_Element()
        else:
            print self.codeInformation()
            raise NotImplementedError('stopping in bad OEF thermal element')
        ###

    def readOEF1_Data(self):
        tfsCode = [self.tableCode,self.formatCode,self.sortCode]
        #self.skipOES_Element() # skipping entire table
        #self.skippedCardsFile.write('skipping atfsCode=%s in %s-OEF\n' %(self.atfsCode,self.tableName))
        
        if self.thermal==1:
            self.readOEF_Thermal()
        else:
            #sys.exit('non-thermal results in OEF not supported')
            self.skipOES_Element() # skipping entire table
        #sys.exit('stopping in OEF')
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

    #def readOEF_2D_3D(self,data):
        #print "read_2D_3D"
        #gridDevice, = unpack('i',data[0:4])
        #grid = (gridDevice-self.deviceCode)/10
        #eType = ''.join(unpack('cccccccc',data[4:12]))
        #(xGrad,yGrad,zGrad,xFlux,yFlux,zFlux) = unpack('ffffff',data[12:36])
        #print "grid=%g dx=%i dy=%i dz=%i rx=%i ry=%i rz=%i" %(grid,xGrad,yGrad,zGrad,xFlux,yFlux,zFlux)
        #return(grid,eType,xGrad,yGrad,zGrad,xFlux,yFlux,zFlux)

    #def readOEF_CHBDY(self,data):
        #print "read_CHBDYx"
        #gridDevice, = unpack('i',data[0:4])
        #grid = (gridDevice-self.deviceCode)/10
        #eType = ''.join(unpack('cccccccc',data[4:12]))
        #(fApplied,freeConv,forceConv,fRad,fTotal) = unpack('fffff',data[12:32])
        #return(grid,eType,fApplied,freeConv,forceConv,fRad,fTotal)

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
        
