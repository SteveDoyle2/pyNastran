import sys
import copy
from struct import unpack

# pyNastran
from pyNastran.op2.resultObjects.op2_Objects import spcForcesObject
from pyNastran.op2.tables.oug.ougv1_Objects import (
     temperatureObject,displacementObject,  # analysisCode=1, sortCode=0
     fluxObject,                            # analysisCode=1, sortCode=3
     nonlinearTemperatureObject,            # analysisCode=10,sortCode=0
     )

from pyNastran.op2.tables.oug.oug_eigenvectors import (
     eigenVectorObject,                     # analysisCode=2, sortCode=0 formatCode   tableCode=7
     complexEigenVectorObject,              # analysis    =5, sortCode=1 formatCode=1 tableCode=7
     realEigenVectorObject,                 # analysisCode=9, sortCode=1 formatCode=1 tableCode=7
     )

class OQG1(object):
    """Table of spc/mpc forces/momenets"""

    def readTable_OQG1(self):
        self.tableName = 'OQG'
        table3 = self.readTable_OQG1_3
        table4Data = self.readOQG1_Data
        self.readResultsTable(table3,table4Data)
        self.deleteAttributes_OQG()

    def deleteAttributes_OQG(self):
        params = ['lsdvm','mode','eigr','modeCycle','freq','dt','lftsfq','thermal','rCode','fCode','numWide','acousticFlag','thermal']
        self.deleteAttributes(params)
    
    def readTable_OQG1_3(self,iTable): # iTable=-3
        bufferWords = self.getMarker()
        if self.makeOp2Debug:
            self.op2Debug.write('bufferWords=%s\n' %(str(bufferWords)))
        #print "2-bufferWords = ",bufferWords,bufferWords*4,'\n'

        data = self.getData(4)
        bufferSize, = unpack('i',data)
        data = self.getData(4*50)
        #print self.printBlock(data)
        
        (three) = self.parseApproachCode(data)

        self.rCode        = self.getValues(data,'i',8)  ## random code
        self.formatCode   = self.getValues(data,'i',9)  ## format code
        self.numWide      = self.getValues(data,'i',10) ## number of words per entry in record; @note is this needed for this table ???
        self.acousticFlag = self.getValues(data,'f',13) ## acoustic pressure flag
        self.thermal      = self.getValues(data,'i',23) ## thermal flag; 1 for heat ransfer, 0 otherwise
        
        self.dataCode = {'analysisCode': self.analysisCode,'deviceCode':self.deviceCode,
                         'randomCode':self.rCode,'formatCode':self.formatCode,
                         'numWide': self.numWide,'acousticFlag':self.acousticFlag,
                         'thermal': self.thermal}

        #self.printBlock(data) # on
        ## assuming tCode=1
        if self.analysisCode==1:   # statics / displacement / heat flux
            self.lsdvmn = self.getValues(data,'i',5) ## load set number
            self.nonlinearFactor = self.lsdvmn
            self.dataCode['lsdvmn'] = self.lsdvmn
        elif self.analysisCode==2: # real eigenvalues
            self.mode      = self.getValues(data,'i',5) ## mode number
            self.eigr      = self.getValues(data,'f',6) ## real eigenvalue
            self.modeCycle = self.getValues(data,'f',7) ## mode or cycle @todo confused on the type - F1???
            self.dataCode['mode'] = self.mode
            self.dataCode['eigr'] = self.eigr
            self.dataCode['modeCycle'] = self.modeCycle
            print "mode(5)=%s eigr(6)=%s modeCycle(7)=%s" %(self.mode,self.eigr,self.modeCycle)
            self.nonlinearFactor = self.mode
        #elif self.analysisCode==3: # differential stiffness
        #    self.lsdvmn = self.getValues(data,'i',5) ## load set number
        #    self.nonlinearFactor = self.lsdvmn
        #    self.dataCode['lsdvmn'] = self.lsdvmn
        #elif self.analysisCode==4: # differential stiffness
        #    self.lsdvmn = self.getValues(data,'i',5) ## load set number
        #    self.nonlinearFactor = self.lsdvmn
        #    self.dataCode['lsdvmn'] = self.lsdvmn
        elif self.analysisCode==5:   # frequency
            self.freq = self.getValues(data,'f',5) ## frequency
            self.nonlinearFactor = self.freq
            self.dataCode['freq'] = self.freq
        elif self.analysisCode==6: # transient
            self.dt = self.getValues(data,'f',5) ## time step
            self.nonlinearFactor = self.dt
            self.dataCode['dt'] = self.dt
            print "DT(5)=%s" %(self.dt)
        elif self.analysisCode==7: # pre-buckling
            self.lsdvmn = self.getValues(data,'i',5) ## load set
            self.nonlinearFactor = self.lsdvmn
            self.dataCode['lsdvmn'] = self.lsdvmn
            print "LSDVMN(5)=%s" %(self.lsdvmn)
        elif self.analysisCode==8: # post-buckling
            self.lsdvmn = self.getValues(data,'i',5) ## mode number
            self.eigr   = self.getValues(data,'f',6) ## real eigenvalue
            self.nonlinearFactor = self.lsdvmn
            self.dataCode['lsdvmn'] = self.lsdvmn
            self.dataCode['eigr'] = self.eigr
            print "LSDVMN(5)=%s  EIGR(6)=%s" %(self.lsdvmn,self.eigr)
        elif self.analysisCode==9: # complex eigenvalues
            self.mode   = self.getValues(data,'i',5) ## mode
            self.eigr   = self.getValues(data,'f',6) ## real eigenvalue
            self.eigi   = self.getValues(data,'f',7) ## imaginary eigenvalue
            self.nonlinearFactor  = self.mode
            self.dataCode['mode'] = self.mode
            self.dataCode['eigr'] = self.eigr
            self.dataCode['eigi'] = self.eigi
            print "mode(5)=%s  eigr(6)=%s  eigi(7)=%s" %(self.mode,self.eigr,self.eigi)
        elif self.analysisCode==10: # nonlinear statics
            self.lftsfq = self.getValues(data,'f',5) ## load step
            self.dataCode['lftsfq'] = self.lftsfq
            self.nonlinearFactor = self.lftsfq
            print "LFTSFQ(5) = %s" %(self.lftsfq)
        elif self.analysisCode==11: # old geometric nonlinear statics
            self.lsdvmn = self.getValues(data,'i',5)
            self.nonlinearFactor = self.lsdvmn
            self.dataCode['lsdvmn'] = self.lsdvmn
            print "LSDVMN(5)=%s" %(self.lsdvmn)
        elif self.analysisCode==12: # contran ? (may appear as aCode=6)  --> straight from DMAP...grrr...
            self.lsdvmn = self.getValues(data,'i',5)
            self.nonlinearFactor = self.lsdvmn
            self.dataCode['lsdvmn'] = self.lsdvmn
            print "LSDVMN(5)=%s" %(self.lsdvmn)
        else:
            raise RuntimeError('invalid analysis code...analysisCode=%s' %(self.analysisCode))
        # tCode=2
        #if self.analysisCode==2: # sort2
        #    self.lsdvmn = self.getValues(data,'i',5)
        
        #print "*iSubcase=%s"%(self.iSubcase)
        #print "analysisCode=%s tableCode=%s thermal=%s" %(self.analysisCode,self.tableCode,self.thermal)
        #print self.codeInformation()

        #self.printBlock(data)
        self.readTitle()

    def readOQG1_Data(self):
        tfsCode = [self.tableCode,self.formatCode,self.sortCode]
        print "self.analysisCode=%s tableCode(1)=%s thermal(23)=%g" %(self.analysisCode,self.tableCode,self.thermal)
        assert self.thermal in [0,1]

        if tfsCode==[3,1,0]:  # SPC Force vector
            self.readOQG1_Data_table3_format1_sort0()
        #elif tfsCode==[3,1,1]:
        #    self.readOQG1_Data_format1_sort1()
        #elif tfsCode==[3,2,0]:
        #    self.readOQG1_Data_format2_sort0()
        #elif tfsCode==[3,2,1]:
        #    self.readOQG1_Data_format2_sort1()
        #elif tfsCode==[3,3,0]:
        #    self.readOQG1_Data_format3_sort0()
        #elif tfsCode==[3,3,1]:
        #    self.readOQG1_Data_format3_sort1()
        
        #elif tfsCode==[39,1,0]: # MPC Forces
        #    self.readOQG1_Data_table39_format1_sort0()
        #elif tfsCode==[39,1,1]:
        #    self.readOQG1_Data_format1_sort1()
        else:
            self.skipOES_Element(None)
            #raise Exception('bad analysis/table/format/sortCode=%s' %(self.atfsCode))
        ###
        #print self.obj

    def readOQG1_Data_format2_sort0(self):
        print 'not supported OQG solution...atfsCode=%s' %(self.atfsCode)
        self.skipOES_Element(None)

    def readOQG1_Data_format3_sort0(self):
        print 'not supported OQG solution...atfsCode=%s' %(self.atfsCode)
        self.skipOES_Element(None)

    def readOQG1_Data_table3_format1_sort0(self):
        if self.thermal==0:
            if self.analysisCode==1: # displacement
                print "isSPCForces"
                self.createTransientObject(self.spcForces,spcForcesObject)
            #elif self.analysisCode==2: # nonlinear static eigenvector
            #    print "isEigenvector"
            #    self.createTransientObject(self.modalSPCForces,eigenVectorObject)
            #elif self.analysisCode==5: # frequency
            #    print "isFrequencyForces"
            #    self.createTransientObject(self.modalSPCForces,eigenVectorObject)
            elif self.analysisCode==6: # transient forces
                print "isTransientSPCForces"
                self.createTransientObject(self.spcForces,spcForcesObject)
            #elif self.analysisCode==10: # nonlinear static displacement
            #    print "isNonlinearStaticDisplacement"
            #    self.createTransientObject(self.spcForces,displacementObject)
            #elif self.analysisCode==11: # Geometric nonlinear statics
            #    print "isNonlinearStaticDisplacement"
            #    self.createTransientObject(self.spcForces,displacementObject)
            else:
                self.skipOES_Element(None)
                #raise Exception('unsupported OQG1 static solution...atfsCode=%s' %(self.atfsCode))
            ###
        elif self.thermal==1:
            #if self.analysisCode==1: # temperature
            #    print "isTemperature"
            #    self.obj = temperatureObject(self.dataCode,self.iSubcase)
            #    self.temperatures[self.iSubcase] = self.obj
            #elif self.analysisCode==6: # transient forces
            #    print "isTransientForces"
            #    self.createTransientObject(self.nonlinearTemperatures,nonlinearTemperatureObject)
            #elif self.analysisCode==10: # nonlinear static displacement
            #    print "isNonlinearStaticTemperatures"
            #    self.createTransientObject(self.nonlinearTemperatures,nonlinearTemperatureObject)
            if 0:
                pass
            else:
                self.skipOES_Element(None)
                print 'unsupported OQG1 thermal solution...atfsCode=%s' %(self.atfsCode)
                #raise Exception('unsupported OQG1 thermal solution...atfsCode=%s' %(self.atfsCode))
            ###
        ###
        else:
            raise Exception('invalid OQG1 thermal flag...not 0 or 1...flag=%s' %(self.thermal))
        ###
        #print "objName = ",self.obj.name()
        self.readScalars8(self.obj)

    def readOQG1_Data_format1_sort1(self):
        if self.thermal==0:
            if self.analysisCode==5: # frequency
                print "isFrequencyForces"
                self.createTransientObject(self.modalSPCForces,complexEigenVectorObject)
            elif self.analysisCode==9: # frequency
                print "isComplexEigenvalueForces"
                self.createTransientObject(self.modalSPCForces,complexEigenVectorObject)
            ###
            else:
                raise Exception('unsupported OQG1 static solution...atfsCode=%s' %(self.atfsCode))
            ###
        elif thermal==1:
            raise Exception('unsupported OQG1 thermal solution...atfsCode=%s' %(self.atfsCode))
        else:
            raise Exception('invalid OQG1 thermal flag...not 0 or 1...flag=%s' %(self.thermal))
        ###
        #print "objName = ",self.obj.name()
        self.readScalars14(self.obj)

    def readOQG1_Data_format2_sort1(self):
        if self.thermal==0:
            if self.analysisCode==5: # frequency
                print "isFrequencyForces"
                self.createTransientObject(self.modalSPCForces,complexEigenVectorObject)
            else:
                raise Exception('unsupported OQG1 static solution...atfsCode=%s' %(self.atfsCode))
            ###
        elif thermal==1:
            raise Exception('unsupported OQG1 thermal solution...atfsCode=%s' %(self.atfsCode))
        else:
            raise Exception('invalid OQG1 thermal flag...not 0 or 1...flag=%s' %(self.thermal))
        ###
        #print "objName = ",self.obj.name()
        self.readScalars14(self.obj) # readImaginary

    def readOQG1_Data_format3_sort1(self):
        if self.thermal==0:
            if self.analysisCode==5: # frequency
                print "isFrequencyForces"
                self.createTransientObject(self.modalSPCForces,complexEigenVectorObject)
            elif self.analysisCode==6: # transient forces
                print "isTransientForces"
                self.createTransientObject(self.spcForces,spcForcesObject)
            elif self.analysisCode==8: # post-buckling forces
                print "isPostBucklingForces"
                self.createTransientObject(self.spcForces,spcForcesObject)
            elif self.analysisCode==11: # Geometric nonlinear statics
                print "isFrequencyForces"
                self.createTransientObject(self.nonlinearForces,eigenVectorObject)
            else:
                raise Exception('unsupported OQG1 static solution...atfsCode=%s' %(self.atfsCode))
            ###
        elif thermal==1:
            raise Exception('unsupported OQG1 thermal solution...atfsCode=%s' %(self.atfsCode))
        ###
        else:
            raise Exception('invalid OQG1 thermal flag...not 0 or 1...flag=%s' %(self.thermal))
        ###
        #print "objName = ",self.obj.name()
        self.readScalars14(self.obj) # readImaginary

