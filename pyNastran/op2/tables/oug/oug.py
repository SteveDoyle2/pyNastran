import sys
import copy
from numpy import array
from struct import unpack

# pyNastran
#from oug_Objects import (
    #fluxObject,                            # analysisCode=1, formatCode=1 sortCode=3
    #nonlinearTemperatureObject,            # analysisCode=10,formatCode=1 sortCode=0 ???
#     )

from oug_displacements import (
     displacementObject,                    # tableCode=1, formatCode=1 sortCode=0
     complexDisplacementObject,             # analysisCode=5  formatCode=3 sortCode=1
     )

from oug_velocities import (                # tableCode=10,formatCode=1 sortCode=0
     velocityObject,
     complexVelocityObject
     )

from oug_accelerations import (             # tableCode=11,formatCode=1 sortCode=0
     accelerationObject,
     complexAccelerationObject
     )

from oug_temperatures import (              # tableCode=1, formatCode=1 sortCode=0
     temperatureObject,
     )

from oug_eigenvectors import (
     eigenVectorObject,                     # analysisCode=2, sortCode=0 formatCode   tableCode=7
     complexEigenVectorObject,              # analysisCode=5, sortCode=1 formatCode=1 tableCode=7
     realEigenVectorObject,                 # analysisCode=9, sortCode=1 formatCode=1 tableCode=7
     )
from pyNastran.op2.op2_helper import polarToRealImag

class OUG(object):
    """Table of displacements/velocities/acceleration/heat flux/temperature"""

    def readTable_OUG(self):
        #self.tableName = 'OUG'
        table3 = self.readTable_OUG_3
        table4Data = self.readOUG_Data
        self.readResultsTable(table3,table4Data)
        self.deleteAttributes_OUG()

    def deleteAttributes_OUG(self):
        params = ['lsdvm','mode','eigr','modeCycle','freq','dt','lftsfq','thermal','randomCode','fCode','numWide','acousticFlag']
        self.deleteAttributes(params)
    
    def addDataParameter(self,data,Name,Type,FieldNum,applyNonlinearFactor=True):
        #self.mode      = self.getValues(data,'i',5) ## mode number
        value = self.getValues(data,Type,FieldNum)
        setattr(self,Name,value)
        self.dataCode[Name] = value
        
        if applyNonlinearFactor:
            self.nonlinearFactor = value
            self.dataCode['nonlinearFactor'] = value
            self.dataCode['name'] = Name
    
    def applyDataCodeValue(self,Name,value):
        self.dataCode[Name] = value
        
    def readTable_OUG_3(self,iTable): # iTable=-3
        bufferWords = self.getMarker()
        if self.makeOp2Debug:
            self.op2Debug.write('bufferWords=%s\n' %(str(bufferWords)))
        #print "2-bufferWords = ",bufferWords,bufferWords*4,'\n'

        data = self.getData(4)
        bufferSize, = unpack('i',data)
        data = self.getData(4*50)
        #print self.printBlock(data)
        
        (three) = self.parseApproachCode(data)

        self.addDataParameter(data,'randomCode',  'i',8,False)   ## random code
        self.addDataParameter(data,'formatCode',  'i',9,False)   ## format code
        self.addDataParameter(data,'numWide',     'i',10,False)  ## number of words per entry in record; @note is this needed for this table ???
        self.addDataParameter(data,'acousticFlag','f',13,False)  ## acoustic pressure flag
        self.addDataParameter(data,'thermal',     'i',23,False)  ## thermal flag; 1 for heat transfer, 0 otherwise
        
        if not self.isSort1():
            raise NotImplementedError('sort2...')
        assert self.isThermal()==False,self.isThermal

        ## assuming tCode=1
        if self.analysisCode==1:   # statics / displacement / heat flux
            self.addDataParameter(data,'lsdvmn',  'i',5,False)   ## load set number
            self.applyDataCodeValue('dataNames',['lsdvmn'])
        elif self.analysisCode==2: # real eigenvalues
            self.addDataParameter(data,'mode',     'i',5)         ## mode number
            self.addDataParameter(data,'eigr',     'f',6,False)   ## real eigenvalue
            self.addDataParameter(data,'modeCycle','i',7,False)   ## mode or cycle @todo confused on the type - F1???
            self.applyDataCodeValue('dataNames',['mode','eigr','modeCycle'])
        #elif self.analysisCode==3: # differential stiffness
            #self.lsdvmn = self.getValues(data,'i',5) ## load set number
            #self.dataCode['lsdvmn'] = self.lsdvmn
        #elif self.analysisCode==4: # differential stiffness
            #self.lsdvmn = self.getValues(data,'i',5) ## load set number
        elif self.analysisCode==5:   # frequency
            self.addDataParameter(data,'freq','f',5)   ## frequency
            self.applyDataCodeValue('dataNames',['freq'])
        elif self.analysisCode==6: # transient
            self.addDataParameter(data,'dt','f',5)   ## time step
        elif self.analysisCode==7: # pre-buckling
            self.addDataParameter(data,'lsdvmn',  'i',5)   ## load set number
            self.applyDataCodeValue('dataNames',['lsdvmn'])
        elif self.analysisCode==8: # post-buckling
            self.addDataParameter(data,'lsdvmn',  'i',5)   ## load set number
            self.addDataParameter(data,'eigr',    'f',6,False)   ## real eigenvalue
            self.applyDataCodeValue('dataNames',['lsdvmn','eigr'])
        elif self.analysisCode==9: # complex eigenvalues
            self.addDataParameter(data,'mode','i',5)   ## mode number
            self.addDataParameter(data,'eigr','f',6,False)   ## real eigenvalue
            self.addDataParameter(data,'eigi','f',7,False)   ## imaginary eigenvalue
            self.applyDataCodeValue('dataNames',['mode','eigr','eigi'])
        elif self.analysisCode==10: # nonlinear statics
            self.addDataParameter(data,'lftsfq','f',5)   ## load step
            self.applyDataCodeValue('dataNames',['lftsfq'])
        elif self.analysisCode==11: # old geometric nonlinear statics
            self.addDataParameter(data,'lsdvmn',  'i',5)   ## load set number
            self.applyDataCodeValue('dataNames',['lsdvmn'])
        elif self.analysisCode==12: # contran ? (may appear as aCode=6)  --> straight from DMAP...grrr...
            self.addDataParameter(data,'lsdvmn',  'i',5)   ## load set number
            self.applyDataCodeValue('dataNames',['lsdvmn'])
        else:
            raise InvalidAnalysisCodeError('invalid analysisCode...analysisCode=%s' %(self.analysisCode))
        # tCode=2
        #if self.analysisCode==2: # sort2
        #    self.lsdvmn = self.getValues(data,'i',5)
        
        #print "*iSubcase=%s"%(self.iSubcase)
        #print "analysisCode=%s tableCode=%s thermal=%s" %(self.analysisCode,self.tableCode,self.thermal)
        #print self.codeInformation()

        #self.printBlock(data)
        self.readTitle()

    def readOUG_Data(self):
        #print "self.analysisCode=%s tableCode(1)=%s thermal(23)=%g" %(self.analysisCode,self.tableCode,self.thermal)
        tfsCode = [self.tableCode,self.formatCode,self.sortCode]
        #print self.dataCode
        #if self.thermal==2:
        #    self.skipOES_Element()
        #print "tfsCode=%s" %(tfsCode)
        # displacement
        if self.tableCode==1:
            self.readOUG_Data_table1()
        #if   tfsCode==[1,1,0]:  ### was on
            #self.readOUG_Data_table1_format1_sort0()
        #elif tfsCode==[1,1,1]:  ### was on
            #self.readOUG_Data_table1_format1_sort1()
        #elif tfsCode==[1,2,0]:  ### was on
            #self.readOUG_Data_table1_format2_sort0()
        #elif tfsCode==[1,2,1]:
        #    self.readOUG_Data_table1_format2_sort1()
        #elif tfsCode==[1,2,2]:
        #    self.readOUG_Data_table1_format2_sort2()
        #elif tfsCode==[1,3,0]:
        #    self.readOUG_Data_table1_format3_sort0()
        #elif tfsCode==[1,3,1]:
        #    self.readOUG_Data_table1_format3_sort1()
        #elif tfsCode==[1,3,2]:
        #    self.readOUG_Data_table1_format3_sort2()

        # modes
        elif self.tableCode==7:
            self.readOUG_Data_table7()
        #elif tfsCode==[7,1,0]: ### was on
            #self.readOUG_Data_table7_format1_sort0()
        #elif tfsCode==[7,1,1]:
        #    self.readOUG_Data_table7_format1_sort1()
        #elif tfsCode==[7,2,0]: ### was on
            #self.readOUG_Data_table7_format2_sort0()
        #elif tfsCode==[7,2,1]:
        #    self.readOUG_Data_table7_format2_sort1()
        #elif tfsCode==[7,2,2]:
        #    self.readOUG_Data_table7_format2_sort2()
        #elif tfsCode==[7,3,0]:
        #    self.readOUG_Data_table7_format3_sort0()
        #elif tfsCode==[7,3,1]:
        #    self.readOUG_Data_table7_format3_sort1()
        #elif tfsCode==[7,3,2]:
        #    self.readOUG_Data_table7_format3_sort2()

        # velocity
        elif self.tableCode==10:
            self.readOUG_Data_table10()
        #elif tfsCode==[10,1,0]:  ### was on
            #self.readOUG_Data_table10_format1_sort0()
        #elif tfsCode==[10,1,1]:
        #    self.readOUG_Data_table10_format1_sort1()
        #elif tfsCode==[10,2,0]:  ### was on
            #self.readOUG_Data_table10_format2_sort0()
        #elif tfsCode==[10,2,1]:
        #    self.readOUG_Data_table10_format2_sort1()
        #elif tfsCode==[10,2,2]:
        #    self.readOUG_Data_table10_format2_sort2()
        #elif tfsCode==[10,3,0]:
        #    self.readOUG_Data_table10_format3_sort0()
        #elif tfsCode==[10,3,1]:
        #    self.readOUG_Data_table10_format3_sort1()
        #elif tfsCode==[10,3,2]:
        #    self.readOUG_Data_table10_format3_sort2()

        # Acceleration vector
        elif self.tableCode==11:
            self.readOUG_Data_table11()
        #elif tfsCode==[11,1,0]:  ### was on
            #self.readOUG_Data_table11_format1_sort0()
        #elif tfsCode==[11,1,1]:
            #self.readOUG_Data_table11_format1_sort1()
        #elif tfsCode==[11,2,0]:  ### was on
            #self.readOUG_Data_table11_format2_sort0()

        #elif tfsCode==[11,2,1]:
        #    self.readOUG_Data_table11_format2_sort1()
        #elif tfsCode==[11,2,2]:
        #    self.readOUG_Data_table11_format2_sort2()
        #elif tfsCode==[11,3,0]:
        #    self.readOUG_Data_table11_format3_sort0()
        #elif tfsCode==[11,3,1]:
        #    self.readOUG_Data_table11_format3_sort1()
        #elif tfsCode==[11,3,2]:
        #    self.readOUG_Data_table11_format3_sort2()

        else:
            #print "***start skipping***"
            #self.log.debug('skipping approach/table/format/sortCode=%s on %s-OUG table' %(self.atfsCode,self.tableName))
            self.skipOES_Element()
            #print "***end skipping***"
            #print self.codeInformation()
            #raise NotImplementedError('bad approach/table/format/sortCode=%s on %s-OUG table' %(self.atfsCode,self.tableName))
        ###
        #print self.obj

    def readOUG_Data_table1_format1_sort0(self): # displacement

        if self.thermal==0:
            #print self.codeInformation()
            if self.analysisCode==1: # displacement
                #print "isDisplacement"
                self.createTransientObject(self.displacements,displacementObject)

            #elif self.analysisCode==5: # frequency displacement
                #print "isFrequencyDisplacement"
                #self.createTransientObject(self.displacements,eigenVectorObject)
                #raise Exception('verify...displacement5 format1')

            elif self.analysisCode==6: # transient displacement
                #print "isTransientDisplacement"
                self.createTransientObject(self.displacements,displacementObject)

            elif self.analysisCode==7: # pre-buckling displacement
                #print "isPreBucklingDisplacement"
                self.createTransientObject(self.displacements,displacementObject)

            elif self.analysisCode==8: # post-buckling eigenvector
                #print "isPostBucklingEigenvector8_1_0"
                self.createTransientObject(self.displacements,eigenVectorObject)

            #elif self.analysisCode==9: # nonlinear static eigenvector
                #print "isComplexEigenvalues"
                #self.createTransientObject(self.complexEigenvalues,eigenVectorObject)
                #raise Exception('verify...displacement9 format1')

            elif self.analysisCode==10: # nonlinear static displacement
                #print "isNonlinearStaticDisplacement"
                self.createTransientObject(self.displacements,displacementObject)
            elif self.analysisCode==11: # nonlinear geometric static displacement
                #print "isNonlinearStaticDisplacement"
                self.createTransientObject(self.displacements,displacementObject)
            else:
                #print self.codeInformation()
                #raise NotImplementedError('unsupported %s-OUG static solution...atfsCode=%s' %(self.tableName,self.atfsCode))
                pass
            ###

        elif self.thermal==1:
            if self.analysisCode==1: # temperature
                #print "isTemperature"
                self.createTransientObject(self.temperatures,temperatureObject)
            elif self.analysisCode==6: # transient temperature
                #print "isTransientTemperature"
                self.createTransientObject(self.temperatures,temperatureObject)
            elif self.analysisCode==10: # nonlinear static displacement
                #print "isNonlinearStaticTemperatures"
                self.createTransientObject(self.temperatures,temperatureObject)
                #self.createTransientObject(self.temperatures,nonlinearTemperatureObject)
            else:
                #print self.codeInformation()
                #raise NotImplementedError('unsupported %s-OUG thermal solution...atfsCode=%s' %(self.tableName,self.atfsCode))
                pass
            ###
        elif self.thermal in [2,4,8]:
            self.readThermal4()
        else:   # self.thermal>1:  ## @warning thermal>0!!!!
            #print self.codeInformation()
            #raise NotImplementedError('invalid %s-OUG thermal flag...not 0 or 1...flag=%s' %(self.tableName,self.thermal))
            pass
        ###
        
        self.readMappedScalarsOut(debug=False) # handles dtMap

        #if self.obj:
        #    self.readScalars8(debug=False)
        #else:
        #    self.skipOES_Element()
        ###
        #print self.obj
        #return

    def readThermal4(self):
        #print self.codeInformation()
        #print self.printBlock(self.data)
        n=0
        nEntries = len(self.data)//32
        for i in range(nEntries):
            eData = self.data[n:n+32]
            out = unpack('iiffffff',eData)
            nid = (out[0]-self.deviceCode)//10

            #print out
            n+=32
            #print "nid = ",nid
        #sys.exit('thermal4...')
    
    def readOUG_Data_table1(self): # new displacement
        isSort1 = self.isSort1()
        if self.numWide==8:  # real/random
            self.createTransientObject(self.displacements,displacementObject) # real
            self.OUG_RealTable()
        elif self.numWide==14:  # real/imaginary or mag/phase
            self.createTransientObject(self.displacements,complexDisplacementObject) # complex
            self.OUG_ComplexTable()
        else:
            raise NotImplementedError('only numWide=8 or 14 is allowed  numWide=%s' %(self.numWide))
        ###

    def readOUG_Data_table7(self): # new eigenvector
        isSort1 = self.isSort1()
        if self.numWide==8:  # real/random
            self.createTransientObject(self.eigenvectors,eigenVectorObject) # real
            self.OUG_RealTable()
        elif self.numWide==14:  # real/imaginary or mag/phase
            self.createTransientObject(self.eigenvectors,complexEigenVectorObject) # complex
            self.OUG_ComplexTable()
        else:
            raise NotImplementedError('only numWide=8 or 14 is allowed  numWide=%s' %(self.numWide))
        ###

    def readOUG_Data_table10(self): # new velocity
        isSort1 = self.isSort1()
        if self.numWide==8:  # real/random
            self.createTransientObject(self.velocities,velocityObject) # real
            self.OUG_RealTable()
        elif self.numWide==14:  # real/imaginary or mag/phase
            self.createTransientObject(self.velocities,complexVelocityObject) # complex
            self.OUG_ComplexTable()
        else:
            raise NotImplementedError('only numWide=8 or 14 is allowed  numWide=%s' %(self.numWide))
        ###

    def readOUG_Data_table11(self): # new acceleration
        isSort1 = self.isSort1()
        if self.numWide==8:  # real/random
            self.createTransientObject(self.accelerations,accelerationObject) # real
            self.OUG_RealTable()
        elif self.numWide==14:  # real/imaginary or mag/phase
            self.createTransientObject(self.accelerations,complexAccelerationObject) # complex
            self.OUG_ComplexTable()
        else:
            raise NotImplementedError('only numWide=8 or 14 is allowed  numWide=%s' %(self.numWide))
        ###

    def getOUG_FormatStart(self):
        """
        Returns an i or an f depending on if it's SORT2 or not.
        Also returns an extraction function that is called on the first argument
        """
        isSort1 = self.isSort1()
        if isSort1:
            #print "SORT1 - %s" %(self.ElementType(self.elementType))
            format1 = 'i' # SORT1
            extract = self.extractSort1
        else: # values from IDENT   #@todo test this...
            #print "SORT2 - %s" %(self.ElementType(self.elementType))
            if self.analysisCode in [1,2,3,4,7,8,11]:
                format1 = 'i' # SORT1
            elif self.analysisCode in [5,6,9,10,12]:
                format1 = 'f' # SORT1
            else:
                raise InvalidAnalysisCodeError('invalid analysisCode...analysisCode=%s' %(self.analysisCode))
            ###

            extract = self.extractSort2
            #eid = self.nonlinearFactor
        return (format1,extract)

    def OUG_RealTable(self):
        deviceCode = self.deviceCode
        dt = self.nonlinearFactor

        (format1,extract) = self.getOUG_FormatStart()
        format1 += 'iffffff'

        while len(self.data)>=32: # 8*4
            eData     = self.data[0:32]
            self.data = self.data[32: ]
            #print "len(data) = ",len(eData)

            out = unpack(format1, eData)
            (eid,gridType,tx,ty,tz,rx,ry,rz) = out
            eid2  = extract(eid,dt)
            #print "eType=%s" %(eType)
            
            dataIn = [eid2,gridType,tx,ty,tz,rx,ry,rz]
            #print "%s" %(self.ElementType(self.elementType)),dataIn
            #eid = self.obj.addNewEid(out)
            self.obj.add(dt,dataIn)
            #print "len(data) = ",len(self.data)
        ###
        self.handleResultsBuffer(self.OUG_RealTable)

    def OUG_ComplexTable(self):
        deviceCode = self.deviceCode
        dt = self.nonlinearFactor

        (format1,extract) = self.getOUG_FormatStart()
        format1 += 'iffffffffffff'
        #print "format1 = ",format1
        isMagnitudePhase = self.isMagnitudePhase()

        while len(self.data)>=56: # 14*4
            eData     = self.data[0:56]
            self.data = self.data[56: ]
            #print "len(data) = ",len(eData)

            out = unpack(format1, eData)
            (eid,gridType,txr,tyr,tzr,rxr,ryr,rzr,
                          txi,tyi,tzi,rxi,ryi,rzi) = out

            if isMagnitudePhase:
                (txr,txi) = polarToRealImag(txr,txi)
                (tyr,tyi) = polarToRealImag(tyr,tyi)
                (tzr,tzi) = polarToRealImag(tzr,tzi)

                (rxr,rxi) = polarToRealImag(rxr,rxi)
                (ryr,ryi) = polarToRealImag(ryr,ryi)
                (rzr,rzi) = polarToRealImag(rzr,rzi)

            eid2  = extract(eid,dt)
            #print "eType=%s" %(eType)
            
            dataIn = [eid2,gridType,txr,tyr,tzr,rxr,ryr,rzr,
                                    txi,tyi,tzi,rxi,ryi,rzi]
            #print "%s" %(self.ElementType(self.elementType)),dataIn
            #eid = self.obj.addNewEid(out)
            self.obj.add(dt,dataIn)
            #print "len(data) = ",len(self.data)
        ###
        self.handleResultsBuffer(self.OUG_ComplexTable)

