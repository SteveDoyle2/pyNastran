import sys
import copy
from numpy import array
from struct import unpack

# pyNastran
from ogf_Objects import gridPointForcesObject,complexGridPointForcesObject

class OGF(object):
    """Table of Grid Point Forces"""

    def readTable_OGF(self):
        table3 = self.readTable_OGF_3
        table4Data = self.readOGF_Data

        self.dtMap = {}
        self.readResultsTable(table3,table4Data)
        del self.dtMap
        self.deleteAttributes_OGF()

    def deleteAttributes_OGF(self):
        params = ['formatCode','appCode','numWide','value1','value2',]
        self.deleteAttributes(params)
    
    #def addDataParameter(self,data,Name,Type,FieldNum,applyNonlinearFactor=True):
    #    #self.mode = self.getValues(data,'i',5) ## mode number
    #    value = self.getValues(data,Type,FieldNum)
    #    setattr(self,Name,value)
    #    self.dataCode[Name] = value
    #    
    #    if applyNonlinearFactor:
    #        self.nonlinearFactor = value
    #        self.dataCode['nonlinearFactor'] = value
    #        self.dataCode['name'] = Name
    
    def applyDataCodeValue(self,Name,value):
        self.dataCode[Name] = value
        
    def readTable_OGF_3(self,iTable): # iTable=-3
        bufferWords = self.getMarker()
        if self.makeOp2Debug:
            self.op2Debug.write('bufferWords=%s\n' %(str(bufferWords)))
        #print "2-bufferWords = ",bufferWords,bufferWords*4,'\n'

        data = self.getData(4)
        bufferSize, = unpack('i',data)
        data = self.getData(4*50)
        #print self.printBlock(data)
        
        (three) = self.parseApproachCode(data)

        self.addDataParameter(data,'formatCode',  'i',9,False)   ## format code
        self.addDataParameter(data,'appCode',     'i',9,False)   ## approach code ???
        self.addDataParameter(data,'numWide',     'i',10,False)  ## number of words per entry in record; @note is this needed for this table ???
        self.addDataParameter(data,'value1',      'i',11,False)  ## Data Value 1
        self.addDataParameter(data,'value2',      'f',12,False)  ## Data Value 2
        
        #self.printBlock(data) # on
        ## assuming tCode=1
        if self.analysisCode==1:   # statics / displacement / heat flux
            self.extractDt = self.extractInt
        elif self.analysisCode==2: # real eigenvalues
            self.extractDt = self.extractInt
        elif self.analysisCode==3: # differential stiffness
            self.extractDt = self.extractInt
        elif self.analysisCode==4: # differential stiffness
            self.extractDt = self.extractInt
        elif self.analysisCode==5:   # frequency
            self.extractDt = self.extractFloat
        elif self.analysisCode==6: # transient
            self.extractDt = self.extractFloat
        elif self.analysisCode==7: # pre-buckling
            self.extractDt = self.extractInt
        elif self.analysisCode==8: # post-buckling
            self.extractDt = self.extractInt
        elif self.analysisCode==9: # complex eigenvalues
            self.extractDt = self.extractInt
        elif self.analysisCode==10: # nonlinear statics
            self.extractDt = self.extractFloat
        elif self.analysisCode==11: # old geometric nonlinear statics
            self.extractDt = self.extractInt
        elif self.analysisCode==12: # contran ? (may appear as aCode=6)  --> straight from DMAP...grrr...
            self.extractDt = self.extractInt
        else:
            raise InvalidAnalysisCodeError('invalid analysisCode...analysisCode=%s' %(self.analysisCode))

        
        #print "*iSubcase=%s"%(self.iSubcase)
        #print "analysisCode=%s tableCode=%s thermal=%s" %(self.analysisCode,self.tableCode,self.thermal)

        #self.printBlock(data)
        self.readTitle()

    def extractDt(self):
        pass
    def extractInt(self):
        return 'i',self.scaleEid
    def extractFloat(self):
        return 'f',self.scaleDt
    def scaleEid(self,eid):
        return (eid-self.deviceCode)//10
    def scaleDt(self,dt):
        return dt

    def readOGF_Data(self):
        #print "self.analysisCode=%s tableCode(1)=%s thermal(23)=%g" %(self.analysisCode,self.tableCode,self.thermal)
        tfsCode = [self.tableCode,self.formatCode,self.sortCode]
        #print self.dataCode
        #if self.thermal==2:
        #    self.skipOES_Element()
        #print "tfsCode=%s" %(tfsCode)
        # displacement
        if   tfsCode==[19,1,0]:
            self.readOGF_Data_table19_format1_sort0()
        elif tfsCode==[19,1,1]:
            self.readOGF_Data_table19_format1_sort1()
        elif tfsCode==[19,2,0]:
            self.readOGF_Data_table19_format2_sort0()
        #elif tfsCode==[19,2,1]:
        #    self.readOGF_Data_table19_format2_sort1()
        #elif tfsCode==[19,2,2]:
        #    self.readOGF_Data_table19_format2_sort2()
        #elif tfsCode==[19,3,0]:
        #    self.readOGF_Data_table19_format3_sort0()
        #elif tfsCode==[19,3,1]:
        #    self.readOGF_Data_table19_format3_sort1()
        #elif tfsCode==[19,3,2]:
        #    self.readOGF_Data_table19_format3_sort2()

        else:
            #print "***start skipping***"
            #self.log.debug('skipping approach/table/format/sortCode=%s on %s-OGF table' %(self.tableName,self.atfsCode))
            #print self.codeInformation()
            #self.skipOES_Element()
            #print "***end skipping***"
            raise NotImplementedError('bad approach/table/format/sortCode=%s on %s-OGF table' %(self.atfsCode,self.tableName,))
        ###
        #print self.obj

    def readOGF_Data_table19_format1_sort0(self): # displacement
        print self.codeInformation()
        if self.numWide==10:
            self.readOGF_numWide10()
        elif self.numWide==16:
            self.readOGF_numWide16()
        else:
            raise NotImplementedError('%s-OGF only supports numWide=10,16' %(self.tableName))

        #sys.exit('stopping in OGF')        
        
        #readCase = True
        #if self.iSubcase in self.expectedTimes and len(self.expectedTimes[self.iSubcase])>0:
        #    readCase = self.updateDtMap()
        
        #if self.obj and readCase:
        #    self.readScalarsOut(debug=False)
        #else:
        #    self.skipOES_Element()
        ###
        #if self.obj:
        #    self.readScalars8(debug=False)
        #else:
        #    self.skipOES_Element()
        ###
        #print self.obj
        #return
    
    def readOGF_numWide10(self):
        (eKey,scaleValue) = self.extractDt()
        Format = eKey+'issssssssffffff'
        while len(self.data)>=40:
            eData     = self.data[0:4*10]
            self.data = self.data[4*10: ]
            out = unpack(Format,eData)
            (ekey,eid,a,b,c,d,e,f,g,h,f1,f2,f3,m1,m2,m3) = out
            ekey = scaleValue(ekey)
            elemName = a+b+c+d+e+f+g+h
            elemName = elemName.strip()
            data = (eid,elemName,f1,f2,f3,m1,m2,m3)
            self.obj.add(data)
            #print "eid/dt/freq=%s eid=%-6s eName=%-8s f1=%g f2=%g f3=%g m1=%g m2=%g m3=%g" %(ekey,eid,elemName,f1,f2,f3,m1,m2,m3)
        #print len(self.data)
        self.handleResultsBuffer(self.readOGF_numWide10)

    def readOGF_numWide16(self):
        (eKey,scaleValue) = self.extractDt()
        Format = eKey+'issssssssffffffffffff'
        while len(self.data)>=64:
            eData     = self.data[0:4*16]
            self.data = self.data[4*16: ]
            out = unpack(Format,eData)
            (ekey,eid,a,b,c,d,e,f,g,h,f1r,f2r,f3r,m1r,m2r,m3r,f1i,f2i,f3i,m1i,m2i,m3i) = out
            ekey = scaleValue(ekey)
            f1i = f1i*1j
            f1i = f2i*1j
            f1i = f3i*1j
            m1i = m1i*1j
            m2i = m2i*1j
            m3i = m3i*1j
            elemName = a+b+c+d+e+f+g+h
            elemName = elemName.strip()
            print "eid/dt/freq=%s eid=%-6s eName=%-8s f1=%s f2=%s f3=%s m1=%s m2=%s m3=%s" %(ekey,eid,elemName,f1r+f1i,f2r+f2i,f3r+f3i,m1r+m1i,m2r+m2i,m3r+m3i)
        self.handleResultsBuffer(self.readOGF_numWide16)

    def readThermal4(self):
        print self.codeInformation()
        print self.printBlock(self.data)
        n=0
        nEntries = len(self.data)//32
        for i in range(nEntries):
            eData = self.data[n:n+32]
            out = unpack('iiffffff',eData)
            nid = (out[0]-self.deviceCode)//10

            print out
            n+=32
            print "nid = ",nid
        #sys.exit('thermal4...')

    def readOGF_Data_table7_format1_sort0(self):  # modes
        #assert self.formatCode==1 # Real
        #assert self.sortCode==0   # Real
        
        if self.thermal==0:
            #print self.codeInformation()
            if self.analysisCode==2: # nonlinear static eigenvector
                #print "isEigenvector2"
                self.createTransientObject(self.eigenvectors,eigenVectorObject)
            elif self.analysisCode==8: # post-buckling eigenvector
                #print "isPostBucklingEigenvector8"
                self.createTransientObject(self.eigenvectors,realEigenVectorObject)
            else:
                raise Exception('unsupported %s-OGF static solution...atfsCode=%s' %(self.tableName,self.atfsCode))
                pass
            ###
        elif self.thermal==1:
            raise Exception('unsupported %s-OGF thermal solution...atfsCode=%s' %(self.tableName,self.atfsCode))
            pass
        else:
            raise Exception('invalid %s-OGF thermal flag...not 0 or 1...flag=%s' %(self.tableName,self.thermal))
            pass
        ###
        if self.obj:
            self.readScalars8(debug=False)
        else:
            self.skipOES_Element()
        #if self.analysisCode not in [2,8]:
        #    raise Exception('check_format1...')
        ###

    def readOGF_Data_table1_format1_sort1(self): # displacement
        #assert self.formatCode==1 # Real
        #assert self.sortCode==1   # Real/Imaginary
        if self.thermal==0:
            if self.analysisCode==5: # complex displacements (real/imaginary)
                #print "isComplexDisplacement"
                self.createTransientObject(self.displacements,complexDisplacementObject)
            #elif self.analysisCode==7: # pre-buckling displacement
                #print "isPreBucklingDisplacement"
                #self.createTransientObject(self.displacements,displacementObject)
            #elif self.analysisCode==9: # nonlinear static eigenvector
                #print "isComplexEigenvalues"
                #self.createTransientObject(self.displacements,eigenVectorObject)
            #elif self.analysisCode==11: # Geometric nonlinear statics
                #print "isNonlinearStaticDisplacement"
                #self.createTransientObject(self.displacements,displacementObject)
            else:
                raise NotImplementedError('unsupported %s-OGF static table1_format1_sort1 solution...atfsCode=%s' %(self.tableName,self.atfsCode))
            ###
        else:
            raise Exception('invalid thermal flag...not 0 or 1...flag=%s' %(self.thermal))
        ###
        #print "objName = ",self.obj.name()
        self.readScalars14(debug=False)
        #print "---OBJ---"
        #print self.obj
        #raise Exception('format1_sort1')
        #return

    def readOGF_Data_table7_format1_sort1(self): # modes
        #assert self.formatCode==1 # Real
        #assert self.sortCode==1   # Real/Imaginary
        #print self.codeInformation()
        if self.thermal==0:
            #if self.analysisCode==5: # frequency displacement
                #print "isFrequencyDisplacement"
                #self.createTransientObject(self.freqDisplacements,eigenVectorObject)
            #elif self.analysisCode==7: # pre-buckling displacement
                #print "isPreBucklingDisplacement"
                #self.createTransientObject(self.displacements,displacementObject)
            if self.analysisCode==9: # nonlinear static eigenvector
                #print "isComplexEigenvalues"
                self.createTransientObject(self.eigenvectors,complexEigenVectorObject)
            #elif self.analysisCode==11: # Geometric nonlinear statics
                #print "isNonlinearStaticDisplacement"
                #self.createTransientObject(self.displacements,displacementObject)
            else:
                raise NotImplementedError('unsupported %s-OGF static table7_format1_sort1 solution...atfsCode=%s' %(self.tableName,self.atfsCode))
            ###
        else:
            raise Exception('invalid thermal flag...not 0 or 1...flag=%s' %(self.thermal))
        ###
        #print "objName = ",self.obj.name()
        self.readScalars14()
        #print self.obj
        #return

    def readOGF_Data_table1_format2_sort0(self): # displacement
        if self.thermal==0:
            if self.analysisCode==6: # transient displacement
                #print "isTransientDisplacement"
                self.createTransientObject(self.displacements,displacementObject)
            else:
                raise NotImplementedError('unsupported %s-OGF static solution...atfsCode=%s' %(self.tableName,self.atfsCode))
            ###
        elif self.thermal==1:
            raise Exception('unsupported %s-OGF thermal solution...atfsCode=%s' %(self.tableName,self.atfsCode))
        else:
            raise Exception('invalid thermal flag...not 0 or 1...flag=%s' %(self.thermal))
        ###
        #self.readScalars8()
        self.readScalarsOut(debug=False)

    def readOGF_Data_table1_format2_sort1(self): # displacement
        #assert self.formatCode==2 # Real/Imaginary
        #assert self.sortCode==1   # Real/Imaginary
        #print self.codeInformation()
        if self.thermal==0:
            if self.analysisCode==5: # frequency displacement
                #print "isFrequencyDisplacement"
                self.createTransientObject(self.displacements,complexDisplacementObject)
                self.readScalarsF14()
            elif self.analysisCode==7: # pre-buckling displacement
                #print "isPreBucklingDisplacement"
                self.createTransientObject(self.displacements,displacementObject)
                self.readScalarsF14()
            #elif self.analysisCode==9 and self.sortCode==1: # nonlinear static eigenvector
                #print "isComplexEigenvalues"
                #self.createTransientObject(self.complexEigenvalues,eigenVectorObject)
                #self.readScalars8()
            else:
                raise NotImplementedError('unsupported %s-OGF static solution...atfsCode=%s' %(self.tableName,self.atfsCode))
            ###
        elif self.thermal==1:
            raise Exception('unsupported %s-OGF thermal solution...atfsCode=%s' %(self.tableName,self.atfsCode))
        else:
            raise Exception('invalid thermal flag...not 0 or 1...flag=%s' %(self.thermal))
        ###
        #print self.obj
        #return

    def readOGF_Data_table1_format3_sort0(self): # displacement
        #assert self.formatCode==3 # Magnitude/Phase
        #assert self.sortCode==0   # Real
        #print self.codeInformation()
        if self.thermal==0:
            raise Exception('unsupported OGF static solution...atfsCode=%s' %(self.atfsCode))
            if self.analysisCode==1: # displacement
                #print "isDisplacement"
                self.createTransientObject(self.displacements,displacementObject)
            elif self.analysisCode==2: # nonlinear static eigenvector
                #print "isEigenvector_3_0"
                self.createTransientObject(self.eigenvectors,eigenVectorObject)
            elif self.analysisCode==5: # frequency displacement
                #print "isFreqDisplacement"
                self.createTransientObject(self.freqDisplacements,displacementObject)
            elif self.analysisCode==6: # transient displacement
                #print "isTransientDisplacement"
                self.createTransientObject(self.displacements,displacementObject)
            elif self.analysisCode==7: # pre-buckling displacement
                #print "isPreBucklingDisplacement"
                self.createTransientObject(self.preBucklingDisplacements,displacementObject)
            elif self.analysisCode==10: # transient velocity
                #print "isTransientVelocity"
                self.createTransientObject(self.velocities,displacementObject,self.dt)
            else:
                raise NotImplementedError('unsupported %s-OGF static solution...atfsCode=%s' %(self.tableName,self.atfsCode))
            ###
        elif self.thermal==1:
            raise Exception('unsupported %s-OGF thermal solution...atfsCode=%s' %(self.tableName,self.atfsCode))
        else:
            raise Exception('invalid thermal flag...not 0 or 1...flag=%s' %(self.thermal))
        ###
        #self.readScalars8()
        self.readScalarsOut(debug=False)

    def readOGF_Data_table1_format3_sort1(self): # displacemnt
        #assert self.formatCode==3 # Magnitude/Phase
        #assert self.sortCode==1   # Imaginary
        if self.thermal==0:
            if self.analysisCode==5: # complex frequency displacement
                #print "isComplexFreqDisplacement"
                self.createTransientObject(self.displacements,complexDisplacementObject)
            #elif self.analysisCode==7: # pre-buckling displacement
                #print "isComlexPreBucklingDisplacement"
                #self.createTransientObject(self.complexDisplacements,displacementObject)
            #elif self.analysisCode==9: # nonlinear static eigenvector
                #print "isComplexEigenvalues"
                #self.createTransientObject(self.complexEigenvalues,eigenVectorObject)
            else:
                raise NotImplementedError('unsupported %s-OGF static solution...atfsCode=%s' %(self.tableName,self.atfsCode))
            ###
        elif self.thermal==1:
            raise Exception('unsupported %s-OGF thermal solution...atfsCode=%s' %(self.tableName,self.atfsCode))
        else:
            raise Exception('invalid thermal flag...not 0 or 1...flag=%s' %(self.thermal))
        ###
        self.readScalars14()

