import sys
import copy
from struct import unpack

# pyNastran
from pyNastran.op2.op2Errors import *
from oee_objects import *

class OEE(object):
    """Table of energy"""

    def readTable_OEE(self):
        table3 = self.readTable_OEE_3
        table4Data = self.readOEE_Data
        self.readResultsTable(table3,table4Data)
        self.deleteAttributes_OEE()

    def deleteAttributes_OEE(self): # no thermal
        params = ['lsdvm','mode','eigr','freq','dt','lftsfq','formatCode','numWide']
        self.deleteAttributes(params)
    
    def readTable_OEE_3(self,iTable): # iTable=-3
        bufferWords = self.getMarker()
        if self.makeOp2Debug:
            self.op2Debug.write('bufferWords=%s\n' %(str(bufferWords)))
        #print "2-bufferWords = ",bufferWords,bufferWords*4,'\n'

        data = self.getData(4)
        bufferSize, = unpack('i',data)
        data = self.getData(4*50)
        #print self.printBlock(data)

        aCode = self.getBlockIntEntry(data,1)
        self.eTotal = self.parseApproachCode(data) # total energy of all elements in iSubcase/mode

        elementName = ''.join(unpack('cccccccc',data[24:32])).strip() ## element name
        if elementName.isalpha():
            self.dataCode['elementName'] = elementName

        self.addDataParameter(data,'loadSet',   'i',8, False)  ## Load set or zero
        self.addDataParameter(data,'formatCode','i',9, False)  ## format code
        self.addDataParameter(data,'numWide',   'i',10,False)  ## number of words per entry in record; @note is this needed for this table ???
        self.addDataParameter(data,'cvalres',   'i',11,False)  ## C
        self.addDataParameter(data,'setID',     'i',13,False)  ## Set identification number Number
        self.addDataParameter(data,'eigenReal', 'i',14,False)  ## Natural eigenvalue - real part
        self.addDataParameter(data,'eigenImag', 'i',15,False)  ## Natural eigenvalue - imaginary part
        self.addDataParameter(data,'freq',      'f',16,False)  ## Natural frequency
        self.addDataParameter(data,'etotpos',   'f',18)        ## Total positive energy
        self.addDataParameter(data,'etotneg',   'f',19,False)  ## Total negative energy

        if not self.isSort1():
            raise NotImplementedError('sort2...')

        #self.printBlock(data) # on
        if self.analysisCode==1:   # statics / displacement / heat flux
            #del self.dataCode['nonlinearFactor']
            self.applyDataCodeValue('dataNames',['lsdvmn'])
            self.setNullNonlinearFactor()
        elif self.analysisCode==2: # real eigenvalues
            self.addDataParameter(data,'mode','i',5)   ## mode number
            self.applyDataCodeValue('dataNames',['mode'])
            #print "mode(5)=%s eigr(6)=%s modeCycle(7)=%s" %(self.mode,self.eigr,self.modeCycle)
        #elif self.analysisCode==3: # differential stiffness
            #self.lsdvmn = self.getValues(data,'i',5) ## load set number
            #self.dataCode['lsdvmn'] = self.lsdvmn
        #elif self.analysisCode==4: # differential stiffness
            #self.lsdvmn = self.getValues(data,'i',5) ## load set number
        elif self.analysisCode==5:   # frequency
            self.addDataParameter(data,'freq2','f',5)   ## frequency
            self.applyDataCodeValue('dataNames',['freq2'])
        elif self.analysisCode==6: # transient
            self.addDataParameter(data,'time','f',5)   ## time step
            self.applyDataCodeValue('dataNames',['time'])
        #elif self.analysisCode==7: # pre-buckling
            #self.applyDataCodeValue('dataNames',['lsdvmn'])
        elif self.analysisCode==8: # post-buckling
            self.addDataParameter(data,'mode','i',5)   ## mode number
            self.applyDataCodeValue('dataNames',['mode'])
        elif self.analysisCode==9: # complex eigenvalues
            self.addDataParameter(data,'mode','i',5)   ## mode number
            self.applyDataCodeValue('dataNames',['mode'])
        elif self.analysisCode==10: # nonlinear statics
            self.addDataParameter(data,'loadFactor','f',5)   ## load factor
            self.applyDataCodeValue('dataNames',['loadFactor'])
        #elif self.analysisCode==11: # old geometric nonlinear statics
            #self.applyDataCodeValue('dataNames',['lsdvmn'])
        elif self.analysisCode==12: # contran ? (may appear as aCode=6)  --> straight from DMAP...grrr...
            self.addDataParameter(data,'time','f',5)   ## time step
            self.applyDataCodeValue('dataNames',['time'])
        else:
            raise InvalidAnalysisCodeError('invalid analysisCode...analysisCode=%s' %(self.analysisCode))
        ###
        
        #print "*iSubcase=%s elementName=|%s|"%(self.iSubcase,self.elementName)
        #print "analysisCode=%s tableCode=%s" %(self.analysisCode,self.tableCode)
        #print self.codeInformation()

        #self.printBlock(data)
        self.readTitle()

    def readOEE_Data(self):
        #print "self.analysisCode=%s tableCode(1)=%s" %(self.analysisCode,self.tableCode)
        tfsCode = [self.tableCode,self.formatCode,self.sortCode]
        
        if self.tableCode==18:
            self.readStrainEnergy_table18()
        else:
            #self.skipOES_Element()
            raise NotImplementedError('unsupported OEE static solution...aftsCode=%s' %(self.atfsCode))
        ###
        #print str(self.obj)

    
    def readStrainEnergy_table18(self):
        self.createTransientObject(self.strainEnergy,StrainEnergyObject)
        if self.numWide==4:
            self.OEE_Strain4()
        elif self.numWide==5:
            self.OEE_Strain5()
        else:   
            raise NotImplementedError()
        #self.readMappedScalarsOut(debug=False) # handles dtMap, not correct...

    def OEE_Strain4(self):
        deviceCode = self.deviceCode
        dt = self.nonlinearFactor

        (format1,extract) = self.getOUG_FormatStart()  ## @todo change to OEE
        format1 += 'fff'

        while len(self.data)>=16: # 4*4
            eData     = self.data[0:16]
            self.data = self.data[16: ]
            #print "len(data) = ",len(eData)

            out = unpack(format1, eData)
            (eid,energy,percent,density) = out
            eid2  = extract(eid,dt)
            #print "eType=%s" %(eType)
            
            dataIn = [eid2,energy,percent,density]
            #print "%s" %(self.ElementType(self.elementType)),dataIn
            #eid = self.obj.addNewEid(out)
            self.obj.add(dt,dataIn)
            #print "len(data) = ",len(self.data)
        ###
        self.handleResultsBuffer(self.OEE_Strain4)
        #print self.strainEnergy

    def OEE_Strain5(self):
        deviceCode = self.deviceCode
        dt = self.nonlinearFactor

        #(format1,extract) = self.getOUG_FormatStart()  ## @todo change to OEE
        format1 = 'ccccccccfff'

        while len(self.data)>=16: # 5*4
            eData     = self.data[0:20]
            self.data = self.data[20: ]
            #print "len(data) = ",len(eData)

            out = unpack(format1, eData)
            (a,b,c,d,e,f,g,h,energy,percent,density) = out
            #print "out = ",out
            word = a+b+c+d+e+f+g+h
            word = word.strip()
            #print "eType=%s" %(eType)
            
            dataIn = [word,energy,percent,density]
            #print "%s" %(self.ElementType(self.elementType)),dataIn
            #eid = self.obj.addNewEid(out)
            self.obj.add(dt,dataIn)
            #print "len(data) = ",len(self.data)
        ###
        self.handleResultsBuffer(self.OEE_Strain5)
        #print self.strainEnergy
