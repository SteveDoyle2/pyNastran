import sys
import copy
from struct import unpack
from pyNastran.op2.op2Errors import *
from oee_objects import *

# pyNastran
#from pyNastran.op2.resultObjects.ougv1_Objects import (
#     temperatureObject,displacementObject,  # analysisCode=1, sortCode=0
#     eigenVectorObject,                     # analysisCode=2, sortCode=0
#     fluxObject,                            # analysisCode=1, sortCode=3
#     nonlinearTemperatureObject,            # analysisCode=10,sortCode=0
#     )

class OEE(object):
    """Table of energy"""

    def readTable_OEE1(self):
        self.tableName = 'OEE'
        #self.staticEnergy = {} # aCode=1 tCode=18 fCode=1 sortCode=0

        table3 = self.readTable_OEE1_3
        table4Data = self.readOEE1_Data
        self.readResultsTable(table3,table4Data)
        self.deleteAttributes_OEE()
        #sys.exit('end of oee')

    def deleteAttributes_OEE(self): # no thermal
        params = ['lsdvm','mode','eigr','freq','dt','lftsfq','formatCode','numWide']
        self.deleteAttributes(params)
    
    def readTable_OEE1_3(self,iTable): # iTable=-3
        bufferWords = self.getMarker()
        if self.makeOp2Debug:
            self.op2Debug.write('bufferWords=%s\n' %(str(bufferWords)))
        #print "2-bufferWords = ",bufferWords,bufferWords*4,'\n'

        data = self.getData(4)
        bufferSize, = unpack('i',data)
        data = self.getData(4*50)
        #print self.printBlock(data)

        
        
        aCode = self.getBlockIntEntry(data,1)
        #print "aCode = ",aCode
        self.eTotal = self.parseApproachCode(data) # total energy of all elements in iSubcase/mode

        #self.dataCode = {'analysisCode': self.analysisCode,'deviceCode':self.deviceCode,
        #                 'loadSet':self.loadSet,'formatCode':self.formatCode,
        #                 'numWide': self.numWide,'cvalres':self.cvalres,
        #                 'esubt': self.esubt,'setID':self.setID,'eigenReal':self.eigenReal,'eigenImag':self.eigenImag,
        #                 'freq':self.freq,'etotpos':self.etotpos,'etotneg':self.etotneg}

        self.dataCode = {'analysisCode': self.analysisCode,'deviceCode':self.deviceCode,}

        self.elementName  = ''.join(unpack('cccccccc',data[24:32])) ## element name
        #elementName1      = self.getValues(data,'cccc',6)  
        #elementName2      = self.getValues(data,'cccc',7)  ## element name
        #self.elementName  = elementName1+elementName2

        self.addDataParameter(data,'loadSet',   'i',8)   ## Load set or zero
        self.addDataParameter(data,'formatCode','i',9)   ## format code
        self.addDataParameter(data,'numWide',   'i',10)  ## number of words per entry in record; @note is this needed for this table ???
        self.addDataParameter(data,'cvalres',   'i',11)  ## C

        self.addDataParameter(data,'setID',    'i',13)  ## Set identification number Number
        self.addDataParameter(data,'eigenReal','i',14)  ## Natural eigenvalue - real part
        self.addDataParameter(data,'eigenImag','i',15)  ## Natural eigenvalue - imaginary part
        self.addDataParameter(data,'freq',     'f',16)  ## Natural frequency
        self.addDataParameter(data,'etotpos',  'f',18)  ## Total positive energy
        self.addDataParameter(data,'etotneg',  'f',19)  ## Total negative energy


        #self.printBlock(data) # on
        if self.analysisCode==1:   # statics / displacement / heat flux
            pass
        elif self.analysisCode==2: # real eigenvalues
            self.addDataParameter(data,'mode',     'i',5)   ## mode number
            self.nonlinearFactor = self.mode
            #print "mode(5)=%s" %(self.mode)
        elif self.analysisCode==3: # differential stiffness
            pass
        elif self.analysisCode==4: # differential stiffness
            pass
        elif self.analysisCode==5:   # frequency
            self.addDataParameter(data,'freq2','f',5)   ## frequency
            self.nonlinearFactor = self.freq2 ## why are there 2 values of freq?

        elif self.analysisCode==6: # transient
            self.time = self.getValues(data,'f',5) ## time step
            self.nonlinearFactor = self.time
            print "time(5)=%s" %(self.time)
        elif self.analysisCode==7: # pre-buckling
            pass
        elif self.analysisCode==8: # post-buckling
            self.addDataParameter(data,'mode','i',5)   ## mode number
            self.nonlinearFactor = self.mode
            print "mode(5)=%s" %(self.mode)
        elif self.analysisCode==9: # complex eigenvalues
            self.addDataParameter(data,'mode','i',5)   ## mode number
            self.nonlinearFactor = self.mode
            #print "mode(5)=%s" %(self.mode)
        elif self.analysisCode==10: # nonlinear statics
            self.loadFactor = self.getValues(data,'f',5) ## load factor
            self.nonlinearFactor = self.loadFactor
            print "loadFactor(5) = %s" %(self.loadFactor)
        elif self.analysisCode==11: # old geometric nonlinear statics
            pass
        elif self.analysisCode==12: # contran ? (may appear as aCode=6)  --> straight from DMAP...grrr...
            self.time = self.getValues(data,'f',5) ## time step
            self.nonlinearFactor = self.time
            print "time(5)=%s" %(self.time)
        else:
            raise RuntimeError('invalid analysis code...analysisCode=%s' %(self.analysisCode))
        ###
        
        #print "*iSubcase=%s elementName=|%s|"%(self.iSubcase,self.elementName)
        #print "analysisCode=%s tableCode=%s" %(self.analysisCode,self.tableCode)
        if self.numWide==5:
            print self.codeInformation()

        #self.printBlock(data)
        self.readTitle()

    def readOEE1_Data(self):
        print "self.analysisCode=%s tableCode(1)=%s" %(self.analysisCode,self.tableCode)
        tfsCode       =                   [self.tableCode,self.formatCode,self.sortCode]
        
        if tfsCode==[18,1,0]:
            self.readStrainEnergy_table18_format1_sort0()
        #elif fsCode==[18,1,1]:
        #    self.readOEE1_Data_format1_sort1()
        #elif fsCode==[18,2,1]:
        #    self.readOEE1_Data_format2_sort1()
        else:
            self.skipOES_Element()
            raise Exception('unsupported OEE1 static solution...aftsCode=%s' %(self.atfsCode))
        ###
        #print str(self.obj)

    def readStrainEnergy_table18_format1_sort0(self):
        assert self.tableCode==18 # Strain Energy
        assert self.formatCode==1 # Real
        assert self.sortCode==0   # Real

        if self.analysisCode==1: # displacement
            print "isStrainEnergy"
            self.obj = StrainEnergyObject(self.dataCode,self.iSubcase,self.nonlinearFactor)
            self.strainEnergy[self.iSubcase] = self.obj
        elif self.analysisCode==2: # buckling modes
            #print "isBucklingStrainEnergy"
            self.createTransientObject(self.strainEnergy,StrainEnergyObject)
        elif self.analysisCode==5: # freq
            print "isFreqStrainEnergy"
            self.createTransientObject(self.strainEnergy,StrainEnergyObject)
        elif self.analysisCode==6: # transient
            print "isTransientStrainEnergy"
            self.createTransientObject(self.strainEnergy,StrainEnergyObject)
        elif self.analysisCode==9: # nonlinear static eigenvector
            print "isComplexStrainEnergy"
            self.createTransientObject(self.strainEnergy,StrainEnergyObject)
        elif self.analysisCode==10: # nonlinear statics
            print "isNonlinearStrainEnergy"
            self.createTransientObject(self.strainEnergy,StrainEnergyObject)
        else:
            raise Exception('bad analysis/table/format/sortCode=%s on OEE table' %(self.atfsCode))
        ###
        self.readScalars4o(debug=True)

