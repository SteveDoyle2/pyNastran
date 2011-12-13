import sys
import copy
from struct import unpack
from pyNastran.op2.resultObjects.op2_Objects import scalarObject

# pyNastran
#from pyNastran.op2.resultObjects.ougv1_Objects import (
#     temperatureObject,displacementObject,  # approachCode=1, sortCode=0
#     eigenVectorObject,                     # approachCode=2, sortCode=0
#     fluxObject,                            # approachCode=1, sortCode=3
#     nonlinearTemperatureObject,            # approachCode=10,sortCode=0
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
        print "aCode = ",aCode
        self.parseApproachCode(data)

        self.eTotal       = self.getValues(data,'f',3)  ## total energy of all elements in iSubcase/mode
        self.elementName  = ''.join(unpack('cccccccc',data[24:32])) ## element name
        #elementName1      = self.getValues(data,'cccc',6)  
        #elementName2      = self.getValues(data,'cccc',7)  ## element name
        #self.elementName  = elementName1+elementName2

        self.loadSet      = self.getValues(data,'i',8)  ## Load set or zero
        self.formatCode   = self.getValues(data,'i',9)  ## format code
        self.numWide      = self.getValues(data,'i',10) ## number of words per entry in record

        self.cvalres      = self.getValues(data,'i',11) ## C
        self.esubt        = self.getValues(data,'i',12) ## Subtotal of Strain Energy in the Set identification number
        self.setID        = self.getValues(data,'i',13) ## Set identification number Number
        self.eigenReal    = self.getValues(data,'f',14) ## Natural eigenvalue - real part
        self.eigenImag    = self.getValues(data,'f',15) ## Natural eigenvalue - imaginary part

        self.freq         = self.getValues(data,'f',16) ## Natural frequency
        self.etotpos      = self.getValues(data,'f',18) ## Total positive energy
        self.etotneg      = self.getValues(data,'f',19) ## Total negative energy
        self.nonlinearFactor = None

        self.dataCode = {'analysisCode': self.approachCode,'deviceCode':self.deviceCode,
                         'loadSet':self.loadSet,'formatCode':self.formatCode,
                         'numWide': self.numWide,'cvalres':self.cvalres,
                         'esubt': self.esubt,'setID':self.setID,'eigenReal':self.eigenReal,'eigenImag':self.eigenImag,
                         'freq':self.freq,'etotpos':self.etotpos,'etotneg':self.etotneg}

        #self.printBlock(data) # on
        if self.approachCode==1:   # statics / displacement / heat flux
            pass
        elif self.approachCode==2: # real eigenvalues
            self.mode      = self.getValues(data,'i',5) ## mode number
            self.nonlinearFactor = self.mode
            print "mode(5)=%s" %(self.mode)
        elif self.approachCode==3: # differential stiffness
            pass
        elif self.approachCode==4: # differential stiffness
            pass
        elif self.approachCode==5:   # frequency
            self.freq2 = self.getValues(data,'f',5) ## frequency
            self.nonlinearFactor = self.freq2 ## why are there 2 values of freq?

        elif self.approachCode==6: # transient
            self.time = self.getValues(data,'f',5) ## time step
            self.nonlinearFactor = self.time
            print "time(5)=%s" %(self.time)
        elif self.approachCode==7: # pre-buckling
            pass
        elif self.approachCode==8: # post-buckling
            self.mode = self.getValues(data,'i',5) ## mode number
            self.nonlinearFactor = self.mode
            print "mode(5)=%s" %(self.mode)
        elif self.approachCode==9: # complex eigenvalues
            self.mode   = self.getValues(data,'i',5) ## mode number
            self.nonlinearFactor = self.mode
            print "mode(5)=%s" %(self.mode)
        elif self.approachCode==10: # nonlinear statics
            self.loadFactor = self.getValues(data,'f',5) ## load factor
            self.nonlinearFactor = self.loadFactor
            print "loadFactor(5) = %s" %(self.loadFactor)
        elif self.approachCode==11: # old geometric nonlinear statics
            pass
        elif self.approachCode==12: # contran ? (may appear as aCode=6)  --> straight from DMAP...grrr...
            self.time = self.getValues(data,'f',5) ## time step
            self.nonlinearFactor = self.time
            print "time(5)=%s" %(self.time)
        else:
            raise RuntimeError('invalid approach code...approachCode=%s' %(self.approachCode))
        ###
        
        print "*iSubcase=%s elementName=|%s|"%(self.iSubcase,self.elementName)
        print "approachCode=%s tableCode=%s" %(self.approachCode,self.tableCode)
        print self.codeInformation()

        #self.printBlock(data)
        self.readTitle()

    def readOEE1_Data(self):
        print "self.approachCode=%s tableCode(1)=%s" %(self.approachCode,self.tableCode)
        self.atfsCode = [self.approachCode,self.tableCode,self.formatCode,self.sortCode]
        tfsCode       =                   [self.tableCode,self.formatCode,self.sortCode]
        
        if tfsCode==[18,1,0]:
            self.readStrainEnergy_format1_sort0()
        #elif fsCode==[18,1,1]:
        #    self.readOEE1_Data_format1_sort1()
        #elif fsCode==[18,2,1]:
        #    self.readOEE1_Data_format2_sort1()
        else:
            self.skipOES_Element(None)
            raise Exception('unsupported OEE1 static solution...aftsCode=%s' %(self.atfsCode))
        ###
        #print str(self.obj)

    def readStrainEnergy_format1_sort0(self):
        assert self.tableCode==18 # Strain Energy
        assert self.formatCode==1 # Real
        assert self.sortCode==0   # Real

        if self.approachCode==1: # displacement
            print "isStrainEnergy"
            self.obj = StrainEnergyObject(self.dataCode,self.iSubcase,self.nonlinearFactor)
            self.strainEnergy[self.iSubcase] = self.obj
        elif self.approachCode==2: # buckling modes
            print "isBucklingStrainEnergy"
            self.createTransientObject(self.strainEnergy,StrainEnergyObject)
        elif self.approachCode==5: # freq
            print "isFreqStrainEnergy"
            self.createTransientObject(self.strainEnergy,StrainEnergyObject)
        elif self.approachCode==6: # transient
            print "isTransientStrainEnergy"
            self.createTransientObject(self.strainEnergy,StrainEnergyObject)
        elif self.approachCode==9: # nonlinear static eigenvector
            print "isComplexStrainEnergy"
            self.createTransientObject(self.strainEnergy,StrainEnergyObject)
        elif self.approachCode==10: # nonlinear statics
            print "isNonlinearStrainEnergy"
            self.createTransientObject(self.strainEnergy,StrainEnergyObject)
        else:
            raise Exception('bad approach/table/format/sortCode=%s on OEE table' %(self.atfsCode))
        ###
        self.readScalars4(self.obj)

class StrainEnergyObject(scalarObject):
    def __init__(self,dataCode,iSubcase,dt=None):
        scalarObject.__init__(self,dataCode,iSubcase)
        self.energy  = {}
        self.percent = {}
        self.density = {}

        if dt is not None:
            self.dt = dt
            self.addNewTransient()
            self.isTransient = True
            self.add       = self.addTransient
        ###

    def updateDt(self,dataCode,dt):
        """
        this method is called if the object
        already exits and a new time step is found
        """
        self.dataCode = dataCode
        self.applyDataCode()
        #assert dt>=0.
        #print "updating dt...dt=%s" %(dt)
        if dt is not None:
            self.dt = dt
            self.addNewTransient()
        ###

    def addNewTransient(self):
        """
        initializes the transient variables
        @note make sure you set self.dt first
        """
        self.energy[self.dt] = {}
        self.percent[self.dt] = {}
        self.density[self.dt] = {}
        
    def add(self,grid,energy,percent,density):
        assert grid not in self.energy
        self.energy[grid]  = energy
        self.percent[grid] = percent
        self.density[grid] = density

    def addTransient(self,grid,energy,percent,density):
        dt = self.dt
        assert grid not in self.energy[dt]
        self.energy[dt][grid]  = energy
        self.percent[dt][grid] = percent
        self.density[dt][grid] = density
    
    def __repr__(self):
        msg  = '---Strain Energy Object---\n'
        msg += "%-14s "*4 %('EID','Energy','PercentTotal','density')+'\n'
        for eid,energy in sorted(self.energy.items()):
            percent = self.percent[eid]
            density = self.density[eid]
            msg += "%-14g"*4 %(eid,energy,percent,density)+'\n'
            
        return msg