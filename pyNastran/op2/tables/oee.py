import sys
import copy
from struct import unpack

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
        #self.staticEnergy = {} # aCode=1 tCode=18 fCode=1 sortCode=0

        table3 = self.readTable_OEE1_3
        table4Data = self.readTable_OEE1_4_Data
        self.readResultsTable(table3,table4Data)
        self.deleteAttributes_OEE()
        #sys.exit('end of oee')

    def deleteAttributes_OEE(self):
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
        
        #self.printBlock(data) # on
        if self.approachCode==1:   # statics / displacement / heat flux
            pass
        elif self.approachCode==2: # real eigenvalues
            self.mode      = self.getValues(data,'i',5) ## mode number
            print "mode(5)=%s" %(self.mode)
        elif self.approachCode==3: # differential stiffness
            pass
        elif self.approachCode==4: # differential stiffness
            pass
        elif self.approachCode==5:   # frequency
            self.freq = self.getValues(data,'f',5) ## frequency

        elif self.approachCode==6: # transient
            self.time = self.getValues(data,'f',5) ## time step
            print "time(5)=%s" %(self.time)
        elif self.approachCode==7: # pre-buckling
            pass
        elif self.approachCode==8: # post-buckling
            self.mode = self.getValues(data,'i',5) ## mode number
            print "mode(5)=%s" %(self.mode)
        elif self.approachCode==9: # complex eigenvalues
            self.mode   = self.getValues(data,'i',5) ## mode number
            print "mode(5)=%s" %(self.mode)
        elif self.approachCode==10: # nonlinear statics
            self.loadFactor = self.getValues(data,'f',5) ## load factor
            print "loadFactor(5) = %s" %(self.loadFactor)
        elif self.approachCode==11: # old geometric nonlinear statics
            pass
        elif self.approachCode==12: # contran ? (may appear as aCode=6)  --> straight from DMAP...grrr...
            self.time = self.getValues(data,'f',5) ## time step
            print "time(5)=%s" %(self.time)
        else:
            raise RuntimeError('invalid approach code...approachCode=%s' %(self.approachCode))
        ###
        
        print "*iSubcase=%s elementName=|%s|"%(self.iSubcase,self.elementName)
        print "approachCode=%s tableCode=%s" %(self.approachCode,self.tableCode)
        print self.codeInformation()

        #self.printBlock(data)
        self.readTitle()

    def readTable_OEE1_4_Data(self,iTable): # iTable=-4
        isTable4Done = False
        isBlockDone  = False

        bufferWords = self.getMarker('OEE1')
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
        print "self.approachCode=%s tableCode(1)=%s" %(self.approachCode,self.tableCode)

        self.readOEE1_Data()
        #print "-------finished OEE1----------"
        return (isTable4Done,isBlockDone)

    def readOEE1_Data(self):
        tfsCode = [self.tableCode,self.formatCode,self.sortCode]
        
        if tfsCode==[18,1,0]:
            self.readStrainEnergy_format1_sort0()
        else:
            self.skipOES_Element(None)
            raise Exception('unsupported OEE1 static solution...')
        ###
        print str(self.obj)
        
        #if fsCode==[1,0]:
        #    self.readOEE1_Data_format1_sort0()
        #elif fsCode==[1,1]:
        #    self.readOEE1_Data_format1_sort1()
        #elif fsCode==[2,1]:
        #    self.readOEE1_Data_format2_sort1()
        #else:
        #    raise Exception('bad formatCode/sortCode')
        ###

    def readStrainEnergy_format1_sort0(self):
        assert self.tableCode==18 # Strain Energy
        assert self.formatCode==1 # Real
        assert self.sortCode==0   # Real

        if self.approachCode==1: # displacement
            print "isStrainEnergy"
            self.obj = StrainEnergyObject(self.iSubcase)
            self.strainEnergy[self.iSubcase] = self.obj
        elif self.approachCode==2: # buckling modes
            print "isBucklingStrainEnergy"
            self.obj = StrainEnergyObject(self.iSubcase)
            self.modesStrainEnergy[self.iSubcase] = self.obj
        elif self.approachCode==9: # nonlinear static eigenvector
            print "isComplexStrainEnergy"
            self.obj = StrainEnergyObject(self.iSubcase)
            self.complexStrainEnergy[self.iSubcase] = self.obj
        else:
            raise Exception('bad tableCode/formatCode/sortCode=%s on OEE table' %(tfsCode))
        ###
        self.readScalars4(self.obj)

class StrainEnergyObject(object):
    def __init__(self,iSubcase):
        self.energy = {}
        self.percent = {}
        self.density = {}
        
    def add(self,grid,energy,percent,density):
        assert grid not in self.energy
        self.energy[grid]  = energy
        self.percent[grid] = percent
        self.density[grid] = density
    
    def __repr__(self):
        msg  = '---Strain Energy Object---\n'
        msg += "%-14s "*4 %('EID','Energy','PercentTotal','density')+'\n'
        for eid,energy in sorted(self.energy.items()):
            percent = self.percent[eid]
            density = self.density[eid]
            msg += "%-14g"*4 %(eid,energy,percent,density)+'\n'
            
        return msg