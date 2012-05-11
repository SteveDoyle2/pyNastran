import sys
import copy
from struct import unpack

# pyNastran
from oqg_spcForces import spcForcesObject,complexSpcForcesObject
from oqg_mpcForces import mpcForcesObject,complexMpcForcesObject

from pyNastran.op2.tables.oug.oug_Objects import (
    #temperatureObject,displacementObject,  # analysisCode=1, sortCode=0
     fluxObject,                            # analysisCode=1, sortCode=3
     #nonlinearTemperatureObject,            # analysisCode=10,sortCode=0
     )

from pyNastran.op2.tables.oug.oug_eigenvectors import (
     eigenVectorObject,                     # analysisCode=2, sortCode=0 formatCode   tableCode=7
     complexEigenVectorObject,              # analysis    =5, sortCode=1 formatCode=1 tableCode=7
     realEigenVectorObject,                 # analysisCode=9, sortCode=1 formatCode=1 tableCode=7
     )

class OQG(object):
    """Table of spc/mpc forces/momenets"""

    def readTable_OQG(self):
        table3 = self.readTable_OQG_3
        table4Data = self.readOQG_Data
        self.readResultsTable(table3,table4Data)
        self.deleteAttributes_OQG()

    def deleteAttributes_OQG(self):
        #print self.obj
        params = ['lsdvm','mode','eigr','modeCycle','freq','dt','lftsfq','thermal','rCode','fCode','numWide','acousticFlag','thermal']
        self.deleteAttributes(params)
    
    def readTable_OQG_3(self,iTable): # iTable=-3
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
        self.addDataParameter(data,'thermal',     'i',23,False)  ## thermal flag; 1 for heat ransfer, 0 otherwise
        
        if not self.isSort1():
            raise NotImplementedError('sort2...')
        #assert self.isThermal()==False,self.thermal

        #self.printBlock(data) # on
        ## assuming tCode=1
        if self.analysisCode==1:   # statics / displacement / heat flux
            self.addDataParameter(data,'lsdvmn',  'i',5,False)   ## load set number
            self.applyDataCodeValue('dataNames',['lsdvmn'])
            self.setNullNonlinearFactor()
        elif self.analysisCode==2: # real eigenvalues
            self.addDataParameter(data,'mode',     'i',5)         ## mode number
            self.addDataParameter(data,'eigr',     'f',6,False)   ## real eigenvalue
            self.addDataParameter(data,'modeCycle','f',7,False)   ## mode or cycle @todo confused on the type - F1???
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
            self.applyDataCodeValue('dataNames',['dt'])
        elif self.analysisCode==7: # pre-buckling
            self.addDataParameter(data,'lsdvmn',  'i',5)   ## load set number
            self.applyDataCodeValue('dataNames',['lsdvmn'])
        elif self.analysisCode==8: # post-buckling
            self.addDataParameter(data,'lsdvmn',  'i',5)         ## load set number
            self.addDataParameter(data,'eigr',    'f',6,False)   ## real eigenvalue
            self.applyDataCodeValue('dataNames',['lsdvmn','eigr'])
        elif self.analysisCode==9: # complex eigenvalues
            self.addDataParameter(data,'mode','i',5)         ## mode number
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


    def readOQG_Data(self):
        tfsCode = [self.tableCode,self.formatCode,self.sortCode]
        #self.skipOES_Element() # skipping entire table
        #return

        #print "self.analysisCode=%s tableCode(1)=%s thermal(23)=%g" %(self.analysisCode,self.tableCode,self.thermal)
        assert self.thermal in [0,1]

        if self.tableCode==3:  # SPC Force vector
            self.readOQG_Data_table3()
        elif self.tableCode==39:  # MPC Force vector
            self.readOQG_Data_table3()
        else:
            self.skipOES_Element()
            #print self.codeInformation()
            #raise NotImplementedError('bad analysis/table/format/sortCode=%s' %(self.atfsCode))
        ###
        #print self.obj

    def readOQG_Data_table3(self): # SPC Forces
        isSort1 = self.isSort1()
        #print self.codeInformation()
        magPhase = self.isMagnitudePhase()
        if magPhase or self.numWide==14:  # real/imaginary or mag/phase
            if self.thermal==0:
                self.createTransientObject(self.spcForces,complexSpcForcesObject) # complex
            else:
                raise NotImplementedError(self.codeInformation())
            self.OUG_ComplexTable()
        elif self.numWide==8:  # real/random
            #if self.thermal==0:
            self.createTransientObject(self.spcForces,spcForcesObject) # real
            #else:
            #    raise NotImplementedError(self.codeInformation())
            self.OUG_RealTable()
        else:
            raise NotImplementedError('only numWide=8 or 14 is allowed  numWide=%s' %(self.numWide))
        ###
        #if self.thermal!=0:
            #print self.obj
            #raise RuntimeError('check the printout for thermal...')

    def readOQG_Data_table39(self): # MPC Forces
        isSort1 = self.isSort1()
        if self.numWide==8:  # real/random
            if self.thermal==0:
                self.createTransientObject(self.mpcForces,mpcForcesObject) # real
            else:
                raise NotImplementedError(self.codeInformation())
            self.OUG_RealTable()
        elif self.numWide==14:  # real/imaginary or mag/phase
            if self.thermal==0:
                self.createTransientObject(self.mpcForces,complexMpcForcesObject) # complex
            else:
                raise NotImplementedError(self.codeInformation())
            self.OUG_ComplexTable()
        else:
            raise NotImplementedError('only numWide=8 or 14 is allowed  numWide=%s' %(self.numWide))
        ###

