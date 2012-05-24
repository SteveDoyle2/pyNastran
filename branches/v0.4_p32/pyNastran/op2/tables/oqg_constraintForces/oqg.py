## GNU Lesser General Public License
## 
## Program pyNastran - a python interface to NASTRAN files
## Copyright (C) 2011-2012  Steven Doyle, Al Danial
## 
## Authors and copyright holders of pyNastran
## Steven Doyle <mesheb82@gmail.com>
## Al Danial    <al.danial@gmail.com>
## 
## This file is part of pyNastran.
## 
## pyNastran is free software: you can redistribute it and/or modify
## it under the terms of the GNU Lesser General Public License as published by
## the Free Software Foundation, either version 3 of the License, or
## (at your option) any later version.
## 
## pyNastran is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
## 
## You should have received a copy of the GNU Lesser General Public License
## along with pyNastran.  If not, see <http://www.gnu.org/licenses/>.
## 
import sys
import copy
from struct import unpack

# pyNastran
from .oqg_spcForces import spcForcesObject,complexSpcForcesObject
from .oqg_mpcForces import mpcForcesObject,complexMpcForcesObject

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
        #sys.exit('stopping in oqg1.py')
        
    
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
        
        #self.printBlock(data) # on
        ## assuming tCode=1
        if self.analysisCode==1:   # statics / displacement / heat flux
            self.addDataParameter(data,'lsdvmn',  'i',5,False)   ## load set number
        elif self.analysisCode==2: # real eigenvalues
            self.addDataParameter(data,'mode',     'i',5)         ## mode number
            self.addDataParameter(data,'eigr',     'f',6,False)   ## real eigenvalue
            self.addDataParameter(data,'modeCycle','f',7,False)   ## mode or cycle @todo confused on the type - F1???
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
        elif self.analysisCode==8: # post-buckling
            self.addDataParameter(data,'lsdvmn',  'i',5)         ## load set number
            self.addDataParameter(data,'eigr',    'f',6,False)   ## real eigenvalue
            self.applyDataCodeValue('dataNames',['lsdvmn','eigr'])
        elif self.analysisCode==9: # complex eigenvalues
            self.addDataParameter(data,'mode','i',5)         ## mode number
            self.addDataParameter(data,'eigr','f',6,False)   ## real eigenvalue
            self.addDataParameter(data,'eigi','f',7,False)   ## imaginary eigenvalue
        elif self.analysisCode==10: # nonlinear statics
            self.addDataParameter(data,'lftsfq','f',5)   ## load step
        elif self.analysisCode==11: # old geometric nonlinear statics
            self.addDataParameter(data,'lsdvmn',  'i',5)   ## load set number
        elif self.analysisCode==12: # contran ? (may appear as aCode=6)  --> straight from DMAP...grrr...
            self.addDataParameter(data,'lsdvmn',  'i',5)   ## load set number
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

        if tfsCode==[3,1,0]:  # SPC Force vector
            self.readOQG_Data_table3_format1_sort0()
        elif tfsCode==[3,1,1]:
            self.readOQG_Data_table3_format1_sort1()
        elif tfsCode==[3,2,0]:
            self.readOQG_Data_table3_format2_sort0()
        #elif tfsCode==[3,2,1]:
            #self.skipOES_Element()
            #self.readOQG_Data_table3_format2_sort1()
        #elif tfsCode==[3,3,0]:
            #self.readOQG_Data_table3_format3_sort0()
            #self.skipOES_Element()
        #elif tfsCode==[3,3,1]:
            #self.skipOES_Element()
            #self.readOQG_Data_table3_format3_sort1()
        
        elif tfsCode==[39,1,0]: # MPC Forces
            self.readOQG_Data_table39_format1_sort0()
        #elif tfsCode==[39,1,1]:
            #self.skipOES_Element()
            #self.readOQG_Data_table39_format1_sort1()
        else:
            self.skipOES_Element()
            #print self.codeInformation()
            #raise NotImplementedError('bad analysis/table/format/sortCode=%s' %(self.atfsCode))
        ###
        #print self.obj

    def readOQG_Data_table39_format1_sort0(self):
        if self.thermal==0:
            if self.analysisCode==1: # static MPC forces
                self.createTransientObject(self.mpcForces,mpcForcesObject)
            elif self.analysisCode==2: # eigenvector
                self.createTransientObject(self.mpcForces,mpcForcesObject)
            elif self.analysisCode==6: # transient MPC forces
                self.createTransientObject(self.mpcForces,mpcForcesObject)
            elif self.analysisCode==10: # nonlinear MPC forces
                self.createTransientObject(self.mpcForces,mpcForcesObject)
            else:
                pass
                #self.skipOES_Element()
                #print self.codeInformation()
                #raise NotImplementedError('unsupported OQG static solution...atfsCode=%s' %(self.atfsCode))
            ###
        elif self.thermal==1:
            #raise NotImplementedError('thermal not supported for MPC forces...atfsCode=%s' %(self.atfsCode))
            if 0:
                pass
            else:
                #self.skipOES_Element()
                self.log.debug('unsupported OQG thermal solution...atfsCode=%s' %(self.atfsCode))
                #print self.codeInformation()
                #raise NotImplementedError('unsupported OQG thermal solution...atfsCode=%s' %(self.atfsCode))
            ###
        ###
        else:
            #raise NotImplementedError('invalid OQG thermal flag...not 0 or 1...flag=%s' %(self.thermal))
            pass
        ###
        #print "objName = ",self.obj.name()
        self.readMappedScalarsOut(debug=False) # handles dtMap

    def readOQG_Data_table3_format2_sort0(self):
        if self.thermal==0:
            #if self.analysisCode==1: # static SPC forces
                #self.createTransientObject(self.spcForces,spcForcesObject)
            #elif self.analysisCode==2: # nonlinear static eigenvector
                #print "isEigenvector"
                #self.createTransientObject(self.spcForces,spcForcesObject)
            #elif self.analysisCode==5: # frequency
                #print "isFrequencyForces"
                #self.createTransientObject(self.modalSPCForces,eigenVectorObject)
            if self.analysisCode==6: # transient SPC forces
                self.createTransientObject(self.spcForces,spcForcesObject)
            #elif self.analysisCode==8: # post-buckling SPC forces
                #self.createTransientObject(self.spcForces,spcForcesObject)
            #elif self.analysisCode==10: # nonlinear static SPC forces
                #self.createTransientObject(self.spcForces,spcForcesObject)
            #elif self.analysisCode==11: # Geometric nonlinear statics
                #print "isNonlinearStaticDisplacement"
                #self.createTransientObject(self.spcForces,spcForcesObject)
            else:
                #self.skipOES_Element()
                pass
                #print self.codeInformation()
                #raise NotImplementedError('unsupported OQG static solution...atfsCode=%s' %(self.atfsCode))
            ###
        elif self.thermal==1:
            #if self.analysisCode==1: # temperature
                #print "isTemperature"
                #self.obj = temperatureObject(self.dataCode,self.iSubcase)
                #self.temperatures[self.iSubcase] = self.obj
            #elif self.analysisCode==6: # transient forces
                #print "isTransientForces"
                #self.createTransientObject(self.nonlinearTemperatures,nonlinearTemperatureObject)
            #elif self.analysisCode==10: # nonlinear static displacement
                #print "isNonlinearStaticTemperatures"
                #self.createTransientObject(self.nonlinearTemperatures,nonlinearTemperatureObject)
            if 0:
                pass
            else:
                #self.skipOES_Element()
                #pass
                self.log.debug('unsupported OQG thermal solution...atfsCode=%s' %(self.atfsCode))
                #raise NotImplementedError('unsupported OQG thermal solution...atfsCode=%s' %(self.atfsCode))
            ###
        ###
        else:
            #raise NotImplementedError('invalid OQG thermal flag...not 0 or 1...flag=%s' %(self.thermal))
            pass
        ###
        self.readMappedScalarsOut(debug=False) # handles dtMap

    def readOQG_Data_table3_format3_sort0(self):
        if self.thermal==0:
            #if self.analysisCode==1: # static SPC forces
                #self.createTransientObject(self.spcForces,spcForcesObject)
            #elif self.analysisCode==2: # nonlinear static eigenvector
                #print "isEigenvector"
                #self.createTransientObject(self.spcForces,spcForcesObject)
            #elif self.analysisCode==5: # frequency
                #print "isFrequencyForces"
                #self.createTransientObject(self.modalSPCForces,eigenVectorObject)
            if self.analysisCode==6: # transient SPC forces
                self.createTransientObject(self.spcForces,spcForcesObject)
            #elif self.analysisCode==8: # post-buckling SPC forces
                #self.createTransientObject(self.spcForces,spcForcesObject)
            #elif self.analysisCode==10: # nonlinear static SPC forces
                #self.createTransientObject(self.spcForces,spcForcesObject)
            #elif self.analysisCode==11: # Geometric nonlinear statics
                #print "isNonlinearStaticDisplacement"
                #self.createTransientObject(self.spcForces,spcForcesObject)
            else:
                #self.skipOES_Element()
                pass
                #print self.codeInformation()
                #raise NotImplementedError('unsupported OQG static solution...atfsCode=%s' %(self.atfsCode))
            ###
        elif self.thermal==1:
            #if self.analysisCode==1: # temperature
                #print "isTemperature"
                #self.obj = temperatureObject(self.dataCode,self.iSubcase)
                #self.temperatures[self.iSubcase] = self.obj
            #elif self.analysisCode==6: # transient forces
                #print "isTransientForces"
                #self.createTransientObject(self.nonlinearTemperatures,nonlinearTemperatureObject)
            #elif self.analysisCode==10: # nonlinear static displacement
                #print "isNonlinearStaticTemperatures"
                #self.createTransientObject(self.nonlinearTemperatures,nonlinearTemperatureObject)
            if 0:
                pass
            else:
                #self.skipOES_Element()
                #pass
                self.log.debug('unsupported OQG thermal solution...atfsCode=%s' %(self.atfsCode))
                #raise NotImplementedError('unsupported OQG thermal solution...atfsCode=%s' %(self.atfsCode))
            ###
        ###
        else:
            #raise NotImplementedError('invalid OQG thermal flag...not 0 or 1...flag=%s' %(self.thermal))
            pass
        ###
        self.readMappedScalarsOut(debug=False) # handles dtMap

    def readOQG_Data_table3_format1_sort0(self):
        if self.thermal==0:
            if self.analysisCode==1: # static SPC forces
                self.createTransientObject(self.spcForces,spcForcesObject)
            elif self.analysisCode==2: # nonlinear static eigenvector
                #print "isEigenvector"
                self.createTransientObject(self.spcForces,spcForcesObject)
            #elif self.analysisCode==5: # frequency
                #print "isFrequencyForces"
                #self.createTransientObject(self.modalSPCForces,eigenVectorObject)
            elif self.analysisCode==6: # transient SPC forces
                self.createTransientObject(self.spcForces,spcForcesObject)
            elif self.analysisCode==8: # post-buckling SPC forces
                self.createTransientObject(self.spcForces,spcForcesObject)
            elif self.analysisCode==10: # nonlinear static SPC forces
                self.createTransientObject(self.spcForces,spcForcesObject)
            elif self.analysisCode==11: # Geometric nonlinear statics
                #print "isNonlinearStaticDisplacement"
                self.createTransientObject(self.spcForces,spcForcesObject)
            else:
                #self.skipOES_Element()
                pass
                #print self.codeInformation()
                #raise NotImplementedError('unsupported OQG static solution...atfsCode=%s' %(self.atfsCode))
            ###
        elif self.thermal==1:
            #if self.analysisCode==1: # temperature
                #print "isTemperature"
                #self.obj = temperatureObject(self.dataCode,self.iSubcase)
                #self.temperatures[self.iSubcase] = self.obj
            #elif self.analysisCode==6: # transient forces
                #print "isTransientForces"
                #self.createTransientObject(self.nonlinearTemperatures,nonlinearTemperatureObject)
            #elif self.analysisCode==10: # nonlinear static displacement
                #print "isNonlinearStaticTemperatures"
                #self.createTransientObject(self.nonlinearTemperatures,nonlinearTemperatureObject)
            if 0:
                pass
            else:
                #self.skipOES_Element()
                pass
                self.log.debug('unsupported OQG thermal solution...atfsCode=%s' %(self.atfsCode))
                #raise NotImplementedError('unsupported OQG thermal solution...atfsCode=%s' %(self.atfsCode))
            ###
        ###
        else:
            #raise NotImplementedError('invalid OQG thermal flag...not 0 or 1...flag=%s' %(self.thermal))
            pass
        ###
        self.readMappedScalarsOut(debug=False) # handles dtMap


    def readOQG_Data_table3_format1_sort1(self):
        if self.thermal==0:
            if self.analysisCode==5: # frequency
                #print "isFrequencyForces"
                self.createTransientObject(self.spcForces,complexSpcForcesObject)
            elif self.analysisCode==9: # frequency
                #print "isComplexEigenvalueForces"
                self.createTransientObject(self.spcForces,complexSpcForcesObject)
            ###
            else:
                #raise NotImplementedError('unsupported OQG static solution...atfsCode=%s' %(self.atfsCode))
                pass
            ###
        elif thermal==1:
            #raise NotImplementedError('unsupported OQG thermal solution...atfsCode=%s' %(self.atfsCode))
            pass
        else:
            #raise NotImplementedError('invalid OQG thermal flag...not 0 or 1...flag=%s' %(self.thermal))
            pass
        ###
        #print "objName = ",self.obj.name()
        if self.obj:
            self.readScalars14(debug=False)
        else:
            self.skipOES_Element()

    def readOQG_Data_format2_sort1(self):
        if self.thermal==0:
            if self.analysisCode==5: # frequency
                #print "isFrequencyForces"
                self.createTransientObject(self.spcForces,complexSpcForcesObject)
            else:
                #raise NotImplementedError('unsupported OQG static solution...atfsCode=%s' %(self.atfsCode))
                pass
            ###
        elif thermal==1:
            #raise NotImplementedError('unsupported OQG thermal solution...atfsCode=%s' %(self.atfsCode))
            pass
        else:
            #raise NotImplementedError('invalid OQG thermal flag...not 0 or 1...flag=%s' %(self.thermal))
            pass
        ###
        #print "objName = ",self.obj.name()
        if self.obj:
            self.readScalars14(debug=True) # readImaginary
        else:
            self.skipOES_Element()

    def readOQG_Data_format3_sort1(self):
        if self.thermal==0:
            if self.analysisCode==5: # frequency
                #print "isFrequencyForces"
                self.createTransientObject(self.modalSPCForces,complexEigenVectorObject)
            elif self.analysisCode==6: # transient forces
                #print "isTransientForces"
                self.createTransientObject(self.spcForces,spcForcesObject)
            elif self.analysisCode==8: # post-buckling forces
                #print "isPostBucklingForces"
                self.createTransientObject(self.spcForces,spcForcesObject)
            elif self.analysisCode==11: # Geometric nonlinear statics
                #print "isFrequencyForces"
                self.createTransientObject(self.nonlinearForces,eigenVectorObject)
            else:
                #raise NotImplementedError('unsupported OQG static solution...atfsCode=%s' %(self.atfsCode))
                pass
            ###
        elif thermal==1:
            #raise NotImplementedError('unsupported OQG thermal solution...atfsCode=%s' %(self.atfsCode))
            pass
        ###
        else:
            #raise NotImplementedError('invalid OQG thermal flag...not 0 or 1...flag=%s' %(self.thermal))
            pass
        ###
        #print "objName = ",self.obj.name()
        if self.obj:
            self.readScalars14(debug=True) # readImaginary
        else:
            self.skipOES_Element()

