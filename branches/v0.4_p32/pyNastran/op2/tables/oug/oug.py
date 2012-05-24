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
     )

from oug_accelerations import (             # tableCode=11,formatCode=1 sortCode=0
     accelerationObject,
     )

from oug_temperatures import (              # tableCode=1, formatCode=1 sortCode=0
     temperatureObject,
     )

from oug_eigenvectors import (
     eigenVectorObject,                     # analysisCode=2, sortCode=0 formatCode   tableCode=7
     complexEigenVectorObject,              # analysisCode=5, sortCode=1 formatCode=1 tableCode=7
     realEigenVectorObject,                 # analysisCode=9, sortCode=1 formatCode=1 tableCode=7
     )

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
        
        #self.printBlock(data) # on
        ## assuming tCode=1
        if self.analysisCode==1:   # statics / displacement / heat flux
            self.addDataParameter(data,'lsdvmn',  'i',5,False)   ## load set number
            self.applyDataCodeValue('dataNames',['lsdvmn'])
        elif self.analysisCode==2: # real eigenvalues
            self.addDataParameter(data,'mode',     'i',5)         ## mode number
            self.addDataParameter(data,'eigr',     'f',6,False)   ## real eigenvalue
            self.addDataParameter(data,'modeCycle','i',7,False)   ## mode or cycle @todo confused on the type - F1???
            self.applyDataCodeValue('dataNames',['mode','eigr','modeCycle'])

            assert self.dataCode['name']=='mode'
            #print "mode(5)=%s eigr(6)=%s modeCycle(7)=%s" %(self.mode,self.eigr,self.modeCycle)
        #elif self.analysisCode==3: # differential stiffness
        #    self.lsdvmn = self.getValues(data,'i',5) ## load set number
        #elif self.analysisCode==4: # differential stiffness
        #    self.lsdvmn = self.getValues(data,'i',5) ## load set number
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
        if   tfsCode==[1,1,0]:
            self.readOUG_Data_table1_format1_sort0()
        elif tfsCode==[1,1,1]:
            self.readOUG_Data_table1_format1_sort1()
        elif tfsCode==[1,2,0]:
            self.readOUG_Data_table1_format2_sort0()
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
        elif tfsCode==[7,1,0]:
            self.readOUG_Data_table7_format1_sort0()
        #elif tfsCode==[7,1,1]:
        #    self.readOUG_Data_table7_format1_sort1()
        elif tfsCode==[7,2,0]:
            self.readOUG_Data_table7_format2_sort0()
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
        elif tfsCode==[10,1,0]:
            self.readOUG_Data_table10_format1_sort0()
        #elif tfsCode==[10,1,1]:
        #    self.readOUG_Data_table10_format1_sort1()
        elif tfsCode==[10,2,0]:
            self.readOUG_Data_table10_format2_sort0()
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
        elif tfsCode==[11,1,0]:
            self.readOUG_Data_table11_format1_sort0()
        #elif tfsCode==[11,1,1]:
        #    self.readOUG_Data_table11_format1_sort1()
        elif tfsCode==[11,2,0]:
            self.readOUG_Data_table11_format2_sort0()
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

    def readOUG_Data_table7_format1_sort0(self):  # modes
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
            elif self.analysisCode==9: # complex eigenvector
                #print "isComplexEigenvector9"
                self.createTransientObject(self.eigenvectors,complexEigenVectorObject)
            else:
                #print self.codeInformation()
                #raise NotImplementedError('unsupported %s-OUG static solution...atfsCode=%s' %(self.tableName,self.atfsCode))
                pass
            ###
        elif self.thermal==1:
            #print self.codeInformation()
            #raise NotImplementedError('unsupported %s-OUG thermal solution...atfsCode=%s' %(self.tableName,self.atfsCode))
            pass
        else:
            #print self.codeInformation()
            #raise NotImplementedError('invalid %s-OUG thermal flag...not 0 or 1...flag=%s' %(self.tableName,self.thermal))
            pass
        ###
        if self.obj:
            self.readScalars8(debug=False)
        else:
            self.skipOES_Element()
        #if self.analysisCode not in [2,8]:
        #    raise NotImplementedError('check_format1...')
        ###

    def readOUG_Data_table10_format1_sort0(self): # velocity
        #assert self.formatCode==1 # Real
        #assert self.sortCode==0   # Real
        if self.thermal==0 or self.thermal>1:  ## @warning dont leave the thermal>0!!!!
            if self.analysisCode==1: # velocity
                #print "isVelocity"
                self.createTransientObject(self.velocities,velocityObject)
            elif self.analysisCode==6: # transient velocity
                #print "isTransientVelocity"
                self.createTransientObject(self.velocities,velocityObject)
            else:
                #print self.codeInformation()
                #raise NotImplementedError('unsupported OUG static solution...atfsCode=%s' %(self.atfsCode))
                pass
            ###
        elif self.thermal==1:
            #print self.codeInformation()
            #raise NotImplementedError('unsupported OUG thermal solution...atfsCode=%s' %(self.atfsCode))
            pass
        else:
            #print self.codeInformation()
            #raise NotImplementedError('invalid OUG thermal flag...not 0 or 1...flag=%s' %(self.thermal))
            pass
        ###
        self.readMappedScalarsOut(debug=False) # handles dtMap

        #if self.obj:
        #    self.readScalarsOut(debug=False)
        #    #self.readScalars8(debug=False)
        #else:
        #    self.skipOES_Element()
        ###

    def readOUG_Data_table10_format2_sort0(self): # velocity
        #assert self.formatCode==1 # Real
        #assert self.sortCode==0   # Real
        if self.thermal==0 or self.thermal>1:  ## @warning dont leave the thermal>0!!!!
            if 0:
                pass
            #if self.analysisCode==1: # velocity
                #print "isVelocity"
                #self.createTransientObject(self.velocities,velocityObject)
            elif self.analysisCode==6: # transient velocity
                self.createTransientObject(self.velocities,velocityObject)
            else:
                #print self.codeInformation()
                #raise NotImplementedError('unsupported OUG static solution...atfsCode=%s' %(self.atfsCode))
                pass
            ###
        elif self.thermal==1:
            #print self.codeInformation()
            #raise NotImplementedError('unsupported OUG thermal solution...atfsCode=%s' %(self.atfsCode))
            pass
        else:
            #print self.codeInformation()
            #raise NotImplementedError('invalid OUG thermal flag...not 0 or 1...flag=%s' %(self.thermal))
            pass
        ###
        self.readMappedScalarsOut(debug=False) # handles dtMap

        #if self.obj:
        #    self.readScalarsOut(debug=False)
        #    #self.readScalars8(debug=False)
        #else:
        #    self.skipOES_Element()
        ###

    def readOUG_Data_table11_format1_sort0(self): # acceleration
        #assert self.formatCode==1 # Real
        #assert self.sortCode==0   # Real
        if self.thermal==0 or self.thermal>1:  ## @warning dont leave the thermal>0!!!!
            if self.analysisCode==1: # acceleration
                #print "isAcceleration"
                self.createTransientObject(self.accelerations,accelerationObject)
            elif self.analysisCode==6: # transient acceleration
                #print "isTransientAcceleration"
                self.createTransientObject(self.accelerations,accelerationObject)
            else:
                #print self.codeInformation()
                #raise NotImplementedError('unsupported %s-OUG static solution...atfsCode=%s' %(self.tableName,self.atfsCode))
                pass
            ###
        elif self.thermal==1:
            #print self.codeInformation()
            #raise NotImplementedError('unsupported %s-OUG thermal solution...atfsCode=%s' %(self.tableName,self.atfsCode))
            pass
        else:
            #print self.codeInformation()
            #raise NotImplementedError('invalid %s-OUG thermal flag...not 0 or 1...flag=%s' %(self.tableName,self.thermal))
            pass
        ###
        self.readMappedScalarsOut(debug=False) # handles dtMap

    def readOUG_Data_table11_format2_sort0(self): # acceleration
        #assert self.formatCode==1 # Real
        #assert self.sortCode==0   # Real
        if self.thermal==0 or self.thermal>1:  ## @warning dont leave the thermal>0!!!!
            if 0:
                pass
            #if self.analysisCode==1: # acceleration
                #print "isAcceleration"
                #self.createTransientObject(self.accelerations,accelerationObject)
            elif self.analysisCode==6: # transient acceleration
                #print "isTransientAcceleration"
                self.createTransientObject(self.accelerations,accelerationObject)
            else:
                #print self.codeInformation()
                #raise NotImplementedError('unsupported %s-OUG static solution...atfsCode=%s' %(self.tableName,self.atfsCode))
                pass
            ###
        elif self.thermal==1:
            #print self.codeInformation()
            #raise NotImplementedError('unsupported %s-OUG thermal solution...atfsCode=%s' %(self.tableName,self.atfsCode))
            pass
        else:
            #print self.codeInformation()
            #raise NotImplementedError('invalid %s-OUG thermal flag...not 0 or 1...flag=%s' %(self.tableName,self.thermal))
            pass
        ###
        self.readMappedScalarsOut(debug=False) # handles dtMap

    def readScalarsX1(n,sFormat,debug):
        #self.readScalarsX1(8,sFormat,debug=False)
        if n==8:
            if   sFormat=='iiffffff':
                self.readScalars8(debug)
            elif sFormat=='fiffffff':
                self.readScalarsF8(debug)
            else:
                #print self.codeInformation()
                #raise RuntimeError('not supported format...')
                pass
            ###
        ###
        else:
            #print self.codeInformation()
            #raise RuntimeError('not supported format...')
            pass
        ###

    def readOUG_Data_table1_format1_sort1(self): # displacement
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
                #print self.codeInformation()
                #raise NotImplementedError('unsupported %s-OUG static table1_format1_sort1 solution...atfsCode=%s' %(self.tableName,self.atfsCode))
                pass
            ###
        else:
            #print self.codeInformation()
            #raise NotImplementedError('invalid thermal flag...not 0 or 1...flag=%s' %(self.thermal))
            pass
        ###
        #print "objName = ",self.obj.name()
        if self.obj:
            self.readScalars14(debug=False)
        else:
            self.skipOES_Element()
        #print "---OBJ---"
        #print self.obj
        #raise Exception('format1_sort1')

    def readOUG_Data_table7_format1_sort1(self): # modes
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
                #print self.codeInformation()
                #raise NotImplementedError('unsupported %s-OUG static table7_format1_sort1 solution...atfsCode=%s' %(self.tableName,self.atfsCode))
                pass
            ###
        else:
            #print self.codeInformation()
            #raise NotImplementedError('invalid thermal flag...not 0 or 1...flag=%s' %(self.thermal))
            pass
        ###
        #print "objName = ",self.obj.name()
        if self.obj:
            self.readScalars14(debug=False)
        else:
            self.skipOES_Element()
        #print self.obj

    def readOUG_Data_table10_format1_sort1(self): # velocity
        #assert self.formatCode==1 # Real
        #assert self.sortCode==1   # Real/Imaginary
        #print self.codeInformation()
        if self.thermal==0:
            if self.analysisCode==5: # frequency velocity
                #print "isFrequencyVelocity"
                self.createTransientObject(self.velocities,complexEigenVectorObject)
            #elif self.analysisCode==7: # pre-buckling displacement
                #print "isPreBucklingDisplacement"
                #self.createTransientObject(self.displacements,displacementObject)
            elif self.analysisCode==9: # nonlinear static eigenvector
                #print "isComplexVelocities"
                self.createTransientObject(self.velocites,complexDisplacementObject)
            #elif self.analysisCode==11: # Geometric nonlinear statics
                #print "isNonlinearStaticDisplacement"
                #self.createTransientObject(self.displacements,displacementObject)
            else:
                #print self.codeInformation()
                #raise NotImplementedError('unsupported %s-OUG static table7_format1_sort1 solution...atfsCode=%s' %(self.tableName,self.atfsCode))
                pass
            ###
        else:
            #print self.codeInformation()
            #raise NotImplementedError('invalid thermal flag...not 0 or 1...flag=%s' %(self.thermal))
            pass
        ###
        #print "objName = ",self.obj.name()
        #print self.obj
        if self.obj:
            self.readScalars14()
        else:
            self.skipOES_Element()

    def readOUG_Data_table1_format2_sort0(self): # displacement
        if self.thermal==0:
            if self.analysisCode==6: # transient displacement
                #print "isTransientDisplacement"
                self.createTransientObject(self.displacements,displacementObject)
            else:
                #print self.codeInformation()
                #raise NotImplementedError('unsupported %s-OUG static solution...atfsCode=%s' %(self.tableName,self.atfsCode))
                pass
            ###
        elif self.thermal==1:
            #print self.codeInformation()
            #raise NotImplementedError('unsupported %s-OUG thermal solution...atfsCode=%s' %(self.tableName,self.atfsCode))
            pass
        else:
            #print self.codeInformation()
            #raise NotImplementedError('invalid thermal flag...not 0 or 1...flag=%s' %(self.thermal))
            pass
        ###
        if self.obj:
            #self.readScalars8()
            self.readScalarsOut(debug=False)
        else:
            self.skipOES_Element()

    def readOUG_Data_table1_format2_sort1(self): # displacement
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
                #print self.codeInformation()
                #raise NotImplementedError('unsupported %s-OUG static solution...atfsCode=%s' %(self.tableName,self.atfsCode))
                pass
            ###
        elif self.thermal==1:
            #print self.codeInformation()
            #raise NotImplementedError('unsupported %s-OUG thermal solution...atfsCode=%s' %(self.tableName,self.atfsCode))
            pass
        else:
            #print self.codeInformation()
            #raise NotImplementedError('invalid thermal flag...not 0 or 1...flag=%s' %(self.thermal))
            pass
        ###
        #print self.obj

    def readOUG_Data_table1_format3_sort0(self): # displacement
        #assert self.formatCode==3 # Magnitude/Phase
        #assert self.sortCode==0   # Real
        #print self.codeInformation()
        if self.thermal==0:
            #print self.codeInformation()
            #raise NotImplementedError('unsupported OUG static solution...atfsCode=%s' %(self.atfsCode))
            
            if 1:
                pass # break
            elif self.analysisCode==1: # displacement
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
                #print self.codeInformation()
                #raise NotImplementedError('unsupported %s-OUG static solution...atfsCode=%s' %(self.tableName,self.atfsCode))
                pass
            ###
        elif self.thermal==1:
            #print self.codeInformation()
            #raise NotImplementedError('unsupported %s-OUG thermal solution...atfsCode=%s' %(self.tableName,self.atfsCode))
            pass
        else:
            #print self.codeInformation()
            #raise NotImplementedError('invalid thermal flag...not 0 or 1...flag=%s' %(self.thermal))
            pass
        ###
        #self.readScalars8()
        self.readMappedScalarsOut(debug=False) # handles dtMap

    def readOUG_Data_table1_format3_sort1(self): # displacemnt
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
                #print self.codeInformation()
                #raise NotImplementedError('unsupported %s-OUG static solution...atfsCode=%s' %(self.tableName,self.atfsCode))
                pass
            ###
        elif self.thermal==1:
            #print self.codeInformation()
            #raise NotImplementedError('unsupported %s-OUG thermal solution...atfsCode=%s' %(self.tableName,self.atfsCode))
            pass
        else:
            #print self.codeInformation()
            #raise NotImplementedError('invalid thermal flag...not 0 or 1...flag=%s' %(self.thermal))
            pass
        ###
        if self.obj:
            self.readScalars14()
        else:
            self.skipOES_Element()
        ###
