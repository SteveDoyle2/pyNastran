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
from pyNastran.op2.op2Errors import *
from pyNastran.op2.tables.oes_stressStrain.elementsStressStrain import ElementsStressStrain
from pyNastran.op2.tables.oes_stressStrain.oes_rods   import rodStressObject,rodStrainObject
from pyNastran.op2.tables.oes_stressStrain.oes_shear  import shearStressObject,shearStrainObject
from pyNastran.op2.tables.oes_stressStrain.oes_bars   import barStressObject,barStrainObject
from pyNastran.op2.tables.oes_stressStrain.oes_beams  import beamStressObject,beamStrainObject

from pyNastran.op2.tables.oes_stressStrain.oes_solids import solidStressObject,solidStrainObject
from pyNastran.op2.tables.oes_stressStrain.oes_plates import plateStressObject,plateStrainObject
from pyNastran.op2.tables.oes_stressStrain.oes_compositePlates import compositePlateStressObject,compositePlateStrainObject
from pyNastran.op2.tables.oes_stressStrain.oes_springs   import celasStressObject,celasStrainObject
from pyNastran.op2.tables.oes_stressStrain.oes_nonlinear import nonlinearRodObject,nonlinearQuadObject,hyperelasticQuadObject
from pyNastran.op2.tables.oes_stressStrain.oes_triax import ctriaxStressObject,ctriaxStrainObject

class OES(ElementsStressStrain):
    """Table of stresses/strains"""

    def readTable_OES(self):
        table3 = self.readTable_OES_3
        table4Data = self.readTable_OES_4_Data
        self.readResultsTable(table3,table4Data,flag=1) # flag=1 defines old style
        self.deleteAttributes_OES()

    def deleteAttributes_OES(self):
        params = ['sCode','elementType','obj','markerStart','loadSet','formatCode','sCode','thermal',
                  'lsdvmn','mode','eign','modeCycle','freq','mode','eigr','eigi','dt']
        self.deleteAttributes(params)

    def readTable_OES_3(self,iTable):
        #print "*iTable3 = ",iTable
        #if 0:
            #markers = self.readMarkers([0,2])
            #print "markers=%s" %(markers)
            #block = self.readBlock()
            #print "block = ",block
            #markers = self.readMarkers([-1,7])
            #print "markers=%s" %(markers)
            #print self.printSection(200)
        
        bufferWords = self.getBufferWords()
        
        data = self.getData(4)
        bufferSize, = unpack('i',data)
        if self.makeOp2Debug:
            self.op2Debug.write('bufferSize=|%s|\n' %(str(bufferSize)))

        data = self.getData(4*50)
        #self.printBlock(data)
        if self.makeOp2Debug:
            self.op2Debug.write('block3header\n')
        
        self.parseApproachCode(data) # 3
        self.addDataParameter(data,'elementType', 'i',3,False)   ## element type
        self.addDataParameter(data,'loadSet',     'i',8,False)   ## load set ID
        self.addDataParameter(data,'formatCode',  'i',9,False)   ## format code
        self.addDataParameter(data,'numWide',     'i',10,False)  ## number of words per entry in record; @note is this needed for this table ???
        self.addDataParameter(data,'sCode',       'i',11,False)  ## stress/strain codes
        self.addDataParameter(data,'thermal',     'i',23,False)  ## thermal flag; 1 for heat ransfer, 0 otherwise

        #print "loadset=%s formatCode=%s numWordsEntry=%s sCode=%s" %(self.loadSet,self.formatCode,self.numWide,self.sCode)
        #print "thermal(23)=%s elementType(3)=%s" %(self.thermal,self.elementType)


        ## assuming tCode=1
        if self.analysisCode==1:   # statics / displacement / heat flux
            self.addDataParameter(data,'lsdvmn','i',5,False)   ## load set number
        elif self.analysisCode==2: # real eigenvalues
            self.addDataParameter(data,'mode',     'i',5)         ## mode number
            self.addDataParameter(data,'eign',     'f',6,False)   ## real eigenvalue
            self.addDataParameter(data,'modeCycle','f',7,False)   ## mode or cycle @todo confused on the type - F1???
            self.applyDataCodeValue('dataNames',['mode','eigr','modeCycle'])
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
            self.addDataParameter(data,'lsdvmn','i',5)   ## load set
            self.applyDataCodeValue('dataNames',['lsdvmn'])
        elif self.analysisCode==8: # post-buckling
            self.addDataParameter(data,'lsdvmn','i',5)       ## mode number
            self.addDataParameter(data,'eigr','f',6,False)   ## real eigenvalue
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
            self.addDataParameter(data,'lsdvmn','i',5)   ## load set number
            self.applyDataCodeValue('dataNames',['lsdvmn'])
        elif self.analysisCode==12: # contran ? (may appear as aCode=6)  --> straight from DMAP...grrr...
            self.addDataParameter(data,'dt','f',5)   ## Time step ??? --> straight from DMAP
            self.applyDataCodeValue('dataNames',['lsdvmn'])
        else:
            raise InvalidAnalysisCodeError('invalid analysisCode...analysisCode=%s' %(self.analysisCode))
        ###
        # tCode=2
        #if self.analysisCode==2: # sort2
        #    self.lsdvmn = self.getValues(data,'i',5)

        self.readTitle()
        #print "n4 = ",self.n

    def parseStressCode(self):
        """
        sCode =  0 -> stressBits = [0,0,0,0,0]
        sCode =  1 -> stressBits = [0,0,0,0,1]
        sCode =  2 -> stressBits = [0,0,0,1,0]
        sCode =  3 -> stressBits = [0,0,0,1,1]
        etc.
        sCode = 32 -> stressBits = [1,1,1,1,1]

        stressBits[0] = 0 -> isMaxShear=True       isVonMises=False
        stressBits[0] = 1 -> isMaxShear=False      isVonMises=True

        stressBits[1] = 0 -> isStress=True         isStrain=False
        stressBits[2] = 0 -> isFiberCurvature=True isFiberDistance=False
        stressBits[3] = 0 -> duplicate of Bit[1] (stress/strain)
        stressBits[4] = 0 -> material coordinate system flag
        """
        bits = [0,0,0,0,0]
        
        sCode = self.sCode
        i=4
        while sCode>0:
            value = sCode%2
            sCode = (sCode - value)//2
            bits[i] = value
            #print "    *bit = ",value
            #print "    sCode = ",sCode
            i-=1
        #bits.reverse()
        #print "stressBits = ",bits
        self.stressBits = bits
        self.dataCode['stressBits'] = self.stressBits

    def readTable_OES_4_Data(self,iTable):
        isTable4Done = False
        isBlockDone  = False
        #print self.printSection(100)

        data = self.getData(16)
        #print self.printBlock(data) # on
        #print "16 block..."
        #self.printBlock(data)
        #self.printBlock(data)
        bufferWords, = unpack('i',data[4:8])
        #print "bufferWords = ",bufferWords
        if self.makeOp2Debug:
            self.op2Debug.write('bufferWords=|%s|\n' %(str(bufferWords)))

        #print "*********************"
        #bufferWords = self.getMarker() # 87 - buffer
        #print "OES4 bufferWords = ",bufferWords,bufferWords*4
        #self.verifyBufferSize(bufferWords)
        
        isBlockDone = not(bufferWords)
        #print "self.firstPass = ",self.firstPass
        
        ## table -4 is done, restarting table -3
        if self.isBufferDone: # table is done when the buffer is done
            isTable4Done = True
            #print "exitA"
            return isTable4Done,isBlockDone
        if bufferWords==0:
            #print "bufferWords 0 - done with Table4"
            isTable4Done = True
            #isBlockDone  = True
            self.printSection(40)
            #print "exitB"
            return isTable4Done,isBlockDone

        self.readElementTable()
        return isTable4Done,isBlockDone

    def readElementTable(self):
        #print "**self.readElementTable"
        #print "*elementType = ",self.elementType
        #print "op2.tell=%s n=%s" %(self.op2.tell(),self.n)
        
        self.rewind(4)
        self.data = self.readBlock()  # 348
        #print "len(self.data) = ",len(self.data)

        if self.makeOp2Debug:
            self.op2Debug.write('reading big data block\n')
        #print self.printBlock(self.data)

        #msg = 'elementType=%s -> %s' %(self.elementType,self.ElementType(self.elementType))
        tfsCode = [self.tableCode,self.formatCode,self.sortCode]

        self.parseStressCode()

        if not self.isValidSubcase():# lets the user skip a certain subcase
            self.log.debug("***skipping table=%s iSubcase=%s" %(self.tableName,self.iSubcase))
            self.skipOES_Element()
        elif self.thermal==0:
            # Stress / Strain
            self.dataCode['elementName'] = self.ElementType(self.elementType)
            if   tfsCode==[5,1,0]:
                self.readOES_Data_format1_sort0()
            #elif tfsCode==[5,1,1]:
                #self.readOES_Data_format1_sort1()
            #elif tfsCode==[5,2,0]:
                #self.readOES_Data_format2_sort0()
            #elif tfsCode==[5,2,1]:
                #self.readOES_Data_format2_sort1()
            #elif tfsCode==[5,3,0]:
                #self.readOES_Data_format3_sort0()
            #elif tfsCode==[5,3,1]:
                #self.readOES_Data_format3_sort1()
            #elif tfsCode==[5,3,2]:
                #self.readOES_Data_format3_sort2()
            else:
                #raise InvalidATFSCodeError('invalid atfsCode=%s' %(self.atfsCode))
                self.skipOES_Element()
                pass
            ###
        elif self.thermal==1:
            self.OES_Thermal()
            #raise NotImplementedError('thermal stress')
        else:
            raise Exception('invalid thermal option...')
        ###
        
        #print self.obj

    def readOES_Data_format1_sort1(self):
        #msg = '%s-OES format3_sort2 elementType=%-3s -> %-6s is not supported - fname=%s\n' %(self.tableName,self.elementType,self.ElementType(self.elementType),self.op2FileName)
        msg = ''
        if 0:
            pass
        else:
            #self.printBlock(self.data[0:100])
            self.skipOES_Element()
            msg = '%s-OES format1_sort1 elementType=%-3s -> %-6s is not supported - fname=%s\n' %(self.tableName,self.elementType,self.ElementType(self.elementType),self.op2FileName)
            self.log.debug(msg)
            #raise RuntimeError(msg)
        ###
        self.skippedCardsFile.write(msg)

    def readOES_Data_format2_sort0(self):
        #msg = '%s-OES format2_sort0 elementType=%-3s -> %-6s is not supported - fname=%s\n' %(self.tableName,self.elementType,self.ElementType(self.elementType),self.op2FileName)
        msg = ''
        if 0:
            pass
        else:
            #self.printBlock(self.data[0:100])
            self.skipOES_Element()
            msg = '%s-OES format2_sort0 elementType=%-3s -> %-6s is not supported - fname=%s\n' %(self.tableName,self.elementType,self.ElementType(self.elementType),self.op2FileName)
            self.log.debug(msg)
            #raise NotImplementedError(msg)
        ###
        self.skippedCardsFile.write(msg)

    def readOES_Data_format2_sort1(self):
        #msg = '%s-OES format2_sort1 elementType=%-3s -> %-6s is not supported - fname=%s\n' %(self.tableName,self.elementType,self.ElementType(self.elementType),self.op2FileName)
        msg = ''
        if 0:
            pass
        #if self.elementType in [1,3,10]: # crod/ctube/conrod
            #self.makeOES_Object(self.rodStress,rodStressObject,
            #                    self.rodStrain,rodStrainObject)
            #self.basicElement()
        #elif self.elementType in [10,11,12,33,74]:
            #raise AddNewElementError('add element=%s' %(self.elementType))
        else:
            #self.printBlock(self.data[0:100])
            self.skipOES_Element()
            msg = '%s-OES format2_sort1 elementType=%-3s -> %-6s is not supported - fname=%s\n' %(self.tableName,self.elementType,self.ElementType(self.elementType),self.op2FileName)
            self.log.debug(msg)
            #raise NotImplementedError(msg)
        ###
        self.skippedCardsFile.write(msg)

    def readOES_Data_format3_sort0(self):
        msg = '%s-OES format3_sort0 elementType=%-3s -> %-6s is not supported - fname=%s\n' %(self.tableName,self.elementType,self.ElementType(self.elementType),self.op2FileName)
        #msg = ''
        if 0:
            pass
        else:
            #self.printBlock(self.data[0:100])
            self.skipOES_Element()
            msg = '%s-OES format3_sort0 elementType=%-3s -> %-6s is not supported - fname=%s\n' %(self.tableName,self.elementType,self.ElementType(self.elementType),self.op2FileName)
            self.log.debug(msg)
            #raise RuntimeError(msg)
        ###
        self.skippedCardsFile.write(msg)

    def readOES_Data_format3_sort1(self):
        #msg = '%s-OES format3_sort1 elementType=%-3s -> %-6s is not supported - fname=%s\n' %(self.tableName,self.elementType,self.ElementType(self.elementType),self.op2FileName)
        msg = ''
        if 0:
            pass
        else:
            #self.printBlock(self.data[0:100])
            self.skipOES_Element()
            msg = '%s-OES format3_sort1 elementType=%-3s -> %-6s is not supported - fname=%s\n' %(self.tableName,self.elementType,self.ElementType(self.elementType),self.op2FileName)
            self.log.debug(msg)
            #raise NotImplementedError(msg)
        ###
        self.skippedCardsFile.write(msg)

    def readOES_Data_format3_sort2(self):
        #msg = '%s-OES format3_sort2 elementType=%-3s -> %-6s is not supported - fname=%s\n' %(self.tableName,self.elementType,self.ElementType(self.elementType),self.op2FileName)
        msg = ''
        if 0:
            pass
        else:
            #self.printBlock(self.data[0:100])
            self.skipOES_Element()
            msg = '%s-OES format3_sort2 elementType=%-3s -> %-6s is not supported - fname=%s\n' %(self.tableName,self.elementType,self.ElementType(self.elementType),self.op2FileName)
            self.log.debug(msg)
            #raise NotImplementedError(msg)
        ###
        self.skippedCardsFile.write(msg)

    def readOES_Data_format1_sort0(self):
        #msg = '%s-OES elementType=%-3s -> %-6s\n' %(self.tableName,self.elementType,self.ElementType(self.elementType))
        msg = ''
        #if self.analysisCode not in [1,6,10]:
            #raise InvalidATFSCodeError('self.atfsCode=%s' %(self.atfsCode))

        readCase = True
        if self.iSubcase in self.expectedTimes and len(self.expectedTimes[self.iSubcase])>0:
            readCase = self.updateDtMap()
        
        if readCase==False:
            self.skipOES_Element()
            return

        #print 'self.elementType  = ',self.elementType
        if self.elementType in [1,3,10]: # crod/ctube/conrod
            #if self.elementType==1:    self.dataCode['elementName'] = 'CROD'
            #if self.elementType==3:    self.dataCode['elementName'] = 'CTUBE'
            #if self.elementType==10:   self.dataCode['elementName'] = 'CONROD'
            
            self.makeOES_Object(self.rodStress,rodStressObject,
                                self.rodStrain,rodStrainObject)
            self.basicElement()
        elif self.elementType == 2:   # cbeam
            self.dataCode['elementName'] = 'CBEAM'
            self.makeOES_Object(self.beamStress,beamStressObject,
                                self.beamStrain,beamStrainObject)
            self.CBEAM_2()

        elif self.elementType in [4]: # cshear
            self.dataCode['elementName'] = 'CSHEAR'
            self.makeOES_Object(self.shearStress,shearStressObject,
                                self.shearStrain,shearStrainObject)
            self.basicElement()
        elif self.elementType in [11,12,13]:   # celas1/celas2/celas3
            #print "    found celas2_12"
            #if   self.elementType==11: self.dataCode['elementName'] = 'CELAS1'
            #elif self.elementType==12: self.dataCode['elementName'] = 'CELAS2'
            #elif self.elementType==13: self.dataCode['elementName'] = 'CELAS3'
            #else:  raise Exception('not implemented error')
            
            self.makeOES_Object(self.celasStress,celasStressObject,
                                self.celasStrain,celasStrainObject)
            self.basicElement()
        elif self.elementType == 34:   # cbar
            self.dataCode['elementName'] = 'CBAR'
            self.makeOES_Object(self.barStress,barStressObject,
                                self.barStrain,barStrainObject)
            self.CBAR_34()

        elif self.elementType==33: # cquad4_33
            self.dataCode['elementName'] = 'CQUAD4'
            self.makeOES_Object(self.plateStress,plateStressObject,
                                self.plateStrain,plateStrainObject)
            self.CQUAD4_33()
            #print self.obj.writeF06(['',''],'PAGE ',1)[0]

        elif self.elementType==53: # ctriax6
            self.dataCode['elementName'] = 'CTRIAX6'
            #print self.codeInformation()
            self.makeOES_Object(self.ctriaxStress,ctriaxStressObject,
                                self.ctriaxStrain,ctriaxStrainObject)
            self.CTRIAX6_53()

        elif self.elementType==74:  # ctria
            self.dataCode['elementName'] = 'CTRIA3'
            self.makeOES_Object(self.plateStress,plateStressObject,
                                self.plateStrain,plateStrainObject)
            self.CTRIA3_74()
        elif self.elementType in [64,144,70,75]: # cquad8/cquad4/ctriar/ctria6
            #print "    found cquad_144"
            if     self.elementType==64:  self.dataCode['elementName'] = 'CQUAD8'
            elif   self.elementType==144: self.dataCode['elementName'] = 'CQUAD4'
            elif   self.elementType==70:  self.dataCode['elementName'] = 'CTRIAR'
            elif   self.elementType==75:  self.dataCode['elementName'] = 'CTRIA6'
            else:  raise NotImplementedError('card not implemented elementType=%s' %(self.elementType))
            self.makeOES_Object(self.plateStress,plateStressObject,
                                self.plateStrain,plateStrainObject)
            self.CQUAD4_144()

        elif self.elementType in [39,67,68]:   # ctetra/chexa/cpenta (linear)
            #print "    found ctetra_39 / hexa_67 / cpenta_68"
            self.makeOES_Object(self.solidStress,solidStressObject,
                                self.solidStrain,solidStrainObject)
            self.CSOLID_67()

        elif self.elementType in [85]:   # ctetra/chexa/cpenta (91,93)  (nonlinear)
            #print "    found ctetra_85 / hexa_93 / cpenta_91"
            self.makeOES_Object(self.solidStress,solidStressObject,
                                self.solidStrain,solidStrainObject)
            self.CSOLID_85()

        elif self.elementType in [87,89,92]:   # CTUBENL, RODNL, CONRODNL
            #print "    found RODNL_89"
            self.makeOES_Object(self.nonlinearRodStress,nonlinearRodObject,
                                self.nonlinearRodStrain,nonlinearRodObject)
            self.RODNL_89_92()

        elif self.elementType in [88,90]: # CTRIA3NL, CQUAD4NL
            #print "cquad4_90"
            self.makeOES_Object(self.nonlinearPlateStress,nonlinearQuadObject,
                                self.nonlinearPlateStrain,nonlinearQuadObject)
            self.CQUAD4NL_90()

        #elif self.elementType in [91]: # CPENTANL
        #    #print "hexa_93"
        #    self.CPENTANL_91()
        #elif self.elementType in [93]: # CHEXANL
        #    #print "hexa_93"
        #    self.CHEXANL_93()

        elif self.elementType in [95,96,97,98]: # CQUAD4, CQUAD8, CTRIA3, CTRIA6 (composite)
            #print "    found a 95/96/97 or 98!"
            self.eid2 = None # stores the previous elementID
            self.makeOES_Object(self.compositePlateStress,compositePlateStressObject,
                                self.compositePlateStrain,compositePlateStrainObject)
            self.CQUAD4_95()
            del self.eid2

        #elif self.elementType in [94]: # CBEAM (nonlinear)
            #print "    found a 94!"
            #self.eid2 = None # stores the previous elementID
            #self.makeOES_Object(self.beamStress,beamStressObject,
            #                    self.beamStrain,beamStrainObject)
            #self.CBEAM_94()
            #raise NotImplementedError('stoping at end of CBEAM_94')
            #del self.eid2

        elif self.elementType in [139]:   # QUAD4FD (hyperelastic)
            #print "    found QUAD4FD_139"
            self.makeOES_Object(self.hyperelasticPlateStress,hyperelasticQuadObject,
                                self.hyperelasticPlateStrain,hyperelasticQuadObject)
            self.QUAD4FD_139()

        #elif self.elementType in [2,53,61,70,86,88,90,94,102,189,232,]:
            #self.skipOES_Element()
            #elementType=53  -> TRIAX6  is not supported
            #elementType=61  -> DUM9    is not supported
            #elementType=70  -> TRIAR   is not supported
            #elementType=86  -> GAPNL   is not supported
            #elementType=88  -> TRIA3NL is not supported
            #elementType=90  -> QUAD4NL is not supported
            #elementType=94  -> BEAMNL  is not supported
            #elementType=102 -> BUSH    is not supported
            #elementType=189 -> VUQUAD  is not supported
            #elementType=232 -> QUADRLC is not supported
        #elif self.elementType in [100]:   # BARS
        #    self.makeOES_Object(self.barStress,barStressObject,
        #                        self.barStrain,barStrainObject)
        #    self.basicElement()
        #elif self.elementType in [75,89,90,92,93]:
        #    msg = '%s-OES format1_sort0 elementType=%-3s -> %s is not supported - fname=%s\n' %(self.tableName,self.elementType,self.ElementType(self.elementType),self.op2FileName)
        #    raise AddNewElementError(msg)
        else:
            #self.printBlock(self.data[0:100])
            self.skipOES_Element()
            msg = '%s-OES format1_sort0 elementType=%-3s -> %s is not supported - fname=%s\n' %(self.tableName,self.elementType,self.ElementType(self.elementType),self.op2FileName.strip())
            self.log.debug(msg)
            self.skippedCardsFile.write(msg)
            #raise NotImplementedError(msg)
        ###
        #elif self.elementType == 1:    # crod     (done)
        #elif self.elementType == 2:    # cbeam    (done)
        #elif self.elementType == 3:    # ctube    (done)
        #elif self.elementType == 4:    # cshear   (done)
        #elif self.elementType == 10:   # conrod   (done)

        #elif self.elementType == 33:   # cquad4_33 (done)
        #elif self.elementType == 34:   # cbar      (done)
        #elif self.elementType == 35:   # cconeax
        #elif self.elementType == 38:   # cgap
        #elif self.elementType == 39:   # ctetra    (done)
        #elif self.elementType == 40:   # cbush1d

        #elif self.elementType == 47:   # caxif2
        #elif self.elementType == 48:   # caxif3
        #elif self.elementType == 49:   # caxif4
        #elif self.elementType == 50:   # cslot3
        #elif self.elementType == 51:   # cslot4
        #elif self.elementType == 53:   # ctriax6
        #elif self.elementType == 55:   # cdum3
        #elif self.elementType == 56:   # cdum4
        #elif self.elementType == 57:   # cdum5
        #elif self.elementType == 58:   # cdum6
        #elif self.elementType == 59:   # cdum7
        #elif self.elementType == 60:   # cdum8/crac2d
        #elif self.elementType == 64:   # cquad8   (done)
        #elif self.elementType == 67:   # chexa    (done)
        #elif self.elementType == 68:   # cpenta   (done)
        #elif self.elementType == 69:   # cbend
        #elif self.elementType == 70:   # ctriar
        #elif self.elementType == 74:   # ctria3  (done)
        #elif self.elementType == 75:   # ctria6  (done)
        #elif self.elementType == 82:   # cquadr
        #elif self.elementType == 95:   # composite CQUAD4    (done)
        #elif self.elementType == 96:   # composite CQUAD8    (done)
        #elif self.elementType == 97:   # composite CTRIA3    (done)
        #elif self.elementType == 98:   # composite CTRIA6    (done)
        #elif self.elementType == 100:  # cbar w/ cbarao or pload1

        #elif self.elementType == 102:  # cbush
        #elif self.elementType == 144:  # cquad_144 - corner stresses (done)


        # rods/bars/beams
        #elif self.elementType == 1:    # crod   (done)
        #elif self.elementType == 2:    # cbeam  (done)
        #elif self.elementType == 3:    # ctube  (done)
        #elif self.elementType == 10:   # conrod (done)
        #elif self.elementType == 34:   # cbar   (done)

        #springs
        #elif self.elementType == 11:   # celas1 (done)
        #elif self.elementType == 12:   # celas2 (done)
        #elif self.elementType == 13:   # celas3 (done)
        
        #plate
        #elif self.elementType == 33:   # cquad_33 (done)
        #elif self.elementType == 74:   # ctria3_74  - (done)
        #elif self.elementType == 75:   # ctria6_75  - (in progress)
        #elif self.elementType == 64:   # cquad8_64  - corner stresses (done)
        #elif self.elementType == 144:  # cquad4_144 - corner stresses (done)

        #solid (???)
        #elif self.elementType == 39:  # ctetra (done)
        #elif self.elementType == 67:  # chexa  (done)
        #elif self.elementType == 68:  # cpenta (done)

        # composite plate
        #elif self.elementType == 95: # CQUAD4 (done)
        #elif self.elementType == 96: # CQUAD8 (done)
        #elif self.elementType == 97: # CTRIA3 (done)
        #elif self.elementType == 98: # CTRIA6 (done)

        # nonlinear
        #elif self.elementType == 85:   # tetra  (nonlinear,not integrated)
        #elif self.elementType == 86:   # gap    (nonlinear)
        #elif self.elementType == 87:   # ctube  (nonlinear, done)
        #elif self.elementType == 88:   # tria3  (nonlinear) - same as quad4
        #elif self.elementType == 89:   # crod   (nonlinear, done)
        #elif self.elementType == 90:   # quad4  (nonlinear, done)
        #elif self.elementType == 91:   # cpenta (nonlinear)
        #elif self.elementType == 92:   # conrod (nonlinear, done)
        #elif self.elementType == 93:   # chexa  (nonlinear,not integrated)
        #elif self.elementType == 94:   # cbeam  (nonlinear)
        
        # acoustic
        #elif self.elementType == 76:   # chexa  (acoustic)
        #elif self.elementType == 77:   # cpenta (acoustic)
        #elif self.elementType == 78:   # ctetra (acoustic)
        #elif self.elementType == 101:  # caabsf (acoustic)

    def isStress(self):
        if self.stressBits[1]==0:
            return True
        return False

    def makeOES_Object(self,stress,stressObject,strain,strainObject):
        """
        Creates a stress/strain object if necessary
        """
        if self.isStress():
            self.createTransientObject(stress,stressObject)
        else:
            self.createTransientObject(strain,strainObject)
        ###
        #print "loading",self.obj.__class__.__name__
        return self.obj
