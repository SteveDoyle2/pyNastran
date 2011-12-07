import sys
import copy
from struct import unpack

# pyNastran
from pyNastran.op2.resultObjects.ougv1_Objects import (
     temperatureObject,displacementObject,  # approachCode=1, sortCode=0
     eigenVectorObject,                     # approachCode=2, sortCode=0
     fluxObject,                            # approachCode=1, sortCode=3
     nonlinearTemperatureObject,            # approachCode=10,sortCode=0
     )

class OUGV1(object):
    """Table of displacements/velocities/acceleration/heat flux/temperature"""

    def readTable_OUG1(self):
        table3 = self.readTable_OUGV1_3
        table4Data = self.readTable_OUGV1_4_Data
        self.readResultsTable(table3,table4Data)
        self.deleteAttributes_OUG()

    def deleteAttributes_OUG(self):
        params = ['lsdvm','mode','eigr','modeCycle','freq','dt','lftsfq','thermal','rCode','fCode','numWide','acousticFlag','thermal']
        self.deleteAttributes(params)
    
    def readTable_OUGV1_3(self,iTable): # iTable=-3
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
        #iSubcase = self.getValues(data,'i',4)

        self.rCode        = self.getValues(data,'i',8)  ## random code
        self.formatCode   = self.getValues(data,'i',9)  ## format code
        self.numWide      = self.getValues(data,'i',10) ## number of words per entry in record; @note is this needed for this table ???
        self.acousticFlag = self.getValues(data,'f',13) ## acoustic pressure flag
        self.thermal      = self.getValues(data,'i',23) ## thermal flag; 1 for heat ransfer, 0 otherwise
        
        #self.printBlock(data) # on
        ## assuming tCode=1
        if self.approachCode==1:   # statics / displacement / heat flux
            self.lsdvmn = self.getValues(data,'i',5) ## load set number
        elif self.approachCode==2: # real eigenvalues
            self.mode      = self.getValues(data,'i',5) ## mode number
            self.eigr      = self.getValues(data,'f',6) ## real eigenvalue
            self.modeCycle = self.getValues(data,'f',7) ## mode or cycle @todo confused on the type - F1???
            print "mode(5)=%s eigr(6)=%s modeCycle(7)=%s" %(self.mode,self.eigr,self.modeCycle)
        elif self.approachCode==3: # differential stiffness
            self.lsdvmn = self.getValues(data,'i',5) ## load set number
        elif self.approachCode==4: # differential stiffness
            self.lsdvmn = self.getValues(data,'i',5) ## load set number
        elif self.approachCode==5:   # frequency
            self.freq = self.getValues(data,'f',5) ## frequency

        elif self.approachCode==6: # transient
            self.dt = self.getValues(data,'f',5) ## time step
            print "DT(5)=%s" %(self.dt)
        elif self.approachCode==7: # pre-buckling
            self.lsdvmn = self.getValues(data,'i',5) ## load set
            print "LSDVMN(5)=%s" %(self.lsdvmn)
        elif self.approachCode==8: # post-buckling
            self.lsdvmn = self.getValues(data,'i',5) ## mode number
            self.eigr   = self.getValues(data,'f',6) ## real eigenvalue
            print "LSDVMN(5)=%s  EIGR(6)=%s" %(self.lsdvmn,self.eigr)
        elif self.approachCode==9: # complex eigenvalues
            self.mode   = self.getValues(data,'i',5) ## mode
            self.eigr   = self.getValues(data,'f',6) ## real eigenvalue
            self.eigi   = self.getValues(data,'f',7) ## imaginary eigenvalue
            print "mode(5)=%s  eigr(6)=%s  eigi(7)=%s" %(self.mode,self.eigr,self.eigi)
        elif self.approachCode==10: # nonlinear statics
            self.lftsfq = self.getValues(data,'f',5) ## load step
            print "LFTSFQ(5) = %s" %(self.lftsfq)
        elif self.approachCode==11: # old geometric nonlinear statics
            self.lsdvmn = self.getValues(data,'i',5)
            print "LSDVMN(5)=%s" %(self.lsdvmn)
        elif self.approachCode==12: # contran ? (may appear as aCode=6)  --> straight from DMAP...grrr...
            self.lsdvmn = self.getValues(data,'i',5)
            print "LSDVMN(5)=%s" %(self.lsdvmn)
        else:
            raise RuntimeError('invalid approach code...approachCode=%s' %(self.approachCode))
        # tCode=2
        #if self.analysisCode==2: # sort2
        #    self.lsdvmn = self.getValues(data,'i',5)
        
        print "*iSubcase=%s"%(self.iSubcase)
        print "approachCode=%s tableCode=%s thermal=%s" %(self.approachCode,self.tableCode,self.thermal)
        print self.codeInformation()

        #self.printBlock(data)
        self.readTitle()

        #return (analysisCode,tableCode,thermal)

        #if self.j==3:
        #    #print str(self.obj)
        #    sys.exit('checkA...j=%s dt=6E-2 dx=%s dtActual=%f' %(self.j,'1.377e+01',self.dt))
        ###

    def readTable_OUGV1_4_Data(self,iTable): # iTable=-4
        isTable4Done = False
        isBlockDone  = False

        bufferWords = self.getMarker('OUGV1')
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
        print "self.approachCode=%s tableCode(1)=%s thermal(23)=%g" %(self.approachCode,self.tableCode,self.thermal)

        self.readOUGV1_Data()
        #print "-------finished OUGV1----------"
        return (isTable4Done,isBlockDone)

    def readOUGV1_Data(self):
        tfsCode = [self.tableCode,self.formatCode,self.sortCode]
        #if self.thermal==2:
        #    self.skipOES_Element(None)
        if 0:
            pass
        else:
            if   tfsCode in [[1,1,0],[7,1,0]]: # displacement/modes
                self.readOUGV1_Data_format1_sort0()
            elif tfsCode in [[1,1,1],[7,1,1]]:
                self.readOUGV1_Data_format1_sort1()
            #elif tfsCode==[1,2,1]:
            #    self.readOUGV1_Data_format2_sort1()
            elif tfsCode==[1,3,1]:
                self.readOUGV1_Data_format3_sort1()
            else:
                raise Exception('bad tableCode/formatCode/sortCode=%s' %(tfsCode))
            ###
        ###

    def readOUGV1_Data_format3_sort1(self):
        assert self.formatCode==3 # Magnitude/Phase
        assert self.sortCode==1   # Imaginary
        if self.thermal==0:
            if self.approachCode==5: # frequency displacement
                print "isFreqDisplacement"
                self.createTransientObject(self.freqDisplacements,displacementObject,self.freq)
                self.freqDisplacements[self.iSubcase] = self.obj
            else:
                raise Exception('unsupported OUGV1 static solution...')
            ###
        elif self.thermal==1:
            raise Exception('unsupported OUGV1 thermal solution...')            
        else:
            raise Exception('invalid thermal flag...not 0 or 1...flag=%s' %(self.thermal))
        ###
        self.readScalars14(self.obj)

    def readOUGV1_Data_format1_sort0(self):
        assert self.formatCode==1 # Real
        assert self.sortCode==0   # Real
        if self.thermal==0:
            if self.approachCode==1: # displacement
                print "isDisplacement"
                self.obj = displacementObject(self.iSubcase)
                self.displacements[self.iSubcase] = self.obj
            elif self.approachCode==1: # spc forces
                print "isForces"
                raise Exception('is this correct???')
                self.obj = spcForcesObject(self.iSubcase)
                self.spcForces[self.iSubcase] = self.obj
                
            elif self.approachCode==2: # nonlinear static eigenvector
                print "isEigenvector"
                #self.obj = eigenVectorObject(self.iSubcase,self.eigr)
                self.createTransientObject(self.eigenvectors,eigenVectorObject,self.eigr)
                self.eigenvectors[self.iSubcase] = self.obj
                #print "****self", type(self.obj)

            elif self.approachCode==5: # frequency displacement
                print "isFrequencyDisplacement"
                self.createTransientObject(self.freqDisplacements,eigenVectorObject,self.freq)
                self.freqDisplacements[self.iSubcase] = self.obj

            elif self.approachCode==6: # transient displacement
                #print "isTransientDisplacement"
                #self.createTransientObject(self.displacements,displacementObject,self.dt)
                #self.displacements[self.iSubcase] = self.obj
                self.obj = None

            elif self.approachCode==7: # pre-buckling displacement
                print "isPreBucklingDisplacement"
                self.createTransientObject(self.preBucklingDisplacements,displacementObject,self.lsdvmn)
                self.preBucklingDisplacements[self.iSubcase] = self.obj

            elif self.approachCode==8: # post-buckling eigenvector
                print "isPostBucklingEigenvector"
                #print "pbe = ",self.postBucklingEigenvector

                #self.obj = eigenVectorObject(self.iSubcase,self.eigr)
                self.createTransientObject(self.postBucklingEigenvector,eigenVectorObject,self.eigr)
                self.postBucklingEigenvector[self.iSubcase] = self.obj
                #print "****self", type(self.obj)

            elif self.approachCode==9: # nonlinear static eigenvector
                print "isComplexEigenvalues"
                #self.obj = eigenVectorObject(self.iSubcase,self.eigr)
                self.createTransientObject(self.complexEigenvalues,eigenVectorObject,(self.mode,self.eigr,self.eigi))
                self.complexEigenvalues[self.iSubcase] = self.obj
                #print "****self", type(self.obj)

            elif self.approachCode==10: # nonlinear static displacement
                print "isNonlinearStaticDisplacement"
                self.createTransientObject(self.nonlinearDisplacements,displacementObject,self.lftsfq)
                self.nonlinearDisplacements[self.iSubcase] = self.obj
            else:
                raise Exception('unsupported OUGV1 static solution...')
            ###

        elif self.thermal==1:
            if self.approachCode==1: # temperature
                print "isTemperature"
                self.temperatures[self.iSubcase] = temperatureObject(self.iSubcase)

            #elif self.approachCode==1 and self.sortCode==1: # heat fluxes
            #    print "isFluxes"
            #    self.obj = fluxObject(self.iSubcase,dt=self.dt)
            #    self.fluxes[self.iSubcase] = self.obj
            elif self.approachCode==6: # transient temperature
                print "isTransientTemperature"
                self.createTransientObject(self.temperatures,temperatureObject,self.dt)
                self.temperatures[self.iSubcase] = self.obj

            elif self.approachCode==10: # nonlinear static displacement
                print "isNonlinearStaticTemperatures"
                self.createTransientObject(self.nonlinearTemperatures,nonlinearTemperatureObject,self.lftsfq)
                self.nonlinearTemperatures[self.iSubcase] = self.obj
            else:
                raise Exception('unsupported OUGV1 thermal solution...')
            ###
        else:
            raise Exception('invalid thermal flag...not 0 or 1...flag=%s' %(self.thermal))
        ###
        if self.obj:
            self.readScalars8(self.obj)
        else:
            self.skipOES_Element(None)
        ###
        #print self.obj
        #return

    def readOUGV1_Data_format1_sort1(self):
        assert self.formatCode==1 # Real
        assert self.sortCode==1   # Real/Imaginary
        if self.thermal==0:
            if self.approachCode==5: # frequency displacement
                print "isFrequencyDisplacement"
                self.createTransientObject(self.freqDisplacements,eigenVectorObject,self.freq)
                self.freqDisplacements[self.iSubcase] = self.obj
            elif self.approachCode==9: # nonlinear static eigenvector
                print "isComplexEigenvalues"
                #self.obj = eigenVectorObject(self.iSubcase,self.eigr)
                self.createTransientObject(self.complexEigenvalues,eigenVectorObject,(self.mode,self.eigr,self.eigi))
                self.complexEigenvalues[self.iSubcase] = self.obj
                #print "****self", type(self.obj)
            else:
                raise Exception('unsupported OUGV1 static format1_sort1 solution...')
            ###
        else:
            raise Exception('invalid thermal flag...not 0 or 1...flag=%s' %(self.thermal))
        ###
        self.readScalars14(self.obj)
        #print self.obj
        #return

    def readOUGV1_Data_format2_sort1(self):
        assert self.formatCode==2 # Real/Imaginary
        assert self.sortCode==1   # Real/Imaginary
        if self.thermal==0:
            if self.approachCode==9 and self.sortCode==1: # nonlinear static eigenvector
                print "isComplexEigenvalues"
                #self.obj = eigenVectorObject(self.iSubcase,self.eigr)
                self.createTransientObject(self.complexEigenvalues,eigenVectorObject,(self.mode,self.eigr,self.eigi))
                self.complexEigenvalues[self.iSubcase] = self.obj
                #print "****self", type(self.obj)
            else:
                raise Exception('unsupported OUGV1 static solution...')
            ###
        else:
            raise Exception('invalid thermal flag...not 0 or 1...flag=%s' %(self.thermal))
        ###
        self.readScalars8(self.obj)
        #print self.obj
        #return


