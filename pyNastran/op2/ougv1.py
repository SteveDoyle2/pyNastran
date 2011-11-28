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
        ## OUGV1
        tableName = self.readTableName(rewind=False) # OUGV1
        self.tableInit(tableName)
        print "tableName = |%r|" %(tableName)

        self.readMarkers([-1,7],'OUGV1')
        ints = self.readIntBlock()
        #print "*ints = ",ints

        self.readMarkers([-2,1,0],'OUGV1') # 7
        bufferWords = self.getMarker()
        self.op2Debug.write('bufferWords=%s\n' %(str(bufferWords)))
        #print "1-bufferWords = ",bufferWords,bufferWords*4
        ints = self.readIntBlock()
        self.op2Debug.write('ints=%s\n' %(str(ints)))
        
        markerA = -4
        markerB = 0

        #self.j = 0
        iTable = -3
        self.readMarkers([iTable,1,0],'OUGV1')
        while [markerA,markerB]!=[0,2]:
            self.readTable_OUGV1_3(iTable)
            isBlockDone = self.readTable_OUGV1_4(iTable-1)
            iTable -= 2

            if isBlockDone:
                #print "***"
                #print "iTable = ",iTable
                #print "$$$$"
                #self.n = self.markerStart
                #self.op2.seek(self.n)
                break
            ###

            n = self.n
            #self.printSection(100)
            self.readMarkers([iTable,1,0],'OUGV1')
            print "i read the markers!!!"
   
            #if self.j==3:
            #    print str(self.obj)
            #    sys.exit('check...j=%s dt=6E-2 dx=%s' %(self.j,'1.377e+01'))
            #self.j+=1
        ###
        self.readMarkers([iTable,1,0],'OUGV1')
        #self.printSection(100)
        #print str(self.obj)
        self.deleteAttributes_OUG()

    def deleteAttributes_OUG(self):
        params = ['lsdvm','mode','eigr','modeCycle','freq','dt','lftsfq','thermal','rCode','fCode','numWide','acousticFlag','thermal']
        self.deleteAttributes(params)
    
    def readTable_OUGV1_3(self,iTable): # iTable=-3
        bufferWords = self.getMarker()
        self.op2Debug.write('bufferWords=%s\n' %(str(bufferWords)))
        #print "2-bufferWords = ",bufferWords,bufferWords*4,'\n'

        data = self.getData(4)
        bufferSize, = unpack('i',data)
        data = self.getData(4*50)
        #print self.printBlock(data)
        
        
        aCode = self.getBlockIntEntry(data,1)
        print "aCode = ",aCode
        (three) = self.parseApproachCode(data)
        #iSubcase = self.getValues(data,'i',4)

        self.rCode   = self.getValues(data,'i',8)  ## random code
        self.fCode   = self.getValues(data,'i',9)  ## format code
        self.numWide = self.getValues(data,'i',10) ## number of words per entry in record; @note is this needed for this table ???
        self.acousticFlag = self.getValues(data,'f',13) ## acoustic pressure flag
        self.thermal      = self.getValues(data,'i',23) ## thermal flag; 1 for heat ransfer, 0 otherwise
        
        self.printBlock(data)
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
            print "LFTSFQ(5)=%s  EIGR(6)=%s  EIGI(7)=%s" %(self.lftsfq,self.eigr,self.eigi)
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

    def readTable_OUGV1_4(self,iTable):
        #self.readMarkers([iTable,1,0])
        markerA = 4
        
        while markerA>None:
            self.markerStart = copy.deepcopy(self.n)
            #self.printSection(180)
            self.readMarkers([iTable,1,0])
            print "starting OUGV1 table 4..."
            (isTable4Done,isBlockDone) = self.readTable_OUGV1_4_Data(iTable)
            if isTable4Done:
                print "done with OUGV1-4"
                self.n = self.markerStart
                self.op2.seek(self.n)
                break
            print "finished reading ougv1 table..."
            markerA = self.getMarker('A',debug=False)
            self.n-=12
            self.op2.seek(self.n)
            
            self.n = self.op2.tell()
            #print "***markerA = ",markerA
            #self.printSection(140)

            #print self.plateStress[self.iSubcase]
            
            iTable-=1
            print "isBlockDone = ",isBlockDone
        ###
        #print "isBlockDone = ",isBlockDone
        return isBlockDone

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
        if self.thermal==0:
            if self.approachCode==1 and self.sortCode==0: # displacement
                print "isDisplacement"
                self.obj = displacementObject(self.iSubcase)
                self.displacements[self.iSubcase] = self.obj
            elif self.approachCode==1 and self.sortCode==1: # spc forces
                print "isForces"
                raise Exception('is this correct???')
                self.obj = spcForcesObject(self.iSubcase)
                self.spcForces[self.iSubcase] = self.obj

            elif self.approachCode==2 and self.sortCode==0: # nonlinear static eigenvector
                print "isEigenvector"
                #self.obj = eigenVectorObject(self.iSubcase,self.eigr)
                self.createTransientObject(self.eigenvectors,eigenVectorObject,self.eigr)
                self.eigenvectors[self.iSubcase] = self.obj
                #print "****self", type(self.obj)

            elif self.approachCode==6 and self.sortCode==0: # transient displacement
                print "isTransientDisplacement"
                self.createTransientObject(self.displacments,displacementObject,self.dt)
                self.displacements[self.iSubcase] = self.obj

            elif self.approachCode==10 and self.sortCode==0: # nonlinear static displacement
                print "isNonlinearStaticDisplacement"
                self.createTransientObject(self.nonlinearDisplacments,displacementObject,self.dt)
                self.nonlinearDisplacements[self.iSubcase] = self.obj

            else:
                raise Exception('unsupported OUGV1 solution...')
            ###

        elif self.thermal==1:
            if self.approachCode==1 and self.sortCode==0: # temperature
                print "isTemperature"
                self.temperatures[self.iSubcase] = temperatureObject(self.iSubcase)

            elif self.approachCode==1 and self.sortCode==1: # heat fluxes
                print "isFluxes"
                self.obj = fluxObject(self.iSubcase,dt=self.dt)
                self.fluxes[self.iSubcase] = self.obj
            elif self.approachCode==6 and self.sortCode==0: # transient temperature
                print "isTransientTemperature"
                self.createTransientObject(self.temperatures,temperatureObject,self.dt)
                self.temperatures[self.iSubcase] = self.obj

            elif self.approachCode==10 and self.sortCode==0: # nonlinear static displacement
                print "isNonlinearStaticTemperatures"
                self.createTransientObject(self.nonlinearTemperatures,nonlinearTemperatureObject,self.lftsfq)
                self.nonlinearTemperatures[self.iSubcase] = self.obj
            else:
                raise Exception('unsupported OUGV1 thermal solution...')
            ###
        else:
            raise Exception('invalid thermal flag...not 0 or 1...flag=%s' %(self.thermal))
        ###
        self.readScalars(self.obj)
        #print self.obj
        
        print "-------finished OUGV1----------"
        return (isTable4Done,isBlockDone)

