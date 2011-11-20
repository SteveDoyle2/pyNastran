import sys
import copy
from struct import unpack

# pyNastran
from ougv1_Objects import (temperatureObject,displacementObject,
                           nonlinearTemperatureObject,
                           fluxObject,nonlinearFluxObject)

class OEF(object):
    """Table of element forces"""
    def readTable_OEF(self):
        ## OEF
        tableName = self.readTableName(rewind=False) # OEF
        print "tableName = |%r|" %(tableName)

        self.readMarkers([-1,7],'OEF')
        ints = self.readIntBlock()
        print "*ints = ",ints

        self.readMarkers([-2,1,0],'OEF') # 7
        bufferWords = self.getMarker()
        print "1-bufferWords = ",bufferWords,bufferWords*4
        ints = self.readIntBlock()
        print "*ints = ",ints
        
        markerA = -4
        markerB = 0

        #self.j = 0
        iTable=-3
        self.readMarkers([iTable,1,0],'OEF')
        while [markerA,markerB]!=[0,2]:
            self.readTable_OEF_3(iTable)
            isBlockDone = self.readTable_OEF_4(iTable-1)
            iTable -= 2

            if isBlockDone:
                print "***"
                print "iTable = ",iTable
                print "$$$$"
                #self.n = self.markerStart
                #self.op2.seek(self.n)
                break
            ###

            n = self.n
            #self.printSection(100)

            self.readMarkers([iTable,1,0],'OEF')
            print "i read the markers!!!"
   
            #if self.j==3:
            #    print str(self.obj)
            #    sys.exit('check...j=%s dt=6E-2 dx=%s' %(self.j,'1.377e+01'))
            #self.j+=1

            #self.printSection(120)
        ###
        self.readMarkers([iTable,1,0],'OEF')
        #self.printSection(100)
        print str(self.obj)

    def readTable_OEF_3(self,iTable): # iTable=-3
        bufferWords = self.getMarker()
        print "2-bufferWords = ",bufferWords,bufferWords*4,'\n'

        data = self.getData(4)
        bufferSize, = unpack('i',data)
        data = self.getData(4*50)
        
        #self.printBlock(data)
        
        
        aCode = self.getBlockIntEntry(data,1)
        print "aCode = ",aCode
        (three) = self.parseApproachCode(data)
        #iSubcase = self.getValues(data,'i',4)

        self.dLoadID  = self.getValues(data,'i',8)  ## dynamic load set ID/random code
        self.fCode    = self.getValues(data,'i',9)  ## format code
        self.numwide  = self.getValues(data,'i',10) ## number of words per entry in record; @note is this needed for this table ???
        self.oCode    = self.getValues(data,'i',11) ## undefined in DMAP...
        self.thermal  = self.getValues(data,'i',23) ## thermal flag; 1 for heat ransfer, 0 otherwise
        print "dLoadID(8)=%s fCode(9)=%s numwde(10)=%s oCode(11)=%s thermal(23)=%s" %(self.dLoadID,self.fCode,self.numwide,self.oCode,self.thermal)
        
        ## assuming tCode=1
        if self.approachCode==1:   # statics
            self.loadID = self.getValues(data,'i',5) ## load set ID number
        elif self.approachCode==2: # normal modes/buckling (real eigenvalues)
            self.mode      = self.getValues(data,'i',5) ## mode number
            self.eign      = self.getValues(data,'f',6) ## eigenvalue
        elif self.approachCode==3: # differential stiffness 0
            self.loadID = self.getValues(data,'i',5) ## load set ID number
        elif self.approachCode==4: # differential stiffness 1
            self.loadID = self.getValues(data,'i',5) ## load set ID number
        elif self.approachCode==5:   # frequency
            self.freq = self.getValues(data,'f',5) ## frequency

        elif self.approachCode==6: # transient
            self.time = self.getValues(data,'f',5) ## time step
            print "TIME(5)=%s" %(self.time)
        elif self.approachCode==7: # pre-buckling
            self.loadID = self.getValues(data,'i',5) ## load set ID number
            print "LOADID(5)=%s" %(self.loadID)
        elif self.approachCode==8: # post-buckling
            self.loadID = self.getValues(data,'i',5) ## load set ID number
            self.eigr   = self.getValues(data,'f',6) ## real eigenvalue
            print "LOADID(5)=%s  EIGR(6)=%s" %(self.loadID,self.eigr)
        elif self.approachCode==9: # complex eigenvalues
            self.mode   = self.getValues(data,'i',5) ## mode
            self.eigr   = self.getValues(data,'f',6) ## real eigenvalue
            self.eigi   = self.getValues(data,'f',7) ## imaginary eigenvalue
            print "MODE(5)=%s  EIGR(6)=%s  EIGI(7)=%s" %(self.mode,self.eigr,self.eigi)
        elif self.approachCode==10: # nonlinear statics
            self.loadStep = self.getValues(data,'f',5) ## load step
            print "loadStep(5) = %s" %(self.loadStep)
        elif self.approachCode==11: # geometric nonlinear statics
            self.loadID = self.getValues(data,'i',5) ## load set ID number
            print "LOADID(5)=%s" %(self.loadID)
        else:
            raise RuntimeError('invalid approach code...approachCode=%s' %(self.approachCode))

        # tCode=2
        #if self.analysisCode==2: # sort2
        #    self.loadID = self.getValues(data,'i',5) ## load set ID number
        
        print "*iSubcase=%s"%(self.iSubcase)
        print "approachCode=%s tableCode=%s thermal=%s" %(self.approachCode,self.tableCode,self.thermal)
        print self.codeInformation(sCode=None,tCode=None,thermal=self.thermal)

        #self.printBlock(data)
        self.readTitle()

        #return (analysisCode,tableCode,thermal)

        #if self.j==3:
        #    #print str(self.obj)
        #    sys.exit('checkA...j=%s dt=6E-2 dx=%s dtActual=%f' %(self.j,'1.377e+01',self.dt))
        ###

    def readTable_OEF_4(self,iTable):
        #self.readMarkers([iTable,1,0])
        markerA = 4
        
        j = 0
        while markerA>None:
            self.markerStart = copy.deepcopy(self.n)
            #self.printSection(180)
            self.readMarkers([iTable,1,0])
            print "starting OEF table 4..."
            isTable4Done,isBlockDone = self.readTable_OEF_4_Data(iTable)
            if isTable4Done:
                print "done with OEF4"
                self.n = self.markerStart
                self.op2.seek(self.n)
                break
            print "finished reading oef table..."
            markerA = self.getMarker('A')
            self.n-=12
            self.op2.seek(self.n)
            
            self.n = self.op2.tell()
            print "***markerA = ",markerA
            #self.printSection(140)

            #print self.plateStress[self.iSubcase]
            
            try:
                del self.analysisCode
                del self.tableCode
                del self.thermal
                #del self.dt
            except:
                pass
            iTable-=1
            if j>10000:
                sys.exit('check...')
            j+=1
            print "isBlockDone = ",isBlockDone
            
        print "isBlockDone = ",isBlockDone
        return isBlockDone

    def readTable_OEF_4_Data(self,iTable): # iTable=-4
        isTable4Done = False
        isBlockDone = False

        bufferWords = self.getMarker('OEF')
        #print len(bufferWords)
        data = self.readBlock()
        #self.printBlock(data)

        if bufferWords==146:  # table -4 is done, restarting table -3
            isTable4Done = True
            return isTable4Done,isBlockDone
        elif bufferWords==0:
            print "bufferWords 0 - done with Table4"
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
            elif self.approachCode==6 and self.sortCode==0: # transient displacement
                print "isTransientDisplacement"

                self.createTransientObject(self.displacments,displacementObject,self.dt)
                self.displacements[self.iSubcase] = self.obj

            elif self.approachCode==10 and self.sortCode==0: # nonlinear static displacement
                print "isNonlinearStaticDisplacement"
                self.createTransientObject(self.nonlinearDisplacments,displacementObject,self.dt)
                self.nonlinearDisplacements[self.iSubcase] = self.obj
            else:
                raise Exception('not supported OEF solution...')
            ###

        elif self.thermal==1:
            if self.approachCode==1 and self.sortCode==0: # temperature
                print "isTemperature"
                raise Exception('verify...')
                self.temperatures[self.iSubcase] = temperatureObject(self.iSubcase)

            elif self.approachCode==1 and self.sortCode==1: # heat fluxes
                print "isFluxes"
                raise Exception('verify...')
                self.obj = fluxObject(self.iSubcase,dt=self.dt)
                self.fluxes[self.iSubcase] = self.obj
            elif self.approachCode==6 and self.sortCode==0: # transient temperature
                print "isTransientTemperature"
                #raise Exception('verify...')
                self.createTransientObject(self.temperatures,temperatureObject,self.dt)
                self.temperatures[self.iSubcase] = self.obj  ## @todo modify the name of this...

            elif self.approachCode==10 and self.sortCode==0: # nonlinear static displacement
                print "isNonlinearStaticTemperatures"
                self.createTransientObject(self.nonlinearFluxes,nonlinearFluxObject,self.lftsfq)
                self.nonlinearFluxes[self.iSubcase] = self.obj
                self.readForcesNonlinear(data,self.obj)
            else:
                raise Exception('not supported OEF solution...')
            ###
        else:
            raise Exception('invalid thermal flag...not 0 or 1...flag=%s' %(self.thermal))
        ###
        
        self.readForces(data,self.obj)
        #print self.obj
        
        print "-------finished OEF----------"
        return (isTable4Done,isBlockDone)

        
    def readForces(self,data,scalarObject):
        #self.printBlock(data[0:self.numwide*4])
        while data:
            #print "len(data) = ",len(data)
            #self.printBlock(data[32:])
            gridDevice, = unpack('i',data[0:4])
            eType = ''.join(unpack('cccccccc',data[4:12]))
            #print "eType = ",eType
            #print "len(data[8:40]"
            if self.numwide in [9,10]:
                (xGrad,yGrad,zGrad,xFlux,yFlux,zFlux) = unpack('ffffff',data[12:36])
            elif self.numwide==8: ## @todo CHBDY - how do i add this to the case...
                (fApplied,freeConv,forceConv,fRad,fTotal) = unpack('fffff',data[12:32])
                sys.stderr.write('skipping CHBDY')
                data = data[self.numwide*4:]
                continue
            else:
                raise Exception('only CBEAM/CBAR/CTUBE/2D/3D elements supported...so no special thermal elements...numwde=%s' %(self.numwide))
            
            #print "gridDevice = ",gridDevice
            #print "deviceCode = ",deviceCode
            grid = (gridDevice-self.deviceCode)/10
            #print "grid=%g dx=%g dy=%g dz=%g rx=%g ry=%g rz=%g" %(grid,xGrad,yGrad,zGrad,xFlux,yFlux,zFlux)
            #print type(scalarObject)
            scalarObject.add(grid,xGrad,yGrad,zGrad,xFlux,yFlux,zFlux)
            data = data[self.numwide*4:]
        ###
        #print self.obj
        #sys.exit('check...')

    def readForcesNonlinear(self,data,scalarObject):
        while data:
            #print "len(data) = ",len(data)
            #self.printBlock(data[32:])
            gridDevice, = unpack('i',data[0:4])
            eType = ''.join(unpack('cccccccc',data[4:12]))
            #print "eType = ",eType
            #print "len(data[8:40]"
            if self.numwide in [9,10]:
                (xGrad,yGrad,zGrad,xFlux,yFlux,zFlux) = unpack('ffffff',data[12:36])
            elif self.numwide==8: ## @todo CHBDY - how do i add this to the case...
                (fApplied,freeConv,forceConv,fRad,fTotal) = unpack('fffff',data[12:32])
                sys.stderr.write('skipping CHBDY')
                data = data[self.numwide*4:]
                continue
            else:
                raise Exception('only CBEAM/CBAR/CTUBE/2D/3D elements supported...so no special thermal elements...numwde=%s' %(self.numwide))
            
            #print "gridDevice = ",gridDevice
            #print "deviceCode = ",deviceCode
            grid = (gridDevice-self.deviceCode)/10
            #print "grid=%g dx=%g dy=%g dz=%g rx=%g ry=%g rz=%g" %(grid,xGrad,yGrad,zGrad,xFlux,yFlux,zFlux)
            print type(scalarObject)
            scalarObject.add(grid,eType,xGrad,yGrad,zGrad,xFlux,yFlux,zFlux)
            data = data[self.numwide*4:]
        ###
        #print self.obj
        #sys.exit('check...')
        
    def isDisplacement(self):
        if self.approachCode==1 and self.thermal==0:
            return True
        return False

    def isTransientDisplacement(self):
        if self.approachCode==6 and self.sortCode==0 and self.thermal==0:
            return True
        return False

    def isTemperature(self):
        if self.approachCode==1 and self.sortCode==0 and self.thermal==1:
            return True
        return False

    def isTransientTemperature(self):
        if self.approachCode==6 and self.sortCode==0 and self.thermal==1:
            return True
        return False

    def isForces(self):
        if(approachCode==1 and self.sortCode==1 and self.thermal==0):
            return True
        return False

    def isFluxes(self):
        if(self.approachCode==1 and self.sortCode==1 and self.thermal==1):
            return True
        return False
