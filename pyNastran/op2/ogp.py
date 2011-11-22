import sys
import copy
from struct import unpack

# pyNastran
#from pyNastran.op2.resultObjects.ougv1_Objects import (
#    temperatureObject,
#    nonlinearTemperatureObject,
#    fluxObject,nonlinearFluxObject)
from pyNastran.op2.resultObjects.opg_Objects import appliedLoadsObject


class OGP(object):
    """Table of element forces"""
    def readTable_OGP1(self):
        ## OGP
        tableName = self.readTableName(rewind=False) # OGP
        print "tableName = |%r|" %(tableName)

        self.readMarkers([-1,7],'OGP')
        ints = self.readIntBlock()
        print "*ints = ",ints

        self.readMarkers([-2,1,0],'OGP') # 7
        bufferWords = self.getMarker()
        print "1-bufferWords = ",bufferWords,bufferWords*4
        ints = self.readIntBlock()
        print "*ints = ",ints
        
        markerA = -4
        markerB = 0

        #self.j = 0
        iTable=-3
        self.readMarkers([iTable,1,0],'OGP')
        while [markerA,markerB]!=[0,2]:
            self.readTable_OGP_3(iTable)
            isBlockDone = self.readTable_OGP_4(iTable-1)
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
            self.readMarkers([iTable,1,0],'OGP')
            print "i read the markers!!!"
   
            #if self.j==3:
            #    print str(self.obj)
            #    sys.exit('check...j=%s dt=6E-2 dx=%s' %(self.j,'1.377e+01'))
            #self.j+=1

            #self.printSection(120)
        ###
        self.readMarkers([iTable,1,0],'OGP')
        #self.printSection(100)
        #print str(self.obj)
        self.deleteAttributes_OGP()

    def deleteAttributes_OGP(self):
        params = ['lsdvm','mode','eigr','eign','eigi','modeCycle','freq','time','lftsfq','dLoadID','fCode','numWide','oCode']
        self.deleteAttributes(params)

    def readTable_OGP_3(self,iTable): # iTable=-3
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
        self.numWide  = self.getValues(data,'i',10) ## number of words per entry in record; @note is this needed for this table ???
        self.oCode    = self.getValues(data,'i',11) ## undefined in DMAP...
        self.thermal  = self.getValues(data,'i',23) ## thermal flag; 1 for heat ransfer, 0 otherwise
        print "dLoadID(8)=%s fCode(9)=%s numWide(10)=%s oCode(11)=%s thermal(23)=%s" %(self.dLoadID,self.fCode,self.numWide,self.oCode,self.thermal)
        
        ## assuming tCode=1
        if self.approachCode==1:   # statics
            self.lsdvmn = self.getValues(data,'i',5) ## load set number
        elif self.approachCode==2: # normal modes/buckling (real eigenvalues)
            self.mode      = self.getValues(data,'i',5) ## mode number
            self.eign      = self.getValues(data,'f',6) ## real eigenvalue
            self.modeCycle = self.getValues(data,'f',7) ## mode or cycle @todo confused on the type ???
        elif self.approachCode==3: # differential stiffness 0
            self.lsdvmn = self.getValues(data,'i',5) ## load set number
        elif self.approachCode==4: # differential stiffness 1
            self.lsdvmn = self.getValues(data,'i',5) ## load set number
        elif self.approachCode==5:   # frequency
            self.freq = self.getValues(data,'f',5) ## frequency

        elif self.approachCode==6: # transient
            self.time = self.getValues(data,'f',5) ## time step
            print "TIME(5)=%s" %(self.time)
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
        #    self.lsdvmn = self.getValues(data,'i',5) ## load set, Mode number
        
        print "*iSubcase=%s"%(self.iSubcase)
        print "approachCode=%s tableCode=%s thermal=%s" %(self.approachCode,self.tableCode,self.thermal)
        print self.codeInformation()

        #self.printBlock(data)
        self.readTitle()

        #if self.j==3:
        #    #print str(self.obj)
        #    sys.exit('checkA...j=%s dt=6E-2 dx=%s dtActual=%f' %(self.j,'1.377e+01',self.dt))
        ###

    def readTable_OGP_4(self,iTable):
        #self.readMarkers([iTable,1,0])
        markerA = 4
        
        j = 0
        while markerA>None:
            self.markerStart = copy.deepcopy(self.n)
            #self.printSection(180)
            self.readMarkers([iTable,1,0])
            print "starting OGP table 4..."
            isTable4Done,isBlockDone = self.readTable_OGP_4_Data(iTable)
            if isTable4Done:
                print "done with OGP4"
                self.n = self.markerStart
                self.op2.seek(self.n)
                break
            print "finished reading ogp table..."
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

    def readTable_OGP_4_Data(self,iTable): # iTable=-4
        isTable4Done = False
        isBlockDone = False

        bufferWords = self.getMarker('OGP')
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
                self.obj = appliedLoadsObject(self.iSubcase)
                self.appliedLoads[self.iSubcase] = self.obj
            elif self.approachCode==1 and self.sortCode==1: # spc forces
                print "isForces"
                raise Exception('not implemented in OGP...')
                self.obj = spcForcesObject(self.iSubcase)
                self.spcForces[self.iSubcase] = self.obj
            elif self.approachCode==6 and self.sortCode==0: # transient displacement
                print "isTransientDisplacement"
                raise Exception('not implemented in OGP...')
                self.createTransientObject(self.displacments,displacementObject,self.time)
                self.displacements[self.iSubcase] = self.obj

            elif self.approachCode==10 and self.sortCode==0: # nonlinear static displacement
                print "isNonlinearStaticDisplacement"
                raise Exception('not implemented in OGP...')
                self.createTransientObject(self.nonlinearDisplacments,displacementObject,self.dt)
                self.nonlinearDisplacements[self.iSubcase] = self.obj
            else:
                raise Exception('not supported OGP solution...')
            ###

        elif self.thermal==1:
            if self.approachCode==1 and self.sortCode==0: # temperature
                print "isTemperature"
                raise Exception('not implemented in OGP...')
                self.temperatures[self.iSubcase] = temperatureObject(self.iSubcase)

            elif self.approachCode==1 and self.sortCode==1: # heat fluxes
                print "isFluxes"
                raise Exception('not implemented in OGP...')
                self.obj = fluxObject(self.iSubcase,dt=self.dt)
                self.fluxes[self.iSubcase] = self.obj
            elif self.approachCode==6 and self.sortCode==0: # transient temperature
                print "isTransientTemperature"
                raise Exception('not implemented in OGP...')
                self.createTransientObject(self.temperatureForces,temperatureObject,self.time)
                self.temperatureForces[self.iSubcase] = self.obj  ## @todo modify the name of this...
                #self.readForces(data,self.obj)

            elif self.approachCode==10 and self.sortCode==0: # nonlinear static displacement
                print "isNonlinearStaticTemperatures"
                raise Exception('not implemented in OGP...')
                self.createTransientObject(self.nonlinearFluxes,nonlinearFluxObject,self.lftsfq)
                self.nonlinearFluxes[self.iSubcase] = self.obj
                self.readForcesNonlinear(data,self.obj)
            else:
                raise Exception('not supported OGP solution...')
            ###
        else:
            raise Exception('invalid thermal flag...not 0 or 1...flag=%s' %(self.thermal))
        ###
        
        self.readOGPForces(data,self.obj)
        #print self.obj
        
        print "-------finished OGP----------"
        return (isTable4Done,isBlockDone)

        
    def readOGPForces(self,data,scalarObject):
        deviceCode = self.deviceCode
        #self.printBlock(data[0:self.numWide*4])
        dn = self.numWide*4
        while data:
            #print "len(data) = ",len(data)
            eData = data[0:dn]
            #self.printBlock(data[:dn])
            gridDevice,eid = unpack('ii',data[0:8])
            nodeID = (gridDevice-deviceCode)/10
            
            source = ''.join(unpack('cccccccc',data[8:16]))
            #self.printBlock(data[16:dn])
            (dx,dy,dz,rx,ry,rz) = unpack('ffffff',data[16:40])
            #print "source = |%s|" %(source)
            #print "len(data[8:40]"
            
            print "nodeID=%s eid=%s source=|%s| dx=%-4i dy=%-4i dz=%-4i rx=%-4i ry=%-4i rz=%-4i" %(nodeID,eid,source,dx,dy,dz,rx,ry,rz)
            scalarObject.add(nodeID,eid,source,dx,dy,dz,rx,ry,rz)
            #print "gridDevice = ",gridDevice
            #print "deviceCode = ",deviceCode
            grid = (gridDevice-self.deviceCode)/10
            #print "grid=%g dx=%g dy=%g dz=%g rx=%g ry=%g rz=%g" %(grid,xGrad,yGrad,zGrad,xFlux,yFlux,zFlux)
            #print type(scalarObject)
            data = data[dn:]
            #sys.exit('asd')
        ###
        print "***********"
        #print self.obj
        #sys.exit('check...')

