import sys
import copy
from struct import unpack

# pyNastran
from op2_Objects import temperatureObject,fluxObject,displacementObject

class OUGV1(object):

    def getValues(self,data,sFormat,iWordStart,iWordStop=None):
        """
        extracts the ith word from the data structure as the provided type
        supports multiple inputs with iWordStop (note this is words, not outputs)
        @warning
            works with nastran syntax, not standard python syntax
            this makes it work with what's documented in the DMAP manual
        """
        if iWordStop==None:
            print "iWordStart=%s data[%s:%s]" %(iWordStart,iWordStart*4,(iWordStart+1)*4)
            ds = data[(iWordStart-1)*4:iWordStart*4]
            return unpack(sFormat,ds)[0]
            
        #print "type(data) = ",type(data)
        ds = data[(iWordStart-1)*4:(iWordStop-1)*4]
        return unpack(sFormat,ds)
        
    def getValues8(self,data,sFormat,iWordStart,iWordStop=None):
        if iWordStop==None:
            ds = data[iWordStart*8:(iWordStart+1)*8]
            return unpack(sFormat,ds)[0]
            
        #print "type(data) = ",type(data)
        ds = data[iWordStart*8:iWordStop*8]
        return unpack(sFormat,ds)
        
    def readTable_OUGV1(self):
        ## OUGV1
        tableName = self.readTableName(rewind=False) # OUGV1
        print "tableName = |%r|" %(tableName)

        self.readMarkers([-1,7],'OUGV1')
        ints = self.readIntBlock()
        print "*ints = ",ints

        self.readMarkers([-2,1,0],'OUGV1') # 7
        bufferWords = self.getMarker()
        print "1-bufferWords = ",bufferWords,bufferWords*4
        ints = self.readIntBlock()
        print "*ints = ",ints
        
        markerA = -4
        markerB = 0

        #self.j = 0
        iTable=-3
        self.readMarkers([iTable,1,0],'OUGV1')
        while [markerA,markerB]!=[0,2]:
            self.readTable_OUGV1_3(iTable)
            isBlockDone = self.readTable_OUGV1_4(iTable-1)
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

            #markerA = self.getMarker('A')
            #markerB = self.getMarker('B')
            #self.n-=24
            #self.op2.seek(self.n)
            #print "markerA=%s markerB=%s" %(markerA,markerB)
            self.readMarkers([iTable,1,0],'OUGV1')
            print "i read the markers!!!"
   
            #if self.j==3:
            #    print str(self.obj)
            #    sys.exit('check...j=%s dt=6E-2 dx=%s' %(self.j,'1.377e+01'))
            #self.j+=1

            #self.printSection(120)
            #break
        self.readMarkers([iTable,1,0],'OUGV1')
        #self.printSection(100)
        print str(self.obj)

        #sys.exit('end of displacementA')

    def readTable_OUGV1_3(self,iTable): # iTable=-3
        bufferWords = self.getMarker()
        print "2-bufferWords = ",bufferWords,bufferWords*4,'\n'

        data = self.getData(4)
        bufferSize, = unpack('i',data)
        print "bufferSize = ",bufferSize
        data = self.getData(4*50)
        
        #print "---dataBlock 200---"
        #self.printBlock(data)
        
        
        aCode = self.getBlockIntEntry(data,1)
        print "aCode = ",aCode
        (self.tableCode,three,self.iSubcase) = self.parseApproachCode(data)
        #iSubcase = self.getValues(data,'i',4)

        self.rCode  = self.getValues(data,'i',8) ## random code
        self.fCode  = self.getValues(data,'i',9) ## format code
        self.numwde = self.getValues(data,'i',10) ## number of words per entry in record; @note is this needed for this table ???
        self.acousticFlag = self.getValues(data,'f',13) ## acoustic pressure flag
        self.thermal      = self.getValues(data,'i',23) ## thermal flag; 1 for heat ransfer, 0 otherwise
        
        ## assuming tCode=1
        if self.approachCode==1:   # statics / displacement / heat flux
            self.lsdvmn = self.getValues(data,'i',5) ## load set number
        elif self.approachCode==2: # real eigenvalues
            self.mode      = self.getValues(data,'i',5) ## mode number
            self.eigr      = self.getValues(data,'f',6) ## real eigenvalue
            self.modeCycle = self.getValues(data,'f',7) ## mode or cycle @todo confused on the type ???
        elif self.approachCode==3: # differential stiffness
            self.lsdvmn = self.getValues(data,'i',5) ## load set number
        elif self.approachCode==4: # differential stiffness
            self.lsdvmn = self.getValues(data,'i',5) ## load set number
        elif self.approachCode==5:   # frequency
            self.freq = self.getValues(data,'f',5) ## frequency

        elif self.approachCode==6: # transient
            self.dt = self.getValues(data,'f',5) ## time step
            print "*****self.dt = ",self.dt
        elif self.approachCode==7: # pre-buckling
            self.lsdvmn = self.getValues(data,'i',5) ## load set
        elif self.approachCode==8: # post-buckling
            self.lsdvmn = self.getValues(data,'i',5) ## mode number
            self.eigr   = self.getValues(data,'f',6) ## real eigenvalue
        elif self.approachCode==9: # complex eigenvalues
            self.mode   = self.getValues(data,'i',5) ## mode
            self.eigr   = self.getValues(data,'f',6) ## real eigenvalue
            self.eigi   = self.getValues(data,'f',7) ## imaginary eigenvalue
        elif self.approachCode==10: # nonlinear statics
            self.lftsfq = self.getValues(data,'i',5) ## load step
        elif self.approachCode==11: # old geometric nonlinear statics
            self.lsdvmn = self.getValues(data,'i',5)
        elif self.approachCode==12: # contran ? (may appear as aCode=6)  --> straight from DMAP...grrr...
            self.lsdvmn = self.getValues(data,'i',5)
        else:
            raise RuntimeError('invalid approach code...approachCode=%s' %(self.approachCode))
        # tCode=2
        #if self.analysisCode==2: # sort2
        #    self.lsdvmn = self.getValues(data,'i',5)
        
        print "*iSubcase=%s"%(self.iSubcase)
        print "approachCode=%s tableCode=%s thermal=%s" %(self.approachCode,self.tableCode,self.thermal)
        print self.codeInformation(sCode=None,tCode=None,thermal=self.thermal)

        #self.printBlock(data)

        word = self.readString(384) # titleSubtitleLabel
        #print "word = |%s|" %(word)
        #word = self.readString(4*(63+33)) # title, subtitle, and label
        Title    = word[0:128]
        Subtitle = word[128:256]
        Label    = word[256:]
        #print "Title    %s |%s|" %(len(Title   ),Title)
        #print "Subtitle %s |%s|" %(len(Subtitle),Subtitle)
        #print "Label    %s |%s|" %(len(Label   ),Label)
        print "Title    %s |%s|" %(len(Title   ),Title.strip())
        print "Subtitle %s |%s|" %(len(Subtitle),Subtitle.strip())
        print "Label    %s |%s|" %(len(Label   ),Label.strip())


        self.readHollerith()
        #return (analysisCode,tableCode,thermal)

        #if self.j==3:
        #    #print str(self.obj)
        #    sys.exit('checkA...j=%s dt=6E-2 dx=%s dtActual=%f' %(self.j,'1.377e+01',self.dt))
        ###

    def readTable_OUGV1_4(self,iTable):
        #self.readMarkers([iTable,1,0])
        markerA = 4
        
        j = 0
        while markerA>None:
            self.markerStart = copy.deepcopy(self.n)
            #self.printSection(180)
            self.readMarkers([iTable,1,0])
            print "starting OUGV1 table 4..."
            isTable4Done,isBlockDone = self.readTable_OUGV1_4_Data(iTable)
            if isTable4Done:
                print "done with OUGV14"
                self.n = self.markerStart
                self.op2.seek(self.n)
                break
            print "finished reading ougv1 table..."
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

    def readTable_OUGV1_4_Data(self,iTable): # iTable=-4
        isTable4Done = False
        isBlockDone = False

        bufferWords = self.getMarker('OUGV1')
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
            if self.approachCode==1 and self.tableCode==1: # displacement
                print "isDisplacement"
                self.obj = displacementObject(self.iSubcase)
                self.displacements[self.iSubcase] = self.obj
            elif self.approachCode==1 and self.tableCode==3: # spc forces
                print "isForces"
            elif self.approachCode==6 and self.tableCode==1: # transient displacement
                print "isTransientDisplacement"
                if self.iSubcase in self.displacements:
                    self.obj = self.displacements[self.iSubcase]
                    self.obj.updateDt(self.dt)
                else:
                    self.obj = displacementObject(self.iSubcase,dt=self.dt)
                self.displacements[self.iSubcase] = self.obj
            else:
                raise Exception('not supported OUGV1 solution...')
            ###
        elif self.thermal==1:
            if self.approachCode==1 and self.tableCode==1: # temperature
                print "isTemperature"
                self.temperatures[self.iSubcase] = temperatureObject(self.iSubcase)

            elif self.approachCode==1 and self.tableCode==3: # heat fluxes
                print "isFluxes"
                self.obj = temperatureObject(self.iSubcase,dt=self.dt)
                self.fluxes[self.iSubcase] = self.obj
            elif self.approachCode==6 and self.tableCode==1: # transient temperature
                print "isTransientTemperature"
                if self.iSubcase in self.temperatures:
                    self.obj = self.temperatures[self.iSubcase]
                    self.obj.updateDt(self.dt)
                else:
                    self.obj = temperatureObject(self.iSubcase,dt=self.dt)
                self.temperatures[self.iSubcase] = self.obj
            else:
                raise Exception('not supported OUGV1 solution...')
            ###
        else:
            raise Exception('invalid thermal flag...not 0 or 1...flag=%s' %(self.thermal))
        ###
        self.readScalars(data,self.obj)
        #print self.obj
        
        print "-------finished OUGV1----------"
        return (isTable4Done,isBlockDone)

    def isDisplacement(self):
        if self.approachCode==1 and self.thermal==0:
            return True
        return False

    def isTransientDisplacement(self):
        if self.approachCode==6 and self.tableCode==1 and self.thermal==0:
            return True
        return False

    def isTemperature(self):
        if self.approachCode==1 and self.thermal==1:
            return True
        return False

    def isTransientTemperature(self):
        if self.approachCode==6 and self.tableCode==1 and self.thermal==1:
            return True
        return False

    def isForces(self,tableCode,approachCode,thermal):
        if(approachCode==1 and tableCode==3 and thermal==0):
            return True
        return False

    def isFluxes(self,tableCode,approachCode,thermal):
        if(approachCode==1 and tableCode==3 and thermal==1):
            return True
        return False

    def readTable_OEF1X(self):
        ## OEF1X
        tableName = self.readTableName(rewind=False) # OEF1X
        print "tableName = |%r|" %(tableName)

        self.readMarkers([-1,7])
        ints = self.readIntBlock()
        print "*ints = ",ints

