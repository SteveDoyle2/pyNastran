import sys
from struct import unpack
from op2_Objects import temperatureObject,fluxObject,displacementObject

class OUGV1(object):

    def getValues(self,data,sFormat,iWordStart,iWordStop=None):
        if iWordStop==None:
            ds = data[iWordStart*4:(iWordStart+1)*4]
            return unpack(sFormat,ds)[0]
            
        #print "type(data) = ",type(data)
        ds = data[iWordStart*4:iWordStop*4]
        return unpack(sFormat,ds)
        
    def readTable_OUGV1_3(self,iTable): # iTable=-3
        self.readMarkers([iTable,1,0],'OUGV1')
        bufferWords = self.getMarker()
        print "bufferWords = ",bufferWords,bufferWords*4,'\n'

        data = self.getData(4)
        bufferSize, = unpack('i',data)
        print "bufferSize = ",bufferSize
        data = self.getData(4*50)
        aCode = self.getBlockIntEntry(data,1)
        print "aCode = ",aCode
        (self.tableCode,three,self.iSubcase) = self.parseApproachCode(data)
        #iSubcase = self.getValues(data,'i',4)
        dt       = self.getValues(data,'f',5)
        self.thermal  = self.getValues(data,'i',23)
        print "*iSubcase=%s dt=%s"%(self.iSubcase,dt)
        print "approachCode=%s tableCode=%s thermal=%s" %(self.approachCode,self.tableCode,self.thermal)
        print self.codeInformation(sCode=None,tCode=None,thermal=self.thermal)

        #self.printBlock(data)

        word = self.readString(384)
        print "word = |%s|" %(word)
        self.readHollerith()
        #return (analysisCode,tableCode,thermal)

    def isDisplacement(self):
        if self.approachCode==1:
            return True
        return False

    def isTransientDisplacement(self):
        if self.approachCode==6 and self.tableCode==1:
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

    def readTable_OUGV1_4(self,iTable): # iTable=-4
        self.readMarkers([iTable,1,0]) # iTable=-4
        bufferWords = self.getMarker('OUGV1')
        data = self.readBlock()
        #self.printBlock(data)
        
        if self.isDisplacement():
            self.dispObj = displacementObject(self.iSubcase)
        elif self.isTransientDisplacement():
            self.dispObj = displacementObject(self.iSubcase)
        #elif self.isForces():
            
        else:
            raise Exception('not supported OUGV1 solution...')
        self.readScalars(data,self.dispObj)

    def readTable_OUGV1(self):
        ## OUGV1
        word = self.readTableName(rewind=False) # OUGV1
        print "word = |%r|" %(word)

        self.readMarkers([-1,7],'OUGV1')
        ints = self.readIntBlock()
        print "*ints = ",ints

        self.readMarkers([-2,1,0],'OUGV1') # 7
        bufferWords = self.getMarker()
        print "bufferWords = ",bufferWords,bufferWords*4
        ints = self.readIntBlock()
        print "*ints = ",ints
        
        iTable=-3

        markerA = -4
        markerB = 0

        while [markerA,markerB]!=[0,2]:
            self.readTable_OUGV1_3(iTable)
            self.readTable_OUGV1_4(iTable-1)
            iTable -= 2
        
            print str(self.dispObj)

            n = self.n
            self.readMarkers([iTable,1,0],'OUGV1')

            markerA = self.getMarker('A')
            markerB = self.getMarker('B')
            self.n-=24
            self.op2.seek(self.n)
            print "markerA=%s markerB=%s" %(markerA,markerB)
            

            self.printSection(120)
            #break

        #sys.exit('end of displacementA')
        return



        #self.readMarkers([-5,1,0])
        #sys.exit('end of displacement')

        bufferWords = self.getMarker()
        data = self.readBlock()

        sys.exit('end of displacementB')
        self.readMarkers([-6,1,0])
        bufferWords = self.getMarker()
        #self.printSection(100)
        data = self.readBlock()
        #sys.exit('ougv1 - stop')
        self.readMarkers([-7,1,0])
        bufferWords = self.getMarker()
        data = self.readBlock()

        self.readMarkers([-8,1,0])
        bufferWords = self.getMarker()
        data = self.readBlock()

        self.readMarkers([-9,1,0])
        bufferWords = self.getMarker()
        data = self.readBlock()

        self.readMarkers([-10,1,0])
        bufferWords = self.getMarker()
        data = self.readBlock()

        self.readMarkers([-11,1,0])
        bufferWords = self.getMarker()
        data = self.readBlock()

        self.readMarkers([-12,1,0])
        bufferWords = self.getMarker()
        data = self.readBlock()

        self.readMarkers([-13,1,0])
        bufferWords = self.getMarker()
        data = self.readBlock()

        self.readMarkers([-14,1,0])
        bufferWords = self.getMarker()
        data = self.readBlock()

        self.readMarkers([-15,1,0])
        bufferWords = self.getMarker()
        data = self.readBlock()

        self.readMarkers([-16,1,0])
        bufferWords = self.getMarker()
        data = self.readBlock()

        self.readMarkers([-17,1,0])
        bufferWords = self.getMarker()
        bufferWords = self.getMarker()
        data = self.readBlock()



        self.readMarkers([-1,7])
        ints = self.readIntBlock()

        self.readMarkers([-2,1,0])
        bufferWords = self.getMarker()
        data = self.readBlock()

        self.readMarkers([-3,1,0])
        bufferWords = self.getMarker()
        data = self.readBlock()

        self.readMarkers([-4,1,0])
        bufferWords = self.getMarker()
        data = self.readBlock()


        self.readMarkers([-5,1,0])
        bufferWords = self.getMarker()
        data = self.readBlock()


        self.readMarkers([-6,1,0])
        bufferWords = self.getMarker()
        data = self.readBlock()


        self.readMarkers([-7,1,0])
        bufferWords = self.getMarker()
        data = self.readBlock()


        self.readMarkers([-8,1,0])
        bufferWords = self.getMarker()
        data = self.readBlock()

        self.readMarkers([-9,1,0])
        bufferWords = self.getMarker()
        data = self.readBlock()

        self.readMarkers([-10,1,0])
        bufferWords = self.getMarker()
        data = self.readBlock()

        self.readMarkers([-11,1,0])
        bufferWords = self.getMarker()
        data = self.readBlock()


        self.readMarkers([-12,1,0])
        bufferWords = self.getMarker()
        data = self.readBlock()

        self.readMarkers([-13,1,0])
        bufferWords = self.getMarker()
        data = self.readBlock()


        self.readMarkers([-14,1,0])
        bufferWords = self.getMarker()
        data = self.readBlock()

        self.readMarkers([-15,1,0])
        bufferWords = self.getMarker()
        data = self.readBlock()

        self.readMarkers([-16,1,0])
        bufferWords = self.getMarker()
        data = self.readBlock()

        self.readMarkers([-17,1,0])
        bufferWords = self.getMarker()
        data = self.readBlock()

        self.readMarkers([-18,1,0])
        bufferWords = self.getMarker()
        eid = unpack('i',data[0:4])

        data = self.readBlock()
        print "---------"
        self.printBlock(data)
        
        
        

        #self.readMarkers([-1,7])
        #bufferWords = self.getMarker()
        #data = self.readBlock()
        #self.printSection(160)


        #assert self.op2.tell()==4780,self.op2.tell()
        #sys.exit('end of displacements')
        #self.readTable_OEF1X()
        sys.exit('end of ougv1')

    def readTable_OEF1X(self):
        ## OEF1X
        word = self.readTableName(rewind=False) # OEF1X
        print "word = |%r|" %(word)

        self.readMarkers([-1,7])
        ints = self.readIntBlock()
        print "*ints = ",ints

