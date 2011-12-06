import os
import sys
import struct
from struct import unpack

#from pyNastran.op2.op2Errors import *
from pyNastran.bdf.cards.nodes import GRID
from pyNastran.bdf.cards.coordinateSystems import CORD1R,CORD2R,CORD2C,CORD3G #CORD1C,CORD1S,CORD2S


class Geometry1(object):

    def readDesvar(self,data):
        out = unpack('iiccccccccfff',data[0:28])
        (idvid,dvid,alabel1,bLabel1,cLabel1,dLabel1,alabel2,bLabel2,cLabel2,dLabel2,vmin,vmax,delx) = out
        label1 = alabel1+bLabel1+cLabel1+dLabel1
        label2 = alabel2+bLabel2+cLabel2+dLabel2
        #print "ivid=%s dvid=%s label1=%s label2=%s vmax=%g vmin=%g delx=%g" %(idvid,dvid,label1,label2,vmin,vmax,delx)
        #self.readData(8)

    def readTable_DesTab(self):
        tableName = self.readTableName(rewind=False) # DESTAB
        print "tableName = |%r|" %(tableName)

        self.readMarkers([-1,7],'DESTAB')
        ints = self.readIntBlock()
        #print "*ints = ",ints

        self.readMarkers([-2,1,0],'DESTAB')
        bufferWords = self.getMarker()
        print "bufferWords = ",bufferWords
        word = self.readStringBlock() # DESTAB
        #print "word = ",word

        iTable = -3
        imax   = -241
        while iTable > imax:
            #self.printSection(80)
            self.readMarkers([iTable,1,0],'DESTAB')
            bufferWords = self.getMarker()
            #print "bufferWords = ",bufferWords
            data = self.readBlock()
            #print "len(data) = ",len(data)
            self.readDesvar(data)
            iTable -= 1


        #self.op2Debug.write('bufferWords=%s\n' %(str(bufferWords)))
        #print "1-bufferWords = ",bufferWords,bufferWords*4

        print self.printSection(300)
        sys.exit('asdf')

    def readTable_Geom1(self):
        self.iTableMap = {
                            #(1701,17,6):    self.readCord1C, # record 1
                            (1801,18,5):     self.readCord1R, # record 2
                            #(1901,19,7):    self.readCord1S, # record 3
                            (2001,20,9):     self.readCord2C, # record 4
                            (2101,21,8):     self.readCord2R, # record 5
                            #(2201,22,10):   self.readCord2S, # record 6
                            (14301,143,651): self.readCord3G, # record 7
                            (4501,45,1):     self.readGrid,   # record 17 - slow
                            (5301,53,4):     self.readSEQGP,  # record 27 - not done

                            (2201,22, 10):    self.readFake,
                            (6101,61,388):    self.readFake,
                         }
        self.readRecordTable('GEOM1')

    def readSEQGP(self,data):
        pass

    def readCord1C(self,data):
        """
        (1701,17,6) - the marker for Record 1
        """
        print "reading CORD1C"
        while len(data)>=24: # 6*4
            eData = data[:24]
            data  = data[24:]
            (cid,one,two,g1,g2,g3) = unpack('iiiiii',eData)
            dataIn = [cid,g1,g2,g3]
            coord = CORD1R(None,dataIn)
            self.addCoord(coord)
        ###

    def readCord1R(self,data):
        """
        (1801,18,5) - the marker for Record 2
        """
        print "reading CORD1R"
        while len(data)>=24: # 6*4
            eData = data[:24]
            data  = data[24:]
            (cid,one,one,g1,g2,g3) = unpack('iiiiii',eData)
            dataIn = [cid,g1,g2,g3]
            coord = CORD1R(None,dataIn)
            self.addCoord(coord)
        ###

    def readCord1S(self,data):
        """
        (1901,19,7) - the marker for Record 3
        """
        print "reading CORD1S"
        while len(data)>=24: # 6*4
            eData = data[:24]
            data  = data[24:]
            (cid,three,one,g1,g2,g3) = unpack('iiiiii',eData)
            dataIn = [cid,g1,g2,g3]
            coord = CORD1S(None,dataIn)
            self.addCoord(coord,allowOverwrites=True)
        ###

    def readCord2C(self,data):
        """
        (2001,20,9) - the marker for Record 4
        """
        print "reading CORD2C"
        while len(data)>=52: # 13*4
            eData = data[:52]
            data  = data[52:]
            (cid,two,two,rid,a1,a2,a3,b1,b2,b3,c1,c2,c3) = unpack('iiiifffffffff',eData)
            #print "cid=%s two=%s two=%s rid=%s a1=%s a2=%s a3=%s b1=%s b2=%s b3=%s c1=%s c2=%s c3=%s" %(cid,two,two,rid,a1,a2,a3,b1,b2,b3,c1,c2,c3)
            dataIn = [cid,rid,a1,a2,a3,b1,b2,b3,c1,c2,c3]
            coord = CORD2C(None,dataIn)
            self.addCoord(coord,allowOverwrites=True)
        ###

    def readCord2R(self,data):
        """
        (2101,21,8) - the marker for Record 5
        """
        print "reading CORD2R"
        while len(data)>=52: # 13*4
            eData = data[:52]
            data  = data[52:]
            (cid,one,two,rid,a1,a2,a3,b1,b2,b3,c1,c2,c3) = unpack('iiiifffffffff',eData)
            dataIn = [cid,rid,a1,a2,a3,b1,b2,b3,c1,c2,c3]
            #print "cid=%s one=%s two=%s rid=%s a1=%s a2=%s a3=%s b1=%s b2=%s b3=%s c1=%s c2=%s c3=%s" %(cid,one,two,rid,a1,a2,a3,b1,b2,b3,c1,c2,c3)
            coord = CORD2R(None,dataIn)
            self.addCoord(coord,allowOverwrites=True)
        ###

    def readCord2S(self,data):
        """
        (2201,22,10) - the marker for Record 6
        """
        print "reading CORD2S"
        while len(data)>=52: # 13*4
            eData = data[:52]
            data  = data[52:]
            (cid,sixty5,eight,rid,a1,a2,a3,b1,b2,b3,c1,c2,c3) = unpack('iiiifffffffff',eData)
            #print "cid=%s sixty5=%s eight=%s rid=%s a1=%s a2=%s a3=%s b1=%s b2=%s b3=%s c1=%s c2=%s c3=%s" %(cid,sixty5,eight,rid,a1,a2,a3,b1,b2,b3,c1,c2,c3)
            dataIn = [cid,rid,a1,a2,a3,b1,b2,b3,c1,c2,c3]
            coord = CORD2S(dataIn)
            self.addCoord(coord,allowOverwrites=True)
        ###

    def readCord3G(self,data):
        """
        (14301,143,651) - the marker for Record 7
        @todo isnt this a CORD3G, not a CORD3R ???
        """
        print "reading CORD3G"
        while len(data)>=16: # 4*4
            eData = data[:16]
            data  = data[16:]
            (cid,n1,n2,n3) = unpack('iiii',eData)
            dataIn = [cid,n1,n2,n3]
            coord = CORD3G(None,dataIn)
            self.addCoord(coord,allowOverwrites=True)
        ###

    def readGrid(self,data):  # 21.8 sec, 18.9
        """(4501,45,1) - the marker for Record 17"""
        print "reading GRID"
        #return
        
        n=0
        nEntries = len(data)/32
        for i in range(nEntries):
            eData = data[n:n+32]
            out = unpack('iifffiii',eData)

            (nID,cp,x1,x2,x3,cd,ps,seid) = out
            if cd>=0 and nID<10000000:
                node = GRID(None,out)
                self.addNode(node)
            else:
                print "*nID=%s cp=%s x1=%-5.2f x2=%-5.2f x3=%-5.2f cd=%-2s ps=%s seid=%s" %(nID,cp,x1,x2,x3,cd,ps,seid)
            ###
            #print str(grid)[:-1]
            n+=32
        ###
        data = data[n:]
        #assert len(data)==0,'len(data)!=0   len(data)=%s' %(len(data))
        #print "len(data) = ",len(data)

# integrated into readGeom1 function
#-------------------------------------
# not integrated

