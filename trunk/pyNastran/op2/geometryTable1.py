import os
import sys
import struct
from struct import unpack

#from pyNastran.op2.op2Errors import *
from pyNastran.bdf.cards.nodes import GRID

class GeomObj(object):
    def __init__(self):
        pass
    def geomFunc(self,data):
        pass

class Geometry1(object):

    def readTable_Geom1(self):
        self.nodes  = {}
        self.coords = {}

        tableName = self.readTableName(rewind=False) # GEOM1
        self.tableInit(tableName)
        print "*tableName = |%r|" %(tableName)

        self.readMarkers([-1,7])
        fields = self.readIntBlock()
        #print "fields = ",fields

        self.readMarkers([-2,1,0]) # 2
        bufferWords = self.getMarker()
        #print "bufferWords = ",bufferWords,bufferWords*4
        word = self.readStringBlock()

        self.iTableMap = {
                            (1701,17,6):     self.readCord1C, # record 1
                            (1801,18,5):     self.readCord1R, # record 2
                            (1901,19,7):     self.readCord1S, # record 3
                            (2001,20,9):     self.readCord2C, # record 4
                            (2101,21,8):     self.readCord2R, # record 5
                            (2201,22,10):    self.readCord2S, # record 6
                            (14301,143,651): self.readCord3R, # record 7
                            (4501,45,1):     self.readNodes,  # record 17
                         }
        iTable = -3
        while 1:  ## @todo could this cause an infinite loop...i dont this so...
            (tableName,isNextTable,isNextSubTable) = self.readGeomSubTable(iTable)
        
            if self.checkForNextTable():
                #sys.exit('end of geom1')
                return
            iTable -= 1
        ###
        sys.exit('end of geom1-huh')

    def readCord1C(self,data):
        """
        (1701,17,6) - the marker for Record 1
        """
        print "reading CORD1C"
        while len(data)>=24: # 6*4
            eData = data[:24]
            data  = data[24:]
            (cid,one,two,g1,g2,g3) = unpack('iiiiii',eData)
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
            #print "cid=%s two=%s two=%s rid=%s a1=%s a2=%s a3=%s b1=%s b2=%s b3=%s c1=%s c2=%s c3=%s" %(cid,one,two,rid,a1,a2,a3,b1,b2,b3,c1,c2,c3)
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
            #print "cid=%s one=%s two=%s rid=%s a1=%s a2=%s a3=%s b1=%s b2=%s b3=%s c1=%s c2=%s c3=%s" %(cid,one,two,rid,a1,a2,a3,b1,b2,b3,c1,c2,c3)
        ###

    def readCord2S(self,data):
        """(2201,22,10) - the marker for Record 6"""
        print "reading CORD2S"
        while len(data)>=52: # 13*4
            eData = data[:52]
            data  = data[52:]
            (cid,sixty5,eight,rid,a1,a2,a3,b1,b2,b3,c1,c2,c3) = unpack('iiiifffffffff',eData)
            #print "cid=%s sixty5=%s eight=%s rid=%s a1=%s a2=%s a3=%s b1=%s b2=%s b3=%s c1=%s c2=%s c3=%s" %(cid,sixty5,eight,rid,a1,a2,a3,b1,b2,b3,c1,c2,c3)
        ###

    def readCord3R(self,data):
        """
        (14301,143,651) - the marker for Record 7
        @todo isnt this a CORD3G, not a CORD3R ???
        """
        print "reading CORD3R"
        while len(data)>=16: # 4*4
            eData = data[:16]
            data  = data[16:]
            (cid,n1,n2,n3) = unpack('iiii',eData)
        ###

    def readNodes(self,data):
        """(4501,45,1) - the marker for Record 17"""
        print "reading NODES"
        while len(data)>=32: # 8*4
            eData = data[:32]
            data  = data[32:]
            out = unpack('iifffiii',eData)

            (nID,cp,x1,x2,x3,cd,ps,seid) = out
            #print "nID=%s cp=%s x1=%s x2=%s x3=%s cd=%s ps=%s seid=%s" %(nID,cp,x1,x2,x3,cd,ps,seid)
            node = GRID(None,out)
            self.addNode(node)
            
            #print str(grid)[:-1]
        ###
        #print "len(data) = ",len(data)
        

# GEOM1
#----
    def readGeomSubTable(self,iTable):
        i=0
        isNextTable=False
        isNextSubTable=False
        self.readMarkers([iTable,1,0])
        #print self.iTableMap

        tableName = self.readTableName(rewind=True,stopOnFailure=False)
        if tableName:
            print "**tableName = |%r|" %(tableName)
            return tableName,isNextTable,isNextSubTable

        data = ''
        isTableActive=False
        while isNextSubTable==False and isNextTable==False:
            #print self.printSection(100)
            marker = self.getMarker()
            #print "marker = ",marker
            if marker<0:
                msg = 'marker is less than 0...'
                raise Exception(msg)
            data += self.readBlock()
            if not isTableActive:
                tableType = unpack('iii',data[:12])
                data = data[12:]

            #print "iTable=%s lenGeomData=%s" %(iTable,len(data))

            if tableType in self.iTableMap:
                #print "reading  iTable=%-3s with tableType=%s" %(iTable,tableType)
                self.iTableMap[tableType](data)
            else:
                print "skipping iTable=%-3s with tableType=%s" %(iTable,tableType)

            #self.op2Debug.write('ints = %s\n' %(str(ints)))

            isNextTable = self.checkForNextTable()
            isNextSubTable = self.checkForNextSubTable(iTable-1)
            #print "i=%s tell=%s isNextTable=%s isNextSubTable=%s" %(i,self.op2.tell(),isNextTable,isNextSubTable)
            #if i==13:
            #    sys.exit('stopA')
            i+=1
            isTableActive=True
        ### while
        
        #print "exiting the geom sub table"
        return (tableName,isNextTable,isNextSubTable)

