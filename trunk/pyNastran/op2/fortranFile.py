#import os
#import sys

import struct
from struct import unpack,pack

class FortranFile(object):
    def __init__(self):
        self.endian = '>'
    
    def setEndian(self,endian):
        self.endian = endian

    def readHollerith(self):
        self.skip(4)  # weird hollerith

    def readHeader(self):
        data = self.op2.read(12)
        ints = self.getInts(data)
        print "header ints = %s" %(repr(ints))
        self.n += 12
        assert ints[0]==ints[2]==4,"header ints = (%s, %2s, %s)" %(ints)
        return ints[1]

    def readString(self,nData):
        data = self.op2.read(nData)
        string = ''.join(self.getStrings(data))
        self.n += nData
        return string

    def readString(self,nData):
        data = self.op2.read(nData)
        self.n += nData
        words = self.getStrings(data)
        return ''.join(words)

    def readInts(self,nInts):
        nData = 4*nInts
        #print "nData = ",nData
        data = self.op2.read(nData)

        iFormat = 'i'*nInts
        ints = unpack(iFormat,data)
        self.n+=nData
        return ints

    def readDoubles(self,nData):
        data = self.op2.read(nData)
        self.n += nData
        doubles = self.getDoubles(data)
        return doubles

    def readFloats(self,nData):
        data = self.op2.read(nData)
        self.n += nData
        floats = self.getFloats(data)
        return floats

    def getStrings(self,data):
        n = len(data)
        sFormat = 's'*n
        strings = unpack(sFormat,data)
        return strings

    def getInts(self,data):
        n = len(data)
        nInts = n/4
        #print "nInts = ",nInts
        iFormat = 'i'*nInts
        ints = unpack(iFormat,data[:nInts*4])
        return ints

    def getLongs(self,data):
        n = len(data)
        nLongs = n/4
        #print "nLongs = ",nLongs
        #a = pack('l',200)
        #print "len(a) = ",len(a)
        
        LFormat = 'l'*nLongs
        longs = unpack(LFormat,data[:nLongs*4])
        return longs

    def getFloats(self,data):
        n = len(data)
        nFloats = n/4
        fFormat = 'f'*nFloats
        ints = unpack(fFormat,data[:nFloats*4])
        return ints

    def getDoubles(self,data):
        n = len(data)
        nDoubles = n/8
        dFormat = 'd'*nDoubles
        ints = unpack(dFormat,data[:nDoubles*8])
        return ints
    
    def printBlock(self,data):
        ints    = self.getInts(data)
        #longs   = self.getLongs(data)
        floats  = self.getFloats(data)
        doubles = self.getDoubles(data)
        strings = self.getStrings(data)
        print "ints    = ",ints
        #print "longs   = ",longs
        print "floats  = ",floats
        print "doubles = ",doubles
        print "strings = |%r|" %(''.join(strings))

    def getData(self,n):
        data = self.op2.read(n)
        self.n+=n
        return data

    def getBlockIntEntry(self,data,n):
        data2 = data[4*(n-1):4*(n-1)+4]
        return struct.unpack('i',data2)[0]
        
    def printSection(self,nBytes):
        """
        prints data, but doesnt move the cursor
        """
        data = self.op2.read(nBytes)
        self.printBlock(data)
        self.op2.seek(self.n)
    
    def skip(self,n):
        #print "\n--SKIP--"
        #print "tell = ",self.op2.tell()
        #print "n = ",n
        #print "self.n = ",self.n
        self.n += n
        self.op2.seek(self.n)
        #print "*tell = ",self.op2.tell()
        #print "\n"

    def scan(self,n):
        data = self.op2.read(n)
        self.n+=n

    def getTableCode(self):
        tableCode = self.readHeader()
        return tableCode

    def getMarker(self):
        tableCode = self.readHeader()
        return tableCode
        
    def readMarkers(self,markers):
        for marker in markers:
            tableCode = self.readHeader()
            assert marker==tableCode
        ###

    def getNMarkers(self,nMarkers,rewind=False):
        markers = []
        for iMarker in range(nMarkers):
            tableCode = self.readHeader()
            markers.append(tableCode)
        ###
        if rewind:
            self.n -= 12*nMarkers
            self.op2.seek(self.n)
            
        return markers

    def isTableDone(self,expectedMarkers):
        markers = self.getNMarkers(len(expectedMarkers),rewind=True)
        print "getMarkers = ",markers

        if markers==[-1,7]:
            return True
        elif markers==expectedMarkers:
            return False
        else:
            raise RuntimeError('this should never happen...invalid markers...expected=%s markers=%s' %(expectedMarkers,markers))

    def goto(self,n):
        self.op2.seek(n)

    def readStringBlock(self):
        data = self.readBlock()
        nLetters = len(data)
        letters = unpack('s'*nLetters,data)
        word = ''.join(letters)
        #print "word = |%s|" %(word)
        return word

    def readBlock(self):
        data = self.op2.read(4)
        nValues, = unpack('i',data)
        self.n+=4
        data = self.op2.read(nValues)
        self.n+=nValues+4
        self.goto(self.n)
        return data

    def readIntBlock(self):
        data = self.readBlock()
        nInts = len(data)/4
        print "**nInts = ",nInts
        ints = unpack('i'*nInts,data)
        return ints

    def readFloatBlock(self):
        data = self.readBlock()
        nFloats = len(data)/4
        floats = unpack('f'*nFloats,data)
        return floats

    def readFloatBlock(self):
        data = self.readBlock()
        nDoubles = len(data)/8
        doubles = unpack('d'*nDoubles,data)
        return doubles

    def rewind(self,n):
        """doesnt support a full rewind, only a partial"""
        self.n -= n
        self.op2.seek(self.n)

    def readTableName(self,rewind=True):
        n = self.n
        #print ""
        self.readMarkers([0,2])
        word = self.readStringBlock()  # GEOM1
        #print "*word = |%r|" %(word)

        #print "n      = ",n
        #print "self.n = ",self.n
        #print "op2.tell = ",self.op2.tell()
        #print "******"
        if rewind:
            self.n = n
            self.op2.seek(n)
        #print "n      = ",n
        #print "self.n = ",self.n
        #print "op2.tell = ",self.op2.tell()
        return word.strip()

    def skipNextTable(self,bufferSize=10000):
        word = self.readTableName(rewind=False) # GEOM1
        print "skippingTable |%s|" %(word)

        self.readMarkers([-1,7])
        
        dataPack = (4,1,4,  4,0,4,  4,0,4)  # marks the end of the table
        binaryData = pack('iiiiiiiii',*dataPack)
        #[1,0,0] marks the end of the table
        
        i = 0
        error = 50
        n = self.n
        endIndex = -1
        data = "dummy"
        while endIndex== -1 and len(data)>0:
            data     = self.op2.read(bufferSize+error)
            endIndex = data.find(binaryData)
            self.op2.seek(n+i*bufferSize)
            i+=1
        
        assert endIndex>0,'couldnt find the end of the table'
        self.n = self.n+endIndex+36 # 36 so it gets to the end of the table markersNext=[0] or [2]
        self.op2.seek(self.n)
        "self.op2.tell() = ",self.op2.tell()

        marker = self.getMarker()
        if marker==2:
            isAnotherTable = True
        else:# marker=0
            isAnotherTable = False
        ###
        #print "isAnotherTable = ",isAnotherTable
        self.n -= 24  # subtract off the header [0,2] or [0,0]
        self.op2.seek(self.n)
        return isAnotherTable

    def hasMoreTables(self):
        print self.printSection(120)
        try:
            marker1 = self.getMarker()
            marker2 = self.getMarker()
            marker = [marker1,marker2]
            print "marker = ",marker
            if marker==[0,2]:
                isAnotherTable = True
            else:# marker=0
                isAnotherTable = False
            ###
            #print "isAnotherTable = ",isAnotherTable
            self.n -= 24  # subtract off the header [0,2] or [0,0]
            self.op2.seek(self.n)
        except IndexError:
            isAnotherTable = False
        return isAnotherTable
    