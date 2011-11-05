#import os
#import sys

import struct
from struct import unpack

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
        print "header ints = ",ints
        self.n += 12
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
        iFormat = 'i'*nInts
        ints = unpack(iFormat,data[:nInts*4])
        return ints

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
        floats  = self.getFloats(data)
        doubles = self.getDoubles(data)
        strings = self.getStrings(data)
        print "ints    = ",ints
        print "floats  = ",floats
        print "doubles = ",doubles
        print "strings = |%r|" %(''.join(strings))

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

    def startTable(self,markers):
        self.readMarkers(markers)
        word = self.readStringBlock()
        return word

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

    def goto(self,n):
        self.op2.seek(n)

    def readStringBlock(self):
        data = self.readBlock()
        nLetters = len(data)
        letters = unpack('s'*nLetters,data)
        word = ''.join(letters)
        print "word = |%s|" %(word)
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
