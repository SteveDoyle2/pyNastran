#import os
#import sys

import struct
from struct import unpack

class FortranFile(object):
    def __init__(self):
        self.endian = '>'
    
    def readHeader(self):
        data = self.op2.read(12)
        ints = self.getInts(data)
        print "header ints = ",ints
        self.n += 12
        return ints[1]

    def getStrings(self,data):
        n = len(data)
        sFormat = 's'*n
        strings = unpack(sFormat,data)
        return strings

    def readInts(self,nInts):
        nData = 4*nInts
        #print "nData = ",nData
        data = self.op2.read(nData)

        iFormat = 'i'*nInts
        ints = unpack(iFormat,data)
        self.n+=nData
        return ints

    def getInts(self,data):
        n = len(data)
        nInts = (n/4)
        iFormat = 'i'*nInts
        ints = unpack(iFormat,data[:nInts*4])
        return ints

    def getFloats(self,data):
        n = len(data)
        nFloats = (n/4)
        iFormat = 'f'*nFloats
        ints = unpack(iFormat,data[:nFloats*4])
        return ints

    def printSection(self,nBytes):
        """
        prints data, but doesnt move the cursor
        """
        data = self.op2.read(nBytes)
        ints    = self.getInts(data)
        floats  = self.getFloats(data)
        strings = self.getStrings(data)
        print "ints    = ",ints
        print "floats  = ",floats
        print "strings = |%r|" %(''.join(strings))
        self.op2.seek(self.n)
    
    def readString(self,nData):
        data = self.op2.read(nData)
        string = ''.join(self.getStrings(data))
        self.n += nData
        return string

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
        
    def readMarkers(self,markers):
        for marker in markers:
            tableCode = self.readHeader()
            assert marker==tableCode
        ###

    def readFloats(self,nData):
        data = self.op2.read(nData)
        self.n += nData
        floats = self.getFloats(data)
        return floats

    def readString(self,nData):
        data = self.op2.read(nData)
        self.n += nData
        words = self.getStrings(data)
        return ''.join(words)

    def readStringBlock(self):
        data = self.op2.read(4)
        nValues, = unpack('i',data)
        
        self.n+=4
        word = self.readString(nValues)
        print "word = |%s|" %(word)
        self.skip(4)
        return word

    def readIntBlock(self):
        data = self.op2.read(4)
        self.n+=4

        nValues, = unpack('i',data)
        #print "nValues = ",nValues/4
        ints = self.readInts(nValues/4)
        self.skip(4)
        #print "ints = ",ints
        return ints
        
