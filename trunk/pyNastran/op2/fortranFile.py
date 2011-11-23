#import os

import sys
from struct import unpack,pack

#pyNastran
from op2Errors import *

class FortranFile(object):
    def __init__(self):
        self.endian = '>'
        self.bufferSize = 65535
    
    def setEndian(self,endian):
        self.endian = endian

    def readHollerith(self):
        self.skip(4)  # weird hollerith
        
    def readHeader(self,expected,debug=True):
        """
        a header is defined as (4,i,4), where i is an integer
        """
        #self.printSection(60)
        #data = self.op2.read(12)
        ints = self.readFullIntBlock()
        #print "header ints = %s" %(repr(ints))
        #self.n += 12*4
        
        if len(ints)==5:
            print "   might be a buffer block..."
            ints = self.readFullIntBlock()
            if debug:
                self.op2Debug.write('bufferBlcok = |%s|\n' %(str(ints)))
            
        elif len(ints)==0:
            return None
        
        assert ints[0]==ints[2]==4,"header ints=(%s) expected=%s\n" %(str(ints[0:5]),expected)
        if debug:
            self.op2Debug.write('[4,%s,4]\n' %(ints[1]))
        return ints[1]

    def readString(self,nData):
        """
        reads nCharacters that are assumed to be a string
        """
        data = self.op2.read(nData)
        string = ''.join(self.getStrings(data))
        self.op2Debug.write('|%s|\n' %(str(string)))
        self.n += nData
        return string

    #def readString(self,nData):
    #    data = self.op2.read(nData)
    #    self.n += nData
    #    words = self.getStrings(data)
    #    return ''.join(words)

    def readInts(self,nInts,debug=True):
        """
        reads nIntegers
        """
        nData = 4*nInts
        #print "nData = ",nData
        data = self.op2.read(nData)

        iFormat = 'i'*nInts
        ints = unpack(iFormat,data)
        if debug:
            self.op2Debug.write('|%s|\n' %(str(ints)))

        self.n+=nData
        return ints

    def readDoubles(self,nData,debug=True):
        """
        reads nDoubles
        """
        data = self.op2.read(nData)
        self.n += nData
        doubles = self.getDoubles(data)
        if debug:
            self.op2Debug.write('|%s|\n' %(str(doubles)))
        return doubles

    def readFloats(self,nData,debug=True):
        """
        reads nFloats
        """
        data = self.op2.read(nData)
        self.n += nData
        floats = self.getFloats(data)
        if debug:
            self.op2Debug.write('|%s|\n' %(str(floats)))
        return floats

    def getStrings(self,data):
        """
        unpacks a data set into a series of characters
        """
        n = len(data)
        sFormat = 's'*n
        strings = unpack(sFormat,data)
        #self.op2Debug.write('|%s|\n' %(str(strings)))
        return strings

    def getInts(self,data,debug=True):
        """
        unpacks a data set into a series of ints
        """
        n = len(data)
        nInts = n/4
        #print "nInts = ",nInts
        iFormat = 'i'*nInts
        ints = unpack(iFormat,data[:nInts*4])
        if debug:
            self.op2Debug.write('|%s|\n' %(str(ints)))
        return ints

    def getLongs(self,data):
        """
        unpacks a data set into a series of longs
        """
        n = len(data)
        nLongs = n/4
        #print "nLongs = ",nLongs
        #a = pack('l',200)
        #print "len(a) = ",len(a)
        
        LFormat = 'l'*nLongs
        longs = unpack(LFormat,data[:nLongs*4])
        return longs

    def getFloats(self,data):
        """
        unpacks a data set into a series of floats
        """
        n = len(data)
        nFloats = n/4
        fFormat = 'f'*nFloats
        ints = unpack(fFormat,data[:nFloats*4])
        return ints

    def getDoubles(self,data):
        """
        unpacks a data set into a series of doubles
        """
        n = len(data)
        nDoubles = n/8
        dFormat = 'd'*nDoubles
        ints = unpack(dFormat,data[:nDoubles*8])
        return ints
    
    def printBlock(self,data):
        """
        prints a data set in int/float/double/string format to
        determine table info.  doesn't move cursor.
        """
        ints    = self.getInts(data)
        #longs   = self.getLongs(data)
        floats  = self.getFloats(data)
        doubles = self.getDoubles(data)
        strings = self.getStrings(data)
        print "ints    = ",ints
        #print "longs   = ",longs
        print "floats  = ",floats
        #print "doubles = ",doubles
        print "strings = |%r|" %(''.join(strings))
        print "nWords = ",len(data)/4

    def getData(self,n):
        """
        gets a data set of length N
        """
        assert n>0
        #assert self.op2.tell()==self.n,'tell=%s n=%s' %(self.op2.tell(),self.n)
        data = self.op2.read(n)
        self.n+=n
        #print "n =",n
        #assert self.op2.tell()==self.n,'tell=%s n=%s' %(self.op2.tell(),self.n)
        return data

    def getBlockIntEntry(self,data,n):
        """
        given a data set, grabs the nth word and casts it as an integer
        """
        data2 = data[4*(n-1):4*(n-1)+4]
        return unpack('i',data2)[0]

    def printSection(self,nBytes):
        """
        prints data, but doesn't move the cursor
        """
        data = self.op2.read(nBytes)
        self.printBlock(data)
        self.op2.seek(self.n)
    
    def skip(self,n):
        """skips nBits"""
        self.op2Debug.write('skipping\n')
        self.printSection(4)
        self.op2Debug.write('skipped\n')
        #print "\n--SKIP--"
        #print "tell = ",self.op2.tell()
        #print "n = ",n
        #print "self.n = ",self.n
        self.n += n
        self.op2.seek(self.n)
        #print "*tell = ",self.op2.tell()
        #print "\n"

    def scan(self,n):
        """same as skip, but actually reads the data instead of using seek"""
        self.op2Debug.write('skipping %s\n' %(n))
        self.printSection(4)
        self.op2Debug.write('skipped\n')
        data = self.op2.read(n)
        self.n+=n

    def getTableCode(self,expected=None,debug=True):
        tableCode = self.readHeader(expected,debug)
        return tableCode

    def getMarker(self,expected=None,debug=True):
        tableCode = self.readHeader(expected,debug)
        return tableCode
        
    def readMarker(self,expected=None):
        return self.getMarker(expected)
        
    def readMarkers(self,markers,tableName=None,debug=False):
        """
        reads a set of predefined markers e.g. [-4,1,0]
        and makes sure it is correct
        """
        for marker in markers:
            tableCode = self.readHeader(marker,debug)
            if tableCode==None:
                return
            assert marker==tableCode,'tableName=%s found=%s expected=%s leftover=%s' %(tableName,tableCode,marker,self.printSection(40))
        ###
        msg = ''
        for i in markers:
            msg += '[4,'+str(i)+',4] + '
        self.op2Debug.write(msg[:-3]+'\n')
        print "@markers = ",markers
        print ""

    def getNMarkers(self,nMarkers,rewind=False):
        """gets the next N markers, verifies they're correct"""
        markers = []
        for iMarker in range(nMarkers):
            tableCode = self.readHeader(None)
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
            #sys.exit(expectedMarkers)
            return False
        else:
            raise RuntimeError('this should never happen...invalid markers...expected=%s markers=%s' %(expectedMarkers,markers))

    def goto(self,n):
        #print "goto n = |%s|" %(n)
        assert n>0
        self.op2.seek(n)

    def readBlock(self):
        """
        reads a fortran formatted data block
        nWords  data1 data2 data3 nWords
        """
        data = self.op2.read(4)
        if len(data)==0:
            raise EndOfFileError("data=('')")
        nValues, = unpack('i',data)
        self.n+=4
        data = self.op2.read(nValues)
        self.n+=nValues+4
        self.goto(self.n)
        return data

    def readFullBlock(self):
        """
        reads a fortran formatted data block
        nWords  data1 data2 data3 nWords
        includes nWords in the output
        """
        data = self.op2.read(4)
        nValues, = unpack('i',data)
        self.n+=4
        data = self.op2.read(nValues)
        self.n+=nValues+4
        self.goto(self.n)

    def readFullIntBlock(self):
        """
        reads a fortran formatted block
        assumes that the data is made up of integers only
        """
        """
        reads a fortran formatted data block
        nWords  data1 data2 data3 nWords
        includes nWords in the output
        """
        data = self.op2.read(4)
        if len(data)==0:
            print "found the end of the file..."
            return []
        nValues, = unpack('i',data)
        self.n+=4
        data = self.op2.read(nValues)
        self.n+=nValues+4
        self.goto(self.n)

        #nInts = len(data)/4
        nInts = len(data)/4
        #print "**nInts = ",nInts
        ints = unpack('i'*nInts,data)
        return [nValues]+list(ints)+[nValues]

    def readStringBlock(self,debug=True):
        """
        reads a fortran formatted block
        assumes that the data is made up of characters only
        """
        data = self.readBlock()
        nLetters = len(data)
        letters = unpack('s'*nLetters,data)  ## @todo should this be c instead of s???
        word = ''.join(letters)
        #print "word = |%s|" %(word)
        print "nLetters=%s word=|%s|" %(nLetters,word)
        if debug:
            self.op2Debug.write('|%s|\n' %(str(word)))
        return word

    def readIntBlock(self):
        """
        reads a fortran formatted block
        assumes that the data is made up of integers only
        """
        data = self.readBlock()
        nInts = len(data)/4
        #print "**nInts = ",nInts
        ints = unpack('i'*nInts,data)
        return ints

    def readFloatBlock(self):
        """
        reads a fortran formatted block
        assumes that the data is made up of floats only
        """
        data = self.readBlock()
        nFloats = len(data)/4
        floats = unpack('f'*nFloats,data)
        return floats

    def readDoubleBlock(self):
        """
        reads a fortran formatted block
        assumes that the data is made up of doubles only
        """
        data = self.readBlock()
        nDoubles = len(data)/8
        doubles = unpack('d'*nDoubles,data)
        return doubles

    def rewind(self,n):
        """
        rewinds the file nBytes
        @warning
            doesnt support a full rewind, only a partial
        """
        self.n -= n
        self.op2.seek(self.n)

    def readTableName(self,rewind=True,stopOnFailure=True):
        """
        peeks into a table to check it's name
        """
        debug = True
        if rewind:
            debug = False
        n = self.n
        try:
            #print ""
            self.readMarkers([0,2],debug)
            word = self.readStringBlock(debug)
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
        except:
            if rewind and not stopOnFailure:
                self.n = n
                self.op2.seek(n)
                return
            raise
            ###
        ###

    def skipNextTable(self,bufferSize=10000):
        """
        skips a table
        @todo fix bugs
        """
        tableName = self.readTableName(rewind=False) # GEOM1
        self.tableInit(tableName)
        print "skippingTable |%s|" %(tableName)
        print "self.n = ",self.n

        self.readMarkers([-1,7],tableName)
        
        dataPack = (4,1,4,  4,0,4,  4,0,4)  # marks the end of the table
        binaryData = pack('iiiiiiiii',*dataPack)
        #[1,0,0] marks the end of the table
        
        i = 0
        error = 80
        n = self.n
        endIndex = -1
        data = "dummy"
        while endIndex== -1 and len(data)>0:
            data     = self.op2.read(bufferSize+error)
            endIndex = data.find(binaryData)
            
            self.op2.seek(n+i*bufferSize)
            i+=1
        
        #print "i = ",i
        assert endIndex>0,'couldnt find the end of the table'
        self.n = self.n+(i-1)*bufferSize + endIndex # 36 so it gets to the end of the table markersNext=[0] or [2]
        
        n = self.n
        self.n += 36  ## @todo sometimes this is needed
        #ints = self.readIntBlock()
        #print "*?*ints = ",ints
        #if len(ints)==0:
        #    pass

        self.op2.seek(self.n)
        "self.op2.tell() = ",self.op2.tell()
        #self.printSection(200)
        marker = self.getMarker()
        print "marker = ",marker
        if marker==2:
            isAnotherTable = True
        else:# marker=0
            isAnotherTable = False
        isAnotherTable = True
        ###
        #print "isAnotherTable = ",isAnotherTable
        self.n -= 24  # subtract off the header [0,2] or [0,0]
        self.op2.seek(self.n)
        print "self.n = ",self.n
        print "---table %s is skipped---" %(tableName)
        
        return isAnotherTable

    def hasMoreTables(self):
        #print self.printSection(120)
        try:
            marker1 = self.getMarker('[4,0,4]')
            marker2 = self.getMarker('[4,0,4] or [4,2,4]')

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
    