from __future__ import (nested_scopes, generators, division, absolute_import,
                        print_function, unicode_literals)
#import sys
from struct import unpack, pack



class FortranFile(object):
    def __init__(self):
        ## the endian of the processor (typically '<' for Windows/Linux/Mac,
        ##                              '>' for old HPCs)
        ## currently does nothing
        self.endian = '<'
        ## currently does nothing
        self.bufferSize = 65535

    def setEndian(self, endian='<'):
        """
        sets the endian
        @todo hasnt been implemented
        """
        self.endian = endian

    def readHollerith(self):
        """
        doesnt really read a hollerith, it's an integer
        of value=528 which corresponds to the length of
        iTable=3
        """
        self.skip(4)

    def readHeader(self, expected=None, debug=True):
        """
        a header is defined as (4,i,4), where i is an integer
        """
        #print(self.printSection(60))
        #data = self.op2.read(12)
        ints = self.readFullIntBlock()
        #print("header ints = %s" % (repr(ints)))
        #self.n += 12*4

        if len(ints) == 5:
            #print "   buffer block..."
            ## flag to help know if a buffer was found
            self.hasBuffer = True
            ints = self.readFullIntBlock()
            if debug and self.makeOp2Debug:
                self.op2Debug.write('bufferBlock = |%s|\n' % (str(ints)))

        elif len(ints) == 0:
            return None

        if not(ints[0] == ints[2] == 4):
            msg = "pyNastran reading failed because an improperly formatted (or unsupported) table is in the OP2.\n"
            msg += "If you remove the offending table (make sure you're using PARAM,POST,-1 first) the code should work.\n"
            msg += "header ints=(%s) expected=%s\n" % (
                str(ints[0:5]), expected)
            msg += 'tableName=|%s|' % (self.tableName)
            raise SyntaxError("Invalid Marker: %s", msg)

        #if ints[1]==2:  # buffer???
        #    ints = self.readFullIntBlock()
        #    print("bufferInts1=%s" % (ints))
        #    ints = self.readFullIntBlock()
        #    print("bufferInts2=%s" % (ints))
        #print("marker=%s" % (ints[1]))
        if debug and self.makeOp2Debug:
            self.op2Debug.write('[4,%s,4]\n' % (ints[1]))
        return ints[1]

    def readString(self, nData):
        """
        reads nCharacters that are assumed to be a string
        """
        data = self.op2.read(nData)
        string = self.getStrings(data)
        if self.makeOp2Debug:
            self.op2Debug.write('|%s|\n' % (str(string)))
        self.n += nData
        return string

    #def readString(self, nData):
    #    data = self.op2.read(nData)
    #    self.n += nData
    #    words = self.getStrings(data)
    #    return ''.join(words)

    def readInts(self, nInts, debug=True):
        """
        reads a list of nIntegers
        @param self the object pointer
        @param nInts the number of ints to read
        @param debug for developer: debug combined with makeOp2Debug
        """
        nData = 4 * nInts
        #print "nData = ",nData
        data = self.op2.read(nData)

        iFormat = str(nInts) + 'i'
        iFormat = bytes(iFormat)
        ints = unpack(iFormat, data)
        if debug and self.makeOp2Debug:
            self.op2Debug.write('|%s|\n' % (str(ints)))

        self.n += nData
        return ints

    def readDoubles(self, nData, debug=True):
        """
        reads a list of nDoubles
        @param self the object pointer
        @param nData the number of doubles to read
        @param debug for developer: debug combined with makeOp2Debug
        """
        data = self.op2.read(nData)
        self.n += nData
        doubles = self.getDoubles(data)
        if debug and self.makeOp2Debug:
            self.op2Debug.write('|%s|\n' % (str(doubles)))
        return doubles

    def readFloats(self, nData, debug=True):
        """
        reads nFloats
        """
        data = self.op2.read(nData)
        self.n += nData
        floats = self.getFloats(data)
        if debug and self.makeOp2Debug:
            self.op2Debug.write('|%s|\n' % (str(floats)))
        return floats

    def getStrings(self, data):
        """
        unpacks a data set into a series of characters
        """
        n = len(data)
        iFormat = str(n) + 's'
        iFormat = bytes(iFormat)
        strings, = unpack(iFormat, data)
        return strings  # .encode('utf-8')

    def getStrings2(self, data, endian):
        """
        unpacks a data set into a series of characters
        """
        n = len(data)
        iFormat = endian + str(n) + 's'
        iFormat = bytes(iFormat)
        strings = unpack(iFormat, data)
        return strings

    def getInts(self, data, debug=True):
        """
        unpacks a data set into a series of ints
        """
        n = len(data)
        nInts = n // 4
        #print "nInts = ",nInts
        iFormat = str(nInts) + 'i'
        iFormat = bytes(iFormat)
        ints = unpack(iFormat, data[:nInts * 4])
        return ints

    def getInts2(self, data, endian, debug=True):
        """
        unpacks a data set into a series of ints
        """
        n = len(data)
        nInts = n // 4
        iFormat = endian + str(nInts) + 'i'
        iFormat = bytes(iFormat)
        ints = unpack(iFormat, data[:nInts * 4])
        return ints

    def getLongs(self, data):
        """
        unpacks a data set into a series of longs
        """
        n = len(data)
        nLongs = n // 4
        #print "nLongs = ", nLongs
        #a = pack('l', 200)
        #print "len(a) = ", len(a)

        iFormat = str(nLongs) + 'l'
        iFormat = bytes(iFormat)
        longs = unpack(iFormat, data[:nLongs * 4])
        return longs

    def getFloats(self, data):
        """
        unpacks a data set into a series of floats
        """
        n = len(data)
        nFloats = n // 4
        iFormat = str(nFloats) + 'f'
        iFormat = bytes(iFormat)
        ints = unpack(iFormat, data[:nFloats * 4])
        return ints

    def getFloats2(self, data, endian):
        """
        unpacks a data set into a series of floats
        """
        n = len(data)
        nFloats = n // 4
        iFormat = endian + str(nFloats) + 'f'
        iFormat = bytes(iFormat)
        ints = unpack(iFormat, data[:nFloats * 4])
        return ints

    def getDoubles(self, data):
        """
        unpacks a data set into a series of doubles
        """
        n = len(data)
        nDoubles = n // 8
        iFormat = str(nDoubles) + 'd'
        iFormat = bytes(iFormat)
        ints = unpack(iFormat, data[:nDoubles * 8])
        return ints

    def printBlock(self, data, nMax=200):
        """
        prints a data set in int/float/double/string format to
        determine table info.  doesn't move cursor.
        @note this is a great function for debugging
        """
        data2 = data
        nData = len(data)
        #if nData>nMax:
            #print("oops!  too much data...limiting to %s bytes.  nData=%s" % (nMax,nData))
            #data2 = data[:nMax]

        msg = ''
        ints = self.getInts(data2)
        #longs   = self.getLongs(data2)
        floats = self.getFloats(data2)
        #doubles = self.getDoubles(data2)
        strings = self.getStrings(data2)
        msg += "n       = %s\n" % (self.n)
        msg += "ints    = %s\n" % (str(ints))
        #msg += "longs  = %s\n" % (longs)
        msg += "floats  = %s\n" % (str(floats))
        #msg += "doubles = %s\n" % (str(doubles))
        msg += "strings = |%r|\n" % (strings)
        msg += "nWords  = %s\n" % (len(data) // 4)
        #msg += "tell    = %s\n" % (self.op2.tell())
        return msg

    def printBlock2(self, data, endian):
        """
        prints a data set in int/float/double/string format to
        determine table info.  doesn't move cursor.
        @note this is a great function for debugging
        """
        msg = ''
        ints = self.getInts2(data, endian)
        #longs   = self.getLongs(data)
        floats = self.getFloats2(data, endian)
        #doubles = self.getDoubles(data)
        strings = self.getStrings2(data, endian)
        msg += "ints    = %s\n" % (str(ints))
        #msg += "longs  = %s\n" % (longs)
        msg += "floats  = %s\n" % (str(floats))
        #msg += "doubles = %s\n" % (doubles)
        msg += "strings = |b%r|\n" % (''.join(strings))
        msg += "nWords  = %s\n" % (len(data) // 4)
        #msg += "tell    = %s\n" % (self.op2.tell())
        return msg

    def getData(self, n):
        """
        gets a data set of length N
        """
        if n <= 0:
            raise RuntimeError('Zero Buffer Error')

        #assert self.op2.tell()==self.n,'tell=%s n=%s' % (self.op2.tell(),self.n)
        data = self.op2.read(n)
        self.n += n
        #print "n =",n
        #assert self.op2.tell()==self.n,'tell=%s n=%s' % (self.op2.tell(),self.n)
        return data

    def readData(self, n):
        return self.getData(n)

    def getBlockIntEntry(self, data, n):
        """
        given a data set, grabs the nth word and casts it as an integer
        """
        data2 = data[4 * (n - 1):4 * (n - 1) + 4]
        iFormat = 'i'
        iFormat = bytes(iFormat)
        return unpack(iFormat, data2)[0]

    def printSection(self, nBytes):
        """
        prints data, but doesn't move the cursor
        @param self the object pointer
        @param nBytes the number of bytes to print the data specs on
        @retval msg ints/floats/strings of the next nBytes (handles poorly sized nBytes; uncrashable :) )
        @note this the BEST function when adding new cards/tables/debugging
        """
        data = self.op2.read(nBytes)
        msg = self.printBlock(data)
        self.op2.seek(self.n)
        return msg

    def printSection2(self, nBytes, endian):
        """
        prints data, but doesn't move the cursor
        @param self the object pointer
        @param nBytes the number of bytes to print the data specs on
        @retval msg ints/floats/strings of the next nBytes (handles poorly sized nBytes; uncrashable :) )
        @note this the BEST function when adding new cards/tables/debugging
        """
        data = self.op2.read(nBytes)
        msg = self.printBlock2(data, endian)
        self.op2.seek(self.n)
        return msg

    def skip(self, n):
        """skips nBits"""
        self.n += n
        self.op2.seek(self.n)

    def scan(self, n):
        """same as skip, but actually reads the data instead of using seek"""
        self.op2.read(n)
        self.n += n

    def getTableCode(self, expected=None, debug=True):
        tableCode = self.readHeader(expected, debug)
        return tableCode

    def getMarker(self, expected=None, debug=True):
        tableCode = self.readHeader(expected, debug)
        return tableCode

    def readMarker(self, expected=None):
        return self.getMarker(expected)

    def readMarkers(self, markers, tableName=None, debug=False, printErrorOnFailure=True):
        """
        Reads a set of predefined markers e.g. [-3,1,0]
        and makes sure it is correct.

        A marker (e.g. a -3) is a series of 3 integers [4,-3,4].  Typically 3
        markers are put together (e.g. [-3,1,0]) such that the integers are
        [4,-3,4, 4,1,4, 4,0,4] to mark important parts of the table.

        Markers will "increment" during table reading, such that the first marker
        is [-1,1,0], then [-2,1,0], etc.  Tables will end (or really the next table starts)
        when a [-1,1,0] or a [0,1,0] marker is found.

        # Verify the following statement...
        Occassionally, buffer markers will be embedded inside the
        marker [-3,1,0], (e.g. [4,2^16,4] <- the default BUFFSIZE), which make
        reading the marker more difficult.
        """
        #print("markers = ",markers)
        foundMarkers = []
        for marker in markers:
            tableCode = self.readHeader(marker, debug)
            #if tableCode==2:
            #    tableCode = self.readHeader(marker,debug)

            #print("tableCode=",tableCode)
            if tableCode is None:
                return
            if marker != tableCode:
                msg = ''
                if printErrorOnFailure:
                    msg = '\nmarkers=%s foundMarkers=%s\n' % (
                        markers, foundMarkers)
                    msg += 'tableName=%s found=%s expected=%s leftover=%s' % (tableName, tableCode, marker, self.printSection(40))
                    #print(msg)
                raise SyntaxError("Invalid Markers: %s", msg)
            foundMarkers.append(marker)
            ###
        ###
        msg = ''
        for i in markers:
            msg += '[4,' + str(i) + ',4] + '
        if self.makeOp2Debug:
            self.op2Debug.write(msg[:-3] + '\n')
        if debug:
            self.log.debug("@markers = %s" % (markers))
            self.log.debug("")
        ###

    def getNMarkers(self, nMarkers, rewind=False):
        """gets the next N markers, verifies they're correct"""
        markers = []
        for iMarker in xrange(nMarkers):
            tableCode = self.readHeader(None)
            markers.append(tableCode)
        ###
        if rewind:
            self.n -= 12 * nMarkers
            self.op2.seek(self.n)

        return markers

    def isTableDone(self, expectedMarkers):
        markers = self.getNMarkers(len(expectedMarkers), rewind=True)
        #print "getMarkers = ",markers

        if markers == [-1, 7]:
            return True
        elif markers == expectedMarkers:
            #sys.exit(expectedMarkers)
            return False
        else:
            raise RuntimeError('this should never happen...invalid markers...expected=%s markers=%s' % (expectedMarkers, markers))

    def goto(self, n):
        """
        jumps to position n in the file
        @param self the object pointer
        @param n the position to goto
        @note n>0
        """
        #print "goto n = |%s|" % (n)
        assert n > 0
        self.n = n
        self.op2.seek(n)

    def readBlock(self):
        """
        reads a fortran formatted data block
        nWords  data1 data2 data3 nWords
        """
        data = self.op2.read(4)
        if len(data) == 0:
            raise EOFError("data=('')")

        iFormat = 'i'
        iFormat = bytes(iFormat)
        nValues, = unpack(iFormat, data)
        self.n += 4
        data = self.op2.read(nValues)
        self.n += nValues + 4
        self.goto(self.n)
        return data

    def readFullBlock(self):
        """
        reads a fortran formatted data block
        nWords  data1 data2 data3 nWords
        includes nWords in the output
        """
        data = self.op2.read(4)
        iFormat = 'i'
        iFormat = bytes(iFormat)
        nValues, = unpack(iFormat, data)
        self.n += 4
        data = self.op2.read(nValues)
        self.n += nValues + 4
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
        if len(data) == 0:
            self.log.debug("found the end of the file...")
            return []
        iFormat = 'i'
        iFormat = bytes(iFormat)
        nValues, = unpack(iFormat, data)
        self.n += 4
        data = self.op2.read(nValues)
        self.n += nValues + 4
        self.goto(self.n)

        nInts = len(data) // 4
        iFormat = str(nInts) + 'i'
        iFormat = bytes(iFormat)
        ints = unpack(iFormat, data)
        return [nValues] + list(ints) + [nValues]

    def readStringBlock(self, debug=True):
        """
        reads a fortran formatted block
        assumes that the data is made up of characters only
        """
        data = self.readBlock()
        nLetters = len(data)
        iFormat = str(nLetters) + 's'
        iFormat = bytes(iFormat)
        word, = unpack(iFormat, data)

        #print "word = |%s|" % (word)
        #print "nLetters=%s word=|%s|" % (nLetters,word)
        if debug and self.makeOp2Debug:
            self.op2Debug.write('|%s|\n' % (str(word)))
        return word

    def readIntBlock(self):
        """
        reads a fortran formatted block
        assumes that the data is made up of integers only
        """
        data = self.readBlock()
        nInts = len(data) // 4
        iFormat = str(nInts) + 'i'
        iFormat = bytes(iFormat)
        ints = unpack(iFormat, data)
        return ints

    def readFloatBlock(self):
        """
        reads a fortran formatted block
        assumes that the data is made up of floats only
        """
        data = self.readBlock()
        nFloats = len(data) // 4
        iFormat = str(nFloats) + 'f'
        iFormat = bytes(iFormat)
        floats = unpack(iFormat, data)
        return floats

    def readDoubleBlock(self):
        """
        reads a fortran formatted block
        assumes that the data is made up of doubles only
        """
        data = self.readBlock()
        nDoubles = len(data) // 8
        iFormat = str(nDoubles) + 'd'
        iFormat = bytes(iFormat)
        doubles = unpack(iFormat, data)
        return doubles

    def rewind(self, n):
        """
        rewinds the file nBytes
        @warning
            doesnt support a full rewind, only a partial
        """
        self.n -= n
        self.op2.seek(self.n)

    def readTableName(self, rewind=True, debug=True, stopOnFailure=True):
        """
        peeks into a table to check it's name
        """
        #debug = True
        if rewind:
            debug = False
        n = self.n
        try:
            #print ""
            self.readMarkers([0, 2], debug)
            word = self.readStringBlock(debug)
            #print("*word = |%r|" % (word))

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
            tableName = word.strip()
            return tableName.decode('utf-8')
        except:
            if rewind and not stopOnFailure:
                self.n = n
                self.op2.seek(n)
                return
            raise
            ###
        ###

    def skipNextTable(self, bufferSize=10000):
        """
        skips a table
        @todo fix bugs
        """
        tableName = self.readTableName(rewind=False)  # GEOM1
        self.tableInit(tableName)
        self.log.debug("skippingTable |%s|" % (tableName))
        self.log.debug("self.n = %s" % (self.n))

        self.readMarkers([-1, 7], tableName)

        dataPack = (4, 1, 4, 4, 0, 4, 4, 0, 4)  # marks the end of the table
        binaryData = pack('9i', *dataPack)
        #[1,0,0] marks the end of the table

        i = 0
        error = 80
        n = self.n
        endIndex = -1
        data = "dummy"
        while endIndex == -1 and len(data) > 0:
            data = self.op2.read(bufferSize + error)
            endIndex = data.find(binaryData)

            self.op2.seek(n + i * bufferSize)
            i += 1

        #print "i = ",i
        assert endIndex > 0, 'couldnt find the end of the table'
        self.n = self.n + (i - 1) * bufferSize + endIndex  # 36 so it gets to the end of the table markersNext=[0] or [2]

        n = self.n
        self.n += 36  ## @todo sometimes this is needed
        #ints = self.readIntBlock()
        #print "*?*ints = ",ints
        #if len(ints)==0:
        #    pass

        self.op2.seek(self.n)
        self.log.debug("self.op2.tell() = %s" % (self.op2.tell()))
        #self.printSection(200)
        marker = self.getMarker()
        #print "marker = ",marker
        if marker == 2:
            isAnotherTable = True
        else:  # marker=0
            isAnotherTable = False
        isAnotherTable = True
        ###
        #print "isAnotherTable = ",isAnotherTable
        self.n -= 24  # subtract off the header [0,2] or [0,0]
        self.op2.seek(self.n)
        self.log.debug("self.n = %s" % (self.n))
        self.log.debug("---table %s is skipped---" % (tableName))

        return isAnotherTable

    def hasMoreTables(self):
        #print self.printSection(120)
        try:
            marker1 = self.getMarker('[4,0,4]')
            marker2 = self.getMarker('[4,0,4] or [4,2,4]')

            marker = [marker1, marker2]
            #print "marker = ",marker
            if marker == [0, 2]:
                isAnotherTable = True
            else:  # marker=0
                isAnotherTable = False
            ###
            #print "isAnotherTable = ",isAnotherTable
            self.n -= 24  # subtract off the header [0,2] or [0,0]
            self.op2.seek(self.n)
        except IndexError:
            isAnotherTable = False
        return isAnotherTable
