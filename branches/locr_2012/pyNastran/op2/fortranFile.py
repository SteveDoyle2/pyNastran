#pylint:  disable=C0103
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
        self.buffer_size = 65535
        self.table_name = None
        self.op2 = None
        self.make_op2_debug = False
        self.op2Debug = None
        self.log = None
        self.n = 0

    def set_endian(self, endian='<'):
        """
        Sets the endian
        @todo hasnt been implemented
        """
        self.endian = endian

    def read_hollerith(self):
        """
        Doesnt really read a hollerith, it's an integer
        of value=528 which corresponds to the length of
        iTable=3
        """
        self.skip(4)

    def read_header(self, expected=None, debug=True):
        """
        A header is defined as (4,i,4), where i is an integer
        """
        #print(self.print_section(60))
        #data = self.op2.read(12)
        ints = self.read_full_int_block()
        #print("header ints = %s" % (repr(ints)))
        #self.n += 12*4

        if len(ints) == 5:
            #print "   buffer block..."
            ## flag to help know if a buffer was found
            self.hasBuffer = True
            ints = self.read_full_int_block()
            if debug and self.make_op2_debug:
                self.op2Debug.write('bufferBlock = |%s|\n' % str(ints))

        elif len(ints) == 0:
            return None

        if not(ints[0] == ints[2] == 4):
            msg = ("pyNastran reading failed because an improperly formatted "
                   "(or unsupported) table is in the OP2.\n"
                   "If you remove the offending table (make sure you're using "
                   "PARAM,POST,-1 first) the code should work.\n"
                   "header ints=(%s) expected=%s\n" % (str(ints[0:5]),
                                 expected))
            msg += 'table_name=|%s|' % self.table_name
            raise SyntaxError("Invalid Marker: %s" % msg)

        #if ints[1]==2:  # buffer???
        #    ints = self.read_full_int_block()
        #    print("bufferInts1=%s" % (ints))
        #    ints = self.read_full_int_block()
        #    print("bufferInts2=%s" % (ints))
        #print("marker=%s" % (ints[1]))
        if debug and self.make_op2_debug:
            self.op2Debug.write('[4,%s,4]\n' % ints[1])
        return ints[1]

    def read_string(self, nData):
        """
        Reads nCharacters that are assumed to be a string
        """
        data = self.op2.read(nData)
        string = self.get_strings(data)
        if self.make_op2_debug:
            self.op2Debug.write('|%s|\n' % str(string))
        self.n += nData
        return string

    #def read_string(self, nData):
    #    data = self.op2.read(nData)
    #    self.n += nData
    #    words = self.getStrings(data)
    #    return ''.join(words)

    def read_ints(self, nInts, debug=True):
        """
        Reads a list of nIntegers
        @param self the object pointer
        @param nInts the number of ints to read
        @param debug for developer: debug combined with make_op2_debug
        """
        nData = 4 * nInts
        #print "nData = ",nData
        data = self.op2.read(nData)

        iFormat = str(nInts) + 'i'
        iFormat = bytes(iFormat)
        ints = unpack(iFormat, data)
        if debug and self.make_op2_debug:
            self.op2Debug.write('|%s|\n' % str(ints))
        self.n += nData
        return ints

    def read_doubles(self, nData, debug=True):
        """
        Reads a list of nDoubles
        @param self the object pointer
        @param nData the number of doubles to read
        @param debug for developer: debug combined with make_op2_debug
        """
        data = self.op2.read(nData)
        self.n += nData
        doubles = self.get_doubles(data)
        if debug and self.make_op2_debug:
            self.op2Debug.write('|%s|\n' % str(doubles))
        return doubles

    def read_floats(self, nData, debug=True):
        """
        Reads nFloats
        """
        data = self.op2.read(nData)
        self.n += nData
        floats = self.get_floats(data)
        if debug and self.make_op2_debug:
            self.op2Debug.write('|%s|\n' % str(floats))
        return floats

    def get_strings(self, data):
        """
        Unpacks a data set into a series of characters
        """
        n = len(data)
        iFormat = str(n) + 's'
        iFormat = bytes(iFormat)
        strings, = unpack(iFormat, data)
        return strings  # .encode('utf-8')

    def get_strings2(self, data, endian):
        """
        Unpacks a data set into a series of characters
        """
        n = len(data)
        iFormat = endian + str(n) + 's'
        iFormat = bytes(iFormat)
        strings = unpack(iFormat, data)
        return strings

    def get_ints(self, data):
        """
        Unpacks a data set into a series of ints
        """
        n = len(data)
        nInts = n // 4
        #print "nInts = ",nInts
        iFormat = str(nInts) + 'i'
        iFormat = bytes(iFormat)
        ints = unpack(iFormat, data[:nInts * 4])
        return ints

    def get_ints2(self, data, endian):
        """
        Unpacks a data set into a series of ints
        """
        n = len(data)
        nInts = n // 4
        iFormat = endian + str(nInts) + 'i'
        iFormat = bytes(iFormat)
        ints = unpack(iFormat, data[:nInts * 4])
        return ints

    def get_longs(self, data):
        """
        Unpacks a data set into a series of longs
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

    def get_floats(self, data):
        """
        Unpacks a data set into a series of floats
        """
        n = len(data)
        nFloats = n // 4
        iFormat = str(nFloats) + 'f'
        iFormat = bytes(iFormat)
        ints = unpack(iFormat, data[:nFloats * 4])
        return ints

    def get_floats2(self, data, endian):
        """
        Unpacks a data set into a series of floats
        """
        n = len(data)
        nFloats = n // 4
        iFormat = endian + str(nFloats) + 'f'
        iFormat = bytes(iFormat)
        ints = unpack(iFormat, data[:nFloats * 4])
        return ints

    def get_doubles(self, data):
        """
        Unpacks a data set into a series of doubles
        """
        n = len(data)
        nDoubles = n // 8
        iFormat = str(nDoubles) + 'd'
        iFormat = bytes(iFormat)
        ints = unpack(iFormat, data[:nDoubles * 8])
        return ints

    def print_block(self, data, nMax=200):
        """
        Prints a data set in int/float/double/string format to
        determine table info.  doesn't move cursor.
        @note this is a great function for debugging
        """
        data2 = data
        nData = len(data)
        #if nData>nMax:
            #print("oops!  too much data...limiting to %s bytes.  nData=%s" % (nMax,nData))
            #data2 = data[:nMax]

        msg = ''
        ints = self.get_ints(data2)
        #longs   = self.getLongs(data2)
        floats = self.get_floats(data2)
        #doubles = self.getDoubles(data2)
        strings = self.get_strings(data2)
        msg += "n       = %s\n" % (self.n)
        msg += "ints    = %s\n" % (str(ints))
        #msg += "longs  = %s\n" % (longs)
        msg += "floats  = %s\n" % (str(floats))
        #msg += "doubles = %s\n" % (str(doubles))
        msg += "strings = |%r|\n" % (strings)
        msg += "nWords  = %s\n" % (len(data) // 4)
        #msg += "tell    = %s\n" % (self.op2.tell())
        return msg

    def print_block2(self, data, endian):
        """
        Prints a data set in int/float/double/string format to
        determine table info.  doesn't move cursor.
        @note this is a great function for debugging
        """
        msg = ''
        ints = self.get_ints2(data, endian)
        #longs   = self.getLongs(data)
        floats = self.get_floats2(data, endian)
        #doubles = self.getDoubles(data)
        strings = self.get_strings2(data, endian)
        msg += "ints    = %s\n" % (str(ints))
        #msg += "longs  = %s\n" % (longs)
        msg += "floats  = %s\n" % (str(floats))
        #msg += "doubles = %s\n" % (doubles)
        msg += "strings = |b%r|\n" % (''.join(strings))
        msg += "nWords  = %s\n" % (len(data) // 4)
        #msg += "tell    = %s\n" % (self.op2.tell())
        return msg

    def get_data(self, n):
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

    def read_data(self, n):
        return self.get_data(n)

    def get_block_int_entry(self, data, n):
        """
        given a data set, grabs the nth word and casts it as an integer
        """
        data2 = data[4 * (n - 1):4 * (n - 1) + 4]
        iFormat = 'i'
        iFormat = bytes(iFormat)
        return unpack(iFormat, data2)[0]

    def print_section(self, nBytes):
        """
        Prints data, but doesn't move the cursor
        @param self the object pointer
        @param nBytes the number of bytes to print the data specs on
        @retval msg ints/floats/strings of the next nBytes
          (handles poorly sized nBytes; uncrashable :) )
        @note this the BEST function when adding new cards/tables/debugging
        """
        data = self.op2.read(nBytes)
        msg = self.print_block(data)
        self.op2.seek(self.n)
        return msg

    def print_section2(self, nBytes, endian):
        """
        Prints data, but doesn't move the cursor
        @param self the object pointer
        @param nBytes the number of bytes to print the data specs on
        @retval msg ints/floats/strings of the next nBytes
         (handles poorly sized nBytes; uncrashable :) )
        @note this the BEST function when adding new cards/tables/debugging
        """
        data = self.op2.read(nBytes)
        msg = self.print_block2(data, endian)
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

    def get_table_code(self, expected=None, debug=True):
        tablecode = self.read_header(expected, debug)
        return tablecode

    def get_marker(self, expected=None, debug=True):
        tablecode = self.read_header(expected, debug)
        return tablecode

    def read_marker(self, expected=None):
        return self.get_marker(expected)

    def read_markers(self, markers, table_name=None, debug=False,
                     printErrorOnFailure=True):
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
            table_code = self.read_header(marker, debug)
            #if table_code==2:
            #    table_code = self.read_header(marker,debug)

            #print("table_code=",table_code)
            if table_code is None:
                return
            if marker != table_code:
                msg = ''
                if printErrorOnFailure:
                    msg = '\nmarkers=%s foundMarkers=%s\n' % (markers,
                                                              foundMarkers)
                    msg += 'table_name=%s found=%s expected=%s leftover=%s' % (table_name, table_code, marker, self.print_section(40))
                    #print(msg)
                raise SyntaxError("Invalid Markers: %s" % msg)
            foundMarkers.append(marker)

        msg = ''
        for i in markers:
            msg += '[4,' + str(i) + ',4] + '
        if self.make_op2_debug:
            self.op2Debug.write(msg[:-3] + '\n')
        if debug:
            self.log.debug("@markers = %s" % (markers))
            self.log.debug("")

    def getNMarkers(self, nMarkers, rewind=False):
        """gets the next N markers, verifies they're correct"""
        markers = []
        for iMarker in xrange(nMarkers):
            table_code = self.read_header(None)
            markers.append(table_code)

        if rewind:
            self.n -= 12 * nMarkers
            self.op2.seek(self.n)

        return markers

    def is_table_done(self, expectedMarkers):
        markers = self.getNMarkers(len(expectedMarkers), rewind=True)
        #print "get_markers = ",markers

        if markers == [-1, 7]:
            return True
        elif markers == expectedMarkers:
            #sys.exit(expectedMarkers)
            return False
        else:
            msg = ('this should never happen...invalid markers...'
                   'expected=%s markers=%s' % (expectedMarkers, markers))
            raise RuntimeError(msg)

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

    def read_block(self):
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

    def read_full_block(self):
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

    def read_full_int_block(self):
        """
        Reads a fortran formatted block
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

    def read_string_block(self, debug=True):
        """
        reads a fortran formatted block
        assumes that the data is made up of characters only
        """
        data = self.read_block()
        nLetters = len(data)
        iFormat = str(nLetters) + 's'
        iFormat = bytes(iFormat)
        word, = unpack(iFormat, data)

        #print "word = |%s|" % (word)
        #print "nLetters=%s word=|%s|" % (nLetters,word)
        if debug and self.make_op2_debug:
            self.op2Debug.write('|%s|\n' % (str(word)))
        return word

    def read_int_block(self):
        """
        Reads a fortran formatted block
        assumes that the data is made up of integers only
        """
        data = self.read_block()
        nInts = len(data) // 4
        iFormat = str(nInts) + 'i'
        iFormat = bytes(iFormat)
        ints = unpack(iFormat, data)
        return ints

    def read_float_block(self):
        """
        Reads a fortran formatted block
        assumes that the data is made up of floats only
        """
        data = self.read_block()
        nFloats = len(data) // 4
        iFormat = str(nFloats) + 'f'
        iFormat = bytes(iFormat)
        floats = unpack(iFormat, data)
        return floats

    def read_double_block(self):
        """
        Reads a fortran formatted block
        assumes that the data is made up of doubles only
        """
        data = self.read_block()
        nDoubles = len(data) // 8
        iFormat = str(nDoubles) + 'd'
        iFormat = bytes(iFormat)
        doubles = unpack(iFormat, data)
        return doubles

    def rewind(self, n):
        """
        Rewinds the file nBytes
        @warning
            doesnt support a full rewind, only a partial
        """
        self.n -= n
        self.op2.seek(self.n)

    def read_table_name(self, rewind=True, debug=True, stopOnFailure=True):
        """
        Peeks into a table to check it's name
        """
        #debug = True
        if rewind:
            debug = False
        n = self.n
        try:
            #print ""
            self.read_markers([0, 2], debug)
            word = self.read_string_block(debug)
            #print("*word = |%r|" % word)

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
            table_name = word.strip()
            return table_name.decode('utf-8')
        except:
            if rewind and not stopOnFailure:
                self.n = n
                self.op2.seek(n)
                return
            raise

    def skip_next_table(self, buffer_size=10000):
        """
        skips a table
        @todo fix bugs
        """
        table_name = self.read_table_name(rewind=False)  # GEOM1
        self.table_init(table_name)
        self.log.debug("skippingTable |%s|" % table_name)
        self.log.debug("self.n = %s" % self.n)

        self.read_markers([-1, 7], table_name)

        dataPack = (4, 1, 4, 4, 0, 4, 4, 0, 4)  # marks the end of the table
        binaryData = pack('9i', *dataPack)
        #[1,0,0] marks the end of the table

        i = 0
        error = 80
        n = self.n
        endIndex = -1
        data = "dummy"
        while endIndex == -1 and len(data) > 0:
            data = self.op2.read(buffer_size + error)
            endIndex = data.find(binaryData)

            self.op2.seek(n + i * buffer_size)
            i += 1

        #print "i = ",i
        assert endIndex > 0, 'couldnt find the end of the table'

        # 36 so it gets to the end of the table markersNext=[0] or [2]
        self.n = self.n + (i - 1) * buffer_size + endIndex
        n = self.n
        self.n += 36  # TODO sometimes this is needed
        #ints = self.read_int_block()
        #print "*?*ints = ",ints
        #if len(ints)==0:
        #    pass

        self.op2.seek(self.n)
        self.log.debug("self.op2.tell() = %s" % self.op2.tell())
        #self.print_section(200)
        marker = self.get_marker()
        #print "marker = ",marker
        if marker == 2:
            isAnotherTable = True
        else:  # marker=0
            isAnotherTable = False
        isAnotherTable = True

        #print "isAnotherTable = ",isAnotherTable
        self.n -= 24  # subtract off the header [0,2] or [0,0]
        self.op2.seek(self.n)
        self.log.debug("self.n = %s" % (self.n))
        self.log.debug("---table %s is skipped---" % table_name)

        return isAnotherTable

    def has_more_tables(self):
        #print self.print_section(120)
        try:
            marker1 = self.get_marker('[4,0,4]')
            marker2 = self.get_marker('[4,0,4] or [4,2,4]')

            marker = [marker1, marker2]
            #print "marker = ",marker
            if marker == [0, 2]:
                isAnotherTable = True
            else:  # marker=0
                isAnotherTable = False

            #print "isAnotherTable = ",isAnotherTable
            self.n -= 24  # subtract off the header [0,2] or [0,0]
            self.op2.seek(self.n)
        except IndexError:
            isAnotherTable = False
        return isAnotherTable
