import sys
from struct import unpack

from pyNastran.op2.tables.geom.geom1 import Geometry1
from pyNastran.op2.tables.geom.geom2 import Geometry2
from pyNastran.op2.tables.geom.geom3 import Geometry3
from pyNastran.op2.tables.geom.geom4 import Geometry4
from pyNastran.op2.tables.ept import EPT
from pyNastran.op2.tables.mpt import MPT
from pyNastran.op2.tables.geom.dynamics import DYNAMICS
from pyNastran.op2.tables.geom.dit import DIT


class GeomObj(object):
    def __init__(self):
        pass

    def geomFunc(self, data):
        pass


class GeometryTables(Geometry1, Geometry2, Geometry3, Geometry4, EPT, MPT, DIT,
                     DYNAMICS):

    def readRecordTable(self, expectedTableName):
        """
        .. note:: assumes self.iTableMap has already been set
        """
        if self.make_geom == False:
            self.iTableMap = {}

        table_name = self.read_table_name(rewind=False)  # GEOM1
        self._table_init(table_name)

        #print "*table_name = |%r|" %(table_name)

        self.read_markers([-1, 7])
        fields = self.read_int_block()
        #print "fields = ",fields

        self.read_markers([-2, 1, 0])  # 2
        buffer_words = self.get_marker()
        #print "buffer_words = ",buffer_words,buffer_words*4
        word = self.read_string_block()

        iTable = -3
        while 1:  # TODO could this cause an infinite loop...i dont this so...
            (table_name, isNextTable, isNextSubTable,
                isFileDone) = self.readGeomSubTable(iTable)

            if self.checkForNextTable() or isFileDone:
                #sys.exit('end of geom1')
                return
            iTable -= 1
        sys.exit('end of %s-this should never happen...' % (expectedTableName))

    def checkForNextTable(self):
        foundTable = False
        #print "---checking---"
        word = self.read_table_name(rewind=True, debug=False,
                                    stopOnFailure=False)
        if word is not None:
            foundTable = True
        #print '---checked---'
        #print "geomWord = ",word
        return foundTable

    def checkFileDone(self, n):
        isFileDone = False
        #print "tell = ",self.op2.tell()

        nOld = self.op2.tell()
        try:
            #print self.print_section(60)
            self.read_markers([n, 1, 0])
            markerA = self.get_marker()
            markerB = self.get_marker()
            #print "markerA=%s markerB=%s" %(markerA,markerB)
            #self.read_markers([0,0])
            #print "subtable :) = ",foundSubTable
            if [markerA, markerB] == [0, 0]:
                isFileDone = True
        except:
            pass
        self.n = nOld
        self.op2.seek(self.n)
        #print "isFileDone = ",isFileDone
        return isFileDone

    def checkForNextSubTable(self, n):
        foundSubTable = False
        #print "tell = ",self.op2.tell()

        try:
            nOld = self.op2.tell()
            self.read_markers([n, 1, 0])
            foundSubTable = True
            #print "subtable :) = ",foundSubTable
        except SyntaxError:
            foundSubTable = False
        self.n = nOld
        self.op2.seek(self.n)

        return foundSubTable

    def readGeomSubTable(self, iTable):
        i = 0
        isNextTable = False
        isNextSubTable = False
        self.read_markers([iTable, 1, 0])
        #print self.iTableMap

        table_name = self.read_table_name(rewind=True, stopOnFailure=False)
        if table_name:
            #print "**table_name = |%r|" %(table_name)
            return table_name, isNextTable, isNextSubTable, False

        data = b''
        isTableActive = False
        while isNextSubTable == False and isNextTable == False:
            #print self.print_section(200)
            marker = self.get_marker()
            #print "marker = ",marker
            if marker < 0:
                msg = 'marker is less than 0...'
                raise Exception(msg)
            data += self.read_block()
            if not isTableActive:
                tableType = unpack('iii', data[:12])
                data = data[12:]

            #print "iTable=%s lenGeomData=%s" %(iTable,len(data))

            if tableType in self.iTableMap:
                #print "reading  iTable=%-3s with tableType=%s" %(iTable,tableType)
                self.iTableMap[tableType](data)
            else:
                if not(tableType[0] == tableType[1] == tableType[2]):
                    msg = "skipping %s iTable=%-3s with tableType=%s" % (
                        self.table_name, iTable, tableType)
                    #self.skippedCardsFile.write(msg+'\n')
                    if self.make_geom:
                        self.log.debug(msg)
            #self.op2Debug.write('ints = %s\n' %(str(ints)))

            isNextTable = self.checkForNextTable()
            isNextSubTable = self.checkForNextSubTable(iTable - 1)
            isFileDone = self.checkFileDone(iTable - 1)
            if isFileDone:
                isNextTable = True

            #print "i=%s tell=%s isNextTable=%s isNextSubTable=%s" %(i,self.op2.tell(),isNextTable,isNextSubTable)
            #if i==13:
            #    sys.exit('stopA')
            i += 1
            isTableActive = True

        #print "exiting the geom sub table"
        return (table_name, isNextTable, isNextSubTable, isFileDone)

    def readTable_PCOMPTS(self):
        #self.iTableMap = {
        #                 }
        #self.readRecordTable('PCOMPTS')

        table_name = self.read_table_name(rewind=False)  # PCOMP
        self._table_init(table_name)
        self.read_markers([-1, 7])
        ints = self.read_int_block()  # ??? ints
        #print(ints)
        #data = self.read_block()
        #print self.print_block(data)
        #print("fields = ",fields)

        #-------------------------------------------
        self.read_markers([-2, 1, 0])  # 2
        self.read_markers([2])  # 2
        strings = self.read_string_block()  # IPCOMPT
        #print(strings)

        #-------------------------------------------
        #print "3"
        iTable = -3
        while 1:
            self.read_markers([iTable, 1, 0])  # 3
            n = self.op2.tell()
            try:
                buffer_words = self.get_marker()
                if buffer_words == 0:
                    self.goto(n)
                    #print "returning from table=-3"
                    return
                elif buffer_words < 0:
                    self.goto(n)
                else:
                    #print "buffer_words = ",buffer_words,buffer_words*4
                    data = self.get_data(4)
                    buffer_size, = unpack('i', data)

                    #print "buffer_size = ",buffer_size
                    data = self.get_data(buffer_words * 4)
                    data = self.get_data(4)
            except:
                raise RuntimeError('error in iTable=% of %s...' %
                                   (self.iTable, self.table_name))
            iTable -= 1

    def readTable_SDF(self):
        table_name = self.read_table_name(rewind=False)  # SDF
        self._table_init(table_name)
        self.read_markers([-1, 7])
        ints = self.read_int_block()  # ??? ints
        #print ints

        #-------------------------------------------

        iTable = -2

        #print "iTable = ",iTable
        self.read_markers([iTable, 1, 0])  # 2
        buffer_words = self.get_marker()
        #print "buffer_words = ",buffer_words
        data = self.get_data(4)
        buffer_size, = unpack('i', data)
        data = self.get_data(buffer_words * 4)
        #print self.print_block(data)
        data = self.get_data(4)
        iTable -= 1

        self.read_markers([iTable, 1, 1])  # 3
        buffer_words = self.get_marker()  # 12
        data = self.get_data(4)
        buffer_size, = unpack('i', data)  # 52
        data = self.get_data(buffer_size)
        #print self.print_block(data)
        data = self.get_data(4)
        iTable -= 1

        self.read_markers([iTable, 1, 0])  # 4

        #-------------------------------------------
        #print self.print_section(240)
        #sys.exit('SDF...')
    def readTable_CASECC(self):
        table_name = self.read_table_name(rewind=False)  # CASECC
        #print '*table_name = ',table_name
        self._table_init(table_name)
        self.read_markers([-1, 7])
        data = self.get_data(4)
        buffer_size, = unpack('i', data)
        #print "buffer_size = ",buffer_size
        data = self.get_data(buffer_size)
        #print self.print_block(data)
        data = self.get_data(4)
        #print "---------------"

        self.read_markers([-2, 1, 0])
        buffer_words = self.get_marker()
        data = self.get_data(4)
        buffer_size, = unpack('i', data)
        #print "buffer_size = ",buffer_size
        data = self.get_data(buffer_size)
        #print self.print_block(data)
        data = self.get_data(4)

        #data = self.read_block()
        #print self.print_block(data)
        #print "---------------"
        self.read_markers([-3, 1, 0])
        buffer_words = self.get_marker()
        data = self.get_data(4)
        buffer_size, = unpack('i', data)
        #print "buffer_words = ",buffer_words
        #print "buffer_size = ",buffer_size
        data = self.get_data(buffer_size)
        #print self.print_block(data)
        data = self.get_data(4)

        #print "---------------"
        self.read_markers([-4, 1, 0])

        #print(self.print_section(240))
        sys.exit('CASECC...')

    def readTable_OMM2(self):
        #-------------------------------------------
        table_name = self.read_table_name(rewind=False)  # PCOMP
        self._table_init(table_name)
        self.read_markers([-1, 7])
        ints = self.read_int_block()  # ??? ints
        #print ints

        #-------------------------------------------
        iTable = -2
        while 1:
            self.read_markers([iTable, 1, 0])  # 2
            try:
                n = self.op2.tell()
                #print "iTable = ",iTable
                buffer_words = self.get_marker()
                #print "buffer_words = ",buffer_words
                data = self.get_data(4)
                buffer_size, = unpack('i', data)
                #print "buffer_size = ",buffer_size
                data = self.get_data(buffer_words * 4)
                #print self.print_block(data)
                data = self.get_data(4)
                iTable -= 1
            except:
                self.goto(n)
                break

        #-------------------------------------------
        #print self.print_section(400)
        #sys.exit('OMM2...stop...')

    def readTable_DUMMY_GEOM(self, table_name):
        self.iTableMap = {}
        self.readRecordTable(table_name)
