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
        @note assumes self.iTableMap has already been set
        """
        if self.makeGeom == False:
            self.iTableMap = {}

        tableName = self.readTableName(rewind=False)  # GEOM1
        self.tableInit(tableName)
        #print "*tableName = |%r|" %(tableName)

        self.readMarkers([-1, 7])
        fields = self.readIntBlock()
        #print "fields = ",fields

        self.readMarkers([-2, 1, 0])  # 2
        bufferWords = self.getMarker()
        #print "bufferWords = ",bufferWords,bufferWords*4
        word = self.readStringBlock()

        iTable = -3
        while 1:  ## @todo could this cause an infinite loop...i dont this so...
            (tableName, isNextTable, isNextSubTable,
                isFileDone) = self.readGeomSubTable(iTable)

            if self.checkForNextTable() or isFileDone:
                #sys.exit('end of geom1')
                return
            iTable -= 1
        ###
        sys.exit('end of %s-this should never happen...' % (expectedTableName))

    def checkForNextTable(self):
        foundTable = False
        #print "---checking---"
        word = self.readTableName(
            rewind=True, debug=False, stopOnFailure=False)
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
            #print self.printSection(60)
            self.readMarkers([n, 1, 0])
            markerA = self.getMarker()
            markerB = self.getMarker()
            #print "markerA=%s markerB=%s" %(markerA,markerB)
            #self.readMarkers([0,0])
            #print "subtable :) = ",foundSubTable
            if [markerA, markerB] == [0, 0]:
                isFileDone = True
            ###
        except:
            pass
        ###
        self.n = nOld
        self.op2.seek(self.n)
        #print "isFileDone = ",isFileDone
        return isFileDone

    def checkForNextSubTable(self, n):
        foundSubTable = False
        #print "tell = ",self.op2.tell()

        try:
            nOld = self.op2.tell()
            self.readMarkers([n, 1, 0])
            foundSubTable = True
            #print "subtable :) = ",foundSubTable
        except SyntaxError:
            foundSubTable = False
        ###
        self.n = nOld
        self.op2.seek(self.n)

        return foundSubTable

    def readGeomSubTable(self, iTable):
        i = 0
        isNextTable = False
        isNextSubTable = False
        self.readMarkers([iTable, 1, 0])
        #print self.iTableMap

        tableName = self.readTableName(rewind=True, stopOnFailure=False)
        if tableName:
            #print "**tableName = |%r|" %(tableName)
            return tableName, isNextTable, isNextSubTable, False

        data = b''
        isTableActive = False
        while isNextSubTable == False and isNextTable == False:
            #print self.printSection(200)
            marker = self.getMarker()
            #print "marker = ",marker
            if marker < 0:
                msg = 'marker is less than 0...'
                raise Exception(msg)
            data += self.readBlock()
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
                        self.tableName, iTable, tableType)
                    #self.skippedCardsFile.write(msg+'\n')
                    if self.makeGeom:
                        self.log.debug(msg)
                    ###
                ###
            ###
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
        ### while

        #print "exiting the geom sub table"
        return (tableName, isNextTable, isNextSubTable, isFileDone)

    def readTable_PCOMPTS(self):
        #self.iTableMap = {
        #                 }
        #self.readRecordTable('PCOMPTS')

        tableName = self.readTableName(rewind=False)  # PCOMP
        self.tableInit(tableName)
        self.readMarkers([-1, 7])
        ints = self.readIntBlock()  # ??? ints
        #print(ints)
        #data = self.readBlock()
        #print self.printBlock(data)
        #print("fields = ",fields)

        #-------------------------------------------
        self.readMarkers([-2, 1, 0])  # 2
        self.readMarkers([2])  # 2
        strings = self.readStringBlock()  # IPCOMPT
        #print(strings)

        #-------------------------------------------
        #print "3"
        iTable = -3
        while 1:
            self.readMarkers([iTable, 1, 0])  # 3
            n = self.op2.tell()
            try:
                bufferWords = self.getMarker()
                if bufferWords == 0:
                    self.goto(n)
                    #print "returning from table=-3"
                    return
                elif bufferWords < 0:
                    self.goto(n)
                else:
                    #print "bufferWords = ",bufferWords,bufferWords*4
                    data = self.getData(4)
                    bufferSize, = unpack('i', data)

                    #print "bufferSize = ",bufferSize
                    data = self.getData(bufferWords * 4)
                    data = self.getData(4)
                ###
            except:
                raise RuntimeError('error in iTable=% of %s...' %
                                   (self.iTable, self.tableName))
            ###
            iTable -= 1
            ###
        ###

    def readTable_SDF(self):
        tableName = self.readTableName(rewind=False)  # SDF
        self.tableInit(tableName)
        self.readMarkers([-1, 7])
        ints = self.readIntBlock()  # ??? ints
        #print ints

        #-------------------------------------------

        iTable = -2

        #print "iTable = ",iTable
        self.readMarkers([iTable, 1, 0])  # 2
        bufferWords = self.getMarker()
        #print "bufferWords = ",bufferWords
        data = self.getData(4)
        bufferSize, = unpack('i', data)
        data = self.getData(bufferWords * 4)
        #print self.printBlock(data)
        data = self.getData(4)
        iTable -= 1

        self.readMarkers([iTable, 1, 1])  # 3
        bufferWords = self.getMarker()  # 12
        data = self.getData(4)
        bufferSize, = unpack('i', data)  # 52
        data = self.getData(bufferSize)
        #print self.printBlock(data)
        data = self.getData(4)
        iTable -= 1

        self.readMarkers([iTable, 1, 0])  # 4

        #-------------------------------------------
        #print self.printSection(240)
        #sys.exit('SDF...')
    def readTable_CASECC(self):
        tableName = self.readTableName(rewind=False)  # CASECC
        #print '*tableName = ',tableName
        self.tableInit(tableName)
        self.readMarkers([-1, 7])
        data = self.getData(4)
        bufferSize, = unpack('i', data)
        #print "bufferSize = ",bufferSize
        data = self.getData(bufferSize)
        #print self.printBlock(data)
        data = self.getData(4)
        #print "---------------"

        self.readMarkers([-2, 1, 0])
        bufferWords = self.getMarker()
        data = self.getData(4)
        bufferSize, = unpack('i', data)
        #print "bufferSize = ",bufferSize
        data = self.getData(bufferSize)
        #print self.printBlock(data)
        data = self.getData(4)

        #data = self.readBlock()
        #print self.printBlock(data)
        #print "---------------"
        self.readMarkers([-3, 1, 0])
        bufferWords = self.getMarker()
        data = self.getData(4)
        bufferSize, = unpack('i', data)
        #print "bufferWords = ",bufferWords
        #print "bufferSize = ",bufferSize
        data = self.getData(bufferSize)
        #print self.printBlock(data)
        data = self.getData(4)

        #print "---------------"
        self.readMarkers([-4, 1, 0])

        #print(self.printSection(240))
        sys.exit('CASECC...')

    def readTable_OMM2(self):
        #-------------------------------------------
        tableName = self.readTableName(rewind=False)  # PCOMP
        self.tableInit(tableName)
        self.readMarkers([-1, 7])
        ints = self.readIntBlock()  # ??? ints
        #print ints

        #-------------------------------------------
        iTable = -2
        while 1:
            self.readMarkers([iTable, 1, 0])  # 2
            try:
                n = self.op2.tell()
                #print "iTable = ",iTable
                bufferWords = self.getMarker()
                #print "bufferWords = ",bufferWords
                data = self.getData(4)
                bufferSize, = unpack('i', data)
                #print "bufferSize = ",bufferSize
                data = self.getData(bufferWords * 4)
                #print self.printBlock(data)
                data = self.getData(4)
                iTable -= 1
            except:
                self.goto(n)
                break
            ###
        ###

        #-------------------------------------------
        #print self.printSection(400)
        #sys.exit('OMM2...stop...')

    def readTable_DUMMY_GEOM(self, tableName):
        self.iTableMap = {}
        self.readRecordTable(tableName)
