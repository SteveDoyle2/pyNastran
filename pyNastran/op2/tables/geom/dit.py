import sys
from struct import unpack

from pyNastran.bdf.cards.aero import GUST
from pyNastran.bdf.cards.tables import (TABLED1, TABLED2, TABLED3, TABLEM1,
                                        TABLEM2, TABLEM3, TABLEM4)


class DIT(object):
    def readTable_DIT(self):
        self.iTableMap = {
            (1005, 10, 174): self.readGust,     # record 1
            (1105, 11, 133): self.readTableD1,  # record 4
            (1205, 12, 134): self.readTableD2,  # record 5
            (1305, 13, 140): self.readTableD3,  # record 6

            #(105,1,93): self.readTableM1, # record 9
            (205, 2, 94): self.readTableM2,  # record 10
            (305, 3, 95): self.readTableM3,  # record 11
            #(405,4,96): self.readTableM4, # record 12

            (15, 21, 162): self.readFake,
            (56, 26, 303): self.readFake,
            (3105, 31, 97): self.readFake,  # record 13 - TABLES1
        }
        self.readRecordTable('DIT')

    def readGust(self, data):
        """
        GUST(1005,10,174)    - the marker for Record 1
        """
        #print "reading GUST"
        n = 0
        nEntries = len(data) // 20  # 5*4
        for i in xrange(nEntries):
            eData = data[n:n + 20]
            out = unpack('iifff', eData)
            (sid, dload, wg, x0, V) = out
            gust = GUST(None, out)
            self.addGUST(gust)
            n += 20
        ###
        data = data[n:]

#TABDMP1
#TABLE3D

    def readTableD1(self, data):
        """
        TABLED1(1105,11,133) - the marker for Record 4
        """
        self.skippedCardsFile.write('skipping TABLED1 in DIT\n')
        return
        #print "reading TABLED1"
        func = TABLED1
        self.readTable1(func, data)

    def readTable1(self, func, data):
        n = 0
        #nEntries = len(data)//40 # 10*4
        while len(data) >= 40:
            eData = data[n:n + 40]
            out = unpack('iiiiiiiiff', eData)
            (sid, codeX, codeY, a, a, a, a, a, x, y) = out
            dataIn = [sid, codeX, codeY]
            n += 40
            while 1:
                (xInt, yInt) = unpack('ii', data[n:n + 8])
                (x, y) = unpack('ff', data[n:n + 8])

                n += 8
                if [xInt, yInt] != [-1, -1]:
                    break
                else:
                    dataIn += [x, y]
                ###
            dataIn += [x, y]
            table = func(None, out)
            self.addTable(table)
        ###
        data = data[n:]

    def readTableD2(self, data):
        """
        TABLED2(1205,12,134) - the marker for Record 5
        """
        #print "reading TABLED2"
        func = TABLED2
        self.readTable2(func, data)

    def readTable2(self, func, data):
        n = 0
        while len(data) >= 40:
            eData = data[n:n + 40]
            out = unpack('ifiiiiiiff', eData)
            (sid, x1, a, a, a, a, a, a, x, y) = out
            dataIn = [sid, x1]
            n += 40
            while 1:
                (xInt, yInt) = unpack('ii', data[n:n + 8])
                (x, y) = unpack('ff', data[n:n + 8])

                n += 8
                if [xInt, yInt] != [-1, -1]:
                    break
                else:
                    dataIn += [x, y]
                ###
            dataIn += [x, y]
            table = func(None, out)
            self.addTable(table)
        ###
        data = data[n:]

    def readTableD3(self, data):
        """
        TABLED3(1305,13,140) - the marker for Record 6
        """
        #print "reading TABLED3"
        func = TABLED3
        self.readTable3(func, data)

    def readTable3(self, func, data):
        n = 0
        while len(data) >= 40:
            eData = data[n:n + 40]
            out = unpack('iffiiiiiff', eData)
            (sid, x1, x2, a, a, a, a, a, x, y) = out
            dataIn = [sid, x1, x2]
            n += 40
            while 1:
                (xInt, yInt) = unpack('ii', data[n:n + 8])
                (x, y) = unpack('ff', data[n:n + 8])

                n += 8
                if [xInt, yInt] != [-1, -1]:
                    break
                else:
                    dataIn += [x, y]
                ###
            dataIn += [x, y]
            table = func(None, out)
            self.addTable(table)
        ###
        data = data[n:]

#TABLEDR

    def readTableM1(self, data):
        """
        TABLEM1(105,1,93) - the marker for Record 9
        """
        self.skippedCardsFile.write('skipping TABLEM1 in DIT\n')
        return
        #print "reading TABLED1"
        func = TABLEM1
        self.readTable1(func, data)

    def readTableM2(self, data):
        """
        TABLEM2(205,2,94) - the marker for Record 10
        """
        #print "reading TABLEM2"
        func = TABLEM2
        self.readTable2(func, data)

    def readTableM3(self, data):
        """
        TABLEM3(305,3,95) - the marker for Record 11
        """
        #print "reading TABLED3"
        func = TABLEM3
        self.readTable3(func, data)

    def readTableM4(self, data):
        """
        TABLEM4(405,4,96) - the marker for Record 12
        """
        #print "reading TABLED3"
        func = TABLEM4
        self.readTable4(func, data)


#TABLES1
#TABLEST
#TABRND1
#TABRNDG
