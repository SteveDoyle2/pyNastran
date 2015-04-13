from six.moves import range
from struct import unpack

from pyNastran.bdf.cards.aero import GUST
from pyNastran.bdf.cards.bdf_tables import (TABLED1, TABLED2, TABLED3, TABLEM1,
                                        TABLEM2, TABLEM3, TABLEM4)


class DIT(object):

    def _read_dit_4(self, data):
        return self._read_geom_4(self._dit_map, data)

    def __init__(self):
        self._dit_map = {
            (1005, 10, 174): self.readGust,     # record 1
            (1105, 11, 133): self.readTableD1,  # record 4
            (1205, 12, 134): self.readTableD2,  # record 5
            (1305, 13, 140): self.readTableD3,  # record 6

            #(105,1,93): self.readTableM1, # record 9
            (205, 2, 94): self.readTableM2,  # record 10
            (305, 3, 95): self.readTableM3,  # record 11
            #(405,4,96): self.readTableM4, # record 12

            (15, 21, 162): self._readFake,
            (56, 26, 303): self._readFake,
            (3105, 31, 97): self._readFake,  # record 13 - TABLES1
        }

    def readGust(self, data, n):
        """
        GUST(1005,10,174)    - the marker for Record 1
        """
        #print("reading GUST")
        nEntries = (len(data) - n) // 20  # 5*4
        for i in range(nEntries):
            eData = data[n:n + 20]
            out = unpack('ii3f', eData)
            (sid, dload, wg, x0, V) = out
            gust = GUST(None, out)
            self.add_GUST(gust)
            n += 20
        return n

#TABDMP1
#TABLE3D

    def readTableD1(self, data, n):
        """
        TABLED1(1105,11,133) - the marker for Record 4
        """
        #self.skippedCardsFile.write('skipping TABLED1 in DIT\n')
        return
        #print("reading TABLED1")
        func = TABLED1
        n = self.readTable1(func, data, n)
        return n

    def readTable1(self, func, data, n):
        #nEntries = len(data)//40 # 10*4
        n = 0
        ndata = len(data)
        while ndata - n >= 40:
            eData = data[n:n + 40]
            out = unpack('8iff', eData)
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

            dataIn += [x, y]
            table = func(None, out)
            self.add_table(table)
        return n

    def readTableD2(self, data, n):
        """
        TABLED2(1205,12,134) - the marker for Record 5
        """
        #print("reading TABLED2")
        func = TABLED2
        n = self.readTable2(func, data)
        return n

    def readTable2(self, func, data):
        n = 0
        return len(data)
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
            dataIn += [x, y]
            table = func(None, out)
            self.add_table(table)
        return len(data)

    def readTableD3(self, data, n):
        """
        TABLED3(1305,13,140) - the marker for Record 6
        """
        #print("reading TABLED3")
        func = TABLED3
        n = self.readTable3(func, data)
        return n

    def readTable3(self, func, data):
        n = 0
        ndata = len(data)
        while ndata - n >= 40:
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
            dataIn += [x, y]
            table = func(None, out)
            self.add_table(table)
        return len(data)

#TABLEDR

    def readTableM1(self, data, n):
        """
        TABLEM1(105,1,93) - the marker for Record 9
        """
        self.skippedCardsFile.write('skipping TABLEM1 in DIT\n')
        return
        #print("reading TABLED1")
        func = TABLEM1
        n = self.readTable1(func, data)
        return n

    def readTableM2(self, data, n):
        """
        TABLEM2(205,2,94) - the marker for Record 10
        """
        #print("reading TABLEM2")
        func = TABLEM2
        n = self.readTable2(func, data)
        return n

    def readTableM3(self, data, n):
        """
        TABLEM3(305,3,95) - the marker for Record 11
        """
        #print("reading TABLED3")
        func = TABLEM3
        n = self.readTable3(func, data)
        return n

    def readTableM4(self, data, n):
        """
        TABLEM4(405,4,96) - the marker for Record 12
        """
        #print("reading TABLED3")
        func = TABLEM4
        n = self.readTable4(func, data)
        return n

#TABLES1
#TABLEST
#TABRND1
#TABRNDG
