from struct import unpack
from six.moves import range

from pyNastran.bdf.cards.aero import GUST
from pyNastran.bdf.cards.bdf_tables import (TABLED1, TABLED2, TABLED3, TABLEM1,
                                            TABLEM2, TABLEM3, TABLEM4)
from pyNastran.op2.tables.geom.geom_common import GeomCommon


class DIT(GeomCommon):
    """defines methods for reading op2 tables"""

    def _read_dit_4(self, data, ndata):
        return self._read_geom_4(self._dit_map, data, ndata)

    def __init__(self):
        GeomCommon.__init__(self)
        self._dit_map = {
            (1005, 10, 174): self._read_gust,     # record 1
            (1105, 11, 133): self._read_tabled1,  # record 4
            (1205, 12, 134): self._read_tabled2,  # record 5
            (1305, 13, 140): self._read_tabled3,  # record 6

            #(105,1,93): self.readTableM1, # record 9
            (205, 2, 94): self._read_tablem2,  # record 10
            (305, 3, 95): self._read_tablem3,  # record 11
            #(405,4,96): self._read_tablem4, # record 12

            (15, 21, 162): self._read_fake,
            (56, 26, 303): self._read_fake,
            (3105, 31, 97): self._read_fake,  # record 13 - TABLES1
        }

    def _read_gust(self, data, n):
        """
        GUST(1005,10,174)    - the marker for Record 1
        """
        nentries = (len(data) - n) // 20  # 5*4
        for i in range(nentries):
            edata = data[n:n + 20]
            out = unpack('ii3f', edata)
            # (sid, dload, wg, x0, V) = out
            gust = GUST.add_op2_data(out)
            self.add_gust(gust)
            n += 20
        return n

#TABDMP1
#TABLE3D

    def _read_tabled1(self, data, n):
        """
        TABLED1(1105,11,133) - the marker for Record 4
        """
        #self.skippedCardsFile.write('skipping TABLED1 in DIT\n')
        return
        #print("reading TABLED1")
        cls = TABLED1
        n = self.read_table1(cls, data, n)
        return n

    def _read_table1(self, cls, data, n):
        #nentries = len(data)//40 # 10*4
        n = 0
        ndata = len(data)
        while ndata - n >= 40:
            edata = data[n:n + 40]
            out = unpack('8iff', edata)
            (sid, code_x, code_y, a, a, a, a, a, x, y) = out
            data_in = [sid, code_x, code_y]
            n += 40
            while 1:
                (xint, yint) = unpack('ii', data[n:n + 8])
                (x, y) = unpack('ff', data[n:n + 8])

                n += 8
                if [xint, yint] != [-1, -1]:
                    break
                else:
                    data_in += [x, y]

            data_in += [x, y]
            table = cls.add_op2_data(out)
            self.add_table(table)
        return n

    def _read_tabled2(self, data, n):
        """
        TABLED2(1205,12,134) - the marker for Record 5
        """
        cls = TABLED2
        n = self._read_table2(cls, data)
        return n

    def _read_table2(self, cls, data):
        n = 0
        return len(data)
        while len(data) >= 40:
            edata = data[n:n + 40]
            out = unpack('ifiiiiiiff', edata)
            (sid, x1, a, a, a, a, a, a, x, y) = out
            data_in = [sid, x1]
            n += 40
            while 1:
                (xint, yint) = unpack('ii', data[n:n + 8])
                (x, y) = unpack('ff', data[n:n + 8])

                n += 8
                if [xint, yint] != [-1, -1]:
                    break
                else:
                    data_in += [x, y]
            data_in += [x, y]
            table = cls.add_op2_data(out)
            self.add_table(table)
        return len(data)

    def _read_tabled3(self, data, n):
        """
        TABLED3(1305,13,140) - the marker for Record 6
        """
        cls = TABLED3
        n = self._read_table3(cls, data)
        return n

    def _read_table3(self, cls, data):
        n = 0
        ndata = len(data)
        while ndata - n >= 40:
            edata = data[n:n + 40]
            out = unpack('iffiiiiiff', edata)
            (sid, x1, x2, a, a, a, a, a, x, y) = out
            data_in = [sid, x1, x2]
            n += 40
            while 1:
                (xint, yint) = unpack('ii', data[n:n + 8])
                (x, y) = unpack('ff', data[n:n + 8])

                n += 8
                if [xint, yint] != [-1, -1]:
                    break
                else:
                    data_in += [x, y]
            data_in += [x, y]
            table = cls.add_op2_data(out)
            self.add_table(table)
        return len(data)

#TABLEDR

    def _read_tablem1(self, data, n):
        """
        TABLEM1(105,1,93) - the marker for Record 9
        """
        self.skippedCardsFile.write('skipping TABLEM1 in DIT\n')
        return
        cls = TABLEM1
        n = self._read_table1(cls, data)
        return n

    def _read_tablem2(self, data, n):
        """
        TABLEM2(205,2,94) - the marker for Record 10
        """
        cls = TABLEM2
        n = self._read_table2(cls, data)
        return n

    def _read_tablem3(self, data, n):
        """
        TABLEM3(305,3,95) - the marker for Record 11
        """
        cls = TABLEM3
        n = self._read_table3(cls, data)
        return n

    def _read_tablem4(self, data, n):
        """
        TABLEM4(405,4,96) - the marker for Record 12
        """
        cls = TABLEM4
        n = self._read_table4(cls, data)
        return n

#TABLES1
#TABLEST
#TABRND1
#TABRNDG
