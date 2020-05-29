# pylint: disable=C0103
"""
defines readers for BDF objects in the OP2 DIT/DITS table
"""
from struct import Struct, error as struct_error
import numpy as np

from pyNastran.bdf.cards.aero.dynamic_loads import GUST
from pyNastran.bdf.cards.bdf_tables import (TABLED1, TABLED2, TABLED3, TABLED4,
                                            TABLEM1, TABLEM2, TABLEM3, TABLEM4,
                                            TABRND1, TABDMP1, TABLES1)
from pyNastran.op2.tables.geom.geom_common import GeomCommon
from pyNastran.op2.op2_interface.op2_reader import mapfmt, reshape_bytes_block


class DIT(GeomCommon):
    """defines methods for reading op2 tables"""

    def _read_dit_4(self, data, ndata):
        return self._read_geom_4(self._dit_map, data, ndata)

    def __init__(self):
        GeomCommon.__init__(self)
        self._dit_map = {
            (1005, 10, 174): ['GUST', self._read_gust],     # record 1
            (1105, 11, 133): ['TABLED1', self._read_tabled1],  # record 4
            (1205, 12, 134): ['TABLED2', self._read_tabled2],  # record 5
            (1305, 13, 140): ['TABLED3', self._read_tabled3],  # record 6
            (1405, 14, 141) : ['TABLED4', self._read_tabled4],  # record 7-MSC
            #(4201, 42, 648) : ['TABLEDR', self._read_tabledr], # record 8-MSC

            (105, 1, 93): ['TABLEM1', self._read_tablem1], # record 9
            (205, 2, 94): ['TABLEM2', self._read_tablem2],  # record 10
            (305, 3, 95): ['TABLEM3', self._read_tablem3],  # record 11
            (405, 4, 96): ['TABLEM4', self._read_tablem4], # record 12

            (55, 25, 191) : ['TABRND1', self._read_tabrnd1],
            (15, 21, 162): ['TABDMP1', self._read_tabdmp1],   # NX
            (56, 26, 303): ['TABRNDG', self._read_tabrndg],   # NX
            (3105, 31, 97): ['TABLES1', self._read_tables1],  # record 13 - TABLES1 (NX)
            (4000, 40, 460) : ['TABLE3D', self._read_fake],

            # F:\work\pyNastran\examples\Dropbox\move_tpl\htab11.op2
            (14705, 147, 618) : ['TABLEHT', self._read_fake],
            (14605, 146, 617) : ['TABLEH1', self._read_fake],

            # F:\work\pyNastran\examples\Dropbox\move_tpl\n10640b.op2
            (1905, 19, 178) : ['TABLEST', self._read_fake],
        }

    def _read_tabdmp1(self, data: bytes, n: int) -> int:
        """
        TABDMP1(15, 21, 162)

        1 ID    I  Table identification number
        9 F     RS Natural frequency
        10 G    RS Damping
        Words 9 through 10 repeat until (-1,-1) occurs
        """
        #nfields = (ndata - n) // 4

        datan = data[n:]
        ints = np.frombuffer(datan, self.idtype8)
        floats = np.frombuffer(datan, self.fdtype8)
        iminus1_delta = get_iend_from_ints(ints)
        istart = 0
        nentries = 0
        for iend in iminus1_delta:
            #datai = data[n+istart*4 : n+iend*4]
            tid = ints[istart].copy()
            deltai = iend - istart - 8 # subtract 2 for sid, global scale
            assert deltai % 2 == 0, (self.show_data(data[n+istart*4 : n+iend*4], 'if'))

            xy = floats[istart+8:iend].reshape(deltai//2, 2).copy()
            x = xy[:, 0]
            y = xy[:, 1]
            table = TABDMP1(tid, x, y, Type='G')
            if tid in self.tables_sdamping:
                assert table == self.tables_sdamping[tid]
            else:
                self._add_table_sdamping_object(table)
            istart = iend + 2
            nentries += 1
        self.increase_card_count('TABDMP1', nentries)
        return len(data)

    def _read_tabrndg(self, data: bytes, n: int) -> int:
        """
        TABRNDG(56, 26, 303)
        Power spectral density for gust loads in aeroelastic analysis.

        1 ID        I   Table identification number
        2 TYPE      I   Power spectral density type
        3 LU        RS  Scale of turbulence divided by velocity
        4 WG        RS  Root-mean-square gust velocity
        5 UNDEF(4) none Not used
        Words 1 through 8 repeat until (-1,-1) occurs

        """
        ndata = len(data)# - n
        assert ndata == 52, ndata
        struct_2i2f4i = Struct('2i2f4i')
        #struct_ff = Struct('ff')
        #struct_2i = self.struct_2i
        while ndata - n >= 32:
            edata = data[n:n + 32]
            out = struct_2i2f4i.unpack(edata)
            (tid, table_type, lu, wg, unused_dunno_a,
             unused_dunno_b, unused_dunno_c, unused_dunno_d) = out
            if tid > 100000000:
                tid = -(tid - 100000000)
            n += 32
            self.add_tabrndg(tid, table_type, lu, wg, comment='')
            #nentries += 1
        #self.increase_card_count('TABRNDG', nentries)
        n += 8  #  for the (-1,-1)
        return n

    def _read_tables1(self, data: bytes, n: int) -> int:
        """TABLES1(3105, 31, 97)"""
        n = self._read_table1(TABLES1, self.tables, self._add_table_object, data, n, 'TABLES1',
                              add_codes=False)
        return n

    def _read_gust(self, data: bytes, n: int) -> int:
        """
        GUST(1005,10,174) - the marker for Record 1
        """
        nentries = (len(data) - n) // 20  # 5*4
        struct_2i3f = Struct('ii3f')
        for unused_i in range(nentries):
            edata = data[n:n + 20]
            out = struct_2i3f.unpack(edata)
            # (sid, dload, wg, x0, V) = out
            gust = GUST.add_op2_data(out)
            self._add_gust_object(gust)
            n += 20
        return n

#TABDMP1
#TABLE3D

    def _read_tabled1(self, data: bytes, n: int) -> int:
        """
        TABLED1(1105,11,133) - the marker for Record 4
        """
        n = self._read_table1(TABLED1, self.tables_d, self._add_tabled_object, data, n, 'TABLED1')
        return n

    def _read_table1(self, cls, slot, add_method, data, n, table_name, add_codes=True):
        nentries = 0
        ndata = len(data)
        if self.size == 4:
            struct_8i2f = Struct('8iff')
            struct_ff = Struct('ff')
            struct_2i = self.struct_2i
            ntotal1 = 40
            ntotal2 = 8
        else:
            struct_8i2f = Struct('8qdd')
            struct_ff = Struct('dd')
            struct_2i = self.struct_2q
            ntotal1 = 80
            ntotal2 = 16

        while ndata - n >= ntotal1:
            edata = data[n:n + ntotal1]
            out = struct_8i2f.unpack(edata)
            (tid, code_x, code_y, unused_a, unused_b, unused_c, unused_d, unused_e,
             x, y) = out
            if tid > 100000000:
                tid = -(tid - 100000000)
            if add_codes:
                data_in = [tid, code_x, code_y, x, y]
            else:
                data_in = [tid, x, y]

            n += ntotal1
            while 1:
                (xint, yint) = struct_2i.unpack(data[n:n + ntotal2])
                (x, y) = struct_ff.unpack(data[n:n + ntotal2])

                n += ntotal2
                if [xint, yint] == [-1, -1]:
                    break
                else:
                    data_in += [x, y]

            #print('data_in =', data_in)
            table = cls.add_op2_data(data_in)
            if tid in slot:
                assert table == slot[tid]
            else:
                add_method(table)
            nentries += 1
        self.increase_card_count(table_name, nentries)
        return n

    def _read_tabled2(self, data: bytes, n: int) -> int:
        """
        TABLED2(1205,12,134) - the marker for Record 5
        """
        n = self._read_table2(TABLED2, self.tables_d, self._add_tabled_object, data, n, 'TABLED2')
        return n

    def _read_table2(self, cls, slot, add_method, data, n, table_name):
        """
        1 ID    I  Table identification number
        2 X1    RS X-axis shift
        3 FLAG  I  Extrapolation on/off flag
        4 UNDEF I  None
        9  X RS X  value
        10 Y RS Y  value
        Words 9 through 10 repeat until (-1,-1) occurs
        """
        ndata = len(data)
        nentries = 0
        if self.size == 4:
            struct1 = Struct('ifiiiiiiff')
            struct_ff = Struct('ff')
            struct_2i = self.struct_2i
            ntotal1 = 40
            ntotal2 = 8
        else:
            struct1 = Struct('qdqqqqqqdd')
            struct_ff = Struct('dd')
            struct_2i = self.struct_2q
            ntotal1 = 80
            ntotal2 = 16
        while n < ndata:
            edata = data[n:n + ntotal1]
            out = struct1.unpack(edata)
            (tid, x1, unused_a, unused_b, unused_c, unused_d, unused_e, unused_f,
             x, y) = out
            data_in = [tid, x1, x, y]
            n += ntotal1
            while 1:
                (xint, yint) = struct_2i.unpack(data[n:n + ntotal2])
                (x, y) = struct_ff.unpack(data[n:n + ntotal2])

                n += ntotal2
                if [xint, yint] == [-1, -1]:
                    break
                else:
                    data_in += [x, y]
            table = cls.add_op2_data(data_in)
            add_method(table)
            nentries += 1
        self.increase_card_count(table_name, nentries)
        return n

    def _read_tabled3(self, data: bytes, n: int) -> int:
        """
        TABLED3(1305,13,140) - the marker for Record 6
        """
        n = self._read_table3(TABLED3, self.tables_d, self._add_tabled_object, data, n, 'TABLED3')
        return n

    def _read_tabled4(self, data: bytes, n: int) -> int:
        """
        TABLED4 - the marker for Record 7
        """
        n = self._read_table4(TABLED4, self.tables_d, self._add_tabled_object, data, n, 'TABLED4')
        return n

#TABLEDR

    def _read_tablem1(self, data: bytes, n: int) -> int:
        """
        TABLEM1(105,1,93) - the marker for Record 9
        """
        n = self._read_table1(TABLEM1, self.tables_m, self._add_tablem_object, data, n, 'TABLEM1')
        return n

    def _read_tablem2(self, data: bytes, n: int) -> int:
        """
        TABLEM2(205,2,94) - the marker for Record 10
        """
        n = self._read_table2(TABLEM2, self.tables_m, self._add_tablem_object, data, n, 'TABLEM2')
        return n

    def _read_tablem3(self, data: bytes, n: int) -> int:
        """
        TABLEM3(305,3,95) - the marker for Record 11
        """
        n = self._read_table3(TABLEM3, self.tables_m, self._add_tablem_object, data, n, 'TABLEM3')
        return n

    def _read_tablem4(self, data: bytes, n: int) -> int:
        """
        TABLEM4(405,4,96) - the marker for Record 12
        """
        n = self._read_table4(TABLEM4, self.tables_m, self._add_tablem_object, data, n, 'TABLEM4')
        return n

    def _read_table3(self, cls, slot, add_method, data, n, table_name):
        nentries = 0
        ndata = len(data)
        if self.size == 4:
            ntotal1 = 40
            ntotal2 = 8
            struct1 = Struct('iffiiiiiff')
            struct_2i = self.struct_2i
            struct_ff = Struct('ff')
        else:
            ntotal1 = 80
            ntotal2 = 16
            struct1 = Struct('qddqqqqqdd')
            struct_2i = self.struct_2q
            struct_ff = Struct('dd')

        while ndata - n >= ntotal1:
            edata = data[n:n + ntotal1]
            out = struct1.unpack(edata)
            (tid, x1, x2, unused_a, unused_b, unused_c, unused_d, unused_e,
             x, y) = out
            data_in = [tid, x1, x2, x, y]
            n += ntotal1
            while 1:
                (xint, yint) = struct_2i.unpack(data[n:n + ntotal2])
                (x, y) = struct_ff.unpack(data[n:n + ntotal2])

                n += ntotal2
                if [xint, yint] == [-1, -1]:
                    break
                else:
                    data_in += [x, y]
            table = cls.add_op2_data(data_in)
            add_method(table)
            nentries += 1
        self.increase_card_count(table_name, nentries)
        return n

    def _read_table4(self, cls, slot, add_method, data, n, table_name):
        """
        1 ID I Table identification number
        2 X1 RS X-axis shift
        3 X2 RS X-axis normalization
        4 X3 RS X value when x is less than X3
        5 X4 RS X value when x is greater than X4
        6 UNDEF(3 ) None
        9 A RS
        Word 9 repeats until End of Record (-1)
        """
        n0 = n
        nentries = 0
        ndata = len(data)
        size = self.size
        struct1 = Struct(mapfmt(self._endian + b'i 4f 3i f i', size))
        struct_i = self.struct_i if size == 4 else self.struct_q
        struct_f = Struct(self._endian + b'f') if size == 4 else Struct(self._endian + b'd')
        ntotal1 = 40 * self.factor
        ntotal2 = 36 * self.factor
        try:
            while ndata - n >= ntotal1:
                edata = data[n:n + ntotal1]
                out = struct1.unpack(edata)
                (tid, x1, x2, x3, x4, unused_a, unused_b, unused_c, x, test_minus1) = out
                data_in = [tid, x1, x2, x3, x4, x]
                n += ntotal2
                if test_minus1 == -1:
                    n += size
                else:
                    while 1:
                        xint, = struct_i.unpack(data[n:n + size])
                        x, = struct_f.unpack(data[n:n + size])

                        n += size
                        if xint == -1:
                            break
                        else:
                            data_in.append(x)
                table = cls.add_op2_data(data_in)
                add_method(table)
                nentries += 1
        except struct_error:
            self.log.error('failed parsing %s' % table_name)
            self.show_data(data[n0:], 'if')
            self.show_data(edata, 'if')
            #n = n0 + ndata
            raise
        self.increase_card_count(table_name, nentries)
        return n

#TABLEST
    def _read_tabrnd1(self, data: bytes, n: int) -> int:
        """
        TABRND1(55,25,191)

        1 ID    I Table identification number
        2 CODEX I Type of interpolation for the x-axis
        3 CODEY I Type of interpolation for the y-axis
        4 UNDEF(5) None
        9 F    RS Frequency
        10 G   RS Power spectral density
        Words 9 through 10 repeat until (-1,-1) occurs
        """
        #nfields = (ndata - n) // 4

        datan = data[n:]
        ints = np.frombuffer(datan, self.idtype).copy()
        floats = np.frombuffer(datan, self.fdtype)
        iminus1_delta = get_iend_from_ints(ints)
        istart = 0
        nentries = 0
        for iend in iminus1_delta:
            #datai = data[n+istart*4 : n+iend*4]
            tid = ints[istart]
            codex = ints[istart + 1]
            codey = ints[istart + 2]
            #print('  sid=%s global_scale=%s' % (sid, global_scale))
            deltai = iend - istart - 8 # subtract 2 for sid, global scale
            assert deltai % 2 == 0, (self.show_data(data[n+istart*4 : n+iend*4], 'if'))

            xy = floats[istart+8:iend].reshape(deltai//2, 2).copy()
            x = xy[:, 0]
            y = xy[:, 1]
            if codex == 0:
                xaxis = 'LINEAR'
            elif codex == 1:
                xaxis = 'LOG'
            else:
                raise NotImplementedError(codex) # LOG
            if codey == 0:
                yaxis = 'LINEAR'
            elif codey == 1:
                yaxis = 'LOG'
            else:
                raise NotImplementedError(codey) # LOG

            table = TABRND1(tid, x, y, xaxis=xaxis, yaxis=yaxis)
            self._add_random_table_object(table)
            istart = iend + 2
            nentries += 1
        self.increase_card_count('TABRND1', nentries)
        return len(data)

def get_iend_from_ints(ints):
    """
    istart = iend + 2
    for (-1, -1) flag to end a table
    """
    #debug = True
    iminus1 = np.where(ints == -1)[0]
    delta = iminus1[1:] - iminus1[:-1]
    #if debug:
        #print("ints =", ints)
        #print("iminus1 =", iminus1)
        #print("delta =", delta)
    idelta = np.where(delta == 1)[0]
    #if debug:
        #print("idelta =", idelta)
        #print('delta[idelta] =', delta[idelta])
    iminus1_delta = iminus1[idelta]
    #if debug:
        #print("iminus1_delta =", iminus1_delta)
    return iminus1_delta
