"""
Defines various tables that don't fit in other sections:
  - MinorTables
    - _read_ibulk(self)
    - _read_fol(self)
    - _read_gpl(self)
    - _read_extdb(self)
    - _skip_pcompts(self)
    - _read_pcompts(self)
    - _read_meff(self)
    - _read_intmod(self)
    - _read_frl(self)
    - _read_sdf(self)
    - _read_aemonpt(self)
    - _read_monitor(self)
    - _read_r1tabrg(self, data, ndata)
    - _read_hisadd(self)
    - _skip_matrix_mat(self)
    - _get_matrix_row_fmt_nterms_nfloats(self, nvalues, tout)
    - _read_matrix(self, table_name)
    - _read_matpool_matrix(self)
    - _read_matrix_mat(self)

  - grids_comp_array_to_index(grids1, comps1, grids2, comps2,
                              make_matrix_symmetric)
"""

from __future__ import print_function
from struct import unpack
from six import b
import numpy as np
import scipy
from pyNastran.op2.tables.matrix import Matrix
from pyNastran.op2.tables.design_response import WeightResponse, FlutterResponse, Convergence
from pyNastran.op2.op2_interface.op2_common import OP2Common


class MinorTables(OP2Common):
    """reads various tables that don't fit into a larger category"""
    def __init__(self):
        OP2Common.__init__(self)

    def _read_ibulk(self):
        self.table_name = self._read_table_name(rewind=False)
        self.log.debug('table_name = %r' % self.table_name)
        if self.is_debug_file:
            self.binary_debug.write('_read_geom_table - %s\n' % self.table_name)
        self.read_markers([-1])
        if self.is_debug_file:
            self.binary_debug.write('---markers = [-1]---\n')
        data = self._read_record()

        markers = self.get_nmarkers(1, rewind=True)
        marker = -2
        while 1:
            self.read_markers([marker, 1, 0])
            nfields = self.get_marker1(rewind=True)
            if nfields > 0:
                data = self._read_record()
                #self.show_data(data, types='s', endian=None)
            elif nfields == 0:
                #self.show_ndata(100, types='ifs')
                break
            else:
                raise RuntimeError('nfields=%s' % nfields)
            marker -= 1
        marker_end = self.get_marker1(rewind=False)

    def _read_fol(self):
        """
        Reads the FOL table
        Frequency response frequency output list

        +------+---------+-------+-----------------+
        | Word |  Name   | Type  |   Description   |
        +======+=========+=======+=================+
        |  1   | NAME(2) | CHAR4 | Data block name |
        +------+---------+-------+-----------------+
        |  3   |  FREQ   |  RS   |   Frequency     |
        +------+---------+-------+-----------------+
        | Word 3 repeats until End of Record       |
        +------------------------------------------+

        +------+----------+------+-----------------------------+
        | Word |  Name    | Type |   Description               |
        +======+==========+======+=============================+
        |  1   |  WORD1   |  I   | Number of frequencies       |
        +------+----------+------+-----------------------------+
        |  2   |  WORD2   |  I   | Frequency set record number |
        +------+----------+------+-----------------------------+
        |  3   |  WORD3   |  I   | Number of loads             |
        +------+----------+------+-----------------------------+
        |  4   | UNDEF(3) | None | Not used                    |
        +------+----------+------+-----------------------------+
        """
        self.log.debug("table_name = %r" % self.table_name)
        self.table_name = self._read_table_name(rewind=False)
        self.read_markers([-1])
        data = self._read_record()
        self.read_markers([-2, 1, 0])
        data = self._read_record()
        ndata = len(data)
        subtable_name_raw, = unpack(b(self._endian + '8s'), data[:8])
        subtable_name = subtable_name_raw.strip()
        assert subtable_name == b'FOL', 'subtable_name=%r' % subtable_name

        nfloats = (ndata - 8) // 4
        assert nfloats * 4 == (ndata - 8)
        fmt = b(self._endian + '%sf' % nfloats)
        freqs = np.array(list(unpack(fmt, data[8:])), dtype='float32')
        self._frequencies = freqs
        if self.is_debug_file:
            self.binary_debug.write('  recordi = [%r, freqs]\n'  % (subtable_name_raw))
            self.binary_debug.write('  subtable_name=%r\n' % subtable_name)
            self.binary_debug.write('  freqs = %s' % freqs)
        self._read_subtables()

    def _read_gpl(self):
        """reads the GPL table (grid point list?)"""
        self.table_name = self._read_table_name(rewind=False)
        self.log.debug('table_name = %r' % self.table_name)
        if self.is_debug_file:
            self.binary_debug.write('_read_geom_table - %s\n' % self.table_name)
        self.read_markers([-1])
        if self.is_debug_file:
            self.binary_debug.write('---markers = [-1]---\n')
        data = self._read_record()

        markers = self.get_nmarkers(1, rewind=True)
        if self.is_debug_file:
            self.binary_debug.write('---marker0 = %s---\n' % markers)
        n = -2
        while markers[0] != 0:
            self.read_markers([n, 1, 0])
            if self.is_debug_file:
                self.binary_debug.write('---markers = [%i, 1, 0]---\n' % n)

            markers = self.get_nmarkers(1, rewind=True)
            if markers[0] == 0:
                markers = self.get_nmarkers(1, rewind=False)
                break
            data = self._read_record()
            #self.show_data(data, 'i')
            n -= 1
            markers = self.get_nmarkers(1, rewind=True)

    def _read_extdb(self):
        r"""
        fails if a streaming block:
         - nx_spike\extse04c_0.op2
        """
        self.table_name = self._read_table_name(rewind=False)
        self.log.debug('table_name = %r' % self.table_name)
        if self.is_debug_file:
            self.binary_debug.write('_read_geom_table - %s\n' % self.table_name)
        self.read_markers([-1])
        if self.is_debug_file:
            self.binary_debug.write('---markers = [-1]---\n')
        data = self._read_record()

        markers = self.get_nmarkers(1, rewind=True)
        if self.is_debug_file:
            self.binary_debug.write('---marker0 = %s---\n' % markers)

        marker = -2
        while 1:
            try:
                self.read_markers([marker, 1, 0])
            except FortranMarkerError:
                self.show_ndata(100)
                raise
            nfields = self.get_marker1(rewind=True)
            if nfields > 0:
                data = self._read_record()
                #self.show_data(data, types='s', endian=None)
            elif nfields == 0:
                #self.show_ndata(100, types='ifs')
                break
            else:
                raise RuntimeError('nfields=%s' % nfields)
            marker -= 1
        marker_end = self.get_marker1(rewind=False)

    def _skip_pcompts(self):
        """
        Reads the PCOMPTS table (poorly).
        The PCOMPTS table stores information about the PCOMP cards???
        """
        self.log.debug("table_name = %r" % self.table_name)
        table_name = self._read_table_name(rewind=False)

        self.read_markers([-1])
        data = self._skip_record()

        self.read_markers([-2, 1, 0])
        data = self._skip_record()
        #table_name, = self.struct_8s.unpack(data)

        self.read_markers([-3, 1, 0])
        markers = self.get_nmarkers(1, rewind=True)
        if markers != [-4]:
            data = self._skip_record()

        self.read_markers([-4, 1, 0])
        markers = self.get_nmarkers(1, rewind=True)
        if markers != [0]:
            data = self._skip_record()
        else:
            self.read_markers([0])
            return

        self.read_markers([-5, 1, 0])
        data = self._skip_record()

        self.read_markers([-6, 1, 0])
        self.read_markers([0])

    def _read_pcompts(self):
        """
        Reads the PCOMPTS table (poorly).
        The PCOMPTS table stores information about the PCOMP cards???
        """
        self._skip_pcompts()
        return
        #if self.read_mode == 1:
            #return
        #self.log.debug("table_name = %r" % self.table_name)
        #table_name = self._read_table_name(rewind=False)

        #self.read_markers([-1])
        #data = self._read_record()

        #self.read_markers([-2, 1, 0])
        #data = self._read_record()
        #table_name, = self.struct_8s.unpack(data)
        ##print "table_name = %r" % table_name

        #self.read_markers([-3, 1, 0])
        #markers = self.get_nmarkers(1, rewind=True)
        #if markers != [-4]:
            #data = self._read_record()

        #self.read_markers([-4, 1, 0])
        #markers = self.get_nmarkers(1, rewind=True)
        #if markers != [0]:
            #data = self._read_record()
        #else:
            #self.read_markers([0])
            #return

        #self.read_markers([-5, 1, 0])
        #data = self._read_record()

        #self.read_markers([-6, 1, 0])
        #self.read_markers([0])

    def _read_meff(self):
        self.table_name = self._read_table_name(rewind=False)
        self.log.debug('table_name = %r' % self.table_name)
        if self.is_debug_file:
            self.binary_debug.write('_read_geom_table - %s\n' % self.table_name)
        self.read_markers([-1])
        if self.is_debug_file:
            self.binary_debug.write('---markers = [-1]---\n')
        data = self._read_record()

        markers = self.get_nmarkers(1, rewind=True)
        if self.is_debug_file:
            self.binary_debug.write('---marker0 = %s---\n' % markers)
        self.read_markers([-2, 1, 0])
        data = self._read_record()

        for n in [-3, -4, -5, -6, -7, -8]:
            self.read_markers([n, 1, 1])
            markers = self.get_nmarkers(1, rewind=False)
            #print('markers =', markers)
            nbytes = markers[0]*4 + 12
            data = self.f.read(nbytes)
            self.n += nbytes
        n = -9
        self.read_markers([n, 1, 0, 0])

    def _read_intmod(self):
        """reads the INTMOD table"""
        self.table_name = self._read_table_name(rewind=False)
        #self.log.debug('table_name = %r' % self.table_name)
        if self.is_debug_file:
            self.binary_debug.write('_read_geom_table - %s\n' % self.table_name)
        self.read_markers([-1])
        if self.is_debug_file:
            self.binary_debug.write('---markers = [-1]---\n')
        data = self._read_record()
        #print('intmod data1')
        #self.show_data(data)

        markers = self.get_nmarkers(1, rewind=True)
        if self.is_debug_file:
            self.binary_debug.write('---marker0 = %s---\n' % markers)
        self.read_markers([-2, 1, 0])
        data = self._read_record()
        #print('intmod data2')
        #self.show_data(data)

        for n in [-3, -4, -5, -6, -7, -8,]:
            self.read_markers([n, 1, 1])
            markers = self.get_nmarkers(1, rewind=False)
            #print('markers =', markers)
            nbytes = markers[0]*4 + 12
            data = self.f.read(nbytes)
            #print('intmod data%i' % n)
            #self.show_data(data)
            self.n += nbytes

        n = -9
        self.read_markers([n, 1, 0, 0])
        #self.show(50)
        #raise NotImplementedError(self.table_name)

    def _read_frl(self):
        #self.log.debug("table_name = %r" % self.table_name)
        #self.table_name = self._read_table_name(rewind=False)
        #self.read_markers([-1])
        #data = self._read_record()

        #self.read_markers([-2, 1, 0])
        #data = self._read_record()
        self._skip_table(self.table_name)

    def _read_sdf(self):
        self.log.debug("table_name = %r" % self.table_name)
        self.table_name = self._read_table_name(rewind=False)
        self.read_markers([-1])
        data = self._read_record()

        self.read_markers([-2, 1, 0])
        data, ndata = self._read_record_ndata()
        if ndata == 16:
            subtable_name, dummy_a, dummy_b = unpack(b(self._endian + '8sii'), data)
            if self.is_debug_file:
                self.binary_debug.write('  recordi = [%r, %i, %i]\n'  % (
                    subtable_name, dummy_a, dummy_b))
                self.binary_debug.write('  subtable_name=%r\n' % subtable_name)
                assert dummy_a == 170, dummy_a
                assert dummy_b == 170, dummy_b
        else:
            strings, ints, floats = self.show_data(data)
            msg = 'len(data) = %i\n' % ndata
            msg += 'strings  = %r\n' % strings
            msg += 'ints     = %r\n' % str(ints)
            msg += 'floats   = %r' % str(floats)
            raise NotImplementedError(msg)

        self.read_markers([-3, 1, 1])

        markers0 = self.get_nmarkers(1, rewind=False)
        record = self.read_block()

        #data = self._read_record()
        self.read_markers([-4, 1, 0, 0])
        #self._read_subtables()

    def _read_aemonpt(self):
        """reads the AEMONPT table"""
        #self.log.debug("table_name = %r" % self.table_name)
        table_name = self._read_table_name(rewind=False)

        #print('-----------------------')
        #print('record 1')
        self.read_markers([-1])
        data = self._read_record()
        #self.show_data(data)
        if self.read_mode == 2:
            a, bi, c, d, e, f, g = unpack(b('%s7i'% self._endian), data)
            assert a == 101, a
            assert bi == 1, bi
            assert c == 27, c
            assert d == 1, d
            assert e == 6, e
            assert f == 0, f
            assert g == 0, g
        #print('-----------------------')
        #print('record 2')
        self.read_markers([-2, 1, 0])
        data = self._read_record()
        #if self.read_mode == 2:
        word, = unpack(b('%s8s' % self._endian), data)
        assert word == b'AECFMON ', word
        #self.show_data(data)
        #print('-----------------------')
        #print('record 3')
        self.read_markers([-3, 1, 0])
        data = self._read_record()
        #self.show_data(data)

        if self.read_mode == 2:
            ndata = len(data)
            assert ndata == 108, ndata
            n = 8 + 56 + 20 + 12 + 12
            aero, name, comps, cp, bi, c, d, coeff, word, e, f, g = unpack(b('8s 56s 5i 4s 8s 3i'), data[:n])
            print('aero=%r' % aero)
            print('name=%r' % name)
            print('comps=%r cp=%s b,c,d=(%s, %s, %s)' % (comps, cp, bi, c, d))
            print('coeff=%r' % coeff)
            print('word=%r (e, f, g)=(%s, %s, %s)' % (word, e, f, g)) # (1, 2, 0)
            assert cp == 2, cp
            assert bi == 0, bi
            assert c == 0, c
            assert d == 0, d
            assert e == 1, e
            assert f == 2, f
            assert g == 0, g

        #print('-----------------------')
        #print('record 4')
        self.read_markers([-4, 1, 0])
        #data = self._read_record()
        #self.show_data(data)
        self.read_markers([0])

        #print('-----------------------')
        #print('end')
        #self.show(200)
        #aaa

    def _read_monitor(self):
        """reads the MONITOR table"""
        self.log.debug("table_name = %r" % self.table_name)
        table_name = self._read_table_name(rewind=False)

        #print('-----------------------')
        #print('record 1')
        self.read_markers([-1])
        data = self._read_record()
        #self.show_data(data)
        if self.read_mode == 2:
            a, bi, c, d, e, f, g = unpack(b('%s7i' % self._endian), data)
            assert a == 101, a
            assert bi == 1, bi
            assert c == 27, c
            assert d == 0, d
            assert e == 6, e
            assert f == 0, f
            assert g == 0, g
        #print('-----------------------')
        #print('record 2')
        self.read_markers([-2, 1, 0])
        data = self._read_record()
        if self.read_mode == 2:
            word, = unpack(b('%s8s' % self._endian), data)
            assert word == b'STCFMON ', word
        #self.show_data(data)
        #print('-----------------------')
        #print('record 3')
        self.read_markers([-3, 1, 0])
        data = self._read_record()
        #self.show_data(data[96:108])

        if self.read_mode == 2:
            ndata = len(data)
            assert ndata == 108, ndata
            aero, name, comps, cp, x, y, z, coeff, word, column, cd, ind_dof = unpack(b'8s 56s 2i 3f 4s 8s 3i', data[:108])
            print('aero=%r' % aero)
            print('name=%r' % name)
            print('comps=%s cp=%s (x, y, z)=(%s, %s, %s)' % (comps, cp, x, y, z))
            print('coeff=%r' % coeff)
            print('word=%r (column, cd, ind_dof)=(%s, %s, %s)' % (word, column, cd, ind_dof))
            assert cp == 2, cp
            assert x == 0.0, x
            assert y == 0.0, y
            assert d == 0.0, z
            assert column == 1, column
            assert cd == 2, cd
            assert ind_dof == 0, ind_dof
            self.monitor_data = [{
                'name' : name,
                'cp' : cp,
                'cd' : cd,
                'xyz' : [x, y, z],
                'comps' : comps,
            }]

        #print('-----------------------')
        #print('record 4')
        self.read_markers([-4, 1, 0])
        #data = self._read_record()
        #self.show_data(data)
        self.read_markers([0])

        #print('-----------------------')
        #print('end')
        #self.show(200)
        #aaa

    def _read_r1tabrg(self, data, ndata):
        """
        Design Responses:
          - Weight
          - Flutter Speed
          - Stress
          - Strain
          - Displacement
        """
        if self._table4_count == 0:
            self._count += 1
        self._table4_count += 1

        #if self._table4_count == 0:
            #self._count += 1
        #self._table4_count += 1

        if self.read_mode == 1:
            assert data is not None, data
            assert len(data) > 12, len(data)
            Type, = unpack(self._endian + 'i', data[8:12])
            #assert Type in [1, 6, 10, 84], Type
            if Type == 1:
                if self.weight_response is None:
                    self.weight_response = WeightResponse()
                else:
                    self.weight_response.n += 1
            elif Type == 4:
                #TYPE =4 EIGN or FREQ
                #8 MODE I Mode number
                #9 APRX I Approximation code
                pass
            elif Type == 5:
                #TYPE =5 DISP
                #8 COMP I Displacement component
                #9 UNDEF None
                #10 GRID I Grid identification number
                pass
            elif Type == 15:
                # CEIG
                #8 MODE I Mode number
                #9 ICODE I 1: Real component or 2: Imaginary component
                pass
            elif Type == 84:
                if self.flutter_response is None:
                    self.flutter_response = FlutterResponse()
                else:
                    self.flutter_response.n += 1
            return ndata
            #else: # response not added...
                #pass

        read_r1tabrg = True
        if read_r1tabrg:
            #self.show_data(data, types='ifs', endian=None)
            out = unpack(self._endian + 'iii 8s iiii i iiiii', data)
            # per the R1TAB DMAP page:
            #   all indicies are downshift by 1
            #   indices above out[3] are off by +2 because of the 2 field response_label
            internal_id = out[0]
            dresp_id = out[1]
            response_type = out[2]
            response_label = out[3].strip()
            # -1 for 2 field wide response_label
            region = out[4]
            subcase = out[5]
            type_flag = out[12]  # no meaning per MSC DMAP 2005
            seid = out[13]

            if response_type == 1:
                #                             -----  WEIGHT RESPONSE  -----
                # ---------------------------------------------------------------------------------
                #  INTERNAL  DRESP1  RESPONSE  ROW  COLUMN  LOWER     INPUT      OUTPUT     UPPER
                #     ID       ID     LABEL     ID    ID    BOUND     VALUE       VALUE     BOUND
                # ---------------------------------------------------------------------------------
                #       1       1    WEIGHT     3     3       N/A   2.9861E+05  2.9852E+05   N/A
                #(1, 1,    1, 'WEIGHT  ', 0, 1011, 3, 3, 0, 0, 0, 0, 0, 0)
                #(1, 1000, 1, 'W       ', 0, 1,    3, 3, 0, 0, 0, 0, 0, 0)
                #print(out)
                #row_id = out[4]

                # these should be blank?
                row_id = out[6]
                column_id = out[7]
                seid_weight = out[8]

                assert np.abs(out[8:-1]).sum() == 0.0, 'out=%s 8=%s' % (out, out[8:-1])
                assert out[-1] in [0, 1, 2, 3, 4, 5], out
                #dunno_8 = out[8]
                #dunno_9 = out[9]
                #dunno_10 = out[10]
                #dunno_11 = out[11]
                #dunno_12 = out[12]
                #dunno_13 = out[13]
                #msg = ('WEIGHT - response_type=%r response_label=%r row_id=%r column_id=%r '
                       #'6=%r 7=%r 8=%r 9=%r 10=%r 11=%r 12=%r 13=%r' % (
                           #response_type, response_label, row_id, column_id,
                           #dunno_6, dunno_7, dunno_8, dunno_9, dunno_10, dunno_11, dunno_12, dunno_13))
                #out = unpack(self._endian + 'iii 8s iiff f fffff', data)
                #print(out)
                msg = 'WEIGHT - label=%r region=%s subcase=%s row_id=%r column_id=%r' % (
                    response_label, region, subcase, row_id, column_id)
                self.weight_response.append(internal_id, dresp_id, response_label, region,
                                            subcase, type_flag, seid,
                                            row_id, column_id)
                #print(msg)
                #self.log.debug(msg)
            elif response_type == 5:  # DISP
                # out = (1, 101, 5, 'DISP1   ', 101, 1, 3, 0, 1, 0, 0, 0, 0, 0)

                #print(out[6:])
                # (3,   0,  1,    0,   0,   0,   0,   0)
                # (???, NA, comp, ???, ???, ???, ???, ???)
                pass
            elif response_type == 6:  # STRESS
                #                              -----   STRESS RESPONSES   -----
                #  -------------------------------------------------------------------------------------------
                #   INTERNAL  DRESP1  RESPONSE  ELEMENT   VIEW   COMPONENT  LOWER   INPUT    OUTPUT    UPPER
                #      ID       ID     LABEL       ID    ELM ID     NO.     BOUND   VALUE     VALUE    BOUND
                #  -------------------------------------------------------------------------------------------
                #         21      209  S09L      144747             17       N/A   4.85E+04  5.00E+04  5.00E+04
                # (21, 209, 6, 'S09L    ', 30, 1011, 17, 0, 144747, 0, 0, 0, 0, 0)
                stress_code = out[6]
                pid = out[8]
                msg = ('STRESS - response_type=%r label=%r region=%s subcase=%s '
                       'stress_code=%s pid=%s' % (
                           response_type, response_label, region, subcase,
                           stress_code, pid))

            #elif response_type == 5:  # DISP
                #pass
            #elif response_type == 7:  # STRAIN
                #pass
            elif response_type == 10:  # CSTRESS
                stress_code = out[6]
                ply = out[7]
                pid = out[8]  # is this element id?
                msg = 'CSTRESS - label=%r region=%s subcase=%s stress_code=%s ply=%s pid=%s' % (
                    response_label, region, subcase, stress_code, ply, pid)
                #print(msg)
            #elif response_type == 10:  # CSTRAIN
                #pass
            elif response_type == 24:  # FRSTRE
                #8 ICODE I Stress item code
                #9 UNDEF None
                #10 ELID I Element identification number
                #11 FREQ RS Frequency
                #12 IFLAG I Integrated response flag. See Remark 20 of
                #DRESP1.
                #Value is -1 to -6, for SUM, AVG, SSQ,
                pass
            elif response_type == 28:  # RMSACCL
                #8 COMP I RMS Acceleration component
                #9 RANDPS I RANDPS entry identification number
                #10 GRID I Grid identification number
                #11 DMFREQ RS Dummy frequency for internal use
                pass
            elif response_type == 84:  # FLUTTER  (iii, label, mode, (Ma, V, rho), flutter_id, fff)
                out = unpack(self._endian + 'iii 8s iii fff i fff', data)
                mode = out[6]
                mach = out[7]
                velocity = out[8]
                density = out[9]
                flutter_id = out[10]
                msg = ('FLUTTER - _count=%s label=%r region=%s subcase=%s mode=%s '
                       'mach=%s velocity=%s density=%s flutter_id=%s' % (
                           self._count, response_label, region, subcase, mode,
                           mach, velocity, density, flutter_id))
                self.flutter_response.append(internal_id, dresp_id, response_label, region,
                                             subcase, type_flag, seid,
                                             mode, mach, velocity, density, flutter_id)
                #print(msg)
                #self.log.debug(msg)
            else:
                self.log.debug('R1TABRG response response_type=%s not supported' % response_type)
                #raise NotImplementedError(response_type)
            assert len(out) == 14, len(out)
        #self.response1_table[self._count] = out
        return ndata

    def _read_hisadd(self):
        """optimization history (SOL200) table"""
        self.table_name = self._read_table_name(rewind=False)

        if self.read_mode == 1:
            self.read_markers([-1])
            self._skip_record()
            self.read_markers([-2, 1, 0])
            self._skip_record()
            self.read_markers([-3, 1, 0])

            if self.convergence_data is None:
                data = self._read_record()
                ndvs = len(data) // 4 - 7
                self.convergence_data = Convergence(ndvs)
            else:
                self._skip_record()
                self.convergence_data.n += 1

            self.read_markers([-4, 1, 0, 0])
            return

        if self.is_debug_file:
            self.binary_debug.write('_read_geom_table - %s\n' % self.table_name)
        #self.log.info('----marker1----')
        self.read_markers([-1])
        if self.is_debug_file:
            self.binary_debug.write('---markers = [-1]---\n')
        data = self._read_record()  # ()102, 303, 0, 0, 0, 0, 0) date???
        #print('hisadd data1')
        #self.show_data(data)

        #self.log.info('----marker2----')
        markers = self.get_nmarkers(1, rewind=True)
        if self.is_debug_file:
            self.binary_debug.write('---marker0 = %s---\n' % markers)
        self.read_markers([-2, 1, 0])
        data = self._read_record()  # ('HISADD', )
        #print('hisadd data2')
        #self.show_data(data)

        #self.log.info('----marker3----')
        self.read_markers([-3, 1, 0])
        data = self._read_record()

        (design_iter, iconvergence, conv_result, obj_intial, obj_final,
         constraint_max, row_constraint_max) = unpack(b(self._endian + '3i3fi'), data[:28])
        if iconvergence == 1:
            iconvergence = 'soft'
        elif iconvergence == 2:
            iconvergence = 'hard'
        elif iconvergence == 6:
            self.log.warning('HISADD iconverge=6')
            iconvergence = '???'
        else:
            msg = 'iconvergence=%s\n' % iconvergence
            self.show_data(data, types='ifs', endian=None)
            raise NotImplementedError(msg)

        if conv_result == 0:
            conv_result = 'no'
        elif conv_result == 1:
            conv_result = 'soft'
        elif conv_result == 2:
            conv_result = 'hard'
        elif conv_result in [3, 4]:
            self.log.warning('HISADD conv_result=%s' % conv_result)
            # not sure why this happens, but the field is wrong
            # it seems to apply to one step before this one
            conv_result = 'best_design'
        else:
            self.log.debug('design_iter=%s iconvergence=%s conv_result=%s obj_intial=%s '
                           'obj_final=%s constraint_max=%s row_constraint_max=%s' % (
                               design_iter, iconvergence, conv_result, obj_intial,
                               obj_final, constraint_max, row_constraint_max))
            raise NotImplementedError('conv_result=%s' % conv_result)
        #self.log.debug('design_iter=%s iconvergence=%s conv_result=%s obj_intial=%s '
                       #'obj_final=%s constraint_max=%s row_constraint_max=%s' % (
                           #design_iter, iconvergence, conv_result, obj_intial,
                           #obj_final, constraint_max, row_constraint_max))

        ndvs = len(data) // 4 - 7
        desvar_values = unpack('%sf' % ndvs, data[28:])

        self.convergence_data.append(design_iter, iconvergence, conv_result, obj_intial,
                                     obj_final, constraint_max, row_constraint_max, desvar_values)
        self.read_markers([-4, 1, 0, 0])

    def _skip_matrix_mat(self):
        table_name = self._read_table_name(rewind=False, stop_on_failure=True)
        self.read_markers([-1])
        data = self._skip_record()

        self.read_markers([-2, 1, 0])
        data = self._skip_record()

        itable = -3
        niter = 0
        niter_max = 100000000

        jj = 1
        while niter < niter_max:
            #nvalues = self.get_marker1(rewind=True)
            #print('nvalues4a =', nvalues)
            self.read_markers([itable, 1])
            one = self.get_marker1(rewind=False)

            if one:  # if keep going
                nvalues = self.get_marker1(rewind=True)
                while nvalues >= 0:
                    nvalues = self.get_marker1(rewind=False)
                    data = self._skip_block()
                    nvalues = self.get_marker1(rewind=True)
                jj += 1
            else:
                nvalues = self.get_marker1(rewind=False)
                assert nvalues == 0, nvalues
                return
            itable -= 1
            niter += 1
        raise RuntimeError('this should never happen; n=%s' % niter_max)

    def _get_matrix_row_fmt_nterms_nfloats(self, nvalues, tout):
        """
        +------+---------------------------+
        | Type | Meaning                   |
        +------+---------------------------+
        |  1   | Real, single precision    |
        |  2   | Real, double precision    |
        |  3   | Complex, single precision |
        |  4   | Complex, double precision |
        +------+---------------------------+
        """
        if tout == 1:
            nfloats = nvalues
            nterms = nvalues
            fmt = self._endian + 'i %if' % nfloats
        elif tout == 2:
            nfloats = nvalues // 2
            nterms = nvalues // 2
            fmt = self._endian + 'i %id' % nfloats
        elif tout == 3:
            nfloats = nvalues
            nterms = nvalues // 2
            fmt = self._endian + 'i %if' % nfloats
        elif tout == 4:
            nfloats = nvalues // 2
            nterms = nvalues // 4
            fmt = self._endian + 'i %id' % nfloats
        else:
            raise RuntimeError('tout = %s' % tout)
        return fmt, nfloats, nterms


    def _read_matrix(self, table_name):
        """
        general method for reading matrices and MATPOOL matrices

        .. todo:: Doesn't support checking matrices vs. MATPOOLs
        .. todo:: MATPOOLs are disabled because they're not parsed properly
        """
        i = self.f.tell()
        # if we skip on read_mode=1, we don't get debugging
        # if we just use read_mode=2, some tests fail
        #
        #if self.read_mode == 1:
        if self.read_mode == 2 and not self.debug_file:
            try:
                self._skip_matrix_mat()  # doesn't work for matpools
            except:
                self._goto(i)
                self._skip_table(table_name)
            return

        #enable_matpool = True
        #if enable_matpool:
        try:
            self._read_matrix_mat()
        except:
            # read matpool matrix
            self._goto(i)
            try:
                self._read_matpool_matrix()
            except:
                raise
                #self._goto(i)
                #self._skip_table(self.table_name)
        #else:
            #try:
                #self._read_matrix_mat()
            #except:
                #self._goto(i)
                #self._skip_table(self.table_name)

    def _read_matpool_matrix(self):
        """
        Reads a MATPOOL matrix

        MATPOOL matrices are always sparse

        +------+-----------------+
        | Form | Meaning         |
        +======+=================+
        |  1   | Square          |
        |  2   | Rectangular     |
        |  6   | Symmetric       |
        |  9   | Pseudo identity |
        +------+-----------------+
        """
        table_name = self._read_table_name(rewind=False, stop_on_failure=True)
        utable_name = table_name.decode('utf-8')
        self.read_markers([-1])
        data = self._read_record()

        self.read_markers([-2, 1, 0])
        data = self._read_record()

        self.read_markers([-3, 1, 0])
        data = self._read_record()

        nvalues = len(data) // 4
        assert len(data) % 4 == 0, len(data) / 4.

        header = unpack('3i 8s 7i', data[:48]) # 48=4*12
        assert header[:3] == (114, 1, 120), 'header[:3]=%s header=%s' % (header[:3], header)

        # ncols_gset is needed for form=9
        #  list of header values:
        #    4:5   matrix name
        #    6     placeholder
        #    7     matrix shape (1 = square, 2 or 9 = rectangular, 6 = symmetric)
        #    8     input type flag (1 = single, 2 = double, 3 = complex single,
        #                           4 = complex double)
        #    9     output type flag (0 = precision set by system cell,
        #                            1 = single, 2 = double, 3 = complex single,
        #                            4 = complex double)
        #   10     complex flag (0 = real/imaginary, >0 = magnitude/phase)
        #   11     placeholder
        #   12     number of columns in the G set (only necessary for matrix
        #                                          shape 9)
        matrix_name, junk1, matrix_shape, tin, tout, is_phase, junk2, ncols_gset = header[3:]
        matrix_name = matrix_name.strip()

        #self.log.debug('matrix_name=%s, junk1=%s, matrix_shape=%s, tin=%s, tout=%s, '
                       #'is_phase=%s, junk2=%s, ncols_gset=%s' % (
                           #matrix_name, junk1, matrix_shape, tin, tout,
                           #is_phase, junk2, ncols_gset))

        is_complex = False
        if tin > 2 or tout > 2:
            is_complex = True
            assert is_phase == 0, 'is_phase=%s' % is_phase
            imags = []

        if tout == 1:
            dtype = 'float32'
            fdtype = self.fdtype
        elif tout == 2:
            dtype = 'float64'
            fdtype = self.double_dtype
        elif tout == 3:
            dtype = 'complex64'
            fdtype = self.fdtype
        elif tout == 4:
            dtype = 'complex128'
            fdtype = self.double_dtype
        else:
            dtype = '???'
            msg = ('matrix_name=%s, junk1=%s, matrix_shape=%s, tin=%s, tout=%s, '
                   'is_phase=%s, junk2=%s, ncols_gset=%s' % (
                       matrix_name, junk1, matrix_shape, tin, tout,
                       is_phase, junk2, ncols_gset))
            self.log.warning(msg)
            raise RuntimeError(msg)

        is_symmetric = matrix_shape == 6
        is_phase_flag = is_phase > 0

        if tout in [1, 3]:
            # works for float32, complex64
            ints = np.fromstring(data[48:], dtype=self.idtype)
            floats = np.fromstring(data[48:], dtype=self.fdtype)
            temp_ints = ints
        else:
            # works for float64, complex128
            temp_ints = np.fromstring(data[48:], dtype=self.idtype)

        # find the first index with ()-1,-1)
        iminus1 = np.where(temp_ints[:-1] == -1)[0]
        double_minus1 = (iminus1[:-1] + 1 == iminus1[1:])[:-1]

        # the field after our stop
        # we'll handle the off by 1 later with arange
        istop = iminus1[:-2][double_minus1]

        # 2 fields after is the start position
        # add on a 0 to the beginning to account for the starting position
        # istart defines icol
        istart = np.hstack([0, istop[:-1] + 2])

        col_nids_short = temp_ints[istart]
        col_dofs_short = temp_ints[istart+1]
        #nj2 = len(istart)  ## TODO: why is this wrong???

        row_nids = []
        row_dofs = []
        col_nids = []
        col_dofs = []
        reals = []
        imags = []
        for col_nidi, col_dofi, istarti, istopi in zip(
            col_nids_short, col_dofs_short, istart + 2, istop):

            ## TODO: preallocate arrays
            imag = None
            # The float32/complex64 blocks are really simple
            # because we can just use the data is from the temp_ints block.
            # We calculate istart/istop and directly access the float data.
            #
            # The float64/complex128 blocks are more complicated.
            # It's easier to just use istart/istop to calculate datai
            # case that, and then slice it.
            #
            # In all cases, we use istart/istop to calculate row/col nid/dof
            # because they are always int32.  Only the real/imag types change.
            if dtype == 'float32':
                irow = np.arange(istarti, istopi-1, step=3, dtype='int32')
                real = floats[irow + 2]
            elif dtype == 'complex64':
                irow = np.arange(istarti, istopi-1, step=4, dtype='int32')
                real = floats[irow + 2]
                imag = floats[irow + 3]

            elif dtype == 'float64':
                datai = data[48+(istarti*4) : 48+(istopi*4)]
                irow = np.arange(istarti, istopi-1, step=4, dtype='int32')
                assert len(datai) % 8 == 0, len(datai) / 8
                real = np.fromstring(datai, dtype=fdtype)[1::2]

            elif dtype == 'complex128':
                datai = data[48+(istarti*4) : 48+(istopi*4)]

                # iword
                # -----
                #   0    1    3     5   <---- iword
                #   1    1    2     2   <---- nwords
                # (nid, dof, real, imag)
                irow = np.arange(istarti, istopi-1, step=6, dtype='int32')
                assert len(datai) % 8 == 0, len(datai) / 8
                floats = np.fromstring(datai, dtype=fdtype)

                # ndoubles
                # --------
                #  <---0--->   1     2    <----- iword
                #      1       1     1    <----- nwords
                # (nid, dof, real, imag)
                real = floats[1::3]
                imag = floats[2::3]
            else:
                msg = '%s is not supported' % dtype
                self.log.error(msg)
                raise RuntimeError(msg)

            if len(irow) != len(real):
                msg = 'nrow=%s nreal=%s nimag=%s' % (len(irow), len(real), len(imag))
                raise RuntimeError(msg)

            # the row index; [1, 2, ..., 43]
            row_nid = temp_ints[irow]

            # the dof; [0, 0, ..., 0.]
            row_dof = temp_ints[irow + 1]
            urow_dof = np.unique(row_dof)
            for udofi in urow_dof:
                if udofi not in [0, 1, 2, 3, 4, 5, 6]:
                    msg = 'udofi=%s is invalid; must be in [0, 1, 2, 3, 4, 5, 6]; dofs=%s' % (
                        udofi, np.asarray(urow_dof, dtype='int32').tolist())
                    raise ValueError(msg)

            ni = len(irow)
            col_nid = np.ones(ni, dtype='int32') * col_nidi
            col_dof = np.ones(ni, dtype='int32') * col_dofi

            row_nids.append(row_nid)
            row_dofs.append(row_dof)
            col_nids.append(col_nid)
            col_dofs.append(col_dof)
            reals.append(real)
            imags.append(imag)

        row_nids_array = np.hstack(row_nids)
        row_dofs_array = np.hstack(row_dofs)

        col_nids_array = np.hstack(col_nids)
        col_dofs_array = np.hstack(col_dofs)
        real_array = np.hstack(reals)
        if is_complex:
            complex_array = np.hstack(imags)
            assert len(real_array) == len(complex_array)
            real_imag_array = real_array + 1.j * complex_array
        else:
            real_imag_array = real_array

        # TODO: this is way slower than it should be
        #       because we didn't preallocate the data and the
        #       horrific grids_comp_array_to_index function
        grids1 = col_nids_array
        comps1 = col_dofs_array

        grids2 = row_nids_array
        comps2 = row_dofs_array
        assert len(grids1) == len(comps1), 'ngrids1=%s ncomps1=%s' % (len(grids1), len(comps1))
        assert len(grids1) == len(grids2), 'ngrids1=%s ngrids2=%s' % (len(grids1), len(grids2))
        assert len(comps1) == len(comps2), 'ncomps1=%s ncomps2=%s' % (len(comps1), len(comps2))

        apply_symmetry = True
        make_matrix_symmetric = apply_symmetry and matrix_shape == 'symmetric'
        j1, j2, nj1, nj2, nj = grids_comp_array_to_index(
            grids1, comps1, grids2, comps2, make_matrix_symmetric)
        assert len(j1) == len(j2), 'nj1=%s nj2=%s' % (len(j1), len(j2))
        assert len(grids1) == len(real_array), 'ngrids1=%s nreals=%s' % (len(j1), len(real_array))

        # not 100% on these, they might be flipped
        #ncols = len(np.unique(j1))
        #mrows = len(np.unique(j2))

        if is_symmetric:
            mrows = nj
            ncols = nj
            #print('  j1 =', j1)
            #print('  j2 =', j2)
        else:
            ncols = nj1
            mrows = nj2

        try:
            matrix = scipy.sparse.coo_matrix(
                (real_imag_array, (j2, j1)),
                shape=(mrows, ncols), dtype=dtype)
        except ValueError:
            msg = 'Passed all the checks; cannot build MATPOOL sparse matrix...\n'
            spaces = '                                          '
            msg += '%sname=%s dtype=%s nrows=%s ncols=%s nj1=%s nj2=%s nj=%s' % (
                spaces, table_name, dtype, mrows, ncols, nj1, nj2, nj)
            self.log.error(msg)
            raise


        # enforce symmetry if necessary
        if make_matrix_symmetric:
            # get the upper and lower triangular matrices
            upper_tri = scipy.sparse.triu(matrix)
            lower_tri = scipy.sparse.tril(matrix)

            # extracts a [1, 2, 3, ..., n] off the diagonal of the matrix
            # diagonal_array = diagional(upper_tri)
            #
            # make it a diagonal matrix
            # diagi = diags(diagonal_array)
            diagi = scipy.sparse.diags(scipy.sparse.diagional(upper_tri))

            # Check to see which triangle is populated.
            # If they both are, make sure they're equal
            # or average them and throw a warning
            lnnz = (lower_tri - diagi).nnz
            unnz = (upper_tri - diagi).nnz
            assert isinstance(lnnz, int), type(lnnz)
            assert isinstance(unnz, int), type(unnz)

            # both upper and lower triangle are populated
            if lnnz > 0 and unnz > 0:
                upper_tri_t = upper_tri.T
                if lower_tri == upper_tri_t:
                    matrix = upper_tri + upper_tri_t - diagi
                else:
                    self.log.warning(
                        'Matrix marked as symmetric does not contain '
                        'symmetric data.  Data will be symmetrized.')
                    matrix = (matrix + matrix.T) / 2.
            elif lnnz > 0:
                #  lower triangle is populated
                matrix = lower_tri + lower_tri.T - diagi
            elif unnz > 0:
                #  upper triangle is populated
                matrix = upper_tri + upper_tri_t - diagi
            else:
                # matrix is diagonal (or null)
                matrix = diagi
            data = matrix

            # matrix is symmetric, but is not stored as symmetric
            matrix_shape = 'rectangular'

        m = Matrix(table_name, is_matpool=True, form=matrix_shape)
        m.data = matrix
        m.col_nid = col_nids_array
        m.col_dof = col_dofs_array
        m.row_nid = row_nids_array
        m.row_dof = row_dofs_array
        m.form = matrix_shape
        self.matrices[utable_name] = m
        self.log.debug(m)

        self.read_markers([-4, 1, 0])
        data = self._read_record()

        if len(data) == 12:
            self.read_markers([-5, 1, 0, 0])
            return
        raise RuntimeError('failed on read_matpool_matrix')

    def _read_matrix_mat(self):
        """
        Matrix Trailer:
        +------+---------------------------------------------------+
        | Word | Contents                                          |
        +======+===================================================+
        |  1   | Number of columns in matrix                       |
        |  2   | Number of rows in matrix                          |
        |  3   | Form of the matrix                                |
        |  4   | Type of matrix                                    |
        |  5   | Largest number of nonzero words among all columns |
        |  6   | Density of the matrix multiplied by 10000         |
        |  7   | Size in blocks                                    |
        |  8   | Maximum string length over all strings            |
        |  9   | Number of strings                                 |
        |  10  | Average bandwidth                                 |
        |  11  | Maximum bandwidth                                 |
        |  12  | Number of null columns                            |
        +------+---------------------------------------------------+

        +------+--------------------------------+
        | Form | Meaning                        |
        +======+================================+
        |  1   | Square                         |
        |  2   | Rectangular                    |
        |  3   | Diagonal                       |
        |  4   | Lower triangular factor        |
        |  5   | Upper triangular factor        |
        |  6   | Symmetric                      |
        |  8   | Identity                       |
        |  9   | Pseudo identity                |
        |  10  | Cholesky factor                |
        |  11  | Trapezoidal factor             |
        |  13  | Sparse lower triangular factor |
        |  15  | Sparse upper triangular factor |
        +------+--------------------------------+

        +------+---------------------------+
        | Type | Meaning                   |
        +======+===========================+
        |  1   | Real, single precision    |
        |  2   | Real, double precision    |
        |  3   | Complex, single precision |
        |  4   | Complex, double precision |
        +------+---------------------------+
        """
        allowed_forms = [1, 2, 3, 4, 5, 6, 8, 9, 10, 11, 13, 15]
        #self.log.debug('----------------------------------------------------------------')
        table_name = self._read_table_name(rewind=False, stop_on_failure=True)
        self.read_markers([-1])
        data = self._read_record()

        # old-bad
        #matrix_num, form, mrows, ncols, tout, nvalues, g = unpack(self._endian + '7i', data)

        #           good   good   good  good  ???    ???
        matrix_num, ncols, mrows, form, tout, nvalues, g = unpack(self._endian + '7i', data)
        #print('g =', g)

        m = Matrix(table_name, form=form)
        self.matrices[table_name.decode('utf-8')] = m

        # matrix_num is a counter (101, 102, 103, ...)
        # 101 will be the first matrix 'A' (matrix_num=101),
        # then we'll read a new matrix 'B' (matrix_num=102),
        # etc.
        #
        # the matrix is Mrows x Ncols
        #
        # it has nvalues in it
        #
        # tout is the precision of the matrix
        # 0 - set precision by cell
        # 1 - real, single precision (float32)
        # 2 - real, double precision (float64)
        # 3 - complex, single precision (complex64)
        # 4 - complex, double precision (complex128)

        # form (bad name)
        # 1 - column matrix
        # 2 - factor matrix
        # 3 - factor matrix
        if tout == 1:
            dtype = 'float32'
        elif tout == 2:
            dtype = 'float64'
        elif tout == 3:
            dtype = 'complex64'
        elif tout == 4:
            dtype = 'complex128'
        else:
            dtype = '???'
            msg = ('unexpected tout for %s: matrix_num=%s form=%s '
                   'mrows=%s ncols=%s tout=%s nvalues=%s g=%s'  % (
                       table_name, matrix_num, form, mrows, ncols, tout, nvalues, g))
            self.log.warning(msg)
            raise RuntimeError(msg)

        if form == 1:
            if ncols != mrows:
                self.log.warning('unexpected size for %s; form=%s mrows=%s ncols=%s' % (
                    table_name, form, mrows, ncols))
        elif form not in allowed_forms:
            self.log.error('name=%r matrix_num=%s form=%s mrows=%s '
                           'ncols=%s tout=%s nvalues=%s g=%s' % (
                               table_name, matrix_num, form, mrows, ncols, tout, nvalues, g))
            raise RuntimeError('form=%s; allowed=%s' % (form, allowed_forms))
        #self.log.info('name=%r matrix_num=%s form=%s mrows=%s ncols=%s tout=%s nvalues=%s g=%s' % (
            #table_name, matrix_num, form, mrows, ncols, tout, nvalues, g))

        self.read_markers([-2, 1, 0])
        data = self._read_record()

        if len(data) == 16:
            name, ai, bi = unpack(self._endian + '8s 2i', data)
            assert ai == 170, ai
            assert bi == 170, bi
        else:
            self.log.warning('unexpected matrix length=%s' % len(data))
            self.log.warning(self.show_data(data, types='if'))

        itable = -3
        j = None

        niter = 0
        niter_max = 100000000

        GCi = []
        GCj = []
        reals = []
        jj = 1
        while niter < niter_max:
            #nvalues = self.get_marker1(rewind=True)
            self.read_markers([itable, 1])
            one = self.get_marker1(rewind=False)

            if one:  # if keep going
                nvalues = self.get_marker1(rewind=True)

                while nvalues >= 0:
                    nvalues = self.get_marker1(rewind=False)
                    fmt, nfloats, nterms = self._get_matrix_row_fmt_nterms_nfloats(nvalues, tout)
                    GCjj = [jj] * nterms
                    GCj += GCjj

                    #-----------
                    data = self.read_block()
                    #self.show_data(data)
                    #print('***itable=%s nvalues=%s fmt=%r' % (itable, nvalues, fmt))
                    out = unpack(fmt, data)
                    ii = out[0]
                    values = out[1:]

                    GCii = list(range(ii, ii + nterms))
                    GCi += GCii
                    reals += values
                    nvalues = self.get_marker1(rewind=True)
                    if self.debug_file:
                        self.binary_debug.write('  GCi = %s\n' % GCii)
                        self.binary_debug.write('  GCj = %s\n' % GCjj)
                        self.binary_debug.write('  reals/imags = %s\n' % str(values))
                assert len(GCi) == len(GCj), 'nGCi=%s nGCj=%s' % (len(GCi), len(GCj))
                if tout in [1, 2]:
                    assert len(GCi) == len(reals), 'tout=%s nGCi=%s nReals=%s' % (tout, len(GCi), len(reals))
                else:
                    assert len(GCi)*2 == len(reals), 'tout=%s nGCi=%s nReals=%s' % (tout, len(GCi)*2, len(reals))
                jj += 1
            else:
                nvalues = self.get_marker1(rewind=False)
                assert nvalues == 0, nvalues
                # print('nvalues =', nvalues)
                # print('returning...')

                #assert max(GCi) <= mrows, 'GCi=%s GCj=%s mrows=%s' % (GCi, GCj, mrows)
                #assert max(GCj) <= ncols, 'GCi=%s GCj=%s ncols=%s' % (GCi, GCj, ncols)
                GCi = np.array(GCi, dtype='int32') - 1
                GCj = np.array(GCj, dtype='int32') - 1
                try:
                    # we subtract 1 to the indicides to account for Fortran
                    #    huh??? we don't...
                    if dtype == '???':
                        matrix = None
                        self.log.warning('what is the dtype?')
                    elif tout in [1, 2]:
                        real_array = np.array(reals, dtype=dtype)
                        matrix = scipy.sparse.coo_matrix(
                            (real_array, (GCi, GCj)),
                            shape=(mrows, ncols), dtype=dtype)
                        matrix = matrix.todense()
                        #self.log.info('created %s' % self.table_name)
                    elif tout in [3, 4]:
                        real_array = np.array(reals, dtype=dtype)
                        nvalues_matrix = real_array.shape[0] // 2
                        real_complex = real_array.reshape((nvalues_matrix, 2))
                        real_imag = real_complex[:, 0] + real_complex[:, 1]*1j
                        if self.binary_debug:
                            #self.binary_debug.write('reals = %s' % real_complex[:, 0])
                            #self.binary_debug.write('imags = %s' % real_complex[:, 1])
                            self.binary_debug.write('real_imag = %s' % real_imag)
                        matrix = scipy.sparse.coo_matrix(
                            (real_imag, (GCi, GCj)),
                            shape=(mrows, ncols), dtype=dtype)
                        msg = 'created %s...verify the complex matrix' % self.table_name
                        self.log.warning(msg)
                        #raise RuntimeError(msg)
                    else:
                        raise RuntimeError('this should never happen')
                except ValueError:
                    self.log.warning('shape=(%s, %s)' % (mrows, ncols))
                    self.log.warning('cant make a coo/sparse matrix...trying dense')

                    if dtype == '???':
                        matrix = None
                        self.log.warning('what is the dtype?')
                    else:
                        real_array = np.array(reals, dtype=dtype)
                        self.log.debug('shape=%s mrows=%s ncols=%s' % (
                            str(real_array.shape), mrows, ncols))
                        if len(reals) == mrows * ncols:
                            real_array = real_array.reshape(mrows, ncols)
                            self.log.info('created %s' % self.table_name)
                        else:
                            self.log.warning('cant reshape because invalid sizes : created %s' %
                                             self.table_name)

                        matrix = real_array

                m.data = matrix
                if matrix is not None:
                    self.matrices[table_name.decode('utf-8')] = m
                #nvalues = self.get_marker1(rewind=True)
                return
            itable -= 1
            niter += 1
        raise RuntimeError('MaxIteration: this should never happen; n=%s' % niter_max)

def grids_comp_array_to_index(grids1, comps1, grids2, comps2,
                              make_matrix_symmetric):
    """maps the dofs"""
    #from pyNastran.utils.mathematics import unique2d
    ai = np.vstack([grids1, comps1]).T
    bi = np.vstack([grids2, comps2]).T
    #print('grids2 =', grids2)
    #print('comps2 =', comps2)
    from itertools import count
    #c = np.vstack([a, b])
    #assert c.shape[1] == 2, c.shape
    #urows = unique2d(c)
    #urows = c

    nid_comp_to_dof_index = {}
    j = 0
    a_keys = set()
    for nid_dof in ai:
        #nid_dof = (int(nid), int(dof))
        nid_dof = tuple(nid_dof)
        if nid_dof not in a_keys:
            a_keys.add(nid_dof)
            if nid_dof not in nid_comp_to_dof_index:
                nid_comp_to_dof_index[nid_dof] = j
                j += 1
    nja = len(a_keys)
    del a_keys

    b_keys = set()
    for nid_dof in bi:
        nid_dof = tuple(nid_dof)
        if nid_dof not in b_keys:
            b_keys.add(nid_dof)
        if nid_dof not in nid_comp_to_dof_index:
            nid_comp_to_dof_index[nid_dof] = j
            j += 1
    njb = len(b_keys)
    del b_keys


    nj = len(nid_comp_to_dof_index)
    if make_matrix_symmetric:
        ja = np.zeros(nj, dtype='int32')
        for i, nid_dof in zip(count(), ai):
            j[i] = nid_comp_to_dof_index[tuple(nid_dof)]
        for i, nid_dof in zip(count(), bi):
            j[i] = nid_comp_to_dof_index[tuple(nid_dof)]
        return j, j, nj, nj, nj
    else:
        ja = np.zeros(grids1.shape, dtype='int32')
        for i, nid_dof in zip(count(), ai.tolist()):
            ja[i] = nid_comp_to_dof_index[tuple(nid_dof)]

        jb = np.zeros(grids2.shape, dtype='int32')
        for i, nid_dof in zip(count(), bi.tolist()):
            jb[i] = nid_comp_to_dof_index[tuple(nid_dof)]

        return ja, jb, nja, njb, nj


