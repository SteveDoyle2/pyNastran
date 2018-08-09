"""
Defines various tables that don't fit in other sections:
  - OP2Reader
    - read_cmodeext(self)
    - read_cmodeext_helper(self)
    - read_aemonpt(self)
    - read_monitor(self)
    - read_r1tabrg(self, data, ndata)
    - read_hisadd(self)

    - read_cstm()
    - read_dit()
    - read_extdb()
    - read_fol()
    - read_frl()
    - read_gpl()
    - read_ibulk()
    - read_intmod()
    - read_meff()
    - read_omm2()
    - read_sdf()
    - read_tol()
    - _skip_pcompts(self)
    - _read_pcompts(self)

  - Others
    - read_markers()
    - _skip_subtables()
    - _skip_table_helper()
    - _skip_matrix_mat()
   - _print_month(month, day, year, zero, one)
"""
from __future__ import print_function
#from copy import deepcopy
from struct import unpack
from six import b #iteritems
import numpy as np
from pyNastran.f06.errors import FatalError
from pyNastran.op2.errors import FortranMarkerError #, SortCodeError
from pyNastran.op2.tables.design_response import (
    WeightResponse, StressResponse, StrainResponse, ForceResponse,
    FlutterResponse, Convergence)


class OP2Reader(object):
    """Stores methods that aren't useful to an end user"""
    def __init__(self, op2):
        self.op2 = op2
        self.mapped_tables = {
            b'GPL' : self.read_gpl,
            #b'MEFF' : self.read_meff,
            b'INTMOD' : self.read_intmod,
            b'HISADD' : self.read_hisadd,
            b'EXTDB' : self.read_extdb,
            b'OMM2' : self.read_omm2,
            b'TOL' : self.read_tol,
            b'PCOMPTS' : self._read_pcompts,
            b'MONITOR' : self._read_monitor,
            b'AEMONPT' : self._read_aemonpt,
            b'FOL' : self.read_fol,
            b'SDF' : self.read_sdf,
            b'IBULK' : self.read_ibulk,
            b'CDDATA' : self.read_ibulk,
            b'CMODEXT' : self._read_cmodext,
            b'CSTM' : self.read_cstm,
        }
        #self.op2_skip = OP2Skip(op2)

    def _read_aemonpt(self):
        """reads the AEMONPT table"""
        op2 = self.op2
        #self.log.debug("table_name = %r" % op2.table_name)
        table_name = self.op2_reader._read_table_name(rewind=False)

        #print('-----------------------')
        #print('record 1')
        self.read_markers([-1])
        data = self._read_record()
        #self.show_data(data)
        if op2.read_mode == 2:
            a, bi, c, d, e, f, g = unpack(self._endian + b'7i', data)
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
        #if op2.read_mode == 2:
        word, = unpack(self._endian + b'8s', data)
        assert word == b'AECFMON ', word
        #op2.show_data(data)
        #print('-----------------------')
        #print('record 3')
        self.read_markers([-3, 1, 0])
        data = self._read_record()
        #op2.show_data(data)

        if op2.read_mode == 2:
            ndata = len(data)
            assert ndata == 108, ndata
            n = 8 + 56 + 20 + 12 + 12
            out = unpack(self._endian + b'8s 56s 5i 4s 8s 3i', data[:n])
            (aero, name, comps, cp, bi, c, d, coeff, word, e, f, g) = out
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
        #op2.show_data(data)
        self.read_markers([0])

        #print('-----------------------')
        #print('end')
        #op2.show(200)
        #aaa

    def _read_monitor(self):
        """reads the MONITOR table"""
        op2 = self.op2
        self.log.debug("table_name = %r" % op2.table_name)
        table_name = self._read_table_name(rewind=False)

        #print('-----------------------')
        #print('record 1')
        self.read_markers([-1])
        data = self._read_record()
        #self.show_data(data)
        if op2.read_mode == 2:
            a, bi, c, d, e, f, g = unpack(self._endian + b'7i', data)
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
        if op2.read_mode == 2:
            word, = op2.struct_8s.unpack(data)
            assert word == b'STCFMON ', word
        #op2.show_data(data)
        #print('-----------------------')
        #print('record 3')
        self.read_markers([-3, 1, 0])
        data = self._read_record()
        #self.show_data(data[96:108])

        if op2.read_mode == 2:
            ndata = len(data)
            assert ndata == 108, ndata
            (aero, name, comps, cp, x, y, z, coeff, word, column, cd,
             ind_dof) = unpack(self._endian + b'8s 56s 2i 3f 4s 8s 3i', data[:108])
            #print('aero=%r' % aero)
            #print('name=%r' % name)
            #print('comps=%s cp=%s (x, y, z)=(%s, %s, %s)' % (comps, cp, x, y, z))
            #print('coeff=%r' % coeff)
            #print('word=%r (column, cd, ind_dof)=(%s, %s, %s)' % (word, column, cd, ind_dof))
            assert cp == 2, cp
            assert x == 0.0, x
            assert y == 0.0, y
            assert d == 0.0, z
            assert column == 1, column
            assert cd == 2, cd
            assert ind_dof == 0, ind_dof
            op2.monitor_data = [{
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
        #op2.show_data(data)
        self.read_markers([0])

        #print('-----------------------')
        #print('end')
        #op2.show(200)

    def _read_cmodext(self):
        r"""
        fails if a streaming block???:
         - nx_spike\mnf16_0.op2
        """
        op2 = self.op2
        op2.table_name = self._read_table_name(rewind=False)
        self.log.debug('table_name = %r' % op2.table_name)
        if self.is_debug_file:
            self.binary_debug.write('_read_geom_table - %s\n' % op2.table_name)
        self.read_markers([-1])
        if self.is_debug_file:
            self.binary_debug.write('---markers = [-1]---\n')
        data = self._read_record()

        markers = self.get_nmarkers(1, rewind=True)
        if self.is_debug_file:
            self.binary_debug.write('---marker0 = %s---\n' % markers)

        marker = -2
        markers = self.read_markers([marker, 1, 0])

        data = self._read_record()
        table_name, oneseventy_a, oneseventy_b = unpack('8sii', data)
        assert oneseventy_a == 170, oneseventy_a
        assert oneseventy_b == 170, oneseventy_b
        print('170*4 =', 170*4)
        #op2.show_data(data)
        marker -= 1
        marker = self._read_cmodext_helper(marker) # -3
        marker = self._read_cmodext_helper(marker)
        marker = self._read_cmodext_helper(marker)
        marker = self._read_cmodext_helper(marker)
        marker = self._read_cmodext_helper(marker)
        print('table8')
        marker = self._read_cmodext_helper(marker, debug=True)
        op2.show_ndata(100)

    def _read_cmodext_helper(self, marker_orig, debug=False):
        op2 = self.op2
        marker = marker_orig
        #markers = self.read_nmarkers([marker, 1, 1]) # -3

        if debug:
            op2.show_ndata(100)
        markers = self.get_nmarkers(3, rewind=False)
        assert markers == [marker_orig, 1, 1], markers
        print('markers =', markers)

        #marker = self.get_nmarkers(1, rewind=False, macro_rewind=False)[0]
        val_old = 0
        if debug:
            print('-----------------------------')
        i = 0
        #icheck = 7
        while 1:
            #print('i = %i' % i)
            marker = self.get_nmarkers(1, rewind=False, macro_rewind=False)[0]
            if marker != 6:
                print('marker = %s' % marker)

            assert marker == 6, marker
            data = self.read_block()
            val = unpack('i', data[:4])[0]
            if debug:
                print('val=%s delta=%s' % (val, val - val_old))
                op2.show_data(data, types='ifs')
            assert len(data) > 4
            #print('i=%s val=%s delta=%s' % (i, val, val - val_old))
            val_old = val

            marker2 = self.get_nmarkers(1, rewind=True, macro_rewind=False)[0]
            #print(marker2)
            if marker2 == 696:
                break
            i += 1
        if debug:
            print('----------------------------------------')

        marker = self.get_nmarkers(1, rewind=False, macro_rewind=False)[0]
        if debug:
            print('****marker = %s' % marker)
        assert marker == 696, marker
        data = self.read_block()
        #self.show_data(data)

        marker = self.get_nmarkers(1, rewind=True, macro_rewind=False)[0]
        assert marker == (marker_orig - 1), marker

        if debug:
            op2.show_ndata(200)
        return marker

        #data = self._read_record()
        #marker -= 1
        #op2.show_ndata(100)

        ##marker -= 1
        ##marker_end = op2.get_marker1(rewind=False)

    def read_cstm(self):
        """
        Reads the CSTM table, which defines the transform from global to basic.

        Returns 14-column matrix 2-d array of the CSTM data:
        ::
          [
           [ id1 type xo yo zo T(1,1:3) T(2,1:3) T(3,1:3) ]
           [ id2 type xo yo zo T(1,1:3) T(2,1:3) T(3,1:3) ]
           ...
          ]

        T is transformation from local to basic for the coordinate system.
        """
        op2 = self.op2
        table_name = self._read_table_name(rewind=False)
        self.read_markers([-1])
        data = self._read_record()
        assert len(data) == 28, len(data)

        self.read_markers([-2, 1, 0])
        data = self._read_record()
        assert len(data) == 8, len(data)

        self.read_markers([-3, 1, 0])
        data = self._read_record()
        ints = np.frombuffer(data, dtype='int32')
        floats = np.frombuffer(data, dtype='float32')
        nints = len(ints)
        assert nints % 14 == 0, 'nints=%s' % (nints)
        ncstm = nints // 14
        ints = ints.reshape(ncstm, 14)[:, :2]
        floats = floats.reshape(ncstm, 14)[:, 2:]
        #assert ncstm == 1, 'ncoords = %s' % ncstm
        #print(self.coords)
        coord_type_map = {
            1 : 'CORD2R',
            2 : '???',
        }
        for i, coord in enumerate(ints):
            cid = ints[i, 0]
            coord_type_int = ints[i, 1]
        if coord_type_int in coord_type_map:
            coord_type = coord_type_map[coord_type_int]
        else:  # pragma: no cover
            msg = 'cid=%s coord_type_int=%s is not supported\n' % (cid, coord_type_int)
            if hasattr(self, 'coords'):
                print(op2.coords)
            raise RuntimeError(msg)

        #print(self.coords)
        #print('cid = ', cid)
        #print('coord_type = ', coord_type)
        #print('myints =', ints)
        #print('floats =', floats)

        self.read_markers([-4, 1, 0, 0])

    def _read_dit(self):
        """
        Reads the DIT table (poorly).
        The DIT table stores information about table cards
        (e.g. TABLED1, TABLEM1).

        """
        op2 = self.op2
        table_name = self._read_table_name(rewind=False)
        self.read_markers([-1])
        data = self._read_record()

        self.read_markers([-2, 1, 0])
        data = self._read_record()
        table_name, = op2.struct_8s.unpack(data)

        self.read_markers([-3, 1, 0])
        data = self._read_record()

        self.read_markers([-4, 1, 0])
        data = self._read_record()

        self.read_markers([-5, 1, 0])

        itable = -6
        while 1:
            markers = self.get_nmarkers(1, rewind=True)
            if markers == [0]:
                break
            data = self._read_record()
            self.read_markers([itable, 1, 0])
            itable -= 1

        #self.show(100)
        self.read_markers([0])

    def read_extdb(self):
        r"""
        fails if a streaming block:
         - nx_spike\extse04c_0.op2
        """
        op2 = self.op2
        op2.table_name = self._read_table_name(rewind=False)
        self.log.debug('table_name = %r' % op2.table_name)
        if self.is_debug_file:
            self.binary_debug.write('_read_geom_table - %s\n' % op2.table_name)
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
                op2.show_ndata(100)
                raise
            nfields = op2.get_marker1(rewind=True)
            if nfields > 0:
                data = self._read_record()
                #self.show_data(data, types='s', endian=None)
            elif nfields == 0:
                #self.show_ndata(100, types='ifs')
                break
            else:
                raise RuntimeError('nfields=%s' % nfields)
            marker -= 1
        marker_end = op2.get_marker1(rewind=False)

    def read_fol(self):
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
        op2 = self.op2
        op2.log.debug("table_name = %r" % op2.table_name)
        op2.table_name = self._read_table_name(rewind=False)
        self.read_markers([-1])
        data = self._read_record()
        self.read_markers([-2, 1, 0])
        data = self._read_record()
        ndata = len(data)
        subtable_name_raw, = op2.struct_8s.unpack(data[:8])
        subtable_name = subtable_name_raw.strip()
        assert subtable_name == b'FOL', 'subtable_name=%r' % subtable_name

        nfloats = (ndata - 8) // 4
        assert nfloats * 4 == (ndata - 8)
        fmt = b(op2._uendian + '%sf' % nfloats)
        freqs = np.array(list(unpack(fmt, data[8:])), dtype='float32')
        op2._frequencies = freqs
        if self.is_debug_file:
            self.binary_debug.write('  recordi = [%r, freqs]\n'  % (subtable_name_raw))
            self.binary_debug.write('  subtable_name=%r\n' % subtable_name)
            self.binary_debug.write('  freqs = %s' % freqs)
        op2._read_subtables()

    def read_frl(self):
        """reads the FRL (Frequency Response List) table"""
        op2 = self.op2
        #op2.log.debug("table_name = %r" % op2.table_name)
        #op2.table_name = self._read_table_name(rewind=False)
        #self.read_markers([-1])
        #data = self._read_record()

        #self.read_markers([-2, 1, 0])
        #data = self._read_record()
        op2._skip_table(op2.table_name)

    def read_gpl(self):
        """reads the GPL table (grid point list?)"""
        op2 = self.op2
        op2.table_name = self._read_table_name(rewind=False)
        self.log.debug('table_name = %r' % op2.table_name)
        if self.is_debug_file:
            self.binary_debug.write('_read_geom_table - %s\n' % op2.table_name)
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
                markers = op2.get_nmarkers(1, rewind=False)
                break
            data = self._read_record()
            #op2.show_data(data, 'i')
            n -= 1
            markers = self.get_nmarkers(1, rewind=True)

    def read_hisadd(self):
        """optimization history (SOL200) table"""
        op2 = self.op2
        op2.table_name = self._read_table_name(rewind=False)

        if op2.read_mode == 1:
            self.read_markers([-1])
            self._skip_record()
            self.read_markers([-2, 1, 0])
            self._skip_record()
            self.read_markers([-3, 1, 0])

            if op2.responses.convergence_data is None:
                data = self._read_record()
                ndvs = len(data) // 4 - 7
                op2.responses.convergence_data = Convergence(ndvs)
            else:
                self._skip_record()
                op2.responses.convergence_data.n += 1

            self.read_markers([-4, 1, 0, 0])
            return

        if self.is_debug_file:
            self.binary_debug.write('_read_geom_table - %s\n' % op2.table_name)
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
        #op2.show_data(data)

        #self.log.info('----marker3----')
        self.read_markers([-3, 1, 0])
        data = self._read_record()

        (design_iter, iconvergence, conv_result, obj_intial, obj_final,
         constraint_max, row_constraint_max) = unpack(self._endian + b'3i3fi', data[:28])
        if iconvergence == 1:
            iconvergence = 'soft'
        elif iconvergence == 2:
            iconvergence = 'hard'
        elif iconvergence == 6:
            self.log.warning('HISADD iconverge=6')
            iconvergence = '???'
        else:
            msg = 'iconvergence=%s\n' % iconvergence
            op2.show_data(data, types='ifs', endian=None)
            raise NotImplementedError(msg)

        if conv_result == 0:
            conv_result = 'no'
        elif conv_result == 1:
            conv_result = 'soft'
        elif conv_result == 2:
            conv_result = 'hard'
        elif conv_result in [3, 4]:
            #self.log.warning('HISADD conv_result=%s' % conv_result)
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

        op2.responses.convergence_data.append(
            design_iter, iconvergence, conv_result, obj_intial,
            obj_final, constraint_max, row_constraint_max, desvar_values)
        self.read_markers([-4, 1, 0, 0])

    def read_ibulk(self):
        op2 = self.op2
        op2.table_name = self._read_table_name(rewind=False)
        op2.log.debug('table_name = %r' % op2.table_name)
        if self.is_debug_file:
            self.binary_debug.write('_read_geom_table - %s\n' % op2.table_name)
        self.read_markers([-1])
        if self.is_debug_file:
            self.binary_debug.write('---markers = [-1]---\n')
        data = self._read_record()

        markers = self.get_nmarkers(1, rewind=True)
        marker = -2
        while 1:
            self.read_markers([marker, 1, 0])
            nfields = op2.get_marker1(rewind=True)
            if nfields > 0:
                data = self._read_record()
                #self.show_data(data, types='s', endian=None)
            elif nfields == 0:
                #op2.show_ndata(100, types='ifs')
                break
            else:
                raise RuntimeError('nfields=%s' % nfields)
            marker -= 1
        marker_end = op2.get_marker1(rewind=False)

    def read_omm2(self):
        """reads the OMM2 table"""
        op2 = self.op2
        op2.log.debug("table_name = %r" % op2.table_name)
        op2.table_name = self._read_table_name(rewind=False)
        self.read_markers([-1])
        data = self._read_record()

        self.read_markers([-2, 1, 0])
        data = self._read_record()
        if len(data) == 28:
            subtable_name, month, day, year, zero, one = unpack(self._endian + b'8s5i', data)
            if self.is_debug_file:
                self.binary_debug.write('  recordi = [%r, %i, %i, %i, %i, %i]\n'  % (
                    subtable_name, month, day, year, zero, one))
                self.binary_debug.write('  subtable_name=%r\n' % subtable_name)
            self._print_month(month, day, year, zero, one)
        else:
            raise NotImplementedError(op2.show_data(data))
        op2._read_subtables()

    def _read_pcompts(self):
        """
        Reads the PCOMPTS table (poorly).
        The PCOMPTS table stores information about the PCOMP cards???
        """
        self._skip_pcompts()
        return
        #if op2.read_mode == 1:
            #return
        #op2.log.debug("table_name = %r" % op2.table_name)
        #table_name = self._read_table_name(rewind=False)

        #self.read_markers([-1])
        #data = self._read_record()

        #self.read_markers([-2, 1, 0])
        #data = self._read_record()
        #table_name, = op2.struct_8s.unpack(data)
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

    def _skip_pcompts(self):
        """
        Reads the PCOMPTS table (poorly).
        The PCOMPTS table stores information about the PCOMP cards???
        """
        op2 = self.op2
        op2.log.debug("table_name = %r" % op2.table_name)
        table_name = self._read_table_name(rewind=False)

        self.read_markers([-1])
        data = self._skip_record()

        self.read_markers([-2, 1, 0])
        data = self._skip_record()
        #table_name, = op2.struct_8s.unpack(data)

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

    def read_meff(self):
        """reads the MEFF table"""
        op2 = self.op2
        op2.table_name = self._read_table_name(rewind=False)
        op2.log.debug('table_name = %r' % op2.table_name)
        if self.is_debug_file:
            self.binary_debug.write('_read_geom_table - %s\n' % op2.table_name)
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
            data = op2.f.read(nbytes)
            op2.n += nbytes
        n = -9
        self.read_markers([n, 1, 0, 0])

    def read_intmod(self):
        """reads the INTMOD table"""
        op2 = self.op2
        op2.table_name = self._read_table_name(rewind=False)
        #op2.log.debug('table_name = %r' % op2.table_name)
        if self.is_debug_file:
            self.binary_debug.write('_read_geom_table - %s\n' % op2.table_name)
        self.read_markers([-1])
        if self.is_debug_file:
            self.binary_debug.write('---markers = [-1]---\n')
        data = self._read_record()
        #print('intmod data1')
        #op2.show_data(data)

        markers = self.get_nmarkers(1, rewind=True)
        if self.is_debug_file:
            self.binary_debug.write('---marker0 = %s---\n' % markers)
        self.read_markers([-2, 1, 0])
        data = self._read_record()
        #print('intmod data2')
        #op2.show_data(data)

        for n in [-3, -4, -5, -6, -7, -8,]:
            self.read_markers([n, 1, 1])
            markers = self.get_nmarkers(1, rewind=False)
            #print('markers =', markers)
            nbytes = markers[0]*4 + 12
            data = op2.f.read(nbytes)
            #print('intmod data%i' % n)
            #op2.show_data(data)
            op2.n += nbytes

        n = -9
        self.read_markers([n, 1, 0, 0])
        #op2.show(50)
        #raise NotImplementedError(op2.table_name)

    def read_r1tabrg(self, data, ndata):
        """
        Design Responses:
          - Weight
          - Flutter Speed
          - Stress
          - Strain
          - Displacement
        """
        op2 = self.op2
        responses = op2.responses
        if op2._table4_count == 0:
            op2._count += 1
        op2._table4_count += 1

        #if op2._table4_count == 0:
            #op2._count += 1
        #op2._table4_count += 1

        if op2.read_mode == 1:
            assert data is not None, data
            assert len(data) > 12, len(data)
            response_type, = op2.struct_i.unpack(data[8:12])
            #assert response_type in [1, 6, 10, 84], response_type
            if response_type == 1:
                if responses.weight_response is None:
                    responses.weight_response = WeightResponse()
                else:
                    responses.weight_response.n += 1
            elif response_type == 4:
                #TYPE =4 EIGN or FREQ
                #8 MODE I Mode number
                #9 APRX I Approximation code
                pass
            elif response_type == 5:
                #TYPE =5 DISP
                #8 COMP I Displacement component
                #9 UNDEF None
                #10 GRID I Grid identification number
                pass
            elif response_type == 6:
                if responses.stress_response is None:
                    responses.stress_response = StressResponse()
                else:
                    responses.stress_response.n += 1
            elif response_type == 7:
                if responses.strain_response is None:
                    responses.strain_response = StrainResponse()
                else:
                    responses.strain_response.n += 1

            elif response_type == 8:
                if responses.force_response is None:
                    responses.force_response = ForceResponse()
                else:
                    responses.force_response.n += 1
            elif response_type == 15:
                # CEIG
                #8 MODE I Mode number
                #9 ICODE I 1: Real component or 2: Imaginary component
                pass
            elif response_type == 84:
                if responses.flutter_response is None:
                    responses.flutter_response = FlutterResponse()
                else:
                    responses.flutter_response.n += 1
            return ndata
            #else: # response not added...
                #pass

        read_r1tabrg = True
        if read_r1tabrg:
            #self.show_data(data, types='ifs', endian=None)
            out = unpack(self._endian + b'iii 8s iiii i iiiii', data)
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
                responses.weight_response.add_from_op2(out, self.log)
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
                #msg = ('STRESS - response_type=%r label=%r region=%s subcase=%s '
                       #'stress_code=%s pid=%s' % (
                           #response_type, response_label, region, subcase,
                           #stress_code, pid))
                responses.stress_response.append(
                    internal_id, dresp_id, response_label, region,
                    subcase, type_flag, seid,
                    stress_code, pid)

            elif response_type == 7:  # STRAIN
                strain_code = out[6]
                pid = out[8]
                responses.strain_response.append(
                    internal_id, dresp_id, response_label, region,
                    subcase, type_flag, seid,
                    strain_code, pid)
            elif response_type == 8:  # FORCE
                #print('internal_id=%s dresp_id=%s response_type=%s response_label=%s region=%s subcase=%s type_flag=%s seid=%s' % (
                    #internal_id, dresp_id, response_type, response_label, region, subcase, type_flag, seid
                #))
                force_code = out[6]
                pid = out[8]
                #msg = 'FORCE - label=%r region=%s subcase=%s force_code=%s pid=%s' % (
                    #response_label, region, subcase, force_code, pid)
                #print(msg)
                #print(out)
                responses.force_response.append(
                    internal_id, dresp_id, response_label, region,
                    subcase, type_flag, seid,
                    force_code, pid)

            elif response_type == 10:  # CSTRESS
                stress_code = out[6]
                ply = out[7]
                pid = out[8]  # is this element id?
                #msg = 'CSTRESS - label=%r region=%s subcase=%s stress_code=%s ply=%s pid=%s' % (
                    #response_label, region, subcase, stress_code, ply, pid)
                #print(msg)
            #elif response_type == 10:  # CSTRAIN
                #pass
            elif response_type == 24:  # FRSTRE
                #8 ICODE I Stress item code
                #9 UNDEF None
                #10 ELID I Element identification number
                #11 FREQ RS Frequency
                #12 IFLAG I Integrated response flag. See Remark 20 of DRESP1.
                #Value is -1 to -6, for SUM, AVG, SSQ,
                pass
            elif response_type == 28:  # RMSACCL
                #8 COMP I RMS Acceleration component
                #9 RANDPS I RANDPS entry identification number
                #10 GRID I Grid identification number
                #11 DMFREQ RS Dummy frequency for internal use
                pass
            elif response_type == 84:
                # FLUTTER  (iii, label, mode, (Ma, V, rho), flutter_id, fff)
                out = unpack(self._endian + b'iii 8s iii fff i fff', data)
                mode = out[6]
                mach = out[7]
                velocity = out[8]
                density = out[9]
                flutter_id = out[10]
                msg = ('FLUTTER - _count=%s label=%r region=%s subcase=%s mode=%s '
                       'mach=%s velocity=%s density=%s flutter_id=%s' % (
                           op2._count, response_label, region, subcase, mode,
                           mach, velocity, density, flutter_id))
                responses.flutter_response.append(
                    internal_id, dresp_id, response_label, region,
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

    def read_sdf(self):
        """reads the SDF table"""
        op2 = self.op2
        op2.log.debug("table_name = %r" % op2.table_name)
        op2.table_name = self._read_table_name(rewind=False)
        self.read_markers([-1])
        data = self._read_record()

        self.read_markers([-2, 1, 0])
        data, ndata = self._read_record_ndata()
        if ndata == 16:
            subtable_name, dummy_a, dummy_b = unpack(self._endian + b'8sii', data)
            if self.is_debug_file:
                self.binary_debug.write('  recordi = [%r, %i, %i]\n'  % (
                    subtable_name, dummy_a, dummy_b))
                self.binary_debug.write('  subtable_name=%r\n' % subtable_name)
                assert dummy_a == 170, dummy_a
                assert dummy_b == 170, dummy_b
        else:
            strings, ints, floats = op2.show_data(data)
            msg = 'len(data) = %i\n' % ndata
            msg += 'strings  = %r\n' % strings
            msg += 'ints     = %r\n' % str(ints)
            msg += 'floats   = %r' % str(floats)
            raise NotImplementedError(msg)

        self.read_markers([-3, 1, 1])

        markers0 = self.get_nmarkers(1, rewind=False)
        record = op2.read_block()

        #data = self._read_record()
        self.read_markers([-4, 1, 0, 0])
        #self._read_subtables()

    def read_tol(self):
        """
        This is probably broken for MSC Nastran

        TOL
        ---
        -2 - nitimes?
        -3 - list of times?

        """
        table_name = self._read_table_name(rewind=False, stop_on_failure=True)
        self.read_markers([-1])
        data = self._read_record()
        #op2.show_data(data)

        self.read_markers([-2, 1, 0])
        #op2.show_ndata(440, types='if')
        data = self._read_record()
        #print('----')
        self.read_markers([-3, 1, 0])
        #op2.show_ndata(440, types='if')
        #print('----')
        self.read_markers([0])
        #data = self._read_record()


        #op2.show_ndata(440, types='ifs')

        #op2.show_data(data)
        #aaaa

    def _skip_matrix_mat(self):
        """
        Reads a matrix in "standard" form.

        See
        ---
        read_matrix_mat
        """
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

    #---------------------------------------------------------------------------
    def get_nmarkers(self, n, rewind=True, macro_rewind=False):
        """
        Gets n markers, so if n=2, it will get 2 markers.

        Parameters
        ----------
        n : int
            number of markers to get
        rewind : bool
            should the file be returned to the starting point

        Returns
        -------
        markers : List[int]
            list of [1, 2, 3, ...] markers

        """
        op2 = self.op2
        ni = op2.n
        markers = []
        for i in range(n):
            data = self.read_block()
            marker, = op2.struct_i.unpack(data)
            markers.append(marker)
        if rewind:
            op2.n = ni
            op2.f.seek(op2.n)
            #for i in range(n):
                #self.binary_debug.write('get_nmarkers- [4, %i, 4]; macro_rewind=%s\n' % (
                    #i, macro_rewind or rewind))
        else:
            #if not macro_rewind:
            if self.is_debug_file:
                for i in range(n):
                    self.binary_debug.write('get_nmarkers- [4, %i, 4]; macro_rewind=%s\n' % (
                        i, macro_rewind or rewind))
        return markers

    def read_markers(self, markers, macro_rewind=True):
        """
        Gets specified markers, where a marker has the form of [4, value, 4].
        The "marker" corresponds to the value, so 3 markers takes up 9 integers.
        These are used to indicate position in the file as well as the number
        of bytes to read.

        Because we're checking the markers vs. what we expect, we just throw
        the data away.

        Parameters
        ----------
        markers : List[int]
            markers to get; markers = [-10, 1]

        Raises
        ------
        FortranMarkerError
            if the expected table number is not found

        """
        op2 = self.op2
        for i, marker in enumerate(markers):
            data = op2.read_block()
            imarker, = op2.struct_i.unpack(data)
            if marker != imarker:
                import os
                msg = 'marker=%r imarker=%r; markers=%s; i=%s; table_name=%r; iloc=%s/%s' % (
                    marker, imarker, markers, i, op2.table_name,
                    op2.f.tell(), os.path.getsize(op2.op2_filename))
                raise FortranMarkerError(msg)
            if self.is_debug_file:
                self.binary_debug.write('  read_markers -> [4, %i, 4]\n' % marker)

    def _skip_table(self, table_name):
        """bypasses the next table as quickly as possible"""
        if table_name in ['DIT', 'DITS']:  # tables
            self._read_dit()
        elif table_name in ['PCOMPTS']:
            self._read_pcompts()
        else:
            self._skip_table_helper()

    def _print_month(self, month, day, year, zero, one):
        """
        Creates the self.date attribute from the 2-digit year.

        Parameters
        ----------
        month : int
            the month (integer <= 12)
        day :  int
            the day (integer <= 31)
        year : int
            the day (integer <= 99)
        zero : int
            a dummy integer (???)
        one : int
            a dummy integer (???)

        """
        month, day, year = self.op2._set_op2_date(month, day, year)

        #self.log.debug("%s/%s/%4i zero=%s one=%s" % (month, day, year, zero, one))
        #if self.is_debug_file:
        if self.is_debug_file:
            self.binary_debug.write('  [subtable_name, month=%i, day=%i, year=%i, '
                                    'zero=%i, one=%i]\n\n' % (month, day, year, zero, one))
        #assert zero == 0, zero  # is this the RTABLE indicator???
        assert one in [0, 1], one  # 0, 50

    #----------------------------------------------------------------------------------------
    def _read_record(self, stream=False, debug=True, macro_rewind=False):
        """interface to the op2 object"""
        return self.op2._read_record(stream=stream, debug=debug, macro_rewind=macro_rewind)

    def _read_record_ndata(self, stream=False, debug=True, macro_rewind=False):
        """reads a record and the length of the record"""
        op2 = self.op2
        markers0 = self.get_nmarkers(1, rewind=False, macro_rewind=macro_rewind)
        if self.is_debug_file and debug:
            self.binary_debug.write('read_record - marker = [4, %i, 4]; macro_rewind=%s\n' % (
                markers0[0], macro_rewind))
        record, nrecord = self._read_block_ndata()

        if self.is_debug_file and debug:
            msg = 'read_record - record = [%i, recordi, %i]; macro_rewind=%s\n' % (
                nrecord, nrecord, macro_rewind)
            self.binary_debug.write(msg)
        if markers0[0]*4 != len(record):
            raise FortranMarkerError('markers0=%s*4 len(record)=%s; table_name=%r' % (
                markers0[0]*4, len(record), op2.table_name))

        markers1 = self.get_nmarkers(1, rewind=True)
        if markers1[0] > 0:
            nloop = 0
            records = [record]
            while markers1[0] > 0:
                markers1 = self.get_nmarkers(1, rewind=False)
                if self.is_debug_file and debug:
                    self.binary_debug.write('read_record - markers1 = [4, %i, 4]\n' % markers1[0])
                recordi, nrecordi = self._read_block_ndata()
                nrecord += nrecordi
                records.append(recordi)
                #record += recordi
                markers1 = self.get_nmarkers(1, rewind=True)
                if self.is_debug_file and debug:
                    self.binary_debug.write('read_record - markers1 = [4, %i, 4]\n' % markers1[0])
                nloop += 1

            # if nloop == 0:
                # record = records[0]
            # elif nloop == 1:
                # record = records[0] + records[1]
            # else:
            record = b''.join(records)
        return record, nrecord

    def _read_block_ndata(self):
        """
        Reads a block following a pattern of:
            [nbytes, data, nbytes]

        Returns
        -------
        data : bytes
            the data in binary
        ndata : int
            len(data)

        """
        op2 = self.op2
        data = op2.f.read(4)
        ndata, = op2.struct_i.unpack(data)

        data_out = op2.f.read(ndata)
        data = op2.f.read(4)
        op2.n += 8 + ndata
        return data_out, ndata

    #------------------------------------------------------------------
    def _read_table_name(self, rewind=False, stop_on_failure=True):
        """Reads the next OP2 table name (e.g. OUG1, OES1X1)"""
        op2 = self.op2
        table_name = None
        data = None
        if self.is_debug_file:
            self.binary_debug.write('_read_table_name - rewind=%s\n' % rewind)
        ni = op2.n
        structi = op2.struct_8s
        if stop_on_failure:
            data = self._read_record(debug=False, macro_rewind=rewind)
            table_name, = structi.unpack(data)
            if self.is_debug_file and not rewind:
                self.binary_debug.write('marker = [4, 2, 4]\n')
                self.binary_debug.write('table_header = [8, %r, 8]\n\n' % table_name)
            table_name = table_name.strip()
        else:
            try:
                data = self._read_record(macro_rewind=rewind)
                table_name, = structi.unpack(data)
                table_name = table_name.strip()
            except:
                # we're done reading
                op2.n = ni
                op2.f.seek(op2.n)

                try:
                    # we have a trailing 0 marker
                    self.read_markers([0], macro_rewind=rewind)
                except:
                    # if we hit this block, we have a FATAL error
                    if not op2._nastran_format.lower().startswith('imat') and op2.post != -4:
                        op2.f.seek(op2.n)
                        op2.show(1000)
                        raise FatalError('There was a Nastran FATAL Error.  '
                                         'Check the F06.\nlast table=%r; post=%s' % (
                                             op2.table_name, op2.post))
                table_name = None

                # we're done reading, so we're going to ignore the rewind
                rewind = False

        if rewind:
            op2.n = ni
            op2.f.seek(op2.n)
        return table_name


    def read_block(self):
        """interface to the op2 object"""
        return self.op2.read_block()

    def get_marker1(self, rewind=True, macro_rewind=False):
        """
        Gets 1 marker
        See get_n_markers(...)

        Parameters
        ----------
        rewind : bool
            should the file be returned to the starting point
        macro_rewind : bool
            ???

        Returns
        -------
        markers : int
            a single marker

        """
        op2 = self.op2
        ni = op2.n
        #markers = []
        data = self.read_block()
        marker, = op2.struct_i.unpack(data)
        if rewind:
            op2.n = ni
            op2.f.seek(op2.n)
        else:
            if self.is_debug_file:
                msg = 'get_marker - [4, %i, 4]; macro_rewind=%s\n' % (
                    marker, macro_rewind or rewind)
                self.binary_debug.write(msg)
        return marker

    #------------------------------------------------------------------
    @property
    def is_debug_file(self):
        """interface to the op2 object"""
        return self.op2.is_debug_file

    @property
    def binary_debug(self):
        """interface to the op2 object"""
        return self.op2.binary_debug

    @property
    def log(self):
        """interface to the op2 object"""
        return self.op2.log

    @property
    def _endian(self):
        """interface to the op2 object"""
        return self.op2._endian

    #------------------------------------------------------------------
    # skip methods

#class OP2Skip(object):
    #def __init__(self, op2):
        #self.op2 = op2

    #@property
    #def is_debug_file(self):
        #return self.op2.is_debug_file

    #@property
    #def binary_debug(self):
        #return self.op2.binary_debug

    #@property
    #def log(self):
        #return self.op2.log

    def _skip_table_helper(self):
        """
        Skips the majority of geometry/result tables as they follow a very standard format.
        Other tables don't follow this format.

        """
        op2 = self.op2
        op2.table_name = self._read_table_name(rewind=False)
        if self.is_debug_file:
            self.binary_debug.write('skipping table...%r\n' % op2.table_name)
        self.log.warning('    skipping table_helper = %s' % op2.table_name)

        self.read_markers([-1])
        data = self._skip_record()
        self.read_markers([-2, 1, 0])
        data = self._skip_record()
        self._skip_subtables()

    def _skip_subtables(self):
        """skips a set of subtables"""
        op2 = self.op2
        op2.isubtable = -3
        self.read_markers([-3, 1, 0])

        markers = self.get_nmarkers(1, rewind=True)
        while markers[0] != 0:
            unused_data = self._skip_record()
            if self.is_debug_file:
                self.log.debug("skipping table_name = %r" % op2.table_name)
            #if len(data) == 584:
                #self._parse_results_table3(data)
            #else:
                #data = self._parse_results_table4(data)

            op2.isubtable -= 1
            self.read_markers([op2.isubtable, 1, 0])
            markers = self.get_nmarkers(1, rewind=True)
        self.read_markers([0])

    def _skip_record(self):
        """
        the skip version of ``_read_record``

        Returns
        -------
        record : None
            a record of None indicates a skipped block

        """
        op2 = self.op2
        unused_markers0 = self.get_nmarkers(1, rewind=False)
        record = self._skip_block()

        markers1 = self.get_nmarkers(1, rewind=True)
        # handling continuation blocks
        while markers1[0] > 0:
            markers1 = op2.get_nmarkers(1, rewind=False)
            record = self._skip_block()
            markers1 = self.get_nmarkers(1, rewind=True)
        return record

    def _skip_record_ndata(self, stream=False, debug=True, macro_rewind=False):
        """the skip version of ``_read_record_ndata``"""
        op2 = self.op2
        markers0 = self.get_nmarkers(1, rewind=False, macro_rewind=macro_rewind)
        if self.is_debug_file and debug:
            self.binary_debug.write('read_record - marker = [4, %i, 4]; macro_rewind=%s\n' % (
                markers0[0], macro_rewind))
        record, nrecord = self._skip_block_ndata()

        if self.is_debug_file and debug:
            self.binary_debug.write('read_record - record = [%i, recordi, %i]; '
                                    'macro_rewind=%s\n' % (nrecord, nrecord, macro_rewind))
        if markers0[0]*4 != nrecord:
            msg = 'markers0=%s*4 len(record)=%s; table_name=%r' % (
                markers0[0]*4, nrecord, op2.table_name)
            raise FortranMarkerError(msg)

        markers1 = self.get_nmarkers(1, rewind=True)

        if markers1[0] > 0:
            while markers1[0] > 0:
                markers1 = self.get_nmarkers(1, rewind=False)
                if self.is_debug_file and debug:
                    self.binary_debug.write('read_record - markers1 = [4, %i, 4]\n' % markers1[0])
                unused_recordi, nrecordi = self._skip_block_ndata()
                nrecord += nrecordi

                markers1 = self.get_nmarkers(1, rewind=True)
                if self.is_debug_file and debug:
                    self.binary_debug.write('read_record - markers1 = [4, %i, 4]\n' % markers1[0])
        return record, nrecord

    def _skip_block(self):
        """
        Skips a block following a pattern of:
            [nbytes, data, nbytes]

        Returns
        -------
        data :  since data can never be None, a None value
                indicates something bad happened.
        """
        return self._skip_block_ndata()[0]

    def _skip_block_ndata(self):
        """
        Skips a block following a pattern of:
            [nbytes, data, nbytes]

        Returns
        -------
        data :  since data can never be None, a None value
                indicates something bad happened.

        """
        op2 = self.op2
        data = op2.f.read(4)
        ndata, = op2.struct_i.unpack(data)
        op2.n += 8 + ndata
        self._goto(op2.n)
        return None, ndata

    def _goto(self, n):
        """
        Jumps to position n in the file

        Parameters
        ----------
        n : int
            the position to goto

        """
        self.op2.n = n
        self.op2.f.seek(n)


    def is_valid_subcase(self):
        """
        Lets the code check whether or not to read a subcase

        Returns
        -------
        is_valid : bool
            should this subcase defined by self.isubcase be read?

        """
        op2 = self.op2
        if not op2.is_all_subcases:
            if hasattr(op2, 'isubcase') and op2.isubcase in op2.valid_subcases:
                return True
            return False
        return True
