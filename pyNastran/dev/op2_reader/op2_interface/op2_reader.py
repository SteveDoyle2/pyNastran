"""
Defines various tables that don't fit in other sections:
  - OP2Reader
    - read_cmodeext(self)
    - read_cmodeext_helper(self)
    - read_aemonpt(self)
    - read_monitor(self)
    - read_r1tabrg(self, data, ndata)
    - read_hisadd(self)

    - read_cstm(self)
    - read_dit(self)
    - read_extdb(self)
    - read_fol(self)
    - read_frl(self)
    - read_gpl(self)
    - read_ibulk(self)
    - read_intmod(self)
    - read_meff(self)
    - read_omm2(self)
    - read_sdf(self)
    - read_tol(self)
    - _skip_pcompts(self)
    - _read_pcompts(self)

  - Matrix
    - _get_matrix_row_fmt_nterms_nfloats(self, nvalues, tout)
    - _skip_matrix_mat(self)
    - read_matrix(self, table_name)
    - _read_matpool_matrix(self)
    - _read_matrix_mat(self)
    - grids_comp_array_to_index(grids1, comps1, grids2, comps2,
                                make_matrix_symmetric)

  - Others
    - _get_marker_n(self, nmarkers)
    - read_markers(self)
    - _skip_subtables(self)
    - _skip_table_helper(self)
    - _print_month(self, month, day, year, zero, one)
    - read_results_table(self)
"""
from __future__ import print_function, unicode_literals
import sys
from copy import deepcopy
from itertools import count
from struct import unpack, Struct
from six import b, string_types
import numpy as np
from numpy import frombuffer
import scipy  # type: ignore

from pyNastran.utils.numpy_utils import integer_types
from pyNastran.f06.errors import FatalError
from pyNastran.op2.errors import FortranMarkerError, SortCodeError
from pyNastran.op2.tables.design_response import (
    WeightResponse, StressResponse, StrainResponse, ForceResponse,
    FlutterResponse, Convergence)
from pyNastran.op2.tables.matrix import Matrix

from pyNastran.op2.op2_interface.utils import apply_mag_phase, update_label2, get_superelement_adaptivity_index
from pyNastran.op2.op2_interface.op2_codes import determine_sort_bits_meaning

from pyNastran.op2.tables.lama_eigenvalues.lama import LAMA
from pyNastran.op2.tables.oes_stressStrain.oes import OES
from pyNastran.op2.tables.oqg_constraintForces.oqg import OQG
from pyNastran.op2.tables.oee_energy.onr import ONR
from pyNastran.op2.tables.ogf_gridPointForces.ogpf import OGPF

from pyNastran.op2.tables.oef_forces.oef import OEF
#from pyNastran.op2.tables.oes_stressStrain.oesm import OESM
from pyNastran.op2.tables.ogs_grid_point_stresses.ogs import OGS

from pyNastran.op2.tables.opg_appliedLoads.opg import OPG
from pyNastran.op2.tables.oug.oug import OUG
from pyNastran.op2.tables.ogpwg import OGPWG

#class MinorTables(object):
    #def __init__(self, op2_reader):
        #self.op2_reader = op2_reader

from pyNastran.op2.op2_interface.op2_codes import Op2Codes

class OP2Reader(Op2Codes):
    """Stores methods that aren't useful to an end user"""
    def __init__(self, op2):
        Op2Codes.__init__(self)

        # how many optimization passes have there been
        self._count = 0

        self.op2 = op2
        self.table_reader = TableReader(self)
        self.lama = LAMA(self, op2)
        self.oes = OES(self, op2)
        self.oqg = OQG(self, op2)
        self.oug = OUG(self, op2)
        #self.minor_tables = MinorTables(self)

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
            b'FOL' : self.read_fol,  # frequency response list
            b'FRL' : self.read_frl,  # frequency response list
            b'SDF' : self.read_sdf,
            b'IBULK' : self.read_ibulk,
            b'CDDATA' : self.read_ibulk,
            b'CMODEXT' : self._read_cmodext,
            b'CSTM' : self.read_cstm,  # coordinate system transformation matrices
        }
        #self.op2_skip = OP2Skip(op2)
        self.numwide = -1

    @property
    def idtype(self):
        return self.op2.idtype
    @property
    def fdtype(self):
        return self.op2.fdtype

    @property
    def use_vector(self):
        return self.op2.use_vector

    @property
    def result_names(self):
        return self.op2.result_names

    @property
    def encoding(self):
        """interface to the op2 object"""
        return self.op2.encoding

    @property
    def f(self):
        return self.op2.f
    #@f.setter
    #def f(self, f):
        #self.op2.f = f

    @property
    def n(self):
        return self.op2.n
    @n.setter
    def n(self, n):
        self.op2.n = n

    def _function1(self, value):
        """function1(value)"""
        if value // 1000 in [2, 3, 6]:
            return 2
        return 1

    def set_null_nonlinear_factor(self):
        """
        Initializes the nonlinear factor, which lets us know if
        this is a transient solution or not.

        """
        self.nonlinear_factor = None
        self.data_code['nonlinear_factor'] = None

    def _read_title(self, data):
        self._read_title_helper(data)

        if hasattr(self, 'isubcase'):
            op2 = self.op2
            if self.isubcase not in op2.isubcase_name_map:
                # 100 from label
                # 20 from subtitle line
                # 'SUBCASE 2'
                #op2.isubcase_name_map[isubcase] = [self.subtitle, self.label]
                op2.isubcase_name_map[self.isubcase] = [
                    self.subtitle, self.superelement_adaptivity_index,
                    self.analysis_code, self.label]
        else:
            raise  RuntimeError('isubcase is not defined')

        if hasattr(self, 'subtitle') and hasattr(self, 'label'):
            ogs = 0
            if hasattr(self, 'ogs'):
                ogs = self.ogs
            unused_code = (self.isubcase, self.analysis_code,
                           self.superelement_adaptivity_index,
                           self.pval_step, ogs)
            #code = (self.isubcase, self.analysis_code,
                    #self.superelement_adaptivity_index, self.table_name_str)
            #print("code =", code)
            #if code not in op2.labels:
                #op2.subtitles[op2.isubcase].append(self.subtitle)
                #op2.labels[code] = self.label

    def _read_title_helper(self, data):
        assert len(data) == 584, len(data)
        # titleSubtitleLabel
        title, subtitle, label = unpack(self._endian + b'128s128s128s', data[200:])
        self.title = title.decode(self.encoding).strip()
        subtitle = subtitle.decode(self.encoding)

        label = label.decode(self.encoding).strip()
        nlabel = 65
        label2 = label[nlabel:]
        try:
            label2 = update_label2(label2, self.isubcase)
        except AssertionError:
            pass

        assert len(label[:nlabel]) <= nlabel, 'len=%s \nlabel     =%r \nlabel[:%s]=%r' % (len(label), label, nlabel, label[:nlabel])
        assert len(label2) <= 55, 'len=%s label = %r\nlabel[:%s]=%r\nlabel2    =%r' % (len(label2), label, nlabel, label[:nlabel], label2)
        # not done...
        # 65 + 55 = 120 < 128

        self.label = label
        self.pval_step = label2

        #split_label = label.split()
        #if len(split_label) == 2:
            #word, value1 = split_label
            #assert word == 'SUPERELEMENT', 'split_label=%s' % split_label
            #subtitle = '%s; SUPERELEMENT %s' % (subtitle, value1)
            #value1 = int(value1)

            #if superelement_adaptivity_index:
                #superelement_adaptivity_index = '%s; SUPERELEMENT %s' % (
                    #superelement_adaptivity_index, value1)
            #else:
                #superelement_adaptivity_index = 'SUPERELEMENT %ss' % value1
        #elif len(split_label) == 4:
            #word, value1, comma, value2 = split_label
            #assert word == 'SUPERELEMENT', 'split_label=%s' % split_label
            #value1 = int(value1)
            #value2 = int(value2)

            #if superelement_adaptivity_index:
                #superelement_adaptivity_index = '%s; SUPERELEMENT %s,%s' % (
                    #superelement_adaptivity_index, value1, value2)
            #else:
                #superelement_adaptivity_index = 'SUPERELEMENT %s,%s' % (value1, value2)
        #else:
            #raise RuntimeError(split_label)


        nsubtitle_break = 67
        adpativity_index = subtitle[nsubtitle_break:99].strip()
        superelement = subtitle[99:].strip()

        subtitle = subtitle[:nsubtitle_break].strip()
        assert len(superelement) <= 26, 'len=%s superelement=%r' % (len(superelement), superelement)
        superelement = superelement.strip()

        assert len(subtitle) <= 67, 'len=%s subtitle=%r' % (len(subtitle), subtitle)
        superelement_adaptivity_index = get_superelement_adaptivity_index(subtitle, superelement)

        if adpativity_index:
            assert 'ADAPTIVITY INDEX=' in adpativity_index
            # F:\work\pyNastran\examples\Dropbox\move_tpl\pet1018.op2
            #'ADAPTIVITY INDEX=      1'
            split_adpativity_index = adpativity_index.split()
            assert len(split_adpativity_index) == 3, split_adpativity_index
            word1, word2, adpativity_index_value = split_adpativity_index
            assert word1 == 'ADAPTIVITY', 'split_adpativity_index=%s' % split_adpativity_index
            assert word2 == 'INDEX=', 'split_adpativity_index=%s' % split_adpativity_index

            adpativity_index_value = int(adpativity_index_value)
            subtitle = '%s; ADAPTIVITY_INDEX=%s' % (subtitle, adpativity_index_value)
            if superelement_adaptivity_index:
                superelement_adaptivity_index = '%s; ADAPTIVITY_INDEX=%s' % (
                    superelement_adaptivity_index, adpativity_index_value)
            else:
                superelement_adaptivity_index = 'ADAPTIVITY_INDEX=%s' % adpativity_index_value

        self.subtitle = subtitle
        self.superelement_adaptivity_index = superelement_adaptivity_index
        assert len(self.label) <= 124, 'len=%s label=%r' % (len(self.label), self.label)

        #: the subtitle of the subcase
        data_code = self.data_code
        data_code['subtitle'] = self.subtitle

        # the sub-key
        data_code['pval_step'] = self.pval_step
        data_code['superelement_adaptivity_index'] = self.superelement_adaptivity_index

        #: the label of the subcase
        data_code['label'] = self.label
        data_code['title'] = self.title

        if self.is_debug_file:
            self.binary_debug.write(
                '  %-14s = %r\n' * 6 % (
                    'count', self._count,
                    'title', self.title,
                    'subtitle', self.subtitle,
                    'label', self.label,
                    'pval_step', self.pval_step,
                    'superelement_adaptivity_index', self.superelement_adaptivity_index))

    def _write_debug_bits(self):
        """
        s_code =  0 -> stress_bits = [0,0,0,0,0]
        s_code =  1 -> stress_bits = [0,0,0,0,1]
        s_code =  2 -> stress_bits = [0,0,0,1,0]
        s_code =  3 -> stress_bits = [0,0,0,1,1]
        etc.
        s_code = 32 -> stress_bits = [1,1,1,1,1]

        stress_bits[0] = 0 -> isMaxShear=True       isVonMises=False
        stress_bits[0] = 1 -> isMaxShear=False      isVonMises=True

        stress_bits[1] = 0 -> is_stress=True        is_strain=False
        stress_bits[2] = 0 -> isFiberCurvature=True isFiberDistance=False
        stress_bits[3] = 0 -> duplicate of Bit[1] (stress/strain)
        stress_bits[4] = 0 -> material coordinate system flag

        """
        if self.is_debug_file:
            msg = ''
            assert len(self.words) in [0, 28], 'table_name=%r len(self.words)=%s words=%s' % (self.table_name, len(self.words), self.words)
            for i, param in enumerate(self.words):
                if param == 's_code':
                    try:
                        s_word = get_scode_word(self.s_code, self.stress_bits)
                    except AttributeError:
                        raise
                    self.binary_debug.write('  s_code         = %s -> %s\n' % (self.s_code, s_word))
                    self.binary_debug.write('    stress_bits[0] = %i -> is_von_mises    =%-5s vs is_max_shear\n' % (self.stress_bits[0], self.is_von_mises))
                    self.binary_debug.write('    stress_bits[1] = %i -> is_strain       =%-5s vs is_stress\n' % (self.stress_bits[1], self.is_strain))
                    self.binary_debug.write('    stress_bits[2] = %i -> strain_curvature=%-5s vs fiber_dist\n' % (self.stress_bits[2], self.is_curvature))
                    self.binary_debug.write('    stress_bits[3] = %i -> is_strain       =%-5s vs is_stress\n' % (self.stress_bits[3], self.is_strain))
                    self.binary_debug.write('    stress_bits[4] = %i -> material coordinate system flag=%s vs ???\n' % (self.stress_bits[4], self.stress_bits[4]))
                elif param == '???':
                    param = 0
                msg += '%s, ' % param
                if i % 5 == 4:
                    msg += '\n             '

            if hasattr(self, 'format_code'):
                try:
                    is_complex = self.is_complex
                except AssertionError:
                    self.binary_debug.write('\n  ERROR: cannot determine is_complex properly; '
                                            'check_sort_bits!!!\n')
                    is_complex = '???'

                try:
                    is_random = self.is_random
                except AssertionError:
                    is_random = '???'

                try:
                    is_sort1 = self.is_sort1
                except AssertionError:
                    is_sort1 = '???'

                try:
                    is_real = self.is_real
                except AssertionError:
                    is_real = '???'

                if is_complex:
                    msg = '\n  %-14s = %i -> is_mag_phase vs is_real_imag vs. is_random\n' % (
                        'format_code', self.format_code)
                    self.binary_debug.write(msg)
                else:
                    self.binary_debug.write('  %-14s = %i\n' % ('format_code', self.format_code))
                self.binary_debug.write('    sort_bits[0] = %i -> is_random=%s vs mag/phase\n' % (self.sort_bits[0], is_random))
                self.binary_debug.write('    sort_bits[1] = %i -> is_sort1 =%s vs sort2\n' % (self.sort_bits[1], is_sort1))
                self.binary_debug.write('    sort_bits[2] = %i -> is_real  =%s vs real/imag\n' % (self.sort_bits[2], is_real))

                try:
                    sort_method, is_real, is_random = self._table_specs()
                    self.binary_debug.write('    sort_method = %s\n' % sort_method)
                except AssertionError:
                    self.binary_debug.write('    sort_method = ???\n')

                if is_complex:
                    msg = '\n  %-14s = %i -> is_mag_phase vs is_real_imag vs. is_random\n' % (
                        'format_code', self.format_code)
                    self.binary_debug.write(msg)
                else:
                    self.binary_debug.write('  %-14s = %i\n' % ('format_code', self.format_code))

            self.binary_debug.write('  recordi = [%s]\n\n' % msg)

    def _setup_op2_subcase(self, word):
        """
        Parameters
        ----------
        word : str
            displacement
            FLUX
        """
        if self.read_mode == 1:
            op2 = self.op2
            case_control_deck = op2.case_control_deck
            if self.isubcase not in case_control_deck.subcases:
                op2.subcase = case_control_deck.create_new_subcase(self.isubcase)
            else:
                op2.subcase = case_control_deck.subcases[self.isubcase]
            op2.subcase.add_op2_data(self.data_code, word, self.log)

    def read_nastran_version(self):
        """reads the version header"""
        #try:
        op2 = self.op2
        markers = self.get_nmarkers(1, rewind=True)
        #except:
            #self._goto(0)
            #try:
                #self.f.read(4)
            #except:
                #raise FatalError("The OP2 is empty.")
            #raise
        if self.is_debug_file:
            if self.read_mode == 1:
                self.binary_debug.write('read_mode = %s (vectorized; 1st pass)\n' % self.read_mode)
            elif self.read_mode == 2:
                self.binary_debug.write('read_mode = %s (vectorized; 2nd pass)\n' % self.read_mode)

        if markers == [3,]:  # PARAM, POST, -1
            if self.is_debug_file:
                self.binary_debug.write('marker = 3 -> PARAM,POST,-1?\n')
            op2.post = -1
            self.read_markers([3])
            data = self.read_block()   # TODO: is this the date?
            #assert len(data) == 12, len(data)

            self.read_markers([7])
            data = self.read_block()

            if data == b'NASTRAN FORT TAPE ID CODE - ':
                macro_version = 'nastran'
            elif b'IMAT v' in data:
                imat_version = data[6:11].encode('utf8')
                macro_version = 'IMAT %s' % imat_version
            else:
                version_ints = Struct(self._endian + b'7i').unpack(data)
                if version_ints == (1, 2, 3, 4, 5, 6, 7):
                    macro_version = 'MSFC'
                else:
                    self.show_data(data)
                    raise NotImplementedError(data)
                #self.show_data(data)
            if self.is_debug_file:
                self.binary_debug.write('%r\n' % data)
            #print('macro_version = %r' % macro_version)

            data = self._read_record()
            if self.is_debug_file:
                self.binary_debug.write('%r\n' % data)
            version = data.strip()

            if macro_version == 'nastran':
                if version.startswith(b'NX'):
                    op2.set_as_nx()
                    op2.set_table_type()
                elif version.startswith(b'MODEP'):
                    # TODO: why is this separate?
                    # F:\work\pyNastran\pyNastran\master2\pyNastran\bdf\test\nx_spike\out_ac11103.op2
                    op2.set_as_nx()
                    op2.set_table_type()
                elif version.startswith(b'AEROFREQ'):
                    # TODO: why is this separate?
                    # C:\Users\Steve\Dropbox\pyNastran_examples\move_tpl\loadf.op2
                    op2.set_as_msc()
                    op2.set_table_type()
                elif version.startswith(b'AEROTRAN'):
                    # TODO: why is this separate?
                    # C:\Users\Steve\Dropbox\pyNastran_examples\move_tpl\loadf.op2
                    op2.set_as_msc()
                    op2.set_table_type()
                elif version in [b'XXXXXXXX', b'V2005R3B']:
                    op2.set_as_msc()
                    op2.set_table_type()
                elif version == b'OS12.210':
                    op2.set_as_optistruct()
                    op2.set_table_type()
                elif version == b'OS11XXXX':
                    op2.set_as_radioss()
                    op2.set_table_type()
                #elif data[:20] == b'XXXXXXXX20141   0   ':
                    #self.set_as_msc()
                    #self.set_table_type()
                else:
                    raise RuntimeError('unknown version=%r' % version)

            elif macro_version.startswith('IMAT'):
                assert version.startswith(b'ATA'), version
                op2._nastran_format = macro_version
            elif macro_version == 'MSFC':
                op2._nastran_format = macro_version
                #self.show_data(version)
                #assert version.startswith(b'ATA'), version

            if self.is_debug_file:
                self.binary_debug.write(data.decode(self._encoding) + '\n')
            self.read_markers([-1, 0])
        elif markers == [2,]:  # PARAM, POST, -2
            if self.is_debug_file:
                self.binary_debug.write('marker = 2 -> PARAM,POST,-2?\n')
            op2.post = -2
        else:
            raise NotImplementedError(markers)

    def _read_aemonpt(self):
        """reads the AEMONPT table"""
        #self.log.debug("table_name = %r" % op2.table_name)
        unused_table_name = self._read_table_name(rewind=False)

        #print('-----------------------')
        #print('record 1')
        self.read_markers([-1])
        data = self._read_record()
        #self.show_data(data)
        if self.read_mode == 2:
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
        #if self.read_mode == 2:
        word, = unpack(self._endian + b'8s', data)
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
        #self.show_data(data)
        self.read_markers([0])

        #print('-----------------------')
        #print('end')
        #self.show(200)
        #aaa

    def _read_monitor(self):
        """reads the MONITOR table"""
        op2 = self.op2
        self.log.debug("table_name = %r" % op2.table_name)
        unused_table_name = self._read_table_name(rewind=False)

        #print('-----------------------')
        #print('record 1')
        self.read_markers([-1])
        data = self._read_record()
        #self.show_data(data)
        if self.read_mode == 2:
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
        if self.read_mode == 2:
            word, = op2.struct_8s.unpack(data)
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
            (unused_aero, name, comps, cp, x, y, z, unused_coeff, word, column, cd,
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
        #self.show_data(data)
        self.read_markers([0])

        #print('-----------------------')
        #print('end')
        #self.show(200)

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
        unused_table_name, oneseventy_a, oneseventy_b = unpack('8sii', data)
        assert oneseventy_a == 170, oneseventy_a
        assert oneseventy_b == 170, oneseventy_b
        print('170*4 =', 170*4)
        #self.show_data(data)
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
                self.show_data(data, types='ifs')
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
        unused_table_name = self._read_table_name(rewind=False)
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
        for i, unused_coord in enumerate(ints):
            cid = ints[i, 0]
            coord_type_int = ints[i, 1]
        if coord_type_int in coord_type_map:
            unused_coord_type = coord_type_map[coord_type_int]
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
        unused_table_name = self._read_table_name(rewind=False)
        self.read_markers([-1])
        data = self._read_record()

        self.read_markers([-2, 1, 0])
        data = self._read_record()
        unused_table_name, = op2.struct_8s.unpack(data)

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
        unused_data = self._read_record()

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
                unused_data = self._read_record()
                #self.show_data(data, types='s', endian=None)
            elif nfields == 0:
                #self.show_ndata(100, types='ifs')
                break
            else:
                raise RuntimeError('nfields=%s' % nfields)
            marker -= 1
        unused_marker_end = op2.get_marker1(rewind=False)

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
        fmt = b(self._uendian + '%sf' % nfloats)
        freqs = np.array(list(unpack(fmt, data[8:])), dtype='float32')

        if self.read_mode == 2:
            if op2._frequencies is not None and not np.array_equal(freqs, op2._frequencies):
                msg = (
                    'Cannot overwrite op2._frequencies...\n'
                    'op2._frequencies = %s\n'
                    'new_freqs = %s\n' % (op2._frequencies, freqs))
                raise RuntimeError(msg)
            op2._frequencies = freqs
            if self.is_debug_file:
                self.binary_debug.write('  recordi = [%r, freqs]\n'  % (subtable_name_raw))
                self.binary_debug.write('  subtable_name=%r\n' % subtable_name)
                self.binary_debug.write('  freqs = %s' % freqs)
        self._read_subtables()

    def read_frl(self):
        """reads the FRL (Frequency Response List) table"""
        op2 = self.op2
        op2.table_name = self._read_table_name(rewind=False)
        self.read_markers([-1])
        data = self._read_record()
        idata = unpack(self._endian + b'7i', data)
        assert idata[0] == 101, idata
        assert idata[1] == 1, idata
        assert idata[2] == 0, idata
        assert idata[3] == 0, idata
        assert idata[4] == 0, idata
        assert idata[5] == 0, idata
        assert idata[6] == 0, idata
        #print(self.show_data(data))


        self.read_markers([-2, 1, 0])
        data = self._read_record()
        assert len(data) == 12, '\n'.join(str(d) for d in self.show_data(data))

        subtable_name_raw, = op2.struct_8s.unpack(data[:8])
        subtable_name = subtable_name_raw.strip()
        assert subtable_name == b'FRL0', 'subtable_name=%r' % subtable_name

        self.read_markers([-3, 1, 0])
        isubtable = -3
        markers = self.get_nmarkers(1, rewind=True)
        while markers[0] != 0:
            if self.read_mode == 1:
                self._skip_record()
            else:
                data = self._read_record()
                #self.show_data(data)
                freqs = np.frombuffer(data, dtype=op2.fdtype).copy()
                #print('read_mode=%s freqs=%s' % (self.read_mode, freqs.tolist()))
                if op2._frequencies is not None and not np.array_equal(freqs, op2._frequencies):
                    msg = (
                        'Cannot overwrite op2._frequencies...\n'
                        'op2._frequencies = %s\n'
                        'new_freqs = %s\n' % (op2._frequencies, freqs))
                    raise RuntimeError(msg)
                op2._frequencies = freqs

            isubtable -= 1
            self.read_markers([isubtable, 1, 0])
            markers = self.get_nmarkers(1, rewind=True)
        del isubtable
        self.read_markers([0])

    def read_gpl(self):
        """reads the GPL table (grid point list?)"""
        op2 = self.op2
        op2.table_name = self._read_table_name(rewind=False)
        self.log.debug('table_name = %r' % op2.table_name)
        if self.is_debug_file:
            self.binary_debug.write('read_geom_table - %s\n' % op2.table_name)
        self.read_markers([-1])
        if self.is_debug_file:
            self.binary_debug.write('---markers = [-1]---\n')
        unused_data = self._read_record()

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
            unused_data = self._read_record()
            #self.show_data(data, 'i')
            n -= 1
            markers = self.get_nmarkers(1, rewind=True)

    def read_hisadd(self):
        """optimization history (SOL200) table"""
        op2 = self.op2
        op2.table_name = self._read_table_name(rewind=False)

        if self.read_mode == 1:
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
        #self.show_data(data)

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
            self.show_data(data, types='ifs', endian=None)
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
            self.binary_debug.write('read_geom_table - %s\n' % op2.table_name)
        self.read_markers([-1])
        if self.is_debug_file:
            self.binary_debug.write('---markers = [-1]---\n')
        unused_data = self._read_record()

        unused_markers = self.get_nmarkers(1, rewind=True)
        marker = -2
        while 1:
            self.read_markers([marker, 1, 0])
            nfields = self.get_marker1(rewind=True)
            if nfields > 0:
                unused_data = self._read_record()
                #self.show_data(data, types='s', endian=None)
            elif nfields == 0:
                #op2.show_ndata(100, types='ifs')
                break
            else:
                raise RuntimeError('nfields=%s' % nfields)
            marker -= 1
        unused_marker_end = self.get_marker1(rewind=False)

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
            raise NotImplementedError(self.show_data(data))
        self._read_subtables()

    def _read_pcompts(self):
        """
        Reads the PCOMPTS table (poorly).
        The PCOMPTS table stores information about the PCOMP cards???
        """
        self._skip_pcompts()
        return
        #if self.read_mode == 1:
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
        unused_table_name = self._read_table_name(rewind=False)

        self.read_markers([-1])
        unused_data = self._skip_record()

        self.read_markers([-2, 1, 0])
        unused_data = self._skip_record()
        #table_name, = op2.struct_8s.unpack(data)

        self.read_markers([-3, 1, 0])
        markers = self.get_nmarkers(1, rewind=True)
        if markers != [-4]:
            unused_data = self._skip_record()

        self.read_markers([-4, 1, 0])
        markers = self.get_nmarkers(1, rewind=True)
        if markers != [0]:
            unused_data = self._skip_record()
        else:
            self.read_markers([0])
            return

        self.read_markers([-5, 1, 0])
        unused_data = self._skip_record()

        self.read_markers([-6, 1, 0])
        self.read_markers([0])

    def read_meff(self):
        """reads the MEFF table"""
        op2 = self.op2
        op2.table_name = self._read_table_name(rewind=False)
        op2.log.debug('table_name = %r' % op2.table_name)
        if self.is_debug_file:
            self.binary_debug.write('read_geom_table - %s\n' % op2.table_name)
        self.read_markers([-1])
        if self.is_debug_file:
            self.binary_debug.write('---markers = [-1]---\n')
        unused_data = self._read_record()

        markers = self.get_nmarkers(1, rewind=True)
        if self.is_debug_file:
            self.binary_debug.write('---marker0 = %s---\n' % markers)
        self.read_markers([-2, 1, 0])
        unused_data = self._read_record()

        for n in [-3, -4, -5, -6, -7, -8]:
            self.read_markers([n, 1, 1])
            markers = self.get_nmarkers(1, rewind=False)
            #print('markers =', markers)
            nbytes = markers[0]*4 + 12
            unused_data = self.f.read(nbytes)
            op2.n += nbytes
        n = -9
        self.read_markers([n, 1, 0, 0])

    def read_intmod(self):
        """reads the INTMOD table"""
        op2 = self.op2
        op2.table_name = self._read_table_name(rewind=False)
        #op2.log.debug('table_name = %r' % op2.table_name)
        if self.is_debug_file:
            self.binary_debug.write('read_geom_table - %s\n' % op2.table_name)
        self.read_markers([-1])
        if self.is_debug_file:
            self.binary_debug.write('---markers = [-1]---\n')
        unused_data = self._read_record()
        #print('intmod data1')
        #self.show_data(data)

        markers = self.get_nmarkers(1, rewind=True)
        if self.is_debug_file:
            self.binary_debug.write('---marker0 = %s---\n' % markers)
        self.read_markers([-2, 1, 0])
        unused_data = self._read_record()
        #print('intmod data2')
        #self.show_data(data)

        for n in [-3, -4, -5, -6, -7, -8,]:
            self.read_markers([n, 1, 1])
            markers = self.get_nmarkers(1, rewind=False)
            #print('markers =', markers)
            nbytes = markers[0]*4 + 12
            unused_data = self.f.read(nbytes)
            #print('intmod data%i' % n)
            #self.show_data(data)
            op2.n += nbytes

        n = -9
        self.read_markers([n, 1, 0, 0])
        #self.show(50)
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
        if self._table4_count == 0:
            self._count += 1
        self._table4_count += 1

        #if self._table4_count == 0:
            #self._count += 1
        #self._table4_count += 1

        if self.read_mode == 1:
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
                #print('internal_id=%s dresp_id=%s response_type=%s response_label=%s'
                      #' region=%s subcase=%s type_flag=%s seid=%s' % (
                    #internal_id, dresp_id, response_type, response_label,
                    #region, subcase, type_flag, seid
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
                #msg = ('FLUTTER - _count=%s label=%r region=%s subcase=%s mode=%s '
                       #'mach=%s velocity=%s density=%s flutter_id=%s' % (
                           #self._count, response_label, region, subcase, mode,
                           #mach, velocity, density, flutter_id))
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
            strings, ints, floats = self.show_data(data)
            msg = 'len(data) = %i\n' % ndata
            msg += 'strings  = %r\n' % strings
            msg += 'ints     = %r\n' % str(ints)
            msg += 'floats   = %r' % str(floats)
            raise NotImplementedError(msg)

        self.read_markers([-3, 1, 1])

        unused_markers0 = self.get_nmarkers(1, rewind=False)
        unused_record = self.read_block()

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
        unused_table_name = self._read_table_name(rewind=False, stop_on_failure=True)
        self.read_markers([-1])
        unused_data = self._read_record()
        #self.show_data(data)

        self.read_markers([-2, 1, 0])
        #op2.show_ndata(440, types='if')
        unused_data = self._read_record()
        #print('----')
        self.read_markers([-3, 1, 0])
        #op2.show_ndata(440, types='if')
        #print('----')
        self.read_markers([0])
        #data = self._read_record()


        #op2.show_ndata(440, types='ifs')

        #self.show_data(data)
        #aaaa

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
            fmt = b(self._uendian + 'i %if' % nfloats)
        elif tout == 2:
            nfloats = nvalues // 2
            nterms = nvalues // 2
            fmt = b(self._uendian + 'i %id' % nfloats)
        elif tout == 3:
            nfloats = nvalues
            nterms = nvalues // 2
            fmt = b(self._uendian + 'i %if' % nfloats)
        elif tout == 4:
            nfloats = nvalues // 2
            nterms = nvalues // 4
            fmt = b(self._uendian + 'i %id' % nfloats)
        else:
            raise RuntimeError('tout = %s' % tout)
        return fmt, nfloats, nterms

    def _read_matrix_mat(self):
        """
        Reads a matrix in "standard" form.  The forms are::
            standard:
                Return a matrix that looks similar to a matrix found
                in the OP4.  Created by:
                ``ASSIGN output2='model.op2', UNIT=12,UNFORMATTED,DELETE``
                ``OUTPUT2 KGG//0/12``
            matpool:
                Return a matrix that looks similar to a DMIG matrix
                (e.g., it contains the node id and DOF).  Created by:
                ``ASSIGN output2='model.op2', UNIT=12,UNFORMATTED,DELETE``
                ``TODO: add the magic keyword...``
                ``OUTPUT2 KGG//0/12``

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
        op2 = self.op2
        allowed_forms = [1, 2, 3, 4, 5, 6, 8, 9, 10, 11, 13, 15]
        #self.log.debug('----------------------------------------------------------------')
        table_name = self._read_table_name(rewind=False, stop_on_failure=True)
        self.read_markers([-1])
        data = self._read_record()

        # old-bad
        #matrix_num, form, mrows, ncols, tout, nvalues, g = unpack(self._endian + b'7i', data)

        #           good   good   good  good  ???    ???
        matrix_num, ncols, mrows, form, tout, nvalues, g = unpack(self._endian + b'7i', data)
        #print('g =', g)

        m = Matrix(table_name, form=form)
        op2.matrices[table_name.decode('utf-8')] = m

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

        #self.log.error('name=%r matrix_num=%s form=%s mrows=%s '
        #               'ncols=%s tout=%s nvalues=%s g=%s' % (
        #                   table_name, matrix_num, form, mrows, ncols, tout, nvalues, g))
        if form == 1:
            if ncols != mrows:
                self.log.warning('unexpected size for %s; form=%s mrows=%s ncols=%s' % (
                    table_name, form, mrows, ncols))
        elif form not in allowed_forms:
            self.log.error('name=%r matrix_num=%s form=%s mrows=%s '
                           'ncols=%s tout=%s nvalues=%s g=%s' % (
                               table_name, matrix_num, form, mrows, ncols,
                               tout, nvalues, g))
            raise RuntimeError('form=%s; allowed=%s' % (form, allowed_forms))
        self.log.debug('name=%r matrix_num=%s form=%s mrows=%s ncols=%s tout=%s '
                       'nvalues=%s g=%s' % (
                           table_name, matrix_num, form, mrows, ncols, tout, nvalues, g))

        self.read_markers([-2, 1, 0])
        data = self._read_record()

        if len(data) == 16:
            unused_name, ai, bi = unpack(self._endian + b'8s 2i', data)
            assert ai == 170, ai
            assert bi == 170, bi
        else:
            self.log.warning('unexpected matrix length=%s' % len(data))
            self.log.warning(self.show_data(data, types='if'))

        itable = -3
        unused_j = None

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
                    fmt, unused_nfloats, nterms = self._get_matrix_row_fmt_nterms_nfloats(
                        nvalues, tout)
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
                    assert len(GCi) == len(reals), 'tout=%s nGCi=%s nreals=%s' % (tout, len(GCi), len(reals))
                else:
                    assert len(GCi)*2 == len(reals), 'tout=%s nGCi=%s nreals=%s' % (tout, len(GCi)*2, len(reals))
                jj += 1
            else:
                nvalues = self.get_marker1(rewind=False)
                assert nvalues == 0, nvalues

                matrix = self._cast_matrix_mat(GCi, GCj, mrows, ncols, reals, tout, dtype)
                m.data = matrix
                if matrix is not None:
                    op2.matrices[table_name.decode('utf-8')] = m
                #nvalues = self.get_marker1(rewind=True)
                return
            itable -= 1
            niter += 1
        raise RuntimeError('MaxIteration: this should never happen; n=%s' % niter_max)

    def _cast_matrix_mat(self, GCi, GCj, mrows, ncols, reals, tout, dtype):
        """helper method for _read_matrix_mat"""
        op2 = self.op2
        #assert max(GCi) <= mrows, 'GCi=%s GCj=%s mrows=%s' % (GCi, GCj, mrows)
        #assert max(GCj) <= ncols, 'GCi=%s GCj=%s ncols=%s' % (GCi, GCj, ncols)

        # we subtract 1 to the indicides to account for Fortran
        GCi = np.array(GCi, dtype='int32') - 1
        GCj = np.array(GCj, dtype='int32') - 1
        try:
            if dtype == '???':
                matrix = None
                self.log.warning('what is the dtype?')
            elif tout in [1, 2]:
                # real
                real_array = np.array(reals, dtype=dtype)
                matrix = scipy.sparse.coo_matrix(
                    (real_array, (GCi, GCj)),
                    shape=(mrows, ncols), dtype=dtype)
                #self.log.info('created %s (real)' % self.table_name)
            elif tout in [3, 4]:
                # complex
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
                #msg = 'created %s (complex)' % self.table_name
                #self.log.debug(msg)
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
                    self.log.info('created %s' % op2.table_name)
                else:
                    self.log.warning('cant reshape because invalid sizes : created %s' %
                                     op2.table_name)

                matrix = real_array
        return matrix

    def _skip_matrix_mat(self):
        """
        Reads a matrix in "standard" form.

        See also
        --------
        read_matrix_mat
        """
        unused_table_name = self._read_table_name(rewind=False, stop_on_failure=True)
        self.read_markers([-1])
        unused_data = self._skip_record()

        self.read_markers([-2, 1, 0])
        unused_data = self._skip_record()

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
                    unused_data = self._skip_block()
                    nvalues = self.get_marker1(rewind=True)
                jj += 1
            else:
                nvalues = self.get_marker1(rewind=False)
                assert nvalues == 0, nvalues
                return
            itable -= 1
            niter += 1
        raise RuntimeError('this should never happen; n=%s' % niter_max)

    def read_matrix(self, table_name):
        """
        general method for reading matrices and MATPOOL matrices

        .. todo:: Doesn't support checking matrices vs. MATPOOLs
        .. todo:: MATPOOLs are disabled because they're not parsed properly
        """
        op2 = self.op2
        i = self.f.tell()
        # if we skip on read_mode=1, we don't get debugging
        # if we just use read_mode=2, some tests fail
        #
        if self.read_mode == 2 and not self.debug_file:
            try:
                self._skip_matrix_mat()  # doesn't work for matpools
            except MemoryError:
                raise
            except(RuntimeError, AssertionError, ValueError):
                self._goto(i)
                self._skip_table(table_name)
            return

        try:
            self._read_matrix_mat()
        except MemoryError:
            raise
        except(RuntimeError, AssertionError, ValueError):
            # read matpool matrix
            self._goto(i)
            try:
                self._read_matrix_matpool()
            except(RuntimeError, AssertionError, ValueError):
                #raise
                self._goto(i)
                self._skip_table(op2.table_name)

    def _read_matrix_matpool(self):
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
        op2 = self.op2
        table_name = self._read_table_name(rewind=False, stop_on_failure=True)
        utable_name = table_name.decode('utf-8')
        self.read_markers([-1])
        data = self._read_record()

        self.read_markers([-2, 1, 0])
        data = self._read_record()

        self.read_markers([-3, 1, 0])
        data = self._read_record()

        #nvalues = len(data) // 4
        assert len(data) % 4 == 0, len(data) / 4.

        header = unpack(self._endian + b'3i 8s 7i', data[:48]) # 48=4*12
        assert header[:3] == (114, 1, 120), 'header[:3]=%s header=%s' % (header[:3], header)

        # ncols_gset is needed for form=9
        #  list of header values:
        #    4:5   matrix name
        #    6     placeholder
        #    7     matrix shape (1=square, 2 or 9 = rectangular, 6=symmetric)
        #    8     input type flag (1=single, 2=double, 3=complex single,
        #                           4=complex double)
        #    9     output type flag (0=precision set by system cell,
        #                            1=single, 2=double, 3=complex single,
        #                            4=complex double)
        #   10     complex flag (0=real/imaginary, >0=magnitude/phase)
        #   11     placeholder
        #   12     number of columns in the G set
        #          (only necessary for matrix shape 9)
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
            fdtype = op2.fdtype
        elif tout == 2:
            dtype = 'float64'
            fdtype = op2.double_dtype
        elif tout == 3:
            dtype = 'complex64'
            fdtype = op2.fdtype
        elif tout == 4:
            dtype = 'complex128'
            fdtype = op2.double_dtype
        else:
            dtype = '???'
            msg = ('matrix_name=%s, junk1=%s, matrix_shape=%s, tin=%s, tout=%s, '
                   'is_phase=%s, junk2=%s, ncols_gset=%s' % (
                       matrix_name, junk1, matrix_shape, tin, tout,
                       is_phase, junk2, ncols_gset))
            self.log.warning(msg)
            raise RuntimeError(msg)

        is_symmetric = matrix_shape == 6
        #is_phase_flag = is_phase > 0

        if tout in [1, 3]:
            # works for float32, complex64
            ints = np.frombuffer(data[48:], dtype=op2.idtype).copy()
            floats = np.frombuffer(data[48:], dtype=op2.fdtype).copy()
            temp_ints = ints
        else:
            # works for float64, complex128
            temp_ints = np.frombuffer(data[48:], dtype=op2.idtype).copy()

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
                real = np.frombuffer(datai, dtype=fdtype)[1::2].copy()

            elif dtype == 'complex128':
                datai = data[48+(istarti*4) : 48+(istopi*4)]

                # iword
                # -----
                #   0    1    3     5   <---- iword
                #   1    1    2     2   <---- nwords
                # (nid, dof, real, imag)
                irow = np.arange(istarti, istopi-1, step=6, dtype='int32')
                assert len(datai) % 8 == 0, len(datai) / 8
                floats = np.frombuffer(datai, dtype=fdtype).copy()

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

        self._cast_matrix_matpool(utable_name, real_imag_array,
                                  col_nids_array, col_dofs_array,
                                  row_nids_array, row_dofs_array,
                                  matrix_shape, dtype, is_symmetric)

    def _cast_matrix_matpool(self, table_name, real_imag_array,
                             col_nids_array, col_dofs_array,
                             row_nids_array, row_dofs_array,
                             matrix_shape, dtype, is_symmetric):
        """helper method for _read_matpool_matrix"""
        op2 = self.op2
        make_matrix_symmetric = op2.apply_symmetry and matrix_shape == 'symmetric'

        # TODO: this is way slower than it should be
        #       because we didn't preallocate the data and the
        #       grids_comp_array_to_index function needs work
        grids1 = col_nids_array
        comps1 = col_dofs_array
        grids2 = row_nids_array
        comps2 = row_dofs_array
        assert len(grids1) == len(comps1), 'ngrids1=%s ncomps1=%s' % (len(grids1), len(comps1))
        assert len(grids1) == len(grids2), 'ngrids1=%s ngrids2=%s' % (len(grids1), len(grids2))
        assert len(comps1) == len(comps2), 'ncomps1=%s ncomps2=%s' % (len(comps1), len(comps2))

        j1, j2, nj1, nj2, nj = grids_comp_array_to_index(
            grids1, comps1, grids2, comps2, make_matrix_symmetric)
        assert len(j1) == len(j2), 'nj1=%s nj2=%s' % (len(j1), len(j2))
        assert len(grids1) == len(real_imag_array), 'ngrids1=%s nreals=%s' % (len(j1), len(real_imag_array))

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
            # and make it a diagonal matrix
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
                        'Matrix %r marked as symmetric does not contain '
                        'symmetric data.  Data will be symmetrized by averaging.' % table_name)
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
        op2.matrices[table_name] = m
        self.log.debug(m)

        self.read_markers([-4, 1, 0])
        data = self._read_record()

        if len(data) == 12:
            self.read_markers([-5, 1, 0, 0])
            return
        raise RuntimeError('failed on _read_matpool_matrix')

    #---------------------------------------------------------------------------

    def _get_marker_n(self, nmarkers):
        """
        Gets N markers

        A marker is a flag that is used.  It's a series of 3 ints (4, n, 4)
        where n changes from marker to marker.

        Parameters
        ----------
        nmarkers : int
            the number of markers to read

        Returns
        -------
        markers : List[int, int, int]
            a list of nmarker integers

        """
        op2 = self.op2
        markers = []
        struc = Struct('3i')
        for unused_i in range(nmarkers):
            block = self.f.read(12)
            marker = struc.unpack(block)
            markers.append(marker)
        return markers

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
        ni = self.n
        markers = []
        for i in range(n):
            data = self.read_block()
            marker, = op2.struct_i.unpack(data)
            markers.append(marker)
        if rewind:
            self.n = ni
            self.f.seek(self.n)
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
            data = self.read_block()
            imarker, = op2.struct_i.unpack(data)
            if marker != imarker:
                import os
                msg = 'marker=%r imarker=%r; markers=%s; i=%s; table_name=%r; iloc=%s/%s' % (
                    marker, imarker, markers, i, op2.table_name,
                    self.f.tell(), os.path.getsize(op2.op2_filename))
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
        month, day, year = self._set_op2_date(month, day, year)

        #self.log.debug("%s/%s/%4i zero=%s one=%s" % (month, day, year, zero, one))
        #if self.is_debug_file:
        if self.is_debug_file:
            self.binary_debug.write('  [subtable_name, month=%i, day=%i, year=%i, '
                                    'zero=%i, one=%i]\n\n' % (month, day, year, zero, one))
        #assert zero == 0, zero  # is this the RTABLE indicator???
        assert one in [0, 1], one  # 0, 50

    def _set_op2_date(self, month, day, year):
        """sets the date the job was run"""
        date = (month, day, 2000 + year)
        self.op2.date = date
        return date

    #----------------------------------------------------------------------------------------
    def _read_record(self, debug=True, macro_rewind=False):
        """
        Reads a record.

        A record is defined N blocks.  Blocks are split every 2^12 bytes,
        which is an oddity of the OP2, which is a "Fortran formatted" file.
        You can think of a block as a partial result.  A record is a full
        result.

        If a block is small enough, it will fit into 2^12 bytes and a record
        is a block.

        """
        return self._read_record_ndata(debug, macro_rewind)[0]

    def _read_record_ndata(self, debug=True, macro_rewind=False):
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
        data = self.f.read(4)
        ndata, = op2.struct_i.unpack(data)

        data_out = self.f.read(ndata)
        data = self.f.read(4)
        self.n += 8 + ndata
        return data_out, ndata

    #------------------------------------------------------------------
    def _read_table_name(self, rewind=False, stop_on_failure=True):
        """Reads the next OP2 table name (e.g. OUG1, OES1X1)"""
        op2 = self.op2
        table_name = None
        data = None
        if self.is_debug_file:
            self.binary_debug.write('_read_table_name - rewind=%s\n' % rewind)
        ni = self.n
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
                self.n = ni
                self.f.seek(self.n)

                try:
                    # we have a trailing 0 marker
                    self.read_markers([0], macro_rewind=rewind)
                except:
                    # if we hit this block, we have a FATAL error
                    if not op2._nastran_format.lower().startswith('imat') and op2.post != -4:
                        self.f.seek(self.n)
                        self.show(1000)
                        raise FatalError('There was a Nastran FATAL Error.  '
                                         'Check the F06.\nlast table=%r; post=%s' % (
                                             self.table_name, op2.post))
                table_name = None

                # we're done reading, so we're going to ignore the rewind
                rewind = False

        if rewind:
            self.n = ni
            self.f.seek(self.n)
        return table_name

    def read_block(self):
        """
        Reads a block following a pattern of:
            [nbytes, data, nbytes]

        Returns
        -------
        data : bytes
            the data in binary

        See also
        --------
        read_3_blocks - reds 3 blocks
        ??? - reads multi-blocks

        """
        op2 = self.op2
        data = self.f.read(4)
        ndata, = op2.struct_i.unpack(data)

        data_out = self.f.read(ndata)
        data = self.f.read(4)
        self.n += 8 + ndata
        return data_out

    def read_3_blocks(self):
        """
        Reads a block following a pattern of:
            [nbytes, data, nbytes]
            [nbytes, data, nbytes]
            [nbytes, data, nbytes]

        This is intended to be used for reading marker triples

        Returns
        -------
        data : bytes
            the data in binary

        """
        op2 = self.op2
        data_out = b''
        for unused_i in range(3):
            data = self.f.read(4)
            ndata, = op2.struct_i.unpack(data)

            data_out += self.f.read(ndata)
            data = self.f.read(4)
            self.n += 8 + ndata
        return data_out

    def read_3_markers(self, markers, macro_rewind=True):
        """
        Micro-optimizes ``read_markers`` for 3 markers.

        Parameters
        ----------
        markers : List[int, int, int]
            markers to get; markers = [-10, 1, 0]

        Raises
        ------
        FortranMarkerError
            if the expected table number is not found

        """
        op2 = self.op2
        data = self.read_3_blocks()
        unused_imarkers = op2.struct_3i.unpack(data)
        for imarker, marker in enumerate(markers):
            if marker != imarker:
                msg = 'marker=%r imarker=%r; markers=%s; table_name=%r' % (
                    marker, imarker, markers, op2.table_name)
                raise FortranMarkerError(msg)
            if self.is_debug_file:
                self.binary_debug.write('  read_markers -> [4, %i, 4]\n' % marker)

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
        ni = self.n
        #markers = []
        data = self.read_block()
        marker, = op2.struct_i.unpack(data)
        if rewind:
            self.n = ni
            self.f.seek(self.n)
        else:
            if self.is_debug_file:
                msg = 'get_marker - [4, %i, 4]; macro_rewind=%s\n' % (
                    marker, macro_rewind or rewind)
                self.binary_debug.write(msg)
        return marker

    #def _read_subtables(self):
        #"""interface to the op2 object"""
        #return self.op2._read_subtables()

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

    @property
    def _uendian(self):
        """interface to the op2 object"""
        return self.op2._uendian

    @property
    def _encoding(self):
        """interface to the op2 object"""
        return self.op2._encoding

    @property
    def read_mode(self):
        """interface to the op2 object"""
        return self.op2.read_mode

    @property
    def debug_file(self):
        """interface to the op2 object"""
        return self.op2.debug_file

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
        self.table_name = self._read_table_name(rewind=False)
        if self.is_debug_file:
            self.binary_debug.write('skipping table...%r\n' % op2.table_name)
        self.log.warning('    skipping table_helper = %s' % op2.table_name)

        self.read_markers([-1])
        unused_data = self._skip_record()
        self.read_markers([-2, 1, 0])
        unused_data = self._skip_record()
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
        unused_markers0 = self.get_nmarkers(1, rewind=False)
        record = self._skip_block()

        markers1 = self.get_nmarkers(1, rewind=True)
        # handling continuation blocks
        while markers1[0] > 0:
            markers1 = self.get_nmarkers(1, rewind=False)
            record = self._skip_block()
            markers1 = self.get_nmarkers(1, rewind=True)
        return record

    def _skip_record_ndata(self, debug=True, macro_rewind=False):
        """the skip version of ``_read_record_ndata``"""
        op2 = self.op2
        marker0 = self.get_marker1(rewind=False, macro_rewind=macro_rewind)
        if self.is_debug_file and debug:
            self.binary_debug.write('read_record - marker = [4, %i, 4]; macro_rewind=%s\n' % (
                marker0, macro_rewind))
        record, nrecord = self._skip_block_ndata()

        if self.is_debug_file and debug:
            self.binary_debug.write('read_record - record = [%i, recordi, %i]; '
                                    'macro_rewind=%s\n' % (nrecord, nrecord, macro_rewind))
        if marker0*4 != nrecord:
            msg = 'marker0=%s*4 len(record)=%s; table_name=%r' % (
                marker0*4, nrecord, op2.table_name)
            raise FortranMarkerError(msg)

        marker1 = self.get_marker1(rewind=True)

        if marker1 > 0:
            while marker1 > 0:
                marker1 = self.get_marker1(rewind=False)
                if self.is_debug_file and debug:
                    self.binary_debug.write('read_record - marker1 = [4, %i, 4]\n' % marker1)
                unused_recordi, nrecordi = self._skip_block_ndata()
                nrecord += nrecordi

                marker1 = self.get_marker1(rewind=True)
                if self.is_debug_file and debug:
                    self.binary_debug.write('read_record - marker1 = [4, %i, 4]\n' % marker1)
        return record, nrecord

    def _get_record_length(self):
        """
        The record length helps us figure out data block size, which is used
        to quickly size the arrays.  We just need a bit of meta data and can
        jump around quickly.

        Returns
        -------
        record_length : int
            the length of the data block

        """
        op2 = self.op2
        if self.is_debug_file:
            self.binary_debug.write('_get_record_length\n')
        len_record = 0
        n0 = self.n
        markers0 = self.get_nmarkers(1, rewind=False)
        if self.is_debug_file:
            self.binary_debug.write('  markers0=%s\n' % markers0)

        n = self.n
        unused_record = self._skip_block()
        len_record += self.n - n - 8  # -8 is for the block
        if self.is_debug_file:
            self.binary_debug.write('  len_record=%s\n' % len_record)

        markers1 = self.get_nmarkers(1, rewind=True)
        # handling continuation blocks
        while markers1[0] > 0:
            markers1 = self.get_nmarkers(1, rewind=False)
            if self.is_debug_file:
                self.binary_debug.write('  markers1=%s\n' % markers1)
            n = self.n
            unused_record = self._skip_block()
            len_record += self.n - n - 8  # -8 is for the block
            markers1 = self.get_nmarkers(1, rewind=True)
        self._goto(n0)
        return len_record

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
        data = self.f.read(4)
        ndata, = op2.struct_i.unpack(data)
        self.n += 8 + ndata
        self._goto(self.n)
        return None, ndata

    #---------------------------------------------------------------------------
    def _goto(self, n):
        """
        Jumps to position n in the file

        Parameters
        ----------
        n : int
            the position to goto

        """
        self.n = n
        self.f.seek(n)

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
            if hasattr(op2, 'isubcase') and self.isubcase in op2.valid_subcases:
                return True
            return False
        return True

    def read_results_table(self):
        """Reads a results table"""
        #op2 = self.op2
        if self.is_debug_file:
            self.binary_debug.write('read_results_table - %s\n' % self.table_name)
        self.table_name = self._read_table_name(rewind=False)
        self.read_markers([-1])
        if self.is_debug_file:
            self.binary_debug.write('---markers = [-1]---\n')
            #self.binary_debug.write('marker = [4, -1, 4]\n')
        data = self._read_record()

        self.read_markers([-2, 1, 0])
        if self.is_debug_file:
            self.binary_debug.write('---markers = [-2, 1, 0]---\n')
        data, ndata = self._read_record_ndata()
        if ndata == 8:
            subtable_name = op2.struct_8s.unpack(data)
            if self.is_debug_file:
                self.binary_debug.write('  recordi = [%r]\n'  % subtable_name)
                self.binary_debug.write('  subtable_name=%r\n' % subtable_name)
        elif ndata == 28:
            subtable_name, month, day, year, zero, one = unpack(self._endian + b'8s5i', data)
            if self.is_debug_file:
                self.binary_debug.write('  recordi = [%r, %i, %i, %i, %i, %i]\n'  % (
                    subtable_name, month, day, year, zero, one))
                self.binary_debug.write('  subtable_name=%r\n' % subtable_name)
            self._print_month(month, day, year, zero, one)
        elif ndata == 612: # ???
            strings, ints, floats = self.show_data(data)
            msg = 'len(data) = %i\n' % ndata
            #msg += 'strings  = %r\n' % strings
            #msg += 'ints     = %r\n' % str(ints)
            #msg += 'floats   = %r' % str(floats)
            print(msg)
            subtable_name, = op2.struct_8s.unpack(data[:8])
            print('subtable_name = %r' % subtable_name.strip())
        else:
            strings, ints, floats = self.show_data(data)
            msg = 'len(data) = %i\n' % ndata
            msg += 'strings  = %r\n' % strings
            msg += 'ints     = %r\n' % str(ints)
            msg += 'floats   = %r' % str(floats)
            raise NotImplementedError(msg)
        if hasattr(self, 'subtable_name'):
            raise RuntimeError('the file hasnt been cleaned up; subtable_name_old=%s new=%s' % (
                self.subtable_name, subtable_name))
        self.subtable_name = subtable_name
        self._read_subtables()

    def read_geom_table(self):
        """Reads a geometry table"""
        op2 = self.op2
        op2.table_name = self._read_table_name(rewind=False)
        if self.is_debug_file:
            self.binary_debug.write('read_geom_table - %s\n' % op2.table_name)
        self.read_markers([-1])
        data = self._read_record()

        self.read_markers([-2, 1, 0])
        data, ndata = self._read_record_ndata()
        if ndata == 8:
            subtable_name, = op2.struct_8s.unpack(data)
        else:
            strings, ints, floats = self.show_data(data)
            msg = 'Unhandled table length error\n'
            msg += 'table_name = %s\n' % op2.table_name
            msg += 'len(data) = %i\n' % ndata
            msg += 'strings  = %r\n' % strings
            msg += 'ints     = %r\n' % str(ints)
            msg += 'floats   = %r' % str(floats)
            raise NotImplementedError(msg)

        self.subtable_name = subtable_name.rstrip()
        self._read_subtables()

    def _read_subtables(self):
        """reads a series of subtables"""
        # this parameters is used for numpy streaming
        op2 = self.op2
        op2._table4_count = 0
        op2.is_table_1 = True
        op2._data_factor = 1

        #nstart = self.n
        op2.isubtable = -3
        self.read_markers([-3, 1, 0])
        if self.is_debug_file:
            self.binary_debug.write('***isubtable = %i\n' % op2.isubtable)
            self.binary_debug.write('---markers = [-3, 1, 0]---\n')
        table_mapper = op2._get_table_mapper()

        # get the parsing functions (table3_parser, table4_parser)
        # or find out we're going to be skipping the tables
        #
        # table3 - the table with the meta data (e.g. subcase_id, time, is_stress/strain)
        # table4 - the actual results data
        #
        # we indicate table3/4 by isubtable, which starts from -3 (so table3) and counts
        # down (yes down) to 4 to indicate table4.  If we count down again, we end up
        # back at table 3 (with isubtable=-5), which will occur in the case of multiple
        # times/element types/results in a single macro table (e.g. OUG, OES).
        if op2.table_name in table_mapper:
            #if self.read_mode == 2:
                #self.log.debug("table_name = %r" % op2.table_name)
            table3_parser, table4_parser = table_mapper[op2.table_name]
            passer = False
        else:
            if self.read_mode == 2:
                self.log.info("skipping table_name = %r" % op2.table_name)
                    #raise NotImplementedError(op2.table_name)
            table3_parser = None
            table4_parser = None
            passer = True

        # we need to check the marker, so we read it and rewind, so we don't
        # screw up our positioning in the file
        markers = self.get_nmarkers(1, rewind=True)
        if self.is_debug_file:
            self.binary_debug.write('---marker0 = %s---\n' % markers)

        # while the subtables aren't done
        while markers[0] != 0:
            op2.is_start_of_subtable = True
            if self.is_debug_file:
                self.binary_debug.write('***isubtable = %i\n' % op2.isubtable)
            self._read_subtable_3_4(table3_parser, table4_parser, passer)
            #force_table4 = self._read_subtable_3_4(table3_parser, table4_parser, passer)
            op2.isubtable -= 1
            self.read_markers([op2.isubtable, 1, 0])
            markers = self.get_nmarkers(1, rewind=True)
        if self.is_debug_file:
            self.binary_debug.write('breaking on marker=%r\n' % str(markers))

        # we've finished reading all subtables, but have one last marker to read
        self.read_markers([0])
        self._finish()

    def _read_subtable_3_4(self, table3_parser, table4_parser, passer):
        """
        Reads a series of subtable 3/4

        Parameters
        ----------
        table3_parser : bool / function
            None : just to break the code if we're on the array sizing step
            function : the table 3 reading function
        table4_parser :bool / function
            None : just to break the code if we're on the array sizing step
            function : the table 4 reading function
        passer : bool
            flag to see if we're skipping tables

        Returns
        -------
        flag : bool
            True : ???
            False : failed???
            None : passed???

        """
        if self.binary_debug:
            self.binary_debug.write('-' * 60 + '\n')
        # this is the length of the current record inside table3/table4
        record_len = self._get_record_length()
        if self.is_debug_file:
            self.binary_debug.write('record_length = %s\n' % record_len)

        oes_nl = [b'OESNLXD', b'OESNL1X', b'OESNLXR']
        if record_len == 584:  # table3 has a length of 584
            if self.table_name in oes_nl and hasattr(self, 'num_wide') and self.num_wide == 146:
                data_code_old = deepcopy(self.data_code)

            self.data_code = {'_encoding' : self._encoding}
            self.obj = None
            data, ndata = self._read_record_ndata()
            if not passer:
                try:
                    table3_parser(data, ndata)
                except SortCodeError:
                    if self.is_debug_file:
                        self.binary_debug.write('except SortCodeError!\n')
                    if self.table_name in oes_nl:
                        self.data_code = data_code_old
                        for key, value in data_code_old.items():
                            setattr(self, key, value)
                        table4_parser(data, ndata)
                        return False
                    raise RuntimeError(self.code_information())
                #if hasattr(self, 'isubcase'):
                    #print("code = ", self._get_code())
        else:
            if passer or not self.is_valid_subcase():
                data = self._skip_record()
            else:
                if hasattr(self, 'num_wide'):
                    # num_wide is the result size and is usually found in
                    # table3, but some B-list tables don't have it
                    unused_n = self._read_subtable_results(table4_parser, record_len)
                else:
                    data, ndata = self._read_record_ndata()
                    unused_n = table4_parser(data, ndata)
                #del n

    def _read_subtable_results(self, table4_parser, record_len):
        """
        # if reading the data
        # 1 - 1st pass to size the array (vectorized)
        # 2 - 2nd pass to read the data  (vectorized)

        Parameters
        ----------
        table4_parser : function
            the parser function for table 4
        record_len : int
            the length of the record block

        Returns
        -------
        n : None / int
            None : an error occurred or we're in read_mode=1/array sizeing (???)
            int : the number of bytes that have been read

        """
        datai = b''
        n = 0
        if self.read_mode == 2:
            self.ntotal = 0

            data, ndata = self._read_record_ndata()
            n = table4_parser(data, ndata)
            assert isinstance(n, integer_types), self.table_name

            self._reset_vector_counter()

        elif self.read_mode == 1:
            # if we're checking the array size

            #n = self._skip_record()
            #n = table4_parser(datai, 300000)
            if self.table_name in [b'R1TABRG', b'ONRGY1']:
                data, ndata = self._read_record_ndata()
            else:
                data, ndata = self._skip_record_ndata()
            n = table4_parser(data, ndata)
            if not isinstance(n, integer_types):
                msg = 'n is not an integer; table_name=%s n=%s table4_parser=%s' % (
                    self.table_name, n, table4_parser)
                raise TypeError(msg)

            #self._goto(n)
            #n = self._skip_record()

            self._init_vector_counter(record_len)
        else:
            raise RuntimeError(self.read_mode)
        self._cleanup_data_members()
        return n

    def _init_vector_counter(self, record_len):
        """
        Sets the table size

        Parameters
        ----------
        record_len : int
            the length of the record block
        """
        if not(hasattr(self, 'obj') and self.obj is not None):
            return

        if hasattr(self.obj, 'ntimes'):
            if not hasattr(self.obj, '_reset_indices'):
                #methods = '\ndir(obj)=%s' % ', '.join(sorted(dir(self.obj)))
                #msg = 'is %s vectorized because its missing _reset_indices...%s' % (
                    #self.obj.__class__.__name__, methods)
                return None
                #raise RuntimeError(msg)
            self.obj._reset_indices()
            self.obj.ntimes += 1
            ntotal = record_len // (self.num_wide * 4) * self._data_factor

            # this has a problem with XYPLOT data if there is a result
            #    request in the same format (e.g. OESNLXD/OES1X1 tables
            #    if they both have the same element ID)
            #
            #class_name = self.obj.__class__.__name__
            #if class_name == 'RealBush1DStressArray':
                #print('%s.ntotal = %s' % (class_name, ntotal))
                #print('num_wide=%s factor=%s len=%s ntotal=%s' % (
                    #self.num_wide, self._data_factor, record_len, ntotal))
            self.obj.ntotal = ntotal
            self.obj._ntotals.append(ntotal)

            assert isinstance(self.obj.ntotal, integer_types), type(self.obj.ntotal)
        else:
            self.log.warning('obj=%s doesnt have ntimes' % self.obj.__class__.__name__)

    def _cleanup_data_members(self):
        """deletes variables from previous tables"""
        del_words = [
            'words',
            #'Title',
            #'ID',
            'analysis_code',
            #'result_names',
            #'labels',
            #'data_names',
        ]
        msg = ''
        if hasattr(self, 'words'):
            if not len(self.words) in [0, 28]:
                msg = 'table_name=%r len(self.words)=%s words=%s' % (
                    self.table_name, len(self.words), self.words)
                raise RuntimeError(msg)

            for word in self.words:
                if word in ['???', 'Title']:
                    continue
                if not hasattr(self, word):
                    continue
                delattr(self, word)
            self.words = []
        if hasattr(self, 'analysis_code'):
            del self.analysis_code
        #if hasattr(self, 'data_names') and self.data_names is not None:
            #print(object_attributes(self))

        if hasattr(self, 'data_code'):
            del self.data_code

        for word in del_words:
            if hasattr(self, word):
                val = getattr(self, word)
                if isinstance(val, list) and len(val) == 0:
                    continue
                msg += '  %s=%s\n' % (word, val)
        if msg:
            print(object_attributes(self))
            print(msg)

    def show(self, n, types='ifs', endian=None):  # pragma: no cover
        assert self.n == self.f.tell()
        nints = n // 4
        data = self.f.read(4 * nints)
        strings, ints, floats = self.show_data(data, types=types, endian=endian)
        self.f.seek(self.n)
        return strings, ints, floats

    def show_data(self, data, types='ifs', endian=None):  # pragma: no cover
        """
        Shows a data block as various types

        Parameters
        ----------
        data : bytes
            the binary string bytes
        types : str; default='ifs'
            i - int
            f - float
            s - string
            d - double (float64; 8 bytes)
            q - long long (int64; 8 bytes)

            l - long (int; 4 bytes)
            I - unsigned int (int; 4 bytes)
            L - unsigned long (int; 4 bytes)
            Q - unsigned long long (int; 8 bytes)
        endian : str; default=None -> auto determined somewhere else in the code
            the big/little endian {>, <}

        .. warning:: 's' is apparently not Python 3 friendly

        """
        return self._write_data(sys.stdout, data, types=types, endian=endian)

    def _write_data(self, f, data, types='ifs', endian=None):  # pragma: no cover
        """
        Useful function for seeing what's going on locally when debugging.

        Parameters
        ----------
        data : bytes
            the binary string bytes
        types : str; default='ifs'
            i - int
            f - float
            s - string
            d - double (float64; 8 bytes)
            q - long long (int64; 8 bytes)

            l - long (int; 4 bytes)
            I - unsigned int (int; 4 bytes)
            L - unsigned long (int; 4 bytes)
            Q - unsigned long long (int; 8 bytes)
        endian : str; default=None -> auto determined somewhere else in the code
            the big/little endian {>, <}

        """
        n = len(data)
        nints = n // 4
        ndoubles = n // 8
        strings = None
        ints = None
        floats = None
        longs = None

        if endian is None:
            endian = self._uendian
            assert endian is not None, endian

        f.write('\nndata = %s:\n' % n)
        for typei in types:
            assert typei in 'sifdq lIL', 'type=%r is invalid' % typei

        if 's' in types:
            strings = unpack('%s%is' % (endian, n), data)
            f.write("  strings = %s\n" % str(strings))
        if 'i' in types:
            ints = unpack('%s%ii' % (endian, nints), data)
            f.write("  ints    = %s\n" % str(ints))
        if 'f' in types:
            floats = unpack('%s%if' % (endian, nints), data)
            f.write("  floats  = %s\n" % str(floats))
        if 'd' in types:
            doubles = unpack('%s%id' % (endian, ndoubles), data[:ndoubles*8])
            f.write("  doubles (float64) = %s\n" % str(doubles))

        if 'l' in types:
            longs = unpack('%s%il' % (endian, nints), data)
            f.write("  long  = %s\n" % str(longs))
        if 'I' in types:
            ints2 = unpack('%s%iI' % (endian, nints), data)
            f.write("  unsigned int = %s\n" % str(ints2))
        if 'L' in types:
            longs2 = unpack('%s%iL' % (endian, nints), data)
            f.write("  unsigned long = %s\n" % str(longs2))
        if 'q' in types:
            longs = unpack('%s%iq' % (endian, ndoubles), data[:ndoubles*8])
            f.write("  long long (int64) = %s\n" % str(longs))
        f.write('\n')
        return strings, ints, floats

    def show_ndata(self, n, types='ifs'):  # pragma: no cover
        return self._write_ndata(sys.stdout, n, types=types)

    def _write_ndata(self, f, n, types='ifs'):  # pragma: no cover
        """
        Useful function for seeing what's going on locally when debugging.
        """
        nold = self.n
        data = self.f.read(n)
        self.n = nold
        self.f.seek(self.n)
        return self._write_data(f, data, types=types)

    def _finish(self):
        """
        Clears out the data members contained within the self.words variable.
        This prevents mixups when working on the next table, but otherwise
        has no effect.

        """
        if hasattr(self, 'words'):
            for word in self.words:
                if word != '???' and hasattr(self, word):
                    if word not in ['Title', 'reference_point']:
                        delattr(self, word)
        self.obj = None
        if hasattr(self, 'subtable_name'):
            del self.subtable_name

    #-------------------------------------------------------------------------------------------
    # oug-style reader
    def create_transient_object(self, storage_obj, class_obj, is_cid=False, debug=False):
        """
        Creates a transient object (or None if the subcase should be skippied).

        Parameters
        ----------
        storageName : str
            the name of the dictionary to store the object in (e.g. 'displacements')
        class_obj : object()
            the class object to instantiate
        debug : bool
            developer debug

        .. python ::

            slot = self.displacements
            slot_vector = RealDisplacementArray
            self.create_transient_object(slot, slot_vector, is_cid=is_cid)

        .. note:: dt can also be load_step depending on the class

        """
        op2 = self.op2
        assert not isinstance(class_obj, string_types), 'class_obj=%r' % class_obj
        assert class_obj is not None, class_obj
        if debug:
            print("create Transient Object")
            print("***NF = %s" % self.nonlinear_factor)
        #if not hasattr(self, storageName):
            #attrs =  object_attributes(obj, mode="public")
            #msg = 'storage_obj=%r does not exist.\n' % storage_obj
            #msg += 'Attributes = [%s]' , ', %s'.join(attrs)
            #raise RuntimeError(msg)
        #storage_obj = getattr(self, storageName)
        #assert class_obj is not None, 'name=%r has no associated classObject' % storageName

        #self.log.debug('self.table_name=%s isubcase=%s subtitle=%r' % (
            #self.table_name, self.isubcase, self.subtitle.strip()))
        self.data_code['table_name'] = self.table_name.decode(self.encoding)
        assert self.log is not None

        code = self._get_code()
        #print('code =', code)
        if hasattr(self, 'isubcase'):
            if self.code in storage_obj:
                print('made %s' % str(self.code))
                obj = storage_obj[code]
                if self.nonlinear_factor is not None:
                    if obj.nonlinear_factor is None:
                        msg = (
                            'The object is flipping from a static (e.g. preload)\n'
                            'result to a transient/frequency based results\n'
                            '%s -> %s\n' % (obj.nonlinear_factor, self.nonlinear_factor))
                        msg += (
                            'code = (subcase=%s, analysis_code=%s, sort=%s, count=%s, ogs=%s, '
                            'superelement_adaptivity_index=%r pval_step=%r)\n' % tuple(code))
                        msg += '%s\n' % str(obj)
                        msg += '\nIf this isnt correct, check if the data code was applied on the object'
                        raise MultipleSolutionNotImplementedError(msg)
                obj.update_data_code(copy.deepcopy(self.data_code))
            else:
                print('making %s' % class_obj)
                class_obj.is_cid = is_cid
                is_sort1 = self.is_sort1  # uses the sort_bits

                obj = class_obj(self.data_code, is_sort1, self.isubcase, self.nonlinear_factor)
            self.obj = obj
            storage_obj[code] = obj
        else:
            if code in storage_obj:
                self.obj = storage_obj[code]
            else:
                storage_obj[code] = self.obj

    def _create_table_vector(self, result_name, nnodes,
                             slot, slot_vector, is_cid=False):
        assert isinstance(result_name, string_types), result_name
        assert isinstance(slot, dict), slot
        auto_return = False
        print('read_mode=%s %s nnodes=%s' % (self.read_mode, result_name, nnodes))
        self.result_names.add(result_name)
        assert nnodes > 0

        if self.read_mode == 1:
            self.create_transient_object(slot, slot_vector, is_cid=is_cid)
            self.result_names.add(result_name)
            self.obj._nnodes += nnodes
            auto_return = True
        elif self.read_mode == 2:
            self.code = self._get_code()
            obj = slot[self.code]
            print(obj)
            #self.obj.update_data_code(self.data_code)
            obj.build()
            self.obj = obj
        else:
            auto_return = True
        return auto_return

    def _fix_format_code(self, format_code=1):
        """
        Nastran can mess up the format code by using what the user specified,
        which may be wrong.

        For a SOL 101, if the user uses the following in their BDF:
            DISP(PLOT,PHASE)=ALL
        it's wrong, and should be:
            DISP(PLOT,REAL)=ALL

        """
        # we'll probably remove this later because we're fixing
        #it before we get to the object
        return
        # if self.format_code != format_code:
            # self.format_code = format_code
            # self.obj.format_code = format_code
            # self.obj.data_code['format_code'] = format_code

    def create_transient_object(self, storage_obj, class_obj, is_cid=False, debug=False):
        """
        Creates a transient object (or None if the subcase should be skippied).

        Parameters
        ----------
        storageName : str
            the name of the dictionary to store the object in (e.g. 'displacements')
        class_obj : object()
            the class object to instantiate
        debug : bool
            developer debug

        .. python ::

            slot = self.displacements
            slot_vector = RealDisplacementArray
            self.create_transient_object(slot, slot_vector, is_cid=is_cid)

        .. note:: dt can also be load_step depending on the class

        """
        print("create Transient Object")
        assert not isinstance(class_obj, string_types), 'class_obj=%r' % class_obj
        assert class_obj is not None, class_obj
        if debug:
            print("create Transient Object")
            print("***NF = %s" % self.nonlinear_factor)
        #if not hasattr(self, storageName):
            #attrs =  object_attributes(obj, mode="public")
            #msg = 'storage_obj=%r does not exist.\n' % storage_obj
            #msg += 'Attributes = [%s]' , ', %s'.join(attrs)
            #raise RuntimeError(msg)
        #storage_obj = getattr(self, storageName)
        #assert class_obj is not None, 'name=%r has no associated classObject' % storageName

        #self.log.debug('self.table_name=%s isubcase=%s subtitle=%r' % (
            #self.table_name, self.isubcase, self.subtitle.strip()))
        self.data_code['table_name'] = self.table_name.decode(self.encoding)
        assert self.log is not None

        code = self._get_code()
        #print('code =', code)
        if hasattr(self, 'isubcase'):
            if self.code in storage_obj:
                obj = storage_obj[code]
                if self.nonlinear_factor is not None:
                    if obj.nonlinear_factor is None:
                        msg = ('The object is flipping from a static (e.g. preload)\n'
                               'result to a transient/frequency based results\n'
                               '%s -> %s\n' % (obj.nonlinear_factor, self.nonlinear_factor))
                        msg += (
                            'code = (subcase=%s, analysis_code=%s, sort=%s, count=%s, ogs=%s, '
                            'superelement_adaptivity_index=%r pval_step=%r)\n' % tuple(code))
                        msg += '%s\n' % str(obj)
                        msg += '\nIf this isnt correct, check if the data code was applied on the object'
                        raise MultipleSolutionNotImplementedError(msg)
                obj.update_data_code(copy.deepcopy(self.data_code))
            else:
                class_obj.is_cid = is_cid
                is_sort1 = self.is_sort1  # uses the sort_bits
                obj = class_obj(self.data_code, is_sort1, self.isubcase, self.nonlinear_factor)

            self.obj = obj
            storage_obj[code] = obj
        else:
            if code in storage_obj:
                self.obj = storage_obj[code]
            else:
                storage_obj[code] = self.obj

    def _get_code(self):
        code = self.isubcase
        ogs = 0
        if hasattr(self, 'ogs'):
            ogs = self.ogs
        #if self.binary_debug:
            #self.binary_debug.write(self.code_information(include_time=True))

        code = (self.isubcase, self.analysis_code, self._sort_method, self._count, ogs,
                self.superelement_adaptivity_index, self.pval_step)
        #code = (self.isubcase, self.analysis_code, self._sort_method, self._count,
                #self.superelement_adaptivity_index, self.table_name_str)
        #print('%r' % self.subtitle)
        self.code = code
        #self.log.debug('code = %s' % str(self.code))
        return self.code

    @property
    def _sort_method(self):
        try:
            sort_method, is_real, is_random = self._table_specs()
        except:
            sort_method = get_sort_method_from_table_name(self.table_name)
        #is_sort1 = self.table_name.endswith('1')
        #is_sort1 = self.is_sort1  # uses the sort_bits
        assert sort_method in [1, 2], 'sort_method=%r\n%s' % (sort_method, self.code_information())
        return sort_method

    def update_mode_cycle(self, name):
        value = getattr(self, name)
        if value == 0.0:
            #print('table_name=%r mode=%s eigr=%s' % (self.table_name, self.mode, self.eigr))
            value = np.sqrt(np.abs(self.eign)) / (2. * np.pi)
            setattr(self, name, value)
            self.data_code[name] = value

    def _get_code(self):
        code = self.isubcase
        ogs = 0
        if hasattr(self, 'ogs'):
            ogs = self.ogs
        #if self.binary_debug:
            #self.binary_debug.write(self.code_information(include_time=True))

        code = (self.isubcase, self.analysis_code, self._sort_method, self._count, ogs,
                self.superelement_adaptivity_index, self.pval_step)
        #code = (self.isubcase, self.analysis_code, self._sort_method, self._count,
                #self.superelement_adaptivity_index, self.table_name_str)
        #print('%r' % self.subtitle)
        self.code = code
        #self.log.debug('code = %s' % str(self.code))
        return self.code

    def _reset_vector_counter(self):
        """
        if reading the data
        0 - non-vectorized
        1 - 1st pass to size the array (vectorized)
        2 - 2nd pass to read the data  (vectorized)

        vectorized objects are stored as self.obj
        they have obj.itime which is their table3 counter
        """
        if not(hasattr(self, 'obj') and hasattr(self.obj, 'itime')):
            #print('self.obj.name=%r doesnt have itime' % self.obj.__class__.__name__)
            return
        #ntotal = record_len // (self.num_wide * 4) * self._data_factor

        # we reset the itime counter when we fill up the
        # total number of nodes/elements/layers in the
        # result, where ntotal is the critical length of
        # interest.  This let's us start back at the correct
        # spot the next time we read table3
        #
        # For displacements, ntotal=nnodes
        #
        # For a CBAR, it's ntotal=nelements*2, where 2 is
        # the number of nodes; points A/B
        #
        # For a CTRIA3 / linear CQUAD4, it's
        # ntotal=nelements*2, where 2 is the number of
        # layers (top/btm) and we only get a centroidal
        # result.
        #
        # For a CQUAD4 bilinear, it's
        # ntotal=nelements*(nnodes+1)*2, where 2 is the
        # number of layers and nnodes is 4 (we get an extra
        # result at the centroid).
        #
        # For a PCOMP, it's ntotal=sum(nelements*nlayers),
        # where each element can have a different number
        # of layers
        obj = self.obj
        if obj.ntotal == obj.data.shape[1]:
            #if self.table_name_str in ['OESRMS2', 'OESNO2', 'OSTRRMS2', 'OSTRNO2', 'OESATO2']:
                #print('resetting %r indicies; itime=%s; shape=%s' % (
                    #self.obj.class_name, self.obj.itime, self.obj.data.shape))
            obj._reset_indices()
            obj.words = self.words
            obj.itime += 1
        else:
            # This happens when self._data_factor hasn't been reset
            # or is set wrong.
            # can it happen any other time?
            msga = 'self.obj.name=%r has itime' % obj.__class__.__name__
            self.log.debug(msga)
            msgb = 'ntotal=%s shape=%s shape[1]=%s _data_factor=%s\n' % (
                obj.ntotal, str(obj.data.shape),
                obj.data.shape[1], self._data_factor)
            msgb += 'obj._ntotals=%s' % obj._ntotals
            self.log.error(msgb)
            raise RuntimeError(msga + '\n' + msgb)

    def _read_table_vectorized(self, data, ndata, result_name, storage_obj,
                               real_vector, complex_vector,
                               node_elem, random_code=None, is_cid=False):
        """Reads a generalized real/complex SORT1/SORT2 table"""
        return self.table_reader._read_table_vectorized(
            data, ndata, result_name, storage_obj,
            real_vector, complex_vector,
            node_elem, random_code=random_code, is_cid=is_cid)


class TableReader(object):
    def __init__(self, op2_reader):
        self.op2_reader = op2_reader

    @property
    def _endian(self):
        """interface to the op2 object"""
        return self.op2_reader._endian

    @property
    def is_debug_file(self):
        """interface to the op2 object"""
        return self.op2_reader.is_debug_file

    @property
    def binary_debug(self):
        """interface to the op2 object"""
        return self.op2_reader.binary_debug

    @property
    def log(self):
        """interface to the op2 object"""
        return self.op2_reader.log

    def _read_table_vectorized(self, data, ndata, result_name, storage_obj,
                               real_vector, complex_vector,
                               node_elem, random_code=None, is_cid=False):
        """Reads a generalized real/complex SORT1/SORT2 table"""
        op2_reader = self.op2_reader
        assert isinstance(result_name, string_types), 'result_name=%r' % result_name
        assert isinstance(storage_obj, dict), 'storage_obj=%r' % storage_obj
        #print('self.num_wide =', self.num_wide)
        #print('random...%s' % self.isRandomResponse())
        #if not self.isRandomResponse():
        is_vectorized = True
        if op2_reader.format_code == 1 and op2_reader.num_wide == 8:  # real/random
            # real
            nnodes = ndata // 32  # 8*4
            auto_return = op2_reader._create_table_vector(
                result_name, nnodes, storage_obj, real_vector, is_cid=is_cid)
            if auto_return:
                return ndata

            op2_reader._fix_format_code(format_code=1)
            if op2_reader.is_sort1:
                if op2_reader.nonlinear_factor is None:
                    n = self._read_real_table_static(data, is_vectorized, nnodes, result_name, node_elem, is_cid=is_cid)
                else:
                    n = self._read_real_table_sort1(data, is_vectorized, nnodes, result_name, node_elem, is_cid=is_cid)
            else:
                n = self._read_real_table_sort2(data, is_vectorized, nnodes, result_name, node_elem, is_cid=is_cid)
                #n = ndata
                #msg = self.code_information()
                #n = self._not_implemented_or_skip(data, ndata, msg)
        elif op2_reader.format_code in [2, 3] and op2_reader.num_wide == 14:  # real or real/imaginary or mag/phase
            # complex
            nnodes = ndata // 56  # 14*4
            if op2_reader.is_debug_file:
                op2_reader.binary_debug.write('nnodes=%s' % nnodes)
            auto_return = self._create_table_vector(
                result_name, nnodes, storage_obj, complex_vector)
            if auto_return:
                return ndata
            if op2_reader.is_sort1:
                if op2_reader.is_magnitude_phase():
                    n = self._read_complex_table_sort1_mag(data, is_vectorized, nnodes, result_name, node_elem)
                else:
                    n = self._read_complex_table_sort1_imag(data, is_vectorized, nnodes, result_name, node_elem)
            else:
                if op2_reader.is_magnitude_phase():
                    n = self._read_complex_table_sort2_mag(data, is_vectorized, nnodes, result_name, node_elem)
                else:
                    n = self._read_complex_table_sort2_imag(data, is_vectorized, nnodes, result_name, node_elem)
                #msg = self.code_information()
                #n = self._not_implemented_or_skip(data, ndata, msg)
        else:
            #msg = 'COMPLEX/PHASE is included in:\n'
            #msg += '  DISP(PLOT)=ALL\n'
            #msg += '  but the result type is REAL\n'
            msg = op2_reader.code_information()
            n = op2_reader._not_implemented_or_skip(data, ndata, msg)
        #else:
        #msg = 'invalid random_code=%s num_wide=%s' % (random_code, op2_reader.num_wide)
        #n = self._not_implemented_or_skip(data, ndata, msg)
        return n

    def _read_real_table_static(self, data, is_vectorized, nnodes, result_name,
                                flag, is_cid=False):
        """
        With a static (e.g. SOL 101) result, reads a complex OUG-style
        table created by:
              DISP(PLOT,SORT1,REAL) = ALL

        """
        if self.is_debug_file:
            self.binary_debug.write('  _read_real_table_static\n')
        assert flag in ['node', 'elem'], flag
        op2_reader = self.op2_reader
        dt = op2_reader.nonlinear_factor
        assert op2_reader.obj is not None
        obj = op2_reader.obj

        if op2_reader.use_vector and is_vectorized:
            n = nnodes * 4 * 8
            itotal2 = obj.itotal + nnodes
            #print('ndata=%s n=%s nnodes=%s' % (ndata, n, nnodes))
            ints = frombuffer(data, dtype=op2_reader.idtype).reshape(nnodes, 8)
            floats = frombuffer(data, dtype=op2_reader.fdtype).reshape(nnodes, 8)
            obj._times[obj.itime] = dt
            #self.node_gridtype[self.itotal, :] = [node_id, grid_type]
            #self.data[self.itime, self.itotal, :] = [v1, v2, v3, v4, v5, v6]
            #obj.node_gridtype[obj.itotal:itotal2, :] = ints[:, 0:1]
            nids = ints[:, 0] // 10
            assert nids.min() > 0, nids.min()
            obj.node_gridtype[obj.itotal:itotal2, 0] = nids
            obj.node_gridtype[obj.itotal:itotal2, 1] = ints[:, 1].copy()
            obj.data[obj.itime, obj.itotal:itotal2, :] = floats[:, 2:].copy()
            obj.itotal = itotal2
        else:
            n = 0
            dt = None
            s = Struct(self._endian + b'2i6f')
            for inode in range(nnodes):
                out = s.unpack(data[n:n+32])
                (eid_device, grid_type, tx, ty, tz, rx, ry, rz) = out
                eid = eid_device // 10
                if self.is_debug_file:
                    self.binary_debug.write('  %s=%i; %s\n' % (flag, eid, str(out)))
                obj.add_sort1(None, eid, grid_type, tx, ty, tz, rx, ry, rz)
                n += 32
        return n

