from struct import unpack
from struct import error as StructError

from pyNastran.f06.f06 import FatalError
from pyNastran.f06.tables.grid_point_weight import GridPointWeight

from pyNastran.op2.dev.oef import OEF
from pyNastran.op2.dev.oes import OES

from pyNastran.op2.dev.opg import OPG
from pyNastran.op2.dev.oqg import OQG
from pyNastran.op2.dev.oug import OUG
from pyNastran.op2.dev.ogpwg import OGPWG
from pyNastran.op2.dev.results import Results
from pyNastran.op2.dev.fortran_format import FortranFormat

class OP2(OEF, OES, OPG, OQG, OUG, OGPWG, FortranFormat, Results):
    def set_subcases(self, isubcases):
        pass
    def get_op2_stats(self):
        pass
    def write_bdf(self, *args, **kwargs):
        pass
    def write_f06(self, *args, **kwargs):
        pass
    def __init__(self, op2_filename, make_geom=False, debug=False, log=None):
        if op2_filename:
            self.op2_filename = op2_filename

        #self.tables_to_read = []

        #BDF.__init__(self, debug=debug, log=log)
        #F06Writer.__init__(self)

        OEF.__init__(self)
        OES.__init__(self)

        OPG.__init__(self)
        OQG.__init__(self)
        OUG.__init__(self)
        OGPWG.__init__(self)
        FortranFormat.__init__(self)
        Results.__init__(self)

        self.grid_point_weight = GridPointWeight()
        self.debug = True
        if self.debug:
            self.binary_debug = open('debug.out', 'wb')

        self.show_table3_map = [
            #'OUGV1',
            #'OEF1X',
            #'OES1X1',
        ]
        self.show_table4_map = [
            #'OUGV1',
            #'OEF1X',
            #'OES1X1',
        ]

    def debug3(self):
        return True
        if self.debug and self.table_name in self.show_table3_map:
            return True
        return False

    def debug4(self):
        return True
        if self.debug and self.table_name in self.show_table4_map:
            return True
        return False

    def read_op2(self, op2_filename=None):
        if op2_filename:
            self.op2_filename = op2_filename
        self.n = 0
        self.table_name = None
        self.f = open(self.op2_filename, 'rb')

        markers = self.get_nmarkers(1, rewind=True)
        if markers == [3,]:  # PARAM, POST, -2
            self.read_markers([3])
            data = self.read_block()

            self.read_markers([7])
            data = self.read_block()
            #self.show(100)

            data = self._read_record()

            self.read_markers([-1, 0])
        elif markers == [2,]:  # PARAM, POST, -1
            pass
        else:
            raise NotImplementedError(markers)

        #=================
        table_name = self.read_table_name(rewind=True)
        keep_going = True
        table_names = []
        while table_name is not None:
            #print "----------------------------------"
            #print "**", table_name
            table_names.append(table_name)

            if self.debug:
                self.binary_debug.write('-' * 80 + '\n')
                self.binary_debug.write('table_name = %r; f.tell()=%s\n' % (table_name, self.f.tell()))

            #if table_name in self.tables_to_read:
            if 0:
                self.table_name = table_name
                if table_name in ['DIT']:
                    self._read_dit()
                elif table_name in ['PCOMPTS']:
                    self._read_pcompts()
                else:
                    self._skip_table()
            else:
                self.table_name = table_name
                if table_name in ['GEOM1', 'GEOM2', 'GEOM3', 'GEOM4',
                                  'GEOM1N',
                                  'GEOM1OLD',

                                  'EPT',
                                  'MPT', 'MPTS',

                                  'PVT0', 'CASECC',
                                  'EDOM', 'BGPDT', 'OGPFB1',
                                  'DYNAMIC', 'DYNAMICS',
                                  'EQEXIN', 'EQEXINS',
                                  'GPDT', 'ERRORN',
                                  'DESTAB', 'R1TABRG', 'HISADD', 'GPL',

                                   # eigenvalues
                                   'BLAMA', 'LAMA',
                                   # strain energy
                                   'ONRGY1',
                                   # other
                                   'CONTACT', 'VIEWTB', 'OMM2',
                                   # grid point weight
                                   'OGPWG',
                                  ]:
                    self._read_geom_table()  # DIT (agard)
                elif table_name in ['DIT']:
                    self._read_dit()
                elif table_name in ['PCOMPTS']: # blade
                    self._read_pcompts()
                elif table_name in [
                                    # stress
                                    'OES1X1', 'OES1', 'OES1X', 'OES1C', 'OESCP',
                                    'OESNLXR','OESNLXD','OESNLBR','OESTRCP',
                                    'OESNL1X','OESRT',
                                    # strain
                                    'OSTR1X', 'OSTR1C',
                                    # forces
                                    'OEFIT', 'OEF1X', 'OEF1', 'OEF1X',
                                    # spc forces
                                    'OQG1',
                                    # mpc forces
                                    'OQMG1',
                                    # displacement/velocity/acceleration/eigenvector
                                    'OUG1', 'OUGV1', 'BOUGV1',
                                    # applied loads
                                    'OPG1',#'OPG2',

                                    'OUPV1', 'OGS1','OPNL1',
                                    # other
                                    ]:
                    self._read_results_table()
                else:
                    raise NotImplementedError('%r' % table_name)

            table_name = self.read_table_name(rewind=True, stop_on_failure=False)
            #if table_name is None:
                #self.show(100)

        if self.debug:
            self.binary_debug.write('-' * 80 + '\n')
            self.binary_debug.write('f.tell()=%s\ndone...\n' % self.f.tell())

        #self.show_data(data)
        #print "----------------"
        #print "done..."
        return table_names

    def _read_dit(self):
        table_name = self.read_table_name(rewind=False)
        self.read_markers([-1])
        data = self._read_record()

        self.read_markers([-2, 1, 0])
        data = self._read_record()
        table_name, = unpack(b'8s', data)

        self.read_markers([-3, 1, 0])
        data = self._read_record()

        self.read_markers([-4, 1, 0])
        data = self._read_record()

        self.read_markers([-5, 1, 0])
        self.read_markers([0])

    def _read_pcompts(self):
        table_name = self.read_table_name(rewind=False)

        self.read_markers([-1])
        data = self._read_record()

        self.read_markers([-2, 1, 0])
        data = self._read_record()
        table_name, = unpack(b'8s', data)
        #print "table_name = %r" % table_name

        self.read_markers([-3, 1, 0])
        self.read_markers([-4, 1, 0])
        data = self._read_record()

        self.read_markers([-5, 1, 0])
        data = self._read_record()
        self.read_markers([-6, 1, 0])
        self.read_markers([0])
        #self.show(100)

    def read_table_name(self, rewind=False, stop_on_failure=True):
        ni = self.n
        if stop_on_failure:
            data = self._read_record(debug=False)
            table_name, = unpack(b'8s', data)
            if self.debug:
                self.binary_debug.write('marker = [4, 2, 4]\n')
                self.binary_debug.write('table_header = [8, %r, 8]\n\n' % table_name)
            table_name = table_name.strip()
        else:
            try:
                data = self._read_record()
                table_name, = unpack(b'8s', data)
                table_name = table_name.strip()
            except:
                # we're done reading
                self.n = ni
                self.f.seek(self.n)

                try:
                    # we have a trailing 0 marker
                    self.read_markers([0])
                except:
                    # if we hit this block, we have a FATAL error
                    raise FatalError('last table=%r' % self.table_name)
                table_name = None
                rewind = False  # we're done reading, so we're going to ignore the rewind
        #print "table_name1 = %r" % table_name

        if rewind:
            self.n = ni
            self.f.seek(self.n)
        return table_name

    def _skip_table(self):
        if self.debug:
            self.binary_debug.write('skipping table...\n')
        self.table_name = self.read_table_name(rewind=False)
        self.read_markers([-1])
        data = self._skip_record()

        self.read_markers([-2, 1, 0])
        data = self._skip_record()
        self._skip_subtables()

    def _read_geom_table(self):
        self.table_name = self.read_table_name(rewind=False)
        self.read_markers([-1])
        data = self._read_record()
        #self.show_data(data)

        self.read_markers([-2, 1, 0])
        data = self._read_record()
        table_name, = unpack(b'8s', data)
        self._read_subtables()

    def _read_results_table(self):
        if self.debug:
            self.binary_debug.write('read_results_table\n')
        self.table_name = self.read_table_name(rewind=False)
        self.read_markers([-1])
        data = self._read_record()

        self.read_markers([-2, 1, 0])
        data = self._read_record()
        if len(data) == 8:
            subtable_name = unpack(b'8s', data)
            if self.debug:
                self.binary_debug.write('  recordi = [%r]\n'  % subtable_name)
                self.binary_debug.write('  subtable_name=%r\n' % subtable_name)
        elif len(data) == 28:
            subtable_name, month, day, year, zero, one = unpack(b'8s5i', data)
            if self.debug:
                self.binary_debug.write('  recordi = [%r, %i, %i, %i, %i, %i]\n'  % (subtable_name, month, day, year, zero, one))
                self.binary_debug.write('  subtable_name=%r\n' % subtable_name)
            self._print_month(month, day, year, zero, one)
        else:
            raise NotImplementedError(self.show_data(data))

        self._read_subtables()

    def _print_month(self, month, day, year, zero, one):
        self.date = (month, day, 2000 + year)
        #print "%s/%s/%4i" % self.date
        if self.debug:
            self.binary_debug.write('  [subtable_name, month=%i, day=%i, year=%i, zero=%i, one=%i]\n\n' % (month, day, year, zero, one))
        assert zero == 0
        assert one == 1

    def _parse_results_table3(self, data):
        if self.debug:
            self.binary_debug.write('  Table3\n')
        if self.table_name in [# stress
                               'OES1X1', 'OES1', 'OES1X', 'OES1C', 'OESCP',
                               'OESNLXR','OESNLXD','OESNLBR','OESTRCP',
                               'OESNL1X','OESRT',
                               # strain
                               'OSTR1X', 'OSTR1C',]:
            self.read_oes1_3(data)

        elif self.table_name in ['OQG1', 'OQGV1', 'OQP1', 'OQMG1']:
            self.read_oqg1_3(data)
        elif self.table_name in ['OUG1', 'OUGV1']:
            self.read_oug1_3(data)
        elif self.table_name in ['OPG1']:
            self.read_opg1_3(data)
        elif self.table_name in ['OEF1', 'OEFIT', 'OEF1X']:
            self.read_oef1_3(data)
        elif self.table_name in ['OGPWG',]:
            self._read_ogpwg_3(data)
        elif self.table_name in  ['GEOM1', 'GEOM2', 'GEOM3', 'GEOM4',
                                  'GEOM1S',
                                  'GEOM1OLD',

                                  'EPT',
                                  'MPT', 'MPTS',

                                  'PVT0', 'CASECC',
                                  'EDOM', 'BGPDT', 'OGPFB1',
                                  'DYNAMIC', 'DYNAMICS',
                                  'EQEXIN', 'EQEXINS',
                                  'GPDT', 'ERRORN',
                                  'DESTAB', 'R1TABRG', 'HISADD', 'GPL',

                                   # eigenvalues
                                   'BLAMA', 'LAMA',
                                   # strain energy
                                   'ONRGY1',
                                   # other
                                   'CONTACT', 'VIEWTB', 'OMM2',
                                 ]:
            pass
        else:
            raise NotImplementedError(self.table_name)

    def _parse_results_table4(self, data):
        if self.debug:
            self.binary_debug.write('  Table4\n')
        assert len(data) > 0
        if self.table_name in [# stress
                               'OES1X1', 'OES1', 'OES1X', 'OES1C', 'OESCP',
                               'OESNLXR','OESNLXD','OESNLBR','OESTRCP',
                               'OESNL1X','OESRT',
                               # strain
                               'OSTR1X', 'OSTR1C',]:
            self.read_oes1_4(data)

        elif self.table_name in ['OQG1', 'OQGV1', 'OQP1', 'OQMG1']:
            self.read_oqg1_4(data)
        elif self.table_name in ['OUG1', 'OUGV1']:
            self.read_oug1_4(data)
        elif self.table_name in ['OPG1']:
            self.read_opg1_4(data)
        elif self.table_name in ['OEF1', 'OEFIT', 'OEF1X']:
            self.read_oef1_4(data)
        elif self.table_name in ['OGPWG',]:
            self._read_ogpwg_4(data)
        elif self.table_name in  ['GEOM1', 'GEOM2', 'GEOM3', 'GEOM4',
                                  'GEOM1S'
                                  'GEOM1OLD',

                                  'EPT',
                                  'MPT', 'MPTS',

                                  'PVT0', 'CASECC',
                                  'EDOM', 'BGPDT', 'OGPFB1',
                                  'DYNAMIC', 'DYNAMICS',
                                  'EQEXIN', 'EQEXINS',
                                  'GPDT', 'ERRORN',
                                  'DESTAB', 'R1TABRG', 'HISADD', 'GPL',

                                   # eigenvalues
                                   'BLAMA', 'LAMA',
                                   # strain energy
                                   'ONRGY1',
                                   # other
                                   'CONTACT', 'VIEWTB', 'OMM2',
                                 ]:
            pass
        else:
            raise NotImplementedError(self.table_name)

    def finish(self):
        if self.table_name == 'OES1X1':
            self.finish_oes()

    def parse_approach_code(self, data):
        (aCode, tCode, int3, isubcase) = unpack(b'4i', data[:16])
        self.aCode = aCode
        self.tCode = tCode
        self.int3 = int3

        #: the local subcase ID
        self.isubcase = isubcase
        #print("isubcase = %s" %(isubcase))
        #self.subcases.add(self.isubcase)  # set notation

        #: the type of result being processed
        self.table_code = tCode % 1000
        #: used to create sort_bits
        self.sort_code = tCode // 1000
        #: what type of data was saved from the run; used to parse the
        #: approach_code and grid_device.  device_code defines what options
        #: inside a result, STRESS(PLOT,PRINT), are used.
        self.device_code = aCode % 10
        #: what solution was run (e.g. Static/Transient/Modal)
        self.analysis_code = (aCode - self.device_code) // 10

        if self.device_code == 3:
            #sys.stderr.write('The op2 may be inconsistent...\n')
            #sys.stderr.write("  print and plot can cause bad results..."
            #                 "if there's a crash, try plot only\n")
            self.device_code = 1

            #self.log.info('The op2 may be inconsistent...')
            #self.log.info('  print and plot can cause bad results...'
            #              'if there's a crash, try plot only')

        #print('parse_approach_code - aCode=%s tCode=%s int3=%s isubcase=%s' % (aCode, tCode, int3, isubcase))
        #print('                 so - analysis_code=%s device_code=%s table_code=%s sort_code=%s\n' % (self.analysis_code, self.device_code, self.table_code, self.sort_code))
        self._parse_sort_code()

    #===================================
    def _parse_sort_code(self):
        """
        sort_code = 0 -> sort_bits = [0,0,0]
        sort_code = 1 -> sort_bits = [0,0,1]
        sort_code = 2 -> sort_bits = [0,1,0]
        sort_code = 3 -> sort_bits = [0,1,1]
        etc.
        sort_code = 7 -> sort_bits = [1,1,1]

        sort_bits[0] = 0 -> is_sort1=True  isSort2=False
        sort_bits[1] = 0 -> isReal=True   isReal/Imaginary=False
        sort_bits[2] = 0 -> isSorted=True isRandom=False
        """
        bits = [0, 0, 0]
        sort_code = self.sort_code
        i = 2
        while sort_code > 0:
            value = sort_code % 2
            sort_code = (sort_code - value) // 2
            bits[i] = value
            i -= 1
        #: the bytes describe the SORT information
        self.sort_bits = bits
        #print "sort_bits =", bits
        #self.data_code['sort_bits'] = self.sort_bits

    def is_sort1(self):
        if self.sort_bits[0] == 0:
            return True
        return False

    def is_sort2(self):
        return not(self.is_sort1())

    def is_real(self):
        return not(self.is_complex())

    def is_complex(self):
        if self.sort_bits[1] == 1:
            return True
        return False

    def is_random(self):
        if self.sort_bits[1] == 1:
            return True
        return False

    def is_mag_phase(self):
        assert self.format_code in [0, 1], self.format_code
        return bool(self.format_code)

    def is_magnitude_phase(self):
        if self.format_code == 3:
            return True
        return False


if __name__ == '__main__':
    import sys
    o = OP2(sys.argv[1])
    o.read_op2(sys.argv[1])



