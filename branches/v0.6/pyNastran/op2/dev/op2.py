from struct import unpack
from struct import error as StructError

from pyNastran.f06.f06 import FatalError
from pyNastran.f06.tables.grid_point_weight import GridPointWeight

from pyNastran.op2.dev.oef import OEF
from pyNastran.op2.dev.oes import OES
from pyNastran.op2.dev.ogs import OGS

from pyNastran.op2.dev.opg import OPG
from pyNastran.op2.dev.oqg import OQG
from pyNastran.op2.dev.oug import OUG
from pyNastran.op2.dev.ogpwg import OGPWG
from pyNastran.op2.dev.fortran_format import FortranFormat

class OP2(OEF, OES, OGS, OPG, OQG, OUG, OGPWG, FortranFormat):
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
        OGS.__init__(self)

        OPG.__init__(self)
        OQG.__init__(self)
        OUG.__init__(self)
        OGPWG.__init__(self)
        FortranFormat.__init__(self)

        self.grid_point_weight = GridPointWeight()
        self.words = []
        self.debug = True
        if self.debug:
            self.binary_debug = open('debug.out', 'wb')

        self._table_mapper = {
            #'DIT': self.read_dit,
            #=======================
            # OEF
            # internal forces
            'OEFIT' : [self.read_oef1_3, self.read_oef1_4],
            'OEF1X' : [self.read_oef1_3, self.read_oef1_4],
            'OEF1'  : [self.read_oef1_3, self.read_oef1_4],
            'DOEF1' : [self.read_oef1_3, self.read_oef1_4],
            #=======================
            # OQG
            # spc forces
            'OQG1'  : [self.read_oqg1_3, self.read_oqg1_4],
            'OQGV1' : [self.read_oqg1_3, self.read_oqg1_4],
            # mpc forces
            'OQMG1' : [self.read_oqg1_3, self.read_oqg1_4],
            #=======================
            # OPG
            # applied loads
            'OPG1'  : [self.read_opg1_3, self.read_opg1_4],
            'OPGV1' : [self.read_opg1_3, self.read_opg1_4],
            'OPNL1' : [self.read_opg1_3, self.read_opg1_4],

            #=======================
            # OES
            # stress
            'OES1X1'  : [self.read_oes1_3, self.read_oes1_4],
            'OES1'    : [self.read_oes1_3, self.read_oes1_4],
            'OES1X'   : [self.read_oes1_3, self.read_oes1_4],
            'OES1C'   : [self.read_oes1_3, self.read_oes1_4],
            'OESCP'   : [self.read_oes1_3, self.read_oes1_4],
            'OESNLXR' : [self.read_oes1_3, self.read_oes1_4],
            'OESNLXD' : [self.read_oes1_3, self.read_oes1_4],
            'OESNLBR' : [self.read_oes1_3, self.read_oes1_4],
            'OESTRCP' : [self.read_oes1_3, self.read_oes1_4],
            'OESNL1X' : [self.read_oes1_3, self.read_oes1_4],
            'OESRT'   : [self.read_oes1_3, self.read_oes1_4],
            # strain
            'OSTR1X'  : [self.read_oes1_3, self.read_oes1_4],
            'OSTR1C'  : [self.read_oes1_3, self.read_oes1_4],

            #=======================
            # OUG
            # displacement/velocity/acceleration/eigenvector
            'OUG1'   : [self.read_oug1_3, self.read_oug1_4],
            'OUGV1'  : [self.read_oug1_3, self.read_oug1_4],
            'BOUGV1' : [self.read_oug1_3, self.read_oug1_4],
            'OUPV1'  : [self.read_oug1_3, self.read_oug1_4],

            #=======================
            # OGPWG
            # grid point weight
            'OGPWG' : [self._read_ogpwg_3, self._read_ogpwg_4],

            #=======================
            # OGS
            # grid point stresses
            'OGS1' : [self._read_ogs1_3, self._read_ogs1_4],
            #=======================
            # OEE
            #'ONRGY1' : [self._read_oee1_3, self._read_oee1_4],
        }

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
        table_name = self.read_table_name(rewind=True, stop_on_failure=False)
        if table_name is None:
            raise FatalError('no tables exists...')

        keep_going = True
        table_names = []
        while table_name is not None:
            #print "----------------------------------"
            table_names.append(table_name)

            if self.debug:
                print "**", table_name
                self.binary_debug.write('-' * 80 + '\n')
                self.binary_debug.write('table_name = %r; f.tell()=%s\n' % (table_name, self.f.tell()))

            self.table_name = table_name
            #if table_name in self.tables_to_read:
            if 0:
                self._skip_table(table_name)
            else:
                if table_name in ['GEOM1', 'GEOM2', 'GEOM3', 'GEOM4',  # regular
                                  'GEOM1S', 'GEOM2S', 'GEOM3S', 'GEOM4S', # superelements
                                  'GEOM1N',
                                  'GEOM1OLD', 'GEOM2OLD', 'GEOM4OLD',

                                  'EPT', 'EPTS', 'EPTOLD',
                                  'MPT', 'MPTS',

                                  'PVT0', 'CASECC',
                                  'EDOM', 'OGPFB1',
                                  'BGPDT', 'BGPDTOLD',
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
                elif table_name in ['DIT']:  # tables
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
                                    'OEFIT', 'OEF1X', 'OEF1', 'DOEF1',
                                    # spc forces
                                    'OQG1', 'OQGV1',
                                    # mpc forces
                                    'OQMG1',
                                    # ??? forces
                                    'OQP1',
                                    # displacement/velocity/acceleration/eigenvector
                                    'OUG1', 'OUGV1', 'BOUGV1', 'OUPV1',
                                    # applied loads
                                    'OPG1', 'OPGV1', 'OPNL1', #'OPG2',

                                    # grid point stresses
                                    'OGS1',

                                    # other
                                    'OPNL1','OFMPF2M',
                                    'OSMPF2M','OPMPF2M','OLMPF2M', 'OGPMPF2M',

                                    'OAGPSD2', 'OAGCRM2', 'OAGRMS2', 'OAGATO2', 'OAGNO2',
                                    'OESPSD2', 'OESCRM2', 'OESRMS2', 'OESATO2', 'OESNO2',
                                    'OEFPSD2', 'OEFCRM2', 'OEFRMS2', 'OEFATO2', 'OEFNO2',
                                    'OPGPSD2', 'OPGCRM2', 'OPGRMS2', 'OPGATO2', 'OPGNO2', 
                                    'OQGPSD2', 'OQGCRM2', 'OQGRMS2', 'OQGATO2', 'OQGNO2', 
                                    'OQMPSD2', 'OQMCRM2', 'OQMRMS2', 'OQMATO2', 'OQMNO2',
                                    'OSTRPSD2','OSTRCRM2','OSTRRMS2','OSTRATO2','OSTRNO2',
                                    'OUGPSD2', 'OUGCRM2', 'OUGRMS2', 'OUGATO2', 'OUGNO2',
                                    'OVGPSD2', 'OVGCRM2', 'OVGRMS2', 'OVGATO2', 'OVGNO2',
                                    'OCRUG',
                                    'OCRPG', 
                                    'STDISP', 'FOL',
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

    def _skip_table(table_name):
        if table_name in ['DIT']:  # tables
            self._read_dit()
        elif table_name in ['PCOMPTS']:
            self._read_pcompts()
        else:
            self._skip_table_helper()

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

        markers = self.get_nmarkers(1, rewind=True)
        if markers != [0]:
            data = self._read_record()
            self.read_markers([-6, 1, 0])

        markers = self.get_nmarkers(1, rewind=True)
        if markers != [0]:
            data = self._read_record()
            self.read_markers([-7, 1, 0])

        markers = self.get_nmarkers(1, rewind=True)
        if markers != [0]:
            data = self._read_record()
            self.read_markers([-8, 1, 0])

        markers = self.get_nmarkers(1, rewind=True)
        if markers != [0]:
            data = self._read_record()
            self.read_markers([-9, 1, 0])

        #self.show(100)
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
        markers = self.get_nmarkers(1, rewind=True)
        if markers != [-4]:
            data = self._read_record()

        self.read_markers([-4, 1, 0])
        markers = self.get_nmarkers(1, rewind=True)
        if markers != [0]:
            data = self._read_record()
        else:
            self.read_markers([0])
            return
            

        self.read_markers([-5, 1, 0])
        data = self._read_record()

        self.read_markers([-6, 1, 0])
        self.read_markers([0])

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

    def _skip_table_helper(self):
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

    def finish(self):
        for word in self.words:
            if word != '???' and hasattr(self, word):
                delattr(self, word)


if __name__ == '__main__':
    import sys
    o = OP2(sys.argv[1])
    o.read_op2(sys.argv[1])



