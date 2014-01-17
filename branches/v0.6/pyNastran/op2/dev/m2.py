from struct import unpack
from struct import error as StructError
from pyNastran.f06.f06 import FatalError

class OES(object):
    def __init__(self):
        self.oes_mapper = {
            1: 'ROD',
            10: 'CONROD',
            11: 'CELAS1',
            12: 'CELAS2',
            13: 'CELAS3',
            33: 'QUAD4', # 1 node
            34 : 'CBAR',
            39: 'TETRA',
            67: 'HEXA',
            68: 'PENTA',
            74: 'TRIA3', # 1 node
            144: 'QUAD144', # 5 nodes
        }
        pass

    def read_oes1x1_4(self, data):
        #assert self.isubtable == -4, self.isubtable
        self.element_name = self.oes_mapper[self.element_type]

        if self.is_sort1():
            print "element_name =", self.element_name
            n = 0
            if self.element_type in [1, 10]: # CROD
                ntotal = 5 * 4
                nelements = len(data) // ntotal
                for i in xrange(nelements):
                    edata = data[n:n+ntotal]
                    out = unpack('i4f', edata)
                    (eid_device, axial, axial_margin, torsion, torsion_margin) = out
                    eid = (eid_device - self.device_code) // 10
                    print "eid =", eid
                    n += ntotal
            elif self.element_type == 34: # CBAR
                ntotal = 16 * 4
                nelements = len(data) // ntotal
                for i in xrange(nelements):
                    edata = data[n:n+ntotal]
                    out = unpack('i15f', edata)
                    (eid_device, s1a, s2a, s3a, s4a, axial, smaxa, smina, MSt,
                         s1b, s2b, s3b, s4b, smaxb, sminb, MSc) = out
                    eid = (eid_device - self.device_code) // 10
                    n += ntotal
                    print "eid =", eid
            elif self.element_type in [39, 67, 68]: # TETRA, HEXA, PENTA
                if self.element_type == 39: # CTETRA
                    nnodes_expected = 5  # 1 centroid + 4 corner points
                elif self.element_type == 67:
                    nnodes_expected = 9
                elif self.element_type == 68:
                    nnodes_expected = 7
                else:
                    raise NotImplementedError('sort1 Type=%s num=%s' % (self.element_name, self.element_type))

                ntotal = 16 + 84 * nnodes_expected
                #self.show_data(data[:ntotal])
                nelements = len(data) // ntotal
                print "nelements =", nelements
                for i in xrange(nelements):
                    edata = data[n:n+ntotal]
                    eid_device, = unpack('i', data[n:n+4])
                    eid = (eid_device - self.device_code) // 10
                    #out = unpack('ii4si', edata)
                    n += ntotal
                    print "eid =", eid

            #=========================
            elif self.element_type in [33]: # QUAD4-centroidal
                ntotal = 68  # 4*17
                format1 = 'i16f'
                nelements = len(data) // ntotal
                #assert nelements == 2, nelements
                for i in xrange(nelements):
                    edata = data[n:n+ntotal]
                    out = unpack(format1, edata)
                    #if self.make_op2_debug:
                        #self.op2_debug.write('CQUAD4-33A - %s\n' % (str(out)))

                    (eid_device, fd1, sx1, sy1, txy1, angle1, major1, minor1, maxShear1,
                                 fd2, sx2, sy2, txy2, angle2, major2, minor2, maxShear2) = out

                    eid = (eid_device - self.device_code) // 10

                    #print "eid=%i grid=%s fd1=%-3.1f sx1=%i sy1=%i txy1=%i angle1=%i major1=%i minor1=%i vm1=%i" % (eid,'C',fd1,sx1,sy1,txy1,angle1,major1,minor1,maxShear1)
                    #print   "             fd2=%-3.1f sx2=%i sy2=%i txy2=%i angle2=%i major2=%i minor2=%i vm2=%i\n"       % (fd2,sx2,sy2,txy2,angle2,major2,minor2,maxShear2)
                    #print "nNodes = ",nNodes
                    #self.obj.add_new_eid('CQUAD4', dt, eid, 'C', fd1, sx1, sy1,
                    #                   txy1, angle1, major1, minor1, maxShear1)
                    #self.obj.add(dt, eid, 'C', fd2, sx2, sy2, txy2,
                    #             angle2, major2, minor2, maxShear2)
                    #print "eid =", eid
                    n += ntotal

            elif self.element_type in [74]: # TRIA3
                ntotal = 68  # 4*17
                format1 = 'i16f'
                nelements = len(data) // ntotal
                for i in xrange(nelements):
                    edata = data[n:n + ntotal]
                    out = unpack(format1, edata)
                    #if self.make_op2_debug:
                        #self.op2_debug.write('CTRIA3-74 - %s\n' % (str(out)))

                    (eid_device, fd1, sx1, sy1, txy1, angle1, major1, minor1, vm1,
                                 fd2, sx2, sy2, txy2, angle2, major2, minor2, vm2,) = out
                    eid = (eid_device - self.device_code) // 10
                    print "eid =", eid
                    n += ntotal
            elif self.element_type in [144]: # CQUAD4-bilinear
                nnodes_expected = 4 # + 1 centroid
                return
                ntotal = 16 + 84 * nnodes_expected
                assert ntotal == 348, ntotal
                nelements = len(data) // ntotal
                for i in xrange(nelements):
                    edata = data[n:n+ntotal]
                    eid_device, = unpack('i', data[:4])
                    eid = (eid_device - self.device_code) // 10
                    #out = unpack('ii4si', edata)
                    n += ntotal
                    print "eid =", eid
            else:
                raise NotImplementedError('sort1 Type=%s num=%s' % (self.element_name, self.element_type))
        else:
            raise NotImplementedError('sort2 Type=%s num=%s' % (self.element_name, self.element_type))
        assert len(data) % ntotal == 0, '%s n=%s len=%s ntotal=%s' % (self.element_name, len(data) % ntotal, len(data), ntotal)
        print "element_name =", self.element_name
        print "----------------"
        #sys.exit('stopping...')

    def apply_data_code_value(self, name, value):
        pass
    def setNullNonlinearFactor(self):
        pass

    def is_sort1(self):
        if self.sort_bits[0] == 0:
            return True
        return False

    def is_sort2(self):
        return not(self.is_sort1())

    def read_oes1x1_3(self, data):
        self.words = ['aCode',       'tCode',    'element_type', 'isubcase',
                 '???',         '???',      '???',          'load_set'
                 'format_code', 'num_wide', 's_code',       '???',
                 '???',         '???',      '???',          '???',
                 '???',         '???',      '???',          '???',
                 '???',         '???',      '???',          '???',
                 '???', 'Title', 'subtitle', 'label']

        #self.show_data(data)
        self.parse_approach_code(data)  # 3
        ## element type
        self.add_data_parameter(data, 'element_type', 'i', 3, False)
        ## load set ID
        self.add_data_parameter(data, 'load_set', 'i', 8, False)
        ## format code
        self.add_data_parameter(data, 'format_code', 'i', 9, False)
        ## number of words per entry in record
        ## .. note:: is this needed for this table ???
        self.add_data_parameter(data, 'num_wide', 'i', 10, False)
        ## stress/strain codes
        self.add_data_parameter(data, 's_code', 'i', 11, False)
        ## thermal flag; 1 for heat ransfer, 0 otherwise
        self.add_data_parameter(data, 'thermal', 'i', 23, False)

        ## assuming tCode=1
        if self.analysis_code == 1:   # statics / displacement / heat flux
            ## load set number
            self.add_data_parameter(data, 'lsdvmn', 'i', 5, False)
            self.apply_data_code_value('dataNames', ['lsdvmn'])
            self.setNullNonlinearFactor()
        elif self.analysis_code == 2:  # real eigenvalues
            ## mode number
            self.add_data_parameter(data, 'mode', 'i', 5)
            ## real eigenvalue            self.add_data_parameter(data, 'eign', 'f', 6, False)
            ## mode or cycle TODO confused on the type - F1???
            self.add_data_parameter(data, 'mode_cycle', 'f', 7, False)
            self.apply_data_code_value('dataNames', ['mode', 'eigr', 'mode_cycle'])
        #elif self.analysis_code==3: # differential stiffness
            #self.lsdvmn = self.get_values(data,'i',5) ## load set number
            #self.data_code['lsdvmn'] = self.lsdvmn
        #elif self.analysis_code==4: # differential stiffness
        #    self.lsdvmn = self.get_values(data,'i',5) ## load set number
        elif self.analysis_code == 5:   # frequency
            ## frequency
            self.add_data_parameter(data, 'freq', 'f', 5)
            self.apply_data_code_value('dataNames', ['freq'])
        elif self.analysis_code == 6:  # transient
            ## time step
            self.add_data_parameter(data, 'dt', 'f', 5)
            self.apply_data_code_value('dataNames', ['dt'])
        elif self.analysis_code == 7:  # pre-buckling
            ## load set
            self.add_data_parameter(data, 'lsdvmn', 'i', 5)
            self.apply_data_code_value('dataNames', ['lsdvmn'])
        elif self.analysis_code == 8:  # post-buckling
            ## mode number
            self.add_data_parameter(data, 'lsdvmn', 'i', 5)
            self.add_data_parameter(data, 'eigr', 'f', 6, False)  # real eigenvalue
            self.apply_data_code_value('dataNames', ['lsdvmn', 'eigr'])
        elif self.analysis_code == 9:  # complex eigenvalues
            ## mode number
            self.add_data_parameter(data, 'mode', 'i', 5)
            ## real eigenvalue
            self.add_data_parameter(data, 'eigr', 'f', 6, False)
            ## imaginary eigenvalue            self.add_data_parameter(data, 'eigi', 'f', 7, False)
            self.apply_data_code_value('dataNames', ['mode', 'eigr', 'eigi'])
        elif self.analysis_code == 10:  # nonlinear statics
            ## load step
            self.add_data_parameter(data, 'lftsfq', 'f', 5)
            self.apply_data_code_value('dataNames', ['lftsfq'])
        elif self.analysis_code == 11:  # old geometric nonlinear statics
            ## load set number
            self.add_data_parameter(data, 'lsdvmn', 'i', 5)
            self.apply_data_code_value('dataNames', ['lsdvmn'])
        elif self.analysis_code == 12:  # contran ? (may appear as aCode=6)  --> straight from DMAP...grrr...
            ## Time step ??? --> straight from DMAP
            self.add_data_parameter(data, 'dt', 'f', 5)
            self.apply_data_code_value('dataNames', ['dt'])
        else:
            raise RuntimeError('invalid analysis_code...analysis_code=%s' %
                               self.analysis_code)
        # tCode=2
        #if self.analysis_code==2: # sort2
        #    self.lsdvmn = self.get_values(data,'i',5)

        if not self.is_sort1() and not isRelease:
            raise NotImplementedError('sort2...')

        if self.debug:
            self.binary_debug.write('  element_name = %r\n' % self.oes_mapper[self.element_type])
            self.binary_debug.write('  aCode    = %r\n' % self.aCode)
            self.binary_debug.write('  tCode    = %r\n' % self.tCode)
            self.binary_debug.write('  isubcase = %r\n' % self.isubcase)

        #self.show_data(data[200:])
        self.read_title(data)
        if self.element_type not in self.oes_mapper:
            raise NotImplementedError(self.element_type)

        if self.debug:
            msg = ''
            for i, param in enumerate(self.words):
                if param == '???':
                    param = 0
                msg += '%s, ' % param
                if i % 5 == 4:
                    msg += '\n             '

            self.binary_debug.write('  recordi = [%s]\n\n' % msg)
        #sys.exit()

    def read_title(self, data):
        assert len(data) == 584, len(data)
        Title, subtitle, label = unpack('128s128s128s', data[200:])  # titleSubtitleLabel

        self.Title = Title.strip()
        #: the subtitle of the subcase
        self.subtitle = subtitle.strip()
        #: the label of the subcase
        self.label = label.strip()
        if self.debug:
            self.binary_debug.write('  Title    = %r\n' % self.Title)
            self.binary_debug.write('  subtitle = %r\n' % self.subtitle)
            self.binary_debug.write('  label    = %r\n' % self.label)

    def finish_oes(self):
        del self.element_type
        del self.load_set
        del self.format_code
        del self.num_wide
        del self.s_code
        del self.thermal

class OP2(OES):
    def __init__(self):
        #self.tables_to_read = []
        OES.__init__(self)
        self.debug = True
        self.binary_debug = open('debug.out', 'wb')

    def show(self, n):
        assert self.n == self.f.tell()
        nints = n//4
        data = self.f.read(n)
        self.show_data(data)
        self.f.seek(self.n)

    def show_data(self, data):
        n = len(data)
        nints = n // 4
        strings = unpack('%is' % n, data)
        ints    = unpack('%ii' % nints, data)
        floats  = unpack('%if' % nints, data)
        print "strings =", strings
        print "ints    =", ints
        #print "floats  =", floats

    def skip_block(self):
        data = self.f.read(4)
        ndata, = unpack('i', data)
        self.n += 8 + ndata
        self.goto(self.n)
        return None

    def read_block(self):
        data = self.f.read(4)
        ndata, = unpack('i', data)

        data_out = self.f.read(ndata)
        data = self.f.read(4)
        self.n += 8 + ndata
        return data_out

    def read_markers(self, markers):
        for i, marker in enumerate(markers):
            data = self.read_block()
            imarker, = unpack('i', data)
            assert marker == imarker, 'marker=%r imarker=%r; markers=%s i=%s' % (marker, imarker, markers, i)

    def get_nmarkers(self, n, rewind=True):
        ni = self.n
        markers = []
        for i in xrange(n):
            data = self.read_block()
            marker, = unpack('i', data)
            #print "marker=%s" % marker
            markers.append(marker)
        if rewind:
            self.n = ni
            self.f.seek(self.n)
        return markers

    def read_op2(self, op2_filename):
        self.op2_filename = op2_filename
        self.n = 0
        self.table_name = None
        self.f = open(op2_filename, 'rb')

        markers = self.get_nmarkers(1, rewind=True)
        #print "markers =", markers
        if markers == [3,]:
            self.read_markers([3])
            data = self.read_block()

            self.read_markers([7])
            data = self.read_block()
            #self.show(100)

            data = self.read_record()

            self.read_markers([-1, 0])

        #=================
        table_name = self.read_table_name(rewind=True)
        keep_going = True
        table_names = []
        while table_name is not None:
            print "----------------------------------"
            print "**", table_name
            table_names.append(table_name)

            if self.debug:
                self.binary_debug.write('-' * 80 + '\n')
                self.binary_debug.write('table_name = %r; f.tell()=%s\n' % (table_name, self.f.tell()))

            #if table_name in self.tables_to_read:
            if 0:
                self.table_name = table_name
                self._skip_table()
            else:
                self.table_name = table_name
                if table_name in ['GEOM1', 'GEOM2', 'GEOM3', 'GEOM4',
                                  'EPT', 'MPT', 'MPTS', 'PVT0', 'CASECC',
                                  'EDOM', 'BGPDT', 'EQEXINS', 'OGPWG', 'OGPFB1',
                                  'DYNAMICS', 'DESTAB', 'R1TABRG', 'HISADD', 'GPL',
                                  'GPDT']:
                    self.read_geom_table()  # DIT (agard)
                elif table_name in ['DIT']:
                    self.read_dit()
                elif table_name in ['PCOMPTS']: # blade
                    self.read_pcompts()
                elif table_name in ['OQMG1', 'OSTR1X', 'OESNLXR', 'OQG1', 'OUGV1',
                                    'OEF1X', 'OES1X1', 'OPG1', 'OES1C', 'OSTR1C',
                                    'OEFIT', 'BOUGV1', 'OES1', 'OEF1', ]:
                    self.read_results_table()
                else:
                    raise NotImplementedError('%r' % table_name)

            table_name = self.read_table_name(rewind=True, stop_on_failure=False)
            if table_name is None:
                self.show(100)

        if self.debug:
            self.binary_debug.write('-' * 80 + '\n')
            self.binary_debug.write('f.tell()=%s\ndone...\n' % self.f.tell())

        #self.show_data(data)
        print "----------------"
        print "done..."
        return table_names

    def read_dit(self):
        table_name = self.read_table_name(rewind=False)
        self.read_markers([-1])
        data = self.read_record()

        self.read_markers([-2, 1, 0])
        data = self.read_record()
        table_name, = unpack('8s', data)
        #print "table_name = %r" % table_name

        self.read_markers([-3, 1, 0])
        data = self.read_record()

        self.read_markers([-4, 1, 0])
        data = self.read_record()

        self.read_markers([-5, 1, 0])
        self.read_markers([0])

        #self.show(100)
        #sys.exit()

    def read_pcompts(self):
        table_name = self.read_table_name(rewind=False)

        self.read_markers([-1])
        data = self.read_record()

        self.read_markers([-2, 1, 0])
        data = self.read_record()
        table_name, = unpack('8s', data)
        #print "table_name = %r" % table_name

        self.read_markers([-3, 1, 0])
        self.read_markers([-4, 1, 0])
        data = self.read_record()

        self.read_markers([-5, 1, 0])
        data = self.read_record()
        self.read_markers([-6, 1, 0])
        self.read_markers([0])
        #self.show(100)

    def read_table_name(self, rewind=False, stop_on_failure=True):
        ni = self.n
        if stop_on_failure:
            data = self.read_record(debug=False)
            table_name, = unpack('8s', data)
            if self.debug:
                self.binary_debug.write('marker = [4, 2, 4]\n')
                self.binary_debug.write('table_header = [8, %r, 8]\n\n' % table_name)
            table_name = table_name.strip()
        else:
            try:
                data = self.read_record()
                table_name, = unpack('8s', data)
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
        data = self.skip_record()

        self.read_markers([-2, 1, 0])
        data = self.skip_record()
        self.skip_subtables()

    def read_geom_table(self):
        self.table_name = self.read_table_name(rewind=False)
        self.read_markers([-1])
        data = self.read_record()
        #self.show_data(data)

        self.read_markers([-2, 1, 0])
        data = self.read_record()
        table_name, = unpack('8s', data)
        self.read_subtables()

    def read_results_table(self):
        if self.debug:
            self.binary_debug.write('read_results_table\n')
        self.table_name = self.read_table_name(rewind=False)
        self.read_markers([-1])
        data = self.read_record()

        self.read_markers([-2, 1, 0])
        data = self.read_record()

        subtable_name, month, day, year, zero, one = unpack('8s5i', data)

        if self.debug:
            self.binary_debug.write('  recordi = [%r, %i, %i, %i, %i, %i]\n'  % (subtable_name, month, day, year, zero, one))
            self.binary_debug.write('  subtable_name=%r\n' % subtable_name)
        self.print_month(month, day, year, zero, one)
        self.read_subtables()
        #if self.table_name == 'OES1X1':
            #asdfa

    def print_month(self, month, day, year, zero, one):
        self.date = (month, day, 2000 + year)
        #print "%s/%s/%4i" % self.date
        if self.debug:
            self.binary_debug.write('  [subtable_name, month=%i, day=%i, year=%i, zero=%i, one=%i]\n\n' % (month, day, year, zero, one))
        assert zero == 0
        assert one == 1

    def parse_results_table3(self, data):
        if self.table_name in ['OES1X1', 'OES1']:
            self.read_oes1x1_3(data)

    def parse_results_table4(self, data):
        if self.table_name in ['OES1X1', 'OES1']:
            self.read_oes1x1_4(data)

    def skip_subtables(self):
        self.isubtable = -3
        self.read_markers([-3, 1, 0])

        markers = self.get_nmarkers(1, rewind=True)
        while markers[0] != 0:
            data = self.skip_record()
            #if len(data) == 584:
                #self.parse_results_table3(data)
            #else:
                #data = self.parse_results_table4(data)

            self.isubtable -= 1
            self.read_markers([self.isubtable, 1, 0])
            markers = self.get_nmarkers(1, rewind=True)
        self.read_markers([0])

    def read_subtables(self):
        self.isubtable = -3
        self.read_markers([-3, 1, 0])
        #data = self.read_record()
        #self.parse_results_table3(data)

        #self.isubtable -= 1 # -4
        #self.read_markers([self.isubtable, 1, 0])
        markers = self.get_nmarkers(1, rewind=True)
        while markers[0] != 0:
            data = self.read_record()
            if len(data) == 584:
                self.parse_results_table3(data)
            else:
                data = self.parse_results_table4(data)

            self.isubtable -= 1
            self.read_markers([self.isubtable, 1, 0])
            markers = self.get_nmarkers(1, rewind=True)
        self.read_markers([0])
        self.finish()

    def finish(self):
        if self.table_name == 'OES1X1':
            self.finish_oes()

    def read_ougv1(self):
        pass

    def parse_approach_code(self, data):
        (aCode, tCode, int3, isubcase) = unpack(b'iiii', data[:16])
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

    def add_data_parameter(self, data, var_name, Type, slot, debug=False):
        datai = data[4*(slot-1) : 4*(slot)]
        assert len(datai) == 4, len(datai)
        value, = unpack(Type, datai)
        #print "%-12s = %r" % (var_name, value)
        if self.debug:
            self.binary_debug.write('  %-12s = %r\n' % (var_name, value))
        setattr(self, var_name, value)
        self.words[slot-1] = var_name

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
        print "sort_bits =", bits

        #self.data_code['sort_bits'] = self.sort_bits

    def goto(self, n):
        self.n = n
        self.f.seek(n)

    def skip_record(self):
        markers0 = self.get_nmarkers(1, rewind=False)

        #record = self.read_block()
        record = self.skip_block()

        markers1 = self.get_nmarkers(1, rewind=True)
        # handling continuation blocks
        if markers1[0] > 0:
            nloop = 0
            while markers1[0] > 0: #, 'markers0=%s markers1=%s' % (markers0, markers1)
                markers1 = self.get_nmarkers(1, rewind=False)
                record = self.read_block()
                markers1 = self.get_nmarkers(1, rewind=True)

                nloop += 1
            if nloop > 0:
                print "nloop = %s" % nloop

        return record

    def read_record(self, debug=True):
        markers0 = self.get_nmarkers(1, rewind=False)
        if self.debug and debug:
            self.binary_debug.write('marker = [4, %i, 4]\n' % markers0[0])
        #print "markers =", markers0
        record = self.read_block()
        if self.debug and debug:
            nrecord = len(record)
            self.binary_debug.write('record = [%i, recordi, %i]\n' % (nrecord, nrecord))
        #if self.table_name == 'OES1X1':
            #self.show_data(record)
        assert (markers0[0]*4) == len(record), 'markers0=%s*4 len(record)=%s' % (markers0[0]*4, len(record))

        markers1 = self.get_nmarkers(1, rewind=True)
        #print "markers1 =", markers1

        # handling continuation blocks
        if markers1[0] > 0:
            nloop = 0
            records = [record]
            while markers1[0] > 0: #, 'markers0=%s markers1=%s' % (markers0, markers1)
                markers1 = self.get_nmarkers(1, rewind=False)
                record = self.read_block()
                markers1 = self.get_nmarkers(1, rewind=True)

                records.append(record)
                nloop += 1
            if nloop > 0:
                print "nloop = %s" % nloop
                record = ''.join(records)
        return record


if __name__ == '__main__':
    o = OP2()
    import sys
    print o.read_op2(sys.argv[1])



