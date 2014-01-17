from struct import unpack
from struct import error as StructError
from pyNastran.f06.f06 import FatalError
from pyNastran.op2.dev.oes import OES

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
                                  'GPDT', 'LAMA', 'EQEXIN']:
                    self.read_geom_table()  # DIT (agard)
                elif table_name in ['DIT']:
                    self.read_dit()
                elif table_name in ['PCOMPTS']: # blade
                    self.read_pcompts()
                elif table_name in ['OQMG1', 'OSTR1X', 'OESNLXR', 'OQG1', 'OUGV1',
                                    'OEF1X', 'OES1X1', 'OPG1', 'OES1C', 'OSTR1C',
                                    'OEFIT', 'BOUGV1', 'OES1', 'OEF1', 'ONRGY1']:
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



