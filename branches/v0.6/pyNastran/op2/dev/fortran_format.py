from struct import unpack

class FortranFormat(object):
    def __init__(self):
        self.n = 0
        self.f = None
        
    def show(self, n):
        assert self.n == self.f.tell()
        nints = n // 4
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
        print "ints    =", ints, '\n'
        print "floats  =", floats

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
            markers.append(marker)
        if rewind:
            self.n = ni
            self.f.seek(self.n)
        return markers

    def _skip_subtables(self):
        self.isubtable = -3
        self.read_markers([-3, 1, 0])

        markers = self.get_nmarkers(1, rewind=True)
        while markers[0] != 0:
            data = self._skip_record()
            #if len(data) == 584:
                #self._parse_results_table3(data)
            #else:
                #data = self._parse_results_table4(data)

            self.isubtable -= 1
            self.read_markers([self.isubtable, 1, 0])
            markers = self.get_nmarkers(1, rewind=True)
        self.read_markers([0])

    def _read_subtables(self):
        self.isubtable = -3
        self.read_markers([-3, 1, 0])
        #data = self._read_record()
        #self._parse_results_table3(data)

        #self.isubtable -= 1 # -4
        #self.read_markers([self.isubtable, 1, 0])
        markers = self.get_nmarkers(1, rewind=True)
        while markers[0] != 0:
            data = self._read_record()
            if len(data) == 584:
                self._parse_results_table3(data)
            else:
                data = self._parse_results_table4(data)

            self.isubtable -= 1
            self.read_markers([self.isubtable, 1, 0])
            markers = self.get_nmarkers(1, rewind=True)
        self.read_markers([0])
        self.finish()

    def goto(self, n):
        self.n = n
        self.f.seek(n)

    def _skip_record(self):
        markers0 = self.get_nmarkers(1, rewind=False)

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
        return record

    def _read_record(self, debug=True):
        markers0 = self.get_nmarkers(1, rewind=False)
        if self.debug and debug:
            self.binary_debug.write('marker = [4, %i, 4]\n' % markers0[0])
        record = self.read_block()
        if self.debug and debug:
            nrecord = len(record)
            self.binary_debug.write('record = [%i, recordi, %i]\n' % (nrecord, nrecord))
        assert (markers0[0]*4) == len(record), 'markers0=%s*4 len(record)=%s' % (markers0[0]*4, len(record))

        markers1 = self.get_nmarkers(1, rewind=True)

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
                record = ''.join(records)
        return record