from struct import unpack
import copy

class FortranFormat(object):
    def __init__(self):
        self.n = 0
        self.f = None
        self.obj = None
        self.data_code = None
        self.table_name = None
        self.isubcase = None
        self._table_mapper = {}

        #: stores if the user entered [] for iSubcases
        self.isAllSubcases = True
        self.valid_subcases = []

    def show(self, n):
        assert self.n == self.f.tell()
        nints = n // 4
        data = self.f.read(n)
        strings, ints, floats = self.show_data(data)
        self.f.seek(self.n)
        return strings, ints, floats

    def show_data(self, data):
        n = len(data)
        nints = n // 4
        strings = unpack(b'%is' % n, data)
        ints    = unpack(b'%ii' % nints, data)
        floats  = unpack(b'%if' % nints, data)
        print "strings =", strings
        print "ints    =", ints, '\n'
        print "floats  =", floats
        return strings, ints, floats

    def skip_block(self):
        """
        Skips a block following a pattern of:
            [nbytes, data, nbytes]
        :retval data: since data can never be None, a None value
                      indicates something bad happened.
        """
        data = self.f.read(4)
        ndata, = unpack(b'i', data)
        self.n += 8 + ndata
        self.goto(self.n)
        return None

    def read_block(self):
        """
        Reads a block following a pattern of:
            [nbytes, data, nbytes]
        :retval data: the data in binary
        """
        data = self.f.read(4)
        ndata, = unpack(b'i', data)

        data_out = self.f.read(ndata)
        data = self.f.read(4)
        self.n += 8 + ndata
        return data_out

    def read_markers(self, markers):
        """
        Gets specified markers, where a marker has the form of [4, value, 4].
        The "marker" corresponds to the value, so 3 markers takes up 9 integers.
        These are used to indicate position in the file as well as
        the number of bytes to read.

        :param markers: markers to get; markers = [-10, 1]
        """
        for i, marker in enumerate(markers):
            data = self.read_block()
            imarker, = unpack(b'i', data)
            assert marker == imarker, 'marker=%r imarker=%r; markers=%s i=%s' % (marker, imarker, markers, i)

    def get_nmarkers(self, n, rewind=True):
        """
        Gets n markers, so if n=2, it will get 2 markers.

        :param n: number of markers to get
        :param rewind: should the file be returned to the starting point
        :retval markers: list of [1, 2, 3, ...] markers
        """
        ni = self.n
        markers = []
        for i in xrange(n):
            data = self.read_block()
            marker, = unpack(b'i', data)
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
            self.log.debug("skipping table_name = %r" % self.table_name)
            #if len(data) == 584:
                #self._parse_results_table3(data)
            #else:
                #data = self._parse_results_table4(data)

            self.isubtable -= 1
            self.read_markers([self.isubtable, 1, 0])
            markers = self.get_nmarkers(1, rewind=True)
        self.read_markers([0])

    def passer(self, data):
        """
        dummy function used for unsupported tables
        """
        pass

    def _get_table_mapper(self):
        raise NotImplementedError('this should be overwritten')

    def _read_subtables(self):
        self._data_factor = 1
        nstart = self.n
        self.isubtable = -3
        self.read_markers([-3, 1, 0])
        #data = self._read_record()
        #self._parse_results_table3(data)

        #self.isubtable -= 1 # -4
        #self.read_markers([self.isubtable, 1, 0])
        table_mapper = self._get_table_mapper()

        if self.table_name in table_mapper:
            if self.read_mode in [0, 2]:
                self.log.debug("table_name = %r" % self.table_name)
            table3_parser, table4_parser = table_mapper[self.table_name]
            passer = False
        else:
            #raise NotImplementedError(self.table_name)
            if self.read_mode in [0, 2]:
                self.log.debug("skipping table_name = %r" % self.table_name)
            table3_parser = None
            table4_parser = None
            passer = True

        markers = self.get_nmarkers(1, rewind=True)
        while markers[0] != 0:
            record_len = self._get_record_length()
            if record_len == 584:
                self.data_code = {}  #'log': self.log,}  # resets the logger
                self.obj = None
                data = self._read_record()
                if not passer:
                    table3_parser(data)
                    #if hasattr(self, 'isubcase'):
                        #print "code = ", self._get_code()

            else:
                if passer or not self.is_valid_subcase():
                    data = self._skip_record()
                else:
                    if hasattr(self, 'num_wide'):
                        datai = b''
                        if self.read_mode in [0, 2]:
                            self.ntotal = 0
                            #print "read_mode", self.read_mode
                            #print "-----------------------------------------"
                            for data in self._stream_record():
                                data = datai + data
                                n = table4_parser(data)
                                assert isinstance(n, int), self.table_name
                                datai = data[n:]

                            if hasattr(self, 'eid_old'):
                                del self.eid_old

                            if self.read_mode == 2:
                                if hasattr(self, 'obj') and hasattr(self.obj, 'itime'):
                                    #ntotal = record_len // (self.num_wide * 4) * self._data_factor
                                    if self.obj.ntotal == self.obj.data.shape[1]:
                                        self.obj._reset_indices()
                                        self.obj.itime += 1
                                        #self.obj.itotal = 0
                                    else:
                                        print('self.obj.name=%r has itime', self.obj.__class__.__name__)
                                        print('ntotal=%s shape=%s' % (self.obj.ntotal, str(self.obj.data.shape)))
                                #else:
                                    #print('self.obj.name=%r doesnt have itime' % self.obj.__class__.__name__)


                        elif self.read_mode == 1:
                            #n = self._skip_record()
                            #n = table4_parser(datai, 300000)
                            if 1:
                                self.ntotal = 0
                                #n = self.n
                                n = 0
                                for i, data in enumerate(self._stream_record()):
                                    #if i == 0:
                                        data = datai + data
                                        ndata = len(data)
                                        n = table4_parser(data)
                                        assert isinstance(n, int), self.table_name
                                        datai = data[n:]
                                        #if self.obj is not None:
                                            #print "len(datai) =", len(datai), ndata / (self.num_wide * 4.)
                                assert len(datai) == 0, len(datai)
                                #n = record_len
                                #break
                            #self.goto(n)
                            #n = self._skip_record()

                            if hasattr(self, 'obj') and self.obj is not None:
                                if hasattr(self.obj, 'ntimes'):
                                    if not hasattr(self.obj, '_reset_indices'):
                                        methods = '\ndir(obj)=%s' % ', '.join(sorted(dir(self.obj)))
                                        msg = 'is %s vectorized because its missing _reset_indices...%s' % (
                                            self.obj.__class__.__name__, methods)
                                        break
                                        raise RuntimeError(msg)
                                    #print self.freq
                                    self.obj._reset_indices()
                                    #print "isubcase", self.isubcase
                                    #print "record_len", record
                                    # _len
                                    self.obj.ntimes += 1
                                    self.obj.ntotal = record_len // (self.num_wide * 4) * self._data_factor
                                    #self.obj.ntotal = record_len // nwide
                                    #print "ntotal         =", self.ntotal, self.obj.ntotal  #, type(self.obj)
                                    #sys.exit()
                                else:
                                    print('obj=%s doesnt have ntimes' % self.obj.__class__.__name__)
                    else:
                        data = self._read_record()
                        n = table4_parser(data)
                    del n

            self.isubtable -= 1
            self.read_markers([self.isubtable, 1, 0])
            markers = self.get_nmarkers(1, rewind=True)
        self.read_markers([0])
        self.finish()

    def is_valid_subcase(self):
        """
        Lets the code check whether or not to read a subcase

        :param self: the object pointer
        :retval is_valid: should this subcase defined by self.isubcase be read?
        """
        if not self.isAllSubcases:
            if hasattr(self, 'isubcase') and self.isubcase in self.valid_subcases:
                return True
            return False
        return True

    def goto(self, n):
        """
        Jumps to position n in the file

        :param self: the object pointer
        :param n:    the position to goto
        """
        self.n = n
        self.f.seek(n)

    def _get_record_length(self):
        len_record = 0
        n0 = self.n
        markers0 = self.get_nmarkers(1, rewind=False)

        n = self.n
        record = self.skip_block()
        len_record += self.n - n - 8  # -8 is for the block

        markers1 = self.get_nmarkers(1, rewind=True)
        # handling continuation blocks
        if markers1[0] > 0:
            #nloop = 0
            while markers1[0] > 0:
                # 'markers0=%s markers1=%s' % (markers0, markers1)
                markers1 = self.get_nmarkers(1, rewind=False)
                n = self.n
                record = self.skip_block()
                len_record += self.n - n - 8  # -8 is for the block
                markers1 = self.get_nmarkers(1, rewind=True)
        self.goto(n0)
        return len_record

    def _skip_record(self):
        markers0 = self.get_nmarkers(1, rewind=False)
        record = self.skip_block()

        markers1 = self.get_nmarkers(1, rewind=True)
        # handling continuation blocks
        if markers1[0] > 0:
            #nloop = 0
            while markers1[0] > 0:
                # 'markers0=%s markers1=%s' % (markers0, markers1)
                markers1 = self.get_nmarkers(1, rewind=False)
                record = self.read_block()
                markers1 = self.get_nmarkers(1, rewind=True)
                #nloop += 1
        return record

    def _stream_record(self, debug=True):
        markers0 = self.get_nmarkers(1, rewind=False)
        if self.debug and debug:
            self.binary_debug.write('marker = [4, %i, 4]\n' % markers0[0])
        record = self.read_block()

        if self.debug and debug:
            nrecord = len(record)
            self.binary_debug.write('record = [%i, recordi, %i]\n' % (nrecord, nrecord))
        assert (markers0[0]*4) == len(record), 'markers0=%s*4 len(record)=%s' % (markers0[0]*4, len(record))
        yield record

        markers1 = self.get_nmarkers(1, rewind=True)

        # handling continuation blocks
        if markers1[0] > 0:
            nloop = 0
            while markers1[0] > 0:
                # 'markers0=%s markers1=%s' % (markers0, markers1)
                markers1 = self.get_nmarkers(1, rewind=False)

                record = self.read_block()
                yield record
                markers1 = self.get_nmarkers(1, rewind=True)
                nloop += 1

    def _read_record(self, stream=False, debug=True):
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
            while markers1[0] > 0:
                # 'markers0=%s markers1=%s' % (markers0, markers1)
                markers1 = self.get_nmarkers(1, rewind=False)
                record = self.read_block()
                records.append(record)
                markers1 = self.get_nmarkers(1, rewind=True)
                nloop += 1

            if nloop > 0:
                record = ''.join(records)
        return record