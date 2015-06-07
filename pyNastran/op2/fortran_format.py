from __future__ import print_function
from six.moves import range
import sys
from struct import unpack
import copy

class FortranFormat(object):
    def __init__(self):
        """
        :param self:    the OP2 object pointer
        """
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

    def show(self, n, types='ifs'):
        """
        :param self:    the OP2 object pointer
        """
        assert self.n == self.f.tell()
        nints = n // 4
        data = self.f.read(4 * n)
        strings, ints, floats = self.show_data(data, types=types)
        self.f.seek(self.n)
        return strings, ints, floats

    def show_data(self, data, types='ifs'):
        return self.write_data(sys.stdout, data, types=types)

    def write_data(self, f, data, types='ifs'):
        """
        Useful function for seeing what's going on locally when debugging.

        :param self:    the OP2 object pointer
        """
        n = len(data)
        nints = n // 4
        ndoubles = n // 8
        strings = None
        ints = None
        floats = None
        longs = None
        if 's' in types:
            strings = unpack(b'%is' % n, data)
            f.write("strings = %s\n" % str(strings))
        if 'i' in types:
            ints = unpack(b'%ii' % nints, data)
            f.write("ints    = %s\n" % str(ints))
        if 'f' in types:
            floats = unpack(b'%if' % nints, data)
            f.write("floats  = %s\n" % str(floats))

        if 'l' in types:
            longs = unpack(b'%il' % nints, data)
            f.write("long  = %s\n" % str(longs))
        if 'I' in types:
            ints2 = unpack(b'%iI' % nints, data)
            f.write("unsigned int = %s\n" % str(ints2))
        if 'L' in types:
            longs2 = unpack(b'%iL' % nints, data)
            f.write("unsigned long = %s\n" % str(longs2))
        if 'q' in types:
            longs = unpack(b'%iq' % ndoubles, data[:ndoubles*8])
            f.write("long long = %s\n" % str(longs))
        return strings, ints, floats

    def show_ndata(self, n, types='ifs'):
        return self.write_ndata(sys.stdout, n, types=types)

    def write_ndata(self, f, n, types='ifs'):
        """
        Useful function for seeing what's going on locally when debugging.

        :param self:    the OP2 object pointer
        """
        nold = self.n
        data = self.f.read(n)
        self.n = nold
        self.f.seek(self.n)
        return self.write_data(f, data, types=types)

    def skip_block(self):
        """
        Skips a block following a pattern of:
            [nbytes, data, nbytes]

        :param self:    the OP2 object pointer
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

    def read_markers(self, markers, macro_rewind=True):
        """
        Gets specified markers, where a marker has the form of [4, value, 4].
        The "marker" corresponds to the value, so 3 markers takes up 9 integers.
        These are used to indicate position in the file as well as
        the number of bytes to read.

        :param self:    the OP2 object pointer
        :param markers: markers to get; markers = [-10, 1]
        """
        for i, marker in enumerate(markers):
            data = self.read_block()
            imarker, = unpack(b'i', data)
            assert marker == imarker, 'marker=%r imarker=%r; markers=%s i=%s' % (marker, imarker, markers, i)
            self.binary_debug.write('  read_markers -> [4, %i, 4]\n' % marker)

    def get_nmarkers(self, n, rewind=True, macro_rewind=False):
        """
        Gets n markers, so if n=2, it will get 2 markers.

        :param self:    the OP2 object pointer
        :param n:        number of markers to get
        :param rewind:   should the file be returned to the starting point
        :retval markers: list of [1, 2, 3, ...] markers
        """
        ni = self.n
        markers = []
        for i in range(n):
            data = self.read_block()
            marker, = unpack(b'i', data)
            markers.append(marker)
        if rewind:
            self.n = ni
            self.f.seek(self.n)
            #for i in range(n):
                #self.binary_debug.write('get_nmarkers- [4, %i, 4]; macro_rewind=%s\n' % (i, macro_rewind or rewind))
        else:
            #if not macro_rewind:
            for i in range(n):
                self.binary_debug.write('get_nmarkers- [4, %i, 4]; macro_rewind=%s\n' % (i, macro_rewind or rewind))
        return markers

    def _skip_subtables(self):
        """
        :param self:    the OP2 object pointer
        """
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
        :param self:    the OP2 object pointer
        """
        pass

    def _get_table_mapper(self):
        """
        :param self:    the OP2 object pointer
        """
        raise NotImplementedError('this should be overwritten')

    def _read_subtables(self):
        """
        :param self:    the OP2 object pointer
        """
        # this parameters is used for numpy streaming
        self._data_factor = 1

        nstart = self.n
        self.isubtable = -3
        self.read_markers([-3, 1, 0])
        self.binary_debug.write('***isubtable = %i\n' % self.isubtable)
        self.binary_debug.write('---markers = [-3, 1, 0]---\n')
        table_mapper = self._get_table_mapper()

        # get the parsing functions (table3_parser, table4_parser)
        # or find out we're going to be skipping the tables
        #
        # table3 - the table with the meta data (e.g. subcaseID, time, is_stress/strain)
        # table4 - the actual results data
        #
        # we indicate table3/4 by isubtable, which starts from -3 (so table3) and counts
        # down (yes down) to 4 to indicate table4.  If we count down again, we end up
        # back at table 3 (with isubtable=-5), which will occur in the case of multiple
        # times/element types/results in a single macro table (e.g. OUG, OES).
        if self.table_name in table_mapper:
            #if self.read_mode in [0, 2]:
                #self.log.debug("table_name = %r" % self.table_name)

            table3_parser, table4_parser = table_mapper[self.table_name]
            passer = False
        else:
            if self.read_mode in [0, 2]:
                self.log.debug("skipping table_name = %r" % self.table_name)
                #raise NotImplementedError(self.table_name)
            table3_parser = None
            table4_parser = None
            passer = True

        # we need to check the marker, so we read it and rewind, so we don't
        # screw up our positioning in the file
        markers = self.get_nmarkers(1, rewind=True)
        self.binary_debug.write('---marker0 = %s---\n' % markers)

        # while the subtables aren't done
        while markers[0] != 0:
            self.is_start_of_subtable = True
            self.binary_debug.write('***isubtable = %i\n' % self.isubtable)
            self._read_subtable_3_4(table3_parser, table4_parser, passer)
            self.isubtable -= 1
            self.read_markers([self.isubtable, 1, 0])
            markers = self.get_nmarkers(1, rewind=True)

        # we've finished reading all subtables, but have one last marker to read
        self.read_markers([0])
        self.finish()

    def _read_subtable_3_4(self, table3_parser, table4_parser, passer):
        # this is the length of the current record inside table3/table4
        record_len = self._get_record_length()
        self.binary_debug.write('record_length = %s\n' % record_len)
        if record_len == 584:  # table3 has a length of 584
            self.data_code = {}
            self.obj = None
            data = self._read_record()
            if not passer:
                table3_parser(data)
                #if hasattr(self, 'isubcase'):
                    #print("code = ", self._get_code())
        else:
            if passer or not self.is_valid_subcase():
                data = self._skip_record()
            else:
                if hasattr(self, 'num_wide'):
                    # num_wide is the result size and is usually found in
                    # table3, but some B-list tables don't have it
                    n = self._read_subtable_results(table4_parser, record_len)
                else:
                    data = self._read_record()
                    n = table4_parser(data)
                #del n

    def _read_subtable_results(self, table4_parser, record_len):
        """
        # if reading the data
        # 0 - non-vectorized
        # 1 - 1st pass to size the array (vectorized)
        # 2 - 2nd pass to read the data  (vectorized)
        """
        datai = b''
        n = 0
        if self.read_mode in [0, 2]:
            self.ntotal = 0

            # we stream the record because we get it in partial blocks
            for data in self._stream_record():
                data = datai + data

                n = table4_parser(data)
                assert isinstance(n, int), self.table_name
                datai = data[n:]

            # PCOMPs are stupid, so we need an element flag
            if hasattr(self, 'eid_old'):
                del self.eid_old

            # if reading the data
            # 0 - non-vectorized
            # 1 - 1st pass to size the array (vectorized)
            # 2 - 2nd pass to read the data  (vectorized)
            if self.read_mode == 2:
                # vectorized objects are stored as self.obj
                # they have obj.itime which is their table3 counter
                if hasattr(self, 'obj') and hasattr(self.obj, 'itime'):
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
                    if self.obj.ntotal == self.obj.data.shape[1]:
                        self.obj._reset_indices()
                        self.obj.words = self.words
                        self.obj.itime += 1
                    else:
                        print('self.obj.name=%r has itime' % self.obj.__class__.__name__)
                        print('ntotal=%s shape=%s shape[1]=%s _data_factor=%s' % (
                            self.obj.ntotal, str(self.obj.data.shape),
                            self.obj.data.shape[1], self._data_factor))
                #else:
                    #print('self.obj.name=%r doesnt have itime' % self.obj.__class__.__name__)

        elif self.read_mode == 1:
            # if we're checking the array size

            #n = self._skip_record()
            #n = table4_parser(datai, 300000)
            if 1:
                self.ntotal = 0
                #n = self.n
                n = 0
                for i, data in enumerate(self._stream_record()):
                    data = datai + data
                    ndata = len(data)
                    n = table4_parser(data)
                    assert isinstance(n, int), self.table_name
                    datai = data[n:]
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
                        return None
                        raise RuntimeError(msg)
                    self.obj._reset_indices()
                    self.obj.ntimes += 1
                    self.obj.ntotal = record_len // (self.num_wide * 4) * self._data_factor
                else:
                    print('obj=%s doesnt have ntimes' % self.obj.__class__.__name__)
        else:
            raise RuntimeError(self.read_mode)
        from pyNastran.utils import object_attributes
        del_words = [
            'words',
            #'Title',
            #'ID',
            'analysis_code',
            #'result_names',
            #'labels',
            #'dataNames',
        ]
        msg = ''
        if hasattr(self, 'words'):
            for word in self.words:
                if word in ['???', 'Title']:
                    continue
                if not hasattr(self, word):
                    continue
                delattr(self, word)
            self.words = []
        if hasattr(self, 'analysis_code'):
            del self.analysis_code
        #if hasattr(self, 'dataNames') and self.dataNames is not None:
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
        return n

    def is_valid_subcase(self):
        """
        Lets the code check whether or not to read a subcase

        :param self:    the OP2 object pointer
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

        :param self:    the OP2 object pointer
        :param n:    the position to goto
        """
        self.n = n
        self.f.seek(n)

    def _get_record_length(self):
        """
        The record length helps us figure out data block size, which is used
        to quickly size the arrays.  We just need a bit of meta data and can
        jump around quickly.
        """
        self.binary_debug.write('_get_record_length\n')
        len_record = 0
        n0 = self.n
        markers0 = self.get_nmarkers(1, rewind=False)
        self.binary_debug.write('  markers0=%s\n' % markers0)

        n = self.n
        record = self.skip_block()
        len_record += self.n - n - 8  # -8 is for the block
        self.binary_debug.write('  len_record=%s\n' % len_record)

        markers1 = self.get_nmarkers(1, rewind=True)
        # handling continuation blocks
        if markers1[0] > 0:
            while markers1[0] > 0:
                markers1 = self.get_nmarkers(1, rewind=False)
                self.binary_debug.write('  markers1=%s\n' % markers1)
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
            while markers1[0] > 0:
                markers1 = self.get_nmarkers(1, rewind=False)
                record = self.read_block()
                markers1 = self.get_nmarkers(1, rewind=True)
        return record

    def _stream_record(self, debug=True):
        """
        Creates a "for" loop that keeps giving us records until we're done.

        :param self:    the OP2 object pointer
        """
        self.istream = 0
        markers0 = self.get_nmarkers(1, rewind=False)
        if self.debug and debug:
            self.binary_debug.write('_stream_record - marker = [4, %i, 4]\n' % markers0[0])
        record = self.read_block()
        if self.debug and debug:
            nrecord = len(record)
            self.binary_debug.write('_stream_record - record = [%i, recordi, %i]\n' % (nrecord, nrecord))
        assert (markers0[0]*4) == len(record), 'markers0=%s*4 len(record)=%s' % (markers0[0]*4, len(record))
        yield record
        self.istream += 1

        markers1 = self.get_nmarkers(1, rewind=True)
        if self.debug and debug:
            self.binary_debug.write('_stream_record - markers1 = [4, %s, 4]\n' % str(markers1))

        # handling continuation blocks
        if markers1[0] > 0:
            nloop = 0
            while markers1[0] > 0:
                markers1 = self.get_nmarkers(1, rewind=False)
                if self.debug and debug:
                    self.binary_debug.write('_stream_record - markers1 = [4, %s, 4]\n' % str(markers1))
                record = self.read_block()
                yield record
                self.istream += 1

                markers1 = self.get_nmarkers(1, rewind=True)
                nloop += 1

    def _read_record(self, stream=False, debug=True, macro_rewind=False):
        """
        :param self:  the OP2 object pointer
        """
        markers0 = self.get_nmarkers(1, rewind=False, macro_rewind=macro_rewind)
        if self.debug and debug:
            self.binary_debug.write('read_record - marker = [4, %i, 4]; macro_rewind=%s\n' % (markers0[0], macro_rewind))
        record = self.read_block()

        if self.debug and debug:
            nrecord = len(record)
            self.binary_debug.write('read_record - record = [%i, recordi, %i]; macro_rewind=%s\n' % (nrecord, nrecord, macro_rewind))
        assert (markers0[0]*4) == len(record), 'markers0=%s*4 len(record)=%s' % (markers0[0]*4, len(record))

        markers1 = self.get_nmarkers(1, rewind=True)

        # handling continuation blocks
        if markers1[0] > 0:
            nloop = 0
            records = [record]
            while markers1[0] > 0:
                markers1 = self.get_nmarkers(1, rewind=False)
                if self.debug and debug:
                    self.binary_debug.write('read_record - markers1 = [4, %i, 4]\n' % markers1[0])
                record = self.read_block()
                records.append(record)
                markers1 = self.get_nmarkers(1, rewind=True)
                if self.debug and debug:
                    self.binary_debug.write('read_record - markers1 = [4, %i, 4]\n' % markers1[0])
                nloop += 1

            if nloop > 0:
                record = b''.join(records)
        return record
