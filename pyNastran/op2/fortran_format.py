"""
Defines:
 - FortranFormat

"""
from pyNastran.utils import object_attributes
from pyNastran.utils.numpy_utils import integer_types
#from pyNastran.op2.errors import FortranMarkerError, SortCodeError


class FortranFormat:
    """defines basic methods for reading Fortran formatted data files"""
    def __init__(self):
        self.n = 0
        self.f = None
        self.obj = None
        self.table_name = None
        self.isubcase = None
        self.binary_debug = None
        self.read_mode = 1
        self._endian = None
        self._table_mapper = {}
        self._nastran_format = None

        #: stores if the user entered [] for isubcases
        self.is_all_subcases = True
        self.valid_subcases = []
        #self.op2_reader = OP2Reader()
        self.IS_TESTING = True

    def show(self, n: int, types: str='ifs', endian=None, force: bool=False):  # pragma: no cover
        """Shows binary data"""
        return self.op2_reader.show(n, types=types, endian=endian, force=force)

    def show_data(self, data, types: str='ifs', endian=None, force: bool=False):  # pragma: no cover
        """Shows binary data"""
        return self.op2_reader.show_data(data, types=types, endian=endian, force=force)

    def show_ndata(self, n: int, types: str='ifs', force: bool=False):  # pragma: no cover
        self.op2_reader.show_ndata(n, types=types, force=force)

    #def passer(self, data):
        #"""
        #dummy function used for unsupported tables
        #"""
        #pass

    def _get_table_mapper(self):
        raise NotImplementedError('this should be overwritten')

    def _finish(self):
        raise NotImplementedError('overwrite this')

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
        op2_reader = self.op2_reader
        #datai = b''
        n = 0
        if self.read_mode == 2:
            self.ntotal = 0

            data, ndata = op2_reader._read_record_ndata()
            n = table4_parser(data, ndata)
            assert isinstance(n, integer_types), self.table_name

            self._reset_vector_counter()

        elif self.read_mode == 1:
            # if we're checking the array size

            #n = op2_reader._skip_record()
            #n = table4_parser(datai, 300000)
            if self.table_name in {b'R1TABRG', b'ONRGY1', b'PVT', b'PVT0', b'PVTS'}:
                # these tables are always fully parsed
                # PVT/PVTS - we want to know what the PARAM cards are,
                #            so we can determine the NXVER
                data, ndata = op2_reader._read_record_ndata()
            else:
                data, ndata = op2_reader._skip_record_ndata()
            n = table4_parser(data, ndata)
            if not isinstance(n, integer_types):
                msg = 'n is not an integer; table_name=%s n=%s table4_parser=%s' % (
                    self.table_name, n, table4_parser)
                raise TypeError(msg)

            #op2_reader._goto(n)
            #n = op2_reader._skip_record()

            self._init_vector_counter(record_len)
        else:
            raise RuntimeError(self.read_mode)
        self._cleanup_data_members()
        return n

    def _reset_vector_counter(self) -> None:
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
        if self.obj.ntotal == self.obj.data.shape[1] or 1:
            #if self.table_name_str in ['OESRMS2', 'OESNO2', 'OSTRRMS2', 'OSTRNO2', 'OESATO2']:
                #print('resetting %r indicies; itime=%s; shape=%s' % (
                    #self.obj.class_name, self.obj.itime, self.obj.data.shape))
            self.obj._reset_indices()
            self.obj.words = self.words
            self.obj.itime += 1
        else:
            # This happens when self._data_factor hasn't been reset
            # or is set wrong.
            # can it happen any other time?
            #
            # yup, when you have sort2...
            msga = 'self.obj.name=%r has itime' % self.obj.__class__.__name__
            self.log.debug(msga)
            msgb = 'ntotal=%s shape=%s shape[1]=%s _data_factor=%s\n' % (
                self.obj.ntotal, str(self.obj.data.shape),
                self.obj.data.shape[1], self._data_factor)
            msgb += 'obj._ntotals=%s' % self.obj._ntotals
            self.log.error(msgb)
            raise RuntimeError(msga + '\n' + msgb)

    def _init_vector_counter(self, record_len: int) -> None:
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
        return

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
            if len(self.words) not in [0, 28]:
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
        if hasattr(self, 'mode'):
            del self.mode

        for word in del_words:
            if hasattr(self, word):
                val = getattr(self, word)
                if isinstance(val, list) and len(val) == 0:
                    continue
                msg += '  %s=%s\n' % (word, val)
        if msg:
            print(object_attributes(self))
            print(msg)
