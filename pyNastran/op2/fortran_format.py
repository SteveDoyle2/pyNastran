"""
Defines:
 - FortranFormat

"""
from typing import Optional
from pyNastran.utils import object_attributes
from pyNastran.utils.numpy_utils import integer_types
#from pyNastran.op2.errors import FortranMarkerError, SortCodeError
from pyNastran.op2.errors import EmptyRecordError

from pyNastran.op2.tables.oef_forces.oef import OEF
from pyNastran.op2.tables.oes_stressStrain.oes import OES
from pyNastran.op2.tables.ogs_grid_point_stresses.ogs import OGS

from pyNastran.op2.tables.oee_energy.onr import ONR
from pyNastran.op2.tables.ogf_gridPointForces.ogpf import OGPF
#from pyNastran.op2.tables.oes_stressStrain.oesm import OESM

from pyNastran.op2.tables.opg_appliedLoads.opg import OPG
from pyNastran.op2.tables.oqg_constraintForces.oqg import OQG
from pyNastran.op2.tables.oug.oug import OUG
from pyNastran.op2.tables.oug.otemp1 import OTEMP1
from pyNastran.op2.tables.oug.ougpk1 import OUGPK1

from pyNastran.op2.tables.lama_eigenvalues.lama import LAMA
from pyNastran.op2.tables.onmd import ONMD
from pyNastran.op2.tables.opr import OPR
from pyNastran.op2.tables.ogpwg import OGPWG

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
        self._data_factor = 1

        #: stores if the user entered [] for isubcases
        self.is_all_subcases = True
        self.valid_subcases = []
        #self.op2_reader = OP2Reader()
        self.IS_TESTING = False
        self.reader_onmd = ONMD(self)
        self.reader_ogpwg = OGPWG(self)
        self.reader_lama = LAMA(self)

        self.reader_oes = OES(self)
        self.reader_opg = OPG(self)
        self.reader_oef = OEF(self)
        self.reader_oqg = OQG(self)
        self.reader_opr = OPR(self)
        self.reader_ogs = OGS(self)
        self.reader_onr = ONR(self)
        self.reader_ogpf = OGPF(self)

        # OUG - displacement, velocity, acceleration, eigenvector, temperature
        self.reader_oug = OUG(self)
        self.reader_otemp1 = OTEMP1(self)  # Siemens
        self.reader_ougpk1 = OUGPK1(self)  # STK

    def show(self, n: int, types: str='ifs', endian=None, force: bool=False):  # pragma: no cover
        """Shows binary data"""
        return self.op2_reader.show(n, types=types, endian=endian, force=force)

    def show_data(self, data, types: str='ifs', endian=None, force: bool=False):  # pragma: no cover
        """Shows binary data"""
        return self.op2_reader.show_data(data, types=types, endian=endian, force=force)

    def show_ndata(self, n: int, types: str='ifs', force: bool=False, endian=None):  # pragma: no cover
        self.op2_reader.show_ndata(n, types=types, force=force, endian=endian)

    #def passer(self, data):
        #"""
        #dummy function used for unsupported tables
        #"""
        #pass

    def _get_table_mapper(self):
        raise NotImplementedError('this should be overwritten')

    def _finish(self):
        raise NotImplementedError('overwrite this')

    def _read_subtable_results(self, table4_parser, record_len: int) -> Optional[int]:
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
            None : an error occurred or we're in read_mode=1/array sizing (???)
            int : the number of bytes that have been read

        """
        op2_reader = self.op2_reader  # type: OP2Reader
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
            #self.show(100, types='ifs', endian=None, force=False)
            if self.table_name in {b'R1TABRG', b'ONRGY1', b'PVT', b'PVT0', b'PVTS'}:
                # these tables are always fully parsed
                # PVT/PVTS - we want to know what the PARAM cards are,
                #            so we can determine the NXVER
                data, ndata = op2_reader._read_record_ndata()
            else:
                try:
                    data, ndata = op2_reader._skip_record_ndata()
                except EmptyRecordError:
                    self.log.error('error round 2...')
                    raise
                    op2_reader.read_markers([1, 0], macro_rewind=False)
                    marker146 = op2_reader.get_nmarkers4(1, rewind=True)
                    if marker146 == 146:
                        #print('marker146', marker146)
                        return n
                    self._cleanup_data_members()
                    #n = self.n
                    #return n
                    raise

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
            #print(f'self.obj.name={self.obj.__class__.__name__!r} doesnt have itime')
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
        #
        # TODO: the or 1 flag breaks tests...that could be why shapes are bad...
        #
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
            msga = f'self.obj.name={self.obj.__class__.__name__!r} has itime'
            self.log.debug(msga)
            msgb = 'ntotal=%s shape=%s shape[1]=%s _data_factor=%s\n' % (
                self.obj.ntotal, str(self.obj.data.shape),
                self.obj.data.shape[1], self._data_factor)
            msgb += f'obj._ntotals={self.obj._ntotals}'
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
            #ntotal = record_len // (self.num_wide * 4) * self._data_factor
            ntotal = record_len // (self.num_wide * self.size) * self._data_factor
            #if self._data_factor == 2 and self.num_wide == 77:
                #print(f'record_len={record_len} num_wide={self.num_wide}')
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

    def _cleanup_data_members(self) -> None:
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
                msg += f'  {word}={val}\n'
        if msg:
            print(object_attributes(self))
            print(msg)
