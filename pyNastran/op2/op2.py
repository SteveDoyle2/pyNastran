#pylint: disable=C0301,C0103,W0613,C0111,W0612,R0913
"""
Defines the OP2 class.
"""
import os
from struct import unpack, Struct
#from struct import error as StructError

from numpy import array

from pyNastran.f06.errors import FatalError
from pyNastran.f06.tables.grid_point_weight import GridPointWeight

#from pyNastran.op2.tables.geom.geom1 import GEOM1
#from pyNastran.op2.tables.geom.geom2 import GEOM2
#from pyNastran.op2.tables.geom.geom3 import GEOM3
#from pyNastran.op2.tables.geom.geom4 import GEOM4

#from pyNastran.op2.tables.geom.ept import EPT
#from pyNastran.op2.tables.geom.mpt import MPT

#from pyNastran.op2.tables.geom.dit import DIT
#from pyNastran.op2.tables.geom.dynamics import DYNAMICS
#============================

from pyNastran.op2.tables.lama_eigenvalues.lama import LAMA
from pyNastran.op2.tables.oee_energy.onr import ONR
from pyNastran.op2.tables.opg_appliedLoads.ogpf import OGPF

from pyNastran.op2.tables.oef_forces.oef import OEF
from pyNastran.op2.tables.oes_stressStrain.oes import OES
from pyNastran.op2.tables.ogs import OGS

from pyNastran.op2.tables.opg_appliedLoads.opg import OPG
from pyNastran.op2.tables.oqg_constraintForces.oqg import OQG
from pyNastran.op2.tables.oug.oug import OUG
from pyNastran.op2.tables.ogpwg import OGPWG
from pyNastran.op2.fortran_format import FortranFormat

#from pyNastran.bdf.bdf import BDF

from pyNastran.utils import is_binary
from pyNastran.utils.log import get_logger


class TrashWriter(object):
    """
    A dummy file that just trashes all data
    """
    def __init__(self, *args, **kwargs):
        pass
    def open(self, *args, **kwargs):
        pass
    def write(self, *args, **kwargs):
        pass
    def close(self, *args, **kwargs):
        pass

class OP2( #BDF,
          #GEOM1, GEOM2, GEOM3, GEOM4, EPT, MPT, DIT, DYNAMICS,
          LAMA, ONR, OGPF,
          OEF, OES, OGS, OPG, OQG, OUG, OGPWG, FortranFormat):
    """
    Defines an interface for the Nastran OP2 file.
    """
    def set_as_vectorized(self, vectorized=False, ask=False):
        if vectorized is True:
            msg = 'OP2 class doesnt support vectorization.  Use OP2_Vectorized '
            msg += 'from pyNastran.op2.dev_explicit.op2_vectorized instead.'
            raise RuntimeError(msg)
        if ask is True:
            msg = 'OP2 class doesnt support ask.'
            raise RuntimeError(msg)

    def set_subcases(self, subcases=None):
        """
        Allows you to read only the subcases in the list of iSubcases

        :param subcases: list of [subcase1_ID,subcase2_ID]
                         (default=None; all subcases)
        """
        #: stores the set of all subcases that are in the OP2
        self.subcases = set()
        if subcases is None or subcases == []:
            #: stores if the user entered [] for iSubcases
            self.isAllSubcases = True
            self.valid_subcases = []
        else:
            #: should all the subcases be read (default=True)
            self.isAllSubcases = False

            #: the set of valid subcases -> set([1,2,3])
            self.valid_subcases = set(subcases)
        self.log.debug("set_subcases - subcases = %s" % self.valid_subcases)

    def set_transient_times(self, times):  # TODO this name sucks...
        """
        Takes a dictionary of list of times in a transient case and
        gets the output closest to those times.::

          times = {subcaseID_1: [time1, time2],
                   subcaseID_2: [time3, time4]}
        """
        expected_times = {}
        for (isubcase, eTimes) in times.iteritems():
            eTimes = list(times)
            eTimes.sort()
            expected_times[isubcase] = array(eTimes)
        self.expected_times = expected_times

    def __init__(self, make_geom=False,
                 debug=False, log=None, debug_file=None):
        """
        Initializes the OP2 object

        :param make_geom: reads the BDF tables (default=False)
        :param debug: enables the debug log and sets the debug in the logger (default=False)
        :param log: a logging object to write debug messages to
         (.. seealso:: import logging)
        :param debug_file: sets the filename that will be written to (default=None -> no debug)
        """
        assert isinstance(make_geom, bool), 'make_geom=%r' % make_geom
        assert isinstance(debug, bool), 'debug=%r' % debug

        self.make_geom = make_geom

        #BDF.__init__(self, debug=debug, log=log)
        self.log = get_logger(log, 'debug' if debug else 'info')

        #GEOM1.__init__(self)
        #GEOM2.__init__(self)
        #GEOM3.__init__(self)
        #GEOM4.__init__(self)

        #EPT.__init__(self)
        #MPT.__init__(self)
        #DIT.__init__(self)
        #DYNAMICS.__init__(self)

        LAMA.__init__(self)
        ONR.__init__(self)
        OGPF.__init__(self)

        OEF.__init__(self)
        OES.__init__(self)
        OGS.__init__(self)

        OPG.__init__(self)
        OQG.__init__(self)
        OUG.__init__(self)
        OGPWG.__init__(self)
        FortranFormat.__init__(self)

        self.read_mode = 0
        self.is_vectorized = False
        self._close_op2 = True

        self.result_names = set([])

        self.grid_point_weight = GridPointWeight()
        self.words = []
        self.debug = debug
        #self.debug = True
        #self.debug = False
        #debug_file = None
        if debug_file is None:
            self.debug_file = None
        else:
            assert isinstance(debug_file, basestring), debug_file
            self.debug_file = debug_file
        self.make_geom = False

    def _get_table_mapper(self):
        table_mapper = {
            #=======================
            # OEF
            # element forces
            'OEFIT' : [self._read_oef1_3, self._read_oef1_4],
            'OEF1X' : [self._read_oef1_3, self._read_oef1_4],
            'OEF1'  : [self._read_oef1_3, self._read_oef1_4],
            'DOEF1' : [self._read_oef1_3, self._read_oef1_4],
            #=======================
            # OQG
            # spc forces
            'OQG1'  : [self._read_oqg1_3, self._read_oqg1_4],  # spc forces in the nodal frame
            'OQGV1' : [self._read_oqg1_3, self._read_oqg1_4],  # spc forces in the nodal frame
            # mpc forces
            'OQMG1' : [self._read_oqg1_3, self._read_oqg1_4],

            # ???? - passer
            #'OQP1': [self._table_passer, self._table_passer],
            #=======================
            # OPG
            # applied loads
            'OPG1'  : [self._read_opg1_3, self._read_opg1_4],  # applied loads in the nodal frame
            'OPGV1' : [self._read_opg1_3, self._read_opg1_4],
            'OPNL1' : [self._read_opg1_3, self._read_opg1_4],

            # OGPFB1
            # grid point forces
            'OGPFB1' : [self._read_ogpf1_3, self._read_ogpf1_4],  # grid point forces

            # ONR/OEE
            # strain energy density
            'ONRGY1' : [self._read_onr1_3, self._read_onr1_4],  # strain energy density
            #=======================
            # OES
            # stress
            'OES1X1'  : [self._read_oes1_3, self._read_oes1_4],  # stress
            'OES1'    : [self._read_oes1_3, self._read_oes1_4],  # stress
            'OES1X'   : [self._read_oes1_3, self._read_oes1_4],  # stress
            'OES1C'   : [self._read_oes1_3, self._read_oes1_4],  # stress - composite
            'OESCP'   : [self._read_oes1_3, self._read_oes1_4],
            'OESNLXR' : [self._read_oes1_3, self._read_oes1_4],
            'OESNLXD' : [self._read_oes1_3, self._read_oes1_4],
            'OESNLBR' : [self._read_oes1_3, self._read_oes1_4],
            'OESTRCP' : [self._read_oes1_3, self._read_oes1_4],
            'OESNL1X' : [self._read_oes1_3, self._read_oes1_4],
            'OESRT'   : [self._read_oes1_3, self._read_oes1_4],

            # strain
            'OSTR1X'  : [self._read_oes1_3, self._read_oes1_4],  # strain - isotropic
            'OSTR1C'  : [self._read_oes1_3, self._read_oes1_4],  # strain - composite

            #=======================
            # OUG
            # displacement/velocity/acceleration/eigenvector/temperature
            'OUG1'    : [self._read_oug1_3, self._read_oug1_4],  # displacements in nodal frame
            'OUGV1'   : [self._read_oug1_3, self._read_oug1_4],  # displacements in nodal frame
            'BOUGV1'  : [self._read_oug1_3, self._read_oug1_4],  # OUG1 on the boundary???
            'OUGV1PAT': [self._read_oug1_3, self._read_oug1_4],  # OUG1 + coord ID
            'OUPV1'   : [self._read_oug1_3, self._read_oug1_4],

            #=======================
            # OGPWG
            # grid point weight
            'OGPWG'  : [self._read_ogpwg_3, self._read_ogpwg_4],  # grid point weight
            'OGPWGM' : [self._read_ogpwg_3, self._read_ogpwg_4],  # modal? grid point weight

            #=======================
            # OGS
            # grid point stresses
            'OGS1' : [self._read_ogs1_3, self._read_ogs1_4],  # grid point stresses
            #=======================
            # eigenvalues
            'BLAMA': [self._read_buckling_eigenvalue_3, self._read_buckling_eigenvalue_4],  # buckling eigenvalues
            'CLAMA': [self._read_complex_eigenvalue_3,  self._read_complex_eigenvalue_4],   # complex eigenvalues
            'LAMA' : [self._read_real_eigenvalue_3,     self._read_real_eigenvalue_4],      # eigenvalues

            # ===geom passers===
            # geometry
            'GEOM1': [self._table_passer, self._table_passer],
            'GEOM2': [self._table_passer, self._table_passer],
            'GEOM3': [self._table_passer, self._table_passer],
            'GEOM4': [self._table_passer, self._table_passer],

            # superelements
            'GEOM1S': [self._table_passer, self._table_passer],  # GEOMx + superelement
            'GEOM2S': [self._table_passer, self._table_passer],
            'GEOM3S': [self._table_passer, self._table_passer],
            'GEOM4S': [self._table_passer, self._table_passer],

            'GEOM1N': [self._table_passer, self._table_passer],
            'GEOM2N': [self._table_passer, self._table_passer],
            'GEOM3N': [self._table_passer, self._table_passer],
            'GEOM4N': [self._table_passer, self._table_passer],

            'GEOM1OLD': [self._table_passer, self._table_passer],
            'GEOM2OLD': [self._table_passer, self._table_passer],
            'GEOM3OLD': [self._table_passer, self._table_passer],
            'GEOM4OLD': [self._table_passer, self._table_passer],

            'EPT' : [self._table_passer, self._table_passer],  # elements
            'EPTS': [self._table_passer, self._table_passer],  # elements - superelements
            'EPTOLD' : [self._table_passer, self._table_passer],

            'MPT' : [self._table_passer, self._table_passer],  # materials
            'MPTS': [self._table_passer, self._table_passer],  # materials - superelements

            'DYNAMIC': [self._table_passer, self._table_passer],
            'DYNAMICS': [self._table_passer, self._table_passer],
            'DIT': [self._table_passer, self._table_passer],

            # geometry
            #'GEOM1': [self._read_geom1_4, self._read_geom1_4],
            #'GEOM2': [self._read_geom2_4, self._read_geom2_4],
            #'GEOM3': [self._read_geom3_4, self._read_geom3_4],
            #'GEOM4': [self._read_geom4_4, self._read_geom4_4],

            # superelements
            #'GEOM1S': [self._read_geom1_4, self._read_geom1_4],
            #'GEOM2S': [self._read_geom2_4, self._read_geom2_4],
            #'GEOM3S': [self._read_geom3_4, self._read_geom3_4],
            #'GEOM4S': [self._read_geom4_4, self._read_geom4_4],

            #'GEOM1N': [self._read_geom1_4, self._read_geom1_4],
            #'GEOM2N': [self._read_geom2_4, self._read_geom2_4],
            #'GEOM3N': [self._read_geom3_4, self._read_geom3_4],
            #'GEOM4N': [self._read_geom4_4, self._read_geom4_4],

            #'GEOM1OLD': [self._read_geom1_4, self._read_geom1_4],
            #'GEOM2OLD': [self._read_geom2_4, self._read_geom2_4],
            #'GEOM3OLD': [self._read_geom3_4, self._read_geom3_4],
            #'GEOM4OLD': [self._read_geom4_4, self._read_geom4_4],

            #'EPT' : [self._read_ept_4, self._read_ept_4],
            #'EPTS': [self._read_ept_4, self._read_ept_4],
            #'EPTOLD' : [self._read_ept_4, self._read_ept_4],

            #'MPT' : [self._read_mpt_4, self._read_mpt_4],
            #'MPTS': [self._read_mpt_4, self._read_mpt_4],

            #'DYNAMIC': [self._read_dynamics_4, self._read_dynamics_4],
            #'DYNAMICS': [self._read_dynamics_4, self._read_dynamics_4],
            #'DIT': [self._read_dit_4, self._read_dit_4],   # table objects (e.g. TABLED1)

            # ===passers===
            'EQEXIN': [self._table_passer, self._table_passer],
            'EQEXINS': [self._table_passer, self._table_passer],

            'GPDT': [self._table_passer, self._table_passer],     # grid points?
            'BGPDT': [self._table_passer, self._table_passer],
            'BGPDTS': [self._table_passer, self._table_passer],
            'BGPDTOLD': [self._table_passer, self._table_passer],

            'PVT0': [self._table_passer, self._table_passer],
            'DESTAB': [self._table_passer, self._table_passer],
            'STDISP': [self._table_passer, self._table_passer],
            'R1TABRG': [self._table_passer, self._table_passer],
            'CASECC': [self._table_passer, self._table_passer],  # case control deck

            'HISADD': [self._table_passer, self._table_passer],  # optimization history (SOL200)
            'EDTS': [self._table_passer, self._table_passer],
            'FOL': [self._table_passer, self._table_passer],
            'MONITOR': [self._table_passer, self._table_passer],  # monitor points
            'PERF': [self._table_passer, self._table_passer],
            'VIEWTB': [self._table_passer, self._table_passer],   # view elements

            #==================================
            'OUGATO2': [self._table_passer, self._table_passer],
            'OUGCRM2': [self._table_passer, self._table_passer],
            'OUGNO2': [self._table_passer, self._table_passer],
            'OUGPSD2': [self._table_passer, self._table_passer],  # psd
            'OUGRMS2': [self._table_passer, self._table_passer],  # rms

            'OQGATO2': [self._table_passer, self._table_passer],
            'OQGCRM2': [self._table_passer, self._table_passer],

            'OQGNO2': [self._table_passer, self._table_passer],
            'OQGPSD2': [self._table_passer, self._table_passer],
            'OQGRMS2': [self._table_passer, self._table_passer],

            'OFMPF2M': [self._table_passer, self._table_passer],
            'OLMPF2M': [self._table_passer, self._table_passer],
            'OPMPF2M': [self._table_passer, self._table_passer],
            'OSMPF2M': [self._table_passer, self._table_passer],
            'OGPMPF2M': [self._table_passer, self._table_passer],

            'OEFATO2': [self._table_passer, self._table_passer],
            'OEFCRM2': [self._table_passer, self._table_passer],
            'OEFNO2': [self._table_passer, self._table_passer],
            'OEFPSD2': [self._table_passer, self._table_passer],
            'OEFRMS2': [self._table_passer, self._table_passer],

            'OESATO2': [self._table_passer, self._table_passer],
            'OESCRM2': [self._table_passer, self._table_passer],
            'OESNO2': [self._table_passer, self._table_passer],
            'OESPSD2': [self._table_passer, self._table_passer],
            'OESRMS2': [self._table_passer, self._table_passer],

            'OVGATO2': [self._table_passer, self._table_passer],
            'OVGCRM2': [self._table_passer, self._table_passer],
            'OVGNO2': [self._table_passer, self._table_passer],
            'OVGPSD2': [self._table_passer, self._table_passer],
            'OVGRMS2': [self._table_passer, self._table_passer],

            #==================================
            #'GPL': [self._table_passer, self._table_passer],
            'OMM2': [self._table_passer, self._table_passer],
            'ERRORN': [self._table_passer, self._table_passer],
            #==================================

            'OCRPG': [self._table_passer, self._table_passer],
            'OCRUG': [self._table_passer, self._table_passer],

            'EDOM': [self._table_passer, self._table_passer],

            'OAGPSD2': [self._table_passer, self._table_passer],
            'OAGATO2': [self._table_passer, self._table_passer],
            'OAGRMS2': [self._table_passer, self._table_passer],
            'OAGNO2': [self._table_passer, self._table_passer],
            'OAGCRM2': [self._table_passer, self._table_passer],

            'OPGPSD2': [self._table_passer, self._table_passer],
            'OPGATO2': [self._table_passer, self._table_passer],
            'OPGRMS2': [self._table_passer, self._table_passer],
            'OPGNO2': [self._table_passer, self._table_passer],
            'OPGCRM2': [self._table_passer, self._table_passer],

            'OSTRPSD2': [self._table_passer, self._table_passer],
            'OSTRATO2': [self._table_passer, self._table_passer],
            'OSTRRMS2': [self._table_passer, self._table_passer],
            'OSTRNO2': [self._table_passer, self._table_passer],
            'OSTRCRM2': [self._table_passer, self._table_passer],

            'OQMPSD2': [self._table_passer, self._table_passer],
            'OQMATO2': [self._table_passer, self._table_passer],
            'OQMRMS2': [self._table_passer, self._table_passer],
            'OQMNO2': [self._table_passer, self._table_passer],
            'OQMCRM2': [self._table_passer, self._table_passer],

            'AAA': [self._table_passer, self._table_passer],
            'AAA': [self._table_passer, self._table_passer],
            'AAA': [self._table_passer, self._table_passer],
            'AAA': [self._table_passer, self._table_passer],
        }
        return table_mapper

    def _not_available(self, data):
        if len(data) > 0:
            raise RuntimeError('this should never be called...table_name=%r len(data)=%s' % (self.table_name, len(data)))

    def _table_passer(self, data):
        return len(data)

    def readFake(self, data, n):
        return n

    def read_op2(self, op2_filename=None):
        """
        Starts the OP2 file reading

        :param op2_filename: if a string is set, the filename specified in the
            __init__ is overwritten.

            init=None,  op2_filename=None   -> a dialog is popped up  (not implemented; crash)
            init=fname, op2_filename=fname  -> fname is used
        """
        sr = list(self._saved_results)
        sr.sort()
        self.log.debug('_saved_results = %s' % str(sr))
        if op2_filename is None:
            from pyNastran.utils.gui_io import load_file_dialog
            wildcard_wx = "Nastran OP2 (*.op2)|*.op2|" \
                "All files (*.*)|*.*"
            wildcard_qt = "Nastran OP2 (*.op2);;All files (*)"
            title = 'Please select a OP2 to load'
            op2_filename = load_file_dialog(title, wildcard_wx, wildcard_qt)
            assert op2_filename is not None, op2_filename

        if not is_binary(op2_filename):
            if os.path.getsize(op2_filename) == 0:
                raise IOError('op2_filename=%r is empty.' % op2_filename)
            raise IOError('op2_filename=%r is not a binary OP2.' % op2_filename)

        self.log.debug('op2_filename = %r' % op2_filename)
        bdf_extension = '.bdf'
        f06_extension = '.f06'
        (fname, extension) = os.path.splitext(op2_filename)

        self.op2_filename = op2_filename
        self.bdf_filename = fname + bdf_extension
        self.f06_filename = fname + f06_extension

        if self.debug_file is not None:
            #: an ASCII version of the op2 (creates lots of output)
            self.binary_debug = open(self.debug_file, 'wb')
            self.binary_debug.write(self.op2_filename + '\n')
        else:
            self.binary_debug = open(os.devnull, 'wb')  #TemporaryFile()
            self.binary_debug = TrashWriter('debug.out', 'wb')

        #: file index
        self.n = 0
        self.table_name = None

        #: the OP2 file object
        self.f = open(self.op2_filename, 'rb')
        try:
            markers = self.get_nmarkers(1, rewind=True)
        except:
            self.goto(0)
            try:
                self.f.read(4)
            except:
                raise FatalError("The OP2 is empty.")
            raise

        if markers == [3,]:  # PARAM, POST, -1
            self.read_markers([3])
            data = self.read_block()

            self.read_markers([7])
            data = self.read_block()
            #self.show(100)
            data = self._read_record()
            self.read_markers([-1, 0])
        elif markers == [2,]:  # PARAM, POST, -2
            pass
        else:
            raise NotImplementedError(markers)

        #=================
        table_name = self.read_table_name(rewind=True, stop_on_failure=False)
        if table_name is None:
            raise FatalError('no tables exists...')

        table_names = self._read_tables(table_name)
        if self.debug:
            self.binary_debug.write('-' * 80 + '\n')
            self.binary_debug.write('f.tell()=%s\ndone...\n' % self.f.tell())
        self.binary_debug.close()
        if self._close_op2:
            self.f.close()
            del self.binary_debug
            del self.f
        #self.remove_unpickable_data()
        return table_names

    def create_unpickable_data(self):
        raise NotImplementedError()
        #==== not needed ====
        #self.f
        #self.binary_debug

        # needed
        self._geom1_map
        self._geom2_map
        self._geom3_map
        self._geom4_map
        self._dit_map
        self._dynamics_map
        self._ept_map
        self._mpt_map
        self._table_mapper

    def remove_unpickable_data(self):
        del self.f
        del self.binary_debug
        del self._geom1_map
        del self._geom2_map
        del self._geom3_map
        del self._geom4_map
        del self._dit_map
        del self._dynamics_map
        del self._ept_map
        del self._mpt_map
        del self._table_mapper

    def _read_tables(self, table_name):
        """
        Reads all the geometry/result tables.
        The OP2 header is not read by this function.

        :param table_name: the first table's name
        """
        table_names = []
        while table_name is not None:
            #print "----------------------------------"
            table_names.append(table_name)

            if self.debug:
                self.binary_debug.write('-' * 80 + '\n')
                self.binary_debug.write('table_name = %r; f.tell()=%s\n' % (table_name, self.f.tell()))

            self.table_name = table_name
            if 0:
                self._skip_table(table_name)
            else:
                if table_name in ['GEOM1', 'GEOM2', 'GEOM3', 'GEOM4',  # regular
                                  'GEOM1S', 'GEOM2S', 'GEOM3S', 'GEOM4S', # superelements
                                  'GEOM1N',
                                  'GEOM1OLD', 'GEOM2OLD', 'GEOM4OLD',

                                  'EPT', 'EPTS', 'EPTOLD',
                                  'EDTS',
                                  'MPT', 'MPTS',

                                  'PVT0', 'CASECC',
                                  'EDOM', 'OGPFB1',
                                  'BGPDT', 'BGPDTS', 'BGPDTOLD',
                                  'DYNAMIC', 'DYNAMICS',
                                  'EQEXIN', 'EQEXINS',
                                  'GPDT', 'ERRORN',
                                  'DESTAB', 'R1TABRG', 'HISADD',

                                   # eigenvalues
                                   'BLAMA', 'LAMA',
                                   # strain energy
                                   'ONRGY1',
                                   # grid point weight
                                   'OGPWG', 'OGPWGM',

                                   # other
                                   'CONTACT', 'VIEWTB',
                                   'KDICT', 'MONITOR', 'PERF',
                                  ]:
                    self._read_geom_table()  # DIT (agard)
                elif table_name in ['GPL', ]:
                    self._read_gpl()
                elif table_name in ['MEFF', ]:
                    self._read_meff()
                elif table_name in ['OMM2', ]:
                    self._read_omm2()
                elif table_name in ['DIT']:  # tables
                    self._read_dit()
                elif table_name in ['KELM']:
                    self._read_kelm()
                elif table_name in ['PCOMPTS']: # blade
                    self._read_pcompts()
                elif table_name == 'FOL':
                    self._read_fol()
                elif table_name in ['SDF', 'PMRF']:  #, 'PERF'
                    self._read_sdf()
                elif table_name in [
                                    # stress
                                    'OES1X1', 'OES1', 'OES1X', 'OES1C', 'OESCP',
                                    'OESNLXR', 'OESNLXD', 'OESNLBR', 'OESTRCP',
                                    'OESNL1X', 'OESRT',
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
                                    # displacement/velocity/acceleration/eigenvector/temperature
                                    'OUG1', 'OUGV1', 'BOUGV1', 'OUPV1', 'OUGV1PAT',
                                    # applied loads
                                    'OPG1', 'OPGV1', 'OPNL1', #'OPG2',

                                    # grid point stresses
                                    'OGS1',

                                    # other
                                    'OPNL1', 'OFMPF2M',
                                    'OSMPF2M', 'OPMPF2M', 'OLMPF2M', 'OGPMPF2M',

                                    'OAGPSD2', 'OAGCRM2', 'OAGRMS2', 'OAGATO2', 'OAGNO2',
                                    'OESPSD2', 'OESCRM2', 'OESRMS2', 'OESATO2', 'OESNO2',
                                    'OEFPSD2', 'OEFCRM2', 'OEFRMS2', 'OEFATO2', 'OEFNO2',
                                    'OPGPSD2', 'OPGCRM2', 'OPGRMS2', 'OPGATO2', 'OPGNO2',
                                    'OQGPSD2', 'OQGCRM2', 'OQGRMS2', 'OQGATO2', 'OQGNO2',
                                    'OQMPSD2', 'OQMCRM2', 'OQMRMS2', 'OQMATO2', 'OQMNO2',
                                    'OUGPSD2', 'OUGCRM2', 'OUGRMS2', 'OUGATO2', 'OUGNO2',
                                    'OVGPSD2', 'OVGCRM2', 'OVGRMS2', 'OVGATO2', 'OVGNO2',
                                    'OSTRPSD2', 'OSTRCRM2', 'OSTRRMS2', 'OSTRATO2', 'OSTRNO2',
                                    'OCRUG',
                                    'OCRPG',
                                    'STDISP',
                                    ]:
                    self._read_results_table()
                else:
                    raise NotImplementedError('%r' % table_name)

            table_name = self.read_table_name(rewind=True, stop_on_failure=False)
            #if table_name is None:
                #self.show(100)

        #self.show_data(data)
        #print "----------------"
        #print "done..."
        return table_names

    def _skip_table(self, table_name):
        """bypasses the next table as quickly as possible"""
        if table_name in ['DIT']:  # tables
            self._read_dit()
        elif table_name in ['PCOMPTS']:
            self._read_pcompts()
        else:
            self._skip_table_helper()

    def _read_dit(self):
        """
        Reads the DIT table (poorly).
        The DIT table stores information about table cards (e.g. TABLED1, TABLEM1).
        """
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

    def _read_kelm(self):
        """
        ..todo:: this table follows a totally different pattern...
        The KELM table stores information about the K matrix???
        """
        self.log.debug("table_name = %r" % self.table_name)
        table_name = self.read_table_name(rewind=False)
        self.read_markers([-1])
        data = self._read_record()

        self.read_markers([-2, 1, 0])
        data = self._read_record()
        if len(data) == 16:  # KELM
            table_name, dummyA, dummyB = unpack(b'8sii', data)
            assert dummyA == 170, dummyA
            assert dummyB == 170, dummyB
        else:
            raise NotImplementedError(self.show_data(data))

        self.read_markers([-3, 1, 1])

        # data = read_record()
        for i in xrange(170//10):
            self.read_markers([2])
            data = self.read_block()
            #print "i=%s n=%s" % (i, self.n)

        self.read_markers([4])
        data = self.read_block()

        for i in xrange(7):
            self.read_markers([2])
            data = self.read_block()
            #print "i=%s n=%s" % (i, self.n)


        self.read_markers([-4, 1, 1])
        for i in xrange(170//10):
            self.read_markers([2])
            data = self.read_block()
            #print "i=%s n=%s" % (i, self.n)

        self.read_markers([4])
        data = self.read_block()

        for i in xrange(7):
            self.read_markers([2])
            data = self.read_block()
            #print "i=%s n=%s" % (i, self.n)

        self.read_markers([-5, 1, 1])
        self.read_markers([600])
        data = self.read_block()  # 604

        self.read_markers([-6, 1, 1])
        self.read_markers([188])
        data = self.read_block()

        self.read_markers([14])
        data = self.read_block()
        self.read_markers([16])
        data = self.read_block()
        self.read_markers([18])
        data = self.read_block()
        self.read_markers([84])
        data = self.read_block()
        self.read_markers([6])
        data = self.read_block()

        self.read_markers([-7, 1, 1])
        self.read_markers([342])
        data = self.read_block()

        self.read_markers([-8, 1, 1])
        data = self.read_block()
        data = self.read_block()
        while 1:
            n = self.get_nmarkers(1, rewind=True)[0]
            if n not in [2, 4, 6, 8]:
                #print "n =", n
                break
            n = self.get_nmarkers(1, rewind=False)[0]
            #print n
            data = self.read_block()


        i = -9
        while i != -13:
            n = self.get_nmarkers(1, rewind=True)[0]

            self.read_markers([i, 1, 1])
            while 1:
                n = self.get_nmarkers(1, rewind=True)[0]
                if n not in [2, 4, 6, 8]:
                    #print "n =", n
                    break
                n = self.get_nmarkers(1, rewind=False)[0]
                #print n
                data = self.read_block()
            i -= 1

        #print "n=%s" % (self.n)
        #strings, ints, floats = self.show(100)
        pass

    def _skip_pcompts(self):
        """
        Reads the PCOMPTS table (poorly).
        The PCOMPTS table stores information about the PCOMP cards???
        """
        self.log.debug("table_name = %r" % self.table_name)
        table_name = self.read_table_name(rewind=False)

        self.read_markers([-1])
        data = self._skip_record()

        self.read_markers([-2, 1, 0])
        data = self._skip_record()
        #table_name, = unpack(b'8s', data)
        #print "table_name = %r" % table_name

        self.read_markers([-3, 1, 0])
        markers = self.get_nmarkers(1, rewind=True)
        if markers != [-4]:
            data = self._skip_record()

        self.read_markers([-4, 1, 0])
        markers = self.get_nmarkers(1, rewind=True)
        if markers != [0]:
            data = self._skip_record()
        else:
            self.read_markers([0])
            return

        self.read_markers([-5, 1, 0])
        data = self._skip_record()

        self.read_markers([-6, 1, 0])
        self.read_markers([0])

    def _read_pcompts(self):
        """
        Reads the PCOMPTS table (poorly).
        The PCOMPTS table stores information about the PCOMP cards???
        """
        self._skip_pcompts()
        return
        if self.read_mode == 1:
            return
        self.log.debug("table_name = %r" % self.table_name)
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
        """Reads the next OP2 table name (e.g. OUG1, OES1X1)"""
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
        """
        Skips the majority of geometry/result tables as they follow a very standard format.
        Other tables don't follow this format.
        """
        self.table_name = self.read_table_name(rewind=False)
        if self.debug:
            self.binary_debug.write('skipping table...%r\n' % self.table_name)
        self.read_markers([-1])
        data = self._skip_record()
        self.read_markers([-2, 1, 0])
        data = self._skip_record()
        self._skip_subtables()

    def _read_omm2(self):
        """
        :param self:    the OP2 object pointer
        """
        self.log.debug("table_name = %r" % self.table_name)
        self.table_name = self.read_table_name(rewind=False)
        self.read_markers([-1])
        data = self._read_record()

        self.read_markers([-2, 1, 0])
        data = self._read_record()
        if len(data) == 28:
            subtable_name, month, day, year, zero, one = unpack(b'8s5i', data)
            if self.debug:
                self.binary_debug.write('  recordi = [%r, %i, %i, %i, %i, %i]\n'  % (subtable_name, month, day, year, zero, one))
                self.binary_debug.write('  subtable_name=%r\n' % subtable_name)
            self._print_month(month, day, year, zero, one)
        else:
            raise NotImplementedError(self.show_data(data))
        self._read_subtables()

    def _read_fol(self):
        """
        :param self:    the OP2 object pointer
        """
        self.log.debug("table_name = %r" % self.table_name)
        self.table_name = self.read_table_name(rewind=False)
        self.read_markers([-1])
        data = self._read_record()

        self.read_markers([-2, 1, 0])
        data = self._read_record()
        if len(data) == 12:
            subtable_name, double = unpack(b'8sf', data)
            if self.debug:
                self.binary_debug.write('  recordi = [%r, %f]\n'  % (subtable_name, double))
                self.binary_debug.write('  subtable_name=%r\n' % subtable_name)
        else:
            strings, ints, floats = self.show_data(data)
            msg = 'len(data) = %i\n' % len(data)
            msg += 'strings  = %r\n' % strings
            msg += 'ints     = %r\n' % str(ints)
            msg += 'floats   = %r' % str(floats)
            raise NotImplementedError(msg)
        self._read_subtables()

    def _read_gpl(self):
        """
        :param self:    the OP2 object pointer
        """
        self.table_name = self.read_table_name(rewind=False)
        self.log.debug('table_name = %r' % self.table_name)
        if self.debug:
            self.binary_debug.write('_read_geom_table - %s\n' % self.table_name)
        self.read_markers([-1])
        if self.debug:
            self.binary_debug.write('---markers = [-1]---\n')
        data = self._read_record()

        markers = self.get_nmarkers(1, rewind=True)
        self.binary_debug.write('---marker0 = %s---\n' % markers)
        n = -2
        while markers[0] != 0:
            self.read_markers([n, 1, 0])
            if self.debug:
                self.binary_debug.write('---markers = [%i, 1, 0]---\n' % n)

            markers = self.get_nmarkers(1, rewind=True)
            if markers[0] == 0:
                markers = self.get_nmarkers(1, rewind=False)
                break
            data = self._read_record()
            #self.show_data(data, 'i')
            n -= 1
            markers = self.get_nmarkers(1, rewind=True)

    def _read_meff(self):
        """
        :param self:    the OP2 object pointer
        """
        raise NotImplementedError(self.table_name)
        self.table_name = self.read_table_name(rewind=False)
        self.log.debug('table_name = %r' % self.table_name)
        if self.debug:
            self.binary_debug.write('_read_geom_table - %s\n' % self.table_name)
        self.read_markers([-1])
        if self.debug:
            self.binary_debug.write('---markers = [-1]---\n')
        data = self._read_record()

        markers = self.get_nmarkers(1, rewind=True)
        self.binary_debug.write('---marker0 = %s---\n' % markers)
        self.read_markers([-2, 1, 0])
        data = self._read_record()

        for n in [-3, -4, -5, -6, -7, -8, ]:
            self.read_markers([n, 1, 1])
            markers = self.get_nmarkers(1, rewind=False)
            print('markers =', markers)
            self.f.read(markers[0]*4+12)
            self.n += markers[0]*4+12

        n = -9
        self.read_markers([n, 1, 0, 0])
        #data = self._read_record()
        self.show(50)

    def get_marker_n(self, n):
        markers = []
        s = Struct('3i')
        for i in xrange(n):
            block = self.f.read(12)
            marker = s.unpack(block)
            print('markers =', marker)
            markers.append(marker)
        print("markers =", markers)
        return markers

    def _read_geom_table(self):
        """
        Reads a geometry table
        :param self:    the OP2 object pointer
        """
        self.table_name = self.read_table_name(rewind=False)
        if self.debug:
            self.binary_debug.write('_read_geom_table - %s\n' % self.table_name)
        self.read_markers([-1])
        data = self._read_record()

        self.read_markers([-2, 1, 0])
        data = self._read_record()
        if len(data) == 8:
            table_name, = unpack(b'8s', data)
        else:
            strings, ints, floats = self.show_data(data)
            msg = 'len(data) = %i\n' % len(data)
            msg += 'strings  = %r\n' % strings
            msg += 'ints     = %r\n' % str(ints)
            msg += 'floats   = %r' % str(floats)
            raise NotImplementedError(msg)
        self._read_subtables()

    def _read_sdf(self):
        """
        :param self:    the OP2 object pointer
        """
        self.log.debug("table_name = %r" % self.table_name)
        self.table_name = self.read_table_name(rewind=False)
        self.read_markers([-1])
        data = self._read_record()

        self.read_markers([-2, 1, 0])
        data = self._read_record()
        if len(data) == 16:
            subtable_name, dummyA, dummyB = unpack(b'8sii', data)
            if self.debug:
                self.binary_debug.write('  recordi = [%r, %i, %i]\n'  % (subtable_name, dummyA, dummyB))
                self.binary_debug.write('  subtable_name=%r\n' % subtable_name)
                assert dummyA == 170, dummyA
                assert dummyB == 170, dummyB
        else:
            strings, ints, floats = self.show_data(data)
            msg = 'len(data) = %i\n' % len(data)
            msg += 'strings  = %r\n' % strings
            msg += 'ints     = %r\n' % str(ints)
            msg += 'floats   = %r' % str(floats)
            raise NotImplementedError(msg)

        self.read_markers([-3, 1, 1])

        markers0 = self.get_nmarkers(1, rewind=False)
        record = self.read_block()

        #data = self._read_record()
        self.read_markers([-4, 1, 0, 0])
        #self._read_subtables()

    def _read_results_table(self):
        """
        Reads a results table
        """
        if self.debug:
            self.binary_debug.write('read_results_table - %s\n' % self.table_name)
        self.table_name = self.read_table_name(rewind=False)
        self.read_markers([-1])
        if self.debug:
            self.binary_debug.write('---markers = [-1]---\n')
            #self.binary_debug.write('marker = [4, -1, 4]\n')
        data = self._read_record()

        self.read_markers([-2, 1, 0])
        if self.debug:
            self.binary_debug.write('---markers = [-2, 1, 0]---\n')
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
            strings, ints, floats = self.show_data(data)
            msg = 'len(data) = %i\n' % len(data)
            msg += 'strings  = %r\n' % strings
            msg += 'ints     = %r\n' % str(ints)
            msg += 'floats   = %r' % str(floats)
            raise NotImplementedError(msg)
        self._read_subtables()

    def _print_month(self, month, day, year, zero, one):
        """
        Creates the self.date attribute from the 2-digit year.

        :param month: the month (integer <= 12)
        :param day:  the day (integer <= 31)
        :param year: the day (integer <= 99)
        :param zero: a dummy integer (???)
        :param one:  a dummy integer (???)
        """
        month, day, year = self._set_op2_date(month, day, year)

        #self.log.debug("%s/%s/%4i zero=%s one=%s" % (month, day, year, zero, one))
        if self.debug:
            self.binary_debug.write('  [subtable_name, month=%i, day=%i, year=%i, zero=%i, one=%i]\n\n' % (month, day, year, zero, one))
        #assert zero == 0  # is this the RTABLE indicator???
        assert one == 1

    def finish(self):
        """
        Clears out the data members contained within the self.words variable.
        This prevents mixups when working on the next table, but otherwise
        has no effect.
        """
        for word in self.words:
            if word != '???' and hasattr(self, word):
                if word not in ['Title', 'reference_point']:
                    delattr(self, word)
        self.obj = None

    def get_op2_stats(self):
        """
        Gets info about the contents of the different attributes of the
        OP2 class.
        """
        table_types = self.get_table_types()
        msg = []
        for table_type in table_types:
            table = getattr(self, table_type)
            for isubcase, subcase in sorted(table.iteritems()):
                if hasattr(subcase, 'get_stats'):
                    msg.append('op2.%s[%s]\n' % (table_type, isubcase))
                    msg.extend(subcase.get_stats())
                    msg.append('\n')
                else:
                    msg.append('skipping %s op2.%s[%s]\n\n' % (subcase.__class__.__name__, table_type, isubcase))
                    #raise RuntimeError('skipping %s op2.%s[%s]\n\n' % (subcase.__class__.__name__, table_type, isubcase))
        return ''.join(msg)

if __name__ == '__main__':  # pragma: no conver
    from pickle import dumps, dump, load, loads
    txt_filename = 'solid_shell_bar.txt'
    f = open(txt_filename, 'wb')
    op2_filename = 'solid_shell_bar.op2'
    op2 = OP2()
    op2.read_op2(op2_filename)
    print op2.displacements[1]
    dump(op2, f)
    f.close()

    f = open(txt_filename, 'r')
    op2 = load(f)
    f.close()
    print op2.displacements[1]


    #import sys
    #op2_filename = sys.argv[1]

    #o = OP2(op2_filename)
    #o.read_op2(op2_filename)
    #(model, ext) = os.path.splitext(op2_filename)
    #f06_outname = model + '.test_op2.f06'
    #o.write_f06(f06_outname)
