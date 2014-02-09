#pylint: disable=C0301,C0103,W0613,C0111,W0612,R0913
"""
Defines the OP2 class.
"""
import os
from struct import unpack, Struct
#from struct import error as StructError

from numpy import array

from pyNastran.f06.f06 import FatalError
from pyNastran.f06.tables.grid_point_weight import GridPointWeight

from pyNastran.op2.dev_explicit.tables.geom1 import GEOM1
from pyNastran.op2.dev_explicit.tables.geom2 import GEOM2
from pyNastran.op2.dev_explicit.tables.geom3 import GEOM3
from pyNastran.op2.dev_explicit.tables.geom4 import GEOM4

from pyNastran.op2.dev_explicit.tables.ept import EPT
from pyNastran.op2.dev_explicit.tables.mpt import MPT

from pyNastran.op2.dev_explicit.tables.oef import OEF
from pyNastran.op2.dev_explicit.tables.oes import OES
from pyNastran.op2.dev_explicit.tables.ogs import OGS

from pyNastran.op2.dev_explicit.tables.opg import OPG
from pyNastran.op2.dev_explicit.tables.oqg import OQG
from pyNastran.op2.dev_explicit.tables.oug import OUG
from pyNastran.op2.dev_explicit.tables.ogpwg import OGPWG
from pyNastran.op2.dev_explicit.fortran_format import FortranFormat

from pyNastran.bdf.bdf import BDF

from pyNastran.utils import is_binary


class OP2(BDF, GEOM1, GEOM2, GEOM3, GEOM4, EPT, MPT,
          OEF, OES, OGS, OPG, OQG, OUG, OGPWG, FortranFormat):
    """
    Defines an interface for the Nastran OP2 file.
    """
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

    def set_result_types(self, result_types):
        """
        :param result_types: allows a user to reduce the amount of data that's extracted

        Example
        =======
        op2 = OP2(op2_filename)
        op2.set_result_types(['displacements', 'solidStress'])
        """
        table_mapper = {  # very incomplete list
            "geometry" : [],  # includes materials/properties
            "loads_constraints" : [],

            'displacements' : ['OUG1', 'OUGV1', 'BOUGV1', 'OUPV1', ],
            'velocities'    : ['OUG1', 'OUGV1', 'BOUGV1', 'OUPV1', ],
            'accelerations' : ['OUG1', 'OUGV1', 'BOUGV1', 'OUPV1', ],
            'eigenvectors'  : ['OUG1', 'OUGV1', 'BOUGV1', 'OUPV1', ],
            'temperatures'  : ['OUG1', 'OUGV1', 'BOUGV1', 'OUPV1', ],

            # applied loads
            'appliedLoads' : ['OPG1', 'OPGV1', 'OPNL1', ],
            'loadVectors' : ['OPG1', 'OPGV1', 'OPNL1', ],
            'thermalLoadVectors' : ['OPG1', 'OPGV1', 'OPNL1', ],

            'barStress'    : ['OES1X1', 'OES1', 'OES1X', 'OESNLXR', 'OESNLXD', 'OESNLBR', 'OESTRCP', 'OESNL1X', 'OESRT', ],
            'rodStress'    : ['OES1X1', 'OES1', 'OES1X', 'OESNLXR', 'OESNLXD', 'OESNLBR', 'OESTRCP', 'OESNL1X', 'OESRT', ],
            'beamStress'   : ['OES1X1', 'OES1', 'OES1X', 'OESNLXR', 'OESNLXD', 'OESNLBR', 'OESTRCP', 'OESNL1X', 'OESRT', ],
            'solidStress'  : ['OES1X1', 'OES1', 'OES1X', 'OESNLXR', 'OESNLXD', 'OESNLBR', 'OESTRCP', 'OESNL1X', 'OESRT', ],
            'plateStress'  : ['OES1X1', 'OES1', 'OES1X', 'OESNLXR', 'OESNLXD', 'OESNLBR', 'OESTRCP', 'OESNL1X', 'OESRT', ],
            'shearStress'  : ['OES1X1', 'OES1', 'OES1X', 'OESNLXR', 'OESNLXD', 'OESNLBR', 'OESTRCP', 'OESNL1X', 'OESRT', ],
            'compositePlateStress': ['OES1C', 'OESCP'],

            'barStrain'    : ['OSTR1X'],
            'rodStrain'    : ['OSTR1X'],
            'beamStrain'   : ['OSTR1X'],
            'solidStrain'  : ['OSTR1X'],
            'plateStrain'  : ['OSTR1X'],
            'shearStrain'  : ['OSTR1X'],
            'compositePlateStrain': ['OSTR1C'],

            'barForces'    : ['OEFIT', 'OEF1X', 'OEF1', 'DOEF1'],
            'beamForces'   : ['OEFIT', 'OEF1X', 'OEF1', 'DOEF1'],
            'rodForces'    : ['OEFIT', 'OEF1X', 'OEF1', 'DOEF1'],
            'plateForces'  : ['OEFIT', 'OEF1X', 'OEF1', 'DOEF1'],
            'bar100Forces' : ['OEFIT', 'OEF1X', 'OEF1', 'DOEF1'],
            #'solidForces'  : ['OEFIT', 'OEF1X', 'OEF1', 'DOEF1'],

            'grid_point_stresses' : ['OGS1'],

            'grid_point_weight': ['OGPWG'],

            'spcForces' : ['OQG1', 'OQGV1'],
            'mpcForces' : ['OQMG1'],
        }

        table_names_set = set([])
        for result_type in result_types:
            if not hasattr(self, result_type):
                raise RuntimeError('result_type=%r is invalid')
            if not result_type in table_mapper:
                raise NotImplementedError('result_type=%r is not supported by set_result_types')

            itables = table_mapper[result_type]
            for itable in itables:
                table_names_set.add(itable)
        self.set_tables_to_read(table_names_set)  # is this the right function name???

    def set_tables_to_read(self, table_names):
        """
        :param table_names: list/set of tables as strings
        """
        self.tables_to_read = table_names

    def __init__(self, op2_filename, make_geom=False, save_skipped_cards=False, debug=False, log=None):
        self.save_skipped_cards = save_skipped_cards
        if op2_filename:
            self.op2_filename = op2_filename
        self.make_geom = make_geom

        #self.tables_to_read = []
        BDF.__init__(self, debug=debug, log=log)

        GEOM1.__init__(self)
        GEOM2.__init__(self)
        GEOM3.__init__(self)
        GEOM4.__init__(self)

        EPT.__init__(self)
        MPT.__init__(self)

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
            'OEFIT' : [self._read_oef1_3, self._read_oef1_4],
            'OEF1X' : [self._read_oef1_3, self._read_oef1_4],
            'OEF1'  : [self._read_oef1_3, self._read_oef1_4],
            'DOEF1' : [self._read_oef1_3, self._read_oef1_4],
            #=======================
            # OQG
            # spc forces
            'OQG1'  : [self._read_oqg1_3, self._read_oqg1_4],
            'OQGV1' : [self._read_oqg1_3, self._read_oqg1_4],
            # mpc forces
            'OQMG1' : [self._read_oqg1_3, self._read_oqg1_4],
            #=======================
            # OPG
            # applied loads
            'OPG1'  : [self.read_opg1_3, self.read_opg1_4],
            'OPGV1' : [self.read_opg1_3, self.read_opg1_4],
            'OPNL1' : [self.read_opg1_3, self.read_opg1_4],

            # OGPFB1
            # grid point forces
            'OGPFB1' : [self.read_ogpf1_3, self.read_ogpf1_4],

            # ONR
            # grid point forces
            'ONRGY1' : [self.read_onr1_3, self.read_onr1_4],
            #=======================
            # OES
            # stress
            'OES1X1'  : [self._read_oes1_3, self._read_oes1_4],
            'OES1'    : [self._read_oes1_3, self._read_oes1_4],
            'OES1X'   : [self._read_oes1_3, self._read_oes1_4],
            'OES1C'   : [self._read_oes1_3, self._read_oes1_4],
            'OESCP'   : [self._read_oes1_3, self._read_oes1_4],
            'OESNLXR' : [self._read_oes1_3, self._read_oes1_4],
            'OESNLXD' : [self._read_oes1_3, self._read_oes1_4],
            'OESNLBR' : [self._read_oes1_3, self._read_oes1_4],
            'OESTRCP' : [self._read_oes1_3, self._read_oes1_4],
            'OESNL1X' : [self._read_oes1_3, self._read_oes1_4],
            'OESRT'   : [self._read_oes1_3, self._read_oes1_4],
            # strain
            'OSTR1X'  : [self._read_oes1_3, self._read_oes1_4],
            'OSTR1C'  : [self._read_oes1_3, self._read_oes1_4],

            #=======================
            # OUG
            # displacement/velocity/acceleration/eigenvector
            'OUG1'   : [self._read_oug1_3, self._read_oug1_4],
            'OUGV1'  : [self._read_oug1_3, self._read_oug1_4],
            'BOUGV1' : [self._read_oug1_3, self._read_oug1_4],

            'OUPV1'  : [self._read_oug1_3, self._read_oug1_4],

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

            # eigenvalues
            'BLAMA': [self._read_complex_eigenvalue_3, self._read_complex_eigenvalue_4],
            'LAMA': [self._read_real_eigenvalue_3, self._read_real_eigenvalue_4],

            # geometry
            'GEOM1': [self._not_available, self._read_geom1_4],
            'GEOM2': [self._not_available, self._read_geom2_4],
            'GEOM3': [self._not_available, self._read_geom3_4],
            'GEOM4': [self._not_available, self._read_geom4_4],

            # superelements
            'GEOM1S': [self._not_available, self._read_geom1_4],
            'GEOM2S': [self._not_available, self._read_geom2_4],
            'GEOM3S': [self._not_available, self._read_geom3_4],
            'GEOM4S': [self._not_available, self._read_geom4_4],

            'GEOM1N': [self._not_available, self._read_geom1_4],
            #'GEOM2N': [self._not_available, self._read_geom2_4],
            #'GEOM3N': [self._not_available, self._read_geom3_4],
            #'GEOM4N': [self._not_available, self._read_geom4_4],

            'GEOM1OLD': [self._not_available, self._read_geom1_4],
            'GEOM2OLD': [self._not_available, self._read_geom2_4],
            #'GEOM3OLD': [self._not_available, self._read_geom3_4],
            'GEOM4OLD': [self._not_available, self._read_geom4_4],

            'EPT' : [self._not_available, self._read_ept_4],
            'EPTS': [self._not_available, self._read_ept_4],
            'EPTOLD' : [self._not_available, self._read_ept_4],

            'MPT' : [self._not_available, self._read_mpt_4],
            'MPTS': [self._not_available, self._read_mpt_4],
        }
        self.make_geom = False

    def read_onr1_3(self, data):  # TODO: this is wrong...
        self.read_opg1_3(data)

    def read_onr1_4(self, data):
        if self.table_code == 18:  # element strain energy
            assert self.table_name in ['ONRGY1'], 'table_name=%s table_code=%s' % (self.table_name, self.table_code)
            n = self._read_element_strain_energy(data)
        else:
            raise NotImplementedError(self.table_code)
        return n

    def read_ogpf1_3(self, data):
        self.read_opg1_3(data)  # TODO: this is wrong...

    def read_ogpf1_4(self, data):
        if self.table_code == 19:  # grid point force balance
            assert self.table_name in ['OGPFB1'], 'table_name=%s table_code=%s' % (self.table_name, self.table_code)
            n = self._read_grid_point_forces(data)
        else:
            raise NotImplementedError(self.table_code)
        return n

    def _read_element_strain_energy(self, data):
        """
        table_code = 19
        """
        from pyNastran.op2.tables.oee_energy.oee_objects import StrainEnergyObject
        dt = self.nonlinear_factor
        n = 0
        if self.thermal == 0:
            result_name = 'strainEnergy'
            if self.num_wide == 4:
                self.create_transient_object(self.strainEnergy, StrainEnergyObject)
                s = Struct(b'i3f')

                ntotal = 16
                nnodes = len(data) // ntotal
                for i in xrange(nnodes):
                    eData = data[n:n+ntotal]

                    out = s.unpack(eData)
                    (eid_device, energy, percent, density) = out
                    eid = (eid_device - self.device_code) // 10
                    #print "eType=%s" %(eType)

                    dataIn = [eid, energy, percent, density]
                    #print "%s" %(self.get_element_type(self.element_type)),dataIn
                    self.obj.add(dt, dataIn)
                    n += ntotal
            elif self.num_wide == 5:
                self.create_transient_object(self.strainEnergy, StrainEnergyObject)  # why is this not different?
                s = Struct(b'i8s3f')
                ntotal = 20
                nnodes = len(data) // ntotal
                for i in xrange(nnodes):
                    eData = self.data[n:n+20]
                    out = unpack(format1, eData)
                    (word, energy, percent, density) = out
                    #print "out = ",out
                    word = word.strip()
                    #print "eType=%s" %(eType)

                    dataIn = [word, energy, percent, density]
                    #print "%s" %(self.get_element_type(self.element_type)),dataIn
                    #eid = self.obj.add_new_eid(out)
                    self.obj.add(dt, dataIn)
                    n += ntotal
            else:
                raise NotImplementedError('num_wide = %s' % (self.num_wide))

            #complex_obj = complexGridPointForcesObject

            #self._read_table(data, storage_obj, real_obj, complex_obj, 'node')
        #elif self.thermal == 1:
            #result_name = 'thermalLoadVectors'
            #storage_obj = self.thermalLoadVectors
            #real_obj = ThermalLoadVectorObject
            #complex_obj = None
            #self._read_table(data, storage_obj, real_obj, complex_obj, 'node')
        else:
            raise NotImplementedError(self.thermal)
        return n

    def _read_grid_point_forces(self, data):
        """
        table_code = 19
        """
        from pyNastran.op2.tables.ogf_gridPointForces.ogf_Objects import gridPointForcesObject, complexGridPointForcesObject
        dt = self.nonlinear_factor
        n = 0
        if self.thermal == 0:
            result_name = 'gridPointForces'
            if self.num_wide == 10:
                self.create_transient_object(self.gridPointForces, gridPointForcesObject)
                s = Struct(b'ii8s6f')
                ntotal = 40
                nnodes = len(data) // ntotal
                for i in xrange(nnodes):
                    eData = data[n:n+ntotal]
                    out = s.unpack(eData)
                    (eKey, eid, elemName, f1, f2, f3, m1, m2, m3) = out
                    ekey = (eKey - self.device_code) // 10
                    elemName = elemName.strip()
                    #data = (eid,elemName,f1,f2,f3,m1,m2,m3)
                    self.obj.add(dt, eKey, eid, elemName, f1, f2, f3, m1, m2, m3)
                    #print "eid/dt/freq=%s eid=%-6s eName=%-8s f1=%g f2=%g f3=%g m1=%g m2=%g m3=%g" %(ekey,eid,elemName,f1,f2,f3,m1,m2,m3)
                    n += ntotal
            else:
                raise NotImplementedError('num_wide = %s' % (self.num_wide))

            #complex_obj = complexGridPointForcesObject

            #self._read_table(data, storage_obj, real_obj, complex_obj, 'node')
        #elif self.thermal == 1:
            #result_name = 'thermalLoadVectors'
            #storage_obj = self.thermalLoadVectors
            #real_obj = ThermalLoadVectorObject
            #complex_obj = None
            #self._read_table(data, storage_obj, real_obj, complex_obj, 'node')
        else:
            raise NotImplementedError(self.thermal)
        return n

    def _not_available(self, data):
        raise RuntimeError('this should never be called...table_name=%r' % self.table_name)

    def _read_complex_eigenvalue_3(self, data):
        self.show_data(data)
        aaa

    def _read_complex_eigenvalue_4(self, data):
        self.show_data(data)
        bbb

    def _read_real_eigenvalue_3(self, data):
        #self.show_data(data)

        (three) = self.parse_approach_code(data)

        self.add_data_parameter(data, 'seven', 'i', 10, False)  # seven
        ## residual vector augmentation flag
        self.add_data_parameter(data, 'resFlag', 'i', 11, False)
        ## fluid modes Flag
        self.add_data_parameter(data, 'fldFlag', 'i', 12, False)

        #print self.data_code
        #self.add_data_parameter(data,'format_code',  'i',9,False)   ## format code
        #self.add_data_parameter(data,'num_wide',     'i',10,False)  ## number of words per entry in record; .. note:: is this needed for this table ???

        #if self.analysis_code==2: # sort2
        #    self.lsdvmn = self.get_values(data,'i',5)

        #print "*isubcase=%s"%(self.isubcase)
        #print "analysis_code=%s table_code=%s thermal=%s" %(self.analysis_code,self.table_code,self.thermal)

        #self.print_block(data)
        self._read_title(data)


    def _read_real_eigenvalue_4(self, data):
        #self.show_data(data)
        nModes = len(data) // 28
        from pyNastran.op2.tables.lama_eigenvalues.lama_objects import (
            RealEigenvalues, ComplexEigenvalues)

        n = 0
        ntotal = 28
        lama = RealEigenvalues(self.isubcase)
        self.eigenvalues[self.isubcase] = lama
        s = Struct('ii5f')
        for i in xrange(nModes):
            edata = data[n:n+28]
            out = s.unpack(edata)
            if self.debug4():
                self.binary_debug.write('  eigenvalue%s - %s\n' % (i, str(out)))
            #(iMode,order,eigen,omega,freq,mass,stiff) = out
            #(modeNum,extractOrder,eigenvalue,radian,cycle,genM,genK) = line
            #print out
            lama.addF06Line(out)
            n += ntotal
        return n

    def readFake(self, data, n):
        return n

    def _read_geom1_4(self, data):
        self._read_geom_4(self._geom1_map, data)

    def _read_geom2_4(self, data):
        self._read_geom_4(self._geom2_map, data)

    def _read_geom3_4(self, data):
        self._read_geom_4(self._geom3_map, data)

    def _read_geom4_4(self, data):
        self._read_geom_4(self._geom4_map, data)

    def _read_ept_4(self, data):
        self._read_geom_4(self._ept_map, data)
    def _read_mpt_4(self, data):
        self._read_geom_4(self._mpt_map, data)


    def _read_geom_4(self, mapper, data):
        if not self.make_geom:
            return
        n = 0
        keys = unpack('3i', data[n:n+12])
        n += 12
        if len(data) == 12:
            pass
        elif keys in mapper:
            func = mapper[keys]
            if isinstance(func, list):
                name, func = func
                print "found keys=(%5s,%4s,%4s) name=%s" % (keys[0], keys[1], keys[2], name)
            else:
                print "found keys=(%5s,%4s,%4s)" % (keys[0], keys[1], keys[2])
            n = func(data, n)  # gets all the grid/mat cards
        else:
            raise NotImplementedError('keys=%s not found' % str(keys))
        #assert n == len(data), 'n=%s len(data)=%s' % (n, len(data))
        #self.show_data(data[n:])

    def read_op2(self, op2_filename=None):
        """
        Starts the OP2 file reading
        :param op2_filename: if a string is set, the filename specified in the
            __init__ is overwritten.

            init=None,  op2_filename=None   -> a dialog is popped up  (not implemented; crash)
            init=fname, op2_filename=fname  -> fname is used
            init=fname, op2_filename=fname2 -> fname2 is used
            init=None,  op2_filename=fname  -> fname is used
        """
        if op2_filename:
            self.op2_filename = op2_filename
        if self.save_skipped_cards:
            self.skippedCardsFile = open('skippedCards.out', 'a')
        self.n = 0
        self.table_name = None
        if not is_binary(self.op2_filename):
            if os.path.getsize(self.op2_filename) == 0:
                raise IOError('op2_filename=%r is empty.' % self.op2_filename)
            raise IOError('op2_filename=%r is not a binary OP2.' % self.op2_filename)
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

        table_names = self._read_tables(table_name)
        if self.debug:
            self.binary_debug.write('-' * 80 + '\n')
            self.binary_debug.write('f.tell()=%s\ndone...\n' % self.f.tell())
        if self.save_skipped_cards:
            self.skippedCardsFile.close()
        return table_names

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
                                   # grid point weight
                                   'OGPWG',

                                   # other
                                   'CONTACT', 'VIEWTB',
                                   'KDICT', 'MONITOR',
                                  ]:
                    self._read_geom_table()  # DIT (agard)
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
                                    # displacement/velocity/acceleration/eigenvector
                                    'OUG1', 'OUGV1', 'BOUGV1', 'OUPV1',
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
            print "i=%s n=%s" % (i, self.n)

        self.read_markers([4])
        data = self.read_block()

        for i in xrange(7):
            self.read_markers([2])
            data = self.read_block()
            print "i=%s n=%s" % (i, self.n)


        self.read_markers([-4, 1, 1])
        for i in xrange(170//10):
            self.read_markers([2])
            data = self.read_block()
            print "i=%s n=%s" % (i, self.n)

        self.read_markers([4])
        data = self.read_block()

        for i in xrange(7):
            self.read_markers([2])
            data = self.read_block()
            print "i=%s n=%s" % (i, self.n)

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
                print "n =", n
                break
            n = self.get_nmarkers(1, rewind=False)[0]
            print n
            data = self.read_block()


        i = -9
        while i != -13:
            n = self.get_nmarkers(1, rewind=True)[0]

            self.read_markers([i, 1, 1])
            while 1:
                n = self.get_nmarkers(1, rewind=True)[0]
                if n not in [2, 4, 6, 8]:
                    print "n =", n
                    break
                n = self.get_nmarkers(1, rewind=False)[0]
                print n
                data = self.read_block()
            i -= 1

        print "n=%s" % (self.n)
        #strings, ints, floats = self.show(100)

        pass

    def _read_pcompts(self):
        """
        Reads the PCOMPTS table (poorly).
        The PCOMPTS table stores information about the PCOMP cards???
        """
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
        if self.debug:
            self.binary_debug.write('skipping table...\n')
        self.table_name = self.read_table_name(rewind=False)
        self.read_markers([-1])
        data = self._skip_record()

        self.read_markers([-2, 1, 0])
        data = self._skip_record()
        self._skip_subtables()

    def _read_omm2(self):
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

    def _read_geom_table(self):
        """
        Reads a geometry table
        """
        self.table_name = self.read_table_name(rewind=False)
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
        self.show(100)

        #sys.exit()
        #self._read_subtables()

    def _read_results_table(self):
        """
        Reads a results table
        """
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
        table_types = [
            ## OUG - displacement
            'displacements',
            'displacementsPSD',
            'displacementsATO',
            'displacementsRMS',
            'displacementsCRM',
            'displacementsNO',
            'scaledDisplacements',

            ## OUG - temperatures
            'temperatures',

            ## OUG - eigenvectors
            'eigenvectors',

            ## OUG - velocity
            'velocities',

            ## OUG - acceleration
            'accelerations',

            # OQG - spc/mpc forces
            'spcForces',
            'mpcForces',
            'thermalGradientAndFlux',

            ## OGF - grid point forces
            'gridPointForces',

            ## OPG - summation of loads for each element
            'loadVectors',
            'thermalLoadVectors',
            'appliedLoads',
            'forceVectors',

            # OES - tCode=5 thermal=0 s_code=0,1 (stress/strain)
            ## OES - CELAS1/CELAS2/CELAS3/CELAS4 stress
            'celasStress',
            ## OES - CELAS1/CELAS2/CELAS3/CELAS4 strain
            'celasStrain',

            ## OES - isotropic CROD/CONROD/CTUBE stress
            'rodStress',
            'conrodStress',
            'ctubeStress',

            ## OES - isotropic CROD/CONROD/CTUBE strain
            'rodStrain',
            'conrodStrain',
            'ctubeStrain',

            ## OES - isotropic CBAR stress
            'barStress',
            ## OES - isotropic CBAR strain
            'barStrain',
            ## OES - isotropic CBEAM stress
            'beamStress',
            ## OES - isotropic CBEAM strain
            'beamStrain',

            ## OES - isotropic CTRIA3/CQUAD4 stress
            'plateStress',
            # OES - isotropic CTRIA3/CQUAD4 strain
            'plateStrain',
            ## OES - isotropic CTETRA/CHEXA/CPENTA stress
            'solidStress',
            ## OES - isotropic CTETRA/CHEXA/CPENTA strain
            'solidStrain',

            ## OES - CSHEAR stress
            'shearStress',
            ## OES - CSHEAR strain
            'shearStrain',
            ## OES - CEALS1 224, CELAS3 225
            'nonlinearSpringStress',
            ## OES - GAPNL 86
            'nonlinearGapStress',
            ## OES - CBUSH 226
            'nolinearBushStress',
        ]

        table_types += [
            # LAMA
            'eigenvalues',

            # OEF - Forces - tCode=4 thermal=0
            'rodForces',
            'conrodForces',
            'ctubeForces',

            'barForces',
            'bar100Forces',
            'beamForces',
            'bendForces',
            'bushForces',
            'coneAxForces',
            'damperForces',
            'gapForces',
            'plateForces',
            'plateForces2',
            'shearForces',
            'solidPressureForces',
            'springForces',
            'viscForces',

            'force_VU',
            'force_VU_2D',

            #OEF - Fluxes - tCode=4 thermal=1
            'thermalLoad_CONV',
            'thermalLoad_CHBDY',
            'thermalLoad_1D',
            'thermalLoad_2D_3D',
            'thermalLoad_VU',
            'thermalLoad_VU_3D',
            'thermalLoad_VUBeam',
            #self.temperatureForces
        ]
        table_types += [
            ## OES - CTRIAX6
            'ctriaxStress',
            'ctriaxStrain',

            'bushStress',
            'bushStrain',
            'bush1dStressStrain',

            ## OES - nonlinear CROD/CONROD/CTUBE stress
            'nonlinearRodStress',
            'nonlinearRodStrain',

            ## OESNLXR - CTRIA3/CQUAD4 stress
            'nonlinearPlateStress',
            'nonlinearPlateStrain',
            'hyperelasticPlateStress',
            'hyperelasticPlateStrain',

            ## OES - composite CTRIA3/CQUAD4 stress
            'compositePlateStress',
            'compositePlateStrain',

            ## OGS1 - grid point stresses
            'gridPointStresses',        # tCode=26
            'gridPointVolumeStresses',  # tCode=27

            ## OEE - strain energy density
            'strainEnergy',  # tCode=18
        ]
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

    def print_results(self):
        """
        Prints an ASCII summary of the OP2.
        ..warning:: Will take lots of memory for large OP2 files.
        """
        results = [
            # OUG - Displacements/Velocity/Acceleration/Temperature
            self.displacements, self.displacementsPSD,
            self.displacementsATO,
            self.temperatures,
            self.eigenvalues,
            self.eigenvectors,
            self.velocities,
            self.accelerations,

            # OEF - Applied Forces/Temperatures - ???

            # OQG1 - SPC/MPC Forces
            self.spcForces, self.mpcForces,

            # OGF - Grid Point Forces
            self.gridPointForces,

            # OPG - Applied Force/Moment
            self.appliedLoads,
            self.loadVectors, self.thermalLoadVectors,
            self.forceVectors,

            # OES - Stress/Strain
            self.celasStress, self.celasStrain,
            self.rodStress, self.rodStrain,
            self.conrodStress, self.conrodStrain,
            self.nonlinearRodStress, self.nonlinearRodStrain,

            self.barStress, self.barStrain,
            self.beamStress, self.beamStrain,
            self.plateStress, self.plateStrain,
            self.solidStress, self.solidStrain,
            self.compositePlateStress, self.compositePlateStrain,
            self.ctriaxStress, self.ctriaxStrain,

            # OEE - Strain Energy
            self.strainEnergy,
        ]

        msg = '---ALL RESULTS---\n'
        for result in results:
            for (isubcase, res) in sorted(result.iteritems()):
                msg += 'isubcase = %s\n' % isubcase
                try:
                    msg += str(res) + '\n'
                except:
                    print 'failed on %s' % res.__class__.__name__
                    raise
        return msg


if __name__ == '__main__':
    import sys
    op2_filename = sys.argv[1]

    o = OP2(op2_filename)
    o.read_op2(op2_filename)
    (model, ext) = os.path.splitext(op2_filename)
    f06_outname = model + '.test_op2.f06'
    o.write_f06(f06_outname)
