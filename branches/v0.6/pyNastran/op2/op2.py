#pylint: disable=C0103,W0201,W0223,R0901,R0902,R0904
"""
Main OP2 class
"""
from __future__ import (nested_scopes, generators, division, absolute_import,
                        print_function, unicode_literals)
import os
import sys
import warnings
from numpy import array
from struct import unpack

from pyNastran.op2.fortranFile import FortranFile
from pyNastran.op2.op2Codes import Op2Codes

from pyNastran.op2.tables.resultTable import ResultTable
from pyNastran.op2.tables.geom.geometryTables import GeometryTables

from pyNastran.bdf.bdf import BDF
from pyNastran.f06.f06Writer import F06Writer
from pyNastran.utils.gui_io import load_file_dialog


class OP2Deprecated(object):

    def readOP2(self):
        """
        Reads the op2 file

        .. deprecated: will be replaced in version 0.7 with :func:`read_op2`
        """
        warnings.warn('readOP2 has been deprecated; use '
                      'read_op2', DeprecationWarning, stacklevel=2)
        self.read_op2()

    def setSubcases(self, iSubcases):
        """
        .. deprecated: will be replaced in version 0.7 with :func:`set_subcases`
        """
        warnings.warn('setSubcases has been deprecated; use '
                      'set_subcases', DeprecationWarning, stacklevel=2)
        self.set_subcases(iSubcases)

    def setTransientTimes(self, times):
        """
        .. todo:: rename this
        .. deprecated: will be replaced in version 0.7 with :func:`set_transient_times`
        """
        warnings.warn('setTransientTimes has been deprecated; use '
                      'set_transient_times', DeprecationWarning, stacklevel=2)
        self.set_transient_times(times)


class OP2(BDF,
          FortranFile, Op2Codes, GeometryTables, ResultTable, F06Writer,
          OP2Deprecated):

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

    def is_valid_subcase(self):
        """
        Lets the code check whether or not to read a subcase
        """
        if not self.isAllSubcases:
            if self.isubcase in self.valid_subcases:
                return True
            return False
        return True

    def __init__(self, op2FileName=None, make_geom=False, debug=True, log=None):
        """
        Initializes the OP2 object

        :param op2FileName: the file to be parsed (string or None for GUI)
        :param make_geom: reads the BDF tables (default=False)
        :param debug: prints data about how the OP2 was parsed (default=False)
        :param log: a logging object to write debug messages to
         (.. seealso:: import logging)
        """
        BDF.__init__(self, debug=debug, log=log)
        F06Writer.__init__(self)
        self.set_subcases()  # initializes the variables

        if op2FileName is None:
            wildcard_wx = "Nastran OP2 (*.op2)|*.op2|" \
                "All files (*.*)|*.*"
            wildcard_qt = "Nastran OP2 (*.op2);;All files (*)"
            title = 'Please select a OP2 to load'
            op2FileName = load_file_dialog(title, wildcard_wx, wildcard_qt)

        self.log.debug('op2FileName = %s' % op2FileName)
        bdfExtension = '.bdf'
        f06Extension = '.f06'
        (fname, extension) = os.path.splitext(op2FileName)
        self.table_name = 'temp'

        #: should the BDF tables be parsed
        self.make_geom = make_geom

        #: the input OP2 filename
        self.op2FileName = op2FileName

        #: the expected BDF filename (guessed)
        self.bdfFileName = fname + bdfExtension

        #: the expected F06 filename (guessed)
        self.f06FileName = fname + f06Extension
        #print "bdfFileName = ",self.bdfFileName

        #: developer parameter to write the OP2 is ASCII format
        #: to better understand it
        self.make_op2_debug = False

        #: BDF Title
        self.Title = ''
        #: limit output DTs
        self.expected_times = {}
        #self.expected_times = {1:array([0.1,0.12])}

        #: file object containing the skipped cards
        self.skippedCardsFile = open('skippedCards.out', 'a')

        #: the list of supported tables (dont edit this)
        self.tablesToRead = [
            # nodes/geometry/loads/BCs
            'GEOM1', 'GEOM2', 'GEOM3', 'GEOM4',      # regular
            'GEOM1S', 'GEOM2S', 'GEOM3S', 'GEOM4S',  # superelements

            'GEOM1OLD', 'GEOM1N',  # ???
            'EPT', 'MPT',    # properties/materials
            'EPTS', 'MPTS',  # properties/materials - superelements
            'EDTS',          # ???
            'DYNAMIC', 'DYNAMICS',
            'DIT',                           # tables (e.g. TABLED1)
            'LAMA', 'BLAMA',                 # eigenvalues

            'BGPDT', 'BGPDTS',               # boundary grids???
            'EQEXIN', 'EQEXINS', 'PVT0', 'CASECC', 'EDOM', 'CASEXX',
            'DESTAB',                        # design variables
            'OQG1', 'OQGV1',                 # spc forces
            'OQMG1',                         # mpc forces

            'OUGV1', 'OUG1',                 # displacements
            'OGPFB1',                        # grid point forces
            'OGS1',                          # grid point stresses

            'OEF1X', 'DOEF1', 'OEFIT',        # element forces
            'OPG1', 'OPGV1', 'OPNL1',          # applied forces
            'OES1', 'OES1X', 'OES1X1', 'OES1C',  # stress
            'OSTR1C', 'OSTR1X',               # strains
            'OESNLXR', 'OESNLXD', 'OESNL1X', 'OESNLBR',  # nonlinear stress

            'OESCP',                       # cylinder stress???
            'OESTRCP',                     # cylinder strain???
            'OESRT',                       # rotational stress?

            'ONRGY1',  # energy
            'ONRGY2',  # energy (sort2, unsupported)

            'R1TABRG', 'HISADD',  # SOL 200

            # unsupported frequency results
            'OAGPSD2', 'OAGATO2', 'OAGRMS2', 'OAGNO2', 'OAGCRM2',
            'OEFPSD2', 'OEFATO2', 'OEFRMS2', 'OEFNO2', 'OEFCRM2',
            'OESPSD2', 'OESATO2', 'OESRMS2', 'OESNO2', 'OESCRM2',
            'OPGPSD2', 'OPGATO2', 'OPGRMS2', 'OPGNO2', 'OPGCRM2',
            'OQGPSD2', 'OQGATO2', 'OQGRMS2', 'OQGNO2', 'OQGCRM2',

            # supported-ish
            'OQMPSD2', 'OQMATO2', 'OQMRMS2', 'OQMNO2', 'OQMCRM2',
            'OSTRPSD2', 'OSTRATO2', 'OSTRRMS2', 'OSTRNO2', 'OSTRCRM2',
            'OUGPSD2', 'OUGATO2', 'OUGCRM2', 'OUGNO2', 'OUGRMS2',
            'OVGPSD2', 'OVGATO2', 'OVGRMS2', 'OVGNO2', 'OVGCRM2',

            # TODO what do these do???
            'OUPV1',
            'VIEWTB', 'ERRORN',
            'OFMPF2M', 'OSMPF2M', 'OPMPF2M', 'OGPMPF2M', 'OLMPF2M',
            'PCOMPTS',
            'OMM2',
            'EDOM',
            'STDISP', 'SDF', 'MONITOR', 'AEMONPT',
            'OGPWG',  # grid point weight
            'OQP1',
            'OCRUG', 'OCRPG',

            # new
            'AFRF', 'AGRF',
            'PMRF', 'PERF', 'PFRF',
            'FOL',

            'DBCOPT','CONTACT',
        ]

        #: list of OP2 tables that were read
        #: mainly for debugging
        self.tablenames = []

        self.__objects_init__()

    def __objects_init__(self):
        """More variable declarations"""
        #: ESE
        self.eigenvalues = {}

        #: OUG - displacement
        self.displacements = {}           # tCode=1 thermal=0
        self.displacementsPSD = {}        # random
        self.displacementsATO = {}        # random
        self.displacementsRMS = {}        # random
        self.displacementsCRM = {}        # random
        self.displacementsNO = {}        # random
        self.scaledDisplacements = {}     # tCode=1 thermal=8

        #: OUG - temperatures
        self.temperatures = {}           # tCode=1 thermal=1

        #: OUG - eigenvectors
        self.eigenvectors = {}            # tCode=7 thermal=0

        #: OUG - velocity
        self.velocities = {}              # tCode=10 thermal=0

        #: OUG - acceleration
        self.accelerations = {}           # tCode=11 thermal=0

        # OEF - Forces - tCode=4 thermal=0
        # rods
        self.rodForces = {}
        self.conrodForces = {}
        self.ctubeForces = {}

        self.barForces = {}
        self.bar100Forces = {}
        self.beamForces = {}
        self.bendForces = {}
        self.bushForces = {}
        self.coneAxForces = {}
        self.damperForces = {}
        self.gapForces = {}
        self.plateForces = {}
        self.plateForces2 = {}
        self.shearForces = {}
        self.solidPressureForces = {}
        self.springForces = {}
        self.viscForces = {}

        self.force_VU = {}
        self.force_VU_2D = {}

        #OEF - Fluxes - tCode=4 thermal=1
        self.thermalLoad_CONV = {}
        self.thermalLoad_CHBDY = {}
        self.thermalLoad_1D = {}
        self.thermalLoad_2D_3D = {}
        self.thermalLoad_VU = {}
        self.thermalLoad_VU_3D = {}
        self.thermalLoad_VUBeam = {}
        #self.temperatureForces = {}

        # OES - tCode=5 thermal=0 s_code=0,1 (stress/strain)
        #: OES - CELAS1/CELAS2/CELAS3/CELAS4 stress
        self.celasStress = {}
        #: OES - CELAS1/CELAS2/CELAS3/CELAS4 strain
        self.celasStrain = {}

        #: OES - CTRIAX6
        self.ctriaxStress = {}
        self.ctriaxStrain = {}

        #: OES - isotropic CROD/CONROD/CTUBE stress
        self.rodStress = {}
        self.conrodStress = {}
        self.ctubeStress = {}

        #: OES - isotropic CROD/CONROD/CTUBE strain
        self.rodStrain = {}
        self.conrodStrain = {}
        self.ctubeStrain = {}

        #: OES - nonlinear CROD/CONROD/CTUBE stress
        self.nonlinearRodStress = {}
        #: OES - nonlinear CROD/CONROD/CTUBE strain
        self.nonlinearRodStrain = {}
        #: OES - isotropic CBAR stress
        self.barStress = {}
        #: OES - isotropic CBAR strain
        self.barStrain = {}
        #: OES - isotropic CBEAM stress
        self.beamStress = {}
        #: OES - isotropic CBEAM strain
        self.beamStrain = {}
        #: OES - isotropic CBUSH stress
        self.bushStress = {}
        #: OES - isotropic CBUSH strain
        self.bushStrain = {}
         #: OES - isotropic CBUSH1D strain/strain
        self.bush1dStressStrain = {}

        #: OES - isotropic CTRIA3/CQUAD4 stress
        self.plateStress = {}
        #: OES - isotropic CTRIA3/CQUAD4 strain
        self.plateStrain = {}
        #: OESNLXR - CTRIA3/CQUAD4 stress
        self.nonlinearPlateStress = {}
        #: OESNLXR - CTRIA3/CQUAD4 strain
        self.nonlinearPlateStrain = {}
        self.hyperelasticPlateStress = {}
        self.hyperelasticPlateStrain = {}

        #: OES - isotropic CTETRA/CHEXA/CPENTA stress
        self.solidStress = {}
        #: OES - isotropic CTETRA/CHEXA/CPENTA strain
        self.solidStrain = {}
        #: OES - composite CTRIA3/CQUAD4 stress
        self.compositePlateStress = {}
        #: OES - composite CTRIA3/CQUAD4 strain
        self.compositePlateStrain = {}

        #: OES - CSHEAR stress
        self.shearStress = {}
        #: OES - CSHEAR strain
        self.shearStrain = {}

        #: OES - CELAS1 224, CELAS3 225,
        self.nonlinearSpringStress = {}
        #: OES - GAPNL 86
        self.nonlinearGapStress = {}
        #: OES - CBUSH 226
        self.nolinearBushStress = {}

        # OQG - spc/mpc forces
        self.spcForces = {}  # tCode=3?
        self.mpcForces = {}  # tCode=39

        # OQG - thermal forces
        self.thermalGradientAndFlux = {}

        #: OGF - grid point forces
        self.gridPointForces = {}  # tCode=19

        #: OGS1 - grid point stresses
        self.gridPointStresses = {}       # tCode=26
        self.gridPointVolumeStresses = {}  # tCode=27

        #: OPG - summation of loads for each element
        self.loadVectors = {}       # tCode=2  thermal=0
        self.thermalLoadVectors = {}  # tCode=2  thermal=1
        self.appliedLoads = {}       # tCode=19 thermal=0
        self.forceVectors = {}       # tCode=12 thermal=0

        #: OEE - strain energy density
        self.strainEnergy = {}  # tCode=18

    def get_op2_stats(self):
        """
        Gets info about the contents of the different attributes of the
        OP2 class
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
                if hasattr(subcase,'get_stats'):
                    msg.append('op2.%s[%s]\n' % (table_type, isubcase))
                    msg.extend(subcase.get_stats())
                    msg.append('\n')
                else:
                    msg.append('skipping %s op2.%s[%s]\n\n' % (subcase.__class__.__name__, table_type, isubcase))
                    #raise RuntimeError('skipping %s op2.%s[%s]\n\n' % (subcase.__class__.__name__, table_type, isubcase))

        return ''.join(msg)

    def read_tape_code(self):
        """
        Reads the OP2 header.  This table is still very much in development.

        .. todo:: whats in this table?
        """
        #self.print_section(500)
        #sys.exit('op2-readTapeCode')

        self.read_markers([3])
        #print(self.print_section(20))
        ints = self.read_int_block()
        if self.make_op2_debug:
            self.op2Debug.write('%s\n' % str(ints))
        #print "*ints = ",ints
        self.read_markers([7])

        word = self.read_string_block()  # Nastran Fort Tape ID Code -
        #print "word = |%r|" %(word)

        self.read_markers([2])
        ints = self.read_int_block()
        #print "*ints = ",ints

        self.read_markers([-1])

        #data = self.get_data(60)
        #self.print_block(data)

    def read_tape_code_post2(self):
        # PVTO
        #self.read_markers([2])
        #ints = self.read_int_block()
        #print("*ints = ",ints)

        #ints = self.read_block()
        #print("*ints = ",ints)
        #print(self.print_section(40))

        block = self.read_new_block()
        #print("len(block) = %s" %(len(block)))
        if len(block) == 12:  # post = -1
            form = -1
            (month, day, year) = unpack(b'iii', block)
            print("date = %s/%s/%s" %(month, day, 2000+year))

            block = self.read_new_block()
            #print("len(block) = %s" % len(block))
            print(self.print_block(block))

            block = self.read_new_block()
            #print("len(block) = %s" % len(block))
            print(self.print_block(block))

            self.read_markers([-1, 0])

            block = self.read_new_block()
            #print("len(block) = %s" % len(block))
            table_name, = unpack(b'8s', block)
            self.log.debug("table_name=%r self.n=%s" % (table_name, self.n))

            marker = self.read_marker()
            print('marker = %s' % marker)
            #self.read_markers([-1])
            block = self.read_new_block()
            print(self.print_block(block))

            self.read_markers([-2])
            block = self.read_new_block()
            print(self.print_block(block))

            marker = self.read_marker()
            print('marker = %s' % marker)

            #self.read_markers([-1])
            #self.readPartOfNewTable()  # PVTO

            #self.readNewTable()
            #self.readNewTable()
            print(self.print_section(100))

            sys.exit('stoppingD')

        elif len(block) == 8: # post = -2
            self.n = 0
            self.op2.seek(self.n)
            table_name, = unpack(b'8s', block)
            self.log.debug("table_name = %r\n" % table_name)
            #sys.exit('stoppingA')

            block = self.read_new_block()
            #print(self.print_block(block))
            self.read_markers([-1])
            #print("")

            self.read_part_of_new_table()  # PVTO

            block = self.read_new_block()
            print(self.print_block(block))
            self.read_part_of_new_table()  # PRTMAXIM=YES AUTOSPC=NO

            #block = self.readNewBlock()
            #print(self.print_section(200))
            #return

            self.read_new_table()
            self.read_new_table()
            self.read_new_table()
            #print(self.print_section(40))
            sys.exit('stoppingB')

            return


        sys.exit()

        # -----------------------

        foundMoreTables = True
        while foundMoreTables and table_name:
            self.log.debug("table_name=%r self.n=%s" % (table_name, self.n))
            print("-----------------------------")
            try:
                n = self.n
                table_name = self.read_new_table()
            except:
                self.op2.seek(n)
                foundMoreTables = False

        print("**************")
        self.op2.seek(n)
        print(self.print_section(240))
        sys.exit('stoppingC')

    def read_new_table(self):
        self.op2.read(8)
        self.n += 8
        data = self.op2.read(8)
        self.n += 12

        if len(data) < 8:
            return

        self.op2.seek(self.n)
        table_name, = unpack(b'8s', data)
        self.log.debug("table_name = %r\n" % table_name)

        self.read_markers([-1])
        block = self.read_new_block()
        #print(self.print_block(block))

        self.read_part_of_new_table()
        return table_name

    def read_part_of_new_table(self):
        n = -2
        keepGoingOnTable = True
        while keepGoingOnTable:
            print("n = %s" % n)
            try:
                nStar = self.n
                self.read_markers([n, 1, 0])
                block = self.read_new_block()
                #print(self.print_block(block))
                n -= 1
            except:
                self.n = nStar
                self.op2.seek(nStar)
                keepGoingOnTable = False

    def read_new_block(self):
        data = self.op2.read(16)
        #print(self.print_block(data))
        (four, n, four, fourN) = unpack(b'4i', data)
        #print('n = %s' %(n))
        self.n += 16

        #if n > 70000:
            #asf
        data = self.op2.read(n * 4)
        self.n += n * 4 + 4
        #print("self.n = ",self.n)
        self.op2.seek(self.n)

        return data

    def read_op2(self):
        #: the OP2 file object
        self.op2 = open(self.op2FileName, 'rb')
        try:
            if self.make_op2_debug:
                #: a developer debug file (largely unsupported)
                self.op2Debug = open('debug.out', 'wb')
            #: the byte position in the OP2
            self.n = self.op2.tell()

            try:
                #self.readTapeCodePost2()
                self.read_tape_code()
            except:
                #raise
                msg = ('When this happens, the analysis failed or '
                      'the code bombed...check the F06.\n'
                      'If the F06 is OK:\n'
                      '  1.  Make sure you used PARAM,POST,-1 in your '
                      'BDF/DAT/NAS\n'
                      '  2.  Make sure the following line is not included in the BDF\n'
                      "      ASSIGN OUTPUT2 = '%s', UNIT = 12, FORM = FORMATTED\n"
                      '  3.  Run the problem on a different Operating System\n'
                      '  4.  Are you running an OP2? :)  \n'
                      'fname=%s' % (self.op2FileName, self.op2FileName))
                #raise RuntimeError("Tape Code Error: %s" % msg)
                raise

            isAnotherTable = True
            while isAnotherTable:
                self.log.debug('-' * 80)
                try:
                    table_name = self.read_table_name(rewind=True, stopOnFailure=False)
                    self.log.debug("table_name = %r" % table_name)
                except EOFError:  # the isAnotherTable method sucks...
                    isAnotherTable = False
                    self.log.debug("***ok exit, but it could be better...")
                    break
                except SyntaxError:  # Invalid Markers, the isAnotherTable method sucks...
                    isAnotherTable = False
                    self.log.debug("***poor exit, but it worked...")
                    #raise
                    break
                except:
                    raise
                self.log.debug("table_name = %r" % table_name)
                #print("table_name = |%r|" % table_name)

                if table_name is None:
                    break
                else:
                    isAnotherTable = self.read_table(table_name)

            self.log.debug("---end of all tables---")
        except:
            self.op2.close()
            self.skippedCardsFile.close()
            raise
        self.op2.close()
        self.skippedCardsFile.close()

    def read_table(self, table_name):
        if table_name in self.tablesToRead:
            self.table_name = table_name
            self.isRegular = True
            try:
                #print("startTell = %s" %(self.op2.tell()))
                if table_name == 'GEOM1':  # nodes,coords,etc.
                    self.readTable_Geom1()
                elif table_name == 'GEOM1S':  # superelements
                    self.readTable_Geom1S()  #  - nodes,coords,etc.
                elif table_name == 'GEOM2S':  # superelements
                    self.readTable_Geom2S()  #  - elements
                elif table_name == 'GEOM3S':  # superelements
                    self.readTable_Geom3S()  #  - static/thermal loads
                elif table_name == 'GEOM4S':  # superelements
                    self.readTable_Geom4S()  #  - constraints

                #elif table_name=='GEOM1OLD':
                #    self.readTable_Geom1Old()
                #elif table_name=='GEOM1N':
                #    self.readTable_Geom1N()
                elif table_name == 'GEOM2':  # elements
                    self.readTable_Geom2()
                elif table_name == 'GEOM3':  # static/thermal loads
                    self.readTable_Geom3()
                elif table_name == 'GEOM4':  # constraints
                    self.readTable_Geom4()

                elif table_name in ['EPT', 'EPTS']:  # element properties
                    self.readTable_EPT()
                elif table_name in ['MPT', 'MPTS']:  # material properties
                    self.readTable_MPTS()
                elif table_name in ['DYNAMIC', 'DYNAMICS']:  # dyanmic info
                    self.readTable_DYNAMICS()
                elif  table_name in ['DIT']:  # tables...
                    self.readTable_DIT()     # TABLED1/TABLEM1/TABLES1/GUST
                elif table_name in ['LAMA', 'BLAMA']:  # eigenvalue
                    self.readTable_LAMA()

                elif table_name in ['VIEWTB', 'EQEXIN', 'EQEXINS', 'OEFIT',
                                   'GEOM1N', 'OGPWG', 'GEOM1OLD']:
                    self.readTable_DUMMY_GEOM(table_name)
                elif table_name in ['OMM2']:
                    self.readTable_OMM2()
                elif table_name in ['DESTAB']:  # design variable table
                    self.readTable_DesTab()
                elif table_name in ['R1TABRG']:  # not done - response table
                    self.readTable_R1TAB()
                    self.isOptimization = True
                elif table_name in ['HISADD']:  # not done
                    self.readTable_R1TAB()
                    self.isOptimization = True
                elif table_name in ['ERRORN']:  # not done
                    self.readTable_R1TAB()

                elif table_name in ['OPG1', 'OPNL1', 'OGS1', 'OPGV1']:
                    self.readTable_OPG()  # table of applied loads
                elif table_name in ['OGPFB1', ]:
                    self.readTable_OGF()
                elif table_name in ['OCRUG', 'OCRPG']:  # totally guessing...
                    self.readTable_OUG()

                elif table_name in ['OEF1X', 'DOEF1', 'OEFPSD2', 'OEFATO2',
                                   'OEFRMS2', 'OEFNO2', 'OEFCRM2', ]:
                    self.readTable_OEF()  # applied loads
                elif table_name in ['OQG1', 'OQGV1', 'OQP1', ]:  # spc forces
                    self.readTable_OQG()
                elif table_name in ['OQMG1', 'OQMPSD2', 'OQMATO2', 'OQMRMS2',
                                   'OQMNO2', 'OQMCRM2', ]:  # mpc forces
                    #self.readTable_OQG()
                    self.readTable_DUMMY_GEOM(table_name)

                elif table_name in ['OUGV1', 'OUPV1', 'OUG1']:
                    self.readTable_OUG()  # displacements/velocity/acceleration
                elif table_name in ['OUGPSD2', 'OUGATO2', 'OUGRMS2', 'OUGNO2',
                                   'OUGCRM2']:  # OUG tables???
                    self.readTable_OUG2()
                    #self.readTable_OUG()

                elif table_name in ['OES1', 'OES1X', 'OES1X1', 'OSTR1X',
                                   'OES1C', 'OESCP', 'OESRT', 'OESNLXR',
                                   'OESNL1X']:
                    self.readTable_OES()  # stress
                elif table_name in ['OSTR1X', 'OSTR1C', ]:
                    self.readTable_OES()  # strain
                elif table_name in ['OESTRCP', 'OESNLXD', 'OESNLXR', ]:
                    self.readTable_OES()  # ??? stress/strain
                elif table_name in ['OSTRATO2', 'OSTRPSD2', 'OESRMS2', 'OESNO2',
                                   'OESCRM2', 'OSTRRMS2', 'OESRMS2', 'OSTRNO2',
                                   'OESCRM2', 'OSTRCRM2', ]:  # unhandled
                    self.readTable_OES()

                #elif table_name in ['OESNLXD',]: # dont use this, testing only
                    #self.readTable_OES() # NLXD
                #elif table_name in ['OESNLXR',]: # dont use this
                    #self.readTable_OES()  # NLXR

                elif table_name in ['ONRGY1']:  # energy
                    self.readTable_OEE()
                elif table_name in ['ONRGY2']:
                    self.readTable_OEE()

                elif table_name in ['PCOMPTS']:
                    self.readTable_PCOMPTS()
                elif table_name in ['SDF']:  # ???
                    self.readTable_SDF()
                #elif table_name in ['CASECC']:
                    #self.readTable_CASECC()

                # not done
                elif table_name in []:
                    self.readTableB_DUMMY()
                elif table_name in ['MONITOR']:
                    self.read_monitor()
                elif table_name in ['PMRF', 'PERF', 'PFRF',
                                   'AEMONPT', 'FOL', 'AFRF', 'AGRF', ]:
                    self.readTableB_DUMMY()
                #elif table_name in []:
                #    self.readTableB_DUMMY()

                elif table_name in ['STDISP', 'FOL', 'OFMPF2M', 'OSMPF2M',
                                   'OPMPF2M', 'OGPMPF2M', 'OLMPF2M',
                                   'OVGPSD2']:
                    self.readTable_DUMMY_GEOM(table_name)
                elif table_name in ['OVGATO2', 'OVGRMS2', 'OVGNO2']:
                    self.readTable_DUMMY_GEOM(table_name)

                elif table_name in ['OESNLXR', 'OESNL1X', 'OESPSD2', 'OESNLBR',
                                   'OESATO2', ]:
                    self.readTable_DUMMY_GEOM(table_name)
                elif table_name in ['OVGCRM2', 'OAGPSD2', 'OAGATO2', 'OAGRMS2',
                                   'OAGNO2', 'OAGCRM2', 'OPGPSD2', 'OPGPSD2',
                                   'OPGPSD2', 'OPGATO2']:
                    self.readTable_DUMMY_GEOM(table_name)
                elif table_name in ['OPGRMS2', 'OPGNO2', 'OPGCRM2', 'OQGPSD2']:
                    self.readTable_DUMMY_GEOM(table_name)
                elif table_name in ['OQGPSD2', 'OQGATO2', 'OQGRMS2', 'OQGNO2',
                                   'OQGCRM2', 'PVT0', 'CASECC', 'CASEXX', 'EDOM', 'DBCOPT',]:
                    self.readTable_DUMMY_GEOM(table_name)
                elif table_name in ['BGPDT', 'BGPDTS', 'EDTS', 'CONTACT']:
                    self.readTable_DUMMY_GEOM(table_name)
                else:
                    msg = 'unhandled table_name=%r' % table_name
                    raise KeyError(msg)
                #print("endTell   = ",self.op2.tell())
                #print("---isAnotherTable---")
                (is_another_table) = self.has_more_tables()
                #is_another_table = True
            except EOFError:
                is_another_table = False
        else:
            if table_name not in [None]:
                raise NotImplementedError('%r is not supported' % table_name)
            (is_another_table) = self.skip_next_table()
            #return isAnotherTable
        #print(self.print_section(140))
        self.log.debug("*** finished table_name = %r" % table_name)
        return is_another_table

    def read_monitor(self):
        table_name = self.read_table_name(rewind=False)  # LAMA
        self.table_init(table_name)
        #print("tablename1 = |%r|" % table_name)
        #print("tablename2 = |%r|" % self.table_name)

        self.read_markers([-1, 7], 'MONITOR')
        ints = self.read_int_block()
        #print("*ints = ", ints)

        self.read_markers([-2, 1, 0], 'MONITOR')
        buffer_words = self.get_marker()
        #print("buffer_words = ", buffer_words)

        word = self.read_string_block()  # MONITOR
        #print("word = |%s|" % word)

        self.read_markers([-3, 1, 0], 'MONITOR')
        buffer_words = self.get_marker()
        #print("buffer_words = ", buffer_words, buffer_words * 4)

        data = self.op2.read(4)
        data = self.op2.read(buffer_words * 4)
        Format = str(buffer_words*4) + 's'
        Format = bytes(Format)
        word, = unpack(Format, data)
        #print("word = ", word)
        data = self.op2.read(4)
        self.n += buffer_words * 4 + 8


        self.read_markers([-4, 1, 0], 'MONITOR')
        buffer_words = self.get_marker()
        #print("buffer_words = ", buffer_words, buffer_words * 4)
        data = self.op2.read(4)
        data = self.op2.read(buffer_words * 4)
        Format = str(buffer_words * 4) + 's'
        Format = bytes(Format)
        word, = unpack(Format, data)
        #print("word = ", word)
        data = self.op2.read(4)
        self.n += buffer_words * 4 + 8

        self.read_markers([-5, 1, 0], 'MONITOR')

        #word = self.read_string_block()  # MONITOR
        #print("word = |%s|" %(word))

        #print(self.print_section(200))
        #sys.exit()

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
        #print("***sort_code = ",self.sort_code)
        while sort_code > 0:
            value = sort_code % 2
            sort_code = (sort_code - value) // 2
            bits[i] = value
            #print("    *bit = %s" %(value))
            #print("    sort_code = %s" %(sort_code))
            i -= 1
        #print("sort_bits = %s" %(bits))
        #: the bytes describe the SORT information
        self.sort_bits = bits

        self.data_code['sort_bits'] = self.sort_bits

    def parse_approach_code(self, data):
        """
        int3 is the 3rd word in table=-3 and may be
        element_type or something else depending on the table type
        """
        (aCode, tCode, int3, isubcase) = unpack(b'iiii', data[:16])
        #: the local subcase ID
        self.isubcase = isubcase
        #print("isubcase = %s" %(isubcase))
        self.subcases.add(self.isubcase)  # set notation

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

        #: data_code stores the active variables; these pass important
        #: self variables into the result object
        self.data_code = {'analysis_code': self.analysis_code,
                         'device_code': self.device_code,
                         'table_code': self.table_code,
                         'sort_code': self.sort_code,
                         'dt': None,
                         'log': self.log,
                         }
        #print("isubcase = ",self.isubcase)
        self._parse_sort_code()

        #print('aCode(1)=%s analysis_code=%s device_code=%s '
        #      'tCode(2)=%s table_code=%s sort_code=%s isubcase(4)=%s'
        #      %(aCode, self.analysis_code, self.device_code, tCode,
        #        self.table_code, self.sort_code, self.isubcase))
        #self.log.debug(self.print_table_code(self.table_code))
        return (int3)

    def parse_approach_code2(self, data):
        """
        int3 is the 3rd word in table=-3 and may be
        element_type or something else depending on the table type
        """
        (aCode, tCode, int3, ID) = unpack(b'iiii', data[:16])
        #: the local subcase ID
        self.ID = ID
        #print("isubcase = %s" %(isubcase))
        self.subcases.add(self.ID)  # set notation

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

        #: data_code stores the active variables; these pass important
        #: self variables into the result object
        self.data_code = {'analysis_code': self.analysis_code,
                         'device_code': self.device_code,
                         'table_code': self.table_code,
                         'sort_code': self.sort_code,
                         'dt': None,
                         'log': self.log,
                         }
        #print("isubcase = ",self.isubcase)
        self._parse_sort_code()

        #print('aCode(1)=%s analysis_code=%s device_code=%s '
        #      'tCode(2)=%s table_code=%s sort_code=%s isubcase(4)=%s'
        #      %(aCode, self.analysis_code, self.device_code, tCode,
        #        self.table_code, self.sort_code, self.isubcase))
        #self.log.debug(self.print_table_code(self.table_code))
        return (int3)

    def get_values(self, data, iFormat, iWordStart, iWordStop=None):
        """
        Extracts the ith word from the data structure as the provided type
        supports multiple inputs with iWordStop (note this is words, not
        outputs)

        :param self: the object pointer
        :param data: the binary data that is as long as the buffer size
        :param iWordStart: the word to start reading from
        :param iWordStop:  the word to stop reading on (largely unused)
        .. warning:: works with nastran syntax, not standard python syntax
            this makes it work with what's documented in the DMAP manual
        """
        if iWordStop is None:
            #print("iWordStart=%s data[%s:%s]" %(iWordStart,iWordStart*4,
            #                                   (iWordStart+1)*4))
            ds = data[(iWordStart - 1) * 4:iWordStart * 4]
            iFormat = bytes(iFormat)
            return unpack(iFormat, ds)[0]
        #print("type(data) = ",type(data))
        ds = data[(iWordStart - 1) * 4:(iWordStop - 1) * 4]
        return unpack(iFormat, ds)

    def _delete_attributes(self, params):
        """
        Deletes any parameters before going to the next table to avoid messing
        up data
        """
        params += ['data_code', 'device_code', 'analysis_code', 'table_code',
                   'sort_code', 'isubcase', 'data', 'num_wide',
                   'nonlinear_factor', 'obj', 'subtitle', 'label']
        for param in params:
            if hasattr(self, param):
                delattr(self, param)

    def get_buffer_words(self):
        buffer_words = self.get_marker()
        #print("buff_marker = |%s|" % buffer_words)
        #print("buffer_words = ",buffer_words, buffer_words * 4)
        if buffer_words <= 0:
            msg = ('An invalid buffer_size was found...buffer_words=%s '
                   'table_name=%s section=\n%s' % (buffer_words,
                    self.table_name, self.print_section(200)))
            raise BufferError(msg)
        return buffer_words

    def verify_buffer_size(self, buffer_words):
        assert buffer_words > 0, self.print_section(220)

    def read_title(self):
        r"""
        Reads the Title, Subtitle, and Label.

        Puts them in self.iSubcaseNameMap[isubcase] = [Subtitle,Label]
        """
        #: the title of the analysis
        word = self.read_string(384)  # titleSubtitleLabel
        self.Title = word[0:128].strip()
        #: the subtitle of the subcase
        self.subtitle = word[128:256].strip()
        #: the label of the subcase
        self.label = word[256:328].strip()

        # not really a hollerith, just the end of the block (so buffer_words*4)
        self.read_hollerith()

        self.data_code['subtitle'] = self.subtitle
        self.data_code['label'] = self.label

        if hasattr(self,'isubcase'):
            if self.isubcase not in self.iSubcaseNameMap:
                self.iSubcaseNameMap[self.isubcase] = [self.subtitle, self.label]
        else:
            pass

    def table_init(self, word):
        """
        Starts a new table
        """
        #: the local table name
        self.table_name = word.strip()
        #: the names of all the tables
        self.tablenames.append(word)
        #msg = '*' * 20 + word + '*' * 20 + '\n'
        #print(msg)

    def get_table_names_from_op2(self):
        """
        Returns the list of parsed tables
        """
        return self.tablenames

    def print_results(self):
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
                msg += 'isubcase = %s\n' % (isubcase)
                try:
                    msg += str(res) + '\n'
                except:
                    print('failed on %s' % (res.__class__.__name__))
                    raise
        return msg