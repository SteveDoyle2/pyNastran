# pylint: disable=C0103,C0302,R0902,R0904,R0914,W0201,W0611
from __future__ import (nested_scopes, generators, division, absolute_import,
                        print_function, unicode_literals)
import os
import sys
import warnings

from pyNastran.general.general import ListPrint

from .cards.elements.elements import CFAST, CGAP, CRAC2D, CRAC3D
from .cards.properties.properties import (PFAST, PGAP, PLSOLID, PSOLID,
                                          PRAC2D, PRAC3D, PCONEAX)

from .cards.elements.springs import (CELAS1, CELAS2, CELAS3, CELAS4,
                                     SpringElement)
from .cards.properties.springs import PELAS, PELAST

from .cards.elements.solid import (CTETRA4, CTETRA10, CPENTA6, CPENTA15,
                                   CHEXA8, CHEXA20, SolidElement)
from .cards.elements.rigid import (RBAR, RBAR1, RBE1, RBE2, RBE3, RigidElement)

from .cards.elements.shell import (CQUAD, CQUAD4, CQUAD8, CQUADR, CQUADX,
                                   CSHEAR, CTRIA3, CTRIA6, CTRIAX,
                                   CTRIAX6, CTRIAR, ShellElement)
from .cards.properties.shell import PSHELL, PCOMP, PCOMPG, PSHEAR
from .cards.elements.bush import CBUSH, CBUSH1D, CBUSH2D
from .cards.properties.bush import PBUSH, PBUSH1D
from .cards.elements.damper import (CVISC, CDAMP1, CDAMP2, CDAMP3, CDAMP4,
                                    CDAMP5, DamperElement)
from .cards.properties.damper import (PVISC, PDAMP, PDAMP5, PDAMPT)
from .cards.elements.bars import (CROD, CONROD, CTUBE, CBAR, CBEAM, CBEAM3,
                                  CBEND, LineElement, RodElement)
from .cards.properties.bars import (PROD, PTUBE, PBAR, PBARL,
                                    PBEAM, PBEAML, PBCOMP)  # PBEND
from .cards.elements.mass import (CONM1, CONM2, CMASS1, CMASS2, CMASS3, CMASS4,
                                  PointElement, PointMassElement)  # CMASS5
from .cards.properties.mass import (PMASS, NSM)
from .cards.aero import (AEFACT, AELINK, AELIST, AEPARM, AESTAT, AESURF,
                         AESURFS, AERO, AEROS, CSSCHD, CAERO1, CAERO2, CAERO3,
                         CAERO4, CAERO5, FLFACT, FLUTTER, GUST, MKAERO1,
                         MKAERO2, PAERO1, PAERO2, SPLINE1, SPLINE2, SPLINE4,
                         SPLINE5, TRIM)
from .cards.constraints import (SPC, SPCADD, SPCD, SPCAX, SPC1,
                                MPC, MPCADD, SUPORT1, SUPORT,
                                constraintObject2)
from .cards.coordinateSystems import (CORD1R, CORD1C, CORD1S,
                                      CORD2R, CORD2C, CORD2S, CORD3G)
from .cards.dmig import (DEQATN, DMIG, DMI, DMIJ, DMIK, DMIJI, NastranMatrix)
from .cards.dynamic import (FREQ, FREQ1, FREQ2, FREQ4, TSTEP, TSTEPNL, NLPARM)
from .cards.loads.loads import (LSEQ, SLOAD, DLOAD, DAREA, TLOAD1, TLOAD2,
                                RLOAD1, RLOAD2, RANDPS, RFORCE)
from .cards.loads.staticLoads import (LOAD, GRAV, ACCEL1, FORCE,
                                      FORCE1, FORCE2, MOMENT, MOMENT1, MOMENT2,
                                      PLOAD, PLOAD1, PLOAD2, PLOAD4, PLOADX1)

from .cards.materials import (MAT1, MAT2, MAT3, MAT4, MAT5,
                              MAT8, MAT9, MAT10,
                              MATHP, CREEP, MATS1)
from .cards.methods import (EIGB, EIGC, EIGR, EIGP, EIGRL)
from .cards.nodes import GRID, GRDSET, SPOINTs
from .cards.optimization import (DCONSTR, DESVAR, DDVAL, DOPTPRM, DLINK,
                                 DRESP1, DRESP2, DVMREL1, DVPREL1, DVPREL2)
from .cards.params import PARAM
from .cards.sets import (ASET, BSET, CSET, QSET,
                         ASET1, BSET1, CSET1, QSET1,
                         SET1, SET3, SESET, SEQSEP, RADSET)
from .cards.thermal.loads import (QBDY1, QBDY2, QBDY3, QHBDY, TEMP, TEMPD)
from .cards.thermal.thermal import (CHBDYE, CHBDYG, CHBDYP, PCONV, PCONVM,
                                    PHBDY, CONV, RADM, RADBC,)
from .cards.tables import (TABLED1, TABLED2, TABLED3,
                           TABLEM1, TABLEM2, TABLEM3, TABLEM4,
                           TABLES1, TABLEST, TABRND1, TABRNDG, TIC)
from pyNastran.bdf.caseControlDeck import CaseControlDeck
from pyNastran.bdf.bdf_Methods import BDFMethods
from .bdfInterface.getCard import GetMethods
from .bdfInterface.addCard import AddMethods
from .bdfInterface.BDF_Card import BDFCard
from .bdfInterface.bdf_Reader import BDFReader
from .bdfInterface.bdf_writeMesh import WriteMesh
from .bdfInterface.bdf_cardMethods import CardMethods
from .bdfInterface.crossReference import XrefMesh

class BDF(BDFReader, BDFMethods, GetMethods, AddMethods, WriteMesh,
          CardMethods, XrefMesh):
    """
    NASTRAN BDF Reader/Writer/Editor class.
    """
    modelType = 'nastran'
    isStructured = False

    def __init__(self, debug=True, log=None, nCardLinesMax=1000):
        """
        Initializes the BDF object
        @param self
          the object pointer
        @param debug
          used to set the logger if no logger is passed in
        @param log
          a python logging module object
        @param nCardLinesMax
          the number of lines of the longest card in the deck (default=1000)
        """
        ## allows the BDF variables to be scoped properly (i think...)
        BDFReader.__init__(self, debug, log)
        GetMethods.__init__(self)
        AddMethods.__init__(self)
        BDFMethods.__init__(self)
        WriteMesh.__init__(self)
        CardMethods.__init__(self, nCardLinesMax)
        XrefMesh.__init__(self)

        ## useful in debugging errors in input
        self.debug = debug
        self._init_solution()

        ## flag that allows for OpenMDAO-style optimization syntax to be used
        self._is_dynamic_syntax = False
        ## lines that were rejected b/c they were for a card
        ## that isnt supported
        self.rejects = []
        ## cards that were created, but not processed
        self.reject_cards = []
        ## list of execive control deck lines
        self.executive_control_lines = []
        ## list of case control deck lines
        self.case_control_lines = []

        ## the list of possible cards that will be parsed
        self.cardsToRead = set([
            'PARAM',
            'GRID', 'GRDSET', 'SPOINT',  # 'RINGAX',
            'POINT',

            # elements
            'CONM2', 'CMASS1', 'CMASS2', 'CMASS3', 'CMASS4',
            # 'CONM1',
            'CELAS1', 'CELAS2', 'CELAS3', 'CELAS4',
            # 'CELAS5',
            'CBUSH', 'CBUSH1D', 'CBUSH2D',

            'CDAMP1', 'CDAMP2', 'CDAMP3', 'CDAMP4', 'CDAMP5',
            'CFAST',

            'CBAR', 'CROD', 'CTUBE', 'CBEAM', 'CBEAM3', 'CONROD', 'CBEND',
            'CTRIA3', 'CTRIA6', 'CTRIAR', 'CTRIAX', 'CTRIAX6',
            'CQUAD4', 'CQUAD8', 'CQUADR', 'CQUADX', 'CQUAD',
            'CTETRA', 'CPENTA', 'CHEXA',
            'CSHEAR', 'CVISC', 'CRAC2D', 'CRAC3D',
            'CGAP',

            # rigid elements
            'RBAR', 'RBAR1', 'RBE1', 'RBE2', 'RBE3',

            # properties
            'PMASS',
            'PELAS', 'PGAP', 'PFAST',
            'PBUSH', 'PBUSH1D',
            'PDAMP', 'PDAMP5', 'PDAMPT',
            'PROD', 'PBAR', 'PBARL', 'PBEAM', 'PTUBE', 'PBEND', 'PBCOMP',
            'PBEAML', # 'PBEAM3',

            'PSHELL', 'PCOMP', 'PCOMPG', 'PSHEAR',
            'PSOLID', 'PLSOLID', 'PVISC','PRAC2D','PRAC3D',

            # creep materials
            'CREEP',

            # materials
            'MAT1', 'MAT2', 'MAT3', 'MAT8', 'MAT9', 'MAT10', 'MATHP',
            #'MATT1', 'MATT2', 'MATT3', 'MATT4', 'MATT5', 'MATT8', 'MATT9',
            'MATS1',

            # thermal materials
            'MAT4', 'MAT5',

            # spc/mpc constraints
            'SPC', 'SPCADD', 'SPC1', 'SPCD', 'SPCAX',
            'MPC', 'MPCADD',
            'SUPORT', 'SUPORT1',

            # loads
            'LOAD', 'LSEQ', 'RANDPS',
            'DLOAD', 'SLOAD', 'TLOAD1', 'TLOAD2', 'RLOAD1', 'RLOAD2',
            'FORCE', 'FORCE1', 'FORCE2',
            'MOMENT', 'MOMENT1', 'MOMENT2',
            'GRAV', 'ACCEL1',
            'PLOAD', 'PLOAD1', 'PLOAD2', 'PLOAD4',
            'PLOADX1', 'RFORCE',

            # aero cards
            'AERO', 'AEROS', 'GUST', 'FLUTTER', 'FLFACT', 'MKAERO1', 'MKAERO2',
            'AEFACT', 'AELINK', 'AELIST', 'AEPARAM', 'AESTAT', 'AESURF',
            'CAERO1', 'CAERO2',  # 'CAERO3', 'CAERO4', 'CAERO5',
            'PAERO1', 'PAERO2',  # 'PAERO3', 'PAERO4', 'PAERO5',
            'SPLINE1', 'SPLINE2', 'SPLINE4', 'SPLINE5',
            #'SPLINE3', 'SPLINE6', 'SPLINE7',
            'TRIM',

            # coords
            'CORD1R', 'CORD1C', 'CORD1S',
            'CORD2R', 'CORD2C', 'CORD2S',

            # temperature cards
            'TEMP',  # 'TEMPD',
            'QBDY1', 'QBDY2', 'QBDY3', 'QHBDY',
            'CHBDYE', 'CHBDYG', 'CHBDYP',
            'PCONV', 'PCONVM', 'PHBDY',
            'RADBC', 'CONV',  # 'RADM',

            # dynamic cards
            'DAREA', 'NLPARM', 'TSTEP', 'TSTEPNL',

            # frequencies
            'FREQ', 'FREQ1', 'FREQ2',

            # direct matrix input cards
            'DMIG', 'DMIJ', 'DMIJI', 'DMIK', 'DMI',
            'DEQATN',

            # optimization cards
            'DCONSTR', 'DESVAR', 'DDVAL', 'DRESP1', 'DRESP2', 'DVPREL1', 'DVPREL2',
            'DOPTPRM', 'DVMREL1', 'DLINK', 'DRESP3',

            # sets
            'ASET', 'BSET', 'CSET', 'QSET',  # 'USET',
            'ASET1', 'BSET1', 'CSET1', 'QSET1',  # 'USET1',
            'SET1', 'SET3',

            # super-element sets
            'SESET',

            # tables
            #'DTABLE', 'TABLEHT', 'TABRNDG',
            'TABLED1', 'TABLED2', 'TABLED3',  # 'TABLED4',
            'TABLEM1', 'TABLEM2', 'TABLEM3', 'TABLEM4',
            'TABLES1', 'TABLEST',
            'TABRND1', 'TABRNDG',

            # initial conditions - sid (set ID)
            #'TIC',  (in tables.py)

            ## methods - @todo EIGRL not done???
            'EIGB', 'EIGR', 'EIGRL',

            ## cMethods - @todo EIGC not done???
            'EIGC', 'EIGP',

            # other
            'INCLUDE',  # '='
            'ENDDATA',
            ])

        caseControlCards = set(['FREQ', 'GUST', 'MPC', 'SPC', 'NLPARM', 'NSM',
                            'TEMP', 'TSTEPNL', 'INCLUDE'])
        self.uniqueBulkDataCards = self.cardsToRead.difference(
            caseControlCards)

        ## / is the delete from restart card
        self.specialCards = ['DEQATN', '/']
        ## was playing around with an idea...does nothing for now...
        self.cardsToWrite = self.cardsToRead

    def disable_cards(self, cards):
        """
        Method for removing broken cards from the reader
        @param self the object pointer
        @param cards a list/set of cards that should not be read
        """
        disableSet = set(cards)
        self.cardsToRead.difference(disableSet)

    def is_thermal_solution(self):
        """
        @todo implement case control deck checker
        @warning dont use this...returns False
        """
        return False
        #self.caseControlDeck.is_thermal_solution()

    def _init_solution(self):
        """
        creates storage objects for the BDF object
        this would be in the init but doing it this way allows for
        better inheritance

        References:
          1.  http://www.mscsoftware.com/support/library/conf/wuc87/p02387.pdf
        """
        self.bdf_filename = None
        self._auto_reject = False
        self._solmap_to_value = {
            'NONLIN': 101,  # 66 -> 101 per Reference 1
            'SESTATIC': 101,
            'SESTATICS': 101,
            'SEMODES': 103,
            'BUCKLING': 105,
            'SEBUCKL': 105,
            'NLSTATIC': 106,
            'SEDCEIG': 107,
            'SEDFREQ': 108,
            'SEDTRAN': 109,
            'SEMCEIG': 110,
            'SEMFREQ': 111,
            'SEMTRAN': 112,
            'CYCSTATX': 114,
            'CYCMODE': 115,
            'CYCBUCKL': 116,
            'CYCFREQ': 118,
            'NLTRAN': 129,
            'AESTAT': 144,
            'FLUTTR': 145,
            'SEAERO': 146,
            'NLSCSH': 153,
            'NLTCSH': 159,
            'DBTRANS': 190,
            'DESOPT': 200,

            # guessing
            #'CTRAN'     : 115,
            'CFREQ': 118,

            # solution 200 names
            'STATICS': 101,
            'MODES': 103,
            'BUCK': 105,
            'DFREQ': 108,
            'MFREQ': 111,
            'MTRAN': 112,
            'DCEIG': 107,
            'MCEIG': 110,
            #'HEAT'     : None,
            #'STRUCTURE': None,
            #'DIVERGE'  : None,
            'FLUTTER': 145,
            'SAERO': 146,
            }

        self.rsolmap_toStr = {
            66: 'NONLIN',
            101: 'SESTSTATIC',  # linear static
            103: 'SEMODES',  # modal
            105: 'BUCKLING',  # buckling
            106: 'NLSTATIC',  # non-linear static
            107: 'SEDCEIG',  # direct complex frequency response
            108: 'SEDFREQ',  # direct frequency response
            109: 'SEDTRAN',  # direct transient response
            110: 'SEMCEIG',  # modal complex eigenvalue
            111: 'SEMFREQ',  # modal frequency response
            112: 'SEMTRAN',  # modal transient response
            114: 'CYCSTATX',
            115: 'CYCMODE',
            116: 'CYCBUCKL',
            118: 'CYCFREQ',
            129: 'NLTRAN',  # nonlinear transient
            144: 'AESTAT',  # static aeroelastic
            145: 'FLUTTR',  # flutter/aeroservoelastic
            146: 'SEAERO',  # dynamic aeroelastic
            153: 'NLSCSH',  # nonlinear static thermal
            159: 'NLTCSH',  # nonlinear transient thermal
            190: 'DBTRANS',
            200: 'DESOPT',  # optimization
        }
        self._init_structural_defaults()
        self._init_aero_defaults()
        self._init_thermal_defaults()

    def _is_special_card(self, cardName):
        """these cards are listed in the case control and the bulk data deck"""
        if cardName in self.specialCards:
            return True
        return False

    def _init_structural_defaults(self):
        """
        creates storage objects for the BDF object
        this would be in the init but doing it this way allows for
        better inheritance
        """
        ## the analysis type
        self.sol = None
        ## used in solution 600, method
        self.solMethod = None
        ## the line with SOL on it, marks ???
        self.iSolLine = None
        self.caseControlDeck = CaseControlDeck([], self.log)
        #self.executive_control_lines = [self.sol]

        # main structural block
        ## store the PARAM cards
        self.params = {}
        ## stores SPOINT, GRID cards
        self.nodes = {}
        ## stores POINT cards
        self.points = {}
        
        self.spoints = None
        ## stores GRIDSET card
        self.gridSet = None
        ## stores elements (CQUAD4, CTRIA3, CHEXA8, CTETRA4, CROD, CONROD,
        ## etc.)
        self.elements = {}
        ## stores rigid elements (RBE2, RBE3, RJOINT, etc.)
        self.rigidElements = {}
        ## store CMASS1,CMASS2,CMASS3,CMASS4,CMASS5
        #self.massElements = {}
        ## stores LOTS of propeties (PBAR, PBEAM, PSHELL, PCOMP, etc.)
        self.properties = {}
        ## stores MAT1, MAT2, MAT3,...MAT10 (no MAT4, MAT5)
        self.materials = {}
        ## stores MATS1
        self.materialDeps = {}
        ## stores the CREEP card
        self.creepMaterials = {}

        # loads
        ## stores LOAD, FORCE, MOMENT, etc.
        self.loads = {}
        #self.gusts  = {} # Case Control GUST = 100
        #self.random = {} # Case Control RANDOM = 100

        ## stores coordinate systems
        self.coords = {0: CORD2R()}

        # constraints
        ## stores SUPORT1s
        #self.constraints = {} # suport1, anything else???
        self.suports = []  # suport, suport1

        ## stores SPCADD,SPC,SPC1,SPCD,SPCAX
        self.spcObject2 = constraintObject2()
        ## stores MPCADD,MPC
        self.mpcObject2 = constraintObject2()

        self.spcs = {}
        self.spcadds = {}

        self.mpcs = {}
        self.mpcadds = {}

        # dynamic cards
        ## stores DAREA
        self.dareas = {}
        ## stores NLPARM
        self.nlparms = {}
        ## stores TSTEPs
        self.tsteps = {}
        ## stores TSTEPNL
        self.tstepnls = {}

        ## direct matrix input - DMIG
        self.dmis = {}
        self.dmigs = {}
        self.dmijs = {}
        self.dmijis = {}
        self.dmiks = {}
        self.dequations = {}

        ## frequencies
        self.frequencies = {}

        # optimization
        self.dconstrs = {}
        self.desvars = {}
        self.ddvals = {}
        self.dlinks = {}
        self.dresps = {}
        ## stores DVPREL1, DVPREL2...might change to DVxRel
        self.dvprels = {}
        self.dvmrels = {}
        self.doptprm = None

        ## SETx
        self.sets = {}
        self.asets = []
        self.bsets = []
        self.csets = []
        self.qsets = []
        ## SESETx
        self.setsSuper = {}

        ## tables
        self.tables = {}
        ## randomTables
        self.randomTables = {}

        ## EIGB, EIGR, EIGRL methods
        self.methods = {}
        # EIGC, EIGP methods
        self.cMethods = {}

    def _init_aero_defaults(self):
        """
        creates storage objects for the BDF object
        this would be in the init but doing it this way allows for
        better inheritance
        """
        # aero cards
        ## stores CAEROx
        self.caeros = {}
        ## stores PAEROx
        self.paeros = {}
        ## stores AERO
        self.aero = {}
        ## stores AEROS
        self.aeros = {}

        ## stores AEFACT
        self.aefacts = {}
        ## stores AELINK
        self.aelinks = {}
        ## stores AELIST
        self.aelists = {}
        ## stores AEPARAM
        self.aeparams = {}
        ## stores AESURF
        self.aesurfs = {}
        ## stores AESTAT
        self.aestats = {}

        ## stores GUST cards
        self.gusts = {}
        ## stores FLFACT
        self.flfacts = {}  # can this be simplified ???
        ## stores FLUTTER
        self.flutters = {}
        ## mkaeros
        self.mkaeros = []
        ## store SPLINE1,SPLINE2,SPLINE4,SPLINE5
        self.splines = {}
        ## stores TRIM
        self.trims = {}

    def _init_thermal_defaults(self):
        """initializes some bdf parameters"""
        # BCs
        ## stores thermal boundary conditions - CONV,RADBC
        self.bcs = {}  # e.g. RADBC
        ## defines the MAT4, MAT5, MATT4, etc. @todo verify MATT4
        self.thermalMaterials = {}

        # elements
        # see self.elements

        # properties
        ## stores other thermal properties - unused ???
        #self.thermalProperties    = {}
        ## stores PHBDY
        self.phbdys = {}
        ## stores convection properties - PCONV, PCONVM ???
        self.convectionProperties = {}

    def readBDF(self, infilename, includeDir=None, xref=True):
        """
        main read method for the bdf
        @param infilename
          the input bdf
        @param includeDir
          the relative path to any include files (default=None
          if no include files)
        @param xref
          should the bdf be cross referenced (default=True)
        """
        self._set_infile(infilename, includeDir)

        fname = self.print_filename(self.bdf_filename)
        self.log.debug('---starting BDF.readBDF of %s---' % fname)
        #self.log.info('xref=%s' %(xref))
        sys.stdout.flush()

        #self.debug = True
        if self.debug:
            self.log.debug("*BDF.readBDF")
        self._read_executive_control_deck()
        self._read_case_control_deck(self.bdf_filename)
        self._read_bulk_data_deck()

        self.cross_reference(xref=xref)
        if self.debug:
            self.log.debug("***BDF.readBDF")

        self.log.debug('---finished BDF.readBDF of %s---' % fname)
        sys.stdout.flush()

    def readBDF_Punch(self, infilename, includeDir=None, xref=True):
        """
        BDF punch file reader
        @param infilename
          the input bdf
        @param includeDir
          the relative path to any include files (default=None
          if no include files)
        @param xref
          should the bdf be cross referenced (default=True)
        """
        self._set_infile(infilename, includeDir)

        fname = self.print_filename(self.bdf_filename)
        self.log.debug('---starting BDF.readBDF_Punch of %s---' % fname)
        #self.log.info('xref=%s' %(xref))
        sys.stdout.flush()

        #self.debug = True
        if self.debug:
            self.log.debug("*BDF.readBDF_Punch")
        self._read_bulk_data_deck()
        #self.close_file()
        self.cross_reference(xref=xref)
        if self.debug:
            self.log.debug("***BDF.readBDF_Punch")

        self.log.debug('---finished BDF.readBDF_Punch of %s---' % fname)
        sys.stdout.flush()

    def _is_executive_control_deck(self, line):
        """@todo code this..."""
        return True

    def _read_executive_control_deck(self):
        """reads the executive control deck"""
        self.open_file(self.bdf_filename)
        line = ''
        emptyLines = 0
        while emptyLines < 50 and 'CEND' not in line.upper()[:4]:
            line = self.infilesPack[-1].readline()
            line = line.rstrip('\n\r\t ')
            if len(line) > 0:
                self.executive_control_lines.append(line)
            else:
                emptyLines += 1

        if 0:  # old method; breaks DMAP alters
            while len(self._active_filenames) > 0:  # keep going until finished
                lineIn = self.get_next_line()
                #print(lineIn)
                if lineIn is None:  # file was closed and a 2nd readCaseControl
                    return          # was called

                line = lineIn.strip()
                if self.debug:
                    (n) = self.get_line_number()
                    self.log.debug("executiveLine[%s] = |%s|" % (n, line))
                self.executive_control_lines.append(lineIn)
                if 'CEND' in line.upper():
                    break

        self._parse_executive_control_deck()

    def _parse_executive_control_deck(self):
        """extracts the solution from the executive control deck"""
        for (i, eline) in enumerate(self.executive_control_lines):
            #print 'eLine = |%r|' %(eline)
            uline = eline.strip().upper()  # uppercase line
            uline = uline.split('$')[0]
            if 'SOL ' in uline[:4]:
                #print "line = ",uline
                if ',' in uline:
                    sline = uline.split(',')  # SOL 600,method
                    solValue = sline[0]
                    method = sline[1]

                    #print "sline    = |%s|" %(sline)
                    #print "sline2   = |%s|" %(sline2)
                else:
                    solValue = uline
                    method = None

                #print "solValue = |%s|" %(solValue)
                sol = solValue[3:].strip()

                assert self.sol is None, ('cannot overwrite solution existing='
                                          '|SOL %s| new =|%s|' % (self.sol, uline))
                self.iSolLine = i

                try:
                    self.update_solution(sol, method)
                except:
                    msg = ('update_solution failed...line=%s' % uline)
                    raise RuntimeError(msg)

        #print("sol = ", sol)

    def updateSolution(self, sol, method=None):
        """
        @see update_solution
        """
        warnings.warn('updateSolution has been deprecated; use '
                      'update_solution', DeprecationWarning, stacklevel=2)
        self.update_solution(sol, method)

    def update_solution(self, sol, method=None):
        """
        updates the overall solution type (e.g. 101,200,600)
        @param self   the object pointer
        @param sol    the solution type (101,103, etc)
        @param method the solution method (only for SOL=600), default=None
        """
        ## the integer of the solution type (e.g. SOL 101)
        try:
            self.sol = int(sol)
            #print "sol = |%s|" %(sol)
        except ValueError:
            #print "sol = |%r|" %(sol)
            self.sol = self._solmap_to_value[sol]
            #print "sol = ",self.sol

        if self.sol == 600:
            ## solution 600 method modifier
            self.solMethod = method.strip()
            self.log.debug("sol=%s method=%s" % (self.sol, self.solMethod))
        else:  # very common
            self.solMethod = None

        #print "sol=%s method=%s" %(self.sol,self.solMethod)

    def setDynamicSyntax(self, dictOfVars):
        """@see set_dynamic_syntax"""
        warnings.warn('setDynamicSyntax has been deprecated; use '
                      'set_dynamic_syntax', DeprecationWarning, stacklevel=2)
        self.set_dynamic_syntax(dictOfVars)

    def set_dynamic_syntax(self, dictOfVars):
        """
        uses the OpenMDAO syntax of %varName in an embedded BDF to
        update the values for an optimization study.
        Variables should be 7 characters to fit in an 8-character field.
        %varName
        dictOfVars = {'varName': 10}
        """
        self.dictOfVars = {}
        for (key, value) in sorted(dictOfVars.iteritems()):
            assert len(key) <= 7, ('max length for key is 7; '
                                   'len(%s)=%s' % (key, len(key)))
            self.dictOfVars[key.upper()] = value

        self._is_dynamic_syntax = True

    def _is_case_control_deck(self, line):
        """@todo not done..."""
        #print "line = |%r|" %(line)
        lineUpper = line.upper().strip()
        #print "line = |%s|" %(lineUpper)
        if 'CEND' in line.upper():
            raise SyntaxError('invalid Case Control Deck card...CEND...')
        if '=' in lineUpper or ' ' in lineUpper:
            #print "True1"
            return True
        for card in self.uniqueBulkDataCards:
            lenCard = len(card)
            if card in lineUpper[:lenCard]:
                #print "*card = |%s|" %(card)
                #print "False1"
                return False
            #print "card = |%s|" %(card)
        #print "True2"
        return True

    def _read_case_control_deck(self, infilename):
        """
        reads the case control deck
        @note called with recursion if an INCLUDE file is found
        """
        #print "opening |%s|" %(infilename)
        self.open_file(infilename)
        #self.log.info("reading Case Control Deck...")
        line = ''
        #self.caseControlControlLines = []

        i = 0
        while len(self._active_filenames) > 0:  # keep going until finished
        #while 'BEGIN BULK' not in line:
            lineIn = self.get_next_line()
            #print("lineIn = |%r|" % (lineIn))
            #print("lineIn = ", lineIn)
            if lineIn is None:
                return  # file was closed and a 2nd readCaseControl was called
            if not self._is_case_control_deck(lineIn):
                self.linesPack = [lineIn] + self.linesPack
            line = lineIn.strip().split('$')[0].strip()
            lineUpper = line.upper()

            (line, lineUpper) = self._checkForIncludeFile_CaseControlDeck(
                lineIn, line, lineUpper)

            #print "*line = |%s|" %(line)
            if 'BEGIN' in lineUpper and 'BULK' in lineUpper:
                self.log.debug('found the end of the case control deck!')
                #print "breaking"
                break
            if i > 10000:
                raise RuntimeError('there are too many lines in the '
                                   'Case Control Deck; max=10000')
            i += 1
        #self.log.info("finished with Case Control Deck..")
        #print("self.case_controlLines = ",self.case_control_lines)

        #for line in self.case_control_lines:
            #print "** line=|%r|" %(line)

        self.caseControlDeck = CaseControlDeck(self.case_control_lines,
                                               self.log)
        self.caseControlDeck.solmap_toValue = self._solmap_to_value
        self.caseControlDeck.rsolmap_toStr = self.rsolmap_toStr

        #print "done w/ case control..."
        #print '***********************'
        return self.case_control_lines

    def _checkForIncludeFile_CaseControlDeck(self, lineIn, line, lineUpper):
        """
        Special parsing must be done if an INCLUDE file is found.
        @param lineIn
          the current line without any stripping of comments (defined by $)
        @param line
          the current line
        @param lineUpper
          the current line (in all caps)
        @retval line
          the new current line
        @retval lineUpper
          the new current line (in all caps)
        """
        if lineUpper.startswith('INCLUDE'):
            nextLine = self.get_next_line().strip().split('$')[0].strip()
            includeLines = [line]
            #print "^&*1",nextLine
            while '\\' in nextLine or '/' in nextLine:  # more includes
                includeLines.append(nextLine)
                nextLine = self.get_next_line().strip().split('$')[0].strip()
                #print "^&*2",nextLine

            #print "include lines = |%s|" %(includeLines)
            filename = self._get_include_file_name(includeLines)

            self._add_include_file(filename)
            #self.open_file(filename)
            self._read_case_control_deck(filename)
            line = nextLine
            #print "appending |%r|" %(nextLine)
            self.case_control_lines.append(nextLine)
        else:
            #print "appending |%r|" %(lineIn)
            self.case_control_lines.append(lineIn)

        return (line, lineUpper)

    def _get_card_name(self, cardLines):
        """
        Given a list of card lines, determines the cardName.
        @param self      the object pointer
        @param cardLines the list of lines that define the card
        @retval cardName the name of the card
        @note
            Parses the first 8 characters of a card, splits bassed on a comma,
            pulls off any spaces or asterisks and returns what is left.
        """
        #self.log.debug("getting cardName...")
        cardName = cardLines[0][0:8].strip()
        if ',' in cardName:
            cardName = cardName.split(',')[0].strip()

        cardName = cardName.lstrip().rstrip(' *')
        #self.log.debug("_get_card_name cardName=|%s|" %(cardName))
        return cardName

    def _is_reject(self, cardName):
        """Can the card be read"""
        #cardName = self._get_card_name(card)
        if cardName.startswith('='):
            return False
        elif cardName in self.cardsToRead:
            #print "*card = ",card
            #print "RcardName = |%s|" %(cardName)
            return False
        if cardName:
            if cardName not in self.rejectCount:
                self.log.info("RejectCardName = |%s|" % cardName)
                self.rejectCount[cardName] = 0
            self.rejectCount[cardName] += 1
        return True

    def _get_include_file_name(self, cardLines):
        """Parses an INCLUDE file split into multiple lines (as a list).
        @param self
          the BDF object
        @param cardLines
          the list of lines in the include card (all the lines!)
        @param cardName
          INCLUDE or include (needed to strip it off without converting
          the case)
        @retval filename the INCLUDE filename
        """
        cardLines2 = []
        for line in cardLines:
            line = line.strip('\t\r\n ')
            cardLines2.append(line)

        cardLines2[0] = cardLines2[0][7:].strip()  # truncate the cardname
        filename = ''.join(cardLines2)
        filename = filename.strip('"').strip("'")
        if ':' in filename:
            ifilepaths = filename.split(':')
            filename = os.path.join(*ifilepaths)
        #print 'filename = |%s|' % (filename)
        filename = os.path.join(self.includeDir, filename)
        return filename

    def _add_include_file(self, infileName):
        """
        This method must be called before opening an INCLUDE file.
        Identifies the new file as being opened.
        @param self the object pointer
        @param infileName the new INCLUDE file
        @note isOpened[fileName] is really initialized to False
        """
        self.isOpened[infileName] = False

    def _read_bulk_data_deck(self):
        """parses the Bulk Data Deck"""
        if self.debug:
            self.log.debug("*read_bulk_data_deck")
        self.open_file(self.bdf_filename)

        #old_card_obj = BDFCard()
        while len(self._active_filenames) > 0:  # keep going until finished
            ## gets the cardLines
            (rawCard, card, cardName) = self._get_card(debug=False)
            #print "outcard = ",card

            if cardName == 'INCLUDE':
                #print "rawCard = ",rawCard
                #print "card    = ",card
                filename = self._get_include_file_name(rawCard)
                #print 'filename = ', os.path.relpath(filename)
                self._add_include_file(filename)
                self.open_file(filename)
                reject = '$ INCLUDE processed:  %s\n' % (filename)
                self.rejects.append([reject])
                continue
            #elif cardName is None:
                #self.close_file()
                #continue

            passCard = False
            if self._is_special_card(cardName):
                passCard = True
                card = rawCard
                #print(rawCard)
            elif not self._is_reject(cardName):
                #print ""
                #print "not a reject"
                card = self.processCard(card)  # parse the card into fields
                #print "processedCard = ",card
            elif card[0].strip() == '':
                #print "funny strip thing..."
                pass
            else:
                #print "reject!"
                self.rejects.append(card)
                continue
                #print " rejecting card = ",card
                #card = self.processCard(card)

            #print "card2 = ",list_print(card)
            #print "card = ",card
            #cardName = self._get_card_name(card)

            if 'ENDDATA' in cardName:
                #print cardName
                break  # exits while loop
            #self.log.debug('cardName = |%s|' %(cardName))

            #cardObj = BDFCard(card,oldCardObj)
            #if cardName in self.specialCards:
            #    cardObj = card
            #else:
            #    cardObj = BDFCard(card)
            #cardObj = BDFCard(card)

            nCards = 1
            #special = False
            if '=' in cardName:
                nCards = cardName.strip('=()')
                if nCards:
                    nCards = int(nCards)
                else:
                    nCards = 1
                    #special = True
                #print "nCards = ",nCards
                #cardName = oldCardObj.field(0)

            for iCard in xrange(nCards):
                #if special:
                    #print "iCard = ",iCard
                self.add_card(card, cardName, iCard=0, old_card_obj=None)
                #if self.foundEndData:
                #    break

            if self.doneReading or len(self.linesPack[-1]) == 0:
                #print("doneReading=%s len(pack)=%s"
                #    %(self.doneReading, len(self.linesPack[-1])))
                self.close_file()

            #oldCardObj = copy.deepcopy(cardObj) # used for =(*1) stuff
            #print ""

            #print "self.linesPack[-1] = ",len(self.linesPack[-1])
            #print "self._active_filenames = ",self._active_filenames
        ### end of while loop

        #self.debug = True
        if self.debug:
            #for nid,node in self.nodes.iteritems():
                #print node
            #for eid,element in self.elements.iteritems():
                #print element

            #self.log.debug("\n$REJECTS")
            #for reject in self.rejects:
                #print printCard(reject)
                #print ''.join(reject)
            self.log.debug("***read_bulk_data_deck")


    def addCard(self, card, cardName, iCard=0, oldCardObj=None):
        """@see add_card"""
        warnings.warn('addCard has been deprecated; use add_card',
                      DeprecationWarning, stacklevel=2)
        return self.add_card(card, cardName, iCard=iCard)

    def add_card(self, card, cardName, iCard=0, old_card_obj=None):
        """
        adds a card object to the BDF object.
        @param self
          the object pointer
        @param card
          the list of the card fields -> ['GRID',1,2,]
        @param cardName
          the cardName -> 'GRID'
        @param iCard
          used when reading Nastran Free-Format (disabled)
        @param old_card_obj
          the last card object that was returned (type=BDFCard or None;
          default=None)
        @retval cardObject
          the card object representation of card
        @note
          this is a very useful method for interfacing with the code
        @note
           the cardObject is not a card-type object...so not a GRID card
           or CQUAD4 object.  It's a BDFCard Object.  However, you know the
           type (assuming a GRID), so just call the mesh.Node(nid) to get the
           No de object that was just created.
        @warning
          cardObject is not returned
        @todo this method is 600+ lines long...refactor time...
        """
        #print "card = ",card

        #if cardName in self.specialCards:
        #    pass #cardObj = card
        #else:
        cardObj = BDFCard(card, oldCardObj=None)

        #cardObj = BDFCard(card, old_card_obj=None)
        #if self.debug:
        #    self.log.debug(card)
        #    self.log.debug("*old_card_obj = \n%s" %(old_card_obj))
        #    self.log.debug("*cardObj = \n%s" %(cardObj))
        #cardObj.applyOldFields(iCard)

        try:
            if self._auto_reject == True:
                print('rejecting processed %s' % (card))
                self.reject_cards.append(card)
            elif card == [] or cardName == '':
                pass
            elif cardName == 'DMIG':  # not done...
                if cardObj.field(2) == 'UACCEL':  # special DMIG card
                    self.reject_cards.append(card)
                elif cardObj.field(2) == 0:
                    dmig = DMIG(cardObj)
                    self.addDMIG(dmig)
                else:
                    name = cardObj.field(1)
                    dmig = self.dmigs[name]
                    dmig.addColumn(cardObj)

            elif cardName == 'DMIJ':
                if cardObj.field(2) == 0:
                    dmij = DMIJ(cardObj)
                    self.addDMIJ(dmij)
                else:
                    name = cardObj.field(1)
                    dmij = self.dmijs[name]
                    dmij.addColumn(cardObj)
            elif cardName == 'DMIJI':
                if cardObj.field(2) == 0:
                    dmiji = DMIJI(cardObj)
                    self.addDMIJI(dmiji)
                else:
                    name = cardObj.field(1)
                    dmiji = self.dmijis[name]
                    dmiji.addColumn(cardObj)
            elif cardName == 'DMIK':
                if cardObj.field(2) == 0:
                    dmik = DMIK(cardObj)
                    self.addDMIK(dmik)
                else:
                    name = cardObj.field(1)
                    dmik = self.dmiks[name]
                    dmik.addColumn(cardObj)
            elif cardName == 'DMI':
                if cardObj.field(2) == 0:
                    dmi = DMI(cardObj)
                    self.addDMI(dmi)
                else:
                    name = cardObj.field(1)
                    dmi = self.dmis[name]
                    dmi.addColumn(cardObj)

            elif cardName == 'DEQATN':  # buggy for commas
                #print 'DEQATN:  cardObj.card=%s' %(cardObj.card)
                #equation = DEQATN(cardObj)
                #self.addDEQATN(equation)
                self.rejects.append(card)

            elif cardName == 'GRID':
                node = GRID(cardObj)
                self.addNode(node)
            elif cardName == 'GRDSET':
                self.gridSet = GRDSET(cardObj)
            #elif cardName == 'POINT':
            #    point = POINT(cardObj)
            #    self.addPoint(point)

            #elif cardName == 'RINGAX':
            #    node = RINGAX(cardObj)
            #    self.addNode(node)
            elif cardName == 'SPOINT':
                spoint = SPOINTs(cardObj)
                self.addSPoint(spoint)

            elif cardName == 'CQUAD4':
                elem = CQUAD4(cardObj)
                self.addElement(elem)
            elif cardName == 'CQUAD8':
                elem = CQUAD8(cardObj)
                self.addElement(elem)
            elif cardName == 'CQUAD':
                elem = CQUAD(cardObj)
                self.addElement(elem)
            elif cardName == 'CQUADR':
                elem = CQUADR(cardObj)
                self.addElement(elem)
            elif cardName == 'CQUADX':
                elem = CQUADX(cardObj)
                self.addElement(elem)

            elif cardName == 'CTRIA3':
                elem = CTRIA3(cardObj)
                self.addElement(elem)
            elif cardName == 'CTRIA6':
                elem = CTRIA6(cardObj)
                self.addElement(elem)

            elif cardName == 'CTRIAR':
                elem = CTRIAR(cardObj)
                self.addElement(elem)
            elif cardName == 'CTRIAX':
                elem = CTRIAX(cardObj)
                self.addElement(elem)
            elif cardName == 'CTRIAX6':
                elem = CTRIAX6(cardObj)
                self.addElement(elem)

            elif cardName == 'CTETRA':  # 4/10 nodes
                nFields = cardObj.nFields()
                if   nFields == 7:
                    elem = CTETRA4(cardObj)  # 4+3
                else:
                    elem = CTETRA10(cardObj)  # 10+3
                self.addElement(elem)
            elif cardName == 'CHEXA':  # 8/20 nodes
                nFields = cardObj.nFields()
                if   nFields == 11:
                    elem = CHEXA8(cardObj)   # 8+3
                else:
                    elem = CHEXA20(cardObj)  # 20+3
                self.addElement(elem)
            elif cardName == 'CPENTA':  # 6/15 nodes
                nFields = cardObj.nFields()
                if   nFields == 9:
                    elem = CPENTA6(cardObj)   # 6+3
                else:
                    elem = CPENTA15(cardObj)  # 15+3
                self.addElement(elem)

            elif cardName == 'CBAR':
                elem = CBAR(cardObj)
                self.addElement(elem)
            elif cardName == 'CBEAM':
                elem = CBEAM(cardObj)
                self.addElement(elem)
            elif cardName == 'CBEAM3':
                elem = CBEAM3(cardObj)
                self.addElement(elem)
            elif cardName == 'CROD':
                elem = CROD(cardObj)
                self.addElement(elem)
            elif cardName == 'CONROD':
                elem = CONROD(cardObj)
                self.addElement(elem)
            elif cardName == 'CTUBE':
                elem = CTUBE(cardObj)
                self.addElement(elem)
            elif cardName == 'CBEND':
                elem = CBEND(cardObj)
                self.addElement(elem)

            elif cardName == 'CELAS1':
                elem = CELAS1(cardObj)
                self.addElement(elem)
            elif cardName == 'CELAS2':
                (elem) = CELAS2(cardObj)
                self.addElement(elem)
            elif cardName == 'CELAS3':
                (elem) = CELAS3(cardObj)
                self.addElement(elem)
            elif cardName == 'CELAS4':
                (elem) = CELAS4(cardObj)
                self.addElement(elem)

            elif cardName == 'CBUSH':
                (elem) = CBUSH(cardObj)
                self.addDamperElement(elem)
            elif cardName == 'CBUSH1D':
                (elem) = CBUSH1D(cardObj)
                self.addDamperElement(elem)
            elif cardName == 'CFAST':
                (elem) = CFAST(cardObj)
                self.addDamperElement(elem)

            elif cardName == 'CDAMP1':
                (elem) = CDAMP1(cardObj)
                self.addDamperElement(elem)
            elif cardName == 'CDAMP2':
                (elem) = CDAMP2(cardObj)
                self.addDamperElement(elem)
            elif cardName == 'CDAMP3':
                (elem) = CDAMP3(cardObj)
                self.addDamperElement(elem)
            elif cardName == 'CDAMP4':
                (elem) = CDAMP4(cardObj)
                self.addDamperElement(elem)
            elif cardName == 'CDAMP5':
                (elem) = CDAMP5(cardObj)
                self.addDamperElement(elem)

            elif cardName == 'CONM1':
                elem = CONM1(cardObj)
                self.addElement(elem)
            elif cardName == 'CONM2':
                elem = CONM2(cardObj)
                self.addElement(elem)

            elif cardName == 'CMASS1':
                elem = CMASS1(cardObj)
                self.addElement(elem)
            elif cardName == 'CMASS2':
                elem = CMASS2(cardObj)
                self.addElement(elem)
            elif cardName == 'CMASS3':
                elem = CMASS3(cardObj)
                self.addElement(elem)
            elif cardName == 'CMASS4':
                elem = CMASS4(cardObj)
                self.addElement(elem)

            elif cardName == 'CVISC':
                elem = CVISC(cardObj)
                self.addElement(elem)
            elif cardName == 'CSHEAR':
                elem = CSHEAR(cardObj)
                self.addElement(elem)
            elif cardName == 'CGAP':
                elem = CGAP(cardObj)
                self.addElement(elem)

            elif cardName == 'CRAC2D':
                elem = CRAC2D(cardObj)
                self.addElement(elem)
            elif cardName == 'CRAC3D':
                elem = CRAC3D(cardObj)
                self.addElement(elem)

            # rigid elements
            elif cardName == 'RBAR':
                (elem) = RBAR(cardObj)
                self.addRigidElement(elem)
            elif cardName == 'RBAR1':
                (elem) = RBAR1(cardObj)
                self.addRigidElement(elem)
            elif cardName == 'RBE1':
                (elem) = RBE1(cardObj)
                self.addRigidElement(elem)
            elif cardName == 'RBE2':
                (elem) = RBE2(cardObj)
                self.addRigidElement(elem)
            elif cardName == 'RBE3':
                (elem) = RBE3(cardObj)
                self.addRigidElement(elem)

            elif cardName == 'PSHELL':
                prop = PSHELL(cardObj)
                self.addProperty(prop)
            elif cardName == 'PCOMP':
                prop = PCOMP(cardObj)
                self.addProperty(prop)
            elif cardName == 'PCOMPG':  # hasnt been verified
                prop = PCOMPG(cardObj)
                self.addProperty(prop)
            elif cardName == 'PSHEAR':
                prop = PSHEAR(cardObj)
                self.addProperty(prop)

            elif cardName == 'PSOLID':
                prop = PSOLID(cardObj)
                self.addProperty(prop)
            elif cardName == 'PBAR':
                prop = PBAR(cardObj)
                self.addProperty(prop)
            elif cardName == 'PBARL':
                prop = PBARL(cardObj)
                self.addProperty(prop)
            elif cardName == 'PBEAM':
                prop = PBEAM(cardObj)
                self.addProperty(prop)
            elif cardName == 'PBCOMP':
                prop = PBCOMP(cardObj)
                self.addProperty(prop)
            #elif cardName == 'PBEAM3':
            #    prop = PBEAM3(cardObj)
            #    self.addProperty(prop)
            elif cardName == 'PBEAML':   # disabled
                prop = PBEAML(cardObj)
                self.addProperty(prop)
            elif cardName == 'PELAS':
                prop = PELAS(cardObj)
                self.addProperty(prop)
                if cardObj.field(5):
                    prop = PELAS(cardObj, 1)  # makes 2nd PELAS card
                    self.addProperty(prop)
            elif cardName == 'PVISC':
                prop = PVISC(cardObj)
                self.addProperty(prop)
                if cardObj.field(5):
                    prop = PVISC(card, 1)
                    self.addProperty(prop)
            elif cardName == 'PROD':
                prop = PROD(cardObj)
                self.addProperty(prop)
            elif cardName == 'PTUBE':
                prop = PTUBE(cardObj)
                self.addProperty(prop)
            elif cardName == 'PMASS':
                prop = PMASS(cardObj, nOffset=0)
                self.addProperty(prop)

                if cardObj.field(3) is not None:
                    prop = PMASS(cardObj, nOffset=1)
                    self.addProperty(prop)
                if cardObj.field(5) is not None:
                    prop = PMASS(cardObj, nOffset=2)
                    self.addProperty(prop)
                if cardObj.field(7) is not None:
                    prop = PMASS(cardObj, nOffset=3)
                    self.addProperty(prop)

            elif cardName == 'PLSOLID':
                prop = PLSOLID(cardObj)
                self.addProperty(prop)

            elif cardName == 'PBUSH1D':
                prop = PBUSH1D(cardObj)
                self.addProperty(prop)
            #elif cardName == 'PBUSHT':
            #    prop = PBUSHT(cardObj)
            #    self.addProperty(prop)
            elif cardName == 'PBUSH':
                prop = PBUSH(cardObj)
                self.addProperty(prop)
            elif cardName == 'PFAST':
                prop = PFAST(cardObj)
                self.addProperty(prop)

            elif cardName == 'PDAMPT':
                prop = PDAMPT(cardObj)
                self.addProperty(prop)
            elif cardName == 'PDAMP':
                prop = PDAMP(cardObj)
                self.addProperty(prop)
                if cardObj.field(3):
                    prop = PDAMP(cardObj, 1)  # makes 2nd PDAMP card
                    self.addProperty(prop)
                if cardObj.field(5):
                    prop = PDAMP(cardObj, 1)  # makes 3rd PDAMP card
                    self.addProperty(prop)

            elif cardName == 'PDAMP5':
                prop = PDAMP5(cardObj)
                self.addProperty(prop)
            elif cardName == 'PGAP':
                elem = PGAP(cardObj)
                self.addProperty(elem)
            elif cardName == 'PRAC2D':
                prop = PRAC2D(cardObj)
                self.addProperty(prop)
            elif cardName == 'PRAC3D':
                prop = PRAC3D(cardObj)
                self.addProperty(prop)

            elif cardName == 'CREEP':  # hasnt been verified
                creep = CREEP(cardObj)
                # links up to MAT1, MAT2, MAT9 w/ same MID
                self.addCreepMaterial(creep)
            elif cardName == 'MAT1':
                material = MAT1(cardObj)
                self.addMaterial(material)
            elif cardName == 'MAT2':
                material = MAT2(cardObj)
                self.addMaterial(material)
            elif cardName == 'MAT3':
                material = MAT3(cardObj)
                self.addMaterial(material)
            elif cardName == 'MAT4':
                material = MAT4(cardObj)
                self.addThermalMaterial(material)
            elif cardName == 'MAT5':
                material = MAT5(cardObj)
                self.addThermalMaterial(material)
            elif cardName == 'MAT8':  # note there is no MAT6 or MAT7
                material = MAT8(cardObj)
                self.addMaterial(material)
            elif cardName == 'MAT9':
                material = MAT9(cardObj)
                self.addMaterial(material)
            elif cardName == 'MAT10':
                material = MAT10(cardObj)
                self.addMaterial(material)

            elif cardName == 'MATHP':
                material = MATHP(cardObj)
                self.addMaterial(material)

            elif cardName == 'MATS1':
                material = MATS1(cardObj)
                self.addMaterialDependence(material)
            #elif cardName == 'MATT1':
            #    material = MATT1(cardObj)
            #    self.addTempMaterial(material)
            #elif cardName == 'MATT2':
            #    material = MATT2(cardObj)
            #    self.addTempMaterial(material)
            #elif cardName == 'MATT3':
            #    material = MATT3(cardObj)
            #    self.addTempMaterial(material)
            #elif cardName == 'MATT4':
            #    material = MATT4(cardObj)
            #    self.addTempMaterial(material)
            #elif cardName == 'MATT5':
            #    material = MATT5(cardObj)
            #    self.addTempMaterial(material)
            #elif cardName == 'MATT8':
            #    material = MATT8(cardObj)
            #    self.addTempMaterial(material)
            #elif cardName == 'MATT9':
            #    material = MATT9(cardObj)
            #    self.addTempMaterial(material)

            # cards contained within the LOAD card
            elif cardName == 'FORCE':
                force = FORCE(cardObj)
                self.addLoad(force)
            elif cardName == 'FORCE1':
                force = FORCE1(cardObj)
                self.addLoad(force)
            elif cardName == 'FORCE2':
                force = FORCE2(cardObj)
                self.addLoad(force)
            elif cardName == 'MOMENT':
                moment = MOMENT(cardObj)
                self.addLoad(moment)
            elif cardName == 'MOMENT1':
                moment = MOMENT1(cardObj)
                self.addLoad(moment)
            elif cardName == 'MOMENT2':
                moment = MOMENT2(cardObj)
                self.addLoad(moment)
            elif cardName == 'GRAV':
                grav = GRAV(cardObj)
                self.addLoad(grav)
            elif cardName == 'ACCEL1':
                grav = ACCEL1(cardObj)
                self.addLoad(grav)
            elif cardName == 'LOAD':
                load = LOAD(cardObj)
                self.addLoad(load)
            elif cardName == 'PLOAD':
                load = PLOAD(cardObj)
                self.addLoad(load)
            elif cardName == 'PLOAD1':
                load = PLOAD1(cardObj)
                self.addLoad(load)
            elif cardName == 'PLOAD2':
                load = PLOAD2(cardObj)
                self.addLoad(load)
            elif cardName == 'PLOAD4':
                load = PLOAD4(cardObj)
                self.addLoad(load)
            elif cardName == 'PLOADX1':
                load = PLOADX1(cardObj)
                self.addLoad(load)
            elif cardName == 'RFORCE':
                load = RFORCE(cardObj)
                self.addLoad(load)

            elif cardName == 'LSEQ':
                load = LSEQ(cardObj)
                self.addLSeq(load)

            elif cardName == 'DLOAD':
                load = DLOAD(cardObj)
                self.addLoad(load)
            elif cardName == 'SLOAD':
                load = SLOAD(cardObj)
                self.addLoad(load)

            elif cardName == 'TLOAD1':
                load = TLOAD1(cardObj)
                self.addLoad(load)
            elif cardName == 'TLOAD2':
                load = TLOAD2(cardObj)
                self.addLoad(load)
            elif cardName == 'RLOAD1':
                load = RLOAD1(cardObj)
                self.addLoad(load)
            elif cardName == 'RLOAD2':
                load = RLOAD2(cardObj)
                self.addLoad(load)
            elif cardName == 'RANDPS':
                load = RANDPS(cardObj)
                self.addLoad(load)

            # thermal loads
            elif cardName == 'TEMP':
                load = TEMP(cardObj)
                self.addThermalLoad(load)
            #elif cardName == 'TEMPD':
            #    load = TEMPD(cardObj)
            #    self.addThermalLoad(load)
            elif cardName == 'QBDY1':
                load = QBDY1(cardObj)
                self.addThermalLoad(load)
            elif cardName == 'QBDY2':
                load = QBDY2(cardObj)
                self.addThermalLoad(load)
            elif cardName == 'QBDY3':
                load = QBDY3(cardObj)
                self.addThermalLoad(load)
            elif cardName == 'QHBDY':
                load = QHBDY(cardObj)
                self.addThermalLoad(load)

            # thermal elements
            elif cardName == 'CHBDYE':
                element = CHBDYE(cardObj)
                self.addThermalElement(element)
            elif cardName == 'CHBDYG':
                element = CHBDYG(cardObj)
                self.addThermalElement(element)
            elif cardName == 'CHBDYP':
                element = CHBDYP(cardObj)
                self.addThermalElement(element)

            # thermal properties
            elif cardName == 'PCONV':
                prop = PCONV(cardObj)
                self.addConvectionProperty(prop)
            elif cardName == 'PCONVM':
                prop = PCONVM(cardObj)
                self.addConvectionProperty(prop)
            elif cardName == 'PHBDY':
                prop = PHBDY(cardObj)
                self.addPHBDY(prop)

            # thermal BCs
            elif cardName == 'CONV':
                bc = CONV(cardObj)
                self.addThermalBC(bc, bc.eid)
            #elif cardName == 'RADM':
            #    bc = RADM(cardObj)
            #    self.addThermalBC(bc, bc.nodamb)
            elif cardName == 'RADBC':
                bc = RADBC(cardObj)
                self.addThermalBC(bc, bc.nodamb)

            #elif cardName == 'TABLEH1':
            #    load = TABLEH1(cardObj)
            #    self.addTable(load)

            # constraints
            elif cardName == 'MPC':
                constraint = MPC(cardObj)
                self.addConstraint_MPC(constraint)
            elif cardName == 'MPCADD':
                constraint = MPCADD(cardObj)
                #assert not isinstance(constraint,SPCADD)
                self.addConstraint_MPC(constraint)

            # constraints
            elif cardName == 'SPC':
                constraint = SPC(cardObj)
                self.addConstraint_SPC(constraint)
            elif cardName == 'SPC1':
                constraint = SPC1(cardObj)
                self.addConstraint_SPC(constraint)
            elif cardName == 'SPCAX':
                constraint = SPCAX(cardObj)
                self.addConstraint_SPC(constraint)
            elif cardName == 'SPCD':
                constraint = SPCD(cardObj)
                self.addConstraint_SPC(constraint)
            elif cardName == 'SPCADD':
                constraint = SPCADD(cardObj)
                #assert not isinstance(constraint,MPCADD)
                self.addConstraint_SPC(constraint)
            elif cardName == 'SUPORT':  # pseudo-constraint
                suport = SUPORT(cardObj)
                self.addSuport(suport)
            elif cardName == 'SUPORT1':  # pseudo-constraint
                suport1 = SUPORT1(cardObj)
                self.addConstraint(suport1)

            # aero
            elif cardName == 'SPLINE1':
                aero = SPLINE1(cardObj)
                self.addSpline(aero)
            elif cardName == 'SPLINE2':
                aero = SPLINE2(cardObj)
                self.addSpline(aero)
            #elif cardName == 'SPLINE3':
                #aero = SPLINE3(cardObj)
                #self.addSpline(aero)
            elif cardName == 'SPLINE4':
                aero = SPLINE4(cardObj)
                self.addSpline(aero)
            elif cardName == 'SPLINE5':
                aero = SPLINE5(cardObj)
                self.addSpline(aero)

            elif cardName == 'CAERO1':
                aero = CAERO1(cardObj)
                self.addCAero(aero)
            elif cardName == 'CAERO2':
                aero = CAERO2(cardObj)
                self.addCAero(aero)
            #elif cardName == 'CAERO3':
            #    aero = CAERO3(cardObj)
            #    self.addCAero(aero)
            elif cardName == 'PAERO1':
                aero = PAERO1(cardObj)
                self.addPAero(aero)
            elif cardName == 'PAERO2':
                aero = PAERO2(cardObj)
                self.addPAero(aero)
            #elif cardName == 'PAERO3':
            #    aero = PAERO3(cardObj)
            #    self.addPAero(aero)
            elif cardName == 'AERO':
                aero = AERO(cardObj)
                self.addAero(aero)
            elif cardName == 'AEROS':
                aeros = AEROS(cardObj)
                self.addAeros(aeros)

            elif cardName == 'AEFACT':
                aefact = AEFACT(cardObj)
                self.addAEFact(aefact)
            elif cardName == 'AELINK':
                aelink = AELINK(cardObj)
                self.addAELink(aelink)
            elif cardName == 'AELIST':
                aelist = AELIST(cardObj)
                self.addAEList(aelist)
            elif cardName == 'AEPARM':
                aeparm = AEPARM(cardObj)
                self.addAEParam(aeparm)
            elif cardName == 'AESTAT':
                aestat = AESTAT(cardObj)
                self.addAEStat(aestat)
            elif cardName == 'AESURF':
                aesurf = AESURF(cardObj)
                self.addAESurf(aesurf)

            elif cardName == 'TRIM':
                trim = TRIM(cardObj)
                self.addTrim(trim)

            elif cardName == 'FLUTTER':
                flutter = FLUTTER(cardObj)
                self.addFlutter(flutter)
            elif cardName == 'MKAERO1':
                mkaero = MKAERO1(cardObj)
                self.addMKAero(mkaero)
            elif cardName == 'MKAERO2':
                mkaero = MKAERO2(cardObj)
                self.addMKAero(mkaero)

            elif cardName == 'FLFACT':
                flfact = FLFACT(cardObj)
                self.addFLFACT(flfact)
            elif cardName == 'GUST':
                gust = GUST(cardObj)
                self.addGUST(gust)

            # dynamic
            elif cardName == 'DAREA':
                darea = DAREA(cardObj)
                self.addDArea(darea)
                if cardObj.field(5):
                    darea = DAREA(cardObj, 1)
                    self.addDArea(darea)
            elif cardName == 'NLPARM':
                nlparmObj = NLPARM(cardObj)
                self.addNLParm(nlparmObj)
            elif cardName == 'TSTEP':
                tstepObj = TSTEP(cardObj)
                self.addTSTEP(tstepObj)
            elif cardName == 'TSTEPNL':
                tstepnlObj = TSTEPNL(cardObj)
                self.addTSTEPNL(tstepnlObj)

            # frequencies
            elif cardName == 'FREQ':
                freq = FREQ(cardObj)
                self.addFREQ(freq)
            elif cardName == 'FREQ1':
                freq = FREQ1(cardObj)
                self.addFREQ(freq)
            elif cardName == 'FREQ2':
                freq = FREQ2(cardObj)
                self.addFREQ(freq)
            #elif cardName == 'FREQ3':
            #    freq = FREQ3(cardObj)
            #    self.addFREQ(freq)
            #elif cardName == 'FREQ4':
            #    freq = FREQ4(cardObj)
            #    self.addFREQ(freq)
            #elif cardName == 'FREQ5':
            #    freq = FREQ5(cardObj)
            #    self.addFREQ(freq)

            # SETx
            elif cardName == 'ASET':
                set_obj = ASET(cardObj)
                self.addASet(set_obj)
            elif cardName == 'BSET':
                set_obj = BSET(cardObj)
                self.addBSet(set_obj)
            elif cardName == 'CSET':
                set_obj = CSET(cardObj)
                self.addCSet(set_obj)
            elif cardName == 'QSET':
                set_obj = QSET(cardObj)
                self.addQSet(set_obj)
            elif cardName == 'ASET1':
                set_obj = ASET1(cardObj)
                self.addASet(set_obj)
            elif cardName == 'BSET1':
                set_obj = BSET1(cardObj)
                self.addBSet(set_obj)
            elif cardName == 'CSET1':
                set_obj = CSET1(cardObj)
                self.addCSet(set_obj)
            elif cardName == 'QSET1':
                set_obj = QSET1(cardObj)
                self.addQSet(set_obj)

            #elif cardName == 'USET':
            #    set_obj = USET(cardObj)
            #    self.addSet(set_obj)
            #elif cardName == 'USET1':
            #    set_obj = USET1(cardObj)
            #    self.addSet(set_obj)
            elif cardName == 'SET1':
                set_obj = SET1(cardObj)
                self.addSet(set_obj)
            #elif cardName == 'SET2':
            #    set_obj = SET2(cardObj)
            #    self.addSet(set_obj)
            elif cardName == 'SET3':
                set_obj = SET3(cardObj)
                self.addSet(set_obj)

            elif cardName == 'SESET':
                set_obj = SESET(cardObj)
                self.addSetSuper(set_obj)

            # optimization
            elif cardName == 'DCONSTR':
                flutter = DCONSTR(cardObj)
                self.addDConstr(flutter)
            elif cardName == 'DESVAR':
                desvar = DESVAR(cardObj)
                self.addDesvar(desvar)
            elif cardName == 'DDVAL':
                ddval = DDVAL(cardObj)
                self.addDDVal(ddval)
            elif cardName == 'DLINK':
                dlink = DLINK(cardObj)
                self.addDLink(dlink)
            elif cardName == 'DRESP1':
                ddval = DRESP1(cardObj)
                self.addDResp(ddval)
            elif cardName == 'DRESP2':
                ddval = DRESP2(cardObj)
                self.addDResp(ddval)
            elif cardName == 'DVPREL1':
                dvprel = DVPREL1(cardObj)
                self.addDvprel(dvprel)
            elif cardName == 'DVPREL2':
                dvprel = DVPREL2(cardObj)
                self.addDvprel(dvprel)
            elif cardName == 'DVMREL1':
                dvmrel = DVMREL1(cardObj)
                self.addDvmrel(dvmrel)
            #elif cardName == 'DVMREL2':
            #    dvmrel = DVMREL2(cardObj)
            #    self.addDvmrel(dvmrel)
            elif cardName == 'DOPTPRM':
                doptprm = DOPTPRM(cardObj)
                #self.addDoptprm(doptprm)
                self.doptprm = doptprm

            # coordinate systems
            elif cardName == 'CORD2R':
                coord = CORD2R(cardObj)
                self.addCoord(coord)
            elif cardName == 'CORD2C':
                coord = CORD2C(cardObj)
                self.addCoord(coord)
            elif cardName == 'CORD2S':
                coord = CORD2S(cardObj)
                self.addCoord(coord)

            elif cardName == 'CORD1R':
                coord = CORD1R(cardObj)
                self.addCoord(coord)
                if cardObj.field(5):
                    coord = CORD1R(cardObj, nCoord=1)
                    self.addCoord(coord)

            elif cardName == 'CORD1C':
                coord = CORD1C(cardObj)
                self.addCoord(coord)
                if cardObj.field(5):
                    coord = CORD1C(cardObj, nCoord=1)
                    self.addCoord(coord)

            elif cardName == 'CORD1S':
                coord = CORD1S(cardObj)
                self.addCoord(coord)
                if cardObj.field(5):
                    coord = CORD1S(cardObj, nCoord=1)
                    self.addCoord(coord)

            #elif cardName == 'CORD3G':
            #    coord = CORD3G(cardObj)
            #    self.addCoord(coord)

            # Table
            elif cardName == 'TABLED1':
                table = TABLED1(cardObj)
                self.addTable(table)
            elif cardName == 'TABLED2':
                table = TABLED2(cardObj)
                self.addTable(table)
            elif cardName == 'TABLED3':
                table = TABLED3(cardObj)
                self.addTable(table)
            #elif cardName == 'TABLED4':
            #    table = TABLED4(cardObj)
            #    self.addTable(table)

            elif cardName == 'TABLEM1':
                table = TABLEM1(cardObj)
                self.addTable(table)
            elif cardName == 'TABLEM2':
                table = TABLEM2(cardObj)
                self.addTable(table)
            elif cardName == 'TABLEM3':
                table = TABLEM3(cardObj)
                self.addTable(table)
            elif cardName == 'TABLEM4':
                table = TABLEM4(cardObj)
                self.addTable(table)

            elif cardName == 'TABLES1':
                table = TABLES1(cardObj)
                self.addTable(table)
            elif cardName == 'TABLEST':
                table = TABLEST(cardObj)
                self.addTable(table)

            # randomTables
            elif cardName == 'TABRND1':
                table = TABRND1(cardObj)
                self.addRandomTable(table)
            elif cardName == 'TABRNDG':
                table = TABRNDG(cardObj)
                self.addRandomTable(table)

            # Methods
            elif cardName == 'EIGB':
                method = EIGB(cardObj)
                self.addMethod(method)
            elif cardName == 'EIGR':
                method = EIGR(cardObj)
                self.addMethod(method)
            elif cardName == 'EIGRL':
                method = EIGRL(cardObj, sol=self.sol)
                self.addMethod(method)

            # CMethods
            elif cardName == 'EIGC':
                method = EIGC(cardObj)
                self.addCMethod(method)
            elif cardName == 'EIGP':
                method = EIGP(cardObj)
                self.addCMethod(method)

            elif cardName == 'PARAM':
                param = PARAM(cardObj)
                self.addParam(param)
            elif 'ENDDATA' in cardName:
                self.foundEndData = True
                #break
            else:
                ## @warning cards with = signs in them
                ## are not announced when they are rejected
                if '=' not in card[0]:
                    print('rejecting processed %s' % (card))
                self.reject_cards.append(card)
            ###
        except:
            print("cardName = |%r|" % (cardName))
            print("failed! Unreduced Card=%s\n" % (ListPrint(card)))
            print("filename = %s\n" % (self.bdf_filename))
            sys.stdout.flush()
            raise

        return cardObj
