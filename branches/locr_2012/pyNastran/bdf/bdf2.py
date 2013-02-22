# pylint: disable=C0103,C0302,R0901,R0902,R0904,R0912,R0914,R0915,W0201,W0611
"""
Main BDF class
"""
from __future__ import (nested_scopes, generators, division, absolute_import,
                        print_function, unicode_literals)


import sys
import os
#import io
import warnings

from pyNastran.utils import list_print
from pyNastran.utils import object_attributes
from pyNastran.utils.log import get_logger

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
from .bdfInterface.BDF_Card import BDFCard, wipe_empty_fields
#from .bdfInterface.bdf_Reader import BDFReader
from .bdfInterface.bdf_writeMesh import WriteMesh
from .bdfInterface.bdf_cardMethods import interpretValue # CardMethods
from .bdfInterface.crossReference import XrefMesh

class BDF(BDFMethods, GetMethods, AddMethods, WriteMesh, XrefMesh):
    """
    NASTRAN BDF Reader/Writer/Editor class.
    """
    modelType = 'nastran'

    def __init__(self, debug=True, log=None):
        """
        Initializes the BDF object
        @param self the BDF object
        @param debug used to set the logger if no logger is passed in
        @param log a python logging module object
        """
        debug = False
        self.active_file = None
        self.active_files = []
        self.active_filenames = []
        self.used_filenames = []
        self.stored_lines = {}
        self.stored_comments = {}
        ## list of all read in cards - useful in determining if entire BDF
        ## was read & really useful in debugging
        self.cardCount = {}
        ## stores the cardCount of cards that have been rejected
        self.reject_count = {}
        ## was an ENDDATA card found
        self.foundEndData = False

        #BDFReader.__init__(self, debug, log)
        self.relpath = True
        if sys.version_info < (2, 6):
            self.relpath = False
            #raise RuntimeError("must use python 2.6 or greater...version=%s"
            #                   %(str(sys.version_info)))
        self.log = get_logger(log, 'debug' if debug else 'info')

        ## allows the BDF variables to be scoped properly (i think...)
        GetMethods.__init__(self)
        AddMethods.__init__(self)
        BDFMethods.__init__(self)
        WriteMesh.__init__(self)
        XrefMesh.__init__(self)

        ## useful in debugging errors in input
        self.debug = debug
        self._init_solution()

        ## flag that allows for OpenMDAO-style optimization syntax to be used
        self._is_dynamic_syntax = False
        ## lines that were rejected b/c they were for a card that isnt supported
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
            #'PBEAML',  # not fully supported
            # 'PBEAM3',

            'PSHELL', 'PCOMP', 'PCOMPG', 'PSHEAR',
            'PSOLID', 'PLSOLID', 'PVISC', 'PRAC2D', 'PRAC3D',

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
            'DCONSTR', 'DESVAR', 'DDVAL', 'DRESP1', 'DRESP2',
            'DVPREL1', 'DVPREL2',
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

    def print_filename(self, filename):
        """
        Takes a path such as C:/work/fem.bdf and locates the file using
        relative paths.  If it's on another drive, the path is not modified.
        @param self the BDF object
        @param filename a filename string
        @retval filenameString a shortened representation of the filename
        """
        driveLetter = os.path.splitdrive(os.path.abspath(filename))[0]
        if driveLetter == os.path.splitdrive(os.curdir)[0] and self.relpath:
            return os.path.relpath(filename)
        return filename

    def disable_cards(self, cards):
        """
        Method for removing broken cards from the reader
        @param self the BDF object
        @param cards a list/set of cards that should not be read
        """
        disableSet = set(cards)
        self.cardsToRead.difference(disableSet)

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
        Creates storage objects for the BDF object
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

        self.pbusht = {}
        self.pdampt = {}
        self.pelast = {}

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
        Creates storage objects for the BDF object
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
        self.flfacts = {}  ## TODO can this be simplified ???
        ## stores FLUTTER
        self.flutters = {}
        ## mkaeros
        self.mkaeros = []
        ## store SPLINE1,SPLINE2,SPLINE4,SPLINE5
        self.splines = {}
        ## stores TRIM
        self.trims = {}

    def _init_thermal_defaults(self):
        """Initializes some bdf parameters"""
        # BCs
        ## stores thermal boundary conditions - CONV,RADBC
        self.bcs = {}  # e.g. RADBC
        ## defines the MAT4, MAT5, MATT4, etc. TODO verify MATT4
        self.thermalMaterials = {}

        ## stores PHBDY
        self.phbdys = {}
        ## stores convection properties - PCONV, PCONVM ???
        self.convectionProperties = {}

    def readBDF(self, bdf_filename, include_dir=None, xref=True, punch=False):
        """
        Read method for the bdf files
        @param self the BDF object
        @param bdf_filename the input bdf
        @param include_dir the relative path to any include files (default=None
          if no include files)
        @param xref should the bdf be cross referenced (default=True)
        @param punch indicates whether the file is a punch file (default=False)
        """
        ## the active filename (string)
        self.bdf_filename = bdf_filename

        if include_dir is None:
            include_dir = os.path.dirname(bdf_filename)
        ## the directory of the 1st BDF (include BDFs are relative to this one)
        self.include_dir = include_dir

        self.open_file(self.bdf_filename)

        self.log.debug('---starting BDF.readBDF of %s---' % self.bdf_filename)
        #self.log.info('xref=%s' % xref)

        if not punch:
            self._read_executive_control_deck()
            self._read_case_control_deck()

        self._read_bulk_data_deck()
        self.cross_reference(xref=xref)

        self.log.debug('---finished BDF.readBDF of %s---' % self.bdf_filename)

    def _read_executive_control_deck(self):
        """Reads the executive control deck"""
        lineUpper = ''
        while 'CEND' not in lineUpper[:4] and 'BEGIN' not in lineUpper and 'BULK' not in lineUpper:
            line, comment = self.get_next_line(False)
            #print("exec line",line)
            line = line.rstrip('\n\r\t ')
            if len(line) > 0:
                self.executive_control_lines.append(line)
            lineUpper = line.upper()
        self._parse_executive_control_deck()

    def _parse_executive_control_deck(self):
        """
        Extracts the solution from the executive control deck
        """
        for (i, eline) in enumerate(self.executive_control_lines):
            uline = eline.strip().upper()  # uppercase line
            uline = uline.split('$')[0]
            if 'SOL ' in uline[:4]:
                if ',' in uline:
                    sline = uline.split(',')  # SOL 600,method
                    solValue = sline[0]
                    method = sline[1]
                else:
                    solValue = uline
                    method = None

                sol = solValue[3:].strip()
                if self.sol is not None:
                    raise ValueError('cannot overwrite solution existing='
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
        @warning will be removed after v0.7 in favor of update_solution
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
        except ValueError:
            self.sol = self._solmap_to_value[sol]

        if self.sol == 600:
            ## solution 600 method modifier
            self.solMethod = method.strip()
            self.log.debug("sol=%s method=%s" % (self.sol, self.solMethod))
        else:  # very common
            self.solMethod = None

    def setDynamicSyntax(self, dictOfVars):
        """
        @see set_dynamic_syntax
        @warning will be removed after v0.7 in favor of set_dynamic_syntax
        """
        warnings.warn('setDynamicSyntax has been deprecated; use '
                      'set_dynamic_syntax', DeprecationWarning, stacklevel=2)
        self.set_dynamic_syntax(dictOfVars)

    def set_dynamic_syntax(self, dict_of_vars):
        """
        Uses the OpenMDAO syntax of %varName in an embedded BDF to
        update the values for an optimization study.
        Variables should be 7 characters to fit in an 8-character field.
        %varName
        dict_of_vars = {'varName': 10}
        @param self the BDF object
        @param dict_of_vars dictionary of 7 character variable names to map.
        """
        self.dict_of_vars = {}
        for (key, value) in sorted(dict_of_vars.iteritems()):
            assert len(key) <= 7, ('max length for key is 7; '
                                   'len(%s)=%s' % (key, len(key)))
            self.dict_of_vars[key.upper()] = value
        self._is_dynamic_syntax = True

    def _is_case_control_deck(self, line):
        """
        @todo not done...
        """
        lineUpper = line.upper().strip()
        if 'CEND' in line.upper():
            raise SyntaxError('invalid Case Control Deck card...CEND...')
        if '=' in lineUpper or ' ' in lineUpper:
            return True
        for card in self.uniqueBulkDataCards:
            lenCard = len(card)
            if card in lineUpper[:lenCard]:
                return False
        return True

    def _read_case_control_deck(self):
        """
        reads the case control deck
        @note called with recursion if an INCLUDE file is found
        """
        self.log.info("reading Case Control Deck...")
        line = ''
        while len(self.active_filenames) > 0:  # keep going until finished
            lineIn, comment = self.get_next_line(False)
            if lineIn is None:
                return  # file was closed
            line = lineIn.strip().split('$')[0].strip()
            lineUpper = line.upper()
            #print("lineUpper = |%s|" % lineUpper)
            if lineUpper.startswith('INCLUDE'):
                #print("INCLUDE!!!")
                next_line, comment = self.get_next_line(False)
                if next_line:
                    next_line = next_line.strip().split('$')[0].strip()
                else:
                    next_line = ''
                include_lines = [line]
                while '\\' in next_line or '/' in next_line:  # more includes
                    include_lines.append(next_line)
                    next_line, comment = self.get_next_line(False)
                    next_line = next_line.strip().split('$')[0].strip()

                filename = self._get_include_file_name(include_lines)

                self.open_file(filename)
                #line = next_line
                self.case_control_lines.append(next_line)
            else:
                self.case_control_lines.append(lineUpper)

            if 'BEGIN' in lineUpper and 'BULK' in lineUpper:
                self.log.debug('found the end of the case control deck!')
                break
        #self.log.info("finished with Case Control Deck..")

        #for line in self.case_control_lines:
            #print "** line=|%r|" %(line)

        self.caseControlDeck = CaseControlDeck(self.case_control_lines, self.log)
        self.caseControlDeck.solmap_toValue = self._solmap_to_value
        self.caseControlDeck.rsolmap_toStr = self.rsolmap_toStr

        return self.case_control_lines

    def _get_include_file_name(self, cardLines):
        """Parses an INCLUDE file split into multiple lines (as a list).
        @param self the BDF object
        @param cardLines the list of lines in the include card (all the lines!)
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
        filename = os.path.join(self.include_dir, filename)
        return filename

    def is_reject(self, card_name):
        """
        Can the card be read
        @param self the BDF object
        @param card_name the card_name -> 'GRID'
        """
        if card_name.startswith('='):
            return False
        elif card_name in self.cardsToRead:
            return False
        if card_name:
            if card_name not in self.reject_count:
                self.log.info("RejectCardName = |%s|" % card_name)
                self.reject_count[card_name] = 0
            self.reject_count[card_name] += 1
        return True

    def open_file(self, fname):
        """
        Opens the primary bdf/dat file and all subsequent INCLUDE files.
        @param fname:  the name of the bdf/dat file to open
        @return: None
        @note Doesn't allow reuse of the same bdf/dat file twice.
        """
        if fname in self.used_filenames:
            raise RuntimeError('fname=%s has already been opened once' % fname)
        #print("fname!!",fname)
        self.log.info("opening %s" % fname)
        self.used_filenames.append(fname)
        self.active_filenames.append(fname)
        self.active_filename = fname
        f = open(fname, 'rb')
        self.active_files.append(f)
        self.active_file = f
        self.stored_lines[self.active_filename] = []
        self.stored_comments[self.active_filename] = []

    def _close_file(self):
        """
        Closes the primary bdf/dat file and all subsequent INCLUDE files.
        @return: None
        """
        self.log.info("closing %s" % self.active_filename)
        if len(self.active_files) > 1:
            self.active_file.close()

            # pop the closing file
            self.active_files.pop()
            self.active_filenames.pop()

            # set the next file
            self.active_file = self.active_files[-1]
            self.active_filename = self.active_filenames[-1]
        else: # all files are closed, any lines left have been loaded
            self.active_file = None
            self.active_filename = None

    def get_next_line(self, access_stored_lines):
        """
        Gets the next Executive or Case Control Deck line.
        @param access_stored_lines: No purpose
        @return: line, comment
        """
        line = ''
        comment = ''
        try:
            while not line and not self._is_end_of_file():
                line = self.active_file.readline().rstrip('\n\r\t ')
                #print("*line = |%r|" % line)
                if '$' in line:
                    i = line.index('$')
                    comment = line[i:] + '\n'
                    line = line[:i].rstrip('\t ')
                if self._is_end_of_file():
                    self._close_file()
                    if line:
                        return line, comment
                    elif self.active_file is None:
                        return None, None
        except AttributeError:
            line = None
            comment = None
        # except ValueError:
        #     line = None
        #     comment = None
        return line, comment

    def _is_end_of_file(self):
        """
        Has the end of the current file been reached?
        @param self the BDF object
        @return:  True/False
        """
        try:
            fsize = os.fstat(self.active_file.fileno()).st_size
            if self.active_file.tell() == fsize:
                return True
        except:
            print("self.active_file = ", self.active_file)
            #print("self.active_file.tell() = ", self.active_file.tell())
            print(self.active_file.readline())
            #raise IOError('self.active_filename = |%s|' % self.active_filename)
            raise
        return False

    def _clean_comment(self, comment):
        """
        Removes specific pyNastran comment lines so duplicate lines aren't
        created.
        @param self the BDF object
        @param comment the comment to possibly remove
        """
        if comment[:-1] in ['$NODES', '$SPOINTS', '$ELEMENTS',
                       '$PARAMS', '$PROPERTIES', '$ELEMENTS_WITH_PROPERTIES',
                       '$ELEMENTS_WITH_NO_PROPERTIES (PID=0 and unanalyzed properties)',
                       '$UNASSOCIATED_PROPERTIES',
                       '$MATERIALS', '$THERMAL MATERIALS',
                       '$CONSTRAINTS', '$SPCs', '$MPCs', '$RIGID ELEMENTS',
                       '$LOADS', '$AERO', '$AERO CONTROL SURFACES',
                       '$FLUTTER', '$DYNAMIC', '$OPTIMIZATION',
                       '$COORDS', '$THERMAL', '$TABLES', '$RANDOM TABLES',
                       '$SETS', '$REJECTS', '$REJECT_LINES']:
            comment = ''
        return comment

    def get_next_bulk_line(self, access_stored_lines):
        """
        @param self the BDF object
        @param access_stored_lines Should all the stored lines be returned?
          This happens on the first line of a new card.  This keeps all the
          comments for it together.
        """
        comment = ''
        stored_lines = self.stored_lines[self.active_filename]
        if access_stored_lines and self.active_filename and stored_lines:
            if self.stored_comments[self.active_filename]:
                comment = ''.join(self.stored_comments[self.active_filename])
                self.stored_comments[self.active_filename] = []
            else:
                comment = ''
            line = stored_lines.pop(0)
            return line, comment

        try:
            line = self.active_file.readline().rstrip('\n\r\t ')
            if '$' in line:
                i = line.index('$')
                comment = self._clean_comment(line[i:] + '\n')
                line = line[:i].rstrip('\t ')
        except AttributeError:
            line = None
            comment = None
        except ValueError:
            line = None
            comment = None
        return line, comment

    def get_card(self):
        """
        Returns the next Bulk Data Card in the BDF
        @param self the BDF object
        @retval lines the lines of the card
        @retval comment the comment for the card
        @retval cardname the name of the card
        """
        lines = []
        comments = []

        line, comment = self.get_next_bulk_line(True)
        comments.append(comment)
        if line is None:
            return None, None, None

        while line is None or len(line) == 0 or not line[0:1].isalpha():
            line, comment = self.get_next_bulk_line(False)
            comments.append(comment)
            if self._is_end_of_file():
                lines.append(line)
                cardname = self.get_card_name(lines)
                self._close_file()
                return lines, comments, cardname
        lines.append(line)

        while 1:
            line2, comment = self.get_next_bulk_line(False)
            if line2 is None:
                break
            elif len(line2) > 0 and line2[0].isalpha():  # new card
                if self.active_filename is not None:
                    for line in self.stored_lines[self.active_filename]:
                        lines.append(line)
                    self.stored_lines[self.active_filename] = [line2]
                    if comment:
                        self.stored_comments[self.active_filename].append(comment)
                break
            elif len(line2) == 0: # empty line
                self.stored_lines[self.active_filename].append(line2)
                if comment:
                    self.stored_comments[self.active_filename].append(comment)
            else:
                if self.stored_lines[self.active_filename]:
                    lines += self.stored_lines[self.active_filename]
                    comments += self.stored_comments[self.active_filename]
                    self.stored_lines[self.active_filename] = []
                    self.stored_comments[self.active_filename] = []
                lines.append(line2)
                comments += comment

            if self._is_end_of_file():
                lines.append(line)
                cardname = self.get_card_name(lines)
                self._close_file()
                return(lines, comments, cardname)

        cardname = self.get_card_name(lines)
        return lines, comments, cardname

    def get_card_name(self, lines):
        """
        Returns the name of the card defined by the provided lines
        @param self the BDF object
        @param lines the lines of the card
        @retval cardname the name of the card
        """
        card_name = lines[0][:8].rstrip('\t, ').split(',')[0].split('\t')[0].strip('\t ')
        if len(card_name) == 0:
            return None
        assert ' ' not in card_name and len(card_name) > 0, 'card_name=|%s|\nline=|%s| in filename=%s is invalid' % (card_name, lines[0], self.active_filename)
        return card_name.upper()
    
    def _read_bulk_data_deck(self):
        """
        Parses the Bulk Data Deck
        @param self the BDF object
        """
        self.log.debug("*read_bulk_data_deck")
        while len(self.active_filenames) > 0 or self.stored_lines[self.active_filename]:
            (raw_card, comment, card_name) = self.get_card()
            #print("raw_card",raw_card)
            if card_name == 'INCLUDE':
                filename = self._get_include_file_name(raw_card)
                self.open_file(filename)
                reject = '$ INCLUDE processed:  %s\n' % (filename)
                self.rejects.append([reject])
                continue
            elif card_name is None:
                break

            self._increase_card_count(card_name)
            if not self.is_reject(card_name):
                self.add_card(raw_card, card_name, comment)
            else:
                self.rejects.append(raw_card)
                continue

            if 'ENDDATA' in card_name:
                break  # exits while loop
        ### end of while loop

    def _increase_card_count(self, card_name):
        """
        Used for testing to check that the number of cards going in is the
        same as each time the model is read verifies proper writing of cards
        @param self the BDF object
        @param card_name the card_name -> 'GRID'
        @warning this wont guarantee proper reading of cards, but will help
        """
        if card_name == '':  # stupid null case
            return

        if card_name in self.cardCount:
            self.cardCount[card_name] += 1
        else:
            self.cardCount[card_name] = 1

    def add_card(self, card_lines, card_name, comment=''):
        """
        adds a card object to the BDF object.
        @param self the BDF object
        @param card_lines the list of the card fields -> ['GRID',1,2,]
        @param card_name the card_name -> 'GRID'
        @param comment an optional the comment for the card
        @retval card_object the card object representation of card
        @note
          this is a very useful method for interfacing with the code
        @note
           the cardObject is not a card-type object...so not a GRID card
           or CQUAD4 object.  It's a BDFCard Object.  However, you know the
           type (assuming a GRID), so just call the mesh.Node(nid) to get the
           Node object that was just created.
        @warning
          cardObject is not returned
        """
        comment = ''.join(comment)
        if card_name in ['DEQATN']:
            card_obj = card_lines
            card = card_lines
        else:
            fields = to_fields(card_lines, card_name)

            # apply OPENMDAO syntax
            if self._is_dynamic_syntax:
                fields = [self.parse_dynamic_syntax(field) if '%' in field[0:1]
                          else field for field in fields]

            card = wipe_empty_fields([interpretValue(field, fields)
                                      for field in fields])
            card_obj = BDFCard(card)
        
        # function that gets by name the initialized object (from global scope)
        try:
            _get_cls = lambda name: globals()[name](card_obj, comment=comment)
        except Exception as e:
            if not e.args: 
                e.args=('',)
            e.args = (e.args[0] + "\ncard = %s" % card,)+e.args[1:]
        _cls = lambda name: globals()[name]

        if self._auto_reject:
            self.reject_cards.append(card)
            print('rejecting processed %s' % (card))
            return card_obj
        try:
            # cards that have their own method add_CARDNAME to add them
            if card_name in ['LSEQ', 'PHBDY', 'AERO', 'AEROS', 'AEFACT',
              'AELINK', 'AELIST', 'AEPARM', 'AESTAT', 'AESURF', 'TRIM',
              'FLUTTER', 'FLFACT', 'GUST', 'NLPARM', 'TSTEP', 'TSTEPNL',
              'SESET', 'DCONSTR', 'DESVAR', 'DDVAL', 'DLINK', 'PARAM',
              'PDAMPT', 'PELAST', 'PBUSHT']:
                try:
                    getattr(self, 'add_' + card_name)(_get_cls(card_name))
                except Exception as e:
                    if not e.args: 
                        e.args=('',)
                    e.args = (e.args[0] + "\ncard = %s" % card,)+e.args[1:]
                return card_obj

            # dictionary of cards. Key is the name of the function to add the
            # card
            # 'PCOMPG':  # hasnt been verified
            # 'MAT8':  # note there is no MAT6 or MAT7
            _cards = {
             'add_node': ['GRID'],
             'add_element': ['CQUAD4', 'CQUAD8', 'CQUAD', 'CQUADR', 'CQUADX',
                             'CTRIA3', 'CTRIA6', 'CTRIAR', 'CTRIAX', 'CTRIAX6',
                             'CBAR', 'CBEAM', 'CBEAM3', 'CROD', 'CONROD',
                             'CTUBE', 'CBEND', 'CELAS1', 'CELAS2', 'CELAS3',
                             'CELAS4', 'CONM1',  'CONM2', 'CMASS1', 'CMASS2',
                             'CMASS3', 'CMASS4', 'CVISC', 'CSHEAR', 'CGAP',
                             'CRAC2D', 'CRAC3D'],
             'add_damper_element': ['CBUSH', 'CBUSH1D', 'CFAST', 'CDAMP1',
                                    'CDAMP2', 'CDAMP3', 'CDAMP4', 'CDAMP5'],
             'add_rigid_element': ['RBAR', 'RBAR1', 'RBE1', 'RBE2', 'RBE3'],
             'add_property': ['PSHELL', 'PCOMP', 'PCOMPG', 'PSHEAR', 'PSOLID',
                              'PBAR', 'PBARL', 'PBEAM', 'PBCOMP', 'PBEAML',
                              'PROD', 'PTUBE', 'PLSOLID', 'PBUSH1D', 'PBUSH',
                              'PFAST', 'PDAMP5', 'PGAP', 'PRAC2D', 'PRAC3D'],

             # hasnt been verified, links up to MAT1, MAT2, MAT9 w/ same MID
             'add_creep_material': ['CREEP'],
             'add_material': ['MAT1', 'MAT2', 'MAT3', 'MAT8', 'MAT9', 'MAT10',
                              'MATHP'],
             'add_thermal_material': ['MAT4', 'MAT5'],
             'add_material_dependence': ['MATS1'],
             'add_load': ['FORCE', 'FORCE1', 'FORCE2', 'MOMENT', 'MOMENT1',
                          'MOMENT2', 'GRAV', 'ACCEL1', 'LOAD', 'PLOAD',
                              'PLOAD1', 'PLOAD2', 'PLOAD4', 'PLOADX1',
                              'RFORCE', 'DLOAD', 'SLOAD', 'TLOAD1', 'TLOAD2',
                              'RLOAD1', 'RLOAD2', 'RANDPS'],
             'add_thermal_load': ['TEMP', 'QBDY1', 'QBDY2', 'QBDY3', 'QHBDY'],
             'add_thermal_element': ['CHBDYE', 'CHBDYG', 'CHBDYP'],
             'add_convection_property': ['PCONV', 'PCONVM'],
             'add_thermal_BC': ['CONV', 'RADBC'],
             'add_constraint_MPC': ['MPC', 'MPCADD'],
             'add_constraint_SPC': ['SPC', 'SPC1', 'SPCAX', 'SPCD', 'SPCADD'],
             'add_suport': ['SUPORT'],  # pseudo-constraint
             'add_constraint': ['SUPORT1'],  # pseudo-constraint
             'add_SPLINE': ['SPLINE1', 'SPLINE2', 'SPLINE4', 'SPLINE5'],
             'add_CAERO': ['CAERO1', 'CAERO2'],
             'add_PAERO': ['PAERO1', 'PAERO2'],
             'add_MKAERO': ['MKAERO1', 'MKAERO2'],
             'add_FREQ': ['FREQ', 'FREQ1', 'FREQ2'],
             'add_ASET': ['ASET', 'ASET1'], 'add_BSET': ['BSET', 'BSET1'],
             'add_CSET': ['CSET', 'CSET1'], 'add_QSET': ['QSET', 'QSET1'],
             'add_SET': ['SET1', 'SET3'],
             'add_DRESP': ['DRESP1', 'DRESP2'],
             'add_DVPREL': ['DVPREL1', 'DVPREL2'],
             'add_coord': ['CORD2R', 'CORD2C', 'CORD2S'],
             'add_table': ['TABLED1', 'TABLED2', 'TABLED3', 'TABLEM1',
                           'TABLEM2', 'TABLEM3', 'TABLEM4', 'TABLES1',
                           'TABLEST'],
             'add_random_table': ['TABRND1', 'TABRNDG'],
             'add_method': ['EIGB', 'EIGR', 'EIGRL'],
             'add_cmethod': ['EIGC', 'EIGP'],
             'add_DVMREL': ['DVMREL1'],
            }

            for func, names in _cards.iteritems():
                if card_name in names:
                    try:
                        getattr(self, func)(_get_cls(card_name))
                    except Exception as e:
                        if not e.args: 
                            e.args=('',)
                        e.args = (e.args[0] + "\ncard = %s" % card,)+e.args[1:]
                        raise
                    return card_obj

            # card that requires more careful processing, elements
            _dct = {'CTETRA': (7, CTETRA4, CTETRA10), 'CHEXA': (11, CHEXA8,
                    CHEXA20), 'CPENTA': (9, CPENTA6, CPENTA15)}
            if card_name in _dct:
                d = _dct[card_name]
                self.add_element((d[1] if card_obj.nFields() == d[0]
                                       else d[2])(card_obj, comment=comment))
                return card_obj
            
            # dampers
            _dct = {'PELAS': (5,), 'PVISC': (5,), 'PDAMP': (3, 5)}
            if card_name in _dct:
                try:
                    self.add_property(_get_cls(card_name))
                except Exception as e:
                    if not e.args: 
                        e.args=('',)
                    e.args = (e.args[0] + "\ncard = %s" % card,)+e.args[1:]
                    raise
                for i in _dct[card_name]:
                    if card_obj.field(i):
                        self.add_property(_cls(card_name)(card_obj, 1,
                                          comment=comment))
                return card_obj

            if card_name in ['DEQATN']:  # buggy for commas
                #print 'DEQATN:  card_obj.card=%s' %(card_obj.card)
                #self.add_DEQATN(DEQATN(card_obj)) # should be later moved to
                self.rejects.append(card)          # for loop below
            elif card_name == 'GRDSET':
                self.gridSet = GRDSET(card_obj, comment=comment)
            elif card_name == 'DOPTPRM':
                self.doptprm = DOPTPRM(card_obj, comment=comment)

            elif card_name == 'DMIG':  # not done...
                if card_obj.field(2) == 'UACCEL':  # special DMIG card
                    self.reject_cards.append(card)
                elif card_obj.field(2) == 0:
                    self.add_DMIG(DMIG(card_obj, comment=comment))
                else:
                    self.dmigs[card_obj.field(1)].addColumn(card_obj,
                                                            comment=comment)

            elif card_name in ['DMI', 'DMIJ', 'DMIJI', 'DMIK']:
                if card_obj.field(2) == 0:
                    getattr(self, 'add_' + card_name)(_get_cls(card_name),
                                                      comment=comment)
                else:
                    getattr(self, card_name.lower() +
                            's')[card_obj.field(1)].addColumn(card_obj)
            # dynamic
            elif card_name == 'DAREA':
                self.add_DAREA(DAREA(card_obj, comment=comment))
                if card_obj.field(5):
                    self.add_DAREA(DAREA(card_obj, 1, comment=comment))

            elif card_name in ['CORD1R', 'CORD1C', 'CORD1S']:
                self.add_coord(_get_cls(card_name))
                if card_obj.field(5):
                    self.add_coord(_cls(card_name)(card_obj, nCoord=1,
                                                   comment=comment))

            elif card_name == 'PMASS':
                self.add_property(PMASS(card_obj, nOffset=0, comment=comment))
                for (i, j) in enumerate([3, 5, 7]):
                    if card_obj.field(j) is not None:
                        self.add_property(PMASS(card_obj, nOffset=i+1,
                                                comment=comment))

            elif card_name == 'SPOINT':
                self.add_SPOINT(SPOINTs(card_obj, comment=comment))
            elif 'ENDDATA' in card_name:
                self.foundEndData = True
            else:
                ## @warning cards with = signs in them
                ## are not announced when they are rejected
                if '=' not in card[0]:
                    self.log.info('rejecting processed %s' % card)
                self.reject_cards.append(card)
        except Exception as e:
            print(str(e))
            self.log.debug("card_name = |%r|" % (card_name))
            self.log.debug("failed! Unreduced Card=%s\n" % list_print(card) )
            self.log.debug("filename = %s\n" % self.bdf_filename)
            raise

        return card_obj

    def card_stats(self):
        """
        Print statistics for the BDF
        @param self the BDF object
        @note
          if a card is not supported and not added to the proper lists,
          this method will fail
        """
        card_stats = [
            'params', 'nodes', 'points', 'elements', 'rigidElements',
            'properties', 'materials', 'materialDeps', 'creepMaterials',
            'coords', 'mpcs', 'mpcadds',

            # dynamic cards
            'dareas', 'nlparms', 'tsteps', 'tstepnls',

            # direct matrix input - DMIG - dict
            'dmis', 'dmigs', 'dmijs', 'dmijis', 'dmiks',
            'dequations',

            # frequencies - dict
            'frequencies',

            # optimization - dict
            'dconstrs', 'desvars', 'ddvals', 'dlinks', 'dresps',
            'dvprels', 'dvmrels',

            # SESETx - dict
            'setsSuper',

            # tables
            'tables', 'randomTables',

            # methods
            'methods', 'cMethods',

            # aero
            'caeros', 'paeros', 'aero', 'aeros', 'aefacts', 'aelinks',
            'aelists', 'aeparams', 'aesurfs', 'aestats', 'gusts', 'flfacts',
            'flutters', 'splines', 'trims',

            # thermal
            'bcs', 'thermalMaterials', 'phbdys',
            'convectionProperties', ]

        ignored_types = set([
            'spoints', 'spointi',  # singleton
            'gridSet',  # singleton

            'spcs', 'spcadds',

            'suports',  # suport, suport1 - list
            'doptprm',  # singleton

            # SETx - list
            'sets', 'asets', 'bsets', 'csets', 'qsets',
        ])

        ignored_types2 = set([
            'caseControlDeck', 'spcObject2', 'mpcObject2',

            # done
            'sol', 'loads', 'mkaeros',
            'rejects', 'reject_cards',

            # not cards
            'debug',  'executive_control_lines',
            'case_control_lines', 'cardsToRead', 'cardCount',
            'isStructured', 'uniqueBulkDataCards',
            'nCardLinesMax', 'modelType', 'includeDir',
            'cardsToWrite', 'solMethod', 'log', 'doneReading',
            'linesPack', 'lineNumbers', 'iSolLine',
            'reject_count', 'relpath', 'isOpened',
            'foundEndData', 'specialCards',
            'infilesPack'])

        all_params = object_attributes(self)
        # removing variables that are not supported
        for attribute_name in ignored_types.union(ignored_types2):
            try:
                all_params.remove(attribute_name)
                #print('removing attribute_name=%s' % (attribute_name))
            except ValueError:
                pass

        msg = ['---BDF Statistics---']
        # sol
        msg.append('SOL %s\n' % self.sol)

        # loads
        for (lid, loads) in sorted(self.loads.iteritems()):
            msg.append('bdf.loads[%s]' % lid)
            groups = {}
            for load in loads:
                groups[load.type] = groups.get(load.type, 0) + 1
            for name, n in sorted(groups.iteritems()):
                msg.append('  %-8s %s' % (name + ':', n))
            msg.append('')

        #mkaeros
        if self.mkaeros:
            msg.append('bdf:mkaeros')
            msg.append('  %-8s %s' % ('MKAERO:', len(self.mkaeros)))

        for card_group_name in card_stats:
            card_group = getattr(self, card_group_name)
            groups = set([])
            #print("card_group_name = ", card_group_name)
            #if isinstance(card_group, list):
                #print("card_group_name = ", card_group_name)

            for card in card_group.itervalues():
                if isinstance(card, list):
                    for card2 in card:
                        groups.add(card2.type)
                else:
                    groups.add(card.type)

            group_msg = []
            for card_name in sorted(groups):
                try:
                    ncards = self.cardCount[card_name]
                    group_msg.append('  %-8s %s' % (card_name + ':', ncards))
                except KeyError:
                    assert card_name == 'CORD2R'
            if group_msg:
                msg.append('bdf.%s' % (card_group_name))
                msg.append('\n'.join(group_msg))
                msg.append('')

        # rejects
        if self.rejects:
            msg.append('Rejected Cards')
            for name, counter in sorted(self.cardCount.iteritems()):
                if name not in self.cardsToRead:
                    msg.append('  %-8s %s' % (name + ':', counter))
        msg.append('')
        return '\n'.join(msg)


def to_fields(card_lines, card_name):
    """
    Converts a series of lines in a card into string versions of the field.
    Handles large, small, and CSV formatted cards.
    @param lines the lines of the BDF card object
    @param card_name the card_name -> 'GRID'
    @retval fields the string formatted fields of the card
    """
    fields = []
    
    # first line
    line = card_lines.pop(0)
    if '=' in line:
            raise SyntaxError('card_name=%s\nequal signs are not supported...'
                              'line=|%r|' % (card_name, line))
        
    if '\t' in line:
        line = line.expandtabs()
        if ',' in line:
            raise SyntaxError('tabs and commas in the same line are '
                              'not supported...line=|%r|' % line)

    if '*' in line:  # large field
        if ',' in line:  # csv
            new_fields = line[:72].split(',')[:5]
            fields += new_fields
            for i in range(5-len(new_fields)):
                fields.append('')
        else:  # standard
            fields = [line[0:8], line[8:24], line[24:40], line[40:56],
                      line[56:72]]
    else: # small field
        if ',' in line:  # csv
            new_fields = line[:72].split(',')[:8]
            fields += new_fields
            for i in range(9-len(new_fields)):
                fields.append('')
        else:  # standard
            fields += [line[0:8], line[8:16], line[16:24], line[24:32],
                       line[32:40], line[40:48], line[48:56], line[56:64],
                       line[64:72]]
    
    for j, line in enumerate(card_lines): # ccntinuation lines
        #for i, field in enumerate(fields):
        #    if field.strip() == '+':
        #        raise RuntimeError('j=%s field[%s] is a +' % (j,i))

        if '=' in line and card_name != 'EIGRL':
            raise SyntaxError('card_name=%s\nequal signs are not supported...'
                              'line=|%r|' % (card_name, line))
        if '\t' in line:
            line = line.expandtabs()
            if ',' in line:
                raise SyntaxError('tabs and commas in the same line are '
                                  'not supported...line=|%r|' % line)
        
        if '*' in line:  # large field
            if ',' in line:  # csv
                new_fields = line[:72].split(',')[1:5]
                assert len(new_fields) == 4
                fields += new_fields
                for i in range(4-len(new_fields)):
                    fields.append('')
            else:  # standard
                fields += [line[8:24], line[24:40], line[40:56], line[56:72]]
        else: # small field
            if ',' in line:  # csv
                new_fields = line[:72].split(',')[1:8]
                fields += new_fields
                for i in range(8-len(new_fields)):
                    fields.append('')
            else:  # standard
                fields += [line[8:16], line[16:24], line[24:32], line[32:40],
                           line[40:48], line[48:56], line[56:64], line[64:72]]
    return fields


if __name__ == '__main__':
    bdf = BDF()
    infilename = 'f.bdf'
    bdf.readBDF(infilename, include_dir=None, xref=True, punch=False)
    print(bdf)
