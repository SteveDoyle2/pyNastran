# pylint: disable=C0103,C0302,R0901,R0902,R0904,R0912,R0914,R0915,W0201,W0611
"""
Main BDF class
"""
from __future__ import (nested_scopes, generators, division, absolute_import,
                        unicode_literals)

#from __future__ import (nested_scopes, generators, division, absolute_import)
                        #print_function, unicode_literals)

#import io
import os
import sys
import warnings
#from pyNastran.bdf.bdf2 import (to_fields, wipe_empty_fields, interpret_value,
#                                get_include_filename, parse_executive_control_deck)
from pyNastran.bdf.bdf2 import *


from pyNastran.utils import list_print
from pyNastran.utils import object_attributes
from pyNastran.utils.log import get_logger

if 0:
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
                                    constraintObject)
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

    from .caseControlDeck import CaseControlDeck
    from .bdf_Methods import BDFMethods
    from .bdfInterface.getCard import GetMethods
    from .bdfInterface.addCard import AddMethods
    from .bdfInterface.BDF_Card import BDFCard, wipe_empty_fields
    from .bdfInterface.bdf_Reader import print_filename # BDFReader
    from .bdfInterface.bdf_writeMesh import WriteMesh
    from .bdfInterface.bdf_cardMethods import interpret_value
    from .bdfInterface.crossReference import XrefMesh


class BDFDeprecated(object):
    def readBDF(self, bdf_filename, include_dir=None, xref=True, punch=False):
        """
        @see read_bdf
        @warning will be removed after v0.7 in favor of read_bdf
        """
        warnings.warn('readBDF has been deprecated; use '
                      'read_bdf', DeprecationWarning, stacklevel=2)
        self.read_bdf(bdf_filename, include_dir, xref, punch)

    def updateSolution(self, sol, method=None):
        """
        @see update_solution
        @warning will be removed after v0.7 in favor of update_solution
        """
        warnings.warn('updateSolution has been deprecated; use '
                      'update_solution', DeprecationWarning, stacklevel=2)
        self.update_solution(sol, method)

    def setDynamicSyntax(self, dictOfVars):
        """
        @see set_dynamic_syntax
        @warning will be removed after v0.7 in favor of set_dynamic_syntax
        """
        warnings.warn('setDynamicSyntax has been deprecated; use '
                      'set_dynamic_syntax', DeprecationWarning, stacklevel=2)
        self.set_dynamic_syntax(dictOfVars)

    def addCard(self, card, cardName, iCard=0, oldCardObj=None):
        """
        @see add_card
        @warning will be removed after v0.7 in favor of add_card
        """
        warnings.warn('addCard has been deprecated; use add_card',
                      DeprecationWarning, stacklevel=2)
        return self.add_card(card, cardName, icard=iCard,
                             old_card_obj=oldCardObj)

    def disableCards(self, cards):
        """
        @see disable_cards
        @warning will be removed after v0.7 in favor of disable_cards
        """
        warnings.warn('disableCards has been deprecated; use '
                      'disable_cards', DeprecationWarning, stacklevel=2)
        self.disable_cards(cards)
   
class BDF(BDFMethods, GetMethods, AddMethods, WriteMesh, XrefMesh, BDFDeprecated):
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
        # file management parameters
        self.include_dir = ''
        self.active_filename = None
        self.active_filenames = []
        self.used_filenames = []
        self._line_generators = []

        self.stored_Is = []
        self.stored_lines = []
        self.stored_comments = []
        
        ## list of all read in cards - useful in determining if entire BDF
        ## was read & really useful in debugging
        self.card_count = {}
        ## stores the card_count of cards that have been rejected
        self.reject_count = {}
        ## was an ENDDATA card found
        self.foundEndData = False

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

        self.case_control_lines = []
        self.executive_control_lines = []
        self.__init_attributes()

        ## the list of possible cards that will be parsed
        self.cardsToRead = set([
            'PARAM',
            'GRID', 'GRDSET', 'SPOINT',  # 'RINGAX',
            #'POINT',

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
            #'DMIG', 'DMIJ', 'DMIJI', 'DMIK', 'DMI',
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

    def disable_cards(self, cards):
        """
        Method for removing broken cards from the reader
        @param self the BDF object
        @param cards a list/set of cards that should not be read
        """
        disableSet = set(cards)
        self.cardsToRead.difference(disableSet)

    #def _is_special_card(self, cardName):
        #"""These cards are listed in the case control and the bulk data deck"""
        #if cardName in self.specialCards:
            #return True
        #return False

    def __init_attributes(self):
        """
        Creates storage objects for the BDF object.
        This would be in the init but doing it this way allows for better
        inheritance

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

        # ------------------------ structural defaults -----------------------
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
        self.spcObject2 = constraintObject()
        ## stores MPCADD,MPC
        self.mpcObject2 = constraintObject()

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

        # --------------------------- aero defaults --------------------------
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
        self.flfacts = {}  ## @todo can this be simplified ???
        ## stores FLUTTER
        self.flutters = {}
        ## mkaeros
        self.mkaeros = []
        ## store SPLINE1,SPLINE2,SPLINE4,SPLINE5
        self.splines = {}
        ## stores TRIM
        self.trims = {}

        # ------------------------- thermal defaults -------------------------
        # BCs
        ## stores thermal boundary conditions - CONV,RADBC
        self.bcs = {}  # e.g. RADBC
        ## defines the MAT4, MAT5, MATT4, etc.  @todo verify MATT4
        self.thermalMaterials = {}

        ## stores PHBDY
        self.phbdys = {}
        ## stores convection properties - PCONV, PCONVM ???
        self.convectionProperties = {}

    def read_bdf(self, bdf_filename, include_dir=None, xref=True, punch=False):
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

        self._open_file(self.bdf_filename)
        self.log.debug('---starting BDF.read_bdf of %s---' % self.bdf_filename)

        #self.get_line_gen = self.get_next_line()
        self.gen2 = self.get_card()
        #print "gen2 = ", self.gen2

        if not punch:
            self._read_executive_control_deck()
            self._read_case_control_deck()

        self._read_bulk_data_deck()
        self.cross_reference(xref=xref)

        self.log.debug('---finished BDF.read_bdf of %s---' % self.bdf_filename)
    
    def _read_executive_control_deck(self):
        """Reads the executive control deck"""
        lineUpper = ''
        while 'CEND' not in lineUpper[:4] and 'BEGIN' not in lineUpper and 'BULK' not in lineUpper:
            (i, line, comment) = self.get_line_gen.next()
            #print("exec line",line)
            line = line.rstrip('\n\r\t ')
            if len(line) > 0:
                self.executive_control_lines.append(line)
            lineUpper = line.upper()

        if 'CEND' in lineUpper[:4]:
            self.has_case_control_deck = True
        else:
            self.has_case_control_deck = False
            (i, line, comment) = self.get_line_gen.next()   # BEGIN BULK
            #print('execline2 = ',line)

        sol, method, iSolLine = parse_executive_control_deck(self.executive_control_lines)
        self.sol = sol
        self.iSolLine = iSolLine
        try:
            self.update_solution(sol, method)
        except:
            msg = ('update_solution failed...line=%s' % uline)
            raise RuntimeError(msg)

    def update_solution(self, sol, method=None):
        """
        Updates the overall solution type (e.g. 101,200,600)
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

    def set_dynamic_syntax(self, dict_of_vars):
        """
        Uses the OpenMDAO syntax of %varName in an embedded BDF to
        update the values for an optimization study.
        Variables should be 7 characters to fit in an 8-character field.
        Case sensitivity is supported.
        %varName
        dict_of_vars = {'varName': 10}
        @param self the BDF object
        @param dict_of_vars dictionary of 7 character variable names to map.
        """
        self.dict_of_vars = {}
        assert len(dict_of_vars) > 0, 'nvars = %s' % len(dict_of_vars)
        for (key, value) in sorted(dict_of_vars.iteritems()):
            assert len(key) <= 7, ('max length for key is 7; '
                                   'len(%s)=%s' % (key, len(key)))
            assert len(key) <= 1, ('max length for key is 1; '
                                   'len(%s)=%s' % (key, len(key)))
            assert isinstance(key, str), 'key=%s must be a string' % key
            self.dict_of_vars[key] = value
        self._is_dynamic_syntax = True

    def _parse_dynamic_syntax(self, key):
        """
        Applies the dynamic syntax for %varName
        @param self the object pointer
        @param key the uppercased key
        @retval value the dynamic value defined by dictOfVars
        @note
          %varName is actually %VARNAME b/c of auto-uppercasing the string,
          so the setDynamicSyntax method uppercases the key prior to this step.
        @see setDynamicSyntax
        """
        #print "*** valueRaw.lstrip() = |%r|" % valueRaw.lstrip()
        #key = key.lstrip('%%')
        key = key[1:].strip()
        self.log.info("dynamic key = |%r|" % key)
        #self.dict_of_vars = {'P5':0.5,'ONEK':1000.}
        if key not in self.dict_of_vars:
            msg = "key=|%r| not found in keys=%s" % (key,
                                                     self.dict_of_vars.keys())
            raise KeyError(msg)
        return self.dict_of_vars[key]

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
        Reads the case control deck
        @note called with recursion if an INCLUDE file is found
        """
        if not self.has_case_control_deck:
            return
        self.log.info("reading Case Control Deck...")
        line = ''
        while len(self.active_filenames) > 0:  # keep going until finished
            (i, lineIn, comment) = self.get_line_gen.next()
            if lineIn is None:
                return  # file was closed
            line = lineIn.strip().split('$')[0].strip()
            lineUpper = line.upper()
            #print("lineUpper = |%s|" % lineUpper)
            if lineUpper.startswith('INCLUDE'):
                #print("INCLUDE!!!")
                (i, next_line, comment) = self.get_line_gen.next()
                if next_line:
                    next_line = next_line.strip().split('$')[0].strip()
                else:
                    next_line = ''
                include_lines = [line]
                while '\\' in next_line or '/' in next_line:  # more includes
                    include_lines.append(next_line)
                    (i, line_next, comment) = self.get_line_gen.next()
                    next_line = next_line.strip().split('$')[0].strip()

                filename = get_include_filename(include_lines,
                                                 include_dir=self.include_dir)

                self._open_file(filename)
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

    def _open_file(self, bdf_filename):
        """
        Opens the primary bdf/dat file and all subsequent INCLUDE files.
        @param fname:  the name of the bdf/dat file to open
        @return: None
        @note Doesn't allow reuse of the same bdf/dat file twice.
        """
        #print "opening self.active_filename=%s" % bdf_filename
        if bdf_filename in self.used_filenames:
            msg = 'bdf_filename=%s has already been opened once.\nused_filenames=%s' % self.used_filenames
            raise RuntimeError(msg)
        self.log.info("opening %s" % bdf_filename)
        self.active_filename = bdf_filename
        self.used_filenames.append(bdf_filename)
        
        #if len(self.used_filenames) > 1:
        self.active_filenames.append(bdf_filename)
        #print "active_filenames =", self.active_filenames
        self.get_line_gen = self.get_next_line()
        self._line_generators.append(self.get_line_gen)

    def _close_file(self):
        """
        Closes the primary bdf/dat file and all subsequent INCLUDE files.
        @return: None
        """
        self.log.info("closing %s" % self.active_filename)
        #print ""
        #print "closing self.active_filename=%s" % self.active_filename
        #print "*active_filenames =", self.active_filenames
        if self.active_filenames:
            # pop the closing file
            self.active_filename = self.active_filenames.pop()
            #print "the new active file is self.active_filename=%s" % self.active_filename
            self.get_line_gen = self._line_generators.pop()
        else:
            self.active_filename = None
            self.get_line_gen = None
        
    def get_next_line(self):
        """
        Gets the next Executive or Case Control Deck line.
        builds a generator to read the BDF in a much cleaner way
        @return: line, comment
        """
        with open(self.active_filename) as file:
            for n,line in enumerate(file):
                #print "%%%%%%%%%%%%%%%%%%%%%%%%%%%%"
                line = line.rstrip('\t\r\n ')
                if '$' in line:
                    i = line.index('$')
                    comment = line[i:] + '\n'
                    line = line[:i].rstrip('\t ')
                    #print "  **n=%s line=%s comment=%s" % (n, line.strip(), comment)
                else:
                    #print "  **n=%s line=%s" % (n, line.strip())
                    pass
                
                comment = ''
                #print "  stored lines =", self.stored_lines
                #print "  yieldingB new line=|%s|" % line
                #print "############################"
                yield n, line, comment

                #print "  stored lines2 =", self.stored_lines
                while self.stored_lines:
                    comment = ''
                    #print "  trying to store lines...line=%s" % line
                    i2 = self.stored_Is.pop(0)
                    line2 = self.stored_lines.pop(0)
                    if self.stored_comments:
                        comment = ''.join(self.stored_comments)
                        self.stored_comments = []

                    #self.stored_Is.append(n)
                    #self.stored_lines.append(line)
                    #print "  yieldingA line=|%s|" % line2
                    #print "############################"
                    yield i2, line2, comment
                
        file.close()
        #while self.stored_lines:
            #yield -1, self.stored_lines.pop(0), ''
        ###

    def get_next_bulk_line(self, access_stored_lines):
        """
        @param self the BDF object
        @param access_stored_lines Should all the stored lines be returned?
          This happens on the first line of a new card.  This keeps all the
          comments for it together.

        @todo has a problem with comment lines that are embedded in a single line
        TABLED1,1000,,,,,,,,+TAB1 $ THE ONE PER CENT DAMPING LINE
        +TAB1   0.0     0.0     5.12    200.    100.    200.    200.    0.      +TAB2
        $       F0      A0      F1      A1      F2      A2      F3      A3
        +TAB2   1000.   0.      ENDT
        $       F4      A4

        but this should also be supported...
        PBEAM,
        +
        1.0
        where the + line represents a line (a + is not required, hence the difficulty)
        """
        comment = ''

    def get_card(self):
        """
        Returns the next Bulk Data Card in the BDF
        @param self the BDF object
        @retval lines the lines of the card
        @retval comment the comment for the card
        @retval cardname the name of the card
        """
        #print '--------------------------'
        for (i, line, comment) in self.get_line_gen:
            #print "new card...line =|%s|" % line
            Is = []
            lines = []
            comments = []
            #c = ''

            if comment:
                #c=' comment=|%s|' % comment.strip()
                comments.append(comment)
            #print "lineA1 %s line=|%s|%s" % (i, line.strip(), c)

            if line:
                #print "adding line"
                Is.append(i)
                lines.append(line)
            else:
                while len(line)==0:
                    #print "you cant have an empty first line..."
                    (i, line, comment) = self.get_line_gen.next()
                    Is.append(i)
                    lines.append(line)

                    c = ''
                    if comment:
                        #c=' comment=|%s|' % comment.strip()
                        comments.append(comment)
                    #print "lineA2 %s line=|%s|%s" % (i, line.strip(), c)

            #print "lines =",lines

            # get another line
            #print "*lineC %s line=|%s|%s" % (i, line.strip(), c)
            try:
                (i, line, comment) = self.get_line_gen.next()
            except StopIteration:
                yield lines, comments

            #if comment:  c=' comment=|%s|' % comment.strip()
            #print "lineC %s line=|%s|%s" % (i, line.strip(), c)
            #print "lines =",lines
            #print ""
            
            in_loop = False
            while len(line)==0 or line[0] in [' ', '*', '+', ',']:
                in_loop = True
                #print "into the loop!"
                Is.append(i)
                lines.append(line)
                #sys.exit('into the loop')
                if comment:
                    comments.append(comment)
                
                if comment:  c=' comment=|%s|' % comment.strip()
                #print "lineD %s line=|%s|%s" % (i, line.strip(), c)
                (i, line, comment) = self.get_line_gen.next()

            #if not in_loop:
                #self.stored_Is.append(i)
                #self.stored_lines.append(line)
                #self.stored_comments.append(comment)
                #print "non-continuation line = ", line
                #print "lines = ", lines

            #sys.exit('out of the the loop')
            if line[0] not in [' ', '*', '+', ',']:
            #if in_loop:
                #print "stored_lines =", self.stored_lines
                #print "storing line..."
                self.stored_Is.append(i)
                self.stored_lines.append(line)
                if comment:
                    self.stored_comments.append(comment)
            #if comment:
            #    self.stored_comments.append(comment)
            #print "linesE = %s" % lines
            
            #print "lines =",lines
            lines2 = []
            found_lines = False
            for line in lines:
                if line:
                    found_lines = True
                if found_lines:
                    lines2.append(line)
            #print "lines2 =",lines2
            yield lines2, comments
            #print '--------------------------'

        #return line, comment

    def get_card_name(self, lines):
        """
        Returns the name of the card defined by the provided lines
        @param self the BDF object
        @param lines the lines of the card
        @retval cardname the name of the card
        """
        card_name = lines[0][:8].rstrip('\t, ').split(',')[0].split('\t')[0].strip('*\t ')
        if len(card_name) == 0:
            return None
        assert ' ' not in card_name and len(card_name) > 0, 'card_name=|%r|\nline=|%s| in filename=%s is invalid' % (card_name, lines[0], self.active_filename)
        return card_name.upper()

    def _read_bulk_data_deck(self):
        """
        Parses the Bulk Data Deck
        @param self the BDF object
        """
        self.log.debug("*read_bulk_data_deck")
        n = 1
        isDone = False
        for (lines, comments) in self.gen2:
            #print "*asdf*lines =", lines
            if not lines:
                self._close_file()
                if self.active_filename:
                    continue
                break
            
            card_name = self.get_card_name(lines)
            #print "card_name =", card_name
            if card_name == 'INCLUDE':
                bdf_filename = get_include_filename(lines, include_dir=self.include_dir)
                #print "newfname =", newfname
                self._open_file(bdf_filename)
                reject = '$ INCLUDE processed:  %s\n' % bdf_filename
                self.rejects.append([reject])
                continue
            #print "*asdf*comments =", comments
            #print "*asdf*card_name =", card_name
            #fields = to_fields(lines, card_name)
            #card = wipe_empty_fields([interpret_value(field, fields)
            #                          for field in fields])

            self._increase_card_count(card_name)
            if not self.is_reject(card_name):
                self.add_card(lines, card_name)
            else:
                self.rejects.append(lines)
                continue

            if 'ENDDATA' in card_name:
                break  # exits while loop
            
            
            #print "bulk card =", card
            
            isDone = False

            c = ''
            if comments:  c=' comments=%s' % comments
            #print "ncard=%s lines=%s%s" % (n, lines, c)
            n += 1

            #print "**********************"
        ### end of all the files

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

        if card_name in self.card_count:
            self.card_count[card_name] += 1
        else:
            self.card_count[card_name] = 1

    def add_card(self, card_lines, card_name, comment=''):
        """
        Adds a card object to the BDF object.
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
            #print("_is_dynamic_syntax =", self._is_dynamic_syntax)
            if self._is_dynamic_syntax:
                fields = [self._parse_dynamic_syntax(field) if '%' in
                          field[0:1] else field for field in fields]

            card = wipe_empty_fields([interpret_value(field, fields)
                                      for field in fields])
            print "card =", card
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
                    ## PHBDY -> add_PHBDY
                    getattr(self, 'add_' + card_name)(_get_cls(card_name))
                except Exception as e:
                    if not e.args: 
                        e.args = ('',)
                    e.args = (e.args[0] + "\ncard = %s" % card,) + e.args[1:]
                    raise
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
                    getattr(self, 'add_' + card_name)(_get_cls(card_name))
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

            elif card_name == 'CONV':
                bc = CONV(card_obj, comment=comment)
                self.add_thermal_BC(bc, bc.eid)
            #elif card_name == 'RADM':
            #    bc = RADM(card_obj, comment=comment)
            #    self.add_thermal_BC(bc, bc.nodamb)
            elif card_name == 'RADBC':
                bc = RADBC(card_obj, comment=comment)
                self.add_thermal_BC(bc, bc.nodamb)

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
            'case_control_lines', 'cardsToRead', 'card_count',
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
                #print('removing attribute_name=%s' % attribute_name)
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
                    ncards = self.card_count[card_name]
                    group_msg.append('  %-8s %s' % (card_name + ':', ncards))
                except KeyError:
                    assert card_name == 'CORD2R'
            if group_msg:
                msg.append('bdf.%s' % card_group_name)
                msg.append('\n'.join(group_msg))
                msg.append('')

        # rejects
        if self.rejects:
            msg.append('Rejected Cards')
            for name, counter in sorted(self.card_count.iteritems()):
                if name not in self.cardsToRead:
                    msg.append('  %-8s %s' % (name + ':', counter))
        msg.append('')
        return '\n'.join(msg)

def _clean_comment(comment):
    """
    Removes specific pyNastran comment lines so duplicate lines aren't
    created.
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


def to_fields(card_lines, card_name):
    """
    Converts a series of lines in a card into string versions of the field.
    Handles large, small, and CSV formatted cards.
    @param lines the lines of the BDF card object
    @param card_name the card_name -> 'GRID'
    @retval fields the string formatted fields of the card
    """
    fields = []
    #print('---------------')
    # first line
    line = card_lines.pop(0)
    print "first line=|%s|" % line.rstrip()
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
            for i in range(5-len(new_fields)):
                new_fields.append('')
        else:  # standard
            new_fields = [line[0:8], line[8:24], line[24:40], line[40:56],
                          line[56:72]]
        fields += new_fields
        assert len(fields) == 5
    else:  # small field
        print "small line=|%s|" % line.rstrip()
        if ',' in line:  # csv
            new_fields = line[:72].split(',')[:9]
            for i in range(9-len(new_fields)):
                new_fields.append('')
        else:  # standard
            new_fields = [line[0:8], line[8:16], line[16:24], line[24:32],
                          line[32:40], line[40:48], line[48:56], line[56:64],
                          line[64:72]]
        fields += new_fields
        assert len(fields) == 9

    #print("new_fieldsA =",new_fields)

    for j, line in enumerate(card_lines): # continuation lines
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
                for i in range(4-len(new_fields)):
                    new_fields.append('')
            else:  # standard
                new_fields = [line[8:24], line[24:40], line[40:56], line[56:72]]
            assert len(new_fields) == 4
        else:  # small field
            print "small lin2=|%s|" % line.rstrip()
            if ',' in line:  # csv
                new_fields = line[:72].split(',')[1:9]
                for i in range(8-len(new_fields)):
                    new_fields.append('')
            else:  # standard
                new_fields = [line[8:16], line[16:24], line[24:32], line[32:40],
                              line[40:48], line[48:56], line[56:64], line[64:72]]
            assert len(new_fields) == 8, 'nFields=%s new_fields=%s' % (len(new_fields), new_fields)

        fields += new_fields
        #print("new_fieldsB =", new_fields)

    return fields

def get_include_filename(cardLines, include_dir=''):
    """
    Parses an INCLUDE file split into multiple lines (as a list).
    @param cardLines the list of lines in the include card (all the lines!)
    @param include_dir the include directory (default='')
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
    filename = os.path.join(include_dir, filename)
    return filename


def parse_executive_control_deck(executive_control_lines):
    """
    Extracts the solution from the executive control deck
    """
    sol = None
    for (i, eline) in enumerate(executive_control_lines):
        uline = eline.strip().upper()  # uppercase line
        uline = uline.split('$')[0]
        print(uline)
        if 'SOL ' in uline[:4]:
            if ',' in uline:
                sline = uline.split(',')  # SOL 600,method
                solValue = sline[0]
                method = sline[1]
            else:
                solValue = uline
                method = None

            #print("sol =", sol)
            if sol is None:
                sol = solValue[3:].strip()
            else:
                raise ValueError('cannot overwrite solution existing='
                                 '|SOL %s| new =|%s|' % (sol, uline))
            iSolLine = i
    #print("sol = ", sol)
    return sol, method, iSolLine


if __name__ == '__main__':
    bdf = BDF()
    import pyNastran
    pkg_path = pyNastran.__path__[0]
    bdfname = os.path.join(pkg_path, '..', 'models', 'solidBending.bdf')
    #print "bdfname =", bdfname
    bdf.read_bdf(bdfname)
    bdf.write_bdf('fem.out.bdf')

