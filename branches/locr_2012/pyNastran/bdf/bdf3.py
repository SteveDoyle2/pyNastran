#from __future__ import (nested_scopes, generators, division, absolute_import,
#                        print_function, unicode_literals)

#from __future__ import (nested_scopes, generators, division, absolute_import)
                        #print_function, unicode_literals)

import os
import sys
#from pyNastran.bdf.bdf2 import (to_fields, wipe_empty_fields, interpret_value,
#                                get_include_filename, parse_executive_control_deck)
from pyNastran.bdf.bdf2 import *

#import io
import os
import sys
import warnings

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


class BDF3(BDFMethods, GetMethods, AddMethods, WriteMesh, XrefMesh):
    """
    NASTRAN BDF Reader/Writer/Editor class.
    """
    modelType = 'nastran'

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

    def __init__(self, debug=True, log=None):
        # file management parameters
        self.include_dir = ''
        self.active_filename = None
        self.active_filenames = []
        self.all_filenames = []
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
        if include_dir is None:
            include_dir = os.path.dirname(bdf_filename)
        ## the directory of the 1st BDF (include BDFs are relative to this one)
        self.include_dir = include_dir

        #self.log.debug('---starting BDF.read_bdf of %s---' % self.bdf_filename)
        self.open_file(bdf_filename)

        #self.get_line_gen = self.stream_file()
        self.gen2 = self.get_card()
        #print "gen2 = ", self.gen2

        if not punch:
            self._read_executive_control_deck()
            self._read_case_control_deck()
        self.read_bulk_data_deck()
        #self.cross_reference(xref=xref)

        #self.log.debug('---finished BDF.read_bdf of %s---' % self.bdf_filename)
    
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

        #self._parse_executive_control_deck()
        #sol, method, iSolLine = parse_executive_control_deck(self.executive_control_lines)
        #self.sol = sol
        #self.iSolLine = iSolLine
        #try:
            #self.update_solution(sol, method)
        #except:
            #msg = ('update_solution failed...line=%s' % uline)
            #raise RuntimeError(msg)

    #def read_executive_deck(self):
        #self.executive_lines = []
        #(i, line, comment) = self.get_line_gen.next()
        #while 'CEND' not in line:
            #self.executive_lines.append(line)
            #(i, line, comment) = self.get_line_gen.next()
            #print "line =", line
        #self.executive_lines.append(line)
        #print "end of executive control",line


    def _read_case_control_deck(self):
        """
        Reads the case control deck
        @note called with recursion if an INCLUDE file is found
        """
        if not self.has_case_control_deck:
            return
        #self.log.info("reading Case Control Deck...")
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

                self.open_file(filename)
                #line = next_line
                self.case_control_lines.append(next_line)
            else:
                self.case_control_lines.append(lineUpper)

            if 'BEGIN' in lineUpper and 'BULK' in lineUpper:
                #self.log.debug('found the end of the case control deck!')
                break
        #self.log.info("finished with Case Control Deck..")

        #for line in self.case_control_lines:
            #print "** line=|%r|" %(line)

        self.caseControlDeck = CaseControlDeck(self.case_control_lines, self.log)
        self.caseControlDeck.solmap_toValue = self._solmap_to_value
        self.caseControlDeck.rsolmap_toStr = self.rsolmap_toStr

    #def _read_case_control_deck(self):
        #self.case_control_lines = []
        #(i, line, comment) = self.get_line_gen.next()
        #while 'BEGIN BULK' not in line:
            #self.case_control_lines.append(line)
            #(i, line, comment) = self.get_line_gen.next()
            #print "line =", line
        #self.case_control_lines.append(line)
        #print "end of case control", line

    def read_bulk_data_deck(self):
        n = 1
        isDone = False
        for (lines, comments) in self.gen2:
            #print "*asdf*lines =", lines
            if not lines:
                self.close_file()
                if self.active_filename:
                    continue
                break
            
            card_name = self.get_card_name(lines)
            #print "card_name =", card_name
            if card_name == 'INCLUDE':
                new_bdf_filename = get_include_filename(lines, include_dir=self.include_dir)
                #print "newfname =", newfname
                self.open_file(new_bdf_filename)
                #asfd
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
            #if n >3:
                #sys.exit('check')
            
            isDone = False

            c = ''
            if comments:  c=' comments=%s' % comments
            #print "ncard=%s lines=%s%s" % (n, lines, c)
            n += 1

            #print "**********************"
        #bbb
        print "done..."

    def get_card_name(self, lines):
        """
        Returns the name of the card defined by the provided lines
        @param self the BDF object
        @param lines the lines of the card
        @retval cardname the name of the card
        """
        #print "card_name lines=%s" % lines
        card_name = lines[0][:8].rstrip('\t, ').split(',')[0].split('\t')[0].strip('*\t ')
        #print card_name
        if len(card_name) == 0:
            return None
        assert ' ' not in card_name and len(card_name) > 0, 'card_name=|%r|\nline=|%s| in filename=%s is invalid' % (card_name, lines[0], self.active_filename)
        return card_name.upper()

    def get_card(self):
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
            while len(line)==0 or line[0] in [' ', '*', '+']:
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
            if line[0] not in [' ', '*', '+']:
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

    def open_file(self, bdf_filename):
        #print "opening self.active_filename=%s" % bdf_filename
        if bdf_filename in self.all_filenames:
            msg = 'bdf_filename=%s already found in all_filenames=%s' % self.all_filenames
            raise RuntimeError(msg)
        self.active_filename = bdf_filename
        self.all_filenames.append(bdf_filename)
        
        #if len(self.all_filenames) > 1:
        self.active_filenames.append(bdf_filename)
        #print "active_filenames =", self.active_filenames
        self.get_line_gen = self.stream_file()
        self._line_generators.append(self.get_line_gen)

    def close_file(self):
        #print ""
        #print "closing self.active_filename=%s" % self.active_filename
        #print "*active_filenames =", self.active_filenames
        if self.active_filenames:
            self.active_filename = self.active_filenames.pop()
            #print "the new active file is self.active_filename=%s" % self.active_filename
            self.get_line_gen = self._line_generators.pop()
        else:
            self.active_filename = None
            self.get_line_gen = None
        
    def stream_file(self):
        """
        builds a generator to read the BDF in a much cleaner way
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

#lines_required = 100
#gen = get_line()
#chunk = [i for i, j in zip(gen, range(lines_required))]        
                
                

if __name__ == '__main__':
    #print "***********************************"
    bdf = BDF3()
    import pyNastran
    pkg_path = pyNastran.__path__[0]
    bdfname = os.path.join(pkg_path, '..', 'models', 'solidBending.bdf')
    #print "bdfname =", bdfname
    bdf.read_bdf(bdfname)
    sys.exit('exit...')

    cards_exact = [
        ['PARAM','POST'],
        ['PARAM','PRTMAXIM'],
        ['GRID',1],
        ['GRID',2],
        ['GRID',3],
        ['GRID',4],
        ['GRID',5],
        ['GRID',6],
        ['GRID',7],
        ['GRID',8],
        ['GRID',9],
        ['GRID',10],
        ['GRID',11],
        ['GRID',12],
        ['GRID',13],
        ['GRID',14],
        ['GRID',15],
        ['GRID',16],
        ['GRID',17],
        ['GRID',18],
        ['GRID',19],
        ['GRID',20],
        ['GRID',21],
        ['GRID',22],
        ['GRID',23],
        ['GRID',24],
        ['GRID',25],
        ['GRID',26],
        ['GRID',27],
        ['GRID',28],
        ['GRID',29],
        ['GRID',30],
        ['GRID',31],
        ['GRID',32],
        ['GRID',33],
        ['GRID',34],
        ['GRID',35],
        ['GRID',36],
        ['GRID',37],
        ['GRID',38],
        ['GRID',39],
        ['GRID',40],
        ['GRID',41],
        ['GRID',42],
        ['GRID',43],
        ['GRID',44],
        ['GRID',45],
        ['GRID',46],
        ['GRID',47],
        ['GRID',48],
        ['GRID',49],
        ['GRID',50],
        ['GRID',51],
        ['GRID',52],
        ['GRID',53],
        ['GRID',54],
        ['GRID',55],
        ['GRID',56],
        ['GRID',57],
        ['GRID',58],
        ['GRID',59],
        ['GRID',60],
        ['GRID',61],
        ['GRID',62],
        ['GRID',63],
        ['GRID',64],
        ['GRID',65],
        ['GRID',66],
        ['GRID',67],
        ['GRID',68],
        ['GRID',69],
        ['GRID',70],
        ['GRID',71],
        ['GRID',72],
        ['PSOLID',1],
        ['CTETRA',1],
        ['CTETRA',2],
        ['CTETRA',3],
        ['CTETRA',4],
        ['CTETRA',5],
        ['CTETRA',6],
        ['CTETRA',7],
        ['CTETRA',8],
        ['CTETRA',9],
        ['CTETRA',10],
        ['CTETRA',11],
        ['CTETRA',12],
        ['CTETRA',13],
        ['CTETRA',14],
        ['CTETRA',15],
        ['CTETRA',16],
        ['CTETRA',17],
        ['CTETRA',18],
        ['CTETRA',19],
        ['CTETRA',20],
        ['CTETRA',21],
        ['CTETRA',22],
        ['CTETRA',23],
        ['CTETRA',24],
        ['CTETRA',25],
        ['CTETRA',26],
        ['CTETRA',27],
        ['CTETRA',28],
        ['CTETRA',29],
        ['CTETRA',30],
        ['CTETRA',31],
        ['CTETRA',32],
        ['CTETRA',33],
        ['CTETRA',34],
        ['CTETRA',35],
        ['CTETRA',36],
        ['CTETRA',37],
        ['CTETRA',38],
        ['CTETRA',39],
        ['CTETRA',40],
        ['CTETRA',41],
        ['CTETRA',42],
        ['CTETRA',43],
        ['CTETRA',44],
        ['CTETRA',45],
        ['CTETRA',46],
        ['CTETRA',47],
        ['CTETRA',48],
        ['CTETRA',49],
        ['CTETRA',50],
        ['CTETRA',51],
        ['CTETRA',52],
        ['CTETRA',53],
        ['CTETRA',54],
        ['CTETRA',55],
        ['CTETRA',56],
        ['CTETRA',57],
        ['CTETRA',58],
        ['CTETRA',59],
        ['CTETRA',60],
        ['CTETRA',61],
        ['CTETRA',62],
        ['CTETRA',63],
        ['CTETRA',64],
        ['CTETRA',65],
        ['CTETRA',66],
        ['CTETRA',67],
        ['CTETRA',68],
        ['CTETRA',69],
        ['CTETRA',70],
        ['CTETRA',71],
        ['CTETRA',72],
        ['CTETRA',73],
        ['CTETRA',74],
        ['CTETRA',75],
        ['CTETRA',76],
        ['CTETRA',77],
        ['CTETRA',78],
        ['CTETRA',79],
        ['CTETRA',80],
        ['CTETRA',81],
        ['CTETRA',82],
        ['CTETRA',83],
        ['CTETRA',84],
        ['CTETRA',85],
        ['CTETRA',86],
        ['CTETRA',87],
        ['CTETRA',88],
        ['CTETRA',89],
        ['CTETRA',90],
        ['CTETRA',91],
        ['CTETRA',92],
        ['CTETRA',93],
        ['CTETRA',94],
        ['CTETRA',95],
        ['CTETRA',96],
        ['CTETRA',97],
        ['CTETRA',98],
        ['CTETRA',99],
        ['CTETRA',100],
        ['CTETRA',101],
        ['CTETRA',102],
        ['CTETRA',103],
        ['CTETRA',104],
        ['CTETRA',105],
        ['CTETRA',106],
        ['CTETRA',107],
        ['CTETRA',108],
        ['CTETRA',109],
        ['CTETRA',110],
        ['CTETRA',111],
        ['CTETRA',112],
        ['CTETRA',113],
        ['CTETRA',114],
        ['CTETRA',115],
        ['CTETRA',116],
        ['CTETRA',117],
        ['CTETRA',118],
        ['CTETRA',119],
        ['CTETRA',120],
        ['CTETRA',121],
        ['CTETRA',122],
        ['CTETRA',123],
        ['CTETRA',124],
        ['CTETRA',125],
        ['CTETRA',126],
        ['CTETRA',127],
        ['CTETRA',128],
        ['CTETRA',129],
        ['CTETRA',130],
        ['CTETRA',131],
        ['CTETRA',132],
        ['CTETRA',133],
        ['CTETRA',134],
        ['CTETRA',135],
        ['CTETRA',136],
        ['CTETRA',137],
        ['CTETRA',138],
        ['CTETRA',139],
        ['CTETRA',140],
        ['CTETRA',141],
        ['CTETRA',142],
        ['CTETRA',143],
        ['CTETRA',144],
        ['CTETRA',145],
        ['CTETRA',146],
        ['CTETRA',147],
        ['CTETRA',148],
        ['CTETRA',149],
        ['CTETRA',150],
        ['CTETRA',151],
        ['CTETRA',152],
        ['CTETRA',153],
        ['CTETRA',154],
        ['CTETRA',155],
        ['CTETRA',156],
        ['CTETRA',157],
        ['CTETRA',158],
        ['CTETRA',159],
        ['CTETRA',160],
        ['CTETRA',161],
        ['CTETRA',162],
        ['CTETRA',163],
        ['CTETRA',164],
        ['CTETRA',165],
        ['CTETRA',166],
        ['CTETRA',167],
        ['CTETRA',168],
        ['CTETRA',169],
        ['CTETRA',170],
        ['CTETRA',171],
        ['CTETRA',172],
        ['CTETRA',173],
        ['CTETRA',174],
        ['CTETRA',175],
        ['CTETRA',176],
        ['CTETRA',177],
        ['CTETRA',178],
        ['CTETRA',179],
        ['CTETRA',180],
        ['CTETRA',181],
        ['CTETRA',182],
        ['CTETRA',183],
        ['CTETRA',184],
        ['CTETRA',185],
        ['CTETRA',186],
        ['MAT1',1],
        ['FORCE',1],
        ['FORCE',1],
        ['FORCE',1],
        ['FORCE',1],
        ['FORCE',1],
        ['FORCE',1],
        ['FORCE',1],
        ['FORCE',1],
        ['FORCE',1],
        ['FORCE',1],
        ['FORCE',1],
        ['FORCE',1],
        ['FORCE',1],
        ['FORCE',1],
        ['FORCE',1],
        ['FORCE',1],
        ['FORCE',1],
        ['FORCE',1],
        ['FORCE',1],
        ['FORCE',1],
        ['FORCE',1],
        ['FORCE',1],
        ['FORCE',1],
        ['LOAD',2],
        ['SPC1',1],
        ['SPCADD',2],
        ['SPC1',3],
        ['ENDDATA','6A7D6A3B']
    ]

    bdf = BDF3()
    bdf.open_file(bdfname)
    bdf.get_line_gen = bdf.stream_file()
    bdf.gen2 = bdf.get_card()

    bdf._read_executive_control_deck()
    bdf.read_case_control_deck()
    
    for icard, card_exact in enumerate(cards_exact):
        lines, comments = bdf.gen2.next()
        card_name = self.get_card_name(lines)

        #print 'bulk lines =', lines
        fields = to_fields(lines, card_name)
        #interp = [interpret_value(field, fields)
        #                          for field in fields]
        
        card = wipe_empty_fields([interpret_value(field, fields)
                                  for field in fields])
        #print "interp = ", card
        #print "*bulk card =", card[0:2]
        
        #print "*icard=%s card_exact=%s" % (icard, card_exact[0:2])
        #msg = 'card=%s card_exact=%s' % (card, card_exact)
        assert card[0:2] == card_exact[0:2], msg
        #print "------------------------------------------------------"
        #asf
    