# pylint: disable=W0212,C0103,W0633,W0611,W0201,C0301,R0915,R0912
# coding: utf-8
"""
Main BDF class.  Defines:
  - BDF
  - BDFDeprecated (unused for major releases)
"""
from __future__ import (nested_scopes, generators, division, absolute_import,
                        print_function, unicode_literals)
from six import string_types

#from codecs import open as codec_open
import io
import os
import sys
import warnings
import traceback

from pyNastran.utils import (object_attributes, print_bad_path)
from pyNastran.utils.dev import list_print
from pyNastran.utils.log import get_logger

from pyNastran.bdf.bdfInterface.assign_type import (integer,
    integer_or_string, string)

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
from .cards.properties.shell import PSHELL, PCOMP, PCOMPG, PSHEAR, PLPLANE
from .cards.elements.bush import CBUSH, CBUSH1D, CBUSH2D
from .cards.properties.bush import PBUSH, PBUSH1D
from .cards.elements.damper import (CVISC, CDAMP1, CDAMP2, CDAMP3, CDAMP4,
                                    CDAMP5, DamperElement)
from .cards.properties.damper import (PVISC, PDAMP, PDAMP5, PDAMPT)
from .cards.elements.rods import CROD, CONROD, CTUBE, RodElement
from .cards.elements.bars import CBAR, CBEND, LineElement, CBEAM3
from .cards.elements.beam import CBEAM
from .cards.properties.rods import PROD, PTUBE
from .cards.properties.bars import (PBAR, PBARL, )  # PBEND
from .cards.properties.beam import  PBEAM, PBEAML, PBCOMP
from .cards.elements.mass import (CONM1, CONM2, CMASS1, CMASS2, CMASS3, CMASS4,
                                  PointElement, PointMassElement)  # CMASS5
from .cards.properties.mass import (PMASS, NSM)
from .cards.aero import (AEFACT, AELINK, AELIST, AEPARM, AESTAT, AESURF,
                         AESURFS, AERO, AEROS, CSSCHD,
                         CAERO1, CAERO2, CAERO3, CAERO4, CAERO5,
                         PAERO1, PAERO2, PAERO3,
                         FLFACT, FLUTTER, GUST, MKAERO1,
                         MKAERO2, SPLINE1, SPLINE2, SPLINE4,
                         SPLINE5, TRIM)
from .cards.constraints import (SPC, SPCADD, SPCD, SPCAX, SPC1,
                                MPC, MPCADD, SUPORT1, SUPORT,
                                ConstraintObject)
from .cards.coordinateSystems import (CORD1R, CORD1C, CORD1S,
                                      CORD2R, CORD2C, CORD2S, CORD3G)
from .cards.dmig import (DEQATN, DMIG, DMI, DMIJ, DMIK, DMIJI, NastranMatrix)
from .cards.dynamic import (FREQ, FREQ1, FREQ2, FREQ4, TSTEP, TSTEPNL, NLPARM,
                            NLPCI)
from .cards.loads.loads import (LSEQ, SLOAD, DLOAD, DAREA, TLOAD1, TLOAD2,
                                RLOAD1, RLOAD2, RANDPS, RFORCE)
from .cards.loads.staticLoads import (LOAD, GRAV, ACCEL, ACCEL1, FORCE,
                                      FORCE1, FORCE2, MOMENT, MOMENT1, MOMENT2,
                                      PLOAD, PLOAD1, PLOAD2, PLOAD4, PLOADX1)

from .cards.materials import (MAT1, MAT2, MAT3, MAT4, MAT5,
                              MAT8, MAT9, MAT10, MAT11,
                              MATHP, CREEP, EQUIV)
from .cards.material_deps import (MATT1, MATT2, MATT4, MATT5, MATS1)  # TODO: add MATT3, MATT8, MATT9

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
from .cards.tables import (TABLED1, TABLED2, TABLED3, TABLED4,
                           TABLEM1, TABLEM2, TABLEM3, TABLEM4,
                           TABLES1, TABDMP1, TABLEST, TABRND1, TABRNDG, TIC)
from .cards.contact import BCRPARA, BCTADD, BCTSET, BSURF, BSURFS, BCTPARA
from .caseControlDeck import CaseControlDeck
from .bdf_Methods import BDFMethods
from .bdfInterface.getCard import GetMethods
from .bdfInterface.addCard import AddMethods
from .bdfInterface.BDF_Card import BDFCard, wipe_empty_fields
from .bdfInterface.assign_type import interpret_value
from .bdfInterface.bdf_writeMesh import WriteMesh
from .bdfInterface.crossReference import XrefMesh


#class BDFDeprecated(object):
    #def updateSolution(self, sol, method=None):
        #"""
        #.. deprecated: will be replaced in version 0.7 with
        #               :func: update_solution
        #"""
        #warnings.warn('updateSolution has been deprecated; use '
        #              'update_solution', DeprecationWarning, stacklevel=2)
        #self.update_solution(sol, method)
    #def __init__(self):
        #pass

# BDFDeprecated
class BDF(BDFMethods, GetMethods, AddMethods, WriteMesh, XrefMesh):
    """
    NASTRAN BDF Reader/Writer/Editor class.
    """
    #: this is a nastran model
    modelType = 'nastran'

    #: Flips between a dictionary based storage BDF storage method and
    #: a list based method.  Don't modify this.
    _isDict = True

    #: required for sphinx bug
    #: http://stackoverflow.com/questions/11208997/autoclass-and-instance-attributes
    #__slots__ = ['_is_dynamic_syntax']
    def __init__(self, debug=True, log=None):
        """
        Initializes the BDF object

        :param self:  the BDF object
        :param debug: used to set the logger if no logger is passed in
        :param log:   a python logging module object
        """
        assert debug in [True, False], 'debug=%r' % debug

        # file management parameters
        self._ifile = -1
        self.include_dir = ''
        self.active_filename = None
        self.active_filenames = []
        if self._isDict:
            self._stored_Is = {}
            self._stored_lines = {}
            self._stored_comments = {}
            self._line_streams = {}
            self._card_streams = {}
        else:
            self._stored_Is = []
            self._stored_lines = []
            self._stored_comments = []
            self._line_streams = []
            self._card_streams = []
        self._break_comment = None

        self._relpath = True
        if sys.version_info < (2, 6):
            self._relpath = False
        self.log = get_logger(log, 'debug' if debug else 'info')

        #: list of all read in cards - useful in determining if entire BDF
        #: was read & really useful in debugging
        self.card_count = {}
        #: stores the card_count of cards that have been rejected
        self.reject_count = {}

        #: was an ENDDATA card found
        #self.foundEndData = False

        #: allows the BDF variables to be scoped properly (i think...)
        GetMethods.__init__(self)
        AddMethods.__init__(self)
        BDFMethods.__init__(self)
        WriteMesh.__init__(self)
        XrefMesh.__init__(self)

        #: useful in debugging errors in input
        self.debug = debug

        #: flag that allows for OpenMDAO-style optimization syntax to be used
        self._is_dynamic_syntax = False

        #: lines that were rejected b/c they were for a card that isnt supported
        self.rejects = []

        #: cards that were created, but not processed
        self.reject_cards = []

        #: list of execive control deck lines
        self.executive_control_lines = []

        #: list of case control deck lines
        self.case_control_lines = []

        self.__init_attributes()

        #: the list of possible cards that will be parsed
        self.cardsToRead = set([
            'PARAM',
            'GRID', 'GRDSET', 'SPOINT',  # 'RINGAX',
            #'POINT', 'POINTAX', 'RINGAX',

            # mass
            'CONM1', 'CONM2', 'CMASS1', 'CMASS2', 'CMASS3', 'CMASS4',

            # elements
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
            'PELAS', 'PGAP', 'PFAST', 'PLPLANE',
            'PBUSH', 'PBUSH1D',
            'PDAMP', 'PDAMP5', 'PDAMPT',
            'PROD', 'PBAR', 'PBARL', 'PBEAM', 'PTUBE', 'PBEND', 'PBCOMP',
            'PBEAML',  # not fully supported
            # 'PBEAM3',

            'PSHELL', 'PCOMP', 'PCOMPG', 'PSHEAR',
            'PSOLID', 'PLSOLID', 'PVISC', 'PRAC2D', 'PRAC3D',

            # creep materials
            'CREEP',

            # materials
            'MAT1', 'MAT2', 'MAT3', 'MAT8', 'MAT9', 'MAT10', 'MAT11', 'MATHP',
            'MATT1', 'MATT2', 'MATT4', 'MATT5',  #'MATT3', 'MATT8', 'MATT9',
            'MATS1', #'MATS3', 'MATS8',
            # 'MATHE'
            #'EQUIV', # testing only, should never be activated...

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
            'GRAV', 'ACCEL', 'ACCEL1',
            'PLOAD', 'PLOAD1', 'PLOAD2', 'PLOAD4',
            'PLOADX1', 'RFORCE',

            # aero cards
            'AERO', 'AEROS', 'GUST', 'FLUTTER', 'FLFACT', 'MKAERO1', 'MKAERO2',
            'AEFACT', 'AELINK', 'AELIST', 'AEPARAM', 'AESTAT', 'AESURF',
            'CAERO1', 'CAERO2', 'CAERO3', 'CAERO4', # 'CAERO5',
            'PAERO1', 'PAERO2',  'PAERO3', # 'PAERO4', 'PAERO5',
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
            'DAREA', 'NLPARM', 'NLPCI', 'TSTEP', 'TSTEPNL',

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
            'TABLED1', 'TABLED2', 'TABLED3', 'TABLED4',
            'TABLEM1', 'TABLEM2', 'TABLEM3', 'TABLEM4',
            'TABDMP1',
            'TABLES1', 'TABLEST',
            'TABRND1', 'TABRNDG',

            # initial conditions - sid (set ID)
            #'TIC',  (in tables.py)

            #: methods - .. todo:: EIGRL not done???
            'EIGB', 'EIGR', 'EIGRL',

            #: cMethods - .. todo:: EIGC not done???
            'EIGC', 'EIGP',

            #: contact
            'BCTPARA',
            'BCRPARA', 'BCTADD', 'BCTSET', 'BSURF', 'BSURFS',

            # other
            'INCLUDE',  # '='
            'ENDDATA',
        ])

        caseControlCards = set(['FREQ', 'GUST', 'MPC', 'SPC', 'NLPARM', 'NSM',
                                'TEMP', 'TSTEPNL', 'INCLUDE'])
        self.uniqueBulkDataCards = self.cardsToRead.difference(caseControlCards)

        #: / is the delete from restart card
        self.specialCards = ['DEQATN', '/']

    def disable_cards(self, cards):
        """
        Method for removing broken cards from the reader

        :param self:  the BDF object
        :param cards: a list/set of cards that should not be read
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
            #'CTRAN' : 115,
            'CFREQ' : 118,

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
        #: the analysis type
        self.sol = None
        #: used in solution 600, method
        self.solMethod = None
        #: the line with SOL on it, marks ???
        self.iSolLine = None
        self.caseControlDeck = None

        #: store the PARAM cards
        self.params = {}

        # ------------------------------- nodes -------------------------------
        # main structural block
        #: stores SPOINT, GRID cards
        self.nodes = {}
        #: stores POINT cards
        self.points = {}
        #self.grids = {}
        self.spoints = None
        #self.epoints = None
        #: stores GRIDSET card
        self.gridSet = None

        #: stores elements (CQUAD4, CTRIA3, CHEXA8, CTETRA4, CROD, CONROD,
        #: etc.)
        self.elements = {}

        #: stores rigid elements (RBE2, RBE3, RJOINT, etc.)
        self.rigidElements = {}

        #: store CONM1, CONM2, CMASS1,CMASS2, CMASS3, CMASS4, CMASS5
        self.masses = {}
        self.properties_mass = {} # PMASS

        #: stores LOTS of propeties (PBAR, PBEAM, PSHELL, PCOMP, etc.)
        self.properties = {}

        #: stores MAT1, MAT2, MAT3, MAT8, MAT10, MAT11
        self.materials = {}

        #: defines the MAT4, MAT5
        self.thermalMaterials = {}

        #: defines the MATHE, MATHP
        self.hyperelasticMaterials = {}

        #: stores MATSx
        self.MATS1 = {}
        self.MATS3 = {}
        self.MATS8 = {}

        #: stores MATTx
        self.MATT1 = {}
        self.MATT2 = {}
        self.MATT3 = {}
        self.MATT4 = {}
        self.MATT5 = {}
        self.MATT8 = {}
        self.MATT9 = {}

        #: stores the CREEP card
        self.creepMaterials = {}

        # loads
        #: stores LOAD, FORCE, MOMENT, etc.
        self.loads = {}
        #self.gusts  = {} # Case Control GUST = 100
        #self.random = {} # Case Control RANDOM = 100

        #: stores coordinate systems
        self.coords = {0: CORD2R()}

        # --------------------------- constraints ----------------------------
        #: stores SUPORT1s
        #self.constraints = {} # suport1, anything else???
        self.suports = []  # suport, suport1

        #: stores SPCADD,SPC,SPC1,SPCD,SPCAX
        self.spcObject = ConstraintObject()
        #: stores MPCADD,MPC
        self.mpcObject = ConstraintObject()

        self.spcs = {}
        self.spcadds = {}

        self.mpcs = {}
        self.mpcadds = {}

        # --------------------------- dynamic ----------------------------
        #: stores DAREA
        self.dareas = {}

        self.pbusht = {}
        self.pdampt = {}
        self.pelast = {}

        #: frequencies
        self.frequencies = {}

        # ----------------------------------------------------------------
        #: direct matrix input - DMIG
        self.dmis = {}
        self.dmigs = {}
        self.dmijs = {}
        self.dmijis = {}
        self.dmiks = {}
        self.dequations = {}

        # ----------------------------------------------------------------
        #: SETx
        self.sets = {}
        self.asets = []
        self.bsets = []
        self.csets = []
        self.qsets = []
        #: SESETx
        self.setsSuper = {}

        # ----------------------------------------------------------------
        #: tables
        self.tables = {}
        #: randomTables
        self.randomTables = {}

        # ----------------------------------------------------------------
        #: EIGB, EIGR, EIGRL methods
        self.methods = {}
        # EIGC, EIGP methods
        self.cMethods = {}

        # ---------------------------- optimization --------------------------
        # optimization
        self.dconstrs = {}
        self.desvars = {}
        self.ddvals = {}
        self.dlinks = {}
        self.dresps = {}
        #: stores DVPREL1, DVPREL2...might change to DVxRel
        self.dvprels = {}
        self.dvmrels = {}
        self.doptprm = None

        # ------------------------- nonlinear defaults -----------------------
        #: stores NLPCI
        self.nlpcis = {}
        #: stores NLPARM
        self.nlparms = {}
        #: stores TSTEPs
        self.tsteps = {}
        #: stores TSTEPNL
        self.tstepnls = {}
        # --------------------------- aero defaults --------------------------
        # aero cards
        #: stores CAEROx
        self.caeros = {}
        #: stores PAEROx
        self.paeros = {}
        #: stores AERO
        self.aero = {}
        #: stores AEROS
        self.aeros = {}

        #: stores AEFACT
        self.aefacts = {}
        #: stores AELINK
        self.aelinks = {}
        #: stores AELIST
        self.aelists = {}
        #: stores AEPARAM
        self.aeparams = {}
        #: stores AESURF
        self.aesurfs = {}
        #: stores AESTAT
        self.aestats = {}

        #: stores GUST cards
        self.gusts = {}
        #: stores FLFACT
        self.flfacts = {}  #: .. todo:: can this be simplified ???
        #: stores FLUTTER
        self.flutters = {}
        #: mkaeros
        self.mkaeros = []
        #: store SPLINE1,SPLINE2,SPLINE4,SPLINE5
        self.splines = {}
        #: stores TRIM
        self.trims = {}

        # ------------------------- thermal defaults -------------------------
        # BCs
        #: stores thermal boundary conditions - CONV,RADBC
        self.bcs = {}  # e.g. RADBC

        #: stores PHBDY
        self.phbdys = {}
        #: stores convection properties - PCONV, PCONVM ???
        self.convectionProperties = {}

        # -------------------------contact cards-------------------------------
        self.bcrparas = {}
        self.bctadds = {}
        self.bctparas = {}
        self.bctsets = {}
        self.bsurf = {}
        self.bsurfs = {}

    def _verify_bdf(self):
        """
        Cross reference verification method.
        """
        xref = self._xref
        #for key, card in sorted(self.params.iteritems()):
            #try:
            #card._verify(xref)
        for key, card in sorted(self.nodes.iteritems()):
            try:
                card._verify(xref)
            except:
                print(str(card))
                raise
        for key, card in sorted(self.coords.iteritems()):
            try:
                card._verify(xref)
            except:
                print(str(card))
                raise
        for key, card in sorted(self.elements.iteritems()):
            try:
                card._verify(xref)
            except:
                print(str(card))
                raise
        for key, card in sorted(self.properties.iteritems()):
            try:
                card._verify(xref)
            except:
                print(str(card))
                raise
        for key, card in sorted(self.materials.iteritems()):
            try:
                card._verify(xref)
            except:
                print(str(card))
                raise

    def read_bdf(self, bdf_filename=None, include_dir=None,
                 xref=True, punch=False):
        """
        Read method for the bdf files

        :param self:         the BDF object
        :param bdf_filename: the input bdf (default=None; popup a dialog)
        :param include_dir:  the relative path to any include files
                             (default=None if no include files)
        :param xref:  should the bdf be cross referenced (default=True)
        :param punch: indicates whether the file is a punch file (default=False)

        >>> bdf = BDF()
        >>> bdf.read_bdf(bdf_filename, xref=True)
        >>> g1 = bdf.Node(1)
        >>> print(g1.Position())
        [10.0, 12.0, 42.0]
        >>> bdf.write_bdf(bdf_filename2)
        >>> print(bdf.card_stats())
        ---BDF Statistics---
        SOL 101
        bdf.nodes = 20
        bdf.elements = 10
        etc.
        """
        #self._encoding = encoding

        if bdf_filename is None:
            from pyNastran.utils.gui_io import load_file_dialog
            wildcard_wx = "Nastran BDF (*.bdf; *.dat; *.nas; *.pch)|" \
                "*.bdf;*.dat;*.nas;*.pch|" \
                "All files (*.*)|*.*"
            wildcard_qt = "Nastran BDF (*.bdf *.dat *.nas *.pch);;All files (*)"
            title = 'Please select a BDF/DAT/PCH to load'
            bdf_filename = load_file_dialog(title, wildcard_wx, wildcard_qt)
            assert bdf_filename is not None, bdf_filename

        #: the active filename (string)
        self.bdf_filename = bdf_filename
        if include_dir is None:
            include_dir = os.path.dirname(bdf_filename)

        #: the directory of the 1st BDF (include BDFs are relative to this one)
        self.include_dir = include_dir

        if not os.path.exists(bdf_filename):
            msg = 'cannot find bdf_filename=%r\n%s' % (bdf_filename, print_bad_path(bdf_filename))
            raise IOError(msg)
        if bdf_filename.lower().endswith('.pch'):
            punch = True

        try:
            self._open_file(self.bdf_filename)
            self.log.debug('---starting BDF.read_bdf of %s---' % self.bdf_filename)
            if not punch:
                self.log.debug('---reading executive & case control decks---')
                self._read_executive_control_deck()
                self._read_case_control_deck()
            else:
                self.log.debug('---skipping executive & case control decks---')

            self._read_bulk_data_deck()
            self.cross_reference(xref=xref)
            self._xref = xref
            self._cleanup_file_streams()
        except:
            self._cleanup_file_streams()
            raise
        self.log.debug('---finished BDF.read_bdf of %s---' % self.bdf_filename)

    def _cleanup_file_streams(self):
        """
        This function is required to prevent too many files being opened.
        The while loop closes them.
        """
        self._break_comment = False  # speeds up self._get_line()
        while self._get_line():
            pass
        if self._isDict:
            self._stored_Is = {}
            self._stored_lines = {}
            self._stored_comments = {}
            self._line_streams = {}
            self._card_streams = {}
        else:
            self._stored_Is = []
            self._stored_lines = []
            self._stored_comments = []
            self._line_streams = []
            self._card_streams = []
        #del self._break_comment

    def _read_executive_control_deck(self):
        """Reads the executive control deck"""
        self._break_comment = False
        lineUpper = ''
        while('CEND' not in lineUpper[:4] and 'BEGIN' not in lineUpper and
              'BULK' not in lineUpper):
            (i, line, comment) = self._get_line()
            line = line.rstrip('\n\r\t ')

            lineUpper = line.upper()
            if lineUpper == '$EXECUTIVE CONTROL DECK':
                continue  # skip this comment

            if len(line) > 0:
                self.executive_control_lines.append(line)
            lineUpper = lineUpper.split('$')[0]

        if 'CEND' in lineUpper[:4]:
            self.has_case_control_deck = True
        else:
            self.has_case_control_deck = False
            (i, line, comment) = self._get_line()   # BEGIN BULK

        sol, method, iSolLine = parse_executive_control_deck(self.executive_control_lines)
        self.update_solution(sol, method, iSolLine)

    def update_solution(self, sol, method, iSolLine):
        """
        Updates the overall solution type (e.g. 101,200,600)

        :param self:     the object pointer
        :param sol:      the solution type (101,103, etc)
        :param method:   the solution method (only for SOL=600)
        :param iSolLine: the line to put the SOL/method on
        """
        self.iSolLine = iSolLine
        # the integer of the solution type (e.g. SOL 101)
        if sol is None:
            self.sol = None
            self.solMethod = None
            return

        try:
            self.sol = int(sol)
        except ValueError:
            self.sol = self._solmap_to_value[sol]

        if self.sol == 600:
            #: solution 600 method modifier
            self.solMethod = method.strip()
            self.log.debug("sol=%s method=%s" % (self.sol, self.solMethod))
        else:  # very common
            self.solMethod = None

    def set_dynamic_syntax(self, dict_of_vars):
        """
        Uses the OpenMDAO syntax of %varName in an embedded BDF to
        update the values for an optimization study.

        :param self:         the BDF object
        :param dict_of_vars: dictionary of 7 character variable names to map.

        ::

          GRID, 1, %xVar, %yVar, %zVar

        >>> dict_of_vars = {'xVar': 1.0, 'yVar', 2.0, 'zVar':3.0}
        >>> bdf = BDF()
        >>> bdf.set_dynamic_syntax(dict_of_vars)
        >>> bdf,read_bdf(bdf_filename, xref=True)
        >>>

        .. note:: Case sensitivity is supported.
        .. note:: Variables should be 7 characters or less to fit in an
                     8-character field.
        .. warning:: Type matters!
        """
        self.dict_of_vars = {}
        assert len(dict_of_vars) > 0, 'nvars = %s' % len(dict_of_vars)
        for (key, value) in sorted(dict_of_vars.iteritems()):
            assert len(key) <= 7, ('max length for key is 7; '
                                   'len(%s)=%s' % (key, len(key)))
            assert len(key) >= 1, ('min length for key is 1; '
                                   'len(%s)=%s' % (key, len(key)))
            if not isinstance(key, string_types):
                msg = 'key=%r must be a string.  type=%s' % (key, type(key))
                raise TypeError(msg)
            self.dict_of_vars[key] = value
        self._is_dynamic_syntax = True

    def _parse_dynamic_syntax(self, key):
        """
        Applies the dynamic syntax for %varName

        :param self: the BDF object
        :param key:  the uppercased key
        :returns value: the dynamic value defined by dict_of_vars
        .. seealso:: :func: `set_dynamic_syntax`
        """
        key = key[1:].strip()
        self.log.debug("dynamic key = %r" % key)
        #self.dict_of_vars = {'P5':0.5,'ONEK':1000.}
        if key not in self.dict_of_vars:
            msg = "key=%r not found in keys=%s" % (key, self.dict_of_vars.keys())
            raise KeyError(msg)
        return self.dict_of_vars[key]

    def _is_case_control_deck(self, line):
        """
        .. todo:: not done...
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

        :param self: the BDF object
        .. note:: called with recursion if an INCLUDE file is found
        """
        self._break_comment = False
        if not self.has_case_control_deck:
            return
        line = ''
        while self.active_filename:  # keep going until finished
            #lines = []
            (i, lineIn, comment) = self._get_line()
            if lineIn is None:
                return  # file was closed
            line = lineIn.strip().split('$')[0].strip()
            lineUpper = line.upper()

            if lineUpper.startswith('INCLUDE'):
                try:
                    (i, next_line, comment) = self._get_line()
                except:
                    next_line = None

                if next_line:
                    next_line = next_line.strip().split('$')[0].strip()
                else:
                    next_line = ''
                include_lines = [line]
                while '\\' in next_line or '/' in next_line:  # more includes
                    include_lines.append(next_line)
                    # TODO: should this be next_line instead of line_next???
                    (i, line_next, comment) = self._get_line()
                    next_line = next_line.strip().split('$')[0].strip()
                self.case_control_lines.append(next_line)
                filename = get_include_filename(include_lines,
                                                include_dir=self.include_dir)
                self._open_file(filename)
            else:
                self.case_control_lines.append(lineUpper)

            if 'BEGIN' in lineUpper and ('BULK' in lineUpper or 'SUPER' in lineUpper):
                self.log.debug('found the end of the Case Control Deck!')
                break
        self.log.debug("finished with Case Control Deck...")

        #for line in self.case_control_lines:
            #print("** line=%r" % line)

        self.caseControlDeck = CaseControlDeck(self.case_control_lines, self.log)
        self.caseControlDeck.solmap_toValue = self._solmap_to_value
        self.caseControlDeck.rsolmap_toStr = self.rsolmap_toStr

    def is_reject(self, card_name):
        """
        Can the card be read

        :param self: the BDF object
        :param card_name: the card_name -> 'GRID'
        """
        if card_name.startswith('='):
            return False
        elif card_name in self.cardsToRead:
            return False
        if card_name:
            if card_name not in self.reject_count:
                self.log.info("reject card_name = %r" % card_name)
                self.reject_count[card_name] = 0
            self.reject_count[card_name] += 1
        return True

    def _open_file(self, bdf_filename):
        """
        Opens the primary bdf/dat file and all subsequent INCLUDE files.

        :param bdf_filename: the name of the bdf/dat file to open
        :returns: None

        .. note:: Doesn't allow reuse of the same bdf/dat file twice.
        """
        if len(self.active_filenames) > 1:
            bdf_filename = os.path.join(self.include_dir, str(bdf_filename))
        if not os.path.exists(bdf_filename):
            msg = 'No such bdf_filename: %r\n' % bdf_filename
            msg += 'cwd: %r' % os.getcwd()
            raise IOError(msg)

        if bdf_filename in self.active_filenames:
            msg = 'bdf_filename=%s is already active.\nactive_filenames=%s' \
                % (bdf_filename, self.active_filenames)
            raise RuntimeError(msg)
        self.log.info('opening %r' % bdf_filename)

        self._ifile += 1
        self.active_filename = bdf_filename
        self.active_filenames.append(bdf_filename)

        self._stored_Is[self._ifile] = []
        self._stored_lines[self._ifile] = []
        self._stored_comments[self._ifile] = []

        line_gen = self._stream_line()
        if self._isDict:
            self._line_streams[self._ifile] = line_gen
            self._card_streams[self._ifile] = self._stream_card(line_gen)
        else:
            self._line_streams.append(line_gen)
            self._card_streams.append(self._stream_card(line_gen))


    def _close_file(self):
        """
        handles closing the file stream and resetting the active file
        """
        self.log.info('closing %r' % self.active_filename)
        if self._ifile == 0:
            self._ifile = -1
            self.active_filename = None
            return
        if self._isDict:
            del self._stored_Is[self._ifile]
            del self._stored_lines[self._ifile]
            del self._stored_comments[self._ifile]
            del self._line_streams[self._ifile]
            del self._card_streams[self._ifile]
        else:
            self._stored_Is.pop()
            self._stored_lines.pop()
            self._stored_comments.pop()
            self._line_streams.pop()
            self._card_streams.pop()
        self._ifile -= 1

        self.active_filenames.pop()
        self.active_filename = self.active_filenames[-1]

    def _stream_card(self, line_stream):
        """
        Returns the next Bulk Data Card in the BDF

        :param self:        the BDF object
        :param line_stream: the generator for the file
        :returns lines:    the lines of the card
        :returns comment:  the comment for the card
        :returns cardname: the name of the card
        """
        for (i, line, comment) in line_stream:
            #-----------------------------------------------------------------
            # get the first line of the card
            Is = []
            lines = []
            comments = []

            comment = _clean_comment(comment)
            if comment:
                comments.append(comment)

            # If the first line is valid, continue.
            # Otherwise, keep getting lines until one isn't blank.
            if line:
                Is.append(i)
                lines.append(line)
            else:
                while len(line) == 0:
                    # you cant have an empty first line
                    (i, line, comment) = self._get_line()
                    if line:
                        break

                    comment = _clean_comment(comment)
                    if comment:
                        comments.append(comment)
                Is.append(i)
                lines.append(line)
                if comment:
                    comments.append(comment)
            assert len(lines) == 1, lines

            #-----------------------------------------------------------------
            # get another line
            try:
                (i, line, comment) = self._get_line()
            except TypeError:
                lines2 = clean_empty_lines(lines)
                yield lines2, ''.join(comments)

            #-----------------------------------------------------------------
            # We define a continuation by either a regular,
            # large field, small field, tab, or CSV formatted line.
            # Large field - a * is in the first character
            # Small field - a + or ' ' is in the first character
            #               or the line is blank
            # Tab - tab separated value; large or small formatted line
            # CSV - comma separated value; large or small formatted line

            # If the line is a continuation line, keep going.
            #in_loop = False

            Is2 = []
            lines2 = []
            comments2 = []
            while len(line) == 0 or line[0] in [' ', '*', '+', ',', '\t']:
                in_loop = True
                if len(line):
                    if Is2:
                        Is += Is2
                        lines += lines2
                        comments += comments2
                    Is.append(i)
                    lines.append(line)

                    Is2 = []
                    lines2 = []
                    comments2 = []
                    comment = _clean_comment(comment)
                    if comment:
                        comments.append(comment)
                else:
                    Is2.append(i)
                    lines2.append(line)
                    comment = _clean_comment(comment)
                    if comment:
                        comments2.append(comment)

                try:
                    (i, line, comment) = self._get_line()
                except TypeError:
                    lines2 = clean_empty_lines(lines)
                    comment = ''.join(comments+comments2)
                    yield lines2, comment

            # the extra lines we grabbed in the while loop should go on the
            # next card
            if Is2:
                self._stored_Is[self._ifile] = Is2
                self._stored_lines[self._ifile] = lines2
                self._stored_comments[self._ifile] = comments2

            #-----------------------------------------------------------------
            # We maybe got one too many lines
            if line[0] not in [' ', '*', '+', ',', '\t']:
                self._stored_Is[self._ifile].append(i)
                self._stored_lines[self._ifile].append(line)
                comment = _clean_comment(comment)
                if comment:
                    self._stored_comments[self._ifile].append(comment)

            lines2 = clean_empty_lines(lines)
            comment = ''.join(comments)
            yield lines2, comment
        return

    def _get_card_name(self, lines):
        """
        Returns the name of the card defined by the provided lines

        :param self:  the BDF object
        :param lines: the lines of the card
        :returns cardname: the name of the card
        """
        card_name = lines[0][:8].rstrip('\t, ').split(',')[0].split('\t')[0].strip('*\t ')
        if len(card_name) == 0:
            return None
        if ' ' in card_name or len(card_name) == 0:
            msg = 'card_name=%r\nline=%r in filename=%r is invalid' \
                  % (card_name, lines[0], self.active_filename)
            raise RuntimeError(msg)
        return card_name.upper()

    def _read_bulk_data_deck(self):
        """
        Parses the Bulk Data Deck

        :param self: the BDF object
        """
        self.log.debug("reading Bulk Data Deck...")
        self._break_comment = True
        #isEndData = False
        while self.active_filename:
            try:
                (lines, comment) = self._card_streams[self._ifile].next()
            except StopIteration:
                self._close_file()
                continue
            assert len(lines) > 0

            card_name = self._get_card_name(lines)
            if not isinstance(comment, string_types):
                raise TypeError('comment=%s type=%s' % (comment, type(comment)))

            if card_name == 'INCLUDE':
                bdf_filename = get_include_filename(lines, include_dir=self.include_dir)
                self._open_file(bdf_filename)
                reject = '$ INCLUDE processed:  %s\n' % bdf_filename
                if comment:
                    self.rejects.append([comment])
                self.rejects.append([reject])
                continue
            elif 'ENDDATA' in card_name:
                self._increase_card_count(card_name)
                #isEndData = True  # exits while loop
                break

            if not self.is_reject(card_name):
                # card_count is increased in add_card function
                self.add_card(lines, card_name, comment, is_list=False)
            else:
                self._increase_card_count(card_name)
                if comment:
                    self.rejects.append([comment])
                self.rejects.append(lines)

    def _increase_card_count(self, card_name):
        """
        Used for testing to check that the number of cards going in is the
        same as each time the model is read verifies proper writing of cards

        :param self:      the BDF object
        :param card_name: the card_name -> 'GRID'

        >>> bdf.read_bdf(bdf_filename)
        >>> bdf.card_count['GRID']
        50
        """
        if card_name == '':  # stupid null case
            return

        if card_name in self.card_count:
            self.card_count[card_name] += 1
        else:
            self.card_count[card_name] = 1

    def _get_line(self):
        """
        Gets the next line in the BDF from the current or sub-BDF
        """
        try:
            return self._line_streams[self._ifile].next()
        except StopIteration:
            self._close_file()
            return self._get_line()
        except KeyError:
            return

    def _stream_line(self):
        """
        Uses generators to open the file and stream the next line into
        a (line_number, comment, and line).
        """
        with open(self.active_filename, 'r') as f:
            for n, line in enumerate(f):
                line = line.rstrip('\t\r\n ')
                comment = ''
                if self._break_comment and '$' in line:
                    i = line.index('$')
                    comment = line[i:] + '\n'
                    line = line[:i].rstrip('\t ')
                yield n, line, comment

                while self._stored_lines[self._ifile]:
                    comment = ''
                    i2 = self._stored_Is[self._ifile].pop(0)
                    line2 = self._stored_lines[self._ifile].pop(0)
                    if self._stored_comments:
                        comment = ''.join(self._stored_comments[self._ifile])
                        self._stored_comments[self._ifile] = []
                    yield i2, line2, comment

    def process_card(self, card_lines):
        """
        :param self: the BDF object
        """
        card_name = self._get_card_name(card_lines)
        fields = to_fields(card_lines, card_name)
        if self._is_dynamic_syntax:
            fields = [self._parse_dynamic_syntax(field) if '%' in
                      field[0:1] else field for field in fields]
        card = wipe_empty_fields(fields)
        card[0] = card_name
        return card

    def add_card(self, card_lines, card_name, comment='', is_list=True):
        """
        Adds a card object to the BDF object.

        :param self:       the BDF object
        :param card_lines: the list of the card fields
         >>> ['GRID,1,2',]  # (is_list = False)
         >>> ['GRID',1,2,]  # (is_list = True; default)

        :param card_name: the card_name -> 'GRID'
        :param comment:   an optional the comment for the card
        :param is_list:   changes card_lines from a list of lines to
                          a list of fields
        :returns card_object: the card object representation of card

        .. note:: this is a very useful method for interfacing with the code
        .. note:: the cardObject is not a card-type object...so not a GRID
                  card or CQUAD4 object.  It's a BDFCard Object.  However,
                  you know the type (assuming a GRID), so just call the
                  *mesh.Node(nid)* to get the Node object that was just
                  created.
        .. warning:: cardObject is not returned
        """
        card_name = card_name.upper()
        self._increase_card_count(card_name)
        if card_name in ['DEQATN']:
            card_obj = card_lines
            card = card_lines
        else:
            if is_list:
                fields = card_lines
            else:
                fields = to_fields(card_lines, card_name)

            # apply OPENMDAO syntax
            if self._is_dynamic_syntax:
                fields = [self._parse_dynamic_syntax(field) if '%' in
                          field[0:1] else field for field in fields]

                card = wipe_empty_fields([interpret_value(field, fields)
                                          if field is not None
                                          else None for field in fields])
            else:  # leave everything as strings
                card = wipe_empty_fields(fields)
            card_obj = BDFCard(card)

        if self._auto_reject:
            self.reject_cards.append(card)
            print('rejecting processed auto=rejected %s' % card)
            return card_obj

        # function that gets by name the initialized object (from global scope)
        try:
            _get_cls = lambda name: globals()[name](card_obj, comment=comment)
        except Exception as e:
            if not e.args:
                e.args = ('',)
            e.args = ('%s' % e.args[0] + "\ncard = %s" % card,) + e.args[1:]
            raise
        _cls = lambda name: globals()[name]

        try:
            # cards that have their own method add_CARDNAME to add them
            if card_name in ['LSEQ', 'PHBDY', 'AERO', 'AEROS', 'AEFACT',
              'AELINK', 'AELIST', 'AEPARM', 'AESTAT', 'AESURF', 'TRIM',
              'FLUTTER', 'FLFACT', 'GUST', 'NLPARM', 'NLPCI', 'TSTEP',
              'TSTEPNL', 'SESET', 'DCONSTR', 'DESVAR', 'DDVAL', 'DLINK',
              'PARAM', 'PDAMPT', 'PELAST', 'PBUSHT']:
                try:
                    # PHBDY -> add_PHBDY
                    getattr(self, 'add_' + card_name)(_get_cls(card_name))
                except Exception as e:
                    if not e.args:
                        e.args = ('',)
                    e.args = ('%s' % e.args[0] + "\ncard = %s" % card,) + e.args[1:]
                    raise
                return card_obj

            # dictionary of cards. Key is the name of the function to add the
            # card
            # 'PCOMPG':  # hasnt been verified
            # 'MAT8':  # note there is no MAT6 or MAT7
            _cards = {
             'add_node': ['GRID'],
             'add_mass': ['CONM1', 'CONM2', 'CMASS1',
                          'CMASS2', 'CMASS3', 'CMASS4', ],
             'add_element': ['CQUAD4', 'CQUAD8', 'CQUAD', 'CQUADR', 'CQUADX',
                             'CTRIA3', 'CTRIA6', 'CTRIAR', 'CTRIAX', 'CTRIAX6',
                             'CBAR', 'CBEAM', 'CBEAM3', 'CROD', 'CONROD',
                             'CTUBE', 'CBEND', 'CELAS1', 'CELAS2', 'CELAS3',
                             'CELAS4', 'CVISC', 'CSHEAR', 'CGAP',
                             'CRAC2D', 'CRAC3D'],
             'add_damper': ['CBUSH', 'CBUSH1D', 'CFAST', 'CDAMP1',
                            'CDAMP2', 'CDAMP3', 'CDAMP4', 'CDAMP5'],
             'add_rigid_element': ['RBAR', 'RBAR1', 'RBE1', 'RBE2', 'RBE3'],
             'add_property': ['PSHELL', 'PCOMP', 'PCOMPG', 'PSHEAR', 'PSOLID',
                              'PBAR', 'PBARL', 'PBEAM', 'PBCOMP', 'PBEAML',
                              'PROD', 'PTUBE', 'PLSOLID', 'PBUSH1D', 'PBUSH',
                              'PFAST', 'PDAMP5', 'PGAP', 'PRAC2D', 'PRAC3D',
                              'PLPLANE',],

             # hasnt been verified, links up to MAT1, MAT2, MAT9 w/ same MID
             'add_creep_material': ['CREEP'],
             'add_structural_material': ['MAT1', 'MAT2', 'MAT3', 'MAT8',
                                         'MAT9', 'MAT10', 'MAT11',
                                         'EQUIV'],
             'add_hyperelastic_material': ['MATHE', 'MATHP',],
             'add_thermal_material': ['MAT4', 'MAT5'],
             'add_material_dependence': ['MATS1', 'MATS3', 'MATS8',
                                         'MATT1', 'MATT2', 'MATT3', 'MATT4',
                                         'MATT5', 'MATT8', 'MATT9'],
             'add_load': ['FORCE', 'FORCE1', 'FORCE2', 'MOMENT', 'MOMENT1',
                          'MOMENT2', 'GRAV', 'ACCEL', 'ACCEL1', 'LOAD', 'PLOAD',
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
             'add_SPLINE': ['SPLINE1', 'SPLINE2', 'SPLINE3', 'SPLINE4', 'SPLINE5'],
             'add_CAERO': ['CAERO1', 'CAERO2', 'CAERO3', 'CAERO4', 'CAERO5'],
             'add_PAERO': ['PAERO1', 'PAERO2', 'PAERO3', 'PAERO4', 'PAERO5'],
             'add_MKAERO': ['MKAERO1', 'MKAERO2'],
             'add_FREQ': ['FREQ', 'FREQ1', 'FREQ2'],
             'add_ASET': ['ASET', 'ASET1'], 'add_BSET': ['BSET', 'BSET1'],
             'add_CSET': ['CSET', 'CSET1'], 'add_QSET': ['QSET', 'QSET1'],
             'add_SET': ['SET1', 'SET3'],
             'add_DRESP': ['DRESP1', 'DRESP2'],
             'add_DVPREL': ['DVPREL1', 'DVPREL2'],
             'add_coord': ['CORD2R', 'CORD2C', 'CORD2S'],
             'add_table': ['TABLED1', 'TABLED2', 'TABLED3', 'TABLED4',
                           'TABLEM1', 'TABLEM2', 'TABLEM3', 'TABLEM4',
                           'TABLES1', 'TABLEST', 'TABDMP1'],
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
                            e.args = ('',)
                        e.args = ('%s' % e.args[0] + "\ncard = %s" % card,) + e.args[1:]
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
                        e.args = ('',)
                    e.args = ('%s' % e.args[0] + "\ncard = %s" % card,) + e.args[1:]
                    raise
                for i in _dct[card_name]:
                    if card_obj.field(i):
                        self.add_property(_cls(card_name)(card_obj, 1, comment=comment))
                return card_obj

            if card_name in ['DEQATN']:  # buggy for commas
                if comment:
                    self.rejects.append([comment])
                #print 'DEQATN:  card_obj.card=%s' %(card_obj.card)
                #self.add_DEQATN(DEQATN(card_obj)) # should be later moved to
                self.rejects.append(card)          # for loop below
            elif card_name == 'GRDSET':
                self.gridSet = GRDSET(card_obj, comment=comment)
            elif card_name == 'DOPTPRM':
                self.doptprm = DOPTPRM(card_obj, comment=comment)

            elif card_name == 'DMIG':  # not done...
                field2 = integer_or_string(card_obj, 2, 'flag')
                if field2 == 'UACCEL':  # special DMIG card
                    self.reject_cards.append(card)
                elif field2 == 0:
                    self.add_DMIG(DMIG(card_obj, comment=comment))
                else:
                    name = string(card_obj, 1, 'name')
                    try:
                        dmig = self.dmigs[name]
                    except KeyError:
                        msg = 'cannot find DMIG name=%r in names=%s' \
                            % (name, self.dmigs.keys())
                        raise KeyError(msg)
                    dmig.addColumn(card_obj, comment=comment)

            elif card_name in ['DMI', 'DMIJ', 'DMIJI', 'DMIK']:
                field2 = integer(card_obj, 2, 'flag')
                if field2 == 0:
                    getattr(self, 'add_' + card_name)(_get_cls(card_name))
                else:
                    name = string(card_obj, 1, 'name')
                    getattr(self, card_name.lower() + 's')[name].addColumn(card_obj)
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
                self.add_property_mass(PMASS(card_obj, nOffset=0, comment=comment))
                for (i, j) in enumerate([3, 5, 7]):
                    if card_obj.field(j) is not None:
                        self.add_property_mass(PMASS(card_obj, nOffset=i+1,
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

            elif card_name == 'BCRPARA':
                card = BCRPARA(card_obj, comment=comment)
                self.add_BCRPARA(card)
            elif card_name == 'BCTADD':
                card = BCTADD(card_obj, comment=comment)
                self.add_BCTADD(card)
            elif card_name == 'BCTPARA':
                card = BCTPARA(card_obj, comment=comment)
                self.add_BCTPARA(card)
            elif card_name == 'BCTSET':
                card = BCTSET(card_obj, comment=comment, sol=self.sol)
                self.add_BCTSET(card)
            elif card_name == 'BSURF':
                card = BSURF(card_obj, comment=comment)
                self.add_BSURF(card)
            elif card_name == 'BSURFS':
                card = BSURFS(card_obj, comment=comment)
                self.add_BSURFS(card)

            elif card_name == 'SPOINT':
                self.add_SPOINT(SPOINTs(card_obj, comment=comment))
            elif card_name == 'PBEAML':
                prop = PBEAML(card_obj, comment=comment)
                self.add_property(prop)

            elif 'ENDDATA' in card_name:
                raise RuntimeError('this should never happen...')
                #self.foundEndData = True
            else:
                #: ..warning:: cards with = signs in them
                #:             are not announced when they are rejected
                if '=' not in card[0]:
                    self.log.info('rejecting processed equal signed card %s' % card)
                self.reject_cards.append(card)
        except Exception as e:
            print(str(e))
            self.log.debug("card_name = %r" % card_name)
            self.log.debug("failed! Unreduced Card=%s\n" % list_print(card))
            self.log.debug("filename = %r\n" % self.bdf_filename)
            raise
        return card_obj

    def get_bdf_stats(self, return_type='string'):
        """
        Print statistics for the BDF

        :param self: the BDF object
        .. note:: if a card is not supported and not added to the proper
                  lists, this method will fail
        """
        card_stats = [
            'params', 'nodes', 'points', 'elements', 'rigidElements',
            'properties', 'materials', 'creepMaterials',
            'MATT1', 'MATT2', 'MATT3', 'MATT4', 'MATT5', 'MATT8', 'MATT9',
            'MATS1', 'MATS3', 'MATT8',
            'coords', 'mpcs', 'mpcadds',

            # dynamic cards
            'dareas', 'nlparms', 'nlpcis', 'tsteps', 'tstepnls',

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
            'debug', 'executive_control_lines',
            'case_control_lines', 'cardsToRead', 'card_count',
            'isStructured', 'uniqueBulkDataCards',
            'nCardLinesMax', 'modelType', 'includeDir',
            'cardsToWrite', 'solMethod', 'log', 'doneReading',
            'linesPack', 'lineNumbers', 'iSolLine',
            'reject_count', '_relpath', 'isOpened',
            #'foundEndData',
            'specialCards',
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
                    group_msg.append('  %-8s : %s' % (card_name, ncards))
                except KeyError:
                    group_msg.append('  %-8s : ???' % card_name)
                    #assert card_name == 'CORD2R', self.card_count
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
        if return_type == 'string':
            return '\n'.join(msg)
        else:
            return msg


def _clean_comment(comment, end=-1):
    """
    Removes specific pyNastran comment lines so duplicate lines aren't
    created.

    :param comment: the comment to possibly remove
    """
    if comment[:end] in ['$EXECUTIVE CONTROL DECK',
            '$CASE CONTROL DECK',
            '$NODES', '$SPOINTS', '$ELEMENTS',
            '$PARAMS', '$PROPERTIES', '$ELEMENTS_WITH_PROPERTIES',
            '$ELEMENTS_WITH_NO_PROPERTIES (PID=0 and unanalyzed properties)',
            '$UNASSOCIATED_PROPERTIES',
            '$MATERIALS', '$THERMAL MATERIALS',
            '$CONSTRAINTS', '$SPCs', '$MPCs', '$RIGID ELEMENTS',
            '$LOADS', '$AERO', '$AERO CONTROL SURFACES',
            '$FLUTTER', '$DYNAMIC', '$OPTIMIZATION',
            '$COORDS', '$THERMAL', '$TABLES', '$RANDOM TABLES',
            '$SETS', '$CONTACT', '$REJECTS', '$REJECT_LINES',
            '$PROPERTIES_MASS', '$MASSES']:
        comment = ''
    return comment


def to_fields(card_lines, card_name):
    """
    Converts a series of lines in a card into string versions of the field.
    Handles large, small, and CSV formatted cards.

    :param lines:     the lines of the BDF card object
    :param card_name: the card_name -> 'GRID'
    :returns fields:  the string formatted fields of the card
    """
    fields = []
    # first line
    line = card_lines.pop(0)
    if '=' in line:
        raise SyntaxError('card_name=%r\nequal signs are not supported...'
                          'line=%r' % (card_name, line))

    if '\t' in line:
        line = line.expandtabs()
        if ',' in line:
            raise SyntaxError('tabs and commas in the same line are '
                              'not supported...line=%r' % line)

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
        if ',' in line:  # csv
            new_fields = line[:72].split(',')[:9]
            for i in range(9 - len(new_fields)):
                new_fields.append('')
        else:  # standard
            new_fields = [line[0:8], line[8:16], line[16:24], line[24:32],
                          line[32:40], line[40:48], line[48:56], line[56:64],
                          line[64:72]]
        fields += new_fields
        assert len(fields) == 9

    for j, line in enumerate(card_lines): # continuation lines
        #for i, field in enumerate(fields):
        #    if field.strip() == '+':
        #        raise RuntimeError('j=%s field[%s] is a +' % (j,i))

        if '=' in line and card_name != 'EIGRL':
            raise SyntaxError('card_name=%r\nequal signs are not supported...'
                              'line=%r' % (card_name, line))
        if '\t' in line:
            line = line.expandtabs()
            if ',' in line:
                raise SyntaxError('tabs and commas in the same line are '
                                  'not supported...line=%r' % line)

        if '*' in line:  # large field
            if ',' in line:  # csv
                new_fields = line[:72].split(',')[1:5]
                for i in range(4 - len(new_fields)):
                    new_fields.append('')
            else:  # standard
                new_fields = [line[8:24], line[24:40], line[40:56], line[56:72]]
            assert len(new_fields) == 4
        else:  # small field
            if ',' in line:  # csv
                new_fields = line[:72].split(',')[1:9]
                for i in range(8 - len(new_fields)):
                    new_fields.append('')
            else:  # standard
                new_fields = [line[8:16], line[16:24], line[24:32],
                              line[32:40], line[40:48], line[48:56],
                              line[56:64], line[64:72]]
            if len(new_fields) != 8:
                nfields = len(new_fields)
                msg = 'nFields=%s new_fields=%s' % (nfields, new_fields)
                raise RuntimeError(msg)

        fields += new_fields
    return [field.strip() for field in fields]


def get_include_filename(card_lines, include_dir=''):
    """
    Parses an INCLUDE file split into multiple lines (as a list).

    :param card_lines:  the list of lines in the include card (all the lines!)
    :param include_dir: the include directory (default='')
    :returns filename:  the INCLUDE filename
    """
    cardLines2 = []
    for line in card_lines:
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
    method = None
    iSolLine = None
    for (i, eline) in enumerate(executive_control_lines):
        uline = eline.strip().upper()  # uppercase line
        uline = uline.split('$')[0].expandtabs()
        if uline[:4] in ['SOL ']:
            if ',' in uline:
                sline = uline.split(',')  # SOL 600,method
                solValue = sline[0].strip()
                method = sline[1].strip()
            else:
                solValue = uline
                method = None

            if sol is None:
                sol = solValue[3:].strip()
            else:
                raise ValueError('cannot overwrite solution existing='
                                 '|SOL %s| new =|%s|' % (sol, uline))
            iSolLine = i
    return sol, method, iSolLine


def clean_empty_lines(lines):
    """
    removes leading and trailing empty lines
    don't remove internally blank lines
    """
    found_lines = False
    if len(lines) < 2:
        return lines

    for i, line in enumerate(lines):
        if not found_lines and line:
            found_lines = True
            n1 = i
            n2 = i + 1
        elif found_lines and line:
            n2 = i + 1
    lines2 = lines[n1:n2]
    return lines2


if __name__ == '__main__':  # pragma: no cover
    from pyNastran.bdf.test.test_bdf import main
    main()
