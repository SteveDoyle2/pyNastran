# coding: utf-8
# pylint: disable=W0212,W0633,W0611,W0201,C0301,R0915,R0912
"""
Main BDF class.  Defines:
  - BDF
"""
from __future__ import (nested_scopes, generators, division, absolute_import,
                        print_function, unicode_literals)
from six import string_types, iteritems, itervalues
from collections import defaultdict

from codecs import open as codec_open
import io
import os
import sys
import traceback

from numpy import unique

from pyNastran.bdf.utils import _parse_pynastran_header
from pyNastran.utils import object_attributes, print_bad_path
from pyNastran.bdf.utils import (to_fields, get_include_filename,
                                 parse_executive_control_deck,
                                 clean_empty_lines, CardParseSyntaxError)
from pyNastran.bdf.field_writer_8 import print_card_8
from pyNastran.bdf.cards.utils import wipe_empty_fields

from pyNastran.utils import (object_attributes, print_bad_path)
from pyNastran.utils.dev import list_print
from pyNastran.utils.log import get_logger2

from pyNastran.bdf.bdfInterface.assign_type import (integer,
                                                    integer_or_string, string)

from pyNastran.bdf.cards.elements.elements import CFAST, CGAP, CRAC2D, CRAC3D
from pyNastran.bdf.cards.properties.properties import (PFAST, PGAP, PLSOLID, PSOLID,
                                                       PRAC2D, PRAC3D, PCONEAX)

from pyNastran.bdf.cards.elements.springs import (CELAS1, CELAS2, CELAS3, CELAS4,
                                                  SpringElement)
from pyNastran.bdf.cards.properties.springs import PELAS, PELAST

from pyNastran.bdf.cards.elements.solid import (CTETRA4, CTETRA10, CPYRAM5, CPYRAM13,
                                                CPENTA6, CPENTA15,
                                                CHEXA8, CHEXA20, SolidElement)
from pyNastran.bdf.cards.elements.rigid import (RBAR, RBAR1, RBE1, RBE2, RBE3, RigidElement)

from pyNastran.bdf.cards.elements.shell import (CQUAD, CQUAD4, CQUAD8, CQUADR, CQUADX,
                                                CSHEAR, CTRIA3, CTRIA6, CTRIAX,
                                                CTRIAX6, CTRIAR, ShellElement)
from pyNastran.bdf.cards.properties.shell import PSHELL, PCOMP, PCOMPG, PSHEAR, PLPLANE
from pyNastran.bdf.cards.elements.bush import CBUSH, CBUSH1D, CBUSH2D
from pyNastran.bdf.cards.properties.bush import PBUSH, PBUSH1D
from pyNastran.bdf.cards.elements.damper import (CVISC, CDAMP1, CDAMP2, CDAMP3, CDAMP4,
                                                 CDAMP5, DamperElement)
from pyNastran.bdf.cards.properties.damper import (PVISC, PDAMP, PDAMP5, PDAMPT)
from pyNastran.bdf.cards.elements.rods import CROD, CONROD, CTUBE, RodElement
from pyNastran.bdf.cards.elements.bars import CBAR, LineElement, CBEAM3 # CBEND
from pyNastran.bdf.cards.elements.beam import CBEAM
from pyNastran.bdf.cards.properties.rods import PROD, PTUBE
from pyNastran.bdf.cards.properties.bars import (PBAR, PBARL, )  # PBEND
from pyNastran.bdf.cards.properties.beam import  PBEAM, PBEAML, PBCOMP
from pyNastran.bdf.cards.elements.mass import (CONM1, CONM2, CMASS1, CMASS2, CMASS3, CMASS4,
                                               PointElement, PointMassElement)  # CMASS5
from pyNastran.bdf.cards.properties.mass import (PMASS, NSM)
from pyNastran.bdf.cards.aero import (AEFACT, AELINK, AELIST, AEPARM, AESTAT, AESURF,
                                      AESURFS, AERO, AEROS, CSSCHD,
                                      CAERO1, CAERO2, CAERO3, CAERO4, CAERO5,
                                      PAERO1, PAERO2, PAERO3,
                                      FLFACT, FLUTTER, GUST, MKAERO1,
                                      MKAERO2, SPLINE1, SPLINE2, SPLINE4,
                                      SPLINE5, TRIM)
from pyNastran.bdf.cards.constraints import (SPC, SPCADD, SPCAX, SPC1,
                                             MPC, MPCADD, SUPORT1, SUPORT, SESUP,
                                             GMSPC,
                                             ConstraintObject)
from pyNastran.bdf.cards.coordinateSystems import (CORD1R, CORD1C, CORD1S,
                                                   CORD2R, CORD2C, CORD2S, CORD3G,
                                                   GMCORD)
from pyNastran.bdf.cards.dmig import (DEQATN, DMIG, DMI, DMIJ, DMIK, DMIJI, NastranMatrix)
from pyNastran.bdf.cards.dynamic import (FREQ, FREQ1, FREQ2, FREQ4, TSTEP, TSTEPNL, NLPARM,
                                         NLPCI)
from pyNastran.bdf.cards.loads.loads import LSEQ, SLOAD, DAREA, RANDPS, RFORCE, SPCD
from pyNastran.bdf.cards.loads.dloads import DLOAD, TLOAD1, TLOAD2, RLOAD1, RLOAD2
from pyNastran.bdf.cards.loads.staticLoads import (LOAD, GRAV, ACCEL, ACCEL1, FORCE,
                                                   FORCE1, FORCE2, MOMENT, MOMENT1, MOMENT2,
                                                   PLOAD, PLOAD1, PLOAD2, PLOAD4, PLOADX1,
                                                   GMLOAD)

from pyNastran.bdf.cards.materials import (MAT1, MAT2, MAT3, MAT4, MAT5,
                                           MAT8, MAT9, MAT10, MAT11,
                                           MATHP, CREEP, EQUIV)
from pyNastran.bdf.cards.material_deps import MATT1, MATT2, MATT4, MATT5, MATS1  # TODO: add MATT3, MATT8, MATT9

from pyNastran.bdf.cards.methods import EIGB, EIGC, EIGR, EIGP, EIGRL
from pyNastran.bdf.cards.nodes import GRID, GRDSET, SPOINTs
from pyNastran.bdf.cards.optimization import (DCONSTR, DESVAR, DDVAL, DOPTPRM, DLINK,
                                              DRESP1, DRESP2, DVMREL1, DVPREL1, DVPREL2)
from pyNastran.bdf.cards.params import PARAM
from pyNastran.bdf.cards.bdf_sets import (ASET, BSET, CSET, QSET, USET,
                                          ASET1, BSET1, CSET1, QSET1, USET1,
                                          SET1, SET3, RADSET,
                                          SEBSET, SECSET, SEQSET, # SEUSET
                                          SEBSET1, SECSET1, SEQSET1, # SEUSET1
                                          SESET, SEQSEP)
from pyNastran.bdf.cards.thermal.loads import QBDY1, QBDY2, QBDY3, QHBDY, TEMP, TEMPD, QVOL
from pyNastran.bdf.cards.thermal.thermal import (CHBDYE, CHBDYG, CHBDYP, PCONV, PCONVM,
                                                 PHBDY, CONV, RADM, RADBC)
from pyNastran.bdf.cards.bdf_tables import (TABLED1, TABLED2, TABLED3, TABLED4,
                                            TABLEM1, TABLEM2, TABLEM3, TABLEM4,
                                            TABLES1, TABDMP1, TABLEST, TABRND1, TABRNDG, TIC)
from pyNastran.bdf.cards.contact import BCRPARA, BCTADD, BCTSET, BSURF, BSURFS, BCTPARA
from pyNastran.bdf.caseControlDeck import CaseControlDeck
from pyNastran.bdf.bdf_Methods import BDFMethods
from pyNastran.bdf.bdfInterface.getCard import GetMethods
from pyNastran.bdf.bdfInterface.addCard import AddMethods
from pyNastran.bdf.bdfInterface.BDF_Card import BDFCard
from pyNastran.bdf.bdfInterface.assign_type import interpret_value
from pyNastran.bdf.bdfInterface.bdf_writeMesh import WriteMesh
from pyNastran.bdf.bdfInterface.crossReference import XrefMesh, CrossReferenceError
from pyNastran.bdf.bdfInterface.attributes import BDFAttributes


class DuplicateIDsError(RuntimeError):
    pass

class BDF(BDFMethods, GetMethods, AddMethods, WriteMesh, XrefMesh, BDFAttributes):
    """
    NASTRAN BDF Reader/Writer/Editor class.
    """
    #: this is a nastran model
    modelType = 'nastran'

    #: required for sphinx bug
    #: http://stackoverflow.com/questions/11208997/autoclass-and-instance-attributes
    #__slots__ = ['_is_dynamic_syntax']
    def __init__(self, debug=True, log=None):
        """
        Initializes the BDF object

        :param self:  the BDF object
        :param debug: used to set the logger if no logger is passed in
                      True:  logs debug/info/error messages
                      False: logs info/error messages
                      None:  logs error messages
        :param log:   a python logging module object;
                      if log is set, debug is ignored and uses the
                      settings the logging object has
        """
        assert debug in [True, False, None], 'debug=%r' % debug
        self.echo = False

        # file management parameters
        self.active_filenames = []
        self.active_filename = None
        self.include_dir = ''

        self._relpath = True
        if sys.version_info < (2, 6):
            self._relpath = False

        self.log = get_logger2(log, debug)

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
        BDFAttributes.__init__(self)

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
        self.cards_to_read = set([
            'ECHOON', 'ECHOOFF',
            'PARAM',

            ## nodes
            'GRID', 'GRDSET', 'SPOINT',
            #'POINT', 'POINTAX', 'RINGAX', 'EPOINT', 'GRIDG'

            # mass
            'CONM1', 'CONM2', 'CMASS1', 'CMASS2', 'CMASS3', 'CMASS4',

            ## elements
            # springs
            'CELAS1', 'CELAS2', 'CELAS3', 'CELAS4', # 'CELAS5',
            # bushings
            'CBUSH', 'CBUSH1D', 'CBUSH2D',
            # dampers
            'CDAMP1', 'CDAMP2', 'CDAMP3', 'CDAMP4', 'CDAMP5',
            'CFAST',

            'CBAR', 'CROD', 'CTUBE', 'CBEAM', 'CBEAM3', 'CONROD', 'CBEND',
            'CTRIA3', 'CTRIA6', 'CTRIAR', 'CTRIAX', 'CTRIAX6',
            'CQUAD4', 'CQUAD8', 'CQUADR', 'CQUADX', 'CQUAD',
            'CTETRA', 'CPYRAM', 'CPENTA', 'CHEXA',
            'CSHEAR', 'CVISC', 'CRAC2D', 'CRAC3D',
            'CGAP',

            ## rigidElements
            'RBAR', 'RBAR1', 'RBE1', 'RBE2', 'RBE3',

            ## properties
            'PMASS',
            'PELAS', 'PGAP', 'PFAST', 'PLPLANE',
            'PBUSH', 'PBUSH1D',
            'PDAMP', 'PDAMP5',
            'PROD', 'PBAR', 'PBARL', 'PBEAM', 'PTUBE', 'PBEND', 'PBCOMP',
            'PBEAML',  # not fully supported
            # 'PBEAM3',

            'PSHELL', 'PCOMP', 'PCOMPG', 'PSHEAR',
            'PSOLID', 'PLSOLID', 'PVISC', 'PRAC2D', 'PRAC3D',
            # PQUAD4

            ## pdampt
            'PDAMPT',

            ## pelast
            'PELAST',

            ## creepMaterials
            'CREEP',

            ## materials
            'MAT1', 'MAT2', 'MAT3', 'MAT8', 'MAT9', 'MAT10', 'MAT11', 'MATHP',

            ## Material dependence - MATT1/MATT2/etc.
            'MATT1', 'MATT2', 'MATT4', 'MATT5',  #'MATT3', 'MATT8', 'MATT9',
            'MATS1', #'MATS3', 'MATS8',
            # 'MATHE'
            #'EQUIV', # testing only, should never be activated...

            ## thermalMaterials
            'MAT4', 'MAT5',

            ## spcs/spcadds
            'SPC', 'SPCADD', 'SPC1', 'SPCAX',
            'GMSPC',

            ## mpcs/mpcadds
            'MPC', 'MPCADD',

            ## suport/suport1/se_suport
            'SUPORT', 'SUPORT1', 'SESUP',

            ## loads
            'LOAD', 'LSEQ', 'RANDPS',
            'DLOAD', 'SLOAD', 'TLOAD1', 'TLOAD2', 'RLOAD1', 'RLOAD2',
            'FORCE', 'FORCE1', 'FORCE2',
            'MOMENT', 'MOMENT1', 'MOMENT2',
            'GRAV', 'ACCEL', 'ACCEL1',
            'PLOAD', 'PLOAD1', 'PLOAD2', 'PLOAD4',
            'PLOADX1', 'RFORCE',
            'GMLOAD', 'SPCD',
            #thermal
            'QVOL',

            # aero cards
            'AERO',  ## aero
            'AEROS',  ## aeros
            'GUST',  ## gusts
            'FLUTTER',   ## flutters
            'FLFACT',   ## flfacts
            'MKAERO1', 'MKAERO2',  ## mkaeros
            'AEFACT',   ## aefacts
            'AELINK',   ## aelinks
            'AELIST',   ## aelists
            'AEPARAM',   ## aeparams
            'AESTAT',   ## aestats
            'AESURF',  ## aesurfs
            'CAERO1', 'CAERO2', 'CAERO3', 'CAERO4',  ## caeros
            # 'CAERO5',
            'PAERO1', 'PAERO2', 'PAERO3',  ## paeros
            # 'PAERO4', 'PAERO5',
            'SPLINE1', 'SPLINE2', 'SPLINE4', 'SPLINE5',  ## splines
            #'SPLINE3', 'SPLINE6', 'SPLINE7',
            'TRIM',  ## trims
            'CSSCHD', ## csschds

            ## coords
            'CORD1R', 'CORD1C', 'CORD1S',
            'CORD2R', 'CORD2C', 'CORD2S',
            'GMCORD',

            # temperature cards
            'TEMP',  # 'TEMPD',
            'QBDY1', 'QBDY2', 'QBDY3', 'QHBDY',
            'CHBDYE', 'CHBDYG', 'CHBDYP',
            'PCONV', 'PCONVM', 'PHBDY',
            'RADBC', 'CONV',  # 'RADM',

            # ---- dynamic cards ---- #
            'DAREA',
            'NLPARM',  ## nlparms
            'NLPCI',  ## nlpcis
            'TSTEP',  ## tsteps
            'TSTEPNL',  ## tstepnls

            ## frequencies
            'FREQ', 'FREQ1', 'FREQ2',

            # direct matrix input cards
            'DMIG', 'DMIJ', 'DMIJI', 'DMIK', 'DMI',
            'DEQATN',

            # optimization cards
            'DCONSTR', 'DESVAR', 'DDVAL', 'DRESP1', 'DRESP2',
            'DVPREL1', 'DVPREL2',
            'DOPTPRM', 'DVMREL1', 'DLINK', 'DRESP3',
            #'DSCREEN',

            'SET1', 'SET3',  ## sets
            'ASET', 'ASET1',  ## asets
            'BSET', 'BSET1',  ## bsets
            'CSET', 'CSET1',  ## csets
            'QSET', 'QSET1',  ## qsets
            'USET', 'USET1',  ## usets


            # super-element sets
            'SESET',  ## se_sets

            'SEBSET', 'SEBSET1',  ## se_bsets
            'SECSET', 'SECSET1',  ## se_csets
            'SEQSET', 'SEQSET1',  ## se_qsets
            #'SEUSET', 'SEUSET1',  ## se_usets
            'SEQSEP',
            #'RADSET',

            #------------------------------------------------------------------
            ## tables
            #'DTABLE', 'TABLEHT', 'TABRNDG',
            'TABLED1', 'TABLED2', 'TABLED3', 'TABLED4',  # dynamic tables - freq/time loads
            'TABLEM1', 'TABLEM2', 'TABLEM3', 'TABLEM4',  # material tables - temperature

            # nonlinear elastic temperature dependent materials (e.g. creep)
            # sees TABLES1
            'TABLEST',
            # material tables - stress (MATS1, CREEP, MATHP)
            'TABLES1',


            ## modal damping table - tables_sdamping
            'TABDMP1',

            ## randomTables
            # PSD=func(freq); used by RANDPS card
            'TABRND1',
            # gust for aeroelastic response; used by RANDPS card
            'TABRNDG',

            #------------------------------------------------------------------

            # initial conditions - sid (set ID)
            #'TIC',  (in bdf_tables.py)

            #: methods
            'EIGB', 'EIGR', 'EIGRL',

            #: cMethods
            'EIGC', 'EIGP',

            #: contact
            'BCTPARA',  ## bctpara
            'BCRPARA',  ## bcrpara
            'BCTADD',  ## bctadds
            'BCTSET',  ## bctsets
            'BSURF',  ## bsurf
            'BSURFS',  ## bsurfs

            # other
            'INCLUDE',  # '='
            'ENDDATA',
        ])

        caseControlCards = set(['FREQ', 'GUST', 'MPC', 'SPC', 'NLPARM', 'NSM',
                                'TEMP', 'TSTEPNL', 'INCLUDE'])
        self.uniqueBulkDataCards = self.cards_to_read.difference(caseControlCards)

        #: / is the delete from restart card
        self.specialCards = ['DEQATN', '/']

    def disable_cards(self, cards):
        """
        Method for removing broken cards from the reader

        :param self:  the BDF object
        :param cards: a list/set of cards that should not be read
        """
        disableSet = set(cards)
        self.cards_to_read.difference(disableSet)

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

        # ------------------------ bad duplicates ----------------------------
        self._iparse_errors = 0
        self._nparse_errors = 0
        self._stop_on_parsing_error = True
        self._stored_parse_errors = []

        self._duplicate_nodes = []
        self._duplicate_elements = []
        self._duplicate_properties = []
        self._duplicate_materials = []
        self._duplicate_masses = []
        self._duplicate_thermal_materials = []
        self._duplicate_coords = []

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
        self.epoints = None
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
        #: stores LOAD, FORCE, FORCE1, FORCE2, MOMENT, MOMENT1, MOMENT2,
        #: PLOAD, PLOAD2, PLOAD4, SLOAD
        #: GMLOAD, SPCD,
        #: QVOL
        self.loads = {}
        self.tics = {}

        # stores DLOAD entries.
        self.dloads = {}
        # stores RLOAD1, RLOAD2, TLOAD1, TLOAD2, and ACSRCE entries.
        self.dload_entries = {}

        #self.gusts  = {} # Case Control GUST = 100
        #self.random = {} # Case Control RANDOM = 100

        #: stores coordinate systems
        self.coords = {0: CORD2R()}

        # --------------------------- constraints ----------------------------
        #: stores SUPORT1s
        #self.constraints = {} # suport1, anything else???
        self.suport = []
        self.suport1 = {}
        self.se_suport = []

        #: stores SPCADD, SPC, SPC1, SPCAX, GMSPC
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
        #: SETy
        self.sets = {}
        self.asets = []
        self.bsets = []
        self.csets = []
        self.qsets = []
        self.usets = {}

        #: SExSETy
        self.se_bsets = []
        self.se_csets = []
        self.se_qsets = []
        self.se_usets = {}
        self.se_sets = {}

        # ----------------------------------------------------------------
        #: tables
        self.tables = {}
        #: randomTables
        self.randomTables = {}
        #: TABDMP1
        self.tables_sdamping = {}

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
        self.dscreen = {}

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
        #: stores CSSCHD
        self.csschds = {}

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

        # ---------------------------------------------------------------------
        self._type_to_id_map = defaultdict(list)
        self._type_to_slot_map = None
        self._slot_to_type_map = {
            'params' : ['PARAM'],
            'nodes' : ['GRID', 'SPOINT', ], # 'RINGAX',
            'gridSet' : ['GRDSET'],
            #'POINT', 'POINTAX', 'RINGAX',

            # CMASS4 lies in the QRG
            'masses' : ['CONM1', 'CONM2', 'CMASS1', 'CMASS2', 'CMASS3'],

            'elements' : [
                'CELAS1', 'CELAS2', 'CELAS3', 'CELAS4',
                # 'CELAS5',
                'CBUSH', 'CBUSH1D', 'CBUSH2D',

                'CDAMP1', 'CDAMP2', 'CDAMP3', 'CDAMP4', 'CDAMP5',
                'CFAST',

                'CBAR', 'CROD', 'CTUBE', 'CBEAM', 'CBEAM3', 'CONROD', 'CBEND',
                'CTRIA3', 'CTRIA6', 'CTRIAR', 'CTRIAX', 'CTRIAX6',
                'CQUAD4', 'CQUAD8', 'CQUADR', 'CQUADX', 'CQUAD',
                'CTETRA', 'CPYRAM', 'CPENTA', 'CHEXA',
                'CSHEAR', 'CVISC', 'CRAC2D', 'CRAC3D',
                'CGAP',

                # thermal
                'CHBDYE', 'CHBDYG', 'CHBDYP',
            ],
            'rigidElements' : ['RBAR', 'RBAR1', 'RBE1', 'RBE2', 'RBE3',],

            'mass_properties' : ['PMASS'],
            'properties' : [
                'PELAS', 'PGAP', 'PFAST', 'PLPLANE',
                'PBUSH', 'PBUSH1D',
                'PDAMP', 'PDAMP5',
                'PROD', 'PBAR', 'PBARL', 'PBEAM', 'PTUBE', 'PBEND', 'PBCOMP',
                'PBEAML',  # not fully supported
                # 'PBEAM3',

                'PSHELL', 'PCOMP', 'PCOMPG', 'PSHEAR',
                'PSOLID', 'PLSOLID', 'PVISC', 'PRAC2D', 'PRAC3D',
            ],
            'pdampt' : ['PDAMPT',],
            'pelast' : ['PELAST',],

            # materials
            'materials' : ['MAT1', 'MAT2', 'MAT3', 'MAT8', 'MAT9', 'MAT10', 'MAT11'],
            'hyperelasticMaterials' : ['MATHP',],
            'creepMaterials' : ['CREEP'],
            'MATT1' : ['MATT1'],
            'MATT2' : ['MATT2'],
            'MATT3' : ['MATT3'],
            'MATT4' : ['MATT4'], # thermal
            'MATT5' : ['MATT5'], # thermal
            'MATT8' : ['MATT8'],
            'MATT9' : ['MATT9'],
            'MATS1' : ['MATS1'],
            'MATS3' : ['MATS3'],
            'MATS8' : ['MATS8'],

            # 'MATHE'
            #'EQUIV', # testing only, should never be activated...

            # thermal materials
            'thermalMaterials' : ['MAT4', 'MAT5',],

            # spc/mpc constraints - TODO: is this correct?
            'spcs' : ['SPC', 'SPC1', 'SPCAX', 'SPCADD', 'GMSPC'],
            #'spcadds' : ['SPCADD'],
            #'mpcadds' : ['MPCADD'],
            'mpcs' : ['MPC', 'MPCADD'],
            'suport' : ['SUPORT'],
            'suport1' : ['SUPORT1'],
            'se_suport' : ['SESUP'],

            # loads
            'loads' : [
                'LOAD', 'LSEQ', 'RANDPS',
                'FORCE', 'FORCE1', 'FORCE2',
                'MOMENT', 'MOMENT1', 'MOMENT2',
                'GRAV', 'ACCEL', 'ACCEL1',
                'PLOAD', 'PLOAD1', 'PLOAD2', 'PLOAD4',
                'PLOADX1', 'RFORCE', 'SLOAD',
                'GMLOAD', 'SPCD',

                # thermal
                'TEMP', 'QBDY1', 'QBDY2', 'QBDY3', 'QHBDY',
                'QVOL',
                ],
            'dloads' : ['DLOAD', ],
            # stores RLOAD1, RLOAD2, TLOAD1, TLOAD2, and ACSRCE entries.
            'dload_entries' : ['TLOAD1', 'TLOAD2', 'RLOAD1', 'RLOAD2',],

            # aero cards
            'aero' : ['AERO'],
            'aeros' : ['AEROS'],
            'gusts' : ['GUST'],
            'flutters' : ['FLUTTER'],
            'flfacts' : ['FLFACT'],
            'mkaeros' : ['MKAERO1', 'MKAERO2'],
            'aefacts' : ['AEFACT'],
            'aelinks' : ['AELINK'],
            'aelists' : ['AELIST'],
            'aeparams' : ['AEPARM'],
            'aesurfs' : ['AESURF'],
            'aestats' : ['AESTAT'],
            'caeros' : ['CAERO1', 'CAERO2', 'CAERO3', 'CAERO4',], # 'CAERO5',
            'paeros' : ['PAERO1', 'PAERO2', 'PAERO3',], # 'PAERO4', 'PAERO5',
            'splines' : ['SPLINE1', 'SPLINE2', 'SPLINE4', 'SPLINE5',],
            'csschds' : ['CSSCHD',],
            #'SPLINE3', 'SPLINE6', 'SPLINE7',
            'trims' : ['TRIM',],

            # coords
            'coords' : ['CORD1R', 'CORD1C', 'CORD1S',
                        'CORD2R', 'CORD2C', 'CORD2S',
                        'GMCORD'],

            # temperature cards
            # 'TEMPD',

            'phbdys' : ['PHBDY'],
            'convectionProperties' : ['PCONV', 'PCONVM'],

            # stores thermal boundary conditions
            'bcs' : ['CONV', 'RADBC', 'RADM'],


            # dynamic cards
            'dareas' : ['DAREA'],
            'nlparms' : ['NLPARM'],
            'nlpcis' : ['NLPCI'],
            'tsteps' : ['TSTEP'],
            'tstepnls' : ['TSTEPNL'],

            'frequencies' : ['FREQ', 'FREQ1', 'FREQ2'],

            # direct matrix input cards
            'dmigs' : ['DMIG'],
            'dmijs' : ['DMIJ'],
            'dmijis' : ['DMIJI'],
            'dmiks' : ['DMIK'],
            'dmis' : ['DMI'],
            'dequations' : ['DEQATN'],

            # optimization cards
            'dconstrs' : ['DCONSTR'],
            'desvars' : ['DESVAR'],
            'ddvals' : ['DDVAL'],
            'dlinks' : ['DLINK'],
            'dresps' : ['DRESP1', 'DRESP2', 'DRESP3',],
            'dvprels' : ['DVPREL1', 'DVPREL2'],
            'dvmrels' : ['DVMREL1', 'DVMREL2'],
            'doptprm' : ['DOPTPRM'],
            'dscreen' : ['DSCREEN'],


            # sets
            'asets' : ['ASET', 'ASET1'],
            'bsets' : ['BSET', 'BSET1',],
            'qsets' : ['QSET', 'QSET1'],
            'csets' : ['CSET', 'CSET1',],
            'usets' : ['USET', 'USET1',],
            'sets' : ['SET1', 'SET3',],

            # super-element sets
            'se_bsets' : ['SEBSET', 'SEBSET1'],
            'se_csets' : ['SECSET', 'SECSET1',],
            'se_qsets' : ['SEQSET', 'SEQSET1'],
            'se_usets' : ['SEUSET', 'SEQSET1'],
            'se_sets' : ['SESET'],
            # SEBSEP

            'tables' : [
                'DTABLE', 'TABLEHT', 'TABRNDG',
                'TABLED1', 'TABLED2', 'TABLED3', 'TABLED4',
                'TABLEM1', 'TABLEM2', 'TABLEM3', 'TABLEM4',
                'TABLES1', 'TABLEST',
                ],
            'tables_sdamping' : ['TABDMP1'],
            'randomTables' : ['TABRND1', 'TABRNDG',],

            # initial conditions - sid (set ID)
            ##'TIC',  (in bdf_tables.py)

            # methods
            'methods' : ['EIGB', 'EIGR', 'EIGRL',],

            # cMethods
            'cMethods' : ['EIGC', 'EIGP',],

            # contact
            'bctparas' : ['BCTPARA'],
            'bcrparas' : ['BCRPARA'],
            'bctadds' : ['BCTADD'],
            'bctsets' : ['BCTSET'],
            'bsurf' : ['BSURF'],
            'bsurfs' : ['BSURFS'],

            ## other
            #'INCLUDE',  # '='
            #'ENDDATA',
        }

    def set_error_storage(self, nparse_errors=100, stop_on_parsing_error=True,
                          nxref_errors=100, stop_on_xref_error=True):
        """
        Catch parsing errors and store them up to print them out all at once
        (not all errors are caught).

        :param self:                  the BDF object
        :param nparse_errors:         how many parse errors should be stored
                                      (default=0; all=None; no storage=0)
        :param stop_on_parsing_error: should an error be raised if there
                                      are parsing errors (default=True)
        :param nxref_errors:          how many cross-reference errors
                                      should be stored (default=0; all=None; no storage=0)
        :param stop_on_xref_error:    should an error be raised if there
                                      are cross-reference errors (default=True)
        """
        self._nparse_errors = nparse_errors
        self._nxref_errors = nxref_errors
        self._stop_on_parsing_error = stop_on_parsing_error
        self._stop_on_xref_error = stop_on_xref_error

    def read_bdf(self, bdf_filename=None, include_dir=None,
                 xref=True, punch=False, encoding=None):
        """
        Read method for the bdf files

        :param self:         the BDF object
        :param bdf_filename: the input bdf (default=None; popup a dialog)
        :param include_dir:  the relative path to any include files
                             (default=None if no include files)
        :param xref:  should the bdf be cross referenced (default=True)
        :param punch: indicates whether the file is a punch file (default=False)
        :param encoding:  the unicode encoding (default=None; system default)

        .. code-block:: python

            >>> bdf = BDF()
            >>> bdf.read_bdf(bdf_filename, xref=True)
            >>> g1 = bdf.Node(1)
            >>> print(g1.Position())
            [10.0, 12.0, 42.0]
            >>> bdf.write_card(bdf_filename2)
            >>> print(bdf.card_stats())

            ---BDF Statistics---
            SOL 101
            bdf.nodes = 20
            bdf.elements = 10
            etc.
        """
        #self.set_error_storage(nparse_errors=None, stop_on_parsing_error=True,
        #                       nxref_errors=None, stop_on_xref_error=True)
        if encoding is None:
            encoding = sys.getdefaultencoding()
        self._encoding = encoding
        if bdf_filename is None:
            from pyNastran.utils.gui_io import load_file_dialog
            wildcard_wx = "Nastran BDF (*.bdf; *.dat; *.nas; *.pch)|" \
                "*.bdf;*.dat;*.nas;*.pch|" \
                "All files (*.*)|*.*"
            wildcard_qt = "Nastran BDF (*.bdf *.dat *.nas *.pch);;All files (*)"
            title = 'Please select a BDF/DAT/PCH to load'
            bdf_filename, wildcard_level = load_file_dialog(title, wildcard_wx, wildcard_qt)
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

        self.punch = punch
        if bdf_filename.lower().endswith('.pch'):  # .. todo:: should this be removed???
            self.punch = True
        self._parse_primary_file_header(bdf_filename)

        try:
            self.log.debug('---starting BDF.read_bdf of %s---' % self.bdf_filename)
            executive_control_lines, case_control_lines, \
                bulk_data_lines = self._get_lines(self.bdf_filename, self.punch)

            self.case_control_lines = case_control_lines
            self.executive_control_lines = executive_control_lines

            sol, method, iSolLine = parse_executive_control_deck(executive_control_lines)
            self.update_solution(sol, method, iSolLine)

            self.caseControlDeck = CaseControlDeck(self.case_control_lines, self.log)
            self.caseControlDeck.solmap_toValue = self._solmap_to_value
            self.caseControlDeck.rsolmap_toStr = self.rsolmap_toStr

            cards, card_count = self.get_bdf_cards(bulk_data_lines)
            self._parse_cards(cards, card_count)
            self.pop_parse_errors()
            self.cross_reference(xref=xref)
            self._xref = xref
        except:
            raise
        self.log.debug('---finished BDF.read_bdf of %s---' % self.bdf_filename)
        self.pop_xref_errors()

    def pop_parse_errors(self):
        if self._stop_on_parsing_error:
            if self._iparse_errors == 1 and self._nparse_errors == 0:
                raise
            is_error = False
            msg = ''
            if self._duplicate_elements:
                duplicate_eids = [elem.eid for elem in self._duplicate_elements]
                uduplicate_eids = unique(duplicate_eids)
                uduplicate_eids.sort()
                msg += 'self.elements IDs are not unique=%s\n' % uduplicate_eids
                for eid in uduplicate_eids:
                    msg += 'old_element=\n%s\n' % self.elements[eid].print_repr_card()
                    msg += 'new_elements=\n'
                    for elem, eidi in zip(self._duplicate_elements, duplicate_eids):
                        if eidi == eid:
                            msg += elem.print_repr_card()
                    msg += '\n'
                    is_error = True

            if self._duplicate_properties:
                duplicate_pids = [prop.pid for prop in self._duplicate_properties]
                uduplicate_pids = unique(duplicate_pids)
                uduplicate_pids.sort()
                msg += 'self.properties IDs are not unique=%s\n' % uduplicate_pids
                for pid in duplicate_pids:
                    msg += 'old_property=\n%s\n' % self.properties[pid].print_repr_card()
                    msg += 'new_properties=\n'
                    for prop, pidi in zip(self._duplicate_properties, duplicate_pids):
                        if pidi == pid:
                            msg += prop.print_repr_card()
                    msg += '\n'
                    is_error = True

            if self._duplicate_masses:
                duplicate_eids = [elem.eid for elem in self._duplicate_masses]
                uduplicate_eids = unique(duplicate_eids)
                uduplicate_eids.sort()
                msg += 'self.massses IDs are not unique=%s\n' % uduplicate_eids
                for eid in uduplicate_eids:
                    msg += 'old_mass=\n%s\n' % self.masses[eid].print_repr_card()
                    msg += 'new_masses=\n'
                    for elem, eidi in zip(self._duplicate_masses, duplicate_eids):
                        if eidi == eid:
                            msg += elem.print_repr_card()
                    msg += '\n'
                    is_error = True

            if self._duplicate_materials:
                duplicate_mids = [mat.mid for mat in self._duplicate_materials]
                uduplicate_mids = unique(duplicate_mids)
                uduplicate_mids.sort()
                msg += 'self.materials IDs are not unique=%s\n' % uduplicate_mids
                for mid in uduplicate_mids:
                    msg += 'old_material=\n%s\n' % self.materials[mid].print_repr_card()
                    msg += 'new_materials=\n'
                    for mat, midi in zip(self._duplicate_materials, duplicate_mids):
                        if midi == mid:
                            msg += mat.print_repr_card()
                    msg += '\n'
                    is_error = True

            if self._duplicate_thermal_materials:
                duplicate_mids = [mat.mid for mat in self._duplicate_thermal_materials]
                uduplicate_mids = unique(duplicate_mids)
                uduplicate_mids.sort()
                msg += 'self.thermalMaterials IDs are not unique=%s\n' % uduplicate_mids
                for mid in uduplicate_mids:
                    msg += 'old_thermal_material=\n%s\n' % self.thermalMaterials[mid].print_repr_card()
                    msg += 'new_thermal_materials=\n'
                    for mat, midi in zip(self._duplicate_thermal_materials, duplicate_mids):
                        if midi == mid:
                            msg += mat.print_repr_card()
                    msg += '\n'
                    is_error = True

            if self._duplicate_coords:
                duplicate_cids = [coord.cid for coord in self._duplicate_coords]
                uduplicate_cids = unique(duplicate_cids)
                uduplicate_cids.sort()
                msg += 'self.coords IDs are not unique=%s\n' % uduplicate_cids
                for cid in uduplicate_cids:
                    msg += 'old_coord=\n%s\n' % self.coords[cid].print_repr_card()
                    msg += 'new_coords=\n'
                    for coord, cidi in zip(self._duplicate_coords, duplicate_cids):
                        if cidi == cid:
                            msg += coord.print_repr_card()
                    msg += '\n'
                    is_error = True
            if is_error:
                msg = 'There are dupliate cards.\n\n' + msg

            if self._stop_on_xref_error:
                msg += 'There are parsing errors.\n\n'

                for (card, an_error) in self._stored_parse_errors:
                    #msg += '%scard=%s\n' % (an_error[0], card)
                    msg += '%s\n\n'% an_error[0]
                    is_error = True

            if is_error:
                raise DuplicateIDsError(msg.rstrip())

    def pop_xref_errors(self):
        if self._stop_on_xref_error:
            if self._ixref_errors == 1 and self._nxref_errors == 0:
                raise
            is_error = False
            if self._stored_xref_errors:
                msg = 'There are cross-reference errors.\n\n'
                for (card, an_error) in self._stored_xref_errors:
                    msg += '%scard=%s\n' % (an_error[0], card)
                    is_error = True

                if is_error and self._stop_on_xref_error:
                    raise CrossReferenceError(msg.rstrip())

    def get_bdf_cards(self, bulk_data_lines):
        cards = defaultdict(list)
        card_count = defaultdict(int)
        full_comment = ''
        card_lines = []
        old_card_name = None
        backup_comment = ''
        nlines = len(bulk_data_lines)
        for i, line in enumerate(bulk_data_lines):
            #print('    backup=%r' % backup_comment)
            comment = ''
            if '$' in line:
                line, comment = line.split('$', 1)
            card_name = line.split(',', 1)[0].split('\t', 1)[0][:8].rstrip().upper()
            if card_name and card_name[0] not in ['+', '*']:
                if old_card_name:
                    #print('-------')
                    #print('applying %s' % card_name)
                    #print(full_comment)
                    #for linei in card_lines:
                        #print('  %r' % linei)
                    #print('-------')
                    if self.echo:
                        self.log.info('Reading %s:\n' % old_card_name + full_comment + ''.join(card_lines))

                    cards[old_card_name].append([full_comment, card_lines])
                    card_count[old_card_name] += 1
                    card_lines = []
                    full_comment = ''

                    if old_card_name == 'ECHOON':
                        self.echo = True
                    elif old_card_name == 'ECHOOFF':
                        self.echo = False
                old_card_name = card_name.rstrip(' *')
                if old_card_name == 'ENDDATA':
                    self.card_count['ENDDATA'] = 1
                    if nlines - i > 1:
                        self.log.debug('exiting due to ENDDATA found with %i lines left' % (nlines - i - 1))
                    #print('enddata')
                    return cards, card_count
                #print("card_name = %s" % card_name)

            comment = _clean_comment(comment)
            if line:
                card_lines.append(line)
                if backup_comment:
                    if comment:
                        full_comment += backup_comment + '$' + comment + '\n'
                    else:
                        full_comment += backup_comment
                    backup_comment = ''
                elif comment:
                    full_comment += '$' + comment + '\n'
                    backup_comment = ''

            elif comment:
                backup_comment += '$' + comment + '\n'
                #print('add backup=%r' % backup_comment)
            #elif comment:
                #backup_comment += '$' + comment + '\n'

        if card_lines:
            if self.echo:
                self.log.info('Reading %s:\n' % old_card_name + full_comment + ''.join(card_lines))
            #print('end_add %s' % card_lines)
            cards[old_card_name].append([backup_comment + full_comment, card_lines])
            card_count[old_card_name] += 1
        return cards, card_count

    #def write_test(self):
        #f = open('test.bdf', 'wb')
        #for line in executive_control_lines:
            #f.write(line)
        #for line in case_control_lines:
            #f.write(line)
        #for card_name, card in sorted(cards.iteritems()):
            ##print('---%r---' % card_name)
            #if card_name == 'ENDDATA':
                #continue
            #for cardi in card:
                #f.write('$-------------------------\n')
                #if cardi[0]:
                    #f.write(cardi[0])
                #f.write('\n'.join(cardi[1]))
                #f.write('\n')
        #if 'ENDDATA' in cards:
            #for cardi in cards['ENDDATA']:
                #if cardi[0]:
                    #f.write(cardi[0])
                #f.write('\n'.join(cardi[1]))
                #f.write('\n')

    def update_solution(self, sol, method, isol_line):
        """
        Updates the overall solution type (e.g. 101,200,600)

        :param self:      the object pointer
        :param sol:       the solution type (101,103, etc)
        :param method:    the solution method (only for SOL=600)
        :param isol_line: the line to put the SOL/method on
        """
        self.iSolLine = isol_line
        # the integer of the solution type (e.g. SOL 101)
        if sol is None:
            self.sol = None
            self.solMethod = None
            return

        try:
            self.sol = int(sol)
        except ValueError:
            try:
                self.sol = self._solmap_to_value[sol]
            except KeyError:
                self.sol = sol

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

        .. code-block:: python

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
        for (key, value) in sorted(iteritems(dict_of_vars)):
            assert len(key) <= 7, ('max length for key is 7; '
                                   'len(%s)=%s' % (key, len(key)))
            assert len(key) >= 1, ('min length for key is 1; '
                                   'len(%s)=%s' % (key, len(key)))
            if not isinstance(key, string_types):
                msg = 'key=%r must be a string.  type=%s' % (key, type(key))
                raise TypeError(msg)
            self.dict_of_vars[key] = value
        self._is_dynamic_syntax = True

    def is_reject(self, card_name):
        """
        Can the card be read.

        If the card is rejected, it's added to self.reject_count

        :param self: the BDF object
        :param card_name: the card_name -> 'GRID'
        """
        if card_name.startswith('='):
            return False
        elif card_name in self.cards_to_read:
            return False
        if card_name:
            if card_name not in self.reject_count:
                self.reject_count[card_name] = 0
            self.reject_count[card_name] += 1
        return True

    def process_card(self, card_lines):
        """
        Converts card_lines into a card.
        Considers dynamic syntax and removes empty fields

        :param self:        the BDF object
        :param card_lines:  list of strings that represent the card's lines
        :returns fields:    the parsed card's fields
        :returns card_name: the card's name

        .. code-block:: python

        >>> card_lines = ['GRID,1,,1.0,2.0,3.0,,']
        >>> model = BDF()
        >>> fields, card_name = model.process_card(card_lines)
        >>> fields
        ['GRID', '1', '', '1.0', '2.0', '3.0']
        >>> card_name
        'GRID'
        """
        card_name = self._get_card_name(card_lines)
        fields = to_fields(card_lines, card_name)
        if self._is_dynamic_syntax:
            fields = [self._parse_dynamic_syntax(field) if '%' in
                      field[0:1] else field for field in fields]
        card = wipe_empty_fields(fields)
        card[0] = card_name
        return card

    def create_card_object(self, card_lines, card_name, is_list=True):
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
        return card_obj, card

    def add_card(self, card_lines, card_name, comment='', is_list=True):
        """
        Adds a card object to the BDF object.

        :param self:       the BDF object
        :param card_lines: the list of the card fields
        :param card_name: the card_name -> 'GRID'
        :param comment:   an optional the comment for the card
        :param is_list:   changes card_lines from a list of lines to
                          a list of fields
        :returns card_object: the card object representation of card

        .. code-block:: python

          >>> model = BDF()

          # is_list is a somewhat misleading name; is it a list of card_lines?
          # where a card_line is an unparsed string
          >>> card_lines =['GRID,1,2']
          >>> comment = 'this is a comment'
          >>> model.add_card(card_lines, 'GRID', comment, is_list=True)

          # here is_list=False because it's been parsed
         >>> card = ['GRID', 1, 2,]
          >>> model.add_card(card_lines, 'GRID', comment, is_list=False)

        .. note:: this is a very useful method for interfacing with the code
        .. note:: the card_object is not a card-type object...so not a GRID
                  card or CQUAD4 object.  It's a BDFCard Object.  However,
                  you know the type (assuming a GRID), so just call the
                  *mesh.Node(nid)* to get the Node object that was just
                  created.
        """
        card_name = card_name.upper()
        card_obj, card = self.create_card_object(card_lines, card_name, is_list=is_list)

        if self._auto_reject:
            self.reject_cards.append(card)
            print('rejecting processed auto=rejected %s' % card)
            return card_obj

        if card_name == 'ECHOON':
            self.echo = True
            return
        elif card_name == 'ECHOOFF':
            self.echo = False
            return

        if self.echo:
            print(print_card_8(card_obj).rstrip())

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
                             'FLUTTER', 'FLFACT', 'GUST', 'CSSCHD',
                             'NLPARM', 'NLPCI', 'TSTEP',
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
            _cards = {
                'add_node' : ['GRID'],
                'add_mass' : ['CONM1', 'CONM2', 'CMASS1',
                              'CMASS2', 'CMASS3',
                              # CMASS4 - added later because documentation is wrong
                              ],
                'add_element' : ['CQUAD4', 'CQUAD8', 'CQUAD', 'CQUADR', 'CQUADX',
                                 'CTRIA3', 'CTRIA6', 'CTRIAR', 'CTRIAX', 'CTRIAX6',
                                 'CBAR', 'CBEAM', 'CBEAM3', 'CROD', 'CONROD',
                                 'CTUBE', 'CELAS1', 'CELAS2', 'CELAS3', # 'CBEND',
                                 'CELAS4', 'CVISC', 'CSHEAR', 'CGAP',
                                 'CRAC2D', 'CRAC3D'],
                'add_damper' : ['CBUSH', 'CBUSH1D', 'CFAST', 'CDAMP1',
                                'CDAMP2', 'CDAMP3', 'CDAMP5',
                                # CDAMP4 - is added later because the documentation is wrong
                                ],
                'add_rigid_element' : ['RBAR', 'RBAR1', 'RBE1', 'RBE2', 'RBE3'],
                'add_property' : ['PSHELL', 'PCOMP', 'PCOMPG', 'PSHEAR', 'PSOLID',
                                  'PBAR', 'PBARL', 'PBEAM', 'PBCOMP', 'PBEAML',
                                  'PROD', 'PTUBE', 'PLSOLID', 'PBUSH1D', 'PBUSH',
                                  'PFAST', 'PDAMP5', 'PGAP', 'PRAC2D', 'PRAC3D',
                                  'PLPLANE',],

                # hasnt been verified, links up to MAT1, MAT2, MAT9 w/ same MID
                'add_creep_material': ['CREEP'],

                # note there is no MAT6 or MAT7
                'add_structural_material' : ['MAT1', 'MAT2', 'MAT3', 'MAT8',
                                             'MAT9', 'MAT10', 'MAT11',
                                             'EQUIV'],
                'add_hyperelastic_material' : ['MATHE', 'MATHP',],
                'add_thermal_material' : ['MAT4', 'MAT5'],
                'add_material_dependence' : ['MATS1', 'MATS3', 'MATS8',
                                             'MATT1', 'MATT2', 'MATT3', 'MATT4',
                                             'MATT5', 'MATT8', 'MATT9'],

                'add_load' : ['FORCE', 'FORCE1', 'FORCE2', 'MOMENT', 'MOMENT1',
                              'MOMENT2', 'GRAV', 'ACCEL', 'ACCEL1', 'LOAD', 'PLOAD',
                              'PLOAD1', 'PLOAD2', 'PLOAD4', 'PLOADX1',
                              'RFORCE', 'SLOAD', 'RANDPS',
                              'GMLOAD', 'SPCD',
                              #thermal
                              'QVOL',
                              ],
                'add_dload' : ['DLOAD'],
                'add_dload_entry' : ['TLOAD1', 'TLOAD2', 'RLOAD1', 'RLOAD2',],

                'add_thermal_load' : ['TEMP', 'QBDY1', 'QBDY2', 'QBDY3', 'QHBDY'],
                'add_thermal_element' : ['CHBDYE', 'CHBDYG', 'CHBDYP'],
                'add_convection_property' : ['PCONV', 'PCONVM'],
                'add_constraint_MPC' : ['MPC', 'MPCADD'],
                'add_constraint_SPC' : ['SPC', 'SPC1', 'SPCAX', 'SPCADD',
                                        'GMSPC'],
                'add_suport' : ['SUPORT'],  # pseudo-constraint
                #'add_constraint' : ['SUPORT1'],  # pseudo-constraint - TODO: should this be add_suport1???
                'add_suport1' : ['SUPORT1'],  # pseudo-constraint
                'add_SESUP' : ['SESUP'],  # pseudo-constraint
                'add_SPLINE' : ['SPLINE1', 'SPLINE2', 'SPLINE3', 'SPLINE4', 'SPLINE5'],
                'add_CAERO' : ['CAERO1', 'CAERO2', 'CAERO3', 'CAERO4', 'CAERO5'],
                'add_PAERO' : ['PAERO1', 'PAERO2', 'PAERO3', 'PAERO4', 'PAERO5'],
                'add_MKAERO' : ['MKAERO1', 'MKAERO2'],
                'add_FREQ' : ['FREQ', 'FREQ1', 'FREQ2'],
                'add_ASET' : ['ASET', 'ASET1'],
                'add_BSET' : ['BSET', 'BSET1'],
                'add_CSET' : ['CSET', 'CSET1'],
                'add_QSET' : ['QSET', 'QSET1'],
                'add_USET' : ['USET', 'USET1'],
                'add_SET' : ['SET1', 'SET3'],
                'add_SEBSET' : ['SEBSET', 'SEBSET1'],
                'add_SECSET' : ['SECSET', 'SECSET1'],
                'add_SEQSET' : ['SEQSET', 'SEQSET1'],
                'add_SEUSET' : ['SEUSET', 'SEUSET1'],
                'add_DRESP' : ['DRESP1', 'DRESP2'],
                'add_DVPREL' : ['DVPREL1', 'DVPREL2'],
                'add_coord' : ['CORD2R', 'CORD2C', 'CORD2S', 'GMCORD'],
                'add_table' : ['TABLED1', 'TABLED2', 'TABLED3', 'TABLED4',
                               'TABLEM1', 'TABLEM2', 'TABLEM3', 'TABLEM4',
                               'TABLES1', 'TABLEST', 'TABDMP1'],
                'add_table_sdamping' : ['TABDMP1'],
                'add_random_table' : ['TABRND1', 'TABRNDG'],
                'add_method' : ['EIGB', 'EIGR', 'EIGRL'],
                'add_cmethod' : ['EIGC', 'EIGP'],
                'add_DVMREL' : ['DVMREL1'],
            }

            for func, names in iteritems(_cards):
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
            _dct = {
                'CTETRA' : (7, CTETRA4, CTETRA10),
                'CPYRAM' : (8, CPYRAM5, CPYRAM13),
                'CPENTA' : (9, CPENTA6, CPENTA15),
                'CHEXA' : (11, CHEXA8, CHEXA20),
            }
            if card_name in _dct:
                d = _dct[card_name]
                self.add_element((d[1] if card_obj.nfields == d[0]
                                  else d[2])(card_obj, comment=comment))
                return card_obj


            # masses that violate the QRG...no duplicate elements on card...lies
            _dct = {'CMASS4': (5,),}
            if card_name in _dct:
                try:
                    self.add_mass(_get_cls(card_name))
                except Exception as e:
                    if not e.args:
                        e.args = ('',)
                    e.args = ('%s' % e.args[0] + "\ncard = %s" % card,) + e.args[1:]
                    raise
                for i in _dct[card_name]:
                    if card_obj.field(i):
                        self.add_mass(_cls(card_name)(card_obj, 1, comment=comment))
                return card_obj

            # elements that violate the QRG...no duplicate elements on card...lies
            if card_name == 'CDAMP4':
                    self.add_damper(CDAMP4(card_obj, comment=comment))
                    if card_obj.field(5):
                        self.add_damper(CDAMP4(card_obj, 1, comment=comment))
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
                # should be later moved to loop below
                #self.add_DEQATN(DEQATN(card_obj))
                self.rejects.append(card)
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
                    dmig._add_column(card_obj, comment=comment)

            elif card_name in ['DMI', 'DMIJ', 'DMIJI', 'DMIK']:
                field2 = integer(card_obj, 2, 'flag')
                if field2 == 0:
                    getattr(self, 'add_' + card_name)(_get_cls(card_name))
                else:
                    name = string(card_obj, 1, 'name')
                    getattr(self, card_name.lower() + 's')[name]._add_column(card_obj)
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
            #elif card_name == 'EPOINT':
                #self.add_EPOINT(EPOINTs(card_obj, comment=comment))
            elif card_name == 'PBEAML':
                prop = PBEAML(card_obj, comment=comment)
                self.add_property(prop)

            elif 'ENDDATA' in card_name:
                raise RuntimeError('this should never happen...')
                #self.foundEndData = True
            else:
                #: .. warning:: cards with = signs in them
                #:              are not announced when they are rejected
                if '=' not in card[0]:
                    self.log.info('rejecting processed equal signed card %s' % card)
                self.reject_cards.append(card)
        except (SyntaxError, RuntimeError, AssertionError, KeyError, ValueError) as e:
            # NameErrors should be caught
            self._iparse_errors += 1
            var = traceback.format_exception_only(type(e), e)
            self._stored_parse_errors.append((card, var))
            if self._iparse_errors > self._nparse_errors:
                self.pop_parse_errors()
                #print(str(e))
                #self.log.debug("card_name = %r" % card_name)
                #self.log.debug("failed! Unreduced Card=%s\n" % list_print(card))
                #self.log.debug("filename = %r\n" % self.bdf_filename)
                #raise
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
            'suport1',
            'se_sets',
            'se_usets',

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

        # These are ignored because they're lists
        ignored_types = set([
            'spoints', 'spointi',  # singleton
            'gridSet',  # singleton

            'spcs', 'spcadds',

            'suport', 'se_suport', # suport, suport1 - list
            'doptprm',  # singleton

            # SETx - list
            'sets', 'asets', 'bsets', 'csets', 'qsets',
            'se_bsets', 'se_csets', 'se_qsets',
        ])

        ## TODO: why are some of these ignored?
        ignored_types2 = set([
            'caseControlDeck', 'spcObject2', 'mpcObject2',

            # done
            'sol', 'loads', 'mkaeros',
            'rejects', 'reject_cards',

            # not cards
            'debug', 'executive_control_lines',
            'case_control_lines', 'cards_to_read', 'card_count',
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
        for (lid, loads) in sorted(iteritems(self.loads)):
            msg.append('bdf.loads[%s]' % lid)
            groups = {}
            for load in loads:
                groups[load.type] = groups.get(load.type, 0) + 1
            for name, n in sorted(iteritems(groups)):
                msg.append('  %-8s %s' % (name + ':', n))
            msg.append('')

        # dloads
        for (lid, loads) in sorted(iteritems(self.dloads)):
            msg.append('bdf.dloads[%s]' % lid)
            groups = {}
            for load in loads:
                groups[load.type] = groups.get(load.type, 0) + 1
            for name, n in sorted(iteritems(groups)):
                msg.append('  %-8s %s' % (name + ':', n))
            msg.append('')

        for (lid, loads) in sorted(iteritems(self.dload_entries)):
            msg.append('bdf.dload_entries[%s]' % lid)
            groups = {}
            for load in loads:
                groups[load.type] = groups.get(load.type, 0) + 1
            for name, n in sorted(iteritems(groups)):
                msg.append('  %-8s %s' % (name + ':', n))
            msg.append('')

        #mkaeros
        if self.mkaeros:
            msg.append('bdf:mkaeros')
            msg.append('  %-8s %s' % ('MKAERO:', len(self.mkaeros)))

        for card_group_name in card_stats:
            card_group = getattr(self, card_group_name)
            groups = set([])

            if not isinstance(card_group, dict):
                msg = '%s is a %s; not dictionary' % (card_group_name, type(card_group))
                raise RuntimeError(msg)
            for card in itervalues(card_group):
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
            for name, counter in sorted(iteritems(self.card_count)):
                if name not in self.cards_to_read:
                    msg.append('  %-8s %s' % (name + ':', counter))
        msg.append('')
        if return_type == 'string':
            return '\n'.join(msg)
        else:
            return msg


    def _get_card_name(self, lines):
        """
        Returns the name of the card defined by the provided lines

        :param self:  the BDF object
        :param lines: the lines of the card
        :returns card_name: the name of the card
        """
        card_name = lines[0][:8].rstrip('\t, ').split(',')[0].split('\t')[0].strip('*\t ')
        if len(card_name) == 0:
            return None
        if ' ' in card_name or len(card_name) == 0:
            msg = 'card_name=%r\nline=%r in filename=%r is invalid' \
                  % (card_name, lines[0], self.active_filename)
            raise CardParseSyntaxError(msg)
        return card_name.upper()

    def _show_bad_file(self, bdf_filename):
        lines = []
        print('ENCODING - show_bad_file=%r' % self._encoding)
        with codec_open(bdf_filename, 'rU', encoding=self._encoding) as f:
            i = 0
            n = 0
            while 1:
                try:
                    line = f.readline().rstrip()
                except UnicodeDecodeError as e:
                    i0 = max([i-10, 0])
                    for i1, line in enumerate(lines[i0:i]):
                        self.log.error('lines[%i]=%r' % (i0+i1, line))
                    msg = "\n%s encoding error on line=%s of %s; not '%s'" % (self._encoding, i, bdf_filename, self._encoding)
                    raise RuntimeError(msg)
                if line:
                    n = 0
                else:
                    n += 1
                if n == 20:
                    raise RuntimeError('20 blank lines')
                i += 1
                lines.append(line)

    def _get_lines(self, bdf_filename, punch=False):

        with self._open_file(bdf_filename) as f:
            try:
                lines = f.readlines()
            except:
                self._show_bad_file(bdf_filename)

        nlines = len(lines)

        i = 0
        while i < nlines:
            line = lines[i].rstrip('\r\n\t')
            uline = line.upper()
            #print(uline.rstrip())
            if uline.startswith('INCLUDE'):
                j = i + 1
                #print('*** %s' % line)
                #bdf_filename2 = line[7:].strip(" '")
                bdf_filename2 = get_include_filename([line], include_dir=self.include_dir)
                #print('****f = %r' % bdf_filename2)


                with self._open_file(bdf_filename2) as f:
                    #print('f.name = %s' % f.name)
                    lines2 = f.readlines()

                #print('lines2 = %s' % lines2)
                nlines += len(lines2)

                #line2 = lines[j].split('$')
                #if not line2[0].isalpha():
                    #print('** %s' % line2)

                include_comment = '$ INCLUDE processed:  %s' % bdf_filename2
                #for line in lines2:
                    #print("  ?%s" % line.rstrip())
                lines = lines[:i] + [include_comment] + lines2 + lines[j:]
                #for line in lines:
                    #print("  *%s" % line.rstrip())
            i += 1

        executive_control_lines = []
        case_control_lines = []
        bulk_data_lines = []

        if punch:
            flag = 3
        else:
            flag = 1
        for i, line in enumerate(lines):
            if flag == 1:
                #line = line.upper()
                if line.upper().startswith('CEND'):
                    assert flag == 1
                    flag = 2
                executive_control_lines.append(line.rstrip())
            elif flag == 2:
                uline = line.upper()
                if 'BEGIN' in uline and ('BULK' in uline or 'SUPER' in uline):
                    assert flag == 2
                    flag = 3
                case_control_lines.append(line.rstrip())
            else:
                break
        for line in lines[i:]:
            bulk_data_lines.append(line.rstrip())
        del lines
        #for line in bulk_data_lines:
            #print(line)

        # clean comments
        executive_control_lines = [_clean_comment(line) for line in executive_control_lines]
        case_control_lines = [_clean_comment(line) for line in case_control_lines]
        return executive_control_lines, case_control_lines, bulk_data_lines

    def _increase_card_count(self, card_name, n=1):
        """
        Used for testing to check that the number of cards going in is the
        same as each time the model is read verifies proper writing of cards

        :param self:      the BDF object
        :param card_name: the card_name -> 'GRID'
        :param n:         the amount to increment by (default=1)

        >>> bdf.read_bdf(bdf_filename)
        >>> bdf.card_count['GRID']
        50
        """
        if card_name in self.card_count:
            self.card_count[card_name] += n
        else:
            self.card_count[card_name] = n

    def _open_file(self, bdf_filename):
        #import io
        bdf_filename = os.path.join(self.include_dir, str(bdf_filename))
        if not os.path.exists(bdf_filename):
            msg = 'No such bdf_filename: %r\n' % bdf_filename
            msg += 'cwd: %r\n' % os.getcwd()
            msg += print_bad_path(bdf_filename)
            raise IOError(msg)

        if bdf_filename in self.active_filenames:
            msg = 'bdf_filename=%s is already active.\nactive_filenames=%s' \
                % (bdf_filename, self.active_filenames)
            raise RuntimeError(msg)
        self.log.info('opening %r' % bdf_filename)
        self.active_filenames.append(bdf_filename)

        #print('ENCODING - _open_file=%r' % self._encoding)
        f = codec_open(bdf_filename, 'rU', encoding=self._encoding)
        return f

    def _parse_cards(self, cards, card_count):
        #print('card_count = %s' % card_count)
        for card_name, card in sorted(cards.iteritems()):
            #print('---%r---' % card_name)
            if self.is_reject(card_name):
                self.log.info('    rejecting card_name = %s' % card_name)
                for cardi in card:
                    self._increase_card_count(card_name)
                    self.rejects.append([cardi[0]] + cardi[1])
            elif card_name == 'DMIG':
                # the DMIG cards need special handling because there is no
                # guarantee that the header card comes first
                dmig_cards = []

                # read the header cards and save the other cards
                for comment, card_lines in card:
                    card_obj, card_fields = self.create_card_object(card_lines, card_name, is_list=False)
                    field2 = integer_or_string(card_obj, 2, 'flag')
                    if field2 == 0:
                        self.add_DMIG(DMIG(card_obj, comment=comment))
                    else:
                        dmig_cards.append([comment, card_obj])

                # load the matrix entries
                for comment, card_obj in dmig_cards:
                    field2 = integer_or_string(card_obj, 2, 'flag')
                    if field2 == 'UACCEL':  # special DMIG card
                        self.reject_cards.append(card)
                    #elif field2 == 0:
                        #self.add_DMIG(DMIG(card_obj, comment=comment))
                    else:
                        name = string(card_obj, 1, 'name')
                        try:
                            dmig = self.dmigs[name]
                        except KeyError:
                            msg = 'cannot find DMIG name=%r in names=%s' \
                                % (name, self.dmigs.keys())
                            raise KeyError(msg)
                        dmig._add_column(card_obj, comment=comment)

            else:
                for comment, card_lines in card:
                    assert card_lines != []
                    self.add_card(card_lines, card_name, comment=comment, is_list=False)

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

    #def _is_case_control_deck(self, line):
        #line_upper = line.upper().strip()
        #if 'CEND' in line.upper():
            #raise SyntaxError('invalid Case Control Deck card...CEND...')
        #if '=' in line_upper or ' ' in line_upper:
            #return True
        #for card in self.uniqueBulkDataCards:
            #lenCard = len(card)
            #if card in line_upper[:lenCard]:
                #return False
        #return True

    def _parse_primary_file_header(self, bdf_filename):
        """
        Extract encoding, nastran_format, and punch from the primary BDF.

        ..code-block :: python
            $ pyNastran: version=NX
            $ pyNastran: encoding=latin-1
            $ pyNastran: punch=True

        ..warning :: pyNastran lines must be at the top of the file
        """
        f = open(bdf_filename, 'r')
        check_header = True
        while check_header:
            try:
                line = f.readline()
            except:
                break

            if line.startswith('$'):
                key, value = _parse_pynastran_header(line)

                if key:
                    #print('pyNastran key=%s value=%s' % (key, value))
                    if key == 'version':
                        self.nastran_format = value
                    elif key == 'encoding':
                        self._encoding = value
                    elif key == 'punch':
                        self.punch = True if value == 'true' else False
                    else:
                        raise NotImplementedError(key)
                else:
                    break
            else:
                break
        f.close()

    def _verify_bdf(self):
        """
        Cross reference verification method.
        """
        xref = self._xref
        #for key, card in sorted(iteritems(self.params)):
            #card._verify(xref)
        for key, card in sorted(iteritems(self.nodes)):
            try:
                card._verify(xref)
            except:
                print(str(card))
                raise
        for key, card in sorted(iteritems(self.coords)):
            try:
                card._verify(xref)
            except:
                print(str(card))
                raise
        for key, card in sorted(iteritems(self.elements)):
            try:
                card._verify(xref)
            except:
                print(str(card))
                raise
        for key, card in sorted(iteritems(self.properties)):
            try:
                card._verify(xref)
            except:
                print(str(card))
                raise
        for key, card in sorted(iteritems(self.materials)):
            try:
                card._verify(xref)
            except:
                print(str(card))
                raise

IGNORE_COMMENTS = (
    '$EXECUTIVE CONTROL DECK',
    '$CASE CONTROL DECK',
    'NODES', 'SPOINTS', 'ELEMENTS',
    'PARAMS', 'PROPERTIES', 'ELEMENTS_WITH_PROPERTIES',
    'ELEMENTS_WITH_NO_PROPERTIES (PID=0 and unanalyzed properties)',
    'UNASSOCIATED_PROPERTIES',
    'MATERIALS', 'THERMAL MATERIALS',
    'CONSTRAINTS', 'SPCs', 'MPCs', 'RIGID ELEMENTS',
    'LOADS', 'AERO', 'AERO CONTROL SURFACES',
    'FLUTTER', 'DYNAMIC', 'OPTIMIZATION',
    'COORDS', 'THERMAL', 'TABLES', 'RANDOM TABLES',
    'SETS', 'CONTACT', 'REJECTS', 'REJECT_LINES',
    'PROPERTIES_MASS', 'MASSES')

def _clean_comment(comment):
    """
    Removes specific pyNastran comment lines so duplicate lines aren't
    created.

    :param comment: the comment to possibly remove
    """
    if comment == '':
        pass
    elif comment in IGNORE_COMMENTS:
        comment = ''
    elif 'pynastran' in comment.lower():
        comment = ''

    #if comment:
        #print(comment)
    return comment


def main():
    import pyNastran
    pkg_path = pyNastran.__path__[0]
    bdf_filename = os.path.abspath(os.path.join(pkg_path, '..', 'models', 'solid_bending', 'solid_bending.bdf'))
    model = BDF()
    model.read_bdf(bdf_filename, encoding='latin-1')
    node1 = model.nodes[1]

    # decode when we receive, encode on send
    note = 'helló wörld from two'  # must be same encoding as the header (utf-8)
    #note = b'helló wörld from two\n'.decode('utf-8')
    print(note)

    # this will be wrong because it's inconsistent with the header (utf-8)
    #note = b'á'.decode('latin-1')

    # this will work
    note = 'á'

    # so will this
    #note = b'á'.decode('utf-8')

    # The encoding that goes into the comment must be consisent with the local
    # file, so if your print doesn't work right, your comment will be bad too.
    #
    # If the print is correct, assuming all the characters are supported in your
    # desired encoding, it *should* work.
    print(note)

    # Comments are unmodified, so you can inadvertantly add cards/bugs.
    # A comment is a single string where all lines start with $ and end
    # with an endline character.
    node1._comment = '$ ' + note + '\n'

    # in other words, msg is a bad comment:
    msg = '$ line 1\n'
    msg += 'line 2\n'
    msg += '$ line 3\n'
    model.write_bdf('test.bdf')


if __name__ == '__main__':  # pragma: no cover
    from pyNastran.bdf.dev_unicode.test.test_bdf import main
    #main()
