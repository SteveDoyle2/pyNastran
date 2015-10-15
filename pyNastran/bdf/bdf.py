# coding: utf-8
# pylint: disable=W0212,W0633,W0201,C0301,R0915,R0912
  #W0611,
"""
Main BDF class.  Defines:
  - BDF
"""
from __future__ import (nested_scopes, generators, division, absolute_import,
                        print_function, unicode_literals)
from six import string_types, iteritems, itervalues
from collections import defaultdict

from codecs import open as codec_open
# import io
import os
import sys
import traceback

from numpy import unique

from pyNastran.bdf.utils import _parse_pynastran_header, deprecated
from pyNastran.utils import object_attributes, print_bad_path
from pyNastran.bdf.utils import (to_fields, get_include_filename,
                                 parse_executive_control_deck,
                                 CardParseSyntaxError)
from pyNastran.bdf.field_writer_8 import print_card_8
from pyNastran.bdf.cards.utils import wipe_empty_fields

from pyNastran.utils.log import get_logger2

from pyNastran.bdf.bdfInterface.assign_type import (integer,
                                                    integer_or_string, string)

from pyNastran.bdf.cards.elements.elements import CFAST, CGAP, CRAC2D, CRAC3D, PLOTEL
from pyNastran.bdf.cards.properties.properties import (PFAST, PGAP, PLSOLID, PSOLID,
                                                       PRAC2D, PRAC3D, PCONEAX)

from pyNastran.bdf.cards.elements.springs import (CELAS1, CELAS2, CELAS3, CELAS4,)
from pyNastran.bdf.cards.properties.springs import PELAS, PELAST

from pyNastran.bdf.cards.elements.solid import (CTETRA4, CTETRA10, CPYRAM5, CPYRAM13,
                                                CPENTA6, CPENTA15,
                                                CHEXA8, CHEXA20)
from pyNastran.bdf.cards.elements.rigid import RBAR, RBAR1, RBE1, RBE2, RBE3

from pyNastran.bdf.cards.elements.shell import (CQUAD, CQUAD4, CQUAD8, CQUADR, CQUADX,
                                                CSHEAR, CTRIA3, CTRIA6, CTRIAX,
                                                CTRIAX6, CTRIAR)
from pyNastran.bdf.cards.properties.shell import PSHELL, PCOMP, PCOMPG, PSHEAR, PLPLANE
from pyNastran.bdf.cards.elements.bush import CBUSH, CBUSH1D, CBUSH2D
from pyNastran.bdf.cards.properties.bush import PBUSH, PBUSH1D, PBUSHT
from pyNastran.bdf.cards.elements.damper import (CVISC, CDAMP1, CDAMP2, CDAMP3, CDAMP4,
                                                 CDAMP5)
from pyNastran.bdf.cards.properties.damper import PVISC, PDAMP, PDAMP5, PDAMPT
from pyNastran.bdf.cards.elements.rods import CROD, CONROD, CTUBE
from pyNastran.bdf.cards.elements.bars import CBAR, CBEAM3, CBEND
from pyNastran.bdf.cards.elements.beam import CBEAM
from pyNastran.bdf.cards.properties.rods import PROD, PTUBE
from pyNastran.bdf.cards.properties.bars import PBAR, PBARL  # PBEND
from pyNastran.bdf.cards.properties.beam import PBEAM, PBEAML, PBCOMP
from pyNastran.bdf.cards.elements.mass import (CONM1, CONM2, CMASS1, CMASS2, CMASS3, CMASS4,
                                               )  # CMASS5
from pyNastran.bdf.cards.properties.mass import (PMASS, NSM)
from pyNastran.bdf.cards.aero import (AECOMP, AEFACT, AELINK, AELIST, AEPARM, AESTAT,
                                      AESURF, AESURFS, AERO, AEROS, CSSCHD,
                                      CAERO1, CAERO2, CAERO3, CAERO4, CAERO5,
                                      PAERO1, PAERO2, PAERO3, PAERO5, # PAERO4
                                      MONPNT1,
                                      FLFACT, FLUTTER, GUST, MKAERO1,
                                      MKAERO2, SPLINE1, SPLINE2, SPLINE3, SPLINE4,
                                      SPLINE5, TRIM)
from pyNastran.bdf.cards.constraints import (SPC, SPCADD, SPCAX, SPC1,
                                             MPC, MPCADD, SUPORT1, SUPORT, SESUP,
                                             GMSPC, ConstraintObject)
from pyNastran.bdf.cards.coordinateSystems import (CORD1R, CORD1C, CORD1S,
                                                   CORD2R, CORD2C, CORD2S, CORD3G,
                                                   GMCORD)
from pyNastran.bdf.cards.dmig import DMIG, DMI, DMIJ, DMIK, DMIJI
from pyNastran.bdf.cards.deqatn import DEQATN
from pyNastran.bdf.cards.dynamic import (DELAY, FREQ, FREQ1, FREQ2, FREQ4, TSTEP, TSTEPNL,
                                         NLPARM, NLPCI, TF)
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
from pyNastran.bdf.cards.nodes import GRID, GRDSET, SPOINTs, EPOINTs
from pyNastran.bdf.cards.optimization import (DCONADD, DCONSTR, DESVAR, DDVAL, DOPTPRM, DLINK,
                                              DRESP1, DRESP2, DRESP3,
                                              DVMREL1,
                                              DVPREL1, DVPREL2)
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
                                            TABLES1, TABDMP1, TABLEST, TABRND1, TABRNDG, TIC,
                                            DTABLE)
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
from pyNastran.bdf.field_writer_16 import print_field_16

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

        Parameters
        ----------
        self :  BDF()
            the BDF object
        debug : bool/None
            used to set the logger if no logger is passed in
                True:  logs debug/info/error messages
                False: logs info/error messages
                None:  logs error messages
        log : logging module object / None
            if log is set, debug is ignored and uses the
            settings the logging object has
        """
        assert debug in [True, False, None], 'debug=%r' % debug
        self.echo = False

        # file management parameters
        self.active_filenames = []
        self.active_filename = None
        self.include_dir = ''
        self._encoding = None
        self.punch = None

        # this flag will be flipped to True someday (and then removed), but
        # doesn't support 100% of cards yet.  It enables a new method for card
        # parsing.
        #
        # 80.3 seconds -> 67.2 seconds for full_bay model
        # (multiple BDF passes among other things)
        self._fast_add = True

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
            'GRID', 'GRDSET', 'SPOINT', 'EPOINT',
            #'POINT', 'POINTAX', 'RINGAX', 'GRIDG'

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

            ## plotels
            'PLOTEL',

            ## properties
            'PMASS',
            'PELAS', 'PGAP', 'PFAST', 'PLPLANE',
            'PBUSH', 'PBUSH1D',
            'PDAMP', 'PDAMP5',
            'PROD', 'PBAR', 'PBARL', 'PBEAM', 'PTUBE', 'PBCOMP', # 'PBEND',
            'PBEAML',  # not fully supported
            # 'PBEAM3',

            'PSHELL', 'PCOMP', 'PCOMPG', 'PSHEAR',
            'PSOLID', 'PLSOLID', 'PVISC', 'PRAC2D', 'PRAC3D',
            # PQUAD4

            ## pdampt
            'PDAMPT',

            ## pelast
            'PELAST',

            ## pbusht
            'PBUSHT',

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
            'AECOMP',   ## aecomps
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
            'MONPNT1',                                   ## monitor_points
            'SPLINE1', 'SPLINE2', 'SPLINE4', 'SPLINE5',  ## splines
            #'SPLINE3', 'SPLINE6', 'SPLINE7',
            'TRIM',  ## trims
            'CSSCHD', ## csschds

            ## coords
            'CORD1R', 'CORD1C', 'CORD1S',
            'CORD2R', 'CORD2C', 'CORD2S',
            'GMCORD',

            # temperature cards
            'TEMP', 'TEMPD',
            'QBDY1', 'QBDY2', 'QBDY3', 'QHBDY',
            'CHBDYE', 'CHBDYG', 'CHBDYP',
            'PCONV', 'PCONVM', 'PHBDY',
            'RADBC', 'CONV',  # 'RADM',

            # ---- dynamic cards ---- #
            'DAREA',
            'DELAY',  ## delays
            'NLPARM',  ## nlparms
            'NLPCI',  ## nlpcis
            'TSTEP',  ## tsteps
            'TSTEPNL',  ## tstepnls
            'TF',  ## transfer_functions

            ## frequencies
            'FREQ', 'FREQ1', 'FREQ2', #'FREQ4',

            # direct matrix input cards
            'DMIG', 'DMIJ', 'DMIJI', 'DMIK', 'DMI',


            # optimization cards
            'DEQATN', 'DTABLE',
            'DCONSTR', 'DESVAR', 'DDVAL', 'DRESP1', 'DRESP2', 'DRESP3',
            'DVPREL1', 'DVPREL2',
            'DVMREL1',
            'DOPTPRM', 'DLINK', 'DCONADD',
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
            #'TABLEHT', 'TABRNDG',
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

        case_control_cards = set(['FREQ', 'GUST', 'MPC', 'SPC', 'NLPARM', 'NSM',
                                  'TEMP', 'TSTEPNL', 'INCLUDE'])
        self._unique_bulk_data_cards = self.cards_to_read.difference(case_control_cards)

        #: / is the delete from restart card
        self.specialCards = ['DEQATN', '/']
        self._make_card_parser()

    def deprecated(self, old_name, new_name, deprecated_version):
        """deprecates methods"""
        return deprecated(old_name, new_name, deprecated_version, levels=[0, 1, 2])

    def save_object(self, obj_filename='model.obj'):
        """
        ..warning:: doesn't work right
        """
        #import cPickle as pickle
        import pickle
        del self.log
        del self.spcObject
        del self.mpcObject
        self.case_control_lines = str(self.case_control_deck).split('\n')
        del self.case_control_deck
        self.uncross_reference()
        obj_file = open(obj_filename, "w")
        pickle.dump(self, obj_file)
        obj_file.close()

    def load_object(self, obj_filename='model.obj'):
        """
        ..warning:: doesn't work right
        """
        #import cPickle as pickle
        import pickle
        #del self.log
        #del self.spcObject
        #del self.mpcObject
        #lines = print(self.case_control_deck)
        #self.case_control_lines = lines.split('\n')
        #del self.case_control_deck
        #self.uncross_reference()
        #import types
        obj_file = open(obj_filename, "r")
        obj = pickle.load(obj_file)
        obj_file.close()

        keys_to_skip = [
            'case_control_deck', 'caseControlDeck',
            'log', 'mpcObject', 'spcObject',
            'node_ids', 'coord_ids', 'element_ids', 'property_ids',
            'material_ids', 'caero_ids', 'is_long_ids',
            'nnodes', 'ncoords', 'nelements', 'nproperties',
            'nmaterials', 'ncaeros',
        ]
        for key in object_attributes(self, mode="all", keys_to_skip=keys_to_skip):
            if key.startswith('__') and key.endswith('__'):
                continue

            val = getattr(obj, key)
            #print(key)
            #if isinstance(val, types.FunctionType):
                #continue
            setattr(self, key, val)

        self.log.debug('done loading!')

    def disable_cards(self, cards):
        """
        Method for removing broken cards from the reader

        Paramters
        ---------
        self : BDF()
            the BDF object
        cards : List[str]; Set[str]
            a list/set of cards that should not be read

        .. python ::

            bdfModel.disable_cards(['DMIG', 'PCOMP'])
        """
        if isinstance(cards, string_types):
            disable_set = set([cards])
        else:
            disable_set = set(cards)
        self.cards_to_read = self.cards_to_read.difference(disable_set)

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
        self._stop_on_duplicate_error = True
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
        self.case_control_deck = None

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
        #: stores PLOTELs
        self.plotels = {}

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
        self._dmig_temp = defaultdict(list)

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
        self.dconadds = {}
        self.dconstrs = {}
        self.desvars = {}
        self.ddvals = {}
        self.dlinks = {}
        self.dresps = {}

        self.dtable = None
        self.dequations = {}

        #: stores DVPREL1, DVPREL2...might change to DVxRel
        self.dvprels = {}
        self.dvmrels = {}
        self.dvcrels = {}
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
        #: stores TF
        self.transfer_functions = {}
        #: stores DELAY
        self.delays = {}

        # --------------------------- aero defaults --------------------------
        # aero cards
        #: stores CAEROx
        self.caeros = {}
        #: stores PAEROx
        self.paeros = {}
        # stores MONPNT1
        self.monitor_points = {}
        #: stores AERO
        self.aero = {}
        #: stores AEROS
        self.aeros = {}

        #: stores AECOMP
        self.aecomps = {}
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
        #: stores TEMPD
        self.tempds = {}

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
            'nodes' : ['GRID', 'SPOINT', 'EPOINT'], # 'RINGAX',
            'gridSet' : ['GRDSET'],
            #'POINT', 'POINTAX', 'RINGAX',

            # CMASS4 lies in the QRG
            'masses' : ['CONM1', 'CONM2', 'CMASS1', 'CMASS2', 'CMASS3', 'CMASS4'],

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
            'plotels' : ['PLOTEL',],

            'properties_mass' : ['PMASS'],
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
            'pbusht' : ['PBUSHT',],

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
            'aecomps' : ['AECOMP'],
            'aefacts' : ['AEFACT'],
            'aelinks' : ['AELINK'],
            'aelists' : ['AELIST'],
            'aeparams' : ['AEPARM'],
            'aesurfs' : ['AESURF', 'AESURFS'],
            'aestats' : ['AESTAT'],
            'caeros' : ['CAERO1', 'CAERO2', 'CAERO3', 'CAERO4',], # 'CAERO5',
            'paeros' : ['PAERO1', 'PAERO2', 'PAERO3', 'PAERO5'], # 'PAERO4',
            'monitor_points' : ['MONPNT1'],
            'splines' : ['SPLINE1', 'SPLINE2', 'SPLINE4', 'SPLINE5',],
            'csschds' : ['CSSCHD',],
            #'SPLINE3', 'SPLINE6', 'SPLINE7',
            'trims' : ['TRIM',],

            # coords
            'coords' : ['CORD1R', 'CORD1C', 'CORD1S',
                        'CORD2R', 'CORD2C', 'CORD2S',
                        'GMCORD'],

            # temperature cards
            'tempds' : ['TEMPD'],

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
            'transfer_functions' : ['TF'],
            'delays' : ['DELAY'],

            'frequencies' : ['FREQ', 'FREQ1', 'FREQ2', 'FREQ4'],

            # direct matrix input cards
            'dmigs' : ['DMIG'],
            'dmijs' : ['DMIJ'],
            'dmijis' : ['DMIJI'],
            'dmiks' : ['DMIK'],
            'dmis' : ['DMI'],

            # optimzation
            'dequations' : ['DEQATN'],
            'dtable' : ['DTABLE'],
            'dconadds' : ['DCONADD'],
            'dconstrs' : ['DCONSTR'],
            'desvars' : ['DESVAR'],
            'ddvals' : ['DDVAL'],
            'dlinks' : ['DLINK'],
            'dresps' : ['DRESP1', 'DRESP2', 'DRESP3',],
            'dvprels' : ['DVPREL1', 'DVPREL2'],
            'dvmrels' : ['DVMREL1', 'DVMREL2'],
            'dvcrels' : ['DVCREL1', 'DVCREL2'],
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
                'TABLEHT', 'TABRNDG',
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

        Parameters
        ----------
        self : BDF()
            the BDF object
        nparse_errors : int
            how many parse errors should be stored
            (default=0; all=None; no storage=0)
        stop_on_parsing_error : bool
            should an error be raised if there
            are parsing errors (default=True)
        nxref_errors : int
            how many cross-reference errors
            should be stored (default=0; all=None; no storage=0)
        stop_on_xref_error : bool
            should an error be raised if there
            are cross-reference errors (default=True)
        """
        self._nparse_errors = nparse_errors
        self._nxref_errors = nxref_errors
        self._stop_on_parsing_error = stop_on_parsing_error
        self._stop_on_xref_error = stop_on_xref_error

    def read_bdf(self, bdf_filename=None,
                 xref=True, punch=False, encoding=None):
        """
        Read method for the bdf files

        Parameters
        ----------
        self : BDF()
            the BDF object
        bdf_filename : str / None
            the input bdf (default=None; popup a dialog)
        xref :  bool
            should the bdf be cross referenced (default=True)
        punch : bool
            indicates whether the file is a punch file (default=False)
        encoding : str
            the unicode encoding (default=None; system default)

        .. code-block:: python

            >>> bdf = BDF()
            >>> bdf.read_bdf(bdf_filename, xref=True)
            >>> g1 = bdf.Node(1)
            >>> print(g1.get_position())
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

        if not os.path.exists(bdf_filename):
            msg = 'cannot find bdf_filename=%r\n%s' % (bdf_filename, print_bad_path(bdf_filename))
            raise IOError(msg)

        #: the active filename (string)
        self.bdf_filename = bdf_filename

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

            sol, method, isol_line = parse_executive_control_deck(executive_control_lines)
            self.update_solution(sol, method, isol_line)

            self.case_control_deck = CaseControlDeck(self.case_control_lines, self.log)
            self.case_control_deck.solmap_toValue = self._solmap_to_value
            self.case_control_deck.rsolmap_toStr = self.rsolmap_toStr

            cards, card_count = self.get_bdf_cards(bulk_data_lines)
            self._parse_cards(cards, card_count)
            self.pop_parse_errors()
            self.fill_dmigs()

            self.cross_reference(xref=xref)
            self._xref = xref
        except:
            raise
        self.log.debug('---finished BDF.read_bdf of %s---' % self.bdf_filename)
        self.pop_xref_errors()

    def fill_dmigs(self):
        """fills the DMIx cards with the column data that's been stored"""
        for name, card_comments in iteritems(self._dmig_temp):
            card0, comment0 = card_comments[0]
            card_name = card0[0]
            card_name = card_name.rstrip(' *').upper()

            if card_name == 'DMIG':
                # if field2 == 'UACCEL':  # special DMIG card
                card = self.dmigs[name]
            elif card_name == 'DMI':
                card = self.dmis[name]
            elif card_name == 'DMIJ':
                card = self.dmijs[name]
            elif card_name == 'DMIJI':
                card = self.dmijis[name]
            elif card_name == 'DMIK':
                card = self.dmiks[name]
            else:
                raise NotImplementedError(card_name)

            for (card_obj, comment) in card_comments:
                card._add_column(card_obj, comment=comment)

        self._dmig_temp = defaultdict(list)

    def pop_parse_errors(self):
        """raises an error if there are parsing errors"""
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
        """raises an error if there are cross-reference errors"""
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


    def update_solution(self, sol, method, isol_line):
        """
        Updates the overall solution type (e.g. 101,200,600)

        Parameters
        ----------
        self : BDF
            the object pointer
        sol : int
            the solution type (101,103, etc)
        method : str
            the solution method (only for SOL=600)
        isol_line : int
            the line to put the SOL/method on
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

        Parameters
        ----------
        self : BDF()
            the BDF object
        dict_of_vars : dict[str] = int/float/str
            dictionary of 7 character variable names to map.

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

        Parameters
        ----------
        self : BDF()
            the BDF object
        card_name : str
            the card_name -> 'GRID'
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

        Parameters
        ----------
        self : BDF()
            the BDF object
        card_lines : List[str]
            list of strings that represent the card's lines

        Returns
        -------
        fields : list[str]
            the parsed card's fields
        card_name : str
            the card's name

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

    def create_card_object(self, card_lines, card_name, is_list=True, has_none=True):
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
                fields = [print_field_16(self._parse_dynamic_syntax(field)) if '%' in
                          field.strip()[0:1] else print_field_16(field) for field in fields]
                has_none = False

            if has_none:
                card = wipe_empty_fields([print_field_16(field) for field in fields])
            else:
                #card = remove_trailing_fields(fields)
                card = wipe_empty_fields(fields)
            card_obj = BDFCard(card, has_none=False)
        return card_obj, card

    def _make_card_parser(self):
        """creates the card parser variables that are used by add_card"""
        self._card_parser = {
            'GRID' : (GRID, self.add_node),

            'FORCE' : (FORCE, self.add_load),
            'FORCE1' : (FORCE1, self.add_load),
            'FORCE2' : (FORCE2, self.add_load),
            'MOMENT' : (MOMENT, self.add_load),
            'MOMENT1' : (MOMENT1, self.add_load),
            'MOMENT2' : (MOMENT2, self.add_load),

            'GRAV' : (GRAV, self.add_load),
            'ACCEL' : (ACCEL, self.add_load),
            'ACCEL1' : (ACCEL1, self.add_load),
            'LOAD' : (LOAD, self.add_load),
            'PLOAD' : (PLOAD, self.add_load),
            'PLOAD1' : (PLOAD1, self.add_load),
            'PLOAD2' : (PLOAD2, self.add_load),
            'PLOAD4' : (PLOAD4, self.add_load),
            'PLOADX1' : (PLOADX1, self.add_load),
            'RFORCE' : (RFORCE, self.add_load),
            'SLOAD' : (SLOAD, self.add_load),
            'RANDPS' : (RANDPS, self.add_load),
            'GMLOAD' : (GMLOAD, self.add_load),
            'SPCD' : (SPCD, self.add_load),  # enforced displacement
            'QVOL' : (QVOL, self.add_load),  # thermal


            'PLOTEL' : (PLOTEL, self.add_plotel),
            'CQUAD' : (CQUAD, self.add_element),
            'CQUAD4' : (CQUAD4, self.add_element),
            'CQUAD8' : (CQUAD8, self.add_element),
            'CQUADX' : (CQUADX, self.add_element),
            'CQUADR' : (CQUADR, self.add_element),
            'CTRIA3' : (CTRIA3, self.add_element),
            'CTRIA6' : (CTRIA6, self.add_element),
            'CTRIAR' : (CTRIAR, self.add_element),
            'CTRIAX' : (CTRIAX, self.add_element),
            'CTRIAX6' : (CTRIAX6, self.add_element),
            'PCOMP' : (PCOMP, self.add_property),
            'PCOMPG' : (PCOMPG, self.add_property),
            'PLPLANE' : (PLPLANE, self.add_property),
            'PSHELL' : (PSHELL, self.add_property),

            'CSHEAR' : (CSHEAR, self.add_element),
            'PSHEAR' : (PSHEAR, self.add_property),

            # CTETRA - added later
            # CHEXA  - added later
            # CPENTA - added later
            # CPYRAM - added later
            'PSOLID' : (PSOLID, self.add_property),
            'PLSOLID' : (PLSOLID, self.add_property),
            'PCONEAX' : (PCONEAX, self.add_property),

            'CBAR' : (CBAR, self.add_element),
            'PBAR' : (PBAR, self.add_property),
            'PBARL' : (PBARL, self.add_property),

            'CBEAM' : (CBEAM, self.add_element),
            'PBEAM' : (PBEAM, self.add_property),
            'PBEAML' : (PBEAML, self.add_property),
            'PBCOMP' : (PBCOMP, self.add_property),

            'CBEAM3' : (CBEAM3, self.add_element),
            #'PBEAM3' : (PBEAM3, self.add_property),

            'CBEND' : (CBEND, self.add_element),
            #'PBEND' : (PBEND, self.add_property),

            'CBUSH' : (CBUSH, self.add_damper),
            'CBUSH1D' : (CBUSH1D, self.add_damper),
            'CBUSH2D' : (CBUSH2D, self.add_damper),
            'PBUSH' : (PBUSH, self.add_property),
            'PBUSH1D' : (PBUSH1D, self.add_property),

            'CELAS1' : (CELAS1, self.add_element),
            'CELAS2' : (CELAS2, self.add_element),
            'CELAS3' : (CELAS3, self.add_element),
            'CELAS4' : (CELAS4, self.add_element),
            'CVISC' : (CVISC, self.add_element),
            'PELAST' : (PELAST, self.add_PELAST),

            'CDAMP1' : (CDAMP1, self.add_damper),
            'CDAMP2' : (CDAMP2, self.add_damper),
            'CDAMP3' : (CDAMP3, self.add_damper),
            # CDAMP4 added later because the documentation is wrong
            'CDAMP5' : (CDAMP5, self.add_damper),
            'PDAMP5' : (PDAMP5, self.add_property),

            'CONROD' : (CONROD, self.add_element),

            'CROD' : (CROD, self.add_element),
            'PROD' : (PROD, self.add_property),

            'CTUBE' : (CTUBE, self.add_element),
            'PTUBE' : (PTUBE, self.add_property),

            'CFAST' : (CFAST, self.add_damper),
            'PFAST' : (PFAST, self.add_property),

            'CGAP' : (CGAP, self.add_element),
            'PGAP' : (PGAP, self.add_property),

            'CRAC2D' : (CRAC2D, self.add_element),
            'PRAC2D' : (PRAC2D, self.add_property),

            'CRAC3D' : (CRAC3D, self.add_element),
            'PRAC3D' : (PRAC3D, self.add_property),

            'PDAMPT' : (PDAMPT, self.add_PDAMPT),
            'PBUSHT' : (PBUSHT, self.add_PBUSHT),

            'LSEQ' : (LSEQ, self.add_LSEQ),
            'PHBDY' : (PHBDY, self.add_PHBDY),
            'AERO' : (AERO, self.add_AERO),
            'AEROS' : (AEROS, self.add_AEROS),
            'AECOMP' : (AECOMP, self.add_AECOMP),
            'AEFACT' : (AEFACT, self.add_AEFACT),
            'AELINK' : (AELINK, self.add_AELINK),
            'AELIST' : (AELIST, self.add_AELIST),
            'AEPARM' : (AEPARM, self.add_AEPARM),
            'AESTAT' : (AESTAT, self.add_AESTAT),
            'AESURF' : (AESURF, self.add_AESURF),
            'AESURFS' : (AESURFS, self.add_AESURF),

            'TRIM' : (TRIM, self.add_TRIM),
            'FLUTTER' : (FLUTTER, self.add_FLUTTER),
            'FLFACT' : (FLFACT, self.add_FLFACT),
            'GUST' : (GUST, self.add_GUST),
            'CSSCHD' : (CSSCHD, self.add_CSSCHD),
            'NLPARM' : (NLPARM, self.add_NLPARM),
            'NLPCI' : (NLPCI, self.add_NLPCI),
            'TSTEP' : (TSTEP, self.add_TSTEP),
            'TSTEPNL' : (TSTEPNL, self.add_TSTEPNL),

            'SESET' : (SESET, self.add_SESET),
            'DCONSTR' : (DCONSTR, self.add_DCONSTR),
            'DESVAR' : (DESVAR, self.add_DESVAR),
            'DDVAL' : (DDVAL, self.add_DDVAL),
            'DLINK' : (DLINK, self.add_DLINK),
            'PARAM' : (PARAM, self.add_PARAM),

            'TF' : (TF, self.add_TF),
            'DELAY' : (DELAY, self.add_DELAY),
            'DCONADD' : (DCONADD, self.add_DCONADD),

            'RBAR' : (RBAR, self.add_rigid_element),
            'RBAR1' : (RBAR1, self.add_rigid_element),
            'RBE1' : (RBE1, self.add_rigid_element),
            'RBE2' : (RBE2, self.add_rigid_element),
            'RBE3' : (RBE3, self.add_rigid_element),

            # there is no MAT6 or MAT7
            'MAT1' : (MAT1, self.add_structural_material),
            'MAT2' : (MAT2, self.add_structural_material),
            'MAT3' : (MAT3, self.add_structural_material),
            'MAT8' : (MAT8, self.add_structural_material),
            'MAT9' : (MAT9, self.add_structural_material),
            'MAT10' : (MAT10, self.add_structural_material),
            'MAT11' : (MAT11, self.add_structural_material),
            'EQUIV' : (EQUIV, self.add_structural_material),

            #'MATHE' : (MATHE, self.add_hyperelastic_material),
            'MATHP' : (MATHP, self.add_hyperelastic_material),
            'MAT4' : (MAT4, self.add_thermal_material),
            'MAT5' : (MAT5, self.add_thermal_material),

            'MATS1' : (MATS1, self.add_material_dependence),
            #'MATS3' : (MATS3, self.add_material_dependence),
            #'MATS8' : (MATS8, self.add_material_dependence),
            'MATT1' : (MATT1, self.add_material_dependence),
            'MATT2' : (MATT2, self.add_material_dependence),
            #'MATT3' : (MATT3, self.add_material_dependence),
            'MATT4' : (MATT4, self.add_material_dependence),
            'MATT5' : (MATT5, self.add_material_dependence),
            #'MATT8' : (MATT8, self.add_material_dependence),
            #'MATT9' : (MATT9, self.add_material_dependence),

            # hasnt been verified, links up to MAT1, MAT2, MAT9 w/ same MID
            'CREEP' : (CREEP, self.add_creep_material),


            'DLOAD' : (DLOAD, self.add_dload),
            'TLOAD1' : (TLOAD1, self.add_dload_entry),
            'TLOAD2' : (TLOAD2, self.add_dload_entry),
            'RLOAD1' : (RLOAD1, self.add_dload_entry),
            'RLOAD2' : (RLOAD2, self.add_dload_entry),

            'TEMP' : (TEMP, self.add_thermal_load),
            'QBDY1' : (QBDY1, self.add_thermal_load),
            'QBDY2' : (QBDY2, self.add_thermal_load),
            'QBDY3' : (QBDY3, self.add_thermal_load),
            'QHBDY' : (QHBDY, self.add_thermal_load),

            'CHBDYE' : (CHBDYE, self.add_thermal_element),
            'CHBDYG' : (CHBDYG, self.add_thermal_element),
            'CHBDYP' : (CHBDYP, self.add_thermal_element),

            'PCONV' : (PCONV, self.add_convection_property),
            'PCONVM' : (PCONVM, self.add_convection_property),

            'MPC' : (MPC, self.add_constraint_MPC),
            'MPCADD' : (MPCADD, self.add_constraint_MPC),

            'SPC' : (SPC, self.add_constraint_SPC),
            'SPC1' : (SPC1, self.add_constraint_SPC),
            'SPCAX' : (SPCAX, self.add_constraint_SPC),
            'SPCADD' : (SPCADD, self.add_constraint_SPC),
            'GMSPC' : (GMSPC, self.add_constraint_SPC),

            'SUPORT' : (SUPORT, self.add_suport), # pseudo-constraint

            'CAERO1' : (CAERO1, self.add_CAERO),
            'CAERO2' : (CAERO2, self.add_CAERO),
            'CAERO3' : (CAERO3, self.add_CAERO),
            'CAERO4' : (CAERO4, self.add_CAERO),
            'CAERO5' : (CAERO5, self.add_CAERO),

            'PAERO1' : (PAERO1, self.add_PAERO),
            'PAERO2' : (PAERO2, self.add_PAERO),
            'PAERO3' : (PAERO3, self.add_PAERO),
            #'PAERO4' : (PAERO4, self.add_PAERO),
            'PAERO5' : (PAERO5, self.add_PAERO),

            'SPLINE1' : (SPLINE1, self.add_SPLINE),
            'SPLINE2' : (SPLINE2, self.add_SPLINE),
            'SPLINE3' : (SPLINE3, self.add_SPLINE),
            'SPLINE4' : (SPLINE4, self.add_SPLINE),
            'SPLINE5' : (SPLINE5, self.add_SPLINE),

            'MONPNT1' : (MONPNT1, self.add_MONPNT),
            'MKAERO1' : (MKAERO1, self.add_MKAERO),
            'MKAERO2' : (MKAERO2, self.add_MKAERO),

            'SUPORT1' : (SUPORT1, self.add_suport1),  # pseudo-constraint
            #'SESUP' : (SESUP, self.add_SESUP),  # pseudo-constraint

            'FREQ' : (FREQ, self.add_FREQ),
            'FREQ1' : (FREQ1, self.add_FREQ),
            'FREQ2' : (FREQ2, self.add_FREQ),
            'FREQ4' : (FREQ4, self.add_FREQ),

            'ASET' : (ASET, self.add_ASET),
            'ASET1' : (ASET1, self.add_ASET),

            'BSET' : (BSET, self.add_BSET),
            'BSET1' : (BSET1, self.add_BSET),

            'CSET' : (CSET, self.add_CSET),
            'CSET1' : (CSET1, self.add_CSET),

            'QSET' : (QSET, self.add_QSET),
            'QSET1' : (QSET1, self.add_QSET),

            'USET' : (USET, self.add_USET),
            'USET1' : (USET1, self.add_USET),

            'SET1' : (SET1, self.add_SET),
            'SET3' : (SET3, self.add_SET),

            'SEBSET' : (SEBSET, self.add_SEBSET),
            'SEBSET1' : (SEBSET1, self.add_SEBSET),

            'SECSET' : (SECSET, self.add_SECSET),
            'SECSET1' : (SECSET1, self.add_SECSET),

            'SEQSET' : (SEQSET, self.add_SEQSET),
            'SEQSET1' : (SEQSET1, self.add_SEQSET),

            #'SEUSET' : (SEUSET, self.add_SEUSET),
            #'SEUSET1' : (SEUSET1, self.add_SEUSET),

            'DTABLE' : (DTABLE, self.add_DTABLE),

            'DRESP1' : (DRESP1, self.add_DRESP),
            'DRESP2' : (DRESP2, self.add_DRESP), # deqatn
            'DRESP3' : (DRESP3, self.add_DRESP),

            'DVPREL1' : (DVPREL1, self.add_DVPREL),
            'DVPREL2' : (DVPREL2, self.add_DVPREL), # deqatn

            'DVMREL1' : (DVMREL1, self.add_DVMREL),
            #'DVMREL2' : (DVMREL2, self.add_DVMREL), # deqatn
            #DVCREL1
            # DVCREL2 - deqatn

            'CORD2R' : (CORD2R, self.add_coord),
            'CORD2C' : (CORD2C, self.add_coord),
            'CORD2S' : (CORD2S, self.add_coord),
            'GMCORD' : (GMCORD, self.add_coord),

            'TABLED1' : (TABLED1, self.add_table),
            'TABLED2' : (TABLED2, self.add_table),
            'TABLED3' : (TABLED3, self.add_table),
            'TABLED4' : (TABLED4, self.add_table),
            'TABLEM1' : (TABLEM1, self.add_table),
            'TABLEM2' : (TABLEM2, self.add_table),
            'TABLEM3' : (TABLEM3, self.add_table),
            'TABLEM4' : (TABLEM4, self.add_table),

            'TABLES1' : (TABLES1, self.add_table),
            'TABLEST' : (TABLEST, self.add_table),

            'TABDMP1' : (TABDMP1, self.add_table_sdamping),
            'TABRND1' : (TABRND1, self.add_random_table),
            'TABRNDG' : (TABRNDG, self.add_random_table),

            'EIGB' : (EIGB, self.add_method),
            'EIGR' : (EIGR, self.add_method),
            'EIGRL' : (EIGRL, self.add_method),
            'EIGC' : (EIGC, self.add_cmethod),
            'EIGP' : (EIGP, self.add_cmethod),

            'BCRPARA' : (BCRPARA, self.add_BCRPARA),
            'BCTADD' : (BCTADD, self.add_BCTADD),
            'BCTPARA' : (BCTPARA, self.add_BCTPARA),
            'BSURF' : (BSURF, self.add_BSURF),
            'BSURFS' : (BSURFS, self.add_BSURFS),

            'CONM1' : (CONM1, self.add_mass),
            'CONM2' : (CONM2, self.add_mass),
            'CMASS1' : (CMASS1, self.add_mass),
            'CMASS2' : (CMASS2, self.add_mass),
            'CMASS3' : (CMASS3, self.add_mass),
            # CMASS4 - added later because documentation is wrong

            'DOPTPRM' : (DOPTPRM, self._add_doptprm),
            'SPOINT' : (SPOINTs, self.add_SPOINT),
            'EPOINT' : (EPOINTs, self.add_EPOINT),

    #elif card_name == 'BCTSET':
        #card = BCTSET(card_obj, comment=comment, sol=self.sol)
        #self.add_BCTSET(card)
        }
        self._card_parser_b = {
            'CTETRA' : self._prepare_ctetra,
            'CPYRAM' : self._prepare_cpyram,
            'CPENTA' : self._prepare_cpenta,
            'CHEXA' : self._prepare_chexa,

            'CORD1R' : self._prepare_cord1r,
            'CORD1C' : self._prepare_cord1c,
            'CORD1S' : self._prepare_cord1s,
            #'CORD3G' : self._prepare_CORD3G,

            'DAREA' : self._prepare_darea,
            'PMASS' : self._prepare_pmass,
            'CMASS4' : self._prepare_cmass4,
            'CDAMP4' : self._prepare_cdamp4,

            'DMIG' : self._prepare_dmig,
            'DMI' : self._prepare_dmi,
            'DMIJ' : self._prepare_dmij,
            'DMIK' : self._prepare_dmik,
            'DMIJI' : self._prepare_dmiji,

            'DEQATN' : self._prepare_dequatn,

            'PVISC' : self._prepare_pvisc,
            'PELAS' : self._prepare_pelas,
            'PDAMP' : self._prepare_pdamp,

            'TEMPD' : self._prepare_tempd,
            'CONV' : self._prepare_conv,
            'RADM' : self._prepare_radm,
            'RADBC' : self._prepare_radbc,
            'GRDSET' : self._prepare_grdset,

            'BCTSET' : self._prepare_bctset,
        }

    def _prepare_bctset(self, card, card_obj, comment=''):
        card = BCTSET(card_obj, comment=comment, sol=self.sol)
        self.add_BCTSET(card)

    def _prepare_grdset(self, card, card_obj, comment=''):
        self.gridSet = GRDSET(card_obj, comment=comment)

    def _prepare_cdamp4(self, card, card_obj, comment=''):
        self.add_damper(CDAMP4(card_obj, comment=comment))
        if card_obj.field(5):
            self.add_damper(CDAMP4(card_obj, 1, comment=''))
        return card_obj

    def _prepare_conv(self, card, card_obj, comment=''):
        bc = CONV(card_obj, comment=comment)
        self.add_thermal_BC(bc, bc.eid)

    def _prepare_radm(self, card, card_obj, comment=''):
        bc = RADM(card_obj, comment=comment)
        self.add_thermal_BC(bc, bc.radmid)

    def _prepare_radbc(self, card, card_obj, comment=''):
        bc = RADBC(card_obj, comment=comment)
        self.add_thermal_BC(bc, bc.nodamb)

    def _prepare_tempd(self, card, card_obj, comment=''):
        self.add_TEMPD(TEMPD(card_obj, 0, comment=''))
        if card_obj.field(3):
            self.add_TEMPD(TEMPD(card_obj, 1, comment=''))
            if card_obj.field(5):
                self.add_TEMPD(TEMPD(card_obj, 1, comment=''))
                if card_obj.field(7):
                    self.add_TEMPD(TEMPD(card_obj, 1, comment=''))

    def _add_doptprm(self, doptprm, comment=''):
        self.doptprm = doptprm

    def _prepare_dequatn(self, card, card_obj, comment=''):
        if hasattr(self, 'test_deqatn') or 1:
            self.add_DEQATN(DEQATN(card_obj, comment=comment))
        else:
            if comment:
                self.rejects.append([comment])
            self.rejects.append(card)

    def _prepare_dmig(self, card, card_obj, comment=''):
        """adds a DMIG"""
        # not done...
        field2 = integer_or_string(card_obj, 2, 'flag')
        if field2 == 0:
            card = DMIG(card_obj, comment=comment)
            self.add_DMIG(card)
        # elif field2 == 'UACCEL':  # special DMIG card
            # self.reject_cards(card)
        else:
            name = string(card_obj, 1, 'name')
            self._dmig_temp[name].append((card_obj, comment))


    def _prepare_dmix(self, class_obj, add_method, card_obj, comment=''):
        """adds a DMIx"""
        #elif card_name in ['DMI', 'DMIJ', 'DMIJI', 'DMIK']:
        field2 = integer(card_obj, 2, 'flag')
        if field2 == 0:
            add_method(class_obj(card_obj, comment=comment))
        else:
            name = string(card_obj, 1, 'name')
            self._dmig_temp[name].append((card_obj, comment))

    def _prepare_dmi(self, card, card_obj, comment=''):
        """adds a DMI"""
        self._prepare_dmix(DMI, self.add_DMI, card_obj, comment=comment)

    def _prepare_dmij(self, card, card_obj, comment=''):
        """adds a DMIJ"""
        self._prepare_dmix(DMIJ, self.add_DMIJ, card_obj, comment=comment)

    def _prepare_dmik(self, card, card_obj, comment=''):
        """adds a DMIK"""
        self._prepare_dmix(DMIK, self.add_DMIK, card_obj, comment=comment)

    def _prepare_dmiji(self, card, card_obj, comment=''):
        """adds a DMIJI"""
        self._prepare_dmix(DMIJI, self.add_DMIJI, card_obj, comment=comment)

    def _prepare_cmass4(self, card, card_obj, comment=''):
        """adds a CMASS4"""
        class_instance = CMASS4(card_obj, icard=0, comment=comment)
        self.add_mass(class_instance)
        if card_obj.field(5):
            class_instance = CMASS4(card_obj, icard=1, comment=comment)
            self.add_mass(class_instance)

    def _prepare_pelas(self, card, card_obj, comment=''):
        """adds a PELAS"""
        class_instance = PELAS(card_obj, icard=0, comment=comment)
        self.add_property(class_instance)
        if card_obj.field(5):
            class_instance = PELAS(card_obj, icard=1, comment=comment)
            self.add_property(class_instance)

    def _prepare_pvisc(self, card, card_obj, comment=''):
        """adds a PVISC"""
        class_instance = PVISC(card_obj, icard=0, comment=comment)
        self.add_property(class_instance)
        if card_obj.field(5):
            class_instance = PVISC(card_obj, icard=1, comment=comment)
            self.add_property(class_instance)

    def _prepare_pdamp(self, card, card_obj, comment=''):
        """adds a PDAMP"""
        class_instance = PDAMP(card_obj, icard=0, comment=comment)
        self.add_property(class_instance)
        if card_obj.field(3):
            class_instance = PDAMP(card_obj, icard=1, comment=comment)
            self.add_property(class_instance)
        if card_obj.field(5):
            class_instance = PDAMP(card_obj, icard=2, comment=comment)
            self.add_property(class_instance)
        if card_obj.field(7):
            class_instance = PDAMP(card_obj, icard=3, comment=comment)
            self.add_property(class_instance)

    def _prepare_pmass(self, card, card_obj, comment=''):
        """adds a PMASS"""
        card_instance = PMASS(card_obj, icard=0, comment=comment)
        self.add_property_mass(card_instance)
        for (i, j) in enumerate([3, 5, 7]):
            if card_obj.field(j):
                card_instance = PMASS(card_obj, icard=i+1, comment=comment)
                self.add_property_mass(card_instance)

    def _prepare_darea(self, card, card_obj, comment=''):
        """adds a DAREA"""
        class_instance = DAREA(card_obj, comment=comment)
        self.add_DAREA(class_instance)
        if card_obj.field(5):
            class_instance = DAREA(card_obj, icard=1, comment=comment)
            self.add_DAREA(class_instance)

    def _prepare_cord1r(self, card, card_obj, comment=''):
        """adds a CORD1R"""
        class_instance = CORD1R(card_obj, comment=comment)
        self.add_coord(class_instance)
        if card_obj.field(5):
            class_instance = CORD1R(card_obj, icard=1, comment=comment)
            self.add_coord(class_instance)

    def _prepare_cord1c(self, card, card_obj, comment=''):
        """adds a CORD1C"""
        class_instance = CORD1C(card_obj, comment=comment)
        self.add_coord(class_instance)
        if card_obj.field(5):
            class_instance = CORD1C(card_obj, icard=1, comment=comment)
            self.add_coord(class_instance)

    def _prepare_cord1s(self, card, card_obj, comment=''):
        """adds a CORD1S"""
        class_instance = CORD1S(card_obj, comment=comment)
        self.add_coord(class_instance)
        if card_obj.field(5):
            class_instance = CORD1S(card_obj, icard=1, comment=comment)
            self.add_coord(class_instance)

    def _prepare_ctetra(self, card, card_obj, comment=''):
        """adds a CTETRA"""
        card_class = CTETRA4 if card_obj.nfields == 7 else CTETRA10
        class_instance = card_class(card_obj, comment=comment)
        self.add_element(class_instance)

    def _prepare_cpyram(self, card, card_obj, comment=''):
        """adds a CPYRAM"""
        card_class = CPYRAM5 if card_obj.nfields == 8 else CPYRAM13
        class_instance = card_class(card_obj, comment=comment)
        self.add_element(class_instance)

    def _prepare_cpenta(self, card, card_obj, comment=''):
        """adds a CPENTA"""
        card_class = CPENTA6 if card_obj.nfields == 9 else CPENTA15
        class_instance = card_class(card_obj, comment=comment)
        self.add_element(class_instance)

    def _prepare_chexa(self, card, card_obj, comment=''):
        """adds a CHEXA"""
        card_class = CHEXA8 if card_obj.nfields == 11 else CHEXA20
        class_instance = card_class(card_obj, comment=comment)
        self.add_element(class_instance)

    def add_card(self, card_lines, card_name, comment='', is_list=True, has_none=True):
        """
        Adds a card object to the BDF object.

        Parameters
        ----------
        self : BDF
            the BDF object
        card_lines: list[str]
            the list of the card fields
        card_name : str
            the card_name -> 'GRID'
        comment : str
            an optional the comment for the card
        is_list : bool, optional
            False : input is a list of card fields -> ['GRID', 1, 2, 3.0, 4.0, 5.0]
            True :  input is list of card_lines -> ['GRID, 1, 2, 3.0, 4.0, 5.0']

        Returns
        -------
        card_object : BDFCard()
            the card object representation of card

        .. code-block:: python

          >>> model = BDF()

          # is_list is a somewhat misleading name; is it a list of card_lines
          # where a card_line is an unparsed string
          >>> card_lines = ['GRID,1,2']
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
        card_obj, card = self.create_card_object(card_lines, card_name, is_list=is_list, has_none=has_none)

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
        failed = True
        try:
            card_class, add_card_function = self._card_parser[card_name]
            failed = False
        except KeyError:
            try:
                add_card_function = self._card_parser_b[card_name]
            except KeyError:
                #: .. warning:: cards with = signs in them
                #:              are not announced when they are rejected
                #if '=' not in card[0]:
                #    self.log.info('rejecting processed equal signed card %s' % card)
                #self.reject_cards.append(card)
                msg = 'card_name=%r not in card_parser_b'  % (card_name)
                raise KeyError(msg)
                #print(msg)
                #return

            try:
                add_card_function(card, card_obj, comment=comment)
            except (SyntaxError, AssertionError, KeyError, ValueError) as e:
                # WARNING: Don't catch RuntimeErrors or a massive memory leak can occur
                #tpl/cc451.bdf
                raise
                # NameErrors should be caught
                self._iparse_errors += 1
                var = traceback.format_exception_only(type(e), e)
                self._stored_parse_errors.append((card, var))
                if self._iparse_errors > self._nparse_errors:
                    self.pop_parse_errors()
        else:
            assert failed == False, failed
            # KeyError won't trigger this

            try:
                class_instance = card_class(card_obj, comment=comment)
                add_card_function(class_instance)
            except (SyntaxError, AssertionError, ValueError) as e:
                # WARNING: Don't catch RuntimeErrors or a massive memory leak can occur
                #tpl/cc451.bdf
                #raise
                # NameErrors should be caught
                self._iparse_errors += 1
                var = traceback.format_exception_only(type(e), e)
                self._stored_parse_errors.append((card, var))
                if self._iparse_errors > self._nparse_errors:
                    self.pop_parse_errors()
        return card_obj

    def get_bdf_stats(self, return_type='string'):
        """
        Print statistics for the BDF

        Parameters
        ----------
        self : BDF()
            the BDF object

        Returns
        -------
        return_type : str, optional
            the output type ('list', 'string')
                'list' : list of strings
                'string' : single, joined string

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
            'dconadds', 'dconstrs', 'desvars', 'ddvals', 'dlinks', 'dresps',
            'dvcrels', 'dvmrels', 'dvprels',

            # SESETx - dict
            'suport1',
            'se_sets',
            'se_usets',

            # tables
            'tables', 'randomTables',

            # methods
            'methods', 'cMethods',

            # aero
            'caeros', 'paeros', 'aero', 'aeros', 'aecomps', 'aefacts', 'aelinks',
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
            'case_control_deck', 'caseControlDeck',
            'spcObject2', 'mpcObject2',

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

        unsupported_types = ignored_types.union(ignored_types2)
        all_params = object_attributes(self, keys_to_skip=unsupported_types)

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

    def get_displcement_index_transforms(self):
        """
        Get index and transformation matricies for nodes with
        their output in coordinate systems other than the global.
        Used in combination with ``OP2.transform_displacements_to_global``

        Parameters
        ----------
        self : BDF
            BDF object.

        Returns
        ----------
        i_transform : dict{float:ndarray}
            Dictionary from coordinate id to index of the nodes in
            ``self.point_ids`` that their output (`CD`) in that
            coordinate system.
        transforms : dict{float:ndarray}
            Dictionary from coordinate id to 3 x 3 transformation
            matrix for that coordinate system.
        """
        nids_transform = {}
        i_transform = {}
        transforms = {}
        if len(self.coords) < 2:
            return i_transform, transforms
        for nid, node in sorted(iteritems(self.nodes)):
            cid_d = node.Cd()
            if cid_d:
                if cid_d not in nids_transform:
                    nids_transform[cid_d] = []
                nids_transform[cid_d].append(nid)

        nids_all = np.array(sorted(self.point_ids))
        for cid in sorted(nids_transform.keys()):
            nids = np.array(nids_transform[cid])
            i_transform[cid] = np.where(np.in1d(nids_all, nids))[0]
            transforms[cid] = self.coords[cid].beta()
        return i_transform, transforms

    def _get_card_name(self, lines):
        """
        Returns the name of the card defined by the provided lines

        Parameters
        ----------
        self : BDF()
            the BDF object
        lines : list[str]
            the lines of the card

        Returns
        -------
        card_name : str
            the name of the card
        """
        card_name = lines[0][:8].rstrip('\t, ').split(',')[0].split('\t')[0].strip('*\t ')
        if len(card_name) == 0:
            return None
        if ' ' in card_name or len(card_name) == 0:
            msg = 'card_name=%r\nline=%r in filename=%r is invalid' \
                  % (card_name, lines[0], self.active_filename)
            print(msg)
            raise CardParseSyntaxError(msg)
        return card_name.upper()

    def _show_bad_file(self, bdf_filename):
        """
        Prints the 10 lines before the UnicodeDecodeError occurred.

        Paramters
        ---------
        bdf_filename : str
            the filename to print the lines of
        """
        lines = []
        print('ENCODING - show_bad_file=%r' % self._encoding)
        with codec_open(bdf_filename, 'rU', encoding=self._encoding) as bdf_file:
            iline = 0
            nblank = 0
            while 1:
                try:
                    line = bdf_file.readline().rstrip()
                except UnicodeDecodeError as e:
                    i0 = max([iline - 10, 0])
                    self.log.error('filename=%s' % self.bdf_filename)
                    for i1, line in enumerate(lines[i0:iline]):
                        self.log.error('lines[%i]=%r' % (i0 + i1, line))
                    msg = "\n%s encoding error on line=%s of %s; not '%s'" % (self._encoding, iline, bdf_filename, self._encoding)
                    raise RuntimeError(msg)
                if line:
                    nblank = 0
                else:
                    nblank += 1
                if nblank == 20:
                    raise RuntimeError('20 blank lines')
                iline += 1
                lines.append(line)

    def _get_lines(self, bdf_filename, punch=False):
        """
        Opens the bdf and extracts the lines

        Parameters
        ----------

        bdf_filename : str
            the main bdf_filename
        punch : bool, optional
            is this a punch file (default=False; no executive/case control decks)

        Returns
        -------
        executive_control_lines : list[str]
            the executive control deck as a list of strings
        case_control_lines : list[str]
            the case control deck as a list of strings
        bulk_data_lines : list[str]
            the bulk data deck as a list of strings
        """
        #: the directory of the 1st BDF (include BDFs are relative to this one)
        self.include_dir = os.path.dirname(os.path.abspath(bdf_filename))

        with self._open_file(bdf_filename, basename=True) as bdf_file:
            try:
                lines = bdf_file.readlines()
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

                with self._open_file(bdf_filename2, basename=False) as bdf_file:
                    #print('bdf_file.name = %s' % bdf_file.name)
                    lines2 = bdf_file.readlines()

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

        Parameters
        ----------
        self : BDF()
            the BDF object
        card_name : str
            the card_name -> 'GRID'
        n : int, optional
            the amount to increment by (default=1)

        >>> bdf.read_bdf(bdf_filename)
        >>> bdf.card_count['GRID']
        50
        """
        if card_name in self.card_count:
            self.card_count[card_name] += n
        else:
            self.card_count[card_name] = n

    def _open_file(self, bdf_filename, basename=False):
        """
        Opens a new bdf_filename with the proper encoding and include directory

        Parameters
        ----------
        bdf_filename : str
            the filename to open
        basename : bool (default=False)
            should the basename of bdf_filename be appended to the include directory
        """
        if basename:
            bdf_filename_inc = os.path.join(self.include_dir, os.path.basename(bdf_filename))
        else:
            bdf_filename_inc = os.path.join(self.include_dir, bdf_filename)
        if not os.path.exists(bdf_filename_inc):
            msg = 'No such bdf_filename: %r\n' % bdf_filename_inc
            msg += 'cwd: %r\n' % os.getcwd()
            msg += 'include_dir: %r\n' % self.include_dir
            msg += print_bad_path(bdf_filename_inc)
            raise IOError(msg)
        elif bdf_filename_inc.endswith('.op2'):
            raise IOError('Invalid filetype: bdf_filename=%r' % bdf_filename_inc)
        bdf_filename = bdf_filename_inc

        if bdf_filename in self.active_filenames:
            msg = 'bdf_filename=%s is already active.\nactive_filenames=%s' \
                % (bdf_filename, self.active_filenames)
            raise RuntimeError(msg)
        self.log.debug('opening %r' % bdf_filename)
        self.active_filenames.append(bdf_filename)

        #print('ENCODING - _open_file=%r' % self._encoding)
        bdf_file = codec_open(bdf_filename, 'rU', encoding=self._encoding)
        return bdf_file

    def _parse_cards(self, cards, card_count):
        """creates card objects and adds the parsed cards to the deck"""
        #print('card_count = %s' % card_count)
        for card_name, card in sorted(iteritems(cards)):
            #print('---%r---' % card_name)
            if self.is_reject(card_name):
                self.log.info('    rejecting card_name = %s' % card_name)
                for cardi in card:
                    self._increase_card_count(card_name)
                    self.rejects.append([cardi[0]] + cardi[1])
            else:
                for comment, card_lines in card:
                    self.add_card(card_lines, card_name, comment=comment, is_list=False, has_none=False)

    def _parse_dynamic_syntax(self, key):
        """
        Applies the dynamic syntax for %varName

        Parameters
        ----------
        self : BDF()
            the BDF object
        key : str
            the uppercased key

        Returns
        -------
        value : int/float/str
            the dynamic value defined by dict_of_vars
        .. seealso:: :func: `set_dynamic_syntax`
        """
        key = key.strip()[1:]
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

        Parameters
        ----------
        self : BDF()
            the BDF object
        bdf_filename : str
            the input filename

        ..code-block :: python
            $ pyNastran: version=NX
            $ pyNastran: encoding=latin-1
            $ pyNastran: punch=True

        ..warning :: pyNastran lines must be at the top of the file
        """
        bdf_file = open(bdf_filename, 'r')
        check_header = True
        while check_header:
            try:
                line = bdf_file.readline()
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
                    elif key in ['nnodes', 'nelements']:
                        pass
                    else:
                        raise NotImplementedError(key)
                else:
                    break
            else:
                break
        bdf_file.close()

    def _verify_bdf(self):
        """
        Cross reference verification method.
        """
        try:
            xref = self._xref
        except AttributeError:
            xref = True
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

        for key, card in sorted(iteritems(self.dresps)):
            #try:
            card._verify(xref)
            #except:
                #print(str(card))
                #raise

        for key, card in sorted(iteritems(self.dvcrels)):
            try:
                card._verify(xref)
            except:
                print(str(card))
                raise
        for key, card in sorted(iteritems(self.dvmrels)):
            try:
                card._verify(xref)
            except:
                print(str(card))
                raise
        for key, card in sorted(iteritems(self.dvprels)):
            try:
                card._verify(xref)
            except:
                print(str(card))
                raise

IGNORE_COMMENTS = (
    '$EXECUTIVE CONTROL DECK',
    '$CASE CONTROL DECK',
    'NODES', 'SPOINTS', 'EPOINTS', 'ELEMENTS',
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

    Parameters
    ----------
    comment : str
        the comment to possibly remove

    Returns
    -------
    updated_comment : str
        the comment
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
    """
    shows off how unicode works becausee it's overly complicated
    """
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
    from pyNastran.bdf.test.test_bdf import main
    #main()
