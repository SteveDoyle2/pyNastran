# pylint: disable=W0212,C0103,W0633,W0611,W0201,C0301,R0915,R0912
"""
Main BDF class.  Defines:
  - BDF

see https://docs.plm.automation.siemens.com/tdoc/nxnastran/10/help/#uid:index
"""
from __future__ import (nested_scopes, generators, division, absolute_import,
                        print_function, unicode_literals)
import os
import sys
import traceback
from codecs import open as codec_open
from collections import defaultdict
import io

from six import string_types, iteritems

import numpy as np

import pyNastran
from pyNastran.bdf.field_writer_8 import print_card_8
from pyNastran.bdf.field_writer_16 import print_field_16
from pyNastran.bdf.utils import (
    to_fields, get_include_filename,
    parse_executive_control_deck, clean_empty_lines)
from pyNastran.bdf.cards.methods import EIGB, EIGC, EIGR, EIGP, EIGRL

from pyNastran.utils import _filename, print_bad_path
from pyNastran.utils.dev import list_print, object_attributes
from pyNastran.utils.log import get_logger
from pyNastran.utils.gui_io import load_file_dialog

# coords
from .cards.coord import Coord

# nodes
from .cards.nodes.grid import GRID, GRDSET
from .cards.nodes.point import POINT
from .cards.nodes.spoint import SPOINT
from .cards.nodes.epoint import EPOINT
from .cards.nodes.pointax import POINTAX

###################Elements##################
from pyNastran.bdf.dev_vectorized.cards.elements.elements import Elements
from pyNastran.bdf.dev_vectorized.cards.elements.properties import Properties

# bar
from .cards.elements.bar.cbar import CBAR, CBAROR
#from .cards.elements.bar.pbar import PBAR
#from .cards.elements.bar.pbarl import PBARL
from .cards.elements.bar.properties_bar import PropertiesBar

# beams
from .cards.elements.beam.cbeam import CBEAM #, CBEAMOR
from .cards.elements.beam.properties_beam import PropertiesBeam

# bush
from pyNastran.bdf.dev_vectorized.cards.elements.bush.cbush import CBUSH
from pyNastran.bdf.dev_vectorized.cards.elements.bush.pbush import PBUSH

# mass
from pyNastran.bdf.dev_vectorized.cards.elements.mass.mass import Mass

# rigid
from pyNastran.bdf.dev_vectorized.cards.elements.rigid.rbe2 import RBE2
from pyNastran.bdf.dev_vectorized.cards.elements.rigid.rbe3 import RBE3

# rods
from pyNastran.bdf.dev_vectorized.cards.elements.rod.prod import PROD
from pyNastran.bdf.dev_vectorized.cards.elements.rod.crod import CROD
from pyNastran.bdf.dev_vectorized.cards.elements.rod.conrod import CONROD
from pyNastran.bdf.dev_vectorized.cards.elements.rod.ptube import PTUBE
from pyNastran.bdf.dev_vectorized.cards.elements.rod.ctube import CTUBE

# shear
from pyNastran.bdf.dev_vectorized.cards.elements.shear.cshear import CSHEAR
from pyNastran.bdf.dev_vectorized.cards.elements.shear.pshear import PSHEAR

# shell
from pyNastran.bdf.dev_vectorized.cards.elements.shell.elements_shell import ElementsShell
from pyNastran.bdf.dev_vectorized.cards.elements.shell.properties_shell import PropertiesShell

# solid
from pyNastran.bdf.dev_vectorized.cards.elements.solid.elements_solid import ElementsSolid
from pyNastran.bdf.dev_vectorized.cards.elements.solid.properties_solid import PropertiesSolid

# spring
from pyNastran.bdf.dev_vectorized.cards.elements.spring.elements_spring import ElementsSpring
from pyNastran.bdf.dev_vectorized.cards.elements.spring.pelas import PELAS

#===========================
# aero

from pyNastran.bdf.dev_vectorized.cards.aero.caero import CAero
from pyNastran.bdf.dev_vectorized.cards.aero.paero import PAero
from pyNastran.bdf.dev_vectorized.cards.aero.spline1 import SPLINE1

from pyNastran.bdf.dev_vectorized.cards.aero.trim import TRIM
from pyNastran.bdf.dev_vectorized.cards.aero.aero import AERO
from pyNastran.bdf.dev_vectorized.cards.aero.aeros import AEROS
#from pyNastran.bdf.cards.aero import AEFACT, AELIST, AELINK, AEPARM, AESTAT, AESURF, AESURFS
from pyNastran.bdf.cards.aero import AEFACT, AELIST, AELINK, AEPARM, AESTAT, AESURF

#===========================

# materials
from pyNastran.bdf.dev_vectorized.cards.materials.materials import Materials


# loads
from pyNastran.bdf.dev_vectorized.cards.loads.loads import Loads
from pyNastran.bdf.dev_vectorized.cards.loads.temp import TEMPs
#=============================
# dynamic
#from pyNastran.bdf.dev_vectorized.cards.nonlinear.nlpci import NLPCI
#from pyNastran.bdf.dev_vectorized.cards.nonlinear.nlparm import NLPARM
from pyNastran.bdf.cards.dynamic import NLPARM, NLPCI, TSTEP, TSTEPNL

#=============================

# constraints
from .cards.constraints.spc import SPC, get_spc_constraint
from .cards.constraints.spcd import SPCD

# spc
from .cards.constraints.spc1 import SPC1, get_spc1_constraint
from .cards.constraints.spcadd import SPCADD, get_spcadd_constraint

# mpc
from .cards.constraints.mpc import MPC, get_mpc_constraint
#from .cards.constraints.mpcax import MPCAX
from .cards.constraints.mpcadd import MPCADD


#from pyNastran.bdf.dev_vectorized.cards.coordinateSystems import (CORD1R, CORD1C, CORD1S,
                                                                  #CORD2R, CORD2C, CORD2S, CORD3G)
#from .cards.coordinate_systems import (CORD1R, CORD1C, CORD1S,
#                                      CORD2R, CORD2C, CORD2S, CORD3G) old...
from pyNastran.bdf.cards.params import PARAM
from pyNastran.bdf.case_control_deck import CaseControlDeck
from .bdf_methods import BDFMethods
from .bdf_interface.get_methods import GetMethods
from .bdf_interface.add_card import AddCard
from .bdf_interface.utils import wipe_empty_fields
from .bdf_interface.assign_type import interpret_value
from .bdf_interface.write_mesh import WriteMesh
from .bdf_interface.cross_reference import XRefMesh

# old
from pyNastran.bdf.bdf_interface.bdf_card import BDFCard

# sets
from pyNastran.bdf.cards.bdf_sets import SET1, SET3

def class_obj_defaultdict(class_obj, *args, **kwargs):

    class ClassObjDefaultDict(dict):
        def __init__(self):
            dict.__init__(self)

        def __getitem__(self, key):
            #print('getting key=%r' % key)
            #print('dir(d) = %s' % dir(self))
            try:
                val = self[key]
                return val
            except KeyError:
                val = class_obj(self.args) # , self.kwargs
                self[key] = [val]
                return self[key]

        def __setitem__(self, key, value):
            #print('setting key=%r' % key)
            #print('dir(d) =', dir(self))
            try:
                val = self[key].append(value)
            except KeyError:
                self[key] = [class_obj(args)] # , self.kwargs

    return ClassObjDefaultDict


class BDF(BDFMethods, GetMethods, AddCard, WriteMesh, XRefMesh):
    """
    NASTRAN BDF Reader/Writer/Editor class.
    """
    #: this is a nastran model
    model_type = 'nastran'

    #: required for sphinx bug
    #: http://stackoverflow.com/questions/11208997/autoclass-and-instance-attributes
    #__slots__ = ['_is_dynamic_syntax']
    def __init__(self, debug=True, log=None, precision='double'):
        """
        Initializes the BDF object

        Parameters
        debug : bool; default=True
            used to set the logger if no logger is passed in
        log : logger; default=None -> logger
            a python logging module object
        precision : str; default='double'
            'single'/'float32'
            'double'/'float64'
            that is used by all the objects
        """
        assert debug in [True, False], 'debug=%r' % debug

        self.set_precision(precision)

        # file management parameters
        self._ifile = -1
        self.include_dir = ''
        self.active_filename = None
        self.active_filenames = []
        #self.used_filenames = []
        #self._stored_Is = {}
        #self._stored_lines = {}
        #self._stored_comments = {}
        #self._line_streams = {}
        #self._card_streams = {}
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
        self.foundEndData = False

        #: allows the BDF variables to be scoped properly (i think...)
        GetMethods.__init__(self)
        AddCard.__init__(self)
        BDFMethods.__init__(self)
        WriteMesh.__init__(self)
        XRefMesh.__init__(self)

        #: useful in debugging errors in input
        self.debug = debug

        #: flag that allows for OpenMDAO-style optimization syntax to be used
        self._is_dynamic_syntax = False

        self._element_type_to_element_name_mapper = {
            1 : 'CROD',
            2 : 'CBEAM', 34 : 'CBAR',
            4 : 'CSHEAR',
            10 : 'CONROD',
            3 : 'CTUBE',

            11 : 'CELAS1', 12 : 'CELAS2', 13 : 'CELAS3', 14 : 'CELAS4',
            15 : 'CDAMP1', 16 : 'CDAMP2', 17 : 'CDAMP3', 18 : 'CDAMP4',
            19 : 'CVISC', 20 : 'CBUSH', 21 : 'CBUSH1D', 23 : 'CBUSH2D',
            30 : 'CMASS1', 31 : 'CMASS2', 32 : 'CMASS3', 33 : 'CMASS4',

            35 : 'CONM1', # ???
            36 : 'CONM2', # ???

            73 : 'CTRIA3', 74: 'CTRIA6', 75 : 'CTRIAX', 76: 'CTRIAX6',
            144 : 'CQUAD4', 145: 'CQUAD8', 146: 'CQUAD', 147: 'CQUADX',

            60 : 'CTETRA4', 61 : 'CTETRA10',
            62 : 'CPENTA6', 63 : 'CPENTA15',
            64 : 'CHEXA8', 65 : 'CHEXA20',

            # ???
            100 : 'CAERO1'
        }
        self._element_name_to_element_type_mapper = {
            v:k for k, v in iteritems(self._element_type_to_element_name_mapper)}

        #: lines that were rejected b/c they were for a card that isnt supported
        self.rejects = []

        #: cards that were created, but not processed
        self.reject_cards = []

        #: list of execive control deck lines
        self.executive_control_lines = []

        #: list of case control deck lines
        self.case_control_lines = []

        #: should the cards be echoed to the console as they are read; bool
        #: ECHOON and ECHOOFF will toggle this in the BDF
        self.echo = False

        self._is_cards_dict = True

        # file management parameters
        self.active_filenames = []
        self.active_filename = None
        self.include_dir = ''
        self.dumplines = False

        self.__init_attributes()

        #: the list of possible cards that will be parsed
        self.cards_to_read = set([
            'ECHOON', 'ECHOOFF',
            'PARAM',
            'GRID', 'GRDSET', 'SPOINT', 'EPOINT', 'POINT', 'POINTAX', # 'RINGAX',

            # elements
            'CONM1', 'CONM2', 'CMASS1', 'CMASS2', 'CMASS3', 'CMASS4',

            # elements
            'CELAS1', 'CELAS2', 'CELAS3', 'CELAS4',
            # 'CELAS5',
            'CBUSH', 'CBUSH1D', 'CBUSH2D',

            # dampers
            'CDAMP1', 'CDAMP2', 'CDAMP3', 'CDAMP4', 'CDAMP5',
            'CFAST',

            'CBAR', 'CROD', 'CTUBE', 'CBEAM', 'CBEAM3', 'CONROD', 'CBEND',
            'CTRIA3', 'CTRIA6', 'CTRIAR', 'CTRIAX', 'CTRIAX6',
            'CQUAD4', 'CQUAD8', 'CQUADR', 'CQUADX', 'CQUAD',

            'CTETRA', 'CPENTA', 'CHEXA',
            'CTETRA4', 'CPYRAM5', 'CPENTA6', 'CHEXA8',
            'CTETRA10', 'CPYRAM13', 'CPENTA15', 'CHEXA20',

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
            'DLOAD', 'SLOAD', 'ACSRCE', 'TLOAD1', 'TLOAD2', 'RLOAD1', 'RLOAD2',
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
            'TEMP', 'TEMPD', 'TEMPP1',
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

            # optimization cards
            'DEQATN',
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
            'SESET',

            # tables
            #'DTABLE', 'TABLEHT', 'TABRNDG',
            'TABLED1', 'TABLED2', 'TABLED3', 'TABLED4',
            'TABLEM1', 'TABLEM2', 'TABLEM3', 'TABLEM4',
            'TABDMP1',
            'TABLES1', 'TABLEST',
            'TABRND1', 'TABRNDG',
            #------------------------------------------------------------------

            # initial conditions - sid (set ID)
            #'TIC',  (in tables.py)

            #: methods - .. todo:: EIGRL not done???
            'EIGB', 'EIGR', 'EIGRL',

            #: cMethods - .. todo:: EIGC not done???
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
        self.special_cards = ['DEQATN', '/']

    def set_precision(self, precision='double'):
        """
        Sets the float precision.

        :param precision:  string of 'single'/'float32' or
          'double'/'float64' that is used by all the objects
        """
        if precision in ('double', 'float64'):
            self.float = 'float64'
        elif precision == ('single', 'float32'):
            self.float = 'float32'
        else:
            raise NotImplementedError('precision=%r' % precision)

    def disable_cards(self, cards):
        """
        Method for removing broken cards from the reader

        Parameters
        ----------
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

        model = self
        self.grid = GRID(model)
        self.point = POINT(model)
        self.grdset = GRDSET(model)
        self.spoint = SPOINT(model)
        self.epoint = EPOINT(model)
        self.pointax = POINTAX(model)

        self.coords = Coord(model)

        #: stores rigid elements (RBE2, RBE3, RJOINT, etc.)
        #self.rigid_elements = {}
        self.rbe2 = RBE2(model)
        self.rbe3 = RBE3(model)

        #self.properties_spring = PropertiesSpring(model)
        #self.properties_rod = PropertiesRod(v)
        #self.properties_bar = PropertiesBar(model)
        #self.properties_solid = PropertiesSolid(model)

        #self.elements_spring = ElementsSpring(model)
        #self.elements_rod = ElementsRod(model)
        #self.elements_bar = ElementsBar(model)

        # shell
        #: stores PSHELL, PCOMP, PCOMPG
        self.properties_shell = PropertiesShell(model)
        #: stores CTRIA3, CTRIA6, CQUAD4, CQUAD8
        self.elements_shell = ElementsShell(model)

        class DummyCard(object):
            def __init__(self, model):
                pass
            def add(self, card, comment=''):
                pass
        # bush
        self.cbush = CBUSH(model)
        self.cbush1d = DummyCard(model)
        self.cbush2d = DummyCard(model)
        self.pbush = PBUSH(model)

        # shear
        #: stores CSHEAR
        self.cshear = CSHEAR(model)
        #: stores PSHEAR
        self.pshear = PSHEAR(model)

        # spring
        self.elements_spring = ElementsSpring(model)
        self.pelas = PELAS(model)

        # rods
        self.conrod = CONROD(model)
        self.prod = PROD(model)
        self.crod = CROD(model)
        self.ptube = PTUBE(model)
        self.ctube = CTUBE(model)

        # mass
        #: stores CONM1, CONM2, CMASS1, CMASS2, CMASS3, CMASS4, CMASS5, PMASS
        self.mass = Mass(model)

        # bars
        #: stores CBAR
        self.cbar = CBAR(model)
        #: stores PBAR, PBARL
        self.properties_bar = PropertiesBar(model)

        # stores CBAROR
        #self.cbaror = None
        self.cbaror = CBAROR()

        # beams
        #: stores CBEAM
        self.cbeam = CBEAM(model)
        #: stores PBEAM, PBEAML
        self.properties_beam = PropertiesBeam(model)
        # stores CBEAMOR
        self.cbeamor = None

        # disabled 1D
        self.cbeam3 = None
        self.pbeam3 = None
        self.cbend = None
        self.pbend = None

        # solids
        #: stores CTETRA4, CPENTA6, CHEXA8, CTETRA10, CPENTA15, CHEXA20
        self.elements_solid = ElementsSolid(self)
        #: stores PSOLID, PLSOLID
        self.properties_solid = PropertiesSolid(self)

        # control structure
        self.elements = Elements(self)
        self.properties = Properties(self)

        #===================================
        # methods
        self.eigrl = {}
        self.eigb = {}
        self.eigc = {}

        # dynamic
        self.nlpcis = {}
        #===================================
        # vectorized: aero

        #: stores CAERO1, CAERO2, CAERO3, CAERO4, CAERO5
        self.caero = CAero(model)
        #: stores PAERO1, PAERO2, PAERO3, PAERO4, PAERO5
        self.paero = PAero(model)

        self.caero1 = self.caero.caero1
        self.paero1 = self.paero.paero1

        self.spline1 = SPLINE1(model)
        #self.spline2 = SPLINE2(model)
        #self.spline3 = SPLINE3(model)
        #self.spline4 = SPLINE4(model)
        #self.spline5 = SPLINE5(model)

        self.aero = AERO(model)
        self.aeros = AEROS(model)

        #===================================
        #: stores coordinate systems
        #self.coord = Coord(model)

        # loads
        #self.loadcase = LoadCase(model)

        #: stores LOAD, FORCE, MOMENT, etc.
        self.loads = Loads(model)
        self.temps = TEMPs(model)

        #: stores MAT1, MAT2, MAT3,...MAT10 (no MAT4, MAT5)
        self.materials = Materials(model)

        self.__define_unvectorized()

    def __define_unvectorized(self):
        #: stores TRIM
        self.trim = {}
        self.trims = {}
        #===================================
        # optimization
        #self.dconstr = DCONSTR(model)

        #===================================
        #: stores MATS1
        self.materialDeps = {}
        #: stores the CREEP card
        self.creep_materials = {}

        # loads
        #self.gusts  = {} # Case Control GUST = 100
        #self.random = {} # Case Control RANDOM = 100

        # --------------------------- constraints ----------------------------
        #: stores SUPORT1s
        #self.constraints = {} # suport1, anything else???
        self.suports = []  # suport, suport1

        #: stores SPCADD,SPC,SPC1,SPCD,SPCAX
        #self.spcObject = ConstraintObject()
        #: stores MPCADD,MPC
        #self.mpcObject = ConstraintObject()

        # these work
        #self.spc = defaultdict(SPC)
        self.spc = {} #class_obj_defaultdict(SPC, model)
        self.spcd = {} #class_obj_defaultdict(SPCD, model)
        self.spc1 = {} #class_obj_defaultdict(SPC1, model)
        self.spcadd = {}
        self.mpc = {}  # the new form, not added...
        self.mpcadd = {}

        # these don't...
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
        self._dmig_temp = defaultdict(list)
        self.dmis = {}
        self.dmigs = {}
        self.dmijs = {}
        self.dmijis = {}
        self.dmiks = {}
        self.dequations = {}

        # ----------------------------------------------------------------
        #: SETx
        self.set1 = {}
        self.set2 = {}
        self.set3 = {}
        self.aset = {}
        self.bset = {}
        self.cset = {}
        self.qset = {}
        #: SESETx
        self.seset = {}

        # ----------------------------------------------------------------
        #: tables
        self.tables = {}
        #: random_tables
        self.random_tables = {}

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
        self.nlpci = {}
        #: stores NLPARM
        self.nlparm = {}
        self.nlparms = {}
        #: stores TSTEPs
        self.tstep = {}
        self.tsteps = {}
        #: stores TSTEPNL
        self.tstepnl = {}
        self.tstepnls = {}
        # --------------------------- aero defaults --------------------------
        # aero cards
        #: stores CAEROx
        self.caeros = {}
        #: stores PAEROx
        self.paeros = {}
        #: stores AERO
        #self.aero = {}
        #: stores AEROS
        #self.aeros = {}

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

        # ------------------------- thermal defaults -------------------------
        # BCs
        #: stores thermal boundary conditions - CONV,RADBC
        self.bcs = {}  # e.g. RADBC
        #: defines the MAT4, MAT5, MATT4, etc.  .. todo:: verify MATT4
        self.thermal_materials = {}

        #: stores PHBDY
        self.phbdys = {}
        #: stores convection properties - PCONV, PCONVM ???
        self.convection_properties = {}

        # -------------------------contact cards-------------------------------
        self.bcrparas = {}
        self.bctadds = {}
        self.bctparas = {}
        self.bctsets = {}
        self.bsurf = {}
        self.bsurfs = {}

    def read_bdf(self, bdf_filename=None,
                 xref=True, punch=False, encoding=None):
        """
        Read method for the bdf files

        Parameters
        ----------
        bdf_filename : str / None
            the input bdf (default=None; popup a dialog)
        xref :  bool
            should the bdf be cross referenced (default=True)
        punch : bool
            indicates whether the file is a punch file (default=False)
        encoding : str
            the unicode encoding (default=None; system default)

        .. todo:: this is out of date
        >>> bdf = BDF()
        >>> bdf.read_bdf(bdf_filename, xref=True)
        >>> g1 = bdf.Node(1)
        >>> print(g1.get_position())
        [10.0, 12.0, 42.0]
        >>> bdf.write_bdf(bdf_filename2)
        >>> print(bdf.card_stats())
        ---BDF Statistics---
        SOL 101
        bdf.nodes = 20
        bdf.elements = 10
        etc.
        """
        self._read_bdf_helper(bdf_filename, encoding, punch)


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

        if self._is_cards_dict:
            cards, card_count = self.get_bdf_cards_dict(bulk_data_lines)
        else:
            cards, card_count = self.get_bdf_cards(bulk_data_lines)
        self._parse_cards(cards, card_count)

        #if not self.inspect:
            #import time
            #t0 = time.time()
            #bdf_temp = BDF(debug=False)
            #bdf_temp.inspect = True
            #bdf_temp.add_card = bdf_temp._add_card
            #bdf_temp.add_reject = bdf_temp._add_reject
            #bdf_temp.read_bdf(bdf_filename=bdf_filename, include_dir=include_dir, xref=False, punch=punch)
            #print('card_count = %s' % bdf_temp.card_count)
            #self.allocate(bdf_temp.card_count)
            #del bdf_temp
            #print('dt = %s' % (time.time() - t0))

        #try:
            #self._open_file(self.bdf_filename)
            #self.log.debug('---starting BDF.read_bdf of %s---' % self.bdf_filename)
            #if not punch:
                #self.log.debug('---reading executive & case control decks---')
                #self._read_executive_control_deck()
                #self._read_case_control_deck()
            #else:
                #self.log.debug('---skipping executive & case control decks---')

            #self._read_bulk_data_deck()
            #self.build(xref=xref)
            #self._xref = xref
            #self._cleanup_file_streams()
        #except:
            #self._cleanup_file_streams()
            #raise
        self.build(xref=xref)
        self.fill_dmigs()
        self.log.debug('---finished BDF.read_bdf of %s---' % self.bdf_filename)

    def _read_bdf_helper(self, bdf_filename, encoding, punch):
        """creates the file loading if bdf_filename is None"""
        if encoding is None:
            encoding = sys.getdefaultencoding()
        self._encoding = encoding
        if bdf_filename is None:
            #from pyNastran.utils.gui_io import load_file_dialog
            wildcard_wx = "Nastran BDF (*.bdf; *.dat; *.nas; *.pch, *.ecd)|" \
                "*.bdf;*.dat;*.nas;*.pch|" \
                "All files (*.*)|*.*"
            wildcard_qt = "Nastran BDF (*.bdf *.dat *.nas *.pch *.ecd);;All files (*)"
            title = 'Please select a BDF/DAT/PCH/ECD to load'
            bdf_filename = load_file_dialog(title, wildcard_wx, wildcard_qt)[0]
            assert bdf_filename is not None, bdf_filename

        if not os.path.exists(bdf_filename):
            msg = 'cannot find bdf_filename=%r\n%s' % (bdf_filename, print_bad_path(bdf_filename))
            raise IOError(msg)
        if bdf_filename.lower().endswith('.pch'):  # .. todo:: should this be removed???
            punch = True

        #: the active filename (string)
        self.bdf_filename = bdf_filename

        #: is this a punch file (no executive control deck)
        self.punch = punch

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
            card.finalize()

        self._dmig_temp = defaultdict(list)

    def allocate(self, card_count, cards):
        """
        Sets the size of the card objects.

        Parameters
        ----------
        card_count : dict
            dictionary of {card_name : ncards}
            card_name : string
            ncards : int

        .. note::  Sometimes there are cards that have 2 cards per line.
                   Depending on the card (e.g. CELASx)
        .. note::  Some cards (e.g. CTETRA4, CTETRA10) have a distinction
                    despite all being CTETRAs.  Card_count considers
                    CTETRA4s and CTETRA10s.

        .. note::  Not all cards require allocation (e.g. CORD1x), so we
           can be a more lax about how we calculate the card_count.  The
           CORD1x cards can have 2 cards per line, so not needing this
           is convenient.
        .. note::  In this method, we only allocate the macro groups.
           For example, since a CTRIA3 is an element, it is set with the
           other elements.
        """
        self.grid.allocate(card_count)
        self.coords.allocate(card_count=card_count)

        card_splits = [
            # card_name, card_name1, nfields1, card_name2
            ('CTETRA', 'CTETRA4', 7, self.elements_solid.ctetra4, 'CTETRA10', self.elements_solid.ctetra10),
            #('CYPRAM', 'CPYRAM5', self.elements_solid.cpyram5, 8, 'CPYRAM13', self.elements_solid.cpyram12),
            ('CPENTA', 'CPENTA6', self.elements_solid.cpenta6, 9, 'CPENTA15', self.elements_solid.cpenta15),
            ('CHEXA', 'CHEXA8', self.elements_solid.chexa8, 11, 'CHEXA20', self.elements_solid.chexa20),
        ]
        for (card_name, card_name1, slot1, nfields1, card_name2, slot2) in card_splits:
            if card_name in card_count:
                cards1 = [] # CTETRA4, CPENTA6, CHEXA8, CPYRAM5
                cards2 = [] # CTETRA10, CPENTA15, CHEXA20, CPYRAM15
                cardsi = cards[card_name]
                #assert len(cardsi) == 2, cardsi
                is_list = False
                has_none = False
                for comment, card_lines in cardsi:
                    card_obj, card = self.create_card_object(card_lines, card_name,
                                                             is_list=is_list, has_none=has_none)
                    if len(card_obj) == nfields1:
                        cards1.append((comment, card_obj))
                    else:
                        cards2.append((comment, card_obj))

                count1 = len(cards1)
                if count1:
                    card_count[card_name1] = count1
                    #cards[card_name1] = cards1
                    self.log.info(' allocating %s->%s' % (card_name1, count1))
                    slot1.allocate(count1)
                    for comment, card_obj in cards1:
                        self._add_card_object(card_name1, card_obj, comment=comment)

                count2 = len(cards2)
                if count2:
                    slot2.allocate(count2)
                    card_count[card_name1] = count2
                    self.log.info(' allocating %s->%s' % (card_name2, count2))
                    for comment, card_obj in cards2:
                        self._add_card_object(card_name2, card_obj, comment=comment)

                del card_count[card_name]
                del cards[card_name]
        print(cards.keys())




        #if card_name in ['CTETRA', 'CPENTA', 'CHEXA']:
            #card = self.process_card(card_lines)
            #if card_name == 'CTETRA':
                #if len(card) == 7:
                    #card_name = 'CTETRA4'
                #else:
                    #card_name = 'CTETRA10'
            #elif card_name == 'CPENTA':
                #if len(card) == 9:
                    #card_name = 'CPENTA6'
                #else:
                    #card_name = 'CPENTA15'
            #elif card_name == 'CHEXA':
                #if len(card) == 11:
                    #card_name = 'CHEXA8'
                #else:
                    #card_name = 'CHEXA20'
        #self._increase_card_count(card_name)


        self.elements.allocate(card_count)
        self.materials.allocate(card_count)
        self.loads.allocate(card_count)
        self.temps.allocate(card_count)

        if 'CAERO1' in card_count:
            self.caero1.allocate(card_count['CAERO1'])
        if 'PAERO1' in card_count:
            self.paero1.allocate(card_count['PAERO1'])
        if 'SPLINE1' in card_count:
            self.spline1.allocate(card_count['SPLINE1'])
        if 'AERO' in card_count:
            self.aero.allocate(card_count['AERO'])
        if 'AEROS' in card_count:
            self.aeros.allocate(card_count['AEROS'])

    def get_bdf_cards(self, bulk_data_lines):
        """Parses the BDF lines into a list of card_lines"""
        cards = []
        #cards = defaultdict(list)
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
                    if self.echo:
                        self.log.info('Reading %s:\n' %
                                      old_card_name + full_comment + ''.join(card_lines))

                    # old dictionary version
                    # cards[old_card_name].append([full_comment, card_lines])

                    # new list version
                    cards.append([old_card_name, full_comment, card_lines])

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
                        nleftover = nlines - i - 1
                        msg = 'exiting due to ENDDATA found with %i lines left' % nleftover
                        self.log.debug(msg)
                    return cards, card_count
                #print("card_name = %s" % card_name)

            comment = _clean_comment(comment)
            if line.rstrip():
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

            # old dictionary version
            #cards[old_card_name].append([backup_comment + full_comment, card_lines])

            # new list version
            cards.append([old_card_name, backup_comment + full_comment, card_lines])
            card_count[old_card_name] += 1
        return cards, card_count

    def get_bdf_cards_dict(self, bulk_data_lines):
        """Parses the BDF lines into a list of card_lines"""
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
                    if self.echo:
                        self.log.info('Reading %s:\n' %
                                      old_card_name + full_comment + ''.join(card_lines))

                    # old dictionary version
                    cards[old_card_name].append([full_comment, card_lines])

                    # new list version
                    #cards.append([old_card_name, full_comment, card_lines])

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
                        nleftover = nlines - i - 1
                        msg = 'exiting due to ENDDATA found with %i lines left' % nleftover
                        self.log.debug(msg)
                    return cards, card_count
                #print("card_name = %s" % card_name)

            comment = _clean_comment(comment)
            if line.rstrip():
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

            # old dictionary version
            cards[old_card_name].append([backup_comment + full_comment, card_lines])

            # new list version
            #cards.append([old_card_name, backup_comment + full_comment, card_lines])
            card_count[old_card_name] += 1
        return cards, card_count

    def update_solution(self, sol, method, isol_line):
        """
        Updates the overall solution type (e.g. 101,200,600)

        Parameters
        ----------
        sol : int
            the solution type (101, 103, etc)
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

    def _is_case_control_deck(self, line):
        """
        .. todo:: not done...
        """
        line_upper = line.upper().strip()
        if 'CEND' in line.upper():
            raise SyntaxError('invalid Case Control Deck card...CEND...')
        if '=' in line_upper or ' ' in line_upper:
            return True
        for card in self.unique_bulk_data_cards:
            len_card = len(card)
            if card in line_upper[:len_card]:
                return False
        return True

    def is_reject(self, card_name):
        """
        Should the card be read?

        Parameters
        ----------
        card_name : str
            the card_name -> 'GRID'

        Returns
        -------
        is_reject : bool
            True/False
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

    def _get_card_name(self, lines):
        """
        Returns the name of the card defined by the provided lines

        Parameters
        lines : List[str]
            the lines of the card

        Returns
        -------
        cardname : str
            the name of the card
        """
        card_name = lines[0][:8].rstrip('\t, ').split(',')[0].split('\t')[0].strip('*\t ')
        if len(card_name) == 0:
            return None
        if ' ' in card_name or len(card_name) == 0:
            msg = 'card_name=%r\nline=%r in filename=%r is invalid' \
                  % (card_name, lines[0], self.active_filename)
            raise RuntimeError(msg)
        return card_name.upper()

    def _add_reject(self, comment, lines):
        """
        The dummy add_reject method used by self.inspect=False
        """
        pass

    def add_reject(self, comment, lines):
        """
        Appends self.rejects with the list of lines in the card

        Parameters
        ----------
        comment : str
            the card comment
        lines : List[str]
            the lines that are part of the card
        """
        if comment:
            self.rejects.append([comment])
        self.rejects.append(lines)

    def _dump_file(self, bdf_dump_filename, lines, i):
        """
        Writes a BDF up to some failed line index

        Parameters
        ----------
        bdf_dump_filename : str
            the bdf filename to dump
        lines : List[str]
            the entire list of lines
        i : int
            the last index to write
        """
        with codec_open(_filename(bdf_dump_filename),
                        'w', encoding=self._encoding) as crash_file:
            for line in lines[:i]:
                crash_file.write(line)

    def _increase_card_count(self, card_name):
        """
        Used for testing to check that the number of cards going in is the
        same as each time the model is read verifies proper writing of cards

        Parameters
        ----------
        card_name : str
            the card_name -> 'GRID'

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

    def process_card(self, card_lines):
        """
        Parameters
        ----------
        card_lines : List[str]
            the lines that are part of the current card

        Returns
        -------
        card : list[str]
            a list of (typically) string fields that is properly sized

        .. note:: If dynamic syntax is used, then the fields may be
                  ints/floats.

        .. note:: a properly sized card's last field is a value (not None)
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
        """
        Creates a BDFCard object, which is really just a list that
        allows indexing past the last field
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

    def write_sorted_card(self, card, n):
        """
        1-        CELAS2  1       3.      1       1       2       1               7.0
        """
        if self.echo:
            msg = '%21s-%8s' % (n, '')
            for field in card.fields():
                if field is None:
                    field = ''
                msg += '%-8s' % field
            self.f06.write('%-110s\n' % msg)
        #return msg

    def add_card(self, card_lines, card_name, comment='', is_list=True, has_none=True):
        """
        Adds a card object to the BDF object.

        Parameters
        ----------
        card_lines : List[str]
            the list of the card fields

         >>> ['GRID,1,2',]  # (is_list = False)
         >>> ['GRID',1,2,]  # (is_list = True; default)

        card_name : str
            the card_name -> 'GRID'
        comment : str; default=''
            an optional the comment for the card
        is_list : bool; default=True
            changes card_lines from a list of lines to a list of fields

        Returns
        -------
        card_object : BDFCard()
            the card object representation of card

        .. note:: this is a very useful method for interfacing with the code
        .. note:: the cardObject is not a card-type object...so not a
                  GRID card or CQUAD4 object.  It's a BDFCard Object.
                  However, you know the type (assuming a GRID), so just
                  call the *mesh.Node(nid)* to get the Node object that
                  was just created.
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
            assert isinstance(card, list), 'type=%s; card=\n%r' % (type(card), card)
            print('***card', card)
            card_obj = BDFCard(card)

        self._add_card_object(card_name, card_obj, comment=comment)

    def _add_card_object(self, card_name, card_obj, comment=''):
        self.log.info('_add_card_object; card_name=%s' % card_name)
        # function that gets by name the initialized object (from global scope)

        name = card_name # card[0]
        ##print("name = %r" % name)
        #self.write_sorted_card(card_obj, icard)

        if card_name == 'ECHOON':
            self.echo = True
            return
        elif card_name == 'ECHOOFF':
            self.echo = False
            return

        if self.echo:
            print(print_card_8(card_obj).rstrip())

        #if icard % 10000 == 0:
            #self.log.debug('icard = %i' % icard)

        if name == 'PARAM':
            param = PARAM.add_card(card_obj, comment=comment)
            self.add_PARAM(param)
        elif name == 'BCRPARA':
            pass
        elif name == 'BCTADD':
            pass
        elif name == 'BCTSET':
            pass
        elif name == 'BSURFS':
            pass
        elif name == 'BSURF':
            #self.bsurf.add(card_obj, comment=comment)
            pass

        #========================
        # cshear / pshear
        elif name == 'CSHEAR':
            self.cshear.add(card_obj, comment=comment)
        elif name == 'PSHEAR':
            self.pshear.add(card_obj, comment=comment)
        #========================
        # elements_shell
        elif name == 'CTRIA3':
            self.elements_shell.add_ctria3(card_obj, comment=comment)
        elif name == 'CTRIA6':
            self.elements_shell.add_ctria6(card_obj, comment=comment)

        elif name == 'CQUAD':
            self.elements_shell.add_cquad(card_obj, comment=comment)
        elif name == 'CQUADX':
            self.elements_shell.add_cquadx(card_obj, comment=comment)
        elif name == 'CQUADR':
            self.elements_shell.add_cquadr(card_obj, comment=comment)
        elif name == 'CQUAD4':
            self.elements_shell.add_cquad4(card_obj, comment=comment)
        elif name == 'CQUAD8':
            self.elements_shell.add_cquad8(card_obj, comment=comment)
        elif name == 'CQUAD9':
            self.elements_shell.add_cquad9(card_obj, comment=comment)

        # properties shell
        elif name == 'PSHELL':
            self.properties_shell.add_pshell(card_obj, comment=comment)
        elif name == 'PCOMP':
            self.properties_shell.add_pcomp(card_obj, comment=comment)
        elif name == 'PCOMPG':
            self.properties_shell.add_pcompg(card_obj, comment=comment)

        #========================
        # rotation elements/loads
        elif name == 'CTRIAX':
            self.elements_shell.add_ctriax(card_obj, comment=comment)
        elif name == 'CTRIAX6':
            self.elements_shell.add_ctriax6(card_obj, comment=comment)
        #elif name == 'PLOADX1':
            #self.ploadx1.add(card_obj, comment=comment)

        #========================
        # elements_solid
        elif name == 'CTETRA4':
            self.elements_solid.add_ctetra4(card_obj, comment=comment)
        elif name == 'CTETRA10':
            self.elements_solid.add_ctetra10(card_obj, comment=comment)
        elif name == 'CPENTA6':
            self.elements_solid.add_cpenta6(card_obj, comment=comment)
        elif name == 'CPENTA15':
            self.log.debug('len(CPENTA) = %s' % len(card_obj))
            self.elements_solid.add_cpenta15(card_obj, comment=comment)
        elif name == 'CHEXA8':
            self.elements_solid.add_chexa8(card_obj, comment=comment)
        elif name == 'CHEXA20':
            self.log.debug('len(CHEXA) = %s' % len(card_obj))
            self.elements_solid.add_chexa20(card_obj, comment=comment)

        elif name == 'PSOLID':
            self.properties_solid.add_psolid(card_obj, comment=comment)
        elif name == 'PLSOLID':
            self.properties_solid.add_plsolid(card_obj, comment=comment)
        #========================
        # conrod/rod
        elif name == 'CROD':
            self.crod.add(card_obj, comment=comment)
        elif name == 'CONROD':
            self.conrod.add(card_obj, comment=comment)
        elif name == 'PROD':
            self.prod.add(card_obj, comment=comment)

        #========================
        # tube
        elif name == 'CTUBE':
            self.ctube.add(card_obj, comment=comment)
        elif name == 'PTUBE':
            self.ptube.add(card_obj, comment=comment)

        #========================
        # thermal boundary conditions
        elif name == 'CHBDYE':
            pass
        elif name == 'CHBDYG':
            pass
        elif name == 'CHBDYP':
            pass
        elif name == 'PHBDY':
            pass
        #========================
        # convection
        elif name == 'CONV':
            pass
        elif name == 'CONVM':
            pass
        elif name == 'PCONV':
            pass
        elif name == 'PCOMVM':
            pass
        #========================
        # radiation
        elif name == 'RADM':
            pass
        elif name == 'RADBC':
            pass
        #========================
        # thermal load
        elif name == 'QBDY1':
            pass
        elif name == 'QBDY2':
            pass
        #========================
        # design_variables
        elif name == 'DESVAR':
            #self.desvar.add(card_obj, comment=comment)
            pass
        elif name == 'DCONSTR':
            #self.dconstr.add(card_obj, comment=comment)
            pass
        elif name == 'DOPTPRM':
            #self.doptprm.add(card_obj, comment=comment)
            pass
        elif name == 'DVPREL1':
            pass
        elif name == 'DRESP1':
            pass
        elif name == 'DRESP2':
            pass

        #========================
        #aero...
        elif name == 'MONPNT1':
            pass
        elif name == 'AECOMP':
            pass
        elif name == 'AELINK':
            pass
        elif name == 'AESTAT':
            pass
        elif name == 'AESURF':
            pass
        #elif name == 'AESURFS':
            #pass
        elif name == 'AELIST':
            pass
        elif name == 'AERO':
            self.add_AERO(card_obj)
            #self.aero.add(card_obj, comment=comment)
        elif name == 'AEROS':
            #self.aeros.add(card_obj, comment=comment)
            pass

        elif name == 'FLUTTER':
            #self.flutter.add(card_obj, comment=comment)
            pass
        elif name == 'FLFACT':
            #self.flfact.add(card_obj, comment=comment)
            pass

        elif name == 'MKAERO1':
            #self.mkaero1.add(card_obj, comment=comment)
            pass
        #========================
        elif name == 'SET1':
            set1 = SET1(card_obj)
            key = set1.sid
            if key in self.set1:
                set1_base = self.set1[key]
                set1_base.union(set1)
            else:
                self.set1[key] = set1
        elif name == 'SET3':
            set3 = SET3(card_obj)
            key = set3.sid
            if key in self.set3:
                set3_base = self.set3[key]
                set3_base.union(set3)
            else:
                self.set3[key] = set3
        elif name == 'ASET1':
            pass
        elif name == 'QSET1':
            pass
        #========================
        elif name == 'SPLINE1':
            self.spline1.add(card_obj, comment=comment)
        elif name == 'SPLINE2':
            self.spline2.add(card_obj, comment=comment)
        elif name == 'SPLINE3':
            self.spline3.add(card_obj, comment=comment)
        elif name == 'SPLINE4':
            self.spline4.add(card_obj, comment=comment)
        elif name == 'SPLINE5':
            self.spline5.add(card_obj, comment=comment)
        #========================
        # caero/paero
        elif name == 'CAERO1':
            self.caero.add_caero1(card_obj, comment=comment)
        elif name == 'CAERO2':
            self.caero.add_caero2(card_obj, comment=comment)
        elif name == 'CAERO3':
            self.caero.add_caero3(card_obj, comment=comment)
        elif name == 'CAERO4':
            self.caero.add_caero4(card_obj, comment=comment)
        elif name == 'CAERO5':
            self.caero.add_caero5(card_obj, comment=comment)

        elif name == 'PAERO1':
            self.paero.add_paero1(card_obj, comment=comment)
        elif name == 'PAERO2':
            self.paero.add_paero2(card_obj, comment=comment)
        elif name == 'PAERO3':
            self.paero.add_paero3(card_obj, comment=comment)
        elif name == 'PAERO4':
            self.paero.add_paero4(card_obj, comment=comment)
        elif name == 'PAERO5':
            self.paero.add_paero5(card_obj, comment=comment)

        elif name == 'SPLINE1':
            self.spline1.add(card_obj, comment=comment)

        elif name == 'TRIM':
            trim = TRIM(self)
            trim.add(card_obj, comment=comment)
            self.trim[trim.trim_id] = trim
        #========================
        # bushing
        elif name == 'CBUSH':
            self.cbush.add(card_obj, comment=comment)
        elif name == 'CBUSH1D':
            self.cbush1d.add(card_obj, comment=comment)
        elif name == 'CBUSH2D':
            self.cbush2d.add(card_obj, comment=comment)
        elif name == 'CBUSH':
            self.cbush.add(card_obj, comment=comment)
        elif name == 'PBUSH':
            self.pbush.add(card_obj, comment=comment)
        elif name == 'PBUSHT':
            self.pbusht.add(card_obj, comment=comment)

        #========================
        # fast
        elif name == 'PFAST':
            self.pfast.add(card_obj, comment=comment)
        #========================
        # springs
        elif name == 'PELAS':
            self.pelas.add(card_obj, 0, comment)
            self.pelas.add(card_obj, 1, comment)
        elif name == 'CELAS1':
            self.elements_spring.add_celas1(card_obj, comment)
        elif name == 'CELAS2':
            self.elements_spring.add_celas2(card_obj, comment)
        elif name == 'CELAS3':
            self.elements_spring.add_celas3(card_obj, comment)
        elif name == 'CELAS4':
            self.elements_spring.add_celas4(card_obj, comment)
        #========================
        # dampers
        elif name == 'PDAMP':
            self.pdamp.add(card_obj, comment=comment)
        elif name == 'PDAMPT':
            self.pdampt.add(card_obj, comment=comment)
        elif name == 'CDAMP1':
            self.elements_damper.cdamp1.add(card_obj, comment=comment)
        elif name == 'CDAMP2':
            self.elements_damper.cdamp2.add(card_obj, comment=comment)
        elif name == 'CDAMP3':
            self.elements_damper.cdamp3.add(card_obj, comment=comment)
        elif name == 'CDAMP4':
            self.elements_damper.cdamp4.add(card_obj, comment=comment)

        #========================
        # bars
        elif name == 'CBAR':
            self.cbar.add(card_obj, comment=comment)
        elif name == 'CBAROR':
            self.cbaror.add(card_obj, comment=comment)
        elif name == 'PBAR':
            self.properties_bar.add_pbar(card_obj, comment=comment)
        elif name == 'PBARL':
            self.properties_bar.add_pbarl(card_obj, comment=comment)

        # beams
        elif name == 'CBEAM':
            self.cbeam.add(card_obj, comment=comment)
        elif name == 'CBEAMOR':
            self.cbeamor.add(card_obj, comment=comment)
        elif name == 'PBEAM':
            self.properties_beam.add_pbeam(card_obj, comment=comment)
        elif name == 'PBEAML':
            self.properties_beam.add_pbeaml(card_obj, comment=comment)
        elif name == 'PBCOMP':
            self.properties_beam.add_pbcomp(card_obj, comment=comment)

        # beam3
        # bend
        #========================
        # mass
        elif name == 'PMASS':
            self.mass.add_pmass(card_obj, comment=comment)
        elif name == 'CONM1':
            self.mass.add_conm1(card_obj, comment=comment)
        elif name == 'CONM2':
            self.mass.add_conm2(card_obj, comment=comment)
        elif name == 'CMASS1':
            self.mass.add_cmass1(card_obj, comment=comment)
        elif name == 'CMASS2':
            self.mass.add_cmass2(card_obj, comment=comment)
        elif name == 'CMASS3':
            self.mass.add_cmass3(card_obj, comment=comment)
        elif name == 'CMASS4':
            self.mass.add_cmass4(card_obj, comment=comment)
        elif name == 'CMASS5':
            self.mass.add_cmass5(card_obj, comment=comment)
        #========================

        # load combinations
        elif name == 'DLOAD':
            load = DLOAD(self)
            load.add_from_bdf(card_obj, comment=comment)
            self.loads.dload[load.load_id].append(load)
        elif name == 'LSEQ':
            #load = LSEQ(self)
            #load.add_from_bdf(card_obj, comment=comment)
            #self.loads.lseq[load.load_id].append(load)
            pass
        elif name == 'SLOAD':
            load = SLOAD(self)
            load.add_from_bdf(card_obj, comment=comment)
            self.loads.sload[load.load_id].append(load)
            #pass
        elif name == 'LOAD':
            #self.model.log.debug(card_obj)
            self.loads.load.add(card_obj, comment=comment)

        elif name == 'TEMPD':
            self.temps.add_tempd(card_obj, comment=comment)
        elif name == 'TEMP':
            self.temps.add_temp(card_obj, comment=comment)
        elif name == 'TEMPP1':
            self.temps.add_tempp1(card_obj, comment=comment)

        #========================
        # applied loads
        elif name == 'FORCE':
            self.loads.force.add(card_obj, comment=comment)
        elif name == 'FORCE1':
            self.loads.force1.add(card_obj, comment=comment)
        elif name == 'FORCE2':
            self.loads.force2.add(card_obj, comment=comment)
        elif name == 'MOMENT':
            self.loads.moment.add(card_obj, comment=comment)
        elif name == 'MOMENT1':
            self.loads.moment1.add(card_obj, comment=comment)
        elif name == 'MOMENT2':
            self.loads.moment2.add(card_obj, comment=comment)

        elif name == 'GRAV':
            self.loads.grav.add(card_obj, comment=comment)
        #========================
        # pressure loads
        elif name == 'PLOAD':
            self.loads.pload.add(card_obj, comment=comment)
        elif name == 'PLOAD1':
            self.loads.pload1.add(card_obj, comment=comment)
        elif name == 'PLOAD2':
            self.loads.pload2.add(card_obj, comment=comment)
        elif name == 'PLOAD4':
            self.loads.pload4.add(card_obj, comment=comment)
        elif name == 'PLOADX1':
            self.loads.ploadx1.add(card_obj, comment=comment)

        #========================
        # time loads
        elif name == 'TLOAD1':
            self.loads.tload1.add(card_obj, comment=comment)
        elif name == 'TLOAD2':
            self.loads.tload2.add(card_obj, comment=comment)

        # frequency loads
        #elif name == 'ACSRCE':
            #pass
        elif name == 'RLOAD1':
            self.loads.rload1.add(card_obj, comment=comment)
        elif name == 'RLOAD2':
            self.loads.rload2.add(card_obj, comment=comment)

        # other
        elif name == 'RFORCE':
            self.loads.rforce.add(card_obj, comment=comment)
        elif name == 'DAREA':
            self.loads.darea.add(card_obj, comment=comment)


        #========================
        # tables
        elif name == 'TABLED1':
            pass
        elif name == 'TABLED2':
            pass
        elif name == 'TABLED3':
            pass
        elif name == 'TABLED4':
            pass
        elif name == 'TABLED5':
            pass

        elif name == 'TABLEM1':
            pass
        elif name == 'TABLEM2':
            pass
        elif name == 'TABLEM3':
            pass
        elif name == 'TABLEM4':
            pass

        elif name == 'TABLES1':
            #self.materials.add_mats1(card_obj, comment=comment)
            pass
        elif name == 'TABLEST':
            pass
        #========================
        # mpc
        elif name == 'MPCADD':
            constraint_id, node_ids = get_spcadd_constraint(card_obj)
            mpcadd = self.mpcadd.setdefault(constraint_id, MPCADD(self))
            mpcadd.add(constraint_id, node_ids, comment=comment)
        elif name == 'MPC':
            constraint_id, constraint = get_mpc_constraint(card_obj)
            mpc = self.mpc.setdefault(constraint_id, MPC(self))
            mpc.add(constraint_id, constraint, comment=comment)
        #========================
        # spc
        elif name == 'SPCADD':
            constraint_id, node_ids = get_spcadd_constraint(card_obj)
            spcadd = self.spcadd.setdefault(constraint_id, SPCADD(self))
            spcadd.add(constraint_id, node_ids, comment=comment)
        elif name in ['SPC', 'SPCD']:
            for i in [0, 1]:
                constraint_id, node_id, dofs, enforced_motion = get_spc_constraint(card_obj, i)
                if enforced_motion != 0.0:
                    spcd = self.spcd.setdefault(constraint_id, SPCD(self))
                    spcd.add(constraint_id, node_id, dofs, enforced_motion, comment=comment)
                else:
                    spc = self.spc.setdefault(constraint_id, SPC(self))
                    spc.add(constraint_id, node_id, dofs, enforced_motion, comment=comment)
        elif name == 'SPC1':
            constraint_id, dofs, node_ids = get_spc1_constraint(card_obj)
            #print("type(spc1", type(self.spc1))
            spc1 = self.spc1.setdefault(constraint_id, SPC1(self))
            spc1.add(constraint_id, dofs, node_ids, comment=comment)
        #========================
        elif name == 'SUPORT':
            #self.suport.add(card_obj, comment=comment)
            pass
        elif name == 'SUPORT1':
            #self.suport1.add(card_obj, comment=comment)
            pass
        #========================
        # freq
        elif name == 'EIGB':
            card = EIGB.add_card(card_obj, comment=comment)
            #self.eigb[card.sid] = card
            self.methods[card.sid] = card
        elif name == 'EIGC':
            card = EIGC.add_card(card_obj, comment=comment)
            #self.eigc[card.sid] = card
            self.methods[card.sid] = card
        elif name == 'EIGR':
            card = EIGR.add_card(card_obj, comment=comment)
            #self.eigr[card.sid] = card
            self.methods[card.sid] = card
        elif name == 'EIGRL':
            card = EIGRL.add_card(card_obj, comment=comment)
            #self.eigrl[card.sid] = card
            self.methods[card.sid] = card

        #elif name == 'FREQ':
            #pass
        #elif name == 'FREQ1':
            #pass
        #elif name == 'FREQ2':
            #pass

        #========================
        # materials
        elif name == 'MAT1':
            self.materials.add_mat1(card_obj, comment=comment)
            #self.mat1.add(card_obj, comment=comment)
        elif name == 'MATS1':
            self.materials.add_mats1(card_obj, comment=comment)
        elif name == 'MAT4':
            self.materials.add_mat4(card_obj, comment=comment)
        elif name == 'MAT5':
            self.materials.add_mat5(card_obj, comment=comment)
        elif name == 'MAT8':
            self.materials.add_mat8(card_obj, comment=comment)
        elif name == 'MAT10':
            self.materials.add_mat10(card_obj, comment=comment)
        elif name == 'MATHP':
            self.materials.add_mathp(card_obj, comment=comment)
        elif name == 'MATHE':
            self.materials.add_mathe(card_obj, comment=comment)

        #========================
        elif name == 'RBAR':
            pass
        elif name == 'RBE2':
            self.elements.rbe2.add(card_obj, comment=comment)
        elif name == 'RBE3':
            self.elements.rbe3.add(card_obj, comment=comment)
        #========================
        elif name == 'CORD1R':
            self.coords.add_cord1r(card_obj, comment=comment)
        elif name == 'CORD1C':
            self.coords.add_cord1c(card_obj, comment=comment)
        elif name == 'CORD1S':
            self.coords.add_cord1s(card_obj, comment=comment)

        elif name == 'CORD2R':
            self.coords.add_cord2r(card_obj, comment=comment)
        elif name == 'CORD2C':
            self.coords.add_cord2c(card_obj, comment=comment)
        elif name == 'CORD2S':
            self.coords.add_cord2s(card_obj, comment=comment)

        #========================
        elif name == 'ACMODL':
            pass
        #========================
        # nodes
        elif name == 'GRID':
            self.grid.add(card_obj, comment=comment)
        elif name == 'GRDSET':
            self.grdset.add(card_obj, comment=comment)

        elif name == 'POINT':
            self.point.add(card_obj, comment=comment)
        elif name == 'SPOINT':
            self.spoint.add(card_obj, comment=comment)
        elif name == 'EPOINT':
            self.epoint.add(card_obj, comment=comment)
        elif name == 'POINTAX':
            self.pointax.add(card_obj, comment=comment)
        #elif name == 'RINGAX':
            #self.ringax.add(card_obj, comment=comment)
        #========================
        # nonlinear
        elif name == 'TSTEP':
            card = TSTEP.add_card(card_obj, comment=comment)
            self.tstep[tid] = card
        elif name == 'TSTEPNL':
            card = TSTEPNL.add_card(card_obj, comment=comment)
            self.tstep[tid] = card
        elif name == 'NLPARM':
            card = NLPARM.add_card(card_obj, comment=comment)
            self.nlparm[card.nlparm_id] = card
        elif name == 'NLPCI':
            card = NLPCI.add_card(card_obj, comment=comment)
            self.nlpci[card.nlpci_id] = card
        #========================
        else:
            raise NotImplementedError(name)
        return

    def add_AERO(self, card_obj, comment=''):
        assert self.aero is None, self.aero
        aero = AERO.add_card(card_obj, comment=comment)
        #assert key not in self.aero, '\naero=\n%s oldAERO=\n%s' % (
            #aero, self.aero[key])
        self.aero = aero

    def add_AEROS(self, card_obj, comment=''):
        assert self.aeros is None, self.aero
        aero = AEROS.add_card(card_obj, comment=comment)
        #assert key not in self.aeros, '\naeros=\n%s oldAEROS=\n%s' % (
            #aero, self.aeros[key])
        self.aeros = aero

    def add_AEFACT(self, card_obj, comment='', allow_overwrites=False):
        aefact = AEFACT.add_card(card_obj, comment=comment)
        key = aefact.sid
        if key in self.aefacts and not allow_overwrites:
            if not aefact._is_same_card(self.aefacts[key]):
                assert key not in self.aefacts, 'sid=%s\noldAEFACT=\n%snewAEFACT=\n%s' % (key, self.aefacts[key], aefact)
        else:
            assert key > 0, 'sid=%s method=\n%s' % (key, aefact)
            self.aefacts[key] = aefact

    def add_AELIST(self, card_obj, comment=''):
        aelist = AELIST(card_obj, comment=comment)
        key = aelist.sid
        assert key not in self.aelists, '\naelist=\n%s oldAELIST=\n%s' % (
            aelist, self.aelists[key])
        assert key >= 0
        self.aelists[key] = aelist

    def add_AELINK(self, card_obj, comment=''):
        aelink = AELINK.add_card(card_obj, comment=comment)
        key = aelink.id
        assert key >= 0
        if key not in self.aelinks:
            self.aelinks[key] = []
        self.aelinks[key].append(aelink)
        #assert key not in self.aestats,'\naestat=%s oldAESTAT=\n%s' %(aelink, self.aelinks[key])

    def add_AEPARM(self, card_obj, comment=''):
        aeparam = AEPARM.add_card(card_obj, comment=comment)
        key = aeparam.id
        assert key not in self.aeparams, '\naeparam=\n%s oldAESTAT=\n%s' % (
            aeparam, self.aeparams[key])
        assert key >= 0
        self.aeparams[key] = aeparam

    def add_AESTAT(self, card_obj, comment=''):
        aestat = AESTAT.add_card(card_obj, comment=comment)
        key = aestat.id
        assert key not in self.aestats, '\naestat=\n%s oldAESTAT=\n%s' % (
            aestat, self.aestats[key])
        assert key >= 0
        self.aestats[key] = aestat

    def add_AESURF(self, card_obj, comment=''):
        aesurf = AESURF.add_card(card_obj, comment=comment)
        key = aesurf.aesid
        assert key not in self.aesurfs, '\naesurf=\n%s oldAESURF=\n%s' % (
            aesurf, self.aesurfs[key])
        assert key >= 0
        self.aesurfs[key] = aesurf

    def get_bdf_stats(self, return_type='string'):
        """
        Print statistics for the BDF

        Parameters
        ----------
        return_type : str; default='string'
            'string'/'list':
                'string' returns one big string, while 'list' of each line in the string.
                'list' useful for the GUI.

        .. note:: if a card is not supported and not added to the proper
                  lists, this method will fail
        """
        msg = ['---BDF Statistics---']
        # sol
        msg.append('SOL %s\n' % self.sol)
        msg += self.grid.get_stats()

        msg += self.crod.get_stats()
        msg += self.prod.get_stats()
        msg += self.conrod.get_stats()
        msg += self.ctube.get_stats()
        msg += self.ptube.get_stats()

        msg += self.cbar.get_stats()
        msg += self.properties_bar.get_stats()
        #msg += self.pbar.get_stats()
        #msg += self.pbarl.get_stats()

        msg += self.cbeam.get_stats()
        msg += self.properties_beam.get_stats()
        #msg += self.pbeam.get_stats()
        #msg += self.pbeaml.get_stats()

        #msg += self.cbeam3.get_stats()
        #msg += self.pbeam3.get_stats()
        #msg += self.cbend.get_stats()
        #msg += self.pbend.get_stats()

        msg += self.elements_shell.get_stats()
        msg += self.properties_shell.get_stats()

        msg += self.elements_solid.get_stats()
        msg += self.properties_solid.get_stats()

        msg += self.materials.get_stats()

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

    def get_SPCx_ids(self, exclude_spcadd=False):
        """
        Get the SPC/SPCADD/SPC1/SPCAX IDs.

        Parameters
        ----------
        exclude_spcadd : bool; default=False
            you can exclude SPCADD if you just want a list of all the
            SPCs in the model.  For example, apply all the SPCs when
            there is no SPC=N in the case control deck, but you don't
            need to apply SPCADD=N twice.
        """
        spcs = {
            'SPC': self.spc,
            'SPC1': self.spc1,
            #'SPCAX': self.spcax,
            'SPCD': self.spcd,
        }
        if not exclude_spcadd:
            spcs['SPCADD'] = self.spcadd

        spc_ids = []
        for spc_type, spc in iteritems(spcs):
            spc_ids.extend(spc.keys())
        return np.unique(spc_ids)

    def SPC(self, spc_id, resolve=True, used_ids=None):
        """
        Gets all the MPCs that are in:
          - SPCADD
          - SPC
          - SPC1
          - SPCD

        Doesn't get:
          - SPCAX

        Parameters
        ----------
        spc_id : int
            the ID to get
        resolve : bool
            removes the SPCADDs by turning them into the
            cards they point to
        used_ids : list; default=None -> []
            an internal parameter
        """
        if used_ids is None:
            used_ids = []
        spc_out = []
        used_ids.append(spc_id)
        spcs = {
            'SPCADD' : self.spcadd,
            'SPC' : self.spc,
            'SPC1' : self.spc1,
            #'SPCAX' : self.spcax,
            'SPCD' : self.spcd,
        }

        if not resolve:
            for spc_type, spc in iteritems(spcs):
                if spc_id in spc:
                    spc_out.append(out)
            return spc_out

        for spc_type, spc in iteritems(spcs):
            if spc_id in spc:
                out = spc[spc_id]
                if spc_type == 'SPCADD':
                    for spci in out.spc_ids:
                        if spci in used_ids:
                            raise RuntimeError('duplicate SPC id=%i' % spci)
                        spc_outi = self.SPC(spci, resolve, used_ids)
                        for spcii in spc_outi:
                            spc_out.append(spcii)
                else:
                    spc_out.append(out)
        return spc_out

    def MPC(self, mpc_id, resolve=True, used_ids=None):
        """
        Gets all the MPCs that are in:
          - MPCADD
          - MPC

        Doesn't get:
          - MPCAX

        Parameters
        ----------
        mpc_id : int
            the ID to get
        resolve : bool
            removes the MPCADDs by turning them into the cards they
            point to
        used_ids : list; default=None -> []
            an internal parameter;
        """
        if used_ids is None:
            used_ids = []
        mpc_out = []
        mpcs = {
            'MPCADD' : self.mpcadd,
            'MPC' : self.mpc,
            #'MPCAX'  : self.mpcax,
        }
        used_ids.append(mpc_id)

        if not resolve:
            for mpc_type, mpc in iteritems(mpcs):
                if mpc_id in mpc:
                    mpc_out.append(out)
            return mpc_out

        for mpc_type, mpc in iteritems(mpcs):
            if mpc_id in mpc:
                out = mpc[mpc_id]
                if mpc_type == 'MPCADD':
                    for mpci in out.mpc_ids:
                        if mpci in used_ids:
                            raise RuntimeError('duplicate MPC id=%i' % mpci)
                        mpc_outi = self.MPC(mpci, resolve, used_ids)
                        for mpcii in mpc_outi:
                            mpc_out.append(mpcii)
                elif mpc_type in 'MPC':
                    mpc_out.append(out)
                else:
                    mpc_out.append(out)
        return mpc_out

    def _get_mass_types(self):
        """
        Gets the list of element types that have mass
        """
        types = [
            # O-D
            self.mass,

            # 1-D
            self.cbar, self.conrod, self.crod, self.ctube,
            self.cbeam, self.cbeam3, self.cbend,

            # 2-D
            self.elements_shell,

            # 3-D
            self.elements_solid,]
        return reduce_types(types)

    def _get_stiffness_types(self):
        """
        Gets the list of element types that have stiffness
        """
        types = [
            # O-D
            self.elements_spring,
            #celas1, self.celas2, self.celas3, self.celas4,

            # 1-D
            self.cbar, self.conrod, self.crod, self.ctube,
            self.cbeam, self.cbend, self.cbeam3,

            # 2-D
            self.elements_shell,

            # 3-D
            self.elements_solid,]
        return reduce_types(types)

    def _get_damping_types(self):
        """
        Gets the list of element types that have damping
        """
        types = [
            # O-D
            self.cdamp1, self.cdamp2, self.cdamp3, self.cdamp4, self.cdamp5,

            # 1-D
            self.cbar, self.conrod, self.crod, self.ctube,

            # 2-D
            self.elements_shell,

            # 3-D
            self.elements_solid,]
        return reduce_types(types)

    def mass_properties(self, total=False, sym_axis=None, scale=None):
        """
        .. todo:: consider calling BDF.elements.mass_properties
        """
        mass_types = self._get_mass_types()
        massi = []
        for mass_type in mass_types:
            massii = mass_type.get_mass_by_element_id(total=False)
            assert massii is not None, mass_type
            assert not isinstance(massii, float), mass_type
            #print("f massii =", massii)
            massi.extend(massii)

        massi = np.array(massi)
        total = True
        if total:
            mass = massi.sum()
        else:
            mass = massi
        return mass

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
            try:
                line = lines[i].rstrip('\r\n\t')
            except IndexError:
                break
            uline = line.upper()
            if uline.startswith('INCLUDE'):
                j = i + 1
                line_base = line.split('$')[0]
                include_lines = [line_base.strip()]
                # print('----------------------')

                line_base = line_base[8:].strip()
                if line_base.startswith("'") and line_base.endswith("'"):
                    pass
                else:
                    while not line.split('$')[0].endswith("'") and j < nlines:
                        # print('j=%s nlines=%s less?=%s'  % (j, nlines, j < nlines))
                        try:
                            line = lines[j].split('$')[0].strip()
                        except IndexError:
                            # print('bdf_filename=%r' % bdf_filename)
                            crash_name = 'pyNastran_crash.bdf'
                            self._dump_file(crash_name, lines, i+1)
                            msg = 'There was an invalid filename found while parsing (index).\n'
                            msg += 'Check the end of %r\n' % crash_name
                            msg += 'bdf_filename2 = %r' % bdf_filename
                            raise IndexError(msg)
                        # print('endswith_quote=%s; %r' % (line.split('$')[0].strip().endswith(""), line.strip()))
                        include_lines.append(line.strip())
                        j += 1
                    # print('j=%s nlines=%s less?=%s'  % (j, nlines, j < nlines))

                    #print('*** %s' % line)
                    #bdf_filename2 = line[7:].strip(" '")
                    #include_lines = [line] + lines[i+1:j]
                #print(include_lines)
                bdf_filename2 = get_include_filename(include_lines, include_dir=self.include_dir)

                try:
                    self._open_file_checks(bdf_filename2)
                except IOError:
                    crash_name = 'pyNastran_crash.bdf'
                    self._dump_file(crash_name, lines, j)
                    msg = 'There was an invalid filename found while parsing.\n'
                    msg += 'Check the end of %r\n' % crash_name
                    msg += 'bdf_filename2 = %r' % bdf_filename2
                    #msg += 'len(bdf_filename2) = %s' % len(bdf_filename2)
                    raise IOError(msg)

                with self._open_file(bdf_filename2, basename=False) as bdf_file:
                    #print('bdf_file.name = %s' % bdf_file.name)
                    lines2 = bdf_file.readlines()

                #print('lines2 = %s' % lines2)
                nlines += len(lines2)

                #line2 = lines[j].split('$')
                #if not line2[0].isalpha():
                    #print('** %s' % line2)

                include_comment = '\n$ INCLUDE processed:  %s\n' % bdf_filename2
                #for line in lines2:
                    #print("  ?%s" % line.rstrip())
                lines = lines[:i] + [include_comment] + lines2 + lines[j:]
                #for line in lines:
                    #print("  *%s" % line.rstrip())
            i += 1

        if self.dumplines:
            self._dump_file('pyNastran_dump.bdf', lines, i)
        return _lines_to_decks(lines, i, punch)

    def _open_file(self, bdf_filename, basename=False, check=True):
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

        self._validate_open_file(bdf_filename, bdf_filename_inc, check)


        self.log.debug('opening %r' % bdf_filename_inc)
        self.active_filenames.append(bdf_filename_inc)

        #print('ENCODING - _open_file=%r' % self._encoding)
        bdf_file = codec_open(_filename(bdf_filename_inc), 'r', encoding=self._encoding)
        return bdf_file

    def _validate_open_file(self, bdf_filename, bdf_filename_inc, check):
        """
        checks that the file doesn't have obvious errors
         - hasn't been used
         - not a directory
         - is a file

        Parameters
        ----------
        bdf_filename : str
           the current bdf filename
        bdf_filename_inc : str
           the next bdf filename

        Raises
        ------
        RuntimeError : file is active
        IOError : Invalid file type
        """
        if check:
            if not os.path.exists(_filename(bdf_filename_inc)):
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
            elif os.path.isdir(_filename(bdf_filename)):
                current_filename = self.active_filename if len(self.active_filenames) > 0 else 'None'
                raise IOError('Found a directory: bdf_filename=%r\ncurrent_file=%s' % (
                    bdf_filename_inc, current_filename))
            elif not os.path.isfile(_filename(bdf_filename)):
                raise IOError('Not a file: bdf_filename=%r' % bdf_filename)

    def _parse_cards(self, cards, card_count):
        """creates card objects and adds the parsed cards to the deck"""
        #print('card_count = %s' % card_count)
        #print(type(cards))
        if isinstance(cards, dict): # self._is_cards_dict = True
            self.allocate(card_count, cards)
            print('cards.keys() = ', cards.keys())
            for card_name, card in sorted(iteritems(cards)):
                if self.is_reject(card_name):
                    self.log.info('    rejecting card_name = %s' % card_name)
                    self.log.info('    rejecting card = %s' % card)
                    for cardi in card:
                        self._increase_card_count(card_name)
                        print(cardi[0])
                        print(cardi[1])
                        self.rejects.append([cardi[0]] + cardi[1])
                else:
                    self.log.info('    card_name = %s' % card_name)
                    #self.log.info('    card = %s' % card)
                    #self.log.info('    len(card) = %s' % len(card))
                    for comment, card_lines in card:
                        #self.log.info('    card_lines = %s' % card_lines)
                        self.add_card(card_lines, card_name, comment=comment,
                                      is_list=False, has_none=False)
                        #self._add_card_object(card_name, card_lines, comment=comment)
        else:
            asf
            for card in cards:
                card_name, comment, card_lines = card
                if card_name is None:
                    msg = 'card_name = %r\n' % card_name
                    msg += 'card_lines = %s' % card_lines
                    raise RuntimeError(msg)
                if self.is_reject(card_name):
                    if card_name not in self.card_count:
                        if ' ' in card_name:
                            msg = (
                                'No spaces allowed in card name %r.  '
                                'Should this be a comment?\n%s%s' % (
                                    card_name, comment, card_lines))
                            raise RuntimeError(msg)
                        if card_name in ['SUBCASE ', 'CEND']:
                            raise RuntimeError('No executive/case control deck was defined.')
                        self.log.info('    rejecting card_name = %s' % card_name)
                    self._increase_card_count(card_name)
                    self.rejects.append([comment] + card_lines)
                else:
                    self.add_card(card_lines, card_name, comment=comment,
                                  is_list=False, has_none=False)

    def _parse_dynamic_syntax(self, key):
        """
        Applies the dynamic syntax for %varName

        Parameters
        ----------
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

    def _parse_primary_file_header(self, bdf_filename):
        """
        Extract encoding, nastran_format, and punch from the primary BDF.

        Parameters
        ----------
        bdf_filename : str
            the input filename

        ..code-block :: python

            $ pyNastran: version=NX
            $ pyNastran: encoding=latin-1
            $ pyNastran: punch=True
            $ pyNastran: dumplines=True
            $ pyNastran: nnodes=10
            $ pyNastran: nelements=100
            $ pyNastran: skip_cards=PBEAM,CBEAM
            $ pyNastran: units=in,lb,s

        ..warning :: pyNastran lines must be at the top of the file
        """
        with open(bdf_filename, 'r') as bdf_file:
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
                        elif key == 'dumplines':
                            self.dumplines = True if value == 'true' else False
                        elif key == 'skip_cards':
                            cards = {value.strip() for value in value.upper().split(',')}
                            self.cards_to_read = self.cards_to_read - cards
                        elif 'skip ' in key:
                            type_to_skip = key[5:].strip()
                            #values = [int(value) for value in value.upper().split(',')]
                            values = parse_patran_syntax(value)
                            if type_to_skip not in self.object_attributes():
                                raise RuntimeError('%r is an invalid key' % type_to_skip)
                            if type_to_skip not in self.values_to_skip:
                                self.values_to_skip[type_to_skip] = values
                            else:
                                self.values_to_skip[type_to_skip] = np.hstack([
                                    self.values_to_skip[type_to_skip],
                                    values
                                ])
                        #elif key == 'skip_elements'
                        #elif key == 'skip_properties'
                        elif key == 'units':
                            self.units = [value.strip() for value in value.upper().split(',')]
                        else:
                            raise NotImplementedError(key)
                    else:
                        break
                else:
                    break

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
        #if 0:
            #for key, card in sorted(iteritems(self.elements)):
                #try:
                    #card._verify(xref)
                #except:
                    #print(str(card))
                    #raise
            #for key, card in sorted(iteritems(self.properties)):
                #try:
                    #card._verify(xref)
                #except:
                    #print(str(card))
                    #raise
        self.materials._verify(xref)


def _lines_to_decks(lines, i, punch):
    """
    Splits the lines into their deck.
    """
    executive_control_lines = []
    case_control_lines = []
    bulk_data_lines = []

    if punch:
        bulk_data_lines = lines
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

def reduce_types(types):
    types2 = []
    for etype in types:
        if etype is not None and etype.n:
            types2.append(etype)
    return types2

def _clean_comment(comment, end=-1):
    """
    Removes specific pyNastran comment lines so duplicate lines aren't
    created.

    Parameters
    ----------
    comment : str
        the comment to possibly remove
    end : int; default=-1
        currently leftover from the unvectorized version

    Returns
    -------
    updated_comment : str
        the comment
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
                         '$SETS', '$CONTACT', '$REJECTS', '$REJECT_LINES']:
        comment = ''
    return comment

#def print_filename(filename, relpath):
    #"""
    #Takes a path such as C:/work/fem.bdf and locates the file using
    #relative paths.  If it's on another drive, the path is not modified.

    #:param filename: a filename string
    #:returns filename_string: a shortened representation of the filename
    #"""
    #driveLetter = os.path.splitdrive(os.path.abspath(filename))[0]
    #if driveLetter == os.path.splitdrive(os.curdir)[0] and self._relpath:
        #return os.path.relpath(filename)
    #return filename



if __name__ == '__main__':  # pragma: no cover
    bdf = BDF()
    pkg_path = pyNastran.__path__[0]
    bdfname = sys.argv[1]
    #print("bdfname =", bdfname)
    bdf.read_bdf(bdfname)
    bdf.write_bdf('fem.out.bdf')
