"""
Main BDF class
"""
#from __future__ import (nested_scopes, generators, division, absolute_import,
#                        unicode_literals)

from __future__ import (nested_scopes, generators, division, absolute_import,
                        print_function, unicode_literals)

import io
import os
import sys
import warnings
import traceback
from collections import defaultdict

from numpy import unique, array

from pyNastran.bdf.bdf import (to_fields, get_include_filename,
    parse_executive_control_deck, clean_empty_lines)
from pyNastran.bdf.bdf import EIGRL

from pyNastran.utils import list_print, is_string, object_attributes
from pyNastran.utils.log import get_logger
from pyNastran.utils.gui_io import load_file_dialog


from .cards.coord import Coord

# nodes
from .cards.nodes.grid import GRID, GRDSET
from .cards.nodes.point import POINT
from .cards.nodes.spoint import SPOINT
from .cards.nodes.epoint import EPOINT
from .cards.nodes.pointax import POINTAX


# spring
from .cards.elements.spring.elements_spring import ElementsSpring
from .cards.elements.spring.pelas import PELAS

# shell
from .cards.elements.shell.elements_shell import ElementsShell
from .cards.elements.shell.properties_shell import PropertiesShell

# solid
from .cards.elements.solid.elements_solid import ElementsSolid
from .cards.elements.solid.properties_solid import PropertiesSolid

# rods
from .cards.elements.rod.prod import PROD
from .cards.elements.rod.crod import CROD
from .cards.elements.rod.conrod import CONROD

# shear
from .cards.elements.shear.cshear import CSHEAR
from .cards.elements.shear.pshear import PSHEAR

# bar
from .cards.elements.bar.cbar import CBAR #, CBAROR
#from .cards.elements.bar.pbar import PBAR
#from .cards.elements.bar.pbarl import PBARL
from .cards.elements.bar.properties_bar import PropertiesBar

# beams

# mass
from .cards.elements.mass.mass import Mass

#===========================
# aero

from .cards.aero.caero import CAero
from .cards.aero.paero import PAero
from .cards.aero.trim import TRIM
#from .cards.aero.aero import AERO
#from .cards.aero.aeros import AEROS

#===========================

# materials
from .cards.materials.materials import Materials


# loads
from .cards.loads.load import LOAD
from .cards.loads.dload import DLOAD
from .cards.loads.dload import DLOAD as LSEQ
from .cards.loads.dload import DLOAD as SLOAD
from .cards.loads.grav import GRAV
from .cards.loads.force import FORCE
from .cards.loads.moment import MOMENT

#from .cards.loads.force1 import FORCE1
#from .cards.loads.force2 import FORCE2
#from .cards.loads.moment1 import MOMENT1
#from .cards.loads.moment2 import MOMENT2

# ACCEL1
# PLOAD3
# DAREA
# TLOAD1
# TLOAD2
# RLOAD1
# RLOAD2
# RANDPS


from .cards.loads.pload  import PLOAD
from .cards.loads.pload1 import PLOAD1
from .cards.loads.pload2 import PLOAD2
#from .cards.loads.pload3 import PLOAD3
#from .cards.loads.pload4 import PLOAD4

from .cards.loads.ploadx1 import PLOADX1
#from .cards.loads.grav import GRAV

from .cards.loads.rforce import RFORCE
#from .cards.loads.sload import SLOAD

#from .cards.loads.loadcase import LoadCase
#from .cards.loads.loadset import LOADSET


#=============================
# dynamic
from .cards.nonlinear.nlpci import NLPCI
from .cards.nonlinear.nlparm import NLPARM

#=============================

# constraints
from .cards.constraints.spc import SPC, get_spc_constraint
from .cards.constraints.spcd import SPCD

from .cards.constraints.spc1 import SPC1, get_spc1_constraint
from .cards.constraints.spcadd import SPCADD, get_spcadd_constraint
from .cards.constraints.mpcadd import MPCADD

#from .cards.elements.elements import CFAST, CGAP, CRAC2D, CRAC3D
#from .cards.properties.properties import (PFAST, PGAP, PLSOLID, PSOLID,
#                                          PRAC2D, PRAC3D, PCONEAX)

#from .cards.elements.springs import (CELAS3, CELAS4)
#from .cards.properties.springs import PELAS, PELAST

#from .cards.elements.solid import (CTETRA4, CTETRA10, CPENTA6, CPENTA15,
#                                   CHEXA8, CHEXA20, SolidElement)
#from .cards.elements.rigid import (RBAR, RBAR1, RBE1, RBE2, RBE3, RigidElement)
#
#from .cards.elements.shell import (CQUAD, CQUADR, CQUADX,
#                                   CTRIAX, CTRIAR)
#from .cards.properties.shell import PSHELL, PCOMP, PCOMPG, PSHEAR, PLPLANE
#from .cards.elements.bush import CBUSH, CBUSH1D, CBUSH2D
#from .cards.properties.bush import PBUSH, PBUSH1D
#from .cards.elements.damper import (CVISC, CDAMP1, CDAMP2, CDAMP3, CDAMP4,
#                                    CDAMP5, DamperElement)
#from .cards.properties.damper import (PVISC, PDAMP, PDAMP5, PDAMPT)
#from .cards.elements.bars import (CTUBE, CBEAM, CBEAM3,
#                                  CBEND, LineElement, RodElement)
#from .cards.properties.bars import (PTUBE, PBAR, PBARL,
#                                    PBEAM, PBEAML, PBCOMP)  # PBEND
#from .cards.properties.mass import (PMASS, NSM)
#from .cards.aero import (AEFACT, AELINK, AELIST, AEPARM, AESTAT, AESURF,
#                         AESURFS, AERO, AEROS, CSSCHD, CAERO1, CAERO2, CAERO3,
#                         CAERO4, CAERO5, FLFACT, FLUTTER, GUST, MKAERO1,
#                         MKAERO2, PAERO1, PAERO2, SPLINE1, SPLINE2, SPLINE4,
#                         SPLINE5)
#from .cards.constraints import (SPC, SPCADD, SPCD, SPCAX, SPC1,
#                                MPC, MPCADD, SUPORT1, SUPORT,
#                                ConstraintObject)
from .cards.coordinateSystems import (CORD1R, CORD1C, CORD1S,
                                      CORD2R, CORD2C, CORD2S, CORD3G)
#from .cards.dmig import (DEQATN, DMIG, DMI, DMIJ, DMIK, DMIJI, NastranMatrix)
#from .cards.dynamic import (FREQ, FREQ1, FREQ2, FREQ4, TSTEP, TSTEPNL, NLPARM, NLPCI)
#from .cards.loads.loads import (LSEQ, SLOAD, DLOAD, DAREA, TLOAD1, TLOAD2,
#                                RLOAD1, RLOAD2, RANDPS)
#from .cards.loads.staticLoads import (LOAD, ACCEL1,
#                                      FORCE1, FORCE2, MOMENT1, MOMENT2,
#                                      PLOAD4)
#
#from .cards.materials import (MAT2, MAT3, MAT4, MAT5,
#                              MAT8, MAT9, MAT10,
#                              MATHP, CREEP, EQUIV)
#from .cards.methods import (EIGB, EIGC, EIGR, EIGP, EIGRL)
#from .cards.optimization import (DCONSTR, DESVAR, DDVAL, DOPTPRM, DLINK,
#                                 DRESP1, DRESP2, DVMREL1, DVPREL1, DVPREL2)
from .cards.params import PARAM
#from .cards.sets import (ASET, BSET, CSET, QSET,
#                         ASET1, BSET1, CSET1, QSET1,
#                         SET1, SET3, SESET, SEQSEP, RADSET)
#from .cards.thermal.loads import (QBDY1, QBDY2, QBDY3, QHBDY, TEMP, TEMPD)
#from .cards.thermal.thermal import (CHBDYE, CHBDYG, CHBDYP, PCONV, PCONVM,
#                                    PHBDY, CONV, RADM, RADBC,)
#from .cards.tables import (TABLED1, TABLED2, TABLED3,
#                           TABLEM1, TABLEM2, TABLEM3, TABLEM4,
#                           TABLES1, TABLEST, TABRND1, TABRNDG, TIC)
#from .cards.contact import BCRPARA, BCTADD, BCTSET, BSURF, BSURFS
from pyNastran.bdf.caseControlDeck import CaseControlDeck
from .bdf_methods import BDFMethods
from .bdf_interface.get_methods import GetMethods
from .bdf_interface.add_card import AddCard
from .bdf_interface.utils import wipe_empty_fields
from .bdf_interface.assign_type import interpret_value
from .bdf_interface.write_mesh import WriteMesh
from .bdf_interface.cross_reference import XRefMesh

# old
from pyNastran.bdf.bdfInterface.BDF_Card import BDFCard


def class_obj_defaultdict(class_obj, *args, **kwargs):

    class ClassObjDefaultDict(dict):
        def __init__(self):
            dict.__init__(self)

        def __getitem__(self, key):
            #print('getting key=%r' % key)
            #print('dir(d) =', dir(self))
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
    modelType = 'nastran'
    #: Flips between a dictionary based storage BDF storage method and
    #: a list based method.  Don't modify this.
    _isDict = True

    #: required for sphinx bug
    #: http://stackoverflow.com/questions/11208997/autoclass-and-instance-attributes
    #__slots__ = ['_is_dynamic_syntax']
    def __init__(self, debug=True, log=None, precision='double'):
        """
        Initializes the BDF object
        :param self:  the BDF object
        :param debug: used to set the logger if no logger is passed in
        :param log:   a python logging module object
        """
        assert debug in [True, False], 'debug=%r' % debug
        
        if precision == 'double':
            self.float = 'float64'
        elif precision == 'single':
            self.float = 'float32'
        else:
            raise NotImplementedError('precision=%r' % precision)

        # file management parameters
        self._ifile = -1
        self.include_dir = ''
        self.active_filename = None
        self.active_filenames = []
        #self.used_filenames = []
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
        #: lines that were rejected b/c they were for a card that isnt supported
        self.rejects = []
        #: cards that were created, but not processed
        self.reject_cards = []
        #: list of execive control deck lines
        self.executive_control_lines = []
        #: list of case control deck lines
        self.case_control_lines = []
        
        self.echo = False

        self.__init_attributes()

        #: the list of possible cards that will be parsed
        self.cards_to_read = set([
            'PARAM',
            'GRID', 'GRDSET', 'SPOINT', 'EPOINT', 'POINT', 'POINTAX', # 'RINGAX',

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
            'MAT1', 'MAT2', 'MAT3', 'MAT8', 'MAT9', 'MAT10', 'MATHP',
            #'MATT1', 'MATT2', 'MATT3', 'MATT4', 'MATT5', 'MATT8', 'MATT9',
            'MATS1',# 'MATHE'
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
            'GRAV', 'ACCEL1',
            'PLOAD', 'PLOAD1', 'PLOAD2', 'PLOAD4',
            'PLOADX1', 'RFORCE',

            # aero cards
            'AERO', 'AEROS', 'GUST', 'FLUTTER', 'FLFACT', 'MKAERO1', 'MKAERO2',
            'AEFACT', 'AELINK', 'AELIST', 'AEPARAM', 'AESTAT', 'AESURF',
            'CAERO1', 'CAERO2', 'CAERO4', # 'CAERO3', 'CAERO5',
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
            'DAREA', 'NLPARM', 'NLPCI', 'TSTEP', 'TSTEPNL',

            # frequencies
            'FREQ', 'FREQ1', 'FREQ2',

            # direct matrix input cards
            'DMIG', 'DMIJ', 'DMIJI', 'DMIK', #'DMI',
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

            #: methods - .. todo:: EIGRL not done???
            'EIGB', 'EIGR', 'EIGRL',

            #: cMethods - .. todo:: EIGC not done???
            'EIGC', 'EIGP',

            #: contact
            'BCRPARA', 'BCTADD', 'BCTSET', 'BSURF', 'BSURFS',

            # other
            'INCLUDE',  # '='
            'ENDDATA',
        ])

        caseControlCards = set(['FREQ', 'GUST', 'MPC', 'SPC', 'NLPARM', 'NSM',
                                'TEMP', 'TSTEPNL', 'INCLUDE'])
        self.uniqueBulkDataCards = self.cards_to_read.difference(
            caseControlCards)

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
        
        model = self
        self.grid = GRID(model)
        self.point = POINT(model)
        self.grdset = GRDSET(model)
        self.spoint = SPOINT(model)
        self.epoint = EPOINT(model)
        self.pointax = POINTAX(model)
        
        self.coords = {0 : CORD2R() }
        
        #: stores elements (CQUAD4, CTRIA3, CHEXA8, CTETRA4, CROD, CONROD,
        #: etc.)
        self.elements = {}
        #: stores rigid elements (RBE2, RBE3, RJOINT, etc.)
        self.rigidElements = {}
        #: stores LOTS of propeties (PBAR, PBEAM, PSHELL, PCOMP, etc.)
        self.properties = {}
        
        #self.properties_spring = PropertiesSpring(model)
        #self.proeprties_rod = PropertiesRod(v)
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
        
        # mass
        #: stores CONM1, CONM2, CMASS1, CMASS2, CMASS3, CMASS4, CMASS5, PMASS
        self.mass = Mass(model)
        
        # bars
        #: stores CBAR
        self.cbar = CBAR(model)
        #: stores PBAR, PBARL
        self.properties_bar = PropertiesBar(model)

        # beams

        # solids
        #: stores CTETRA4, CPENTA6, CHEXA8, CTETRA10, CPENTA15, CHEXA20
        self.elements_solid = ElementsSolid(self)
        #: stores PSOLID, PLSOLID
        self.properties_solid = PropertiesSolid(self)

        #===================================
        # methods
        self.eigrl = {}
        self.eigb = {}
        self.eigc = {}
        #===================================
        # aero

        #: stores CAERO1, CAERO2, CAERO3, CAERO4, CAERO5
        self.caero = CAero(model)
        #: stores PAERO1, PAERO2, PAERO3, PAERO4, PAERO5
        self.paero = PAero(model)

        #: stores TRIM
        self.trim = {}
        #self.aero = AERO(model)
        #self.aeros = AEROS(model)
        #===================================
        # optimization
        #self.dconstr = DCONSTR(model)
        
        #===================================
        # loads
        #self.loadcase = LoadCase(model)

        self.load = defaultdict(list)
        self.dload = defaultdict(list)
        #self.loadset = LOADSET(model)
        
        self.force = FORCE(model)
        #self.force1 = FORCE1(model)
        #self.force2 = FORCE2(model)
        self.moment = MOMENT(model)
        #self.moment1 = MOMENT1(model)
        #self.moment2 = MOMENT2(model)
        self.grav = GRAV(model)
        self.rforce = RFORCE(model)
        
        self.pload = PLOAD(model)
        self.pload1 = PLOAD1(model)
        self.pload2 = PLOAD2(model)
        #self.pload3 = PLOAD3(model)
        #self.pload4 = PLOAD4(model)
        
        self.ploadx1 = PLOADX1(model)

        #: stores MAT1, MAT2, MAT3,...MAT10 (no MAT4, MAT5)
        #self.materials = {}
        self.materials = Materials(model)
        
        #: stores MATS1
        self.materialDeps = {}
        #: stores the CREEP card
        self.creepMaterials = {}

        # loads
        #: stores LOAD, FORCE, MOMENT, etc.
        #self.loads = {}
        #self.gusts  = {} # Case Control GUST = 100
        #self.random = {} # Case Control RANDOM = 100

        #: stores coordinate systems
        self.coord = Coord(model)

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

        # ------------------------- thermal defaults -------------------------
        # BCs
        #: stores thermal boundary conditions - CONV,RADBC
        self.bcs = {}  # e.g. RADBC
        #: defines the MAT4, MAT5, MATT4, etc.  .. todo:: verify MATT4
        self.thermalMaterials = {}

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
        if 0:
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
        self.materials._verify(xref)

    def read_bdf(self, bdf_filename=None, include_dir=None, xref=True, punch=False):
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
        >>> print g1.Position()
        [10.0, 12.0, 42.0]
        >>> bdf.write_bdf(bdf_filename2)
        >>> print bdf.card_stats()
        ---BDF Statistics---
        SOL 101
        bdf.nodes = 20
        bdf.elements = 10
        etc.
        """
        if bdf_filename is None:
            wildcard_wx = "Nastran BDF (*.bdf; *.dat; *.nas; *.pch)|*.bdf;*.dat;*.nas;*.pch|" \
                "All files (*.*)|*.*"
            wildcard_qt = "Nastran BDF (*.bdf *.dat *.nas *.pch);;All files (*)"
            title = 'Please select a BDF/DAT/PCH to load'
            bdf_filename = load_file_dialog(title, wildcard_wx, wildcard_qt)

        #: the active filename (string)
        self.bdf_filename = bdf_filename
        if include_dir is None:
            include_dir = os.path.dirname(bdf_filename)

        #: the directory of the 1st BDF (include BDFs are relative to this one)
        self.include_dir = include_dir

        #: will this model be cross referenced
        self._xref = xref

        if not os.path.exists(bdf_filename):
            raise IOError('cannot find bdf_filename=%r' % bdf_filename)
        if bdf_filename.lower().endswith('.pch'):
            punch = True
        
        #: is this a punch file (no executive control deck)
        self._punch = punch
        try:
            self._open_file(self.bdf_filename)
            self.log.debug('---starting BDF.read_bdf of %s---' % self.bdf_filename)
            if not punch:
                self.log.debug('---reading executive and case control decks---')
                self._read_executive_control_deck()
                self._read_case_control_deck()
            else:
                self.log.debug('---skipping executive and case control decks---')

            self._read_bulk_data_deck()
            self.cross_reference(xref=xref)
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
        while 'CEND' not in lineUpper[:4] and 'BEGIN' not in lineUpper and 'BULK' not in lineUpper:
            #lines = []
            try:
                (i, line, comment) = self._get_line()
            except TypeError:
                msg = 'Failed getting line.  If this file does not contain an executive control deck, \n'
                if self.include_dir == '':
                    include_dir = ''
                else:
                    include_dir = 'include_dir=%s, ' % self.include_dir
                if self._xref:
                    xref = ''
                else:
                    xref = 'xref=False, '
                msg += 'call read_bdf(bdf_filename=%r, %s%spunch=%s)' % (self.bdf_filename, include_dir, xref, True)
                msg += ' instead.\n'
                raise RuntimeError(msg)
            line = line.rstrip('\n\r\t ')
            #print("line exec = %r" % line)

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
        #self.sol = sol
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
            assert isinstance(key, str), 'key=%s must be a string' % key
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
        #print "*** valueRaw.lstrip() = |%r|" % valueRaw.lstrip()
        #key = key.lstrip('%%')
        key = key[1:].strip()
        self.log.info("dynamic key = |%r|" % key)
        #self.dict_of_vars = {'P5':0.5,'ONEK':1000.}
        if key not in self.dict_of_vars:
            msg = "key=|%r| not found in keys=%s" % (key, self.dict_of_vars.keys())
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
        #print("reading Case Control Deck...")
        if not self.has_case_control_deck:
            return
        line = ''
        while self.active_filename:  # keep going until finished
            #print "top of loop"
            lines = []
            (i, lineIn, comment) = self._get_line()
            if lineIn is None:
                return  # file was closed
            line = lineIn.strip().split('$')[0].strip()
            lineUpper = line.upper()
            #print("lineUpper = %r" % str(lineUpper))
            if lineUpper.startswith('INCLUDE'):
                #print("INCLUDE!!!")
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
                    (i, line_next, comment) = self._get_line()
                    next_line = next_line.strip().split('$')[0].strip()
                self.case_control_lines.append(next_line)
                #except StopIteration:
                    #include_lines = [line]
                filename = get_include_filename(include_lines,
                                                 include_dir=self.include_dir)
                self._open_file(filename)
                #line = next_line
            else:
                self.case_control_lines.append(lineUpper)

            if 'BEGIN' in lineUpper and (('BULK' in lineUpper) or ('SUPER' in lineUpper)):
                self.log.debug('found the end of the Case Control Deck!')
                break
        self.log.info("finished with Case Control Deck...")

        #for line in self.case_control_lines:
            #print "** line=|%r|" %(line)

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
        elif card_name in self.cards_to_read:
            return False
        if card_name:
            if card_name not in self.reject_count:
                self.log.info("reject card_name = |%s|" % card_name)
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
        bdf_filename = os.path.join(self.include_dir, str(bdf_filename))
        if not os.path.exists(bdf_filename):
            msg = 'No such bdf_filename: %r\n' % bdf_filename
            msg += 'cwd: %r' % os.getcwd()
            raise IOError(msg)

        #print "opening self.active_filename=%s" % bdf_filename
        if bdf_filename in self.active_filenames:
            msg = 'bdf_filename=%s is already active.\nactive_filenames=%s' % (bdf_filename, self.active_filenames)
            raise RuntimeError(msg)
        self.log.info('opening %r' % bdf_filename)

        self._ifile += 1
        self.active_filename = bdf_filename
        #self.used_filenames.append(bdf_filename)

        self.active_filenames.append(bdf_filename)

        self._stored_Is[self._ifile] = []
        self._stored_lines[self._ifile] = []
        self._stored_comments[self._ifile] = []

        line_gen = self._stream_line()
        if self._isDict:
            #print "making ifile=%s" % self._ifile
            self._line_streams[self._ifile] = line_gen
            self._card_streams[self._ifile] = self._stream_card(line_gen)
        else:
            self._line_streams.append(line_gen)
            self._card_streams.append(self._stream_card(line_gen))


    def _close_file(self):
        self.log.info('closing %r' % self.active_filename)
        if self._ifile == 0:
            self._ifile = -1
            self.active_filename = None
            return
        if self._isDict:
            #print "killing ifile=%s" % self._ifile
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
        #print "active_filenames =", self.active_filenames
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
        #while 1:
            #-----------------------------------------------------------------
            # get the first line of the card
            #self.log.debug("stored_lines = %s" % self._stored_lines)
            #self.log.debug("new card...line =%r" % line)
            Is = []
            lines = []
            comments = []

            #c = ''
            comment = _clean_comment(comment)
            if comment:
                #c=' comment=|%s|' % comment.strip()
                comments.append(comment)
            #self.log.debug("lineA1 %s line=|%s|%s" % (i, line.strip(), c))

            # If the first line is valid, continue.
            # Otherwise, keep getting lines until one isn't blank.
            if line:
                #self.log.debug("adding line...")
                Is.append(i)
                lines.append(line)
            else:
                #self.log.debug("while block...")
                while len(line)==0:
                    #print "lines while len(line)==0; ",lines
                    #print "you cant have an empty first line..."
                    (i, line, comment) = self._get_line()
                    if line:
                        break

                    #c = ''
                    comment = _clean_comment(comment)
                    if comment:
                        #c=' comment=|%s|' % comment.strip()
                        comments.append(comment)
                    #print "lineA2 %s line=|%s|%s" % (i, line.strip(), c)
                Is.append(i)
                lines.append(line)
                if comment:
                    comments.append(comment)
            assert len(lines) == 1, lines
            #if 'PBARL' in lines[0]:
            #print "===end of block 1 (get first line)===; lines=%r" % lines
            #pass

            #print "lines =",lines

            #-----------------------------------------------------------------
            # get another line
            #self.log.debug("*lineC %s line=|%s|%s" % (i, line.strip(), c))
            try:
                (i, line, comment) = self._get_line()
            except TypeError:
                #print 'type yield...'
                #raise
                lines2 = clean_empty_lines(lines)
                yield lines2, ''.join(comments)
            if comment:  c=' comment=|%s|' % comment.strip()
            #self.log.debug("lineC %s line=|%s|%s" % (i, line.strip(), c))
            #print "lines =",lines
            #print ""

            #-----------------------------------------------------------------
            # We define a continuation by either a regular,
            # large field, small field, tab, or CSV formatted line.
            # Large field - a * is in the first character
            # Small field - a + or ' ' is in the first character
            #               or the line is blank
            # Tab - tab separated value; large or small formatted line
            # CSV - comma separated value; large or small formatted line

            # If the line is a continuation line, keep going.
            in_loop = False
            #print "checkline = %r" % line

            Is2 = []
            lines2 = []
            comments2 = []
            while len(line)==0 or line[0] in [' ', '*', '+', ',', '\t']:
                in_loop = True
                #print "into the loop!"
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

                if comment:  c=' comment=|%s|' % comment.strip()
                #self.log.debug("lineD %s line=|%s|%s" % (i, line.strip(), c))
                try:
                    (i, line, comment) = self._get_line()
                except TypeError:
                    #raise
                    lines2 = clean_empty_lines(lines)
                    comment = ''.join(comments+comments2)
                    #self.log.debug("*yieldingTYPE lines2=%r comments=%s" %(lines2, comments))
                    yield lines2, comment

                #print "len(line)=%s line[0]=%r" % (len(line), line[0])
            #if not in_loop:
                #print "len(line)=%i line=%r" % (len(line), line)

            # the extra lines we grabbed in the while loop should go on the
            # next card
            if Is2:
                #self.log.debug("storing lines2 %s" % lines2)
                self._stored_Is[self._ifile] = Is2
                self._stored_lines[self._ifile] = lines2
                self._stored_comments[self._ifile] = comments2
            #Is2 = []
            #lines2 = []
            #comments2 = []

            c = ''
            if comment:  c =' comment=|%s|' % comment.strip()
            #print "lineD2 %s line=|%s|%s" % (i, line.strip(), c)
            #if 'PBARL' in lines[0]:
            #print '===end of block 2 (done with continuation lines)=== lines=%r' % lines
            #print '===end of block 2 (done with continuation lines)=== lines=%r comments=%s' % (lines, comments)
            #pass

            #if not in_loop:
                #self._stored_Is.append(i)
                #self._stored_lines[self._ifile].append(line)
                #self._stored_comments[self._ifile].append(comment)
                #print "non-continuation line = ", line
                #print "lines = ", lines

            #-----------------------------------------------------------------
            # We maybe got one too many lines
            if line[0] not in [' ', '*', '+', ',', '\t']:
            #if in_loop:
                #print "stored_lines =", self._stored_lines
                #print "storing line..."
                self._stored_Is[self._ifile].append(i)
                self._stored_lines[self._ifile].append(line)
                comment = _clean_comment(comment)
                if comment:
                    self._stored_comments[self._ifile].append(comment)
            #if comment:
            #    self._stored_comments.append(comment)
            #self.log.debug("linesE = %s" % lines)

            #print "lines =",lines

            lines2 = clean_empty_lines(lines)
            comment = ''.join(comments)
            #if 'PBARL' in lines[0]:
            #self.log.debug("*yielding lines2=%r" % lines2)
            #self.log.debug("*yielding lines2=%r comments=%s" %(lines2, comments))
            yield lines2, comment
            #print '--------------------------'
        return
        #print "end of stream card..."

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
        assert ' ' not in card_name and len(card_name) > 0, 'card_name=|%r|\nline=|%s| in filename=%s is invalid' % (card_name, lines[0], self.active_filename)
        return card_name.upper()

    def _read_bulk_data_deck(self):
        """
        Parses the Bulk Data Deck

        :param self: the BDF object
        """
        self.log.info("reading Bulk Data Deck...")
        self._break_comment = True
        n = 1
        #isDone = False
        isEndData = False
        icard = 1
        while self.active_filename: # or self._stored_lines:
            #self.log.debug("self._stored_lines = %s" % self._stored_lines)
            try:
                (lines, comment) = self._card_streams[self._ifile].next()
            except StopIteration:
                self._close_file()
                continue
            assert len(lines) > 0
            #if not lines:
                #print "closing from get_card_name"
                #self._close_file()
                #if self.active_filename:
                    #continue
                #break
            #isEndData = self._process_bulk_card(lines, comment)
            n += 1

            card_name = self._get_card_name(lines)
            assert is_string(comment), type(comment)
            #comment = ''.join(comments)
            #print "*asdf*lines = %r" % lines

            if card_name == 'INCLUDE':
                bdf_filename = get_include_filename(lines, include_dir=self.include_dir)
                #print "newfname =", newfname
                self._open_file(bdf_filename)
                reject = '$ INCLUDE processed:  %s\n' % bdf_filename
                #print reject
                if comment:
                    self.rejects.append([comment])
                self.rejects.append([reject])
                #print "reject return False"
                continue
            elif 'ENDDATA' in card_name:
                isEndData = True  # exits while loop
                #print("found ENDDATA in %r" % self.active_filename)
            #if isEndData:
                break

            self._increase_card_count(card_name)
            #print('card_count =', self.card_count.keys())
            if not self.is_reject(card_name):
                self.add_card(icard, lines, card_name, comment, is_list=False)
                icard += 1
            else:
                #print('card_count =', self.card_count.keys())
                if comment:
                    self.rejects.append([comment])
                self.rejects.append(lines)
                #print "reject return False"

            #print '*'*60
            #print "self.active_filename =", self.active_filename
        #print "normal exit from _read_bulk_data_deck"
        #self._close_file()
        #del self._line_streams

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
        #print "ifile =", self._ifile
        try:
            return self._line_streams[self._ifile].next()
        except StopIteration:
            #print "ifile =", self._ifile
            self._close_file()
            #print "ifile =", self._ifile
            return self._get_line()
            #print "ifile =", self._ifile
            #lines = self._stored_lines[self._ifile]
            #comments = self._stored_comments[self._ifile]
            #return lines, comments
        except KeyError:
            #return None, None, None
            return

    def _stream_line(self):
        with open(self.active_filename, 'r') as f:
            for n,line in enumerate(f):
                line = line.rstrip('\t\r\n ')
                comment = ''
                if self._break_comment and '$' in line:
                    i = line.index('$')
                    comment = line[i:] + '\n'
                    line = line[:i].rstrip('\t ')

                #self.log.debug("  stored lines = %s" % self._stored_lines)
                #self.log.debug("  yieldingB new line=|%s|" % line)
                #self.log.debug("############################")
                yield n, line, comment

                while self._stored_lines[self._ifile]:
                    comment = ''
                    i2 = self._stored_Is[self._ifile].pop(0)
                    line2 = self._stored_lines[self._ifile].pop(0)
                    if self._stored_comments:
                        comment = ''.join(self._stored_comments[self._ifile])
                        self._stored_comments[self._ifile] = []

                    #self.log.debug("  yieldingA line=|%s|" % line2)
                    #self.log.debug("############################")
                    yield i2, line2, comment
        #file.close()

    def process_card(self, card_lines, debug=False):
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

    def add_card(self, icard, card_lines, card_name, comment='', is_list=True):
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
        #self.log.info('card_lines = %s' % card_lines)
        #self.log.info('-'*80)
        #comment = ''.join(comment)
        #self.log.debug("card_name = |%r|" % card_name)
        #print("add_card; card_name=%r is_list=%s" % (card_name, is_list))
        if card_name in ['DEQATN']:
            card_obj = card_lines
            card = card_lines
        else:
            if is_list:
                fields = card_lines
            else:
                fields = to_fields(card_lines, card_name)

            # apply OPENMDAO syntax
            #print("_is_dynamic_syntax =", self._is_dynamic_syntax)
            if self._is_dynamic_syntax:
                fields = [self._parse_dynamic_syntax(field) if '%' in
                          field[0:1] else field for field in fields]

            card = wipe_empty_fields([interpret_value(field, fields) if field is not None
                                      else None for field in fields])
            #print "add_card =", card
            card_obj = BDFCard(card)

        # function that gets by name the initialized object (from global scope)
        
        name = card[0]
        self.write_sorted_card(card_obj, icard)

        
        if icard % 10000 == 0:
            self.log.debug('icard = %i' % icard)
        #print(card_obj)
        #========================
        if name == 'PARAM':
            param = PARAM(card_obj, comment=comment)
            self.add_PARAM(param)
        #========================
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
            #self.elements_shell.add_cquadx(card_obj, comment=comment)
            pass
        elif name == 'CQUADR':
            #self.elements_shell.add_cquadr(card_obj, comment=comment)
            pass
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
            pass


        #========================
        # rotation elements/loads
        elif name == 'CTRIAX':
            self.elements_shell.add_ctriax(card_obj, comment=comment)
        elif name == 'CTRIAX6':
            self.elements_shell.add_ctriax6(card_obj, comment=comment)
        elif name == 'PLOADX1':
            self.ploadx1.add(card_obj, comment=comment)

        #========================
        # elements_solid
        elif name == 'CTETRA':
            if len(card) == 7:
                self.elements_solid.add_ctetra4(card_obj, comment=comment)
            else: # 13
                self.elements_solid.add_ctetra10(card_obj, comment=comment)
        elif name == 'CPENTA':
            if len(card) == 9:
                self.elements_solid.add_cpenta6(card_obj, comment=comment)
            else:
                print('len(CPENTA) =', len(card_obj))
                self.elements_solid.add_cpenta15(card_obj, comment=comment)
        elif name == 'CHEXA':
            if len(card) == 11:
                self.elements_solid.add_chexa8(card_obj, comment=comment)
            else:
                print('len(CHEXA) =', len(card_obj))
                self.elements_solid.add_chexa20(card_obj, comment=comment)

        elif name == 'PSOLID':
            self.properties_solid.add_psolid(card_obj, comment=comment)
        elif name == 'PLSOLID':
            self.properties_solid.add_plsolid(card_obj, comment=comment)
            pass
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
            #self.ctube.add(card_obj, comment=comment)
            pass
        elif name == 'PROD':
            #self.ptube.add(card_obj, comment=comment)
            pass

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
        # applied temperature
        elif name == 'TEMP':
            #self.temp.add(card_obj, comment=comment)
            pass
        elif name == 'TEMPD':
            #self.tempd.add(card_obj, comment=comment)
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
            #self.aero.add(card_obj, comment=comment)
            pass
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
            pass
        elif name == 'ASET1':
            pass
        elif name == 'QSET1':
            pass
        #========================
        elif name == 'SPLINE1':
            pass
        elif name == 'SPLINE2':
            pass
        elif name == 'SPLINE3':
            pass
        elif name == 'SPLINE4':
            pass
        elif name == 'SPLINE5':
            pass
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


        elif name == 'TRIM':
            trim = TRIM(self)
            trim.add(card_obj, comment=comment)
            self.trim[trim.trim_id] = trim
        #========================
        # bushing
        elif name == 'CBUSH1D':
            self.cbush1d.add(card_obj, comment=comment)
            pass
        elif name == 'CBUSH2D':
            self.cbush2d.add(card_obj, comment=comment)
            pass
        elif name == 'CBUSH':
            self.cbush.add(card_obj, comment=comment)
            pass
        elif name == 'PBUSH':
            self.pbush.add(card_obj, comment=comment)
            pass
        elif name == 'PBUSHT':
            self.pbusht.add(card_obj, comment=comment)
            pass

        #========================
        # fast
        elif name == 'PFAST':
            self.pfast.add(card_obj, comment=comment)
            pass
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
        elif name == 'CDAMP':
            self.cdampt.add(card_obj, comment=comment)
        
        #========================
        # bars
        elif name == 'CBAR':
            self.cbar.add(card_obj, comment=comment)
        elif name == 'PBAR':
            self.properties_bar.add_pbar(card_obj, comment=comment)
        elif name == 'PBARL':
            self.properties_bar.add_pbarl(card_obj, comment=comment)

        # beams
        elif name == 'CBEAM':
            #self.cbeam.add(card_obj, comment=comment)
            pass
        elif name == 'PBEAM':
            #self.pbeam.add(card_obj, comment=comment)
            pass
        elif name == 'PBEAML':
            #self.pbeaml.add(card_obj, comment=comment)
            pass

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
            self.dload[load.load_id].append(load)
        elif name == 'LSEQ':
            #load = LSEQ(self)
            #load.add_from_bdf(card_obj, comment=comment)
            #self.lseq[load.load_id].append(load)
            pass
        elif name == 'SLOAD':
            load = SLOAD(self)
            load.add_from_bdf(card_obj, comment=comment)
            self.sload[load.load_id].append(load)
            #pass
        elif name == 'LOAD':
            load = LOAD(self)
            load.add_from_bdf(card_obj, comment=comment)
            self.load[load.load_id].append(load)
            #self.load.add(card_obj, comment=comment)
            pass
        
        #========================
        # applied loads
        elif name == 'FORCE':
            self.force.add(card_obj, comment=comment)
        elif name == 'FORCE1':
            self.force1.add(card_obj, comment=comment)
        elif name == 'FORCE2':
            self.force2.add(card_obj, comment=comment)
        elif name == 'MOMENT':
            self.moment.add(card_obj, comment=comment)
        elif name == 'MOMENT1':
            self.moment1.add(card_obj, comment=comment)
        elif name == 'MOMENT2':
            self.moment2.add(card_obj, comment=comment)

        elif name == 'GRAV':
            self.grav.add(card_obj, comment=comment)
        #========================
        # pressure loads
        elif name == 'PLOAD':
            self.pload.add(card_obj, comment=comment)
        elif name == 'PLOAD1':
            self.pload1.add(card_obj, comment=comment)
        elif name == 'PLOAD2':
            self.pload2.add(card_obj, comment=comment)
        elif name == 'PLOAD4':
            self.pload4.add(card_obj, comment=comment)
            pass
        #elif name == 'PLOADX1':
            #self.ploadx1.add(card_obj, comment=comment)
            pass

        #========================
        # time loads
        elif name == 'TLOAD1':
            #self.tload1.add(card_obj, comment=comment)
            pass
        elif name == 'TLOAD2':
            #self.tload2.add(card_obj, comment=comment)
            pass

        # frequency loads
        elif name == 'RLOAD1':
            #self.rload1.add(card_obj, comment=comment)
            pass
        elif name == 'RLOAD2':
            #self.rload2.add(card_obj, comment=comment)
            pass

        # other
        elif name == 'RFORCE':
            self.rforce.add(card_obj, comment=comment)
        elif name == 'DAREA':
            #self.darea.add(card_obj, comment=comment)
            pass


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
            #self.mpc.add(card_obj, comment=comment)
            pass
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
            card = EIGB(card_obj, comment=comment)
            self.eigb[card.sid] = card
        elif name == 'EIGC':
            card = EIGC(card_obj, comment=comment)
            self.eigc[card.sid] = card
        elif name == 'EIGR':
            card = EIGR(card_obj, comment=comment)
            self.eigr[card.sid] = card
        elif name == 'EIGRL':
            card = EIGRL(card_obj, comment=comment)
            self.eigrl[card.sid] = card

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
            pass
        elif name == 'MAT4':
            #self.materials.add_mat4(card_obj, comment=comment)
            pass
        elif name == 'MAT8':
            #self.materials.add_mat8(card_obj, comment=comment)
            pass
        elif name == 'MAT10':
            self.materials.add_mat10(card_obj, comment=comment)

        #========================
        elif name == 'RBAR':
            pass
        elif name == 'RBE2':
            pass
        elif name == 'RBE3':
            pass
        #========================
        elif name == 'CORD1R':
            self.cord1r.add(card, comment)
        elif name == 'CORD1C':
            self.cord1c.add(card, comment)
        elif name == 'CORD1S':
            self.cord1s.add(card, comment)

        elif name == 'CORD2R':
            self.cord2r.add(card, comment)
        elif name == 'CORD2C':
            self.cord2c.add(card, comment)
        elif name == 'CORD2S':
            self.cord2s.add(card, comment)

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
        elif name == 'RINGAX':
            self.ringax.add(card_obj, comment=comment)
        #========================
        # nonlinear
        elif name == 'TSTEP':
            pass
        elif name == 'NLPARM':
            card = NLPARM()
            card.add(card_obj, comment=comment)
            self.nlparm[card.nlparm_id] = card
        elif name == 'NLPCI':
            card = NLPCI()
            card.add(card_obj, comment=comment)
            self.nlpci[card.nlpci_id] = card
        #========================
        else:
            raise NotImplementedError(name)
        return

    def print_filename(self, filename):
        """
        Takes a path such as C:/work/fem.bdf and locates the file using
        relative paths.  If it's on another drive, the path is not modified.

        :param self:     the BDF object
        :param filename: a filename string
        :returns filename_string: a shortened representation of the filename
        """
        driveLetter = os.path.splitdrive(os.path.abspath(filename))[0]
        if driveLetter == os.path.splitdrive(os.curdir)[0] and self._relpath:
            return os.path.relpath(filename)
        return filename

    def card_stats(self, return_type='string'):
        """
        Print statistics for the BDF

        :param self: the BDF object
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

        msg += self.elements_shell.get_stats()
        msg += self.properties_shell.get_stats()

        msg += self.elements_solid.get_stats()
        msg += self.properties_solid.get_stats()

        msg += self.materials.get_stats()
        
        # rejects
        if self.rejects:
            msg.append('Rejected Cards')
            for name, counter in sorted(self.card_count.iteritems()):
                if name not in self.cards_to_read:
                    msg.append('  %-8s %s' % (name + ':', counter))
        msg.append('')
        if return_type == 'string':
            return '\n'.join(msg)
        else:
            return msg

    def get_SPCx_ids(self, exclude_spcadd=False):
        """
        :param exclude_spcadd: you can exclude SPCADD if you just want a
            list of all the SPCs in the model.  For example, apply all the
            SPCs when there is no SPC=N in the case control deck, but you
            don't need to apply SPCADD=N twice.
        """
        spcs = {
            'SPC': self.spc,
            'SPC1': self.spc1,
            #'SPCAX': self.spcax,
            'SPCD': self.spcd,
        }
        if not exclude_spcadd:
            spcs['SPCADD'] = self.spcadd,

        spc_ids = []
        for spc_type, spc in spcs.iteritems():
            spc_ids.extend(spc.keys())
        return unique(spc_ids)

    def SPC(self, spc_id, used_ids=[]):
        used_ids.append(spc_id)
        spcs = {
            'SPCADD' : self.spcadd,
            'SPC': self.spc,
            'SPC1': self.spc1,
            #'SPCAX': self.spcax,
            'SPCD': self.spcd,
        }
        spc_out = []
        for spc_type, spc in spcs.iteritems():
            if spc_id in spc:
                out = spc[spc_id]
                if spc_type == 'SPCADD':
                    for spci in out.spc_ids:
                        if spci in used_ids:
                            raise RuntimeError('duplicate SPC id=%i' % spci)
                        spc_outi = self.SPC(spci, used_ids)
                        for spcii in spc_outi:
                            spc_out.append(spcii)
                else:
                    spc_out.append(out)
        print("spc_out =", spc_out)
        return spc_out

    def MPC(self, mpc_id, used_ids=[]):
        used_ids.append(mpc_id)
        mpcs = {
            'MPCADD' : self.mpcadd,
            'MPC': self.mpc,
        }
        mpc_out = []
        for mpc_type, mpc in mpcs.iteritems():
            if mpc_id in mpc:
                out = mpc[mpc_id]
                if mpc_type == 'MPCADD':
                    for spci in out.mpc_ids:
                        if mpci in used_ids:
                            raise RuntimeError('duplicate MPC id=%i' % mpci)
                        mpc_outi = self.MPC(mpci, used_ids)
                        for mpcii in mpc_outi:
                            mpc_out.append(mpcii)
                else:
                    mpc_out.append(out)
        print("mpc_out =", mpc_out)
        return mpc_out

    def _get_mass_types(self):
        types = [
            # O-D
            self.mass,

            # 1-D
            self.cbar, self.conrod, self.crod, #self.ctube, self.cbeam,

            # 2-D
            self.elements_shell,

            # 3-D
            self.elements_solid,]
        return reduce_types(types)

    def _get_stiffness_types(self):
        types = [
            # O-D
            self.elements_spring,
            #celas1, self.celas2, self.celas3, self.celas4,

            # 1-D
            self.cbar, self.conrod, self.crod,

            # 2-D
            self.elements_shell,

            # 3-D
            self.elements_solid,]
        return reduce_types(types)

    def _get_damping_types(self):
        types = [
            # O-D
            self.cdamp1, self.cdamp2, self.cdamp3, self.cdamp4, self.cdamp5,

            # 1-D
            self.cbar, self.conrod, self.crod,

            # 2-D
            self.elements_shell,

            # 3-D
            self.elements_solid,]
        return reduce_types(types)

    def mass_properties(self, total=False):
        mass_types = self._get_mass_types()
        massi = []
        for mass_type in mass_types:
            massii = mass_type.get_mass(total=False)
            assert massii is not None, mass_type
            assert not isinstance(massii, float), mass_type
            #print("f massii =", massii)
            massi.extend(massii)

        massi = array(massi)
        total = True
        if total:
            mass = massi.sum()
        else:
            mass = massi
        #print('mass =', mass, type(mass))
        return mass

def reduce_types(types):
    types2 = []
    for etype in types:
        if etype.n:
            types2.append(etype)
    return types2

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
                   '$SETS', '$CONTACT', '$REJECTS', '$REJECT_LINES']:
        comment = ''
    return comment


if __name__ == '__main__':
    bdf = BDF()
    import pyNastran
    pkg_path = pyNastran.__path__[0]
    bdfname = sys.argv[1]
    #print "bdfname =", bdfname
    bdf.read_bdf(bdfname)
    bdf.write_bdf('fem.out.bdf')
