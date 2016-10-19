"""defines the BDF attributes"""
from __future__ import print_function
from collections import defaultdict

from pyNastran.utils import object_attributes, object_methods

# nodes
from pyNastran.bdf.dev_vectorized.cards.coord import Coord
from pyNastran.bdf.dev_vectorized.cards.nodes.grid import GRID, GRDSET
from pyNastran.bdf.dev_vectorized.cards.nodes.spoint import SPOINT
from pyNastran.bdf.dev_vectorized.cards.nodes.point import POINT
from pyNastran.bdf.dev_vectorized.cards.nodes.spoint import SPOINT
from pyNastran.bdf.dev_vectorized.cards.nodes.epoint import EPOINT
from pyNastran.bdf.dev_vectorized.cards.nodes.pointax import POINTAX

# springs
from pyNastran.bdf.dev_vectorized.cards.elements.spring.pelas import PELAS
from pyNastran.bdf.dev_vectorized.cards.elements.spring.celas1 import CELAS1
from pyNastran.bdf.dev_vectorized.cards.elements.spring.celas2 import CELAS2
from pyNastran.bdf.dev_vectorized.cards.elements.spring.celas3 import CELAS3
from pyNastran.bdf.dev_vectorized.cards.elements.spring.celas4 import CELAS4

# rods/tubes
from pyNastran.bdf.dev_vectorized.cards.elements.rod.prod import PROD
from pyNastran.bdf.dev_vectorized.cards.elements.rod.crod import CROD
from pyNastran.bdf.dev_vectorized.cards.elements.rod.conrod import CONROD
from pyNastran.bdf.dev_vectorized.cards.elements.rod.ptube import PTUBE
from pyNastran.bdf.dev_vectorized.cards.elements.rod.ctube import CTUBE

# bars
from pyNastran.bdf.dev_vectorized.cards.elements.bar.cbar import CBAR
#from pyNastran.bdf.dev_vectorized.cards.elements.bar.cbaror import CBAROR
from pyNastran.bdf.dev_vectorized.cards.elements.bar.pbar import PBAR
from pyNastran.bdf.dev_vectorized.cards.elements.bar.pbarl import PBARL
from pyNastran.bdf.dev_vectorized.cards.elements.bar.properties_bar import PropertiesBar

# beams
from pyNastran.bdf.dev_vectorized.cards.elements.beam.cbeam import CBEAM
from pyNastran.bdf.dev_vectorized.cards.elements.beam.pbeam import PBEAM
from pyNastran.bdf.dev_vectorized.cards.elements.beam.pbeaml import PBEAML
from pyNastran.bdf.dev_vectorized.cards.elements.beam.properties_beam import PropertiesBeam

# shear
from pyNastran.bdf.dev_vectorized.cards.elements.shear.cshear import CSHEAR
from pyNastran.bdf.dev_vectorized.cards.elements.shear.pshear import PSHEAR

# shells
from pyNastran.bdf.dev_vectorized.cards.elements.shell.cquad4 import CQUAD4
from pyNastran.bdf.dev_vectorized.cards.elements.shell.ctria3 import CTRIA3
from pyNastran.bdf.dev_vectorized.cards.elements.shell.ctria6 import CTRIA6
from pyNastran.bdf.dev_vectorized.cards.elements.shell.cquad8 import CQUAD8
#from pyNastran.bdf.dev_vectorized.cards.elements.shell.cquad import CQUAD
from pyNastran.bdf.dev_vectorized.cards.elements.shell.pshell import PSHELL
from pyNastran.bdf.dev_vectorized.cards.elements.shell.pcomp import PCOMP
from pyNastran.bdf.dev_vectorized.cards.elements.shell.pcompg import PCOMPG

# solids
from pyNastran.bdf.dev_vectorized.cards.elements.solid.psolid import PSOLID
from pyNastran.bdf.dev_vectorized.cards.elements.solid.ctetra4 import CTETRA4
from pyNastran.bdf.dev_vectorized.cards.elements.solid.ctetra10 import CTETRA10
from pyNastran.bdf.dev_vectorized.cards.elements.solid.cpenta6 import CPENTA6
from pyNastran.bdf.dev_vectorized.cards.elements.solid.cpenta15 import CPENTA15
from pyNastran.bdf.dev_vectorized.cards.elements.solid.chexa8 import CHEXA8
from pyNastran.bdf.dev_vectorized.cards.elements.solid.chexa20 import CHEXA20

# mass
from pyNastran.bdf.dev_vectorized.cards.elements.mass.mass import Mass
from pyNastran.bdf.dev_vectorized.cards.elements.mass.conm1 import CONM1
from pyNastran.bdf.dev_vectorized.cards.elements.mass.conm2 import CONM2

#from pyNastran.bdf.dev_vectorized.cards.elements.mass.pmass import PMASS
#from pyNastran.bdf.dev_vectorized.cards.elements.mass.cmass1 import CMASS1
#from pyNastran.bdf.dev_vectorized.cards.elements.mass.cmass2 import CMASS2
#from pyNastran.bdf.dev_vectorized.cards.elements.mass.cmass3 import CMASS3
#from pyNastran.bdf.dev_vectorized.cards.elements.mass.cmass4 import CMASS4
#from pyNastran.bdf.dev_vectorized.cards.elements.mass.cmass5 import CMASS5


# cards
from pyNastran.bdf.dev_vectorized.cards.elements.elements import Elements
from pyNastran.bdf.dev_vectorized.cards.elements.properties import Properties


# materials
from pyNastran.bdf.dev_vectorized.cards.materials.mat1 import MAT1
#from pyNastran.bdf.dev_vectorized.cards.materials.mat2 import MAT2
from pyNastran.bdf.dev_vectorized.cards.materials.mat8 import MAT8


#------------------------------------------
from pyNastran.bdf.dev_vectorized.cards.elements.shell.elements_shell import ElementsShell
from pyNastran.bdf.dev_vectorized.cards.elements.shell.properties_shell import PropertiesShell

# solid
from pyNastran.bdf.dev_vectorized.cards.elements.solid.elements_solid import ElementsSolid
from pyNastran.bdf.dev_vectorized.cards.elements.solid.properties_solid import PropertiesSolid

# spring
from pyNastran.bdf.dev_vectorized.cards.elements.spring.elements_spring import ElementsSpring


class BDFAttributes(object):
    """defines attributes of the BDF"""
    def __init__(self):
        """creates the attributes for the BDF"""
        self.is_msc = True
        self._is_cards_dict = True
        self.reject_lines = []

        self.set_precision()
        #----------------------------------------
        self.params = {}
        self.grid = GRID(self)
        self.point = POINT(self)
        self.grdset = GRDSET(self)
        self.spoint = SPOINT(self)
        self.epoint = EPOINT(self)
        self.pointax = POINTAX(self)
        self.coords = Coord(self)
        #----------------------------------------

        # springs
        self.pelas = PELAS(self)
        self.celas1 = CELAS1(self)
        self.celas2 = CELAS2(self)
        self.celas3 = CELAS3(self)
        self.celas4 = CELAS4(self)
        self.elements_spring = ElementsSpring(self)

        # rods/tubes
        self.prod = PROD(self)
        self.crod = CROD(self)
        self.conrod = CONROD(self)
        self.ptube = PTUBE(self)
        self.ctube = CTUBE(self)

        # bars
        self.cbar = CBAR(self)
        #self.cbaror = CBAROR(self)
        self.pbar = PBAR(self)
        self.pbarl = PBARL(self)
        self.properties_bar = PropertiesBar(self)

        # beams
        self.cbeam = CBEAM(self)
        self.pbeam = PBEAM(self)
        self.pbeaml = PBEAML(self)
        #: stores PBEAM, PBEAML
        self.properties_beam = PropertiesBeam(self)

        # shear
        #: stores CSHEAR
        self.cshear = CSHEAR(self)
        #: stores PSHEAR
        self.pshear = PSHEAR(self)

        # shells
        self.pshell = PSHELL(self)
        self.pcomp = PCOMP(self)
        self.pcompg = PCOMPG(self)
        self.cquad4 = CQUAD4(self)
        self.ctria3 = CTRIA3(self)
        self.ctria6 = CTRIA6(self)
        self.cquad8 = CQUAD8(self)
        #self.cquad = CQUAD(self)
        #: stores PSHELL, PCOMP, PCOMPG
        self.properties_shell = PropertiesShell(self)
        #: stores CTRIA3, CTRIA6, CQUAD4, CQUAD8
        self.elements_shell = ElementsShell(self)

        # solids
        self.psolid = PSOLID(self)
        self.ctetra4 = CTETRA4(self)
        self.ctetra10 = CTETRA10(self)
        self.cpenta6 = CPENTA6(self)
        self.cpenta15 = CPENTA15(self)
        self.chexa8 = CHEXA8(self)
        self.chexa20 = CHEXA20(self)
        #: stores CTETRA4, CPENTA6, CHEXA8, CTETRA10, CPENTA15, CHEXA20
        self.elements_solid = ElementsSolid(self)
        #: stores PSOLID, PLSOLID
        self.properties_solid = PropertiesSolid(self)

        #----------------------------------------
        # mass
        self.conm1 = CONM1(self)
        self.conm2 = CONM2(self)
        #self.pmass = PMASS(self)
        #self.cmass1 = CMASS1(self)
        #self.cmass2 = CMASS2(self)
        #self.cmass3 = CMASS3(self)
        #self.cmass4 = CMASS4(self)
        #self.cmass5 = CMASS5(self)
        self.pmass = None
        self.cmass1 = None
        self.cmass2 = None
        self.cmass3 = None
        self.cmass4 = None
        self.cmass5 = None
        self.mass = Mass(self)
        #----------------------------------------
        # b-list elements
        self.rbe2 = None
        self.rbe3 = None
        self.cbush = None
        self.pbush = None
        self.cbush1d = None
        self.pbush1d = None
        self.cbush2d = None
        self.pbush2d = None

        #----------------------------------------
        # control structure
        self.elements = Elements(self)
        self.properties = Properties(self)
        #----------------------------------------



        self.mat1 = MAT1(self)
        #self.mat2 = MAT2(self)
        self.mat8 = MAT8(self)

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
        #: random_tables
        self.random_tables = {}
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
        self.dvgrids = {}
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
        self.monitor_points = []

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
        self.aesurf = {}
        #: stores AESURFS
        self.aesurfs = {}
        #: stores AESTAT
        self.aestats = {}
        #: stores CSSCHD
        self.csschds = {}

        #: store SPLINE1,SPLINE2,SPLINE4,SPLINE5
        self.splines = {}

        # ------ SOL 144 ------
        #: stores AEROS
        self.aeros = None

        #: stores TRIM
        self.trims = {}

        #: stores DIVERG
        self.divergs = {}

        # ------ SOL 145 ------
        #: stores AERO
        self.aero = None

        #: stores FLFACT
        self.flfacts = {}  #: .. todo:: can this be simplified ???
        #: stores FLUTTER
        self.flutters = {}
        #: mkaeros
        self.mkaeros = []

        # ------ SOL 146 ------
        #: stores GUST cards
        self.gusts = {}
        # ------------------------- thermal defaults -------------------------
        # BCs
        #: stores thermal boundary conditions - CONV,RADBC
        self.bcs = {}  # e.g. RADBC

        #: stores PHBDY
        self.phbdys = {}
        #: stores convection properties - PCONV, PCONVM ???
        self.convection_properties = {}
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

        self.rsolmap_to_str = {
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

    def set_as_msc(self):
        self._nastran_format = 'msc'
        self.is_nx = False
        self.is_msc = True

    def set_as_nx(self):
        self._nastran_format = 'nx'
        self.is_nx = True
        self.is_msc = False


    def object_attributes(self, mode='public', keys_to_skip=None):
        """
        List the names of attributes of a class as strings. Returns public
        attributes as default.

        Parameters
        ----------
        mode : str
            defines what kind of attributes will be listed
            * 'public' - names that do not begin with underscore
            * 'private' - names that begin with single underscore
            * 'both' - private and public
            * 'all' - all attributes that are defined for the object
        keys_to_skip : List[str]; default=None -> []
            names to not consider to avoid deprecation warnings

        Returns
        -------
        attribute_names : List[str]
            sorted list of the names of attributes of a given type or None
            if the mode is wrong
        """
        if keys_to_skip is None:
            keys_to_skip = []

        my_keys_to_skip = [
            #'case_control_deck',
            'log', 'mpcObject', 'spcObject',
            'node_ids', 'coord_ids', 'element_ids', 'property_ids',
            'material_ids', 'caero_ids', 'is_long_ids',
            'nnodes', 'ncoords', 'nelements', 'nproperties',
            'nmaterials', 'ncaeros',

            'point_ids', 'subcases',
            '_card_parser', '_card_parser_b',
            'object_methods', 'object_attributes',
        ]
        return object_attributes(self, mode=mode, keys_to_skip=keys_to_skip+my_keys_to_skip)

    def object_methods(self, mode='public', keys_to_skip=None):
        """
        List the names of methods of a class as strings. Returns public methods
        as default.

        Parameters
        ----------
        obj : instance
            the object for checking
        mode : str
            defines what kind of methods will be listed
            * "public" - names that do not begin with underscore
            * "private" - names that begin with single underscore
            * "both" - private and public
            * "all" - all methods that are defined for the object
        keys_to_skip : List[str]; default=None -> []
            names to not consider to avoid deprecation warnings

        Returns
        -------
        method : List[str]
            sorted list of the names of methods of a given type
            or None if the mode is wrong
        """
        if keys_to_skip is None:
            keys_to_skip = []
        my_keys_to_skip = []

        my_keys_to_skip = [
            #'case_control_deck',
            'log', #'mpcObject', 'spcObject',
            'node_ids', 'coord_ids', 'element_ids', 'property_ids',
            'material_ids', 'caero_ids', 'is_long_ids',
            'nnodes', 'ncoords', 'nelements', 'nproperties',
            'nmaterials', 'ncaeros',

            'point_ids', 'subcases',
            '_card_parser', '_card_parser_b',
            'object_methods', 'object_attributes',
        ]
        return object_methods(self, mode=mode, keys_to_skip=keys_to_skip+my_keys_to_skip)


    @property
    def rejects(self):
        """access the rejected lines"""
        #: lines that were rejected b/c they were for a card that isnt supported
        return self.reject_lines

    @rejects.setter
    def rejects(self, rejects):
        """set the rejected lines"""
        self.reject_lines = rejects

    @property
    def card_name_to_obj(self):
        """turns the card_name into the object"""
        mapper = {
            'GRID' : self.grid,
            'GRIDSET' : self.grdset,
            'SPOINT' : self.spoint,

            # springs
            'PELAS' : self.pelas,
            'CELAS1' : self.celas1,
            'CELAS2' : self.celas2,
            'CELAS3' : self.celas3,
            'CELAS4' : self.celas4,

            # rods/tube
            'PROD' : self.prod,
            'CROD' : self.crod,
            'CONROD' : self.conrod,
            'PTUBE' : self.ptube,
            'CTUBE' : self.ctube,

            # bars
            'CBAR' : self.cbar,
            #'CBAROR' : self.cbaror,
            'PBAR' : self.pbar,
            'PBARL' : self.pbarl,

            # beams
            'CBEAM' : self.cbeam,
            'PBEAM' : self.pbeam,
            'PBEAML' : self.pbeaml,

            # shells
            'PSHELL' : self.pshell,
            'PCOMP' : self.pcomp,
            'PCOMPG' : self.pcompg,
            'CTRIA3' : self.ctria3,
            'CQUAD4' : self.cquad4,
            'CQUAD8' : self.cquad8,
            #'CQUAD' : self.cquad,

            # mass
            'CONM1' : self.conm1,
            'CONM2' : self.conm2,
            'PMASS' : self.pmass,
            'CMASS1' : self.cmass1,
            'CMASS2' : self.cmass2,
            'CMASS3' : self.cmass3,
            'CMASS4' : self.cmass4,

            'MAT1' : self.mat1,
            #'MAT2' : self.mat2,
            'MAT8' : self.mat8,
            'PSOLID' : self.psolid,
        }
        return mapper

    def set_precision(self, precision='double'):
        """
        Sets the float precision.

        Parameters
        ----------
        precision : str
            string of 'single'/'float32' or 'double'/'float64'
            that is used by all the objects
        """
        if precision in ('double', 'float64'):
            self.float_fmt = 'float64'
        elif precision == ('single', 'float32'):
            self.float_fmt = 'float32'
        else:
            raise NotImplementedError('precision=%r' % precision)

    @property
    def nastran_format(self):
        return self._nastran_format

    @nastran_format.setter
    def nastran_format(self, nastran_format):
        fmt_lower = nastran_format.lower().strip()
        if fmt_lower not in ['nx', 'msc']:
            raise RuntimeError(nastran_format)
        self._nastran_format = fmt_lower

    @property
    def is_long_ids(self):
        if self._nastran_format == 'nx':
            return True
        return False
