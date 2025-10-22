"""defines the BDF attributes"""
from collections import defaultdict

from pyNastran.utils import object_attributes, object_methods

# nodes
from pyNastran.dev.bdf_vectorized.cards.coord import Coord
from pyNastran.dev.bdf_vectorized.cards.nodes.grid import GRID, GRDSET
from pyNastran.dev.bdf_vectorized.cards.nodes.point import POINT
from pyNastran.dev.bdf_vectorized.cards.nodes.spoint import SPOINT
from pyNastran.dev.bdf_vectorized.cards.nodes.epoint import EPOINT
from pyNastran.dev.bdf_vectorized.cards.nodes.pointax import POINTAX

# springs
from pyNastran.dev.bdf_vectorized.cards.elements.spring.pelas import PELAS
from pyNastran.dev.bdf_vectorized.cards.elements.spring.celas1 import CELAS1
from pyNastran.dev.bdf_vectorized.cards.elements.spring.celas2 import CELAS2
from pyNastran.dev.bdf_vectorized.cards.elements.spring.celas3 import CELAS3
from pyNastran.dev.bdf_vectorized.cards.elements.spring.celas4 import CELAS4

from pyNastran.dev.bdf_vectorized.cards.elements.damper.cdamp1 import CDAMP1
from pyNastran.dev.bdf_vectorized.cards.elements.damper.cdamp2 import CDAMP2
from pyNastran.dev.bdf_vectorized.cards.elements.damper.cdamp3 import CDAMP3
from pyNastran.dev.bdf_vectorized.cards.elements.damper.cdamp4 import CDAMP4

# bush
from pyNastran.dev.bdf_vectorized.cards.elements.bush.cbush import CBUSH
from pyNastran.dev.bdf_vectorized.cards.elements.bush.pbush import PBUSH

# rods/tubes
from pyNastran.dev.bdf_vectorized.cards.elements.rod.prod import PROD
from pyNastran.dev.bdf_vectorized.cards.elements.rod.crod import CROD
from pyNastran.dev.bdf_vectorized.cards.elements.rod.conrod import CONROD
from pyNastran.dev.bdf_vectorized.cards.elements.rod.ptube import PTUBE
from pyNastran.dev.bdf_vectorized.cards.elements.rod.ctube import CTUBE

# bars
from pyNastran.dev.bdf_vectorized.cards.elements.bar.cbar import CBAR
#from pyNastran.dev.bdf_vectorized.cards.elements.bar.baror import BAROR
from pyNastran.dev.bdf_vectorized.cards.elements.bar.pbar import PBAR
from pyNastran.dev.bdf_vectorized.cards.elements.bar.pbarl import PBARL
from pyNastran.dev.bdf_vectorized.cards.elements.bar.properties_bar import PropertiesBar

# beams
from pyNastran.dev.bdf_vectorized.cards.elements.beam.cbeam import CBEAM
from pyNastran.dev.bdf_vectorized.cards.elements.beam.pbeam import PBEAM
from pyNastran.dev.bdf_vectorized.cards.elements.beam.pbeaml import PBEAML
from pyNastran.dev.bdf_vectorized.cards.elements.beam.properties_beam import PropertiesBeam

# shear
from pyNastran.dev.bdf_vectorized.cards.elements.shear.cshear import CSHEAR
from pyNastran.dev.bdf_vectorized.cards.elements.shear.pshear import PSHEAR

# shells
from pyNastran.dev.bdf_vectorized.cards.elements.shell.cquad4 import CQUAD4
from pyNastran.dev.bdf_vectorized.cards.elements.shell.ctria3 import CTRIA3
from pyNastran.dev.bdf_vectorized.cards.elements.shell.ctria6 import CTRIA6
from pyNastran.dev.bdf_vectorized.cards.elements.shell.cquad8 import CQUAD8
#from pyNastran.dev.bdf_vectorized.cards.elements.shell.cquad import CQUAD
#from pyNastran.dev.bdf_vectorized.cards.elements.shell.cquadx import CQUADX
#from pyNastran.dev.bdf_vectorized.cards.elements.shell.cquad9 import CQUAD9
#from pyNastran.dev.bdf_vectorized.cards.elements.shell.ctriax import CTRIAX
from pyNastran.dev.bdf_vectorized.cards.elements.shell.ctriax6 import CTRIAX6
from pyNastran.dev.bdf_vectorized.cards.elements.shell.pshell import PSHELL
from pyNastran.dev.bdf_vectorized.cards.elements.shell.pcomp import PCOMP
from pyNastran.dev.bdf_vectorized.cards.elements.shell.pcompg import PCOMPG

# solids
from pyNastran.dev.bdf_vectorized.cards.elements.solid.psolid import PSOLID
from pyNastran.dev.bdf_vectorized.cards.elements.solid.plsolid import PLSOLID
from pyNastran.dev.bdf_vectorized.cards.elements.solid.ctetra4 import CTETRA4
from pyNastran.dev.bdf_vectorized.cards.elements.solid.ctetra10 import CTETRA10
from pyNastran.dev.bdf_vectorized.cards.elements.solid.cpenta6 import CPENTA6
from pyNastran.dev.bdf_vectorized.cards.elements.solid.cpenta15 import CPENTA15
from pyNastran.dev.bdf_vectorized.cards.elements.solid.chexa8 import CHEXA8
from pyNastran.dev.bdf_vectorized.cards.elements.solid.chexa20 import CHEXA20

# mass
from pyNastran.dev.bdf_vectorized.cards.elements.mass.mass import Mass
from pyNastran.dev.bdf_vectorized.cards.elements.mass.conm1 import CONM1
from pyNastran.dev.bdf_vectorized.cards.elements.mass.conm2 import CONM2
#from pyNastran.dev.bdf_vectorized.cards.elements.mass.pmass import PMASS
#from pyNastran.dev.bdf_vectorized.cards.elements.mass.cmass1 import CMASS1
#from pyNastran.dev.bdf_vectorized.cards.elements.mass.cmass2 import CMASS2
#from pyNastran.dev.bdf_vectorized.cards.elements.mass.cmass3 import CMASS3
#from pyNastran.dev.bdf_vectorized.cards.elements.mass.cmass4 import CMASS4
#from pyNastran.dev.bdf_vectorized.cards.elements.mass.cmass5 import CMASS5


# cards
from pyNastran.dev.bdf_vectorized.cards.elements.elements import Elements
from pyNastran.dev.bdf_vectorized.cards.elements.properties import Properties


# materials
from pyNastran.dev.bdf_vectorized.cards.materials.mat1 import MAT1
#from pyNastran.dev.bdf_vectorized.cards.materials.mat2 import MAT2
from pyNastran.dev.bdf_vectorized.cards.materials.mat8 import MAT8
from pyNastran.dev.bdf_vectorized.cards.materials.mats1 import MATS1

from pyNastran.dev.bdf_vectorized.cards.materials.mathp import MATHP
from pyNastran.dev.bdf_vectorized.cards.materials.materials import Materials

# constraints
# spc
#from pyNastran.dev.bdf_vectorized.cards.constraints.spcadd import SPCADD

# mpc
#from pyNastran.dev.bdf_vectorized.cards.constraints.mpc import MPC
#from pyNastran.dev.bdf_vectorized.cards.constraints.mpcax import MPCAX
#from pyNastran.dev.bdf_vectorized.cards.constraints.mpcadd import MPCADD

#------------------------------------------
from pyNastran.dev.bdf_vectorized.cards.elements.shell.elements_shell import ElementsShell
from pyNastran.dev.bdf_vectorized.cards.elements.shell.properties_shell import PropertiesShell

# solid
from pyNastran.dev.bdf_vectorized.cards.elements.solid.elements_solid import ElementsSolid
from pyNastran.dev.bdf_vectorized.cards.elements.solid.properties_solid import PropertiesSolid

# spring
from pyNastran.dev.bdf_vectorized.cards.elements.spring.elements_spring import ElementsSpring

#-------------------------------------------------------------------------------
# loads
from pyNastran.dev.bdf_vectorized.cards.loads.loads import Loads, LOADs
from pyNastran.dev.bdf_vectorized.cards.loads.temp import TEMPs, TEMPP1, TEMP

#from pyNastran.dev.bdf_vectorized.cards.loads.load import LOAD
#from pyNastran.dev.bdf_vectorized.cards.loads.dload import DLOAD
#from pyNastran.dev.bdf_vectorized.cards.loads.dload import DLOAD as LSEQ
#from pyNastran.dev.bdf_vectorized.cards.loads.dload import DLOAD as SLOAD
from pyNastran.dev.bdf_vectorized.cards.loads.static_loads import (
    GRAV, FORCE, MOMENT, FORCE1, FORCE2, MOMENT1, MOMENT2,
    PLOAD, PLOAD1, PLOAD2, PLOAD4
)

# ACCEL1
# PLOAD3
from pyNastran.dev.bdf_vectorized.cards.loads.darea import DAREA

from pyNastran.dev.bdf_vectorized.cards.loads.tload1 import TLOAD1
from pyNastran.dev.bdf_vectorized.cards.loads.tload2 import TLOAD2
from pyNastran.dev.bdf_vectorized.cards.loads.delay import DELAY
from pyNastran.dev.bdf_vectorized.cards.loads.rload1 import RLOAD1
#from pyNastran.dev.bdf_vectorized.cards.loads.rload2 import RLOAD2
from pyNastran.dev.bdf_vectorized.cards.loads.dphase import DPHASE
# RANDPS


from pyNastran.dev.bdf_vectorized.cards.loads.ploadx1 import PLOADX1

from pyNastran.dev.bdf_vectorized.cards.loads.rforce import RFORCE
#from pyNastran.dev.bdf_vectorized.cards.loads.sload import SLOAD

#from pyNastran.dev.bdf_vectorized.cards.loads.loadcase import LoadCase
#from pyNastran.dev.bdf_vectorized.cards.loads.loadset import LOADSET
#-------------------------------------------------------------------------------


class BDFAttributes:
    """defines attributes of the BDF"""

    def __init__(self):
        """creates the attributes for the BDF"""
        self._nastran_format = 'msc'
        self.is_nx = False
        self.is_msc = True
        self.max_int = 100000000

        #----------------------------------------
        self._is_cards_dict = True
        self.punch = None
        self.reject_lines = []

        self.set_precision()
        #----------------------------------------
        self.grid = GRID(self)
        self.grdset = GRDSET(self)
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

        self.celas1 = CDAMP1(self)
        self.celas2 = CDAMP2(self)
        self.celas3 = CDAMP3(self)
        self.celas4 = CDAMP4(self)

        # rods/tubes
        self.prod = PROD(self)
        self.crod = CROD(self)
        self.conrod = CONROD(self)
        self.ptube = PTUBE(self)
        self.ctube = CTUBE(self)

        # bars
        self.cbar = CBAR(self)
        #self.baror = BAROR(self)
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
        self.ctriax6 = CTRIAX6(self)
        #self.cquad = CQUAD(self)
        #: stores PSHELL, PCOMP, PCOMPG
        self.properties_shell = PropertiesShell(self)
        #: stores CTRIA3, CTRIA6, CQUAD4, CQUAD8
        self.elements_shell = ElementsShell(self)

        # solids
        self.psolid = PSOLID(self)
        self.plsolid = PLSOLID(self)
        self.ctetra4 = CTETRA4(self)
        self.ctetra10 = CTETRA10(self)
        self.cpyram5 = None
        self.cpyram13 = None
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
        #self.rbe2 = None
        #self.rbe3 = None
        self.cbush = CBUSH(self)
        self.pbush = PBUSH(self)
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
        self.mats1 = MATS1(self)
        #self.mat2 = MAT2(self)
        #self.mat2 = MAT2(self)
        #self.mat4 = MAT4(self)
        #self.mat5 = MAT5(self)
        self.mat8 = MAT8(self)
        #self.mat10 = MAT10(self)
        #self.mat11 = MAT11(self)
        self.mathp = MATHP(self)

        self.materials = Materials(self)

        # ----------------------------------------------------------------

        self.load = LOADs(self)
        self.dload = LOADs(self)
        #self.dload = defaultdict(list)
        #self.loadset = LOADSET(model)

        self.force = FORCE(self)
        self.force1 = FORCE1(self)
        self.force2 = FORCE2(self)
        self.moment = MOMENT(self)
        self.moment1 = MOMENT1(self)
        self.moment2 = MOMENT2(self)
        self.grav = GRAV(self)
        self.rforce = RFORCE(self)

        self.pload = PLOAD(self)
        self.pload1 = PLOAD1(self)
        self.pload2 = PLOAD2(self)
        #self.pload3 = PLOAD3(self)
        self.pload4 = PLOAD4(self)
        self.ploadx1 = PLOADX1(self)

        self.tload1 = TLOAD1(self)
        self.tload2 = TLOAD2(self)
        self.delay = DELAY(self)

        self.rload1 = RLOAD1(self)
        #self.rload2 = RLOAD2(self)
        self.dphase = DPHASE(self)

        self.darea = DAREA(self)

        #: stores LOAD, FORCE, MOMENT, etc.
        self.loads = Loads(self)
        # ----------------------------------------------------------------
        self.tempp1 = TEMPP1(self)
        self.temp = TEMP(self)
        self.temps = TEMPs(self)

        # ----------------------------------------------------------------
        #self.spc1 = SPC1(self)
        #self.spcadd = SPCADD(self)
        self.spc = {} #class_obj_defaultdict(SPC, model)
        self.spcd = {} #class_obj_defaultdict(SPCD, model)
        self.spc1 = {} #class_obj_defaultdict(SPC1, model)
        self.spcadd = {}
        self.mpc = {}  # the new form, not added...
        self.mpcadd = {}

        # ----------------------------------------------------------------
        #: stores PARAMs
        self.params = {}

        #: stores rigid elements (RBE2, RBE3, RJOINT, etc.)
        self.rigid_elements = {}
        #: stores PLOTELs
        self.plotels = {}

        # --------------------------- dynamic ----------------------------
        #: stores DAREA
        #self.dareas = {}
        #self.dphases = {}

        self.pbusht = {}
        self.pdampt = {}
        self.pelast = {}

        #: frequencies
        self.frequencies = {}

        # ----------------------------------------------------------------
        #: direct matrix input - DMIG
        self.dmi = {}
        self.dmig = {}
        self.dmij = {}
        self.dmiji = {}
        self.dmik = {}
        self._dmig_temp = defaultdict(list)

        # ----------------------------------------------------------------
        #: SETy
        self.sets = {} # SET1, SET3
        self.asets = []
        self.bsets = []
        self.csets = []
        self.qsets = []
        self.omits = []
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
        self.trims: dict[key, TRIM] = {}

        #: stores DIVERG
        self.divergs: dict[key, DIVERG] = {}

        # ------ SOL 145 ------
        #: stores AERO
        self.aero = None

        #: stores FLFACT
        self.flfacts: dict[key, FLFACT] = {}  #: .. todo:: can this be simplified ???
        #: stores FLUTTER
        self.flutters: dict[key, FLUTTER] = {}
        #: mkaeros
        self.mkaeros = []

        # ------ SOL 146 ------
        #: stores GUST cards
        self.gusts: dict[key, GUST] = {}
        # ------------------------- thermal defaults -------------------------
        # BCs
        #: stores thermal boundary conditions - CONV,RADBC
        self.bcs: dict[key, CONV | RADBC] = {}

        #: stores PHBDY
        self.phbdys = {}  # type: dict[key, Any]
        #: stores convection properties - PCONV, PCONVM ???
        self.convection_properties: dict[key, PCONV] = {}
        #: stores TEMPD
        self.tempds: dict[key, TEMPD] = {}

        # -------------------------contact cards-------------------------------
        self.bcrparas: dict[key, BCRPARA] = {}
        self.bctadds: dict[key, BCTADD] = {}
        self.bctparas: dict[key, BCTPARA] = {}
        self.bctsets: dict[key, BCTSET] = {}
        self.bsurf: dict[key, BURF] = {}
        self.bsurfs: dict[key, BSURFS] = {}

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
            'MFREQI': 111,
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


    def object_attributes(self, mode='public', keys_to_skip=None,
                          filter_properties=False):
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
        keys_to_skip : list[str]; default=None -> []
            names to not consider to avoid deprecation warnings

        Returns
        -------
        attribute_names : list[str]
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
            '_card_parser',
            'object_methods', 'object_attributes',
        ]
        return object_attributes(self, mode=mode, keys_to_skip=keys_to_skip+my_keys_to_skip,
                                 filter_properties=filter_properties)

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
        keys_to_skip : list[str]; default=None -> []
            names to not consider to avoid deprecation warnings

        Returns
        -------
        method : list[str]
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
            '_card_parser',
            'object_methods', 'object_attributes',
        ]
        return object_methods(self, mode=mode, keys_to_skip=keys_to_skip+my_keys_to_skip)

    @property
    def subcases(self):
        """gets the subcases"""
        if self.case_control_deck is None:
            return {}
        return self.case_control_deck.subcases

    @property
    def rejects(self):
        """access the rejected lines"""
        #: lines that were rejected b/c they were for a card that isn't supported
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
            'EPOINT' : self.epoint,

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
            #'BAROR' : self.baror,
            'CBAR' : self.cbar,
            'PBAR' : self.pbar,
            'PBARL' : self.pbarl,

            # beams
            #'BEAMOR' : self.beamor,
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

            # solid
            'PSOLID' : self.psolid,
            'PLSOLID' : self.plsolid,

            # mass
            'CONM1' : self.conm1,
            'CONM2' : self.conm2,
            'PMASS' : self.pmass,
            'CMASS1' : self.cmass1,
            'CMASS2' : self.cmass2,
            'CMASS3' : self.cmass3,
            'CMASS4' : self.cmass4,

            'CBUSH' : self.cbush,
            'PBUSH' : self.pbush,

            'CBUSH1D' : self.cbush1d,
            'PBUSH1D' : self.pbush1d,

            'CBUSH2D' : self.cbush2d,
            'PBUSH2D' : self.pbush2d,

            'MAT1' : self.mat1,
            #'MAT2' : self.mat2,
            'MAT8' : self.mat8,

            # loads
            #'LOAD' : self.load,
            'GRAV' : self.grav,

            'FORCE' : self.force,
            'FORCE1' : self.force1,
            'FORCE2' : self.force2,

            'MOMENT' : self.moment,
            'MOMENT1' : self.moment1,
            'MOMENT2' : self.moment2,

            'PLOAD' : self.pload,
            'PLOAD2' : self.pload2,
            'PLOAD4' : self.pload4,
            'PLOADX1' : self.ploadx1,

            'TLOAD1' : self.tload1,
            'TLOAD2' : self.tload2,
            'DELAY' : self.delay,
            'RLOAD1' : self.rload1,
            #'RLOAD2' : self.rload2,
            'DPHASE' : self.dphase,

            # constraints
            #'SPC1' : self.spc1,
            #'SPCADD' : self.spcadd,
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
        else:  # pragma: no cover
            raise NotImplementedError('precision=%r' % precision)

    @property
    def nastran_format(self) -> str:
        return self._nastran_format

    @nastran_format.setter
    def nastran_format(self, nastran_format: str) -> None:
        fmt_lower = nastran_format.lower().strip()
        if fmt_lower not in ['nx', 'msc']:
            raise RuntimeError(nastran_format)
        self._nastran_format = fmt_lower
        map_version(fmt_lower)

    @property
    def is_long_ids(self) -> bool:
        if self._nastran_format == 'nx':
            return True
        return False
