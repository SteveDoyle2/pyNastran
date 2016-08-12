from __future__ import print_function
from collections import defaultdict
from numpy import array

from pyNastran.utils import object_attributes, object_methods
from pyNastran.bdf.utils import deprecated
#from pyNastran.bdf.case_control_deck import CaseControlDeck
from pyNastran.bdf.cards.coordinate_systems import CORD2R
#from pyNastran.bdf.cards.constraints import ConstraintObject

class BDFAttributes(object):

    def __init__(self):
        self.__init_attributes()
        self._is_cards_dict = False
        self.set_as_msc()
        self.units = []

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

    def deprecated(self, old_name, new_name, deprecated_version):
        """deprecates methods"""
        return deprecated(old_name, new_name, deprecated_version, levels=[0, 1, 2])

    def __init_attributes(self):
        """
        Creates storage objects for the BDF object.
        This would be in the init but doing it this way allows for better
        inheritance

        References:
          1.  http://www.mscsoftware.com/support/library/conf/wuc87/p02387.pdf
        """
        self.bdf_filename = None
        self.punch = None
        self._encoding = None

        #: list of execive control deck lines
        self.executive_control_lines = []

        #: list of case control deck lines
        self.case_control_lines = []

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
        self.values_to_skip = {}

        # ------------------------ structural defaults -----------------------
        #: the analysis type
        self._sol = None
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
        self.rigid_elements = {}
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
        self.thermal_materials = {}

        #: defines the MATHE, MATHP
        self.hyperelastic_materials = {}

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
        self.creep_materials = {}

        # loads
        #: stores LOAD, FORCE, FORCE1, FORCE2, MOMENT, MOMENT1, MOMENT2,
        #: PLOAD, PLOAD2, PLOAD4, SLOAD
        #: GMLOAD, SPCD,
        #: QVOL
        self.loads = {}
        self.tics = {}

        # stores DLOAD entries.
        self.dloads = {}
        # stores ACSRCE, RLOAD1, RLOAD2, TLOAD1, TLOAD2, and ACSRCE entries.
        self.dload_entries = {}

        #self.gusts  = {} # Case Control GUST = 100
        #self.random = {} # Case Control RANDOM = 100

        #: stores coordinate systems
        origin = array([0., 0., 0.])
        zaxis = array([0., 0., 1.])
        xzplane = array([1., 0., 0.])
        coord = CORD2R(cid=0, rid=0, origin=origin, zaxis=zaxis, xzplane=xzplane)
        self.coords = {0 : coord}

        # --------------------------- constraints ----------------------------
        #: stores SUPORT1s
        #self.constraints = {} # suport1, anything else???
        self.suport = []
        self.suport1 = {}
        self.se_suport = []

        #: stores SPCADD, SPC, SPC1, SPCAX, GMSPC
        #self.spcObject = ConstraintObject()
        #: stores MPCADD,MPC
        #self.mpcObject = ConstraintObject()

        self.spcs = {}
        self.spcadds = {}

        self.mpcs = {}
        self.mpcadds = {}

        # --------------------------- dynamic ----------------------------
        #: stores DAREA
        self.dareas = {}
        self.dphases = {}

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
                'CTETRA', 'CPYRAM', 'CPENTA', 'CHEXA', 'CIHEX1',
                'CSHEAR', 'CVISC', 'CRAC2D', 'CRAC3D',
                'CGAP',

                # thermal
                'CHBDYE', 'CHBDYG', 'CHBDYP',
            ],
            'rigid_elements' : ['RBAR', 'RBAR1', 'RBE1', 'RBE2', 'RBE3', 'RROD', 'RSPLINE'],
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
                'PIHEX', 'PCOMPS',
            ],
            'pdampt' : ['PDAMPT',],
            'pelast' : ['PELAST',],
            'pbusht' : ['PBUSHT',],

            # materials
            'materials' : ['MAT1', 'MAT2', 'MAT3', 'MAT8', 'MAT9', 'MAT10', 'MAT11'],
            'hyperelastic_materials' : ['MATHP',],
            'creep_materials' : ['CREEP'],
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
            'thermal_materials' : ['MAT4', 'MAT5',],

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
                'PLOADX1', 'RFORCE', 'RFORCE1', 'SLOAD',
                'GMLOAD', 'SPCD', 'LOADCYN',

                # thermal
                'TEMP', 'QBDY1', 'QBDY2', 'QBDY3', 'QHBDY',
                'QVOL',
                ],
            'dloads' : ['DLOAD', ],
            # stores RLOAD1, RLOAD2, TLOAD1, TLOAD2, and ACSRCE entries.
            'dload_entries' : ['ACSRCE', 'TLOAD1', 'TLOAD2', 'RLOAD1', 'RLOAD2',],

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
            'aesurf' : ['AESURF'],
            'aesurfs' : ['AESURFS'],
            'aestats' : ['AESTAT'],
            'caeros' : ['CAERO1', 'CAERO2', 'CAERO3', 'CAERO4', 'CAERO5'],
            'paeros' : ['PAERO1', 'PAERO2', 'PAERO3', 'PAERO4', 'PAERO5'],
            'monitor_points' : ['MONPNT1'],
            'splines' : ['SPLINE1', 'SPLINE2', 'SPLINE4', 'SPLINE5',],
            'csschds' : ['CSSCHD',],
            #'SPLINE3', 'SPLINE6', 'SPLINE7',
            'trims' : ['TRIM',],
            'divergs' : ['DIVERG'],

            # coords
            'coords' : ['CORD1R', 'CORD1C', 'CORD1S',
                        'CORD2R', 'CORD2C', 'CORD2S',
                        'GMCORD'],

            # temperature cards
            'tempds' : ['TEMPD'],

            'phbdys' : ['PHBDY'],
            'convection_properties' : ['PCONV', 'PCONVM'],

            # stores thermal boundary conditions
            'bcs' : ['CONV', 'RADBC', 'RADM'],


            # dynamic cards
            'dareas' : ['DAREA'],
            'dphases' : ['DPHASE'],
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
            'dconstrs' : ['DCONSTR', 'DCONADD'],
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
            'random_tables' : ['TABRND1', 'TABRNDG',],

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

    #@property
    #def caseControlDeck(self):
        #self.deprecated('self.caseControlDeck', 'self.case_control_deck', '0.8')
        #return self.case_control_deck

    #@caseControlDeck.setter
    #def caseControlDeck(self, value):
        #self.deprecated('self.caseControlDeck', 'self.case_control_deck', '0.8')
        #self.case_control_deck = value

    @property
    def sol(self):
        """gets the solution (e.g. 101, 103)"""
        return self._sol

    @sol.setter
    def sol(self, sol):
        """sets the solution (e.g. 101, 103)"""
        self._sol = sol
        if len(self.executive_control_lines) == 0:
            self.executive_control_lines = ['SOL %s' % sol, 'CEND']
            self.iSolLine = 0
        return self._sol

    @property
    def subcases(self):
        """gets the subcases"""
        if self.case_control_deck is None:
            return {}
        return self.case_control_deck.subcases

    @property
    def rejects(self):
        #: lines that were rejected b/c they were for a card that isnt supported
        return self.reject_lines

    @rejects.setter
    def rejects(self, rejects):
        self.reject_lines = rejects

    @property
    def nnodes(self):
        """gets the number of GRIDs"""
        return len(self.nodes)

    @property
    def node_ids(self):
        """gets the GRID ids"""
        return self.nodes.keys()

    @property
    def point_ids(self):
        """gets the GRID, SPOINT, EPOINT ids"""
        if self.spoints is not None and self.epoints is not None:
            return set(self.node_ids) | self.spoints.points | self.epoints.points
        elif self.spoints is not None:
            return set(self.node_ids) | self.spoints.points
        elif self.epoints is not None:
            return set(self.node_ids) | self.epoints.points
        return set(self.node_ids)

    #--------------------
    # Elements CARDS

    @property
    def nelements(self):
        """gets the number of element"""
        return len(self.elements)

    @property
    def element_ids(self):
        """gets the element ids"""
        return self.elements.keys()

    #--------------------
    # Property CARDS

    @property
    def nproperties(self):
        """gets the number of properties"""
        return len(self.properties)

    @property
    def property_ids(self):
        """gets the property ids"""
        return self.properties.keys()

    #--------------------
    # Material CARDS

    @property
    def nmaterials(self):
        """gets the number of materials"""
        return len(self.materials)

    @property
    def material_ids(self):
        """gets the material ids"""
        return self.materials.keys()

    #--------------------
    # Coords CARDS

    @property
    def ncoords(self):
        """gets the number of coordinate systems"""
        return len(self.coords)

    @property
    def coord_ids(self):
        """gets the number of coordinate system ids"""
        return self.coords.keys()

    #--------------------

    @property
    def ncaeros(self):
        """gets the number of CAEROx panels"""
        return len(self.caeros)

    @property
    def caero_ids(self):
        """gets the CAEROx ids"""
        return self.caeros.keys()
