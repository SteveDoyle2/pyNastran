"""defines the BDF attributes"""
from __future__ import print_function, unicode_literals
from collections import defaultdict
from typing import List, Dict, Optional, Any
from numpy import array  # type: ignore

from pyNastran.utils import object_attributes, object_methods
from pyNastran.bdf.bdf_interface.utils import deprecated
#from pyNastran.bdf.case_control_deck import CaseControlDeck
from pyNastran.bdf.cards.coordinate_systems import CORD2R
#from pyNastran.bdf.cards.constraints import ConstraintObject
from pyNastran.bdf.cards.aero.zona import ZONA


class BDFAttributes(object):
    """defines attributes of the BDF"""

    def __init__(self):
        """creates the attributes for the BDF"""
        self.__init_attributes()
        self._is_cards_dict = False

        self.is_nx = False
        self.is_msc = False
        self.is_zona = False
        self.save_file_structure = False
        self.is_superelements = False
        self.set_as_msc()
        self.units = []  # type: List[str]

    def set_as_msc(self):
        self._nastran_format = 'msc'
        self.is_nx = False
        self.is_msc = True
        self.is_zona = False

    def set_as_nx(self):
        self._nastran_format = 'nx'
        self.is_nx = True
        self.is_msc = False
        self.is_zona = False

    def set_as_zona(self):
        self._nastran_format = 'zona'
        self.is_nx = False
        self.is_msc = False
        self.is_zona = True

    def __properties__(self):
        """the list of @property attributes"""
        return ['nastran_format', 'is_long_ids', 'sol', 'subcases',
                'nnodes', 'node_ids', 'point_ids', 'npoints',
                'nelements', 'element_ids', 'nproperties', 'property_ids',
                'nmaterials', 'material_ids', 'ncoords', 'coord_ids',
                'ncaeros', 'caero_ids', 'wtmass', 'is_bdf_vectorized', 'nid_map']

    def object_attributes(self, mode='public', keys_to_skip=None):
        # type: (str, Optional[List[str]]) -> List[str]
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
            '_card_parser', '_card_parser_b', '_card_parser_prepare',
            'object_methods', 'object_attributes',
        ]
        return object_attributes(self, mode=mode, keys_to_skip=keys_to_skip+my_keys_to_skip)

    def object_methods(self, mode='public', keys_to_skip=None):
        # type: (str, Optional[List[str]]) -> List[str]
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
        my_keys_to_skip = []  # type: List[str]

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
        # type: (str, str, str) -> None
        """deprecates methods"""
        return deprecated(old_name, new_name, deprecated_version, levels=[0, 1, 2])

    def clear_attributes(self):
        # type: () -> None
        """removes the attributes from the model"""
        self.__init_attributes()

        self.nodes = {}
        self.loads = {}  # type: Dict[int, List[Any]]
        self.load_combinations = {}  # type: Dict[int, List[Any]]

    def reset_errors(self):
        """removes the errors from the model"""
        self._ixref_errors = 0
        self._stored_xref_errors = []

    def __init_attributes(self):
        # type: () -> None
        """
        Creates storage objects for the BDF object.
        This would be in the init but doing it this way allows for better
        inheritance

        References:
          1.  http://www.mscsoftware.com/support/library/conf/wuc87/p02387.pdf
        """
        self.reset_errors()
        self.bdf_filename = None
        self.punch = None
        self._encoding = None
        self._is_long_ids = False # ids > 8 characters

        #: ignore any ECHOON flags
        self.force_echo_off = True

        #: list of Nastran SYSTEM commands
        self.system_command_lines = []  # type: List[str]

        #: list of execive control deck lines
        self.executive_control_lines = []  # type: List[str]

        #: list of case control deck lines
        self.case_control_lines = []  # type: List[str]

        # dictionary of BDFs
        self.superelement_models = {}
        self.initial_superelement_models = []  # the keys before superelement mirroring

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

        # ------------------------ bad duplicates ----------------------------
        self._iparse_errors = 0
        self._nparse_errors = 0
        self._stop_on_parsing_error = True
        self._stop_on_duplicate_error = True
        self._stored_parse_errors = []  # type: List[str]

        self._duplicate_nodes = []  # type: List[str]
        self._duplicate_elements = []  # type: List[str]
        self._duplicate_properties = []  # type: List[str]
        self._duplicate_materials = []  # type: List[str]
        self._duplicate_masses = []  # type: List[str]
        self._duplicate_thermal_materials = []  # type: List[str]
        self._duplicate_coords = []  # type: List[str]
        self.values_to_skip = {}  # type: Dict[str, List[int]]

        # ------------------------ structural defaults -----------------------
        #: the analysis type
        self._sol = None
        #: used in solution 600, method
        self.sol_method = None
        #: the line with SOL on it, marks ???
        self.sol_iline = None  # type : Optional[int]
        self.case_control_deck = None  # type: Optional[Any]

        #: store the PARAM cards
        self.params = {}  # type: Dict[str, Any]
        # ------------------------------- nodes -------------------------------
        # main structural block
        #: stores POINT cards
        self.points = {}  # type: Dict[int, Any]
        #self.grids = {}

        self.spoints = {}  # type: Dict[int, Any]
        self.epoints = {}  # type: Dict[int, Any]

        #: stores GRIDSET card
        self.grdset = None  # type: Optional[Any]

        #: stores SEQGP cards
        self.seqgp = None  # type: Optional[Any]

        ## stores RINGAX
        self.ringaxs = {}  # type: Dict[int, Any]

        ## stores GRIDB
        self.gridb = {}  # type: Dict[int, Any]

        #: stores elements (CQUAD4, CTRIA3, CHEXA8, CTETRA4, CROD, CONROD,
        #: etc.)
        self.elements = {}  # type: Dict[int, Any]

        #: stores CBARAO, CBEAMAO
        self.ao_element_flags = {}  # type: Dict[int, Any]
        #: stores BAROR
        self.baror = None  # type: Optional[Any]
        #: stores BEAMOR
        self.beamor = None  # type: Optional[Any]
        #: stores SNORM
        self.normals = {}  # type: Dict[int, Any]

        #: stores rigid elements (RBE2, RBE3, RJOINT, etc.)
        self.rigid_elements = {}  # type: Dict[int, Any]
        #: stores PLOTELs
        self.plotels = {}  # type: Optional[Any]

        #: stores CONM1, CONM2, CMASS1,CMASS2, CMASS3, CMASS4, CMASS5
        self.masses = {}  # type: Dict[int, Any]
        #: stores PMASS
        self.properties_mass = {}  # type: Dict[int, Any]

        #: stores NSM, NSM1, NSML, NSML1
        self.nsms = {}  # type: Dict[int, List[Any]]
        #: stores NSMADD
        self.nsmadds = {}  # type: Dict[int, List[Any]]

        #: stores LOTS of propeties (PBAR, PBEAM, PSHELL, PCOMP, etc.)
        self.properties = {}  # type: Dict[int, Any]

        #: stores MAT1, MAT2, MAT3, MAT8, MAT10, MAT11
        self.materials = {}  # type: Dict[int, Any]

        #: defines the MAT4, MAT5
        self.thermal_materials = {}  # type: Dict[int, Any]

        #: defines the MATHE, MATHP
        self.hyperelastic_materials = {}  # type: Dict[int, Any]

        #: stores MATSx
        self.MATS1 = {}  # type: Dict[int, Any]
        self.MATS3 = {}  # type: Dict[int, Any]
        self.MATS8 = {}  # type: Dict[int, Any]

        #: stores MATTx
        self.MATT1 = {}  # type: Dict[int, Any]
        self.MATT2 = {}  # type: Dict[int, Any]
        self.MATT3 = {}  # type: Dict[int, Any]
        self.MATT4 = {}  # type: Dict[int, Any]
        self.MATT5 = {}  # type: Dict[int, Any]
        self.MATT8 = {}  # type: Dict[int, Any]
        self.MATT9 = {}  # type: Dict[int, Any]
        self.nxstrats = {}  # type: Dict[int, Any]

        #: stores the CREEP card
        self.creep_materials = {}  # type: Dict[int, Any]

        self.tics = {}  # type: Optional[Any]

        # stores DLOAD entries.
        self.dloads = {}    # type: Dict[int, Any]
        # stores ACSRCE, RLOAD1, RLOAD2, TLOAD1, TLOAD2, and ACSRCE,
        #        and QVECT entries.
        self.dload_entries = {}    # type: Dict[int, Any]

        #self.gusts = {}  # Case Control GUST = 100
        #self.random = {} # Case Control RANDOM = 100

        #: stores coordinate systems
        origin = array([0., 0., 0.])
        zaxis = array([0., 0., 1.])
        xzplane = array([1., 0., 0.])
        coord = CORD2R(cid=0, rid=0, origin=origin, zaxis=zaxis, xzplane=xzplane)
        self.coords = {0 : coord}   # type: Dict[int, Any]

        # --------------------------- constraints ----------------------------
        #: stores SUPORT1s
        #self.constraints = {} # suport1, anything else???
        self.suport = []  # type: List[Any]
        self.suport1 = {}  # type: Dict[int, Any]
        self.se_suport = []  # type: List[Any]

        #: stores SPC, SPC1, SPCAX, GMSPC
        self.spcs = {}  # type: Dict[int, List[Any]]
        #: stores SPCADD
        self.spcadds = {}  # type: Dict[int, List[Any]]
        self.spcoffs = {}  # type: Dict[int, List[Any]]

        self.mpcs = {}  # type: Dict[int, List[Any]]
        self.mpcadds = {}  # type: Dict[int, List[Any]]

        # --------------------------- dynamic ----------------------------
        #: stores DAREA
        self.dareas = {}   # type: Dict[int, Any]
        self.dphases = {}  # type: Dict[int, Any]

        self.pbusht = {}  # type: Dict[int, Any]
        self.pdampt = {}  # type: Dict[int, Any]
        self.pelast = {}  # type: Dict[int, Any]

        #: frequencies
        self.frequencies = {}  # type: Dict[int, List[Any]]

        # ----------------------------------------------------------------
        #: direct matrix input - DMIG
        self.dmis = {}  # type: Dict[str, Any]
        self.dmigs = {}  # type: Dict[str, Any]
        self.dmijs = {}  # type: Dict[str, Any]
        self.dmijis = {}  # type: Dict[str, Any]
        self.dmiks = {}  # type: Dict[str, Any]
        self.dti = {}  # type: Dict[str, Any]
        self._dmig_temp = defaultdict(list)  # type: Dict[str, List[str]]

        # ----------------------------------------------------------------
        #: SETy
        self.sets = {}  # type: Dict[int, Any]
        self.asets = []  # type: List[Any]
        self.omits = []  # type: List[Any]
        self.bsets = []  # type: List[Any]
        self.csets = []  # type: List[Any]
        self.qsets = []  # type: List[Any]
        self.usets = {}  # type: Dict[str, Any]

        #: SExSETy
        self.se_bsets = []  # type: List[Any]
        self.se_csets = []  # type: List[Any]
        self.se_qsets = []  # type: List[Any]
        self.se_usets = {}  # type: Dict[str, Any]
        self.se_sets = {}  # type: Dict[str, Any]

        # ----------------------------------------------------------------
        #: tables
        # TABLES1, ...
        self.tables = {}  # type: Dict[int, Any]

        # TABLEDx
        self.tables_d = {}  # type: Dict[int, Any]

        # TABLEMx
        self.tables_m = {}  # type: Dict[int, Any]

        #: random_tables
        self.random_tables = {}  # type: Dict[int, Any]
        #: TABDMP1
        self.tables_sdamping = {}  # type: Dict[int, Any]

        # ----------------------------------------------------------------
        #: EIGB, EIGR, EIGRL methods
        self.methods = {}  # type: Dict[int, Any]
        # EIGC, EIGP methods
        self.cMethods = {}  # type: Dict[int, Any]

        # ---------------------------- optimization --------------------------
        # optimization
        self.dconadds = {}  # type: Dict[int, Any]
        self.dconstrs = {}  # type: Dict[int, Any]
        self.desvars = {}  # type: Dict[int, Any]
        self.ddvals = {}  # type: Dict[int, Any]
        self.dlinks = {}  # type: Dict[int, Any]
        self.dresps = {}  # type: Dict[int, Any]

        self.dtable = None  # type: Optional[Any]
        self.dequations = {}  # type: Dict[int, Any]

        #: stores DVPREL1, DVPREL2...might change to DVxRel
        self.dvprels = {}  # type: Dict[int, Any]
        self.dvmrels = {}  # type: Dict[int, Any]
        self.dvcrels = {}  # type: Dict[int, Any]
        self.dvgrids = {}  # type: Dict[int, Any]
        self.doptprm = None  # type: Optional[Any]
        self.dscreen = {}  # type: Dict[int, Any]

        # ------------------------- nonlinear defaults -----------------------
        #: stores NLPCI
        self.nlpcis = {}  # type: Dict[int, Any]
        #: stores NLPARM
        self.nlparms = {}  # type: Dict[int, Any]
        #: stores TSTEPs, TSTEP1s
        self.tsteps = {}  # type: Dict[int, Any]
        #: stores TSTEPNL
        self.tstepnls = {}  # type: Dict[int, Any]
        #: stores TF
        self.transfer_functions = {}  # type: Dict[int, Any]
        #: stores DELAY
        self.delays = {}  # type: Dict[int, Any]

        #: stores ROTORG
        self.rotors = {}  # type: Dict[int, Any]

        # --------------------------- aero defaults --------------------------
        # aero cards
        #: stores CAEROx
        self.caeros = {}  # type: Dict[int, Any]
        #: stores PAEROx
        self.paeros = {}  # type: Dict[int, Any]
        # stores MONPNT1
        self.monitor_points = []  # type: List[Any]

        #: stores AECOMP
        self.aecomps = {}  # type: Dict[int, Any]
        #: stores AEFACT
        self.aefacts = {}  # type: Dict[int, Any]
        #: stores AELINK
        self.aelinks = {}  # type: Dict[int, List[Any]]
        #: stores AELIST
        self.aelists = {}  # type: Dict[int, Any]
        #: stores AEPARAM
        self.aeparams = {}  # type: Dict[int, Any]
        #: stores AESURF
        self.aesurf = {}  # type: Dict[int, Any]
        #: stores AESURFS
        self.aesurfs = {}  # type: Dict[int, Any]
        #: stores AESTAT
        self.aestats = {}  # type: Dict[int, Any]
        #: stores CSSCHD
        self.csschds = {}  # type: Dict[int, Any]

        #: store SPLINE1,SPLINE2,SPLINE4,SPLINE5
        self.splines = {}  # type: Dict[int, Any]
        self.zona = ZONA(self)

        # axisymmetric
        self.axic = None  # type: Optional[Any]
        self.axif = None  # type: Optional[Any]
        self.ringfl = {}  # type: Dict[int, Any]
        self._is_axis_symmetric = False

        # ------ SOL 144 ------
        #: stores AEROS
        self.aeros = None  # type: Optional[Any]

        #: stores TRIM
        self.trims = {}  # type: Dict[int, Any]

        #: stores DIVERG
        self.divergs = {}  # type: Dict[int, Any]

        # ------ SOL 145 ------
        #: stores AERO
        self.aero = None  # type: Optional[Any]

        #: stores FLFACT
        #: .. todo:: can this be simplified ???
        self.flfacts = {}  # type: Dict[int, Any]

        #: stores FLUTTER
        self.flutters = {} # type: Dict[int, Any]

        #: mkaeros
        self.mkaeros = []  # type: List[Any]

        # ------ SOL 146 ------
        #: stores GUST cards
        self.gusts = {}  # type: Dict[int, Any]

        # ------------------------- thermal defaults -------------------------
        # BCs
        #: stores thermal boundary conditions - CONV,RADBC
        self.bcs = {}  # type: Dict[int, Any]

        #: stores PHBDY
        self.phbdys = {}  # type: Dict[int, Any]
        #: stores convection properties - PCONV, PCONVM ???
        self.convection_properties = {}  # type: Dict[int, Any]
        #: stores TEMPD
        self.tempds = {}  # type: Dict[int, Any]

        #: stores VIEW
        self.views = {}
        #: stores VIEW3D
        self.view3ds = {}
        self.radset = None
        self.radcavs = {}
        self.radmtx = {}

        # -------------------------contact cards-------------------------------
        self.bcrparas = {}  # type: Dict[int, Any]
        self.bctadds = {}  # type: Dict[int, Any]
        self.bctparas = {}  # type: Dict[int, Any]
        self.bctsets = {}  # type: Dict[int, Any]
        self.bsurf = {}  # type: Dict[int, Any]
        self.bsurfs = {}  # type: Dict[int, Any]
        self.bconp = {}  # type: Dict[int, Any]
        self.blseg = {}  # type: Dict[int, Any]


        #--------------------------superelements------------------------------
        self.setree = {}  # type: Dict[int, Any]
        self.senqset = {}  # type: Dict[int, Any]
        self.sebulk = {}  # type: Dict[int, Any]
        self.sebndry = {}  # type: Dict[int, Any]
        self.seloc = {}  # type: Dict[int, Any]
        self.sempln = {}  # type: Dict[int, Any]
        self.seconct = {}  # type: Dict[int, Any]
        self.selabel = {}  # type: Dict[int, Any]
        self.seexcld = {}  # type: Dict[int, Any]
        self.seelt = {}  # type: Dict[int, Any]
        self.seload = {}  # type: Dict[int, Any]
        self.csuper = {}  # type: Dict[int, Any]
        self.csupext = {}  # type: Dict[int, Any]

        # ---------------------------------------------------------------------
        self._type_to_id_map = defaultdict(list)  # type: Dict[int, List[Any]]
        self._slot_to_type_map = {
            'params' : ['PARAM'],
            'nodes' : ['GRID', 'SPOINT', 'EPOINT'], # 'RINGAX',
            'points' : ['POINT'],
            'ringaxs' : ['RINGAX', 'POINTAX'],
            'ringfl' : ['RINGFL'],
            'axic' : ['AXIC'],
            'axif' : ['AXIF'],
            'grdset' : ['GRDSET'],
            'gridb' : ['GRIDB'],
            'seqgp' : ['SEQGP'],
            'ao_element_flags' : ['CBARAO'],
            #'POINTAX', 'RINGAX',

            # CMASS4 lies in the QRG
            'masses' : ['CONM1', 'CONM2', 'CMASS1', 'CMASS2', 'CMASS3', 'CMASS4'],

            'elements' : [
                'CELAS1', 'CELAS2', 'CELAS3', 'CELAS4',
                # 'CELAS5',
                'CBUSH', 'CBUSH1D', 'CBUSH2D',

                'CDAMP1', 'CDAMP2', 'CDAMP3', 'CDAMP4', 'CDAMP5',
                'CFAST', 'GENEL',

                'CBAR', 'CROD', 'CTUBE', 'CBEAM', 'CBEAM3', 'CONROD', 'CBEND',
                'CTRIA3', 'CTRIA6', 'CTRIAR',
                'CQUAD4', 'CQUAD8', 'CQUADR', 'CQUAD',
                'CPLSTN3', 'CPLSTN6', 'CPLSTN4', 'CPLSTN8',
                'CPLSTS3', 'CPLSTS6', 'CPLSTS4', 'CPLSTS8',
                'CTRAX3', 'CTRAX6', 'CTRIAX', 'CTRIAX6',
                'CQUADX', 'CQUADX4', 'CQUADX8',

                'CTETRA', 'CPYRAM', 'CPENTA', 'CHEXA', 'CIHEX1', 'CIHEX2',
                'CSHEAR', 'CVISC', 'CRAC2D', 'CRAC3D',
                'CGAP',

                # thermal
                'CHBDYE', 'CHBDYG', 'CHBDYP',
            ],
            'normals' : ['SNORM'],
            'nsms' : ['NSM', 'NSM1', 'NSML', 'NSML1'],
            'nsmadds' : ['NSMADD'],
            'rigid_elements' : ['RBAR', 'RBAR1', 'RBE1', 'RBE2', 'RBE3', 'RROD', 'RSPLINE', 'RSSCON'],
            'plotels' : ['PLOTEL',],

            'properties_mass' : ['PMASS'],
            'properties' : [
                'PELAS', 'PGAP', 'PFAST', 'PLPLANE', 'PPLANE',
                'PBUSH', 'PBUSH1D',
                'PDAMP', 'PDAMP5',
                'PROD', 'PBAR', 'PBARL', 'PBEAM', 'PTUBE', 'PBEND', 'PBCOMP', 'PBRSECT', 'PBMSECT',
                'PBEAML',  # not fully supported
                # 'PBEAM3',

                'PSHELL', 'PCOMP', 'PCOMPG', 'PSHEAR',
                'PSOLID', 'PLSOLID', 'PVISC', 'PRAC2D', 'PRAC3D',
                'PIHEX', 'PCOMPS',
                'PCONEAX',
            ],
            'pdampt' : ['PDAMPT',],
            'pelast' : ['PELAST',],
            'pbusht' : ['PBUSHT',],

            # materials
            'materials' : ['MAT1', 'MAT2', 'MAT3', 'MAT8', 'MAT9', 'MAT10', 'MAT11',
                           'MAT3D', 'MATG'],
            'hyperelastic_materials' : ['MATHE', 'MATHP',],
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
            'nxstrats' : ['NXSTRAT'],

            # 'MATHE'
            #'EQUIV', # testing only, should never be activated...

            # thermal materials
            'thermal_materials' : ['MAT4', 'MAT5',],

            # spc/mpc constraints - TODO: is this correct?
            'spcadds' : ['SPCADD'],
            'spcs' : ['SPC', 'SPC1', 'SPCAX', 'GMSPC'],
            'spcoffs' : ['SPCOFF', 'SPCOFF1'],
            'mpcadds' : ['MPCADD'],
            'mpcs' : ['MPC'],
            'suport' : ['SUPORT'],
            'suport1' : ['SUPORT1'],
            'se_suport' : ['SESUP'],

            'setree' : ['SETREE'],
            'senqset' : ['SENQSET'],
            'sebulk' : ['SEBULK'],
            'sebndry' : ['SEBNDRY'],
            'seloc' : ['SELOC'],
            'sempln' : ['SEMPLN'],
            'seconct' : ['SECONCT'],
            'selabel' : ['SELABEL'],
            'seexcld' : ['SEEXCLD'],
            'seelt' : ['SEELT'],
            'seload' : ['SELOAD'],
            'csuper' : ['CSUPER'],
            'csupext' : ['CSUPEXT'],

            # loads
            'load_combinations' : ['LOAD', 'LSEQ'],
            'loads' : [
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
            'dload_entries' : ['ACSRCE', 'TLOAD1', 'TLOAD2', 'RLOAD1', 'RLOAD2',
                               'QVECT', 'RANDPS', 'RANDT1'],

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
            'caeros' : ['CAERO1', 'CAERO2', 'CAERO3', 'CAERO4', 'CAERO5', 'CAERO7', 'BODY7'],
            'paeros' : ['PAERO1', 'PAERO2', 'PAERO3', 'PAERO4', 'PAERO5', 'SEGMESH'],
            'monitor_points' : ['MONPNT1', 'MONPNT2', 'MONPNT3'],
            'splines' : ['SPLINE1', 'SPLINE2', 'SPLINE3', 'SPLINE4', 'SPLINE5', 'SPLINE6', 'SPLINE7'],
            'panlsts' : ['PANLST1', 'PANLST2', 'PANLST3',],
            'csschds' : ['CSSCHD',],
            #'SPLINE3', 'SPLINE6', 'SPLINE7',
            'trims' : ['TRIM', 'TRIM2'],
            'divergs' : ['DIVERG'],

            # coords
            'coords' : ['CORD1R', 'CORD1C', 'CORD1S',
                        'CORD2R', 'CORD2C', 'CORD2S',
                        'GMCORD', 'ACOORD', 'CORD3G'],

            # temperature cards
            'tempds' : ['TEMPD'],

            'phbdys' : ['PHBDY'],
            'convection_properties' : ['PCONV', 'PCONVM'],

            # stores thermal boundary conditions
            'bcs' : ['CONV', 'CONVM', 'RADBC', 'RADM', 'TEMPBC'],


            # dynamic cards
            'dareas' : ['DAREA'],
            'tics' : ['TIC'],
            'dphases' : ['DPHASE'],
            'nlparms' : ['NLPARM'],
            'nlpcis' : ['NLPCI'],
            'tsteps' : ['TSTEP'],
            'tstepnls' : ['TSTEPNL', 'TSTEP1'],
            'transfer_functions' : ['TF'],
            'delays' : ['DELAY'],
            'rotors' : ['ROTORG', 'ROTORD'],

            'frequencies' : ['FREQ', 'FREQ1', 'FREQ2', 'FREQ3', 'FREQ4', 'FREQ5'],

            # direct matrix input cards
            'dmigs' : ['DMIG'],
            'dmijs' : ['DMIJ'],
            'dmijis' : ['DMIJI'],
            'dmiks' : ['DMIK'],
            'dmis' : ['DMI'],
            'dti' : ['DTI'],

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
            'dvgrids' : ['DVGRID'],
            'doptprm' : ['DOPTPRM'],
            'dscreen' : ['DSCREEN'],


            # sets
            'asets' : ['ASET', 'ASET1'],
            'omits' : [
                #'OMIT',
                'OMIT1'
            ],
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
            'radset' : ['RADSET'],
            'radcavs' : ['RADCAV', 'RADLST'],
            'radmtx' : ['RADMTX'],
            # SEBSEP

            'tables' : [
                'TABLEH1', 'TABLEHT',
                'TABLES1', 'TABLEST',
                ],
            'tables_d' : ['TABLED1', 'TABLED2', 'TABLED3', 'TABLED4', 'TABLED5'],
            'tables_m' : ['TABLEM1', 'TABLEM2', 'TABLEM3', 'TABLEM4'],
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
            'bconp' : ['BCONP'],
            'blseg' : ['BLSEG'],
            'views' : ['VIEW'],
            'view3ds' : ['VIEW3D'],

            ## other
            #'INCLUDE',  # '='
            #'ENDDATA',
        }  # type: Dict[str, List[str]]
        self._type_to_slot_map = self.get_rslot_map()

    @property
    def type_slot_str(self):
        """helper method for printing supported cards"""
        nchars = len('Card Group')

        nchars_cards = 0
        for card_group in self._slot_to_type_map:
            nchars = max(nchars, len(card_group))

        nline = 58
        fmt =     '| %%-%ss | %%-%ss |\n' % (nchars, nline)
        fmt_plus = '+%%-%ss+%%-%ss+\n' % (nchars + 2, nline + 2)

        dash1 = '-' * (nchars + 2)
        dash2 = '-' * (nline + 2)
        dash_plus = fmt_plus % (dash1, dash2)
        html_msg = [
            dash_plus,
            fmt % ('Card Group', 'Cards'),
        ]
        for card_group, card_types in sorted(self._slot_to_type_map.items()):
            valid_cards = [card_type for card_type in card_types
                           if card_type in self.cards_to_read]
            valid_cards.sort()
            if len(valid_cards) == 0:
                continue

            i = 0
            sublines = []
            subline = ''
            while valid_cards:
                card_type = valid_cards.pop(0)

                # the +2 is for the comma and space
                len_card_type = len(card_type) + 2
                nline_new = len(subline) + len_card_type
                if nline_new > nline:
                    sublines.append(subline.rstrip(' '))
                    subline = ''
                subline += '%s, ' % card_type

            if subline:
                sublines.append(subline.rstrip(', '))

            html_msg.append(dash_plus)
            for subline in sublines:
                html_msg.append(fmt % (card_group, subline))
                card_group = ''

        html_msg.append(dash_plus)

        #for card_group, card_types in sorted(self._slot_to_type_map.items()):
            #html_msg.append('| %s | %s |' % (card_group, ', '.join(card_types)))

        #html_msg.append(
            #fmt_plus % ('-'*(nchars + 2), '-'*(nline + 2))
        #)
        msg = ''.join(html_msg)
        return msg

    @property
    def nastran_format(self):
        # type: () -> str
        return self._nastran_format

    @nastran_format.setter
    def nastran_format(self, nastran_format):
        # type: (str) -> None
        fmt_lower = nastran_format.lower().strip()
        if fmt_lower not in ['nx', 'msc', 'zona']:
            raise RuntimeError(nastran_format)
        self._nastran_format = fmt_lower

    @property
    def is_long_ids(self):
        # type: () -> bool
        return self._is_long_ids
        #if self._nastran_format == 'nx' or self._is_long_ids:
            #return True
        #return False

    @property
    def sol(self):
        # type: () -> int
        """gets the solution (e.g. 101, 103)"""
        return self._sol

    @sol.setter
    def sol(self, sol):
        # type: (int) -> int
        """sets the solution (e.g. 101, 103)"""
        self._sol = sol
        if len(self.executive_control_lines) == 0:
            self.executive_control_lines = ['SOL %s' % sol, 'CEND']
            self.sol_iline = 0
        return self._sol

    @property
    def subcases(self):
        # type: () -> Dict[int, Optional[Any]]
        """gets the subcases"""
        if self.case_control_deck is None:
            return {}
        return self.case_control_deck.subcases

    #@property
    #def grids(self):
        #"""might be renaming self.nodes to self.grids"""
        #return self.nodes

    #@property.setter
    #def grids(self, grids):
        #"""might be renaming self.nodes to self.grids"""
        #self.nodes = grids

    @property
    def nnodes(self):
        # type: () -> int
        """gets the number of GRIDs"""
        return len(self.nodes)

    @property
    def node_ids(self):
        """gets the GRID ids"""
        return self.nodes.keys()

    @property
    def point_ids(self):
        """gets the GRID, SPOINT, EPOINT ids"""
        return set(self.node_ids) | set(list(self.spoints.keys())) | set(list(self.epoints.keys()))

    @property
    def npoints(self):
        # type: () -> int
        """gets the number of GRID, SPOINT, EPOINT ids"""
        return len(self.point_ids)

    #--------------------
    # Elements CARDS

    @property
    def nelements(self):
        # type: () -> int
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
        # type: () -> int
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
        # type: () -> int
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
        # type: () -> int
        """gets the number of coordinate systems"""
        return len(self.coords)

    @property
    def coord_ids(self):
        """gets the number of coordinate system ids"""
        return self.coords.keys()

    #--------------------

    @property
    def ncaeros(self):
        # type: () -> int
        """gets the number of CAEROx panels"""
        return len(self.caeros)

    @property
    def caero_ids(self):
        """gets the CAEROx ids"""
        return self.caeros.keys()

    @property
    def wtmass(self):
        """
        Gets the PARAM,WTMASS value, which defines the weight to mass
        conversion factor

        kg -> kg : 1.0
        lb -> slug : 1/32.2
        lb -> slinch : 1/(32.2*12)=1/386.4
        """
        wtmass = 1.0
        if 'WTMASS' in self.params:
            param = self.params['WTMASS']
            wtmass = param.values[0]
        return wtmass

    def set_param(self, key, values):
        # type: (str, Union[int, float, str, List[float]]) -> None
        """sets a param card; creates it if necessary"""
        if isinstance(values, (int, float, str)):
            values = [values]
        key = key.upper()
        if key in self.params:
            param = self.params[key]
            param.update_values(self, *values)
