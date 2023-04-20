"""defines the BDF attributes"""
from __future__ import annotations
from collections import defaultdict
from typing import Optional, Any, Union, TYPE_CHECKING
from numpy import array
from pyNastran.bdf.cards.superelements import SEEXCLD  # type: ignore

from pyNastran.utils import object_attributes, object_methods, deprecated
#from pyNastran.bdf.case_control_deck import CaseControlDeck
from pyNastran.bdf.cards.coordinate_systems import CORD2R
#from pyNastran.bdf.cards.constraints import ConstraintObject
from pyNastran.bdf.cards.aero.zona import ZONA
if TYPE_CHECKING:  # pragma: no cover
    from pyNastran.bdf import (
        # BDF,
        CaseControlDeck,
        #params,
        PARAM, MDLPRM,
        # grids/points
        POINT, SPOINT, EPOINT,
        GRDSET, SEQGP, GRIDB,
        # bar
        BAROR, BEAMOR,
        # plot
        PLOTEL,
        # dynamic
        TSTEP, TSTEP1, TSTEPNL,
        NLPCI, NLPARM,
        TABLES1,
        TABLED1, TABLED2, TABLED3, TABLED4,
        TABLEM1, TABLEM2, TABLEM3, TABLEM4,
        TABDMP1,
        TF, DELAY, #DPHASE,
        # axisymmetric
        RINGAX, CYAX, AXIF, RINGFL, CYJOIN, AXIC,
        # shells
        SNORM,
        #CQUAD4, CQUAD8, CQUADR, CQUAD,
        #CTRIA3, CTRIA6, CTRIAR,
        # solids
        #CTETRA4, CPYRAM5, CPENTA6, CHEXA8,
        #CTETRA10, CPYRAM13, CPENTA15, CHEXA20,
        # loads
        TEMPD,
        # thermal
        #CHBYDP, CHBDYE, CHBDYP,
        PHBDY,
        CONV, PCONV, PCONVM, #CONVM,
        RADCAV, RADMTX, VIEW, VIEW3D,
        RADBC, #TEMPBC,
        # aero
        MONPNT1, MONPNT2, MONPNT3,
        AECOMP, AEFACT, AELINK, AELIST, AEPARM, AESURF, AESURFS, AESTAT,
        AERO, AEROS,
        CAERO1, CAERO2, CAERO3, CAERO4, CAERO5,
        PAERO1, PAERO2, PAERO3, PAERO4, PAERO5,
        SPLINE1, SPLINE2, SPLINE3, SPLINE4, SPLINE5,
        FLUTTER, MKAERO1, MKAERO2, FLFACT,
        TRIM, TRIM2, GUST, GUST2, DIVERG, CSSCHD,
        # roter
        ROTORD, ROTORG,
        # modal
        EIGRL, EIGR, EIGP, EIGC, EIGB, MODTRAK,
        # optimization
        DESVAR, DLINK, TOPVAR, DVGRID,
        DEQATN, DDVAL, DSCREEN,
        DTABLE, DRESP1, DRESP2, DRESP3,
        DVPREL1, DVCREL1, DVMREL1, DVTREL1,
        DVPREL2, DVCREL2, DVMREL2, DVTREL2,
        DCONADD, DCONSTR, DOPTPRM,
        DMNCON, GROUP,
        # contact
        BCPARA, BCBODY, BCTPARAM, BGSET, BCTADD, BSURF, BSURFS, BCONP, BLSEG, BFRIC,
        BCTSET, BCTPARA, BGADD, BCTPARA, BCRPARA,
        # superelements
        SEBULK, SENQSET, SENQSET1, SEBNDRY, RELEASE, SELOC, SEMPLN, SETREE,
        SELABEL, SECONCT, SEEXCLD, SEELT, SELOAD, CSUPER, CSUPEXT,
    )
    #Coord = Union[CORD1R, CORD1C, CORD1S,
    #              CORD2R, CORD2C, CORD2S]
    from pyNastran.bdf.cards.dmig import DMIG, DMI, DMIJ, DMIK, DMIJI
    from pyNastran.bdf.subcase import Subcase

BDF_FORMATS = {'nx', 'msc', 'optistruct', 'zona', 'nasa95', 'mystran'}

class BDFAttributes:
    """defines attributes of the BDF"""

    def __init__(self):
        """creates the attributes for the BDF"""
        self.__init_attributes()
        self._is_cards_dict = False

        self.is_nx = False
        self.is_msc = False
        self.is_mystran = False
        self.is_nasa95 = False
        self.is_zona = False
        self.save_file_structure = False
        self.is_superelements = False
        self.set_as_msc()
        self.units = []  # type: list[str]

    def set_as_msc(self):
        self._nastran_format = 'msc'
        self.is_nx = False
        self.is_msc = True
        self.is_optistruct = False
        self.is_mystran = False
        self.is_nasa95 = False
        self.is_zona = False

    def set_as_nx(self):
        self._nastran_format = 'nx'
        self.is_nx = True
        self.is_msc = False
        self.is_optistruct = False
        self.is_mystran = False
        self.is_nasa95 = False
        self.is_zona = False

    def set_as_optistruct(self):
        self._nastran_format = 'optistruct'
        self.is_nx = False
        self.is_msc = False
        self.is_optistruct = True
        self.is_mystran = False
        self.is_nasa95 = False
        self.is_zona = False

    def set_as_zona(self):
        self._nastran_format = 'zona'
        self.is_nx = False
        self.is_msc = False
        self.is_optistruct = False
        self.is_mystran = False
        self.is_nasa95 = False
        self.is_zona = True

    def set_as_mystran(self):
        self._nastran_format = 'mystran'
        self.is_nx = False
        self.is_msc = False
        self.is_optistruct = False
        self.is_mystran = True
        self.is_nasa95 = False
        self.is_zona = False
        self._update_for_mystran()

    def set_as_nasa95(self):
        self._nastran_format = 'nasa95'
        self.is_nx = False
        self.is_msc = False
        self.is_optistruct = False
        self.is_mystran = False
        self.is_nasa95 = True
        self.is_zona = False
        self._update_for_nasa95()

    def __properties__(self):
        """the list of @property attributes"""
        return ['nastran_format', 'is_long_ids', 'sol', 'subcases',
                'nnodes', 'node_ids', 'point_ids', 'npoints',
                'nelements', 'element_ids', 'nproperties', 'property_ids',
                'nmaterials', 'material_ids', 'ncoords', 'coord_ids',
                'ncaeros', 'caero_ids', 'wtmass', 'is_bdf_vectorized', 'nid_map']

    def object_attributes(self, mode: str='public',
                          keys_to_skip: Optional[list[str]]=None,
                          filter_properties: bool=False) -> list[str]:
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
        filter_properties: bool: default=False
            filters the @property objects

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
            'log',
            'node_ids', 'coord_ids', 'element_ids', 'property_ids',
            'material_ids', 'caero_ids', 'is_long_ids',
            'nnodes', 'ncoords', 'nelements', 'nproperties',
            'nmaterials', 'ncaeros', 'npoints',

            'point_ids', 'subcases',
            '_card_parser', '_card_parser_b', '_card_parser_prepare',
            'object_methods', 'object_attributes',
        ]
        return object_attributes(self, mode=mode, keys_to_skip=keys_to_skip+my_keys_to_skip,
                                 filter_properties=filter_properties)

    def object_methods(self, mode: str='public', keys_to_skip: Optional[list[str]]=None) -> list[str]:
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
        my_keys_to_skip = []  # type: list[str]

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

    def deprecated(self, old_name: str, new_name: str, deprecated_version: str) -> None:
        """deprecates methods"""
        return deprecated(old_name, new_name, deprecated_version, levels=[0, 1, 2])

    def clear_attributes(self) -> None:
        """removes the attributes from the model"""
        self.__init_attributes()

        self.nodes = {}
        self.loads = {}  # type: dict[int, list[Any]]
        self.load_combinations = {}  # type: dict[int, list[Any]]

    def reset_errors(self) -> None:
        """removes the errors from the model"""
        self._ixref_errors = 0
        self._stored_xref_errors = []

    def __init_attributes(self) -> None:
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
        self.system_command_lines = []  # type: list[str]

        #: list of execive control deck lines
        self.executive_control_lines = []  # type: list[str]

        #: list of case control deck lines
        self.case_control_lines = []  # type: list[str]

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
            103: 'SEMODES',   # modal
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
            #187 - Dynamic Design Analysis Method
            190: 'DBTRANS',
            200: 'DESOPT',  # optimization
        }

        # ------------------------ bad duplicates ----------------------------
        self._iparse_errors = 0
        self._nparse_errors = 0
        self._stop_on_parsing_error = True
        self._stop_on_duplicate_error = True
        self._stored_parse_errors = []  # type: list[str]

        self._duplicate_nodes = []  # type: list[str]
        self._duplicate_elements = []  # type: list[str]
        self._duplicate_properties = []  # type: list[str]
        self._duplicate_materials = []  # type: list[str]
        self._duplicate_masses = []  # type: list[str]
        self._duplicate_thermal_materials = []  # type: list[str]
        self._duplicate_coords = []  # type: list[str]
        self.values_to_skip = {}  # type: dict[str, list[int]]

        # ------------------------ structural defaults -----------------------
        #: the analysis type
        self._sol = None
        #: used in solution 600, method
        self.sol_method = None
        #: the line with SOL on it, marks ???
        self.sol_iline = None  # type : Optional[int]
        self.case_control_deck = None  # type: Optional[CaseControlDeck]

        #: store the PARAM cards
        self.params = {}    # type: dict[str, PARAM]
        self.mdlprm = None  # type: MDLPRM
        # ------------------------------- nodes -------------------------------
        # main structural block
        #: stores POINT cards
        self.points = {}  # type: dict[int, POINT]
        #self.grids = {}

        self.spoints = {}  # type: dict[int, SPOINT]
        self.epoints = {}  # type: dict[int, EPOINT]

        #: stores GRIDSET card
        self.grdset = None  # type: Optional[GRDSET]

        #: stores SEQGP cards
        self.seqgp = None  # type: Optional[SEQGP]

        ## stores RINGAX
        self.ringaxs = {}  # type: dict[int, RINGAX]

        ## stores GRIDB
        self.gridb = {}  # type: dict[int, GRIDB]

        #: stores elements (CQUAD4, CTRIA3, CHEXA8, CTETRA4, CROD, CONROD,
        #: etc.)
        self.elements = {}  # type: dict[int, Any]

        #: stores CBARAO, CBEAMAO
        self.ao_element_flags = {}  # type: dict[int, Any]
        #: stores BAROR
        self.baror = None  # type: Optional[BAROR]
        #: stores BEAMOR
        self.beamor = None  # type: Optional[BEAMOR]
        #: stores SNORM
        self.normals = {}  # type: dict[int, SNORM]

        #: stores rigid elements (RBE2, RBE3, RJOINT, etc.)
        self.rigid_elements = {}  # type: dict[int, Any]
        #: stores PLOTELs
        self.plotels = {}  # type: Optional[PLOTEL]

        #: stores CONM1, CONM2, CMASS1,CMASS2, CMASS3, CMASS4, CMASS5
        self.masses = {}  # type: dict[int, Any]
        #: stores PMASS
        self.properties_mass = {}  # type: dict[int, Any]

        #: stores NSM, NSM1, NSML, NSML1
        self.nsms = {}  # type: dict[int, list[Any]]
        #: stores NSMADD
        self.nsmadds = {}  # type: dict[int, list[Any]]

        #: stores LOTS of properties (PBAR, PBEAM, PSHELL, PCOMP, etc.)
        self.properties = {}  # type: dict[int, Any]

        #: stores MAT1, MAT2, MAT3, MAT8, MAT10, MAT11
        self.materials = {}  # type: dict[int, Any]

        #: defines the MAT4, MAT5
        self.thermal_materials = {}  # type: dict[int, Any]

        #: defines the MATHE, MATHP
        self.hyperelastic_materials = {}  # type: dict[int, Any]

        #: stores MATSx
        self.MATS1 = {}  # type: dict[int, Any]
        self.MATS3 = {}  # type: dict[int, Any]
        self.MATS8 = {}  # type: dict[int, Any]

        #: stores MATTx
        self.MATT1 = {}  # type: dict[int, Any]
        self.MATT2 = {}  # type: dict[int, Any]
        self.MATT3 = {}  # type: dict[int, Any]
        self.MATT4 = {}  # type: dict[int, Any]
        self.MATT5 = {}  # type: dict[int, Any]
        self.MATT8 = {}  # type: dict[int, Any]
        self.MATT9 = {}  # type: dict[int, Any]
        self.MATDMG = {}  # type: dict[int, Any]
        self.nxstrats = {}  # type: dict[int, Any]

        #: stores the CREEP card
        self.creep_materials = {}  # type: dict[int, Any]

        self.tics = {}  # type: Optional[Any]

        # stores DLOAD entries.
        self.dloads = {}    # type: dict[int, Any]
        # stores ACSRCE, RLOAD1, RLOAD2, TLOAD1, TLOAD2, and ACSRCE,
        #        and QVECT entries.
        self.dload_entries = {}    # type: dict[int, Any]

        #self.gusts = {}  # Case Control GUST = 100
        #self.random = {} # Case Control RANDOM = 100

        #: stores coordinate systems
        origin = array([0., 0., 0.])
        zaxis = array([0., 0., 1.])
        xzplane = array([1., 0., 0.])
        coord = CORD2R(cid=0, rid=0, origin=origin, zaxis=zaxis, xzplane=xzplane)
        self.coords = {0 : coord}   # type: dict[int, Any]
        self.MATCID = {}

        # --------------------------- constraints ----------------------------
        #: stores SUPORT1s
        #self.constraints = {} # suport1, anything else???
        self.suport = []  # type: list[Any]
        self.suport1 = {}  # type: dict[int, Any]
        self.se_suport = []  # type: list[Any]

        #: stores SPC, SPC1, SPCAX, GMSPC
        self.spcs = {}  # type: dict[int, list[Any]]
        #: stores SPCADD
        self.spcadds = {}  # type: dict[int, list[Any]]
        self.spcoffs = {}  # type: dict[int, list[Any]]

        self.mpcs = {}  # type: dict[int, list[Any]]
        self.mpcadds = {}  # type: dict[int, list[Any]]

        # --------------------------- dynamic ----------------------------
        #: stores DAREA
        self.dareas = {}   # type: dict[int, Any]
        self.dphases = {}  # type: dict[int, Any]

        self.pbusht = {}  # type: dict[int, Any]
        self.pdampt = {}  # type: dict[int, Any]
        self.pelast = {}  # type: dict[int, Any]

        #: frequencies
        self.frequencies = {}  # type: dict[int, list[Any]]

        # ----------------------------------------------------------------
        #: direct matrix input - DMIG
        self.dmi = {}  # type: dict[str, DMI]
        self.dmig = {}  # type: dict[str, DMIG]
        self.dmij = {}  # type: dict[str, DMIJ]
        self.dmiji = {}  # type: dict[str, DMIJI]
        self.dmik = {}  # type: dict[str, DMIK]
        self.dmiax = {}  # type: dict[str, DMIAX]
        self.dti = {}  # type: dict[str, Any]
        self._dmig_temp = defaultdict(list)  # type: dict[str, list[str]]

        # ----------------------------------------------------------------
        #: SETy
        self.sets = {}  # type: dict[int, Any]
        self.asets = []  # type: list[Any]
        self.omits = []  # type: list[Any]
        self.bsets = []  # type: list[Any]
        self.csets = []  # type: list[Any]
        self.qsets = []  # type: list[Any]
        self.usets = {}  # type: dict[str, Any]

        #: SExSETy
        self.se_bsets = []  # type: list[Any]
        self.se_csets = []  # type: list[Any]
        self.se_qsets = []  # type: list[Any]
        self.se_usets = {}  # type: dict[str, Any]
        self.se_sets = {}  # type: dict[str, Any]

        # ----------------------------------------------------------------
        #: parametric
        self.pset = {}
        self.pval = {}
        self.gmcurv = {}
        self.gmsurf = {}
        self.feedge = {}
        self.feface = {}

        # ----------------------------------------------------------------
        #: tables
        # TABLES1, ...
        self.tables = {}  # type: dict[int, TABLES1]

        # TABLEDx
        self.tables_d = {}  # type: dict[int, Union[TABLED1, TABLED2, TABLED3, TABLED4]]

        # TABLEMx
        self.tables_m = {}  # type: dict[int, Union[TABLEM1, TABLEM2, TABLEM3, TABLEM4]]

        #: random_tables
        self.random_tables = {}  # type: dict[int, Any]
        #: TABDMP1
        self.tables_sdamping = {}  # type: dict[int, TABDMP1]

        # ----------------------------------------------------------------
        #: EIGB, EIGR, EIGRL methods
        self.methods = {}  # type: dict[int, Union[EIGR, EIGRL, EIGB]]
        # EIGC, EIGP methods
        self.cMethods = {}  # type: dict[int, Union[EIGC, EIGP]]

        # ---------------------------- optimization --------------------------
        # optimization
        self.dconadds = {}  # type: dict[int, DCONADD]
        self.dconstrs = {}  # type: dict[int, DCONSTR]
        self.desvars = {}  # type: dict[int, DESVAR]
        self.topvar = {}  # type: dict[int, TOPVAR]
        self.ddvals = {}  # type: dict[int, DDVAL]
        self.dlinks = {}  # type: dict[int, DLINK]
        self.dresps = {}  # type: dict[int, Union[DRESP1, DRESP2, DRESP3]]

        self.dtable = None  # type: Optional[DTABLE]
        self.dequations = {}  # type: dict[int, DEQATN]

        #: stores DVPREL1, DVPREL2...might change to DVxRel
        self.dvprels = {}  # type: dict[int, Union[DVPREL1, DVPREL2]]
        self.dvmrels = {}  # type: dict[int, Union[DVMREL1, DVMREL2]]
        self.dvcrels = {}  # type: dict[int, Union[DVCREL1, DVCREL2]]
        self.dvgrids = {}  # type: dict[int, DVGRID]
        self.doptprm = None  # type: Optional[DOPTPRM]
        self.dscreen = {}  # type: dict[int, DSCREEN]

        # nx optimization
        self.group = {}    # type: dict[int, GROUP]
        self.dmncon = {}   # type: dict[int, DMNCON]
        self.dvtrels = {}  # type: dict[int, Union[DVTREL1, DVTREL2]]

        # ------------------------- nonlinear defaults -----------------------
        #: stores NLPCI
        self.nlpcis = {}  # type: dict[int, NLPCI]
        #: stores NLPARM
        self.nlparms = {}  # type: dict[int, NLPARM]
        #: stores TSTEPs, TSTEP1s
        self.tsteps = {}  # type: dict[int, Union[TSTEP, TSTEP1]]
        #: stores TSTEPNL
        self.tstepnls = {}  # type: dict[int, TSTEPNL]
        #: stores TF
        self.transfer_functions = {}  # type: dict[int, TF]
        #: stores DELAY
        self.delays = {}  # type: dict[int, DELAY]

        #: stores ROTORD, ROTORG
        self.rotors = {}  # type: dict[int, Union[ROTORD, ROTORG]]

        # --------------------------- aero defaults --------------------------
        # aero cards
        #: stores CAEROx
        self.caeros = {}  # type: dict[int, Union[CAERO1, CAERO2, CAERO3, CAERO4, CAERO5]]
        #: stores PAEROx
        self.paeros = {}  # type: dict[int, Union[PAERO1, PAERO2, PAERO3, PAERO4, PAERO5]]
        # stores MONPNT1
        self.monitor_points = []  # type: list[Union[MONPNT1, MONPNT2, MONPNT3]]

        #: stores AECOMP
        self.aecomps = {}  # type: dict[int, AECOMP]
        #: stores AEFACT
        self.aefacts = {}  # type: dict[int, AEFACT]
        #: stores AELINK
        self.aelinks = {}  # type: dict[int, list[AELINK]]
        #: stores AELIST
        self.aelists = {}  # type: dict[int, AELIST]
        #: stores AEPARM
        self.aeparams = {}  # type: dict[int, AEPARM]
        #: stores AESURF
        self.aesurf = {}  # type: dict[int, AESURF]
        #: stores AESURFS
        self.aesurfs = {}  # type: dict[int, AESURFS]
        #: stores AESTAT
        self.aestats = {}  # type: dict[int, AESTAT]
        #: stores CSSCHD
        self.csschds = {}  # type: dict[int, CSSCHD]

        #: store SPLINE1,SPLINE2,SPLINE4,SPLINE5
        self.splines = {}  # type: dict[int, Union[SPLINE1, SPLINE2, SPLINE3, SPLINE4, SPLINE5]]
        self.zona = ZONA(self)

        # axisymmetric
        self.axic = None  # type: Optional[AXIC]
        self.axif = None  # type: Optional[AXIF]
        self.ringfl = {}  # type: dict[int, RINGFL]
        self._is_axis_symmetric = False

        # cyclic
        self.cyax = None  # type: Optional[CYAX]
        self.cyjoin = {}  # type: dict[int, CYJOIN]

        self.modtrak = None  # type: Optional[MODTRAK]

        # acoustic
        self.acmodl = None

        # ------ SOL 144 ------
        #: stores AEROS
        self.aeros = None  # type: Optional[AEROS]

        #: stores TRIM, TRIM2
        self.trims = {}  # type: dict[int, Union[TRIM, TRIM2]]

        #: stores DIVERG
        self.divergs = {}  # type: dict[int, DIVERG]

        # ------ SOL 145 ------
        #: stores AERO
        self.aero = None  # type: Optional[AERO]

        #: stores FLFACT
        self.flfacts = {}  # type: dict[int, FLFACT]

        #: stores FLUTTER
        self.flutters = {} # type: dict[int, FLUTTER]

        #: mkaeros
        self.mkaeros = []  # type: list[Union[MKAERO1,MKAERO2]]

        # ------ SOL 146 ------
        #: stores GUST cards
        self.gusts = {}  # type: dict[int, GUST]

        # ------------------------- thermal defaults -------------------------
        # BCs
        #: stores thermal boundary conditions - CONV,RADBC
        self.bcs = {}  # type: dict[int, Union[CONV, RADBC]]

        #: stores PHBDY
        self.phbdys = {}  # type: dict[int, PHBDY]
        #: stores convection properties - PCONV, PCONVM ???
        self.convection_properties = {}  # type: dict[int, Union[PCONV, PCONVM]]
        #: stores TEMPD
        self.tempds = {}  # type: dict[int, TEMPD]

        #: stores VIEW
        self.views = {}  # type: dict[int, VIEW]
        #: stores VIEW3D
        self.view3ds = {}  # type: dict[int, VIEW3D]
        self.radset = None
        self.radcavs = {}  # type: dict[int, RADCAV]
        self.radmtx = {}  # type: dict[int, RADMTX]

        # -------------------------contact cards-------------------------------
        self.bcbodys = {}  # type: dict[int, BCBODY]
        self.bcparas = {}  # type: dict[int, BCPARA]
        self.bcrparas = {}  # type: dict[int, BCRPARA]
        self.bctparas = {}  # type: dict[int, BCTPARA]

        self.bctadds = {}  # type: dict[int, BCTADD]
        self.bctsets = {}  # type: dict[int, BCTSET]
        self.bsurf = {}  # type: dict[int, BSURF]
        self.bsurfs = {}  # type: dict[int, BSURFS]
        self.bconp = {}  # type: dict[int, BCONP]
        self.blseg = {}  # type: dict[int, BLSEG]
        self.bfric = {}  # type: dict[int, BFRIC]
        self.bgadds = {}  # type: dict[int, BGADD]
        self.bgsets = {}  # type: dict[int, BGSET]
        self.bctparms = {}  # type: dict[int, BCTPARAM]

        #--------------------------superelements------------------------------
        self.setree = {}  # type: dict[int, SETREE]
        self.senqset = {}  # type: dict[int, Union[SENQSET, SENQSET1]]
        self.sebulk = {}  # type: dict[int, SEBULK]
        self.sebndry = {}  # type: dict[int, SEBNDRY]
        self.release = {}  # type: dict[int, RELEASE]
        self.seloc = {}  # type: dict[int, SELOC]
        self.sempln = {}  # type: dict[int, SEMPLN]
        self.seconct = {}  # type: dict[int, SECONCT]
        self.selabel = {}  # type: dict[int, SELABEL]
        self.seexcld = {}  # type: dict[int, SEEXCLD]
        self.seelt = {}  # type: dict[int, SEELT]
        self.seload = {}  # type: dict[int, SELOAD]
        self.csuper = {}  # type: dict[int, CSUPER]
        self.csupext = {}  # type: dict[int, CSUPEXT]

        self.bolt = {} # type: dict[int, BOLT]
        self.boltseq = {} # type: dict[int, BOLTSEQ]
        self.boltfor = {} # type: dict[int, BOLTFOR]
        self.boltfrc = {}
        self.boltld = {}
        # ---------------------------------------------------------------------
        self._type_to_id_map = defaultdict(list)  # type: dict[int, list[Any]]
        self._slot_to_type_map = {
            'params' : ['PARAM'],
            'mdlprm': ['MDLPRM'],
            'nodes' : ['GRID', 'SPOINT', 'EPOINT'], # 'RINGAX',
            'points' : ['POINT'],
            'ringaxs' : ['RINGAX', 'POINTAX'],
            'ringfl' : ['RINGFL'],
            'axic' : ['AXIC'],
            'axif' : ['AXIF'],
            'acmodl' : ['ACMODL'],
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
                'CCONEAX',

                'CTETRA', 'CPYRAM', 'CPENTA', 'CHEXA', 'CIHEX1', 'CIHEX2', 'CHEXA1', 'CHEXA2',
                'CSHEAR', 'CVISC', 'CRAC2D', 'CRAC3D',
                'CGAP',

                # thermal
                'CHBDYE', 'CHBDYG', 'CHBDYP',

                # acoustic
                'CHACAB', 'CAABSF', 'CHACBR',
            ],
            'normals' : ['SNORM'],
            'nsms' : ['NSM', 'NSM1', 'NSML', 'NSML1'],
            'nsmadds' : ['NSMADD'],
            'rigid_elements' : ['RBAR', 'RBAR1', 'RBE1', 'RBE2', 'RBE3', 'RROD', 'RSPLINE', 'RSSCON'],
            'plotels' : ['PLOTEL'],

            'properties_mass' : ['PMASS'],
            #'properties_acoustic' : ['PACABS'],
            'properties' : [
                #  acoustic
                'PACABS', 'PAABSF', 'PACBAR',

                # 0d
                'PELAS', 'PGAP', 'PFAST',
                'PBUSH', 'PBUSH1D',
                'PDAMP', 'PDAMP5',

                # 1d
                'PROD', 'PBAR', 'PBARL', 'PBEAM', 'PTUBE', 'PBEND', 'PBCOMP', 'PBRSECT', 'PBMSECT',
                'PBEAML',  # not fully supported
                'PBEAM3',

                # 2d
                'PLPLANE', 'PPLANE',
                'PSHELL', 'PCOMP', 'PCOMPG', 'PSHEAR',
                'PSOLID', 'PLSOLID', 'PVISC', 'PRAC2D', 'PRAC3D',
                'PIHEX', 'PCOMPS', 'PCOMPLS',
                'PCONEAX',
            ],
            'pdampt' : ['PDAMPT'],
            'pelast' : ['PELAST'],
            'pbusht' : ['PBUSHT'],

            # materials
            'materials' : ['MAT1', 'MAT2', 'MAT3', 'MAT8', 'MAT9', 'MAT10', 'MAT11',
                           'MAT3D', 'MATG'],
            'hyperelastic_materials' : ['MATHE', 'MATHP'],
            'creep_materials' : ['CREEP'],
            'MATT1' : ['MATT1'],
            'MATT2' : ['MATT2'],
            'MATT3' : ['MATT3'],
            'MATT4' : ['MATT4'], # thermal
            'MATT5' : ['MATT5'], # thermal
            'MATT8' : ['MATT8'],
            'MATT9' : ['MATT9'],
            'MATS1' : ['MATS1'],
            'MATDMG': ['MATDMG'],
            'MATS3' : ['MATS3'],
            'MATS8' : ['MATS8'],
            'nxstrats' : ['NXSTRAT'],

            # 'MATHE'
            #'EQUIV', # testing only, should never be activated...

            # thermal materials
            'thermal_materials' : ['MAT4', 'MAT5'],

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
            'release' : ['RELEASE'],
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
            'load_combinations' : ['LOAD', 'LSEQ', 'CLOAD'],
            'loads' : [
                'FORCE', 'FORCE1', 'FORCE2',
                'MOMENT', 'MOMENT1', 'MOMENT2',
                'GRAV', 'ACCEL', 'ACCEL1',
                'PLOAD', 'PLOAD1', 'PLOAD2', 'PLOAD4',
                'RFORCE', 'RFORCE1', 'SLOAD',
                'SPCD', 'LOADCYN', 'LOADCYH', 'DEFORM',

                # msgmesh
                #'GMLOAD',

                # thermal
                'TEMP', 'TEMPB3', 'TEMPRB',
                'QBDY1', 'QBDY2', 'QBDY3', 'QHBDY', 'QVOL',

                # axisymmetric
                'PLOADX1', 'FORCEAX', 'PRESAX', 'TEMPAX',
                ],
            'cyjoin' : ['CYJOIN'],
            'cyax' : ['CYAX'],
            'modtrak' : ['MODTRAK'],
            'dloads' : ['DLOAD'],
            # stores RLOAD1, RLOAD2, TLOAD1, TLOAD2, and ACSRCE entries.
            'dload_entries' : ['ACSRCE', 'TLOAD1', 'TLOAD2', 'RLOAD1', 'RLOAD2',
                               'QVECT', 'RANDPS', 'RANDT1'],

            # aero cards
            'aero' : ['AERO'],
            'aeros' : ['AEROS'],
            'gusts' : ['GUST', 'GUST2'],
            'flutters' : ['FLUTTER'],
            'flfacts' : ['FLFACT'],
            'mkaeros' : ['MKAERO1', 'MKAERO2'],
            'aecomps' : ['AECOMP', 'AECOMPL'],
            'aefacts' : ['AEFACT'],
            'aelinks' : ['AELINK'],
            'aelists' : ['AELIST'],
            'aeparams' : ['AEPARM'],
            'aesurf' : ['AESURF'],
            'aesurfs' : ['AESURFS'],
            'aestats' : ['AESTAT'],
            'caeros' : ['CAERO1', 'CAERO2', 'CAERO3', 'CAERO4', 'CAERO5', 'CAERO7', 'BODY7'],
            'paeros' : ['PAERO1', 'PAERO2', 'PAERO3', 'PAERO4', 'PAERO5', 'SEGMESH'],
            'monitor_points' : ['MONPNT1', 'MONPNT2', 'MONPNT3', 'MONDSP1'],
            'splines' : ['SPLINE1', 'SPLINE2', 'SPLINE3', 'SPLINE4', 'SPLINE5', 'SPLINE6', 'SPLINE7'],
            'panlsts' : ['PANLST1', 'PANLST2', 'PANLST3'],
            'csschds' : ['CSSCHD',],
            #'SPLINE3', 'SPLINE6', 'SPLINE7',
            'trims' : ['TRIM', 'TRIM2'],
            'divergs' : ['DIVERG'],

            # coords
            'coords' : ['CORD1R', 'CORD1C', 'CORD1S',
                        'CORD2R', 'CORD2C', 'CORD2S',
                        'GMCORD', 'ACOORD', 'CORD3G'],
            'matcid': ['MATCID'],

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
            'dmig' : ['DMIG'],
            'dmiax' : ['DMIAX'],
            'dmij' : ['DMIJ'],
            'dmiji' : ['DMIJI'],
            'dmik' : ['DMIK'],
            'dmi' : ['DMI'],
            'dti' : ['DTI'],

            # optimzation
            'dequations' : ['DEQATN'],
            'dtable' : ['DTABLE'],
            'dconstrs' : ['DCONSTR', 'DCONADD'],
            'desvars' : ['DESVAR'],
            'topvar' : ['TOPVAR'],
            'ddvals' : ['DDVAL'],
            'dlinks' : ['DLINK'],
            'dresps' : ['DRESP1', 'DRESP2', 'DRESP3'],
            'dvprels' : ['DVPREL1', 'DVPREL2'],
            'dvmrels' : ['DVMREL1', 'DVMREL2'],
            'dvcrels' : ['DVCREL1', 'DVCREL2'],
            'dvgrids' : ['DVGRID'],
            'doptprm' : ['DOPTPRM'],
            'dscreen' : ['DSCREEN'],

            # optimization - nx
            'dmncon' : ['DMNCON'],
            'dvtrels' : ['DVTREL1'],
            'group' : ['GROUP'],

            # sets
            'asets' : ['ASET', 'ASET1'],
            'omits' : ['OMIT', 'OMIT1'],
            'bsets' : ['BSET', 'BSET1'],
            'qsets' : ['QSET', 'QSET1'],
            'csets' : ['CSET', 'CSET1'],
            'usets' : ['USET', 'USET1'],
            'sets' : ['SET1', 'SET3'],

            # super-element sets
            'se_bsets' : ['SEBSET', 'SEBSET1'],
            'se_csets' : ['SECSET', 'SECSET1'],
            'se_qsets' : ['SEQSET', 'SEQSET1'],
            'se_usets' : ['SEUSET', 'SEQSET1'],
            'se_sets' : ['SESET'],
            'radset' : ['RADSET'],
            'radcavs' : ['RADCAV', 'RADLST'],
            'radmtx' : ['RADMTX'],
            # SEBSEP

            # parametric
            'pset' : ['PSET'],
            'pval' : ['PVAL'],
            'gmcurv' : ['GMCURV'],
            'gmsurf' : ['GMSURF'],
            'feedge' : ['FEEDGE'],
            'feface' : ['FEFACE'],

            # tables
            'tables' : [
                'TABLEH1', 'TABLEHT',
                'TABLES1', 'TABLEST',
                ],
            'tables_d' : ['TABLED1', 'TABLED2', 'TABLED3', 'TABLED4', 'TABLED5'],
            'tables_m' : ['TABLEM1', 'TABLEM2', 'TABLEM3', 'TABLEM4'],
            'tables_sdamping' : ['TABDMP1'],
            'random_tables' : ['TABRND1', 'TABRNDG'],

            # initial conditions - sid (set ID)
            ##'TIC',  (in bdf_tables.py)

            # methods
            'methods' : ['EIGB', 'EIGR', 'EIGRL'],

            # cMethods
            'cMethods' : ['EIGC', 'EIGP'],

            # contact
            'bcbodys' : ['BCBODY'],
            'bcparas' : ['BCPARA'],
            'bctparas' : ['BCTPARA'],
            'bcrparas' : ['BCRPARA'],
            'bctparms' : ['BCTPARM'],

            'bctadds' : ['BCTADD'],
            'bctsets' : ['BCTSET'],
            'bgadds' : ['BGADD'],
            'bgsets' : ['BGSET'],
            'bsurf' : ['BSURF'],
            'bsurfs' : ['BSURFS'],
            'bconp' : ['BCONP'],
            'blseg' : ['BLSEG'],
            'bfric' : ['BFRIC'],
            'views' : ['VIEW'],
            'view3ds' : ['VIEW3D'],

            # nx bolts
            'bolt': ['BOLT'],
            'boltld': ['BOLTLD'],
            'boltfor': ['BOLTFOR'],
            'boltfrc': ['BOLTFRC'],
            'boltseq': ['BOLTSEQ'],
            ## other
            #'INCLUDE',  # '='
            #'ENDDATA',
        }  # type: dict[str, list[str]]
        self._type_to_slot_map = self.get_rslot_map()

    @property
    def type_slot_str(self) -> str:
        """helper method for printing supported cards"""
        nchars = len('Card Group')

        #nchars_cards = 0
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

            #i = 0
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
            for isub, subline in enumerate(sublines):
                if isub > 0:  # adds intermediate dash lines
                    html_msg.append(dash_plus)
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
    def nastran_format(self) -> str:
        return self._nastran_format

    @nastran_format.setter
    def nastran_format(self, nastran_format: str) -> None:
        assert isinstance(nastran_format, str), nastran_format
        fmt_lower = nastran_format.lower().strip()
        if fmt_lower not in BDF_FORMATS:
            raise RuntimeError(nastran_format)
        self._nastran_format = fmt_lower

    @property
    def is_long_ids(self) -> bool:
        return self._is_long_ids
        #if self._nastran_format == 'nx' or self._is_long_ids:
            #return True
        #return False

    def _set_punch(self) -> None:
        """updates the punch flag"""
        if self.punch is None:
            # writing a mesh without using read_bdf
            if self.system_command_lines or self.executive_control_lines or self.case_control_deck:
                self.punch = False
            else:
                self.punch = True

    @property
    def sol(self) -> int:
        """gets the solution (e.g. 101, 103)"""
        return self._sol

    @sol.setter
    def sol(self, sol: int) -> int:
        """sets the solution (e.g. 101, 103)"""
        self._sol = sol
        if len(self.executive_control_lines) == 0:
            self.executive_control_lines = ['SOL %s' % sol, 'CEND']
            self.sol_iline = 0
        return self._sol

    @property
    def subcases(self) -> dict[int, Subcase]:
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
    def nnodes(self) -> int:
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
    def npoints(self) -> int:
        """gets the number of GRID, SPOINT, EPOINT ids"""
        return len(self.point_ids)

    #--------------------
    # Elements CARDS

    @property
    def nelements(self) -> int:
        """gets the number of element"""
        return len(self.elements)

    @property
    def element_ids(self):
        """gets the element ids"""
        return self.elements.keys()

    #--------------------
    # Property CARDS

    @property
    def nproperties(self) -> int:
        """gets the number of properties"""
        return len(self.properties)

    @property
    def property_ids(self):
        """gets the property ids"""
        return self.properties.keys()

    #--------------------
    # Material CARDS

    @property
    def nmaterials(self) -> int:
        """gets the number of materials"""
        return len(self.materials)

    @property
    def material_ids(self):
        """gets the material ids"""
        return self.materials.keys()

    #--------------------
    # Coords CARDS

    @property
    def ncoords(self) -> int:
        """gets the number of coordinate systems"""
        return len(self.coords)

    @property
    def coord_ids(self):
        """gets the number of coordinate system ids"""
        return self.coords.keys()

    #--------------------

    @property
    def ncaeros(self) -> int:
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

    def set_param(self, key: str, values: Union[int, float, str, list[float]], comment: str='') -> None:
        """sets a param card; creates it if necessary"""
        if isinstance(values, (int, float, str)):
            values = [values]
        key = key.upper()
        if key in self.params:
            param = self.params[key]
            param.update_values(*values)
        else:
            self.add_param(key, values, comment=comment)

    def get_param(self, key: str, default: Union[int, float, str, list[float]]
                  ) -> Union[int, float, str, list[float]]:
        """gets a param card"""
        key = key.upper()
        if key in self.params:
            param = self.params[key]
            return param.value
        return default

    #--------------------
    # deprecations
    @property
    def dmis(self) -> dict[str, DMI]:
        return self.dmi
    @property
    def dmigs(self) -> dict[str, DMIG]:
        return self.dmig
    @property
    def dmiks(self) -> dict[str, DMIK]:
        return self.dmik
    @property
    def dmijs(self) -> dict[str, DMIJ]:
        return self.dmij
    @property
    def dmijis(self) -> dict[str, DMIJI]:
        return self.dmiji

    @dmis.setter
    def dmis(self, dmi):
        self.dmi = dmi
    @dmigs.setter
    def dmigs(self, dmig):
        self.dmig = dmig
    @dmiks.setter
    def dmiks(self, dmik):
        self.dmik = dmik
    @dmijs.setter
    def dmijs(self, dmij):
        self.dmij = dmij
    @dmijis.setter
    def dmijis(self, dmiji):
        self.dmiji = dmiji
