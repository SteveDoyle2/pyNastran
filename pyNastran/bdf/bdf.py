# coding: utf-8
# pylint: disable=W0201,R0915,R0912
"""
Main BDF class.  Defines:
  - BDF

"""
# see https://docs.plm.automation.siemens.com/tdoc/nxnastran/10/help/#uid:index
import os
import sys
from copy import deepcopy
from io import StringIO, IOBase
import traceback
from collections import defaultdict

from typing import List, Dict, Set, Tuple, Sequence, Optional, Union, Any # , cast
from pickle import load, dump, dumps  # type: ignore

import numpy as np  # type: ignore
from cpylog import get_logger2

from pyNastran.utils import object_attributes, check_path
from .utils import parse_patran_syntax
from .bdf_interface.utils import (
    _parse_pynastran_header, to_fields, parse_executive_control_deck,
    fill_dmigs, _get_card_name, _parse_dynamic_syntax,
)

from .bdf_interface.replication import (
    to_fields_replication, get_nrepeats, int_replication, float_replication,
    _field, repeat_cards)

from .field_writer_8 import print_card_8
from .field_writer_16 import print_card_16, print_field_16

from .cards.base_card import _format_comment
from .cards.utils import wipe_empty_fields

#from .write_path import write_include
from .bdf_interface.assign_type import (integer,
                                        integer_or_string, string)

from .cards.elements.elements import CFAST, CGAP, CRAC2D, CRAC3D, PLOTEL, GENEL
from .cards.properties.properties import PFAST, PGAP, PRAC2D, PRAC3D
from .cards.properties.solid import PLSOLID, PSOLID, PIHEX, PCOMPS
from .cards.cyclic import CYAX, CYJOIN
from .cards.msgmesh import CGEN

from .cards.elements.springs import CELAS1, CELAS2, CELAS3, CELAS4
from .cards.properties.springs import PELAS, PELAST

from .cards.elements.solid import (
    #CTETRA, CPYRAM, CPENTA, CHEXA,
    CIHEX1, CIHEX2,
    CTETRA4, CPYRAM5, CPENTA6, CHEXA8,
    CTETRA10, CPYRAM13, CPENTA15, CHEXA20,
)
from .cards.elements.rigid import RBAR, RBAR1, RBE1, RBE2, RBE3, RROD, RSPLINE, RSSCON

from .cards.axisymmetric.axisymmetric import (
    AXIF, RINGFL,
    AXIC, RINGAX, POINTAX, CCONEAX, PCONEAX, )
from .cards.axisymmetric.loads import PLOADX1, FORCEAX, PRESAX, TEMPAX
from .cards.elements.axisymmetric_shells import (
    CTRAX3, CTRAX6, CTRIAX, CTRIAX6, CQUADX, CQUADX4, CQUADX8)
from pyNastran.bdf.cards.elements.shell import (
    CQUAD, CQUAD4, CQUAD8, CQUADR, CSHEAR,
    CTRIA3, CTRIA6, CTRIAR,
    CPLSTN3, CPLSTN4, CPLSTN6, CPLSTN8,
    CPLSTS3, CPLSTS4, CPLSTS6, CPLSTS8,
    SNORM)
from .cards.properties.shell import PSHELL, PCOMP, PCOMPG, PSHEAR, PLPLANE, PPLANE
from .cards.elements.acoustic import (
    CHACAB, CAABSF, CHACBR, PACABS, PAABSF, PACBAR, ACMODL)
from .cards.elements.bush import CBUSH, CBUSH1D, CBUSH2D
from .cards.properties.bush import PBUSH, PBUSH1D, PBUSHT
from .cards.elements.damper import (CVISC, CDAMP1, CDAMP2, CDAMP3, CDAMP4,
                                    CDAMP5)
from .cards.properties.damper import PVISC, PDAMP, PDAMP5, PDAMPT
from .cards.elements.rods import CROD, CONROD, CTUBE
from .cards.elements.bars import CBAR, BAROR, CBARAO, CBEAM3, CBEND
from .cards.elements.beam import CBEAM, BEAMOR
from .cards.properties.rods import PROD, PTUBE
from .cards.properties.bars import PBAR, PBARL, PBRSECT, PBEND, PBEAM3
from .cards.properties.beam import PBEAM, PBEAML, PBCOMP, PBMSECT
# CMASS5
from .cards.elements.mass import CONM1, CONM2, CMASS1, CMASS2, CMASS3, CMASS4
from .cards.properties.mass import PMASS, NSM, NSM1, NSML, NSML1, NSMADD
from .cards.constraints import (SPC, SPCADD, SPCAX, SPC1, SPCOFF, SPCOFF1,
                                MPC, MPCADD, SUPORT1, SUPORT, SESUP,
                                GMSPC)
from .cards.coordinate_systems import (CORD1R, CORD1C, CORD1S,
                                       CORD2R, CORD2C, CORD2S, #CORD3G,
                                       GMCORD, transform_coords_vectorized,
                                       CORDx)
from .cards.deqatn import DEQATN
from .cards.dynamic import (
    DELAY, DPHASE, FREQ, FREQ1, FREQ2, FREQ3, FREQ4, FREQ5,
    TSTEP, TSTEP1, TSTEPNL, NLPARM, NLPCI, TF, ROTORG, ROTORD, TIC)
from .cards.loads.loads import (
    LSEQ, SLOAD, DAREA, RFORCE, RFORCE1, SPCD, DEFORM, LOADCYN, LOADCYH)
from .cards.loads.dloads import ACSRCE, DLOAD, TLOAD1, TLOAD2, RLOAD1, RLOAD2
from .cards.loads.static_loads import (LOAD, CLOAD, GRAV, ACCEL, ACCEL1, FORCE,
                                       FORCE1, FORCE2, MOMENT, MOMENT1, MOMENT2,
                                       PLOAD, PLOAD1, PLOAD2, PLOAD4,
                                       GMLOAD)
from .cards.loads.random_loads import RANDPS, RANDT1

from .cards.materials import (MAT1, MAT2, MAT3, MAT4, MAT5,
                              MAT8, MAT9, MAT10, MAT11, MAT3D,
                              MATG, MATHE, MATHP, CREEP, EQUIV,
                              NXSTRAT)
from .cards.material_deps import (
    MATT1, MATT2, MATT3, MATT4, MATT5, MATT8, MATT9, MATS1)

from .cards.methods import EIGB, EIGC, EIGR, EIGP, EIGRL, MODTRAK
from .cards.nodes import GRID, GRDSET, SPOINTs, EPOINTs, POINT, SEQGP, GRIDB
from .cards.aero.aero import (
    AECOMP, AECOMPL, AEFACT, AELINK, AELIST, AEPARM, AESURF, AESURFS,
    CAERO1, CAERO2, CAERO3, CAERO4, CAERO5,
    PAERO1, PAERO2, PAERO3, PAERO4, PAERO5,
    MONPNT1, MONPNT2, MONPNT3, MONDSP1,
    SPLINE1, SPLINE2, SPLINE3, SPLINE4, SPLINE5)
from .cards.aero.static_loads import AESTAT, AEROS, CSSCHD, TRIM, TRIM2, DIVERG
from .cards.aero.dynamic_loads import AERO, FLFACT, FLUTTER, GUST, MKAERO1, MKAERO2
from .cards.optimization import (
    DCONADD, DCONSTR, DESVAR, TOPVAR, DDVAL, DOPTPRM, DLINK,
    DRESP1, DRESP2, DRESP3,
    DVCREL1, DVCREL2,
    DVMREL1, DVMREL2,
    DVPREL1, DVPREL2,
    DVGRID, DSCREEN)
from .cards.superelements import (
    RELEASE, SEBNDRY, SEBULK, SECONCT, SEELT, SEEXCLD,
    SELABEL, SELOAD, SELOC, SEMPLN, SENQSET, SETREE,
    CSUPER, CSUPEXT,
)
from .cards.bdf_sets import (
    ASET, BSET, CSET, QSET, USET,
    ASET1, BSET1, CSET1, QSET1, USET1,
    OMIT, OMIT1,
    SET1, SET3,
    SEBSET, SECSET, SEQSET, # SEUSET
    SEBSET1, SECSET1, SEQSET1, # SEUSET1
    SESET, #SEQSEP
    RADSET,
)
from .cards.params import PARAM
from .cards.dmig import DMIG, DMI, DMIJ, DMIK, DMIJI, DMIG_UACCEL, DTI, DMIAX
from .cards.thermal.loads import (QBDY1, QBDY2, QBDY3, QHBDY, TEMP, TEMPD, TEMPB3,
                                  TEMPRB, QVOL, QVECT)
from .cards.thermal.thermal import (CHBDYE, CHBDYG, CHBDYP, PCONV, PCONVM,
                                    PHBDY, CONV, CONVM, TEMPBC)
from .cards.thermal.radiation import RADM, RADBC, RADCAV, RADLST, RADMTX, VIEW, VIEW3D
from .cards.bdf_tables import (TABLED1, TABLED2, TABLED3, TABLED4,
                               TABLEM1, TABLEM2, TABLEM3, TABLEM4,
                               TABLES1, TABDMP1, TABLEST, TABLEHT, TABLEH1,
                               TABRND1, TABRNDG,
                               DTABLE)
from .cards.contact import (
    BCRPARA, BCTADD, BCTSET, BSURF, BSURFS, BCTPARA, BCONP, BLSEG, BFRIC)
from .cards.parametric.geometry import PSET, PVAL, FEEDGE, FEFACE, GMCURV, GMSURF

from .case_control_deck import CaseControlDeck, Subcase
from .bdf_methods import BDFMethods
from .bdf_interface.get_card import GetCard
from .bdf_interface.add_card import AddCards
from .bdf_interface.bdf_card import BDFCard
from .bdf_interface.write_mesh_file import WriteMeshs
from .bdf_interface.uncross_reference import UnXrefMesh
from .bdf_interface.verify_validate import verify_bdf, validate_bdf
from .bdf_interface.stats import get_bdf_stats

from .errors import (CrossReferenceError, DuplicateIDsError,
                                  CardParseSyntaxError, UnsupportedCard, DisabledCardError,
                                  SuperelementFlagError, ReplicationError)
from .bdf_interface.pybdf import (
    BDFInputPy, _clean_comment, _clean_comment_bulk, EXECUTIVE_CASE_SPACES)

#from .bdf_interface.add_card import CARD_MAP


def load_bdf_object(obj_filename:str, xref: bool=True, log=None, debug: bool=True):
    model = BDF(log=log, debug=debug)
    model.load(obj_filename=obj_filename)
    model.cross_reference(xref=xref, xref_nodes=True, xref_elements=True,
                          xref_nodes_with_elements=True,
                          xref_properties=True,
                          xref_masses=True,
                          xref_materials=True,
                          xref_loads=True,
                          xref_constraints=True,
                          xref_aero=True,
                          xref_sets=True,
                          xref_optimization=True)
    return model


class BDF_(BDFMethods, GetCard, AddCards, WriteMeshs, UnXrefMesh):
    """
    Base class for the BDF Reader/Writer/Editor class that's used by the
    main BDF object and (temporarily) the in-development vectorized object
    to keep things working.

    If you add very few methods and attributes to this, you get the ``BDF``
    class.  The point of this class is to break out a attributes, so the
    names (e.g., nodes) can be reused when vectorize the data.

    """
    #: required for sphinx bug
    #: http://stackoverflow.com/questions/11208997/autoclass-and-instance-attributes
    #__slots__ = ['_is_dynamic_syntax']
    def __init__(self, debug: bool=True, log=None, mode: str='msc') -> None:
        # SimpleLogger
        """
        Initializes the BDF_ object

        Parameters
        ----------
        debug : bool/None; default=True
            used to set the logger if no logger is passed in
                True:  logs debug/info/error messages
                False: logs info/error messages
                None:  logs error messages
        log : logging module object / None
            if log is set, debug is ignored and uses the
            settings the logging object has
        mode : str; default='msc'
            the type of Nastran
            valid_modes = {'msc', 'nx'}

        """
        assert debug in [True, False, None], 'debug=%r' % debug
        self.echo = False
        self.read_includes = True

        # file management parameters
        self.active_filenames = []  # type: List[str]
        self.active_filename = None  # type: Optional[str]
        self.include_dir = ''
        self.dumplines = False

        self.log = get_logger2(log, debug)

        # list of all read in cards - useful in determining if entire BDF
        # was read & really useful in debugging
        self.card_count = {}  # type: Dict[str, int]

        # stores the card_count of cards that have been rejected
        self.reject_count = {}  # type: Dict[str, int]

        # allows the BDF variables to be scoped properly (i think...)
        GetCard.__init__(self)
        AddCards.__init__(self)
        BDFMethods.__init__(self)
        WriteMeshs.__init__(self)
        UnXrefMesh.__init__(self)

        # useful in debugging errors in input
        self.debug = debug

        # flag that allows for OpenMDAO-style optimization syntax to be used
        self._is_dynamic_syntax = False

        # lines that were rejected b/c they were for a card that isnt supported
        self.reject_lines = []  # type: List[List[str]]

        # cards that were created, but not processed
        self.reject_cards = []  # type: List[str]

        self.include_filenames = defaultdict(list)
        # self.__init_attributes()

        cards_to_read = [
            '/',
            'ECHOON', 'ECHOOFF',
            'PARAM',

            ## nodes
            'GRID', 'GRDSET', 'SPOINT', 'EPOINT', 'SEQGP', 'GRIDB',

            # points
            'POINT',
            #'GRIDG'

            ## ringfl
            'RINGFL',
            ## ringaxs
            'RINGAX', 'POINTAX',

            ## masses
            'CONM1', 'CONM2',
            'CMASS1', 'CMASS2', 'CMASS3', 'CMASS4',

            ## nsms
            'NSM', 'NSM1', 'NSML', 'NSML1',

            ## nsmadds
            'NSMADD',

            ## elements
            # springs
            'CELAS1', 'CELAS2', 'CELAS3', 'CELAS4', # 'CELAS5',
            # bushings
            'CBUSH', 'CBUSH1D', 'CBUSH2D',
            # dampers
            'CDAMP1', 'CDAMP2', 'CDAMP3', 'CDAMP4', 'CDAMP5',
            'CFAST',

            'CBAR', 'CBARAO', 'BAROR',
            'CROD', 'CTUBE', 'CBEAM', 'CBEAM3', 'CONROD', 'CBEND', 'BEAMOR',
            'CTRIA3', 'CTRIA6', 'CTRIAR',
            'CQUAD4', 'CQUAD8', 'CQUADR', 'CQUAD',
            'CTRAX3', 'CTRAX6', 'CTRIAX', 'CTRIAX6', 'CQUADX', 'CQUADX4', 'CQUADX8',
            'SNORM',

            'CPLSTN3', 'CPLSTN4', 'CPLSTN6', 'CPLSTN8', # plate strain
            'CPLSTS3', 'CPLSTS4', 'CPLSTS6', 'CPLSTS8', # plate stress

            # acoustic
            'CHACAB', 'CAABSF', 'CHACBR',
            'PACABS', 'PAABSF', 'PACBAR', 'ACMODL',

            'CTETRA', 'CPYRAM', 'CPENTA', 'CHEXA',
            'CIHEX1', 'CIHEX2',
            'CSHEAR', 'CVISC', 'CRAC2D', 'CRAC3D',
            'CGAP',
            'GENEL',

            ## rigid_elements
            'RBAR', 'RBAR1', 'RBE1', 'RBE2', 'RBE3', 'RROD', 'RSPLINE', 'RSSCON',

            ## plotels
            'PLOTEL',

            ## properties
            'PMASS',
            'PELAS', 'PGAP', 'PFAST', 'PLPLANE', 'PPLANE',
            'PBUSH', 'PBUSH1D',
            'PDAMP', 'PDAMP5',
            'PROD', 'PBAR', 'PBARL', 'PBEAM', 'PTUBE', 'PBCOMP', 'PBRSECT', 'PBEND',
            'PBEAML', 'PBMSECT', # not fully supported
            'PBEAM3',  # v1.3

            'PSHELL', 'PCOMP', 'PCOMPG', 'PSHEAR',
            'PSOLID', 'PLSOLID', 'PVISC', 'PRAC2D', 'PRAC3D',
            'PIHEX', 'PCOMPS',
            # PQUAD4

            # axixsymmetric
            'CCONEAX', # element
            'PCONEAX', # property
            'AXIC', # axic
            'AXIF', # axif
            'FORCEAX', # loads

            ## pdampt
            'PDAMPT',

            ## pelast
            'PELAST',

            ## pbusht
            'PBUSHT',

            ## creep_materials
            'CREEP',

            ## materials
            'MAT1', 'MAT2', 'MAT3', 'MAT8', 'MAT9', 'MAT10', 'MAT11', 'MAT3D',
            'MATG', 'MATHE', 'MATHP',

            ## Material dependence - MATT1/MATT2/etc.
            'MATT1', 'MATT2', 'MATT3', 'MATT4', 'MATT5', 'MATT8', 'MATT9',
            'MATS1', #'MATS3', 'MATS8',
            # 'MATHE'
            #'EQUIV', # testing only, should never be activated...

            ## nxstrats
            'NXSTRAT',

            ## thermal_materials
            'MAT4', 'MAT5',

            ## spcs
            'SPC', 'SPCADD', 'SPC1', 'SPCAX', 'SPCOFF', 'SPCOFF1',
            'GMSPC',

            ## mpcs
            'MPC', 'MPCADD',

            ## suport/suport1/se_suport
            'SUPORT', 'SUPORT1', 'SESUP',

            ## dloads
            'DLOAD',

            ## dload_entries
            'ACSRCE', 'TLOAD1', 'TLOAD2', 'RLOAD1', 'RLOAD2',
            'QVECT',
            'RANDPS', 'RANDT1', # random

            ## loads
            'LOAD', 'CLOAD', 'LSEQ', 'LOADCYN', 'LOADCYH',
            'SLOAD',
            'FORCE', 'FORCE1', 'FORCE2',
            'MOMENT', 'MOMENT1', 'MOMENT2',
            'GRAV', 'ACCEL', 'ACCEL1',
            'PLOAD', 'PLOAD1', 'PLOAD2', 'PLOAD4',
            'PLOADX1', 'RFORCE', 'RFORCE1',
            'GMLOAD', 'SPCD', 'DEFORM',

            # axisymmetric
            'PRESAX',

            #thermal
            'QVOL',

            # aero cards
            'AERO',  ## aero
            'AEROS',  ## aeros
            'GUST',  ## gusts
            'FLUTTER',   ## flutters
            'FLFACT',   ## flfacts
            'MKAERO1', 'MKAERO2',  ## mkaeros
            'AECOMP', 'AECOMPL',   ## aecomps
            'AEFACT',   ## aefacts
            'AELINK',   ## aelinks
            'AELIST',   ## aelists
            'AEPARM',   ## aeparams
            'AESTAT',   ## aestats
            'AESURF',  ## aesurf
            'AESURFS', ## aesurfs
            'CAERO1', 'CAERO2', 'CAERO3', 'CAERO4', 'CAERO5', ## caeros
            'PAERO1', 'PAERO2', 'PAERO3', 'PAERO4', 'PAERO5', ## paeros

            'MONPNT1', 'MONPNT2', 'MONPNT3', 'MONDSP1', ## monitor_points
            'SPLINE1', 'SPLINE2', 'SPLINE3', 'SPLINE4', 'SPLINE5',  ## splines
            'SPLINE6', 'SPLINE7',
            'TRIM', 'TRIM2',  ## trims
            'CSSCHD',  ## csschds
            'DIVERG',  ## divergs

            ## coords
            'CORD1R', 'CORD1C', 'CORD1S',
            'CORD2R', 'CORD2C', 'CORD2S',
            'GMCORD',

            # temperature cards
            'TEMP', 'TEMPD', 'TEMPB3', 'TEMPAX',
            'QBDY1', 'QBDY2', 'QBDY3', 'QHBDY',
            'CHBDYE', 'CHBDYG', 'CHBDYP',
            'PCONV', 'PCONVM', 'PHBDY',
            'RADBC', 'CONV',
            'RADM', 'VIEW', 'VIEW3D',  # TODO: not validated


            'RADCAV', ## radcavs

            # ---- dynamic cards ---- #
            'DAREA',  ## dareas
            'DPHASE',  ## dphases
            'DELAY',  ## delays
            'NLPARM',  ## nlparms
            'ROTORG', 'ROTORD', ## rotors
            'NLPCI',  ## nlpcis
            'TSTEP',  ## tsteps
            'TSTEPNL', 'TSTEP1',  ## tstepnls
            'TF',  ## transfer_functions
            'TIC', ## initial conditions - sid (set ID)

            ## frequencies
            'FREQ', 'FREQ1', 'FREQ2', 'FREQ3', 'FREQ4', 'FREQ5',

            # direct matrix input cards
            'DMIG', 'DMIJ', 'DMIJI', 'DMIK', 'DMI', 'DTI',
            'DMIAX',

            # optimization cards
            'DEQATN', 'DTABLE',
            'DCONSTR', 'DESVAR', 'TOPVAR', 'DDVAL', 'DRESP1', 'DRESP2', 'DRESP3',
            'DVCREL1', 'DVCREL2',
            'DVPREL1', 'DVPREL2',
            'DVMREL1', 'DVMREL2',
            'DOPTPRM', 'DLINK', 'DCONADD', 'DVGRID',
            'DSCREEN',

            # sets
            'SET1', 'SET3',  ## sets
            'ASET', 'ASET1',  ## asets
            'OMIT', 'OMIT1',  ## omits
            'BSET', 'BSET1',  ## bsets
            'CSET', 'CSET1',  ## csets
            'QSET', 'QSET1',  ## qsets
            'USET', 'USET1',  ## usets

            'RADSET',  # radset

            # superelements
            'SETREE', 'SENQSET', 'SEBULK', 'SEBNDRY', 'SEELT', 'SELOC', 'SEMPLN',
            'SECONCT', 'SELABEL', 'SEEXCLD', 'CSUPER', 'CSUPEXT',
            'SELOAD', 'RELEASE',

            # super-element sets
            'SESET',  ## se_sets

            'SEBSET', 'SEBSET1',  ## se_bsets
            'SECSET', 'SECSET1',  ## se_csets
            'SEQSET', 'SEQSET1',  ## se_qsets
            #'SEUSET', 'SEUSET1',  ## se_usets
            'SEQSEP',

            #------------------------------------------------------------------
            ## parametric
            'PSET', 'PVAL', 'GMCURV', 'GMSURF', 'FEEDGE', 'FEFACE',

            #------------------------------------------------------------------
            ## tables
            'TABLED1', 'TABLED2', 'TABLED3', 'TABLED4',  # dynamic tables - freq/time loads
            'TABLEM1', 'TABLEM2', 'TABLEM3', 'TABLEM4',  # material tables - temperature

            # nonlinear elastic temperature dependent materials (e.g. creep)
            # sees TABLES1
            'TABLEST',
            # material tables - stress (MATS1, CREEP, MATHP)
            'TABLES1',

            ## modal damping table - tables_sdamping
            'TABDMP1',

            ## random_tables
            # PSD=func(freq); used by RANDPS card
            'TABRND1',
            # gust for aeroelastic response; used by RANDPS card
            'TABRNDG',

            # ???
            'TABLEHT', 'TABLEH1',

            #------------------------------------------------------------------
            #: methods
            'EIGB', 'EIGR', 'EIGRL',

            #: cMethods
            'EIGC', 'EIGP',

            # : modtrak
            'MODTRAK',

            #: contact
            'BCTPARA',  ## bctpara
            'BCRPARA',  ## bcrpara
            'BCTADD',  ## bctadds
            'BCTSET',  ## bctsets
            'BSURF',  ## bsurf
            'BSURFS',  ## bsurfs
            'BCONP', ## bconp
            'BLSEG', ## blseg
            'BFRIC', ## bfric

            'TEMPBC',
            #'RADMT',
            'RADLST', 'RADMTX', #'RADBND',
            #'TEMPP1',
            'TEMPRB',
            'CONVM',
            ## ???
            #'PANEL', 'SWLDPRM',
            #'CWELD', 'PWELD', 'PWSEAM', 'CWSEAM', 'CSEAM', 'PSEAM', 'DVSHAP', 'BNDGRID',
            #'CYSYM', 'CYJOIN', 'MODTRAK', 'DSCONS', 'DVAR', 'DVSET', 'DYNRED',
            #'BNDFIX', 'BNDFIX1',
            #'AEFORCE', 'UXVEC', 'GUST2',

            # cyclic
            'CYJOIN', 'CYAX',

            # other
            'INCLUDE',  # '='
            'ENDDATA',
        ]
        set_cards_to_read = set(cards_to_read)
        from collections import Counter
        if len(cards_to_read) != len(set_cards_to_read):  # pragma: no cover
            bad_cards = [key for key, value in Counter(cards_to_read).items()
                         if value > 1]
            raise RuntimeError(bad_cards)

        # the list of possible cards that will be parsed
        self.cards_to_read = set_cards_to_read

        self._xref = False

        #case_control_cards = {'FREQ', 'GUST', 'MPC', 'SPC', 'NLPARM', 'NSM',
                              #'TEMP', 'TSTEPNL', 'INCLUDE'}
        #self._unique_bulk_data_cards = self.cards_to_read.difference(CASE_CONTROL_CARDS)

        #: / is the delete from restart card
        self.special_cards = ['DEQATN', '/']
        self._make_card_parser()

        self._nastran_format = mode
        if mode == 'msc':
            self.set_as_msc()
        elif mode == 'nx':
            self.set_as_nx()
        elif mode == 'zona':
            self.set_as_zona()
        elif mode == 'nasa95':
            self.set_as_nasa95()
        else:  # pragma: no cover
            msg = 'mode=%r is not supported; modes=[msc, nx, zona]' % self._nastran_format
            raise NotImplementedError(msg)

    def __getstate__(self):
        """clears out a few variables in order to pickle the object"""
        # Copy the object's state from self.__dict__ which contains
        # all our instance attributes. Always use the dict.copy()
        # method to avoid modifying the original state.
        state = self.__dict__.copy()
        # Remove the unpicklable entries.
        del state['_card_parser'], state['log']
        if hasattr(self, '_card_parser_b'):
            del state['_card_parser_b']
        if hasattr(self, '_card_parser_prepare'):
            del state['_card_parser_prepare']
        return state

    def get_h5attrs(self) -> List[str]:
        """helper method for dict_to_h5py"""
        attrs = self.object_attributes(mode='both', keys_to_skip=None)
        return attrs

    def export_hdf5_filename(self, hdf5_filename:str) -> None:
        """
        Converts the BDF objects into hdf5 object

        Parameters
        ----------
        hdf5_filename : str
            the path to the hdf5 file

        TODO: doesn't support:
          - BucklingEigenvalues

        """
        import h5py
        with h5py.File(hdf5_filename, 'w') as hdf5_file:
            #self.log.info('starting export_hdf5_file of %r' % hdf5_filename)
            self.export_hdf5_file(hdf5_file)

    def export_hdf5_file(self, hdf5_file, exporter=None) -> None:
        """
        Converts the BDF objects into hdf5 object

        Parameters
        ----------
        hdf5_file : H5File()
            an h5py object
        exporter : HDF5Exporter; default=None
            unused

        """
        from pyNastran.bdf.bdf_interface.hdf5_exporter import export_bdf_to_hdf5_file
        export_bdf_to_hdf5_file(hdf5_file, self)

    def load_hdf5_filename(self, hdf5_filename:str) -> None:
        """
        Loads a BDF object from an hdf5 filename

        Parameters
        ----------
        hdf5_filename : str
            the path to the hdf5 file

        """
        import h5py
        with h5py.File(hdf5_filename, 'r') as hdf5_file:
            #self.log.info('starting load_hdf5_file of %r' % hdf5_filename)
            self.load_hdf5_file(hdf5_file)

    def load_hdf5_file(self, h5_file) -> None:
        """
        Loads a BDF object from an hdf5 object

        Parameters
        ----------
        hdf5_file : H5File()
            an h5py object
        exporter : HDF5Exporter; default=None
            unused

        """
        from pyNastran.bdf.bdf_interface.hdf5_loader import load_bdf_from_hdf5_file
        load_bdf_from_hdf5_file(h5_file, self)

    def saves(self, unxref: bool=True) -> str:
        """Saves a pickled string"""
        if unxref:
            self.uncross_reference()
        return dumps(self)

    def save(self, obj_filename: str='model.obj', unxref: bool=True) -> None:
        """Saves a pickleable object"""
        #del self.log
        #del self._card_parser, self._card_parser_prepare

        #try:
            #del self.log
        #except AttributeError:
            #pass
        #self.case_control_lines = str(self.case_control_deck).split('\n')
        #del self.case_control_deck

        if unxref:
            self.uncross_reference()
        with open(obj_filename, 'wb') as obj_file:
            dump(self, obj_file)

    def load(self, obj_filename: str='model.obj') -> None:
        """Loads a pickleable object"""
        #del self.log
        #lines = print(self.case_control_deck)
        #self.case_control_lines = lines.split('\n')
        #del self.case_control_deck
        #self.uncross_reference()
        #import types
        with open(obj_filename, 'rb') as obj_file:
            obj = load(obj_file)

        # these are properties, functions, etc.
        keys_to_skip = [
            'case_control_deck',
            'log',
            'node_ids', 'coord_ids', 'element_ids', 'property_ids',
            'material_ids', 'caero_ids', 'is_long_ids',
            'nnodes', 'npoints', 'ncoords', 'nelements', 'nproperties',
            'nmaterials', 'ncaeros', 'nid_map',
            'is_bdf_vectorized', 'type_slot_str',
            #'dmigs', 'dmijs', 'dmiks', 'dmijis', 'dtis', 'dmis',

            'point_ids', 'subcases',
            '_card_parser', '_card_parser_b', '_card_parser_prepare',
            'wtmass',
        ]
        for key in object_attributes(self, mode='all', keys_to_skip=keys_to_skip):
            if key.startswith('__') and key.endswith('__'):
                continue

            val = getattr(obj, key)
            #print(key)
            #if isinstance(val, types.FunctionType):
                #continue
            try:
                setattr(self, key, val)
            except AttributeError:  # pragma: no cover
                raise AttributeError('key=%r val=%s\nupdate ~line 860 of bdf.py and '
                                     'add the new key (%s)' % (key, val, key))

        self.case_control_deck = CaseControlDeck(self.case_control_lines, log=self.log)
        #self.log.debug('done loading!')
        for model in self.superelement_models.values():
            model.log = self.log

    def replace_cards(self, replace_model) -> None:
        """
        Replaces the common cards from the current (self) model from the
        ones in the new replace_model.  The intention is that you're
        going to replace things like PSHELLs and DESVARs from a pch file
        in order to update your BDF with the optimized geometry.

        .. todo:: only does a subset of cards.

        Notes
        -----
        loads/spcs (not supported) are tricky because you can't replace
        cards one-to-one...not sure what to do

        """
        for nid, node in replace_model.nodes.items():
            self.nodes[nid] = node
        for eid, elem in self.elements.items():
            self.elements[eid] = elem
        for eid, elem in replace_model.rigid_elements.items():
            self.rigid_elements[eid] = elem
        for pid, prop in replace_model.properties.items():
            self.properties[pid] = prop
        for mid, mat in replace_model.materials.items():
            self.materials[mid] = mat

        for dvid, desvar in replace_model.desvars.items():
            self.desvars[dvid] = desvar
        for dvid, dvprel in replace_model.dvprels.items():
            self.dvprels[dvid] = dvprel
        for dvid, dvmrel in replace_model.dvmrels.items():
            self.dvmrels[dvid] = dvmrel
        for dvid, dvgrid in replace_model.dvgrids.items():
            self.dvgrids[dvid] = dvgrid

    def disable_cards(self, cards: Sequence[str]) -> None:
        """
        Method for removing broken cards from the reader

        Parameters
        ----------
        cards : List[str]; Set[str]
            a list/set of cards that should not be read

        .. python ::

            bdf_model.disable_cards(['DMIG', 'PCOMP'])

        """
        if cards is None:
            return
        elif isinstance(cards, str):
            disable_set = set([cards])
        else:
            disable_set = set(cards)
        self.cards_to_read = self.cards_to_read.difference(disable_set)

    def set_cards(self, cards: Union[List[str],Set[str]]) -> None:
        """
        Method for setting the cards that will be processed

        Parameters
        ----------
        cards : List[str]; Set[str]
            a list/set of cards that should not be read

        .. python ::

            bdf_model.set_cards(['GRID', 'CTRIA3'])

        """
        if cards is None:
            return
        elif isinstance(cards, str):
            enable_set = set([cards])
        else:
            enable_set = set(cards)
        self.cards_to_read = enable_set

    def set_error_storage(self, nparse_errors: int=100, stop_on_parsing_error: bool=True,
                          nxref_errors: int=100, stop_on_xref_error: bool=True) -> None:
        """
        Catch parsing errors and store them up to print them out all at once
        (not all errors are caught).

        Parameters
        ----------
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
        assert isinstance(nparse_errors, int), type(nparse_errors)
        assert isinstance(nxref_errors, int), type(nxref_errors)
        self._nparse_errors = nparse_errors
        self._nxref_errors = nxref_errors
        self._stop_on_parsing_error = stop_on_parsing_error
        self._stop_on_xref_error = stop_on_xref_error

    def validate(self) -> None:
        """runs some checks on the input data beyond just type checking"""
        validate_bdf(self)

    def _verify_bdf(self, xref: Optional[bool]=None) -> None:
        """Cross reference verification method."""
        if xref is None:
            xref = self._xref
        verify_bdf(self, xref)

    def include_zip(self, bdf_filename: Optional[str]=None,
                    encoding: Optional[str]=None,
                    make_ilines: bool=True) -> (List[str], Any):
        """
        Read a bdf without perform any other operation, except (optionally)
        insert the INCLUDE files in the bdf

        Parameters
        ----------
        bdf_filename : str / None
            the input bdf (default=None; popup a dialog)
        encoding : str; default=None -> system default
            the unicode encoding
        make_ilines : bool; default=True
            flag for ilines

        Returns
        -------
        all_lines : List[str]
            all the lines packed into a single line stream
        ilines : (nlines, 2) int ndarray
            if make_ilines = True:
                the [ifile, iline] pair for each line in the file
            if make_ilines = False:
                 ilines = None

        .. note::  Setting read_includes to False is kind of pointless if
                   called directly; it's useful for ``read_bdf``

        """
        punch = False #  doesn't really matter
        read_includes = True
        self._read_bdf_helper(bdf_filename, encoding, punch, read_includes)
        self._parse_primary_file_header(bdf_filename)

        obj = BDFInputPy(self.read_includes, self.dumplines, self._encoding,
                         nastran_format=self.nastran_format,
                         consider_superelements=self.is_superelements,
                         log=self.log, debug=self.debug)
        main_lines = obj.get_main_lines(self.bdf_filename)
        all_lines, ilines = obj.lines_to_deck_lines(main_lines, make_ilines=make_ilines)
        self._set_pybdf_attributes(obj, save_file_structure=False)
        return all_lines, ilines

    def _set_pybdf_attributes(self, obj: BDFInputPy,
                              save_file_structure: bool=False) -> None:
        # type: (BDFInputPy, bool) -> None
        """common method for all functions that use BDFInputPy"""
        # these are include line pairs
        #print(obj.include_lines)
        self.active_filenames = []
        self.reject_lines = []
        self.include_filenames = defaultdict(list)
        for ifile, include_lines_filename_pairs in obj.include_lines.items():
            assert len(include_lines_filename_pairs) > 0, include_lines_filename_pairs
            for include_lines, bdf_filename2 in include_lines_filename_pairs:
                #print(ifile, include_lines)
                self.include_filenames[ifile].append(bdf_filename2)
                if not save_file_structure and not obj.read_includes:
                    self.reject_lines += include_lines
        #print('-------------ssett (end)----------')
        self.active_filenames += obj.active_filenames
        self.active_filename = obj.active_filename
        self.include_dir = obj.include_dir

    def read_bdf(self, bdf_filename: Optional[str]=None,
                 validate: bool=True,
                 xref: bool=True,
                 punch: bool=False,
                 read_includes: bool=True,
                 save_file_structure: bool=False,
                 encoding: Optional[str]=None) -> None:
        """
        Read method for the bdf files

        Parameters
        ----------
        bdf_filename : str / None
            the input bdf (default=None; popup a dialog)
        validate : bool; default=True
            runs various checks on the BDF
        xref :  bool; default=True
            should the bdf be cross referenced
        punch : bool; default=False
            indicates whether the file is a punch file
        read_includes : bool; default=True
            indicates whether INCLUDE files should be read
        save_file_structure : bool; default=False
            enables the ``write_bdfs`` method
        encoding : str; default=None -> system default
            the unicode encoding

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
        self.save_file_structure = save_file_structure
        if bdf_filename and not isinstance(bdf_filename, (StringIO, list)):
            check_path(bdf_filename, 'bdf_filename')
        self._read_bdf_helper(bdf_filename, encoding, punch, read_includes)
        self.log.debug('---starting BDF.read_bdf of %s---' % self.bdf_filename)
        self._parse_primary_file_header(bdf_filename)

        obj = BDFInputPy(self.read_includes, self.dumplines, self._encoding,
                         nastran_format=self.nastran_format,
                         consider_superelements=self.is_superelements,
                         log=self.log, debug=self.debug)
        out = obj.get_lines(bdf_filename, punch=self.punch, make_ilines=True)
        system_lines, executive_control_lines, case_control_lines, bulk_data_lines, bulk_data_ilines, superelement_lines, superelement_ilines = out
        self._set_pybdf_attributes(obj, save_file_structure)

        self.system_command_lines = system_lines
        self.executive_control_lines = executive_control_lines
        self.case_control_lines = case_control_lines

        sol, method, sol_iline = parse_executive_control_deck(executive_control_lines)
        self.update_solution(sol, method, sol_iline)

        self.case_control_deck = CaseControlDeck(case_control_lines, self.log)
        self.case_control_deck.solmap_to_value = self._solmap_to_value
        self.case_control_deck.rsolmap_to_str = self.rsolmap_to_str

        try:
            self._parse_all_cards(bulk_data_lines, bulk_data_ilines)
        except SuperelementFlagError:
            if self.is_superelements:
                raise

            self.clear_attributes()
            self.log.error('Attempting to use is_superelements=True.')
            self.is_superelements = True
            self.read_bdf(bdf_filename=bdf_filename, validate=validate, xref=xref, punch=punch,
                          read_includes=read_includes, save_file_structure=save_file_structure,
                          encoding=encoding)
            return

        if superelement_lines:
            self._add_superelements(superelement_lines, superelement_ilines)

        self.pop_parse_errors()
        fill_dmigs(self)

        if validate:
            self.validate()

        self.cross_reference(xref=xref)
        self._xref = xref

        self.log.debug('---finished BDF.read_bdf of %s---' % self.bdf_filename)

    def _add_superelements(self, superelement_lines: List[str],
                           superelement_ilines: Any) -> None:  # pragma: no cover
        self.log.warning('_add_superelements should be overwritten')

    def _parse_all_cards(self, bulk_data_lines: List[str], bulk_data_ilines: Any) -> None:
        """creates and loads all the cards the bulk data section"""
        cards_list = []
        cards_dict = {}
        if self._is_cards_dict:
            cards_dict, card_count = self.get_bdf_cards_dict(
                bulk_data_lines, bulk_data_ilines)
            #if 0:
                #with open('dump.bdf', 'w') as bdf_file_obj:
                    #bdf_file_obj.write('\n'.join(executive_control_lines))
                    #bdf_file_obj.write(str(case_control_deck))
                    #for cardname, cards in cards.items():
                        #for (comment, cardlines) in cards:
                            ##bdf_file_obj.write(comment + '\n')
                            #bdf_file_obj.write('\n'.join(cardlines) + '\n')
                        #bdf_file_obj.write('\n')
        else:
            cards_list, cards_dict, card_count = self.get_bdf_cards(
                bulk_data_lines, bulk_data_ilines)
            #for card in cards_list:
                #card_name = card[0]
                #if card_name == 'CBAR':
                    #print(card)
        self._parse_cards(cards_list, cards_dict, card_count)

        if self.values_to_skip:
            for key, values in self.values_to_skip.items():
                dict_values = getattr(self, key)
                if not isinstance(dict_values, dict):
                    msg = '%r is an invalid type; only dictionaries are supported' % key
                    raise TypeError(msg)
                for value in values:
                    del dict_values[value]
            # TODO: redo get_card_ids_by_card_types & card_count

    def _read_bdf_helper(self, bdf_filename: str, encoding: str,
                         punch: bool, read_includes: bool):
        """creates the file loading if bdf_filename is None"""
        #self.set_error_storage(nparse_errors=None, stop_on_parsing_error=True,
        #                       nxref_errors=None, stop_on_xref_error=True)
        if encoding is None:
            encoding = sys.getdefaultencoding()
        self._encoding = encoding

        self.read_includes = read_includes
        self.active_filenames = []

        if bdf_filename is None:
            from pyNastran.utils.gui_io import load_file_dialog
            wildcard_wx = "Nastran BDF (*.bdf; *.dat; *.nas; *.pch, *.ecd)|" \
                "*.bdf;*.dat;*.nas;*.pch;*.ecd|" \
                "All files (*.*)|*.*"
            wildcard_qt = "Nastran BDF (*.bdf *.dat *.nas *.pch *.ecd);;All files (*)"
            title = 'Please select a BDF/DAT/PCH/ECD to load'
            bdf_filename = load_file_dialog(title, wildcard_wx, wildcard_qt)[0]
            assert bdf_filename is not None, bdf_filename

        elif isinstance(bdf_filename, str):
            pass
        elif isinstance(bdf_filename, (StringIO, IOBase)):
            self.bdf_filename = bdf_filename
            self.punch = punch
            return
        else:
            raise NotImplementedError(bdf_filename)

        check_path(bdf_filename, 'bdf_filename')
        if bdf_filename.lower().endswith('.pch'):  # .. todo:: should this be removed???
            punch = True

        #: the active filename (string)
        self.bdf_filename = bdf_filename

        #: is this a punch file (no executive control deck)
        self.punch = punch
        assert not bdf_filename.lower().endswith('.op2'), bdf_filename

    def pop_parse_errors(self) -> None:
        """raises an error if there are parsing errors"""
        if self._stop_on_parsing_error:
            if self._iparse_errors == 1 and self._nparse_errors == 0:
                raise
            is_error = False
            msg = ''
            if self._duplicate_elements:
                duplicate_eids = [elem.eid for elem in self._duplicate_elements]
                uduplicate_eids = np.unique(duplicate_eids)
                msg += 'self.elements IDs are not unique=%s\n' % uduplicate_eids
                for eid in uduplicate_eids:
                    msg += 'old_element=\n%s\n' % str(self.elements[eid])
                    msg += 'new_elements=\n'
                    for elem, eidi in zip(self._duplicate_elements, duplicate_eids):
                        if eidi == eid:
                            msg += str(elem)
                    msg += '\n'
                    is_error = True
                    raise DuplicateIDsError(msg)

            if self._duplicate_properties:
                duplicate_pids = [prop.pid for prop in self._duplicate_properties]
                uduplicate_pids = np.unique(duplicate_pids)
                msg += 'self.properties IDs are not unique=%s\n' % uduplicate_pids
                for pid in duplicate_pids:
                    msg += 'old_property=\n%s\n' % str(self.properties[pid])
                    msg += 'new_properties=\n'
                    for prop, pidi in zip(self._duplicate_properties, duplicate_pids):
                        if pidi == pid:
                            msg += str(prop)
                    msg += '\n'
                    is_error = True

            if self._duplicate_masses:
                duplicate_eids = [elem.eid for elem in self._duplicate_masses]
                uduplicate_eids = np.unique(duplicate_eids)
                msg += 'self.massses IDs are not unique=%s\n' % uduplicate_eids
                for eid in uduplicate_eids:
                    msg += 'old_mass=\n%s\n' % str(self.masses[eid])
                    msg += 'new_masses=\n'
                    for elem, eidi in zip(self._duplicate_masses, duplicate_eids):
                        if eidi == eid:
                            msg += str(elem)
                    msg += '\n'
                    is_error = True

            if self._duplicate_materials:
                duplicate_mids = [mat.mid for mat in self._duplicate_materials]
                uduplicate_mids = np.unique(duplicate_mids)
                msg += 'self.materials IDs are not unique=%s\n' % uduplicate_mids
                for mid in uduplicate_mids:
                    msg += 'old_material=\n%s\n' % str(self.materials[mid])
                    msg += 'new_materials=\n'
                    for mat, midi in zip(self._duplicate_materials, duplicate_mids):
                        if midi == mid:
                            msg += str(mat)
                    msg += '\n'
                    is_error = True

            if self._duplicate_thermal_materials:
                duplicate_mids = [mat.mid for mat in self._duplicate_thermal_materials]
                uduplicate_mids = np.unique(duplicate_mids)
                msg += 'self.thermal_materials IDs are not unique=%s\n' % uduplicate_mids
                for mid in uduplicate_mids:
                    msg += 'old_thermal_material=\n%s\n' % str(self.thermal_materials[mid])
                    msg += 'new_thermal_materials=\n'
                    for mat, midi in zip(self._duplicate_thermal_materials, duplicate_mids):
                        if midi == mid:
                            msg += str(mat)
                    msg += '\n'
                    is_error = True

            if self._duplicate_coords:
                duplicate_cids = [coord.cid for coord in self._duplicate_coords]
                uduplicate_cids = np.unique(duplicate_cids)
                msg += 'self.coords IDs are not unique=%s\n' % uduplicate_cids
                for cid in uduplicate_cids:
                    msg += 'old_coord=\n%s\n' % str(self.coords[cid])
                    msg += 'new_coords=\n'
                    for coord, cidi in zip(self._duplicate_coords, duplicate_cids):
                        if cidi == cid:
                            msg += str(coord)
                    msg += '\n'
                    is_error = True

            if is_error:
                msg = 'There are duplicate cards.\n\n' + msg

            if self._stop_on_xref_error:
                msg += 'There are parsing errors.\n\n'
                for (card, an_error) in self._stored_parse_errors:
                    msg += '%scard=%s\n' % (an_error[0], card)
                    msg += 'xref error: %s\n\n'% an_error[0]
                    is_error = True

            if is_error:
                print('%s' % msg)
                raise DuplicateIDsError(msg.rstrip())

    def pop_xref_errors(self) -> None:
        """raises an error if there are cross-reference errors"""
        is_error = False
        if self._stop_on_xref_error:
            if self._ixref_errors == 1 and self._nxref_errors == 0:
                raise
            if self._stored_xref_errors:
                msg = 'There are cross-reference errors.\n\n'
                for (card, an_error) in self._stored_xref_errors:
                    msg += '%scard=%s\n' % (an_error[0], card)
                    is_error = True

                if is_error and self._stop_on_xref_error:
                    raise CrossReferenceError(msg.rstrip())

    def get_bdf_cards(self, bulk_data_lines, bulk_data_ilines=None):
        # type: (List[str], Optional[Any]) -> Any, Any, Any
        """Parses the BDF lines into a list of card_lines"""
        if bulk_data_ilines is None:
            bulk_data_ilines = np.zeros((len(bulk_data_lines), 2), dtype='int32')

        cards_list = []
        cards_dict = defaultdict(list)
        dict_cards = ['BAROR', 'BEAMOR']
        #cards = defaultdict(list)
        card_count = defaultdict(int)
        full_comment = ''
        card_lines = []
        old_ifile_iline = None
        old_card_name = None
        backup_comment = ''
        nlines = len(bulk_data_lines)
        if len(bulk_data_lines) != len(bulk_data_ilines):
            msg = 'len(bulk_data_lines)=%s len(bulk_data_ilines)=%s' % (
                len(bulk_data_lines), len(bulk_data_ilines))
            self.log.warning(msg)

        for iline_bulk, line in enumerate(bulk_data_lines):
            ifile_iline = bulk_data_ilines[iline_bulk, :]
            #print(iline_bulk, ifile_iline, line)
            #print('    backup=%r' % backup_comment)
            comment = ''
            if '$' in line:
                line, comment = line.split('$', 1)
            card_name = line.split(',', 1)[0].split('\t', 1)[0][:8].rstrip().upper()
            if card_name and card_name[0] not in ['+', '*']:
                if old_card_name:
                    if self.echo and not self.force_echo_off:
                        self.log.info('Reading %s:\n' %
                                      old_card_name + full_comment + ''.join(card_lines))

                    # old dictionary version
                    # cards_list[old_card_name].append([full_comment, card_lines, ifile_iline])

                    # new list version
                    #if full_comment:
                        #print('full_comment = ', full_comment)
                    if old_card_name in dict_cards:
                        cards_dict[old_card_name].append([_prep_comment(full_comment),
                                                          card_lines, ifile_iline])
                    else:
                        cards_list.append([old_card_name, _prep_comment(full_comment),
                                           card_lines, old_ifile_iline])

                    card_count[old_card_name] += 1
                    card_lines = []
                    full_comment = ''

                    if old_card_name == 'ECHOON':
                        self.echo = True
                    elif old_card_name == 'ECHOOFF':
                        self.echo = False
                old_ifile_iline = ifile_iline
                old_card_name = card_name.rstrip(' *')

                if old_card_name == 'ENDDATA':
                    self.card_count['ENDDATA'] = 1
                    if nlines - iline_bulk > 1:
                        nleftover = nlines - iline_bulk - 1
                        msg = 'exiting due to ENDDATA found with %i lines left' % nleftover
                        self.log.debug(msg)
                    return cards_list, cards_dict, card_count
                #print("card_name = %s" % card_name)

            comment = _clean_comment(comment)

            #TODO: these additional \n need to be there for rejected cards
            #      but not parsed cards
            if line.rstrip():
                card_lines.append(line)
                if backup_comment:
                    if comment:
                        full_comment += backup_comment + comment + '\n'
                    else:
                        full_comment += backup_comment
                    backup_comment = ''
                elif comment:
                    full_comment += comment + '\n'
                    backup_comment = ''

            elif comment:
                backup_comment += comment + '\n'
            #elif comment:
                #backup_comment += '$' + comment + '\n'


        if card_lines:
            if self.echo and not self.force_echo_off:
                self.log.info('Reading %s:\n' % old_card_name + full_comment + ''.join(card_lines))
            #print('end_add %s' % card_lines)

            # old dictionary version
            #cards[old_card_name].append([backup_comment + full_comment, card_lines])

            # new list version
            #if backup_comment + full_comment:
                #print('backup_comment + full_comment = ', backup_comment + full_comment)
            if old_card_name in dict_cards:
                cards_dict[old_card_name].append([_prep_comment(
                    backup_comment + full_comment), card_lines, ifile_iline])
            else:
                cards_list.append([old_card_name, _prep_comment(
                    backup_comment + full_comment), card_lines, ifile_iline])
            card_count[old_card_name] += 1
        self.echo = False
        return cards_list, cards_dict, card_count

    def get_bdf_cards_dict(self, bulk_data_lines, bulk_data_ilines=None):
        """Parses the BDF lines into a list of card_lines"""
        if bulk_data_ilines is None:
            bulk_data_ilines = np.zeros((len(bulk_data_lines), 2), dtype='int32')

        cards_dict = defaultdict(list)
        card_count = defaultdict(int)
        full_comment = ''
        card_lines = []
        old_card_name = None
        backup_comment = ''
        nlines = len(bulk_data_lines)
        for iline_bulk, line in enumerate(bulk_data_lines):
            ifile_iline = bulk_data_ilines[iline_bulk, :]
            #print('    backup=%r' % backup_comment)
            comment = ''
            if '$' in line:
                line, comment = line.split('$', 1)
            card_name = line.split(',', 1)[0].split('\t', 1)[0][:8].rstrip().upper()
            if card_name and card_name[0] not in ['+', '*']:
                if old_card_name:
                    if self.echo and not self.force_echo_off:
                        self.log.info('Reading %s:\n' %
                                      old_card_name + full_comment + ''.join(card_lines))

                    # old dictionary version
                    cards_dict[old_card_name].append([full_comment, card_lines, ifile_iline])

                    # new list version
                    #cards.append([old_card_name, full_comment, card_lines, ifile_iline])

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
                    if nlines - iline_bulk > 1:
                        nleftover = nlines - iline_bulk - 1
                        msg = 'exiting due to ENDDATA found with %i lines left' % nleftover
                        self.log.debug(msg)
                    return cards_dict, card_count
                #print("card_name = %s" % card_name)

            comment = _clean_comment_bulk(comment)
            if line.rstrip():
                card_lines.append(line)
                if backup_comment:
                    if comment:
                        full_comment += backup_comment + comment + '\n'
                    else:
                        full_comment += backup_comment
                    backup_comment = ''
                elif comment:
                    full_comment += comment + '\n'
                    backup_comment = ''

            elif comment:
                backup_comment += comment + '\n'
                #print('add backup=%r' % backup_comment)
            #elif comment:
                #backup_comment += comment + '\n'

        if card_lines:
            if self.echo and not self.force_echo_off:
                self.log.info('Reading %s:\n' % old_card_name + full_comment + ''.join(card_lines))
            #print('end_add %s' % card_lines)

            # old dictionary version
            cards_dict[old_card_name].append(
                [backup_comment + full_comment, card_lines, ifile_iline])

            # new list version
            #cards.append([old_card_name, backup_comment + full_comment, card_lines])
            card_count[old_card_name] += 1
        return cards_dict, card_count

    def update_solution(self, sol: int,
                        method: Optional[str],
                        sol_iline: int) -> None:
        """
        Updates the overall solution type (e.g. 101,200,600)

        Parameters
        ----------
        sol : int
            the solution type (101, 103, etc)
        method : str
            the solution method (only for SOL=600)
        sol_iline : int
            the line to put the SOL/method on

        """
        self.sol_iline = sol_iline
        # the integer of the solution type (e.g. SOL 101)
        if sol is None:
            self.sol = None
            self.sol_method = None
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
            if method is None:
                method = ''
            self.sol_method = method.strip()
            self.log.debug("sol=%s method=%s" % (self.sol, self.sol_method))
        else:  # very common
            self.sol_method = None

    def update_card(self, card_name: str, icard: int, ifield: int,
                    value: Union[int, float, str]) -> None:
        """
        Updates a Nastran card based on standard Nastran optimization names

        Parameters
        ----------
        card_name : str
            the name of the card
            (e.g. GRID)
        icard : int
            the unique 1-based index identifier for the card
            (e.g. the GRID id)
        ifield : int
            the index on the card
            (e.g. X on GRID card as an integer representing the field number)
        value : varies
            the value to assign

        Returns
        -------
        obj : varies
            the corresponding object
            (e.g. the GRID object)

        # On GRID 100, set Cp (2) to 42
        >>> model.update_card('GRID', 100, 2, 42)

        # On GRID 100, set X (3) to 43.
        >>> model.update_card('GRID', 100, 3, 43.)

        """
        #rslot_map = self.get_rslot_map(reset_type_to_slot_map=False)

        for key in self.card_count:
            assert isinstance(key, str), 'key=%r' % key
            if key not in self._type_to_slot_map:
                msg = 'add %r to self._type_to_slot_map\n%s' % (key, str(self._type_to_slot_map))
                raise RuntimeError(msg)
        #_slot_to_type_map['nodes'] : ['GRID']
        #_type_to_slot_map['GRID'] : ['nodes']

        # get the storage object
        try:
            field_str = self._type_to_slot_map[card_name] # 'nodes'
        except KeyError:
            msg = 'Updating card card_name=%r is not supported\nkeys=%s' % (
                card_name, list(self._type_to_slot_map.keys()))
            raise KeyError(msg)

        objs = getattr(self, field_str) # self.nodes
        # get the specific card
        try:
            obj = objs[icard]
        except KeyError:
            msg = 'Could not find %s ID=%r' % (card_name, icard)
            raise KeyError(msg)

        # update the card
        obj.update_field(ifield, value)
        return obj

    def set_dynamic_syntax(self, dict_of_vars: Dict[str, Union[int, float, str]]) -> None:
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
           >>> bdf.read_bdf(bdf_filename, xref=True)

        Notes
        -----
        Case sensitivity is supported.
        Variables should be 7 characters or less to fit in an
        8-character field.

        .. warning:: Type matters!

        """
        self.dict_of_vars = {}
        assert len(dict_of_vars) > 0, 'nvars = %s' % len(dict_of_vars)
        for (key, value) in sorted(dict_of_vars.items()):
            assert len(key) <= 7, ('max length for key is 7; '
                                   'len(%s)=%s' % (key, len(key)))
            assert len(key) >= 1, ('min length for key is 1; '
                                   'len(%s)=%s' % (key, len(key)))
            if not isinstance(key, str):
                msg = 'key=%r must be a string.  type=%s' % (key, type(key))
                raise TypeError(msg)
            self.dict_of_vars[key] = value
        self._is_dynamic_syntax = True

    def is_reject(self, card_name: str) -> bool:
        """
        Can the card be read.

        If the card is rejected, it's added to self.reject_count

        Parameters
        ----------
        card_name : str
            the card_name -> 'GRID'

        """
        assert '=' not in card_name, card_name
        if card_name.startswith('='):
            return False
        elif card_name in self.cards_to_read:
            return False
        if card_name:
            if card_name not in self.reject_count:
                self.reject_count[card_name] = 0
            self.reject_count[card_name] += 1
        return True

    def _process_card(self, card_lines: List[str]) -> List[str]:
        """
        Converts card_lines into a card.
        Considers dynamic syntax and removes empty fields

        Parameters
        ----------
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
            >>> fields, card_name = model._process_card(card_lines)
            >>> fields
            ['GRID', '1', '', '1.0', '2.0', '3.0']
            >>> card_name
            'GRID'

        """
        card_name = _get_card_name(card_lines, self.active_filename)
        fields = to_fields(card_lines, card_name)
        if self._is_dynamic_syntax:
            fields = [
                print_field_16(_parse_dynamic_syntax(field, self.dict_of_vars, self.log))
                if '%' in field[0:1] else field
                for field in fields]
        card = wipe_empty_fields(fields)
        card[0] = card_name
        return card

    def create_card_object(self, card_lines, card_name, is_list=True, has_none=True):
        """
        Creates a BDFCard object, which is really just a list that
        allows indexing past the last field

        Parameters
        ----------
        card_lines: list[str]
            the list of the card fields
            input is list of card_lines -> ['GRID, 1, 2, 3.0, 4.0, 5.0']
        card_name : str
            the card_name -> 'GRID'
        is_list : bool; default=True
            True : this is a list of fields
            False : this is a list of lines
        has_none : bool; default=True
            can there be trailing Nones in the card data (e.g. ['GRID, 1, 2, 3.0, 4.0, 5.0, '])

        Returns
        -------
        card_object : BDFCard()
            the card object representation of card
        card : list[str]
            the card with empty fields removed

        """
        card_name = card_name.upper()
        self.increase_card_count(card_name)
        if card_name in ['DEQATN', 'PBRSECT', 'PBMSECT', 'GMCURV', 'GMSURF', 'OUTPUT', 'ADAPT']:
            card_obj = card_lines
            card = card_lines
        else:
            if is_list:
                fields = card_lines
            else:
                fields = to_fields(card_lines, card_name)

            # apply OPENMDAO syntax
            if self._is_dynamic_syntax:
                fields = [
                    print_field_16(_parse_dynamic_syntax(field, self.dict_of_vars, self.log))
                    if '%' in field.strip()[0:1] else print_field_16(field)
                    for field in fields]
                has_none = False

            if has_none:
                card = wipe_empty_fields([print_field_16(field) for field in fields])
            else:
                #card = remove_trailing_fields(fields)
                card = wipe_empty_fields(fields)
            card_obj = BDFCard(card, has_none=False)
        return card_obj, card

    def _parse_dynamic_syntax(self, key: str) -> Dict[str, Any]:
        return _parse_dynamic_syntax(key, self.dict_of_vars, self.log)

    def _make_card_parser(self) -> None:
        """creates the card parser variables that are used by add_card"""
        class Crash:
            """class for crashing on specific cards"""
            def __init__(self) -> None:
                """dummy init"""
                pass
            @classmethod
            def add_card(cls, card, comment=''):
                """the method that forces the crash"""
                #raise CardParseSyntaxError(card)
                msg = _format_comment(comment) + str(card)
                raise UnsupportedCard(msg)

        #class CrashIgnore:
            #"""class for crashing on specific cards"""
            #def __init__(self):
                #"""dummy init"""
                #pass
            #@classmethod
            #def add_card(cls, card, comment=''):
                #"""the method that forces the crash"""
                ##raise CardParseSyntaxError(card)
                #msg = _format_comment(comment) + str(card)
                #raise DisabledCardError(msg)

        #: a storage of card_name to (card_class, add_method)
        self._card_parser = {
            #'=' : (Crash, None),
            '/' : (Crash, None),

            #'CGEN' : (CrashIgnore, None),

            'SETREE' : (SETREE, self._add_setree_object),
            'SENQSET' : (SENQSET, self._add_senqset_object),
            'SEBULK' : (SEBULK, self._add_sebulk_object),
            'RELEASE': (RELEASE, self._add_release_object),
            'SEBNDRY' : (SEBNDRY, self._add_sebndry_object),
            'SEELT' : (SEELT, self._add_seelt_object),
            'SELOC' : (SELOC, self._add_seloc_object),
            'SEMPLN' : (SEMPLN, self._add_sempln_object),
            'SECONCT' : (SECONCT, self._add_seconct_object),
            'SELABEL' : (SELABEL, self._add_selabel_object),
            'SEEXCLD' : (SEEXCLD, self._add_seexcld_object),
            'CSUPER' : (CSUPER, self._add_csuper_object),
            'CSUPEXT' : (CSUPEXT, self._add_csupext_object),
            'SELOAD' : (SELOAD, self._add_seload_object),

            'CHACAB': (CHACAB, self._add_element_object),
            'CHACBR': (CHACBR, self._add_element_object),
            'CAABSF': (CAABSF, self._add_element_object),
            'PACABS': (PACABS, self._add_acoustic_property_object),
            'PAABSF': (PAABSF, self._add_acoustic_property_object),
            'PACBAR': (PACBAR, self._add_acoustic_property_object),
            #'PANEL' : (Crash, None),

            'BCONP' : (BCONP, self._add_bconp_object),
            'BLSEG' : (BLSEG, self._add_blseg_object),
            'BFRIC' : (BFRIC, self._add_bfric_object),
            'MODTRAK' : (MODTRAK, self._add_modtrak_object),

            #'BGADD', 'BGSET', 'BOLT', 'BOLTFOR'
            'BGADD' : (Crash, None),
            'BGSET' : (Crash, None),
            'BOLT' : (Crash, None),
            'BOLTFOR' : (Crash, None),

            #'CBEAR', 'PBEAR', 'ROTORB',
            'CBEAR' : (Crash, None),
            'PBEAR' : (Crash, None),
            'ROTORB' : (Crash, None),

            #'SWLDPRM' : (Crash, None),

            #'CWELD' : (Crash, None),
            #'PWELD' : (Crash, None),
            #'PWSEAM' : (Crash, None),
            #'CWSEAM' : (Crash, None),
            #'CSEAM' : (Crash, None),
            #'PSEAM' : (Crash, None),

            #'DVSHAP' : (Crash, None),
            #'BNDGRID' : (Crash, None),

            #'CYSYM' : (Crash, None),
            #'TEMPP1' : (Crash, None),
            #'DSCONS' : (Crash, None),
            #'DVAR' : (Crash, None),
            #'DVSET' : (Crash, None),
            #'DYNRED' : (Crash, None),
            #'BNDFIX' : (Crash, None),
            #'BNDFIX1' : (Crash, None),

            #'AEFORCE' : (Crash, None),
            #'UXVEC' : (Crash, None),
            'GUST2' : (Crash, None),

            #'RADBND' : (Crash, None),


            # nodes
            'GRID' : (GRID, self._add_node_object),
            'SPOINT' : (SPOINTs, self._add_spoint_object),
            'EPOINT' : (EPOINTs, self._add_epoint_object),
            'RINGAX' : (RINGAX, self._add_ringax_object),
            'POINTAX' : (POINTAX, self._add_ringax_object),
            'POINT' : (POINT, self._add_point_object),
            'SEQGP' : (SEQGP, self._add_seqgp_object),
            'GRIDB' : (GRIDB, self._add_gridb_object),

            'PARAM' : (PARAM, self._add_param_object),

            'CORD2R' : (CORD2R, self._add_coord_object),
            'CORD2C' : (CORD2C, self._add_coord_object),
            'CORD2S' : (CORD2S, self._add_coord_object),

            # parametric
            'PSET' : (PSET, self._add_pset),
            'PVAL' : (PVAL, self._add_pval),
            'GMCURV' : (GMCURV, self._add_gmcurv),
            'GMSURF' : (GMSURF, self._add_gmsurf),
            'FEFACE' : (FEFACE, self._add_feface),
            'FEEDGE' : (FEEDGE, self._add_feedge),

            # msgmesh
            'GMCORD' : (GMCORD, self._add_coord_object),
            'CGEN' : (CGEN, self._add_element_object),

            'CONROD' : (CONROD, self._add_element_object),
            'CROD' : (CROD, self._add_element_object),
            'PROD' : (PROD, self._add_property_object),
            'CTUBE' : (CTUBE, self._add_element_object),
            'PTUBE' : (PTUBE, self._add_property_object),

            'BAROR' : (BAROR, self._add_baror_object),
            'CBARAO' : (CBARAO, self._add_ao_object),
            'PBAR' : (PBAR, self._add_property_object),
            'PBARL' : (PBARL, self._add_property_object),
            'PBRSECT' : (PBRSECT, self._add_property_object),

            'BEAMOR' : (BEAMOR, self._add_beamor_object),
            'PBEAM' : (PBEAM, self._add_property_object),
            'PBEAML' : (PBEAML, self._add_property_object),
            'PBCOMP' : (PBCOMP, self._add_property_object),
            'PBMSECT' : (PBMSECT, self._add_property_object),

            'CBEAM3' : (CBEAM3, self._add_element_object),
            'PBEAM3' : (PBEAM3, self._add_property_object),

            'CBEND' : (CBEND, self._add_element_object),
            'PBEND' : (PBEND, self._add_property_object),

            'CTRIA3' : (CTRIA3, self._add_element_object),
            'CQUAD4' : (CQUAD4, self._add_element_object),
            'CQUAD' : (CQUAD, self._add_element_object),
            'CQUAD8' : (CQUAD8, self._add_element_object),
            'CQUADX' : (CQUADX, self._add_element_object),
            'CQUADX4' : (CQUADX4, self._add_element_object),
            'CQUADX8' : (CQUADX8, self._add_element_object),
            'CQUADR' : (CQUADR, self._add_element_object),
            'CTRIA6' : (CTRIA6, self._add_element_object),
            'CTRIAR' : (CTRIAR, self._add_element_object),
            'CTRAX3' : (CTRAX3, self._add_element_object),
            'CTRAX6' : (CTRAX6, self._add_element_object),
            'CTRIAX' : (CTRIAX, self._add_element_object),
            'CTRIAX6' : (CTRIAX6, self._add_element_object),
            'SNORM' : (SNORM, self._add_normal_object),
            'PCOMP' : (PCOMP, self._add_property_object),
            'PCOMPG' : (PCOMPG, self._add_property_object),
            'PSHELL' : (PSHELL, self._add_property_object),
            'PLPLANE' : (PLPLANE, self._add_property_object),
            'CPLSTN3' : (CPLSTN3, self._add_element_object),
            'CPLSTN4' : (CPLSTN4, self._add_element_object),
            'CPLSTN6' : (CPLSTN6, self._add_element_object),
            'CPLSTN8' : (CPLSTN8, self._add_element_object),
            'CPLSTS3' : (CPLSTS3, self._add_element_object),
            'CPLSTS4' : (CPLSTS4, self._add_element_object),
            'CPLSTS6' : (CPLSTS6, self._add_element_object),
            'CPLSTS8' : (CPLSTS8, self._add_element_object),
            'PPLANE' : (PPLANE, self._add_property_object),

            'CSHEAR' : (CSHEAR, self._add_element_object),
            'PSHEAR' : (PSHEAR, self._add_property_object),

            'CIHEX1' : (CIHEX1, self._add_element_object),
            'CIHEX2' : (CIHEX2, self._add_element_object),
            'PIHEX' : (PIHEX, self._add_property_object),
            'PSOLID' : (PSOLID, self._add_property_object),
            'PLSOLID' : (PLSOLID, self._add_property_object),
            'PCOMPS' : (PCOMPS, self._add_property_object),

            'CELAS1' : (CELAS1, self._add_element_object),
            'CELAS2' : (CELAS2, self._add_element_object),
            'CELAS3' : (CELAS3, self._add_element_object),
            'CELAS4' : (CELAS4, self._add_element_object),
            'CVISC' : (CVISC, self._add_element_object),
            'PELAST' : (PELAST, self._add_pelast_object),

            'CDAMP1' : (CDAMP1, self._add_damper_object),
            'CDAMP2' : (CDAMP2, self._add_damper_object),
            'CDAMP3' : (CDAMP3, self._add_damper_object),
            # CDAMP4 added later because the documentation is wrong
            'CDAMP5' : (CDAMP5, self._add_damper_object),
            'PDAMP5' : (PDAMP5, self._add_property_object),

            'CFAST' : (CFAST, self._add_damper_object),
            'PFAST' : (PFAST, self._add_property_object),

            'CGAP' : (CGAP, self._add_element_object),
            'PGAP' : (PGAP, self._add_property_object),

            'CBUSH' : (CBUSH, self._add_damper_object),
            'CBUSH1D' : (CBUSH1D, self._add_damper_object),
            'CBUSH2D' : (CBUSH2D, self._add_damper_object),
            'PBUSH' : (PBUSH, self._add_property_object),
            'PBUSH1D' : (PBUSH1D, self._add_property_object),

            'CRAC2D' : (CRAC2D, self._add_element_object),
            'PRAC2D' : (PRAC2D, self._add_property_object),

            'CRAC3D' : (CRAC3D, self._add_element_object),
            'PRAC3D' : (PRAC3D, self._add_property_object),

            'GENEL' : (GENEL, self._add_element_object),


            'PDAMPT' : (PDAMPT, self._add_pdampt_object),
            'PBUSHT' : (PBUSHT, self._add_pbusht_object),

            'CCONEAX' : (CCONEAX, self._add_element_object),
            'PCONEAX' : (PCONEAX, self._add_property_object),
            'AXIC' : (AXIC, self._add_axic_object),
            'AXIF' : (AXIF, self._add_axif_object),
            'CYAX' : (CYAX, self._add_cyax_object),

            'RBAR' : (RBAR, self._add_rigid_element_object),
            'RBAR1' : (RBAR1, self._add_rigid_element_object),
            'RBE1' : (RBE1, self._add_rigid_element_object),
            'RBE2' : (RBE2, self._add_rigid_element_object),
            'RBE3' : (RBE3, self._add_rigid_element_object),
            'RROD' : (RROD, self._add_rigid_element_object),
            'RSPLINE' : (RSPLINE, self._add_rigid_element_object),
            'RSSCON' : (RSSCON, self._add_rigid_element_object),


            ## there is no MAT6 or MAT7
            'MAT1' : (MAT1, self._add_structural_material_object),
            'MAT2' : (MAT2, self._add_structural_material_object),
            'MAT3' : (MAT3, self._add_structural_material_object),
            'MAT8' : (MAT8, self._add_structural_material_object),
            'MAT9' : (MAT9, self._add_structural_material_object),
            'MAT10' : (MAT10, self._add_structural_material_object),
            'MAT11' : (MAT11, self._add_structural_material_object),
            'MAT3D' : (MAT3D, self._add_structural_material_object),
            'EQUIV' : (EQUIV, self._add_structural_material_object),
            'MATG' : (MATG, self._add_structural_material_object),

            'MATHE' : (MATHE, self._add_hyperelastic_material_object),
            'MATHP' : (MATHP, self._add_hyperelastic_material_object),
            'MAT4' : (MAT4, self._add_thermal_material_object),
            'MAT5' : (MAT5, self._add_thermal_material_object),

            'MATS1' : (MATS1, self._add_material_dependence_object),
            #'MATS3' : (MATS3, self._add_material_dependence_object),
            #'MATS8' : (MATS8, self._add_material_dependence_object),
            'MATT1' : (MATT1, self._add_material_dependence_object),
            'MATT2' : (MATT2, self._add_material_dependence_object),
            'MATT3' : (MATT3, self._add_material_dependence_object),
            'MATT4' : (MATT4, self._add_material_dependence_object),
            'MATT5' : (MATT5, self._add_material_dependence_object),
            'MATT8' : (MATT8, self._add_material_dependence_object),
            'MATT9' : (MATT9, self._add_material_dependence_object),
            'NXSTRAT' : (NXSTRAT, self._add_nxstrat_object),

            # hasnt been verified, links up to MAT1, MAT2, MAT9 w/ same MID
            'CREEP' : (CREEP, self._add_creep_material_object),

            'NSMADD' : (NSMADD, self._add_nsmadd_object),
            'NSM1' : (NSM1, self._add_nsm_object),
            'NSML1' : (NSML1, self._add_nsm_object),

            'CONM1' : (CONM1, self._add_mass_object),
            'CONM2' : (CONM2, self._add_mass_object),
            'CMASS1' : (CMASS1, self._add_mass_object),
            'CMASS2' : (CMASS2, self._add_mass_object),
            'CMASS3' : (CMASS3, self._add_mass_object),
            # CMASS4 - added later because documentation is wrong

            'MPC' : (MPC, self._add_constraint_mpc_object),
            'MPCADD' : (MPCADD, self._add_constraint_mpcadd_object),

            'SPC' : (SPC, self._add_constraint_spc_object),
            'SPC1' : (SPC1, self._add_constraint_spc_object),
            'SPCOFF' : (SPCOFF, self._add_constraint_spcoff_object),
            'SPCOFF1' : (SPCOFF1, self._add_constraint_spcoff_object),
            'SPCAX' : (SPCAX, self._add_constraint_spc_object),
            'SPCADD' : (SPCADD, self._add_constraint_spcadd_object),
            'GMSPC' : (GMSPC, self._add_constraint_spc_object),

            'SESUP' : (SESUP, self._add_sesuport_object), # pseudo-constraint
            'SUPORT' : (SUPORT, self._add_suport_object), # pseudo-constraint
            'SUPORT1' : (SUPORT1, self._add_suport1_object),  # pseudo-constraint

            'FORCE' : (FORCE, self._add_load_object),
            'FORCE1' : (FORCE1, self._add_load_object),
            'FORCE2' : (FORCE2, self._add_load_object),
            'MOMENT' : (MOMENT, self._add_load_object),
            'MOMENT1' : (MOMENT1, self._add_load_object),
            'MOMENT2' : (MOMENT2, self._add_load_object),

            'LSEQ' : (LSEQ, self._add_lseq_object),
            'LOAD' : (LOAD, self._add_load_combination_object),
            'CLOAD' : (CLOAD, self._add_load_combination_object),
            'LOADCYN' : (LOADCYN, self._add_load_object),
            'LOADCYH' : (LOADCYH, self._add_load_object),

            # basic static loads
            'GRAV' : (GRAV, self._add_load_object),
            'ACCEL' : (ACCEL, self._add_load_object),
            'ACCEL1' : (ACCEL1, self._add_load_object),
            'PLOAD' : (PLOAD, self._add_load_object),
            'PLOAD1' : (PLOAD1, self._add_load_object),
            'PLOAD2' : (PLOAD2, self._add_load_object),
            'PLOAD4' : (PLOAD4, self._add_load_object),
            'RFORCE' : (RFORCE, self._add_load_object),
            'RFORCE1' : (RFORCE1, self._add_load_object),
            'SLOAD' : (SLOAD, self._add_load_object),
            'GMLOAD' : (GMLOAD, self._add_load_object),
            'SPCD' : (SPCD, self._add_load_object),  # enforced displacement
            'QVOL' : (QVOL, self._add_load_object),  # thermal

            # axisymmetric loads
            'FORCEAX' : (FORCEAX, self._add_load_object),
            'PLOADX1' : (PLOADX1, self._add_load_object),
            'PRESAX' : (PRESAX, self._add_load_object),  # axisymmetric

            'DLOAD' : (DLOAD, self._add_dload_object),

            'ACSRCE' : (ACSRCE, self._add_dload_entry),
            'TLOAD1' : (TLOAD1, self._add_dload_entry),
            'TLOAD2' : (TLOAD2, self._add_dload_entry),
            'RLOAD1' : (RLOAD1, self._add_dload_entry),
            'RLOAD2' : (RLOAD2, self._add_dload_entry),
            'RANDPS' : (RANDPS, self._add_dload_entry), # random
            'RANDT1' : (RANDT1, self._add_dload_entry), # random
            'QVECT' : (QVECT, self._add_dload_entry),

            'FREQ' : (FREQ, self._add_freq_object),
            'FREQ1' : (FREQ1, self._add_freq_object),
            'FREQ2' : (FREQ2, self._add_freq_object),
            'FREQ3' : (FREQ3, self._add_freq_object),
            'FREQ4' : (FREQ4, self._add_freq_object),
            'FREQ5' : (FREQ5, self._add_freq_object),

            'DOPTPRM' : (DOPTPRM, self._add_doptprm_object),
            'DESVAR' : (DESVAR, self._add_desvar_object),
            'TOPVAR' : (TOPVAR, self._add_topvar_object),
            # BCTSET

            'TEMPRB' : (TEMPRB, self._add_thermal_load_object),
            'TEMP' : (TEMP, self._add_thermal_load_object),
            'TEMPB3' : (TEMPB3, self._add_thermal_load_object),
            'QBDY1' : (QBDY1, self._add_thermal_load_object),
            'QBDY2' : (QBDY2, self._add_thermal_load_object),
            'QBDY3' : (QBDY3, self._add_thermal_load_object),
            'QHBDY' : (QHBDY, self._add_thermal_load_object),
            'PHBDY' : (PHBDY, self._add_phbdy_object),

            'CHBDYE' : (CHBDYE, self._add_thermal_element_object),
            'CHBDYG' : (CHBDYG, self._add_thermal_element_object),
            'CHBDYP' : (CHBDYP, self._add_thermal_element_object),
            'PCONV' : (PCONV, self._add_convection_property_object),
            'PCONVM' : (PCONVM, self._add_convection_property_object),

            'VIEW' : (VIEW, self._add_view_object),
            'VIEW3D' : (VIEW3D, self._add_view3d_object),

            # aero
            'AECOMP' : (AECOMP, self._add_aecomp_object),
            'AECOMPL' : (AECOMPL, self._add_aecomp_object),
            'AEFACT' : (AEFACT, self._add_aefact_object),
            'AELINK' : (AELINK, self._add_aelink_object),
            'AELIST' : (AELIST, self._add_aelist_object),
            'AEPARM' : (AEPARM, self._add_aeparm_object),
            'AESTAT' : (AESTAT, self._add_aestat_object),
            'AESURF' : (AESURF, self._add_aesurf_object),
            'AESURFS' : (AESURFS, self._add_aesurfs_object),

            'CAERO1' : (CAERO1, self._add_caero_object),
            'CAERO2' : (CAERO2, self._add_caero_object),
            'CAERO3' : (CAERO3, self._add_caero_object),
            'CAERO4' : (CAERO4, self._add_caero_object),
            'CAERO5' : (CAERO5, self._add_caero_object),

            'PAERO1' : (PAERO1, self._add_paero_object),
            'PAERO2' : (PAERO2, self._add_paero_object),
            'PAERO3' : (PAERO3, self._add_paero_object),
            'PAERO4' : (PAERO4, self._add_paero_object),
            'PAERO5' : (PAERO5, self._add_paero_object),

            'SPLINE1' : (SPLINE1, self._add_spline_object),
            'SPLINE2' : (SPLINE2, self._add_spline_object),
            'SPLINE3' : (SPLINE3, self._add_spline_object),
            'SPLINE4' : (SPLINE4, self._add_spline_object),
            'SPLINE5' : (SPLINE5, self._add_spline_object),

            # SOL 144
            'AEROS' : (AEROS, self._add_aeros_object),
            'TRIM' : (TRIM, self._add_trim_object),
            'TRIM2' : (TRIM2, self._add_trim_object),
            'DIVERG' : (DIVERG, self._add_diverg_object),

            # SOL 145
            'AERO' : (AERO, self._add_aero_object),
            'FLUTTER' : (FLUTTER, self._add_flutter_object),
            'FLFACT' : (FLFACT, self._add_flfact_object),
            'MKAERO1' : (MKAERO1, self._add_mkaero_object),
            'MKAERO2' : (MKAERO2, self._add_mkaero_object),

            'GUST' : (GUST, self._add_gust_object),
            'CSSCHD' : (CSSCHD, self._add_csschd_object),
            'MONPNT1' : (MONPNT1, self._add_monpnt_object),
            'MONPNT2' : (MONPNT2, self._add_monpnt_object),
            'MONPNT3' : (MONPNT3, self._add_monpnt_object),
            'MONDSP1' : (MONDSP1, self._add_monpnt_object),

            'NLPARM' : (NLPARM, self._add_nlparm_object),
            'NLPCI' : (NLPCI, self._add_nlpci_object),
            'TSTEP' : (TSTEP, self._add_tstep_object),
            'TSTEP1' : (TSTEP1, self._add_tstepnl_object),
            'TSTEPNL' : (TSTEPNL, self._add_tstepnl_object),

            'TF' : (TF, self._add_tf_object),
            'TIC' : (TIC, self._add_tic_object),

            'DCONADD' : (DCONADD, self._add_dconstr_object),
            'DCONSTR' : (DCONSTR, self._add_dconstr_object),
            'DDVAL' : (DDVAL, self._add_ddval_object),
            'DLINK' : (DLINK, self._add_dlink_object),
            'DSCREEN' : (DSCREEN, self._add_dscreen_object),

            'DTABLE' : (DTABLE, self._add_dtable_object),
            'DRESP1' : (DRESP1, self._add_dresp_object), # dresps
            'DRESP2' : (DRESP2, self._add_dresp_object),
            'DRESP3' : (DRESP3, self._add_dresp_object),
            'DVCREL1' : (DVCREL1, self._add_dvcrel_object), # dvcrels
            'DVCREL2' : (DVCREL2, self._add_dvcrel_object),
            'DVPREL1' : (DVPREL1, self._add_dvprel_object), # dvprels
            'DVPREL2' : (DVPREL2, self._add_dvprel_object),
            'DVMREL1' : (DVMREL1, self._add_dvmrel_object), # ddvmrels
            'DVMREL2' : (DVMREL2, self._add_dvmrel_object),
            'DVGRID' : (DVGRID, self._add_dvgrid_object), # dvgrids

            # tables
            'TABLES1' : (TABLES1, self._add_table_object),
            'TABLEST' : (TABLEST, self._add_table_object),
            'TABLEHT' : (TABLEHT, self._add_table_object),
            'TABLEH1' : (TABLEH1, self._add_table_object),

            # dynamic tables
            'TABLED1' : (TABLED1, self._add_tabled_object),
            'TABLED2' : (TABLED2, self._add_tabled_object),
            'TABLED3' : (TABLED3, self._add_tabled_object),
            'TABLED4' : (TABLED4, self._add_tabled_object),

            # material tables
            'TABLEM1' : (TABLEM1, self._add_tablem_object),
            'TABLEM2' : (TABLEM2, self._add_tablem_object),
            'TABLEM3' : (TABLEM3, self._add_tablem_object),
            'TABLEM4' : (TABLEM4, self._add_tablem_object),

            # other tables
            'TABDMP1' : (TABDMP1, self._add_table_sdamping_object),
            'TABRND1' : (TABRND1, self._add_random_table_object),
            'TABRNDG' : (TABRNDG, self._add_random_table_object),

            'EIGB' : (EIGB, self._add_method_object),
            'EIGR' : (EIGR, self._add_method_object),
            'EIGRL' : (EIGRL, self._add_method_object),
            'EIGC' : (EIGC, self._add_cmethod_object),
            'EIGP' : (EIGP, self._add_cmethod_object),

            'BCRPARA' : (BCRPARA, self._add_bcrpara_object),
            'BCTADD' : (BCTADD, self._add_bctadd_object),
            'BCTPARA' : (BCTPARA, self._add_bctpara_object),
            'BSURF' : (BSURF, self._add_bsurf_object),
            'BSURFS' : (BSURFS, self._add_bsurfs_object),
            # 'BGADD', 'BGSET', 'BOUTPUT', 'BOLT', 'BOLTFOR', 'BOLTFRC',
            'BGADD': (Crash, None),
            'BGSET': (Crash, None),
            'BOUTPUT': (Crash, None),
            'BOLT': (Crash, None),
            'BOLTFOR': (Crash, None),
            'BOLTFRC': (Crash, None),

            'RADCAV' : (RADCAV, self._add_radcav_object), #
            #'RADLST' : (RADLST, self._add_radcav_object), # TestOP2.test_bdf_op2_thermal_02
            #'RADMTX' : (RADMTX, self._add_radmtx_object), # TestOP2.test_bdf_op2_thermal_02
            #'RADMT' : (Crash, None),

            'ASET' : (ASET, self._add_aset_object),
            'ASET1' : (ASET1, self._add_aset_object),

            'BSET' : (BSET, self._add_bset_object),
            'BSET1' : (BSET1, self._add_bset_object),

            'CSET' : (CSET, self._add_cset_object),
            'CSET1' : (CSET1, self._add_cset_object),

            'QSET' : (QSET, self._add_qset_object),
            'QSET1' : (QSET1, self._add_qset_object),

            'USET' : (USET, self._add_uset_object),
            'USET1' : (USET1, self._add_uset_object),

            'OMIT' : (OMIT, self._add_omit_object),
            'OMIT1' : (OMIT1, self._add_omit_object),

            'SET1' : (SET1, self._add_set_object),
            'SET3' : (SET3, self._add_set_object),

            # radset
            'RADSET' : (RADSET, self._add_radset_object),

            'SESET' : (SESET, self._add_seset_object),

            'SEBSET' : (SEBSET, self._add_sebset_object),
            'SEBSET1' : (SEBSET1, self._add_sebset_object),

            'SECSET' : (SECSET, self._add_secset_object),
            'SECSET1' : (SECSET1, self._add_secset_object),

            'SEQSET' : (SEQSET, self._add_seqset_object),
            'SEQSET1' : (SEQSET1, self._add_seqset_object),

            #'SESUP' : (SESUP, self._add_sesup_object),  # pseudo-constraint

            #'SEUSET' : (SEUSET, self._add_seuset_object),
            #'SEUSET1' : (SEUSET1, self._add_seuset_object),

            # BCTSET
            'ROTORG' : (ROTORG, self._add_rotor_object),
            'ROTORD' : (ROTORD, self._add_rotor_object),

            'DAREA' : (DAREA, self._add_darea_object),
            'DPHASE' : (DPHASE, self._add_dphase_object),
            'DELAY' : (DELAY, self._add_delay_object),

            'CYJOIN' : (CYJOIN, self._add_cyjoin_object),
        }

        self._card_parser_prepare = {
            'PLOTEL': self._prepare_plotel,
            'CBAR' : self._prepare_cbar,
            'CBEAM' : self._prepare_cbeam,
            'CTETRA' : self._prepare_ctetra,
            'CPENTA' : self._prepare_cpenta,
            'CHEXA' : self._prepare_chexa,
            'CPYRAM' : self._prepare_cpyram,

            'CORD1R' : self._prepare_cord1r,
            'CORD1C' : self._prepare_cord1c,
            'CORD1S' : self._prepare_cord1s,
            #'CORD3G' : self._prepare_CORD3G,

            'PMASS' : self._prepare_pmass,
            'CMASS4' : self._prepare_cmass4,
            'CDAMP4' : self._prepare_cdamp4,
            'DEFORM' : self._prepare_deform,  # enforced displacement

            'DTI' : self._prepare_dti,
            'DMIG' : self._prepare_dmig,
            'DMIAX' : self._prepare_dmiax,
            'DMI' : self._prepare_dmi,
            'DMIJ' : self._prepare_dmij,
            'DMIK' : self._prepare_dmik,
            'DMIJI' : self._prepare_dmiji,
            'RINGFL' : self._prepare_ringfl,

            'DEQATN' : self._prepare_dequatn,

            'NSM' : self._prepare_nsm,
            'NSML' : self._prepare_nsml,
            'PVISC' : self._prepare_pvisc,
            'PELAS' : self._prepare_pelas,
            'PDAMP' : self._prepare_pdamp,

            'TEMPAX' : self._prepare_tempax,
            'TEMPD' : self._prepare_tempd,
            'TEMPBC' : self._prepare_tempbc,
            'CONVM' : self._prepare_convm,
            'CONV' : self._prepare_conv,
            'RADM' : self._prepare_radm,
            'RADBC' : self._prepare_radbc,
            # GRDSET-will be last card to update from _card_parser_prepare
            'GRDSET' : self._prepare_grdset,

            'BCTSET' : self._prepare_bctset,
            'ACMODL' : self._prepare_acmodl,
        }

        new_reject_method = False
        if new_reject_method:
            self.cards_to_read.remove('CONM2')
            # bwb_saero: 6.37-6.57 before using

            all_cards = set(self._card_parser.keys() + self._card_parser_prepare.keys())
            reject_cards = all_cards.difference(self.cards_to_read)
            #print('reject_cards =', reject_cards)
            for card_name in reject_cards:
                if card_name in self._card_parser:
                    del self._card_parser[card_name]
                self._card_parser_prepare[card_name] = self._reject_card_obj2

    def _reject_card_obj2(self, unused_card_name, card_obj):
        """rejects a card object"""
        self.reject_cards.append(card_obj)

    def reject_card_lines(self, card_name: str, card_lines: List[str],
                          show_log: bool=True, comment: str='') -> None:
        """rejects a card"""
        if card_name.isdigit():
            # TODO: this should technically work (I think), but it's a problem
            #       for the code
            #
            # prevents:
            # spc1,100,456,10013832,10013833,10013830,10013831,10013836,10013837,
            # 10013834,10013835,10013838,10013839,10014508,10008937,10008936,10008935,
            msg = 'card_name=%r was misparsed...\ncard_lines=%s' % (
                card_name, card_lines)
            raise RuntimeError(msg)

        if card_name not in self.card_count:
            _check_for_spaces(card_name, card_lines, comment, self.log)
            if card_name == '\ufeff':
                self.log.warning('    rejecting card_name = %r' % card_name)
                self.log.warning('    comment:\n')
                print(comment)
                self.log.warning('    lines:\n')
                for line in card_lines:
                    print(line)
            elif show_log:
                self.log.info('    rejecting card_name = %r' % card_name)
            assert isinstance(show_log, bool), show_log
        self.increase_card_count(card_name)
        self.reject_lines.append([_format_comment(comment)] + card_lines)

    def _write_reject_message(self, card_name, unused_card_obj, comment=''):
        """common method to not write duplicate reject card names"""
        if card_name not in self.card_count:
            #if ' ' in card_name:
                #_check_for_spaces(card_name, card_lines, comment, self.log)
            self.log.info('    rejecting card_name = %s' % card_name)

    def _prepare_plotel(self, unused_card: List[str], card_obj: BDFCard, comment='') -> None:
        """adds a PLOTEL"""
        #['PLOTEL', '3101', '3101', '3102', None, '3102', '3102', '3103']
        plotels = [PLOTEL.add_card(card_obj, 0, comment=comment)]
        if card_obj.field(5):  # eid
            plotels.append(PLOTEL.add_card(card_obj, 1, comment=''))
        for plotel in plotels:
            self._add_plotel_object(plotel)
        return plotels

    def _prepare_cbar(self, unused_card: List[str], card_obj: BDFCard, comment='') -> None:
        """adds a CBAR"""
        elem = CBAR.add_card(card_obj, baror=self.baror, comment=comment)
        self._add_element_object(elem)
        return elem

    def _prepare_cbeam(self, unused_card: List[str], card_obj: BDFCard, comment='') -> None:
        """adds a CBEAM"""
        elem = CBEAM.add_card(card_obj, beamor=self.beamor, comment=comment)
        self._add_element_object(elem)
        return elem

    def _prepare_ctetra(self, unused_card: List[str], card_obj: BDFCard, comment='') -> None:
        """adds a CTETRA4/CTETRA10"""
        if len(card_obj) == 7:
            elem = CTETRA4.add_card(card_obj, comment=comment)
        else:
            elem = CTETRA10.add_card(card_obj, comment=comment)
        self._add_element_object(elem)
        return elem

    def _prepare_cpyram(self, unused_card: List[str], card_obj: BDFCard, comment='') -> None:
        """adds a CPYRAM5/CPYRAM13"""
        if len(card_obj) == 8:
            elem = CPYRAM5.add_card(card_obj, comment=comment)
        else:
            elem = CPYRAM13.add_card(card_obj, comment=comment)
        self._add_element_object(elem)
        return elem

    def _prepare_cpenta(self, unused_card: List[str], card_obj: BDFCard, comment='') -> None:
        """adds a CPENTA6/CPENTA15"""
        if len(card_obj) == 9:
            elem = CPENTA6.add_card(card_obj, comment=comment)
        else:
            elem = CPENTA15.add_card(card_obj, comment=comment)
        self._add_element_object(elem)
        return elem

    def _prepare_chexa(self, unused_card: List[str], card_obj: BDFCard, comment='') -> None:
        """adds a CHEXA8/CHEXA20"""
        if len(card_obj) == 11:
            elem = CHEXA8.add_card(card_obj, comment=comment)
        else:
            elem = CHEXA20.add_card(card_obj, comment=comment)
        self._add_element_object(elem)
        return elem

    def _prepare_bctset(self, unused_card: List[str], card_obj: BDFCard, comment='') -> None:
        """adds a BCTSET"""
        bctset = BCTSET.add_card(card_obj, comment=comment, sol=self.sol)
        self._add_bctset_object(bctset)
        return bctset

    def _prepare_grdset(self, unused_card: List[str], card_obj: BDFCard, comment='') -> None:
        """adds a GRDSET"""
        self.grdset = GRDSET.add_card(card_obj, comment=comment)
        return self.grdset

    def _prepare_cdamp4(self, unused_card: List[str], card_obj: BDFCard, comment='') -> None:
        """adds a CDAMP4"""
        dampers = [CDAMP4.add_card(card_obj, comment=comment)]
        if card_obj.field(5):
            dampers.append(CDAMP4.add_card(card_obj, 1, comment=''))
        for damper in dampers:
            self._add_damper_object(damper)
        return dampers

    def _prepare_deform(self, unused_card: List[str], card_obj: BDFCard, comment='') -> None:
        """adds a DEFORM"""
        loads = [DEFORM.add_card(card_obj, comment=comment)]
        if card_obj.field(4):
            loads.append(DEFORM.add_card(card_obj, 1, comment=comment))
        if card_obj.field(6):
            loads.append(DEFORM.add_card(card_obj, 2, comment=comment))
        for loadi in loads:
            self._add_load_object(loadi)
        return loads

    def _prepare_tempbc(self, unused_card: List[str], card_obj: BDFCard, comment='') -> None:
        """adds a TEMPBC"""
        boundary_condition = TEMPBC.add_card(card_obj, comment=comment)
        self._add_thermal_bc_object(boundary_condition, boundary_condition.eid)
        return boundary_condition

    def _prepare_convm(self, unused_card: List[str], card_obj: BDFCard, comment='') -> None:
        """adds a CONVM"""
        boundary_condition = CONVM.add_card(card_obj, comment=comment)
        self._add_thermal_bc_object(boundary_condition, boundary_condition.eid)
        return boundary_condition

    def _prepare_conv(self, unused_card: List[str], card_obj: BDFCard, comment='') -> None:
        """adds a CONV"""
        boundary_condition = CONV.add_card(card_obj, comment=comment)
        self._add_thermal_bc_object(boundary_condition, boundary_condition.eid)
        return boundary_condition

    def _prepare_radm(self, unused_card: List[str], card_obj: BDFCard, comment='') -> None:
        """adds a RADM"""
        boundary_condition = RADM.add_card(card_obj, comment=comment)
        self._add_thermal_bc_object(boundary_condition, boundary_condition.radmid)
        return boundary_condition

    def _prepare_radbc(self, unused_card: List[str], card_obj: BDFCard, comment='') -> None:
        """adds a RADBC"""
        boundary_condition = RADBC.add_card(card_obj, comment=comment)
        self._add_thermal_bc_object(boundary_condition, boundary_condition.nodamb)
        return boundary_condition

    def _prepare_tempd(self, unused_card: List[str], card_obj: BDFCard, comment='') -> None:
        """adds a TEMPD"""
        tempds = [TEMPD.add_card(card_obj, 0, comment=comment)]
        if card_obj.field(3):
            tempds.append(TEMPD.add_card(card_obj, 1, comment=''))
            if card_obj.field(5):
                tempds.append(TEMPD.add_card(card_obj, 2, comment=''))
                if card_obj.field(7):
                    tempds.append(TEMPD.add_card(card_obj, 3, comment=''))
        for tempd in tempds:
            self._add_tempd_object(tempd)
        return tempds

    def _prepare_tempax(self, unused_card: List[str], card_obj: BDFCard, comment='') -> None:
        """adds a TEMPAX"""
        tempaxs = [TEMPAX.add_card(card_obj, 0, comment=comment)]
        if card_obj.field(5):
            tempaxs.append(TEMPAX.add_card(card_obj, 1, comment=''))
        for tempax in tempaxs:
            self._add_load_object(tempax)
        return tempaxs

    def _prepare_dequatn(self, unused_card: List[str], card_obj: BDFCard, comment='') -> None:
        """adds a DEQATN"""
        deqatn = DEQATN.add_card(card_obj, comment=comment)
        self._add_deqatn_object(deqatn)
        return deqatn

    def _prepare_dti(self, unused_card_name, card_obj, comment=''):
        """adds a DTI"""
        #name = string(card_obj, 1, 'name')
        dti = DTI.add_card(card_obj, comment=comment)
        self._add_dti_object(dti)
        return dti

    def _prepare_dmig(self, unused_card: List[str], card_obj: BDFCard, comment='') -> None:
        """adds a DMIG"""
        name = string(card_obj, 1, 'name')
        field2 = integer_or_string(card_obj, 2, 'flag')

        if name == 'UACCEL':  # special DMIG card
            if field2 == 0:
                dmig = DMIG_UACCEL.add_card(card_obj, comment=comment)
                self._add_dmig_object(dmig)
            else:
                dmig = -1
                self._dmig_temp[name].append((card_obj, comment))
        else:
            field2 = integer_or_string(card_obj, 2, 'flag')
            if field2 == 0:
                dmig = DMIG.add_card(card_obj, comment=comment)
                self._add_dmig_object(dmig)
            else:
                dmig = -1
                self._dmig_temp[name].append((card_obj, comment))
        return dmig

    def _prepare_dmix(self, class_obj, add_method, card_obj, comment=''):
        """adds a DMI, DMIJ, DMIJI, or DMIK"""
        field2 = integer(card_obj, 2, 'flag')
        if field2 == 0:
            dmix = class_obj.add_card(card_obj, comment=comment)
            add_method(dmix)
        else:
            dmix = -1
            name = string(card_obj, 1, 'name')
            self._dmig_temp[name].append((card_obj, comment))
        return dmix

    def _prepare_dmiax(self, unused_card: List[str], card_obj: BDFCard, comment='') -> None:
        """adds a DMIAX"""
        return self._prepare_dmix(DMIAX, self._add_dmiax_object, card_obj, comment=comment)

    def _prepare_dmi(self, unused_card: List[str], card_obj: BDFCard, comment='') -> None:
        """adds a DMI"""
        return self._prepare_dmix(DMI, self._add_dmi_object, card_obj, comment=comment)

    def _prepare_dmij(self, unused_card: List[str], card_obj: BDFCard, comment='') -> None:
        """adds a DMIJ"""
        return self._prepare_dmix(DMIJ, self._add_dmij_object, card_obj, comment=comment)

    def _prepare_dmik(self, unused_card: List[str], card_obj: BDFCard, comment='') -> None:
        """adds a DMIK"""
        return self._prepare_dmix(DMIK, self._add_dmik_object, card_obj, comment=comment)

    def _prepare_dmiji(self, unused_card: List[str], card_obj: BDFCard, comment='') -> None:
        """adds a DMIJI"""
        return self._prepare_dmix(DMIJI, self._add_dmiji_object, card_obj, comment=comment)

    def _prepare_cmass4(self, unused_card: List[str], card_obj: BDFCard, comment='') -> None:
        """adds a CMASS4"""
        elements = [CMASS4.add_card(card_obj, icard=0, comment=comment)]
        if card_obj.field(5):
            elements.append(CMASS4.add_card(card_obj, icard=1, comment=comment))
        for elem in elements:
            self._add_mass_object(elem)
        return elements

    def _prepare_pelas(self, unused_card: List[str], card_obj: BDFCard, comment='') -> None:
        """adds a PELAS"""
        properties = [PELAS.add_card(card_obj, icard=0, comment=comment)]
        if card_obj.field(5):
            properties.append(PELAS.add_card(card_obj, icard=1, comment=comment))
        for prop in properties:
            self._add_property_object(prop)
        return properties

    def _prepare_nsm(self, unused_card: List[str], card_obj: BDFCard, comment='') -> None:
        """adds an NSM"""
        nfields = len(card_obj)
        ncards = (nfields - 3) // 2
        nextra = (nfields - 3) % 2
        assert nextra == 0, 'NSM error; nfields=%s must have an odd number of fields\ncard=%s' % (
            nfields, card_obj)

        nsms = []
        for icard in range(ncards):
            nsms.append(NSM.add_card(card_obj, icard, comment=comment))
        for nsm in nsms:
            self._add_nsm_object(nsm)
        return nsms

    def _prepare_nsml(self, unused_card: List[str], card_obj: BDFCard, comment='') -> None:
        """adds an NSML"""
        nfields = len(card_obj)
        ncards = (nfields - 3) // 2
        nextra = (nfields - 3) % 2
        assert nextra == 0, 'NSML error; nfields=%s must have an odd number of fields\ncard=%s' % (
            nfields, card_obj)
        nsms = []
        for icard in range(ncards):
            nsms.append(NSML.add_card(card_obj, icard, comment=comment))
        for nsm in nsms:
            self._add_nsm_object(nsm)
        return nsms

    def _prepare_pvisc(self, unused_card: List[str], card_obj: BDFCard, comment='') -> None:
        """adds a PVISC"""
        properties = [PVISC.add_card(card_obj, icard=0, comment=comment)]
        if card_obj.field(5):
            properties.append(PVISC.add_card(card_obj, icard=1, comment=comment))
        for prop in properties:
            self._add_property_object(prop)
        return properties

    def _prepare_ringfl(self, unused_card: List[str], card_obj: BDFCard, comment='') -> None:
        """adds a RINGFL"""
        rings = [RINGFL.add_card(card_obj, icard=0, comment=comment)]
        if card_obj.field(3):
            rings.append(RINGFL.add_card(card_obj, icard=1, comment=comment))
        if card_obj.field(5):
            rings.append(RINGFL.add_card(card_obj, icard=2, comment=comment))
        for ring in rings:
            self._add_ringfl_object(ring)
        return rings

    def _prepare_pdamp(self, unused_card: List[str], card_obj: BDFCard, comment='') -> None:
        """adds a PDAMP"""
        properties = [PDAMP.add_card(card_obj, icard=0, comment=comment)]
        if card_obj.field(3):
            properties.append(PDAMP.add_card(card_obj, icard=1, comment=comment))
        if card_obj.field(5):
            properties.append(PDAMP.add_card(card_obj, icard=2, comment=comment))
        if card_obj.field(7):
            properties.append(PDAMP.add_card(card_obj, icard=3, comment=comment))
        for prop in properties:
            self._add_property_object(prop)
        return properties

    def _prepare_pmass(self, unused_card: List[str], card_obj: BDFCard, comment='') -> None:
        """adds a PMASS"""
        properties = [PMASS.add_card(card_obj, icard=0, comment=comment)]
        for (i, j) in enumerate([3, 5, 7]):
            if card_obj.field(j):
                properties.append(PMASS.add_card(card_obj, icard=i+1, comment=comment))
        for prop in properties:
            self._add_property_mass_object(prop)
        return properties

    def _prepare_cord1r(self, unused_card: List[str], card_obj: BDFCard, comment='') -> None:
        """adds a CORD1R"""
        coords = [CORD1R.add_card(card_obj, comment=comment)]
        if card_obj.field(5):
            coords.append(CORD1R.add_card(card_obj, icard=1, comment=comment))
        for coord in coords:
            self._add_coord_object(coord)
        return coords

    def _prepare_cord1c(self, unused_card: List[str], card_obj: BDFCard, comment='') -> None:
        """adds a CORD1C"""
        coords = [CORD1C.add_card(card_obj, comment=comment)]
        if card_obj.field(5):
            coords.append(CORD1C.add_card(card_obj, icard=1, comment=comment))
        for coord in coords:
            self._add_coord_object(coord)
        return coords

    def _prepare_cord1s(self, unused_card: List[str], card_obj: BDFCard, comment='') -> None:
        """adds a CORD1S"""
        coords = [CORD1S.add_card(card_obj, comment=comment)]
        if card_obj.field(5):
            coords.append(CORD1S.add_card(card_obj, icard=1, comment=comment))
        for coord in coords:
            self._add_coord_object(coord)
        return coords

    def _prepare_acmodl(self, unused_card: List[str], card_obj: BDFCard, comment='') -> None:
        acmodl = ACMODL.add_card(card_obj, self._nastran_format, comment=comment)
        self._add_acmodl_object(acmodl)

    def add_card_ifile(self, ifile: int, card_lines: List[str], card_name: str,
                       comment: str='', is_list: bool=True, has_none: bool=True) -> Any:
        """Same as ``add_card`` except it has an ifile parameter"""
        assert isinstance(ifile, (int, np.int32)), 'ifile=%s type=%s' % (ifile, type(ifile))
        card_name = card_name.upper()
        card_obj, unused_card = self.create_card_object(
            card_lines, card_name,
            is_list=is_list, has_none=has_none)
        self._add_card_helper_ifile(ifile, card_obj, card_name, card_name, comment)
        return card_obj

    def add_card(self, card_lines: List[str], card_name: str,
                 comment: str='', ifile=None, is_list: bool=True, has_none: bool=True) -> Any:
        """
        Adds a card object to the BDF object.

        Parameters
        ----------
        card_lines: list[str]
            the list of the card fields
        card_name : str
            the card_name -> 'GRID'
        comment : str
            an optional the comment for the card
        is_list : bool, optional
            False : input is a list of card fields -> ['GRID', 1, None, 3.0, 4.0, 5.0]
            True :  input is list of card_lines -> ['GRID, 1,, 3.0, 4.0, 5.0']
        has_none : bool; default=True
            can there be trailing Nones in the card data (e.g. ['GRID', 1, 2, 3.0, 4.0, 5.0, None])
            can there be trailing Nones in the card data (e.g. ['GRID, 1, 2, 3.0, 4.0, 5.0, '])

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

           # here is_list=False because it's been parsed
           # Note the None at the end of the 1st line, which is there
           #      because the CONM2 card has a blank field.
           #      It must be there.
           # We also set i32 on the 2nd line, so it will default to 0.0
           >>> card = [
                   'CONM2', eid, nid, cid, mass, x1, x2, x3, None,
                            i11, i21, i22, i31, None, i33,
               ]
           >>> model.add_card(card_lines, 'CONM2', comment, is_list=False)

           # here's an alternate approach for the CONM2
           # we use Nastran's CSV format
           # There are many blank fields, but it's parsed exactly like a
           # standard CONM2.
           >>> card = [
                   'CONM2,1,2,3,10.0',
                   ',1.0,,5.0'
               ]
           >>> model.add_card(card_lines, 'CONM2', comment, is_list=True)

        .. note:: this is a very useful method for interfacing with the code
        .. note:: the card_object is not a card-type object...so not a GRID
                  card or CQUAD4 object.  It's a BDFCard Object.  However,
                  you know the type (assuming a GRID), so just call the
                  *mesh.Node(nid)* to get the Node object that was just
                  created.

        """
        card_name = card_name.upper()
        card_obj, unused_card = self.create_card_object(
            card_lines, card_name,
            is_list=is_list, has_none=has_none)
        self._add_card_helper(card_obj, card_name, card_name, ifile, comment)
        return card_obj

    def add_card_fields(self, card_lines, card_name, comment='', has_none=True):
        """
        Adds a card object to the BDF object.

        Parameters
        ----------
        card_lines: list[str]
            the list of the card fields
            input is a list of card fields -> ['GRID', 1, 2, 3.0, 4.0, 5.0]
        card_name : str
            the card_name -> 'GRID'
        comment : str
            an optional the comment for the card
        has_none : bool; default=True
            can there be trailing Nones in the card data (e.g. ['GRID', 1, 2, 3.0, 4.0, 5.0, None])

        Returns
        -------
        card_object : BDFCard()
            the card object representation of card

        """
        card_name = card_name.upper()
        card_obj, card = self.create_card_object(card_lines, card_name,
                                                 is_list=True, has_none=has_none)
        self._add_card_helper(card_obj, card, card_name, comment)

    def add_card_lines(self, card_lines, card_name, comment='', has_none=True):
        """
        Adds a card object to the BDF object.

        Parameters
        ----------
        card_lines: list[str]
            the list of the card fields
            input is list of card_lines -> ['GRID, 1, 2, 3.0, 4.0, 5.0']
        card_name : str
            the card_name -> 'GRID'
        comment : str; default=''
            an optional the comment for the card
        has_none : bool; default=True
            can there be trailing Nones in the card data (e.g. ['GRID, 1, 2, 3.0, 4.0, 5.0, '])

        """
        card_name = card_name.upper()
        card_obj, card = self.create_card_object(card_lines, card_name,
                                                 is_list=False, has_none=has_none)
        self._add_card_helper(card_obj, card, card_name, comment)

    def get_xyz_in_coord_no_xref(self, cid=0, fdtype='float64', sort_ids=True):
        """see get_xyz_in_coord"""
        npoints, nids, all_nodes = self._get_npoints_nids_allnids()

        xyz_cid0 = np.zeros((npoints, 3), dtype=fdtype)
        if cid == 0:
            for i, nid in enumerate(nids):
                node = self.nodes[nid]
                xyz = node.get_position_no_xref(self)
                xyz_cid0[i, :] = xyz
        else:
            for i, nid in enumerate(nids):
                node = self.nodes[nid]
                xyz = node.get_position_wrt_no_xref(self, cid)
                xyz_cid0[i, :] = xyz
        if sort_ids:
            isort = np.argsort(all_nodes)
            xyz_cid0 = xyz_cid0[isort, :]
        return all_nodes, xyz_cid0

    def _get_npoints_nids_allnids(self):
        """helper method for get_xyz_in_coord"""
        if self.is_bdf_vectorized:
            return self.nodes.nids

        nnodes = len(self.nodes)
        nspoints = 0
        nepoints = 0
        nids = list(self.node_ids)
        all_nodes = list(self.node_ids)
        if self.spoints:
            spoints = list(self.spoints)
            nspoints = len(spoints)
            all_nodes += spoints
        if self.epoints:
            epoints = list(self.epoints)
            nepoints = len(epoints)
            all_nodes += epoints
        #self.log.debug('all_nodes = %s' % all_nodes)

        npoints = nnodes + nspoints + nepoints
        if len(all_nodes) != npoints:
            msg = 'len(unique(all_nodes))=%s npoints=%s\n' % (len(np.unique(all_nodes)), npoints)
            msg += 'npoints = nnodes+nspoints+nepoints = %s + %s + %s\n' % (
                nnodes, nspoints, nepoints)
            msg += 'all_nodes=%s' % (all_nodes)
            raise RuntimeError(msg)
        if npoints == 0:
            msg = 'nnodes=%s nspoints=%s nepoints=%s' % (nnodes, nspoints, nepoints)
            raise ValueError(msg)
        return npoints, nids, all_nodes

    def get_xyz_in_coord(self, cid=0, fdtype='float64', sort_ids=True):
        """
        Gets the xyz points (including SPOINTS) in the desired coordinate frame

        Parameters
        ----------
        cid : int; default=0
            the desired coordinate system
        fdtype : str; default='float64'
            the data type of the xyz coordinates
        sort_ids : bool; default=True
            sort the ids

        Returns
        -------
        xyz : (n, 3) ndarray
            the xyz points in the cid coordinate frame

        """
        #return self.get_displacement_index_xyz_cp_cd(cid=cid, fdtype=dtype)[2]
        npoints, nids, all_nodes = self._get_npoints_nids_allnids()
        xyz_cid0 = np.zeros((npoints, 3), dtype=fdtype)
        if cid == 0:
            for i, nid in enumerate(nids):
                node = self.nodes[nid]
                xyz = node.get_position()
                xyz_cid0[i, :] = xyz
        else:
            for i, nid in enumerate(nids):
                node = self.nodes[nid]
                xyz = node.get_position_wrt(self, cid)
                xyz_cid0[i, :] = xyz
        if sort_ids:
            isort = np.argsort(all_nodes)
            xyz_cid0 = xyz_cid0[isort, :]
        return xyz_cid0

    def _add_card_helper_ifile(self, ifile, card_obj, card, card_name, comment=''):
        """See ``_add_card_helper``"""
        if card_name == 'ECHOON':
            self.echo = True
            return
        elif card_name == 'ECHOOFF':
            self.echo = False
            return

        if self.echo and not self.force_echo_off:
            _echo_card(card, card_obj)

        if card_name in self._card_parser:
            card_class, add_card_function = self._card_parser[card_name]
            try:
                class_instance = card_class.add_card(card_obj, comment=comment)
                class_instance.ifile = ifile
                add_card_function(class_instance)
            except TypeError:
                # this should never be turned on, but is useful for testing
                msg = 'problem adding %s' % card_obj
                print(msg)
                raise
            except (SyntaxError, AssertionError, KeyError, ValueError) as exception:
                self._iparse_errors += 1
                self.log.error(card_obj)
                var = traceback.format_exception_only(type(exception), exception)
                self._stored_parse_errors.append((card, var))
                if self._iparse_errors > self._nparse_errors:
                    self.pop_parse_errors()

        elif card_name in self._card_parser_prepare:
            add_card_function = self._card_parser_prepare[card_name]
            try:
                obj = add_card_function(card, card_obj, comment=comment)
            except (SyntaxError, AssertionError, KeyError, ValueError) as exception:
                #raise
                self._iparse_errors += 1
                self.log.error(card_obj)
                var = traceback.format_exception_only(type(exception), exception)
                self._stored_parse_errors.append((card, var))
                if self._iparse_errors > self._nparse_errors:
                    self.pop_parse_errors()

            if obj is None:
                print(add_card_function)
                print(card)
                print(card_obj)
                raise RuntimeError('_prepare_%s needs to implement obj' % card_name.lower())
            elif isinstance(obj, list):
                for obji in obj:
                    obji.ifile = ifile
            elif obj == -1:
                pass
            else:
                obj.ifile = ifile

        else:
            self.reject_cards.append(card_obj)

    def _add_card_helper(self, card_obj: BDFCard, card: List[str],
                         card_name: str, ifile: int,
                         comment: str='') -> None:
        """
        Adds a card object to the BDF object.

        Parameters
        ----------
        card_object : BDFCard()
            the card object representation of card
        card : List[str]
            the fields of the card object; used for rejection and special cards
        card_name : str
            the card_name -> 'GRID'
        ifile : int
            the file number
        comment : str
            an optional the comment for the card

        """
        if card_name == 'ECHOON':
            self.echo = True
            return
        elif card_name == 'ECHOOFF':
            self.echo = False
            return

        if self.echo and not self.force_echo_off:
            _echo_card(card, card_obj)

        if card_name in self._card_parser:
            card_class, add_card_function = self._card_parser[card_name]
            try:
                class_instance = card_class.add_card(card_obj, comment=comment)
                add_card_function(class_instance)
            except TypeError:
                # this should never be turned on, but is useful for testing
                msg = 'problem adding %s' % card_obj
                print(msg)
                raise
                #raise TypeError(msg)
            except (SyntaxError, AssertionError, KeyError, ValueError) as exception:
                print('problem adding %s' % card_obj)
                #msg = 'problem adding card from %s' % self.active_filenames[ifile]
                #raise
                # WARNING: Don't catch RuntimeErrors or a massive memory leak can occur
                #tpl/cc451.bdf
                #raise
                # NameErrors should be caught
                self._iparse_errors += 1
                #self.log.error(str(card_obj))
                var = traceback.format_exception_only(type(exception), exception)
                self._stored_parse_errors.append((card, var))
                if self._iparse_errors > self._nparse_errors:
                    self.pop_parse_errors()
                #raise
            #except AssertionError as exception:
                #self.log.error(str(card_obj))

        elif card_name in self._card_parser_prepare:
            add_card_function = self._card_parser_prepare[card_name]
            try:
                add_card_function(card, card_obj, comment=comment)
            except (SyntaxError, AssertionError, KeyError, ValueError) as exception:
                #raise
                # WARNING: Don't catch RuntimeErrors or a massive memory leak can occur
                #tpl/cc451.bdf
                #raise
                # NameErrors should be caught
                self._iparse_errors += 1
                self.log.error(str(card_obj))
                var = traceback.format_exception_only(type(exception), exception)
                self._stored_parse_errors.append((card, var))
                if self._iparse_errors > self._nparse_errors:
                    self.pop_parse_errors()
            #except AssertionError as exception:
                #self.log.error(str(card_obj))
                #raise
        else:
            #raise RuntimeError(card_obj)
            self.reject_cards.append(card_obj)

    def get_bdf_stats(self, return_type: str='string') -> Union[str, List[str]]:
        """
        Print statistics for the BDF

        Parameters
        ----------
        return_type : str (default='string')
            the output type ('list', 'string')
                'list' : list of strings
                'string' : single, joined string

        Returns
        -------
        return_data : str, optional
            the output data

        .. note:: if a card is not supported and not added to the proper
                  lists, this method will fail

        .. todo:: RBE3s from OP2s can show up as ???s

        """
        return get_bdf_stats(self, return_type=return_type)

    def get_displacement_index_xyz_cp_cd(self, fdtype: str='float64', idtype: str='int32',
                                         sort_ids: bool=True) -> Any:
        """
        Get index and transformation matricies for nodes with
        their output in coordinate systems other than the global.
        Used in combination with ``OP2.transform_displacements_to_global``

        Parameters
        ----------
        fdtype : str; default='float64'
            the type of xyz_cp
        idtype : str; default='int32'
            the type of nid_cp_cd
        sort_ids : bool; default=True
            sort the ids

        Returns
        -------
        icd_transform : dict{int cd : (n,) int ndarray}
            Dictionary from coordinate id to index of the nodes in
            ``self.point_ids`` that their output (`CD`) in that
            coordinate system.
        icp_transform : dict{int cp : (n,) int ndarray}
            Dictionary from coordinate id to index of the nodes in
            ``self.point_ids`` that their input (`CP`) in that
            coordinate system.
        xyz_cp : (n, 3) float ndarray
            points in the CP coordinate system
        nid_cp_cd : (n, 3) int ndarray
            node id, CP, CD for each node

        Examples
        --------
        # assume GRID 1 has a CD=10, CP=0
        # assume GRID 2 has a CD=10, CP=0
        # assume GRID 5 has a CD=50, CP=0
        >>> model.point_ids
        [1, 2, 5]
        >>> out = model.get_displacement_index_xyz_cp_cd()
        >>> icd_transform, icp_transform, xyz_cp, nid_cp_cd = out
        >>> nid_cp_cd
        [
           [1, 0, 10],
           [2, 0, 10],
           [5, 0, 50],
        ]
        >>> icd_transform[10]
        [0, 1]

        >>> icd_transform[50]
        [2]

        """
        nnodes = len(self.nodes)
        nspoints = 0
        nepoints = 0
        spoints = None
        epoints = None
        nrings = len(self.ringaxs)
        if self.spoints:
            spoints = list(self.spoints)
            nspoints = len(spoints)
        if self.epoints:
            epoints = list(self.epoints)
            nepoints = len(epoints)

        ngridb = len(self.gridb)
        if nnodes + nspoints + nepoints + ngridb + nrings == 0:
            msg = 'nnodes=%s nspoints=%s nepoints=%s nrings=%s' % (
                nnodes, nspoints, nepoints, nrings)
            raise ValueError(msg)

        if idtype == 'int32':
            try:
                out = _set_nodes(self, spoints, epoints,
                                 nnodes, nspoints, nepoints, ngridb,
                                 idtype, fdtype)
            except OverflowError:
                out = _set_nodes(self, spoints, epoints,
                                 nnodes, nspoints, nepoints, ngridb,
                                 'int64', fdtype)
        else:
            out = _set_nodes(self, spoints, epoints,
                             nnodes, nspoints, nepoints, ngridb,
                             idtype, fdtype)
        nid_cp_cd, xyz_cp, nids_cd_transform, nids_cp_transform = out

        if sort_ids:
            nids = nid_cp_cd[:, 0]
            isort = nids.argsort()
            nid_cp_cd = nid_cp_cd[isort, :]
            xyz_cp = xyz_cp[isort, :]

        icp_transform = {}
        icd_transform = {}
        nids_all = nid_cp_cd[:, 0]

        # get the indicies of the xyz array where the nodes that
        # need to be transformed are
        for cd, nids in sorted(nids_cd_transform.items()):
            if cd in [0, -1]:
                continue
            nids = np.array(nids)
            icd_transform[cd] = np.where(np.in1d(nids_all, nids))[0]

        for cp, nids in sorted(nids_cp_transform.items()):
            if cp in [-1]:
                continue
            nids = np.array(nids)
            icp_transform[cp] = np.where(np.in1d(nids_all, nids))[0]
        return icd_transform, icp_transform, xyz_cp, nid_cp_cd

    def get_xyz_in_coord_array(self, cid: int=0,
                               fdtype: str='float64', idtype: str='int32') -> Tuple[Any, Any, Any, Dict[int, Any], Dict[int, Any]]:
        """
        Gets the xyzs as an array in an arbitrary coordinate system

        Parameters
        ----------
        fdtype : str; default='float64'
            the type of xyz_cp
        idtype : str; default='int32'
            the type of nid_cp_cd
        cid : int; default=0
            the coordinate system to get xyz in

        Returns
        -------
        nid_cp_cd : (n, 3) int ndarray
            node id, CP, CD for each node
        xyz_cid : (n, 3) float ndarray
            points in the CID coordinate system
        xyz_cp : (n, 3) float ndarray
            points in the CP coordinate system
        icd_transform : dict{int cd : (n,) int ndarray}
            Dictionary from coordinate id to index of the nodes in
            ``self.point_ids`` that their output (`CD`) in that
            coordinate system.
        icp_transform : dict{int cp : (n,) int ndarray}
            Dictionary from coordinate id to index of the nodes in
            ``self.point_ids`` that their input (`CP`) in that
            coordinate system.

        .. todo:: how are SPOINTs/EPOINTs identified?

        Examples
        --------
        >>> out = model.get_xyz_in_coord_array(cid=0, fdtype='float64', idtype='int32')
        >>> nid_cp_cd, xyz_cid, xyz_cp, icd_transform, icp_transform = out
        """
        icd_transform, icp_transform, xyz_cp, nid_cp_cd = self.get_displacement_index_xyz_cp_cd(
            fdtype=fdtype, idtype=idtype, sort_ids=True)
        nids = nid_cp_cd[:, 0]
        xyz_cid = self.transform_xyzcp_to_xyz_cid(xyz_cp, nids, icp_transform,
                                                  cid=cid, in_place=False, atol=1e-6)
        return nid_cp_cd, xyz_cid, xyz_cp, icd_transform, icp_transform

    def transform_xyzcp_to_xyz_cid(self, xyz_cp: Any, nids: Any, icp_transform: Any,
                                   cid: int=0, in_place: bool=False, atol: float=1e-6) -> Any:
        """
        Vectorized method for calculating node locations in an arbitrary
        coordinate system.

        Parameters
        ----------
        xyz_cp : (n, 3) float ndarray
            points in the CP coordinate system
        nids : (n, ) int ndarray
            the GRID/SPOINT/EPOINT ids corresponding to xyz_cp
        icp_transform : dict{int cp : (n,) int ndarray}
            Dictionary from coordinate id to index of the nodes in
            ``self.point_ids`` that their input (`CP`) in that
            coordinate system.
        cid : int; default=0
            the coordinate system to get xyz in
        in_place : bool, default=False
            If true the original xyz_cp is modified, otherwise a
            new one is created.

        Returns
        -------
        xyz_cid : (n, 3) float ndarray
            points in the CID coordinate system

        Examples
        --------
        # assume GRID 1 has a CD=10, CP=0
        # assume GRID 2 has a CD=10, CP=0
        # assume GRID 5 has a CD=50, CP=1
        >>> model.point_ids
        [1, 2, 5]
        >>> out = model.get_displacement_index_xyz_cp_cd()
        >>> icd_transform, icp_transform, xyz_cp, nid_cp_cd = out
        >>> nids = nid_cp_cd[:, 0]
        >>> xyz_cid0 = model.transform_xyzcp_to_xyz_cid(
                xyz_cp, nids, icp_transform,
                cid=0)
        >>> xyz_cid1 = model.transform_xyzcp_to_xyz_cid(
                xyz_cp, nids, icp_transform,
                cid=1)

        """
        #F:\work\pyNastran\examples\femap_examples\Support\nast\tpl\heli112em7.dat
        if self.is_bdf_vectorized:
            # this is used when xref=False (only for vectorized=True)
            # we now require nids, where the other approach
            # (the one with xref=True) does not
            in_place = False
            cps_to_check = list(self.coords.keys())
        else:
            # this requires xref
            #cps_to_check = list(icp_transform.keys())
            # xref allows in_place=True

            # this is more general and slightly slower
            # requires in_place=False???
            cps_to_check = list(self.coords.keys())
        cps_to_check.sort()
        assert 0 in cps_to_check, cps_to_check

        nodes = self.nodes
        coord2 = self.coords[cid]
        #assert in_place is False, 'in_place=%s' % in_place
        if in_place:
            xyz_cid0 = xyz_cp
        else:
            xyz_cid0 = np.copy(xyz_cp)

        do_checks = False
        xyz_cid0_correct = None
        if do_checks:
            # transform the grids to the global coordinate system
            xyz_cid0_correct = self.get_xyz_in_coord(fdtype=xyz_cid0.dtype, cid=0)

        #cps_to_check = list(icp_transform.keys())
        #ncoords_to_setup = len(icp_transform)
        ncoords_to_setup = len(cps_to_check)
        nids_checked = []
        while ncoords_to_setup > 0:
            #print('--------------------------------------------------------------------------')
            #print('ncoords_to_setup = ', ncoords_to_setup)
            ncoords_to_setup = 0
            nids_checkedi, cps_checked, cps_to_check = self._transform(
                cps_to_check, icp_transform,
                nids, xyz_cp, xyz_cid0, xyz_cid0_correct,
                in_place, do_checks)

            if cps_to_check:
                nids_checked.append(nids_checkedi)
                #print("nids_checkedi =", nids_checkedi)
                out = self._get_coords_to_update(cps_to_check, cps_checked, nids_checked)
                unused_ncoords_to_setup, cord1s_to_update, cord2s_to_update, nids_checked = out
                #print('CPs not handled=%s\n  cord1s_to_update=%s\n  cord2s_to_update=%s' % (
                    #cps_to_check, cord1s_to_update, cord2s_to_update))

                for cp in cord2s_to_update:
                    coord = self.coords[cp]
                    coord.rid_ref = self.coords[coord.rid]
                    coord.setup_no_xref(self)

                for cp in cord1s_to_update:
                    coord = self.coords[cp]
                    nid1, nid2, nid3 = coord.node_ids
                    if self.is_bdf_vectorized or 1:
                        i1, i2, i3 = np.searchsorted(nids, coord.node_ids)
                        assert nids[i1] == nid1
                        assert nids[i2] == nid2
                        assert nids[i3] == nid3
                        coord.e1 = xyz_cid0[i1, :] #: the origin in the local frame
                        coord.e2 = xyz_cid0[i2, :] #: a point on the z-axis
                        coord.e3 = xyz_cid0[i3, :] #: a point on the xz-plane
                    else:
                        g1_ref = nodes[nid1]
                        g2_ref = nodes[nid2]
                        g3_ref = nodes[nid3]
                        coord.e1 = g1_ref.get_position() #: the origin in the local frame
                        coord.e2 = g2_ref.get_position() #: a point on the z-axis
                        coord.e3 = g3_ref.get_position() #: a point on the xz-plane
                    coord.setup_no_xref(self)
                    #coord.rid_ref = self.coords[coord.rid]
                    #coord.setup_no_xref(self)

                ncoords_to_setup = len(cord1s_to_update) + len(cord2s_to_update)
            #print('ncoords_next = ', ncoords_to_setup)
        #print('--------------------------------------------------------------------------')
        #print('ncoords_to_setup = ', ncoords_to_setup)

        #if ncoords == 0:
        if cps_to_check:
            msg = 'CPs not handled=%s cord1s_to_update=%s cord2s_to_update=%s\n' % (
                cps_to_check, cord1s_to_update, cord2s_to_update)
            for cp in cps_to_check:
                coord = self.coords[cp]
                msg += coord.rstrip() + '\n'
                if coord.type in ['CORD2R', 'CORD2C', 'CORD2S']:
                    rid_ref = self.coords[coord.rid]
                    msg += rid_ref.rstrip() + '\n'
                    msg += '  rid=%r origin=%s\n\n' % (coord.rid, rid_ref.origin)
                else:
                    nid1, nid2, nid3 = coord.node_ids
                    #coord.e1 = xyz_cid0[i1, :] #: the origin in the local frame
                    #coord.e2 = xyz_cid0[i2, :] #: a point on the z-axis
                    #coord.e3 = xyz_cid0[i3, :] #: a point on the xz-plane
                    if self.is_bdf_vectorized:
                        i1, i2, i3 = np.searchsorted(nids, coord.node_ids)
                        cp1 = nodes.cp[i1]
                        cp2 = nodes.cp[i2]
                        cp3 = nodes.cp[i3]
                    else:
                        cp1 = nodes[nid1].cp
                        cp2 = nodes[nid2].cp
                        cp3 = nodes[nid3].cp
                    msg += '  g1=%s xyz=%s cp=%s\n' % (nid1, coord.e1, cp1)
                    msg += '  g2=%s xyz=%s cp=%s\n' % (nid2, coord.e2, cp2)
                    msg += '  g3=%s xyz=%s cp=%s\n' % (nid3, coord.e3, cp3)
                    #break
            raise RuntimeError(msg)

        if do_checks and not np.allclose(xyz_cid0, xyz_cid0_correct, atol=atol):
            #np.array_equal(xyz_cid, xyz_cid_alt):
            out = self.get_displacement_index_xyz_cp_cd(fdtype=xyz_cid0.dtype, sort_ids=True)
            unused_icd_transform, icp_transform, xyz_cp, nid_cp_cd = out
            msg = ('xyz_cid0:\n%s\n'
                   'xyz_cid0_correct:\n%s\n'
                   'nid_cp_cd:\n%s\n'
                   'xyz_cp:\n%s'% (xyz_cid0, xyz_cid0_correct, nid_cp_cd, xyz_cp))
            raise ValueError(msg)

        if cid == 0:
            return xyz_cid0

        # transform the grids to the local coordinate system
        #is_beta = np.diagonal(beta2).min() != 1.
        #is_origin = np.abs(coord2.origin).max() != 0.
        #if is_beta and is_origin:
        xyz_cid = coord2.transform_node_to_local_array(xyz_cid0)
        #xyz_cid = coord2.xyz_to_coord_array(np.dot(xyz_cid0 - coord2.origin, beta2.T))
        #elif is_beta:
            #xyz_cid = coord2.xyz_to_coord_array(xyz_cid0 @ beta2.T)
        #else:
            #xyz_cid = coord2.xyz_to_coord_array(xyz_cid0 - coord2.origin)

        if atol is not None:
            xyz_cid_correct = self.get_xyz_in_coord(cid=cid)
            if not np.allclose(xyz_cid, xyz_cid_correct, atol=atol):
                #np.array_equal(xyz_cid, xyz_cid_correct):
                msg = ('xyz_cid:\n%s\n'
                       'xyz_cid_correct:\n%s'% (xyz_cid, xyz_cid_correct))
                raise ValueError(msg)
        return xyz_cid

    def _transform(self, cps_to_check0, icp_transform,
                   nids, xyz_cp, xyz_cid0, xyz_cid0_correct,
                   unused_in_place, do_checks):
        """
        Transforms coordinates in a vectorized way
        Helper method for ``transform_xyzcp_to_xyz_cid``

        Parameters
        ----------
        cps_to_check0 : List[int]
            the Cps to check
        icp_transform : dict{int cp : (n,) int ndarray}
            Dictionary from coordinate id to index of the nodes in
            ``self.point_ids`` that their input (`CP`) in that
            coordinate system.
        nids : (n, ) int ndarray
            the GRID/SPOINT/EPOINT ids corresponding to xyz_cp
        xyz_cp : (n, 3) float ndarray
            points in the CP coordinate system
        xyz_cid : (n, 3) float ndarray
            points in the CID coordinate system
        xyz_cid_correct : (n, 3) float ndarray
            points in the CID coordinate system
        unused_in_place : bool, default=False
            If true the original xyz_cp is modified, otherwise a
            new one is created.
        do_checks : bool; default=False
            internal value for testing
            True : makes use of xyz_cid_correct
            False : xyz_cid_correct is unused

        Returns
        -------
        nids_checked : (nnodes_checked,) int ndarray
           the node ids that were checked
        cps_checked : List[int]
            the Cps that were checked
        cps_to_check : List[int]
            the Cps that are unreferenceable given the current information

        """
        nids_checked, cps_checked, cps_to_check = transform_coords_vectorized(
            cps_to_check0, icp_transform,
            nids, xyz_cp, xyz_cid0, xyz_cid0_correct,
            self.coords, do_checks)
        return nids_checked, cps_checked, cps_to_check

    def _get_coords_to_update(self, cps_to_check, cps_checked, nids_checked):
        """helper method for ``transform_xyzcp_to_xyz_cid``"""
        cord1s_to_update_temp = []
        cord2s_to_update = []
        for cp in sorted(cps_to_check):
            coord = self.coords[cp]
            if coord.type in ['CORD2R', 'CORD2C', 'CORD2S']:
                if coord.rid in cps_checked:
                    cord2s_to_update.append(cp)
            elif coord.type in ['CORD1R', 'CORD1C', 'CORD1S']:
                cord1s_to_update_temp.append(cp)
            else:
                raise NotImplementedError(coord.rstrip())

        cord1s_to_update = set()
        if cord1s_to_update_temp:
            if len(nids_checked) == 0:
                raise RuntimeError('len(nids_checked)=0...this shouldnt happen.')
            elif len(nids_checked) == 1:
                pass
            else:
                nids_checked = [np.hstack(nids_checked)]

            nids_checkedi = nids_checked[0]
            if len(nids_checkedi) == 0:
                #print("no cord1s to check...")
                cord1s_to_update = []
            else:
                #print('nids_checked = ', nids_checkedi)
                for cp in cord1s_to_update_temp:
                    coord = self.coords[cp]
                    nids = coord.node_ids
                    #print('cp=%s nids=%s' % (cp, nids))
                    for nid in nids:
                        if nid not in nids_checkedi:
                            #print('  nid=%s break...' % nid)
                            break
                    else:
                        #print('  passed')
                        # all nids passed
                        cord1s_to_update.add(cp)
                cord1s_to_update = list(cord1s_to_update)
                cord1s_to_update.sort()

        ncoords = len(cord1s_to_update) + len(cord2s_to_update)
        #if ncoords == 0:
            #msg = 'CPs not handled=%s cord1s_to_update=%s cord2s_to_update=%s\n' % (
                #cps_to_check, cord1s_to_update, cord2s_to_update)
            #for cp in (cord1s_to_update + cord2s_to_update):
                #msg += str(cp)
            #raise RuntimeError(msg)
        return ncoords, cord1s_to_update, cord2s_to_update, nids_checked

    @property
    def is_bdf_vectorized(self):
        """Returns False for the ``BDF`` class"""
        return hasattr(self, 'grid')

    def get_displacement_index(self) -> Tuple[Any, Any, Dict[int, Any]]:
        """
        Get index and transformation matricies for nodes with
        their output in coordinate systems other than the global.
        Used in combination with ``OP2.transform_displacements_to_global``

        Returns
        -------
        nids_all : (nnodes,) int ndarray
            the GRID/SPOINT/EPOINT ids
        nids_transform : dict[cd] : (nnodesi,) int ndarray
            the indicies in nids_all that correspond to cd > 0
            cd : int
                the CD coordinate system
            nnodesi : int
                nnodesi <= nnodes
        icd_transform : dict{int cid : (n,) int ndarray}
            Dictionary from coordinate id to index of the nodes in
            ``self.point_ids`` that their output (`CD`) in that
            coordinate system.

        Examples
        --------
        # assume GRID 1 has a CD=10
        # assume GRID 2 has a CD=10
        # assume GRID 5 has a CD=50
        >>> model.point_ids
        [1, 2, 5]
        >>> icd_transform = model.get_displacement_index()
        >>> icd_transform[10]
        [0, 1]

        >>> icd_transform[50]
        [2]

        """
        nids_transform = defaultdict(list)
        icd_transform = {}
        if len(self.coords) == 1:  # was ncoords > 2; changed b/c seems dangerous
            return icd_transform

        for nid, node in sorted(self.nodes.items()):
            cid_d = node.Cd()
            if cid_d:
                nids_transform[cid_d].append(nid)

        nids_all = np.array(sorted(self.point_ids))
        for cid in sorted(nids_transform.keys()):
            nids = np.array(nids_transform[cid])
            icd_transform[cid] = np.where(np.in1d(nids_all, nids))[0]
        return nids_all, nids_transform, icd_transform

    def increase_card_count(self, card_name: str, count_num: int=1) -> None:
        """
        Used for testing to check that the number of cards going in is the
        same as each time the model is read verifies proper writing of cards

        Parameters
        ----------
        card_name : str
            the card_name -> 'GRID'
        count_num : int, optional
            the amount to increment by (default=1)

        >>> bdf.read_bdf(bdf_filename)
        >>> bdf.card_count['GRID']
        50

        """
        assert '=' not in card_name, card_name
        if card_name in self.card_count:
            self.card_count[card_name] += count_num
        else:
            self.card_count[card_name] = count_num

    def _old_card_fields(self, card_lines: List[str], card_name: str,
                         log: Any,
                         is_list: bool=False, has_none: bool=True,
                         is_dynamic_syntax: bool=False) -> BDFCard:
        """replication helper"""
        if is_list:
            fields = card_lines
        else:
            fields = to_fields(card_lines, card_name)

        # apply OPENMDAO syntax
        if is_dynamic_syntax:
            #_parse_dynamic_syntax(key, dict_of_vars, log)
            fields = [print_field_16(_parse_dynamic_syntax(field, self.dict_of_vars, log))
                      if '%' in field.strip()[0:1] else print_field_16(field)
                      for field in fields]
            has_none = False

        if has_none:
            card = wipe_empty_fields([print_field_16(field) for field in fields])
        else:
            card = wipe_empty_fields(fields)
        card_obj = BDFCard(card, has_none=False)
        return card_obj

    def _expand_replication(self, card_name, icard, cards_list, card_lines_new, dig=True):
        """replication helper"""
        #dig_str = '  ' if dig is False else ''
        #print(dig_str, '-----------************---------')
        #print(dig_str, '--dig=%s--' % dig)
        #print(dig_str, 'card_lines_new=%s' % card_lines_new)
        card = []
        cards = []
        card_lines_old = cards_list[icard-1][2]

        is_star_lines = any('*' in line for line in card_lines_old)
        if is_star_lines:
            #new_fields = to_fields_replication(card_lines_old)
            old_card = to_fields_replication(card_lines_old)
        else:
            #old_card, unused_card = self.create_card_object(
                #card_lines_old, card_name,
                #is_list=False, has_none=True)
            #print(card_lines_old)
            old_card = self._old_card_fields(card_lines_old, card_name, self.log,
                                             is_list=False, has_none=True,
                                             is_dynamic_syntax=self._is_dynamic_syntax)
            #print(old_card)
            #assert '=' not in card_name

        nlines = len(card_lines_new)
        #print(dig_str, "card_lines_new =", card_lines_new)
        new_card = to_fields_replication(card_lines_new)
        assert len(card_lines_new) == nlines, card_lines_new

        #print(dig_str, 'old_card =', old_card)
        #print(dig_str, 'card_name = %r' % card_name)
        old_card_real = None
        if old_card[0] == '=':
            #print(dig_str, 'A!!!')
            if dig is False:
                raise ReplicationError(f'dig=False...old_card=\n{old_card}')

            cards2 = self._expand_replication(
                card_name, icard-1, cards_list, card_lines_old, dig=False)
            assert len(cards2) == 1, 'cards2=%s; ncards=%s' % (cards2, len(cards2))
            #print(dig_str, 'cards_equal =', cards2)
            old_card_fields = cards2[0]
            old_card_real = old_card
            #print(dig_str, 'old_card_fields =', old_card_fields)
            #print(dig_str, 'old_card_real =', old_card_real)
            #print(dig_str, 'card_lines_old =', card_lines_old)
            old_card = self._old_card_fields(old_card_fields, card_name, self.log,
                                             is_list=True, has_none=True,
                                             is_dynamic_syntax=self._is_dynamic_syntax)
        elif '=' in card_name:
            #print(dig_str, 'B!!!')
            #print(dig_str, 'old_card =', old_card)
            #print(dig_str, 'card_lines_new =', card_lines_new)
            #print(dig_str, 'card_name = %r' % card_name)

            # good
            #new_card = [u'=3']
            #old_card = [u'CQUAD4', u'64', u'1', u'88', u'89', u'101', u'100']
            #old_card_real = [u'=', u'*1', u'=', u'*1', u'*1', u'*1', u'*1']

            # bad
            #card_name = u'=(7)'
            #new_card = [u'=(7)', u'*(10)', u'=', u'=', u'=', u'*(1.0)']
            #old_card = [u'grid', u'1001', None, u'0.', u'0.', u'0.']
            #old_card_real = [u'grid', u'1001', None, u'0.', u'0.', u'0.']
            old_card_real = new_card
        #else:
            #print('old_card[0] %r' % old_card[0])

        #print(dig_str, "old_card =", old_card)
        #print(dig_str, "new_card =", new_card)
        for ifield, field in enumerate(new_card):
            if field is None:
                field2 = old_card.field(ifield)
                #print(' %i: %r -> %r' % (ifield, field, field2))
                #assert field2 is None, 'field=%s field2=%s' % (field, field2)
                card.append(field2)
                continue

            #if field == '':
                #pass

            field = field.strip()
            if field == '=':
                field2 = old_card[ifield]
                #field2 = old_card.field(ifield)
            elif field == '==':
                # just append the remaining fields
                card.extend(old_card[ifield:])
                #print(dig_str, ' %i : extending %s' % (ifield, old_card[ifield:]))
                #print(dig_str, ' break _expand_replication...')
                break
            elif '=' in field:
                # =4
                assert ifield == 0, 'ifield=%s field=%r new_card=%s' % (ifield, field, new_card)
                nrepeats = get_nrepeats(field, old_card, new_card)
                if old_card_real is None:
                    #old_card_real = old_card
                    msg = (
                        'Invalid Replication Syntax (continuations arent supported)\n'
                        'old:\n%s\n'
                        'new:\n%s'
                        % (old_card, new_card))
                    raise RuntimeError(msg)

                #new_card = [u'=3']
                #old_card = [u'CQUAD4', u'64', u'1', u'88', u'89', u'101', u'100']
                #old_card_real = [u'=', u'*1', u'=', u'*1', u'*1', u'*1', u'*1']

                #print('---')
                #print(dig_str, "nrepeats =", nrepeats)
                new_card[0] = '='
                #print(dig_str, "new_card =", new_card)
                #print(dig_str, "old_card =", old_card)
                #print(dig_str, "old_card_real =", old_card_real)
                for unused_irepeat in range(nrepeats):
                    repeated_cards = repeat_cards(old_card, old_card_real)
                    if len(repeated_cards) != 1:
                        for repeated_card in repeated_cards:
                            print("  repeated_card =", repeated_card)
                        raise RuntimeError('too many repeated cards')
                    #for repeated_card in repeated_cards:
                        #print("  repeated_card =", repeated_card)
                    repeated_card = repeated_cards[0]
                    cards.append(repeated_card)
                    old_card = repeated_card
                    #print(dig_str, "  repeated_card =", repeated_card)
                #print(dig_str, 'breaking...')
                return cards

            elif '*' in field:
                # this is an increment, not multiplication...
                old_field = _field(old_card, ifield)
                assert old_field is not None, f'old_card:{old_card}\nnew_card:\n{new_card}'
                try:
                    if '.' in field:
                        field2 = float_replication(field, old_field)
                    else:
                        field2 = int_replication(field, old_field)
                except:
                    self.log.error(f'old_card:{old_card}\nnew_card:\n{new_card}')
                    raise
            else:
                assert '(' not in field, 'field=%r' % field
                assert '*' not in field, 'field=%r' % field
                assert '=' not in field, 'field=%r' % field
                field2 = field
            #print(dig_str, ' %i: %r -> %r' % (ifield, field, field2))
            card.append(field2)
        if card:
            cards.append(card)
            #print(dig_str, 'card_expanded = %s' % card)
        else:  # pragma: no cover
            raise RuntimeError(card)
        return cards

    def _parse_cards(self, cards_list, cards_dict, card_count):
        # type: (List[List[str]], Dict[str, List[str]], Dict[str, int]) -> None
        """creates card objects and adds the parsed cards to the deck"""
        # we don't want replication markers in the card_count
        card_names_to_remove = (card_name for card_name in list(card_count.keys())
                                if '=' in card_name)
        for card_name in card_names_to_remove:
            del card_count[card_name]

        self.echo = False
        if cards_dict: # self._is_cards_dict = True
            self._parse_cards_dict(cards_dict)

        if cards_list:
            # this is the block that actually runs
            self._parse_cards_list(cards_list)

    def _parse_cards_dict(self, cards_dict: Dict[str, List[str]]) -> None:
        """parses the cards that are in dictionary format"""
        if self.save_file_structure:
            raise NotImplementedError('save_file_structure=True is not supported\n%s' % (
                list(cards_dict.keys())))

        for card_name, cards in sorted(cards_dict.items()):
            if self.is_reject(card_name):
                self.log.info('    rejecting card_name = %s' % card_name)
                for comment, card_lines, unused_ifile_iline in cards:
                    self.increase_card_count(card_name)
                    self.reject_lines.append([_format_comment(comment)] + card_lines)
            else:
                for comment, card_lines, (ifile, unused_iline) in cards:
                    self.add_card(card_lines, card_name, comment=comment, ifile=ifile,
                                  is_list=False, has_none=False)

    def _parse_cards_list(self, cards_list):
        """parses the cards that are in list format"""
        save_file_structure = self.save_file_structure
        if save_file_structure:
            for icard, card in enumerate(cards_list):
                card_name, comment, card_lines, (ifile, unused_iline) = card
                if card_name is None:
                    msg = 'card_name = %r\n' % card_name
                    msg += 'card_lines = %s' % card_lines
                    raise RuntimeError(msg)

                if '=' in card_name:
                    #print(card)
                    try:
                        replicated_cards = self._expand_replication(
                            card_name, icard, cards_list, card_lines)
                    except ReplicationError:
                        self.log.error('failed to expand %s\n%s' % (card_name, ''.join(card_lines)))
                        raise

                    _check_replicated_cards(replicated_cards)
                    for replicated_card in replicated_cards:
                        self.add_card_ifile(ifile, replicated_card, replicated_card[0],
                                            comment=comment, is_list=True, has_none=True)
                    continue

                if self.is_reject(card_name):  # pragma: no cover
                    msg = "save_file_structure=True doesn't support %s" % card_name
                    raise NotImplementedError(msg)
                    #self.reject_card_lines(card_name, card_lines, comment)
                else:
                    self.add_card_ifile(ifile, card_lines, card_name, comment=comment,
                                        is_list=False, has_none=False)

        else:
            for icard, card in enumerate(cards_list):
                card_name, comment, card_lines, (ifile, unused_iline) = card
                #print(unused_iline, card_lines[0])
                if card_name is None:
                    msg = 'card_name = %r\n' % card_name
                    msg += 'card_lines = %s' % card_lines
                    raise RuntimeError(msg)

                if '=' in card_name:
                    #print(card)
                    try:
                        replicated_cards = self._expand_replication(
                            card_name, icard, cards_list, card_lines)
                    except ReplicationError:
                        self.log.error('failed to expand %s\n%s' % (card_name, ''.join(card_lines)))
                        raise

                    _check_replicated_cards(replicated_cards)
                    for replicated_card in replicated_cards:
                        self.add_card(replicated_card, replicated_card[0], comment=comment,
                                      is_list=True, has_none=True)
                    continue

                if self.is_reject(card_name):
                    self.reject_card_lines(card_name, card_lines, comment=comment)
                else:
                    self.add_card(card_lines, card_name, comment=comment, ifile=ifile,
                                  is_list=False, has_none=False)

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

    def _parse_primary_file_header(self, bdf_filename: Union[str, StringIO]) -> None:
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
        try:
            with open(bdf_filename, 'r') as bdf_file:
                lines = bdf_file.readlines()
        except UnicodeDecodeError:
            with open(bdf_filename, 'r', errors='ignore') as bdf_file:
                line = bdf_file.readline()
                iline = 1
                while line:
                    self.log.debug(f'Line {iline}: {line.strip()}')
                    line = bdf_file.readline()
                    iline += 1
        except (AttributeError, TypeError) as error:
            if hasattr(bdf_filename, 'read') and hasattr(bdf_filename, 'write'):
                lines = bdf_filename.readlines()
                bdf_filename.seek(0)  # need to rewind the buffer!
            else:
                raise error

        self._check_pynastran_header(lines, check_header=True)
        if self.nastran_format == 'zona':
            self.zona.update_for_zona()

    def _check_pynastran_header(self, lines: List[str], check_header: bool=True) -> None:
        """updates the $pyNastran: key=value variables"""
        if not check_header:
            return
        for line in lines:
            if not line.startswith('$'):
                break

            key, value = _parse_pynastran_header(line)
            if not key:
                break

            # key/value are lowercase
            if key == 'version':
                assert value.lower() in ['msc', 'nx', 'zona', 'nasa95'], f'version={value!r} is not supported'
                self.nastran_format = value
            elif key == 'encoding':
                self._encoding = value
            elif key == 'punch':
                self.punch = _bool(value)
            elif key in ['nnodes', 'nelements']:
                pass
            elif key == 'dumplines':
                self.dumplines = _bool(value)
            elif key == 'is_superelements':
                self.is_superelements = _bool(value)
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
            elif key in ['code-block', 'code_block']:
                value = line.split('=', 1)[1]
                if not hasattr(self, 'code_block'):
                    self.code_block = ''
                    char0 = value.lstrip()[0]
                    indent = value.index(char0)
                self.code_block += value[indent:]
            else:
                raise NotImplementedError(key)

        if hasattr(self, 'code_block'):
            exec(self.code_block)

#---------------------------------------------------------------------------------------------------
    # HDF5
    def _read_bdf_cards(self, bdf_filename: Optional[str]=None,
                        punch: bool=False,
                        read_includes: bool=True, encoding: Optional[str]=None) -> None:
        """
        Read method for the bdf files

        Parameters
        ----------
        bdf_filename : str / None
            the input bdf (default=None; popup a dialog)
        punch : bool; default=False
            indicates whether the file is a punch file
        read_includes : bool; default=True
            indicates whether INCLUDE files should be read
        encoding : str; default=None -> system default
            the unicode encoding

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
        self._is_cards_dict = True

        self._read_bdf_helper(bdf_filename, encoding, punch, read_includes)
        self.log.debug('---starting BDF.read_bdf of %s---' % self.bdf_filename)
        self._parse_primary_file_header(bdf_filename)

        obj = BDFInputPy(self.read_includes, self.dumplines, self._encoding,
                         nastran_format=self.nastran_format,
                         log=self.log, debug=self.debug)
        out = obj.get_lines(bdf_filename, punch=self.punch, make_ilines=True)
        (system_lines, executive_control_lines, case_control_lines,
         bulk_data_lines, bulk_data_ilines,
         superelement_lines, superelement_ilines) = out
        self._set_pybdf_attributes(obj, save_file_structure=False)

        self.system_command_lines = system_lines
        self.executive_control_lines = executive_control_lines
        self.case_control_lines = case_control_lines

        sol, method, sol_iline = parse_executive_control_deck(executive_control_lines)
        self.update_solution(sol, method, sol_iline)

        #self._is_cards_dict = True
        if self._is_cards_dict:
            cards, card_count = self.get_bdf_cards_dict(bulk_data_lines, bulk_data_ilines)
        cards_out = self._parse_cards_hdf5(cards, card_count)
        assert isinstance(cards_out, dict), cards_out

        self.case_control_deck = CaseControlDeck(case_control_lines, log=self.log)
        self.case_control_deck.solmap_to_value = self._solmap_to_value
        self.case_control_deck.rsolmap_to_str = self.rsolmap_to_str

        return cards_out

    def create_subcases(self, subcase_ids: Union[int, List[int], None]=None) -> Dict[int, Subcase]:
        """creates a series of subcases"""
        if subcase_ids is None:
            subcase_ids = []
        elif isinstance(subcase_ids, int):
            subcase_ids = [subcase_ids]

        if self.case_control_deck is None:
            self.case_control_deck = CaseControlDeck([], log=self.log)

        subcases = {}
        for subcase_id in subcase_ids:
            subcase = self.case_control_deck.create_new_subcase(subcase_id)
            subcases[subcase_id] = subcase
        return subcases

    def _parse_cards_hdf5(self, cards: Dict[str, Any], unused_card_count):
        """creates card objects and adds the parsed cards to the deck"""
        self.echo = False
        cards_out = {}
        for card_name, card in sorted(cards.items()):
            cards_list = []
            cards_out[card_name] = cards_list
            if self.is_reject(card_name):
                self.log.info('    rejecting card_name = %s' % card_name)
                for comment, card_lines, unused_ifile_iline in card:
                    self.increase_card_count(card_name)
                    self.reject_lines.append([_format_comment(comment)] + card_lines)
            else:
                for comment, card_lines, unused_ifile_iline in card:
                    class_instance = self._add_card_hdf5(card_lines, card_name, comment=comment,
                                                         is_list=False, has_none=False)
                    cards_list.append(class_instance)
        return cards_out


    def _add_card_hdf5(self, card_lines: List[str], card_name: str,
                       comment='', is_list: bool=True, has_none: bool=True) -> Any:
        """
        Creates a BaseCard object that will be used to simplify HDF5 adding.

        Parameters
        ----------
        card_lines: list[str]
            the list of the card fields
        card_name : str
            the card_name -> 'GRID'
        comment : str
            an optional the comment for the card
        is_list : bool, optional
            False : input is a list of card fields -> ['GRID', 1, None, 3.0, 4.0, 5.0]
            True :  input is list of card_lines -> ['GRID, 1,, 3.0, 4.0, 5.0']
        has_none : bool; default=True
            can there be trailing Nones in the card data (e.g. ['GRID', 1, 2, 3.0, 4.0, 5.0, None])
            can there be trailing Nones in the card data (e.g. ['GRID, 1, 2, 3.0, 4.0, 5.0, '])

        Returns
        -------
        class_instance : BaseCard()
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

           # here is_list=False because it's been parsed
           # Note the None at the end of the 1st line, which is there
           #      because the CONM2 card has a blank field.
           #      It must be there.
           # We also set i32 on the 2nd line, so it will default to 0.0
           >>> card = [
                   'CONM2', eid, nid, cid, mass, x1, x2, x3, None,
                            i11, i21, i22, i31, None, i33,
               ]
           >>> model.add_card(card_lines, 'CONM2', comment, is_list=False)

           # here's an alternate approach for the CONM2
           # we use Nastran's CSV format
           # There are many blank fields, but it's parsed exactly like a
           # standard CONM2.
           >>> card = [
                   'CONM2,1,2,3,10.0',
                   ',1.0,,5.0'
               ]
           >>> model.add_card(card_lines, 'CONM2', comment, is_list=True)

        Notes
        -----
        This is a very useful method for interfacing with the code.

        The card_object is not a card-type object...so not a GRID
        card or CQUAD4 object.  It's a BDFCard Object.  However,
        you know the type (assuming a GRID), so just call the
        *mesh.Node(nid)* to get the Node object that was just
        created.

        """
        card_name = card_name.upper()
        card_obj, unused_card = self.create_card_object(
            card_lines, card_name,
            is_list=is_list, has_none=has_none)
        class_instance = self._add_card_helper_hdf5(card_obj, card_name, card_name, comment)
        return class_instance

    def _add_card_helper_hdf5(self, card_obj, card, card_name, comment=''):
        # type: (BDFCard, List[str], str, str) -> None
        """
        Adds a card object to the BDF object.

        Parameters
        ----------
        card_object : BDFCard()
            the card object representation of card
        card : List[str]
            the fields of the card object; used for rejection and special cards
        card_name : str
            the card_name -> 'GRID'
        comment : str
            an optional the comment for the card

        Returns
        -------
        class_instance : obj or None
            obj : GRID, CQUAD4, ...
            None : ECHOON, ECHOOFF, reject

        """
        if card_name == 'ECHOON':
            self.echo = True
            return None
        elif card_name == 'ECHOOFF':
            self.echo = False
            return None

        if self.echo and not self.force_echo_off:
            try:
                print(print_card_8(card_obj).rstrip())
            except:
                if card in ['DEQATN']:
                    print(str(card_obj).rstrip())
                else:
                    print(print_card_16(card_obj).rstrip())

        if card_name in self._card_parser:
            card_class, add_card_function = self._card_parser[card_name]

            # simplified, so no error catching
            class_instance = card_class.add_card(card_obj, comment=comment)
            #add_card_function(class_instance)

        elif card_name in self._card_parser_prepare:
            add_card_function = self._card_parser_prepare[card_name]
            # simplified, so no error catching
            class_instance = add_card_function(card, card_obj, comment=comment)

        else:
            self.reject_cards.append(card_obj)
            class_instance = None
        return class_instance

class BDF(BDF_):
    """NASTRAN BDF Reader/Writer/Editor class."""
    _properties = ['is_bdf_vectorized', 'nid_map', 'wtmass', 'type_slot_str'] + [
        'nastran_format', 'is_long_ids', 'sol', 'subcases',
        'nnodes', 'node_ids', 'point_ids', 'npoints',
        'nelements', 'element_ids', 'nproperties', 'property_ids',
        'nmaterials', 'material_ids', 'ncoords', 'coord_ids',
        'ncaeros', 'caero_ids', 'wtmass', 'is_bdf_vectorized', 'nid_map',
        #'dmigs', 'dmijs', 'dmiks', 'dmijis', 'dtis', 'dmis',
    ]

    def __init__(self, debug: Optional[bool]=True, log: Any=None, mode: str='msc') -> None:
        """
        Initializes the BDF object

        Parameters
        ----------
        debug : bool/None; default=True
            used to set the logger if no logger is passed in
                True:  logs debug/info/warning/error messages
                False: logs info/warning/error messages
                None:  logs warning/error messages
        log : logging module object / None
            if log is set, debug is ignored and uses the
            settings the logging object has
        mode : str; default='msc'
            the type of Nastran
            valid_modes = {'msc', 'nx'}

        """
        BDF_.__init__(self, debug=debug, log=log, mode=mode)
        #: stores SPOINT, GRID cards
        self.nodes = {}  # type: Dict[int, Any]

        # loads
        #: stores LOAD, FORCE, FORCE1, FORCE2, MOMENT, MOMENT1, MOMENT2,
        #: PLOAD, PLOAD2, PLOAD4, SLOAD
        #: GMLOAD, SPCD, DEFORM,
        #: QVOL
        self.loads = {}  # type: Dict[int, List[Any]]
        self.load_combinations = {}  # type: Dict[int, List[Any]]

    def __deepcopy__(self, memo: Dict[str, Any]):
        """performs a deepcopy"""
        #newone = type(self)()
        #newone.__dict__.update(self.__dict__)
        #return newone
        cls = self.__class__
        result = cls.__new__(cls)
        memo[id(self)] = result
        for key, value in self.__dict__.items():
            setattr(result, key, deepcopy(value, memo))
        return result

    def __copy__(self):
        """performs a copy"""
        newone = type(self)()
        newone.__dict__.update(self.__dict__)
        return newone

    def _add_superelements(self, superelement_lines, superelement_ilines):
        for superelement_id, superelement_line in sorted(superelement_lines.items()):
            assert isinstance(superelement_line, list), superelement_line

            # hack to get rid of extra 'BEGIN SUPER=2' lines
            iminus = 0
            for line in superelement_line:
                uline = line.upper()
                if not uline.startswith('BEGIN '):
                    break
                iminus += 1

            nlines = len(superelement_line) - iminus
            model = BDF()
            model.active_filenames = self.active_filenames
            model.log = self.log
            model.punch = True
            #model.nastran_format = ''
            superelement_ilines = np.zeros((nlines, 2), dtype='int32')  ## TODO: calculate this
            model._parse_all_cards(superelement_line[iminus:], superelement_ilines)
            self.superelement_models[superelement_id] = model
            self.initial_superelement_models.append(superelement_id)


def _echo_card(card, card_obj):
    """echos a card"""
    try:
        print(print_card_8(card_obj).rstrip())
    except:
        if card in ['DEQATN']:
            print(str(card_obj).rstrip())
        else:
            print(print_card_16(card_obj).rstrip())

def read_bdf(bdf_filename: Optional[str]=None, validate: bool=True, xref: bool=True, punch: bool=False,
             save_file_structure: bool=False,
             skip_cards: Optional[List[str]]=None,
             read_cards: Optional[List[str]]=None,
             encoding: Optional[str]=None,
             log=None,
             debug: bool=True, mode: str='msc') -> BDF:
    # Optional[SimpleLogger]
    """
    Creates the BDF object

    Parameters
    ----------
    bdf_filename : str (default=None -> popup)
        the bdf filename
    debug : bool/None
        used to set the logger if no logger is passed in
            True:  logs debug/info/error messages
            False: logs info/error messages
            None:  logs error messages
    log : logging module object / None
        if log is set, debug is ignored and uses the
        settings the logging object has
    validate : bool; default=True
        runs various checks on the BDF
    xref :  bool; default=True
        should the bdf be cross referenced
    punch : bool; default=False
        indicates whether the file is a punch file
    save_file_structure : bool; default=False
        enables the ``write_bdfs`` method
    skip_cards : List[str]; default=None
        None : include all cards
        list of cards to skip
    read_cards : List[str]; default=None
        None : include all cards
        list of cards to read (all the cards)
    encoding : str; default=None -> system default
        the unicode encoding
    mode : str; default='msc'
        the type of Nastran
        valid_modes = {'msc', 'nx'}

    Returns
    -------
    model : BDF()
        an BDF object

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

    .. note :: this method will change in order to return an object that
               does not have so many methods
    .. todo:: finish this

    """
    model = BDF(log=log, debug=debug, mode=mode)
    if read_cards and skip_cards:
        msg = 'read_cards=%s skip_cards=%s cannot be used at the same time'
        raise NotImplementedError(msg)
    if skip_cards:
        model.disable_cards(skip_cards)
    elif read_cards:
        model.set_cards(read_cards)

    if bdf_filename and not isinstance(bdf_filename, StringIO):
        check_path(bdf_filename, 'bdf_filename')
    model.read_bdf(bdf_filename=bdf_filename, validate=validate,
                   xref=xref, punch=punch, read_includes=True,
                   save_file_structure=save_file_structure,
                   encoding=encoding)

    #if 0:
        ### TODO: remove all the extra methods

        #keys_to_suppress = []
        #method_names = model.object_methods(keys_to_skip=keys_to_suppress)

        #methods_to_remove = [
            #'_process_card', 'read_bdf', 'disable_cards', 'set_dynamic_syntax',
            #'create_card_object', 'create_card_object_fields', 'create_card_object_list',

            #'add_AECOMP', 'add_AEFACT', 'add_AELINK', 'add_AELIST', 'add_AEPARM', 'add_AERO',
            #'add_AEROS', 'add_AESTAT', 'add_AESURF', 'add_ASET', 'add_BCRPARA', 'add_BCTADD',
            #'add_BCTPARA', 'add_BCTSET', 'add_BSET', 'add_BSURF', 'add_BSURFS', 'add_CAERO',
            #'add_DIVERG',
            #'add_CSET', 'add_CSSCHD', 'add_DAREA', 'add_DCONADD', 'add_DCONSTR', 'add_DDVAL',
            #'add_DELAY', 'add_DEQATN', 'add_DESVAR', 'add_DLINK', 'add_DMI', 'add_DMIG',
            #'add_DMIJ', 'add_DMIJI', 'add_DMIK', 'add_DPHASE', 'add_DRESP', 'add_DTABLE',
            #'add_DVMREL', 'add_DVPREL', 'add_EPOINT', 'add_FLFACT', 'add_FLUTTER', 'add_FREQ',
            #'add_GUST', 'add_LSEQ', 'add_MKAERO', 'add_MONPNT', 'add_NLPARM', 'add_NLPCI',
            #'add_PAERO', 'add_PARAM', 'add_PBUSHT', 'add_PDAMPT', 'add_PELAST', 'add_PHBDY',
            #'add_QSET', 'add_SEBSET', 'add_SECSET', 'add_SEQSET', 'add_SESET', 'add_SET',
            #'add_SEUSET', 'add_SPLINE', 'add_spoint', 'add_tempd', 'add_TF', 'add_TRIM',
            #'add_TSTEP', 'add_TSTEPNL', 'add_USET',

            #'add_card', 'add_card_fields', 'add_card_lines', 'add_cmethod', 'add_constraint',
            #'add_constraint_MPC', 'add_constraint_MPCADD',
            #'add_constraint_SPC', 'add_constraint_SPCADD',
            #'add_convection_property', 'add_coord', 'add_creep_material', 'add_damper',
            #'add_dload', '_add_dload_entry', 'add_element', 'add_hyperelastic_material',
            #'add_load', 'add_mass', 'add_material_dependence', 'add_method', 'add_node',
            #'add_plotel', 'add_property', 'add_property_mass', 'add_random_table',
            #'add_rigid_element', 'add_structural_material', 'add_suport', 'add_suport1',
            #'add_table', 'add_table_sdamping', 'add_thermal_BC', 'add_thermal_element',
            #'add_thermal_load', 'add_thermal_material',

            #'set_as_msc',
            #'set_as_nx',

            #'pop_parse_errors',
            ##'pop_xref_errors',
            #'set_error_storage',
            #'is_reject',
        #]
        #for method_name in method_names:
            #if method_name not in methods_to_remove + keys_to_suppress:
                ##print(method_name)
                #pass
            #else:
                ### TODO: doesn't work...
                ##delattr(model, method_name)
                #pass
        #model.get_bdf_stats()
    return model


def _prep_comment(comment):
    return comment.rstrip()
    #print('comment = %r' % comment)
    #comment = '  this\n  is\n  a comment\n'
    #print(comment.rstrip('\n').split('\n'))
    #sline = [comment[1:] if len(comment) and comment[0] == ' ' else comment
             #for comment in comment.rstrip().split('\n')]
    #print('sline = ', sline)

def _check_for_spaces(card_name: str, card_lines: List[str], comment: str, log: Any) -> None:
    if ' ' in card_name:
        if card_name.startswith(EXECUTIVE_CASE_SPACES):  # TODO verify upper
            msg = (
                'No spaces allowed in card name %r.\n'
                'Did you mean to call read_bdf(punch=False) instead of '
                'read_bdf(punch=True)?\n%s' % (
                    card_name, card_lines))
            raise RuntimeError(msg)
        elif card_name.startswith('BEGIN '):
            uline = card_lines[0].upper()
            if 'SUPER' in uline:
                msg = (
                    'Misindentified Superelement section.  Use:\n'
                    '$ pyNastran: is_superelements=True\n')
            else:
                msg = (
                    'Is there a second BEGIN BULK in your deck?\n'
                    'Another possibility is that punch=True and there is a '
                    'BEGIN BULK in your deck.\n')
            msg += '%s\n' % card_lines
            log.error(msg)
            raise SuperelementFlagError(msg)
        else:
            msg = (
                'No spaces allowed in card name %r.\n'
                'Should this be a comment?\n%s%s' % (
                    card_name, comment, card_lines))
        raise RuntimeError(msg)

    if card_name in ['SUBCASE ', 'CEND']:
        raise RuntimeError('No executive/case control deck was defined.')

def _check_replicated_cards(replicated_cards):
    """helper method for ``parse_cards_list``"""
    replicated_card_old = []
    try:
        for replicated_card in replicated_cards:
            assert replicated_card != replicated_card_old
            replicated_card_old = replicated_card
    except AssertionError:
        #print('card_list = %s' % card_list)
        #print('card_lines = %s' % card_lines)
        replicated_card_old = []
        for replicated_card in replicated_cards:
            #print('adding ', replicated_card)
            assert replicated_card != replicated_card_old
            replicated_card_old = replicated_card
        raise


def _set_nodes(model: BDF,
               spoints, epoints,
               nnodes: int, nspoints: int, nepoints: int, ngridb: int,
               idtype, fdtype):
    """helper method for ``get_displacement_index_xyz_cp_cd``"""
    i = 0
    nids_cd_transform = defaultdict(list)  # type: Dict[int, np.ndarray]
    nids_cp_transform = defaultdict(list)  # type: Dict[int, np.ndarray]
    nxyz = nnodes + nspoints + nepoints + ngridb
    xyz_cp = np.zeros((nxyz, 3), dtype=fdtype)
    nid_cp_cd = np.zeros((nxyz, 3), dtype=idtype)
    for nid, node in sorted(model.nodes.items()):
        cd = node.Cd()
        cp = node.Cp()
        nids_cp_transform[cp].append(nid)
        nids_cd_transform[cd].append(nid)
        nid_cp_cd[i, :] = [nid, cp, cd]
        xyz_cp[i, :] = node.xyz
        i += 1
    if nspoints:
        for nid in sorted(spoints):
            nid_cp_cd[i, 0] = nid
            i += 1
    if nepoints:
        for nid in sorted(epoints):
            nid_cp_cd[i, 0] = nid
            i += 1
    if ngridb:
        for nid, node in sorted(model.gridb.items()):
            phi = node.phi
            cd = node.cd
            ringfl_id = node.ringfl
            ringfl = model.ringfl[ringfl_id]
            x = ringfl.xa  ## TODO: what about xb?
            y = phi  ## TODO: is this really phi?
            z = 0.
            axif = model.axif
            cp = axif.cid
            nid_cp_cd[i, :] = [nid, cp, cd]
            xyz_cp[i, :] = [x, y, z]
            i += 1
    return nid_cp_cd, xyz_cp, nids_cd_transform, nids_cp_transform

def _bool(value):
    """casts a lower string to a booean"""
    return True if value == 'true' else False

def main():  # pragma: no cover
    """shows off how unicode works becausee it's overly complicated"""
    import pyNastran
    pkg_path = pyNastran.__path__[0]
    bdf_filename = os.path.abspath(os.path.join(
        pkg_path, '..', 'models', 'solid_bending', 'solid_bending.bdf'))
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
    node1.comment = '$ ' + note + '\n'

    # in other words, msg is a bad comment:
    msg = '$ line 1\n'
    msg += 'line 2\n'
    msg += '$ line 3\n'
    print(msg)
    model.write_bdf('test.bdf')


#if __name__ == '__main__':  # pragma: no cover
    #from pyNastran.bdf.test.test_bdf import main
    #main()
