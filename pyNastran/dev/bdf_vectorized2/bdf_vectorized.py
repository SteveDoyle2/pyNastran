# coding: utf-8
# pylint: disable=W0201,R0915,R0912
"""
Main BDF class.  Defines:
  - BDF

"""
# TABLE3D TID X0 Y0 Z0 F0
# X1 Y1 Z1 F1 X2 Y2 Z2 F2
# X3 Y3 Z3 F3 X4 Y4 Z4 F4
# -etc.- ENDT

# see https://docs.plm.automation.siemens.com/tdoc/nxnastran/10/help/#uid:index
from __future__ import annotations

import sys
from io import StringIO, IOBase
from pathlib import PurePath
import traceback
from collections import defaultdict
from typing import (
    Sequence, TextIO, Optional, Any, TYPE_CHECKING)
from pickle import load, dump, dumps  # type: ignore

import numpy as np  # type: ignore
from cpylog import get_logger2, __version__ as CPYLOG_VERSION

from pyNastran.bdf.bdf import (LOAD, _bool, _check_replicated_cards,
                               _get_coords_to_update, map_update, map_version,
                               _prep_comment, _echo_card)
#from pyNastran.bdf.bdf import BDF as BDF_, LOAD
from .cards.nodes import GRIDv, Nodes
from .cards.elements.elements import Elements

from .cards.elements.springs import (
    CELAS1, CELAS2, CELAS3, CELAS4, Springs)
from .cards.elements.dampers import (
    CDAMP1, CDAMP2, CDAMP3, CDAMP4, CDAMP5, CVISCv, PLOTELv, Dampers)
from .cards.elements.bush import (
    CBUSHv, Bushes)

from .cards.elements.rods import (
    CONRODv, CRODv, CTUBEv, Rods)
from .cards.elements.masses import (
    CONM2v, Masses)
from .cards.elements.bars import CBARv, Bars
from .cards.elements.beams import CBEAMv, Beams
from .cards.elements.shears import CSHEARv, Shears
from .cards.elements.shells import (
    CTRIA3v, CTRIA6v, CTRIARv, CQUAD4v, CQUAD8v, CQUADv, CQUADRv, Shells)
from .cards.elements.solids import (
    CTETRA4v, CPENTA6v, CHEXA8v, CPYRAM5v,
    CTETRA10v, CPENTA15v, CHEXA20v, CPYRAM13v, Solids)

from .cards.loads.loads import (
    Loads, FORCEv, FORCE1v, FORCE2v,
    MOMENTv, MOMENT1v, MOMENT2v,
    SLOADv, SPCDv, GRAVv, LSEQv)
from .cards.loads.pressure_loads import (
    PLOADv, PLOAD1v, PLOAD2v, PLOAD4v)
from .cards.loads.thermal_loads import (
    TEMPv, TEMPDv)

# ------------------------------------------------------------------------------
from pyNastran.utils import object_attributes, check_path
from pyNastran.bdf.utils import parse_patran_syntax
from pyNastran.bdf.bdf_interface.utils import (
    _parse_pynastran_header, to_fields, parse_executive_control_deck,
    fill_dmigs, _get_card_name, _parse_dynamic_syntax
)
from pyNastran.bdf.bdf_interface.add_card import CARD_MAP
from pyNastran.bdf.bdf_interface.replication import (
    to_fields_replication, get_nrepeats, int_replication, float_replication,
    _field, repeat_cards)

#from pyNastran.bdf.field_writer_8 import print_card_8
from pyNastran.bdf.field_writer_16 import print_field_16 # print_card_16,

from pyNastran.bdf.cards.base_card import _format_comment
from pyNastran.bdf.cards.utils import wipe_empty_fields

#from .write_path import write_include
from pyNastran.bdf.bdf_interface.assign_type import (integer,
                                        integer_or_string, string)

from pyNastran.bdf.cards.elements.elements import PLOTEL # , CFAST, CGAP, CRAC2D, CRAC3D, GENEL
from pyNastran.bdf.cards.properties.properties import PFAST, PGAP, PRAC2D, PRAC3D
from pyNastran.bdf.cards.properties.solid import PLSOLID, PSOLID, PIHEX, PCOMPS
#from pyNastran.bdf.cards.cyclic import CYAX, CYJOIN
#from pyNastran.bdf.cards.msgmesh import CGEN

#from pyNastran.bdf.cards.elements.springs import CELAS1, CELAS2, CELAS3, CELAS4
from pyNastran.bdf.cards.properties.springs import PELAS, PELAST

#from pyNastran.bdf.cards.elements.solid import (
    #CTETRA, CPYRAM, CPENTA, CHEXA,
    #CIHEX1, CIHEX2, CHEXA1, CHEXA2,
    #CTETRA4, CPYRAM5, CPENTA6, CHEXA8,
    #CTETRA10, CPYRAM13, CPENTA15, CHEXA20,
#)
from pyNastran.bdf.cards.elements.rigid import RBAR, RBAR1, RBE1, RBE2, RBE3, RROD, RSPLINE, RSSCON

#from pyNastran.bdf.cards.axisymmetric.axisymmetric import (
    #AXIF, RINGFL,
    #AXIC, RINGAX, POINTAX, CCONEAX, PCONEAX, )
#from pyNastran.bdf.cards.axisymmetric.loads import PLOADX1, FORCEAX, PRESAX, TEMPAX
#from pyNastran.bdf.cards.elements.axisymmetric_shells import (
    #CTRAX3, CTRAX6, CTRIAX, CTRIAX6, CQUADX, CQUADX4, CQUADX8)
#from pyNastran.bdf.cards.elements.shell import (
    #CQUAD, CQUAD4, CQUAD8, CQUADR, CSHEAR,
    #CTRIA3, CTRIA6, CTRIAR,
    #CPLSTN3, CPLSTN4, CPLSTN6, CPLSTN8,
    #CPLSTS3, CPLSTS4, CPLSTS6, CPLSTS8,
    #SNORM,)
#from pyNastran.bdf.cards.elements.shell_nasa95 import (
    #CTRSHL, CQUAD1, PQUAD1)

from pyNastran.bdf.cards.properties.shell import PSHELL, PCOMP, PCOMPG, PSHEAR, PLPLANE, PPLANE, PTRSHL
from pyNastran.bdf.cards.elements.acoustic import (
    CHACAB, CAABSF, CHACBR, PACABS, PAABSF, PACBAR, ACMODL)
#from pyNastran.bdf.cards.elements.bush import CBUSH, CBUSH1D, CBUSH2D
from pyNastran.bdf.cards.properties.bush import PBUSH, PBUSH1D, PBUSHT, PBUSH_OPTISTRUCT
#from pyNastran.bdf.cards.elements.damper import (CVISC, CDAMP1, CDAMP2, CDAMP3, CDAMP4,
                                    #CDAMP5)
from pyNastran.bdf.cards.properties.damper import PVISC, PDAMP, PDAMP5, PDAMPT
#from pyNastran.bdf.cards.elements.rods import CROD, CONROD, CTUBE
#from pyNastran.bdf.cards.elements.bars import CBAR, BAROR, CBARAO, CBEAM3, CBEND
#from pyNastran.bdf.cards.elements.beam import CBEAM, BEAMOR
from pyNastran.bdf.cards.properties.rods import PROD, PTUBE
from pyNastran.bdf.cards.properties.bars import PBAR, PBARL, PBRSECT, PBEND, PBEAM3
from pyNastran.bdf.cards.properties.beam import PBEAM, PBEAML, PBCOMP, PBMSECT
## CMASS5
#from pyNastran.bdf.cards.elements.mass import CONM1, CONM2, CMASS1, CMASS2, CMASS3, CMASS4
#from pyNastran.bdf.cards.properties.mass import PMASS, NSM, NSM1, NSML, NSML1, NSMADD
from pyNastran.bdf.cards.constraints import (SPC, SPCADD, SPCAX, SPC1, SPCOFF, SPCOFF1,
                                             MPC, MPCADD, SUPORT1, SUPORT, SESUP,
                                             GMSPC)
from pyNastran.bdf.cards.coordinate_systems import (CORD1R, CORD1C, CORD1S,
                                                    CORD2R, CORD2C, CORD2S, #CORD3G,
                                                    transform_coords_vectorized,
                                                    Coord)
#from pyNastran.bdf.cards.msgmesh import CGEN, GMCORD, GMLOAD
from pyNastran.bdf.cards.deqatn import DEQATN
from pyNastran.bdf.cards.dynamic import (
    DELAY, DPHASE, FREQ, FREQ1, FREQ2, FREQ3, FREQ4, FREQ5,
    TSTEP, TSTEP1, TSTEPNL, NLPARM, NLPCI, TF, ROTORG, ROTORD, TIC)
from pyNastran.bdf.cards.loads.loads import (
    LSEQ, SLOAD, DAREA, RFORCE, RFORCE1, SPCD, DEFORM, LOADCYN, LOADCYH)
from pyNastran.bdf.cards.loads.dloads import ACSRCE, DLOAD, TLOAD1, TLOAD2, RLOAD1, RLOAD2
from pyNastran.bdf.cards.loads.static_loads import (LOAD, CLOAD) # ACCEL, ACCEL1,
                                       #GMLOAD)
from pyNastran.bdf.cards.loads.random_loads import RANDPS, RANDT1

from pyNastran.bdf.cards.materials import (MAT1, MAT2, MAT3, MAT4, MAT5,
                                           MAT8, MAT9, MAT10, MAT11, MAT3D,
                                           MATG, MATHE, MATHP, MATEV,
                                           CREEP, EQUIV, NXSTRAT)
from pyNastran.bdf.cards.material_deps import (
    MATT1, MATT2, MATT3, MATT4, MATT5, MATT8, MATT9, MATS1)

from pyNastran.bdf.cards.methods import EIGB, EIGC, EIGR, EIGP, EIGRL, MODTRAK
from pyNastran.bdf.cards.nodes import GRDSET # , GRID, SPOINTs, EPOINTs, POINT, SEQGP, GRIDB
from pyNastran.bdf.cards.aero.aero import (
    AECOMP, AECOMPL, AEFACT, AELINK, AELIST, AEPARM, AESURF, AESURFS,
    CAERO1, CAERO2, CAERO3, CAERO4, CAERO5,
    PAERO1, PAERO2, PAERO3, PAERO4, PAERO5,
    MONPNT1, MONPNT2, MONPNT3, MONDSP1,
    SPLINE1, SPLINE2, SPLINE3, SPLINE4, SPLINE5)
from pyNastran.bdf.cards.aero.static_loads import AESTAT, AEROS, CSSCHD, TRIM, TRIM2, DIVERG
from pyNastran.bdf.cards.aero.dynamic_loads import AERO, FLFACT, FLUTTER, GUST, MKAERO1, MKAERO2
from pyNastran.bdf.cards.optimization import (
    DCONADD, DCONSTR, DESVAR, TOPVAR, DDVAL, DOPTPRM, DLINK,
    DRESP1, DRESP2, DRESP3,
    DVCREL1, DVCREL2,
    DVMREL1, DVMREL2,
    DVPREL1, DVPREL2,
    DVGRID, DSCREEN)
#from pyNastran.bdf.cards.superelements import (
    #RELEASE, SEBNDRY, SEBULK, SECONCT, SEELT, SEEXCLD,
    #SELABEL, SELOAD, SELOC, SEMPLN, SENQSET, SETREE,
    #CSUPER, CSUPEXT,
#)
from pyNastran.bdf.cards.bdf_sets import (
    ASET, BSET, CSET, QSET, USET,
    ASET1, BSET1, CSET1, QSET1, USET1,
    OMIT, OMIT1,
    SET1, SET3,
    SEBSET, SECSET, SEQSET, # SEUSET
    SEBSET1, SECSET1, SEQSET1, # SEUSET1
    SESET, #SEQSEP
    RADSET,
)
from pyNastran.bdf.cards.params import PARAM, PARAM_MYSTRAN, PARAM_NASA95
from pyNastran.bdf.cards.dmig import DMIG, DMI, DMIJ, DMIK, DMIJI, DMIG_UACCEL, DTI, DMIAX
#from pyNastran.bdf.cards.thermal.loads import (QBDY1, QBDY2, QBDY3, QHBDY, TEMP, TEMPD, TEMPB3,
                                  #TEMPRB, QVOL, QVECT)
from pyNastran.bdf.cards.thermal.thermal import (CHBDYE, CHBDYG, CHBDYP, PCONV, PCONVM,
                                    PHBDY, CONV, CONVM, TEMPBC)
from pyNastran.bdf.cards.thermal.radiation import RADM, RADBC, RADCAV, RADLST, RADMTX, VIEW, VIEW3D
from pyNastran.bdf.cards.bdf_tables import (TABLED1, TABLED2, TABLED3, TABLED4,
                               TABLEM1, TABLEM2, TABLEM3, TABLEM4,
                               TABLES1, TABDMP1, TABLEST, TABLEHT, TABLEH1,
                               TABRND1, TABRNDG,
                               DTABLE)
from pyNastran.bdf.cards.contact import (
    BCRPARA, BCTADD, BCTSET, BSURF, BSURFS, BCPARA, BCTPARA, BCONP, BLSEG, BFRIC,
    BCTPARM, BGADD, BGSET, BCBODY)
#from pyNastran.bdf.cards.parametric.geometry import PSET, PVAL, FEEDGE, FEFACE, GMCURV, GMSURF

from pyNastran.bdf.case_control_deck import CaseControlDeck, Subcase
from pyNastran.bdf.bdf_methods import BDFMethods
from pyNastran.bdf.bdf_interface.get_card import GetCard
from pyNastran.bdf.bdf_interface.add_card import AddCards
from pyNastran.bdf.bdf_interface.bdf_card import BDFCard
from pyNastran.bdf.bdf_interface.write_mesh_file import WriteMeshs
from pyNastran.bdf.bdf_interface.uncross_reference import UnXrefMesh
from pyNastran.bdf.bdf_interface.verify_validate import verify_bdf, validate_bdf
from pyNastran.bdf.bdf_interface.stats import get_bdf_stats

from pyNastran.bdf.errors import (CrossReferenceError, DuplicateIDsError,
                                  CardParseSyntaxError, UnsupportedCard, DisabledCardError,
                                  SuperelementFlagError, ReplicationError)
from pyNastran.bdf.bdf_interface.pybdf import (
    BDFInputPy, _clean_comment, _clean_comment_bulk,
    _check_for_spaces, add_superelements_from_deck_lines)

#from .bdf_interface.add_card import CARD_MAP
if TYPE_CHECKING:  # pragma: no cover
    from cpylog import SimpleLogger

SOL_700 = {
    ## Explicit Nonlinear (SOL 700)
    'ABINFL', 'AIRBAG', 'ATBACC', 'ATBJNT', 'ATBSEG',
    'BARRIER',
    'BCBODY', 'BCBODY1', 'BCBOX', 'BCELIPS', 'BCGRID', 'BCMATL', 'BCONECT',
    'BCONPRG', 'BCONPRP', 'BCPROP', 'BCSEG', 'BCTABL1', 'BCTABLE',
    'BIAS', 'BJOIN', 'BSURF',
    'CDAMP1D', 'CDAMP2D', 'CELAS1D', 'CELAS2D', 'CMARKB2', 'CMARKN1', 'COHFRIC',
    'COMPUDS', 'CORD3R',
    'COUCOHF', 'COUOPT', 'COUP1FL', 'COUPINT', 'COUPLE',
    'CSPR', 'CYLINDR', 'DETSPH', 'DYFSISW', 'DYPARAM',
    'EOSDEF', 'EOSGAM', 'EOSGRUN', 'EOSIG', 'EOSJWL',
    'EOSMG', 'EOSNA', 'EOSPOL', 'EOSUDS',
    'EOSTAIT', 'EULFOR', 'EULFOR1', 'EULFREG',
    'FAILJC', 'FAILMPS', 'FAILUDS', 'FFCONTR',
    'FLOW', 'FLOWC', 'FLOWDEF', 'FLOWT', 'FLOWUDS', 'FORCE2', 'FORCUDS',
    'GBAG', 'GBAGCOU',
    'HEATLOS', 'HGSUPPR', 'HTRCONV', 'HTRRAD', 'HYDSTAT',
    'INFLCG', 'INFLFRC', 'INFLGAS', 'INFLHB', 'INFLTNK', 'INFLTR', 'INITGAS', 'LEAKAGE',
    'MATBV', 'MATDEUL', 'MATEP', 'MATF', 'MATFAB', 'MATHE', 'MATORT', 'MATRIG', 'MATVE',
    'MESH', 'MOMENT2', 'NLOUTUD',
    'PBEAML', 'PBELT', 'PCOMPA', 'PELAS1', 'PERMEAB', 'PERMGBG', 'PEULER', 'PEULER1', 'PMARKER', 'PMINC',
    'PORFCPL', 'PORFGBG', 'PORFLOW', 'PORFLWT', 'PORHOLE', 'PORHYDS', 'PORUDS',
    'PSHELL1', 'PSPRMAT', 'PVISC1', 'RBJOINT', 'RELEX',
    'SHREL', 'SHRPOL', 'SHRUDS',
    'SPHERE', 'SURFINI',
    'TABLUDS', 'TIC3', 'TICEL', 'TICEUDS', 'TICEUL1', 'TICREG', 'TICVAL',
    'TODYNA',
    'WALL',
    'YLDHY', 'YLDJC', 'YLDMC', 'YLDMSS', 'YLDPOL', 'YLDRPL',
    'YLDSG', 'YLDTM', 'YLDUDS', 'YLDVM', 'YLDZA',
}
MISSING_CARDS = {
    'CGEN', 'GMSPC', 'GMCURV', 'GMLOAD', 'FEFACE', 'GMSURF', 'GMINTS', 'PVAL', 'PINTS',
    'EGRID', 'ADAPT', 'GRIDG', 'MESHOPT', 'OUTPUT', 'OUTRCV', 'SPCG', 'GRIDU',
    'GMINTC', 'PINTC', 'GMBNDS', 'GMBC', 'GMCONV', 'PGEN', 'GMQVOL',
    'CNGRNT',

    # slot
    'GRIDF',
    'CSLOT3', 'CSLOT4', 'CAXIF2', 'CAXIF3', 'AXSLOT', 'SLBDY',

    'EQUIV', 'EXTRN', 'DSCONS', 'DVGEOM', 'DVAR', 'DVSET',
    'ADUM1', 'ADUM8', 'ADUM9', 'GRIDS',

    'CFLUID2', 'CFLUID3', 'CFLUID4', 'FSLIST', 'BNDGRID', 'BDYLIST', 'PRESPT',
    'FREEPT', 'FLSYM',
    # ----------------------------
    'RJOINT', 'RTRPLT', 'RTRPLT1', 'MDLPRM', 'DYNRED',

    ## fatigue
    'FTGDEF', 'FTGPARM', 'FTGEVNT', 'FTGLOAD', 'FTGSEQ',
    'MATFTG', 'PFTG', 'TABLFTG', 'UDNAME', 'SET4',
    'TOPSTR', 'TOPDMG',

    ## acoustic
    'CACINF3', 'CACINF4',

    ## Non Linear (SOL 400)
    'BCBODY', 'BSURF', 'BCTABLE', 'BCPARA', 'PSLDN1',


    ## Implicit Nonlinear (SOL 600) - marc
    'PARAMARC', 'MARCIN', 'MARCOUT',
    'RESTART',
    'MATD001', 'MATD003', 'MATD005','MATD006', 'MATD007', 'MATD009',
    'MATD012', 'MATD013', 'MATD014', 'MATD015', 'MATD018', 'MATD019',
    'MATD020', 'MATD022', 'MATD024', 'MATD026', 'MATD027', 'MATD028',
    'MATD030', 'MATD031', 'MATD034', 'MATD036',
    'MATD054', 'MATD055', 'MATD057', 'MATD059',
    'MATD062', 'MATD063', 'MATD064', 'MATD077',
    'MATD080', 'MATD081', 'MATD127', 'MATD181',
    'MATD20M',
    'MATD2AN', 'MATD2OR',
    'MATORT', 'MATTORT',
    'MATEP', 'MATTEP',
    'MATHE', 'MATTHE',
    'MATVP',
    'MATSMA',
    'MATVE', 'MATTVE',
    'MATG', 'MATTG'
    'MATF',
    'MATVB',
    'MATHED',
    'COHSEIV',
    'DEACTEL', 'ACTIVAT',
    # properties
    'PCOMPF', 'MSTACK', 'GASKET',
    # control
    'NLAUTO', 'NLDAMP', 'NLSTRAT',
    # solid -> shell
    'CSSHLH', 'CSSHLP',
    # solid element connector
    'CSSHL', 'PSSHL',
    # contact bodies
    'BCBODY', 'BCHANGE',
    'GMNURB', 'BSURF', 'BCBOX', 'BCPROP', 'BCMATL',
    # contact parameters
    'BCONTACT', 'BCPARA',
    # contact table
    'BCTABLE',
    # contact movement
    'BCMOVE',
    # thermal contact
    'MPHEAT', 'NLHEAT', 'MCHSTAT', 'MINSTAT',
    'MHEATSHL', 'MTHERM',

    ## bolts
    'MBOLT', 'MBOLTUS', 'BOLT', 'BOLTFRC', 'BOLTFOR',

    ## uds
    'PORUDS', 'YLDUDS', 'SHRUDS', 'FAILUDS', 'COMPUDS',
    'FLOWUDS','BCONUDS', 'ELEMUDS', 'UDSESV', 'MATUDS',
    'TICEUDS', 'EOSUDS', 'TABLUDS', 'GENUDS', 'ENTUDS',

    # explosives
    'EXPLSV', 'PLBLAST',

    # rotor
    'ROTOR', 'ROTORB', 'ROTPARM', 'ROTORAX', 'ROTORSE',

    # elements/properties
    'CELAS1D',
    'PRODN1',
    'CBARG',
    'PBEMN1',
    'PSHELL1', 'PSHL3D',
    'PSLDN1', 'PSHELLD', 'PSOLIDD',
    'PSHLN1', 'PSHLN2',
    'PCOMPG1', 'PCOMPLS', 'PLCOMP',
    'PBUSH2D', 'PBSH2DT',
    'CYSYM', 'CSEAM', 'PSEAM', 'CWSEAM', 'PWSEAM',
    'PACINF', 'PAXSYMH',
    'CBEAR', 'PBEAR',
    'CSPOT', 'CFILLET',
    'CBUTT',
    'PCOHE',
    'CWELD', 'PWELD',


    ## rigid_elements
    'RSPINT', 'RSPINR', 'RBE2GS',

    ## materials
    'MATF', 'MATEP', 'MATVE', 'MCOHE', 'MATM',
    'MATDT01', 'MATDIGI', 'MATUSR', 'MATTC',
    'MATORT', 'MATTORT', 'MATTHE', 'MATPLCY',
    'MATSMA', 'MAT8A', 'MATTEP',
    'MATPOR', 'MATDMG',
    'MAT2F', 'MAT8F',

    ## loads
    'FORCDST', 'PLOADX', 'PLOADB3', 'PLOADG', 'QBDY4', 'SLOADN1',
    'TEMPP1', 'TTEMP', 'RGYRO', 'PLOADE1', 'RCROSSC', 'LOADCYT',
    'TEMPN1',

    ## boundaries
    'BNDFREE', 'BNDFRE1',
    'BNDFIX', 'BNDFIX1',
    'SPCR',

    ## aero
    'UXVEC', 'GUST2', 'RVDOF', 'RVDOF1', 'AEFORCE', 'AEPRESS',

    ## acoustics
    'MICPNT',

    ## coords
    'CORD1RX', 'CORD3G',

    ## brakes
    'BRKSQL',

    ## d2r
    'D2R0000', 'D2RAUTO', 'D2RINER',

    'MATD016', 'MATD029', 'MATD053', 'MATD066', 'MATD067', 'MATD069',
    'MATD070', 'MATD078', 'MATD093', 'MATD094', 'MATD095', 'MATD097',
    'MATD098', 'MATD099',
    'MATDS01', 'MATDS02','MATDS03','MATDS04', 'MATDS05','MATDS06','MATDS07','MATDS08',
    'MATDS13','MATDS14','MATDS15',
    'MATDSW1', 'MATDSW2', 'MATDSW3', 'MATDSW4', 'MATDSW5',

    # ???
    'MONGRP', 'THPAD', 'TABLEDR',
    'DMIGOUT', 'MPCOUT', 'BCBZIER', 'BCPATCH', 'MDBULK', 'BCSCAP',
    'PRJCON',
    'TABLRPC', 'CIFQDX', 'NLADAPT', 'EOSTAB', 'EOSTABC',
    'TABLED5',
    'TIMNAT', 'ROTHYBD',
    'GRIA',
    'PANEL', 'TRMCPL',
    'DAMPING', 'DAMPGBL',
    'MONCARL',
    'PLCYISO', 'PLCYKIN', 'PLCYRUP',
    'BCTABLE', 'BCTABL1', 'BCAUTOP', 'BCBODY1', 'BCPROP', 'BCPROPS', 'BCBMRAD',
    'BGPARM', 'BCHANGE', 'BOUTPUT', 'TCNTPRM', 'BSQUEAL',
    'DELETE', 'RENAME',
    'ITER', 'GROUP', 'LIST',
    'MATRIG',
    'DVSHAP', 'DVBSHAP',
    'DISTORT', 'UNGLUE',
    'NLSTEP', 'NLAUTO', 'NLMOPTS', 'NLOUT', 'ERPPNL',
    'RBAX3D', 'ACPEMCP',
    'FTGDEF', 'CAMPBLL', 'UNBALNC', 'ECHO',
    'CINTC', 'GMBNDC',
    'SET4', 'PFTG', 'FTGPARM', 'UDNAME', 'FTGSEQ',
    'ELIST', 'MFLUID',
    'TEMPF', 'TEMPG', 'HADAPTL', 'HADACRI',
    'CIFPENT',
    'ELAR2', 'EBDSET', 'TMCPARA',
    'PEULER1', 'MATDEUL', 'EOSPOL', 'SHREL', 'YLDMC', 'PMINC', 'COUPLE', 'COUCOHF',
    'COHFRIC', 'MESH', 'TICVAL', 'TICEUL1', 'TICREG', 'CYLINDR',
    'SWLDPRM', 'PLOTOPT', 'PLOTE', 'PLOTG',
    'FTGEVNT', 'RADBND', 'RADMT', 'NOLIN1', 'NOLIN2', 'NOLIN3', 'NOLIN4',
    'NLFREQ1', 'NLRGAP',
    'ACADAPT', 'ACORDER', 'AMLREG'
    'NLCNTL', 'NLCNTLG', 'NLCNTL2', 'NLSTRAT', 'NLRSFD', 'NLHARM',
    'DMRLAW',
    'TABLE3D', 'TABL3D0',
    'CONTRLT',
    'MPOINT',
    'CFTUBE', 'CHBDY', 'MONSUM', 'TMPSET', 'HYBDAMP',
    'TOMVAR', 'STOCHAS', 'NHRMPRM', 'MONSUM', 'CYLINDER', 'EOSGAM',
    'SPHERE', 'BJOIN', 'PMARKER', 'TICEUL', 'PEULER', 'BCBOX',
    'NLOUTUD', 'YLDVM', 'MATBV', 'SEDRSP2', 'SEDRSP3', 'SEDLINK', 'MFUN',
    'MTAB', 'TIC3', 'TABSCTL', 'SPLINEX', 'BCONECT', 'BCONPRG',
    'BCSURF', 'BCGRID', 'PAXISYM', 'BEADVAR', 'BLDOUT', 'GRIDMOD',
    'IPSTRAIN', 'RESTART', 'TICD',
    'RCROSS', 'CYSUP', 'FBODYLD', 'BDYOR', 'RANDVAR',
    'L16MOD', 'SECTAX',
    'ACIFPRM', 'BCBDPRP', 'MDMIOUT', 'MPROCS', 'FBODYSB',
    'CCRSFIL', 'COMBWLD',
    'CIFHEX',
    'FSICTRL', 'FLOW', 'FLOWDEF', 'BCSEG', 'CMARKN1', 'LEAKAGE', 'TICEL',
    'EOSMG', 'YLDJC', 'MAT1A', 'HGSUPPR',
    'EOSTAIT', 'EOSJWL', 'SURFINI', 'WALL', 'MPCY',
    'TIMNVH', 'CIFQUAD', 'VCCT',
    'DYTIMHS', 'PRESTRS', 'SEQROUT', 'BCRIGID', 'BCMOVE', 'ISTRESS',
    'WETLOAD', 'WETELMG', 'COUP1FL', 'PORFLOW', 'PERMEAB', 'ENDDYNA',
    'COUOPT', 'COUPLE1', 'COUP1INT', 'FAILJC',
    'DETSPH', 'BIAS', 'PORFCPL', 'RELAX', 'FLOW', 'ISTRSSH',
    'NTHICK', 'TIM2PSD', 'WETSURF', 'WETELME', 'BCNURBS',
    'BCTRIM', 'IMPGEOM', 'IMPCASE', 'SPCD2', 'SPRBCK', 'BCRGSRF', 'BCNURB2',
    'CAXISYM', 'RADC', 'VIEWEX',
    'ACLOAD', 'PBARN1'
    , 'TABL3D1', 'AEGRID', 'AEQUAD4', 'SPBLND1', 'CONCTL',
    'MAT1F', 'HYDROS', 'HYDROC', 'DVLREL1', 'DTABLE2', 'DVPSURF', 'PFASTT',
    'FRFRELS', 'FRFCONN', 'FRFXIT1', 'FBALOAD', 'MATS8', 'METADATA',
    'PRIM1', 'PRIM7', 'CONV3', 'GRIDA', 'MAT10F', 'RADCOL', 'SPLINRB',
    'TABL3D2', 'DYMAT24', 'FBAPHAS', 'FRFFLEX', 'FBADLAY', 'FRFXIT',
    'FRFCOMP', 'FRFSPC1', 'PCOMPFQ', 'PDISTB', 'MASSSET', 'PSLDN2', 'MATSORT',
    'NLHEATC', 'MDRBE2', 'MDRBE3', 'MDEXCLD', 'MDWELD', 'MDBCNCT', 'MDCONCT',
    'MDBCTB1', 'MDMPC', 'MDRROD', 'ACCSSPT', 'MDLOC', 'MDMPLN', 'MDMOVE', 'MDWELD',
    'DEFUSET', 'MDFAST', 'MDBOLT', 'MDSEAM', 'ALIASM', 'POSTBUK', 'MATDB01',
    'PBEAM71', 'PSHEARN', 'MATD010', 'PBDISCR', 'PBELTD', 'CORD3RX', 'RBE2A',
    'ACCMETR', 'CBELT', 'RBJSTIF', 'MDRJNT', 'BCPFLG', 'MDLABEL', 'MDBNDRY',
    'CHEXCZ', 'CPENTCZ', 'BEDGE', 'DVEREL1', 'DMNCON', 'DVTREL1', 'NLCNTL',
    'MATCRP', 'TLOAD3', 'MATSR', 'DTEMP',
    'MDDMIG', 'MDTRAN', 'MDROT1', 'MDROT2', 'MDMIR1', 'MDMIR2', 'MDMIAUX',
    'MPCREEP',
    'COHESIV', 'CSSHLM', 'PBMARB6', 'PBMNUM6', 'DMIGROT', 'GRNDSPR', 'SUPORT6',
    'MARPRN', 'MATNLE2', 'MATNLE3', 'MATNLE4', 'MATNLE5', 'MATNLE6', 'MATPDR',
    'MGRSPR', 'MIXTURE', 'MNF600', 'MT16SEL', 'MTABRV', 'NLBSH3D', 'MONCNCM',
    'FREQV', 'VATVFS', 'PMIC', 'MAT10C', 'ALOAD', 'ELAR', 'ATVBULK', 'AMLREG',
    'ATVFS', 'BOLTLD', 'BCTPAR2', 'MATFT', 'PLOTEL4', 'CYCADD', 'MATT11'
}


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
    def __init__(self, debug: str | bool | None,
                 log: Optional[SimpleLogger]=None,
                 mode: str='msc') -> None:
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
            valid_modes = {'msc', 'nx', 'mystran', 'nasa95', 'zona'}

        """
        assert debug in [True, False, None], f'debug={debug!r}'
        self.echo = False
        self.read_includes = True

        # file management parameters
        self.active_filenames: list[str] = []
        self.active_filename: Optional[str] = None
        self.include_dir = ''
        self.dumplines = False

        log_args = {} if CPYLOG_VERSION <= '1.5.0' else {'nlevels': 3}
        self.log = get_logger2(log=log, debug=debug, **log_args)

        # list of all read in cards - useful in determining if entire BDF
        # was read & really useful in debugging
        self.card_count: dict[str, int] = {}

        # stores the card_count of cards that have been rejected
        self.reject_count: dict[str, int] = {}

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

        # lines that were rejected b/c they were for a card that isn't supported
        self.reject_lines: list[list[str]] = []

        # cards that were created, but not processed
        self.reject_cards: list[str] = []

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
            'CTRSHL', 'CQUAD1',
            'SNORM',

            'CPLSTN3', 'CPLSTN4', 'CPLSTN6', 'CPLSTN8', # plate strain
            'CPLSTS3', 'CPLSTS4', 'CPLSTS6', 'CPLSTS8', # plate stress

            # acoustic
            'CHACAB', 'CAABSF', 'CHACBR',
            'PACABS', 'PAABSF', 'PACBAR', 'ACMODL',

            'CTETRA', 'CPYRAM', 'CPENTA', 'CHEXA',
            'CIHEX1', 'CIHEX2', 'CHEXA1', 'CHEXA2',
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

            'PSHELL', 'PCOMP', 'PCOMPG', 'PSHEAR', 'PTRSHL', 'PQUAD1',
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
            'MATG', 'MATHE', 'MATHP', 'MATEV',

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
            'SPC', 'SPCADD', 'SPC1', 'SPCOFF', 'SPCOFF1',  # 'SPCAX', removed
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
            'BCBODY',  ## bcbody
            'BCPARA',  ## bcpara
            'BCTPARA',  ## bctpara
            'BCRPARA',  ## bcrpara
            'BCTPARM', ## bctparm
            'BGADD',  ## bgadds
            'BGSET',  ## bgsets
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
        map_version(self, mode)

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

    #def get_h5attrs(self) -> list[str]:
        #"""helper method for dict_to_h5py"""
        #attrs = self.object_attributes(mode='both', keys_to_skip=None)
        #return attrs

    #def export_hdf5_filename(self, hdf5_filename: str) -> None:
        #"""
        #Converts the BDF objects into hdf5 object

        #Parameters
        #----------
        #hdf5_filename : str
            #the path to the hdf5 file

        #TODO: doesn't support:
          #- BucklingEigenvalues

        #"""
        #import h5py
        #with h5py.File(hdf5_filename, 'w') as hdf5_file:
            #self.export_hdf5_file(hdf5_file)

    #def export_hdf5_file(self, hdf5_file, exporter=None) -> None:
        #"""
        #Converts the BDF objects into hdf5 object

        #Parameters
        #----------
        #hdf5_file : H5File()
            #an h5py object
        #exporter : HDF5Exporter; default=None
            #unused

        #"""
        #from pyNastran.bdf.bdf_interface.hdf5_exporter import export_bdf_to_hdf5_file
        #export_bdf_to_hdf5_file(hdf5_file, self)

    #def load_hdf5_filename(self, hdf5_filename: str) -> None:
        #"""
        #Loads a BDF object from an hdf5 filename

        #Parameters
        #----------
        #hdf5_filename : str
            #the path to the hdf5 file

        #"""
        #import h5py
        #with h5py.File(hdf5_filename, 'r') as hdf5_file:
            #self.load_hdf5_file(hdf5_file)

    #def load_hdf5_file(self, h5_file) -> None:
        #"""
        #Loads a BDF object from an hdf5 object

        #Parameters
        #----------
        #hdf5_file : H5File()
            #an h5py object
        #exporter : HDF5Exporter; default=None
            #unused

        #"""
        #from pyNastran.bdf.bdf_interface.hdf5_loader import load_bdf_from_hdf5_file
        #load_bdf_from_hdf5_file(h5_file, self)

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
            'type_slot_str',
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
        cards : list[str]; set[str]
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

    def enable_cards(self, cards: list[str] | set[str]) -> None:
        """
        Method for setting the cards that will be processed

        Parameters
        ----------
        cards : list[str]; set[str]
            a list/set of cards that should not be read

        .. python ::

            bdf_model.enable_cards(['GRID', 'CTRIA3'])

        """
        if cards is None:
            return
        elif isinstance(cards, str):
            enable_set = set([cards])
        else:
            enable_set = set(cards)
        self.cards_to_read = enable_set

    def set_error_storage(self, nparse_errors: int=100,
                          stop_on_parsing_error: bool=True,
                          nxref_errors: int=100,
                          stop_on_xref_error: bool=True) -> None:
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
        self.xref_obj.set_error_storage(
            nparse_errors=nparse_errors,
            stop_on_parsing_error=stop_on_parsing_error,
            nxref_errors=nxref_errors,
            stop_on_xref_error=stop_on_xref_error
        )
        assert isinstance(nparse_errors, int), type(nparse_errors)
        #assert isinstance(nxref_errors, int), type(nxref_errors)
        self._nparse_errors = nparse_errors
        #self._nxref_errors = nxref_errors
        self._stop_on_parsing_error = stop_on_parsing_error
        #self._stop_on_xref_error = stop_on_xref_error

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
                    make_ilines: bool=True) -> (list[str], Any):
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
        all_lines : list[str]
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
        if not make_ilines:
            ilines = None
        self._set_pybdf_attributes(obj, save_file_structure=False)
        return all_lines, ilines

    def _set_pybdf_attributes(self, obj: BDFInputPy,
                              save_file_structure: bool=False) -> None:
        """common method for all functions that use BDFInputPy"""
        # these are include line pairs
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
        self.log.debug(f'---starting BDF.read_bdf of {self.bdf_filename}---')
        self._parse_primary_file_header(bdf_filename)

        obj = BDFInputPy(self.read_includes, self.dumplines, self._encoding,
                         nastran_format=self.nastran_format,
                         consider_superelements=self.is_superelements,
                         log=self.log, debug=self.debug)
        obj.use_new_parser = self.use_new_deck_parser

        out = obj.get_lines(bdf_filename, punch=self.punch, make_ilines=True)
        (system_lines,
         executive_control_lines,
         case_control_lines,
         bulk_data_lines, bulk_data_ilines,
         additional_deck_lines, additional_deck_lines) = out
        self._set_pybdf_attributes(obj, save_file_structure)

        self.system_command_lines = system_lines
        self.executive_control_lines = executive_control_lines
        self.case_control_lines = case_control_lines

        sol, method, sol_iline, app = parse_executive_control_deck(executive_control_lines)
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

        if additional_deck_lines:
            add_superelements_from_deck_lines(self, BDF, additional_deck_lines)

        self.pop_parse_errors()
        fill_dmigs(self)

        if validate:
            self.validate()

        self.cross_reference(xref=xref)
        self._xref = xref

        self.log.debug('---finished BDF.read_bdf of %s---' % self.bdf_filename)
        if self.reject_cards:
            rejected_cards = list(set([card[0] for card in self.reject_cards]))
            rejected_cards.sort()
            self.log.warning(f'rejecting {rejected_cards}')
        #for card_type in self.reject_cards:

    def _parse_all_cards(self, bulk_data_lines: list[str], bulk_data_ilines: Any) -> None:
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
                    msg = f'{key!r} is an invalid type; only dictionaries are supported'
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
        xref_obj = self.xref_obj
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

            if xref_obj._stop_on_xref_error:
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
        self.xref_obj.pop_xref_errors()

    def get_bdf_cards(self, bulk_data_lines: list[str],
                      bulk_data_ilines: Optional[Any]=None) -> tuple[Any, Any, Any]:
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
                        msg = f'exiting due to ENDDATA found with {nleftover:d} lines left'
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
            self.log.debug(f'sol={self.sol} method={self.sol_method!r}')
        else:  # very common
            self.sol_method = None

    def update_card(self, card_name: str, icard: int, ifield: int,
                    value: int | float | str) -> None:
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
            assert isinstance(key, str), f'key={key!r}'
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

    def set_dynamic_syntax(self, dict_of_vars: dict[str, int | float | str]) -> None:
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
        assert len(dict_of_vars) > 0, f'nvars = {len(dict_of_vars):d}'
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

    def _process_card(self, card_lines: list[str]) -> list[str]:
        """
        Converts card_lines into a card.
        Considers dynamic syntax and removes empty fields

        Parameters
        ----------
        card_lines : list[str]
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

    def _parse_dynamic_syntax(self, key: str) -> dict[str, Any]:
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
        add_methods = self._add_methods
        self._card_parser = {
            ##'=' : (Crash, None),
            #'/' : (Crash, None),

            ##'CGEN' : (CrashIgnore, None),

            #'SETREE' : (SETREE, self._add_setree_object),
            #'SENQSET' : (SENQSET, self._add_senqset_object),
            #'SEBULK' : (SEBULK, self._add_sebulk_object),
            #'RELEASE': (RELEASE, self._add_release_object),
            #'SEBNDRY' : (SEBNDRY, self._add_sebndry_object),
            #'SEELT' : (SEELT, self._add_seelt_object),
            #'SELOC' : (SELOC, self._add_seloc_object),
            #'SEMPLN' : (SEMPLN, self._add_sempln_object),
            #'SECONCT' : (SECONCT, self._add_seconct_object),
            #'SELABEL' : (SELABEL, self._add_selabel_object),
            #'SEEXCLD' : (SEEXCLD, self._add_seexcld_object),
            #'CSUPER' : (CSUPER, self._add_csuper_object),
            #'CSUPEXT' : (CSUPEXT, self._add_csupext_object),
            #'SELOAD' : (SELOAD, self._add_seload_object),

            #'CHACAB': (CHACAB, self._add_element_object),
            #'CHACBR': (CHACBR, self._add_element_object),
            #'CAABSF': (CAABSF, self._add_element_object),
            #'PACABS': (PACABS, self._add_acoustic_property_object),
            #'PAABSF': (PAABSF, self._add_acoustic_property_object),
            #'PACBAR': (PACBAR, self._add_acoustic_property_object),
            ##'PANEL' : (Crash, None),

            #'BCONP' : (BCONP, self._add_bconp_object),
            #'BLSEG' : (BLSEG, self._add_blseg_object),
            #'BFRIC' : (BFRIC, self._add_bfric_object),
            #'MODTRAK' : (MODTRAK, self._add_modtrak_object),

            ##  nx contact
            #'BCPARA' : (BCPARA, self._add_bcpara_object),
            #'BCTPARM' : (BCTPARM, self._add_bctparam_object),
            #'BGADD' : (BGADD, self._add_bgadd_object),
            #'BGSET' : (BGSET, self._add_bgset_object),
            #'BCBODY' : (BCBODY, self._add_bcbody_object),

            ## 'BOLT', 'BOLTFOR'
            #'BOLT' : (Crash, None),
            #'BOLTFOR' : (Crash, None),

            ##'CBEAR', 'PBEAR', 'ROTORB',
            #'CBEAR' : (Crash, None),
            #'PBEAR' : (Crash, None),
            #'ROTORB' : (Crash, None),

            ##'SWLDPRM' : (Crash, None),

            ##'CWELD' : (Crash, None),
            ##'PWELD' : (Crash, None),
            ##'PWSEAM' : (Crash, None),
            ##'CWSEAM' : (Crash, None),
            ##'CSEAM' : (Crash, None),
            ##'PSEAM' : (Crash, None),

            ##'DVSHAP' : (Crash, None),
            ##'BNDGRID' : (Crash, None),

            ##'CYSYM' : (Crash, None),
            ##'TEMPP1' : (Crash, None),
            ##'DSCONS' : (Crash, None),
            ##'DVAR' : (Crash, None),
            ##'DVSET' : (Crash, None),
            ##'DYNRED' : (Crash, None),
            ##'BNDFIX' : (Crash, None),
            ##'BNDFIX1' : (Crash, None),

            ##'AEFORCE' : (Crash, None),
            ##'UXVEC' : (Crash, None),
            #'GUST2' : (Crash, None),

            ##'RADBND' : (Crash, None),


            ## nodes
            #'GRID' : (GRID, self._add_node_object),
            #'SPOINT' : (SPOINTs, self._add_spoint_object),
            #'EPOINT' : (EPOINTs, self._add_epoint_object),
            #'RINGAX' : (RINGAX, self._add_ringax_object),
            #'POINTAX' : (POINTAX, self._add_ringax_object),
            #'POINT' : (POINT, self._add_point_object),
            #'SEQGP' : (SEQGP, self._add_seqgp_object),
            #'GRIDB' : (GRIDB, self._add_gridb_object),

            'PARAM' : (PARAM, add_methods._add_param_object),

            'CORD2R' : (CORD2R, add_methods._add_coord_object),
            'CORD2C' : (CORD2C, add_methods._add_coord_object),
            'CORD2S' : (CORD2S, add_methods._add_coord_object),

            ## parametric
            #'PSET' : (PSET, self._add_pset),
            #'PVAL' : (PVAL, self._add_pval),
            #'GMCURV' : (GMCURV, self._add_gmcurv),
            #'GMSURF' : (GMSURF, self._add_gmsurf),
            #'FEFACE' : (FEFACE, self._add_feface),
            #'FEEDGE' : (FEEDGE, self._add_feedge),

            ## msgmesh
            #'GMCORD' : (GMCORD, self._add_coord_object),
            #'CGEN' : (CGEN, self._add_element_object),

            #'CONROD' : (CONROD, self._add_element_object),
            #'CROD' : (CROD, self._add_element_object),
            'PROD' : (PROD, add_methods._add_property_object),
            #'CTUBE' : (CTUBE, self._add_element_object),
            'PTUBE' : (PTUBE, add_methods._add_property_object),

            #'BAROR' : (BAROR, self._add_baror_object),
            #'CBARAO' : (CBARAO, self._add_ao_object),
            'PBAR' : (PBAR, add_methods._add_property_object),
            'PBARL' : (PBARL, add_methods._add_property_object),
            #'PBRSECT' : (PBRSECT, self._add_property_object),

            #'BEAMOR' : (BEAMOR, self._add_beamor_object),
            'PBEAM' : (PBEAM, add_methods._add_property_object),
            'PBEAML' : (PBEAML, add_methods._add_property_object),
            #'PBCOMP' : (PBCOMP, self._add_property_object),
            #'PBMSECT' : (PBMSECT, self._add_property_object),

            #'CBEAM3' : (CBEAM3, self._add_element_object),
            #'PBEAM3' : (PBEAM3, self._add_property_object),

            #'CBEND' : (CBEND, self._add_element_object),
            'PBEND' : (PBEND, add_methods._add_property_object),

            # ------------------------------------------------------------------

            #'CTRSHL' : (CTRSHL, self._add_element_object),  # nasa95
            #'CTRIA3' : (CTRIA3, self._add_element_object),
            #'CQUAD1' : (CQUAD1, self._add_element_object),  # nasa95
            #'CQUAD4' : (CQUAD4, self._add_element_object),
            #'CQUAD' : (CQUAD, self._add_element_object),
            #'CQUAD8' : (CQUAD8, self._add_element_object),
            #'CQUADX' : (CQUADX, self._add_element_object),
            #'CQUADX4' : (CQUADX4, self._add_element_object),
            #'CQUADX8' : (CQUADX8, self._add_element_object),
            #'CQUADR' : (CQUADR, self._add_element_object),
            #'CTRIA6' : (CTRIA6, self._add_element_object),
            #'CTRIAR' : (CTRIAR, self._add_element_object),
            #'CTRAX3' : (CTRAX3, self._add_element_object),
            #'CTRAX6' : (CTRAX6, self._add_element_object),
            #'CTRIAX' : (CTRIAX, self._add_element_object),
            #'CTRIAX6' : (CTRIAX6, self._add_element_object),
            #'SNORM' : (SNORM, self._add_normal_object),
            'PCOMP' : (PCOMP, add_methods._add_property_object),
            'PCOMPG' : (PCOMPG, add_methods._add_property_object),
            'PSHELL' : (PSHELL, add_methods._add_property_object),
            #'PTRSHL' : (PTRSHL, self._add_property_object),
            #'PQUAD1' : (PQUAD1, self._add_property_object),
            'PLPLANE' : (PLPLANE, add_methods._add_property_object),
            #'CPLSTN3' : (CPLSTN3, self._add_element_object),
            #'CPLSTN4' : (CPLSTN4, self._add_element_object),
            #'CPLSTN6' : (CPLSTN6, self._add_element_object),
            #'CPLSTN8' : (CPLSTN8, self._add_element_object),
            #'CPLSTS3' : (CPLSTS3, self._add_element_object),
            #'CPLSTS4' : (CPLSTS4, self._add_element_object),
            #'CPLSTS6' : (CPLSTS6, self._add_element_object),
            #'CPLSTS8' : (CPLSTS8, self._add_element_object),
            #'PPLANE' : (PPLANE, self._add_property_object),

            #'CSHEAR' : (CSHEAR, self._add_element_object),
            'PSHEAR' : (PSHEAR, add_methods._add_property_object),

            ## nastran95
            #'CIHEX1' : (CIHEX1, self._add_element_object),
            #'CIHEX2' : (CIHEX2, self._add_element_object),
            #'CHEXA1' : (CHEXA1, self._add_element_object),
            #'CHEXA2' : (CHEXA2, self._add_element_object),
            #'PIHEX' : (PIHEX, self._add_property_object),

            ## msc/nx
            'PSOLID' : (PSOLID, add_methods._add_property_object),
            'PLSOLID' : (PLSOLID, add_methods._add_property_object),
            #'PCOMPS' : (PCOMPS, self._add_property_object),

            #'CELAS1' : (CELAS1, self._add_element_object),
            #'CELAS2' : (CELAS2, self._add_element_object),
            #'CELAS3' : (CELAS3, self._add_element_object),
            #'CELAS4' : (CELAS4, self._add_element_object),
            #'CVISC' : (CVISC, self._add_element_object),
            #'PELAST' : (PELAST, self._add_pelast_object),

            #'CDAMP1' : (CDAMP1, self._add_damper_object),
            #'CDAMP2' : (CDAMP2, self._add_damper_object),
            #'CDAMP3' : (CDAMP3, self._add_damper_object),
            ## CDAMP4 added later because the documentation is wrong
            #'CDAMP5' : (CDAMP5, self._add_damper_object),
            #'PDAMP5' : (PDAMP5, self._add_property_object),

            #'CFAST' : (CFAST, self._add_damper_object),
            #'PFAST' : (PFAST, self._add_property_object),

            #'CGAP' : (CGAP, self._add_element_object),
            #'PGAP' : (PGAP, self._add_property_object),

            #'CBUSH' : (CBUSH, self._add_damper_object),
            #'CBUSH1D' : (CBUSH1D, self._add_damper_object),
            #'CBUSH2D' : (CBUSH2D, self._add_damper_object),
            'PBUSH' : (PBUSH, add_methods._add_property_object),
            'PBUSH1D' : (PBUSH1D, add_methods._add_property_object),

            #'CRAC2D' : (CRAC2D, self._add_element_object),
            #'PRAC2D' : (PRAC2D, self._add_property_object),

            #'CRAC3D' : (CRAC3D, self._add_element_object),
            #'PRAC3D' : (PRAC3D, self._add_property_object),

            #'GENEL' : (GENEL, self._add_element_object),


            #'PDAMPT' : (PDAMPT, self._add_pdampt_object),
            #'PBUSHT' : (PBUSHT, self._add_pbusht_object),

            #'CCONEAX' : (CCONEAX, self._add_element_object),
            #'PCONEAX' : (PCONEAX, self._add_property_object),
            #'AXIC' : (AXIC, self._add_axic_object),
            #'AXIF' : (AXIF, self._add_axif_object),
            #'CYAX' : (CYAX, self._add_cyax_object),

            'RBAR' : (RBAR, add_methods._add_rigid_element_object),
            'RBAR1' : (RBAR1, add_methods._add_rigid_element_object),
            'RBE1' : (RBE1, add_methods._add_rigid_element_object),
            'RBE2' : (RBE2, add_methods._add_rigid_element_object),
            'RBE3' : (RBE3, add_methods._add_rigid_element_object),
            'RROD' : (RROD, add_methods._add_rigid_element_object),
            'RSPLINE' : (RSPLINE, add_methods._add_rigid_element_object),
            'RSSCON' : (RSSCON, add_methods._add_rigid_element_object),


            ## there is no MAT6 or MAT7
            'MAT1' : (MAT1, add_methods._add_structural_material_object),
            'MAT2' : (MAT2, add_methods._add_structural_material_object),
            'MAT3' : (MAT3, add_methods._add_structural_material_object),
            'MAT8' : (MAT8, add_methods._add_structural_material_object),
            'MAT9' : (MAT9, add_methods._add_structural_material_object),
            'MAT10' : (MAT10, add_methods._add_structural_material_object),
            'MAT11' : (MAT11, add_methods._add_structural_material_object),
            'MAT3D' : (MAT3D, add_methods._add_structural_material_object),
            'EQUIV' : (EQUIV, add_methods._add_structural_material_object),
            'MATG' : (MATG, add_methods._add_structural_material_object),

            'MATHE' : (MATHE, add_methods._add_hyperelastic_material_object),
            'MATHP' : (MATHP, add_methods._add_hyperelastic_material_object),
            'MATEV' : (MATEV, add_methods._add_structural_material_object),
            'MAT4' : (MAT4, add_methods._add_thermal_material_object),
            'MAT5' : (MAT5, add_methods._add_thermal_material_object),

            'MATS1' : (MATS1, add_methods._add_material_dependence_object),
            #'MATS3' : (MATS3, add_methods._add_material_dependence_object),
            #'MATS8' : (MATS8, add_methods._add_material_dependence_object),
            'MATT1' : (MATT1, add_methods._add_material_dependence_object),
            'MATT2' : (MATT2, add_methods._add_material_dependence_object),
            'MATT3' : (MATT3, add_methods._add_material_dependence_object),
            'MATT4' : (MATT4, add_methods._add_material_dependence_object),
            'MATT5' : (MATT5, add_methods._add_material_dependence_object),
            'MATT8' : (MATT8, add_methods._add_material_dependence_object),
            'MATT9' : (MATT9, add_methods._add_material_dependence_object),
            'NXSTRAT' : (NXSTRAT, add_methods._add_nxstrat_object),

            ## hasn't been verified, links up to MAT1, MAT2, MAT9 w/ same MID
            'CREEP' : (CREEP, add_methods._add_creep_material_object),

            #'NSMADD' : (NSMADD, self._add_nsmadd_object),
            #'NSM1' : (NSM1, self._add_nsm_object),
            #'NSML1' : (NSML1, self._add_nsm_object),

            #'CONM1' : (CONM1, self._add_mass_object),
            #'CONM2' : (CONM2, self._add_mass_object),
            #'CMASS1' : (CMASS1, self._add_mass_object),
            #'CMASS2' : (CMASS2, self._add_mass_object),
            #'CMASS3' : (CMASS3, self._add_mass_object),
            ## CMASS4 - added later because documentation is wrong

            'MPC' : (MPC, add_methods._add_constraint_mpc_object),
            'MPCADD' : (MPCADD, add_methods._add_constraint_mpcadd_object),

            'SPC' : (SPC, add_methods._add_constraint_spc_object),
            'SPC1' : (SPC1, add_methods._add_constraint_spc_object),
            'SPCOFF' : (SPCOFF, add_methods._add_constraint_spcoff_object),
            'SPCOFF1' : (SPCOFF1, add_methods._add_constraint_spcoff_object),
            'SPCAX' : (SPCAX, add_methods._add_constraint_spc_object),
            'SPCADD' : (SPCADD, add_methods._add_constraint_spcadd_object),
            #'GMSPC' : (GMSPC, add_methods._add_constraint_spc_object),

            'SESUP' : (SESUP, add_methods._add_sesuport_object), # pseudo-constraint
            'SUPORT' : (SUPORT, add_methods._add_suport_object), # pseudo-constraint
            'SUPORT1' : (SUPORT1, add_methods._add_suport1_object),  # pseudo-constraint

            #'LSEQ' : (LSEQ, add_methods._add_lseq_object),
            #'LOAD' : (LOAD, add_methods._add_load_combination_object),
            'CLOAD' : (CLOAD, add_methods._add_load_combination_object),
            #'LOADCYN' : (LOADCYN, add_methods._add_load_object),
            'LOADCYH' : (LOADCYH, add_methods._add_load_object),

            ## basic static loads
            #'GRAV' : (GRAV, self._add_load_object),
            #'ACCEL' : (ACCEL, self._add_load_object),
            #'ACCEL1' : (ACCEL1, self._add_load_object),
            #'RFORCE' : (RFORCE, self._add_load_object),
            #'RFORCE1' : (RFORCE1, self._add_load_object),
            #'SLOAD' : (SLOAD, self._add_load_object),
            #'GMLOAD' : (GMLOAD, self._add_load_object),
            #'SPCD' : (SPCD, self._add_load_object),  # enforced displacement
            #'QVOL' : (QVOL, self._add_load_object),  # thermal

            ## axisymmetric loads
            #'FORCEAX' : (FORCEAX, self._add_load_object),
            #'PLOADX1' : (PLOADX1, self._add_load_object),
            #'PRESAX' : (PRESAX, self._add_load_object),  # axisymmetric

            #'DLOAD' : (DLOAD, self._add_dload_object),

            'ACSRCE' : (ACSRCE, add_methods._add_dload_entry),
            'TLOAD1' : (TLOAD1, add_methods._add_dload_entry),
            'TLOAD2' : (TLOAD2, add_methods._add_dload_entry),
            'RLOAD1' : (RLOAD1, add_methods._add_dload_entry),
            'RLOAD2' : (RLOAD2, add_methods._add_dload_entry),
            'RANDPS' : (RANDPS, add_methods._add_dload_entry), # random
            'RANDT1' : (RANDT1, add_methods._add_dload_entry), # random
            #'QVECT' : (QVECT, add_methods._add_dload_entry),

            'FREQ' : (FREQ, add_methods._add_freq_object),
            'FREQ1' : (FREQ1, add_methods._add_freq_object),
            'FREQ2' : (FREQ2, add_methods._add_freq_object),
            'FREQ3' : (FREQ3, add_methods._add_freq_object),
            'FREQ4' : (FREQ4, add_methods._add_freq_object),
            'FREQ5' : (FREQ5, add_methods._add_freq_object),

            'DOPTPRM' : (DOPTPRM, add_methods._add_doptprm_object),
            'DESVAR' : (DESVAR, add_methods._add_desvar_object),
            'TOPVAR' : (TOPVAR, add_methods._add_topvar_object),
            ## BCTSET

            #'TEMPRB' : (TEMPRB, self._add_thermal_load_object),
            #'TEMP' : (TEMP, self._add_thermal_load_object),
            #'TEMPB3' : (TEMPB3, self._add_thermal_load_object),
            #'QBDY1' : (QBDY1, self._add_thermal_load_object),
            #'QBDY2' : (QBDY2, self._add_thermal_load_object),
            #'QBDY3' : (QBDY3, self._add_thermal_load_object),
            #'QHBDY' : (QHBDY, self._add_thermal_load_object),
            #'PHBDY' : (PHBDY, self._add_phbdy_object),

            #'CHBDYE' : (CHBDYE, self._add_thermal_element_object),
            #'CHBDYG' : (CHBDYG, self._add_thermal_element_object),
            #'CHBDYP' : (CHBDYP, self._add_thermal_element_object),
            #'PCONV' : (PCONV, self._add_convection_property_object),
            #'PCONVM' : (PCONVM, self._add_convection_property_object),

            #'VIEW' : (VIEW, self._add_view_object),
            #'VIEW3D' : (VIEW3D, self._add_view3d_object),

            ## aero
            'AECOMP' : (AECOMP, add_methods._add_aecomp_object),
            'AECOMPL' : (AECOMPL, add_methods._add_aecomp_object),
            'AEFACT' : (AEFACT, add_methods._add_aefact_object),
            'AELINK' : (AELINK, add_methods._add_aelink_object),
            'AELIST' : (AELIST, add_methods._add_aelist_object),
            'AEPARM' : (AEPARM, add_methods._add_aeparm_object),
            'AESTAT' : (AESTAT, add_methods._add_aestat_object),
            'AESURF' : (AESURF, add_methods._add_aesurf_object),
            'AESURFS' : (AESURFS, add_methods._add_aesurfs_object),

            'CAERO1' : (CAERO1, add_methods._add_caero_object),
            'CAERO2' : (CAERO2, add_methods._add_caero_object),
            'CAERO3' : (CAERO3, add_methods._add_caero_object),
            'CAERO4' : (CAERO4, add_methods._add_caero_object),
            'CAERO5' : (CAERO5, add_methods._add_caero_object),

            'PAERO1' : (PAERO1, add_methods._add_paero_object),
            'PAERO2' : (PAERO2, add_methods._add_paero_object),
            'PAERO3' : (PAERO3, add_methods._add_paero_object),
            'PAERO4' : (PAERO4, add_methods._add_paero_object),
            'PAERO5' : (PAERO5, add_methods._add_paero_object),

            'SPLINE1' : (SPLINE1, add_methods._add_spline_object),
            'SPLINE2' : (SPLINE2, add_methods._add_spline_object),
            'SPLINE3' : (SPLINE3, add_methods._add_spline_object),
            'SPLINE4' : (SPLINE4, add_methods._add_spline_object),
            'SPLINE5' : (SPLINE5, add_methods._add_spline_object),

            ## SOL 144
            'AEROS' : (AEROS, add_methods._add_aeros_object),
            'TRIM' : (TRIM, add_methods._add_trim_object),
            'TRIM2' : (TRIM2, add_methods._add_trim_object),
            'DIVERG' : (DIVERG, add_methods._add_diverg_object),

            ## SOL 145
            'AERO' : (AERO, add_methods._add_aero_object),
            'FLUTTER' : (FLUTTER, add_methods._add_flutter_object),
            'FLFACT' : (FLFACT, add_methods._add_flfact_object),
            'MKAERO1' : (MKAERO1, add_methods._add_mkaero_object),
            'MKAERO2' : (MKAERO2, add_methods._add_mkaero_object),

            'GUST' : (GUST, add_methods._add_gust_object),
            'CSSCHD' : (CSSCHD, add_methods._add_csschd_object),
            'MONPNT1' : (MONPNT1, add_methods._add_monpnt_object),
            'MONPNT2' : (MONPNT2, add_methods._add_monpnt_object),
            'MONPNT3' : (MONPNT3, add_methods._add_monpnt_object),
            'MONDSP1' : (MONDSP1, add_methods._add_monpnt_object),

            'NLPARM' : (NLPARM, add_methods._add_nlparm_object),
            'NLPCI' : (NLPCI, add_methods._add_nlpci_object),
            'TSTEP' : (TSTEP, add_methods._add_tstep_object),
            'TSTEP1' : (TSTEP1, add_methods._add_tstepnl_object),
            'TSTEPNL' : (TSTEPNL, add_methods._add_tstepnl_object),

            'TF' : (TF, add_methods._add_tf_object),
            'TIC' : (TIC, add_methods._add_tic_object),

            'DCONADD' : (DCONADD, add_methods._add_dconstr_object),
            'DCONSTR' : (DCONSTR, add_methods._add_dconstr_object),
            'DDVAL' : (DDVAL, add_methods._add_ddval_object),
            'DLINK' : (DLINK, add_methods._add_dlink_object),
            'DSCREEN' : (DSCREEN, add_methods._add_dscreen_object),

            'DTABLE' : (DTABLE, add_methods._add_dtable_object),
            'DRESP1' : (DRESP1, add_methods._add_dresp_object), # dresps
            'DRESP2' : (DRESP2, add_methods._add_dresp_object),
            'DRESP3' : (DRESP3, add_methods._add_dresp_object),
            'DVCREL1' : (DVCREL1, add_methods._add_dvcrel_object), # dvcrels
            'DVCREL2' : (DVCREL2, add_methods._add_dvcrel_object),
            'DVPREL1' : (DVPREL1, add_methods._add_dvprel_object), # dvprels
            'DVPREL2' : (DVPREL2, add_methods._add_dvprel_object),
            'DVMREL1' : (DVMREL1, add_methods._add_dvmrel_object), # ddvmrels
            'DVMREL2' : (DVMREL2, add_methods._add_dvmrel_object),
            'DVGRID' : (DVGRID, add_methods._add_dvgrid_object), # dvgrids

            ## tables
            'TABLES1' : (TABLES1, add_methods._add_table_object),
            'TABLEST' : (TABLEST, add_methods._add_table_object),
            'TABLEHT' : (TABLEHT, add_methods._add_table_object),
            'TABLEH1' : (TABLEH1, add_methods._add_table_object),

            ## dynamic tables
            'TABLED1' : (TABLED1, add_methods._add_tabled_object),
            'TABLED2' : (TABLED2, add_methods._add_tabled_object),
            'TABLED3' : (TABLED3, add_methods._add_tabled_object),
            'TABLED4' : (TABLED4, add_methods._add_tabled_object),

            ## material tables
            'TABLEM1' : (TABLEM1, add_methods._add_tablem_object),
            'TABLEM2' : (TABLEM2, add_methods._add_tablem_object),
            'TABLEM3' : (TABLEM3, add_methods._add_tablem_object),
            'TABLEM4' : (TABLEM4, add_methods._add_tablem_object),

            ## other tables
            'TABDMP1' : (TABDMP1, add_methods._add_table_sdamping_object),
            'TABRND1' : (TABRND1, add_methods._add_random_table_object),
            'TABRNDG' : (TABRNDG, add_methods._add_random_table_object),

            'EIGB' : (EIGB, add_methods._add_method_object),
            'EIGR' : (EIGR, add_methods._add_method_object),
            'EIGRL' : (EIGRL, add_methods._add_method_object),
            'EIGC' : (EIGC, add_methods._add_cmethod_object),
            'EIGP' : (EIGP, add_methods._add_cmethod_object),

            'BCRPARA' : (BCRPARA, add_methods._add_bcrpara_object),
            'BCTADD' : (BCTADD, add_methods._add_bctadd_object),
            'BCTPARA' : (BCTPARA, add_methods._add_bctpara_object),
            'BSURF' : (BSURF, add_methods._add_bsurf_object),
            'BSURFS' : (BSURFS, add_methods._add_bsurfs_object),
            ## 'BOUTPUT', 'BOLT', 'BOLTFOR', 'BOLTFRC',
            #'BOUTPUT': (Crash, None),
            #'BOLT': (Crash, None),
            #'BOLTFOR': (Crash, None),
            #'BOLTFRC': (Crash, None),

            #'RADCAV' : (RADCAV, self._add_radcav_object), #
            ##'RADLST' : (RADLST, add_methods._add_radcav_object), # TestOP2.test_bdf_op2_thermal_02
            ##'RADMTX' : (RADMTX, add_methods._add_radmtx_object), # TestOP2.test_bdf_op2_thermal_02
            ##'RADMT' : (Crash, None),

            'ASET' : (ASET, add_methods._add_aset_object),
            'ASET1' : (ASET1, add_methods._add_aset_object),

            'BSET' : (BSET, add_methods._add_bset_object),
            'BSET1' : (BSET1, add_methods._add_bset_object),

            'CSET' : (CSET, add_methods._add_cset_object),
            'CSET1' : (CSET1, add_methods._add_cset_object),

            'QSET' : (QSET, add_methods._add_qset_object),
            'QSET1' : (QSET1, add_methods._add_qset_object),

            'USET' : (USET, add_methods._add_uset_object),
            'USET1' : (USET1, add_methods._add_uset_object),

            'OMIT' : (OMIT, add_methods._add_omit_object),
            'OMIT1' : (OMIT1, add_methods._add_omit_object),

            'SET1' : (SET1, add_methods._add_set_object),
            'SET3' : (SET3, add_methods._add_set_object),

            ## radset
            #'RADSET' : (RADSET, self._add_radset_object),

            'SESET' : (SESET, add_methods._add_seset_object),

            'SEBSET' : (SEBSET, add_methods._add_sebset_object),
            'SEBSET1' : (SEBSET1, add_methods._add_sebset_object),

            'SECSET' : (SECSET, add_methods._add_secset_object),
            'SECSET1' : (SECSET1, add_methods._add_secset_object),

            'SEQSET' : (SEQSET, add_methods._add_seqset_object),
            'SEQSET1' : (SEQSET1, add_methods._add_seqset_object),

            ##'SESUP' : (SESUP, self._add_sesup_object),  # pseudo-constraint

            ##'SEUSET' : (SEUSET, add_methods._add_seuset_object),
            ##'SEUSET1' : (SEUSET1, add_methods._add_seuset_object),

            ## BCTSET
            'ROTORG' : (ROTORG, add_methods._add_rotor_object),
            'ROTORD' : (ROTORD, add_methods._add_rotor_object),

            #'DAREA' : (DAREA, self._add_darea_object),
            #'DPHASE' : (DPHASE, self._add_dphase_object),
            #'DELAY' : (DELAY, self._add_delay_object),

            #'CYJOIN' : (CYJOIN, self._add_cyjoin_object),
        }
        self._card_parser_prepare = {
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
            'TEMPBC' : self._prepare_tempbc,
            'CONVM' : self._prepare_convm,
            'CONV' : self._prepare_conv,
            'RADM' : self._prepare_radm,
            'RADBC' : self._prepare_radbc,
            # GRDSET-will be last card to update from _card_parser_prepare
            'GRDSET' : self._prepare_grdset,

            'BCTSET' : self._prepare_bctset,
            'ACMODL' : self._prepare_acmodl,

            # ------------------------------------------------------------------
            # vectorized
            # 'CONM1' : self._prepare_conm1
            'CONM2' : self._prepare_conm2,

            #'CMASS1' : self._prepare_cmass1,
            #'CMASS2' : self._prepare_cmass2,
            #'CMASS3' : self._prepare_cmass3,
            #'CMASS4' : self._prepare_cmass4,

            'CELAS1' : self._prepare_celas1,
            'CELAS2' : self._prepare_celas2,
            'CELAS3' : self._prepare_celas3,
            'CELAS4' : self._prepare_celas4,

            # ------------------------------------------------------------------
            'CDAMP1' : self._prepare_cdamp1,
            'CDAMP2' : self._prepare_cdamp2,
            'CDAMP3' : self._prepare_cdamp3,
            'CDAMP4' : self._prepare_cdamp4, # no (???)

            #'CBUSH1D' : self._prepare_cbush1d,
            #'CBUSH2D' : self._prepare_cbush2d,
            'CBUSH' : self._prepare_cbush,
            'CVISC':  self._prepare_cvisc,
            'PLOTEL':  self._prepare_plotel,
            # ------------------------------------------------------------------

            'CONROD' : self._prepare_conrod,
            'CROD':  self._prepare_crod,
            'CTUBE':  self._prepare_ctube,
            # ------------------------------------------------------------------

            'CBAR':  self._prepare_cbar,
            'CBEAM':  self._prepare_cbeam,
            # ------------------------------------------------------------------
            'GRID' : self._prepare_grid,
            # --------------------------------------
            'CTRIA3' : self._prepare_ctria3,
            'CTRIA6' : self._prepare_ctria6,
            'CTRIAR' : self._prepare_ctriar,
            'CQUAD4' : self._prepare_cquad4,
            'CQUAD8' : self._prepare_cquad8,

            'CQUAD' : self._prepare_cquad,
            'CQUADR' : self._prepare_cquadr,
            # shears
            'CSHEAR' : self._prepare_cshear,

            # -------------------------------------------
            # static loads
            'GRAV' : self._prepare_grav,
            'PLOAD' : self._prepare_pload,
            'PLOAD1' : self._prepare_pload1,
            'PLOAD2' : self._prepare_pload2,
            'PLOAD4' : self._prepare_pload4,
            'FORCE' : self._prepare_force,
            'FORCE1' : self._prepare_force1,
            'FORCE2' : self._prepare_force2,
            'MOMENT' : self._prepare_moment,
            'MOMENT1' : self._prepare_moment1,
            'MOMENT2' : self._prepare_moment2,
            # -------------------------------------------
            # dynamic loads
            'LOAD': self._prepare_load,
            'LSEQ': self._prepare_lseq,
            'SLOAD': self._prepare_sload,
            'SPCD': self._prepare_spcd,

            'TEMP': self._prepare_temp,
            'TEMPD': self._prepare_tempd,
            'PLOADX1': self._prepare_ploadx1,
            'ACCEL': self._prepare_accel,
            'ACCEL1': self._prepare_accel1,

            'QVOL': self._prepare_qvol,
            'QHBDY': self._prepare_qhbdy,
            'RFORCE': self._prepare_rforce,
            'RFORCE1': self._prepare_rforce1,
            'QBDY1': self._prepare_qbdy1,
            'QBDY2': self._prepare_qbdy2,
            'QBDY3': self._prepare_qbdy3,

            #'GMLOAD': self._prepare_gmload,
            'LOADCYN': self._prepare_loadcyn,
            'PRESAX': self._prepare_presax,
        }
        set_prep = set(self._card_parser_prepare)
        set_parse = set(self._card_parser)
        duplicate_names = list(set_prep.intersection(set_parse))
        duplicate_names.sort()
        assert len(duplicate_names) == 0, duplicate_names

        #new_reject_method = False
        #if new_reject_method:
            #self.cards_to_read.remove('CONM2')
            ## bwb_saero: 6.37-6.57 before using

            #all_cards = set(self._card_parser.keys() + self._card_parser_prepare.keys())
            #reject_cards = all_cards.difference(self.cards_to_read)
            ##print('reject_cards =', reject_cards)
            #for card_name in reject_cards:
                #if card_name in self._card_parser:
                    #del self._card_parser[card_name]
                #self._card_parser_prepare[card_name] = self._reject_card_obj2

    #def _reject_card_obj2(self, unused_card_name, card_obj):
        #"""rejects a card object"""
        #self.reject_cards.append(card_obj)

    def reject_card_lines(self, card_name: str, card_lines: list[str],
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

        #elif card_name in SOL_700 or card_name in MISSING_CARDS:
            #self.log.warning(f'    rejecting card_name = {card_name!r}')
            #return

        if card_name not in self.card_count:
            _check_for_spaces(card_name, card_lines, comment, self.log)
            #raise RuntimeError(card_name)
            if card_name == '\ufeff':
                self.log.warning(f'    rejecting card_name = {card_name!r}')
                self.log.warning('    comment:\n')
                print(comment)
                self.log.warning('    lines:\n')
                for line in card_lines:
                    print(line)
            elif show_log:
                self.log.info(f'    rejecting card_name = {card_name!r}')
            assert isinstance(show_log, bool), show_log
        self.increase_card_count(card_name)
        self.reject_lines.append([_format_comment(comment)] + card_lines)

    def _write_reject_message(self, card_name, unused_card_obj, comment=''):
        """common method to not write duplicate reject card names"""
        if card_name not in self.card_count:
            #if ' ' in card_name:
                #_check_for_spaces(card_name, card_lines, comment, self.log)
            self.log.info('    rejecting card_name = %s' % card_name)

    def _prepare_plotel(self, unused_card: list[str], card_obj: BDFCard, comment='') -> None:
        """adds a PLOTEL"""
        #['PLOTEL', '3101', '3101', '3102', None, '3102', '3102', '3103']
        plotels = [PLOTEL.add_card(card_obj, 0, comment=comment)]
        if card_obj.field(5):  # eid
            plotels.append(PLOTEL.add_card(card_obj, 1, comment=''))
        for plotel in plotels:
            self._add_methods._add_plotel_object(plotel)
        return plotels

    #def _prepare_cbar(self, unused_card: list[str], card_obj: BDFCard, comment='') -> None:
        #"""adds a CBAR"""
        #elem = CBAR.add_card(card_obj, baror=self.baror, comment=comment)
        #self._add_element_object(elem)
        #return elem

    def _prepare_cbeam(self, unused_card: list[str], card_obj: BDFCard, comment='') -> None:
        """adds a CBEAM"""
        elem = CBEAM.add_card(card_obj, beamor=self.beamor, comment=comment)
        self._add_methods._add_element_object(elem)
        return elem

    def _prepare_ctetra(self, unused_card: list[str], card_obj: BDFCard, comment='') -> None:
        """adds a CTETRA4/CTETRA10"""
        if len(card_obj) == 7:
            elem = CTETRA4.add_card(card_obj, comment=comment)
        else:
            elem = CTETRA10.add_card(card_obj, comment=comment)
        self._add_methods._add_element_object(elem)
        return elem

    def _prepare_cpyram(self, unused_card: list[str], card_obj: BDFCard, comment='') -> None:
        """adds a CPYRAM5/CPYRAM13"""
        if len(card_obj) == 8:
            elem = CPYRAM5.add_card(card_obj, comment=comment)
        else:
            elem = CPYRAM13.add_card(card_obj, comment=comment)
        self._add_methods._add_element_object(elem)
        return elem

    def _prepare_cpenta(self, unused_card: list[str], card_obj: BDFCard, comment='') -> None:
        """adds a CPENTA6/CPENTA15"""
        if len(card_obj) == 9:
            elem = CPENTA6.add_card(card_obj, comment=comment)
        else:
            elem = CPENTA15.add_card(card_obj, comment=comment)
        self._add_methods._add_element_object(elem)
        return elem

    def _prepare_chexa(self, unused_card: list[str], card_obj: BDFCard, comment='') -> None:
        """adds a CHEXA8/CHEXA20"""
        if len(card_obj) == 11:
            elem = CHEXA8.add_card(card_obj, comment=comment)
        else:
            elem = CHEXA20.add_card(card_obj, comment=comment)
        self._add_methods._add_element_object(elem)
        return elem

    def _prepare_bctset(self, unused_card: list[str], card_obj: BDFCard, comment='') -> None:
        """adds a BCTSET"""
        bctset = BCTSET.add_card(card_obj, comment=comment, sol=self.sol)
        self._add_methods._add_bctset_object(bctset)
        return bctset

    def _prepare_grdset(self, unused_card: list[str], card_obj: BDFCard, comment='') -> None:
        """adds a GRDSET"""
        self.grdset = GRDSET.add_card(card_obj, comment=comment)
        return self.grdset

    def _prepare_cdamp4(self, unused_card: list[str], card_obj: BDFCard, comment='') -> None:
        """adds a CDAMP4"""
        dampers = [CDAMP4.add_card(card_obj, comment=comment)]
        if card_obj.field(5):
            dampers.append(CDAMP4.add_card(card_obj, 1, comment=''))
        for damper in dampers:
            self._add_methods._add_damper_object(damper)
        return dampers

    def _prepare_deform(self, unused_card: list[str], card_obj: BDFCard, comment='') -> None:
        """adds a DEFORM"""
        loads = [DEFORM.add_card(card_obj, comment=comment)]
        if card_obj.field(4):
            loads.append(DEFORM.add_card(card_obj, 1, comment=comment))
        if card_obj.field(6):
            loads.append(DEFORM.add_card(card_obj, 2, comment=comment))
        for loadi in loads:
            self._add_methods._add_load_object(loadi)
        return loads

    def _prepare_tempbc(self, unused_card: list[str], card_obj: BDFCard, comment='') -> None:
        """adds a TEMPBC"""
        boundary_condition = TEMPBC.add_card(card_obj, comment=comment)
        self._add_methods._add_thermal_bc_object(boundary_condition, boundary_condition.eid)
        return boundary_condition

    def _prepare_convm(self, unused_card: list[str], card_obj: BDFCard, comment='') -> None:
        """adds a CONVM"""
        boundary_condition = CONVM.add_card(card_obj, comment=comment)
        self._add_methods._add_thermal_bc_object(boundary_condition, boundary_condition.eid)
        return boundary_condition

    def _prepare_conv(self, unused_card: list[str], card_obj: BDFCard, comment='') -> None:
        """adds a CONV"""
        boundary_condition = CONV.add_card(card_obj, comment=comment)
        self._add_methods._add_thermal_bc_object(boundary_condition, boundary_condition.eid)
        return boundary_condition

    def _prepare_radm(self, unused_card: list[str], card_obj: BDFCard, comment='') -> None:
        """adds a RADM"""
        boundary_condition = RADM.add_card(card_obj, comment=comment)
        self._add_methods._add_thermal_bc_object(boundary_condition, boundary_condition.radmid)
        return boundary_condition

    def _prepare_radbc(self, unused_card: list[str], card_obj: BDFCard, comment='') -> None:
        """adds a RADBC"""
        boundary_condition = RADBC.add_card(card_obj, comment=comment)
        self._add_methods._add_thermal_bc_object(boundary_condition, boundary_condition.nodamb)
        return boundary_condition

    def _prepare_tempd(self, unused_card: list[str], card_obj: BDFCard, comment='') -> None:
        """adds a TEMPD"""
        tempds = [TEMPD.add_card(card_obj, 0, comment=comment)]
        if card_obj.field(3):
            tempds.append(TEMPD.add_card(card_obj, 1, comment=''))
            if card_obj.field(5):
                tempds.append(TEMPD.add_card(card_obj, 2, comment=''))
                if card_obj.field(7):
                    tempds.append(TEMPD.add_card(card_obj, 3, comment=''))
        for tempd in tempds:
            self._add_methods._add_tempd_object(tempd)
        return tempds

    def _prepare_tempax(self, unused_card: list[str], card_obj: BDFCard, comment='') -> None:
        """adds a TEMPAX"""
        tempaxs = [TEMPAX.add_card(card_obj, 0, comment=comment)]
        if card_obj.field(5):
            tempaxs.append(TEMPAX.add_card(card_obj, 1, comment=''))
        for tempax in tempaxs:
            self._add_methods._add_load_object(tempax)
        return tempaxs

    def _prepare_dequatn(self, unused_card: list[str], card_obj: BDFCard, comment='') -> None:
        """adds a DEQATN"""
        deqatn = DEQATN.add_card(card_obj, comment=comment)
        self._add_methods._add_deqatn_object(deqatn)
        return deqatn

    def _prepare_dti(self, unused_card_name, card_obj, comment=''):
        """adds a DTI"""
        #name = string(card_obj, 1, 'name')
        dti = DTI.add_card(card_obj, comment=comment)
        self._add_methods._add_dti_object(dti)
        return dti

    def _prepare_dmig(self, unused_card: list[str], card_obj: BDFCard, comment='') -> None:
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

    def _prepare_dmiax(self, unused_card: list[str], card_obj: BDFCard, comment='') -> None:
        """adds a DMIAX"""
        return self._prepare_dmix(DMIAX, self._add_dmiax_object, card_obj, comment=comment)

    def _prepare_dmi(self, unused_card: list[str], card_obj: BDFCard, comment='') -> None:
        """adds a DMI"""
        return self._prepare_dmix(DMI, self._add_dmi_object, card_obj, comment=comment)

    def _prepare_dmij(self, unused_card: list[str], card_obj: BDFCard, comment='') -> None:
        """adds a DMIJ"""
        return self._prepare_dmix(DMIJ, self._add_dmij_object, card_obj, comment=comment)

    def _prepare_dmik(self, unused_card: list[str], card_obj: BDFCard, comment='') -> None:
        """adds a DMIK"""
        return self._prepare_dmix(DMIK, self._add_dmik_object, card_obj, comment=comment)

    def _prepare_dmiji(self, unused_card: list[str], card_obj: BDFCard, comment='') -> None:
        """adds a DMIJI"""
        return self._prepare_dmix(DMIJI, self._add_dmiji_object, card_obj, comment=comment)

    def _prepare_cmass4(self, unused_card: list[str], card_obj: BDFCard, comment='') -> None:
        """adds a CMASS4"""
        elements = [CMASS4.add_card(card_obj, icard=0, comment=comment)]
        if card_obj.field(5):
            elements.append(CMASS4.add_card(card_obj, icard=1, comment=comment))
        for elem in elements:
            self._add_methods._add_mass_object(elem)
        return elements

    def _prepare_pelas(self, unused_card: list[str], card_obj: BDFCard, comment='') -> None:
        """adds a PELAS"""
        properties = [PELAS.add_card(card_obj, icard=0, comment=comment)]
        if card_obj.field(5):
            properties.append(PELAS.add_card(card_obj, icard=1, comment=comment))
        for prop in properties:
            self._add_methods._add_property_object(prop)
        return properties

    def _prepare_nsm(self, unused_card: list[str], card_obj: BDFCard, comment='') -> None:
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
            self._add_methods._add_nsm_object(nsm)
        return nsms

    def _prepare_nsml(self, unused_card: list[str], card_obj: BDFCard, comment='') -> None:
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
            self._add_methods._add_nsm_object(nsm)
        return nsms

    def _prepare_pvisc(self, unused_card: list[str], card_obj: BDFCard, comment='') -> None:
        """adds a PVISC"""
        properties = [PVISC.add_card(card_obj, icard=0, comment=comment)]
        if card_obj.field(5):
            properties.append(PVISC.add_card(card_obj, icard=1, comment=comment))
        for prop in properties:
            self._add_methods._add_property_object(prop)
        return properties

    def _prepare_ringfl(self, unused_card: list[str], card_obj: BDFCard, comment='') -> None:
        """adds a RINGFL"""
        rings = [RINGFL.add_card(card_obj, icard=0, comment=comment)]
        if card_obj.field(3):
            rings.append(RINGFL.add_card(card_obj, icard=1, comment=comment))
        if card_obj.field(5):
            rings.append(RINGFL.add_card(card_obj, icard=2, comment=comment))
        for ring in rings:
            self._add_methods._add_ringfl_object(ring)
        return rings

    def _prepare_pdamp(self, unused_card: list[str], card_obj: BDFCard, comment='') -> None:
        """adds a PDAMP"""
        properties = [PDAMP.add_card(card_obj, icard=0, comment=comment)]
        if card_obj.field(3):
            properties.append(PDAMP.add_card(card_obj, icard=1, comment=comment))
        if card_obj.field(5):
            properties.append(PDAMP.add_card(card_obj, icard=2, comment=comment))
        if card_obj.field(7):
            properties.append(PDAMP.add_card(card_obj, icard=3, comment=comment))
        for prop in properties:
            self._add_methods._add_property_object(prop)
        return properties

    def _prepare_pmass(self, unused_card: list[str], card_obj: BDFCard, comment='') -> None:
        """adds a PMASS"""
        properties = [PMASS.add_card(card_obj, icard=0, comment=comment)]
        for (i, j) in enumerate([3, 5, 7]):
            if card_obj.field(j):
                properties.append(PMASS.add_card(card_obj, icard=i+1, comment=comment))
        for prop in properties:
            self._add_methods._add_property_mass_object(prop)
        return properties

    def _prepare_cord1r(self, unused_card: list[str], card_obj: BDFCard, comment='') -> None:
        """adds a CORD1R"""
        coords = [CORD1R.add_card(card_obj, comment=comment)]
        if card_obj.field(5):
            coords.append(CORD1R.add_card(card_obj, icard=1, comment=comment))
        for coord in coords:
            self._add_methods._add_coord_object(coord)
        return coords

    def _prepare_cord1c(self, unused_card: list[str], card_obj: BDFCard, comment='') -> None:
        """adds a CORD1C"""
        coords = [CORD1C.add_card(card_obj, comment=comment)]
        if card_obj.field(5):
            coords.append(CORD1C.add_card(card_obj, icard=1, comment=comment))
        for coord in coords:
            self._add_methods._add_coord_object(coord)
        return coords

    def _prepare_cord1s(self, unused_card: list[str], card_obj: BDFCard, comment='') -> None:
        """adds a CORD1S"""
        coords = [CORD1S.add_card(card_obj, comment=comment)]
        if card_obj.field(5):
            coords.append(CORD1S.add_card(card_obj, icard=1, comment=comment))
        for coord in coords:
            self._add_methods._add_coord_object(coord)
        return coords

    def _prepare_acmodl(self, unused_card: list[str], card_obj: BDFCard, comment='') -> None:
        acmodl = ACMODL.add_card(card_obj, self._nastran_format, comment=comment)
        self._add_methods._add_acmodl_object(acmodl)

    def add_card_ifile(self, ifile: int, card_lines: list[str], card_name: str,
                       comment: str='', is_list: bool=True, has_none: bool=True) -> Any:
        """Same as ``add_card`` except it has an ifile parameter"""
        assert isinstance(ifile, (int, np.int32)), 'ifile=%s type=%s' % (ifile, type(ifile))
        card_name = card_name.upper()
        card_obj, unused_card = self.create_card_object(
            card_lines, card_name,
            is_list=is_list, has_none=has_none)
        self._add_card_helper_ifile(ifile, card_obj, card_name, card_name, comment)
        return card_obj

    def add_card(self, card_lines: list[str], card_name: str,
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

    def _get_npoints_nids_allnids(self):
        """helper method for get_xyz_in_coord"""
        return self.nodes.nids

        #nnodes = len(self.nodes)
        #nspoints = 0
        #nepoints = 0
        #nids = list(self.node_ids)
        #all_nodes = list(self.node_ids)
        #if self.spoints:
            #spoints = list(self.spoints)
            #nspoints = len(spoints)
            #all_nodes += spoints
        #if self.epoints:
            #epoints = list(self.epoints)
            #nepoints = len(epoints)
            #all_nodes += epoints
        ##self.log.debug('all_nodes = %s' % all_nodes)

        #npoints = nnodes + nspoints + nepoints
        #if len(all_nodes) != npoints:
            #msg = 'len(unique(all_nodes))=%s npoints=%s\n' % (len(np.unique(all_nodes)), npoints)
            #msg += 'npoints = nnodes+nspoints+nepoints = %s + %s + %s\n' % (
                #nnodes, nspoints, nepoints)
            #msg += 'all_nodes=%s' % (all_nodes)
            #raise RuntimeError(msg)
        #if npoints == 0:
            #msg = 'nnodes=%s nspoints=%s nepoints=%s' % (nnodes, nspoints, nepoints)
            #raise ValueError(msg)
        #return npoints, nids, all_nodes

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
                msg = f'problem adding {card_obj}'
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
                raise RuntimeError(f'_prepare_{card_name.lower()} needs to implement obj')
            elif isinstance(obj, list):
                for obji in obj:
                    obji.ifile = ifile
            elif obj == -1:
                pass
            else:
                obj.ifile = ifile

        else:
            self.reject_cards.append(card_obj)

    def _add_card_helper(self, card_obj: BDFCard, card: list[str],
                         card_name: str, ifile: int,
                         comment: str='') -> None:
        """
        Adds a card object to the BDF object.

        Parameters
        ----------
        card_object : BDFCard()
            the card object representation of card
        card : list[str]
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
                print('problem adding %s' % card_obj)
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

    def get_bdf_stats(self, return_type: str='string') -> str | list[str]:
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
        if self.spoints:
            spoints = list(self.spoints)
            nspoints = len(spoints)
        if self.epoints:
            epoints = list(self.epoints)
            nepoints = len(epoints)

        if nnodes + nspoints + nepoints == 0:
            msg = 'nnodes=%s nspoints=%s nepoints=%s' % (
                nnodes, nspoints, nepoints)
            raise ValueError(msg)

        if idtype == 'int32':
            try:
                out = _set_nodes(self, spoints, epoints,
                                 nnodes, nspoints, nepoints,
                                 idtype, fdtype)
            except OverflowError:
                out = _set_nodes(self, spoints, epoints,
                                 nnodes, nspoints, nepoints,
                                 'int64', fdtype)
        else:
            out = _set_nodes(self, spoints, epoints,
                             nnodes, nspoints, nepoints,
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
            icd_transform[cd] = np.where(np.isin(nids_all, nids))[0]

        for cp, nids in sorted(nids_cp_transform.items()):
            if cp in [-1]:
                continue
            nids = np.array(nids)
            icp_transform[cp] = np.where(np.isin(nids_all, nids))[0]
        return icd_transform, icp_transform, xyz_cp, nid_cp_cd

    def get_xyz_in_coord_array(self, cid: int=0,
                               fdtype: str='float64', idtype: str='int32') -> tuple[Any, Any, Any,
                                                                                    dict[int, Any], dict[int, Any]]:
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
        # this is used when xref=False (only for vectorized=True)
        # we now require nids, where the other approach
        # (the one with xref=True) does not
        in_place = False
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
                out = _get_coords_to_update(
                    self.coords, cps_to_check, cps_checked, nids_checked)
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
                    i1, i2, i3 = np.searchsorted(nids, coord.node_ids)
                    assert nids[i1] == nid1
                    assert nids[i2] == nid2
                    assert nids[i3] == nid3
                    coord.e1 = xyz_cid0[i1, :] #: the origin in the local frame
                    coord.e2 = xyz_cid0[i2, :] #: a point on the z-axis
                    coord.e3 = xyz_cid0[i3, :] #: a point on the xz-plane
                    #else:
                        #g1_ref = nodes[nid1]
                        #g2_ref = nodes[nid2]
                        #g3_ref = nodes[nid3]
                        #coord.e1 = g1_ref.get_position() #: the origin in the local frame
                        #coord.e2 = g2_ref.get_position() #: a point on the z-axis
                        #coord.e3 = g3_ref.get_position() #: a point on the xz-plane
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
                    msg += f'  rid={coord.rid!r} origin={rid_ref.origin}\n\n'
                else:
                    nid1, nid2, nid3 = coord.node_ids
                    #coord.e1 = xyz_cid0[i1, :] #: the origin in the local frame
                    #coord.e2 = xyz_cid0[i2, :] #: a point on the z-axis
                    #coord.e3 = xyz_cid0[i3, :] #: a point on the xz-plane
                    #if self.is_bdf_vectorized:
                    i1, i2, i3 = np.searchsorted(nids, coord.node_ids)
                    cp1 = nodes.cp[i1]
                    cp2 = nodes.cp[i2]
                    cp3 = nodes.cp[i3]
                    #else:
                        #cp1 = nodes[nid1].cp
                        #cp2 = nodes[nid2].cp
                        #cp3 = nodes[nid3].cp
                    msg += f'  g1={nid1} xyz={coord.e1} cp={cp1}\n'
                    msg += f'  g2={nid2} xyz={coord.e2} cp={cp2}\n'
                    msg += f'  g3={nid3} xyz={coord.e3} cp={cp3}\n'
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
        cps_to_check0 : list[int]
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
        cps_checked : list[int]
            the Cps that were checked
        cps_to_check : list[int]
            the Cps that are unreferenceable given the current information

        """
        nids_checked, cps_checked, cps_to_check = transform_coords_vectorized(
            cps_to_check0, icp_transform,
            nids, xyz_cp, xyz_cid0, xyz_cid0_correct,
            self.coords, do_checks)
        return nids_checked, cps_checked, cps_to_check

    #@property
    #def is_bdf_vectorized(self):
        #"""Returns False for the ``BDF`` class"""
        #return hasattr(self, 'grid')

    def get_displacement_index(self) -> tuple[Any, Any, dict[int, Any]]:
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
            icd_transform[cid] = np.where(np.isin(nids_all, nids))[0]
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
        if card_name in self.card_count:
            self.card_count[card_name] += count_num
        else:
            assert '=' not in card_name, card_name
            self.card_count[card_name] = count_num

    def _old_card_fields(self, card_lines: list[str], card_name: str,
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
            assert len(cards2) == 1, f'cards2={cards2}; ncards={len(cards2)}'
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
                assert ifield == 0, f'ifield={ifield} field={field!r} new_card={new_card}'
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
                except Exception:
                    self.log.error(f'old_card:{old_card}\nnew_card:\n{new_card}')
                    raise
            else:
                assert '(' not in field, f'field={field!r}'
                assert '*' not in field, f'field={field!r}'
                assert '=' not in field, f'field={field!r}'
                field2 = field
            #print(dig_str, ' %i: %r -> %r' % (ifield, field, field2))
            card.append(field2)
        if card:
            cards.append(card)
            #print(dig_str, 'card_expanded = %s' % card)
        else:  # pragma: no cover
            raise RuntimeError(card)
        return cards

    def _parse_cards(self, cards_list: list[list[str]],
                     cards_dict: dict[str, list[str]],
                     card_count: dict[str, int]) -> None:
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

    def _parse_cards_dict(self, cards_dict: dict[str, list[str]]) -> None:
        """parses the cards that are in dictionary format"""
        if self.save_file_structure:
            raise NotImplementedError('save_file_structure=True is not supported\n%s' % (
                list(cards_dict.keys())))

        for card_name, cards in sorted(cards_dict.items()):
            if self.is_reject(card_name):
                self.log.info(f'    rejecting card_name = {card_name}')
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
                    msg = f'card_name = {card_name!r}\n'
                    msg += f'card_lines = {card_lines}'
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
                    msg = f"save_file_structure=True doesn't support {card_name}"
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
                    msg = f'card_name = {card_name!r}\n'
                    msg += f'card_lines = {card_lines}'
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

    def _parse_primary_file_header(self, bdf_filename: str | StringIO) -> None:
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
        if isinstance(bdf_filename, (str, PurePath)):
            try:
                with open(bdf_filename, 'r') as bdf_file:
                    lines = bdf_file.readlines()
            except UnicodeDecodeError:
                with open(bdf_filename, 'r', errors='replace') as bdf_file:
                    line = 'temp'
                    lines = []
                    while line:
                        line = bdf_file.readline()
                        lines.append(line)
                        try:
                            # try to force a crash
                            unused_bytes_line = line.encode('ascii')  # TODO: use the encoding
                        except UnicodeEncodeError:
                            break
                    n = 20
                    i0 = len(lines) - n
                    for i, line in enumerate(lines[-n:-1]):
                        self.log.debug(f'Line {i0+i}: {line.strip()!r}')
                    self.log.error(f'Line {i0+i+1}: {lines[-1].strip()!r}')
                    raise
        else:
            # StringIO
            if hasattr(bdf_filename, 'read') and hasattr(bdf_filename, 'write'):
                lines = bdf_filename.readlines()
                bdf_filename.seek(0)  # need to rewind the buffer!

        self._check_pynastran_header(lines, check_header=True)
        map_update(self, self.nastran_format)

    def _update_for_nastran(self):
        """updates for msc/nx/optistruct"""
        # TODO: undo the changes for zona
        card_parser = self._card_parser
        CARD_MAP['PARAM'] = PARAM
        card_parser['PARAM'] = (PARAM, self._add_methods._add_param_object)
        self.add_param = self._add_param_nastran

    def _update_for_optistruct(self):
        """updates for mystran"""
        self._update_for_nastran() # copies this...
        card_parser = self._card_parser
        CARD_MAP['PBUSH_OPTISTRUCT'] = PBUSH_OPTISTRUCT
        card_parser['PBUSH'] = (PBUSH_OPTISTRUCT, self._add_methods._add_property_object)

    def _update_for_mystran(self):
        """updates for mystran"""
        card_parser = self._card_parser
        CARD_MAP['PARAM'] = PARAM_MYSTRAN
        card_parser['PARAM'] = (PARAM_MYSTRAN, self._add_methods._add_param_object)
        self.add_param = self._add_param_mystran

    def _update_for_nasa95(self):
        """updates for nasa95"""
        CARD_MAP['PARAM'] = PARAM_NASA95
        card_parser = self._card_parser
        card_parser['PARAM'] = (PARAM_NASA95, self._add_methods._add_param_object)
        self.add_param = self._add_param_nasa95

    def _check_pynastran_header(self, lines: list[str],
                                check_header: bool=True) -> None:
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
                assert value.lower() in ['msc', 'nx', 'optistruct', 'zona', 'nasa95', 'mystran'], f'version={value!r} is not supported'
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
                    raise RuntimeError(f'{type_to_skip!r} is an invalid key')
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
            exec(self.code_block, sys._getframe().f_locals)

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
        self.log.debug(f'---starting BDF.read_bdf of {self.bdf_filename}---')
        self._parse_primary_file_header(bdf_filename)

        obj = BDFInputPy(self.read_includes, self.dumplines, self._encoding,
                         nastran_format=self.nastran_format,
                         log=self.log, debug=self.debug)
        obj.use_new_parser = self.use_new_deck_parser

        out = obj.get_lines(bdf_filename, punch=self.punch, make_ilines=True)
        (system_lines, executive_control_lines, case_control_lines,
         bulk_data_lines, bulk_data_ilines,
         superelement_lines, superelement_ilines) = out
        self._set_pybdf_attributes(obj, save_file_structure=False)

        self.system_command_lines = system_lines
        self.executive_control_lines = executive_control_lines
        self.case_control_lines = case_control_lines

        sol, method, sol_iline, app = parse_executive_control_deck(executive_control_lines)
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

    def create_subcases(self, subcase_ids: int | list[int] | None=None) -> dict[int, Subcase]:
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

    #def _parse_cards_hdf5(self, cards: dict[str, Any], unused_card_count):
        #"""creates card objects and adds the parsed cards to the deck"""
        #self.echo = False
        #cards_out = {}
        #for card_name, card in sorted(cards.items()):
            #cards_list = []
            #cards_out[card_name] = cards_list
            #if self.is_reject(card_name):
                #self.log.info(f'    rejecting card_name = {card_name}')
                #for comment, card_lines, unused_ifile_iline in card:
                    #self.increase_card_count(card_name)
                    #self.reject_lines.append([_format_comment(comment)] + card_lines)
            #else:
                #for comment, card_lines, unused_ifile_iline in card:
                    #class_instance = self._add_card_hdf5(card_lines, card_name, comment=comment,
                                                         #is_list=False, has_none=False)
                    #cards_list.append(class_instance)
        #return cards_out


    #def _add_card_hdf5(self, card_lines: list[str], card_name: str,
                       #comment='', is_list: bool=True, has_none: bool=True) -> Any:
        #"""
        #Creates a BaseCard object that will be used to simplify HDF5 adding.

        #Parameters
        #----------
        #card_lines: list[str]
            #the list of the card fields
        #card_name : str
            #the card_name -> 'GRID'
        #comment : str
            #an optional the comment for the card
        #is_list : bool, optional
            #False : input is a list of card fields -> ['GRID', 1, None, 3.0, 4.0, 5.0]
            #True :  input is list of card_lines -> ['GRID, 1,, 3.0, 4.0, 5.0']
        #has_none : bool; default=True
            #can there be trailing Nones in the card data (e.g. ['GRID', 1, 2, 3.0, 4.0, 5.0, None])
            #can there be trailing Nones in the card data (e.g. ['GRID, 1, 2, 3.0, 4.0, 5.0, '])

        #Returns
        #-------
        #class_instance : BaseCard()
            #the card object representation of card

        #.. code-block:: python

           #>>> model = BDF()

           ## is_list is a somewhat misleading name; is it a list of card_lines
           ## where a card_line is an unparsed string
           #>>> card_lines = ['GRID,1,2']
           #>>> comment = 'this is a comment'
           #>>> model.add_card(card_lines, 'GRID', comment, is_list=True)

           ## here is_list=False because it's been parsed
           #>>> card = ['GRID', 1, 2,]
           #>>> model.add_card(card_lines, 'GRID', comment, is_list=False)

           ## here is_list=False because it's been parsed
           ## Note the None at the end of the 1st line, which is there
           ##      because the CONM2 card has a blank field.
           ##      It must be there.
           ## We also set i32 on the 2nd line, so it will default to 0.0
           #>>> card = [
                   #'CONM2', eid, nid, cid, mass, x1, x2, x3, None,
                            #i11, i21, i22, i31, None, i33,
               #]
           #>>> model.add_card(card_lines, 'CONM2', comment, is_list=False)

           ## here's an alternate approach for the CONM2
           ## we use Nastran's CSV format
           ## There are many blank fields, but it's parsed exactly like a
           ## standard CONM2.
           #>>> card = [
                   #'CONM2,1,2,3,10.0',
                   #',1.0,,5.0'
               #]
           #>>> model.add_card(card_lines, 'CONM2', comment, is_list=True)

        #Notes
        #-----
        #This is a very useful method for interfacing with the code.

        #The card_object is not a card-type object...so not a GRID
        #card or CQUAD4 object.  It's a BDFCard Object.  However,
        #you know the type (assuming a GRID), so just call the
        #*mesh.Node(nid)* to get the Node object that was just
        #created.

        #"""
        #card_name = card_name.upper()
        #card_obj, unused_card = self.create_card_object(
            #card_lines, card_name,
            #is_list=is_list, has_none=has_none)
        #class_instance = self._add_card_helper_hdf5(card_obj, card_name, card_name, comment)
        #return class_instance

    #def _add_card_helper_hdf5(self, card_obj: BDFCard,
                              #card: list[str],
                              #card_name: str, comment: str='') -> None:
        #"""
        #Adds a card object to the BDF object.

        #Parameters
        #----------
        #card_object : BDFCard()
            #the card object representation of card
        #card : list[str]
            #the fields of the card object; used for rejection and special cards
        #card_name : str
            #the card_name -> 'GRID'
        #comment : str
            #an optional the comment for the card

        #Returns
        #-------
        #class_instance : obj or None
            #obj : GRID, CQUAD4, ...
            #None : ECHOON, ECHOOFF, reject

        #"""
        #if card_name == 'ECHOON':
            #self.echo = True
            #return None
        #elif card_name == 'ECHOOFF':
            #self.echo = False
            #return None

        #if self.echo and not self.force_echo_off:
            #try:
                #print(print_card_8(card_obj).rstrip())
            #except Exception:
                #if card in ['DEQATN']:
                    #print(str(card_obj).rstrip())
                #else:
                    #print(print_card_16(card_obj).rstrip())

        #if card_name in self._card_parser:
            #card_class, add_card_function = self._card_parser[card_name]

            ## simplified, so no error catching
            #class_instance = card_class.add_card(card_obj, comment=comment)
            ##add_card_function(class_instance)

        #elif card_name in self._card_parser_prepare:
            #add_card_function = self._card_parser_prepare[card_name]
            ## simplified, so no error catching
            #class_instance = add_card_function(card, card_obj, comment=comment)

        #else:
            #self.reject_cards.append(card_obj)
            #class_instance = None
        #return class_instance


class BDF(BDF_):
    """
    Uses the BDF class, but overwrites a few classes/methods

    Vectorized:
     - GRID
    """
    is_bdf_vectorized = True
    def __init__(self, debug: Optional[bool]=True,
                 log: Optional[SimpleLogger]=None,
                 mode: str='msc'):
        BDF_.__init__(self, debug=debug, log=log, mode=mode)
        self.use_new_deck_parser = True
        #super(BDF, self).__init__(debug=debug, log=log, mode=mode)
        #self._grids_temp = []

        model = self
        self.grid = GRIDv(model)
        self.nodes = Nodes(model)

        self.celas1 = CELAS1(model)
        self.celas2 = CELAS2(model)
        self.celas3 = CELAS3(model)
        self.celas4 = CELAS4(model)
        self.springs = Springs(model)

        self.cdamp1 = CDAMP1(model)
        self.cdamp2 = CDAMP2(model)
        self.cdamp3 = CDAMP3(model)
        self.cdamp4 = CDAMP4(model)
        #self.cdamp5 = CDAMP5(model)    # TODO: temp
        self.cvisc = CVISCv(model)    # this is in dampers right now
        self.plotel = PLOTELv(model)  # this is in dampers right now
        self.dampers = Dampers(model)

        self.cbush = CBUSHv(model)
        self.bushes = Bushes(model)

        #self.conm1 = CONM1v(model)
        self.conm2 = CONM2v(model)
        #self.cmass1 = CMASS1v(model)
        #self.cmass2 = CMASS2v(model)
        #self.cmass3 = CMASS3v(model)
        #self.cmass4 = CMASS4v(model)
        self.masses2 = Masses(model)

        self.crod = CRODv(model)
        self.conrod = CONRODv(model)
        self.ctube = CTUBEv(model)
        self.rods = Rods(model)

        self.cbar = CBARv(model)
        self.bars = Bars(model)

        self.cbeam = CBEAMv(model)
        self.beams = Beams(model)

        self.ctria3 = CTRIA3v(model)
        self.cquad4 = CQUAD4v(model)
        self.ctria6 = CTRIA6v(model)
        self.cquad8 = CQUAD8v(model)
        self.cquad = CQUADv(model)
        self.cquadr = CQUADRv(model)
        self.ctriar = CTRIARv(model)

        self.shells = Shells(model)
        #self.pshell = PSHELLv(model)  # TODO: temp

        self.cshear = CSHEARv(model)
        self.shears = Shears(model)

        #self.ctriax = CTRIA3v(model)   # TODO: temp
        #self.cquadx = CTRIA3v(model)   # TODO: temp
        #self.ctriax6 = CTRIA3v(model)  # TODO: temp
        #self.cquadx8 = CTRIA3v(model)  # TODO: temp

        self.ctetra4 = CTETRA4v(model)
        self.ctetra10 = CTETRA10v(model)
        self.chexa8 = CHEXA8v(model)
        self.chexa20 = CHEXA20v(model)
        self.cpenta6 = CPENTA6v(model)
        self.cpenta15 = CPENTA15v(model)
        self.cpyram5 = CPYRAM5v(model)
        self.cpyram13 = CPYRAM13v(model)
        self.solids = Solids(model)

        self.elements2 = Elements(model)  # TODO: change this name


        self.sload = SLOADv(model)
        self.grav = GRAVv(model)
        self.force = FORCEv(model)
        self.force1 = FORCE1v(model)
        self.force2 = FORCE2v(model)
        self.pload = PLOADv(model)
        self.pload1 = PLOAD1v(model)
        self.pload2 = PLOAD2v(model)
        self.pload4 = PLOAD4v(model)

        self.moment = MOMENTv(model)
        self.moment1 = MOMENT1v(model)
        self.moment2 = MOMENT2v(model)

        self.spcd = SPCDv(model)
        self.temp = TEMPv(model)
        self.tempd = TEMPDv(model)

        self.load_combinations = {}
        #def lseqi():
            #return LSEQv(model)
        #self.lseqs = defaultdict(lseqi)
        self.lseq = LSEQv(model)
        self.loads = Loads(model)

        self._update_card_parser()

    def clear_attributes(self) -> None:
        """removes the attributes from the model"""
        self.log.info('clearing vectorized BDF model')
        #self.__init_attributes()
        model = self
        #self.superelement_models = {}
        self.grid = GRIDv(model)
        self.nodes = Nodes(model)

        self.celas1 = CELAS1(model)
        self.celas2 = CELAS2(model)
        self.celas3 = CELAS3(model)
        self.celas4 = CELAS4(model)
        self.springs = Springs(model)

        self.cdamp1 = CDAMP1(model)
        self.cdamp2 = CDAMP2(model)
        self.cdamp3 = CDAMP3(model)
        self.cdamp4 = CDAMP4(model)
        #self.cdamp5 = CDAMP5(model)    # TODO: temp
        self.cvisc = CVISCv(model)    # this is in dampers right now
        self.plotel = PLOTELv(model)  # this is in dampers right now
        self.dampers = Dampers(model)

        self.cbush = CBUSHv(model)
        self.bushes = Bushes(model)

        #self.conm1 = CONM1v(model)
        self.conm2 = CONM2v(model)
        #self.cmass1 = CMASS1v(model)
        #self.cmass2 = CMASS2v(model)
        #self.cmass3 = CMASS3v(model)
        #self.cmass4 = CMASS4v(model)
        self.masses2 = Masses(model)

        self.crod = CRODv(model)
        self.conrod = CONRODv(model)
        self.ctube = CTUBEv(model)
        self.rods = Rods(model)

        self.cbar = CBARv(model)
        self.bars = Bars(model)

        self.cbeam = CBEAMv(model)
        self.beams = Beams(model)

        self.ctria3 = CTRIA3v(model)
        self.cquad4 = CQUAD4v(model)
        self.ctria6 = CTRIA6v(model)
        self.cquad8 = CQUAD8v(model)
        self.cquad = CQUADv(model)
        self.cquadr = CQUADRv(model)
        self.ctriar = CTRIARv(model)

        self.shells = Shells(model)
        #self.pshell = PSHELLv(model)  # TODO: temp

        self.cshear = CSHEARv(model)
        self.shears = Shears(model)

        #self.ctriax = CTRIA3v(model)   # TODO: temp
        #self.cquadx = CTRIA3v(model)   # TODO: temp
        #self.ctriax6 = CTRIA3v(model)  # TODO: temp
        #self.cquadx8 = CTRIA3v(model)  # TODO: temp

        self.ctetra4 = CTETRA4v(model)
        self.ctetra10 = CTETRA10v(model)
        self.chexa8 = CHEXA8v(model)
        self.chexa20 = CHEXA20v(model)
        self.cpenta6 = CPENTA6v(model)
        self.cpenta15 = CPENTA15v(model)
        self.cpyram5 = CPYRAM5v(model)
        self.cpyram13 = CPYRAM13v(model)
        self.solids = Solids(model)

        self.elements2 = Elements(model)  # TODO: change this name


        self.sload = SLOADv(model)
        self.grav = GRAVv(model)
        self.force = FORCEv(model)
        self.force1 = FORCE1v(model)
        self.force2 = FORCE2v(model)
        self.pload = PLOADv(model)
        self.pload1 = PLOAD1v(model)
        self.pload2 = PLOAD2v(model)
        self.pload4 = PLOAD4v(model)

        self.moment = MOMENTv(model)
        self.moment1 = MOMENT1v(model)
        self.moment2 = MOMENT2v(model)

        self.spcd = SPCDv(model)
        self.temp = TEMPv(model)
        self.tempd = TEMPDv(model)

        self.load_combinations = {}
        #def lseqi():
            #return LSEQv(model)
        #self.lseqs = defaultdict(lseqi)
        self.lseq = LSEQv(model)
        self.loads = Loads(model)
        #add_superelements_from_deck_lines(self, BDF, additional_deck_lines)

    def uncross_reference(self) -> None:
        pass
    def _prepare_grid(self, card, card_obj, comment=''):
        self.grid.add_card(card_obj, comment=comment)

    def _add_node_object(self, node, allow_overwrites=False):
        raise AttributeError("'BDF' object has no attribute '_add_node_object'")
    def _add_element_object(self, elem, allow_overwrites=False):
        key = elem.eid
        assert key > 0, 'eid=%s must be positive; elem=\n%s' % (key, elem)
        #if key in self.shells.eids:
        if key in self.elements and not allow_overwrites:
            if not elem == self.elements[key]:
                self._duplicate_elements.append(elem)
                if self._stop_on_duplicate_error:
                    self.pop_parse_errors()
        else:
            self.elements[key] = elem
            self._type_to_id_map[elem.type].append(key)

    #def _prepare_conm1(self, card, card_obj, comment=''):
        #self.conm1.add_card(card_obj, comment=comment)
    def _prepare_conm2(self, card, card_obj, comment=''):
        self.conm2.add_card(card_obj, comment=comment)
    #def _prepare_cmass1(self, card, card_obj, comment=''):
        #self.cmass1.add_card(card_obj, comment=comment)
    #def _prepare_cmass2(self, card, card_obj, comment=''):
        #self.cmass2.add_card(card_obj, comment=comment)
    #def _prepare_cmass3(self, card, card_obj, comment=''):
        #self.cmass3.add_card(card_obj, comment=comment)
    #def _prepare_cmass4(self, card, card_obj, comment=''):
        #self.cmass4.add_card(card_obj, comment=comment)

    def _prepare_celas1(self, card, card_obj, comment=''):
        self.celas1.add_card(card_obj, comment=comment)
    def _prepare_celas2(self, card, card_obj, comment=''):
        self.celas2.add_card(card_obj, comment=comment)
    def _prepare_celas3(self, card, card_obj, comment=''):
        self.celas3.add_card(card_obj, comment=comment)
    def _prepare_celas4(self, card, card_obj, comment=''):
        self.celas4.add_card(card_obj, comment=comment)

    def _prepare_cdamp1(self, card, card_obj, comment=''):
        self.cdamp1.add_card(card_obj, comment=comment)
    def _prepare_cdamp2(self, card, card_obj, comment=''):
        self.cdamp2.add_card(card_obj, comment=comment)
    def _prepare_cdamp3(self, card, card_obj, comment=''):
        self.cdamp3.add_card(card_obj, comment=comment)
    def _prepare_cdamp4(self, card, card_obj, comment=''):
        self.cdamp4.add_card(card_obj, comment=comment)
    #def _prepare_cdamp5(self, card, card_obj, comment=''):
        #self.cdamp5.add_card(card_obj, comment=comment)

    def _prepare_cvisc(self, card, card_obj, comment=''):
        self.cvisc.add_card(card_obj, comment=comment)
    def _prepare_plotel(self, card, card_obj, comment=''):
        self.plotel.add_card(card_obj, comment=comment)
    def _prepare_cbush(self, card, card_obj, comment=''):
        self.cbush.add_card(card_obj, comment=comment)
    #def _prepare_cbush1d(self, card, card_obj, comment=''):
        #self.cbush1d.add_card(card_obj, comment=comment)
    #def _prepare_cbush2d(self, card, card_obj, comment=''):
        #self.cbush2d.add_card(card_obj, comment=comment)

    def _prepare_conrod(self, card, card_obj, comment=''):
        self.conrod.add_card(card_obj, comment=comment)
    def _prepare_crod(self, card, card_obj, comment=''):
        self.crod.add_card(card_obj, comment=comment)
    def _prepare_ctube(self, card, card_obj, comment=''):
        self.ctube.add_card(card_obj, comment=comment)

    def _prepare_cbar(self, card, card_obj, comment=''):
        self.cbar.add_card(card_obj, comment=comment)
    def _prepare_cbeam(self, card, card_obj, comment=''):
        self.cbeam.add_card(card_obj, comment=comment)

    def _prepare_cquad4(self, card, card_obj, comment=''):
        """adds a CQUAD4"""
        self.cquad4.add_card(card_obj, comment=comment)
    def _prepare_ctria3(self, card, card_obj, comment=''):
        """adds a CTRIA3"""
        self.ctria3.add_card(card_obj, comment=comment)
    def _prepare_ctria6(self, card, card_obj, comment=''):
        """adds a CTRIA6"""
        self.ctria6.add_card(card_obj, comment=comment)
    def _prepare_cquad8(self, card, card_obj, comment=''):
        """adds a CQUAD8"""
        self.cquad8.add_card(card_obj, comment=comment)
    def _prepare_cquad(self, card, card_obj, comment=''):
        """adds a CQUAD"""
        self.cquad.add_card(card_obj, comment=comment)
    def _prepare_cquadr(self, card, card_obj, comment=''):
        """adds a CQUADR"""
        self.cquadr.add_card(card_obj, comment=comment)
    def _prepare_ctriar(self, card, card_obj, comment=''):
        """adds a CTRIAR"""
        self.ctriar.add_card(card_obj, comment=comment)

    def _prepare_cshear(self, card, card_obj, comment=''):
        """adds a CSHEAR"""
        self.cshear.add_card(card_obj, comment=comment)

    def _prepare_ctetra(self, card, card_obj, comment=''):
        """adds a CTETRA4/CTETRA10"""
        if len(card_obj) == 7:
            self.ctetra4.add_card(card_obj, comment=comment)
        else:
            self.ctetra10.add_card(card_obj, comment=comment)

    def _prepare_cpenta(self, card, card_obj, comment=''):
        """adds a CPENTA6/CPENTA15"""
        if len(card_obj) == 9:
            self.cpenta6.add_card(card_obj, comment=comment)
        else:
            self.cpenta15.add_card(card_obj, comment=comment)

    def _prepare_chexa(self, card, card_obj, comment=''):
        """adds a CHEXA8/CHEXA20"""
        if len(card_obj) == 11:
            self.chexa8.add_card(card_obj, comment=comment)
        else:
            self.chexa20.add_card(card_obj, comment=comment)

    def _prepare_cpyram(self, card, card_obj, comment=''):
        """adds a CPYRAM5/CPYRAM13"""
        if len(card_obj) == 8:
            self.cpyram5.add_card(card_obj, comment=comment)
        else:
            self.cpyram13.add_card(card_obj, comment=comment)

    def _prepare_load(self, card, card_obj, comment=''):
        load = LOAD.add_card(card_obj, comment=comment)
        key = load.sid
        assert key not in self.load_combinations
        self.load_combinations[key] = load

    def _prepare_lseq(self, card, card_obj, comment=''):
        #self.lseq.add_card(card_obj, comment=comment)
        #lseq = LSEQ.add_card(card_obj, comment=comment)
        #key = lseq.sid
        #assert key not in self.lseqs
        self.lseq.add_card(card_obj, comment=comment)

    def _prepare_grav(self, card, card_obj, comment=''):
        self.grav.add_card(card_obj, comment=comment)
    def _prepare_accel(self, card, card_obj, comment=''):
        if self.card_count['ACCEL'] == 1:
            self.log.warning('skipping %s' % str(card))
        #self.accel.add_card(card_obj, comment=comment)
    def _prepare_accel1(self, card, card_obj, comment=''):
        if self.card_count['ACCEL1'] == 1:
            self.log.warning('skipping %s' % str(card))
        #self.accel1.add_card(card_obj, comment=comment)
    def _prepare_sload(self, card, card_obj, comment=''):
        self.sload.add_card(card_obj, comment=comment)
    def _prepare_force(self, card, card_obj, comment=''):
        self.force.add_card(card_obj, comment=comment)
    def _prepare_force1(self, card, card_obj, comment=''):
        self.force1.add_card(card_obj, comment=comment)
    def _prepare_force2(self, card, card_obj, comment=''):
        self.force2.add_card(card_obj, comment=comment)
    def _prepare_moment(self, card, card_obj, comment=''):
        self.moment.add_card(card_obj, comment=comment)
    def _prepare_moment1(self, card, card_obj, comment=''):
        self.moment1.add_card(card_obj, comment=comment)
    def _prepare_moment2(self, card, card_obj, comment=''):
        self.moment2.add_card(card_obj, comment=comment)

    def _prepare_pload(self, card, card_obj, comment=''):
        self.pload.add_card(card_obj, comment=comment)
    def _prepare_pload1(self, card, card_obj, comment=''):
        self.pload1.add_card(card_obj, comment=comment)
    def _prepare_pload2(self, card, card_obj, comment=''):
        self.pload2.add_card(card_obj, comment=comment)
    def _prepare_pload4(self, card, card_obj, comment=''):
        self.pload4.add_card(card_obj, comment=comment)

    def _prepare_spcd(self, card, card_obj, comment=''):
        self.spcd.add_card(card_obj, comment=comment)
    def _prepare_temp(self, card, card_obj, comment=''):
        self.temp.add_card(card_obj, comment=comment)
    def _prepare_tempd(self, card, card_obj, comment=''):
        self.tempd.add_card(card_obj, comment=comment)

    def _prepare_ploadx1(self, card, card_obj, comment=''):
        if self.card_count['PLOADX1'] == 1:
            self.log.warning('skipping %s' % str(card))
    def _prepare_qvol(self, card, card_obj, comment=''):
        if self.card_count['QVOL'] == 1:
            self.log.warning('skipping %s' % str(card))
    def _prepare_qhbdy(self, card, card_obj, comment=''):
        if self.card_count['QHBDY'] == 1:
            self.log.warning('skipping %s' % str(card))
    def _prepare_rforce(self, card, card_obj, comment=''):
        if self.card_count['RFORCE'] == 1:
            self.log.warning('skipping %s' % str(card))
    def _prepare_rforce1(self, card, card_obj, comment=''):
        if self.card_count['RFORCE1'] == 1:
            self.log.warning('skipping %s' % str(card))
    def _prepare_qbdy1(self, card, card_obj, comment=''):
        if self.card_count['QBDY1'] == 1:
            self.log.warning('skipping %s' % str(card))
    def _prepare_qbdy2(self, card, card_obj, comment=''):
        if self.card_count['QBDY2'] == 1:
            self.log.warning('skipping %s' % str(card))
    def _prepare_qbdy3(self, card, card_obj, comment=''):
        if self.card_count['QBDY3'] == 1:
            self.log.warning('skipping %s' % str(card))
    #def _prepare_gmload(self, card, card_obj, comment=''):
        #if self.card_count['GMLOAD'] == 1:
            #self.log.warning('skipping %s' % str(card))
    def _prepare_loadcyn(self, card, card_obj, comment=''):
        if self.card_count['LOADCYN'] == 1:
            self.log.warning('skipping %s' % str(card))
    def _prepare_presax(self, card, card_obj, comment=''):
        if self.card_count['PRESAX'] == 1:
            self.log.warning('skipping %s' % str(card))

    def _update_card_parser(self):
        pass


    #def add_grid(self, nid, xyz, cp=0, cd=0, ps='', seid=0, comment=''):
        #pass
    #def cross_reference(self, xref=True, xref_nodes=True, xref_elements=True,
                       #xref_nodes_with_elements=True,
                       #xref_properties=True,
                       #xref_masses=True,
                       #xref_materials=True,
                       #xref_loads=True,
                       #xref_constraints=True,
                       #xref_aero=True, xref_sets=True,
                       #xref_optimization=True):
        #self.grid.cross_reference(self)

#def _add_node_object(self, node: Any, allow_overwrites: bool=False) -> None:
    #"""adds a GRID card"""
    #key = node.nid
    #if key in self.nodes and not allow_overwrites:
        #if not node == self.nodes[key]:
            #assert node.nid not in self.nodes, 'nid=%s\nold_node=\n%snew_node=\n%s' % (node.nid, self.nodes[key], node)
        #else:
            ##print('GRID was duplicated...nid=%s; node=\n%s' % (key, node))
            #pass
    #else:
        #assert key > 0, 'nid=%s node=%s' % (key, node)
        #self.node.add_[key] = node
        #self._type_to_id_map[node.type].append(key)

    @property
    def grids(self):
        raise NotImplementedError('grids is removed; use grid')
        #return self.grid

    @grids.setter
    def grids(self, grid):
        key = grid.nid
        self.grid.update(grid)
        #self.grid.add_grid(grid.nid, grid.xyz, cp=grid.cp, cd=grid.cd,
                           #ps=grid.ps, seid=grid.seid, comment=grid.comment)

    #def xyz_cid0(self):
        #return self.xyz_cid0


    #def _write_header(self, bdf_file, encoding):
        #"""
        #Writes the executive and case control decks.
        #"""
        #if self.punch is None:
            ## writing a mesh without using read_bdf
            #if self.executive_control_lines or self.case_control_deck:
                #self.punch = False
            #else:
                #self.punch = True

        #if self.nastran_format:
            #bdf_file.write('$pyNastran: version=%s\n' % self.nastran_format)
            #bdf_file.write('$pyNastran: punch=%s\n' % self.punch)
            #bdf_file.write('$pyNastran: encoding=%s\n' % encoding)
            #bdf_file.write('$pyNastran: nnodes=%s\n' % self.grid.n)
            #bdf_file.write('$pyNastran: nelements=%s\n' % len(self.elements))

        #if not self.punch:
            #self._write_executive_control_deck(bdf_file)
            #self._write_case_control_deck(bdf_file)

    #def _write_nodes(self, bdf_file: TextIO, size: int=8, is_double: bool=False) -> None:
        #"""
        #Writes the NODE-type cards
        #"""
        #BDF_._write_nodes(self, bdf_file, size=size, is_double=is_double)

    def _write_grids(self, bdf_file: TextIO,
                     size: int=8, is_double: bool=False,
                     is_long_ids: Optional[bool]=None) -> None:
        """Writes the GRID-type cards"""
        self.nodes.write_card(size=size, is_double=is_double, bdf_file=bdf_file)

    def _write_elements_interspersed(self, bdf_file: TextIO,
                                     size: int=8, is_double: bool=False,
                                     is_long_ids: Optional[bool]=None) -> None:
        """spoofed method"""
        self._write_elements(bdf_file, size=size, is_double=is_double, is_long_ids=is_long_ids)
        self._write_properties(bdf_file, size=size, is_double=is_double, is_long_ids=is_long_ids)

    def _write_elements(self, bdf_file: TextIO, size: int=8, is_double: bool=False,
                        is_long_ids: Optional[bool]=None) -> None:
        """Writes the elements in a sorted order"""
        if self.elements:
            bdf_file.write('$ELEMENTS\n')
            if self.is_long_ids:
                for (eid, element) in sorted(self.elements.items()):
                    bdf_file.write(element.write_card_16(is_double))
            else:
                for (eid, element) in sorted(self.elements.items()):
                    try:
                        bdf_file.write(element.write_card(size, is_double))
                    except Exception:
                        print('failed printing element...'
                              'type=%s eid=%s' % (element.type, eid))
                        raise
        self.elements2.write_card(size, is_double, bdf_file)
        #bdf_file.write(self.shells.write_card(size, is_double))
        #bdf_file.write(self.solids.write_card(size, is_double))
        if self.ao_element_flags:
            for (eid, element) in sorted(self.ao_element_flags.items()):
                bdf_file.write(element.write_card(size, is_double))
        self._write_nsm(bdf_file, size, is_double, is_long_ids=is_long_ids)

    #def _write_loads(self, bdf_file: TextIO, size: int=8, is_double: bool=False) -> None:
        #"""Writes the loads in a sorted order"""
        #BDF_._write_loads(self, bdf_file, size=size, is_double=is_double)
        ##for key, loadi in sorted(self.loads):
            ##bdf_file.write(loadi.write_card(size=size, is_double=is_double))

    def _write_loads(self, bdf_file: TextIO, size: int=8, is_double: bool=False,
                     is_long_ids: Optional[bool]=None) -> None:
        """Writes the load cards sorted by ID"""
        if self.loads or self.tempds:
            #msg = ['$LOADS\n']
            try:
                self.loads.write_card(size=size, is_double=is_double, bdf_file=bdf_file)
            except TypeError:
                print(self.loads)
                raise

            try:
                for key, load_combinations in sorted(self.load_combinations.items()):
                    bdf_file.write(load_combinations.write_card(size=size, is_double=is_double))
                    #for load_combination in load_combinations:
                        #bdf_file.write(load_combination.write_card(size=size, is_double=is_double))
            except(TypeError, AttributeError):
                print(self.load_combinations)
                raise

        assert len(self.tempds) == 0, self.tempds
        self._write_dloads(bdf_file, size=size, is_double=is_double, is_long_ids=is_long_ids)

    def get_displacement_index_xyz_cp_cd(self, fdtype: str='float64', idtype: str='int32',
                                         sort_ids: bool=True) -> Any:
        """
        Get index and transformation matricies for nodes with
        their output in coordinate systems other than the global.
        Used in combination with ``OP2.transform_displacements_to_global``

        Parameters
        ----------
        fdtype : str
            the type of xyz_cp
        idtype : str
            the type of nid_cp_cd

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
        return self.nodes.get_displacement_index_xyz_cp_cd(
            fdtype=fdtype, idtype=idtype)

    def get_node_index(self, nids, allow0=False):
        """maps the requested nodes to their desired index in the array"""
        i = self.nodes.get_node_index(nids, allow0=allow0)
        return i

    def validate(self):
        pass
    #def get_bdf_stats(self, return_type='string'):
        #pass

def read_bdf(bdf_filename: Optional[str]=None, validate: bool=True, xref: bool=True, punch: bool=False,
             skip_cards: list[str]=None,
             encoding: Optional[str]=None, log: Optional[SimpleLogger]=None, debug: bool=True,
             mode: str='msc'):
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
    skip_cards : list[str]; default=None
        None : include all cards
        list of cards to skip
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
    if skip_cards:
        model.disable_cards(skip_cards)
    model.read_bdf(bdf_filename=bdf_filename, validate=validate,
                   xref=False, punch=punch, read_includes=True, encoding=encoding)
    model.elements2.make_current()
    return model
