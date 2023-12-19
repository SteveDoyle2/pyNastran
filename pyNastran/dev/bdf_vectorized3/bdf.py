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
import os
import sys
from copy import deepcopy
from collections import Counter
from io import StringIO, IOBase
from pathlib import PurePath
from functools import partial
from collections import defaultdict
import traceback

from typing import (
    Sequence, Optional, Union, Any, TYPE_CHECKING)
from pickle import load, dump, dumps  # type: ignore

import numpy as np  # type: ignore
from cpylog import get_logger2

from pyNastran.utils import object_attributes, check_path, PathLike
from pyNastran.bdf.utils import parse_patran_syntax
from pyNastran.bdf.bdf_interface.utils import (
    _parse_pynastran_header, to_fields, parse_executive_control_deck,
    fill_dmigs, _get_card_name, _parse_dynamic_syntax,
)
from pyNastran.bdf.bdf_interface.add_card import CARD_MAP
from pyNastran.bdf.bdf_interface.replication import (
    to_fields_replication, get_nrepeats, int_replication, float_replication,
    _field, repeat_cards)

from pyNastran.bdf.field_writer_8 import print_card_8
from pyNastran.bdf.field_writer_16 import print_card_16, print_field_16

from pyNastran.bdf.cards.base_card import _format_comment
from pyNastran.bdf.cards.utils import wipe_empty_fields

#from .write_path import write_include
from pyNastran.bdf.bdf_interface.assign_type import (
    integer, integer_or_string, string)
from pyNastran.bdf.bdf_interface.model_group import ModelGroup
from pyNastran.bdf.cards.coordinate_systems import transform_coords_vectorized

#-------------------------------
from pyNastran.dev.bdf_vectorized3.bdf_interface.h5_pytables.h5_geometry import read_h5_geometry
from .cards.elements.bar import BAROR
from .cards.elements.beam import BEAMOR

#from pyNastran.bdf.cards.elements.elements import CRAC2D, CRAC3D
#from pyNastran.bdf.cards.properties.properties import PRAC2D, PRAC3D
#from pyNastran.bdf.cards.properties.solid import PIHEX
#from pyNastran.bdf.cards.cyclic import CYAX, CYJOIN
#from pyNastran.bdf.cards.msgmesh import CGEN

#from pyNastran.bdf.cards.elements.solid import (
    #CIHEX1, CIHEX2, CHEXA1, CHEXA2,
#)
#from pyNastran.bdf.cards.elements.rigid import RSPLINE, RSSCON

#from pyNastran.bdf.cards.axisymmetric.axisymmetric import (
    #AXIF, RINGFL,
    #AXIC, RINGAX, POINTAX, CCONEAX, PCONEAX)
#from pyNastran.bdf.cards.axisymmetric.loads import PLOADX1, FORCEAX, PRESAX, TEMPAX
#from pyNastran.bdf.cards.elements.shell_nasa95 import (
    #CTRSHL, CQUAD1, PQUAD1)

#from .cards.properties.shell import PTRSHL
#from .cards.elements.acoustic import (
    #CHACAB, CHACBR, PACABS, PAABSF, PACBAR, ACMODL)
#from .cards.elements.bush import CBUSH2D
#from .cards.properties.bush import PBUSH_OPTISTRUCT
#from .cards.elements.bars import CBEAM3, CBEND
#from .cards.elements.beam import CBEAM, BEAMOR
#from .cards.properties.bars import PBRSECT, PBEAM3
#from .cards.properties.beam import PBCOMP, PBMSECT
## CMASS5
#from pyNastran.bdf.cards.constraints import SPCAX, SESUP, GMSPC
#from .cards.coordinate_systems import (#CORD3G,
                                       #transform_coords_vectorized,
                                       #CORDx)
#from .cards.coordinate_systems.msgmesh import CGEN, GMCORD, GMLOAD
from .cards.deqatn import DEQATN
from pyNastran.bdf.cards.dynamic import (
    FREQ, FREQ1, FREQ2, FREQ3, FREQ4, FREQ5,
    TSTEP, TSTEP1, TSTEPNL, NLPARM, NLPCI, ROTORG, ROTORD,
)
#from .cards.loads.loads import (
    #LOADCYN, LOADCYH)
#from .cards.loads.dloads import ACSRCE
#from .cards.loads.static_loads import CLOAD

#from .cards.loads.random_loads import RANDT1

#from .cards.materials import (MAT3D,
                              #MATG, MATHE, MATEV,
                              #CREEP, EQUIV)
from pyNastran.bdf.cards.materials import NXSTRAT
#from .cards.material_deps import (
    #MATT2, MATT3, MATT4, MATT5, MATT9, MATS1)

from pyNastran.bdf.cards.methods import EIGB, EIGC, EIGR, EIGP, EIGRL
from .cards.grid import GRDSET
#from pyNastran.bdf.cards.nodes import GRDSET # SEQGP, GRIDB
from pyNastran.bdf.cards.aero.aero import MONDSP1
from pyNastran.bdf.cards.aero.static_loads import AEROS, DIVERG
from pyNastran.bdf.cards.aero.dynamic_loads import AERO, MKAERO1, MKAERO2
from pyNastran.bdf.cards.optimization import DOPTPRM
    #TOPVAR,
    #DRESP3,
    #DSCREEN)
#from .cards.optimization_nx import (
    #DVTREL1, GROUP, DMNCON,
#)
#from .cards.superelements import (
    #RELEASE, SEBNDRY, SEBULK, SECONCT, SEELT, SEEXCLD,
    #SELABEL, SELOAD, SELOC, SEMPLN, SENQSET, SETREE,
    #CSUPER, CSUPEXT,
#)
#from .cards.bdf_sets import (
    #SET2,
    #SEUSET, SEUSET1
    #SEQSEP
#)
from pyNastran.bdf.cards.params import PARAM, PARAM_MYSTRAN, PARAM_NASA95, MDLPRM
from pyNastran.bdf.cards.dmig import DMIG, DMI, DMIJ, DMIK, DMIJI, DMIG_UACCEL, DTI, DTI_UNITS, DMIAX
#from .bdf.cards.elements.thermal import BDYOR
#from .cards.thermal.loads import (TEMPB3, TEMPRB)
from pyNastran.dev.bdf_vectorized3.cards.elements.thermal import (
    CHBDYE, CHBDYG, CHBDYP, PCONV, PCONVM,
    PHBDY, CONV, CONVM, # TEMPBC,
)
#from .cards.thermal.radiation import RADM, RADBC, RADCAV, RADLST, RADMTX, VIEW, VIEW3D
from pyNastran.bdf.cards.bdf_tables import (
    TABLED1, TABLED2, TABLED3, TABLED4,
    TABLEM1, TABLEM2, TABLEM3, TABLEM4,
    TABLES1, # TABDMP1,
    TABLEST, TABLEH1, # TABLEHT,
    #TABRND1, TABRNDG,
    DTABLE,
)
from pyNastran.bdf.cards.contact import (
    BCRPARA, BCPARA, BCTPARA, BLSEG, BCTPARM, BCBODY)
#from .cards.parametric.geometry import PSET, PVAL, FEEDGE, FEFACE, GMCURV, GMSURF

from pyNastran.bdf.case_control_deck import CaseControlDeck, Subcase

from pyNastran.dev.bdf_vectorized3.bdf_interface.bdf_attributes import BDFAttributes
from pyNastran.dev.bdf_vectorized3.bdf_interface.add_card import AddCards
from pyNastran.dev.bdf_vectorized3.bdf_interface.write_mesh import WriteMesh
from pyNastran.bdf.bdf_interface.bdf_card import BDFCard
#from pyNastran.bdf.bdf_interface.write_mesh_file import WriteMeshs
#from pyNastran.bdf.bdf_interface.uncross_reference import UnXrefMesh
#from pyNastran.bdf.bdf_interface.verify_validate import verify_bdf, validate_bdf
#from pyNastran.bdf.bdf_interface.stats import get_bdf_stats

from pyNastran.bdf.errors import (CrossReferenceError, DuplicateIDsError,
                                  CardParseSyntaxError, UnsupportedCard, DisabledCardError,
                                  SuperelementFlagError, ReplicationError)
from pyNastran.bdf.bdf_interface.pybdf import (
    BDFInputPy, _clean_comment, _clean_comment_bulk, EXECUTIVE_CASE_SPACES)

#from .bdf_interface.add_card import CARD_MAP
if TYPE_CHECKING:  # pragma: no cover
    from cpylog import SimpleLogger

REMOVED_CARDS = {
    'ADAPT',
    'PVAL', 'GMCURV', 'GMSURF', 'FEEDGE', 'FEFACE', 'GMSPC', 'GMLOAD',
    #'OUTPUT'
    #'OUTRCV'
    'GMBNDS', 'GMINTS', 'PINTS',
    'GMBNDC', 'GMINTC', 'PINTC'
    'CGEN', 'EGRID', 'GRIDG', 'SPCG',
}

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
    # msgmesh
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
    'RJOINT', 'RTRPLT', 'RTRPLT1', 'DYNRED',

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
    'MATPOR', # 'MATDMG',
    'MAT2F', 'MAT8F',

    ## loads
    'FORCDST', 'PLOADX', 'PLOADB3', 'PLOADG', 'QBDY4', 'SLOADN1',
    'TEMPP1', 'TTEMP', 'RGYRO', 'PLOADE1', 'RCROSSC', 'LOADCYT',
    'TEMPN1',

    ## boundaries
    'BNDFREE', 'BNDFREE1',
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
    'ACLOAD', 'PBARN1',
    'TABL3D1', 'AEGRID', 'AEQUAD4', 'SPBLND1', 'CONCTL',
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
OBJ_CARDS = {
    'PARAM', 'MDLPRM', 'TSTEP', 'TSTEP1', 'TSTEPNL',
    'NLPCI', 'NLPARM',
    'EIGR', 'EIGRL', 'EIGB',
    'EIGC', 'EIGP', 'NXSTRAT',
    'FREQ', 'FREQ1', 'FREQ2', 'FREQ3', 'FREQ4', 'FREQ5',
    # aero
    'AERO', 'AEROS',
    'MKAERO1', 'MKAERO2',
    # contact
    'BCTPARA', 'BCTPARM',
    # optimization
    'DOPTPRM',
    # tables
    'TABLED1', 'TABLED2', 'TABLED3', 'TABLED4', 'TABLED5',
    'TABLEM1', 'TABLEM2', 'TABLEM3', 'TABLEM4',
    'TABLES1', 'TABLEST', 'TABLEH1', 'TABLEHT',
    'DTABLE',
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


#class BDF_(BDFMethods, GetCard, AddCards, WriteMeshs, UnXrefMesh):
class BDF(AddCards, WriteMesh): # BDFAttributes
    """NASTRAN BDF Reader/Writer/Editor class."""
    _properties = ['nid_map', 'wtmass', 'type_slot_str'] + [
        'nastran_format', 'is_long_ids', 'sol', 'subcases',
        'nnodes', 'node_ids', 'point_ids', 'npoints',
        'nelements', 'element_ids', 'nproperties', 'property_ids',
        'nmaterials', 'material_ids', 'ncoords', 'coord_ids',
        'ncaeros', 'caero_ids', 'wtmass', 'nid_map',
        #'dmigs', 'dmijs', 'dmiks', 'dmijis', 'dtis', 'dmis',
    ]

    #: required for sphinx bug
    #: http://stackoverflow.com/questions/11208997/autoclass-and-instance-attributes
    #__slots__ = ['_is_dynamic_syntax']
    def __init__(self, debug: Union[str, bool, None]=True,
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
        self.log = get_logger2(log=log, debug=debug)
        super().__init__()
        #AddCards.__init__(self)
        #WriteMesh.__init__(self)
        #BDFAttributes.__init__(self)
        assert debug in [True, False, None], f'debug={debug!r}'
        self.echo = False
        self.read_includes = True
        self._remove_disabled_cards = False

        self.is_lax_parser = False

        # file management parameters
        self.active_filenames: list[str] = []
        self.active_filename: Optional[str] = ''
        self.include_dir = ''
        self.dumplines = False

        # list of all read in cards - useful in determining if entire BDF
        # was read & really useful in debugging
        self.card_count: dict[str, int] = {}

        # stores the card_count of cards that have been rejected
        self.reject_count: dict[str, int] = {}

        # allows the BDF variables to be scoped properly (i think...)
        #GetCard.__init__(self)
        AddCards.__init__(self)
        #BDFMethods.__init__(self)
        WriteMesh.__init__(self)
        #UnXrefMesh.__init__(self)

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
            'PARAM', 'MDLPRM',

            ## nodes
            'GRID', 'GRDSET', 'SPOINT',
            'EPOINT', # 'SEQGP', 'GRIDB',

            # points
            'POINT',
            #'GRIDG'

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
            'CBUSH', 'CBUSH1D', # 'CBUSH2D',
            # dampers
            'CDAMP1', 'CDAMP2', 'CDAMP3', 'CDAMP4', 'CDAMP5', 'CVISC',
            #  misc
            'CFAST', 'CGAP', 'CSHEAR', 'GENEL',

            'CBAR', 'CBARAO', 'BAROR',
            'CROD', 'CTUBE', 'CONROD',
            'CBEAM', 'CBEAM3', 'CBEND', 'BEAMOR',
            'CTRIA3', 'CTRIA6', 'CTRIAR',
            'CQUAD4', 'CQUAD8', 'CQUADR', 'CQUAD',
            'CTRAX3', 'CTRAX6', 'CTRIAX', 'CTRIAX6', 'CQUADX', 'CQUADX4', 'CQUADX8',
            #'SNORM',

            'CPLSTN3', 'CPLSTN4', 'CPLSTN6', 'CPLSTN8', # plate strain
            'CPLSTS3', 'CPLSTS4', 'CPLSTS6', 'CPLSTS8', # plate stress

            # solids
            'CTETRA', 'CPYRAM', 'CPENTA', 'CHEXA',
            'CHEXCZ',

            # acoustic
            #'CHACAB',
            'CAABSF', # 'CHACBR',
            #'PACABS',
            'PAABSF', #'PACBAR', 'ACMODL',

            # crack
            #'CRAC2D', 'CRAC3D',

            ## rigid_elements
            'RBAR', 'RBAR1',
            'RBE1', 'RBE2', 'RBE3', 'RROD', # 'RSPLINE', 'RSSCON',

            ## plotels
            'PLOTEL', 'PLOTEL3', 'PLOTEL4', 'PLOTEL6', 'PLOTEL8',

            ## properties
            'PMASS',
            'PELAS', 'PGAP', 'PFAST',
            'PLPLANE', 'PPLANE',
            'PBUSH', 'PBUSH1D',
            'PDAMP', 'PDAMP5',
            'PROD', 'PBAR', 'PBARL', 'PTUBE',
            'PBEAM', 'PBEAML', 'PBCOMP', # 'PBRSECT',
            'PBEND',
            #'PBMSECT', # not fully supported
            #'PBEAM3',  # v1.3

            'PSHELL', 'PCOMP', 'PCOMPG',
            'PSHEAR', 'PSHLN1', 'PSHLN2',
            'PSOLID', 'PLSOLID', 'PVISC', # 'PRAC2D', 'PRAC3D',
            'PCOMPS', 'PCOMPLS',

            ## pdampt
            'PDAMPT',

            ## pelast
            'PELAST',

            ## pbusht
            'PBUSHT',

            ## creep_materials
            #'CREEP',

            ## materials
            'MAT1', 'MAT2', 'MAT3',
            'MAT8', 'MAT9',
            'MAT10', 'MAT11', # 'MAT3D',
            #'MATG',
            #'MATHE',
            'MATHP', 'MATORT', 'MAT10C', # 'MATEV',

            ## Material dependence - MATT1/MATT2/etc.
            'MATT1', 'MATT2', # 'MATT3', 'MATT4', 'MATT5',
            'MATT8', 'MATT9',
            'MATS1', #'MATS3',
            'MATS8',
            #'EQUIV', # testing only, should never be activated...

            ## nxstrats
            'NXSTRAT',

            ## thermal_materials
            'MAT4', 'MAT5',

            ## spcs
            'SPC', 'SPCADD', 'SPC1',
            'SPCOFF', 'SPCOFF1',
            'BNDFIX', 'BNDFIX1',
            'BNDFREE', 'BNDFREE1',

            ## mpcs
            'MPC', 'MPCADD',

            ## suport/suport1/se_suport
            'SUPORT', 'SUPORT1',  # suport
            #'SESUP',

            ## dloads
            'DLOAD',

            ## dload_entries
            'TLOAD1', 'TLOAD2', 'RLOAD1', 'RLOAD2',
            'QVECT',
            'ACSRCE', 'RANDPS', # 'RANDT1', # random

            ## loads
            'LOAD', # 'CLOAD',
            'LSEQ', # 'LOADCYN', 'LOADCYH',
            'SLOAD',
            'FORCE', 'FORCE1', 'FORCE2',
            'MOMENT', 'MOMENT1', 'MOMENT2',
            'GRAV', 'ACCEL', 'ACCEL1',
            'PLOAD', 'PLOAD1', 'PLOAD2', 'PLOAD4',
            'RFORCE', 'RFORCE1',
            'SPCD', 'DEFORM',

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
            'CAERO1', 'CAERO2', 'CAERO3', 'CAERO4', 'CAERO5', 'CAERO7', ## caeros
            'PAERO1', 'PAERO2', 'PAERO3', 'PAERO4', 'PAERO5', ## paeros
            #'AEFORCE', 'UXVEC', 'GUST2',

            'SPLINE1', 'SPLINE2', 'SPLINE3', 'SPLINE4', 'SPLINE5',  ## splines
            #'SPLINE6', 'SPLINE7',
            'TRIM', # 'TRIM2',  ## trims
            'CSSCHD',  ## csschds
            #'DIVERG',  ## divergs

            'MONPNT1', 'MONPNT2', 'MONPNT3', # 'MONDSP1', ## monitor_points

            ## coords
            'CORD1R', 'CORD1C', 'CORD1S',
            'CORD2R', 'CORD2C', 'CORD2S',

            # temperature cards
            'TEMP', 'TEMPD', # 'TEMPB3', 'TEMPAX',
            'QBDY1', 'QBDY2', 'QBDY3', 'QHBDY',
            'CHBDYE', 'CHBDYG', 'CHBDYP', 'BDYOR',
            'PCONV', 'PCONVM', 'PHBDY',
            'RADBC', 'CONV',
            'RADM', # 'VIEW', 'VIEW3D',  # TODO: not validated

            #'RADCAV', ## radcavs

            # ---- dynamic cards ---- #
            'DAREA',  ## darea
            'DPHASE',  ## dphase
            'DELAY',  ## delay
            'NLPARM',  ## nlparms
            #'ROTORG', 'ROTORD', ## rotors
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
            'DCONSTR', 'DESVAR', 'DDVAL', 'DRESP1', 'DRESP2', # 'DRESP3', 'TOPVAR',
            'DVCREL1', 'DVCREL2',
            'DVPREL1', 'DVPREL2',
            'DVMREL1', 'DVMREL2',
            'DOPTPRM',
            'DLINK', 'DCONADD', 'DVGRID',
            'DSCREEN',

            # nx optimization
            #'DVTREL1', 'GROUP', 'DMNCON',

            # sets
            'SET1', 'SET3',  ## sets
            #'SET2',

            'ASET', 'ASET1',  ## aset
            'OMIT', 'OMIT1',  ## omit
            'BSET', 'BSET1',  ## bset
            'CSET', 'CSET1',  ## cset
            'QSET', 'QSET1',  ## qset
            'USET', 'USET1',  ## uset

            'RADSET',  # radset

            # superelements
            #'SETREE', 'SENQSET', 'SEBULK', 'SEBNDRY', 'SEELT', 'SELOC', 'SEMPLN',
            #'SECONCT', 'SELABEL', 'SEEXCLD', 'CSUPER', 'CSUPEXT',
            #'SELOAD',

            # super-element sets
            'SESET',  ## se_set

            'SEBSET', 'SEBSET1',  ## se_bset
            'SECSET', 'SECSET1',  ## se_cset
            'SEQSET', 'SEQSET1',  ## se_qset
            #'SEUSET', 'SEUSET1',  ## se_uset
            'RELEASE',
            #'SEQSEP',

            #------------------------------------------------------------------
            # axixsymmetric
            #'CCONEAX', # element
            #'PCONEAX', # property
            #'AXIC', # axic
            #'AXIF', # axif

            #'FORCEAX', # loads
            #'PRESAX', # loads
            #'PLOADX1',
            #'SPCAX', # spcs

            #------------------------------------------------------------------
            ## parametric
            #'PSET', 'PVAL', 'GMCURV', 'GMSURF', 'FEEDGE', 'FEFACE',
            #'GMSPC',  # spcs

            # msgmesh
            #'GMLOAD',  # loads
            #'GMCORD',  # coords

            #  nastran 95
            #'CIHEX1', 'CIHEX2', 'CHEXA1', 'CHEXA2',
            #'CTRSHL', 'CQUAD1',
            #'PTRSHL', 'PQUAD1',
            #'PIHEX', # PQUAD4

            # ringfl
            #'RINGFL',
            # ringaxs
            #'RINGAX', 'POINTAX',
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
            #'TABDMP1',

            ## random_tables
            # PSD=func(freq); used by RANDPS card
            #'TABRND1',
            # gust for aeroelastic response; used by RANDPS card
            #'TABRNDG',

            # ???
            #'TABLEHT',
            'TABLEH1',

            #------------------------------------------------------------------
            #: methods
            'EIGB', 'EIGR', 'EIGRL',

            #: cMethods
            'EIGC', 'EIGP',

            # : modtrak
            'MODTRAK',

            #: contact
            #'BCBODY',  ## bcbody
            #'BCPARA',  ## bcpara

            # nx-contact
            #-----------
            # glue...
            'BGADD',  ## bgadds
            'BGSET',  ## bgset  - glue set (points to BSURF/BSURFS/BCPROP/BCPROPS)
            'BEDGE',

            # not glue...
            'BCRPARA',  ## bcrpara
            'BCTPARM', ## bctparm
            'BCTADD',  ## bctadd
            'BCTSET',   ## bctset  - contact_id to source/target_ids
            'BSURF',    ## bsurf   - source/target_id to shell eid
            'BSURFS',   ## bsurfs  - source/target_id to solid eid
            'BCPROP',   ## bcprop  - source/target_id to shell pid
            'BCPROPS',  ## bcprops - source/target_id to solid pid
            'BCONP', ## bconp
            'BLSEG', ## blseg - glue edge
            'BFRIC', ## bfric

            #'TEMPBC',
            #'RADMT',
            #'RADLST', 'RADMTX', #'RADBND',
            #'TEMPP1',
            #'TEMPRB',
            'CONVM',
            ## ???
            #'PANEL', 'SWLDPRM',
            #'CWELD', 'PWELD',
            #'CWELD', 'PWELD', 'PWSEAM', 'CWSEAM', 'CSEAM', 'PSEAM', 'DVSHAP', 'BNDGRID',
            #'CYSYM', 'CYJOIN',
            # 'DSCONS', 'DVAR', 'DVSET', 'DYNRED',

            # cyclic
            #'CYJOIN', 'CYAX',

            # other
            'INCLUDE',  # '='
            'ENDDATA',
        ]
        set_cards_to_read = set(cards_to_read)
        if len(cards_to_read) != len(set_cards_to_read):  # pragma: no cover
            bad_cards = [key for key, value in Counter(cards_to_read).items()
                         if value > 1]
            raise RuntimeError(f'duplicate cards in cards_to_read={bad_cards}')

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

    def get_h5attrs(self) -> list[str]:
        """helper method for dict_to_h5py"""
        attrs = self.object_attributes(mode='both', keys_to_skip=None)
        return attrs

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
        #try:
            #with h5py.File(hdf5_filename, 'w') as hdf5_file:
                ##self.log.info('starting export_hdf5_file of %r' % hdf5_filename)
                #self.export_hdf5_file(hdf5_file)
        #except OSError:
            #self.log.error(f'failed to export {hdf5_filename!r}')
            #raise

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
            ##self.log.info('starting load_hdf5_file of %r' % hdf5_filename)
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

    def save_obj(self, obj_filename: str='model.obj', unxref: bool=True) -> None:
        """Saves a pickleable object"""
        #del self.log
        #del self._card_parser, self._card_parser_prepare

        #try:
            #del self.log
        #except AttributeError:
            #pass
        #self.case_control_lines = str(self.case_control_deck).split('\n')
        #del self.case_control_deck

        with open(obj_filename, 'wb') as obj_file:
            dump(self, obj_file)

    def load_obj(self, obj_filename: str='model.obj') -> None:
        """Loads a pickleable object"""
        #del self.log
        #lines = print(self.case_control_deck)
        #self.case_control_lines = lines.split('\n')
        #del self.case_control_deck
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
            #'dmig', 'dmij', 'dmik', 'dmiji', 'dti', 'dmi',

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
                raise AttributeError(f'key={key!r} val={val}\nupdate ~line 1050 of bdf.py and '
                                     f'add the new key ({key})')

        self.case_control_deck = CaseControlDeck(self.case_control_lines, log=self.log)
        #self.log.debug('done loading!')
        for model in self.superelement_models.values():
            model.log = self.log

    #def replace_cards(self, replace_model) -> None:
        #"""
        #Replaces the common cards from the current (self) model from the
        #ones in the new replace_model.  The intention is that you're
        #going to replace things like PSHELLs and DESVARs from a pch file
        #in order to update your BDF with the optimized geometry.

        #.. todo:: only does a subset of cards.

        #Notes
        #-----
        #loads/spcs (not supported) are tricky because you can't replace
        #cards one-to-one...not sure what to do

        #"""
        #for nid, node in replace_model.nodes.items():
            #self.nodes[nid] = node
        #for eid, elem in self.elements.items():
            #self.elements[eid] = elem
        #for eid, elem in replace_model.rigid_elements.items():
            #self.rigid_elements[eid] = elem
        #for pid, prop in replace_model.properties.items():
            #self.properties[pid] = prop
        #for mid, mat in replace_model.materials.items():
            #self.materials[mid] = mat

        #for dvid, desvar in replace_model.desvars.items():
            #self.desvars[dvid] = desvar
        #for dvid, dvprel in replace_model.dvprels.items():
            #self.dvprels[dvid] = dvprel
        #for dvid, dvmrel in replace_model.dvmrels.items():
            #self.dvmrels[dvid] = dvmrel
        #for dvid, dvgrid in replace_model.dvgrids.items():
            #self.dvgrids[dvid] = dvgrid

    def disable_cards(self, cards: Sequence[str]) -> None:
        """
        Method for removing broken cards from the reader

        Parameters
        ----------
        cards : list[str]; Set[str]
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

    def enable_cards(self, cards: Sequence[str]) -> None:
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
        assert isinstance(nparse_errors, int), type(nparse_errors)
        assert isinstance(nxref_errors, int), type(nxref_errors)
        self._nparse_errors = nparse_errors
        self._nxref_errors = nxref_errors
        self._stop_on_parsing_error = stop_on_parsing_error
        self._stop_on_xref_error = stop_on_xref_error

    def validate(self) -> None:
        """runs some checks on the input data beyond just type checking"""
        pass
        #self.setup()
        #validate_bdf(self)
        self.log.warning('no validate')
    def _get_rigid(self):
        self.log.warning('no _get_rigid')
    def get_rigid_elements_with_node_ids(self, node_ids: list[int]):
        self.log.warning('no get_rigid_elements_with_node_ids')
    #def get_rigid_elements_with_node_ids(self, node_ids: list[int]):
        #self.log.warning('no get_rigid_elements_with_node_ids')
    #def get_rigid_elements_with_node_ids(self, node_ids: list[int]):
        #self.log.warning('no get_rigid_elements_with_node_ids')
    #def get_rigid_elements_with_node_ids(self, node_ids: list[int]):
        #self.log.warning('no get_rigid_elements_with_node_ids')

    def _verify_bdf(self, xref: Optional[bool]=None) -> None:
        """Cross reference verification method."""
        self.log.warning('skipping _verify_bdf')
        return
        if xref is None:
            xref = self._xref
        verify_bdf(self, xref)

    def include_zip(self, bdf_filename: Optional[str]=None,
                    encoding: Optional[str]=None,
                    make_ilines: bool=True) -> tuple[list[str], Any]:
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
        self._set_pybdf_attributes(obj, save_file_structure=False)
        return all_lines, ilines

    def _set_pybdf_attributes(self, obj: BDFInputPy,
                              save_file_structure: bool=False) -> None:
        """common method for all functions that use BDFInputPy"""
        # these are include line pairs
        #print(obj.include_lines)
        self.active_filenames = []
        self.reject_lines = []
        include_filenames = defaultdict(list)
        for ifile, include_lines_filename_pairs in obj.include_lines.items():
            assert len(include_lines_filename_pairs) > 0, include_lines_filename_pairs
            for include_lines, bdf_filename2 in include_lines_filename_pairs:
                #print(ifile, include_lines)
                include_filenames [ifile].append(bdf_filename2)
                if not save_file_structure and not obj.read_includes:
                    self.reject_lines += include_lines

        self.include_filenames: dict[int, list[str]] = dict(include_filenames)
        #print('-------------ssett (end)----------')
        self.active_filenames += obj.active_filenames
        self.active_filename = obj.active_filename
        self.include_dir = obj.include_dir

    def read_h5(self, h5_filename: PathLike):
        read_h5_geometry(self, h5_filename)

    def read_bdf(self, bdf_filename: Optional[PathLike]=None,
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
        out = obj.get_lines(bdf_filename, punch=self.punch, make_ilines=True)
        (system_lines,
         executive_control_lines,
         case_control_lines,
         bulk_data_lines, bulk_data_ilines,
         superelement_lines, superelement_ilines) = out
        self._set_pybdf_attributes(obj, save_file_structure)

        #assert system_lines == [], system_lines
        #assert executive_control_lines == [], executive_control_lines
        #assert case_control_lines == [], case_control_lines
        self.system_command_lines = system_lines
        self.executive_control_lines = executive_control_lines
        self.case_control_lines = case_control_lines

        sol, method, sol_iline, app = parse_executive_control_deck(executive_control_lines)
        self.app = app
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

        if self._remove_disabled_cards:
            all_cards = set(self.card_count.keys())
            union_cards = all_cards.intersection(REMOVED_CARDS)
            if union_cards:
                raise DisabledCardError(f'the following cards have been removed: {list(union_cards)}')

        self.cross_reference(run_geom_check=xref)
        self._xref = xref

        self.log.debug('---finished BDF.read_bdf of %s---' % self.bdf_filename)

    def _add_superelements(self, superelement_lines: list[str],
                           superelement_ilines: Any) -> None:  # pragma: no cover
        self.log.warning('_add_superelements should be overwritten')

    def _parse_all_cards(self, bulk_data_lines: list[str], bulk_data_ilines: Any) -> None:
        """creates and loads all the cards the bulk data section"""
        strict = True
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
        self._parse_cards(cards_list, cards_dict, card_count, strict=strict)

        if self.values_to_skip:
            for key, values in self.values_to_skip.items():
                dict_values = getattr(self, key)
                if not isinstance(dict_values, dict):
                    msg = f'{key!r} is an invalid type; only dictionaries are supported'
                    raise TypeError(msg)
                for value in values:
                    del dict_values[value]
            # TODO: redo get_card_ids_by_card_types & card_count

    def _read_bdf_helper(self, bdf_filename: Optional[PathLike], encoding: str,
                         punch: bool, read_includes: bool):
        """creates the file loading if bdf_filename is None"""
        #self.set_error_storage(nparse_errors=None, stop_on_parsing_error=True,
        #                       nxref_errors=None, stop_on_xref_error=True)
        if encoding is None:
            encoding = sys.getdefaultencoding()
        self._encoding = encoding

        self.read_includes = read_includes
        self.active_filenames = []

        # turns bdf_filename -> bdf_filename2
        if bdf_filename is None:
            from pyNastran.utils.gui_io import load_file_dialog
            wildcard_wx = "Nastran BDF (*.bdf; *.dat; *.nas; *.pch, *.ecd)|" \
                "*.bdf;*.dat;*.nas;*.pch;*.ecd|" \
                "All files (*.*)|*.*"
            wildcard_qt = "Nastran BDF (*.bdf *.dat *.nas *.pch *.ecd);;All files (*)"
            title = 'Please select a BDF/DAT/PCH/ECD to load'
            bdf_filename2 = load_file_dialog(title, wildcard_wx, wildcard_qt)[0]
            assert bdf_filename2 is not None, bdf_filename2

        elif isinstance(bdf_filename, (str, PurePath)):
            bdf_filename2 = bdf_filename
        elif isinstance(bdf_filename, (StringIO, IOBase)):
            self.bdf_filename = bdf_filename
            self.punch = punch
            return
        else:
            raise NotImplementedError(bdf_filename)

        #-------------------------------
        check_path(bdf_filename2, 'bdf_filename')
        ext = os.path.splitext(bdf_filename2)[1]
        if ext == '.pch':  # .. todo:: should this be removed???
            punch = True

        #: the active filename (string)
        self.bdf_filename = bdf_filename2

        #: is this a punch file (no executive control deck)
        self.punch = punch
        assert ext != '.op2', bdf_filename2

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

    def get_bdf_cards(self, bulk_data_lines: list[str],
                      bulk_data_ilines: Optional[Any]=None) -> tuple[Any, Any, Any]:
        """Parses the BDF lines into a list of card_lines"""
        if bulk_data_ilines is None:
            bulk_data_ilines = np.zeros((len(bulk_data_lines), 2), dtype='int32')

        cards_list: list[Any] = []
        cards_dict = defaultdict(list)  # type: dict[str, list[Any]]
        dict_cards = ['BAROR', 'BEAMOR']
        #cards = defaultdict(list)
        card_count = defaultdict(int)  # dict[str, int]
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
                strip_comment = comment.strip()
                if strip_comment.lower().startswith('group:'):
                    #'group: name="ULFuseCanardAtch MainFuseStruct Fixed Gridpoints"; nodes=1'
                    strip_comment2 = strip_comment.split(':', 1)[1].strip()
                    if ';' in strip_comment2:
                        group = ModelGroup.create_from_line(strip_comment2)
                        name = group.name
                        if name in self.model_groups:
                            og_group = self.model_groups[name]
                            #print(og_group)
                            og_group.union(group)
                            #print('->', og_group)
                            del og_group
                            continue
                        self.model_groups[name] = group
                    else:
                        self.log.warning(f'unknown group={strip_comment}')

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

    def update_card(self, card_name: str, ids: np.ndarray, ifield: int,
                    values: np.ndarray) -> None:
        """
        Updates a Nastran card based on standard Nastran optimization names

        Parameters
        ----------
        card_name : str
            the name of the card
            (e.g. GRID)
        ids : (n, ) int np.ndarray
            the unique 1-based index identifier for the card
            (e.g. the GRID id)
        ifield : int
            the index on the card
            (e.g. X on GRID card as an integer representing the field number)
        value : (n, ) int/float/str np.ndarray
            the value to assign (float/int/str)

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
        card = getattr(self, card_name.lower())

        #ids = [1, 2, 3]  # -> GRID 1, 2, 3
        #ifield = 2   # CP
        #ifield = 3   # X
        #ifield = 4   # Y
        #ifield = 5   # Z
        #values = [0.1, 0.2, 0.3]
        card.update_field(ifield, ids, values)
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
        if '=' in card_name:
            raise ReplicationError('unparsed replication format')

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
        if card_name in ['DEQATN', 'PBRSECT', 'PBMSECT', 'GMCURV', 'GMSURF', 'OUTPUT', 'ADAPT',
                         'MONDSP1']:
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
            #'=' : (Crash, None),
            '/' : (Crash, None),

            #'SETREE' : (SETREE, add_methods._add_setree_object),
            #'SENQSET' : (SENQSET, add_methods._add_senqset_object),
            #'SEBULK' : (SEBULK, add_methods._add_sebulk_object),
            #'RELEASE': (RELEASE, add_methods._add_release_object),
            #'SEBNDRY' : (SEBNDRY, add_methods._add_sebndry_object),
            #'SEELT' : (SEELT, add_methods._add_seelt_object),
            #'SELOC' : (SELOC, add_methods._add_seloc_object),
            #'SEMPLN' : (SEMPLN, add_methods._add_sempln_object),
            #'SECONCT' : (SECONCT, add_methods._add_seconct_object),
            #'SELABEL' : (SELABEL, add_methods._add_selabel_object),
            #'SEEXCLD' : (SEEXCLD, add_methods._add_seexcld_object),
            #'CSUPER' : (CSUPER, add_methods._add_csuper_object),
            #'CSUPEXT' : (CSUPEXT, add_methods._add_csupext_object),
            #'SELOAD' : (SELOAD, add_methods._add_seload_object),

            #'PANEL' : (Crash, None),

            # nx contact
            #'BCPARA' : (BCPARA, add_methods._add_bcpara_object),
            'BCTPARM' : (BCTPARM, add_methods._add_bctparm_object),
            #'BCBODY' : (BCBODY, add_methods._add_bcbody_object),

            # nx bolts
            'BOLT' : (Crash, None),
            'BOLTFOR' : (Crash, None),
            'BOLTLD' : (Crash, None),

            #'CBEAR', 'PBEAR', 'ROTORB',
            #'CBEAR' : (Crash, None),
            #'PBEAR' : (Crash, None),
            #'ROTORB' : (Crash, None),

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

            #'AEFORCE' : (Crash, None),
            #'UXVEC' : (Crash, None),
            #'GUST2' : (Crash, None),

            #'RADBND' : (Crash, None),

            # nodes
            #'RINGAX' : (RINGAX, add_methods._add_ringax_object),
            #'POINTAX' : (POINTAX, add_methods._add_ringax_object),
            #'POINT' : (POINT, add_methods._add_point_object),
            #'SEQGP' : (SEQGP, add_methods._add_seqgp_object),
            #'GRIDB' : (GRIDB, add_methods._add_gridb_object),

            'PARAM' : (PARAM, add_methods._add_param_object),
            'MDLPRM' : (MDLPRM, add_methods._add_mdlprm_object),

            # parametric
            #'PSET' : (PSET, add_methods._add_pset),
            #'PVAL' : (PVAL, add_methods._add_pval),
            #'GMCURV' : (GMCURV, add_methods._add_gmcurv),
            #'GMSURF' : (GMSURF, add_methods._add_gmsurf),
            #'FEFACE' : (FEFACE, add_methods._add_feface),
            #'FEEDGE' : (FEEDGE, add_methods._add_feedge),

            # msgmesh
            #'CGEN' : (CrashIgnore, None),
            #'GMCORD' : (GMCORD, add_methods._add_coord_object), # coords
            #'CGEN' : (CGEN, add_methods._add_element_object),   # elements
            #'GMLOAD' : (GMLOAD, add_methods._add_load_object),  # basic loads

            #'CBEAM3' : (CBEAM3, add_methods._add_element_object),
            #'PBEAM3' : (PBEAM3, add_methods._add_property_object),

            # nastran95
            #'CTRSHL' : (CTRSHL, add_methods._add_element_object),
            #'CQUAD1' : (CQUAD1, add_methods._add_element_object),
            #'PTRSHL' : (PTRSHL, add_methods._add_property_object),
            #'PQUAD1' : (PQUAD1, add_methods._add_property_object),
            #
            #'CIHEX1' : (CIHEX1, add_methods._add_element_object),
            #'CIHEX2' : (CIHEX2, add_methods._add_element_object),
            #'CHEXA1' : (CHEXA1, add_methods._add_element_object),
            #'CHEXA2' : (CHEXA2, add_methods._add_element_object),
            #'PIHEX' : (PIHEX, add_methods._add_property_object),

            #'CRAC2D' : (CRAC2D, add_methods._add_element_object),
            #'PRAC2D' : (PRAC2D, add_methods._add_property_object),

            #'CRAC3D' : (CRAC3D, add_methods._add_element_object),
            #'PRAC3D' : (PRAC3D, add_methods._add_property_object),

            #'CCONEAX' : (CCONEAX, add_methods._add_element_object),
            #'PCONEAX' : (PCONEAX, add_methods._add_property_object),
            #'AXIC' : (AXIC, add_methods._add_axic_object),
            #'AXIF' : (AXIF, add_methods._add_axif_object),
            #'CYAX' : (CYAX, add_methods._add_cyax_object),

            #'MAT3D' : (MAT3D, add_methods._add_structural_material_object),
            #'MATEV' : (MATEV, add_methods._add_structural_material_object),
            #'EQUIV' : (EQUIV, add_methods._add_structural_material_object),
            #'MATG' : (MATG, add_methods._add_structural_material_object),

            #'MATEV' : (MATEV, add_methods._add_structural_material_object),
            'NXSTRAT' : (NXSTRAT, add_methods._add_nxstrat_object),

            # hasn't been verified, links up to MAT1, MAT2, MAT9 w/ same MID
            #'CREEP' : (CREEP, add_methods._add_creep_material_object),

            #'SPCAX' : (SPCAX, add_methods._add_constraint_spc_object),
            ## parametric
            #'GMSPC' : (GMSPC, add_methods._add_constraint_spc_object),

            #'SESUP' : (SESUP, add_methods._add_sesuport_object), # pseudo-constraint

            #'CLOAD' : (CLOAD, add_methods._add_load_combination_object),
            #'LOADCYN' : (LOADCYN, add_methods._add_load_object),
            #'LOADCYH' : (LOADCYH, add_methods._add_load_object),

            'FREQ' : (FREQ, add_methods._add_freq_object),
            'FREQ1' : (FREQ1, add_methods._add_freq_object),
            'FREQ2' : (FREQ2, add_methods._add_freq_object),
            'FREQ3' : (FREQ3, add_methods._add_freq_object),
            'FREQ4' : (FREQ4, add_methods._add_freq_object),
            'FREQ5' : (FREQ5, add_methods._add_freq_object),

            'DOPTPRM' : (DOPTPRM, add_methods._add_doptprm_object),
            #'TOPVAR' : (TOPVAR, add_methods._add_topvar_object),

            #'PCONV' : (PCONV, add_methods._add_convection_property_object),
            #'PCONVM' : (PCONVM, add_methods._add_convection_property_object),

            #'VIEW' : (VIEW, add_methods._add_view_object),
            #'VIEW3D' : (VIEW3D, add_methods._add_view3d_object),

            # SOL 144
            'AEROS' : (AEROS, add_methods._add_aeros_object),
            'DIVERG' : (DIVERG, add_methods._add_diverg_object),

            # SOL 145
            'AERO' : (AERO, add_methods._add_aero_object),
            'MKAERO1' : (MKAERO1, add_methods._add_mkaero_object),
            'MKAERO2' : (MKAERO2, add_methods._add_mkaero_object),

            'NLPARM' : (NLPARM, add_methods._add_nlparm_object),
            'NLPCI' : (NLPCI, add_methods._add_nlpci_object),
            'TSTEP' : (TSTEP, add_methods._add_tstep_object),
            'TSTEP1' : (TSTEP1, add_methods._add_tstepnl_object),
            'TSTEPNL' : (TSTEPNL, add_methods._add_tstepnl_object),

            # nx_opt
            'DTABLE' : (DTABLE, add_methods._add_dtable_object),
            #'DVTREL1' : (DVTREL1, add_methods._add_dvtrel_object), # dvtrels
            #'GROUP' : (GROUP, add_methods._add_group_object), # group
            #'DMNCON' : (DMNCON, add_methods._add_dmncon_object), # dmncon

            # tables
            'TABLES1' : (TABLES1, add_methods._add_table_object),
            'TABLEST' : (TABLEST, add_methods._add_table_object),
            #'TABLEHT' : (TABLEHT, add_methods._add_table_object),
            'TABLEH1' : (TABLEH1, add_methods._add_table_object),

            # dynamic tables
            'TABLED1' : (TABLED1, add_methods._add_tabled_object),
            'TABLED2' : (TABLED2, add_methods._add_tabled_object),
            'TABLED3' : (TABLED3, add_methods._add_tabled_object),
            'TABLED4' : (TABLED4, add_methods._add_tabled_object),

            # material tables
            'TABLEM1' : (TABLEM1, add_methods._add_tablem_object),
            'TABLEM2' : (TABLEM2, add_methods._add_tablem_object),
            'TABLEM3' : (TABLEM3, add_methods._add_tablem_object),
            'TABLEM4' : (TABLEM4, add_methods._add_tablem_object),

            # other tables
            #'TABDMP1' : (TABDMP1, add_methods._add_table_sdamping_object),
            #'TABRND1' : (TABRND1, add_methods._add_random_table_object),
            #'TABRNDG' : (TABRNDG, add_methods._add_random_table_object),

            'EIGB' : (EIGB, add_methods._add_method_object),
            'EIGR' : (EIGR, add_methods._add_method_object),
            'EIGRL' : (EIGRL, add_methods._add_method_object),
            'EIGC' : (EIGC, add_methods._add_cmethod_object),
            'EIGP' : (EIGP, add_methods._add_cmethod_object),

            # 'BOUTPUT', 'BOLT', 'BOLTFOR', 'BOLTFRC',
            'BOUTPUT': (Crash, None),
            #'BOLT': (Crash, None),
            #'BOLTFOR': (Crash, None),
            'BOLTFRC': (Crash, None),
            'MBOLTUS': (Crash, None)

            #'RADCAV' : (RADCAV, add_methods._add_radcav_object), #
            #'RADLST' : (RADLST, add_methods._add_radcav_object), # TestOP2.test_bdf_op2_thermal_02
            #'RADMTX' : (RADMTX, add_methods._add_radmtx_object), # TestOP2.test_bdf_op2_thermal_02
            #'RADMT' : (Crash, None),

            #'SESUP' : (SESUP, add_methods._add_sesup_object),  # pseudo-constraint

            #'SEUSET' : (SEUSET, add_methods._add_seuset_object),
            #'SEUSET1' : (SEUSET1, add_methods._add_seuset_object),

            #'ROTORG' : (ROTORG, add_methods._add_rotor_object),
            #'ROTORD' : (ROTORD, add_methods._add_rotor_object),

            #'CYJOIN' : (CYJOIN, add_methods._add_cyjoin_object),
        }

        self._card_parser_prepare = {
            'PLOTEL': partial(self._prepare_card, self.plotel),
            'PLOTEL3': partial(self._prepare_card, self.plotel3),
            'PLOTEL4': partial(self._prepare_card, self.plotel4),
            'PLOTEL6': partial(self._prepare_card, self.plotel6),
            'PLOTEL8': partial(self._prepare_card, self.plotel8),

            'GRID': partial(self._prepare_card, self.grid),
            'SPOINT': partial(self._prepare_card, self.spoint),
            'EPOINT': partial(self._prepare_card, self.epoint),
            'CORD1R': self._prepare_cord1r,
            'CORD1C': self._prepare_cord1c,
            'CORD1S': self._prepare_cord1s,
            'CORD2R': self._prepare_cord2r,
            'CORD2C': self._prepare_cord2c,
            'CORD2S': self._prepare_cord2s,

            'MODTRAK' : partial(self._prepare_card, self.modtrak),

            # sets
            'SET1': partial(self._prepare_card, self.set1),
            #'SET2': partial(self._prepare_card, self.set2),
            'SET3': partial(self._prepare_card, self.set3),

            # spring
            'CELAS1' : partial(self._prepare_card, self.celas1),
            'CELAS2' : partial(self._prepare_card, self.celas2),
            'CELAS3' : partial(self._prepare_card, self.celas3),
            'CELAS4' : partial(self._prepare_card, self.celas4),
            'PELAS' : partial(self._prepare_card, self.pelas),
            'PELAST' : partial(self._prepare_card, self.pelast),

            # spring
            'CDAMP1' : partial(self._prepare_card, self.cdamp1),
            'CDAMP2' : partial(self._prepare_card, self.cdamp2),
            'CDAMP3' : partial(self._prepare_card, self.cdamp3),
            'CDAMP4' : partial(self._prepare_card, self.cdamp4),
            'CDAMP5' : partial(self._prepare_card, self.cdamp5),
            'PDAMP' : partial(self._prepare_card, self.pdamp),
            'PDAMP5': partial(self._prepare_card, self.pdamp5),
            'PDAMPT': partial(self._prepare_card, self.pdampt),

            # visc
            'CVISC' : partial(self._prepare_card, self.cvisc),
            'PVISC' : partial(self._prepare_card, self.pvisc),

            # bush
            'CBUSH' : partial(self._prepare_card, self.cbush),
            'PBUSH' : partial(self._prepare_card, self.pbush),
            'PBUSHT' : partial(self._prepare_card, self.pbusht),

            'CBUSH1D' : partial(self._prepare_card, self.cbush1d),
            'PBUSH1D' : partial(self._prepare_card, self.pbush1d),
            #'CBUSH2D' : partial(self._prepare_card, self.cbush2d),
            #'PBUSH2D' : partial(self._prepare_card, self.pbush2d),

            # gap
            'CGAP' : partial(self._prepare_card, self.cgap),
            'PGAP' : partial(self._prepare_card, self.pgap),

            # fast
            'CFAST' : partial(self._prepare_card, self.cfast),
            'PFAST' : partial(self._prepare_card, self.pfast),

            'GENEL' : partial(self._prepare_card, self.genel),

            # mass
            'PMASS' : partial(self._prepare_card, self.pmass),
            'CMASS1' : partial(self._prepare_card, self.cmass1),
            'CMASS2' : partial(self._prepare_card, self.cmass2),
            'CMASS3' : partial(self._prepare_card, self.cmass3),
            'CMASS4' : partial(self._prepare_card, self.cmass4),
            'CONM1' : partial(self._prepare_card, self.conm1),
            'CONM2' : partial(self._prepare_card, self.conm2),

            # nonstructural mass
            'NSMADD' : partial(self._prepare_card, self.nsmadd),
            'NSM' : partial(self._prepare_card, self.nsm),
            'NSM1' : partial(self._prepare_card, self.nsm1),
            # lumped
            'NSML' : partial(self._prepare_card, self.nsml),
            'NSML1' : partial(self._prepare_card, self.nsml1),

            # rod
            'CONROD' : partial(self._prepare_card, self.conrod),
            'CROD' : partial(self._prepare_card, self.crod),
            'PROD' : partial(self._prepare_card, self.prod),

            # tube
            'CTUBE' : partial(self._prepare_card, self.ctube),
            'PTUBE' : partial(self._prepare_card, self.ptube),

            # bar
            'CBAR': partial(self._prepare_card, self.cbar),
            'PBAR': partial(self._prepare_card, self.pbar),
            'PBARL': partial(self._prepare_card, self.pbarl),
            'PBRSECT': partial(self._prepare_card, self.pbrsect),
            'CBARAO' : partial(self._prepare_card, self.cbarao),

            # beam
            'CBEAM': partial(self._prepare_card, self.cbeam),
            'PBEAM': partial(self._prepare_card, self.pbeam),
            'PBEAML': partial(self._prepare_card, self.pbeaml),
            'PBCOMP' : partial(self._prepare_card, self.pbcomp),
            #'PBMSECT' : (PBMSECT, add_methods._add_property_object),
            'POINT' : partial(self._prepare_card, self.point),

            # bend
            'CBEND' : partial(self._prepare_card, self.cbend),
            'PBEND' : partial(self._prepare_card, self.pbend),

            # shear
            'CSHEAR': partial(self._prepare_card, self.cshear),
            'PSHEAR': partial(self._prepare_card, self.pshear),

            # shells
            #'SNORM': partial(self._prepare_card, self.snorm),
            'CQUAD4': partial(self._prepare_card, self.cquad4),
            'CTRIA3': partial(self._prepare_card, self.ctria3),
            'CQUAD8': partial(self._prepare_card, self.cquad8),
            'CQUADR': partial(self._prepare_card, self.cquadr),
            'CTRIA6': partial(self._prepare_card, self.ctria6),
            'CTRIAR': partial(self._prepare_card, self.ctriar),
            'CQUAD': partial(self._prepare_card, self.cquad),
            'PSHELL': partial(self._prepare_card, self.pshell),
            'PCOMP': partial(self._prepare_card, self.pcomp),
            'PCOMPG': partial(self._prepare_card, self.pcompg),
            'PSHLN1' : partial(self._prepare_card, self.pshln1),
            'PSHLN2' : partial(self._prepare_card, self.pshln2),

            # planar shells
            'PLPLANE': partial(self._prepare_card, self.plplane),

            # plate stress
            'PPLANE' : partial(self._prepare_card, self.pplane),
            'CPLSTS3' : partial(self._prepare_card, self.cplsts3),
            'CPLSTS4' : partial(self._prepare_card, self.cplsts4),
            'CPLSTS6' : partial(self._prepare_card, self.cplsts6),
            'CPLSTS8' : partial(self._prepare_card, self.cplsts8),

            # plate strain
            # what are the properties?
            'CPLSTN3' : partial(self._prepare_card, self.cplstn3),
            'CPLSTN4' : partial(self._prepare_card, self.cplstn4),
            'CPLSTN6' : partial(self._prepare_card, self.cplstn6),
            'CPLSTN8' : partial(self._prepare_card, self.cplstn8),

            # axisymmetric shells
            'CQUADX': partial(self._prepare_card, self.cquadx),
            'CQUADX4' : partial(self._prepare_card, self.cquadx4),
            'CQUADX8' : partial(self._prepare_card, self.cquadx8),
            'CTRIAX': partial(self._prepare_card, self.ctriax),
            'CTRIAX6' : partial(self._prepare_card, self.ctriax6),
            'CTRAX3' : partial(self._prepare_card, self.ctrax3),
            'CTRAX6' : partial(self._prepare_card, self.ctrax6),

            # solids
            'CTETRA' : partial(self._prepare_card, self.ctetra),
            'CPENTA' : partial(self._prepare_card, self.cpenta),
            'CHEXA' : partial(self._prepare_card, self.chexa),
            'CPYRAM' : partial(self._prepare_card, self.cpyram),
            'CHEXCZ' : partial(self._prepare_card, self.chexcz),
            'PSOLID': partial(self._prepare_card, self.psolid),
            'PLSOLID': partial(self._prepare_card, self.plsolid),
            'PCOMPS': partial(self._prepare_card, self.pcomps),
            'PCOMPLS': partial(self._prepare_card, self.pcompls),

            # acoustic shells?
            'CAABSF': partial(self._prepare_card, self.caabsf),
            'PAABSF': partial(self._prepare_card, self.paabsf),

            #'CHACAB' : partial(self._prepare_card, self.chacab),
            #'CHACBR' : partial(self._prepare_card, self.chacbr),
            #'PACABS': (PACABS, add_methods._add_acoustic_property_object),
            #'PACBAR': (PACBAR, add_methods._add_acoustic_property_object),

            # thermal elements
            'CHBDYE': partial(self._prepare_card, self.chbdye),
            'CHBDYG': partial(self._prepare_card, self.chbdyg),
            'CHBDYP': partial(self._prepare_card, self.chbdyp),
            'PHBDY': partial(self._prepare_card, self.phbdy),
            #'TEMPRB' : (TEMPRB, add_methods._add_thermal_load_object),
            #'TEMPB3' : (TEMPB3, add_methods._add_thermal_load_object),

            'CONV': partial(self._prepare_card, self.conv),
            'PCONV': partial(self._prepare_card, self.pconv),
            'CONVM': partial(self._prepare_card, self.convm),
            'PCONVM': partial(self._prepare_card, self.pconvm),

            # static loads
            'LOAD' : partial(self._prepare_card, self.load),
            'PLOAD' : partial(self._prepare_card, self.pload),
            'PLOAD1' : partial(self._prepare_card, self.pload1),
            'PLOAD2' : partial(self._prepare_card, self.pload2),
            'PLOAD4' : partial(self._prepare_card, self.pload4),
            'FORCE' : partial(self._prepare_card, self.force),
            'FORCE1' : partial(self._prepare_card, self.force1),
            'FORCE2' : partial(self._prepare_card, self.force2),
            'MOMENT' : partial(self._prepare_card, self.moment),
            'MOMENT1' : partial(self._prepare_card, self.moment1),
            'MOMENT2' : partial(self._prepare_card, self.moment2),
            'GRAV' : partial(self._prepare_card, self.grav),
            'ACCEL' : partial(self._prepare_card, self.accel),
            'ACCEL1' : partial(self._prepare_card, self.accel1),
            'SLOAD': partial(self._prepare_card, self.sload),
            'TEMP': partial(self._prepare_card, self.temp),
            'TEMPD': partial(self._prepare_card, self.tempd),
            #'DTEMP': partial(self._prepare_card, self.dtemp),
            'LSEQ': partial(self._prepare_card, self.lseq),
            'SPCD' : partial(self._prepare_card, self.spcd), # enforced displacement; load
            'DEFORM' : partial(self._prepare_card, self.deform),  # enforced displacement
            'RFORCE' : partial(self._prepare_card, self.rforce),
            'RFORCE1' : partial(self._prepare_card, self.rforce1),

            # axisymmetric loads
            #'PLOADX1' : partial(self._prepare_card, self.ploadx1),
            #'FORCEAX' : (FORCEAX, add_methods._add_load_object),
            #'PRESAX' : (PRESAX, add_methods._add_load_object),  # axisymmetric

            # thermal loads
            'QBDY1': partial(self._prepare_card, self.qbdy1),
            'QBDY2': partial(self._prepare_card, self.qbdy2),
            'QBDY3': partial(self._prepare_card, self.qbdy3),
            'QVOL': partial(self._prepare_card, self.qvol),
            'TEMPBC': partial(self._prepare_card, self.tempbc),
            'RADBC': partial(self._prepare_card, self.radbc),
            'RADM': partial(self._prepare_card, self.radm),

            # radiation
            'RADSET' : partial(self._prepare_card, self.radset),

            # dynamic loads
            'DLOAD' : partial(self._prepare_card, self.dload),
            'DAREA' : partial(self._prepare_card, self.darea),
            'TLOAD1' : partial(self._prepare_card, self.tload1),
            'TLOAD2' : partial(self._prepare_card, self.tload2),
            'RLOAD1' : partial(self._prepare_card, self.rload1),
            'RLOAD2' : partial(self._prepare_card, self.rload2),
            'TIC' : partial(self._prepare_card, self.tic),
            'TF' : partial(self._prepare_card, self.tf),
            'QVECT': partial(self._prepare_card, self.qvect),
            'QHBDY': partial(self._prepare_card, self.qhbdy),

            # offset
            'DPHASE': partial(self._prepare_card, self.dphase),
            'DELAY': partial(self._prepare_card, self.delay),

            # random laods
            'RANDPS': partial(self._prepare_card, self.randps),
            'ACSRCE': partial(self._prepare_card, self.acsrce),
            #'RANDT1' : (RANDT1, add_methods._add_dload_entry), # random

            # materials
            'MAT1' : partial(self._prepare_card, self.mat1),
            'MAT2' : partial(self._prepare_card, self.mat2),
            'MAT3' : partial(self._prepare_card, self.mat3),
            ## there is no MAT6 or MAT7
            'MAT8' : partial(self._prepare_card, self.mat8),
            'MAT9' : partial(self._prepare_card, self.mat9),
            'MAT10' : partial(self._prepare_card, self.mat10),
            'MAT11' : partial(self._prepare_card, self.mat11),
            'MAT10C' : partial(self._prepare_card, self.mat10c),
            'MATORT' : partial(self._prepare_card, self.matort),

            'MATS1' : partial(self._prepare_card, self.mats1),
            #'MATS3' : partial(self._prepare_card, self.mats3),
            #'MATS8' : partial(self._prepare_card, self.mats8),

            'MATT1' : partial(self._prepare_card, self.matt1),
            #'MATT2' : (MATT2, add_methods._add_material_dependence_object),
            #'MATT3' : (MATT3, add_methods._add_material_dependence_object),
            #'MATT4' : (MATT4, add_methods._add_material_dependence_object),
            #'MATT5' : (MATT5, add_methods._add_material_dependence_object),
            'MATT8' : partial(self._prepare_card, self.matt8),
            'MATT9' : partial(self._prepare_card, self.matt9),

            'MATHE' : partial(self._prepare_card, self.mathe),  # MOONEY only; no OGDEN, ABOYCE, ...
            'MATHP' : partial(self._prepare_card, self.mathp),

            # thermal materials
            'MAT4' : partial(self._prepare_card, self.mat4),
            'MAT5' : partial(self._prepare_card, self.mat5),

            # rigid elements
            'RBAR': partial(self._prepare_card, self.rbar),
            'RBAR1': partial(self._prepare_card, self.rbar1),
            'RBE1': partial(self._prepare_card, self.rbe1),
            'RBE2': partial(self._prepare_card, self.rbe2),
            'RBE3': partial(self._prepare_card, self.rbe3),
            'RROD': partial(self._prepare_card, self.rrod),
            #'RSSCON': partial(self._prepare_card, self.rsscon),  # not supported
            #'RSPLINE' : (RSPLINE, add_methods._add_rigid_element_object),

            # mpcs
            'MPC': partial(self._prepare_card, self.mpc),
            'MPCADD' : partial(self._prepare_card, self.mpcadd),

            'SPC': partial(self._prepare_card, self.spc),
            'SPC1': partial(self._prepare_card, self.spc1),
            'SPCADD' : partial(self._prepare_card, self.spcadd),
            'SPCOFF' : partial(self._prepare_card_by_method, self.spcoff.add_set_card),
            'SPCOFF1' : partial(self._prepare_card_by_method, self.spcoff.add_set1_card),
            'BNDFIX' : partial(self._prepare_card_by_method, self.bndfix.add_set_card),
            'BNDFIX1' : partial(self._prepare_card_by_method, self.bndfix.add_set1_card),
            'BNDFREE' : partial(self._prepare_card_by_method, self.bndfree.add_set_card),
            'BNDFREE1' : partial(self._prepare_card_by_method, self.bndfree.add_set1_card),

            # nx-contact
            'BSURF' : partial(self._prepare_card_by_method, self.bsurf.add_card),
            'BSURFS' : partial(self._prepare_card_by_method, self.bsurfs.add_card),
            'BCPROP' : partial(self._prepare_card_by_method, self.bcprop.add_card),
            'BCPROPS' : partial(self._prepare_card_by_method, self.bcprops.add_card),

            # nx glue contact
            'BGSET' : partial(self._prepare_card_by_method, self.bgset.add_card),
            'BGADD' : partial(self._prepare_card_by_method, self.bgadd.add_card),
            'BEDGE' : partial(self._prepare_card, self.bedge),

            # nx general contact
            'BCTSET' : partial(self._prepare_card_by_method, self.bctset.add_card),
            'BCTADD' : partial(self._prepare_card_by_method, self.bctadd.add_card),

            # ??? contact
            #'BCPARA' : (BCPARA, add_methods._add_bcpara_object),
            #'BCBODY' : (BCBODY, add_methods._add_bcbody_object),
            'BCONP' : partial(self._prepare_card, self.bconp),
            'BFRIC' : partial(self._prepare_card, self.bfric),
            'BLSEG' : partial(self._prepare_card, self.blseg),
            'BCRPARA' : partial(self._prepare_card, self.bcrpara),
            'BCTPARA' : (BCTPARA, add_methods._add_bctpara_object),

            # pseudo-constraint
            'SUPORT': partial(self._prepare_card_by_method, self.suport.add_set_card),
            'SUPORT1': partial(self._prepare_card_by_method, self.suport.add_set1_card),

            'ASET': partial(self._prepare_card_by_method, self.aset.add_set_card),
            'ASET1': partial(self._prepare_card_by_method, self.aset.add_set1_card),

            'BSET': partial(self._prepare_card_by_method, self.bset.add_set_card),
            'BSET1': partial(self._prepare_card_by_method, self.bset.add_set1_card),

            'CSET': partial(self._prepare_card_by_method, self.cset.add_set_card),
            'CSET1': partial(self._prepare_card_by_method, self.cset.add_set1_card),

            'QSET': partial(self._prepare_card_by_method, self.qset.add_set_card),
            'QSET1': partial(self._prepare_card_by_method, self.qset.add_set1_card),

            'OMIT': partial(self._prepare_card_by_method, self.omit.add_set_card),
            'OMIT1': partial(self._prepare_card_by_method, self.omit.add_set1_card),

            'USET': partial(self._prepare_card_by_method, self.uset.add_set_card),
            'USET1': partial(self._prepare_card_by_method, self.uset.add_set1_card),

            'SESET': partial(self._prepare_card_by_method, self.seset.add_card),
            'SEBSET': partial(self._prepare_card_by_method, self.sebset.add_set_card),
            'SECSET': partial(self._prepare_card_by_method, self.secset.add_set_card),
            'SEQSET': partial(self._prepare_card_by_method, self.seqset.add_set_card),

            'SEBSET1': partial(self._prepare_card_by_method, self.sebset.add_set1_card),
            'SECSET1': partial(self._prepare_card_by_method, self.secset.add_set1_card),
            'SEQSET1': partial(self._prepare_card_by_method, self.seqset.add_set1_card),

            'RELEASE': partial(self._prepare_card_by_method, self.release.add_card),

            # aero
            'CAERO1': partial(self._prepare_card, self.caero1),
            'CAERO2': partial(self._prepare_card, self.caero2),
            'CAERO3': partial(self._prepare_card, self.caero3),
            'CAERO4': partial(self._prepare_card, self.caero4),
            'CAERO5': partial(self._prepare_card, self.caero5),
            'CAERO7': partial(self._prepare_card, self.caero7), # zona

            'AESURF': partial(self._prepare_card, self.aesurf),
            'AESURFS' : partial(self._prepare_card, self.aesurfs),
            'AEPARM': partial(self._prepare_card, self.aeparm),
            'AESTAT': partial(self._prepare_card, self.aestat),

            'PAERO1': partial(self._prepare_card, self.paero1),
            'PAERO2': partial(self._prepare_card, self.paero2),
            'PAERO3': partial(self._prepare_card, self.paero3),
            'PAERO4': partial(self._prepare_card, self.paero4),
            'PAERO5': partial(self._prepare_card, self.paero5),

            'SPLINE1': partial(self._prepare_card, self.spline1),
            'SPLINE2': partial(self._prepare_card, self.spline2),
            'SPLINE3': partial(self._prepare_card, self.spline3),
            'SPLINE4': partial(self._prepare_card, self.spline4),
            'SPLINE5': partial(self._prepare_card, self.spline5),

            'AELIST': partial(self._prepare_card, self.aelist),    # aesurf boxes
            'AEFACT' : partial(self._prepare_card, self.aefact),   # caero paneling
            'AELINK': partial(self._prepare_card, self.aelink),    # control surface linkage
            'AECOMP': partial(self._prepare_card, self.aecomp),    # monpntx set1/aelist ids
            'AECOMPL': partial(self._prepare_card, self.aecompl),  # groups AECOMP/AECOMPLs

            # aero loads------------------------
            # sol 144 - trim/divergence
            'CSSCHD': partial(self._prepare_card, self.csschd),
            'TRIM': partial(self._prepare_card, self.trim),
            #'TRIM2': partial(self._prepare_card, self.trim2),
            # diverg

            # sol 145 - flutter
            'FLFACT' : partial(self._prepare_card, self.flfact),
            'FLUTTER' : partial(self._prepare_card, self.flutter),

            # sol 146 - gust
            'GUST' : partial(self._prepare_card, self.gust),

            'MONPNT1' : partial(self._prepare_card, self.monpnt1),
            'MONPNT2' : partial(self._prepare_card, self.monpnt2),
            'MONPNT3' : partial(self._prepare_card, self.monpnt3),
            #'MONDSP1' : (MONDSP1, add_methods._add_monpnt_object),

            # optimization
            'DESVAR': partial(self._prepare_card, self.desvar),
            'DLINK': partial(self._prepare_card, self.dlink),
            'DVGRID': partial(self._prepare_card, self.dvgrid),
            'DRESP1': partial(self._prepare_card, self.dresp1),
            'DRESP2': partial(self._prepare_card, self.dresp2),
            'DCONSTR': partial(self._prepare_card, self.dconstr),
            'DVPREL1': partial(self._prepare_card, self.dvprel1),
            'DVPREL2': partial(self._prepare_card, self.dvprel2),
            'DVCREL1': partial(self._prepare_card, self.dvcrel1),
            'DVCREL2': partial(self._prepare_card, self.dvcrel2),
            'DVMREL1': partial(self._prepare_card, self.dvmrel1),
            'DVMREL2': partial(self._prepare_card, self.dvmrel2),

            'DCONADD': partial(self._prepare_card, self.dconadd),
            'DDVAL': partial(self._prepare_card, self.ddval),
            'DSCREEN': partial(self._prepare_card, self.dscreen),

            #'DTABLE' : (DTABLE, add_methods._add_dtable_object),
            #'DRESP3' : (DRESP3, add_methods._add_dresp_object),

            #'CORD3G' : self._prepare_CORD3G,

            'DTI' : self._prepare_dti,
            'DMIG' : self._prepare_dmig,
            'DMIAX' : self._prepare_dmiax,
            'DMI' : self._prepare_dmi,
            'DMIJ' : self._prepare_dmij,
            'DMIK' : self._prepare_dmik,
            'DMIJI' : self._prepare_dmiji,
            #'RINGFL' : self._prepare_ringfl,

            'DEQATN' : self._prepare_dequatn,

            #'TEMPAX' : self._prepare_tempax,
            #'CONVM' : self._prepare_convm,
            #'CONV' : self._prepare_conv,
            # GRDSET-will be last card to update from _card_parser_prepare
            'GRDSET' : self._prepare_grdset,
            'BAROR' : self._prepare_baror,
            'BEAMOR' : self._prepare_beamor,
            #'BDYOR': self._prepare_bdyor,

            #'ACMODL' : self._prepare_acmodl,
        }
        ncards_supported = len(self._card_parser_prepare) + len(self._card_parser)
        #print(f'ncards = {ncards_supported}')
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

    def _prepare_cord1r(self, unused_card: list[str], card_obj: BDFCard, comment='') -> int:
        """adds a CORD1R"""
        out = self.coord.add_cord1r_bdf(card_obj, icard=0, comment=comment)
        if card_obj.field(5):
            out = self.coord.add_cord1r_bdf(card_obj, icard=1, comment=comment)
        return out
    def _prepare_cord1c(self, unused_card: list[str], card_obj: BDFCard, comment='') -> int:
        """adds a CORD1C"""
        out = self.coord.add_cord1c_bdf(card_obj, icard=0, comment=comment)
        if card_obj.field(5):
            out = self.coord.add_cord1c_bdf(card_obj, icard=1, comment=comment)
        return out
    def _prepare_cord1s(self, unused_card: list[str], card_obj: BDFCard, comment='') -> int:
        """adds a CORD1S"""
        out = self.coord.add_cord1s_bdf(card_obj, icard=0, comment=comment)
        if card_obj.field(5):
            out = self.coord.add_cord1s_bdf(card_obj, icard=1, comment=comment)
        return out

    def _prepare_cord2r(self, unused_card: list[str], card_obj: BDFCard, comment='') -> int:
        """adds a CORD2R"""
        out = self.coord.add_cord2r_bdf(card_obj, comment=comment)
        return out
    def _prepare_cord2c(self, unused_card: list[str], card_obj: BDFCard, comment='') -> int:
        """adds a CORD2C"""
        out = self.coord.add_cord2c_bdf(card_obj, comment=comment)
        return out
    def _prepare_cord2s(self, unused_card: list[str], card_obj: BDFCard, comment='') -> int:
        """adds a CORD2S"""
        out = self.coord.add_cord2s_bdf(card_obj, comment=comment)
        return out

    def _prepare_card(self, cls, unused_card: list[str], card_obj: BDFCard, comment='') -> None:
        """adds a CBAR"""
        i = cls.add_card(card_obj, comment=comment)
        return i

    def _prepare_card_by_method(self, func: Callable, unused_card: list[str], card_obj: BDFCard, comment='') -> None:
        """adds a CBAR"""
        i = func(card_obj, comment=comment)
        return i

    def _prepare_grdset(self, unused_card: list[str], card_obj: BDFCard, comment='') -> None:
        """adds a GRDSET"""
        assert self.grdset is None, self.grdset
        self.grdset = GRDSET.add_card(card_obj, comment=comment)
        return self.grdset
    def _prepare_baror(self, unused_card: list[str], card_obj: BDFCard, comment='') -> None:
        """adds a BAROR"""
        assert self.baror is None, self.baror
        self.baror = BAROR.add_card(self, card_obj, comment=comment)
        return self.baror
    def _prepare_beamor(self, unused_card: list[str], card_obj: BDFCard, comment='') -> None:
        """adds a BEAMOR"""
        assert self.beamor is None, self.beamor
        self.beamor = BEAMOR.add_card(card_obj, comment=comment)
        return self.beamor
    def _prepare_bdyor(self, unused_card: list[str], card_obj: BDFCard, comment='') -> None:
        """adds a BYDOR"""
        assert self.bdyor is None, self.bdyor
        self.bdyor = BDYOR.add_card(card_obj, comment=comment)
        return self.bdyor

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

    def _prepare_dti(self, unused_card_name, card_obj, comment='') -> DTI:
        """adds a DTI"""
        name = string(card_obj, 1, 'name')
        if name == 'UNITS':
            dti = DTI_UNITS.add_card(card_obj, comment=comment)
        else:
            dti = DTI.add_card(card_obj, comment=comment)
        self._add_methods._add_dti_object(dti)
        return dti

    def _prepare_dmig(self, unused_card: list[str], card_obj: BDFCard, comment='') -> DMIG:
        """adds a DMIG"""
        name = string(card_obj, 1, 'name')
        field2 = integer_or_string(card_obj, 2, 'flag')

        if name == 'UACCEL':  # special DMIG card
            if field2 == 0:
                dmig = DMIG_UACCEL.add_card(card_obj, comment=comment)
                self._add_methods._add_dmig_object(dmig)
            else:
                dmig = -1
                self._dmig_temp[name].append((card_obj, comment))
        else:
            field2 = integer_or_string(card_obj, 2, 'flag')
            if field2 == 0:
                dmig = DMIG.add_card(card_obj, comment=comment)
                self._add_methods._add_dmig_object(dmig)
            else:
                dmig = -1
                self._dmig_temp[name].append((card_obj, comment))
        return dmig

    def _prepare_dmix(self, class_obj, add_method, card_obj, comment='') -> Union[DMI, DMIJ, DMIJI, DMIK]:
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

    def _prepare_dmiax(self, unused_card: list[str], card_obj: BDFCard, comment='') -> DMIAX:
        """adds a DMIAX"""
        return self._prepare_dmix(DMIAX, self._add_methods._add_dmiax_object, card_obj, comment=comment)

    def _prepare_dmi(self, unused_card: list[str], card_obj: BDFCard, comment='') -> DMI:
        """adds a DMI"""
        return self._prepare_dmix(DMI, self._add_methods._add_dmi_object, card_obj, comment=comment)

    def _prepare_dmij(self, unused_card: list[str], card_obj: BDFCard, comment='') -> DMIJ:
        """adds a DMIJ"""
        return self._prepare_dmix(DMIJ, self._add_methods._add_dmij_object, card_obj, comment=comment)

    def _prepare_dmik(self, unused_card: list[str], card_obj: BDFCard, comment='') -> DMIK:
        """adds a DMIK"""
        return self._prepare_dmix(DMIK, self._add_methods._add_dmik_object, card_obj, comment=comment)

    def _prepare_dmiji(self, unused_card: list[str], card_obj: BDFCard, comment='') -> DMIJI:
        """adds a DMIJI"""
        return self._prepare_dmix(DMIJI, self._add_methods._add_dmiji_object, card_obj, comment=comment)

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
                 comment: str='', ifile=None,
                 is_list: bool=True, has_none: bool=True) -> Any:
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
            can there be trailing Nones in the card data (e.g. ['GRID', 1, 2, 3.0, 4.0, 5.0, '])

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
        self._add_card_helper(card_obj, card_name, card_name, ifile,
                              comment=comment)
        return card_obj

    def add_card_lax(self, card_lines: list[str], card_name: str,
                     comment: str='', ifile=None, is_list: bool=True, has_none: bool=True) -> Any:
        """see ``add_card``"""
        card_name = card_name.upper()
        #if card_name not in self.card_count:
            #print(card_name)
        card_obj, unused_card = self.create_card_object(
            card_lines, card_name,
            is_list=is_list, has_none=has_none)
        self._add_card_helper_lax(card_obj, card_name, card_name, ifile,
                                  comment=comment)
        return card_obj

    def add_card_fields(self, card_lines: list[str], card_name: str,
                        comment: str='', has_none: bool=True):
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

    def add_card_lines(self, card_lines: list[str], card_name: str, comment: str='', has_none=True):
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

    def _get_npoints_nids_allnids(self) -> np.ndarray:
        """helper method for get_xyz_in_coord"""
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

    def get_xyz_in_coord(self, cid: int=0, fdtype: str='float64', sort_ids: bool=True):
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

        if card_name in self._crash_cards:
            raise DisabledCardError(f'stopping on card={card_name}')

        if card_name in self._card_parser:
            #self.log.debug(card_name)
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
            #self.log.debug(card_name)
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

    def _add_card_helper_lax(self, card_obj: BDFCard, card: list[str],
                             card_name: str, ifile: int,
                             comment: str='') -> None:
        #if card_name not in ['GRID', 'CQUAD4', 'CTRIA3']:
            #print(card_obj)

        if card_name == 'ECHOON':
            self.echo = True
            return
        elif card_name == 'ECHOOFF':
            self.echo = False
            return

        if self.echo and not self.force_echo_off:
            _echo_card(card, card_obj)

        if card_name in self._crash_cards:
            raise DisabledCardError(f'stopping on card={card_name}')

        if card_name in self._card_parser:

            card_class, add_card_function = self._card_parser[card_name]
            if hasattr(card_class, 'add_card_lax'):
                class_instance = card_class.add_card_lax(card_obj, comment=comment)
            else:

                class_instance = card_class.add_card(card_obj, comment=comment)
            add_card_function(class_instance)

        elif card_name in self._card_parser_prepare:
            add_card_function = self._card_parser_prepare[card_name]
            add_card_function(card, card_obj, comment=comment)
        else:
            #raise RuntimeError(card_obj)
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
                if self.is_lax_parser and hasattr(card_class, 'add_card_lax'):
                    class_instance = card_class.add_card_lax(card_obj, comment=comment)
                else:
                    class_instance = card_class.add_card(card_obj, comment=comment)
                card_idi = add_card_function(class_instance)
                if card_name not in OBJ_CARDS:
                    if not isinstance(card_idi, int):
                        msg = (f'card_name={card_name!r} card_idi={card_idi}; '
                               '\nvalue should be an int or defined in OBJ_CARDS')
                        raise TypeError(msg)

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
                raise
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
            self.log.warning(f'rejecting {card_obj}')
            #raise RuntimeError(card_obj)
            self.reject_cards.append(card_obj)

    def get_bdf_stats(self, return_type: str='string') -> Union[str, list[str]]:
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
        msg = ''
        if hasattr(self, 'sol'):
            msg += f'sol = {self.sol}\n'
        for card, count in self.card_count.items():
            msg += f'{card}: {count}\n'
        return msg
        return get_bdf_stats(self, return_type=return_type)

    def get_displacement_index_xyz_cp_cd(self, fdtype: str='float64',
                                         idtype: str='int32',
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
        nnodes = self.grid.n # + self.spoint.n
        spoints = None
        epoints = None
        nrings = len(self.ringaxs)
        nspoints = self.spoint.n
        nepoints = 0
        assert nnodes + nspoints > 0, (nnodes, nspoints)

        #nepoints = self.epoint.n
        if nspoints:
            spoints = list(self.spoint.ids)
        if nepoints:
            epoints = list(self.epoint.ids)

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
                               fdtype: str='float64',
                               idtype: str='int32') -> tuple[np.ndarray, np.ndarray, np.ndarray,
                                                             dict[int, np.ndarray], dict[int, np.ndarray]]:
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
        if cid == 0:
            #xyz_cid = self.grid.xyz_cid0()
            xyz_cid = xyz_cp
        else:
            xyz_cid = self.transform_xyzcp_to_xyz_cid(xyz_cp, nids, icp_transform,
                                                      cid=cid, in_place=False, atol=1e-6)
        return nid_cp_cd, xyz_cid, xyz_cp, icd_transform, icp_transform

    def transform_xyzcp_to_xyz_cid(self, xyz_cp: np.ndarray,
                                   nids: np.ndarray,
                                   icp_transform: dict[int, np.ndarray],
                                   cid: int=0,
                                   in_place: bool=False,
                                   atol: float=1e-6) -> np.ndarray:
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
        asdf
        #F:\work\pyNastran\examples\femap_examples\Support\nast\tpl\heli112em7.dat
        # this is used when xref=False (only for vectorized=True)
        # we now require nids, where the other approach
        # (the one with xref=True) does not
        in_place = False

        cps_to_check = list(self.coord.coord_id)
        cps_to_check.sort()
        assert 0 in cps_to_check, cps_to_check

        nodes = self.grid
        coord2 = self.coord.slice_card_by_id(cid)
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
                    i1, i2, i3 = np.searchsorted(nids, coord.node_ids)
                    cp1 = nodes.cp[i1]
                    cp2 = nodes.cp[i2]
                    cp3 = nodes.cp[i3]
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

    def _old_card_fields(self, card_lines: list[str], card_name: str,
                         log: SimpleLogger,
                         is_list: bool=False, has_none: bool=True,
                         is_dynamic_syntax: bool=False) -> BDFCard:
        """replication helper"""
        if is_list:
            fields = card_lines
        else:
            fields = to_fields(card_lines, card_name)

        # apply OPENMDAO syntax
        if is_dynamic_syntax:
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

    def _expand_replication(self, card_name: str, icard: int,
                            cards_list, card_lines_new, dig: bool=True):
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
                     card_count: dict[str, int],
                     strict: bool=True) -> None:
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
            self._parse_cards_list(cards_list, strict=strict)

    def _parse_cards_dict(self, cards_dict: dict[str, list[str]]) -> None:
        """parses the cards that are in dictionary format"""
        if self.save_file_structure:
            raise NotImplementedError('save_file_structure=True is not supported\n%s' % (
                list(cards_dict.keys())))

        for card_name, cards in sorted(cards_dict.items()):
            #if card_name == 'GRID':
                #print(card_name, len(cards))
            try:
                is_reject = self.is_reject(card_name)
            except ReplicationError as error:
                card_strs = [f'{icard}: {str(card)}' for icard, (comment, card, ifile) in enumerate(cards)]
                msg = '\n'.join(card_strs)
                raise ReplicationError('Unparsable Replication:\n' + msg) from error

            if is_reject:
                self.log.info(f'    rejecting card_name = {card_name}')
                for comment, card_lines, unused_ifile_iline in cards:
                    self.increase_card_count(card_name)
                    if card_name in self._crash_cards:
                        raise DisabledCardError(f'stopping on card={card_name}')
                    self.reject_lines.append([_format_comment(comment)] + card_lines)
            else:
                for comment, card_lines, (ifile, unused_iline) in cards:
                    card_idi = self.add_card(card_lines, card_name, comment=comment, ifile=ifile,
                                             is_list=False, has_none=False)

    def _parse_cards_list(self, cards_list: list[str], strict: bool=True):
        """parses the cards that are in list format"""
        add_card = self.add_card if strict else self.add_card_lax
        del strict

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
                        add_card(replicated_card, replicated_card[0], comment=comment,
                                 is_list=True, has_none=True)
                    continue

                if self.is_reject(card_name):
                    try:
                        self.reject_card_lines(card_name, card_lines, comment=comment)
                    except:
                        card_old = cards_list[icard-1]
                        old_card_name, old_comment, old_card_lines, unused_ifile_iline = card_old
                        self.log.error('Last card was:\n%s' % '\n'.join(old_card_lines))
                        print('Last card was:\n%s' % '\n'.join(old_card_lines))
                        raise
                else:
                    add_card(card_lines, card_name, comment=comment, ifile=ifile,
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
                    i = 0
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

    def _check_pynastran_header(self, lines: list[str], check_header: bool=True) -> None:
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
            exec(self.code_block)

#---------------------------------------------------------------------------------------------------
    # HDF5
    def _read_bdf_cards(self, bdf_filename: Optional[str]=None,
                        punch: bool=False,
                        read_includes: bool=True, encoding: Optional[str]=None) -> dict[str, list[Any]]:
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
        out = obj.get_lines(bdf_filename, punch=self.punch, make_ilines=True)
        (system_lines, executive_control_lines, case_control_lines,
         bulk_data_lines, bulk_data_ilines,
         superelement_lines, superelement_ilines) = out
        self._set_pybdf_attributes(obj, save_file_structure=False)

        self.system_command_lines = system_lines
        self.executive_control_lines = executive_control_lines
        self.case_control_lines = case_control_lines

        sol, method, sol_iline, app = parse_executive_control_deck(executive_control_lines)
        self.app = app
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

    def create_subcases(self, subcase_ids: Union[int, list[int], None]=None) -> dict[int, Subcase]:
        """creates a series of subcases"""
        if subcase_ids is None:
            subcase_ids = []
        elif isinstance(subcase_ids, int):
            subcase_ids = [subcase_ids]

        self.punch = False
        if self.case_control_deck is None:
            self.case_control_deck = CaseControlDeck([], log=self.log)

        subcases = {}
        for subcase_id in subcase_ids:
            subcase = self.case_control_deck.create_new_subcase(subcase_id)
            subcases[subcase_id] = subcase
        return subcases

    def _parse_cards_hdf5(self, cards: dict[str, Any], unused_card_count) -> dict[str, list[Any]]:
        """creates card objects and adds the parsed cards to the deck"""
        self.echo = False
        cards_out = {}
        for card_name, card in sorted(cards.items()):
            cards_list = []
            cards_out[card_name] = cards_list
            if self.is_reject(card_name):
                self.log.info(f'    rejecting card_name = {card_name}')
                for comment, card_lines, unused_ifile_iline in card:
                    self.increase_card_count(card_name)
                    self.reject_lines.append([_format_comment(comment)] + card_lines)
            else:
                for comment, card_lines, unused_ifile_iline in card:
                    class_instance = self._add_card_hdf5(card_lines, card_name, comment=comment,
                                                         is_list=False, has_none=False)
                    cards_list.append(class_instance)
        return cards_out

    def _add_card_hdf5(self, card_lines: list[str], card_name: str,
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

    def _add_card_helper_hdf5(self, card_obj: BDFCard,
                              card: list[str],
                              card_name: str, comment: str='') -> None:
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
            except Exception:
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

    def __deepcopy__(self, memo: dict[str, Any]):
        """performs a deepcopy"""
        #newone = type(self)()
        #newone.__dict__.update(self.__dict__)
        #return newone
        cls = self.__class__
        result = cls.__new__(cls)
        memo[id(self)] = result
        bad_attrs = []
        for key, value in self.__dict__.items():
            try:
                memo2 = deepcopy(value, memo)
            except SyntaxError:
                if isinstance(value, dict):
                    for keyi, valuei in value.items():
                        self.log.warn(valuei)
                        #print(valuei.object_attributes())
                        #break
                raise
                bad_attrs.append(key)
            setattr(result, key, memo2)
        if bad_attrs:
            raise RuntimeError(f'failed copying {bad_attrs}')
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

    def _add_disabled_cards(self):
        self._remove_disabled_cards = False
        self.cards_to_read.update(REMOVED_CARDS)  # add

    def __repr__(self) -> str:
        """prints the card count of the model"""
        msg = 'BDF:\n'
        cards = sorted(self._cards_to_setup, key=lambda card: card.type)
        #msg = '\n'.join([card for card in cards if card.n])
        for card in cards:
            if card.n == 0:
                continue
            msg += str(card) + '\n'

        card_count = {}
        singletons = {'AERO', 'AEROS', 'DTABLE', 'DOPTPRM'}
        for card in self._cards_not_setup:
            if card is None:
                continue
            elif isinstance(card, list):
                card_counti = defaultdict(int)
                for cardi in card:
                    card_counti[cardi.type] += 1
                card_counti = dict(card_counti)
                for name, counti in sorted(card_count.items()):
                    card_count[name] = counti

            elif isinstance(card, dict):
                card_counti = defaultdict(int)
                for cardi in card.values():
                    card_counti[cardi.type] += 1
                card_counti = dict(card_counti)
                for name, counti in sorted(card_count.items()):
                    card_count[name] = counti

            elif card.type in singletons:
                card_count[card.type] = 1
            else:
                raise RuntimeError(card)

        if len(card_count):
            msg = 'Vectorized ' + msg + 'Unvectorized:\n'
            for name, counti in card_count.items():
                msg += f'{name}: {counti}\n'
        return msg

def _echo_card(card, card_obj):
    """echos a card"""
    try:
        print(print_card_8(card_obj).rstrip())
    except Exception:
        if card in ['DEQATN']:
            print(str(card_obj).rstrip())
        else:
            print(print_card_16(card_obj).rstrip())

def read_bdf(bdf_filename: Optional[str]=None, validate: bool=True, xref: bool=True, punch: bool=False,
             save_file_structure: bool=False,
             skip_cards: Optional[list[str]]=None,
             read_cards: Optional[list[str]]=None,
             encoding: Optional[str]=None,
             log: Optional[SimpleLogger]=None,
             debug: bool=True, mode: str='msc') -> BDF:
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
    skip_cards : list[str]; default=None
        None : include all cards
        list of cards to skip
    read_cards : list[str]; default=None
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
        model.enable_cards(read_cards)

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

            #'add_AELIST', 'add_AEPARM',
            #'add_DIVERG',
            #'add_FLFACT', ,
            #'add_GUST', 'add_NLPARM', 'add_NLPCI',
            #'add_PARAM', 'add_PHBDY',
            #'add_SET', 'add_SEUSET',

            #'add_card', 'add_card_fields', 'add_card_lines', 'add_cmethod', 'add_constraint',
            #'add_convection_property', 'add_creep_material',
            #'add_material_dependence', 'add_method',
            #'add_random_table',
            #'add_table', 'add_table_sdamping', 'add_thermal_BC', 'add_thermal_element',

            #'set_as_msc',
            #'set_as_nx',

            #'pop_parse_errors',
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

def _check_for_spaces(card_name: str, card_lines: list[str], comment: str, log: SimpleLogger) -> None:
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
    nids_cd_transform = defaultdict(list)  # type: dict[int, np.ndarray]
    nids_cp_transform = defaultdict(list)  # type: dict[int, np.ndarray]
    nxyz = nnodes + nspoints + nepoints + ngridb
    xyz_cp = np.zeros((nxyz, 3), dtype=fdtype)
    nid_cp_cd = np.zeros((nxyz, 3), dtype=idtype)

    if nnodes:
        nids_cp_transform, nids_cd_transform, nid_cp_cd, xyz_cp = model.grid.get_nid_cp_cd_xyz_transforms()

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


#def _get_coords_to_update(coords: dict[int, Union[CORD1R, CORD1C, CORD1S,
                                                  #CORD2R, CORD2C, CORD2S]],
                          #cps_to_check: list[int],
                          #cps_checked: list[int],
                          #nids_checked: list[int]) -> tuple[int, list[int], list[int], list[int]]:
    #"""helper method for ``transform_xyzcp_to_xyz_cid``"""
    #cord1s_to_update_temp = []
    #cord2s_to_update_list = []
    #for cp in sorted(cps_to_check):
        #coord = coords[cp]
        #if coord.type in {'CORD2R', 'CORD2C', 'CORD2S'}:
            #if coord.rid in cps_checked:
                #cord2s_to_update_list.append(cp)
        #elif coord.type in {'CORD1R', 'CORD1C', 'CORD1S'}:
            #cord1s_to_update_temp.append(cp)
        #else:
            #raise NotImplementedError(coord.rstrip())

    #cord1s_to_update_list = []
    #if cord1s_to_update_temp:
        #cord1s_to_update = set()
        #if len(nids_checked) == 0:
            #raise RuntimeError('len(nids_checked)=0...this should not happen.')
        #elif len(nids_checked) == 1:
            #pass
        #else:
            #nids_checked = [np.hstack(nids_checked)]

        #nids_checkedi = nids_checked[0]
        #if len(nids_checkedi) != 0:
            ## check the CORD1x cards
            ##
            ##print('nids_checked = ', nids_checkedi)
            #for cp in cord1s_to_update_temp:
                #coord = coords[cp]
                #nids = coord.node_ids
                ##print('cp=%s nids=%s' % (cp, nids))
                #for nid in nids:
                    #if nid not in nids_checkedi:
                        ##print('  nid=%s break...' % nid)
                        #break
                #else:
                    ##print('  passed')
                    ## all nids passed
                    #cord1s_to_update.add(cp)
            #cord1s_to_update_list = list(cord1s_to_update)
            #cord1s_to_update_list.sort()

    #ncoords = len(cord1s_to_update_list) + len(cord2s_to_update_list)
    ##if ncoords == 0:
        ##msg = 'CPs not handled=%s cord1s_to_update=%s cord2s_to_update=%s\n' % (
            ##cps_to_check, cord1s_to_update, cord2s_to_update)
        ##for cp in (cord1s_to_update + cord2s_to_update):
            ##msg += str(cp)
        ##raise RuntimeError(msg)
    #return ncoords, cord1s_to_update_list, cord2s_to_update_list, nids_checked

class Zona():
    def __init__(self):
        pass
    def update_for_zona(self):
        pass
def map_version(fem: BDF, version: str):
    fem.nastran_format = version
    fem.zona = Zona()

    version_map = {
        'msc': fem.set_as_msc,
        'nx': fem.set_as_nx,
        'optistruct': fem.set_as_optistruct,
        'mystran': fem.set_as_mystran,
        'nasa95': fem.set_as_nasa95,
        'zona': fem.set_as_zona,
    }
    try:
        func = version_map[version]
    except KeyError: # msc, nx, zona, nasa95, mystran
        fmts = ', '.join(version_map)
        msg = f'mode={version!r} is not supported; modes=[{fmts}]'
        raise RuntimeError(msg)
    func()

def map_update(fem: BDF, version: str):
    #if self.nastran_format == 'zona':
        #self.zona.update_for_zona()
    #elif self.nastran_format == 'mystran':
        #self._update_for_mystran()
    #elif self.nastran_format == 'nasa95':
        #self._update_for_nasa95()
    #else:
        # msc / nx / optistruct
        #self._update_for_nastran()

    version_map = {
        'msc': fem._update_for_nastran,
        'nx': fem._update_for_nastran,
        'optistruct': fem._update_for_optistruct,
        'mystran': fem._update_for_mystran,
        'nasa95': fem._update_for_nasa95,
        'zona': fem.zona.update_for_zona,
    }
    try:
        func = version_map[version]
    except KeyError:
        msg = f'mode={version!r} is not supported; modes=[msc, nx, optistruct, zona, nasa95, mystran]'
        raise RuntimeError(msg)
    func()


#if mode == 'msc':
    #self.set_as_msc()
#elif mode == 'nx':
    #self.set_as_nx()
#elif mode == 'nasa95':
    #self.set_as_nasa95()
#elif mode == 'mystran':
    #self.set_as_mystran()
#elif mode == 'zona':
    #self.set_as_zona()
#else:  # pragma: no cover
    #msg = f'mode={self._nastran_format!r} is not supported; modes=[msc, nx, zona, nasa95, mystran]'
    #raise NotImplementedError(msg)


def main():  # pragma: no cover
    """shows off how unicode works because it's overly complicated"""
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

    # The encoding that goes into the comment must be consistent with the local
    # file, so if your print doesn't work right, your comment will be bad too.
    #
    # If the print is correct, assuming all the characters are supported in your
    # desired encoding, it *should* work.
    print(note)

    # Comments are unmodified, so you can inadvertently add cards/bugs.
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
