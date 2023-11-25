from __future__ import annotations
from collections import defaultdict
from itertools import zip_longest
import numpy as np
from typing import TYPE_CHECKING, Set, Optional, Any

from pyNastran.bdf.cards.params import PARAM, MDLPRM
from pyNastran.bdf.cards.dmig import DMI, DMIG, DMIG_UACCEL, DMIAX, DMIJ, DMIJI, DMIK
#from pyNastran.bdf.cards.coordinate_systems import CORD2R
from pyNastran.dev.bdf_vectorized3.cards.bdf_sets import (
    ASET, BSET, CSET, QSET, OMIT, USET, SUPORT,
    SEBSET, SECSET, SEQSET, SESET, RELEASE,
    SET1, SET3, RADSET) # , SET2
from pyNastran.dev.bdf_vectorized3.cards.grid import GRID, EPOINT, SPOINT, GRDSET, POINT
from pyNastran.dev.bdf_vectorized3.cards.elements.rod import CROD, PROD, CONROD, CTUBE, PTUBE
from pyNastran.dev.bdf_vectorized3.cards.elements.bar import BAROR, CBAR, CBARAO, PBAR, PBARL, PBRSECT
from pyNastran.dev.bdf_vectorized3.cards.elements.bush import (
    CBUSH, PBUSH, PBUSHT, CBUSH1D, PBUSH1D, CBUSH2D, PBUSH2D)
from pyNastran.dev.bdf_vectorized3.cards.elements.fast import CFAST, PFAST
from pyNastran.dev.bdf_vectorized3.cards.elements.genel import GENEL
from pyNastran.dev.bdf_vectorized3.cards.elements.spring import CELAS1, CELAS2, CELAS3, CELAS4, PELAS, PELAST
from pyNastran.dev.bdf_vectorized3.cards.elements.damper import (
    CDAMP1, CDAMP2, CDAMP3, CDAMP4, CDAMP5,
    PDAMP, PDAMP5, PDAMPT, CVISC, PVISC, CGAP, PGAP)
from pyNastran.dev.bdf_vectorized3.cards.elements.beam import (
    CBEAM, PBEAM, PBEAML, PBCOMP, # , PBMSECT
    CBEND, PBEND)
from pyNastran.dev.bdf_vectorized3.cards.elements.shear import CSHEAR, PSHEAR
from pyNastran.dev.bdf_vectorized3.cards.elements.shell import (
    CQUAD4, CTRIA3, CQUAD8, CTRIA6, CTRIAR, CQUADR, CQUAD,
    # SNORM,
    CAABSF, # acoustic shells
)
from pyNastran.dev.bdf_vectorized3.cards.elements.shell_properties import (
    PSHELL, PCOMP, PCOMPG,
    PLPLANE, PSHLN1, PSHLN2,
    PAABSF,
)
from pyNastran.dev.bdf_vectorized3.cards.elements.shell_axi import (
    CTRIAX, CTRIAX6,
    CQUADX, CQUADX4, CQUADX8,
    CTRAX3, CTRAX6)

from pyNastran.dev.bdf_vectorized3.cards.elements.plate_stress_strain import (
    PPLANE,
    CPLSTS3, CPLSTS4, CPLSTS6, CPLSTS8,
    CPLSTN3, CPLSTN4, CPLSTN6, CPLSTN8,
)
from pyNastran.dev.bdf_vectorized3.cards.elements.solid import (
    CTETRA, CHEXA, CPENTA, CPYRAM,
    PSOLID, PLSOLID, PCOMPS, PCOMPLS,
    CHEXCZ,
    #CHACAB, CHACBR,
)
from pyNastran.dev.bdf_vectorized3.cards.elements.mass import CONM1, CONM2
from pyNastran.dev.bdf_vectorized3.cards.elements.cmass import PMASS, CMASS1, CMASS2, CMASS3, CMASS4
from pyNastran.dev.bdf_vectorized3.cards.elements.nsm import NSMADD, NSM, NSM1, NSML, NSML1
#from pyNastran.dev.bdf_vectorized3.cards.elements.thermal import CHBDYE, CHBDYP, CHBDYG, CONV, PCONV, CONVM, PCONVM, PHBDY
from pyNastran.dev.bdf_vectorized3.cards.elements.plot import PLOTEL, PLOTEL3, PLOTEL4, PLOTEL6, PLOTEL8

from pyNastran.dev.bdf_vectorized3.cards.loads.static_loads import (
    LOAD, SLOAD,
    FORCE, FORCE1, FORCE2,
    MOMENT, MOMENT1, MOMENT2,
    TEMP, TEMPD, # DTEMP, DTEMP,
    SPCD, DEFORM,
    RFORCE, RFORCE1,
    GRAV, ACCEL, ACCEL1,
)
from pyNastran.dev.bdf_vectorized3.cards.loads.static_pressure_loads import (
    PLOAD, PLOAD1, PLOAD2, PLOAD4, # PLOADX1,
)
from pyNastran.dev.bdf_vectorized3.cards.loads.types import Loads as StaticLoad

from pyNastran.dev.bdf_vectorized3.cards.loads.dynamic_loads import (
    LSEQ, DLOAD,
    DAREA, TLOAD1, TLOAD2, RLOAD1, RLOAD2, TIC, QVECT,
    RANDPS, DELAY, DPHASE, TF, ACSRCE,
)

from pyNastran.dev.bdf_vectorized3.cards.loads.thermal_loads import (
    QHBDY, QBDY1, QBDY2, QBDY3, QVOL, TEMPBC, RADBC, RADM)

from pyNastran.bdf.cards.materials import NXSTRAT

from pyNastran.dev.bdf_vectorized3.cards.materials import (
    MAT1,
    MAT2, MAT3, MAT4, MAT5,
    MAT8, MAT9, MAT10, MAT11,
    MAT10C, MATORT, MATHP, MATHE,
)
from pyNastran.dev.bdf_vectorized3.cards.materials_dep import (
    MATS1, # MATS3,
    #MATS8, # MSC
    MATT1, MATT2, MATT8, MATT9)

from pyNastran.dev.bdf_vectorized3.cards.contact import (
    BSURF, BSURFS, BCPROP, BCPROPS,
    BGSET, BCTSET, BGADD, BCTADD,
    BCONP, BFRIC, BLSEG, BCRPARA, BEDGE)
from pyNastran.dev.bdf_vectorized3.cards.coord import COORD
from pyNastran.dev.bdf_vectorized3.cards.constraints import (
    SPC, SPC1, SPCADD,
    SPCOFF,  # SPCOFF1,
    BNDFIX,  # BNDFIX1,
    BNDFREE, # BNDFREE1
    MPC, MPCADD)
from pyNastran.dev.bdf_vectorized3.cards.elements.rigid import (
    RBAR, RBAR1,
    RROD,
    RBE1, RBE2, RBE3,
    #RSSCON,
)
from pyNastran.dev.bdf_vectorized3.cards.monitor import MONPNT1, MONPNT2, MONPNT3
#from pyNastran.dev.bdf_vectorized3.cards.aero.aero import (
    #CAERO1, CAERO2, CAERO3, CAERO4, CAERO5, CAERO7,
    #PAERO1, PAERO2, PAERO3, PAERO4, PAERO5,
    #SPLINE1, SPLINE2, SPLINE3, SPLINE4, SPLINE5,
    #AECOMP, AECOMPL, AELIST, AEFACT, FLFACT, AEPARM, AELINK, AESTAT,
    #GUST, AESURF, AESURFS, CSSCHD, TRIM)
from pyNastran.dev.bdf_vectorized3.cards.optimization import (
    DESVAR, DLINK, DVGRID,
    DRESP1, DRESP2, DCONSTR,
    DVPREL1, DVPREL2,
    DVMREL1, DVMREL2,
    DVCREL1, DVCREL2, DCONADD,
    DSCREEN, # DDVAL,
    MODTRAK,
)
from .loads_summation import (
    get_static_loads_by_subcase_id,
    get_reduced_static_load,
    sum_forces_moments)
from .breakdowns import (
    get_mass_breakdown, get_length_breakdown,
    get_area_breakdown, get_volume_breakdown,
    NO_LENGTH, NO_AREA, NO_VOLUME, NO_MASS,
)


if TYPE_CHECKING:  # pragma: no cover
    #from pyNastran.dev.bdf_vectorized3.bdf import PARAM, MDLPRM, FLUTTER
    from pyNastran.bdf.case_control_deck import CaseControlDeck

class BDFAttributes:
    def __init__(self):
        # basic settings
        self.filter_midside_nodes = True

        # DESVAR xinit being out of range of [xlb, xub] causes nastran to crash
        # just push it to the xlb/xub
        self.apply_clip_to_desvar_range = True

        # run the area, length, volume, mass, inertia, and quality methods
        self.run_testing_checks = False

        self.fdtype = 'float64'
        self.idtype = 'int32'
        self.punch = None
        self._encoding = None
        self.save_file_structure = False

        #: ignore any ECHOON flags
        self.force_echo_off = True


        self.initial_superelement_models = []
        self._crash_cards = {
            'CGEN', 'ADAPT', 'FEEDGE', 'FEFACE',
            'GMCURV', 'GMLOAD', 'GMSPC', 'GMSURF', 'PVAL', }


        # ------------------------ structural defaults -----------------------
        #: the analysis type
        self._sol = None
        #: used in solution 600, method
        self.sol_method = None
        #: the line with SOL on it, marks ???
        self.sol_iline: Optional[int] = None
        self.case_control_deck: Optional[CaseControlDeck] = None
        self.app = ''

        # ---------------------------------------------------------------------
        #: store the PARAM cards
        self.params: dict[str, PARAM] = {}
        self.mdlprm: Optional[MDLPRM] = None

        self.grdset: Optional[GRDSET] = None # GRDSET(cp=0, cd=0, ps=0, seid=0, comment='')
        self.grid = GRID(self)
        self.spoint = SPOINT(self)
        self.epoint = EPOINT(self)
        self.point = POINT(self)
        self.coord = COORD(self)

        # plot
        self.plotel = PLOTEL(self)
        self.plotel3 = PLOTEL3(self)
        self.plotel4 = PLOTEL4(self)
        self.plotel6 = PLOTEL6(self)
        self.plotel8 = PLOTEL8(self)

        # spring
        self.celas1 = CELAS1(self)
        self.celas2 = CELAS2(self)
        self.celas3 = CELAS3(self)
        self.celas4 = CELAS4(self)
        self.pelas = PELAS(self)
        self.pelast = PELAST(self)

        # damper
        self.cdamp1 = CDAMP1(self)
        self.cdamp2 = CDAMP2(self)
        self.cdamp3 = CDAMP3(self)
        self.cdamp4 = CDAMP4(self)
        self.cdamp5 = CDAMP5(self)
        self.pdamp = PDAMP(self)
        self.pdamp5 = PDAMP5(self)
        self.pdampt = PDAMPT(self)

        # contact
        self.bsurf = BSURF(self)      # source/target_id to shell eid
        self.bsurfs = BSURFS(self)    # source/target_id to solid eid
        self.bcprop = BCPROP(self)    # source/target_id to shell pid
        self.bcprops = BCPROPS(self)  # source/target_id to solid pid
        self.blseg = BLSEG(self)
        self.bcrpara = BCRPARA(self)
        self.bedge = BEDGE(self)
        self.bctpara = {}
        self.bctparm = {}

        self.bconp = BCONP(self)      # sideline contact
        self.bfric = BFRIC(self)      # friction for BCONP

        self.bgset = BGSET(self)    # glue_id to source/target_id
        self.bctset = BCTSET(self)  # contact_id to source/target_ids
        self.bgadd = BGADD(self)    # add glue_ids
        self.bctadd = BCTADD(self)  # add contact_ids

        # sets
        self.suport = SUPORT(self)
        self.aset = ASET(self)  # ASET, ASET1
        self.bset = BSET(self)  # BSET, BSET1
        self.cset = CSET(self)  # CSET, CSET1
        self.qset = QSET(self)  # QSET, QSET1
        self.omit = OMIT(self)  # OMIT, OMIT1
        self.uset = USET(self)  # USET, USET1

        self.seset = SESET(self)
        self.sebset = SEBSET(self)
        self.secset = SECSET(self)
        self.seqset = SEQSET(self)
        self.release = RELEASE(self)

        self.set1 = SET1(self)
        #self.set2 = SET2(self)
        self.set3 = SET3(self)

        self.nxstrats: dict[int, NXSTRAT] = {}

        # monitor
        self.monpnt1 = MONPNT1(self)
        self.monpnt2 = MONPNT2(self)
        self.monpnt3 = MONPNT3(self)

        # aero geometry
        #self.caero1 = CAERO1(self)
        #self.caero2 = CAERO2(self)
        #self.caero3 = CAERO3(self)
        #self.caero4 = CAERO4(self)
        #self.caero5 = CAERO5(self)
        #self.caero7 = CAERO7(self)  # zona

        #self.paero1 = PAERO1(self)
        #self.paero2 = PAERO2(self)
        #self.paero3 = PAERO3(self)
        #self.paero4 = PAERO4(self)
        #self.paero5 = PAERO5(self)

        #self.spline1 = SPLINE1(self)
        #self.spline2 = SPLINE2(self)
        #self.spline3 = SPLINE3(self)
        #self.spline4 = SPLINE4(self)
        #self.spline5 = SPLINE5(self)  #  faked

        #self.aefact = AEFACT(self)  # caero fractions
        #self.aelist = AELIST(self)  # box ids for control surfaces
        #self.aecomp = AECOMP(self)  # for monitor points
        #self.aecompl = AECOMPL(self)  # for monitor points
        #self.aeparm = AEPARM(self)
        #self.aelink = AELINK(self)  # links control surfaces
        #self.aestat = AESTAT(self)  # degrees of freedom

        # flutter
        #self.flfact = FLFACT(self)  # Mach, Vel, rho for FLUTTER

        #  gust
        #self.gust = GUST(self)

        #-------------------------------------------------
        # visc
        self.cvisc = CVISC(self)
        self.pvisc = PVISC(self)

        # gap
        self.pgap = PGAP(self)
        self.cgap = CGAP(self)

        # rod
        self.crod = CROD(self)
        self.prod = PROD(self)
        self.conrod = CONROD(self)

        # tube
        self.ctube = CTUBE(self)
        self.ptube = PTUBE(self)

        # bush
        self.cbush = CBUSH(self)
        self.pbush = PBUSH(self)
        self.pbusht = PBUSHT(self)

        self.cbush1d = CBUSH1D(self)
        self.pbush1d = PBUSH1D(self)
        #self.cbush2d = CBUSH2D(self)
        #self.pbush2d = PBUSH2D(self)

        # fast
        self.cfast = CFAST(self)
        self.pfast = PFAST(self)

        # genel
        self.genel = GENEL(self)

        self.baror: Optional[BAROR] = None
        self.cbar = CBAR(self)
        self.cbarao = CBARAO(self)
        self.pbar = PBAR(self)
        self.pbarl = PBARL(self)
        self.pbrsect = PBRSECT(self)
        #self.bar_properties = [self.pbar, self.pbarl, self.pbrsect]

        self.beamor = None
        self.cbeam = CBEAM(self)
        self.pbeam = PBEAM(self)
        self.pbeaml = PBEAML(self)
        self.pbcomp = PBCOMP(self)
        #self.beam_properties = [self.pbeaml]

        self.cbend = CBEND(self)
        self.pbend = PBEND(self)

        self.cshear = CSHEAR(self)
        self.pshear = PSHEAR(self)

        #self.snorm = SNORM(self)
        self.ctria3 = CTRIA3(self)
        self.cquad4 = CQUAD4(self)
        self.ctria6 = CTRIA6(self)
        self.cquad8 = CQUAD8(self)
        self.ctriar = CTRIAR(self)
        self.cquadr = CQUADR(self)
        self.cquad = CQUAD(self)
        self.pshell = PSHELL(self)
        self.pcomp = PCOMP(self)
        self.pcompg = PCOMPG(self)
        self.pshln1 = PSHLN1(self)
        self.pshln2 = PSHLN2(self)
        #self.shell_properties = [self.pshell, self.pcomp, self.pcompg]

        # planar shells
        self.plplane = PLPLANE(self)
        self.pplane = PPLANE(self)

        # plane stress
        # what are the properties?
        self.cplsts3 = CPLSTS3(self)
        self.cplsts4 = CPLSTS4(self)
        self.cplsts6 = CPLSTS6(self)
        self.cplsts8 = CPLSTS8(self)

        # plate strain
        # what are the properties?
        self.cplstn3 = CPLSTN3(self)
        self.cplstn4 = CPLSTN4(self)
        self.cplstn6 = CPLSTN6(self)
        self.cplstn8 = CPLSTN8(self)

        # axisymmetric shells
        self.cquadx = CQUADX(self)
        self.cquadx4 = CQUADX4(self)
        self.cquadx8 = CQUADX8(self)
        self.ctriax = CTRIAX(self)
        self.ctriax6 = CTRIAX6(self)

        self.ctrax3 = CTRAX3(self)
        self.ctrax6 = CTRAX6(self)

        # acoustic shells?
        self.caabsf = CAABSF(self)
        self.paabsf = PAABSF(self)

        # solid elements
        self.ctetra = CTETRA(self)
        self.chexa = CHEXA(self)
        self.cpenta = CPENTA(self)
        self.cpyram = CPYRAM(self)
        self.psolid = PSOLID(self)
        self.plsolid = PLSOLID(self)
        self.pcomps = PCOMPS(self)
        self.pcompls = PCOMPLS(self)
        #self.solid_properties = []

        self.chexcz = CHEXCZ(self)

        # mass
        self.pmass = PMASS(self)
        self.cmass1 = CMASS1(self)
        self.cmass2 = CMASS2(self)
        self.cmass3 = CMASS3(self)
        self.cmass4 = CMASS4(self)
        self.conm1 = CONM1(self)
        self.conm2 = CONM2(self)

        # nonstructural mass
        self.nsmadd = NSMADD(self)
        self.nsm = NSM(self)
        self.nsm1 = NSM1(self)
        # lumped
        self.nsml = NSML(self)
        self.nsml1 = NSML1(self)

        # thermal
        #self.bdyor = None
        #self.chbdye = CHBDYE(self)
        #self.chbdyp = CHBDYP(self)
        #self.chbdyg = CHBDYG(self)
        #self.phbdy = PHBDY(self)

        #self.chacab = CHACAB(self)
        #self.chacbr = CHACBR(self)

        # thermal boundary conditions
        #self.conv = CONV(self)
        #self.pconv = PCONV(self)
        #self.convm = CONVM(self)
        #self.pconvm = PCONVM(self)

        # loads
        self.load = LOAD(self)
        self.sload = SLOAD(self)
        self.force = FORCE(self)
        self.force1 = FORCE1(self)
        self.force2 = FORCE2(self)
        self.moment = MOMENT(self)
        self.moment1 = MOMENT1(self)
        self.moment2 = MOMENT2(self)
        self.pload = PLOAD(self)
        self.pload1 = PLOAD1(self)
        self.pload2 = PLOAD2(self)
        self.pload4 = PLOAD4(self)
        self.grav = GRAV(self)
        self.accel = ACCEL(self)
        self.accel1 = ACCEL1(self)
        self.temp = TEMP(self)
        self.tempd = TEMPD(self)  # default temp
        #self.dtemp = DTEMP(self)  # has nodes
        self.lseq = LSEQ(self)    # static load sequence
        self.spcd = SPCD(self)   # enforced displacment; load
        self.deform = DEFORM(self) # encforced displacement; load
        self.rforce = RFORCE(self)    # rotational force
        self.rforce1 = RFORCE1(self)  # rotational force

        # axisymmetric loads
        #self.ploadx1 = PLOADX1(self)

        # other thermal...
        #self.dtemp = DTEMP(self)
        self.qvect = QVECT(self)
        self.qhbdy = QHBDY(self)
        self.qbdy1 = QBDY1(self)
        self.qbdy2 = QBDY2(self)
        self.qbdy3 = QBDY3(self)
        self.qvol = QVOL(self)
        self.tempbc = TEMPBC(self)
        self.radbc = RADBC(self)
        self.radm = RADM(self)

        # radiation
        self.radset = RADSET(self)

        # dynamic loads
        self.dload = DLOAD(self)
        self.darea = DAREA(self)
        self.tload1 = TLOAD1(self)
        self.tload2 = TLOAD2(self)
        self.rload1 = RLOAD1(self)
        self.rload2 = RLOAD2(self)
        self.tic = TIC(self)  # initial conditions
        self.tf = TF(self)  # transfer function

        # dynamic load offset
        self.dphase = DPHASE(self)
        self.delay = DELAY(self)

        # random loads
        self.randps = RANDPS(self)
        self.acsrce = ACSRCE(self)

        self.mat1 = MAT1(self)
        self.mat2 = MAT2(self)
        self.mat3 = MAT3(self)
        self.mat4 = MAT4(self)
        self.mat5 = MAT5(self)
        self.mat8 = MAT8(self)
        self.mat9 = MAT9(self)
        self.mat10 = MAT10(self)
        self.mat11 = MAT11(self)
        self.mat10c = MAT10C(self)
        self.matort = MATORT(self)
        self.mathe = MATHE(self)
        self.mathp = MATHP(self)

        self.mats1 = MATS1(self)
        #self.mats3 = MATS3(self)
        #self.mats8 = MATS8(self)
        self.matt1 = MATT1(self)
        self.matt2 = MATT2(self)
        self.matt8 = MATT8(self)
        self.matt9 = MATT9(self)

        self.spc = SPC(self)
        self.spc1 = SPC1(self)
        self.spcoff = SPCOFF(self)  # SPCOFF1
        self.bndfix = BNDFIX(self)  # BNDFIX1
        self.bndfree = BNDFREE(self)  # BNDFREE1
        self.spcadd = SPCADD(self)
        self.mpc = MPC(self)
        self.mpcadd = MPCADD(self)

        # rigid elements
        self.rbar = RBAR(self)
        self.rbar1 = RBAR1(self)
        self.rbe1 = RBE1(self)
        self.rbe2 = RBE2(self)
        self.rbe3 = RBE3(self)
        self.rrod = RROD(self)
        #self.rsscon = RSSCON(self)  # not supported

        # modes
        self.methods = {}
        self.cMethods = {}

        # nonlinear
        self.nlparms = {}

        # transient
        self.tsteps = {}

        # frequency solutions
        self.frequencies = {}

        # nonlinear transient
        self.tstepnls = {}
        self.nlpcis = {}

        # optimization
        self.modtrak = MODTRAK(self)
        self.desvar = DESVAR(self)
        self.dlink = DLINK(self)
        self.dvgrid = DVGRID(self)
        self.dresp1 = DRESP1(self)
        self.dresp2 = DRESP2(self)
        self.dconstr = DCONSTR(self)

        self.dvprel1 = DVPREL1(self)
        self.dvprel2 = DVPREL2(self)
        self.dvmrel1 = DVMREL1(self)
        self.dvmrel2 = DVMREL2(self)
        self.dvcrel1 = DVCREL1(self)
        self.dvcrel2 = DVCREL2(self)

        self.dconadd = DCONADD(self)
        #self.ddval = DDVAL(self)
        self.dscreen = DSCREEN(self)

        # ---------------------------------------------------
        ## unoptimized
        # aero model
        #self.aefacts = {}

        #self.aesurf = AESURF(self)
        #self.aesurfs = AESURFS(self)
        #self.csschd = CSSCHD(self)
        #self.trim = TRIM(self)
        #self.trim2 = TRIM(self)

        # static aero
        self.aeros = None
        self.trims = {}
        self.divergs = {}

        #  flutter
        self.aero = None
        self.flutters = {}
        self.mkaeros = []

        # ----------------------------------------------------------------
        # matrices
        #: direct matrix input - DMIG
        self.dmi: dict[str, DMI]= {}
        self.dmig: dict[str, DMIG | DMIG_UACCEL] = {}
        self.dmij: dict[str, DMIJ] = {}
        self.dmiji: dict[str, DMIJI] = {}
        self.dmik: dict[str, DMIK] = {}
        self.dmiax: dict[str, DMIAX] = {}
        self.dti: dict[str, DMI] = {}
        self._dmig_temp = defaultdict(list)  # type: dict[str, list[str]]
        # ----------------------------------------
        #self.suport1 = {}
        #self.suport = []
        #self.suport = SUPORT(self)

        self.system_command_lines: list[str] = []
        self.executive_control_lines: list[str] = []
        self.superelement_models = {}

        #origin = np.array([0., 0., 0.])
        #zaxis = np.array([0., 0., 1.])
        #xzplane = np.array([1., 0., 0.])
        #coord = CORD2R(cid=0, rid=0, origin=origin, zaxis=zaxis, xzplane=xzplane)
        #self.coords = {0 : coord}   # type: dict[int, Any]

        ## old style
        self.ringaxs = {}
        self.gridb = {}
        self.mdlprm = None

        # optimization
        self.doptprm = None

    def aero_coord(self) -> int:
        """gets the aerodynamic coordinate system"""
        aero = self.aero
        aeros = self.aeros
        if aero is None and aeros is None:
            msg = 'neither AERO nor AEROS cards exist'
            #raise RuntimeError(msg)
            self.log.warning(msg)
            return None

        if aero is not None and aeros is not None:
            acsid_aero = aero.Acsid()
            acsid_aeros = aeros.Acsid()
            assert acsid_aero == acsid_aeros, f'AERO acsid={acsid_aero}, AEROS acsid={acsid_aeros}'
            coord = acsid_aeros
        elif aero is not None:
            coord = aero.Acsid()
        elif aeros is not None:
            coord = aeros.Acsid()
        return coord

    @property
    def plot_element_cards(self) -> list[Any]:
        elements = [
            self.plotel,
            self.plotel3, self.plotel4,
            self.plotel6, self.plotel8,
        ]
        return elements

    @property
    def spring_element_cards(self) -> list[Any]:
        elements = [
            self.celas1, self.celas2, self.celas3, self.celas4,
        ]
        return elements

    @property
    def damper_element_cards(self) -> list[Any]:
        elements = [
            self.cdamp1, self.cdamp2, self.cdamp3, self.cdamp4, self.cdamp5,
        ]
        return elements

    @property
    def shell_element_cards(self) -> list[Any]:
        elements = [
            self.ctria3, self.cquad4, self.ctria6, self.cquad8,
            self.ctriar, self.cquadr, self.cquad,
        ]
        return elements

    @property
    def solid_element_cards(self) -> list[Any]:
        elements = [
            self.ctetra, self.cpenta, self.chexa, self.cpyram,
            self.chexcz,
        ]
        return elements

    @property
    def element_cards(self) -> list[Any]:
        axisymmetric_element_cards = [
            #self.ctriax,
            #self.cquadx, self.cquadx4, self.cquadx8,
            #self.ctriax, self.ctriax6,
            #self.ctrax3, self.ctrax6,
        ]
        acoustic_element_cards = [
            #self.chacab, self.chacbr,
        ]
        elements = self.spring_element_cards + self.damper_element_cards + [
            self.cvisc, self.cgap,
            self.cbush, self.cbush1d, # self.cbush2d,
            self.cfast,
            self.crod, self.conrod, self.ctube,
            self.cbar,
            self.cbeam,
            self.cbend,
            self.cshear,
            self.caabsf, # acoustic shells
            self.genel,
        ] + self.shell_element_cards + self.solid_element_cards + \
        axisymmetric_element_cards + [
            self.conm1, self.conm2,
            self.cmass1, self.cmass2, self.cmass3, self.cmass4,
            self.cplsts3, self.cplsts4, self.cplsts6, self.cplsts8,
            self.cplstn3, self.cplstn4, self.cplstn6, self.cplstn8,
        ] + acoustic_element_cards
        return elements

    @property
    def bar_property_cards(self) -> list[Any]:
        properties = [
            self.pbar, self.pbarl, self.pbrsect,
        ]
        return properties

    @property
    def beam_property_cards(self) -> list[Any]:
        properties = [
            self.pbeam, self.pbeaml, self.pbcomp,
        ]
        return properties

    @property
    def line_property_cards(self) -> list[Any]:
        properties = self.bar_property_cards + self.beam_property_cards
        return properties

    @property
    def shell_property_cards(self) -> list[Any]:
        properties = [
            self.pshell, self.pcomp, self.pcompg,
            self.plplane, self.pplane,
            # self.pshln1, self.pshln2,
        ]
        return properties
    @property
    def nonlinear_shell_property_cards(self) -> list[Any]:
        properties = [
            self.pshln1, self.pshln2,
        ]
        return properties

    @property
    def solid_property_cards(self) -> list[Any]:
        properties = [
            self.psolid, self.plsolid, self.pcomps, self.pcompls,
        ]
        return properties

    @property
    def property_cards(self) -> list[Any]:
        acoustic_property_cards = [self.paabsf]
        zero_d_properties = [
            self.pelas, self.pelast,
            self.pdamp, self.pdampt,
            self.pbush, self.pbusht,
            self.pbush1d, # self.pbush2d,
            self.pfast,
            self.pvisc, self.pgap,
            self.pmass,
        ]
        line_properties = [
            self.prod, self.ptube, self.pbend,
        ] + self.bar_property_cards + self.beam_property_cards
        shell_properties = [self.pshear] + \
            self.shell_property_cards + self.nonlinear_shell_property_cards
        properties = zero_d_properties + line_properties + shell_properties + \
            self.solid_property_cards + acoustic_property_cards
        return properties

    @property
    def material_cards(self) -> list[Any]:
        materials = (
            self.structural_material_cards + self.thermal_material_cards +
            self.hyperelastic_material_cards + self.thermal_material_dep_cards)
        return materials

    @property
    def structural_material_cards(self) -> list[Any]:
        materials = [
            self.mat1, self.mat2, self.mat3,
            self.mat8, self.mat9, self.mat10, self.mat11,
            self.mat10c, self.matort,
        ]
        return materials

    @property
    def thermal_material_dep_cards(self) -> list[Any]:
        materials = [
            self.matt1, self.matt8, self.matt9,
            self.mats1,
        ]
        return materials

    @property
    def thermal_material_cards(self) -> list[Any]:
        materials = [
            self.mat4, self.mat5,
        ]
        return materials

    @property
    def hyperelastic_material_cards(self) -> list[Any]:
        materials = [
            self.mathe,
            self.mathp,
        ]
        return materials

    @property
    def optimization_cards(self) -> list[Any]:
        optimization = [
            self.desvar, self.dlink, self.dvgrid,
            self.dresp1, self.dresp2, self.dconstr, self.dconadd,
            self.dvprel1, self.dvprel2,
            self.dvmrel1, self.dvmrel2,
            self.dvcrel1, self.dvcrel2,
            self.dscreen,
            self.modtrak,
        ]
        return optimization

    @property
    def load_cards(self) -> list[Any]:
        loads = [
            self.load, self.lseq,
            self.force, self.force1, self.force2,
            self.moment, self.moment1, self.moment2,
            self.pload, self.pload1, self.pload2, self.pload4,
            self.grav, self.accel, self.accel1,
            self.sload,
            self.temp, self.tempd,
            #self.dtemp, # has nodes
            self.qhbdy, self.qbdy1, self.qbdy2, self.qbdy3,
            self.qvol, self.qvect,
            self.spcd, self.deform,
            self.rforce, self.rforce1,
            #self.ploadx1,
        ]
        return loads

    @property
    def dynamic_load_cards(self) -> list[Any]:
        loads = [
            self.dload,
            self.darea,
            self.tload1, self.tload2, # self.tload1,
            self.rload1, self.rload2, # self.rloadex

            # random loads
            self.randps, self.acsrce,
        ]
        return loads

    @property
    def dynamic_cards(self) -> list[Any]:
        cards = [
            self.tic, self.delay, self.dphase, self.tf,
        ]
        return cards

    @property
    def rigid_element_cards(self) -> list[Any]:
        rigid_elements = [
            self.rbar, self.rbar1,
            self.rrod,
            self.rbe1, self.rbe2, self.rbe3,
            #self.rsscon,
        ]
        return rigid_elements

    @property
    def nonstructural_mass_cards(self) -> list[Any]:
        cards = [
            self.nsmadd,
            self.nsm, self.nsm1,
            self.nsml, self.nsml1, # lumped
        ]
        return cards

    @property
    def thermal_element_cards(self) -> list[Any]:
        thermal_elements = [
            #self.chbdye, self.chbdyg, self.chbdyp, self.phbdy,
        ]
        return thermal_elements

    @property
    def thermal_boundary_condition_cards(self) -> list[Any]:
        boundary_conditions = [
            #self.conv, self.pconv,
            #self.convm, self.pconvm,
            #self.tempbc,
            self.radbc, self.radm, self.radset,
        ]
        return boundary_conditions

    @property
    def set_cards(self) -> list[Any]:
        """handles ASET/ASET1, ..."""
        sets = [
            self.set1, # self.set2,
            self.set3,
            self.aset, self.bset, self.cset, self.qset,
            self.omit, self.uset,
        ]
        return sets

    @property
    def seset_cards(self) -> list[Any]:
        """handles SEBSET/SEBSET1, ..."""
        sesets = [
            self.sebset, self.secset, self.seqset, self.release,
            self.seset,
            #self.seuset,
        ]
        return sesets

    @property
    def aero_elements(self) -> list[Any]:
        elements = [
            #self.caero1, self.caero2, self.caero3, self.caero4, self.caero5,
            #self.caero7,
        ]
        return elements
    @property
    def aero_properties(self) -> list[Any]:
        properties = [
            #self.paero1, self.paero2, self.paero3, self.paero4, self.paero5,
        ]
        return properties
    @property
    def aero_splines(self) -> list[Any]:
        splines = [
            #self.spline1, self.spline2, self.spline3,
            #self.spline4, self.spline5,
        ]
        return splines

    @property
    def aero_loads(self) -> list[Any]:
        loads = [
            #self.gust, self.csschd, self.trim, self.trim2, # self.diverg,
        ]
        return loads

    @property
    def aero_objects(self) -> list[Any]:
        aero = [
            #self.aecomp, self.aecompl, self.aesurf, self.aesurfs, self.aestat,
            #self.aelist, self.aeparm, self.aelink,
            #self.aefact, self.flfact,
            ] + self.aero_elements + \
            self.aero_properties + self.aero_splines + self.monitor_point_cards + \
            self.aero_loads
        return aero

    @property
    def monitor_point_cards(self) -> list[Any]:
        monitor_points = [
            self.monpnt1, self.monpnt2, self.monpnt3,
        ]
        return monitor_points

    # ------------------------------------------------------------------------
    # constraints
    @property
    def spc_cards(self) -> list[Any]:
        return [
            self.spc, self.spc1, self.spcadd, self.spcoff, self.bndfix, self.bndfree,
        ]
    @property
    def mpc_cards(self) -> list[Any]:
        return [
            self.mpc, self.mpcadd, # self.mpcax,
        ]
    @property
    def nsm_cards(self) -> list[Any]:
        return [
            self.nsm, self.nsm1, self.nsml, self.nsml1, self.nsmadd,
        ]

    @property
    def contact_cards(self) -> list[Any]:
        return [self.bsurf, self.bsurfs, self.bcprop, self.bcprops,
                self.bgset, self.bctset,
                self.bgadd, self.bctadd,
                self.bconp, self.bfric, self.blseg, self.bcrpara, self.bedge]

    @property
    def _cards_to_setup(self) -> list[Any]:
        cards = [
            self.grid, self.spoint, self.epoint, self.point,
            self.coord,
            #self.snorm,
            self.suport, # self.suport1,
            self.cbarao,
        ] + self.spc_cards + self.mpc_cards + self.element_cards + self.rigid_element_cards + \
        self.property_cards + self.material_cards + self.optimization_cards + \
        self.load_cards + self.dynamic_load_cards + \
        self.plot_element_cards + self.thermal_element_cards + \
        self.thermal_boundary_condition_cards + self.set_cards + self.seset_cards + \
        self.aero_objects + self.dynamic_cards + self.nonstructural_mass_cards + self.contact_cards
        #for i, card in enumerate(cards):
            #assert isinstance(card.type, str), card
            #print(card.type)
        return cards
    # ------------------------------------------------------------------------
    @property
    def spring_element_ids(self) -> np.ndarray:
        elements = self.spring_element_cards
        element_ids = check_element_ids(elements)
        return element_ids

    @property
    def damper_element_ids(self) -> np.ndarray:
        elements = self.damper_element_cards
        element_ids = check_element_ids(elements)
        return element_ids

    @property
    def shell_element_ids(self) -> np.ndarray:
        elements = self.shell_element_cards
        element_ids = check_element_ids(elements)
        return element_ids
    @property
    def solid_element_ids(self) -> np.ndarray:
        elements = self.solid_element_cards
        element_ids = check_element_ids(elements)
        return element_ids
    # ------------------------------------------------------------------------
    @property
    def bar_property_ids(self) -> np.ndarray:
        properties = self.bar_property_cards
        property_ids = check_property_ids(properties, 'Bar ')
        return property_ids

    @property
    def beam_property_ids(self) -> np.ndarray:
        properties = self.beam_property_cards
        property_ids = check_property_ids(properties, 'Beam ')
        return property_ids

    @property
    def shell_property_ids(self) -> np.ndarray:
        properties = self.shell_property_cards
        property_ids = check_property_ids(properties, 'Shell ')
        return property_ids

    @property
    def solid_property_ids(self) -> np.ndarray:
        properties = self.solid_property_cards
        property_ids = check_property_ids(properties, 'Solid ')
        return property_ids

    @property
    def spline_ids(self) -> np.ndarray:
        list_spline_ids = [caero.spline_id for caero in self.aero_splines]
        if len(list_spline_ids) == 0:
            return np.array([], dtype='int32')
        spline_ids = np.hstack(list_spline_ids)
        uspline_ids = np.unique(spline_ids)
        if len(spline_ids) != len(uspline_ids):
            count = np.bincount(spline_ids)
            count_where = np.where(count > 1)[0]
            count2 = count[count_where]
            msg = ''
            raise RuntimeError(f'Duplicate {msg}SPLINEx IDs\n'
                               f'spline_ids={count_where}\ncount={count2}')
        return uspline_ids

    @property
    def caero_ids(self) -> np.ndarray:
        list_caero_ids = [caero.element_id for caero in self.aero_elements]
        if len(list_caero_ids) == 0:
            return np.array([], dtype='int32')
        caero_ids = np.hstack(list_caero_ids)
        ucaero_ids = np.unique(caero_ids)
        if len(caero_ids) != len(ucaero_ids):
            count = np.bincount(caero_ids)
            count_where = np.where(count > 1)[0]
            count2 = count[count_where]
            msg = ''
            raise RuntimeError(f'Duplicate {msg}CAEROx IDs\n'
                               f'caero_ids={count_where}\ncount={count2}')
        return ucaero_ids
    # ------------------------------------------------------------------------

    def set_as_msc(self) -> None:
        self._nastran_format = 'msc'
        self.is_nx = False
        self.is_msc = True
        self.is_optistruct = False
        self.is_mystran = False
        self.is_nasa95 = False
        self.is_zona = False

    def set_as_nx(self) -> None:
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

    def set_as_mystran(self) -> None:
        self._nastran_format = 'mystran'
        self.is_nx = False
        self.is_msc = False
        self.is_optistruct = False
        self.is_mystran = True
        self.is_nasa95 = False
        self.is_zona = False
        self._update_for_mystran()

    def set_as_nasa95(self) -> None:
        self._nastran_format = 'nasa95'
        self.is_nx = False
        self.is_msc = False
        self.is_optistruct = False
        self.is_mystran = False
        self.is_nasa95 = True
        self.is_zona = False
        self._update_for_nasa95()

    # ------------------------------------------------------------------------
    # optimization
    @property
    def dresp_ids(self) -> np.ndarray:
        return np.hstack([self.dresp1.dresp_id, self.dresp2.dresp_id])

    @property
    def is_thermal(self) -> bool:
        if self.app == 'HEAT':
            return True
        is_thermal = False
        subcases = self.subcases
        for subcase_id, subcase in subcases.items():
            if 'ANALYSIS' in subcase:
                value, options = subcase['ANALYSIS']
                if value in {'HEAT', 'HSTAT', 'HTRAN'}:
                    is_thermal = True
        return is_thermal

    def self_static_loads_by_subcase_id(self, subcase_ids: list[int]=None) -> dict[int, Any]:
        res = get_static_loads_by_subcase_id(self, subcase_ids=subcase_ids)
        return res

    def get_reduced_static_load(self, load: Optional[StaticLoad]=None):
        reduced_loads = get_reduced_static_load(self, load=load)
        return reduced_loads

    def sum_forces_moments(self) -> tuple[dict[int, np.ndarray],
                                          dict[int, np.ndarray],]:
        loads_by_load_id, loads_by_subcase_id = sum_forces_moments(self, [0., 0., 0.])
        return loads_by_load_id, loads_by_subcase_id

    def get_mass_breakdown(self, property_ids=None,
                           stop_if_no_mass: bool=True) -> tuple[dict[int, float], dict[str, float]]:
        """
        Returns
        -------
        pids_to_mass : dict {int : float, ...}
            Map from property id to mass.
        mass_type_to_mass : dict {str : float, ...}
            Map from mass id to mass for mass elements.

        """
        pid_to_mass, mass_type_to_mass = get_mass_breakdown(
            self,
            stop_if_no_mass=stop_if_no_mass)
        return pid_to_mass, mass_type_to_mass

    def get_length_breakdown(self, property_ids=None,
                             stop_if_no_length: bool=True):
        """
        Gets a breakdown of the length by property region.

        Returns
        -------
        pids_to_length : dict[int pid] : float length
            the pid to length dictionary

        TODO: What about CONRODs?

        """
        pids_to_length = get_length_breakdown(
            self,
            property_ids=property_ids,
            stop_if_no_length=stop_if_no_length)
        return pids_to_length

    def get_area_breakdown(self,
                           property_ids=None,
                           stop_if_no_area: bool=True,
                           sum_bar_area: bool=True) -> dict[int, float]:
        """
        Gets a breakdown of the area by property region.

        Parameters
        ----------
        property_ids : list[int] / int
            list of property ID
        stop_if_no_area : bool; default=True
            prevents crashing if there are no elements
        sum_bar_area : bool; default=True
            sum the areas for CBAR/CBEAM/CROD/CONROD/CTUBE elements
            True : get the area of the model by property id (e.g., A=A_pid*Nelements)
            False : only get the cross sectional properties (e.g., A=A_pid)

        Returns
        -------
        pids_to_area : dict[int pid] : float area
            the pid to area dictionary

        TODO: What about CONRODs?
            #'PBRSECT', 'PBCOMP', 'PBMSECT', 'PBEAM3', 'PBEND', 'PIHEX', 'PCOMPS',

        """
        pids_to_area = get_area_breakdown(
            self,
            property_ids=property_ids,
            stop_if_no_area=stop_if_no_area,
            sum_bar_area=sum_bar_area)
        return pids_to_area

    def get_volume_breakdown(self,
                             property_ids=None,
                             stop_if_no_volume=True,
                             ) -> dict[int, float]:
        """
        Gets a breakdown of the volume by property region.

        Parameters
        ----------
        property_ids : list[int] / int
            list of property ID
        stop_if_no_volume : bool; default=True
            prevents crashing if there are no elements

        TODO: What about CONRODs?
        #'PBRSECT',
        #'PBCOMP',
        #'PBMSECT',
        #'PBEAM3',
        #'PBEND',
        #'PIHEX',

        """
        pids_to_volume = get_volume_breakdown(
            self,
            property_ids=property_ids,
            stop_if_no_volume=stop_if_no_volume)
        return pids_to_volume

    def get_element_property_ids(self):
        NO_PROPERTY = {
            'CONROD', 'CONM1', 'CONM2',
        }
        element_id = []
        property_id = []
        for card in self.elements:
            if card.n == 0: # or card.type in NO_MASS:
                continue
            element_id.append(card.element_id)
            if card.type in NO_PROPERTY:
                nelement = len(card.element_id)
                property_id.append(np.full(nelement, 0, dtype='int32'))
                continue
            property_id.append(card.property_id)
        element_ids = np.hstack(element_id)
        property_ids = np.hstack(property_id)
        return element_ids, property_ids

    def get_mass_breakdown_by_property_id_by_material_id(self):
        """TODO: not done"""
        #property_id = []
        #material_id = []
        #mass = []
        NO_DETAILED_MASS = {
            'CROD', 'CBAR', 'CBEAM',
            'CTETRA', 'CPENTA', 'CHEXA', 'CPYRAM',
        }
        NO_MATERIAL_MASS = {
            'CONM1', 'CONM2', 'CMASS1', 'CMASS2', 'CMASS3', 'CMASS4',
        }

        for card in self.elements:
            if card.n == 0 or (card_type := card.type) in NO_MASS:
                continue
            nelement = len(card.element_id)
            if card_type in NO_MATERIAL_MASS:
                property_id = np.zeros(nelement, dtype='float64')
                material_id = np.zeros(nelement, dtype='float64')
                mass = card.mass()
            elif card_type in NO_DETAILED_MASS:  # solids
                property_id = card.property_id
                material_id = card.mass_material_id()
                mass = card.mass()
            else:
                # shells
                element_id, property_id, material_id = card.mass_material_id()

            #card.get_detailed_mass()
            #if card.type in
        x = 1
        #raise RuntimeError('get_mass_breakdown_by_property_id_by_material_id')
        return x

    def nonstructural_mass(self, element_id=None) -> float:
        self.nsmadd.get_nsms_by_nsm_id()

        #self.nsm.get_nsms_by_nsm_id()
        #self.nsm1.get_nsms_by_nsm_id()
        #self.nsml.get_nsms_by_nsm_id()
        #self.nsml1.get_nsms_by_nsm_id()

    def mass(self) -> tuple[np.ndarray, np.ndarray]:
        """
        TODO: limit element ids somehow...by id???
        """
        log = self.log
        element_ids_all = []
        masses = []
        #mass = 0.
        for card in self.element_cards:
            if card.n == 0 or card.type in NO_MASS:
                continue
            element_ids_all.append(card.element_id)
            massi = card.mass()
            massi_sum = massi.sum()
            if np.isnan(massi_sum):
                log.error(f'{card.type} has nan mass; mass={massi}')
                raise RuntimeError(f'{card.type} has nan mass; mass={massi}')
            #log.debug(f'{card.type}; mass={massi_sum}')
            masses.append(massi)

        if len(masses) == 0:
            element_ids = np.array([], dtype='int32')
            mass = np.array([], dtype='float64')
            log.error(f'No elements found; mass={mass}')
            return element_ids, mass

        masses = np.hstack(masses)
        element_ids = np.hstack(element_ids_all)
        return element_ids, masses

    def mass_sum(self, element_id: Optional[np.ndarray]=None) -> float:
        eids, massv = self.mass()
        if element_id is None:
            mass = massv.sum()
        else:
            element_id = np.atleast_1d(element_id).astype('int32')
            isort = np.argsort(eids)
            eids2 = eids[isort]
            mass2 = massv[isort]
            i = np.searchsorted(eids2, element_id)
            assert np.array_equal(eids2[i], element_id)
            mass = mass2[i]
        return mass

    def inertia_sum(self, element_id: Optional[np.ndarray]=None) -> tuple[float, np.ndarray, np.ndarray]:
        """
        mass moment of inertia
        """
        log = self.log
        element_ids_all = []
        inertias = []
        total_mass = 0.
        mass_cg = np.zeros(3, dtype='float64')
        for card in self.element_cards:
            if card.n == 0 or card.type in NO_MASS:
                continue
            card.center_of_mass()
            #if card.type == 'CONM2':
                #continue

            element_ids_all.append(card.element_id)
            #inertiai = card.inertia()
            massi = card.mass()
            massi_sum = massi.sum()
            assert massi.shape == (card.n, ), massi.shape

            if hasattr(card, 'center_of_mass'):
                centroid = card.center_of_mass()
                assert centroid.shape == (card.n, 3), centroid.shape
            else:
                centroid = card.centroid()
                assert centroid.shape == (card.n, 3), centroid.shape

            total_mass += massi_sum
            if np.any(np.isnan(centroid)):
                log.error(f'{card.type} has nan centroid; centroid={centroid}')
                raise RuntimeError(f'{card.type} has nan centroid; centroid={centroid}')
            #mass += massi.sum()
            #log.debug(f'{card.type}; mass={massi_sum}')
            mass_centroid = centroid * massi[:, np.newaxis]
            mass_cg += mass_centroid.sum(axis=0)
            inertias.append((massi, centroid))

        if total_mass == 0.:
            cg = np.full(3, np.nan, dtype='float64')
            inertia = np.full(6, np.nan, dtype='float64')
            log.error('no elements with mass...inertia is nan')
            return total_mass, cg, inertia

        cg = mass_cg / total_mass
        Ixx = 0.
        Iyy = 0.
        Izz = 0.
        Ixy = 0.
        Ixz = 0.
        Iyz = 0.
        for (massi, centroid) in inertias:
            dx = centroid[:, 0] - cg[0]
            dy = centroid[:, 1] - cg[1]
            dz = centroid[:, 2] - cg[2]
            Ixx += (massi * (dy ** 2 + dz ** 2)).sum()
            Iyy += (massi * (dx ** 2 + dz ** 2)).sum()
            Izz += (massi * (dx ** 2 + dy ** 2)).sum()
            Ixy += (massi * (dx * dy)).sum()
            Ixz += (massi * (dx * dz)).sum()
            Iyz += (massi * (dy * dz)).sum()
        inertia = np.array([Ixx, Iyy, Izz, Ixy, Ixz, Iyz], dtype='float64')

        #if fem1.wtmass != 1.0:
            #print('weight = %s' % (mass1 / fem1.wtmass))
        log.debug(f'mass = {total_mass}')
        log.debug(f'cg   = {cg}')
        #print('Ixx=%s, Iyy=%s, Izz=%s \nIxy=%s, Ixz=%s, Iyz=%s' % tuple(inertia1))
        log.debug(f'Ixx={Ixx:g} Iyy={Iyy:g} Izz={Izz:g}\nIxy={Ixy:g} Ixz={Ixz:g} Iyz={Iyz:g}')
        #if element_id is None:
        return total_mass, cg, inertia

    def inertia(self) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
        """
        [mass, cg, mass moment of inertia]
        """
        log = self.log
        element_ids_all = []
        #inertias = []
        #mass = 0.
        mass_cg = np.zeros(3, dtype='float64')
        masses = []
        centroids = []
        for card in self.element_cards:
            if card.n == 0 or card.type in NO_MASS:
                continue

            element_ids_all.append(card.element_id)
            massi = card.mass()
            masses.append(massi)
            centroidi = card.centroid()
            if np.any(np.isnan(centroidi)):
                log.error(f'{card.type} has nan centroid; centroid={centroidi}')
                raise RuntimeError(f'{card.type} has nan centroid; centroid={centroidi}')
            centroids.append(centroidi)

        if len(masses) == 0:
            element_id = np.array([], dtype='int32')
            mass = np.array([], dtype='float64')
            cg = np.zeros((0, 3), dtype='float64')
            inertia = np.zeros((0, 6), dtype='float64')
            log.error('no elements with mass/inertia')
            return element_id, mass, cg, inertia

        element_id = np.hstack(element_ids_all)
        mass = np.hstack(masses)
        centroid = np.vstack(centroids)
        #abs_mass = np.abs(mass).sum()
        neids = len(element_id)
        #if abs_mass == 0.:
            #assert len(element_id) > 0, element_id
            #cg = np.full(3, np.nan, dtype='float64')
            #inertia = np.full(6, np.nan, dtype='float64')
            #log.error('no elements with mass...inertia is nan')
            #return element_id, abs_mass, cg, inertia
        mass_cg = mass[:, None] * centroid
        imass = (mass != 0)
        cg = np.full(centroid.shape, np.nan, dtype=centroid.dtype)
        cg[imass] = mass_cg[imass, :] / mass[imass, np.newaxis]

        #cg = mass_cg.sum(axis=0) / mass.sum()
        #assert len(cg) == 3, cg
        assert cg.shape == (neids, 3), cg.shape

        dxyz = centroid - cg
        dx = dxyz[:, 0]
        dy = dxyz[:, 1]
        dz = dxyz[:, 2]
        Ixx = mass * (dy ** 2 + dz ** 2)
        Iyy = mass * (dx ** 2 + dz ** 2)
        Izz = mass * (dx ** 2 + dy ** 2)
        Ixy = mass * (dx * dy)
        Ixz = mass * (dx * dz)
        Iyz = mass * (dy * dz)
        inertia = np.stack([Ixx, Iyy, Izz, Ixy, Ixz, Iyz], axis=1, out=None)

        nrows = len(mass)
        assert inertia.shape == (nrows, 6), inertia.shape
        return element_id, mass, centroid, inertia

    def length(self) -> float:
        NAN_LENGTHS = {'CBUSH', 'CBUSH1D'}
        length = 0.
        for card in self.element_cards:
            if card.n == 0 or card.type in NO_LENGTH:
                continue
            lengthi = card.length()
            if card.type not in NAN_LENGTHS and np.any(np.isnan(lengthi)):
                self.log.error(f'{card.type} has nan length; length={lengthi}')
                raise RuntimeError(f'{card.type} has nan length; length={lengthi}')
            length += lengthi.sum()
        return length

    def area(self) -> float:
        area = 0.
        for card in self.element_cards:
            if card.n == 0 or card.type in NO_AREA:
                continue
            areai = card.area()
            if np.any(np.isnan(areai)):
                self.log.error(f'{card.type} has nan area; area={areai}')
                raise RuntimeError(f'{card.type} has nan area; area={areai}')
            area += areai.sum()
        return area

    def volume(self) -> float:
        volume = 0.
        for card in self.element_cards:
            if card.n == 0 or card.type in NO_VOLUME:
                continue
            volumei = card.volume()
            if np.any(np.isnan(volumei)):
                self.log.error(f'{card.type} has nan volume; volume={volumei}')
                raise RuntimeError(f'{card.type} has nan volume; volume={volumei}')
            volume += volumei.sum()
        return volume

    def quality(self, cards_to_read: Optional[Set[str]]=None):
        """
        cards_to_read is intended for the gui, which doesn't support
        every card and also masses, which don't have results and so are not
        really elements

        Parameters
        ----------
        cards_to_read : set[str]
            cards that won't be included in element_ids
        """
        # list of elements that have no quality results
        NO_QUALITY = {
            'CELAS1', 'CELAS2', 'CELAS3', 'CELAS4',
            'CDAMP1', 'CDAMP2', 'CDAMP3', 'CDAMP4',
            'CMASS1', 'CMASS2', 'CMASS3', 'CMASS4', 'CONM1', 'CONM2',
            'CROD', 'CONROD', 'CTUBE', 'CBAR', 'CBEAM',
            'CVISC', 'CGAP',
            'CBUSH', 'CBUSH1D', 'CBUSH2D',
            'CSHEAR',
        }
        nelements = 0
        elements2 = [card for card in self.element_cards if card.n > 0]
        for card in elements2:
            if cards_to_read is not None and card.type not in cards_to_read:
                continue
            #if card.type in NOT_ELEMENT:
                #continue
            #if card.n == 0 or card.type in NO_QUALITY:
                #continue
            # all elements have a quality value; it's just None
            #print(f'adding {card.type} to quality; n={card.n}')
            nelements += card.n

        #area = np.full(nelements, np.nan, dtype='float64')
        element_ids = np.full(nelements, -1, dtype='int32')
        taper_ratio = np.full(nelements, np.nan, dtype='float64')
        area_ratio = np.full(nelements, np.nan, dtype='float64')
        max_skew = np.full(nelements, np.nan, dtype='float64')
        aspect_ratio = np.full(nelements, np.nan, dtype='float64')
        min_theta = np.full(nelements, np.nan, dtype='float64')
        max_theta = np.full(nelements, np.nan, dtype='float64')
        dideal_theta = np.full(nelements, np.nan, dtype='float64')
        min_edge_length = np.full(nelements, np.nan, dtype='float64')
        max_warp = np.full(nelements, np.nan, dtype='float64')
        icard_type = np.full(nelements, np.nan, dtype='float64')

        i0 = 0
        icard = 0
        for card in elements2:
            if cards_to_read is not None and card.type not in cards_to_read:
                continue
            n = card.n
            #if n > 0:
                #print(f'adding {card.type} to quality with {n} elements')

            if n == 0:
                continue

            #print(f'card_type={card.type}; icard_type={icard}; [{i0}:{i0+n}]; n={n} eids={card.element_id}')
            element_ids[i0:i0+n] = card.element_id
            if card.type in NO_QUALITY:
                icard += 1
                i0 += n
                continue
            if not hasattr(card, 'quality'):
                icard += 1
                self.log.warning(f'{card.type} has no quality() method')
                i0 += n
                continue
            qualityi = card.quality()
            (areai, taper_ratioi, area_ratioi, max_skewi, aspect_ratioi,
             min_thetai, max_thetai, dideal_thetai, min_edge_lengthi, max_warpi) = qualityi
            taper_ratio[i0:i0+n] = taper_ratioi
            area_ratio[i0:i0+n] = area_ratioi
            max_skew[i0:i0+n] = max_skewi
            aspect_ratio[i0:i0+n] = aspect_ratioi
            min_theta[i0:i0+n] = min_thetai
            max_theta[i0:i0+n] = max_thetai
            dideal_theta[i0:i0+n] = dideal_thetai
            min_edge_length[i0:i0+n] = min_edge_lengthi
            max_warp[i0:i0+n] = max_warpi
            icard_type[i0:i0+n] = icard
            i0 += n
            icard += 1
            #if i0 + n >= 46934:
                #print('*********')

            #if np.any(np.isnan(areai)):
                #self.log.error(f'{card.type} has nan area; area={areai}')
                #raise RuntimeError(f'{card.type} has nan area; area={areai}')
            #area += areai.sum()
        out = (
            element_ids, taper_ratio, area_ratio, max_skew, aspect_ratio,
            min_theta, max_theta, dideal_theta, min_edge_length, max_warp,
            #icard_type,
        )
        return out

    def get_mklist(self) -> np.ndarray:
        """gets the MKLIST vector from MKAERO1/MKAERO2"""
        mklist = []
        mkarray = np.array([])
        for mkaero in self.mkaeros:
            mklist += mkaero.mklist()
        if mklist:
            mkarray = np.hstack([mklist])
            #new_array = [tuple(row) for row in mkarray]
            #unique_pairs = np.lib.arraysetops.unique(new_array, axis=0).tolist()
        return mkarray

    # FLUTTER CARDS

    #def FLFACT(self, sid: int, msg: str='') -> FLFACT:
        #"""gets an FLFACT"""
        #try:
            #return self.flfacts[sid]
        #except KeyError:
            #raise KeyError('sid=%s not found%s.  Allowed FLFACTs=%s'
                           #% (sid, msg, _unique_keys(self.flfacts)))

    def Flutter(self, fid: int, msg: str='') -> FLUTTER:
        """gets a FLUTTER"""
        try:
            return self.flutters[fid]
        except KeyError:
            raise KeyError('fid=%s not found%s.  Allowed FLUTTERs=%s'
                           % (fid, msg, _unique_keys(self.flutters)))

    def Element(self, element_id: int) -> list[Element]:
        elements = []
        for card in self.element_cards:
            if card.n == 0:
                continue
            if element_id in card.element_id:
                elem = card.slice_card_by_element_id(element_id)
                elements.append(elem)
        return elements

    def Property(self, property_id: int) -> list[Property]:
        props = []
        for card in self.property_cards:
            if card.n == 0:
                continue
            if property_id in card.property_id:
                prop = card.slice_card_by_property_id(property_id)
                props.append(prop)
        return props

    def Material(self, material_id: int) -> list[Material]:
        materials = []
        for card in self.material_cards:
            if card.n == 0:
                continue
            if material_id in card.material_id:
                mat = card.slice_card_by_material_id(material_id)
                materials.append(mat)
        return materials

def _unique_keys(mydict: dict[int, Any]) -> np.ndarray:
    """helper method"""
    return np.unique(list(mydict.keys()))

def check_element_ids(elements: list[Any]) -> np.ndarray:
    elements2 = [card for card in elements if card.n > 0]
    list_element_ids2 = [card.element_id for card in elements2]

    if len(list_element_ids2) == 0:
        return np.array([], dtype=elements[0].element_id.dtype)
    element_ids = np.hstack(list_element_ids2)
    uelement_ids = np.unique(element_ids)
    if len(element_ids) != len(uelement_ids):
        #assert len(element_ids) == len(uelement_ids), '{}'
        count_where, count, count2, str_cards = _get_duplicate_cards(
            elements2, element_ids, list_element_ids2)
        msg = ''
        raise RuntimeError(f'Duplicate {msg}Element IDs\n'
                           f'element_ids={count_where}\ncount={count2}\n{str_cards}')
    return uelement_ids


def check_property_ids(properties: list[Any], msg: str='') -> np.ndarray:
    properties2 = [card for card in properties if card.n > 0]
    list_property_ids2 = [card.property_id for card in properties2]
    if len(list_property_ids2) == 0:
        return np.array([], dtype=properties[0].property_id.dtype)

    property_ids = np.hstack(list_property_ids2)
    uproperty_ids = np.unique(property_ids)
    if len(property_ids) != len(uproperty_ids):
        msg2 = ''
        for card in properties2:
            msg2 += card.write()
        print(msg2)
        count_where, count, count2, str_cards = _get_duplicate_cards(
            properties2, property_ids, list_property_ids2)
        raise RuntimeError(f'Duplicate {msg}Property IDs\n'
                           f'property_ids={count_where}\ncount={count2}\n{str_cards}')
        #assert len(property_ids) == len(uproperty_ids)
    return uproperty_ids

def _get_duplicate_cards(properties: list[Any],
                         property_ids: np.ndarray,
                         list_property_ids: list[np.ndarray],
                         ) -> tuple[np.ndarray, np.ndarray, np.ndarray, list[str]]:
    count = np.bincount(property_ids)
    count_where = np.where(count > 1)[0]
    count2 = count[count_where]
    #print(len(properties))
    #print('--------------------------')
    assert len(properties) == len(list_property_ids)
    cards = []
    for card, ids in zip_longest(properties, list_property_ids):
        if card.n == 0:
            continue
        common_ids = np.intersect1d(ids, count_where)
        #log = card.model.log
        #log.info(f'processing {card.type!r}; common_ids={common_ids}')
        if len(common_ids):
            #index = card.index(common_ids)
            #new_card = card.slice_card_by_index(index)
            new_card = card.slice_card_by_id(common_ids)
            if new_card.write() == '':
                raise RuntimeError(f'card={card.type!r} has a write problem')
            #print()
            cards.append(new_card)
        #else:
            #print('nope...')

    str_cards = ''.join((card.write() for card in cards))
    return count_where, count, count2, str_cards
