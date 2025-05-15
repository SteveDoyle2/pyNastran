acoustic = {
    # loads
    'ACLOAD', 'ACSRCE',

    # element - property
    'CHACBR', 'PACBAR',
    'CHACAB', 'PACABS',
    'CACINF3', 'CACINF4', 'PACINF',
    'CSLOT3', 'CSLOT4', 'AXSLOT', 'GRIDS',
    'CAABSF', 'PAABSF',
}

aero_geom = {
    # aero cards
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
}
static_aero = {
    'AEROS',  ## aeros
    'TRIM', 'TRIM2',  ## trims
    'CSSCHD',  ## csschds
    'DIVERG',  ## divergs
}
dynamic_aero = {
    'AERO',  ## aero
    'GUST',  ## gusts
    'FLUTTER',   ## flutters
    'FLFACT',   ## flfacts
    'MKAERO1', 'MKAERO2',  ## mkaeros
}

basic = {
    '/',
    'ECHOON', 'ECHOOFF',
    'PARAM',

    'GRID', 'GRDSET', 'SPOINT',

    ## coords
    'CORD1R', 'CORD1C', 'CORD1S',
    'CORD2R', 'CORD2C', 'CORD2S',
    #'GMCORD',
}
old_axisymmetric_structure = {
    # axixsymmetric
    #'CCONEAX', # element - removed
    #'PCONEAX', # property - removed
    #'AXIC', # axic - removed
    #'AXIF', # axif - removed
    #'FORCEAX', # loads - removed
}

axisymmetric_structure = {
    'CTRAX3', 'CTRAX6', 'CTRIAX', 'CTRIAX6', 'CQUADX', 'CQUADX4', 'CQUADX8',
}
axisymmetric_loads = {}
structure = {
    ## masses
    'CONM1', 'CONM2',
    'CMASS1', 'CMASS2', 'CMASS3', 'CMASS4',

    ## elements
    # springs
    'CELAS1', 'CELAS2', 'CELAS3', 'CELAS4', # 'CELAS5',
    # bushings
    'CBUSH', 'CBUSH1D', 'CBUSH2D',
    # dampers
    'CDAMP1', 'CDAMP2', 'CDAMP3', 'CDAMP4', 'CDAMP5',
    # fasteners
    'CFAST', 'CWELD',

    'CBAR', 'CBARAO', 'BAROR',
    'CROD', 'CTUBE', 'CBEAM', 'CBEAM3', 'CONROD', 'CBEND', 'BEAMOR',
    'CTRIA3', 'CTRIA6', 'CTRIAR',
    'CQUAD4', 'CQUAD8', 'CQUADR', 'CQUAD',
    'SNORM',

    'CTETRA', 'CPYRAM', 'CPENTA', 'CHEXA',
    'CSHEAR', 'CVISC', 'CRAC2D', 'CRAC3D',
    'CGAP',
    'GENEL',

    ## properties
    'PMASS',
    'PELAS', 'PGAP', 'PFAST', 'PWELD', 'PLPLANE', 'PPLANE',
    'PBUSH', 'PBUSH1D',
    'PDAMP', 'PDAMP5',
    'PROD', 'PBAR', 'PBARL', 'PBEAM', 'PTUBE', 'PBCOMP', 'PBRSECT', 'PBEND',
    'PBEAML', 'PBMSECT', # not fully supported
    'PBEAM3',  # v1.3

    'PSHELL', 'PCOMP', 'PCOMPG', 'PSHEAR', 'PTRSHL',
    'PSOLID', 'PLSOLID', 'PVISC', 'PRAC2D', 'PRAC3D',
    'PCOMPS', 'PCOMPLS',

    ## pdampt
    'PDAMPT',

    ## pelast
    'PELAST',

    ## pbusht
    'PBUSHT',
}
#structure_loads = {

#}
modes = {
    #: methods
    'EIGB', 'EIGR', 'EIGRL',

    #: cMethods
    'EIGC', 'EIGP',
}
materials = {
    ## materials
    'MAT1', 'MAT2', 'MAT3', 'MAT8', 'MAT9', 'MAT10', 'MAT11', 'MAT3D',
    'MATG', 'MATHE', 'MATHP', 'MATEV',

    # 'MATHE'
    #'EQUIV', # testing only, should never be activated...
    ## thermal_materials
    'MAT4', 'MAT5',
}
nonlinear_materials = {
    ## creep_materials
    'CREEP',

    ## Material dependence - MATT1/MATT2/etc.
    'MATT1', 'MATT2', 'MATT3', 'MATT4', 'MATT5', 'MATT8', 'MATT9', 'MATT11',
    'MATS1', #'MATS3', 'MATS8',

    ## tables
    'TABLEM1', 'TABLEM2', 'TABLEM3', 'TABLEM4',  # material tables - temperature

    # nonlinear elastic temperature dependent materials (e.g. creep)
    # sees TABLES1
    'TABLEST',
    # material tables - stress (MATS1, CREEP, MATHP)
    'TABLES1',
}

dynamic_loads = {
    ## dloads
    'DLOAD',

    ## dload_entries
    'ACSRCE', 'TLOAD1', 'TLOAD2', 'RLOAD1', 'RLOAD2',
    'QVECT',
    'RANDPS', 'RANDT1', # random

    # ---- dynamic cards ---- #
    'DAREA',  ## dareas
    'DPHASE',  ## dphases
    'TF',  ## transfer_functions
    'TIC', ## initial conditions - sid (set ID)
}
transient = {
    'DELAY',  ## delays
    'NLPARM',  ## nlparms
    'NLPCI',  ## nlpcis
    'TSTEP',  ## tsteps
    'TSTEPNL', 'TSTEP1',  ## tstepnls
}
frequency = {
    ## frequencies
    'FREQ', 'FREQ1', 'FREQ2', 'FREQ3', 'FREQ4', 'FREQ5',
}

static_loads = {
    ## loads
    'LOAD', 'CLOAD', 'LSEQ', 'LOADCYN', 'LOADCYH',
    'SLOAD',
    'FORCE', 'FORCE1', 'FORCE2',
    'MOMENT', 'MOMENT1', 'MOMENT2',
    'GRAV', 'ACCEL', 'ACCEL1',
    'PLOAD', 'PLOAD1', 'PLOAD2', 'PLOAD4',
    'PLOADX1', 'RFORCE', 'RFORCE1',
    'SPCD', 'DEFORM',

    # msgmesh
    #'GMLOAD',

    # axisymmetric
    'PRESAX',

    #thermal
    'QVOL',
}
static_thermal_loads = {
    # temperature cards
    'TEMP', 'TEMPD', 'TEMPB3', 'TEMPAX',
}
dynamic_thermal_loads = {
    'QBDY1', 'QBDY2', 'QBDY3', 'QHBDY',
    'CHBDYE', 'CHBDYG', 'CHBDYP',
    'PCONV', 'PCONVM', 'PHBDY',
    'RADBC', 'CONV',
    'RADM', 'VIEW', 'VIEW3D',  # TODO: not validated
    'RADCAV', ## radcavs
}
optimization = {
    # optimization cards
    'DEQATN', 'DTABLE',
    'DCONSTR', 'DESVAR', 'TOPVAR', 'DDVAL', 'DRESP1', 'DRESP2', 'DRESP3',
    'DVCREL1', 'DVCREL2',
    'DVPREL1', 'DVPREL2',
    'DVMREL1', 'DVMREL2',
    'DOPTPRM', 'DLINK', 'DCONADD', 'DVGRID',
    'DSCREEN',
}

superelements = {
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
}

sol_101 = basic + structure + materials + static_loads + static_thermal_loads # statics
sol_103 = basic + structure + materials + modes # modes
sol_105 = basic + structure + materials + static_loads # buckling
sol_144 = basic + structure + materials + aero_geom + static_aero + static_loads # static aero
sol_145 = basic + structure + materials + aero_geom + dynamic_aero + modes # flutter
sol_146 = basic + structure + materials + aero_geom + dynamic_aero + modes + frequency + static_loads + dynamic_loads  # gust
sol_153 = basic + structure + materials # + static_loads
sol_200 = optimization + sol_144 + sol_145 + sol_146

cards_to_read = [
    ## nodes
    'EPOINT', 'SEQGP', 'GRIDB',

    # points
    'POINT',
    #'GRIDG'

    ## ringfl
    'RINGFL',
    ## ringaxs
    'RINGAX', 'POINTAX',

    ## nsms
    'NSM', 'NSM1', 'NSML', 'NSML1',

    ## nsmadds
    'NSMADD',

    #'CTRSHL',

    'CPLSTN3', 'CPLSTN4', 'CPLSTN6', 'CPLSTN8', # plate strain
    'CPLSTS3', 'CPLSTS4', 'CPLSTS6', 'CPLSTS8', # plate stress

    # acoustic
    'CHACAB', 'CAABSF', 'CHACBR',
    'PACABS', 'PAABSF', 'PACBAR', 'ACMODL',

    ## rigid_elements
    'RBAR', 'RBAR1', 'RBE1', 'RBE2', 'RBE3', 'RROD', 'RSPLINE', 'RSSCON',

    ## plotels
    'PLOTEL',

    ## nxstrats
    'NXSTRAT',

    ## spcs
    'SPC', 'SPCADD', 'SPC1', 'SPCAX', 'SPCOFF', 'SPCOFF1',
    'GMSPC',

    ## mpcs
    'MPC', 'MPCADD',

    ## suport/suport1/se_suport
    'SUPORT', 'SUPORT1', 'SESUP',

    'ROTORG', 'ROTORD', ## rotors

    # direct matrix input cards
    'DMIG', 'DMIJ', 'DMIJI', 'DMIK', 'DMI', 'DTI',
    'DMIAX',

    # sets
    'SET1', 'SET3',  ## sets
    'ASET', 'ASET1',  ## asets
    'OMIT', 'OMIT1',  ## omits
    'BSET', 'BSET1',  ## bsets
    'CSET', 'CSET1',  ## csets
    'QSET', 'QSET1',  ## qsets
    'USET', 'USET1',  ## usets

    'RADSET',  # radset

    #------------------------------------------------------------------
    ## parametric
    'PSET', 'PVAL', 'GMCURV', 'GMSURF', 'FEEDGE', 'FEFACE',

    #------------------------------------------------------------------
    ## tables
    'TABLED1', 'TABLED2', 'TABLED3', 'TABLED4',  # dynamic tables - freq/time loads

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
