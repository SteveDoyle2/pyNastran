"""defines a list of cards used by the bdf reader"""

#CASE_CONTROL_CARDS = CASE_CONTROL_INT_CARDS.union(CASE_CONTROL_PLOTTABLE_TYPES)

# case control cards that are also bulk data cards
FLAGGED_CARDS = {
    # of the form 'LOAD = 5', so 'PARAM,POST,-1' doesn't count
    'LOAD', 'SPC', 'FREQ', 'MPC',  # case control + bulk data cards
    'FORCE', 'TRIM', 'DESVAR', 'TSTEP', 'TSTEPNL', 'NSM', 'CLOAD', 'SUPORT1',
    'CSSCHD', 'SDAMPING', 'DLOAD', 'TRIM',
    'SUPORT', # short for SUPORT1
    'ACCEL',  # short for ACCELERATION
    # 'PARAM', # equals sign is problematic
}

# we're going to say that all these cards are are BDF cards, not that they are
# all read
BULK_DATA_CARDS = {
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
    'CPLSTN3', 'CPLSTN6', 'CPLSTN4', 'CPLSTN8',
    'CPLSTS3', 'CPLSTS6', 'CPLSTS4', 'CPLSTS8',
    'CTRAX3', 'CTRAX6', 'CTRIAX', 'CTRIAX6', 'CQUADX', 'CQUADX4', 'CQUADX8',
    'SNORM',

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
    'EQUIV',

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
    'LOAD', 'LSEQ', 'LOADCYN',
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

    'MONPNT1', 'MONPNT2', 'MONPNT3',  ## monitor_points
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
    'RADM', 'VIEW', 'VIEW3D',


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

    # optimization cards
    'DEQATN', 'DTABLE',
    'DCONSTR', 'DESVAR', 'DDVAL', 'DRESP1', 'DRESP2', 'DRESP3',
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
    'SELOAD',

    # super-element sets
    'SESET',  ## se_sets

    'SEBSET', 'SEBSET1',  ## se_bsets
    'SECSET', 'SECSET1',  ## se_csets
    'SEQSET', 'SEQSET1',  ## se_qsets
    #'SEUSET', 'SEUSET1',  ## se_usets
    'SEQSEP',

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

    #: contact
    'BCTPARA',  ## bctpara
    'BCRPARA',  ## bcrpara
    'BCTADD',  ## bctadds
    'BCTSET',  ## bctsets
    'BSURF',  ## bsurf
    'BSURFS',  ## bsurfs
    'BCONP', ## bconp
    'BLSEG', ## blseg
    'BFRIC',


    'TEMPBC',
    'RADMT',
    'RADLST', 'RADMTX', 'RADBND',
    'TEMPP1',
    'TEMPRB',
    'CONVM',

    # cyclic
    'CYJOIN', 'CYSYM',

    # acoustic
    'CHACAB', 'PACABS',

    ## ???
    'ACMODL', 'PANEL', 'SWLDPRM',
    'CWELD', 'PWELD', 'PWSEAM', 'CWSEAM', 'CSEAM', 'PSEAM', 'DVSHAP', 'BNDGRID',
    'MODTRAK', 'DSCONS', 'DVAR', 'DVSET', 'DYNRED',
    'BNDFIX', 'BNDFIX1',
    'AEFORCE', 'UXVEC', 'GUST2',

    # other
    'INCLUDE',  # '='
    'ENDDATA',
}

NASA_CARDS = {
    # axisymmetric
    'CCONEAX', 'AXIC', 'MOMAX', 'PCONEAX', 'POINTAX', 'RINGAX',
    'CTRAPAX', 'PRESAX', 'PTRIAAX', 'CTRIAAX',

    # fluid
    'AXIF', 'CFLUID2', 'CFLUID3', 'CFLUID4', 'RINGFL',
    'BDYLIST', 'FSLIST', 'FLSYM', 'GRIDB',
    'PRESPT',

    # slot
    'GRIDF', 'CAXIF2', 'CAXIF3', 'CAXIF4',
    'GRIDS', 'AXSLOT', 'CSLOT3', 'CSLOT4',
    'SLBDY',

    # free surface fluid
    'MATF', 'CFFREE', 'CFLSTR', 'CFHEX2', 'CFWEDGE',

    # 0d
    'PELAS', 'CELAS2', 'CELAS3', 'CELAS4',
    'PDAMP', 'CDAMP3',
    'PMASS', 'CMASS3', 'CMASS4',
    'CONM2',


    # rod
    'CROD', 'CONROD', 'CBAR', 'BAROR',
    'CNGRNT', # ???

    'CSHEAR',
    'CTRIM6',
    'CTRIARG',

    'CTRIA1', 'CTRIA2', 'CQUAD1', 'CQUAD2',
    'CTRMEM',
    'CQDMEM', 'CQDMEM1', 'CQDMEM2',
    'CTRIA3', 'CQUAD4',

    'CHEXA1', 'CHEXA2', 'CIHEX1', 'CIHEX2', 'CIHEX3',
    'CWEDGE', 'CTETRA',
    'GENEL', 'DMI', 'DMIG', 'DTI', 'PLOTEL',

    # rigid
    'CRIGD1', 'CRIGD2', 'CRIGD3', 'CRIGDR',
    'PBAR', 'PROD', 'PSHEAR', 'PTUBE',
    'PTRIM6',
    'PTRIA1', 'PTRIA2', 'PQUAD1', 'PQUAD2',
    'PSHELL', 'PCOMP', 'PCOMP1', 'PCOMP2',

    'PTRMEM',
    'PQDMEM', 'PQDMEM1', 'PQDMEM2',
    'PIHEX',

    # no idea...
    'CTRSHL', 'PTRSHL',
    'CELBOW', 'PELBOW',
    'CIS2D8', 'PIS2D8',
    'CTRPLT1', 'PTRPLT1',
    'CTORDRG', 'PTORDRG',
    'CPSE2', 'CPSE3', 'PPSE',
    'CTRAPRG', # ???
    'SPCFLD', 'REMFLUX', 'CEMLOOP', 'GEMLOOP', 'BFIELD',

    'GRDSET', 'GRID', 'SPOINT', 'EPOINT', 'OMIT', 'OMIT1', 'PARAM',
    'CORD2R', 'CORD2C', 'CORD2S',

    # materials
    'MAT1', 'MAT2', 'MAT4', 'MAT6', 'MAT8',
    'MATT1', 'MATT4', 'MATS1',
    'TABLEM1', 'TABLEM3',
    'TABLED1',
    'TABLES1',
    'TABDMP1',
    'TABRND1',

    # static loads
    'FORCE', 'FORCE1', 'MOMENT', 'PLOAD2', 'PLOAD3', 'PLOAD4',
    'RFORCE', 'GRAV', 'DEFORM',
    # thermal load
    'TEMP', 'TEMPD', 'TEMPRB', 'TEMPP1',
    'CHBDY', 'PHBDY', 'QBDY1', 'QVOL', 'QVECT',
    'RADLST', 'RADMTX',

    # dynamic loads
    'RLOAD1',
    'TLOAD1', 'TLOAD2',
    'RANDPS', 'RANDT1',
    'DLOAD', 'LOADC',
    'DAREA', 'DAREAS',
    'DELAY', 'DPHASE',
    'FREQ', 'FREQ1', 'FREQ2',
    'PLFACT', # ???

    # aero
    'SET1',
    'AEFACT',
    'MKAERO1', 'MKAERO2',
    'AERO',
    'CAERO1',
    'SPLINE1', 'SPLINE2',
    'FLUTTER', 'FLFACT',
    'STREAML1', 'STREAML2',

    # constraints
    'SPC', 'SPC1', 'SPCD', 'SPCS1', 'SPCADD',
    'MPC', 'MPCADD', 'MPCS',
    'TICS', 'TSTEP',
    'SEQGP',
    'EIGR', 'EIGB', 'EIGC', 'EIGP',
    'TIC', 'TF', 'SUPORT',

    # cyclic
    'CYJOIN',

    # ???
    'BDYC', 'BDYS1', 'GTRAN', 'TRANS',
    'ENDDATA',
}
