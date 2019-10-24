from collections import Counter
import numpy as np

NX_GEOM_TABLES = [
    b'CASECC',
    b'PVT', b'PVT0', b'PVTS',
    #b'GPLS',
    b'LAMA', b'CLAMA',
    b'OGPWG',
    b'EDT', b'EDTS',
    b'CONTACT', b'CONTACTS', # surface contact definition
    b'GEOM1', b'GEOM2', b'GEOM3', b'GEOM4', b'EPT', b'MPT', b'DYNAMIC', b'DIT', b'EDOM',
    b'GEOM1S', b'GEOM2S', b'GEOM3S', b'GEOM4S', b'EPTS', b'MPTS', b'DYNAMICS',
    b'GEOM1VU', b'GEOM2VU',
    b'GEOM1N',
    b'VIEWTB',
    b'R1TABRG',
    b'ERRORN',
]

NX_MATRIX_TABLES = [
    b'XSOP2DIR',
    b'RADEFMP', # Modal Effective Inertia Matrix - Modal Matrix (per Vibrata)

    # hasn't been validated
    #b'RAFGEN', # Load Set Modal Forces  - Modal generalized force vectors  (per Vibrata)
    #b'RADAMPZ',
    #b'RADAMPG',

    b'EFMFSMS', b'EFMASSS', b'RBMASSS', b'EFMFACS', b'MPFACS', b'MEFMASS', b'MEFWTS',

    # hasn't been validated
    #b'K4HH', b'KELMP', b'MELMP',

    # not-MATPOOL
    # hasn't been validated
    #b'DELTAK', b'DELTAM', b'RBM0', b'DELTAM0',

    # MATPOOL
    # hasn't been validated
    #b'MRGGT', b'UEXPT',

    # MATRIX/MATPOOL - testing-remove this
    # hasn't been validated
    #b'PATRN', b'IDENT', b'RANDM', b'CMPLX',
    #b'MPATRN', b'MIDENT', b'MRANDM', b'MCMPLX',

    b'MATPOOL',
    ##b'KELM',

    # hasn't been validated
    b'MELM', # b'BELM',

    b'BHH', b'KHH',
    b'DSCM2',
]

NX_EXTRA_TABLES = [
    # geometry, but buggy in the geometry block...
    b'ICASE',

    # geometry
    b'BGPDTVU', # basic grid point defintion table for a superelement and related to geometry with view-grids added
    b'DESCYC',
    b'DSCMCOL', # design sensitivity parameters
    b'DBCOPT',  # design optimization history for post-processing
    #--------------------------------------------------------------------------

    # RADx...
    b'RADCONS', b'RADEFFM', b'RADEATC',

    # stress
    b'OES1', b'OES1X', b'OES1X1', b'OES1C',
    b'OES2',
    b'OESNLXR', b'OESNLBR', b'OESNLXD',  # nonlinear
    b'OESNLXR2', b'OESNLBR2',
    b'OESCP', # ???
    b'OESRT', # ???
    b'OESXRMS1', # random with RMS von mises stress

    # strain
    b'OSTR1', b'OSTR1X', b'OSTR1C',
    b'OSTR2',
    b'OESTRCP', # ???

    # contact
    b'OSPDSI1', b'OSPDSI2', # intial separation distance
    b'OSPDS1', b'OSPDS2',   # final separation distance
    b'OBC1', b'OBC2',       # contact pressures and tractions at grid points

    # glue
    b'OBG1', # glue normal and tangential tractions at grid points in basic coordinate system

    # RMAXMIN - Defines parameters to output the minimum, maximum, absolute
    #           value maximum, average, and RMS value of stress, force, and
    #           displacement results for SOLs 101, 109, and 112.
    b'OUGV1MX', # max displacement?
    b'OES1MX',  # max stress?
    b'OEF1MX',  # max force?
    b'OSMPF2M',
    b'OFMPF2M',
    b'OPMPF2M',
    b'OLMPF2M',
    b'OEKE1',
]

NX_RESULT_TABLES = [
    # displacements, velocity, acceleration
    # BOUGV1 - G-set results (displacement, velocity, acceleration, eigenvector)
    #          in the global (CD) frame
    # OUGVi  - Displacements in the global (CD) coordinate system
    b'BOUGV1',  # G-set results boundary in the basic (cid=0) frame
    b'OUGV1', b'OUGV2',
    b'OUGATO2', b'OUGCRM2', b'OUGPSD2',
    b'OUGNO1', b'OUGRMS1',

    # eigenvectors
    b'OPHIG', # Eigenvectors in the basic (cid=0) coordinate system

    # eigenvectors
    b'BOPHIG',  # basic (cid=0) frame

    # temperature
    b'TOUGV1',

    #------------------------
    # solution set
    # OUXYi - ??? set in ??? frame
    #------------------------

    # spc forces
    # OQGx   - SPC forces in the G-set
    #        - can also be MPC forces, but generally not anymore for NX Nastran
    # OQGCFx - ???
    # OQGGFx - ???
    b'OQG1', b'OQG2',
    b'OQGCF1', b'OQGCF2', # ???
    b'OQGGF1', b'OQGGF2', # ???

    # mpc forces
    # OQMGx - MPC forces in the G-set
    b'OQMG1', b'OQMG2',

    # load vector
    # OPGi  - G-set load vectors in the global (CD) frame
    # OPNLi - Nonlinear loads in for the h-set or d-set.
    b'OPG1', b'OPG2',
    b'OPNL1', b'OPNL2',
    b'OPGNO1',
    b'OPGRMS1',

    # Grid point stresses
    b'OGS1',

    # strain energy
    b'ONRGY1', b'ONRGY2', b'ONRGY',

    # failure indicies
    b'OEFIT',

    #-----------------------
    # OESVM1  - OES Table of           element stresses
    # OESVM1C - OES Table of composite element stresses
    #           for frequency response analysis that includes von Mises stress
    #           output in SORT1 format.
    b'OESVM1', b'OESVM1C',
    b'OESVM2',
    b'OSTRVM1', b'OSTRVM1C',
    b'OSTRVM2',

    b'OES2C', b'OSTR2C',

    # hasn't been validated
    b'OESPSD2C', b'OSTPSD2C',
    b'OSTRRMS1', b'OSTRMS1C',
    b'OSTRNO1', b'OSTNO1C',

    # sol 401?
    b'OESNL2',
    b'OSTR1IN',  # OES output table of initial strains at corner grids in the basic coordinate system
    b'OSTR1G',   # Table of total strain at Gauss points in SORT1 format
    b'OSTR1PLG', # Table of plastic strain at Gauss points in SORT1 format
    b'OSTR1THG', # Table of thermal strain at Gauss points in SORT1 format
    b'OSTR1ELG', # Table of elastic strain at Gauss points in SORT1 format
    b'OSTR1TH',  # Table of thermal strain in SORT1 format
    b'OSTR1EL',  # Table of elastic strain in SORT1 format
    b'OSTR1ING', # OES output table of initial strains at corner Gauss points in the basic coordinate system

    b'OES1G',   # Grid point stress or strain table in SORT1 format and interpolated from the centroidal stress table, OES1M.

    #----------------------
    # hasn't been validated...
    b'MDICT', b'BDICT', b'KDICTP', b'MDICTP',

    #----------------------
    # forces
    # OEF1X - Element forces with intermediate (CBAR and CBEAM) station forces
    #         and forces on nonlinear elements
    # OEFx  - Element forces for shells/solids/rods
    b'OEF1X',
    b'OEF1', b'OEF2',

    # heat flux
    b'HOEF1',


    # ---------------------------------------------
    # nx2019.2

    # geometry
    #b'GPDTS',

    # results - supported
    b'OPHSA',   # Displacement output table in SORT1
    b'OUXY1',   # Displacements in SORT1 format for h-set or d-set.
    b'OUXY2',   # Displacements in SORT2 format for h-set or d-set.
    b'OTEMP1',  # Grid point temperature output

    b'LAMAS', # Normal modes eigenvalue summary table for the structural portion of the model
    b'LAMAF', # Normal modes eigenvalue summary table for the fluid portion of the model

    b'BOPHIGF',  # Eigenvectors in the basic coordinate system for the fluid portion of the model.

    # hasn't been validated
    #b'BOPHIGS',  # Eigenvectors in the basic coordinate system for the structural portion of the model.

    # Grid Point Forces - SORT1/SORT2
    b'OGPFB1', b'OGPFB2',

    # ---------------
    # new results
    b'PSDF',    # Power spectral density table.

    # random stress
    b'OESXRM1C', # Table of composite element RMS stresses in SORT1 format for random analysis that includes von Mises stress output.
    b'OESXNO1C',
    b'OESXNO1',

    # ---------------
    # results - unsupported
    b'GPLS',  # needs to be here to prevent a crash
    b'TRMBU',  # Transfomration matrices from undeformed to basic
    b'TRMBD',  # Transformation matrices from deformed to basic
    #b'PVTS', # PVT0?
    b'OEFMXORD', # List of element IDs with maximum frequency and element order diagnostics for FEMAO solution
    b'OBCKL',    # Load factor vs. cumulative arc-length in SORT2 format
    b'ONMD',     # Normalized material density for topology optimization output
    b'OACPERF',  # Performance data that indicates computation time in seconds and memory consumed in GB per frequency per subcase for FEMAO analysis.
    b'OGSTR1',   # Grid point strains of superelement
    b'OSHT1',   # Shell element thickness results (created by SHELLTHK)
    b'OEFIIS',  # Inter-laminar shear failure indices.

    # bolt
    b'OBOLT1',  # Bolt output data block

    # damage
    b'OELAR',    # Element status (active or inactive)
    b'OJINT',    # J-integral for a crack defined by CRAKTP.
    b'ODAMGPFE', # Damage energy for ply failure
    b'ODAMGPFD', # Damage values for ply failure
    b'ODAMGPFS', # Damage status for ply failure
    b'ODAMGPFR', # Crack density for ply failure EUD model from SOL 401. Crack density at corner grids on the middle of plies. The values are unitless

    # contact / glue
    b'OSLIDEG1', # Glue slide distance output
    b'OCONST1',  # Contact status in SORT1 format
    b'OSLIDE1',  # Incremental and total slide output for contact/glue.

    b'OCPSDF',  # Cross-power-spectral-density functions.
    b'OCCORF',  # Cross-correlation functions.
    b'OCPSDFC', # Cross-power spectral density functions for composites.
    b'OCCORFC', # Cross-correlation functions for composites.

    b'OEDE1', # Elemental energy loss.

    # grid point pressure
    b'OPRNO1',   # SORT1 - NO
    b'OPRRMS1',  # SORT1 - RMS
    b'OPRPSD2',  # SORT2 - PSD
    b'OPRATO2',  # SORT2 - AUTO
    b'OPRCRM2',  # SORT2 - CRMS

    # modal contribution
    b'OUGMC1',  # Modal contributions for displacements, velocities, accelerations.
    b'OUGMC2',  # Modal contributions for displacements, velocities, accelerations.
    b'OQGMC1',  # Modal contributions of single point constraint forces - SORT1
    b'OQGMC2',  # Modal contributions of single point constraint forces - SORT2
    b'OESMC1',  # Element stress modal contributions - SORT1
    b'OESMC2',  # Element stress modal contributions - SORT2
    b'OSTRMC1', # Modal contributions of element strains - SORT1
    b'OSTRMC2', # Modal contributions of element strains - SORT2
    b'OEFMC1',  # Modal contributions of element forces - SORT1
    b'OEFMC2',  # Modal contributions of element forces - SORT2

    # modal strain energy
    b'OMSEC1',  # Constant modal strain energy - SORT1
    b'OMSEC2',  # Constant modal strain energy - SORT2
    b'OMECON1', # Constant total modal energies - SORT1
    b'OMECON2', # Constant total modal energies - SORT2
    b'OMEOSC1', # Oscillating total modal energies - SORT1
    b'OMEOSC2', # Oscillating total modal energies - SORT2
    b'OMKEC1',  # Constant modal kinetic energies - SORT1
    b'OMKEC2',  # Constant modal kinetic energies - SORT2
    b'OMKEO2',  # Oscillating modal kinetic energies - SORT2
    b'OMSEO1',  # Oscillating modal strain energies - SORT2
    b'OMSEO2',  # Oscillating modal strain energies - SORT2
    b'OMKEO1',  # Oscillating modal kinetic energies - SORT1

    # radiated power
    b'OERP',     # Equivalent radiated power output.
    b'OERPEL1',  # Element equivalent radiated power (element output)
    b'OERPEL2',  # Element equivalent radiated power output.

    # acoustic
    b'OUGPC1',  # Table of panel contributions - SORT1
    b'OUGPC2',  # Table of panel contributions - SORT2
    b'OUGF1',   # Acoustic pressures at microphone points in SORT1 format
    b'OUGF2',   # Acoustic pressures at microphone points in SORT2 format
    b'OUGGC1',  # Table of grid contributions - SORT1
    b'OUGGC2',  # Table of grid contributions - SORT2
    b'OUGRC1',  # Reciprocal panel contributions - SORT1
    b'OUGRC2',  # Reciprocal panel contributions - SORT2
    b'BOUGF1',  # Acoustic pressures at microphone points in SORT1 format - basic frame

    # acoustic acceleration
    b'OACCQ',    # Acoustic coupling quality
    b'OACINT1',  # Acoustic intensities at microphone points - SORT1
    b'OACINT2',  # Acoustic intensities at microphone points - SORT2
    b'OACVELO1', # Acoustic velocities at microphone points - SORT1
    b'OACVELO2', # Acoustic velocities at microphone points - SORT2
    b'OACPWR2',  # Acoustic power for AML regions and GROUPs of 2D elements - SORT2
    b'OACPWRI2', # Acoustic incident power - SORT2
    b'OACPWRT2', # Transmitted acoustic power for AML regions and GROUPs of 2D elements - SORT2
    b'OACTRLS2', # Acoustic transmission loss - SORT2

    # random acoustic
    b'OAPPSD2', # Acoustic power for the PSD function - SORT2
]

if len(NX_RESULT_TABLES) != len(np.unique(NX_RESULT_TABLES)):  # pragma: no cover
    counter = Counter(NX_RESULT_TABLES)
    _MSG = 'Invalid count:\n'
    for key, cvaluei in counter.items():
        if cvaluei != 1:
            _MSG += '%s = %s\n' % (key, cvaluei)
    raise RuntimeError(_MSG)

NX_RESULT_TABLES += NX_EXTRA_TABLES

NX_TABLE_CONTENT = {
    # nx 8.5
    0: '',
    1: 'OUG - Displacement vector',
    2: 'OPG - Load vector',
    3: 'OQG - SPC/MPC Force vector',
    4: 'OEF - Element force/flux',
    5: 'OES - Element stress/strain',
    6: 'LAMA - Eigenvalue summary',
    7: 'OUG - Eigenvector',
    8: 'Grid Point Singularity Table (obsolete)',
    9: 'OEIGS - Eigenvalue analysis summary',
    10: 'OUG - Velocity vector',
    11: 'OUG - Acceleration vector',
    12: 'OPG - Nonlinear force vector',
    13: 'OGPWG - Grid point weight generator',
    14: 'OUG - Eigenvector (solution set)',
    15: 'OUG - Displacement vector (solution set)',
    16: 'OUG - Velocity vector (solution set)',
    17: 'OUG - Acceleration vector (solution set)',
    18: 'OEE - Element strain energy',
    19: 'OGF - Grid point force balance',
    20: 'OES - Stresses at grid points',
    21: 'OES - Strain/curvature at grid points',
    22: 'OELOF1 - Element internal forces/moments',
    23: 'OELOP1 - Summation of element oriented forces on adjacent elements',
    24: 'OEP - Element pressures',
    25: 'OEF - Composite failure indices',
    26: 'OGS - Grid point stresses (surface)',
    27: 'OGS - Grid point stresses (volume - direct)',
    28: 'OGS - Grid point stresses (volume - princial)',
    29: 'OGS - Element stress discontinuities (surface)',
    30: 'OGS - Element stress discontinuities (volume - direct)',
    31: 'OGS - Element stress discontinuities (volume - princial)',
    32: 'OGS - Grid point stress discontinuities (surface)',
    33: 'OGS - Grid point stress discontinuities (volume - direct)',
    34: 'OGS - Grid point stress discontinuities (volume - princial)',
    35: 'OGS - Grid point stresses (plane strain)',
    36: 'OEE - Element kinetic energy',
    37: 'OEE - Element energy loss',

    38 : 'OMSEC - Constant modal strain energy',
    39 : 'OMSED - Oscillating modal strain energy',
    40 : 'OMKEC - Constant modal kinetic energy',
    41 : 'OMKED - Oscillating modal kinetic energy',
    42 : 'OMECON - Constant total modal energy',
    43 : 'OMEOSC - Oscillating total modal energy',
    44 : 'OUGMC - Displacement/velocity/acceleration modal contributions',
    45 : 'OEFMC - Element force modal contributions',
    46 : 'OESMC - Element stress modal contributions',
    47 : 'OSTRMC - Element strain modal contributions',
    48 : 'OQGMC - SPC force modal contributions',
    49 : 'OUGPC - Panel contributions',
    50 : 'OUGGC - Grid contributions',
    51 : 'OUGRC - Reciprocal panel contributions',
    #
    53 : 'OACVELO - Acoustic velocity',
    54 : 'OACINT - Acoustic intensity',
    55 : 'OACPWR - Acoustic power',
    56 : 'OACPWRI - Acoustic incident power',
    57 : 'OACPWRT - Acoustic transmitted power',
    58 : 'OACTRLS - Acoustic transmission loss',
    #
    61 : 'OGK - Gasket Element Results',
    62 : 'OBC - Contact Pressure and Traction Results',
    63 : 'OQG - Contact Force Results',
    64 : 'OSPDSI - Contact Separation Distance - Initial',
    65 : 'OSPDS - Contact Separation Distance',
    66 : 'OBG - Glue force results (normal and in-plane tractions)',
    #67 : 'OQG - Glue force results',
    68 : 'ELRSCALV - Tosca normalized material properties',
    69 : 'OERP - Element equivalent radiated power (panel output)',
    70 : 'OERPEL - Element equivalent radiated power (element output)',
    71 : 'Reserved for FE-Design',
    72 : 'OTEMP - Grid point temperature output',
    73 : 'JINT - Crack front J-integral output',
    74 : 'SLIDE - Contact/Glue slide output',
    75 : 'CONSTAT - Contact status output',
    76 : 'OERR - Error estimator output',
    77 : 'OPRESS - Grid point pressure output',
    78 : 'STATE - Variables output',
    79 : 'INITSTR - Initial strain output',
    80 : 'OBOLT - Bolt preload output',

    81 : 'OCKGAP1 - Opening gap values for chocking elements',
    82 : 'ODAMGCZD - Damage values for cohesive elements',
    83 : 'ODAMGCZR - Relative displacements for cohesive elements',
    84 : 'ODAMGCZT - Tractions for cohesive elements',
    85 : 'ODAMGPFD - Damage values for PFA',
    86 : 'ODAMGPFE - Damage energy for PFA',
    87 : 'ODAMGPFS - Damage status for PFA',
    88 : 'ODAMGPFR - Damage crack density for PFA / EUD',
    89 : '??? - Composite strength ratios',
    90 : 'TRMBD - Transformation matrices from deformed to basic',
    91 : 'TRMBU - Transfomration matrices from undeformed to basic',
    92 : 'ONMD - Normalized material density for topology optimization output',
    93 : 'OBCKL - SORT2 output for Load Factor versus Cummulative Arc-length from a SOL 401 arc-length solution',
    #
    # nx 2019.2
    #
    804 : 'OEFRMS1 - ???',
    805 : 'OESXRMS1 - element RMS stresses for random analysis that includes von Mises stress output.',
    905 : 'OESXNO1C - Cumulative Root Mean Square output',
}
