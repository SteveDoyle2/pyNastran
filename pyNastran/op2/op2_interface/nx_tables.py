from collections import Counter
from typing import List, Dict
import numpy as np

NX_VERSIONS = ['8.0', '8.5', '9.1', '10.1', '11.0', '11.0.1', '12.0', '12.0.2', '2019.2']

NX_ELEMENTS = {
    0 : 'GRID',
    1 : 'CROD',
    2 : 'CBEAM',
    3 : 'CTUBE',
    4 : 'CSHEAR',
    5 : 'FORCEi/MOMENTi follower',
    6 : 'CTRIA1-old',
    7 : 'PLOAD4 follower',
    8 : 'FLOAD/PLOAD2 follower',
    9 : 'PLOAD/PLOAD2',
    10 : 'CONROD',
    11 : 'CELAS1',
    12 : 'CELAS2',
    13 : 'CELAS3',
    14 : 'CELAS4',
    15 : 'AEROT3',
    16 : 'AEROBEAM',
    17 : 'CTRIA2-old',
    18 : 'CQUAD2-old',
    19 : 'CQUAD1-old',
    20 : 'CDAMP1',
    21 : 'CDAMP2',
    22 : 'CDAMP3',
    23 : 'CDAMP4',
    24 : 'CVISC',
    25 : 'CMASS1',
    26 : 'CMASS2',
    27 : 'CMASS3',
    28 : 'CMASS4',
    29 : 'CONM1',

    30 : 'CONM2',
    31 : 'PLOTEL',
    32 : None,
    33 : 'CQUAD4',
    34 : 'CBAR-34',
    35 : 'CCONE',
    36 : 'CTRIARG-old',
    37 : 'CTRAPRG-old',
    38 : 'CGAP',
    39 : 'CTETRA',

    40 : 'CBUSH1D',
    41 : 'CHEXA1-old',
    42 : 'CHEXA2-old',
    43 : 'CFLUID2',
    44 : 'CFLUID3',
    45 : 'CFLUID4',
    46 : 'FLMASS',
    47 : 'AXIF2',
    48 : 'AXIF3',
    49 : 'AXIF4',

    50 : 'SLOT3',
    51 : 'SLOT4',
    52 : 'CHBDYG/CHBDYP',
    53 : 'CTRIAX6',
    54 : 'TRIM6-old',
    55 : 'CDUM3',
    56 : 'CDUM4',
    57 : 'CDUM5',
    58 : 'CDUM6',
    59 : 'CDUM7',

    60 : 'CDUM8',
    61 : 'CDUM9',
    62 : 'CQDMEM1-old',
    63 : 'CQDMEM2-old',
    64 : 'CQUAD8',
    65 : 'CHEX8-old',
    66 : 'CHEX20-old',
    67 : 'CHEXA',
    68 : 'CPENTA',
    69 : 'CBEND',

    70 : 'CTRIAR',
    71 : '',
    72 : 'AEROQ4',
    73 : 'CFTUBE-old',
    74 : 'CTRIA3',
    75 : 'CTRIA6',
    76 : 'CHEXPR',
    77 : 'CPENPR',
    78 : 'CTETPR',
    79 : 'CPYRAM',
    80 : '',
    81 : '',
    82 : 'CQUADR',
    83 : 'HACAB',
    84 : 'HACBR',
    85 : 'TETRA-nonlinear',
    86 : 'GAP-nonlinear',
    87 : 'TUBE-nonlinear',
    88 : 'TRIA3-nonlinear',
    89 : 'ROD-nonlinear',
    90 : 'QUAD4-nonlinear',
    91 : 'PENTA-nonlinear',
    92 : 'CONROD-nonlinear',
    93 : 'HEXA-nonlinear',
    94 : 'BEAM-nonlinear',
    95 : 'QUAD4-nonlinear',
    96 : 'QUAD8-nonlinear',
    97 : 'TRIA3-nonlinear',
    98 : 'TRIA6-nonlinear',
    99 : '',
    100 : 'CBAR-100',
    101 : 'AABSF',
    102 : 'CBUSH',
    103 : 'CQUADP',
    104 : 'CTRIAP',
    105 : 'CBEAMP',
    106 : 'CDAMP5',
    107 : 'CHBDYE',
    108 : 'CHBDYG',
    109 : 'CHBDYP',

    110 : 'CONV',
    111 : 'CONVM',
    112 : 'QBDY3',
    113 : 'QVECT',
    114 : 'QVOL',
    115 : 'RADBC',
    116 : 'SLIF1D',
    117 : 'CWELDC',  # unlisted in the main table, used in OEF table
    118 : 'CWELDP',  # unlisted in the main table, used in OEF table
    119 : 'CFAST',  # unlisted in the main table, used in OEF table
    120 : '',

    121 : '',
    122 : '',
    123 : '',
    124 : '',
    125 : '',
    126 : '',
    127 : 'CQUAD',  # unlisted in the main table, used in OEF table
    128 : 'CQUADX',  # unlisted in the main table, used in OEF table
    129 : 'RELUC',  # unlisted in the main table, used in OEF table
    130 : 'RES',  # unlisted in the main table, used in OEF table

    131 : 'TETRAE',  # unlisted in the main table, used in OEF table
    132 : 'CTRIA',  # unlisted in the main table, used in OEF table
    133 : 'CTRIAX',  # unlisted in the main table, used in OEF table
    134 : 'LINEOB',  # unlisted in the main table, used in OEF table
    135 : 'LINXOB',  # unlisted in the main table, used in OEF table
    136 : 'QUADOB',  # unlisted in the main table, used in OEF table
    137 : 'TRIAOB',  # unlisted in the main table, used in OEF table
    138 : 'LINX',  # unlisted in the main table, used in OEF table
    139 : 'CQUAD4FD',
    140 : 'CHEXA8FD',

    141 : 'CHEXAP',
    142 : 'CPENTAP',
    143 : 'CTETRAP',
    144 : 'CQUAD144',
    145 : 'VUHEXA',
    146 : 'VUPENTA',
    147 : 'VUTETRA',
    148 : 'HEXAM',  # unlisted in the main table, used in OEF table
    149 : 'PENTAM',  # unlisted in the main table, used in OEF table
    150 : 'TETRAM',  # unlisted in the main table, used in OEF table

    151 : 'QUADM',  # unlisted in the main table, used in OEF table
    152 : 'TRIAM',  # unlisted in the main table, used in OEF table
    153 : 'QUADXM',  # unlisted in the main table, used in OEF table
    154 : 'TRIAXM',  # unlisted in the main table, used in OEF table
    155 : 'QUADPW',  # unlisted in the main table, used in OEF table
    156 : 'TRIAPW',  # unlisted in the main table, used in OEF table
    157 : 'LINEPW',  # unlisted in the main table, used in OEF table
    158 : 'QUADOBM',  # unlisted in the main table, used in OEF table
    159 : 'TRIAOBM',  # unlisted in the main table, used in OEF table
    160 : 'CPENTA6FD',

    161 : 'CTETRA4FD',
    162 : 'CTRIA3FD',
    163 : 'CHEXAFD',
    164 : 'CQUADFD',
    165 : 'CPENTAFD',
    166 : 'CTETRAFD',
    167 : 'CTRIAFD',
    168 : 'CTRIAX3FD',
    169 : 'CTRIAXFD',
    170 : 'CQUADX4FD',

    171 : 'CQUADXFD',
    172 : '',
    173 : '',
    174 : 'LINEOBM',  # unlisted in the main table, used in OEF table
    175 : 'LINXOBM',  # unlisted in the main table, used in OEF table
    176 : 'QUADWGM',  # unlisted in the main table, used in OEF table
    177 : 'TRIAWGM',  # unlisted in the main table, used in OEF table
    178 : 'QUADIB',  # unlisted in the main table, used in OEF table
    179 : 'TRIAIB',  # unlisted in the main table, used in OEF table
    180 : 'LINEIB',  # unlisted in the main table, used in OEF table

    181 : 'LINXIB',  # unlisted in the main table, used in OEF table
    182 : 'QUADIBM',  # unlisted in the main table, used in OEF table
    183 : 'TRIAIBM',  # unlisted in the main table, used in OEF table
    184 : 'LINEIBM',  # unlisted in the main table, used in OEF table
    185 : 'LINXIBM',  # unlisted in the main table, used in OEF table
    186 : 'QUADPWM',  # unlisted in the main table, used in OEF table
    187 : 'TRIAPWM',  # unlisted in the main table, used in OEF table
    188 : 'LINEPWM',  # unlisted in the main table, used in OEF table
    189 : 'VUQUAD',
    190 : 'VUTRIA',

    191 : 'VUBEAM',
    192 : 'CVINT',
    193 : 'QUADFR',  # unlisted in the main table, used in OEF table
    194 : 'TRIAFR',  # unlisted in the main table, used in OEF table
    195 : 'LINEFR',  # unlisted in the main table, used in OEF table
    196 : 'LINXFR',  # unlisted in the main table, used in OEF table
    197 : 'SFINT',  # TODO: GMINTS-OEF??
    198 : 'CNVPEL',
    199 : 'VUHBDY',
    200 : 'CWELD',

    201 : 'CQUAD4FD',
    202 : 'CHEXA8FD',
    203 : 'SLIF1D',
    204 : 'CPENTA6FD',
    205 : 'CTETRA4FD',
    206 : 'CTRIA3FD',
    207 : 'CHEXAFD',
    208 : 'CQUADFD',
    209 : 'CPENTAFD',
    210 : 'CTETRAFD',

    211 : 'CTRIAFD',
    212 : 'CTRIAX3FD',
    213 : 'CTRIAXFD',
    214 : 'CQUADX4FD',
    215 : 'CQUADXFD',
    216 : 'CTETRA4FD',
    217 : 'CTRIA3FD',
    218 : 'CHEXAFD',
    219 : 'CQUADFD',
    220 : 'CPENTAFD',

    221 : 'CTETRAFD',
    222 : 'CTRIAX3FD',
    223 : 'CQUADXFD',
    224 : 'CELAS1',
    225 : 'CELAS3',
    226 : 'CBUSH',
    227 : 'CTRIAR',
    228 : 'CQUADR',
    229 : '',
    230 : '',

    231 : '',
    232 : 'CQUADR-composite',
    233 : 'CTRIAR-composite',
    234 : '',
    235 : '',
    236 : '',
    237 : '',
    238 : '',
    239 : '',

    # per nx 8.5 QRG
    240 : 'CTRIA6',
    241 : 'CQUAD8',
    242 : 'CTRAX3',
    243 : 'CQUADX4',
    244 : 'CTRAX6',
    245 : 'CQUADX8',
    246 : 'CTRAX3',
    247 : 'CQUADX4',
    248 : 'CTRAX6',
    249 : 'CQUADX8',
    255 : 'CPYRAM',
    256 : 'CPYRAM',
    257 : 'CPYRAMFD',
    258 : 'CPYRAMFD',
    259 : 'CTRAX3FD',

    261 : 'CTRAX3FD',
    262 : 'CQUADX4FD',
    263 : 'CTRAX6FD',
    264 : 'QUADX3FD',

    266 : 'QUADX8FD',

    269 : 'CHEXAL',
    270 : 'CPENTAL',
    271 : 'CPLSTN3',
    272 : 'CPLSTN4',
    273 : 'CPLSTN6',
    274 : 'CPLSTN3',
    275 : 'CPLSTS3',
    276 : 'CPLSTS4',
    277 : 'CPLSTS6',
    278 : 'CPLSTS8',
    # 279
    280 : 'CBEAR',
    281 : 'CPLSTN3',
    282 : 'CPLSTN4',
    283 : 'CPLSTN6',
    284 : 'CPLSTN8',
    285 : 'CPLSTS3',
    287 : 'CPLSTS6',
    288 : 'CPLSTS8',
    289 : 'CPLSTN4',
    290 : 'CPLSTS4',
    291 : 'CPLSTN3',
    292 : 'CPLSTN4',
    293 : 'CPLSTN6',
    294 : 'CPLSTS8',
    295 : 'CPLSTS3',
    296 : 'CPLSTS4',
    297 : 'CPLSTS6',
    298 : 'CPLSTS8',

    300 : 'HEXA',
    301 : 'PENTA',
    302 : 'TETRA',
    303 : 'PYRAM',

    304 : 'HEXAL',
    305 : 'PENTAL',
    306 : 'HEXALN',
    307 : 'CPENTALN',

    312 : 'TRAX3',
    313 : 'QUADX4',
    314 : 'TRAX6',
    315 : 'QUADX8',

    316 : 'PLSTN3',
    317 : 'PLSTN4',
    318 : 'PLSTN6',
    319 : 'PLSTN8',

    320 : 'PLSTS3',
    321 : 'PLSTS4',
    322 : 'PLSTS6',
    323 : 'PLSTS8',

    328 : 'GPLSTN3',
    329 : 'GPLSTN4',
    330 : 'GPLSTN6',
    331 : 'GPLSTN8',

    337 : 'CHOCK3',
    338 : 'CHOCK4',
    339 : 'CHOCK6',
    340 : 'CHOCK8',
    # SOL 401
    341 : 'CTRIA3 SOL 401',
    342 : 'CQUAD4 SOL 401',
    343 : 'CTRIA6 SOL 401',
    344 : 'CQUAD8 SOL 401',
    345 : 'CTRIAR SOL 401',
    346 : 'CQUADR SOL 401',
    347 : 'CBAR SOL 401',
    348 : 'CBEAM SOL 401',
    349 : 'CBUSH1D SOL 401',
    350 : 'CELAS1 SOL 401',
    351 : 'CELAS2 SOL 401',
    352 : 'CBUSH SOL 401',

    #ELTYPE=267 - Composite HEXA element (CHEXAL)
    #ELTYPE=268 - Composite PENTA element (CPENTAL)
    355 : 'Composite triangular shell element (CTRIA6); SOL 402?',
    356 : 'Composite quadrilateral shell element (CQUAD8); SOL 402?',
    357 : 'Composite triangular shell element (CTRIAR); SOL 402?',
    358 : 'Composite quadrilateral shell element (CQUADR); SOL 402?',

    # SOL 402
    363 : 'CROD SOL 402',

    400 : 'CELAS1 - Basic System',
    401 : 'CELAS2 - Basic System',
    402 : 'CDAMP1 - Basic System',
    403 : 'CDAMP2 - Basic System',
    404 : 'CBUSH1D - Basic System',
    405 : 'CBUSH - Basic System',
    406 : 'CVISC - Basic System',
}

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
    b'GEOM1EXA', b'GEOM2EXA', b'GEOM4EXA',
    b'VIEWTB',
    b'R1TABRG',
    b'ERRORN',
    b'BGPDTVU', # basic grid point defintion table for a superelement and related to geometry with view-grids added

    # ???
    b'GEOM1ATV', b'GEOM2ATV', b'EPTATV', b'PTMIC', b'ATVMAP',
]  # type: List[bytes]

NX_MATRIX_TABLES = [
    b'ATV',
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
]  # type: List[bytes]

NX_EXTRA_TABLES = [
    # geometry, but buggy in the geometry block...
    b'ICASE',

    # geometry
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
]  # type: List[bytes]

NX_RESULT_TABLES = [
    # ???
    b'OSTR1THC',
    b'OSTR1PLC',
    b'OSTR1CRC',
    b'OSTR1PL',
    b'OSTR1CR',
    b'OEFIIP',
    b'OESRIP',
    b'OESRIS',
    b'ODELBGPD',
    b'ODAMGCZT',
    b'ODAMGCZR',
    b'ODAMGCZD',
    b'XCASECC',
    b'RST',

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
    b'OESCRM2C',
    b'OSTCRM2C',

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

    b'OGK1', # gasket
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
    67 : 'OQG - Glue force results ???',
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
    504 : b'OEFCRM1 - Cumulative Root Mean Square output',
    604 : 'OEFPSD2 - ???',
    605 : 'OSTPSD2C - ???',
    804 : 'OEFRMS1 - ???',
    805 : 'OESXRMS1 - element RMS stresses for random analysis that includes von Mises stress output.',
    904 : 'OEFNO1 - Cumulative Root Mean Square output',
    905 : 'OESXNO1C - Cumulative Root Mean Square output',
}  # type: Dict[int, str]

NX_OEF_REAL_MAPPER = {
    1: 3,    # CROD
    2: 1 + (10 - 1) * 11,  # CBEAM
    3: 3,    # CTUBE
    4: 17,   # CSHEAR
    10: 3,    # CONROD
    11: 2,    # CELAS1
    12: 2,    # CELAS2
    13: 2,    # CELAS3
    14: 2,    # CELAS4

    20: 2,    # CDAMP1
    21: 2,    # CDAMP2
    22: 2,    # CDAMP3
    23: 2,    # CDAMP4
    24: 3,    # CVISC

    33: 9,    # CQUAD4
    34: 9,    # CBAR
    35: 7,    # CCONEAX
    38: 9,    # CGAP
    40: 8,    # CBUSH1D ???
    64: 2 + (11 - 2) * 5,  # CQUAD8
    69: 1 + (8 - 1) * 2,  # CBEND
    70: 2 + (11 - 2) * 4,  # CTRIAR
    74: 9,    # CTRIA3
    75: 2 + (11 - 2) * 4,  # CTRIA6

    #76:  16,   # Acoustic Velocity/Pressure CHEXA ???
    76: None,  # dummy so it doesnt go into the real results
    77: 10,   # Acoustic Velocity/Pressure CPENTA
    78: 10,   # Acoustic Velocity/Pressure CTETRA

    82: 2 + (11 - 2) * 5,  # CQUADR
    95: 9,    # composite CQUAD4 ???
    96: 9,    # composite CQUAD8 ???
    97: 9,    # composite CTRIA3 ???
    98: 9,    # composite CTRIA6 ???
    100: 8,    # CBAR-100
    102: 7,    # CBUSH
    144: 2 + (11 - 2) * 5,  # bilinear CQUAD4
    189: 6 + (19 - 6) * 4,  # VUQUAD
    190: 6 + (19 - 6) * 3,  # VUTRIA
    191: 4 + (12 - 4) * 2,  # VUBEAM
    200: 9,    # CWELD
    232: 9,    # composite CQUADR ???
    233: 9,    # composite TRIAR ???
    235: 9,    # punch CQUADR...num_wide in DMAP is wrong...left out first entry...
    236: 8,    # punch CTRIAR
}
NX_OEF_IMAG_MAPPER = {
    1: 5,    # CROD
    2: 1 + (17 - 1) * 11,  # CBEAM
    3: 5,     # CTUBE
    4: 33,    # CSHEAR
    10: 5,    # CONROD

    11: 3,    # CELAS1
    12: 3,    # CELAS2
    13: 3,    # CELAS3
    14: 3,    # CELAS4

    20: 3,    # CDAMP1
    21: 3,    # CDAMP2
    22: 3,    # CDAMP3
    23: 3,    # CDAMP4
    24: 5,    # CVISC

    33: 17,   # CQUAD4-centroid
    34: 17,   # CBAR-34
    35: 7,    # CCONEAX # needed to not crash the code...
    38: 9,    # CGAP
    40: 8,    # CBUSH1D ???

    64: 2 + (19 - 2) * 5,  # CQUAD8
    69: 1 + (14 - 1) * 2,  # CBEND
    70: 2 + (19 - 2) * 4,  # CTRIAR
    74: 17,   # CTRIA3
    75: 2 + (19 - 2) * 4,  # CTRIA6

    76: 16,   # Acoustic Velocity/Pressure CHEXA_PR
    77: 16,   # Acoustic Velocity/Pressure CPENTA_PR
    78: 16,   # Acoustic Velocity/Pressure CTETRA_PR

    82: 2 + (19 - 2) * 5,  # CQUADR
    95: 9,    # composite CQUAD4 ???
    96: 9,    # composite CQUAD8 ???
    97: 9,    # composite CTRIA3 ???
    98: 9,    # composite CTRIA6 ???
    100: 14,   # BARS
    102: 13,   # CBUSH

    144: 2 + (19 - 2) * 5,  # CQUAD4-bilinear
    189: 6 + (31 - 6) * 4,  # VUQUAD
    190: 6 + (31 - 6) * 3,  # VUTRIA
    191: 4 + (18 - 4) * 2,  # VUBEAM
    200: 17,   # CWELD
    232: 9,    # composite CQUADR ???
    233: 9,    # composite TRIAR ???
    235: 17,   # punch CQUADR...num_wide in DMAP is wrong...left out first entry...
    236: 16,   # punch CTRIAR
}
