from collections import Counter
from typing import List, Dict
import numpy as np

MSC_ELEMENTS = {
    None: '',
    0: 'GRID',
    1: 'CROD',
    2: 'CBEAM',
    3: 'CTUBE',
    4: 'CSHEAR',
    5: 'FORMON12-FORCEi/MOMENTi follower',
    6: 'FORCE',
    7: 'PLOAD4',
    8: 'PLOADX1',
    9: 'PLOAD/PLOAD2',
    10: 'CONROD',
    11: 'CELAS1',
    12: 'CELAS2',
    13: 'CELAS3',
    14: 'CELAS4',
    15: 'AEROT3',
    16: 'AEROBEAM',
    17: 'CTRIAX',
    18: 'CQUADX',
    19: 'CQUAD1-old',
    20: 'CDAMP1',
    21: 'CDAMP2',
    22: 'CDAMP3',
    23: 'CDAMP4',
    24: 'CVISC',
    25: 'CMASS1',
    26: 'CMASS2',
    27: 'CMASS3',
    28: 'MASS4',
    29: 'CONM1',

    30: 'CONM2',
    31: 'PLOTEL',
    32: None,
    33: 'CQUAD4',
    34: 'CBAR-34',
    35: 'CCONEAX',
    36: 'CTRIARG-old',
    37: 'CTRAPRG-old',
    38: 'CGAP',
    39: 'CTETRA',

    40: 'CBUSH1D',
    41: 'CHEXA1-old',
    42: 'CHEXA2-old',
    43: 'CFLUID2',
    44: 'CFLUID3',
    45: 'CFLUID4',
    46: 'FLMASS',
    47: 'AXIF2',
    48: 'AXIF3',
    49: 'AXIF4',

    50: 'SLOT3',
    51: 'SLOT4',
    52: 'CHBDYG/CHBDYP',
    53: 'CTRIAX6',
    54: 'TRIM6-old',
    55: 'CDUM3',
    56: 'CDUM4',
    57: 'CDUM5',
    58: 'CDUM6',
    59: 'CDUM7',

    60: 'CDUM8/CCRAC2D',
    61: 'CDUM9/CCRAC3D',
    62: 'CQDMEM1',
    63: 'IFQUAD',
    64: 'CQUAD8',
    65: 'IFHEX',
    66: 'IFPENT',
    67: 'CHEXA',
    68: 'CPENTA',
    69: 'CBEND',

    70: 'CTRIAR',
    71: None,
    72: 'AEROQ4',
    73: 'IFQDX',
    74: 'CTRIA3',
    75: 'CTRIA6',
    76: 'HEXPR',
    77: 'PENPR',
    78: 'TETPR',
    79: None,
    80: None,
    81: None,
    82: 'QUADR',
    83: 'HACAB',
    84: 'HACBR',
    85: 'TETRANL-nonlinear',
    86: 'GAPNL-nonlinear',
    87: 'TUBENL-nonlinear',
    88: 'TRIA3NL-nonlinear',
    89: 'RODNL-nonlinear',
    90: 'QUAD4NL-nonlinear',
    91: 'PENTANL-nonlinear',
    92: 'CONRODNL-nonlinear',
    93: 'HEXANL-nonlinear',
    94: 'BEAMNL-nonlinear',
    95: 'QUAD4LC-composite',
    96: 'QUAD8LC-composite',
    97: 'TRIA3LC-composite',
    98: 'TRIA6LC-composite',
    99: None,
    100: 'BAR-100',
    101: 'AABSF',
    102: 'CBUSH',
    103: 'QUADP',
    104: 'TRIAP',
    105: 'BEAMP',
    106: 'DAMP5',
    107: 'CHBDYE',
    108: 'CHBDYG',
    109: 'CHBDYP',

    110: 'CONV',
    111: 'CONVM',
    112: 'QBDY3',
    113: 'QVECT',
    114: 'QVOL',
    115: 'RADBC',
    116: 'SLIF1D',
    117: 'WELDC',
    118: 'WELDP',
    119: 'SEAM',
    120: 'GENEL',

    121: 'DMIG',
    122: 'DIEL-old',
    123: 'HEXAE-old',
    124: 'IND-old',
    125: 'LINE-old',
    126: 'FASTP',
    127: 'CQUAD-old',
    128: 'CQUADX-old',
    129: 'RELUC-old',
    130: 'RES-old',

    131: 'CTETRAE-old',
    132: 'CTRIA-old',
    133: 'CTRIAX-old',
    134: 'LINEOB-old',
    135: 'LINXOB-old',
    136: 'QUADOB-old',
    137: 'TRIAOB-old',
    138: 'LINEX-old',
    139: 'QUAD4FD',
    140: 'HEXA8FD',

    141: 'HEXAP',
    142: 'PENTAP',
    143: 'TETRAP',
    144: 'QUAD144',
    145: 'VUHEXA',
    146: 'VUPENTA',
    147: 'VUTETRA',
    148: 'HEXAM-old',
    149: 'PENTAM-old',
    150: 'TETRAM-old',
    151: 'QUAD-old',
    152: 'TRIAM-old',
    153: 'QUADXM-old',
    154: 'TRIAXM-old',
    155: 'RADINT',
    156: 'BUSH2D',
    157: 'LINEPW-old',
    158: 'QUADOBM-old',
    159: 'SEAMP',
    160: 'PENTA6FD',
    161: 'TETRA4FD',
    162: 'TRIA3FD',
    163: 'HEXAFD',
    164: 'QUADFD',
    165: 'PENTAFD',
    166: 'TETRAFD',
    167: 'TRIAFD',
    168: 'TRIAX3FD',
    169: 'TRIAXFD',
    170: 'QUADX4FD',
    171: 'QUADXFD',
    172: 'QUADRNL',
    173: 'TRIARNL',
    174: 'LINEOBM-old',
    175: 'LINXOBM-old',
    176: 'QUADWGM-old',
    177: 'TRIAWGM-old',
    178: 'QUADIB-old',
    179: 'TRIAIB-old',
    180: 'LINEIB-old',
    181: 'LINEXIB-old',
    182: 'QUADIBM-old',
    183: 'TRIAIBM-old',
    184: 'CBEAM3',
    185: 'LINXIBM-old',
    186: 'QUADPWM-old',
    187: 'TRIAPWM-old',
    188: 'LINEPWM-old',
    189: 'VUQUAD',
    190: 'VUTRIA',
    191: 'VUBEAM',
    192: 'CVINT',
    193: 'QUADFR-old',
    194: 'TRIAFR-old',
    195: 'LINEFR-old',
    196: 'LINXFR-old',
    197: 'SFINT',
    198: 'CNVPEL',
    199: 'VUHBDY',
    200: 'WELD',
    201: 'QUAD4FD',
    202: 'HEXA8FD',
    203: 'SLIF1D',
    204: 'PENTA6FD',
    205: 'TETRA4FD',
    206: 'TRIA3FD',
    207: 'HEXAFD',
    208: 'QUADFD',
    209: 'PENTAFD',
    210: 'TETRAFD',
    211: 'TRIAFD',
    212: 'TRIAX3FD',
    213: 'TRIAXFD',
    214: 'QUADX4FD',
    215: 'QUADXFD',
    216: 'TETRA4FD',
    217: 'TRIA3FD',
    218: 'HEXAFD',
    219: 'QUADFD',
    220: 'PENTAFD',
    221: 'TETRAFD',
    222: 'TRIAX3FD',
    223: 'QUADXFD',
    224: 'ELAS1',
    225: 'ELAS3',
    226: 'BUSH',
    227: 'RBAR',
    228: 'RBE1',
    229: 'RBE3',
    230: 'RJOINT',
    231: 'RROD',
    232: 'QUADRLC',
    233: 'TRIARLC',
    234: None,
    235: 'CQUADR',
    236: 'CTRIAR',
    237 : 'CTRIAR',
    238 : 'CBAR',
    239 : 'CEAM',
    240 : 'CBAR',
    241 : 'AXISYM',
}

MSC_GEOM_TABLES = [
    # GEOM2 - Table of Bulk Data entry images related to element connectivity andscalar points
    # GEOM4 - Table of Bulk Data entry images related to constraints, degree-of-freedom membership and rigid element connectivity.
    b'GEOM1', b'GEOM2', b'GEOM3', b'GEOM4',  # regular
    b'GEOM1S', b'GEOM2S', b'GEOM3S', b'GEOM4S', # superelements
    b'GEOM1N', b'GEOM1VU', b'GEOM2VU',
    b'GEOM1OLD', b'GEOM2OLD', b'GEOM4OLD',

    b'EPT', b'EPTS', b'EPTOLD',
    b'EDTS',
    b'MPT', b'MPTS',
    b'AXIC',

    b'DIT', b'DITS',

    b'PVT', b'PVT0', b'CASECC',
    b'EDOM',
    b'DYNAMIC', b'DYNAMICS',

    b'ERRORN',
    b'DESTAB', b'R1TABRG',

    # eigenvalues
    b'BLAMA', b'LAMA', b'CLAMA',  #CLAMA is new

    # grid point weight
    b'OGPWG', b'OGPWGM',

    # other
    b'CONTACT', b'VIEWTB',
    b'KDICT',
    #b'MDICTP' where does this go?

    # aero?
    #b'MONITOR',
    b'CASEXX',
]

MSC_MATRIX_TABLES = [
    #b'TOLD',
    b'SDT', #b'STDISP',
    b'TOLB2', b'ADSPT', #b'MONITOR',
    b'PMRT', b'PFRT', b'PGRT', # b'AEMONPT',
    b'AFRT', b'AGRT',
    b'QHHA',

    b'A', b'SOLVE,', b'UMERGE,', b'AA', b'AAP', b'ADELUF', b'ADELUS', b'ADELX',
    b'ADJG', b'ADJGT', b'ADRDUG', b'AEDBUXV', b'AEDW', b'AEFRC', b'AEIDW',
    b'AEIPRE', b'AEPRE', b'AG', b'AGD', b'AGG', b'AGX', b'AH', b'AJJT', b'AM2',
    b'AM3', b'ANORM', b'APART', b'APIMAT', b'APIMATT', b'APL', b'APPLOD',
    b'APU', b'ARVEC', b'AUG1', b'B', b'SOLVE,', b'B2DD', b'B2GG', b'B2PP',
    b'BAA', b'BACK', b'BANDPV', b'BASVEC', b'BASVEC0', b'BCONXI', b'BCONXT',
    b'BDD', b'BDIAG', b'BFEFE', b'BFHH', b'BHH', b'BHH1', b'BKK', b'BP', b'BPP',
    b'BRDD', b'BXX', b'BUX', b'C', b'CDELB', b'CDELK', b'CDELM', b'CFSAB',
    b'CLAMMAT', b'CLFMAT', b'CMAT', b'CMBXPHG', b'CMSQE', b'CMSTQE', b'CNVTST',
    b'CON', b'CONS1T', b'CONSBL', b'CONTVDIF', b'COORD', b'COORDO', b'CPH1',
    b'CPH2', b'CPHP', b'CPHX', b'CPHL', b'CVAL', b'CVALO', b'CVAL', b'CVALR',
    b'CVALRG', b'CVECT', b'D', b'D1JE', b'D1JK', b'D2JE', b'D2JK', b'DAR',
    b'DBUG', b'DCLDXT', b'DELB1', b'DELBSH', b'DELCE', b'DELDV', b'DELF1',
    b'DELFL', b'DELGM', b'DELGS', b'DELS', b'DELS1', b'DELTGM', b'DELVS',
    b'DELWS', b'DELX', b'DELX1', b'DESVCP', b'DESVEC', b'DESVECP', b'DJX',
    b'DM', b'DPHG', b'DPLDXI', b'DPLDXT', b'DRDUG', b'DRDUGM', b'DSCM',
    b'DSCM2', b'DSCMG', b'DSCMR', b'DSDIV', b'DSEGM', b'DSESM', b'DSTABR',
    b'DSTABU', b'DUGNI', b'DUX', b'DXDXI', b'DXDXIT', b'E', b'EFMASMTT',
    b'EFMMCOL', b'EFMMAT', b'EGK', b'EGM', b'EGTX', b'EGX', b'EMAT', b'EMM',
    b'ENEMAT', b'ENFLODB', b'ENFLODK', b'ENFLODM', b'ENFMOTN', b'ERHM',
    b'EUHM', b'EXCITEFX', b'EXCITF', b'EXCITP', b'F', b'F2J', b'FFAJ', b'FGNL',
    b'FMPF', b'FN', b'FOLMAT', b'FORE', b'FREQMASS', b'FRMDS', b'GC', b'GDGK',
    b'GDKI', b'GDKSKS', b'GEG', b'GLBRSP', b'GLBRSPDS', b'GM', b'GMD', b'GMNE',
    b'GMS', b'GOA', b'GOD', b'GPFMAT', b'GPGK', b'GPKH', b'GPIK', b'GPKE',
    b'GPMPF', b'GRDRM', b'GS', b'HMKT', b'IFD', b'IFG', b'IFP', b'IFS', b'IFST',
    b'IMAT', b'IMATG', b'K2DD', b'K2GG', b'K2PP', b'K4AA', b'K4KK', b'K4XX',
    b'KAA', b'KAAL', b'KDD', b'KDICTDS', b'KDICTX', b'KFHH', b'KFS', b'KGG',
    b'KGG1', b'KGGNL', b'KGGNL1', b'KGGT', b'KHH', b'KHH1', b'KKK', b'KLL',
    b'KLR', b'KMM', b'KNN', b'KOO', b'KPP', b'KRDD', b'KRFGG', b'KRR', b'KRZX',
    b'KSAZX', b'KSGG', b'KSS', b'KTTP', b'KTTS', b'KUX', b'KXWAA', b'KXX',
    b'LAJJT', b'LAM1DD', b'LAMAM', b'LAMMAT', b'LCPHL', b'LCPHP', b'LCPHX',
    b'LMPF', b'LSCM', b'LSEQ', b'LTF', b'M2DD', b'M2GG', b'M2PP', b'MA', b'MAA',
    b'MABXWGG', b'MAT', b'MAT1', b'MAT1N', b'MAT2', b'MAT2N', b'MATS', b'MATM',
    b'MBSP', b'MCHI', b'MCHI2', b'MDD', b'MDUGNI', b'MEA', b'MEF', b'MEM', b'MES',
    b'MEW', b'MFEFE', b'MFHH', b'MGG', b'MGGCOMB', b'MHH', b'MHH1', b'MI', b'MKK',
    b'MKNRGY', b'MLAM', b'MLAM2', b'MLL', b'MLR', b'MMP', b'MNRGYMTF', b'MOA',
    b'MOO', b'MPJN2O', b'MPP', b'MQG', b'MR', b'MRR', b'MSNRGY', b'MUG', b'MUGNI',
    b'MULNT', b'MUPN', b'MUX', b'MXWAA', b'MXX', b'MZZ', b'OTMT', b'P2G', b'PA',
    b'PBYG', b'PC1', b'PD', b'PDF', b'PDT', b'PDT1', b'PFP', b'PG', b'PG1',
    b'PGG', b'PGRV', b'PGT', b'PGUP', b'PGVST', b'PHA', b'PHA1', b'PHAREF1',
    b'PHASH2', b'PHDFH', b'PHDH', b'PHF', b'PHF1', b'PHG', b'PHG1', b'PHGREF',
    b'PHGREF1', b'PHT', b'PHX', b'PHXL', b'PHZ', b'PJ', b'PKF', b'PKYG', b'PL',
    b'PLI', b'PMPF', b'PMYG', b'PNL', b'PNLT', b'PO', b'POI', b'PPF', b'PPL',
    b'PPLT', b'PPT', b'PRBDOFS', b'PROPI', b'PROPO', b'PS', b'PSF', b'PSI',
    b'PST', b'PUG', b'PUGD', b'PUGS', b'PX', b'PXA', b'PXF', b'PXT', b'PXTDV',
    b'PXT1', b'PZ', b'QG', b'QHH', b'QHHL', b'QHJ', b'QHJK', b'QHJL', b'QKH',
    b'QKHL', b'QLL', b'QMG', b'QMPF', b'QPF', b'QR', b'QXX', b'R', b'R1VAL',
    b'R1VALO', b'R1VALR', b'R1VALRG', b'R2VAL', b'R2VALO', b'R2VALR',
    b'R2VALRG', b'R3VAL', b'R3VALO', b'R3VALR', b'R3VALRG', b'RBF', b'RECM',
    b'RDG', b'RESMATFT', b'RESMAX', b'RESMAX0', b'RGG', b'RHMCF', b'RMAT',
    b'RMATG', b'RMG', b'RMG1', b'RMPTQM', b'RMSVAL', b'RMSVALR', b'RMSVLR',
    b'RPH', b'RPV', b'RPX', b'RQA', b'RSPTQS', b'RSTAB', b'RUG', b'RUL', b'RUO',
    b'SCLFMAT', b'SEQMAP', b'SHPVEC', b'SKJ', b'SLIST', b'SMPF', b'SNORMM',
    b'SORTBOOL', b'SRKS', b'SRKT', b'SVEC', b'SYSE', b'TR', b'TRX', b'UA',
    b'UACCE', b'UAJJT', b'UAM1DD', b'UD', b'UD1', b'UDISP', b'UE', b'UG', b'UGD',
    b'UGDS', b'UGDS1', b'UGG', b'UGNI', b'UGNT', b'UGT', b'MATMOD', b'UGX',
    b'UGX1', b'UH', b'UHF', b'UHFF', b'UHFS', b'UI', b'UL', b'ULNT', b'UNITDISP',
    b'UO', b'UOO', b'UPF', b'UPNL0', b'UPNT', b'UTF', b'UVELO', b'UX', b'UXDIFV',
    b'UXF', b'UXR', b'UXT', b'UXT1', b'UXU', b'UXV', b'UXVBRL', b'UXVF', b'UXVP',
    b'UXVST', b'UXVW', b'VA', b'VG', b'VGD', b'WGTM', b'WJ', b'WRJVBRL',
    b'WSKJF', b'SOLVE,', b'XAA', b'XD', b'XDD', b'XDICT', b'XDICTB', b'XDICTDS',
    b'XDICTX', b'XG', b'XGG', b'XH', b'XINIT', b'XJJ', b'XO', b'XORTH', b'XP',
    b'XPP', b'SOLVIT', b'XSF', b'XSS', b'XZ', b'YACCE', b'YPF', b'YPO', b'YPT',
    b'YS', b'YS0', b'YSD', b'YVELO', b'Z1ZX', b'ZZX',

    # not sure - per BAH_Plane_cont_gust.f06 (MONITOR point deck)
    b'PMRF', b'PERF', b'PFRF', b'PGRF', b'AFRF', b'AGRF', b'MP3F',
] # type: List[bytes]

MSC_RESULT_TABLES = [b'ASSIG', b'ASEPS'] + [
    # Qualifier Information Table - which is what???
    b'QUALINFO',

    # new
    b'TOLD',
    b'RAPCONS', b'RAQCONS', b'RADCONS', b'RASCONS', b'RAFCONS', b'RAECONS',
    b'RANCONS', b'RAGCONS', b'RADEFFM', b'RAPEATC', b'RAQEATC', b'RADEATC',
    b'RASEATC', b'RAFEATC', b'RAEEATC', b'RANEATC', b'RAGEATC',

    # other
    b'MKLIST',

    # stress
    b'OES1X1', b'OES1', b'OES1X', b'OES1C', b'OESCP',
    b'OESNLXR', b'OESNLXD', b'OESNLBR', b'OESTRCP',
    b'OESNL1X', b'OESRT',
    #----------------------
    # strain
    b'OSTR1', b'OSTR1X', b'OSTR1C',

    #----------------------
    # forces
    # OEF1  - Element forces (linear elements only)
    # HOEF1 - Element heat flux
    # OEF1X - Element forces with intermediate (CBAR and CBEAM) station forces
    #         and forces on nonlinear elements
    # DOEF1 - Scaled Response Spectra
    b'OEFIT', b'OEF1X', b'OEF1', b'DOEF1',
    b'OEFITSTN', # output.op2


    # Table of Max values?
    # Table of RMS values?
    b'OEF1MX', b'OUGV1MX',

    #----------------------
    # spc forces - gset - sort 1
    b'OQG1', b'OQGV1',
    # mpc forces - gset - sort 1
    b'OQMG1', b'OQMG2',
    # ??? forces
    b'OQP1',

    #----------------------
    # displacement/velocity/acceleration/eigenvector/temperature
    # OUPV1 - Scaled Response Spectra - displacements
    b'OUG1', b'OAG1',
    b'OUGV1', b'BOUGV1', b'OUGV1PAT',
    b'OUPV1',
    b'OUG1F',  #  I think this is

    # OUGV1PAT - Displacements in the basic coordinate system
    # OUGV1  - Output (O) Displacements (U) in the global/g-set (G)
    #          coordinate system in vector (V) format and SORT1
    # BOUGV1 - Displacements in the basic coordinate system
    # BOPHIG - Eigenvectors in the basic coordinate system
    # ROUGV1 - Relative OUGV1
    # TOUGV1 - Temperature OUGV1
    b'ROUGV1', b'TOUGV1', b'RSOUGV1', b'RESOES1', b'RESEF1',

    #----------------------
    # applied loads
    # OPG1 - Applied static loads
    b'OPNL1', # nonlinear applied loads - sort 1
    b'OPG1', # applied loads - gset? - sort 1
    b'OPGV1',
    b'OPG2', # applied loads - sort 2 - v0.8

    # grid point stresses
    b'OGS1', # grid point stresses/strains - sort 1

    #-----------------------------------------------------
    # random analysis
    # OP2 tables:
    #   CRM - cumulative root mean square
    #   PSD - power spectral density function
    #   RMS - root mean square
    #   NO - number of zero crossings???
    #   ATO -autocorrelation???

    # QRG:
    #   RALL - PSDF, ATOC, CRMS
    #   ATOC - autocorrelation
    #   CRMS - cumulative root mean square
    #   PSDF - power spectral density function

    # msc displacement/velocity/acceleration
    # nx displacement
    b'OUGATO1', b'OUGCRM1', b'OUGPSD1', b'OUGRMS1', b'OUGNO1',
    b'OUGATO2', b'OUGCRM2', b'OUGPSD2', b'OUGRMS2', b'OUGNO2',

    # nx velocity
    b'OVGATO1', b'OVGCRM1', b'OVGPSD1', b'OVGRMS1', b'OVGNO1',
    b'OVGATO2', b'OVGCRM2', b'OVGPSD2', b'OVGRMS2', b'OVGNO2',

    # nx acceleration
    b'OAGATO1', b'OAGCRM1', b'OAGPSD1', b'OAGRMS1', b'OAGNO1',
    b'OAGATO2', b'OAGCRM2', b'OAGPSD2', b'OAGRMS2', b'OAGNO2',

    # msc spc/mpc forces
    # nx spc forces
    b'OQGATO1', b'OQGCRM1', b'OQGPSD1', b'OQGRMS1', b'OQGNO1',
    b'OQGATO2', b'OQGCRM2', b'OQGPSD2', b'OQGRMS2', b'OQGNO2',

    # nx mpc forces
    b'OQMATO2', b'OQMCRM2', b'OQMPSD2', b'OQMRMS2', b'OQMNO2',

    # stress
    b'OESATO1', b'OESCRM1', b'OESPSD1', b'OESRMS1', b'OESNO1', b'OESXRMS1',
    b'OESATO2', b'OESCRM2', b'OESPSD2', b'OESRMS2', b'OESNO2',

    # load vector ???
    b'OPGATO1', b'OPGCRM1', b'OPGPSD1', b'OPGRMS1', b'OPGNO1',
    b'OPGATO2', b'OPGCRM2', b'OPGPSD2', b'OPGRMS2', b'OPGNO2',

    #------------------
    # strain
    b'OSTRATO1', b'OSTRCRM1', b'OSTRPSD1', b'OSTRRMS1', b'OSTRNO1',
    b'OSTRATO2', b'OSTRCRM2', b'OSTRPSD2', b'OSTRRMS2', b'OSTRNO2',

    # force
    b'OEFATO1', b'OEFCRM1', b'OEFPSD1', b'OEFRMS1', b'OEFNO1',
    b'OEFATO2', b'OEFCRM2', b'OEFPSD2', b'OEFRMS2', b'OEFNO2',
    #b'OEFPSD2', b'OEFCRM2', b'OEFRMS2', b'OEFATO2', b'OEFNO2',

    #-----------------------------------------------------
    # other
    b'OFMPF2M',
    b'OSMPF2M', b'OPMPF2M', b'OLMPF2M', b'OGPMPF2M',

    b'OCRUG',
    b'OCRPG',

    b'STDISP', b'AEDISP', #b'TOLB2',

    # autoskip
    #b'MATPOOL',
    b'CSTM', b'CSTMS',
    b'BOPHIG',
    b'BOPG1',
    b'HOEF1',

    #------------------------------------------
    # strain energy
    b'ONRGY1',
    b'ONRGY2',

    #------------------------------------------
    # new skipped - v0.8

    # buggy
    b'IBULK',
    #b'FRL',  # frequency response list
    b'TOL',
    b'DSCM2', # normalized design sensitivity coeff. matrix

    # dont seem to crash
    b'DESCYC',
    b'DBCOPT', # design optimization history for post-processing
    b'PVT',    # table containing PARAM data
    b'XSOP2DIR',
    b'ONRGY',
    #b'CLAMA',
    b'DSCMCOL', # design sensitivity parameters
    b'CONTACTS',
    b'EDT', # aero and element deformations
    b'EXTDB',


    b'OQG2', # single point forces - sort 2 - v0.8
    b'OBC1', b'OBC2', # contact pressures and tractions at grid points
    b'OBG1', # glue normal and tangential tractions at grid points in basic coordinate system
    b'OES2', # element stresses - sort 2 - v0.8
    b'OEF2', # element forces - sort 2 - v0.8
    b'OUGV2', # absolute displacements/velocity/acceleration - sort 2

    # contact
    b'OSPDSI1', # intial separation distance
    b'OSPDS1',  # final separation distance
    b'OQGCF1', b'OQGCF2', # contact force at grid point
    b'OQGGF1', b'OQGGF2', # glue forces in grid point basic coordinate system

    b'OSPDSI2',
    b'OSPDS2',
    b'OSTR2',
    b'OESNLXR2',

    b'CMODEXT',
    b'ROUGV2',  # relative displacement
    b'CDDATA',  # cambpell diagram table
    b'OEKE1',
    b'OES1MX', # extreme stresses?
    b'OESNLBR2',
    b'BGPDTVU', # basic grid point defintion table for a superelement and related to geometry with view-grids added

    b'OUG2T',
    b'AEMONPT',
    #b'KDICT',
]
if len(MSC_RESULT_TABLES) != len(np.unique(MSC_RESULT_TABLES)):  # pragma: no cover
    counter = Counter(MSC_RESULT_TABLES)
    _MSG = 'Invalid count:\n'
    for key, cvaluei in counter.items():
        if cvaluei != 1:
            _MSG += '%s = %s\n' % (key, cvaluei)
    raise RuntimeError(_MSG)

MSC_TABLE_CONTENT = {
    # dmap 2014
    0: '',
    1: 'OUG - Displacement vector',
    2: 'OPG - Load vector',
    3: 'OQG - SPC Force vector',
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
    17: 'OUG - Acceleration vector (solutin set)',
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
    38: 'OMM - MaxMin summary',
    39: 'OQG - MPC forces',
    40: 'OGPKE - Grip point kinetic energy',
    #51: 'OFMPF2M - ???',
    #52: 'OSMPF2M - ???',
    #53: 'OPMPF2M - ???',
    #54: 'OLMPF2M - ???',
    #55: 'OGMPF2M - ???',
    89: 'OESRT - Ply strength ratios',
} # type: Dict[int, str]
MSC_OEF_REAL_MAPPER = {
    1: 3,   # CROD
    2: 1 + (10 - 1) * 11,  # CBEAM
    3: 3,   # CTUBE
    4: 17,  # CSHEAR
    10: 3,  # CONROD
    11: 2,  # CELAS1
    12: 2,  # CELAS2
    13: 2,  # CELAS3
    14: 2,  # CELAS4

    20: 2,  # CDAMP1
    21: 2,  # CDAMP2
    22: 2,  # CDAMP3
    23: 2,  # CDAMP4
    24: 3,  # CVISC

    33: 9,  # CQUAD4-centroid
    34: 9,  # CBAR-34
    35: 7,  # CCONEAX
    38: 9,  # CGAP
    #40: 8,  # CBUSH1D

    64: 2 + (11 - 2) * 5,  # CQUAD8
    #67: 14,  # CHEXA-acoustic
    #68: 14,  # CPENTA-acoustic
    69: 1 + (8 - 1) * 2,   # CBEND
    70: 2 + (11 - 2) * 4,  # CTRIAR
    74: 9,  # CTRIA3
    75: 2 + (11 - 2) * 4,  # CTRIA6
    #76: 10,  # CHEXA-acoustic/pressure
    77: 10,  # CPENTA-acoustic/pressure
    78: 10,  # CTETRA-acoustic/pressure

    82: 2 + (11 - 2) * 5,  # CQUADR
    #93: 14, # CHEXA-acoustic
    95: 9, # CQUAD4-composite
    96: 9, # CQUAD8-composite
    97: 9, # CTRIA3-composite
    98: 9, # CTRIA6-composite

    100: 8, # CBAR-100
    102: 7, # CBUSH
    #107: 8, # CHBDYE
    #108: 8, # CHBDYG
    #109: 8, # CHBDYP
    #110: 4, # CONV
    #111: 4, # CONVM
    117: 9, # WELDC
    118: 9, # WELDP
    119: 9, # CSEAM

    144: 2 + (11 - 2) * 5,  # CQUAD4
    189: 6 + (11 - 6) * 4,  # VUQUAD
}
MSC_OEF_IMAG_MAPPER = {
    1: 5,  # CROD
    2: 1 + (17 - 1) * 11,  # CBEAM
    3: 5,   # CTUBE
    4: 33,  # CSHEAR
    10: 5,  # CONROD

    11: 3,  # CELAS1
    12: 3,  # CELAS2
    13: 3,  # CELAS3
    14: 3,  # CELAS4

    20: 3,  # CDAMP1
    21: 3,  # CDAMP2
    22: 3,  # CDAMP3
    23: 3,  # CDAMP4
    24: 5,  # CVISC

    33: 17,  # CQUAD4-centroid
    34: 17,  # CBAR-34
    38: 9,   # CGAP
    #40: 8,  # CBUSH1D

    64: 2 + (19 - 2) * 5,  # CQUAD8
    #67: 14,  # CHEXA-acoustic
    #68: 14,  # CPENTA-acoustic
    69: 1 + (14 - 1) * 2,  # CBEND
    70: 2 + (19 - 2) * 4,  # CTRIAR
    74: 17,   # CTRIA3
    75: 2 + (19 - 2) * 4,  # CTRIA6

    #76: 16,  # CHEXA-acoustic velocity/pressure
    77: 16,  # CPENTA-acoustic velocity/pressure
    78: 16,  # CTETRA-acoustic velocity/pressure

    82: 2 + (19 - 2) * 5,  # CQUADR
    #93: 14, # CHEXA-acoustic

    100: 14, # CBAR-100
    102: 13, # CBUSH
    #107: 8, # CHBDYE
    #108: 8, # CHBDYG
    #109: 8, # CHBDYP
    #110: 8, # CONV
    #111: 4, # CONVM
    117: 17, # WELDC
    118: 17, # WELDP
    119: 17, # CSEAM

    144: 2 + (19 - 2) * 5,  # CQUAD4-bilinear
    189: 6 + (31 - 6) * 4,  # VUQUAD
}
