#pylint: disable=W0613,W0612,R0913
"""
Defines the OP2 class.
"""
from __future__ import (nested_scopes, generators, division, absolute_import,
                        print_function, unicode_literals)
import os
#import sys
from struct import unpack, Struct
from collections import Counter

from six import string_types, iteritems, PY2, PY3, b
from six.moves import range

from numpy import array
import numpy as np

from pyNastran import is_release
from pyNastran.f06.errors import FatalError
from pyNastran.op2.errors import SortCodeError, DeviceCodeError, FortranMarkerError
from pyNastran.op2.tables.grid_point_weight import GridPointWeight

#============================

from pyNastran.op2.tables.lama_eigenvalues.lama import LAMA
from pyNastran.op2.tables.oee_energy.onr import ONR
from pyNastran.op2.tables.ogf_gridPointForces.ogpf import OGPF

from pyNastran.op2.tables.oef_forces.oef import OEF
from pyNastran.op2.tables.oes_stressStrain.oes import OES
#from pyNastran.op2.tables.oes_stressStrain.oesm import OESM
from pyNastran.op2.tables.ogs import OGS

from pyNastran.op2.tables.opg_appliedLoads.opg import OPG
from pyNastran.op2.tables.oqg_constraintForces.oqg import OQG
from pyNastran.op2.tables.oug.oug import OUG
from pyNastran.op2.tables.ogpwg import OGPWG
from pyNastran.op2.tables.minor_tables import MinorTables
from pyNastran.op2.fortran_format import FortranFormat

from pyNastran.utils import is_binary_file
from pyNastran.utils.log import get_logger

GEOM_TABLES = [
    # GEOM2 - Table of Bulk Data entry images related to element connectivity andscalar points
    # GEOM4 - Table of Bulk Data entry images related to constraints, degree-of-freedom membership and rigid element connectivity.
    b'GEOM1', b'GEOM2', b'GEOM3', b'GEOM4',  # regular
    b'GEOM1S', b'GEOM2S', b'GEOM3S', b'GEOM4S', # superelements
    b'GEOM1N', b'GEOM1VU', b'GEOM2VU',
    b'GEOM1OLD', b'GEOM2OLD', b'GEOM4OLD',

    b'EPT', b'EPTS', b'EPTOLD',
    b'EDTS',
    b'MPT', b'MPTS',

    b'DIT', b'DITS',

    b'PVT0', b'CASECC',
    b'EDOM', b'OGPFB1',
    # GPDT  - Grid point definition table
    # BGPDT - Basic grid point definition table.
    b'GPDT', b'BGPDT', b'BGPDTS', b'BGPDTOLD',
    b'DYNAMIC', b'DYNAMICS',

    # EQEXIN - quivalence table between external and internal grid/scalaridentification numbers.
    b'EQEXIN', b'EQEXINS',
    b'ERRORN',
    b'DESTAB', b'R1TABRG',# b'HISADD',

    # eigenvalues
    b'BLAMA', b'LAMA', b'CLAMA',  #CLAMA is new
    # strain energy
    b'ONRGY1',
    # grid point weight
    b'OGPWG', b'OGPWGM',

    # other
    b'CONTACT', b'VIEWTB',
    #b'KDICT',

    # aero?
    #b'MONITOR',
]

NX_RESULT_TABLES = [
    # OESVM1  - OES Table of           element stresses for frequency response analysis that includes von Mises stress output in SORT1 format.
    # OESVM1C - OES Table of composite element stresses for frequency response analysis that includes von Mises stress output in SORT1 format.
    b'OESVM1', b'OSTRVM1',
    b'OESVM1C', b'OSTRVM1C',

    b'OESPSD2C', b'OSTPSD2C',
    b'OSTRRMS1', b'OSTRMS1C',
    b'OSTRNO1', b'OSTNO1C',
    b'MDICT', b'BDICT', b'KDICTP', b'MDICTP',

    #----------------------
    # displacement/velocity/acceleration/eigenvector/temperature
    # OUGV1  - Displacements in the global coordinate system
    b'OUGV1',
    #----------------------
    # mpc forces - gset - sort 1
    b'OQMG1',

    #----------------------
    # forces
    # OEF1X - Element forces with intermediate (CBAR and CBEAM) station forces
    #         and forces on nonlinear elements
    b'OEF1X',

    b'OES2C',
    b'OSTR2C',
    b'OSTRRMS1',
    b'OSTRMS1C',
    b'OSTRNO1',
    b'OSTNO1C',
    b'OESPSD2C',
    b'OSTPSD2C',

    b'OPHIG', # Eigenvectors in the basic coordinate system.
]

MSC_RESULT_TABLES = [b'ASSIG', b'ASEPS'] + [
    # new
    b'TOLD',
    b'RAPCONS', b'RAQCONS', b'RADCONS', b'RASCONS', b'RAFCONS', b'RAECONS',
    b'RANCONS', b'RAGCONS', b'RADEFFM', b'RAPEATC', b'RAQEATC', b'RADEATC',
    b'RASEATC', b'RAFEATC', b'RAEEATC', b'RANEATC', b'RAGEATC',

    # stress
    b'OES1X1', b'OES1', b'OES1X', b'OES1C', b'OESCP',
    b'OESNLXR', b'OESNLXD', b'OESNLBR', b'OESTRCP',
    b'OESNL1X', b'OESRT',
    #----------------------
    # strain
    b'OSTR1X', b'OSTR1C',

    #----------------------
    # forces
    # OEF1  - Element forces (linear elements only)
    # HOEF1 - Element heat flux
    # OEF1X - Element forces with intermediate (CBAR and CBEAM) station forces
    #         and forces on nonlinear elements
    # DOEF1 - Scaled Response Spectra
    b'OEFIT', b'OEF1X', b'OEF1', b'DOEF1',


    # Table of Max values?
    # Table of RMS values?
    b'OEF1MX', b'OUGV1MX',

    #----------------------
    # spc forces - gset - sort 1
    b'OQG1', b'OQGV1',
    # mpc forces - gset - sort 1
    b'OQMG1',
    # ??? forces
    b'OQP1',

    #----------------------
    # displacement/velocity/acceleration/eigenvector/temperature
    # OUPV1 - Scaled Response Spectra - displacements
    b'OUG1', b'OAG1',
    b'OUGV1', b'BOUGV1', b'OUPV1', b'OUGV1PAT',

    # OUGV1PAT - Displacements in the basic coordinate system
    # OUGV1  - Displacements in the global coordinate system
    # BOUGV1 - Displacements in the basic coordinate system
    # BOPHIG - Eigenvectors in the basic coordinate system
    # ROUGV1 - Relative OUGV1
    # TOUGV1 - Temperature OUGV1
    b'ROUGV1', b'TOUGV1', b'RSOUGV1', b'RESOES1', b'RESEF1',

    #----------------------
    # applied loads
    # OPG1 - Applied static loads
    b'OPNL1', # nonlinear applied loads - sort 1
    b'OPG1', b'OPGV1', # applied loads - gset? - sort 1
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
    b'CSTM',
    b'AXIC',
    b'BOPHIG',
    b'HOEF1',
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
]

if len(MSC_RESULT_TABLES) != len(np.unique(MSC_RESULT_TABLES)):
    counter = Counter(MSC_RESULT_TABLES)
    _MSG = 'Invalid count:\n'
    for key, value in counter.items():
        if value != 1:
            _MSG += '%s = %s\n' % (key, value)
    raise RuntimeError(_MSG)


NX_MATRIX_TABLES = [
    b'RADEFMP', # Modal Effective Inertia Matrix - Modal Matrix (per Vibrata)
    b'RAFGEN', # Load Set Modal Forces  - Modal generalized force vectors  (per Vibrata)
    b'RADAMPZ',
    b'RADAMPG',
    b'EFMFSMS', b'EFMASSS', b'RBMASSS', b'EFMFACS', b'MPFACS', b'MEFMASS', b'MEFWTS',
    b'K4HH', b'KELMP', b'MELMP',

    # not-MATPOOL
    b'DELTAK', b'DELTAM', b'RBM0', b'DELTAM0',

    # MATPOOL
    b'MRGGT', b'UEXPT',

    # MATRIX/MATPOOL - testing-remove this
    b'PATRN', b'IDENT', b'RANDM', b'CMPLX',
    b'MPATRN', b'MIDENT', b'MRANDM', b'MCMPLX',
]


MSC_MATRIX_TABLES = [
    #b'TOLD',
    b'SDT', #b'STDISP',
    b'TOLB2', b'ADSPT', #b'MONITOR',
    b'PMRT', b'PFRT', b'PGRT', # b'AEMONPT',
    b'AFRT', b'AGRT',

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
]
AUTODESK_MATRIX_TABLES = [
    #b'BELM',
    #b'KELM',
    #b'MELM',
]
# this will be split later
RESULT_TABLES = NX_RESULT_TABLES + MSC_RESULT_TABLES
MATRIX_TABLES = NX_MATRIX_TABLES + MSC_MATRIX_TABLES + AUTODESK_MATRIX_TABLES + [b'MEFF',]


class OP2_Scalar(LAMA, ONR, OGPF,
                 OEF, OES, OGS, OPG, OQG, OUG, OGPWG, MinorTables, FortranFormat):
    """
    Defines an interface for the Nastran OP2 file.
    """
    @property
    def total_effective_mass_matrix(self):
        """6x6 matrix"""
        return self.matrices['EFMFSMS']

    @property
    def effective_mass_matrix(self):
        """6x6 matrix"""
        return self.matrices['EFMASSS']

    @property
    def rigid_body_mass_matrix(self):
        """6x6 matrix"""
        return self.matrices['RBMASS']

    @property
    def modal_effective_mass_fraction(self):
        """6xnmodes matrix"""
        return self.matrices['EFMFACS']#.dataframe

    @property
    def modal_participation_factors(self):
        """6xnmodes matrix"""
        return self.matrices['MPFACS']#.dataframe

    @property
    def modal_effective_mass(self):
        """6xnmodes matrix"""
        return self.matrices['MEFMASS']#.dataframe

    @property
    def modal_effective_weight(self):
        """6xnmodes matrix"""
        return self.matrices['MEFWTS']#.dataframe

    def set_as_nx(self):
        self.is_nx = True
        self.is_msc = False
        self.is_optistruct = False
        self.is_radioss = False
        self._nastran_format = 'msc'

    def set_as_msc(self):
        self.is_nx = False
        self.is_msc = True
        self.is_optistruct = False
        self.is_radioss = False
        self._nastran_format = 'nx'

    def set_as_optistruct(self):
        self.is_nx = False
        self.is_msc = False
        self.is_optistruct = True
        self.is_radioss = False
        self._nastran_format = 'optistruct'

    def set_as_radioss(self):
        self.is_nx = False
        self.is_msc = False
        self.is_optistruct = False
        self.is_radioss = True
        self._nastran_format = 'radioss'

    def __init__(self, debug=False, log=None, debug_file=None):
        """
        Initializes the OP2_Scalar object

        Parameters
        ----------
        debug : bool; default=False
            enables the debug log and sets the debug in the logger
        log : Log()
            a logging object to write debug messages to
            (.. seealso:: import logging)
        debug_file : str; default=None (No debug)
            sets the filename that will be written to
        """
        assert isinstance(debug, bool), 'debug=%r' % debug

        self.log = get_logger(log, 'debug' if debug else 'info')
        self._count = 0
        self.op2_filename = None
        self.bdf_filename = None
        self.f06_filename = None
        self._encoding = 'utf8'

        LAMA.__init__(self)
        ONR.__init__(self)
        OGPF.__init__(self)

        OEF.__init__(self)
        OES.__init__(self)
        #OESM.__init__(self)
        OGS.__init__(self)

        OPG.__init__(self)
        OQG.__init__(self)
        OUG.__init__(self)
        OGPWG.__init__(self)
        MinorTables.__init__(self)
        FortranFormat.__init__(self)

        self.is_vectorized = False
        self._close_op2 = True

        self.result_names = set([])

        self.grid_point_weight = GridPointWeight()
        self.words = []
        self.debug = debug
        self._last_comment = None
        #self.debug = True
        #self.debug = False
        #debug_file = None
        if debug_file is None:
            self.debug_file = None
        else:
            assert isinstance(debug_file, string_types), debug_file
            self.debug_file = debug_file

    def set_as_vectorized(self, vectorized=False, ask=False):
        """don't call this...testing"""
        if vectorized is True:
            msg = 'OP2_Scalar class doesnt support vectorization.  Use OP2 '
            msg += 'from pyNastran.op2.op2 instead.'
            raise RuntimeError(msg)
        if ask is True:
            msg = 'OP2_Scalar class doesnt support ask.'
            raise RuntimeError(msg)

    def set_subcases(self, subcases=None):
        """
        Allows you to read only the subcases in the list of iSubcases

        Parameters
        ----------
        subcases : List[int, ...] / int; default=None->all subcases
            list of [subcase1_ID,subcase2_ID]
        """
        #: stores the set of all subcases that are in the OP2
        #self.subcases = set([])
        if subcases is None or subcases == []:
            #: stores if the user entered [] for iSubcases
            self.is_all_subcases = True
            self.valid_subcases = []
        else:
            #: should all the subcases be read (default=True)
            self.is_all_subcases = False

            if isinstance(subcases, int):
                subcases = [subcases]

            #: the set of valid subcases -> set([1,2,3])
            self.valid_subcases = set(subcases)
        self.log.debug("set_subcases - subcases = %s" % self.valid_subcases)

    def set_transient_times(self, times):  # TODO this name sucks...
        """
        Takes a dictionary of list of times in a transient case and
        gets the output closest to those times.

        .. code-block:: python
          times = {subcase_id_1: [time1, time2],
                   subcase_id_2: [time3, time4]}

        .. warning:: I'm not sure this still works...
        """
        expected_times = {}
        for (isubcase, etimes) in iteritems(times):
            etimes = list(times)
            etimes.sort()
            expected_times[isubcase] = array(etimes)
        self.expected_times = expected_times

    def _get_table_mapper(self):
        """gets the dictionary of function3 / function4"""
        table_mapper = {

            # per NX
            b'OESVM1' : [self._read_oes1_3, self._read_oes1_4],    # isat_random
            b'OESVM1C' : [self._read_oes1_3, self._read_oes1_4],   # isat_random
            b'OSTRVM1' : [self._read_oes1_3, self._read_ostr1_4],   # isat_random
            b'OSTRVM1C' : [self._read_oes1_3, self._read_ostr1_4],  # isat_random

            b'OSTR2' : [self._read_oes2_3, self._read_oes2_4],
            b'OES2C' : [self._read_oes2_3, self._read_oes2_4],
            b'OSTR2C' : [self._read_oes2_3, self._read_oes2_4],

            # MSC TABLES

            # common tables

            # unorganized
            b'RADCONS': [self._read_oug1_3, self._read_oug_4],     # Displacement Constraint Mode (OUG)
            b'RADEFFM': [self._read_oug1_3, self._read_oug_4], # Displacement Effective Inertia Mode (OUG)
            b'RADEATC': [self._read_oug1_3, self._read_oug_4], # Displacement Equivalent Inertia Attachment mode (OUG)

            # broken - isat_launch_100hz.op2 - wrong numwide
            #b'RAQCONS': [self._read_oqg1_3, self._read_oqg_4], # Constraint mode MPC force table (OQG)
            #b'RAQEATC': [self._read_oqg1_3, self._read_oqg_4], # Attachment mode MPC force table (OQG)
            b'RAQCONS': [self._table_passer, self._table_passer], # temporary
            b'RAQEATC': [self._table_passer, self._table_passer], # temporary

            #b'RAFCONS': [self._read_oef1_3, self._read_oef1_4], # Element Force Constraint Mode (OEF)
            #b'RAFEATC': [self._read_oef1_3, self._read_oef1_4], # Element Force Equivalent Inertia Attachment mode (OEF)
            b'RAFCONS': [self._table_passer, self._table_passer], # temporary
            b'RAFEATC': [self._table_passer, self._table_passer], # temporary

            #b'RAGCONS': [self._read_oef1_3, self._read_oef1_4], # Grid Point Forces Constraint Mode (OGPFB)
            #b'RAGEATC': [self._table_passer, self._table_passer], # Grid Point Forces Equivalent Inertia Attachment mode (OEF)
            b'RAGCONS': [self._table_passer, self._table_passer], # temporary
            b'RAGEATC': [self._table_passer, self._table_passer], # temporary

            b'RAPCONS': [self._table_passer, self._table_passer], # Constraint mode ply stress table (OES)
            b'RAPEATC': [self._table_passer, self._table_passer], # Attachment mode ply stress table (OES)

            #b'RASCONS': [self._read_oes1_3, self._read_oes1_4], # Stress Constraint Mode (OES)
            #b'RASEATC': [self._read_oes1_3, self._read_oes1_4], # Stress Equivalent Inertia Attachment mode (OES)
            b'RASCONS': [self._table_passer, self._table_passer], # temporary
            b'RASEATC': [self._table_passer, self._table_passer], # temporary

            #b'RAEEATC': [self._table_passer, self._table_passer], # Strain Equivalent Inertia Attachment mode (OES)
            #b'RAECONS': [self._read_oes1_3, self._read_oes1_4], # Strain Constraint Mode (OSTR)
            b'RAEEATC': [self._table_passer, self._table_passer], # temporary
            b'RAECONS': [self._table_passer, self._table_passer], # temporary

            b'RANEATC': [self._table_passer, self._table_passer], # Strain Energy Equivalent Inertia Attachment mode (ORGY1)
            b'RANCONS': [self._table_passer, self._table_passer], # Constraint mode element strain energy table (ORGY1)


            b'R1TABRG': [self._table_passer, self._read_r1tabrg],
            #b'TOL': [self._table_passer, self._table_passer],

            b'MATPOOL': [self._table_passer, self._table_passer], # DMIG bulk data entries
            b'CSTM':    [self._table_passer, self._table_passer],
            b'AXIC':    [self._table_passer, self._table_passer],
            b'ONRGY2':  [self._table_passer, self._table_passer],

            b'RSOUGV1': [self._table_passer, self._table_passer],
            b'RESOES1': [self._table_passer, self._table_passer],
            b'RESEF1' : [self._table_passer, self._table_passer],
            b'DESCYC' : [self._table_passer, self._table_passer],
            #b'AEMONPT' : [self._read_aemonpt_3, self._read_aemonpt_4],
            #=======================
            # OEF
            # element forces
            b'OEFIT' : [self._read_oef1_3, self._read_oef1_4],  # failure indices
            b'OEF1X' : [self._read_oef1_3, self._read_oef1_4],  # element forces at intermediate stations
            b'OEF1'  : [self._read_oef1_3, self._read_oef1_4],  # element forces or heat flux
            b'HOEF1':  [self._read_oef1_3, self._read_oef1_4], # element heat flux
            b'DOEF1' : [self._read_oef1_3, self._read_oef1_4],  # scaled response spectra - forces

            # off force
            b'OEF2'    : [self._table_passer, self._table_passer],  # element forces or heat flux
            #=======================
            # OQG
            # spc forces
            # OQG1/OQGV1 - spc forces in the nodal frame
            # OQP1 - scaled response spectra - spc-forces
            b'OQG1' : [self._read_oqg1_3, self._read_oqg_4],
            b'OQG2' : [self._read_oqg2_3, self._read_oqg_4],

            b'OQGV1' : [self._read_oqg1_3, self._read_oqg_4],
            b'OQGV2' : [self._read_oqg2_3, self._read_oqg_4],

            b'OQP1' : [self._read_oqg1_3, self._read_oqg_4],
            b'OQP2' : [self._read_oqg2_3, self._read_oqg_4],

            # SPC/MPC tables depending on table_code
            # SPC - NX/MSC
            # MPC - MSC
            b'OQGATO1' : [self._read_oqg1_3, self._read_oqg_4],
            b'OQGCRM1' : [self._read_oqg1_3, self._read_oqg_4],
            b'OQGPSD1' : [self._read_oqg1_3, self._read_oqg_4],
            b'OQGRMS1' : [self._read_oqg1_3, self._read_oqg_4],
            b'OQGNO1'  : [self._read_oqg1_3, self._read_oqg_4],

            b'OQGATO2' : [self._read_oqg2_3, self._read_oqg_4],
            b'OQGCRM2' : [self._read_oqg2_3, self._read_oqg_4],
            b'OQGPSD2' : [self._read_oqg2_3, self._read_oqg_4],
            b'OQGRMS2' : [self._read_oqg2_3, self._read_oqg_4],
            b'OQGNO2'  : [self._read_oqg2_3, self._read_oqg_4],

            #=======================
            # MPC Forces
            # these are NX tables

            # OQGM1 - mpc forces in the nodal frame
            b'OQMG1'   : [self._read_oqg1_3, self._read_oqg_mpc_forces],
            b'OQMATO1' : [self._read_oqg1_3, self._read_oqg_mpc_ato],
            b'OQMCRM1' : [self._read_oqg1_3, self._read_oqg_mpc_crm],
            b'OQMPSD1' : [self._read_oqg1_3, self._read_oqg_mpc_psd],
            b'OQMRMS1' : [self._read_oqg1_3, self._read_oqg_mpc_rms],
            b'OQMNO1'  : [self._read_oqg1_3, self._read_oqg_mpc_no],

            #b'OQMG2'   : [self._read_oqg1_3, self._read_oqg_mpc_forces],
            b'OQMATO2' : [self._read_oqg2_3, self._read_oqg_mpc_ato],
            b'OQMCRM2' : [self._read_oqg2_3, self._read_oqg_mpc_crm],
            b'OQMPSD2' : [self._read_oqg2_3, self._read_oqg_mpc_psd],
            b'OQMRMS2' : [self._read_oqg2_3, self._read_oqg_mpc_rms],
            b'OQMNO2'  : [self._read_oqg2_3, self._read_oqg_mpc_no],

            #=======================
            # OPG
            # applied loads
            b'OPG1'  : [self._read_opg1_3, self._read_opg1_4],  # applied loads in the nodal frame
            b'OPG2' : [self._table_passer, self._table_passer],

            b'OPGV1' : [self._read_opg1_3, self._read_opg1_4],  # solution set applied loads?
            b'OPNL1' : [self._read_opg1_3, self._read_opg1_4],  # nonlinear loads

            b'OPGATO1' : [self._table_passer, self._table_passer],
            b'OPGCRM1' : [self._table_passer, self._table_passer],
            b'OPGPSD1' : [self._table_passer, self._table_passer],
            b'OPGRMS1' : [self._table_passer, self._table_passer],
            b'OPGNO1'  : [self._table_passer, self._table_passer],

            b'OPGATO2' : [self._table_passer, self._table_passer],
            b'OPGCRM2' : [self._table_passer, self._table_passer],
            b'OPGPSD2' : [self._table_passer, self._table_passer],
            b'OPGRMS2' : [self._table_passer, self._table_passer],
            b'OPGNO2'  : [self._table_passer, self._table_passer],
            #=======================
            # OGPFB1
            # grid point forces
            b'OGPFB1' : [self._read_ogpf1_3, self._read_ogpf1_4],  # grid point forces

            #=======================
            # ONR/OEE
            # strain energy density
            b'ONRGY'  : [self._read_onr1_3, self._read_onr1_4],
            b'ONRGY1' : [self._read_onr1_3, self._read_onr1_4],  # strain energy density
            #=======================
            # OES
            # stress
            # OES1C - Table of composite element stresses or strains in SORT1 format
            # OESRT - Table of composite element ply strength ratio. Output by SDRCOMP
            b'OES1X1' : [self._read_oes1_3, self._read_oes1_4],  # stress - nonlinear elements
            b'OES1'   : [self._read_oes1_3, self._read_oes1_4],  # stress - linear only
            b'OES1X'  : [self._read_oes1_3, self._read_oes1_4],  # element stresses at intermediate stations & nonlinear stresses
            b'OES1C'  : [self._read_oes1_3, self._read_oes1_4],  # stress - composite
            b'OESCP'  : [self._read_oes1_3, self._read_oes1_4],  # stress - nonlinear???
            b'OESRT'  : [self._read_oes1_3, self._read_oes1_4], # ply strength ratio

            # special nonlinear tables
            # OESNLBR - Slideline stresses
            # OESNLXD - Nonlinear transient stresses
            # OESNLXR - Nonlinear stress
            #           Table of nonlinear element stresses in SORT1 format and appended for all subcases

            b'OESNLXR' : [self._read_oes1_3, self._read_oes1_4],  # nonlinear stresses
            b'OESNLXD' : [self._read_oes1_3, self._read_oes1_4],  # nonlinear transient stresses
            b'OESNLBR' : [self._read_oes1_3, self._read_oes1_4],
            b'OESNL1X' : [self._read_oes1_3, self._read_oes1_4],

            #b'OESNLXR' : [self._table_crasher, self._table_crasher],  # nonlinear stresses
            #b'OESNLXD' : [self._table_crasher, self._table_crasher],  # nonlinear transient stresses
            #b'OESNLBR' : [self._table_crasher, self._table_crasher],
            #b'OESNL1X' : [self._table_crasher, self._table_crasher],

            b'OESNLXR2' : [self._table_passer, self._table_passer],
            b'OESNLBR2' : [self._table_passer, self._table_passer],

            # off stress
            b'OES2'    : [self._table_passer, self._table_passer],  # stress - linear only
            b'OESPSD2C' : [self._table_passer, self._table_passer],
            #b'OESATO2' : [self._table_passer, self._table_passer],
            #b'OESRMS2' : [self._table_passer, self._table_passer],
            #b'OESNO2'  : [self._table_passer, self._table_passer],
            #b'OESCRM2' : [self._table_passer, self._table_passer],
            #=======================
            # strain
            b'OSTR1X'  : [self._read_oes1_3, self._read_ostr1_4],  # strain - isotropic
            b'OSTR1C'  : [self._read_oes1_3, self._read_ostr1_4],  # strain - composite
            b'OESTRCP' : [self._read_oes1_3, self._read_ostr1_4],

            # off strain
            b'OSTRATO1' : [self._table_passer, self._table_passer],
            b'OSTRCRM1' : [self._table_passer, self._table_passer],
            b'OSTRPSD1' : [self._table_passer, self._table_passer],
            b'OSTRRMS1' : [self._table_passer, self._table_passer], # isat_random
            b'OSTRNO1' : [self._table_passer, self._table_passer],  # isat_random

            b'OSTRATO2' : [self._table_passer, self._table_passer],
            b'OSTRCRM2' : [self._table_passer, self._table_passer],
            b'OSTRPSD2' : [self._table_passer, self._table_passer],
            b'OSTRRMS2' : [self._table_passer, self._table_passer],
            b'OSTRNO2'  : [self._table_passer, self._table_passer],

            b'OSTRMS1C' : [self._table_passer, self._table_passer], # isat_random
            b'OSTNO1C' : [self._table_passer, self._table_passer],  # isat_random
            #b'OSTRRMS1' : [self._read_oes1_3, self._read_oes1_4], # isat_random
            #b'OSTRNO1' : [self._read_oes1_3, self._read_oes1_4],  # isat_random
            #b'OSTRMS1C' : [self._read_oes1_3, self._read_oes1_4], # isat_random
            #b'OSTRMS1C' : [self._read_oes1_3, self._read_oes1_4], # isat_random
            #b'OSTNO1C' : [self._read_oes1_3, self._read_oes1_4],  # isat_random

            b'OSTPSD2C' : [self._table_passer, self._table_passer],
            #=======================
            # OUG
            # displacement/velocity/acceleration/eigenvector/temperature
            b'OUG1'    : [self._read_oug1_3, self._read_oug_4],  # displacements in nodal frame
            b'OAG1'    : [self._read_oug1_3, self._read_oug_4],  # accelerations in nodal frame

            b'OUGV1'   : [self._read_oug1_3, self._read_oug_4],  # displacements in nodal frame
            b'BOUGV1'  : [self._read_oug1_3, self._read_oug_4],  # OUG1 on the boundary???
            b'OUGV1PAT': [self._read_oug1_3, self._read_oug_4],  # OUG1 + coord ID
            b'OUPV1'   : [self._read_oug1_3, self._read_oug_4],  # scaled response spectra - displacement
            b'TOUGV1'  : [self._read_oug1_3, self._read_oug_4],  # grid point temperature
            b'ROUGV1'  : [self._read_oug1_3, self._read_oug_4], # relative OUG

            b'OPHIG' :  [self._read_oug1_3, self._read_oug_4],  # eigenvectors in basic coordinate system
            b'BOPHIG':  [self._read_oug1_3, self._read_oug_4],  # eigenvectors in basic coordinate system

            b'OUGV2'   : [self._read_oug2_3, self._read_oug_4],  # displacements in nodal frame
            b'ROUGV2'  : [self._read_oug2_3, self._read_oug_4], # relative OUG

            b'OUGATO1' : [self._read_oug1_3, self._read_oug_ato],
            b'OUGCRM1' : [self._read_oug1_3, self._read_oug_crm],
            b'OUGPSD1' : [self._read_oug1_3, self._read_oug_psd],
            b'OUGRMS1' : [self._read_oug1_3, self._read_oug_rms],
            b'OUGNO1'  : [self._read_oug1_3, self._read_oug_no],

            b'OUGATO2' : [self._read_oug2_3, self._read_oug_ato],
            b'OUGCRM2' : [self._read_oug2_3, self._read_oug_crm],
            b'OUGPSD2' : [self._read_oug2_3, self._read_oug_psd],
            b'OUGRMS2' : [self._read_oug2_3, self._read_oug_rms],
            b'OUGNO2'  : [self._read_oug2_3, self._read_oug_no],

            #=======================
            # extreme values of the respective table
            b'OUGV1MX' : [self._table_passer, self._table_passer],
            b'OEF1MX' : [self._table_passer, self._table_passer],
            b'OES1MX' : [self._table_passer, self._table_passer],

            #=======================
            # contact
            b'OQGCF1' : [self._table_passer, self._table_passer], # Contact force at grid point.
            b'OQGCF2' : [self._table_passer, self._table_passer], # Contact force at grid point.

            b'OSPDS1' : [self._table_passer, self._table_passer],  # Final separation distance.
            b'OSPDS2' : [self._table_passer, self._table_passer],

            b'OSPDSI1' : [self._table_passer, self._table_passer], # Initial separation distance.
            b'OSPDSI2' : [self._table_passer, self._table_passer], # Output contact separation distance results.

            b'OBC1' : [self._table_passer, self._table_passer],
            b'OBC2' : [self._table_passer, self._table_passer], # Contact pressures and tractions at grid points.

            b'OBG1' : [self._table_passer, self._table_passer], # Glue normal and tangential tractions at grid point in basic coordinate system

            b'OQGGF1' : [self._table_passer, self._table_passer], # Glue forces at grid point in basic coordinate system
            b'OQGGF2' : [self._table_passer, self._table_passer],
            #=======================
            # OGPWG
            # grid point weight
            b'OGPWG'  : [self._read_ogpwg_3, self._read_ogpwg_4],  # grid point weight
            b'OGPWGM' : [self._read_ogpwg_3, self._read_ogpwg_4],  # modal? grid point weight

            #=======================
            # OGS
            # grid point stresses
            b'OGS1' : [self._read_ogs1_3, self._read_ogs1_4],  # grid point stresses
            #=======================
            # eigenvalues
            b'BLAMA' : [self._read_buckling_eigenvalue_3, self._read_buckling_eigenvalue_4], # buckling eigenvalues
            b'CLAMA' : [self._read_complex_eigenvalue_3, self._read_complex_eigenvalue_4],   # complex eigenvalues
            b'LAMA'  : [self._read_real_eigenvalue_3, self._read_real_eigenvalue_4],         # eigenvalues

            # ===geom passers===
            # geometry
            b'GEOM1' : [self._table_passer, self._table_passer], # GEOM1-Geometry-related bulk data
            b'GEOM2' : [self._table_passer, self._table_passer], # GEOM2-element connectivity and SPOINT-related data
            b'GEOM3' : [self._table_passer, self._table_passer], # GEOM3-Static and thermal loads
            b'GEOM4' : [self._table_passer, self._table_passer], # GEOM4-constraints, DOF membership entries, MPC, and R-type element data

            # superelements
            b'GEOM1S' : [self._table_passer, self._table_passer],  # GEOMx + superelement
            b'GEOM2S' : [self._table_passer, self._table_passer],
            b'GEOM3S' : [self._table_passer, self._table_passer],
            b'GEOM4S' : [self._table_passer, self._table_passer],

            b'GEOM1VU' : [self._table_passer, self._table_passer],
            b'GEOM2VU' : [self._table_passer, self._table_passer],
            b'BGPDTVU' : [self._table_passer, self._table_passer],

            b'GEOM1N' : [self._table_passer, self._table_passer],
            b'GEOM2N' : [self._table_passer, self._table_passer],
            b'GEOM3N' : [self._table_passer, self._table_passer],
            b'GEOM4N' : [self._table_passer, self._table_passer],

            b'GEOM1OLD' : [self._table_passer, self._table_passer],
            b'GEOM2OLD' : [self._table_passer, self._table_passer],
            b'GEOM3OLD' : [self._table_passer, self._table_passer],
            b'GEOM4OLD' : [self._table_passer, self._table_passer],

            b'EPT' : [self._table_passer, self._table_passer],  # elements
            b'EPTS' : [self._table_passer, self._table_passer],  # elements - superelements
            b'EPTOLD' : [self._table_passer, self._table_passer],

            b'MPT' : [self._table_passer, self._table_passer],  # materials
            b'MPTS' : [self._table_passer, self._table_passer],  # materials - superelements

            b'DYNAMIC' : [self._table_passer, self._table_passer],
            b'DYNAMICS' : [self._table_passer, self._table_passer],
            b'DIT' : [self._table_passer, self._table_passer],
            b'DITS' : [self._table_passer, self._table_passer],

            # geometry
            #b'GEOM1': [self._read_geom1_4, self._read_geom1_4],
            #b'GEOM2': [self._read_geom2_4, self._read_geom2_4],
            #b'GEOM3': [self._read_geom3_4, self._read_geom3_4],
            #b'GEOM4': [self._read_geom4_4, self._read_geom4_4],

            # superelements
            #b'GEOM1S': [self._read_geom1_4, self._read_geom1_4],
            #b'GEOM2S': [self._read_geom2_4, self._read_geom2_4],
            #b'GEOM3S': [self._read_geom3_4, self._read_geom3_4],
            #b'GEOM4S': [self._read_geom4_4, self._read_geom4_4],

            #b'GEOM1N': [self._read_geom1_4, self._read_geom1_4],
            #b'GEOM2N': [self._read_geom2_4, self._read_geom2_4],
            #b'GEOM3N': [self._read_geom3_4, self._read_geom3_4],
            #b'GEOM4N': [self._read_geom4_4, self._read_geom4_4],

            #b'GEOM1OLD': [self._read_geom1_4, self._read_geom1_4],
            #b'GEOM2OLD': [self._read_geom2_4, self._read_geom2_4],
            #b'GEOM3OLD': [self._read_geom3_4, self._read_geom3_4],
            #b'GEOM4OLD': [self._read_geom4_4, self._read_geom4_4],

            #b'EPT' : [self._read_ept_4, self._read_ept_4],
            #b'EPTS': [self._read_ept_4, self._read_ept_4],
            #b'EPTOLD' : [self._read_ept_4, self._read_ept_4],

            #b'MPT' : [self._read_mpt_4, self._read_mpt_4],
            #b'MPTS': [self._read_mpt_4, self._read_mpt_4],

            #b'DYNAMIC': [self._read_dynamics_4, self._read_dynamics_4],
            #b'DYNAMICS': [self._read_dynamics_4, self._read_dynamics_4],
            #b'DIT': [self._read_dit_4, self._read_dit_4],   # table objects (e.g. TABLED1)

            # ===passers===
            b'EQEXIN': [self._table_passer, self._table_passer],
            b'EQEXINS': [self._table_passer, self._table_passer],

            b'GPDT' : [self._table_passer, self._table_passer],     # grid points?
            b'BGPDT' : [self._table_passer, self._table_passer],    # basic grid point defintion table
            b'BGPDTS' : [self._table_passer, self._table_passer],
            b'BGPDTOLD' : [self._table_passer, self._table_passer],

            b'PVT' : [self._table_passer, self._table_passer], # PVT - Parameter Variable Table
            b'PVT0' : [self._table_passer, self._table_passer],  # user parameter value table
            b'DESTAB' : [self._table_passer, self._table_passer],
            b'TOLD' : [self._table_passer, self._table_passer],
            b'CASECC' : [self._table_passer, self._table_passer],  # case control deck

            b'STDISP' : [self._table_passer, self._table_passer], # matrix?
            b'AEDISP' : [self._table_passer, self._table_passer], # matrix?
            #b'TOLB2' : [self._table_passer, self._table_passer], # matrix?

            # EDT - element deformation, aerodynamics, p-element, divergence analysis,
            #       and iterative solver input (includes SET1 entries)
            b'EDT' : [self._table_passer, self._table_passer],
            b'EDTS' : [self._table_passer, self._table_passer],

            b'FOL' : [self._table_passer, self._table_passer],
            #b'MONITOR' : [self._read_monitor_3, self._read_monitor_4],  # monitor points
            b'PERF' : [self._table_passer, self._table_passer],
            b'VIEWTB' : [self._table_passer, self._table_passer],   # view elements

            # DSCMCOL - Correlation table for normalized design sensitivity coefficient matrix.
            #           Output by DSTAP2.
            # DBCOPT - Design optimization history table for
            b'CONTACT' : [self._table_passer, self._table_passer],
            b'CONTACTS' : [self._table_passer, self._table_passer],
            b'OEKE1' : [self._table_passer, self._table_passer],
            b'DSCMCOL' : [self._table_passer, self._table_passer],
            b'DBCOPT' : [self._table_passer, self._table_passer],
            #b'FRL0': [self._table_passer, self._table_passer],  # frequency response list

            #==================================
            # new
            b'OFMPF2M' : [self._table_passer, self._table_passer],
            b'OLMPF2M' : [self._table_passer, self._table_passer],
            b'OPMPF2M' : [self._table_passer, self._table_passer],
            b'OSMPF2M' : [self._table_passer, self._table_passer],
            b'OGPMPF2M' : [self._table_passer, self._table_passer],

            # velocity
            b'OVGATO1' : [self._read_oug1_3, self._read_oug_ato],
            b'OVGCRM1' : [self._read_oug1_3, self._read_oug_crm],
            b'OVGPSD1' : [self._read_oug1_3, self._read_oug_psd],
            b'OVGRMS1' : [self._read_oug1_3, self._read_oug_rms],
            b'OVGNO1'  : [self._read_oug1_3, self._read_oug_no],

            b'OVGATO2' : [self._read_oug2_3, self._read_oug_ato],
            b'OVGCRM2' : [self._read_oug2_3, self._read_oug_crm],
            b'OVGPSD2' : [self._read_oug2_3, self._read_oug_psd],
            b'OVGRMS2' : [self._read_oug2_3, self._read_oug_rms],
            b'OVGNO2'  : [self._read_oug2_3, self._read_oug_no],

            #==================================
            #b'GPL': [self._table_passer, self._table_passer],
            b'OMM2' : [self._table_passer, self._table_passer],
            b'ERRORN' : [self._table_passer, self._table_passer],  # p-element error summary table
            #==================================

            b'OCRPG' : [self._table_passer, self._table_passer],
            b'OCRUG' : [self._table_passer, self._table_passer],

            b'EDOM' : [self._table_passer, self._table_passer],
            b'OUG2T' : [self._table_passer, self._table_passer],

            # acceleration
            b'OAGATO1' : [self._read_oug1_3, self._read_oug_ato],
            b'OAGCRM1' : [self._read_oug1_3, self._read_oug_crm],
            b'OAGPSD1' : [self._read_oug1_3, self._read_oug_psd],
            b'OAGRMS1' : [self._read_oug1_3, self._read_oug_rms],
            b'OAGNO1'  : [self._read_oug1_3, self._read_oug_no],

            b'OAGATO2' : [self._read_oug2_3, self._read_oug_ato],
            b'OAGCRM2' : [self._read_oug2_3, self._read_oug_crm],
            b'OAGPSD2' : [self._read_oug2_3, self._read_oug_psd],
            b'OAGRMS2' : [self._read_oug2_3, self._read_oug_rms],
            b'OAGNO2'  : [self._read_oug2_3, self._read_oug_no],

            # stress
            b'OESATO1' : [self._table_passer, self._table_passer],
            b'OESCRM1' : [self._table_passer, self._table_passer],
            b'OESPSD1' : [self._table_passer, self._table_passer],
            b'OESRMS1' : [self._table_passer, self._table_passer],
            b'OESNO1'  : [self._table_passer, self._table_passer],
            b'OESXRMS1' : [self._table_passer, self._table_passer],

            b'OESATO2' : [self._table_passer, self._table_passer],
            b'OESCRM2' : [self._table_passer, self._table_passer],
            b'OESPSD2' : [self._table_passer, self._table_passer],
            b'OESRMS2' : [self._table_passer, self._table_passer],
            b'OESNO2'  : [self._table_passer, self._table_passer],

            # force
            b'OEFATO1' : [self._read_oef1_3, self._read_oef1_4],
            b'OEFCRM1' : [self._read_oef1_3, self._read_oef1_4],
            b'OEFPSD1' : [self._read_oef1_3, self._read_oef1_4],
            b'OEFRMS1' : [self._read_oef1_3, self._read_oef1_4],
            b'OEFNO1'  : [self._read_oef1_3, self._read_oef1_4],

            b'OEFATO2' : [self._table_passer, self._table_passer],
            b'OEFCRM2' : [self._table_passer, self._table_passer],
            b'OEFPSD2' : [self._table_passer, self._table_passer],
            b'OEFRMS2' : [self._table_passer, self._table_passer],
            b'OEFNO2'  : [self._table_passer, self._table_passer],
        }
        return table_mapper

    def _not_available(self, data, ndata):
        """testing function"""
        if ndata > 0:
            raise RuntimeError('this should never be called...'
                               'table_name=%r len(data)=%s' % (self.table_name, ndata))

    def _table_crasher(self, data, ndata):
        """auto-table crasher"""
        if self.is_debug_file:
            self.binary_debug.write('  crashing table = %s\n' % self.table_name)
            raise NotImplementedError(self.table_name)
        return ndata

    def _table_passer(self, data, ndata):
        """auto-table skipper"""
        if self.is_debug_file:
            self.binary_debug.write('  skipping table = %s\n' % self.table_name)
        return ndata

    def _validate_op2_filename(self, op2_filename):
        """
        Pops a GUI if the op2_filename hasn't been set.

        Parameters
        ----------
        op2_filename : str
            the filename to check (None -> gui)

        Returns
        -------
        op2_filename : str
            a valid file string
        """
        if op2_filename is None:
            from pyNastran.utils.gui_io import load_file_dialog
            wildcard_wx = "Nastran OP2 (*.op2)|*.op2|" \
                "All files (*.*)|*.*"
            wildcard_qt = "Nastran OP2 (*.op2);;All files (*)"
            title = 'Please select a OP2 to load'
            op2_filename, wildcard_level = load_file_dialog(
                title, wildcard_wx, wildcard_qt, dirname='')
            assert op2_filename is not None, op2_filename
        return op2_filename

    def _create_binary_debug(self):
        """
        Instatiates the ``self.binary_debug`` variable/file
        """
        if hasattr(self, 'binary_debug') and self.binary_debug is not None:
            self.binary_debug.close()
            del self.binary_debug

        wb = 'w'
        if PY2:
            wb = 'wb'

        if self.debug_file is not None:
            #: an ASCII version of the op2 (creates lots of output)
            self.log.debug('debug_file = %s' % self.debug_file)
            self.binary_debug = open(self.debug_file, wb)
            self.binary_debug.write(self.op2_filename + '\n')
            self.is_debug_file = True
        else:
            self.is_debug_file = False

    def read_op2(self, op2_filename=None, combine=False):
        """
        Starts the OP2 file reading

        Parameters
        ----------
        op2_filename : str
            the op2 file

        +--------------+-----------------------+
        | op2_filename | Description           |
        +--------------+-----------------------+
        |     None     | a dialog is popped up |
        +--------------+-----------------------+
        |    string    | the path is used      |
        +--------------+-----------------------+
        """
        self._count = 0
        if self.read_mode == 1:
            #sr = list(self._results.saved)
            #sr.sort()
            #self.log.debug('_results.saved = %s' % str(sr))
            #self.log.info('_results.saved = %s' % str(sr))
            pass

        if self.read_mode != 2:
            op2_filename = self._validate_op2_filename(op2_filename)
            self.log.info('op2_filename = %r' % op2_filename)
            if not is_binary_file(op2_filename):
                if os.path.getsize(op2_filename) == 0:
                    raise IOError('op2_filename=%r is empty.' % op2_filename)
                raise IOError('op2_filename=%r is not a binary OP2.' % op2_filename)

        bdf_extension = '.bdf'
        f06_extension = '.f06'
        (fname, extension) = os.path.splitext(op2_filename)

        self.op2_filename = op2_filename
        self.bdf_filename = fname + bdf_extension
        self.f06_filename = fname + f06_extension

        self._create_binary_debug()

        #: file index
        self.n = 0
        self.table_name = None

        if not hasattr(self, 'f') or self.f is None:
            #: the OP2 file object
            self.f = open(self.op2_filename, 'rb')
            self._endian = None
            flag_data = self.f.read(20)
            self.f.seek(0)

            if unpack(b'>5i', flag_data)[0] == 4:
                self._endian = '>'
            elif unpack(b'<5i', flag_data)[0] == 4:
                self._endian = '<'
            #elif unpack(b'<ii', flag_data)[0] == 4:
                #self._endian = '<'
            else:
                # Matrices from test show
                # (24, 10, 10, 6, 2) before the Matrix Name...
                #self.show_data(flag_data, types='iqlfsld', endian='<')
                #print('----------')
                #self.show_data(flag_data, types='iqlfsld', endian='>')
                raise FatalError('cannot determine endian')
            if PY2:
                self._endian = b(self._endian)
        else:
            self._goto(self.n)


        if self.read_mode == 1:
            self._set_structs()

        #try:
        markers = self.get_nmarkers(1, rewind=True)
        #except:
            #self._goto(0)
            #try:
                #self.f.read(4)
            #except:
                #raise FatalError("The OP2 is empty.")
            #raise
        if self.is_debug_file:
            if self.read_mode == 1:
                self.binary_debug.write('read_mode = %s (vectorized; 1st pass)\n' % self.read_mode)
            elif self.read_mode == 2:
                self.binary_debug.write('read_mode = %s (vectorized; 2nd pass)\n' % self.read_mode)

        if markers == [3,]:  # PARAM, POST, -1
            if self.is_debug_file:
                self.binary_debug.write('marker = 3 -> PARAM,POST,-1?\n')
            self.post = -1
            self.read_markers([3])
            data = self.read_block()
            #assert len(data) == 12, len(data)

            self.read_markers([7])
            data = self.read_block()
            assert data == b'NASTRAN FORT TAPE ID CODE - ', '%r' % data
            if self.is_debug_file:
                self.binary_debug.write('%r\n' % data)

            data = self._read_record()
            if self.is_debug_file:
                self.binary_debug.write('%r\n' % data)
            version = data.strip()
            if version.startswith(b'NX'):
                self.set_as_nx()
                self.set_table_type()
            elif version.startswith(b'MODEP'):
                # TODO: why is this separate?
                # F:\work\pyNastran\pyNastran\master2\pyNastran\bdf\test\nx_spike\out_ac11103.op2
                self.set_as_nx()
                self.set_table_type()
            elif version.startswith(b'AEROFREQ'):
                # TODO: why is this separate?
                # C:\Users\Steve\Dropbox\pyNastran_examples\move_tpl\loadf.op2
                self.set_as_msc()
                self.set_table_type()
            elif version.startswith(b'AEROTRAN'):
                # TODO: why is this separate?
                # C:\Users\Steve\Dropbox\pyNastran_examples\move_tpl\loadf.op2
                self.set_as_msc()
                self.set_table_type()
            elif version in [b'XXXXXXXX', b'V2005R3B']:
                self.set_as_msc()
                self.set_table_type()
            elif version == b'OS12.210':
                self.set_as_optistruct()
                self.set_table_type()
            elif version == b'OS11XXXX':
                self.set_as_radioss()
                self.set_table_type()
            else:
                raise RuntimeError('unknown version=%r' % version)

            if self.is_debug_file:
                self.binary_debug.write(data.decode(self._encoding) + '\n')
            self.read_markers([-1, 0])
        elif markers == [2,]:  # PARAM, POST, -2
            if self.is_debug_file:
                self.binary_debug.write('marker = 2 -> PARAM,POST,-2?\n')
            self.post = -2
        else:
            raise NotImplementedError(markers)

        #=================
        table_name = self._read_table_name(rewind=True, stop_on_failure=False)
        if table_name is None:
            raise FatalError('There was a Nastran FATAL Error.  Check the F06.\nNo tables exist...')

        self._make_tables()
        table_names = self._read_tables(table_name)
        if self.is_debug_file:
            self.binary_debug.write('-' * 80 + '\n')
            self.binary_debug.write('f.tell()=%s\ndone...\n' % self.f.tell())
            self.binary_debug.close()
        if self._close_op2:
            self.f.close()
            del self.binary_debug
            del self.f
        #self.remove_unpickable_data()
        return table_names

    #def create_unpickable_data(self):
        #raise NotImplementedError()
        ##==== not needed ====
        ##self.f
        ##self.binary_debug

        ## needed
        #self._geom1_map
        #self._geom2_map
        #self._geom3_map
        #self._geom4_map
        #self._dit_map
        #self._dynamics_map
        #self._ept_map
        #self._mpt_map
        #self._table_mapper

    #def remove_unpickable_data(self):
        #del self.f
        #del self.binary_debug
        #del self._geom1_map
        #del self._geom2_map
        #del self._geom3_map
        #del self._geom4_map
        #del self._dit_map
        #del self._dynamics_map
        #del self._ept_map
        #del self._mpt_map
        #del self._table_mapper

    def _make_tables(self):
        return
        #global RESULT_TABLES, NX_RESULT_TABLES, MSC_RESULT_TABLES
        #table_mapper = self._get_table_mapper()
        #RESULT_TABLES = table_mapper.keys()

    def _read_tables(self, table_name):
        """
        Reads all the geometry/result tables.
        The OP2 header is not read by this function.

        Parameters
        ----------
        table_name : bytes str
            the first table's name
        """
        table_names = []
        while table_name is not None:
            table_names.append(table_name)

            if self.is_debug_file:
                self.binary_debug.write('-' * 80 + '\n')
                self.binary_debug.write('table_name = %r\n' % (table_name))

            if is_release:
                self.log.debug('  table_name=%r' % table_name)

            self.table_name = table_name
            #if 0:
                #self._skip_table(table_name)
            #else:
            if table_name in GEOM_TABLES:
                self._read_geom_table()  # DIT (agard)
            elif table_name == b'GPL':
                self._read_gpl()
            #elif table_name == b'MEFF':
                #self._read_meff()
            elif table_name == b'INTMOD':
                self._read_intmod()
            elif table_name == b'HISADD':
                self._read_hisadd()
            elif table_name == b'FRL':  # frequency response list
                self._skip_table(self.table_name)
            elif table_name == b'EXTDB':
                self._read_extdb()
            elif table_name == b'OMM2':
                self._read_omm2()
            elif table_name == b'TOL':
                self._read_tol()
            elif table_name == b'PCOMPTS': # blade
                self._read_pcompts()
            elif table_name == b'MONITOR':
                self._read_monitor()
            elif table_name == b'AEMONPT':
                self._read_aemonpt()
            elif table_name == b'FOL':
                self._read_fol()
            elif table_name == b'SDF':
                self._read_sdf()
            elif table_name in [b'IBULK', b'CDDATA']:
                self._read_ibulk()
            elif table_name == b'CMODEXT':
                self._read_cmodext()

            elif table_name in MATRIX_TABLES:
                self._read_matrix(table_name)
            elif table_name in RESULT_TABLES:
                self._read_results_table()
            elif self.skip_undefined_matrices:
                self._read_matrix(table_name)
            elif table_name.strip() in self.additional_matrices:
                self._read_matrix(table_name)
            else:
                msg = 'geom/results split: %r\n\n' % table_name
                msg += 'If you the table is a result table, see:\n'
                msg += '  model.set_additional_result_tables_to_read(methods_dict)\n'
                msg += "  methods_dict = {\n"
                msg += "      b'OUGV1' : [method3, method4],\n"
                msg += "      b'OES1X1' : False,\n"
                msg += '  }\n'
                msg += 'If you have matrices that you want to read, see:\n'
                msg += '  model.set_additional_matrices_to_read(matrices)'
                msg += '  matrices = {\n'
                msg += "      b'BHH' : True,\n"
                msg += "      b'KHH' : False,\n"
                msg += '  }\n'
                raise NotImplementedError(msg)

            table_name = self._read_table_name(rewind=True, stop_on_failure=False)
        return table_names

    def _read_tol(self):
        """
        This is probably broken for MSC Nastran

        TOL
        ---
        -2 - nitimes?
        -3 - list of times?
        """
        table_name = self._read_table_name(rewind=False, stop_on_failure=True)
        self.read_markers([-1])
        data = self._read_record()
        #self.show_data(data)

        self.read_markers([-2, 1, 0])
        #self.show_ndata(440, types='if')
        data = self._read_record()
        #print('----')
        self.read_markers([-3, 1, 0])
        #self.show_ndata(440, types='if')
        #print('----')
        self.read_markers([0])
        #data = self._read_record()


        #self.show_ndata(440, types='ifs')

        #self.show_data(data)
        #aaaa

    def _skip_table(self, table_name):
        """bypasses the next table as quickly as possible"""
        if table_name in ['DIT', 'DITS']:  # tables
            self._read_dit()
        elif table_name in ['PCOMPTS']:
            self._read_pcompts()
        else:
            self._skip_table_helper()

    def _read_dit(self):
        """
        Reads the DIT table (poorly).
        The DIT table stores information about table cards
        (e.g. TABLED1, TABLEM1).
        """
        table_name = self._read_table_name(rewind=False)
        self.read_markers([-1])
        data = self._read_record()

        self.read_markers([-2, 1, 0])
        data = self._read_record()
        table_name, = self.struct_8s.unpack(data)

        self.read_markers([-3, 1, 0])
        data = self._read_record()

        self.read_markers([-4, 1, 0])
        data = self._read_record()

        self.read_markers([-5, 1, 0])

        itable = -6
        while 1:
            markers = self.get_nmarkers(1, rewind=True)
            if markers == [0]:
                break
            data = self._read_record()
            self.read_markers([itable, 1, 0])
            itable -= 1

        #self.show(100)
        self.read_markers([0])

    def _read_table_name(self, rewind=False, stop_on_failure=True):
        """Reads the next OP2 table name (e.g. OUG1, OES1X1)"""
        table_name = None
        data = None
        if self.is_debug_file:
            self.binary_debug.write('_read_table_name - rewind=%s\n' % rewind)
        ni = self.n
        s = self.struct_8s
        if stop_on_failure:
            data = self._read_record(debug=False, macro_rewind=rewind)
            table_name, = s.unpack(data)
            if self.is_debug_file and not rewind:
                self.binary_debug.write('marker = [4, 2, 4]\n')
                self.binary_debug.write('table_header = [8, %r, 8]\n\n' % table_name)
            table_name = table_name.strip()
        else:
            try:
                data = self._read_record(macro_rewind=rewind)
                table_name, = s.unpack(data)
                table_name = table_name.strip()
            except:
                # we're done reading
                self.n = ni
                self.f.seek(self.n)

                try:
                    # we have a trailing 0 marker
                    self.read_markers([0], macro_rewind=rewind)
                except:
                    # if we hit this block, we have a FATAL error
                    raise FatalError('There was a Nastran FATAL Error.  '
                                     'Check the F06.\nlast table=%r' % self.table_name)
                table_name = None

                # we're done reading, so we're going to ignore the rewind
                rewind = False

        if rewind:
            self.n = ni
            self.f.seek(self.n)
        return table_name

    def set_additional_result_tables_to_read(self, matrices):
        """
        Adds methods to read additional result tables.
        This is expected to really only be used for skipping
        unsupported tables or disabling enabled tables that are
        buggy (e.g., OUGV1).

        Parameters
        ----------
        matrices : Dict[bytes] = varies
            a dictionary of key=name, value=list[method3, method4]/False,
            False : skips a table
                applies self._table_passer to method3 and method4
            method3 : function
                function to read table 3 results (e.g., metadata)
            method4 : function
                function to read table 4 results (e.g., the actual results)
        """
        global NX_RESULT_TABLES
        global MSC_RESULT_TABLES
        global RESULT_TABLES
        failed_keys = []
        keys = list(matrices.keys())
        for _key in keys:
            if PY3:
                if not isinstance(_key, bytes):
                    failed_keys.append(_key)
            if self.is_nx:
                NX_RESULT_TABLES.append(_key)
            else:
                MSC_RESULT_TABLES.append(_key)
        if failed_keys:
            failed_keys_str = [str(_key) for _key in failed_keys]
            raise TypeError('[%s] must be bytes' % ', '. join(failed_keys_str))
        RESULT_TABLES = NX_RESULT_TABLES + MSC_RESULT_TABLES

        #RESULT_TABLES.sort()
        #assert 'OESXRMS1' in RESULT_TABLES, RESULT_TABLES

        table_mapper = self._get_table_mapper()
        #is_added = False
        def func():
            """overloaded version of _get_table_mapper"""
            #if is_added:
                #return table_mapper
            for _key, methods in iteritems(matrices):
                if methods is False:
                    table_mapper[_key] = [self._table_passer, self._table_passer]
                else:
                    assert len(methods) == 2, methods
                    table_mapper[_key] = methods
            #is_added = True
            return table_mapper
        self._get_table_mapper = func

    def set_additional_matrices_to_read(self, matrices):
        """
        Matrices (e.g., KHH) can be sparse or dense.

        Parameters
        ----------
        matrices : Dict[bytes] = bool
            a dictionary of key=name, value=True/False,
            where True/False indicates the matrix should be read

        .. note:: If you use an already defined table (e.g. KHH), it
                  will be ignored.  If the table you requested doesn't
                  exist, there will be no effect.
        .. note:: Do not use this for result tables like OUGV1, which
                  store results like displacement.  Those are not matrices.
                  Matrices are things like DMIGs.
        """
        self.additional_matrices = matrices
        if PY2:
            self.additional_matrices = matrices
        else:
            self.additional_matrices = {}
            for matrix_name, matrix in iteritems(matrices):
                self.additional_matrices[b(matrix_name)] = matrix

    def _skip_table_helper(self):
        """
        Skips the majority of geometry/result tables as they follow a very standard format.
        Other tables don't follow this format.
        """
        self.table_name = self._read_table_name(rewind=False)
        if self.is_debug_file:
            self.binary_debug.write('skipping table...%r\n' % self.table_name)
        self.read_markers([-1])
        data = self._skip_record()
        self.read_markers([-2, 1, 0])
        data = self._skip_record()
        self._skip_subtables()

    def _read_omm2(self):
        """reads the OMM2 table"""
        self.log.debug("table_name = %r" % self.table_name)
        self.table_name = self._read_table_name(rewind=False)
        self.read_markers([-1])
        data = self._read_record()

        self.read_markers([-2, 1, 0])
        data = self._read_record()
        if len(data) == 28:
            subtable_name, month, day, year, zero, one = unpack(b(self._endian + '8s5i'), data)
            if self.is_debug_file:
                self.binary_debug.write('  recordi = [%r, %i, %i, %i, %i, %i]\n'  % (
                    subtable_name, month, day, year, zero, one))
                self.binary_debug.write('  subtable_name=%r\n' % subtable_name)
            self._print_month(month, day, year, zero, one)
        else:
            raise NotImplementedError(self.show_data(data))
        self._read_subtables()

    def _read_cmodext(self):
        r"""
        fails if a streaming block???:
         - nx_spike\mnf16_0.op2
        """
        self.table_name = self._read_table_name(rewind=False)
        self.log.debug('table_name = %r' % self.table_name)
        if self.is_debug_file:
            self.binary_debug.write('_read_geom_table - %s\n' % self.table_name)
        self.read_markers([-1])
        if self.is_debug_file:
            self.binary_debug.write('---markers = [-1]---\n')
        data = self._read_record()

        markers = self.get_nmarkers(1, rewind=True)
        if self.is_debug_file:
            self.binary_debug.write('---marker0 = %s---\n' % markers)

        marker = -2
        markers = self.read_markers([marker, 1, 0])

        data = self._read_record()
        table_name, oneseventy_a, oneseventy_b = unpack('8sii', data)
        assert oneseventy_a == 170, oneseventy_a
        assert oneseventy_b == 170, oneseventy_b
        print('170*4 =', 170*4)
        #self.show_data(data)
        marker -= 1
        marker = self._read_cmodext_helper(marker) # -3
        marker = self._read_cmodext_helper(marker)
        marker = self._read_cmodext_helper(marker)
        marker = self._read_cmodext_helper(marker)
        marker = self._read_cmodext_helper(marker)
        print('table8')
        marker = self._read_cmodext_helper(marker, debug=True)
        self.show_ndata(100)

    def _read_cmodext_helper(self, marker_orig, debug=False):
        marker = marker_orig
        #markers = self.read_nmarkers([marker, 1, 1]) # -3

        if debug:
            self.show_ndata(100)
        markers = self.get_nmarkers(3, rewind=False)
        assert markers == [marker_orig, 1, 1], markers
        print('markers =', markers)

        #marker = self.get_nmarkers(1, rewind=False, macro_rewind=False)[0]
        val_old = 0
        if debug:
            print('-----------------------------')
        i = 0
        #icheck = 7
        while 1:
            #print('i = %i' % i)
            marker = self.get_nmarkers(1, rewind=False, macro_rewind=False)[0]
            if marker != 6:
                print('marker = %s' % marker)

            assert marker == 6, marker
            data = self.read_block()
            val = unpack('i', data[:4])[0]
            if debug:
                print('val=%s delta=%s' % (val, val - val_old))
                self.show_data(data, types='ifs')
            assert len(data) > 4
            #print('i=%s val=%s delta=%s' % (i, val, val - val_old))
            val_old = val

            marker2 = self.get_nmarkers(1, rewind=True, macro_rewind=False)[0]
            #print(marker2)
            if marker2 == 696:
                break
            i += 1
        if debug:
            print('----------------------------------------')

        marker = self.get_nmarkers(1, rewind=False, macro_rewind=False)[0]
        if debug:
            print('****marker = %s' % marker)
        assert marker == 696, marker
        data = self.read_block()
        #self.show_data(data)

        marker = self.get_nmarkers(1, rewind=True, macro_rewind=False)[0]
        assert marker == (marker_orig - 1), marker

        if debug:
            self.show_ndata(200)
        return marker

        #data = self._read_record()
        #marker -= 1



        #self.show_ndata(100)

        ##marker -= 1
        ##marker_end = self.get_marker1(rewind=False)
        #asdf

    def _get_marker_n(self, nmarkers):
        """
        Gets N markers

        A marker is a flag that is used.  It's a series of 3 ints (4, n, 4)
        where n changes from marker to marker.

        Parameters
        ----------
        nmarkers : int
            the number of markers to read

        Returns
        -------
        markers : List[int, int, int]
            a list of nmarker integers
        """
        markers = []
        struc = Struct('3i')
        for i in range(nmarkers):
            block = self.f.read(12)
            marker = struc.unpack(block)
            markers.append(marker)
        return markers

    def _read_geom_table(self):
        """
        Reads a geometry table
        """
        self.table_name = self._read_table_name(rewind=False)
        if self.is_debug_file:
            self.binary_debug.write('_read_geom_table - %s\n' % self.table_name)
        self.read_markers([-1])
        data = self._read_record()

        self.read_markers([-2, 1, 0])
        data = self._read_record()
        if len(data) == 8:
            subtable_name, = self.struct_8s.unpack(data)
        else:
            strings, ints, floats = self.show_data(data)
            msg = 'Unhandled table length error\n'
            msg += 'table_name = %s\n' % self.table_name
            msg += 'len(data) = %i\n' % len(data)
            msg += 'strings  = %r\n' % strings
            msg += 'ints     = %r\n' % str(ints)
            msg += 'floats   = %r' % str(floats)
            raise NotImplementedError(msg)

        self.subtable_name = subtable_name.rstrip()
        self._read_subtables()

    def _read_results_table(self):
        """
        Reads a results table
        """
        if self.is_debug_file:
            self.binary_debug.write('read_results_table - %s\n' % self.table_name)
        self.table_name = self._read_table_name(rewind=False)
        self.read_markers([-1])
        if self.is_debug_file:
            self.binary_debug.write('---markers = [-1]---\n')
            #self.binary_debug.write('marker = [4, -1, 4]\n')
        data = self._read_record()

        self.read_markers([-2, 1, 0])
        if self.is_debug_file:
            self.binary_debug.write('---markers = [-2, 1, 0]---\n')
        data, ndata = self._read_record_ndata()
        if ndata == 8:
            subtable_name = self.struct_8s.unpack(data)
            if self.is_debug_file:
                self.binary_debug.write('  recordi = [%r]\n'  % subtable_name)
                self.binary_debug.write('  subtable_name=%r\n' % subtable_name)
        elif ndata == 28:
            subtable_name, month, day, year, zero, one = unpack(b(self._endian + '8s5i'), data)
            if self.is_debug_file:
                self.binary_debug.write('  recordi = [%r, %i, %i, %i, %i, %i]\n'  % (
                    subtable_name, month, day, year, zero, one))
                self.binary_debug.write('  subtable_name=%r\n' % subtable_name)
            self._print_month(month, day, year, zero, one)
        elif ndata == 612: # ???
            strings, ints, floats = self.show_data(data)
            msg = 'len(data) = %i\n' % ndata
            #msg += 'strings  = %r\n' % strings
            #msg += 'ints     = %r\n' % str(ints)
            #msg += 'floats   = %r' % str(floats)
            print(msg)
            subtable_name, = unpack(b(self._endian) + '8s', data[:8])
            print('subtable_name = %r' % subtable_name.strip())
        else:
            strings, ints, floats = self.show_data(data)
            msg = 'len(data) = %i\n' % ndata
            msg += 'strings  = %r\n' % strings
            msg += 'ints     = %r\n' % str(ints)
            msg += 'floats   = %r' % str(floats)
            raise NotImplementedError(msg)
        if hasattr(self, 'subtable_name'):
            raise RuntimeError('the file hasnt been cleaned up; subtable_name_old=%s new=%s' % (
                self.subtable_name, subtable_name))
        self.subtable_name = subtable_name
        self._read_subtables()

    def _print_month(self, month, day, year, zero, one):
        """
        Creates the self.date attribute from the 2-digit year.

        Parameters
        ----------
        month : int
            the month (integer <= 12)
        day :  int
            the day (integer <= 31)
        year : int
            the day (integer <= 99)
        zero : int
            a dummy integer (???)
        one : int
            a dummy integer (???)
        """
        month, day, year = self._set_op2_date(month, day, year)

        #self.log.debug("%s/%s/%4i zero=%s one=%s" % (month, day, year, zero, one))
        #if self.is_debug_file:
        if self.is_debug_file:
            self.binary_debug.write('  [subtable_name, month=%i, day=%i, year=%i, '
                                    'zero=%i, one=%i]\n\n' % (month, day, year, zero, one))
        #assert zero == 0, zero  # is this the RTABLE indicator???
        assert one in [0, 1], one  # 0, 50

    def finish(self):
        """
        Clears out the data members contained within the self.words variable.
        This prevents mixups when working on the next table, but otherwise
        has no effect.
        """
        for word in self.words:
            if word != '???' and hasattr(self, word):
                if word not in ['Title', 'reference_point']:
                    delattr(self, word)
        self.obj = None
        if hasattr(self, 'subtable_name'):
            del self.subtable_name


def main():  # pragma: no cover
    """testing pickling"""
    from pickle import dumps, dump, load, loads
    txt_filename = 'solid_shell_bar.txt'
    pickle_file = open(txt_filename, 'wb')
    op2_filename = 'solid_shell_bar.op2'
    op2 = OP2_Scalar()
    op2.read_op2(op2_filename)
    print(op2.displacements[1])
    dump(op2, pickle_file)
    pickle_file.close()

    pickle_file = open(txt_filename, 'r')
    op2 = load(pickle_file)
    pickle_file.close()
    print(op2.displacements[1])


    #import sys
    #op2_filename = sys.argv[1]

    #o = OP2_Scalar()
    #o.read_op2(op2_filename)
    #(model, ext) = os.path.splitext(op2_filename)
    #f06_outname = model + '.test_op2.f06'
    #o.write_f06(f06_outname)


if __name__ == '__main__':  # pragma: no conver
    main()
