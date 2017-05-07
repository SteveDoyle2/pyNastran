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
import scipy

from pyNastran import is_release
from pyNastran.f06.errors import FatalError
from pyNastran.op2.errors import SortCodeError, DeviceCodeError, FortranMarkerError
from pyNastran.op2.tables.grid_point_weight import GridPointWeight
from pyNastran.op2.tables.matrix import Matrix

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
from pyNastran.op2.fortran_format import FortranFormat

from pyNastran.utils import is_binary_file
from pyNastran.utils.log import get_logger
from pyNastran.op2.tables.design_response import WeightResponse, FlutterResponse, Convergence

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
    b'KDICT',

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
    b'OUG1', b'OUGV1', b'BOUGV1', b'OUPV1', b'OUGV1PAT',

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
    b'MATPOOL',
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
    b'BELM', b'KELM', b'MELM',
]
# this will be split later
RESULT_TABLES = NX_RESULT_TABLES + MSC_RESULT_TABLES
MATRIX_TABLES = NX_MATRIX_TABLES + MSC_MATRIX_TABLES + AUTODESK_MATRIX_TABLES + [b'MEFF',]


class OP2_Scalar(LAMA, ONR, OGPF,
                 OEF, OES, OGS, OPG, OQG, OUG, OGPWG, FortranFormat):
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
        self._nastran_format = 'msc'

    def set_as_msc(self):
        self.is_nx = False
        self.is_msc = True
        self.is_optistruct = False
        self._nastran_format = 'nx'

    def set_as_optistruct(self):
        self.is_nx = False
        self.is_msc = False
        self.is_optistruct = True
        self._nastran_format = 'optistruct'

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
            self.isAllSubcases = True
            self.valid_subcases = []
        else:
            #: should all the subcases be read (default=True)
            self.isAllSubcases = False

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

    def _read_aemonpt(self):
        """reads the AEMONPT table"""
        #self.log.debug("table_name = %r" % self.table_name)
        table_name = self._read_table_name(rewind=False)

        #print('-----------------------')
        #print('record 1')
        self.read_markers([-1])
        data = self._read_record()
        #self.show_data(data)
        if self.read_mode == 2:
            a, bi, c, d, e, f, g = unpack(b('%s7i'% self._endian), data)
            assert a == 101, a
            assert bi == 1, bi
            assert c == 27, c
            assert d == 1, d
            assert e == 6, e
            assert f == 0, f
            assert g == 0, g
        #print('-----------------------')
        #print('record 2')
        self.read_markers([-2, 1, 0])
        data = self._read_record()
        #if self.read_mode == 2:
        word, = unpack(b('%s8s' % self._endian), data)
        assert word == b'AECFMON ', word
        #self.show_data(data)
        #print('-----------------------')
        #print('record 3')
        self.read_markers([-3, 1, 0])
        data = self._read_record()
        #self.show_data(data)

        if self.read_mode == 2:
            ndata = len(data)
            assert ndata == 108, ndata
            n = 8 + 56 + 20 + 12 + 12
            aero, name, comps, cp, bi, c, d, coeff, word, e, f, g = unpack(b('8s 56s 5i 4s 8s 3i'), data[:n])
            print('aero=%r' % aero)
            print('name=%r' % name)
            print('comps=%r cp=%s b,c,d=(%s, %s, %s)' % (comps, cp, bi, c, d))
            print('coeff=%r' % coeff)
            print('word=%r (e, f, g)=(%s, %s, %s)' % (word, e, f, g)) # (1, 2, 0)
            assert cp == 2, cp
            assert bi == 0, bi
            assert c == 0, c
            assert d == 0, d
            assert e == 1, e
            assert f == 2, f
            assert g == 0, g

        #print('-----------------------')
        #print('record 4')
        self.read_markers([-4, 1, 0])
        #data = self._read_record()
        #self.show_data(data)
        self.read_markers([0])

        #print('-----------------------')
        #print('end')
        #self.show(200)
        #aaa

    def _read_monitor(self):
        """reads the MONITOR table"""
        self.log.debug("table_name = %r" % self.table_name)
        table_name = self._read_table_name(rewind=False)

        #print('-----------------------')
        #print('record 1')
        self.read_markers([-1])
        data = self._read_record()
        #self.show_data(data)
        if self.read_mode == 2:
            a, bi, c, d, e, f, g = unpack(b('%s7i' % self._endian), data)
            assert a == 101, a
            assert bi == 1, bi
            assert c == 27, c
            assert d == 0, d
            assert e == 6, e
            assert f == 0, f
            assert g == 0, g
        #print('-----------------------')
        #print('record 2')
        self.read_markers([-2, 1, 0])
        data = self._read_record()
        if self.read_mode == 2:
            word, = unpack(b('%s8s' % self._endian), data)
            assert word == b'STCFMON ', word
        #self.show_data(data)
        #print('-----------------------')
        #print('record 3')
        self.read_markers([-3, 1, 0])
        data = self._read_record()
        #self.show_data(data[96:108])

        if self.read_mode == 2:
            ndata = len(data)
            assert ndata == 108, ndata
            aero, name, comps, cp, x, y, z, coeff, word, column, cd, ind_dof = unpack(b'8s 56s 2i 3f 4s 8s 3i', data[:108])
            print('aero=%r' % aero)
            print('name=%r' % name)
            print('comps=%s cp=%s (x, y, z)=(%s, %s, %s)' % (comps, cp, x, y, z))
            print('coeff=%r' % coeff)
            print('word=%r (column, cd, ind_dof)=(%s, %s, %s)' % (word, column, cd, ind_dof))
            assert cp == 2, cp
            assert x == 0.0, x
            assert y == 0.0, y
            assert d == 0.0, z
            assert column == 1, column
            assert cd == 2, cd
            assert ind_dof == 0, ind_dof
            self.monitor_data = [{
                'name' : name,
                'cp' : cp,
                'cd' : cd,
                'xyz' : [x, y, z],
                'comps' : comps,
            }]

        #print('-----------------------')
        #print('record 4')
        self.read_markers([-4, 1, 0])
        #data = self._read_record()
        #self.show_data(data)
        self.read_markers([0])

        #print('-----------------------')
        #print('end')
        #self.show(200)
        #aaa

    def _table_passer(self, data, ndata):
        """auto-table skipper"""
        if self.is_debug_file:
            self.binary_debug.write('  skipping table = %s\n' % self.table_name)
        return ndata

    def _read_r1tabrg(self, data, ndata):
        """
        Design Responses:
          - Weight
          - Flutter Speed
          - Stress
          - Strain
          - Displacement
        """
        if self._table4_count == 0:
            self._count += 1
        self._table4_count += 1

        #if self._table4_count == 0:
            #self._count += 1
        #self._table4_count += 1

        if self.read_mode == 1:
            assert data is not None, data
            assert len(data) > 12, len(data)
            Type, = unpack(self._endian + 'i', data[8:12])
            #assert Type in [1, 6, 10, 84], Type
            if Type == 1:
                if self.weight_response is None:
                    self.weight_response = WeightResponse()
                else:
                    self.weight_response.n += 1
            elif Type == 4:
                #TYPE =4 EIGN or FREQ
                #8 MODE I Mode number
                #9 APRX I Approximation code
                pass
            elif Type == 5:
                #TYPE =5 DISP
                #8 COMP I Displacement component
                #9 UNDEF None
                #10 GRID I Grid identification number
                pass
            elif Type == 15:
                # CEIG
                #8 MODE I Mode number
                #9 ICODE I 1: Real component or 2: Imaginary component
                pass
            elif Type == 84:
                if self.flutter_response is None:
                    self.flutter_response = FlutterResponse()
                else:
                    self.flutter_response.n += 1
            return ndata
            #else: # response not added...
                #pass

        read_r1tabrg = True
        if read_r1tabrg:
            #self.show_data(data, types='ifs', endian=None)
            out = unpack(self._endian + 'iii 8s iiii i iiiii', data)
            # per the R1TAB DMAP page:
            #   all indicies are downshift by 1
            #   indices above out[3] are off by +2 because of the 2 field response_label
            internal_id = out[0]
            dresp_id = out[1]
            response_type = out[2]
            response_label = out[3].strip()
            # -1 for 2 field wide response_label
            region = out[4]
            subcase = out[5]
            type_flag = out[12]  # no meaning per MSC DMAP 2005
            seid = out[13]

            if response_type == 1:
                #                             -----  WEIGHT RESPONSE  -----
                # ---------------------------------------------------------------------------------
                #  INTERNAL  DRESP1  RESPONSE  ROW  COLUMN  LOWER     INPUT      OUTPUT     UPPER
                #     ID       ID     LABEL     ID    ID    BOUND     VALUE       VALUE     BOUND
                # ---------------------------------------------------------------------------------
                #       1       1    WEIGHT     3     3       N/A   2.9861E+05  2.9852E+05   N/A
                #(1, 1,    1, 'WEIGHT  ', 0, 1011, 3, 3, 0, 0, 0, 0, 0, 0)
                #(1, 1000, 1, 'W       ', 0, 1,    3, 3, 0, 0, 0, 0, 0, 0)
                #print(out)
                #row_id = out[4]

                # these should be blank?
                row_id = out[6]
                column_id = out[7]
                seid_weight = out[8]

                assert np.abs(out[8:-1]).sum() == 0.0, 'out=%s 8=%s' % (out, out[8:-1])
                assert out[-1] in [0, 1, 2, 3, 4, 5], out
                #dunno_8 = out[8]
                #dunno_9 = out[9]
                #dunno_10 = out[10]
                #dunno_11 = out[11]
                #dunno_12 = out[12]
                #dunno_13 = out[13]
                #msg = ('WEIGHT - response_type=%r response_label=%r row_id=%r column_id=%r '
                       #'6=%r 7=%r 8=%r 9=%r 10=%r 11=%r 12=%r 13=%r' % (
                           #response_type, response_label, row_id, column_id,
                           #dunno_6, dunno_7, dunno_8, dunno_9, dunno_10, dunno_11, dunno_12, dunno_13))
                #out = unpack(self._endian + 'iii 8s iiff f fffff', data)
                #print(out)
                msg = 'WEIGHT - label=%r region=%s subcase=%s row_id=%r column_id=%r' % (
                    response_label, region, subcase, row_id, column_id)
                self.weight_response.append(internal_id, dresp_id, response_label, region,
                                            subcase, type_flag, seid,
                                            row_id, column_id)
                #print(msg)
                #self.log.debug(msg)
            elif response_type == 5:  # DISP
                # out = (1, 101, 5, 'DISP1   ', 101, 1, 3, 0, 1, 0, 0, 0, 0, 0)

                #print(out[6:])
                # (3,   0,  1,    0,   0,   0,   0,   0)
                # (???, NA, comp, ???, ???, ???, ???, ???)
                pass
            elif response_type == 6:  # STRESS
                #                              -----   STRESS RESPONSES   -----
                #  -------------------------------------------------------------------------------------------
                #   INTERNAL  DRESP1  RESPONSE  ELEMENT   VIEW   COMPONENT  LOWER   INPUT    OUTPUT    UPPER
                #      ID       ID     LABEL       ID    ELM ID     NO.     BOUND   VALUE     VALUE    BOUND
                #  -------------------------------------------------------------------------------------------
                #         21      209  S09L      144747             17       N/A   4.85E+04  5.00E+04  5.00E+04
                # (21, 209, 6, 'S09L    ', 30, 1011, 17, 0, 144747, 0, 0, 0, 0, 0)
                stress_code = out[6]
                pid = out[8]
                msg = ('STRESS - response_type=%r label=%r region=%s subcase=%s '
                       'stress_code=%s pid=%s' % (
                           response_type, response_label, region, subcase,
                           stress_code, pid))

            #elif response_type == 5:  # DISP
                #pass
            #elif response_type == 7:  # STRAIN
                #pass
            elif response_type == 10:  # CSTRESS
                stress_code = out[6]
                ply = out[7]
                pid = out[8]  # is this element id?
                msg = 'CSTRESS - label=%r region=%s subcase=%s stress_code=%s ply=%s pid=%s' % (
                    response_label, region, subcase, stress_code, ply, pid)
                #print(msg)
            #elif response_type == 10:  # CSTRAIN
                #pass
            elif response_type == 24:  # FRSTRE
                #8 ICODE I Stress item code
                #9 UNDEF None
                #10 ELID I Element identification number
                #11 FREQ RS Frequency
                #12 IFLAG I Integrated response flag. See Remark 20 of
                #DRESP1.
                #Value is -1 to -6, for SUM, AVG, SSQ,
                pass
            elif response_type == 28:  # RMSACCL
                #8 COMP I RMS Acceleration component
                #9 RANDPS I RANDPS entry identification number
                #10 GRID I Grid identification number
                #11 DMFREQ RS Dummy frequency for internal use
                pass
            elif response_type == 84:  # FLUTTER  (iii, label, mode, (Ma, V, rho), flutter_id, fff)
                out = unpack(self._endian + 'iii 8s iii fff i fff', data)
                mode = out[6]
                mach = out[7]
                velocity = out[8]
                density = out[9]
                flutter_id = out[10]
                msg = ('FLUTTER - _count=%s label=%r region=%s subcase=%s mode=%s '
                       'mach=%s velocity=%s density=%s flutter_id=%s' % (
                           self._count, response_label, region, subcase, mode,
                           mach, velocity, density, flutter_id))
                self.flutter_response.append(internal_id, dresp_id, response_label, region,
                                             subcase, type_flag, seid,
                                             mode, mach, velocity, density, flutter_id)
                #print(msg)
                #self.log.debug(msg)
            else:
                self.log.debug('R1TABRG response response_type=%s not supported' % response_type)
                #raise NotImplementedError(response_type)
            assert len(out) == 14, len(out)
        #self.response1_table[self._count] = out
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
            #elif table_name in [b'DIT', b'DITS']:  # tables
                #self._read_dit()
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
            elif table_name in [b'SDF']:
                self._read_sdf()
            elif table_name in [b'IBULK', b'CDDATA']:
                self._read_ibulk()
            elif table_name in [b'CMODEXT']:
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

    def _skip_matrix_mat(self):
        table_name = self._read_table_name(rewind=False, stop_on_failure=True)
        self.read_markers([-1])
        data = self._skip_record()

        self.read_markers([-2, 1, 0])
        data = self._skip_record()

        itable = -3
        niter = 0
        niter_max = 100000000

        jj = 1
        while niter < niter_max:
            #nvalues = self.get_marker1(rewind=True)
            #print('nvalues4a =', nvalues)
            self.read_markers([itable, 1])
            one = self.get_marker1(rewind=False)

            if one:  # if keep going
                nvalues = self.get_marker1(rewind=True)
                while nvalues >= 0:
                    nvalues = self.get_marker1(rewind=False)
                    data = self._skip_block()
                    nvalues = self.get_marker1(rewind=True)
                jj += 1
            else:
                nvalues = self.get_marker1(rewind=False)
                assert nvalues == 0, nvalues
                return
            itable -= 1
            niter += 1
        raise RuntimeError('this should never happen; n=%s' % niter_max)

    def _get_matrix_row_fmt_nterms_nfloats(self, nvalues, tout):
        """
        +------+---------------------------+
        | Type | Meaning                   |
        +------+---------------------------+
        |  1   | Real, single precision    |
        |  2   | Real, double precision    |
        |  3   | Complex, single precision |
        |  4   | Complex, double precision |
        +------+---------------------------+
        """
        if tout == 1:
            nfloats = nvalues
            nterms = nvalues
            fmt = self._endian + 'i %if' % nfloats
        elif tout == 2:
            nfloats = nvalues // 2
            nterms = nvalues // 2
            fmt = self._endian + 'i %id' % nfloats
        elif tout == 3:
            nfloats = nvalues
            nterms = nvalues // 2
            fmt = self._endian + 'i %if' % nfloats
        elif tout == 4:
            nfloats = nvalues // 2
            nterms = nvalues // 4
            fmt = self._endian + 'i %id' % nfloats
        else:
            raise RuntimeError('tout = %s' % tout)
        return fmt, nfloats, nterms


    def _read_matrix(self, table_name):
        """
        general method for reading matrices and MATPOOL matrices

        .. todo:: Doesn't support checking matrices vs. MATPOOLs
        .. todo:: MATPOOLs are disabled because they're not parsed properly
        """
        i = self.f.tell()
        # if we skip on read_mode=1, we don't get debugging
        # if we just use read_mode=2, some tests fail
        #
        #if self.read_mode == 1:
        if self.read_mode == 2 and not self.debug_file:
            try:
                self._skip_matrix_mat()  # doesn't work for matpools
            except:
                self._goto(i)
                self._skip_table(table_name)
            return

        #enable_matpool = True
        #if enable_matpool:
        try:
            self._read_matrix_mat()
        except:
            # read matpool matrix
            self._goto(i)
            try:
                self._read_matpool_matrix()
            except:
                raise
                #self._goto(i)
                #self._skip_table(self.table_name)
        #else:
            #try:
                #self._read_matrix_mat()
            #except:
                #self._goto(i)
                #self._skip_table(self.table_name)

    def _read_matpool_matrix(self):
        """
        Reads a MATPOOL matrix

        MATPOOL matrices are always sparse

        +------+-----------------+
        | Form | Meaning         |
        +======+=================+
        |  1   | Square          |
        |  2   | Rectangular     |
        |  6   | Symmetric       |
        |  9   | Pseudo identity |
        +------+-----------------+
        """
        table_name = self._read_table_name(rewind=False, stop_on_failure=True)
        utable_name = table_name.decode('utf-8')
        self.read_markers([-1])
        data = self._read_record()

        self.read_markers([-2, 1, 0])
        data = self._read_record()

        self.read_markers([-3, 1, 0])
        data = self._read_record()

        nvalues = len(data) // 4
        assert len(data) % 4 == 0, len(data) / 4.

        header = unpack('3i 8s 7i', data[:48]) # 48=4*12
        assert header[:3] == (114, 1, 120), 'header[:3]=%s header=%s' % (header[:3], header)

        # ncols_gset is needed for form=9
        #  list of header values:
        #    4:5   matrix name
        #    6     placeholder
        #    7     matrix shape (1 = square, 2 or 9 = rectangular, 6 = symmetric)
        #    8     input type flag (1 = single, 2 = double, 3 = complex single,
        #                           4 = complex double)
        #    9     output type flag (0 = precision set by system cell,
        #                            1 = single, 2 = double, 3 = complex single,
        #                            4 = complex double)
        #   10     complex flag (0 = real/imaginary, >0 = magnitude/phase)
        #   11     placeholder
        #   12     number of columns in the G set (only necessary for matrix
        #                                          shape 9)
        matrix_name, junk1, matrix_shape, tin, tout, is_phase, junk2, ncols_gset = header[3:]
        matrix_name = matrix_name.strip()

        #self.log.debug('matrix_name=%s, junk1=%s, matrix_shape=%s, tin=%s, tout=%s, '
                       #'is_phase=%s, junk2=%s, ncols_gset=%s' % (
                           #matrix_name, junk1, matrix_shape, tin, tout,
                           #is_phase, junk2, ncols_gset))

        is_complex = False
        if tin > 2 or tout > 2:
            is_complex = True
            assert is_phase == 0, 'is_phase=%s' % is_phase
            imags = []

        if tout == 1:
            dtype = 'float32'
            fdtype = self.fdtype
        elif tout == 2:
            dtype = 'float64'
            fdtype = self.double_dtype
        elif tout == 3:
            dtype = 'complex64'
            fdtype = self.fdtype
        elif tout == 4:
            dtype = 'complex128'
            fdtype = self.double_dtype
        else:
            dtype = '???'
            msg = ('matrix_name=%s, junk1=%s, matrix_shape=%s, tin=%s, tout=%s, '
                   'is_phase=%s, junk2=%s, ncols_gset=%s' % (
                       matrix_name, junk1, matrix_shape, tin, tout,
                       is_phase, junk2, ncols_gset))
            self.log.warning(msg)
            raise RuntimeError(msg)

        is_symmetric = matrix_shape == 6
        is_phase_flag = is_phase > 0

        if tout in [1, 3]:
            # works for float32, complex64
            ints = np.fromstring(data[48:], dtype=self.idtype)
            floats = np.fromstring(data[48:], dtype=self.fdtype)
            temp_ints = ints
        else:
            # works for float64, complex128
            temp_ints = np.fromstring(data[48:], dtype=self.idtype)

        # find the first index with ()-1,-1)
        iminus1 = np.where(temp_ints[:-1] == -1)[0]
        double_minus1 = (iminus1[:-1] + 1 == iminus1[1:])[:-1]

        # the field after our stop
        # we'll handle the off by 1 later with arange
        istop = iminus1[:-2][double_minus1]

        # 2 fields after is the start position
        # add on a 0 to the beginning to account for the starting position
        # istart defines icol
        istart = np.hstack([0, istop[:-1] + 2])

        col_nids_short = temp_ints[istart]
        col_dofs_short = temp_ints[istart+1]
        #nj2 = len(istart)  ## TODO: why is this wrong???

        row_nids = []
        row_dofs = []
        col_nids = []
        col_dofs = []
        reals = []
        imags = []
        for col_nidi, col_dofi, istarti, istopi in zip(
            col_nids_short, col_dofs_short, istart + 2, istop):

            ## TODO: preallocate arrays
            imag = None
            # The float32/complex64 blocks are really simple
            # because we can just use the data is from the temp_ints block.
            # We calculate istart/istop and directly access the float data.
            #
            # The float64/complex128 blocks are more complicated.
            # It's easier to just use istart/istop to calculate datai
            # case that, and then slice it.
            #
            # In all cases, we use istart/istop to calculate row/col nid/dof
            # because they are always int32.  Only the real/imag types change.
            if dtype == 'float32':
                irow = np.arange(istarti, istopi-1, step=3, dtype='int32')
                real = floats[irow + 2]
            elif dtype == 'complex64':
                irow = np.arange(istarti, istopi-1, step=4, dtype='int32')
                real = floats[irow + 2]
                imag = floats[irow + 3]

            elif dtype == 'float64':
                datai = data[48+(istarti*4) : 48+(istopi*4)]
                irow = np.arange(istarti, istopi-1, step=4, dtype='int32')
                assert len(datai) % 8 == 0, len(datai) / 8
                real = np.fromstring(datai, dtype=fdtype)[1::2]

            elif dtype == 'complex128':
                datai = data[48+(istarti*4) : 48+(istopi*4)]

                # iword
                # -----
                #   0    1    3     5   <---- iword
                #   1    1    2     2   <---- nwords
                # (nid, dof, real, imag)
                irow = np.arange(istarti, istopi-1, step=6, dtype='int32')
                assert len(datai) % 8 == 0, len(datai) / 8
                floats = np.fromstring(datai, dtype=fdtype)

                # ndoubles
                # --------
                #  <---0--->   1     2    <----- iword
                #      1       1     1    <----- nwords
                # (nid, dof, real, imag)
                real = floats[1::3]
                imag = floats[2::3]
            else:
                msg = '%s is not supported' % dtype
                self.log.error(msg)
                raise RuntimeError(msg)

            if len(irow) != len(real):
                msg = 'nrow=%s nreal=%s nimag=%s' % (len(irow), len(real), len(imag))
                raise RuntimeError(msg)

            # the row index; [1, 2, ..., 43]
            row_nid = temp_ints[irow]

            # the dof; [0, 0, ..., 0.]
            row_dof = temp_ints[irow + 1]
            urow_dof = np.unique(row_dof)
            for udofi in urow_dof:
                if udofi not in [0, 1, 2, 3, 4, 5, 6]:
                    msg = 'udofi=%s is invalid; must be in [0, 1, 2, 3, 4, 5, 6]; dofs=%s' % (
                        udofi, np.asarray(urow_dof, dtype='int32').tolist())
                    raise ValueError(msg)

            ni = len(irow)
            col_nid = np.ones(ni, dtype='int32') * col_nidi
            col_dof = np.ones(ni, dtype='int32') * col_dofi

            row_nids.append(row_nid)
            row_dofs.append(row_dof)
            col_nids.append(col_nid)
            col_dofs.append(col_dof)
            reals.append(real)
            imags.append(imag)

        row_nids_array = np.hstack(row_nids)
        row_dofs_array = np.hstack(row_dofs)

        col_nids_array = np.hstack(col_nids)
        col_dofs_array = np.hstack(col_dofs)
        real_array = np.hstack(reals)
        if is_complex:
            complex_array = np.hstack(imags)
            assert len(real_array) == len(complex_array)
            real_imag_array = real_array + 1.j * complex_array
        else:
            real_imag_array = real_array

        # TODO: this is way slower than it should be
        #       because we didn't preallocate the data and the
        #       horrific grids_comp_array_to_index function
        grids1 = col_nids_array
        comps1 = col_dofs_array

        grids2 = row_nids_array
        comps2 = row_dofs_array
        assert len(grids1) == len(comps1), 'ngrids1=%s ncomps1=%s' % (len(grids1), len(comps1))
        assert len(grids1) == len(grids2), 'ngrids1=%s ngrids2=%s' % (len(grids1), len(grids2))
        assert len(comps1) == len(comps2), 'ncomps1=%s ncomps2=%s' % (len(comps1), len(comps2))

        apply_symmetry = True
        make_matrix_symmetric = apply_symmetry and matrix_shape == 'symmetric'
        j1, j2, nj1, nj2, nj = grids_comp_array_to_index(
            grids1, comps1, grids2, comps2, make_matrix_symmetric)
        assert len(j1) == len(j2), 'nj1=%s nj2=%s' % (len(j1), len(j2))
        assert len(grids1) == len(real_array), 'ngrids1=%s nreals=%s' % (len(j1), len(real_array))

        # not 100% on these, they might be flipped
        #ncols = len(np.unique(j1))
        #mrows = len(np.unique(j2))

        if is_symmetric:
            mrows = nj
            ncols = nj
            #print('  j1 =', j1)
            #print('  j2 =', j2)
        else:
            ncols = nj1
            mrows = nj2

        try:
            matrix = scipy.sparse.coo_matrix(
                (real_imag_array, (j2, j1)),
                shape=(mrows, ncols), dtype=dtype)
        except ValueError:
            msg = 'Passed all the checks; cannot build MATPOOL sparse matrix...\n'
            spaces = '                                          '
            msg += '%sname=%s dtype=%s nrows=%s ncols=%s nj1=%s nj2=%s nj=%s' % (
                spaces, table_name, dtype, mrows, ncols, nj1, nj2, nj)
            self.log.error(msg)
            raise


        # enforce symmetry if necessary
        if make_matrix_symmetric:
            # get the upper and lower triangular matrices
            upper_tri = scipy.sparse.triu(matrix)
            lower_tri = scipy.sparse.tril(matrix)

            # extracts a [1, 2, 3, ..., n] off the diagonal of the matrix
            # diagonal_array = diagional(upper_tri)
            #
            # make it a diagonal matrix
            # diagi = diags(diagonal_array)
            diagi = scipy.sparse.diags(scipy.sparse.diagional(upper_tri))

            # Check to see which triangle is populated.
            # If they both are, make sure they're equal
            # or average them and throw a warning
            lnnz = (lower_tri - diagi).nnz
            unnz = (upper_tri - diagi).nnz
            assert isinstance(lnnz, int), type(lnnz)
            assert isinstance(unnz, int), type(unnz)

            # both upper and lower triangle are populated
            if lnnz > 0 and unnz > 0:
                upper_tri_t = upper_tri.T
                if lower_tri == upper_tri_t:
                    matrix = upper_tri + upper_tri_t - diagi
                else:
                    self.log.warning(
                        'Matrix marked as symmetric does not contain '
                        'symmetric data.  Data will be symmetrized.')
                    matrix = (matrix + matrix.T) / 2.
            elif lnnz > 0:
                #  lower triangle is populated
                matrix = lower_tri + lower_tri.T - diagi
            elif unnz > 0:
                #  upper triangle is populated
                matrix = upper_tri + upper_tri_t - diagi
            else:
                # matrix is diagonal (or null)
                matrix = diagi
            data = matrix

            # matrix is symmetric, but is not stored as symmetric
            matrix_shape = 'rectangular'

        m = Matrix(table_name, is_matpool=True, form=matrix_shape)
        m.data = matrix
        m.col_nid = col_nids_array
        m.col_dof = col_dofs_array
        m.row_nid = row_nids_array
        m.row_dof = row_dofs_array
        m.form = matrix_shape
        self.matrices[utable_name] = m
        self.log.debug(m)

        self.read_markers([-4, 1, 0])
        data = self._read_record()

        if len(data) == 12:
            self.read_markers([-5, 1, 0, 0])
            return
        raise RuntimeError('failed on read_matpool_matrix')

    def _read_matrix_mat(self):
        """
        Matrix Trailer:
        +------+---------------------------------------------------+
        | Word | Contents                                          |
        +======+===================================================+
        |  1   | Number of columns in matrix                       |
        |  2   | Number of rows in matrix                          |
        |  3   | Form of the matrix                                |
        |  4   | Type of matrix                                    |
        |  5   | Largest number of nonzero words among all columns |
        |  6   | Density of the matrix multiplied by 10000         |
        |  7   | Size in blocks                                    |
        |  8   | Maximum string length over all strings            |
        |  9   | Number of strings                                 |
        |  10  | Average bandwidth                                 |
        |  11  | Maximum bandwidth                                 |
        |  12  | Number of null columns                            |
        +------+---------------------------------------------------+

        +------+--------------------------------+
        | Form | Meaning                        |
        +======+================================+
        |  1   | Square                         |
        |  2   | Rectangular                    |
        |  3   | Diagonal                       |
        |  4   | Lower triangular factor        |
        |  5   | Upper triangular factor        |
        |  6   | Symmetric                      |
        |  8   | Identity                       |
        |  9   | Pseudo identity                |
        |  10  | Cholesky factor                |
        |  11  | Trapezoidal factor             |
        |  13  | Sparse lower triangular factor |
        |  15  | Sparse upper triangular factor |
        +------+--------------------------------+

        +------+---------------------------+
        | Type | Meaning                   |
        +======+===========================+
        |  1   | Real, single precision    |
        |  2   | Real, double precision    |
        |  3   | Complex, single precision |
        |  4   | Complex, double precision |
        +------+---------------------------+
        """
        allowed_forms = [1, 2, 3, 4, 5, 6, 8, 9, 10, 11, 13, 15]
        #self.log.debug('----------------------------------------------------------------')
        table_name = self._read_table_name(rewind=False, stop_on_failure=True)
        self.read_markers([-1])
        data = self._read_record()

        # old-bad
        #matrix_num, form, mrows, ncols, tout, nvalues, g = unpack(self._endian + '7i', data)

        #           good   good   good  good  ???    ???
        matrix_num, ncols, mrows, form, tout, nvalues, g = unpack(self._endian + '7i', data)
        #print('g =', g)

        m = Matrix(table_name, form=form)
        self.matrices[table_name.decode('utf-8')] = m

        # matrix_num is a counter (101, 102, 103, ...)
        # 101 will be the first matrix 'A' (matrix_num=101),
        # then we'll read a new matrix 'B' (matrix_num=102),
        # etc.
        #
        # the matrix is Mrows x Ncols
        #
        # it has nvalues in it
        #
        # tout is the precision of the matrix
        # 0 - set precision by cell
        # 1 - real, single precision (float32)
        # 2 - real, double precision (float64)
        # 3 - complex, single precision (complex64)
        # 4 - complex, double precision (complex128)

        # form (bad name)
        # 1 - column matrix
        # 2 - factor matrix
        # 3 - factor matrix
        if tout == 1:
            dtype = 'float32'
        elif tout == 2:
            dtype = 'float64'
        elif tout == 3:
            dtype = 'complex64'
        elif tout == 4:
            dtype = 'complex128'
        else:
            dtype = '???'
            msg = ('unexpected tout for %s: matrix_num=%s form=%s '
                   'mrows=%s ncols=%s tout=%s nvalues=%s g=%s'  % (
                       table_name, matrix_num, form, mrows, ncols, tout, nvalues, g))
            self.log.warning(msg)
            raise RuntimeError(msg)

        if form == 1:
            if ncols != mrows:
                self.log.warning('unexpected size for %s; form=%s mrows=%s ncols=%s' % (
                    table_name, form, mrows, ncols))
        elif form not in allowed_forms:
            self.log.error('name=%r matrix_num=%s form=%s mrows=%s '
                           'ncols=%s tout=%s nvalues=%s g=%s' % (
                               table_name, matrix_num, form, mrows, ncols, tout, nvalues, g))
            raise RuntimeError('form=%s; allowed=%s' % (form, allowed_forms))
        #self.log.info('name=%r matrix_num=%s form=%s mrows=%s ncols=%s tout=%s nvalues=%s g=%s' % (
            #table_name, matrix_num, form, mrows, ncols, tout, nvalues, g))

        self.read_markers([-2, 1, 0])
        data = self._read_record()

        if len(data) == 16:
            name, ai, bi = unpack(self._endian + '8s 2i', data)
            assert ai == 170, ai
            assert bi == 170, bi
        else:
            self.log.warning('unexpected matrix length=%s' % len(data))
            self.log.warning(self.show_data(data, types='if'))

        itable = -3
        j = None

        niter = 0
        niter_max = 100000000

        GCi = []
        GCj = []
        reals = []
        jj = 1
        while niter < niter_max:
            #nvalues = self.get_marker1(rewind=True)
            self.read_markers([itable, 1])
            one = self.get_marker1(rewind=False)

            if one:  # if keep going
                nvalues = self.get_marker1(rewind=True)

                while nvalues >= 0:
                    nvalues = self.get_marker1(rewind=False)
                    fmt, nfloats, nterms = self._get_matrix_row_fmt_nterms_nfloats(nvalues, tout)
                    GCjj = [jj] * nterms
                    GCj += GCjj

                    #-----------
                    data = self.read_block()
                    #self.show_data(data)
                    #print('***itable=%s nvalues=%s fmt=%r' % (itable, nvalues, fmt))
                    out = unpack(fmt, data)
                    ii = out[0]
                    values = out[1:]

                    GCii = list(range(ii, ii + nterms))
                    GCi += GCii
                    reals += values
                    nvalues = self.get_marker1(rewind=True)
                    if self.debug_file:
                        self.binary_debug.write('  GCi = %s\n' % GCii)
                        self.binary_debug.write('  GCj = %s\n' % GCjj)
                        self.binary_debug.write('  reals/imags = %s\n' % str(values))
                assert len(GCi) == len(GCj), 'nGCi=%s nGCj=%s' % (len(GCi), len(GCj))
                if tout in [1, 2]:
                    assert len(GCi) == len(reals), 'tout=%s nGCi=%s nReals=%s' % (tout, len(GCi), len(reals))
                else:
                    assert len(GCi)*2 == len(reals), 'tout=%s nGCi=%s nReals=%s' % (tout, len(GCi)*2, len(reals))
                jj += 1
            else:
                nvalues = self.get_marker1(rewind=False)
                assert nvalues == 0, nvalues
                # print('nvalues =', nvalues)
                # print('returning...')

                #assert max(GCi) <= mrows, 'GCi=%s GCj=%s mrows=%s' % (GCi, GCj, mrows)
                #assert max(GCj) <= ncols, 'GCi=%s GCj=%s ncols=%s' % (GCi, GCj, ncols)
                GCi = array(GCi, dtype='int32') - 1
                GCj = array(GCj, dtype='int32') - 1
                try:
                    # we subtract 1 to the indicides to account for Fortran
                    #    huh??? we don't...
                    if dtype == '???':
                        matrix = None
                        self.log.warning('what is the dtype?')
                    elif tout in [1, 2]:
                        real_array = np.array(reals, dtype=dtype)
                        matrix = scipy.sparse.coo_matrix(
                            (real_array, (GCi, GCj)),
                            shape=(mrows, ncols), dtype=dtype)
                        matrix = matrix.todense()
                        #self.log.info('created %s' % self.table_name)
                    elif tout in [3, 4]:
                        real_array = np.array(reals, dtype=dtype)
                        nvalues_matrix = real_array.shape[0] // 2
                        real_complex = real_array.reshape((nvalues_matrix, 2))
                        real_imag = real_complex[:, 0] + real_complex[:, 1]*1j
                        if self.binary_debug:
                            #self.binary_debug.write('reals = %s' % real_complex[:, 0])
                            #self.binary_debug.write('imags = %s' % real_complex[:, 1])
                            self.binary_debug.write('real_imag = %s' % real_imag)
                        matrix = scipy.sparse.coo_matrix(
                            (real_imag, (GCi, GCj)),
                            shape=(mrows, ncols), dtype=dtype)
                        msg = 'created %s...verify the complex matrix' % self.table_name
                        self.log.warning(msg)
                        #raise RuntimeError(msg)
                    else:
                        raise RuntimeError('this should never happen')
                except ValueError:
                    self.log.warning('shape=(%s, %s)' % (mrows, ncols))
                    self.log.warning('cant make a coo/sparse matrix...trying dense')

                    if dtype == '???':
                        matrix = None
                        self.log.warning('what is the dtype?')
                    else:
                        real_array = np.array(reals, dtype=dtype)
                        self.log.debug('shape=%s mrows=%s ncols=%s' % (
                            str(real_array.shape), mrows, ncols))
                        if len(reals) == mrows * ncols:
                            real_array = real_array.reshape(mrows, ncols)
                            self.log.info('created %s' % self.table_name)
                        else:
                            self.log.warning('cant reshape because invalid sizes : created %s' %
                                             self.table_name)

                        matrix = real_array

                m.data = matrix
                if matrix is not None:
                    self.matrices[table_name.decode('utf-8')] = m
                #nvalues = self.get_marker1(rewind=True)
                return
            itable -= 1
            niter += 1
        raise RuntimeError('MaxIteration: this should never happen; n=%s' % niter_max)

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

    def _skip_pcompts(self):
        """
        Reads the PCOMPTS table (poorly).
        The PCOMPTS table stores information about the PCOMP cards???
        """
        self.log.debug("table_name = %r" % self.table_name)
        table_name = self._read_table_name(rewind=False)

        self.read_markers([-1])
        data = self._skip_record()

        self.read_markers([-2, 1, 0])
        data = self._skip_record()
        #table_name, = self.struct_8s.unpack(data)

        self.read_markers([-3, 1, 0])
        markers = self.get_nmarkers(1, rewind=True)
        if markers != [-4]:
            data = self._skip_record()

        self.read_markers([-4, 1, 0])
        markers = self.get_nmarkers(1, rewind=True)
        if markers != [0]:
            data = self._skip_record()
        else:
            self.read_markers([0])
            return

        self.read_markers([-5, 1, 0])
        data = self._skip_record()

        self.read_markers([-6, 1, 0])
        self.read_markers([0])

    def _read_pcompts(self):
        """
        Reads the PCOMPTS table (poorly).
        The PCOMPTS table stores information about the PCOMP cards???
        """
        self._skip_pcompts()
        return
        #if self.read_mode == 1:
            #return
        #self.log.debug("table_name = %r" % self.table_name)
        #table_name = self._read_table_name(rewind=False)

        #self.read_markers([-1])
        #data = self._read_record()

        #self.read_markers([-2, 1, 0])
        #data = self._read_record()
        #table_name, = self.struct_8s.unpack(data)
        ##print "table_name = %r" % table_name

        #self.read_markers([-3, 1, 0])
        #markers = self.get_nmarkers(1, rewind=True)
        #if markers != [-4]:
            #data = self._read_record()

        #self.read_markers([-4, 1, 0])
        #markers = self.get_nmarkers(1, rewind=True)
        #if markers != [0]:
            #data = self._read_record()
        #else:
            #self.read_markers([0])
            #return

        #self.read_markers([-5, 1, 0])
        #data = self._read_record()

        #self.read_markers([-6, 1, 0])
        #self.read_markers([0])

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

    def _read_fol(self):
        """
        Reads the FOL table
        Frequency response frequency output list

        +------+---------+-------+-----------------+
        | Word |  Name   | Type  |   Description   |
        +======+=========+=======+=================+
        |  1   | NAME(2) | CHAR4 | Data block name |
        +------+---------+-------+-----------------+
        |  3   |  FREQ   |  RS   |   Frequency     |
        +------+---------+-------+-----------------+
        | Word 3 repeats until End of Record       |
        +------------------------------------------+

        +------+----------+------+-----------------------------+
        | Word |  Name    | Type |   Description               |
        +======+==========+======+=============================+
        |  1   |  WORD1   |  I   | Number of frequencies       |
        +------+----------+------+-----------------------------+
        |  2   |  WORD2   |  I   | Frequency set record number |
        +------+----------+------+-----------------------------+
        |  3   |  WORD3   |  I   | Number of loads             |
        +------+----------+------+-----------------------------+
        |  4   | UNDEF(3) | None | Not used                    |
        +------+----------+------+-----------------------------+
        """
        self.log.debug("table_name = %r" % self.table_name)
        self.table_name = self._read_table_name(rewind=False)
        self.read_markers([-1])
        data = self._read_record()
        self.read_markers([-2, 1, 0])
        data = self._read_record()
        ndata = len(data)
        subtable_name_raw, = unpack(b(self._endian + '8s'), data[:8])
        subtable_name = subtable_name_raw.strip()
        assert subtable_name == b'FOL', 'subtable_name=%r' % subtable_name

        nfloats = (ndata - 8) // 4
        assert nfloats * 4 == (ndata - 8)
        fmt = b(self._endian + '%sf' % nfloats)
        freqs = np.array(list(unpack(fmt, data[8:])), dtype='float32')
        self._frequencies = freqs
        if self.is_debug_file:
            self.binary_debug.write('  recordi = [%r, freqs]\n'  % (subtable_name_raw))
            self.binary_debug.write('  subtable_name=%r\n' % subtable_name)
            self.binary_debug.write('  freqs = %s' % freqs)
        self._read_subtables()

    def _read_gpl(self):
        """reads the GPL table (grid point list?)"""
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
        n = -2
        while markers[0] != 0:
            self.read_markers([n, 1, 0])
            if self.is_debug_file:
                self.binary_debug.write('---markers = [%i, 1, 0]---\n' % n)

            markers = self.get_nmarkers(1, rewind=True)
            if markers[0] == 0:
                markers = self.get_nmarkers(1, rewind=False)
                break
            data = self._read_record()
            #self.show_data(data, 'i')
            n -= 1
            markers = self.get_nmarkers(1, rewind=True)

    def _read_extdb(self):
        r"""
        fails if a streaming block:
         - nx_spike\extse04c_0.op2
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
        while 1:
            try:
                self.read_markers([marker, 1, 0])
            except FortranMarkerError:
                self.show_ndata(100)
                raise
            nfields = self.get_marker1(rewind=True)
            if nfields > 0:
                data = self._read_record()
                #self.show_data(data, types='s', endian=None)
            elif nfields == 0:
                #self.show_ndata(100, types='ifs')
                break
            else:
                raise RuntimeError('nfields=%s' % nfields)
            marker -= 1
        marker_end = self.get_marker1(rewind=False)

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

    def _read_ibulk(self):
        self.table_name = self._read_table_name(rewind=False)
        self.log.debug('table_name = %r' % self.table_name)
        if self.is_debug_file:
            self.binary_debug.write('_read_geom_table - %s\n' % self.table_name)
        self.read_markers([-1])
        if self.is_debug_file:
            self.binary_debug.write('---markers = [-1]---\n')
        data = self._read_record()

        markers = self.get_nmarkers(1, rewind=True)
        marker = -2
        while 1:
            self.read_markers([marker, 1, 0])
            nfields = self.get_marker1(rewind=True)
            if nfields > 0:
                data = self._read_record()
                #self.show_data(data, types='s', endian=None)
            elif nfields == 0:
                #self.show_ndata(100, types='ifs')
                break
            else:
                raise RuntimeError('nfields=%s' % nfields)
            marker -= 1
        marker_end = self.get_marker1(rewind=False)

    def _read_meff(self):
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
        self.read_markers([-2, 1, 0])
        data = self._read_record()

        for n in [-3, -4, -5, -6, -7, -8]:
            self.read_markers([n, 1, 1])
            markers = self.get_nmarkers(1, rewind=False)
            #print('markers =', markers)
            nbytes = markers[0]*4 + 12
            data = self.f.read(nbytes)
            self.n += nbytes
        n = -9
        self.read_markers([n, 1, 0, 0])

    def _read_intmod(self):
        """reads the INTMOD table"""
        self.table_name = self._read_table_name(rewind=False)
        #self.log.debug('table_name = %r' % self.table_name)
        if self.is_debug_file:
            self.binary_debug.write('_read_geom_table - %s\n' % self.table_name)
        self.read_markers([-1])
        if self.is_debug_file:
            self.binary_debug.write('---markers = [-1]---\n')
        data = self._read_record()
        #print('intmod data1')
        #self.show_data(data)

        markers = self.get_nmarkers(1, rewind=True)
        if self.is_debug_file:
            self.binary_debug.write('---marker0 = %s---\n' % markers)
        self.read_markers([-2, 1, 0])
        data = self._read_record()
        #print('intmod data2')
        #self.show_data(data)

        for n in [-3, -4, -5, -6, -7, -8,]:
            self.read_markers([n, 1, 1])
            markers = self.get_nmarkers(1, rewind=False)
            #print('markers =', markers)
            nbytes = markers[0]*4 + 12
            data = self.f.read(nbytes)
            #print('intmod data%i' % n)
            #self.show_data(data)
            self.n += nbytes

        n = -9
        self.read_markers([n, 1, 0, 0])
        #self.show(50)
        #raise NotImplementedError(self.table_name)


    def _read_hisadd(self):
        """optimization history (SOL200) table"""
        self.table_name = self._read_table_name(rewind=False)

        if self.read_mode == 1:
            self.read_markers([-1])
            self._skip_record()
            self.read_markers([-2, 1, 0])
            self._skip_record()
            self.read_markers([-3, 1, 0])

            if self.convergence_data is None:
                data = self._read_record()
                ndvs = len(data) // 4 - 7
                self.convergence_data = Convergence(ndvs)
            else:
                self._skip_record()
                self.convergence_data.n += 1

            self.read_markers([-4, 1, 0, 0])
            return

        if self.is_debug_file:
            self.binary_debug.write('_read_geom_table - %s\n' % self.table_name)
        #self.log.info('----marker1----')
        self.read_markers([-1])
        if self.is_debug_file:
            self.binary_debug.write('---markers = [-1]---\n')
        data = self._read_record()  # ()102, 303, 0, 0, 0, 0, 0) date???
        #print('hisadd data1')
        #self.show_data(data)

        #self.log.info('----marker2----')
        markers = self.get_nmarkers(1, rewind=True)
        if self.is_debug_file:
            self.binary_debug.write('---marker0 = %s---\n' % markers)
        self.read_markers([-2, 1, 0])
        data = self._read_record()  # ('HISADD', )
        #print('hisadd data2')
        #self.show_data(data)

        #self.log.info('----marker3----')
        self.read_markers([-3, 1, 0])
        data = self._read_record()

        (design_iter, iconvergence, conv_result, obj_intial, obj_final,
         constraint_max, row_constraint_max) = unpack(b(self._endian + '3i3fi'), data[:28])
        if iconvergence == 1:
            iconvergence = 'soft'
        elif iconvergence == 2:
            iconvergence = 'hard'
        elif iconvergence == 6:
            self.log.warning('HISADD iconverge=6')
            iconvergence = '???'
        else:
            msg = 'iconvergence=%s\n' % iconvergence
            self.show_data(data, types='ifs', endian=None)
            raise NotImplementedError(msg)

        if conv_result == 0:
            conv_result = 'no'
        elif conv_result == 1:
            conv_result = 'soft'
        elif conv_result == 2:
            conv_result = 'hard'
        elif conv_result in [3, 4]:
            self.log.warning('HISADD conv_result=%s' % conv_result)
            # not sure why this happens, but the field is wrong
            # it seems to apply to one step before this one
            conv_result = 'best_design'
        else:
            self.log.debug('design_iter=%s iconvergence=%s conv_result=%s obj_intial=%s '
                           'obj_final=%s constraint_max=%s row_constraint_max=%s' % (
                               design_iter, iconvergence, conv_result, obj_intial,
                               obj_final, constraint_max, row_constraint_max))
            raise NotImplementedError('conv_result=%s' % conv_result)
        #self.log.debug('design_iter=%s iconvergence=%s conv_result=%s obj_intial=%s '
                       #'obj_final=%s constraint_max=%s row_constraint_max=%s' % (
                           #design_iter, iconvergence, conv_result, obj_intial,
                           #obj_final, constraint_max, row_constraint_max))

        ndvs = len(data) // 4 - 7
        desvar_values = unpack('%sf' % ndvs, data[28:])

        self.convergence_data.append(design_iter, iconvergence, conv_result, obj_intial,
                                     obj_final, constraint_max, row_constraint_max, desvar_values)
        self.read_markers([-4, 1, 0, 0])


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

    def _read_frl(self):
        #self.log.debug("table_name = %r" % self.table_name)
        #self.table_name = self._read_table_name(rewind=False)
        #self.read_markers([-1])
        #data = self._read_record()

        #self.read_markers([-2, 1, 0])
        #data = self._read_record()
        self._skip_table(self.table_name)


    def _read_sdf(self):
        self.log.debug("table_name = %r" % self.table_name)
        self.table_name = self._read_table_name(rewind=False)
        self.read_markers([-1])
        data = self._read_record()

        self.read_markers([-2, 1, 0])
        data, ndata = self._read_record_ndata()
        if ndata == 16:
            subtable_name, dummy_a, dummy_b = unpack(b(self._endian + '8sii'), data)
            if self.is_debug_file:
                self.binary_debug.write('  recordi = [%r, %i, %i]\n'  % (
                    subtable_name, dummy_a, dummy_b))
                self.binary_debug.write('  subtable_name=%r\n' % subtable_name)
                assert dummy_a == 170, dummy_a
                assert dummy_b == 170, dummy_b
        else:
            strings, ints, floats = self.show_data(data)
            msg = 'len(data) = %i\n' % ndata
            msg += 'strings  = %r\n' % strings
            msg += 'ints     = %r\n' % str(ints)
            msg += 'floats   = %r' % str(floats)
            raise NotImplementedError(msg)

        self.read_markers([-3, 1, 1])

        markers0 = self.get_nmarkers(1, rewind=False)
        record = self.read_block()

        #data = self._read_record()
        self.read_markers([-4, 1, 0, 0])
        #self._read_subtables()

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

def grids_comp_array_to_index(grids1, comps1, grids2, comps2,
                              make_matrix_symmetric):
    """maps the dofs"""
    #from pyNastran.utils.mathematics import unique2d
    ai = np.vstack([grids1, comps1]).T
    bi = np.vstack([grids2, comps2]).T
    #print('grids2 =', grids2)
    #print('comps2 =', comps2)
    from itertools import count
    #c = np.vstack([a, b])
    #assert c.shape[1] == 2, c.shape
    #urows = unique2d(c)
    #urows = c

    nid_comp_to_dof_index = {}
    j = 0
    a_keys = set()
    for nid_dof in ai:
        #nid_dof = (int(nid), int(dof))
        nid_dof = tuple(nid_dof)
        if nid_dof not in a_keys:
            a_keys.add(nid_dof)
            if nid_dof not in nid_comp_to_dof_index:
                nid_comp_to_dof_index[nid_dof] = j
                j += 1
    nja = len(a_keys)
    del a_keys

    b_keys = set()
    for nid_dof in bi:
        nid_dof = tuple(nid_dof)
        if nid_dof not in b_keys:
            b_keys.add(nid_dof)
        if nid_dof not in nid_comp_to_dof_index:
            nid_comp_to_dof_index[nid_dof] = j
            j += 1
    njb = len(b_keys)
    del b_keys


    nj = len(nid_comp_to_dof_index)
    if make_matrix_symmetric:
        ja = np.zeros(nj, dtype='int32')
        for i, nid_dof in zip(count(), ai):
            j[i] = nid_comp_to_dof_index[tuple(nid_dof)]
        for i, nid_dof in zip(count(), bi):
            j[i] = nid_comp_to_dof_index[tuple(nid_dof)]
        return j, j, nj, nj, nj
    else:
        ja = np.zeros(grids1.shape, dtype='int32')
        for i, nid_dof in zip(count(), ai.tolist()):
            ja[i] = nid_comp_to_dof_index[tuple(nid_dof)]

        jb = np.zeros(grids2.shape, dtype='int32')
        for i, nid_dof in zip(count(), bi.tolist()):
            jb[i] = nid_comp_to_dof_index[tuple(nid_dof)]

        return ja, jb, nja, njb, nj


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
