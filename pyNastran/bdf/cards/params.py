"""
defines the following card:
 - PARAM
"""
# pylint: disable=C0103,R0902,R0904,R0914
import numpy as np
from pyNastran.bdf import MAX_INT
from pyNastran.bdf.cards.base_card import BaseCard
from pyNastran.bdf.bdf_interface.bdf_card import BDFCard
from pyNastran.bdf.bdf_interface.assign_type import (
    integer, double, integer_or_blank, double_or_blank, string, string_or_blank,
    integer_double_string_or_blank, blank)
from pyNastran.bdf.field_writer_8 import print_card_8
from pyNastran.bdf.field_writer_16 import print_card_16
from pyNastran.utils.numpy_utils import integer_types, float_types
from typing import Union

#float_words_1 = [
    #b'K6ROT', b'WTMASS', b'SNORM', b'PATVER', b'MAXRATIO', b'EPSHT',
    #b'SIGMA', b'TABS']

SMALL_FIELD_PARAMS = [
    'ACOUT', 'ACOWEAK', 'ACSYM', 'ADJMETH', 'AESMAXIT', 'AESMETH', 'ADSTAT',
    'MAXLINES'] #+ INT_WORDS_1 + STR_WORDS_1


# per NX 11 QRG
temperatures = ['', '-C', '-F', '-K', 'R']
systems = [
    ('N', 'M'), ('N', 'MM'), ('MN', 'MM'), ('CN', 'CM'), ('KGF', 'M'), ('KGF', 'MM'),
    ('LBF', 'FT'), ('LBF', 'IN'), ('PDL', 'FT'),
]
unit_systems = ['']
for system in systems:
    force, length = system
    for temperature in temperatures:
        unit_system = f'{force}-{length}{temperature}'
        unit_systems.append(unit_system)
del temperatures, systems, system, temperature, force, length

PARAMS = (
    # A
    ('ACSYM', 'YES', ['YES', 'NO']),
    ('APDCON', 1.0),
    ('ADSTAT', 'YES', ['YES', 'NO']),
    ('AESDISC', 1e-8),
    ('AESMAXIT', 15),
    ('AESMETH', 'SELECT', ['SELECT', 'AUTO', 'DIRECT', 'RITZ', 'ITER']),
    ('AESRNDM', 2),
    ('AESTOL', 1e-8),
    ('AFNORM', 1.0),
    ('AFZERO', 1.0),
    ('AGGPCH', 'NO', ['YES', 'NO']),
    ('ALPHA1', (0., 0.)),
    ('ALPHA2', (0., 0.)),
    ('ALTRED', 'NO', ['YES', 'NO']),
    ('AMPCZ', 1e-6),
    ('ASCOUP', 'YES', ['YES', 'NO']),
    ('ASING', 0),
    ('AUNITS', 1.0),
    ('AUTOADJ', 'YES', ['YES', 'NO']),
    ('AUTOMPC', 'NO', ['YES', 'NO']),
    ('AUTOSEEL', 'NO', ['YES', 'NO']),
    #('AUTOSPC', 'NO', ['YES', 'NO']),  # varies by solution
    ('AUTOSPCR', 'NO', ['YES', 'NO']),
    ('AUTOSPRT', 'YES', ['YES', 'NO']),
    # B
    ('BAILOUT', 0),
    ('BETA', (1./3., 0.)),
    ('BIGER', 0.),
    ('BIGER1', 0.),
    ('BIGER2', 0.),
    ('BOLTFACT', 1.E7),
    ('BSHDAMP', 'DIFF', ['SAME', 'DIFF']),
    ('BUCKLE', -1),
    # C
    ('CA1', (1., 0.)),
    ('CA2', (1., 0.)),
    ('CB1', (1., 0.)),
    ('CB2', (1., 0.)),
    ('CDIF', 'YES', ['YES', 'NO']),
    ('CDITER', 0),
    ('CDPCH', 'NO', ['YES', 'NO']),
    ('CDPRT', 'YES', ['YES', 'NO']),
    ('CHECKOUT', 'NO', ['YES', 'NO']),
    ('CK1', (1., 0.)),
    ('CK2', (1., 0.)),
    ('CK3', (1., 0.)),
    ('CK41', (1., 0.)),
    ('CK42', (1., 0.)),
    ('CLOSE', 1.),
    ('CM1', (1., 0.)),
    ('CM2', (1., 0.)),
    # CNSTRT
    ('CNTASET', 'NO', ['YES', 'NO']),
    ('COLPHEXA', 'NO', ['YES', 'NO']),
    #('COMPMATT', 'NO', ['YES', 'NO', 'NONSMEAR']), # MSC only: 'NONSMEAR'
    ('CONFAC', 1E-5),
    ('CORROPT', 'NO', ['YES', 'NO']),
    ('COUPMASS', -1),
    ('CP1', (1., 0.)),
    ('CP2', (1., 0.)),
    ('CURV', -1),
    ('CURVPLOT', -1),
    # D
    #('DBALL', 'DBAL', ['DBAL']),
    #('DBCCONV', 'XL', ['XL']),
    ('DBCDIAG', 0),
    ('DBCOVWRT', 'YES', ['YES', 'NO']),
    ('DBDICT', -1),
    #('DBDN', 'DBAL', ['DBAL']),
    ('DBDRPRJ', 0),
    ('DBDRVER', 0),
    #('DBRCV', 'DBAL', ['DBAL']),
    #('DBUP', 'DBAL', ['DBAL']),
    ('DDRMM', 0),
    ('DESPCH', 0),
    ('DESPCH1', 6),
    ('DFREQ', 1E-5),
    ('DIGITS', 15),
    ('DOF123', 0, [0, 1]),
    ('DOPT', 0),
    ('DPEPS', 1E-4),
    ('DPREONLY', 'NO', ['YES', 'NO']),
    ('DSNOKD', 0.),
    ('DSZERO', 0.),
    ('DYNSPCF', 'NEW', ['NEW', 'OLD']),
    # E
    ('ELITASPC', 'NO', ['YES', 'NO']),
    ('EPPRT', 1E-8),
    ('EPZERO', 1E-8),
    ('ERROR', -1),
    ('EST', 2),
    ('EXTBEMI', 0, [0, 1]),
    ('EXTBEMO', 0, [0, 1]),
    ('EXTDR', 'NO', ['NO']), # missing values?
    ('EXTDROUT', 'NO', ['NO', 'MATRIXDB', 'DMIGDB', 'DMIGOP2']),  # missing values?
    ('EXTDRUNT', 31),
    ('EXTOUT', 'NO', ['NO', 'MATRIXDB', 'DMIGDB', 'DMIGOP2', 'DMIGPCH']),
    ('EXTRCV', 0),
    ('EXTUNIT', 30),
    # F
    ('F56', 'NO', ['YES', 'NO']),
    ('FACTOR', 10000),
    ('FCTC1', 1.), ('FCTC2', 1.), ('FCTC3', 1.), ('FCTC4', 1.), ('FCTC5', 1.),
    ('FCTM1', 1.), ('FCTM2', 1.), ('FCTM3', 1.), ('FCTM4', 1.), ('FCTM5', 1.),
    ('FCTZ1', 1.), ('FCTZ2', 1.), ('FCTZ3', 1.), ('FCTZ4', 1.), ('FCTZ5', 1.),
    ('FIXEDB', 0, [-2, -1, 0, 1]),
    ('FKSYMFAC', 1.),
    ('FLEXINCR', 'NO', ['YES', 'NO']),
    ('FLUIDMP', 0),
    ('FLUIDSE', 0),
    ('FOLLOWK', 'YES', ['YES', 'NO']),
    ('FRQDEPO', 'NO', ['YES', 'NO']),
    ('FRRU', 'YES', ['YES', 'NO']),
    ('FRUMIN', 500),
    ('FZERO', 1E-3),
    # G
    ('G', 0.),
    ('GFL', 0.),
    # ('Gi')
    ('GDAMPF', 'YES', ['YES', 'NO']),
    ('GEOMU', 40),
    ('GPECT', -1, [-1, 1]),
    ('GRDPNT', -1),
    ('GUSTAERO', 1, [1, -1]),
    # H
    ('HEATSTAT', 'NO', ['YES', 'NO']),
    ('HFREQ', 1E30),
    ('HFREQFL', 1E30),
    # I
    #('IFP', '')
    ('IFTM', 0, [0, 1, 2]),
    ('INEMETH', 1, [0, 1]),
    ('INP4FMT', 32, [32, 64]),
    ('INREL', 0, [-2, -1, 0]),
    ('INRELM', 0, [0, -1]),
    ('IRES', -1, [1, -1]),
    ('ITAPE', -1),
    # K
    ('K6ROT', 100.),
    ('KDAMP', 1, [-1, 1]),
    ('KDAMPFL', 1, [-1, 1]),
    ('KDIAG', -1.0),
    ('KGGCPCH', 0, [0, 1]),
    ('KGGLPCH', 0, [0, 1]),

    ('OUGCORD', '', ['', 'GLOBAL', 'BASIC']),
    ('UNITSYS', '', unit_systems),
    # L
)
string_params = {}
int_params = {}
int_params_allowed = {}
float_params = {}
float2_params = {}
for param in PARAMS:
    key = param[0]
    value = param[1]
    nparam = len(param)
    if isinstance(value, str):
        allowed = param[2]
        string_params[key] = (value, allowed)
    elif isinstance(value, int):
        if nparam == 2:
            int_params[key] = value
        else: # 3
            allowed = param[2]
            int_params_allowed[key] = (value, allowed)
    elif isinstance(value, float):
        float_params[key] = value
    else: # two values
        nvalues = len(value)
        assert nvalues == 2, param
        if isinstance(value[0], float):
            float2_params[key] = value
        else:
            raise RuntimeError(param)
del allowed

STR_WORDS_1 = {
    'POSTEXT', 'PRTMAXIM', 'AUTOSPC', 'OGEOM', 'PRGPST',
    'RESVEC', 'RESVINER', 'OGPS', 'OIBULK', 'OMACHPR',
    'UNITSYS', 'F56', 'OUGCORD', 'OGEM', 'EXTSEOUT',

    'AUTOSPC', 'CDIF', 'EXTDROUT', 'OGEOM', 'OMID', 'PRTMAXIM', 'PRGPST',
    'POSTEXT', 'RESVEC', 'SUPAERO', 'SHLDAMP', 'ZROCMAS',
    'ALTRED', 'RSCON', 'RESVINER', 'MESH', 'SKINOUT', 'VUPENTA',
    'DYNSEN', 'PRTGPL', 'PRTEQXIN', 'PRTGPDT', 'PRTCSTM', 'PRTBGPDT', 'PRTGPTT',
    'PRTMGG', 'PRTPG', 'FOLLOWK', 'DBCCONV', 'COMPMATT', 'ADB', 'HEATSTAT',
    'AUTOSPCR', 'CHECKOUT', 'EPSILONT', 'RMS', 'RESPATH', 'AUTOMSET', 'COMPMATT',
    'F56', 'METHFL', 'DMIGNRG', 'METHT', 'AEDB', 'SAVEOFP', 'DBDN', 'DBCONV',
    'ASCOUP', 'DBUP', 'VUHEXA', 'FLEXINCR', 'SRCOMPS', 'VUTETRA', 'DYNSPCF',
    'SRCOMPS', 'SERST', 'PGRPST', 'CDPRT', 'CFDIAGP', 'CWDIAGP', 'ENFMOTN', 'RMSINT',
    'DEBUG', 'CHKOUT', 'DBALL', 'OIBULK', 'RMXTRAN', 'CDPCH', 'AUTOADJ', 'SEMAP',
    'SOFTEXIT', 'SM', 'RESFLEX', 'SENSUOO', 'OGEM', 'XBYMODE', 'STRESS', 'OPTION',
    'AOTOSPC', 'ARBMAS', 'ARBMSS', 'OMACHPR', 'AUTOMPC', 'BSHDAMP', 'ASCII4',
    'SPCSTR', 'AUTOSPRT', 'PBRPROP', 'ELITASPC', 'SECOMB', 'CNTASET', 'WMODAL',
    'RESVALT', 'PRTRESLT', 'ITFPRNT', 'AMLS', 'PRGPOST', 'SDAMPUP', 'COLPHEXA',
    'ELEMITER', 'ROTSYNC', 'MECHFIX', 'CTYPE', 'SESDAMP', 'SYNCDAMP',
}
INT_WORDS_1 = {
    'ACTH', 'ADMPOST', 'ALPHA', 'ALTSHAPE', 'AMGOK', 'APPEN', 'AROWS',
    'ASING', 'BAILOUT', 'BAND', 'BUCKLE', 'CDITER', 'CHKMASS',
    'CHKSTIF', 'CNSTRT', 'COMPARE', 'COUPMASS', 'CURV', 'CURVPLOT',
    'CYCIO', 'CYCSEQ', 'DBCDIAG', 'DBCLEAN', 'DBDICT', 'DBDRNL',
    'DBDRNL', 'DBDROPT', 'DBNBLKS', 'DBRNL', 'DBSORT', 'DCOMP',
    'DDRMM', 'DEBUGPRT', 'DESPCH', 'DESPCH1', 'DESPCHL', 'DIA600',
    'DLOAD', 'DPHFLG', 'DRMH', 'DUMMY', 'DYENDTIM', 'DYLDKND',
    'DYMATS1', 'DYNAMES', 'EFWGT', 'EIGVECDS', 'EIGVECDS', 'ERROR',
    'EST', 'EXTBEMI', 'EXTBEMO', 'FIXEDB', 'FLUIDMP', 'FLUIDSE',
    'FLYLOAD', 'GEOMPLT', 'GEOMU', 'GPECT', 'GPFORCE', 'GPWG', 'GRDEQ',
    'GRDPNT', 'GUSTAERO', 'GYROAVG', 'GYROFATL', 'I', 'I', 'IEXT',
    'IFTM', 'IFTM', 'INFILE', 'INREL', 'INRLM', 'INT1', 'INT2', 'INT3',
    'INT4', 'IPAD', 'IPU', 'IRES', 'ITSMAX', 'IUNIT4', 'J', 'KDAMP',
    'KDAMPFL', 'KEPRT', 'KGGCPCH', 'LAMA', 'LANGLE', 'LGDIS', 'LGDISP',
    'LGSTRN', 'LMODES', 'LMODESFL', 'LOADU', 'LOOPID', 'LSTRN',
    'MAKEMAT', 'MARCASUM', 'MARCAUCU', 'MARCAUDM', 'MARCAUNW',
    'MARCAUST', 'MARCAUTN', 'MARCAUTO', 'MARCAUTR', 'MARCCENT',
    'MARCEKND', 'MARCGAUS', 'MARCNMPC', 'MARCONLY', 'MARCOPT',
    'MARCPOS', 'MARCPOST', 'MARCPRN', 'MARCPTH', 'MARCRBE2', 'MARCRBE3',
    'MARCREVR', 'MARCRIGD', 'MARCRUN', 'MARCSETT', 'MARCSINC',
    'MARCSLHT', 'MARCSUMY', 'MARCT19', 'MARCTABL', 'MARCTNSF',
    'MARCTNSF', 'MARCTOL', 'MARCVERS', 'MATFILE', 'MATNL', 'MAXIT',
    'MAXITER', 'MAXLINES', 'MESHG', 'METHCMRS', 'MODACC', 'MODTRK',
    'MPCX', 'MPTDUMP', 'MRALIAS', 'MRFOLOW1', 'MRFOLOW3', 'MRFOLOW4',
    'MRORINTS', 'MROUTLAY', 'MRTIMING', 'NASPRT', 'NBRUPT', 'NEWSEQ',
    'NEWSET', 'NLAYERS', 'NLDISP', 'NLPACK', 'NLTOL', 'NMLOOP', 'NOAP',
    'NOBAND', 'NOCOMP', 'NOCOMPS', 'NODATA', 'NOELOF', 'NOELOP',
    'NOFISR', 'NOGPF', 'NOINT', 'NONCUP', 'NOPRT', 'NOREAL', 'NQSET',
    'NT', 'NTOL', 'NUMOUT', 'NUMOUT1', 'OCORD', 'OELMSET', 'OGRAV',
    'OGRDOPT', 'OGRDSET', 'OLDSEQ', 'OMAXR', 'OMSGLVL', 'OP2FMT',
    'OPCHSET', 'OPGTKG', 'OPPHIPA', 'OPTEXIT', 'ORDER', 'ORIGID',
    'OSWELM', 'OSWPPT', 'OSWPTT', 'OUNIT2', 'OUNIT2M', 'P1', 'PANELMP',
    'PARAM', 'PARTFAC', 'PDRMSG', 'PEDGEP', 'PH2OUT', 'PIVOT', 'PLTMSG',
    'PNCHDB', 'PNCHDB', 'POSTU', 'PROUT', 'PVALDR', 'PVALINIT', 'QRES',
    'R', 'RANREAL', 'RBMODES', 'RFORCE', 'RMXCRT', 'RMXTRN', 'ROTCSV',
    'ROTGPF', 'RSOPT', 'RSPECTRA', 'S1', 'S1A', 'S1AG', 'S1AM', 'S1G',
    'S1M', 'SCRSPEC', 'SELOCAL', 'SEMAPPRT', 'SEOP2CV', 'SEP1XOVR',
    'SEQOUT', 'SESEF', 'SIGNB', 'SIZE', 'SMALL', 'SNORMPRT', 'SOLADJC',
    'STRUCTMP', 'SUBID', 'SUBSKP', 'SUPER', 'TDAMP', 'TEMPMATE',
    'TESTNEG', 'TESTNO', 'TORSIN', 'TSTATIC', 'UNUSED', 'USETPRT',
    'USETSEL', 'VMOPT', 'XFLAG', 'XYUNIT', 'OPPHIB',
}
FLOAT_PARAMS = {
    'AUNITS', 'CHIMERA', 'EPSHT', 'EPPRT', 'G', 'GFL', 'K6ROT', 'MAXRATIO',
    'PATVER', 'PRPHIVZ', 'SIGMA', 'SNORM', 'TABS', 'TINY', 'W3', 'W4', 'WTMASS',
    'HFREQ', 'CONFAC', 'VREF', 'CFRANDEL', 'DSNOKD', 'DFREQ', 'TESTSE', # 'TESTNEG',
    'KDIAG', 'NDAMP', 'M', 'Q', 'TABST', 'TOLRSC', 'DPEPS', 'PREFDB', 'MACH',
    'LFREQ', 'EPZER0', 'W4FL', 'MARCAUT0', 'DSZERO', 'BETAXX', 'STIME',
    'CWRANDEL', 'MARCAUTT', 'AMORT', 'SEAU', 'TOLLRSC', 'EXTRCV', 'PRPA', 'LFREQFL',
    'MARCAUTS', 'MARCFRIC', 'SOFFSET', 'LMFACT', 'PENFN', 'CLOSE', 'BIGER', 'WR3',
    'HFREQFL', 'EPZERO', 'PRPJ', 'FRIC', 'LAMBDAS', 'MARCSCL', 'MARCAUTB',
    'PATVERS', 'RPM', 'W3FL', 'ROTOMRF', 'MARCAUTM', 'TOL', 'MXLAGM1', 'MARCTVL',
    'BETA', 'FKSYMFAC', 'BOLTFACT', 'ROTCMRF', 'LMSCAL', 'RESVRAT', 'SMALLQ',
    'LAMLIM', 'MARCAUTX', 'MARCRSCL', 'BIGER1', 'DYCOWPRD', 'DYDTOUT', 'FZERO',
    'BIGER2', 'DYCOWPRP', 'MTRFMAX', 'REAL1', 'REAL2', 'REAL3',
}

FLOAT2_PARAMS = {'BETA1', 'CA1', 'CA2'}

STR_WORDS_1.update(set(string_params.keys()))
INT_WORDS_1.update(set(int_params.keys()))
FLOAT_PARAMS.update(set(float_params.keys()))
FLOAT2_PARAMS.update(set(float2_params.keys()))

INT_STR_WORDS_1 = INT_WORDS_1 | STR_WORDS_1

class PARAM(BaseCard):
    type = 'PARAM'
    _field_map = {1: 'key'}

    def _update_field_helper(self, n, value):
        if n - 2 >= 0:
            try:
                self.values[n - 2] = value
            except IndexError:
                msg = 'Field %r=%r is an invalid %s entry for key=%r.' % (
                    n, value, self.type, self.key.upper())
                raise IndexError(msg)
        else:
            raise KeyError('Field %r=%r is an invalid %s entry.' % (n, value, self.type))

    @classmethod
    def _init_from_empty(cls):
        key = 'POST'
        values = -1
        return PARAM(key, values, comment='')

    def __init__(self, key, values, comment=''):
        """
        Creates a PARAM card

        Parameters
        ----------
        key : str
            the name of the PARAM
        values : int/float/str/list
            varies depending on the type of PARAM
        comment : str; default=''
            a comment for the card

        """
        if comment:
            self.comment = comment
        self.key = key
        if isinstance(values, list):
            pass
        elif isinstance(values, (integer_types, float_types, str)):
            values = [values]
        self.values = values
        if isinstance(self.values, tuple) or isinstance(self.values[0], (list, tuple)):
            raise TypeError((key, self.values))
        #assert not isinstance(values, tuple), values
        #if isinstance(values, list):
            #assert not isinstance(values[0], tuple), values

    @classmethod
    def add_card(cls, card, comment=''):
        """
        Adds a PARAM card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card

        """
        key = string(card, 1, 'key')

        n = 1
        value = None
        if key == 'ACOUT':
            value = string_or_blank(card, 2, 'value', default='PEAK')
        elif key == 'ACOWEAK':
            value = string_or_blank(card, 2, 'value', default='NO')
        elif key == 'ADJMETH':
            value = integer_or_blank(card, 2, 'value', default=0)
        elif key == 'ADPCON':
            value = double_or_blank(card, 2, 'value', default=1.0)
        #elif key == 'ADMPOST':
            #value = string_or_blank(card, 2, 'value', default=0) ## TODO: 0 is not a string
        elif key == 'ADSDISC':
            value = double_or_blank(card, 2, 'value', default=1e-8)
        elif key == 'AESMAXIT':
            value = integer_or_blank(card, 2, 'value', default=15)
        elif key == 'AESMETH':
            value = string_or_blank(card, 2, 'value', 'SELECT')
            assert value in ['SELECT', 'AUTO', 'DIRECT', 'RITZ', 'ITER'], 'value=%s' % value
        elif key == 'AESTOL':
            value = double_or_blank(card, 2, 'value', default=1e-10)
        elif key in ['ALPHA1FL', 'ALPHA2FL']:  # check alpha1/alpha1FL
            value1 = double_or_blank(card, 2, 'value1', default=0.0)
            value2 = double_or_blank(card, 3, 'value2', default=0.0)
            n = 2
        elif key == 'COMPMATT':
            #('COMPMATT', 'NO', ['YES', 'NO', 'NONSMEAR']), # MSC only: 'NONSMEAR'
            value = string_or_blank(card, 2, 'value1', 'NO')
            if value == 'NONS':  # assume
                value = 'NONSMEAR'
            if value in {'SMEAR', 'SMEARED'}:  # assume
                value = 'YES'
            assert value in {'YES', 'NO', 'NONSMEAR'}, 'value=%r' % value

        elif key == 'POST':
            value = integer_or_blank(card, 2, 'value', default=1)
        elif key == 'UNITSYS':
            value = string_or_blank(card, 2, 'value', default='')
        elif key == 'MDOF':
            # NX  - integer
            # MSC - string
            value = integer_string_or_blank(card, 2, 'value', default=None)

        #-------------------------------------------------------------
        # strings; has defaults
        elif key in string_params:
            default, allowed_values = string_params[key]
            value = string_or_blank(card, 2, 'value', default=default)
            assert value in allowed_values, f'value={value} allowed={allowed_values}'

        # ints; has defaults
        elif key in int_params:
            default=int_params[key]
            value = integer_or_blank(card, 2, 'value', default=default)
        elif key in int_params_allowed:
            default, allowed_values = int_params_allowed[key]
            value = integer_or_blank(card, 2, 'value', default=default)
            assert value in allowed_values, f'value={value} allowed={allowed_values}'

        # floats; has defaults
        elif key in float_params:
            default=float_params[key]
            value = double_or_blank(card, 2, 'value', default=default)
        elif key in float2_params:
            defaults = float2_params[key]
            value = double_or_blank(card, 2, 'value', default=defaults[0])
            value = double_or_blank(card, 2, 'value', default=defaults[1])
            n = 2

        # unchecked catch all
        elif key in STR_WORDS_1:
            value = string(card, 2, 'value')
        elif key in INT_WORDS_1:
            value = integer(card, 2, 'value')
        elif key in FLOAT_PARAMS:
            value = double(card, 2, 'value')
        elif key in FLOAT2_PARAMS:
            value1 = double(card, 2, 'value1')
            value2 = double(card, 3, 'value2')
            values = [value1, value2]
            n = 2

        #-------------------------------------------------------------
        else:
            #raise NotImplementedError(card)
            n = 2
            value1 = integer_double_string_or_blank(card, 2, 'value1')
            value2 = integer_double_string_or_blank(card, 3, 'value2')
            if value2 is None:
                value = value1
                n = 1

        if value is None:
            # n=2 or blank
            if isinstance(value1, str):
                assert ' ' not in value1, f'PARAM value1={value1!r}'
            if isinstance(value2, str):
                assert ' ' not in value2, f'PARAM value2={value2!r}'
            values = [value1, value2]
        else:
            # n=1
            if isinstance(value, str):
                assert ' ' not in value, f'PARAM value={value!r}'
            values = [value]

        if n == 1:
            assert len(card) <= 3, f'len(PARAM card)={len(card):d} card={card!r}'
        else:
            assert len(card) <= 4, f'len(PARAM card)={len(card):d} card={card!r}'
        return PARAM(key, values, comment=comment)

    def update_values(self, value1=None, value2=None):
        """
        Updates value1 and value2.  Performs type checking based on the PARAM
        type after setting any default value(s).

        Parameters
        ----------
        value1 : varies; default=None
            the main value
        value2 : varies; default=None
            optional value

        If you want to access the data directly, use:
        >>>  param = bdf.params['POST']
        >>> param.values[0] = -1  # value1
        >>> param.values[1] = 3   # value2
        >>>

        .. note::  Most PARAM cards only have one value.  Some have two.

        """
        if self.key == 'ACOUT':
            if value1 is None:
                value1 = 'PEAK'
            if not isinstance(value1, str):
                msg = f'key={self.key} value1={value1!r} must be a string.'
                raise TypeError(msg)

        elif self.key == 'ACOWEAK':
            if value1 is None:
                value1 = 'NO'
            if not isinstance(value1, str):
                msg = f'key={self.key} value1={value1!r} must be a string.'
                raise TypeError(msg)

        elif self.key == 'ACSYM':
            if value1 is None:
                value1 = 'YES'
            if not isinstance(value1, str):
                msg = f'key={self.key} value1={value1!r} must be a string.'
                raise TypeError(msg)

        elif self.key == 'ADJMETH':
            if value1 is None:
                value1 = 0
            if not isinstance(value1, int):
                msg = f'key={self.key} value1={value1!r} must be an integer.'
                raise TypeError(msg)

        #elif self.key == 'ADMPOST': ## TODO: 0 is not a string
            #value = string_or_blank(card, 2, 'value', 0)

        elif self.key == 'ADSTAT':
            if value1 is None:
                value1 = 'YES'
            if not isinstance(value1, str):
                msg = f'key={self.key} value1={value1!r} must be a string.'
                raise TypeError(msg)

        elif self.key in ['ALPHA1', 'ALPHA2', 'ALPHA1FL', 'ALPHA2FL']:
            if value1 is None:
                value1 = 0.0
            if value2 is None:
                value2 = 0.0
            if not isinstance(value1, float):
                msg = f'key={self.key} value1={value1!r} must be a float.'
                raise TypeError(msg)
            if isinstance(value2, float):
                msg = f'key={self.key} value2={value2!r} must be a float.'
                raise TypeError(msg)

        elif self.key in ['CB1', 'CB2', 'CK1', 'CK2', 'CK3', 'CM1', 'CM2', 'CP1', 'CP2']:
            if value1 is None:
                value1 = 1.0
            if value2 is None:
                value2 = 0.0
            if not isinstance(value1, float):
                msg = f'key={self.key} value1={value1!r} must be a float.'
                raise TypeError(msg)
            if isinstance(value2, float):
                msg = f'key={self.key} value2={value2!r} must be a float.'
                raise TypeError(msg)

        else:
            if not isinstance(value1, (int, float, str)):
                msg = f'key={self.key} value1={value1!r} must be an integer, float, or string.'
                raise TypeError(msg)

        self.values = [value1]
        if value2 is not None:
            self.values.append(value2)

    def raw_fields(self):
        list_fields = ['PARAM', self.key] + self.values
        return list_fields

    def repr_fields(self):
        return self.raw_fields()

    def write_card(self, size: int=8, is_double: bool=False) -> str:
        card = self.raw_fields()
        if isinstance(self.values[0], int) and self.values[0] > MAX_INT:
            return self.comment + print_card_16(['PARAM', self.key, self.values[0]])

        if self.key in INT_STR_WORDS_1:
            return '%sPARAM   %8s%8s\n' % (
                self.comment, self.key, self.values[0])

        if size == 8 or self.key in SMALL_FIELD_PARAMS:
            return self.comment + print_card_8(card)
        return self.comment + print_card_16(card)


class PARAM_NASA95(BaseCard):
    @classmethod
    def _init_from_empty(cls):
        key = 'POST'
        values = -1
        return PARAM_NASA95(key, values, comment='')

    def __init__(self, key, values, comment=''):
        """
        Creates a PARAM card

        Parameters
        ----------
        key : str
            the name of the PARAM
        values : int/float/str/list
            varies depending on the type of PARAM
        comment : str; default=''
            a comment for the card

        """
        if comment:
            self.comment = comment
        self.key = key
        if isinstance(values, (int, float, str)):
            values = [values]
        self.values = values

    @classmethod
    def add_card(cls, card, comment=''):
        """
        Adds a PARAM card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card

        """
        key = string(card, 1, 'key')

        if key == 'AUTOSPC':
            values = [
                integer_or_blank(card, 2, 'value', 0)
            ]
        elif key == 'BETA':
            values = [
                double_or_blank(card, 2, 'value', 0.55)
            ]
        elif key == 'BETAD':
            values = [
                integer_or_blank(card, 2, 'value', 4)
            ]
        elif key == 'COUPMASS':
            values = [
                integer_or_blank(card, 2, 'value', 1)
            ]
        # CTYPE
        elif key == 'CYCIO':
            values = [
                integer_or_blank(card, 2, 'value', 1)
            ]
        elif key == 'CYCSEQ':
            values = [
                integer_or_blank(card, 2, 'value', -1)
            ]
        elif key == 'EPSHT':
            values = [
                double_or_blank(card, 2, 'value', 0.001)
            ]
        elif key == 'ESPIO':
            values = [
                double_or_blank(card, 2, 'value', 1e-5)
            ]
        elif key == 'G':
            values = [
                double(card, 2, 'G')
            ]
        else:
            ifield = 2
            #                 x
            #[PARAM, AUTOSPC, 2]
            values = []
            while ifield < len(card):
                value = integer_double_string_or_blank(card, ifield, 'field_{i}')
                values.append(value)
                ifield += 1
        return PARAM_NASA95(key, values, comment=comment)

    def raw_fields(self):
        list_fields = ['PARAM', self.key] + self.values
        return list_fields

    def repr_fields(self):
        return self.raw_fields()

    def write_card(self, size: int=8, is_double: bool=False) -> str:
        card = self.raw_fields()
        if size == 8:
            return self.comment + print_card_8(card)
        return self.comment + print_card_16(card)

class PARAM_MYSTRAN(BaseCard):
    """
    Parameter Name DataType Function of Parameter

    ARP_TOL; real; Default=1x10^-6
        Tolerance to use in Lanczos eigenvalue extraction method for convergence
    ART_KED
        (for diff stiffness – not fully implemented)
        char; field 3: ART_KED, default=N. If Y add artificial stiff to diag of KED stiff matrix
        char; field 4: ART_TRAN_MASS: value for translation degrees of freedom, default 1x10^-6
        char; field 5: ART_ROT_MASS: value for translation degrees of freedom, default 1x10^-6
    ART_MASS
        char; field 3: ART_MASS, default=N. If Y add artificial mass to diag of MGG mass matrix
        real; Field 4: ART_TRAN_MASS: value for translation degrees of freedom, default 1x10^-6
        real; Field 5: ART_ROT_MASS:  value for translation degrees of freedom, default 1x10^-6

    AUTOSPC
        ['PARAM', 'AUTOSPC', 'Y', '1.0E-9', None, 'Y', 'Y']
        Char  Field 3: AUTOSPC value, default=Y (AUTOSPC), N turns AUTOSPC off.
        Real  Field 4: AUTOSPC_RAT, default=1x10^-6 (see Section 3.4.1.1)
        Int   Field 5: AUTOSPC_NSET, default=1 (see Section 3.4.1.1)
        Char  Field 6: AUTOSPC_INFO, default=N. If Y then print messages about the AUTOSPC's
        Char  Field 7: AUTOSPC_SPCF, default=N. If Y print AUTOSPC forces of constraint
    BAILOUT Int; Default=1
        If > 0 quit if a singularity in decomposing a matrix is detected.
        If <= 0 do not quit
    CBMIN3 Real; Default=2.0
        The constant CB used in tuning the shear correction factor in
        Ref 3 for the TRIA3 plate element.
    CBMIN4 real; default=3.6
       The constant CB used in tuning the shear correction factor in
       Ref 4 for the QUAD4 plate element (QUAD4TYP='MIN4'). See Ref 4
    CBMIN4T real; default=3.6
       The constant CB used in tuning the shear correction factor in
       Ref 4 for the QUAD4 plate element (QUAD4TYP='MIN4T').
    CHKGRDS
    char; default=Y.
        If N do not check that all grids for all elements exist
    CUSERIN
       char; If this parameter is present, Bulk Data entries for Craig-Bampton (CB)
             reduced models will be written to the F06 file as a CUSERIN element
             (including grids, coordinate systems, etc)
       int;  field 3: element ID, default=9999999
       int;  field 4: property ID, default=9999999
       int;  field 5: starting index for the SPOINT's to represent modes of the CB model, default=1001
       int;  field 6: IN4 file number that goes on the PUSERIN entry for this CUSERIN element, default=9999999
       char; field 7: Set-ID for the CUSERIN element (typically the 'R', or boundary, set), default is blank field
       int;  field 8: Format for how to write the component numbers
                     (1 thru 6) for each grid of the CUSERIN element.
                     If 0, write them in compact form (e.g. 1356). If > 0
                     write them in expanded form (1 3 56), default=0
    DARPACK int; default=2
        how many extra modes to find above EIG_N2 on the EIGRL entry.
        These few highest mode are not used due to difficulty with
        getting good GP force balance.
    EIGESTL int; default=5000
        For eigenvalue problems by the Lanczos method, if the number of
        L-set DOF's exceed EIGESTL the method for specifying the search
        range will be changed from F1 to F2 to N (see EIGRL Bulk Data
        entry) to avoid excessive run times (since the code to estimate
        the number of eigens in the F1 to F2 range can be excessive).
    EIGNORM2 char; default=N
        If 'Y' then eigenvectors will be renormalized a last time by
        multiplying by a set of scale factors (1 per eigenvector)
        supplied in a file with the same name as the input file and
        extension 'EIN' (if it exists)
    ELFORCEN; char; default=GLOBAL
        If and nodal forces have been requested in Case Control and:
          ELFORCEN=GLOBAL, they will be output in the global coordinate system.
          ELFORCEN=BASIC, they will be output in the basic coordinate systeml.
          ELFORCEN=LOCAL, they will be output in the local element coordinate system.
    EPSERR char; default=Y.
        If N, do not calculate the NASTRAN like 'epsilon error estimate'
    EPSIL Real
        There are 3 EPSIL(i) values each of which requires a separate
        PAPAM EPSIL Bulk Data entry with the index (i) in field 3 and
        EPSIL(i) value in field 4.  These are small numbers used in
        MYSTRAN for the purposes indicated below:
          1) EPSIL(1) (default=1x10^-15) is used in MYSTRAN such that, in
                      any real number comparisons, any real number whose
                      absolute magnitude is less than EPSIL(1) is
                      considered to be zero. If no PARAM EPSIL 1 entry is
                      in the data file then this value is reset (from the
                      default) in LINK1 to a value based on machine
                      precision calculated using LAPACK BLAS function
                      DLAMCH. If the user has a PARAM EPSIL 1 entry, this
                      value will be used for EPSIL(1) instead of the
                      LAPACK machine precision.
          2) Currently not used
          3) EPSIL(3) is used in the Inv erse Power method of eigenvalue
                      extraction to test convergence of an eigenvalue. The
                      default value (% change) is 1x10^-5 %
          4) EPSIL(4) is used to calculate the maximum warp for quadrilateral
                      plate elements, above which a warning message will be
                      written. This maximum warp is EPSIL(2) times the
                      average length of the quadrilateral's two diagonals.
                      The default for EPSIL(2) is 1.x10^-1.
          5) EPSIL(5) (default 1.x10^-6) is used in BAR and ROD margin of
                      safety calculations. If a stress magnitude is less
                      than EPSIL(5) a 1.x10^10 margin of safety will
                      printed out for that stress (in other words, an
                      infinite margin of safety)
          6) EPSIL(6) (default 1.x10^-15) is used in BAR margin of safety calculations
    EQCHECK
        int;  field 3: default=0 (basic origin) or reference grid to use
                       in calculating the rigid body displacement matrix
                       for the equilibrium check
        int;  field 4: If nonzero, do equilibrium check on the G-set
        int;  field 5: If nonzero, do equilibrium check on the N-set
        int;  field 6: If nonzero, do equilibrium check on the F-set
        int;  field 7: If nonzero, do equilibrium check on the A-set
        int;  field 8: If nonzero, do equilibrium check on the L-set
        real; field 9: EQCHK_TINY, default=1x10^-5. I Do not print grid forces smaller than this
        char; field 10: default=N. If Y, normalize the grid forces on diagonal stiffness

        The value in fields 4-8 can be:
        1: print loads due to rigid body displacements
        2: print strain energy due to rigid body displacements
        3: print both
    GRDPNT Int default=-1.
       If not -1 then the value is interpreted as a grid number
       If GRDPNT != 0, calculate total mass properties of the model
       relative to the basic coordinate system origin or relative
       to the specified grid.
    GRIDSEQ
        char; field 3: GRIDSEQ value (default=BANDIT). Other values are
                       GRID and INPUT. BANDIT is automatic grid sequencing.
                       GRID is sequencing in grid ID numerical order.
                       INPUT is sequencing in the grid input order.
        char; field 4: SEQQUIT, default=N. If Y, then quit in the sequence
                       processor if BANDIT did not run correctly.
        char; field 5: SEQPRT, default=N. If Y, print SEQGP card images
                       generated by BANDIT to the F06 output file
    HEXAXIS; Char
        'SIDE12', use side 1-2 as the local elem x axis.
        'SPLITD' (default), use angle that splits the 2 diags to define the elem x axis
    IORQ1M; Int; default=2
        Gaussian integration order for membrane direct stress terms for
        the QUAD4, QUAD4K quadrilateral elements
    IORQ1S; Int; default=1
        Gaussian integration order for membrane shear stress terms for all quad elements
    IORQ1B; Int; default=2
        Gaussian integration order for bending stress terms for the QUAD4K element
    IORQ2B; Int; default=2
        Gaussian integration order for bending stress terms for the QUAD4 element
    IORQ2T; Int; default=3
        Gaussian integration order for transverse shear stress terms for the QUAD4 element
    ITMAX; Int; default=5
        Max number of iterations in refining the solution when parameter UREFINE=Y
    KLLRAT; Char; default=Y
        to tell whether to calc ratio of max/min KLL diagonal terms
    KOORAT; Char; default=Y
        to tell whether to calc ratio of max/min KOO diagonal terms
    LANCMETH; Char
        Procedure to use for Lanczos eigenvalue extraction (ARPACK or TRLan)
    MATSPARS; Char
        If=Y (default), use sparse matrix routines for add/multiply in all matrix operations.
        If N, use full matrix add/multiply (not recommended)
    MAXRATIO; Real; Default=1X10^7
        Max value of matrix diagonal to factor diagonal before messages
        are written and BAILOUT tested for aborting run
    MEFMCORD; Int; default=0
        The coordinate system in which to calculate modal mass and participation factors
    MEFMLOC; Char
        Reference location for calculating modal effective mass in Craig-Bampton (SOL 31) analyses.
        This only affects the rotational modal effective masses.
        Field 3 can be GRDPNT, GRID or CG:
           If field 3=GRDPNT (default): ref point is the same as the one for PARAM GRDPNT
           If field 3=CG: use the model center of gravity as the reference point
           If field 3=GRID: use the grid point number in field 4 as the reference point
        Field 4: MEFMGRID (grid to use when field 3 is GRID)
    MEMAFAC; Int; default=0.9
        Factor to multiply the size request of memory to be allocated when
        looping to find an allowable amount of memory to allocate. Used
        when the initial request for memory (in subrs ESP or EMP) cannot
        be met and we know that the request is conservative.
    MIN4TRED; Char; default=STC
        Defines the method for how the 5th node of the MIN4T element is
        reduced out (to get a 4 node quad element). STC (default) is static
        condensation. B%$ (not implemented as of Version 3.0) uses a method
        developed by the element author (see Users Reference manual)
    MKLFACij; Char; default=INDEF
        Matrix type for use in decomposing matrices in various subroutines
        in MYSTRAN when PARAM SOLLIB is IntMKL'
        MKLFAC21 is for use in subr REDUCE_KAA_TO_KFF
        MKLFAC31 is for use in subr LINK3
        MKLFAC41 is for use in subr EIG_INV_PWR
        MKLFAC42 is for use in subr EIG_LANCZOS_ARPACK
        MKLFAC61 is for use in subr CALC_KRRcb
        MKLFAC62 is for use in subr SOLVE_DLR
        MKLFAC63 is for use in subr SOLVE_PHIZL1
    MKLMATST; Char; default=NONSYM
        Matrix structure to use when PARAM SOLLIB=IntMKL.
        Values can be NONSYM or SYM
    MKLSTATS; char; default=N
        If Y write stats on matrix decomposition when PARAM SOLLIB=IntMKL
    MPFOUT; char
        (1) '6' (default) indicates to output modal participation
                factors (MPF) relative to the 6 DOF's at grid MEFMGRID (see PARAM MEFMLOC)
        (2) 'R' indicates to output MPF's for all of the R-set DOF's individually
    MXALLOCA; int; default=10
        Max number of attempts to allow when trying to allocate memory in subroutine ALLOCATE_STF_ARRAYS
    MXITERI; int; default=50
        Max number of iterations to use in the Inverse Power eigenvalue extraction method
    MXITERL; int; default=50
        Max number of iterations to use in the Lanczos eigenvalue extraction method
    OTMSKIP; int
        Number of lines to skip between segments of OTM text file descriptors
    PBARLDEC; int; default=5
        Number of decimal digits when writing PBAR equivalents for PBARL entry real data
    PBARLSHR; char; default=Y
        Include K1, K2 for PBAR equiv to PBARL BAR properties
    PCHSPC1
        char; Field 3: PCHSPC1 value (default=N, do not punch SPC1 card
                       images for constraints generated by the AUTOSPC
                       feature, use Y to punch these)
        int;  Field 4: SPC1SID value (default=9999999, the set ID to put
                       on the SPC1 card images)
        char; Field 5: SPC1QUIT value (default=N, do not stop after SPC's
                       are punched, or Y to stop processing)
    PCMPTSTM
        Real; Factor to multiply composite ply thickness for effective shear thickness
    PCOMPEQ; int; default=0
        Indicator to write equiv PSHELL, MAT2 to F06 for PCOMP's.
        If > 0, write the equivalent PSHELL amd MAT2 Bulk Data entries
        for the PCOMP. If > 1 also write the data in a format with a
        greater number of digits of accuracy.
    POST; int
        If=-1 then write FEMAP neutral file for post processing of MYSTRAN outputs
    PRTBASIC; int
        If=1 print grid coordinates in the basic coordinate system
    PRTCGLTM; int
        If=1 print CB matrix for C.G. LTM loads
    PRTCONN; int
        If=1, print table of elements connected to each grid
    PRTCORD; int
        If PRTCORD=1 print coordinate system transformation data
    PRTDISP; int
        PRTDISP(I), I=1-4 go in fields 3-6 of the PARAM PRTDISP entry
        that prints displacement matrices for G, N, F, and/or A-sets:
          V1=PRTDISP(1)=1 print UG
          V2=PRTDISP(2)=1 or 3 print UN, 2 or 3 print UM
          V3=PRTDISP(3)=1 or 3 print UF, 2 or 3 print US
          V4=PRTDISP(4)=1 or 3 print UA, 2 or 3 print UO
          V5=PRTDISP(5)=1 or 3 print UL, 2 or 3 print UR
    PRTDLR; int
        If=1, the DLR matrix will be printed
    PRTDOF; int
        If PRTDOF=1 or 3 print TDOF table, in grid point ID numerical
        order, which gives a list of the degree of freedom numbers for
        each displacement set (size is number of degrees of freedom x number
        of displacement sets).
        If PRTDOF=2 or 3 print TDOF table, in degree of freedom numerical
        order, which gives a list of the degree of freedom numbers for each
        displacement set (size is number of degrees of freedom x number
        of displacement sets)
    PRTFOR; int
        PRTFOR(I), I=1-4 go in fields 3-6 of the PARAM PRTFOR entry that
        prints sparse force matrices for G, N, F, and/or A-sets:
          V1=PRTFOR(1)=1 print sparse PG
          V2=PRTFOR(2)=1 or 3 print sparse PN, 2 or 3 print PM
          V3=PRTFOR(3)=1 or 3 print sparse PF, 2 or 3 print PS
          V4=PRTFOR(4)=1 or 3 print sparse PA, 2 or 3 print PO
          V5=PRTFOR(5)=1 or 3 print sparse PL, 2 or 3 print PR
    PRTGMN; int
        If PRTGMN=1, print GMN matrix
    PRTGOA; int
        If PRTGOA=1, print GOA matrix

    PRTHMN; int
        If=1 print HMN constraint matrix
    PRTIFLTM; int
        If=1 print CB matrix for Interface Forces LTM
    PRTKXX; int
        If=1 print CB matrix KXX
    PRTMASSD; int
        Same as PRTMASS, except only print diagonal terms
    PRTMASS; int
        PRTMASS(I), I=1-4 go in fields 3-6 of the PARAM PRTMASS entry that
        prints sparse mass matrices for G, N, F, and/or A-sets:
          V1=PRTMASS(1)=1 print sparse MGG
          V2=PRTMASS(2)=1 or 3 print sparse MNN, 2 or 3 print MNM, MMM
          V3=PRTMASS(3)=1 or 3 print sparse MFF, 2 or 3 print MFS, MSS
          V4=PRTMASS(4)=1 or 3 print sparse MAA, 2 or 3 print MAO, MOO
          V5=PRTMASS(5)=1 or 3 print sparse MLL, 2 or 3 print MLR, MRR
    PRTMXX; int
        If=1 print CB matrix MXX
    PRTOU4; int
        If > 0 write all OU4 (OUTPUT4) matrices to F06 file
    PRTPHIXA; int
        If=1 print CB matrix PHIXA
    PRTPHIZL; int
        If=1 print CB matrix PHIZL
    PRTPSET; int
        If > 0 print the OUTPUT4 matrix partitioning vector sets
    PRTQSYS; int
        If=1 print matrix QSYS
    PRTRMG; int
        If PRTRMG=1 or 3, print constraint matrix RMG
        If PRTRMG=2 or 3, print partitions RMN and RMM of constraint matrix RMG
    PRTSCP; int;If PRTSCP=1 print data generated in the subcase processor
    PRTSTIFD; int
        Same as PRTSTIFF, except only print diagonal terms
    PRTSTIFF; int; Defaults=0 for PRTSTIFF(I), I=1-4 which go in
        fields 3-6 of the PARAM PRTSTIFF entry that prints sparse
        stiffness matrices for G, N, F, and/or A-sets:
          V1=PRTSTIFF(1)=1 print sparse KGG
          V2=PRTSTIFF(2)=1 or 3 print sparse KNN, 2 or 3 print KNM, KMM
          V3=PRTSTIFF(3)=1 or 3 print sparse KFF, 2 or 3 print KFS, KSS
          V4=PRTSTIFF(4)=1 or 3 print sparse KAA, 2 or 3 print KAO, KOO
          V5=PRTSTIFF(5)=1 or 3 print sparse KLL, 2 or 3 print KLR, KRR
    PRTTSET; int
        If PRTSET=1 print TSET table which gives the character name of
        the displacement sets that each degree of freedom belongs to
        (size is number of grids x 6)
    PRTUO0; int
        If=1 print UO0
    PRTUSET; int
        If > 0 print the user defined set (U1 or U2) definitions
    PRTYS; int
        If=1 print matrix YS
    Q4SURFIT; int; default=6
        Polynomial order for the surface fit of QUAD4 stress/strain when
        stresses are requested for other than corner locations
    QUAD4TYP; char
        'MIN4T' ! Which element to use in MYSTRAN as the QUAD4 element
        'MIN4T (default)': Use Tessler's MIN4T element made up of 4 MIN3 triangles
        'MIN4 ': Use Tessler's MIN4 element
    QUADAXIS; char; default='SIDE12'
        This determines how the quad element local x axis is defined.
        'SIDE12' means that the axis between grids 1 and 2 of the quad define the local x axis.
        'SPLITD' means that the axis is defined as the direction that splits the angle between the quad diagonals

    RCONDK; char
        If RCONDK=Y, then LAPACK calculates the condition number of the A-set stiffness matrix.
        This is required if LAPACK error bounds on the A-set displacement solution are desired.
        This can require significant solution time.
    RELINK3; char
        'Y' or 'N' to specify whether to rerun LINK3 and also LINK5 in a restart
    SETLKTK
        int; Field 3: SETLKTK value. default=0. Method to estimate number of
                      nonzeros in G-set stiffness matrix so array can be allocated.
            (1) If SETLKTK=0, estimate LTERM_KGG based on full element stiffness
                                matrices unconnected (most conservative but not time consuming).
            (2) If SETLKTK=1, estimate LTERM_KGG based on KGG bandwidth.
            (3) If SETLKTK=2, estimate LTERM_KGG based on KGG density of nonzero terms
            (4) If SETLKTK=3, estimate LTERM_KGG based on actual element stiffness matrices unconnected.
            (5) f SETLKTK=4, estimate LTERM_KGG on value input by user in
                 field 5 of the PARAM SETLKT entry (PARAM USR_LTERM_KGG).
        char; field 4: ESP0_PAUSE value (default=N, do not pause after subr
                       ESP0 to let user input LTERM_KGG, or pause if=Y
        int;  field 5: User input value of LTERM_KGG
    SETLKTM
        Same as SETLKTK but for the G-set mass matrix. Only the values
        for SETLKTM=1, 3, 4 are available
    SHRFXFAC; real; default=1x10^6
        Factor used to adjust transverse shear stiffness when user has
        indicated zero shear flexibility for shell elements. The shear
        stiffness will be reset from infinite (zero flexibility) to
        SHRFXFAC times the average of the bending stiffnesses in the 2 planes
    SKIPMKGG; char; default=N
        'Y', 'N' indicator to say whether to skip calculation of MGG,
        KGG in which case MGG, KGG will be read from previously generated,
        and saved, files (LINK1L for KGG, LINK1R for MGG)
    SOLLIB; char; default=IntMKL
        Denotes which library to use for matrix decomposition and equation solution. Options are:
        IntMKL: Intel Math Kernel Library (matrices stored in sparse form)
        LAPACK (matrices stored in band form)
        YaleSMP: (matrices stored in sparse form) – not fully implemented in MYSTRAN
    SORT_MAX; int; default=5
        Max number of times to run algorithm when sorting arrays before fatal message.
    SPARSTOR; char; default=SYM
        If SYM, symmetric matrices are stored with only the terms on and above the diagonal.
        If NONSYM all terms are stored. SYM requires less disk storage but NONSYM can save
        significant time in sparse matrix partitioning and multiply operations.
    STR_CID; int; default=-1.
        Indicator for the coordinate system to use ID for elem stress, strain and emgineering force output:
        -1 is local element coordinate system (default)
        0 is basic coordinate system
        j (any other integer) is a defined coordinate system for output
    SUPINFO; char; default=Y
        Indicator of whether some information messages should be suppressed in the F06 output file.
          N indicates to suppress
          Y indicates to not suppress messages in the file.
    SUPWARN; char; default=Y
       Indicator of whether warning messages should be suppressed in the F06 output file.
        N indicates to suppress, Y indicates to not suppress messages in the file.
    THRESHK; real; default=0.1
        User defined value for the threshold in deciding whether to
        equilibrate the A-set stiffness matrix in LAPACK subroutine
        DLAQSB. Default value 0.1, LAPACK suggests

    TINY; Real
        Do not print matrix values whose absolute value is less than this parameter value
    TRLLOHI; int
        For TRLan eigen extraction (default=-1) - which end of spectrum to compute:
        < 0, the smallest eigenvalues
        = 0, whichever converges first
        > 0, the largest eigenvalues
    TRLMXLAN; int
        For TRLan eigen extraction (default=7) - Max num Lanczos basis vectors to be used
        (If user enters a Bulk Data PARAM TRLMXLAN then internal parameter
        USER_TRLMXLAN is set to 'Y')
    TRLMXMV; int
        For TRLan eigen extraction (default=-2000) - Max number of
        matrix-vector multiplications allowed
    TRLREST; int
       For TRLan eigen extraction (default=1) - Index of restarting schemes
    TRLTOL; real
        For TRLan eigen extraction (default=1.4901D-8) - Eigenpair is
        declared converged if its residual norm is < tol*||OP||
    TRLVERB; int
         For TRLan eigen extraction (default=0) - Level of output data written by TRLan
    TSTM_DEF; real; default=5/6 = 0.833333
        Value for TS/TM on PSHELL Bulk data entry when that field on the
        PSHELL is blank
    USETSTR; char
        Requests output of the internal sequence order for displacement
        sets (e.g. G-set, etc). See section 3.6 for a discussion of
        displacement sets. In addition to the sets in section 3.7, the
        user displacement sets U1 and U2 (see Bulk Data entry USET and
        USET1) can also have the internal sort order output to the F06
        file. As an example, to obtain a row oriented tabular output
        of the internal sort order for the R-set, include the Bulk data entry:
            PARAM, USETSTR, R
    USR_JCT; int
        User supplied value for JCT - used in shell sort subroutines.
        If USR_JCT=0, internal values for JCT will be used in the shell sort.
    WINAMEM; real; default=2.0 GB
        Max memory Windows allows for any array. If it is exceeded, a
        message is printed out and execution is aborted. This is used to
        avoid a failure which aborts MYSTRAN catastrophically (due to a
        system fault).
    WTMASS; real; default=1.0
        Multiplier for mass matrix after the model total mass is output
        in the Grid Point Weight Generator (GPWG). This allows user to
        input mass terms as weight to get model mass properties in weight
        units and then to convert back to mass units after the GPWG has
        run. For example, if the model units are lb-sec2/inch for mass
        and inches for length and the input data file has lb for 'mass'
        (read weight), then 1/386, or 0.002591 would be the value for
        WTMASS needed to convert the 'mass' matrix from weight units to
        mass units.

    Note
    ----
    Default values of parameters are:
        N for Char
        0 for Int
        0.0 for real

    """
    type = 'PARAM'
    _field_map = {1: 'key'}

    def _update_field_helper(self, n, value):
        if n - 2 >= 0:
            try:
                self.values[n - 2] = value
            except IndexError:
                msg = 'Field %r=%r is an invalid %s entry for key=%r.' % (
                    n, value, self.type, self.key.upper())
                raise IndexError(msg)
        else:
            raise KeyError('Field %r=%r is an invalid %s entry.' % (n, value, self.type))

    @classmethod
    def _init_from_empty(cls):
        key = 'POST'
        values = -1
        return PARAM(key, values, comment='')

    def __init__(self, key, values, comment=''):
        """
        Creates a PARAM card

        Parameters
        ----------
        key : str
            the name of the PARAM
        values : int/float/str/list
            varies depending on the type of PARAM
        comment : str; default=''
            a comment for the card

        """
        if comment:
            self.comment = comment
        self.key = key
        if isinstance(values, (int, float, str)):
            values = [values]
        self.values = values

    @classmethod
    def add_card(cls, card, comment=''):
        """
        Adds a PARAM card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card

        """
        key = string(card, 1, 'key')

        if key == 'ARP_TOL':
            values = [
                double_or_blank(card, 2, 'value', 1e-6)
            ]
        elif key == 'AUTOSPC':
            # ['PARAM', 'AUTOSPC', 'Y', '1.0E-9', None, 'Y', 'Y']
            # enabled, tol, nset, info, print_
            enabled = string_or_blank(card, 2, 'ENABLED', 'Y')

            tol = double_or_blank(card, 3, 'TOL', 1e-6)
            nset = integer_or_blank(card, 4, 'NSET', 1)
            info = string_or_blank(card, 5, 'INFO', 'N')
            print_ = string_or_blank(card, 6, 'PRINT', 'N')

            assert enabled in {'Y', 'N'}, enabled
            assert info in {'Y', 'N'}, info
            assert print_ in {'Y', 'N'}, print_
            values = [
                enabled, tol, nset, info, print_
            ]
        elif key == 'ART_MASS':
            # ['PARAM', 'ART_MASS', 'Y', '0.', '1.E-6']
            values = [
                string_or_blank(card, 2, 'is_mass', 'N'),
                double_or_blank(card, 3, 'trans_mass', 1e-6),
                double_or_blank(card, 4, 'rot_mass', 1e-6),
            ]
        elif key == 'EIGNORM2':
            #['PARAM', 'EIGNORM2', 'Y']
            values = [string_or_blank(card, 2, 'enabled', 'N')]
        elif key == 'EQCHECK':
            #['PARAM', 'EQCHECK', '0', '3', '3', None, None, None, '-1.E10']
            values = [
                integer_or_blank(card, 2, 'refgrid', 0),
                integer_or_blank(card, 3, 'g-set', 0),
                integer_or_blank(card, 4, 'n-set', 0),
                integer_or_blank(card, 5, 'f-set', 0),
                integer_or_blank(card, 6, 'a-set', 0),
                integer_or_blank(card, 7, 'l-set', 0),
                double_or_blank(card, 8, 'grid_forces_tol', 1e-5),
                string_or_blank(card, 9, 'norm_forces', 'N'),
            ]
        elif key == 'GRDPNT':
            #['PARAM', 'GRDPNT', '0']
            values = [
                integer_or_blank(card, 2, 'GRDPNT', -1),
            ]
        elif key == 'MEFMLOC':
            # ['PARAM', 'MEFMLOC', 'CG']
            values = [string_or_blank(card, 2, 'enabled', 'GRDPNT')]
            assert values[0] in ['GRDPNT', 'CG', 'GRID'], values
        elif key == 'POST':
            #['PARAM', 'POST', '-1']
            values = [
                integer_or_blank(card, 2, 'POST', -1),
            ]
        elif key == 'RCONDK':
            #['PARAM', 'RCONDK', 'Y']
            values = [
                string_or_blank(card, 2, 'RCONDK', 'N'),
            ]
        elif key == 'SOLLIB':
            # ['PARAM', 'SOLLIB', 'LAPACK']
            values = [
                string_or_blank(card, 2, 'library', 'INTELMKL'),
            ]
            #SOLLIB; char; default=IntMKL
                #Denotes which library to use for matrix decomposition and equation solution. Options are:
                #IntMKL: Intel Math Kernel Library (matrices stored in sparse form)
                #LAPACK (matrices stored in band form)
                #YaleSMP: (matrices stored in sparse form) – not fully implemented in MYSTRAN
            assert values[0] in ['INTMKL', 'LAPACK', 'YALESMP'], values
        elif key == 'SPARSTOR':
            values = [
                string_or_blank(card, 2, 'storage', 'SYM'),
            ]
            assert values[0] in {'SYM', 'NONSYM'}, values
        elif key == 'WTMASS':
            #['PARAM', 'WTMASS', '.002591']
            values = [
                double_or_blank(card, 2, 'WTMASS', 1.0),
            ]

        # basic print flags
        elif key == 'PRTBASIC':
            #['PARAM', 'PRTBASIC', '1']
            values = [
                integer_or_blank(card, 2, 'PRTBASIC', 0),
            ]
            assert values[0] in [0, 1], values
        elif key == 'PRTCGLTM':
            #['PARAM', 'PRTCGLTM', '1']
            values = [
                integer_or_blank(card, 2, 'PRTCGLTM', 0),
            ]
            assert values[0] in [0, 1], values
        elif key == 'PRTDLR':
            #['PARAM', 'PRTDLR', '1']
            values = [
                integer_or_blank(card, 2, 'PRTDLR', 0),
            ]
            assert values[0] in [0, 1], values
        elif key == 'PRTGMN':
            #['PARAM', 'PRTGMN', '1']
            values = [
                integer_or_blank(card, 2, 'PRTGMN', 0),
            ]
            assert values[0] in [0, 1], values
        elif key == 'PRTGOA':
            #['PARAM', 'PRTGOA', '1']
            values = [
                integer_or_blank(card, 2, 'PRTGMN', 0),
            ]
            assert values[0] in [0, 1], values
        elif key == 'PRTSCP':
            #['PARAM', 'PRTSCP', '1']
            values = [
                integer_or_blank(card, 2, 'PRTSCP', 0),
            ]
            assert values[0] in [0, 1], values
        elif key == 'PRTTSET':
            #['PARAM', 'PRTTSET', '1']
            values = [
                integer_or_blank(card, 2, 'PRTTSET', 0),
            ]
            assert values[0] in [0, 1], values

        elif key == 'PRTHMN':
            #['PARAM', 'PRTHMN', '1']
            values = [
                integer_or_blank(card, 2, 'PRTHMN', 0),
            ]
            assert values[0] in [0, 1], values
        elif key == 'PRTIFLTM':
            #['PARAM', 'PRTIFLTM', '1']
            values = [
                integer_or_blank(card, 2, 'PRTIFLTM', 0),
            ]
            assert values[0] in [0, 1], values
        elif key == 'PRTKXX':
            #['PARAM', 'PRTKXX', '1']
            values = [
                integer_or_blank(card, 2, 'PRTKXX', 0),
            ]
            assert values[0] in [0, 1], values
        elif key == 'PRTMXX':
            #['PARAM', 'PRTMXX', '1']
            values = [
                integer_or_blank(card, 2, 'PRTMXX', 0),
            ]
            assert values[0] in [0, 1], values
        elif key == 'PRTPHIXA':
            #['PARAM', 'PRTPHIXA', '1']
            values = [
                integer_or_blank(card, 2, 'PRTPHIXA', 0),
            ]
            assert values[0] in [0, 1], values
        elif key == 'PRTPHIZL':
            #['PARAM', 'PRTPHIZL', '1']
            values = [
                integer_or_blank(card, 2, 'PRTPHIZL', 0),
            ]
            assert values[0] in [0, 1], values
        elif key == 'PRTQSYS':
            #['PARAM', 'PRTQSYS', '1']
            values = [
                integer_or_blank(card, 2, 'PRTQSYS', 0),
            ]
            assert values[0] in [0, 1], values
        elif key == 'PRTUO0':
            #['PARAM', 'PRTUO0', '1']
            values = [
                integer_or_blank(card, 2, 'PRTUO0', 0),
            ]
            assert values[0] in [0, 1], values
        elif key == 'PRTYS':
            #['PARAM', 'PRTYS', '1']
            values = [
                integer_or_blank(card, 2, 'PRTUO0', 0),
            ]
            assert values[0] in [0, 1], values
        # ----------------------------------
        # not basic print flags
        elif key == 'PRTCORD':
            #['PARAM', 'PRTCORD', '2']
            values = [
                integer_or_blank(card, 2, 'PRTCORD', 0),
            ]
            assert values[0] in [0, 1, 2], values # TODO: why is 2 allowed?
        elif key == 'PRTDOF':
            #['PARAM', 'PRTDOF', '1']
            values = [
                integer_or_blank(card, 2, 'PRTDOF', 0),
            ]
            assert values[0] in [0, 1, 2, 3], values
        elif key == 'PRTRMG':
            #['PARAM', 'PRTRMG', '3']
            values = [
                integer_or_blank(card, 2, 'PRTRMG', 0),
            ]
            assert values[0] in [0, 1, 2, 3], values


        elif key == 'PRTSTIFF':
            # ['PARAM', 'PRTSTIFF', '1', '3', '3', '3', '3']
            values = [
                integer(card, 2, 'PRTSTIFF-2'),
                integer(card, 3, 'PRTSTIFF-3'),
                integer(card, 4, 'PRTSTIFF-4'),
                integer(card, 5, 'PRTSTIFF-5'),
                integer(card, 6, 'PRTSTIFF-6'),
            ]
            for i, value in enumerate(values):
                assert value in [1, 3], f'i={i} values={values}'
        elif key == 'PRTSTIFD':
            # ['PARAM', 'PRTSTIFD', '1', '3', '3', '3', '3']
            values = [
                integer(card, 2, 'PRTSTIFD-2'),
                integer(card, 3, 'PRTSTIFD-3'),
                integer(card, 4, 'PRTSTIFD-4'),
                integer(card, 5, 'PRTSTIFD-5'),
                integer(card, 6, 'PRTSTIFD-6'),
            ]
            for i, value in enumerate(values):
                assert value in [1, 3], f'i={i} values={values}'
        elif key == 'PRTFOR':
            # ['PARAM', 'PRTFOR', '1', '3', '3', '3', '3']
            values = [
                integer(card, 2, 'PRTFOR-2'),
                integer(card, 3, 'PRTFOR-3'),
                integer(card, 4, 'PRTFOR-4'),
                integer(card, 5, 'PRTFOR-5'),
                integer(card, 6, 'PRTFOR-6'),
            ]
            for i, value in enumerate(values):
                assert value in [1, 3], f'i={i} values={values}'
        elif key == 'PRTMASS':
            # ['PARAM', 'PRTMASS', '1', '3', '3', '3', '3']
            values = [
                integer(card, 2, 'PRTMASS-2'),
                integer(card, 3, 'PRTMASS-3'),
                integer(card, 4, 'PRTMASS-4'),
                integer(card, 5, 'PRTMASS-5'),
                integer(card, 6, 'PRTMASS-6'),
            ]
            for i, value in enumerate(values):
                assert value in [1, 3], f'i={i} values={values}'
        else:
            raise NotImplementedError(card)

        return PARAM_MYSTRAN(key, values, comment=comment)

    def raw_fields(self):
        list_fields = ['PARAM', self.key] + self.values
        return list_fields

    def repr_fields(self):
        return self.raw_fields()

    def write_card(self, size: int=8, is_double: bool=False) -> str:
        card = self.raw_fields()
        if size == 8:
            return self.comment + print_card_8(card)
        return self.comment + print_card_16(card)


MDLPRM_INT_KEYS_1 = {
    'BRTOBM', 'BUSHRT', 'DCFLTEXP', 'GEV1417', 'GNLSTN',
    'HDF5', 'H5GM34', 'H5INFO', 'H5MDL', 'H5MTX', 'H5NORDOF', 'H5XHH',
    'IGNSHBDN', 'INTOUT', 'LMT2MPC', 'MLTSPLIN', 'MPCF129', 'NLDIFF',
    'NONUPIV', 'NSGRDS4', 'PRDIDPVT', 'PRDITRFN', # 'PIVTHRSH',
    'PRDMTYPE', 'PRDOOC', 'PRDWMTCH', 'QR6ROT', 'QRSHEAR',
    'RDBOTH', 'RELAXF', 'REUPSE', 'RMRBE3RT', 'RSTIGNDP', 'SHRTOQ4',
    'STREQCNT', 'TWBRBML',

    # undefined in MSC
    'RBEDOF', 'NLDIAG', 'ITRFMT', 'NLSPCD', 'MRCONV', 'LA3FLG', 'TIMADJ',
}
MDLPRM_STR_KEYS_1 = {'COMPN1', 'SHEARP', 'OFFDEF',
                     'PRTELAS', 'PRTFAST', 'PRTMASS', 'PRTSEAM', 'PRTWELD'}
MDLPRM_FLOAT_KEYS_1 = {
    'DBCTOLE', 'DELELAS', 'DELFAST', 'DELMASS', 'DELSEAM', 'DELWELD',
    'PEXTS4', 'PIVTHRSH', 'SPBLNDX', 'PEXTS4', }
MDLPRM_KEYS = MDLPRM_INT_KEYS_1 | MDLPRM_STR_KEYS_1

class MDLPRM(BaseCard):
    """MSC Nastran card"""
    type = 'MDLPRM'
    _field_map = {1: 'key'}

    @classmethod
    def _init_from_empty(cls):
        key = 'POST'
        values = -1
        return MDLPRM(key, values, comment='')

    def __init__(self, mdlprm_dict: dict[str, Union[int, float]], comment=''):
        """
        Creates a MDLPRM card

        Parameters
        ----------
        mdlprm_dict : dict[name, value]
            name : str
                the name of the MDLPRM
            value: int/float
                varies depending on the type of MDLPRM
        comment : str; default=''
            a comment for the card

        """
        if comment:
            self.comment = comment
        self.mdlprm_dict = mdlprm_dict

        assert isinstance(mdlprm_dict, dict), mdlprm_dict
        for key, value in sorted(mdlprm_dict.items()):
            if key in MDLPRM_INT_KEYS_1:
                assert isinstance(value, integer_types), f'MDLPRM key={key!r} value={value!r} must be an integer'
            elif key in MDLPRM_STR_KEYS_1:
                assert isinstance(value, str), f'MDLPRM key={key!r} value={value!r} must be an integer'
            elif key in MDLPRM_FLOAT_KEYS_1:
                assert isinstance(value, float_types), f'MDLPRM key={key!r} value={value!r} must be an float'
            else:
                raise RuntimeError(f'MDLPRM key={key!r} value={value!r} is not supported')
        assert len(mdlprm_dict) > 0, mdlprm_dict

    @classmethod
    def add_card(cls, card: BDFCard, comment=''):
        """
        Adds a MDLPRM card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card

        """
        mdlprm_dict = {}
        ifield = 1
        while ifield < len(card):
            key = string_or_blank(card, ifield, 'key')
            if key is None:
                value = blank(card, ifield+1, 'blank_value')
                ifield += 2
                continue
            elif key in MDLPRM_INT_KEYS_1:
                value = integer_or_blank(card, ifield+1, 'value')
            elif key in MDLPRM_STR_KEYS_1:
                value = string(card, ifield+1, key)
            elif key in MDLPRM_FLOAT_KEYS_1:
                value = double(card, ifield+1, key)
            else:
                raise RuntimeError(f'MDLPRM key={key!r} is not supported; value={card.field(ifield+1)}')
            mdlprm_dict[key] = value
            ifield += 2
        obj = MDLPRM(mdlprm_dict, comment=comment)
        return obj

    #@classmethod
    def export_to_hdf5(self, h5_file, model, desvar_ids):
        """exports the elements in a vectorized way"""
        #comments = []
        encoding = model.get_encoding()
        keys = []
        flags = []
        value_ints = []
        value_floats = []
        value_strs = []
        for key, value in sorted(self.mdlprm_dict.items()):
            value_int = 0
            value_float = np.nan
            value_str = b''
            if isinstance(value, integer_types):
                value_int = value
                flag = 1
            elif isinstance(value, float_types):
                value_float = value
                flag = 2
            elif isinstance(value, str):
                value_str = value.encode(encoding)
                flag = 3
            else:
                raise NotImplementedError((key, value))
            keys.append(key.encode(encoding))
            flags.append(flag)
            value_ints.append(value_int)
            value_floats.append(value_float)
            value_strs.append(value_str)
        h5_file.create_dataset('keys', data=keys)
        h5_file.create_dataset('ifs_flag', data=flags)
        h5_file.create_dataset('value_ints', data=value_ints)
        h5_file.create_dataset('value_floats', data=value_floats)
        h5_file.create_dataset('value_strs', data=value_strs)

    def load_hdf5_file(self, hdf5_file, encoding) -> None:
        from pyNastran.utils.dict_to_h5py import _cast

        keys = list(hdf5_file.keys())
        assert len(keys) == 5, keys
        ifs_flag = _cast(hdf5_file['ifs_flag'])
        keys = _cast(hdf5_file['keys'])
        value_ints = _cast(hdf5_file['value_ints'])
        value_floats = _cast(hdf5_file['value_floats'])
        value_strs = _cast(hdf5_file['value_strs'])
        for key, flag, xint, xfloat, xstr in zip(keys, ifs_flag, value_ints, value_floats, value_strs):
            key_str = key.decode(encoding)
            if flag == 1:
                value = xint
            elif flag == 2:
                value = xfloat
            elif flag == 3:
                value = xstr.decode(encoding)
                assert isinstance(value, str), value
            else:  # pragma: no cover
                raise RuntimeError((key_str, flag, xint, xfloat, xstr))
            #print(key_str, value)
            self.mdlprm_dict[key_str] = value
        #print(str(self))
        return

    def raw_fields(self):
        list_fields = ['MDLPRM']
        for key, value in self.mdlprm_dict.items():
            list_fields += [key, value]
        assert len(list_fields) > 1, self.mdlprm_dict
        return list_fields

    def repr_fields(self):
        return self.raw_fields()

    def write_card(self, size: int=8, is_double: bool=False) -> str:
        card = self.raw_fields()
        if size == 8:
            return self.comment + print_card_8(card)
        return self.comment + print_card_16(card)

