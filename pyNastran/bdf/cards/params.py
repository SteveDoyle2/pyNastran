"""
defines the following card:
 - PARAM
"""
# pylint: disable=C0103,R0902,R0904,R0914
from pyNastran.bdf.cards.base_card import BaseCard
from pyNastran.bdf.bdf_interface.bdf_card import BDFCard
from pyNastran.bdf.bdf_interface.assign_type import (
    integer, double, integer_or_blank, double_or_blank, string, string_or_blank,
    integer_double_string_or_blank)
from pyNastran.bdf.field_writer_8 import print_card_8
from pyNastran.bdf.field_writer_16 import print_card_16

#float_words_1 = [
    #b'K6ROT', b'WTMASS', b'SNORM', b'PATVER', b'MAXRATIO', b'EPSHT',
    #b'SIGMA', b'TABS']

SMALL_FIELD_PARAMS = [
    'ACOUT', 'ACOWEAK', 'ACSYM', 'ADJMETH', 'AESMAXIT', 'AESMETH', 'ADSTAT',
    'MAXLINES'] #+ INT_WORDS_1 + STR_WORDS_1


# per NX 11 QRG
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
    'ALTRED', 'OUGCORD', 'RSCON', 'RESVINER', 'MESH', 'SKINOUT', 'VUPENTA',
    'DYNSEN', 'PRTGPL', 'PRTEQXIN', 'PRTGPDT', 'PRTCSTM', 'PRTBGPDT', 'PRTGPTT',
    'PRTMGG', 'PRTPG', 'FOLLOWK', 'DBCCONV', 'COMPMATT', 'ADB', 'HEATSTAT',
    'AUTOSPCR', 'CHECKOUT', 'EPSILONT', 'RMS', 'RESPATH', 'AUTOMSET', 'COMPMATT',
    'F56', 'METHFL', 'DMIGNRG', 'METHT', 'AEDB', 'SAVEOFP', 'DBDN', 'DBCONV',
    'ASCOUP', 'DBUP', 'VUHEXA', 'FLEXINCR', 'SRCOMPS', 'VUTETRA', 'DYNSPCF',
    'SRCOMPS', 'SERST', 'PGRPST', 'CDPRT', 'CFDIAGP', 'CWDIAGP', 'ENFMOTN', 'RMSINT',
    'DEBUG', 'CHKOUT', 'DBALL', 'OIBULK', 'RMXTRAN', 'CDPCH', 'AUTOADJ', 'SEMAP',
    'SOFTEXIT', 'SM', 'RESFLEX', 'SENSUOO', 'OGEM', 'XBYMODE', 'STRESS', 'OPTION',
    'AOTOSPC', 'ARBMAS', 'ARBMSS', 'OMACHPR', 'AUTOMPC', 'BSHDAMP', 'MDOF', 'ASCII4',
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
        values : int/float/str/List
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

        n = 1
        value = None
        if key == 'ACOUT':
            value = string_or_blank(card, 2, 'value', 'PEAK')
        elif key == 'ACOWEAK':
            value = string_or_blank(card, 2, 'value', 'NO')
        elif key == 'ADJMETH':
            value = integer_or_blank(card, 2, 'value', 0)
        elif key == 'ADPCON':
            value = double_or_blank(card, 2, 'value', 1.0)
        #elif key == 'ADMPOST':
            #value = string_or_blank(card, 2, 'value', 0) ## TODO: 0 is not a string
        elif key == 'ADSDISC':
            value = double_or_blank(card, 2, 'value', 1e-8)
        elif key == 'AESMAXIT':
            value = integer_or_blank(card, 2, 'value', 15)
        elif key == 'AESMETH':
            value = string_or_blank(card, 2, 'value', 'SELECT')
            assert value in ['SELECT', 'AUTO', 'DIRECT', 'RITZ', 'ITER'], 'value=%s' % value
        elif key == 'AESTOL':
            value = double_or_blank(card, 2, 'value', 1e-10)
        elif key in ['ALPHA1FL', 'ALPHA2FL']:  # check alpha1/alpha1FL
            value1 = double_or_blank(card, 2, 'value1', 0.0)
            value2 = double_or_blank(card, 3, 'value2', 0.0)
            n = 2
        elif key == 'COMPMATT':
            #('COMPMATT', 'NO', ['YES', 'NO', 'NONSMEAR']), # MSC only: 'NONSMEAR'
            value = string_or_blank(card, 2, 'value1', 'NO')
            if value == 'NONS':  # assume
                value = 'NONSMEAR'
            if value == 'SMEAR':  # assume
                value = 'YES'
            assert value in ['YES', 'NO', 'NONSMEAR'], 'value=%r' % value

        elif key == 'POST':
            value = integer_or_blank(card, 2, 'value', 1)
        elif key == 'UNITSYS':
            value = string(card, 2, 'value')

        #-------------------------------------------------------------
        # strings; has defaults
        elif key in string_params:
            default, allowed_values = string_params[key]
            value = string_or_blank(card, 2, 'value', default=default)
            assert value in allowed_values, 'value=%s allowed=%s' % (value, allowed_values)

        # ints; has defaults
        elif key in int_params:
            default = int_params[key]
            value = integer_or_blank(card, 2, 'value', default=default)
        elif key in int_params_allowed:
            default, allowed_values = int_params_allowed[key]
            value = integer_or_blank(card, 2, 'value', default=default)
            assert value in allowed_values, 'value=%s allowed=%s' % (value, allowed_values)

        # floats; has defaults
        elif key in float_params:
            default = float_params[key]
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
                assert ' ' not in value1, 'PARAM value1=%r' % value1
            if isinstance(value2, str):
                assert ' ' not in value2, 'PARAM value2=%r' % value2
            values = [value1, value2]
        else:
            # n=1
            if isinstance(value, str):
                assert ' ' not in value, 'PARAM value=%r' % value
            values = [value]

        if n == 1:
            assert len(card) <= 3, 'len(PARAM card)=%i card=%r' % (len(card), card)
        else:
            assert len(card) <= 4, 'len(PARAM card)=%i card=%r' % (len(card), card)
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
                msg = 'key=%s value1=%r must be an string.' % (self.key, value1)
                raise TypeError(msg)

        elif self.key == 'ACOWEAK':
            if value1 is None:
                value1 = 'NO'
            if not isinstance(value1, str):
                msg = 'key=%s value1=%r must be an string.' % (self.key, value1)
                raise TypeError(msg)

        elif self.key == 'ACSYM':
            if value1 is None:
                value1 = 'YES'
            if not isinstance(value1, str):
                msg = 'key=%s value1=%r must be an string.' % (self.key, value1)
                raise TypeError(msg)

        elif self.key == 'ADJMETH':
            if value1 is None:
                value1 = 0
            if not isinstance(value1, int):
                msg = 'key=%s value1=%r must be an integer.' % (self.key, value1)
                raise TypeError(msg)

        #elif self.key == 'ADMPOST': ## TODO: 0 is not a string
            #value = string_or_blank(card, 2, 'value', 0)

        elif self.key == 'ADSTAT':
            if value1 is None:
                value1 = 'YES'
            if not isinstance(value1, str):
                msg = 'key=%s value1=%r must be an string.' % (self.key, value1)
                raise TypeError(msg)

        elif self.key in ['ALPHA1', 'ALPHA2', 'ALPHA1FL', 'ALPHA2FL']:
            if value1 is None:
                value1 = 0.0
            if value2 is None:
                value2 = 0.0
            if not isinstance(value1, float):
                msg = 'key=%s value1=%r must be an float.' % (self.key, value1)
                raise TypeError(msg)
            if isinstance(value2, float):
                msg = 'key=%s value2=%r must be an float.' % (self.key, value2)
                raise TypeError(msg)

        elif self.key in ['CB1', 'CB2', 'CK1', 'CK2', 'CK3', 'CM1', 'CM2', 'CP1', 'CP2']:
            if value1 is None:
                value1 = 1.0
            if value2 is None:
                value2 = 0.0
            if not isinstance(value1, float):
                msg = 'key=%s value1=%r must be an float.' % (self.key, value1)
                raise TypeError(msg)
            if isinstance(value2, float):
                msg = 'key=%s value2=%r must be an float.' % (self.key, value2)
                raise TypeError(msg)

        else:
            if not isinstance(value1, (int, float, str)):
                msg = 'key=%s value1=%r must be an integer, float, or string.' % (self.key, value1)
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
        if self.key in INT_STR_WORDS_1:
            return '%sPARAM   %8s%8s\n' % (
                self.comment, self.key, self.values[0])

        if size == 8 or self.key in SMALL_FIELD_PARAMS:
            return self.comment + print_card_8(card)
        return self.comment + print_card_16(card)
