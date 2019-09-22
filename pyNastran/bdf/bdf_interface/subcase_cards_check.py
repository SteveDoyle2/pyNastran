from typing import List, Dict, Tuple, Union, Set, Any
from pyNastran.bdf.bdf_interface.subcase_cards import CheckCard


class DISPLACEMENT(CheckCard):
    """
    DISPLACEMENT=5
    DISPLACEMENT(REAL)=ALL
    DISPLACEMENT(SORT2, PUNCH, REAL)=ALL

    DISPLACEMENT=5
    DISPLACEMENTS(REAL)=ALL
    DISPLACEMENT(SORT2, PUNCH, REAL)=ALL
    DISPLACEMENT(SORT2, PRINT, PSDF, CRMS, RPUNCH)=20
    DISPLACEMENT(PRINT, RALL, NORPRINT)=ALL
    DISP(T1=1.-3, T3=1.-2) = ALL
    DISP(TM=1.-3, PRINT,PLOT) = ALL
    DISP(TM=1.-3,PRINT,PLOT,SORT2) = 20

    """
    type = 'DISPLACEMENT'
    alternate_names = {'DISP',
                       'PRESS', 'PRESSURE'}
    allowed_keys = {
        'SORT1', 'SORT2',
        'PRINT', 'PUNCH', 'PLOT',
        'REAL', 'IMAG', 'PHASE',
        'ABS', 'REL', 'PSDF', 'ATOC', 'CRMS', 'RMS', 'RALL', 'RPRINT',
        'NOPRINT', 'RPUNCH',
        'TM', 'T1', 'T2', 'T3',
        'RM', 'R1', 'R2', 'R3',
    }
    allowed_strings = {'ALL', 'NONE'}
    duplicate_names = {
        'PHAS' : 'PHASE',
    }
    allowed_values = {}  # type: Dict[str, Union[str, int]]
    allow_ints = True

    def __init__(self, key, value, options):
        CheckCard.__init__(self, key, value, options)

class VELOCITY(CheckCard):
    """
    VELOCITY=5
    VELOCITY(REAL)=ALL
    VELOCITY(SORT2, PUNCH, REAL)=ALL

    """
    type = 'VELOCITY'
    alternate_names = {'VELO'}

    allowed_keys = {
        'SORT1', 'SORT2', 'PRINT', 'PUNCH', 'PLOT', 'REAL', 'IMAG', 'PHASE',
        'ABS', 'REL', 'PSDF', 'ATOC', 'CRMS', 'RMS', 'RALL', 'RPRINT',
        'NOPRINT', 'RPUNCH',
    }
    allowed_strings = {'ALL', 'NONE'}
    duplicate_names = {
        'PHAS' : 'PHASE',
    }
    allowed_values = {}  # type: Dict[str, Union[str, int]]
    allow_ints = True

    def __init__(self, key, value, options):
        CheckCard.__init__(self, key, value, options)

class ACCELERATION(CheckCard):
    """
    ACCELERATION=5
    ACCELERATION(REAL)=ALL
    ACCELERATION(SORT2, PUNCH, REAL)=ALL

    """
    type = 'ACCELERATION'
    alternate_names = {'ACCE'}
    allowed_keys = {
        'SORT1', 'SORT2', 'PRINT', 'PUNCH', 'PLOT', 'REAL', 'IMAG', 'PHASE',
        'ABS', 'REL', 'PSDF', 'ATOC', 'CRMS', 'RMS', 'RALL', 'RPRINT',
        'NOPRINT', 'RPUNCH',
    }
    allowed_strings = {'ALL', 'NONE'}
    allowed_values = {}  # type: Dict[str, Union[str, int]]
    allow_ints = True

    def __init__(self, key, value, options):
        CheckCard.__init__(self, key, value, options)

class VECTOR(CheckCard):
    """
    VECTOR=5
    VECTOR(REAL)=ALL
    VECTOR(SORT2, PUNCH, REAL)=ALL

    """
    type = 'VECTOR'
    alternate_names = {'VECT'}
    allowed_keys = {
        'SORT1', 'SORT2', 'PRINT', 'PUNCH', 'PLOT', 'REAL', 'IMAG', 'PHASE',
        #'ABS', 'REL', 'PSDF', 'ATOC', 'CRMS', 'RMS', 'RALL', 'RPRINT',
        #'NOPRINT', 'RPUNCH',
    }
    allowed_strings = {'ALL', 'NONE'}
    duplicate_names = {
        'PHAS' : 'PHASE',
    }
    allowed_values = {}  # type: Dict[str, Union[str, int]]
    allow_ints = True

    def __init__(self, key, value, options):
        CheckCard.__init__(self, key, value, options)

class SDISPLACEMENT(CheckCard):
    """
    SDISPLACEMENT=5
    SDISPLACEMENT(REAL)=ALL
    SDISPLACEMENT(SORT2, PUNCH, REAL)=ALL

    """
    type = 'SDISPLACEMENT'
    alternate_names = {'SDISP'}
    allowed_keys = {
        'SORT1', 'SORT2', 'PRINT', 'PUNCH', 'PLOT', 'REAL', 'IMAG', 'PHASE',
        #'ABS', 'REL', 'PSDF', 'ATOC', 'CRMS', 'RMS', 'RALL', 'RPRINT',
        #'NOPRINT', 'RPUNCH',
    }
    allowed_strings = {'ALL', 'NONE'}
    allowed_values = {}  # type: Dict[str, Union[str, int]]
    allow_ints = True

    def __init__(self, key, value, options):
        CheckCard.__init__(self, key, value, options)


class SVELOCITY(CheckCard):
    """
    SVELOCITY=5
    SVELOCITY(REAL)=ALL
    SVELOCITY(SORT2, PUNCH, REAL)=ALL

    """
    type = 'SVELOCITY'
    alternate_names = {'SVELO'}
    allowed_keys = {
        'SORT1', 'SORT2', 'PRINT', 'PUNCH', 'PLOT', 'REAL', 'IMAG', 'PHASE',
        #'ABS', 'REL', 'PSDF', 'ATOC', 'CRMS', 'RMS', 'RALL', 'RPRINT',
        #'NOPRINT', 'RPUNCH',
    }
    allowed_strings = {'ALL', 'NONE'}
    allowed_values = {}  # type: Dict[str, Union[str, int]]
    allow_ints = True

    def __init__(self, key, value, options):
        CheckCard.__init__(self, key, value, options)


class SACCELERATION(CheckCard):
    """
    SACCELERATION=5
    SACCELERATION(REAL)=ALL
    SACCELERATION(SORT2, PUNCH, REAL)=ALL

    """
    type = 'SACCELERATION'
    alternate_names = {'SACCE', 'SACCEL'}
    allowed_keys = {
        'SORT1', 'SORT2', 'PRINT', 'PUNCH', 'PLOT', 'REAL', 'IMAG', 'PHASE',
        #'ABS', 'REL', 'PSDF', 'ATOC', 'CRMS', 'RMS', 'RALL', 'RPRINT',
        #'NOPRINT', 'RPUNCH',
    }
    allowed_strings = {'ALL', 'NONE'}
    allowed_values = {}  # type: Dict[str, Union[str, int]]
    allow_ints = True

    def __init__(self, key, value, options):
        CheckCard.__init__(self, key, value, options)

class SVECTOR(CheckCard):
    """
    SVECTOR=5
    SVECTOR(REAL)=ALL
    SVECTOR(SORT2, PUNCH, REAL)=ALL

    """
    type = 'SVECTOR'
    alternate_names = {'SVEC'}
    allowed_keys = {
        'SORT1', 'SORT2', 'PRINT', 'PUNCH', 'PLOT', 'REAL', 'IMAG', 'PHASE',
        #'ABS', 'REL', 'PSDF', 'ATOC', 'CRMS', 'RMS', 'RALL', 'RPRINT',
        #'NOPRINT', 'RPUNCH',
    }
    allowed_strings = {'ALL', 'NONE'}
    allowed_values = {}  # type: Dict[str, Union[str, int]]
    allow_ints = True

    def __init__(self, key, value, options):
        CheckCard.__init__(self, key, value, options)

class SPCFORCES(CheckCard):
    """
    SPCFORCES=5
    SPCFORCES(REAL)=ALL
    SPCFORCES(SORT2, PUNCH, REAL)=ALL

    """
    type = 'SPCFORCES'
    alternate_names = {'SPCF'}
    allowed_keys = {
        'SORT1', 'SORT2', 'PRINT', 'PUNCH', 'PLOT', 'REAL', 'IMAG', 'PHASE',
        'ABS', 'REL', 'PSDF', 'ATOC', 'CRMS', 'RMS', 'RALL', 'RPRINT',
        'NOPRINT', 'RPUNCH',
    }
    allowed_strings = {'ALL', 'NONE'}
    allowed_values = {}  # type: Dict[str, Union[str, int]]
    allow_ints = True

    def __init__(self, key, value, options):
        CheckCard.__init__(self, key, value, options)

class MPCFORCES(CheckCard):
    """
    MPCFORCES=5
    MPCFORCES(REAL)=ALL
    MPCFORCES(SORT2, PUNCH, REAL)=ALL

    """
    type = 'MPCFORCES'
    alternate_names = {'MPCF'}
    allowed_keys = {
        'SORT1', 'SORT2', 'PRINT', 'PUNCH', 'PLOT', 'REAL', 'IMAG', 'PHASE',
        'ABS', 'REL', 'PSDF', 'ATOC', 'CRMS', 'RMS', 'RALL', 'RPRINT',
        'NOPRINT', 'RPUNCH',
    }
    allowed_strings = {'ALL', 'NONE'}
    allowed_values = {}  # type: Dict[str, Union[str, int]]
    allow_ints = True

    def __init__(self, key, value, options):
        CheckCard.__init__(self, key, value, options)

class NLLOAD(CheckCard):
    """
    NLLOAD(PRINT, PUNCH)=ALL
    NLLOAD(PRINT, PUNCH)=N
    NLLOAD(PRINT, PUNCH)=NONE

    """
    type = 'NLLOAD'
    #short_name = type
    allowed_keys = {'PRINT', 'PUNCH', 'PLOT'}
    allowed_strings = {'ALL', 'NONE'}
    allowed_values = {}  # type: Dict[str, Union[str, int]]
    allow_ints = True

    def __init__(self, key, value, options):
        CheckCard.__init__(self, key, value, options)

class NLSTRESS(CheckCard):
    """
    NLSTRESS=5
    NLSTRESS (SORT1,PRINT,PUNCH,PHASE)=15
    NLSTRESS(PLOT)=ALL

    """
    type = 'NLSTRESS'
    #short_name = type
    allowed_keys = {'SORT1', 'SORT2', 'PRINT', 'PUNCH', 'PLOT'}
    allowed_strings = {'ALL', 'NONE'}
    allowed_values = {}  # type: Dict[str, Union[str, int]]
    allow_ints = True

    def __init__(self, key, value, options):
        CheckCard.__init__(self, key, value, options)


class NOUTPUT(CheckCard):
    """
    NOUTPUT (R)=ALL
    NOUTPUT (2)=5
    NOUTPUT (4,L)=10
    """
    # doesn't work right...
    type = 'NOUTPUT'
    allowed_keys = {'L', 'R'}
    allowed_strings = {'ALL'}
    #allowed_values = {}
    allow_ints = True

    #def __init__(self, key, value, options):
        #super(NLLOAD, self).__init__(key, value, options)

class OLOAD(CheckCard):
    """
    OLOAD=ALL
    OLOAD(SORT1,PHASE)=5

    """
    type = 'OLOAD'
    alternate_names = {'OLOA'}
    allowed_keys = {'SORT1', 'SORT2', 'PRINT', 'PUNCH', 'PLOT',
                    'REAL', 'IMAG', 'PHASE', 'PSDF', 'ATOC', 'CRMS',
                    'RMS', 'RALL', 'RPRINT', 'NORPRINT', 'RPUNCH'}
    allowed_strings = {'ALL', 'NONE'}
    allowed_values = {}  # type: Dict[str, Union[str, int]]
    allow_ints = True

    def __init__(self, key, value, options):
        CheckCard.__init__(self, key, value, options)

class PRESSURE(CheckCard):
    """
    PRESSURE=5
    PRESSURE(REAL)=ALL
    PRESSURE(SORT2,PUNCH,REAL)=ALL

    """
    type = 'PRESSURE'
    allowed_keys = {'PRINT', 'PUNCH', 'PLOT'
                    'SORT1', 'SORT2',
                    'REAL'}
    allowed_strings = {'ALL'} # , 'NONE'}
    #allowed_values = {}  # type: Dict[str, Union[str, int]]
    alternate_names = {'PRES'}
    allow_ints = True

    def __init__(self, key, value, options):
        CheckCard.__init__(self, key, value, options)

class OPRESS(CheckCard):
    """
    OPRESS=ALL
    OPRESS(PRINT,PUNCH)=17

    """
    type = 'OPRESS'
    allowed_keys = {'PRINT', 'PUNCH'}
    allowed_strings = {'ALL', 'NONE'}
    allowed_values = {}  # type: Dict[str, Union[str, int]]
    allow_ints = True

    def __init__(self, key, value, options):
        CheckCard.__init__(self, key, value, options)

class OTEMP(CheckCard):
    """
    OTEMP=ALL
    OTEMP(PRINT,PUNCH)=17

    """
    type = 'OTEMP'
    allowed_keys = {'PRINT', 'PUNCH'}
    allowed_strings = {'ALL', 'NONE'}
    allowed_values = {}  # type: Dict[str, Union[str, int]]
    allow_ints = True

    def __init__(self, key, value, options):
        CheckCard.__init__(self, key, value, options)


#----------------------------------
# didn't check...
class STRESS(CheckCard):
    """
    STRESS=5
    STRESS(REAL)=ALL
    STRESS(SORT2, PUNCH, REAL)=ALL

    """
    type = 'STRESS'
    alternate_names = {'STRE', 'STRES',
                       'ELSTR', 'ELSTRESS'}
    allowed_keys = {
        'SORT1', 'SORT2',
        'PRINT', 'PUNCH', 'PLOT',
        'REAL', 'IMAG', 'PHASE',
        'CORNER', 'BILIN', 'CENTER', 'CUBIC', 'SGAGE', 'CORN',
        'FIBER', 'STRCUR',
        'VONM', 'VONMISES',
        'MAXSHEAR', 'SHEAR', 'MAXS',
        'PSDF', 'RMS', 'RALL',
        'CPLYBMT', 'CPLYBT', 'CPLYMID',
        #'ABS', 'REL', 'PSDF', 'ATOC', 'CRMS', 'RMS', 'RALL', 'RPRINT',
        #'NOPRINT', 'RPUNCH',
    }
    allowed_strings = {'ALL', 'NONE'}
    allowed_values = {}  # type: Dict[str, Union[str, int]]
    allow_ints = True

    def __init__(self, key, value, options):
        CheckCard.__init__(self, key, value, options)

class STRAIN(CheckCard):
    """
    STRAIN=5
    STRAIN(REAL)=ALL
    STRAIN(SORT2, PUNCH, REAL)=ALL

    """
    type = 'STRAIN'
    alternate_names = {'STRA', 'STRAI'}
    allowed_keys = {
        'SORT1', 'SORT2',
        'PRINT', 'PUNCH', 'PLOT',
        'REAL', 'IMAG', 'PHASE',
        'CORNER', 'BILIN', 'CENTER', 'CUBIC', 'CORN',
        'FIBER', 'STRCUR',
        'VONM', 'VONMISES',
        'MAXSHEAR', 'SHEAR', 'MAXS',
        #'ABS', 'REL', 'PSDF', 'ATOC', 'CRMS', 'RMS', 'RALL', 'RPRINT',
        #'NOPRINT', 'RPUNCH',
    }
    allowed_strings = {'ALL', 'NONE'}
    allowed_values = {}  # type: Dict[str, Union[str, int]]
    allow_ints = True

    def __init__(self, key, value, options):
        CheckCard.__init__(self, key, value, options)

class FORCE(CheckCard):
    """
    FORCE=5
    FORCE(REAL)=ALL
    FORCE(SORT2, PUNCH, REAL)=ALL

    """
    type = 'FORCE'
    alternate_names = {'FORC', 'ELFO', 'ELFORCE'}
    allowed_keys = {
        'SORT1', 'SORT2',
        'PRINT', 'PUNCH', 'PLOT',
        'REAL', 'IMAG', 'PHASE',
        'CORNER', 'BILIN', 'CUBIC', 'CENTER',
        'SGAGE',
        'RALL', 'RPRINT',
        #'ABS', 'REL', 'PSDF', 'ATOC', 'CRMS', 'RMS', 'RALL', 'RPRINT',
        #'NOPRINT', 'RPUNCH',
    }
    allowed_strings = {'ALL', 'NONE'}
    allowed_values = {}  # type: Dict[str, Union[str, int]]
    allow_ints = True

    def __init__(self, key, value, options):
        CheckCard.__init__(self, key, value, options)

class GPFORCE(CheckCard):
    """
    GPFORCE=5
    GPFORCE(REAL)=ALL
    GPFORCE(SORT2, PUNCH, REAL)=ALL

    """
    type = 'GPFORCE'
    alternate_names = {'GPFO'}
    allowed_keys = {
        'SORT1', 'SORT2', 'PRINT', 'PUNCH', 'PLOT', 'REAL', 'IMAG', 'PHASE',
        #'ABS', 'REL', 'PSDF', 'ATOC', 'CRMS', 'RMS', 'RALL', 'RPRINT',
        #'NOPRINT', 'RPUNCH',
    }
    allowed_strings = {'ALL', 'NONE'}
    allowed_values = {}  # type: Dict[str, Union[str, int]]
    allow_ints = True

    def __init__(self, key, value, options):
        CheckCard.__init__(self, key, value, options)

class EDE(CheckCard):
    """
    EDE=5
    EDE(AMPL) = ALL

    """
    type = 'EDE'
    allowed_keys = {
        'IMAG', 'PHASE',
        'PEAK',
        'AMPL', 'AVGAMP',
        #'SORT1', 'SORT2', 'PRINT', 'PUNCH', 'PLOT', 'REAL', 'IMAG', 'PHASE',
        #'PEAK',
        #'NOPRINT',
        #'ABS', 'REL', 'PSDF', 'ATOC', 'CRMS', 'RMS', 'RALL', 'RPRINT',
        #'NOPRINT', 'RPUNCH',
    }
    allowed_strings = {'ALL', } # 'NONE'
    allowed_values = {}  # type: Dict[str, Union[str, int]]
    allow_ints = True

    def __init__(self, key, value, options):
        CheckCard.__init__(self, key, value, options)

class ESE(CheckCard):
    """
    ESE=5
    ESE(REAL)=ALL
    ESE(SORT2, PUNCH, REAL)=ALL

    """
    type = 'ESE'
    allowed_keys = {
        'SORT1', 'SORT2', 'PRINT', 'PUNCH', 'PLOT', 'REAL', 'IMAG', 'PHASE',
        'PEAK',
        'NOPRINT',
        #'ABS', 'REL', 'PSDF', 'ATOC', 'CRMS', 'RMS', 'RALL', 'RPRINT',
        #'NOPRINT', 'RPUNCH',
    }
    allowed_strings = {'ALL', 'NONE'}
    allowed_values = {}  # type: Dict[str, Union[str, int]]
    allow_ints = True

    def __init__(self, key, value, options):
        CheckCard.__init__(self, key, value, options)

class GPSTRESS(CheckCard):
    """
    GPSTRESS=5
    GPSTRESS(REAL)=ALL
    GPSTRESS(SORT2, PUNCH, REAL)=ALL

    """
    type = 'GPSTRESS'
    alternate_names = {'GPST'}
    allowed_keys = {
        'SORT1', 'SORT2', 'PRINT', 'PUNCH', 'PLOT', 'REAL', 'IMAG', 'PHASE',
        #'ABS', 'REL', 'PSDF', 'ATOC', 'CRMS', 'RMS', 'RALL', 'RPRINT',
        #'NOPRINT', 'RPUNCH',
    }
    allowed_strings = {'ALL', 'NONE'}
    allowed_values = {}  # type: Dict[str, Union[str, int]]
    allow_ints = True

    def __init__(self, key, value, options):
        CheckCard.__init__(self, key, value, options)


class GPSTRAIN(CheckCard):
    """
    GPSTRAIN=5
    GPSTRAIN(REAL)=ALL
    GPSTRAIN(SORT2, PUNCH, REAL)=ALL

    """
    type = 'GPSTRAIN'
    alternate_names = {'GPST'}
    allowed_keys = {
        'PLOT', 'PRINT', 'PUNCH',
        'FIBER',
        #'SORT1', 'SORT2', 'PRINT', 'PUNCH', 'PLOT', 'REAL', 'IMAG', 'PHASE',
        #'ABS', 'REL', 'PSDF', 'ATOC', 'CRMS', 'RMS', 'RALL', 'RPRINT',
        #'NOPRINT', 'RPUNCH',
    }
    allowed_strings = {'ALL', 'NONE'}
    allowed_values = {}  # type: Dict[str, Union[str, int]]
    allow_ints = True

    def __init__(self, key, value, options):
        CheckCard.__init__(self, key, value, options)

class STRFIELD(CheckCard):
    """
    STRFIELD=5
    STRFIELD(REAL)=ALL
    STRFIELD(SORT2, PUNCH, REAL)=ALL

    """
    type = 'STRFIELD'
    alternate_names = {'STRFIEL'}
    allowed_keys = {
        'SORT1', 'SORT2', 'PRINT', 'PUNCH', 'PLOT', 'REAL', 'IMAG', 'PHASE',
        #'ABS', 'REL', 'PSDF', 'ATOC', 'CRMS', 'RMS', 'RALL', 'RPRINT',
        #'NOPRINT', 'RPUNCH',
    }
    allowed_strings = {'ALL', 'NONE'}
    allowed_values = {}  # type: Dict[str, Union[str, int]]
    allow_ints = True

    def __init__(self, key, value, options):
        CheckCard.__init__(self, key, value, options)

class GPSDCON(CheckCard):
    """
    GPSDCON=5
    GPSDCON(REAL)=ALL
    GPSDCON(SORT2, PUNCH, REAL)=ALL

    """
    type = 'GPSDCON'
    #alternate_names = {'GPSD'}
    allowed_keys = {
        'SORT1', 'SORT2', 'PRINT', 'PUNCH', 'PLOT', 'REAL', 'IMAG', 'PHASE',
        #'ABS', 'REL', 'PSDF', 'ATOC', 'CRMS', 'RMS', 'RALL', 'RPRINT',
        #'NOPRINT', 'RPUNCH',
    }
    allowed_strings = {'ALL', 'NONE'}
    allowed_values = {}  # type: Dict[str, Union[str, int]]
    allow_ints = True

    def __init__(self, key, value, options):
        CheckCard.__init__(self, key, value, options)


class ELSDCON(CheckCard):
    """
    ELSDCON=5
    ELSDCON(REAL)=ALL
    ELSDCON(SORT2, PUNCH, REAL)=ALL

    """
    type = 'ELSDCON'
    alternate_names = {'ELSD'}
    allowed_keys = {
        'SORT1', 'SORT2', 'PRINT', 'PUNCH', 'PLOT', 'REAL', 'IMAG', 'PHASE',
        #'ABS', 'REL', 'PSDF', 'ATOC', 'CRMS', 'RMS', 'RALL', 'RPRINT',
        #'NOPRINT', 'RPUNCH',
    }
    allowed_strings = {'ALL', 'NONE'}
    allowed_values = {}  # type: Dict[str, Union[str, int]]
    allow_ints = True

    def __init__(self, key, value, options):
        CheckCard.__init__(self, key, value, options)

class OUTPUT(CheckCard):
    """
    OUTPUT
    OUTPUT(PLOT)
    OUTPUT(XYOUT)

    """
    # checked
    type = 'OUTPUT'
    allowed_keys = {
        'PLOT', 'POST', 'XYOUT', 'XYPLOT', 'CARDS',
    }
    allow_ints = True
    allow_equals = False

    def __init__(self, key, value, options):
        CheckCard.__init__(self, key, value, options)


class ENTHALPY(CheckCard):
    """
    ENTHALPY(PLOT) = ALL

    """
    type = 'ENTHALPY'
    allowed_keys = {
        'PLOT', 'PRINT',
        'SORT1',
    }
    allowed_strings = {'ALL'} # , 'NONE'}
    allow_ints = True

    def __init__(self, key, value, options):
        CheckCard.__init__(self, key, value, options)

class THERMAL(CheckCard):
    """
    THERMAL(PLOT) = ALL

    """
    type = 'THERMAL'
    allowed_keys = {
        'PLOT', 'PRINT', 'PUNCH',
        'SORT1', 'SORT2',
    }
    allowed_strings = {'ALL'} # , 'NONE'}
    allow_ints = True

    def __init__(self, key, value, options):
        CheckCard.__init__(self, key, value, options)

class TEMPERATURE(CheckCard):
    """
    TEMPERATURE(INITIAL) = 8
    TEMPERATURE(MATERIAL) = 8
    TEMPERATURE(LOAD) = 8
    TEMPERATURE(BOTH) = 8

    """
    type = 'TEMPERATURE'
    alternate_names = {'TEMP'}
    allowed_keys = {
        # estimate is in MSC examples, but is not documented
        'INITIAL', 'LOAD', 'MATERIAL', 'BOTH', 'ESTIMATE',
    }
    duplicate_names = {
        'INIT' : 'INITIAL',
        'MATE' : 'MATERIAL',
        'ESTI' : 'ESTIMATE',
    }
    #allowed_strings = {'ALL'} # , 'NONE'}
    allow_ints = True

    def __init__(self, key, value, options):
        CheckCard.__init__(self, key, value, options)

class FLUX(CheckCard):
    """
    FLUX(PLOT) = ALL

    """
    type = 'FLUX'
    allowed_keys = {
        'PLOT', 'PRINT',
        'SORT1', 'SORT2',
    }
    allowed_strings = {'ALL'} # , 'NONE'}
    allow_ints = True

    def __init__(self, key, value, options):
        CheckCard.__init__(self, key, value, options)

#class PRESSURE(CheckCard):
    #"""
    #PRESSURE(PLOT,PRINT) = ALL

    #"""
    #type = 'PRESSURE'
    ##short_name = 'PRESS'
    #alternate_names = {'PRESS'}

    #allowed_keys = {
        #'PLOT', 'PRINT',
    #}
    #allowed_strings = {'ALL'} # , 'NONE'}
    ##allow_ints = True

    #def __init__(self, key, value, options):
        #CheckCard.__init__(self, key, value, options)

class RESVEC(CheckCard):
    """
    RESVEC(APPLOD,RVDOF) = YES

    """
    type = 'RESVEC'
    short_name = 'RESVEC'
    allowed_keys = {
        'APPLOD', 'RVDOF', 'INRLOD',
        'NOAPPL', 'NOINER', 'NOINRL', 'NORVDO', 'NODYNRSP',
    }
    allowed_strings = {'YES', 'NO', 'NOCOMPONENT', 'NOSYSTEM', 'SYSTEM'} # , 'NONE'}
    allow_ints = True

    def __init__(self, key, value, options):
        CheckCard.__init__(self, key, value, options)

class ELSUM(CheckCard):
    """
    ELSUM(BOTH) = ALL
    ELSUM(EID,PIDSUM,NSMCONT) = ALL
    LSUM(EID,PID,BOTH,EIDSUM,PIDSUM,NSMCONT) = {ALL, NONE, N}
    """
    type = 'ELSUM'
    #short_name = 'RESVEC'
    allowed_keys = {'BOTH', 'EID', 'PID', 'BOTM', 'PIDSUM', 'NSMCONT',
                    # why is PRINT in the examples?
                    'PRINT'}
    allowed_strings = {'ALL', 'NONE'}
    allow_ints = True

    def __init__(self, key, value, options):
        CheckCard.__init__(self, key, value, options)



class DATAREC(CheckCard):
    """
    DATAREC = 11
    DATAREC(SORT1) = 11
    DATAREC(SORT2) = 11
    """
    type = 'DATAREC'
    allowed_keys = {'SORT1', 'SORT2'}
    allow_ints = True
    def __init__(self, key, value, options):
        CheckCard.__init__(self, key, value, options)

class BOUTPUT(CheckCard):
    """
    BOUTPUT(SORT2) = 2
    """
    type = 'BOUTPUT'
    allowed_keys = {'SORT1', 'SORT2'}
    allow_ints = True
    def __init__(self, key, value, options):
        CheckCard.__init__(self, key, value, options)

class BGRESULTS(CheckCard):
    """
    BGRESULTS(TRACTION,FORCE)=ALL
    """
    type = 'BGRESULTS'
    allowed_keys = {'TRACTION', 'FORCE'}
    #allow_ints = True
    #def __init__(self, key, value, options):
        #CheckCard.__init__(self, key, value, options)


class SENSITY(CheckCard):
    """
    SENSITY(STORE) = ALL

    """
    type = 'SENSITY'
    allowed_keys = {'STORE', }
    allowed_strings = {'ALL'} # , 'NONE'}
    #allowed_values = {}  # type: Dict[str, Union[str, int]]
    #alternate_names = {'PRES'}
    #allow_ints = True

    #def __init__(self, key, value, options):
        #CheckCard.__init__(self, key, value, options)

CHECK_CARDS = [
    DATAREC, DISPLACEMENT, VELOCITY, ACCELERATION, VECTOR,
    SDISPLACEMENT, SVELOCITY, SACCELERATION, SVECTOR,
    NLLOAD, NLSTRESS, OLOAD, OPRESS,
    OTEMP, SPCFORCES, MPCFORCES,
    EDE, ESE, STRESS, STRAIN, FORCE, GPFORCE, GPSTRESS, GPSTRAIN,
    ENTHALPY, THERMAL, FLUX, TEMPERATURE,
    STRFIELD, GPSDCON, ELSDCON,
    NOUTPUT, OUTPUT, RESVEC, ELSUM,
    BOUTPUT, BGRESULTS, PRESSURE, SENSITY,

]  # type: List[Any]

#-------------------------------------------------------------------------------
CHECK_CARD_DICT = {card.type : card for card in CHECK_CARDS} # type: Dict[str, str]
for card in CHECK_CARDS:
    if hasattr(card, 'alternate_names'):
        for name in card.alternate_names:
            CHECK_CARD_DICT[name] = card


# CHECK_CARD_NAMES = tuple([card.short_name for card in CHECK_CARDS])  # type: Tuple[str]
CHECK_CARD_NAMES = tuple(CHECK_CARD_DICT.keys())
