from __future__ import print_function, unicode_literals
from typing import List, Dict, Tuple, Union, Set, Any
from six import string_types, text_type
from pyNastran.bdf.bdf_interface.subcase_utils import write_set

class CaseControlCard(object):
    """basic card similar to the BaseCard class for the BDF"""
    def __iter__(self):
        """temporary method to emulate the old list access style"""
        value = self
        options = None
        param_type = 'OBJ-type'
        return iter([value, options, param_type])

#-------------------------------------------------------------------------------
class IntCard(CaseControlCard):
    """
    interface for cards of the form:
       NAME = 10

    """
    type = 'IntCard'
    def __init__(self, value):
        """
        Creates an IntCard

        Parameters
        ----------
        value : int
            the value for the card

        """
        super(IntCard, self).__init__()
        self.value = int(value)

    def __iter__(self):
        """temporary method to emulate the old list access style"""
        value = self
        options = []
        #param_type = 'STRESS-type'
        param_type = 'OBJ-type'
        return iter([value, options, param_type])

    @classmethod
    def add_from_case_control(cls, line, line_upper, lines, i):
        """
        Creates a card from the Case Control Deck

        Parameters
        ----------
        line : str
            the line of the card
        line_upper : str
            unused
        lines : List[str]
            unused
        i : int
            unused

        """
        value = line.split('=')[1]
        return cls(value)

    def __repr__(self):
        """writes a card"""
        return '%s = %i\n' % (self.type, self.value)

    def write(self, spaces):
        """writes a card with spaces"""
        return spaces + str(self)

class IntStrCard(IntCard):
    """
    interface for cards of the form:
       NAME = 10
       NAME = ALL

    """
    type = 'IntStrCard'
    allowed_strings = set([]) # type: Set[str]
    def __init__(self, value):
        """
        Creates an IntStrCard

        Parameters
        ----------
        value : int/str
            the value for the card

        """
        super(IntStrCard, self).__init__()
        try:
            self.value = int(value)
        except ValueError:
            if value not in self.allowed_strings:
                msg = 'value=%r not in [%s]' % (
                    self.value, ', '.join(self.allowed_strings))
                raise ValueError(msg)

    def __repr__(self):
        """writes a card"""
        return '%s = %s\n' % (self.type, self.value)


class ADACT(IntStrCard):
    type = 'ADACT'
    allowed_strings = {'ALL', 'NONE'}
    def __init__(self, value):
        super(ADACT, self).__init__(value)

class AEROF(IntStrCard):
    type = 'AEROF'
    allowed_strings = {'ALL'}
    def __init__(self, value):
        super(AEROF, self).__init__(value)

class APRES(IntStrCard):
    type = 'APRES'
    allowed_strings = {'ALL'}
    def __init__(self, value):
        super(APRES, self).__init__(value)

class GPRSORT(IntStrCard):
    type = 'GPRSORT'
    allowed_strings = {'ALL'}
    def __init__(self, value):
        super(GPRSORT, self).__init__(value)

class GPSDCON(IntStrCard):
    type = 'GPSDCON'
    allowed_strings = {'ALL'}
    def __init__(self, value):
        super(GPSDCON, self).__init__(value)

class HARMONICS(IntStrCard):
    type = 'HARMONICS'
    allowed_strings = {'ALL', 'NONE'}
    def __init__(self, value):
        super(HARMONICS, self).__init__(value)

class OFREQUENCY(IntStrCard):
    type = 'OFREQUENCY'
    allowed_strings = {'ALL'}
    def __init__(self, value):
        super(OFREQUENCY, self).__init__(value)

class OMODES(IntStrCard):
    type = 'OMODES'
    allowed_strings = {'ALL'}
    def __init__(self, value):
        super(OMODES, self).__init__(value)

INTSTR_CARDS = [
    ADACT, AEROF, APRES, GPRSORT, GPSDCON, HARMONICS, OFREQUENCY, OMODES]
INTSTR_CARD_DICT = {card.type : card for card in INTSTR_CARDS}
INTSTR_CARD_NAMES = tuple([card.type for card in INTSTR_CARDS])

#--------------------------------------------------------------
class AUXMODEL(IntCard):
    type = 'AUXMODEL'
    def __init__(self, value):
        super(AUXMODEL, self).__init__(value)

class BC(IntCard):
    type = 'BC'
    def __init__(self, value):
        super(BC, self).__init__(value)

class BCSET(IntCard):
    type = 'BCSET'
    def __init__(self, value):
        super(BCSET, self).__init__(value)

class BGSET(IntCard):
    type = 'BGSET'
    def __init__(self, value):
        super(BGSET, self).__init__(value)

class BOLTLD(IntCard):
    type = 'BOLTLD'
    def __init__(self, value):
        super(BOLTLD, self).__init__(value)

class CLOAD(IntCard):
    type = 'CLOAD'
    def __init__(self, value):
        super(CLOAD, self).__init__(value)

class CMETHOD(IntCard):
    type = 'CMETHOD'
    def __init__(self, value):
        super(CMETHOD, self).__init__(value)

class CSSCHD(IntCard):
    type = 'CSSCHD'
    def __init__(self, value):
        super(CSSCHD, self).__init__(value)

class DEFORM(IntCard):
    type = 'DEFORM'
    def __init__(self, value):
        super(DEFORM, self).__init__(value)

class DESGLB(IntCard):
    type = 'DESGLB'
    def __init__(self, value):
        super(DESGLB, self).__init__(value)

class DESSUB(IntCard):
    type = 'DESSUB'
    def __init__(self, value):
        super(DESSUB, self).__init__(value)

class DIVERG(IntCard):
    type = 'DIVERG'
    def __init__(self, value):
        super(DIVERG, self).__init__(value)

class DLOAD(IntCard):
    type = 'DLOAD'
    def __init__(self, value):
        super(DLOAD, self).__init__(value)

class DRSPAN(IntCard):
    type = 'DRSPAN'
    def __init__(self, value):
        super(DRSPAN, self).__init__(value)

class DTEMP(IntCard):
    type = 'DTEMP'
    def __init__(self, value):
        super(DTEMP, self).__init__(value)

class EBDSET(IntCard):
    type = 'EBDSET'
    def __init__(self, value):
        super(EBDSET, self).__init__(value)

class FMETHOD(IntCard):
    type = 'FMETHOD'
    def __init__(self, value):
        super(FMETHOD, self).__init__(value)

class FREQUENCY(IntCard):
    type = 'FREQUENCY'
    def __init__(self, value):
        super(FREQUENCY, self).__init__(value)

class GUST(IntCard):
    type = 'GUST'
    def __init__(self, value):
        super(GUST, self).__init__(value)

class LINE(IntCard):
    type = 'LINE'
    def __init__(self, value):
        super(LINE, self).__init__(value)

class LOAD(IntCard):
    type = 'LOAD'
    def __init__(self, value):
        super(LOAD, self).__init__(value)

class LOADSET(IntCard):
    type = 'LOADSET'
    def __init__(self, value):
        super(LOADSET, self).__init__(value)

class MAXLINES(IntCard):
    type = 'MAXLINES'
    def __init__(self, value):
        super(MAXLINES, self).__init__(value)

class MFLUID(IntCard):
    type = 'MFLUID'
    def __init__(self, value):
        super(MFLUID, self).__init__(value)

class MODES(IntCard):
    type = 'MODES'
    def __init__(self, value):
        super(MODES, self).__init__(value)

class MODTRAK(IntCard):
    type = 'MODTRAK'
    def __init__(self, value):
        super(MODTRAK, self).__init__(value)

class MPC(IntCard):
    type = 'MPC'
    def __init__(self, value):
        super(MPC, self).__init__(value)

class NLCNTL(IntCard):
    type = 'NLCNTL'
    def __init__(self, value):
        super(NLCNTL, self).__init__(value)

class NLPARM(IntCard):
    type = 'NLPARM'
    def __init__(self, value):
        super(NLPARM, self).__init__(value)

class NONLINEAR(IntCard):
    type = 'NONLINEAR'
    def __init__(self, value):
        super(NONLINEAR, self).__init__(value)

class NSM(IntCard):
    type = 'NSM'
    def __init__(self, value):
        super(NSM, self).__init__(value)

class OUTRCV(IntCard):
    type = 'OUTRCV'
    def __init__(self, value):
        super(OUTRCV, self).__init__(value)

class PARTN(IntCard):
    type = 'PARTN'
    def __init__(self, value):
        super(PARTN, self).__init__(value)

class REPCASE(IntCard):
    type = 'REPCASE'
    def __init__(self, value):
        super(REPCASE, self).__init__(value)

class RSMETHOD(IntCard):
    type = 'RSMETHOD'
    def __init__(self, value):
        super(RSMETHOD, self).__init__(value)

class SEFINAL(IntCard):
    type = 'SEFINAL'
    def __init__(self, value):
        super(SEFINAL, self).__init__(value)

class SMETHOD(IntCard):
    type = 'SMETHOD'
    def __init__(self, value):
        super(SMETHOD, self).__init__(value)

class SPC(IntCard):
    type = 'SPC'
    def __init__(self, value):
        super(SPC, self).__init__(value)

class SUPORT1(IntCard):
    type = 'SUPORT1'
    def __init__(self, value):
        super(SUPORT1, self).__init__(value)

class SYM(IntCard):
    type = 'SYM'
    def __init__(self, value):
        super(SYM, self).__init__(value)

class SYMCOM(IntCard):
    type = 'SYMCOM'
    def __init__(self, value):
        super(SYMCOM, self).__init__(value)

class TRIM(IntCard):
    type = 'TRIM'
    def __init__(self, value):
        super(TRIM, self).__init__(value)

class TSTEP(IntCard):
    type = 'TSTEP'
    def __init__(self, value):
        super(TSTEP, self).__init__(value)

class TSTEPNL(IntCard):
    type = 'TSTEPNL'
    def __init__(self, value):
        super(TSTEPNL, self).__init__(value)

class TSTRU(IntCard):
    type = 'TSTRU'
    def __init__(self, value):
        super(TSTRU, self).__init__(value)


INT_CARDS = [
    AUXMODEL, BC, BCSET, BGSET, BOLTLD, CLOAD, CMETHOD, CSSCHD,
    DEFORM, DESGLB, DESSUB, DIVERG, DLOAD, DRSPAN, DESSUB, DTEMP,
    EBDSET, FMETHOD, FREQUENCY, GUST, LINE, LOAD, LOADSET, MAXLINES,
    MFLUID, MODES, MODTRAK, MPC, NLCNTL, NLPARM, NONLINEAR, NSM,
    OUTRCV, PARTN, REPCASE, RSMETHOD, SEFINAL, SMETHOD, SPC,
    SUPORT1, SYM, SYMCOM, TRIM, TSTEP, TSTEPNL, TSTRU,
]
INT_CARD_DICT = {card.type : card for card in INT_CARDS}
INT_CARD_NAMES = tuple([card.type for card in INT_CARDS])

#-------------------------------------------------------------------------------

class StringCard(CaseControlCard):
    type = 'StringCard'
    allowed_values = [] # type: List[str]
    def __init__(self, value):
        super(StringCard, self).__init__()
        self.value = value.strip()
        self.validate()

    @classmethod
    def add_from_case_control(cls, line, line_upper, lines, i):
        """add method used by the CaseControl class"""
        value = line_upper.split('=')[1]
        return cls(value)

    def validate(self):
        if self.value not in self.allowed_values:
            msg = 'value=%r not in [%s]' % (self.value, ', '.join(self.allowed_values))
            raise ValueError(msg)

    def __repr__(self):
        """writes a card"""
        return '%s = %s\n' % (self.type, self.value)

    def write(self, spaces):
        return spaces + str(self)


class AESYMXY(StringCard):
    type = 'AESYMXY'
    allowed_values = ['SYMMETRIC', 'ANTISYMMETRIC', 'ASYMMTRIC']
    def __init__(self, value):
        super(AESYMXY, self).__init__(value)

class AESYMXZ(StringCard):
    type = 'AESYMXZ'
    allowed_values = ['SYMMETRIC', 'ANTISYMMETRIC', 'ASYMMTRIC']
    def __init__(self, value):
        super(AESYMXZ, self).__init__(value)

class AXISYMMETRIC(StringCard):
    type = 'AXISYMMETRIC'
    allowed_values = ['SINE', 'COSINE', 'FLUID']
    def __init__(self, value):
        super(AXISYMMETRIC, self).__init__(value)

class DSYM(StringCard):
    type = 'DSYM'
    allowed_values = ['S', 'A', 'SS', 'SA', 'AS', 'AA']
    def __init__(self, value):
        super(DSYM, self).__init__(value)

class SEQDEP(StringCard):
    type = 'SEQDEP'
    allowed_values = ['YES', 'NO']
    def __init__(self, value):
        super(SEQDEP, self).__init__(value)

STR_CARDS = [AESYMXY, AESYMXZ, AXISYMMETRIC, DSYM, SEQDEP]
STR_CARD_DICT = {card.type : card for card in STR_CARDS}
STR_CARD_NAMES = tuple([card.type for card in STR_CARDS])

#-------------------------------------------------------------------------------

class SET(CaseControlCard):
    type = 'SET'
    def __init__(self, set_id, values):
        super(SET, self).__init__()
        self.set_id = int(set_id)

        #values2 = expand_thru_case_control(values)
        self.value = values

    @property
    def key(self):
        """temporary method to emulate the old key attribute"""
        return '%s %s' % (self.type, self.set_id)

    def __iter__(self):
        """temporary method to emulate the old list access style"""
        value = self
        options = None
        param_type = 'OBJ-type'
        return iter([value, options, param_type])

    @classmethod
    def add_from_case_control(cls, line_upper, lines, i):
        """add method used by the CaseControl class"""
        line = lines[i]
        sline = line_upper.split('=')
        assert len(sline) == 2, sline

        key, value = sline
        try:
            (key, set_id) = key.split()
        except:
            raise RuntimeError(key)

        assert key.upper() == key, key
        uused_options = int(set_id)

        #if self.debug:
            #self.log.debug('SET-type key=%r set_id=%r' % (key, set_id))
        fivalues = value.rstrip(' ,').split(',')  # float/int values

        #: .. todo:: should be more efficient multiline reader...
        # read more lines....
        if line[-1].strip() == ',':
            i += 1
            #print("rawSETLine = %r" % (lines[i]))
            while 1:
                if lines[i].strip()[-1] == ',':
                    fivalues += lines[i][:-1].split(',')
                else:  # last case
                    fivalues += lines[i].split(',')
                    #print("fivalues last = i=%s %r" % (i, lines[i]))
                    i += 1
                    break
                i += 1
        #print("len(fivalues) = %s" % len(fivalues))
        return cls(set_id, fivalues)

    def write(self, spaces):
        """
        writes
        SET 80 = 3926, 3927, 3928, 4141, 4142, 4143, 4356, 4357, 4358, 4571,
             4572, 4573, 3323 THRU 3462, 3464 THRU 3603, 3605 THRU 3683,
             3910 THRU 3921, 4125 THRU 4136, 4340 THRU 4351

        """
        return write_set(self.set_id, self.value, spaces=spaces)

    def __repr__(self):
        """see `write`"""
        return write_set(self.set_id, self.value)

class SETMC(SET):
    """
    SETMC 121 = ACCE/99(T3),1200(T1),1399(R2)
    SETMC 222 = STRESS/134(22)
    SETMC 343 = ACCE/99(T3),1200(T1),1399(R2),STRESS/134(22)
    SETMC 122 = DISP/45(T1) 45(T2) 45(T3),
                 38(T1) 38(T2) 38(T3),
            VELO/45(T1) 45(T2) 45(T3),
                 38(T1) 38(T2) 38(T3),
            ACCE/45(T1) 45(T2) 45(T3),
                 38(T1) 38(T2) 38(T3)

    """
    type = 'SETMC'
    def __init__(self, set_id, values):
        super(SETMC, self).__init__(set_id, values)

class CheckCard(CaseControlCard):
    """
    Creates a card that validates the input

    GROUNDCHECK=YES
    GROUNDCHECK(GRID=12,SET=(G,N,A),THRESH=1.E-5,DATAREC=YES)=YES
    GROUNDCHECK(SET=ALL)=YES

    WEIGHTCHECK=YES
    WEIGHTCHECK(GRID=12,SET=(G,N,A),MASS)=YES
    WEIGHTCHECK(SET=ALL)=YES

    """
    type = 'CheckCard'
    allowed_keys = set([])  # type: Set[str]

    # key:(type, allowed_values)
    allowed_values = {}  # type: Dict[str, Union[float, str]]
    allowed_strings = set([]) # type: Set[str]
    duplicate_names = {} # type: Dict[Any, Any]
    allow_ints = False

    def export_to_hdf5(self, hdf5_file, encoding):
        #print(hdf5_file)
        #print('values* =', self.value)
        #print('options* =', self.options)
        if isinstance(self.options, list):
            options_bytes = [
                option.encode(encoding) if isinstance(option, text_type) else option
                for option in self.options]
            #print('optins =', options_bytes)
            hdf5_file.create_dataset('options', data=options_bytes)
        else:
            raise NotImplementedError(self.options)
        #else:
            #sub_group.create_dataset('options', data=self.options)

        if isinstance(self.data, list):
            data_group = hdf5_file.create_group('data')
            keys = []
            values = []
            for (key, value) in self.data:
                keys.append(key)
                values.append(value)
            #print('keys = ', keys)
            #print('values = ', values)
            keys_bytes = [
                key.encode(encoding) if isinstance(key, text_type) else key
                for key in keys]
            values_bytes = [
                value.encode(encoding) if isinstance(value, text_type) else value
                for value in values]
            data_group.create_dataset('keys', data=keys_bytes)
            data_group.create_dataset('values', data=values_bytes)
            #hdf5_file.create_dataset('data', data=data_bytes)
        else:
            raise NotImplementedError(self.data)

        hdf5_file.create_dataset('key', data=self.key)
        hdf5_file.create_dataset('value', data=self.value)
        #hdf5_file.create_dataset('options', data=self.options)
        #hdf5_file.create_dataset('data', data=self.data)

    def __init__(self, key, value, options):
        """
        Creates a card of the form:
            key(options) = value

        Parameters
        ----------
        key : str
            the name of the card
        value : List[str]
            the options
        value : str
            the response value

        """
        super(CheckCard, self).__init__()
        self.key = key
        self.options = options
        self.data = []
        for key_value in options:
            if key_value.upper().startswith('SET'):
                key = self._parse_set(key_value, options)
            else:
                key = self._parse(key_value, options)

            if key not in self.allowed_keys:
                msg = '%s: key=%r allowed_keys=[%r]' % (
                    self.type, key, ', '.join(self.allowed_keys))
                raise KeyError(msg)

        if isinstance(value, string_types):
            value = value.strip().upper()
        if self.allow_ints:
            try:
                value = int(value)
            except ValueError:
                if value not in self.allowed_strings:
                    msg = '%s: value=%r not in [%s]' % (
                        self.type, value, ', '.join(self.allowed_strings))
                    raise ValueError(msg)
        else:
            if value not in self.allowed_strings:
                msg = '%s: value=%r not in [%s]' % (
                    self.type, value, ', '.join(self.allowed_strings))
                raise ValueError(msg)
        self.value = value

    def _parse(self, key_value, options):
        if '=' in key_value:
            key, valuei = key_value.split('=')
            key = key.strip()
            valuei = valuei.strip()

            if key in self.duplicate_names:
                key = self.duplicate_names[key]

            if key in self.allowed_values:
                key_type, allowed_values = self.allowed_values[key]
                try:
                    valuei = key_type(valuei)
                except ValueError:
                    msg = 'cannot make %r a %s in %r' % (valuei, key_type, key_value)
                    raise ValueError(msg)
                except TypeError:
                    msg = 'cannot make %r a %s in %r' % (valuei, key_type, key_value)
                    raise TypeError(msg)

                # parse the value
                # SET=(G,N,A)
                if allowed_values is not None:
                    try:
                        sline = valuei.strip('(,)').split(',')
                    except AttributeError:
                        msg = 'cannot make %r a %s in %r of the form SET=(G,N,A)' % (
                            valuei, key_type, key_value)
                        raise ValueError(msg)

                    for val in sline:
                        if val not in allowed_values:
                            msg = '%s: key=%r value=%r allowed_values=[%r]' % (
                                self.type, key, val, ', '.join(allowed_values))
                            msg += '\noptions = %r' % options
                            raise ValueError(msg)

            key = key.upper()
            if isinstance(valuei, string_types):
                valuei = valuei.upper()
            self.data.append((key, valuei))
        else:
            key = key_value.upper()
            self.data.append((key, None))
        return key

    def _parse_set(self, key_value, options):
        """SET=(G,N,N+AUTOSPC,F,A)"""
        if '=' in key_value:
            key, valuei = key_value.split('=')
            key = key.strip()
            valuei = valuei.strip()

            if key in self.duplicate_names:
                key = self.duplicate_names[key]

            if key in self.allowed_values:
                key_type, allowed_values = self.allowed_values[key]
                try:
                    valuei = key_type(valuei)
                except ValueError:
                    msg = 'cannot make %r a %s in %r' % (valuei, key_type, key_value)
                    raise ValueError(msg)
                except TypeError:
                    msg = 'cannot make %r a %s in %r' % (valuei, key_type, key_value)
                    raise TypeError(msg)

                # parse the value
                # SET=(G,N,A)
                if allowed_values is not None:
                    try:
                        sline = valuei.strip('(,)').split(',')
                    except AttributeError:
                        msg = 'cannot make %r a %s in %r of the form SET=(G,N,A)' % (
                            valuei, key_type, key_value)
                        raise ValueError(msg)

                    for val in sline:
                        if '+' not in val or 'ALL' in val.upper():
                            # typical case
                            if val not in allowed_values:
                                msg = '%s: key=%r value=%r allowed_values=[%r]' % (
                                    self.type, key, val, ', '.join(allowed_values))
                                msg += '\noptions = %r' % options
                                raise ValueError(msg)
                        else:
                            vals = val.split('+')
                            # N+AUTOSPC
                            for vali in vals:
                                if vali not in allowed_values:
                                    msg = '%s: key=%r value=%r allowed_values=[%r]' % (
                                        self.type, key, val, ', '.join(allowed_values))
                                    msg += '\noptions = %r' % options
                                    raise ValueError(msg)

            self.data.append((key, valuei))
        else:
            key = key_value
            self.data.append((key, None))

        key = key.upper()
        return key

    @classmethod
    def add_from_case_control(cls, line, line_upper, lines, i):
        """add method used by the CaseControl class"""
        equals_count = line.count('=')
        if equals_count == 1:
            #GROUNDCHECK=YES
            #WEIGHTCHECK=YES

            if '=' in line:
                (key, value) = line_upper.strip().split('=')
            else:
                msg = 'expected item of form "name = value"   line=%r' % line.strip()
                raise RuntimeError(msg)

            key = key.strip().upper()
            value = value.strip()
            #if self.debug:
                #self.log.debug("key=%r value=%r" % (key, value))
            #param_type = 'STRESS-type'
            assert key.upper() == key, key

            if '(' in key:  # comma may be in line - STRESS-type
                #param_type = 'STRESS-type'
                sline = key.strip(')').split('(')
                key = sline[0]
                options = sline[1].split(',')

                # handle TEMPERATURE(INITIAL) and TEMPERATURE(LOAD) cards
                if key in ['TEMPERATURE', 'TEMP']:
                    option = options[0]
                    if option == '':
                        option = 'BOTH'
                    key = 'TEMPERATURE'
                    options = [option]
            else:
                # DISPLACEMENT = ALL
                options = []

        elif equals_count >= 2 and '(' in line:
            #GROUNDCHECK(PRINT,SET=(G,N,N+AUTOSPC,F,A),DATAREC=NO)=YES
            #WEIGHTCHECK(PRINT,SET=(G,N,F,A),CGI=NO,WEIGHT)=YES
            assert len(lines) == 1, lines
            line = lines[0]
            try:
                key, value_options = line.split('(', 1)
                #GROUNDCHECK, PRINT,SET=(G,N,N+AUTOSPC,F,A),DATAREC=NO)=YES
                #WEIGHTCHECK, PRINT,SET=(G,N,F,A),CGI=NO,WEIGHT)=YES
            except ValueError:
                msg = 'Expected a "(", but did not find one.\n'
                msg += 'Looking for something of the form:\n'
                msg += '   GROUNDCHECK(PRINT,SET=(G,N,N+AUTOSPC,F,A),DATAREC=NO)=YES\n'
                msg += '%r' % line
                raise ValueError(msg)

            try:
                options_paren, value = value_options.rsplit('=', 1)
                #'GROUNDCHECK', 'PRINT,SET=(G,N,N+AUTOSPC,F,A),DATAREC=NO)', 'YES'
                #'WEIGHTCHECK', 'PRINT,SET=(G,N,F,A),CGI=NO,WEIGHT)',        'YES'
            except ValueError:
                msg = 'Expected a "=", but did not find one.\n'
                msg += 'Looking for something of the form:\n'
                msg += '   GROUNDCHECK(PRINT,SET=(G,N,N+AUTOSPC,F,A),DATAREC=NO)=YES\n'
                msg += 'value_options=%r\n' % value_options
                msg += '%r' % line
                raise ValueError(msg)
            options_paren = options_paren.strip()

            value = value.strip()
            if value.isdigit():
                value = int(value)
            if not options_paren.endswith(')'):
                raise RuntimeError(line)
            str_options = options_paren[:-1]
            #'GROUNDCHECK', 'PRINT,SET=(G,N,N+AUTOSPC,F,A),DATAREC=NO', 'YES'
            #'WEIGHTCHECK', 'PRINT,SET=(G,N,F,A),CGI=NO,WEIGHT',        'YES'

            if '(' in str_options:
                options = split_by_mixed_commas_parentheses(str_options)
            else:
                options = str_options.split(',')
            #param_type = 'STRESS-type'
            key = key.upper()
        #elif equals_count == 2:
            #GROUNDCHECK(SET=ALL)=YES
            #WEIGHTCHECK(SET=ALL, PRINT, THRESH=0.01, DATAREC=NO)=YES
        else:
            raise RuntimeError('equals_count=%s; line = %r' % (equals_count, line))
        return cls(key, value, options)

    def write(self, spaces):
        msg = spaces + str(self)
        return msg

    def __repr__(self):
        """writes a card"""
        msg = '%s' % self.type
        if self.data:
            msg += '('
            for key, value in self.data:
                if value is None:
                    msg += '%s, ' % key
                else:
                    msg += '%s=%s, ' % (key, value)
            msg = msg.strip(', ') + ') = %s' % self.value
        else:
            msg += ' = %s' % self.value
        return msg + '\n'

def split_by_mixed_commas_parentheses(str_options):
    # type: (str) -> List[str]
    """
    Excessively complicated function to split something excessively
    complicated.  Thankfully, it only has one set of parentheses
    and no nested blocks.

    Parameters
    ----------
    str_options : str
        a nasty section of a case control line
        'PRINT,SET=(G,N,N+AUTOSPC,F,A),DATAREC=NO'
        'PRINT,SET=(G,N,F,A),CGI=NO,WEIGHT'

    Returns
    -------
    options : List[str]
        something that's actually parseable
        ['PRINT', 'SET=(G,N,N+AUTOSPC,F,A)', 'DATAREC=NO']
        ['PRINT', 'SET=(G,N,F,A)',           'CGI=NO',    'WEIGHT']

    """
    options_start = []
    options_end = []
    options_start_new = []  # type: List[str]
    options_end_new = []  # type: List[str]

    # search for ',' until one is '(' closer to the beginning
    # of the string; put it in options_start
    icomma = str_options.index(',')
    iparen = str_options.index('(')
    #print('icomma=%s iparen=%s' % (icomma, iparen))
    while icomma < iparen:
        base, str_options = str_options.split(',', 1)
        str_options = str_options.strip()
        icomma = str_options.index(',')
        iparen = str_options.index('(')
        options_start.append(base.strip())
        #print('  icomma=%s iparen=%s' % (icomma, iparen))
        #print('  options_start=%s' % options_start)

    # search for ',' until one is ')' closer to the end
    # of the string; put it in options_end
    icomma = str_options.rindex(',')
    iparen = str_options.rindex(')')
    #print('icomma=%s iparen=%s' % (icomma, iparen))
    while icomma > iparen:
        str_options, end = str_options.rsplit(')', 1)
        str_options = str_options.strip() + ')'
        iparen = str_options.rindex(')')
        if ',' in str_options:
            icomma = str_options.rindex(',')
        else:
            icomma = -1
        options_end.append(end.strip(' ,'))
        #print('  icomma=%s iparen=%s' % (icomma, iparen))
        #print('  options_end=%s' % options_end[::-1])

    #print()
    #print('options_start=%s' % options_start)
    #print('options_end=%s' % options_end)
    #print('leftover = %r' % str_options)

    # clean up the block and make sure we didn't mess up parsing the line
    for option in options_start:
        assert '(' not in option, option
        assert ')' not in option, option
        options_start_new += [optioni.strip() for optioni in option.split(',')]

    # we created options_end from right to left, so we need to reverse it
    for option in options_end[::-1]:
        assert '(' not in option, option
        assert ')' not in option, option
        options_end_new += [optioni.strip() for optioni in option.split(',')]

    options = options_start_new + [str_options] + options_end_new
    return options

class GROUNDCHECK(CheckCard):
    """
    GROUNDCHECK=YES
    GROUNDCHECK(GRID=12,SET=(G,N,A),THRESH=1.E-5,DATAREC=YES)=YES

    """
    type = 'GROUNDCHECK'
    allowed_keys = {'GRID', 'SET', 'PRINT', 'NOPRINT', 'THRESH', 'DATAREC', 'RTHRESH'}
    allowed_strings = {'YES'}
    allowed_values = {
        'CGI' : (str, ['YES', 'NO']),
        'SET' : (str, ['G', 'N', 'AUTOSPC', 'F', 'A', 'ALL']),
        'THRESH' : (float, None),
        'DATAREC' : (str, ['YES', 'NO']),
        'RTHRESH' : (float, None),
        'GRID' : (int, None),
    }

    def __init__(self, key, value, options):
        CheckCard.__init__(self, key, value, options)

class WEIGHTCHECK(CheckCard):
    """
    WEIGHTCHECK=YES
    WEIGHTCHECK(GRID=12,SET=(G,N,A),MASS)=YES

    """
    type = 'WEIGHTCHECK'
    allowed_keys = {'GRID', 'SET', 'PRINT', 'NOPRINT', 'CGI', 'MASS', 'WEIGHT'}
    allowed_strings = {'YES'}
    allowed_values = {
        'CGI': (str, ['YES', 'NO']),
        'SET': (str, ['G', 'N', 'AUTOSPC', 'F', 'A', 'V', 'ALL']),
        'GRID' : (int, None),
    }

    def __init__(self, key, value, options):
        CheckCard.__init__(self, key, value, options)

class MODCON(CheckCard):
    """
    MODCON=123
    MODCON(SORT1,PHASE,PRINT,PUNCH,BOTH,TOPS=5)=ALL

    """
    type = 'MODCON'
    allowed_keys = {'SORT1', 'SORT2', 'REAL', 'IMAG', 'PHASE', 'PRINT', 'NOPRINT',
                    'PUNCH', 'ABS', 'NORM', 'BOTH', 'TOPS', 'TOPF', 'SOLUTION',
                    'PANELMC'}
    duplicate_names = {
        'TOP' : 'TOPS',
        'SOLU' : 'SOLUTION',
        'PANE' : 'PANELMC',
    }
    allowed_strings = {'ALL', 'NONE'}
    allow_ints = True

    allowed_values = {
        'TOPS': (int, None),
        'TOPF': (int, None),
        'SOLUTION' : (int, None),  ## TODO: is this right???
    }

    def __init__(self, key, value, options):
        CheckCard.__init__(self, key, value, options)

class EXTSEOUT(CaseControlCard):
    """
    EXTSEOUT
    EXTSEOUT(ASMBULK,EXTID=100)
    EXTSEOUT(ASMBULK,EXTBULK,EXTID=200)
    EXTSEOUT(EXTBULK,EXTID=300)
    EXTSEOUT(DMIGDB)
    EXTSEOUT(ASMBULK,EXTID=400,DMIGOP2=21)
    EXTSEOUT(EXTID=500,DMIGPCH)
    EXTSEOUT(ASMBULK,EXTBULK,EXTID=500,DMIGSFIX=XSE500,DMIGPCH)
    EXTSEOUT(ASMBULK,EXTBULK,EXTID=500,DMIGSFIX=EXTID,DMIGPCH)
    EXTSEOUT(STIF,MASS,DAMP,EXTID=600,ASMBULK,EXTBULK,MATDB)
    EXTSEOUT(STIF,MASS,DAMP,GEOM,EXTID=600)

    """
    type = 'EXTSEOUT'
    allowed_keys = {'EXTID', 'ASMBULK', 'EXTBULK', 'MATDB', 'MATRIXDB',
                    'GEOM', 'DMIGSFIX', 'DMIGDB',
                    'STIFF', 'STIFFNESS', 'MASS',
                    'DAMP', 'DAMPING', 'K4DAMP',
                    'LOADS',
                    'DMIGOP2', 'DMIGPCH',
                    'MATOP4', 'MATRIXOP4'}

    def __init__(self, data):
        super(EXTSEOUT, self).__init__()
        self.data = data

    def export_to_hdf5(self, hdf5_file, encoding):
        if isinstance(self.data, list):
            data_group = hdf5_file.create_group('data')
            keys = []
            values = []
            for (key, value) in self.data:
                keys.append(key)
                values.append(value)
            #print('keys = ', keys)
            #print('values = ', values)
            keys_bytes = [
                key.encode(encoding) if isinstance(key, text_type) else key
                for key in keys]
            values_bytes = [
                value.encode(encoding) if isinstance(value, text_type) else value
                for value in values]
            data_group.create_dataset('keys', data=keys_bytes)

            if None in values_bytes:
                value_group = data_group.create_group('values')
                for i, value in enumerate(values):
                    if value is None:
                        continue
                    value_group.create_dataset(str(i), data=value)
            else:
                data_group.create_dataset('values', data=values_bytes)
            #hdf5_file.create_dataset('data', data=data_bytes)
        else:
            raise NotImplementedError(self.data)

    @classmethod
    def add_from_case_control(cls, line):
        """add method used by the CaseControl class"""
        data_list = []
        if '(' not in line:
            assert line == 'EXTSEOUT', line
        else:
            assert line.startswith('EXTSEOUT('), line
            assert line.endswith(')'), line
            data = line[9:-1].split(',')
            #print('data EXTSEOUT =', data)
            for key_value in data:
                key_value = key_value.strip()
                if '=' in key_value:
                    key, value = key_value.split('=')
                    key = cls._update_key(key)
                    value = value.strip()

                    data_list.append((key, value))
                else:
                    key = cls._update_key(key_value)
                    data_list.append((key, None))

                if key not in cls.allowed_keys:
                    msg = 'EXTSEOUT: key=%r allowed_keys=[%s]' % (key, ', '.join(cls.allowed_keys))
                    raise KeyError(msg)
        return EXTSEOUT(data_list)

    @staticmethod
    def _update_key(key):
        """
        STIFFNESS, DAMPING, K4DAMP, and LOADS may be abbreviated to STIF,
        DAMP, K4DA, and LOAD, respectively.

        """
        key = key.strip()
        if key == 'STIF':
            key = 'STIFFNESS'
        elif key == 'DAMP':
            key = 'DAMPING'
        elif key == 'K4DA':
            key = 'K4DAMP'
        elif key == 'LOAD':
            key = 'LOADS'
        return key

    def write(self, spaces):
        msg = spaces + str(self)
        return msg

    def __repr__(self):
        """writes a card"""
        msg = 'EXTSEOUT'
        if self.data:
            msg += '('
            for key, value in self.data:
                if value is None:
                    msg += '%s, ' % key
                else:
                    msg += '%s=%s, ' % (key, value)
            msg = msg.strip(', ') + ')'
        return msg + '\n'

#-------------------------------------------------------------------------------
class DISPLACEMENT(CheckCard):
    """
    DISPLACEMENT=5
    DISPLACEMENT(REAL)=ALL
    DISPLACEMENT(SORT2, PUNCH, REAL)=ALL

    """
    type = 'DISPLACEMENT'
    short_name = 'DISP'
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

class VELOCITY(CheckCard):
    """
    VELOCITY=5
    VELOCITY(REAL)=ALL
    VELOCITY(SORT2, PUNCH, REAL)=ALL

    """
    type = 'VELOCITY'
    short_name = 'VELO'

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

class ACCELERATION(CheckCard):
    """
    ACCELERATION=5
    ACCELERATION(REAL)=ALL
    ACCELERATION(SORT2, PUNCH, REAL)=ALL

    """
    type = 'ACCELERATION'
    short_name = 'ACCE'
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

class SPCFORCES(CheckCard):
    """
    SPCFORCES=5
    SPCFORCES(REAL)=ALL
    SPCFORCES(SORT2, PUNCH, REAL)=ALL

    """
    type = 'SPCFORCES'
    short_name = 'SPCF'
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
    short_name = 'MPCF'
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
    short_name = type
    allowed_keys = {'PRINT', 'PUNCH'}
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
    short_name = type
    allowed_keys = {'SORT1', 'SORT2', 'PRINT', 'PUNCH', 'PLOT'}
    allowed_strings = {'ALL', 'NONE'}
    allowed_values = {}  # type: Dict[str, Union[str, int]]
    allow_ints = True

    def __init__(self, key, value, options):
        CheckCard.__init__(self, key, value, options)

#class NOUTPUT(CheckCard):
    #"""
    #NOUTPUT (R)=ALL
    #NOUTPUT (2)=5
    #NOUTPUT (4,L)=10
    #"""
    #type = 'NOUTPUT'
    #allowed_keys = {'SORT1', 'SORT2', 'PRINT', 'PUNCH', 'PLOT']
    #allowed_strings = {'ALL'}
    #allowed_values = {}

    #def __init__(self, key, value, options):
        #super(NLLOAD, self).__init__(key, value, options)

class OLOAD(CheckCard):
    """
    OLOAD=ALL
    OLOAD(SORT1,PHASE)=5

    """
    type = 'OLOAD'
    short_name = type
    allowed_keys = {'SORT1', 'SORT2', 'PRINT', 'PUNCH', 'PLOT',
                    'REAL', 'IMAG', 'PHASE', 'PSDF', 'ATOC', 'CRMS',
                    'RMS', 'RALL', 'RPRINT', 'NORPRINT', 'RPUNCH'}
    allowed_strings = {'ALL', 'NONE'}
    allowed_values = {}  # type: Dict[str, Union[str, int]]
    allow_ints = True

    def __init__(self, key, value, options):
        CheckCard.__init__(self, key, value, options)


class OPRESS(CheckCard):
    """
    OPRESS=ALL
    OPRESS(PRINT,PUNCH)=17

    """
    type = 'OPRESS'
    short_name = type
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
    short_name = type
    allowed_keys = {'PRINT', 'PUNCH'}
    allowed_strings = {'ALL', 'NONE'}
    allowed_values = {}  # type: Dict[str, Union[str, int]]
    allow_ints = True

    def __init__(self, key, value, options):
        CheckCard.__init__(self, key, value, options)

CHECK_CARDS = [
    DISPLACEMENT, VELOCITY, ACCELERATION, NLLOAD, NLSTRESS, OLOAD, OPRESS,
    OTEMP, SPCFORCES, MPCFORCES,
]  # type: List[Any]
CHECK_CARD_DICT = {card.type : card for card in CHECK_CARDS} # type: Dict[str, str]
CHECK_CARD_NAMES = tuple([card.short_name for card in CHECK_CARDS])  # type: Tuple[str]

#-------------------------------------------------------------------------------

CLASS_MAP = {
    'GROUNDCHECK' : GROUNDCHECK,
    'EXTSEOUT' : EXTSEOUT,
    'WEIGHTCHECK' : WEIGHTCHECK,
    'MODCON' : MODCON,
    'SET' : SET,
    'SETMC' : SETMC,
}
