from typing import List, Dict, Union, Set, Any
from pyNastran.bdf.bdf_interface.subcase_utils import write_set
from pyNastran.utils.numpy_utils import bytes_type

def decode_bytes_list(bytes_list, encoding):
    out = [bytes_str.decode(encoding) if isinstance(bytes_str, bytes_type) else bytes_str
           for bytes_str in bytes_list]
    return out

class CaseControlCard:
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
        value = line_upper.split('=')[1]
        try:
            out = cls(value)
        except ValueError:
            print(line)
            raise
        return out

    def export_to_hdf5(self, h5_file, encoding):
        h5_file.create_dataset('value', data=self.value)

    @classmethod
    def load_hdf5(cls, h5_file, encoding):
        from pyNastran.utils.dict_to_h5py import _cast
        value = h5_file['value']
        value2 = _cast(value)
        return cls(value2), []

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
        #super(IntStrCard, self).__init__()
        try:
            self.value = int(value)
        except ValueError:
            value = value.strip()
            if value not in self.allowed_strings:
                msg = 'value=%r not in [%s]' % (
                    value, ', '.join(self.allowed_strings))
                raise ValueError(msg)
            self.value = value

    def export_to_hdf5(self, h5_file, encoding):
        value_bytes = self.value.encode(encoding) if isinstance(self.value, str) else self.value
        #sub_group = h5_file.create_group(self.type)
        h5_file.create_dataset('value', data=value_bytes)

    @classmethod
    def load_hdf5(cls, h5_file, encoding):
        from pyNastran.utils.dict_to_h5py import _cast
        value = h5_file['value']

        casted_value = _cast(value)
        if isinstance(casted_value, int):
            value2 = casted_value
        else:
            value2 = casted_value.decode(encoding)# if isinstance(value, bytes) else value
        return cls(value2), []

    def __repr__(self):
        """writes a card"""
        return '%s = %s\n' % (self.type, self.value)


class ADACT(IntStrCard):
    type = 'ADACT'
    allowed_strings = {'ALL', 'NONE'}
    def __init__(self, value):
        super().__init__(value)

class AEROF(IntStrCard):
    type = 'AEROF'
    allowed_strings = {'ALL'}
    def __init__(self, value):
        super().__init__(value)

class APRES(IntStrCard):
    type = 'APRES'
    allowed_strings = {'ALL'}
    def __init__(self, value):
        super().__init__(value)

class GPRSORT(IntStrCard):
    type = 'GPRSORT'
    allowed_strings = {'ALL'}
    def __init__(self, value):
        super().__init__(value)

class GPSDCON(IntStrCard):
    type = 'GPSDCON'
    allowed_strings = {'ALL'}
    def __init__(self, value):
        super().__init__(value)

class HARMONICS(IntStrCard):
    type = 'HARMONICS'
    allowed_strings = {'ALL', 'NONE'}
    def __init__(self, value):
        super().__init__(value)

class OFREQUENCY(IntStrCard):
    type = 'OFREQUENCY'
    alternate_names = {'OFREQ'}
    allowed_strings = {'ALL'}
    def __init__(self, value):
        super().__init__(value)

class OMODES(IntStrCard):
    type = 'OMODES'
    allowed_strings = {'ALL'}
    def __init__(self, value):
        super().__init__(value)

class SUPER(IntStrCard):
    type = 'SUPER'
    allowed_strings = {'ALL'}
    #def __init__(self, value):
        #super().__init__(value)

#----------------------------

class GPSTRESS(IntStrCard):
    type = 'GPSTRESS'
    allowed_strings = {'ALL'}
    def __init__(self, value):
        super().__init__(value)

class SEALL(IntStrCard):
    type = 'SEALL'
    allowed_strings = {'ALL'}
    def __init__(self, value):
        super().__init__(value)

class SEDR(IntStrCard):
    type = 'SEDR'
    allowed_strings = {'ALL'}
    def __init__(self, value):
        super().__init__(value)

class GPKE(IntStrCard):
    type = 'GPKE'
    allowed_strings = {'ALL'}
    def __init__(self, value):
        super().__init__(value)

INTSTR_CARDS = [
    ADACT, AEROF, APRES, GPRSORT, GPSDCON, HARMONICS, OFREQUENCY, OMODES,
    SUPER, SEALL, SEDR,
] + [GPSTRESS, GPKE, ]
INTSTR_CARD_DICT = {card.type : card for card in INTSTR_CARDS}
INTSTR_CARD_NAMES = tuple([card.type for card in INTSTR_CARDS])

#-------------------------------------------------------------------------------

class StringCard(CaseControlCard):
    type = 'StringCard'
    allowed_values = [] # type: List[str]
    def __init__(self, value, validate=True):
        super(StringCard, self).__init__()
        self.value = value.strip()
        if validate:
            self.validate()

    @classmethod
    def add_from_case_control(cls, line, line_upper, lines, i):
        """add method used by the CaseControl class"""
        value = line_upper.split('=')[1]
        return cls(value)

    def validate(self):
        if self.value not in self.allowed_values:
            msg = '%s: value=%r not in [%s]' % (
                self.type, self.value, ', '.join(self.allowed_values))
            raise ValueError(msg)

    def __repr__(self):
        """writes a card"""
        return '%s = %s\n' % (self.type, self.value)

    def write(self, spaces):
        return spaces + str(self)

    def export_to_hdf5(self, h5_file, encoding):
        value_bytes = self.value.encode(encoding)
        #sub_group = h5_file.create_group(self.type)
        h5_file.create_dataset('value', data=value_bytes)

    @classmethod
    def load_hdf5(cls, h5_file, encoding):
        from pyNastran.utils.dict_to_h5py import _cast
        value = h5_file['value']
        try:
            value2 = _cast(value).decode(encoding)
        except AttributeError:
            print(cls.type, _cast(value))
            raise
        return cls(value2), []

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
        unused_options = int(set_id)

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

    # the allowed value for the key, options, value approach
    allowed_strings = set([]) # type: Set[str]

    # maps something like INIT to INITIAL
    duplicate_names = {} # type: Dict[Any, Any]

    # enables values as integers instead of just strings
    allow_ints = False

    # 'OUTPUT(PLOT)' instead of 'OUTPUT(PLOT)=YES'
    allow_equals = True

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

            if key in self.duplicate_names:
                key = self.duplicate_names[key]

            if key not in self.allowed_keys:
                msg = '%s: key=%r allowed_keys=[%r]' % (
                    self.type, key, ', '.join(self.allowed_keys))
                raise KeyError(msg)

        if isinstance(value, str):
            value = value.strip().upper()
        if self.allow_equals:
            if self.allow_ints:
                try:
                    value = int(value)
                except ValueError:
                    if value not in self.allowed_strings:
                        msg = '%s: value=%r not in [%s]; options=%s' % (
                            self.type, value, ', '.join(self.allowed_strings), options)
                        raise ValueError(msg)
            else:
                if value not in self.allowed_strings:
                    msg = '%s: value=%r not in [%s]; options=%s' % (
                        self.type, value, ', '.join(self.allowed_strings), options)
                    raise ValueError(msg)
        else:
            assert value is None, value
        self.value = value

    def _parse(self, key_value, options):
        if '=' in key_value:
            assert self.allow_equals is True, key_value
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
            if isinstance(valuei, str):
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
        if cls.allow_equals:
            if equals_count == 1:
                line_upper = line_upper.replace(' ', '')
                # GROUNDCHECK=YES
                # WEIGHTCHECK=YES
                key, value, options = cls._parse_single_equals(line, line_upper, lines, i)

            elif equals_count >= 2 and '(' in line:
                #GROUNDCHECK(PRINT,SET=(G,N,N+AUTOSPC,F,A),DATAREC=NO)=YES
                #WEIGHTCHECK(PRINT,SET=(G,N,F,A),CGI=NO,WEIGHT)=YES
                key, value, options = cls._parse_multi_equals(line, line_upper, lines, i)

            #elif equals_count == 2:
                #GROUNDCHECK(SET=ALL)=YES
                #WEIGHTCHECK(SET=ALL, PRINT, THRESH=0.01, DATAREC=NO)=YES
            else:
                raise RuntimeError('equals_count=%s; line = %r' % (equals_count, line))
        else:
            value = None
            if '(' in line_upper:
                (class_name, options_str) = line_upper.strip(')').split('(')
                options = options_str.split(',')
            else:
                class_name = line_upper
                assert class_name == 'OUTPUT', class_name
                options = []
            key = class_name
            #print(f'options_str = {options_str!r}')
        return cls(key, value, options)

    @classmethod
    def _parse_single_equals(cls, line, line_upper, lines, i):
        """
        GROUNDCHECK=YES
        WEIGHTCHECK=YES
        """
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
        return key, value, options

    @classmethod
    def _parse_multi_equals(cls, line, line_upper, lines, i):
        """
        #GROUNDCHECK(PRINT,SET=(G,N,N+AUTOSPC,F,A),DATAREC=NO)=YES
        #WEIGHTCHECK(PRINT,SET=(G,N,F,A),CGI=NO,WEIGHT)=YES
        """
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
        return key, value, options

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

def split_by_mixed_commas_parentheses(str_options: str) -> List[str]:
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

    def export_to_hdf5(self, hdf5_file, encoding):
        export_to_hdf5_check(self, hdf5_file, encoding)


class MEFFMASS(CheckCard):
    """
    MEFFMASS
    MEFFMASS(GRID=12,SUMMARY,PARTFAC)
    MEFFMASS(PLOT,ALL,THRESH=0.001)=YES

    """
    type = 'MEFFMASS'
    allowed_keys = {
        'PRINT', 'PLOT', 'PUNCH',
        'MINT1', 'MINT2', 'MINT3', 'MAXIT',
        'THRESH', 'GRID',
        'SUMMARY', 'PARTFAC', 'MEFFM', 'MEFFW', 'FRACSUM', 'ALL'}
    allowed_strings = {'YES', 'NO'}
    alternate_names = {'MEFF'}
    allowed_values = {
        'GRID' : (int, None),
        'MINT1' : (int, None),
        'MINT2' : (int, None),
        'MINT3' : (int, None),
        'MAXIT' : (int, None),
        'THRESH' : (float, None),
    }  # type: Dict[str, Union[str, int]]
    #alternate_names = {'PRES'}
    #allow_ints = True

    #def __init__(self, key, value, options):
        #CheckCard.__init__(self, key, value, options)

    def export_to_hdf5(self, hdf5_file, encoding):
        export_to_hdf5_check(self, hdf5_file, encoding)


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

    def export_to_hdf5(self, hdf5_file, encoding):
        export_to_hdf5_check(self, hdf5_file, encoding)

class DSAPRT(CheckCard):
    """
    DSAPRT(END=SENS)=ALL
    DSAPRT(FORMATTED,EXPORT)
    DSAPRT(FORMATTED,START=FIRST,BY=3,END=LAST)=101
    DSAPRT(UNFORMATTED,START=FIRST)
    DSAPRT(UNFORMATTED,EXPORT)
    DSAPRT(FORMATTED,END=4)=ALL
    DSAPRT(UNFORMATTED,END=SENS)=ALL
    DSAPRT(NOPRINT, EXPORT)
    """
    # not done...
    type = 'DSAPRT'
    allowed_keys = {
        'FORMATTED', 'UNFORMATTED', 'EXPORT',
        'START', 'BY', 'END',}
    allowed_strings = {'ALL'}
    allowed_values = {
        'START' : (str, ['FIRST']),
        'BY' : (int, None),
        'END' : (str, ['SENS', 'LAST']),
    }

    def __init__(self, key, value, options):
        CheckCard.__init__(self, key, value, options)

    def export_to_hdf5(self, hdf5_file, encoding):
        export_to_hdf5_check(self, hdf5_file, encoding)


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
            values = []
            keys = []
            keys_none = []
            for (key, value) in self.data:
                if value is None:
                    keys_none.append(key)
                else:
                    keys.append(key)
                    values.append(value)

            if keys_none:
                keys_none_bytes = [key.encode(encoding) for key in keys_none]
                data_group.create_dataset('keys_none', data=keys_none_bytes)
            if keys:
                keys_bytes = [key.encode(encoding) for key in keys]
                values_bytes = [value.encode(encoding) if isinstance(value, str) else value
                                for value in values]
                data_group.create_dataset('keys', data=keys_bytes)
                data_group.create_dataset('values', data=values_bytes)
        else:
            raise NotImplementedError(self.data)

    @classmethod
    def load_hdf5(self, subgroup, encoding):
        """loads EXTSEOUT from an h5py HDF5 file"""
        from pyNastran.utils.dict_to_h5py import _cast
        for key in subgroup.keys():
            subgroupi = subgroup[key]
            if key == 'data':
                data_keys = _cast(subgroupi['keys'])
                data_keys = decode_bytes_list(data_keys, encoding)

                keys_none = subgroupi['keys_none']
                keys_none = decode_bytes_list(keys_none, encoding)

                data_values = _cast(subgroupi['values'])
                data_values = decode_bytes_list(data_values, encoding)

                keys = data_keys + keys_none
                values = data_values + [None] * len(keys_none)
                data = [(key, value) for (key, value) in zip(keys, values)]
            else:
                raise NotImplementedError(key)
        return EXTSEOUT(data), []

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

class VOLUME(CaseControlCard):
    """
    VOLUME 21 SET 2
    VOLUME id SET sid, [PRINCIPAL, DIRECT STRESS] [SYSTEM {ELEMENT, CORD cid, BASIC}]
    """
    type = 'VOLUME'

    def __init__(self, data):
        super().__init__()
        self.data = data

    #def export_to_hdf5(self, hdf5_file, encoding):
        #if isinstance(self.data, list):
            #data_group = hdf5_file.create_group('data')
            #keys = []
            #values = []
            #for (key, value) in self.data:
                #keys.append(key)
                #values.append(value)
            ##print('keys = ', keys)
            ##print('values = ', values)
            #keys_bytes = [
                #key.encode(encoding) if isinstance(key, str) else key
                #for key in keys]
            #values_bytes = [
                #value.encode(encoding) if isinstance(value, str) else value
                #for value in values]
            #data_group.create_dataset('keys', data=keys_bytes)

            #if None in values_bytes:
                #value_group = data_group.create_group('values')
                #for i, value in enumerate(values):
                    #if value is None:
                        #continue
                    #value_group.create_dataset(str(i), data=value)
            #else:
                #data_group.create_dataset('values', data=values_bytes)
            ##hdf5_file.create_dataset('data', data=data_bytes)
        #else:
            #raise NotImplementedError(self.data)

    @classmethod
    def add_from_case_control(cls, line):
        """add method used by the CaseControl class"""
        sline = line.split()

        i = 0
        data = {}
        while i < len(sline):
            word = sline[i]
            if word == 'VOLUME':
                value = sline[i+1]
                data[word] = int(value)
                i += 2
            elif word == 'SET':
                value = sline[i+1]
                data[word] = int(value)
                i += 2
            elif word == 'DIRECT':
                # this is confusing...
                value = 'NONE'
                if i+1 < len(sline):
                    #print(sline)
                    value = sline[i+1]
                #data[word] = None
                data[word] = value
                assert value in ['NONE'], 'DIRECT value=%r' % value
                i += 2
            else:
                raise RuntimeError(f'VOLUME: {word}; {sline}\n{line}')
        return VOLUME(data)

    def write(self, spaces):
        msg = spaces + str(self)
        return msg

    def __repr__(self):
        """writes a card"""
        msg = 'VOLUME %i' % self.data['VOLUME']
        for key, value in sorted(self.data.items()):
            if key == 'VOLUME':
                continue
            msg += ' %s %s' % (key, value)
        return msg + '\n'

class SURFACE(CaseControlCard):
    """
    SURFACE 10 SET 9 NORMAL X3
    SURFACE 41 SET 42 FIBRE ALL NORMAL Z
    VOLUME id SET sid, [PRINCIPAL, DIRECT STRESS] [SYSTEM {ELEMENT, CORD cid, BASIC}]

    """
    type = 'VOLUME'

    def __init__(self, data):
        super().__init__()
        self.data = data

    #def export_to_hdf5(self, hdf5_file, encoding):
        #if isinstance(self.data, list):
            #data_group = hdf5_file.create_group('data')
            #keys = []
            #values = []
            #for (key, value) in self.data:
                #keys.append(key)
                #values.append(value)
            ##print('keys = ', keys)
            ##print('values = ', values)
            #keys_bytes = [
                #key.encode(encoding) if isinstance(key, str) else key
                #for key in keys]
            #values_bytes = [
                #value.encode(encoding) if isinstance(value, str) else value
                #for value in values]
            #data_group.create_dataset('keys', data=keys_bytes)

            #if None in values_bytes:
                #value_group = data_group.create_group('values')
                #for i, value in enumerate(values):
                    #if value is None:
                        #continue
                    #value_group.create_dataset(str(i), data=value)
            #else:
                #data_group.create_dataset('values', data=values_bytes)
            ##hdf5_file.create_dataset('data', data=data_bytes)
        #else:
            #raise NotImplementedError(self.data)

    @classmethod
    def add_from_case_control(cls, line):
        """add method used by the CaseControl class"""
        sline = line.split()

        i = 0
        data = {}
        while i < len(sline):
            word = sline[i]
            if word == 'SURFACE':
                value = sline[i+1]
                data[word] = int(value)
                i += 2
            elif word == 'SET':
                value = sline[i+1]
                data[word] = int(value)
                i += 2
            elif word == 'FIBRE':
                value = sline[i+1]
                data[word] = value
                assert value in  ['ALL'], 'SURFACE: %r=%r' % (word, value)
                i += 2
            elif word == 'NORMAL':
                value = sline[i+1]
                data[word] = value
                assert value in  ['Z', 'X1', 'X2', 'X3'], 'SURFACE: %r=%r' % (word, value)
                i += 2
            else:
                raise RuntimeError(word)
        return SURFACE(data)

    def write(self, spaces):
        msg = spaces + str(self)
        return msg

    def __repr__(self):
        """writes a card"""
        msg = 'SURFACE %i' % self.data['SURFACE']
        for key, value in sorted(self.data.items()):
            if key == 'SURFACE':
                continue
            msg += ' %s %s' % (key, value)
        return msg + '\n'


class CSCALE(CaseControlCard):
    """
    CSCALE 1.3

    """
    type = 'CSCALE'

    def __init__(self, data):
        super().__init__()
        self.data = data

    #def export_to_hdf5(self, hdf5_file, encoding):
        #if isinstance(self.data, list):
            #data_group = hdf5_file.create_group('data')
            #keys = []
            #values = []
            #for (key, value) in self.data:
                #keys.append(key)
                #values.append(value)
            ##print('keys = ', keys)
            ##print('values = ', values)
            #keys_bytes = [
                #key.encode(encoding) if isinstance(key, str) else key
                #for key in keys]
            #values_bytes = [
                #value.encode(encoding) if isinstance(value, str) else value
                #for value in values]
            #data_group.create_dataset('keys', data=keys_bytes)

            #if None in values_bytes:
                #value_group = data_group.create_group('values')
                #for i, value in enumerate(values):
                    #if value is None:
                        #continue
                    #value_group.create_dataset(str(i), data=value)
            #else:
                #data_group.create_dataset('values', data=values_bytes)
            ##hdf5_file.create_dataset('data', data=data_bytes)
        #else:
            #raise NotImplementedError(self.data)

    @classmethod
    def add_from_case_control(cls, line):
        """add method used by the CaseControl class"""
        sline = line.split()
        value = float(sline[1])
        return CSCALE(value)

    def write(self, spaces):
        msg = spaces + str(self)
        return msg

    def __repr__(self):
        """writes a card"""
        msg = 'CSCALE %s\n' % self.value
        return msg

def export_to_hdf5_check(self, hdf5_file, encoding):
    #print(hdf5_file)
    #print('values* =', self.value)
    #print('options* =', self.options)

    if isinstance(self.options, list):
        options_bytes = [
            option.encode(encoding) if isinstance(option, str) else option
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
            key.encode(encoding) if isinstance(key, str) else key
            for key in keys]
        values_bytes = [
            value.encode(encoding) if isinstance(value, str) else value
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


#-------------------------------------------------------------------------------

