from typing import Union, Any
from pyNastran.utils.numpy_utils import bytes_type
from .subcase_base import CaseControlCard

from .int_str_cards import (
    AEROF, APRES, GPKE, GPRSORT, GPSDCON, GPSTRESS,
    HARMONICS, OFREQUENCY, OMODES, SEALL, SEDR, SUPER)
from .matrix import (
    A2GG,
    B2GG, M2GG, K2GG, P2G, K42GG,
    B2PP, M2PP, K2PP, )
from .sets import SET, SETMC
#from pyNastran.bdf.field_writer_8 import print_float_8

def encode_str_list(strings: list[str], encoding: str) -> list[bytes]:
    return [stri.encode(encoding) for stri in strings]

def encode_str_value_list(strings: list[Union[int, float, str]],
                          encoding: str) -> list[Union[int, float, bytes]]:
    values_bytes = [stri.encode(encoding) if isinstance(stri, str) else stri
                    for stri in strings]
    return values_bytes

def decode_bytes_list(bytes_list: list[bytes], encoding: str) -> list[str]:
    out = [bytes_str.decode(encoding) if isinstance(bytes_str, bytes_type) else bytes_str
           for bytes_str in bytes_list]
    return out

#-------------------------------------------------------------------------------

class ECHO(CaseControlCard):
    """
    ECHO = NONE
    ECHO = NOSORT
    ECHO = BOTH

    ECHO = SORT(PARAM,EIGC,EIGRL,FREQ,DESVAR,DCONSTR,DRESP1,DRESP2,DEQATN,DVPREL1)
    ECHO = PUNCH,SORT(MAT1,PARAM)
    ECHO = SORT(EXCEPT DMI,DMIG)
    ECHO = PUNCH(BSTBULK)
    ECHO = PUNCH(NEWBULK)
    ECHO = SORT,PUNCH(BSTBULK)
    """
    type = 'ECHO'
    def __init__(self, value: str):
        self.value = value

    def __repr__(self) -> str:
        """writes a card"""
        rows = self._repr_rows()
        return '\n'.join(rows)

    def write(self, spaces: str):
        """writes a card with spaces"""
        rows = self._repr_rows()
        out = spaces + ('\n' + spaces).join(rows) + '\n'
        return out

    def _repr_rows(self) -> list[str]:
        """writes a card"""
        target = 68
        #['SORT', ['EXCEPT DMI', 'DMIG'], 'PUNCH', ['BSTBULK']]
        rows = []
        msg = 'ECHO = '
        base = '       '
        for valuei in self.value:
            if isinstance(valuei, str):
                new = valuei + ','
                if len(msg) + len(new) < target:
                    msg += new
                else:
                    rows.append(msg)
                    msg = base + new
            elif isinstance(valuei, list):
                #msg =  + '(%s),' % ','.join(valuei)  #  too long
                msg = msg[:-1] + '('
                for valueii in valuei:
                    new = valueii + ','
                    if len(msg) + len(new) < target:
                        msg += new
                    else:
                        rows.append(msg)
                        msg = base + new
                msg = msg[:-1] + '),'
            else:
                raise TypeError(self.value)

        if msg:
            rows.append(msg[:-1])
        else:
            raise RuntimeError(msg)

        for row in rows:
            assert len(row) < target, f'row={row!r}'
        return rows

    @classmethod
    def add_from_case_control(cls, line: str):
        """
        Creates a card from the Case Control Deck

        Parameters
        ----------
        line : str
            the line of the card

        """
        value = line.split('=')[1].strip() #  NONE
        if value in {'NONE', 'NOSORT', 'BOTH', 'UNSORT'}:
            # ECHO = NONE
            # ECHO = NOSORT
            # ECHO = BOTH
            return cls([value])
        options = []
        _set_options_from_line(line, value, options)

        options2 = _strip_echo_options(options)
        return cls(options2)

    def export_to_hdf5(self, h5_file, encoding: str) -> None:
        values = []
        for value in self.value:
            if isinstance(value, str):
                values.append(value)
            else:
                valuei = '(%s)' % ','.join(value)
                values.append(valuei)
        values_bytes = encode_str_list(values, encoding)
        h5_file.create_dataset('value', data=values_bytes)

    @classmethod
    def load_hdf5(cls, h5_file, encoding: str) -> Any:
        from pyNastran.utils.dict_to_h5py import _cast
        values_bytes = _cast(h5_file['value'])
        values_str = decode_bytes_list(values_bytes, encoding)

        values = []
        for value in values_str:
            if value.startswith('('):
                #print(values_str)
                #print(value)
                #['SORT', '(PARAM)']
                pass
                values.append(value)
                #asdf
            else:
                values.append(value)

        return cls(values)

def _set_options_from_line(line: str, value: str, options: list[str]) -> list[str]:
    """ECHO = PUNCH,SORT(MAT1,PARAM)"""
    is_comma = ',' in value
    is_paren = '(' in value
    if is_comma and is_paren:
        icomma, iparen = _get_icomma_iparen(value)
        #if icomma == -1 and iparen == -1:
            #raise NotImplementedError(f'ECHO is_comma/is_paren; line={line!r}')

        if iparen < icomma:
            # ECHO = SORT(PARAM,EIGC,EIGRL,FREQ,DESVAR,DCONSTR,DRESP1,DRESP2,DEQATN,DVPREL1)
            # ECHO = SORT(EXCEPT DMI,DMIG)
            # ECHO = UNSORT / SORT (EXCEPT DMI, DMIG), PUNCH(BSTBULK/NEWBULK)
            # ECHO = SORT(EXCEPT DMI,DMIG),PUNCH(BSTBULK)
            base, paren_open = value.split('(', 1) # [:iparen], value[iparen+1:]
            options.append(base)
            paren, end = paren_open.split(')', 1)
            sparen = paren.split(',')
            options.append(sparen)
            end = end.strip(' ,')
            if end == '':
                return

            _set_options_from_line(line, end, options)
        else:
            # ECHO = SORT,PUNCH(BSTBULK)
            #'PUNCH,SORT(MAT1,PARAM)'
            base, paren = value.split(',', 1)
            options.append(base)
            _set_options_from_line(line, paren, options)
    elif is_comma:
        sline = value.split(',', 1)
        options += sline
    elif is_paren:
        # ECHO = PUNCH(BSTBULK)
        # ECHO = PUNCH(NEWBULK)
        value2, options_ = value.split('(', 1)
        value2 = value2.strip()
        options_ = options_.strip()
        assert options_[-1] == ')', options_
        inner = options_[:-1]
        options += [value2, [inner]]
    else:
        # single value
        options.append(value)
    return options

def _get_icomma_iparen(value: str) -> tuple[int, int]:
    #nvalue = len(value)
    icomma = value.find(',')
    iparen = value.find('(')
    #if icomma == nvalue:
        #icomma = -1
    #if iparen == nvalue:
        #iparen = -1
    return icomma, iparen

def _strip_echo_options(options: list[str]) -> list[str]:
    options2 = []
    for option in options:
        if isinstance(option, str):
            options2.append(option)
        elif isinstance(option, list):
            option2 = [optioni.strip() for optioni in option]
            options2.append(option2)
        else:
            raise NotImplementedError(option)
    return options2
#-------------------------------------------------------------------------------

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
    allowed_keys: set[str] = set([])

    # key:(type, allowed_values)
    allowed_values: dict[str, Union[float, str]] = {}

    # the allowed value for the key, options, value approach
    allowed_strings: set[str] = set([])

    # maps something like INIT to INITIAL
    duplicate_names: dict[Any, Any] = {}

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
        value : list[str]
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

def split_by_mixed_commas_parentheses(str_options: str) -> list[str]:
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
    options : list[str]
        something that's actually parseable
        ['PRINT', 'SET=(G,N,N+AUTOSPC,F,A)', 'DATAREC=NO']
        ['PRINT', 'SET=(G,N,F,A)',           'CGI=NO',    'WEIGHT']

    """
    options_start = []
    options_end = []
    options_start_new = []  # type: list[str]
    options_end_new = []  # type: list[str]

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

    @classmethod
    def load_hdf5(cls, h5_file, encoding: str) -> Any:
        from pyNastran.utils.dict_to_h5py import _cast
        #'data', 'key', 'options', 'param_type', 'value'
        #print(h5_file, h5_file.keys())
        values = _cast(h5_file['value'])
        options_bytes = _cast(h5_file['options'])
        #data = _cast(h5_file['data'])
        key = _cast(h5_file['key'])
        #b'GROUNDCHECK' b'YES' [b'SET=(G,F)', b'DATAREC=NO', b'RTHRESH=.8']
        #'options'
        value = values.decode(encoding)
        options = [option.decode(encoding) for option in options_bytes]
        print(key, values, options, ) # data
        #'data'
        return GROUNDCHECK('GROUNDCHECK', value, options)


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
    }  # type: dict[str, Union[str, int]]
    #alternate_names = {'PRES'}
    #allow_ints = True

    def __init__(self):
        super(CheckCard, self).__init__()

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
                    'MATOP4', 'MATRIXOP4',
                    'MODACC'}

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
                keys_none_bytes = encode_str_value_list(keys_none, encoding)
                data_group.create_dataset('keys_none', data=keys_none_bytes)
            if keys:
                keys_bytes = encode_str_list(keys, encoding)
                values_bytes = encode_str_value_list(values, encoding)
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
                subgroupi_keys = list(subgroupi)
                if len(subgroupi_keys) > 0:
                    print(subgroupi, subgroupi_keys)
                    data_keys_group = subgroupi['keys']
                    data_keys = _cast(data_keys_group)
                    data_keys = decode_bytes_list(data_keys, encoding)

                    keys_none = subgroupi['keys_none']
                    keys_none = decode_bytes_list(keys_none, encoding)

                    data_values = _cast(subgroupi['values'])
                    data_values = decode_bytes_list(data_values, encoding)

                    keys = data_keys + keys_none
                    values = data_values + [None] * len(keys_none)
                    data = [(key, value) for (key, value) in zip(keys, values)]
                else:
                    data = []
            elif key == 'param_type':
                pass
            else:
                raise NotImplementedError(key)
        return EXTSEOUT(data)

    @classmethod
    def add_from_case_control(cls, line: str):
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

MATRIX_MAP = {
    #A2GG=ADMIG
    #A2GG=ADMIG1, ADMIG2, ADMIG3
    #A2GG=1.25*ADMIG1, 1.0*ADMIG2, 0.75*ADMIG3
    #SET 100=A1, A2
    #A2GG=100
    'A2GG': A2GG,

    #B2GG=BDMIG
    #B2GG=BDMIG1, BDMIG2, BDMIG3
    #B2GG=1.25*BDMIG1, 1.0*BDMIG2, 0.75*BDMIG3
    #SET 100=B1, B2
    #B2GG=100
    'B2GG': B2GG,

    #M2GG=MDMIG
    #M2GG=MDMIG1, MDMIG2, MDMIG3
    #M2GG=1.25*MDMIG1, 1.0*MDMIG2, 0.75*MDMIG3
    #SET 100=M1, M2
    #M2GG=100
    'M2GG': M2GG,

    #K2GG=KDMIG
    #K2GG=KDMIG1, KDMIG2, KDMIG3
    #K2GG=1.25*KDMIG1, 1.0*KDMIG2, 0.75*KDMIG3
    #SET 100=K1, K2
    #K2GG=100
    'K2GG': K2GG,

    #P2G=LDMIG
    #P2G=L1, L2, L3
    #SET 100=LDMIG, L1, L8
    #P2G=100
    #P2G=1.25*L1, 1.0*L2, 0.75*L3
    'P2G': P2G,
    # -----------------------------
    #B2PP=BDMIG
    #B2PP=BDMIG1, BDMIG2, BDMIG3
    #B2PP=1.25*BDMIG1, 1.0*BDMIG2, 0.75*BDMIG3
    #B2PP=(1.25,0.5)*BDMIG1, (1.0,0.0)*BDMIG2, (0.75,-2.2)*BDMIG3
    #SET 100=B1, B2
    #B2PP=100
    'B2PP': B2PP,

    #M2PP=MDMIG
    #M2PP=MDMIG1, MDMIG2, MDMIG3
    #M2PP=1.25*MDMIG1, 1.0*MDMIG2, 0.75*MDMIG3
    #M2PP=(1.25,0.5)*MDMIG1, (1.0,0.0)*MDMIG2, (0.75,-2.2)*MDMIG3
    #SET 100=M1, M2
    #M2PP=100
    'M2PP': M2PP,

    #K2PP=KDMIG
    #K2PP=KDMIG1, KDMIG2, KDMIG3
    #K2PP=1.25*KDMIG1, 1.0*KDMIG2, 0.75*KDMIG3
    #K2PP=(1.25,0.5)*KDMIG1, (1.0,0.0)*KDMIG2, (0.75,-2.2)*KDMIG3
    #SET 100=K1
    'K2PP': K2PP,

    # --------------------------------
    #K42GG=KDMIG
    #K42GG=KDMIG1, KDMIG2, KDMIG3
    #K42GG=1.25*KDMIG1, 1.0*KDMIG2, 0.75*KDMIG3
    #SET 100=K1, K2
    #K42GG=100
    'K42GG': K42GG,
}

OBJ_MAP = {
    'EXTSEOUT': EXTSEOUT,
    'GROUNDCHECK': GROUNDCHECK,
    'MEFFMASS': MEFFMASS,
    'VOLUME': VOLUME,
    'SURFACE': SURFACE,

    'AEROF': AEROF,
    'APRES': APRES,
    'GPKE': GPKE,
    'GPRSORT': GPRSORT,
    'GPSDCON': GPSDCON,
    'GPSTRESS': GPSTRESS,
    'HARMONICS': HARMONICS,
    'MODCON': MODCON,
    'OFREQUENCY': OFREQUENCY,
    'OMODES': OMODES,
    'SEALL': SEALL,
    'SEDR': SEDR,
    'SET': SET,
    'SETMC': SETMC,
    'SUPER': SUPER,
    'ECHO': ECHO,
}
OBJ_MAP.update(MATRIX_MAP)
