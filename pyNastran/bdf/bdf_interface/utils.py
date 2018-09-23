"""
Defines various utilities for BDF parsing including:
 - to_fields
"""
from __future__ import print_function, unicode_literals
import os
import sys
import inspect
import warnings
from six import StringIO, string_types
from typing import List, Union, Dict, Tuple, Optional

import pyNastran
from pyNastran.bdf.errors import CardParseSyntaxError
from pyNastran.bdf.cards.utils import wipe_empty_fields


_REMOVED_LINES = [
    '$EXECUTIVE CONTROL DECK',
    '$CASE CONTROL DECK',
    '$NODES', '$SPOINTS', '$ELEMENTS',
    '$PARAMS', '$PROPERTIES', '$ELEMENTS_WITH_PROPERTIES',
    '$ELEMENTS_WITH_NO_PROPERTIES (PID=0 and unanalyzed properties)',
    '$UNASSOCIATED_PROPERTIES',
    '$MATERIALS', '$THERMAL MATERIALS',
    '$CONSTRAINTS', '$SPCs', '$MPCs', '$RIGID ELEMENTS',
    '$LOADS', '$AERO', '$AERO CONTROL SURFACES',
    '$STATIC AERO', '$FLUTTER', '$GUST',
    '$DYNAMIC', '$OPTIMIZATION',
    '$COORDS', '$THERMAL', '$TABLES', '$RANDOM TABLES',
    '$SETS', '$CONTACT', '$REJECTS', '$REJECT_LINES',
    '$PROPERTIES_MASS', '$MASSES',
]
EXPECTED_HEADER_KEYS_CHECK = ['version', 'encoding', 'punch', 'nnodes', 'nelements', 'dumplines']
EXPECTED_HEADER_KEYS_NO_CHECK = ['skip_cards', 'units']


def _to_fields_mntpnt1(card_lines):
    # type: (List[str]) -> List[str]
    assert len(card_lines) == 2, card_lines
    line1, line2 = card_lines

    base = line1[:24]
    unused_comment = line1[24:]  # len=56 max
    assert ',' not in base, base
    assert '\t' not in base, base

    assert ',' not in line2, card_lines
    assert '\t' not in line2, card_lines
    assert '*' not in line2, card_lines
    fields = [
        line1[0:8],
        line1[8:16], line1[16:24], line1[24:32], line1[32:40], line1[40:48],
        line1[48:56], line1[56:64], line1[64:72],

        line2[8:16], line2[16:24], line2[24:32], line2[32:40], line2[40:48],
        line2[48:56], line2[56:64], line2[64:72],
    ]
    return fields

#def to_fields_long(card_lines, card_name):
    ## type: (List[str], str) -> List[str]
    #"""
    #Converts a series of lines in a card into string versions of the field.
    #Handles large, small, and CSV formatted cards.

    #Doesn't consider Nastran's rule about 72 character width fields,
    #which is nice when you have poorly formatted BDFs.

    #Parameters
    #----------
    #lines : List[str]
        #the lines of the BDF card object
    #card_name : str
        #the card_name -> 'GRID'

    #Returns
    #-------
    #fields : List[str]
        #the string formatted fields of the card

    #.. warning:: this function is used by the reader and isn't intended
                 #to be called by a separate process

    #.. code-block:: python

      #>>> card_lines = ['GRID,1,,1.0,2.0,3.0']
      #>>> card_name = 'GRID'
      #>>> fields = to_fields_long(lines, card_name)
      #>>> fields
      #['GRID', '1', '', '1.0', '2.0', '3.0']
    #"""
    #fields = []  # type: List[str]

    #if card_name == 'MONPNT1':
        #return _to_fields_mntpnt1(card_lines)

    ## first line
    #line = card_lines.pop(0)
    #if '=' in line:
        #msg = 'card_name=%r\nequal signs are not supported...line=%r' % (card_name, line)
        #raise CardParseSyntaxError(msg)

    #if '\t' in line:
        #line = line.expandtabs()
        #if ',' in line:
            #msg = 'tabs and commas in the same line are not supported...line=%r' % line
            #raise CardParseSyntaxError(msg)

    #if '*' in line:  # large field
        #if ',' in line:  # csv
            #new_fields = line.split(',')[:5]
            #for i in range(5 - len(new_fields)):
                #new_fields.append('')
        #else:  # standard
            #new_fields = [line[0:8], line[8:24], line[24:40], line[40:56],
                          #line[56:72]]
        #fields += new_fields
        #assert len(fields) == 5, fields
    #else:  # small field
        #if ',' in line:  # csv
            #new_fields = line.split(',')[:9]
            #for i in range(9 - len(new_fields)):
                #new_fields.append('')
        #else:  # standard
            #new_fields = [line[0:8], line[8:16], line[16:24], line[24:32],
                          #line[32:40], line[40:48], line[48:56], line[56:64],
                          #line[64:72]]
        #fields += new_fields
        #assert len(fields) == 9, fields

    #for j, line in enumerate(card_lines): # continuation lines
        #if '=' in line and card_name != 'EIGRL':
            #msg = 'card_name=%r\nequal signs are not supported...line=%r' % (card_name, line)
            #raise CardParseSyntaxError(msg)
        #if '\t' in line:
            #line = line.expandtabs()
            #if ',' in line:
                #msg = 'tabs and commas in the same line are not supported...line=%r' % line
                #raise CardParseSyntaxError(msg)

        #if '*' in line:  # large field
            #if ',' in line:  # csv
                #new_fields = line.split(',')[1:5]
                #for i in range(4 - len(new_fields)):
                    #new_fields.append('')
            #else:  # standard
                #new_fields = [line[8:24], line[24:40], line[40:56], line[56:72]]
            #assert len(new_fields) == 4, new_fields
        #else:  # small field
            #if ',' in line:  # csv
                #new_fields = line.split(',')[1:9]
                #for i in range(8 - len(new_fields)):
                    #new_fields.append('')
            #else:  # standard
                #new_fields = [line[8:16], line[16:24], line[24:32],
                              #line[32:40], line[40:48], line[48:56],
                              #line[56:64], line[64:72]]
            #if len(new_fields) != 8:
                #nfields = len(new_fields)
                #msg = 'nFields=%s new_fields=%s' % (nfields, new_fields)
                #raise RuntimeError(msg)

        #fields += new_fields
    #return fields #[field.strip() for field in fields]

def to_fields(card_lines, card_name):
    # type: (List[str], str) -> List[str]
    """
    Converts a series of lines in a card into string versions of the field.
    Handles large, small, and CSV formatted cards.

    Parameters
    ----------
    lines : List[str]
        the lines of the BDF card object
    card_name : str
        the card_name -> 'GRID'

    Returns
    -------
    fields : List[str]
        the string formatted fields of the card

    .. warning:: this function is used by the reader and isn't intended
                 to be called by a separate process

    .. code-block:: python

      >>> card_lines = ['GRID,1,,1.0,2.0,3.0']
      >>> card_name = 'GRID'
      >>> fields = to_fields(lines, card_name)
      >>> fields
      ['GRID', '1', '', '1.0', '2.0', '3.0']

    """
    fields = []  # type: List[str]

    if card_name == 'MONPNT1':
        return _to_fields_mntpnt1(card_lines)

    # first line
    line = card_lines[0]
    if '=' in line:
        msg = 'card_name=%r\nequal signs are not supported...line=%r' % (card_name, line)
        raise CardParseSyntaxError(msg)

    if '\t' in line:
        line = expand_tabs(line)

    if '*' in line:  # large field
        if ',' in line:  # csv
            new_fields = line.split(',')[:5]
            for unused_i in range(5 - len(new_fields)):
                new_fields.append('')
        else:  # standard
            new_fields = [line[0:8], line[8:24], line[24:40], line[40:56],
                          line[56:72]]
        fields += new_fields
        assert len(fields) == 5, fields
    else:  # small field
        if ',' in line:  # csv
            new_fields = line.split(',')[:9]
            for unused_i in range(9 - len(new_fields)):
                new_fields.append('')
        else:  # standard
            new_fields = [line[0:8], line[8:16], line[16:24], line[24:32],
                          line[32:40], line[40:48], line[48:56], line[56:64],
                          line[64:72]]
        fields += new_fields
        assert len(fields) == 9, fields

    for line in card_lines[1:]: # continuation lines
        if '=' in line and card_name != 'EIGRL':
            msg = 'card_name=%r\nequal signs are not supported...\nline=%r' % (card_name, line)
            raise CardParseSyntaxError(msg)
        if '\t' in line:
            line = expand_tabs(line)

        if '*' in line:  # large field
            if ',' in line:  # csv
                new_fields = line.split(',')[1:5]
                for unused_i in range(4 - len(new_fields)):
                    new_fields.append('')
            else:  # standard
                new_fields = [line[8:24], line[24:40], line[40:56], line[56:72]]
            assert len(new_fields) == 4, new_fields
        else:  # small field
            if ',' in line:  # csv
                new_fields = line.split(',')[1:9]
                for unused_i in range(8 - len(new_fields)):
                    new_fields.append('')
            else:  # standard
                new_fields = [line[8:16], line[16:24], line[24:32],
                              line[32:40], line[40:48], line[48:56],
                              line[56:64], line[64:72]]
            if len(new_fields) != 8:
                nfields = len(new_fields)
                msg = 'nfields=%s new_fields=%s' % (nfields, new_fields)
                raise RuntimeError(msg)

        fields += new_fields
    return fields

def expand_tabs(line):
    """expands the tabs; breaks if you mix commas and tabs"""
    line = line.expandtabs()
    if ',' in line:
        line = line.replace('\t', '')
        msg = 'tabs and commas in the same line are not supported...\nline=%r' % line
        raise CardParseSyntaxError(msg)
    return line

def to_fields_replication(card_lines):
    """
    Converts a series of lines in a card into string versions of the field.
    Handles large, small, and CSV formatted cards.  Same as to_fields, but
    uses a different method to determine if it's large or small field format.

    Parameters
    ----------
    lines : List[str]
        the lines of the BDF card object

    Returns
    -------
    fields : List[str]
        the string formatted fields of the card

    .. warning:: this function is used by the reader and isn't intended
                 to be called by a separate process

    .. code-block:: python

      >>> card_lines = ['GRID,1,,1.0,2.0,3.0']
      >>> card_name = 'GRID'
      >>> fields = to_fields_replication(lines)
      >>> fields
      ['GRID', '1', '', '1.0', '2.0', '3.0']

    """
    #print('to_fields_replicationA =', card_lines)
    #line0 = card_lines[0]
    fields = []
    for iline, line in enumerate(card_lines):
        line = line.rstrip()

        if '\t' in line:
            line = expand_tabs(line)
        #print('  line = %r' % line)
        if ',' in line:
            assert '\t' not in line, '%r' % line
            sline = line.split(',')
            #print('sline = ', iline, sline)
            if iline == 0:
                card_name = sline[0]
                #assert card_name != '=', card_lines
                fields.append(card_name)

            if '*' in sline[0]:
                fields.extend(sline[1:5])
                for unused_i in range(5 - len(sline)):
                    fields.append('')
            else:
                fields.extend(sline[1:9])
                for unused_i in range(9 - len(sline)):
                    fields.append('')
        else:
            assert '\t' not in line, line
            if iline == 0:
                card_name = line[0:8]
                assert card_name != '=', card_lines
                fields.append(card_name)

            if '*' in card_name:
                fields.extend([line[8:24], line[24:40], line[40:56],
                               line[56:72]])
            else:
                fields.extend([line[8:16], line[16:24], line[24:32],
                               line[32:40], line[40:48], line[48:56], line[56:64],
                               line[64:72]])
    assert '*' not in card_name, card_lines
    wiped_fields = wipe_empty_fields(fields)
    for field in fields:
        sfield = field.strip()
        #while '= ' in sfield:
            #sfield = field.replace('= ','=')
            #print('sfield=%r' % sfield)
        if ' ' in sfield:
            raise RuntimeError('field=%r has embedded blanks\nfields=%s' % (sfield, fields))
    return wiped_fields

def get_nrepeats(field, old_card, new_card):
    """=4, =(11)"""
    msg = 'field=%r; expected =(1), =2, ...\nold_card=%s\nnew_card=%s' % (
        field, old_card, new_card)
    assert field[0] == '=', msg
    assert '.' not in field, msg
    assert '*' not in field, msg
    if '(' in field or ')' in field:
        assert field[1] == '(', msg
        assert field[-1] == ')', msg
        fieldi = field[2:-1]
        assert '(' not in fieldi, msg
        assert ')' not in fieldi, msg
    else:
        assert '(' not in field, msg
        assert ')' not in field, msg
        fieldi = field[1:]
    nrepeats = int(fieldi)
    return nrepeats

def float_replication(field, old_field):
    """*4., *(11.5)"""
    msg = 'field=%r; expected *(1.), *2., ..., *11.' % field
    assert field[0] == '*', msg
    assert '.' in field, msg
    assert len(field) >= 3, msg
    if '(' in field:
        assert field[1] == '(', msg
        assert field[-1] == ')', msg
        fieldi = field[2:-1]
        assert '(' not in fieldi, msg
        assert ')' not in fieldi, msg
    else:
        assert '(' not in field, msg
        assert ')' not in field, msg
        fieldi = field[1:]
    nfloat = float(fieldi)
    field2 = nfloat + float(old_field)
    return field2

def int_replication(field, old_field):
    """*4, *(11)"""
    msg = 'field=%r; expected *(1), *2, ..., *11' % field
    assert field[0] == '*', msg
    assert '.' not in field, msg
    if '(' in field:
        assert field[1] == '(', msg
        assert field[-1] == ')', msg
        fieldi = field[2:-1]
        assert '(' not in fieldi, msg
        assert ')' not in fieldi, msg
    else:
        assert '(' not in field, msg
        assert ')' not in field, msg
        fieldi = field[1:]

    assert '*' not in fieldi, msg
    assert len(field) >= 2, msg
    #assert len(field) == 2, 'field=%r; expected *1, *2, *3, ..., *9' % field
    # integer
    nint = int(fieldi)
    field2 = nint + int(old_field)
    return field2

def parse_executive_control_deck(executive_control_lines):
    """Extracts the solution from the executive control deck"""
    sol = None
    method = None
    sol_iline = None
    for (i, eline) in enumerate(executive_control_lines):
        uline = eline.strip().upper()  # uppercase line
        uline = uline.split('$')[0].expandtabs()
        if uline[:4] in ['SOL ']:
            if ',' in uline:
                sline = uline.split(',')  # SOL 600,method
                sol_value = sline[0].strip()
                method = sline[1].strip()
            else:
                sol_value = uline
                method = None

            if sol is None:
                sol = sol_value[3:].strip()
            else:
                raise ValueError('cannot overwrite solution existing='
                                 '|SOL %s| new =|%s|' % (sol, uline))
            sol_iline = i
    return sol, method, sol_iline


def _parse_pynastran_header(line):
    # type: (str) -> Tuple[Optional[str], Optional[str]]
    """
    Parameters
    ----------
    line : str
        the line to parse (e.g., '$ pyNastran: version=NX')

    Returns
    -------
    key : str / None
        the key for the parameters
        str : valid (e.g., 'version')
        None : invalid
    value : str / None
        the key for the parameters
        str : valid (e.g., 'NX')
        None : invalid

    Search for data of the form:
        ..code-block :: python
            $ pyNastran: version=NX
            $ pyNastran: encoding=latin-1
            $ pyNastran: punch=True
            $ pyNastran: dumplines=True
            $ pyNastran: nnodes=10
            $ pyNastran: nelements=100
            $ pyNastran: skip_cards=PBEAM,CBEAM
            $ pyNastran: units=in,lb,s
            $ pyNastran: skip elements=12345,6,7,8

    If we find:
        ..code-block :: python
            $$ pyNastran: version=NX

    or a line without a valid pyNastran flag, we'll stop reading,
    even a valid header statement is on the following line.

    """
    lline = line[1:].lower().strip()
    if len(lline) == 0 or lline[0] == '$':
        key = None
        value = None
    elif 'pynastran' in lline:
        base, word = lline.split(':', 1)
        if base.strip() != 'pynastran':
            msg = 'unrecognized pyNastran marker\n'
            msg += 'line=%r' % line
            raise SyntaxError(msg)
        key, value = word.strip().split('=')
        key = key.strip()
        value = value.strip()
        if key in EXPECTED_HEADER_KEYS_CHECK:
            assert ' ' not in value, 'value=%r' % value
        elif key in EXPECTED_HEADER_KEYS_NO_CHECK:
            pass
        elif 'skip ' in key:
            pass
        else:
            msg = '\nunrecognized pyNastran key=%r type(key)=%s\n' % (key, type(key))
            msg += 'line=%r\n' % line
            msg += 'expected_keys = [%s]\n' % ', '.join(EXPECTED_HEADER_KEYS_CHECK + EXPECTED_HEADER_KEYS_NO_CHECK)
            msg += 'type(key0) = %s' % type(EXPECTED_HEADER_KEYS_CHECK[0])
            print(msg)
            raise KeyError(msg)
    else:
        key = None
        value = None
    return key, value


#def clean_empty_lines(lines):
    ## type: (List[str]) -> List[str]
    #"""
    #Removes leading and trailing empty lines
    #don't remove internally blank lines
    #"""
    #found_lines = False
    #if len(lines) < 2:
        #return lines

    #for i, line in enumerate(lines):
        #if not found_lines and line:
            #found_lines = True
            #n1 = i
            #n2 = i + 1
        #elif found_lines and line:
            #n2 = i + 1
    #lines2 = lines[n1:n2]
    #return lines2


def print_filename(filename, relpath=True):
    # type: (str, str) -> str
    """
    Takes a path such as C:/work/fem.bdf and locates the file using
    relative paths.  If it's on another drive, the path is not modified.

    Parameters
    ----------
    filename : str
        a filename string

    Returns
    -------
    filename_string : str
        a shortened representation of the filename

    """
    if isinstance(filename, StringIO):
        return '<StringIO>'
    drive_letter = os.path.splitdrive(os.path.abspath(filename))[0]
    if drive_letter == os.path.splitdrive(os.curdir)[0] and relpath:
        return os.path.relpath(filename)
    return filename

def deprecated(old_name, new_name, deprecated_version, levels=None):
    # type: (str, str, str, Optional[List[int]]) -> None
    """
    Throws a deprecation message and crashes if past a specific version.

    Parameters
    ----------
    old_name : str
        the old function name
    new_name : str
        the new function name
    deprecated_version : float
        the version the method was first deprecated in
    levels : List[int]
        the deprecation levels to show
        [1, 2, 3] shows 3 levels up from this function

    TODO: turn this into a decorator?

    """
    assert isinstance(deprecated_version, string_types), type(deprecated_version)
    assert isinstance(levels, list), type(levels)
    assert old_name != new_name, "'%s' and '%s' are the same..." % (old_name, new_name)

    version = pyNastran.__version__.split('_')[0]
    dep_ver_tuple = tuple([int(i) for i in deprecated_version.split('.')])
    ver_tuple = tuple([int(i) for i in version.split('.')[:2]])

    #new_line = ''
    msg = "'%s' was deprecated in v%s (current=%s)" % (
        old_name, deprecated_version, version)
    if new_name:
        msg += "; replace it with '%s'\n" % new_name

    for level in levels:
        # jump to get out of the inspection code
        frame = sys._getframe(3 + level)
        line_no = frame.f_lineno
        code = frame.f_code
        try:
            #filename = os.path.basename(frame.f_globals['__file__'])
            filename = os.path.basename(inspect.getfile(code))
        except:
            print(code)
            raise

        source_lines, line_no0 = inspect.getsourcelines(code)
        delta_nlines = line_no - line_no0
        try:
            line = source_lines[delta_nlines]
        except:
            break
        msg += '  %-25s:%-4s %s\n' % (filename, str(line_no) + ';', line.strip())

    if ver_tuple > dep_ver_tuple: # or 'id' in msg:
        # fail
        raise NotImplementedError(msg)
    else:
        warnings.warn(msg, DeprecationWarning)
