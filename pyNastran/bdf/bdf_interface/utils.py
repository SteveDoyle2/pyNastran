"""
Defines various utilities for BDF parsing including:
 - to_fields

"""
from __future__ import annotations
import os
from io import StringIO
from collections import defaultdict
from typing import List, Dict, Tuple, Optional, Any, TYPE_CHECKING

#import pyNastran
from pyNastran.bdf.errors import CardParseSyntaxError
if TYPE_CHECKING:  # pragma: no cover
    from pyNastran.bdf.bdf import BDF

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
EXPECTED_HEADER_KEYS_CHECK = [
    'version', 'encoding', 'nnodes', 'nelements',
    'punch', 'dumplines', 'is_superelements', # booleans
]
EXPECTED_HEADER_KEYS_NO_CHECK = ['skip_cards', 'units', 'code-block']


def _to_fields_mntpnt1(card_lines: List[str]) -> List[str]:
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


def to_fields(card_lines: List[str], card_name: str) -> List[str]:
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

def expand_tabs(line: str) -> str:
    """expands the tabs; breaks if you mix commas and tabs"""
    line = line.expandtabs()
    if ',' in line:
        line = line.replace('\t', '')
        msg = 'tabs and commas in the same line are not supported...\nline=%r' % line
        raise CardParseSyntaxError(msg)
    return line

def parse_executive_control_deck(
        executive_control_lines: List[str]) -> Tuple[Optional[int], Optional[str], Optional[int]]:
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


def _parse_pynastran_header(line: str) -> Tuple[Optional[str], Optional[str]]:
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
        try:
            key, value = word.strip().split('=', 1)
        except ValueError:
            msg = (
                'expected header of the form:\n'
                '$ pyNastran: version=NX\n'
                '$ pyNastran: encoding=latin-1\n'
                '$ pyNastran: punch=True\n'
                '$ pyNastran: dumplines=True\n'
                '$ pyNastran: nnodes=10\n'
                '$ pyNastran: nelements=100\n'
                '$ pyNastran: skip_cards=PBEAM,CBEAM\n'
                '$ pyNastran: units=in,lb,s\n'
                '$ pyNastran: skip elements=12345,6,7,8\n'
            )
            raise SyntaxError(msg)
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
            msg += 'expected_keys = [%s]\n' % ', '.join(
                EXPECTED_HEADER_KEYS_CHECK + EXPECTED_HEADER_KEYS_NO_CHECK)
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


def print_filename(filename: str, relpath: bool=True) -> str:
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


def _parse_dynamic_syntax(key: str,
                          dict_of_vars: Dict[str, Any],
                          log: Any) -> Dict[str, Any]:
    """
    Applies the dynamic syntax for %varName

    Parameters
    ----------
    key : str
        the uppercased key

    Returns
    -------
    value : int/float/str
        the dynamic value defined by dict_of_vars

    .. seealso:: :func: `set_dynamic_syntax`

    """

    key = key.strip()[1:]
    log.debug("dynamic key = %r" % key)
    #dict_of_vars = {'P5':0.5,'ONEK':1000.}
    if key not in dict_of_vars:
        msg = "key=%r not found in keys=%s" % (key, dict_of_vars.keys())
        raise KeyError(msg)
    return dict_of_vars[key]

def _get_card_name(lines: List[str], active_filename: str) -> Optional[str]:
    """
    Returns the name of the card defined by the provided lines

    Parameters
    ----------
    lines : list[str]
        the lines of the card

    Returns
    -------
    card_name : str
        the name of the card

    """
    card_name = lines[0][:8].rstrip('\t, ').split(',')[0].split('\t')[0].strip('*\t ')
    if len(card_name) == 0:
        return None
    if ' ' in card_name or len(card_name) == 0:
        msg = 'card_name=%r\nline=%r in filename=%r is invalid' \
              % (card_name, lines[0], active_filename)
        print(msg)
        raise CardParseSyntaxError(msg)
    return card_name.upper()

def fill_dmigs(model: BDF) -> None:
    """fills the DMIx cards with the column data that's been stored"""
    for name, card_comments in model._dmig_temp.items():
        card0, unused_comment0 = card_comments[0]
        card_name = card0[0]
        card_name = card_name.rstrip(' *').upper()

        if card_name == 'DMIG':
            # if field2 == 'UACCEL':  # special DMIG card
            card = model.dmig[name]
        elif card_name == 'DMI':
            card = model.dmi[name]
        elif card_name == 'DMIJ':
            card = model.dmij[name]
        elif card_name == 'DMIJI':
            card = model.dmiji[name]
        elif card_name == 'DMIK':
            card = model.dmik[name]
        elif card_name == 'DMIAX':
            card = model.dmiax[name]
        else:  # pragma: no cover
            raise NotImplementedError(card_name)

        for (card_obj, comment) in card_comments:
            card._add_column(card_obj, comment=comment)
        card.finalize()

    # empty the _dmig_temp variable
    model._dmig_temp = defaultdict(list)
