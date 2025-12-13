"""
Defines various utilities for BDF parsing including:
 - to_fields

"""
from __future__ import annotations
import os
from io import StringIO
import warnings
from collections import defaultdict
from typing import Optional, Any, TYPE_CHECKING

from pyNastran.bdf.errors import CardParseSyntaxError
if TYPE_CHECKING:  # pragma: no cover
    from pyNastran.utils import PathLike
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
    'punch', 'dumplines', 'is_superelements',  # booleans
]
EXPECTED_HEADER_KEYS_NO_CHECK = ['skip_cards', 'units', 'code-block']


def _to_fields_mntpnt1(card_lines: list[str]) -> list[str]:
    """splits a MONPNT1"""
    assert len(card_lines) == 2, card_lines
    line1, line2 = card_lines

    fields = _monpnt_line1_fields(line1)

    #assert ',' not in line2, card_lines
    #assert '\t' not in line2, card_lines
    #assert '*' not in line2, card_lines
    fields += _monpnt_line2_to_fields(line2)
    return fields


def _to_fields_mntpnt3(card_lines: list[str]) -> list[str]:
    assert len(card_lines) in {2, 3}, card_lines
    line1 = card_lines[0]

    fields = _monpnt_line1_fields(line1)

    #assert ',' not in line2, card_lines
    #assert '\t' not in line2, card_lines
    #assert '*' not in line2, card_lines
    card_lines_end = card_lines[1:]
    for line in card_lines_end:
        fields += _monpnt_line2_to_fields(line)
    return fields

#def _to_fields_monsumt(card_lines: list[str]) -> list[str]:
    #assert len(card_lines) == 3, card_lines
    #line1, line2, line3 = card_lines

    #fields = _monpnt_line1_fields(line1)

    ##assert ',' not in line2, card_lines
    ##assert '\t' not in line2, card_lines
    ##assert '*' not in line2, card_lines
    #fields += _monpnt_line2_to_fields(line2)
    #fields += _monpnt_line2_to_fields(line3)
    #return fields


def _monpnt_line1_fields(line1: str) -> list[str]:
    """splits the first line of a MONPNT1/MONPNT3"""

    if '\t' in line1:
        line1 = line1.expandtabs()
        assert ',' not in line1[:16], line1

    label = line1[:24]
    unused_comment = line1[24:]  # len=56 max
    #assert ',' not in label, f'base={label!r}'
    #assert '\t' not in label, f'base={label!r}'

    fields = [
        line1[0:8],
        line1[8:16], line1[16:24], line1[24:32], line1[32:40], line1[40:48],
        line1[48:56], line1[56:64], line1[64:72],
    ]
    return fields


def _monpnt_line2_to_fields(line2: str) -> list[str]:
    """splits a MONPNT3"""
    if '\t' in line2:
        line2 = line2.expandtabs()
        assert ',' not in line2[:16], line2

    if ',' in line2:
        # drop off the first field of row2
        fields = line2.split(',')[1:]
    else:
        fields = [
            line2[8:16], line2[16:24], line2[24:32], line2[32:40], line2[40:48],
            line2[48:56], line2[56:64], line2[64:72],
        ]
    return fields


def _to_fields_micpnt(card_lines: list[str]) -> list[str]:
    """splits a MICPNT"""
    line1 = card_lines[0]
    fields = _expand_2_values_name(line1)
    return fields


def _to_fields_trimfnc(card_lines: list[str]) -> list[str]:
    """splits a TRIMFNC"""
    assert len(card_lines) == 1, card_lines
    line = card_lines[0]
    if '\t' in line:
        line = expand_tabs(line)

    new_fields = [line[0:8], line[8:16], line[16:24], line[24:32],
                  line[32:40], line[40:48], line[48:56], line[56:72]]
    return new_fields

def _to_fields_extfile(card_lines: list[str]) -> list[str]:
    """splits a EXTFILE"""
    assert len(card_lines) == 1, card_lines
    line = card_lines[0]
    assert '\t' not in line, line
    # if '\t' in line:
    #     line = expand_tabs(line)
    new_fields = [line[0:8], line[8:16], line[16:72]]
    return new_fields

def _to_fields_dollar(card_lines: list[str], card_name: str) -> list[str]:
    return _to_fields_standard(card_lines, card_name)

def _to_fields_amlreg(card_lines: list[str]) -> list[str]:
    """splits an AMLREG"""
    line1 = card_lines[0]
    line2 = card_lines[1]
    fields = _expand_2_values_name(line1)
    fields += _monpnt_line2_to_fields(line2)
    return fields


def _to_fields_group(card_lines: list[str]) -> list[str]:
    """
    splits an GROUP


    +-------+------+------------------------------------+
    | GROUP | 10   | Assembly AA4                       |
    +-------+------+------------------------------------+
    |       | META | 100 RPM                            |
    +-------+------+------------------------------------+
    |       | META | Optionally continue the meta data  |
    +-------+------+-----+------+-----+-----+---+---+---+
    |       | GRID |  1  |  2   |  3  |  4  | 5 | 6 | 7 |
    +-------+------+-----+------+-----+-----+---+---+---+
    |       |      |  8  |      |     |     |   |   |   |
    +-------+------+-----+------+-----+-----+---+---+---+
    |       | GRID | 10  | THRU | 20  |     |   |   |   |
    +-------+------+-----+------+-----+-----+---+---+---+
    |       | GRID | 100 | THRU | 200 |     |   |   |   |
    +-------+------+-----+------+-----+-----+---+---+---+
    |       | GRID | 341 | THRU | 360 |  BY | 2 |   |   |
    +-------+------+-----+------+-----+-----+---+---+---+
    |       | ELEM | 30  | THRU | 40  |     |   |   |   |
    +-------+------+-----+------+-----+-----+---+---+---+
    |       | PROP | ALL |      |     |     |   |   |   |
    +-------+------+-----+------+-----+-----+---+---+---+

    # TODO: buggy
    GROUP   991001  Description                                             +
    +       GRID    291     293     301     302     303     304     305     +
    +               306     1301    68408   68409   69144   69145   69146
    """
    # for iline, line in enumerate(card_lines):
    #     print(f'GROUP {iline}: {line!r}')

    fields = _monpnt_line1_fields(card_lines[0])
    for j, line in enumerate(card_lines[1:]):
        line_strip = line.lstrip('+ ')
        line = line.rstrip()
        if line_strip.startswith('META'):
            fieldsi = [
                line[8:16], line[16:24], line[24:32], line[32:40], line[40:48],
                line[48:56], line[56:64], line[64:72],
            ]
        elif line_strip.startswith(('GRID', 'ELEM', 'PROP')):
            fieldsi = [
                line[8:16], line[16:24], line[24:32], line[32:40], line[40:48],
                line[48:56], line[56:64], line[64:72],
            ]
            # fieldsi = _expand_2_values_name(line)
            # #print(fieldsi)
        else:
            fieldsi = [
                line[8:16], line[16:24], line[24:32], line[32:40], line[40:48],
                line[48:56], line[56:64], line[64:72],
            ]
            # msg = f'line:\n{line!r}\nline_strip:\n{line_strip!r}'
            # raise NotImplementedError(msg)
        # print(f'GROUP fields[{j}]: fields={fieldsi}')
        fields.extend(fieldsi)
    # print(f'GROUP: fields={fields}')
    return fields


def _expand_2_values_name(line1: str) -> list[str]:
    """helper for AMLREG, MICPNT"""
    line1 = line1.rstrip()
    if '\t' in line1:
        line1 = line1.expandtabs()
        assert ',' not in line1[:16], line1
    while ',' in line1[:24]:
        line1 = line1.replace(',', '\t', 1)
        line1 = line1.expandtabs()
    label = line1[24:72]
    unused_comment = line1[24:]  # len=56 max
    #assert ',' not in label, f'base={label!r}'
    #assert '\t' not in label, f'base={label!r}'

    fields = [
        line1[0:8],
        line1[8:16], line1[16:24],

        label,
        #line1[24:32], line1[32:40], line1[40:48],
        #line1[48:56], line1[56:64], line1[64:72],
    ]
    return fields


def _to_fields_set1(card_lines: list[str], card_name: str) -> list[str]:
    fields = []
    throw_length_warning = False
    for iline, line in enumerate(card_lines):
        if '\t' in line:
            line = expand_tabs(line)
        # print(f'line = {line}')
        line = line.rstrip('\n\r,')
        # print(f'line2 = {line}')
        if '*' in line:  # large field
            if ',' in line:  # csv
                new_fields = line.split(',')  # [:5]
                new_fields = [field.strip() for field in new_fields if field.strip()]
                if iline > 0:
                    new_fields = [''] + new_fields
                assert len(new_fields) <= 5, new_fields
                # for unused_i in range(5 - len(new_fields)):
                #     new_fields.append('')
                assert len(new_fields) == 5, new_fields
            else:  # standard
                new_fields = [line[0:8], line[8:24], line[24:40], line[40:56],
                              line[56:72]]
                end = line[72:].rstrip('+ ')
                assert len(end) == 0, line
        else:  # small field
            length_max = 9  # if iline == 0 else 8
            if ',' in line:  # csv
                new_fields = line.split(',')  # [:9]
                new_fields = [field.strip() for field in new_fields if field.strip()]
                if iline > 0:
                    new_fields = [''] + new_fields
                if len(new_fields) >= length_max:
                    throw_length_warning = True
                for unused_i in range(9 - len(new_fields)):
                    new_fields.append('')
                # assert len(new_fields) == 9, f'{new_fields}; {len(new_fields)}'
            else:  # standard
                new_fields = [line[0:8], line[8:16], line[16:24], line[24:32],
                              line[32:40], line[40:48], line[48:56], line[56:64],
                              line[64:72]]
                end = line[72:].rstrip('+ ')
                assert len(end) == 0, line

        # if ',' in line:
        #     sline = line.split(',')
        # elif '*' in line:
        #     sline = line.split()
        # else:
        #     sline = line.split()
        #
        # sline2 = [field.strip() for field in sline if field.strip()]
        # for field in sline2:
        #     assert len(field) < 8, fields
        #
        # print(sline2)
        fields += new_fields
    warnings.warn(f'SET1 was too long; {fields}')
    return fields


def to_fields_line0(card_line: str, card_name: str) -> list[str]:
    """
    Converts the first line of a card into the string versions of the line.
    Handles large, small, and CSV formatted cards.
    """
    # first line
    line = card_line.rstrip()
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
            assert len(new_fields) == 5, new_fields
        else:  # standard
            new_fields = [line[0:8], line[8:24], line[24:40], line[40:56],
                          line[56:72]]
    else:  # small field
        if ',' in line:  # csv
            new_fields = line.split(',')[:9]
            for unused_i in range(9 - len(new_fields)):
                new_fields.append('')
            assert len(new_fields) == 9, new_fields
        else:  # standard
            new_fields = [line[0:8], line[8:16], line[16:24], line[24:32],
                          line[32:40], line[40:48], line[48:56], line[56:64],
                          line[64:72]]
    return new_fields


def to_fields(card_lines: list[str], card_name: str,
              allow_tabs: bool=True) -> list[str]:
    """
    Converts a series of lines in a card into string versions of the field.
    Handles large, small, and CSV formatted cards.

    Parameters
    ----------
    card_lines : list[str]
        the lines of the BDF card object
    card_name : str
        the card_name -> 'GRID'
    allow_tabs : bool; default=True
        tabs are ok

    Returns
    -------
    fields : list[str]
        the string formatted fields of the card

    .. warning:: this function is used by the reader and isn't intended
                 to be called by a separate process

    .. code-block:: python

       >>> card_lines = ['GRID,1,,1.0,2.0,3.0']
       >>> card_name = 'GRID'
       >>> fields = to_fields(card_lines, card_name)
       >>> fields
       ['GRID', '1', '', '1.0', '2.0', '3.0']

    """
    if not allow_tabs:
        joined_lines = '\n'.join(card_lines)
        if not allow_tabs and '\t' in joined_lines:
            raise RuntimeError(f'There are tabs in:\n{joined_lines}')

    if card_name in ['MONPNT1', 'MONDSP1']:
        return _to_fields_mntpnt1(card_lines)
    elif card_name in ['MONPNT3', 'MONSUMT']:
        return _to_fields_mntpnt3(card_lines)
    elif card_name == 'GROUP':
        return _to_fields_group(card_lines)
    elif card_name == 'AMLREG':
        return _to_fields_amlreg(card_lines)
    elif card_name == 'MICPNT':
        return _to_fields_micpnt(card_lines)
    elif card_name == 'TRIMFNC':
        return _to_fields_trimfnc(card_lines)
    elif card_name == 'EXTFILE':
        return _to_fields_extfile(card_lines)
    elif card_name == 'MKAEROZ':
        return _to_fields_dollar(card_lines, card_name)
    #elif card_name == 'SET1':
        #return _to_fields_set1(card_lines, card_name)
    return _to_fields_standard(card_lines, card_name)


def _to_fields_standard(card_lines: list[str], card_name: str) -> list[str]:
    fields: list[str] = []
    # first line
    line = card_lines[0].rstrip()
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
            assert len(new_fields) == 5, new_fields
        else:  # standard
            new_fields = [line[0:8], line[8:24], line[24:40], line[40:56],
                          line[56:72]]
        fields += new_fields
    else:  # small field
        if ',' in line:  # csv
            new_fields = line.split(',')[:9]
            for unused_i in range(9 - len(new_fields)):
                new_fields.append('')
            assert len(new_fields) == 9, new_fields
        else:  # standard
            new_fields = [line[0:8], line[8:16], line[16:24], line[24:32],
                          line[32:40], line[40:48], line[48:56], line[56:64],
                          line[64:72]]
        fields += new_fields

    for line in card_lines[1:]:  # continuation lines
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
                assert len(new_fields) == 4, new_fields
            else:  # standard
                new_fields = [line[8:24], line[24:40], line[40:56], line[56:72]]
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
        msg = f'tabs and commas in the same line are not supported...\nline={line!r}'
        raise CardParseSyntaxError(msg)
    return line


def parse_executive_control_deck(
        executive_control_lines: list[str]) -> tuple[Optional[int], Optional[str], Optional[int], str]:
    """Extracts the solution from the executive control deck"""
    sol = None
    method = None
    sol_iline = None
    app = ''
    for (i, eline) in enumerate(executive_control_lines):
        uline = eline.strip().upper()  # uppercase line
        uline = uline.split('$')[0].expandtabs()
        if uline[:4] == 'SOL ':
            if ',' in uline:
                sline = uline.split(',')  # SOL 600,method
                sol_value = sline[0].strip()
                method = sline[1].strip()
            else:
                sol_value = uline
                method = None

            if sol is None:
                sol = sol_value[3:].strip(' \t=')
                if ',' not in sol:
                    try:
                        # SOL 101
                        sol = int(sol)
                    except ValueError:
                        # SOL SESTATIC
                        pass
            else:
                raise ValueError('cannot overwrite solution existing='
                                 f'|SOL {sol}| new={uline!r}')
            sol_iline = i
        elif uline.startswith('APP '):
            #print('uline = %r' % uline)
            sline = uline.strip().split()
            assert len(sline) == 2, sline
            app = sline[1]
            assert app in {'HEAT', 'DISP', 'COUPLED', 'DISPLACEMENT', 'DMAP'}, f'uline={uline!r}'
    return sol, method, sol_iline, app


def _parse_pynastran_header(line: str) -> tuple[Optional[str], Optional[str]]:
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


# def clean_empty_lines(lines: list[str]) -> list[str]:
#     """
#     Removes leading and trailing empty lines
#     don't remove internally blank lines
#     """
#     found_lines = False
#     if len(lines) < 2:
#         return lines
#
#     for i, line in enumerate(lines):
#         if not found_lines and line:
#             found_lines = True
#             n1 = i
#             n2 = i + 1
#         elif found_lines and line:
#             n2 = i + 1
#     lines2 = lines[n1:n2]
#     return lines2


def print_filename(filename: PathLike | StringIO,
                   relpath: bool=True) -> str:
    """
    Takes a path such as C:/work/fem.bdf and locates the file using
    relative paths.  If it's on another drive, the path is not modified.

    Parameters
    ----------
    filename : str
        a filename string
    relpath: bool; default=True
       should the relative path be returned

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
                          dict_of_vars: dict[str, Any],
                          log: Any) -> dict[str, Any]:
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


def _get_card_name(lines: list[str], active_filename: str) -> Optional[str]:
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

        # dmig, dmi, dmij, dmiji, dmik, dmiax
        try:
            card_dict = getattr(model, card_name.lower())
        except AttributeError:
            raise AttributeError(f'No slot to store {card_name.lower()}')

        # find the header for the MGG matrix (or other)
        try:
            card = card_dict[name]
        except KeyError:
            raise KeyError(f'No matrix header for {card_name}: {name}')

        # if 0:
        #     if card_name == 'DMIG':
        #         # if field2 == 'UACCEL':  # special DMIG card
        #         card = model.dmig[name]
        #     elif card_name == 'DMI':
        #         card = model.dmi[name]
        #     elif card_name == 'DMIJ':
        #         card = model.dmij[name]
        #     elif card_name == 'DMIJI':
        #         card = model.dmiji[name]
        #     elif card_name == 'DMIK':
        #         card = model.dmik[name]
        #     elif card_name == 'DMIAX':
        #         card = model.dmiax[name]
        #     else:  # pragma: no cover
        #         raise NotImplementedError(card_name)

        for (card_obj, comment) in card_comments:
            card._add_column(card_obj, comment=comment)
        card.finalize()

    # empty the _dmig_temp variable
    model._dmig_temp = defaultdict(list)


def _prep_comment(comment: str) -> str:
    """cleans up the comment"""
    return comment.rstrip()
    # print('comment = %r' % comment)
    # comment = '  this\n  is\n  a comment\n'
    # print(comment.rstrip('\n').split('\n'))
    # sline = [comment[1:] if len(comment) and comment[0] == ' ' else comment
    #          for comment in comment.rstrip().split('\n')]
    # print('sline = ', sline)
