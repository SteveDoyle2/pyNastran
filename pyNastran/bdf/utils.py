"""
Defines various utilities for BDF parsing including:
 - to_fields
 - parse_patran_syntax
 - parse_patran_syntax_dict
 - Position
 - PositionWRT
 - TransformLoadWRT
"""
from __future__ import print_function, unicode_literals
import os
import sys
import inspect
import warnings
from copy import deepcopy
from six import iteritems, StringIO, string_types

import numpy as np
from numpy import unique, cross, dot, array

import pyNastran
from pyNastran.bdf.errors import CardParseSyntaxError
from pyNastran.bdf.cards.collpase_card import collapse_colon_packs


is_windows = 'nt' in os.name
#is_linux = 'posix' in os.name
#is_mac = 'darwin' in os.name

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


def _clean_comment(comment, end=-1):
    """
    Removes specific pyNastran comment lines so duplicate lines aren't
    created.

    Parameters
    ----------
    comment : str
         the comment to possibly remove
    end : int; default=-1
        lets you remove trailing characters (e.g. a ``\n``)
    """
    raise RuntimeError('is this used...')
    if comment[:end] in _REMOVED_LINES:
        comment = ''
    elif 'pynastran' in comment.lower():
        comment = ''
    return comment


def _to_fields_mntpnt1(card_lines):
    assert len(card_lines) == 2, card_lines
    line1, line2 = card_lines

    base = line1[:24]
    comment = line1[24:]  # len=56 max
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

def to_fields_long(card_lines, card_name):
    """
    Converts a series of lines in a card into string versions of the field.
    Handles large, small, and CSV formatted cards.

    Doesn't consider Nastran's rule about 72 character width fields,
    which is nice when you have poorly formatted BDFs.

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
      >>> fields = to_fields_long(lines, card_name)
      >>> fields
      ['GRID', '1', '', '1.0', '2.0', '3.0']
    """
    fields = []

    if card_name == 'MONPNT1':
        return _to_fields_mntpnt1(card_lines)

    # first line
    line = card_lines.pop(0)
    if '=' in line:
        msg = 'card_name=%r\nequal signs are not supported...line=%r' % (card_name, line)
        raise CardParseSyntaxError(msg)

    if '\t' in line:
        line = line.expandtabs()
        if ',' in line:
            msg = 'tabs and commas in the same line are not supported...line=%r' % line
            raise CardParseSyntaxError(msg)

    if '*' in line:  # large field
        if ',' in line:  # csv
            new_fields = line.split(',')[:5]
            for i in range(5 - len(new_fields)):
                new_fields.append('')
        else:  # standard
            new_fields = [line[0:8], line[8:24], line[24:40], line[40:56],
                          line[56:72]]
        fields += new_fields
        assert len(fields) == 5, fields
    else:  # small field
        if ',' in line:  # csv
            new_fields = line.split(',')[:9]
            for i in range(9 - len(new_fields)):
                new_fields.append('')
        else:  # standard
            new_fields = [line[0:8], line[8:16], line[16:24], line[24:32],
                          line[32:40], line[40:48], line[48:56], line[56:64],
                          line[64:72]]
        fields += new_fields
        assert len(fields) == 9, fields

    for j, line in enumerate(card_lines): # continuation lines
        if '=' in line and card_name != 'EIGRL':
            msg = 'card_name=%r\nequal signs are not supported...line=%r' % (card_name, line)
            raise CardParseSyntaxError(msg)
        if '\t' in line:
            line = line.expandtabs()
            if ',' in line:
                msg = 'tabs and commas in the same line are not supported...line=%r' % line
                raise CardParseSyntaxError(msg)

        if '*' in line:  # large field
            if ',' in line:  # csv
                new_fields = line.split(',')[1:5]
                for i in range(4 - len(new_fields)):
                    new_fields.append('')
            else:  # standard
                new_fields = [line[8:24], line[24:40], line[40:56], line[56:72]]
            assert len(new_fields) == 4, new_fields
        else:  # small field
            if ',' in line:  # csv
                new_fields = line.split(',')[1:9]
                for i in range(8 - len(new_fields)):
                    new_fields.append('')
            else:  # standard
                new_fields = [line[8:16], line[16:24], line[24:32],
                              line[32:40], line[40:48], line[48:56],
                              line[56:64], line[64:72]]
            if len(new_fields) != 8:
                nfields = len(new_fields)
                msg = 'nFields=%s new_fields=%s' % (nfields, new_fields)
                raise RuntimeError(msg)

        fields += new_fields
    return fields #[field.strip() for field in fields]

def to_fields(card_lines, card_name):
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
    fields = []

    if card_name == 'MONPNT1':
        return _to_fields_mntpnt1(card_lines)

    # first line
    line = card_lines.pop(0)
    if '=' in line:
        msg = 'card_name=%r\nequal signs are not supported...line=%r' % (card_name, line)
        raise CardParseSyntaxError(msg)

    if '\t' in line:
        line = line.expandtabs()
        if ',' in line:
            msg = 'tabs and commas in the same line are not supported...line=%r' % line
            raise CardParseSyntaxError(msg)

    if '*' in line:  # large field
        if ',' in line:  # csv
            new_fields = line.split(',')[:5]
            for i in range(5 - len(new_fields)):
                new_fields.append('')
        else:  # standard
            new_fields = [line[0:8], line[8:24], line[24:40], line[40:56],
                          line[56:72]]
        fields += new_fields
        assert len(fields) == 5, fields
    else:  # small field
        if ',' in line:  # csv
            new_fields = line.split(',')[:9]
            for i in range(9 - len(new_fields)):
                new_fields.append('')
        else:  # standard
            new_fields = [line[0:8], line[8:16], line[16:24], line[24:32],
                          line[32:40], line[40:48], line[48:56], line[56:64],
                          line[64:72]]
        fields += new_fields
        assert len(fields) == 9, fields

    for j, line in enumerate(card_lines): # continuation lines
        if '=' in line and card_name != 'EIGRL':
            msg = 'card_name=%r\nequal signs are not supported...line=%r' % (card_name, line)
            raise CardParseSyntaxError(msg)
        if '\t' in line:
            line = line.expandtabs()
            if ',' in line:
                msg = 'tabs and commas in the same line are not supported...line=%r' % line
                raise CardParseSyntaxError(msg)

        if '*' in line:  # large field
            if ',' in line:  # csv
                new_fields = line.split(',')[1:5]
                for i in range(4 - len(new_fields)):
                    new_fields.append('')
            else:  # standard
                new_fields = [line[8:24], line[24:40], line[40:56], line[56:72]]
            assert len(new_fields) == 4, new_fields
        else:  # small field
            if ',' in line:  # csv
                new_fields = line.split(',')[1:9]
                for i in range(8 - len(new_fields)):
                    new_fields.append('')
            else:  # standard
                new_fields = [line[8:16], line[16:24], line[24:32],
                              line[32:40], line[40:48], line[48:56],
                              line[56:64], line[64:72]]
            if len(new_fields) != 8:
                nfields = len(new_fields)
                msg = 'nFields=%s new_fields=%s' % (nfields, new_fields)
                raise RuntimeError(msg)

        fields += new_fields
    return fields #[field.strip() for field in fields]


def get_include_filename(card_lines, include_dir=''):
    """
    Parses an INCLUDE file split into multiple lines (as a list).

    Parameters
    ----------
    card_lines : List[str]
        the list of lines in the include card (all the lines!)
    include_dir : str; default=''
        the include directory

    Returns
    -------
    filename : str
        the INCLUDE filename
    """
    card_lines2 = []
    for line in card_lines:
        line = line.strip('\t\r\n ')
        card_lines2.append(line)

    card_lines2[0] = card_lines2[0][7:].strip()  # truncate the cardname
    filename = ''.join(card_lines2)
    filename = filename.strip('"').strip("'")
    if ':' in filename:
        ifilepaths = filename.split(':')
        filename = os.path.join(*ifilepaths)
    if is_windows:
        filename = os.path.join(include_dir, filename).replace('/', '\\')
    else:
        filename = os.path.join(include_dir, filename).replace('\\', '/')
    return filename


def parse_executive_control_deck(executive_control_lines):
    """
    Extracts the solution from the executive control deck
    """
    sol = None
    method = None
    isol_line = None
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
            isol_line = i
    return sol, method, isol_line


def _parse_pynastran_header(line):
    """
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


def clean_empty_lines(lines):
    """
    Removes leading and trailing empty lines
    don't remove internally blank lines
    """
    found_lines = False
    if len(lines) < 2:
        return lines

    for i, line in enumerate(lines):
        if not found_lines and line:
            found_lines = True
            n1 = i
            n2 = i + 1
        elif found_lines and line:
            n2 = i + 1
    lines2 = lines[n1:n2]
    return lines2


def print_filename(filename, relpath):
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


def parse_patran_syntax(node_sets, pound=None):
    """
    Parses Patran's syntax for compressing nodes/elements

    Parameters
    ----------
    node_sets : str
        the node_set to parse
    pound : int / str
        value : the pound value (e.g. # in 1:#, which means all)

    Returns
    -------
    nodes : List[int]
        the integer values

    Patran has a short syntax of the form:

      +------------+----------------------+
      |  String    | Output               |
      +------------+----------------------+
      |"1 2 3"     | [1, 2, 3]            |
      +------------+----------------------+
      |"5:10"      | [5, 6, 7, 8, 9, 10]  |
      +------------+----------------------+
      |"12:20:2"   | [12, 14, 16, 18, 20] |
      +------------+----------------------+

    Example 1
    ----------
    .. code-block:: python

      >>> node_sets = "1 2 3 5:10 12:20:2"
      >>> data = parse_patran_syntax(node_sets)
      >>> data
      data = [1, 2, 3, 5, 6, 7, 8, 9, 10, 12, 14, 16, 18, 20]

    Example 2
    ----------
    .. code-block:: python

      >>> node_sets = "1 2 3:#"
      >>> data = parse_patran_syntax(node_sets, pound=10)
      >>> data
      data = [1, 2, 3, 5, 6, 7, 8, 9, 10]


    .. warning::  Don't include the n/node or e/element or any other
                  identifier, just a string of "1 2 3 5:10 12:20:2".
                  Use parse_patran_syntax_dict to consider the identifier.
    """
    assert isinstance(node_sets, string_types), type(node_sets)
    if pound is not None:
        assert isinstance(pound, (string_types, int)), type(pound)
        node_sets = node_sets.replace('#', str(pound).strip())
    if len(node_sets) == 0:
        return array([], dtype='int32')

    snodes = node_sets.split()
    nodes = []
    for snode in snodes:
        if ':' in snode:
            ssnode = snode.split(':')
            if len(ssnode) == 2:
                nmin = int(ssnode[0])
                nmax = int(ssnode[1])
                new_set = list(range(nmin, nmax + 1))
            elif len(ssnode) == 3:
                nmin = int(ssnode[0])
                nmax = int(ssnode[1])
                delta = int(ssnode[2])
                nmin, nmax = min([nmin, nmax]), max([nmin, nmax])
                if delta > 0:
                    new_set = list(range(nmin, nmax + 1, delta))
                else:
                    new_set = list(range(nmin, nmax + 1, -delta))
            else:
                raise NotImplementedError(snode)
            nodes += new_set
        else:
            nodes.append(int(snode))
    return unique(nodes)

def write_patran_syntax_dict(dict_sets):
    msg = ''
    for key, dict_set in iteritems(dict_sets):
        singles, doubles = collapse_colon_packs(dict_set)
        double_list = ['%s:%s ' % (double[0], double[2])
                       if len(double) == 3 else '%s:%s:%s ' % (double[0], double[2], double[4])
                       for double in doubles]
        double_str = ''.join(double_list)
        msg += '%s %s %s' % (key,
                             ' '.join(str(single) for single in singles),
                             double_str,
                            )
    assert '%' not in msg, msg
    return msg


def parse_patran_syntax_dict(node_sets, pound_dict=None):
    """
    Parses Patran's syntax for compressing nodes/elements

    Parameters
    ----------
    node_sets : str
        the node_set to parse
    pound_dict : List[str] : int
        key : the string
        value : the pound value (e.g. 1:#)

    Returns
    -------
    nodes : Dict[str] = List[int]
        str : the key
        values : the integer values for that key

    Example 1
    ---------
    .. code-block:: python

      node_sets = "e 1:3 n 2:6:2 Node 10:13"
      data = parse_patran_syntax_dict(node_sets)
      data = {
          'e'    : [1, 2, 3],
          'n'    : [2, 4, 6],
          'Node' : [10, 11, 12, 13],
      }

    Example 2
    ---------
    .. code-block:: python

      node_sets = "e 1:3 n 2:6:2 Node 10:#"

      # a pound character will be set to 20, but only for 'Node', but not 'n'
      # so define it twice if needed
      pounds = {'Node' : 20}
      data = parse_patran_syntax_dict(node_sets, pounds=pounds)
      data = {
          'e'    : [1, 2, 3],
          'n'    : [2, 4, 6],
          'Node' : [10, 11, 12, 13],
      }

    .. note:: an identifier (e.g. "e") must be used.
              Use parse_patran_syntax to skip the identifier.
    .. warning:: case sensitive
    """
    data = {}
    try:
        snodes = node_sets.split()
    except AttributeError:
        print('node_sets =', node_sets, type(node_sets))
        raise
    except TypeError:
        print('node_sets =', node_sets, type(node_sets))
        raise

    if pound_dict is None:
        pound_dict = {}

    key = None
    for snode in snodes:
        if ':' in snode:
            ssnode = snode.split(':')
            if len(ssnode) == 2:
                if ssnode[0].isdigit():
                    nmin = int(ssnode[0])
                else:
                    raise NotImplementedError('ssnode=%s must be int,int' % ssnode)

                if ssnode[1].isdigit():
                    nmax = int(ssnode[1])
                elif ssnode[1] == '#' and key in pound_dict:
                    nmax = int(pound_dict[key])
                else:
                    raise NotImplementedError('ssnode=%s must be int,int' % ssnode)

                new_set = list(range(nmin, nmax + 1))
            elif len(ssnode) == 3:
                if ssnode[0].isdigit():
                    nmin = int(ssnode[0])
                else:
                    raise NotImplementedError('ssnode=%s must be int,int,int' % ssnode)

                if ssnode[1].isdigit():
                    nmax = int(ssnode[1])
                elif ssnode[1] == '#' and key in pound_dict:
                    nmax = int(pound_dict[key])
                else:
                    raise NotImplementedError('ssnode=%s must be int,int,int' % ssnode)

                delta = int(ssnode[2])
                nmin, nmax = min([nmin, nmax]), max([nmin, nmax])
                if delta > 0:
                    new_set = list(range(nmin, nmax + 1, delta))
                else:
                    new_set = list(range(nmin, nmax + 1, -delta))
            else:
                raise NotImplementedError(snode)
            if key is None:
                msg = 'data must be of the form "Node 10:13", not "10:13"\n'
                msg += 'new_set=%s' % array(new_set, dtype='int32')
                raise SyntaxError(msg)
            data[key] += new_set
        else:
            if snode.isdigit():
                data[key].append(int(snode))
            else:
                key = snode
                if key is None:
                    msg = 'data must be of the form "Node 10:13", not "10:13"'
                    raise SyntaxError(msg)
                if key not in data:
                    data[key] = []
    for key, ints in iteritems(data):
        data[key] = unique(ints)
    return data


def Position(xyz, cid, model, is_cid_int=True):
    """
    Gets the point in the global XYZ coordinate system.

    Parameters
    ----------
    xyz : (3,) ndarray
        the position of the GRID in an arbitrary coordinate system
    cid : int
        the coordinate ID for xyz
    model : BDF()
        the BDF model object

    Returns
    -------
    xyz2 : (3,) ndarray
        the position of the GRID in an arbitrary coordinate system
    """
    if is_cid_int:
        cp_ref = model.Coord(cid)
    else:
        cp_ref = cid
    xyz2 = cp_ref.transform_node_to_global(xyz)
    return xyz2


def TransformLoadWRT(F, M, cid, cid_new, model, is_cid_int=True):
    """
    Transforms a force/moment from an arbitrary coordinate system to another
    coordinate system.

    Parameters
    ----------
    Fxyz : (3, ) float ndarray
        the force in an arbitrary coordinate system
    Mxyz : (3, ) float ndarray
        the moment in an arbitrary coordinate system
    cid : int
        the coordinate ID for xyz
    cid_new : int
        the desired coordinate ID
    model : BDF()
        the BDF model object
    is_cid_int : bool
        is cid/cid_new an integer or a Coord object

    Returns
    -------
    Fxyz_local : (3, ) float ndarray
        the force in an arbitrary coordinate system
    Mxyz_local : (3, ) float ndarray
        the force in an arbitrary coordinate system
    """
    if cid == cid_new: # same coordinate system
        return F, M

    # find the vector r for doing:
    #     M = r x F
    if is_cid_int:
        cp_ref = model.Coord(cid)
        coord_to_ref = model.Coord(cid_new)
    else:
        cp_ref = cid
        coord_to_ref = cid_new
    r = cp_ref.origin - coord_to_ref.origin

    # change R-theta-z to xyz
    Fxyz_local_1 = cp_ref.coord_to_xyz(F)
    Mxyz_local_1 = cp_ref.coord_to_xyz(M)

    # pGlobal = pLocal1 * beta1 + porigin1
    # pGlobal = pLocal2 * beta2 + porigin2
    # pLocal1 * beta1 + porigin1 = pLocal2 * beta2 + porigin2
    # plocal1 * beta1 + porigin1 - porigin2 = plocal2 * beta2
    # (plocal1 * beta1 + porigin1 - porigin2) * beta2.T = plocal2
    #
    # origin transforms only apply to nodes, so...
    # Fglobal = Flocal1 * beta1
    # Flocal2 = (Flocal1 * beta1) * beta2.T

    Fxyz_global = dot(Fxyz_local_1, cp_ref.beta())
    Fxyz_local_2 = dot(dot(Fxyz_local_1, cp_ref.beta()), coord_to_ref.beta().T)

    # find the moment about the new origin due to the force
    Mxyz_global = cross(r, Fxyz_global)
    dMxyz_local_2 = cross(r, Fxyz_local_2)
    Mxyz_local_2 = Mxyz_local_1 + dMxyz_local_2

    # rotate the delta moment into the local frame
    M_local = coord_to_ref.xyz_to_coord(Mxyz_local_2)

    return Fxyz_local_2, Mxyz_local_2


def PositionWRT(xyz, cid, cid_new, model, is_cid_int=True):
    """
    Gets the location of the GRID which started in some arbitrary system and
    returns it in the desired coordinate system

    Parameters
    ----------
    xyz : (3, ) float ndarray
        the position of the GRID in an arbitrary coordinate system
    cid : int
        the coordinate ID for xyz
    cid_new : int
        the desired coordinate ID
    model : BDF()
        the BDF model object
    is_cid_int : bool
        is cid/cid_new an integer or a Coord object

    Returns
    -------
    xyz_local : (3, ) float ndarray
        the position of the GRID in an arbitrary coordinate system
    """
    if cid == cid_new: # same coordinate system
        return xyz

    if is_cid_int:
        cp_ref = model.Coord(cid)
        coord_to_ref = model.Coord(cid_new)
    else:
        cp_ref = cid
        coord_to_ref = cid_new

    if 0:
        # pGlobal = pLocal1 * beta1 + porigin1
        # pGlobal = pLocal2 * beta2 + porigin2
        # pLocal1 * beta1 + porigin1 = pLocal2 * beta2 + porigin2
        # plocal1 * beta1 + porigin1 - porigin2 = plocal2 * beta2
        # (plocal1 * beta1 + porigin1 - porigin2) * beta2.T = plocal2

        # convert R-Theta-Z_1 to xyz_1
        p1_local = cp_ref.coord_to_xyz(xyz)

        # transform xyz_1 to xyz_2
        p2_local = dot(
            dot(p1_local, cp_ref.beta()) + cp_ref.origin - coord_to_ref.origin,
            coord_to_ref.beta().T)

        # convert xyz_2 to R-Theta-Z_2
        xyz_local = coord_to_ref.xyz_to_coord(p2_local)
    else:
        # converting the xyz point arbitrary->global
        xyz_global = cp_ref.transform_node_to_global(xyz)

        # a matrix global->local matrix is found
        matrix = coord_to_ref.beta()
        xyz_local = coord_to_ref.transformToLocal(xyz_global, matrix)
    return xyz_local



def deprecated(old_name, new_name, deprecated_version, levels=None):
    """
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

    new_line = ''
    if new_name:
        new_line = "; replace it with '%s'\n" % new_name
    msg = "'%s' was deprecated in v%s%s" % (old_name, deprecated_version, new_line)

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
        print(code)

        source_lines, line_no0 = inspect.getsourcelines(code)
        di = line_no - line_no0
        try:
            line = source_lines[di]
        except:
            break
        msg += '  %-25s lineNo=%-4s %s\n' % (filename, str(line_no) + ';', line.strip())

    if ver_tuple > dep_ver_tuple:
        # fail
        raise NotImplementedError(msg)
    else:
        warnings.warn(msg, DeprecationWarning)

def split_eids_along_nids(model, eids, nids):
    """
    Dissassociate a list of elements along a list of nodes.

    The expected use of this function is that you have two bodies that
    are incorrectly equivalenced and you would like to create duplicate
    nodes at the same location and associate the new nodes with one half
    of the elements.

    Pick the nodes along the line and the elements along one side of the line.

    Parameters
    ----------
    model : BDF()
        the BDF model
    eids : list/tuple
        element ids to disassociate
    nids : list/tuple
        node ids to disassociate

    Implicitly returns model with additional nodes.

    .. note :: xref should be set to False for this function.
    """
    #assert model.xref == False, model.xref
    nid = max(model.nodes.keys()) + 1

    nid_map = {}
    for nidi in nids:
        node = model.nodes[nidi]
        node2 = deepcopy(node)
        node2.nid = nid
        model.nodes[nid] = node2
        nid_map[nidi] = nid
        nid += 1

    for eid in eids:
        nodes = []
        elem = model.elements[eid]
        for nidi in elem.nodes:
            if nidi in nid_map:
                nodes.append(nid_map[nidi])
            else:
                nodes.append(nidi)
            assert len(np.unique(nodes)) == len(nodes), 'nodes=%s' % nodes
        elem.nodes = nodes
