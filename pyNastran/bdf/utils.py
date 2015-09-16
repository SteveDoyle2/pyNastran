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
from six import iteritems, StringIO, string_types
import os
from numpy import unique, cross, dot, array

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
    '$FLUTTER', '$DYNAMIC', '$OPTIMIZATION',
    '$COORDS', '$THERMAL', '$TABLES', '$RANDOM TABLES',
    '$SETS', '$CONTACT', '$REJECTS', '$REJECT_LINES',
    '$PROPERTIES_MASS', '$MASSES',
]

def _clean_comment(comment, end=-1):
    """
    Removes specific pyNastran comment lines so duplicate lines aren't
    created.

    :param comment: the comment to possibly remove
    :param end: lets you remove trailing characters (e.g. a ``\n``; default=-1)
    """
    if comment[:end] in _REMOVED_LINES:
        comment = ''
    elif 'pynastran' in line.lower():
        comment = ''
    return comment


class CardParseSyntaxError(SyntaxError):
    """
    Class that is used for testing.
    Users should just treat this as a SyntaxError.
    """
    pass


def _to_fields_mntpnt1(card_lines, card_name):
    print(card_lines)
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

def to_fields(card_lines, card_name):
    """
    Converts a series of lines in a card into string versions of the field.
    Handles large, small, and CSV formatted cards.

    :param lines:     the lines of the BDF card object
    :param card_name: the card_name -> 'GRID'
    :returns fields:  the string formatted fields of the card

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
        return _to_fields_mntpnt1(card_lines, card_name)

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
            new_fields = line[:72].split(',')[:5]
            for i in range(5 - len(new_fields)):
                new_fields.append('')
        else:  # standard
            new_fields = [line[0:8], line[8:24], line[24:40], line[40:56],
                          line[56:72]]
        fields += new_fields
        assert len(fields) == 5, fields
    else:  # small field
        if ',' in line:  # csv
            new_fields = line[:72].split(',')[:9]
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
                new_fields = line[:72].split(',')[1:5]
                for i in range(4 - len(new_fields)):
                    new_fields.append('')
            else:  # standard
                new_fields = [line[8:24], line[24:40], line[40:56], line[56:72]]
            assert len(new_fields) == 4, new_fields
        else:  # small field
            if ',' in line:  # csv
                new_fields = line[:72].split(',')[1:9]
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

    :param card_lines:  the list of lines in the include card (all the lines!)
    :param include_dir: the include directory (default='')
    :returns filename:  the INCLUDE filename
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
    filename = os.path.join(include_dir, filename)
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
    lline = line[1:].lower()
    if 'pynastran' in lline:
        base, word = lline.split(':')
        if base.strip() != 'pynastran':
            msg = 'unrecognized pyNastran marker\n'
            msg += 'line=%r' % line
            raise SyntaxError(msg)
        key, value = word.strip().split('=')
        key = key.strip()
        value = value.strip()
        if key in ['version', 'encoding', 'punch']:
            pass
        else:
            msg = '\nunrecognized pyNastran key=%r\n' % key
            msg += 'line=%r' % line
            print(msg)
            raise KeyError(msg)

        assert ' ' not in value, 'value=%r' % value
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

    :param filename: a filename string
    :returns filename_string: a shortened representation of the filename
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

    :param node_sets: the node_set to parse
    :returns nodes: list of integers

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

    Example
    -------
    .. code-block:: python

      >>> node_sets = "1 2 3 5:10 12:20:2"
      >>> data = parse_patran_syntax(node_sets)
      >>> data
      data = [1, 2, 3, 5, 6, 7, 8, 9, 10, 12, 14, 16, 18, 20]

    .. warning::
      Don't include the n/node or e/element or any other identifier,
      just a string of "1 2 3 5:10 12:20:2".
      Use parse_patran_syntax_dict to consider the identifier.
    """
    assert isinstance(node_sets, string_types), type(node_sets)
    if pound is not None:
        assert isinstance(pound, (str, int)), type(pound)
        node_sets = node_sets.replace('#', str(pound).strip())

    snodes = node_sets.split()
    nodes = []
    for snode in snodes:
        if ':' in snode:
            ssnode = snode.split(':')
            if len(ssnode) == 2:
                nmin, nmax = int(ssnode[0]), int(ssnode[1])
                new_set = list(range(nmin, nmax + 1))
            elif len(ssnode) == 3:
                nmin, nmax, delta = int(ssnode[0]), int(ssnode[1]), int(ssnode[2])
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


def parse_patran_syntax_dict(node_sets):
    """
    Parses Patran's syntax for compressing nodes/elements

    :param node_sets: the node_set to parse
    :returns nodes: list of integers

    .. code-block:: python

      node_sets = "e 1:3 n 2:6:2 Node 10:13"
      data = parse_patran_syntax_dict(node_sets)
      data = {
          'e'    : [1, 2, 3],
          'n'    : [2, 4, 6],
          'Node' : [10, 11, 12, 13],
      }

    .. note::
      the identifier (e.g. "e") must be used.  Use parse_patran_syntax to skip the identifier.
    .. note::
      doesn't support "1:#"
    """
    data = {}
    try:
        snodes = node_sets.split()
    except TypeError:
        print('node_sets =', node_sets, type(node_sets))
        raise
    key = None
    for snode in snodes:
        if ':' in snode:
            ssnode = snode.split(':')
            if len(ssnode) == 2:
                nmin, nmax = int(ssnode[0]), int(ssnode[1])
                new_set = list(range(nmin, nmax + 1))
            elif len(ssnode) == 3:
                nmin, nmax, delta = int(ssnode[0]), int(ssnode[1]), int(ssnode[2])
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

    :param xyz:    the position of the GRID in an arbitrary
                   coordinate system
    :type xyz:     TYPE = NDARRAY.  SIZE=(3,)
    :param cid:    the coordinate ID for xyz
    :type cid:     int
    :param model:  the BDF model object
    :type model:   BDF()

    :returns xyz2:  the position of the GRID in an arbitrary coordinate system
    :type xyz2:     TYPE = NDARRAY.  SIZE=(3,)
    """
    if is_cid_int:
        cp = model.Coord(cid)
    else:
        cp = cid
    xyz2, matrix = cp.transformToGlobal(xyz)
    return xyz2


def TransformLoadWRT(F, M, cid, cid_new, model, is_cid_int=True):
    """
    Transforms a force/moment from an arbitrary coordinate system to another
    coordinate system.

    :param Fxyz:     the force in an arbitrary coordinate system
    :type Fxyz:      TYPE = NDARRAY.  SIZE=(3,)
    :param Mxyz:     the moment in an arbitrary coordinate system
    :type Mxyz:      TYPE = NDARRAY.  SIZE=(3,)
    :param cid:      the coordinate ID for xyz
    :type cid:       int
    :param cid_new:  the desired coordinate ID
    :type cid_new:   int
    :param model:    the BDF model object
    :type model:     BDF()
    :param is_cid_int:  is cid/cid_new an integer or a Coord object
    :type is_cid_int:  bool

    :returns Fxyz_local:  the force in an arbitrary coordinate system
    :type Fxyz_local:     TYPE = NDARRAY.  SIZE=(3,)
    :returns Mxyz_local:  the force in an arbitrary coordinate system
    :type MxyMxyz_local:  TYPE = NDARRAY.  SIZE=(3,)
    """
    if cid == cid_new: # same coordinate system
        return F, M

    # find the vector r for doing:
    #     M = r x F
    if is_cid_int:
        cp = model.Coord(cid)
        coord_to = model.Coord(cid_new)
    else:
        cp = cid
        coord_to = cid_new
    r = cp.origin - coord_to.origin

    # change R-theta-z to xyz
    Fxyz_local_1 = cp.coordToXYZ(F)
    Mxyz_local_1 = cp.coordToXYZ(M)

    # pGlobal = pLocal1 * beta1 + porigin1
    # pGlobal = pLocal2 * beta2 + porigin2
    # pLocal1 * beta1 + porigin1 = pLocal2 * beta2 + porigin2
    # plocal1 * beta1 + porigin1 - porigin2 = plocal2 * beta2
    # (plocal1 * beta1 + porigin1 - porigin2) * beta2.T = plocal2
    #
    # origin transforms only apply to nodes, so...
    # Fglobal = Flocal1 * beta1
    # Flocal2 = (Flocal1 * beta1) * beta2.T

    Fxyz_global = dot(Fxyz_local_1, cp.beta())
    Fxyz_local_2 = dot(dot(Fxyz_local_1, cp.beta()), coord_to.beta().T)

    # find the moment about the new origin due to the force
    Mxyz_global = cross(r, Fxyz_global)
    dMxyz_local_2 = cross(r, Fxyz_local_2)
    Mxyz_local_2 = Mxyz_local_1 + dMxyz_local_2

    # rotate the delta moment into the local frame
    M_local = coord_to.XYZtoCoord(Mxyz_local_2)

    return Fxyz_local_2, Mxyz_local_2


def PositionWRT(xyz, cid, cid_new, model, is_cid_int=True):
    """
    Gets the location of the GRID which started in some arbitrary system and
    returns it in the desired coordinate system

    :param xyz:      the position of the GRID in an arbitrary
                     coordinate system
    :type xyz:       TYPE = NDARRAY.  SIZE=(3,)
    :param cid:      the coordinate ID for xyz
    :type cid:       int
    :param cid_new:  the desired coordinate ID
    :type cid_new:   int
    :param model:    the BDF model object
    :type model:     BDF()

    :returns xyz_local:  the position of the GRID in an arbitrary coordinate system
    :type xyz_local:     TYPE = NDARRAY.  SIZE=(3,)
    """
    if cid == cid_new: # same coordinate system
        return xyz

    if is_cid_int:
        cp = model.Coord(cid)
        coord_to = model.Coord(cid_new)
    else:
        cp = cid
        coord_to = cid_new

    if 0:
        # pGlobal = pLocal1 * beta1 + porigin1
        # pGlobal = pLocal2 * beta2 + porigin2
        # pLocal1 * beta1 + porigin1 = pLocal2 * beta2 + porigin2
        # plocal1 * beta1 + porigin1 - porigin2 = plocal2 * beta2
        # (plocal1 * beta1 + porigin1 - porigin2) * beta2.T = plocal2

        # convert R-Theta-Z_1 to xyz_1
        p1_local = cp.coordToXYZ(xyz)

        # transform xyz_1 to xyz_2
        p2_local = dot(dot(p1_local, cp.beta()) + cp.origin - coord_to.origin, coord_to.beta().T)

        # convert xyz_2 to R-Theta-Z_2
        xyz_local = coord_to.XYZtoCoord(p2_local)
    else:
        # converting the xyz point arbitrary->global
        xyz_global, matrix_dum = cp.transformToGlobal(xyz)

        # a matrix global->local matrix is found
        matrix = coord_to.beta()
        xyz_local = coord_to.transformToLocal(xyz_global, matrix)
    return xyz_local
