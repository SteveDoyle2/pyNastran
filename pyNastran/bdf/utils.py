from six import iteritems
import os
from numpy import unique, cross, dot, array


def _clean_comment(comment, end=-1):
    """
    Removes specific pyNastran comment lines so duplicate lines aren't
    created.

    :param comment: the comment to possibly remove
    """
    if comment[:end] in ['$EXECUTIVE CONTROL DECK',
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
            '$PROPERTIES_MASS', '$MASSES']:
        comment = ''
    return comment


class CardParseSyntaxError(SyntaxError):
    pass


def to_fields(card_lines, card_name):
    """
    Converts a series of lines in a card into string versions of the field.
    Handles large, small, and CSV formatted cards.

    :param lines:     the lines of the BDF card object
    :param card_name: the card_name -> 'GRID'
    :returns fields:  the string formatted fields of the card
    """
    fields = []
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
        assert len(fields) == 5
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
        assert len(fields) == 9

    for j, line in enumerate(card_lines): # continuation lines
        #for i, field in enumerate(fields):
        #    if field.strip() == '+':
        #        raise RuntimeError('j=%s field[%s] is a +' % (j,i))

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
            assert len(new_fields) == 4
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
    return [field.strip() for field in fields]


def get_include_filename(card_lines, include_dir=''):
    """
    Parses an INCLUDE file split into multiple lines (as a list).

    :param card_lines:  the list of lines in the include card (all the lines!)
    :param include_dir: the include directory (default='')
    :returns filename:  the INCLUDE filename
    """
    cardLines2 = []
    for line in card_lines:
        line = line.strip('\t\r\n ')
        cardLines2.append(line)

    cardLines2[0] = cardLines2[0][7:].strip()  # truncate the cardname
    filename = ''.join(cardLines2)
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
    iSolLine = None
    for (i, eline) in enumerate(executive_control_lines):
        uline = eline.strip().upper()  # uppercase line
        uline = uline.split('$')[0].expandtabs()
        if uline[:4] in ['SOL ']:
            if ',' in uline:
                sline = uline.split(',')  # SOL 600,method
                solValue = sline[0].strip()
                method = sline[1].strip()
            else:
                solValue = uline
                method = None

            if sol is None:
                sol = solValue[3:].strip()
            else:
                raise ValueError('cannot overwrite solution existing='
                                 '|SOL %s| new =|%s|' % (sol, uline))
            iSolLine = i
    return sol, method, iSolLine


def clean_empty_lines(lines):
    """
    removes leading and trailing empty lines
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
    driveLetter = os.path.splitdrive(os.path.abspath(filename))[0]
    if driveLetter == os.path.splitdrive(os.curdir)[0] and relpath:
        return os.path.relpath(filename)
    return filename


def parse_patran_syntax(node_sets):
    """
    Parses Patran's syntax for compressing nodes/elements

    :param node_sets: the node_set to parse
    :returns nodes: list of integers

    Patran has a short syntax of the form:

      +------------+--------------------+
      |  String    | Output             |
      +------------+--------------------+
      |"1 2 3"   | [1, 2, 3]            |
      |"5:10"    | [5, 6, 7, 8, 9, 10]  |
      |"12:20:2" | [12, 14, 16, 18, 20] |
      +------------+--------------------+

    Example
    -------
    node_sets = "1 2 3 5:10 12:20:2"
    data = parse_patran_syntax(node_sets)

    :warning:
      Don't include the n/node or e/element or any other identifier,
      just a string of "1 2 3 5:10 12:20:2".
      Use parse_patran_syntax_dict to consider the identifier.
    :note:
      doesn't support "1:#"
    """
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

                nmin, nmax, delta = int(ssnode[0]), int(ssnode[1]), int(ssnode[2])
                if delta > 0:
                    new_set = list(range(nmin, nmax + 1, delta))
                else:
                    new_set = list(range(nmax, nmin + 1, -delta))
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

    Example
    -------
    node_sets = "e 1:3 n 2:6:2 Node 10:13"
    data = parse_patran_syntax_dict(node_sets)
    data = {
        'e'    : [1, 2, 3],
        'n'    : [2, 4, 6],
        'Node' : [10, 11, 12, 13],
    }
    :note:
      the identifier (e.g. "e") must be used.  Use parse_patran_syntax to skip the identifier.
    :note:
      doesn't support "1:#"
    """
    data = {}
    snodes = node_sets.split()
    key = None
    for snode in snodes:
        if ':' in snode:
            ssnode = snode.split(':')
            if len(ssnode) == 2:
                nmin, nmax = int(ssnode[0]), int(ssnode[1])
                new_set = list(range(nmin, nmax + 1))
            elif len(ssnode) == 3:
                nmin, nmax, delta = int(ssnode[0]), int(ssnode[1]), int(ssnode[2])
                if delta > 0:
                    new_set = list(range(nmin, nmax + 1, delta))
                else:
                    new_set = list(range(nmax, nmin + 1, -delta))
            else:
                raise NotImplementedError(snode)
            if key is None:
                msg =  'data must be of the form "Node 10:13", not "10:13"\n'
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
        coordB = model.Coord(cid_new)
    else:
        cp = cid
        coordB = cid_new
    r = cp.origin - coordB.origin

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
    Fxyz_local_2 = dot(dot(Fxyz_local_1, cp.beta()), coordB.beta().T)

    # find the moment about the new origin due to the force
    Mxyz_global = cross(r, Fxyz_global)
    dMxyz_local_2 = cross(r, Fxyz_local_2)
    Mxyz_local_2 = Mxyz_local_1 + dMxyz_local_2

    # rotate the delta moment into the local frame
    M_local = coordB.XYZtoCoord(Mxyz_local_2)

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
        coordB = model.Coord(cid_new)
    else:
        cp = cid
        coordB = cid_new

    if 0:
        # pGlobal = pLocal1 * beta1 + porigin1
        # pGlobal = pLocal2 * beta2 + porigin2
        # pLocal1 * beta1 + porigin1 = pLocal2 * beta2 + porigin2
        # plocal1 * beta1 + porigin1 - porigin2 = plocal2 * beta2
        # (plocal1 * beta1 + porigin1 - porigin2) * beta2.T = plocal2

        # convert R-Theta-Z_1 to xyz_1
        p1_local = cp.coordToXYZ(xyz)

        # transform xyz_1 to xyz_2
        p2_local = dot(dot(p1_local, cp.beta()) + cp.origin - coordB.origin, coordB.beta().T)

        # convert xyz_2 to R-Theta-Z_2
        xyz_local = coordB.XYZtoCoord(p2_local)
    else:
        # converting the xyz point arbitrary->global
        xyz_global, matrixDum = cp.transformToGlobal(xyz)

        # a matrix global->local matrix is found
        matrix = coordB.beta()
        xyz_local = coordB.transformToLocal(xyz_global, matrix)
    return xyz_local
