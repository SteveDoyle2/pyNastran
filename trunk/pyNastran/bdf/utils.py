from six import iteritems
import os
from numpy import unique, cross, dot

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
    data = parse_patran_syntax_dict(node_sets)

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
                new_set = list(range(nmin, nmax+1))
            elif len(ssnode) == 3:
                nmin, nmax, delta = int(ssnode[0]), int(ssnode[1]), int(ssnode[2])
                new_set = list(range(nmin, nmax+1, delta))
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
    #from collections import defaultdict
    #data = defaultdict(list)
    data = {}
    snodes = node_sets.split()
    key = None
    for snode in snodes:
        if ':' in snode:
            ssnode = snode.split(':')
            if len(ssnode) == 2:
                nmin, nmax = int(ssnode[0]), int(ssnode[1])
                new_set = list(range(nmin, nmax+1))
            elif len(ssnode) == 3:
                nmin, nmax, delta = int(ssnode[0]), int(ssnode[1]), int(ssnode[2])
                new_set = list(range(nmin, nmax+1, delta))
            else:
                raise NotImplementedError(snode)
            data[key] += new_set
        else:
            if snode.isdigit():
                data[key].append(int(snode))
            else:
                key = snode
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

    Fxyz_global  = dot(Fxyz_local_1, cp.beta())
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
