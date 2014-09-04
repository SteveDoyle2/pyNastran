from numpy import unique, cross

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
    for key, ints in data.iteritems():
        data[key] = unique(ints)
    return data

def Position(xyz, cid, model):
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
    cp = model.Coord(cid)
    xyz2, matrix = cp.transformToGlobal(xyz)
    return xyz2

def TransformLoadWRT(Fxyz, Mxyz, cid, cid_new, model):
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

    :returns Fxyz_local:  the force in an arbitrary coordinate system
    :type Fxyz_local:     TYPE = NDARRAY.  SIZE=(3,)
    :returns Mxyz_local:  the force in an arbitrary coordinate system
    :type MxyMxyz_local:  TYPE = NDARRAY.  SIZE=(3,)
    """
    if cid == cid_new: # same coordinate system
        return Fxyz, Mxyz

    # find the vector r for doing:
    #     M = r x F
    cp = model.Coord(cid)
    coordB = model.Coord(cid_new)
    r = cp.origin - coordB.origin

    Fxyz_global = Position(Fxyz, cid, model)
    Fxyz_local = PositionWRT(Fxyz_global, 0, cid_new, model)
    Mxyz_local = PositionWRT(Mxyz, cid, cid_new, model)

    # find the moment about the new origin due to the force
    Mxyz_global_delta += cross(r, Fxyz_global)

    # rotate the delta moment into the local frame
    Mxyz_local_delta = PositionWRT(Mxyz_global_delta, cid, cid_new, model)

    # add the delta load
    Mxyz_local += Mxyz_local_delta

    return Fxyz_local, Mxyz_local

def PositionWRT(xyz, cid, cid_new, model):
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

    # converting the xyz point arbitrary->global
    cp = model.Coord(cid)
    xyz_global, matrixDum = cp.transformToGlobal(xyz)
    coordB = model.Coord(cid_new)

    # a matrix global->local matrix is found
    matrix = coordB.beta()
    xyz_local = coordB.transformToLocal(xyz_global, matrix)
    return xyz_local
