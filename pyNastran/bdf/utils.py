"""
Defines various utilities including:
 - parse_patran_syntax
 - parse_patran_syntax_dict
 - Position
 - PositionWRT
 - TransformLoadWRT
"""
from __future__ import print_function, unicode_literals
from copy import deepcopy
from six import string_types
from typing import List, Union, Dict, Tuple, Optional

import numpy as np  # type: ignore
from numpy import unique, cross, dot, array  # type: ignore

from pyNastran.bdf.cards.collpase_card import collapse_colon_packs
from pyNastran.bdf.bdf_interface.utils import deprecated
from pyNastran.utils.numpy_utils import integer_types


def parse_patran_syntax(node_sets, pound=None):
    # type: (str, Optional[int]) -> np.ndarray
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

    Examples
    --------
    **Example 1**

    >>> node_sets = "1 2 3 5:10 12:20:2"
    >>> data = parse_patran_syntax(node_sets)
    >>> data
    data = [1, 2, 3, 5, 6, 7, 8, 9, 10, 12, 14, 16, 18, 20]

    **Example 2**

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
        assert isinstance(pound, (string_types, integer_types)), type(pound)
        node_sets = node_sets.replace('#', str(pound).strip())
    if len(node_sets) == 0:
        return array([], dtype='int32')

    snodes = node_sets.split()
    nodes = []  # type: List[int]
    for snode in snodes:
        _apply_comma_colon_int_node(nodes, snode)
    return unique(nodes)

def _apply_comma_colon_int_node(nodes, snode):
    """helper method for parse_patran_syntax"""
    if ',' in snode:
        comma_split_node = snode.split(',')
        for comma_node in comma_split_node:
            _apply_comma_colon_int_node(nodes, comma_node)
    elif ':' in snode:
        new_set = _apply_colon_set(snode)
        nodes += new_set
    else:
        nodes.append(int(snode))

def _apply_colon_set(snode):
    """helper method for parse_patran_syntax"""
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
    return new_set

def write_patran_syntax_dict(dict_sets):
    # type: (Dict[str, np.ndarray]) -> str
    """
    writes partran syntax

    Parameters
    ----------
    dict_sets : Dict[str] = List[int]
        str : the key
        values : the integer values for that key

    Returns
    -------
    node_sets : str
        the node_set to parse

    See ``parse_patran_syntax_dict`` for explanation of usage

    """
    msg = ''
    for key, dict_set in sorted(dict_sets.items()):
        singles, doubles = collapse_colon_packs(dict_set, thru_split=4)
        double_list = ('%s:%s' % (double[0], double[2])
                       if len(double) == 3 else '%s:%s:%s' % (double[0], double[2], double[4])
                       for double in doubles)
        double_str = ' '.join(double_list)
        msg += '%s %s %s ' % (
            key,
            ' '.join(str(single) for single in singles),
            double_str,
        )
    assert '%' not in msg, msg
    return msg.strip().replace('  ', ' ')


def parse_patran_syntax_dict(node_sets, pound_dict=None, msg=''):
    # type: (str, Dict[str, Optional[int]], str) -> Dict[str, np.ndarray]
    """
    Parses Patran's syntax for compressing nodes/elements

    Parameters
    ----------
    node_sets : str
        the node_set to parse
    pound_dict : List[str] : int
        key : the string
        value : the pound value (e.g. 1:#)
    msg : str
        error message; currently unused

    Returns
    -------
    nodes : Dict[str] = List[int]
        str : the key
        values : the integer values for that key

    Examples
    --------
    **Example 1**

    >>> node_sets = "e 1:3 n 2:6:2 Node 10:13"
    >>> data = parse_patran_syntax_dict(node_sets)
    >>> data = {
          'e'    : [1, 2, 3],
          'n'    : [2, 4, 6],
          'Node' : [10, 11, 12, 13],
    }


    **Example 2**

    >>> node_sets = "e 1:3 n 2:6:2 Node 10:#"

    # a pound character will be set to 20, but only for 'Node', but not
    # 'n' so define it twice if needed
    >>> pounds = {'Node' : 20}
    >>> data = parse_patran_syntax_dict(node_sets, pounds=pounds)
    >>> data = {
          'e'    : [1, 2, 3],
          'n'    : [2, 4, 6],
          'Node' : [10, 11, 12, 13],
      }

    Notes
    -----
    An identifier (e.g. "e") must be used.
    Use parse_patran_syntax to skip the identifier.

    .. warning:: case sensitive

    """
    data = {}  # type: Dict[str, List[int]]
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
    for key, ints in data.items():
        data[key] = unique(ints)
    return data


def parse_patran_syntax_dict_map(node_sets, type_map, msg=''):
    # type: (str, Dict[str, str], str) -> Dict[str, np.ndarray]
    """
    Parses Patran's syntax for compressing nodes/elements

    Parameters
    ----------
    node_sets : str
        the node_set to parse
    type_map : dict[key_in] : key_out
        key_in : str
            the name of the input string
        key_out : str
            the name of the out string
    #pound_dict : List[str] : int
        #key : the string
        #value : the pound value (e.g. 1:#)
    msg : str
        error message; currently unused

    Returns
    -------
    nodes : Dict[str] = List[int]
        str : the key
        values : the integer values for that key

    Examples
    --------
    **Example 1**
    .. code-block:: python

      # we drop the coordinate systems because we didn't request them
      # (coord is not referenced)
      #
      >>> node_sets = "e 1:3 n 2:6:2 Node 10:13 N 15 coord 1:10"
      >>> type_map = {
          'n' : 'Node',
          'Node' : 'Node',
          'e' : 'Element',
          'Elm' : 'Element',
          'Element' : 'Element',
      }

      **Example 2**
      >>> data = parse_patran_syntax_dict(node_sets, type_map)
      >>> data = {
          'Element' : [1, 2, 3],
          'Node' : [2, 4, 6, 10, 11, 12, 13, 15],
      }

    .. todo:: doesn't support pound_dict

    """
    # makes it so we can pass in 'N' and 'n' and still get 'Node' out
    update_type_map = {}  # type: Dict[str, str]
    for key, value in type_map.items():
        if key in update_type_map:
            assert update_type_map[key] == value
        update_type_map[key.upper()] = value

    dict_in = parse_patran_syntax_dict(node_sets.upper(), pound_dict=None)
    dict_temp = {}  # type: Dict[str, np.ndarray]
    for key_in, value in sorted(dict_in.items()):
        key_in2 = key_in.upper()
        if key_in2 in update_type_map:
            key_out = update_type_map[key_in2]
            #print('key_in=%r key_out=%r' % (key_in, key_out))
            if key_out in dict_temp:
                dict_temp[key_out].append(value)
            else:
                dict_temp[key_out] = [value]
        else:
            print('skipping key=%r while parsing %s' % (key_in, msg))

    dict_out = {}  # type: Dict[str, np.ndarray]
    for key, value_list in dict_temp.items():
        if len(value_list) == 1:
            value = value_list[0]
        else:
            value = np.hstack(value_list)
            value.sort()
        dict_out[key] = value
    return dict_out


def Position(xyz, cid, model, is_cid_int=None):
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
    is_cid_int : bool
        is cid/cid_new an integer or a Coord object (deprecated)

    Returns
    -------
    xyz2 : (3,) ndarray
        the position of the GRID in an arbitrary coordinate system

    """
    if is_cid_int is not None:  # pragma: no cover
        deprecated('Position(xyz, cid, model, is_cid_int=%s)' % is_cid_int,
                   'Position(xyz, cid, model)', '1.2',
                   levels=[-1])

    cp_ref = _coord(model, cid)
    xyz2 = cp_ref.transform_node_to_global(xyz)
    return xyz2

def TransformLoadWRT(F, M, cid, cid_new, model, is_cid_int=None):
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
        is cid/cid_new an integer or a Coord object (deprecated)

    Returns
    -------
    Fxyz_local : (3, ) float ndarray
        the force in an arbitrary coordinate system
    Mxyz_local : (3, ) float ndarray
        the force in an arbitrary coordinate system

    """
    if is_cid_int is not None:  # pragma: no cover
        deprecated('TransformLoadWRT(F, M, cid, cid_new, model, is_cid_int=%s)' % is_cid_int,
                   'TransformLoadWRT(F, M, cid, cid_new, model)', '1.2',
                   levels=[-1])

    if cid == cid_new: # same coordinate system
        return F, M

    # find the vector r for doing:
    #     M = r x F
    cp_ref = _coord(model, cid)
    coord_to_ref = _coord(model, cid_new)
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


def PositionWRT(xyz, cid, cid_new, model, is_cid_int=None):
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
        is cid/cid_new an integer or a Coord object (deprecated)

    Returns
    -------
    xyz_local : (3, ) float ndarray
        the position of the GRID in an arbitrary coordinate system

    """
    if is_cid_int is not None:  # pragma: no cover
        deprecated('PositionWRT(xyz, cid, cid_new, model, is_cid_int=%s)' % is_cid_int,
                   'PositionWRT(xyz, cid, cid_new, model)', '1.2',
                   levels=[-1])

    if cid == cid_new: # same coordinate system
        return xyz

    cp_ref = _coord(model, cid)
    coord_to_ref = _coord(model, cid_new)

    if 0:  # pragma: no cover
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

        # now converting it to the output coordinate system
        xyz_local = coord_to_ref.transform_node_to_local(xyz_global)

    return xyz_local


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

    Notes
    -----
    xref should be set to False for this function.

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

def _coord(model, cid):
    """helper method"""
    if isinstance(cid, integer_types):
        cp_ref = model.Coord(cid)
    else:
        cp_ref = cid
    return cp_ref
