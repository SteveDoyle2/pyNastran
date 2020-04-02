"""
Defines various utilities including:
 - parse_patran_syntax
 - parse_patran_syntax_dict

"""
from typing import List, Dict, Optional
import numpy as np  # type: ignore

from pyNastran.bdf.cards.collpase_card import collapse_colon_packs
from pyNastran.utils.numpy_utils import integer_types


def parse_patran_syntax(node_sets: str, pound: Optional[int]=None) -> np.ndarray:
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
    assert isinstance(node_sets, str), type(node_sets)
    if pound is not None:
        assert isinstance(pound, (str, integer_types)), type(pound)
        node_sets = node_sets.replace('#', str(pound).strip())
    if len(node_sets) == 0:
        return np.array([], dtype='int32')

    snodes = node_sets.split()
    nodes = []  # type: List[int]
    for snode in snodes:
        _apply_comma_colon_int_node(nodes, snode)
    return np.unique(nodes)

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

def _apply_colon_set(snode: str) -> List[int]:
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

def write_patran_syntax_dict(dict_sets: Dict[str, np.ndarray]) -> str:
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


def parse_patran_syntax_dict(node_sets: str, pound_dict: Dict[str, Optional[int]]=None,
                             msg: str='') -> Dict[str, np.ndarray]:
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
                msg += 'new_set=%s' % np.array(new_set, dtype='int32')
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
        data[key] = np.unique(ints)
    return data


def parse_patran_syntax_dict_map(node_sets: str,
                                 type_map: Dict[str, str],
                                 msg: str='') -> Dict[str, np.ndarray]:
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
