"""
Defines various utilities including:
 - parse_patran_syntax
 - parse_patran_syntax_dict
 - Position
 - PositionWRT
 - transform_load

"""
from __future__ import annotations
import os
from copy import deepcopy
from typing import TYPE_CHECKING
import numpy as np  # type: ignore

from pyNastran.utils import deprecated
from pyNastran.utils.numpy_utils import integer_types
from pyNastran.bdf.patran_utils.colon_syntax import (
    parse_patran_syntax, parse_patran_syntax_dict, parse_patran_syntax_dict_map,
    write_patran_syntax_dict)  # pragma: disable=unused-import
from pyNastran.utils import print_bad_path, PathLike

if TYPE_CHECKING:  # pragma: no cover
    from pyNastran.nptyping_interface import NDArray3float
    from pyNastran.bdf.bdf import BDF
    from pyNastran.bdf.cards.coordinate_systems import Coord


def split_comment_to_femap_comment(comment: str,
                                   ) -> tuple[bool, str, tuple[str, int, str]]:
    """
    Parameters
    ----------
    comment : str
        examples are as follows:
        ['$ Femap Region 12345 : Wing NSM']
        ['$ Femap Property 100 : Wing Skin 20 Plies\n$\n$ Femap Layup 101 : 20 Ply\n']
        ['$$ Femap Material 202 : Steel:42\n']
        ['$ Femap Property 8000007 : Aileron, Steel Pin dia=.375', '$ Femap PropShape 8000007 : 5,0,0.1875,0.,0.,0.,0.,0.', '$ Femap PropMethod 8000007 : 5,0,1,0.', '$ Femap PropOrient 8000007 : 5,0,0.,1.,2.,3.,4.,-1.,0.,0.']

    Returns
    -------
    is_passed : bool
        was the comment split successfully
    base_comment
        TODO: maybe remove this?
        empty if successful
        otherwise, just comment
    word_id_name : tuple[str, int, str]
        word : str
            'Femap Property'
        id : int
            100
        name : str
            'Wing Skin 20 Plies'
    """
    comment = comment.rstrip()
    if 'Femap' not in comment:
        return False, comment, ('', -1, '')

    lines = comment.split('\n')
    lines = [line.rstrip('$\t\n\r ') for line in lines if line.rstrip('$\t\n\r ') and ':' in line]
    if len(lines) == 1:
        line = lines[0]
        word, idi, name = _femap_comment_to_sline(line)
    else:
        for line in lines:
            out = _femap_comment_to_sline(line)
        return False, comment, ('', -1, '')
    return True, '', (word, idi, name)

    #print(f'femap_word_id = {femap_word_id!r}')

def _femap_comment_to_sline(line: str) -> tuple[str, int, str]:
    # print(line)
    try:
        femap_word_id, name = line.split(':', 1)
    except ValueError:
        print(line)
        raise
    name = name.strip()
    word_id = femap_word_id.split('Femap')[1]
    word, id_str = word_id.strip().split()
    word = word.strip()
    idi = int(id_str)
    return word, idi, name

def parse_femap_syntax_copy(filename_lines: PathLike | list[str],
                            combine_rows: bool=False) -> dict[int, np.ndarray]:
    """Parses the following syntax from FEMAP (use Copy vs. Copy-As-List):


    lines = ['7203615,7203654,1 7203990,7204010,1 7204032,7204050,1',
             '7203594,7203614,1 7203655,7203675,1 7203969,7203989,1 7204011,7204031,1',]

    Returns
    -------
    out : dict[key, values]
        key : int
            the line number
        values : (n,) int np.ndarray
            sorted ids

    .. note:: assume_unique_row avoids combining the "groups"
    """
    if isinstance(filename_lines, list):
        lines = filename_lines
    else:
        assert os.path.exists(filename_lines), print_bad_path(filename_lines)
        with open(filename_lines, 'r') as prop_file:
            lines = prop_file.readlines()

    lines2 = []
    for line in lines:
        line = line.split('#')[0].strip()
        if line:
            lines2.append(line)
    assert len(lines2) > 0, lines2

    row_values = {}
    for i, line in enumerate(lines2):
        # split line and remove blank entries
        sline = [val.strip() for val in line.split(' ') if val.strip()]

        values = []
        row_values[i] = values
        for pair in sline:
            if ',' in pair:
                # 7203615,7203654,1
                sline = pair.split(',')
                assert ' ' not in sline, sline
                if len(sline) == 1:
                    values.append(int(sline[0]))
                elif len(sline) == 3:
                    start, stop, step = [int(val) for val in sline]
                    valuesi = range(start, stop+step, step)
                    values.extend(valuesi)
                else:
                    raise NotImplementedError(f'line={line!r} pair={pair!r} sline={sline}')

    out = {key: np.unique(values) for key, values in row_values.items()}
    if combine_rows:
        values = [values for key, values in out.items()]
        out2 = {
            0: np.unique(np.hstack(values)),
        }
        return out2
    return out

def parse_femap_syntax(lines: list[str]) -> np.ndarray:
    """Parses the following syntax from FEMAP:

    Add            1646           0           1
    Add            1422        1502           1
    Add            1505        1645           1

    Returns
    -------
    out : (n,) int np.ndarray
        sorted ids

    .. note:: A list of lines is expected
    """
    assert isinstance(lines, list), lines
    values = []
    for line in lines:
        line = line.strip()
        if len(line) == 0:
            continue

        assert '\n' not in line, line
        sline = line.split()
        assert len(sline) == 4, sline
        assert sline[0] == 'Add', sline
        word, start, stop, step = sline
        istart = int(start)
        if stop == '0':
            values.append(istart)
        elif step == '1':
            istop = int(stop)
            valuesi = np.arange(istart, istop+1)
            values.append(valuesi)
        else:
            istop = int(stop)
            istep = int(step)
            valuesi = np.arange(istart, istop+1, istep)
            values.append(valuesi)

    values2 = np.unique(np.hstack(values))
    return values2


def get_femap_property_comments_dict(data_dict):
    return _get_femap_comments_dict(data_dict, word='Femap Property')


def get_femap_material_comments_dict(data_dict):
    return _get_femap_comments_dict(data_dict, word='Femap Material')


def _get_femap_comments_dict(data_dict, word: str):
    word = word.lower()
    comment_dict = {}
    for pid, prop in data_dict.items():
        lines = prop.comment.split('\n')
        commenti = ''
        for line in lines:
            if word in line.lower():
                base, commenti = line.split(':', 1)
                break
        comment_dict[pid] = commenti.strip()
    return comment_dict


def Position(xyz: NDArray3float, cid: int, model: BDF) -> np.ndarray:
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
    cp_ref = _coord(model, cid)
    xyz2 = cp_ref.transform_node_to_global(xyz)
    return xyz2


def transform_load(F: np.ndarray, M: np.ndarray,
                   cid: int, cid_new: int, model: BDF) -> tuple[np.ndarray, np.ndarray]:
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

    Returns
    -------
    Fxyz_local : (3, ) float ndarray
        the force in an arbitrary coordinate system
    Mxyz_local : (3, ) float ndarray
        the force in an arbitrary coordinate system

    """
    if cid == cid_new:  # same coordinate system
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

    # rotate force and moment from coord1 to coord2
    beta1 = cp_ref.beta()
    beta2 = coord_to_ref.beta()
    Fxyz_local_2 = (Fxyz_local_1 @ beta1) @ beta2.T
    Mxyz_local_2 = (Mxyz_local_1 @ beta1) @ beta2.T

    # moment contribution from force at new origin: dM = r x F
    # r and F must be in the same frame — use coord2's local frame
    r_local_2 = r @ beta2.T
    dMxyz_local_2 = np.cross(r_local_2, Fxyz_local_2)
    Mxyz_local_2 = Mxyz_local_2 + dMxyz_local_2

    return Fxyz_local_2, Mxyz_local_2


def PositionWRT(xyz: NDArray3float, cid: int, cid_new: int, model: BDF) -> NDArray3float:
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

    Returns
    -------
    xyz_local : (3, ) float ndarray
        the position of the GRID in an arbitrary coordinate system

    """
    if cid == cid_new:  # same coordinate system
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
        p2_local = np.dot(
            (p1_local @ cp_ref.beta()) + cp_ref.origin - coord_to_ref.origin,
            coord_to_ref.beta().T)

        # convert xyz_2 to R-Theta-Z_2
        xyz_local = coord_to_ref.xyz_to_coord(p2_local)
    else:
        # converting the xyz point arbitrary->global
        xyz_global = cp_ref.transform_node_to_global(xyz)

        # now converting it to the output coordinate system
        xyz_local = coord_to_ref.transform_node_to_local(xyz_global)

    return xyz_local


def get_xyz_cid0_dict(model: BDF,
                      xyz_cid0: dict[int, NDArray3float]=None) -> dict[int, NDArray3float]:
    """
    helper method

    Parameters
    ----------
    model : BDF()
        a BDF object
    xyz_cid0 : None / dict[int] = (3, ) ndarray
        the nodes in the global coordinate system

    Returns
    -------
    xyz_cid0_dict
    """
    if xyz_cid0 is None:
        xyz = {}
        for nid, node in model.nodes.items():
            xyz[nid] = node.get_position()
    else:
        xyz = xyz_cid0
    return xyz


def split_eids_along_nids(model: BDF, eids: list[int], nids: list[int]) -> None:
    """
    Disassociate a list of elements along a list of nodes.

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


def _coord(model: BDF, cid: int) -> Coord:
    """helper method"""
    if isinstance(cid, integer_types):
        cp_ref = model.Coord(cid)
    else:
        cp_ref = cid
    return cp_ref
