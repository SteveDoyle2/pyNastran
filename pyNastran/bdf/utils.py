"""
Defines various utilities including:
 - parse_patran_syntax
 - parse_patran_syntax_dict
 - Position
 - PositionWRT
 - transform_load

"""
from __future__ import annotations
from copy import deepcopy
from typing import List, Dict, TYPE_CHECKING
import numpy as np  # type: ignore
from numpy import cross, dot  # type: ignore

from pyNastran.utils import deprecated
from pyNastran.utils.numpy_utils import integer_types
from pyNastran.bdf.patran_utils.colon_syntax import (
    parse_patran_syntax, parse_patran_syntax_dict, parse_patran_syntax_dict_map,
    write_patran_syntax_dict)  # pragma: disable=unused-import
if TYPE_CHECKING:  # pragma: no cover
    from pyNastran.nptyping import NDArray3float
    from pyNastran.bdf.bdf import BDF
    from pyNastran.bdf.cards.coordinate_systems import Coord



def Position(xyz: NDArray3float, cid: int, model: BDF):
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


def TransformLoadWRT(F, M, cid, cid_new, model):
    deprecated('TransformLoadWRT', 'transform_load', '1.3', levels=[0, 1, 2])
    return transform_load(F, M, cid, cid_new, model)

def transform_load(F, M, cid: int, cid_new: int, model: BDF):
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
    unused_Mxyz_global = cross(r, Fxyz_global)
    dMxyz_local_2 = cross(r, Fxyz_local_2)
    Mxyz_local_2 = Mxyz_local_1 + dMxyz_local_2

    # rotate the delta moment into the local frame
    unused_M_local = coord_to_ref.xyz_to_coord(Mxyz_local_2)

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


def get_xyz_cid0_dict(model: BDF,
                      xyz_cid0: Dict[int, NDArray3float]=None) -> Dict[int, NDArray3float]:
    """
    helper method

    Parameters
    ----------
    model : BDF()
        a BDF object
    xyz_cid0 : None / Dict[int] = (3, ) ndarray
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

def split_eids_along_nids(model: BDF, eids: List[int], nids: List[int]) -> None:
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

def _coord(model: BDF, cid: int) -> Coord:
    """helper method"""
    if isinstance(cid, integer_types):
        cp_ref = model.Coord(cid)
    else:
        cp_ref = cid
    return cp_ref
