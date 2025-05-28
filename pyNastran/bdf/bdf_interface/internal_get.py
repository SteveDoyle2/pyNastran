from __future__ import annotations
from typing import Optional, TYPE_CHECKING

if TYPE_CHECKING:  # pragma: no cover
    from pyNastran.bdf.cards.nodes import GRID
    from pyNastran.bdf.cards.coordinate_systems import Coord
    from pyNastran.bdf.cards.aero.aero import (
        CAEROs, PAEROs, AELIST, AEFACT)
    from pyNastran.bdf.cards.bdf_sets import SET1, SET2
    from pyNastran.bdf.cards.bdf_tables import TABLEDs
    from pyNastran.bdf.cards.optimization_nx import GROUP


def node_id(node_ref: GRID, nid: int) -> int:
    if node_ref is not None:
        return node_ref.nid
    return nid


def coord_id_blank(coord_ref: Coord, cid: Optional[int]) -> int:
    assert cid is None or cid >= 0, cid
    if coord_ref is not None:
        return coord_ref.cid
    return cid

def coord_id_negative(coord_ref: Coord, cid: int) -> int:
    assert cid >= -1, cid
    if coord_ref is not None:
        return coord_ref.cid
    return cid

def coord_id(coord_ref: Coord, cid: int) -> int:
    if coord_ref is not None:
        return coord_ref.cid
    return cid


# def property_id(prop_ref, pid: int) -> int:
#     if prop_ref is not None:
#         return prop_ref.pid
#     return pid


def element_id(elem_ref, eid: int) -> int:
    if elem_ref is not None:
        return elem_ref.eid
    return eid


def material_id(mid_ref, mid: int) -> int:
    if mid_ref is None:
        return mid
    return mid_ref.mid



def table_id(table_ref: TABLEDs, tid: int) -> int:
    assert isinstance(tid, int), tid
    if table_ref is None:
        return tid
    return table_ref.tid


def paero_id(prop_ref: PAEROs, pid: int) -> int:
    if prop_ref is not None:
        return prop_ref.pid
    return pid


def caero_id(caero_ref: CAEROs, caero: int) -> int:
    if caero_ref is not None:
        return caero_ref.eid
    return caero


def set_id(setg_ref: SET1 | SET2, setg: int) -> int:
    if setg_ref is not None:
        return setg_ref.sid
    return setg


def aelist_id(aelist_ref: AELIST, aelist: int) -> int:
    if aelist_ref is not None:
        return aelist_ref.sid
    return aelist


def aefact_id(aefact_ref: AEFACT, aefact: int) -> int:
    if aefact_ref is not None:
        return aefact_ref.sid
    return aefact


def set_group(set_group_ref: GROUP | SET1,
              set_group_id: int) -> int:
    if hasattr(set_group_ref, 'group_id'):
        return set_group_ref.group_id
    if hasattr(set_group_ref, 'sid'):
        return set_group_ref.sid
    return set_group_id
