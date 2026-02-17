from __future__ import annotations
from typing import Any, TYPE_CHECKING
import numpy as np
if TYPE_CHECKING:
    from pyNastran.dev.bdf_vectorized3.bdf import BDF


def recast_data(idtype: str,
                fdtype: str) -> tuple[str, str]:
    idtype = 'int32'
    fdtype = 'float32'
    return idtype, fdtype

def get_ieids_eids(model: BDF, etype: str, eids_str: str,
                   idtype: str='int32') -> tuple[int, Any, Any, Any]:
    """helper for the stress/strain/force/displacment recovery"""
    # eids = np.array(model._type_to_id_map[etype], dtype=idtype)
    elem = getattr(model, etype.lower())
    if len(elem) == 0:
        return 0, None, None

    eids = elem.element_id
    if eids_str == 'ALL':
        neids = len(eids)
        ieids = np.arange(neids, dtype=idtype)
    else:
        ieids = np.searchsorted(eids_str, eids)
        neids = len(ieids)
    return neids, ieids, eids


def get_element(model: BDF, element_name: str,
                ieids: Optional[np.ndarray],
                eids: np.ndarray):
    elem = getattr(model, element_name.lower())
    if ieids is not None:
        elem = elem.slice_card_by_element_id(eids)
    return elem
