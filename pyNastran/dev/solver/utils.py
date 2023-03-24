from __future__ import annotations
from typing import Any, TYPE_CHECKING
import numpy as np

if TYPE_CHECKING:  # pragma: no cover
    from pyNastran.bdf.bdf import BDF



DOF_MAP = dict[tuple[int, int], int]

def get_ieids_eids(model: BDF, etype: str, eids_str,
                   idtype: str='int32') -> tuple[int, Any, Any, Any]:
    """helper for the stress/strain/force/displacment recovery"""
    eids = np.array(model._type_to_id_map[etype], dtype=idtype)
    if len(eids) == 0:
        return 0, None, None

    if eids_str == 'ALL':
        neids = len(eids)
        ieids = np.arange(neids, dtype=idtype)
    else:
        ieids = np.searchsorted(eids_str, eids)
        neids = len(ieids)
    return neids, ieids, eids


def lambda1d(dxyz, debug=True):
    """
    ::
      3d  [l,m,n,0,0,0]  2x6
          [0,0,0,l,m,n]
    """
    #xyz1 = model.Node(n1).get_position()
    #xyz2 = model.Node(n2).get_position()
    #v1 = xyz2 - xyz1
    if debug:
        print("v1=%s" % dxyz)
    n = np.linalg.norm(dxyz)
    if n == 0:
        raise ZeroDivisionError(dxyz)

    (l, m, n) = dxyz / n
    Lambda = np.zeros((2, 6), dtype='float64')
    Lambda[0, 0] = Lambda[1, 3] = l
    Lambda[0, 1] = Lambda[1, 4] = m
    Lambda[0, 2] = Lambda[1, 5] = n

    #debug = True
    if debug:
        print("Lambda = \n" + str(Lambda))
    return Lambda
