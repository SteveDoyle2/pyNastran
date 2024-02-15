from typing import Callable
from collections import defaultdict
import numpy as np

from pyNastran.femutils.utils import abs_nan_min_max, safe_nanstd


def nan_difference(x: np.ndarray, axis: int) -> np.ndarray:
    diff = np.nanmax(x, axis=axis) - np.nanmin(x, axis=axis)
    return diff


derivation_map: dict[str, Callable[[np.ndarray, int], np.ndarray]]= {
    'Absolute Max': abs_nan_min_max,
    'Mean': np.nanmean,
    'Max': np.nanmax,
    'Min': np.nanmin,
    'Difference': nan_difference,
    'Std. Dev.': safe_nanstd,
}

def abs_max_scalar(x: np.ndarray) -> np.ndarray:
    mini = np.nanmin(x)
    maxi = np.nanmax(x)
    if np.abs(mini) > np.abs(maxi):
        return mini
    return maxi

def difference_scalar(x: np.ndarray) -> np.ndarray:
    out = np.nanmax(x) - np.nanmin(x)
    return out
def safe_nanstdi(x: np.ndarray) -> np.ndarray:
    return safe_nanstd(x, axis=None)

nodal_combine_map: dict[str, Callable[[np.ndarray], np.ndarray]]= {
    'Absolute Max': abs_max_scalar,
    'Mean': np.nanmean,
    'Max': np.nanmax,
    'Min': np.nanmin,
    'Difference': difference_scalar,  # nan-subtract
    'Std. Dev.': safe_nanstdi,
}

def nodal_average(nodal_combine_func: Callable[[np.ndarray], np.ndarray],
                  element_node: np.ndarray,
                  data: np.ndarray,
                  nids: np.ndarray,
                  #inid: np.ndarray,
                  nid_to_inid_map: dict[int, int]) -> np.ndarray:
    data_dict = defaultdict(list)
    nnode = len(nids)

    data2 = np.full(nnode, np.nan, dtype=data.dtype)
    for (eid, nid), datai in zip(element_node, data):
        data_dict[nid].append(datai)
    for nid, datasi in data_dict.items():
        collapsed_value = nodal_combine_func(datasi)
        inidi = nid_to_inid_map[nid]
        data2[inidi] = collapsed_value
    return data2
