from __future__ import annotations
import numpy as np
integer_types = (int, )
from typing import TYPE_CHECKING
if TYPE_CHECKING:  # pragma: no cover
    from pyNastran.converters.avl.surface import Surface
    from pyNastran.converters.avl.body import Body

def save_wing_elements(isurface: int,
                       point: np.ndarray, element: np.ndarray,
                       xyz_scale: np.ndarray,
                       dxyz: np.ndarray,
                       nodes: list[np.ndarray], quad_elements: list[np.ndarray],
                       surfaces: list[Union[Surface, Body]],
                       is_cs_list, is_cs,
                       ipoint: int) -> int:
    npoint = point.shape[0]
    nelement = element.shape[0]

    # scaling is applied before translation
    point[:, 0] *= xyz_scale[0]
    point[:, 1] *= xyz_scale[1]
    point[:, 2] *= xyz_scale[2]
    #print(point)
    point += dxyz

    surfaces.append(np.ones(nelement, dtype='int32') * isurface)
    assert len(is_cs) == nelement
    is_cs_list.append(is_cs)

    nodes.append(point)
    quad_elements.append(ipoint + element)

    ipoint += npoint
    return ipoint

def get_spacing(nchord: int, k: float) -> np.ndarray:
    """
    Parameters
    ----------
    nchord : int
        number of chordwise points
    k : spacing
        factor

    Returns
    -------
    x' : (n,) ndarray
        ranges from [0, 1]
    """
    xeq = np.linspace(0., 1., num=nchord+1,
                      endpoint=True, retstep=False, dtype=None)
    if k in {-3.0, 0.0, 3.0}:
        return xeq

    theta = xeq * np.pi
    xcos = 1. / 2 * (1 - np.cos(theta))
    xsine = 1. / 2 * np.sin(theta/2.)
    keq = 0.
    kcos = 0.
    ksin = 0.
    kabs = abs(k)
    if 0.0 <= kabs <= 1.0:
        keq = 1.0 - kabs
        kcos = kabs
    elif kabs <= 2.0:
        kcos = kabs - 1.0
        ksin = 2.0 - kabs
    elif kabs <= 3.0:
        ksin = kabs - 2.0
        keq = 3.0 - kabs
    else:
        raise RuntimeError(f'kabs={kabs} and must be <= 3.0')
    xprime = keq * xeq + ksin * xsine + kcos * xcos
    return xprime

def save_wing_elements(isurface: int,
                       point: np.ndarray, element: np.ndarray,
                       xyz_scale: np.ndarray,
                       dxyz: np.ndarray,
                       nodes: list[np.ndarray], quad_elements: list[np.ndarray],
                       surfaces: list[Union[Surface, Body]],
                       is_cs_list, is_cs,
                       ipoint: int) -> int:
    npoint = point.shape[0]
    nelement = element.shape[0]

    # scaling is applied before translation
    point[:, 0] *= xyz_scale[0]
    point[:, 1] *= xyz_scale[1]
    point[:, 2] *= xyz_scale[2]
    #print(point)
    point += dxyz

    surfaces.append(np.ones(nelement, dtype='int32') * isurface)
    assert len(is_cs) == nelement
    is_cs_list.append(is_cs)

    nodes.append(point)
    quad_elements.append(ipoint + element)

    ipoint += npoint
    return ipoint
