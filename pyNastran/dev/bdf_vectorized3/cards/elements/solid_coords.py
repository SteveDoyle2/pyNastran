"""Vectorized element and material coordinate systems for solid elements.

Reference: Nastran Quick Reference Guide (CHEXA, CPENTA, CTETRA element descriptions)
"""
from __future__ import annotations
from typing import TYPE_CHECKING
import numpy as np

if TYPE_CHECKING:
    from pyNastran.dev.bdf_vectorized3.bdf import BDF


def chexa_element_coordinate_system(
        n1: np.ndarray, n2: np.ndarray, n3: np.ndarray, n4: np.ndarray,
        n5: np.ndarray, n6: np.ndarray, n7: np.ndarray, n8: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Compute the element coordinate system for CHEXA elements.

    Returns (centroid, xe, ye, ze) where each is (nelements, 3).
    xe is along the 1-2/4-3/5-6/8-7 face direction.
    """
    centroid = (n1 + n2 + n3 + n4 + n5 + n6 + n7 + n8) / 8.
    xe = ((n2 + n3 + n6 + n7) - (n1 + n4 + n8 + n5)) / 4.
    xe_norm = np.linalg.norm(xe, axis=1, keepdims=True)
    xe /= xe_norm

    v = ((n3 + n7 + n8 + n4) - (n1 + n2 + n6 + n5)) / 4.
    ze = np.cross(xe, v, axis=1)
    ze_norm = np.linalg.norm(ze, axis=1, keepdims=True)
    ze /= ze_norm

    ye = np.cross(ze, xe, axis=1)
    ye_norm = np.linalg.norm(ye, axis=1, keepdims=True)
    ye /= ye_norm
    return centroid, xe, ye, ze


def cpenta_element_coordinate_system(
        n1: np.ndarray, n2: np.ndarray, n3: np.ndarray,
        n4: np.ndarray, n5: np.ndarray, n6: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Compute the element coordinate system for CPENTA elements.

    Returns (centroid, xe, ye, ze) where each is (nelements, 3).
    """
    centroid = (n1 + n2 + n3 + n4 + n5 + n6) / 6.
    origin = (n1 + n4) / 2.
    xe = (n2 + n3 + n5 + n6) / 4. - origin
    xe_norm = np.linalg.norm(xe, axis=1, keepdims=True)
    xe /= xe_norm

    v = ((n1 + n3 + n4 + n6) - (n1 + n2 + n4 + n5)) / 4.
    ze = np.cross(xe, v, axis=1)
    ze_norm = np.linalg.norm(ze, axis=1, keepdims=True)
    ze /= ze_norm

    ye = np.cross(ze, xe, axis=1)
    ye_norm = np.linalg.norm(ye, axis=1, keepdims=True)
    ye /= ye_norm
    return centroid, xe, ye, ze


def ctetra_element_coordinate_system(
        n1: np.ndarray, n2: np.ndarray,
        n3: np.ndarray, n4: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Compute the element coordinate system for CTETRA elements.

    Returns (centroid, xe, ye, ze) where each is (nelements, 3).
    """
    centroid = (n1 + n2 + n3 + n4) / 4.
    xe = (n2 + n3 + n4) / 3. - n1
    xe_norm = np.linalg.norm(xe, axis=1, keepdims=True)
    xe /= xe_norm

    v = ((n1 + n3 + n4) - (n1 + n2 + n4)) / 3.
    ze = np.cross(xe, v, axis=1)
    ze_norm = np.linalg.norm(ze, axis=1, keepdims=True)
    ze /= ze_norm

    ye = np.cross(ze, xe, axis=1)
    ye_norm = np.linalg.norm(ye, axis=1, keepdims=True)
    ye /= ye_norm
    return centroid, xe, ye, ze


def cpyram_element_coordinate_system(
        n1: np.ndarray, n2: np.ndarray, n3: np.ndarray,
        n4: np.ndarray, n5: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Compute the element coordinate system for CPYRAM elements.

    Returns (centroid, xe, ye, ze) where each is (nelements, 3).
    xe is along the 1-2/4-3 direction, ze from base toward apex.
    """
    centroid = (n1 + n2 + n3 + n4 + n5) / 5.
    xe = ((n2 + n3) - (n1 + n4)) / 2.
    xe_norm = np.linalg.norm(xe, axis=1, keepdims=True)
    xe /= xe_norm

    base_centroid = (n1 + n2 + n3 + n4) / 4.
    v = n5 - base_centroid
    ze = np.cross(xe, np.cross(v, xe, axis=1), axis=1)
    ze_norm = np.linalg.norm(ze, axis=1, keepdims=True)
    ze /= ze_norm

    ye = np.cross(ze, xe, axis=1)
    ye_norm = np.linalg.norm(ye, axis=1, keepdims=True)
    ye /= ye_norm
    return centroid, xe, ye, ze


def solid_material_coordinate_system(
        model: BDF,
        element_id: np.ndarray,
        property_id: np.ndarray,
        centroid: np.ndarray,
        xe: np.ndarray,
        ye: np.ndarray,
        ze: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Transform element coordinate system to material coordinate system
    based on PSOLID CORDM field.

    CORDM values:
        -1: use basic coordinate system (default for blank)
         0: use basic coordinate system
        >0: use referenced CORDx coordinate system

    Returns (xm, ym, zm) in global coordinates, each (nelements, 3).
    When CORDM specifies a coordinate system, the material axes are the
    axes of that system evaluated at the element centroid.
    When CORDM is 0 or -1, the element coordinate system IS the material system.
    """
    nelements = len(element_id)
    xm = xe.copy()
    ym = ye.copy()
    zm = ze.copy()

    # look up PSOLID coord_id for each element
    psolid = model.psolid
    if psolid.n == 0:
        return xm, ym, zm

    iprop = np.searchsorted(psolid.property_id, property_id)
    iprop = np.clip(iprop, 0, psolid.n - 1)
    valid = (psolid.property_id[iprop] == property_id)

    if not np.any(valid):
        return xm, ym, zm

    cordm = psolid.coord_id[iprop]

    # CORDM <= 0 means use element coordinate system (already set)
    has_cordm = valid & (cordm > 0)
    if not np.any(has_cordm):
        return xm, ym, zm

    # For elements with CORDM > 0, get the coordinate system axes
    coord = model.coord
    unique_cids = np.unique(cordm[has_cordm])
    for cid in unique_cids:
        imask = has_cordm & (cordm == cid)
        if not np.any(imask):
            continue
        # get the transformation matrix for this coord system
        T = coord.xyz_to_global_transform[cid]
        xm[imask] = T[0, :]
        ym[imask] = T[1, :]
        zm[imask] = T[2, :]

    return xm, ym, zm
