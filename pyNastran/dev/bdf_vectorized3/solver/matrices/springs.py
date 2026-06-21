from __future__ import annotations
import copy
from itertools import count
from typing import Any, TYPE_CHECKING

import numpy as np
from scipy.sparse import dok_matrix

from ..utils import DOF_MAP


if TYPE_CHECKING:  # pragma: no cover
    from pyNastran.nptyping_interface import NDArrayNNfloat
    from pyNastran.dev.bdf_vectorized3.bdf import BDF
    from pyNastran.dev.bdf_vectorized3.bdf_interface.bdf_attributes import (
        CELAS1, CELAS2, CELAS3, CELAS4,)


def _build_kbb_celas1(model: BDF, Kbb: dok_matrix,
                      dof_map: DOF_MAP) -> None:
    """fill the CELAS1 Kbb matrix"""
    celas = model.celas1
    pelas = model.pelas
    if celas.n == 0:
        return celas.n
    pelas = pelas.slice_card_by_id(celas.property_id, assume_sorted=True)
    nids1 = celas.nodes[:, 0]
    nids2 = celas.nodes[:, 1]
    c1s = celas.components[:, 0]
    c2s = celas.components[:, 1]
    # pelas = model.celas2
    ks = pelas.k
    for nid1, nid2, c1, c2, ki in zip(nids1, nids2, c1s, c2s, ks):
        i = dof_map[(nid1, c1)]
        j = dof_map[(nid2, c2)]
        k = ki * np.array([
            [1, -1,],
            [-1, 1],
        ])
        ibe = [
            (i, 0),
            (j, 1),
        ]
        for ib1, ie1 in ibe:
            for ib2, ie2 in ibe:
                Kbb[ib1, ib2] += k[ie1, ie2]
    return celas.n


def _build_kbb_celas2(model: BDF, Kbb: dok_matrix, dof_map: DOF_MAP) -> int:
    """fill the CELAS2 Kbb matrix"""
    celas = model.celas2
    if celas.n == 0:
        return celas.n
    pelas = celas
    nids1 = celas.nodes[:, 0]
    nids2 = celas.nodes[:, 1]
    c1s = celas.components[:, 0]
    c2s = celas.components[:, 1]
    pelas = model.celas2
    ks = pelas.k
    k_unscaled = np.array([
        [1, -1,],
        [-1, 1],
    ])
    for nid1, nid2, c1, c2, ki in zip(nids1, nids2, c1s, c2s, ks):
        i = dof_map[(nid1, c1)]
        j = dof_map[(nid2, c2)]
        k = ki * k_unscaled
        ibe = [
            (i, 0),
            (j, 1),
        ]
        for ib1, ie1 in ibe:
            for ib2, ie2 in ibe:
                Kbb[ib1, ib2] += k[ie1, ie2]
    return celas.n


def _build_kbb_celas3(model: BDF, Kbb: dok_matrix, dof_map: DOF_MAP) -> None:
    """fill the CELAS3 Kbb matrix"""
    element = model.celas3
    nelement = len(element)
    if nelement == 0:
        return nelement
    eids = element.element_id
    pelas = model.pelas.slice_card_by_id(element.property_id, assume_sorted=True)
    nids1 = element.spoints[:, 0]
    nids2 = element.spoints[:, 1]

    ks = pelas.k
    Ke = np.full((nelement, 2, 2), np.nan, dtype="float64")
    for ielement, nid1, nid2, ki in zip(count(), nids1, nids2, ks):
        # i = dof_map[(nid1, c1)]
        # j = dof_map[(nid2, c2)]
        # for eid in eids:
        # elem = model.elements[eid]
        # ki = elem.K()
        # print(elem, ki)
        # print(elem.get_stats())
        ke = _build_kbbi_celas34(Kbb, dof_map, nid1, nid2, ki)
        Ke[ielement, :, :] = ke
    return nelement


def _build_kbb_celas4(model: BDF, Kbb: dok_matrix, dof_map: DOF_MAP) -> int:
    """fill the CELAS4 Kbb matrix"""
    element = model.celas4
    nelement = len(element)
    if nelement == 0:
        return nelement

    eids = element.element_id
    ks = element.k
    nids1 = element.spoints[:, 0]
    nids2 = element.spoints[:, 1]

    Ke = np.full((nelement, 2, 2), np.nan, dtype="float64")
    for ielement, eid, nid1, nid2, ki in zip(count(), eids, nids1, nids2, ks):
        ke = _build_kbbi_celas34(Kbb, dof_map, nid1, nid2, ki)
        Ke[ielement, :, :] = ke
    return nelement


def _build_kbbi_celas12(Kbb: dok_matrix, dof_map: DOF_MAP,
                        elem: CELAS1 | CELAS2, ki: float) -> np.ndarray:
    """fill the CELASx Kbb matrix"""
    nid1, nid2 = elem.nodes
    c1, c2 = elem.c1, elem.c2
    i = dof_map[(nid1, c1)]
    j = dof_map[(nid2, c2)]
    ke = ki * np.array([
        [1, -1,],
        [-1, 1],
    ])
    ibe = [
        (i, 0),
        (j, 1),
    ]
    for ib1, ie1 in ibe:
        for ib2, ie2 in ibe:
            Kbb[ib1, ib2] += ke[ie1, ie2]
    # Kbb[j, i] += ki
    # Kbb[i, j] += ki
    # del i, j, ki, nid1, nid2, c1, c2
    return ke


def _build_kbbi_celas34(Kbb: dok_matrix, dof_map: DOF_MAP,
                        nid1: int, nid2: int, ki: float) -> np.ndarray:
    """fill the CELASx Kbb matrix"""
    # print(dof_map)
    i = dof_map[(nid1, 0)]
    j = dof_map[(nid2, 0)]
    ke = ki * np.array([
        [1, -1,],
        [-1, 1],
    ])
    ibe = [
        (i, 0),
        (j, 1),
    ]
    for ib1, ie1 in ibe:
        for ib2, ie2 in ibe:
            Kbb[ib1, ib2] += ke[ie1, ie2]
    # Kbb[j, i] += ki
    # Kbb[i, j] += ki
    # del i, j, ki, nid1, nid2, c1, c2
    return ke

