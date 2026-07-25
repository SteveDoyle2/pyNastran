from __future__ import annotations
from typing import Optional, TYPE_CHECKING
import numpy as np
from pyNastran.utils.numpy_utils import integer_types
from pyNastran.bdf.mesh_utils.mass_properties import _get_sym_axis
from pyNastran.dev.bdf_vectorized3.bdf_interface.breakdowns import (
    NO_MASS)

if TYPE_CHECKING:
    from pyNastran.dev.bdf_vectorized3.bdf import BDF


def mass_properties(model: BDF, nsm_id: int = 0,
                    element_id: np.ndarray | list[int] | None = None,
                    reference_point: Optional[np.ndarray] = None,
                    inertia_reference: str = "cg",
                    sym_axis: str = "",
                    scale: Optional[float] = None,
                    include_base_mass: bool = True) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    reference_xyz, coord2, is_cg = _update_reference_point(
        model, reference_point, inertia_reference)

    log = model.log
    element_ids_all = []
    # inertias = []
    # mass = 0.
    mass_cg = np.zeros(3, dtype='float64')
    masses = []
    centroids = []
    element_cards = [card for card in model.element_cards
                     if card.n > 0 and card.type not in NO_MASS]

    if nsm_id == 0:
        element_id_nsm = np.array([], dtype='int32')
        mass_nsm = np.array([], dtype='float64')
        centroid_nsm = np.zeros((0, 3), dtype='float64')
        inertia_nsm = np.zeros((0, 6), dtype='float64')
    else:
        nsm_cards = get_nsm_cards(model, nsm_id)
        for nsm_card in nsm_cards:
            #print(nsm_card.get_stats())
            element_id_nsm, mass_nsm, centroid_nsm, inertia_nsm = nsm_card.inertia(element_id=element_id)
            neidi = len(element_id_nsm)
            assert len(mass_nsm) == neidi, (mass_nsm, neidi)
            assert centroid_nsm.shape == (neidi, 3), f'neidi={neidi}, centroid_nsm.shape={centroid_nsm.shape}'
            assert inertia_nsm.shape == (neidi, 6), f'neidi={neidi}, inertia_nsm.shape={inertia_nsm.shape}'
            element_ids_all.append(element_id_nsm)
            masses.append(mass_nsm)
            centroids.append(centroid_nsm)

    if include_base_mass:
        for card in element_cards:
            if element_id:
                eids_common = np.intersect1d(element_id, card.element_id)
                # print(f'eids_common = {eids_common}')
                if len(eids_common) == 0:
                    continue
                card = card.slice_card_by_id(eids_common)

            element_ids_all.append(card.element_id)
            massi = card.mass()
            masses.append(massi)
            centroidi = card.centroid()
            if np.any(np.isnan(centroidi)):
                log.error(f'{card.type} has nan centroid; centroid={centroidi}')
                raise RuntimeError(f'{card.type} has nan centroid; centroid={centroidi}')
            centroids.append(centroidi)

    if len(masses) == 0:
        element_id = np.array([], dtype='int32')
        mass = np.array([], dtype='float64')
        cg = np.zeros((0, 3), dtype='float64')
        inertia = np.zeros((0, 6), dtype='float64')
        log.error('no elements with mass/inertia')
        return element_id, mass, cg, inertia
        # return element_id_nsm, mass_nsm, centroid_nsm, inertia_nsm

    element_id_out = np.hstack(element_ids_all)
    mass = np.hstack(masses)
    centroid = np.vstack(centroids)

    # Find unique keys and map them to clean 0, 1, 2... indices
    # element_ids_unique, inverse_indices = np.unique(element_id_out, return_inverse=True)

    # Sum using mapped indices
    # sums = np.bincount(inverse_indices, weights=mass)

    # abs_mass = np.abs(mass).sum()
    neids = len(element_id_out)
    # if abs_mass == 0.:
    # assert len(element_id) > 0, element_id
    # cg = np.full(3, np.nan, dtype='float64')
    # inertia = np.full(6, np.nan, dtype='float64')
    # log.error('no elements with mass...inertia is nan')
    # return element_id, abs_mass, cg, inertia
    mass_cg = mass[:, None] * centroid
    imass = (mass != 0)
    cg = np.full(centroid.shape, np.nan, dtype=centroid.dtype)
    cg[imass] = mass_cg[imass, :] / mass[imass, np.newaxis]

    # cg = mass_cg.sum(axis=0) / mass.sum()
    # assert len(cg) == 3, cg
    assert cg.shape == (neids, 3), f'neid={neids}, cg.shape={cg.shape}'

    if is_cg:
        # only transform if we're calculating the inertia about the cg
        #     xyz_ref = reference_xyz
        #     xyz_ref2 = cg
        #     inertia = transform_inertia(
        #         mass, cg, xyz_ref, xyz_ref2, inertia, coord1=coord1, coord2=coord2
        #     )
        dxyz = centroid - cg
    else:
        dxyz = centroid - reference_point
    dx = dxyz[:, 0]
    dy = dxyz[:, 1]
    dz = dxyz[:, 2]
    Ixx = mass * (dy ** 2 + dz ** 2)
    Iyy = mass * (dx ** 2 + dz ** 2)
    Izz = mass * (dx ** 2 + dy ** 2)
    Ixy = mass * (dx * dy)
    Ixz = mass * (dx * dz)
    Iyz = mass * (dy * dz)
    inertia = np.stack([Ixx, Iyy, Izz, Ixy, Ixz, Iyz], axis=1, out=None)

    nrows = len(mass)
    assert inertia.shape == (nrows, 6), inertia.shape
    mass.sum()
    mass, cg, inertia = _apply_mass_symmetry(
        model, sym_axis, scale, mass, cg, inertia)
    return element_id_out, mass, centroid, inertia


def get_nsm_cards(model: BDF, nsm_id: int) -> list:
    nsm_ids = model.nsmadd.get_reduced_nsms(stop_on_failure=False)
    log = model.log
    if nsm_id in nsm_ids:
        log.debug('method A')
        log.warning(f'nsm_ids = {nsm_ids}')
        nsm_cards = []
        for cards in nsm_ids[nsm_id]:
            if isinstance(cards, list):
                nsm_cards.extend(cards)
            else:
                raise TypeError(cards)
        # print(f'nsm_cards = {nsm_cards}')
        # raise NotImplementedError(f'nsm_id={nsm_id} not implemented')
    else:
        log.debug('method B')
        nsm_cards = []
        for nsm in [model.nsm, model.nsm1, model.nsml1, model.nsml]:
            if nsm_id not in nsm.nsm_id:
                continue
            nsm_card = nsm.slice_card_by_id(nsm_id)
            nsm_cards.append(nsm_card)
        # log.warning(f'cards = {cards}')
    # assert len(nsm_cards) == 1, nsm_cards
    # assert nsm_cards[0].n == 1, nsm_cards
    return nsm_cards


def _update_reference_point(
    model: BDF, reference_point: np.ndarray, inertia_reference: str = "cg"
) -> tuple[np.ndarray, CORD2R, bool]:
    """helper method for handling reference point"""
    inertia_reference = inertia_reference.lower()
    if inertia_reference == "cg":
        is_cg = True  # nastran-style inertia is always about the cg
    elif inertia_reference == "ref":
        is_cg = False  # inertia is about the reference point
    else:
        raise ValueError("inertia_reference=%r and must be 'cg' or 'ref'" % inertia_reference)

    cid = 0
    if reference_point is None:
        reference_xyz = np.array([0.0, 0.0, 0.0])
    elif isinstance(reference_point, integer_types):
        nid_ref = model.grid.slice_card_by_node_id(reference_point)
        reference_xyz = nid_ref.xyz_cid0().squeeze()
        cid = nid_ref.cd.squeeze()
        coord = model.coord.slice_card_by_id(cid)
        assert coord is not None, nid_ref.get_stats()
        assert coord.type in {"CORD1R", "CORD2R"}, coord
    else:
        # TODO: this method doesn't support coord
        reference_xyz = np.asarray(reference_point, dtype="float64")
        if len(reference_xyz.shape) != 1 or len(reference_xyz) != 3:
            msg = (
                "reference_point=%r and must be None, "
                "a list of 3 floats, or an integer (node id)" % reference_point
            )
            raise ValueError(msg)
    return reference_xyz, cid, is_cg

def _apply_mass_symmetry(
    model: BDF,
    sym_axis: str,
    scale: float,
    mass: float,
    cg: np.ndarray,
    inertia: np.ndarray,
) -> tuple[float, np.ndarray, np.ndarray]:
    """
    Scales the mass & moment of inertia based on the symmetry axes
    and the PARAM WTMASS card

    """
    sym_axis_set = _get_sym_axis(model, sym_axis)
    del sym_axis

    if sym_axis_set:
        # either we figured sym_axis out from the AERO cards or the user told us
        model.log.debug(f"Mass/MOI sym_axis = {sym_axis_set!r}")

        if "xz" in sym_axis_set:
            # y inertias are 0
            cg[1] = 0.0
            mass *= 2.0
            inertia[0] *= 2.0
            inertia[1] *= 2.0
            inertia[2] *= 2.0
            inertia[3] *= 0.0  # Ixy
            inertia[4] *= 2.0  # Ixz; no y
            inertia[5] *= 0.0  # Iyz

        if "xy" in sym_axis_set:
            # z inertias are 0
            cg[2] = 0.0
            mass *= 2.0
            inertia[0] *= 2.0
            inertia[1] *= 2.0
            inertia[2] *= 2.0
            inertia[3] *= 2.0  # Ixy; no z
            inertia[4] *= 0.0  # Ixz
            inertia[5] *= 0.0  # Iyz

        if "yz" in sym_axis_set:
            # x inertias are 0
            cg[0] = 0.0
            mass *= 2.0
            inertia[0] *= 2.0
            inertia[1] *= 2.0
            inertia[2] *= 2.0
            inertia[3] *= 0.0  # Ixy
            inertia[4] *= 0.0  # Ixz
            inertia[5] *= 2.0  # Iyz; no x

    wtmass = model.wtmass
    if scale is None:
        scale = wtmass
        if scale != 1.0:
            model.log.info(f"WTMASS scale = {scale!r}")
    mass *= scale
    inertia *= scale
    return mass, cg, inertia
