"""
Cut a BDF shell model at spanwise stations and generate ZAERO
BODY7 / PBODY7 / SEGMESH / AEFACT cards from the concave hull
of each cross-section.
"""
from __future__ import annotations
import copy
from pathlib import Path
from itertools import count
from typing import Any, Optional, TYPE_CHECKING

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

from cpylog import SimpleLogger
from pyNastran.utils import PathLike
from pyNastran.bdf.bdf import BDF, read_bdf
from pyNastran.bdf.field_writer_8 import print_card_8
from pyNastran.bdf.cards.coordinate_systems import CORD2R

from pyNastran.utils.concave_hull import (
    concave_hull_2d, order_hull_by_angle,)
from pyNastran.bdf.mesh_utils.cut.cut_model_by_plane import (
    cut_face_model_by_coord, _setup_faces,)

#if TYPE_CHECKING:
#    pass

Rods = tuple[np.ndarray, np.ndarray, np.ndarray]


def get_cut_points(rods: Rods) -> np.ndarray:
    """Extract the unique 3D points from the cut rods.

    Parameters
    ----------
    rods : Rods
        (rod_eid_nodes, rod_nids, rod_xyzs) from cut_face_model_by_coord

    Returns
    -------
    points : (npoints, 3) float ndarray
        unique xyz coordinates of the cut
    """
    rod_eid_nodes, rod_nids, rod_xyzs = rods
    n1 = rod_eid_nodes[:, 1]
    n2 = rod_eid_nodes[:, 2]
    inid1 = np.searchsorted(rod_nids, n1)
    inid2 = np.searchsorted(rod_nids, n2)
    xyz1 = rod_xyzs[inid1, :]
    xyz2 = rod_xyzs[inid2, :]
    all_xyz = np.vstack([xyz1, xyz2])
    _, idx = np.unique(np.round(all_xyz, decimals=10), axis=0, return_index=True)
    return all_xyz[np.sort(idx)]


def cut_and_generate_body7(
    bdf_filename: PathLike | BDF,
    normal_plane: np.ndarray,
    log: SimpleLogger,
    stations: list[float] | np.ndarray,
    coords: list[CORD2R],
    body_id: int = 1,
    label: str = "BODY",
    acoord: int = 0,
    nradial: int = 10,
    alpha: float = 0.0,
    segmesh_id_start: int = 100,
    aefact_id_start: int = 1000,
    pbody7_id: int = 0,
    include_lines: bool = False,
    include_solids: bool = False,
    face_data: Optional[Any] = None,
    plane_atol: float = 1e-5,
    debug_vectorize: bool = True,
    debug_v3: bool = False,
    stop_on_failure: bool = False,
    zaero_model: BDF = None,
    output_filename: PathLike = "",) -> dict[str, Any]:
    """
    Cut a shell model at stations and produce
    ZAERO BODY7/SEGMESH/AEFACT cards.

    Parameters
    ----------
    bdf_filename : PathLike or BDF
        the BDF model (path or object)
    normal_plane : (3,) float ndarray
        the plane normal for the cutting tool
    log : SimpleLogger
        logger
    stations : list[float] or ndarray
        y-locations to cut at (in the coord system of each coord)
    coords : list[CORD2R]
        one coordinate system per station defining the cut plane
    body_id : int
        BODY7 element id
    label : str
        BODY7 label (up to 8 chars)
    acoord : int
        ACOORD id for body centerline orientation (0 = basic)
    nradial : int
        number of circumferential points for each SEGMESH station
    alpha : float
        concave hull alpha parameter (0 = convex hull)
    segmesh_id_start : int
        starting SEGMESH id
    aefact_id_start : int
        starting AEFACT id
    pbody7_id : int
        PBODY7 id (0 = no PBODY7)
    include_lines : bool
        include line elements in face setup
    include_solids : bool
        include solid elements in face setup
    face_data : optional
        pre-computed face data from _setup_faces
    plane_atol : float
        plane intersection tolerance
    debug_vectorize : bool
        use vectorized cutting
    debug_v3 : bool
        use v3 cutting
    stop_on_failure : bool
        raise on failed cuts
    output_filename : PathLike
        if provided, write the ZAERO cards to this file

    Returns
    -------
    result : dict
        stations_found : (n,) float ndarray
            stations where cuts succeeded
        hulls : list[ndarray]
            the concave hull (y, z) points per station
        centroids : (n, 3) float ndarray
            centroid of each cut
        model: BDF
            the zaero cards
    """
    if isinstance(bdf_filename, (str, Path)):
        skip_cards = [
            'DMIG', 'DMIJ', 'DMIK', 'CBAR', 'CBEAM',
            "FORCE", "FORCE1", "FORCE2",
            "MOMENT", "MOMENT1", "MOMENT2",
            "PLOAD", "PLOAD1", "PLOAD2", "PLOAD4",
            'CELAS1', 'CELAS2', 'CELAS3', 'CELAS4',
            "CDAMP1", "CDAMP2", "CDAMP3", "CDAMP4",
            "CBUSH", "CGAP",
            "RBE1", "RBE2", "RBE3", "RBAR", "RBAR1", "RSPLINE",
            "CONM2", "CONM1",
            "CMASS1", "CMASS2", "CMASS3", "CMASS4",
        ]
        model = read_bdf(bdf_filename, log=log, skip_cards=skip_cards)
        model_static = copy.deepcopy(model)
    else:
        model = bdf_filename
        model_static = bdf_filename

    if zaero_model is None:
        zaero_model = BDF()
    zaero_model.set_as_zaero()
    zaero_obj = zaero_model.zaero

    if face_data is None:
        _log, *face_data = _setup_faces(
            model,
            include_lines=include_lines,
            include_solids=include_solids,)

    assert len(stations) == len(coords), (len(stations), len(coords))
    assert len(stations) > 0, "no stations provided"

    cut_results: list[dict[str, Any]] = []

    for icut, dy, coord in zip(count(), stations, coords):
        model.coords[coord.cid] = coord
        nodal_result = None
        try:
            out = cut_face_model_by_coord(
                model_static,
                coord,
                nodal_result,
                plane_atol=plane_atol,
                skip_cleanup=True,
                csv_filename="",
                plane_bdf_filename1="",
                plane_bdf_filename2="",
                plane_y_offset=dy,
                face_data=face_data,
                debug_vectorize=debug_vectorize,
                debug_v3=debug_v3,
                stop_on_failure=stop_on_failure,
            )
        except (RuntimeError, PermissionError):
            log.warning(f"cut failed at station {dy}")
            continue

        found_cut, _geom, _results, rods = out
        if not found_cut:
            log.debug(f"no cut found at station {dy}")
            continue

        points_3d = get_cut_points(rods)
        if len(points_3d) < 3:
            log.warning(f"fewer than 3 points at station {dy}, skipping")
            continue

        centroid_3d = points_3d.mean(axis=0)

        result = {
            "station": dy,
            "points_3d": points_3d,
            "centroid_3d": centroid_3d,
            "icut": icut,
        }
        cut_results.append(result)

    if not cut_results:
        raise RuntimeError("no valid cuts found at any station")

    log.info(f"{len(cut_results)} cuts succeeded out of {len(stations)} stations")

    stations_found = np.array([r["station"] for r in cut_results])
    centroids = np.array([r["centroid_3d"] for r in cut_results])

    hulls = []
    hull_yz_list: list[np.ndarray] = []
    cut_points_2d_list: list[np.ndarray] = []
    cambers: list[float] = []

    for result in cut_results:
        pts = result["points_3d"]
        # cut points are in local coords: col 0 = x (in-plane),
        # col 1 = y ≈ 0 (cut normal), col 2 = z (in-plane)
        pts_xz = pts[:, [0, 2]]
        cut_points_2d_list.append(pts_xz)

        hull = concave_hull_2d(pts_xz, alpha=alpha)
        hull = order_hull_by_angle(hull)
        hulls.append(hull)

        hull_open = hull[:-1]
        y_center = hull_open[:, 0].mean()
        z_center = hull_open[:, 1].mean()
        cambers.append(z_center)

        npts = len(hull_open)
        if npts != nradial:
            angles = np.linspace(0, 2 * np.pi, nradial, endpoint=False)
            hull_centered = hull_open - np.array([y_center, z_center])
            hull_angles = np.arctan2(hull_centered[:, 1], hull_centered[:, 0])
            sort_idx = np.argsort(hull_angles)
            hull_angles_sorted = hull_angles[sort_idx]
            hull_centered_sorted = hull_centered[sort_idx]

            hull_angles_ext = np.concatenate([
                hull_angles_sorted - 2 * np.pi,
                hull_angles_sorted,
                hull_angles_sorted + 2 * np.pi,
            ])
            hull_y_ext = np.concatenate([
                hull_centered_sorted[:, 0],
                hull_centered_sorted[:, 0],
                hull_centered_sorted[:, 0],
            ])
            hull_z_ext = np.concatenate([
                hull_centered_sorted[:, 1],
                hull_centered_sorted[:, 1],
                hull_centered_sorted[:, 1],
            ])

            interp_y = np.interp(angles, hull_angles_ext, hull_y_ext)
            interp_z = np.interp(angles, hull_angles_ext, hull_z_ext)

            resampled_y = interp_y + y_center
            resampled_z = interp_z + z_center
        else:
            resampled_y = hull_open[:, 0]
            resampled_z = hull_open[:, 1]

        hull_yz_list.append(np.column_stack([resampled_y, resampled_z]))

    nstations = len(cut_results)
    nseg = 1
    segmesh_id = segmesh_id_start
    aefact_id = aefact_id_start

    aefact_idy_ids: list[int] = []
    aefact_idz_ids: list[int] = []
    for i, yz_pts in enumerate(hull_yz_list):
        idy = aefact_id
        idz = aefact_id + 1

        y_vals = yz_pts[:, 0].tolist()
        z_vals = yz_pts[:, 1].tolist()

        zaero_model.add_aefact(idy, y_vals)
        zaero_model.add_aefact(idz, z_vals)

        aefact_idy_ids.append(idy)
        aefact_idz_ids.append(idz)
        aefact_id += 2

    itypes: list[int] = []
    x_list: list[float] = []
    cam_list: list[float] = []
    yr_list: list[float] = []
    zr_list: list[float] = []
    idy_list: list[int | None] = []
    idz_list: list[int | None] = []

    for i in range(nstations):
        itypes.append(3)
        x_list.append(float(stations_found[i]))
        cam_list.append(float(cambers[i]))
        yr_list.append(0.0)
        zr_list.append(0.0)
        idy_list.append(aefact_idy_ids[i])
        idz_list.append(aefact_idz_ids[i])

    nose_radius = None
    iaxis = 1
    zaero_obj.add_segmesh(
        segmesh_id, nstations, nradial, nose_radius, iaxis, itypes,
        x_list, cam_list, yr_list, zr_list, idy_list, idz_list)

    idmeshes = [segmesh_id]
    # eid: int, label: str, pid: int,
    # nseg: int, idmeshes: list[int], acoord: int = 0, comment

    comment = " ZAERO BODY7 cards generated by pyNastran cut_body7"
    body7 = zaero_model.zaero.add_body7(
        body_id, label, pbody7_id,
        nseg, idmeshes, acoord=acoord, comment=comment)

    if pbody7_id:
        zaero_model.zaero.add_pbody7(
            pbody7_id, 1, [-0.2, 1.3, 1.1, 0.0, 0.0, 0])

    if output_filename:
        zaero_model.write_bdf(output_filename)
        log.info(f"wrote ZAERO body cards to {str(output_filename)}")

    zaero_model.cross_reference()
    if nstations > 1:
        points, elems = body7.get_points_elements_3d()

    out = {
        "stations_found": stations_found,
        "hulls": hulls,
        "centroids": centroids,
        "hull_yz_resampled": hull_yz_list,
        "cut_points_2d": cut_points_2d_list,
        'model': zaero_model,
    }
    return out
