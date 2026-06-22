import os
import copy
from datetime import date
from itertools import count
from typing import TextIO, Any

import numpy as np
import scipy as sp
from scipy.sparse import csc_matrix, lil_matrix, dok_matrix, issparse
from scipy.sparse.linalg import ArpackNoConvergence

from cpylog import SimpleLogger

from pyNastran.utils.solver_utils import get_solver
import pyNastran
from pyNastran.utils.numpy_utils import integer_types  # , float_types
from pyNastran.dev.bdf_vectorized3.bdf import BDF, Subcase

from pyNastran.f06.errors import FatalError

from pyNastran.dev.bdf_vectorized3.solver.elements.beam import (
    beam_pg_distributed, beam_pg_point,
    consistent_mass, lumped_mass,
    beam_transforms,)
#-----------------------------------------------


def apply_pload1(model, load, scale, Fb, dof_map, log):
    """Apply PLOAD1 (distributed/concentrated loads on CBAR/CBEAM) to Fb."""

    load_eids = load.element_id
    eids_bar  = np.intersect1d(model.cbar.element_id, load_eids)
    eids_beam = np.intersect1d(model.cbeam.element_id, load_eids)
    all_eids = np.union1d(eids_bar, eids_beam)
    missing_eids = np.setdiff1d(load_eids, all_eids)
    if len(missing_eids):
        log.warning(f"  PLOAD1: element {missing_eids} not found")

    ibar = np.array([i for i, eid in enumerate(load.element_id)
                     if eid in eids_bar])
    ibeam = np.array([i for i, eid in enumerate(load.element_id)
                      if eid in eids_beam])

    if len(ibar):
        i = ibar
        eid = load.element_id[i]
        elem = model.cbar.slice_card_by_id(eid)
        set_Fb(model, elem, load, ibar, Fb, dof_map, scale=scale)

    if len(ibeam):
        i = ibeam
        eid = load.element_id[i]
        elem = model.cbeam.slice_card_by_id(eid)
        set_Fb(model, elem, load, ibeam, Fb, dof_map, scale=scale)
    return

#   for i_load in range(load.n):
#       eid = load.element_id[i_load]
#       load_type = str(load.load_type[i_load]).strip()
#       scale_type = str(load.scale[i_load]).strip()
#       x1, x2 = load.x[i_load]
#       p1, p2 = load.pressure[i_load]
#
#       # Find the element (CBAR or CBEAM)
#       # TODO: split this...
#       elem = None
#       for e in [model.cbar, model.cbeam]:
#           if e.n == 0:
#               continue
#           idx = np.where(e.element_id == eid)[0]
#           if len(idx) > 0:
#               elem = e.slice_card_by_index(idx)
#               break
#       if elem is None:
#           log.warning(f"  PLOAD1: element {eid} not found")
#           continue
#
#       # Get element properties
#       xyz1, xyz2 = elem.get_xyz()
#       xyz1i = xyz1[0]
#       xyz2i = xyz2[0]
#       L = np.linalg.norm(xyz2i - xyz1i)
#       v, ihat, yhat, zhat, wa, wb = elem.get_axes(xyz1, xyz2)
#       Teb = beam_transform(ihat[0], yhat[0], zhat[0])
#
#       e_g_nus = elem.e_g_nu()
#       E_val, G_val, nu_val = e_g_nus[0]
#       inertia = elem.inertia()
#       I1, I2, I12, J_val = inertia[0]
#       area = elem.area()[0]
#       k_arr = elem.k()
#       k1, k2 = k_arr[0]
#
#       # Convert scale type to actual x positions
#       if scale_type in ("FR", "FRPR"):
#           xa = x1 * L
#           xb = x2 * L
#       else:
#           xa = x1
#           xb = x2
#
#       # Determine load direction in element local coords
#       # FX/FY/FZ = basic coord; FXE/FYE/FZE = element coord
#       is_element = load_type.endswith("E")
#       base_type = load_type.rstrip("E")
#
#       qx = qy = qz = 0.0
#       qx_end = qy_end = qz_end = 0.0
#       is_force = base_type.startswith("F")
#
#       if is_force:
#           direction = base_type[1]  # X, Y, or Z
#           if is_element:
#               if direction == "X":
#                   qx, qx_end = p1, p2
#               elif direction == "Y":
#                   qy, qy_end = p1, p2
#               else:
#                   qz, qz_end = p1, p2
#           else:
#               # Basic coord -> element local via rotation
#               T = np.vstack([ihat[0], yhat[0], zhat[0]])
#               if direction == "X":
#                   q_basic = np.array([1.0, 0.0, 0.0])
#               elif direction == "Y":
#                   q_basic = np.array([0.0, 1.0, 0.0])
#               else:
#                   q_basic = np.array([0.0, 0.0, 1.0])
#               assert T.shape == (3, 3), T.shape
#               q_local = T @ q_basic
#               qx, qy, qz = p1 * q_local
#               qx_end, qy_end, qz_end = p2 * q_local
#       else:
#           # Moment loads (MX, MY, MZ) — not yet supported
#           log.warning(f"  PLOAD1: moment load type {load_type} not yet supported")
#           continue
#
#       # Point load vs distributed
#       is_point = abs(xa - xb) < 1e-14 * L
#       if is_point:
#           fe = beam_pg_point(
#               qx * scale, qy * scale, qz * scale,
#               xa, L, E_val, G_val, area, I1, I2, k1, k2,
#           )
#       else:
#           fe = beam_pg_distributed(
#               qx * scale,
#               qy * scale,
#               qz * scale,
#               L, E_val, G_val, area, I1, I2, k1, k2,
#               x_start=xa, x_end=xb,
#               qx_end=qx_end * scale,
#               qy_end=qy_end * scale,
#               qz_end=qz_end * scale,
#           )
#
#       # Transform to basic and assemble
#       PG = Teb.T @ fe
#
#       nid1, nid2 = elem.nodes[0]
#       gi1 = dof_map[(nid1, 1)]
#       gi2 = dof_map[(nid2, 1)]
#       Fb[gi1:gi1 + 6] += PG[:6]
#       Fb[gi2:gi2 + 6] += PG[6:]


def coord_transforms(ihat, jhat, khat):
    neid = len(ihat)
    
    T2 = np.full((neid, 3, 3), np.nan, dtype='float32')
    T2[:, 0, :] = ihat
    T2[:, 1, :] = jhat
    T2[:, 2, :] = khat
    #T3 = np.vstack([ihat, jhat, normal])

    T = np.full((neid, 3, 3), np.nan, dtype=ihat.dtype)
    for i, ihati, jhati, khati in zip(count(), ihat, jhat, khat):
        Ti = np.vstack([ihati, jhati, khati])
        T[i, :, :] = Ti
    assert np.allclose(T, T2)
    return T


def set_Fb(model: BDF,
           elem, load, ibar,
           Fb, dof_map, scale: float=1.0):
    xyz1, xyz2 = elem.get_xyz()
    LAIJEG = elem.stiffness_info()
    # columns: [length, area, I1, I2, I12, J, k1, k2, E, G]
    assert LAIJEG.shape[1] == 10, LAIJEG.shape
    L = LAIJEG[:, 0]
    A = LAIJEG[:, 1]
    I = LAIJEG[:, [2, 3, 4]]
    J = LAIJEG[:, 5]
    k1 = LAIJEG[:, 6] 
    k2 = LAIJEG[:, 7]
    if elem.type == 'CBAR':
        E = LAIJEG[:, 8]
        G = LAIJEG[:, 9]
    elif elem.type == 'CBEAM':
        s1 = LAIJEG[:, 8]
        s2 = LAIJEG[:, 9]
        E = LAIJEG[:, 10]
        G = LAIJEG[:, 11]

    # Get element properties
    xyz1, xyz2 = elem.get_xyz()
    iload = np.arange(load.n)[ibar]

    v, ihat, yhat, zhat, wa, wb = elem.get_axes(xyz1, xyz2)

    T = coord_transforms(ihat, yhat, zhat)
    Teb = beam_transforms(ihat, yhat, zhat)

    for (i, (nid1, nid2), Li, Ai, (i1i, i2i, i12ii), Ji,
            k1i, k2i, Tebi, Ti, Ei, Gi) in zip(
            iload, elem.nodes,
            L, A, I, J, k1, k2, Teb, T, E, G):
        eid = load.element_id[i]
        load_type = load.load_type[i]
        scale_type = load.scale[i]
        x1, x2 = load.x[i]
        p1, p2 = load.pressure[i]

        # Convert scale type to actual x positions
        if scale_type in ("FR", "FRPR"):
            xa = x1 * Li
            xb = x2 * Li
        else:
            xa = x1
            xb = x2

        # Determine load direction in element local coords
        # FX/FY/FZ = basic coord; FXE/FYE/FZE = element coord
        is_element = load_type.endswith("E")
        base_type = load_type.rstrip("E")

        qx = qy = qz = 0.0
        qx_end = qy_end = qz_end = 0.0
        is_force = base_type.startswith("F")

        if is_force:
            direction = base_type[1]  # X, Y, or Z
            if is_element:
                if direction == "X":
                    qx, qx_end = p1, p2
                elif direction == "Y":
                    qy, qy_end = p1, p2
                elif direction == 'Z':
                    qz, qz_end = p1, p2
                else:
                    raise RuntimeError(load_type)
            else:
                # Basic coord -> element local via rotation
                #T = Tebi
                if direction == "X":
                    q_basic = np.array([1.0, 0.0, 0.0])
                elif direction == "Y":
                    q_basic = np.array([0.0, 1.0, 0.0])
                elif direction == 'Z':
                    q_basic = np.array([0.0, 0.0, 1.0])
                else:
                    raise RuntimeError(load_type)
                assert Ti.shape == (3, 3), Ti.shape
                q_local = Ti @ q_basic
                qx, qy, qz = p1 * q_local
                qx_end, qy_end, qz_end = p2 * q_local
        else:
            # Moment loads (MX, MY, MZ) — not yet supported
            NotImplementedError(f"  PLOAD1: moment load type {load_type} not yet supported")

        # Point load vs distributed
        is_point = abs(xa - xb) < 1e-14 * Li
        if is_point:
            fe = beam_pg_point(
                qx * scale, qy * scale, qz * scale,
                xa, Li, Ei, Gi, Ai, i1i, i2i, k1i, k2i,
            )
        else:
            fe = beam_pg_distributed(
                qx * scale, qy * scale, qz * scale,
                Li, Ei, Gi, Ai, i1i, i2i, k1i, k2i,
                x_start=xa, x_end=xb,
                qx_end=qx_end * scale, qy_end=qy_end * scale, qz_end=qz_end * scale,
            )

        # Transform to basic and assemble
        PG = Tebi.T @ fe

        gi1 = dof_map[(nid1, 1)]
        gi2 = dof_map[(nid2, 1)]
        Fb[gi1:gi1 + 6] += PG[:6]
        Fb[gi2:gi2 + 6] += PG[6:]
