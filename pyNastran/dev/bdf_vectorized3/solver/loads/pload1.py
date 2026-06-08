import os
import copy
from datetime import date
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
    beam_transform,)
#-----------------------------------------------


def apply_pload1(model, load, scale, Fb, dof_map, log):
    """Apply PLOAD1 (distributed/concentrated loads on CBAR/CBEAM) to Fb."""
    for i_load in range(load.n):
        eid = load.element_id[i_load]
        load_type = str(load.load_type[i_load]).strip()
        scale_type = str(load.scale[i_load]).strip()
        x1, x2 = load.x[i_load]
        p1, p2 = load.pressure[i_load]

        # Find the element (CBAR or CBEAM)
        elem = None
        for e in [model.cbar, model.cbeam]:
            if e.n == 0:
                continue
            idx = np.where(e.element_id == eid)[0]
            if len(idx) > 0:
                elem = e.slice_card_by_index(idx)
                break
        if elem is None:
            log.warning(f"  PLOAD1: element {eid} not found")
            continue

        # Get element properties
        xyz1, xyz2 = elem.get_xyz()
        xyz1i = xyz1[0]
        xyz2i = xyz2[0]
        L = np.linalg.norm(xyz2i - xyz1i)
        v, ihat, yhat, zhat, wa, wb = elem.get_axes(xyz1, xyz2)
        e_g_nus = elem.e_g_nu()
        E_val, G_val, nu_val = e_g_nus[0]
        inertia = elem.inertia()
        I1, I2, I12, J_val = inertia[0]
        area = elem.area()[0]
        k_arr = elem.k()
        k1, k2 = k_arr[0]

        # Convert scale type to actual x positions
        if scale_type in ("FR", "FRPR"):
            xa = x1 * L
            xb = x2 * L
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
                else:
                    qz, qz_end = p1, p2
            else:
                # Basic coord -> element local via rotation
                T = np.vstack([ihat[0], yhat[0], zhat[0]])
                if direction == "X":
                    q_basic = np.array([1.0, 0.0, 0.0])
                elif direction == "Y":
                    q_basic = np.array([0.0, 1.0, 0.0])
                else:
                    q_basic = np.array([0.0, 0.0, 1.0])
                q_local = T @ q_basic
                qx, qy, qz = p1 * q_local
                qx_end, qy_end, qz_end = p2 * q_local
        else:
            # Moment loads (MX, MY, MZ) — not yet supported
            log.warning(f"  PLOAD1: moment load type {load_type} not yet supported")
            continue

        # Point load vs distributed
        is_point = abs(xa - xb) < 1e-14 * L
        if is_point:
            fe = beam_pg_point(
                qx * scale, qy * scale, qz * scale,
                xa, L, E_val, G_val, area, I1, I2, k1, k2,
            )
        else:
            fe = beam_pg_distributed(
                qx * scale, qy * scale, qz * scale,
                L, E_val, G_val, area, I1, I2, k1, k2,
                x_start=xa, x_end=xb,
                qx_end=qx_end * scale, qy_end=qy_end * scale, qz_end=qz_end * scale,
            )

        # Transform to basic and assemble
        Teb = beam_transform(ihat[0], yhat[0], zhat[0])
        PG = Teb.T @ fe

        nid1, nid2 = elem.nodes[0]
        gi1 = dof_map[(nid1, 1)]
        gi2 = dof_map[(nid2, 1)]
        Fb[gi1:gi1 + 6] += PG[:6]
        Fb[gi2:gi2 + 6] += PG[6:]
