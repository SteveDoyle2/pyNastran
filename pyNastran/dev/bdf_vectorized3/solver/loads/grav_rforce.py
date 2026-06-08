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


def apply_grav(model, load, scale, Fb, dof_map, ndof, log):
    """Apply GRAV load (gravity) to Fb using grid-level mass approach.

    Following MYSTRAN's approach: F_i = M_grid_i * acceleration.
    The assembled mass matrix is used to compute forces at each grid.
    """
    for i_load in range(load.n):
        sid = load.load_id[i_load]
        N_vec = load.N[i_load]  # direction vector (3,)
        mag = load.scale_factor[i_load]  # scale factor (acceleration magnitude)
        cid = load.coord_id[i_load]

        # Acceleration vector in basic frame
        norm_N = np.linalg.norm(N_vec)
        if norm_N == 0.0:
            continue
        accel_basic = mag * N_vec / norm_N

        # Transform from CID to basic if needed
        if cid != 0:
            coord = model.coord
            beta = coord.xyz_to_global_transform[cid]
            accel_basic = beta.T @ accel_basic

        accel_basic *= scale

        log.debug(f"  GRAV: accel_basic={accel_basic}")

        # Build lumped mass matrix and apply F = M * a at each grid
        # Use the assembled Mbb to get per-grid mass (6x6 diagonal blocks)
        # More efficient: iterate grids, get diagonal mass, apply accel
        grid = model.grid
        for nid in grid.node_id:
            nid_int = int(nid)
            i1 = dof_map[(nid_int, 1)]
            # Get translational mass from Mbb diagonal (will be assembled later)
            # Instead, apply accel to the force vector; the mass contribution
            # is handled by assembling F = rho*V*a for each element

        # Grid-level approach: we need the assembled mass matrix.
        # Build a temporary Mbb for this calculation.
        # For efficiency, use the same Mbb that will be built for the analysis.
        # The caller (build_Fb) is called AFTER Mbb is built in SOL 101.
        # However, build_Fb is called before Mbb->Mgg transform in the flow.
        # We must build our own mass here.
        temp_subcase = Subcase(id=0)
        Mbb_temp = build_Mbb(model, temp_subcase, dof_map, ndof, fdtype="float64")

        # Apply F = M * a for each grid (translational DOFs only)
        for nid in grid.node_id:
            nid_int = int(nid)
            i1 = dof_map[(nid_int, 1)]
            # Extract the 3x3 translational mass block for this grid
            m_diag = np.array([
                Mbb_temp[i1, i1],
                Mbb_temp[i1 + 1, i1 + 1],
                Mbb_temp[i1 + 2, i1 + 2],
            ])
            Fb[i1: i1 + 3] += m_diag * accel_basic


def apply_rforce(model, load, scale, Fb, dof_map, ndof, log):
    """Apply RFORCE load (centrifugal/rotational force) to Fb.

    Following MYSTRAN's approach:
    F_i = M_i * [omega x (omega x r_i) + alpha x r_i]
    where r_i = position of grid i relative to the reference grid.
    """
    for i_load in range(load.n):
        nid_ref = load.node_id[i_load]
        cid = load.coord_id[i_load]
        a_scale = load.scale_factor[i_load]
        r1, r2, r3 = load.r123[i_load]
        method = load.method[i_load] if hasattr(load, 'method') else 1
        racc = load.racc[i_load] if hasattr(load, 'racc') else 0.0

        # Angular velocity vector
        omega = a_scale * np.array([r1, r2, r3])

        # Transform from CID to basic if needed
        if cid != 0:
            coord = model.coord
            beta = coord.xyz_to_global_transform[cid]
            omega = beta.T @ omega

        # Angular acceleration (for tangential component)
        alpha = racc * np.array([r1, r2, r3])
        if cid != 0:
            alpha = beta.T @ alpha

        omega *= scale
        alpha *= scale

        log.debug(f"  RFORCE: omega={omega}, alpha={alpha}, nid_ref={nid_ref}")

        # Reference point position
        grid = model.grid
        all_nids = grid.node_id
        xyz_cid0 = grid.xyz_cid0()

        if nid_ref > 0:
            iref = np.searchsorted(all_nids, nid_ref)
            xyz_ref = xyz_cid0[iref]
        else:
            xyz_ref = np.zeros(3)

        # Build temporary Mbb for force computation
        temp_subcase = Subcase(id=0)
        Mbb_temp = build_Mbb(model, temp_subcase, dof_map, ndof, fdtype="float64")

        # Apply F = M * a_centrifugal at each grid
        for nid, xyz_i in zip(all_nids, xyz_cid0):
            nid_int = int(nid)
            i1 = dof_map[(nid_int, 1)]

            # Relative position
            r_i = xyz_i - xyz_ref

            # Centripetal: omega x (omega x r_i)
            omega_cross_r = np.cross(omega, r_i)
            accel_centripetal = np.cross(omega, omega_cross_r)

            # Tangential: alpha x r_i
            accel_tangential = np.cross(alpha, r_i)

            # Total acceleration (negative for inertial force direction in Nastran)
            accel_total = -(accel_centripetal + accel_tangential)

            # Mass at this grid
            m_diag = np.array([
                Mbb_temp[i1, i1],
                Mbb_temp[i1 + 1, i1 + 1],
                Mbb_temp[i1 + 2, i1 + 2],
            ])
            Fb[i1: i1 + 3] += m_diag * accel_total
