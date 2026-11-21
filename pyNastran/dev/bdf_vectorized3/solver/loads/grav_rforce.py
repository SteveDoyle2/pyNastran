import os
import copy
from datetime import date
from typing import TextIO, Any

import numpy as np
import scipy as sp
from scipy import sparse
from scipy.sparse import csc_matrix, lil_matrix, dok_matrix, issparse
from scipy.sparse.linalg import ArpackNoConvergence

from cpylog import SimpleLogger

from pyNastran.utils.solver_utils import get_solver
import pyNastran
from pyNastran.utils.numpy_utils import integer_types  # , float_types
from pyNastran.dev.bdf_vectorized3.bdf import BDF, Subcase
#from pyNastran.dev.bdf_vectorized3.cards.load.gravity import GRAV
from pyNastran.f06.errors import FatalError

from pyNastran.dev.bdf_vectorized3.solver.elements.beam import (
    beam_pg_distributed, beam_pg_point,
    consistent_mass, lumped_mass,)
from pyNastran.bdf.mesh_utils.loads import rforce_alpha_operator
DOF_MAP = dict[tuple[int, int], int]
Array = np.ndarray


def apply_grav(model: BDF,
               load: 'GRAV',
               scale: float,
               Fb: np.ndarray,
               dof_map: DOF_MAP, ndof: int,
               log: SimpleLogger,
               Mbb: Array) -> Array:
    """Apply GRAV load (gravity) to Fb using grid-level mass approach.

    Following MYSTRAN's approach: F_i = M_grid_i * acceleration.
    The assembled mass matrix is used to compute forces at each grid.
    """
    assert Mbb is not None, Mbb

    coord = model.coord
    betas = {}
    for cid in np.unique(load.coord_id):
        betas[cid] = coord.xyz_to_global_transform[cid]
            
    for i_load in range(load.n):
        sid = load.load_id[i_load]
        N_vec = load.N[i_load]  # direction vector (3,)
        if load.type == 'GRAV':
            mag = load.scale[i_load]  # accel
        else:
            # # scale factor (acceleration)
            mag = load.scale_factor[i_load]
            asdf
        cid = load.coord_id[i_load]

        # Acceleration vector in basic frame
        norm_N = np.linalg.norm(N_vec)
        if norm_N == 0.0:
            continue
        # QRG GRAV remark 1: a = A * N, so |a| = A * |N| (N is not normalized)
        accel_basic = mag * N_vec

        # Transform from CID to basic if needed
        if cid != 0:
            beta = betas[cid]
            accel_basic = beta.T @ accel_basic

        accel_basic *= scale

        log.debug(f"  GRAV: accel_basic={accel_basic}")

        # P = Mbb @ D a: every grid translates by the same acceleration.
        # Using the full Mbb (not just its diagonal) keeps the m*e x a moment
        # of offset masses (e.g. CONM2 X1-X3), so the gravity resultant acts
        # through the true cg.
        uaccel = np.zeros(ndof, dtype='float64')
        for nid in model.grid.node_id:
            i1 = dof_map[(nid, 1)]
            uaccel[i1: i1 + 3] = accel_basic
        Fb += Mbb @ uaccel


def apply_rforce(model: BDF,
                 load: 'RFORCE',
                 scale: float,
                 Fb: np.ndarray,
                 dof_map: DOF_MAP, ndof: int,
                 log: SimpleLogger,
                 xyz_cid0: np.ndarray | None = None,
                 Mbb: Array | None = None) -> None:
    """Apply RFORCE load (centrifugal/angular-acceleration force) to Fb.

    Per the QRG (RFORCE remarks 1-2)::

        F_i = -[m]_i [omega x (omega x (r_i - r_a)) + alpha x (r_i - r_a)]
        omega = 2*pi * A * R       (R is NOT normalized)
        alpha = 2*pi * RACC * R

    where ``[m]_i`` is the 3x3 translational mass at grid i and ``r_a`` is
    grid G (G=0 is the origin of the CID system).

    Notes
    -----
    - The angular-acceleration (RACC) part uses ``rforce_alpha_operator``,
      which reproduces NX OLOAD output for METHOD=1 (offset CONM2s and rotary
      inertia included, the way NX does it) and METHOD=2 (translational mass
      only).
    - The centrifugal (A) part uses only the translational mass at each grid.
      TODO: METHOD=1 CONM2 offset terms for the centrifugal (A) part.
    - The rotation vector is rotated from CID to basic by the coordinate
      system axes.

    """
    assert Mbb is not None, Mbb
    grid = model.grid
    all_nids = grid.node_id
    if xyz_cid0 is None:
        xyz_cid0 = grid.xyz_cid0()

    # per-grid 6x6 mass blocks; translational mass is the T1-T3 diagonal
    i1 = np.array([dof_map[(nid, 1)] for nid in all_nids])
    local = np.arange(6)
    rows = (i1[:, None, None] + local[None, :, None]) * np.ones((1, 1, 6), dtype=int)
    cols = (i1[:, None, None] + local[None, None, :]) * np.ones((1, 6, 1), dtype=int)
    mass_blocks = np.asarray(
        sparse.csr_matrix(Mbb)[rows.ravel(), cols.ravel()]).reshape(len(i1), 6, 6)
    m_trans = np.diagonal(mass_blocks[:, :3, :3], axis1=1, axis2=2)

    coord = model.coord
    for i_load in range(load.n):
        nid_ref = load.node_id[i_load]
        cid = load.coord_id[i_load]
        r123 = load.r[i_load]

        # rev/time -> rad/time; R is not normalized (QRG remark 2)
        omega = 2 * np.pi * load.scale[i_load] * r123
        alpha = 2 * np.pi * load.racc[i_load] * r123

        xyz_ref = np.zeros(3)
        if cid != 0:
            beta = coord.xyz_to_global_transform[cid]
            omega = beta.T @ omega
            alpha = beta.T @ alpha
            if nid_ref == 0:
                # The QRG G field says G=0 is the origin of the CID system,
                # while remark 4 says the basic origin.  This follows the
                # field definition.
                icoord = np.where(coord.coord_id == cid)[0][0]
                xyz_ref = coord.origin[icoord]
        if nid_ref > 0:
            iref = np.searchsorted(all_nids, nid_ref)
            xyz_ref = xyz_cid0[iref]
        log.debug(f"  RFORCE: omega={omega}, alpha={alpha}, nid_ref={nid_ref}")

        dxyz = xyz_cid0 - xyz_ref[np.newaxis, :]
        # centrifugal: -m omega x (omega x r)
        centrifugal = np.cross(omega[np.newaxis, :], np.cross(omega[np.newaxis, :], dxyz))
        loads = np.zeros((len(i1), 6), dtype='float64')
        loads[:, :3] = -m_trans * centrifugal
        # angular acceleration, per the METHOD (see rforce_alpha_operator)
        method = int(load.method[i_load])
        loads += rforce_alpha_operator(mass_blocks, dxyz, method) @ alpha

        # A LOAD scale multiplies the load vector, not omega (the centrifugal
        # term is quadratic in omega; QRG RFORCE remarks 11 and 19)
        for j in range(6):
            Fb[i1 + j] += scale * loads[:, j]
