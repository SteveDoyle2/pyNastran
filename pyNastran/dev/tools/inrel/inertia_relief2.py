"""
Inertia relief on a Nastran model (principal-axis formulation).

Purpose
-------
Given a set of lumped masses and an unbalanced set of applied loads, compute
the d'Alembert (inertial) loads that place the free-free structure in
equilibrium, so that a static solve is well posed without artificial
constraints.  The caller adds the returned deltas to the applied loads::

    force_out  = force  + dforce
    moment_out = moment + dmoment
    sum(force_out) == 0
    sum(moment_out) + sum(r × force_out, about any point) == 0

Method
------
The rigid-body equations of motion for a free body about the cg are::

    sum(F) = m * a                 (uncoupled, easy)
    sum(M) = [I] * alpha           ([I] is full -> coupled)

``[I]`` generally has non-zero products of inertia, so ``alpha`` requires a
3x3 solve.  This module instead rotates into the **principal axis frame**,
where the cross-inertia terms vanish and the moment equation decouples into
three independent scalar equations::

    [I]        = [S] [I_principal] [S].T       (symmetric eigendecomposition)
    M_t        = transform(M)
    alpha_t[i] = M_t[i] / I_principal[i]        (no matrix solve needed)

``numpy.linalg.eigh`` supplies ``S`` (columns are the eigenvectors) and is used
rather than a general eigensolver because ``[I]`` is symmetric by construction.
The loads are transformed in, the inertial relief is computed axis-by-axis, and
the resulting deltas are transformed back to the global frame.

Sign convention
---------------
Per the Nastran QRG for the CONM2, the 6x6 mass matrix about the cg is::

    [M  0  0                ]
    |0  M  0                |
    |0  0  M                |
    |0  0  0   I11 -I12 -I13|
    |0  0  0  -I12  I22 -I23|
    [0  0  0  -I13 -I23  I33]

so the assembled inertia *tensor* carries **negative** off-diagonal terms,
while the packed 6-term array ``[Ixx, Iyy, Izz, Ixy, Ixz, Iyz]`` used
throughout pyNastran stores the products of inertia as the **positive**
integrals ``Ixy = sum(m*dx*dy)``.  The tensor is therefore assembled as::

    [ Ixx  -Ixy  -Ixz]
    [-Ixy   Iyy  -Iyz]
    [-Ixz  -Iyz   Izz]

which is the form ``_mass_properties_total`` builds.  This is the physically
correct ``[I] = sum m ((r.r) I3 - r r.T)``; verified against that closed form.

CONM2 card values follow QRG CONM2 remark 5: ``I21 = int(rho x1 x2 dV)`` (a
positive integral) and Nastran supplies the negative sign when it forms the
matrix (``-I21`` off-diagonal).  A CONM2's I21/I31/I32 therefore map directly
to Ixy/Ixz/Iyz in the packed array here, and to ``-I21`` etc. in the tensor.
``pyNastran.bdf.mesh_utils.mass_properties`` uses the same packed convention
(fixed 2026-11-17; it previously packed a CONM2's own products as ``-I21``).

Assumptions
-----------
- **Rigid body.**  Only the six rigid-body accelerations are solved for;
  structural flexibility does not feed back into the load distribution.  This
  is the classic "rigid-body inertia relief", not the flexible/residual-vector
  form (Nastran SOL 101 ``PARAM,INREL,-1`` / ``-2`` support both).
- **Quasi-static.**  Velocity-dependent terms are dropped.  There is no
  centrifugal ``omega x (omega x r)`` or Coriolis contribution, and no
  gyroscopic ``omega x [I] omega`` term in the moment balance.  Valid for a
  body released from rest; **not** valid for a spinning/maneuvering vehicle
  where ``omega`` is significant.
- **Small angles / no large rotation.**  The inertia tensor is evaluated once
  in the reference configuration and held fixed.
- **Lumped masses.**  Mass is a discrete set of points; the self-inertia of
  each point is optional and passed separately.
- **Free-free.**  The body is assumed fully unconstrained.  Any SPC/SUPORT in
  the source model is ignored here.

Limitations
-----------
- Degenerate mass distributions produce a singular inertia tensor.  A
  collinear (1-D) set of masses has one zero principal inertia, so the
  rotation about the bar axis is undetermined.  That axis receives no
  rotational relief; if the applied load carries a moment about it, a
  ``RuntimeWarning`` is issued and the moment is left in the residual.
  See the guard in ``inertia_relief``.

Defects fixed 2026-09-21
------------------------
The rotated test case (``test_inrel_bar_force_linear_rotated``) was previously
disabled and failing with a residual ``sum(moment_out)`` of 177.9.  Four
independent root causes were found and corrected:

1. **Forward/back transform mismatch.**  ``eigh`` returns eigenvectors as
   *columns*, so ``S.T @ [I] @ S`` is the diagonal form (confirmed: off-diagonal
   residual 8e-15, versus 1.0e+01 for ``S @ [I] @ S.T``).  Vectors must
   therefore map in as ``S.T @ v`` and back as ``S @ v_t``.  The code forward-
   transformed with ``S`` and returned with ``S.T`` -- the opposite pairing.
2. **Diagonal-only d'Alembert moment.**  The per-node restoring moment was
   built as three uncoupled terms ``-alpha[i] * m * (d_j^2 + d_k^2)``, i.e. the
   diagonal of the inertia tensor only.  The correct distribution is the full
   tensor form ``-m * ((alpha x d) x d)``.  On the rotated bar the full form
   balances to 7e-15 where the diagonal form left 4.77.
3. **Wrong principal inertia.**  ``inertia_final`` was taken from
   ``transform_inertia(..., coord2=CORD2R(...))``, which performs a
   parallel-axis shift to a new reference point rather than a diagonalization.
   On the rotated bar it returned ``[5803.6, 6249.8, 446.7]`` where the true
   principal inertias are ``[0, 6250, 6250]``.  Now taken as ``diag(Mtt)``.
4. **Degenerate axis not detected.**  A collinear bar has a principal inertia
   that is zero in exact arithmetic but ``~1e-12`` out of ``eigh``, so the
   exact ``!= 0`` guard did not exclude it and produced a huge bogus ``alpha``.
   Now excluded by a tolerance relative to the largest principal inertia.

A fifth issue was fixed while validating: per-node **self-inertia** is supplied
in the global frame but was contracted directly with the principal-frame
``alpha``, mixing frames.  Each node's tensor is now rotated by ``S.T I S``
first.  This left a residual of 0.27 on a random cloud with self-inertia.

None of these were visible in the three original tests, which use a bar along
+x where the principal axes coincide with the global axes, ``S`` is the
identity, and all products of inertia vanish.

Validation: force and moment both balance to ~1e-13 over 300 random 3-D mass
clouds, and on collinear, planar, axisymmetric (degenerate eigenvalue pair),
coincident-node, large-magnitude, and non-zero-self-inertia cases.

Remaining limitations
---------------------
- ``_mass_properties_total`` skips diagonalization when
  ``epsilon/delta <= 0.001`` (the Nastran GPWG check), so very small products
  of inertia are intentionally treated as zero.

Rotational relief as forces (2026-11-17)
----------------------------------------
The angular-acceleration relief on each lumped mass is returned as the
tangential nodal *force* ``-m (alpha x d)`` in ``dforce``, consistent with
``inertia_loads.inertia_loads``.  It was previously returned as the
statically-equivalent nodal *moment* ``-m d x (alpha x d)``, which grids with
no rotational stiffness (solids, rods) cannot react.  ``dmoment`` now holds
only the self-inertia term ``-[I0] alpha``.

Interface
---------
Array-only: nothing here imports ``pyNastran.bdf``.  Model adapters build the
``(mass, xyz, inertia)`` arrays and call ``inertia_relief``; see
``inertia_relief_vectorized.py`` for the ``bdf_vectorized3`` adapter.

- ``inertia_relief`` / ``inertia_relief_6x6``: per-grid relief loads
  (point masses / 6x6 lumped blocks).
- ``inertia_relief_grav_rforce``: the GRAV acceleration and RFORCE angular
  acceleration that balance a load, so a few cards replace per-grid
  FORCE/MOMENT cards.  Validated against NX 2412: the balanced LOAD left
  1e-3 of SPC force at a single constrained grid vs 3430 unbalanced.
- ``pyNastran.bdf.mesh_utils.loads.rforce_alpha_operator``: Nastran's RFORCE
  angular-acceleration load (METHOD=1 measured from NX OLOAD; METHOD=2 per
  the QRG).  It lives in the non-dev ``loads.py`` so ``sum_forces_moments``
  can use it too.

"""

from __future__ import annotations
import warnings

import numpy as np
from scipy import sparse

from pyNastran.bdf.mesh_utils.loads import rforce_alpha_operator


def _mass_properties_total(
    mass: np.ndarray, xyz: np.ndarray, inertia: np.ndarray, debug: bool = False
) -> tuple[float, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """
    Assembles the total mass, cg, and inertia tensor about the cg, and finds
    the rotation into the principal axis frame.

    The per-node self-inertia is summed and the parallel-axis ``m*r^2`` terms
    are added, giving the tensor::

        [ Ixx  -Ixy  -Ixz]
        [-Ixy   Iyy  -Iyz]
        [-Ixz  -Iyz   Izz]

    Diagonalization is skipped when the body is already close to principal:
    the test is ``norm(off-diagonals) / norm(diagonal) > 0.001``, which is the
    Nastran GPWG epsilon/delta check, and emits the matching UWM 3042 text when
    it trips.  Below that threshold ``S`` is returned as the identity and the
    small products of inertia are intentionally treated as zero, consistent
    with Nastran.

    Parameters
    ----------
    mass : (nnode,) float ndarray
        mass per node (already scaled by wtmass)
    xyz : (nnode,3) float ndarray
        node locations in the global frame
    inertia : (nnode,6) float ndarray
        per-node self-inertia [Ixx, Iyy, Izz, Ixy, Ixz, Iyz]
    debug : bool; default=False
        prints the assembled inertia tensor

    Returns
    -------
    mass_total : float
        the summed mass
    cg_total : (3,) float ndarray
        the mass-weighted centroid in the global frame
    inertia_total : (6,) float ndarray
        the principal-frame inertia, repacked as
        [Ixx, Iyy, Izz, Ixy, Ixz, Iyz] with the products negated back to the
        positive-integral convention
    Mtt : (3,3) float ndarray
        the inertia tensor in the principal frame (diagonal, up to the
        epsilon/delta threshold)
    S : (3,3) float ndarray
        the eigenvector matrix from ``np.linalg.eigh``; **columns** are the
        principal directions, so ``S.T @ I @ S`` is the diagonal form

    Notes
    -----
    - The epsilon/delta threshold is intentional (Nastran GPWG behavior).
      Because ``inertia_relief`` only uses ``diag(Mtt)``, products of inertia
      below the threshold leave a small moment residual (order 1e-3 of the
      applied moment).
    - ``np.linalg.eigh`` sorts eigenvalues ascending, so the principal axis
      ordering is by magnitude and does **not** track the global x/y/z order.
      A near-axisymmetric body has a degenerate eigenvalue pair whose
      eigenvectors are arbitrary within the plane; the resulting ``S`` is
      then not unique and can flip between runs or between load cases.

    """
    mass_total = mass.sum()

    cg_total = (xyz * mass[:, np.newaxis]).sum(axis=0) / mass_total
    assert len(cg_total) == 3, cg_total
    inertia_totali = inertia.sum(axis=0)
    assert len(inertia_totali) == 6, inertia_totali

    dx = xyz[:, 0] - cg_total[np.newaxis, 0]
    dy = xyz[:, 1] - cg_total[np.newaxis, 1]
    dz = xyz[:, 2] - cg_total[np.newaxis, 2]
    # Ixx, Iyy, Izz, Ixy, Ixz, Iyz = inertia
    ixx = inertia_totali[0] + (mass * (dy**2 + dz**2)).sum()
    iyy = inertia_totali[1] + (mass * (dx**2 + dz**2)).sum()
    izz = inertia_totali[2] + (mass * (dx**2 + dy**2)).sum()
    ixy = inertia_totali[3] + (mass * (dx * dy)).sum()
    ixz = inertia_totali[4] + (mass * (dx * dz)).sum()
    iyz = inertia_totali[5] + (mass * (dy * dz)).sum()
    inertia_mat = np.array(
        [
            [ixx, -ixy, -ixz],
            [-ixy, iyy, -iyz],
            [-ixz, -iyz, izz],
        ]
    )
    if debug:  # pragma: no cover
        print(f"inertia_mat:\n{inertia_mat}")
    Mtt_ = inertia_mat
    # Nastran GPWG epsilon/delta check; the 0.001 threshold is intentional
    delta = np.linalg.norm(np.diag(Mtt_))
    epsilon = np.linalg.norm([Mtt_[0, 1], Mtt_[0, 2], Mtt_[1, 2]])
    if epsilon / delta > 0.001:
        unused_eigvals, S = np.linalg.eigh(Mtt_)
        # eigh returns an orthogonal matrix but not necessarily a proper
        # rotation (det may be -1).  Cross products are pseudovectors, so
        # (S.T a) × (S.T b) = det(S) · S.T(a × b).  An improper rotation
        # silently flips the sign of all moment-arm cross products, which
        # breaks callers that compute moments about different reference
        # points.  Flipping one column makes det(S) = +1 without changing
        # the eigenvalue decomposition (the eigenvector is still valid, just
        # with the opposite sign).
        if np.linalg.det(S) < 0:
            S[:, 0] = -S[:, 0]
        msg = (
            "*** USER WARNING MESSAGE 3042 MODULE = GPWG\n"
            f"INCONSISTENT SCALAR MASSES HAVE BEEN USED. EPSILON/DELTA = {epsilon / delta:.7E}\n"
        )
        print(msg)
        Mtt = S.T @ Mtt_ @ S  # Mt
    else:
        S = np.eye(3, dtype=Mtt_.dtype)
        Mtt = Mtt_.copy()
    inertia_total = np.array(
        [
            Mtt[0, 0],
            Mtt[1, 1],
            Mtt[2, 2],
            -Mtt[0, 1],
            -Mtt[0, 2],
            -Mtt[1, 2],
        ]
    )
    assert Mtt.shape == (3, 3), Mtt.shape
    assert inertia_total.shape == (6,), inertia_total.shape
    return mass_total, cg_total, inertia_total, Mtt, S


def inertia_relief(
    weight: np.ndarray,
    xyz: np.ndarray,
    inertia: np.ndarray,
    force: np.ndarray,
    moment: np.ndarray,
    wtmass: float = 1.0,
    debug: bool = False,
) -> tuple[np.ndarray, np.ndarray]:
    """
    Computes the d'Alembert inertial loads that balance an unbalanced set of
    applied loads on a free-free body.

    Works in the principal-axis frame so that the moment equation decouples
    (``M = I*alpha`` per axis, no 3x3 solve).  See the module docstring for the
    full method, assumptions, and known defects.

    Parameters
    ----------
    weight:  (nnode,) float ndarray
        weight or mass per node
        units of mass or weight (depending on wtmass)
        NOTE: the shape is (nnode,), not (nnode,3) as previously documented;
        it is broadcast against the (nnode,3) force array internally.
    xyz : (nnode,3) float ndarray
        xyz locations of each mass in the global frame
    inertia : (nnode,6) float ndarray
        self-inertia (don't include mr^2)
        packed as [Ixx, Iyy, Izz, Ixy, Ixz, Iyz] with products of inertia
        stored as the positive integrals sum(m*dx*dy)
        units of W*L^2 or M*L^2 depending on wtmass
    force : (nnode,3) float ndarray
        existing forces in the global frame
    moment : (nnode,3) float ndarray
        existing moments in the global frame
    wtmass : float; default=1.0
        weight to mass conversion
        1/g when the caller passes weight rather than mass (PARAM,WTMASS)
    debug : bool; default=False
        prints the assembled inertia tensor

    Returns
    -------
    dforce : (nnode,3) float ndarray
        the inertial force increment, in the global frame:
        ``-m * (a + alpha x d)`` with ``d`` the offset from the cg
    dmoment : (nnode,3) float ndarray
        the inertial moment increment, in the global frame:
        ``-[I0] alpha`` from each node's self-inertia only (zero for pure
        point masses)

    The caller is expected to sum these onto the applied loads::

        force_out  = force  + dforce      # sum(force_out) ~= 0
        moment_out = moment + dmoment     # pure nodal moments

    Equilibrium check (about any reference point ``p``)::

        sum(force_out) ~= 0
        sum(moment_out) + sum(r_p × force_out) ~= 0

    where ``r_p = xyz - p``.  ``dmoment`` does **not** include ``r × dforce``
    couples, so callers compute force-arm moments themselves (as section cuts
    and monitor stations naturally do).

    Assumptions and limitations
    ---------------------------
    - Rigid-body, quasi-static: no omega-dependent (centrifugal, Coriolis,
      gyroscopic) terms.  See the module docstring.
    - ``inertia`` is the *self*-inertia of each lumped mass.  The parallel-axis
      ``m*r^2`` contribution is added internally by ``_mass_properties_total``;
      passing it in ``inertia`` double-counts it.
    - A degenerate mass distribution (collinear masses, or a single point) has
      a zero principal inertia about at least one axis.  That axis is detected
      by a relative tolerance and simply receives no rotational relief, which
      is the physically correct result -- an axis with no inertia cannot
      absorb a moment.  If the applied load *does* carry a moment about such
      an axis the problem is ill-posed; a ``RuntimeWarning`` is issued and the
      moment is left in the residual rather than being silently absorbed.

    """
    mass = weight * wtmass
    inertia = inertia * wtmass
    mass_total, cg_total, inertia_total, Mtt, S = _mass_properties_total(
        mass, xyz, inertia, debug=debug
    )
    assert inertia_total.shape == (6,), inertia_total.shape

    dxyz = xyz - cg_total[np.newaxis, :]

    # Transform forces/moments/offsets into the principal inertia frame.
    #
    # np.linalg.eigh returns the eigenvectors as the COLUMNS of S, so the
    # diagonal form is S.T @ [I] @ S. To stay consistent, a vector maps into
    # the principal frame as v_t = S.T @ v and back out as v = S @ v_t.
    # (Using S in and S.T out -- the opposite pairing -- silently gives the
    # wrong alpha whenever the body is not already aligned with the global
    # axes; that is why the rotated-bar test used to fail.)
    fxyzt = np.einsum("nij,nj->ni", S.T[np.newaxis, :], force)
    mxyzt = np.einsum("nij,nj->ni", S.T[np.newaxis, :], moment)
    dxyzt = np.einsum("nij,nj->ni", S.T[np.newaxis, :], dxyz)

    shape = fxyzt.shape
    dtype = fxyzt.dtype

    assert cg_total.shape == (3,), cg_total.shape
    assert inertia_total.shape == (6,), inertia_total.shape

    # The principal inertias are simply the diagonal of Mtt = S.T @ [I] @ S,
    # which _mass_properties_total already computed.
    #
    # This previously went through transform_inertia(..., coord2=CORD2R(...))
    # built from the eigenvectors, but that performs a parallel-axis shift to
    # a new reference point -- a different operation -- and returned a wrong,
    # non-diagonal inertia (e.g. [5803.6, 6249.8, 446.7] where the true
    # principal values are [0, 6250, 6250]). It also tripped the
    # "coord2 xform not validated yet" warning.
    inertia_final = np.diag(Mtt).copy()
    assert len(inertia_final) == 3, inertia_final

    fxyzt_total = fxyzt.sum(axis=0)
    assert len(fxyzt_total) == 3, fxyzt_total
    # F = m*a
    accelt = fxyzt_total / mass_total

    fxyzt_delta = -mass[:, np.newaxis] * accelt[np.newaxis, :]
    fxyzt_out = fxyzt + fxyzt_delta

    # moment of the linearly-relieved applied forces about the cg
    mxyzt_df = np.cross(dxyzt, fxyzt_out, axis=1)

    assert fxyzt.shape == mxyzt.shape, (fxyzt.shape, mxyzt.shape)

    mxyzt_total = (mxyzt + mxyzt_df).sum(axis=0)

    # Don't divide by ~0.
    #
    # A degenerate (e.g. collinear) mass distribution has a principal inertia
    # that is zero in exact arithmetic but comes out of eigh as ~1e-12, so an
    # exact "!= 0" test does NOT exclude it and would produce a huge bogus
    # alpha. Use a tolerance relative to the largest principal inertia.
    inertia_max = max(inertia_final.max(), 1e-300)
    ipos = np.where(inertia_final > 1e-8 * inertia_max)[0]

    # M = I * alpha, one scalar divide per principal axis (this is the whole
    # point of working in the principal frame). Axes with no inertia get no
    # rotational relief.
    alphat_moment = np.zeros(3, dtype="float64")
    alphat_moment[ipos] = mxyzt_total[ipos] / inertia_final[ipos]

    # A skipped axis that still carries a moment is not physically solvable:
    # there is no inertia to react it, so that moment stays in the residual
    # and the caller's "balanced" load set is not actually balanced. Warn
    # rather than raise -- the other two axes are still perfectly valid.
    izero = np.setdiff1d(np.arange(3), ipos)
    if len(izero) and np.abs(mxyzt_total[izero]).max() > 1e-8 * max(
        np.abs(mxyzt_total).max(), 1e-300
    ):
        warnings.warn(
            f"degenerate mass distribution: principal inertia(s) {izero} are "
            f"~0 (inertia={inertia_final}), but the applied moment about them "
            f"is {mxyzt_total[izero]}. No rotational relief is possible about "
            "a zero-inertia axis, so that moment is left unbalanced.",
            RuntimeWarning,
            stacklevel=2,
        )

    # d'Alembert force from the rigid-body angular acceleration.
    #
    # A point mass at d from the cg has the tangential acceleration alpha x d,
    # so its inertial load is the FORCE
    #     -m * (alpha x d)
    # applied at the node. These sum to zero (sum(m*d) = 0 about the cg) and
    # their moment about the cg is -sum(m d x (alpha x d)) = -[I_pp] alpha,
    # the full tensor (products of inertia included).
    #
    # This was previously returned as the equivalent nodal MOMENT
    # -m * d x (alpha x d). That is statically equivalent for the free body,
    # but a point mass has no rotary inertia of its own, and grids with no
    # rotational stiffness (solids, rods) cannot react a nodal moment, so the
    # relief would be lost there.
    fxyzt_alpha = -mass[:, np.newaxis] * np.cross(alphat_moment[np.newaxis, :], dxyzt)
    fxyzt_delta = fxyzt_delta + fxyzt_alpha

    # Only the self-inertia of each lumped mass produces a nodal moment.
    mxyzt_dm = np.zeros(shape, dtype=dtype)

    # Self-inertia of each lumped mass (its own I0*alpha, no parallel axis).
    #
    # The per-node self-inertia is given in the GLOBAL frame, but alpha lives
    # in the principal frame, so each node's tensor has to be rotated before
    # it can be contracted with alpha. Applying the raw diagonal
    # (inertia[:, :3] * alpha) mixes frames and leaves a residual whenever the
    # self-inertia is non-zero and S is not the identity.
    if np.abs(inertia).max() > 0.0:
        ixx_n, iyy_n, izz_n, ixy_n, ixz_n, iyz_n = inertia.T
        # assemble (nnode,3,3) with the negative-off-diagonal tensor convention
        imat_n = np.empty((len(mass), 3, 3), dtype="float64")
        imat_n[:, 0, 0] = ixx_n
        imat_n[:, 1, 1] = iyy_n
        imat_n[:, 2, 2] = izz_n
        imat_n[:, 0, 1] = imat_n[:, 1, 0] = -ixy_n
        imat_n[:, 0, 2] = imat_n[:, 2, 0] = -ixz_n
        imat_n[:, 1, 2] = imat_n[:, 2, 1] = -iyz_n
        # rotate into the principal frame: I_t = S.T @ I @ S
        imat_t = np.einsum("ij,njk,kl->nil", S.T, imat_n, S)
        mxyzt_dm = mxyzt_dm - np.einsum("nij,j->ni", imat_t, alphat_moment)

    # mxyzt_df is the per-node r×F couple about the CG — needed internally
    # to solve for alpha, but NOT part of the returned dmoment.  Including
    # it would bake force-arm moments into the nodal moment correction,
    # which forces every caller that computes moments about a different
    # point (section cuts, monitor stations) to reverse-engineer and
    # subtract it.  Returning only the self-inertia I0*alpha correction
    # (mxyzt_dm) keeps the output physically clean:
    #   sum(moment + dmoment) + sum(r × (force + dforce)) == 0
    # about any reference point.
    # Transform back out of the principal frame.
    # Forward was v_t = S.T @ v, so the inverse is v = S @ v_t.
    fxyz_delta_out = np.einsum("nij,nj->ni", S[np.newaxis, :], fxyzt_delta)
    mxyz_delta_out = np.einsum("nij,nj->ni", S[np.newaxis, :], mxyzt_dm)

    return fxyz_delta_out, mxyz_delta_out


def lumped_blocks_to_point_masses(
    mass_blocks: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Converts per-grid 6x6 rigid mass blocks into equivalent point masses.

    A 6x6 rigid-body mass block about a grid is exactly equivalent to a single
    mass ``m`` located at an offset ``e`` from the grid with a self-inertia
    ``I0`` about that offset point::

        [ m I3       m skew(e).T             ]
        [ m skew(e)  I0 + m((e.e) I3 - e e.T)]

    so a CONM2 offset (or several masses lumped at one grid) is represented
    without loss for rigid-body purposes.

    Parameters
    ----------
    mass_blocks : (ngrid,6,6) float ndarray
        the diagonal 6x6 block of the basic-frame mass matrix for each grid,
        ordered [T1, T2, T3, R1, R2, R3]

    Returns
    -------
    mass : (ngrid,) float ndarray
        the translational mass; the mean of the T1/T2/T3 diagonal
    offset : (ngrid,3) float ndarray
        the offset ``e`` from the grid to the mass cg (0 where mass is 0)
    inertia : (ngrid,6) float ndarray
        self-inertia about the mass cg, packed as
        [Ixx, Iyy, Izz, Ixy, Ixz, Iyz] (positive-integral products)

    Notes
    -----
    - A block whose T1/T2/T3 masses differ (e.g. a CMASS on one DOF) has no
      single point-mass equivalent; the mean is used.  Callers should check
      ``mass_blocks[:, :3, :3]`` if they care.
    - The rotational block is the tensor, so its off-diagonals are
      ``-I21``-style terms (QRG CONM2 remark 5); they are negated when packed
      so ``inertia`` holds the positive integrals ``Ixy = int(rho x y dV)``.
    - The first moment ``m e`` is read from the lower-left (``R x T``) block
      using the CONM2 sign convention in ``mass_matrix._add_conm2``:
      ``M[3,1] = -m e3, M[3,2] = m e2, M[4,2] = -m e1``.

    """
    nblock = mass_blocks.shape[0]
    assert mass_blocks.shape == (nblock, 6, 6), mass_blocks.shape
    mass = (mass_blocks[:, 0, 0] + mass_blocks[:, 1, 1] + mass_blocks[:, 2, 2]) / 3.0

    # first moment m*e from the symmetric halves of the R x T coupling block
    rt = 0.5 * (mass_blocks[:, 3:, :3] + np.swapaxes(mass_blocks[:, :3, 3:], 1, 2))
    first_moment = np.column_stack([-rt[:, 1, 2], rt[:, 0, 2], -rt[:, 0, 1]])
    offset = np.zeros((nblock, 3), dtype="float64")
    imass = mass != 0.0
    offset[imass] = first_moment[imass] / mass[imass, np.newaxis]

    # remove the parallel-axis term to get the inertia about the mass cg
    e = offset
    ee = np.einsum("ni,nj->nij", e, e)
    parallel = mass[:, np.newaxis, np.newaxis] * (
        np.einsum("ni,ni->n", e, e)[:, np.newaxis, np.newaxis] * np.eye(3)[np.newaxis, :, :] - ee
    )
    rr = mass_blocks[:, 3:, 3:]
    i0 = 0.5 * (rr + np.swapaxes(rr, 1, 2)) - parallel
    inertia = np.column_stack(
        [i0[:, 0, 0], i0[:, 1, 1], i0[:, 2, 2], -i0[:, 0, 1], -i0[:, 0, 2], -i0[:, 1, 2]]
    )
    return mass, offset, inertia


def inertia_relief_6x6(
    mass_blocks: np.ndarray,
    xyz: np.ndarray,
    force: np.ndarray,
    moment: np.ndarray,
    debug: bool = False,
) -> tuple[np.ndarray, np.ndarray]:
    """
    Inertia relief from per-grid 6x6 lumped mass blocks.

    This is the model-agnostic entry point for a lumped mass matrix (e.g. a
    block-diagonal MGG in the basic frame).  Each block is converted to a
    point mass at its own cg (``lumped_blocks_to_point_masses``), the applied
    loads are moved to those points, ``inertia_relief`` is solved, and the
    relief is moved back to the grids::

        M_p = M_g - e x F         (grid -> mass cg)
        dM_g = dM_p + e x dF      (mass cg -> grid)

    Parameters
    ----------
    mass_blocks : (ngrid,6,6) float ndarray
        per-grid 6x6 mass blocks in the basic frame (already scaled by wtmass)
    xyz : (ngrid,3) float ndarray
        grid locations in the basic frame
    force : (ngrid,3) float ndarray
        applied grid forces in the basic frame
    moment : (ngrid,3) float ndarray
        applied grid moments in the basic frame
    debug : bool; default=False
        passed to ``inertia_relief``

    Returns
    -------
    dforce : (ngrid,3) float ndarray
        the inertial force increment at each grid, basic frame
    dmoment : (ngrid,3) float ndarray
        the inertial moment increment at each grid, basic frame

    Notes
    -----
    - Only the diagonal 6x6 blocks are used; mass coupling *between* grids
      (e.g. a CMASS2 connecting two grids) is ignored.
      TODO: fold off-block-diagonal mass terms into the rigid-body properties.

    """
    mass, offset, inertia = lumped_blocks_to_point_masses(mass_blocks)
    xyz_mass = xyz + offset
    moment_mass = moment - np.cross(offset, force)
    dforce, dmoment_mass = inertia_relief(
        mass, xyz_mass, inertia, force, moment_mass, wtmass=1.0, debug=debug
    )
    dmoment = dmoment_mass + np.cross(offset, dforce)
    return dforce, dmoment


def _skew_stack(xyz: np.ndarray) -> np.ndarray:
    """(n,3) -> (n,3,3) skew matrices with ``skew(r) @ v == cross(r, v)``."""
    skew = np.zeros((len(xyz), 3, 3), dtype="float64")
    x, y, z = xyz[:, 0], xyz[:, 1], xyz[:, 2]
    skew[:, 0, 1], skew[:, 0, 2] = -z, y
    skew[:, 1, 0], skew[:, 1, 2] = z, -x
    skew[:, 2, 0], skew[:, 2, 1] = -y, x
    return skew


def inertia_relief_grav_rforce(
    mbb: np.ndarray,
    xyz: np.ndarray,
    force: np.ndarray,
    moment: np.ndarray,
    rforce_method: int = 1,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Finds a GRAV acceleration and an RFORCE angular acceleration that balance
    an applied load, so a few cards replace per-grid FORCE/MOMENT cards.

    Each card's load is modeled the way Nastran applies it:

    - GRAV (QRG GRAV remark 1): the load is ``+[M] D_t g`` with the full mass
      matrix, so CONM2 offsets and coupled mass contribute moments.
    - RFORCE with A=0 about the basic origin: the per-grid load Nastran
      applies for the given METHOD (see ``rforce_alpha_operator``).

    The six unknowns ``(g, alpha)`` are solved so the resultant of
    applied + GRAV + RFORCE is zero about the basic origin::

        [sum F_grav(g) + sum F_rforce(alpha)] = -sum F_applied
        [sum M_grav(g) + sum M_rforce(alpha)] = -sum M_applied

    Parameters
    ----------
    mbb : (6*ngrid, 6*ngrid) float ndarray or sparse matrix
        basic-frame mass matrix (already scaled by WTMASS), 6 DOFs per grid
        in ``xyz`` order
    xyz : (ngrid,3) float ndarray
        grid locations in the basic frame
    force : (ngrid,3) float ndarray
        applied grid forces in the basic frame
    moment : (ngrid,3) float ndarray
        applied grid moments in the basic frame
    rforce_method : int; default=1
        the RFORCE METHOD field (1 or 2); NX rejects METHOD=2 when any CONM2
        has an offset

    Returns
    -------
    grav_accel : (3,) float ndarray
        the GRAV acceleration vector ``A * N`` in the basic frame
    alpha : (3,) float ndarray
        the RFORCE angular acceleration (rad/time^2) about the basic origin;
        RACC = ``|alpha| / (2 pi)`` with ``R = alpha / |alpha|``
    residual : (6,) float ndarray
        the remaining [force, moment] resultant about the origin; ~0 unless
        the mass distribution is degenerate (e.g., collinear)

    Notes
    -----
    - The per-grid loads differ from ``inertia_relief`` at CONM2s with offsets
      (and, for METHOD=2, rotary inertia), because Nastran's RFORCE models
      them differently; the resultant is still exactly balanced because the
      solve uses Nastran's RFORCE load.
    - A degenerate (e.g., collinear) mass distribution can't react a moment
      about the zero-inertia axis.  The system is solved by least squares and
      a ``RuntimeWarning`` reports the residual.

    """
    ngrid = len(xyz)
    assert xyz.shape == (ngrid, 3), xyz.shape
    assert force.shape == (ngrid, 3), force.shape
    assert moment.shape == (ngrid, 3), moment.shape
    assert mbb.shape == (6 * ngrid, 6 * ngrid), (mbb.shape, ngrid)
    skew_r = _skew_stack(xyz)

    # GRAV: grid loads = Mbb @ (unit translation in direction j), j = 1..3
    dt = np.zeros((6 * ngrid, 3), dtype="float64")
    for j in range(3):
        dt[j::6, j] = 1.0
    grav_loads = np.asarray(mbb @ dt).reshape(ngrid, 6, 3)
    grav_f = grav_loads[:, :3, :]  # (ngrid, 3, 3): force per unit g
    grav_m = grav_loads[:, 3:, :]
    a_grav = np.vstack(
        [
            grav_f.sum(axis=0),
            (grav_m + skew_r @ grav_f).sum(axis=0),
        ]
    )

    # RFORCE: per-grid load per unit alpha, as Nastran applies it
    mass_blocks, unused_coupling = get_mass_blocks(sparse.csr_matrix(mbb))
    rforce_op = rforce_alpha_operator(mass_blocks, xyz, rforce_method)
    rforce_f = rforce_op[:, :3, :]
    rforce_m = rforce_op[:, 3:, :]
    a_rforce = np.vstack(
        [
            rforce_f.sum(axis=0),
            (rforce_m + skew_r @ rforce_f).sum(axis=0),
        ]
    )

    amat = np.hstack([a_grav, a_rforce])
    applied = np.hstack(
        [
            force.sum(axis=0),
            moment.sum(axis=0) + np.cross(xyz, force).sum(axis=0),
        ]
    )
    solution, *unused = np.linalg.lstsq(amat, -applied, rcond=1e-12)
    residual = applied + amat @ solution
    if np.abs(residual).max() > 1e-8 * max(np.abs(applied).max(), 1e-300):
        warnings.warn(
            f"GRAV/RFORCE cannot fully balance the load (degenerate mass "
            f"distribution?); residual [F, M] = {residual}",
            RuntimeWarning,
            stacklevel=2,
        )
    return solution[:3], solution[3:], residual


def get_mass_blocks(mbb: sparse.spmatrix) -> tuple[np.ndarray, float]:
    """
    Extracts the diagonal 6x6 block of each grid from a g-set mass matrix.

    Parameters
    ----------
    mbb : (6*ngrid, 6*ngrid) sparse matrix
        the mass matrix with 6 DOFs per grid in [T1, T2, T3, R1, R2, R3] order

    Returns
    -------
    mass_blocks : (ngrid,6,6) float ndarray
        the per-grid diagonal blocks
    coupling_ratio : float
        ``sum(|off-block terms|) / sum(|all terms|)``; 0.0 for a purely
        block-diagonal (lumped) matrix

    """
    ndof = mbb.shape[0]
    assert ndof % 6 == 0, ndof
    ngrid = ndof // 6
    coo = sparse.coo_matrix(mbb)
    igrid = coo.row // 6
    in_block = igrid == coo.col // 6

    mass_blocks = np.zeros((ngrid, 6, 6), dtype="float64")
    np.add.at(
        mass_blocks,
        (igrid[in_block], coo.row[in_block] % 6, coo.col[in_block] % 6),
        coo.data[in_block],
    )
    abs_total = np.abs(coo.data).sum()
    coupling_ratio = np.abs(coo.data[~in_block]).sum() / abs_total if abs_total > 0.0 else 0.0
    return mass_blocks, float(coupling_ratio)
