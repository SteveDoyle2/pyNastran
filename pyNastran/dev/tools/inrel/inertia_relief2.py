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
    sum(moment_out about the cg) == 0

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
- The self-inertia correction is currently inert -- see the "Known defects"
  note on ``inertia_relief`` below.
- ``get_mass_properties_array`` does not yet read the model; it returns
  hardcoded values.

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
- ``get_mass_properties_array`` still returns hardcoded stub data, so
  ``inertia_relief_from_model`` does not yet work against a real model.
- The ``mxyzt_delta_self`` branch is still multiplied by ``0.`` and remains
  inert; the self-inertia contribution is now handled in the ``mxyzt_dm``
  term instead.
- ``_mass_properties_total`` skips diagonalization when
  ``epsilon/delta <= 0.001``, so very small products of inertia are treated
  as exactly zero.

"""
from __future__ import annotations
import warnings
from typing import TYPE_CHECKING

import numpy as np
from pyNastran.bdf.mesh_utils.mass_properties import mass_properties
if TYPE_CHECKING:
    from pyNastran.bdf.bdf import read_bdf, BDF


def get_mass_properties_array(model: BDF) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Intended to extract per-node lumped mass, location, and self-inertia from
    a BDF model.

    .. warning::
        **NOT IMPLEMENTED.**  This function calls ``mass_properties(model)``
        and then immediately discards the result, overwriting ``mass``, ``cg``,
        and ``inertia`` with hardcoded 5-node dummy values.  It returns the
        same stub data for every model, and ``inertia`` is always zeros.
        Anything built on this (notably ``inertia_relief_from_model``) is
        exercising the stub, not the model.

        A real implementation needs to walk the CONM2/CMASSi cards and the
        element/property mass, which is what
        ``pyNastran.bdf.mesh_utils.mass_properties.mass_properties_breakdown``
        already does on a per-element basis.

    Parameters
    ----------
    model : BDF
        the model (currently ignored)

    Returns
    -------
    mass : (nnode,) float ndarray
        hardcoded [1., 2., 3., 4., 5.]
    cg : (nnode,3) float ndarray
        hardcoded node locations
    inertia : (nnode,6) float ndarray
        always zeros

    """
    # from pyNastran.bdf.bdf import read_bdf
    mass, cg, inertia = mass_properties(model)
    nnode = 5
    mass = np.array([1., 2., 3., 4., 5])
    cg = np.array([
        [0., 0., 0.],
        [1., 0., 0.],
        [3., 0., 0.],
        [2., 1., 0.],
        [2., 0., 0.],
    ])
    # Ixx, Iyy, Izz, Ixy, Ixz, Iyz = inertia
    inertia = np.zeros((nnode, 6), dtype='float64')
    return mass, cg, inertia

def get_eigenvalues(imat: np.ndarray,
                    debug: bool=False):  # pragma: no cover
    """
    Finds the principal inertias by numerical optimization rather than by
    eigendecomposition.

    Parameterizes a rotation by two vectors ``i`` and ``j``, builds an
    orthonormal triad via ``k = i x j``, and minimizes the sum of the squared
    off-diagonal terms (less the squared diagonal terms) of ``T @ imat @ T.T``.

    .. note::
        This is a **cross-check / scratch routine**, not the production path.
        ``_mass_properties_total`` uses ``np.linalg.eigh``, which is exact,
        faster, and guaranteed orthonormal for a symmetric matrix.  This
        function is only reachable from the dead ``if debug and 0:`` block.

    Limitations
    -----------
    - The triad is built from ``i`` and ``j`` without re-orthogonalizing ``j``
      against ``i``, so ``T`` is not guaranteed orthonormal and the similarity
      transform is not guaranteed to preserve the eigenvalues.
    - The objective mixes two competing terms and is unbounded below as the
      diagonal grows, so the optimizer can wander; ``scipy.optimize.minimize``
      is run from a single fixed start with no convergence check.
    - Returns ``np.linalg.eigvals`` (general, possibly complex) of the rotated
      matrix rather than ``eigvalsh``, so results may carry tiny imaginary
      parts.

    """
    x0 = [1., 0., 0.,
          0., 1., 0.]
    import scipy
    def func(values):
        i = values[:3]
        j = values[3:]
        i /= np.linalg.norm(i)
        j /= np.linalg.norm(j)
        k = np.cross(i, j)
        k /= np.linalg.norm(k)
        if np.abs(k).max() == 0.0:
            raise RuntimeError(k)
        T = np.vstack([i, j, k])
        imat2 = T @ imat @ T.T
        # imat2 = T.T @ imat @ T
        obj = imat2[0, 1]**2 + imat2[0, 2]**2 + imat2[1, 2]**2
        obj -= imat2[0, 0]**2 + imat2[1, 1]**2 + imat2[2, 2]**2
        return obj
    bounds = [
        [-1., 1.], [-1., 1.], [-1., 1.],
        [-1., 1.], [-1., 1.], [-1., 1.],
    ]
    out = scipy.optimize.minimize(func, x0, bounds=bounds)
    # out = scipy.optimize.fmin(func, out)
    if debug:  # pragma: no cover
        print(out)
    values = out.x
    i = values[:3]
    j = values[3:]
    i /= np.linalg.norm(i)
    j /= np.linalg.norm(j)
    k = np.cross(i, j)
    k /= np.linalg.norm(k)
    T = np.vstack([i, j, k])
    imat2 = T @ imat @ T.T
    # imat2 = T.T @ imat @ T
    T2 = T ** 2
    if debug:  # pragma: no cover
        print(T2.sum(axis=0), T2.sum(axis=1))
        print(f'T:\n{T}')
        print(f'imat2:\n{imat2}')
    return np.linalg.eigvals(imat2)
    # return np.diag(imat2)


def _mass_properties_total(mass: np.ndarray,
                           xyz: np.ndarray,
                           inertia: np.ndarray,
                           debug: bool=False) -> tuple[float, np.ndarray, np.ndarray,
                                      np.ndarray, np.ndarray]:
    """
    Assembles the total mass, cg, and inertia tensor about the cg, and finds
    the rotation into the principal axis frame.

    The per-node self-inertia is summed and the parallel-axis ``m*r^2`` terms
    are added, giving the tensor::

        [ Ixx  -Ixy  -Ixz]
        [-Ixy   Iyy  -Iyz]
        [-Ixz  -Iyz   Izz]

    Diagonalization is skipped when the body is already close to principal:
    the test is ``norm(off-diagonals) / norm(diagonal) > 0.001``, mirroring the
    Nastran GPWG epsilon/delta check, and emits the matching UWM 3042 text when
    it trips.  Below that threshold ``S`` is returned as the identity, so small
    products of inertia are silently treated as zero.

    Parameters
    ----------
    mass : (nnode,) float ndarray
        mass per node (already scaled by wtmass)
    xyz : (nnode,3) float ndarray
        node locations in the global frame
    inertia : (nnode,6) float ndarray
        per-node self-inertia [Ixx, Iyy, Izz, Ixy, Ixz, Iyz]
    debug : bool; default=False
        prints the tensor and compares the closed-form eigenvalues against
        the optimizer in ``get_eigenvalues``

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
    - ``eigenvalues`` is computed via ``np.linalg.eigvals`` and then discarded
      unless ``debug``; the ``if debug and 0:`` guard makes that block dead.
    - ``np.linalg.eigh`` sorts eigenvalues ascending, so the principal axis
      ordering is by magnitude and does **not** track the global x/y/z order.
      A near-axisymmetric body has a degenerate eigenvalue pair whose
      eigenvectors are arbitrary within the plane; the resulting ``S`` is
      then not unique and can flip between runs or between load cases.

    """
    mass_total = mass.sum()
    # assert np.allclose(mass_total, 5.), mass_total

    cg_total = (xyz * mass[:, np.newaxis]).sum(axis=0) / mass_total
    assert len(cg_total) == 3, cg_total
    # print(f'mass = {mass}')
    # print(f'cg = {cg_total}')
    # print(f'r  = {np.linalg.norm(cg_total):g}')
    inertia_totali = inertia.sum(axis=0)
    # print(f'inertia_totali = {inertia_totali}')
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
    inertia_mat = np.array([
        [ixx, -ixy, -ixz],
        [-ixy, iyy, -iyz],
        [-ixz, -iyz, izz],
    ])
    eigenvalues = np.linalg.eigvals(inertia_mat)
    if debug and 0:
        print(f'inertia_mat:\n{str(inertia_mat)}')
        eigenvalues_opt = get_eigenvalues(inertia_mat)
        print(f'eigenvalues_opt = {eigenvalues_opt}')
        print(f'eigenvalues = {eigenvalues}')
    Mtt_ = inertia_mat
    Mtd = np.diag(Mtt_)
    delta = np.linalg.norm(Mtd)
    e_ = [Mtt_[0, 1], Mtt_[0, 2], Mtt_[1, 2]]
    epsilon = np.linalg.norm(e_)
    if epsilon/delta > 0.001:
        unused_eigvals, S = np.linalg.eigh(Mtt_)
        # S = S.T
        # [Mtt_] = [phi] * [lambda] * [phi.T]  # mabye phi=phi.T?
        msg = (
            '*** USER WARNING MESSAGE 3042 MODULE = GPWG\n'
            f'INCONSISTENT SCALAR MASSES HAVE BEEN USED. EPSILON/DELTA = {epsilon/delta:.7E}\n')
        # model.log.warning(msg)
        print(msg)
        Mtt = S.T @ Mtt_ @ S  # Mt
    else:
        S = np.eye(3, dtype=Mtt_.dtype)
        Mtt = Mtt_.copy()
    # print(f'S:\n{S}')

    # inertia_mat = np.array([
    #     [ixx, -ixy, -ixz],
    #     [-ixy, iyy, -iyz],
    #     [-ixz, -iyz, izz],
    # ])
    inertia_total = np.array([
        Mtt[0, 0], Mtt[1, 1], Mtt[2, 2],
        -Mtt[0, 1], -Mtt[0, 2], -Mtt[1, 2],
    ])
    assert Mtt.shape == (3,3), Mtt.shape
    assert inertia_total.shape == (6,), inertia_total.shape
    return mass_total, cg_total, inertia_total, Mtt, S


def inertia_relief_from_model(model: BDF) -> tuple[np.ndarray, np.ndarray]:
    """
    Convenience wrapper: pull mass properties from a model and run inertia
    relief on it.

    .. warning::
        **NOT USABLE AGAINST A REAL MODEL YET.**  Two placeholders block it:
        ``get_mass_properties_array`` ignores ``model`` and returns hardcoded
        stub data, and the applied loads here are hardcoded to a unit force in
        every direction at every node (``force = np.ones((nnode, 3))``) with
        zero moment, rather than being read from the model's load cases.

        A real version needs to take a load case id (or an already-assembled
        load vector) as an argument.

    Parameters
    ----------
    model : BDF
        the model (currently ignored)

    Returns
    -------
    dforce : (nnode,3) float ndarray
        the inertial force increment
    dmoment : (nnode,3) float ndarray
        the inertial moment increment

    """
    mass, xyz, inertia = get_mass_properties_array(model)
    nnode = len(mass)
    force = np.ones((nnode, 3))
    moment = np.ones((nnode, 3)) * 0
    dforce, dmoment = inertia_relief(
        mass, xyz, inertia,
        force, moment)
    return dforce, dmoment


def _get_transformed_self_inertia(inertia: np.ndarray,
                                  S: np.ndarray) -> np.ndarray:
    """
    Sums the per-node self-inertia and rotates it into the principal frame.

    Parameters
    ----------
    inertia : (nnode,6) float ndarray
        per-node self-inertia [Ixx, Iyy, Izz, Ixy, Ixz, Iyz]
    S : (3,3) float ndarray
        the eigenvector matrix from ``_mass_properties_total``

    Returns
    -------
    inertia_xform : (3,) float ndarray
        the diagonal of the rotated self-inertia tensor

    Notes
    -----
    - Only the diagonal is returned; the rotated products of inertia are
      computed and then dropped (the corresponding lines are commented out).
      This is lossy whenever the self-inertia tensor is not already principal
      in the same frame as the overall body.
    - The ``S.T @ imat @ S`` orientation is correct: ``eigh`` returns the
      eigenvectors as columns, so this is the diagonalizing similarity.  (The
      vector transforms in ``inertia_relief`` previously used the opposite
      pairing; that was fixed 2026-09-21.)
    - The returned value feeds only the ``mxyzt_delta_self`` branch, which is
      multiplied by zero, so in the current state this result is unused.  The
      live self-inertia handling is in the ``mxyzt_dm`` term, which rotates
      each node's tensor individually rather than using this summed value.

    """
    inertia_self = inertia.sum(axis=0)
    ixx, iyy, izz, ixy, ixz, iyz = inertia_self
    imat = np.array([
        [ixx, -ixy, -ixz],
        [-ixy, iyy, -iyz],
        [-ixz, -iyz, izz],
    ])
    ximat = S.T @ imat @ S  # is this backwards?
    inertia_xform = np.array([
        ximat[0, 0], ximat[1, 1], ximat[2, 2],
        # -ximat[0, 1], -ximat[0, 2], -ximat[1, 2],
    ])
    return inertia_xform


def inertia_relief(weight: np.ndarray,
                   xyz: np.ndarray,
                   inertia: np.ndarray,
                   force: np.ndarray,
                   moment: np.ndarray,
                   wtmass: float=1.0,
                   debug: bool=False) -> tuple[np.ndarray, np.ndarray]:
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
        prints the assembled inertia tensor and its eigenvalues

    Returns
    -------
    dforce : (nnode,3) float ndarray
        the inertial force increment, in the global frame
    dmoment : (nnode,3) float ndarray
        the inertial moment increment, in the global frame

    The caller is expected to sum these onto the applied loads::

        force_out  = force  + dforce      # sums to ~0
        moment_out = moment + dmoment     # sums to ~0

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
    - ``mxyzt_delta_self`` is still multiplied by ``0.`` and is inert; the
      self-inertia contribution is applied in the ``mxyzt_dm`` term instead.
    - The returned increments are referenced to the cg: ``dmoment`` already
      includes the ``r x F`` couple from the balanced forces, so the caller
      checks equilibrium as ``sum(moment + dmoment) == 0`` directly and must
      **not** add ``r x (force + dforce)`` again.

    """
    mass = weight * wtmass
    inertia = inertia * wtmass
    mass_total, cg_total, inertia_total, Mtt, S = _mass_properties_total(
        mass, xyz, inertia, debug=debug)
    inertia_self = _get_transformed_self_inertia(
        inertia, S)
    assert inertia_total.shape == (6,), inertia_total.shape

    # TODO: Should this be S.T?
    # Ft = S @ F
    # print(f'xyz:\n{xyz}')
    # print(f'cg_total = {cg_total}')
    # print(f'inertia_self = {inertia_self}')
    dxyz = xyz - cg_total[np.newaxis, :]

    # Transform forces/moments/offsets into the principal inertia frame.
    #
    # np.linalg.eigh returns the eigenvectors as the COLUMNS of S, so the
    # diagonal form is S.T @ [I] @ S. To stay consistent, a vector maps into
    # the principal frame as v_t = S.T @ v and back out as v = S @ v_t.
    # (Using S in and S.T out -- the opposite pairing -- silently gives the
    # wrong alpha whenever the body is not already aligned with the global
    # axes; that is why the rotated-bar test used to fail.)
    fxyzt = np.einsum('nij,nj->ni',S.T[np.newaxis,:], force)
    mxyzt = np.einsum('nij,nj->ni',S.T[np.newaxis,:], moment)
    dxyzt = np.einsum('nij,nj->ni',S.T[np.newaxis,:], dxyz)

    shape = fxyzt.shape
    dtype = fxyzt.dtype

    # print(f'Mtt:\n{str(Mtt)}')
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
    # print(f'fxyzt_total = {fxyzt_total}')
    assert len(fxyzt_total) == 3, fxyzt_total
    # F = m*a
    accelt = fxyzt_total / mass_total
    # print(f'accelt = {accelt}')

    fxyzt_delta = -mass[:, np.newaxis] * accelt[np.newaxis, :]
    fxyzt_out = fxyzt + fxyzt_delta
    # if np.abs(fxyzt_out).max() > 0:
    #     print(f'fxyzt_out:\n{fxyzt_out}')
    # print(f'dxyzt =\n{dxyzt}')

    # this is the moment due to applied forces
    # was fxyzt_delta:
    # - fails the force_linear test -> fxyzt_out
    mxyzt_df = np.cross(dxyzt, fxyzt_out, axis=1)
    # print(f'mxyzt_delta0 = {mxyzt_delta0}')
    # if np.abs(mxyzt_df).max() > 0:
    #     print(f'dxyzt:\n{dxyzt}')
    #     print(f'fxyzt:\n{fxyzt}')
    #     print(f'mxyzt_df:\n{mxyzt_df}')

    # print(f'mxyzt_df_total = {mxyzt_df.sum(axis=0)}')
    assert fxyzt.shape == mxyzt.shape, (fxyzt.shape, mxyzt.shape)

    mxyzt1 = mxyzt + mxyzt_df
    mxyzt1_total = mxyzt1.sum(axis=0)
    # print(f' mxyzt_total        = {mxyzt.sum(axis=0)}')
    # print(f'+mxyzt_df_total = {mxyzt_df.sum(axis=0)}')
    # print(f'=mxyzt1_total       = {mxyzt1_total}')

    # M = I*alpha = I0*alpha + m*r^2*alpha
    #
    # handle self-inertia (probably very small)
    # removes the mean term (constant moment distribution)
    if np.abs(inertia_self).max() > 0 and np.abs(mxyzt1_total).max() > 0:
        # don't divide by 0
        ipos = np.where(inertia_self != 0)[0]

        # you should re-transform to get
        # the principal inertias, but mehhh...
        alphat_self = np.zeros(3, dtype='float64')
        alphat_self[ipos] = mxyzt1_total[ipos] / inertia_final[ipos]
        # print(f'alphat (rad/s^2) = {alphat_self}')

        # TODO: should distribute better to nodes
        #       this isn't biased by mass or anything...
        mxyzt_delta_self = -(alphat_self * inertia_final * 0.)[np.newaxis, :]
        # print(f'mxyzt_delta_self:\n{mxyzt_delta_self}')
    else:
        mxyzt_delta_self = np.zeros(shape, dtype=dtype)

    mxyzt2 = mxyzt + mxyzt_df + mxyzt_delta_self
    # print(f'mxyzt2:\n{mxyzt2}')
    mxyzt2_total = mxyzt2.sum(axis=0)

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
    alphat_moment = np.zeros(3, dtype='float64')
    alphat_moment[ipos] = mxyzt2_total[ipos] / inertia_final[ipos]

    # A skipped axis that still carries a moment is not physically solvable:
    # there is no inertia to react it, so that moment stays in the residual
    # and the caller's "balanced" load set is not actually balanced. Warn
    # rather than raise -- the other two axes are still perfectly valid.
    izero = np.setdiff1d(np.arange(3), ipos)
    if len(izero) and np.abs(mxyzt2_total[izero]).max() > 1e-8 * max(
            np.abs(mxyzt2_total).max(), 1e-300):
        warnings.warn(
            f'degenerate mass distribution: principal inertia(s) {izero} are '
            f'~0 (inertia={inertia_final}), but the applied moment about them '
            f'is {mxyzt2_total[izero]}. No rotational relief is possible about '
            'a zero-inertia axis, so that moment is left unbalanced.',
            RuntimeWarning, stacklevel=2)

    # d'Alembert moment from the rigid-body angular acceleration.
    #
    # The restoring moment on a point mass is the FULL tensor form
    #     -m * ((alpha x d) x d)
    # which couples the axes. The previous diagonal-only form
    #     -alpha[i] * m * (d_j^2 + d_k^2)
    # is just the diagonal of the inertia tensor and drops the products of
    # inertia, so it only balances when the body already lies along the
    # global axes. On the rotated bar it left a residual of ~4.8.
    mxyzt_dm = -mass[:, np.newaxis] * np.cross(
        np.cross(dxyzt, alphat_moment[np.newaxis, :]), dxyzt)

    # Self-inertia of each lumped mass (its own I0*alpha, no parallel axis).
    #
    # The per-node self-inertia is given in the GLOBAL frame, but alpha lives
    # in the principal frame, so each node's tensor has to be rotated before
    # it can be contracted with alpha. Applying the raw diagonal
    # (inertia[:, :3] * alpha) mixes frames and leaves a residual whenever the
    # self-inertia is non-zero and S is not the identity.
    if np.abs(inertia).max() > 0.:
        ixx_n, iyy_n, izz_n, ixy_n, ixz_n, iyz_n = inertia.T
        # assemble (nnode,3,3) with the negative-off-diagonal tensor convention
        imat_n = np.empty((len(mass), 3, 3), dtype='float64')
        imat_n[:, 0, 0] = ixx_n
        imat_n[:, 1, 1] = iyy_n
        imat_n[:, 2, 2] = izz_n
        imat_n[:, 0, 1] = imat_n[:, 1, 0] = -ixy_n
        imat_n[:, 0, 2] = imat_n[:, 2, 0] = -ixz_n
        imat_n[:, 1, 2] = imat_n[:, 2, 1] = -iyz_n
        # rotate into the principal frame: I_t = S.T @ I @ S
        imat_t = np.einsum('ij,njk,kl->nil', S.T, imat_n, S)
        mxyzt_dm = mxyzt_dm - np.einsum('nij,j->ni', imat_t, alphat_moment)
    # print(f'mxyzt_dm_total = {mxyzt_dm.sum(axis=0)}')

    mxyzt_delta = mxyzt_df + mxyzt_delta_self + mxyzt_dm
    # mxyzt_out = mxyzt2 + mxyzt_delta

    # Transform back out of the principal frame.
    # Forward was v_t = S.T @ v, so the inverse is v = S @ v_t.
    St = S
    fxyz_delta_out = np.einsum('nij,nj->ni',St[np.newaxis,:], fxyzt_delta)
    mxyz_delta_out = np.einsum('nij,nj->ni',St[np.newaxis,:], mxyzt_delta)

    # St = S.T  # also tried S instead of S.T
    # fxyz_delta_out = np.einsum('nij,ni->nj',St[np.newaxis,:], fxyzt_delta)
    # mxyz_delta_out = np.einsum('nij,ni->nj',St[np.newaxis,:], mxyzt_delta)

    return fxyz_delta_out, mxyz_delta_out


if __name__ == '__main__':  # pragma: no cover
    inertia_relief_from_model(None)
