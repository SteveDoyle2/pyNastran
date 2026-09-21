"""
Forces and moments on lumped masses due to a rigid-body acceleration field.

Purpose
-------
Given a set of point masses and a prescribed rigid-body motion state
(translational ``accel``, angular acceleration ``alpha``, and angular rate
``omega`` about a reference point), compute the d'Alembert load carried by each
mass::

    r_i = cg_i - ref_xyz
    a_i = accel + alpha x r_i + omega x (omega x r_i)
    F_i = -m_i * a_i
    M_i = -([I_i] alpha + omega x ([I_i] omega))

This is the forward problem.  ``inertia_relief2.inertia_relief`` solves the
inverse one: given unbalanced applied loads, find the accelerations.  The two
are consistent for ``omega = 0``; that routine has no rate terms, so it cannot
reproduce a steady-roll or steady-turn load.

Interface
---------
Deliberately array-only.  A lumped mass is a mass, a cg, and a self-inertia;
nothing here knows about ``BDF``, CONM2, or any other card.  The caller
assembles the arrays.  To get loads on specific nodes, pass only those rows.

Sign convention
---------------
Returned loads are **d'Alembert (inertial) loads** -- the load the accelerating
mass applies to the structure, ``F = -m*a``.  This matches the sign returned by
``inertia_relief2.inertia_relief``.  Negate for the reaction that the structure
applies to the mass.

The packed ``(nnode,6)`` self-inertia is ``[Ixx, Iyy, Izz, Ixy, Ixz, Iyz]`` with
the products of inertia stored as the **positive** integrals ``Ixy = sum(m dx dy)``
(the pyNastran/CONM2 convention), so the physical tensor negates the
off-diagonals::

    [ Ixx  -Ixy  -Ixz]
    [-Ixy   Iyy  -Iyz]
    [-Ixz  -Iyz   Izz]

Rate terms
----------
Steady rolls and steady turns are supported: pass ``omega``.  Two distinct
velocity-dependent terms appear, and they are *not* the same thing as an engine
gyroscopic load.

- **Centrifugal**, ``omega x (omega x r_i)``, in the force.  This is the one
  that matters for a roll or a turn -- it throws every off-axis mass outward
  from the rotation axis, and it is quadratic in rate, so it dominates at high
  roll rate even though ``alpha`` is zero in the steady case.
- **Rate-coupling**, ``omega x ([I_i] omega)``, in the moment.  This is the
  body's own inertia resisting a reorientation of its rate vector.  It is zero
  whenever ``omega`` is along a principal axis of ``[I_i]``, which is why a
  pure roll of an axisymmetric store shows nothing here but a coupled
  roll-yaw does.

**Not included: engine/rotor gyroscopic loads.**  A spinning rotor carries an
angular momentum ``[I_rotor] omega_spin`` of its own, about its own shaft, that
is not represented by a point mass's ``[I_i]``.  The resulting
``omega_body x h_rotor`` couple is a separate load and is deliberately out of
scope.  If it is ever wanted it belongs as an explicit per-node angular-momentum
input, not as a side effect of the mass properties.

Assumptions
-----------
- **Rigid body.**  Only the six rigid-body freedoms; structural flexibility
  does not feed back into the load distribution.
- **Instantaneous state.**  ``omega`` and ``alpha`` describe the motion at one
  instant.  Nothing here integrates, so a "steady" roll or turn means the state
  is steady at the evaluated instant; the caller supplies the trimmed rates.
- **No Coriolis.**  ``2 omega x v_rel`` vanishes because the masses are rigidly
  attached -- there is no velocity relative to the rotating frame.  This is a
  consequence of the rigid-body assumption, not an extra approximation.
- **Small angles.**  ``cg`` and ``inertia`` are evaluated once in the reference
  configuration and held fixed.
- ``inertia`` is the **self**-inertia of each lumped mass about its own cg.
  The parallel-axis ``m*r^2`` contribution is already carried by the
  ``alpha x r`` term in the force; adding it to ``inertia`` double-counts it.
- All arrays are in one consistent (basic) frame.  No coordinate transforms
  are performed here.

"""
from __future__ import annotations

import numpy as np


def inertia_loads(mass: np.ndarray,
                  cg: np.ndarray,
                  inertia: np.ndarray,
                  accel: np.ndarray,
                  alpha: np.ndarray,
                  ref_xyz: np.ndarray,
                  omega: np.ndarray | None = None,
                  xyz: np.ndarray | None = None,
                  wtmass: float = 1.0) -> tuple[np.ndarray, np.ndarray]:
    """
    Computes the per-mass inertial force and moment for a rigid-body motion
    state.

    Parameters
    ----------
    mass : (nnode,) float ndarray
        mass (or weight, if ``wtmass`` is used) of each lumped mass
    cg : (nnode,3) float ndarray
        location of each mass's center of gravity.  For a CONM2 this is the
        grid location *plus* the X1-X3 offset, not the grid location.
    inertia : (nnode,6) float ndarray
        self-inertia about each mass's own cg, packed as
        ``[Ixx, Iyy, Izz, Ixy, Ixz, Iyz]`` with products of inertia as positive
        integrals.  Pass zeros for masses with no rotary inertia.
        Units of M*L^2 (or W*L^2, scaled by ``wtmass``).
    accel : (3,) float ndarray
        translational acceleration **of the reference point**, not of the cg.
        Units of L/T^2.  For a gravity/load-factor case pass e.g.
        ``[0, 0, -n_z*g]``.
    alpha : (3,) float ndarray
        angular acceleration about ``ref_xyz``.  Units of 1/T^2.
    ref_xyz : (3,) float ndarray
        the point that ``accel`` is measured at and that ``alpha`` rotates
        about.  Required, and not defaulted to the cg on purpose: an
        accelerometer-measured ``accel`` is referenced to the instrument
        station, and silently re-referencing it to the cg would change the
        answer by ``alpha x (cg - station)``.
    omega : (3,) float ndarray; default=None -> zeros
        angular rate about ``ref_xyz``.  Units of 1/T.  Supply this for a
        steady roll (``omega = [p, 0, 0]``, ``alpha = 0``) or a steady turn.
        Adds the centrifugal force and rate-coupling moment described in the
        module docstring.  Leave as None for a pure maneuver-acceleration case.
    xyz : (nnode,3) float ndarray; default=None -> cg
        where to report the moment.  ``None`` reports at each mass's own cg
        (no offset couple).  Pass the grid locations to report at the grids,
        which adds the ``(cg - xyz) x F`` couple from the CONM2 offset.
    wtmass : float; default=1.0
        weight-to-mass conversion (``PARAM,WTMASS``).  Applied to both ``mass``
        and ``inertia``, since it scales the whole 6x6 mass matrix.  Use 1/g
        when passing weight.

    Returns
    -------
    force : (nnode,3) float ndarray
        inertial force on each mass
    moment : (nnode,3) float ndarray
        inertial moment on each mass, about ``xyz`` (or the cg if ``xyz`` is
        None).  Contains the self-inertia and rate-coupling terms plus, if
        ``xyz`` was given, the offset couple.  It does **not** include
        ``(cg - ref_xyz) x F``; that couple belongs to the resultant about
        ``ref_xyz``, not to the local load.  See ``resultant``.

    Examples
    --------
    A 2 g down load on a single 10-unit mass, no rotation::

        force, moment = inertia_loads(
            np.array([10.]), np.array([[1., 0., 0.]]), np.zeros((1, 6)),
            np.array([0., 0., -2.]), np.zeros(3), np.zeros(3))
        # force -> [[0., 0., 20.]]

    A steady roll at rate p about the x axis, no other motion::

        force, moment = inertia_loads(
            mass, cg, inertia, np.zeros(3), np.zeros(3), ref_xyz,
            omega=np.array([p, 0., 0.]))

    """
    mass = np.asarray(mass, dtype='float64')
    cg = np.asarray(cg, dtype='float64')
    inertia = np.asarray(inertia, dtype='float64')
    accel = np.asarray(accel, dtype='float64').ravel()
    alpha = np.asarray(alpha, dtype='float64').ravel()
    ref_xyz = np.asarray(ref_xyz, dtype='float64').ravel()
    omega = (np.zeros(3, dtype='float64') if omega is None else
             np.asarray(omega, dtype='float64').ravel())

    nnode = len(mass)
    if mass.ndim != 1:
        raise ValueError(f'mass must be (nnode,); got {mass.shape}')
    if cg.shape != (nnode, 3):
        raise ValueError(f'cg must be ({nnode},3); got {cg.shape}')
    if inertia.shape != (nnode, 6):
        raise ValueError(f'inertia must be ({nnode},6); got {inertia.shape}')
    for name, vector in (('accel', accel), ('alpha', alpha),
                         ('ref_xyz', ref_xyz), ('omega', omega)):
        if vector.shape != (3,):
            raise ValueError(f'{name} must be (3,); got {vector.shape}')

    # scale to mass units; wtmass multiplies the entire 6x6 mass matrix
    massi = mass * wtmass
    inertiai = inertia * wtmass

    radius = cg - ref_xyz[np.newaxis, :]
    # old
    # a_i = accel + alpha x r_i + omega x (omega x r_i)
    #radius = cg - ref_xyz[np.newaxis, :]
    #accel_node = accel[np.newaxis, :] + np.cross(alpha[np.newaxis, :], radius)
    # centrifugal; zero when omega is zero, so no branch is needed
    #accel_node += np.cross(omega[np.newaxis, :],
    #                       np.cross(omega[np.newaxis, :], radius))
    #
    # d'Alembert force
    #force = -massi[:, np.newaxis] * accel_node

    # a_i = accel + alpha x r_i + omega x (omega x r_i)
    #     = accel + (skew(alpha) + skew(omega) skew(omega)) r_i
    # Both rate terms are linear in r_i, so they fold into one 3x3 and the
    # whole acceleration field is a single (nnode,3) @ (3,3) matmul.  The
    # centrifugal block is zero when omega is zero, so no branch is needed.
    skew_omega = _skew(omega)
    amat = _skew(alpha) + skew_omega @ skew_omega

    # d'Alembert force, built in place to avoid (nnode,3) temporaries
    force = radius @ amat.T
    force += accel[np.newaxis, :]
    force *= -massi[:, np.newaxis]

    # M_i = -([I_i] alpha + omega x ([I_i] omega))
    # [I_i] v is linear in the packed 6-vector, so this is a (nnode,6) @ (6,3)
    # matmul and the (nnode,3,3) tensor stack never has to be built.
    # The second term is the rate coupling; it vanishes when omega is along a
    # principal axis of [I_i].
    #
    # old
    # tensor stack carries the negated off-diagonals
    #imat = _inertia_tensor(inertiai)
    # rate coupling; vanishes when omega is along a principal axis of [I_i]
    #moment = -np.einsum('nij,j->ni', imat, alpha)
    #moment -= np.cross(omega[np.newaxis, :],
    #                   np.einsum('nij,j->ni', imat, omega))
    #
    bmat = _packed_inertia_dot(alpha) + skew_omega @ _packed_inertia_dot(omega)
    moment = inertiai @ -bmat.T

    if xyz is not None:
        xyz = np.asarray(xyz, dtype='float64')
        if xyz.shape != (nnode, 3):
            raise ValueError(f'xyz must be ({nnode},3); got {xyz.shape}')
        # move the force from the cg to the reporting point
        moment = moment + np.cross(cg - xyz, force)

    return force, moment


def resultant(force: np.ndarray,
              moment: np.ndarray,
              xyz: np.ndarray,
              ref_xyz: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """
    Sums a set of per-node loads into a single force and moment about
    ``ref_xyz``.

    Parameters
    ----------
    force, moment : (nnode,3) float ndarray
        the per-node loads, as returned by ``inertia_loads``
    xyz : (nnode,3) float ndarray
        where each load acts -- the same array passed to ``inertia_loads`` as
        ``xyz``, or ``cg`` if that argument was left as None
    ref_xyz : (3,) float ndarray
        the point to sum moments about

    Returns
    -------
    force_total : (3,) float ndarray
    moment_total : (3,) float ndarray
        includes the ``r x F`` transport couple

    """
    force = np.asarray(force, dtype='float64')
    moment = np.asarray(moment, dtype='float64')
    xyz = np.asarray(xyz, dtype='float64')
    ref_xyz = np.asarray(ref_xyz, dtype='float64').ravel()

    force_total = force.sum(axis=0)
    radius = xyz - ref_xyz[np.newaxis, :]

    # sum(r x F) componentwise as six dot products.  np.cross would build an
    # (nnode,3) temporary only to reduce it away immediately.
    #
    # old
    # moment_total = moment.sum(axis=0) + np.cross(radius, force).sum(axis=0)
    #
    rx, ry, rz = radius[:, 0], radius[:, 1], radius[:, 2]
    fx, fy, fz = force[:, 0], force[:, 1], force[:, 2]
    moment_total = moment.sum(axis=0) + np.array([
        ry @ fz - rz @ fy,
        rz @ fx - rx @ fz,
        rx @ fy - ry @ fx,
    ])
    return force_total, moment_total


def expected_resultant(mass: np.ndarray,
                       cg: np.ndarray,
                       inertia: np.ndarray,
                       accel: np.ndarray,
                       alpha: np.ndarray,
                       ref_xyz: np.ndarray,
                       omega: np.ndarray | None = None,
                       wtmass: float = 1.0) -> tuple[np.ndarray, np.ndarray]:
    """
    Closed-form rigid-body resultant about ``ref_xyz``, for checking
    ``inertia_loads`` independently of how the load was distributed::

        d = cg_total - ref_xyz
        sum(F) = -m_total * (accel + alpha x d + omega x (omega x d))
        sum(M) = -m_total * d x accel
                 - [I_total] alpha
                 - omega x ([I_total] omega)

    where ``[I_total]`` is taken about ``ref_xyz`` (self-inertia plus the
    parallel-axis terms).  The two ``omega`` terms follow from the vector
    identities ``sum m r x (omega x (omega x r)) = omega x ([I] omega)`` and
    ``sum m r x (alpha x r) = [I] alpha``, with ``[I]`` the parallel-axis
    tensor -- verified numerically against a finite-difference of the angular
    momentum, not just asserted.

    This is not a re-derivation of the same algebra: it sums the rigid-body
    equations of motion directly, without ever forming a per-node load, so it
    catches a mis-distributed or double-counted term.

    Returns
    -------
    force_total : (3,) float ndarray
    moment_total : (3,) float ndarray

    """
    mass = np.asarray(mass, dtype='float64')
    cg = np.asarray(cg, dtype='float64')
    inertia = np.asarray(inertia, dtype='float64')
    accel = np.asarray(accel, dtype='float64').ravel()
    alpha = np.asarray(alpha, dtype='float64').ravel()
    ref_xyz = np.asarray(ref_xyz, dtype='float64').ravel()
    omega = (np.zeros(3, dtype='float64') if omega is None else
             np.asarray(omega, dtype='float64').ravel())

    massi = mass * wtmass
    inertiai = inertia * wtmass

    mass_total = massi.sum()
    radius = cg - ref_xyz[np.newaxis, :]

    if mass_total == 0.0:
        cg_radius = np.zeros(3, dtype='float64')
    else:
        cg_radius = (massi[:, np.newaxis] * radius).sum(axis=0) / mass_total

    force_total = -mass_total * (
        accel
        + np.cross(alpha, cg_radius)
        + np.cross(omega, np.cross(omega, cg_radius)))

    # [I] about ref = sum( I_self + m ((r.r) I3 - r r^T) )
    # The self term sums in the packed (6,) form first -- summing 6 columns is
    # far cheaper than building and reducing an (nnode,3,3) stack -- and the
    # parallel-axis term is a weighted Gram matrix, i.e. one (3,nnode)@(nnode,3)
    #
    # old
    #imat = _inertia_tensor(inertiai).sum(axis=0)
    #rdotr = np.einsum('ni,ni->n', radius, radius)
    #
    # matmul rather than a three-index einsum.
    inertia_self = inertiai.sum(axis=0)
    ixx, iyy, izz, ixy, ixz, iyz = inertia_self
    imat = np.array([
        [ixx, -ixy, -ixz],
        [-ixy, iyy, -iyz],
        [-ixz, -iyz, izz],
    ], dtype='float64')

    mass_radius = massi[:, np.newaxis] * radius
    eye = np.eye(3, dtype='float64')

    # old
    #parallel = np.einsum('n,n,ij->ij', massi, rdotr, eye)
    #parallel -= np.einsum('n,ni,nj->ij', massi, radius, radius)

    parallel = eye * (mass_radius * radius).sum()
    parallel -= mass_radius.T @ radius
    imat = imat + parallel

    moment_total = (-mass_total * np.cross(cg_radius, accel)
                    - imat @ alpha
                    - np.cross(omega, imat @ omega))
    return force_total, moment_total


def _skew(vector: np.ndarray) -> np.ndarray:
    """
    The 3x3 skew-symmetric matrix with ``_skew(a) @ b == cross(a, b)``.

    Lets a cross product against many vectors be written as one matmul, which
    is several times faster than ``np.cross`` on large arrays because it hands
    the work to BLAS instead of building intermediates.

    """
    x, y, z = vector
    return np.array([
        [0.0, -z, y],
        [z, 0.0, -x],
        [-y, x, 0.0],
    ], dtype='float64')


def _packed_inertia_dot(vector: np.ndarray) -> np.ndarray:
    """
    The 3x6 matrix ``B`` with ``B @ inertia_packed == [I] @ vector``.

    ``[I] v`` is linear in both ``[I]`` and ``v``, so contracting the packed
    6-term inertia with a fixed vector can be written as a small constant
    matrix.  That turns a per-node ``(nnode,3,3)`` tensor build plus an einsum
    into one ``(nnode,6) @ (6,3)`` matmul.

    The packed order is ``[Ixx, Iyy, Izz, Ixy, Ixz, Iyz]`` with the products of
    inertia stored as positive integrals, so the columns for ``Ixy``, ``Ixz``,
    and ``Iyz`` carry the minus signs of the physical tensor.

    """
    x, y, z = vector
    #                Ixx  Iyy  Izz  Ixy  Ixz  Iyz
    return np.array([
        [x, 0.0, 0.0, -y, -z, 0.0],
        [0.0, y, 0.0, -x, 0.0, -z],
        [0.0, 0.0, z, 0.0, -x, -y],
    ], dtype='float64')


def _inertia_tensor(inertia: np.ndarray) -> np.ndarray:
    """
    Builds an (nnode,3,3) tensor stack from the packed (nnode,6) array.

    The packed form stores products of inertia as positive integrals, so the
    off-diagonals are negated here.  See the module docstring.

    .. note::
        The hot paths in ``inertia_loads`` and ``expected_resultant`` avoid
        this: materializing ``(nnode,3,3)`` costs 9 floats per node where 6
        suffice, and it is the single most expensive step at large ``nnode``.
        Kept because it is the clearest statement of the sign convention, and
        the tests check the fast paths against it.

    """
    ixx = inertia[:, 0]
    iyy = inertia[:, 1]
    izz = inertia[:, 2]
    ixy = inertia[:, 3]
    ixz = inertia[:, 4]
    iyz = inertia[:, 5]

    nnode = len(inertia)
    imat = np.empty((nnode, 3, 3), dtype='float64')
    imat[:, 0, 0] = ixx
    imat[:, 1, 1] = iyy
    imat[:, 2, 2] = izz
    imat[:, 0, 1] = imat[:, 1, 0] = -ixy
    imat[:, 0, 2] = imat[:, 2, 0] = -ixz
    imat[:, 1, 2] = imat[:, 2, 1] = -iyz
    return imat
