"""
Inertia relief on a 1-D load distribution, with plotting (station formulation).

Purpose
-------
A beam/station-level tool, separate from ``inertia_relief2.py``.  Where that
module works on a 3-D cloud of lumped masses and returns load increments for a
Nastran model, this one takes a **running load distribution along a station
axis** (the classic shear/moment-diagram layout for a fuselage or wing) and
produces the balanced distribution plus matplotlib diagnostics.

Scope -- axis-aligned models only
--------------------------------
**This module is only valid for models whose station axis and load directions
are aligned with the global axes.**  Moments are formed by an elementwise
multiply (``force * dxyz_cg``) rather than a cross product, so the coupling
between off-axis components is simply absent.  That happens to give the right
answer for the axis-aligned station layout this tool was written for -- a
fuselage along x loaded in z, say -- because each arm/force pair lies on
distinct single axes that are handled separately.  It is wrong for a rotated
or genuinely 3-D load; see defect 1 below for a worked example.

If the model is not axis-aligned, use ``inertia_relief2.py`` instead: it works
on a 3-D cloud of lumped masses, transforms into the principal inertia frame,
and is covered by ``test_inrel.py``.

Two parallel implementations are provided:

- ``inertia_relief1`` -- scalar/1-D.  ``x`` is a station array and the loads
  are scalars along it.  This is the more complete of the two: it returns the
  eight load arrays and drives the plotting.
- ``inertia_relief3`` -- vector/3-D.  ``xyz`` is an ``(n,3)`` station array
  and the loads are ``(n,3)``.  Solves the coupled ``M = [I] alpha`` with a
  full 3x3 ``np.linalg.solve`` instead of a scalar divide.

Method
------
Both follow the same two-pass sequence:

1. **Linear relief.**  ``a = sum(F) / sum(m)``, then apply ``-m_i * a`` at each
   station.  This zeroes the net force.
2. **Angular relief.**  ``alpha = [I]^-1 sum(M)`` about the cg, then apply
   ``-I_i * alpha`` at each station.  This zeroes the net moment.

The angular increment is converted back to an equivalent *force* by dividing
by the moment arm (``inertial_angular_force = inertial_angular_moment / dxyz_cg``)
so that it can be superposed on the running load.

Assumptions
-----------
- **Rigid body, quasi-static.**  Same as ``inertia_relief2``: no centrifugal,
  Coriolis, or gyroscopic terms; the body is released from rest.
- **Station-ordered data.**  The integrated plots assume ``x`` is monotonic.
  ``fintegrate`` inspects only ``x[0]`` vs ``x[-1]`` to pick a sign, so a
  non-monotonic station array integrates incorrectly without warning.
- **Mass and load are already discretized** into per-station lumps.  These are
  summed directly, so the result depends on the station spacing the caller
  chose; the integrated plots use trapezoidal integration but the balance
  itself does not.
- **Inertia defaults to ``m*d^2``** about the cg when ``inertia_cg`` is not
  supplied -- a point-mass idealization that ignores each station's own
  self-inertia.

Limitations and known defects (verified by execution, 2026-09-21)
-----------------------------------------------------------------
This is ``dev/`` scratch code.  There is no test file covering this module
(``test_inrel.py`` imports only from ``inertia_relief2``), and several of the
issues below are hard errors rather than accuracy concerns:

1. **Moments are formed by elementwise multiply, not a cross product.**
   ``exterior_moment = exterior_force * dxyz_cg`` (and the ``inertial_*``
   equivalents) multiply component-by-component.  The correct rigid-body
   moment is ``d x F``.  These differ entirely: for ``F = [0,0,10]`` at
   ``d = [5,0,0]`` the code yields ``[0,0,0]`` -- the moment is silently lost --
   where ``cross(d,F)`` gives ``[0,-50,0]``.  ``inertia_relief3`` is therefore
   wrong for any genuinely 3-D load; it only coincides with the right answer
   in the 1-D case where the arm and the force are along different single axes
   handled separately.
2. **The inertia tensor is assembled with the wrong sign.**  ``Imat`` (and the
   per-station ``Imati``) use ``+Ixy, +Ixz, +Iyz`` off-diagonal, but the
   physical tensor ``sum m ((r.r) I3 - r r.T)`` requires them negative -- which
   is what ``inertia_relief2`` correctly does, and what the Nastran CONM2
   convention states.  On a random 3-D mass cloud this changes ``alpha`` by up
   to ~95%.
3. **(Not a defect -- by design.)**  The ``tol = 1e-16`` added to the diagonal
   of ``Imat``/``Imati`` exists to handle a degenerate mass distribution where
   a principal inertia is *exactly* zero -- e.g. a collinear (1-D) set of
   masses along x has ``Ixx = 0``, making the tensor singular and
   ``np.linalg.solve`` raise ``LinAlgError``.  The regularizer is deliberately
   chosen far below roundoff so it cannot contaminate a well-conditioned
   solve (verified: it changes ``alpha`` by 0.0 relative on a random 3-D mass
   cloud), while still making the degenerate case invertible.

   The behavior it produces is the physically sensible one.  For a
   *consistent* load -- no applied moment about the degenerate axis -- the
   degenerate component evaluates to ``0/tol = 0``, i.e. no rotational relief
   about an axis with no inertia, which is correct.  For an *inconsistent*
   load -- a non-zero moment about an axis with zero inertia, which has no
   physical solution -- it yields a very large ``alpha`` (``5/1e-16 = 5e16``)
   that makes the ill-posed load case obvious in the output rather than
   silently absorbing it.

   The one improvement worth making is diagnostic, not numerical: the large
   ``alpha`` is currently unlabelled, so a caller has to recognize the
   magnitude themselves.  Emitting a warning when
   ``abs(M[i]) > 0`` and ``I[i] <= tol`` would name the degenerate axis
   explicitly.
4. **``g`` is accepted but never used.**  Both functions do
   ``if g is not None: accel = g``, and ``accel`` is then never read; the code
   always proceeds with ``exterior_accel = sum(F)/sum(m)``.  Passing ``g`` has
   no effect on the result, contrary to the docstrings.
5. **``_get_inertia3`` raises when ``inertia_cg`` is supplied.**  ``dxyz_cg`` is
   only assigned in the ``else`` branch, so the user-supplied-inertia path exits
   with ``UnboundLocalError: local variable 'dxyz_cg' referenced before
   assignment``.  Confirmed by execution.
6. **``inertia_relief3`` returns ``None``.**  It prints a summary and exits;
   unlike ``inertia_relief1`` it never returns the balanced loads, and it
   never plots, so ``show``/``plot_differential``/``plot_integrated``/
   ``locations``/``case`` are all accepted and ignored.
7. **``plot_inertia`` is 1-D only.**  Its ``accel``/``alpha`` are formatted with
   ``f'{accel:.6f}'``, which raises ``TypeError`` on the ``(3,)`` arrays that
   ``inertia_relief3`` holds.  It is only callable from ``inertia_relief1``.
8. ``inertia_relief1``'s signature documents ``x_cg`` but ``test_inertia1``
   calls it with ``xyz_cg=``, which would raise ``TypeError``.  The in-file
   test functions are not run by the test suite and are stale.

"""
from typing import Optional
import numpy as np
from cpylog import SimpleLogger
from scipy.integrate import cumulative_trapezoid
import matplotlib.pyplot as plt

def inertia_relief3(xyz: np.ndarray,
                    exterior_force: np.ndarray,
                    mass: np.ndarray,
                    g: Optional[float]=None,
                    xyz_cg: Optional[float]=None,
                    inertia_cg: Optional[np.ndarray]=None,
                    include_linear_inertia: bool=True,
                    include_angular_inertia: bool=True,
                    mass_units: str='slinch',
                    length_units: str='in',
                    force_units: str='lbf',
                    accel_units: Optional[str]=None,
                    inertia_units: Optional[str]=None,
                    case: str='',
                    show: bool=True,
                    locations: dict[str, float]=None,
                    plot_differential: bool=True,
                    plot_integrated: bool=True):
    """
    3-D (vector) inertia relief on a station-based load distribution.

    Solves the coupled ``M = [I] alpha`` with a full 3x3 solve, unlike the
    scalar ``inertia_relief1``.

    .. warning::
        **This function is not trustworthy in its current state.**  It builds
        moments with an elementwise multiply instead of a cross product,
        assembles the inertia tensor with the wrong off-diagonal sign, and
        returns ``None``.  See defects 1, 2, 3, 5 and 6 in the module
        docstring.  Prefer ``inertia_relief2.inertia_relief`` for 3-D work.

    Parameters
    ----------
    xyz : (n, 3) float array
        the xyz station
    xyz_cg : (3, ) float array; default=None -> mass-weighted centroid
        the center of gravity of the vehicle
    exterior_force : (n, 3) float array
        differential (per-station) applied force
    mass : (n, ) float array
        differential mass
        mass has units of slinch, slug, kg, lbm
    inertia_cg : (n, 6) float array; default=None -> mass * d^2
        differential inertia about the cg, packed
        [Ixx, Iyy, Izz, Ixy, Ixz, Iyz]
        inertia has units of slinch-in^2, slug-ft^2, kg-m^2
        NOTE: supplying this currently raises UnboundLocalError (defect 5).
    g : float; default=None -> sum(F)/m
        the gravity to use for inertia relief
        probably should be positive (a=F/m) to apply inertia relief
        in the opposite direction
        slinch: 12*32.174=386.088
        NOTE: accepted but never applied (defect 4).
    include_linear_inertia / include_angular_inertia : bool; default=True
        toggles the -m*a and -I*alpha passes independently
    mass_units / length_units / force_units: str
        used for prints
    accel_units / inertia_units : Optional[str]
        automatically calculated; used for prints
    locations : dict[str, float] | None
        locations to identify key locations (e.g., nose, wing LE, tail)
        accepted but unused -- this function does not plot
    case / show / plot_differential / plot_integrated
        accepted but unused -- this function does not plot

    Returns
    -------
    None
        Unlike ``inertia_relief1``, the balanced loads are not returned; they
        are only summarized to the log (defect 6).

    """
    log = SimpleLogger()
    if locations is None:
        locations = {}
    assert isinstance(locations, dict), locations
    assert isinstance(xyz, np.ndarray), type(xyz)
    assert isinstance(mass, np.ndarray), type(mass)
    assert isinstance(exterior_force, np.ndarray), type(exterior_force)
    nx = len(mass)

    moment_units, accel_units, alpha_units, inertia_units = _get_units(
        mass_units, length_units, force_units,
        accel_units=accel_units,
        inertia_units=inertia_units)

    sum_mass = mass.sum()
    xyz_cgi, dxyz_cg, inertia_cgi, sum_inertia_cg = _get_mass_cg_inertia3(
        xyz, xyz_cg, mass, sum_mass,
        inertia_cg=inertia_cg)
    del inertia_cg, xyz_cg

    exterior_moment = exterior_force * dxyz_cg  # fintegrate(exterior_force, x=x)
    sum_exterior_force = exterior_force.sum(axis=0)
    sum_exterior_moment = exterior_moment.sum(axis=0)
    assert len(sum_exterior_force) == 3, sum_exterior_force
    assert len(sum_exterior_moment) == 3, sum_exterior_moment

    #  F = m * a
    #  a = F / m
    exterior_accel = sum_exterior_force / sum_mass
    assert len(exterior_accel) == 3, exterior_accel

    # Mcg = Icg * alpha_cg
    # rxx, ryy, rzz, rxy, rxz, ryz
    #print(inertia_cgi.shape)
    #--------------------------------------------
    ixx = inertia_cgi[:, 0]
    iyy = inertia_cgi[:, 1]
    izz = inertia_cgi[:, 2]
    ixy = inertia_cgi[:, 3]
    ixz = inertia_cgi[:, 4]
    iyz = inertia_cgi[:, 5]
    #ixx, iyy, izz, ixy, ixz, iyz = inertia_cgi
    Ixx, Iyy, Izz, Ixy, Ixz, Iyz = sum_inertia_cg

    # Adding a small amount of noise allows us to invert Imat.
    #
    # This handles the degenerate case where a principal inertia is exactly
    # zero -- e.g. a collinear set of masses along x gives Ixx = 0, so the
    # tensor is singular and np.linalg.solve raises LinAlgError.
    #
    # tol is deliberately far below roundoff (eps*I >> tol for any realistic
    # inertia), so it cannot perturb a well-conditioned solve. For a load with
    # no moment about the degenerate axis it gives 0/tol = 0 (no relief about
    # an axis with no inertia -- correct); for a moment about a zero-inertia
    # axis, which is not physically solvable, it gives a very large alpha that
    # makes the ill-posed case visible.
    #
    # DO NOT "fix" this by removing it -- the singular solve comes straight back.
    tol = 1e-16
    Imat = np.array([
        [Ixx, Ixy, Ixz],
        [Ixy, Iyy, Iyz],
        [Ixz, Iyz, Izz],
    ], dtype='float64') + tol * np.eye(3)

    # adding a small amount of noise allows us to invert Imati
    Imati = np.zeros((nx, 3, 3), dtype='float64')
    Imati[:, 0, 0] = ixx + tol
    Imati[:, 0, 1] = Imati[:, 1, 0] = ixy
    Imati[:, 0, 2] = Imati[:, 2, 0] = ixz
    Imati[:, 1, 1] = iyy + tol
    Imati[:, 1, 2] = Imati[:, 2, 1] = iyz
    Imati[:, 2, 2] = izz + tol
    #print(Imat)
    assert Imat.shape == (3, 3), Imat.shape
    #--------------------------------------------
    # [Mx]    [Ixx, Ixy, Ixz] [alpha_x]
    # [My]  = [Ixy, Iyy, Iyz] [alpha_y]
    # [Mz]    [Ixz, Iyz, Izz] [alpha_z]
    #
    # [Ixx, Ixy, Ixz]^-1 [Mx]   [alpha_x]
    # [Ixy, Iyy, Iyz]    [My] = [alpha_y]
    # [Ixz, Iyz, Izz]    [Mz]   [alpha_z]
    sum_exterior_moment = sum_exterior_moment.reshape(3, 1)
    exterior_alpha_cg = np.linalg.solve(Imat, sum_exterior_moment) #/ sum_inertia_cg
    sum_exterior_moment = sum_exterior_moment.flatten()
    exterior_alpha_cg = exterior_alpha_cg.flatten()
    assert len(exterior_alpha_cg) == 3, exterior_alpha_cg

    if g is not None:
        accel = g

    if include_linear_inertia:
        # Fi = mi * a
        inertial_linear_force = -mass[:, np.newaxis] * exterior_accel
        #inertial_linear_moment = fintegrate(inertial_linear_force, x=x)
        inertial_linear_moment = inertial_linear_force * dxyz_cg
    else:
        inertial_linear_force = np.zeros(exterior_force.shape, dtype=exterior_force.dtype)
        inertial_linear_moment = np.zeros(exterior_force.shape, dtype=exterior_force.dtype)
    assert inertial_linear_force.shape == (nx, 3), inertial_linear_force.shape
    assert inertial_linear_moment.shape == (nx, 3), inertial_linear_moment.shape

    #--------
    # Mx_cg            alpha_cg
    # My_cg = [Imat] @ alpha_cg
    # Mz_cg            alpha_cg
    linear_force = exterior_force + inertial_linear_force
    linear_moment = exterior_moment + inertial_linear_moment
    sum_linear_force = linear_force.sum(axis=0)
    sum_linear_moment = linear_moment.sum(axis=0)
    assert len(sum_linear_force) == 3, sum_linear_force
    assert len(sum_linear_moment) == 3, sum_linear_moment

    #linear_alpha_cg = sum_linear_moment / sum_inertia_cg
    linear_alpha_cg = np.linalg.solve(Imat, sum_linear_moment).flatten() #/ sum_inertia_cg

    inertial_angular_force = np.zeros(exterior_force.shape, dtype=exterior_force.dtype)
    log.info(f'dxyz_cg={type(dxyz_cg)} dxyz_cg.shape={str(dxyz_cg.shape)}')
    if include_angular_inertia:
        # Mcgi = -Icgi * alpha_cg
        icg_positive = np.abs(dxyz_cg) > 0
        #icg_positive = dxcg.abs() > 0
        assert linear_alpha_cg.shape == (3,), linear_alpha_cg.shape
        #inertial_angular_moment = -inertia_cgi * linear_alpha_cg
        inertial_angular_moment = -Imati @ linear_alpha_cg
        inertial_angular_force[icg_positive] = inertial_angular_moment[icg_positive] / dxyz_cg[icg_positive]
    else:
        inertial_angular_moment = np.zeros(exterior_force.shape, dtype=exterior_force.dtype)

    #-------------------------------------------------------------------------------------
    total_force = exterior_force + inertial_linear_force + inertial_angular_force
    total_moment = exterior_moment + inertial_linear_moment + inertial_angular_moment
    total_force_sum = total_force.sum(axis=0)
    total_moment_sum = total_moment.sum(axis=0)

    total_accel = total_force_sum / sum_mass
    #total_alpha_cg = total_moment.sum(axis=0) / sum_inertia_cg
    IImat = Imat[np.newaxis, :, :]
    #print(IImat.shape)
    #print(total_moment_sum.shape)
    total_alpha_cg = np.linalg.solve(Imat, total_moment_sum) #/ sum_inertia_cg

    assert len(total_force) == nx
    assert len(total_moment) == nx

    assert len(total_accel) == 3, total_accel
    assert len(total_alpha_cg) == 3, total_alpha_cg
    # print(f'  total_accel = {total_accel} {accel_units}')
    # print(f'  total_alpha_cg = {total_alpha_cg} {accel_units}')

    data_balanced = {
        f'xyz_cg ({length_units})': xyz_cgi,
        f'mass ({mass_units})': sum_mass,
        f'inertia_cg ({inertia_units})': sum_inertia_cg,
        f'exterior_force ({force_units})': sum_exterior_force,
        f'exterior_moment_cg ({moment_units})': sum_exterior_moment,

        f'total_force ({force_units})': total_force_sum,
        f'total_moment_cg ({moment_units})': total_moment_sum,

        f'exterior_accel ({accel_units})': exterior_accel,
        f'exterior_alpha_cg ({alpha_units})': exterior_alpha_cg,
        f'total_alpha_cg ({alpha_units})': total_alpha_cg,
        f'total_accel ({accel_units})': total_accel,
    }
    msg = _write_summary(data_balanced, indent='  ')
    log.info('Summary Balanced Loads')
    log.info('\n' + msg)

    print('done')


def inertia_relief1(x: np.ndarray,
                    exterior_force: np.ndarray,
                    mass: np.ndarray,
                    g: Optional[float]=None,
                    x_cg: Optional[float]=None,
                    inertia_cg: Optional[np.ndarray]=None,
                    include_linear_inertia: bool=True,
                    include_angular_inertia: bool=True,
                    mass_units: str='slinch',
                    length_units: str='in',
                    force_units: str='lbf',
                    accel_units: Optional[str]=None,
                    inertia_units: Optional[str]=None,
                    case: str='',
                    show: bool=True,
                    locations: dict[str, float]=None,
                    plot_differential: bool=True,
                    plot_integrated: bool=True):
    """
    1-D (scalar) inertia relief on a station-based running load.

    This is the more complete of the two implementations here: it returns the
    balanced load arrays and drives the shear/moment plotting.  Because
    everything is scalar along a single station axis, the cross-product and
    inertia-tensor sign defects that affect ``inertia_relief3`` do not arise;
    the moment arm ``dxcg`` is a signed scalar and ``inertia_cg`` is a scalar
    per station.

    .. note::
        **Axis-aligned models only.**  The scalar formulation presumes the
        station axis is a global axis and that the load acts transverse to it,
        which is what makes the elementwise moment ``force * dxcg`` equivalent
        to a cross product here.  There is no way to express a rotated or
        genuinely 3-D load in this signature, and no check that the caller's
        model satisfies the assumption.  For anything off-axis use
        ``inertia_relief2.inertia_relief``.

    Parameters
    ----------
    x : (n, ) float array
        the station coordinate; assumed monotonic (see ``fintegrate``)
    x_cg : float; default=None -> mass-weighted centroid
        the center of gravity of the vehicle
    exterior_force : (n, ) float array
        differential (per-station) applied force
    mass : (n, ) float array
        differential mass
        mass has units of slinch, slug, kg, lbm
    inertia_cg : (n, ) float array; default=None -> mass * dxcg^2
        differential inertia about the cg
        inertia has units of slinch-in^2, slug-ft^2, kg-m^2
    g : float; default=None -> sum(F)/m
        the gravity to use for inertia relief
        probably should be positive (a=F/m) to apply inertia relief
        in the opposite direction
        slinch: 12*32.174=386.088
        NOTE: accepted but never applied (defect 4 in the module docstring).
    include_linear_inertia / include_angular_inertia : bool; default=True
        toggles the -m*a and -I*alpha passes independently
    mass_units / length_units / force_units: str
        used for prints
    accel_units / inertia_units : Optional[str]
        automatically calculated; used for prints
    locations : dict[str, float] | None
        locations to identify key locations (e.g., nose, wing LE, tail)
        drawn as labelled vertical lines on every subplot
    case : str; default=''
        a label appended to the figure titles
    show : bool; default=True
        calls plt.show() before returning
    plot_differential / plot_integrated : bool; default=True
        emit the per-station and the cumulative (shear/moment diagram) figures

    Returns
    -------
    out : tuple of 8 (n,) float arrays
        (exterior_force, inertial_linear_force, inertial_angular_force,
         total_force,
         exterior_moment, inertial_linear_moment, inertial_angular_moment,
         total_moment)

    Notes
    -----
    - The angular pass converts the inertial moment back to an equivalent
      force by dividing by the moment arm, guarded by ``abs(dxcg) > 0`` so the
      station sitting exactly at the cg contributes no angular force.  A
      station very close to (but not exactly at) the cg divides by a near-zero
      arm and produces a large spurious force spike; the guard is an exact
      comparison, not a tolerance.
    - ``sum_inertia_cg`` is a plain scalar sum, so ``exterior_alpha_cg`` is a
      simple divide; there is no matrix solve and no principal-axis rotation
      in the 1-D formulation.

    """
    if locations is None:
        locations = {}
    assert isinstance(locations, dict), locations
    assert isinstance(x, np.ndarray), type(x)
    assert isinstance(mass, np.ndarray), type(mass)
    assert isinstance(exterior_force, np.ndarray), type(exterior_force)
    assert x.ndim == 1, x
    assert mass.ndim == 1, mass.shape

    moment_units, accel_units, alpha_units, inertia_units = _get_units(
        mass_units, length_units, force_units,
        accel_units=accel_units,
        inertia_units=inertia_units)

    sum_mass = mass.sum()
    xcgi, dxcg, inertia_cgi, sum_inertia_cg = _get_mass_cg_inertia1(
        x, x_cg, mass, sum_mass,
        inertia_cg=inertia_cg)

    exterior_moment = exterior_force * dxcg # fintegrate(exterior_force, x=x)
    sum_exterior_force = exterior_force.sum()
    sum_exterior_moment = exterior_moment.sum()

    #  F = m * a
    #  a = F / m
    exterior_accel = sum_exterior_force / sum_mass

    # Mcg = Icg * alpha_cg
    exterior_alpha_cg = sum_exterior_moment / sum_inertia_cg

    data = {
        f'xcg ({length_units})': xcgi,
        f'mass ({mass_units})': sum_mass,
        f'exterior_accel ({accel_units})': exterior_accel,
        f'g ({accel_units})': g,
        f'inertia_cg ({inertia_units})': sum_inertia_cg,
        f'exterior_alpha_cg ({alpha_units})': exterior_alpha_cg,
    }
    msg = _write_summary(data, indent='  ')

    word = ''
    if case:
        word = f' for {case}'
    print(f'Input Summary (about CG){word}')
    print(msg)

    if g is not None:
        accel = g

    if include_linear_inertia:
        # Fi = mi * a
        inertial_linear_force = -mass * exterior_accel
        #inertial_linear_moment = fintegrate(inertial_linear_force, x=x)
        inertial_linear_moment = inertial_linear_force * dxcg
    else:
        inertial_linear_force = np.zeros(exterior_force.shape, dtype=exterior_force.dtype)
        inertial_linear_moment = np.zeros(exterior_force.shape, dtype=exterior_force.dtype)

    # Mcg = Icg * alpha_cg
    linear_force = exterior_force + inertial_linear_force
    linear_moment = exterior_moment + inertial_linear_moment
    sum_linear_force = linear_force.sum()
    sum_linear_moment = linear_moment.sum()
    linear_alpha_cg = sum_linear_moment / sum_inertia_cg

    #---------------------------------------------------------------
    data = {
        #f'xcg ({length_units})': xcgi,
        #f'mass ({mass_units})': total_mass,
        #f'accel ({accel_units})': accel,
        f'F_linear ({force_units})': sum_linear_force,
        f'M_linear ({moment_units})': sum_linear_moment,
        f'inertia_cg ({inertia_units})': sum_inertia_cg,
        f'exterior_alpha_cg ({alpha_units})': exterior_alpha_cg,
        f'linear_alpha_cg ({alpha_units})': linear_alpha_cg,
    }
    msg = _write_summary(data, indent='  ')


    print('Summary after Linear Inertia Applied')
    print(msg)
    inertial_angular_force = np.zeros(exterior_force.shape, dtype=exterior_force.dtype)
    print('dxcg =', type(dxcg), dxcg)
    if include_angular_inertia:
        # Mcgi = -Icgi * alpha_cg
        icg_positive = np.abs(dxcg) > 0
        #icg_positive = dxcg.abs() > 0
        inertial_angular_moment = -inertia_cgi * linear_alpha_cg
        inertial_angular_force[icg_positive] = inertial_angular_moment[icg_positive] / dxcg[icg_positive]
    else:
        inertial_angular_moment = np.zeros(exterior_force.shape, dtype=exterior_force.dtype)

    #-------------------------------------------------------------------------------------
    total_force = exterior_force + inertial_linear_force + inertial_angular_force
    total_moment = exterior_moment + inertial_linear_moment + inertial_angular_moment
    total_accel = total_force.sum() / sum_mass
    total_alpha_cg = total_moment.sum() / sum_inertia_cg
    nx = len(x)
    assert len(total_force) == nx
    assert len(total_moment) == nx

    data_linear = {
        f'xcg ({length_units})': xcgi,
        f'mass ({mass_units})': sum_mass,
        f'exterior_force ({force_units})': sum_exterior_force,
        f'exterior_accel ({accel_units})': exterior_accel,
        f'total_accel ({accel_units})': total_accel,
    }
    msg = _write_summary(data_linear, indent='  ')
    print('Linear Summary')
    print(msg)

    data_angle = {
        f'xcg ({length_units})': xcgi,
        f'inertia_cg ({inertia_units})': sum_inertia_cg,
        f'exterior_moment_cg ({moment_units})': sum_exterior_moment,
        f'exterior_alpha_cg ({alpha_units})': exterior_alpha_cg,
        f'total_alpha_cg ({alpha_units})': total_alpha_cg,
    }
    msg = _write_summary(data_angle, indent='  ')
    print('Angular Summary')
    print(msg)

    data_balanced = {
        f'xcg ({length_units})': xcgi,
        f'mass ({mass_units})': sum_mass,
        f'inertia_cg ({inertia_units})': sum_inertia_cg,
        f'exterior_force ({force_units})': sum_exterior_force,
        f'exterior_moment_cg ({moment_units})': sum_exterior_moment,

        f'total_force ({force_units})': total_force.sum(),
        f'total_moment_cg ({moment_units})': total_moment.sum(),

        f'exterior_accel ({accel_units})': exterior_accel,
        f'exterior_alpha_cg ({alpha_units})': exterior_alpha_cg,
        f'total_alpha_cg ({alpha_units})': total_alpha_cg,
        f'total_accel ({accel_units})': total_accel,
    }
    msg = _write_summary(data_balanced, indent='  ')
    print('Summary Balanced Loads')
    print(msg)
    #print(f'  xcg  ({length_units}) = {xcgi:.3f}')
    #print(f'  inertia_cg ({inertia_units}) = {total_inertia_cg:.3f}')
    #print(f'  exterior_moment_cg ({force_units}) = {total_exterior_moment:.3f}')
    #print(f'  alpha_cg ({accel_units}) = {total_alpha_cg:.3f}')
    out = (
        exterior_force,  inertial_linear_force, inertial_angular_force, total_force,
        exterior_moment, inertial_linear_moment, inertial_angular_moment, total_moment,
    )
    if plot_differential:
        plot_inertia(x,
                     exterior_force, inertial_linear_force, inertial_angular_force, total_force,
                     exterior_moment, inertial_linear_moment, inertial_angular_moment, total_moment,
                     locations,
                     exterior_accel, exterior_alpha_cg,
                     accel_units=accel_units,
                     integrate=False,
                     case=case)
    if plot_integrated:
        plot_inertia(x,
                     exterior_force, inertial_linear_force, inertial_angular_force, total_force,
                     exterior_moment, inertial_linear_moment, inertial_angular_moment, total_moment,
                     locations,
                     exterior_accel, exterior_alpha_cg,
                     accel_units=accel_units,
                     integrate=True,
                     case=case)
    if show:
        plt.show()
    return out

def _get_units(mass_units: str, length_units: str,
               force_units: str,
               accel_units: Optional[str]=None,
               inertia_units: Optional[str]=None,
               ) -> tuple[str, str, str, str]:
    if inertia_units is None:
        inertia_units = f'{mass_units}-{length_units}^2'
    if accel_units is None:
        accel_units = f'{length_units}/s^2'
    moment_units = f'{length_units}-{force_units}'
    alpha_units = 'rad/s^2'
    return moment_units, accel_units, alpha_units, inertia_units

def _get_mass_cg_inertia1(x: np.ndarray,
                          x_cg: np.ndarray,
                          mass: np.ndarray,
                          sum_mass: float,
                          inertia_cg=None):
    assert x.ndim == 1, x
    assert mass.ndim == 1, mass.shape
    x_cgi = _get_x_cg(x, x_cg, mass, sum_mass)
    dx_cg, inertia_cgi = _get_inertia1(x, x_cgi, mass, inertia_cg)
    del inertia_cg
    sum_inertia_cg = inertia_cgi.sum()
    return x_cgi, dx_cg, inertia_cgi, sum_inertia_cg

def _get_mass_cg_inertia3(xyz: np.ndarray,
                          xyz_cg: np.ndarray,
                          mass: np.ndarray,
                          sum_mass: float,
                          inertia_cg=None):
    assert xyz.ndim == 2, xyz
    assert mass.ndim == 1, mass.shape
    xyz_cgi = _get_xyz_cg(xyz, xyz_cg, mass, sum_mass)
    assert len(xyz_cgi) == 3, xyz_cgi
    dxyz_cg, inertia_cgi = _get_inertia3(xyz, xyz_cgi, mass, inertia_cg)

    sum_inertia_cg = inertia_cgi.sum(axis=0)
    Ixx, Iyy, Izz, Ixy, Ixz, Iyz = sum_inertia_cg
    assert len(sum_inertia_cg) == 6, sum_inertia_cg
    return xyz_cgi, dxyz_cg, inertia_cgi, sum_inertia_cg


def _get_x_cg(x: np.ndarray,
              x_cg: Optional[float],
              mass: np.ndarray,
              sum_mass: float) -> float:
    if x_cg is None:
        x_cgi = (x * mass).sum() / sum_mass
    else:
        assert isinstance(x_cg, float), x_cg
        x_cgi = x_cg
    return x_cgi

def _get_xyz_cg(xyz: np.ndarray,
                xyz_cg: Optional[float],
                mass: np.ndarray,
                sum_mass: float) -> float:
    if xyz_cg is None:
        xyz_cgi = (xyz * mass[:, np.newaxis]).sum(axis=0) / sum_mass
    else:
        assert isinstance(xyz_cg, np.ndarray), xyz_cg
        xyz_cgi = xyz_cg
    assert len(xyz_cgi) == 3, xyz_cgi
    return xyz_cgi

def _get_inertia1(x: np.ndarray, xcg: float,
                  mass: np.ndarray,
                  inertia_cg: Optional[np.ndarray]) -> tuple[np.ndarray, np.ndarray]:
    dxcg = x - xcg
    if inertia_cg is not None:
        inertia_cgi = inertia_cg
    else:
        inertia_cgi = mass * dxcg ** 2
    return dxcg, inertia_cgi

def _get_inertia3(xyz: np.ndarray,
                  xyz_cg: np.ndarray,
                  mass: np.ndarray,
                  inertia_cg: Optional[np.ndarray]) -> tuple[np.ndarray, np.ndarray]:
    """
    Builds the per-station offset from the cg and the per-station inertia.

    When ``inertia_cg`` is None the inertia is the point-mass parallel-axis
    form ``m * [dy^2+dz^2, dx^2+dz^2, dx^2+dy^2, dx*dy, dx*dz, dy*dz]``, i.e.
    each station's own self-inertia is neglected.  Products of inertia are
    returned as the positive integrals.

    .. warning::
        **Passing ``inertia_cg`` raises ``UnboundLocalError``.**  ``dxyz_cg`` is
        only assigned in the ``else`` branch, but is returned unconditionally.
        Confirmed by execution.  The fix is to hoist the ``dxyz_cg`` assignment
        above the branch.

    Returns
    -------
    dxyz_cg : (n,3) float ndarray
        station offsets from the cg
    inertia_cgi : (n,6) float ndarray
        per-station inertia [Ixx, Iyy, Izz, Ixy, Ixz, Iyz]

    """
    if inertia_cg is not None:
        inertia_cgi = inertia_cg
    else:
        dxyz_cg = xyz - xyz_cg[np.newaxis, :]
        dx = dxyz_cg[:, 0]
        dy = dxyz_cg[:, 1]
        dz = dxyz_cg[:, 2]
        dx2 = dx * dx
        dy2 = dy * dy
        dz2 = dz * dz
        rxx = dy2 + dz2
        ryy = dx2 + dz2
        rzz = dx2 + dy2
        rxy = dx * dy
        ryz = dy * dz
        rxz = dx * dz
        radius_gyration2 = np.column_stack([rxx, ryy, rzz, rxy, rxz, ryz])
        inertia_cgi = mass[:, np.newaxis] * radius_gyration2
    return dxyz_cg, inertia_cgi

def _write_summary(data: dict[str, float], indent: str='  ') -> str:
    """
    Formats a name -> value dict into an aligned text block for logging.

    Floats print with 3 decimals; everything else (notably numpy arrays) falls
    back to ``str``.  Keys are expected to carry their own units, e.g.
    ``'mass (slinch)'``.

    """
    len_max_name = max((len(key) for key in data))
    msg = ''
    for key, value in data.items():
        if isinstance(value, float):
            fmt = indent + '%%-%ds = %%.3f\n' % len_max_name
        else:
            fmt = indent + '%%-%ds = %%s\n' % len_max_name
        msgi = fmt % (key, value)
        msg += msgi
    return msg

def fintegrate(y, x=None) -> np.ndarray:
    """
    Cumulative trapezoidal integration along a station axis, prepending 0 so
    the result is the same length as the input.

    Used to turn the differential (per-station) loads into the running
    shear/moment diagrams.

    Parameters
    ----------
    y : (n,) float ndarray
        the differential quantity
    x : (n,) float ndarray
        the station coordinate

    Returns
    -------
    xy : (n,) float ndarray
        the cumulative integral, starting at 0

    Notes
    -----
    - The descending-station case is handled by negating ``x`` (via a sign
      check on the endpoints only), which flips the sign of the integral.  A
      **non-monotonic** ``x`` passes this check and integrates incorrectly
      without any warning.
    - ``x`` is keyword-with-default ``None`` but is dereferenced immediately as
      ``x[0]``, so omitting it raises ``TypeError``; it is effectively required.

    """
    if x[0] > x[-1]:
        sign = -1
    else:
        sign = 1
    xy = np.hstack([0., cumulative_trapezoid(y, x=sign*x)])
    return xy

def plot_inertia(x,
                 exterior_force, inertial_linear_force, inertial_angular_force, total_force,
                 exterior_moment, inertial_linear_moment, inertial_angular_moment, total_moment,
                 locations: dict[str, float],
                 accel: float, alpha: float,
                 accel_units:str='in/s^2',
                 alpha_units:str='rad/s^2',
                 length_units:str='in', force_units: str='lbf',
                 integrate: bool=True,
                 case: str='') -> None:
    """
    Plots the load balance as a 2x2 grid of matplotlib axes.

    Layout::

        [0,0] applied force, linear inertial force, and their sum (F-ma)
        [0,1] F-ma, angular inertial force, and the total (should be ~0 net)
        [1,0] applied moment, linear inertial moment, and their sum
        [1,1] the above, the angular inertial moment, and the total

    Parameters
    ----------
    x : (n,) float ndarray
        the station coordinate
    exterior_*, inertial_linear_*, inertial_angular_*, total_* : (n,) float
        the eight load arrays returned by ``inertia_relief1``
    locations : dict[str, float]
        labelled vertical reference lines (nose, wing LE, tail, ...)
    accel : float
        rigid-body linear acceleration, shown in the force subplot title
    alpha : float
        rigid-body angular acceleration, shown in the moment subplot title
    integrate : bool; default=True
        if True, cumulatively integrate first, giving shear/moment diagrams;
        if False, plot the differential (per-station) loads
    case : str; default=''
        appended to the figure suptitle

    Notes
    -----
    - **1-D only.**  ``accel`` and ``alpha`` are formatted with ``:.6f``, which
      raises ``TypeError`` if handed the ``(3,)`` arrays that
      ``inertia_relief3`` works with.
    - ``length_units``/``force_units`` are local defaults here and are not
      passed down by ``inertia_relief1``, so the axis labels read ``in``/``lbf``
      regardless of the units the caller actually used.  The x-label is
      hardcoded to ``'Station (in)'``.
    - ``total_force`` is plotted as-is while the other three force curves are
      integrated when ``integrate=True``, so the ``[0,1]`` subplot mixes an
      integrated and a differential quantity.
    - Creates a figure but does not call ``show()``; the caller does.

    """
    fig = plt.figure()
    axes = fig.subplots(nrows=2, ncols=2)
    moment_units = f'{length_units}-{force_units}'

    word = ''
    if case:
        word = f': {case}'
    if integrate:
        #dword = ''
        exterior_force = fintegrate(exterior_force, x=x)
        inertial_linear_force = fintegrate(inertial_linear_force, x=x)
        inertial_angular_force = fintegrate(inertial_angular_force, x=x)
        exterior_moment = fintegrate(exterior_moment, x=x)
        inertial_linear_moment = fintegrate(inertial_linear_moment, x=x)
        inertial_angular_moment = fintegrate(inertial_angular_moment, x=x)
        total_moment = fintegrate(total_moment, x=x)
        fig.suptitle(f'Integrated Loads{word}')
    else:
        fig.suptitle(f'Differential Loads{word}')
        #dword = 'd'

    fma = exterior_force + inertial_linear_force
    mia = exterior_moment + inertial_linear_moment

    # force
    axes[0, 0].set_title(f'accel {accel:.6f} ({accel_units})')
    axes[0, 0].set_ylabel(f'Force ({force_units})')
    axes[0, 0].plot(x, exterior_force, label=f'F = exterior [{exterior_force.min():.0f}, {exterior_force.max():.0f}]', linestyle='--', )
    axes[0, 0].plot(x, inertial_linear_force, label=f'-ma = linear [{inertial_linear_force.min():.0f}, {inertial_linear_force.max():.0f}]', linestyle='--', )
    axes[0, 0].plot(x, fma, label=f'F-ma [{fma.min():.0f}, {fma.max():.0f}]')

    axes[0, 1].set_ylabel(f'Force ({force_units})')
    axes[0, 1].plot(x, fma, label=f'F-ma [{fma.min():.0f}, {fma.max():.0f}]', linestyle='--', )
    axes[0, 1].plot(x, inertial_angular_force, label=f'-I*alpha/dx = angular [{inertial_angular_force.min():.0f}, {inertial_angular_force.max():.0f}]', linestyle='--', )
    axes[0, 1].plot(x, total_force, label=f'total [{total_force.min():.0f}, {total_force.max():.0f}]')

    #----------------
    # moment
    axes[1, 0].set_title(f'alpha {alpha:.6f} ({alpha_units})')
    axes[1, 0].set_ylabel(f'Moment ({moment_units})')
    axes[1, 0].plot(x, exterior_moment, label=f'M = exterior [{exterior_moment.min():.3g}, {exterior_moment.max():.3g}]', linestyle='--', )
    axes[1, 0].plot(x, inertial_linear_moment, label=f'-I*alpha = linear [{inertial_linear_moment.min():.3g}, {inertial_linear_moment.max():.3g}]', linestyle='--', )
    axes[1, 0].plot(x, mia, label=f'M-I*alpha [{mia.min():.3g}, {mia.max():.3g}]')

    axes[1, 1].set_ylabel(f'Moment ({moment_units})')
    axes[1, 1].plot(x, mia, label=f'(F-ma)*dx [{mia.min():.3g}, {mia.max():.3g}]', linestyle='--', )
    axes[1, 1].plot(x, inertial_angular_moment, label=f'-I*alpha = angular [{inertial_angular_moment.min():.3g}, {inertial_angular_moment.max():.3g}]', linestyle='--', )
    axes[1, 1].plot(x, total_moment, label=f'total [{total_moment.min():.3g}, {total_moment.max():.3g}]')

    for ax in axes.ravel():
        for key, value in locations.items():
            ax.axvline(value, label=key, linestyle='--', color='k')
        ax.set_xlabel('Station (in)')
        ax.grid(True)
        ax.legend()

def test_inertia1():
    x = np.linspace(0., 10., num=51)
    mass = np.ones(x.shape, x.dtype)
    mass = np.ones(x.shape, x.dtype) + 0.1 * x
    #exterior_force = np.zeros(x.shape, x.dtype)
    #exterior_force[0] = 10.
    #exterior_force[1] = 5.
    #exterior_force = np.sin(x)
    exterior_force = mass * 3.
    #exterior_force = 2 * x

    inertia_relief1(
        x, exterior_force, mass, # inertia_cg=inertia_cg,
        xyz_cg=None,
        include_linear_inertia=True, include_angular_inertia=True,
        mass_units='slinch', length_units='in', force_units='lbf',
        plot_differential=True,
        plot_integrated=True,
        #integrate=False,
    )
    # inertia_relief(
    #     x, exterior_force, mass,
    #     xcg=None,
    #     include_linear_inertia=True, include_angular_inertia=True,
    #     mass_units='slinch', length_units='in', force_units='lbf',
    #     plot_differential=True,
    #     plot_integrated=True,
    #     #integrate=True,
    # )
    plt.show()

def test_inertia3():
    nx = 51
    xyz = np.zeros((nx, 3), dtype='float64')
    x = np.linspace(0., 10., num=51)
    xyz[:, 0] = x
    mass = np.ones(x.shape, dtype=x.dtype)
    mass = np.ones(x.shape, dtype=x.dtype) + 0.1 * x

    exterior_force = np.zeros((nx, 3), dtype='float64')
    exterior_force[:, 0] = mass * 3.
    inertia_relief3(xyz, exterior_force, mass)

if __name__ == '__main__':  # pragma: no cover
    #test_inertia1()
    test_inertia3()
