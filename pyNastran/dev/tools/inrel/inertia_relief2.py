"""
The goal of this file is to do inertia relief
on a nastran model.

Per Nastran QRG for the CONM2, the mass matrix about the
CG is:
    [M  0  0                ]
M = [0  M  0                |
    |0  0  M                |
    |0  0  0   I11 -I12 -I13|
    |0  0  0  -I12  I22 -I23|
    [0  0  0  -I13 -I23  I33]
So the inertia matrix has negative signs
on the off-diagonals.
"""
from __future__ import annotations
from typing import TYPE_CHECKING

import numpy as np
from pyNastran.bdf.cards.coordinate_systems import CORD2R
from pyNastran.bdf.mesh_utils.mass_properties import (
    mass_properties, transform_inertia)
if TYPE_CHECKING:
    from pyNastran.bdf.bdf import read_bdf, BDF


def get_mass_properties_array(model: BDF) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
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

def get_eigenvalues(imat: np.ndarray):  # pragma: no cover
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
    T2 = T**2
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

    Parameters
    ----------
    weight:  (nnode,3) float ndarray
        weight or mass oer node
        units of mass or weight (depending on wtmass)
    xyz : (nnode,3) float ndarray
        xyz locations of each mass in the global frame
    inertia : (nnode,6) float ndarray
        self-inertia (don't include mr^2)
        units of W*L^2 or M*L^2 depending on wtmass
    force : (nnode,3) float ndarray
        existing forces in the global frame
    moment : (nnode,3) float ndarray
        existing moments in the global frame
    wtmass : float; default=1.0
        weight to mass conversion
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

    # transform forces/moment to principal inertial frame
    fxyzt = np.einsum('nij,nj->ni',S[np.newaxis,:], force)
    mxyzt = np.einsum('nij,nj->ni',S[np.newaxis,:], moment)
    dxyzt = np.einsum('nij,nj->ni',S[np.newaxis,:], dxyz)

    shape = fxyzt.shape
    dtype = fxyzt.dtype

    # print(f'Mtt:\n{str(Mtt)}')
    i = S[:, 0]
    # j = S[:, 1]
    k = S[:, 2]
    # i = S[0, :]
    # k = S[2, :]
    inertia_ref1 = inertia_total

    xyz_ref1 = cg_total
    xyz_ref2 = xyz_ref1
    origin = xyz_ref1
    # i0 = np.array([1., 0., 0.])
    # k0 = np.array([0., 0., 1.])
    # coord1 = CORD2R(1, origin, origin+k0, origin+i0)
    coord2 = CORD2R(1, origin, origin+k, origin+i)
    assert cg_total.shape == (3,), cg_total.shape
    assert xyz_ref1.shape == (3,), xyz_ref1.shape
    assert xyz_ref2.shape == (3,), xyz_ref2.shape
    assert inertia_ref1.shape == (6,), inertia_ref1.shape
    # print(f'inertia_ref1 = {inertia_ref1}')
    inertia2 = transform_inertia(
        mass_total,
        # cg_total, xyz_ref1, xyz_ref2,
        np.zeros(3), np.zeros(3), np.zeros(3),
        inertia_ref1,
        coord2=coord2,
    )
    inertia_final = inertia2[:3]
    assert len(inertia_final) == 3, inertia_final

    # print(f'Mtt:\n{str(Mtt)}')
    # print(f'inertia2:\n{str(inertia2)}')
    # assert np.allclose(inertia2, Mtt)

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

    assert np.abs(fxyzt_out).max() < 1., fxyzt_out
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

    dxyzt_sq = dxyzt ** 2
    dxt_sq = dxyzt_sq[:, 0]
    dyt_sq = dxyzt_sq[:, 1]
    dzt_sq = dxyzt_sq[:, 2]

    # don't divide by 0
    ipos = np.where(inertia_final != 0)[0]

    # M = I * alpha
    alphat_moment = np.zeros(3, dtype='float64')
    alphat_moment[ipos] = mxyzt2_total[ipos] / inertia_final[ipos]
    mxt_dm = -alphat_moment[0] * mass * (dyt_sq + dzt_sq)
    myt_dm = -alphat_moment[1] * mass * (dxt_sq + dzt_sq)
    mzt_dm = -alphat_moment[2] * mass * (dxt_sq + dyt_sq)
    mxyzt_dm = np.column_stack([mxt_dm, myt_dm, mzt_dm])
    # print(f'mxyzt_dm_total = {mxyzt_dm.sum(axis=0)}')

    mxyzt_delta = mxyzt_df + mxyzt_delta_self + mxyzt_dm
    # mxyzt_out = mxyzt2 + mxyzt_delta

    # wrong
    St = S.T  # also tried S instead of S.T
    fxyz_delta_out = np.einsum('nij,nj->ni',St[np.newaxis,:], fxyzt_delta)
    mxyz_delta_out = np.einsum('nij,nj->ni',St[np.newaxis,:], mxyzt_delta)

    # St = S.T  # also tried S instead of S.T
    # fxyz_delta_out = np.einsum('nij,ni->nj',St[np.newaxis,:], fxyzt_delta)
    # mxyz_delta_out = np.einsum('nij,ni->nj',St[np.newaxis,:], mxyzt_delta)

    return fxyz_delta_out, mxyz_delta_out


# if __name__ == '__main__':
#     inertia_relief_from_model(None)
