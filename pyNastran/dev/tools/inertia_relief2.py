"""
The goal of this file is to do inertia relief
on a nastran model.
"""
from __future__ import annotations
from typing import TYPE_CHECKING

import numpy as np
from pyNastran.bdf.cards.coordinate_systems import CORD2R
from pyNastran.bdf.mesh_utils.mass_properties import transform_inertia
if TYPE_CHECKING:
    from pyNastran.bdf.bdf import read_bdf, BDF


def get_mass_properties_array(model: BDF) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    # from pyNastran.bdf.bdf import read_bdf
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

def _mass_properties_total(mass: float,
                           xyz: np.ndarray,
                           inertia: np.ndarray,
                           ) -> tuple[float, np.ndarray, np.ndarray]:
    mass_total = mass.sum()
    # assert np.allclose(mass_total, 5.), mass_total

    cg_total = (xyz * mass[:, np.newaxis]).sum(axis=0) / mass_total
    assert len(cg_total) == 3, cg_total
    # print(f'mass = {mass}')
    # print(f'cg = {cg_total}')
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
        [ixx, ixy, ixz],
        [ixy, iyy, iyz],
        [ixz, iyz, izz],
    ])
    print(f'inertia_mat:\n{str(inertia_mat)}')

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
    else:
        S = np.eye(3, dtype=Mtt_.dtype)
    Mtt = S.T @ Mtt_ @ S # Mt


    # inertia_mat = np.array([
    #     [ixx, ixy, ixz],
    #     [ixy, iyy, iyz],
    #     [ixz, iyz, izz],
    # ])
    inertia_total = np.array([
        Mtt[0, 0], Mtt[1, 1], Mtt[2, 2],
        Mtt[0, 1], Mtt[0, 2], Mtt[1, 2],
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


def inertia_relief(mass: np.ndarray,
                   xyz: np.ndarray,
                   inertia: np.ndarray,
                   force: np.ndarray,
                   moment: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    mass_total, cg_total, inertia_total, Mtt, S = _mass_properties_total(
        mass, xyz, inertia)
    assert inertia_total.shape == (6,), inertia_total.shape

    # TODO: Should this be S.T?
    # Ft = S @ F
    print(f'xyz:\n{xyz}')
    print(f'cg_total = {cg_total}')
    dxyz = xyz - cg_total[np.newaxis, :]

    # transform forces/moment to principal inertial frame
    fxyzt = np.einsum('nij,nj->ni',S[np.newaxis,:], force)
    mxyzt = np.einsum('nij,nj->ni',S[np.newaxis,:], moment)
    dxyzt = np.einsum('nij,nj->ni',S[np.newaxis,:], dxyz)

    print(f'Mtt:\n{str(Mtt)}')
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
    print(f'inertia_ref1 = {inertia_ref1}')
    inertia2 = transform_inertia(
        mass_total,
        # cg_total, xyz_ref1, xyz_ref2,
        np.zeros(3), np.zeros(3), np.zeros(3),
        inertia_ref1,
        coord2=coord2,
    )
    inertia_final = inertia2[:3]
    assert len(inertia_final) == 3, inertia_final

    print(f'Mtt:\n{str(Mtt)}')
    print(f'inertia2:\n{str(inertia2)}')
    # assert np.allclose(inertia2, Mtt)

    fxyzt_total = fxyzt.sum(axis=0)
    print(f'fxyzt_total = {fxyzt_total}')
    assert len(fxyzt_total) == 3, fxyzt_total
    # F = m*a
    accelt = fxyzt_total / mass_total

    fxyzt_delta = -mass[:, np.newaxis] * accelt[np.newaxis, :]
    print(f'dxyzt =\n{dxyzt}')
    mxyzt_delta0 = np.cross(dxyzt, fxyzt_delta, axis=1)
    # print(f'mxyzt_delta0 = {mxyzt_delta0}')
    print(f'mxyzt_delta0_total = {mxyzt_delta0.sum(axis=0)}')
    assert fxyzt.shape == mxyzt.shape, (fxyzt.shape, mxyzt.shape)

    fxyzt_out = fxyzt + fxyzt_delta
    assert np.abs(fxyzt_out).max() < 1., fxyzt_out
    mxyzt1 = mxyzt + mxyzt_delta0
    mxyzt1_total = mxyzt1.sum(axis=0)
    print(f' mxyzt_total        = {mxyzt.sum(axis=0)}')
    print(f'+mxyzt_delta0_total = {mxyzt_delta0.sum(axis=0)}')
    print(f'=mxyzt1_total       = {mxyzt1_total}')

    # don't divide by 0
    ipos = np.where(inertia_final != 0)[0]

    # M = I*alpha = I0*alpha + m*r^2*alpha
    #
    # handle self-inertia (probably very small)
    alphat1 = np.zeros(3, dtype='float64')
    alphat1[ipos] = mxyzt1_total[ipos] / inertia_final[ipos]

    mxyzt_delta1 = -(alphat1 * inertia_final)[np.newaxis, :]
    mxyzt2 = mxyzt + mxyzt_delta0 + mxyzt_delta1
    mxyzt2_total = mxyzt2.sum(axis=0)

    dxyzt_sq = dxyzt ** 2
    dxt_sq = dxyzt_sq[:, 0]
    dyt_sq = dxyzt_sq[:, 1]
    dzt_sq = dxyzt_sq[:, 2]

    # M = I * alpha
    alphat2 = np.zeros(3, dtype='float64')
    alphat2[ipos] = mxyzt2_total[ipos] / inertia_final[ipos]
    mxt_delta2 = -alphat2[0] * mass * (dyt_sq + dzt_sq)
    myt_delta2 = -alphat2[1] * mass * (dxt_sq + dzt_sq)
    mzt_delta2 = -alphat2[2] * mass * (dxt_sq + dyt_sq)
    mxyzt_delta2 = np.column_stack([mxt_delta2, myt_delta2, mzt_delta2])

    mxyzt_delta = mxyzt_delta1 + mxyzt_delta2
    mxyzt_out = mxyzt2 + mxyzt_delta

    St = S.T
    fxyz_delta_out = np.einsum('nij,nj->ni',St[np.newaxis,:], fxyzt_delta)
    mxyz_delta_out = np.einsum('nij,nj->ni',St[np.newaxis,:], mxyzt_delta)
    return fxyz_delta_out, mxyz_delta_out


def test_bar():
    nnode = 11
    x = np.linspace(0., 100., num=nnode)
    xyz = np.zeros((nnode, 3), dtype='float64')
    xyz[:, 0] = x

    mass_total = 11.
    mass = np.ones(nnode, dtype='float64') * mass_total / nnode
    mass[0] = 5.
    fx = 0.1 * mass
    fy = 0.0 * mass
    fz = 0.5 * mass
    force = np.column_stack([fx, fy, fz])
    moment = force * 0.0
    inertia = np.zeros((nnode, 6), dtype='float64')
    dforce, dmoment = inertia_relief(
        mass, xyz, inertia,
        force, moment)
    force_out = force + dforce
    moment_out = moment + dmoment
    print('df:\n', dforce)
    print('dm:\n', dmoment)
    # print('force_out:\n', force_out)


if __name__ == '__main__':
    test_bar()
    # inertia_relief_from_model(None)
