from __future__ import annotations
from typing import Union, TYPE_CHECKING
import numpy as np

from .solid_volume import volume_chexa, volume_cpenta
from pyNastran.dev.bdf_vectorized3.cards.base_card import (
    searchsorted_filter)

if TYPE_CHECKING:
    from pyNastran.dev.bdf_vectorized3.bdf import BDF
    #from pyNastran.bdf.cards.materials import MAT1, MAT8
    from pyNastran.dev.bdf_vectorized3.cards.grid import GRID
    from .shell_properties import PCOMP, PSHELL, PLPLANE


def tri_area(grid: GRID, nodes: np.ndarray) -> np.ndarray:
    xyz = grid.xyz_cid0()
    nid = grid.node_id
    nelements, nnodes = nodes.shape
    assert nnodes == 3, nodes.shape
    inode = np.searchsorted(nid, nodes)
    actual_nodes = nid[inode]
    assert np.array_equal(actual_nodes, nodes)
    in1 = inode[:, 0]
    in2 = inode[:, 1]
    in3 = inode[:, 2]
    xyz1 = xyz[in1, :]
    xyz2 = xyz[in2, :]
    xyz3 = xyz[in3, :]
    a = xyz1 - xyz2
    b = xyz1 - xyz3
    area = 0.5 * np.linalg.norm(np.cross(a, b), axis=1)
    assert len(area) == nelements
    return area

def tri_centroid(grid: GRID, nodes: np.ndarray) -> np.ndarray:
    xyz = grid.xyz_cid0()
    nid = grid.node_id
    inode = np.searchsorted(nid, nodes)
    assert np.array_equal(nid[inode], nodes)
    in1 = inode[:, 0]
    in2 = inode[:, 1]
    in3 = inode[:, 2]
    xyz1 = xyz[in1, :]
    xyz2 = xyz[in2, :]
    xyz3 = xyz[in3, :]
    centroid = (xyz1 + xyz2 + xyz3) / 3.
    return centroid

def tri_area_centroid_normal(grid: GRID, nodes: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    xyz = grid.xyz_cid0()
    nid = grid.node_id
    nelements, nnodes = nodes.shape
    assert nnodes == 3, nodes.shape
    inode = np.searchsorted(nid, nodes)
    actual_nodes = nid[inode]
    assert np.array_equal(actual_nodes, nodes)
    in1 = inode[:, 0]
    in2 = inode[:, 1]
    in3 = inode[:, 2]
    xyz1 = xyz[in1, :]
    xyz2 = xyz[in2, :]
    xyz3 = xyz[in3, :]
    a = xyz1 - xyz2
    b = xyz1 - xyz3

    normal = np.cross(a, b)
    assert normal.shape[0] == nelements

    norm = np.linalg.norm(normal, axis=1)
    area = 0.5 * norm
    assert len(area) == nelements

    centroid = (xyz1 + xyz2 + xyz3) / 3.
    assert centroid.shape[0] == nelements

    assert normal.shape == (nelements, 3)
    ipos = norm > 0.
    ibad = ~ipos
    if ibad.sum() == 0:
        unit_normal = normal / norm[:, np.newaxis]
    else:
        unit_normal = np.full(normal.shape, np.nan, normal.dtype)
        unit_normal[ipos, :] = normal[ipos, :] / norm[ipos, np.newaxis]
        #print('norm =', norm)
        #print('unit_normal =', unit_normal)
        print('a[bad,:] =', a[ibad, :])
        print('  n1 =', nodes[ibad, 0])
        print('  n2 =', nodes[ibad, 1])
        print('  n3 =', nodes[ibad, 2])
        print('  xyz1 =', xyz1[ibad, :])
        print('  xyz2 =', xyz2[ibad, :])
        print('  xyz3 =', xyz3[ibad, :])
        print('b[bad,:] =', b[ibad, :])

    assert unit_normal.shape == (nelements, 3)
    return area, centroid, unit_normal

def quad_area(grid: GRID, nodes: np.ndarray) -> np.ndarray:
    nelements, nnodes = nodes.shape
    assert nnodes == 4, nnodes
    nid = grid.node_id
    xyz = grid.xyz_cid0()
    inode = np.searchsorted(nid, nodes)
    assert np.array_equal(nid[inode], nodes)
    in1 = inode[:, 0]
    in2 = inode[:, 1]
    in3 = inode[:, 2]
    in4 = inode[:, 3]
    xyz1 = xyz[in1, :]
    xyz2 = xyz[in2, :]
    xyz3 = xyz[in3, :]
    xyz4 = xyz[in4, :]
    #area = 0.5 * norm(cross(n3-n1, n4-n2))
    a = xyz3 - xyz1
    b = xyz4 - xyz2
    area = 0.5 * np.linalg.norm(np.cross(a, b), axis=1)
    assert len(area) == nelements
    return area

def quad_centroid(grid: GRID, nodes: np.ndarray) -> np.ndarray:
    nelements, nnodes = nodes.shape
    assert nnodes == 4, nnodes
    xyz = grid.xyz_cid0()
    nid = grid.node_id
    inode = np.searchsorted(nid, nodes)
    assert np.array_equal(nid[inode], nodes)
    in1 = inode[:, 0]
    in2 = inode[:, 1]
    in3 = inode[:, 2]
    in4 = inode[:, 3]
    xyz1 = xyz[in1, :]
    xyz2 = xyz[in2, :]
    xyz3 = xyz[in3, :]
    xyz4 = xyz[in4, :]
    centroid = (xyz1 + xyz2 + xyz3 + xyz4) / 4.
    assert centroid.shape[0] == nelements
    return centroid

def quad_area_centroid_normal(grid: GRID, nodes: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    nelements, nnodes = nodes.shape
    assert nnodes == 4, nnodes
    nid = grid.node_id
    xyz = grid.xyz_cid0()
    inode = np.searchsorted(nid, nodes)
    assert np.array_equal(nid[inode], nodes)
    in1 = inode[:, 0]
    in2 = inode[:, 1]
    in3 = inode[:, 2]
    in4 = inode[:, 3]
    xyz1 = xyz[in1, :]
    xyz2 = xyz[in2, :]
    xyz3 = xyz[in3, :]
    xyz4 = xyz[in4, :]
    #area = 0.5 * norm(cross(n3-n1, n4-n2))
    a = xyz3 - xyz1
    b = xyz4 - xyz2
    normal = np.cross(a, b)
    assert normal.shape[0] == nelements

    norm = np.linalg.norm(normal, axis=1)
    area = 0.5 * norm
    assert len(area) == nelements

    centroid = (xyz1 + xyz2 + xyz3 + xyz4) / 4.
    assert centroid.shape[0] == nelements

    assert normal.shape == (nelements, 3)
    unit_normal = normal / norm[:, np.newaxis]

    assert unit_normal.shape == (nelements, 3)
    return area, centroid, unit_normal


def tri_volume(grid: GRID,
               nodes: np.ndarray,
               dthickness: np.ndarray) -> np.ndarray:
    """
    Exact quad volume considering differential thicknesses

    Parameters
    ----------
    dthickness : (nelement, 3) float ndarray
        differential thickness on the quad
    """

    nelements, nnodes = nodes.shape
    assert nnodes == 3, nnodes
    nid = grid.node_id
    xyz = grid.xyz_cid0()
    inode = np.searchsorted(nid, nodes)
    assert np.array_equal(nid[inode], nodes)
    in1 = inode[:, 0]
    in2 = inode[:, 1]
    in3 = inode[:, 2]
    xyz1 = xyz[in1, :]
    xyz2 = xyz[in2, :]
    xyz3 = xyz[in3, :]

    a = xyz2 - xyz1
    b = xyz3 - xyz1
    normal = np.cross(a, b)
    assert normal.shape[0] == nelements

    norm = np.linalg.norm(normal, axis=1)
    assert normal.shape == (nelements, 3)
    unit_normal = normal / norm[:, np.newaxis]

    iequal = (dthickness.max(axis=0) == dthickness.min(axis=0))
    inot_equal = ~iequal
    nequal = iequal.sum()
    nnot_equal = inot_equal.sum()
    volume = np.full(nelements, np.nan, dtype='float64')
    if nequal:
        thicknessi = dthickness[iequal, 0]
        areai = 0.5 * norm[iequal]
        volumei = areai * thicknessi
        volume[iequal] = volumei

    if nnot_equal:
        unit_normali = unit_normal[inot_equal, :]
        xyz4 = xyz1[inot_equal, :] + unit_normali * dthickness[inot_equal, 0]
        xyz5 = xyz2[inot_equal, :] + unit_normali * dthickness[inot_equal, 1]
        xyz6 = xyz3[inot_equal, :] + unit_normali * dthickness[inot_equal, 2]
        volumei = volume_cpenta(
            xyz1[inot_equal, :], xyz2[inot_equal, :], xyz3[inot_equal, :],
            xyz4, xyz5, xyz6)
        volume[iequal] = volumei
    return volume

def quad_volume(grid: GRID,
                nodes: np.ndarray,
                dthickness: np.ndarray) -> np.ndarray:
    """
    Exact quad volume considering differential thicknesses

    Parameters
    ----------
    dthickness : (nelement, 4) float ndarray
        differential thickness on the quad
    """

    nelements, nnodes = nodes.shape
    assert nnodes == 4, nnodes
    nid = grid.node_id
    xyz = grid.xyz_cid0()
    inode = np.searchsorted(nid, nodes)
    assert np.array_equal(nid[inode], nodes)
    in1 = inode[:, 0]
    in2 = inode[:, 1]
    in3 = inode[:, 2]
    in4 = inode[:, 3]
    xyz1 = xyz[in1, :]
    xyz2 = xyz[in2, :]
    xyz3 = xyz[in3, :]
    xyz4 = xyz[in4, :]

    a = xyz3 - xyz1
    b = xyz4 - xyz2
    normal = np.cross(a, b)
    assert normal.shape[0] == nelements

    norm = np.linalg.norm(normal, axis=1)
    assert normal.shape == (nelements, 3)
    unit_normal = normal / norm[:, np.newaxis]

    iequal = (dthickness.max(axis=0) == dthickness.min(axis=0))
    inot_equal = ~iequal
    nequal = iequal.sum()
    nnot_equal = inot_equal.sum()
    volume = np.full(nelements, np.nan, dtype='float64')
    if nequal:
        thicknessi = dthickness[iequal, 0]
        areai = 0.5 * norm[iequal]
        volumei = areai * thicknessi
        volume[iequal] = volumei

    if nnot_equal:
        unit_normali = unit_normal[inot_equal, :]
        xyz5 = xyz1[inot_equal, :] + unit_normali * dthickness[inot_equal, 0]
        xyz6 = xyz2[inot_equal, :] + unit_normali * dthickness[inot_equal, 1]
        xyz7 = xyz3[inot_equal, :] + unit_normali * dthickness[inot_equal, 2]
        xyz8 = xyz4[inot_equal, :] + unit_normali * dthickness[inot_equal, 3]
        volumei = volume_chexa(xyz1[inot_equal, :], xyz2[inot_equal, :],
                               xyz3[inot_equal, :], xyz4[inot_equal, :],
                               xyz5, xyz6, xyz7, xyz8)
        volume[iequal] = volumei
    return volume




def shell_thickness(model: BDF,
                    tflag: np.ndarray,
                    T: np.ndarray,
                    property_id: np.ndarray,
                    allowed_properties: list[Union[PCOMP, PSHELL, PLPLANE]]) -> np.ndarray:
    log = model.log
    thickness = np.full(len(property_id), np.nan, dtype='float64')
    assert len(allowed_properties) > 0, allowed_properties
    for prop in allowed_properties:
        ilookup, iall = searchsorted_filter(prop.property_id, property_id)
        if len(iall) == 0:
            continue
        ti = prop.total_thickness()
        ti_all = ti[iall]

        # set the thickness, even if it's nan
        thickness[ilookup] = ti_all
        if prop.type != 'PSHELL':
            continue

        inan = np.isnan(ti_all)
        if not np.any(inan):
            continue

        #print('inan', inan)
        tflag_nan = tflag[ilookup]
        t_nan = T[ilookup, :]
        #print('tflag_nan', tflag_nan)
        #print('t_nan', t_nan)

        i0 = (tflag_nan == 0)
        i1 = ~i0
        if i0.sum():
            mean_thickness = t_nan[i0, :].mean(axis=1)
            ti_all[i0] = mean_thickness
        if i1.sum():
            scale = t_nan[i1, :].mean(axis=1)
            ti_all[i1] *= scale
            raise RuntimeError('tflag=1')
        #log.error(ti_all)
        # tflag=0: Thickness of element at grid points G1 through G4
        # TFLAG=1:Tthickness becomes a product of Ti and the thickness
        # on the PSHELL card. Ti is ignored for hyperelastic elements.
        # See Remark 6. (Real > 0.0 or blank. See Remark 4 for the default.)
        inan = np.isnan(ti_all)
        if np.any(inan):
            msg = (
                f'tflag={tflag_nan[inan]}\n'
                f'T={ti_all[inan, :]}')
            log.error(msg)
            raise RuntimeError(msg)
        thickness[ilookup] = ti_all

    inan = np.isnan(thickness)
    if np.any(inan):
        msg = f'Thickness has nan\nt={thickness}'
        log.error(msg)
        raise RuntimeError(msg)
    return thickness

def shell_total_thickness(property_id: np.ndarray,
                          allowed_properties: list[Union[PCOMP, PSHELL, PLPLANE]],
                          ) -> np.ndarray:
    thickness = np.full(len(property_id), np.nan, dtype='float64')
    assert len(allowed_properties) > 0, allowed_properties
    for prop in allowed_properties:
        ilookup, iall = searchsorted_filter(prop.property_id, property_id)
        if len(iall) == 0:
            continue
        thicknessi = prop.total_thickness()
        thickness[ilookup] = thicknessi[iall]
    return thickness

def shell_nonstructural_mass(property_id: np.ndarray,
                             allowed_properties: list[Union[PCOMP, PSHELL, PLPLANE]],
                             ) -> np.ndarray:
    nsm = np.full(len(property_id), np.nan, dtype='float64')
    assert len(allowed_properties) > 0, allowed_properties
    for prop in allowed_properties:
        ilookup, iall = searchsorted_filter(prop.property_id, property_id)
        if len(iall) == 0:
            continue
        nsmi = prop.nsm
        nsm[ilookup] = nsmi[iall]
    return nsm

def shell_mass_per_area(model: BDF,
                        tflag: np.ndarray,
                        T: np.ndarray,
                        property_id: np.ndarray,
                        allowed_properties: list[Union[PCOMP, PSHELL, PLPLANE]],
                        ) -> np.ndarray:
    nelement = len(property_id)
    assert nelement > 0, property_id
    mass_per_area = np.full(nelement, np.nan, dtype='float64')
    assert len(allowed_properties) > 0, allowed_properties

    for prop in allowed_properties:
        ilookup, iall = searchsorted_filter(prop.property_id, property_id)
        if len(iall) == 0:
            continue

        if prop.type in {'PCOMP', 'PLPLANE'}:
            mass_per_areai = prop.mass_per_area()
            mass_per_areai_all = mass_per_areai[iall]
            mass_per_area[ilookup] = mass_per_areai_all
            continue

        assert prop.type == 'PSHELL', prop.type
        nsm, rho, ti = prop.nsm_rho_thickness()
        nsm_all = nsm[iall]
        rho_all = rho[iall]
        ti_all = ti[iall]
        mass_per_areai_all = nsm_all + rho_all * ti_all
        mass_per_area[ilookup] = mass_per_areai_all

        inan = np.isnan(ti_all)
        if not np.any(inan):
            continue

        #print('inan', inan)
        tflag_nan = tflag[ilookup]
        t_nan = T[ilookup, :]
        #print('tflag_nan', tflag_nan)
        #print('t_nan', t_nan)

        i0 = (tflag_nan == 0)
        i1 = ~i0
        if i0.sum():
            mean_thickness = t_nan[i0, :].mean(axis=1)
            ti_all[i0] = mean_thickness
        if i1.sum():
            scale = t_nan[i1, :].mean(axis=1)
            ti_all[i1] *= scale
            raise RuntimeError('tflag=1')
        #log.error(ti_all)
        # tflag=0: Thickness of element at grid points G1 through G4
        # TFLAG=1:Tthickness becomes a product of Ti and the thickness
        # on the PSHELL card. Ti is ignored for hyperelastic elements.
        # See Remark 6. (Real > 0.0 or blank. See Remark 4 for the default.)
        inan = np.isnan(ti_all)
        if np.any(inan):
            msg = (
                f'tflag={tflag_nan[inan]}\n'
                f'T={ti_all[inan, :]}')
            model.log.error(msg)
            raise RuntimeError(msg)

        mass_per_area_all = nsm_all + rho_all * ti_all
        #nsm_all = nsm[iall]
        #rho_all = rho[iall]
        #ti_all = ti[iall]
        mass_per_area[ilookup] = mass_per_area_all

        inan = np.isnan(ti_all)
        if np.any(inan):
            msg = f'Thickness has nan\nt={ti_all}'
            model.log.error(msg)
            raise RuntimeError(msg)
    assert nelement > 0, nelement
    assert len(mass_per_area) == nelement, mass_per_area
    return mass_per_area


def shell_mass_per_area_breakdown(model: BDF,
                                  tflag: np.ndarray,
                                  T: np.ndarray,
                                  property_id: np.ndarray,
                                  allowed_properties: list[Union[PCOMP, PSHELL, PLPLANE]],
                                  ) -> np.ndarray:
    """
    PCOMP:    [nsm, nan, nan, mass_per_area]
    PLPLANE:  [nan, nan, nan, mass_per_area]
    PSHELL:   [nsm, rho, t,   mass_per_area]
    """
    nelement = len(property_id)
    assert nelement > 0, property_id
    mass_per_area_breakdown = np.full((nelement, 4), np.nan, dtype='float64')
    assert len(allowed_properties) > 0, allowed_properties

    for prop in allowed_properties:
        ilookup, iall = searchsorted_filter(prop.property_id, property_id)
        if len(iall) == 0:
            continue

        if prop.type in {'PCOMP', 'PLPLANE'}:
            if prop.type == 'PCOMP':
                breakdowni = prop.mass_per_area_breakdown() # nsm, mass_per_area
                assert breakdowni.shape[1] == 2, breakdowni.shape
                mass_per_area_breakdown[ilookup, 0] = breakdowni[iall, 0]
                mass_per_area_breakdown[ilookup, 3] = breakdowni[iall, 1]
            elif prop.type == 'PLPLANE':
                mass_per_areai = prop.mass_per_area()
                mass_per_areai_all = mass_per_areai[iall]
                mass_per_area_breakdown[ilookup, 3] = mass_per_areai_all
            else:  ## pragma: no cover
                raise RuntimeError(prop.type)
            #mass_per_area_breakdown[ilookup, 0] = 0.
            #mass_per_area_breakdown[ilookup, 1] = 0.
            #mass_per_area_breakdown[ilookup, 2] = 0.
            continue

        assert prop.type == 'PSHELL', prop.type
        nsm, rho, ti = prop.nsm_rho_thickness()
        nsm_all = nsm[iall]
        rho_all = rho[iall]
        ti_all = ti[iall]
        mass_per_areai_all = nsm_all + rho_all * ti_all
        mass_per_area_breakdown[ilookup, 0] = mass_per_areai_all
        mass_per_area_breakdown[ilookup, 1] = rho_all
        mass_per_area_breakdown[ilookup, 2] = ti_all
        mass_per_area_breakdown[ilookup, 3] = mass_per_areai_all

        inan = np.isnan(ti_all)
        if not np.any(inan):
            continue

        #print('inan', inan)
        tflag_nan = tflag[ilookup]
        t_nan = T[ilookup, :]
        #print('tflag_nan', tflag_nan)
        #print('t_nan', t_nan)

        i0 = (tflag_nan == 0)
        i1 = ~i0
        if i0.sum():
            mean_thickness = t_nan[i0, :].mean(axis=1)
            ti_all[i0] = mean_thickness
        if i1.sum():
            scale = t_nan[i1, :].mean(axis=1)
            ti_all[i1] *= scale
            raise RuntimeError('tflag=1')
        #log.error(ti_all)
        # tflag=0: Thickness of element at grid points G1 through G4
        # TFLAG=1:Tthickness becomes a product of Ti and the thickness
        # on the PSHELL card. Ti is ignored for hyperelastic elements.
        # See Remark 6. (Real > 0.0 or blank. See Remark 4 for the default.)
        inan = np.isnan(ti_all)
        if np.any(inan):
            msg = (
                f'tflag={tflag_nan[inan]}\n'
                f'T={ti_all[inan, :]}')
            model.log.error(msg)
            raise RuntimeError(msg)

        mass_per_area_all = nsm_all + rho_all * ti_all
        #nsm_all = nsm[iall]
        #rho_all = rho[iall]
        #ti_all = ti[iall]
        mass_per_area_breakdown[ilookup, 0] = nsm_all
        mass_per_area_breakdown[ilookup, 1] = rho_all
        mass_per_area_breakdown[ilookup, 2] = ti_all
        mass_per_area_breakdown[ilookup, 3] = mass_per_area_all

        inan = np.isnan(ti_all)
        if np.any(inan):
            msg = f'Thickness has nan\nt={ti_all}'
            model.log.error(msg)
            raise RuntimeError(msg)
    assert nelement > 0, nelement
    assert len(mass_per_area_breakdown) == nelement, mass_per_area_breakdown
    return mass_per_area_breakdown
