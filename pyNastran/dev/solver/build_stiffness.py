from __future__ import annotations
from itertools import count
from typing import Tuple, Any, TYPE_CHECKING
#from typing import List, Dict, Tuple, Union, Any

import numpy as np
import scipy.sparse as sci_sparse

from pyNastran.bdf.cards.elements.shell import transform_shell_material_coordinate_system
from .utils import lambda1d
#from pyNastran.bdf.cards.elements.bars import get_bar_vector, get_bar_yz_transform
if TYPE_CHECKING:  # pragma: no cover
    from pyNastran.bdf.bdf import BDF, CQUAD4, CBAR # , CBEAM
    #from pyNastran.bdf.cards.elements.shell import CQUAD4



def build_Kbb(model: BDF, dof_map, ndof, idtype='int32', fdtype='float32') -> Tuple[np.array, Any]:
    """[K] = d{P}/dx"""
    Kbb = np.zeros((ndof, ndof), dtype=fdtype)
    Kbbs = sci_sparse.dok_matrix((ndof, ndof), dtype=fdtype)
    #print(dof_map)

    #_get_loadid_ndof(model, subcase_id)
    out = model.get_xyz_in_coord_array(cid=0, fdtype=fdtype, idtype=idtype)
    nid_cp_cd, xyz_cid0, xyz_cp, unused_icd_transform, unused_icp_transform = out
    all_nids = nid_cp_cd[:, 0]
    del xyz_cp, nid_cp_cd

    nelements = 0
    nelements += _build_kbb_celas1(model, Kbb, Kbbs, dof_map)
    nelements += _build_kbb_celas2(model, Kbb, Kbbs, dof_map)
    nelements += _build_kbb_conrod(model, Kbb, Kbbs, dof_map)
    nelements += _build_kbb_crod(model, Kbb, Kbbs, dof_map)
    nelements += _build_kbb_ctube(model, Kbb, Kbbs, dof_map)
    nelements += _build_kbb_cbar(model, Kbb, Kbbs, dof_map)
    nelements += _build_kbb_cbeam(model, Kbb, Kbbs, dof_map,
                                  all_nids, xyz_cid0, idtype='int32', fdtype='float64')
    nelements += _build_kbb_cquad4(model, Kbb, Kbbs, dof_map,
                                   all_nids, xyz_cid0, idtype='int32', fdtype='float64')
    assert nelements > 0, nelements
    Kbbs2 = Kbbs.tocsc()
    Kbb2 = Kbbs2.toarray()
    error = np.linalg.norm(Kbb - Kbb2)
    if error > 1e-12:
        model.log.warning(f'error = {error}')
    return Kbb, Kbbs2


def _build_kbb_celas1(model: BDF, Kbb, Kbbs, dof_map):
    """fill the CELAS1 Kbb matrix"""
    eids = model._type_to_id_map['CELAS1']
    for eid in eids:
        elem = model.elements[eid]
        ki = elem.K()
        #print(elem, ki)
        #print(elem.get_stats())
        _build_kbbi_celas12(Kbb, Kbbs, dof_map, elem, ki)
    return len(eids)

def _build_kbb_celas2(model: BDF, Kbb, Kbbs, dof_map):
    """fill the CELAS2 Kbb matrix"""
    eids = model._type_to_id_map['CELAS2']
    #celas3s = model._type_to_id_map['CELAS3']
    #celas4s = model._type_to_id_map['CELAS4']
    for eid in eids:
        elem = model.elements[eid]
        ki = elem.K()
        #print(elem, ki)
        #print(elem.get_stats())
        _build_kbbi_celas12(Kbb, Kbbs, dof_map, elem, ki)
    return len(eids)

def _build_kbbi_celas12(Kbb, Kbbs, dof_map, elem, ki):
    """fill the CELASx Kbb matrix"""
    nid1, nid2 = elem.nodes
    c1, c2 = elem.c1, elem.c2
    i = dof_map[(nid1, c1)]
    j = dof_map[(nid2, c2)]
    k = ki * np.array([[1, -1,],
                       [-1, 1]])
    ibe = [
        (i, 0),
        (j, 1),
    ]
    for ib1, ie1 in ibe:
        for ib2, ie2 in ibe:
            Kbb[ib1, ib2] += k[ie1, ie2]
            Kbbs[ib1, ib2] += k[ie1, ie2]
    #Kbb[j, i] += ki
    #Kbb[i, j] += ki
    #del i, j, ki, nid1, nid2, c1, c2

def _build_kbb_cbar(model, Kbb, Kbbs, dof_map, fdtype='float64'):
    """fill the CBAR Kbb matrix using an Euler-Bernoulli beam"""
    eids = model._type_to_id_map['CBAR']
    nelements = len(eids)
    if nelements == 0:
        return nelements

    for eid in eids:
        elem = model.elements[eid]  # type: CBAR
        nid1, nid2 = elem.nodes
        pid_ref = elem.pid_ref
        mat = pid_ref.mid_ref

        #is_passed, (wa, wb, ihat, jhat, khat) = elem.get_axes(model)
        #T = np.vstack([ihat, jhat, khat])
        #z = np.zeros((3, 3), dtype='float64')
        prop = elem.pid_ref
        mat = prop.mid_ref
        I1 = prop.I11()
        I2 = prop.I22()
        unused_I12 = prop.I12()
        #J = prop.J()
        #E = mat.E()
        #G = mat.G()
        z = np.zeros((3, 3), dtype='float64')
        T = z
        unused_Teb = np.block([
            [T, z],
            [z, T],
        ])
        is_passed, (wa, wb, ihat, jhat, khat) = elem.get_axes(model)
        print(wa, wb)
        xyz1 = elem.nodes_ref[0].get_position() + wa
        xyz2 = elem.nodes_ref[1].get_position() + wb
        dxyz = xyz2 - xyz1
        L = np.linalg.norm(dxyz)
        pid_ref = elem.pid_ref
        mat = pid_ref.mid_ref
        T = np.vstack([ihat, jhat, khat])
        z = np.zeros((3, 3), dtype=fdtype)
        Teb = np.block([
            [T, z, z, z],
            [z, T, z, z],
            [z, z, T, z],
            [z, z, z, T]
        ])

        Ke = _beami_stiffness(L, pid_ref, mat, I1, I2)
        K = Teb.T @ Ke @ Teb
        dofs = [
            (nid1, 1), (nid1, 2), (nid1, 3),
            (nid1, 4), (nid1, 5), (nid1, 6),

            (nid2, 1), (nid2, 2), (nid2, 3),
            (nid2, 4), (nid2, 5), (nid2, 6),
        ]
        n_ijv = [
            dof_map[(nid1, 1)], dof_map[(nid1, 2)], dof_map[(nid1, 3)],
            dof_map[(nid1, 4)], dof_map[(nid1, 5)], dof_map[(nid1, 6)],

            dof_map[(nid2, 1)], dof_map[(nid2, 2)], dof_map[(nid2, 3)],
            dof_map[(nid2, 4)], dof_map[(nid2, 5)], dof_map[(nid2, 6)],
        ]
        for unused_dof1, i1 in zip(dofs, n_ijv):
            for unused_dof2, i2 in zip(dofs, n_ijv):
                ki = K[i1, i2]
                if abs(ki) > 0.:
                    Kbb[i1, i2] = ki
                    Kbbs[i1, i2] = ki
    return nelements

def _build_kbb_crod(model, Kbb, Kbbs, dof_map):
    """fill the CROD Kbb matrix"""
    eids = model._type_to_id_map['CROD']
    for eid in eids:
        elem = model.elements[eid]
        pid_ref = elem.pid_ref
        mat = pid_ref.mid_ref
        _build_kbbi_conrod_crod(Kbb, Kbbs, dof_map, elem, mat)
    return len(eids)

def _build_kbb_ctube(model: BDF, Kbb, Kbbs, dof_map):
    """fill the CTUBE Kbb matrix"""
    ctubes = model._type_to_id_map['CTUBE']
    for eid in ctubes:
        elem = model.elements[eid]
        pid_ref = elem.pid_ref
        mat = pid_ref.mid_ref
        _build_kbbi_conrod_crod(Kbb, Kbbs, dof_map, elem, mat)
    return len(ctubes)

def _build_kbb_conrod(model: BDF, Kbb, Kbbs, dof_map):
    """fill the CONROD Kbb matrix"""
    eids = model._type_to_id_map['CONROD']
    for eid in eids:
        elem = model.elements[eid]
        mat = elem.mid_ref
        _build_kbbi_conrod_crod(Kbb, Kbbs, dof_map, elem, mat)
    return len(eids)

def _build_kbbi_conrod_crod(Kbb, Kbbs, dof_map, elem, mat, fdtype='float64'):
    """fill the ith rod Kbb matrix"""
    nid1, nid2 = elem.nodes
    #mat = elem.mid_ref
    xyz1 = elem.nodes_ref[0].get_position()
    xyz2 = elem.nodes_ref[1].get_position()
    dxyz12 = xyz1 - xyz2
    L = np.linalg.norm(dxyz12)
    E = mat.E
    G = mat.G()
    J = elem.J()
    A = elem.Area()
    E = elem.E()
    #L = elem.Length()
    k_axial = A * E / L
    k_torsion = G * J / L

    assert isinstance(k_axial, float), k_axial
    assert isinstance(k_torsion, float), k_torsion
    #Kbb[i, i] += ki[0, 0]
    #Kbb[i, j] += ki[0, 1]
    #Kbb[j, i] = ki[1, 0]
    #Kbb[j, j] = ki[1, 1]
    k = np.array([[1., -1.],
                  [-1., 1.]])  # 1D rod
    Lambda = lambda1d(dxyz12, debug=False)
    K = Lambda.T @ k @ Lambda
    #i11 = dof_map[(n1, 1)]
    #i12 = dof_map[(n1, 2)]

    #i21 = dof_map[(n2, 1)]
    #i22 = dof_map[(n2, 2)]

    nki, nkj = K.shape
    K2 = np.zeros((nki*2, nkj*2), dtype=fdtype)

    i1 = 0
    i2 = 3 # dof_map[(n1, 2)]
    if k_torsion == 0.0: # axial; 2D or 3D
        K2 = K * k_axial
        n_ijv = [
            # axial
            dof_map[(nid1, 1)], dof_map[(nid1, 2)], dof_map[(nid1, 3)],
            dof_map[(nid2, 1)], dof_map[(nid2, 2)], dof_map[(nid2, 3)],
        ]
        dofs = np.array([
            i1, i1+1, i1+2,
            i2, i2+1, i2+2,
        ], dtype='int32')
    elif k_axial == 0.0: # torsion; assume 3D
        K2 = K * k_torsion
        n_ijv = [
            # torsion
            dof_map[(nid1, 4)], dof_map[(nid1, 5)], dof_map[(nid1, 6)],
            dof_map[(nid2, 4)], dof_map[(nid2, 5)], dof_map[(nid2, 6)],
        ]
        dofs = np.array([
            i1, i1+1, i1+2,
            i2, i2+1, i2+2,
        ], dtype='int32')
    else:  # axial + torsion; assume 3D
        # u1fx, u1fy, u1fz, u2fx, u2fy, u2fz
        K2[:nki, :nki] = K * k_axial

        # u1mx, u1my, u1mz, u2mx, u2my, u2mz
        K2[nki:, nki:] = K * k_torsion

        dofs = np.array([
            i1, i1+1, i1+2,
            i2, i2+1, i2+2,

            i1+3, i1+4, i1+5,
            i2+3, i2+4, i2+5,
        ], dtype='int32')
        n_ijv = [
            # axial
            dof_map[(nid1, 1)], dof_map[(nid1, 2)], dof_map[(nid1, 3)],
            dof_map[(nid2, 1)], dof_map[(nid2, 2)], dof_map[(nid2, 3)],

            # torsion
            dof_map[(nid1, 4)], dof_map[(nid1, 5)], dof_map[(nid1, 6)],
            dof_map[(nid2, 4)], dof_map[(nid2, 5)], dof_map[(nid2, 6)],
        ]
    for dof1, i1 in zip(dofs, n_ijv):
        for dof2, i2 in zip(dofs, n_ijv):
            ki = K2[dof1, dof2]
            if abs(ki) > 0.:
                #print(nij1, nij2, f'({i1}, {i2});', (dof1, dof2), ki)
                Kbb[i1, i2] = ki
                Kbbs[i1, i2] = ki
        #print(K2)
    #print(Kbb)
    return

def _build_kbb_cquad4(model: BDF, Kbb, Kbbs, dof_map,
                      all_nids, xyz_cid0, idtype='int32', fdtype='float64'):
    """fill the CQUAD4 Kbb matrix

    https://www.sciencedirect.com/topics/engineering/quadrilateral-element
    """
    eids = np.array(model._type_to_id_map['CQUAD4'], dtype=idtype)
    nelements = len(eids)
    if nelements == 0:
        return nelements

    eids.sort()
    ncoords = len(model.coords)
    cids = np.zeros(ncoords, dtype='int32')
    coords = np.zeros((ncoords, 3, 3), dtype='float64')
    iaxes = np.zeros((ncoords, 3), dtype='float64')
    for icid, (cid, coord) in zip(count(), sorted(model.coords.items())):
        cids[icid] = cid
        iaxes[icid, :] = coord.i
        coords[icid, :, :] = coord.beta()
        icid += 1

    theta_mcid = []
    normal = []

    nids = np.zeros((nelements, 4), dtype='int32')
    for i, eid in enumerate(eids):
        elem = model.elements[eid]  # type: CQUAD4
        theta_mcidi = elem.theta_mcid
        nids[i, :] = elem.nodes
        # nids.append(elem.nodes)
        theta_mcid.append(theta_mcidi)

    inids = np.searchsorted(all_nids, nids.ravel()).reshape(nelements, 4)
    p1 = xyz_cid0[inids[:, 0], :]
    p2 = xyz_cid0[inids[:, 1], :]
    p3 = xyz_cid0[inids[:, 2], :]
    p4 = xyz_cid0[inids[:, 3], :]


    # normal is correct; matters for +rotation and offsets
    v13 = p1 - p3
    v24 = p2 - p4
    normal = np.cross(v13, v24)
    ni = np.linalg.norm(normal, axis=1)
    if np.allclose(ni, 0.):
        ibad = np.where(ni == 0.)[0]
        bad_eids = eids[ibad]
        raise RuntimeError(f'elements={bad_eids.tolist()} have invalid normal vectors')
    normal /= ni[:, np.newaxis]
    area = 0.5 * ni
    assert len(area) == nelements, 'len(area)=%s nelements=%s' % (len(area), nelements)
    #centroid = (p1 + p2 + p3 + p4) / 4.

    #inids = np.searchsorted(all_nids, nids.ravel()).reshape(nelementsi, 2)
    #p1 = xyz_cid0[inids[:, 0], :]
    #p2 = xyz_cid0[inids[:, 1], :]
    iaxes = p2 - p1
    #length = np.linalg.norm(iaxes, axis=1)
    #centroid = (p1 + p2) / 2.


    # e2i is (npids,3,3)
    # e2 is (nelementsi,3,3)
    #e2i = np.array(e2_dict['shell'], dtype='float64')
    #e2 = e2i[ipids, :, :]
    #assert e2.shape[0] == nelementsi

    # [T^T][e][T]
    #theta_mcid = theta_mcid_dict[etype]
    #telem = breakdown_material_coordinate_system(cids, iaxes, theta_mcid, normal, p1, p2)
    #exx = eyy = ezz = 1.
    #exx, eyy, ezz = _breakdown_transform_shell(e2, telem)
    #mpa = mass_per_area[ipids]
    #npa = nsm_per_area[ipids]
    #mass = mpa * area
    #nsm = npa * area

    # assume the panel is square to calculate w; then multiply by t to get tw
    #assert len(area) > 0, area
    #assert len(thickness) > 0, thickness
    #tw = thickness[ipids] * np.sqrt(area)

#if etype in ['CQUAD4', 'CTRIA3', 'CQUAD8', 'CTRIA6', 'CTRIAR', 'CQUADR', 'CQUAD', 'CTRIAX']:
    #nids_dict[etype].append(elem.nodes)
    #pids_dict[etype].append(elem.pid)
    #thetai = elem.theta_mcid
    #theta_mcid_dict[etype].append(thetai)

    #print(elem, ki)
    #print(elem.get_stats())
    #_build_kbbi_celas12(Kbb, Kbbs, dof_map, elem, ki)
    T = transform_shell_material_coordinate_system(cids, iaxes, theta_mcid, normal, p1, p2,
                                                   idtype=idtype, fdtype=fdtype)
    # tet = np.einsum('nij,njk->nik', telem, et)

    bad_jacobians = []
    for eid, p1i, p2i, p3i, p4i, Ti in zip(eids, p1, p2, p3, p4, T):
        # (4,3) = (4,3) x (3,3)
        xy1 = Ti @ p1i
        xy2 = Ti @ p2i
        xy3 = Ti @ p3i
        xy4 = Ti @ p4i
        #https://math.stackexchange.com/questions/2430691/jacobian-determinant-for-bi-linear-quadrilaterals
        # x, zeta direction = 1 - 2
        # y, eta direction =  2 - 3
        print(xy1)
        print(xy2)
        print(xy3)
        print(xy4)
        x1 = xy1[0]
        x2 = xy2[0]
        x3 = xy3[0]
        x4 = xy4[0]

        y1 = xy1[1]
        y2 = xy2[1]
        y3 = xy3[1]
        y4 = xy4[1]

        #    ^ eta, y
        #    |
        #    |
        # 4-----3
        # |     |
        # |     |---> zeta, x
        # |     |
        # 1-----2
        dx_deta12 = (x2 - x1) / 2.
        dx_deta34 = (x3 - x4) / 2.
        dx_deta = (dx_deta12 + dx_deta34) / 2.

        dy_deta12 = (y2 - y1) / 2.
        dy_deta34 = (y3 - y4) / 2.
        dy_deta = (dy_deta12 + dy_deta34) / 2.

        dx_dzeta14 = (x4 - x1) / 2.
        dx_dzeta23 = (x3 - x2) / 2.
        dx_dzeta = (dx_dzeta14 + dx_dzeta23) / 2.

        dy_dzeta14 = (y4 - y1) / 2.
        dy_dzeta23 = (y3 - y2) / 2.
        dy_dzeta = (dy_dzeta14 + dy_dzeta23) / 2.
        jmat = np.array([
            [dx_deta, dy_deta],
            [dx_dzeta, dy_dzeta],
        ])
        # [du_dzeta]  = [dx_dzeta, dy_dzeta] [du_dx]
        # [du_deta ]    [dx_ zeta, dy_deta ] [du_dy]
        jacobian = np.linalg.det(jmat)
        # K = [N]^T[C][N] * |J|
        #   where C = [A], 2[B], [D] matrices
        #N = [
            #[1., x1, y1, x1*y1, 0., 0., 0., 0.],
            #[0., 0., 0., 0., 1., x2, y2, x2*y2],
        #]
        elem = model.elements[eid]
        nid1, nid2, nid3, nid4 = elem.nodes
        pid_ref = elem.pid_ref
        #materials = []
        ptype = pid_ref.type
        if ptype == 'PSHELL':
            A, B, D = pid_ref.get_individual_ABD_matrices()
            str((A, B, D))
        elif ptype == 'PCOMP':
            A, B, D = pid_ref.get_individual_ABD_matrices()
            str((A, B, D))
        else:
            raise NotImplementedError(pid_ref)

        i1 = dof_map[(nid1, 1)]
        i2 = dof_map[(nid2, 1)]
        i3 = dof_map[(nid3, 1)]
        i4 = dof_map[(nid4, 1)]
        dofs = np.array([
            i1, i1+1,
            i2, i2+1,
            i3, i3+1,
            i4, i4+1,
            #i1+3, i1+4,
            #i2+3, i2+4,
            #i3+3, i3+4,
            #i4+3, i4+4,
        ], dtype='int32')
        n_ijv = [
            # axial
            dof_map[(nid1, 1)], dof_map[(nid1, 2)],
            dof_map[(nid2, 1)], dof_map[(nid2, 2)],
            dof_map[(nid3, 1)], dof_map[(nid3, 2)],
            dof_map[(nid4, 1)], dof_map[(nid4, 2)],

            # torsion
            #dof_map[(nid1, 4)], dof_map[(nid1, 5)]],
            #dof_map[(nid2, 4)], dof_map[(nid2, 5)]],
            #dof_map[(nid3, 4)], dof_map[(nid3, 5)]],
            #dof_map[(nid4, 4)], dof_map[(nid4, 5)]],
        ]
        Ki = A
        #for C in [A, 2*B, D]:
        for C in [2*B, D]:
            Ki += C
            # K += (N.T @ C @ N) * jacobian
        Ki *= jacobian
        model.log.debug(f'Ki:\n{Ki}')
        for unused_dof1, i1 in zip(dofs, n_ijv):
            for unused_dof2, i2 in zip(dofs, n_ijv):
                pass
                # ki = K2[dof1, dof2]
                #Kbbs[i1, i2] += cj
        #print(xy1, xy2, xy3, xy4)

        # TODO: The jacobian ratio is the ratio between the min/max values of the
        #       jacobians for the 4 gauss points.
        #       This is a bandaid...
        jacobian2 = np.linalg.det(jmat / np.abs(jmat).max())
        jratio = jacobian2
        if abs(jacobian) > 1:
            jratio = 1 / jacobian2

        #jratio = jacobians.min() / jacobians.max()
        if not(0.1 <= jratio <= 10.):
            model.log.error(f'eid={eid}; |J|={jacobian:.3f}; Jratio={jratio:.3f} J=\n{jmat}')
            bad_jacobians.append(eid)

    if bad_jacobians:
        raise RuntimeError(f'elements={bad_jacobians} have invalid jacobians')
    return nelements

def _build_kbb_cbeam(model: BDF, Kbb, Kbbs, dof_map,
                     all_nids, xyz_cid0, idtype='int32', fdtype='float64'):
    """TODO: Timoshenko beam"""
    str(all_nids)
    str(xyz_cid0)
    eids = np.array(model._type_to_id_map['CBEAM'], dtype=idtype)
    nelements = len(eids)
    if nelements == 0:
        return nelements

    for eid in eids:
        elem = model.elements[eid]
        nid1, nid2 = elem.nodes
        xyz1 = elem.nodes_ref[0].get_position()
        xyz2 = elem.nodes_ref[1].get_position()
        dxyz = xyz2 - xyz1
        L = np.linalg.norm(dxyz)
        pid_ref = elem.pid_ref
        mat = pid_ref.mid_ref
        is_passed, (wa, wb, ihat, jhat, khat) = elem.get_axes(model)
        T = np.vstack([ihat, jhat, khat])
        z = np.zeros((3, 3), dtype=fdtype)
        Teb = np.block([
            [T, z, z, z],
            [z, T, z, z],
            [z, z, T, z],
            [z, z, z, T]
        ])
        Iy = pid_ref.I11()
        Iz = pid_ref.I22()
        Ke = _beami_stiffness(L, pid_ref, mat, Iy, Iz)
        K = Teb.T @ Ke @ Teb
        dofs = [
            (nid1, 1), (nid1, 2), (nid1, 3),
            (nid1, 4), (nid1, 5), (nid1, 6),

            (nid2, 1), (nid2, 2), (nid2, 3),
            (nid2, 4), (nid2, 5), (nid2, 6),
        ]
        n_ijv = [
            dof_map[(nid1, 1)], dof_map[(nid1, 2)], dof_map[(nid1, 3)],
            dof_map[(nid1, 4)], dof_map[(nid1, 5)], dof_map[(nid1, 6)],

            dof_map[(nid2, 1)], dof_map[(nid2, 2)], dof_map[(nid2, 3)],
            dof_map[(nid2, 4)], dof_map[(nid2, 5)], dof_map[(nid2, 6)],
        ]
        for unused_dof1, i1 in zip(dofs, n_ijv):
            for unused_dof2, i2 in zip(dofs, n_ijv):
                ki = K[i1, i2]
                if abs(ki) > 0.:
                    Kbb[i1, i2] = ki
                    Kbbs[i1, i2] = ki
    return nelements

def _beami_stiffness(L, prop, mat, Iy, Iz):
    """gets the ith Euler-Bernoulli beam stiffness"""
    E = mat.E()
    G = mat.G()
    A = prop.Area()
    J = prop.J()

    kaxial = E * A / L
    ktorsion = G * J / L
    ky = E * Iy / L
    kz = E * Iz / L

    K = np.zeros((12, 12), dtype='float64')
    # axial
    K[0, 0] = K[6, 6] = kaxial
    K[6, 0] = K[0, 6] = -kaxial

    # torsion
    K[3, 3] = K[9, 9] = ktorsion
    K[9, 3] = K[3, 9] = -ktorsion

    #Fx - 0, 6
    #Fy - 1, 7**
    #Fz - 2, 8
    #Mx - 3, 9
    #My - 4, 10
    #Mz - 5, 11**
    # bending z (Fy/1/7, Mz/5/11)
    #      1     5       7   11
    # 1  [12  & 6L   & -12 & 6L
    # 5  [6L  & 4L^2 & -6L & 2L^2
    # 7  [-12 &-6L   &  12 & -6L
    # 11 [6L  & 2L^2 & -6L & 4L^2
    K[1, 1] = K[7, 7] = 12 * kz
    K[1, 7] = K[1, 7] = -12 * kz
    K[1, 5] = K[5, 1] = K[11, 1] = K[1, 11] = 6 * L * kz

    K[5, 7] = K[7, 5] = K[7, 11] = K[11, 7] = -6 * L * kz
    K[5, 11] = K[11, 5] = 2 * L * L * kz
    K[5, 5] = K[11, 11] = 4 * L * L * kz

    #Fx - 0, 6
    #Fy - 1, 7
    #Fz - 2, 8**
    #Mx - 3, 9
    #My - 4, 10**
    #Mz - 5, 11
    # bending y (Fz/2/8, My/4/10)
    #      2     4       8   10
    # 2  [12  & 6L   & -12 & 6L
    # 4  [6L  & 4L^2 & -6L & 2L^2
    # 8  [-12 &-6L   &  12 & -6L
    # 10 [6L  & 2L^2 & -6L & 4L^2
    K[2, 2] = K[8, 8] = 12 * ky
    K[2, 8] = K[2, 8] = -12 * ky
    K[2, 4] = K[4, 2] = K[10, 2] = K[2, 10] = 6 * L * ky

    K[4, 8] = K[8, 4] = K[8, 10] = K[10, 8] = -6 * L * ky
    K[4, 10] = K[10, 4] = 2 * L * L * ky
    K[4, 4] = K[10, 10] = 4 * L * L * ky
    return K
