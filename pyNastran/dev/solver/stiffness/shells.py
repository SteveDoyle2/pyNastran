from __future__ import annotations
from itertools import count
from typing import TYPE_CHECKING

import numpy as np
#import scipy.sparse as sci_sparse

from pyNastran.bdf.cards.elements.shell import transform_shell_material_coordinate_system
from ..utils import DOF_MAP
#from pyNastran.bdf.cards.elements.bars import get_bar_vector, get_bar_yz_transform
if TYPE_CHECKING:  # pragma: no cover
    from pyNastran.nptyping_interface import NDArrayN3float, NDArrayNNfloat
    from pyNastran.bdf.bdf import (
        BDF,
        CQUAD4, CQUAD8, MAT1
    )
    #from pyNastran.bdf.cards.elements.shell import CQUAD4

def build_kbb_cquad4(model: BDF,
                     Kbb,
                     dof_map: DOF_MAP,
                     all_nids, xyz_cid0: NDArrayN3float, idtype='int32', fdtype='float64') -> int:
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
    bad_jacobians2 = []
    sqrt3 = 1 / np.sqrt(3)
    zs_etas = [(-sqrt3, -sqrt3), (sqrt3, -sqrt3), (-sqrt3, sqrt3), (sqrt3, sqrt3)]

    for eid, p1i, p2i, p3i, p4i, Ti in zip(eids, p1, p2, p3, p4, T):
        # (4,3) = (4,3) x (3,3)
        xy1 = Ti @ p1i
        xy2 = Ti @ p2i
        xy3 = Ti @ p3i
        xy4 = Ti @ p4i
        #https://math.stackexchange.com/questions/2430691/jacobian-determinant-for-bi-linear-quadrilaterals
        # x, zeta direction = 1 - 2
        # y, eta direction =  2 - 3
        #print(xy1)
        #print(xy2)
        #print(xy3)
        #print(xy4)
        x1 = xy1[0]
        x2 = xy2[0]
        x3 = xy3[0]
        x4 = xy4[0]

        y1 = xy1[1]
        y2 = xy2[1]
        y3 = xy3[1]
        y4 = xy4[1]

        A0 = ((y4 - y2) * (x3 - x1) - (y3 - y1) * (x4 - x2)) / 8
        A1 = ((y3 - y4) * (x2 - x1) - (y2 - y1) * (x3 - x4)) / 8
        A2 = ((y4 - y1) * (x3 - x2) - (y3 - y2) * (x4 - x1)) / 8

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

        elem = model.elements[eid]
        nid1, nid2, nid3, nid4 = elem.nodes
        pid_ref = elem.pid_ref
        #materials = []
        ptype = pid_ref.type
        if ptype == 'PSHELL':
            A, Bmat, D = pid_ref.get_individual_ABD_matrices()
            str((A, Bmat, D))
        elif ptype == 'PCOMP':
            A, Bmat, D = pid_ref.get_individual_ABD_matrices()
            str((A, Bmat, D))
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
        ], dtype='int32')
        n_ijv = [
            dof_map[(nid1, 1)], dof_map[(nid1, 2)],
            dof_map[(nid2, 1)], dof_map[(nid2, 2)],
            dof_map[(nid3, 1)], dof_map[(nid3, 2)],
            dof_map[(nid4, 1)], dof_map[(nid4, 2)],
        ]

        jacobian2 = []
        Ki = np.zeros((8, 8), dtype='float64')
        for zi, etai in zs_etas:
            jacobian2i = A0 + A1 * zi + A2 * etai
            jacobian2.append(jacobian2i)
            N1x = N2x = etai - 1
            N3x = N4x = etai + 1
            N1y = N4y = zi - 1
            N2y = N3y = zi + 1
            B = np.array([
                [N1x, 0, N2x, 0, N3x, 0, N4x, 0],
                [0, N1y, 0, N2y, 0, N3y, 0, N4y],
                [N1y, N1x, N2y, N2x, N3y, N3x, N4y, N4x],
            ])

            # K = [B]^T[C][B] * |J|
            #   where C = [A], 2[B], [D] matrices
            #N = [
                #[1., x1, y1, x1*y1, 0., 0., 0., 0.],
                #[0., 0., 0., 0., 1., x2, y2, x2*y2],
            #]
            #model.log.debug(f'B {B.shape}:\n{B}')
            Ki += (B.T @ A @ B)
            for C in [2*Bmat, D]:
                Ki += (B.T @ C @ B)
            Ki *= jacobian2i

        if np.abs(Ki).sum() == 0.0:
            if pid_ref.type == 'PSHELL':
                model.log.error(f'K=0; eid={eid} ptype={pid_ref.type} mid1={pid_ref.mid1} mid2={pid_ref.mid2} '
                                f'mid3={pid_ref.mid3} mid4={pid_ref.mid4}')
            else:
                model.log.error(f'K=0; eid={eid} ptype={pid_ref.type} mids={pid_ref.mids}')
            #model.log.debug(f'A {A.shape}:\n{A}')
            #model.log.debug(f'B {Bmat.shape}:\n{Bmat}')
            #model.log.debug(f'D {D.shape}:\n{D}')
            continue

        #model.log.debug(f'Ki {Ki.shape}:\n{Ki}')
        for idof, unused_dof1, i1 in zip(count(), dofs, n_ijv):
            for jdof, unused_dof2, i2 in zip(count(), dofs, n_ijv):
                ki = Ki[idof, jdof]
                Kbb[i1, i2] += ki
        #print(xy1, xy2, xy3, xy4)

        # TODO: The jacobian ratio is the ratio between the min/max values of the
        #       jacobians for the 4 gauss points.
        #       This is a bandaid...
        jacobian3 = np.linalg.det(jmat / np.abs(jmat).max())
        jratio = jacobian3
        #if abs(jacobian) > 1:
            #jratio = 1 / jacobian2

        #jratio = jacobians.min() / jacobians.max()
        jratio2 = max(jacobian2) / min(jacobian2)
        if not(0.1 <= jratio <= 10.):
            model.log.error(f'eid={eid}; |J|={jacobian:.3f}; |J2|={jacobian2}; Jratio={jratio2:.3f} J=\n{jmat}')
            bad_jacobians.append(eid)
        else:
            model.log.debug(f'eid={eid}; |J|={jacobian:.3f}; |J2|={jacobian2}; Jratio={jratio2:.3f} J=\n{jmat}')

    if bad_jacobians:
        raise RuntimeError(f'elements={bad_jacobians} have invalid jacobians')
    return nelements

def build_kbb_cquad8(model: BDF,
                     Kbb,
                     dof_map: DOF_MAP,
                     all_nids, xyz_cid0: NDArrayN3float, idtype='int32', fdtype='float64') -> int:
    """fill the CQUAD8 Kbb matrix

    https://www.sciencedirect.com/topics/engineering/quadrilateral-element
    """
    eids = np.array(model._type_to_id_map['CQUAD8'], dtype=idtype)
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

    nids = np.zeros((nelements, 8), dtype='int32')
    for i, eid in enumerate(eids):
        elem = model.elements[eid]  # type: CQUAD4
        theta_mcidi = elem.theta_mcid
        nids[i, :] = elem.nodes
        # nids.append(elem.nodes)
        theta_mcid.append(theta_mcidi)

    inids = np.searchsorted(all_nids, nids.ravel()).reshape(nelements, 8)
    p1 = xyz_cid0[inids[:, 0], :]
    p2 = xyz_cid0[inids[:, 1], :]
    p3 = xyz_cid0[inids[:, 2], :]
    p4 = xyz_cid0[inids[:, 3], :]
    p5 = xyz_cid0[inids[:, 4], :]
    p6 = xyz_cid0[inids[:, 5], :]
    p7 = xyz_cid0[inids[:, 6], :]
    p8 = xyz_cid0[inids[:, 7], :]


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

    iaxes = p2 - p1
    T = transform_shell_material_coordinate_system(cids, iaxes, theta_mcid, normal, p1, p2,
                                                   idtype=idtype, fdtype=fdtype)
    # tet = np.einsum('nij,njk->nik', telem, et)

    bad_jacobians = []
    for eid, p1i, p2i, p3i, p4i, p5i, p6i, p7i, p8i, Ti in zip(eids, p1, p2, p3, p4, p5, p6, p7, p8, T):
        # (4,3) = (4,3) x (3,3)
        xy1 = Ti @ p1i
        xy2 = Ti @ p2i
        xy3 = Ti @ p3i
        xy4 = Ti @ p4i

        xy5 = Ti @ p5i
        xy6 = Ti @ p6i
        xy7 = Ti @ p7i
        xy8 = Ti @ p8i
        #https://math.stackexchange.com/questions/2430691/jacobian-determinant-for-bi-linear-quadrilaterals
        # x, zeta direction = 1 - 2
        # y, eta direction =  2 - 3
        #print(xy1)
        #print(xy2)
        #print(xy3)
        #print(xy4)
        x1 = xy1[0]
        x2 = xy2[0]
        x3 = xy3[0]
        x4 = xy4[0]
        x5 = xy5[0]
        x6 = xy6[0]
        x7 = xy7[0]
        x8 = xy8[0]

        y1 = xy1[1]
        y2 = xy2[1]
        y3 = xy3[1]
        y4 = xy4[1]
        y5 = xy5[1]
        y6 = xy6[1]
        y7 = xy7[1]
        y8 = xy8[1]

        #    ^ eta, y
        #    |
        #    |
        # 4--7--3
        # |     |
        # 8     6---> zeta, x
        # |     |
        # 1--5--2

        # file:///C:/Users/sdoyle/Downloads/FEM_1_9_8node_2D.pdf
        xi = [x1, x2, x3, x4, x5, x6, x7, x8]
        yi = [y1, y2, y3, y4, y5, y6, y7, y8]
        s3 = 1 / np.sqrt(3)
        xyi = [(-s3, -s3), (s3, -s3), (-s3, s3), (s3, s3), ]  # fixme
        for x, y in xyi:
            #print(x, y)
            N1x = (2 * x + y) * (y - 1) / 4
            N1y = (x - 1) * (x + 2 * y) / 4

            N2x = x * (y - 1)
            N2y = x**2 - 1 / 2

            N3x = (2 * x - y) * (y - 1) / 4
            N3y = (x + 1) * (x - 2 * y) / 4

            N4x = 1/2 - y**2
            N4y = -y * (x + 1)

            N5x = -1 * (2 * x + y) * (y + 1) / 4
            N5y = -1 * (x + 1) * (x + 2 * y) / 4

            N6x = -x * (y+1)
            N6y = 1 / 2 - x**2

            N7x = (-2 * x + y) * (y + 1) / 4
            N7y = (-x + 2 * y) * (x - 1) / 4

            N8x = y ** 2 - 1 / 2
            N8y = y * (x - 1)
            #N = np.array([
                #[N1, 0, N2, 0, N3, 0, N4, 0, N5, 0, N6, 0, N7, 0, N8, 0],
                #[0, N1, 0, N2, 0, N3, 0, N4, 0, N5, 0, N6, 0, N7, 0, N8],
            #])
            B = np.array([
                [N1x, 0, N2x, 0, N3x, 0, N4x, 0, N5x, 0, N6x, 0, N7x, 0, N8x, 0],
                [0, N1y, 0, N2y, 0, N3y, 0, N4y, 0, N5y, 0, N6y, 0, N7y, 0, N8y],
                [N1y, N1x, N2y, N2x, N3y, N3x, N4y, N4x, N5y, N5x, N6y, N6x, N7y, N7x, N8y, N8x],
            ])
            #Ni = [N1, N2, N3, N4, N5, N6, N7, N8]
            Nxi = [N1x, N2x, N3x, N4x, N5x, N6x, N7x, N8x]
            Nyi = [N1y, N2y, N3y, N4y, N5y, N6y, N7y, N8y]

            nxsum = 0.
            nysum = 0.
            for Nx, Ny in zip(Nxi, Nyi):
                jmat = np.array([
                    [Nx * x, Nx * y],
                    [Ny * x, Ny * y],
                ])
                #jmat_inv = np.linalg.inv(jmat)
                jacobian = np.linalg.det(jmat)
                #model.log.error(f'eid={eid}; |J|={jacobian:.3f} jmat={jmat}')
                nxsum += Nx
                nysum += Ny
            #print(f'nxsum={nxsum} nysum={nysum}')
        # K = [N]^T[C][N] * |J|
        #   where C = [A], 2[B], [D] matrices
        #N = [
            #[1., x1, y1, x1*y1, 0., 0., 0., 0.],
            #[0., 0., 0., 0., 1., x2, y2, x2*y2],
        #]
        elem = model.elements[eid]
        nid1, nid2, nid3, nid4, nid5, nid6, nid7, nid8 = elem.nodes
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
            dof_map[(nid1, 1)], dof_map[(nid1, 2)],
            dof_map[(nid2, 1)], dof_map[(nid2, 2)],
            dof_map[(nid3, 1)], dof_map[(nid3, 2)],
            dof_map[(nid4, 1)], dof_map[(nid4, 2)],
            dof_map[(nid5, 1)], dof_map[(nid5, 2)],
            dof_map[(nid6, 1)], dof_map[(nid6, 2)],
            dof_map[(nid7, 1)], dof_map[(nid7, 2)],
            dof_map[(nid8, 1)], dof_map[(nid8, 2)],
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
