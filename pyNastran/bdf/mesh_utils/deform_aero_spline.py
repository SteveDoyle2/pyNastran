import os
import copy
from pathlib import Path
from itertools import count
from typing import Optional
import numpy as np
from pyNastran.utils import PathLike
from pyNastran.bdf.bdf import read_bdf, BDF
from pyNastran.op2.op2 import read_op2, OP2


def _get_nids_xyz_cid0(model: BDF,
                       nids: Optional[np.ndarray]=None,
                       xyz_cid0: Optional[np.ndarray]=None) -> tuple[np.ndarray, np.ndarray]:
    """helper method"""
    if model is not None:
        _npoints, _nids, nids = model._get_npoints_nids_allnids()
        nids = np.array(nids)
        nids.sort()
        del _npoints, _nids
        xyz_cid0 = model.get_xyz_in_coord()
    else:
        assert nids is not None, nids
        assert xyz_cid0 is not None, xyz_cid0
    assert len(nids) == len(xyz_cid0)
    return nids, xyz_cid0


def _get_discrete_aero_mesh(structure_model: Optional[BDF],
                            aero_model: Optional[BDF]) -> tuple[BDF, np.ndarray, int]:
    """
    A discrete aero_mesh is one that has unique nodes for every element
    regardless of the connectivity. This simplifies deflection mapping.
    """
    if aero_model is None:
        pid = 1
        eid = 1
        mid = 1
        t = 0.1
        aero_model = BDF()
        aero_model.add_pshell(pid, mid, t)
        aero_model.add_mat1(mid, 3.0e7, None, 0.3)

        naero_boxes_all = 0
        aero_xyz_cid0_point_list = []
        nid0 = 1
        for caero_id, caero in sorted(structure_model.caeros.items()):
            caero_type = caero.type
            if caero_type == 'CAERO1':
                points, elements = caero.panel_points_elements()
                naero_boxes_all += len(elements)
                for elem in elements:
                    p1 = points[elem[0], :]
                    p2 = points[elem[1], :]
                    p3 = points[elem[2], :]
                    p4 = points[elem[3], :]
                    aero_xyz_cid0_point_list.append(p1)
                    aero_xyz_cid0_point_list.append(p2)
                    aero_xyz_cid0_point_list.append(p3)
                    aero_xyz_cid0_point_list.append(p4)
                    aero_model.add_cquad4(eid, pid, nid0+elem)
                    eid += 1
            #elif caero_type == 'CAERO2':
            else:
                raise NotImplementedError(caero)
            nid0 += len(points)
        aero_xyz_cid0_point = np.array(aero_xyz_cid0_point_list, dtype='float64')
        naero_nodes_all = naero_boxes_all * 4
        assert len(aero_xyz_cid0_point) == naero_nodes_all

        for inid, xyz in enumerate(aero_xyz_cid0_point):
            aero_model.add_grid(inid+1, xyz)
    else:
        # get the entire caero box mesh
        aero_xyz_cid0 = aero_model.get_xyz_in_coord()
        # aero_disp_out_cid0 = np.zeros(aero_xyz_cid0.shape, dtype='float64')
        naero_nodes_all = len(aero_xyz_cid0)
        naero_boxes_all = len(aero_model.elements)
        assert naero_nodes_all > 0, naero_nodes_all
        assert naero_boxes_all > 0, naero_boxes_all
        # del naero_nodes_all

        # split aero model into individual panels.
        # Makes the deflections show up on the individual panels
        # instead of linking them.
        # aero_xyz_cid0_point = np.zeros((naero_boxes_all*4, 3))
        aero_xyz_cid0_point_list = []
        for eid, elem in aero_model.elements.items():
            for node in elem.nodes_ref:
                xyzi = node.get_position()
                aero_xyz_cid0_point_list.append(xyzi)
        aero_xyz_cid0_point = np.array(aero_xyz_cid0_point_list, dtype='float64')

    naero_nodes = len(aero_model.nodes)
    assert naero_nodes > 0, naero_nodes
    return aero_model, aero_xyz_cid0_point, naero_boxes_all


def _get_all_aero_boxs(model: BDF) -> np.ndarray:
    # TODO: assumes splines and caeros are sorted in the same way
    #       what does this impact if that's not true?
    all_aero_box_list = []
    for eid, spline in model.splines.items():
        boxs = spline.aero_element_ids
        all_aero_box_list.append(boxs)
    all_aero_box = np.hstack(all_aero_box_list)
    all_aero_box.sort()
    return all_aero_box


def deform_aero_spline_from_files(bdf_filename: PathLike,
                                  op2_filename: PathLike,
                                  bdf_filename_out: PathLike='',
                                  op2_filename_out: PathLike=''):
    dirname = Path(os.path.dirname(bdf_filename))
    if bdf_filename_out in {None, ''}:
        bdf_filename_out = dirname / 'aero_model_splined.bdf'
    if op2_filename_out in {None, ''}:
        op2_filename_out = dirname / 'aero_model_splined.op2'

    structure_model = read_bdf(bdf_filename, debug=False)
    nids, xyz_cid0 = _get_nids_xyz_cid0(
        structure_model, nids=None, xyz_cid0=None)
    aero_model, aero_xyz_cid0_point, naero_boxes_all = _get_discrete_aero_mesh(
        structure_model, aero_model=None)
    naero_node = len(aero_model.nodes)
    assert naero_node > 0, naero_node
    aero_model.write_bdf(bdf_filename_out)
    # -------------------------------------------------
    results_model = read_op2(op2_filename, debug=False)
    results_model_out = OP2(debug=False, mode='msc')

    for case_key, disp in results_model.displacements.items():
        case = copy.deepcopy(disp)
        _fill_disp_case(structure_model, disp, case, naero_node)
        results_model_out.displacements[case_key] = case

    for case_key, eigenvectors in results_model.eigenvectors.items():
        case = copy.deepcopy(eigenvectors)
        _fill_disp_case(structure_model, eigenvectors, case, naero_node)
        results_model_out.eigenvectors[case_key] = case

    # results_model_out._nastran_format = 'msc'
    results_model_out.write_op2(op2_filename_out,
                                nastran_format='msc')
    return

def _fill_disp_case(structure_model: BDF,
                    disp, case,
                    naero_node: int):
    ntime = case.data.shape[0]
    data_out = np.zeros((ntime, naero_node, 6), dtype='float32')
    data_out2 = data_out[:, :, :3]
    for i in range(ntime):
        datai = disp.data[i, :, :]
        datai2 = deform_aero_spline(
            structure_model,
            aero_model=None,
            nids=None,
            xyz_cid0=None,
            displacement0=datai)
        assert datai2.shape[0] == naero_node, (datai2.shape, naero_node)
        #naero_node, three = datai2.shape
        # assert three == 3, three
        data_out2[i, :, :] = datai2
    assert data_out.shape == (ntime, naero_node, 6), (data_out.shape, (ntime, naero_node, 6))

    # save the ata
    node = np.arange(1, naero_node+1, dtype='int32')
    gridtype = np.zeros(naero_node, dtype='int32')
    case.data = data_out
    case.node_gridtype = np.column_stack([node, gridtype])
    return case

def deform_aero_spline(structure_model: BDF,
                       aero_model: Optional[BDF]=None,
                       nids: Optional[np.ndarray]=None,
                       xyz_cid0: Optional[np.ndarray]=None,
                       displacement0: Optional[np.ndarray]=None):
    """
    Parameters
    ----------
    structure_model : BDF()
        the structural object for the structural mesh/deflections
    aero_model : BDF()
        the model object for the aero mesh
        expected model consists of CQUAD4s and GRIDs and that's it
    nids : (nnodes,) int array
        the node ids
    xyz_cid0 : (nnodes,3) float array
        the xyz locations in the basic frame (cid=0)
        corresponds to the structural nodes
    displacement0 : (nnodes, 6) float array
        the displacements in the basic frame (cid=0)
        corresponds to the structural nodes

    Returns
    -------

    TODO: is the aero mesh connected or unconnected?
    """
    nids, xyz_cid0 = _get_nids_xyz_cid0(structure_model, nids, xyz_cid0)
    assert np.array_equal(nids, np.unique(nids))

    aero_model, aero_xyz_cid0_point, naero_boxes_all = _get_discrete_aero_mesh(
        structure_model, aero_model)
    log = aero_model.log
    maxs = aero_xyz_cid0_point.max(axis=0)
    log.info(f'aero_xyz_cid0_point.max() = {maxs}')

    if displacement0 is None:
        # dummy displacment - parabola in 1, scaled to 1
        nnids = len(nids)
        assert xyz_cid0.shape == (nnids, 3), (nnids, xyz_cid0.shape)
        displacement0 = np.zeros((nnids, 6), dtype='float64')
        y = xyz_cid0[:, 2]
        ymax = np.abs(y).max()
        disp = (xyz_cid0[:, 2] / ymax) ** 2
        displacement0[:, 2] = disp
    assert len(displacement0.shape) == 2, displacement0.shape
    assert displacement0.shape[1] == 6, displacement0.shape

    all_aero_box = _get_all_aero_boxs(structure_model)
    assert aero_xyz_cid0_point.shape == (naero_boxes_all * 4, 3), aero_xyz_cid0_point.shape
    aero_xyz_cid0_point = _deform_aero_spline(
        structure_model, nids, xyz_cid0,
        displacement0,
        aero_model,
        aero_xyz_cid0_point,
        all_aero_box)
    return aero_xyz_cid0_point


def _deform_aero_spline(structure_model: BDF,
                        nids: np.ndarray,
                        xyz_cid0: np.ndarray,
                        displacement0: np.ndarray,
                        aero_model: BDF,
                        aero_xyz_cid0_point: np.ndarray,
                        all_aero_box: np.ndarray) -> np.ndarray:
    for eid, spline in sorted(structure_model.splines.items()):
        comment = f'{spline.type} eid={eid}'
        cid = max(aero_model.coords) + 1

        # get the indices for the aero boxes
        boxs = spline.aero_element_ids
        nboxi = len(boxs)
        ibox = np.searchsorted(all_aero_box, boxs)
        naero_nodesi = nboxi * 4

        # get the structural deflections for the spline points
        set_ids = np.unique(spline.setg_ref.ids)
        nstructure = len(set_ids)
        missing = np.setdiff1d(set_ids, nids)
        #log.info(f'set_ids = {set_ids}')
        #log.info(f'missing = {missing}')

        inid = np.searchsorted(nids, set_ids)
        assert len(missing) == 0, missing
        #log.info(f'inid = {inid}; type={str(type(inid))}')
        #log.info(f'nids = {nids}; type={str(type(nids))}')
        mapped_nids = nids[inid]
        #log.info(f'mapped_nids = {mapped_nids}')
        assert np.allclose(set_ids, mapped_nids), 'missing nodes'

        xyzi_cid0 = xyz_cid0[inid, :]
        disp0 = displacement0[inid, :3]
        assert disp0.shape == (nstructure, 3), disp0.shape

        # get the aero spline xyz locations
        # preallocate
        aero_disp_spline_cid0 = np.zeros((naero_nodesi, 3), dtype='float64')

        # split aero model into individual panels
        aero_xyz_cid0_panels = np.zeros((naero_nodesi, 3))
        j0 = 0
        j1 = 4
        # nboxi = len(ibox)
        # iibox = np.zeros((nboxi, 4), dtype='int32')
        # iibox[:, 0] = 4*ibox
        # iibox[:, 1] = 4*ibox + 1
        # iibox[:, 2] = 4*ibox + 2
        # iibox[:, 3] = 4*ibox + 3
        # iibox2 = iibox.ravel()
        # aero_xyz_cid0_panels2 = aero_xyz_cid0_point[iibox2, :]
        for iboxi in ibox:
            i0 = 4 * iboxi
            i1 = 4 * (iboxi + 1)
            aero_xyz_cid0_panels[j0:j1, :] = aero_xyz_cid0_point[i0:i1]
            j0 += 4
            j1 += 4

        caero = spline.caero_ref
        spline_type = spline.type
        caero_type = caero.type
        if spline_type == 'SPLINE1':
            if caero_type == 'CAERO1':
                # make the local coordinate system
                p1, p2, p3, p4 = caero.get_points()
                assert len(p1) == 3, p1
                xaxis = p2 - p1
                origin = (p1 + p2 + p3 + p4) / 4.
                assert len(origin) == 3, origin
                normal = np.cross(p3 - p1, p4 - p2)
                znorm = np.linalg.norm(normal)
                xzplane = origin + xaxis
                assert znorm > 0, znorm
                zaxis = normal / znorm
                coord = aero_model.add_cord2r(cid, origin, zaxis, xzplane, comment=comment)

                # local to global
                xform = coord.beta()

                # TODO: verify the transform
                # disp' = xform.T @ disp
                # transform from the basic frame (cid=0) to the local frame
                xyz_local = np.einsum('jk,ij->ik', xform.T, xyzi_cid0)
                txyz_disp_local = np.einsum('jk,ij->ik', xform.T, disp0)
                xyz_aero_local = np.einsum('jk,ij->ik', xform.T, aero_xyz_cid0_panels)

                # to local
                # xyz_local = xyz @ xform.T

                # maxs = aero_xyz_cid0_point.max(axis=0)
                # log.info(f'aero_xyz_cid0_panels.max() = {maxs}')
                # maxs = xyz_aero_local.max(axis=0)
                # log.info(f'xyz_aero_local.max() = {maxs}')

                # grab the xy location and the z deflection
                tz_disp_local = txyz_disp_local[:, 2]
                assert xyz_local.shape == (nstructure, 3)
                assert len(tz_disp_local) == nstructure

                # TODO: verify the transform (should be flipped vs. the first xform)
                tz_spline_local_deflected = _map_spline1_deflections(
                    xyz_local, tz_disp_local, xyz_aero_local)
                naero_node = len(xyz_aero_local)

                txyz_spline_local_deflected = np.zeros((naero_node, 3))
                txyz_spline_local_deflected[:, 2] = tz_spline_local_deflected
                txyz_spline_basic_deflected = np.einsum('jk,ij->ik', xform, txyz_spline_local_deflected)
                aero_disp_spline_cid0 += txyz_spline_basic_deflected

            elif caero_type == 'CAERO2':
                points = self.caero_ref.get_points()
                p1 = points[0, :]
                p2 = points[1, :]
                xaxis = p2 - p1
                xnorm = np.linalg.norm(xaxis)
                assert xnorm > 0, xnorm
                xaxis /= xnorm
                origin = p1
                zaxis = self.caero_ref.coord_ref.z
                xzplane = origin + xaxis
                coord = aero_model.add_cord2r(cid, origin, zaxis, xzplane, comment=comment)
                xform = coord.beta()

                # transform from the basic frame (cid=0) to the local frame
                xyz_local = np.einsum('jk,ij->ik', xform.T, xyzi_cid0)
                txyz_disp_local = np.einsum('jk,ij->ik', xform.T, disp0)
                ty_disp_local = txyz_disp_local[:, 1]
                tz_disp_local = txyz_disp_local[:, 2]
                raise NotImplementedError((spline_type, caero_type))
            else:
                raise NotImplementedError((spline_type, caero_type))
        else:
            raise NotImplementedError((spline_type, caero_type))

        # for each spline...
        for iboxi, xyz_aeroi in zip(ibox, aero_disp_spline_cid0):
            i0 = 4 * iboxi
            i1 = 4 * (iboxi + 1)
            aero_xyz_cid0_point[i0:i1, :] = xyz_aeroi
    assert len(aero_xyz_cid0_point.shape) == 2, aero_xyz_cid0_point.shape
    assert aero_xyz_cid0_point.shape[1] == 3, aero_xyz_cid0_point.shape
    return aero_xyz_cid0_point


def _map_spline1_deflections(xyz_structure: np.ndarray,
                             w_structure: np.ndarray,
                             xyz_aero: np.ndarray) -> np.ndarray:
    """
    {w_aero} = [xK] [C]^-1 {w_structure}

    Eq 1-56
    -------
    w = [C] {P}

    Eq 1-41
    -------
    {w_structure} = [C]{P}
    so:
    {P} = [C]^-1 {w_structure}

    Eq 1-42
    -------
               [1   x1A y1A K1A,1, K2A,2 ... K1nA,n]        {0 }
    {w_aero} = [1   x2A y2A K2A,1, K2A,2 ... K2nA,n] [C]^-1 {0 }
               [... ... ... ...    ...   ... ...   ]        {0 } = [xK] [C]^-1 {w_structure}
               [1   xnA ynA KnA,1, K2A,2 ... KnnA,n]        {w1}
                                                            {w2}
                                                            {wn}
    so:
    {w_aero} = [xK] [C]^-1 {w_structure} = [xK] {Cws} = [xK] {P}

    """
    # nids, xyz_local, tz_disp_local, xyz_aero_local
    delta_max = np.abs(w_structure).max()
    #print(f'delta_max = {delta_max}')

    if delta_max == 0.0:
        naero = len(xyz_aero)
        w_aero = np.zeros(naero, dtype='float64')
    else:
        #print(f'w_structure = {w_structure.tolist()}')
        C = get_c_matrix(xyz_structure)
        #print('nids =', nids)
        #print('deflections_structure =', deflections_structure.shape)
        #w_structure = deflections_structure[:, 2]

        w_structure2 = np.hstack([[0., 0., 0.], w_structure])
        w_aero = get_w_aero(C, w_structure2, xyz_structure, xyz_aero)
    return w_aero


def get_c_matrix(xyz_structure: np.ndarray) -> np.ndarray:
    """
    Parameters
    ----------
    xyz_structure

    Returns
    -------
    C

    """
    d = 1.
    pid16 = np.pi * d * 16.

    nnodes = len(xyz_structure)
    #print(f'nids = {nids}; n={nnodes}')
    #print(f'xyz_structure.shape = {str(xyz_structure.shape)}')
    i = 3
    c_array1 = np.zeros((3+nnodes, 3+nnodes), dtype='float64')
    c_array2 = np.zeros((3+nnodes, 3+nnodes), dtype='float64')
    for xyz_cidi in xyz_structure:
        xi, yi, zi = xyz_cidi

        c_array1[0, i] = 1.
        c_array1[1, i] = xi
        c_array1[2, i] = yi

        c_array1[i, 0] = 1.
        c_array1[i, 1] = xi
        c_array1[i, 2] = yi

        c_array2[0, i] = 1.
        c_array2[1, i] = xi
        c_array2[2, i] = yi

        c_array2[i, 0] = 1.
        c_array2[i, 1] = xi
        c_array2[i, 2] = yi

        xj = xyz_structure[:, 0]
        yj = xyz_structure[:, 1]
        rij2 = (xi-xj) ** 2. + (yi-yj) ** 2
        # we'll take the log of 0 if we don't do this...
        rij2[i-3] = 1.
        kij = rij2 * np.log(rij2) / pid16
        js = np.arange(3, 3+nnodes)
        assert len(kij) == len(js)
        """
        [[   0.       0.       0.       1.       1.       1.       1.   ]
         [   0.       0.       0.      74.396   26.338  -39.578  -69.958]
         [   0.       0.       0.      92.576   86.725   90.842   84.774]
         [   1.      74.396   92.576    0.     361.816 2448.467 4135.879]
         [   1.      26.338   86.725  361.816    0.     727.253 1685.963]
         [   1.     -39.578   90.842 2448.467  727.253    0.     131.111]
         [   1.     -69.958   84.774 4135.879 1685.963  131.111    0.   ]]
        """
        c_array1[i, js] = kij

        j = 3
        for xyz_cidj in xyz_structure:
            xj, yj, zj = xyz_cidj
            if i == j:
                c_array2[i, j] = 0.
            else:
                rij2 = (xi-xj)**2. + (yi-yj)**2  # Rij^2
                if rij2 == 0.:
                    c_array2[i, j] = 0.
                else:
                    # log=natural log
                    kij = rij2 * np.log(rij2) / pid16
                    c_array2[i, j] = kij
                    #msg = "i=%s j=%s xi=%s xj=%s yi=%s yj=%s Rij2=%s Kij=%s" %(i,j,xi,xj,yi,yj,rij2,kij)
                    #assert isinstance(kij, float64), msg
            j += 1
        i += 1
    if not np.allclose(c_array1, c_array2):  # pragma: no cover
        print(c_array1)
        print(c_array2)
        for i, r1, r2 in zip(count(), c_array1, c_array2):
            assert np.allclose(r1, r2), f'{i}\n{r1}\n{r2}'
        raise RuntimeError('bad [C]')
    det_c = np.linalg.det(c_array2)
    assert abs(det_c) > 0, det_c
    return c_array2


def get_w_aero(C: np.ndarray,
               w_structure: np.ndarray,
               xyz_structure: np.ndarray,
               xyz_aero: np.ndarray) -> np.ndarray:
    #print(f'C:\n{C}')
    #print(f'w_structure = {w_structure}')
    #print(f'C.shape = {str(C.shape)}, w.shape={str(w_structure.shape)}')
    cws = np.linalg.inv(C) @ w_structure  # Cws matrix, P matrix
    #P = solve(C, wS)
    #C*P=wS
    #P = C^-1*wS

    w_aero = get_xk(cws, xyz_structure, xyz_aero)
    #{w_aero} = [xK] [C] {w_structure}
    return w_aero


def get_xk(cws: np.ndarray,
           xyz_structure: np.ndarray,
           xyz_aero: np.ndarray) -> dict[int, np.ndarray]:
    """
    {w_aero} = [xK] [C] {w_structure}

    Parameters
    ----------
    cws : np.ndarray
        [C] matrix
    xyz_structure : np.ndarray
        the structural xyz locations (typically structural deflection nodes)
    xyz_aero : np.ndarray
        the aero xyz locations (typically CAERO1 nodes)

    Returns
    -------

    """
    d = 1.
    pid16 = np.pi * d * 16.

    naero = len(xyz_aero)
    nstr = len(xyz_structure)
    w_aero = np.zeros(naero, dtype='float64')

    for I, aero_xyzi in enumerate(xyz_aero):
        xk = np.zeros(nstr+3, dtype='float64')

        xa, ya, za = aero_xyzi
        xk[0] = 1.
        xk[1] = xa
        xk[2] = ya

        j = 3
        for str_xyzi in xyz_structure:
            (xs, ys, zs) = str_xyzi

            rij2 = (xa-xs)**2. + (ya-ys)**2  # rij^2
            if rij2 == 0.:
                xk[j] = 0.
            else:
                # log=natural log
                kij = rij2 * np.log(rij2) / pid16
                xk[j] = kij
            j += 1

        wai = xk @ cws
        w_aero[I] = wai  # [0, 0]
        #print "w[%s]=%s" % (I, wi[0,0])
    return w_aero
