from itertools import count
from typing import TYPE_CHECKING
import numpy as np

from pyNastran.op2.op2 import OP2

if TYPE_CHECKING:  # pragma: no cover
    from cpylog import SimpleLogger
    from pyNastran.bdf.bdf import BDF

def get_transform(T1, Te):
    """
    [T_el] = [T1]^T [T]
    [T_element_in_local] = [T_local_in_basic]^T [T_elemental_in_basic]
    [T_element_in_local] = [T_basic_in_local]   [T_elemental_in_basic]
    """
    T = T1.T @ Te # TODO: is this the right order?; I think so...
    return T

def _transform_strain(T, exxii, eyyii, ezzii, exyii, eyzii, exzii):
    strain = np.array([
        [exxii, exyii, exzii],
        [exyii, eyyii, eyzii],
        [exzii, eyzii, ezzii],
    ])

    # TODO: is it T.t @ inner @ T
    #       or    T   @ inner @ T.T
    #
    # TODO: is it Ri  @ inner @ R
    #       or    R   @ inner @ Ri
    #
    # TODO: is the T / R order right?
    strain4 = T.T @ (R @ strain @ Ri) @ T
    #print(T)

    # TODO: or is it another strain?
    strain_cid = strain4

    exxiit = strain_cid[0, 0]
    eyyiit = strain_cid[1, 1]
    ezziit = strain_cid[2, 2]
    exyiit = strain_cid[0, 1]
    eyziit = strain_cid[1, 2]
    exziit = strain_cid[0, 2]
    return exxiit, eyyiit, ezziit, exyiit, eyziit, exziit


def __transform_solids(model: OP2):  # pragma: no cover
    """http://web.mit.edu/course/3/3.11/www/modules/trans.pdf"""
    R = np.array([
        [1., 0., 0.],
        [0., 1., 0.],
        [0., 0., 2.],
    ])
    Ri = np.array([
        [1., 0., 0.],
        [0., 1., 0.],
        [0., 0., 0.5],
    ])
    thetad = 20.
    theta = np.radians(thetad)
    s = np.sin(theta)
    c = np.cos(theta)
    sc = s * c
    c2 = c ** 2
    s2 = s ** 2
    oxx = 1.
    oyy = 2.
    ozz = 3.
    txy = 1.
    txz = 0.
    tyz = 0.
    Ar = np.array([
        [c2, s2, 2. * sc],
        [s2, c2, -2. * sc],
        [-sc, sc, c2 - s2],
    ])
    """
    {ox'         {ox}
    {oy'  = [Ar] {oy}
    {txy'        {txy}
    """
    from pyNastran.bdf.bdf import BDF
    bdf_model = BDF()
    bdf_model.add_grid(1, [1., 0., 0.], cp=0, cd=0, ps='', seid=0, comment='')
    bdf_model.add_grid(2, [1., 1., 0.], cp=0, cd=0, ps='', seid=0, comment='')
    bdf_model.add_grid(3, [0., 1., 0.], cp=0, cd=0, ps='', seid=0, comment='')
    bdf_model.add_grid(4, [0., 0., 1.], cp=0, cd=0, ps='', seid=0, comment='')
    ctetra = bdf_model.add_ctetra(1, 1, [1, 2, 3, 4],)
    bdf_model.add_psolid(1, 1, cordm=0)
    E = 3.0E7
    G = None
    nu = 0.3
    bdf_model.add_mat1(1, E, G, nu)
    bdf_model.cross_reference()

    # this is ACTUALLY the element coordinate system
    centroid, xe, ye, ze = ctetra.material_coordinate_system()
    T = np.vstack([xe, ye, ze]) # Te

    #  we're going to transform the Te

    stress = np.array([
        [oxx, txy, txz],
        [txy, oyy, tyz],
        [txz, tyz, ozz],
    ])
    #stress2 = Ar @ stress
    #strain2 = (R @ A @ Ri) @ strain

    #  which is it?
    stress3 = T.T @ stress @ T
    stress3t = T @ stress @ T.T

    #  this is a test that these are the same...
    stress4 = R @ stress @ Ri
    print(stress)
    print(stress4)
    print('------------')

    strain3 = T.T @ strain @ T
    #strain3t = T.T @ strain @ T
    #strain4 = R @ strain3 @ Ri

    # is this strain3 or strain3t; is it (R @ strain3x @ Ri) or (Ri @ strain3x @ R)
    strain4 = R @ strain3 @ Ri
    #print(T)
    #print(T @ T.T)
    #print(stress2)
    print(stress3)
    print(stress3t)
    x = 1

    R = np.array([
        [1., 0., 0.],
        [0., 1., 0.],
        [0., 0., 2.],
    ])
    Ri = np.array([
        [1., 0., 0.],
        [0., 1., 0.],
        [0., 0., 0.5],
    ])
    nodes_desired = [213972, 213973, 213974, 213975, 213980, 213982, 213989, 213990,
                     213998, 213999, 214420, 214431, 214457, 214458, 214459, 214460]


def transform_solids(bdf_model: BDF, op2_model: OP2, cid: int):
    """
    http://web.mit.edu/course/3/3.11/www/modules/trans.pdf

    [stress_out] = [T_out] [stress_0] [T_out]^T
    [T_out]^T [stress_0] [T_out] = [stress_0] [T_out]
    [stress_out] = [T_out] [T_in]^T [stress_in] [T_in] [T_out]^T
    [stress_out] = [T] [stress_in] [T]^T
    [T] = [T_out] [T_in]^T
    """
    Tout = np.eye(3, dtype='float64')
    if cid != [-1, 0]:
        coord_out = bdf_model.coords[cid]
        assert coord_out.type in ['CORD2R', 'CORD1R'], coord_out
        Tout = coord_out.beta()

    #coord = model.coords[1]
    #T1 = coord.beta()


    # TODO: should be this...
    #strain_obj = model.op2_results.strain.ctetra_strain[1]

    # TODO: all we have for now
    #stress_obj = op2_model.ctetra_stress[1]
    result_types = ['ctetra_stress', 'cpenta_stress', 'cpyram_stress', 'chexa_stress']
    for res_type in result_types:
        res_dict = getattr(op2_model, res_type)
        for subcase, stress_obj in res_dict.items():
            _transform_solid_stress_obj(bdf_model, stress_obj, Tout)


def _transform_solid_stress_obj(bdf_model: BDF, stress_obj, Tout):
   #['oxx', 'oyy', 'ozz', 'txy', 'tyz', 'txz', 'omax', 'omid', 'omin', 'von_mises']
    data = stress_obj.data
    nmodes = data.shape[0]
    if stress_obj.is_stress:
        oxx = data[:, :, 0]
        oyy = data[:, :, 1]
        ozz = data[:, :, 2]
        txy = data[:, :, 3]
        tyz = data[:, :, 4]
        txz = data[:, :, 5]
    else:
        exx = data[:, :, 0]
        eyy = data[:, :, 1]
        ezz = data[:, :, 2]
        exy = data[:, :, 3] / 2
        eyz = data[:, :, 4] / 2
        exz = data[:, :, 5] / 2

    nnodes = 5 # CTETRA4 / CTETRA10
    eids = stress_obj.element_node[:, 0]
    neids = len(eids) // nnodes
    nodes = stress_obj.element_node[:, 1].reshape(neids, nnodes)
    ueids = np.unique(eids)

    eids_cids = stress_obj.element_cid
    cids = eids_cids[:, 1]
    ucids = np.unique(cids)

    for ucid in ucids:
        if ucid == 0:
            continue
        if ucid == cid:
            continue

        ieids = np.where(cids == ucid)[0]
        ueids = eids_cids[ieids, 0]
        nodes_xyz = {}
        for nid in np.unique(nodes[ieids, 1:].ravel()):
            node = bdf_model.nodes[nid] # type: GRID
            nodes_xyz[nid] = node.get_position_wrt(bdf_model, ucid)

        for eid, ieid in zip(ueids, ieids):
            i0 = ieid * nnodes
            i1 = i0 + nnodes
            nodes_eid = nodes[ieid, 1:]
            assert len(nodes_eid) == nnodes - 1
            e_nodes = np.vstack([nodes_xyz[nid] for nid in nodes_eid])
            avg_node = e_nodes.mean(axis=0)
            assert len(avg_node) == 3

            if ucid == -1:
                element = bdf_model.elements[eid]
                centroid, xe, ye, ze = element.material_coordinate_system()
                Te = np.vstack([xe, ye, ze]) # Te
            else:
                coord_in = bdf_model.coords[ucid]
                Tin = coord_in.beta()
                if coord_in.type in ['CORD2R', 'CORD1R']:
                    pass
                #elif coord_in.type in ['CORD2C', 'CORD1C']:
                    #thetad = avg_node[0]
                    #print(avg_node)
                    #theta = np.radians(thetad)
                    #s = np.sin(theta)
                    #c = np.cos(theta)
                    #sc = s * c
                    #c2 = c ** 2
                    #s2 = s ** 2
                    #Ar = np.array([
                        #[c2, s2, 2. * sc],
                        #[s2, c2, -2. * sc],
                        #[-sc, sc, c2 - s2],
                    #])
                    #Ar = np.array([
                        #[c, s, 0.],
                        #[-s, c, 0.],
                        #[0., 0., 1.],
                    #])
                    #Tin2 = Tin @ Ar.T
                else:
                    raise NotImplementedError(coord_in)
            T = Tout @ Tin.T

            for itime in range(nmodes):
                for ielem in range(i0, i1):
                    #exx = data[itime, ielem, 0]
                    #eyy = data[itime, ielem, 1]
                    #ezz = data[itime, ielem, 2]
                    #exy_2 = data[itime, ielem, 3] / 2
                    #eyz_2 = data[itime, ielem, 4] / 2
                    #exz_2 = data[itime, ielem, 5] / 2
                    #strain = np.array([
                        #[exx, exy_2, exz_2],
                        #[exy_2, eyy, eyz_2],
                        #[exz_2, eyz_2, ezz],
                    #])

                    oxx = data[itime, ielem, 0]
                    oyy = data[itime, ielem, 1]
                    ozz = data[itime, ielem, 2]
                    txy = data[itime, ielem, 3]
                    tyz = data[itime, ielem, 4]
                    txz = data[itime, ielem, 5]
                    stress = np.array([
                        [oxx, txy, txz],
                        [txy, oyy, tyz],
                        [txz, tyz, ozz],
                    ])
                    #[stress_out] = [T] [stress_in] [T]^T
                    T11 = T[0, 0]
                    T22 = T[1, 1]
                    #T33 = T[2, 2]
                    T12 = T[0, 1]
                    T13 = T[0, 2]
                    T32 = T23 = T[1, 2]

                    T21 = T[1, 0]
                    T31 = T[2, 0]
                    oxx2 = (oxx*T11**2 + oyy*T21**2 + ozz*T31**2 + 2*txy*T11*T12 + 2*txz*T11*T13 + 2*tyz*T21*T31)
                    oyy2 = (oxx*T12**2 + oyy*T22**2 + ozz*T32**2 + 2*txy*T11*T12 + 2*txz*T11*T13 + 2*tyz*T21*T31)
                    ozz2 = (oxx*T11**2 + oyy*T21**2 + ozz*T31**2 + 2*txy*T11*T12 + 2*txz*T11*T13 + 2*tyz*T21*T31)
                    #oxx2 = (oxx*T11**2 + oyy*T21**2 + ozz*T31**2 + 2*txy*T11*T12 + 2*txz*T11*T13 + 2*tyz*T21*T31)
                    stress2 = T @ stress @ T.T
                    print(eid)
                    print(stress)
                    print(stress2)
                    #ss

    inid0 = 0
    for ieid, eid in enumerate(ueids):

        # ------------------------------------------
        # this is ACTUALLY the element coordinate system
        ctetra = bdf_model.elements[eid]
        #print(ctetra)
        T = get_transform(T1, Te)

        for unused_irange in range(5):
            exxi = exx[:, inid0]
            eyyi = eyy[:, inid0]
            ezzi = ezz[:, inid0]
            exyi = exy[:, inid0]
            eyzi = eyz[:, inid0]
            exzi = exz[:, inid0]

            #  mode loop
            for imode, exxii, eyyii, ezzii, exyii, eyzii, exzii in zip(count(), exxi, eyyi, ezzi, exyi, eyzi, exzi):
                exxiit, eyyiit, ezziit, exyiit, eyziit, exziit = _transform_strain(
                    T, exxii, eyyii, ezzii, exyii, eyzii, exzii)

                #  save_op2 method
                data[imode, inid0, :6] = exxiit, eyyiit, ezziit, exyiit, eyziit, exziit
            inid0 += 1

    #op2_filename_out = os.path.join(dirname, f'xform.op2')
    #op2_model.write_op2(op2_filename_out, post=-1, endian=b'<', skips=None, nastran_format='nx')


def main():  # pragma: no cover
    op2_filename = r'C:\NASA\m4\formats\git\pyNastran\models\solid_bending\solid_bending_coord1.op2'
    from pyNastran.op2.op2 import read_op2
    #from pyNastran.op2.op2_geom import read_op2_geom
    model = read_op2(op2_filename, load_geometry=True)
    cid = 0
    transform_solids(model, model, cid)

if __name__ == '__main__':  # pragma: no cover
    main()
