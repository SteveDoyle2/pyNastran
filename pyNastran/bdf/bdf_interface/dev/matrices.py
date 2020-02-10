# pylint: disable=C0103
"""
defines:
  - make_gpwg(Mgg, reference_point, xyz_cid0, log)

"""
import numpy as np
import scipy as sp
from pyNastran.bdf.mesh_utils.mass_properties import get_sub_eids

def _lambda_1d(v1):
    """
    ::
      3d  [l,m,n,0,0,0]  2x6
          [0,0,0,l,m,n]
    """
    #xyz1 = model.Node(n1).get_position()
    #xyz2 = model.Node(n2).get_position()
    #v1 = xyz2 - xyz1
    n = np.linalg.norm(v1)
    if n == 0:
        raise ZeroDivisionError(v1)
    v1 = v1 / n
    (l, m, n) = v1
    Lambda = np.zeros((2, 6), 'd')
    Lambda[0, 0] = Lambda[1, 3] = l
    Lambda[0, 1] = Lambda[1, 4] = m
    Lambda[0, 2] = Lambda[1, 5] = n
    return Lambda

def triple(A, B):
    """
    Performs the following matrix triple product:

        [C] = [A][B][A]

    .. todo:: not validated
    """
    return np.einsum('ia,aj,ka->ijk', A, B, A)

def make_mass_matrix(model, reference_point, fdtype='float64', idtype='int32'):
    """
    Performs an accurate mass calculation

    ..todo:: not anywhere close to being done
    ..todo:: doesn't support SPOINTs/EPOINTs
    """
    unused_icd_transform, icp_transform, xyz_cp, nid_cp_cd = model.get_displacement_index_xyz_cp_cd(
        fdtype=fdtype, idtype=idtype, sort_ids=True)
    xyz_cid0 = model.transform_xyzcp_to_xyz_cid(
        xyz_cp, icp_transform, cid=0, in_place=False, atol=1e-6)

    nids = nid_cp_cd[:, 0]
    cps = nid_cp_cd[:, 1]

    components = {}
    i = 0
    inids = {}
    for j, nid in enumerate(nids):
        inids[nid] = j
        components[(nid, 1)] = i
        components[(nid, 2)] = i + 1
        components[(nid, 3)] = i + 2
        components[(nid, 4)] = i + 3
        components[(nid, 5)] = i + 4
        components[(nid, 6)] = i + 5
        i += 6
    nrows = len(components)

    mass = sp.sparse.dok_matrix((nrows, nrows), dtype=np.float64)

    no_mass = [
        'CELAS1', 'CELAS2', 'CELAS3', 'CELAS4', #'CLEAS5',
        'CDAMP1', 'CDAMP2', 'CDAMP3', 'CDAMP4', 'CDAMP5',
        'CBUSH', 'CBUSH1D', 'CBUSH2D', 'CVISC', 'CGAP', # is this right?
        'CFAST',
        'CRAC2D', 'CRAC3D',

        'CSSCHD', 'CAERO1', 'CAERO2', 'CAERO3', 'CAERO4', 'CAERO5',
        'CBARAO', 'CORD1R', 'CORD2R', 'CORD1C', 'CORD2C', 'CORD1S', 'CORD2S',
        'CORD3G', 'CONV', 'CONVM', 'CSET', 'CSET1', 'CLOAD',
        'CHBDYG', 'CHBDYE', 'CHBDYP',
    ]
    all_eids = np.array(list(model.elements.keys()), dtype='int32')

    #etypes_skipped = set()
    for etype, eids in model._type_to_id_map.items():
        if etype in no_mass:
            continue

        if etype in ['CROD', 'CONROD']:
            eids2 = get_sub_eids(all_eids, eids, etype)

            # lumped
            mass_mat = np.ones((2, 2), dtype=fdtype)
            mass_mat[0, 0] = mass_mat[2, 2] = 1.
            #mass_mat[2, 2] = mass_mat[5, 5] = 0.

            #
            #mi = (rho * A * L + nsm) / 6.
            #m = array([[2., 1.],
                       #[1., 2.]])  # 1D rod consistent
            #m = array([[1., 0.],
                       #[0., 1.]])  # 1D rod lumped

            for eid in eids2:
                elem = model.elements[eid]
                n1, n2 = elem.node_ids

                i1, i2, i3 = components[(n1, 1)], components[(n1, 2)], components[(n1, 3)]
                j1, j2, j3 = components[(n2, 1)], components[(n2, 2)], components[(n2, 3)]

                inid1 = inids[n1]
                inid2 = inids[n2]
                v1 = xyz_cid0[inid2, :] - xyz_cid0[inid1, :]
                length = np.linalg.norm(v1)
                mpl = elem.MassPerLength()
                massi = mpl * length / 2.
                Lambda = _lambda_1d(v1)
                mass_mat2 = (Lambda.T @ mass_mat @ Lambda) * massi
                assert mass_mat2.shape == (6, 6), mass_mat2
                mass[i1, i1] = mass_mat2[0, 0]
                mass[i2, i2] = mass_mat2[1, 1]
                mass[i3, i3] = mass_mat2[2, 2]

                mass[j1, j1] = mass_mat2[3, 3]
                mass[j2, j2] = mass_mat2[4, 4]
                mass[j3, j3] = mass_mat2[5, 5]

                #centroid = (xyz[n1] + xyz[n2]) / 2.
                #mass = _increment_inertia(centroid, reference_point, m, mass, cg, I)
        elif etype == 'CONM2':
            mass_mat = np.zeros((6, 6), dtype=fdtype)
            eids2 = get_sub_eids(all_eids, eids, etype)
            for eid in eids2:
                elem = model.masses[eid]
                massi = elem.Mass()
                unused_rx, unused_ry, unused_rz = elem.X
                mass_mat[0, 0] = massi
                mass_mat[1, 1] = massi
                mass_mat[2, 2] = massi
                #mass_mat[3, 3] = i11
                #mass_mat[3, 4] = mass_mat[4, 3] = -i12
                #mass_mat[3, 5] = mass_mat[5, 3] = -i13
                #mass_mat[4, 4] = i22
                #mass_mat[4, 5] = mass_mat[5, 4] = -i23
                #mass_mat[5, 5] = i33
                mass_mat[3:5, 3:5] = elem.Inertia()

                i1, i2, i3 = components[(n1, 1)], components[(n1, 2)], components[(n1, 3)]
                i4, i5, i6 = components[(n1, 4)], components[(n1, 5)], components[(n1, 6)]
                j1, j2, j3 = components[(n2, 1)], components[(n2, 2)], components[(n2, 3)]
                j4, j5, j6 = components[(n2, 4)], components[(n2, 5)], components[(n2, 6)]
                mass[i1, j1] = mass_mat[0, 0]
                mass[i2, j2] = mass_mat[1, 1]
                mass[i3, j3] = mass_mat[2, 2]

                mass[i4, j4] = mass_mat[3, 3]
                mass[i4, j5] = mass[i5, j4] = mass_mat[3, 4]
                mass[i4, j6] = mass[i6, j4] = mass_mat[3, 5]

                mass[i5, j5] = mass_mat[4, 4]
                mass[i5, j6] = mass[i6, j5] = mass_mat[4, 5]
                mass[i6, j6] = mass_mat[5, 5]

        else:
            pass

    Mgg = mass
    return make_gpwg(Mgg, reference_point, xyz_cid0, cps, model.coords, model.log)

def make_gpwg(Mgg, reference_point, xyz_cid0, grid_cps, coords, log):
    """
    Calculates the Grid Point Weight Generator (GPWG) table.

    Parameters
    ----------
    reference_point : (3, ) float ndarray
        the reference point
    grid_point : int
        0->origin, x>0, that grid point
    Mgg : (N, N) matrix
        the mass matrix
    xyz_cid0 : (ngrids, 3) float ndarray
        the xyz coordinates of the grids
    grid_cps : (ngrids, ) int ndarray
        array of cp values corresponding to xyz_cid0
    coords : dict[cp] : Coord()
        dict of cp values corresponding to the Cp coordinate systems
    log : logger()
        logging object

    Returns
    -------
    Mo : (6, 6) float ndarray
        the rigid body mass matrix in the basic coordinate system
    S : (3, 3) float ndarray
        the scalar partition matrix (also known as the principal mass axes)
    mass : (3, ) float ndarray
        the mass in the 3 pricincipal (basic) axes
    cg : (3, 3) float ndarray
        the cg in the 3 principal (basic) axes
    II : (3, 3) float ndarray
        inertias relative to the CG
        also called I(S)
    IQ : (3, ) float ndarray
        principal moments of inertia about the CG
        also called I(Q)
    Q : (3, 3) float ndarray
        the coordinate transformation between the S axes and the Q axes

    .. todo:: doesn't consider SPOINTs/EPOINTs
    .. todo:: hasn't been tested
    """
    nnodes = xyz_cid0.shape[0]
    D = np.zeros((nnodes*6, 6), dtype='float32')

    # we subtract ref point so as to not change xyz_cid0
    for i, node in enumerate(xyz_cid0 - reference_point):
        r1, r2, r3 = node
        j = i * 6
        Tr = np.array([[0., r3, -r2],
                       [-r3, 0., r1],
                       [r2, -r1, 0.]], dtype='float32')
        #print('Tr[%i]=\n%s\n' % (i+1, Tr))

        cp = grid_cps[i]
        Ti = coords[cp].beta()
        if not np.array_equal(Ti, np.eye(3)):
            log.info('Ti[%i]=\n%s\n' % (i+1, Ti))
        TiT = Ti.T
        d = np.zeros((6, 6), dtype='float32')
        d[:3, :3] = TiT
        d[3:, 3:] = TiT
        d[:3, 3:] = TiT @ Tr
        D[j:j+6, :] = d

    Mo = np.zeros((6, 6), dtype='float32')
    #print('D=\n%s\n' % D)
    # translati

    Mo = triple(D, Mgg)
    log.info('Mgg=\n%s\n' % Mgg)
    log.info('Mo=\n%s\n' % Mo)

    # t-translation; r-rotation
    Mt_bar = Mo[:3, :3]
    Mtr_bar = Mo[:3, 3:]
    #Mrt_bar = Mo[3:, :3]
    Mr_bar = Mo[3:, 3:]

    #print('dinner =', diag(Mt_bar))
    delta = np.linalg.norm(np.diag(Mt_bar))
    #print('einner =', Mt_bar - diag(Mt_bar))
    epsilon = np.linalg.norm([
        Mt_bar[0, 1],
        Mt_bar[0, 2],
        Mt_bar[1, 2],
    ])
    if epsilon / delta > 0.001:
        # user warning 3042
        pass

    log.info('Mt_bar (correct) =\n%s\n' % Mt_bar)
    log.info('delta=%s' % delta)
    log.info('epsilon=%s' % epsilon)
    log.info('e/d=%s\n' % (epsilon / delta))

    # hermitian eigenvectors
    omega, S = np.linalg.eigh(Mt_bar)
    log.info('omega=%s' % omega)
    log.info('S (right, but not correct order) =\n%s\n' % S)

    Mt = triple(S, Mt_bar)
    Mtr = triple(S, Mtr_bar)
    Mr = triple(S, Mr_bar)

    # 4. determine the principal axis & cg in the principal mass axis system
    # eq G-18
    Mx = Mt[0, 0]
    My = Mt[1, 1]
    Mz = Mt[2, 2]
    mass = np.diag(Mt)
    log.info('mass = %s' % mass)
    #if min(mass) == 0.:
        #raise RuntimeError('mass = %s' % mass)
    cg = np.array([
        [Mtr[0, 0], -Mtr[0, 2], Mtr[0, 1]],
        [Mtr[1, 2], Mtr[1, 1], -Mtr[1, 0]],
        [-Mtr[2, 1], Mtr[2, 0], Mtr[2, 2]],
    ], dtype='float32')
    if mass[0] != 0.:
        cg[0, :] /= Mx
    if mass[1] != 0.:
        cg[1, :] /= My
    if mass[2] != 0.:
        cg[2, :] /= Mz
    #cg = nan_to_num(cg)

    log.info('cg=\n%s\n' % cg)
    #xx = cg[0, 0]
    yx = cg[0, 1]
    zx = cg[0, 2]

    xy = cg[1, 0]
    #yy = cg[1, 1]
    zy = cg[1, 2]

    xz = cg[2, 0]
    yz = cg[2, 1]
    #zz = cg[2, 2]
    I11 = Mr[0, 0] - My * zy ** 2 - Mz * yz ** 2
    I21 = I12 = -Mr[0, 1] - Mz * xz * yz
    I13 = I31 = -Mr[0, 2] - My * xy * zy
    I22 = Mr[1, 1] - Mz * xz ** 2 - Mx * zx ** 2
    I23 = -Mr[1, 2] - Mx * yx * zx
    I32 = I23
    I33 = Mr[2, 2] - Mx * yx ** 2 - My * xy ** 2
    II = np.array([
        [I11, I12, I13],
        [I21, I22, I13],
        [I31, I32, I33],
    ], dtype='float32')
    II = np.nan_to_num(II)

    log.info('I(S)=\n%s\n' % II)


    # 6. Reverse the sign of the off diagonal terms
    np.fill_diagonal(-II, np.diag(II))
    #print('I~=\n%s\n' % II)
    if np.nan in II:
        Q = np.zeros((3, 3), dtype='float32')
    else:
        omegaQ, Q = np.linalg.eig(II)
    #i = argsort(omegaQ)
    log.info('omegaQ = %s' % omegaQ)
    log.info('Q -> wrong =\n%s\n' % Q)
    IQ = triple(Q, II)
    #print('I(Q) -> wrong =\n%s\n' % IQ)

    return Mo, S, mass, cg, II, IQ, Q

def get_Ajj(model, xyz=None):
    """not finished"""
    if xyz is None:
        xyz = {}
        for nid, node in model.nodes.items():
            xyz[nid] = node.get_position()
    for unused_caero_id, caero in model.caeros.items():
        unused_centroids = caero.get_centroids()

    for unused_spline_id, spline in model.splines.items():
        unused_spline_nodes = spline.spline_nodes

    Ajj = None
    return Ajj
