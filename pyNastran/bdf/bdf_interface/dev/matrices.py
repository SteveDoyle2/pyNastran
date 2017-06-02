# pylint: disable=C0103
"""
defines:
  - make_gpwg(Mgg, reference_point, xyz_cid0, log)
"""
from six import iteritems


def make_gpwg(Mgg, reference_point, xyz_cid0, log):
    """
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
    coord_cps : (ncoords, ) int ndarray
        array of cp values corresponding to the coordinate systems
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
    D = zeros((nnodes*6, 6), dtype='float32')

    # we subtract ref point so as to not change xyz_cid0
    for i, node in enumerate(xyz_cid0 - reference_point):
        r1, r2, r3 = node
        j = i * 6
        Tr = array([[0., r3, -r2],
                    [-r3, 0., r1],
                    [r2, -r1, 0.]], dtype='float32')
        #print('Tr[%i]=\n%s\n' % (i+1, Tr))

        cp = grid_cps[i]
        Ti = coords[cp].beta()
        if not array_equal(Ti, eye(3)):
            log.info('Ti[%i]=\n%s\n' % (i+1, Ti))
        TiT = Ti.T
        d = zeros((6, 6), dtype='float32')
        d[:3, :3] = TiT
        d[3:, 3:] = TiT
        d[:3, 3:] = dot(TiT, Tr)
        D[j:j+6, :] = d

    Mo = zeros((6, 6), dtype='float32')
    #print('D=\n%s\n' % D)
    # translati

    Mo = triple(D, Mgg)
    log.info('Mgg=\n%s\n' % Mgg)
    log.info('Mo=\n%s\n' % Mo)

    # t-translation; r-rotation
    Mt_bar = Mo[:3, :3]
    Mtr_bar = Mo[:3, 3:]
    Mrt_bar = Mo[3:, :3]
    Mr_bar = Mo[3:, 3:]

    #print('dinner =', diag(Mt_bar))
    delta = norm(diag(Mt_bar))
    #print('einner =', Mt_bar - diag(Mt_bar))
    epsilon = norm([
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
    omega, S = eigh(Mt_bar)
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
    mass = diag(Mt)
    log.info('mass = %s' % mass)
    #if min(mass) == 0.:
        #raise RuntimeError('mass = %s' % mass)
    cg = array([
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
    xx = cg[0, 0]
    yx = cg[0, 1]
    zx = cg[0, 2]

    xy = cg[1, 0]
    yy = cg[1, 1]
    zy = cg[1, 2]

    xz = cg[2, 0]
    yz = cg[2, 1]
    zz = cg[2, 2]
    I11 = Mr[0, 0] - My * zy ** 2 - Mz * yz ** 2
    I21 = I12 = -Mr[0, 1] - Mz * xz * yz
    I13 = I31 = -Mr[0, 2] - My * xy * zy
    I22 = Mr[1, 1] - Mz * xz ** 2 - Mx * zx ** 2
    I32 = I23 = -Mr[1, 2] - Mx * yx * zx
    I33 = Mr[2, 2] - Mx * yx ** 2 - My * xy ** 2
    II = array([
        [I11, I12, I13],
        [I21, I22, I13],
        [I31, I32, I33],
        ], dtype='float32')
    II = nan_to_num(II)
    log.info('I(S)=\n%s\n' % II)


    # 6. Reverse the sign of the off diagonal terms
    fill_diagonal(-II, diag(II))
    #print('I~=\n%s\n' % II)
    if nan in II:
        Q = zeros((3, 3), dtype='float32')
    else:
        omegaQ, Q = eig(II)
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
        for nid, node in iteritems(model.nodes):
            xyz[nid] = node.get_position()
    for caero_id, caero in iteritems(model.caeros):
        c = caero.get_centroids()

    for spline_id, spline in iteritems(model.splines):
        s = spline.spline_nodes

    Ajj = None
    return Ajj
