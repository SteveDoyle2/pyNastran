"""
defines:
  - make_gpwg(Mgg, reference_point, xyz_cid0, log)

Na: a-set
Ns: s-set
No: omit-set
mass Na x Na mass matrix
stiffness Na x Na stiffness matrix
viscous Na x Na viscous damping matrix
hysteretic Na x Na hysteretic damping matrix
constraint No x Na back expansion matrix

mass MXX, MXXGEXT, MCGG, MRGG, MAAX, M100, M10
stiffness KXX, KXXGEXT, KCGG, KRGG, KAAX, K100, K10
constraint GO, GOGADGX, UEXP, UEXPT, MUG1, MUG1T

"""
from __future__ import annotations
from itertools import count
from typing import TYPE_CHECKING
import numpy as np
import scipy as sp

from pyNastran.bdf.mesh_utils.mass_properties import get_sub_eids
if TYPE_CHECKING:
    from pyNastran.bdf.bdf import BDF, CONM2, CMASS2
    from pyNastran.bdf.cards.coordinate_systems import Coord
    from cpylog import SimpleLogger


def make_mass_matrix(model: BDF,
                     fdtype: str='float64',
                     idtype: str='int32'):
    """
    Performs an accurate mass calculation

    ..todo:: not anywhere close to being done
    ..todo:: doesn't support SPOINTs/EPOINTs
    """
    log = model.log
    unused_icd_transform, icp_transform, xyz_cp, nid_cp_cd = model.get_displacement_index_xyz_cp_cd(
        fdtype=fdtype, idtype=idtype, sort_ids=True)

    nids = nid_cp_cd[:, 0]
    cps = nid_cp_cd[:, 1]
    cds = nid_cp_cd[:, 2]
    xyz_cid0 = model.transform_xyzcp_to_xyz_cid(
        xyz_cp, nids, icp_transform, cid=0, in_place=False, atol=1e-6)

    spoints = np.array(list(model.spoints.keys()), dtype=idtype)
    inids, components = _get_dofs_from_nid_cp_cd(spoints, nids)
    nrows = len(components)

    Mgg = sp.sparse.dok_matrix((nrows, nrows), dtype=np.float64)

    mass_cards = {
        'CROD', 'CONROD', 'CTUBE',
        'CBAR', 'CBEAM', 'CBEND', 'CBEAM3',
        'CMASS1', 'CMASS2', 'CMASS3', 'CMASS4',

        'CTRIA3', 'CTRIA6', 'CQUAD4', 'CQUAD8', 'CQUAD',

        'CTETRA', 'CPENTA', 'CHEXA', 'CPYRAM',
    }
    no_mass_cards = {
        'GRID', 'SPOINT', 'EPOINT',
        'CORD1R', 'CORD1C', 'CORD1S',
        'CORD2R', 'CORD2C', 'CORD2S',
        'EIGRL', 'EIGR', 'EIGB', 'EIGC',

        'FORCE', 'FORCE1', 'FORCE2',
        'MOMENT', 'MOMENT1', 'MOMENT2',
        'PLOAD', 'PLOAD1', 'PLOAD2', 'PLOAD4',

        'PELAS', 'PELAST', 'PDAMP', 'PDAMPT', 'PBUSH', 'PBUSHT',
        'PROD', 'PTUBE', 'PBAR', 'PBARL', 'PBEAM', 'PBEAML', 'PBCOMP',
        'PSHEAR', 'PSHELL', 'PCOMP', 'PCOMPG',
        'PSOLID', 'PLSOLID',

        'MAT1', 'MAT2', 'MAT4', 'MAT5', 'MAT8', 'MAT9', 'MAT10', 'MAT11',
        'MATT1', 'MATT2', 'MATT4', 'MATT5',
        'MATS1',

        'CELAS1', 'CELAS2', 'CELAS3', 'CELAS4', #'CLEAS5',
        'CDAMP1', 'CDAMP2', 'CDAMP3', 'CDAMP4', 'CDAMP5',
        'CBUSH', 'CBUSH1D', 'CBUSH2D', 'CVISC', 'CGAP', # is this right?
        'CFAST',
        'CRAC2D', 'CRAC3D',

        # aero
        'CAERO1', 'CAERO2', 'CAERO3', 'CAERO4', 'CAERO5',
        'PAERO1', 'PAERO2', 'PAERO3', 'PAERO4', 'PAERO5',
        'SPLINE1', 'SPLINE2', 'SPLINE3', 'SPLINE4', 'SPLINE5', 'SPLINE7',
        'AERO', 'AEROS', 'FLFACT', 'FLUTTER', 'TRIM', 'TRIM2', 'CSSCHD',

        'CBARAO',
        'CORD3G', 'CONV', 'CONVM', 'CSET', 'CSET1', 'CLOAD',
        'CHBDYG', 'CHBDYE', 'CHBDYP',
    }
    all_eids = np.hstack([
        list(model.elements),
        list(model.masses),
    ])
    ueids = np.unique(all_eids)
    assert len(all_eids) == len(ueids), 'duplicate element/mass ids'
    all_eids = ueids
    if len(all_eids) == 0:
        raise RuntimeError(f'No elements; model.elements={model.elements}')

    #etypes_skipped = set()
    for etype, eids in model._type_to_id_map.items():
        if etype in no_mass_cards:
            continue
        assert etype in mass_cards, etype

        #print(all_eids, eids, etype)
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
                Mgg[i1, i1] = mass_mat2[0, 0]
                Mgg[i2, i2] = mass_mat2[1, 1]
                Mgg[i3, i3] = mass_mat2[2, 2]

                Mgg[j1, j1] = mass_mat2[3, 3]
                Mgg[j2, j2] = mass_mat2[4, 4]
                Mgg[j3, j3] = mass_mat2[5, 5]

                #centroid = (xyz[n1] + xyz[n2]) / 2.
                #mass = _increment_inertia(centroid, reference_point, m, mass, cg, I)
        elif etype == 'CMASS2':
            eids2 = get_sub_eids(all_eids, eids, etype)
            for eid in eids2:
                elem = model.masses[eid] # type: CMASS2
                n1, n2 = elem.nodes
                c1, c2 = elem.c1, elem.c2
                assert n1 >= 0, n1
                assert n2 >= 0, n2
                i1, i2 = components[(n1, c1)], components[(n2, c2)]
                Mgg[i1, i2] = elem.mass

        elif etype == 'CONM2':
            mass_mat = np.zeros((6, 6), dtype=fdtype)
            eids2 = get_sub_eids(all_eids, eids, etype)
            for eid in eids2:
                elem = model.masses[eid] # type: CONM2
                n1 = elem.nid
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
                inertia = elem.Inertia()
                #inertia = np.ones((3, 3))
                mass_mat[3:, 3:] = inertia
                #print(mass_mat)

                i1, i2, i3 = components[(n1, 1)], components[(n1, 2)], components[(n1, 3)]
                i4, i5, i6 = components[(n1, 4)], components[(n1, 5)], components[(n1, 6)]
                Mgg[i1, i1] += mass_mat[0, 0]
                Mgg[i2, i2] += mass_mat[1, 1]
                Mgg[i3, i3] += mass_mat[2, 2]

                Mgg[i4, i4] += mass_mat[3, 3]
                Mgg[i4, i5] += mass_mat[3, 4]
                Mgg[i5, i4] += mass_mat[3, 4]
                Mgg[i4, i6] += mass_mat[3, 5]
                Mgg[i6, i4] += mass_mat[3, 5]

                Mgg[i5, i5] += mass_mat[4, 4]
                Mgg[i5, i6] += mass_mat[4, 5]
                Mgg[i6, i5] += mass_mat[4, 5]

                Mgg[i6, i6] += mass_mat[5, 5]

        else:
            log.warning(f'skipping etype={etype}')

    return Mgg, nids, cps, cds, xyz_cid0, spoints

def make_reduced_mass_matrix(model: BDF,
                             reference_point: np.ndarray,
                             fdtype: str='float64',
                             idtype: str='int32'):
    """
    Performs an accurate RB mass calculation

    ..todo:: not anywhere close to being done
    ..todo:: doesn't support SPOINTs/EPOINTs
    """
    Mgg, nids, cps, cds, xyz_cid0, spoints = make_mass_matrix(
        model, fdtype=fdtype, idtype=idtype)
    log = model.log
    D, Mo, S, mass, cg, II, IQ, Q = make_gpwg(
        Mgg, reference_point,
        nids, cps, cds, xyz_cid0,
        spoints, model.coords, log)
    return Mgg, D, Mo, S, mass, cg, II, IQ, Q

def make_gpwg(Mgg: sp.sparse.dok_matrix,
              reference_point: np.ndarray,
              grid_nids: np.ndarray,
              grid_cps: np.ndarray,
              grid_cds,
              xyz_cid0: np.ndarray,
              spoints: np.ndarray,
              coords: dict[int, Coord],
              log) -> tuple[str, ]:
    """
    Calculates the Grid Point Weight Generator (GPWG) table.

    The GPWG is performed on the g-size mass matrix, which
    is the mass matrix prior to the processing of the
    rigid elements, MPCs, and SPCs.

    The mass at scalar points and fluid-related masses are
    not included in the GPWG calculation.

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
    fdtype = 'float32'
    #grid_nids: np.ndarray,
                          #grid_cps: np.ndarray,
                          #xyz_cid0: np.ndarray,
                          #coords: dict[int, Coord],
                          #spoints: np.ndarray,
                          #reference_point: np.ndarray,
                          #log: SimpleLogger,

    D = get_6x6_rb_matrix(
        grid_nids, grid_cps, grid_cds, xyz_cid0,
        coords, spoints, reference_point, log,
        fdtype=fdtype)
    Mo, S, mass, cg, II, IQ, Q = _D_to_inertia(
        Mgg, D, log, fdtype=fdtype)
    return D, Mo, S, mass, cg, II, IQ, Q


def _D_to_inertia(Mgg: sp.sparse.dok_matrix,
                  D: np.ndarray,
                  log: SimpleLogger,
                  fdtype: str='float32') -> tuple[
                      np.ndarray,  # Mo
                      np.ndarray,  # S
                      np.ndarray,  # mass
                      np.ndarray,  # cg
                      np.ndarray,  # II
                      np.ndarray,  # IQ
                      np.ndarray,  # Q
                      ]:
    Mo = np.zeros((6, 6), dtype=fdtype)
    #print('D=\n%s\n' % D)
    # translati

    #assert D.shape[1] == Mgg.shape[0], f'D.shape={D.shape} Mgg.shape={Mgg.shape}'
    assert Mgg.shape[1] == D.shape[0], f'D.shape={D.shape} Mgg.shape={Mgg.shape}'
    try:
        Mo = triple(D, Mgg)
    except ValueError:
        raise ValueError(f'D.shape={D.shape} Mgg.shape={Mgg.shape}\n [Mo] = [D][Mgg][D]')
    log.info('Mgg=\n%s\n' % Mgg)
    log.info('Mo=\n%s\n' % Mo)

    # t-translation; r-rotation
    Mtt_bar = Mo[:3, :3]
    Mtr_bar = Mo[:3, 3:]
    #Mrt_bar = Mo[3:, :3]
    Mrr_bar = Mo[3:, 3:]

    #print('dinner =', diag(Mtt_bar))
    delta = np.linalg.norm(np.diag(Mtt_bar))
    #print('einner =', Mt_bar - diag(Mtt_bar))
    epsilon = np.linalg.norm([
        Mtt_bar[0, 1],
        Mtt_bar[0, 2],
        Mtt_bar[1, 2],
    ])


    log.info('Mt_bar (correct) =\n%s\n' % Mtt_bar)
    log.info(f'delta   = {delta:g}')
    log.info(f'epsilon = {epsilon:g}')
    log.info(f'e/d     = {epsilon / delta:g}\n')
    if epsilon / delta > 0.001:
        # user warning 3042
        unused_eigvals, S = np.linalg.eigh(Mtt_bar)
        msg = (
            '*** USER WARNING MESSAGE 3042 MODULE = GPWG\n'
            f'INCONSISTENT SCALAR MASSES HAVE BEEN USED. EPSILON/DELTA = {epsilon/delta:.7E}\n')
        log.warning(msg)

        Mtt = triple(S, Mtt_bar)  # translational
        Mtr = triple(S, Mtr_bar)
        Mrr = triple(S, Mrr_bar)  # rotational
    else:
        S = np.eye(3, dtype=Mtt_bar.dtype)

        Mtt = Mtt_bar  # translational
        Mtr = Mtr_bar
        Mrr = Mrr_bar  # rotational

    #log.info('omega=%s' % omega)
    log.info('S (right, but not correct order) =\n%s\n' % S)

    # 4. determine the principal axis & cg in the principal mass axis system
    # eq G-18
    Mtx = Mtt[0, 0]
    Mty = Mtt[1, 1]
    Mtz = Mtt[2, 2]
    mass = np.diag(Mtt)
    log.info('mass = %s' % mass)
    #if min(mass) == 0.:
        #raise RuntimeError('mass = %s' % mass)
    cg = np.array([
        [ Mtr[0, 0], -Mtr[0, 2],  Mtr[0, 1]],
        [ Mtr[1, 2],  Mtr[1, 1], -Mtr[1, 0]],
        [-Mtr[2, 1],  Mtr[2, 0],  Mtr[2, 2]],
    ], dtype=fdtype)
    if mass[0] != 0.:
        cg[0, :] /= Mtx
    if mass[1] != 0.:
        cg[1, :] /= Mty
    if mass[2] != 0.:
        cg[2, :] /= Mtz
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

    # 5. Following the center of gravity calculation is the calculation to
    # determine the moments of inertia for the center of gravity with
    # respect to the principal mass axes as shown in Equation J-19.
    I11 = Mrr[0, 0] - Mty * zy ** 2 - Mtz * yz ** 2
    I21 = I12 = -Mrr[0, 1] - Mtz * xz * yz
    I13 = I31 = -Mrr[0, 2] - Mty * xy * zy
    I22 = Mrr[1, 1] - Mtz * xz ** 2 - Mtx * zx ** 2
    I23 = -Mrr[1, 2] - Mtx * yx * zx
    I32 = I23
    I33 = Mrr[2, 2] - Mtx * yx ** 2 - Mty * xy ** 2

    II = np.array([
        [I11, I12, I13],
        [I21, I22, I13],
        [I31, I32, I33],
    ], dtype=fdtype)
    II = np.nan_to_num(II)

    log.info('I(S)=\n%s\n' % II)

    # 6. Reverse the sign of the off diagonal terms
    np.fill_diagonal(-II, np.diag(II))
    #print('I~=\n%s\n' % II)
    if np.nan in II:
        Q = np.zeros((3, 3), dtype=fdtype)
    else:
        omegaQ, Q = np.linalg.eig(II)
    #i = argsort(omegaQ)
    log.info('omegaQ = %s' % omegaQ)
    log.info('Q -> wrong =\n%s\n' % Q)
    IQ = triple(Q, II)
    #print('I(Q) -> wrong =\n%s\n' % IQ)
    return Mo, S, mass, cg, II, IQ, Q

def triple(A, B):
    """
    Performs the following matrix triple product:

        [C] = [A][B][A]

    .. todo:: not validated
    """
    #assert A.shape[1] == B.shape[0]
    assert B.shape[1] == A.shape[0]
    #einsum_check = np.einsum('ia,aj,ka->ijk', A, B, A)
    C = A.T @ B @ A
    #assert np.allclose(C, einsum_check)
    return C


def _get_dofs_from_nid_cp_cd(spoints: np.ndarray,
                             nids: np.ndarray,
                             exclude_spoints: bool=False) -> tuple[dict[int, int],
                                                                   dict[tuple[int, int]]]:
    """gets the GRID/SPOINT components"""
    inids = {}
    components = {}
    i = 0
    for j, nid in enumerate(nids):
        if nid in spoints:
            if exclude_spoints:
                continue
            components[(nid, 0)] = i
            i += 1
        else:
            components[(nid, 1)] = i
            components[(nid, 2)] = i + 1
            components[(nid, 3)] = i + 2
            components[(nid, 4)] = i + 3
            components[(nid, 5)] = i + 4
            components[(nid, 6)] = i + 5
            i += 6
        inids[nid] = j
    return inids, components

def get_6x6_rb_matrix(grid_nids: np.ndarray,
                      grid_cps: np.ndarray,
                      grid_cds: np.ndarray,
                      xyz_cid0: np.ndarray,
                      coords: dict[int, Coord],
                      spoints: np.ndarray,
                      reference_point: np.ndarray,
                      log: SimpleLogger,
                      fdtype: str='float32') -> np.ndarray:
    """reduces the (N,N) mass matrix (MGG)"""
    isin = np.isin(grid_nids, spoints)
    nnodes = xyz_cid0.shape[0]
    nspoints = len(spoints)
    nd_dof = nnodes*6 - nspoints*5
    D = np.zeros((nd_dof, 6), dtype=fdtype)

    # we subtract ref point so as to not change xyz_cid0
    for i, isini, nid, xyz in zip(count(), isin, grid_nids, xyz_cid0 - reference_point):
        if isini:
            print(f'skipping {nid}')
            continue
        r1, r2, r3 = xyz
        j = i * 6
        Tr = np.array([[0., r3, -r2],
                       [-r3, 0., r1],
                       [r2, -r1, 0.]], dtype=fdtype)
        #print('Tr[%i]=\n%s\n' % (i+1, Tr))

        cd = grid_cds[i]
        Ti = coords[cd].beta().T
        if not np.array_equal(Ti, np.eye(3)):
            log.info(f'cd={cd:d}; Ti[{i+1:d}]=\n{Ti}\n')
        TiT = Ti.T
        d = np.zeros((6, 6), dtype=fdtype)
        d[:3, :3] = TiT
        d[3:, 3:] = TiT
        d[:3, 3:] = TiT @ Tr
        D[j:j+6, :] = d
    return D


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
