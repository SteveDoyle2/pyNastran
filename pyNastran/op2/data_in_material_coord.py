"""
Defines:
 - data_in_material_coord(bdf, op2, in_place=False)

"""
from __future__ import annotations
import copy
from itertools import count
from typing import TYPE_CHECKING

import numpy as np
from numpy import cos, sin, cross
from numpy.linalg import norm  # type: ignore

from pyNastran.utils.numpy_utils import integer_types

if TYPE_CHECKING:  # pragma: no cover
    from pyNastran.bdf.bdf import BDF, GRID
    from pyNastran.op2.op2 import OP2

force_vectors = ['cquad4_force', 'cquad8_force', 'cquadr_force',
                 'ctria3_force', 'ctria6_force', 'ctriar_force']
stress_vectors = ['cquad4_stress', 'cquad8_stress', 'cquadr_stress',
                  'ctria3_stress', 'ctria6_stress', 'ctriar_stress']
strain_vectors = ['cquad4_strain', 'cquad8_strain', 'cquadr_strain',
                  'ctria3_strain', 'ctria6_strain', 'ctriar_strain']


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

def transf_Mohr(Sxx, Syy, Sxy, thetarad):
    """Mohr's Circle-based Plane Stress Transformation

    Parameters
    ----------
    Sxx, Syy, Sxy : array-like
        Sigma_xx, Sigma_yy, Sigma_xy stresses.
    thetarad : array-like
        Array with angles for which the stresses should be transformed.

    Returns
    -------
    Sxx_theta, Syy_theta, Sxy_theta : np.ndarray
        Transformed stresses.

    """
    Sxx = np.asarray(Sxx)
    Syy = np.asarray(Syy)
    Sxy = np.asarray(Sxy)
    thetarad = np.asarray(thetarad)
    Scenter = (Sxx + Syy)/2.
    R = np.sqrt((Sxx - Scenter)**2 + Sxy**2)
    thetarad_Mohr = np.arctan2(-Sxy, Sxx - Scenter) + 2*thetarad
    cos_Mohr = cos(thetarad_Mohr)
    Sxx_theta = Scenter + R*cos_Mohr
    Syy_theta = Scenter - R*cos_Mohr
    Sxy_theta = -R*sin(thetarad_Mohr)
    return Sxx_theta, Syy_theta, Sxy_theta


def thetadeg_to_principal(Sxx, Syy, Sxy):
    """Calculate the angle to the principal plane stress state

    Parameters
    ----------
    Sxx, Syy, Sxy : array-like
        Sigma_xx, Sigma_yy, Sigma_xy stresses.

    Returns
    -------
    thetadeg : np.ndarray
        Array with angles for which the given stresses are transformed to the
        principal stress state.

    """
    Scenter = (Sxx + Syy)/2.
    thetarad = np.arctan2(Sxy, Scenter - Syy)
    return np.rad2deg(thetarad)/2.


def get_eids_from_op2_vector(vector):
    """Obtain the element ids for a given op2 vector

    Parameters
    ----------
    vector : op2 vector
        An op2 vector obtained, for example, doing::

            vector = op2.cquad4_force[1]
            vector = op2.cquad8_stress[1]
            vector = op2.ctriar_force[1]
            vector = op2.ctria3_stress[1]

    """
    eids = getattr(vector, 'element', None)
    if eids is None:
        eids = vector.element_node[:, 0]
    return eids


def is_mcid(elem):
    """
    Determines if the element uses theta or the mcid (projected material coordinate system)

    Parameters
    ----------
    elem : varies
        an element object
        CQUAD4, CQUAD8, CQUADR
        CTRIA3, CTRIA6, CTRIAR

    Returns
    -------
    is_mcid : bool
        the projected material coordinate system is used

    """
    theta_mcid = getattr(elem, 'theta_mcid', None)
    return isinstance(theta_mcid, integer_types)


def check_theta(elem) -> float:
    theta = getattr(elem, 'theta_mcid', None)
    if theta is None:
        theta = 0.
    elif isinstance(theta, float):
        pass
    elif isinstance(theta, integer_types):
        raise ValueError('MCID is accepted by this function')
    return theta

def angle2vec(v1, v2):
    """
    Using the definition of the dot product to get the angle

    v1 o v2 = |v1| * |v2| * cos(theta)
    theta = np.arccos( (v1 o v2) / (|v1|*|v2|))

    """
    denom = norm(v1, axis=1) * norm(v2, axis=1)
    return np.arccos((v1 * v2).sum(axis=1) / denom)


def calc_imat(normals, csysi):
    """
    Calculates the i vector in the material coordinate system.

    j = k x ihat
    jhat = j / |j|
    i = jhat x k

    Notes
    -----
    i is not a unit vector because k (the element normal)
    is not a unit vector.
    """
    jmat = cross(normals, csysi) # k x i
    jmat /= norm(jmat)
    imat = cross(jmat, normals)
    return imat


def data_in_material_coord(bdf: BDF, op2: OP2, in_place: bool=False) -> OP2:
    """Convert OP2 2D element outputs to material coordinates

    Nastran allows the use of 'PARAM,OMID,YES' to print 2D element forces,
    stresses and strains based on the material direction. However, the
    conversion only takes place in the F06 output file, whereas the OP2 output
    file remains in the element coordinate system.

    This function converts the 2D element vectors to the material OP2
    similarly to most of the post-processing tools (Patran, Femap, HyperView,
    etc). It handles both 2D elements with MCID or THETA.

    Parameters
    ----------
    bdf : :class:`.BDF` object
        A :class:`.BDF` object that corresponds to the 'op2'.
    op2 : :class:`.OP2` object
        A :class:`.OP2` object that corresponds to the 'bdf'.
    in_place : bool; default=False
        If true the original op2 object is modified, otherwise a new one
        is created.

    Returns
    -------
    op2_new : :class:`.OP2` object
        A :class:`.OP2` object with the abovementioned changes.

    .. warning ::  doesn't handle composite stresses/strains/forces
    .. warning ::  doesn't handle solid stresses/strains/forces (e.g. MAT11)
    .. warning ::  zeros out data for CQUAD8s

    """
    if in_place:
        op2_new = op2
    else:
        op2_new = copy.deepcopy(op2)

    eids = np.array(list(bdf.elements.keys()))
    elems = np.array(list(bdf.elements.values()))
    mcid = np.array([is_mcid(e) for e in elems])
    elemsmcid = elems[mcid]
    elemstheta = elems[~mcid]

    thetadeg = np.zeros(elems.shape)
    thetadeg[~mcid] = np.array([check_theta(e) for e in elemstheta])
    thetarad = np.deg2rad(thetadeg)

    #NOTE separating quad types to get vectorizable "corner"
    quad_types = ['CQUAD4', 'CQUAD8', 'CQUADR']
    for quad_type in quad_types:
        # elems with THETA
        thisquad = np.array([quad_type in e.type for e in elemstheta])
        if np.any(thisquad):
            quadelems = elemstheta[thisquad]
            corner = np.array([e.get_node_positions() for e in quadelems])
            g1 = corner[:, 0, :]
            g2 = corner[:, 1, :]
            g3 = corner[:, 2, :]
            g4 = corner[:, 3, :]
            betarad = angle2vec(g3 - g1, g2 - g1)
            gammarad = angle2vec(g4 - g2, g1 - g2)
            alpharad = (betarad + gammarad) / 2.
            tmp = thetarad[~mcid]
            tmp[thisquad] -= betarad
            tmp[thisquad] += alpharad
            thetarad[~mcid] = tmp

        # elems with MCID
        thisquad = np.array([quad_type in e.type for e in elemsmcid])
        if np.any(thisquad):
            quadelems = elemsmcid[thisquad]
            corner = np.array([e.get_node_positions() for e in quadelems])
            g1 = corner[:, 0, :]
            g2 = corner[:, 1, :]
            g3 = corner[:, 2, :]
            g4 = corner[:, 3, :]
            normals = np.array([e.Normal() for e in quadelems])
            csysi = np.array([bdf.coords[e.theta_mcid].i for e in quadelems])
            imat = calc_imat(normals, csysi)
            tmp = thetarad[mcid]
            tmp[thisquad] = angle2vec(g2 - g1, imat)
            # getting sign of THETA
            check_normal = cross(g2 - g1, imat)
            tmp[thisquad] *= np.sign((check_normal * normals).sum(axis=1))
            betarad = angle2vec(g3 - g1, g2 - g1)
            gammarad = angle2vec(g4 - g2, g1 - g2)
            alpharad = (betarad + gammarad) / 2.
            tmp[thisquad] -= betarad
            tmp[thisquad] += alpharad
            thetarad[mcid] = tmp

    tria_types = ['CTRIA3', 'CTRIA6', 'CTRIAR']
    for tria_type in tria_types:
        # elems with MCID
        thistria = np.array([tria_type in e.type for e in elemsmcid])
        if np.any(thistria):
            triaelems = elemsmcid[thistria]
            corner = np.array([e.get_node_positions() for e in triaelems])
            g1 = corner[:, 0, :]
            g2 = corner[:, 1, :]
            g3 = corner[:, 2, :]
            normals = np.array([e.Normal() for e in triaelems])
            csysi = np.array([bdf.coords[e.theta_mcid].i for e in triaelems])
            imat = calc_imat(normals, csysi)
            tmp = thetarad[mcid]
            tmp[thistria] = angle2vec(g2 - g1, imat)
            # getting sign of THETA
            check_normal = cross(g2 - g1, imat)
            tmp[thistria] *= np.sign((check_normal * normals).sum(axis=1))
            thetarad[mcid] = tmp

    thetarad = dict([[eid, theta] for eid, theta in zip(eids, thetarad)])

    for vecname in force_vectors:
        op2_vectors = getattr(op2, vecname)
        new_vectors = getattr(op2_new, vecname)
        for subcase, vector in op2_vectors.items():
            new_vector = new_vectors[subcase]
            veceids = get_eids_from_op2_vector(vector)
            #NOTE assuming thetarad=0 for elements that exist in the op2 but
            #     not in the supplied bdf file
            vecthetarad = np.array([thetarad.get(eid, 0.0) for eid in veceids])

            if veceids.shape[0] == vector.data.shape[1] // 5:
                steps = [5, 5, 5, 5, 5]
            else: # assuming  always that veceids.shape[0] == vector.data.shape[1]
                steps = [1]

            for start, step in enumerate(steps):
                slicei = slice(start, vector.data.shape[1], step)
                # membrane terms
                Sxx = vector.data[:, slicei, 0]
                Syy = vector.data[:, slicei, 1]
                Sxy = vector.data[:, slicei, 2]
                if vector.data.dtype == np.complex64 or vector.data.dtype == np.complex128:
                    Sxx_theta_real, Syy_theta_real, Sxy_theta_real = transf_Mohr(
                        Sxx.real, Syy.real, Sxy.real, vecthetarad)
                    new_vector.data[:, slicei, 0].real = Sxx_theta_real
                    new_vector.data[:, slicei, 1].real = Syy_theta_real
                    new_vector.data[:, slicei, 2].real = Sxy_theta_real
                    Sxx_theta_imag, Syy_theta_imag, Sxy_theta_imag = transf_Mohr(
                        Sxx.imag, Syy.imag, Sxy.imag, vecthetarad)
                    new_vector.data[:, slicei, 0].imag = Sxx_theta_imag
                    new_vector.data[:, slicei, 1].imag = Syy_theta_imag
                    new_vector.data[:, slicei, 2].imag = Sxy_theta_imag
                else:
                    Sxx_theta, Syy_theta, Sxy_theta = transf_Mohr(Sxx, Syy, Sxy, vecthetarad)
                    new_vector.data[:, slicei, 0] = Sxx_theta
                    new_vector.data[:, slicei, 1] = Syy_theta
                    new_vector.data[:, slicei, 2] = Sxy_theta

                # bending terms
                Sxx = vector.data[:, slicei, 3]
                Syy = vector.data[:, slicei, 4]
                Sxy = vector.data[:, slicei, 5]
                if vector.data.dtype == np.complex64 or vector.data.dtype == np.complex128:
                    Sxx_theta_real, Syy_theta_real, Sxy_theta_real = transf_Mohr(
                        Sxx.real, Syy.real, Sxy.real, vecthetarad)
                    new_vector.data[:, slicei, 3].real = Sxx_theta_real
                    new_vector.data[:, slicei, 4].real = Syy_theta_real
                    new_vector.data[:, slicei, 5].real = Sxy_theta_real
                    Sxx_theta_imag, Syy_theta_imag, Sxy_theta_imag = transf_Mohr(
                        Sxx.imag, Syy.imag, Sxy.imag, vecthetarad)
                    new_vector.data[:, slicei, 3].imag = Sxx_theta_imag
                    new_vector.data[:, slicei, 4].imag = Syy_theta_imag
                    new_vector.data[:, slicei, 5].imag = Sxy_theta_imag

                else:
                    Sxx_theta, Syy_theta, Sxy_theta = transf_Mohr(Sxx, Syy, Sxy, vecthetarad)
                    new_vector.data[:, slicei, 3] = Sxx_theta
                    new_vector.data[:, slicei, 4] = Syy_theta
                    new_vector.data[:, slicei, 5] = Sxy_theta

                # transverse terms
                Qx = vector.data[:, slicei, 6]
                Qy = vector.data[:, slicei, 7]
                if vector.data.dtype == np.complex64 or vector.data.dtype == np.complex128:
                    Qx_new_real = cos(vecthetarad)*Qx.real + sin(vecthetarad)*Qy.real
                    Qy_new_real = -sin(vecthetarad)*Qx.real + cos(vecthetarad)*Qy.real
                    new_vector.data[:, slicei, 6].real = Qx_new_real
                    new_vector.data[:, slicei, 7].real = Qy_new_real
                    Qx_new_imag = cos(vecthetarad)*Qx.imag + sin(vecthetarad)*Qy.imag
                    Qy_new_imag = -sin(vecthetarad)*Qx.imag + cos(vecthetarad)*Qy.imag
                    new_vector.data[:, slicei, 6].imag = Qx_new_imag
                    new_vector.data[:, slicei, 7].imag = Qy_new_imag
                else:
                    Qx_new = cos(vecthetarad)*Qx + sin(vecthetarad)*Qy
                    Qy_new = -sin(vecthetarad)*Qx + cos(vecthetarad)*Qy
                    new_vector.data[:, slicei, 6] = Qx_new
                    new_vector.data[:, slicei, 7] = Qy_new

            #TODO implement transformation for corner nodes
            #     for now we just zero the wrong values
            if 'quad8' in vecname:
                for j in [1, 2, 3, 4]:
                    new_vector.data[:, j, :] = 0

    for vecname in stress_vectors:
        op2_vectors = getattr(op2, vecname)
        new_vectors = getattr(op2_new, vecname)
        for subcase, vector in op2_vectors.items():
            new_vector = new_vectors[subcase]
            veceids = get_eids_from_op2_vector(vector)
            check = veceids != 0
            veceids = veceids[check]
            #NOTE assuming thetarad=0 for elements that exist in the op2 but
            #     not in the supplied bdf file
            vecthetarad = np.array([thetarad.get(eid, 0.0) for eid in veceids])

            # bottom and top in-plane stresses
            if vector.data.shape[2] > 3:
                Sxx = vector.data[:, :, 1][:, check]
                Syy = vector.data[:, :, 2][:, check]
                Sxy = vector.data[:, :, 3][:, check]
            else:
                Sxx = vector.data[:, :, 0][:, check]
                Syy = vector.data[:, :, 1][:, check]
                Sxy = vector.data[:, :, 2][:, check]
            if vector.data.dtype == np.complex64 or vector.data.dtype == np.complex128:
                Sxx_theta_real, Syy_theta_real, Sxy_theta_real = transf_Mohr(Sxx.real, Syy.real, Sxy.real, vecthetarad)
                Sxx_theta_imag, Syy_theta_imag, Sxy_theta_imag = transf_Mohr(Sxx.imag, Syy.imag, Sxy.imag, vecthetarad)
                tmp = np.zeros_like(new_vector.data[:, :, 0][:, check])
                tmp.real = Sxx_theta_real
                tmp.imag = Sxx_theta_imag
                new_vector.data[:, :, 0][:, check] = tmp
                tmp.real = Syy_theta_real
                tmp.imag = Syy_theta_imag
                new_vector.data[:, :, 1][:, check] = tmp
                tmp.real = Sxy_theta_real
                tmp.imag = Sxy_theta_imag
                new_vector.data[:, :, 2][:, check] = tmp
            else:
                Sxx_theta, Syy_theta, Sxy_theta = transf_Mohr(Sxx, Syy, Sxy, vecthetarad)
                new_vector.data[:, :, 1][:, check] = Sxx_theta
                new_vector.data[:, :, 2][:, check] = Syy_theta
                new_vector.data[:, :, 3][:, check] = Sxy_theta
                thetadeg_new = thetadeg_to_principal(Sxx_theta, Syy_theta, Sxy_theta)
                new_vector.data[:, :, 4][:, check] = thetadeg_new

            #TODO implement transformation for corner nodes
            #     for now we just zero the wrong values
            if 'quad8' in vecname:
                for i in [2, 3, 4, 5, 6, 7, 8, 9]:
                    new_vector.data[:, i, :] = 0

    for vecname in strain_vectors:
        op2_vectors = getattr(op2, vecname)
        new_vectors = getattr(op2_new, vecname)
        for subcase, vector in op2_vectors.items():
            new_vector = new_vectors[subcase]
            veceids = get_eids_from_op2_vector(vector)
            check = veceids != 0
            veceids = veceids[check]
            #NOTE assuming thetarad=0 for elements that exist in the op2 but
            #     not in the supplied bdf file
            vecthetarad = np.array([thetarad.get(eid, 0.0) for eid in veceids])

            # bottom and top in-plane strains
            if vector.data.shape[2] > 3:
                exx = vector.data[:, :, 1][:, check]
                eyy = vector.data[:, :, 2][:, check]
                exy = vector.data[:, :, 3][:, check] / 2.
            else:
                exx = vector.data[:, :, 0][:, check]
                eyy = vector.data[:, :, 1][:, check]
                exy = vector.data[:, :, 2][:, check] / 2.
            if vector.data.dtype == np.complex64 or vector.data.dtype == np.complex128:
                exx_theta_real, eyy_theta_real, exy_theta_real = transf_Mohr(exx.real, eyy.real, exy.real, vecthetarad)
                exx_theta_imag, eyy_theta_imag, exy_theta_imag = transf_Mohr(exx.imag, eyy.imag, exy.imag, vecthetarad)
                tmp = np.zeros_like(new_vector.data[:, :, 0][:, check])
                tmp.real = exx_theta_real
                tmp.imag = exx_theta_imag
                new_vector.data[:, :, 0][:, check] = tmp
                tmp.real = eyy_theta_real
                tmp.imag = eyy_theta_imag
                new_vector.data[:, :, 1][:, check] = tmp
                tmp.real = exy_theta_real * 2.
                tmp.imag = exy_theta_imag * 2.
                new_vector.data[:, :, 2][:, check] = tmp
            else:
                exx_theta, eyy_theta, exy_theta = transf_Mohr(exx, eyy, exy, vecthetarad)
                thetadeg_new = thetadeg_to_principal(exx_theta, eyy_theta, exy_theta)
                new_vector.data[:, :, 1][:, check] = exx_theta
                new_vector.data[:, :, 2][:, check] = eyy_theta
                new_vector.data[:, :, 3][:, check] = exy_theta * 2.
                new_vector.data[:, :, 4][:, check] = thetadeg_new

            #TODO implement transformation for corner nodes
            #     for now we just zero the wrong values
            if 'quad8' in vecname:
                for i in [2, 3, 4, 5, 6, 7, 8, 9]:
                    new_vector.data[:, i, :] = 0
    return op2_new

def main():  # pragma: no cover
    op2_filename = r'C:\NASA\m4\formats\git\pyNastran\models\solid_bending\solid_bending_coord1.op2'
    from pyNastran.op2.op2 import read_op2
    #from pyNastran.op2.op2_geom import read_op2_geom
    model = read_op2(op2_filename, load_geometry=True)
    cid = 0
    transform_solids(model, model, cid)

if __name__ == '__main__':  # pragma: no cover
    main()
