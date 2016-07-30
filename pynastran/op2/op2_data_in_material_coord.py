import os
import copy

import numpy as np
from numpy import cos, sin, cross
from numpy.linalg import norm

from pyNastran.op2.op2 import OP2
from pyNastran.bdf.bdf import BDF
from pyNastran.bdf.test.bdf_unit_tests import Tester
from pyNastran.utils import integer_types

force_vectors = ['cquad4_force', 'cquad8_force', 'cquadr_force',
        'ctria3_force', 'ctria6_force', 'ctriar_force']
stress_vectors = ['cquad4_stress', 'cquad8_stress', 'cquadr_stress',
        'ctria3_stress', 'ctria6_stress', 'ctriar_stress']
strain_vectors = ['cquad4_strain', 'cquad8_strain', 'cquadr_strain',
        'ctria3_strain', 'ctria6_strain', 'ctriar_strain']


def transf_Mohr(Sxx, Syy, Sxy, thetarad):
    """Mohr's Circle-based Plane Stress Transformation

    Parameters
    ----------
    Sxx, Syy, Sxy : array-like
        Sigma_xx, Sigma_yy, Sigma_xy stresses.
    thetarad : array-like
        Array with angles for wich the stresses should be transformed.

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

            vector = op2.cqua4_force[1]
            vector = op2.cqua8_stress[1]
            vector = op2.ctriar_force[1]
            vector = op2.ctria3_stress[1]

    """
    eids = getattr(vector, 'element', None)
    if eids is None:
        eids = vector.element_node[:, 0]
    return eids


def is_mcid(elem):
    thetaMcid = getattr(elem, 'thetaMcid', None)
    if isinstance(thetaMcid, integer_types):
        return True
    else:
        return False


def check_theta(elem):
    theta = getattr(elem, 'thetaMcid', None)
    if theta is None:
        return 0.
    elif isinstance(theta, float):
        return theta
    elif isinstance(theta, integer_types):
        raise ValueError('MCID is accepted by this function')


def angle2vec(v1, v2):
    den = norm(v1, axis=1) * norm(v2, axis=1)
    return np.arccos((v1 * v2).sum(axis=1) / den)


def calc_imat(normals, csysi):
    jmat = cross(normals, csysi) # k x i
    jmat /= norm(jmat)
    imat = cross(jmat, normals)
    return imat


def data_in_material_coord(bdf, op2, in_place=False):
    """Convert OP2 2D element outputs to material coordinates

    Nastran allows the use of 'PARAM,OMID,YES' to print 2D element forces,
    stresses and strains based on the material direction. However, the
    convertion only takes place in the F06 output file, whereas the OP2 output
    file remains in the element coordinate system.

    This function converts the 2D element vectors to the material OP2
    similarly to most of the post-processing tools (Patran, Femap, HyperView,
    etc).

    Parameters
    ----------
    bdf : :class:`.BDF` object
        A :class:`.BDF` object that corresponds to the 'op2'.
    op2 : :class:`.OP2` object
        A :class:`.OP2` object that corresponds to the 'bdf'.
    in_place : bool, optional
        If true the original op2 object is modified, otherwise a new one is
        created.

    Returns
    -------
    op2_new : :class:`.OP2` object
        A :class:`.OP2` object with the abovementioned changes.

    Warning
    -------
     - doesn't handle composite stresses/strains/forces
     - doesn't handle solid stresses/strains/forces (e.g. MAT11)
     - zeros out data for CQUAD8s
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
            csysi = np.array([bdf.coords[e.thetaMcid].i for e in quadelems])
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
            csysi = np.array([bdf.coords[e.thetaMcid].i for e in triaelems])
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
            vecthetarad = np.array([thetarad[eid] for eid in veceids])

            # membrane terms
            Sxx = vector.data[:, :, 0]
            Syy = vector.data[:, :, 1]
            Sxy = vector.data[:, :, 2]
            Sxx_theta, Syy_theta, Sxy_theta = transf_Mohr(Sxx, Syy, Sxy, vecthetarad)
            new_vector.data[:, :, 0] = Sxx_theta
            new_vector.data[:, :, 1] = Syy_theta
            new_vector.data[:, :, 2] = Sxy_theta

            # bending terms
            Sxx = vector.data[:, :, 3]
            Syy = vector.data[:, :, 4]
            Sxy = vector.data[:, :, 5]
            Sxx_theta, Syy_theta, Sxy_theta = transf_Mohr(Sxx, Syy, Sxy, vecthetarad)
            new_vector.data[:, :, 3] = Sxx_theta
            new_vector.data[:, :, 4] = Syy_theta
            new_vector.data[:, :, 5] = Sxy_theta

            # transverse terms
            Qx = vector.data[:, :, 6]
            Qy = vector.data[:, :, 7]
            Qx_new = cos(vecthetarad)*Qx + sin(vecthetarad)*Qy
            Qy_new = -sin(vecthetarad)*Qx + cos(vecthetarad)*Qy
            new_vector.data[:, :, 6] = Qx_new
            new_vector.data[:, :, 7] = Qy_new

            #TODO implement transformation for corner nodes
            #     for now we just zero the wrong values
            if 'quad8' in vecname:
                for i in [1, 2, 3, 4]:
                    new_vector.data[:, i, :] = 0


    for vecname in stress_vectors:
        op2_vectors = getattr(op2, vecname)
        new_vectors = getattr(op2_new, vecname)
        for subcase, vector in op2_vectors.items():
            new_vector = new_vectors[subcase]
            veceids = get_eids_from_op2_vector(vector)
            check = veceids != 0
            veceids = veceids[check]
            vecthetarad = np.array([thetarad[eid] for eid in veceids])

            # bottom and top in-plane stresses
            Sxx = vector.data[:, :, 1][:, check]
            Syy = vector.data[:, :, 2][:, check]
            Sxy = vector.data[:, :, 3][:, check]
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
            vecthetarad = np.array([thetarad[eid] for eid in veceids])

            # bottom and top in-plane strains
            exx = vector.data[:, :, 1][:, check]
            eyy = vector.data[:, :, 2][:, check]
            exy = vector.data[:, :, 3][:, check] / 2.
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

