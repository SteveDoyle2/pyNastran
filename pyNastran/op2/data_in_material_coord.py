import os
import copy

import numpy as np
from numpy import cos, sin
from numpy.linalg import norm

from pyNastran.op2.op2 import OP2
from pyNastran.bdf.bdf import BDF
from pyNastran.bdf.test.bdf_unit_tests import Tester

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
    R = ((Sxx - Scenter)**2 + Sxy**2)**0.5
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


def check_theta(theta):
    if theta is None:
        return 0.
    elif isinstance(theta, int):
        raise NotImplementedError('MCID currently not supported')
    elif isinstance(theta, float):
        return theta


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

    """
    if in_place:
        op2_new = op2
    else:
        op2_new = copy.deepcopy(op2)

    eids = np.array(list(bdf.elements.keys()))
    elems = np.array(list(bdf.elements.values()))
    thetadeg = np.array([check_theta(e.thetaMcid) for e in elems])
    thetarad = np.deg2rad(thetadeg)

    quad_types = ['CQUAD4', 'CQUAD8', 'CQUADR']
    for quad_type in quad_types:
        isquad = np.array([quad_type in e.type for e in elems])
        if not np.any(isquad):
            continue

        quadelems = elems[isquad]

        coords = np.array([e.get_node_positions() for e in quadelems])
        g1 = coords[:, 0, :]
        g2 = coords[:, 1, :]
        g3 = coords[:, 2, :]
        g4 = coords[:, 3, :]

        v1 = g3 - g1
        v2 = g2 - g1
        den = norm(v1, axis=1) * norm(v2, axis=1)
        betarad = np.arccos((v1 * v2).sum(axis=1) / den)

        v1 = g4 - g2
        v2 = g1 - g2
        den = norm(v1, axis=1) * norm(v2, axis=1)
        gammarad = np.arccos((v1 * v2).sum(axis=1) / den)

        alpharad = (betarad + gammarad) / 2.
        thetarad[isquad] -= betarad
        thetarad[isquad] += alpharad

    thetarad = dict([[eid, theta] for eid, theta in zip(eids, thetarad)])
    thick_array = dict([[eid, e.pid.t] for eid, e, in zip(eids, elems)])

    for vecname in force_vectors:
        op2_vectors = getattr(op2, vecname)
        new_vectors = getattr(op2_new, vecname)
        for subcase, vector in op2_vectors.items():
            new_vector = new_vectors[subcase]
            veceids = get_eids_from_op2_vector(vector)
            vecthick = np.array([thick_array[eid] for eid in veceids])
            vecthetarad = np.array([thetarad[eid] for eid in veceids])

            # membrane terms
            Sxx = vector.data[:, :, 0] / vecthick
            Syy = vector.data[:, :, 1] / vecthick
            Sxy = vector.data[:, :, 2] / vecthick
            Sxx_theta, Syy_theta, Sxy_theta = transf_Mohr(Sxx, Syy, Sxy, vecthetarad)
            new_vector.data[:, :, 0] = Sxx_theta * vecthick
            new_vector.data[:, :, 1] = Syy_theta * vecthick
            new_vector.data[:, :, 2] = Sxy_theta * vecthick

            # bending terms
            Sxx = vector.data[:, :, 3] / vecthick
            Syy = vector.data[:, :, 4] / vecthick
            Sxy = vector.data[:, :, 5] / vecthick
            Sxx_theta, Syy_theta, Sxy_theta = transf_Mohr(Sxx, Syy, Sxy, vecthetarad)
            new_vector.data[:, :, 3] = Sxx_theta * vecthick
            new_vector.data[:, :, 4] = Syy_theta * vecthick
            new_vector.data[:, :, 5] = Sxy_theta * vecthick

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

