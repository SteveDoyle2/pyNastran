import os
import copy

import numpy as np
from numpy import cos, sin

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


def data_in_material_coord(bdf, op2, in_place=False):
    """

    Parameters
    ----------
    in_place : bool, optional
        If true the original op2 object is modified, otherwise a new one is
        created.
    """
    if in_place:
        op2_new = op2
    else:
        op2_new = copy.deepcopy(op2)
    for vecname in force_vectors:
        op2_vectors = getattr(op2, vecname)
        new_vectors = getattr(op2_new, vecname)
        for subcase, vector in op2_vectors.items():
            new_vector = new_vectors[subcase]
            eids = get_eids_from_op2_vector(vector)
            thetadeg = np.array([bdf.elements[eid].thetaMcid for eid in eids])
            thetarad = np.deg2rad(thetadeg)
            thick_array = np.array([bdf.elements[eid].pid.t for eid in eids])

            # membrane terms
            Sxx = vector.data[:, :, 0] / thick_array
            Syy = vector.data[:, :, 1] / thick_array
            Sxy = vector.data[:, :, 2] / thick_array
            Sxx_theta, Syy_theta, Sxy_theta = transf_Mohr(Sxx, Syy, Sxy, thetarad)
            new_vector.data[:, :, 0] = Sxx_theta * thick_array
            new_vector.data[:, :, 1] = Syy_theta * thick_array
            new_vector.data[:, :, 2] = Sxy_theta * thick_array

            # bending terms
            Sxx = vector.data[:, :, 3] / thick_array
            Syy = vector.data[:, :, 4] / thick_array
            Sxy = vector.data[:, :, 5] / thick_array
            Sxx_theta, Syy_theta, Sxy_theta = transf_Mohr(Sxx, Syy, Sxy, thetarad)
            new_vector.data[:, :, 3] = Sxx_theta * thick_array
            new_vector.data[:, :, 4] = Syy_theta * thick_array
            new_vector.data[:, :, 5] = Sxy_theta * thick_array

            # transverse terms
            Qx = vector.data[:, :, 6]
            Qy = vector.data[:, :, 7]
            Qx_new = cos(thetarad)*Qx + sin(thetarad)*Qy
            Qy_new = -sin(thetarad)*Qx + cos(thetarad)*Qy
            new_vector.data[:, :, 6] = Qx_new
            new_vector.data[:, :, 7] = Qy_new

    for vecname in stress_vectors:
        op2_vectors = getattr(op2, vecname)
        new_vectors = getattr(op2_new, vecname)
        for subcase, vector in op2_vectors.items():
            new_vector = new_vectors[subcase]
            eids = get_eids_from_op2_vector(vector)
            check = eids != 0
            eids = eids[check]
            thetadeg = np.array([bdf.elements[eid].thetaMcid for eid in eids])
            thetarad = np.deg2rad(thetadeg)

            # bottom and top in-plane stresses
            Sxx = vector.data[:, :, 1][:, check]
            Syy = vector.data[:, :, 2][:, check]
            Sxy = vector.data[:, :, 3][:, check]
            Sxx_theta, Syy_theta, Sxy_theta = transf_Mohr(Sxx, Syy, Sxy, thetarad)
            new_vector.data[:, :, 1][:, check] = Sxx_theta
            new_vector.data[:, :, 2][:, check] = Syy_theta
            new_vector.data[:, :, 3][:, check] = Sxy_theta
            thetadeg_new = thetadeg_to_principal(Sxx_theta, Syy_theta, Sxy_theta)
            new_vector.data[:, :, 4][:, check] = thetadeg_new
            if True:
                # TODO no need to recompute, keep this only in test routines
                Sxx_theta, Syy_theta, Sxy_theta = transf_Mohr(Sxx_theta, Syy_theta,
                        Sxy_theta, np.deg2rad(thetadeg_new))
                new_vector.data[:, :, 5][:, check] = Sxx_theta
                new_vector.data[:, :, 6][:, check] = Syy_theta


    if 0:
        #FIXME cannot figure out how to do the transformations for strain
        for vecname in strain_vectors:
            op2_vectors = getattr(op2, vecname)
            new_vectors = getattr(op2_new, vecname)
            for subcase, vector in op2_vectors.items():
                new_vector = new_vectors[subcase]
                eids = get_eids_from_op2_vector(vector)
                check = eids != 0
                eids = eids[check]
                modE = np.array([bdf.elements[eid].pid.mid.e for eid in eids])
                modG = np.array([bdf.elements[eid].pid.mid.g for eid in eids])
                nu = np.array([bdf.elements[eid].pid.mid.nu for eid in eids])
                thetadeg = np.array([bdf.elements[eid].thetaMcid for eid in eids])
                thetarad = np.deg2rad(thetadeg)

                # bottom and top in-plane strains
                Sxx = vector.data[:, :, 1][:, check] * modE
                Syy = vector.data[:, :, 2][:, check] * modE
                Sxy = vector.data[:, :, 3][:, check] * modG
                Sxx_theta, Syy_theta, Sxy_theta = transf_Mohr(Sxx, Syy, Sxy, thetarad)
                thetadeg_new = thetadeg_to_principal(Sxx_theta, Syy_theta, Sxy_theta)
                new_vector.data[:, :, 1][:, check] = Sxx_theta / modE
                new_vector.data[:, :, 2][:, check] = Syy_theta / modE
                new_vector.data[:, :, 3][:, check] = Sxy_theta / modG
                new_vector.data[:, :, 4][:, check] = thetadeg_new

    return op2_new

