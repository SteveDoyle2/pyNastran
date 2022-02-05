"""
Defines:
 - data_in_material_coord(bdf, op2, in_place=False)

"""
from __future__ import annotations
import copy
from typing import Tuple, Dict, Union, TYPE_CHECKING

import numpy as np
from numpy import cos, sin, cross
from numpy.linalg import norm  # type: ignore

from pyNastran.utils.numpy_utils import integer_types

if TYPE_CHECKING:  # pragma: no cover
    from cpylog import SimpleLogger
    from pyNastran.bdf.bdf import BDF
    from pyNastran.op2.op2 import OP2
    from pyNastran.op2.tables.oes_stressStrain.real.oes_plates import RealPlateStressArray, RealPlateStrainArray
    from pyNastran.op2.tables.oef_forces.oef_force_objects import RealPlateForceArray, RealPlateBilinearForceArray
    #from pyNastran.op2.tables.oef_forces.oef_complex_force_objects import ComplexPlateForceArray, ComplexPlate2ForceArray

force_vectors = ['cquad4_force', 'cquad8_force', 'cquadr_force',
                 'ctria3_force', 'ctria6_force', 'ctriar_force']
stress_vectors = ['cquad4_stress', 'cquad8_stress', 'cquadr_stress',
                  'ctria3_stress', 'ctria6_stress', 'ctriar_stress']
strain_vectors = ['cquad4_strain', 'cquad8_strain', 'cquadr_strain',
                  'ctria3_strain', 'ctria6_strain', 'ctriar_strain']

def transf_Mohr(Sxx: np.ndarray,
                Syy: np.ndarray,
                Sxy: np.ndarray,
                theta_rad: np.ndarray) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
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
    theta_rad = np.asarray(theta_rad)
    Scenter = (Sxx + Syy) / 2.
    R = np.sqrt((Sxx - Scenter)**2 + Sxy**2)
    theta_rad_Mohr = np.arctan2(-Sxy, Sxx - Scenter) + 2*theta_rad
    cos_Mohr = cos(theta_rad_Mohr)
    Sxx_theta = Scenter + R*cos_Mohr
    Syy_theta = Scenter - R*cos_Mohr
    Sxy_theta = -R*sin(theta_rad_Mohr)
    return Sxx_theta, Syy_theta, Sxy_theta


def theta_deg_to_principal(Sxx: np.ndarray,
                           Syy: np.ndarray,
                           Sxy: np.ndarray) -> np.ndarray:
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
    Scenter = (Sxx + Syy) / 2.
    thetarad = np.arctan2(Sxy, Scenter - Syy)
    return np.rad2deg(thetarad) / 2.


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

def angle2vec(v1: np.ndarray, v2: np.ndarray) -> np.ndarray:
    """
    Using the definition of the dot product to get the angle

    v1 o v2 = |v1| * |v2| * cos(theta)
    theta = np.arccos( (v1 o v2) / (|v1|*|v2|))

    """
    denom = norm(v1, axis=1) * norm(v2, axis=1)
    theta = np.arccos((v1 * v2).sum(axis=1) / denom)
    return theta


def calc_imat(normals: np.ndarray, csysi: np.ndarray) -> np.ndarray:
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


def data_in_material_coord(bdf: BDF, op2: OP2,
                           in_place: bool=False,
                           debug: bool=False) -> OP2:
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
    log = bdf.log

    #eid_to_theta_rad = get_eid_to_theta_rad(bdf, debug)
    eid_to_theta_rad = get_eid_to_theta_rad2(bdf, debug)
    #assert len(eid_to_theta_rad) == len(eid_to_theta_rad2)
    #is_failed = False
    #for eid, theta1 in eid_to_theta_rad.items():
        #theta2 = eid_to_theta_rad2[eid]
        #theta1_deg = np.degrees(theta1)
        #theta2_deg = np.degrees(theta2)
        #if not np.allclose(theta1_deg, theta2_deg):
            #elem = bdf.elements[eid]
            #log.warning(f'eid={eid} {elem.type} theta/mcid={elem.theta_mcid} theta1={theta1_deg} theta2={theta2_deg}')
            #is_failed = True
    #if is_failed:
        #x = 1

    if in_place:
        op2_new = op2
    else:
        op2_new = copy.deepcopy(op2)
        op2_new.log = log

    for vec_name in force_vectors:
        op2_vectors = getattr(op2, vec_name)
        new_vectors = getattr(op2_new, vec_name)
        for subcase, vector in op2_vectors.items():
            new_vector = new_vectors[subcase]
            _transform_shell_force(
                vec_name,
                vector, new_vector,
                eid_to_theta_rad, log)
            if new_vector.data_frame is not None:
                new_vector.build_dataframe()

    for vec_name in stress_vectors:
        op2_vectors = getattr(op2, vec_name)
        new_vectors = getattr(op2_new, vec_name)
        for subcase, vector in op2_vectors.items():
            new_vector = new_vectors[subcase]
            _transform_shell_stress(
                vec_name,
                vector, new_vector,
                eid_to_theta_rad, log)
            if new_vector.data_frame is not None:
                new_vector.build_dataframe()

    for vec_name in strain_vectors:
        op2_vectors = getattr(op2, vec_name)
        new_vectors = getattr(op2_new, vec_name)
        for subcase, vector in op2_vectors.items():
            new_vector = new_vectors[subcase]
            _transform_shell_strain(
                vec_name,
                vector, new_vector,
                eid_to_theta_rad, log)
            if new_vector.data_frame is not None:
                new_vector.build_dataframe()
    return op2_new

def get_eid_to_theta_rad(bdf: BDF, debug: bool) -> Dict[int, float]:
    eids = np.array(list(bdf.elements.keys()))
    elems = np.array(list(bdf.elements.values()))
    mcid = np.array([is_mcid(e) for e in elems])
    elems_mcid = elems[mcid]
    elems_theta = elems[~mcid]

    theta_deg = np.full(elems.shape, np.nan, dtype='float64')
    theta_deg[~mcid] = np.array([check_theta(e) for e in elems_theta])

    theta_rad = np.deg2rad(theta_deg)
    log = bdf.log
    if debug:
        log.info(f'eids = {eids}')
        log.info(f'mcid = {mcid}')
        log.info(f'theta_deg = {theta_deg}')

    #NOTE separating quad types to get vectorizable "corner"
    cquad_types = ['CQUAD4', 'CQUAD8', 'CQUADR']
    for cquad_type in cquad_types:
        #if cquad_type not in bdf.card_count:
            #continue

        # elems with THETA
        this_quad = np.array([cquad_type == e.type for e in elems_theta])
        if not np.any(this_quad):
            continue
        quad_elems = elems_theta[this_quad]
        corner = np.array([e.get_node_positions() for e in quad_elems])
        g1 = corner[:, 0, :]
        g2 = corner[:, 1, :]
        g3 = corner[:, 2, :]
        g4 = corner[:, 3, :]
        beta_rad = angle2vec(g3 - g1, g2 - g1)
        gamma_rad = angle2vec(g4 - g2, g1 - g2)
        alpha_rad = (beta_rad + gamma_rad) / 2.
        tmp = theta_rad[~mcid]
        tmp[this_quad] += alpha_rad - beta_rad
        theta_rad[~mcid] = tmp

        # elems with MCID
        this_quad = np.array([cquad_type in e.type for e in elems_mcid])
        if not np.any(this_quad):
            continue
        quad_elems = elems_mcid[this_quad]
        corner = np.array([e.get_node_positions() for e in quad_elems])
        g1 = corner[:, 0, :]
        g2 = corner[:, 1, :]
        g3 = corner[:, 2, :]
        g4 = corner[:, 3, :]
        normals2 = np.array([e.Normal() for e in quad_elems])
        normals = _get_normal(g3 - g1, g4 - g2)
        assert np.allclose(normals, normals2)

        csysi = np.array([bdf.coords[e.theta_mcid].i for e in quad_elems])
        imat = calc_imat(normals, csysi)
        tmp = theta_rad[mcid]
        tmp[this_quad] = angle2vec(g2 - g1, imat)
        # getting sign of THETA
        check_normal = cross(g2 - g1, imat)
        tmp[this_quad] *= np.sign((check_normal * normals).sum(axis=1))
        beta_rad = angle2vec(g3 - g1, g2 - g1)
        gamma_rad = angle2vec(g4 - g2, g1 - g2)
        alpha_rad = (beta_rad + gamma_rad) / 2.
        tmp[this_quad] += alpha_rad - beta_rad
        theta_rad[mcid] = tmp

    ctria_types = ['CTRIA3', 'CTRIA6', 'CTRIAR']
    for ctria_type in ctria_types:
        #if ctria_type not in bdf.card_count:
            #continue
        # elems with MCID
        this_tria = np.array([ctria_type == e.type for e in elems_mcid])
        if debug:
            log.debug(f'found {ctria_type} with thistria={this_tria}')

        if not np.any(this_tria):
            continue
        tria_elems = elems_mcid[this_tria]
        #eids = [elem.eid for elem in triaelems]

        # corner: (nelments, 6, 3) for a CTRIA6
        corner = np.array([e.get_node_positions() for e in tria_elems])
        g1 = corner[:, 0, :]
        g2 = corner[:, 1, :]
        g3 = corner[:, 2, :]
        normals2 = np.array([e.Normal() for e in tria_elems])
        normals = _get_normal(g2 - g1, g3 - g1)
        assert np.allclose(normals, normals2)

        csysi = np.array([bdf.coords[e.theta_mcid].i for e in tria_elems])
        imat = calc_imat(normals, csysi)
        tmp = theta_rad[mcid]
        tmp[this_tria] = angle2vec(g2 - g1, imat)

        # getting sign of THETA
        check_normal = cross(g2 - g1, imat)
        tmp[this_tria] *= np.sign((check_normal * normals).sum(axis=1))
        theta_rad[mcid] = tmp

    eid_to_theta_rad = dict([[eid, theta] for eid, theta in zip(eids, theta_rad)])
    return eid_to_theta_rad

def _get_tri_nodes(nids: np.ndarray,
                   xyz_cid0: np.ndarray,
                   element_nodes: np.ndarray) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    inode = np.searchsorted(nids, element_nodes)
    g1 = xyz_cid0[inode[:, 0], :]
    g2 = xyz_cid0[inode[:, 1], :]
    g3 = xyz_cid0[inode[:, 2], :]
    return g1, g2, g3

def _get_quad_nodes(nids: np.ndarray,
                    xyz_cid0: np.ndarray,
                    element_nodes: np.ndarray) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    inode = np.searchsorted(nids, element_nodes)
    g1 = xyz_cid0[inode[:, 0], :]
    g2 = xyz_cid0[inode[:, 1], :]
    g3 = xyz_cid0[inode[:, 2], :]
    g4 = xyz_cid0[inode[:, 3], :]
    return g1, g2, g3, g4

def get_eid_to_theta_rad2(model: BDF, debug: bool) -> Dict[int, float]:
    coords = model.coords
    nid_cp_cd, xyz_cid0, *unused = model.get_xyz_in_coord_array(
        cid=0, fdtype='float64', idtype='int32')
    nids = nid_cp_cd[:, 0]
    #eids = []
    #theta_deg = []

    #eid_theta_ = []
    #val_theta_ = []
    eid_tri_mcid_ = []
    val_tri_mcid_ = []
    eid_quad_mcid_ = []
    val_quad_mcid_ = []

    eid_tri_theta_ = []
    val_tri_theta_ = []
    eid_quad_theta_ = []
    val_quad_theta_ = []

    eid_to_theta_rad = {}
    #etype_to_mcid = {}

    node_tri_mcid_ = []
    node_quad_mcid_ = []

    node_tri_theta_ = []
    node_quad_theta_ = []
    for eid, elem in model.elements.items():
        etype = elem.type
        if etype in {'CQUAD4', 'CQUAD8', 'CQUADR'}:
            theta_mcid = elem.theta_mcid
            if isinstance(theta_mcid, float):
                #eids.append(eid)
                #theta_deg.append(theta_mcid)
                eid_quad_theta_.append(eid)
                val_quad_theta_.append(theta_mcid)
                node_quad_theta_.append(elem.nodes[:4])
                continue
            eid_quad_mcid_.append(eid)
            val_quad_mcid_.append(theta_mcid)
            node_quad_mcid_.append(elem.nodes[:4])
        elif etype in {'CTRIA3', 'CTRIA6', 'CTRIAR'}:
            theta_mcid = elem.theta_mcid
            if isinstance(theta_mcid, float):
                #eids.append(eid)
                #theta_deg.append(theta_mcid)
                eid_tri_theta_.append(eid)
                val_tri_theta_.append(theta_mcid)
                node_tri_theta_.append(elem.nodes[:3])
                continue
            eid_tri_mcid_.append(eid)
            val_tri_mcid_.append(theta_mcid)
            node_tri_mcid_.append(elem.nodes[:3])
        continue
    #theta_rad = np.radians(theta_deg)

    # put in dictionary form- maybe change later
    #eid_to_theta_rad = {eid: theta for eid, theta in zip(eids, theta_rad)}


    if len(eid_tri_theta_):
        # the triangle theta is defined from the line g1-g2
        eid_tri_theta = np.array(eid_tri_theta_, dtype='int32')
        val_tri_theta_deg = np.array(val_tri_theta_, dtype='float64')
        #node_tri_theta = np.array(node_tri_theta_, dtype='int32')
        tri_theta_rad = np.radians(val_tri_theta_deg)

        for eid, theta in zip(eid_tri_theta, tri_theta_rad):
            eid_to_theta_rad[eid] = theta

    if len(eid_quad_theta_):
        # the quad theta hs
        # For CQUAD4 elements which are not p-elements and not hyperelastic,
        # the reference coordinate system is the default for output is the
        # element coordinate system.
        eid_quad_theta = np.array(eid_quad_theta_, dtype='int32')
        val_quad_theta_deg = np.array(val_quad_theta_, dtype='float64')
        node_quad_theta = np.array(node_quad_theta_, dtype='int32')
        quad_theta_rad = np.radians(val_quad_theta_deg)

        g1, g2, g3, g4 = _get_quad_nodes(nids, xyz_cid0, node_quad_theta)
        g21 = g2 - g1
        g12 = -g21

        beta_rad = angle2vec(g3 - g1, g21)
        gamma_rad = angle2vec(g4 - g2, g12)
        alpha_rad = (beta_rad + gamma_rad) / 2.
        theta_radi = alpha_rad - beta_rad
        for eid, theta in zip(eid_quad_theta, quad_theta_rad+theta_radi):
            eid_to_theta_rad[eid] = theta


    if len(eid_quad_mcid_):
        eid_quad_mcid = np.array(eid_quad_mcid_, dtype='int32')
        val_quad_mcid = np.array(val_quad_mcid_, dtype='int32')
        node_quad_mcid = np.array(node_quad_mcid_, dtype='int32')

        g1, g2, g3, g4 = _get_quad_nodes(nids, xyz_cid0, node_quad_mcid)
        g42 = g4 - g2
        g31 = g3 - g1
        g21 = g2 - g1
        g12 = -g21
        normals = _get_normal(g31, g42)

        csysi = np.array([coords[mcid].i for mcid in val_quad_mcid])
        imat = calc_imat(normals, csysi)

        theta_radi = angle2vec(g21, imat)
        # getting sign of THETA
        check_normal = cross(g21, imat)
        theta_radi *= np.sign((check_normal * normals).sum(axis=1))
        beta_rad = angle2vec(g31, g21)
        gamma_rad = angle2vec(g42, g12)
        alpha_rad = (beta_rad + gamma_rad) / 2.
        theta_radi += alpha_rad - beta_rad
        for eid, theta in zip(eid_quad_mcid, theta_radi):
            eid_to_theta_rad[eid] = theta

    if len(eid_tri_mcid_):
        eid_tri_mcid = np.array(eid_tri_mcid_, dtype='int32')
        val_tri_mcid = np.array(val_tri_mcid_, dtype='int32')
        node_tri_mcid = np.array(node_tri_mcid_, dtype='int32')
        #eids.append(eid_tri_mcid)

        g1, g2, g3 = _get_tri_nodes(nids, xyz_cid0, node_tri_mcid)
        g21 = g2 - g1
        normals = _get_normal(g21, g3 - g1)

        csysi = np.array([coords[mcid].i for mcid in val_tri_mcid])
        imat = calc_imat(normals, csysi)

        theta_radi = angle2vec(g21, imat)
        # getting sign of THETA
        check_normal = cross(g21, imat)
        theta_radi *= np.sign((check_normal * normals).sum(axis=1))
        for eid, theta in zip(eid_tri_mcid, theta_radi):
            eid_to_theta_rad[eid] = theta
    return eid_to_theta_rad



def _get_normal(v1: np.ndarray, v2: np.ndarray) -> np.ndarray:
    """Cross product and divide by the magnitude so it's a unit normal"""
    normals = np.cross(v1, v2)
    normals /= np.linalg.norm(normals, axis=1)[:, np.newaxis]
    return normals

def _transform_shell_force(vec_name: str,
                           vector: Union[RealPlateForceArray, RealPlateBilinearForceArray],
                           new_vector: Union[RealPlateForceArray, RealPlateBilinearForceArray],
                           eid_to_theta_rad: Dict[int, float],
                           log: SimpleLogger):
    vec_eids = get_eids_from_op2_vector(vector)
    #NOTE assuming thetarad=0 for elements that exist in the op2 but
    #     not in the supplied bdf file
    vec_theta_rad = np.array([eid_to_theta_rad.get(eid, 0.0) for eid in vec_eids])

    if vec_eids.shape[0] == vector.data.shape[1] // 5:
        #log.debug(f'A: vec_eids.shape={vec_eids.shape} vector.data.shape={vector.data.shape}')
        steps = [5, 5, 5, 5, 5]
    else: # assuming  always that veceids.shape[0] == vector.data.shape[1]
        #log.debug(f'B: vec_eids.shape={vec_eids.shape} vector.data.shape={vector.data.shape}')
        steps = [1]

    for start, step in enumerate(steps):
        slicei = slice(start, vector.data.shape[1], step)
        # membrane terms
        Sxx = vector.data[:, slicei, 0]
        Syy = vector.data[:, slicei, 1]
        Sxy = vector.data[:, slicei, 2]
        if vector.data.dtype == np.complex64 or vector.data.dtype == np.complex128:
            Sxx_theta_real, Syy_theta_real, Sxy_theta_real = transf_Mohr(
                Sxx.real, Syy.real, Sxy.real, vec_theta_rad)
            new_vector.data[:, slicei, 0].real = Sxx_theta_real
            new_vector.data[:, slicei, 1].real = Syy_theta_real
            new_vector.data[:, slicei, 2].real = Sxy_theta_real
            Sxx_theta_imag, Syy_theta_imag, Sxy_theta_imag = transf_Mohr(
                Sxx.imag, Syy.imag, Sxy.imag, vec_theta_rad)
            new_vector.data[:, slicei, 0].imag = Sxx_theta_imag
            new_vector.data[:, slicei, 1].imag = Syy_theta_imag
            new_vector.data[:, slicei, 2].imag = Sxy_theta_imag
        else:
            Sxx_theta, Syy_theta, Sxy_theta = transf_Mohr(Sxx, Syy, Sxy, vec_theta_rad)
            new_vector.data[:, slicei, 0] = Sxx_theta
            new_vector.data[:, slicei, 1] = Syy_theta
            new_vector.data[:, slicei, 2] = Sxy_theta

        # bending terms
        Sxx = vector.data[:, slicei, 3]
        Syy = vector.data[:, slicei, 4]
        Sxy = vector.data[:, slicei, 5]
        if vector.data.dtype == np.complex64 or vector.data.dtype == np.complex128:
            Sxx_theta_real, Syy_theta_real, Sxy_theta_real = transf_Mohr(
                Sxx.real, Syy.real, Sxy.real, vec_theta_rad)
            new_vector.data[:, slicei, 3].real = Sxx_theta_real
            new_vector.data[:, slicei, 4].real = Syy_theta_real
            new_vector.data[:, slicei, 5].real = Sxy_theta_real
            Sxx_theta_imag, Syy_theta_imag, Sxy_theta_imag = transf_Mohr(
                Sxx.imag, Syy.imag, Sxy.imag, vec_theta_rad)
            new_vector.data[:, slicei, 3].imag = Sxx_theta_imag
            new_vector.data[:, slicei, 4].imag = Syy_theta_imag
            new_vector.data[:, slicei, 5].imag = Sxy_theta_imag

        else:
            Sxx_theta, Syy_theta, Sxy_theta = transf_Mohr(Sxx, Syy, Sxy, vec_theta_rad)
            new_vector.data[:, slicei, 3] = Sxx_theta
            new_vector.data[:, slicei, 4] = Syy_theta
            new_vector.data[:, slicei, 5] = Sxy_theta

        # transverse terms
        Qx = vector.data[:, slicei, 6]
        Qy = vector.data[:, slicei, 7]
        cos_theta = cos(vec_theta_rad)
        sin_theta = sin(vec_theta_rad)
        if vector.data.dtype == np.complex64 or vector.data.dtype == np.complex128:
            Qx_new_real = cos_theta*Qx.real + sin_theta*Qy.real
            Qy_new_real = -sin_theta*Qx.real + cos_theta*Qy.real
            new_vector.data[:, slicei, 6].real = Qx_new_real
            new_vector.data[:, slicei, 7].real = Qy_new_real

            Qx_new_imag = cos_theta*Qx.imag + sin_theta*Qy.imag
            Qy_new_imag = -sin_theta*Qx.imag + cos_theta*Qy.imag
            new_vector.data[:, slicei, 6].imag = Qx_new_imag
            new_vector.data[:, slicei, 7].imag = Qy_new_imag
        else:
            Qx_new = cos_theta*Qx + sin_theta*Qy
            Qy_new = -sin_theta*Qx + cos_theta*Qy
            new_vector.data[:, slicei, 6] = Qx_new
            new_vector.data[:, slicei, 7] = Qy_new

    #TODO implement transformation for corner nodes
    #     for now we just zero the wrong values
    if 'quad8' in vec_name:
        for j in [1, 2, 3, 4]:
            new_vector.data[:, j, :] = 0

def _transform_shell_stress(
        vec_name: str,
        vector: RealPlateStressArray,
        new_vector: RealPlateStressArray,
        eid_to_theta_rad: Dict[int, float],
        log: SimpleLogger):
    vec_eids = get_eids_from_op2_vector(vector)
    check = (vec_eids != 0)
    if np.any(~check):
        log.warning(f'why are there 0 element ids?  vec_eids={vec_eids}')
    vec_eids = vec_eids[check]
    #NOTE assuming thetarad=0 for elements that exist in the op2 but
    #     not in the supplied bdf file
    vec_theta_rad = np.array([eid_to_theta_rad.get(eid, 0.0) for eid in vec_eids])

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
        Sxx_theta_real, Syy_theta_real, Sxy_theta_real = transf_Mohr(Sxx.real, Syy.real, Sxy.real, vec_theta_rad)
        Sxx_theta_imag, Syy_theta_imag, Sxy_theta_imag = transf_Mohr(Sxx.imag, Syy.imag, Sxy.imag, vec_theta_rad)
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
        Sxx_theta, Syy_theta, Sxy_theta = transf_Mohr(Sxx, Syy, Sxy, vec_theta_rad)
        new_vector.data[:, :, 1][:, check] = Sxx_theta
        new_vector.data[:, :, 2][:, check] = Syy_theta
        new_vector.data[:, :, 3][:, check] = Sxy_theta
        theta_deg_new = theta_deg_to_principal(Sxx_theta, Syy_theta, Sxy_theta)
        new_vector.data[:, :, 4][:, check] = theta_deg_new

    #TODO implement transformation for corner nodes
    #     for now we just zero the wrong values
    if 'quad8' in vec_name:
        for i in [2, 3, 4, 5, 6, 7, 8, 9]:
            new_vector.data[:, i, :] = 0

def _transform_shell_strain(
        vec_name: str,
        vector: RealPlateStrainArray,
        new_vector: RealPlateStrainArray,
        eid_to_theta_rad: Dict[int, float],
        log: SimpleLogger):
    vec_eids = get_eids_from_op2_vector(vector)
    check = (vec_eids != 0)
    if np.any(~check):
        log.warning(f'why are there 0 element ids?  vec_eids={vec_eids}')

    vec_eids = vec_eids[check]
    #NOTE assuming thetarad=0 for elements that exist in the op2 but
    #     not in the supplied bdf file
    theta_rad = np.array([eid_to_theta_rad.get(eid, 0.0) for eid in vec_eids])

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
        exx_theta_real, eyy_theta_real, exy_theta_real = transf_Mohr(exx.real, eyy.real, exy.real, theta_rad)
        exx_theta_imag, eyy_theta_imag, exy_theta_imag = transf_Mohr(exx.imag, eyy.imag, exy.imag, theta_rad)
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
        exx_theta, eyy_theta, exy_theta = transf_Mohr(exx, eyy, exy, theta_rad)
        theta_deg_new = theta_deg_to_principal(exx_theta, eyy_theta, exy_theta) # TODO: can we just add dtheta?
        new_vector.data[:, :, 1][:, check] = exx_theta
        new_vector.data[:, :, 2][:, check] = eyy_theta
        new_vector.data[:, :, 3][:, check] = exy_theta * 2.
        new_vector.data[:, :, 4][:, check] = theta_deg_new

    #TODO implement transformation for corner nodes
    #     for now we just zero the wrong values
    if 'quad8' in vec_name:
        for i in [2, 3, 4, 5, 6, 7, 8, 9]:
            new_vector.data[:, i, :] = 0
