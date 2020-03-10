# coding: utf-8
# pylint: disable=C0103
"""
Defines:
  - mass_poperties
      get the mass & moment of inertia of the model

"""
from itertools import count
from collections import defaultdict
from numpy import array, cross, dot
from numpy.linalg import norm  # type: ignore
import numpy as np
#from pyNastran.bdf.cards.materials import get_mat_props_S
from pyNastran.utils.numpy_utils import integer_types
from pyNastran.utils.mathematics import integrate_positive_unit_line

NO_MASS = {
    # has mass
    'CCONEAX',

    # no mass
    'GRID', 'PARAM', 'FORCE', 'FORCE1', 'FORCE2', 'MOMENT1', 'MOMENT2', 'LOAD',
    'DVPREL1', 'DVPREL2', 'DVCREL1', 'DVCREL2', 'DVMREL1', 'DVMREL2', 'DCONSTR', 'DESVAR',
    'DEQATN', 'DRESP1', 'DRESP2', 'DRESP3',
    'SPC', 'SPC1', 'SPCADD', 'MPC', 'MPCADD',
    'MAT1', 'MAT2', 'MAT4', 'MAT5', 'MAT8', 'MAT10', 'MAT11', 'MAT3D', 'CREEP',
    'MATT1', 'MATT3',

    'PELAS', 'PVISC', 'PBUSH1D',
    'CELAS1', 'CELAS2', 'CELAS3', 'CELAS4', #'CLEAS5',
    'CDAMP1', 'CDAMP2', 'CDAMP3', 'CDAMP4', 'CDAMP5',
    'PDAMP', 'PGAP',
    'CBUSH', 'CBUSH1D', 'CBUSH2D', 'CVISC', 'CGAP', # is this right?
    'CFAST', 'GENEL',
    'CRAC2D', 'CRAC3D',

    'CSSCHD', 'CAERO1', 'CAERO2', 'CAERO3', 'CAERO4', 'CAERO5',
    'CBARAO', 'CORD1R', 'CORD2R', 'CORD1C', 'CORD2C', 'CORD1S', 'CORD2S',
    'CORD3G', 'CONV', 'CONVM', 'CLOAD',
    'CHBDYG', 'CHBDYE', 'CHBDYP', 'TEMP', 'TEMPD', 'QVECT',

    'CTRAX3', 'CTRAX6', 'CQUADX8', 'CQUADX4',
    'CPLSTN3', 'CPLSTN6', 'CPLSTN4', 'CPLSTN8',

    'ASET', 'ASET1', 'BSET', 'BSET1', 'CSET', 'CSET1',
    'QSET', 'QSET1', 'USET', 'USET1', 'OMIT', 'OMIT1',

    'DLOAD', 'TLOAD1', 'PLOAD', 'PLOAD2', 'PLOAD4',
    'TSTEP', 'TSTEPNL', 'TABLED1', 'TABLED2', 'TABLED3', 'TABLED4',
    'TABLEM1', 'TABLEM2', 'TABLEM3', 'TABLEM4', 'TABLEST',

    # aero
    'MONPNT1', 'MONPNT2', 'MONPNT3', 'MONDSP1',
    'AERO', 'AEROS',
    'CAERO1', 'CAERO2', 'CAERO3', 'CAERO4', 'CAERO5',
    'SPLINE1', 'SPLINE2', 'SPLINE3', 'SPLINE4', 'SPLINE5', 'SPLINE6', 'SPLINE7',
    'AEPARM', 'AEFACT', 'AESURF', 'AESURFS', 'AELINK',
    'CSSCHD', 'TRIM', 'TRIM2', 'DIVERG', 'FLUTTER', 'GUST',

    'DVPREL1', 'DVPREL2', 'DVMREL1', 'DVMREL2', 'DVCREL1', 'DVCREL2',
    'DESVAR', 'DCONADD', 'DRESP1', 'DRESP2', 'DRESP3', 'DEQATN', 'DSCREEN',
    'SUPORT', 'SUPORT1',
    'CYJOIN',

    # acoustic
    'CHACAB', 'CAABSF',
}

def transform_inertia(mass, xyz_cg, xyz_ref, xyz_ref2, I_ref):
    """
    Transforms mass moment of inertia using parallel-axis theorem.

    Parameters
    ----------
    mass : float
        the mass
    xyz_cg : (3, ) float ndarray
        the CG location
    xyz_ref : (3, ) float ndarray
        the original reference location
    xyz_ref2 : (3, ) float ndarray
        the new reference location
    I_ref : (6, ) float ndarray
        the mass moment of inertias about the original reference point
        [Ixx, Iyy, Izz, Ixy, Ixz, Iyz]

    Returns
    -------
    I_new : (6, ) float ndarray
        the mass moment of inertias about the new reference point
        [Ixx, Iyy, Izz, Ixy, Ixz, Iyz]

    """
    xcg, ycg, zcg = xyz_cg
    xref, yref, zref = xyz_ref
    xref2, yref2, zref2 = xyz_ref2

    dx1 = xcg - xref
    dy1 = ycg - yref
    dz1 = zcg - zref

    dx2 = xref2 - xcg
    dy2 = yref2 - ycg
    dz2 = zref2 - zcg

    # consistent with mass_properties, not CONM2
    Ixx_ref, Iyy_ref, Izz_ref, Ixy_ref, Ixz_ref, Iyz_ref = I_ref
    dx = dx1**2 - dx2**2
    dy = dy1**2 - dy2**2
    dz = dz1**2 - dz2**2
    Ixx2 = Ixx_ref - mass * (dy + dz)
    Iyy2 = Iyy_ref - mass * (dx + dz)
    Izz2 = Izz_ref - mass * (dx + dy)
    Ixy2 = Ixy_ref - mass * (dx1 * dy1 - dx2 * dy2)
    Ixz2 = Ixz_ref - mass * (dx1 * dz1 - dx2 * dz2)
    Iyz2 = Iyz_ref - mass * (dy1 * dz1 - dy2 * dz2)
    I_new = np.array([Ixx2, Iyy2, Izz2, Ixy2, Ixz2, Iyz2])
    #print('  Iref = %s' % str(I_ref))
    #print('  Inew = %s' % str(I_new))
    return I_new

def _mass_properties_elements_init(model, element_ids, mass_ids):
    """helper method"""
    # if neither element_id nor mass_ids are specified, use everything
    if isinstance(element_ids, integer_types):
        element_ids = [element_ids]
    if isinstance(mass_ids, integer_types):
        mass_ids = [mass_ids]

    if element_ids is None and mass_ids is None:
        elements = model.elements.values()
        masses = model.masses.values()
        element_ids = [elem.eid for elem in elements]
        mass_ids = [elem.eid for elem in masses]

    # if either element_id or mass_ids are specified and the other is not, use only the
    # specified ids
    #
    # TODO: If eids are requested, but don't exist, no warning is thrown.
    #       Decide if this is the desired behavior.
    else:
        if element_ids is None:
            element_ids = []
            elements = []
        else:
            assert len(model.elements) > 0
            elements = [element for eid, element in model.elements.items() if eid in element_ids]

        if mass_ids is None:
            mass_ids = []
            masses = []
        else:
            assert len(model.masses) > 0
            masses = [mass for eid, mass in model.masses.items() if eid in mass_ids]
    assert element_ids is not None, element_ids
    assert mass_ids is not None, mass_ids
    return element_ids, elements, mass_ids, masses

def mass_properties(model, element_ids=None, mass_ids=None,
                    reference_point=None,
                    sym_axis=None, scale=None, inertia_reference='cg'):
    """
    Calculates mass properties in the global system about the
    reference point.

    Parameters
    ----------
    model : BDF()
        a BDF object
    elements : List[int]; ndarray
        the element ids to consider
    masses : List[int]; ndarray
        the mass ids to consider
    reference_point : (3, ) ndarray; default = <0,0,0>.
        an array that defines the origin of the frame.
    inertia_reference : str; default='cg'
        'cg' : inertia is taken about the cg
        'ref' : inertia is about the reference point

    Returns
    -------
    mass : float
        the mass of the model
    cg : (3, ) float NDARRAY
        the cg of the model as an array.
    I : (6, ) float NDARRAY
        moment of inertia array([Ixx, Iyy, Izz, Ixy, Ixz, Iyz])

    .. seealso:: model.mass_properties

    """
    reference_point, is_cg = _update_reference_point(
        model, reference_point, inertia_reference)
    element_ids, elements, mass_ids, masses = _mass_properties_elements_init(
        model, element_ids, mass_ids)
    mass, cg, I = _mass_properties(
        model, elements, masses,
        reference_point, is_cg)
    mass, cg, I = _apply_mass_symmetry(model, sym_axis, scale, mass, cg, I)
    return mass, cg, I

def _update_reference_point(model, reference_point, inertia_reference='cg'):
    """helper method for handling reference point"""
    inertia_reference = inertia_reference.lower()
    if inertia_reference == 'cg':
        is_cg = True  # nastran-style inertia is always about the cg
    elif inertia_reference == 'ref':
        is_cg = False # inertia is about the reference point
    else:
        raise ValueError("inertia_reference=%r and must be 'cg' or 'ref'" % inertia_reference)

    if reference_point is None:
        reference_point = np.array([0., 0., 0.])
    elif isinstance(reference_point, integer_types):
        reference_point = model.nodes[reference_point].get_position()
    else:
        reference_point = np.asarray(reference_point, dtype='float64')
        if len(reference_point.shape) != 1 or len(reference_point) != 3:
            msg = ("reference_point=%r and must be None, "
                   "a list of 3 floats, or an integer (node id)"  % reference_point)
            raise ValueError(msg)
    return reference_point, is_cg

def mass_properties_no_xref(model, element_ids=None, mass_ids=None,
                            reference_point=None,
                            sym_axis=None, scale=None, inertia_reference='cg',):
    """see model.mass_properties_no_xref"""
    reference_point, is_cg = _update_reference_point(
        model, reference_point, inertia_reference)
    element_ids, elements, mass_ids, masses = _mass_properties_elements_init(
        model, element_ids, mass_ids)
    #nelements = len(elements) + len(masses)

    mass, cg, I = _mass_properties_no_xref(
        model, elements, masses,
        reference_point, is_cg)

    mass, cg, I = _apply_mass_symmetry(model, sym_axis, scale, mass, cg, I)
    return mass, cg, I

def _mass_properties(model, elements, masses, reference_point, is_cg):
    """helper method for ``mass_properties``"""
    mass = 0.
    cg = array([0., 0., 0.])
    inertia = array([0., 0., 0., 0., 0., 0., ])
    no_mass = NO_MASS
    for pack in [elements, masses]:
        for element in pack:
            if element.type == 'CBEAM':
                mass = _get_cbeam_mass_no_nsm(model, element, mass, cg, inertia, reference_point)
                continue

            try:
                p = element.center_of_mass()  # was Centroid()
            except AttributeError:
                if element.type in no_mass:
                    continue
                model.log.error(element.rstrip())
                raise

            try:
                m = element.Mass()
                #print('eid=%s type=%s mass=%s'  %(element.eid, element.type, m))
            except:
                #raise
                if element.type in no_mass:
                    continue
                # PLPLANE
                if element.pid_ref.type == 'PSHELL':
                    model.log.warning('p=%s reference_point=%s type(reference_point)=%s' % (
                        p, reference_point, type(reference_point)))
                    raise
                model.log.warning("could not get the inertia for element/property\n%s%s" % (
                    element, element.pid_ref))
                continue
            mass += m
            cg += m * p
            (x, y, z) = p - reference_point
            x2 = x * x
            y2 = y * y
            z2 = z * z
            inertia[0] += m * (y2 + z2)  # Ixx
            inertia[1] += m * (x2 + z2)  # Iyy
            inertia[2] += m * (x2 + y2)  # Izz
            inertia[3] += m * x * y      # Ixy
            inertia[4] += m * x * z      # Ixz
            inertia[5] += m * y * z      # Iyz

    if mass:
        cg /= mass

    # only transform if we're calculating the inertia about the cg
    if is_cg:
        xyz_ref = reference_point
        xyz_ref2 = cg
        inertia = transform_inertia(mass, cg, xyz_ref, xyz_ref2, inertia)
    return mass, cg, inertia

def _mass_properties_no_xref(model, elements, masses, reference_point, is_cg):  # pragma: no cover
    """
    Calculates mass properties in the global system about the
    reference point.

    Parameters
    ----------
    model : BDF()
        a BDF object
    elements : List[int]; ndarray
        the element ids to consider
    masses : List[int]; ndarray
        the mass ids to consider
    reference_point : (3, ) ndarray; default = <0,0,0>.
        an array that defines the origin of the frame.
    is_cg : bool
        is the reference point the CG

    Returns
    -------
    mass : float
        the mass of the model
    cg : (3, ) float NDARRAY
        the cg of the model as an array.
    I : (6, ) float NDARRAY
        moment of inertia array([Ixx, Iyy, Izz, Ixy, Ixz, Iyz])

    .. seealso:: mass_properties(...)

    """
    mass = 0.
    cg = array([0., 0., 0.])
    I = array([0., 0., 0., 0., 0., 0., ])
    for pack in [elements, masses]:
        for element in pack:
            try:
                p = element.Centroid_no_xref(model)
            except:
                #continue
                raise

            try:
                m = element.Mass_no_xref(model)
            except:
                # PLPLANE
                pid_ref = model.Property(element.pid)
                if pid_ref.type == 'PSHELL':
                    model.log.warning('p=%s reference_point=%s type(reference_point)=%s' % (
                        p, reference_point, type(reference_point)))
                    raise
                model.log.warning("could not get the inertia for element/property\n%s%s" % (
                    element, element.pid_ref))
                continue
            (x, y, z) = p - reference_point
            x2 = x * x
            y2 = y * y
            z2 = z * z
            I[0] += m * (y2 + z2)  # Ixx
            I[1] += m * (x2 + z2)  # Iyy
            I[2] += m * (x2 + y2)  # Izz
            I[3] += m * x * y      # Ixy
            I[4] += m * x * z      # Ixz
            I[5] += m * y * z      # Iyz
            mass += m
            cg += m * p

    if mass:
        cg /= mass

    # only transform if we're calculating the inertia about the cg
    if is_cg:
        xyz_ref = reference_point
        xyz_ref2 = cg
        I = transform_inertia(mass, cg, xyz_ref, xyz_ref2, I)
    return mass, cg, I

def _increment_inertia(centroid, reference_point, m, mass, cg, I):
    """helper method"""
    if m == 0.:
        return mass
    (x, y, z) = centroid - reference_point
    x2 = x * x
    y2 = y * y
    z2 = z * z
    I[0] += m * (y2 + z2)  # Ixx
    I[1] += m * (x2 + z2)  # Iyy
    I[2] += m * (x2 + y2)  # Izz
    I[3] += m * x * y      # Ixy
    I[4] += m * x * z      # Ixz
    I[5] += m * y * z      # Iyz
    mass += m
    cg += m * centroid
    return mass

def mass_properties_nsm(model, element_ids=None, mass_ids=None, nsm_id=None,
                        reference_point=None,
                        sym_axis=None, scale=None, inertia_reference='cg',
                        xyz_cid0_dict=None, debug=False):
    """
    Calculates mass properties in the global system about the
    reference point.  Considers NSM, NSM1, NSML, NSML1.

    Parameters
    ----------
    model : BDF()
        a BDF object
    element_ids : list[int]; (n, ) ndarray, optional
        An array of element ids.
    mass_ids : list[int]; (n, ) ndarray, optional
        An array of mass ids.
    nsm_id : int
        the NSM id to consider
    reference_point : ndarray/int, optional
        type : ndarray
            An array that defines the origin of the frame.
            default = <0,0,0>.
        type : int
            the node id
    sym_axis : str, optional
        The axis to which the model is symmetric.
        If AERO cards are used, this can be left blank.
        allowed_values = 'no', x', 'y', 'z', 'xy', 'yz', 'xz', 'xyz'
    scale : float, optional
        The WTMASS scaling value.
        default=None -> PARAM, WTMASS is used
        float > 0.0
    inertia_reference : str; default='cg'
        'cg' : inertia is taken about the cg
        'ref' : inertia is about the reference point
    xyz_cid0_dict : dict[nid] : xyz; default=None -> auto-calculate
        mapping of the node id to the global position
    debug : bool; default=False
        developer debug; may be removed in the future

    Returns
    -------
    mass : float
        The mass of the model.
    cg : ndarray
        The cg of the model as an array.
    inertia : ndarray
        Moment of inertia array([Ixx, Iyy, Izz, Ixy, Ixz, Iyz]).

    inertia = mass * centroid * centroid

    .. math:: I_{xx} = m (dy^2 + dz^2)

    .. math:: I_{yz} = -m * dy * dz

    where:

    .. math:: dx = x_{element} - x_{ref}

    .. seealso:: http://en.wikipedia.org/wiki/Moment_of_inertia#Moment_of_inertia_tensor

    .. note::
       This doesn't use the mass matrix formulation like Nastran.
       It assumes m*r^2 is the dominant term.
       If you're trying to get the mass of a single element, it
       will be wrong, but for real models will be correct.

    Examples
    --------
    **mass properties of entire structure**

    >>> mass, cg, inertia = model.mass_properties()
    >>> Ixx, Iyy, Izz, Ixy, Ixz, Iyz = inertia


    **mass properties of model based on Property ID**
    >>> pids = list(model.pids.keys())
    >>> pid_eids = model.get_element_ids_dict_with_pids(pids)
    >>> for pid, eids in sorted(pid_eids.items()):
    >>>     mass, cg, I = mass_properties(model, element_ids=eids)

    Warnings
    --------
     - If eids are requested, but don't exist, no warning is thrown.
       Decide if this is the desired behavior.
     - If the NSMx ALL option is used, the mass from all elements
       will be considered, even if not included in the element set

    """
    # TODO: check CG for F:\work\pyNastran\examples\Dropbox\move_tpl\ac11102g.bdf
    reference_point, is_cg = _update_reference_point(
        model, reference_point, inertia_reference)

    if xyz_cid0_dict is None:
        xyz = {}
        for nid, node in model.nodes.items():
            xyz[nid] = node.get_position()
    else:
        xyz = xyz_cid0_dict

    element_ids, unused_elements, mass_ids, unused_masses = _mass_properties_elements_init(
        model, element_ids, mass_ids)

    mass = 0.
    cg = array([0., 0., 0.])
    inertia = array([0., 0., 0., 0., 0., 0., ])

    idtype = model._upcast_int_dtype(dtype='int32')
    all_eids = np.array(list(model.elements.keys()), dtype=idtype)
    all_eids.sort()

    all_mass_ids = np.array(list(model.masses.keys()), dtype=idtype)
    all_mass_ids.sort()

    #element_nsms, property_nsms = _get_nsm_data(model, nsm_id, debug=debug)
    #def _increment_inertia0(centroid, reference_point, m, mass, cg, I):
        #"""helper method"""
        #(x, y, z) = centroid - reference_point
        #mass += m
        #cg += m * centroid
        #return mass

    etypes_skipped = set()
    #eid_areas = defaultdict(list)
    area_eids_pids = defaultdict(list)
    areas = defaultdict(list)
    nsm_centroids_area = defaultdict(list)

    length_eids_pids = defaultdict(list)
    nsm_centroids_length = defaultdict(list)
    lengths = defaultdict(list)

    no_mass = NO_MASS
    for etype, eids in model._type_to_id_map.items():
        #assert isinstance(eids, list), f'etype={etype} eids={eids} type={type(eids)}'
        if etype in no_mass or len(eids) == 0:
            continue
        #assert isinstance(eids, list), 'etype=%r eids=%s'%  (etype, eids)
        mass, cg, inertia = _get_mass_nsm(
            model, element_ids, mass_ids,
            all_eids, all_mass_ids, etypes_skipped,
            etype, eids, xyz,
            length_eids_pids, nsm_centroids_length, lengths,
            area_eids_pids, nsm_centroids_area, areas,
            mass, cg, inertia, reference_point)

    model_eids = np.array(list(model.elements.keys()), dtype=idtype)
    model_pids = np.array(list(model.properties.keys()), dtype=idtype)
    if debug:  # pragma: no cover
        model.log.debug('model_pids = %s' % model_pids)

    mass = _apply_nsm(model, nsm_id,
                      model_eids, model_pids,
                      area_eids_pids, areas, nsm_centroids_area,
                      length_eids_pids, lengths, nsm_centroids_length,
                      mass, cg, inertia, reference_point, debug=debug)
    assert mass is not None
    if mass:
        cg /= mass
    # Ixx, Iyy, Izz, Ixy, Ixz, Iyz = inertia

    # only transform if we're calculating the inertia about the cg
    if is_cg:
        xyz_ref = reference_point
        xyz_ref2 = cg
        inertia = transform_inertia(mass, cg, xyz_ref, xyz_ref2, inertia)

    mass, cg, inertia = _apply_mass_symmetry(model, sym_axis, scale, mass, cg, inertia)
    return mass, cg, inertia


def get_sub_eids(all_eids, eids, etype):
    """supports limiting the element/mass ids"""
    eids = np.array(eids)
    ieids = np.searchsorted(all_eids, eids)
    try:
        eids2 = eids[all_eids[ieids] == eids]
    except IndexError:
        print(etype, all_eids, ieids, eids)
        raise
    return eids2

def _get_mass_nsm(model, element_ids, mass_ids,
                  all_eids, all_mass_ids, etypes_skipped,
                  etype, eids, xyz,
                  length_eids_pids, nsm_centroids_length, lengths,
                  area_eids_pids, nsm_centroids_area, areas,
                  mass, cg, I, reference_point):
    """helper method for ``mass_properties_nsm``"""
    element_ids = set(element_ids)
    mass_ids = set(mass_ids)
    if etype in ['CROD', 'CONROD']:
        eids2 = get_sub_eids(all_eids, eids, etype)
        for eid in eids2:
            elem = model.elements[eid]
            n1, n2 = elem.node_ids
            length = norm(xyz[n2] - xyz[n1])
            centroid = (xyz[n1] + xyz[n2]) / 2.
            mpl = elem.MassPerLength()
            if elem.type == 'CONROD':
                #nsm = property_nsms[nsm_id]['CONROD'][eid] + element_nsms[nsm_id][eid]
                length_eids_pids['CONROD'].append((eid, -42))  # faked number
                lengths['CONROD'].append(length)
                nsm_centroids_length['CONROD'].append(centroid)
            else:
                pid = elem.pid
                #nsm = property_nsms[nsm_id]['PROD'][pid] + element_nsms[nsm_id][eid]
                length_eids_pids['PROD'].append((eid, pid))
                nsm_centroids_length['PROD'].append(centroid)
                lengths['PROD'].append(length)
            #m = (mpl + nsm) * length
            massi = mpl * length
            #if massi != elem.Mass():  # pragma: no cover
                #msg = 'mass_new=%s mass_old=%s\n%s' % (massi, elem.Mass(), str(elem))
                #raise RuntimeError(msg)
            if eid in element_ids:
                mass = _increment_inertia(centroid, reference_point, massi, mass, cg, I)
    elif etype == 'CTUBE':
        eids2 = get_sub_eids(all_eids, eids, etype)
        for eid in eids2:
            elem = model.elements[eid]
            pid = elem.pid
            n1, n2 = elem.node_ids
            length = norm(xyz[n2] - xyz[n1])
            centroid = (xyz[n1] + xyz[n2]) / 2.
            mpl = elem.pid_ref.MassPerLength()
            length_eids_pids['PTUBE'].append((eid, pid))
            lengths['PTUBE'].append(length)
            #nsm = property_nsms[nsm_id]['PTUBE'][pid] + element_nsms[nsm_id][eid]
            #m = (mpl + nsm) * length
            massi = mpl * length
            #if massi != elem.Mass():  # pragma: no cover
                #msg = 'mass_new=%s mass_old=%s\n%s' % (massi, elem.Mass(), str(elem))
                #raise RuntimeError(msg)
            if eid in element_ids:
                mass = _increment_inertia(centroid, reference_point, massi, mass, cg, I)
    elif etype == 'CBAR':
        eids2 = get_sub_eids(all_eids, eids, etype)
        for eid in eids2:
            elem = model.elements[eid]
            pid = elem.pid
            n1, n2 = elem.node_ids
            centroid = (xyz[n1] + xyz[n2]) / 2.
            length = norm(xyz[n2] - xyz[n1])
            mpl = elem.pid_ref.MassPerLength()
            length_eids_pids['PBAR'].append((eid, pid))
            lengths['PBAR'].append(length)
            nsm_centroids_length['PBAR'].append(centroid)
            #nsm = property_nsms[nsm_id]['PBAR'][pid] + element_nsms[nsm_id][eid]
            #m = (mpl + nsm) * length
            massi = mpl * length
            #if massi != elem.Mass() or not np.array_equal(centroid, elem.Centroid()):  # pragma: no cover
                #msg = 'mass_new=%s mass_old=%s\n' % (massi, elem.Mass())
                #msg += 'centroid_new=%s centroid_old=%s\n%s' % (
                    #str(centroid), str(elem.Centroid()), str(elem))
                #raise RuntimeError(msg)
            if eid in element_ids:
                mass = _increment_inertia(centroid, reference_point, massi, mass, cg, I)

    elif etype == 'CBEAM':
        mass = _get_cbeam_mass(
            model, xyz, element_ids, all_eids,
            length_eids_pids, lengths, nsm_centroids_length,
            eids, mass, cg, I, reference_point)

    elif etype in ['CTRIA3', 'CTRIA6', 'CTRIAR']:
        mass = _get_tri_mass(
            model, xyz, element_ids, all_eids,
            area_eids_pids, areas, nsm_centroids_area,
            eids, mass, cg, I, reference_point)

    elif etype in ['CQUAD4', 'CQUAD8', 'CQUADR']:
        mass = _get_quad_mass(
            model, xyz, element_ids, all_eids,
            area_eids_pids, areas, nsm_centroids_area,
            eids, mass, cg, I, reference_point)

    elif etype == 'CQUAD':
        eids2 = get_sub_eids(all_eids, eids, etype)
        for eid in eids2:
            elem = model.elements[eid]
            n1, n2, n3, n4 = elem.node_ids[:4]
            prop = elem.pid_ref
            centroid = (xyz[n1] + xyz[n2] + xyz[n3] + xyz[n4]) / 4.
            area = 0.5 * norm(cross(xyz[n3] - xyz[n1], xyz[n4] - xyz[n2]))

            if prop.type == 'PSHELL':
                t = prop.Thickness()
                mpa = prop.nsm + prop.Rho() * t
            elif prop.type in ['PCOMP', 'PCOMPG']:
                mpa = prop.get_mass_per_area()
            elif prop.type == 'PLPLANE':
                continue
                #raise NotImplementedError(prop.type)
            else:
                raise NotImplementedError(prop.type)
            m = area * mpa
            if m != elem.Mass() or not np.array_equal(centroid, elem.Centroid()):  # pragma: no cover
                msg = 'mass_new=%s mass_old=%s\n' % (m, elem.Mass())
                msg += 'centroid_new=%s centroid_old=%s\n%s' % (
                    str(centroid), str(elem.Centroid()), str(elem))
                raise RuntimeError(msg)
            if eid in element_ids:
                mass = _increment_inertia(centroid, reference_point, m, mass, cg, I)

    elif etype == 'CSHEAR':
        eids2 = get_sub_eids(all_eids, eids, etype)
        for eid in eids2:
            elem = model.elements[eid]
            n1, n2, n3, n4 = elem.node_ids
            prop = elem.pid_ref
            pid = elem.pid
            centroid = (xyz[n1] + xyz[n2] + xyz[n3] + xyz[n4]) / 4.
            area = 0.5 * norm(cross(xyz[n3] - xyz[n1], xyz[n4] - xyz[n2]))
            mpa = prop.MassPerArea()

            area_eids_pids['PSHEAR'].append((eid, pid))
            areas['PSHEAR'].append(area)
            nsm_centroids_area['PSHEAR'].append(centroid)

            #nsm = property_nsms[nsm_id]['PSHEAR'][pid] + element_nsms[nsm_id][eid]
            #m = area * (mpa + nsm)
            m = area * mpa
            #if m != elem.Mass() or not np.array_equal(centroid, elem.Centroid()):  # pragma: no cover
                #msg = 'mass_new=%s mass_old=%s\n' % (m, elem.Mass())
                #msg += 'centroid_new=%s centroid_old=%s\n%s' % (
                    #str(centroid), str(elem.Centroid()), str(elem))
                #raise RuntimeError(msg)
            if eid in element_ids:
                mass = _increment_inertia(centroid, reference_point, m, mass, cg, I)
    elif etype in ['CONM1', 'CONM2', 'CMASS1', 'CMASS2', 'CMASS3', 'CMASS4']:
        eids2 = get_sub_eids(all_mass_ids, eids, etype)
        for eid in eids2:
            elem = model.masses[eid]
            m = elem.Mass()
            centroid = elem.Centroid()
            if eid in mass_ids:
                mass = _increment_inertia(centroid, reference_point, m, mass, cg, I)
    elif etype == 'CTETRA':
        eids2 = get_sub_eids(all_eids, eids, etype)
        for eid in eids2:
            elem = model.elements[eid]
            n1, n2, n3, n4 = elem.node_ids[:4]
            centroid = (xyz[n1] + xyz[n2] + xyz[n3] + xyz[n4]) / 4.
            #V = -dot(n1 - n4, cross(n2 - n4, n3 - n4)) / 6.
            volume = -dot(xyz[n1] - xyz[n4], cross(xyz[n2] - xyz[n4], xyz[n3] - xyz[n4])) / 6.
            m = elem.Rho() * volume
            #if m != elem.Mass() or not np.array_equal(centroid, elem.Centroid()):  # pragma: no cover
                #msg = 'mass_new=%s mass_old=%s\n' % (m, elem.Mass())
                #msg += 'centroid_new=%s centroid_old=%s\n%s' % (
                    #str(centroid), str(elem.Centroid()), str(elem))
                #raise RuntimeError(msg)
            if eid in element_ids:
                mass = _increment_inertia(centroid, reference_point, m, mass, cg, I)

    elif etype == 'CPYRAM':
        eids2 = get_sub_eids(all_eids, eids, etype)
        for eid in eids2:
            elem = model.elements[eid]
            n1, n2, n3, n4, n5 = elem.node_ids[:5]
            centroid1 = (xyz[n1] + xyz[n2] + xyz[n3] + xyz[n4]) / 4.
            area1 = 0.5 * norm(cross(xyz[n3]-xyz[n1], xyz[n4]-xyz[n2]))
            centroid5 = xyz[n5]

            #V = (l * w) * h / 3
            #V = A * h / 3
            centroid = (centroid1 + centroid5) / 2.

            #(n1, n2, n3, n4, n5) = self.get_node_positions()
            #area1, c1 = area_centroid(n1, n2, n3, n4)
            #volume = area1 / 3. * norm(c1 - n5)
            volume = area1 / 3. * norm(centroid1 - centroid5)
            m = elem.Rho() * volume
            if m != elem.Mass() or not np.array_equal(centroid, elem.Centroid()):  # pragma: no cover
                msg = 'mass_new=%s mass_old=%s\n' % (m, elem.Mass())
                msg += 'centroid_new=%s centroid_old=%s\n%s' % (
                    str(centroid), str(elem.Centroid()), str(elem))
                raise RuntimeError(msg)
            #print('*eid=%s type=%s mass=%s rho=%s V=%s' % (
                #elem.eid, 'CPYRAM', m, elem.Rho(), volume))
            if eid in element_ids:
                mass = _increment_inertia(centroid, reference_point, m, mass, cg, I)

    elif etype == 'CPENTA':
        eids2 = get_sub_eids(all_eids, eids, etype)
        for eid in eids2:
            elem = model.elements[eid]
            n1, n2, n3, n4, n5, n6 = elem.node_ids[:6]
            area1 = 0.5 * norm(cross(xyz[n3] - xyz[n1], xyz[n2] - xyz[n1]))
            area2 = 0.5 * norm(cross(xyz[n6] - xyz[n4], xyz[n5] - xyz[n4]))
            centroid1 = (xyz[n1] + xyz[n2] + xyz[n3]) / 3.
            centroid2 = (xyz[n4] + xyz[n5] + xyz[n6]) / 3.
            centroid = (centroid1 + centroid2) / 2.
            volume = (area1 + area2) / 2. * norm(centroid1 - centroid2)
            m = elem.Rho() * volume
            #if m != elem.Mass() or not np.array_equal(centroid, elem.Centroid()):  # pragma: no cover
                #msg = 'mass_new=%s mass_old=%s\n' % (m, elem.Mass())
                #msg += 'centroid_new=%s centroid_old=%s\n%s' % (
                    #str(centroid), str(elem.Centroid()), str(elem))
                #raise RuntimeError(msg)
            #print('*eid=%s type=%s mass=%s rho=%s V=%s' % (
                #elem.eid, 'CPENTA', m, elem.Rho(), volume))
            if eid in element_ids:
                mass = _increment_inertia(centroid, reference_point, m, mass, cg, I)

    elif etype == 'CHEXA':
        eids2 = get_sub_eids(all_eids, eids, etype)
        for eid in eids2:
            elem = model.elements[eid]
            n1, n2, n3, n4, n5, n6, n7, n8 = elem.node_ids[:8]
            #(A1, c1) = area_centroid(n1, n2, n3, n4)
            centroid1 = (xyz[n1] + xyz[n2] + xyz[n3] + xyz[n4]) / 4.
            area1 = 0.5 * norm(cross(xyz[n3] - xyz[n1], xyz[n4] - xyz[n2]))
            #(A2, c2) = area_centroid(n5, n6, n7, n8)
            centroid2 = (xyz[n5] + xyz[n6] + xyz[n7] + xyz[n8]) / 4.
            area2 = 0.5 * norm(cross(xyz[n7] - xyz[n5], xyz[n8] - xyz[n6]))

            volume = (area1 + area2) / 2. * norm(centroid1 - centroid2)
            m = elem.Rho() * volume
            centroid = (centroid1 + centroid2) / 2.
            #if m != elem.Mass() or not np.array_equal(centroid, elem.Centroid()):  # pragma: no cover
                #msg = 'mass_new=%s mass_old=%s\n' % (m, elem.Mass())
                #msg = 'centroid_new=%s centroid_old=%s\n%s' % (
                    #str(centroid), str(elem.Centroid()), str(elem))
                #raise RuntimeError(msg)
            #print('*centroid1=%s centroid2=%s' % (str(centroid1), str(centroid2)))
            #print('*area1=%s area2=%s length=%s' % (area1, area2, norm(centroid1 - centroid2)))
            #print('*eid=%s type=%s mass=%s rho=%s V=%s' % (
                #elem.eid, 'CHEXA', m, elem.Rho(), volume))
            if eid in element_ids:
                mass = _increment_inertia(centroid, reference_point, m, mass, cg, I)

    elif etype == 'CBEND':
        model.log.info('elem.type=%s mass is innaccurate' % etype)
        #nsm = property_nsms[nsm_id]['PBEND'][pid] + element_nsms[nsm_id][eid]
        eids2 = get_sub_eids(all_eids, eids, etype)
        for eid in eids2:
            elem = model.elements[eid]
            m = elem.Mass()
            centroid = elem.Centroid()
            if eid in element_ids:
                mass = _increment_inertia(centroid, reference_point, m, mass, cg, I)

    elif etype in ['CQUADX']:
        pass
    elif etype in ['CTRIAX', 'CTRIAX6']:
        mass = _mass_catch_all(model, etype, etypes_skipped,
                               element_ids, all_eids, eids,
                               mass, cg, I, reference_point)
    elif etype in ['CSUPER', 'CSUPEXT']:
        pass
    elif etype.startswith('C'):
        model.log.warning('etype=%r should be explicit' % etype)
        #raise RuntimeError('etype=%r should be explicit' % etype) ## TODO: this is temporary
        mass = _mass_catch_all(model, etype, etypes_skipped,
                               element_ids, all_eids, eids,
                               mass, cg, I, reference_point)

    #property_nsms[nsm_id][nsm.nsm_type][nsm_idi]
    #for nsm_id, prop_types in sorted(property_nsms.items()):
        #for prop_type, prop_id_to_val in sorted(prop_types.items()):
            #for pid, val in sorted(prop_id_to_val.items()):
        #TODO: CRAC2D mass not supported...how does this work???
        #      I know it's an "area" element similar to a CQUAD4
        #TODO: CCONEAX mass not supported...how does this work???
        #TODO: CBEND mass not supported...how do I calculate the length?

        #area_eids['PSHELL'].append(eid)
        #areas['PSHELL'].append(area)
    return mass, cg, I

def _mass_catch_all(model, etype, etypes_skipped,
                    element_ids, all_eids, eids,
                    mass, cg, I, reference_point):
    """helper method for ``get_mass_new``"""
    eids2 = get_sub_eids(all_eids, eids, etype)
    for eid in eids2:
        elem = model.elements[eid]
        #if elem.pid_ref.type in ['PPLANE']:
        try:
            m = elem.Mass()
        except:
            model.log.error('etype = %r' % etype)
            model.log.error(elem)
            model.log.error(elem.pid_ref)
            raise
        centroid = elem.Centroid()
        if m > 0.0:
            model.log.info('elem.type=%r is not supported in new '
                           'mass properties method' % elem.type)
            if eid in element_ids:
                mass = _increment_inertia(centroid, reference_point, m, mass, cg, I)
        elif etype not in etypes_skipped:
            model.log.info('elem.type=%s doesnt have mass' % elem.type)
            etypes_skipped.add(etype)
    return mass

def _get_cbeam_mass(model, xyz, element_ids, all_eids,
                    length_eids_pids, lengths, nsm_centroids_length,
                    eids, mass, cg, inertia, reference_point):
    """helper method for ``get_mass_new``"""
    eids2 = get_sub_eids(all_eids, eids, 'CBEAM')
    for eid in eids2:
        elem = model.elements[eid]
        prop = elem.pid_ref
        pid = elem.pid
        n1, n2 = elem.node_ids
        xyz1 = xyz[n1]
        xyz2 = xyz[n2]
        centroid = (xyz1 + xyz2) / 2.
        length = norm(xyz2 - xyz1)

        is_failed, out = elem.get_axes(model)
        if is_failed:
            model.log.error(out)
            raise RuntimeError(out)
        wa, wb, _ihat, jhat, khat = out
        p1 = xyz1 + wa
        p2 = xyz2 + wb
        if prop.type == 'PBEAM':
            rho = prop.Rho()

            # we don't call the MassPerLength method so we can put the NSM centroid
            # on a different axis (the PBEAM is weird)
            mass_per_lengths = []
            nsm_per_lengths = []
            for (area, nsm) in zip(prop.A, prop.nsm):
                mass_per_lengths.append(area * rho)
                nsm_per_lengths.append(nsm)
            mass_per_length = integrate_positive_unit_line(prop.xxb, mass_per_lengths)
            nsm_per_length = integrate_positive_unit_line(prop.xxb, nsm_per_lengths)
            nsm_n1 = (p1 + jhat * prop.m1a + khat * prop.m2a)
            nsm_n2 = (p2 + jhat * prop.m1b + khat * prop.m2b)
            nsm_centroid = (nsm_n1 + nsm_n2) / 2.
            #if nsm != 0.:
                #p1_nsm = p1 + prop.ma
                #p2_nsm = p2 + prop.mb
        elif prop.type == 'PBEAML':
            mass_per_lengths = prop.get_mass_per_lengths()
            #mass_per_length = prop.MassPerLength() # includes simplified nsm

            # m1a, m1b, m2a, m2b=0.
            nsm_centroid = (p1 + p2) / 2.

            # mass_per_length already includes nsm
            mass_per_length = integrate_positive_unit_line(prop.xxb, mass_per_lengths)
            nsm_per_length = 0.

            #nsm_centroid = np.zeros(3) # TODO: what is this...
            #nsm = prop.nsm[0] * length # TODO: simplified
        elif prop.type == 'PBCOMP':
            mass_per_length = prop.MassPerLength()
            nsm_per_length = prop.nsm
            nsm_n1 = (p1 + jhat * prop.m1 + khat * prop.m2)
            nsm_n2 = (p2 + jhat * prop.m1 + khat * prop.m2)
            nsm_centroid = (nsm_n1 + nsm_n2) / 2.
        elif prop.type == 'PBMSECT':
            continue
            #mass_per_length = prop.MassPerLength()
            #m = mass_per_length * length
            #nsm = prop.nsm
        else:  # pragma: no cover
            raise NotImplementedError(prop.type)

        #mpl = elem.pid_ref.MassPerLength()
        #m = mpl * length

        length_eids_pids['PBEAM'].append((eid, pid))
        lengths['PBEAM'].append(length)
        nsm_centroids_length['PBEAM'].append(nsm_centroid)
        m = mass_per_length * length
        nsm = nsm_per_length * length
        if (m + nsm) != elem.Mass() or not np.array_equal(centroid, elem.Centroid()):  # pragma: no cover
            msg = 'CBEAM; eid=%s; %s pid=%s; m/L=%s nsm/L=%s; length=%s\n' % (
                eid, pid, prop.type, mass_per_length, nsm_per_length, length)
            msg += 'mass_new=%s mass_old=%s\n' % (m, elem.Mass())
            msg += 'centroid_new=%s centroid_old=%s\n%s' % (
                str(centroid), str(elem.Centroid()), str(elem))
            raise RuntimeError(msg)

        if eid not in element_ids:
            continue
        #nsm = (nsm_per_length + nsmi) * length
        (x, y, z) = centroid - reference_point
        (xm, ym, zm) = nsm_centroid - reference_point
        x2 = x * x
        y2 = y * y
        z2 = z * z
        xm2 = xm * xm
        ym2 = ym * ym
        zm2 = zm * zm

        # Ixx, Iyy, Izz, Ixy, Ixz, Iyz
        inertia[0] += m * (y2 + z2) + nsm * (ym2 + zm2)
        inertia[1] += m * (x2 + z2) + nsm * (xm2 + zm2)
        inertia[2] += m * (x2 + y2) + nsm * (xm2 + ym2)
        inertia[3] += m * x * y + nsm * xm * ym
        inertia[4] += m * x * z + nsm * xm * zm
        inertia[5] += m * y * z + nsm * ym * zm
        massi = m + nsm
        mass += massi
        cg += m * centroid + nsm * nsm_centroid
        #print('length=%s mass=%s mass_per_length=%s nsm_per_length=%s m=%s nsm=%s centroid=%s nsm_centroid=%s' % (
            #length, mass, mass_per_length, nsm_per_length, m, nsm, centroid, nsm_centroid))
        if massi != elem.Mass():  # pragma: no cover
            msg = 'mass_new=%s mass_old=%s\n' % (massi, elem.Mass())
            msg += 'centroid_new=%s centroid_old=%s\n%s' % (
                str(centroid), str(elem.Centroid()), str(elem))
            raise RuntimeError(msg)
    return mass

def _get_cbeam_mass_no_nsm(model, elem, mass, cg, inertia, reference_point):
    """helper method for mass_properties"""
    prop = elem.pid_ref
    xyz1, xyz2 = elem.get_node_positions()
    centroid = (xyz1 + xyz2) / 2.
    length = norm(xyz2 - xyz1)

    is_failed, out = elem.get_axes(model)
    if is_failed:
        model.log.error(str(out))
        raise RuntimeError(out)

    wa, wb, _ihat, jhat, khat = out
    p1 = xyz1 + wa
    p2 = xyz2 + wb
    if prop.type == 'PBEAM':
        rho = prop.Rho()
        # we don't call the MassPerLength method so we can put the NSM centroid
        # on a different axis (the PBEAM is weird)
        mass_per_lengths = []
        nsm_per_lengths = []
        for (area, nsm) in zip(prop.A, prop.nsm):
            mass_per_lengths.append(area * rho)
            nsm_per_lengths.append(nsm)
        mass_per_length = integrate_positive_unit_line(prop.xxb, mass_per_lengths)
        nsm_per_length = integrate_positive_unit_line(prop.xxb, nsm_per_lengths)
        #print('nsm/Ls=%s nsm/L=%s' % (nsm_per_lengths, nsm_per_length))
        #print('mass/Ls=%s mass/L=%s' % (mass_per_lengths, mass_per_length))
        nsm_n1 = (p1 + jhat * prop.m1a + khat * prop.m2a)
        nsm_n2 = (p2 + jhat * prop.m1b + khat * prop.m2b)
        nsm_centroid = (nsm_n1 + nsm_n2) / 2.

    elif prop.type == 'PBEAML':
        mass_per_lengths = prop.get_mass_per_lengths()
        #mass_per_length = prop.MassPerLength() # includes simplified nsm

        # m1a, m1b, m2a, m2b=0.
        nsm_centroid = (p1 + p2) / 2.

        # mass_per_length already includes nsm
        mass_per_length = integrate_positive_unit_line(prop.xxb, mass_per_lengths)
        nsm_per_length = 0.

        #nsm_centroid = np.zeros(3) # TODO: what is this...
        #nsm = prop.nsm[0] * length # TODO: simplified
    elif prop.type == 'PBCOMP':
        mass_per_length = prop.MassPerLength()
        nsm_per_length = prop.nsm
        nsm_n1 = (p1 + jhat * prop.m1 + khat * prop.m2)
        nsm_n2 = (p2 + jhat * prop.m1 + khat * prop.m2)
        nsm_centroid = (nsm_n1 + nsm_n2) / 2.
    elif prop.type == 'PBMSECT':
        return mass
        #mass_per_length = prop.MassPerLength()
        #m = mass_per_length * length
        #nsm = prop.nsm
    else:  # pragma: no cover
        raise NotImplementedError(prop.type)

    m = mass_per_length * length
    nsm = nsm_per_length * length
    if (m + nsm) != elem.Mass() or not np.array_equal(centroid, elem.Centroid()):  # pragma: no cover
        msg = 'CBEAM; eid=%s; %s pid=%s; m/L=%s nsm/L=%s; length=%s\n' % (
            elem.eid, elem.pid, prop.type, mass_per_length, nsm_per_length, length)
        msg += 'mass_new=%s mass_old=%s\n' % (m, elem.Mass())
        msg += 'centroid_new=%s centroid_old=%s\n%s' % (
            str(centroid), str(elem.Centroid()), str(elem))
        raise RuntimeError(msg)

    #nsm = (nsm_per_length + nsmi) * length
    (x, y, z) = centroid - reference_point
    (xm, ym, zm) = nsm_centroid - reference_point
    x2 = x * x
    y2 = y * y
    z2 = z * z
    xm2 = xm * xm
    ym2 = ym * ym
    zm2 = zm * zm

    # Ixx, Iyy, Izz, Ixy, Ixz, Iyz
    inertia[0] += m * (y2 + z2) + nsm * (ym2 + zm2)
    inertia[1] += m * (x2 + z2) + nsm * (xm2 + zm2)
    inertia[2] += m * (x2 + y2) + nsm * (xm2 + ym2)
    inertia[3] += m * x * y + nsm * xm * ym
    inertia[4] += m * x * z + nsm * xm * zm
    inertia[5] += m * y * z + nsm * ym * zm
    massi = m + nsm
    mass += massi
    cg += m * centroid + nsm * nsm_centroid
    return mass

def _get_tri_mass(model, xyz, element_ids, all_eids,
                  area_eids_pids, areas, nsm_centroids_area,
                  eids, mass, cg, inertia, reference_point):
    """helper method for ``get_mass_new``"""
    eids2 = get_sub_eids(all_eids, eids, 'tri')
    for eid in eids2:
        elem = model.elements[eid]
        n1, n2, n3 = elem.node_ids[:3]
        prop = elem.pid_ref
        pid = elem.pid
        centroid = (xyz[n1] + xyz[n2] + xyz[n3]) / 3.
        area = 0.5 * norm(cross(xyz[n1] - xyz[n2], xyz[n1] - xyz[n3]))
        #areas_prop[pid] += area
        if prop.type == 'PSHELL':
            tflag = elem.tflag
            ti = prop.Thickness()
            if tflag == 0:
                # absolute
                t1 = ti if elem.T1 is None else elem.T1
                t2 = ti if elem.T2 is None else elem.T2
                t3 = ti if elem.T3 is None else elem.T3
            elif tflag == 1:
                # relative
                t1 = ti if elem.T1 is None else elem.T1 * ti
                t2 = ti if elem.T2 is None else elem.T2 * ti
                t3 = ti if elem.T3 is None else elem.T3 * ti
            else:  # pragma: no cover
                raise RuntimeError('tflag=%r' % tflag)
            assert t1 + t2 + t3 > 0., 't1=%s t2=%s t3=%s' % (t1, t2, t3)
            t = (t1 + t2 + t3) / 3.

            # m/A = rho * A * t + nsm
            #mass_per_area = elem.nsm + rho * elem.t

            #areas_prop[pid] += area
            mpa = prop.nsm + prop.Rho() * t
            #mpa = elem.pid_ref.MassPerArea()
        elif prop.type in ['PCOMP', 'PCOMPG']:
            # TODO: do PCOMPs even support differential thickness?
            #       I don't think so...

            #rho_t = prop.get_rho_t()
            #nsm = prop.nsm
            #rho_t = [mat.Rho() * t for (mat, t) in zip(prop.mids_ref, prop.ts)]
            #mpa = sum(rho_t) + nsm

            # works for PCOMP
            # F:\Program Files\Siemens\NXNastran\nxn10p1\nxn10p1\nast\tpl\cqr3compbuck.dat
            mpa = prop.get_mass_per_area()
        elif prop.type in ['PLPLANE', 'PPLANE']:
            continue
        else:
            raise NotImplementedError(prop.type)
        area_eids_pids['PSHELL'].append((eid, pid))
        areas['PSHELL'].append(area)

        nsm_centroids_area['PSHELL'].append(centroid)
        #nsm = property_nsms[nsm_id]['PSHELL'][pid] + element_nsms[nsm_id][eid]
        #m = area * (mpa + nsm)
        massi = area * mpa
        if not np.array_equal(centroid, elem.Centroid()):  # pragma: no cover
            msg = 'centroid_new=%s centroid_old=%s\n%s' % (
                str(centroid), str(elem.Centroid()), str(elem))
            raise RuntimeError(msg)
        if eid in element_ids:
            mass = _increment_inertia(centroid, reference_point, massi, mass, cg, inertia)
    return mass

def _get_quad_mass(model, xyz, element_ids, all_eids,
                   area_eids_pids, areas, nsm_centroids_area,
                   eids, mass, cg, inertia, reference_point):
    """helper method for ``get_mass_new``"""
    eids2 = get_sub_eids(all_eids, eids, 'quad')
    for eid in eids2:
        elem = model.elements[eid]
        n1, n2, n3, n4 = elem.node_ids[:4]
        prop = elem.pid_ref
        pid = prop.pid
        centroid = (xyz[n1] + xyz[n2] + xyz[n3] + xyz[n4]) / 4.
        area = 0.5 * norm(cross(xyz[n3] - xyz[n1], xyz[n4] - xyz[n2]))

        if prop.type == 'PSHELL':
            tflag = elem.tflag
            ti = prop.Thickness()
            if tflag == 0:
                # absolute
                t1 = ti if elem.T1 is None else elem.T1
                t2 = ti if elem.T2 is None else elem.T2
                t3 = ti if elem.T3 is None else elem.T3
                t4 = ti if elem.T4 is None else elem.T4
            elif tflag == 1:
                # relative
                t1 = ti if elem.T1 is None else elem.T1 * ti
                t2 = ti if elem.T2 is None else elem.T2 * ti
                t3 = ti if elem.T3 is None else elem.T3 * ti
                t4 = ti if elem.T4 is None else elem.T4 * ti
            else:  # pragma: no cover
                raise RuntimeError('tflag=%r' % tflag)
            assert t1 + t2 + t3 + t4 > 0., 't1=%s t2=%s t3=%s t4=%s' % (t1, t2, t3, t4)
            t = (t1 + t2 + t3 + t4) / 4.

            # m/A = rho * A * t + nsm
            #mass_per_area = model.nsm + rho * model.t

            #areas_prop[pid] += area
            mpa = prop.nsm + prop.Rho() * t
            #mpa = elem.pid_ref.MassPerArea()
            #m = mpa * area
        elif prop.type in ['PCOMP', 'PCOMPG']:
            # PCOMP, PCOMPG
            #rho_t = prop.get_rho_t()
            #nsm = prop.nsm
            #rho_t = [mat.Rho() * t for (mat, t) in zip(prop.mids_ref, prop.ts)]
            #mpa = sum(rho_t) + nsm
            mpa = prop.get_mass_per_area()
        elif prop.type in ['PLPLANE', 'PPLANE']:
            continue
            #raise NotImplementedError(prop.type)
        else:
            raise NotImplementedError(prop.type)
        area_eids_pids['PSHELL'].append((eid, pid))
        areas['PSHELL'].append(area)
        nsm_centroids_area['PSHELL'].append(centroid)
        #nsm = property_nsms[nsm_id]['PSHELL'][pid] + element_nsms[nsm_id][eid]
        #m = area * (mpa + nsm)
        m = area * mpa
        if not np.array_equal(centroid, elem.Centroid()):  # pragma: no cover
            msg = 'centroid_new=%s centroid_old=%s\n%s' % (
                str(centroid), str(elem.Centroid()), str(elem))
            raise RuntimeError(msg)
        #elem = model.elements[eid]

        #mass_expected = elem.Mass()
        #if not np.allclose(m, mass_expected):  # pragma: no cover
            #msg = 'massi=%s expected=%s' % (m, mass_expected)
            #for node in elem.nodes_ref:
                #node.comment = ''
            #elem.comment = ''
            #prop.comment = ''
            #prop.mid1_ref.comment = ''
            #msg += elem.get_stats()
            #msg += prop.get_stats()
            #msg += prop.mid_ref.get_stats()
            #raise RuntimeError(msg)
        #print('eid=%s type=%s mass=%s; area=%s mpa=%s'  % (elem.eid, elem.type, m, area, mpa))
        if eid in element_ids:
            mass = _increment_inertia(centroid, reference_point, m, mass, cg, inertia)
    return mass

def _setup_apply_nsm(area_eids_pids, areas, nsm_centroids_area,
                     length_eids_pids, lengths, nsm_centroids_length):
    """
    Sets up the non-structural mass processing

    Parameters
    ----------
    area_eids_pids : ???
        ???
    areas : ???
        ???
    nsm_centroids_area : ???
        ???
    length_eids_pids : ???
        ???
    lengths : ???
        ???
    nsm_centroids_length : ???
        ???

    Returns
    -------
    all_eids_pids : ???
        ???
    area_length : ???
        ???
    is_area : (nelements, ) bool ndarray
        is this an area element
    nsm_centroids : (nelements, 3 ) float ndarray
        the centroids for the elements

    """
    nsm_centroids = []
    all_eids_pids = []
    area_length = []
    is_area = []
    #is_data = False
    #print(areas)
    for ptype, eids_pids in area_eids_pids.items():
        areasi = np.array(areas[ptype], dtype='float64')
        area_eids_pids[ptype] = np.array(eids_pids, dtype='int32')
        areas[ptype] = areasi
        assert len(areasi) > 0, areas
        all_eids_pids += eids_pids
        nsm_centroidsi = np.array(nsm_centroids_area[ptype])
        nsm_centroids.append(nsm_centroidsi)
        assert len(eids_pids) == len(nsm_centroids_area[ptype]), ptype
        area_length.append(areasi)
        is_area += [True] * len(areasi)
        #is_data = True
        nsm_centroids_area[ptype] = nsm_centroidsi

    for ptype, eids_pids in length_eids_pids.items():
        lengthsi = np.array(lengths[ptype], dtype='float64')
        length_eids_pids[ptype] = np.array(eids_pids, dtype='int32')
        lengths[ptype] = lengthsi
        assert len(lengthsi) > 0, lengthsi
        all_eids_pids += eids_pids
        nsm_centroidsi = np.array(nsm_centroids_length[ptype])
        nsm_centroids.append(nsm_centroidsi)
        assert len(eids_pids) == len(nsm_centroids_length[ptype]), ptype
        area_length.append(lengthsi)
        is_area += [False] * len(lengthsi)
        #is_data = True
        nsm_centroids_length[ptype] = nsm_centroidsi

    all_eids_pids = np.array(all_eids_pids, dtype='int32')
    nelements = len(is_area)
    if nelements == 0:
        return all_eids_pids, area_length, is_area, nsm_centroids

    isort = np.argsort(all_eids_pids[:, 0])
    all_eids_pids = all_eids_pids[isort, :]
    area_length = np.hstack(area_length)[isort]

    is_area_array = np.array(is_area, dtype='bool')[isort]
    nsm_centroids = np.vstack(nsm_centroids)[isort]
    return all_eids_pids, area_length, is_area_array, nsm_centroids

def _combine_prop_weighted_area_length_simple(model, eids, area, centroids,
                                              nsm_value, reference_point, mass, cg, I,
                                              is_area, divide_by_sum,
                                              debug=True):
    """
    Calculates the contribution of NSM/NSML/NSM1 cards on mass properties.
    Area/Length are abstracted, so if you have a shell element, the area_sum
    is an area_sum, whereas if you have a line element, it's a lenght_sum.

    The standard (non-NSM card) shell mass is calculated as:
        mass = A * (nsm_per_unit_area + t * rho)

    For an NSM/NSM1 card, we simply modify nsm_per_unit_area either on a
    per element basis or property basis.

    The following NSM/NSML cards area identical (for any number of elements) ::

        NSML, 42, ELEMENT, 100, 0.1
        NSM,  42, ELEMENT, 100, 0.1

    However, the following NSM/NSML property cards are different:
        NSML, 42, PSHELL, 100, 0.1
        NSM,  42, PSHELL, 100, 0.1

    The NSM card defines nsm_unit_area,  while NSML defines
    total nsm  or ``A*nsm_per_unit_area``.

    """
    if debug:  # pragma: no cover
        model.log.debug('_combine_weighted_area_length_simple')
    assert nsm_value is not None
    if len(area) != len(centroids):
        msg = 'len(area)=%s len(centroids)=%s ' % (len(area), len(centroids))
        raise RuntimeError(msg)

    if is_area:
        word = 'area'
    else:
        word = 'length'
    if divide_by_sum:
        area_sum = area.sum()
        area = area / area_sum
        if debug:
            model.log.debug('dividing by %s=%s' % (word, area_sum))

    for eid, areai, centroid in zip(eids, area, centroids):
        #area_sum += areai
        m = nsm_value * areai
        if debug:  # pragma: no cover
            model.log.debug('  eid=%s %si=%s nsm_value=%s mass=%s %s=%s' % (
                eid, word, areai, nsm_value, m, word, areai))
        #elem = model.elements[eid]
        #assert np.allclose(m, elem.Mass()), elem.get_stats()
        mass = _increment_inertia(centroid, reference_point, m, mass, cg, I)
    if debug:  # pragma: no cover
        model.log.debug('mass = %s' % mass)
    return mass

def _combine_prop_weighted_area_length(
        model, areas_ipids, nsm_centroidsi, is_area, area_sum,
        nsm_value, reference_point, mass, cg, I, debug=True):
    """
    Calculates the contribution of NSML cards on mass properties.
    Area/Length are abstracted, so if you have a shell element, the area_sum
    is an area_sum, whereas if you have a line element, it's a lenght_sum.

    The NSML mass contribution is calculated as a distrubted mass:
        mass = nsm_per_unit_area * areai / area_sum

    """
    assert area_sum is not None
    assert nsm_value is not None
    if is_area:
        word = 'area'
    else:
        word = 'length'

    for (area, ipid) in areas_ipids:
        if debug:  # pragma: no cover
            model.log.debug("  nsm_centroidsi = %s" % nsm_centroidsi)
        centroids = nsm_centroidsi[ipid, :]

        area2 = area / area_sum
        for areai, centroid in zip(area2, centroids):
            #print('  areai=%s area_sum=%s nsm_value=%s' % (areai*area_sum, area_sum, nsm_value))
            m = nsm_value * areai
            if debug:  # pragma: no cover
                model.log.debug('  %si=%s %s_sum=%s nsm_value=%s mass=%s' % (
                    word, areai*area_sum, word, area_sum, nsm_value, m))
            #assert np.allclose(m, elem.Mass()), elem.get_stats()
            mass = _increment_inertia(centroid, reference_point, m, mass, cg, I)
    return mass

def _apply_nsm(model, nsm_id,
               unused_model_eids, unused_model_pids,
               area_eids_pids, areas, nsm_centroids_area,
               length_eids_pids, lengths, nsm_centroids_length,
               mass, cg, I, reference_point, debug=False):
    """
    Applies NSM cards to the mass, cg, and inertia.

    Parameters
    ----------
    model : BDF()
        a BDF object
    nsm_id : int
        the NSM id to consider
    reference_point : ndarray/int, optional
        type : ndarray
            An array that defines the origin of the frame.
            default = <0,0,0>.
        type : int
            the node id
    area_eids_pids : Dict[etype_ptype] = eids_pids
        etype_ptype : str
            the element or property type (e.g., CQUAD4, PSHELL)
        eids_pids : (neids, 2) ndarray
            (element_id, property_id) for area elements
    areas : Dict[etype_ptype] = area
        area : (nelements, ) float ndarray
            area for area elements
    nsm_centroids_area : Dict[etype_ptype] = centroids
        centroids : (nelements, 3) float ndarray
            centroids for area elements
    length_eids_pids : Dict[etype_ptype] = eids_pids
        etype_ptype : str
            the element or property type (e.g., CQUAD4, PSHELL)
        eids_pids : (neids, 2) ndarray
            (element_id, property_id) for area elements
    lengths :  Dict[etype_ptype] = length
        length : (nelements, ) float ndarray
            lengths for length elements (e.g., CBAR, PBEAM, CONROD)
    nsm_centroids_length : Dict[etype_ptype] = centroids
        centroids : (nelements, 3) float ndarray
            centroids for length elements
    mass : float
        The mass of the model
    cg : ndarray
        The cg of the model as an array
    I : ndarray
        Moment of inertia array([Ixx, Iyy, Izz, Ixy, Ixz, Iyz])

    Returns
    -------
    mass : float
        The mass of the model

    per MSC QRG 2018.0.1: Undefined property/element IDs are ignored.
    TODO: support ALL

    """
    if not nsm_id:
        return mass

    #print(length_eids_pids)
    #print(lengths)
    nsms = model.get_reduced_nsms(nsm_id, consider_nsmadd=True, stop_on_failure=True)
    if debug:  # pragma: no cover
        for nsm in nsms:
            model.log.debug(nsm.rstrip())

    nsm_type_map = {
        'PSHELL' : 'PSHELL',
        'PCOMP' : 'PSHELL',
        'PCOMPG' : 'PSHELL',

        'PBAR' : 'PBAR',
        'PBARL' : 'PBAR',

        'PBEAM' : 'PBEAM',
        'PBEAML' : 'PBEAM',
        'PBCOMP' : 'PBEAM',

        'PROD' : 'PROD',
        'PBEND' : 'PBEND',
        'PSHEAR' : 'PSHEAR',
        'PTUBE' : 'PTUBE',
        'PCONEAX' : 'PCONEAX',
        'PRAC2D' : 'PRAC2D',
        'CONROD' : 'CONROD',
        'ELEMENT' : 'ELEMENT',
    }
    #all_eid_nsms = []
    if len(nsms) == 0:
        model.log.warning('no nsm...')
        return mass

    all_eids_pids, area_length, is_area_array, nsm_centroids = _setup_apply_nsm(
        area_eids_pids, areas, nsm_centroids_area,
        length_eids_pids, lengths, nsm_centroids_length)

    nelements = len(is_area_array)
    if nelements == 0:
        model.log.debug('  skipping NSM=%s calc because there are no elements\n' % nsm_id)
        return mass

    #print('all_eids_pids =', all_eids_pids)
    #print('area_length =', area_length)
    #print('is_area_array =', is_area_array)
    neids = all_eids_pids.shape[0]
    assert neids == len(area_length), 'len(eids)=%s len(area_length)=%s' % (neids, len(area_length))
    assert neids == len(is_area_array), 'len(eids)=%s len(area_length)=%s' % (neids, len(is_area_array))
    #area_length = area_length[isort]
    #is_area = is_area[isort]

    assert isinstance(is_area_array, np.ndarray), type(is_area_array)
    for nsm in nsms:
        nsm_value = nsm.value
        nsm_type = nsm_type_map[nsm.nsm_type]
        if debug:  # pragma: no cover
            model.log.debug('-' * 80)
            model.log.debug(nsm)
            model.log.debug("nsm_type=%r value=%s" % (nsm_type, nsm_value))

        divide_by_sum = False
        if nsm.type in ['NSML1', 'NSML']:
            divide_by_sum = True

        if nsm.type == 'NSML1':
            if nsm_type == 'PSHELL': # area
                mass = _get_nsml1_prop(
                    model, nsm, nsm_type, nsm_value,
                    area_eids_pids, areas, nsm_centroids_area,
                    mass, cg, I, reference_point, is_area=True, debug=debug)
            elif nsm_type in ['PBAR', 'PBEAM', 'PROD', 'PTUBE']:
                mass = _get_nsml1_prop(
                    model, nsm, nsm_type, nsm_value,
                    length_eids_pids, lengths, nsm_centroids_length,
                    mass, cg, I, reference_point, is_area=False, debug=debug)
            elif nsm_type in ['ELEMENT', 'CONROD']:
                if len(nsm.ids) == 1 and nsm.ids[0] == 'ALL':
                    if nsm_type == 'CONROD':
                        nsm_ids = model._type_to_id_map[nsm_type]
                    else:
                        nsm_ids = list(model.elements.keys())
                else:
                    nsm_ids = nsm.ids
                #eids_pids = length_eids_pids[nsm_type]
                #print('nsm_type=%s eids_pids=%s' % (nsm_type, eids_pids))
                #if len(eids_pids) == 0:
                    #model.log.warning('  *skipping because there are no elements '
                                      #'associated with:\n%s' % str(nsm))
                    #continue

                mass = _nsm1_element(
                    model, nsm, nsm_ids,
                    all_eids_pids, area_length, nsm_centroids,
                    mass, cg, I, reference_point, is_area_array,
                    divide_by_sum, debug=debug)
            else:
                raise NotImplementedError(nsm_type)

        elif nsm.type in ['NSM1', 'NSML', 'NSM']:
            if nsm_type == 'PSHELL': # area
                pids = nsm.ids
                if debug:  # pragma: no cover
                    model.log.debug(nsm.rstrip())
                eids_pids = area_eids_pids[nsm_type]
                if len(eids_pids) == 0:
                    model.log.warning('  *skipping because there are no elements '
                                      'associated with:\n%s' % str(nsm))
                    continue

                area_all = areas[nsm_type]
                all_eids = eids_pids[:, 0]
                all_pids = eids_pids[:, 1]

                is_area = True
                if len(pids) == 1 and pids[0] == 'ALL':
                    #model.log.warning('  *skipping %s/PSHELL/ALL\n%s' % (nsm.type, str(nsm)))
                    centroidsi = nsm_centroids_area[nsm_type]
                    mass = _combine_prop_weighted_area_length_simple(
                        model, all_eids, area_all, centroidsi,
                        nsm_value, reference_point, mass, cg, I,
                        is_area, divide_by_sum,
                        debug=debug)
                else:
                    for pid in pids:
                        assert isinstance(pid, int), 'pid=%s type=%s' % (pid, type(pid))
                        ieidsi = np.where(all_pids == pid)
                        eidsi = all_eids[ieidsi]
                        centroidsi = nsm_centroids_area[nsm_type][ieidsi]
                        areasi = area_all[ieidsi]
                        if len(centroidsi) != len(eidsi):
                            msg = 'ncentroids=%s neids=%s' % (len(centroidsi), len(eidsi))
                            raise RuntimeError(msg)

                        if debug:  # pragma: no cover
                            #print('eids = %s' % all_eids)
                            model.log.debug('  eidsi = %s' % eidsi)
                            model.log.debug('  nsm_centroids_area = %s' % centroidsi)
                            model.log.debug('  centroidsi = %s' % centroidsi)

                        mass = _combine_prop_weighted_area_length_simple(
                            model, eidsi, areasi, centroidsi,
                            nsm_value, reference_point, mass, cg, I,
                            is_area, divide_by_sum,
                            debug=debug)
            elif nsm_type in ['PBAR', 'PBEAM', 'PROD', 'PTUBE']:
                length_all = np.array(lengths[nsm_type])
                pids = nsm.ids
                eids_pids = length_eids_pids[nsm_type]
                if len(eids_pids) == 0:
                    model.log.debug('  *skipping because there are no elements'
                                    ' associated with:\n%s' % str(nsm))
                    continue

                length_all = np.array(lengths[nsm_type])
                all_eids = eids_pids[:, 0]
                all_pids = eids_pids[:, 1]
                is_area = False

                nsm_centroidsi = nsm_centroids_length[nsm_type]
                if len(pids) == 1 and pids[0] == 'ALL':
                    lengthsi = length_all
                    centroidsi = nsm_centroidsi
                    mass = _combine_prop_weighted_area_length_simple(
                        model, all_eids, lengthsi, centroidsi,
                        nsm_value, reference_point, mass, cg, I,
                        is_area, divide_by_sum,
                        debug=debug)
                else:
                    for pid in pids:
                        assert isinstance(pid, int), 'pid=%s type=%s' % (pid, type(pid))
                        ieidsi = np.where(all_pids == pid)
                        eidsi = all_eids[ieidsi]
                        centroidsi = nsm_centroidsi[ieidsi]
                        lengthsi = length_all[ieidsi]
                        if len(centroidsi) != len(eidsi):
                            msg = 'ncentroids=%s neids=%s' % (len(centroidsi), len(eidsi))
                            raise RuntimeError(msg)

                        if debug:  # pragma: no cover
                            model.log.debug('  eidsi = %s' % eidsi)
                            model.log.debug('  nsm_centroids_lengthi = %s' % centroidsi)
                            model.log.debug('  centroidsi = %s' % centroidsi)

                        mass = _combine_prop_weighted_area_length_simple(
                            model, eidsi, lengthsi, centroidsi,
                            nsm_value, reference_point, mass, cg, I,
                            is_area, divide_by_sum,
                            debug=debug)
            elif nsm_type in ['ELEMENT', 'CONROD']:
                if len(nsm.ids) == 1 and nsm.ids[0] == 'ALL':
                    if nsm_type == 'CONROD':
                        nsm_ids = model._type_to_id_map[nsm_type]
                    else:
                        nsm_ids = list(model.elements.keys())
                else:
                    nsm_ids = nsm.ids

                mass = _nsm1_element(
                    model, nsm, nsm_ids,
                    all_eids_pids, area_length, nsm_centroids,
                    mass, cg, I, reference_point, is_area_array,
                    divide_by_sum, debug=debug)
            else:
                raise NotImplementedError(nsm_type)
        else:
            model.log.warning('skipping %s\n%s' % (nsm.type, str(nsm)))


    #print('area:')
    #for ptype, eids_pids in sorted(area_eids_pids.items()):
        #eids = eids_pids[:, 0]
        #pids = eids_pids[:, 1]
        #area = np.array(areas[ptype])
        #ieids = np.argsort(eids)
        #eids_sorted = eids[ieids]
        #area_sorted = area[ieids]
        #print('  ', ptype, eids_sorted, area_sorted)

    #print('length:')
    #for ptype, length_eid in sorted(length_eids_pids.items()):
        #eids = np.array(length_eid, dtype='int32')
        #length = np.array(lengths[ptype])
        #print('  ', ptype, eids, length)
    return mass

def _get_nsml1_prop(
        model, nsm, nsm_type, nsm_value,
        area_eids_pids, areas, nsm_centroids_area,
        mass, cg, I, reference_point, is_area=True, debug=True):
    if is_area:
        word = 'area'
        key_map = {
            'PSHELL' : ['PSHELL', 'PCOMP', 'PCOMPG'],
            'PCOMP' : ['PSHELL', 'PCOMP', 'PCOMPG'],
            'PCOMPG' : ['PSHELL', 'PCOMP', 'PCOMPG'],
        }
    else:
        word = 'length'
        key_map = {
            'PBAR' : ['PBAR', 'PBARL'],
            'PBEAM' : ['PBEAM', 'PBEAML', 'PBCOMP'],
            'PROD' : ['PROD'],
        }
    keys = key_map[nsm_type]

    eids_pids = area_eids_pids[nsm_type]
    area_all = areas[nsm_type]
    if len(eids_pids) == 0:
        model.log.warning('  skipping because there are no elements'
                          ' associated with:\n%s' % str(nsm))
        return mass

    #all_eids = eids_pids[:, 0]
    all_pids = eids_pids[:, 1]
    if len(nsm.ids) == 1 and nsm.ids[0] == 'ALL':
        ids = np.hstack([model._type_to_id_map[key]
                         for key in keys if key in model._type_to_id_map])
    else:
        assert 'ALL' not in nsm.ids, str(nsm)
        ids = np.array(nsm.ids, dtype='int32')

    pids_to_apply = all_pids
    upids = np.unique(pids_to_apply)
    pids_to_apply = np.intersect1d(upids, ids)

    if debug:  # pragma: no cover
        model.log.debug("  all_pids = %s" % all_pids)
        model.log.debug("  nsm_pids = %s" % ids)
        model.log.debug("  pids_to_apply = %s" % pids_to_apply)
    assert len(pids_to_apply) > 0, pids_to_apply

    area_sum = 0.
    areas_ipids = []
    for upid in pids_to_apply:
        ipid = np.where(all_pids == upid)[0]
        eids = eids_pids[ipid, 0]
        area = area_all[ipid]

        #eids_actual = eids[ipid]
        #area_actual = area[ipid]
        if debug:  # pragma: no cover
            model.log.debug('  eids = %s' % eids)
            model.log.debug('  %s = %s' % (word, area))
        area_sum += area.sum()
        areas_ipids.append((area, ipid))

    if debug:  # pragma: no cover
        model.log.debug("%s_sum = %s" % (word, area_sum))
    nsm_centroidsi = nsm_centroids_area[nsm_type]

    mass = _combine_prop_weighted_area_length(
        model, areas_ipids, nsm_centroidsi,
        is_area, area_sum,
        nsm_value, reference_point, mass, cg, I,
        debug=debug)
    return mass

def _nsm1_element(model, nsm, nsm_ids, all_eids_pids, area_length, nsm_centroids,
                  mass, cg, I, reference_point, is_area_array,
                  divide_by_sum, debug=False):
    """calculates the mass of an NSM1 element"""
    nsm_value = nsm.value
    #model.log.warning('  *skipping NSM1/ELEMENT\n%s' % str(nsm))
    #print(nsm.rstrip())
    eids = all_eids_pids[:, 0]
    ids = np.array(nsm_ids, dtype='int32')
    if debug:  # pragma: no cover
        model.log.debug('  ids = %s' % ids)
        model.log.debug('  eids = %s' % eids)
        model.log.debug('  is_area_array = %s' % is_area_array)

    isort = np.searchsorted(eids, ids)

    eids_pids = all_eids_pids
    #area = area_length[is_area]
    #length = area_length[~is_area]
    if debug:  # pragma: no cover
        model.log.debug('  area_length = %s' % area_length)
    #print('  area =', area)
    #print('  length =', length)
    #area_length =
    if debug:  # pragma: no cover
        model.log.debug(nsm.nsm_type)
        model.log.debug(eids_pids)
    #pids = all_eids_pids[:, 1]

    #print('  pids =', pids)
    #print('  area =', area)
    #print('  length =', length)

    ipassed = isort != len(eids)
    if debug:  # pragma: no cover
        model.log.debug('  isort_1 = %s' % isort)
        model.log.debug('  ipassed = %s' % ipassed)

    isort = isort[ipassed]
    if len(isort) == 0:
        model.log.warning('  *no ids found')
        return mass

    if debug:  # pragma: no cover
        model.log.debug('  isort_2 =' % isort)
        model.log.debug('  eids[isort] =' % eids[isort])

    assert len(eids) == len(area_length), 'len(eids)=%s len(area_length)=%s' % (len(eids), len(area_length))
    assert len(eids) == len(is_area_array), 'len(eids)=%s len(area_length)=%s' % (len(eids), len(is_area_array))
    iwhere = eids[isort] == ids

    eids_actual = eids[isort][iwhere]
    area_length_actual = area_length[isort][iwhere]
    is_area_actual = is_area_array[isort][iwhere]
    if len(np.unique(is_area_actual)) != 1:
        #model.log.error('  eids = %s' % eids_actual)
        #model.log.error('  area_length_actual = %s' % area_length_actual)
        #model.log.error('  is_area_actual = %s' % is_area_actual)
        msg = 'Mixed Line/Area element types for:\n%s' % str(nsm)
        for eid in eids_actual:
            msg += str(model.elements[eid])
        raise RuntimeError(msg)

    is_area = is_area_actual[0]
    if debug:  # pragma: no cover
        #print('  is_area_actual =', is_area_actual)
        if is_area:
            model.log.debug('  area!')
        else:
            model.log.debug('  length!')
    nsm_centroid = nsm_centroids[isort, :][iwhere, :]

    # this is area or length depending on eid type
    cgi = area_length_actual[:, np.newaxis] * nsm_centroid

    if debug:  # pragma: no cover
        model.log.debug('  nsm_centroid = %s' % nsm_centroid)
        model.log.debug('  area_length_actual = %s' % area_length_actual)
        model.log.debug('  cgi = %s' % cgi)

    #if is_area:
        #word = 'area'
    #else:
        #word = 'length'

    if divide_by_sum:
        area_sum = area_length_actual.sum()
        #area_sum_str = '%s_sum=%s ' % (word, area_sum)
        area_length_actual2 = area_length_actual / area_sum
    else:
        #area_sum = 1.
        #area_sum_str = ''
        area_length_actual2 = area_length_actual

    for eid, area_lengthi, centroid in zip(eids_actual, area_length_actual2, nsm_centroid):
        massi = nsm_value * area_lengthi
        #if debug:  # pragma: no cover
            #print('  eid=%s %si=%s %snsm_value=%s mass=%s' % (
                #eid, word, area_lengthi, area_sum_str, nsm_value, massi))
        mass = _increment_inertia(centroid, reference_point, massi, mass, cg, I)
    return mass

def _get_sym_axis(model, sym_axis):
    """update the sym_axis"""
    if isinstance(sym_axis, str):
        sym_axis_set = {sym_axis.lower()}
    elif isinstance(sym_axis, (list, tuple)):
        # basically overwrite the existing values on the AERO/AEROS card
        sym_axis_set = {sym_axisi.lower() for sym_axisi in sym_axis}
    else:
        sym_axis_set = set()

        # The symmetry flags on the AERO/AEROS must be the same, so
        # it doesn't matter which we one pick.  However, they might
        # not both be defined.
        #
        # Anti-symmetry refers to load, not geometry.  Geometry is
        # always symmetric.
        #
        if model.aero is not None:
            if model.aero.is_symmetric_xy or model.aero.is_anti_symmetric_xy:
                sym_axis_set.add('xy')
            if model.aero.is_symmetric_xz or model.aero.is_anti_symmetric_xz:
                sym_axis_set.add('xz')

        if model.aeros is not None:
            if model.aeros.is_symmetric_xy or model.aeros.is_anti_symmetric_xy:
                sym_axis_set.add('xy')
            if model.aeros.is_symmetric_xz or model.aeros.is_anti_symmetric_xz:
                sym_axis_set.add('xz')

    is_no = 'no' in sym_axis_set
    if is_no and len(sym_axis_set) > 1:
        raise RuntimeError('no can only be used by itself; sym_axis=%s' % (
            str(list(sym_axis_set))))
    for sym_axisi in sym_axis_set:
        if sym_axisi.lower() not in ['no', 'xy', 'yz', 'xz']:
            msg = 'sym_axis=%r is invalid; sym_axisi=%r; allowed=[no, xy, yz, xz]' % (
                sym_axis, sym_axisi)
            raise RuntimeError(msg)
    return list(sym_axis_set)

def _apply_mass_symmetry(model, sym_axis, scale, mass, cg, inertia):
    """
    Scales the mass & moement of inertia based on the symmetry axes
    and the PARAM WTMASS card

    """
    sym_axis = _get_sym_axis(model, sym_axis)

    if sym_axis:
        # either we figured sym_axis out from the AERO cards or the user told us
        model.log.debug('Mass/MOI sym_axis = %r' % sym_axis)

        if 'xz' in sym_axis:
            # y inertias are 0
            cg[1] = 0.0
            mass *= 2.0
            inertia[0] *= 2.0
            inertia[1] *= 2.0
            inertia[2] *= 2.0
            inertia[3] *= 0.0  # Ixy
            inertia[4] *= 2.0  # Ixz; no y
            inertia[5] *= 0.0  # Iyz

        if 'xy' in sym_axis:
            # z inertias are 0
            cg[2] = 0.0
            mass *= 2.0
            inertia[0] *= 2.0
            inertia[1] *= 2.0
            inertia[2] *= 2.0
            inertia[3] *= 2.0  # Ixy; no z
            inertia[4] *= 0.0  # Ixz
            inertia[5] *= 0.0  # Iyz

        if 'yz' in sym_axis:
            # x inertias are 0
            cg[0] = 0.0
            mass *= 2.0
            inertia[0] *= 2.0
            inertia[1] *= 2.0
            inertia[2] *= 2.0
            inertia[3] *= 0.0  # Ixy
            inertia[4] *= 0.0  # Ixz
            inertia[5] *= 2.0  # Iyz; no x

    wtmass = model.wtmass
    if scale is None:
        scale = wtmass
        if scale != 1.0:
            model.log.info('WTMASS scale = %r' % scale)
    mass *= scale
    inertia *= scale
    return (mass, cg, inertia)

#-------------------------------------------------------------------------------
def mass_properties_breakdown(model, element_ids=None, mass_ids=None, nsm_id=None,
                              reference_point=None,
                              sym_axis=None, scale=None, inertia_reference='cg',
                              xyz_cid0_dict=None, debug=False):
    """Gets an incomplete breakdown the mass properties on a per element basis"""
    reference_point, is_cg = _update_reference_point(
        model, reference_point, inertia_reference)
    #print('is_cg =', is_cg)

    #mass, cg, inertia
    #[Ixx, Iyy, Izz, Ixy, Ixz, Iyz]
    out = model.get_xyz_in_coord_array(cid=0, fdtype='float64', idtype='int32')
    nid_cp_cd, xyz_cid0, xyz_cp, unused_icd_transform, unused_icp_transform = out
    del xyz_cp
    xyz_mean = xyz_cid0.mean(axis=0)
    assert len(xyz_mean) == 3, xyz_mean.shape
    reference_point = np.array([xyz_mean[0], 0., 0.], dtype='float64')
    all_nids = nid_cp_cd[:, 0]

    ncoords = len(model.coords)
    cids = np.zeros(ncoords, dtype='int32')
    coords = np.zeros((ncoords, 3, 3), dtype='float64')
    iaxes = np.zeros((ncoords, 3), dtype='float64')
    for icid, (cid, coord) in zip(count(), sorted(model.coords.items())):
        cids[icid] = cid
        iaxes[icid, :] = coord.i
        coords[icid, :, :] = coord.beta()
        icid += 1


    skip_elements = {
        'CELAS1', 'CELAS2', 'CELAS3', 'CELAS4',
        'CDAMP1', 'CDAMP2', 'CDAMP3', 'CDAMP4',
        #'CHBDYP', 'CHBDYG', 'CRAC2D', 'CRAC3D',
    }
    not_finished_elements = {
        'CBEND',
    }
    skip_elements.update(not_finished_elements)
    eids_dict = defaultdict(list)
    nids_dict = defaultdict(list)
    pids_dict = defaultdict(list)
    mids_dict = defaultdict(list)
    theta_mcid_dict = defaultdict(list)
    g0_dict = defaultdict(list)
    x_dict = defaultdict(list)
    offt_dict = defaultdict(list)
    #element_nsm_dict = defaultdict(list)

    skipped_etypes = set()
    nan = np.full(3, np.nan, dtype='float64')
    for eid, elem in sorted(model.elements.items()):
        etype = elem.type
        if etype in NO_MASS:
            continue
        if etype in ['CQUAD4', 'CTRIA3', 'CQUAD8', 'CTRIA6', 'CTRIAR', 'CQUADR', 'CQUAD', 'CTRIAX']:
            nids_dict[etype].append(elem.nodes)
            pids_dict[etype].append(elem.pid)
            thetai = elem.theta_mcid
            theta_mcid_dict[etype].append(thetai)
        elif etype == 'CTRIAX6':
            nids_dict[etype].append(elem.nodes)
            pids_dict[etype].append(elem.pid)
        elif etype == 'CBEAM':
            if elem.bit is not None:
                continue
            nids_dict[etype].append(elem.nodes)
            pids_dict[etype].append(elem.pid)
            g0_dict[etype].append(-1 if elem.g0 is None else elem.g0)
            x_dict[etype].append(nan if elem.g0 is not None else elem.x)
            offt_dict[etype].append(elem.offt)
        elif etype == 'CBAR':
            nids_dict[etype].append(elem.nodes)
            pids_dict[etype].append(elem.pid)
            g0_dict[etype].append(-1 if elem.g0 is None else elem.g0)
            x_dict[etype].append(nan if elem.g0 is not None else elem.x)
            offt_dict[etype].append(elem.offt)

        elif etype in ['CTETRA', 'CPENTA', 'CHEXA', 'CPYRAM']:
            nids_dict[etype].append(elem.nodes)
            pids_dict[etype].append(elem.pid)
        elif etype in ['CROD', 'CTUBE']:
            nids_dict[etype].append(elem.nodes)
            pids_dict[etype].append(elem.pid)
        elif etype == 'CONROD':
            nids_dict[etype].append(elem.nodes)
            mids_dict[etype].append(elem.mid)
            #element_nsm_dict[etype].append(elem.nsm)
        elif etype in 'CSHEAR':
            nids_dict[etype].append(elem.nodes)
            pids_dict[etype].append(elem.pid)
        #elif etype in ['??']:
            #lazy_eids.append(elem.eid)
            #lazy_mass.append(elem.Mass())
            #lazy_cg.append(elem.center_of_mass())
        else:
            skipped_etypes.add(etype)
            continue
        eids_dict[etype].append(eid)

    nmasses = len(model.masses)
    eids_mass = []
    # [mass + nsm, mass, nsm], [x, y, z], [Ixx, Iyy, Izz, Ixy, Ixz, Iyz]
    data_mass = np.zeros((nmasses, 12), dtype='float64')
    ieid = 0
    for eid, elem in sorted(model.masses.items()):
        etype = elem.type
        if etype in NO_MASS:
            continue

        if etype in ['CONM1', 'CONM2', 'CMASS1', 'CMASS2', 'CMASS3', 'CMASS4']:
            eids_mass.append(eid)
            m = elem.Mass()
            centroid = elem.Centroid()
            data_mass[ieid, [0, 1, 2]] = [m, m, 0.]
            data_mass[ieid, 3:6] = centroid
            ieid += 1
        else:
            skipped_etypes.add(etype)
    if len(eids_mass):
        eids_mass = np.hstack(eids_mass)

    if skipped_etypes:
        model.log.warning('skipping mass_properties_breakdown for %s' % skipped_etypes)
    #---------------------------------------------------------------------------
    all_eids = []
    for eids in eids_dict.values():
        all_eids.append(eids)
    if len(all_eids):
        all_eids = np.hstack(all_eids)
        all_eids.sort()
    nelements = len(all_eids)
    if nelements + nmasses == 0:
        return 0., [0., 0., 0], [0., 0., 0., 0., 0., 0.], None, None, None

    all_mids = []
    for mids in mids_dict.values():
        all_mids.append(mids)
    if len(all_mids):
        all_mids = np.hstack(all_mids)
        all_mids.sort()

    dicts = _breakdown_property_dicts(model)
    (pids_per_length_dict, mass_per_length_dict, nsm_per_length_dict,
     pids_per_area_dict, mass_per_area_dict, nsm_per_area_dict, thickness_dict,
     pids_per_volume_dict, mass_per_volume_dict,
     #e2_dict, e3_dict,
     ) = dicts

    #xaxis = np.array([1., 0., 0.])
    #yaxis = np.array([0., 1., 0.])
    #zaxis = np.array([0., 0., 1.])

    def _get_mass(etype, eids, nids):
        nelementsi = nids.shape[0]
        #print(etype, nelementsi, eids)
        #Ax = 0.
        #Ay = 0.
        #Az = 0.
        if etype == 'CROD':
            pids = pids_dict[etype]
            assert len(pids) > 0, pids
            all_pids = np.array(pids_per_length_dict['PROD'], dtype='int32')
            mass_per_length = np.array(mass_per_length_dict['PROD'], dtype='int32')
            nsm_per_length = np.array(nsm_per_length_dict['PROD'], dtype='int32')
            assert len(nsm_per_length) > 0, nsm_per_length

            ipids = np.searchsorted(all_pids, pids)
            inids = np.searchsorted(all_nids, nids.ravel()).reshape(nelementsi, 2)
            p1 = xyz_cid0[inids[:, 0], :]
            p2 = xyz_cid0[inids[:, 1], :]
            iaxis = p2 - p1
            length = np.linalg.norm(iaxis, axis=1)
            centroid = (p1 + p2) / 2.
            assert len(length) == nelementsi, 'len(length)=%s nelementsi=%s' % (len(length), nelementsi)

            #e2i = np.array(e2_dict['rod'], dtype='float64')
            #e2 = e2i[ipids, :, :]
            #exx = eyy = ezz = 1. / e2[:, 0, 0]

            mpl = mass_per_length[ipids]
            npl = nsm_per_length[ipids]
            mass = mpl * length
            nsm = npl * length
        elif etype == 'CONROD':
            mids = mids_dict[etype]
            assert len(mids) > 0, mids
            mass_per_length = np.array(mass_per_length_dict['CONROD'], dtype='int32')
            nsm_per_length = np.array(nsm_per_length_dict['CONROD'], dtype='int32')
            #assert len(mass_per_length) > 0, mass_per_length
            #assert len(nsm_per_length) > 0, nsm_per_length

            #imids = np.searchsorted(all_mids, mids)
            inids = np.searchsorted(all_nids, nids.ravel()).reshape(nelementsi, 2)
            p1 = xyz_cid0[inids[:, 0], :]
            p2 = xyz_cid0[inids[:, 1], :]
            iaxis = p2 - p1
            length = np.linalg.norm(iaxis, axis=1)
            centroid = (p1 + p2) / 2.
            assert len(length) == nelementsi, 'len(length)=%s nelementsi=%s' % (len(length), nelementsi)

            #e2i = np.array(e2_dict['rod'], dtype='float64')
            #e2 = e2i[imids, :, :]
            #exx = eyy = ezz = 1. / e2[:, 0, 0]
            #exx = eyy = ezz = 1.

            mpl = npl = 0.
            #mpl = mass_per_length[imids]
            #npl = nsm_per_length[imids]
            mass = mpl * length
            nsm = npl * length
        elif etype == 'CTUBE':
            pids = pids_dict[etype]
            assert len(pids) > 0, pids
            all_pids = np.array(pids_per_length_dict['PTUBE'], dtype='int32')
            mass_per_length = np.array(mass_per_length_dict['PTUBE'], dtype='int32')
            nsm_per_length = np.array(nsm_per_length_dict['PTUBE'], dtype='int32')
            assert len(nsm_per_length) > 0, nsm_per_length

            ipids = np.searchsorted(all_pids, pids)
            inids = np.searchsorted(all_nids, nids.ravel()).reshape(nelementsi, 2)
            p1 = xyz_cid0[inids[:, 0], :]
            p2 = xyz_cid0[inids[:, 1], :]
            iaxis = p2 - p1
            length = np.linalg.norm(iaxis, axis=1)
            centroid = (p1 + p2) / 2.
            assert len(length) == nelementsi, 'len(length)=%s nelementsi=%s' % (len(length), nelementsi)
            #e2i = np.array(e2_dict['tube'], dtype='float64')
            #e2 = e2i[ipids, :, :]
            #exx = eyy = ezz = 1. / e2[:, 0, 0]

            mpl = mass_per_length[ipids]
            npl = nsm_per_length[ipids]
            mass = mpl * length
            nsm = npl * length
        elif etype == 'CBAR':
            pids = pids_dict[etype]
            assert len(pids) > 0, pids
            all_pids = np.array(pids_per_length_dict['bar'], dtype='int32')
            #g0 = np.array(g0_dict['CBAR'], dtype='int32')
            #offt = np.array(offt_dict['CBAR'], dtype='|S3')
            x = np.vstack(x_dict['CBAR'])

            #jaxis = _bar_axes(all_nids, xyz_cid0, g0, offt, x, nelementsi)

            mass_per_length = np.array(mass_per_length_dict['bar'])
            nsm_per_length = np.array(nsm_per_length_dict['bar'])
            assert len(nsm_per_length) > 0, nsm_per_length

            ipids = np.searchsorted(all_pids, pids)
            inids = np.searchsorted(all_nids, nids.ravel()).reshape(nelementsi, 2)
            p1 = xyz_cid0[inids[:, 0], :]
            p2 = xyz_cid0[inids[:, 1], :]
            iaxis = p2 - p1
            length = np.linalg.norm(iaxis, axis=1)
            centroid = (p1 + p2) / 2.
            assert len(length) == nelementsi, 'len(length)=%s nelementsi=%s' % (len(length), nelementsi)

            #zaxis = np.cross(iaxis, jaxis)
            #telem = np.zeros((nelementsi, 3, 3), dtype='float64')
            #telem[:, 0, :] = iaxis
            #telem[:, 1, :] = jaxis
            #telem[:, 2, :] = zaxis

            # e2i is (npids,3,3)
            # e2 is (nelementsi,3,3)
            #e2i = np.array(e2_dict['bar'], dtype='float64')
            #e2 = e2i[ipids, :, :]
            #assert e2.shape[0] == nelementsi

            #telem = transform_shell_material_coordinate_system(
                #cids, iaxes, theta_mcid, normal, p1, p2)
            # [T^T][e][T]
            #print(e2.shape, telem.shape)
            #exx, eyy, ezz = _breakdown_transform_shell(e2, telem)

            mpl = mass_per_length[ipids]
            npl = nsm_per_length[ipids]
            mass = mpl * length
            nsm = npl * length
        elif etype == 'CBEAM':
            pids = pids_dict[etype]
            assert len(pids) > 0, pids
            all_pids = np.array(pids_per_length_dict['beam'], dtype='int32')
            #g0 = np.array(g0_dict['CBEAM'], dtype='int32')
            #offt = np.array(offt_dict['CBEAM'], dtype='|S3')
            x = np.vstack(x_dict['CBEAM'])

            #jaxis = _bar_axes(all_nids, xyz_cid0, g0, offt, x, nelementsi)

            mass_per_length = np.array(mass_per_length_dict['beam'])
            nsm_per_length = np.array(nsm_per_length_dict['beam'])
            #print(nsm_per_length_dict)
            assert len(nsm_per_length) > 0, nsm_per_length
            ipids = np.searchsorted(all_pids, pids)
            #print(all_pids, pids, ipids)
            inids = np.searchsorted(all_nids, nids.ravel()).reshape(nelementsi, 2)
            p1 = xyz_cid0[inids[:, 0], :]
            p2 = xyz_cid0[inids[:, 1], :]
            iaxis = p2 - p1
            length = np.linalg.norm(iaxis, axis=1)
            centroid = (p1 + p2) / 2.
            assert len(length) == nelementsi, 'len(length)=%s nelementsi=%s' % (len(length), nelementsi)

            #zaxis = np.cross(iaxis, jaxis)
            #telem = np.zeros((nelementsi, 3, 3), dtype='float64')
            #telem[:, 0, :] = iaxis
            #telem[:, 1, :] = jaxis
            #telem[:, 2, :] = zaxis

            # e2i is (npids,3,3)
            # e2 is (nelementsi,3,3)
            #print(e2_dict)
            #e2i = np.array(e2_dict['beam'], dtype='float64')
            #assert len(e2i) > 0, e2_dict
            #e2 = e2i[ipids, :, :]
            #assert e2.shape[0] == nelementsi


            #telem = transform_shell_material_coordinate_system(
                #cids, iaxes, theta_mcid, normal, p1, p2)
            # [T^T][e][T]
            #print(e2.shape, telem.shape)
            #exx, eyy, ezz = _breakdown_transform_shell(e2, telem)

            mpl = mass_per_length[ipids]
            npl = nsm_per_length[ipids]
            mass = mpl * length
            nsm = npl * length

        elif etype in ['CTRIA3', 'CTRIA6', 'CTRIAR', ]:
            centroid, mass, nsm = _breakdown_tri(
                xyz_cid0, nids, nelementsi, etype,
                all_nids,
                pids_dict, pids_per_area_dict,
                mass_per_area_dict, nsm_per_area_dict, thickness_dict)
        elif etype == 'CSHEAR':
            centroid, mass, nsm = _breakdown_cshear(
                xyz_cid0, nids, nelementsi, etype,
                all_nids,
                pids_dict, pids_per_area_dict,
                mass_per_area_dict, nsm_per_area_dict, thickness_dict)
        elif etype in ['CQUAD4', 'CQUAD8', 'CQUADR', 'CQUAD']:
            centroid, mass, nsm = _breakdown_quad(
                xyz_cid0, nids, nelementsi, etype,
                all_nids,
                pids_dict, pids_per_area_dict,
                mass_per_area_dict, nsm_per_area_dict, thickness_dict)

        elif etype == 'CTETRA':
            centroid, mass, nsm = _breakdown_ctetra(
                xyz_cid0, nids, nelementsi, etype,
                all_nids,
                pids_dict, mass_per_volume_dict)
        elif etype == 'CHEXA':
            centroid, mass, nsm = _breakdown_chexa(
                xyz_cid0, nids, nelementsi, etype,
                all_nids,
                pids_dict, mass_per_volume_dict)
        elif etype == 'CPENTA':
            centroid, mass, nsm = _breakdown_cpenta(
                xyz_cid0, nids, nelementsi, etype,
                all_nids,
                pids_dict, mass_per_volume_dict)
        elif etype == 'CPYRAM':
            centroid, mass, nsm = _breakdown_cpyram(
                xyz_cid0, nids, nelementsi, etype,
                all_nids,
                pids_dict, mass_per_volume_dict)
        else:
            model.log.warning('skipping mass_properties_breakdown for %s' % etype)
            return None, None, None
        return mass, nsm, centroid


    data = np.zeros((nelements, 12), dtype='float64')
    for etype, nids_list in nids_dict.items():
        eids = np.hstack(eids_dict[etype])
        ieids = np.searchsorted(all_eids, eids)
        try:
            nids = np.vstack(nids_list)
        except ValueError as e:
            msg = f'error exctracting mass for {etype}; nids:\n'
            for nidsi in nids_list:
                msg += f'  {nidsi}; n={len(nidsi)}\n'
            raise ValueError(msg) from e

        mass, nsm, centroid = _get_mass(etype, eids, nids)
        if mass is None:
            continue

        # [mass + nsm, mass, nsm], [x, y, z], [Ixx, Iyy, Izz, Ixy, Ixz, Iyz]
        data[ieids, 0] = mass + nsm # total mass
        data[ieids, 1] = mass
        data[ieids, 2] = nsm

        # 0         1     2    3  4  5  6     7    8    9    10   11   12
        #total_mass, mass, nsm, x, y, z, [Ixx, Iyy, Izz, Ixy, Ixz, Iyz]
        data[ieids, 3:6] = centroid
        #data[ieids, 6:11] = inertia
        #data[ieids, 12] = Ax
        #data[ieids, 13] = Ay
        #data[ieids, 14] = Az

    if nmasses and nelements:
        # [mass + nsm, mass, nsm], [x, y, z], [Ixx, Iyy, Izz, Ixy, Ixz, Iyz]
        data = np.vstack([data, data_mass])
    elif nmasses:
        data = data_mass
    #x2 = x * x
    #y2 = y * y
    #z2 = z * z
    #I[0] += m * (y2 + z2)  # Ixx
    #I[1] += m * (x2 + z2)  # Iyy
    #I[2] += m * (x2 + y2)  # Izz
    #I[3] += m * x * y      # Ixy
    #I[4] += m * x * z      # Ixz
    #I[5] += m * y * z      # Iyz
    ix = 3
    iy = 4
    iz = 5
    total_mass = data[:, 0]
    x = data[:, ix]
    y = data[:, iy]
    z = data[:, iz]
    x2 = x ** 2
    y2 = y ** 2
    z2 = z ** 2

    #  mass moi
    data[:, 6] = total_mass * (y2 + z2) # ixx
    data[:, 7] = total_mass * (x2 + z2) # iyy
    data[:, 8] = total_mass * (x2 + y2) # izz
    data[:, 9] = total_mass * (x * y) # ixy
    data[:, 10] = total_mass * (x * z) # ixz
    data[:, 11] = total_mass * (y * z) # iyz

    mass = data[:, :3]
    cg = data[:, 3:6]
    inertia = data[:, 6:12]
    #axyz = data[:, 12:15]
    assert mass.shape[1] == 3
    assert cg.shape[1] == 3
    assert inertia.shape[1] == 6

    #  area
    #ax = axyz[:, 0]
    #ay = axyz[:, 1]
    #az = axyz[:, 2]

    # only transform if we're calculating the inertia about the cg
    total_mass_overall = total_mass.sum()
    inertia_overall = inertia.sum(axis=0)
    assert len(inertia_overall) == 6, len(inertia_overall)

    cg_overall = np.zeros(3, dtype='float64')
    if total_mass_overall == 0.:
        return total_mass_overall, cg_overall, inertia_overall, mass, cg, inertia
        #raise RuntimeError('total_mass_overall=%s' % total_mass_overall)
    cg_overall = (total_mass[:, np.newaxis] * cg).sum(axis=0) / total_mass_overall

    #if is_cg:
        #xyz_ref = reference_point
        #xyz_ref2 = (cg[:, 0], cg[:, 1], cg[:, 2])
        #inertia2 = (inertia[:, 0], inertia[:, 1], inertia[:, 2],
                    #inertia[:, 3], inertia[:, 4], inertia[:, 5], )
        #print('cg_overall =', cg_overall)
        #inertia = transform_inertia(mass, cg_overall, xyz_ref, xyz_ref2, inertia2)

    make_plot = False
    if make_plot:
        plot_inertia(total_mass, cg, inertia)
    return total_mass_overall, cg_overall, inertia_overall, mass, cg, inertia

def _breakdown_tri(xyz_cid0, nids, nelementsi, etype,
                   all_nids,
                   pids_dict, pids_per_area_dict,
                   mass_per_area_dict, nsm_per_area_dict, thickness_dict):
    # no offsets
    nids2 = nids[:, :3]
    pids = np.array(pids_dict[etype], dtype='int32')
    assert len(pids) > 0, pids
    all_pids = np.array(pids_per_area_dict['shell'], dtype='int32')
    mass_per_area = np.array(mass_per_area_dict['shell'])
    nsm_per_area = np.array(nsm_per_area_dict['shell'])
    thickness = np.array(thickness_dict['shell'])
    assert len(mass_per_area) > 0, mass_per_area_dict
    assert len(thickness) > 0, thickness

    ipids = np.searchsorted(all_pids, pids)
    inids = np.searchsorted(all_nids, nids2.ravel()).reshape(nelementsi, 3)
    p1 = xyz_cid0[inids[:, 0], :]
    p2 = xyz_cid0[inids[:, 1], :]
    p3 = xyz_cid0[inids[:, 2], :]

    # normal is correct; matters for +rotation and offsets
    v12 = p1 - p2
    v13 = p1 - p3
    normal = np.cross(v12, v13)
    ni = np.linalg.norm(normal, axis=1)
    normal /= ni[:, np.newaxis]
    area = 0.5 * ni
    assert len(area) == nelementsi, 'len(area)=%s nelementsi=%s' % (len(area), nelementsi)
    centroid = (p1 + p2 + p3) / 3.

    # e2i is (npids,3,3)
    # e2 is (nelementsi,3,3)
    #e2i = np.array(e2_dict['shell'], dtype='float64')
    #e2 = e2i[ipids, :, :]
    #assert e2.shape[0] == nelementsi

    # [T^T][e][T]
    #theta_mcid = theta_mcid_dict[etype]
    #telem = transform_shell_material_coordinate_system(
        #cids, iaxes, theta_mcid, normal, p1, p2)
    #exx, eyy, ezz = _breakdown_transform_shell(e2, telem)
    #exx = eyy = ezz = 1.
    mpa = mass_per_area[ipids]
    npa = nsm_per_area[ipids]
    mass = mpa * area
    nsm = npa * area
    return centroid, mass, nsm

def _breakdown_quad(xyz_cid0, nids, nelementsi, etype,
                    all_nids,
                    pids_dict, pids_per_area_dict,
                    mass_per_area_dict, nsm_per_area_dict, thickness_dict):
    # no offsets
    nids2 = nids[:, :4]
    assert  nids2.shape[0] == nelementsi, nids2.shape
    pids = np.array(pids_dict[etype], dtype='int32')
    assert len(pids) > 0, pids
    all_pids = np.array(pids_per_area_dict['shell'])
    mass_per_area = np.array(mass_per_area_dict['shell'])
    nsm_per_area = np.array(nsm_per_area_dict['shell'])
    thickness = np.array(thickness_dict['shell'])
    assert len(mass_per_area) > 0, mass_per_area_dict
    assert len(thickness) > 0, thickness

    ipids = np.searchsorted(all_pids, pids)
    inids = np.searchsorted(all_nids, nids2.ravel()).reshape(nelementsi, 4)
    p1 = xyz_cid0[inids[:, 0], :]
    p2 = xyz_cid0[inids[:, 1], :]
    p3 = xyz_cid0[inids[:, 2], :]
    p4 = xyz_cid0[inids[:, 3], :]

    # normal is correct; matters for +rotation and offsets
    v13 = p1 - p3
    v24 = p2 - p4
    normal = np.cross(v13, v24)
    ni = np.linalg.norm(normal, axis=1)
    normal /= ni[:, np.newaxis]
    area = 0.5 * ni
    assert len(area) == nelementsi, 'len(area)=%s nelementsi=%s' % (len(area), nelementsi)
    centroid = (p1 + p2 + p3 + p4) / 4.

    # e2i is (npids,3,3)
    # e2 is (nelementsi,3,3)
    #e2i = np.array(e2_dict['shell'], dtype='float64')
    #e2 = e2i[ipids, :, :]
    #assert e2.shape[0] == nelementsi

    # [T^T][e][T]
    #theta_mcid = theta_mcid_dict[etype]
    #telem = transform_shell_material_coordinate_system(
        #cids, iaxes, theta_mcid, normal, p1, p2)
    #exx = eyy = ezz = 1.
    #exx, eyy, ezz = _breakdown_transform_shell(e2, telem)
    mpa = mass_per_area[ipids]
    npa = nsm_per_area[ipids]
    mass = mpa * area
    nsm = npa * area

    # assume the panel is square to calculate w; then multiply by t to get tw
    assert len(area) > 0, area
    assert len(thickness) > 0, thickness
    #tw = thickness[ipids] * np.sqrt(area)
    #Ax = tw * norm(cross(xaxis, normal), axis=1)
    #Ay = tw * norm(cross(yaxis, normal), axis=1)
    #Az = tw * norm(cross(zaxis, normal), axis=1)
    return centroid, mass, nsm

def _breakdown_cshear(xyz_cid0, nids, nelementsi, etype,
                      all_nids,
                      pids_dict, pids_per_area_dict,
                      mass_per_area_dict, nsm_per_area_dict, thickness_dict):
    pids = np.array(pids_dict[etype], dtype='int32')
    assert len(pids) > 0, pids
    all_pids = np.array(pids_per_area_dict['shear'])
    mass_per_area = np.array(mass_per_area_dict['shear'])
    nsm_per_area = np.array(nsm_per_area_dict['shear'])
    thickness = np.array(thickness_dict['shear'])
    assert len(mass_per_area) > 0, mass_per_area_dict
    ipids = np.searchsorted(all_pids, pids)
    inids = np.searchsorted(all_nids, nids.ravel()).reshape(nelementsi, 4)
    p1 = xyz_cid0[inids[:, 0], :]
    p2 = xyz_cid0[inids[:, 1], :]
    p3 = xyz_cid0[inids[:, 2], :]
    p4 = xyz_cid0[inids[:, 3], :]

    # normal is correct; matters for +rotation and offsets
    v13 = p1 - p3
    v24 = p2 - p4
    normal = np.cross(v13, v24)
    ni = np.linalg.norm(normal, axis=1)
    normal /= ni[:, np.newaxis]
    area = 0.5 * ni
    assert len(area) == nelementsi, 'len(area)=%s nelementsi=%s' % (len(area), nelementsi)
    centroid = (p1 + p2 + p3 + p4) / 4.
    #n = _normal(n1 - n3, n2 - n4)

    # e2i is (npids,3,3)
    # e2 is (nelementsi,3,3)
    #e2i = np.array(e2_dict['shear'], dtype='float64')
    #e2 = e2i[ipids, :, :]
    #assert e2.shape[0] == nelementsi

    # [T^T][e][T]
    #theta_mcid = theta_mcid_dict[etype]
    #telem = transform_shell_material_coordinate_system(
        #cids, iaxes, theta_mcid, normal, p1, p2)
    #exx = eyy = ezz = 1.
    #exx, eyy, ezz = _breakdown_transform_shell(e2, telem)
    mpa = mass_per_area[ipids]
    npa = nsm_per_area[ipids]
    mass = mpa * area
    nsm = npa * area

    # assume the panel is square to calculate w; then multiply by t to get tw
    tw = thickness[ipids] * np.sqrt(area)
    assert len(tw) > 0, tw
    #Ax = tw * norm(cross(xaxis, normal), axis=1)
    #Ay = tw * norm(cross(yaxis, normal), axis=1)
    #Az = tw * norm(cross(zaxis, normal), axis=1)
    return centroid, mass, nsm

def _breakdown_ctetra(xyz_cid0, nids, nelementsi, etype,
                      all_nids,
                      pids_dict, mass_per_volume_dict):
    nids2 = nids[:, :4]
    pids = np.array(pids_dict[etype], dtype='int32')
    assert len(pids) > 0, pids
    rho = np.array(mass_per_volume_dict['PSOLID']) # rho
    assert  len(rho) > 0., rho

    inids = np.searchsorted(all_nids, nids2.ravel()).reshape(nelementsi, 4)
    p1 = xyz_cid0[inids[:, 0], :]
    p2 = xyz_cid0[inids[:, 1], :]
    p3 = xyz_cid0[inids[:, 2], :]
    p4 = xyz_cid0[inids[:, 3], :]
    centroid = (p1 + p2 + p3 + p4) / 4.
    a = p1 - p4
    b = cross(p2 - p4, p3 - p4)
    #volume = -dot(a, b) / 6.
    #volume = -np.tensordot(a, b, axes=1)
    volume = -np.einsum('ij,ij->i', a, b) / 6.
    mass = rho * volume
    nsm = np.zeros(nelementsi, dtype='float64')
    #exx = eyy = ezz = 1.
    #e2 = None
    return centroid, mass, nsm

def _breakdown_chexa(xyz_cid0, nids, nelementsi, etype,
                     all_nids,
                     pids_dict, mass_per_volume_dict):
    nids2 = nids[:, :8]
    pids = np.array(pids_dict[etype], dtype='int32')
    assert len(pids) > 0, pids
    rho = np.array(mass_per_volume_dict['PSOLID']) # rho
    assert  len(rho) > 0., rho

    inids = np.searchsorted(all_nids, nids2.ravel()).reshape(nelementsi, 8)
    p1 = xyz_cid0[inids[:, 0], :]
    p2 = xyz_cid0[inids[:, 1], :]
    p3 = xyz_cid0[inids[:, 2], :]
    p4 = xyz_cid0[inids[:, 3], :]
    p5 = xyz_cid0[inids[:, 4], :]
    p6 = xyz_cid0[inids[:, 5], :]
    p7 = xyz_cid0[inids[:, 6], :]
    p8 = xyz_cid0[inids[:, 7], :]

    v13 = p1 - p3
    v24 = p2 - p4
    normal = np.cross(v13, v24)
    ni = np.linalg.norm(normal, axis=1)
    c1 = (p1 + p2 + p3 + p4) / 4.
    area1 = 0.5 * ni

    v57 = p5 - p7
    v68 = p6 - p7
    normal = np.cross(v57, v68)
    ni = np.linalg.norm(normal, axis=1)
    c2 = (p5 + p6 + p7 + p8) / 4.
    area2 = 0.5 * ni

    centroid = (c1 + c2) / 2.
    volume = (area1 + area2) / 2. * norm(c1 - c2, axis=1)
    assert len(volume) == nelementsi, len(volume)

    mass = rho * volume
    nsm = np.zeros(nelementsi, dtype='float64')
    #exx = eyy = ezz = 1.
    #e2 = None
    return centroid, mass, nsm

def _breakdown_cpenta(xyz_cid0, nids, nelementsi, etype,
                      all_nids,
                      pids_dict, mass_per_volume_dict):
    nids2 = nids[:, :6]
    pids = np.array(pids_dict[etype], dtype='int32')
    assert len(pids) > 0, pids
    rho = np.array(mass_per_volume_dict['PSOLID']) # rho
    assert  len(rho) > 0., rho

    inids = np.searchsorted(all_nids, nids2.ravel()).reshape(nelementsi, 6)
    p1 = xyz_cid0[inids[:, 0], :]
    p2 = xyz_cid0[inids[:, 1], :]
    p3 = xyz_cid0[inids[:, 2], :]
    p4 = xyz_cid0[inids[:, 3], :]
    p5 = xyz_cid0[inids[:, 4], :]
    p6 = xyz_cid0[inids[:, 5], :]

    c1 = (p1 + p2 + p3) / 3.
    c2 = (p4 + p5 + p6) / 3.
    centroid = (c1 + c2) / 2.

    area1 = 0.5 * norm(cross(p3 - p1, p2 - p1), axis=1)
    area2 = 0.5 * norm(cross(p6 - p4, p5 - p4), axis=1)
    volume = (area1 + area2) / 2. * norm(c1 - c2, axis=1)
    volume = np.abs(volume)

    mass = rho * volume
    nsm = np.zeros(nelementsi, dtype='float64')
    #exx = eyy = ezz = 1.
    #e2 = None
    return centroid, mass, nsm

def _breakdown_cpyram(xyz_cid0, nids, nelementsi, etype,
                      all_nids,
                      pids_dict, mass_per_volume_dict):
    nids2 = nids[:, :5]
    pids = np.array(pids_dict[etype], dtype='int32')
    assert len(pids) > 0, pids
    rho = np.array(mass_per_volume_dict['PSOLID']) # rho
    assert  len(rho) > 0., rho

    inids = np.searchsorted(all_nids, nids2.ravel()).reshape(nelementsi, 5)
    p1 = xyz_cid0[inids[:, 0], :]
    p2 = xyz_cid0[inids[:, 1], :]
    p3 = xyz_cid0[inids[:, 2], :]
    p4 = xyz_cid0[inids[:, 3], :]
    p5 = xyz_cid0[inids[:, 4], :]

    v13 = p1 - p3
    v24 = p2 - p4
    normal = np.cross(v13, v24)
    ni = np.linalg.norm(normal, axis=1)
    c1 = (p1 + p2 + p3 + p4) / 4.
    area1 = 0.5 * ni

    centroid = (c1 + p5) / 2.

    volume = area1 / 3. * norm(c1 - p5, axis=1)
    #volume = np.abs(volume)
    #return abs(volume)
    mass = rho * volume
    nsm = np.zeros(nelementsi, dtype='float64')
    #exx = eyy = ezz = 1.
    #e2 = None
    return centroid, mass, nsm

def plot_inertia(total_mass, cg, inertia):  # pragma: no cover
    ycg = cg[:, 1]
    ixx = inertia[:, 0]
    iyy = inertia[:, 1]
    izz = inertia[:, 2]
    isort = np.argsort(ycg)
    ycg2 = ycg[isort]

    mass2 = total_mass[isort]
    ixx2 = ixx[isort]
    iyy2 = iyy[isort]
    izz2 = izz[isort]

    # cumulative moi; we integrate from the right side (the iyy2[::-1]
    # and then reverse it with cumsum(...)[::-1]
    cmass = np.cumsum(mass2[::-1])[::-1]
    cixx = np.cumsum(ixx2[::-1])[::-1]
    ciyy = np.cumsum(iyy2[::-1])[::-1]
    cizz = np.cumsum(izz2[::-1])[::-1]

    #import matplotlib
    #matplotlib.use('Qt5cairo')
    import matplotlib.pyplot as plt
    figsize = (8., 6.)
    dpi = 100
    fig = plt.figure(figsize=figsize, dpi=dpi)
    ax = fig.gca()
    ax.plot(ycg2, cmass / cmass[0], label='mass=%g' % cmass[0])
    ax.plot(ycg2, cixx / cixx[0], label='Ixx=%.3e' % cixx[0])
    ax.plot(ycg2, ciyy / ciyy[0], label='Iyy=%.3e' % ciyy[0])
    ax.plot(ycg2, cizz / cizz[0], label='Izz=%.3e' % cizz[0])
    ax.legend()
    ax.set_xlabel('y')
    ax.set_ylabel('Mass Moment of Inertia')
    ax.set_title('Mass Moment of Inertia vs Span')
    ax.grid(True)
    fig.savefig('moi vs span.png')
    #plt.show()

def _breakdown_property_dicts(model):
    """helper method"""
    pids_per_length_dict = defaultdict(list)
    mass_per_length_dict = defaultdict(list)
    nsm_per_length_dict = defaultdict(list)

    pids_per_area_dict = defaultdict(list)
    mass_per_area_dict = defaultdict(list)
    nsm_per_area_dict = defaultdict(list)
    thickness_dict = defaultdict(list)

    pids_per_volume_dict = defaultdict(list)
    mass_per_volume_dict = defaultdict(list)
    #e2_dict = defaultdict(list)
    #e3_dict = defaultdict(list)
    for pid, prop in sorted(model.properties.items()):
        ptype = prop.type
        if ptype in NO_MASS:
            continue

        if ptype in 'PROD':
            pids_per_length_dict[ptype].append(pid)
            mid_ref = prop.mid_ref
            #Sei2, unused_Sei3 = get_mat_props_S(mid_ref)
            #e2_dict['rod'].append(Sei2)
            nsm_per_length_dict[ptype].append(prop.nsm)
            mass_per_length_dict[ptype].append(mid_ref.rho * prop.A)
        elif ptype in 'PTUBE':
            pids_per_length_dict[ptype].append(pid)
            mid_ref = prop.mid_ref
            #Sei2, unused_Sei3 = get_mat_props_S(mid_ref)
            #e2_dict['tube'].append(Sei2)
            nsm_per_length_dict[ptype].append(prop.nsm)
            mass_per_length_dict[ptype].append(mid_ref.rho * prop.Area())

        elif ptype == 'PBAR':
            pids_per_length_dict['bar'].append(pid)
            mid_ref = prop.mid_ref
            rhoi = prop.Rho()
            areai = prop.Area()
            nsmi = prop.Nsm()
            #Sei2, unused_Sei3 = get_mat_props_S(mid_ref)
            #e2_dict['bar'].append(Sei2)
            nsm_per_length_dict['bar'].append(nsmi)
            mass_per_length_dict['bar'].append(areai * rhoi)
        elif ptype == 'PBARL':
            pids_per_length_dict['bar'].append(pid)
            mid_ref = prop.mid_ref
            rhoi = prop.Rho()
            areai = prop.Area()
            nsmi = prop.Nsm()
            #Sei2, unused_Sei3 = get_mat_props_S(mid_ref)
            #e2_dict['bar'].append(Sei2)
            nsm_per_length_dict['bar'].append(nsmi)
            mass_per_length_dict['bar'].append(areai * rhoi)

        elif ptype == 'PBEAM':
            pids_per_length_dict['beam'].append(pid)
            mid_ref = prop.mid_ref
            rhoi = prop.Rho()
            areai = prop.Area()
            nsmi = prop.Nsm()
            #Sei2, unused_Sei3 = get_mat_props_S(mid_ref)
            #e2_dict['beam'].append(Sei2)
            nsm_per_length_dict['beam'].append(nsmi)
            mass_per_length_dict['beam'].append(areai * rhoi)
        elif ptype == 'PBEAML':
            pids_per_length_dict['beam'].append(pid)
            mid_ref = prop.mid_ref
            rhoi = prop.Rho()
            areai = prop.Area()
            nsmi = prop.Nsm()
            #Sei2, unused_Sei3 = get_mat_props_S(mid_ref)
            #e2_dict['beam'].append(Sei2)
            nsm_per_length_dict['beam'].append(nsmi)
            mass_per_length_dict['beam'].append(areai * rhoi)
        elif prop.type == 'PBCOMP':
            pids_per_length_dict['beam'].append(pid)
            mid_ref = prop.mid_ref
            mass_per_length = prop.MassPerLength()
            nsm_per_length = prop.nsm
            #nsm_n1 = (p1 + jhat * prop.m1 + khat * prop.m2)
            #nsm_n2 = (p2 + jhat * prop.m1 + khat * prop.m2)
            #nsm_centroid = (nsm_n1 + nsm_n2) / 2.
            #Sei2, unused_Sei3 = get_mat_props_S(mid_ref)
            #e2_dict['beam'].append(Sei2)
            nsm_per_length_dict['beam'].append(nsm_per_length)
            mass_per_length_dict['beam'].append(mass_per_length)

        elif ptype == 'PSHELL':
            pids_per_area_dict['shell'].append(pid)
            mid_ref = prop.mid_ref
            rhoi = mid_ref.Rho()
            #ei2, ei3 = get_mat_props_S(mid_ref)

            #e2_dict['shell'].append(ei2)
            #e3_dict['shell'].append(ei3)

            # doesn't consider tflag
            #thickness = self.Thickness(tflag=tflag, tscales=tscales)
            thickness = prop.t

            mass_per_area_dict['shell'].append(rhoi * thickness)
            nsm_per_area_dict['shell'].append(prop.nsm)
            thickness_dict['shell'].append(thickness)
        elif ptype == 'PSHEAR':
            pids_per_area_dict['shear'].append(pid)
            mid_ref = prop.mid_ref
            rhoi = mid_ref.Rho()
            #ei2, ei3 = get_mat_props_S(mid_ref)
            #e2_dict['shear'].append(ei2)
            #e3_dict['shear'].append(ei3)

            # doesn't consider tflag
            #thickness = self.Thickness(tflag=tflag, tscales=tscales)
            thickness = prop.t

            thickness_dict['shear'].append(thickness)
            mass_per_area_dict['shear'].append(rhoi * thickness)
            nsm_per_area_dict['shear'].append(prop.nsm)
        elif ptype == 'PLPLANE':
            pids_per_area_dict['shell'].append(pid)
            mid_ref = prop.mid_ref
            rhoi = mid_ref.Rho()
            #ei2, ei3 = get_mat_props_S(mid_ref)
            #e2_dict['shell'].append(ei2)
            #e3_dict['shell'].append(ei3)

            # doesn't consider tflag
            #thickness = self.Thickness(tflag=tflag, tscales=tscales)
            thickness = 0.

            thickness_dict['shell'].append(thickness)
            mass_per_area_dict['shell'].append(rhoi * thickness)
            model.log.warning('assuming nsm for PLPLANE is 0.0')
            nsm_per_area_dict['shell'].append(0.0)
            # nsm_per_area_dict['shell'].append(prop.nsm)
        elif ptype == 'PPLANE':
            pids_per_area_dict['shell'].append(pid)
            mid_ref = prop.mid_ref
            rhoi = mid_ref.Rho()
            #ei2, ei3 = get_mat_props_S(mid_ref)
            #e2_dict['shell'].append(ei2)
            #e3_dict['shell'].append(ei3)

            thickness = prop.Thickness()
            thickness_dict['shell'].append(thickness)
            mass_per_area_dict['shell'].append(rhoi * thickness)
            nsm_per_area_dict['shell'].append(prop.nsm)
        elif ptype == 'PCOMP':
            pids_per_area_dict['shell'].append(pid)
            mids_ref = prop.mids_ref

            #rhoi = [mid_ref.Rho() for mid_ref in mids_ref]
            ksym = 2 if prop.is_symmetrical else 1
            mpai = [ksym * mid_ref.Rho() * t for mid_ref, t in zip(mids_ref, prop.thicknesses)]
            thickness = ksym * sum(prop.thicknesses)
            #nlayers = ksym * len(prop.thicknesses)
            #ti = np.zeros((nlayers, 3, 3), dtype='float64')
            #if prop.is_symmetrical:
                #theta = np.radians(np.hstack([
                    #prop.thetas[::-1],
                    #prop.thetas,
                #]))
            #else:
                #theta = np.radians(prop.thetas)

            #st = np.sin(theta)
            #ct = np.cos(theta)
            #ti[:, 0, 0] = ct
            #ti[:, 1, 1] = ct
            #ti[:, 0, 1] = st
            #ti[:, 1, 0] = -st
            #ti[:, 2, 2] = 1.
            # T = ti @ ti
            #T = np.einsum('ijk,ikj->ijk', ti, ti)

            # [o] = [Q][e]
            # [e] = [S][o]
            # the material matrix [S] rotated into the element frame
            # T = [t][t]^T
            # [Sbar] = [T_inv][S][T_inv^T]
            #S02 = np.zeros((3, 3), dtype='float64')
            #for mid_ref, sti, cti in zip(prop.mids_ref, st, ct):
                #ti = np.array([
                    #[cti, sti, 0.],
                    #[-sti, cti, 0.],
                    #[0., 0., 1.],
                #])
                # T = [t][t]^T
                #T = ti.dot(ti.T)
                #print(ti)
                #print(T)
                #Tinv = np.linalg.inv(T)
                #Sei2, unused_Sei3 = get_mat_props_S(mid_ref)

                # [Sbar] = [T_inv][S][T_inv^T]
                #S02 += Tinv.dot(Sei2.dot(Tinv.T))

                #Q2 = np.linalg.inv(Sei2)
                #Q3 = np.linalg.inv(Sei3)

            # doesn't consider tflag
            #thickness = self.Thickness(tflag=tflag, tscales=tscales)
            #thickness = prop.thickness

            #e2_dict['shell'].append(S02)
            mass_per_area_dict['shell'].append(sum(mpai))
            nsm_per_area_dict['shell'].append(prop.nsm)
            thickness_dict['shell'].append(thickness)

        elif ptype in ['PSOLID', 'PIHEX']:
            pids_per_volume_dict[ptype].append(pid)
            mid_ref = prop.mid_ref
            rho = mid_ref.rho
            assert  rho >= 0., rho
            mass_per_volume_dict[ptype].append(rho)
        else:
            model.log.warning('skipping mass_properties_breakdown for %s' % ptype)

    #for mid, material in sorted(model.materials.items()):
        #mtype = material.type
        #if mtype == 'MAT1':
            #rho = material.rho
        #elif mtype == 'MAT2':
            #rho = material.rho
        #elif mtype == 'MAT8':
            #rho = material.rho
        #else:
            #model.log.warning('skipping mass_properties_breakdown for %s' % mtype)
    #---------------------------------------------------------------------------
    dicts = (
        pids_per_length_dict, mass_per_length_dict, nsm_per_length_dict,
        pids_per_area_dict, mass_per_area_dict, nsm_per_area_dict, thickness_dict,
        pids_per_volume_dict, mass_per_volume_dict,
        #e2_dict, e3_dict,
    )
    return dicts

def _bar_axes(all_nids, xyz_cid0, g0, offt, x, nelements):
    ix = np.where(g0 == -1)[0]
    ig0 = np.where(g0 != -1)[0]

    g0 = g0[ig0]
    x = x[ix, :]
    jaxis = np.full((nelements, 3), np.nan, dtype='float64')
    if len(g0):
        offti = offt[ig0]
        nofft = len(offti)

        iggg = np.where(b'GGG' == offti)[0]
        iig0 = np.searchsorted(all_nids, g0)

        #print('iggg =', iggg, len(iggg))
        #print('x.shape', x.shape)
        #x[iggg, :]
        #jaxis[iggg, :]
        x = xyz_cid0[iig0, :]
        nggg = len(iggg)
        nofft_found = nggg
        if nggg:
            jaxis[iggg, :] = x
        if nofft != nofft_found:
            raise RuntimeError(offti)

    if len(x):
        offti = offt[ix]
        nofft = len(offti)
        #print('g0 =', g0)
        #print('x =', x)
        #print('offt =', offt)
        iggg = np.where(b'GGG' == offti)[0]
        #print('iggg =', iggg, len(iggg))
        #print('x.shape', x.shape)
        #x[iggg, :]
        #jaxis[iggg, :]
        nggg = len(iggg)
        nofft_found = nggg
        if nggg:
            jaxis[iggg, :] = x[iggg, :]
        if nofft != nofft_found:
            raise RuntimeError(offti)
    #print('xmax', x.max())
    assert np.abs(x).max() > 0., 'x has nan values...only supports GGG'
    return jaxis


#def _breakdown_transform_shell(e2, telem):
    #"""helper method for the breakdown"""
    # https://www.brown.edu/Departments/Engineering/Courses/EN224/anis_general/anis_general.htm
    #et = np.einsum('nij,njk->nik', e2, telem)
    #tet = np.einsum('nij,njk->nik', telem, et)
    #print(tet[0, :, :])

    #exx_inv = tet[:, 0, 0]
    #eyy_inv = tet[:, 1, 1]
    #ezz_inv = tet[:, 2, 2]
    #exx = np.zeros(exx_inv.shape, dtype=exx_inv.dtype)
    #eyy = np.zeros(eyy_inv.shape, dtype=eyy_inv.dtype)
    #ezz = np.zeros(ezz_inv.shape, dtype=ezz_inv.dtype)
    #assert abs(exx_inv.max()) >= 0., exx_inv
    #iexx = np.where(exx_inv > 0)[0]
    #ieyy = np.where(eyy_inv > 0)[0]
    #iezz = np.where(ezz_inv > 0)[0]
    #exx[iexx] = 1 / exx_inv[iexx]
    #eyy[ieyy] = 1 / eyy_inv[ieyy]
    #ezz[iezz] = 1 / ezz_inv[iezz]
    #print(exx.tolist())
    #return exx, eyy, ezz
