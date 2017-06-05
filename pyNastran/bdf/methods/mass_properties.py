"""
Defines:
  - mass_poperties
      get the mass & moment of inertia of the model
"""
from __future__ import print_function
from six import string_types, iteritems
from numpy import array, cross, dot
from numpy.linalg import norm
import numpy as np
from pyNastran.utils.mathematics import integrate_positive_unit_line


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
    print('dx1 = <%s, %s, %s>' % (dx1, dy1, dz1))
    print('dx2 = <%s, %s, %s>' % (dx2, dy2, dz2))

    # consistent with mass_properties, not CONM2
    print('I_ref =', I_ref)
    Ixx_ref, Iyy_ref, Izz_ref, Ixy_ref, Ixz_ref, Iyz_ref = I_ref
    Ixx2 = Ixx_ref - mass * (dx1**2 - dx2**2)
    Iyy2 = Iyy_ref - mass * (dy1**2 - dy2**2)
    Izz2 = Izz_ref - mass * (dz1**2 - dz2**2)
    Ixy2 = Ixy_ref - mass * (dx1 * dy1 - dx2 * dy2)
    Ixz2 = Ixz_ref - mass * (dx1 * dz1 - dx2 * dz2)
    Iyz2 = Iyz_ref - mass * (dy1 * dz1 - dy2 * dz2)
    I_new = np.array([Ixx2, Iyy2, Izz2, Ixy2, Ixz2, Iyz2])
    return I_new

def _mass_properties_elements_init(model, element_ids, mass_ids):
    """helper method"""
    # if neither element_id nor mass_ids are specified, use everything
    if element_ids is None and mass_ids is None:
        elements = model.elements.values()
        masses = model.masses.values()

    # if either element_id or mass_ids are specified and the other is not, use only the
    # specified ids
    else:
        if element_ids is None:
            elements = []
        else:
            elements = [element for eid, element in model.elements.items() if eid in element_ids]
        if mass_ids is None:
            masses = []
        else:
            masses = [mass for eid, mass in model.masses.items() if eid in mass_ids]
    return elements, masses

def _mass_properties(model, elements, masses, reference_point):  # pragma: no cover
    """
    Caclulates mass properties in the global system about the
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
    #Ixx Iyy Izz, Ixy, Ixz Iyz
    # precompute the CG location and make it the reference point
    I = array([0., 0., 0., 0., 0., 0., ])
    cg = array([0., 0., 0.])
    if isinstance(reference_point, string_types):
        if reference_point == 'cg':
            mass = 0.
            for pack in [elements, masses]:
                for element in pack:
                    try:
                        p = element.Centroid()
                        m = element.Mass()
                        mass += m
                        cg += m * p
                    except:
                        pass
            if mass == 0.0:
                return mass, cg, I

            reference_point = cg / mass
        else:
            # reference_point = [0.,0.,0.] or user-defined array
            pass

    mass = 0.
    cg = array([0., 0., 0.])
    for pack in [elements, masses]:
        for element in pack:
            try:
                p = element.Centroid()
            except:
                continue

            try:
                m = element.Mass()
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
            except:
                # PLPLANE
                if element.pid_ref.type == 'PSHELL':
                    model.log.warning('p=%s reference_point=%s type(reference_point)=%s' % (
                        p, reference_point, type(reference_point)))
                    raise
                model.log.warning("could not get the inertia for element/property\n%s%s" % (
                    element, element.pid_ref))

                continue

    if mass:
        cg /= mass
    return (mass, cg, I)

def _mass_properties_no_xref(model, elements, masses, reference_point):  # pragma: no cover
    """
    Caclulates mass properties in the global system about the
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

    Returns
    -------
    mass : float
        the mass of the model
    cg : (3, ) float NDARRAY
        the cg of the model as an array.
    I : (6, ) float NDARRAY
        moment of inertia array([Ixx, Iyy, Izz, Ixy, Ixz, Iyz])

    .. seealso:: self.mass_properties
    """
    #Ixx Iyy Izz, Ixy, Ixz Iyz
    # precompute the CG location and make it the reference point
    I = array([0., 0., 0., 0., 0., 0., ])
    cg = array([0., 0., 0.])
    if isinstance(reference_point, string_types):
        if reference_point == 'cg':
            mass = 0.
            for pack in [elements, masses]:
                for element in pack:
                    try:
                        p = element.Centroid_no_xref(model)
                        m = element.Mass_no_xref(model)
                        mass += m
                        cg += m * p
                    except:
                        #pass
                        raise
            if mass == 0.0:
                return mass, cg, I

            reference_point = cg / mass
        else:
            # reference_point = [0.,0.,0.] or user-defined array
            pass

    mass = 0.
    cg = array([0., 0., 0.])
    for pack in [elements, masses]:
        for element in pack:
            try:
                p = element.Centroid_no_xref(model)
            except:
                #continue
                raise

            try:
                m = element.Mass_no_xref(model)
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

    if mass:
        cg /= mass
    return (mass, cg, I)

def _mass_properties_new(model, element_ids=None, mass_ids=None, reference_point=None,
                         sym_axis=None, scale=None, xyz_cid0=None):  # pragma: no cover
    """
    half implemented, not tested, should be faster someday...
    don't use this

    Caclulates mass properties in the global system about the
    reference point.

    Parameters
    ----------
    model : BDF()
        a BDF object
    element_ids : list[int]; (n, ) ndarray, optional
        An array of element ids.
    mass_ids : list[int]; (n, ) ndarray, optional
        An array of mass ids.
    reference_point : ndarray/str/int, optional
        type : ndarray
            An array that defines the origin of the frame.
            default = <0,0,0>.
        type : str
            'cg' is the only allowed string
        type : int
            the node id
    sym_axis : str, optional
        The axis to which the model is symmetric. If AERO cards are used, this can be left blank
        allowed_values = 'no', x', 'y', 'z', 'xy', 'yz', 'xz', 'xyz'
    scale : float, optional
        The WTMASS scaling value.
        default=None -> PARAM, WTMASS is used
        float > 0.0
    xyz_cid0 : dict[nid] : xyz; default=None -> auto-calculate
        mapping of the node id to the global position

    Returns
    -------
    mass : float
        The mass of the model.
    cg : ndarray
        The cg of the model as an array.
    I : ndarray
        Moment of inertia array([Ixx, Iyy, Izz, Ixy, Ixz, Iyz]).

    I = mass * centroid * centroid

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

    Example 1
    ---------
    # mass properties of entire structure
    mass, cg, I = model.mass_properties()
    Ixx, Iyy, Izz, Ixy, Ixz, Iyz = I


    Example 2
    ---------
    # mass properties of model based on Property ID
    pids = list(model.pids.keys())
    pid_eids = model.get_element_ids_dict_with_pids(pids)

    for pid, eids in sorted(iteritems(pid_eids)):
        mass, cg, I = mass_properties(model, element_ids=eids)
    """
    #if reference_point is None:
    reference_point = array([0., 0., 0.])

    if xyz_cid0 is None:
        xyz = {}
        for nid, node in iteritems(model.nodes):
            xyz[nid] = node.get_position()
    else:
        xyz = xyz_cid0

    elements, masses = _mass_properties_elements_init(model, element_ids, mass_ids)

    #mass = 0.
    #cg = array([0., 0., 0.])
    #I = array([0., 0., 0., 0., 0., 0., ])
    #if isinstance(reference_point, string_types):
        #if reference_point == 'cg':
            #mass = 0.
            #for pack in [elements, masses]:
                #for element in pack:
                    #try:
                        #p = element.Centroid()
                        #m = element.Mass()
                        #mass += m
                        #cg += m * p
                    #except:
                        #pass
            #if mass == 0.0:
                #return mass, cg, I

            #reference_point = cg / mass
        #else:
            ## reference_point = [0.,0.,0.] or user-defined array
            #pass

    mass = 0.
    cg = array([0., 0., 0.])
    I = array([0., 0., 0., 0., 0., 0., ])

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

        'CTRAX3', 'CTRAX6', 'CQUADX8', 'CQUADX4',
        'CPLSTN3', 'CPLSTN6', 'CPLSTN4', 'CPLSTN8',
    ]
    all_eids = np.array(list(model.elements.keys()), dtype='int32')
    all_eids.sort()

    all_mass_ids = np.array(list(model.masses.keys()), dtype='int32')
    all_mass_ids.sort()

    #def _increment_inertia0(centroid, reference_point, m, mass, cg, I):
        #"""helper method"""
        #(x, y, z) = centroid - reference_point
        #mass += m
        #cg += m * centroid
        #return mass

    def _increment_inertia(centroid, reference_point, m, mass, cg, I):
        """helper method"""
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

    def get_sub_eids(all_eids, eids):
        """supports limiting the element/mass ids"""
        eids = np.array(eids)
        ieids = np.searchsorted(all_eids, eids)
        eids2 = eids[all_eids[ieids] == eids]
        return eids2

    etypes_skipped = set([])
    for etype, eids in iteritems(model._type_to_id_map):
        if etype in no_mass:
            continue
        elif etype in ['CROD', 'CONROD']:
            eids2 = get_sub_eids(all_eids, eids)
            for eid in eids2:
                elem = model.elements[eid]
                n1, n2 = elem.node_ids
                length = norm(xyz[n2] - xyz[n1])
                centroid = (xyz[n1] + xyz[n2]) / 2.
                mpl = elem.MassPerLength()
                m = mpl * length
                mass = _increment_inertia(centroid, reference_point, m, mass, cg, I)
        elif etype == 'CTUBE':
            eids2 = get_sub_eids(all_eids, eids)
            for eid in eids2:
                elem = model.elements[eid]
                n1, n2 = elem.node_ids
                length = norm(xyz[n2] - xyz[n1])
                centroid = (xyz[n1] + xyz[n2]) / 2.
                mpl = elem.pid_ref.MassPerLength()
                m = mpl * length
                mass = _increment_inertia(centroid, reference_point, m, mass, cg, I)
        elif etype == 'CBAR':
            eids2 = get_sub_eids(all_eids, eids)
            for eid in eids2:
                elem = model.elements[eid]
                n1, n2 = elem.node_ids
                centroid = (xyz[n1] + xyz[n2]) / 2.
                length = norm(xyz[n2] - xyz[n1])
                mpl = elem.pid_ref.MassPerLength()
                m = mpl * length
                mass = _increment_inertia(centroid, reference_point, m, mass, cg, I)
        elif etype == 'CBEAM':
            eids2 = get_sub_eids(all_eids, eids)
            for eid in eids2:
                elem = model.elements[eid]
                prop = elem.pid_ref
                n1, n2 = elem.node_ids
                node1 = xyz[n1]
                node2 = xyz[n2]
                centroid = (node1 + node2) / 2.
                length = norm(node2 - node1)
                #cda = model.nodes[n1].cid_ref
                #cdb = model.nodes[n2].cid_ref

                is_failed, wa, wb, _ihat, jhat, khat = elem.get_axes(model)
                p1 = node1 + wa
                p2 = node2 + wb
                if prop.type == 'PBEAM':
                    rho = prop.Rho()
                    mass_per_lengths = []
                    nsm_per_lengths = []
                    for (area, nsm) in zip(prop.A, prop.nsm):
                        mass_per_lengths.append(area * rho)
                        nsm_per_lengths.append(nsm)
                    mass_per_length = integrate_positive_unit_line(prop.xxb, mass_per_lengths)
                    nsm_per_length = integrate_positive_unit_line(prop.xxb, nsm_per_lengths)
                    #nsm = np.mean(prop.nsm)
                    nsm = nsm_per_length * length
                    m = mass_per_length * length
                    nsm_n1 = (p1 + jhat * prop.m1a + khat * prop.m2a)
                    nsm_n2 = (p2 + jhat * prop.m1b + khat * prop.m2b)
                    nsm_centroid = (nsm_n1 + nsm_n2) / 2.
                    #if nsm != 0.:
                        #p1_nsm = p1 + prop.ma
                        #p2_nsm = p2 + prop.mb
                elif prop.type == 'PBEAML':
                    #mass_per_lengths = elem.get_mass_per_lengths()
                    mass_per_length = prop.MassPerLength() # includes simplified nsm
                    m = mass_per_length * length
                    nsm_centroid = np.zeros(3)
                    nsm = prop.nsm[0]  # TODO: simplified
                elif prop.type == 'PBCOMP':
                    mass_per_length = prop.MassPerLength()
                    m = mass_per_length * length
                    nsm = prop.nsm
                    nsm_n1 = (p1 + jhat * prop.m1 + khat * prop.m2)
                    nsm_n2 = (p2 + jhat * prop.m1 + khat * prop.m2)
                    nsm_centroid = (nsm_n1 + nsm_n2) / 2.
                else:
                    raise NotImplementedError(prop.type)

                #mpl = elem.pid_ref.MassPerLength()
                #m = mpl * length

                (x, y, z) = centroid - reference_point
                (xm, ym, zm) = nsm_centroid - reference_point
                x2 = x * x
                y2 = y * y
                z2 = z * z
                xm2 = xm * xm
                ym2 = ym * ym
                zm2 = zm * zm

                # Ixx, Iyy, Izz, Ixy, Ixz, Iyz
                I[0] += m * (y2 + z2) + nsm * (ym2 + zm2)
                I[1] += m * (x2 + z2) + nsm * (xm2 + zm2)
                I[2] += m * (x2 + y2) + nsm * (xm2 + ym2)
                I[3] += m * x * y + nsm * xm * ym
                I[4] += m * x * z + nsm * xm * zm
                I[5] += m * y * z + nsm * ym * zm
                mass += m + nsm
                cg += m * centroid + nsm * nsm_centroid

        elif etype in ['CTRIA3', 'CTRIA6', 'CTRIAR']:
            eids2 = get_sub_eids(all_eids, eids)
            for eid in eids2:
                elem = model.elements[eid]
                n1, n2, n3 = elem.node_ids[:3]
                prop = elem.pid_ref
                centroid = (xyz[n1] + xyz[n2] + xyz[n3]) / 3.
                area = 0.5 * norm(cross(xyz[n1] - xyz[n2], xyz[n1] - xyz[n3]))
                if prop.type == 'PSHELL':
                    tflag = elem.tflag
                    ti = prop.Thickness()
                    if tflag == 0:
                        # absolute
                        t1 = elem.T1 if elem.T1 else ti
                        t2 = elem.T2 if elem.T2 else ti
                        t3 = elem.T3 if elem.T3 else ti
                    elif tflag == 1:
                        # relative
                        t1 = elem.T1 * ti if elem.T1 else ti
                        t2 = elem.T2 * ti if elem.T2 else ti
                        t3 = elem.T3 * ti if elem.T3 else ti
                    else:
                        raise RuntimeError('tflag=%r' % tflag)
                    assert t1 + t2 + t3 > 0., 't1=%s t2=%s t3=%s' % (t1, t2, t3)
                    t = (t1 + t2 + t3) / 3.

                    # m/A = rho * A * t + nsm
                    #mass_per_area = elem.nsm + rho * elem.t

                    mpa = prop.nsm + prop.Rho() * t
                    #mpa = elem.pid_ref.MassPerArea()
                    m = mpa * area
                elif prop.type in ['PCOMP', 'PCOMPG']:
                    # PCOMP, PCOMPG
                    #rho_t = prop.get_rho_t()
                    #nsm = prop.nsm
                    #rho_t = [mat.Rho() * t for (mat, t) in zip(prop.mids_ref, prop.ts)]
                    #mpa = sum(rho_t) + nsm

                    # works for PCOMP
                    # F:\Program Files\Siemens\NXNastran\nxn10p1\nxn10p1\nast\tpl\cqr3compbuck.dat
                    mpa = prop.get_mass_per_area()
                elif prop.type == 'PLPLANE':
                    continue
                else:
                    raise NotImplementedError(prop.type)
                m = area * mpa
                mass = _increment_inertia(centroid, reference_point, m, mass, cg, I)
        elif etype in ['CQUAD4', 'CQUAD8', 'CQUADR']:
            eids2 = get_sub_eids(all_eids, eids)
            for eid in eids2:
                elem = model.elements[eid]
                n1, n2, n3, n4 = elem.node_ids[:4]
                prop = elem.pid_ref
                centroid = (xyz[n1] + xyz[n2] + xyz[n3] + xyz[n4]) / 4.
                area = 0.5 * norm(cross(xyz[n3] - xyz[n1], xyz[n4] - xyz[n2]))

                if prop.type == 'PSHELL':
                    tflag = elem.tflag
                    ti = prop.Thickness()
                    if tflag == 0:
                        # absolute
                        t1 = elem.T1 if elem.T1 else ti
                        t2 = elem.T2 if elem.T2 else ti
                        t3 = elem.T3 if elem.T3 else ti
                        t4 = elem.T4 if elem.T4 else ti
                    elif tflag == 1:
                        # relative
                        t1 = elem.T1 * ti if elem.T1 else ti
                        t2 = elem.T2 * ti if elem.T2 else ti
                        t3 = elem.T3 * ti if elem.T3 else ti
                        t4 = elem.T4 * ti if elem.T4 else ti
                    else:
                        raise RuntimeError('tflag=%r' % tflag)
                    assert t1 + t2 + t3 + t4 > 0., 't1=%s t2=%s t3=%s t4=%s' % (t1, t2, t3, t4)
                    t = (t1 + t2 + t3 + t4) / 4.

                    # m/A = rho * A * t + nsm
                    #mass_per_area = model.nsm + rho * model.t

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
                elif prop.type == 'PLPLANE':
                    continue
                    #raise NotImplementedError(prop.type)
                else:
                    raise NotImplementedError(prop.type)
                m = area * mpa
                mass = _increment_inertia(centroid, reference_point, m, mass, cg, I)
        elif etype == 'CQUAD':
            eids2 = get_sub_eids(all_eids, eids)
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
                mass = _increment_inertia(centroid, reference_point, m, mass, cg, I)

        elif etype == 'CSHEAR':
            eids2 = get_sub_eids(all_eids, eids)
            for eid in eids2:
                elem = model.elements[eid]
                n1, n2, n3, n4 = elem.node_ids
                prop = elem.pid_ref
                centroid = (xyz[n1] + xyz[n2] + xyz[n3] + xyz[n4]) / 4.
                area = 0.5 * norm(cross(xyz[n3] - xyz[n1], xyz[n4] - xyz[n2]))
                mpa = prop.MassPerArea()
                m = area * mpa
                mass = _increment_inertia(centroid, reference_point, m, mass, cg, I)
        elif etype in ['CONM1', 'CONM2', 'CMASS1', 'CMASS2', 'CMASS3', 'CMASS4']:
            eids2 = get_sub_eids(all_mass_ids, eids)
            for eid in eids2:
                elem = model.masses[eid]
                m = elem.Mass()
                centroid = elem.Centroid()
                mass = _increment_inertia(centroid, reference_point, m, mass, cg, I)
        elif etype == 'CTETRA':
            eids2 = get_sub_eids(all_eids, eids)
            for eid in eids2:
                elem = model.elements[eid]
                n1, n2, n3, n4 = elem.node_ids[:4]
                centroid = (xyz[n1] + xyz[n2] + xyz[n3] + xyz[n4]) / 4.
                #V = -dot(n1 - n4, cross(n2 - n4, n3 - n4)) / 6.
                volume = -dot(xyz[n1] - xyz[n4], cross(xyz[n2] - xyz[n4], xyz[n3] - xyz[n4])) / 6.
                m = elem.Rho() * volume
                mass = _increment_inertia(centroid, reference_point, m, mass, cg, I)
        elif etype == 'CPYRAM':
            eids2 = get_sub_eids(all_eids, eids)
            for eid in eids2:
                elem = model.elements[eid]
                n1, n2, n3, n4, n5 = elem.node_ids[:5]
                centroid1 = (xyz[n1] + xyz[n2] + xyz[n3] + xyz[n4]) / 4.
                area1 = 0.5 * norm(cross(xyz[n3]-xyz[n1], xyz[n4]-xyz[n2]))
                centroid5 = xyz[n5]
                centroid = (centroid1 + centroid5) / 2.
                volume = area1 / 3. * norm(centroid1 - centroid5)
                m = elem.Rho() * volume
                mass = _increment_inertia(centroid, reference_point, m, mass, cg, I)
        elif etype == 'CPENTA':
            eids2 = get_sub_eids(all_eids, eids)
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
                mass = _increment_inertia(centroid, reference_point, m, mass, cg, I)

        elif etype == 'CHEXA':
            eids2 = get_sub_eids(all_eids, eids)
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
                mass = _increment_inertia(centroid, reference_point, m, mass, cg, I)

        elif etype in no_mass:
            continue
        elif etype == 'CBEND':
            model.log.info('elem.type=%s doesnt have mass' % etype)
            continue
        elif etype.startswith('C'):
            eids2 = get_sub_eids(all_eids, eids)
            for eid in eids2:
                elem = model.elements[eid]

                #if elem.pid_ref.type in ['PPLANE']:
                try:
                    m = elem.Mass()
                except:
                    model.log.error('etype = %r' % etype)
                    print(elem)
                    print(elem.pid_ref)
                    raise
                centroid = elem.Centroid()
                if m > 0.0:
                    model.log.info('elem.type=%s is not supported in new '
                                   'mass properties method' % elem.type)
                    mass = _increment_inertia(centroid, reference_point, m, mass, cg, I)
                elif etype not in etypes_skipped:
                    model.log.info('elem.type=%s doesnt have mass' % elem.type)
                    etypes_skipped.add(etype)

    if mass:
        cg /= mass
    mass, cg, I = _apply_mass_symmetry(model, sym_axis, scale, mass, cg, I)
    # Ixx, Iyy, Izz, Ixy, Ixz, Iyz = I
    return mass, cg, I

def _apply_mass_symmetry(model, sym_axis, scale, mass, cg, I):
    """
    Scales the mass & moement of inertia based on the symmetry axes
    and the PARAM WTMASS card
    """
    if isinstance(sym_axis, string_types):
        sym_axis = [sym_axis]
    elif isinstance(sym_axis, (list, tuple)):
        # basically overwrite the existing values on the AERO/AEROS card
        pass
    else:
        sym_axis = []

        # The symmetry flags on the AERO/AEROS must be the same, so
        # it doesn't matter which we one pick.  However, they might
        # not both be defined.
        if model.aero is not None:
            if model.aero.is_symmetric_xy():
                sym_axis.append('xy')
            if model.aero.is_symmetric_xz():
                sym_axis.append('xz')
            if model.aero.is_anti_symmetric_xy():
                sym_axis.append('xy')
                #raise NotImplementedError('%s load is anti-symmetric about the XY plane'
                                          #% str(aero))
            if model.aero.is_anti_symmetric_xz():
                #raise NotImplementedError('%s load is anti-symmetric about the XZ plane'
                                          #% str(aero))
                sym_axis.append('xz')

        if model.aeros is not None:
            if model.aeros.is_symmetric_xy():
                sym_axis.append('xy')
            if model.aeros.is_symmetric_xz():
                sym_axis.append('xz')
            if model.aeros.is_anti_symmetric_xy():
                sym_axis.append('xy')
                #raise NotImplementedError('%s load is anti-symmetric about the XY plane'
                                          #% str(aeros))
            if model.aeros.is_anti_symmetric_xz():
                sym_axis.append('xz')
                #raise NotImplementedError('%s load is anti-symmetric about the XZ plane' %
                                          #str(aeros))

    sym_axis = list(set(sym_axis))
    short_sym_axis = [sym_axisi.lower() for sym_axisi in sym_axis]
    is_no = 'no' in short_sym_axis
    if is_no and len(short_sym_axis) > 1:
        raise RuntimeError('no can only be used by itself; sym_axis=%s' % (str(sym_axis)))
    for sym_axisi in sym_axis:
        if sym_axisi.lower() not in ['no', 'xy', 'yz', 'xz']:
            msg = 'sym_axis=%r is invalid; sym_axisi=%r; allowed=[no, xy, yz, xz]' % (
                sym_axis, sym_axisi)
            raise RuntimeError(msg)

    if sym_axis:
        # either we figured sym_axis out from the AERO cards or the user told us
        model.log.debug('Mass/MOI sym_axis = %r' % sym_axis)

        if 'xz' in sym_axis:
            # y intertias are 0
            cg[1] = 0.0
            mass *= 2.0
            I[0] *= 2.0
            I[1] *= 2.0
            I[2] *= 2.0
            I[3] *= 0.0  # Ixy
            I[4] *= 2.0  # Ixz; no y
            I[5] *= 0.0  # Iyz

        if 'xy' in sym_axis:
            # z intertias are 0
            cg[2] = 0.0
            mass *= 2.0
            I[0] *= 2.0
            I[1] *= 2.0
            I[2] *= 2.0
            I[3] *= 2.0  # Ixy; no z
            I[4] *= 0.0  # Ixz
            I[5] *= 0.0  # Iyz

        if 'yz' in sym_axis:
            # x intertias are 0
            cg[0] = 0.0
            mass *= 2.0
            I[0] *= 2.0
            I[1] *= 2.0
            I[2] *= 2.0
            I[3] *= 0.0  # Ixy
            I[4] *= 0.0  # Ixz
            I[5] *= 2.0  # Iyz; no x

    if scale is None and 'WTMASS' in model.params:
        param = model.params['WTMASS']
        #assert isinstance(param, PARAM), 'param=%s' % param
        scale = param.values[0]
        if scale != 1.0:
            model.log.info('WTMASS scale = %r' % scale)
    elif scale is None:
        scale = 1.0
    mass *= scale
    I *= scale
    return (mass, cg, I)
