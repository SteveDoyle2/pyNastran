"""
Defines:
  - sum_forces_moments
      find the net force/moment on the model
  - sum_forces_moments_elements
      find the net force/moment on the model for a subset of elements
"""
from __future__ import print_function
import numpy as np
from numpy import array, cross, allclose, mean
from numpy.linalg import norm  # type: ignore
from pyNastran.utils.numpy_utils import integer_types
from pyNastran.bdf.cards.loads.static_loads import update_pload4_vector_for_surf

def sum_forces_moments(model, p0, loadcase_id, include_grav=False, xyz_cid0=None):
    """
    Sums applied forces & moments about a reference point p0 for all
    load cases.

    Considers:
      - FORCE, FORCE1, FORCE2
      - MOMENT, MOMENT1, MOMENT2
      - PLOAD, PLOAD2, PLOAD4
      - LOAD

    Parameters
    ----------
    model : BDF()
        a BDF object
    p0 : NUMPY.NDARRAY shape=(3,) or integer (node ID)
        the reference point
    loadcase_id : int
        the LOAD=ID to analyze
    include_grav : bool; default=False
        includes gravity in the summation (not supported)
    xyz_cid0 : None / Dict[int] = (3, ) ndarray
        the nodes in the global coordinate system

    Returns
    -------
    forces : NUMPY.NDARRAY shape=(3,)
        the forces
    moments : NUMPY.NDARRAY shape=(3,)
        the moments

    .. warning:: not full validated
    .. todo:: It's super slow for cid != 0.   We can speed this up a lot
             if we calculate the normal, area, centroid based on
             precomputed node locations.

    Pressure acts in the normal direction per model/real/loads.bdf and loads.f06

    """
    if not isinstance(loadcase_id, integer_types):
        raise RuntimeError('loadcase_id must be an integer; loadcase_id=%r' % loadcase_id)

    p = _get_load_summation_point(model, p0, cid=0)
    loads, scale_factors, unused_is_grav = model.get_reduced_loads(
        loadcase_id, skip_scale_factor0=True)

    F = array([0., 0., 0.])
    M = array([0., 0., 0.])
    xyz = _get_xyz_cid0_dict(model, xyz_cid0=xyz_cid0)

    unsupported_types = set([])
    for load, scale in zip(loads, scale_factors):
        #if load.type not in ['FORCE1']:
            #continue
        if load.type == 'FORCE':
            if load.Cid() != 0:
                cp_ref = load.cid_ref
                #from pyNastran.bdf.bdf import CORD2R
                #cp_ref = CORD2R()
                f = load.mag * cp_ref.transform_vector_to_global(load.xyz) * scale
            else:
                f = load.mag * load.xyz * scale

            node = model.Node(load.node_id)
            r = xyz[node.nid] - p
            m = cross(r, f)
            F += f
            M += m
        elif load.type == 'FORCE1':
            f = load.mag * load.xyz * scale
            node = model.Node(load.node_id)
            r = xyz[node.nid] - p
            m = cross(r, f)
            F += f
            M += m
        elif load.type == 'FORCE2':
            f = load.mag * load.xyz * scale
            node = model.Node(load.node_id)
            r = xyz[node.nid] - p
            m = cross(r, f)
            F += f
            M += m
        elif load.type == 'MOMENT':
            if load.Cid() != 0:
                cp = load.cid_ref
                #from pyNastran.bdf.bdf import CORD2R
                #cp = CORD2R()
                m = load.mag * cp.transform_vector_to_global(load.xyz) * scale
            else:
                m = load.mag * load.xyz * scale
            M += m
        elif load.type == 'MOMENT1':
            m = load.mag * load.xyz * scale
            M += m
        elif load.type == 'MOMENT2':
            m = load.mag * load.xyz * scale
            M += m

        elif load.type == 'PLOAD':
            nodes = load.node_ids
            nnodes = len(nodes)
            if nnodes == 3:
                n1, n2, n3 = xyz[nodes[0]], xyz[nodes[1]], xyz[nodes[2]]
                axb = cross(n1 - n2, n1 - n3)
                centroid = (n1 + n2 + n3) / 3.
            elif nnodes == 4:
                n1, n2, n3, n4 = xyz[nodes[0]], xyz[nodes[1]], xyz[nodes[2]], xyz[nodes[3]]
                axb = cross(n1 - n3, n2 - n4)
                centroid = (n1 + n2 + n3 + n4) / 4.
            else:
                msg = 'invalid number of nodes on PLOAD card; nodes=%s' % str(nodes)
                raise RuntimeError(msg)

            area, normal = _get_area_normal(axb, nodes, xyz)
            r = centroid - p
            f = load.pressure * area * normal * scale
            m = cross(r, f)

            F += f
            M += m

        elif load.type == 'PLOAD1':
            _pload1_total(model, loadcase_id, load, scale, xyz, F, M, p)

        elif load.type == 'PLOAD2':
            pressure = load.pressure * scale
            for eid in load.element_ids:
                elem = model.elements[eid]
                if elem.type in ['CTRIA3', 'CQUAD4', 'CSHEAR', 'CQUADR', 'CTRIAR']:
                    n = elem.Normal()
                    area = elem.Area()
                    f = pressure * n * area
                    r = elem.Centroid() - p
                    m = cross(r, f)
                    F += f
                    M += m
                else:
                    model.log.warning('case=%s etype=%r loadtype=%r not supported' % (
                        loadcase_id, elem.type, load.type))
        elif load.type == 'PLOAD4':
            _pload4_total(loadcase_id, load, scale, xyz, F, M, p)

        elif load.type == 'GRAV':
            if include_grav:  # this will be super slow
                g = load.GravityVector() * scale
                for eid, elem in model.elements.items():
                    centroid = elem.Centroid()
                    mass = elem.Mass()
                    r = centroid - p
                    f = mass * g
                    m = cross(r, f)
                    F += f
                    M += m
        else:
            # we collect them so we only get one print
            unsupported_types.add(load.type)

    for load_type in unsupported_types:
        model.log.debug('case=%s loadtype=%r not supported' % (loadcase_id, load_type))
    return (F, M)

def _pload1_total(model, loadcase_id, load, scale, xyz, F, M, p):
    """helper method for ``sum_forces_moments``"""
    elem = load.eid_ref
    if elem.type in ['CBAR', 'CBEAM']:
        _pload1_bar_beam(model, loadcase_id, load, elem, scale, xyz, F, M, p)
    elif elem.type == 'CBEND':
        model.log.warning('case=%s etype=%r loadtype=%r not supported' % (
            loadcase_id, elem.type, load.type))
    else:
        raise RuntimeError('element.type=%r is not a CBAR, CBEAM, or CBEND' % elem.type)
    return

def _pload1_elements(model, loadcase_id, load, scale, eids, xyz, F, M, p):
    """helper method for ``sum_forces_moments_elements``"""
    #elem = model.elements[load.eid]
    elem = load.eid_ref
    if elem.eid not in eids:
        return
    if elem.type in ['CBAR', 'CBEAM']:
        _pload1_bar_beam(model, loadcase_id, load, elem, scale, xyz, F, M, p)
    elif elem.type == 'CBEND':
        model.log.warning('case=%s etype=%r loadtype=%r not supported' % (
            loadcase_id, elem.type, load.type))
    else:
        raise RuntimeError('element.type=%r is not a CBAR, CBEAM, or CBEND' % elem.type)
    return

def _pload1_bar_beam(model, unused_loadcase_id, load, elem, scale, xyz, F, M, p):
    """
    helper method for ``sum_forces_moments`` and ``sum_forces_moments_elements``
    """
    p1 = load.p1 * scale
    p2 = load.p2 * scale

    nodes = elem.node_ids
    n1, n2 = xyz[nodes[0]], xyz[nodes[1]]
    n1 += elem.wa
    n2 += elem.wb

    bar_vector = n2 - n1
    L = norm(bar_vector)
    try:
        Ldir = bar_vector / L
    except:
        msg = 'Length=0.0; nid1=%s nid2=%s\n' % (nodes[0], nodes[1])
        msg += '%s%s' % (str(elem.nodes[0]), str(elem.nodes[1]))
        raise FloatingPointError(msg)

    if load.scale == 'FR':  # x1, x2 are fractional lengths
        x1 = load.x1
        x2 = load.x2
        #compute_fx = False
    elif load.scale == 'LE': # x1, x2 are actual lengths
        x1 = load.x1 / L
        x2 = load.x2 / L
    elif load.scale == 'LEPR':
        model.log.warning('PLOAD1: LEPR continue')
        return
        #msg = 'scale=%r is not supported.  Use "FR", "LE".' % load.scale
        #raise NotImplementedError(msg)
    elif load.scale == 'FRPR':
        model.log.warning('PLOAD1: FRPR continue')
        return
        #msg = 'scale=%r is not supported.  Use "FR", "LE".' % load.scale
        #raise NotImplementedError(msg)
    else:
        msg = 'PLOAD1: scale=%r is not supported.  Use "FR", "LE".' % load.scale
        raise NotImplementedError(msg)

    # FY - force in basic coordinate system
    # FR - fractional;
    assert x1 <= x2, 'x1=%s x2=%s' % (x1, x2)
    if x1 != x2:
        # continue
        if not load.type in ['FX', 'FY', 'FZ']:
            model.log.warning('PLOAD1 x1 != x2 continue; x1=%s x2=%s; scale=%r\n%s%s'% (
                x1, x2, load.scale, str(elem), str(load)))
            return
        model.log.warning('check this...PLOAD1 x1 != x2; x1=%s x2=%s; scale=%r\n%s%s'% (
            x1, x2, load.scale, str(elem), str(load)))

        # y = (y2-y1)/(x2-x1)*(x-x1) + y1
        # y = (y2-y1) * (x-x1)/(x2-x1) + y1
        # y = y2*(x-x1)/(x2-x1) + y1*(1-(x-x1)/(x2-x1))
        # y = y2 * r + y1 * (1-r)
        # r = (x-x1)/(x2-x1)
        #
        # y = y2 * r + y1 - y1 * r
        # yi = y2 * ri + y1 * x + y1 * ri
        # yi = y2 * ri + y1 * (x2-x1) + y1 * ri
        #
        # ri = integral(r)
        # ri = 1/(x2-x1) * (0.5) * (x1-x2)**2
        #
        # yi = integral(y)
        # yi = y2 * ri + y1 * (x2-x1) + y1 * ri
        # ri = 1./(x2-x1) * (0.5) * (x1-x2)**2
        # y1 = p1
        # y2 = p2
        # yi = y2 * ri + y1 * (x2-x1) + y1 * ri
        # F = yi
        if allclose(p1, -p2):
            Ftotal = p1
            x = (x1 + x2) / 2.
        else:
            Ftotal = L * (x2-x1) * (p1 + p2)/2.
            Mx = L * p1 * (x2-x1)/2. + L * (p2-p1) * (2./3. * x2 + 1./3. * x1)
            x = Mx / Ftotal
        model.log.info('L=%s x1=%s x2=%s p1/L=%s p2/L=%s Ftotal=%s Mtotal=%s x=%s' % (
            L, x1, x2, p1, p2, Ftotal, Mx, x))

        unused_i = Ldir
        if load.Type in ['FX', 'FY', 'FZ']:
            r = (1. - x) * n1 + x * n2
            #print('r=%s n1=%s n2=%s' % (r, n1, n2))
            if load.Type == 'FX':
                force_dir = array([1., 0., 0.])
            elif load.Type == 'FY':
                force_dir = array([0., 1., 0.])
            elif load.Type == 'FZ':
                force_dir = array([0., 0., 1.])
            else:
                raise NotImplementedError('Type=%r is not supported.  '
                                          'Use "FX", "FY", "FZ".' % load.Type)

        Fi = Ftotal * force_dir
        Mi = cross(r - p, force_dir * Ftotal)
        F += Fi
        M += Mi
        model.log.info('Fi=%s Mi=%s x=%s' % (Fi, Mi, x))
    else:
        _bar_eq_pload1(load, elem, xyz, Ldir,
                       n1, n2,
                       x1, x2,
                       p1, p2,
                       F, M, p)
    return

def sum_forces_moments_elements(model, p0, loadcase_id, eids, nids,
                                include_grav=False, xyz_cid0=None):
    """
    Sum the forces/moments based on a list of nodes and elements.

    Parameters
    ----------
    model : BDF()
        a BDF object
    eids : List[int]
        the list of elements to include (e.g. the loads due to a PLOAD4)
    nids : List[int]
        the list of nodes to include (e.g. the loads due to a FORCE card)
    p0 : int; (3,) ndarray
       the point to sum moments about
       type = int
           sum moments about the specified grid point
       type = (3, ) ndarray/list (e.g. [10., 20., 30]):
           the x, y, z location in the global frame
    loadcase_id : int
        the LOAD=ID to analyze
    include_grav : bool; default=False
        includes gravity in the summation (not supported)
    xyz_cid0 : None / Dict[int] = (3, ) ndarray
        the nodes in the global coordinate system

    Returns
    -------
    forces : NUMPY.NDARRAY shape=(3,)
        the forces
    moments : NUMPY.NDARRAY shape=(3,)
        the moments

    Nodal Types  : FORCE, FORCE1, FORCE2,
                   MOMENT, MOMENT1, MOMENT2,
                   PLOAD
    Element Types: PLOAD1, PLOAD2, PLOAD4, GRAV

    If you have a CQUAD4 (eid=3) with a PLOAD4 (sid=3) and a FORCE
    card (nid=5) acting on it, you can incldue the PLOAD4, but
    not the FORCE card by using:

    For just pressure:

    .. code-block:: python

      eids = [3]
      nids = []

    For just force:

    .. code-block:: python

      eids = []
      nids = [5]

    or both:

    .. code-block:: python

      eids = [3]
      nids = [5]

    .. note:: If you split the model into sections and sum the loads
              on each section, you may not get the same result as
              if you summed the loads on the total model.  This is
              due to the fact that nodal loads on the boundary are
              double/triple/etc. counted depending on how many breaks
              you have.

    .. todo:: not done...

    """
    if not isinstance(loadcase_id, integer_types):
        raise RuntimeError('loadcase_id must be an integer; loadcase_id=%r' % loadcase_id)
    p = _get_load_summation_point(model, p0, cid=0)

    if eids is None:
        eids = list(model.element_ids)
    if nids is None:
        nids = list(model.node_ids)

    #for (key, load_case) in model.loads.items():
        #if key != loadcase_id:
            #continue

    loads, scale_factors, unused_is_grav = model.get_reduced_loads(
        loadcase_id, skip_scale_factor0=True)

    F = array([0., 0., 0.])
    M = array([0., 0., 0.])

    xyz = _get_xyz_cid0_dict(model, xyz_cid0)

    unsupported_types = set([])
    shell_elements = [
        'CTRIA3', 'CQUAD4', 'CTRIAR', 'CQUADR',
        'CTRIA6', 'CQUAD8', 'CQUAD', 'CSHEAR']
    for load, scale in zip(loads, scale_factors):
        #if load.type not in ['FORCE1']:
            #continue
        #print(load.type)
        if load.type == 'FORCE':
            if load.node_id not in nids:
                continue
            if load.Cid() != 0:
                cp_ref = load.cid_ref
                #from pyNastran.bdf.bdf import CORD2R
                #cp = CORD2R()
                f = load.mag * cp_ref.transform_vector_to_global(load.xyz) * scale
            else:
                f = load.mag * load.xyz * scale

            node = model.Node(load.node_id)
            r = xyz[node.nid] - p
            m = cross(r, f)
            F += f
            M += m

        elif load.type == 'FORCE1':
            not_found_nid = False
            for nid in load.node_ids:
                if nid not in nids:
                    not_found_nid = True
                    break
            if not_found_nid:
                continue

            f = load.mag * load.xyz * scale
            node = model.Node(load.node_id)
            r = xyz[node.nid] - p
            m = cross(r, f)
            F += f
            M += m
        elif load.type == 'FORCE2':
            not_found_nid = False
            for nid in load.node_ids:
                if nid not in nids:
                    not_found_nid = True
                    break
            if not_found_nid:
                continue

            f = load.mag * load.xyz * scale
            node = model.Node(load.node_id)
            r = xyz[node.nid] - p
            m = cross(r, f)
            F += f
            M += m
        elif load.type == 'MOMENT':
            not_found_nid = False
            for nid in load.node_ids:
                if nid not in nids:
                    not_found_nid = True
                    break
            if not_found_nid:
                continue

            if load.Cid() != 0:
                cp_ref = load.cid_ref
                m = cp_ref.transform_vector_to_global(load.xyz)
            else:
                m = load.xyz
            M += load.mag * m * scale
        elif load.type == 'MOMENT1':
            not_found_nid = False
            for nid in load.node_ids:
                if nid not in nids:
                    not_found_nid = True
                    break
            if not_found_nid:
                continue
            m = load.mag * load.xyz * scale
            M += m
        elif load.type == 'MOMENT2':
            not_found_nid = False
            for nid in load.node_ids:
                if nid not in nids:
                    not_found_nid = True
                    break
            if not_found_nid:
                continue
            m = load.mag * load.xyz * scale
            M += m

        elif load.type == 'PLOAD':
            nodes = load.node_ids
            nnodes = len(nodes)
            nodesi = 0
            if nnodes == 3:
                n1, n2, n3 = xyz[nodes[0]], xyz[nodes[1]], xyz[nodes[2]]
                axb = cross(n1 - n2, n1 - n3)
                centroid = (n1 + n2 + n3) / 3.

            elif nnodes == 4:
                n1, n2, n3, n4 = xyz[nodes[0]], xyz[nodes[1]], xyz[nodes[2]], xyz[nodes[3]]
                axb = cross(n1 - n3, n2 - n4)
                centroid = (n1 + n2 + n3 + n4) / 4.
                if nodes[3] in nids:
                    nodesi += 1
            else:
                raise RuntimeError('invalid number of nodes on PLOAD card; '
                                   'nodes=%s' % str(nodes))
            if nodes[0] in nids:
                nodesi += 1
            if nodes[1] in nids:
                nodesi += 1
            if nodes[2] in nids:
                nodesi += 1

            area, normal = _get_area_normal(axb, nodes, xyz)
            r = centroid - p
            f = load.pressure * area * normal * scale
            m = cross(r, f)

            node_scale = nodesi / float(nnodes)
            F += f * node_scale
            M += m * node_scale

        elif load.type == 'PLOAD1':
            _pload1_elements(model, loadcase_id, load, scale, eids, xyz, F, M, p)

        elif load.type == 'PLOAD2':
            pressure = load.pressure * scale
            for eid in load.element_ids:
                if eid not in eids:
                    continue
                elem = model.elements[eid]
                if elem.type in shell_elements:
                    normal = elem.Normal()
                    area = elem.Area()
                    f = pressure * normal * area
                    r = elem.Centroid() - p
                    m = cross(r, f)
                    F += f
                    M += m
                else:
                    #model.log.warning('case=%s etype=%r loadtype=%r not supported' % (
                        #loadcase_id, elem.type, load.type))
                    raise NotImplementedError('case=%s etype=%r loadtype=%r not supported' % (
                        loadcase_id, elem.type, load.type))
        elif load.type == 'PLOAD4':
            _pload4_elements(loadcase_id, load, scale, eids, xyz, F, M, p)

        elif load.type == 'GRAV':
            if include_grav:  # this will be super slow
                g = load.GravityVector() * scale
                for eid, elem in model.elements.items():
                    if eid not in eids:
                        continue
                    centroid = elem.Centroid()
                    mass = elem.Mass()
                    r = centroid - p
                    f = mass * g
                    m = cross(r, f)
                    F += f
                    M += m
        else:
            # we collect them so we only get one print
            unsupported_types.add(load.type)

    for loadtype in unsupported_types:
        model.log.warning('case=%s loadtype=%r not supported' % (loadcase_id, loadtype))
    #model.log.info("case=%s F=%s M=%s\n" % (loadcase_id, F, M))
    return (F, M)


def _bar_eq_pload1(load, elem, xyz, Ldir,
                   n1, n2,
                   x1, x2,
                   p1, unused_p2,
                   F, M, p):
    """helper for ``_elements_pload1`` and ``_elementi_pload1``"""
    v = elem.get_orientation_vector(xyz)
    i = Ldir
    ki = cross(i, v)
    k = ki / norm(ki)
    j = cross(k, i)

    if load.Type in ['FX', 'FY', 'FZ']:
        #deltaL = n2 - n1
        r = (1 - x1) * n1 + x1 * n2
        #print('    r =', r)
        #print('    n1 =', n1)
        #print('    n2 =', n2)
        #print('    x1 =', x1)
        #print('    1-x1 =', 1-x1)
        #print('    deltaL =', deltaL)
        if load.Type == 'FX' and x1 == x2:
            force_dir = array([1., 0., 0.])
        elif load.Type == 'FY' and x1 == x2:
            force_dir = array([0., 1., 0.])
        elif load.Type == 'FZ' and x1 == x2:
            force_dir = array([0., 0., 1.])
        F += p1 * force_dir
        M += cross(r - p, F)
    elif load.Type in ['MX', 'MY', 'MZ']:
        if load.Type == 'MX' and x1 == x2:
            moment_dir = array([1., 0., 0.])
        elif load.Type == 'MY' and x1 == x2:
            moment_dir = array([0., 1., 0.])
        elif load.Type == 'MZ' and x1 == x2:
            moment_dir = array([0., 0., 1.])
        M += p1 * moment_dir
    elif load.Type in ['FXE', 'FYE', 'FZE']:
        r = (1 - x1) * n1 + x1 * n2
        #print('\n    r =', r)
        #print('    n1 =', n1)
        #print('    n2 =', n2)
        #print('    x1 =', x1)
        #print('    1-x1 =', 1-x1)
        #print('    i    =', i)
        #print('    j    =', j)
        #print('    k    =', k)
        if load.Type == 'FXE' and x1 == x2:
            force_dir = i
        elif load.Type == 'FYE' and x1 == x2:
            force_dir = j
        elif load.Type == 'FZE' and x1 == x2:
            force_dir = k
        #print('    force_dir =', force_dir, load.Type)
        try:
            F += p1 * force_dir
        except FloatingPointError:
            msg = 'eid = %s\n' % elem.eid
            msg += 'i = %s\n' % Ldir
            msg += 'force_dir = %s\n' % force_dir
            msg += 'load = \n%s' % str(load)
            raise FloatingPointError(msg)
        M += cross(r - p, F)
        del force_dir

    elif load.Type in ['MXE', 'MYE', 'MZE']:
        if load.Type == 'MXE' and x1 == x2:
            moment_dir = i
        elif load.Type == 'MYE' and x1 == x2:
            moment_dir = j
        elif load.Type == 'MZE' and x1 == x2:
            moment_dir = k
        try:
            M += p1 * moment_dir
        except FloatingPointError:
            msg = 'eid = %s\n' % elem.eid
            msg += 'moment_dir = %s\n' % moment_dir
            msg += 'load = \n%s' % str(load)
            raise FloatingPointError(msg)
        del moment_dir
    else:
        raise NotImplementedError(
            'Type=%r is not supported.\n'
            'Use [FX, FXE, FY, FYE, FZ, FZE,\n'
            '     MX, MXE, MY, MYE, MZ, MZE]' % load.Type)
    return


def _pload4_total(loadcase_id, load, scale, xyz, F, M, p):
    """helper method for ``sum_forces_moments``"""
    assert load.line_load_dir == 'NORM', 'line_load_dir = %s' % (load.line_load_dir)
    for elem in load.eids_ref:
        fi, mi = _pload4_helper(loadcase_id, load, scale, elem, xyz, p)
        F += fi
        M += mi
    return

def _pload4_elements(loadcase_id, load, scale, eids, xyz, F, M, p):
    """helper method for ``sum_forces_moments_elements``"""
    assert load.line_load_dir == 'NORM', 'line_load_dir = %s' % (load.line_load_dir)
    for elem in load.eids_ref:
        eid = elem.eid
        if eid not in eids:
            continue
        fi, mi = _pload4_helper(loadcase_id, load, scale, elem, xyz, p)
        F += fi
        M += mi
    return

def _get_pload4_area_centroid_normal_nface(loadcase_id, load, elem, xyz):
    """gets the nodes, area, centroid, normal, and nface"""
    etype = elem.type
    if etype in ['CTRIA3', 'CTRIA6', 'CTRIAR',]:
        # triangles
        nodes = elem.node_ids
        n1, n2, n3 = xyz[nodes[0]], xyz[nodes[1]], xyz[nodes[2]]
        axb = cross(n1 - n2, n1 - n3)
        area, normal = _get_area_normal(axb, nodes, xyz)
        centroid = (n1 + n2 + n3) / 3.
        nface = 3
    elif etype in ['CQUAD4', 'CQUAD8', 'CQUAD', 'CQUADR', 'CSHEAR']:
        # quads
        nodes = elem.node_ids
        n1, n2, n3, n4 = xyz[nodes[0]], xyz[nodes[1]], xyz[nodes[2]], xyz[nodes[3]]
        axb = cross(n1 - n3, n2 - n4)
        area, normal = _get_area_normal(axb, nodes, xyz)
        centroid = (n1 + n2 + n3 + n4) / 4.
        nface = 4

    elif etype == 'CTETRA':
        nodes = None
        face_acn = elem.get_face_area_centroid_normal(load.g1_ref.nid, load.g34_ref.nid)
        unused_face, area, centroid, normal = face_acn
        nface = 3

    elif etype == 'CHEXA':
        nodes = None
        face_acn = elem.get_face_area_centroid_normal(load.g34_ref.nid, load.g1_ref.nid)
        # TODO: backwards?
        #face_acn= elem.get_face_area_centroid_normal(load.g1_ref.nid, load.g34_ref.nid)
        unused_face, area, centroid, normal = face_acn
        nface = 4

    elif etype == 'CPENTA':
        nodes = None
        g1 = load.g1_ref.nid
        if load.g34 is None:
            face_acn = elem.get_face_area_centroid_normal(g1)
            nface = 3
        else:
            face_acn = elem.get_face_area_centroid_normal(g1, load.g34_ref.nid)
            nface = 4
        unused_face, area, centroid, normal = face_acn
    else:
        eid = elem.eid
        msg = 'PLOAd4: case=%s eid=%s etype=%r loadtype=%r not supported\n%s%s' % (
            loadcase_id, eid, etype, load.type, str(load), str(elem))
        raise NotImplementedError(msg)
    return nodes, area, centroid, normal, nface

def _pload4_helper(loadcase_id, load, scale, elem, xyz, p):
    """gets the contribution for a single PLOAD4 element"""
    #eid = elem.eid
    nodes, area, centroid, normal, nface = _get_pload4_area_centroid_normal_nface(
        loadcase_id, load, elem, xyz)

    pressures = load.pressures[:nface]
    assert len(pressures) == nface
    if load.surf_or_line == 'SURF':
        pressure = _mean_pressure_on_pload4(pressures, load, elem)
        cid = load.Cid()
        normal = update_pload4_vector_for_surf(load, normal, cid)

        r = centroid - p
        fi = pressure * area * normal * scale
        #load.cid_ref.transform_to_global()
        mi = cross(r, fi)

    elif load.surf_or_line == 'LINE':
        fi, mi = _pload4_helper_line(load, elem, scale, pressures, nodes, xyz, p)
    else:  # pragma: no cover
        msg = 'surf_or_line=%r on PLOAD4 is not supported\n%s' % (
            load.surf_or_line, str(load))
        raise NotImplementedError(msg)
    return fi, mi

def _get_area_normal(axb, nodes, xyz):
    """gets the area/normal vector"""
    nunit = norm(axb)
    area = 0.5 * nunit
    try:
        normal = axb / nunit
    except FloatingPointError:
        msg = ''
        for i, nid in enumerate(nodes):
            msg += 'nid%i=%i node=%s\n' % (i+1, nid, xyz[nodes[i]])
        msg += 'a x b = %s\n' % axb
        msg += 'nunit = %s\n' % nunit
        raise FloatingPointError(msg)
    return area, normal

def _mean_pressure_on_pload4(pressures, load, unused_elem):
    """gets the mean pressure"""
    if min(pressures) != max(pressures):
        pressure = mean(pressures)
        #msg = ('%s%s\npressure.min=%s != pressure.max=%s using average of %%s; '
               #'load=%s eid=%%s'  % (str(load), str(elem), min(pressures),
                                     #max(pressures), load.sid))
        #print(msg % (pressure, eid))
    else:
        pressure = load.pressures[0]
    return pressure

def _get_xyz_cid0_dict(model, xyz_cid0=None):
    """
    helper method

    Parameters
    ----------
    model : BDF()
        a BDF object
    xyz_cid0 : None / Dict[int] = (3, ) ndarray
        the nodes in the global coordinate system

    """
    if xyz_cid0 is None:
        xyz = {}
        for nid, node in model.nodes.items():
            xyz[nid] = node.get_position()
    else:
        xyz = xyz_cid0
    return xyz

def _get_load_summation_point(model, p0, cid=0):
    """
    helper method

    Parameters
    ----------
    model : BDF()
        a BDF object
    p0 : NUMPY.NDARRAY shape=(3,) or integer (node ID)
        the reference point

    """
    if isinstance(p0, integer_types):
        if cid == 0:
            p = model.nodes[p0].get_position()
        else:
            p = model.nodes[p0].get_position_wrt(model, cid)
    else:
        p = array(p0)
    return p

def _pload4_helper_line(load, elem, scale, pressures, nodes, xyz, p):
    # this is pressure per unit length?
    # edge_length * thickness I assume?
    fi = np.zeros(3)
    mi = np.zeros(3)
    edges = []
    if len(pressures) == 4:
        p1, p2, p3, p4 = pressures
        if p1 or p2:
            edges.append((0, 1))
        if p2 or p3:
            edges.append((1, 2))
        if p3 or p4:
            edges.append((2, 3))
        if p4 or p1:
            edges.append((3, 1))
    elif len(pressures) == 3:
        p1, p2, p3, p4 = pressures
        if p1 or p2:
            edges.append((0, 1))
        if p2 or p3:
            edges.append((1, 2))
        if p3 or p1:
            edges.append((2, 1))
    else:  # pragma: no cover
        raise NotImplementedError(pressures)

    thickness = elem.Thickness()
    load_dir = load.nvector / np.linalg.norm(load.nvector)
    for edge in edges:
        inode1, inode2 = edge
        ixyz1 = nodes[inode1]
        ixyz2 = nodes[inode2]
        xyz1 = xyz[ixyz1]
        xyz2 = xyz[ixyz2]
        p1 = pressures[inode1]
        p2 = pressures[inode2]
        area_edge = thickness * np.linalg.norm(xyz2 - xyz1)
        if p1 is None:
            p1 = p2
        if p2 is None:
            p2 = p1

        centroid1 = (xyz2 + xyz1) / 2.
        if p1 > p2:
            dp = p1 - p2
            pnominal = p2
            centroid2 = (2*xyz1 + xyz2) / 3.
        else:
            dp = p2 - p1
            centroid2 = (2*xyz2 + xyz1) / 3.
            pnominal = p1

        r1 = centroid1 - p
        r2 = centroid2 - p
        f1 = pnominal * area_edge * load_dir * scale
        f2 = dp * area_edge * load_dir * scale
        m1 = cross(r1, f1)
        m2 = cross(r2, f2)
        fi += f1 + f2
        mi += m1 + m2

    cid = load.Cid()
    if cid != 0:
        raise NotImplementedError('cid=%r nvector=%s on a PLOAD4 is not supported\n%s' % (
            cid, load.nvector, str(load)))
    return fi, mi
