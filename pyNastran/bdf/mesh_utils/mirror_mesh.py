# coding: utf-8
"""
This file defines:
  - model, nid_offset, eid_offset = bdf_mirror(bdf_filename, plane='xz')
  - model, nid_offset, eid_offset = write_bdf_symmetric(
        bdf_filename, out_filename=None, encoding=None,
        size=8, is_double=False,
        enddata=None, close=True, plane='xz')
  - model = make_symmetric_model(
        bdf_filename, plane='xz', zero_tol=1e-12,
        log=None, debug=True)

"""
from __future__ import (nested_scopes, generators, division, absolute_import,
                        print_function, unicode_literals)
from warnings import warn

import numpy as np

from pyNastran.bdf.cards.loads.static_loads import (
    FORCE, FORCE1, FORCE2, MOMENT, MOMENT1, MOMENT2,
    PLOAD, PLOAD2, PLOAD4)
from pyNastran.bdf.cards.thermal.loads import TEMP, QVOL, QBDY1, QBDY2, QBDY3, QHBDY
from pyNastran.bdf.cards.aero.aero import CAERO1, SPLINE1
from pyNastran.bdf.cards.bdf_sets import SET1 #, SET3
from pyNastran.bdf.bdf import BDF, read_bdf

def get_model(bdf_filename, log=None, debug=True):
    """helper method"""
    if isinstance(bdf_filename, BDF):
        model = bdf_filename
    else:
        # str, StringIO
        model = read_bdf(bdf_filename, validate=True, xref=True,
                         punch=False, skip_cards=None,
                         read_cards=None,
                         encoding=None, log=log,
                         debug=debug, mode='msc')
    return model

def bdf_mirror_plane(bdf_filename, plane, mirror_model=None,
                     log=None, debug=True, use_nid_offset=True):
    """mirrors a model about an arbitrary plane"""
    model = get_model(bdf_filename, log=log, debug=debug)
    if mirror_model is None:
        mirror_model = BDF(debug=debug, log=log, mode='msc')

    nid_offset, plane = _mirror_nodes_plane(model, mirror_model, plane,
                                            use_nid_offset=use_nid_offset)
    eid_offset = _mirror_elements(model, mirror_model, nid_offset)
    #_mirror_loads(model, nid_offset, eid_offset)
    return model, mirror_model, nid_offset, eid_offset


def bdf_mirror(bdf_filename, plane='xz', log=None, debug=True):
    """
    Mirrors the model about the symmetry plane

    Parameters
    ----------
    bdf_filename : str / BDF()
        str : the bdf filename
        BDF : the BDF model object
    plane : str; {'xy', 'yz', 'xz'}; default='xz'
        the plane to mirror about
        xz : +y/-y
        yz : +x/-x
        xy : +z/-z

    Returns
    -------
    model : BDF()
        BDF : the BDF model object
    nid_offset : int
        the offset node id
    eid_offset : int
        the offset element id

    """
    model = get_model(bdf_filename, log=log, debug=debug)
    mirror_model = model
    nid_offset, plane = _mirror_nodes(model, plane=plane)
    eid_offset = _mirror_elements(model, mirror_model, nid_offset)
    _mirror_loads(model, nid_offset, eid_offset)
    _mirror_aero(model, nid_offset, plane=plane)
    return model, nid_offset, eid_offset

def write_bdf_symmetric(bdf_filename, out_filename=None, encoding=None,
                        size=8, is_double=False,
                        enddata=None, close=True, plane='xz', log=None):
    """
    Mirrors the model about the symmetry plane

    Parameters
    ----------
    bdf_filename : str / BDF()
        str : the bdf filename
        BDF : the BDF model object
    out_filename : varies; default=None
        str        - the name to call the output bdf
        file       - a file object
        StringIO() - a StringIO object
        None       - pops a dialog
    encoding : str; default=None -> system specified encoding
        the unicode encoding
        latin1, and utf8 are generally good options
    size : int; {8, 16}
        the field size
    is_double : bool; default=False
        False : small field
        True : large field
    enddata : bool; default=None
        bool - enable/disable writing ENDDATA
        None - depends on input BDF
    close : bool; default=True
        should the output file be closed
    plane : str; {'xy', 'yz', 'xz'}; default='xz'
        the plane to mirror about
        xz : +y/-y
        yz : +x/-x
        xy : +z/-z

    Returns
    -------
    model : BDF()
        BDF : the BDF model object
    nid_offset : int
        the offset node id
    eid_offset : int
        the offset element id

    Notes
    -----
    Updates the BDF object to be symmetric
      - see bdf_mirror if you don't want to write the model
    Doesn't equivalence nodes on the centerline.

    Considers
     - nodes : GRID
     - elements, rigid_elements, mass_elements : see ``_mirror_elements``
     - loads : see ``_mirror_loads``
     - aero cards : see ``_mirror_aero``

    """
    model, nid_offset, eid_offset = bdf_mirror(bdf_filename, plane=plane, log=None)
    model.write_bdf(out_filename=out_filename, encoding=encoding,
                    size=size, is_double=is_double,
                    interspersed=False, enddata=enddata, close=close)
    return model, nid_offset, eid_offset

def _mirror_nodes(model, plane='xz'):
    """
    Mirrors the GRIDs

    .. warning:: doesn't consider coordinate systems;
                  it could, but you'd need 20 new coordinate systems
    .. warning:: doesn't mirror SPOINTs, EPOINTs

    """
    nid_offset = 0
    iy, plane = _plane_to_iy(plane)
    if model.nodes:
        nid_offset = max(model.node_ids)
        for (nid, node) in sorted(model.nodes.items()):
            xyz = node.get_position()
            nid2 = nid + nid_offset
            xyz2 = xyz.copy()
            xyz2[iy] *= -1.
            model.add_grid(nid2, xyz2, cp=0, cd=node.cd, ps=node.ps, seid=node.seid)
    return nid_offset, plane

def _mirror_nodes_plane(model, mirror_model, plane, use_nid_offset=True):
    """
    Mirrors the GRIDs about an arbitrary plane

    .. warning:: doesn't consider coordinate systems;
                  it could, but you'd need 20 new coordinate systems
    .. warning:: doesn't mirror SPOINTs, EPOINTs

    https://mathinsight.org/distance_point_plane

    """
    nid_offset = 0

    if model.nodes:
        all_nodes, xyz_cid0 = model.get_xyz_in_coord_no_xref(cid=0, fdtype='float64', sort_ids=True)
        unused_cid = max(model.coords) + 1
        origin = plane[0, :]

        # just a triangle's normal vector
        n1 = plane[0, :]
        n2 = plane[1, :]
        n3 = plane[2, :]
        normal = np.cross(n2 - n1, n3 - n1)
        normal /= np.linalg.norm(normal)

        #cord2r = model.add_cord2r(cid, plane[0, :], plane[1, :], plane[2, :])
        #del model.coords[cid]
        #print(cord2r)

        #origin = cord2r.origin
        #normal = cord2r.beta()[2, :]
        #print('normal =', normal)
        vector = xyz_cid0 - origin
        assert xyz_cid0.shape == vector.shape, 'xyz_cid0.shape=%s; vector.shape=%s' % (xyz_cid0.shape, vector.shape)
        v_dot_n = vector * normal[np.newaxis, :]
        assert v_dot_n.shape == vector.shape, 'v_dot_n.shape=%s; vector.shape=%s' % (v_dot_n.shape, vector.shape)
        distance = np.linalg.norm(v_dot_n, axis=1)
        assert v_dot_n.shape[0] == len(distance), 'v_dot_n.shape=%s; distance.shape=%s' % (v_dot_n.shape, distance.shape)

        # we're some distance from the plane, but we don't know the
        # direction, so we take the max distance from the plane and
        # project it in the +normal direction and then check that
        # distance in comparison to the known distance
        #
        max_distance = distance.max()
        imax = np.where(distance == max_distance)[0][0]
        distance0 = distance[imax]
        xyz0 = xyz_cid0[imax, :] + distance0 * normal
        v0_dot_n = xyz0 * normal
        distance_plus = np.linalg.norm(v0_dot_n)

        if distance_plus > 1.1*distance0:
            xyz_cid0_2 = xyz_cid0 - 2 * distance[:, np.newaxis] * normal[np.newaxis, :]
        else:
            xyz_cid0_2 = xyz_cid0 + 2 * distance[:, np.newaxis] * normal[np.newaxis, :]

        if use_nid_offset:
            nid_offset = max(all_nodes)
        for nid, xyz2 in zip(all_nodes, xyz_cid0_2):
            node = model.nodes[nid]
            nid2 = nid + nid_offset
            mirror_model.add_grid(nid2, xyz2, cp=0, cd=node.cd, ps=node.ps, seid=node.seid)
    return nid_offset, plane

def _plane_to_iy(plane):
    """gets the index fo the mirror plane"""
    plane = plane.strip().lower()
    if plane == 'yz':
        iy = 0
    if plane == 'xz':
        iy = 1
    elif plane == 'xy':
        iy = 2
    else:  # pragma: no cover
        raise NotImplementedError("plane=%r and must be 'yz', 'xz', or 'xy'." % plane)
    return iy, plane

def _mirror_elements(model, mirror_model, nid_offset, use_eid_offset=True):
    """
    Mirrors the elements

    elements:
       0d : CELAS1, CELAS2, CELAS3, CELAS4, CDAMP1, CDAMP2, CDAMP3, CDAMP4, CDAMP5
            CFAST, CBUSH, CBUSH1D
       2d : CTRIA3, CQUAD4, CTRIA6, CQUAD8, CQUAD, CTRIAR, CQUADR
       3d : ???
       missing : CVISC, CTRIAX, CTRIAX6, CQUADX, CQUADX8, CCONEAX
    rigid_elements:
       loaded: RBE2, RBE3
       missing: RBAR, RBAR1
    mass_elements:
       loaded: CONM2
       missing CONM1, CMASS1, CMASS2, CMASS3, CMASS4

    Notes
    -----
    Doesn't handle CBAR/CBEAM offsets
    Doesn't handle CBEAM SPOINTs

    """
    eid_max_elements = 0
    eid_max_masses = 0
    eid_max_rigid = 0
    if use_eid_offset:
        if model.elements:
            eid_max_elements = max(model.elements.keys())
        if model.masses:
            eid_max_masses = max(model.masses.keys())
        if model.rigid_elements:
            eid_max_rigid = max(model.rigid_elements.keys())
    eid_offset = max(eid_max_elements, eid_max_masses, eid_max_rigid)

    if model.elements:
        shells = {'CTRIA3', 'CQUAD4', 'CTRIA6', 'CQUAD8', 'CQUAD', 'CTRIAR', 'CQUADR'}
        rods = {'CROD', 'CONROD', 'CTUBE'}
        spring_dampers = {
            'CELAS1', 'CELAS2', 'CELAS3', 'CELAS4',
            'CDAMP1', 'CDAMP2', 'CDAMP3', 'CDAMP4', 'CDAMP5',
        }

        eid_offset = max(model.elements)
        def _set_nodes(element, nodes):
            try:
                element.nodes = nodes
            except AttributeError:
                print(element.get_stats())
                print(nodes)
                raise

        for eid, element in sorted(model.elements.items()):
            etype = element.type
            if etype in ['CHBDYG', 'CHBDYE']:
                continue

            nodes = element.node_ids
            #try:
                #nodes = [node_id + nid_offset for node_id in nodes]
            #except TypeError:
                #msg = 'cannot mirror %r eid=%s because None exists in nodes=%s' % (
                    #element.type, eid, nodes)
                #model.log.warning(msg)
                #continue

            eid_mirror = eid + eid_offset
            fields = element.repr_fields()
            fields[1] = eid_mirror
            mirror_model.add_card(fields, etype)
            element2 = mirror_model.elements[eid_mirror]

            if etype in shells:
                nodes = [node_id + nid_offset for node_id in nodes]
                element2.nodes = nodes
                element.flip_normal() # nodes = nodes[::-1]
            elif etype in rods:
                nodes = [node_id + nid_offset for node_id in nodes]
                element2.nodes = nodes
            elif etype in ['CBAR', 'CBEAM']:
                nodes = [node_id + nid_offset for node_id in nodes]
                element2.nodes = nodes
                g0 = element2.g0 + nid_offset if element2.g0 is not None else None
                element2.g0 = g0

            elif etype == 'CGAP':
                #nodes = [node_id + nid_offset if node_id is not None else None
                         #for node_id in nodes]
                #_set_nodes(element2, nodes)
                ga = element2.ga + nid_offset if element2.ga is not None else None
                gb = element2.gb + nid_offset if element2.gb is not None else None
                g0 = element2.g0 + nid_offset if element2.g0 is not None else None
                element2.ga = ga
                element2.gb = gb
                element2.g0 = g0
            elif etype == 'CBUSH1D':
                ga = element2.ga + nid_offset if element2.ga is not None else None
                gb = element2.gb + nid_offset if element2.gb is not None else None
                element2.ga = ga
                element2.gb = gb
            elif etype == 'CFAST':
                ga = element2.ga + nid_offset if element2.ga is not None else None
                gb = element2.gb + nid_offset if element2.gb is not None else None
                gs = element2.gs + nid_offset if element2.gs is not None else None
                element2.ga = ga
                element2.gb = gb
                element2.gs = gs

            elif etype == 'CCONEAX':
                pass
            elif etype in [spring_dampers]:
                nodes = [node_id + nid_offset if node_id is not None else None for node_id in nodes]
                _set_nodes(element2, nodes)
                #print(nodes)
                #element2.nodes = nodes
            elif etype == 'GENEL':
                element2.ul = element2.ul + nid_offset
                element2.ud = element2.ud + nid_offset
            else:
                try:
                    element2.nodes = nodes
                except AttributeError:  # pragma: no cover
                    print(element.get_stats())
                    raise

    if model.masses:
        for eid, element in sorted(model.masses.items()):
            eid_mirror = eid + eid_offset

            #print(element.get_stats())
            if element.type == 'CONM2':
                old_nid = element.nid
                new_nid = old_nid + nid_offset
                mass = element.mass
                cid = element.cid
                X = element.X
                I = element.I
                mirror_model.add_conm2(eid_mirror, new_nid, mass, cid=cid, X=X, I=I, comment='')
            else:  # pragma: no cover
                #print(element.get_stats())
                mirror_model.log.warning('skipping mass element:\n%s' % str(element))

    if model.rigid_elements:
        for eid, rigid_element in sorted(model.rigid_elements.items()):
            if rigid_element.type == 'RBE2':
                Gmi_node_ids = rigid_element.Gmi_node_ids
                Gn = rigid_element.Gn()
                Gijs = None
                ref_grid_id = None
            elif rigid_element.type == 'RBE3':
                Gmi_node_ids = rigid_element.Gmi_node_ids
                Gijs = rigid_element.Gijs
                ref_grid_id = rigid_element.ref_grid_id
                Gn = None
            else:
                msg = '_write_elements_symmetric: %s not implemented' % rigid_element.type
                warn(msg)
                continue
                #raise NotImplementedError(msg)

            Gmi_node_ids_mirror = [node_id + nid_offset for node_id in Gmi_node_ids]
            if Gn:
                Gn_mirror = Gn + nid_offset
            if Gijs:
                Gijs_mirror = [[node_id + nid_offset for node_id in nodes] for nodes in Gijs]
            if ref_grid_id:
                ref_grid_id_mirror = ref_grid_id + nid_offset

            eid_mirror = eid + eid_offset
            if rigid_element.type == 'RBE2':
                rigid_element2 = mirror_model.add_rbe2(eid_mirror, Gn_mirror, rigid_element.cm,
                                                       Gmi_node_ids_mirror)
            elif rigid_element.type == 'RBE3':
                rigid_element2 = mirror_model.add_rbe3(
                    eid_mirror, ref_grid_id_mirror, rigid_element.refc, rigid_element.weights,
                    rigid_element.comps, Gijs_mirror
                )
            else:  # pragma: no cover
                mirror_model.log.warning('skipping:\n%s' % str(rigid_element))

    return eid_offset

def _mirror_loads(model, nid_offset=0, eid_offset=0):
    """
    Mirrors the loads.  A mirrored force acts in the same direction.

    Considers:
     - PLOAD4
        - no coordinate systems (assumes cid=0)
     - TEMP, QHBDY, QBDY1, QBDY2, QBDY3

    """
    for unused_load_id, loads in model.loads.items():
        for load in loads:
            loads_new = []
            load_type = load.type
            if load_type == 'PLOAD4':
                g1 = None
                g34 = None
                if load.g1 is not None:
                    g1 = load.g1 + nid_offset
                if load.g34 is not None:
                    g34 = load.g34 + nid_offset

                eids = [eid + eid_offset for eid in load.eids]
                load = PLOAD4(
                    load.sid, eids, load.pressures, g1, g34,
                    cid=load.cid, nvector=load.nvector,
                    surf_or_line=load.surf_or_line,
                    line_load_dir=load.line_load_dir, comment='')
                loads_new.append(load)

            elif load_type == 'FORCE':
                load = FORCE(load.sid, load.node + nid_offset, load.mag, load.xyz,
                              cid=load.cid, comment='')
                loads_new.append(load)
            elif load_type == 'FORCE1':
                load = FORCE1(load.sid, load.node + nid_offset, load.mag,
                              load.g1 + nid_offset, load.g2 + nid_offset, comment='')
                loads_new.append(load)
            elif load_type == 'FORCE2':
                load = FORCE2(load.sid, load.node + nid_offset, load.mag,
                              load.g1 + nid_offset, load.g2 + nid_offset,
                              load.g3 + nid_offset, load.g4 + nid_offset, comment='')
                loads_new.append(load)

            elif load_type == 'MOMENT':
                load = MOMENT(load.sid, load.node + nid_offset, load.mag, load.xyz,
                              cid=load.cid, comment='')
                loads_new.append(load)
            elif load_type == 'MOMENT1':
                load = MOMENT1(load.sid, load.node + nid_offset, load.mag,
                               load.g1 + nid_offset, load.g2 + nid_offset, comment='')
                loads_new.append(load)
            elif load_type == 'MOMENT2':
                load = MOMENT2(load.sid, load.node + nid_offset, load.mag,
                               load.g1 + nid_offset, load.g2 + nid_offset,
                               load.g3 + nid_offset, load.g4 + nid_offset, comment='')
                loads_new.append(load)

            elif load_type == 'PLOAD':
                nodes = [nid + nid_offset for nid in load.nodes]
                load = PLOAD(load.sid, load.pressure, nodes, comment='')
                loads_new.append(load)
            elif load_type == 'PLOAD2':
                eids = [eid + eid_offset for eid in load.eids]
                load = PLOAD2(load.sid, load.pressure, eids, comment='')
                loads_new.append(load)

            elif load_type == 'QVOL':
                elements = [eid + eid_offset for eid in load.elements]
                load = QVOL(load.sid, load.qvol, nid_offset + load.control_point, elements)
                loads_new.append(load)
            elif load_type == 'QHBDY':
                grids = [nid + nid_offset for nid in load.grids]
                load = QHBDY(load.sid, load.flag, load.q0, grids, af=load.af)
                loads_new.append(load)
            elif load_type == 'QBDY1':
                eids = [eid + eid_offset for eid in load.eids]
                load = QBDY1(load.sid, load.qflux, eids)
                loads_new.append(load)
            elif load_type == 'QBDY2':
                load = QBDY2(load.sid, load.eid + eid_offset, load.qfluxs, comment='')
                loads_new.append(load)
            elif load_type == 'QBDY3':
                eids = [eid + eid_offset for eid in load.eids]
                load = QBDY3(load.sid, load.q0, load.cntrlnd + nid_offset, eids)
                loads_new.append(load)

            elif load_type == 'TEMP':
                temperatures = {}
                for nid, temp in load.temperatures.items():
                    temperatures[nid + nid_offset] = temp
                load = TEMP(load.sid, temperatures)
                loads_new.append(load)
            elif load_type == 'GRAV':
                pass
            else:  # pragma: no cover
                model.log.warning('skipping:\n%s' % load.rstrip())
        if loads_new:
            loads += loads_new

def _mirror_aero(model, nid_offset, plane):
    """
    Mirrors the aero elements

    Considers:
     - AEROS
      - doesn't consider sideslip
     - CAERO1
      - doesn't consider sideslip
      - doesn't consider lspan/lchord
     - SPLINE1
     - SET1

    """
    if model.aeros is not None:
        aeros = model.aeros
        aeros.sref *= 2.
        if plane == 'xz':
            aeros.sym_xz = 0
        elif plane == 'yz':
            aeros.sym_yz = 0
        else:
            model.log.error('not mirroring plane %r; only xz, yz' % plane)

    caero_id_offset = 0
    if len(model.caeros):
        caero_id_max = max(model.caero_ids)
        caero_id_offset = np.max(model.caeros[caero_id_max].box_ids.flat)

        caeros = []
        for unused_caero_id, caero in model.caeros.items():
            if caero.type == 'CAERO1':
                assert caero.lspan == 0, caero
                assert caero.lchord == 0, caero
                lchord = caero.lchord
                nchord = caero.nchord
                lspan = caero.lspan
                nspan = caero.nspan
                p1 = caero.p1.copy()
                p1[1] *= -1.
                x12 = caero.x12
                p4 = caero.p4.copy()
                p4[1] *= -1.
                x43 = caero.x43
                eid2 = caero.eid + caero_id_offset
                caero_new = CAERO1(eid2, caero.pid, caero.igroup,
                                   p1, x12, p4, x43,
                                   cp=caero.cp, nspan=nspan,
                                   lspan=lspan, nchord=nchord, lchord=lchord,
                                   comment='')

                # we flip the normal so if we ever use W2GJ it's going to be consistent
                caero_new.flip_normal()
                caeros.append(caero_new)
                #print(caero)
            else:  # pragma: no cover
                model.log.error('skipping (only support CAERO1):\n%s' % caero.rstrip())

        for caero in caeros:
            model._add_caero_object(caero)

    nsplines = len(model.splines)
    sets_max = max(model.sets) if len(model.sets) else 0
    if caero_id_offset == 0 and nsplines:
        model.log.error("cant mirror splines because CAEROs don't exist...")
    elif nsplines and sets_max == 0:
        model.log.error("cant mirror splines because SET1/3 don't exist...")
    elif nsplines:
        splines = []
        spline_sets_to_duplicate = []
        spline_max = max(model.splines)
        for unused_spline_id, spline in model.splines.items():
            if spline.type == 'SPLINE1':
                #spline = SPLINE1(eid, caero, box1, box2, setg)

                eid = spline.eid + spline_max
                caero = spline.caero + caero_id_offset
                method = spline.method
                usage = spline.usage
                box1 = spline.box1 + caero_id_offset
                box2 = spline.box2 + caero_id_offset
                setg = spline.setg + sets_max
                dz = spline.dz
                melements = spline.melements
                nelements = spline.nelements
                spline_new = SPLINE1(eid, caero, box1, box2, setg, dz=dz,
                                     method=method, usage=usage,
                                     nelements=nelements, melements=melements, comment='')
                splines.append(spline_new)
                spline_sets_to_duplicate.append(spline.setg)
            else:  # pragma: no cover
                model.log.error('skipping (only support SPLINE1):\n%s' % spline.rstrip())

        #print("spline_sets_to_duplicate =", spline_sets_to_duplicate)
        msg = ', which is required to mirror:\n%s' % spline.rstrip()

        sets_to_add = []
        for set_id in spline_sets_to_duplicate:
            set_card = model.Set(set_id, msg=msg)
            if set_card.type == 'SET1':
                sid = set_card.sid + sets_max
                ids = [nid + nid_offset for nid in set_card.ids]
                is_skin = set_card.is_skin
                set_card = SET1(sid, ids, is_skin=is_skin, comment='')
                sets_to_add.append(set_card)
            else:  # pragma: no cover
                model.log.error('skipping (only support SET1):\n%s' % set_card.rstrip())

        for spline in splines:
            model._add_spline_object(spline)
        for set_card in sets_to_add:
            model._add_set_object(set_card)

    model.pop_parse_errors()

def make_symmetric_model(bdf_filename, plane='xz', zero_tol=1e-12, log=None, debug=True):
    """
    Makes a symmetric model from a full model

    Parameters
    ----------
    bdf_filename : str / BDF()
        str : the bdf filename
        BDF : the BDF model object
    plane : str; {'xy', 'yz', 'xz'}; default='xz'
        the plane to mirror about
        xz : +y/-y
        yz : +x/-x
        xy : +z/-z
    zaero_tol : float; default=1e-12
        the symmetry plane tolerance

    Returns
    -------
    model : BDF()
        BDF : the BDF model object

    ## TODO: doesn't handle elements straddling the centerline

    """
    model = get_model(bdf_filename, log=log, debug=debug)
    iy, plane = _plane_to_iy(plane)
    nids_to_remove = []
    eids_to_remove = []
    caero_ids_to_remove = []
    zero = -zero_tol
    for eid, elem in model.elements.items():
        xyz = elem.Centroid()

        if xyz[iy] < zero:
            eids_to_remove.append(eid)

    for nid, node in model.nodes.items():
        xyz = node.get_position()
        if xyz[iy] < zero:
            nids_to_remove.append(nid)

    for nid in nids_to_remove:
        del model.nodes[nid]

    for eid in eids_to_remove:
        del model.elements[eid]

    for caero_id, caero in model.caeros.items():
        if caero.type == 'CAERO1':
            p1, p2, p3, p4 = caero.get_points()
            #print(caero)
            if p1[iy] <= zero and p4[iy] <= zero:
                #print('p1=%s p4=%s' % (p1, p4))
                caero_ids_to_remove.append(caero_id)
            elif p1[iy] < zero:
                p1[iy] = 0.
                caero.set_points([p1, p2, p3, p4])
            elif p4[iy] < zero:
                p4[iy] = 0.
                caero.set_points([p1, p2, p3, p4])
        elif caero.type == 'CAERO2':
            # TODO: a CAERO2 can't be half symmetric...can it?
            # TODO: it can be skewed though...
            p1, p2 = caero.get_points()
            if p1[iy] <= zero and p2[iy] <= zero:
                #print('p1=%s p4=%s' % (p1, p4))
                caero_ids_to_remove.append(caero_id)
        else:  # pragma: no cover
            raise NotImplementedError(caero)

    for caero_id in caero_ids_to_remove:
        del model.caeros[caero_id]

    #print('nids_to_remove =', nids_to_remove)
    for unused_spline_id, spline in model.splines.items():
        caero = spline.caero
        #setg = spline.setg
        #print('caero = ', caero)
        nids = spline.setg_ref.ids  # list
        #spline.uncross_reference()

        #i = 0
        nids = list(set(nids) - set(nids_to_remove))
        nids.sort()
        spline.setg_ref.ids_ref = None
        spline.setg_ref.ids = nids

    plane_to_labels_keep_map = {
        'yz' : ['URDD4', 'URDD2', 'URDD3', 'SIDES', 'YAW'], # yz
        'xz' : ['URDD1', 'URDD5', 'URDD3', 'PITCH', 'ANGLEA'], # xz plane
        'xy' : ['URDD1', 'URDD2', 'URDD6', 'ROLL'], # xy plane
    }

    all_labels = {
        'URDD4', 'URDD2', 'URDD3', 'SIDES', 'YAW',
        'URDD1', 'URDD5', 'URDD3', 'PITCH', 'ANGLEA',
        'URDD1', 'URDD2', 'URDD6', 'ROLL',
    }
    labels_to_keep = plane_to_labels_keep_map[plane]
    labels_to_remove = [label for label in all_labels if label not in labels_to_keep]

    #print('labels_to_remove =', labels_to_remove)
    for aestat_id in list(model.aestats.keys()):
        aestat = model.aestats[aestat_id]
        if aestat.label in labels_to_remove:
            del model.aestats[aestat_id]

    for unused_trim_id, trim in model.trims.items():
        labels = trim.labels
        ilabels_to_remove = [labels.index(label) for label in labels_to_remove
                             if label in labels]
        #print("ilabels_to_remove =", ilabels_to_remove)
        trim.uxz = [trim.uxs[ilabel] for ilabel in ilabels_to_remove]
        trim.labels = [trim.labels[ilabel] for ilabel in ilabels_to_remove]
    return model
