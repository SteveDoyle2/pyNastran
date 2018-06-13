# coding: utf-8
"""
This file defines:
  - nid_offset, eid_offset = bdf_mirror(model, plane='xz')
  - nid_offset, eid_offset = write_bdf_symmetric(
        model, out_filename=None, encoding=None,
        size=8, is_double=False,
        enddata=None, close=True, plane='xz')
"""
from __future__ import (nested_scopes, generators, division, absolute_import,
                        print_function, unicode_literals)
from six import iteritems


def bdf_mirror(model, plane='xz'):
    """
    Mirrors the model

    Parameters
    ----------
    model : BDF()
        the BDF model object
    plane : str; {'xy', 'yz', 'xz'}; default='xz'
        the plane to mirror about
        xz : +y/-y
        yz : +x/-x
        xy : +z/-z

    Returns
    -------
    nid_offset : int
        the offset node id
    eid_offset : int
        the offset element id
    """
    nid_offset = _mirror_nodes(model, plane=plane)
    eid_offset = _mirror_elements(model, nid_offset)
    caero_id_map = _mirror_aero(model, nid_offset, eid_offset, plane=plane)
    return nid_offset, eid_offset, caero_id_map

def write_bdf_symmetric(model, out_filename=None, encoding=None,
                        size=8, is_double=False,
                        enddata=None, close=True, plane='xz'):
    """
    Writes the BDF as a symmetric model.
    Considers shell and line elements.
    Updates the BDF object to be symmetric
      - see bdf_mirror if you don't want to write the model

    Doesn't equivalence nodes on the centerline.
    Doesn't handle solid elements.
    Doesn't consider mass elements.
    Doesn't consider rigid elements.
    Doesn't handle loads.
    Doesn't handle bar/beam offsets.

    Parameters
    ----------
    model : BDF()
        the BDF model object
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
    nid_offset : int
        the offset node id
    eid_offset : int
        the offset element id
    """
    #model.write_caero_model()
    nid_offset, eid_offset, caero_id_map = bdf_mirror(model, plane=plane)

    model.write_bdf(out_filename=out_filename, encoding=encoding,
                    size=size, is_double=is_double,
                    interspersed=False, enddata=enddata, close=close)
    return nid_offset, eid_offset

def _mirror_nodes(model, plane='xz'):
    """
    Mirrors the GRIDs

    .. warning:: doesn't consider coordinate systems;
                  it could, but you'd need 20 new coordinate systems
    .. warning:: doesn't mirror SPOINTs, EPOINTs
    """
    nid_offset = max(model.node_ids)

    iy = _plane_to_iy(plane)
    if model.nodes:
        for (nid, node) in sorted(iteritems(model.nodes)):
            xyz = node.get_position()
            nid2 = nid + nid_offset
            xyz2 = xyz.copy()
            xyz2[iy] *= -1.
            model.add_grid(nid2, xyz2, cp=0, cd=node.cd, ps=node.ps, seid=node.seid)
    return nid_offset

def _plane_to_iy(plane):
    """gets the index fo the mirror plane"""
    plane = plane.strip().lower()
    if plane == 'yz':
        iy = 0
    if plane == 'xz':
        iy = 1
    elif plane == 'xy':
        iy = 2
    else:
        raise NotImplementedError(plane)
    return iy

def _mirror_elements(model, nid_offset):
    """
    Mirrors the elements
    """
    eid_max_elements = 0
    eid_max_rigid = 0
    if model.elements:
        eid_max_elements = max(model.elements.keys())
    if model.rigid_elements:
        eid_max_rigid = max(model.rigid_elements.keys())
    eid_offset = max(eid_max_elements, eid_max_rigid)
    if model.elements:
        for eid, element in sorted(iteritems(model.elements)):
            etype = element.type
            if etype in ['CHBDYG', 'CHBDYE']:
                continue

            nodes = element.node_ids
            try:
                nodes = [node_id + nid_offset for node_id in nodes]
            except TypeError:
                msg = 'cannot mirror %r eid=%s because None exists in nodes=%s' % (
                    element.type, eid, nodes)
                model.log.warning(msg)
                continue

            eid_mirror = eid + eid_offset
            fields = element.repr_fields()
            fields[1] = eid_mirror
            model.add_card(fields, etype)
            element2 = model.elements[eid_mirror]

            if etype in ['CTRIA3', 'CQUAD4']:
                element2.nodes = nodes
                element.flip_normal() # nodes = nodes[::-1]
            else:
                try:
                    element2.nodes = nodes
                except AttributeError:
                    print(element.get_stats())
                    raise

    if model.rigid_elements:
        for eid, rigid_element in sorted(iteritems(model.rigid_elements)):
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
                msg = '_write_elements_symmetric: %s not implimented' % rigid_element.type
                raise NotImplementedError(msg)

            Gmi_node_ids_mirror = [node_id + nid_offset for node_id in Gmi_node_ids]
            if Gn:
                Gn_mirror = Gn + nid_offset
            if Gijs:
                Gijs_mirror = [[node_id + nid_offset for node_id in nodes] for nodes in Gijs]
            if ref_grid_id:
                ref_grid_id_mirror = ref_grid_id + nid_offset

            eid_mirror = eid + eid_offset
            if rigid_element.type == 'RBE2':
                rigid_element2 = model.add_rbe2(eid_mirror, Gn_mirror, rigid_element.cm,
                                                Gmi_node_ids_mirror)
            elif rigid_element.type == 'RBE3':
                rigid_element2 = model.add_rbe3(
                    eid_mirror, ref_grid_id_mirror, rigid_element.refc, rigid_element.weights,
                    rigid_element.comps, Gijs_mirror
                )

    return eid_offset

def _mirror_aero(model, nid_offset, eid_offset, plane='xz'):
    """
    Writes the aero in a sorted order
    """
    caero_id_map = {}
    if model.caeros:
        caero_id_max = max(model.caero_ids)
        caero_id_offset = np.max(model.caeros[caero_id_max].box_ids.flat)

        set_id_offset = max(model.sets.keys())
        spline_id_offset = max(model.splines.keys())

        iy = _plane_to_iy(plane)

        if model.caeros:
            for (caero_id, caero) in sorted(iteritems(model.caeros)):
                if caero.type == 'CAERO1':
                    # reverse points to maintain normal
                    p1_mirror = caero.p4
                    x12_mirror = caero.x43
                    p4_mirror = caero.p1
                    x43_mirror = caero.x12

                    p1_mirror[iy] *= -1.
                    p4_mirror[iy] *= -1.

                    from pyNastran.bdf.bdf import BDF
                    isinstance(model, BDF)

                    caero_id_mirror = caero_id + caero_id_offset
                    caero_id_map[caero_id] = caero_id_mirror

                    caero_mirror = model.add_caero1(
                        caero_id_mirror, caero.pid, caero.igid, p1_mirror, x12_mirror, p4_mirror,
                        x43_mirror, cp=caero.cp, nspan=caero.nspan, lspan=caero.lspan,
                        nchord=caero.nchord, lchord=caero.lchord
                    )
                else:
                    model.log.warning('%s not supported in bdf_mirror:\n%s' % (caero.type, caero))

        if model.splines:
            for (spline_id, spline) in sorted(iteritems(model.splines)):
                if spline.type == 'SPLINE1':
                    set_ref = spline.setg_ref
                    if set_ref is None:
                        set_ref = model.sets[spline.setg]
                    set_ids_mirror = [nid + nid_offset for nid in set_ref.ids]
                    sid_mirror = set_ref.sid + set_id_offset

                    model.add_set1(sid_mirror, set_ids_mirror)

                    box1_mirror = spline.box1 + caero_id_offset
                    box2_mirror = spline.box2 + caero_id_offset
                    caero_mirror = spline.CAero() + caero_id_offset
                    eid_mirror = spline.eid + spline_id_offset
                    setg_mirror = sid_mirror

                    model.add_spline1(eid_mirror, caero_mirror, box1_mirror, box2_mirror, setg_mirror,
                                      dz=spline.dz)
                else:
                    model.log.warning('%s not supported in bdf_mirror:\n%s' % (spline.type, spline))

    return caero_id_map
