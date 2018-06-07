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
from pyNastran.bdf.cards.loads.static_loads import PLOAD4

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
    _mirror_loads(model, nid_offset, eid_offset)
    return nid_offset, eid_offset

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
    nid_offset, eid_offset = bdf_mirror(model, plane=plane)

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
        raise NotImplementedError("plane=%r and must be 'yz', 'xz', or 'xy'." % plane)
    return iy

def _mirror_elements(model, nid_offset):
    """Mirrors the elements"""
    eid_offset = max(model.elements.keys())
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

            eid2 = eid + eid_offset
            fields = element.repr_fields()
            fields[1] = eid2
            model.add_card(fields, etype)
            element2 = model.elements[eid2]

            if etype in ['CTRIA3', 'CQUAD4']:
                element2.nodes = nodes
                element.flip_normal() # nodes = nodes[::-1]
            else:
                try:
                    element2.nodes = nodes
                except AttributeError:
                    print(element.get_stats())
                    raise
    return eid_offset

def _mirror_loads(model, nid_offset=0, eid_offset=0):
    """
    Mirrors the loads

    Considers:
     - PLOAD4
        - no coordinate systems (assumes cid=0)
    """
    for load_id, loads in iteritems(model.loads):
        nloads = len(loads)
        for iload, load in enumerate(loads):
            # TODO: a super hack due to us changing the length of loads
            #       I'm confused why this is required though as I thought
            #       iterators can't change...
            if iload == nloads:
                break
            load_type = load.type
            if load.type == 'PLOAD4':
                g1 = None
                g34 = None
                if load.g1 is not None:
                    g1 = load.g1 + nid_offset
                if load.g34 is not None:
                    g34 = load.g34 + nid_offset
                print(load.rstrip())

                eids = [eid + eid_offset for eid in load.eids]
                load = model.add_pload4(
                    load.sid, eids, load.pressures, g1, g34,
                    cid=load.cid, nvector=load.nvector,
                    surf_or_line=load.surf_or_line,
                    line_load_dir=load.line_load_dir, comment='')
            else:
                model.log.warning('skipping:\n%s' % load.rstrip())

def make_symmetric_model(model, iy=1, zero_tol=1e-12):
    """
    Makes a symmetric model

    ## TODO: doesn't handle elements straddling the centerline
    """
    nids_to_remove = []
    eids_to_remove = []
    caero_ids_to_remove = []
    zero = -zero_tol
    for eid, elem in iteritems(model.elements):
        xyz = elem.Centroid()

        if xyz[iy] < zero:
            eids_to_remove.append(eid)

    for nid, node in iteritems(model.nodes):
        xyz = node.get_position()
        if xyz[iy] < zero:
            nids_to_remove.append(nid)

    for nid in nids_to_remove:
        print('**', nid)
        del model.nodes[nid]

    for eid in eids_to_remove:
        del model.elements[eid]

    for caero_id, caero in iteritems(model.caeros):
        if caero.type == 'CAERO1':
            p1, p2, p3, p4 = caero.get_points()
            #print(caero)
            if p1[iy] <= zero and p4[iy] <= zero:
                #print('p1=%s p4=%s' % (p1, p4))
                caero_ids_to_remove(caero_id)
            elif p1[iy] < zero:
                p1[iy] = 0.
                caero.set_points([p1, p2, p3, p4])
            elif p4[iy] < zero:
                p4[iy] = 0.
                caero.set_points([p1, p2, p3, p4])
        elif caero.type == 'CAERO2':
            # TODO: a CAERO2 can't be half symmetric...
            # TODO: it can be skewed though...
            p1, p2 = caero.get_points()
            if p1[iy] <= zero and p2[iy] <= zero:
                #print('p1=%s p4=%s' % (p1, p4))
                caero_ids_to_remove(caero_id)
        else:
            raise NotImplementedError(caero)

    for caero_id in caero_ids_to_remove:
        del model.caeros[caero_id]

    print('nids_to_remove =', nids_to_remove)
    for spline_id, spline in iteritems(model.splines):
        caero = spline.caero
        setg = spline.setg
        #print('caero = ', caero)
        nids = spline.setg_ref.ids  # list
        #spline.uncross_reference()

        i = 0
        nids = list(set(nids) - set(nids_to_remove))
        nids.sort()
        spline.setg_ref.ids_ref = None
        spline.setg_ref.ids = nids

    plane_to_labels_keep_map = {
        0 : ['URDD4', 'URDD2', 'URDD3', 'SIDES', 'YAW'], # yz
        1 : ['URDD1', 'URDD5', 'URDD3', 'PITCH', 'ANGLEA'], # xz plane
        2 : ['URDD1', 'URDD2', 'URDD6', 'ROLL'], # xy plane
    }

    all_labels = [
        'URDD4', 'URDD2', 'URDD3', 'SIDES', 'YAW',
        'URDD1', 'URDD5', 'URDD3', 'PITCH', 'ANGLEA',
        'URDD1', 'URDD2', 'URDD6', 'ROLL',
    ]
    labels_to_keep = plane_to_labels_keep_map[iy]
    labels_to_remove = [label for label in all_labels if label not in labels_to_keep]

    print('labels_to_remove =', labels_to_remove)
    for aestat_id in list(model.aestats.keys()):
        aestat = model.aestats[aestat_id]
        if aestat.label in labels_to_remove:
            del model.aestats[aestat_id]

    for trid_id, trim in iteritems(model.trims):
        labels = trim.labels
        ilabels_to_remove = [labels.index(label) for label in labels_to_remove
                             if label in labels]
        print("ilabels_to_remove =", ilabels_to_remove)
        trim.uxz = [trim.uxs[ilabel] for ilabel in ilabels_to_remove]
        trim.labels = [trim.labels[ilabel] for ilabel in ilabels_to_remove]
