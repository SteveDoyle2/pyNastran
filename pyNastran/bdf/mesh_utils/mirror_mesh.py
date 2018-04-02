# coding: utf-8
"""
This file defines:
  - write_bdf_symmetric(model, out_filename=None, encoding=None,
                        size=8, is_double=False,
                        enddata=None, close=True, plane='xz')
"""
from __future__ import (nested_scopes, generators, division, absolute_import,
                        print_function, unicode_literals)
from six import iteritems


def bdf_mirror(model, plane='xz'):
    """
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
        raise NotImplementedError(plane)
    return iy

def _mirror_elements(model, nid_offset):
    """
    Mirrors the elements
    """
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
