# coding: utf-8
"""
This file defines:
  - write_bdf_symmetric(model, out_filename=None, encoding=None,
                        size=8, is_double=False,
                        enddata=None, close=True, plane='xz')
"""
from __future__ import (nested_scopes, generators, division, absolute_import,
                        print_function, unicode_literals)
import sys
from codecs import open
from copy import copy
from six import PY2, iteritems

import numpy as np

from pyNastran.bdf.field_writer_8 import print_card_8
from pyNastran.bdf.field_writer_16 import print_card_16
from pyNastran.bdf.field_writer_double import print_card_double
from pyNastran.bdf.cards.nodes import write_xpoints

def write_bdf_symmetric(model, out_filename=None, encoding=None,
                        size=8, is_double=False,
                        enddata=None, close=True, plane='xz'):
    """
    Writes the BDF as a symmetric model.

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
    """
    interspersed = False
    #model.write_caero_model()
    out_filename = model._output_helper(out_filename,
                                        interspersed, size, is_double)
    if encoding is not None:
        pass
    else:
        encoding = model._encoding
        if encoding is None:
            encoding = sys.getdefaultencoding()
    #assert encoding.lower() in ['ascii', 'latin1', 'utf8'], encoding

    if hasattr(out_filename, 'read') and hasattr(out_filename, 'write'):
        bdf_file = out_filename
    else:
        if PY2:
            wb = 'wb'
        else:
            wb = 'w'
        bdf_file = open(out_filename, wb, encoding=encoding)
    model._write_header(bdf_file, encoding)
    model._write_params(bdf_file, size, is_double)
    nid_offset = _write_nodes_symmetric(model, bdf_file, size, is_double, plane=plane)

    if interspersed:
        raise RuntimeError(interspersed)
        #model._write_elements_properties(bdf_file, size, is_double)
    else:
        _write_elements_symmetric(model, bdf_file, size, is_double)
        model._write_properties(bdf_file, size, is_double)
    model._write_materials(bdf_file, size, is_double)

    model._write_masses(bdf_file, size, is_double)
    caero_id_map = _write_aero_symmetric(model, bdf_file, size=size, is_double=is_double,
                                         plane=plane)
    model._write_common(bdf_file, size, is_double)
    if (enddata is None and 'ENDDATA' in model.card_count) or enddata:
        bdf_file.write('ENDDATA\n')
    if close:
        bdf_file.close()

    return nid_offset, caero_id_map

def _get_iy(plane):
    plane = plane.strip().lower()
    if plane == 'yz':
        iy = 0
    elif plane == 'xz':
        iy = 1
    elif plane == 'xy':
        iy = 2
    else:
        raise NotImplementedError(plane)
    return iy

def _write_nodes_symmetric(model, bdf_file, size=8, is_double=False, plane='xz'):
    """
    Writes the NODE-type cards

    .. warning:: doesn't consider coordinate systems;
                  it could, but you'd need 20 new coordinate systems
    .. warning:: doesn't mirror SPOINTs, EPOINTs
    """
    if model.spoints:
        msg = []
        msg.append('$SPOINTS\n')
        msg.append(write_xpoints('SPOINT', model.spoints.keys()))
        bdf_file.write(''.join(msg))
    if model.epoints:
        msg = []
        msg.append('$EPOINTS\n')
        msg.append(write_xpoints('EPOINT', model.epoints.keys()))
        bdf_file.write(''.join(msg))

    iy = _get_iy(plane) + 3

    if model.nodes:
        msg = []
        msg.append('$NODES\n')
        if model.grdset:
            msg.append(model.grdset.print_card(size))

        nid_offset = max(model.nodes.keys())
        if model.is_long_ids:
            print_card_long = print_card_double if is_double else print_card_16
            for (unused_nid, node) in sorted(iteritems(model.nodes)):
                repr_fields = node.repr_fields()
                msg.append(print_card_long(repr_fields))
                repr_fields = node.repr_fields()
                # grid, nid, cp, x, y, z
                repr_fields[1] += nid_offset
                repr_fields[iy] *= -1.0
                msg.append(print_card_long(repr_fields))
        else:
            for (unused_nid, node) in sorted(iteritems(model.nodes)):
                repr_fields = node.repr_fields()
                msg.append(print_card_8(repr_fields))
                repr_fields = node.repr_fields()
                # grid, nid, cp, x, y, z
                repr_fields[1] += nid_offset
                repr_fields[iy] *= -1.0
                msg.append(print_card_8(repr_fields))
        bdf_file.write(''.join(msg))
    #if 0:  # not finished
        #model._write_nodes_associated(bdf_file, size, is_double)
        return nid_offset

def _write_elements_symmetric(model, bdf_file, size=8, is_double=False):
    """
    Writes the elements in a sorted order
    """
    nid_offset = max(model.nodes.keys())
    eid_offset = max(max(model.elements.keys()), max(model.rigid_elements.keys()))
    if model.elements:
        bdf_file.write('$ELEMENTS\n')
        if model.is_long_ids:
            for (eid, element) in sorted(iteritems(model.elements)):
                nodes = element.node_ids
                bdf_file.write(element.write_card_16(is_double))
                element.eid += eid_offset
                nodes = [node_id + nid_offset for node_id in nodes]
                element.nodes = nodes
                bdf_file.write(element.write_card_16(is_double))
        else:
            for (eid, element) in sorted(iteritems(model.elements)):
                nodes = element.node_ids
                bdf_file.write(element.write_card(size, is_double))
                try:
                    nodes = [node_id + nid_offset for node_id in nodes]
                except TypeError:
                    msg = 'cannot mirror %r because None exists in nodes=%s' % (
                        element.type, nodes)
                    model.log.warning(msg)
                    continue

                if element.type in ['CTRIA3', 'CQUAD4']:
                    nodes = nodes[::-1]

                element.uncross_reference()
                try:
                    element.nodes = nodes
                except AttributeError:
                    msg = 'cannot mirror %r because it doesnt have nodes...' % element.type
                    model.log.warning(msg)
                    continue
                element.eid += eid_offset
                bdf_file.write(element.write_card(size, is_double))
    if model.rigid_elements:
        for (eid, rigid_element) in sorted(iteritems(model.rigid_elements)):
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

            bdf_file.write(rigid_element.write_card(size, is_double))

            Gmi_node_ids = [node_id + nid_offset for node_id in Gmi_node_ids]
            if Gn:
                Gn += nid_offset
            if Gijs:
                Gijs = [[node_id + nid_offset for node_id in nodes] for nodes in Gijs]
            if ref_grid_id:
                ref_grid_id += nid_offset

            rigid_element.uncross_reference()
            if rigid_element.type == 'RBE2':
                rigid_element.Gmi = Gmi_node_ids
                rigid_element.gn = Gn
            elif rigid_element.type == 'RBE3':
                rigid_element.Gmi = Gmi_node_ids
                rigid_element.Gijs = Gijs
                rigid_element.refgrid = ref_grid_id

            rigid_element.eid += eid_offset
            bdf_file.write(rigid_element.write_card(size, is_double))

def _write_aero_symmetric(model, bdf_file, size=8, is_double=False, plane='xz'):
    """
    Writes the aero in a sorted order
    """
    caero_id_max = max(model.caero_ids)
    caero_id_offset = np.max(model.caeros[caero_id_max].box_ids.flat)

    nid_offset = max(model.nodes.keys())
    set_id_offset = max(model.sets.keys())
    spline_id_offset = max(model.splines.keys())
    eid_offset = max(max(model.elements.keys()), max(model.rigid_elements.keys()))

    iy = _get_iy(plane)

    caero_id_map = {}
    if model.caeros:
        bdf_file.write('$AERO\n')
        for (caero_id, caero) in sorted(iteritems(model.caeros)):
            bdf_file.write(caero.write_card(size, is_double))

            # reverse points to maintain normal
            p1_mirror = caero.p4
            x12_mirror = caero.x43
            p4_mirror = caero.p1
            x43_mirror = caero.x12

            p1_mirror[iy] *= -1.
            p4_mirror[iy] *= -1.

            caero_mirror = copy(caero)
            caero_mirror.p1 = p1_mirror
            caero_mirror.p4 = p4_mirror
            caero_mirror.x12 = x12_mirror
            caero_mirror.x43 = x43_mirror
            caero_mirror.eid += caero_id_offset
            caero_id_map[caero_id] = caero_mirror.eid

            bdf_file.write(caero_mirror.write_card(size, is_double))

    # paero not mirrored
    for (unused_id, paero) in sorted(iteritems(model.paeros)):
        bdf_file.write(paero.write_card(size, is_double))

    if model.splines:
        bdf_file.write('$SPLINES\n')
        for (spline_id, spline) in sorted(iteritems(model.splines)):
            bdf_file.write(spline.write_card(size, is_double))
            #bdf_file.write(spline.setg_ref.write_card(size, is_double))

            set_copy = copy(spline.setg_ref)
            if set_copy is None:
                set_copy = copy(model.sets[spline.setg])
            set_copy.ids = [nid + nid_offset for nid in set_copy.ids]
            set_copy.ids_ref = None
            set_copy.sid += set_id_offset

            spline.box1 += caero_id_offset
            spline.box2 += caero_id_offset
            spline.caero = spline.CAero() + caero_id_offset
            spline.caero_ref = None
            spline.eid += spline_id_offset
            spline.setg = set_copy.sid
            spline.setg_ref = None

            bdf_file.write(spline.write_card(size, is_double))
            bdf_file.write(set_copy.write_card(size, is_double)) # non mirrored sets writen elsewhere

    # paero not mirrored
    for monitor_point in model.monitor_points:
        bdf_file.write(monitor_point.write_card(size, is_double))

    return caero_id_map
