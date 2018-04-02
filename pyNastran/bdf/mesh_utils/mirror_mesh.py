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
from six import iteritems

from pyNastran.bdf.field_writer_8 import print_card_8
from pyNastran.bdf.field_writer_16 import print_card_16
from pyNastran.bdf.field_writer_double import print_card_double
from pyNastran.bdf.cards.nodes import write_xpoints

def write_bdf_symmetric(model, out_filename=None, encoding=None,
                        size=8, is_double=False,
                        enddata=None, close=True, plane='xz'):
    """
    Writes the BDF as a symmetric model.
    Considers shell and line elements.

    Does not equivalence nodes on the centerline.
    Does not handle coordinate systems.
    Does not handle solid elements.
    Does not handle loads.

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
    """
    interspersed = False
    #model.write_caero_model()
    out_filename = model._output_helper(
        out_filename, interspersed, size, is_double)
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
        bdf_file = open(out_filename, 'w', encoding=encoding)
    model._write_header(bdf_file, encoding)
    model._write_params(bdf_file, size, is_double)
    nid_offset = _write_nodes_symmetric(model, bdf_file, size, is_double, plane=plane)

    if interspersed:
        raise RuntimeError(interspersed)
        #model._write_elements_properties(bdf_file, size, is_double)
    else:
        _write_elements_symmetric(model, bdf_file, nid_offset,
                                  size=size, is_double=is_double)
        model._write_properties(bdf_file, size, is_double)
    model._write_materials(bdf_file, size, is_double)

    model._write_masses(bdf_file, size, is_double)
    model._write_common(bdf_file, size, is_double)
    if (enddata is None and 'ENDDATA' in model.card_count) or enddata:
        bdf_file.write('ENDDATA\n')
    if close:
        bdf_file.close()

def _write_nodes_symmetric(model, bdf_file, size=8, is_double=False, plane='xz'):
    """
    Writes the NODE-type cards

    .. warning:: doesn't consider coordinate systems;
                  it could, but you'd need 20 new coordinate systems
    .. warning:: doesn't mirror SPOINTs, EPOINTs
    """
    nid_offset = 0
    if model.spoints:
        msg = []
        msg.append('$SPOINTS\n')
        spoints = model.spoints.keys()
        msg.append(write_xpoints('SPOINT', spoints))
        bdf_file.write(''.join(msg))
        nid_offset = max(spoints)

    if model.epoints:
        msg = []
        msg.append('$EPOINTS\n')
        epoints = model.epoints.keys()
        msg.append(write_xpoints('EPOINT', epoints))
        bdf_file.write(''.join(msg))
        nid_offset = max(nid_offset, max(epoints))

    iy = _plane_to_iy(plane) + 3
    if model.nodes:
        msg = []
        msg.append('$NODES\n')
        if model.grdset:
            msg.append(model.grdset.print_card(size))

        nids = model.nodes.keys()
        nid_offset = max(nid_offset, max(nids))
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

def _write_elements_symmetric(model, bdf_file, nid_offset, size=8, is_double=False):
    """
    Writes the elements in a sorted order
    """
    eid_offset = max(model.elements.keys())
    if model.elements:
        bdf_file.write('$ELEMENTS\n')
        if model.is_long_ids:
            for unused_eid, element in sorted(iteritems(model.elements)):
                bdf_file.write(element.write_card_16(is_double))
                is_passed = _update_element(element, nid_offset, eid_offset, model.log)
                if is_passed:
                    bdf_file.write(element.write_card_16(is_double))
        else:
            for unused_eid, element in sorted(iteritems(model.elements)):
                bdf_file.write(element.write_card(size, is_double))
                is_passed = _update_element(element, nid_offset, eid_offset, model.log)
                if is_passed:
                    bdf_file.write(element.write_card(size, is_double))

def _update_element(element, nid_offset, eid_offset, log):
    """mirrors the element"""
    is_passed = False
    nodes = element.node_ids
    try:
        nodes = [node_id + nid_offset for node_id in nodes]
    except TypeError:
        msg = 'cannot mirror %r because None exists in nodes=%s' % (
            element.type, nodes)
        log.warning(msg)
        return is_passed

    element.uncross_reference()
    if element.type in ['CTRIA3', 'CQUAD4']:
        element.flip_normal() # nodes = nodes[::-1]
    #elif element.type in ['CTRIA3', 'CQUAD4']:
        #nodes = nodes[::-1]
    try:
        element.nodes = nodes
    except AttributeError:
        msg = 'cannot mirror %r because it doesnt have nodes...' % element.type
        log.warning(msg)
        return is_passed
    element.eid += eid_offset
    is_passed = True
    return is_passed
