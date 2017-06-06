# coding: utf-8
"""
This file defines:
  - WriteMeshes
"""
from __future__ import (nested_scopes, generators, division, absolute_import,
                        print_function, unicode_literals)
import sys
from codecs import open
from six import PY2, iteritems

from pyNastran.bdf.field_writer_8 import print_card_8
from pyNastran.bdf.field_writer_16 import print_card_16
from pyNastran.bdf.field_writer_double import print_card_double
from pyNastran.bdf.bdf_interface.write_mesh import WriteMesh


class WriteMeshes(WriteMesh):
    """
    Defines methods for writing cards

    Major methods:
      - model.write_bdf_symmetric(...)
    """
    def __init__(self):
        """creates methods for writing cards"""
        WriteMesh.__init__(self)

    def write_bdf_symmetric(self, out_filename=None, encoding=None,
                            size=8, is_double=False,
                            enddata=None, close=True, plane='xz'):
        """
        Writes the BDF as a symmetric model.

        Parameters
        ----------
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
        #self.write_caero_model()
        out_filename = self._output_helper(out_filename,
                                           interspersed, size, is_double)
        if encoding is not None:
            pass
        else:
            encoding = self._encoding
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
        self._write_header(bdf_file, encoding)
        self._write_params(bdf_file, size, is_double)
        self._write_nodes_symmetric(bdf_file, size, is_double, plane=plane)

        if interspersed:
            raise RuntimeError(interspersed)
            #self._write_elements_properties(bdf_file, size, is_double)
        else:
            self._write_elements_symmetric(bdf_file, size, is_double)
            self._write_properties(bdf_file, size, is_double)
        self._write_materials(bdf_file, size, is_double)

        self._write_masses(bdf_file, size, is_double)
        self._write_common(bdf_file, size, is_double)
        if (enddata is None and 'ENDDATA' in self.card_count) or enddata:
            bdf_file.write('ENDDATA\n')
        if close:
            bdf_file.close()

    def _write_nodes_symmetric(self, bdf_file, size=8, is_double=False, plane='xz'):
        """
        Writes the NODE-type cards

        .. warning:: doesn't consider coordinate systems;
                      it could, but you'd need 20 new coordinate systems
        .. warning:: doesn't mirror SPOINTs, EPOINTs
        """
        if self.spoints:
            msg = []
            msg.append('$SPOINTS\n')
            msg.append(self.spoints.write_card(size, is_double))
            bdf_file.write(''.join(msg))
        if self.epoints:
            msg = []
            msg.append('$EPOINTS\n')
            msg.append(self.epoints.write_card(size, is_double))
            bdf_file.write(''.join(msg))

        plane = plane.strip().lower()
        if plane == 'xz':
            iy = 4
        elif plane == 'xy':
            iy = 5
        elif plane == 'yz':
            iy = 3
        else:
            raise NotImplementedError(plane)
        if self.nodes:
            msg = []
            msg.append('$NODES\n')
            if self.grdset:
                msg.append(self.grdset.print_card(size))

            nid_offset = max(self.nodes.keys())
            if self.is_long_ids:
                print_card_long = print_card_double if is_double else print_card_16
                for (unused_nid, node) in sorted(iteritems(self.nodes)):
                    repr_fields = node.repr_fields()
                    msg.append(print_card_long(repr_fields))
                    repr_fields = node.repr_fields()
                    # grid, nid, cp, x, y, z
                    repr_fields[1] += nid_offset
                    repr_fields[iy] *= -1.0
                    msg.append(print_card_long(repr_fields))
            else:
                for (unused_nid, node) in sorted(iteritems(self.nodes)):
                    repr_fields = node.repr_fields()
                    msg.append(print_card_8(repr_fields))
                    repr_fields = node.repr_fields()
                    # grid, nid, cp, x, y, z
                    repr_fields[1] += nid_offset
                    repr_fields[iy] *= -1.0
                    msg.append(print_card_8(repr_fields))
            bdf_file.write(''.join(msg))
        #if 0:  # not finished
            #self._write_nodes_associated(bdf_file, size, is_double)

    def _write_elements_symmetric(self, bdf_file, size=8, is_double=False):
        """
        Writes the elements in a sorted order
        """
        nid_offset = max(self.nodes.keys())
        eid_offset = max(self.elements.keys())
        if self.elements:
            bdf_file.write('$ELEMENTS\n')
            if self.is_long_ids:
                for (eid, element) in sorted(iteritems(self.elements)):
                    nodes = element.node_ids
                    bdf_file.write(element.write_card_16(is_double))
                    element.eid += eid_offset
                    nodes = [node_id + nid_offset for node_id in nodes]
                    element.nodes = nodes
                    bdf_file.write(element.write_card_16(is_double))
            else:
                for (eid, element) in sorted(iteritems(self.elements)):
                    nodes = element.node_ids
                    bdf_file.write(element.write_card(size, is_double))
                    try:
                        nodes = [node_id + nid_offset for node_id in nodes]
                    except TypeError:
                        msg = 'cannot mirror %r because None exists in nodes=%s' % (
                            element.type, nodes)
                        self.log.warning(msg)
                        continue

                    if element.type in ['CTRIA3', 'CQUAD4']:
                        nodes = nodes[::-1]
                    try:
                        element.nodes = nodes
                    except AttributeError:
                        msg = 'cannot mirror %r because it doesnt have nodes...' % element.type
                        self.log.warning(msg)
                        continue
                    element.eid += eid_offset
                    bdf_file.write(element.write_card(size, is_double))

