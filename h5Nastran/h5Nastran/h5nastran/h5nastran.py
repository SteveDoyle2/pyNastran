from __future__ import print_function, absolute_import

from collections import OrderedDict

import numpy as np
import tables
from six import iteritems

from ._result import H5NastranResult
from ._bdf import H5NastranBDF


class H5Nastran(H5NastranBDF, H5NastranResult):
    def __init__(self, h5filename, mode='r', nastran_type=None, nastran_version=None, in_memory=False):
        super(H5Nastran, self).__init__(h5filename, mode=mode, nastran_type=nastran_type,
                                        nastran_version=nastran_version, in_memory=in_memory)

        if mode == 'w':
            self._write_info()
        else:
            self._update()

    def visualize(self):
        from ..gui.visualization import to_vtk

        if self.bdf is None:
            self.load_bdf()

        vtk_data = to_vtk(self.bdf)
        vtk_data.visualize()

    def _write_info(self):
        import pyNastran

        info = 'h5Nastran version %s\nPowered by pyNastran version %s' % (self.h5n_version_str, pyNastran.__version__)

        self.h5f.create_array(self.table_paths.about_path, self.table_paths.about_table, obj=info.encode(),
                              title='h5Nastran Info', createparents=True)

        versioning_dtype = np.dtype([
            ('H5NASTRAN_VERSION_STR', 'S8'),
            ('H5NASTRAN_VERSION', '<i8', (3,)),
            ('NASTRAN_TYPE', 'S8'),
            ('NASTRAN_VERSION', '<i8', (3,))
        ])

        format = tables.descr_from_dtype(versioning_dtype)[0]

        self.h5f.create_table(self.table_paths.versioning_path, self.table_paths.versioning_table, format,
                              'VERSIONING', expectedrows=1, createparents=True)

        table = self.h5f.get_node(self.table_paths.versioning)
        data = np.zeros(1, dtype=versioning_dtype)

        data['H5NASTRAN_VERSION_STR'][0] = self.h5n_version_str
        data['H5NASTRAN_VERSION'][0] = self.h5n_version

        nastran_type = self.nastran_type

        if nastran_type is None:
            nastran_type = ''

        data['NASTRAN_TYPE'][0] = nastran_type
        data['NASTRAN_VERSION'][0] = self.nastran_version

        table.append(data)
        self.defaults.save(self)
        self.h5f.flush()

    def _update(self):
        self.nastran.update()
        self.defaults.load(self)

    def element_search(self, elem_types=None, elem_pids=None, elem_ids=None, box=None, partial_fit=True):
        class Dummy(object):
            def __contains__(self, val):
                return True

            def __le__(self, val):
                return True

            def __ge__(self, val):
                return True

        dummy = Dummy()

        if elem_types is None:
            elem_types = dummy

        if elem_pids is None:
            elem_pids = dummy

        if elem_ids is None:
            elem_ids = dummy

        if box is None:
            min_x = dummy
            max_x = dummy
            min_y = dummy
            max_y = dummy
            min_z = dummy
            max_z = dummy
        else:
            x, y, z = box

            if x is None:
                min_x = dummy
                max_x = dummy
            else:
                min_x, max_x = x
            if y is None:
                min_y = dummy
                max_y = dummy
            else:
                min_y, max_y = y
            if z is None:
                min_z = dummy
                max_z = dummy
            else:
                min_z, max_z = z

            if min_x is None:
                min_x = dummy
            if min_y is None:
                min_y = dummy
            if min_z is None:
                min_z = dummy
            if max_x is None:
                max_x = dummy
            if max_y is None:
                max_y = dummy
            if max_z is None:
                max_z = dummy

        elms = []
        bdf = self.bdf

        for eid, elm in iteritems(bdf.elements):
            pid = elm.pid
            etype = elm.type

            if etype not in elem_types:
                continue
            if pid not in elem_pids:
                continue
            if eid not in elem_ids:
                continue

            nodes = elm.node_ids
            in_box = set()
            bdf_nodes = bdf.nodes

            for nid in nodes:
                node = bdf_nodes[nid]
                x, y, z = node.get_position()
                if min_x <= x <= max_x and min_y <= y <= max_y and min_z <= z <= max_z:
                    in_box.add(True)
                else:
                    in_box.add(False)

            if partial_fit:
                if True in in_box:
                    elms.append(elm)
            else:
                if False not in in_box:
                    elms.append(elm)

        return elms

