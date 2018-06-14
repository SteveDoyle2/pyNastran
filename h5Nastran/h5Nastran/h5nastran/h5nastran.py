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
        defaults = self.h5f.get_node(self.table_paths.defaults).read()
        self.defaults.load(self)
