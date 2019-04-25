from __future__ import print_function, absolute_import

from h5Nastran.h5.nastran import NastranNode
from h5Nastran.defaults import Defaults
from h5Nastran.table_paths import TablePaths

import tables


class H5NastranBase(object):
    h5n_version_str = '0.1.0a0'
    h5n_version = (0, 1, 0)
    nastran_type = None
    nastran_version = (0, 0, 0)
    default_driver = None

    def __init__(self, h5filename, mode='r', nastran_type=None, nastran_version=None, in_memory=False):
        self.file_mode = ''
        self.h5filename = ''
        self.h5f = None

        self._open_file(h5filename, mode, nastran_type, nastran_version, in_memory)

        self._card_tables = {}
        self._result_tables = {}
        
        self.nastran = NastranNode(self)
        
        # this node will take care of all custom data in h5 file
        # self.h5nastran = H5NastranNode(self)

        self.table_paths = TablePaths()
        self.defaults = Defaults()

    def close(self):
        self.h5f.close()
        
    def path(self):
        return ['']
    
    def _open_file(self, h5filename, mode, nastran_type, nastran_version, in_memory):
        model = self.file_mode 
        
        if mode == 'r':
            assert nastran_type is None and nastran_version is None
        else:
            if nastran_type is not None:
                self.nastran_type = nastran_type

            if nastran_version is not None:
                assert isinstance(nastran_version, tuple) and len(nastran_version) == 3
                self.nastran_version = nastran_version

        if in_memory:
            driver = 'H5FD_CORE'
            driver_core_backing_store = 0
        else:
            driver = self.default_driver
            driver_core_backing_store = 1

        filters = tables.Filters(complib='zlib', complevel=5)
        self.h5f = tables.open_file(h5filename, mode=mode, filters=filters, driver=driver,
                                    driver_core_backing_store=driver_core_backing_store)
