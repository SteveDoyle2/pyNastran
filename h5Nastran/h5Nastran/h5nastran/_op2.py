from __future__ import print_function, absolute_import

from collections import OrderedDict

from ._result_base import H5NastranResultBase

# from h5Nastran.post_process.result_readers.punch import PunchReader

import numpy as np
import tables
from six import iteritems


class H5NastranResultOP2(H5NastranResultBase):
    def __init__(self, *args, **kwargs):
        super(H5NastranResultOP2, self).__init__(*args, **kwargs)

    def load_op2(self, filename):
        if self._bdf is None:
            raise Exception('BDF must be loaded first!')

        if self._punch is not None:
            raise Exception('Punch file has already been loaded.  Cannot load op2 file after punch.')

        if self._f06 is not None:
            raise Exception('F06 file has already been loaded.  Cannot load op2 file after f06.')

        self._op2 = filename

        # TODO: need OP2Reader
        reader = OP2Reader(filename)
        reader.register_callback(self._load_op2_table)
        reader.read()

        self.h5f.flush()

        for table in self._tables:
            table.finalize()

        self._tables.clear()
        self._write_unsupported_tables()

    def _load_op2_table(self, table_data):
        result_type = table_data.result_type

        table = self._result_tables.get(result_type, None)

        if table is None:
            return self._unsupported_table(table_data)

        table.write_op2_data(table_data)

        self._tables.add(table)
