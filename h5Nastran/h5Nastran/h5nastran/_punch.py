from __future__ import print_function, absolute_import

from collections import OrderedDict

from ._result_base import H5NastranResultBase

from h5Nastran.post_process.result_readers.punch import PunchReader

import numpy as np
import tables
from six import iteritems


class H5NastranResultPunch(H5NastranResultBase):
    def __init__(self, *args, **kwargs):
        super(H5NastranResultPunch, self).__init__(*args, **kwargs)
        
    def load_punch(self, filename):
        if self._bdf is None:
            raise Exception('BDF must be loaded first!')

        if self._f06 is not None:
            raise Exception('F06 has already been loaded.  Cannot load punch file after f06.')

        self._punch = filename
        self._punch_subcase_ids.clear()

        reader = PunchReader(filename)
        reader.register_callback(self._load_punch_table)
        reader.read()

        self.h5f.flush()

        for table in self._tables:
            table.finalize()

        self._tables.clear()
        self._write_unsupported_tables()
        self._punch_finalize()

    def _punch_finalize(self):
        dtype = np.dtype([('SUBCASE_ID', '<i8'), ('LOAD_FACTOR', '<f8'), ('DOMAIN_ID', '<i8')])
        format = tables.descr_from_dtype(dtype)[0]

        self.h5f.create_table(self.table_paths.subcase_path, self.table_paths.subcase_table, format,
                              'SUBCASES', expectedrows=len(self._punch_subcase_ids), createparents=True)

        table = self.h5f.get_node(self.table_paths.subcase)

        data = np.zeros(len(self._punch_subcase_ids), dtype=dtype)
        subcase_id = data['SUBCASE_ID']
        load_factor = data['LOAD_FACTOR']
        domain_id = data['DOMAIN_ID']

        for key, domain_id_ in iteritems(self._punch_subcase_ids):
            index = domain_id_ - 1
            subcase_id_, load_factor_ = key
            subcase_id[index] = subcase_id_
            load_factor[index] = load_factor_
            domain_id[index] = domain_id_

        table.append(data)

        self.h5f.flush()

    def _load_punch_table(self, table_data):
        key = table_data.header.subcase_id_num, table_data.header.load_factor

        if key not in self._punch_subcase_ids:
            self._punch_subcase_ids[key] = len(self._punch_subcase_ids) + 1

        results_type = table_data.header.results_type_basic

        table = self._result_tables.get(results_type, None)

        if table is None:
            return self._unsupported_table(table_data)

        table.write_punch_data(table_data)

        self._tables.add(table)
