from __future__ import print_function, absolute_import

from collections import OrderedDict

import numpy as np
import tables
from six import iteritems

from ._base import H5NastranBase


class H5NastranResultBase(H5NastranBase):
    def __init__(self, *args, **kwargs):
        super(H5NastranResultBase, self).__init__(*args, **kwargs)
        self._tables = set()
        self._unsupported_tables = set()

        self._punch_subcase_ids = OrderedDict()
        self._element_results_tables = {}

        self._punch = None
        self._f06 = None
        self._op2 = None
        
    def register_result_table(self, result_table):
        result_type = result_table.result_type

        if isinstance(result_type, str):
            result_type = [result_type]

        for _result_type in result_type:
            assert _result_type not in self._result_tables

            if 'ELEMENT' in _result_type:
                tmp = _result_type.split()
                try:
                    int(tmp[-3])
                    del tmp[-2]
                except ValueError:
                    pass
                _result_type = ' '.join(tmp)

            self._result_tables[_result_type] = result_table

    # TODO: remove element type from results tables, since it's hard to know what they should be
    # item code is good enough
    # this means remove the QUAD4 in ELEMENT STRESSES QUAD4 33 REAL
    def _find_element_result_table(self, results_type):
        if not results_type.startswith('ELEMENT'):
            return None

        tmp = results_type.split(' ')
        del tmp[2]
        results_type = ' '.join(tmp)

        if len(self._element_results_tables) == 0:
            keys = self._result_tables.keys()
            for key in keys:
                if not key.startswith('ELEMENT'):
                    continue
                tmp = key.split(' ')
                del tmp[2]
                tmp = ' '.join(tmp)
                self._element_results_tables[tmp] = key

        _results_type = self._element_results_tables.get(results_type, None)

        return self._result_tables.get(_results_type, None)

    def _unsupported_table(self, table_data):
        print('Unsupported table %s' % table_data.header.results_type)
        self._unsupported_tables.add(table_data.header.results_type)

    def _write_unsupported_tables(self):
        headers = list(sorted(self._unsupported_tables))
        data = np.array(headers, dtype='S256')

        self.h5f.create_array(self.table_paths.unsupported_result_tables_path,
                              self.table_paths.unsupported_result_tables_table, obj=data,
                              title='UNSUPPORTED RESULT TABLES', createparents=True)
