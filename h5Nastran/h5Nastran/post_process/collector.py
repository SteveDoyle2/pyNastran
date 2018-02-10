from __future__ import print_function, absolute_import

from collections import OrderedDict

from copy import deepcopy

from six.moves import range
from six import iteritems

from typing import List, Dict

from h5Nastran import H5Nastran

from pandas import DataFrame
import numpy as np


class CollectorSubcase(object):
    def __init__(self, file_index=-1, loadcase=-1, factor=-1.):
        self.file_index = file_index
        self.loadcase = loadcase
        self.factor = factor


class CollectorLoadcase(object):
    def __init__(self):
        self.id = -1
        self.description = ''

        self.subcases = []
        """:type: list[CollectorSubcase]"""

    def add_subcase(self, file_index, loadcase, factor):
        self.subcases.append(CollectorSubcase(file_index, loadcase, factor))

    def check_file_index(self):
        try:
            first_index = self.subcases[0].file_index
        except IndexError:
            return True

        for i in range(1, len(self.subcases)):
            if self.subcases[i].file_index != first_index:
                raise ValueError('File indices for subcase must be uniform!')
            
    def model_subcases(self):
        subcases = []
        for subcase in self.subcases:
            subcases.append((subcase.file_index, subcase.loadcase))
        return subcases


class H5NastranCollector(object):
    def __init__(self, filename=None):
        self.filename = None

        self.cases = []  # type: List[CollectorLoadcase]
        self.model_ids = []  # type: List[str]
        self.filenames = []  # type: List[str]
        self.h5n = []  # type: List[H5Nastran]

        self.by_case_id = {}  # type: Dict[int, int]
        self.by_case_description = {}  # type: Dict[str, int]
        
        self.loadcase_selection = []  # type: List[int]
        
        self.model_case_selection = []  # type: List[List[int]]
        self._model_case_selection = []  # type: List[Dict[int, int]]

        if filename is not None:
            self.load_file(filename)

    def __del__(self):
        self.close()
            
    def combine(self, data):
        # type: (List[DataFrame]) -> DataFrame
        
        # if isinstance(data, DataFrame):
        #     data = [data]
            
        assert len(data) == len(self.h5n)

        # FIXME: why is dataframe not using correct dtype?
        _result = np.zeros(len(self.cases), dtype=data[0].the_dtype)
        data_cols = [str(_) for _ in data[0].data_cols]
        
        case_ids = []
        
        for i in range(len(self.cases)):
            case = self.cases[i]
            case_ids.append(case.id)
            for subcase in case.subcases:
                file_index, loadcase, factor = subcase.file_index, subcase.loadcase, subcase.factor
                # TODO: need to handle loadcase selection better
                j = loadcase
                
                _data = data[file_index]
                
                for data_col in data_cols:
                    # print(data_col, i, j)
                    # print(_result[data_col][i])
                    # print(data[file_index][data_col][j])
                    _result[data_col][i] += factor * _data[data_col][j]
                
        _result['LOADCASE'] = case_ids
                
        return data[0].from_records(_result)

    def _set_filenames(self, filenames):
        del self.filenames[:]
        self.filenames[:] = filenames[:]
        
    def open(self):
        if len(self.h5n) == len(self.filenames):
            # TODO: need something better than this
            return
        
        dbs = []
        
        for fname in self.filenames:
            dbs.append(H5Nastran(fname, 'r'))
            
        del self.h5n[:]
        self.h5n.extend(dbs)
        
    def close(self):
        for db in self.h5n:
            db.close()
        del self.h5n[:]

    def load_file(self, filename):
        with open(filename, 'r') as f:
            lines = f.read().split('\n')

        self.close()

        self.filename = filename

        h5files = OrderedDict()
        subcases = OrderedDict()
        loadcases = OrderedDict()

        self.by_case_description.clear()
        self.by_case_id.clear()

        del self.model_ids[:]
        del self.filenames[:]
        del self.cases[:]

        status = None

        for line in lines:
            if line.startswith('#'):
                continue

            if line.strip() == '':
                continue

            _line = line.lower()

            if _line.startswith('database'):
                status = 0
                continue
            elif _line.startswith('subcase'):
                status = 1
                continue
            elif _line.startswith('loadcase'):
                status = 2
                continue

            if status == 0:
                data = line.split(' ')
                if len(data) > 2:
                    raise ValueError('Cannot have spaces in file path! %s' % line[1:])
                h5files[data[0]] = data[1]

            elif status == 1:
                data = line.split(' ')
                subcases[data[0]] = [data[1], int(data[2])]

            elif status == 2:
                data = line.split(' ')
                case_id = data[0]
                case_description = data[1]
                eqn_ = data[2:]

                eqn = []

                i = 0
                while True:
                    factor = float(eqn_[i])
                    i += 1
                    case = eqn_[i]
                    tmp = subcases[case]
                    eqn.append([factor, tmp[0], tmp[1]])
                    i += 1
                    if i >= len(eqn_):
                        break

                loadcases[case_id] = [case_description, eqn]

        h5keys = list(h5files.keys())

        for key in h5keys:
            self.model_ids.append(key)

        data = OrderedDict()

        for key, item in iteritems(loadcases):
            tmp = deepcopy(item[1])

            for tmp2 in tmp:
                tmp2[1] = h5keys.index(tmp2[1])
                tmp2[2] = int(tmp2[2]) - 1

            data[key] = [item[0], tmp]

        filenames = [h5files[key] for key in h5keys]
        self._set_filenames(filenames)

        case_index = 0
        for key, item in iteritems(data):
            case = CollectorLoadcase()

            case.id = int(key)
            case.description = item[0]

            lcs = item[1]

            for i in range(len(lcs)):
                lc = lcs[i]
                case.add_subcase(lc[1], lc[2], lc[0])

            self.cases.append(case)
            case_index += 1

        self.prepare_cases()
        # self.check_file_indices()

    def prepare_cases(self):
        self.by_case_description.clear()
        self.by_case_id.clear()

        for i in range(len(self.cases)):
            case = self.cases[i]
            self.by_case_id[case.id] = i
            self.by_case_description[case.description] = i
            
    def set_selection_by_case_id(self, case_ids):        
        loadcase_selection = []
        
        for case_id in case_ids:
            loadcase_selection.append(self.by_case_id[case_id])
            
        del self.loadcase_selection[:]
        self.loadcase_selection.extend(list(sorted(loadcase_selection)))

        self._model_selections()
            
    def set_section_by_case_description(self, case_descriptions):
        loadcase_selection = []

        for case_descr in case_descriptions:
            loadcase_selection.append(self.by_case_description[case_descr])

        del self.loadcase_selection[:]
        self.loadcase_selection.extend(list(sorted(loadcase_selection)))
        
        self._model_selections()
        
    def _model_selections(self):
        model_selections = []
        _model_selections = []
        
        for i in range(len(self.model_ids)):
            model_selections.append([])
            _model_selections.append({})
            
        for lc in self.loadcase_selection:
            subcases = self.cases[lc].model_subcases()
            
            for file_index, file_lc in subcases:
                model_selections[file_index].append(file_lc)

        del self.model_case_selection[:]
        
        for i in range(len(model_selections)):
            self.model_case_selection.append(list(sorted(model_selections[i])))
            
            model_selection = self.model_case_selection[i]
            
            for j in range(len(model_selection)):
                _model_selections[i][model_selection[j]] = j
                
        del self._model_case_selection[:]
        self._model_case_selection.extend(_model_selections)

    def check_file_indices(self):
        for case in self.cases:
            case.check_file_index()
