from __future__ import print_function, absolute_import
from six import iteritems
from six.moves import range

import numpy as np
import pandas as pd
from typing import List, Tuple, Iterable, Union, Dict

from h5Nastran.result.result_table import ResultTableData


_Vector = Tuple[float, float, float]
_Matrix = Tuple[_Vector, _Vector, _Vector]

Vector = Union[_Vector, np.ndarray]
Matrix = Union[_Matrix, np.ndarray]


gpf_sum_type = np.dtype(
    [
        ('DETAIL', 'U32'),
        ('LOADCASE', np.int64),
        ('F1', np.float64),
        ('F2', np.float64),
        ('F3', np.float64),
        ('M1', np.float64),
        ('M2', np.float64),
        ('M3', np.float64)
    ]
)


class GridPointForceSummationData(ResultTableData):
    data_cols = pd.Index(['F1', 'F2', 'F3', 'M1', 'M2', 'M3'])
    data_group_by = ['DETAIL', 'LOADCASE']


class GridPointForceSummationCalculator(object):
    def __init__(self):
        self.data = None
        """:type: np.ndarray"""

        self.grid_pos = {}  # type: Dict[int, Vector]

        self._elements = {}
        self._nodes = {}
        self._cases = {}

        self._loadcases = []

    def set_data(self, data, grid_pos):
        self.data = data
        self.grid_pos = grid_pos
        self._build()

    def _build(self):
        self._elements.clear()
        self._nodes.clear()
        self._cases.clear()

        eids = self.data['EID']
        nids = self.data['ID']
        cids = self.data['DOMAIN_ID']

        elements = self._elements
        nodes = self._nodes
        cases = self._cases

        lcs = set()

        for i in range(eids.size):
            try:
                elements[eids[i]].append(i)
            except KeyError:
                elements[eids[i]] = [i]

            try:
                nodes[nids[i]].append(i)
            except KeyError:
                nodes[nids[i]] = [i]

            try:
                cases[cids[i]].append(i)
            except KeyError:
                cases[cids[i]] = [i]

            lcs.add(cids[i])

        del self._loadcases[:]
        self._loadcases.extend(sorted(lcs))

    def _indices(self, elements, nodes, loadcases=()):
        ei = []

        for eid in elements:
            try:
                ei.extend(self._elements[eid])
            except KeyError:
                pass

        ni = []

        for nid in nodes:
            try:
                ni.extend(self._nodes[nid])
            except KeyError:
                pass

        lci = []

        for lcid in loadcases:
            try:
                lci.extend(self._cases[lcid])
            except KeyError:
                pass

        ei = set(ei)
        ni = set(ni)

        if len(loadcases) == 0:
            return sorted(ei.intersection(ni))
        else:
            lci = set(lci)
            return sorted(ei.intersection(ni).intersection(lci))

    def sum(self, detail_id, nodes, elements, refpoint, coord, loadcases=(), load_factors=None):
        # type: (str, Iterable[int], Iterable[int], Vector, Matrix, Iterable[int], Vector) -> GridPointForceSummationData
        
        loadcases = sorted(loadcases)

        if len(loadcases) == 0:
            loadcases = self._loadcases

        result = np.zeros(len(loadcases), dtype=gpf_sum_type)

        _indices = set(self._indices(elements, nodes))

        if len(_indices) == 0:
            result['DETAIL'] = detail_id
            result['LDCASE'] = loadcases
            return GridPointForceSummationData.from_records(result)

        for i in range(len(loadcases)):
            lc = loadcases[i]

            casei = self._cases[lc]

            indices = sorted(_indices.intersection(casei))

            if len(indices) == 0:
                tmp = result[i]
                tmp[0] = detail_id
                tmp[1] = lc
                continue

            data = self.data.take(indices)
            result[i] = self._sum(data, detail_id, lc, refpoint, coord)

        if load_factors is not None:
            result['F1'] *= load_factors[0]
            result['F2'] *= load_factors[1]
            result['F3'] *= load_factors[2]
            result['M1'] *= load_factors[3]
            result['M2'] *= load_factors[4]
            result['M3'] *= load_factors[5]

        # print('numpy', result.dtype)

        result = GridPointForceSummationData.from_records(result)

        # print('dataframe', result.dtypes)

        return result

    def _sum(self, data, detail_id, lcid, refpoint, coord):

        fx = 0.
        fy = 0.
        fz = 0.
        mx = 0.
        my = 0.
        mz = 0.

        nids = data['ID'].values
        FX = data['F1'].values
        FY = data['F2'].values
        FZ = data['F3'].values
        MX = data['M1'].values
        MY = data['M2'].values
        MZ = data['M3'].values

        # print(FX.keys())

        grid_pos = self.grid_pos

        for i in range(data.shape[0]):
            _fx = FX[i]
            _fy = FY[i]
            _fz = FZ[i]
            _mx = MX[i]
            _my = MY[i]
            _mz = MZ[i]

            fx += _fx
            fy += _fy
            fz += _fz

            nid = nids[i]

            pos = grid_pos[nid]

            mx += _mx + _fz * (pos[1] - refpoint[1]) - _fy * (pos[2] - refpoint[2])
            my += _my + _fx * (pos[2] - refpoint[2]) - _fz * (pos[0] - refpoint[0])
            mz += _mz + _fy * (pos[0] - refpoint[0]) - _fx * (pos[1] - refpoint[1])

        _fx = np.dot([fx, fy, fz], coord[0])
        _fy = np.dot([fx, fy, fz], coord[1])
        _fz = np.dot([fx, fy, fz], coord[2])
        _mx = np.dot([mx, my, mz], coord[0])
        _my = np.dot([mx, my, mz], coord[1])
        _mz = np.dot([mx, my, mz], coord[2])

        return detail_id, lcid, _fx, _fy, _fz, _mx, _my, _mz

# 
# from zlib import compress as compress_, decompress as decompress_
# 
# 
# def decompress(compressed_data):
#     return decompress_(compressed_data, -15)
# 
# 
# def compress(uncompressed_data, compression_level=6):
#     return compress_(uncompressed_data, compression_level)[2:-4]


