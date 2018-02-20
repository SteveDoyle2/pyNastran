from __future__ import print_function, absolute_import
from six.moves import range

from copy import deepcopy
from typing import Iterable, List, Tuple, Union, Dict, Any

from ...h5nastran import H5Nastran
from ..collector import H5NastranCollector
from .grid_point_force_summation_data import GridPointForceSummationCalculator, Vector, Matrix, GridPointForceSummationData

import numpy as np


class GridPointForceSummationRequest(object):
    def __init__(self, detail_id, nids, eids, refpoint, coordinate_system=None, use_spc_forces=False, use_mpc_forces=False, use_app_forces=False):
        self.detail_id = detail_id
        self.nids = list(nids)
        self.eids = list(eids)
        self.refpoint = refpoint  # type: Vector
        self.coordinate_system = coordinate_system  # type: Matrix

        self.use_spc_forces = use_spc_forces
        self.use_mpc_forces = use_mpc_forces
        self.use_app_forces = use_app_forces
        
        self.data = None  # type: GridPointForceSummationData

    def copy(self):
        return deepcopy(self)

    def run(self, calculator, loadcases=(), load_factors=None):
        # type: (GridPointForceSummationCalculator, Iterable[int], Vector) -> None

        if self.coordinate_system is None:
            self.coordinate_system = np.array([
                [1., 0., 0.],
                [0., 1., 0.],
                [0., 0., 1.]
            ])

        _eids = []
        if self.use_spc_forces:
            _eids.append(-1)
        if self.use_mpc_forces:
            _eids.append(-2)
        if self.use_app_forces:
            _eids.append(-3)

        eids = self.eids + _eids
        
        self.data = -1 * calculator.sum(self.detail_id, self.nids, eids, self.refpoint, self.coordinate_system, loadcases, load_factors)


class GridPointForceSummationRequestList(object):
    def __init__(self):
        self.requests = []  # type: List[GridPointForceSummationRequest]
        self.calculator = GridPointForceSummationCalculator()
        self.db = None  # type: H5Nastran
        self.domain_ids = []  # type: Iterable[int]
        self.grid_pos = None
        self.load_factors = None  # type: Vector

    def set_db(self, db):
        # type: (Union[H5Nastran, H5NastranCollector]) -> None
        self.db = db
        # self.grid_pos = self.db.input.node.grid.grid_in_basic_dict()

    def add_request(self, request):
        # type: (GridPointForceSummationRequest) -> None
        self.requests.append(request)

    def elements(self):
        # type: () -> List[int]
        elements = []
        for request in self.requests:
            elements.extend(request.eids)
        return sorted(list(set(elements)))

    def nodes(self):
        # type: () -> List[int]
        nodes = []
        for request in self.requests:
            nodes.extend(request.nids)
        return sorted(list(set(nodes)))

    def run(self):
        # type: () -> None
        if isinstance(self.db, H5Nastran):
            self._run_h5n()
        elif isinstance(self.db, H5NastranCollector):
            self._run_collector()

    def _run_h5n(self):
        all_nodes = self.nodes()
        grid_forces = self.db.result.nodal.grid_force.search(all_nodes, self.domain_ids, True)

        if self.grid_pos is None:
            self.grid_pos = self.db.input.node.grid.grid_in_basic_dict()

        self.calculator.set_data(grid_forces, self.grid_pos)

        for request in self.requests:
            request.run(self.calculator, load_factors=self.load_factors)

    def _run_collector(self):
        db = self.db  # type: H5NastranCollector
        db.open()

        if len(self.domain_ids) == 0:
            subcase_ids = list(range(1, len(db.cases)+1))
        else:
            subcase_ids = self.domain_ids

        db.set_selection_by_case_id(subcase_ids)

        results = []

        for i in range(len(self.requests)):
            results.append([])

        all_nodes = self.nodes()

        for i in range(len(db.h5n)):
            h5n = db.h5n[i]
            # domains = db.model_case_selection[i]
            # TODO: need to handle loadcase selection better
            grid_forces = h5n.result.nodal.grid_force.search(all_nodes, (), True)

            grid_pos = h5n.input.node.grid.grid_in_basic_dict()

            self.calculator.set_data(grid_forces, grid_pos)

            j = 0
            for request in self.requests:
                request.run(self.calculator, load_factors=self.load_factors)
                results[j].append(request.data)
                j += 1

        for i in range(len(results)):
            self.requests[i].data = db.combine(results[i])
            self.requests[i].data['DETAIL'] = self.requests[i].detail_id
