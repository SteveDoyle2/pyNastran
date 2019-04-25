from __future__ import print_function, absolute_import

from six.moves import range

from h5Nastran.h5nastrannode import H5NastranNode
from ..result_table import ResultTable, TableDef, DataGetter


class Nodal(H5NastranNode):
    def __init__(self, h5n, result):
        self._h5n = h5n
        self._result = result

        self.applied_loads = AppliedLoad(self._h5n, self)
        self.displacement = Displacement(self._h5n, self)
        self.grid_force = GridForce(self._h5n, self)
        self.mpc_force = MPCForce(self._h5n, self)
        self.spc_force = SPCForce(self._h5n, self)

    def path(self):
        return self._h5n.path() + ['NODAL']


########################################################################################################################


class Displacement(ResultTable):
    result_type = 'DISPLACEMENTS REAL'
    table_def = TableDef.create('/NASTRAN/RESULT/NODAL/DISPLACEMENT', result_type)


########################################################################################################################


def _validator(data):
    if data[1] == b'':
        data[1] = 0

    # this is done so that when doing grid point force summation, only nids and eids need
    # to be considered
    eids = {
        b'F-OF-SPC': -1,
        b'F-OF-MPC': -2,
        b'APP-LOAD': -3,
        b'*TOTALS*': -4
    }

    data[1] = eids.get(data[2].strip(), data[1])

    return data


class GridForce(ResultTable):
    result_type = 'GRID POINT FORCE BALANCE REAL'
    table_def = TableDef.create('/NASTRAN/RESULT/NODAL/GRID_FORCE', result_type,
                                indices=DataGetter(indices=[0, 2, 3, 5, 6, 7, 8, 9, 10]),
                                validator=_validator
                                )
    
    def search(self, data_ids, domains=(), convert_to_basic=False):
        _data = super(GridForce, self).search(data_ids, domains)
        
        if not convert_to_basic:
            return _data
        
        _nid = _data.loc[:, 'ID']
        _f1 = _data.loc[:, 'F1']
        _f2 = _data.loc[:, 'F2']
        _f3 = _data.loc[:, 'F3']
        _m1 = _data.loc[:, 'M1']
        _m2 = _data.loc[:, 'M2']
        _m3 = _data.loc[:, 'M3']

        data = _data.copy()
        f1 = data.loc[:, 'F1']
        f2 = data.loc[:, 'F2']
        f3 = data.loc[:, 'F3']
        m1 = data.loc[:, 'M1']
        m2 = data.loc[:, 'M2']
        m3 = data.loc[:, 'M3']

        get_grid = self._h5n.nastran.input.node.grid.get_grid
        vector_to_basic = self._h5n.nastran.input.coordinate_system.h5n_transformation.vector_to_basic

        for i in range(_nid.shape[0]):
            cd = get_grid(_nid[i])[3]

            f = [_f1[i], _f2[i], _f3[i]]
            m =[_m1[i], _m2[i], _m3[i]]

            if cd != 0:
                f = vector_to_basic(f, cd)
                m = vector_to_basic(m, cd)

            f1[i] = f[0]
            f2[i] = f[1]
            f3[i] = f[2]
            m1[i] = m[0]
            m2[i] = m[1]
            m3[i] = m[2]
            
        return data
            

########################################################################################################################


class MPCForce(ResultTable):
    result_type = 'MPCF REAL'
    table_def = TableDef.create('/NASTRAN/RESULT/NODAL/MPC_FORCE', result_type)

########################################################################################################################


class AppliedLoad(ResultTable):
    result_type = 'OLOADS REAL'
    table_def = TableDef.create('/NASTRAN/RESULT/NODAL/APPLIED_LOAD', result_type)

########################################################################################################################


class SPCForce(ResultTable):
    result_type = 'SPCF REAL'
    table_def = TableDef.create('/NASTRAN/RESULT/NODAL/SPC_FORCE', result_type)

########################################################################################################################
