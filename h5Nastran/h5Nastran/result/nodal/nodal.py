from __future__ import print_function, absolute_import
from six import iteritems, itervalues
from six.moves import range

import tables
import numpy as np


from ..result_table import ResultTable, TableDef, DataGetter


class Nodal(object):
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
    """
    <dataset name="DISPLACEMENT" version="1">
        <field name="ID" type="integer" description="Grid identifier"/>
        <field name="X"  type="double" description="X component"/>
        <field name="Y"  type="double" description="Y component"/>
        <field name="Z"  type="double" description="Z component"/>
        <field name="RX" type="double" description="RX component"/>
        <field name="RY" type="double" description="RY component"/>
        <field name="RZ" type="double" description="RZ component"/>
        <field name="DOMAIN_ID" type="integer" description="Domain identifier"/>
    </dataset>
    """

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
    """
    <dataset name="GRID_FORCE">
        <field name="ID" type="integer"/>
        <field name="EID" type="integer"/>
        <field name="ELNAME" type="character" size="8"/>
        <field name="F1" type="double"/>
        <field name="F2" type="double"/>
        <field name="F3" type="double"/>
        <field name="M1" type="double"/>
        <field name="M2" type="double"/>
        <field name="M3" type="double"/>
        <field name="DOMAIN_ID" type="integer"/>
    </dataset>
    """

    result_type = 'GRID POINT FORCE BALANCE REAL'
    table_def = TableDef.create('/NASTRAN/RESULT/NODAL/GRID_FORCE', result_type,
                                indices=DataGetter(indices=[0, 2, 3, 5, 6, 7, 8, 9, 10]),
                                validator=_validator
                                )

########################################################################################################################


class MPCForce(ResultTable):
    """
    <dataset name="MPC_FORCE" version="1">
        <field name="ID" type="integer" description="Grid identifier"/>
        <field name="X"  type="double" description="X component"/>
        <field name="Y"  type="double" description="Y component"/>
        <field name="Z"  type="double" description="Z component"/>
        <field name="RX" type="double" description="RX component"/>
        <field name="RY" type="double" description="RY component"/>
        <field name="RZ" type="double" description="RZ component"/>
        <field name="DOMAIN_ID" type="integer" description="Domain identifier"/>
    </dataset>
    """

    result_type = 'MPCF REAL'
    table_def = TableDef.create('/NASTRAN/RESULT/NODAL/MPC_FORCE', result_type)

########################################################################################################################


class AppliedLoad(ResultTable):
    """
    <dataset name="APPLIED_LOAD" version="1">
        <field name="ID" type="integer" description="Grid identifier"/>
        <field name="X"  type="double" description="X component"/>
        <field name="Y"  type="double" description="Y component"/>
        <field name="Z"  type="double" description="Z component"/>
        <field name="RX" type="double" description="RX component"/>
        <field name="RY" type="double" description="RY component"/>
        <field name="RZ" type="double" description="RZ component"/>
        <field name="DOMAIN_ID" type="integer" description="Domain identifier"/>
    </dataset>
    """

    result_type = 'OLOADS REAL'
    table_def = TableDef.create('/NASTRAN/RESULT/NODAL/APPLIED_LOAD', result_type)

########################################################################################################################


class SPCForce(ResultTable):
    """
    <dataset name="SPC_FORCE" version="1">
        <field name="ID" type="integer" description="Grid identifier"/>
        <field name="X"  type="double" description="X component"/>
        <field name="Y"  type="double" description="Y component"/>
        <field name="Z"  type="double" description="Z component"/>
        <field name="RX" type="double" description="RX component"/>
        <field name="RY" type="double" description="RY component"/>
        <field name="RZ" type="double" description="RZ component"/>
        <field name="DOMAIN_ID" type="integer" description="Domain identifier"/>
    </dataset>
    """

    result_type = 'SPCF REAL'
    table_def = TableDef.create('/NASTRAN/RESULT/NODAL/SPC_FORCE', result_type)

########################################################################################################################