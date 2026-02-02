"""
defines:
 - TableObject
 - RealTableArray
 - ComplexTableArray

these are used by:
 - RealDisplacementArray
 - RealVelocityArray
 - RealAccelerationArray
 - RealEigenvaluesArray
 - RealSPCForcesArray
 - RealMPCForcesArray
 - RealAppliedLoadsArray

 - ComplexDisplacementArray
 - ComplexVelocityArray
 - ComplexAccelerationArray
 - ComplexEigenvaluesArray
 - ComplexSPCForcesArray
 - ComplexMPCForcesArray
 - ComplexAppliedLoadsArray

"""
from __future__ import annotations
import copy
from struct import Struct, pack
from itertools import count
import warnings
from typing import TextIO, Optional

import numpy as np

from pyNastran.bdf import MAX_32_BIT_INT
from pyNastran.op2.result_objects.op2_objects import (
    ScalarObject, get_sort_node_sizes, set_as_sort1,
    NULL_GRIDTYPE, SORT1_TABLES, SORT2_TABLES,
    # GRID_TYPE_INT_TO_STR,
    recast_gridtype_as_string,
)
from pyNastran.op2.result_objects.utils_pandas import build_dataframe_transient_header
from pyNastran.f06.f06_formatting import (
    write_floats_13e, write_floats_13e_long,
    write_imag_floats_13e, write_float_12e)
from pyNastran.op2.errors import SixtyFourBitError
from pyNastran.op2.op2_helper import polar_to_real_imag, real_imag_to_mag_phase
from pyNastran.op2.op2_interface.write_utils import set_table3_field, view_dtype, view_idtype_as_fdtype
from pyNastran.utils.numpy_utils import integer_types, float_types, integer_float_types
from pyNastran.op2.writer.utils import fix_table3_types
from pyNastran.op2.tables.oes_stressStrain.real.oes_objects import (
    set_static_case, set_modal_case, set_transient_case,
    set_freq_case, set_complex_modes_case)
from pyNastran.op2.result_objects.op2_objects import combination_inplace

table_name_to_table_code = {
    # displacement (msc/nx)
    'OUGV1' : 1,
    'BOUGV1' : 1,

    # eigenvector
    # 7

    # velocity
    'OVG1': 10,

    # acceleration
    'OAG1': 11,

    # load vector (msc/nx)
    'OPG1' : 2,
    'BOPG1' : 2,
    #'BOPHIG1' : 5, # ???

    # spc/mpc forces
    'OQG1' : 3,
}
def append_sort1_sort2(data1, data2, to_sort1=True):
    """
    data1 : (ntimes, nnids, 6)
    data2 : (nnids, ntimes, 6)
    """
    assert len(data1.shape) == 3, data1.shape
    assert len(data2.shape) == 3, data2.shape
    ntimes1, nnids1 = data1.shape[:2]
    nnids2, ntimes2 = data2.shape[:2]
    unused_ntimes = ntimes1 + ntimes2
    unused_nnids = nnids1 + nnids2
    assert ntimes1 == ntimes2
    if to_sort1:
        out = np.hstack([
            data1,
            np.swapaxes(data2, 0, 1),])
    else:
        out = np.hstack([
            np.swapaxes(data1, 0, 1),
            data2,])
    return out

def oug_data_code(table_name,
                  is_real: bool=False, is_complex:bool=False,
                  is_sort1=True, is_random=False,
                  random_code=0, title='', subtitle='', label='', is_msc=True):
    sort1_sort_bit = 0 if is_sort1 else 1
    random_sort_bit = 1 if is_random else 0
    sort_method = 1 if is_sort1 else 2
    #if format_code == 1:
        #format_word = "Real"
    #elif format_code == 2:
        #format_word = "Real/Imaginary"
    #elif format_code == 3:
        #format_word = "Magnitude/Phase"
    #DEVICE_CODE_MAP = {
        #1 : "Print",
        #2 : "Plot",
        #3 : "Print and Plot",
        #4 : "Punch",
        #5 : "Print and Punch",
        #6 : "Plot and Punch",
        #7 : "Print, Plot, and Punch",
    #}

    if is_real:
        num_wide = 8
    elif is_complex:
        num_wide = 14

    table_code = table_name_to_table_code[table_name]
    sort_code = 1 # TODO: what should this be???

    #table_code = tCode % 1000
    #sort_code = tCode // 1000
    tCode = table_code * 1000 + sort_code

    device_code = 2  # Plot
    data_code = {
        'nonlinear_factor': None,
        'sort_bits': [0, sort1_sort_bit, random_sort_bit], # real, sort1, random
        'sort_method' : sort_method,
        'is_msc': is_msc,
        'format_code': 1, # real
        'table_code': table_code,
        'tCode': tCode,
        'table_name': table_name, ## TODO: should this be a string?
        'device_code' : device_code,
        'random_code' : random_code,
        'thermal': 0,
        'title' : title,
        'subtitle': subtitle,
        'label': label,
        'num_wide': num_wide,
    }
    return data_code


class TableArray(ScalarObject):  # displacement style table
    """
    Base class for:
     - RealTableArray
     - ComplexTableArray
    """
    def __init__(self, data_code, is_sort1, isubcase, dt):
        self.nonlinear_factor = np.nan
        #self.table_name = None
        #self.approach_code = None
        #self.analysis_code = None

        # no double inheritance
        ScalarObject.__init__(self, data_code, isubcase, apply_data_code=True)
        str(self.is_sort1)
        str(self.is_sort2)
        #self.dt = dt

        #self.code = [self.format_code, self.sort_code, self.s_code]

        #self.ntimes = 0  # or frequency/mode
        self.ntotal = 0
        self._nnodes = 0  # result specific

    def __eq__(self, table):  # pragma: no cover
        return self.assert_equal(table)

    def assert_equal(self, table, rtol=1.e-5, atol=1.e-8):
        self._eq_header(table)
        assert self.is_sort1 == table.is_sort1
        #print(self.node_gridtype)
        #print(table.node_gridtype)
        if not np.array_equal(self.node_gridtype, table.node_gridtype):
            assert self.node_gridtype.shape == table.node_gridtype.shape, 'shape=%s table.shape=%s' % (self.node_gridtype.shape, table.node_gridtype.shape)
            msg = 'table_name=%r class_name=%s\n' % (self.table_name, self.__class__.__name__)
            msg += '%s\n' % str(self.code_information())
            msg += 'nid_gridtype:\n'
            msg += 'gridtype.shape=%s table.gridtype.shape=%s\n' % (str(self.node_gridtype.shape), str(table.node_gridtype.shape))
            for (nid, grid_type), (nid2, grid_type2) in zip(self.node_gridtype, table.node_gridtype):
                msg += '(%s, %s)    (%s, %s)\n' % (nid, grid_type, nid2, grid_type2)
            print(msg)
            raise ValueError(msg)
        if not np.array_equal(self.data, table.data):
            msg = 'table_name=%r class_name=%s\n' % (self.table_name, self.__class__.__name__)
            msg += '%s\n' % str(self.code_information())
            ntimes = self.data.shape[0]

            i = 0
            atols = []
            rtols = []
            if self.is_sort1:
                for itime in range(ntimes):
                    msg += '(nid, gridtype); itime=%s\n' % itime
                    msg += '(tx, ty, tz, rx, ry, rz)\n'
                    for inid, nid_gridtype, in enumerate(self.node_gridtype):
                        (nid, grid_type) = nid_gridtype
                        t1 = self.data[itime, inid, :]
                        t2 = table.data[itime, inid, :]
                        if not np.allclose(t1, t2, rtol=rtol, atol=atol):
                        #if not np.array_equal(t1, t2):
                            inonzero = np.where(t1 != 0.)[0]
                            atoli = np.abs(t2 - t1).max()
                            rtoli = np.abs(t2[inonzero] / t1[inonzero]).max()

                            (tx1, ty1, tz1, rx1, ry1, rz1) = t1
                            (tx2, ty2, tz2, rx2, ry2, rz2) = t2
                            msg += '(%s, %s)\n  (%s, %s, %s, %s, %s, %s)\n  (%s, %s, %s, %s, %s, %s)\n' % (
                                nid, grid_type,
                                tx1, ty1, tz1, rx1, ry1, rz1,
                                tx2, ty2, tz2, rx2, ry2, rz2)
                            i += 1
                            atols.append(atoli)
                            rtols.append(rtoli)
                        if i > 10:
                            msg += 'atol.max() = %s\n' % max(atols)
                            msg += 'rtol.max() = %s\n' % max(rtols)
                            print(msg)
                            raise ValueError(msg)
            else:
                raise NotImplementedError(self.is_sort2)
            if i > 0:
                msg += 'atol.max() = %s\n' % max(atols)
                msg += 'rtol.max() = %s\n' % max(rtols)
                print(msg)
                raise ValueError(msg)
        return True

    def combine(self, result, is_sort1=True):
        #print("combine; result=%s" % result)
        assert self.is_sort1 != result.is_sort1
        assert self.nonlinear_factor is not None
        assert result.nonlinear_factor is not None
        # self.ntimes += result.ntimes
        self.ntotal += result.data.shape[0]
        self.data = append_sort1_sort2(self.data, result.data)
        #print(self._times)
        #print(result._times)
        # self._times = hstack([self._times, result._times])
        self.node_gridtype = np.vstack([self.node_gridtype, result.node_gridtype])
        #print('%s' % ''.join(self.get_stats()))

    def _get_msgs(self, is_mag_phase):
        raise NotImplementedError()

    def data_type(self):
        raise NotImplementedError()

    def slice_data(self, slice_nodes: np.ndarray) -> int:
        assert slice_nodes is not None, slice_nodes
        nodes = self.node_gridtype[:, 0]
        common_nodes = np.intersect1d(nodes, slice_nodes)
        icommon = np.searchsorted(nodes, common_nodes)
        ntotal = len(icommon)
        self.node_gridtype = self.node_gridtype[icommon, :]
        self.data = self.data[:, icommon, :]
        self.ntotal = ntotal
        return ntotal

    def get_stats(self, short: bool=False) -> list[str]:
        if not self.is_built:
            return [
                '<%s>; table_name=%r\n' % (self.__class__.__name__, self.table_name),
                f'  ntimes: {self.ntimes:d}\n',
                f'  ntotal: {self.ntotal:d}\n',
            ]
        #ngrids = len(self.gridTypes)
        if short:
            return self._get_stats_short()
        msg = []

        unused_ntimesi, ntotal = self.data.shape[:2]
        ntimes = len(self._times)
        nnodes = self.node_gridtype.shape[0]

        nmajor = self.ntimes
        nminor = self.ntotal
        if self.is_sort1:
            assert nmajor == ntimes, 'ntimes=%s expected=%s' % (nmajor, ntimes)
            assert nminor == ntotal, 'ntotal=%s expected=%s' % (nminor, nnodes)
        else:
            if not nmajor == nnodes:
                msgi = 'nnodes=%s expected=%s' % (nmajor, nnodes)
                warnings.warn(msgi)
                msg.append('  WARNING: ' + msgi + '\n')
            assert nminor == ntotal, 'ntotal=%s expected=%s' % (nminor, ntimes)

        msg.append('  isubcase = %s\n' % self.isubcase)
        if self.nonlinear_factor not in (None, np.nan):  # transient
            msg.append('  type=%s ntimes=%s nnodes=%s, table_name=%s\n'
                       % (self.__class__.__name__, ntimes, nnodes, self.table_name))
        else:
            msg.append('  type=%s nnodes=%s, table_name=%s\n'
                       % (self.__class__.__name__, nnodes, self.table_name))
        headers = ', '.join(self._get_headers())
        #msg.append('  data: [%s] shape=%s dtype=%s\n'
                   #% (headers, [int(i) for i in self.data.shape], self.data.dtype))
        msg.append('  data: [%s] shape=%s dtype=%s\n'
                   % (headers,
                      [int(i) for i in self.data.shape], self.data.dtype))
        msg.append(f'  node_gridtype.shape = {self.node_gridtype.shape}\n')
        #msg.append('  gridTypes\n  ')
        msg += self.get_data_code()
        return msg

    @property
    def headers(self) -> list[str]:
        return ['t1', 't2', 't3', 'r1', 'r2', 'r3']

    def _reset_indices(self) -> None:
        self.itotal = 0

    def build(self):
        """sizes the vectorized attributes of the TableArray"""
        #print('_nnodes=%s ntimes=%s sort1?=%s ntotal=%s -> _nnodes=%s' % (self._nnodes, self.ntimes, self.is_sort1,
                                                                          #self.ntotal, self._nnodes // self.ntimes))

        # we have a SORT1 data array that will be (ntimes, nnodes, 6)
        # we start by sizing the total number of entries (_nnodes = ntimes * nnodes)
        # we also keep track of the number of times
        # then we compute nnodes
        #
        # for sort1, we just use what was discussed above
        # for sort2, we flip nnodes and ntimes
        #
        # note that in both cases, ntotal is the major dimension:
        #  - SORT1 - ntimes
        #  - SORT2 - nnodes
        #print('ntotal=%s ntimes=%s _nnodes=%s' % (self.ntotal, self.ntimes, self._nnodes))
        #print('ntotal=%s ntimes=%s _nnodes=%s\n' % (self.ntotal, self.ntimes, self._nnodes))
        #if self.ntimes > 1000:
        #    raise RuntimeError(self.ntimes)

        self.itime = 0
        self.itotal = 0

        ntimes, nnodes, ntotal = get_sort_node_sizes(self, debug=False)
        #if self.is_sort1:
            #ntimes = self.ntimes
            #nnodes = self.ntotal
            #ntotal = self.ntotal
            #nx = ntimes
            #ny = nnodes
            #print("SORT1 ntimes=%s nnodes=%s" % (ntimes, nnodes))
        #elif self.is_sort2:
            ## flip this to sort1
            #ntimes = self.ntotal
            #nnodes = self.ntimes
            #ntotal = nnodes
            #nx = ntimes
            #ny = nnodes
            #print("***SORT2 ntotal=%s nnodes=%s ntimes=%s" % (ntotal, nnodes, ntimes))
        #else:
            #raise RuntimeError('expected sort1/sort2\n%s' % self.code_information())
        #self.build_data(ntimes, nnodes, ntotal, self._times_dtype)
        #print(self.class_name, self.analysis_fmt)
        self.build_data(ntimes, nnodes, ntotal, self.analysis_fmt)

    def build_data(self, ntimes, nnodes, ntotal, float_fmt: str):
        """actually performs the build step"""
        self.ntimes = ntimes
        self._nnodes = nnodes
        self.ntotal = ntotal

        _times = np.zeros(ntimes, dtype=float_fmt)
        int_fmt = 'int32' if self.size == 4 else 'int64'
        node_gridtype = np.zeros((nnodes, 2), dtype=int_fmt)

        #[t1, t2, t3, r1, r2, r3]
        data = np.zeros((ntimes, nnodes, 6), self.data_type())
        if self.load_as_h5:
            group = self._get_result_group()
            self._times = group.create_dataset('_times', data=_times)
            self.node_gridtype = group.create_dataset('node_gridtype', data=node_gridtype)
            self.data = group.create_dataset('data', data=data)
        else:
            self._times = _times
            self.node_gridtype = node_gridtype
            self.data = data
        #print('ntimes=%s nnodes=%s; nx=%s ny=%s; ntotal=%s' % (
            #ntimes, nnodes, nx, ny, self.ntotal))

    def build_dataframe(self):
        """creates a pandas dataframe

        works: 0.24.2
        broken: 0.25.0
        """
        import pandas as pd

        headers = self.get_headers()
        # headers = [0, 1, 2, 3, 4, 5]
        # node_gridtype = [self.node_gridtype[:, 0], self.gridtype_str]
        # letter_dims = [
        #     ('G', 6),
        #     ('E', 1),
        #     ('S', 1),
        #     ('H', 6),
        #     ('L', 6),
        # ]
        ntimes, nnodes = self.data.shape[:2]

        # ugridtype_str = np.unique(self.gridtype_str)
        if self.nonlinear_factor not in (None, np.nan):
            # if not self.is_sort1:
            #     print("skipping %s because it's not SORT1" % self.class_name)
            #     return
            column_names, column_values = build_dataframe_transient_header(self)
            #if is_v25:
            #  we start out like this...
            #
            # Mode                             1                 2                   3
            # EigenvalueReal               -0.0              -0.0                -0.0
            # EigenvalueImag          -0.463393          0.463393           -1.705689
            # Damping                        0.0               0.0                 0.0
            # NodeID Type Item
            # 1      G    t1      (0.6558146+0j)    (0.6558146+0j)       (1.034078+0j)
            #             t2                  0j                0j                  0j
            #             t3                  0j                0j                  0j
            #             r1                  0j                0j                  0j
            #             r2                  0j                0j                  0j
            #             r3                  0j                0j                  0j
            #  ...
            #
            # then we call pandas_extract_rows to make it this...
            #
            # Mode                        1              2              3
            # EigenvalueReal           -0.0           -0.0           -0.0
            # EigenvalueImag      -0.463393       0.463393      -1.705689
            # Damping                   0.0            0.0            0.0
            # NodeID Item
            # 1      t1       0.655815+0.0j  0.655815+0.0j  1.034078+0.0j
            #        t2            0.0+0.0j       0.0+0.0j       0.0+0.0j
            #        t3            0.0+0.0j       0.0+0.0j       0.0+0.0j
            #        r1            0.0+0.0j       0.0+0.0j       0.0+0.0j
            #        r2            0.0+0.0j       0.0+0.0j       0.0+0.0j
            #        r3            0.0+0.0j       0.0+0.0j       0.0+0.0j
            # 2      t1       0.999141+0.0j  0.999141+0.0j -0.282216+0.0j
            #        t2            0.0+0.0j       0.0+0.0j       0.0+0.0j
            #        t3            0.0+0.0j       0.0+0.0j       0.0+0.0j
            #        r1            0.0+0.0j       0.0+0.0j       0.0+0.0j
            #        r2            0.0+0.0j       0.0+0.0j       0.0+0.0j
            #        r3            0.0+0.0j       0.0+0.0j       0.0+0.0j
            # 1001   S        0.000859+0.0j  0.000859+0.0j -0.003323+0.0j
            try:
                columns = pd.MultiIndex.from_arrays(column_values, names=column_names)
            except ValueError:
                msg = ''
                for name, values in zip(column_names, column_values):
                    msg += f'{name!r}: {values}; n={len(values)}\n'
                raise ValueError(msg)

            gridtype_str = self.gridtype_str
            ugridtype_str = np.unique(gridtype_str)
            if len(ugridtype_str) == 1 and gridtype_str[0] in ['S', 'M', 'E']:
                nnodes = self.node_gridtype.shape[0]
                node_gridtype = [self.node_gridtype[:, 0], [gridtype_str[0]] * nnodes]

                names = ['NodeID', 'Item']
                index = pd.MultiIndex.from_arrays(node_gridtype, names=names)
                A = self.data[:, :, 0].T
                data_frame = pd.DataFrame(A, columns=columns, index=index)
            else:
                node_gridtype_item = []
                node_ids = self.node_gridtype[:, 0]
                for nid, gridtype in zip(node_ids, gridtype_str):
                    node_gridtype_item.extend([[nid, gridtype, 't1']])
                    node_gridtype_item.extend([[nid, gridtype, 't2']])
                    node_gridtype_item.extend([[nid, gridtype, 't3']])
                    node_gridtype_item.extend([[nid, gridtype, 'r1']])
                    node_gridtype_item.extend([[nid, gridtype, 'r2']])
                    node_gridtype_item.extend([[nid, gridtype, 'r3']])

                names = ['NodeID', 'Type', 'Item']
                index = pd.MultiIndex.from_tuples(node_gridtype_item, names=names)
                A = self.data.reshape(ntimes, nnodes*6).T
                try:
                    data_frame = pd.DataFrame(A, columns=columns, index=index)
                except ValueError:  # pragma: no cover
                    print(f'data.shape={self.data.shape} A.shape={A.shape} '
                          f'ntimes={ntimes} nnodes*6={nnodes*6} ngrids={len(node_ids)}\n'
                          f'column_names={column_names} column_values={column_values} _times={self._times}')
                    raise
                #print(data_frame.to_string())
                data_frame = pandas_extract_rows(data_frame, ugridtype_str, ['NodeID', 'Item'])

            #elif is_v25 and 0:  # pragma: no cover
                #                                                                           t1        t2
                # itime Mode  EigenvalueReal EigenvalueImag Damping  NodeID  Type
                # 0      1      G    -0.0         -0.463393     0.0  1           0.655815+0.0j  0.0+0.0j
                #               G                                    2           0.999141+0.0j  0.0+0.0j
                #               G                                    3                1.0+0.0j  0.0+0.0j
                #               S                                    1001        0.000859+0.0j  0.0+0.0j
                # 1      1      G    -0.0          0.463393     0.0  2           0.655815+0.0j  0.0+0.0j
                #               G                                    2           0.999141+0.0j  0.0+0.0j
                #               G                                    3                1.0+0.0j  0.0+0.0j
                #               S                                    1001        0.000859+0.0j  0.0+0.0j
                # 2      1      G    -0.0        -1.705689      0.0  3           1.034078+0.0j  0.0+0.0j
                #               G                                    2          -0.282216+0.0j  0.0+0.0j
                #               G                                    3          -0.285539+0.0j  0.0+0.0j
                #               S                                    1001       -0.003323+0.0j  0.0+0.0j
                #time_node_gridtype = []
                #from itertools import count
                #for itime in range(ntimes):
                    #column_values2 = [column_value[itime] for column_value in column_values]
                    #for nid, gridtype in zip(self.node_gridtype[:, 0], self.gridtype_str):
                        #time_node_gridtype.append([itime] + column_values2 + [nid, gridtype])

                #names = ['itime'] + column_names + ['NodeID', 'Type']
                #index = pd.MultiIndex.from_tuples(time_node_gridtype, names=names)
                #A = self.data.reshape(ntimes*nnodes, 6)
                #data_frame = pd.DataFrame(A, columns=headers, index=index)
                ##print(self.data_frame.index.names)
                ##data_frame = pandas_extract_rows(self.data_frame, ugridtype_str)
                #print(data_frame)
            #elif is_v25 and 0:  # pragma: no cover
                #node_gridtype2 = []
                #NodeID Type             t1        t2        ...
                #1      G     0.655815+0.0j  0.0+0.0j  0.0+0.0j  0.0+0.0j  0.0+0.0j  0.0+0.0j
                #2      G     0.999141+0.0j  0.0+0.0j  0.0+0.0j  0.0+0.0j  0.0+0.0j  0.0+0.0j
                #3      G          1.0+0.0j  0.0+0.0j  0.0+0.0j  0.0+0.0j  0.0+0.0j  0.0+0.0j
                #1001   S     0.000859+0.0j  0.0+0.0j  0.0+0.0j  0.0+0.0j  0.0+0.0j  0.0+0.0j
                #1      G     0.655815+0.0j  0.0+0.0j  0.0+0.0j  0.0+0.0j  0.0+0.0j  0.0+0.0j
                #2      G     0.999141+0.0j  0.0+0.0j  0.0+0.0j  0.0+0.0j  0.0+0.0j  0.0+0.0j
                #3      G          1.0+0.0j  0.0+0.0j  0.0+0.0j  0.0+0.0j  0.0+0.0j  0.0+0.0j
                #1001   S     0.000859+0.0j  0.0+0.0j  0.0+0.0j  0.0+0.0j  0.0+0.0j  0.0+0.0j
                #1      G     1.034078+0.0j  0.0+0.0j  0.0+0.0j  0.0+0.0j  0.0+0.0j  0.0+0.0j
                #2      G    -0.282216+0.0j  0.0+0.0j  0.0+0.0j  0.0+0.0j  0.0+0.0j  0.0+0.0j
                #3      G    -0.285539+0.0j  0.0+0.0j  0.0+0.0j  0.0+0.0j  0.0+0.0j  0.0+0.0j
                #1001   S    -0.003323+0.0j  0.0+0.0j  0.0+0.0j  0.0+0.0j  0.0+0.0j  0.0+0.0j
                #1      G     1.034078+0.0j  0.0+0.0j  0.0+0.0j  0.0+0.0j  0.0+0.0j  0.0+0.0j
                #2      G    -0.282216+0.0j  0.0+0.0j  0.0+0.0j  0.0+0.0j  0.0+0.0j  0.0+0.0j
                #3      G    -0.285539+0.0j  0.0+0.0j  0.0+0.0j  0.0+0.0j  0.0+0.0j  0.0+0.0j
                #1001   S    -0.003323+0.0j  0.0+0.0j  0.0+0.0j  0.0+0.0j  0.0+0.0j  0.0+0.0j
                #1      G    -0.001818+0.0j  0.0+0.0j  0.0+0.0j  0.0+0.0j  0.0+0.0j  0.0+0.0j
                #2      G    -0.124197+0.0j  0.0+0.0j  0.0+0.0j  0.0+0.0j  0.0+0.0j  0.0+0.0j
                #3      G     0.625574+0.0j  0.0+0.0j  0.0+0.0j  0.0+0.0j  0.0+0.0j  0.0+0.0j
                #1001   S     0.749771+0.0j  0.0+0.0j  0.0+0.0j  0.0+0.0j  0.0+0.0j  0.0+0.0j
                #1      G     0.001011+0.0j  0.0+0.0j  0.0+0.0j  0.0+0.0j  0.0+0.0j  0.0+0.0j
                #2      G    -0.200504+0.0j  0.0+0.0j  0.0+0.0j  0.0+0.0j  0.0+0.0j  0.0+0.0j
                #3      G          1.0+0.0j  0.0+0.0j  0.0+0.0j  0.0+0.0j  0.0+0.0j  0.0+0.0j
                #1001   S     1.200504+0.0j  0.0+0.0j  0.0+0.0j  0.0+0.0j  0.0+0.0j  0.0+0.0j
                #1      G     0.001011+0.0j  0.0+0.0j  0.0+0.0j  0.0+0.0j  0.0+0.0j  0.0+0.0j
                #2      G    -0.200504+0.0j  0.0+0.0j  0.0+0.0j  0.0+0.0j  0.0+0.0j  0.0+0.0j
                #3      G          1.0+0.0j  0.0+0.0j  0.0+0.0j  0.0+0.0j  0.0+0.0j  0.0+0.0j
                #1001   S     1.200504+0.0j  0.0+0.0j  0.0+0.0j  0.0+0.0j  0.0+0.0j  0.0+0.0j
                #for itime in range(ntimes):
                    #column_values2 = [column_value[itime] for column_value in column_values]
                    #for nid, gridtype in zip(self.node_gridtype[:, 0], self.gridtype_str):
                        #node_gridtype2.append([itime, nid, gridtype])

                #names = ['itime', 'NodeID', 'Type']
                #index = pd.MultiIndex.from_tuples(node_gridtype2, names=names)
                #A = self.data.reshape(ntimes*nnodes, 6)
                #data_frame = pd.DataFrame(A, columns=headers, index=index)
                ##print(data_frame.index.names)
                ##data_frame = pandas_extract_rows(data_frame, ugridtype_str)
                #print(data_frame)

            #elif is_v25 and 0:  # pragma: no cover
                #                t1        t2        t3        r1        r2        r3
                # 0   0.655815+0.0j  0.0+0.0j  0.0+0.0j  0.0+0.0j  0.0+0.0j  0.0+0.0j
                # 1   0.999141+0.0j  0.0+0.0j  0.0+0.0j  0.0+0.0j  0.0+0.0j  0.0+0.0j
                # 2        1.0+0.0j  0.0+0.0j  0.0+0.0j  0.0+0.0j  0.0+0.0j  0.0+0.0j
                # 3   0.000859+0.0j  0.0+0.0j  0.0+0.0j  0.0+0.0j  0.0+0.0j  0.0+0.0j
                # 4   0.655815+0.0j  0.0+0.0j  0.0+0.0j  0.0+0.0j  0.0+0.0j  0.0+0.0j
                # 5   0.999141+0.0j  0.0+0.0j  0.0+0.0j  0.0+0.0j  0.0+0.0j  0.0+0.0j
                # 6        1.0+0.0j  0.0+0.0j  0.0+0.0j  0.0+0.0j  0.0+0.0j  0.0+0.0j
                #index = pd.MultiIndex.from_arrays(node_gridtype, names=['NodeID', 'Type'])
                #A = self.data.reshape(nnodes*ntimes, 6)
                #data_frame = pd.DataFrame(A, columns=headers)
                ##data_frame = pd.DataFrame(A, columns=headers, index=index)  # doesn't work
                # doesn't turn into workable table
            #else:
                # old
                # Mode                             1                 2                   3
                # EigenvalueReal               -0.0              -0.0                -0.0
                # EigenvalueImag          -0.463393          0.463393           -1.705689
                # Damping                        0.0               0.0                 0.0
                # NodeID Type Item
                # 1      G    t1      (0.6558146+0j)    (0.6558146+0j)       (1.034078+0j)
                #             t2                  0j                0j                  0j
                #             t3                  0j                0j                  0j
                #             r1                  0j                0j                  0j
                #             r2                  0j                0j                  0j
                #             r3                  0j                0j                  0j

                #   mode    1      2    3
                #   freq    1.0  2.0  3.0
                # nodeid
                #  1  item  1.0  2.0  3.0
                #     t1    etc.
                #     t2
                #     t3
                #     ...
                #  2
                #     t1
                #     t2
                #     t3
                #     ...

                #data_frame = pd.Panel(self.data, items=column_values,
                                      #major_axis=node_gridtype, minor_axis=headers).to_frame()  # to_xarray()
                #data_frame.columns.names = column_names
                #data_frame.index.names = ['NodeID', 'Type', 'Item']
                #print(column_names)
                #print(data_frame)
                #print(self.data_frame.index.names)
                #data_frame = pandas_extract_rows(data_frame, ugridtype_str, ['NodeID', 'Item'])
            #print(data_frame)

        else:
            #self.data_frame = pd.Panel(self.data[0, :, :], major_axis=node_gridtype, minor_axis=headers).to_frame()
            #self.data_frame.columns.names = ['Static']
            #self.data_frame.index.names = ['NodeID', 'Type', 'Item']

            #    NodeID Type      t1           t2            t3   r1   r2   r3
            # 0        1    G     0.0     0.000000  0.000000e+00  0.0  0.0  0.0
            # 1        2    G     0.0     0.000000  0.000000e+00  0.0  0.0  0.0
            # 2        3    G     0.0     0.000000  0.000000e+00  0.0  0.0  0.0
            #self.data_frame = pd.DataFrame(self.data[0], columns=headers, index=self.node_gridtype)

            nnode = self.node_gridtype.shape[0]
            # gridtype_str1 = self.gridtype_str
            gridtype_str = get_gridtype_str(self, self.node_gridtype)
            # print('gridtype_str1 = ', gridtype_str1, gridtype_str1.shape)
            # print('gridtype_str2 = ', gridtype_str2, gridtype_str2.shape)
            data_dict = {
                'NodeID': self.node_gridtype[:, 0],
                'Type': gridtype_str,
            }
            for i, header in enumerate(headers):
                datai = self.data[0, :, i]
                assert len(datai) == nnode, (len(datai), nnode)
                data_dict[header] = datai
            data_frame = pd.DataFrame(data_dict)
        #print(data_frame)
        self.data_frame = data_frame

    def finalize(self) -> None:
        """
        Calls any OP2 objects that need to do any post matrix calcs
        """
        self.set_as_sort1()
        #self.gridtype_str2 = np.chararray(nnodes, unicode=True)
        self.gridtype_str = get_gridtype_str(self, self.node_gridtype)
        #print(self.gridtype_str, self.gridtype_str.dtype)
        #print(self.gridtype_str2, self.gridtype_str2.dtype)
        #asdf
        #del self.itotal, self.itime

    def set_as_sort1(self):
        """changes the table into SORT1"""
        #print(self.class_name, self._times.dtype)
        set_as_sort1(self)

    def add_sort1(self, dt, node_id, grid_type, v1, v2, v3, v4, v5, v6):
        """unvectorized method for adding SORT1 transient data"""
        assert self.is_sort1, self
        assert isinstance(node_id, int) and node_id > 0, 'dt=%s node_id=%s' % (dt, node_id)
        if grid_type in NULL_GRIDTYPE:
            grid_type = -1
        # itotal - the node number
        # itime - the time/frequency step

        # the times/freqs
        self._times[self.itime] = dt
        self.node_gridtype[self.itotal, :] = [node_id, grid_type]
        self.data[self.itime, self.itotal, :] = [v1, v2, v3, v4, v5, v6]
        self.itotal += 1

    def add_sort2(self, dt, node_id, grid_type, v1, v2, v3, v4, v5, v6):
        assert self.is_sort2, self.sort_method
        if grid_type in NULL_GRIDTYPE:
            grid_type = -1
        #self
        #if node_id < 1:
            #msg = self.code_information()
            #msg += "(%s, %s) dt=%g node_id=%s v1=%g v2=%g v3=%g" % (
                #self.itotal, self.itime, dt, node_id, v1, v2, v3)
            ##msg += "                    v4=%g v5=%g v6=%g" % (v4, v5, v6)
            #raise RuntimeError(msg)
        #print(msg)
        self._times[self.itotal] = dt

        # itotal - the time/frequency step
        # itime - the node number
        #print('itime=%s' % self.itime)
        itime = self.itotal
        inode = self.itime
        #print(f'dt={dt} nid={node_id} '
              #f'itime={itime}/{self.ntotal}={self.data.shape[0]} '
              #f'inode={self.itime}/{self.ntimes}={self.data.shape[1]}')
        self.node_gridtype[inode, :] = [node_id, grid_type]
        self.data[itime, inode, :] = [v1, v2, v3, v4, v5, v6]
        #self.ntimes /=
        self.itotal += 1
        #self.itime += 1

    def _write_table_3(self, op2_file, fascii, new_result, itable=-3, itime=0):
        import inspect
        frame = inspect.currentframe()
        call_frame = inspect.getouterframes(frame, 2)
        fascii.write('%s.write_table_3: %s\n' % (self.__class__.__name__, call_frame[1][3]))

        if new_result and itable != -3:
            header = [
                4, 146, 4,
            ]
        else:
            header = [
                4, itable, 4,
                4, 1, 4,
                4, 0, 4,
                4, 146, 4,
            ]
        op2_file.write(pack(b'%ii' % len(header), *header))
        fascii.write('table_3_header = %s\n' % header)
        #op2_file.write(pack('12i', *[4, itable, 4,
                              #4, 1, 4,
                              #4, 0, 4,
                              #4, 146, 4,
                              #]))

        approach_code = self.approach_code
        table_code = self.table_code
        isubcase = self.isubcase
        random_code = self.random_code
        format_code = 1
        num_wide = self.num_wide
        acoustic_flag = self.acoustic_flag if hasattr(self, 'acoustic_flag') else 0
        thermal = self.thermal
        title = b'%-128s' % self.title.encode('ascii')
        subtitle = b'%-128s' % self.subtitle.encode('ascii')  # missing superelement_adaptivity_index
        label = b'%-128s' % self.label.encode('ascii')
        oCode = 0

        ftable3 = b'i' * 50 + b'128s 128s 128s'
        field6 = 0
        field7 = 0

        if isinstance(acoustic_flag, float_types):
            ftable3 = set_table3_field(ftable3, 12, b'f') # field 11

        #print(self.get_stats())
        if self.analysis_code == 1:
            #if hasattr(self, 'lsdvmns'):
            field5 = self.lsdvmns[itime]
            #else:
                #field5 = self.dts[itime]
                #assert isinstance(field5, float_types), type(field5)
                #ftable3 = set_table3_field(ftable3, 5, b'f') # field 5

        elif self.analysis_code == 2:
            field5 = self.modes[itime]
            field6 = self.eigns[itime]
            field7 = self.mode_cycles[itime]
            assert isinstance(field6, float_types), f'field6={field6} type={type(field6)}'
            assert isinstance(field7, float_types), f'field5={field5} field6={field6} field7={field7} type={type(field7)}'
            ftable3 = set_table3_field(ftable3, 6, b'f') # field 6
            ftable3 = set_table3_field(ftable3, 7, b'f') # field 7
        elif self.analysis_code == 5:
            field5 = self.freqs[itime]
            assert isinstance(field5, float_types), f'field5={field5} type={type(field5)}'
            ftable3 = set_table3_field(ftable3, 5, b'f') # field 5
        elif self.analysis_code == 6:
            if hasattr(self, 'dts'):
                field5 = self.dts[itime]
                #assert isinstance(field5, float), type(field5)
            else:
                field5 = self.times[itime]
                #assert isinstance(field5, float), type(field5)
            ftable3 = set_table3_field(ftable3, 5, b'f') # field 5
        elif self.analysis_code == 7:  # pre-buckling
            field5 = self.lsdvmns[itime] # load set number
        elif self.analysis_code == 8:  # post-buckling
            field5 = self.lsdvmns[itime] # load set number
            field5 = int(field5)
            if hasattr(self, 'eigns'):
                field6 = self.eigns[itime]
            elif hasattr(self, 'eigrs'):
                field6 = self.eigrs[itime]
            else:  # pragma: no cover
                raise NotImplementedError('cant find eigns or eigrs on analysis_code=8')
            assert isinstance(field6, float_types), f'field6={field6} type={type(field6)}'
            ftable3 = set_table3_field(ftable3, 6, b'f') # field 6
        elif self.analysis_code == 9:  # complex eigenvalues
            field5 = self.modes[itime]
            if hasattr(self, 'eigns'):
                field6 = self.eigns[itime]
                ftable3 = set_table3_field(ftable3, 6, b'f') # field 6
            #elif hasattr(self, 'eigrs'):
                #field6 = self.eigrs[itime]
                #ftable3 = set_table3_field(ftable3, 6, b'f') # field 6
            #else:  # pragma: no cover
                #raise NotImplementedError('cant find eigns or eigrs on analysis_code=9')
            field7 = self.eigis[itime]
            ftable3 = set_table3_field(ftable3, 7, b'f') # field 7
        elif self.analysis_code == 10:  # nonlinear statics
            field5 = self.lftsfqs[itime]
            ftable3 = set_table3_field(ftable3, 5, b'f') # field 5; load step
        elif self.analysis_code == 11:  # old geometric nonlinear statics
            field5 = self.lsdvmns[itime] # load set number
        else:
            raise NotImplementedError(self.analysis_code)

        table3 = [
            approach_code, table_code, 0, isubcase, field5,
            field6, field7, random_code, format_code, num_wide,
            oCode, acoustic_flag, 0, 0, 0,
            0, 0, 0, 0, 0,
            0, 0, thermal, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            0, 0, 0, 0,
            title, subtitle, label,
        ]
        table3_names = [
            'approach_code', 'table_code', '0', 'isubcase', 'field5',
            'field6', 'field7', 'random_code', 'format_code', 'num_wide',
            'oCode', 'acoustic_flag', 0, 0, 0,
            0, 0, 0, 0, 0,
            0, 0, 'thermal', 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            0, 0, 0, 0,
            title, subtitle, label,
        ]
        assert table3[22] == thermal

        n = 0
        for i, val, ftable3i in zip(count(), table3, ftable3.decode('ascii')):
            assert val is not None, 'i=%s val=%s ftable3i=%s\n%s' % (i, val, ftable3i, self.get_stats())
            if isinstance(val, integer_types):
                n += 4
                assert ftable3i == 'i', 'analysis_code=%r i=%r val=%r type=%r -> name=%r' % (self.analysis_code, i, val, ftable3i, table3_names[i])
            elif isinstance(val, float_types):
                n += 4
                assert ftable3i == 'f', 'analysis_code=%r i=%r val=%r type=%r' % (self.analysis_code, i, val, ftable3i)
            else:
                n += len(val)
        assert n == 584, n
        data = [584] + table3 + [584]
        fmt = b'i' + ftable3 + b'i'

        #op2_file.write(pack(fascii, '%s header 3c' % self.table_name, fmt, data))
        fascii.write('%s header 3c = %s\n' % (self.table_name, data))

        #j = 7
        #print(ftable3[:j])
        #print(table3[:j])
        #pack(ftable3[:j], *table3[:j])
        op2_file.write(pack(fmt, *data))


class RealTableArray(TableArray):
    """
    displacement style table
    """
    def __init__(self, data_code, is_sort1, isubcase, dt):
        TableArray.__init__(self, data_code, is_sort1, isubcase, dt)

    def set_as_static_case(self):
        analysis_code = 1 # static
        device_code = 2  # Plot
        approach_code = analysis_code * 10 + device_code

        self.table_code = table_name_to_table_code[self.table_name_str]
        self.nonlinear_factor = None
        self.data_code['lsdvmns'] = [0] # TODO: ???
        self.data_code['data_names'] = []
        self.data_code['analysis_code'] = analysis_code
        self.data_code['approach_code'] = approach_code
        self.analysis_code = analysis_code
        self.approach_code = approach_code
        self.data_names = []
        self.lsdvmns = [0]
        self._times = [None]

    @classmethod
    def add_static_case(cls, table_name, node_gridtype, data, isubcase,
                        is_sort1=True, is_random=False, is_msc=True,
                        random_code=0, title='', subtitle='', label=''):
        data_code = oug_data_code(table_name,
                                  is_real=True,
                                  is_sort1=is_sort1, is_random=is_random,
                                  random_code=random_code,
                                  title=title, subtitle=subtitle, label=label,
                                  is_msc=is_msc)

        obj = set_static_case(cls, is_sort1, isubcase, data_code,
                              set_real_table, (node_gridtype, data))

        #data_code['lsdvmns'] = [0] # TODO: ???
        #data_code['data_names'] = []

        #times = [None]
        #obj = _set_class(cls, data_code, is_sort1, isubcase, times,
                         #node_gridtype, data)
        return obj

    @classmethod
    def add_transient_case(cls, table_name, node_gridtype, data, isubcase,
                           times,
                           is_sort1=True, is_random=False, is_msc=True,
                           random_code=0, title='', subtitle='', label=''):
        data_code = oug_data_code(table_name,
                                  is_real=True,
                                  is_sort1=is_sort1, is_random=is_random,
                                  random_code=random_code, title=title, subtitle=subtitle, label=label,
                                  is_msc=is_msc)
        obj = set_transient_case(cls, is_sort1, isubcase, data_code,
                                 set_real_table, (node_gridtype, data), times)
        return obj

    @classmethod
    def add_modal_case(cls, table_name, node_gridtype, data, isubcase,
                       modes, eigenvalues, mode_cycles,
                       is_sort1=True, is_random=False, is_msc=True,
                       random_code=0, title='', subtitle='', label=''):

        #elif self.analysis_code == 2:  # real eigenvalues
            ## mode number
            #self.mode = self.add_data_parameter(data, 'mode', b'i', 5)
            ## eigenvalue
            #self.eign = self.add_data_parameter(data, 'eign', b'f', 6, False)
            ## mode or cycle .. todo:: confused on the type - F1???
            #self.mode_cycle = self.add_data_parameter(data, 'mode_cycle', b'i', 7, False)
            #self.update_mode_cycle('mode_cycle')
            #self.data_names = self.apply_data_code_value('data_names', ['mode', 'eign', 'mode_cycle'])

        data_code = oug_data_code(table_name,
                                  is_real=True,
                                  is_sort1=is_sort1, is_random=is_random,
                                  random_code=random_code, title=title, subtitle=subtitle, label=label,
                                  is_msc=is_msc)
        obj = set_modal_case(cls, is_sort1, isubcase, data_code,
                             set_real_table, (node_gridtype, data),
                             modes, eigenvalues, mode_cycles)
        #data_code['modes'] = modes
        #data_code['eigns'] = eigenvalues
        #data_code['mode_cycles'] = mode_cycles
        #data_code['data_names'] = ['modes', 'eigns', 'mode_cycles']

        #obj = set_table_class(cls, data_code, is_sort1, isubcase, modes,
                              #node_gridtype, data)
        #obj.modes = modes
        #obj.eigns = eigenvalues
        #obj.mode_cycles = mode_cycles
        return obj


    def __pos__(self) -> RealTableArray:
        """positive; +a"""
        return self

    def __neg__(self) -> RealTableArray:
        """negative; -a"""
        new_table = copy.deepcopy(self)
        new_table.data *= -1.0
        return new_table

    def __add__(self, table: RealTableArray) -> RealTableArray:
        """a + b"""
        if isinstance(table, RealTableArray):
            self._check_math(table)
            new_data = self.data + table.data
        elif isinstance(table, (integer_types, float_types)):
            new_data = self.data + table
        else:
            raise TypeError(table)
        new_table = copy.deepcopy(self)
        new_table.data = new_data
        return new_table

    # __radd__: reverse order adding (b+a)
    def __iadd__(self, table: RealTableArray) -> RealTableArray:
        """inplace adding; a += b"""
        self._check_math(table)
        self.data += table.data
        return self

    def __sub__(self, table: RealTableArray) -> RealTableArray:
        """a - b"""
        if isinstance(table, RealTableArray):
            self._check_math(table)
            new_data = self.data - table.data
        elif isinstance(table, (integer_types, float_types)):
            new_data = self.data - table
        else:
            raise TypeError(table)
        new_table = copy.deepcopy(self)
        new_table.data = new_data
        return new_table

    def __mul__(self, table: RealTableArray) -> RealTableArray:
        """a * b"""
        if isinstance(table, RealTableArray):
            self._check_math(table)
            new_data = self.data * table.data
        elif isinstance(table, (integer_types, float_types)):
            new_data = self.data * table
        else:
            raise TypeError(table)
        new_table = copy.deepcopy(self)
        new_table.data = new_data
        return new_table

    def __imul__(self, table: RealTableArray) -> RealTableArray:
        """a *= b"""
        if isinstance(table, RealTableArray):
            self._check_math(table)
            self.data *= table.data
        elif isinstance(table, (integer_types, float_types)):
            self.data *= table
        else:
            raise TypeError(table)
        return self

    def __truediv__(self, table: RealTableArray) -> RealTableArray:
        """a / b"""
        if isinstance(table, RealTableArray):
            self._check_math(table)
            new_data = self.data / table.data
        elif isinstance(table, (integer_types, float_types)):
            new_data = self.data / table
        else:
            raise TypeError(table)
        new_table = copy.deepcopy(self)
        new_table.data = new_data
        return new_table

    def linear_combination(self, factor: integer_float_types,
                           data: Optional[np.ndarray]=None,
                           update: bool=True) -> None:
        combination_inplace(self.data, data, factor)
        # if update:
        #     self.update_data_components()

    def update_data_components(self):
        return

    def _check_math(self, table: RealTableArray) -> None:
        """verifies that the shapes are the same"""
        assert self.ntimes == table.ntimes, f'ntimes={self.ntimes} table.times={table.ntimes}'
        assert self.ntotal == table.ntotal, f'ntotal={self.ntotal} table.ntotal={table.ntotal}'
        assert self.node_gridtype.shape == table.node_gridtype.shape, f'node_gridtype.shape={self.node_gridtype.shape} table.node_gridtype.shape={table.node_gridtype.shape}'
        assert self.data.shape == table.data.shape, f'data.shape={self.data.shape} table.data.shape={table.data.shape}'

    @property
    def is_real(self) -> bool:
        return True

    @property
    def is_complex(self) -> bool:
        return False

    def data_type(self) -> str:
        return 'float32'

    def write_op2(self, op2_file, fascii, itable: int, new_result,
                  date, is_mag_phase: bool=False, endian: str='>'):
        """writes an OP2"""
        import inspect
        allowed_tables = [
            'OUGV1', 'BOUGV1',
            'OPHIG', 'BOPHIG',
            'OUPV1', 'OUXY1', # solution set
            'OQP1', 'OQMG1', 'OQG1', 'OQGV1', 'OPNL1',
            'OPG1', 'BOPG1',
            'OPHSA',
            'OPGV1',
            'OUGATO1', 'OUGCRM1', 'OUGNO1', 'OUGPSD1', 'OUGRMS1', # disp/vel/acc/eigenvector
            'OVGATO1', 'OVGCRM1', 'OVGNO1',
            'OAGATO1', 'OAGCRM1', 'OAGNO1', 'OAGPSD1', 'OAGRMS1', # acceleration
                                  'OPGNO1',            'OPGRMS1', # load vector
            'OPRATO1', 'OPRCRM1', 'OPRNO1', 'OPRPSD1',            # pressure
                                  'OQGNO1', 'OQGPSD1', 'OQGRMS1',
                       'OQMCRM1',           'OQMPSD1', # 'OQMRMS1',#'OQMNO1',
            'OCRPG', 'OCRUG', 'OUG1', 'OVG1', 'OAG1',
            'OUGV1PAT',
            'OUGF1',
            'OQGCF1', 'OQGGF1',
            'RADCONS', 'RADEATC', 'RADEFFM',
        ]
        assert self.table_name in allowed_tables, self.table_name

        frame = inspect.currentframe()
        call_frame = inspect.getouterframes(frame, 2)
        fascii.write('%s.write_op2: %s\n' % (self.__class__.__name__, call_frame[1][3]))

        if itable == -1:
            self._write_table_header(op2_file, fascii, date)
            itable = -3

        #print('nonlinear_factor =', self.nonlinear_factor)
        if not self.is_sort1:
            raise NotImplementedError('SORT2')
        #op2_format = endian + b'2i6f'

        node = self.node_gridtype[:, 0]
        gridtype = self.node_gridtype[:, 1]

        # table 4 info
        #ntimes = self.data.shape[0]
        nnodes = self.data.shape[1]
        nnodes_device = self.node_gridtype[:, 0] * 10 + self.device_code

        nnodes = len(node)
        max_id = node.max()
        if max_id > MAX_32_BIT_INT:
            raise SixtyFourBitError(f'64-bit OP2 writing is not supported; max id={max_id}')

        fdtype = self.data.dtype
        if self.size == fdtype.itemsize:
            pass
        else:
            # warnings.warn(f'downcasting {self.class_name}...this is buggy')
            #idtype = np.int32(1)
            fdtype = np.float32(1.0)

        nodedevice_gridtype = np.column_stack([nnodes_device, gridtype])
        node_gridtype_floats = view_idtype_as_fdtype(nodedevice_gridtype, fdtype)

        #print(node_gridtype_floats)
        #node_gridtype_floats = nodedevice_gridtype.view(fdtype) # .reshape(nnodes, 2)

        #(2+6) => (node_id, gridtypei, t1i, t2i, t3i, r1i, r2i, r3i)
        ntotal = nnodes * (2 + 6)

        #print('shape = %s' % str(self.data.shape))
        #assert nnodes > 1, nnodes
        assert ntotal > 1, ntotal

        unused_device_code = self.device_code
        fascii.write('  ntimes = %s\n' % self.ntimes)

        #fmt = '%2i %6f'
        #print('ntotal=%s' % (ntotal))
        for itime in range(self.ntimes):
            self._write_table_3(op2_file, fascii, new_result, itable, itime)

            # record 4
            itable -= 1
            header = [4, itable, 4,
                      4, 1, 4,
                      4, 0, 4,
                      4, ntotal, 4,
                      4*ntotal]
            op2_file.write(pack(b'%ii' % len(header), *header))
            fascii.write('r4 [4, 0, 4]\n')
            fascii.write(f'r4 [4, {itable}, 4]\n')
            fascii.write(f'r4 [4, {4*ntotal:d}, 4]\n')

            datai = view_dtype(self.data[itime, :, :], fdtype)
            node_gridtype_data = np.hstack([node_gridtype_floats, datai])
            op2_file.write(node_gridtype_data)
            assert ntotal == node_gridtype_data.size
            #print(self.data.shape, ntotal)

            itable -= 1
            header = [4 * ntotal,]
            op2_file.write(pack(b'i', *header))
            fascii.write('footer = %s\n' % header)
            new_result = False
        return itable

    def write_csv(self, csv_file: TextIO,
                  is_exponent_format: bool=False,
                  is_mag_phase: bool=False, is_sort1: bool=True,
                  write_header: bool=True):
        """
        Displacement Table
        ------------------
        Flag, SubcaseID,  iTime, NID,       dx,      dy,       dz,      rx,       ry,      rz,  cd,  PointType
        1,            1,  0,     101, 0.014159, 0.03448, 0.019135, 0.00637, 0.008042, 0.00762,   0,  1
        uses cd=-1 for unknown cd

        """
        name = str(self.__class__.__name__)
        if write_header:
            csv_file.write('%s\n' % name)
            headers = ['Flag', 'Subcase', 'iTime', 'Node', ] + self.headers + ['cd', 'PointType']
            csv_file.write('# ' + ','.join(headers) + '\n')
        node = self.node_gridtype[:, 0]
        gridtype = self.node_gridtype[:, 1]

        #unused_times = self._times
        isubcase = self.isubcase

        flag_map = {
            'RealDisplacementArray': 1,
            'RealVelocityArray': 2,
            'RealAccelerationArray': 3,
            'RealEigenvectorArray': 4,

            'RealLoadVectorArray': 5,
            #'RealAppliedLoadsArray': 5,
            'RealSPCForcesArray': 6,
            'RealMPCForcesArray': 7,
            #'Temperature' : 8,

            #'HeatFlux' : 9,
            'RealTemperatureGradientAndFluxArray': 9,
        }
        flag = flag_map[name]

        # sort1 as sort1
        assert is_sort1 is True, is_sort1
        nid_len = '%d' % len(str(node.max()))
        cd = -1
        for itime in range(self.ntimes):
            #dt = self._times[itime]
            t1 = self.data[itime, :, 0]
            t2 = self.data[itime, :, 1]
            t3 = self.data[itime, :, 2]
            r1 = self.data[itime, :, 3]
            r2 = self.data[itime, :, 4]
            r3 = self.data[itime, :, 5]
            for node_id, gridtypei, t1i, t2i, t3i, r1i, r2i, r3i in zip(node, gridtype, t1, t2, t3, r1, r2, r3):
                #unused_sgridtype = self.recast_gridtype_as_string(gridtypei)
                assert is_exponent_format
                if is_exponent_format:
                    vals2 = write_floats_13e_long([t1i, t2i, t3i, r1i, r2i, r3i])
                    (t1i, t2i, t3i, r1i, r2i, r3i) = vals2

                csv_file.write(f'{flag}, {isubcase}, {itime}, {node_id:{nid_len}d}, '
                               f'{t1i}, {t2i}, {t3i}, {r1i}, {r2i}, {r3i}, {cd}, {gridtypei}\n')
        return

    def write_frd(self, frd_file: TextIO,
                  is_exponent_format: bool=False,
                  is_mag_phase: bool=False, is_sort1: bool=True,
                  write_header: bool=True):
        """
        https://web.mit.edu/calculix_v2.7/CalculiX/cgx_2.7/doc/cgx/node174.html

        Nodal Results Block
        Purpose: Stores values on node positions

        1. Record:
        Format:(1X,' 100','C',6A1,E12.5,I12,20A1,I2,I5,10A1,I2)
        Values: KEY,CODE,SETNAME,VALUE,NUMNOD,TEXT,ICTYPE,NUMSTP,ANALYS,
                FORMAT
        Where: KEY    = 100
               CODE   = C
               SETNAME= Name (not used)
               VALUE  = Could be frequency, time or any numerical value
               NUMNOD = Number of nodes in this nodal results block
               TEXT   = Any text
               ICTYPE = Analysis type
                        0  static
                        1  time step
                        2  frequency
                        3  load step
                        4  user named
               NUMSTP = Step number
               ANALYS = Type of analysis (description)
               FORMAT = Format indicator
                        0  short format
                        1  long format
                        2  binary format

        #--------------------------------------------------------------------------
        2. Record:
        Format:(1X,I2,2X,8A1,2I5)
        Values: KEY, NAME, NCOMPS, IRTYPE
        Where: KEY    = -4
               NAME   = Dataset name to be used in the menu
               NCOMPS = Number of entities
               IRTYPE = 1  Nodal data, material independent
                        2  Nodal data, material dependant
                        3  Element data at nodes (not used)

        #--------------------------------------------------------------------------
        3. Type of Record:
        Format:(1X,I2,2X,8A1,5I5,8A1)
        Values: KEY, NAME, MENU, ICTYPE, ICIND1, ICIND2, IEXIST, ICNAME
        Where: KEY    = -5
               NAME   = Entity name to be used in the menu for this comp.
               MENU   = 1
               ICTYPE = Type of entity
                        1  scalar
                        2  vector with 3 components
                        4  matrix
                       12  vector with 3 amplitudes and 3 phase-angles in
                           degree
               ICIND1 = sub-component index or row number
               ICIND2 = column number for ICTYPE=4
               IEXIST = 0  data are provided
                        1  data are to be calculated by predefined
                           functions (not used)
                        2  as 0 but earmarked
               ICNAME = Name of the predefined calculation (not used)
                        ALL  calculate the total displacement if ICTYPE=2
        This record must be repeated for each entity.

        4. Type of Record:  (not used)
        This record will be necessary in combination with the request for
        predefined calculations. This type of record is not allowed in
        combination with binary coding of data.
        Format:(1X,I2,2I5,20I3)
        Values: KEY,IRECTY,NUMCPS,(LSTCPS(I),I=1,NUMCPS)
        Where: KEY    = -6
               IRECTY = Record variant identification number
               NUMCPS = Number of components
               LSTCPS = For each variant component, the position of the
                        corresponding component in attribute definition
        #--------------------------------------------------------------------------

        5. Type of Record:
        The following records are data-records and the format is repeated
        for each node.

        In case of material independent data

        - ascii coding:
        Following records (ascci, FORMAT=0 | 1):
         Short Format:(1X,I2,I5,6E12.5)
         Long Format:(1X,I2,I10,6E12.5)
         Values: KEY, NODE, XX..
         Where: KEY  = -1 if its the first line of data for a given node
                       -2 if its a continuation line
               NODE  = node number or blank if KEY=-2
               XX..  = data

        - binary coding:
         Following records (ascci, FORMAT=2):
         (int,NCOMPS*float)
         int and float are ansi-c data-types
         Values: NODE, XX..
         Where:
               NODE   = node number or blank if KEY=-2
               XX..   = data

        In case of material dependant data
        REMARK: Implemented only for NMATS=1
        - first line:
        Short Format:(1X,I2,4I5)
        Long Format:(1X,I2,I10,3I5)
        Values: KEY, NODENR, NMATS
        Where: KEY    = -1
               NODENR = Node number
               NMATS  = Number of different materials at this node(unused)
        - second and following lines:
        Short Format:(1X,I2,I5,6E12.5)
        Long Format:(1X,I2,I10,6E12.5)
        Values: KEY, MAT, XX, YY, ZZ, XY, YZ, ZX ..
        Where: KEY    = -2
               MAT    = material-property-number if KEY=-2 (unused)
               XX..   = data


        Last Record (only FORMAT=0 | 1 (ascii), omitted for FORMAT=2):
        Format:(1X,'-3')
        Values: KEY
        Displacement Table
        ------------------
        Flag, SubcaseID,  iTime, NID,       dx,      dy,       dz,      rx,       ry,      rz,  cd,  PointType
        1,            1,  0,     101, 0.014159, 0.03448, 0.019135, 0.00637, 0.008042, 0.00762,   0,  1
        uses cd=-1 for unknown cd

        """
        name = str(self.__class__.__name__)
        name_map = {
            'RealDisplacementArray': 'Displacement',
            'RealSPCForcesArray': 'SPC_Force',
            'RealMPCForcesArray': 'MPC_Force',
        }
        name = name_map[name]
        #if write_header:
            #csv_file.write('%s\n' % name)
            #headers = ['Flag', 'Subcase', 'iTime', 'Node', ] + self.headers + ['cd', 'PointType']
            #csv_file.write('# ' + ','.join(headers) + '\n')

        node = self.node_gridtype[:, 0]
        #gridtype = self.node_gridtype[:, 1]

        #unused_times = self._times
        #isubcase = self.isubcase

        # sort1 as sort1
        assert is_sort1 is True, is_sort1
        #nid_len = '%d' % len(str(node.max()))
        #cd = -1
        num_step, nnode = self.data.shape[:2]
        for itime in range(self.ntimes):
            #dt = self._times[itime]
            #1. Record:
            #Format:(1X,  ' 100',    'C',      6A1, E12.5,    I12, 20A1,I2,I5,10A1,I2)
            #Values:       KEY,     CODE,  SETNAME, VALUE, NUMNOD, TEXT,ICTYPE,NUMSTP,ANALYS,
            #        FORMAT
            #Where: KEY    = 100
            #       CODE   = C
            #       SETNAME= Name (not used)
            #       VALUE  = Could be frequency, time or any numerical value
            #       NUMNOD = Number of nodes in this nodal results block
            #       TEXT   = Any text
            #       ICTYPE = Analysis type
            #                0  static
            #                1  time step
            #                2  frequency
            #                3  load step
            #                4  user named
            #       NUMSTP = Step number
            #       ANALYS = Type of analysis (description)
            #       FORMAT = Format indicator
            #                0  short format
            #                1  long format
            #                2  binary format
            # basically table 3
            key = 100

            text = name
            text_str = f'{text:<20s}'
            assert len(text_str) == 20, len(text_str)
            value = self._times[0]
            map_analysis = {
                # displacement -> ???
                1: 1,
            }
            analysis = map_analysis[self.analysis_code]

            code = 'C'
            ic_type = '0'
            data_format = 2
            set_name = 'SET_NAME'
            num_mod = nnode
            frd_file.write(f'RECORD 1: {key} {code} {set_name} {value} {num_mod} {text} {ic_type} {num_step} {analysis} {data_format}\n')

            # 2. Record:
            # Format:(1X, I2, 2X,8A1,2I5)
            # Values:    KEY, NAME, NCOMPS, IRTYPE
            # Where: KEY    = -4
            #        NAME   = Dataset name to be used in the menu
            #        NCOMPS = Number of entities
            #        IRTYPE = 1  Nodal data, material independent
            #                 2  Nodal data, material dependant
            #                 3  Element data at nodes (not used)
            key = -4
            name = f'disp{itime:d}'
            ncomps = 6
            ir_type = 1
            frd_file.write(f'RECORD 2: {key} {name} {ncomps} {ir_type}\n')

            # 3. Type of Record:
            # Format:(1X,I2,2X,8A1,5I5,8A1)
            # Values: KEY, NAME, MENU, ICTYPE, ICIND1, ICIND2, IEXIST, ICNAME
            # Where: KEY    = -5
            #        NAME   = Entity name to be used in the menu for this comp.
            #        MENU   = 1
            #        ICTYPE = Type of entity
            #                 1  scalar
            #                 2  vector with 3 components
            #                 4  matrix
            #                12  vector with 3 amplitudes and 3 phase-angles in
            #                    degree
            #        ICIND1 = sub-component index or row number
            #        ICIND2 = column number for ICTYPE=4
            #        IEXIST = 0  data are provided
            #                 1  data are to be calculated by predefined
            #                    functions (not used)
            #                 2  as 0 but earmarked
            #        ICNAME = Name of the predefined calculation (not used)
            #                 ALL  calculate the total displacement if ICTYPE=2
            # This record must be repeated for each entity
            key = -5
            name = 'Translation'
            menu = 1
            ic_type = 2
            iexist = 0
            ic_name = 'Translation'
            icind1 = itime
            icind2 = 0
            frd_file.write(f'RECORD 3: {key} {name} {menu} {ic_type} {icind1} {icind2} {iexist} {ic_name}\n')

            # 5. Type of Record:
            # The following records are data-records and the format is repeated
            # for each node.
            #
            # In case of material independent data
            #
            # - ascii coding:
            # Following records (ascii, FORMAT=0 | 1):
            #  Short Format:(1X,I2,I5,6E12.5)
            #  Long Format:(1X,I2,I10,6E12.5)
            #  Values: KEY, NODE, XX..
            #  Where: KEY  = -1 if its the first line of data for a given node
            #                -2 if its a continuation line
            #        NODE  = node number or blank if KEY=-2
            #        XX..  = data
            #
            # - binary coding:
            #  Following records (ascci, FORMAT=2):
            #  (int,NCOMPS*float)
            #  int and float are ansi-c data-types
            #  Values: NODE, XX..
            #  Where:
            #        NODE   = node number or blank if KEY=-2
            #        XX..   = data

            t1 = self.data[itime, :, 0]
            t2 = self.data[itime, :, 1]
            t3 = self.data[itime, :, 2]
            r1 = self.data[itime, :, 3]
            r2 = self.data[itime, :, 4]
            r3 = self.data[itime, :, 5]
            for node_id, t1i, t2i, t3i in zip(node, t1, t2, t3):
                frd_file.write(f'RECORD 5: 1 -1 {node_id} {t1i:12.5E} {t2i:12.5E} {t3i:12.5E}\n')

            # This record must be repeated for each entity
            #key = -5
            name = 'Rotation'
            #menu = 1
            #ic_type = 2
            #iexist = 0
            ic_name = 'Rotation'
            frd_file.write(f'RECORD 3:  {key} {name} {menu} {ic_type} {icind1} {icind2} {iexist} {ic_name}\n')
            for node_id, r1i, r2i, r3i in zip(node, r1, r2, r3):
                frd_file.write(f'RECORD 5: 1 -1 {node_id} {r1i:12.5E} {r2i:12.5E} {r3i:12.5E}\n')
        return

    def _write_f06_block(self, words, header, page_stamp, page_num, f06_file: TextIO,
                         write_words,
                         is_mag_phase: bool=False, is_sort1: bool=True):
        if write_words:
            words += [' \n', '      POINT ID.   TYPE          T1             T2             T3             R1             R2             R3\n']
        #words += self.getTableMarker()
        f06_file.write(''.join(header + words))

        node = self.node_gridtype[:, 0]
        gridtype = self.node_gridtype[:, 1]
        t1 = self.data[0, :, 0]
        t2 = self.data[0, :, 1]
        t3 = self.data[0, :, 2]
        r1 = self.data[0, :, 3]
        r2 = self.data[0, :, 4]
        r3 = self.data[0, :, 5]
        for node_id, gridtypei, t1i, t2i, t3i, r1i, r2i, r3i in zip(node, gridtype, t1, t2, t3, r1, r2, r3):
            sgridtype = recast_gridtype_as_string(self, gridtypei)
            vals = [t1i, t2i, t3i, r1i, r2i, r3i]
            vals2 = write_floats_13e(vals)
            (dx, dy, dz, rx, ry, rz) = vals2
            f06_file.write('%14i %6s     %-13s  %-13s  %-13s  %-13s  %-13s  %s\n' % (
                node_id, sgridtype, dx, dy, dz, rx, ry, rz))
        f06_file.write(page_stamp % page_num)
        return page_num

    def _write_sort1_as_sort2(self, f06_file: TextIO, page_num, page_stamp, header, words):
        nodes = self.node_gridtype[:, 0]
        gridtypes = self.node_gridtype[:, 1]
        times = self._times

        for inode, (node_id, gridtypei) in enumerate(zip(nodes, gridtypes)):
            t1 = self.data[:, inode, 0].ravel()
            t2 = self.data[:, inode, 1].ravel()
            t3 = self.data[:, inode, 2].ravel()
            r1 = self.data[:, inode, 3].ravel()
            r2 = self.data[:, inode, 4].ravel()
            r3 = self.data[:, inode, 5].ravel()

            header[1] = ' POINT-ID = %10i\n' % node_id
            f06_file.write(''.join(header + words))
            for dt, t1i, t2i, t3i, r1i, r2i, r3i in zip(times, t1, t2, t3, r1, r2, r3):
                sgridtype = recast_gridtype_as_string(self, gridtypei)
                vals = [t1i, t2i, t3i, r1i, r2i, r3i]
                vals2 = write_floats_13e(vals)
                (dx, dy, dz, rx, ry, rz) = vals2
                if sgridtype in ['G', 'H', 'L']:
                    f06_file.write('%14s %6s     %-13s  %-13s  %-13s  %-13s  %-13s  %s\n' % (
                        write_float_12e(dt), sgridtype, dx, dy, dz, rx, ry, rz))
                elif sgridtype in ['S', 'M', 'E']:
                    f06_file.write('%14s %6s     %s\n' % (node_id, sgridtype, dx))
                else:  # pragma: no cover
                    raise NotImplementedError(sgridtype)
            f06_file.write(page_stamp % page_num)
            page_num += 1
        return page_num

    def _write_sort1_as_sort1(self, f06_file: TextIO, page_num, page_stamp, header, words):
        nodes = self.node_gridtype[:, 0]
        gridtypes = self.node_gridtype[:, 1]
        unused_times = self._times

        for itime in range(self.ntimes):
            dt = self._times[itime]
            t1 = self.data[itime, :, 0]
            t2 = self.data[itime, :, 1]
            t3 = self.data[itime, :, 2]
            r1 = self.data[itime, :, 3]
            r2 = self.data[itime, :, 4]
            r3 = self.data[itime, :, 5]

            if isinstance(dt, float_types):
                header[1] = ' %s = %10.4E\n' % (self.data_code['name'], dt)
            else:
                header[1] = ' %s = %10i\n' % (self.data_code['name'], dt)
            f06_file.write(''.join(header + words))
            for node_id, gridtypei, t1i, t2i, t3i, r1i, r2i, r3i in zip(nodes, gridtypes, t1, t2, t3, r1, r2, r3):
                sgridtype = recast_gridtype_as_string(self, gridtypei)
                vals = [t1i, t2i, t3i, r1i, r2i, r3i]
                vals2 = write_floats_13e(vals)
                (dx, dy, dz, rx, ry, rz) = vals2
                if sgridtype in ['G', 'H', 'L']:
                    f06_file.write('%14i %6s     %-13s  %-13s  %-13s  %-13s  %-13s  %s\n' % (
                        node_id, sgridtype, dx, dy, dz, rx, ry, rz))
                elif sgridtype in ['S', 'M', 'E']:
                    f06_file.write('%14i %6s     %s\n' % (node_id, sgridtype, dx))
                else:  # pragma: no cover
                    raise NotImplementedError(f'node_id={node_id} sgridtype={sgridtype} vals={vals2}')
            f06_file.write(page_stamp % page_num)
            page_num += 1
        return page_num

    #def _write_sort2_as_sort2(self, f06_file, page_num, page_stamp, header, words):
        #nodes = self.node_gridtype[:, 0]
        #gridtypes = self.node_gridtype[:, 1]
        #times = self._times
        #for inode, (node_id, gridtypei) in enumerate(zip(nodes, gridtypes)):
            #t1 = self.data[inode, :, 0]
            #t2 = self.data[inode, :, 1]
            #t3 = self.data[inode, :, 2]
            #r1 = self.data[inode, :, 3]
            #r2 = self.data[inode, :, 4]
            #r3 = self.data[inode, :, 5]

            #header[1] = ' POINT-ID = %10i\n' % node_id
            #f06_file.write(''.join(header + words))
            #for dt, t1i, t2i, t3i, r1i, r2i, r3i in zip(times, t1, t2, t3, r1, r2, r3):
                #sgridtype = self.recast_gridtype_as_string(gridtypei)
                #vals = [t1i, t2i, t3i, r1i, r2i, r3i]
                #vals2 = write_floats_13e(vals)
                #(dx, dy, dz, rx, ry, rz) = vals2
                #if sgridtype in ['G', 'H', 'L']:
                    #f06_file.write('%14s %6s     %-13s  %-13s  %-13s  %-13s  %-13s  %s\n' % (
                        #write_float_12e(dt), sgridtype, dx, dy, dz, rx, ry, rz))
                #elif sgridtype == 'S':
                    #f06_file.write('%14s %6s     %s\n' % (node_id, sgridtype, dx))
                #else:
                    #raise NotImplementedError(sgridtype)
            #f06_file.write(page_stamp % page_num)
            #page_num += 1
        #return page_num

    def _write_f06_transient_block(self, words, header, page_stamp, page_num, f06_file: TextIO, write_words,
                                   is_mag_phase=False, is_sort1=True):
        if write_words:
            words += [' \n', '      POINT ID.   TYPE          T1             T2             T3             R1             R2             R3\n']
        #words += self.getTableMarker()

        if not len(header) >= 3:
            header.append('')

        is_sort2 = not is_sort1
        if self.is_sort1 or self.nonlinear_factor in (None, np.nan):
            if is_sort2 and self.nonlinear_factor is not None:
                page_num = self._write_sort1_as_sort2(f06_file, page_num, page_stamp, header, words)
            else:
                page_num = self._write_sort1_as_sort1(f06_file, page_num, page_stamp, header, words)
        else:
            return page_num - 1
            #raise NotImplementedError('SORT2')
            #page_num = self._write_sort2_as_sort2(f06_file, page_num, page_stamp, header, words)
        return page_num - 1

    def extract_xyplot(self, node_ids, index):
        node_ids = np.asarray(node_ids, dtype='int32')
        i = index - 1
        assert index in [1, 2, 3, 4, 5, 6], index
        nids = self.node_gridtype[:, 0]
        inids = np.searchsorted(nids, node_ids)
        assert all(nids[inids] == node_ids), 'nids=%s expected=%s; all=%s'  % (nids[inids], node_ids, nids)
        return self.data[:, inids, i]

def set_real_table(cls, data_code, is_sort1, isubcase,
                   node_gridtype, data, times):
#def set_real_table(cls, data_code, is_sort1, isubcase, times,
                   #node_gridtype, data):
    assert node_gridtype.ndim == 2, node_gridtype.shape
    assert data.ndim == 3, data.shape
    dt = times[0]
    ntimes = data.shape[0]
    nnodes = data.shape[1]
    data_code['sort_code'] = 0
    obj = cls(data_code, is_sort1, isubcase, dt)
    obj.node_gridtype = node_gridtype
    obj.data = data

    obj.ntimes = ntimes
    obj.ntotal = nnodes
    obj._times = times
    return obj


class ComplexTableArray(TableArray):
    """
    complex displacement style table
    """
    def __init__(self, data_code, is_sort1, isubcase, dt):
        TableArray.__init__(self, data_code, is_sort1, isubcase, dt)

    @classmethod
    def add_freq_case(cls, table_name: str,
                      node_gridtype: np.ndarray,
                      data: np.ndarray,
                      isubcase: int,
                      freqs: np.ndarray,
                      is_sort1=True, is_random=False, is_msc=True,
                      random_code=0, title='', subtitle='', label=''):

        table_name = 'OUGV1'
        #analysis_code = 5 # freq
        data_code = oug_data_code(table_name,
                                  is_complex=True,
                                  is_sort1=is_sort1, is_random=is_random,
                                  random_code=random_code, title=title, subtitle=subtitle, label=label,
                                  is_msc=is_msc)
        #data_code['modes'] = modes
        #data_code['eigns'] = eigenvalues
        #data_code['mode_cycles'] = mode_cycles
        #data_code['data_names'] = ['freq']
        #data_code['name'] = 'FREQ'

        #ntimes = data.shape[0]
        #nnodes = data.shape[1]
        #dt = freqs[0]
        #obj = cls(data_code, is_sort1, isubcase, dt)
        #obj.node_gridtype = node_gridtype
        #obj.data = data

        #obj.freqs = freqs

        #obj.ntimes = ntimes
        #obj.ntotal = nnodes
        #obj._times = freqs

        obj = set_freq_case(cls, is_sort1, isubcase, data_code,
                            set_complex_table, (node_gridtype, data),
                            freqs)
        return obj

    @classmethod
    def add_complex_modes_case(cls, table_name, element, data, isubcase,
                               modes, eigrs, eigis,
                               is_sort1=True, is_random=False, is_msc=True,
                               random_code=0, title='', subtitle='', label=''):
        #data_code = cls._add_case(
            #table_name, isubcase,
            #is_sort1, is_random, is_msc,
            #random_code, title, subtitle, label)
        data_code = oug_data_code(table_name,
                                  is_complex=True,
                                  is_sort1=is_sort1, is_random=is_random,
                                  random_code=random_code, title=title, subtitle=subtitle, label=label,
                                  is_msc=is_msc)

        obj = set_complex_modes_case(cls, is_sort1, isubcase, data_code,
                                     set_complex_table, (element, data),
                                     modes, eigrs, eigis)
        return obj

    def extract_xyplot(self, node_ids, index, index_str):
        j = index_str_to_axis(index_str)

        node_ids = np.asarray(node_ids, dtype='int32')
        i = index - 1
        assert index in [1, 2, 3, 4, 5, 6,
                         7, 8, 9, 10, 11, 12], index
        nids = self.node_gridtype[:, 0]
        inids = np.searchsorted(nids, node_ids)
        assert all(nids[inids] == node_ids), 'nids=%s expected=%s; all=%s'  % (nids[inids], node_ids, nids)
        if j == 1:
            # real
            return self.data[:, inids, i].real
        elif j == 2:
            # imag
            return self.data[:, inids, i].imag
        elif j == 3:
            # mag
            return np.abs(self.data[:, inids, i])
        elif j == 4:
            # phase
            return np.angle(self.data[:, inids, i])
        else:
            raise RuntimeError()

    @property
    def is_real(self) -> bool:
        return False

    @property
    def is_complex(self) -> bool:
        return True

    def data_type(self) -> str:
        return 'complex64'

    #def _write_f06_block(self, words, header, page_stamp, page_num, f06_file, is_mag_phase):
        #self._write_f06_transient_block(words, header, page_stamp, page_num, f06_file, is_mag_phase, is_sort1)

    def _write_f06_transient_block(self, words, header, page_stamp, page_num, f06_file,
                                   is_mag_phase, is_sort1):
        if is_mag_phase:
            words += ['                                                         (MAGNITUDE/PHASE)\n', ]
        else:
            words += ['                                                          (REAL/IMAGINARY)\n', ]

        if not len(header) >= 3:
            header.append('')

        #is_sort1_table = self.is_sort1
        is_sort1_table = SORT1_TABLES # self.table_name[-1] == '1'
        if is_sort1_table:
            assert self.table_name in SORT1_TABLES, self.table_name
            if is_sort1:
                words += [' \n', '      POINT ID.   TYPE          T1             T2             T3             R1             R2             R3\n']
                page_num = self.write_sort1_as_sort1(f06_file, page_num, page_stamp, header, words, is_mag_phase)
            else:
                words += [' \n', '      FREQUENCY   TYPE          T1             T2             T3             R1             R2             R3\n']
                page_num = self.write_sort1_as_sort2(f06_file, page_num, page_stamp, header, words, is_mag_phase)
        else:
            assert self.table_name in SORT2_TABLES, self.table_name
            words += [' \n', '      FREQUENCY   TYPE          T1             T2             T3             R1             R2             R3\n']
            page_num = self.write_sort2_as_sort2(f06_file, page_num, page_stamp, header, words, is_mag_phase)
        return page_num - 1

    def write_sort1_as_sort1(self, f06_file, page_num, page_stamp, header, words, is_mag_phase):
        assert self.ntimes == len(self._times), 'ntimes=%s len(self._times)=%s' % (self.ntimes, self._times)
        words_orig = copy.deepcopy(words)

        for itime, dt in enumerate(self._times):
            if hasattr(self, 'eigrs'):
                words = copy.deepcopy(words_orig)
                eigr = self.eigrs[itime]
                eigi = self.eigis[itime]
                eigr = 0. if eigr == 0 else eigr
                eigi = 0. if eigi == 0 else eigi
                if '%' in words[0]:
                    try:
                        words[0] = words[0] % (eigr, eigi)
                    except TypeError:
                        print('words =', words)
                        raise

                if '%' in words[1]:
                    try:
                        words[1] = words[1] % (itime + 1)
                    except TypeError:
                        print('words =', words)
                        raise

            node = self.node_gridtype[:, 0]
            gridtype = self.node_gridtype[:, 1]
            t1 = self.data[itime, :, 0]
            t2 = self.data[itime, :, 1]
            t3 = self.data[itime, :, 2]
            r1 = self.data[itime, :, 3]
            r2 = self.data[itime, :, 4]
            r3 = self.data[itime, :, 5]

            header[2] = ' %s = %10.4E\n' % (self.data_code['name'], dt)
            f06_file.write(''.join(header + words))
            for node_id, gridtypei, t1i, t2i, t3i, r1i, r2i, r3i in zip(node, gridtype, t1, t2, t3, r1, r2, r3):
                sgridtype = recast_gridtype_as_string(self, gridtypei)
                vals = [t1i, t2i, t3i, r1i, r2i, r3i]
                vals2 = write_imag_floats_13e(vals, is_mag_phase)
                [dxr, dyr, dzr, rxr, ryr, rzr,
                 dxi, dyi, dzi, rxi, ryi, rzi] = vals2
                if sgridtype in ['G', 'H']:
                    f06_file.write('0 %12i %6s     %-13s  %-13s  %-13s  %-13s  %-13s  %-s\n'
                                   '  %12s %6s     %-13s  %-13s  %-13s  %-13s  %-13s  %-s\n' % (
                                       node_id, sgridtype, dxr, dyr, dzr, rxr, ryr, rzr,
                                       '', '', dxi, dyi, dzi, rxi, ryi, rzi))
                elif sgridtype in ['S', 'M', 'E']:
                    f06_file.write('0 %12i %6s     %-13s\n'
                                   '  %12s %6s     %-13s\n' % (node_id, sgridtype, dxr, '', '', dxi))
                else:  # pragma: no cover
                    raise NotImplementedError(f'node_id={node_id} sgridtype={sgridtype} vals={vals2}')
            f06_file.write(page_stamp % page_num)
            page_num += 1
        return page_num

    def write_sort1_as_sort2(self, f06_file: TextIO, page_num: int,
                             page_stamp, header, words, is_mag_phase):
        node = self.node_gridtype[:, 0]
        gridtype = self.node_gridtype[:, 1]

        times = self._times
        # print(self.data.shape)
        for inode, (node_id, gridtypei) in enumerate(zip(node, gridtype)):
            # TODO: for SORT1 pretending to be SORT2
            #t1 = self.data[:, inode, 0].ravel()
            t1 = self.data[:, inode, 0].ravel()
            t2 = self.data[:, inode, 1].ravel()
            t3 = self.data[:, inode, 2].ravel()
            r1 = self.data[:, inode, 3].ravel()
            r2 = self.data[:, inode, 4].ravel()
            r3 = self.data[:, inode, 5].ravel()
            if len(r3) != len(times):
                raise RuntimeError('len(d)=%s len(times)=%s' % (len(r3), len(times)))

            header[2] = ' POINT-ID = %10i\n' % node_id
            f06_file.write(''.join(header + words))
            for dt, t1i, t2i, t3i, r1i, r2i, r3i in zip(times, t1, t2, t3, r1, r2, r3):
                sgridtype = recast_gridtype_as_string(self, gridtypei)
                vals = [t1i, t2i, t3i, r1i, r2i, r3i]
                vals2 = write_imag_floats_13e(vals, is_mag_phase)
                [dxr, dyr, dzr, rxr, ryr, rzr,
                 dxi, dyi, dzi, rxi, ryi, rzi] = vals2
                sdt = write_float_12e(dt)
                #if not is_all_zeros:
                if sgridtype == 'G':
                    f06_file.write('0 %12s %6s     %-13s  %-13s  %-13s  %-13s  %-13s  %-s\n'
                                   '  %13s %6s     %-13s  %-13s  %-13s  %-13s  %-13s  %-s\n' % (
                                       sdt, sgridtype, dxr, dyr, dzr, rxr, ryr, rzr,
                                       '', '', dxi, dyi, dzi, rxi, ryi, rzi))
                elif sgridtype in ['S', 'M', 'E']:
                    f06_file.write('0 %12s %6s     %-13s\n'
                                   '  %12s %6s     %-13s\n' % (sdt, sgridtype, dxr, '', '', dxi))
                else:
                    msg = 'nid=%s dt=%s type=%s dx=%s dy=%s dz=%s rx=%s ry=%s rz=%s' % (
                        node_id, dt, sgridtype, t1i, t2i, t3i, r1i, r2i, r3i)
                    raise NotImplementedError(msg)
            f06_file.write(page_stamp % page_num)
            page_num += 1
        return page_num

    #def write_sort2_as_sort2(self, f06_file, page_num, page_stamp, header, words, is_mag_phase):
        #"""TODO: not validated"""
        #node = self.node_gridtype[:, 0]
        #gridtype = self.node_gridtype[:, 1]

        #times = self._times
        ## print(self.data.shape)
        #for inode, (node_id, gridtypei) in enumerate(zip(node, gridtype)):
            ## TODO: for SORT1 pretending to be SORT2
            ##t1 = self.data[:, inode, 0].ravel()
            #t1 = self.data[:, inode, 0].ravel()
            #t2 = self.data[:, inode, 1].ravel()
            #t3 = self.data[:, inode, 2].ravel()
            #r1 = self.data[:, inode, 3].ravel()
            #r2 = self.data[:, inode, 4].ravel()
            #r3 = self.data[:, inode, 5].ravel()
            #if len(r3) != len(times):
                #raise RuntimeError('len(d)=%s len(times)=%s' % (len(r3), len(times)))

            #header[2] = ' POINT-ID = %10i\n' % node_id
            #f06_file.write(''.join(header + words))
            #for dt, t1i, t2i, t3i, r1i, r2i, r3i in zip(times, t1, t2, t3, r1, r2, r3):
                #sgridtype = self.recast_gridtype_as_string(gridtypei)
                #vals = [t1i, t2i, t3i, r1i, r2i, r3i]
                #vals2 = write_imag_floats_13e(vals, is_mag_phase)
                #[dxr, dyr, dzr, rxr, ryr, rzr,
                 #dxi, dyi, dzi, rxi, ryi, rzi] = vals2
                #sdt = write_float_12e(dt)
                ##if not is_all_zeros:
                #if sgridtype == 'G':
                    #f06_file.write('0 %12s %6s     %-13s  %-13s  %-13s  %-13s  %-13s  %-s\n'
                                   #'  %13s %6s     %-13s  %-13s  %-13s  %-13s  %-13s  %-s\n' % (
                                       #sdt, sgridtype, dxr, dyr, dzr, rxr, ryr, rzr,
                                       #'', '', dxi, dyi, dzi, rxi, ryi, rzi))
                #elif sgridtype in ['S', 'M', 'E']:
                    #f06_file.write('0 %12s %6s     %-13s\n'
                                   #'  %12s %6s     %-13s\n' % (sdt, sgridtype, dxr, '', '', dxi))
                #else:
                    #msg = 'nid=%s dt=%s type=%s dx=%s dy=%s dz=%s rx=%s ry=%s rz=%s' % (
                        #node_id, dt, sgridtype, t1i, t2i, t3i, r1i, r2i, r3i)
                    #raise NotImplementedError(msg)
            #f06_file.write(page_stamp % page_num)
            #page_num += 1
        #return page_num

    def write_op2(self, op2_file: BinaryIO, fascii: TextIO, itable: int, new_result,
                  date, is_mag_phase=False, endian: str='>'):
        """writes an OP2"""
        import inspect
        allowed_tables = [
            'OUGV1', 'BOUGV1',
            'OQG1', 'OQMG1',
            'OPG1',
            'OUXY1',
            'OUG1F',
            'OUG1',
            'OUGF1', 'BOUGF1',
            'OAG1', 'OVG1',
        ]
        assert self.table_name in allowed_tables, self.table_name

        frame = inspect.currentframe()
        call_frame = inspect.getouterframes(frame, 2)
        fascii.write('%s.write_op2: %s\n' % (self.__class__.__name__, call_frame[1][3]))

        if itable == -1:
            self._write_table_header(op2_file, fascii, date)
            itable = -3

        #print('nonlinear_factor =', self.nonlinear_factor)
        if self.is_sort1:
            op2_format = endian + b'2i 12f'
        else:
            raise NotImplementedError('SORT2')
        s = Struct(op2_format)

        node = self.node_gridtype[:, 0]
        max_id = node.max()
        if max_id > 99999999:
            raise SixtyFourBitError(f'64-bit OP2 writing is not supported; max id={max_id}')

        gridtype = self.node_gridtype[:, 1]
        #format_table4_1 = Struct(self._endian + b'15i')
        #format_table4_2 = Struct(self._endian + b'3i')

        # table 4 info
        #ntimes = self.data.shape[0]
        nnodes = self.data.shape[1]
        nnodes_device = self.node_gridtype[:, 0] * 10 + self.device_code

        #fdtype = self.data.real.dtype
        #fdtype = 'float32'  # we don't support complex data yet
        #nodedevice_gridtype = np.column_stack([nnodes_device, gridtype])
        #node_gridtype_bytes = nodedevice_gridtype.tobytes()
        #node_gridtype_floats = np.frombuffer(node_gridtype_bytes,
                                             #dtype=fdtype).reshape(nnodes, 2)

        #(2+6) => (node_id, gridtypei, t1i, t2i, t3i, r1i, r2i, r3i)
        ntotal = nnodes * (2 + 12)

        #print('shape = %s' % str(self.data.shape))
        assert nnodes >= 1, nnodes
        assert ntotal > 1, ntotal

        unused_device_code = self.device_code
        fascii.write('  ntimes = %s\n' % self.ntimes)

        #fmt = '%2i %6f'
        #print('ntotal=%s' % (ntotal))
        for itime in range(self.ntimes):
            self._write_table_3(op2_file, fascii, new_result, itable, itime)

            # record 4
            itable -= 1
            header = [4, itable, 4,
                      4, 1, 4,
                      4, 0, 4,
                      4, ntotal, 4,
                      4*ntotal]
            op2_file.write(pack(b'%ii' % len(header), *header))
            fascii.write('r4 [4, 0, 4]\n')
            fascii.write(f'r4 [4, {itable}, 4]\n')
            fascii.write(f'r4 [4, {4*ntotal:d}, 4]\n')

            #datai = self.data[itime, :, :]
            #node_gridtype_data = np.hstack([
                #node_gridtype_floats,
                #datai.real.astype('float32'),
                #datai.imag.astype('float32'),
            #])
            #node_gridtype_data_bytes = node_gridtype_data.tobytes()
            #op2_file.write(node_gridtype_data_bytes)

            t1 = self.data[itime, :, 0]
            t2 = self.data[itime, :, 1]
            t3 = self.data[itime, :, 2]
            r1 = self.data[itime, :, 3]
            r2 = self.data[itime, :, 4]
            r3 = self.data[itime, :, 5]

            for node_id, gridtypei, t1i, t2i, t3i, r1i, r2i, r3i in zip(nnodes_device, gridtype, t1, t2, t3, r1, r2, r3):
                data = [node_id, gridtypei,
                        t1i.real, t2i.real, t3i.real, r1i.real, r2i.real, r3i.real,
                        t1i.imag, t2i.imag, t3i.imag, r1i.imag, r2i.imag, r3i.imag]
                fascii.write('  nid, grid_type, dx, dy, dz, rx, ry, rz = %s\n' % data)
                op2_file.write(s.pack(*data))

            itable -= 1
            header = [4 * ntotal,]
            op2_file.write(pack(b'i', *header))
            fascii.write('footer = %s\n' % header)
            new_result = False
            #header = [
                #4, itable, 4,
                #4, 1, 4,
                #4, 0, 4,
            #]
            #op2_file.write(pack(b'%ii' % len(header), *header))
            #fascii.write('footer2 = %s\n' % header)
        return itable

    #def write_sort2_as_sort2(self, f06_file, page_num, page_stamp, header, words, is_mag_phase):
        #node = self.node_gridtype[:, 0]
        #gridtype = self.node_gridtype[:, 1]

        #times = self._times
        #for inode, (node_id, gridtypei) in enumerate(zip(node, gridtype)):
            ## TODO: for SORT1 pretending to be SORT2
            ##t1 = self.data[:, inode, 0].ravel()
            #t1 = self.data[inode, :, 0]
            #t2 = self.data[inode, :, 1]
            #t3 = self.data[inode, :, 2]
            #r1 = self.data[inode, :, 3]
            #r2 = self.data[inode, :, 4]
            #r3 = self.data[inode, :, 5]
            #if len(r3) != len(times):
                #raise RuntimeError('len(d)=%s len(times)=%s' % (len(r3), len(times)))

            #header[2] = ' POINT-ID = %10i\n' % node_id
            #f06_file.write(''.join(header + words))
            #for dt, t1i, t2i, t3i, r1i, r2i, r3i in zip(times, t1, t2, t3, r1, r2, r3):
                #sgridtype = self.recast_gridtype_as_string(gridtypei)
                #vals = [t1i, t2i, t3i, r1i, r2i, r3i]
                #vals2 = write_imag_floats_13e(vals, is_mag_phase)
                #[dxr, dyr, dzr, rxr, ryr, rzr,
                 #dxi, dyi, dzi, rxi, ryi, rzi] = vals2
                #sdt = write_float_12e(dt)
                ##if not is_all_zeros:
                #if sgridtype == 'G':
                    #f06_file.write('0 %12s %6s     %-13s  %-13s  %-13s  %-13s  %-13s  %-s\n'
                                   #'  %13s %6s     %-13s  %-13s  %-13s  %-13s  %-13s  %-s\n' % (
                                       #sdt, sgridtype, dxr, dyr, dzr, rxr, ryr, rzr,
                                       #'', '', dxi, dyi, dzi, rxi, ryi, rzi))
                #elif sgridtype == 'S':
                    #f06_file.write('0 %12s %6s     %-13s\n'
                                   #'  %12s %6s     %-13s\n' % (sdt, sgridtype, dxr, '', '', dxi))
                #else:
                    #msg = 'nid=%s dt=%s type=%s dx=%s dy=%s dz=%s rx=%s ry=%s rz=%s' % (
                        #node_id, dt, sgridtype, t1i, t2i, t3i, r1i, r2i, r3i)
                    #raise NotImplementedError(msg)
            #f06_file.write(page_stamp % page_num)
            #page_num += 1
        #return page_num

    def plot_freq(self, nids: np.ndarray, dof: int,
                  mag_unit='', scale: float=1.0,
                  labels: Optional[list[str]]=None,
                  yscale: str='log', ifig: int=1):
        if labels is None:
            nnodes = len(nids)
            labels = [None] * nnodes
        import matplotlib.pyplot as plt
        fig = plt.figure(ifig)
        nrows = 2
        ncols = 1
        ax1, ax2 = fig.subplots(nrows, ncols)

        all_nids = self.node_gridtype[:, 0]
        inids = np.searchsorted(all_nids, nids)
        freq = self._times
        idof = dof - 1
        data = self.data[:, inids, :]  # (freq, inid, idof)
        mag, phase = real_imag_to_mag_phase(data)

        assert len(inids) == len(nids)
        tag = ' (excitation)'
        for i, nid, label in zip(count(), nids, labels):
            # ax1.loglog(freq, mag[:, i, idof], label=f'N={nid} dof={dof}')
            # ax2.semilogx(freq, phase[:, i, idof])

            ax1.plot(freq, scale*mag[:, i, idof], label=f'N={nid} dof={dof}{tag} {label}')
            #ax1.plot(freq, mag[:, i, idof], label=f'N={nid} dof={dof}{tag}')
            ax2.plot(freq, phase[:, i, idof])
            tag = ''

        ax1.set_yscale(yscale)
        ax1.legend()
        ax1.grid()
        ax2.grid()
        ax1.set_xlabel('Frequency [Hz]')
        ax2.set_xlabel('Frequency [Hz]')
        if mag_unit:
            ax1.set_ylabel(f'Magnitude [{mag_unit}]')
        else:
            ax1.set_ylabel(f'Magnitude')
        ax2.set_ylabel('Phase [deg]')
        return fig, (ax1, ax2)

    def plot_tf(self, nids_in: int | list[int],
                nids_out: int | list[int],
                dof: int,
                yscale: str='log', ifig: int=1):
        if isinstance(nids_out, integer_types):
            nids_out = [nids_out]
        nout = len(nids_out)
        if isinstance(nids_in, integer_types):
            nids_in = [nids_in] * nout

        import matplotlib.pyplot as plt
        fig = plt.figure(ifig)
        nrows = 2
        ncols = 1
        ax1, ax2 = fig.subplots(nrows, ncols)

        #nids = [nid_in] + nids_out
        all_nids = self.node_gridtype[:, 0]
        inids_in = np.searchsorted(all_nids, nids_in)
        inids_out = np.searchsorted(all_nids, nids_out)
        freq = self._times
        idof = dof - 1
        data_in = self.data[:, inids_in, :]  # (freq, inid, idof)
        data_out = self.data[:, inids_out, :]  # (freq, inid, idof)
        assert data_in.shape == data_out.shape

        nfreq = data_in.shape[0]
        tf_freq_domain = data_out / data_in
        mag, phase = real_imag_to_mag_phase(tf_freq_domain)
        assert mag.shape == (nfreq, nout, 6), mag.shape

        # nid_in, nid_out = nids
        # inid_in, inid_out = inids
        for inid, nid_out, nid_in in zip(count(), nids_out, nids_in):
            ax1.plot(freq, mag[:, inid, idof], label=f'TF({nid_out}/{nid_in}); dof={dof}')
            # ax1.loglog(freq, mag[:, idof], label=f'TF; dof={dof}')
            # ax2.semilogx(freq, phase[:, idof])
            #ax1.plot(freq, mag[:, idof], label=f'TF; dof={dof}')
            ax2.plot(freq, phase[:, inid, idof])
        ax1.set_yscale(yscale)

        ax1.legend()
        ax1.grid()
        ax2.grid()
        ax1.set_xlabel('Frequency [Hz]')
        ax2.set_xlabel('Frequency [Hz]')
        ax1.set_ylabel('Magnitude')
        ax2.set_ylabel('Phase [deg]')
        return fig, (ax1, ax2)


def index_str_to_axis(index_str: str) -> int:
    index_str = index_str.lower().strip()
    if index_str in ['real', 'r']:
        j = 1
    elif index_str in ['imag', 'i']:
        j = 2
    elif index_str in ['mag', 'magnitude', 'm']:
        j = 3
    elif index_str in ['phase', 'p']:
        j = 4
    else:  # pragma: no cover
        raise ValueError(f'index_str={index_str!r}')
    return j

def set_complex_table(cls, data_code, is_sort1: bool, isubcase: int,
                      node_gridtype: np.ndarray, data: np.ndarray, times: np.ndarray):
    assert node_gridtype.ndim == 2, node_gridtype.shape
    ntimes = data.shape[0]
    nnodes = data.shape[1]
    dt = times[0]
    obj = cls(data_code, is_sort1, isubcase, dt)
    obj.node_gridtype = node_gridtype
    obj.data = data

    ntimes = data.shape[0]
    nnodes = data.shape[1]
    #dt = freqs[0]

    obj.ntimes = ntimes
    obj.ntotal = nnodes
    obj.nelements = nnodes
    obj._times = times
    return obj

def pandas_extract_rows(data_frame, ugridtype_str: np.ndarray, index_names: list[str]):
    """removes the t2-t6 for S and E points"""
    import pandas as pd
    letter_dims = [
        ('G', 6),
        ('E', 1),
        ('S', 1),
        ('H', 6),
        ('L', 6),
    ]
    cat_keys = []
    for (letter, dim) in letter_dims:
        if letter not in ugridtype_str:
            continue
        if dim == 1:
            # Note that I'm only keeping every 6th row
            eig = data_frame.xs(letter, level=1).iloc[0::6]
            eig = eig.reset_index()
            #print(eig.columns)
            #print(eig)
            #item = eig.loc[:, 1]
            #item = eig.loc[:, 'Item']
            #print(dir(eig))
            #print(eig.loc)
            #item = eig['Item']
            #print(item)
            #
            # It looks like the error message was changed between pandas 1.5.3 and 2.0.
            # Same issue though.
            #
            # This may be a solution:
            # https://stackoverflow.com/questions/38663150/pivot-table-error1-ndim-categorical-are-not-supported-at-this-time

            # code
            #eig_replace = eig.replace({'Item' : {'t1' : letter}}).set_index(index_names)
            #eig = eig_replace.set_index(index_names)
            #
            # pandas 1.5.3
            # Freq NodeID Item 9.999999747378752e-06      10.0      20.0      30.0      40.0
            # 0       100   t1              0.0+0.0j  0.0+0.0j  0.0+0.0j  0.0+0.0j  0.0+0.0j
            # 1       101   t1              0.0+0.0j  0.0+0.0j  0.0+0.0j  0.0+0.0j  0.0+0.0j
            # 2       102   t1              0.0+0.0j  0.0+0.0j  0.0+0.0j  0.0+0.0j  0.0+0.0j
            # eig_replace:
            #     Freq NodeID Item 9.999999747378752e-06      10.0      20.0      30.0      40.0
            # 0       100   t1              0.0+0.0j  0.0+0.0j  0.0+0.0j  0.0+0.0j  0.0+0.0j
            # 1       101   t1              0.0+0.0j  0.0+0.0j  0.0+0.0j  0.0+0.0j  0.0+0.0j
            # 2       102   t1              0.0+0.0j  0.0+0.0j  0.0+0.0j  0.0+0.0j  0.0+0.0j
            # NotImplementedError: > 1 ndim Categorical are not supported at this time
            #
            # pandas=2.1.1
            # Freq NodeID Item 9.999999747378752e-06      10.0      20.0      30.0      40.0
            # 0       100   t1              0.0+0.0j  0.0+0.0j  0.0+0.0j  0.0+0.0j  0.0+0.0j
            # 1       101   t1              0.0+0.0j  0.0+0.0j  0.0+0.0j  0.0+0.0j  0.0+0.0j
            # 2       102   t1              0.0+0.0j  0.0+0.0j  0.0+0.0j  0.0+0.0j  0.0+0.0j
            # eig_replace:
            #     Freq NodeID Item 9.999999747378752e-06      10.0      20.0      30.0      40.0
            # 0       100   t1              0.0+0.0j  0.0+0.0j  0.0+0.0j  0.0+0.0j  0.0+0.0j
            # 1       101   t1              0.0+0.0j  0.0+0.0j  0.0+0.0j  0.0+0.0j  0.0+0.0j
            # 2       102   t1              0.0+0.0j  0.0+0.0j  0.0+0.0j  0.0+0.0j  0.0+0.0j
            #ValueError: Length mismatch: Expected axis has 3 elements, new values have 1 elements
            #
            eig_replace = eig.replace({'Item' : {'t1' : letter}})
            if pd.__version__ < '2.0':
                try:
                    eig = eig_replace.set_index(index_names)
                except (TypeError, NotImplementedError):
                    print(f'skipping pandas cleanup due to issue with complex {letter} points')
                    return data_frame
            else:
                try:
                    eig = eig_replace.set_index(index_names)
                except (ValueError):
                    print(f'skipping pandas cleanup due to issue with complex {letter} points')
                    return data_frame

        elif dim == 6:
            eig = data_frame.xs(letter, level=1)
        else:
            raise RuntimeError(dim)
        #log.info('eig = %s' % eig)
        cat_keys.append(eig)
    data_frame = pd.concat(cat_keys)
    return data_frame

#class StaticArrayNode(RealTableArray):
    #def __init__(self, data_code, is_sort1, isubcase, dt):
        #RealTableArray.__init__(self, data_code, is_sort1, isubcase, dt)
    #@property
    #def node_ids(self):
        #return self.node_gridtype[:, 0]


#class StaticArrayElement(RealTableArray):
    #def __init__(self, data_code, is_sort1, isubcase, dt):
        #RealTableArray.__init__(self, data_code, is_sort1, isubcase, dt)
    #@property
    #def element_ids(self):
        #return self.node_gridtype[:, 0]


#class TimeArrayNodeSort1(RealTableArray):
    #def __init__(self, data_code, is_sort1, isubcase, dt):
        #RealTableArray.__init__(self, data_code, is_sort1, isubcase, dt)
    #@property
    #def times(self):
        #return self._times
    #@property
    #def node_ids(self):
        #return self.node_gridtype[:, 0]


#class TimeArrayElementSort1(RealTableArray):
    #def __init__(self, data_code, is_sort1, isubcase, dt):
        #RealTableArray.__init__(self, data_code, is_sort1, isubcase, dt)
    #@property
    #def times(self):
        #return self._times
    #@property
    #def element_ids(self):
        #return self.node_gridtype[:, 0]


#class TimeArrayNodeSort2(RealTableArray):
    #def __init__(self, data_code, is_sort1, isubcase, dt):
        #RealTableArray.__init__(self, data_code, is_sort1, isubcase, dt)
    #@property
    #def times(self):
        #return self._times
    #@property
    #def node_ids(self):
        #return self.node_gridtype[:, 0]


#class TimeArrayElementSort2(RealTableArray):
    #def __init__(self, data_code, is_sort1, isubcase, dt):
        #RealTableArray.__init__(self, data_code, is_sort1, isubcase, dt)
    #@property
    #def times(self):
        #return self._times
    #@property
    #def element_ids(self):
        #return self.node_gridtype[:, 0]


#class FrequencyArrayNodeSort2(ComplexTableArray):
    #def __init__(self, data_code, is_sort1, isubcase, dt):
        #ComplexTableArray.__init__(self, data_code, is_sort1, isubcase, dt)
    #@property
    #def frequencies(self):
        #return self._times
    #@property
    #def node_ids(self):
        #return self.node_gridtype[:, 0]


#class FrequencyArrayElementSort2(ComplexTableArray):
    #def __init__(self, data_code, is_sort1, isubcase, dt):
        #ComplexTableArray.__init__(self, data_code, is_sort1, isubcase, dt)
    #@property
    #def frequencies(self):
        #return self._times
    #@property
    #def node_ids(self):
        #return self.node_gridtype[:, 0]


def get_gridtype_str(self, node_gridtype: np.ndarray) -> np.ndarray:
    gridtypes = node_gridtype[:, 1]
    nnodes = len(gridtypes)
    gridtype_str = np.zeros(nnodes, dtype='U1')
    ugridtypes = np.unique(gridtypes)
    for ugridtype in ugridtypes:
        i = np.where(gridtypes == ugridtype)
        gridtype_str[i] = recast_gridtype_as_string(self, ugridtype)
        # gridtype_str2[i] = recast_gridtype_as_string(self, ugridtype)
    return gridtype_str
