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
import copy
from struct import Struct, pack
import warnings
from typing import List

import numpy as np
#from numpy import float32

from pyNastran.op2.result_objects.op2_objects import ScalarObject
from pyNastran.f06.f06_formatting import write_floats_13e, write_imag_floats_13e, write_float_12e
from pyNastran.op2.errors import SixtyFourBitError
from pyNastran.op2.op2_interface.write_utils import set_table3_field

float_types = (float, np.float32)
integer_types = (int, np.int32)

SORT2_TABLE_NAME_MAP = {
    # sort2_name : sort1_name
    'OUGATO2' : 'OUGATO1',
    'OUGCRM2' : 'OUGCRM1',
    'OUGNO2' : 'OUGNO1',
    'OUGPSD2' : 'OUGPSD1',
    'OUGRMS2' : 'OUGRMS1',

    'OVGATO2' : 'OVGATO1',
    'OVGCRM2' : 'OVGCRM1',
    'OVGNO2' : 'OVGNO1',
    'OVGPSD2' : 'OVGPSD1',
    'OVGRMS2' : 'OVGRMS1',

    'OAGATO2' : 'OAGATO1',
    'OAGCRM2' : 'OAGCRM1',
    'OAGNO2' : 'OAGNO1',
    'OAGPSD2' : 'OAGPSD1',
    'OAGRMS2' : 'OAGRMS1',

    'OQGATO2' : 'OQGATO1',
    'OQGCRM2' : 'OQGCRM1',
    'OQGNO2' : 'OQGNO1',
    'OQGPSD2' : 'OQGPSD1',
    'OQGRMS2' : 'OQGRMS1',

    'OQMATO2' : 'OQMATO1',
    'OQMCRM2' : 'OQMCRM1',
    'OQMNO2' : 'OQMNO1',
    'OQMPSD2' : 'OQMPSD1',
    'OQMRMS2' : 'OQMRMS1',

    'OPGATO2' : 'OPGATO1',
    'OPGCRM2' : 'OPGCRM1',
    'OPGNO2' : 'OPGNO1',
    'OPGPSD2' : 'OPGPSD1',
    'OPGRMS2' : 'OPGRMS1',

    #'OUG2' : 'OUG1',
    'OUGV2' : 'OUGV1',
    'OQG2' : 'OQG1',
    'OQMG2' : 'OQMG1',
    'OPG2' : 'OPG1',
    'OPNL2' : 'OPNL1',
    'OUXY2' : 'OUXY1',
    'OQGGF2' : 'OQGGF1',
    'OQGCF2' : 'OQGCF1',
    'OUGF2' : 'OUGF1',
}
SORT1_TABLES = list(SORT2_TABLE_NAME_MAP.values())
SORT1_TABLES.extend([
    'BOUGV1',
    'OUG1F',
    'BOUGF1',
])
SORT2_TABLES = list(SORT2_TABLE_NAME_MAP.keys())

table_name_to_table_code = {
    # displacement (msc/nx)
    'OUGV1' : 1,
    'BOUGV1' : 1,
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

def oug_data_code(table_name, analysis_code,
                  is_sort1=True, is_random=False,
                  random_code=0, title='', subtitle='', label='', is_msc=True):
    sort1_sort_bit = 0 if is_sort1 else 1
    random_sort_bit = 1 if is_random else 0
    sort_method = 1 if is_sort1 else 2
    assert analysis_code != 0, analysis_code
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

    table_code = table_name_to_table_code[table_name]
    sort_code = 1 # TODO: what should this be???

    #table_code = tCode % 1000
    #sort_code = tCode // 1000
    tCode = table_code * 1000 + sort_code

    device_code = 2  # Plot
    approach_code = analysis_code * 10 + device_code
    #print(f'approach_code={approach_code} analysis_code={analysis_code} device_code={device_code}')
    data_code = {
        'nonlinear_factor': None,
        'approach_code' : approach_code,
        'analysis_code' : analysis_code,
        'sort_bits': [0, sort1_sort_bit, random_sort_bit], # real, sort1, random
        'sort_method' : sort_method,
        'is_msc': is_msc,
        #'is_nasa95': is_nasa95,
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
        'num_wide' : 8, # displacement-style table
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

    def get_stats(self, short=False) -> List[str]:
        if not self.is_built:
            return [
                '<%s>; table_name=%r\n' % (self.__class__.__name__, self.table_name),
                '  ntimes: %i\n' % self.ntimes,
                '  ntotal: %i\n' % self.ntotal,
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
        msg.append('  node_gridtype.shape = %s\n' % str(self.node_gridtype.shape).replace('L', ''))
        #msg.append('  gridTypes\n  ')
        msg += self.get_data_code()
        return msg

    @property
    def headers(self) -> List[str]:
        return ['t1', 't2', 't3', 'r1', 'r2', 'r3']

    def _get_headers(self) -> List[str]:
        return self.headers

    def get_headers(self) -> List[str]:
        return self._get_headers()

    def _reset_indices(self):
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
        self._nnodes //= self.ntimes
        #print('ntotal=%s ntimes=%s _nnodes=%s\n' % (self.ntotal, self.ntimes, self._nnodes))
        #if self.ntimes > 1000:
        #    raise RuntimeError(self.ntimes)

        self.itime = 0
        self.itotal = 0
        self.is_built = True

        if self.is_sort1:
            ntimes = self.ntimes
            nnodes = self.ntotal
            ntotal = self.ntotal
            nx = ntimes
            ny = nnodes
            #print("SORT1 ntimes=%s nnodes=%s" % (ntimes, nnodes))
        elif self.is_sort2:
            # flip this to sort1
            ntimes = self.ntotal
            nnodes = self.ntimes
            ntotal = nnodes
            nx = ntimes
            ny = nnodes
            #print("***SORT2 ntotal=%s nnodes=%s ntimes=%s" % (ntotal, nnodes, ntimes))
        else:
            raise RuntimeError('expected sort1/sort2\n%s' % self.code_information())
        self.build_data(ntimes, nnodes, ntotal, nx, ny, self._times_dtype)

    def build_data(self, ntimes, nnodes, ntotal, nx, ny, float_fmt):
        """actually performs the build step"""
        self.ntimes = ntimes
        self._nnodes = nnodes
        self.ntotal = ntotal

        _times = np.zeros(ntimes, dtype=float_fmt)
        int_fmt = 'int32' if self.size == 4 else 'int64'
        node_gridtype = np.zeros((nnodes, 2), dtype=int_fmt)

        #[t1, t2, t3, r1, r2, r3]
        data = np.zeros((nx, ny, 6), self.data_type())
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
        #is_v25 = pd.__version__ >= '0.25'

        headers = self.get_headers()
        #headers = [0, 1, 2, 3, 4, 5]
        #node_gridtype = [self.node_gridtype[:, 0], self.gridtype_str]

        #letter_dims = [
            #('G', 6),
            #('E', 1),
            #('S', 1),
            #('H', 6),
            #('L', 6),
        #]
        ntimes, nnodes = self.data.shape[:2]

        ugridtype_str = np.unique(self.gridtype_str)
        if self.nonlinear_factor not in (None, np.nan):
            #if not self.is_sort1:
                #print("skipping %s becuase it's not SORT1" % self.class_name)
                #return
            column_names, column_values = self._build_dataframe_transient_header()
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

            columns = pd.MultiIndex.from_arrays(column_values, names=column_names)

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
            self.data_frame = data_frame
            #print(self.data_frame)

        else:
            #self.data_frame = pd.Panel(self.data[0, :, :], major_axis=node_gridtype, minor_axis=headers).to_frame()
            #self.data_frame.columns.names = ['Static']
            #self.data_frame.index.names = ['NodeID', 'Type', 'Item']

            #    NodeID Type      t1           t2            t3   r1   r2   r3
            # 0        1    G     0.0     0.000000  0.000000e+00  0.0  0.0  0.0
            # 1        2    G     0.0     0.000000  0.000000e+00  0.0  0.0  0.0
            # 2        3    G     0.0     0.000000  0.000000e+00  0.0  0.0  0.0
            #self.data_frame = pd.DataFrame(self.data[0], columns=headers, index=self.node_gridtype)

            df1 = pd.DataFrame(self.node_gridtype[:, 0])
            df1.columns = ['NodeID']
            df2 = pd.DataFrame(self.gridtype_str)
            df2.columns = ['Type']
            df3 = pd.DataFrame(self.data[0])
            df3.columns = headers
            self.data_frame = df1.join([df2, df3])

            #df1 = pd.DataFrame(self.node_gridtype)
            #df1.columns = ['NodeID', 'Type']
            #df2 = pd.DataFrame(self.data[0])
            #df2.columns = headers
            #self.data_frame = df1.join([df2])
        #print(self.data_frame)

    def finalize(self):
        """
        Calls any OP2 objects that need to do any post matrix calcs
        """
        self.set_as_sort1()
        gridtypes = self.node_gridtype[:, 1]
        nnodes = len(gridtypes)
        self.gridtype_str = np.chararray((nnodes), unicode=True)
        ugridtypes = np.unique(gridtypes)
        for ugridtype in ugridtypes:
            i = np.where(gridtypes == ugridtype)
            self.gridtype_str[i] = self.recast_gridtype_as_string(ugridtype)
        #del self.itotal, self.itime

    def set_as_sort1(self):
        """changes the table into SORT1"""
        #if not self.table_name != 'OQMRMS1':
            #return
        if self.is_sort1:
            return
        #print('set_as_sort1: table_name=%r' % self.table_name)
        try:
            analysis_method = self.analysis_method
        except AttributeError:
            print(self.code_information())
            raise
        #print(self.get_stats())
        #print(self.node_gridtype)
        #print(self.data.shape)
        self.sort_method = 1
        self.sort_bits[1] = 0
        bit0, bit1, bit2 = self.sort_bits
        self.table_name = SORT2_TABLE_NAME_MAP[self.table_name]
        self.sort_code = bit0 + 2*bit1 + 4*bit2
        #print(self.code_information())
        assert self.is_sort1
        if analysis_method != 'N/A':
            self.data_names[0] = analysis_method
            #print(self.table_name_str, analysis_method, self._times)
            setattr(self, self.analysis_method + 's', self._times)
        del self.analysis_method

    def add_sort1(self, dt, node_id, grid_type, v1, v2, v3, v4, v5, v6):
        """unvectorized method for adding SORT1 transient data"""
        assert isinstance(node_id, int) and node_id > 0, 'dt=%s node_id=%s' % (dt, node_id)
        # itotal - the node number
        # itime - the time/frequency step

        # the times/freqs
        self._times[self.itime] = dt
        self.node_gridtype[self.itotal, :] = [node_id, grid_type]
        self.data[self.itime, self.itotal, :] = [v1, v2, v3, v4, v5, v6]
        self.itotal += 1

    def add_sort2(self, dt, node_id, grid_type, v1, v2, v3, v4, v5, v6):
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
        self.node_gridtype[self.itime, :] = [node_id, grid_type]
        self.data[self.itotal, self.itime, :] = [v1, v2, v3, v4, v5, v6]

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
        assert table3[22] == thermal

        n = 0
        from itertools import count
        for i, val, ftable3i in zip(count(), table3, ftable3.decode('ascii')):
            assert val is not None, 'i=%s val=%s ftable3i=%s\n%s' % (i, val, ftable3i, self.get_stats())
            if isinstance(val, integer_types):
                n += 4
                assert ftable3i == 'i', 'analysis_code=%s i=%s val=%s type=%s' % (self.analysis_code, i, val, ftable3i)
            elif isinstance(val, float_types):
                n += 4
                assert ftable3i == 'f', 'analysis_code=%s i=%s val=%s type=%s' % (self.analysis_code, i, val, ftable3i)
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

        table_name = table_name
        analysis_code = 1 # static
        data_code = oug_data_code(table_name, analysis_code,
                                  is_sort1=is_sort1, is_random=is_random,
                                  random_code=random_code,
                                  title=title, subtitle=subtitle, label=label,
                                  is_msc=is_msc)
        data_code['lsdvmns'] = [0] # TODO: ???
        data_code['data_names'] = []

        ntimes = data.shape[0]
        nnodes = data.shape[1]
        dt = None
        obj = cls(data_code, is_sort1, isubcase, dt)
        obj.node_gridtype = node_gridtype
        obj.data = data

        obj.ntimes = ntimes
        obj.ntotal = nnodes
        obj._times = [None]
        obj.is_built = True
        return obj

    @classmethod
    def add_transient_case(cls, table_name, node_gridtype, data, isubcase,
                           times,
                           is_sort1=True, is_random=False, is_msc=True,
                           random_code=0, title='', subtitle='', label=''):

        analysis_code = 6 # transient
        data_code = oug_data_code(table_name, analysis_code,
                                  is_sort1=is_sort1, is_random=is_random,
                                  random_code=random_code, title=title, subtitle=subtitle, label=label,
                                  is_msc=is_msc)
        data_code['data_names'] = ['dt']

        ntimes = data.shape[0]
        nnodes = data.shape[1]
        dt = times[0]
        obj = cls(data_code, is_sort1, isubcase, dt)
        obj.node_gridtype = node_gridtype
        obj.data = data

        obj.ntimes = ntimes
        obj.ntotal = nnodes
        obj.dts = times
        obj._times = times
        obj.is_built = True
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

        analysis_code = 2 # modal
        data_code = oug_data_code(table_name, analysis_code,
                                  is_sort1=is_sort1, is_random=is_random,
                                  random_code=random_code, title=title, subtitle=subtitle, label=label,
                                  is_msc=is_msc)
        #data_code['modes'] = modes
        #data_code['eigns'] = eigenvalues
        #data_code['mode_cycles'] = mode_cycles
        data_code['data_names'] = ['modes', 'eigns', 'mode_cycles']

        ntimes = data.shape[0]
        nnodes = data.shape[1]
        dt = modes[0]
        obj = cls(data_code, is_sort1, isubcase, dt)
        obj.node_gridtype = node_gridtype
        obj.data = data

        obj.modes = modes
        obj.eigns = eigenvalues
        obj.mode_cycles = mode_cycles

        obj.ntimes = ntimes
        obj.ntotal = nnodes
        obj._times = modes
        obj.is_built = True
        return obj

    @property
    def is_real(self):
        return True

    @property
    def is_complex(self):
        return False

    def data_type(self):
        return 'float32'

    def write_op2(self, op2_file, fascii, itable, new_result,
                  date, is_mag_phase=False, endian='>'):
        """writes an OP2"""
        import inspect
        allowed_tables = [
            'OUGV1', 'BOUGV1', 'BOPHIG', 'BOPG1',
            'OUPV1', 'OUXY1', # solution set
            'OQP1', 'OQMG1', 'OQG1', 'OQGV1', 'OPNL1',
            'OPG1', 'OPGV1',
                       'OUGCRM1', 'OUGNO1', 'OUGPSD1', 'OUGRMS1', # disp/vel/acc/eigenvector
            'OAGATO1', 'OAGCRM1', 'OAGNO1', 'OAGPSD1', 'OAGRMS1', # acceleration
                                  'OPGNO1',            'OPGRMS1', # load vector
            'OQGPSD1',
            'OCRPG', 'OCRUG', 'OUG1',
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
        if self.is_sort1:
            op2_format = endian + b'2i6f'
        else:
            raise NotImplementedError('SORT2')
        s = Struct(op2_format)

        node = self.node_gridtype[:, 0]
        gridtype = self.node_gridtype[:, 1]
        max_id = node.max()
        if max_id > 99999999:
            raise SixtyFourBitError(f'64-bit OP2 writing is not supported; max id={max_id}')

        #format_table4_1 = Struct(self._endian + b'15i')
        #format_table4_2 = Struct(self._endian + b'3i')

        # table 4 info
        #ntimes = self.data.shape[0]
        nnodes = self.data.shape[1]
        nnodes_device = self.node_gridtype[:, 0] * 10 + self.device_code

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
            fascii.write('r4 [4, %s, 4]\n' % (itable))
            fascii.write('r4 [4, %i, 4]\n' % (4*ntotal))

            t1 = self.data[itime, :, 0]
            t2 = self.data[itime, :, 1]
            t3 = self.data[itime, :, 2]
            r1 = self.data[itime, :, 3]
            r2 = self.data[itime, :, 4]
            r3 = self.data[itime, :, 5]

            for node_id, gridtypei, t1i, t2i, t3i, r1i, r2i, r3i in zip(nnodes_device, gridtype, t1, t2, t3, r1, r2, r3):
                data = [node_id, gridtypei, t1i, t2i, t3i, r1i, r2i, r3i]
                fascii.write('  nid, grid_type, dx, dy, dz, rx, ry, rz = %s\n' % data)
                op2_file.write(s.pack(*data))

            itable -= 1
            header = [4 * ntotal,]
            op2_file.write(pack(b'i', *header))
            fascii.write('footer = %s\n' % header)
            new_result = False
        return itable

    def write_csv(self, csv_file, is_mag_phase=False):
        name = str(self.__class__.__name__)
        csv_file.write('%s\n' % name)
        headers = ['Node', 'GridType'] + self.headers
        csv_file.write('%s,' * len(headers) % tuple(headers) + '\n')
        node = self.node_gridtype[:, 0]
        gridtype = self.node_gridtype[:, 1]
        itime = 0
        unused_times = self._times

        # sort1 as sort1
        for itime in range(self.ntimes):
            #dt = self._times[itime]
            t1 = self.data[itime, :, 0]
            t2 = self.data[itime, :, 1]
            t3 = self.data[itime, :, 2]
            r1 = self.data[itime, :, 3]
            r2 = self.data[itime, :, 4]
            r3 = self.data[itime, :, 5]
            for node_id, gridtypei, t1i, t2i, t3i, r1i, r2i, r3i in zip(node, gridtype, t1, t2, t3, r1, r2, r3):
                unused_sgridtype = self.recast_gridtype_as_string(gridtypei)
                csv_file.write('%s,' * 9 % (itime, node_id, gridtypei, t1i, t2i, t3i, r1i, r2i, r3i))
                csv_file.write('\n')
        return

    def _write_f06_block(self, words, header, page_stamp, page_num, f06_file, write_words,
                         is_mag_phase=False, is_sort1=True):
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
            sgridtype = self.recast_gridtype_as_string(gridtypei)
            vals = [t1i, t2i, t3i, r1i, r2i, r3i]
            vals2 = write_floats_13e(vals)
            (dx, dy, dz, rx, ry, rz) = vals2
            f06_file.write('%14i %6s     %-13s  %-13s  %-13s  %-13s  %-13s  %s\n' % (
                node_id, sgridtype, dx, dy, dz, rx, ry, rz))
        f06_file.write(page_stamp % page_num)
        return page_num

    def _write_sort1_as_sort2(self, f06_file, page_num, page_stamp, header, words):
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
                sgridtype = self.recast_gridtype_as_string(gridtypei)
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

    def _write_sort1_as_sort1(self, f06_file, page_num, page_stamp, header, words):
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
                sgridtype = self.recast_gridtype_as_string(gridtypei)
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

    def _write_f06_transient_block(self, words, header, page_stamp, page_num, f06_file, write_words,
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


class ComplexTableArray(TableArray):
    """
    complex displacement style table
    """
    def __init__(self, data_code, is_sort1, isubcase, dt):
        TableArray.__init__(self, data_code, is_sort1, isubcase, dt)

    @classmethod
    def add_freq_case(cls, table_name, node_gridtype, data, isubcase,
                      freqs,
                      is_sort1=True, is_random=False, is_msc=True,
                      random_code=0, title='', subtitle='', label=''):

        table_name = 'OUGV1'
        analysis_code = 5 # freq
        data_code = oug_data_code(table_name, analysis_code,
                                  is_sort1=is_sort1, is_random=is_random,
                                  random_code=random_code, title=title, subtitle=subtitle, label=label,
                                  is_msc=is_msc)
        #data_code['modes'] = modes
        #data_code['eigns'] = eigenvalues
        #data_code['mode_cycles'] = mode_cycles
        data_code['data_names'] = ['freq']
        data_code['name'] = 'FREQ'

        ntimes = data.shape[0]
        nnodes = data.shape[1]
        dt = freqs[0]
        obj = cls(data_code, is_sort1, isubcase, dt)
        obj.node_gridtype = node_gridtype
        obj.data = data

        obj.freqs = freqs

        obj.ntimes = ntimes
        obj.ntotal = nnodes
        obj._times = freqs
        obj.is_built = True
        return obj

    def extract_xyplot(self, node_ids, index, index_str):
        index_str = index_str.lower().strip()
        if index_str in ['real', 'r']:
            j = 1
        elif index_str in ['imag', 'i']:
            j = 2
        elif index_str in ['mag', 'magnitude', 'm']:
            j = 3
        elif index_str in ['phase', 'p']:
            j = 4
        else:
            raise ValueError('index_str=%r' % index_str)

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
    def is_real(self):
        return False

    @property
    def is_complex(self):
        return True

    def data_type(self):
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

                if '%' in words[0]:
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
                sgridtype = self.recast_gridtype_as_string(gridtypei)
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

    def write_sort1_as_sort2(self, f06_file, page_num, page_stamp, header, words, is_mag_phase):
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
                sgridtype = self.recast_gridtype_as_string(gridtypei)
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

    def write_op2(self, op2_file, fascii, itable, new_result,
                  date, is_mag_phase=False, endian='>'):
        """writes an OP2"""
        import inspect
        allowed_tables = [
            'OUGV1', 'BOUGV1',
            'OQG1', 'OQMG1',
            'OPG1',
            'OUXY1',
            'OUGF1',
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
            fascii.write('r4 [4, %s, 4]\n' % (itable))
            fascii.write('r4 [4, %i, 4]\n' % (4*ntotal))

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

def pandas_extract_rows(data_frame, ugridtype_str, index_names):
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
            try:
                eig = eig.replace({'Item' : {'t1' : letter}}).set_index(index_names)
            except (TypeError, NotImplementedError):
                print(f'skipping pandas cleanup due to issue with complex {letter} points')
                return data_frame
                #continue
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
