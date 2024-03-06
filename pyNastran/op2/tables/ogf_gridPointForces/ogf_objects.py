from __future__ import annotations
import sys
import inspect
import warnings
from copy import deepcopy
from struct import Struct, pack
from typing import TextIO, Optional, cast, TYPE_CHECKING

import numpy as np
in1d = np.in1d
#in1d = np.in1d if hasattr(np, 'in1d') else getattr(np, 'in')
#from numpy import zeros, unique, array_equal, empty

from cpylog import SimpleLogger

from pyNastran.op2.result_objects.op2_objects import BaseElement, get_times_dtype
from pyNastran.f06.f06_formatting import (
    write_floats_13e, write_floats_13e_long,
    _eigenvalue_header, write_imag_floats_13e)
from pyNastran.op2.vector_utils import (
    transform_force_moment, transform_force_moment_sum, sortedsum1d)
from pyNastran.utils.numpy_utils import integer_types, float_types
from pyNastran.op2.op2_interface.write_utils import set_table3_field
from pyNastran.op2.writer.utils import fix_table3_types

if TYPE_CHECKING:  # pragma: no cover
    from pyNastran.nptyping_interface import (
        NDArrayN3float, NDArray3float, NDArrayN2int, NDArrayNint, NDArrayNfloat)
    from pyNastran.bdf.bdf import BDF, CORD, SimpleLogger


class GridPointForces(BaseElement):
    def __init__(self, data_code, is_sort1, isubcase):
        BaseElement.__init__(self, data_code, isubcase, apply_data_code=True)
        #self.code = [self.format_code, self.sort_code, self.s_code]

        #self.ntimes = 0  # or frequency/mode
        self.ntotal = 0
        self.itotal = 0

    def _write_table_3(self, op2_file, op2_ascii, new_result, itable, itime): #, itable=-3, itime=0):
        frame = inspect.currentframe()
        call_frame = inspect.getouterframes(frame, 2)
        op2_ascii.write('%s.write_table_3: %s\n' % (self.__class__.__name__, call_frame[1][3]))

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
        op2_ascii.write('table_3_header = %s\n' % header)

        approach_code = self.approach_code
        table_code = self.table_code
        isubcase = self.isubcase
        element_type = 0 #self.element_type
        #[
            #'aCode', 'tCode', 'element_type', 'isubcase',
            #'???', '???', '???', 'load_set'
            #'format_code', 'num_wide', 's_code', '???',
            #'???', '???', '???', '???',
            #'???', '???', '???', '???',
            #'???', '???', '???', '???',
            #'???', 'Title', 'subtitle', 'label']
        #random_code = self.random_code
        format_code = self.format_code
        s_code = 0 # self.s_code
        num_wide = self.num_wide
        acoustic_flag = 0
        thermal = 0
        title = b'%-128s' % self.title.encode('ascii')
        subtitle = b'%-128s' % self.subtitle.encode('ascii')
        label = b'%-128s' % self.label.encode('ascii')
        #oCode = 0
        load_set = 0
        #print(self.code_information())

        ftable3 = b'i' * 50 + b'128s 128s 128s'
        field6 = 0
        field7 = 0
        if self.analysis_code == 1:
            field5 = self.lsdvmns[itime]
        elif self.analysis_code == 2:
            ## mode number
            ## mode or cycle .. todo:: confused on the type - F1???
            #self.mode2 = self.add_data_parameter(data, 'mode2', b'i', 7, False)
            #self.cycle = self.add_data_parameter(data, 'cycle', b'f', 7, False)

            field5 = self.modes[itime]
            field6 = self.eigns[itime]
            field7 = self.cycles[itime]
            assert isinstance(field6, float), type(field6)
            assert isinstance(field7, float), type(field7)
            ftable3 = set_table3_field(ftable3, 6, b'f') # field 6
            ftable3 = set_table3_field(ftable3, 7, b'f') # field 7
        elif self.analysis_code == 5:
            field5 = self.freqs[itime]
            ftable3 = set_table3_field(ftable3, 5, b'f') # field 5
        elif self.analysis_code == 6:
            field5 = self.times[itime]
            ftable3 = set_table3_field(ftable3, 5, b'f') # field 5
        elif self.analysis_code == 10:  # nonlinear statics
            field5 = self.lftsfqs[itime]
            ftable3 = set_table3_field(ftable3, 5, b'f') # field 5; load step
        else:
            raise NotImplementedError(self.analysis_code)

        table3 = [
            approach_code, table_code, element_type, isubcase, field5,
            field6, field7, load_set, format_code, num_wide,
            s_code, acoustic_flag, 0, 0, 0,
            0, 0, 0, 0, 0,
            0, 0, thermal, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            0, 0, 0, 0,
            title, subtitle, label,
        ]
        assert table3[22] == thermal

        table3 = fix_table3_types(table3, size=4)
        data = [584] + table3 + [584]
        fmt = b'i' + ftable3 + b'i'
        #print(fmt)
        #print(data)
        #f.write(pack(fascii, '%s header 3c' % self.table_name, fmt, data))
        op2_ascii.write('%s header 3c = %s\n' % (self.table_name, data))
        op2_file.write(pack(fmt, *data))


class RealGridPointForcesArray(GridPointForces):
    """
                                       G R I D   P O I N T   F O R C E   B A L A N C E
       POINT-ID  ELEMENT-ID   SOURCE        T1       T2    T3            R1   R2   R3
    0     13683        3736  TRIAX6    4.996584E+00  0.0   1.203093E+02  0.0  0.0  0.0
          13683        3737  TRIAX6   -4.996584E+00  0.0  -1.203093E+02  0.0  0.0  0.0
          13683              *TOTALS*  6.366463E-12  0.0  -1.364242E-12  0.0  0.0  0.0

    """
    def __init__(self, data_code, is_sort1, isubcase, dt):
        GridPointForces.__init__(self, data_code, is_sort1, isubcase)
        #self.code = [self.format_code, self.sort_code, self.s_code]

        #self.ntimes = 0  # or frequency/mode
        self.ntotal = 0
        self.itotal = 0

        # do the element_names/node_element vectors change with the time step
        self.is_unique = False
        self.element_names = []

        #self.ielement = 0
        #self.nelements = 0  # result specific
        #self.nnodes = None

    def finalize(self) -> None:
        """required so the OP2 writer works..."""
        self.format_code = 1

    def __pos__(self) -> RealGridPointForcesArray:
        """positive; +a"""
        return self

    def __neg__(self) -> RealGridPointForcesArray:
        """negative; -a"""
        new_table = deepcopy(self)
        new_table.data *= -1.0
        return new_table

    def __add__(self, table: RealGridPointForcesArray) -> RealGridPointForcesArray:
        """a + b"""
        if isinstance(table, RealGridPointForcesArray):
            self._check_math(table)
            new_data = self.data + table.data
            #_fix_min_max(new_data)
        elif isinstance(table, (integer_types, RealGridPointForcesArray)):
            new_data = self.data + table
        else:
            raise TypeError(table)
        new_table = deepcopy(self)
        new_table.data = new_data
        return new_table
    # __radd__: reverse order adding (b+a)
    def __iadd__(self, table: RealGridPointForcesArray) -> RealGridPointForcesArray:
        """inplace adding; a += b"""
        if isinstance(table, RealGridPointForcesArray):
            self._check_math(table)
            self.data += table.data
            #_fix_min_max(self.data)
        elif isinstance(table, (integer_types, float_types)):
            self.data -= table
        else:
            raise TypeError(table)
        return self

    def __sub__(self, table: RealGridPointForcesArray):
        """a - b"""
        if isinstance(table, RealGridPointForcesArray):
            self._check_math(table)
            new_data = self.data - table.data
        elif isinstance(table, (integer_types, float_types)):
            new_data = self.data - table
        else:
            raise TypeError(table)
        new_table = deepcopy(self)
        new_table.data = new_data
        return new_table

    def __mul__(self, table: RealGridPointForcesArray):
        """a * b"""
        if isinstance(table, RealGridPointForcesArray):
            self._check_math(table)
            new_data = self.data * table.data
        elif isinstance(table, (integer_types, float_types)):
            new_data = self.data * table
        else:
            raise TypeError(table)
        new_table = deepcopy(self)
        new_table.data = new_data
        return new_table

    def __truediv__(self, table: RealGridPointForcesArray):
        """a / b"""
        if isinstance(table, RealGridPointForcesArray):
            self._check_math(table)
            new_data = self.data / table.data
        elif isinstance(table, (integer_types, float_types)):
            new_data = self.data / table
        else:
            raise TypeError(table)
        new_table = deepcopy(self)
        new_table.data = new_data
        return new_table

    def _check_math(self, table: RealGridPointForcesArray) -> None:
        """verifies that the shapes are the same"""
        assert self.ntimes == table.ntimes, f'ntimes={self.ntimes} table.times={table.ntimes}'
        assert self.ntotal == table.ntotal, f'ntotal={self.ntotal} table.ntotal={table.ntotal}'
        assert self.node_element.shape == table.node_element.shape, f'node_element.shape={self.node_element.shape} table.element_node.shape={table.node_element.shape}'
        assert self.element_names.shape == table.element_names.shape, f'element_names.shape={self.element_names.shape} table.element_names.shape={table.element_names.shape}'
        assert self.data.shape == table.data.shape, f'data.shape={self.data.shape} table.data.shape={table.data.shape}'

    @property
    def is_real(self) -> bool:
        return True

    @property
    def is_complex(self) -> bool:
        return False

    def _reset_indices(self) -> None:
        self.itotal = 0
        #self.ielement = 0

    @property
    def element_name(self) -> str:
        headers = [name.strip() for name in np.unique(self.element_names) if name.strip()]
        #headers = np.unique(self.element_names)
        return str(', '.join(headers))

    def build(self) -> None:
        """sizes the vectorized attributes of the RealGridPointForcesArray"""
        #print("self.ielement = %s" % self.ielement)
        #print('ntimes=%s nelements=%s ntotal=%s' % (self.ntimes, self.nelements, self.ntotal))
        assert self.ntimes > 0, 'ntimes=%s' % self.ntimes
        #assert self.nelements > 0, 'nelements=%s' % self.nelements
        assert self.ntotal > 0, 'ntotal=%s' % self.ntotal
        #if self.ntotal != max(self._ntotals) or self.ntotal != min(self._ntotals):
            #raise ValueError('RealGridPointForcesArray: ntotal=%s _ntotals=%s' % (
                #self.ntotal, self._ntotals))

        self.is_unique = False
        if self.ntotal != min(self._ntotals) or 1:
            self.ntotal = max(self._ntotals)
            self.is_unique = True
        #self.names = []

        #self.nnodes = nnodes_per_element
        #self.nelements //= nnodes_per_element
        self.itime = 0
        #self.ielement = 0
        self.itotal = 0
        #self.ntimes = 0
        #self.nelements = 0

        #print("***name=%s ntimes=%s ntotal=%s" % (
            #self.element_names, self.ntimes, self.ntotal))
        dtype, idtype, fdtype = get_times_dtype(self.nonlinear_factor, self.size, self.analysis_fmt)
        self._times = np.zeros(self.ntimes, dtype=self.analysis_fmt)

        assert self.ntotal < 2147483647, self.ntotal # max int
        if self.is_unique:
            assert isinstance(self.ntotal, integer_types), 'ntotal=%r type=%s' % (self.ntotal, type(self.ntotal))
            self.node_element = np.zeros((self.ntimes, self.ntotal, 2), dtype=idtype)
            self.element_names = np.empty((self.ntimes, self.ntotal), dtype='U8')
        else:
            self.node_element = np.zeros((self.ntotal, 2), dtype=idtype)
            self.element_names = np.empty(self.ntotal, dtype='U8')

        #[t1, t2, t3, r1, r2, r3]
        self.data = np.zeros((self.ntimes, self.ntotal, 6), dtype=fdtype)

    def build_dataframe(self) -> None:
        """
        major-axis - the axis

        mode              1     2   3
        freq              1.0   2.0 3.0
        nodeID ElementID Item
        1      2         T1
                         T2
                         ...

        major_axis / top = [
            [1, 2, 3],
            [1.0, 2.0, 3.0]
        ]
        minor_axis / headers = [T1, T2, T3, R1, R2, R3]
        name = mode
        """
        #name = self.name
        headers = self.get_headers()
        if self.is_unique:
            data_frame = self._build_unique_dataframe(headers)
        else:
            import pandas as pd
            node_element = [self.node_element[:, 0], self.node_element[:, 1]]
            if self.nonlinear_factor not in (None, np.nan):
                column_names, column_values = self._build_dataframe_transient_header()
                data_frame = pd.Panel(
                    self.data, items=column_values,
                    major_axis=node_element, minor_axis=headers).to_frame()
                data_frame.columns.names = column_names
                data_frame.index.names = ['NodeID', 'ElementID', 'Item']
            else:
                data_frame = pd.Panel(
                    self.data,
                    major_axis=node_element, minor_axis=headers).to_frame()
                data_frame.columns.names = ['Static']
                data_frame.index.names = ['NodeID', 'ElementID', 'Item']
            #print(self.data_frame)
        self.data_frame = data_frame

    def _build_unique_dataframe(self, headers: list[str]):
        import pandas as pd
        ntimes = self.data.shape[0]
        nnodes = self.data.shape[1]
        #nvalues = ntimes * nnodes
        node_element = self.node_element.reshape((ntimes * nnodes, 2))
        if self.nonlinear_factor not in (None, np.nan):
            column_names, column_values = self._build_dataframe_transient_header()
            #column_names = column_names[0]
            #column_values = column_values[0]

            column_values2 = []
            for value in column_values:
                values2 = []
                for valuei in value:
                    values = np.ones(nnodes) * valuei
                    values2.append(values)
                values3 = np.vstack(values2).ravel()
                column_values2.append(values3)
            df1 = pd.DataFrame(column_values2).T
            df1.columns = column_names
            #df1.columns.names = column_names
            #self.data_frame.columns.names = column_names

            df2 = pd.DataFrame(node_element)
            df2.columns = ['NodeID', 'ElementID']
            df3 = pd.DataFrame(self.element_names.ravel())
            df3.columns = ['ElementType']

            dfs = [df2, df3]
            for i, header in enumerate(headers):
                df = pd.DataFrame(self.data[:, :, i].ravel())
                df.columns = [header]
                dfs.append(df)
            data_frame = df1.join(dfs)
            #print(self.data_frame)
        else:
            df1 = pd.DataFrame(node_element)
            df1.columns = ['NodeID', 'ElementID']
            df2 = pd.DataFrame(self.element_names[0, :])
            df2.columns = ['ElementType']
            df3 = pd.DataFrame(self.data[0])
            df3.columns = headers
            data_frame = df1.join([df2, df3])
            #print(data_frame)
        return data_frame

    def __eq__(self, table) -> bool:  # pragma: no cover
        is_valid = False
        try:
            self.assert_equal(table)
            is_valid = True
        except Exception:
            pass
        return is_valid

    def assert_equal(self, table, rtol=1.e-5, atol=1.e-8):
        self._eq_header(table)
        assert self.is_sort1 == table.is_sort1
        if not np.array_equal(self.node_element, table.node_element) and np.array_equal(self.element_names, table.element_names):
            assert self.node_element.shape == table.node_element.shape, 'node_element shape=%s table.shape=%s' % (self.node_element.shape, table.node_element.shape)
            assert self.element_names.shape == table.element_names.shape, 'element_names shape=%s table.shape=%s' % (self.element_names.shape, table.element_names.shape)

            msg = 'table_name=%r class_name=%s\n' % (self.table_name, self.__class__.__name__)
            msg += '%s\n' % str(self.code_information())
            msg += '(Nid, Eid, EName)\n'
            for (nid1, eid1), ename1, (nid2, eid2), ename2 in zip(self.node_element, self.element_names,
                                                                  table.element_names, table.element_names):
                msg += '(%s, %s, %s)    (%s, %s, %s)\n' % (nid1, eid1, ename1, nid2, eid2, ename2)
            print(msg)
            raise ValueError(msg)

        atols = []
        rtols = []
        if self.is_unique:
            # node_element varies with time
            if not np.array_equal(self.data, table.data):
                msg = 'table_name=%r class_name=%s\n' % (self.table_name, self.__class__.__name__)
                msg += '%s\n' % str(self.code_information())
                i = 0
                for itime in range(self.ntimes):
                    #print('node_element = ', self.node_element)
                    #print('shape = ', self.node_element.shape)
                    msg += '#i, Nid, Eid, Name (itime=%s)\n' % itime
                    for ie, e in enumerate(self.node_element[itime, :, :]):
                        #print('e = ', e)
                        (nid, eid) = e
                        ename1 = self.element_names[itime, ie]
                        ename2 = self.element_names[itime, ie]
                        t1 = self.data[itime, ie, :]
                        t2 = table.data[itime, ie, :]

                        if not np.allclose(t1, t2, rtol=rtol, atol=atol):
                            (t11, t21, t31, r11, r21, r31) = t1
                            (t12, t22, t32, r12, r22, r32) = t2
                            inonzero = np.where(t1 != 0.)[0]
                            atoli = np.abs(t2 - t1).max()
                            rtoli = np.abs(t2[inonzero] / t1[inonzero]).max()

                            pre_msg = '(%s, %s, %s, %s)    ' % (ie, nid, eid, ename1)
                            msg += '%s(%s, %s, %s, %s, %s, %s)\n%s(%s, %s, %s, %s, %s, %s)\n' % (
                                pre_msg,
                                t11, t21, t31, r11, r21, r31,
                                ' ' * len(pre_msg),
                                t12, t22, t32, r12, r22, r32)
                            i += 1
                            atols.append(atoli)
                            rtols.append(rtoli)
                            if i > 30:
                                print(atols)
                                msg += 'atol.max() = %s\n' % max(atols)
                                msg += 'rtol.max() = %s\n' % max(rtols)
                                print(msg)
                                raise ValueError(msg)
                    #print(msg)
                    if i > 0:
                        msg += 'atol.max() = %s\n' % max(atols)
                        msg += 'rtol.max() = %s\n' % max(rtols)
                        #print(msg)
                        raise ValueError(msg)
        else:
            # node_element does not vary with time
            if not np.array_equal(self.data, table.data):
                msg = 'table_name=%r class_name=%s\n' % (self.table_name, self.__class__.__name__)
                msg += '%s\n' % str(self.code_information())
                i = 0
                for itime in range(self.ntimes):
                    for ie, e in enumerate(self.node_element):
                        (nid, eid) = e
                        ename1 = self.element_names[ie]
                        ename2 = self.element_names[ie]
                        t1 = self.data[itime, ie, :]
                        t2 = table.data[itime, ie, :]
                        (t11, t21, t31, r11, r21, r31) = t1
                        (t12, t22, t32, r12, r22, r32) = t2

                        if not np.allclose(t1, t2, rtol=rtol, atol=atol):
                            pre_msg = '(%s, %s, %s)    ' % (nid, eid, ename1)
                            msg += '%s(%s, %s, %s, %s, %s, %s)\n%s(%s, %s, %s, %s, %s, %s)\n' % (
                                pre_msg,
                                t11, t21, t31, r11, r21, r31,
                                ' ' * len(pre_msg),
                                t12, t22, t32, r12, r22, r32
                            )
                            i += 1
                            if i > 10:
                                msg += 'atol.max() = %s\n' % max(atols)
                                msg += 'rtol.max() = %s\n' % max(rtols)
                                print(msg)
                                raise ValueError(msg)
                    if i > 0:
                        msg += 'atol.max() = %s\n' % max(atols)
                        msg += 'rtol.max() = %s\n' % max(rtols)
                        #print(msg)
                        raise ValueError(msg)
        return True

    def extract_freebody_loads(self, eids: NDArrayNint,
                               coord_out: CORD,
                               coords: dict[int, CORD],
                               nid_cd: NDArrayN2int,
                               icd_transform: dict[int, NDArrayNint],
                               itime: int=0, debug: bool=False,
                               log: Optional[SimpleLogger]=None):
        """
        Extracts Patran-style freebody loads.  Freebody loads are the
        external loads.

        Parameters
        ----------
        eids : (Nelements, ) int ndarray
            all the elements to consider
        coord_out : CORD2R()
            the output coordinate system
        coords : dict[int] = CORDx
            all the coordinate systems
            key : int
            value : CORDx
        nid_cd : (Nnodes, 2) int ndarray
            the (BDF.point_ids, cd) array
        icd_transform : dict[cd] = (Nnodesi, ) int ndarray
            the mapping for nid_cd
        itime : int; default=0
            the time to extract loads for
        debug : bool; default=False
            debugging flag
        log : log; default=None
            a log object that gets used when debug=True

        Returns
        -------
        force_out : (Nnodes, 3) float ndarray
            the ith float components in the coord_out coordinate frame
        moment_out : (Nnodes, 3) float ndarray
            the ith moment components about the summation point in the
            coord_out coordinate frame

        .. todo:: doesn't seem to handle cylindrical/spherical systems
        .. warning:: not done
        """
        nid_cd = _get_nid_cd_from_nid_cp_cd(nid_cd)
        eids = np.asarray(eids)
        #nids = np.asarray(nids)

        # todo handle multiple values for itime
        gpforce_nids = self.node_element[itime, :, 0]
        gpforce_eids = self.node_element[itime, :, 1]
        # TODO: remove 0s in gpforce_nids/gpforce_eids to handle transient results
        #       be careful of the sum row

        assert isinstance(eids[0], integer_types), type(eids[0])

        is_in = in1d(gpforce_eids, eids, assume_unique=False)
        irange = np.arange(len(gpforce_nids), dtype='int32')[is_in]
        nids = gpforce_nids[irange]

        if debug and log is not None:
            log.debug('gpforce_eids = %s' % gpforce_eids[is_in])
            log.debug('nids = %s' % gpforce_nids[irange])
            log.debug('eids = %s' % gpforce_eids[irange])

        try:
            is_in3 = in1d(nid_cd[:, 0], nids, assume_unique=False)
        except IndexError:
            msg = 'nids_cd=%s nids=%s' % (nid_cd, nids)
            raise IndexError(msg)

        force_global = self.data[itime, irange, :3]
        moment_global = self.data[itime, irange, 3:]

        #data_global = sortedsum1d(nids, self.data[itime, irange, :], axis=0)
        #force_global2 = data_global[:, :3]
        #moment_global2 = data_global[:, 3:]
        force_global = sortedsum1d(nids, force_global)
        moment_global = sortedsum1d(nids, moment_global)
        #print(force_global)
        #print(force_global2)
        #assert np.array_equal(force_global, force_global2)
        #assert np.array_equal(moment_global, moment_global2)

        force, moment = transform_force_moment(
            force_global, moment_global,
            coord_out, coords, nid_cd[is_in3, :],
            icd_transform,
            xyz_cid0=None, summation_point_cid0=None,
            consider_rxf=False,
            debug=debug, log=log)
        return force, moment

    def extract_interface_loads(self,
                                nids: NDArrayNint,
                                eids: NDArrayNint,
                                coord_out: CORD,
                                coords: dict[int, CORD],
                                nid_cd: NDArrayN2int,
                                icd_transform: dict[int, CORD],
                                xyz_cid0: NDArrayN3float,
                                summation_point: Optional[NDArray3float]=None,
                                consider_rxf: bool=True,
                                itime: int=0,
                                assume_sorted: bool=False,
                                stop_on_nan: bool=False,
                                debug: bool=False,
                                log: Optional[SimpleLogger]=None,
                                idtype: str='int32') -> tuple[NDArray3float, NDArray3float]:
        """
        Extracts Patran-style interface loads.  Interface loads are the
        internal loads at a cut.

        Parameters
        ----------
        nids : (Nnodes, ) int ndarray
            all the nodes to consider; must be sorted
        eids : (Nelements, ) int ndarray
            all the elements to consider; must be sorted
        coord_out : CORD2R()
            the output coordinate system
        coords : dict[int] = CORDx
            all the coordinate systems
            key : int
            value : CORDx
        nid_cd : (Nnodes, 2) int ndarray
            the (BDF.point_ids, cd) array
        icd_transform : dict[cd] = (Nnodesi, ) int ndarray
            the mapping for nid_cd
        xyz_cid0 : (nnodes + nspoints + nepoints, 3) ndarray
            the grid locations in coordinate system 0
        summation_point0 : (3, ) float ndarray; default=None
            None : no load summation
            array : the summation point in the global frame
        consider_rxf : bool; default=True
            considers the r x F term
        itime : int; default=0
            the time to extract loads for
        debug : bool; default=False
            debugging flag
        logger : logger; default=None
            a logger object that gets used when debug=True
        assume_sorted : bool; default=False
            sorts the nodes/elements if they're not
            is sorting required?

        Returns
        -------
        force_out_sum : (3, ) float ndarray
            the sum of forces in the coord_out coordinate frame
        moment_out_sum : (3, ) float ndarray
            the sum of moments about the summation point in the
            coord_out coordinate frame

        .. todo:: doesn't seem to handle cylindrical/spherical systems
        .. todo:: Add support for:
                  2D output style:
                    - This would allow for shell problems to have loads
                      applied in the plane of the shells
                    - This would require normals
                  1D output style:
                    - Make loads in the direction of the element
                  This process can't be done for 0D or 3D elements

        """
        # force_out : (Nnodes, 3) float ndarray
        #     the ith float components in the coord_out coordinate frame
        # moment_out : (Nnodes, 3) float ndarray
        #     the ith moment components about the summation point in the
        #     coord_out coordinate frame
        _check_array(nids, 'int', 1)
        _check_array(nid_cd, 'int', 2)
        _check_array(eids, 'int', 1)
        _check_array(xyz_cid0, 'float', 2)

        unused_force_out, unused_moment_out, force_out_sum, moment_out_sum = self._extract_interface_loads(
            nids, eids, coord_out, coords, nid_cd, icd_transform, xyz_cid0,
            summation_point=summation_point, consider_rxf=consider_rxf,
            itime=itime, assume_sorted=assume_sorted,
            stop_on_nan=stop_on_nan,
            debug=debug, log=log, idtype=idtype)
        return force_out_sum, moment_out_sum

    def _extract_interface_loads(self,
                                nids: NDArrayNint,
                                eids: NDArrayNint,
                                coord_out: CORD,
                                coords: dict[int, CORD],
                                nid_cd: NDArrayN2int,
                                icd_transform: dict[int, CORD],
                                xyz_cid0: NDArrayN3float,
                                summation_point: Optional[NDArray3float]=None,
                                consider_rxf: bool=True,
                                itime: int=0,
                                assume_sorted: bool=False,
                                debug: bool=False,
                                stop_on_nan: bool=False,
                                log: Optional[SimpleLogger]=None,
                                idtype: str='int32') -> tuple[NDArrayN3float, NDArrayN3float, NDArray3float, NDArray3float]:
        """see ``extract_interface_loads``"""
        fdtype = self.data.dtype
        nid_cd = _get_nid_cd_from_nid_cp_cd(nid_cd)
        if summation_point is None:
            summation_point = np.array([0., 0., 0.])
        else:
            summation_point = np.asarray(summation_point)
        #assert coord_in.Type == 'R', 'Only rectangular coordinate systems are supported; coord_in=\n%s' % str(coord_in)
        #assert coord_out.Type == 'R', 'Only rectangular coordinate systems are supported; coord_out=\n%s' % str(coord_out)
        assert eids is not None, eids
        assert nids is not None, nids
        assert isinstance(itime, integer_types), type(itime)
        eids = np.asarray(eids)
        nids = np.asarray(nids)
        if not assume_sorted:
            eids.sort()
            nids.sort()
        out = self._extract_interface_loads_helper(
            nids, eids, log, idtype, fdtype,
            itime=itime, debug=debug)
        gpforce_nids, gpforce_eids, irange, force_out, moment_out = out
        if len(irange) == 0:
            force_out_sum = np.full(3, np.nan, dtype=force_out.dtype)
            moment_out_sum = np.full(3, np.nan, dtype=force_out.dtype)
            if stop_on_nan:
                raise RuntimeError(f'no nodes/elements in intersect; summation_point={summation_point}')

            return force_out, moment_out, force_out_sum, moment_out_sum

        if debug and log is not None:
            f06_filename = f'grid_point_forcesi_itime={itime}.debug.f06'
            with open(f06_filename, 'w') as f06_file:
                self.write_f06_time(f06_file, itime=itime, i=irange)

            #log.debug('gpforce_eids =' % gpforce_eids[is_in])
            log = cast(SimpleLogger, log)
            log.debug('gpforce_nids = %s' % gpforce_nids[irange])
            log.debug('gpforce_eids = %s' % gpforce_eids[irange])

        # get analysis coordinate systems
        try:
            is_in_cd = in1d(nid_cd[:, 0], nids, assume_unique=False)
        except IndexError:
            msg = 'nids_cd=%s nids=%s' % (nid_cd, nids)
            raise IndexError(msg)

        nid_cd_used = nid_cd[is_in_cd, :]
        nids_used = nid_cd_used[:, 0]
        gp_nids_used = gpforce_nids[irange]
        isort = np.searchsorted(nids_used, gp_nids_used)

        force_global = self.data[itime, irange, :3]
        moment_global = self.data[itime, irange, 3:]

        # update only the nodes that are in the CD coordinate systems (is_in_cd)
        #force_out_sum = np.full(3, np.nan, dtype=fdtype)
        #moment_out_sum = np.full(3, np.nan, dtype=fdtype)
        nid_cd_used_sorted = nid_cd[is_in_cd, :][isort]
        xyz_cid0_used_sorted = xyz_cid0[is_in_cd, :][isort]
        force_out, moment_out, force_out_sum, moment_out_sum = transform_force_moment_sum(
            force_global, moment_global,
            coord_out, coords, nid_cd_used_sorted,
            icd_transform,
            xyz_cid0_used_sorted, summation_point_cid0=summation_point,
            consider_rxf=consider_rxf,
            debug=debug, log=log)

        if not np.isfinite(force_out_sum[0]) and stop_on_nan:
            raise RuntimeError(f'no nodes/elements in intersect; summation_point={summation_point}')

        return force_out, moment_out, force_out_sum, moment_out_sum

    def _extract_interface_loads_helper(
        self, nids: NDArrayNint, eids: NDArrayNint,
        log: SimpleLogger,
        idtype: str,
        fdtype: str,
        itime:int=0,
        debug: bool=False) -> tuple[NDArrayNint, NDArrayNint,
                                    NDArrayN3float, NDArrayN3float]:
        """
        Returns
        -------
        gpforce_nids : (nnode, ) int array
            all the nodes
        gpforce_eids : (nnode, ) int array
            all the elements
        irange : (nnode_filtered, ) int array
            the index into GPFORCE data of the rows to include
        force_out : (nstations, 3) float array
            the forces
        moment_out : (nstations, 3) float array
            the moments

        """
        force_out = np.full((0, 3), np.nan, dtype=fdtype)
        moment_out = np.full((0, 3), np.nan, dtype=fdtype)
        #force_out_sum = np.full(3, np.nan, dtype=fdtype)
        #moment_out_sum = np.full(3, np.nan, dtype=fdtype)
        # TODO: Handle multiple values for itime
        #       Is this even possible?
        gpforce_nids = self.node_element[itime, :, 0]
        gpforce_eids = self.node_element[itime, :, 1]
        # TODO: Remove 0s in gpforce_nids/gpforce_eids to handle transient results
        #       be careful of the sum row.
        #       Do I even need to do this?

        assert isinstance(eids[0], integer_types), type(eids[0])
        assert isinstance(nids[0], integer_types), type(nids[0])
        # filter out rows not in the node set
        is_in = in1d(gpforce_nids, nids, assume_unique=False)
        if not np.any(is_in):
            msg = 'no nodes found\n'
            if log:
                log.warning(msg)
            else:
                warnings.warn(msg)
            irange = np.array([], dtype='int32')
            return gpforce_nids, gpforce_eids, irange, force_out, moment_out

        # filter out rows not in the element set
        is_in2 = in1d(gpforce_eids[is_in], eids, assume_unique=False)
        if not np.any(is_in2):
            msg = 'no elements found\n'
            log.warning(msg)
            irange = np.array([], dtype='int32')
            return gpforce_nids, gpforce_eids, irange, force_out, moment_out

        igpforce_nids = np.arange(len(gpforce_nids), dtype=idtype)
        irange = igpforce_nids[is_in][is_in2]
        if irange.size == 0:
            msg = (
                'no nodes/elements found\n'
                f'eids={eids}\n'
                f'gpforce_nids={gpforce_nids}\n'
                f'gpforce_eids={gpforce_eids}\n'
                f'gpforce_nids_found={gpforce_nids[is_in][is_in2]}\n'
                f'gpforce_eids_found={gpforce_eids[is_in][is_in2]}\n'
            )
            log.warning(msg)
            #raise RuntimeError(msg)
        return gpforce_nids, gpforce_eids, irange, force_out, moment_out

    def find_centroid_of_load(self, f, m):  # pragma: no cover
        """
        Mx = ry*Fz - rz*Fy
        My = rz*Fx - rx*Fz
        Mz = rx*Fy - ry*Fx

        {M} = [F]{r}
        [F] = [
            [  0, -Fy, Fz],
            [-Fz,   0, Fx],
            [ Fy, -Fx,  0],
        ]
        {r} = [F]^-1 {M}

        When the determinant of [F] is nonzero (2D):
           Life is easy

        When the determinant of [F] is zero:
        When Fx != Fy != Fz and they don't equal 0
        there are 3 solutions:
            where M=[0, 0, 0]

        det([F]) = 0:
           [F]{x} = [lambda]{x}
        where one of the eigenvalues is 0? (the trivial case)
        and

        However, [F] is singular, so let rx=0:
        Mx = ry*Fz - rz*Fy
        My = rz*Fx
        Mz = -ry*Fx
        let Fx=0, so ry, rz != 0, but My=Mz=0
        -> ry = rz*Fy/Fz
        let rz = 1
        -> ry = Fy/Fz
        <0, Fy/Fz, 1>
        """
        raise NotImplementedError()


    def shear_moment_diagram(self, # model,
                             nids: np.ndarray,
                             xyz_cid0: np.ndarray,
                             nid_cd: NDArrayN2int,
                             icd_transform: dict[int, NDArrayNint],

                             eids: np.ndarray,
                             element_centroids_cid0: NDArrayN3float,
                             stations: NDArrayNfloat,
                             coords: dict[int, CORD],
                             coord_out: CORD,
                             iaxis_march: Optional[NDArray3float]=None,
                             itime: int=0,
                             icoord: int=None,
                             idir: int=0,
                             nodes_tol: Optional[float]=None,
                             stop_on_nan: bool=False,
                             debug: bool=False,
                             log: Optional[SimpleLogger]=None) -> tuple[NDArray3float,
                                                                        NDArray3float,
                                                                        dict[int, CORD],
                                                                        NDArrayNint, NDArrayNint]:
        """
        Computes a series of forces/moments at various stations along a
        structure.

        Parameters
        ----------
        xyz_cid0 : (Nnodes, 3) float ndarray
            all the nodes in the model xyz position in the global frame
        eids : (Nelements, ) int ndarray
            an array of element ids to consider
        nids : (Nnodes, ) int ndarray
            an array of node ids corresponding to xyz_cid0
        icd_transform : dict[cd] = (Nnodesi, ) int ndarray
            the mapping for nid_cd
        element_centroids_cid0 : (Nelements, 3) float ndarray
            an array of element centroids corresponding to eids
        coords : dict[int] = CORDx
            all the coordinate systems
            key : int
            value : CORDx
        nid_cd : (Nnodes, 2) int ndarray
            the (BDF.point_ids, cd) array
        stations : (nstations, ) float ndarray
            The station to sum forces/moments about, where a station is
            the value of coord_out in the idir (e.g., for station=1, cid=0,
            idir=0, then x=1).
            It is necessary that the spacing is constant.
            Stations should be sorted (negative to positive), but it not
            necessary (it makes it easier to see a monotonic change in
            the number of nodes/elements).
            Try to avoid picking exactly on symmetry planes/boundaries
            of elements or nodes.
        coord_out : CORD2R()
            the output coordinate system
        iaxis_march : (3,) float narray; default=None -> coord_out.i
            the normalized x-axis that defines the direction to march
        icoord : int; default=None -> coord_out+1
            the starting index for the coordinate systems that will be created
            and placed in new_coords; useful for debugging
        nodes_tol : float; default=None -> dstation
            the tolerance bending the plane to pull nodes from

        Returns
        -------
        force_sum / moment_sum : (nstations, 3) float ndarray
            the forces/moments at the station
        new_coords: dict[int, CORD2R]
            the station march coords starting from icoord
        nelems, nnodes: (nstations,) int ndarray
            the number of elements/nodes included in the summation

        Notes
        -----
        1.  Clip elements based on centroid.
            Elements that are less than the ith station are kept.
        2.  Get the nodes on the opposite side of the clip plane
            and the ones within nodes_tol of the plane.
        3.  Extract the interface loads and sum them about the
            summation point.

        Examples
        --------
        Imagine a swept aircraft wing.  Define a coordinate system
        in the primary direction of the sweep.  Note that station 0
        doesn't have to be perfectly at the root of the wing.

        Create stations from this point.

        .. todo:: Not Tested...Does 3b work?  Can 3a give the right answer?

        """
        _check_array(nids, 'int', 1)
        _check_array(nid_cd, 'int', 2)
        _check_array(eids, 'int', 1)
        _check_array(xyz_cid0, 'float', 2)
        _check_array(element_centroids_cid0, 'float', 2)
        _check_array(stations, 'float', 1)
        assert len(nids) == nid_cd.shape[0]
        assert len(nids) == xyz_cid0.shape[0]
        assert len(eids) == element_centroids_cid0.shape[0]
        assert isinstance(coords, dict), coords
        assert isinstance(icd_transform, dict), icd_transform
        if icoord is None:
            icoord = max(coords) + 1

        nid_cd = _get_nid_cd_from_nid_cp_cd(nid_cd)
        idtype = nid_cd.dtype
        fdtype = xyz_cid0.dtype

        if iaxis_march is None:
            iaxis_march = deepcopy(coord_out.i)

        eids = np.asarray(eids, dtype=idtype)
        nids = np.asarray(nids, dtype=idtype)
        assert len(eids.shape) == 1, eids.shape
        assert len(nids.shape) == 1, nids.shape
        assert len(stations.shape) == 1, stations.shape
        nstations = len(stations)
        assert coord_out.type in ['CORD2R', 'CORD1R'], coord_out.type
        #assert coord_march.type in ['CORD2R', 'CORD1R'], coord_march.type
        #i_axis_march = deepcopy(coord_march.i)
        #del iaxis_march
        assert np.array_equal(nids, np.unique(nids))
        assert np.array_equal(eids, np.unique(eids))
        # ----------------------------------------------------------------------

        element_centroids_coord = coord_out.transform_node_to_local_array(element_centroids_cid0)
        xyz_coord = coord_out.transform_node_to_local_array(xyz_cid0)

        # we can hardcode this cause we transformed it into the output frame
        # which is similar to the output coordinate frame
        # minus sign because we march into the page for axial
        #
        x_elem_centroid = element_centroids_coord[:, idir]
        x_coord = xyz_coord[:, idir]

        eids = np.unique(eids)
        force_sum = np.full((nstations, 3), np.nan, dtype=fdtype)
        moment_sum = np.full((nstations, 3), np.nan, dtype=fdtype)

        offset = np.zeros(3, dtype=fdtype)
        new_coords = {}
        nelems_list = []
        nnodes_list = []

        # calculate the value of dstation_x for the 1st station
        doffset = (stations[1] - stations[0]) * iaxis_march
        dsummation_point = coord_out.origin + doffset
        dsummation_point_coord = coord_out.transform_node_to_local(dsummation_point)
        dstation = dsummation_point_coord[idir]
        if nodes_tol is None:
            nodes_tol = dstation
        #assert nodes_tol >= 0., nodes_tol

        for istation, station in enumerate(stations):
            # summation point creation
            # in the global coordinate system
            offset = station * iaxis_march

            summation_point = coord_out.origin + offset  # in basic frame
            summation_point_coord = coord_out.transform_node_to_local(summation_point)
            station_x = summation_point_coord[idir]
            # we're picking the elements on one side of the centroid
            # and nodes on the other side

            # Calculate the nodes on the boundary.
            #
            # If we make a cutting plane and find all the nodes on
            # one side of the cutting plane, we can take all the
            # nodes within some tolerance of the station direction and
            # find the free nodes
            #
            ielem = np.where(x_elem_centroid <= station_x)[0]
            nelem = len(ielem)
            nelems_list.append(nelem)

            # We want to do this similarly for the elements.  Ideally,
            # we want a single line of elements (and not include extra)
            # to reduce the size of the array sooner.
            #
            # We'll include a few extra nodes (with dstation) to include
            # all the nodes for the boundary elements.  The extra nodes
            # will have +forces and -forces, so they'll cancel.  The alternative
            # is that we miss some the boundary nodes and get no force/moment.
            #
            jnode = np.where(x_coord >= (station_x-nodes_tol))[0]
            nnode = len(jnode)
            nnodes_list.append(nnode)

            coord_save = deepcopy(coord_out)
            coord_save.cid = icoord
            coord_save.translate(offset)

            assert np.allclose(coord_save.e1, summation_point)
            new_coords[icoord] = deepcopy(coord_save)
            #station_x_old = station_x

            # we'd break if we knew the user was traveling in the
            # "correct" direction, but we don't
            icoord += 1
            if nelem == 0:
                continue
            if nnode == 0:
                continue

            eidsi = eids[ielem]
            nidsj = nids[jnode]
            assert len(eidsi) == len(np.unique(eidsi))
            assert len(nidsj) == len(np.unique(nidsj))
            if 0: # pragma: no cover
                # I don't think this will work...
                forcei, momenti = self.extract_freebody_loads(
                    eidsi,
                    coord_out, coords, nid_cd, icd_transform,
                    # xyz_cid0, summation_point,
                    itime=itime, debug=debug, log=log)

                force_sum[istation, :] = forcei.sum(axis=0)
                # TODO: extract_freebody_loads doesn't sum forces/moments
                #       sum loads about summation point
                moment_sum[istation, :] = momenti.sum(axis=0)
            else:
                #print(f'eids={eidsi}; nids={nidsj}')
                force_sumi, moment_sumi = self.extract_interface_loads(
                    nidsj, eidsi,
                    coord_out, coords, nid_cd, icd_transform,
                    xyz_cid0, summation_point, assume_sorted=True,
                    itime=itime, stop_on_nan=stop_on_nan,
                    debug=debug, log=log)
                #log.info(f'neids={len(ielem):d} nnodes={len(jnode):d} station={station:g}; '
                         #f'force={force_sumi} moment={moment_sumi}')
                if not np.isfinite(force_sumi[0]):
                    if stop_on_nan:
                        raise RuntimeError(f'station={station} is nan; force_sum={force_sumi}')
                    continue
                force_sum[istation, :] = force_sumi
                moment_sum[istation, :] = moment_sumi
        nelems = np.array(nelems_list, dtype=idtype)
        nnodes = np.array(nnodes_list, dtype=idtype)
        return force_sum, moment_sum, new_coords, nelems, nnodes

    def add_sort1(self, dt, node_id, eid, ename, t1, t2, t3, r1, r2, r3):
        """unvectorized method for adding SORT1 transient data"""
        assert self.sort_method == 1, self
        assert eid is not None, eid
        #print(self.code_information())
        #assert isinstance(eid, integer_types) and eid > 0, 'dt=%s eid=%s' % (dt, eid)
        assert isinstance(node_id, int), node_id
        self._times[self.itime] = dt

        if self.is_unique:
            self.node_element[self.itime, self.itotal, :] = [node_id, eid]
            self.element_names[self.itime, self.itotal] = ename
        else:
            self.node_element[self.itotal, :] = [node_id, eid]
            self.element_names[self.itotal] = ename

        #self.node_element[self.itotal, :] = [eid, node_id]
        #self.element_names[self.itotal] = ename
        self.data[self.itime, self.itotal, :] = [t1, t2, t3, r1, r2, r3]
        self.itotal += 1

    def get_stats(self, short: bool=False) -> list[str]:
        if not self.is_built:
            return [
                f'<{self.__class__.__name__}>; table_name={self.table_name!r}\n',
                f'  ntimes: {self.ntimes:d}\n',
                f'  ntotal: {self.ntotal:d}\n',
                f'  _ntotals: {self._ntotals}\n',
            ]

        #nelements = self.nelements
        ntimes = self.ntimes
        #nnodes = self.nnodes
        ntotal = self.ntotal
        nelements = np.unique(self.node_element[:, 1]).size

        msg = []
        if self.nonlinear_factor not in (None, np.nan):  # transient
            msgi = '  type=%s ntimes=%i nelements=%i ntotal=%i\n' % (
                self.__class__.__name__, ntimes, nelements, ntotal)
            if self.ntotal != min(self._ntotals):
                msgi += f'  _ntotals={self._ntotals}\n'
            ntimes_word = 'ntimes'
        else:
            msgi = '  type=%s nelements=%i total=%i\n' % (
                self.__class__.__name__, nelements, ntotal)
            if self.ntotal != min(self._ntotals):
                msgi += f'  _ntotals={self._ntotals}\n'
            ntimes_word = '1'
        msg.append(msgi)
        headers = self.get_headers()
        n = len(headers)

        #element_names = [name.strip() for name in unique(self.element_names)]
        msg.append('  data: [%s, ntotal, %i] where %i=[%s]\n' % (
            ntimes_word, n, n, ', '.join(headers)))
        msg.append('  data.shape=%s\n' % str(self.data.shape))
        msg.append(f'  element type: {self.element_name}\n')
        msg += self.get_data_code()
        return msg

    #def get_element_index(self, eids):
        #itot = searchsorted(eids, self.node_element[:, 0])
        #return itot

    #def eid_to_element_node_index(self, eids):
        #ind = ravel([searchsorted(self.node_element[:, 0] == eid) for eid in eids])
        #return ind

    def write_csv(self, csv_file: TextIO,
                  is_exponent_format: bool=False,
                  is_mag_phase: bool=False, is_sort1: bool=True,
                  write_header: bool=True):
        """
        Grid Point Forces Table
        -----------------------
        Flag, NID, SubcaseID, iTime, EID, TYPE,    Fx,      Fy,      Fz,      Mx,      My,      Mz
        13,   101,         1,     1, 0,   APPLIED, 30.9864, 19.7278, 70.2515, 53.3872, 80.9687, 77.4302
        13,   101,         1,     1, 301, RBE3,    41.9012, 53.6651, 0.09483, 76.041,  67.506,  98.0225
        13,   101,         1,     1, 0,   SPC,     71.6306, 97.0527, 89.8733, 5.89262, 61.0523, 48.9043
        13,   101,         1,     1, 301, CTRIA3,  84.5273, 69.36,   92.3295, 52.7074, 77.9904, 68.905
        13,   101,         1,     1, 302, CQUAD4,  97.7843, 11.7545, 99.3901, 44.9476, 70.818,  7.47876

        """
        if write_header:
            name = str(self.__class__.__name__)
            csv_file.write('# %s\n' % name)
            headers = ['Flag', 'NID', 'SubcaseID', 'iTime', 'Nid', 'Eid', 'EName', 'T1', 'T2', 'T3', 'R1', 'R2', 'R3']
            csv_file.write('# ' + ','.join(headers) + '\n')

        flag = 13
        isubcase = self.isubcase

        assert self.is_unique, self.is_unique
        # sort1 as sort1

        nids = self.node_element[:, :, 0]
        eids = self.node_element[:, :, 1]
        nid_len = '%d' % len(str(nids.max()))
        eid_len = '%d' % len(str(eids.max()))

        for itime in range(self.ntimes):
            #dt = self._times[itime]
            t1 = self.data[itime, :, 0]
            t2 = self.data[itime, :, 1]
            t3 = self.data[itime, :, 2]
            r1 = self.data[itime, :, 3]
            r2 = self.data[itime, :, 4]
            r3 = self.data[itime, :, 5]

            nids = self.node_element[itime, :, 0]
            eids = self.node_element[itime, :, 1]
            enames = self.element_names[itime, :]
            #ntotal = self._ntotals[itime]

            for (nid, eid, ename, t1i, t2i, t3i, r1i, r2i, r3i) in zip(
                nids, eids, enames, t1, t2, t3, r1, r2, r3):

                if is_exponent_format:
                    vals2 = write_floats_13e_long([t1i, t2i, t3i, r1i, r2i, r3i])
                    [t1i, t2i, t3i, r1i, r2i, r3i] = vals2

                csv_file.write(f'{flag}, {isubcase}, {itime}, {nid:{nid_len}d}, {eid:{eid_len}d}, {ename.strip():>8s}, '
                               f'{t1i}, {t2i}, {t3i}, {r1i}, {r2i}, {r3i}\n')
        return

    def write_f06(self, f06_file: TextIO, header=None, page_stamp='PAGE %s',
                  page_num: int=1, is_mag_phase: bool=False, is_sort1: bool=True):
        if header is None:
            header = []
        msg = self._get_f06_msg()

        ntimes = self.data.shape[0]
        if self.is_unique:
            for itime in range(ntimes):
                dt = self._times[itime]
                header = _eigenvalue_header(self, header, itime, ntimes, dt)
                f06_file.write(''.join(header + msg))

                #print("self.data.shape=%s itime=%s ieids=%s" % (
                    #str(self.data.shape), itime, str(ieids)))
                #[t1, t2, t3, r1, r2, r3]
                t1 = self.data[itime, :, 0]
                t2 = self.data[itime, :, 1]
                t3 = self.data[itime, :, 2]
                r1 = self.data[itime, :, 3]
                r2 = self.data[itime, :, 4]
                r3 = self.data[itime, :, 5]

                nids = self.node_element[itime, :, 0]
                eids = self.node_element[itime, :, 1]
                enames = self.element_names[itime, :]

                zero = ' '
                ntotal = self._ntotals[itime]
                #print(self._ntotals)
                assert len(eids) == len(nids)
                assert len(enames) == len(nids), 'enames=%s nnids=%s' % (len(enames), len(nids))
                assert len(t1) == len(nids)
                assert len(t2) == len(nids)
                assert len(t3) == len(nids)
                assert len(r1) == len(nids)
                assert len(r2) == len(nids)
                assert len(nids) == ntotal, 'len(nids)=%s ntotal=%s' % (len(nids), ntotal)

                for (i, nid, eid, ename, t1i, t2i, t3i, r1i, r2i, r3i) in zip(
                     range(ntotal), nids, eids, enames, t1, t2, t3, r1, r2, r3):

                    #print(nid, eid, ename, t1i)
                    vals = [t1i, t2i, t3i, r1i, r2i, r3i]
                    vals2 = write_floats_13e(vals)
                    [f1, f2, f3, m1, m2, m3] = vals2
                    if eid == 0:
                        f06_file.write('   %8s    %10s    %s      %-13s  %-13s  %-13s  %-13s  %-13s  %s\n' % (
                            nid, eid, ename, f1, f2, f3, m1, m2, m3))
                        zero = '0'
                    else:
                        f06_file.write('%s  %8s    %10s    %s      %-13s  %-13s  %-13s  %-13s  %-13s  %s\n' % (
                            zero, nid, eid, ename, f1, f2, f3, m1, m2, m3))
                        zero = ' '
                f06_file.write(page_stamp % page_num)
                page_num += 1
        else:
            nids = self.node_element[:, 0]
            eids = self.node_element[:, 1]
            enames = self.element_names

            for itime in range(ntimes):
                dt = self._times[itime]
                header = _eigenvalue_header(self, header, itime, ntimes, dt)
                f06_file.write(''.join(header + msg))

                #print("self.data.shape=%s itime=%s ieids=%s" % (
                    #str(self.data.shape), itime, str(ieids)))

                #[t1, t2, t3, r1, r2, r3]
                t1 = self.data[itime, :, 0]
                t2 = self.data[itime, :, 1]
                t3 = self.data[itime, :, 2]
                r1 = self.data[itime, :, 3]
                r2 = self.data[itime, :, 4]
                r3 = self.data[itime, :, 5]

                zero = ' '
                for (nid, eid, ename, t1i, t2i, t3i, r1i, r2i, r3i) in zip(
                     nids, eids, enames, t1, t2, t3, r1, r2, r3):

                    vals = [t1i, t2i, t3i, r1i, r2i, r3i]
                    vals2 = write_floats_13e(vals)
                    [f1, f2, f3, m1, m2, m3] = vals2
                    if eid == 0:
                        f06_file.write('   %8s    %10s    %s      %-13s  %-13s  %-13s  %-13s  %-13s  %s\n' % (
                            nid, eid, ename, f1, f2, f3, m1, m2, m3))
                        zero = '0'
                    else:
                        f06_file.write('%s  %8s    %10s    %s      %-13s  %-13s  %-13s  %-13s  %-13s  %s\n' % (
                            zero, nid, eid, ename, f1, f2, f3, m1, m2, m3))
                        zero = ' '
                f06_file.write(page_stamp % page_num)
                page_num += 1
        return page_num - 1

    def write_f06_time(self, f06_file, itime=0, i=None, header=None, page_num=1, page_stamp=''):
        if header is None:
            header = []
        dt = self._times[itime]
        msg = self._get_f06_msg()

        ntimes = self.data.shape[0]
        header = _eigenvalue_header(self, header, itime, ntimes, dt)
        f06_file.write(''.join(header + msg))

        #print("self.data.shape=%s itime=%s ieids=%s" % (
            #str(self.data.shape), itime, str(ieids)))
        #[t1, t2, t3, r1, r2, r3]
        t1 = self.data[itime, i, 0]
        t2 = self.data[itime, i, 1]
        t3 = self.data[itime, i, 2]
        r1 = self.data[itime, i, 3]
        r2 = self.data[itime, i, 4]
        r3 = self.data[itime, i, 5]

        nids = self.node_element[itime, i, 0]
        eids = self.node_element[itime, i, 1]
        enames = self.element_names[itime, i]

        zero = ' '
        ntotal = self._ntotals[itime]
        #print(self._ntotals)
        assert len(eids) == len(nids)
        assert len(enames) == len(nids), 'enames=%s nnids=%s' % (len(enames), len(nids))
        assert len(t1) == len(nids)
        assert len(t2) == len(nids)
        assert len(t3) == len(nids)
        assert len(r1) == len(nids)
        assert len(r2) == len(nids)
        assert len(nids) <= ntotal, 'len(nids)=%s ntotal=%s' % (len(nids), ntotal)

        for (i, nid, eid, ename, t1i, t2i, t3i, r1i, r2i, r3i) in zip(
             range(ntotal), nids, eids, enames, t1, t2, t3, r1, r2, r3):

            #print(nid, eid, ename, t1i)
            vals = [t1i, t2i, t3i, r1i, r2i, r3i]
            vals2 = write_floats_13e(vals)
            [f1, f2, f3, m1, m2, m3] = vals2
            if eid == 0:
                f06_file.write('   %8s    %10s    %s      %-13s  %-13s  %-13s  %-13s  %-13s  %s\n' % (
                    nid, eid, ename, f1, f2, f3, m1, m2, m3))
                zero = '0'
            else:
                f06_file.write('%s  %8s    %10s    %s      %-13s  %-13s  %-13s  %-13s  %-13s  %s\n' % (
                    zero, nid, eid, ename, f1, f2, f3, m1, m2, m3))
                zero = ' '
        #f.write(page_stamp % page_num)
        page_num += 1

    def _get_f06_msg(self):
        msg = [
            '                                          G R I D   P O I N T   F O R C E   B A L A N C E\n',
            ' \n',
            '   POINT-ID    ELEMENT-ID     SOURCE             T1             T2             T3             R1             R2             R3\n',
           #'0     13683          3736    TRIAX6         4.996584E+00   0.0            1.203093E+02   0.0            0.0            0.0'
           #'      13683          3737    TRIAX6        -4.996584E+00   0.0           -1.203093E+02   0.0            0.0            0.0'
           #'      13683                  *TOTALS*       6.366463E-12   0.0           -1.364242E-12   0.0            0.0            0.0'
        ]
        return msg

    def get_headers(self) -> list[str]:
        headers = ['f1', 'f2', 'f3', 'm1', 'm2', 'm3']
        return headers

    def write_op2(self, op2_file, op2_ascii, itable, new_result,
                  date, is_mag_phase=False, endian='>'):
        """writes an OP2"""
        frame = inspect.currentframe()
        call_frame = inspect.getouterframes(frame, 2)
        op2_ascii.write(f'{self.__class__.__name__}.write_op2: {call_frame[1][3]}\n')

        if itable == -1:
            self._write_table_header(op2_file, op2_ascii, date, subtable_name_default=b'OGPFB1  ',
                                     include_date=False)
            itable = -3

        #if isinstance(self.nonlinear_factor, float):
            #op2_format = '%sif' % (7 * self.ntimes)
            #raise NotImplementedError()
        #else:
            #op2_format = 'i21f'
        #s = Struct(op2_format)

        # table 4 info
        #ntimes = self.data.shape[0]
        #nnodes = self.data.shape[1]

        # 21 = 1 node, 3 principal, 6 components, 9 vectors, 2 p/ovm
        #ntotal = ((nnodes * 21) + 1) + (nelements * 4)

        ntotali = self.num_wide

        #print('shape = %s' % str(self.data.shape))
        #assert self.ntimes == 1, self.ntimes

        #device_code = self.device_code
        op2_ascii.write(f'  ntimes = {self.ntimes}\n')

        #fmt = '%2i %6f'
        #print('ntotal=%s' % (ntotal))
        #assert ntotal == 193, ntotal

        if self.is_sort1:
            struct1 = Struct(endian + b'2i 8s 6f')
        else:
            raise NotImplementedError('SORT2')

        for itime in range(self.ntimes):
            self._write_table_3(op2_file, op2_ascii, new_result, itable, itime)

            # record 4
            #print('stress itable = %s' % itable)
            itable -= 1

            t1 = self.data[itime, :, 0]
            t2 = self.data[itime, :, 1]
            t3 = self.data[itime, :, 2]
            r1 = self.data[itime, :, 3]
            r2 = self.data[itime, :, 4]
            r3 = self.data[itime, :, 5]

            nids = self.node_element[itime, :, 0]
            eids = self.node_element[itime, :, 1]
            enames = self.element_names[itime, :]

            nids_device = nids * 10 + self.device_code
            assert nids.min() > 0, nids.min()
            nnodes = len(nids)

            ntotal = ntotali * nnodes
            header = [4, itable, 4,
                      4, 1, 4,
                      4, 0, 4,
                      4, ntotal, 4,
                      4 * ntotal]
            op2_file.write(pack('%ii' % len(header), *header))
            op2_ascii.write('r4 [4, 0, 4]\n')
            op2_ascii.write(f'r4 [4, {itable:d}, 4]\n')
            op2_ascii.write(f'r4 [4, {4 * ntotal:d}, 4]\n')

            #zero = ' '
            ntotali = self._ntotals[itime]
            #print(self._ntotals)
            assert len(eids) == len(nids)
            assert len(enames) == len(nids), 'enames=%s nnids=%s' % (len(enames), len(nids))
            assert len(t1) == len(nids)
            assert len(t2) == len(nids)
            assert len(t3) == len(nids)
            assert len(r1) == len(nids)
            assert len(r2) == len(nids)
            assert len(nids) <= ntotali, 'len(nids)=%s ntotali=%s' % (len(nids), ntotali)

            for (i, nid, eid, ename, t1i, t2i, t3i, r1i, r2i, r3i) in zip(
                 range(ntotali), nids_device, eids, enames, t1, t2, t3, r1, r2, r3):

                #print(nid, eid, ename, t1i)
                data = [nid, eid, ename.encode('ascii'), t1i, t2i, t3i, r1i, r2i, r3i]
                #print('  nid=%s eid=%s data=%s' % (nid, eid, str(data[2:])))
                op2_ascii.write('  nid=%-3s eid=%-3s data=%s\n' % (nid, eid, str(data[2:])))
                op2_file.write(struct1.pack(*data))

            itable -= 1
            header = [4 * ntotal,]
            op2_file.write(pack('i', *header))
            op2_ascii.write('footer = %s\n' % header)
            new_result = False
        return itable


class ComplexGridPointForcesArray(GridPointForces):
    def __init__(self, data_code, is_sort1, isubcase, dt):
        GridPointForces.__init__(self, data_code, is_sort1, isubcase)
        #self.code = [self.format_code, self.sort_code, self.s_code]

        # do the element_names/node_element vectors change with the time step
        self.is_unique = False

        #self.ielement = 0
        #self.nelements = 0  # result specific
        #self.nnodes = None

    @property
    def is_real(self) -> bool:
        return False

    @property
    def is_complex(self) -> bool:
        return True

    def _reset_indices(self) -> None:
        self.itotal = 0
        #self.ielement = 0

    @property
    def element_name(self):
        headers = [name.strip() for name in np.unique(self.element_names)]
        #headers = np.unique(self.element_names)
        return str(', '.join(headers))

    def build(self):
        """sizes the vectorized attributes of the ComplexGridPointForcesArray"""
        #print('ntimes=%s nelements=%s ntotal=%s' % (
            #self.ntimes, self.nelements, self.ntotal))
        #self.ntotal += 5  # TODO: remove
        #print("self.ntotal = %s" % self.ntotal)
        #print("self.itotal = %s" % self.itotal)
        #print("self._ntotals = %s" % self._ntotals)

        #self.is_unique = False
        if self.ntotal != max(self._ntotals) or 1:
            self.ntotal = max(self._ntotals)
            self.is_unique = True

        assert self.ntimes > 0, 'ntimes=%s' % self.ntimes
        #assert self.nelements > 0, 'nelements=%s' % self.nelements
        assert self.ntotal > 0, 'ntotal=%s' % self.ntotal
        #self.names = []

        #self.nnodes = nnodes_per_element
        #self.nelements //= nnodes_per_element
        self.itime = 0
        #self.ielement = 0
        self.itotal = 0
        #self.ntimes = 0
        #self.nelements = 0

        #print("***name=%s type=%s nnodes_per_element=%s ntimes=%s nelements=%s ntotal=%s" % (
            #self.element_names, self.element_type, nnodes_per_element,
            #self.ntimes, self.nelements, self.ntotal))
        dtype, idtype, fdtype = get_times_dtype(self.nonlinear_factor, self.size, self.analysis_fmt)

        self._times = np.zeros(self.ntimes, dtype=self.analysis_fmt)

        if self.is_unique:
            self.node_element = np.zeros((self.ntimes, self.ntotal, 2), dtype=idtype)
            self.element_names = np.empty((self.ntimes, self.ntotal), dtype='U8')
        else:
            self.node_element = np.zeros((self.ntotal, 2), dtype=idtype)
            self.element_names = np.empty(self.ntotal, dtype='U8')
        #[t1, t2, t3, r1, r2, r3]
        self.data = np.zeros((self.ntimes, self.ntotal, 6), dtype='complex64')

    def build_dataframe(self):
        """
        major-axis - the axis

        mode              1     2   3
        freq              1.0   2.0 3.0
        nodeID ElementID Item
        1      2         T1
                         T2
                         ...

        major_axis / top = [
            [1, 2, 3],
            [1.0, 2.0, 3.0]
        ]
        minor_axis / headers = [T1, T2, T3, R1, R2, R3]
        name = mode
        """
        import pandas as pd
        headers = self.get_headers()
        #name = self.name
        if self.is_unique:
            ntimes = self.data.shape[0]
            nnodes = self.data.shape[1]
            #nvalues = ntimes * nnodes
            node_element = self.node_element.reshape((ntimes * nnodes, 2))
            df1 = pd.DataFrame(node_element)
            df1.columns = ['NodeID', 'ElementID']
            df2 = pd.DataFrame(self.element_names[0, :])
            df2.columns = ['ElementType']
            df3 = pd.DataFrame(self.data[0])
            df3.columns = headers
            data_frame = df1.join([df2, df3])
            #print(data_frame)
        else:
            node_element = [self.node_element[:, 0], self.node_element[:, 1]]
            if self.nonlinear_factor not in (None, np.nan):
                column_names, column_values = self._build_dataframe_transient_header()
                data_frame = pd.Panel(
                    self.data, items=column_values,
                    major_axis=node_element, minor_axis=headers).to_frame()
                data_frame.columns.names = column_names
                data_frame.index.names = ['NodeID', 'ElementID', 'Item']
            else:
                data_frame = pd.Panel(
                    self.data,
                    major_axis=node_element, minor_axis=headers).to_frame()
                data_frame.columns.names = ['Static']
                data_frame.index.names = ['NodeID', 'ElementID', 'Item']
            #print(self.data_frame)
        self.data_frame = data_frame

    def _build_dataframe(self):
        """::
        major-axis - the axis

        mode              1     2   3
        freq              1.0   2.0 3.0
        nodeID ElementID Item
        1      2         T1
                         T2
                         ...

        major_axis / top = [
            [1, 2, 3],
            [1.0, 2.0, 3.0]
        ]
        minor_axis / headers = [T1, T2, T3, R1, R2, R3]
        name = mode
        """
        headers = self.get_headers()
        import pandas as pd
        column_names, column_values = self._build_dataframe_transient_header()
        if self.is_unique:
            #node_element = [self.node_element[:, 0], self.node_element[:, 1]]
            ntimes = self.data.shape[0]
            nnodes = self.data.shape[1]
            node_element_temp = self.node_element.reshape((ntimes * nnodes, 2))
            node_element = [node_element_temp[:, 0], node_element_temp[:, 1]]
            self.data_frame = pd.Panel(
                self.data, items=column_values,
                major_axis=node_element, minor_axis=headers).to_frame()
            self.data_frame.columns.names = column_names
            self.data_frame.index.names = ['NodeID', 'ElementID', 'Item']
        else:
            node_element = [self.node_element[:, 0], self.node_element[:, 1]]
            #print('column_names =', column_names)
            #for name, values in zip(column_names, column_values):
                #print('  %s = %s' % (name, values))
            self.data_frame = pd.Panel(
                self.data, items=column_values,
                major_axis=node_element, minor_axis=headers).to_frame()
            self.data_frame.columns.names = column_names
            self.data_frame.index.names = ['NodeID', 'ElementID', 'Item']

    def __eq__(self, table):  # pragma: no cover
        self._eq_header(table)
        assert self.is_sort1 == table.is_sort1
        if not np.array_equal(self.node_element, table.node_element) and np.array_equal(self.element_names, table.element_names):
            assert self.node_element.shape == table.node_element.shape, 'node_element shape=%s table.shape=%s' % (self.node_element.shape, table.node_element.shape)
            assert self.element_names.shape == table.element_names.shape, 'element_names shape=%s table.shape=%s' % (self.element_names.shape, table.element_names.shape)

            msg = 'table_name=%r class_name=%s\n' % (self.table_name, self.__class__.__name__)
            msg += '%s\n' % str(self.code_information())
            msg += '(Eid, Nid, EName)\n'
            for (nid1, eid1), ename1, (nid2, eid2), ename2 in zip(self.node_element, self.element_names,
                                                                  table.element_names, table.element_names):
                msg += '(%s, %s, %s)    (%s, %s, %s)\n' % (nid1, eid1, ename1, nid2, eid2, ename2)
            print(msg)
            raise ValueError(msg)
        if self.is_unique:
            if not np.array_equal(self.data, table.data):
                msg = 'table_name=%r class_name=%s\n' % (self.table_name, self.__class__.__name__)
                msg += '%s\n' % str(self.code_information())
                i = 0
                for itime in range(self.ntimes):
                    #print('is_unique =', self.is_unique)
                    sys.stdout.flush()
                    for ie, e in enumerate(self.node_element[itime, :, :]):
                        print(e)
                        (eid, nid) = e
                        ename1 = self.element_names[itime, ie]
                        ename2 = self.element_names[itime, ie]
                        t1 = self.data[itime, ie, :]
                        t2 = table.data[itime, ie, :]
                        (t11, t21, t31, r11, r21, r31) = t1
                        (t12, t22, t32, r12, r22, r32) = t2

                        if not np.array_equal(t1, t2):
                            msg += '(%s, %s, %s)    (%s, %s, %s, %s, %s, %s)  (%s, %s, %s, %s, %s, %s)\n' % (
                                eid, nid, ename1,
                                t11, t21, t31, r11, r21, r31,
                                t12, t22, t32, r12, r22, r32)
                            i += 1
                            if i > 10:
                                print(msg)
                                raise ValueError(msg)
                    #print(msg)
                    if i > 0:
                        raise ValueError(msg)
        else:
            if not np.array_equal(self.data, table.data):
                msg = 'table_name=%r class_name=%s\n' % (self.table_name, self.__class__.__name__)
                msg += '%s\n' % str(self.code_information())
                i = 0
                for itime in range(self.ntimes):
                    for ie, e in enumerate(self.node_element):
                        (eid, nid) = e
                        ename1 = self.element_names[ie]
                        ename2 = self.element_names[ie]
                        t1 = self.data[itime, ie, :]
                        t2 = table.data[itime, ie, :]
                        (t11, t21, t31, r11, r21, r31) = t1
                        (t12, t22, t32, r12, r22, r32) = t2

                        if not np.array_equal(t1, t2):
                            msg += '(%s, %s, %s)    (%s, %s, %s, %s, %s, %s)  (%s, %s, %s, %s, %s, %s)\n' % (
                                eid, nid, ename1,
                                t12, t22, t32, r12, r22, r32,
                                t12, t22, t32, r12, r22, r32)
                            i += 1
                            if i > 10:
                                print(msg)
                                raise ValueError(msg)
                    #print(msg)
                    if i > 0:
                        raise ValueError(msg)
        return True

    def add_sort1(self, dt, node_id, eid, ename, t1, t2, t3, r1, r2, r3):
        """unvectorized method for adding SORT1 transient data"""
        assert self.sort_method == 1, self
        assert eid is not None, eid
        #assert isinstance(eid, integer_types) and eid > 0, 'dt=%s eid=%s' % (dt, eid)
        assert isinstance(node_id, int), node_id

        self._times[self.itime] = dt
        if self.is_unique:
            self.node_element[self.itime, self.itotal, :] = [node_id, eid]
            self.element_names[self.itime, self.itotal] = ename
        else:
            self.node_element[self.itotal, :] = [node_id, eid]
            self.element_names[self.itotal] = ename
        set_3d_data(self.data, self.itime, self.itotal, [t1, t2, t3, r1, r2, r3], self.size)
        self.itotal += 1

    def get_stats(self, short: bool=False) -> list[str]:
        if not self.is_built:
            return [
                f'<{self.__class__.__name__}>; table_name={self.table_name!r}\n',
                f'  ntimes: {self.ntimes:d}\n',
                f'  ntotal: {self.ntotal:d}\n',
            ]

        #nelements = self.nelements
        ntimes = self.ntimes
        #nnodes = self.nnodes
        ntotal = self.ntotal
        nelements = np.unique(self.node_element[:, 1]).size

        msg = []
        if self.nonlinear_factor not in (None, np.nan):  # transient
            msgi = '  type=%s ntimes=%i nelements=%i ntotal=%i\n' % (
                self.__class__.__name__, ntimes, nelements, ntotal)
            ntimes_word = 'ntimes'
        else:
            msgi = '  type=%s nelements=%i total=%i\n' % (
                self.__class__.__name__, nelements, ntotal)
            ntimes_word = '1'
        msg.append(msgi)
        headers = self.get_headers()
        n = len(headers)

        #element_names = [name.strip() for name in unique(self.element_names)]
        msg.append('  data: [%s, ntotal, %i] where %i=[%s]\n' % (
            ntimes_word, n, n, ', '.join(headers)))
        msg.append(f'  node_element.shape={self.node_element.shape}\n')
        msg.append(f'  element_names.shape={self.element_names.shape}\n')
        msg.append(f'  data.shape={self.data.shape}\n')
        msg.append(f'  element type: {self.element_name}\n')
        msg += self.get_data_code()
        return msg

    #def get_element_index(self, eids):
        #itot = searchsorted(eids, self.node_element[:, 0])
        #return itot

    #def eid_to_element_node_index(self, eids):
        #ind = ravel([searchsorted(self.node_element[:, 0] == eid) for eid in eids])
        #return ind

    def write_f06(self, f06_file, header=None, page_stamp='PAGE %s',
                  page_num: int=1, is_mag_phase: bool=False, is_sort1: bool=True):
        if header is None:
            header = []
        msg = self._get_f06_msg(is_mag_phase=is_mag_phase, is_sort1=is_sort1)

        ntimes = self.data.shape[0]
        if self.is_unique:
            for itime in range(ntimes):
                dt = self._times[itime]
                header = _eigenvalue_header(self, header, itime, ntimes, dt)
                f06_file.write(''.join(header + msg))

                #print("self.data.shape=%s itime=%s ieids=%s" % (
                    #str(self.data.shape), itime, str(ieids)))

                #[t1, t2, t3, r1, r2, r3]
                t1 = self.data[itime, :, 0]
                t2 = self.data[itime, :, 1]
                t3 = self.data[itime, :, 2]
                r1 = self.data[itime, :, 3]
                r2 = self.data[itime, :, 4]
                r3 = self.data[itime, :, 5]

                eids = self.node_element[itime, :, 1]
                nids = self.node_element[itime, :, 0]
                enames = self.element_names[itime, :]

                zero = ' '
                ntotal = self._ntotals[itime]
                for (unused_i, nid, eid, ename, t1i, t2i, t3i, r1i, r2i, r3i) in zip(
                         range(ntotal), nids, eids, enames, t1, t2, t3, r1, r2, r3):

                    vals = [t1i, t2i, t3i, r1i, r2i, r3i]
                    vals2 = write_imag_floats_13e(vals, is_mag_phase)
                    [f1r, f2r, f3r, m1r, m2r, m3r, f1i, f2i, f3i, m1i, m2i, m3i] = vals2
                    if eid == 0:
                        f06_file.write(
                            '   %8s    %10s    %8s      %-13s  %-13s  %-13s  %-13s  %-13s  %s\n'
                            '   %8s    %10s    %8s      %-13s  %-13s  %-13s  %-13s  %-13s  %s\n' % (
                                nid, eid, ename, f1r, f2r, f3r, m1r, m2r, m3r,
                                '', '', '', f1i, f2i, f3i, m1i, m2i, m3i,
                        ))
                        zero = '0'
                    else:
                        f06_file.write(
                            '%s  %8s    %10s    %8s      %-13s  %-13s  %-13s  %-13s  %-13s  %s\n'
                            '   %8s    %10s    %8s      %-13s  %-13s  %-13s  %-13s  %-13s  %s\n' % (
                                zero, nid, eid, ename, f1r, f2r, f3r, m1r, m2r, m3r,
                                '', '', '', f1i, f2i, f3i, m1i, m2i, m3i,))
                        zero = ' '
                f06_file.write(page_stamp % page_num)
                page_num += 1
        else:
            eids = self.node_element[:, 1]
            nids = self.node_element[:, 0]
            enames = self.element_names
            for itime in range(ntimes):
                dt = self._times[itime]
                header = _eigenvalue_header(self, header, itime, ntimes, dt)
                f06_file.write(''.join(header + msg))

                #print("self.data.shape=%s itime=%s ieids=%s" % (
                    #str(self.data.shape), itime, str(ieids)))

                #[t1, t2, t3, r1, r2, r3]
                t1 = self.data[itime, :, 0]
                t2 = self.data[itime, :, 1]
                t3 = self.data[itime, :, 2]
                r1 = self.data[itime, :, 3]
                r2 = self.data[itime, :, 4]
                r3 = self.data[itime, :, 5]

                zero = ' '
                for (nid, eid, ename, t1i, t2i, t3i, r1i, r2i, r3i) in zip(
                     nids, eids, enames, t1, t2, t3, r1, r2, r3):

                    vals = [t1i, t2i, t3i, r1i, r2i, r3i]
                    vals2 = write_imag_floats_13e(vals, is_mag_phase)
                    [f1r, f2r, f3r, m1r, m2r, m3r, f1i, f2i, f3i, m1i, m2i, m3i] = vals2
                    if eid == 0:
                        f06_file.write(
                            '   %8s    %10s    %8s      %-13s  %-13s  %-13s  %-13s  %-13s  %s\n'
                            '   %8s    %10s    %8s      %-13s  %-13s  %-13s  %-13s  %-13s  %s\n' % (
                                nid, eid, ename, f1r, f2r, f3r, m1r, m2r, m3r,
                                '', '', '', f1i, f2i, f3i, m1i, m2i, m3i,
                        ))
                        zero = '0'
                    else:
                        f06_file.write(
                            '%s  %8s    %10s    %8s      %-13s  %-13s  %-13s  %-13s  %-13s  %s\n'
                            '   %8s    %10s    %8s      %-13s  %-13s  %-13s  %-13s  %-13s  %s\n' % (
                                zero, nid, eid, ename, f1r, f2r, f3r, m1r, m2r, m3r,
                                '', '', '', f1i, f2i, f3i, m1i, m2i, m3i,))
                        zero = ' '
                f06_file.write(page_stamp % page_num)
                page_num += 1
                eids = self.node_element[:, 1]
                nids = self.node_element[:, 0]
                enames = self.element_names

                for itime in range(ntimes):
                    dt = self._times[itime]
                    header = _eigenvalue_header(self, header, itime, ntimes, dt)
                    f06_file.write(''.join(header + msg))

                    #print("self.data.shape=%s itime=%s ieids=%s" % (
                        #str(self.data.shape), itime, str(ieids)))

                    #[t1, t2, t3, r1, r2, r3]
                    t1 = self.data[itime, :, 0]
                    t2 = self.data[itime, :, 1]
                    t3 = self.data[itime, :, 2]
                    r1 = self.data[itime, :, 3]
                    r2 = self.data[itime, :, 4]
                    r3 = self.data[itime, :, 5]

                    zero = ' '
                    for (nid, eid, ename, t1i, t2i, t3i, r1i, r2i, r3i) in zip(
                         nids, eids, enames, t1, t2, t3, r1, r2, r3):

                        vals = [t1i, t2i, t3i, r1i, r2i, r3i]
                        vals2 = write_imag_floats_13e(vals, is_mag_phase)
                        [f1r, f2r, f3r, m1r, m2r, m3r, f1i, f2i, f3i, m1i, m2i, m3i] = vals2
                        if eid == 0:
                            f06_file.write(
                                '   %8s    %10s    %8s      %-13s  %-13s  %-13s  %-13s  %-13s  %s\n'
                                '   %8s    %10s    %8s      %-13s  %-13s  %-13s  %-13s  %-13s  %s\n' % (
                                    nid, eid, ename, f1r, f2r, f3r, m1r, m2r, m3r,
                                    '', '', '', f1i, f2i, f3i, m1i, m2i, m3i,
                            ))
                            zero = '0'
                        else:
                            f06_file.write(
                                '%s  %8s    %10s    %8s      %-13s  %-13s  %-13s  %-13s  %-13s  %s\n'
                                '   %8s    %10s    %8s      %-13s  %-13s  %-13s  %-13s  %-13s  %s\n' % (
                                    zero, nid, eid, ename, f1r, f2r, f3r, m1r, m2r, m3r,
                                    '', '', '', f1i, f2i, f3i, m1i, m2i, m3i,))
                            zero = ' '
                    f06_file.write(page_stamp % page_num)
                    page_num += 1
        return page_num - 1

    def _get_f06_msg(self, is_mag_phase=True, is_sort1=True):
        msg = [
            '                                  C O M P L E X   G R I D   P O I N T   F O R C E   B A L A N C E\n',
            #sort,
            #' \n'
            #'   POINT-ID    ELEMENT-ID     SOURCE             T1             T2             T3             R1             R2             R3\n'
            #'         19            21    TRIA3          0.0            0.0            0.0            0.0            0.0            0.0'
            #'                                            0.0            0.0            0.0            0.0            0.0            0.0'
        ]
        if is_mag_phase:
            msg += ['                                                          (REAL/IMAGINARY)\n \n']
            raise NotImplementedError('mag/phase')
        else:
            msg += ['                                                          (REAL/IMAGINARY)\n \n']

        if is_sort1:
            #msg += ['   FREQ        ELEMENT-ID     SOURCE             T1             T2             T3             R1             R2             R3\n']
            msg += ['   POINT-ID    ELEMENT-ID     SOURCE             T1             T2             T3             R1             R2             R3\n']
        else:
            # TODO: get this right
            msg += ['   POINT-ID    ELEMENT-ID     SOURCE             T1             T2             T3             R1             R2             R3\n']

        return msg

    def get_headers(self) -> list[str]:
        headers = ['f1', 'f2', 'f3', 'm1', 'm2', 'm3']
        return headers

    def write_op2(self, op2_file, op2_ascii, itable, new_result,
                  date, is_mag_phase=False, endian='>'):
        """writes an OP2"""
        import inspect
        from struct import Struct, pack
        frame = inspect.currentframe()
        call_frame = inspect.getouterframes(frame, 2)
        op2_ascii.write(f'{self.__class__.__name__}.write_op2: {call_frame[1][3]}\n')

        if itable == -1:
            self._write_table_header(op2_file, op2_ascii, date)
            itable = -3

        # table 4 info
        #ntimes = self.data.shape[0]
        #nnodes = self.data.shape[1]

        # 21 = 1 node, 3 principal, 6 components, 9 vectors, 2 p/ovm
        #ntotal = ((nnodes * 21) + 1) + (nelements * 4)

        ntotali = self.num_wide

        #print('shape = %s' % str(self.data.shape))

        #device_code = self.device_code
        op2_ascii.write(f'  ntimes = {self.ntimes}\n')

        #fmt = '%2i %6f'
        #print('ntotal=%s' % (ntotal))

        if self.is_sort1:
            struct1 = Struct(endian + b'2i 8s 12f')
        else:
            raise NotImplementedError('SORT2')

        for itime in range(self.ntimes):
            self._write_table_3(op2_file, op2_ascii, new_result, itable, itime)

            # record 4
            itable -= 1

            nids_all = self.node_element[itime, :, 0]
            inids = np.where(nids_all > 0)[0]
            nids = nids_all[inids]
            eids = self.node_element[itime, inids, 1]
            enames = self.element_names[itime, inids]

            t1 = self.data[itime, inids, 0]
            t2 = self.data[itime, inids, 1]
            t3 = self.data[itime, inids, 2]
            r1 = self.data[itime, inids, 3]
            r2 = self.data[itime, inids, 4]
            r3 = self.data[itime, inids, 5]

            nids_device = nids * 10 + self.device_code
            assert nids.min() > 0, nids.min()
            nnodes = len(nids)

            ntotal = ntotali * nnodes
            header = [4, itable, 4,
                      4, 1, 4,
                      4, 0, 4,
                      4, ntotal, 4,
                      4 * ntotal]
            op2_file.write(pack('%ii' % len(header), *header))
            op2_ascii.write('r4 [4, 0, 4]\n')
            op2_ascii.write(f'r4 [4, {itable:d}, 4]\n')
            op2_ascii.write(f'r4 [4, {4 * ntotal:d}, 4]\n')

            #zero = ' '
            ntotal = self._ntotals[itime]
            #print(self._ntotals)
            assert len(eids) == len(nids)
            assert len(enames) == len(nids), 'enames=%s nnids=%s' % (len(enames), len(nids))
            assert len(t1) == len(nids)
            assert len(t2) == len(nids)
            assert len(t3) == len(nids)
            assert len(r1) == len(nids)
            assert len(r2) == len(nids)
            assert len(nids) <= ntotal, 'len(nids)=%s ntotal=%s' % (len(nids), ntotal)

            for (i, nid, eid, ename, t1i, t2i, t3i, r1i, r2i, r3i) in zip(
                 range(ntotal), nids_device, eids, enames, t1, t2, t3, r1, r2, r3):

                #print(nid, eid, ename, t1i)
                data = [nid, eid, ename.encode('ascii'),
                        t1i.real, t2i.real, t3i.real, r1i.real, r2i.real, r3i.real,
                        t1i.imag, t2i.imag, t3i.imag, r1i.imag, r2i.imag, r3i.imag]
                #print('  nid=%s eid=%s data=%s' % (nid, eid, str(data[2:])))
                op2_ascii.write('  nid=%-3s eid=%-3s data=%s\n' % (nid, eid, str(data[2:])))
                op2_file.write(struct1.pack(*data))
                assert len(data) + 1 == self.num_wide

            itable -= 1
            header = [4 * ntotal,]
            op2_file.write(pack('i', *header))
            op2_ascii.write('footer = %s\n' % header)
            new_result = False
        return itable

def _get_nid_cd_from_nid_cp_cd(nid_cd: np.ndarray) -> np.ndarray:
    if nid_cd.shape[1] == 3:
        assert isinstance(nid_cd[0, :].tolist()[0], int), nid_cd
        nid_cp_cd = nid_cd
        nid_cd = nid_cp_cd[:, [0, 2]]
    assert nid_cd.shape[1] == 2, nid_cd
    return nid_cd

def _check_array(x, dtype, dim: int, msg=''):
    if isinstance(x, np.ndarray):
        assert len(x.shape) == dim, x.shape
        if dtype == 'float':
            assert x.dtype.name in ['float32', 'float64'], x.dtype.name
        elif dtype == 'int':
            assert x.dtype.name in ['int32', 'int64'], x.dtype.name
        else:  # pragma: no cover
            raise NotImplementedError(dtype)
    else:
        assert isinstance(x, (list, tuple)), x
        if dim == 1:
            val = x[0]
        elif dim == 2:
            val = x[0][0]
        else:
            raise NotImplementedError(dim)
        if dtype == 'float':
            assert isinstance(val, (float, np.float32, np.float64)), x
        elif dtype == 'int':
            assert isinstance(val, (int, np.int32, np.int64)), x
        else:  # pragma: no cover
            raise NotImplementedError(dtype)
    return

def set_3d_data(data: np.ndarray, itime: int, itotal: int, values: list[float], size: int):
    """annoying way to handle underflow"""
    if size == 4:
        #MIN_FLOAT32 = np.finfo(np.float32).min
        datai = np.array(values)
        try:
            datai = datai.astype(data.dtype)
            data[itime, itotal, :] = datai
        except FloatingPointError:
            if data.dtype.name in {'float32', 'float64'}:
                for i, dataii in enumerate(datai):
                    try:
                        data[itime, itotal, i] = dataii
                    except FloatingPointError:
                        data[itime, itotal, i] = 0.
            else:
                for i, dataii in enumerate(datai):
                    try:
                        data.real[itime, itotal, i] = dataii.real
                    except FloatingPointError:
                        data.real[itime, itotal, i] = 0.

                    try:
                        data.imag[itime, itotal, i] = dataii.imag
                    except FloatingPointError:
                        data.imag[itime, itotal, i] = 0.
    else:
        data[itime, itotal, :] = values
