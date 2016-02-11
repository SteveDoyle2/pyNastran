from __future__ import print_function
from six import iteritems
import numpy as np
from numpy import zeros, unique, array_equal, empty
from pyNastran.op2.result_objects.op2_objects import ScalarObject
from pyNastran.f06.f06_formatting import write_floats_13e, _eigenvalue_header, write_imag_floats_13e
from pyNastran.op2.vector_utils import transform_force_from_local_to_global, transform_force_from_global_to_local

try:
    import pandas as pd
except ImportError:
    pass


class RealGridPointForcesArray(ScalarObject):
    def __init__(self, data_code, is_sort1, isubcase, dt):
        ScalarObject.__init__(self, data_code, isubcase, apply_data_code=True)
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

    def is_real(self):
        return True

    def is_complex(self):
        return False

    def _reset_indices(self):
        self.itotal = 0
        #self.ielement = 0

    @property
    def element_name(self):
        headers = [name.strip() for name in unique(self.element_names) if name.strip()]
        #headers = unique(self.element_names)
        return str(', '.join(headers))

    def build(self):
        #print("self.ielement = %s" % self.ielement)
        #print('ntimes=%s nelements=%s ntotal=%s' % (self.ntimes, self.nelements, self.ntotal))
        if self.is_built:
            return

        assert self.ntimes > 0, 'ntimes=%s' % self.ntimes
        #assert self.nelements > 0, 'nelements=%s' % self.nelements
        assert self.ntotal > 0, 'ntotal=%s' % self.ntotal
        #if self.ntotal != max(self._ntotals) or self.ntotal != min(self._ntotals):
            #raise ValueError('RealGridPointForcesArray: ntotal=%s _ntotals=%s' % (self.ntotal, self._ntotals))

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
        self.is_built = True

        #print("***name=%s ntimes=%s ntotal=%s" % (
            #self.element_names, self.ntimes, self.ntotal))
        dtype = 'float32'
        if isinstance(self.nonlinear_factor, int):
            dtype = 'int32'

        self._times = zeros(self.ntimes, dtype=dtype)

        if self.is_unique:
            #raise NotImplementedError('not unique')
            self.node_element = zeros((self.ntimes, self.ntotal, 2), dtype='int32')
            self.element_names = empty((self.ntimes, self.ntotal), dtype='U8')
        else:
            self.node_element = zeros((self.ntotal, 2), dtype='int32')
            self.element_names = empty(self.ntotal, dtype='U8')

        #[t1, t2, t3, r1, r2, r3]
        self.data = zeros((self.ntimes, self.ntotal, 6), dtype='float32')

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
        headers = self.get_headers()
        #name = self.name
        if self.is_unique:
            ntimes = self.data.shape[0]
            nnodes = self.data.shape[1]
            nvalues = ntimes * nnodes
            node_element = self.node_element.reshape((ntimes * nnodes, 2))
            if self.nonlinear_factor is not None:
                column_names, column_values = self._build_dataframe_transient_header()

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

                df2 = pd.DataFrame(node_element)
                df2.columns = ['NodeID', 'ElementID']
                df3 = pd.DataFrame(self.element_names.ravel())
                df3.columns = ['ElementType']

                dfs = [df2, df3]
                for i, header in enumerate(headers):
                    df = pd.DataFrame(self.data[:, :, i].ravel())
                    df.columns = [header]
                    dfs.append(df)
                self.data_frame = df1.join(dfs)
                #print(self.data_frame)
            else:
                df1 = pd.DataFrame(node_element)
                df1.columns = ['NodeID', 'ElementID']
                df2 = pd.DataFrame(self.element_names[0, :])
                df2.columns = ['ElementType']
                df3 = pd.DataFrame(self.data[0])
                df3.columns = headers
                self.data_frame = df1.join([df2, df3])
                #print(self.data_frame)
        else:
            node_element = [self.node_element[:, 0], self.node_element[:, 1]]
            if self.nonlinear_factor is not None:
                column_names, column_values = self._build_dataframe_transient_header()
                self.data_frame = pd.Panel(self.data, items=column_values, major_axis=node_element, minor_axis=headers).to_frame()
                self.data_frame.columns.names = column_names
                self.data_frame.index.names = ['NodeID', 'ElementID', 'Item']
            else:
                self.data_frame = pd.Panel(self.data, major_axis=node_element, minor_axis=headers).to_frame()
                self.data_frame.columns.names = ['Static']
                self.data_frame.index.names = ['NodeID', 'ElementID', 'Item']
            #print(self.data_frame)

    def __eq__(self, table):
        assert self.is_sort1() == table.is_sort1()
        assert self.nonlinear_factor == table.nonlinear_factor
        assert self.ntotal == table.ntotal
        assert self.table_name == table.table_name, 'table_name=%r table.table_name=%r' % (self.table_name, table.table_name)
        assert self.approach_code == table.approach_code
        if self.nonlinear_factor is not None:
            assert np.array_equal(self._times, table._times), 'class_name=%s times=%s table.times=%s' % (
                self.class_name, self._times, table._times)
        if not np.array_equal(self.node_element, table.node_element) and array_equal(self.element_names, table.element_names):
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

        if self.is_unique:
            if not np.array_equal(self.data, table.data):
                msg = 'table_name=%r class_name=%s\n' % (self.table_name, self.__class__.__name__)
                msg += '%s\n' % str(self.code_information())
                i = 0
                for itime in range(self.ntimes):
                    for ie, e in enumerate(self.node_element):
                        (nid, eid) = e
                        ename1 = self.element_names[itime, ie]
                        ename2 = self.element_names[itime, ie]
                        t1 = self.data[itime, ie, :]
                        t2 = table.data[itime, ie, :]
                        (t11, t21, t31, r11, r21, r31) = t1
                        (t12, t22, t32, r12, r22, r32) = t2

                        if not np.array_equal(t1, t2):
                            msg += '(%s, %s, %s)    (%s, %s, %s, %s, %s, %s)  (%s, %s, %s, %s, %s, %s)\n' % (
                                nid, eid, ename1,
                                t12, t22, t32, r12, r22, r32,
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
                        (nid, eid) = e
                        ename1 = self.element_names[ie]
                        ename2 = self.element_names[ie]
                        t1 = self.data[itime, ie, :]
                        t2 = table.data[itime, ie, :]
                        (t11, t21, t31, r11, r21, r31) = t1
                        (t12, t22, t32, r12, r22, r32) = t2

                        if not np.array_equal(t1, t2):
                            msg += '(%s, %s, %s)    (%s, %s, %s, %s, %s, %s)  (%s, %s, %s, %s, %s, %s)\n' % (
                                nid, eid, ename1,
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

    def extract_gpforcei(self, panel_nids, panel_eids, coord_in, coord_out,
                         nid_cd, i_transform, beta_transforms, itime=0):
        """
        Parameters
        ----------
        panel_nids : (Nn, ) int ndarray
            all the nodes to consider
        panel_eids : (Ne, ) int ndarray
            all the elements to consider
        coord_in : CORD2R()
            the input coordinate system
        coord_out : CORD2R()
            the output coordinate system
        nid_cd : (M, 2) int ndarray
            the (BDF.point_ids, cd) array
        i_transform : dict[cd] = (Mi, ) intndarray
            the mapping for nid_cd
        beta_transforms : dict[cd] = (3, 3) float ndarray
            the mapping for nid_cd
        itime : int; default=0
            the time to extract loads for

        .. warning :: the function signature will change...
        .. todo :: sum of moments about a point must have an rxF term to get the
                   same value as Patran.
        .. todo :: doesn't support transient/frequency/modal based results
        """
        panel_eids = asarray(panel_eids)
        panel_nids = asarray(panel_nids)

        # todo handle multiple values for itime
        gpforce_nids = self.node_element[itime, :, 0]
        gpforce_eids = self.node_element[itime, :, 1]
        # TODO: remove 0s in gpforce_nids/gpforce_eids to handle transient results
        #       be careful of the sum row

        assert isinstance(panel_eids[0], int), type(panel_eids[0])
        assert isinstance(panel_nids[0], int), type(panel_nids[0])
        is_in = in1d(gpforce_nids, panel_nids, assume_unique=False)
        is_in2 = in1d(gpforce_eids[is_in], panel_eids, assume_unique=False)
        irange = arange(len(gpforce_nids), dtype='int32')[is_in][is_in2]
        force_local = self.data[itime, irange, :]

        # TODO: smash force_local from nelements per node to 1 line per node
        force_global = transform_force_from_local_to_global(force_local, gpforce_nids, nid_cd,
                                                            i_transform, beta_transforms)
        force_in_global = force_global[:, :3]
        total_force_global = force_in_global.sum(axis=0)

        total_force_local = transform_force_from_global_to_local(total_force_global, coord_in, coord_out)

        return total_force_global, total_force_local

    def add(self, dt, node_id, eid, ename, t1, t2, t3, r1, r2, r3):
        assert isinstance(node_id, int), node_id
        self.add_sort1(dt, eid, node_id, ename, t1, t2, t3, r1, r2, r3)

    def add_sort1(self, dt, node_id, eid, ename, t1, t2, t3, r1, r2, r3):
        assert eid is not None, eid
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

    def get_stats(self):
        if not self.is_built:
            return [
                '<%s>\n' % self.__class__.__name__,
                '  ntimes: %i\n' % self.ntimes,
                '  ntotal: %i\n' % self.ntotal,
                '  _ntotals: %s\n' % self._ntotals,
            ]

        #nelements = self.nelements
        ntimes = self.ntimes
        #nnodes = self.nnodes
        ntotal = self.ntotal
        nelements = unique(self.node_element[:, 1]).size

        msg = []
        if self.nonlinear_factor is not None:  # transient
            msgi = '  type=%s ntimes=%i nelements=%i ntotal=%i\n' % (
                self.__class__.__name__, ntimes, nelements, ntotal)
            if self.ntotal != min(self._ntotals):
                msgi += '  _ntotals=%s\n' % (self._ntotals)
            ntimes_word = 'ntimes'
        else:
            msgi = '  type=%s nelements=%i total=%i\n' % (
                self.__class__.__name__, nelements, ntotal)
            if self.ntotal != min(self._ntotals):
                msgi += '  _ntotals=%s\n' % (self._ntotals)
            ntimes_word = 1
        msg.append(msgi)
        headers = self.get_headers()
        n = len(headers)

        #element_names = [name.strip() for name in unique(self.element_names)]
        msg.append('  data: [%s, ntotal, %i] where %i=[%s]\n' % (ntimes_word, n, n,
                                                                 ', '.join(headers)))
        msg.append('  data.shape=%s\n' % str(self.data.shape))
        msg.append('  element type: %s\n  ' % self.element_name)
        msg += self.get_data_code()
        return msg

    #def get_element_index(self, eids):
        #itot = searchsorted(eids, self.node_element[:, 0])
        #return itot

    #def eid_to_element_node_index(self, eids):
        #ind = ravel([searchsorted(self.node_element[:, 0] == eid) for eid in eids])
        #return ind

    def write_f06(self, f, header=None, page_stamp='PAGE %s', page_num=1, is_mag_phase=False, is_sort1=True):
        if header is None:
            header = []
        msg = self._get_f06_msg()

        ntimes = self.data.shape[0]
        if self.is_unique:
            #print('RealGridPointForcesArray.write_f06 with is_unique=True is not supported')

            for itime in range(ntimes):
                dt = self._times[itime]
                header = _eigenvalue_header(self, header, itime, ntimes, dt)
                f.write(''.join(header + msg))

                #print("self.data.shape=%s itime=%s ieids=%s" % (str(self.data.shape), itime, str(ieids)))

                #[t1, t2, t3, r1, r2, r3]
                t1 = self.data[itime, :, 0]
                t2 = self.data[itime, :, 1]
                t3 = self.data[itime, :, 2]
                r1 = self.data[itime, :, 3]
                r2 = self.data[itime, :, 4]
                r3 = self.data[itime, :, 5]

                nids = self.node_element[itime, 0]
                eids = self.node_element[itime, 1]
                enames = self.element_names[itime, :]

                zero = ' '
                ntotal = self._ntotals[itime]
                for (i, nid, eid, ename, t1i, t2i, t3i, r1i, r2i, r3i) in zip(
                     range(ntotal), nids, eids, enames, t1, t2, t3, r1, r2, r3):

                    vals = [t1i, t2i, t3i, r1i, r2i, r3i]
                    vals2 = write_floats_13e(vals)
                    [f1, f2, f3, m1, m2, m3] = vals2
                    if eid == 0:
                        f.write('   %8s    %10s    %s      %-13s  %-13s  %-13s  %-13s  %-13s  %s\n' % (
                                nid, eid, ename, f1, f2, f3, m1, m2, m3))
                        zero = '0'
                    else:
                        f.write('%s  %8s    %10s    %s      %-13s  %-13s  %-13s  %-13s  %-13s  %s\n' % (
                                zero, nid, eid, ename, f1, f2, f3, m1, m2, m3))
                        zero = ' '
                f.write(page_stamp % page_num)
                page_num += 1
        else:
            nids = self.node_element[:, 0]
            eids = self.node_element[:, 1]
            enames = self.element_names

            for itime in range(ntimes):
                dt = self._times[itime]
                header = _eigenvalue_header(self, header, itime, ntimes, dt)
                f.write(''.join(header + msg))

                #print("self.data.shape=%s itime=%s ieids=%s" % (str(self.data.shape), itime, str(ieids)))

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
                        f.write('   %8s    %10s    %s      %-13s  %-13s  %-13s  %-13s  %-13s  %s\n' % (
                                nid, eid, ename, f1, f2, f3, m1, m2, m3))
                        zero = '0'
                    else:
                        f.write('%s  %8s    %10s    %s      %-13s  %-13s  %-13s  %-13s  %-13s  %s\n' % (
                                zero, nid, eid, ename, f1, f2, f3, m1, m2, m3))
                        zero = ' '
                f.write(page_stamp % page_num)
                page_num += 1
        return page_num - 1

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

    def get_headers(self):
        headers = ['f1', 'f2', 'f3', 'm1', 'm2', 'm3']
        return headers


class ComplexGridPointForcesArray(ScalarObject):
    def __init__(self, data_code, is_sort1, isubcase, dt):
        ScalarObject.__init__(self, data_code, isubcase, apply_data_code=True)
        #self.code = [self.format_code, self.sort_code, self.s_code]

        #self.ntimes = 0  # or frequency/mode
        self.ntotal = 0
        self.itotal = 0

        # do the element_names/node_element vectors change with the time step
        self.is_unique = False

        #self.ielement = 0
        #self.nelements = 0  # result specific
        #self.nnodes = None

    def is_real(self):
        return False

    def is_complex(self):
        return True

    def _reset_indices(self):
        self.itotal = 0
        #self.ielement = 0

    @property
    def element_name(self):
        headers = [name.strip() for name in unique(self.element_names)]
        #headers = unique(self.element_names)
        return str(', '.join(headers))

    def build(self):
        #print('ntimes=%s nelements=%s ntotal=%s' % (self.ntimes, self.nelements, self.ntotal))
        if self.is_built:
            return
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
        self.is_built = True

        #print("***name=%s type=%s nnodes_per_element=%s ntimes=%s nelements=%s ntotal=%s" % (
            #self.element_names, self.element_type, nnodes_per_element, self.ntimes, self.nelements, self.ntotal))
        dtype = 'float32'
        if isinstance(self.nonlinear_factor, int):
            dtype = 'int32'

        self._times = zeros(self.ntimes, dtype=dtype)

        if self.is_unique:
            self.node_element = zeros((self.ntimes, self.ntotal, 2), dtype='int32')
            self.element_names = empty((self.ntimes, self.ntotal), dtype='U8')
        else:
            self.node_element = zeros((self.ntotal, 2), dtype='int32')
            self.element_names = empty(self.ntotal, dtype='U8')
        #[t1, t2, t3, r1, r2, r3]
        self.data = zeros((self.ntimes, self.ntotal, 6), dtype='complex64')

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
        headers = self.get_headers()
        #name = self.name
        if self.is_unique:
            ntimes = self.data.shape[0]
            nnodes = self.data.shape[1]
            nvalues = ntimes * nnodes
            node_element = self.node_element.reshape((ntimes * nnodes, 2))
            df1 = pd.DataFrame(node_element)
            df1.columns = ['NodeID', 'ElementID']
            df2 = pd.DataFrame(self.element_names[0, :])
            df2.columns = ['ElementType']
            df3 = pd.DataFrame(self.data[0])
            df3.columns = headers
            self.data_frame = df1.join([df2, df3])
            #print(self.data_frame)
        else:
            node_element = [self.node_element[:, 0], self.node_element[:, 1]]
            if self.nonlinear_factor is not None:
                column_names, column_values = self._build_dataframe_transient_header()
                self.data_frame = pd.Panel(self.data, items=column_values, major_axis=node_element, minor_axis=headers).to_frame()
                self.data_frame.columns.names = column_names
                self.data_frame.index.names = ['NodeID', 'ElementID', 'Item']
            else:
                self.data_frame = pd.Panel(self.data, major_axis=node_element, minor_axis=headers).to_frame()
                self.data_frame.columns.names = ['Static']
                self.data_frame.index.names = ['NodeID', 'ElementID', 'Item']
            #print(self.data_frame)

    def _build_dataframe(self):
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
        headers = self.get_headers()
        column_names, column_values = self._build_dataframe_transient_header()
        if self.is_unique:
            #node_element = [self.node_element[:, 0], self.node_element[:, 1]]
            ntimes = self.data.shape[0]
            nnodes = self.data.shape[1]
            node_element_temp = self.node_element.reshape((ntimes * nnodes, 2))
            node_element = [node_element_temp[:, 0], node_element_temp[:, 1]]
            self.data_frame = pd.Panel(self.data, items=column_values, major_axis=node_element, minor_axis=headers).to_frame()
            self.data_frame.columns.names = column_names
            self.data_frame.index.names = ['NodeID', 'ElementID', 'Item']
        else:
            node_element = [self.node_element[:, 0], self.node_element[:, 1]]
            #print('column_names =', column_names)
            #for name, values in zip(column_names, column_values):
                #print('  %s = %s' % (name, values))
            self.data_frame = pd.Panel(self.data, items=column_values, major_axis=node_element, minor_axis=headers).to_frame()
            self.data_frame.columns.names = column_names
            self.data_frame.index.names = ['NodeID', 'ElementID', 'Item']

    def __eq__(self, table):
        assert self.is_sort1() == table.is_sort1()
        assert self.nonlinear_factor == table.nonlinear_factor
        assert self.ntotal == table.ntotal
        assert self.table_name == table.table_name, 'table_name=%r table.table_name=%r' % (self.table_name, table.table_name)
        assert self.approach_code == table.approach_code
        if self.nonlinear_factor is not None:
            assert np.array_equal(self._times, table._times), 'class_name=%s times=%s table.times=%s' % (
                self.class_name, self._times, table._times)
        if not np.array_equal(self.node_element, table.node_element) and array_equal(self.element_names, table.element_names):
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
                    for ie, e in enumerate(self.node_element):
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
                                t12, t22, t32, r12, r22, r32,
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
        assert eid is not None, eid
        assert isinstance(node_id, int), node_id

        if self.is_unique:
            self.node_element[self.itime, self.itotal, :] = [node_id, eid]
            self.element_names[self.itime, self.itotal] = ename
        else:
            self.node_element[self.itotal, :] = [node_id, eid]
            self.element_names[self.itotal] = ename
        self.data[self.itime, self.itotal, :] = [t1, t2, t3, r1, r2, r3]
        self.itotal += 1

    def get_stats(self):
        if not self.is_built:
            return [
                '<%s>\n' % self.__class__.__name__,
                '  ntimes: %i\n' % self.ntimes,
                '  ntotal: %i\n' % self.ntotal,
            ]

        #nelements = self.nelements
        ntimes = self.ntimes
        #nnodes = self.nnodes
        ntotal = self.ntotal
        nelements = unique(self.node_element[:, 1]).size

        msg = []
        if self.nonlinear_factor is not None:  # transient
            msgi = '  type=%s ntimes=%i nelements=%i ntotal=%i\n' % (
                self.__class__.__name__, ntimes, nelements, ntotal)
            ntimes_word = 'ntimes'
        else:
            msgi = '  type=%s nelements=%i total=%i\n' % (
                self.__class__.__name__, nelements, ntotal)
            ntimes_word = 1
        msg.append(msgi)
        headers = self.get_headers()
        n = len(headers)

        #element_names = [name.strip() for name in unique(self.element_names)]
        msg.append('  data: [%s, ntotal, %i] where %i=[%s]\n' % (ntimes_word, n, n,
                                                                 ', '.join(headers)))
        msg.append('  data.shape=%s\n' % str(self.data.shape))
        msg.append('  element type: %s\n  ' % self.element_name)
        msg += self.get_data_code()
        return msg

    #def get_element_index(self, eids):
        #itot = searchsorted(eids, self.node_element[:, 0])
        #return itot

    #def eid_to_element_node_index(self, eids):
        #ind = ravel([searchsorted(self.node_element[:, 0] == eid) for eid in eids])
        #return ind

    def write_f06(self, f, header=None, page_stamp='PAGE %s', page_num=1, is_mag_phase=False, is_sort1=True):
        if header is None:
            header = []
        msg = self._get_f06_msg(is_mag_phase=is_mag_phase, is_sort1=is_sort1)

        ntimes = self.data.shape[0]

        if self.is_unique:
            for itime in range(ntimes):
                dt = self._times[itime]
                header = _eigenvalue_header(self, header, itime, ntimes, dt)
                f.write(''.join(header + msg))

                #print("self.data.shape=%s itime=%s ieids=%s" % (str(self.data.shape), itime, str(ieids)))

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
                for (i, nid, eid, ename, t1i, t2i, t3i, r1i, r2i, r3i) in zip(
                     range(ntotal), nids, eids, enames, t1, t2, t3, r1, r2, r3):

                    vals = [t1i, t2i, t3i, r1i, r2i, r3i]
                    vals2 = write_imag_floats_13e(vals, is_mag_phase)
                    [f1r, f2r, f3r, m1r, m2r, m3r, f1i, f2i, f3i, m1i, m2i, m3i] = vals2
                    if eid == 0:
                        f.write('   %8s    %10s    %8s      %-13s  %-13s  %-13s  %-13s  %-13s  %s\n'
                                '   %8s    %10s    %8s      %-13s  %-13s  %-13s  %-13s  %-13s  %s\n' % (
                                nid, eid, ename, f1r, f2r, f3r, m1r, m2r, m3r,
                                '', '', '', f1i, f2i, f3i, m1i, m2i, m3i,
                        ))
                        zero = '0'
                    else:
                        f.write('%s  %8s    %10s    %8s      %-13s  %-13s  %-13s  %-13s  %-13s  %s\n'
                                '   %8s    %10s    %8s      %-13s  %-13s  %-13s  %-13s  %-13s  %s\n' % (
                                zero, nid, eid, ename, f1r, f2r, f3r, m1r, m2r, m3r,
                                '', '', '', f1i, f2i, f3i, m1i, m2i, m3i,))
                        zero = ' '
                f.write(page_stamp % page_num)
                page_num += 1
        else:
            eids = self.node_element[:, 1]
            nids = self.node_element[:, 0]
            enames = self.element_names
            for itime in range(ntimes):
                dt = self._times[itime]
                header = _eigenvalue_header(self, header, itime, ntimes, dt)
                f.write(''.join(header + msg))

                #print("self.data.shape=%s itime=%s ieids=%s" % (str(self.data.shape), itime, str(ieids)))

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
                        f.write('   %8s    %10s    %8s      %-13s  %-13s  %-13s  %-13s  %-13s  %s\n'
                                '   %8s    %10s    %8s      %-13s  %-13s  %-13s  %-13s  %-13s  %s\n' % (
                                nid, eid, ename, f1r, f2r, f3r, m1r, m2r, m3r,
                                '', '', '', f1i, f2i, f3i, m1i, m2i, m3i,
                        ))
                        zero = '0'
                    else:
                        f.write('%s  %8s    %10s    %8s      %-13s  %-13s  %-13s  %-13s  %-13s  %s\n'
                                '   %8s    %10s    %8s      %-13s  %-13s  %-13s  %-13s  %-13s  %s\n' % (
                                zero, nid, eid, ename, f1r, f2r, f3r, m1r, m2r, m3r,
                                '', '', '', f1i, f2i, f3i, m1i, m2i, m3i,))
                        zero = ' '
                f.write(page_stamp % page_num)
                page_num += 1
                eids = self.node_element[:, 1]
                nids = self.node_element[:, 0]
                enames = self.element_names

                for itime in range(ntimes):
                    dt = self._times[itime]
                    header = _eigenvalue_header(self, header, itime, ntimes, dt)
                    f.write(''.join(header + msg))

                    #print("self.data.shape=%s itime=%s ieids=%s" % (str(self.data.shape), itime, str(ieids)))

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
                            f.write('   %8s    %10s    %8s      %-13s  %-13s  %-13s  %-13s  %-13s  %s\n'
                                    '   %8s    %10s    %8s      %-13s  %-13s  %-13s  %-13s  %-13s  %s\n' % (
                                    nid, eid, ename, f1r, f2r, f3r, m1r, m2r, m3r,
                                    '', '', '', f1i, f2i, f3i, m1i, m2i, m3i,
                            ))
                            zero = '0'
                        else:
                            f.write('%s  %8s    %10s    %8s      %-13s  %-13s  %-13s  %-13s  %-13s  %s\n'
                                    '   %8s    %10s    %8s      %-13s  %-13s  %-13s  %-13s  %-13s  %s\n' % (
                                    zero, nid, eid, ename, f1r, f2r, f3r, m1r, m2r, m3r,
                                    '', '', '', f1i, f2i, f3i, m1i, m2i, m3i,))
                            zero = ' '
                    f.write(page_stamp % page_num)
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
            mag_phase
        else:
            msg += ['                                                          (REAL/IMAGINARY)\n \n']

        if is_sort1:
            #msg += ['   FREQ        ELEMENT-ID     SOURCE             T1             T2             T3             R1             R2             R3\n']
            msg += ['   POINT-ID    ELEMENT-ID     SOURCE             T1             T2             T3             R1             R2             R3\n']
        else:
            # TODO: get this right
            msg += ['   POINT-ID    ELEMENT-ID     SOURCE             T1             T2             T3             R1             R2             R3\n']

        return msg

    def get_headers(self):
        headers = ['f1', 'f2', 'f3', 'm1', 'm2', 'm3']
        return headers
