from six import iteritems
from numpy import array, zeros, unique, array_equal, empty
from struct import pack
from pyNastran.op2.resultObjects.op2_Objects import ScalarObject
from pyNastran.f06.f06_formatting import writeFloats13E, _eigenvalue_header, writeImagFloats13E



class RealGridPointForcesArray(ScalarObject):
    def __init__(self, data_code, is_sort1, isubcase, dt):
        ScalarObject.__init__(self, data_code, isubcase, apply_data_code=True)
        #self.code = [self.format_code, self.sort_code, self.s_code]

        #self.ntimes = 0  # or frequency/mode
        self.ntotal = 0

        # do the element_names/node_element vectors change with the time step
        self.is_unique = False

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
        headers = [name.strip() for name in unique(self.element_names)]
        #headers = unique(self.element_names)
        return str(b', '.join(headers))

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
        if self.ntotal != min(self._ntotals):
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

        #print("***name=%s type=%s nnodes_per_element=%s ntimes=%s nelements=%s ntotal=%s" % (
            #self.element_names, self.element_type, nnodes_per_element, self.ntimes, self.nelements, self.ntotal))
        dtype = 'float32'
        if isinstance(self.nonlinear_factor, int):
            dtype = 'int32'

        self._times = zeros(self.ntimes, dtype=dtype)

        if self.is_unique:
            self.node_element = zeros((self.ntimes, self.ntotal, 2), dtype='int32')
            self.element_names = empty(self.ntimes, self.ntotal, dtype='S8')
        else:
            self.node_element = zeros((self.ntotal, 2), dtype='int32')
            self.element_names = empty(self.ntotal, dtype='S8')

        #[t1, t2, t3, r1, r2, r3]
        self.data = zeros((self.ntimes, self.ntotal, 6), dtype='float32')

    def __eq__(self, table):
        assert self.is_sort1() == table.is_sort1()
        assert self.nonlinear_factor == table.nonlinear_factor
        assert self.ntotal == table.ntotal
        assert self.table_name == table.table_name, 'table_name=%r table.table_name=%r' % (self.table_name, table.table_name)
        assert self.approach_code == table.approach_code
        if not array_equal(self.node_element, table.node_element) and array_equal(self.element_names, table.element_names):
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
        if not array_equal(self.data, table.data):
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

                    if not array_equal(t1, t2):
                        msg += '(%s, %s, %s)    (%s, %s, %s, %s, %s, %s)  (%s, %s, %s, %s, %s, %s)\n' % (
                            eid, nid, ename,
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

    def add(self, dt, node_id, eid, ename, t1, t2, t3, r1, r2, r3):
        assert isinstance(node_id, int), node_id
        self.add_sort1(dt, eid, node_id, ename, t1, t2, t3, r1, r2, r3)

    def add_sort1(self, dt, node_id, eid, ename, t1, t2, t3, r1, r2, r3):
        assert eid is not None, eid
        assert isinstance(node_id, int), node_id
        self.node_element[self.itotal, :] = [eid, node_id]
        self.element_names[self.itotal] = ename
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

    def write_f06(self, header, page_stamp, page_num=1, f=None, is_mag_phase=False, is_sort1=True):
        msg = self._get_f06_msg()

        ntimes = self.data.shape[0]

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
                (vals2, is_all_zeros) = writeFloats13E(vals)
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
        return str(b', '.join(headers))

    def build(self):
        #print('ntimes=%s nelements=%s ntotal=%s' % (self.ntimes, self.nelements, self.ntotal))
        if self.is_built:
            return
        #self.ntotal += 5  # TODO: remove
        #print("self.ntotal = %s" % self.ntotal)
        #print("self.itotal = %s" % self.itotal)
        #print("self._ntotals = %s" % self._ntotals)

        self.is_unique = False
        if self.ntotal != min(self._ntotals):
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
            self.element_names = empty(self.ntimes, self.ntotal, dtype='S8')
        else:
            self.node_element = zeros((self.ntotal, 2), dtype='int32')
            self.element_names = empty(self.ntotal, dtype='S8')
        #[t1, t2, t3, r1, r2, r3]
        self.data = zeros((self.ntimes, self.ntotal, 6), dtype='complex64')

    def __eq__(self, table):
        assert self.is_sort1() == table.is_sort1()
        assert self.nonlinear_factor == table.nonlinear_factor
        assert self.ntotal == table.ntotal
        assert self.table_name == table.table_name, 'table_name=%r table.table_name=%r' % (self.table_name, table.table_name)
        assert self.approach_code == table.approach_code
        if not array_equal(self.node_element, table.node_element) and array_equal(self.element_names, table.element_names):
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
        if not array_equal(self.data, table.data):
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

                    if not array_equal(t1, t2):
                        msg += '(%s, %s, %s)    (%s, %s, %s, %s, %s, %s)  (%s, %s, %s, %s, %s, %s)\n' % (
                            eid, nid, ename,
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

    def add(self, dt, node_id, eid, ename, t1, t2, t3, r1, r2, r3):
        assert isinstance(node_id, int), node_id
        self.add_sort1(dt, eid, node_id, ename, t1, t2, t3, r1, r2, r3)

    def add_sort1(self, dt, node_id, eid, ename, t1, t2, t3, r1, r2, r3):
        assert eid is not None, eid
        assert isinstance(node_id, int), node_id
        self.node_element[self.itotal, :] = [eid, node_id]
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

    def write_f06(self, header, page_stamp, page_num=1, f=None, is_mag_phase=False, is_sort1=True):
        msg = self._get_f06_msg(is_mag_phase=is_mag_phase, is_sort1=is_sort1)

        ntimes = self.data.shape[0]

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
                (vals2, is_all_zeros) = writeImagFloats13E(vals, is_mag_phase)
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


class RealGridPointForces(ScalarObject):
    def __init__(self, data_code, is_sort1, isubcase, dt):
        ScalarObject.__init__(self, data_code, isubcase)
        self.force_moment = {}
        self.elemName = {}
        self.eids = {}

        self.dt = dt
        if is_sort1:
            if dt is not None:
                self.add = self.add_sort1
        else:
            assert dt is not None, 'is_sort1=%s' % is_sort1
            self.add = self.add_sort2

    def get_stats(self):
        nelements = len(self.eids)

        msg = self.get_data_code()
        if self.dt is not None:  # transient
            ntimes = len(self.force_moment)
            msg.append('  type=%s ntimes=%s nelements=%s\n'
                       % (self.__class__.__name__, ntimes, nelements))
        else:
            msg.append('  type=%s nelements=%s\n' % (self.__class__.__name__,
                                                     nelements))
        msg.append('  force_moment, elemName, eids\n')
        return msg

    def add_new_transient(self, dt):  # eKey
        """initializes the transient variables"""
        self.force_moment[dt] = {}

        # these can shockingly be different sizes, so we can't save space
        self.elemName[dt] = {}
        self.eids[dt] = {}

    def add(self, dt, eKey, eid, elemName, f1, f2, f3, m1, m2, m3):
        if eKey not in self.force_moment:
            self.eids[eKey] = []
            self.force_moment[eKey] = []
            self.elemName[eKey] = []
        self.force_moment[eKey].append(array([f1, f2, f3, m1, m2, m3], dtype='float32'))  # Fx,Fy,Fz
        self.elemName[eKey].append(elemName)
        self.eids[eKey].append(eid)

    def add_f06_data(self, dt, data):
        if dt is None:
            for (nid, eid, source, f1, f2, f3, m1, m2, m3) in data:
                if nid not in self.force_moment:
                    self.eids[nid] = []
                    self.force_moment[nid] = []
                    self.elemName[nid] = []
                self.force_moment[nid].append(array([f1, f2, f3, m1, m2, m3], dtype='float32'))
                self.elemName[nid].append(source)
                self.eids[nid].append(eid)

        else:
            if dt not in self.force_moment:
                self.force_moment[dt] = {}
            for (nid, eid, source, f1, f2, f3, m1, m2, m3) in data:
                if nid not in self.force_moment[dt]:
                    self.eids[nid] = []
                    self.force_moment[dt][nid] = []
                    self.elemName[dt][nid] = []
                    self.eids[dt][nid] = []
                self.force_moment[dt][nid].append(array([f1, f2, f3, m1, m2, m3], dtype='float32'))
                self.elemName[dt][nid].append(source)
                self.eids[dt][nid].append(eid)

    def add_sort1(self, dt, ekey, eid, elemName, f1, f2, f3, m1, m2, m3):
        if dt not in self.force_moment:
            #print "new transient"
            self.add_new_transient(dt)

        #print "%s=%s ekey=%s eid=%s elemName=%s f1=%s" %(self.data_code['name'], dt, ekey, eid, elemName, f1)
        if ekey not in self.force_moment[dt]:
            self.eids[dt][ekey] = []
            self.force_moment[dt][ekey] = []
            self.elemName[dt][ekey] = []
        self.force_moment[dt][ekey].append(array([f1, f2, f3, m1, m2, m3], dtype='float32'))
        self.elemName[dt][ekey].append(elemName)
        self.eids[dt][ekey].append(eid)

    def update_dt(self, data_code, freq):
        self.data_code = data_code
        self.apply_data_code()
        if freq is not None:
            self.log.debug("updating %s...%s=%s  isubcase=%s" % (self.data_code['name'], self.data_code['name'], freq, self.isubcase))
            self.dt = dt
            self.add_new_transient()

    def delete_transient(self, dt):
        del self.force_moment[dt]
        del self.elemName[dt]
        del self.eids[dt]

    def get_transients(self):
        k = self.force_moment.keys()
        k.sort()
        return k

    #def cleanupObj(self):
        #k = self.elemName.keys()
        #self.elemName = self.elemName[k[0]]
        #self.eids = self.eids[k[0]]

    def write_f06(self, header, page_stamp, page_num=1, f=None, is_mag_phase=False, is_sort1=True):
        if self.nonlinear_factor is not None:
            return self._write_f06_transient(header, page_stamp, page_num, f, is_mag_phase=is_mag_phase, is_sort1=is_sort1)

        msg = header + ['                                          G R I D   P O I N T   F O R C E   B A L A N C E\n',
                        ' \n',
                        '   POINT-ID    ELEMENT-ID     SOURCE             T1             T2             T3             R1             R2             R3\n', ]
              #'0     13683          3736    TRIAX6         4.996584E+00   0.0            1.203093E+02   0.0            0.0            0.0'
              #'      13683          3737    TRIAX6        -4.996584E+00   0.0           -1.203093E+02   0.0            0.0            0.0'
              #'      13683                  *TOTALS*       6.366463E-12   0.0           -1.364242E-12   0.0            0.0            0.0'
        zero = ' '
        for eKey, Force in sorted(iteritems(self.force_moment)):
            for iLoad, force in enumerate(Force):
                (f1, f2, f3, m1, m2, m3) = force
                elemName = self.elemName[eKey][iLoad]
                eid = self.eids[eKey][iLoad]
                vals = [f1, f2, f3, m1, m2, m3]
                (vals2, is_all_zeros) = writeFloats13E(vals)
                [f1, f2, f3, m1, m2, m3] = vals2
                if eid == 0:
                    eid = ''
                msg.append('%s  %8s    %10s    %-8s      %-13s  %-13s  %-13s  %-13s  %-13s  %s\n' % (zero, eKey, eid, elemName,
                                                                                                     f1, f2, f3, m1, m2, m3))
                zero = ' '
            zero = '0'
        msg.append(page_stamp % page_num)
        f.write(''.join(msg))
        return page_num

    def _write_f06_transient(self, header, page_stamp, page_num=1, f=None,
                             is_mag_phase=False, is_sort1=True):
        msg_pack = ['                                          G R I D   P O I N T   F O R C E   B A L A N C E\n',
                    ' \n',
                    '   POINT-ID    ELEMENT-ID     SOURCE             T1             T2             T3             R1             R2             R3\n', ]
              #'0     13683          3736    TRIAX6         4.996584E+00   0.0            1.203093E+02   0.0            0.0            0.0'
              #'      13683          3737    TRIAX6        -4.996584E+00   0.0           -1.203093E+02   0.0            0.0            0.0'
              #'      13683                  *TOTALS*       6.366463E-12   0.0           -1.364242E-12   0.0            0.0            0.0'
        ntimes = len(self.force_moment)
        itime = 0
        msg = ['']
        for dt, Forces in sorted(iteritems(self.force_moment)):
            zero = ' '
            header = _eigenvalue_header(self, header, itime, ntimes, dt)
            msg += header + msg_pack
            for eKey, Force in sorted(iteritems(Forces)):
                for iLoad, force in enumerate(Force):
                    (f1, f2, f3, m1, m2, m3) = force
                    elemName = self.elemName[dt][eKey][iLoad]
                    eid = self.eids[dt][eKey][iLoad]

                    vals = [f1, f2, f3, m1, m2, m3]
                    (vals2, is_all_zeros) = writeFloats13E(vals)
                    [f1, f2, f3, m1, m2, m3] = vals2
                    if eid == 0:
                        eid = ''
                    msg.append('%s  %8s    %10s    %-8s      %-13s  %-13s  %-13s  %-13s  %-13s  %s\n' % (zero, eKey, eid, elemName,
                                                                                                        f1, f2, f3, m1, m2, m3))
                    zero = ' '
                zero = '0'

            msg.append(page_stamp % page_num)
            f.write(''.join(msg))
            msg = ['']
            itime += 1
            page_num += 1
        return page_num - 1


class ComplexGridPointForces(ScalarObject):
    def __init__(self, data_code, is_sort1, isubcase, freq=None):
        ScalarObject.__init__(self, data_code, isubcase)
        raise NotImplementedError()

    def get_stats(self):
        nelements = len(self.eids)

        msg = self.get_data_code()
        if self.dt is not None:  # transient
            ntimes = len(self.forces)
            msg.append('  type=%s ntimes=%s nelements=%s\n'
                       % (self.__class__.__name__, ntimes, nelements))
        else:
            msg.append('  type=%s nelements=%s\n' % (self.__class__.__name__,
                                                     nelements))
        msg.append('  forces, moments, elemName, eids\n')
        return msg
