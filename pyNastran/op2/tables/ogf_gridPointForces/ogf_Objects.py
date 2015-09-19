from six import iteritems
from numpy import array
from struct import pack
from pyNastran.op2.resultObjects.op2_Objects import ScalarObject
from pyNastran.f06.f06_formatting import writeFloats13E, _eigenvalue_header


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
            self.add = self.addSort2

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
