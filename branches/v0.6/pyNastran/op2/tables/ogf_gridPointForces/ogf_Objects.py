from struct import pack
from pyNastran.op2.resultObjects.op2_Objects import scalarObject
from pyNastran.f06.f06_formatting import writeFloats13E


class gridPointForcesObject(scalarObject):

    def __init__(self, data_code, is_sort1, isubcase, dt=None):
        scalarObject.__init__(self, data_code, isubcase)
        self.forces = {}
        self.moments = {}
        self.elemName = {}
        self.eids = {}

        self.dt = dt
        if is_sort1:
            if dt is not None:
                self.add = self.add_sort1
        else:
            assert dt is not None
            self.add = self.addSort2

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

    def add_new_transient(self, dt):  # eKey
        """initializes the transient variables"""
        self.forces[dt] = {}
        self.moments[dt] = {}
        #self.elemName[dt] = {}
        #self.eids[dt] = {}
        #self.elemName[eKey] = []
        #self.eids[eKey] = []
        self.elemName = {}
        self.eids = {}

    def add(self, dt, eKey, eid, elemName, f1, f2, f3, m1, m2, m3):
        if eKey not in self.forces:
            self.eids[eKey] = []
            self.forces[eKey] = []
            self.moments[eKey] = []
            self.elemName[eKey] = []
        self.forces[eKey].append([f1, f2, f3])  # Fx,Fy,Fz
        self.moments[eKey].append([m1, m2, m3])  # Mx,My,Mz
        self.elemName[eKey].append(elemName)
        self.eids[eKey].append(eid)

    def add_sort1(self, dt, eKey, eid, elemName, f1, f2, f3, m1, m2, m3):
        if dt not in self.forces:
            #print "new transient"
            self.add_new_transient(dt)

        #print "%s=%s eKey=%s eid=%s elemName=%s f1=%s" %(self.data_code['name'],dt,eKey,eid,elemName,f1)
        if eKey not in self.forces[dt]:
            self.eids[eKey] = []
            self.forces[dt][eKey] = []
            self.moments[dt][eKey] = []
            self.elemName[eKey] = []
        self.forces[dt][eKey].append([f1, f2, f3])  # Mx,My,Mz
        self.moments[dt][eKey].append([m1, m2, m3])  # Fx,Fy,Fz
        self.elemName[eKey].append(elemName)
        self.eids[eKey].append(eid)

    def update_dt(self, data_code, freq):
        self.data_code = data_code
        self.apply_data_code()
        if freq is not None:
            self.log.debug("updating %s...%s=%s  isubcase=%s" % (self.data_code['name'], self.data_code['name'], freq, self.isubcase))
            self.dt = dt
            self.add_new_transient()

    def delete_transient(self, dt):
        del self.forces[dt]
        del self.moments[dt]
        del self.elemName[dt]

    def get_transients(self):
        k = self.forces.keys()
        k.sort()
        return k

    #def cleanupObj(self):
        #k = self.elemName.keys()
        #self.elemName = self.elemName[k[0]]
        #self.eids = self.eids[k[0]]

    def write_f06(self, header, pageStamp, pageNum=1, f=None):
        if self.nonlinear_factor is not None:
            return self._write_f06_transient(header, pageStamp, pageNum, f)

        msg = header + ['                                          G R I D   P O I N T   F O R C E   B A L A N C E\n',
                        ' \n',
                        '   POINT-ID    ELEMENT-ID     SOURCE             T1             T2             T3             R1             R2             R3\n', ]
              #'0     13683          3736    TRIAX6         4.996584E+00   0.0            1.203093E+02   0.0            0.0            0.0'
              #'      13683          3737    TRIAX6        -4.996584E+00   0.0           -1.203093E+02   0.0            0.0            0.0'
              #'      13683                  *TOTALS*       6.366463E-12   0.0           -1.364242E-12   0.0            0.0            0.0'
        zero = ' '
        for eKey, Force in sorted(self.forces.iteritems()):
            for iLoad, force in enumerate(Force):
                (f1, f2, f3) = force
                (m1, m2, m3) = self.moments[eKey][iLoad]
                (elemName) = self.elemName[eKey][iLoad]
                eid = self.eids[eKey][iLoad]
                vals = [f1, f2, f3, m1, m2, m3]
                (vals2, isAllZeros) = writeFloats13E(vals)
                [f1, f2, f3, m1, m2, m3] = vals2
                if eid == 0:
                    eid = ''
                msg.append('%s  %8s    %10s    %-8s      %s  %s  %s  %s  %s  %-s\n' % (zero, eKey, eid, elemName, f1, f2, f3, m1, m2, m3.rstrip()))
                zero = ' '
            zero = '0'

        msg.append(pageStamp + str(pageNum) + '\n')
        if f is not None:
            f.write(''.join(msg))
            msg = ['']
        return (''.join(msg), pageNum)

    def _write_f06_transient(self, header, pageStamp, pageNum=1, f=None):
        msg = header + ['                                          G R I D   P O I N T   F O R C E   B A L A N C E\n',
                        ' \n',
                        '   POINT-ID    ELEMENT-ID     SOURCE             T1             T2             T3             R1             R2             R3\n', ]
              #'0     13683          3736    TRIAX6         4.996584E+00   0.0            1.203093E+02   0.0            0.0            0.0'
              #'      13683          3737    TRIAX6        -4.996584E+00   0.0           -1.203093E+02   0.0            0.0            0.0'
              #'      13683                  *TOTALS*       6.366463E-12   0.0           -1.364242E-12   0.0            0.0            0.0'
        for dt, Forces in sorted(self.forces.iteritems()):
            zero = ' '
            for eKey, Force in sorted(Forces.iteritems()):
                for iLoad, force in enumerate(Force):
                    (f1, f2, f3) = force
                    (m1, m2, m3) = self.moments[dt][eKey][iLoad]
                    (elemName) = self.elemName[eKey][iLoad]
                    eid = self.eids[eKey][iLoad]

                    vals = [f1, f2, f3, m1, m2, m3]
                    (vals2, isAllZeros) = writeFloats13E(vals)
                    [f1, f2, f3, m1, m2, m3] = vals2
                    if eid == 0:
                        eid = ''

                    msg.append('%s  %8s    %10s    %8s      %s  %s  %s  %s  %s  %-s\n' % (zero, eKey, eid, elemName, f1, f2, f3, m1, m2, m3))
                    zero = ' '
                zero = '0'

            msg.append(pageStamp + str(pageNum) + '\n')
            if f is not None:
                f.write(''.join(msg))
                msg = ['']
            pageNum += 1
        return (''.join(msg), pageNum - 1)


class complexGridPointForcesObject(scalarObject):
    def __init__(self, data_code, is_sort1, isubcase, freq=None):
        scalarObject.__init__(self, data_code, isubcase)
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

