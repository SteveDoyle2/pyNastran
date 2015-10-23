class ComplexCBarForce(ScalarObject):  # 34-CBAR
    def __init__(self, data_code, is_sort1, isubcase, dt):
        ScalarObject.__init__(self, data_code, isubcase)
        self.bendingMomentA = {}
        self.bendingMomentB = {}
        self.shear = {}
        self.axial = {}
        self.torque = {}
        self.dt = dt

    def get_stats(self):
        msg = ['  '] + self.get_data_code()
        if self.dt is not None:  # transient
            ntimes = len(self.torque)
            time0 = get_key0(self.torque)
            nelements = len(self.torque[time0])
            msg.append('  type=%s ntimes=%s nelements=%s\n'
                       % (self.__class__.__name__, ntimes, nelements))
        else:
            nelements = len(self.torque)
            msg.append('  type=%s nelements=%s\n' % (self.__class__.__name__,
                                                     nelements))
        msg.append('  bendingMomentA, bendingMomentB, shes, axial, torque\n')
        return msg

    def add_new_transient(self, dt):
        self.dt = dt
        self.bendingMomentA[dt] = {}
        self.bendingMomentB[dt] = {}
        self.shear[dt] = {}
        self.axial[dt] = {}
        self.torque[dt] = {}

    def add_sort1(self, dt, eid, bm1a, bm2a, bm1b, bm2b, ts1, ts2, af, trq):
        if dt not in self.axial:
            self.add_new_transient(dt)
        self.bendingMomentA[dt][eid] = [bm1a, bm2a]
        self.bendingMomentB[dt][eid] = [bm1b, bm2b]
        self.shear[dt][eid] = [ts1, ts2]
        self.axial[dt][eid] = af
        self.torque[dt][eid] = trq

    def add_sort2(self, eid, dt, bm1a, bm2a, bm1b, bm2b, ts1, ts2, af, trq):
        if dt not in self.axial:
            self.add_new_transient(dt)
        self.bendingMomentA[dt][eid] = [bm1a, bm2a]
        self.bendingMomentB[dt][eid] = [bm1b, bm2b]
        self.shear[dt][eid] = [ts1, ts2]
        self.axial[dt][eid] = af
        self.torque[dt][eid] = trq

    def _write_f06_transient(self, header, page_stamp, page_num=1, f=None, is_mag_phase=False, is_sort1=True):
        words = ['                             C O M P L E X   F O R C E S   I N   B A R   E L E M E N T S   ( C B A R )\n',
                 '                                                          (REAL/IMAGINARY)\n',
                 '0    ELEMENT         BEND-MOMENT END-A            BEND-MOMENT END-B                - SHEAR -               AXIAL\n',
                 '       ID.         PLANE 1       PLANE 2        PLANE 1       PLANE 2        PLANE 1       PLANE 2         FORCE         TORQUE\n']

        msg = []
        for dt, bendA in sorted(iteritems(self.bendingMomentA)):
            header[1] = ' %s = %10.4E\n' % (self.data_code['name'], dt)
            msg += header + words
            for eid in sorted(bendA):
                bm1a, bm2a = self.bendingMomentA[dt][eid]
                bm1b, bm2b = self.bendingMomentB[dt][eid]
                ts1, ts2 = self.shear[dt][eid]
                af = self.axial[dt][eid]
                trq = self.torque[dt][eid]
                (vals2, is_all_zeros) = writeImagFloats13E([bm1a, bm2a, bm1b, bm2b, ts1, ts2, af, trq], is_mag_phase)
                [bm1ar, bm2ar, bm1br, bm2br, ts1r, ts2r, afr, trqr,
                 bm1ai, bm2ai, bm1bi, bm2bi, ts1i, ts2i, afi, trqi] = vals2
                msg.append('     %8i    %-13s %-13s  %-13s %-13s  %-13s %-13s  %-13s  %s\n' % (eid, bm1ar, bm2ar, bm1br, bm2br, ts1r, ts2r, afr, trqr))
                msg.append('     %8s    %-13s %-13s  %-13s %-13s  %-13s %-13s  %-13s  %s\n' % ('', bm1ai, bm2ai, bm1bi, bm2bi, ts1i, ts2i, afi, trqi))
            msg.append(page_stamp % page_num)
            f.write(''.join(msg))
            msg = ['']
            page_num += 1
        return page_num - 1

    def write_f06(self, header, page_stamp, page_num=1, f=None, is_mag_phase=False, is_sort1=True):
        if self.nonlinear_factor is not None:
            return self._write_f06_transient(header, page_stamp, page_num, f, is_mag_phase=is_mag_phase, is_sort1=is_sort1)
        raise RuntimeError('is this used...')
        msg = header + ['                             C O M P L E X   F O R C E S   I N   B A R   E L E M E N T S   ( C B A R )\n',
                        '                                                          (REAL/IMAGINARY)\n',
                        '0    ELEMENT         BEND-MOMENT END-A            BEND-MOMENT END-B                - SHEAR -               AXIAL\n',
                        '       ID.         PLANE 1       PLANE 2        PLANE 1       PLANE 2        PLANE 1       PLANE 2         FORCE         TORQUE\n']

        for eid in sorted(self.bendingMomentA):
            bm1a, bm2a = self.bendingMomentA[eid]
            bm1b, bm2b = self.bendingMomentB[eid]
            ts1, ts2 = self.shear[eid]
            af = self.axial[eid]
            trq = self.torque[eid]
            (vals2, is_all_zeros) = writeImagFloats13E([bm1a, bm2a, bm1b, bm2b, ts1, ts2, af, trq], is_mag_phase)
            [bm1ar, bm2ar, bm1br, bm2br, ts1r, ts2r, afr, trqr,
             bm1ai, bm2ai, bm1bi, bm2bi, ts1i, ts2i, afi, trqi] = vals2
            msg.append('     %8i    %-13s %-13s  %-13s %-13s  %-13s %-13s  %-13s  %s\n' % (eid, bm1ar, bm2ar, bm1br, bm2br, ts1r, ts2r, afr, trqr))
            msg.append('     %8s    %-13s %-13s  %-13s %-13s  %-13s %-13s  %-13s  %s\n' % ('', bm1ai, bm2ai, bm1bi, bm2bi, ts1i, ts2i, afi, trqi))
        msg.append(page_stamp % page_num)
        f.write(''.join(msg))
        return page_num


class ComplexCBushForce(ScalarObject):  # 102-CBUSH
    def __init__(self, data_code, is_sort1, isubcase, dt):
        ScalarObject.__init__(self, data_code, isubcase)
        self.force = {}
        self.moment = {}

        self.dt = dt
        if is_sort1:
            if dt is not None:
                self.add = self.add_sort1
        else:
            assert dt is not None
            self.add = self.add_sort2

    def get_stats(self):
        msg = ['  '] + self.get_data_code()
        if self.dt is not None:  # transient
            ntimes = len(self.force)
            time0 = get_key0(self.force)
            nelements = len(self.force[time0])
            msg.append('  type=%s ntimes=%s nelements=%s\n'
                       % (self.__class__.__name__, ntimes, nelements))
        else:
            nelements = len(self.force)
            msg.append('  type=%s nelements=%s\n' % (self.__class__.__name__,
                                                     nelements))
        msg.append('  force, moment\n')
        return msg

    def add_new_transient(self, dt):
        self.dt = dt
        self.force[dt] = {}
        self.moment[dt] = {}

    def add_sort1(self, dt, eid, fx, fy, fz, mx, my, mz):
        if dt not in self.force:
            self.add_new_transient(dt)
        self.force[dt][eid] = [fx, fy, fz]
        self.moment[dt][eid] = [mx, my, mz]

    def add_sort2(self, eid, dt, fx, fy, fz, mx, my, mz):
        if dt not in self.force:
            self.add_new_transient(dt)
        self.force[dt][eid] = [fx, fy, fz]
        self.moment[dt][eid] = [mx, my, mz]


