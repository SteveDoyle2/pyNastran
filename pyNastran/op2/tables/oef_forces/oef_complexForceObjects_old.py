from six import iteritems
from pyNastran.op2.resultObjects.op2_Objects import ScalarObject
from pyNastran.f06.f06_formatting import write_imag_floats_13e, get_key0, write_floats_13e

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
                vals2 = write_imag_floats_13e([bm1a, bm2a, bm1b, bm2b, ts1, ts2, af, trq], is_mag_phase)
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
            vals2 = write_imag_floats_13e([bm1a, bm2a, bm1b, bm2b, ts1, ts2, af, trq], is_mag_phase)
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


class ComplexRodForce(ScalarObject):  # 1-ROD, 3-TUBE, 10-CONROD
    def __init__(self, data_code, is_sort1, isubcase, dt):
        ScalarObject.__init__(self, data_code, isubcase)
        self.axial_force = {}
        self.torque = {}

        self.dt = dt
        if is_sort1:
            if dt is not None:
                self.add = self.add_sort1
        else:
            assert dt is not None
            self.add = self.add_sort2

    def add_new_transient(self, dt):
        self.dt = dt
        self.axial_force[dt] = {}
        self.torque[dt] = {}

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
        msg.append(' axial_force, torque\n')
        return msg

    def add(self, dt, data):
        [eid, axial_force, torque] = data
        self.axial_force[eid] = axial_force
        self.torque[eid] = torque

    def add_sort1(self, dt, data):
        [eid, axial_force, torque] = data
        if dt not in self.axial_force:
            self.add_new_transient(dt)
        self.axial_force[dt][eid] = axial_force
        self.torque[dt][eid] = torque

    def add_sort2(self, eid, data):
        [dt, axial_force, torque] = data
        if dt not in self.axial_force:
            self.add_new_transient(dt)
        self.axial_force[dt][eid] = axial_force
        self.torque[dt][eid] = torque


class ComplexCBeamForce(ScalarObject):  # 2-CBEAM
    def __init__(self, data_code, is_sort1, isubcase, dt):
        ScalarObject.__init__(self, data_code, isubcase)
        self.bendingMoment = {}
        self.shear = {}
        self.axial = {}
        self.totalTorque = {}
        self.warpingTorque = {}

        self.dt = dt
        if is_sort1:
            if dt is not None:
                self.add_new_element = self.addNewElementSort1
                self.add = self.add_sort1
        else:
            assert dt is not None
            self.add_new_element = self.addNewElementSort2
            self.add = self.add_sort2

    def get_stats(self):
        msg = ['  ']
        msg += self.get_data_code()
        if self.dt is not None:  # transient
            ntimes = len(self.shear)
            time0 = get_key0(self.shear)
            nelements = len(self.shear[time0])
            msg.append('  type=%s ntimes=%s nelements=%s\n'
                       % (self.__class__.__name__, ntimes, nelements))
        else:
            nelements = len(self.shear)
            msg.append('  type=%s nelements=%s\n' % (self.__class__.__name__,
                                                     nelements))
        msg.append('  bendingMoment, shear, axial, totalTorque, '
                   'warpingTorque\n')
        return msg

    def add_new_transient(self, dt):
        self.dt = dt
        self.bendingMoment[dt] = {}
        self.shear[dt] = {}
        self.axial[dt] = {}
        self.totalTorque[dt] = {}
        self.warpingTorque[dt] = {}

    def add_new_element(self, dt, data):
        [eid, nid, sd, bm1, bm2, ts1, ts2, af, ttrq, wtrq] = data
        #print "CBEAM addnew",data
        self.bendingMoment[eid] = {sd: [bm1, bm2]}
        self.shear[eid] = {sd: [ts1, ts2]}
        self.axial[eid] = {sd: af}
        self.totalTorque[eid] = {sd: ttrq}
        self.warpingTorque[eid] = {sd: wtrq}

    def add(self, dt, data):
        [eid, nid, sd, bm1, bm2, ts1, ts2, af, ttrq, wtrq] = data
        #print "CBEAM add   ",data
        self.bendingMoment[eid][sd] = [bm1, bm2]
        self.shear[eid][sd] = [ts1, ts2]
        self.axial[eid][sd] = af
        self.totalTorque[eid][sd] = ttrq
        self.warpingTorque[eid][sd] = wtrq

    def addNewElementSort1(self, dt, data):
        [eid, nid, sd, bm1, bm2, ts1, ts2, af, ttrq, wtrq] = data
        self._fillNewObject(
            dt, eid, nid, sd, bm1, bm2, ts1, ts2, af, ttrq, wtrq)

    def add_sort1(self, dt, data):
        [eid, nid, sd, bm1, bm2, ts1, ts2, af, ttrq, wtrq] = data
        self._fill_object(dt, eid, nid, sd, bm1, bm2, ts1, ts2, af, ttrq, wtrq)

    def addNewElementSort2(self, eid, data):
        [dt, nid, sd, bm1, bm2, ts1, ts2, af, ttrq, wtrq] = data
        self._fillNewObject(
            dt, eid, nid, sd, bm1, bm2, ts1, ts2, af, ttrq, wtrq)

    def add_sort2(self, eid, data):
        [dt, nid, sd, bm1, bm2, ts1, ts2, af, ttrq, wtrq] = data
        self._fill_object(dt, eid, nid, sd, bm1, bm2, ts1, ts2, af, ttrq, wtrq)

    def _fill_object(self, dt, eid, nid, sd, bm1, bm2, ts1, ts2, af, ttrq, wtrq):
        #if dt not in self.axial:
            #self.add_new_transient(dt)
        self.bendingMoment[dt][eid][sd] = [bm1, bm2]
        self.shear[dt][eid][sd] = [ts1, ts2]
        self.axial[dt][eid][sd] = af
        self.totalTorque[dt][eid][sd] = ttrq
        self.warpingTorque[dt][eid][sd] = wtrq

    def _fillNewObject(self, dt, eid, nid, sd, bm1, bm2, ts1, ts2, af, ttrq, wtrq):
        if dt not in self.axial:
            self.add_new_transient(dt)
        self.bendingMoment[dt][eid] = {sd: [bm1, bm2]}
        self.shear[dt][eid] = {sd: [ts1, ts2]}
        self.axial[dt][eid] = {sd: af}
        self.totalTorque[dt][eid] = {sd: ttrq}
        self.warpingTorque[dt][eid] = {sd: wtrq}

class ComplexCShearForce(ScalarObject):  # 4-CSHEAR
    def __init__(self, data_code, is_sort1, isubcase, dt):
        ScalarObject.__init__(self, data_code, isubcase)
        self.force41 = {}
        self.force14 = {}
        self.force21 = {}
        self.force12 = {}
        self.force32 = {}
        self.force23 = {}
        self.force43 = {}
        self.force34 = {}
        self.kickForce1 = {}
        self.kickForce2 = {}
        self.kickForce3 = {}
        self.kickForce4 = {}
        self.shear12 = {}
        self.shear23 = {}
        self.shear34 = {}
        self.shear41 = {}

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
            ntimes = len(self.shear12)
            time0 = get_key0(self.shear12)
            nelements = len(self.shear12[time0])
            msg.append('  type=%s ntimes=%s nelements=%s\n'
                       % (self.__class__.__name__, ntimes, nelements))
        else:
            nelements = len(self.shear12)
            msg.append('  type=%s nelements=%s\n' % (self.__class__.__name__,
                                                     nelements))
        msg.append('  force41, force14, force21, force12, force32, force23, '
                   '  force 43, force34, kickForce1, kickForce2, kickForce3, '
                   '  kickForce4, shear12, shear23, shear34, shear41\n')
        return msg

    def add_new_transient(self, dt):
        self.force41[dt] = {}
        self.force14[dt] = {}
        self.force21[dt] = {}
        self.force12[dt] = {}
        self.force32[dt] = {}
        self.force23[dt] = {}
        self.force43[dt] = {}
        self.force34[dt] = {}
        self.kickForce1[dt] = {}
        self.kickForce2[dt] = {}
        self.kickForce3[dt] = {}
        self.kickForce4[dt] = {}
        self.shear12[dt] = {}
        self.shear23[dt] = {}
        self.shear34[dt] = {}
        self.shear41[dt] = {}

    def add(self, dt, data):
        [eid, f41, f21, f12, f32, f23, f43, f34, f14,
         kf1, s12, kf2, s23, kf3, s34, kf4, s41] = data
        #self.eType[eid] = eType
        self.force41[eid] = f41
        self.force14[eid] = f14
        self.force21[eid] = f21
        self.force12[eid] = f12
        self.force32[eid] = f32
        self.force23[eid] = f23
        self.force43[eid] = f43
        self.force34[eid] = f34
        self.kickForce1[eid] = kf1
        self.kickForce2[eid] = kf2
        self.kickForce3[eid] = kf3
        self.kickForce4[eid] = kf4
        self.shear12[eid] = s12
        self.shear23[eid] = s23
        self.shear34[eid] = s34
        self.shear41[eid] = s41

    def add_sort1(self, dt, data):
        [eid, f41, f21, f12, f32, f23, f43, f34, f14,
         kf1, s12, kf2, s23, kf3, s34, kf4, s41] = data
        self._fill_object(dt, eid, f41, f21, f12, f32, f23, f43, f34, f14,
                          kf1, s12, kf2, s23, kf3, s34, kf4, s41)

    def add_sort2(self, eid, data):
        [dt, f41, f21, f12, f32, f23, f43, f34, f14,
         kf1, s12, kf2, s23, kf3, s34, kf4, s41] = data
        self._fill_object(dt, eid, f41, f21, f12, f32, f23, f43, f34, f14,
                          kf1, s12, kf2, s23, kf3, s34, kf4, s41)

    def _fill_object(self, dt, eid, f41, f21, f12, f32, f23, f43, f34, f14,
                     kf1, s12, kf2, s23, kf3, s34, kf4, s41):
        if dt not in self.force41:
            self.add_new_transient(dt)
        self.force41[dt][eid] = f41
        self.force14[dt][eid] = f14
        self.force21[dt][eid] = f21
        self.force12[dt][eid] = f12
        self.force32[dt][eid] = f32
        self.force23[dt][eid] = f23
        self.force43[dt][eid] = f43
        self.force34[dt][eid] = f34
        self.kickForce1[dt][eid] = kf1
        self.kickForce2[dt][eid] = kf2
        self.kickForce3[dt][eid] = kf3
        self.kickForce4[dt][eid] = kf4
        self.shear12[dt][eid] = s12
        self.shear23[dt][eid] = s23
        self.shear34[dt][eid] = s34
        self.shear41[dt][eid] = s41


class ComplexSpringForce(ScalarObject):  # 11-CELAS1,12-CELAS2,13-CELAS3, 14-CELAS4
    def __init__(self, data_code, is_sort1, isubcase, dt):
        ScalarObject.__init__(self, data_code, isubcase)
        self.force = {}
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
        msg.append('  force\n')
        return msg

    def add_new_transient(self, dt):
        self.dt = dt
        self.force[dt] = {}

    def add(self, dt, data):
        [eid, force] = data
        self.force[eid] = force

    def add_sort1(self, dt, data):
        [eid, force] = data
        if dt not in self.force:
            self.add_new_transient(dt)
        self.force[dt][eid] = force

    def add_sort2(self, eid, data):
        [dt, force] = data
        if dt not in self.force:
            self.add_new_transient(dt)
        self.force[dt][eid] = force

    def write_f06(self, header, page_stamp, page_num=1, f=None, is_mag_phase=False, is_sort1=True):
        if self.nonlinear_factor is not None:
            return self._write_f06_transient(header, page_stamp, page_num, f, is_mag_phase=is_mag_phase, is_sort1=is_sort1)
        msg = header + ['                         C O M P L E X   F O R C E S   I N   S C A L A R   S P R I N G S   ( C E L A S 4 )\n',
                        '                                                          (REAL/IMAGINARY)\n',
                        ' \n',
                        '            FREQUENCY                    FORCE                        FREQUENCY                    FORCE\n']
        forces = []
        elements = []
        line = '   '
        for eid, force in sorted(iteritems(self.force)):
            elements.append(eid)
            forces.append(force)
            #pack.append(eid)
            #pack.append(f)
            line += '%-13s  %-13s / %-13s     ' % (eid, force.real, force.imag)
            if len(forces) == 3:
                msg.append(line.rstrip() + '\n')
        if forces:
            msg.append(line.rstrip() + '\n')
        msg.append(page_stamp % page_num)

        f.write(''.join(msg))
        return page_num

    def _write_f06_transient(self, header, page_stamp, page_num=1, f=None, is_mag_phase=False, is_sort1=True):
        words = ['                         C O M P L E X   F O R C E S   I N   S C A L A R   S P R I N G S   ( C E L A S 4 )\n',
                 '                                                          (REAL/IMAGINARY)\n',
                 ' \n',
                 '                ELEMENT                                                   ELEMENT\n',
                 '                  ID.                    FORCE                              ID.                    FORCE\n']
        msg = []
        for dt, Force in sorted(self.force.items()):
            header[1] = ' %s = %10.4E\n' % (self.data_code['name'], dt)
            msg += header + words
            #packs = []
            forces = []
            elements = []
            line = ''
            for eid, force in sorted(Force.items()):
                elements.append(eid)
                forces.append(force)
                #pack.append(eid)
                #pack.append(f)
                [forceReal, forceImag] = write_floats_13e([force.real, force.imag])

                line += '          %13s      %-13s / %-13s' % (eid, forceReal, forceImag)
                if len(forces) == 2:
                    msg.append(line.rstrip() + '\n')
                    line = ''
                    forces = []
            if forces:
                msg.append(line.rstrip() + '\n')
            msg.append(page_stamp % page_num)

            f.write(''.join(msg))
            msg = ['']
            page_num += 1
        return page_num - 1


class ComplexDamperForce(ScalarObject):  # 20-CDAMP1,21-CDAMP2,22-CDAMP3,23-CDAMP4
    def __init__(self, data_code, is_sort1, isubcase, dt):
        ScalarObject.__init__(self, data_code, isubcase)
        self.force = {}

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
        msg.append('  force\n')
        return msg

    def add_new_transient(self, dt):
        self.dt = dt
        self.force[dt] = {}

    def add(self, dt, data):
        [eid, force] = data
        self.force[eid] = force

    def add_sort1(self, dt, data):
        [eid, force] = data
        if dt not in self.force:
            self.add_new_transient(dt)
        self.force[dt][eid] = force

    def add_sort2(self, eid, data):
        [dt, force] = data
        if dt not in self.force:
            self.add_new_transient(dt)
        self.force[dt][eid] = force


class ComplexViscForce(ScalarObject):  # 24-CVISC
    def __init__(self, data_code, is_sort1, isubcase, dt):
        ScalarObject.__init__(self, data_code, isubcase)
        self.axial_force = {}
        self.torque = {}

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
            ntimes = len(self.torque)
            time0 = get_key0(self.torque)
            nelements = len(self.torque[time0])
            msg.append('  type=%s ntimes=%s nelements=%s\n'
                       % (self.__class__.__name__, ntimes, nelements))
        else:
            nelements = len(self.torque)
            msg.append('  type=%s nelements=%s\n' % (self.__class__.__name__,
                                                     nelements))
        msg.append(' axial_force, torque\n')
        return msg

    def add_new_transient(self, dt):
        self.dt = dt
        self.axial_force[dt] = {}
        self.torque[dt] = {}

    def add(self, dt, eid, axial_force, torque):
        self.axial_force[eid] = axial_force
        self.torque[eid] = torque

    def add_sort1(self, dt, eid, axial_force, torque):
        if dt not in self.axial_force:
            self.add_new_transient(dt)
        self.axial_force[dt][eid] = axial_force
        self.torque[dt][eid] = torque

    def add_sort2(self, eid, dt, axial_force, torque):
        if dt not in self.axial_force:
            self.add_new_transient(dt)
        self.axial_force[dt][eid] = axial_force
        self.torque[dt][eid] = torque


class ComplexPlateForce(ScalarObject):  # 33-CQUAD4, 74-CTRIA3
    def __init__(self, data_code, is_sort1, isubcase, dt):
        ScalarObject.__init__(self, data_code, isubcase)
        self.mx = {}
        self.my = {}
        self.mxy = {}
        self.bmx = {}
        self.bmy = {}
        self.bmxy = {}
        self.tx = {}
        self.ty = {}

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
            ntimes = len(self.mx)
            time0 = get_key0(self.mx)
            nelements = len(self.mx[time0])
            msg.append('  type=%s ntimes=%s nelements=%s\n'
                       % (self.__class__.__name__, ntimes, nelements))
        else:
            nelements = len(self.mx)
            msg.append('  type=%s nelements=%s\n' % (self.__class__.__name__,
                                                     nelements))
        msg.append('  mx, my, mxy, bmx, bmy, bmxy, tx, ty\n')
        return msg

    def add_new_transient(self, dt):
        self.dt = dt
        self.mx[dt] = {}
        self.my[dt] = {}
        self.mxy[dt] = {}
        self.bmx[dt] = {}
        self.bmy[dt] = {}
        self.bmxy[dt] = {}
        self.tx[dt] = {}
        self.ty[dt] = {}

    def add(self, dt, eid, mx, my, mxy, bmx, bmy, bmxy, tx, ty):
        self.mx[eid] = mx
        self.my[eid] = my
        self.mxy[eid] = mxy
        self.bmx[eid] = bmx
        self.bmy[eid] = bmy
        self.bmxy[eid] = bmxy
        self.tx[eid] = tx
        self.ty[eid] = ty

    def add_sort1(self, dt, eid, mx, my, mxy, bmx, bmy, bmxy, tx, ty):
        if dt not in self.mx:
            self.add_new_transient(dt)
        self.mx[dt][eid] = mx
        self.my[dt][eid] = my
        self.mxy[dt][eid] = mxy
        self.bmx[dt][eid] = bmx
        self.bmy[dt][eid] = bmy
        self.bmxy[dt][eid] = bmxy
        self.tx[dt][eid] = tx
        self.ty[dt][eid] = ty

    def add_sort2(self, eid, dt, mx, my, mxy, bmx, bmy, bmxy, tx, ty):
        if dt not in self.mx:
            self.add_new_transient(dt)
        self.mx[dt][eid] = mx
        self.my[dt][eid] = my
        self.mxy[dt][eid] = mxy
        self.bmx[dt][eid] = bmx
        self.bmy[dt][eid] = bmy
        self.bmxy[dt][eid] = bmxy
        self.tx[dt][eid] = tx
        self.ty[dt][eid] = ty

    def _get_plate_msg(self, is_mag_phase, is_sort1):
        loads = ['    ELEMENT                - MEMBRANE  FORCES -                        - BENDING MOMENTS -               - TRANSVERSE SHEAR FORCES -\n'
                 '      ID              FX            FY            FXY             MX            MY            MXY             QX            QY\n',]
        if is_mag_phase:
            mag_real = ['                                                         (MAGNITUDE/PHASE)\n \n']
        else:
            mag_real = ['                                                          (REAL/IMAGINARY)\n', ' \n']

        cquad4_bilinear = ['                  C O M P L E X   F O R C E S   I N   Q U A D R I L A T E R A L   E L E M E N T S   ( Q U A D 4 )\n']  # good
        cquad4_linear = ['                  C O M P L E X   F O R C E S   I N   Q U A D R I L A T E R A L   E L E M E N T S   ( Q U A D 4 )\n']  # good
        ctria3 = ['                     C O M P L E X   F O R C E S   I N   T R I A N G U L A R   E L E M E N T S   ( T R I A 3 )\n']  # good
        cquad8 = ['                  C O M P L E X   F O R C E S   I N   Q U A D R I L A T E R A L   E L E M E N T S   ( Q U A D 8 )\n']
        cquadr = ['                  C O M P L E X   F O R C E S   I N   Q U A D R I L A T E R A L   E L E M E N T S   ( Q U A D R )\n']
        ctria6 = ['                     C O M P L E X   F O R C E S   I N   T R I A N G U L A R   E L E M E N T S   ( T R I A 6 )\n']
        ctriar = ['                     C O M P L E X   F O R C E S   I N   T R I A N G U L A R   E L E M E N T S   ( T R I A R )\n']

        is_bilinear = False
        if self.element_type == 144: # CQUAD4
            msg = cquad4_linear + mag_real + loads
        elif self.element_type == 33: # CQUAD4
            msg = cquad4_bilinear + mag_real + loads
        elif self.element_type == 64:  #CQUAD8
            msg = cquad8 + mag_real + loads
        elif self.element_type == 82:  # CQUADR
            msg = cquadr + mag_real + loads
        elif self.element_type == 74: # CTRIA3
            msg = ctria3 + mag_real + loads
        elif self.element_type == 75:  # CTRIA6
            msg = ctria6 + mag_real + loads
        elif self.element_type == 70:  # CTRIAR
            msg = ctriar + mag_real + loads
        else:
            raise NotImplementedError('name=%r type=%s' % (self.element_name, self.element_type))
        return msg

    def _write_f06_transient(self, header, page_stamp, page_num=1, f=None, is_mag_phase=False, is_sort1=True):
        msg_pack = self._get_plate_msg(is_mag_phase, is_sort1)
        dts = list(self.mx.keys())
        dts.sort()
        ntimes = len(dts)
        dt0 = dts[0]
        eids = sorted(self.mx[dt0])

        name = self.data_code['name']
        for dt in dts:
            header[1] = ' %s = %10.4E\n' % (name, dt)
            #print(header)
            #print(msg_pack)
            msg = header + msg_pack
            f.write(''.join(msg))
            for eid in eids:
                mx = self.mx[dt][eid]
                my = self.my[dt][eid]
                mxy = self.mxy[dt][eid]

                bmx = self.bmx[dt][eid]
                bmy = self.bmy[dt][eid]
                bmxy = self.bmxy[dt][eid]

                tx = self.tx[dt][eid]
                ty = self.ty[dt][eid]

                [mxr, myr, mxyr, bmxr, bmyr, bmxyr, txr, tyr,
                 mxi, myi, mxyi, bmxi, bmyi, bmxyi, txi, tyi] = write_imag_floats_13e([mx, my, mxy, bmx, bmy, bmxy, tx, ty], is_mag_phase)

                #"""
                    #ELEMENT                - MEMBRANE  FORCES -                        - BENDING MOMENTS -               - TRANSVERSE SHEAR FORCES -
                      #ID              FX            FY            FXY             MX            MY            MXY             QX            QY
                #0       564       1.543439E+03  7.311177E+02  1.322702E+02    1.080178E+00  1.699104E+00  2.618547E-01    3.877034E+01  4.518554E+00
                                  #358.3129      358.0245      177.5593        177.5292      178.2112        0.0907        358.1465      179.4567
                #"""
                #                fx     fy     fxy     mx     my     mxy    qx      qy
                msg = '0  %8i   %-13s  %-13s  %-13s   %-13s  %-13s  %-13s   %-13s  %s\n' % (eid, mxr, myr, mxyr, bmxr, bmyr, bmxyr, txr, tyr)
                msg += '   %8s   %-13s  %-13s  %-13s   %-13s  %-13s  %-13s   %-13s  %s\n' % ('', mxi, myi, mxyi, bmxi, bmyi, bmxyi, txi, tyi)
                f.write(msg)
            msg = page_stamp % page_num
            f.write(msg)
            page_num += 1
        return page_num -1




class ComplexSolidPressureForce(ScalarObject):  # 76-CHEXA_PR,77-PENTA_PR,78-TETRA_PR
    def __init__(self, data_code, is_sort1, isubcase, dt):
        ScalarObject.__init__(self, data_code, isubcase)
        self.acceleration = {}
        self.velocity = {}
        self.pressure = {}

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
            ntimes = len(self.acceleration)
            time0 = get_key0(self.acceleration)
            nelements = len(self.acceleration[time0])
            msg.append('  type=%s ntimes=%s nelements=%s\n'
                       % (self.__class__.__name__, ntimes, nelements))
        else:
            nelements = len(self.acceleration)
            msg.append('  type=%s nelements=%s\n' % (self.__class__.__name__,
                                                     nelements))
        msg.append('  acceleration, velocity, pressure\n')
        return msg

    def add_new_transient(self, dt):
        self.dt = dt
        self.acceleration[dt] = {}
        self.velocity[dt] = {}
        self.pressure[dt] = {}

    def add(self, dt, data):
        [eid, ename, ax, ay, az, vx, vy, vz, pressure] = data
        self.acceleration[eid] = [ax, ay, az]
        self.velocity[eid] = [vx, vy, vz]
        self.pressure[eid] = pressure

    def add_sort1(self, dt, data):
        [eid, ename, ax, ay, az, vx, vy, vz, pressure] = data
        if dt not in self.acceleration:
            self.add_new_transient(dt)
        self.acceleration[dt][eid] = [ax, ay, az]
        self.velocity[dt][eid] = [vx, vy, vz]
        self.pressure[dt][eid] = pressure

    def add_sort2(self, eid, data):
        [dt, ename, ax, ay, az, vx, vy, vz, pressure] = data
        if dt not in self.acceleration:
            self.add_new_transient(dt)
        self.acceleration[dt][eid] = [ax, ay, az]
        self.velocity[dt][eid] = [vx, vy, vz]
        self.pressure[dt][eid] = pressure
