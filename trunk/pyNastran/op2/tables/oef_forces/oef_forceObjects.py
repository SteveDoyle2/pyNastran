from __future__ import (nested_scopes, generators, division, absolute_import,
                        print_function, unicode_literals)
from pyNastran.op2.resultObjects.op2_Objects import ScalarObject
from pyNastran.f06.f06_formatting import writeFloats13E


class RealRodForce(ScalarObject):  # 1-ROD, 3-TUBE, 10-CONROD
    def __init__(self, data_code, is_sort1, isubcase, dt, read_mode):
        ScalarObject.__init__(self, data_code, isubcase, read_mode)
        #self.eType = {}
        self.axialForce = {}
        self.torque = {}

        self.dt = dt
        if is_sort1:
            if dt is not None:
                self.add = self.add_sort1
        else:
            assert dt is not None
            self.add = self.addSort2

    def add_new_transient(self, dt):
        self.dt = dt
        self.axialForce[dt] = {}
        self.torque[dt] = {}

    def get_stats(self):
        msg = self.get_data_code()
        if self.dt is not None:  # transient
            ntimes = len(self.torque)
            time0 = self.torque.keys()[0]
            nelements = len(self.torque[time0])
            msg.append('  type=%s ntimes=%s nelements=%s\n'
                       % (self.__class__.__name__, ntimes, nelements))
        else:
            nelements = len(self.torque)
            msg.append('  type=%s nelements=%s\n' % (self.__class__.__name__,
                                                     nelements))
        msg.append('  axialForce, torque\n')
        return msg

    def add(self, dt, data):
        [eid, axialForce, torque] = data

        #self.eType[eid] = eType
        self.axialForce[eid] = axialForce
        self.torque[eid] = torque

    def add_sort1(self, dt, data):
        [eid, axialForce, torque] = data
        if dt not in self.axialForce:
            self.add_new_transient(dt)

        #self.eType[eid] = eType
        self.axialForce[dt][eid] = axialForce
        self.torque[dt][eid] = torque

    def addSort2(self, eid, data):
        [dt, axialForce, torque] = data
        if dt not in self.axialForce:
            self.add_new_transient(dt)

        #self.eType[eid] = eType
        self.axialForce[dt][eid] = axialForce
        self.torque[dt][eid] = torque

    def __repr__(self):
        return str(self.axialForce)


class RealCBeamForce(ScalarObject):  # 2-CBEAM
    def __init__(self, data_code, is_sort1, isubcase, dt, read_mode):
        ScalarObject.__init__(self, data_code, isubcase, read_mode)
        #self.eType = {}
        self.nodes = {}
        self.bendingMoment = {}
        self.shear = {}
        self.axial = {}
        self.totalTorque = {}
        self.warpingTorque = {}

        self.dt = dt
        if is_sort1:
            if dt is not None:
                self.addNewElement = self.addNewElementSort1
                self.add = self.add_sort1
        else:
            assert dt is not None
            self.addNewElement = self.addNewElementSort2
            self.add = self.addSort2

    def get_stats(self):
        msg = self.get_data_code()
        if self.dt is not None:  # transient
            ntimes = len(self.shear)
            time0 = self.shear.keys()[0]
            nelements = len(self.shear[time0])
            msg.append('  type=%s ntimes=%s nelements=%s\n'
                       % (self.__class__.__name__, ntimes, nelements))
        else:
            nelements = len(self.shear)
            msg.append('  type=%s nelements=%s\n' % (self.__class__.__name__,
                                                     nelements))
        msg.append('  nodes, bendingMoment, shear, axial, totalTorque, '
                   'warpingTorque\n')
        return msg

    def add_new_transient(self, dt):
        self.dt = dt
        self.bendingMoment[dt] = {}
        self.shear[dt] = {}
        self.axial[dt] = {}
        self.totalTorque[dt] = {}
        self.warpingTorque[dt] = {}

    def addNewElement(self, dt, data):
        [eid, nid, sd, bm1, bm2, ts1, ts2, af, ttrq, wtrq] = data
        #print "CBEAM addnew",data
        #self.eType[eid] = eType
        self.nodes[eid] = {sd: [nid]}
        self.bendingMoment[eid] = {sd: [bm1, bm2]}
        self.shear[eid] = {sd: [ts1, ts2]}
        self.axial[eid] = {sd: af}
        self.totalTorque[eid] = {sd: ttrq}
        self.warpingTorque[eid] = {sd: wtrq}

    def add(self, dt, data):
        [eid, nid, sd, bm1, bm2, ts1, ts2, af, ttrq, wtrq] = data
        #print "CBEAM add   ",data

        #self.eType[eid] = eType
        self.nodes[eid][sd] = nid
        self.bendingMoment[eid][sd] = [bm1, bm2]
        self.shear[eid][sd] = [ts1, ts2]
        self.axial[eid][sd] = af
        self.totalTorque[eid][sd] = ttrq
        self.warpingTorque[eid][sd] = wtrq

    def addNewElementSort1(self, dt, data):
        [eid, nid, sd, bm1, bm2, ts1, ts2, af, ttrq, wtrq] = data
        self._fillObjectNew(
            dt, eid, nid, sd, bm1, bm2, ts1, ts2, af, ttrq, wtrq)

    def add_sort1(self, dt, data):
        [eid, nid, sd, bm1, bm2, ts1, ts2, af, ttrq, wtrq] = data
        self._fillObject(dt, eid, nid, sd, bm1, bm2, ts1, ts2, af, ttrq, wtrq)

    def addNewElementSort2(self, eid, data):
        [dt, nid, sd, bm1, bm2, ts1, ts2, af, ttrq, wtrq] = data
        self._fillObjectNew(
            dt, eid, nid, sd, bm1, bm2, ts1, ts2, af, ttrq, wtrq)

    def addSort2(self, eid, data):
        [dt, nid, sd, bm1, bm2, ts1, ts2, af, ttrq, wtrq] = data
        self._fillObject(dt, eid, nid, sd, bm1, bm2, ts1, ts2, af, ttrq, wtrq)

    def _fillObject(self, dt, eid, nid, sd, bm1, bm2, ts1, ts2, af, ttrq, wtrq):
        #if dt not in self.axial:
            #self.add_new_transient(dt)
        #self.eType[eid] = eType
        self.nodes[eid][sd] = nid
        self.bendingMoment[dt][eid][sd] = [bm1, bm2]
        self.shear[dt][eid][sd] = [ts1, ts2]
        self.axial[dt][eid][sd] = af
        self.totalTorque[dt][eid][sd] = ttrq
        self.warpingTorque[dt][eid][sd] = wtrq

    def _fillObjectNew(self, dt, eid, nid, sd, bm1, bm2, ts1, ts2, af, ttrq, wtrq):
        if dt not in self.axial:
            self.add_new_transient(dt)
        #self.eType[eid] = eType
        self.nodes[eid] = {sd: nid}
        self.bendingMoment[dt][eid] = {sd: [bm1, bm2]}
        self.shear[dt][eid] = {sd: [ts1, ts2]}
        self.axial[dt][eid] = {sd: af}
        self.totalTorque[dt][eid] = {sd: ttrq}
        self.warpingTorque[dt][eid] = {sd: wtrq}

    def _write_f06_transient(self, header, pageStamp, f, pageNum=1, is_mag_phase=False):
        words = ['                                 F O R C E S   I N   B E A M   E L E M E N T S        ( C B E A M )\n',
                 '                    STAT DIST/   - BENDING MOMENTS -            - WEB  SHEARS -           AXIAL          TOTAL          WARPING\n',
                 '   ELEMENT-ID  GRID   LENGTH    PLANE 1       PLANE 2        PLANE 1       PLANE 2        FORCE          TORQUE         TORQUE\n']

        msg = []
        for dt, bms in sorted(self.bendingMoment.iteritems()):
            header[1] = ' %s = %10.4E\n' % (self.data_code['name'], dt)
            msg += header + words
            for eid, bm in sorted(bms.iteritems()):
                for sd in sorted(bm):
                    nid = self.nodes[eid][sd]
                    bm1, bm2 = self.bendingMoment[dt][eid][sd]
                    ts1, ts2 = self.shear[dt][eid][sd]
                    af = self.axial[dt][eid][sd]
                    ttrq = self.totalTorque[dt][eid][sd]
                    wtrq = self.warpingTorque[dt][eid][sd]
                    (vals2, isAllZeros) = writeFloats13E([bm1, bm2, ts1, ts2, af, ttrq, wtrq])
                    [bm1, bm2, ts1, ts2, af, ttrq, wtrq] = vals2

                    if sd == 0.:
                        msg.append('0  %8i\n' % (eid))

                    # TODO store grid ID
                    msg.append('           %8i   %.3f   %13s %13s  %13s %13s  %13s  %13s  %-s\n' % (nid, sd, bm1, bm2, ts1, ts2, af, ttrq, wtrq))
            msg.append(pageStamp + str(pageNum) + '\n')
            f.write(''.join(msg))
            msg = ['']
            pageNum += 1
        return pageNum - 1

    def write_f06(self, header, pageStamp, f, pageNum=1, is_mag_phase=False):
        if self.nonlinear_factor is not None:
            return self._write_f06_transient(header, pageStamp, f, pageNum)
        msg = header + ['                                 F O R C E S   I N   B E A M   E L E M E N T S        ( C B E A M )\n',
                        '                    STAT DIST/   - BENDING MOMENTS -            - WEB  SHEARS -           AXIAL          TOTAL          WARPING\n',
                        '   ELEMENT-ID  GRID   LENGTH    PLANE 1       PLANE 2        PLANE 1       PLANE 2        FORCE          TORQUE         TORQUE\n']
        for eid, bm in sorted(self.bendingMoment.iteritems()):
            for sd in sorted(bm):
                bm1, bm2 = self.bendingMoment[eid][sd]
                ts1, ts2 = self.shear[eid][sd]
                af = self.axial[eid][sd]
                ttrq = self.totalTorque[eid][sd]
                wtrq = self.warpingTorque[eid][sd]
                (vals2, isAllZeros) = writeFloats13E([bm1, bm2, ts1, ts2, af, ttrq, wtrq])
                [bm1, bm2, ts1, ts2, af, ttrq, wtrq] = vals2
                msg.append('0  %8i\n' % (eid))
                msg.append('           %8i   %.3f   %13s %13s  %13s %13s  %13s  %13s  %-s\n' % (eid, sd, bm1, bm2, ts1, ts2, af, ttrq, wtrq))

        msg.append(pageStamp + str(pageNum) + '\n')
        f.write(''.join(msg))
        return pageNum

    def __repr__(self):
        #return str(self.axial)
        return self.write_f06(['', ''], 'PAGE', 1)


class RealCShearForce(ScalarObject):  # 4-CSHEAR
    def __init__(self, data_code, is_sort1, isubcase, dt, read_mode):
        ScalarObject.__init__(self, data_code, isubcase, read_mode)
        #self.eType = {}
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
            self.add = self.addSort2

    def get_stats(self):
        msg = self.get_data_code()
        if self.dt is not None:  # transient
            ntimes = len(self.shear12)
            time0 = self.shear12.keys()[0]
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
        self.dt = dt
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
        self._fillObject(dt, eid, f41, f21, f12, f32, f23, f43, f34, f14,
                         kf1, s12, kf2, s23, kf3, s34, kf4, s41)

    def addSort2(self, eid, data):
        [dt, f41, f21, f12, f32, f23, f43, f34, f14,
            kf1, s12, kf2, s23, kf3, s34, kf4, s41] = data

        self._fillObject(dt, eid, f41, f21, f12, f32, f23, f43, f34, f14,
                         kf1, s12, kf2, s23, kf3, s34, kf4, s41)

    def _fillObject(self, dt, eid, f41, f21, f12, f32, f23, f43, f34, f14,
                    kf1, s12, kf2, s23, kf3, s34, kf4, s41):
        if dt not in self.force41:
            self.add_new_transient(dt)
        #self.eType[eid] = eType
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

    def __repr__(self):
        return str(self.force41)


class RealSpringForce(ScalarObject):  # 11-CELAS1,12-CELAS2,13-CELAS3, 14-CELAS4
    def __init__(self, data_code, is_sort1, isubcase, dt, read_mode):
        ScalarObject.__init__(self, data_code, isubcase, read_mode)
        #self.eType = {}
        self.force = {}

        self.dt = dt
        if is_sort1:
            if dt is not None:
                self.add = self.add_sort1
        else:
            assert dt is not None
            self.add = self.addSort2

    def get_stats(self):
        msg = self.get_data_code()
        if self.dt is not None:  # transient
            ntimes = len(self.force)
            time0 = self.force.keys()[0]
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

        #self.eType[eid] = eType
        self.force[eid] = force

    def add_sort1(self, dt, data):
        [eid, force] = data
        if dt not in self.force:
            self.add_new_transient(dt)

        #self.eType[eid] = eType
        self.force[dt][eid] = force

    def addSort2(self, eid, data):
        [dt, force] = data
        if dt not in self.force:
            self.add_new_transient(dt)

        #self.eType[eid] = eType
        self.force[dt][eid] = force

    def _write_f06_transient(self, header, pageStamp, f, pageNum=1, is_mag_phase=False):
        words = ['                              F O R C E S   I N   S C A L A R   S P R I N G S        ( C E L A S 2 )\n',
                 ' \n',
                 '        TIME          FORCE              TIME          FORCE              TIME          FORCE              TIME          FORCE\n']
        msg = []
        for dt, Force in sorted(self.force.iteritems()):
            header[1] = ' %s = %10.4E\n' % (self.data_code['name'], dt)
            msg += header + words

            #packs = []
            forces = []
            elements = []
            line = '   '
            for eid, force in sorted(Force.iteritems()):
                elements.append(eid)
                forces.append(force)
                #pack.append(eid)
                #pack.append(f)
                line += '%13s  %13s     ' % (eid, f)
                if len(forces) == 3:
                    msg.append(line.rstrip() + '\n')

            if forces:
                msg.append(line.rstrip() + '\n')
            msg.append(pageStamp + str(pageNum) + '\n')
            f.write(''.join(msg))
            msg = ['']
            pageNum += 1

        return pageNum - 1

    def write_f06(self, header, pageStamp, f, pageNum=1, is_mag_phase=False):
        if self.nonlinear_factor is not None:
            return self._write_f06_transient(header, pageStamp, f, pageNum)
        msg = header + ['                              F O R C E S   I N   S C A L A R   S P R I N G S        ( C E L A S 2 )\n',
                        ' \n',
                        '        TIME          FORCE              TIME          FORCE              TIME          FORCE              TIME          FORCE\n']
        #packs = []
        forces = []
        elements = []
        line = '   '
        for eid, force in sorted(self.force.iteritems()):
            elements.append(eid)
            forces.append(force)
            #pack.append(eid)
            #pack.append(f)
            line += '%13s  %13s     ' % (eid, force)
            if len(forces) == 3:
                msg.append(line.rstrip() + '\n')

        if forces:
            msg.append(line.rstrip() + '\n')
        msg.append(pageStamp + str(pageNum) + '\n')

        f.write(''.join(msg))
        return pageNum

    def __repr__(self):
        return str(self.force)


class RealDamperForce(ScalarObject):  # 20-CDAMP1,21-CDAMP2,22-CDAMP3,23-CDAMP4
    def __init__(self, data_code, is_sort1, isubcase, dt, read_mode):
        ScalarObject.__init__(self, data_code, isubcase, read_mode)
        #self.eType = {}
        self.force = {}

        self.dt = dt
        if is_sort1:
            if dt is not None:
                self.add = self.add_sort1
        else:
            assert dt is not None
            self.add = self.addSort2

    def get_stats(self):
        msg = self.get_data_code()
        if self.dt is not None:  # transient
            ntimes = len(self.force)
            time0 = self.force.keys()[0]
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

        #self.eType[eid] = eType
        self.force[eid] = force

    def add_sort1(self, dt, data):
        [eid, force] = data
        if dt not in self.force:
            self.add_new_transient(dt)

        #self.eType[eid] = eType
        self.force[dt][eid] = force

    def addSort2(self, eid, data):
        [dt, force] = data
        if dt not in self.force:
            self.add_new_transient(dt)

        #self.eType[eid] = eType
        self.force[dt][eid] = force

    def _write_f06_transient(self, header, pageStamp, f, pageNum=1, is_mag_phase=False):
        words = ['                              F O R C E S   I N   S C A L A R   S P R I N G S        ( C E L A S 2 )\n',
                 ' \n',
                 '        TIME          FORCE              TIME          FORCE              TIME          FORCE              TIME          FORCE\n']
        msg = []
        for dt, Force in sorted(self.force.iteritems()):
            header[1] = ' %s = %10.4E\n' % (self.data_code['name'], dt)
            msg += header + words

            #packs = []
            forces = []
            elements = []
            line = '   '
            for eid, force in sorted(Force.iteritems()):
                elements.append(eid)
                forces.append(force)
                #pack.append(eid)
                #pack.append(f)
                line += '%13s  %13s     ' % (eid, f)
                if len(forces) == 3:
                    msg.append(line.rstrip() + '\n')

            if forces:
                msg.append(line.rstrip() + '\n')
            msg.append(pageStamp + str(pageNum) + '\n')
            f.write(''.join(msg))
            msg = ['']
            pageNum += 1

        return pageNum - 1

    def write_f06(self, header, pageStamp, f, pageNum=1, is_mag_phase=False):
        if self.nonlinear_factor is not None:
            return self._write_f06_transient(header, pageStamp, f, pageNum)
        msg = header + ['                              F O R C E S   I N   S C A L A R   S P R I N G S        ( C E L A S 2 )\n',
                        ' \n',
                        '        TIME          FORCE              TIME          FORCE              TIME          FORCE              TIME          FORCE\n']
        #packs = []
        forces = []
        elements = []
        line = '   '
        for eid, force in sorted(self.force.iteritems()):
            elements.append(eid)
            forces.append(force)
            #pack.append(eid)
            #pack.append(f)
            line += '%13s  %13s     ' % (eid, force)
            if len(forces) == 3:
                msg.append(line.rstrip() + '\n')

        if forces:
            msg.append(line.rstrip() + '\n')
        msg.append(pageStamp + str(pageNum) + '\n')

        f.write(''.join(msg))
        return pageNum

    def __repr__(self):
        return str(self.force)


class RealViscForce(ScalarObject):  # 24-CVISC
    def __init__(self, data_code, is_sort1, isubcase, dt, read_mode):
        ScalarObject.__init__(self, data_code, isubcase, read_mode)
        #self.eType = {}
        self.axialForce = {}
        self.torque = {}

        self.dt = dt
        if is_sort1:
            if dt is not None:
                self.add = self.add_sort1
        else:
            assert dt is not None
            self.add = self.addSort2

    def get_stats(self):
        msg = self.get_data_code()
        if self.dt is not None:  # transient
            ntimes = len(self.torque)
            time0 = self.torque.keys()[0]
            nelements = len(self.torque[time0])
            msg.append('  type=%s ntimes=%s nelements=%s\n'
                       % (self.__class__.__name__, ntimes, nelements))
        else:
            nelements = len(self.torque)
            msg.append('  type=%s nelements=%s\n' % (self.__class__.__name__,
                                                     nelements))
        msg.append('  axialForce, torque\n')
        return msg

    def add_new_transient(self, dt):
        self.dt = dt
        self.axialForce[dt] = {}
        self.torque[dt] = {}

    def add(self, dt, data):
        [eid, axialForce, torque] = data

        #self.eType[eid] = eType
        self.axialForce[eid] = axialForce
        self.torque[eid] = torque

    def add_sort1(self, dt, data):
        [eid, axialForce, torque] = data
        if dt not in self.axialForce:
            self.add_new_transient(dt)

        #self.eType[eid] = eType
        self.axialForce[dt][eid] = axialForce
        self.torque[dt][eid] = torque

    def addSort2(self, eid, data):
        [dt, axialForce, torque] = data
        if dt not in self.axialForce:
            self.add_new_transient(dt)

        #self.eType[eid] = eType
        self.axialForce[dt][eid] = axialForce
        self.torque[dt][eid] = torque

    def __repr__(self):
        return str(self.axialForce)


class RealPlateForce(ScalarObject):  # 33-CQUAD4, 74-CTRIA3
    def __init__(self, data_code, is_sort1, isubcase, dt, read_mode):
        ScalarObject.__init__(self, data_code, isubcase, read_mode)
        #self.eType = {}
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
            self.add = self.addSort2

    def get_stats(self):
        msg = self.get_data_code()
        if self.dt is not None:  # transient
            ntimes = len(self.mx)
            time0 = self.mx.keys()[0]
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

    def add(self, dt, data):
        [eid, mx, my, mxy, bmx, bmy, bmxy, tx, ty] = data

        #self.eType[eid] = eType
        self.mx[eid] = mx
        self.my[eid] = my
        self.mxy[eid] = mxy
        self.bmx[eid] = bmx
        self.bmy[eid] = bmy
        self.bmxy[eid] = bmxy
        self.tx[eid] = tx
        self.ty[eid] = ty

    def add_sort1(self, dt, data):
        [eid, mx, my, mxy, bmx, bmy, bmxy, tx, ty] = data
        if dt not in self.mx:
            self.add_new_transient(dt)

        #self.eType[eid] = eType
        self.mx[dt][eid] = mx
        self.my[dt][eid] = my
        self.mxy[dt][eid] = mxy
        self.bmx[dt][eid] = bmx
        self.bmy[dt][eid] = bmy
        self.bmxy[dt][eid] = bmxy
        self.tx[dt][eid] = tx
        self.ty[dt][eid] = ty

    def addSort2(self, eid, data):
        [dt, mx, my, mxy, bmx, bmy, bmxy, tx, ty] = data
        if dt not in self.mx:
            self.add_new_transient(dt)

        #self.eType[eid] = eType
        self.mx[dt][eid] = mx
        self.my[dt][eid] = my
        self.mxy[dt][eid] = mxy
        self.bmx[dt][eid] = bmx
        self.bmy[dt][eid] = bmy
        self.bmxy[dt][eid] = bmxy
        self.tx[dt][eid] = tx
        self.ty[dt][eid] = ty

    def __repr__(self):
        return str(self.mx)


class RealPlate2Force(ScalarObject):  # 64-CQUAD8, 75-CTRIA6, 82-CQUADR
    def __init__(self, data_code, is_sort1, isubcase, dt, read_mode):
        ScalarObject.__init__(self, data_code, isubcase, read_mode)
        #self.eType = {}
        self.term = {}
        self.ngrids = {}
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
                self.addNewElement = self.addNewElementSort1
                self.add = self.add_sort1
        else:
            assert dt is not None
            self.addNewElement = self.addNewElementSort2
            self.add = self.addSort2

    def get_stats(self):
        msg = self.get_data_code()
        if self.dt is not None:  # transient
            ntimes = len(self.mx)
            if ntimes == 0:
                nelements = 0
            else:
                time0 = self.mx.keys()[0]
                nelements = len(self.mx[time0])
            msg.append('  type=%s ntimes=%s nelements=%s\n'
                       % (self.__class__.__name__, ntimes, nelements))
        else:
            nelements = len(self.mx)
            msg.append('  type=%s nelements=%s\n' % (self.__class__.__name__,
                                                     nelements))
        msg.append('  term, ngrids, mx, my, mxy, bmx, bmy, bmxy, tx, ty\n')
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

    def addNewElement(self, eid, dt, data):
        #print "eid = ",eid
        [term, nid, mx, my, mxy, bmx, bmy, bmxy, tx, ty] = data

        #self.eType[eid] = eType
        self.term[eid] = term
        self.ngrids[eid] = nid

        self.mx[eid] = [mx]
        self.my[eid] = [my]
        self.mxy[eid] = [mxy]
        self.bmx[eid] = [bmx]
        self.bmy[eid] = [bmy]
        self.bmxy[eid] = [bmxy]
        self.tx[eid] = [tx]
        self.ty[eid] = [ty]

    def add(self, eid, dt, data):
        [nid, mx, my, mxy, bmx, bmy, bmxy, tx, ty] = data

        #self.eType[eid] = eType
        #print "mx = ",self.mx,mx
        self.mx[eid].append(mx)
        self.my[eid].append(my)
        self.mxy[eid].append(mxy)
        self.bmx[eid].append(bmx)
        self.bmy[eid].append(bmy)
        self.bmxy[eid].append(bmxy)
        self.tx[eid].append(tx)
        self.ty[eid].append(ty)

    def addNewElementSort1(self, eid, dt, data):
        [term, nid, mx, my, mxy, bmx, bmy, bmxy, tx, ty] = data
        if dt not in self.mx:
            self.add_new_transient(dt)

        #self.eType[eid] = eType
        self.term[eid] = term
        self.ngrids[eid] = nid
        self.mx[dt][eid] = [mx]
        self.my[dt][eid] = [my]
        self.mxy[dt][eid] = [mxy]
        self.bmx[dt][eid] = [bmx]
        self.bmy[dt][eid] = [bmy]
        self.bmxy[dt][eid] = [bmxy]
        self.tx[dt][eid] = [tx]
        self.ty[dt][eid] = [ty]

    def add_sort1(self, eid, dt, data):
        [nid, mx, my, mxy, bmx, bmy, bmxy, tx, ty] = data
        if dt not in self.mx:
            self.add_new_transient(dt)

        #self.eType[eid] = eType
        self.mx[dt][eid].append(mx)
        self.my[dt][eid].append(my)
        self.mxy[dt][eid].append(mxy)
        self.bmx[dt][eid].append(bmx)
        self.bmy[dt][eid].append(bmy)
        self.bmxy[dt][eid].append(bmxy)
        self.tx[dt][eid].append(tx)
        self.ty[dt][eid].append(ty)

    def addNewElementSort2(self, dt, eid, data):
        [term, nid, mx, my, mxy, bmx, bmy, bmxy, tx, ty] = data
        if dt not in self.mx:
            self.add_new_transient(dt)

        #self.eType[eid] = eType
        self.term[eid] = term
        self.ngrids[eid] = nid

        self.mx[dt][eid] = [mx]
        self.my[dt][eid] = [my]
        self.mxy[dt][eid] = [mxy]
        self.bmx[dt][eid] = [bmx]
        self.bmy[dt][eid] = [bmy]
        self.bmxy[dt][eid] = [bmxy]
        self.tx[dt][eid] = [tx]
        self.ty[dt][eid] = [ty]

    def addSort2(self, dt, eid, data):
        [nid, mx, my, mxy, bmx, bmy, bmxy, tx, ty] = data
        if dt not in self.mx:
            self.add_new_transient(dt)

        #self.eType[eid] = eType
        self.mx[dt][eid].append(mx)
        self.my[dt][eid].append(my)
        self.mxy[dt][eid].append(mxy)
        self.bmx[dt][eid].append(bmx)
        self.bmy[dt][eid].append(bmy)
        self.bmxy[dt][eid].append(bmxy)
        self.tx[dt][eid].append(tx)
        self.ty[dt][eid].append(ty)

    def __repr__(self):
        return str(self.mx)


class RealCBarForce(ScalarObject):  # 34-CBAR
    def __init__(self, data_code, is_sort1, isubcase, dt, read_mode):
        ScalarObject.__init__(self, data_code, isubcase, read_mode)
        #self.eType = {}
        self.shape = {}
        self._inode_start = None
        self._inode_end = None
        self.data = None

        self.bendingMomentA = {}
        self.bendingMomentB = {}
        self.shear = {}
        self.axial = {}
        self.torque = {}

        self.dt = dt
        if is_sort1:
            if dt is not None:
                self.add = self.add_sort1
        else:
            assert dt is not None
            self.add = self.addSort2

    def _increase_size(self, dt, nelements):
        #self.shape += 1
        if dt in self.shape:  # default dictionary
            self.shape[dt] += nelements
        else:
            self.shape[dt] = nelements
            #self.n_nonlinear += 1

    def _get_shape(self):
        ndt = len(self.shape)
        dts = self.shape.keys()
        shape0 = dts[0]
        nelements = self.shape[shape0]
        #print "ndt=%s nelements=%s dts=%s" % (ndt, nelements, str(dts))
        return ndt, nelements, dts

    def get_stats(self):
        msg = self.get_data_code()
        ndt, nelements, dts = self._get_shape()
        if self.dt is not None:  # transient
            ntimes = len(self.torque)
            name = self.data_code['name']
            dt_string = name + ', '
            time0 = self.torque.keys()[0]
            nelements = len(self.torque[time0])
            msg.append('  type=%s n%ss=%s nelements=%s\n'
                       % (self.__class__.__name__, name, ndt, nelements))
        else:
            nelements = len(self.torque)
            dt_string = ''
            msg.append('  real type=%s nelements=%s\n' % (self.__class__.__name__,
                                                     nelements))
        msg.append('  data: index  : %selement_id\n' % dt_string)
        msg.append('  data: results: bendingMomentA, bendingMomentB, shear, axial, torque\n')
        return msg

    def add_new_transient(self, dt):
        self.dt = dt
        self.bendingMomentA[dt] = {}
        self.bendingMomentB[dt] = {}
        self.shear[dt] = {}
        self.axial[dt] = {}
        self.torque[dt] = {}

    def add(self, dt, data):
        [eid, bm1a, bm2a, bm1b, bm2b, ts1, ts2, af, trq] = data

        #self.eType[eid] = eType
        self.bendingMomentA[eid] = [bm1a, bm2a]
        self.bendingMomentB[eid] = [bm1b, bm2b]
        self.shear[eid] = [ts1, ts2]
        self.axial[eid] = af
        self.torque[eid] = trq

    def add_sort1(self, dt, data):
        [eid, bm1a, bm2a, bm1b, bm2b, ts1, ts2, af, trq] = data
        if dt not in self.axial:
            self.add_new_transient(dt)

        #self.eType[eid] = eType
        self.bendingMomentA[dt][eid] = [bm1a, bm2a]
        self.bendingMomentB[dt][eid] = [bm1b, bm2b]
        self.shear[dt][eid] = [ts1, ts2]
        self.axial[dt][eid] = af
        self.torque[dt][eid] = trq

    def addSort2(self, eid, data):
        [dt, bm1a, bm2a, bm1b, bm2b, ts1, ts2, af, trq] = data
        if dt not in self.axial:
            self.add_new_transient(dt)

        #self.eType[eid] = eType
        self.bendingMomentA[dt][eid] = [bm1a, bm2a]
        self.bendingMomentB[dt][eid] = [bm1b, bm2b]
        self.shear[dt][eid] = [ts1, ts2]
        self.axial[dt][eid] = af
        self.torque[dt][eid] = trq

    def _write_f06_transient(self, header, pageStamp, f, pageNum=1, is_mag_phase=False):
        words = ['                                 F O R C E S   I N   B A R   E L E M E N T S         ( C B A R )\n',
                 '0    ELEMENT         BEND-MOMENT END-A            BEND-MOMENT END-B                - SHEAR -               AXIAL\n',
                 '       ID.         PLANE 1       PLANE 2        PLANE 1       PLANE 2        PLANE 1       PLANE 2         FORCE         TORQUE\n']
        msg = []
        for dt, bm in sorted(self.bendingMomentA.iteritems()):
            header[1] = ' %s = %10.4E\n' % (self.data_code['name'], dt)
            msg += header + words
            for eid in sorted(bm):
                bm1a, bm2a = self.bendingMomentA[dt][eid]
                bm1b, bm2b = self.bendingMomentB[dt][eid]
                ts1, ts2 = self.shear[dt][eid]
                af = self.axial[dt][eid]
                trq = self.torque[dt][eid]
                (vals2, isAllZeros) = writeFloats13E([bm1a, bm2a, bm1b, bm2b, ts1, ts2, af, trq])
                [bm1a, bm2a, bm1b, bm2b, ts1, ts2, af, trq] = vals2
                msg.append('      %8i    %13s %13s  %13s %13s  %13s %13s  %13s  %-s\n' % (eid, bm1a, bm2a, bm1b, bm2b, ts1, ts2, af, trq))
#            1     2.504029E+06  9.728743E+06   5.088001E+05  1.976808E+06   1.995229E+06  7.751935E+06  -3.684978E-07  -1.180941E-07

            msg.append(pageStamp + str(pageNum) + '\n')
            f.write(''.join(msg))
            msg = ['']
            pageNum += 1
        return pageNum - 1

    def __repr__(self):
        return self._get_stats()


class RealCBar100Force(ScalarObject):  # 100-CBAR
    def __init__(self, data_code, is_sort1, isubcase, dt, read_mode):
        ScalarObject.__init__(self, data_code, isubcase, read_mode)
        #self.eType = {}
        self.bendingMoment = {}
        self.shear = {}
        self.axial = {}
        self.torque = {}

        self.dt = dt
        if is_sort1:
            if dt is not None:
                self.add = self.add_sort1
        else:
            assert dt is not None
            self.add = self.addSort2

    def get_stats(self):
        msg = self.get_data_code()
        if self.dt is not None:  # transient
            ntimes = len(self.torque)
            time0 = self.torque.keys()[0]
            nelements = len(self.torque[time0])
            msg.append('  type=%s ntimes=%s nelements=%s\n'
                       % (self.__class__.__name__, ntimes, nelements))
        else:
            nelements = len(self.torque)
            msg.append('  type=%s nelements=%s\n' % (self.__class__.__name__,
                                                     nelements))
        msg.append('  bendingMoment, shear, axial, torque\n')
        return msg

    def add_new_transient(self, dt):
        self.dt = dt
        self.bendingMoment[dt] = {}
        self.shear[dt] = {}
        self.axial[dt] = {}
        self.torque[dt] = {}

    def add(self, dt, data):
        [eid, sd, bm1, bm2, ts1, ts2, af, trq] = data

        #self.eType[eid] = eType
        self.bendingMoment[eid] = [bm1, bm2]
        self.shear[eid] = [ts1, ts2]
        self.axial[eid] = af
        self.torque[eid] = trq

    def add_sort1(self, dt, data):
        [eid, sd, bm1, bm2, ts1, ts2, af, trq] = data
        if dt not in self.axial:
            self.add_new_transient(dt)

        #self.eType[eid] = eType
        self.bendingMoment[dt][eid] = [bm1, bm2]
        self.shear[dt][eid] = [ts1, ts2]
        self.axial[dt][eid] = af
        self.torque[dt][eid] = trq

    def addSort2(self, eid, data):
        [dt, sd, bm1, bm2, ts1, ts2, af, trq] = data
        if dt not in self.axial:
            self.add_new_transient(dt)

        #self.eType[eid] = eType
        self.bendingMoment[dt][eid] = [bm1, bm2]
        self.shear[dt][eid] = [ts1, ts2]
        self.axial[dt][eid] = af
        self.torque[dt][eid] = trq

    def __repr__(self):
        return str(self.axial)


class RealConeAxForce(ScalarObject):  # 35-CCONEAX
    def __init__(self, data_code, is_sort1, isubcase, dt, read_mode):
        ScalarObject.__init__(self, data_code, isubcase, read_mode)
        #self.eType = {}
        self.hopa = {}
        self.bmu = {}
        self.bmv = {}
        self.tm = {}
        self.su = {}
        self.sv = {}

        self.dt = dt
        if is_sort1:
            if dt is not None:
                self.add = self.add_sort1
        else:
            assert dt is not None
            self.add = self.addSort2

    def get_stats(self):
        msg = self.get_data_code()
        if self.dt is not None:  # transient
            ntimes = len(self.hopa)
            time0 = self.hopa.keys()[0]
            nelements = len(self.hopa[time0])
            msg.append('  type=%s ntimes=%s nelements=%s\n'
                       % (self.__class__.__name__, ntimes, nelements))
        else:
            nelements = len(self.hopa)
            msg.append('  type=%s nelements=%s\n' % (self.__class__.__name__,
                                                     nelements))
        msg.append('  hopa, bmu, bmv, tm, su, sv\n')
        return msg

    def add_new_transient(self, dt):
        self.dt = dt
        self.hopa[dt] = {}
        self.bmu[dt] = {}
        self.bmv[dt] = {}
        self.tm[dt] = {}
        self.su[dt] = {}
        self.sv[dt] = {}

    def add(self, dt, data):
        [eid, hopa, bmu, bmv, tm, su, sv] = data

        #self.eType[eid] = eType
        self.hopa[eid] = hopa
        self.bmu[eid] = bmu
        self.bmv[eid] = bmv
        self.tm[eid] = tm
        self.su[eid] = su
        self.sv[eid] = sv

    def add_sort1(self, dt, data):
        [eid, hopa, bmu, bmv, tm, su, sv] = data
        if dt not in self.hopa:
            self.add_new_transient(dt)

        #self.eType[eid] = eType
        self.hopa[dt][eid] = hopa
        self.bmu[dt][eid] = bmu
        self.bmv[dt][eid] = bmv
        self.tm[dt][eid] = tm
        self.su[dt][eid] = su
        self.sv[dt][eid] = sv

    def addSort2(self, eid, data):
        [dt, hopa, bmu, bmv, tm, su, sv] = data
        if dt not in self.hopa:
            self.add_new_transient(dt)

        #self.eType[eid] = eType
        self.hopa[dt][eid] = hopa
        self.bmu[dt][eid] = bmu
        self.bmv[dt][eid] = bmv
        self.tm[dt][eid] = tm
        self.su[dt][eid] = su
        self.sv[dt][eid] = sv

    def __repr__(self):
        return str(self.hopa)


class RealCGapForce(ScalarObject):  # 38-CGAP
    def __init__(self, data_code, is_sort1, isubcase, dt, read_mode):
        ScalarObject.__init__(self, data_code, isubcase, read_mode)
        #self.eType = {}
        self.fx = {}
        self.sfy = {}
        self.sfz = {}
        self.u = {}
        self.v = {}
        self.w = {}
        self.sv = {}
        self.sw = {}

        self.dt = dt
        if is_sort1:
            if dt is not None:
                self.add = self.add_sort1
        else:
            assert dt is not None
            self.add = self.addSort2

    def get_stats(self):
        msg = self.get_data_code()
        if self.dt is not None:  # transient
            ntimes = len(self.fx)
            time0 = self.fx.keys()[0]
            nelements = len(self.fx[time0])
            msg.append('  type=%s ntimes=%s nelements=%s\n'
                       % (self.__class__.__name__, ntimes, nelements))
        else:
            nelements = len(self.fx)
            msg.append('  type=%s nelements=%s\n' % (self.__class__.__name__,
                                                     nelements))
        msg.append('  fx, sfy, sfz, u, v, w, sv, sw\n')
        return msg

    def add_new_transient(self, dt):
        self.dt = dt
        self.fx[dt] = {}
        self.sfy[dt] = {}
        self.sfz[dt] = {}
        self.u[dt] = {}
        self.v[dt] = {}
        self.w[dt] = {}
        self.sv[dt] = {}
        self.sw[dt] = {}

    def add(self, dt, data):
        [eid, fx, sfy, sfz, u, v, w, sv, sw] = data

        #self.eType[eid] = eType
        self.fx[eid] = fx
        self.sfy[eid] = sfy
        self.sfz[eid] = sfz
        self.u[eid] = u
        self.v[eid] = v
        self.w[eid] = w
        self.sv[eid] = sv
        self.sw[eid] = sw

    def add_sort1(self, dt, data):
        [eid, fx, sfy, sfz, u, v, w, sv, sw] = data
        if dt not in self.fx:
            self.add_new_transient(dt)

        #self.eType[eid] = eType
        self.fx[dt][eid] = fx
        self.sfy[dt][eid] = sfy
        self.sfz[dt][eid] = sfz
        self.u[dt][eid] = u
        self.v[dt][eid] = v
        self.w[dt][eid] = w
        self.sv[dt][eid] = sv
        self.sw[dt][eid] = sw

    def addSort2(self, eid, data):
        [dt, fx, sfy, sfz, u, v, w, sv, sw] = data
        if dt not in self.fx:
            self.add_new_transient(dt)

        #self.eType[eid] = eType
        self.fx[dt][eid] = fx
        self.sfy[dt][eid] = sfy
        self.sfz[dt][eid] = sfz
        self.u[dt][eid] = u
        self.v[dt][eid] = v
        self.w[dt][eid] = w
        self.sv[dt][eid] = sv
        self.sw[dt][eid] = sw

    def __repr__(self):
        return str(self.fx)


class RealBendForce(ScalarObject):  # 69-CBEND
    def __init__(self, data_code, is_sort1, isubcase, dt, read_mode):
        ScalarObject.__init__(self, data_code, isubcase, read_mode)
        #self.eType = {}
        self.nodeIDs = {}
        self.bendingMoment1 = {}
        self.bendingMoment2 = {}
        self.shearPlane1 = {}
        self.shearPlane2 = {}
        self.axial = {}
        self.torque = {}

        self.dt = dt
        if is_sort1:
            if dt is not None:
                self.add = self.add_sort1
        else:
            assert dt is not None
            self.add = self.addSort2

    def get_stats(self):
        msg = self.get_data_code()
        if self.dt is not None:  # transient
            ntimes = len(self.torque)
            time0 = self.torque.keys()[0]
            nelements = len(self.torque[time0])
            msg.append('  type=%s ntimes=%s nelements=%s\n'
                       % (self.__class__.__name__, ntimes, nelements))
        else:
            nelements = len(self.torque)
            msg.append('  type=%s nelements=%s\n' % (self.__class__.__name__,
                                                     nelements))
        msg.append('  nodeIDs, bendingMoment1, bendingMoment2, '
                   'shearPlate1, shearPlate2, axial, torque\n')
        return msg

    def add_new_transient(self, dt):
        self.dt = dt
        self.bendingMoment1[dt] = {}
        self.bendingMoment2[dt] = {}
        self.shearPlane1[dt] = {}
        self.shearPlane2[dt] = {}
        self.axial[dt] = {}
        self.torque[dt] = {}

    def add(self, dt, data):
        [eid, nidA, bm1A, bm2A, sp1A, sp2A, axialA, torqueA,
         nidB, bm1B, bm2B, sp1B, sp2B, axialB, torqueB] = data

        #self.eType[eid] = eType
        self.nodeIDs[eid] = [nidA, nidB]
        self.bendingMoment1[eid] = [bm1A, bm1B]
        self.bendingMoment2[eid] = [bm2A, bm2B]
        self.shearPlane1[eid] = [sp1A, sp1B]
        self.shearPlane2[eid] = [sp2A, sp2B]
        self.axial[eid] = [axialA, axialB]
        self.torque[eid] = [torqueA, torqueB]

    def add_sort1(self, dt, data):
        [eid, nidA, bm1A, bm2A, sp1A, sp2A, axialA, torqueA,
         nidB, bm1B, bm2B, sp1B, sp2B, axialB, torqueB] = data
        self._fillObject(
            dt, eid, nidA, bm1A, bm2A, sp1A, sp2A, axialA, torqueA,
            nidB, bm1B, bm2B, sp1B, sp2B, axialB, torqueB)

    def addSort2(self, eid, data):
        [dt, nidA, bm1A, bm2A, sp1A, sp2A, axialA, torqueA,
            nidB, bm1B, bm2B, sp1B, sp2B, axialB, torqueB] = data
        self._fillObject(dt, eid, nidA, bm1A, bm2A, sp1A, sp2A, axialA, torqueA,
                                  nidB, bm1B, bm2B, sp1B, sp2B, axialB, torqueB)

    def _fillObject(self, dt, eid, nidA, bm1A, bm2A, sp1A, sp2A, axialA, torqueA,
                                   nidB, bm1B, bm2B, sp1B, sp2B, axialB, torqueB):
        if dt not in self.axial:
            self.add_new_transient(dt)

        #self.eType[eid] = eType
        self.nodeIDs[eid] = [nidA, nidB]
        self.bendingMoment1[dt][eid] = [bm1A, bm1B]
        self.bendingMoment2[dt][eid] = [bm2A, bm2B]
        self.shearPlane1[dt][eid] = [sp1A, sp1B]
        self.shearPlane2[dt][eid] = [sp2A, sp2B]
        self.axial[dt][eid] = [axialA, axialB]
        self.torque[dt][eid] = [torqueA, torqueB]

    def __repr__(self):
        return str(self.axial)


class RealPentaPressureForce(ScalarObject):  # 77-PENTA_PR,78-TETRA_PR
    def __init__(self, data_code, is_sort1, isubcase, dt, read_mode):
        ScalarObject.__init__(self, data_code, isubcase, read_mode)
        #self.eType = {}
        self.acceleration = {}
        self.velocity = {}
        self.pressure = {}

        self.dt = dt
        if is_sort1:
            if dt is not None:
                self.add = self.add_sort1
        else:
            assert dt is not None
            self.add = self.addSort2

    def get_stats(self):
        msg = self.get_data_code()
        if self.dt is not None:  # transient
            ntimes = len(self.acceleration)
            time0 = self.acceleration.keys()[0]
            nelements = len(self.acceleration[time0])
            msg.append('  type=%s ntimes=%s nelements=%s\n'
                       % (self.__class__.__name__, ntimes, nelements))
        else:
            nelements = len(self.force)
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
        [eid, eName, ax, ay, az, vx, vy, vz, pressure] = data

        #self.eType[eid] = eType
        self.acceleration[eid] = [ax, ay, az]
        self.velocity[eid] = [vx, vy, vz]
        self.pressure[eid] = pressure

    def add_sort1(self, dt, data):
        [eid, eName, ax, ay, az, vx, vy, vz, pressure] = data
        if dt not in self.acceleration:
            self.add_new_transient(dt)

        #self.eType[eid] = eType
        self.acceleration[dt][eid] = [ax, ay, az]
        self.velocity[dt][eid] = [vx, vy, vz]
        self.pressure[dt][eid] = pressure

    def addSort2(self, eid, data):
        [dt, eName, ax, ay, az, vx, vy, vz, pressure] = data
        if dt not in self.acceleration:
            self.add_new_transient(dt)

        #self.eType[eid] = eType
        self.acceleration[dt][eid] = [ax, ay, az]
        self.velocity[dt][eid] = [vx, vy, vz]
        self.pressure[dt][eid] = pressure

    def write_f06(self, header, pageStamp, f, pageNum=1, is_mag_phase=False):
        #words = ['                                   P E A K   A C C E L E R A T I O N S   A N D   P R E S S U R E S\n',
        #         ' \n',
        #         '    TIME         EL-TYPE             X-ACCELERATION            Y-ACCELERATION            Z-ACCELERATION            PRESSURE (DB)\n']
        if self.nonlinear_factor is not None:
            return self._write_f06_transient(header, pageStamp, f, pageNum)
        return 'RealPentaPressureForce write_f06 not implemented...\n'
        #raise NotImplementedError()

    def _write_f06_transient(self, header, pageStamp, f, pageNum=1, is_mag_phase=False):
        words = ['                                   P E A K   A C C E L E R A T I O N S   A N D   P R E S S U R E S\n',
                 ' \n',
                 '    TIME         EL-TYPE             X-ACCELERATION            Y-ACCELERATION            Z-ACCELERATION            PRESSURE (DB)\n']
        msg = []
        for dt, acc in sorted(self.acceleration.iteritems()):
            header[1] = ' %s = %10.4E\n' % (self.data_code['name'], dt)
            msg += header + words
            for eid in sorted(acc):
                ax, ay, az = self.acceleration[dt][eid]
                vx, vy, vz = self.velocity[dt][eid]
                pressure = self.pressure[dt][eid]
                vals = [ax, ay, az, pressure]
                (vals2, isAllZeros) = writeFloats13E(vals)
                [ax, ay, az, pressure] = vals2
                eType = 'PENPR'
                msg.append('0%13s    %5s               %13s             %13s             %13s             %-s\n' % (eid, eType, ax, ay, az, pressure))
            msg.append(pageStamp + str(pageNum) + '\n')
            f.write(''.join(msg))
            msg = ['']
            pageNum += 1
        return pageNum - 1

    def __repr__(self):
        return str(self.acceleration)


class RealCBushForce(ScalarObject):  # 102-CBUSH
    def __init__(self, data_code, is_sort1, isubcase, dt, read_mode):
        ScalarObject.__init__(self, data_code, isubcase, read_mode)
        #self.eType = {}
        self.force = {}
        self.moment = {}

        self.dt = dt
        if is_sort1:
            if dt is not None:
                self.add = self.add_sort1
        else:
            assert dt is not None
            self.add = self.addSort2

    def get_stats(self):
        msg = self.get_data_code()
        if self.dt is not None:  # transient
            ntimes = len(self.force)
            time0 = self.force.keys()[0]
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

    def add(self, dt, data):
        [eid, fx, fy, fz, mx, my, mz] = data

        #self.eType[eid] = eType
        self.force[eid] = [fx, fy, fz]
        self.moment[eid] = [mx, my, mz]

    def add_sort1(self, dt, data):
        [eid, fx, fy, fz, mx, my, mz] = data
        if dt not in self.force:
            self.add_new_transient(dt)

        #self.eType[eid] = eType
        self.force[dt][eid] = [fx, fy, fz]
        self.moment[dt][eid] = [mx, my, mz]

    def addSort2(self, eid, data):
        [dt, fx, fy, fz, mx, my, mz] = data
        if dt not in self.force:
            self.add_new_transient(dt)

        #self.eType[eid] = eType
        self.force[dt][eid] = [fx, fy, fz]
        self.moment[dt][eid] = [mx, my, mz]

    def __repr__(self):
        return str(self.force)


class RealForce_VU(ScalarObject):  # 191-VUBEAM
    def __init__(self, data_code, is_sort1, isubcase, dt, read_mode):
        ScalarObject.__init__(self, data_code, isubcase, read_mode)
        #self.eType = {}
        self.parent = {}
        self.coord = {}
        self.icord = {}

        self.forceX = {}
        self.shearY = {}
        self.shearZ = {}
        self.torsion = {}
        self.bendingY = {}
        self.bendingZ = {}

        # TODO if dt=None, handle SORT1 case
        self.dt = dt
        if is_sort1:
            if dt is not None:
                self.add = self.add_sort1
        else:
            assert dt is not None
            self.add = self.addSort2

    def get_stats(self):
        msg = self.get_data_code()
        nelements = len(self.coord)
        if self.dt is not None:  # transient
            ntimes = len(self.forceX)
            msg.append('  type=%s ntimes=%s nelements=%s\n'
                       % (self.__class__.__name__, ntimes, nelements))
        else:
            msg.append('  type=%s nelements=%s\n' % (self.__class__.__name__,
                                                     nelements))
        msg.append('  parent, coord, icord, forceX, shearY, shearZ, torsion, '
                   'bendingY, bendingZ\n')
        return msg

    def add_new_transient(self, dt):
        self.dt = dt
        self.forceX[dt] = {}
        self.shearY[dt] = {}
        self.shearZ[dt] = {}
        self.torsion[dt] = {}
        self.bendingY[dt] = {}
        self.bendingZ[dt] = {}

    def add(self, nNodes, dt, data):
        [eid, parent, coord, icord, forces] = data
        self.parent[eid] = parent
        self.coord[eid] = coord
        self.icord[eid] = icord
        #self.eType[eid]    = eType

        self.forceX[eid] = {}
        self.shearY[eid] = {}
        self.shearZ[eid] = {}
        self.torsion[eid] = {}
        self.bendingY[eid] = {}
        self.bendingZ[eid] = {}

        for force in forces:
            [nid, posit, forceX, shearY, shearZ, torsion,
                bendingY, bendingZ] = force
            self.forceX[eid][nid] = forceX
            self.shearY[eid][nid] = shearY
            self.shearZ[eid][nid] = shearZ
            self.torsion[eid][nid] = torsion
            self.bendingY[eid][nid] = bendingY
            self.bendingZ[eid][nid] = bendingZ

    def add_sort1(self, nNodes, dt, data):
        [eid, parent, coord, icord, forces] = data
        if dt not in self.forceX:
            self.add_new_transient(dt)
        self.parent[eid] = parent
        self.coord[eid] = coord
        self.icord[eid] = icord
        #self.eType[eid]    = eType

        self.forceX[dt][eid] = {}
        self.shearY[dt][eid] = {}
        self.shearZ[dt][eid] = {}
        self.torsion[dt][eid] = {}
        self.bendingY[dt][eid] = {}
        self.bendingZ[dt][eid] = {}

        for force in forces:
            [nid, posit, forceX, shearY, shearZ, torsion,
                bendingY, bendingZ] = force
            self.forceX[dt][eid][nid] = forceX
            self.shearY[dt][eid][nid] = shearY
            self.shearZ[dt][eid][nid] = shearZ
            self.torsion[dt][eid][nid] = torsion
            self.bendingY[dt][eid][nid] = bendingY
            self.bendingZ[dt][eid][nid] = bendingZ

    def addSort2(self, nNodes, eid, data):
        [dt, parent, coord, icord, forces] = data
        if dt not in self.forceX:
            self.add_new_transient(dt)
        self.parent[eid] = parent
        self.coord[eid] = coord
        self.icord[eid] = icord
        #self.eType[eid]    = eType

        self.forceX[dt][eid] = {}
        self.shearY[dt][eid] = {}
        self.shearZ[dt][eid] = {}
        self.torsion[dt][eid] = {}
        self.bendingY[dt][eid] = {}
        self.bendingZ[dt][eid] = {}
        for force in forces:
            [nid, posit, forceX, shearY, shearZ, torsion,
                bendingY, bendingZ] = force
            self.forceX[dt][eid][nid] = forceX
            self.shearY[dt][eid][nid] = shearY
            self.shearZ[dt][eid][nid] = shearZ
            self.torsion[dt][eid][nid] = torsion
            self.bendingY[dt][eid][nid] = bendingY
            self.bendingZ[dt][eid][nid] = bendingZ

    def __repr__(self):
        return str(self.forceX)


class RealForce_VU_2D(ScalarObject):  # 190-VUTRIA # 189-VUQUAD
    def __init__(self, data_code, is_sort1, isubcase, dt, read_mode):
        ScalarObject.__init__(self, data_code, isubcase, read_mode)
        #self.eType = {}
        self.parent = {}
        self.coord = {}
        self.icord = {}
        self.theta = {}

        self.membraneX = {}
        self.membraneY = {}
        self.membraneXY = {}
        self.bendingX = {}
        self.bendingY = {}
        self.bendingXY = {}
        self.shearYZ = {}
        self.shearXZ = {}

        # TODO if dt=None, handle SORT1 case
        self.dt = dt
        if is_sort1:
            if dt is not None:
                self.add = self.add_sort1
        else:
            assert dt is not None
            self.add = self.addSort2

    def get_stats(self):
        msg = self.get_data_code()
        nelements = len(self.coord)
        if self.dt is not None:  # transient
            ntimes = len(self.membraneX)
            msg.append('  type=%s ntimes=%s nelements=%s\n'
                       % (self.__class__.__name__, ntimes, nelements))
        else:
            msg.append('  type=%s nelements=%s\n' % (self.__class__.__name__,
                                                     nelements))
        msg.append('  parent, coord, icord, theta, membraneX, membraneY, '
                   'membraneXY, bendingX, bendingY, bendingXY, '
                   'shearYZ, shearXZ\n')
        return msg

    def add_new_transient(self, dt):
        self.membraneX[dt] = {}
        self.membraneY[dt] = {}
        self.membraneXY[dt] = {}
        self.bendingX[dt] = {}
        self.bendingY[dt] = {}
        self.bendingXY[dt] = {}
        self.shearYZ[dt] = {}
        self.shearXZ[dt] = {}

    def add(self, nNodes, dt, data):
        [eid, parent, coord, icord, theta, forces] = data
        self.parent[eid] = parent
        self.coord[eid] = coord
        self.icord[eid] = icord
        self.theta[eid] = theta
        #self.eType[eid]    = eType

        self.membraneX[eid] = {}
        self.membraneY[eid] = {}
        self.membraneXY[eid] = {}
        self.bendingX[eid] = {}
        self.bendingY[eid] = {}
        self.bendingXY[eid] = {}
        self.shearYZ[eid] = {}
        self.shearXZ[eid] = {}

        for force in forces:
            [nid, membraneX, membraneY, membraneXY, bendingX,
                bendingY, bendingXY, shearYZ, shearXZ] = force
            self.membraneX[eid][nid] = membraneX
            self.membraneY[eid][nid] = membraneY
            self.membraneXY[eid][nid] = membraneXY
            self.bendingX[eid][nid] = bendingX
            self.bendingY[eid][nid] = bendingY
            self.bendingXY[eid][nid] = bendingXY
            self.shearYZ[eid][nid] = shearYZ
            self.shearXZ[eid][nid] = shearXZ

    def add_sort1(self, nNodes, dt, data):
        [eid, parent, coord, icord, theta, forces] = data
        self._fillObject(dt, eid, parent, coord, icord, theta, forces)

    def addSort2(self, nNodes, eid, data):
        [dt, parent, coord, icord, theta, forces] = data
        self._fillObject(dt, eid, parent, coord, icord, theta, forces)

    def _fillObject(self, dt, eid, parent, coord, icord, theta, forces):
        if dt not in self.membraneX:
            self.add_new_transient(dt)
        self.parent[eid] = parent
        self.coord[eid] = coord
        self.icord[eid] = icord
        self.theta[eid] = theta
        #self.eType[eid]    = eType

        self.membraneX[dt][eid] = {}
        self.membraneY[dt][eid] = {}
        self.membraneXY[dt][eid] = {}
        self.bendingX[dt][eid] = {}
        self.bendingY[dt][eid] = {}
        self.bendingXY[dt][eid] = {}
        self.shearYZ[dt][eid] = {}
        self.shearXZ[dt][eid] = {}

        for force in forces:
            [nid, membraneX, membraneY, membraneXY, bendingX,
                bendingY, bendingXY, shearYZ, shearXZ] = force
            self.membraneX[dt][eid][nid] = membraneX
            self.membraneY[dt][eid][nid] = membraneY
            self.membraneXY[dt][eid][nid] = membraneXY
            self.bendingX[dt][eid][nid] = bendingX
            self.bendingY[dt][eid][nid] = bendingY
            self.bendingXY[dt][eid][nid] = bendingXY
            self.shearYZ[dt][eid][nid] = shearYZ
            self.shearXZ[dt][eid][nid] = shearXZ

    def __repr__(self):
        return str(self.membraneX)
