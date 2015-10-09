from __future__ import (nested_scopes, generators, division, absolute_import,
                        print_function, unicode_literals)
from six import iteritems
from pyNastran.op2.resultObjects.op2_Objects import ScalarObject
from pyNastran.f06.f06_formatting import writeFloats13E, writeImagFloats13E, get_key0, write_float_12E
from numpy import zeros

class ComplexRodForce(ScalarObject):  # 1-ROD, 3-TUBE, 10-CONROD
    def __init__(self, data_code, is_sort1, isubcase, dt):
        ScalarObject.__init__(self, data_code, isubcase)
        self.axialForce = {}
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
        self.axialForce[dt] = {}
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
        msg.append('  axialForce, torque\n')
        return msg

    def add(self, dt, data):
        [eid, axialForce, torque] = data
        self.axialForce[eid] = axialForce
        self.torque[eid] = torque

    def add_sort1(self, dt, data):
        [eid, axialForce, torque] = data
        if dt not in self.axialForce:
            self.add_new_transient(dt)
        self.axialForce[dt][eid] = axialForce
        self.torque[dt][eid] = torque

    def add_sort2(self, eid, data):
        [dt, axialForce, torque] = data
        if dt not in self.axialForce:
            self.add_new_transient(dt)
        self.axialForce[dt][eid] = axialForce
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
        self._fillObject(dt, eid, nid, sd, bm1, bm2, ts1, ts2, af, ttrq, wtrq)

    def addNewElementSort2(self, eid, data):
        [dt, nid, sd, bm1, bm2, ts1, ts2, af, ttrq, wtrq] = data
        self._fillNewObject(
            dt, eid, nid, sd, bm1, bm2, ts1, ts2, af, ttrq, wtrq)

    def add_sort2(self, eid, data):
        [dt, nid, sd, bm1, bm2, ts1, ts2, af, ttrq, wtrq] = data
        self._fillObject(dt, eid, nid, sd, bm1, bm2, ts1, ts2, af, ttrq, wtrq)

    def _fillObject(self, dt, eid, nid, sd, bm1, bm2, ts1, ts2, af, ttrq, wtrq):
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
        self._fillObject(dt, eid, f41, f21, f12, f32, f23, f43, f34, f14,
                         kf1, s12, kf2, s23, kf3, s34, kf4, s41)

    def add_sort2(self, eid, data):
        [dt, f41, f21, f12, f32, f23, f43, f34, f14,
         kf1, s12, kf2, s23, kf3, s34, kf4, s41] = data
        self._fillObject(dt, eid, f41, f21, f12, f32, f23, f43, f34, f14,
                         kf1, s12, kf2, s23, kf3, s34, kf4, s41)

    def _fillObject(self, dt, eid, f41, f21, f12, f32, f23, f43, f34, f14,
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
        #packs = []
        forces = []
        elements = []
        line = '   '
        for eid, force in sorted(self.force.items()):
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
                ([forceReal, forceImag], is_all_zeros) = writeFloats13E([force.real, force.imag])

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
        self.axialForce = {}
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
        msg.append('  axialForce, torque\n')
        return msg

    def add_new_transient(self, dt):
        self.dt = dt
        self.axialForce[dt] = {}
        self.torque[dt] = {}

    def add(self, dt, data):
        [eid, axialForce, torque] = data
        self.axialForce[eid] = axialForce
        self.torque[eid] = torque

    def add_sort1(self, dt, data):
        [eid, axialForce, torque] = data
        if dt not in self.axialForce:
            self.add_new_transient(dt)
        self.axialForce[dt][eid] = axialForce
        self.torque[dt][eid] = torque

    def add_sort2(self, eid, data):
        [dt, axialForce, torque] = data
        if dt not in self.axialForce:
            self.add_new_transient(dt)
        self.axialForce[dt][eid] = axialForce
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

                ([mxr, myr, mxyr, bmxr, bmyr, bmxyr, txr, tyr,
                  mxi, myi, mxyi, bmxi, bmyi, bmxyi, txi, tyi], is_all_zeros) = writeImagFloats13E([mx, my, mxy, bmx, bmy, bmxy, tx, ty], is_mag_phase)

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

class ComplexPlate2Force(ScalarObject):  # 64-CQUAD8, 75-CTRIA6, 82-CQUADR
    def __init__(self, data_code, is_sort1, isubcase, dt):
        ScalarObject.__init__(self, data_code, isubcase)
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
                self.add_new_element = self.addNewElementSort1
                self.add = self.add_sort1
        else:
            assert dt is not None
            self.add_new_element = self.addNewElementSort2
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

    def add_new_element(self, eid, dt, data):
        [term, nid, mx, my, mxy, bmx, bmy, bmxy, tx, ty] = data
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

    def add_sort2(self, dt, eid, data):
        [nid, mx, my, mxy, bmx, bmy, bmxy, tx, ty] = data
        if dt not in self.mx:
            self.add_new_transient(dt)
        self.mx[dt][eid].append(mx)
        self.my[dt][eid].append(my)
        self.mxy[dt][eid].append(mxy)
        self.bmx[dt][eid].append(bmx)
        self.bmy[dt][eid].append(bmy)
        self.bmxy[dt][eid].append(bmxy)
        self.tx[dt][eid].append(tx)
        self.ty[dt][eid].append(ty)

    def write_f06(self, header, page_stamp, page_num=1, f=None, is_mag_phase=False, is_sort1=True):
        if self.nonlinear_factor is not None:
            return self._write_f06_transient(header, page_stamp, page_num, f, is_mag_phase=is_mag_phase, is_sort1=is_sort1)
        f.write('%s write_f06 not implemented...\n' % self.__class__.__name__)
        #raise NotImplementedError()
        #words = ['                                             A C C E L E R A T I O N   V E C T O R\n',
        #       ' \n',
        #       '      POINT ID.   TYPE          T1             T2             T3             R1             R2             R3\n']
        #words += self.getTableMarker()
        #return self._writeF06Block(words, header, page_stamp, page_num, f)

    def _write_f06_transient(self, header, page_stamp, page_num=1, f=None, is_mag_phase=False, is_sort1=True):
        f.write('%s _write_f06_transient not implemented...\n' % self.__class__.__name__)
        return page_num
        #raise NotImplementedError()


class ComplexCBarForceArray(ScalarObject):
    def get_headers(self):
        headers = ['bendingMomentA', 'bendingMomentB', 'shear', 'axial', 'torque', ]
        return headers

    def __init__(self, data_code, is_sort1, isubcase, dt):
        #ForceObject.__init__(self, data_code, isubcase)
        ScalarObject.__init__(self, data_code, isubcase)

        #self.eType = {}
        self.result_flag = 0
        #self.code = [self.format_code, self.sort_code, self.s_code]

        #self.ntimes = 0  # or frequency/mode
        #self.ntotal = 0
        self.itime = 0
        self.nelements = 0  # result specific
        self.element_type = 'CBAR'
        #self.cid = {}  # gridGauss

        if is_sort1:
            #sort1
            pass
        else:
            raise NotImplementedError('SORT2')

    def _reset_indices(self):
        self.itotal = 0
        self.ielement = 0

    def build(self):
        #print('ntimes=%s nelements=%s ntotal=%s subtitle=%s' % (self.ntimes, self.nelements, self.ntotal, self.subtitle))
        if self.is_built:
            return
        nnodes = 1

        #self.names = []
        #self.nelements //= nnodes
        self.nelements /= self.ntimes
        #self.ntotal //= self.ntimes
        self.itime = 0
        self.ielement = 0
        self.itotal = 0
        self.is_built = True
        #print('ntotal=%s ntimes=%s nelements=%s' % (self.ntotal, self.ntimes, self.nelements))

        #print("ntimes=%s nelements=%s ntotal=%s" % (self.ntimes, self.nelements, self.ntotal))
        self._times = zeros(self.ntimes, 'float32')
        self.element = zeros(self.ntotal, 'int32')

        # the number is messed up because of the offset for the element's properties

        if not self.nelements * nnodes == self.ntotal:
            msg = 'ntimes=%s nelements=%s nnodes=%s ne*nn=%s ntotal=%s' % (self.ntimes,
                                                                           self.nelements, nnodes,
                                                                           self.nelements * nnodes,
                                                                           self.ntotal)
            raise RuntimeError(msg)
        #[bm1a, bm2a, bm1b, bm2b, ts1, ts2, af, trq]
        self.data = zeros((self.ntimes, self.ntotal, 8), 'complex64')

    def add_sort1(self, dt, eid, bm1a, bm2a, bm1b, bm2b, ts1, ts2, af, trq):
        #[bendingMomentA, bendingMomentB, shear, axial, torque]
        self._times[self.itime] = dt
        self.data[self.itime, self.itotal, :] = [bm1a, bm2a, bm1b, bm2b, ts1, ts2, af, trq]
        self.element[self.itotal] = eid
        self.itotal += 1

    def get_stats(self):
        if not self.is_built:
            return [
                '<%s>\n' % self.__class__.__name__,
                '  ntimes: %i\n' % self.ntimes,
                '  ntotal: %i\n' % self.ntotal,
            ]

        nelements = self.nelements
        ntimes = self.ntimes
        #ntotal = self.ntotal
        msg = []

        if self.nonlinear_factor is not None:  # transient
            msg.append('  type=%s ntimes=%i nelements=%i\n'
                       % (self.__class__.__name__, ntimes, nelements))
        else:
            msg.append('  type=%s nelements=%i\n' % (self.__class__.__name__, nelements))
        msg.append('  eType, cid\n')
        msg.append('  data: [ntimes, nelements, 8] where 8=[%s]\n' % str(', '.join(self.get_headers())))
        msg.append('  data.shape = %s\n' % str(self.data.shape).replace('L', ''))
        msg.append('  is_sort1=%s is_sort2=%s\n' % (self.is_sort1(), self.is_sort2()))
        msg.append('  CBAR\n  ')
        msg += self.get_data_code()
        return msg

    def write_f06(self, header, page_stamp, page_num=1, f=None, is_mag_phase=False, is_sort1=True):
        #msg_temp, nnodes = get_f06_header(self, is_mag_phase, is_sort1)


        #is_sort1 = False
        if is_mag_phase:
            mag_phase = '                                                          (MAGNITUDE/PHASE)\n \n'
        else:
            mag_phase = '                                                          (REAL/IMAGINARY)\n \n'


        name = self.data_code['name']
        if name == 'freq':
            name = 'FREQUENCY'
        #else: # mode
            #raise RuntimeError(name)

        if is_sort1:
            line1 = '0    ELEMENT         BEND-MOMENT-END-A            BEND-MOMENT-END-B                  SHEAR\n'
            line2 = '       ID.         PLANE 1       PLANE 2        PLANE 1       PLANE 2        PLANE 1       PLANE 2        FORCE          TORQUE\n'
        else:
            line1 = '                    BEND-MOMENT-END-A            BEND-MOMENT-END-B                  SHEAR\n'
            line2 = '   %16s       PLANE 1       PLANE 2        PLANE 1       PLANE 2        PLANE 1       PLANE 2        FORCE          TORQUE\n' % name

        # force
        msg_temp = header + [
            '                             C O M P L E X   F O R C E S   I N   B A R   E L E M E N T S   ( C B A R )\n',
            mag_phase,
            ' ',
            line1,
            line2,
        ]
        if self.is_sort1():
            assert self.is_sort1() == True, str(self)
            if is_sort1:
                page_num = self._write_sort1_as_sort1(f, page_num, page_stamp, header, msg_temp, is_mag_phase)
            else:
                self._write_sort1_as_sort2(f, page_num, page_stamp, header, msg_temp, is_mag_phase)
        else:
            assert self.is_sort1() == True, str(self)
        return page_num - 1

    def _write_sort1_as_sort1(self, f, page_num, page_stamp, header, msg_temp, is_mag_phase):
        eids = self.element
        times = self._times
        ntimes = self.data.shape[0]
        for itime in range(ntimes):
            dt = self._times[itime]
            dt_line = ' %14s = %12.5E\n' % (self.data_code['name'], dt)
            header[1] = dt_line
            msg = header + msg_temp
            f.write(''.join(msg))

            # TODO: can I get this without a reshape?
            #bm1a, bm2a, bm1b, bm2b, ts1, ts2, af, trq
            assert self.is_sort1() == True, str(self)
            bm1a = self.data[itime, :, 0]
            bm2a = self.data[itime, :, 1]
            bm1b = self.data[itime, :, 2]
            bm2b = self.data[itime, :, 3]
            ts1 = self.data[itime, :, 4]
            ts2 = self.data[itime, :, 5]
            af = self.data[itime, :, 6]
            trq = self.data[itime, :, 7]

            for eid, bm1ai, bm2ai, bm1bi, bm2bi, ts1i, ts2i, afi, trqi in zip(eids, bm1a, bm2a, bm1b, bm2b, ts1, ts2, af, trq):
                vals = (bm1ai, bm2ai, bm1bi, bm2bi, ts1i, ts2i, afi, trqi)
                (vals2, is_all_zeros) = writeImagFloats13E(vals, is_mag_phase)
                (bm1air, bm2air, bm1bir, bm2bir, ts1ir, ts2ir, afir, trqir,
                 bm1aii, bm2aii, bm1bii, bm2bii, ts1ii, ts2ii, afii, trqii) = vals2

                f.write('0%16i   %-13s  %-13s  %-13s  %-13s  %-13s  %-13s  %-13s  %s\n'
                        ' %14s   %-13s  %-13s  %-13s  %-13s  %-13s  %-13s  %-13s  %s\n' % (
                            eid, bm1air, bm2air, bm1bir, bm2bir, ts1ir, ts2ir, afir, trqir,
                            '', bm1aii, bm2aii, bm1bii, bm2bii, ts1ii, ts2ii, afii, trqii))
            f.write(page_stamp % page_num)
            page_num += 1
        return page_num

    def _write_sort1_as_sort2(self, f, page_num, page_stamp, header, msg_temp, is_mag_phase):
        eids = self.element
        times = self._times
        ntimes = self.data.shape[0]
        for ieid, eid in enumerate(eids):
            eid_line = ' ELEMENT-ID = %s' % (eid)
            header[1] = eid_line
            msg = header + msg_temp
            f.write(''.join(msg))

            # TODO: can I get this without a reshape?
            #bm1a, bm2a, bm1b, bm2b, ts1, ts2, af, trq
            bm1a = self.data[:, ieid, 0]
            bm2a = self.data[:, ieid, 1]
            bm1b = self.data[:, ieid, 2]
            bm2b = self.data[:, ieid, 3]
            ts1 = self.data[:, ieid, 4]
            ts2 = self.data[:, ieid, 5]
            af = self.data[:, ieid, 6]
            trq = self.data[:, ieid, 7]

            for dt, bm1ai, bm2ai, bm1bi, bm2bi, ts1i, ts2i, afi, trqi in zip(times, bm1a, bm2a, bm1b, bm2b, ts1, ts2, af, trq):
                vals = (bm1ai, bm2ai, bm1bi, bm2bi, ts1i, ts2i, afi, trqi)
                (vals2, is_all_zeros) = writeImagFloats13E(vals, is_mag_phase)
                (bm1air, bm2air, bm1bir, bm2bir, ts1ir, ts2ir, afir, trqir,
                 bm1aii, bm2aii, bm1bii, bm2bii, ts1ii, ts2ii, afii, trqii) = vals2

                f.write('0%16s  %-13s  %-13s  %-13s  %-13s  %-13s  %-13s  %-13s  %s\n'
                        ' %15s  %-13s  %-13s  %-13s  %-13s  %-13s  %-13s  %-13s  %s\n' % (
                            write_float_12E(dt),
                            bm1air, bm2air, bm1bir, bm2bir, ts1ir, ts2ir, afir, trqir,
                            '', bm1aii, bm2aii, bm1bii, bm2bii, ts1ii, ts2ii, afii, trqii))
            f.write(page_stamp % page_num)
            page_num += 1
        return page_num

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


class ComplexBendForce(ScalarObject):  # 69-CBEND
    def __init__(self, data_code, is_sort1, isubcase, dt):
        ScalarObject.__init__(self, data_code, isubcase)
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

    def add_sort2(self, eid, data):
        [dt, nidA, bm1A, bm2A, sp1A, sp2A, axialA, torqueA,
            nidB, bm1B, bm2B, sp1B, sp2B, axialB, torqueB] = data
        self._fillObject(
            dt, eid, nidA, bm1A, bm2A, sp1A, sp2A, axialA, torqueA,
            nidB, bm1B, bm2B, sp1B, sp2B, axialB, torqueB)

    def _fillObject(
        self, dt, eid, nidA, bm1A, bm2A, sp1A, sp2A, axialA, torqueA,
        nidB, bm1B, bm2B, sp1B, sp2B, axialB, torqueB):
        if dt not in self.axial:
            self.add_new_transient(dt)
        self.nodeIDs[eid] = [nidA, nidB]
        self.bendingMoment1[dt][eid] = [bm1A, bm1B]
        self.bendingMoment2[dt][eid] = [bm2A, bm2B]
        self.shearPlane1[dt][eid] = [sp1A, sp1B]
        self.shearPlane2[dt][eid] = [sp2A, sp2B]
        self.axial[dt][eid] = [axialA, axialB]
        self.torque[dt][eid] = [torqueA, torqueB]


class ComplexPentaPressureForce(ScalarObject):  # 76-CHEXA_PR,77-PENTA_PR,78-TETRA_PR
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
        [eid, eName, ax, ay, az, vx, vy, vz, pressure] = data
        self.acceleration[eid] = [ax, ay, az]
        self.velocity[eid] = [vx, vy, vz]
        self.pressure[eid] = pressure

    def add_sort1(self, dt, data):
        [eid, eName, ax, ay, az, vx, vy, vz, pressure] = data
        if dt not in self.acceleration:
            self.add_new_transient(dt)
        self.acceleration[dt][eid] = [ax, ay, az]
        self.velocity[dt][eid] = [vx, vy, vz]
        self.pressure[dt][eid] = pressure

    def add_sort2(self, eid, data):
        [dt, eName, ax, ay, az, vx, vy, vz, pressure] = data
        if dt not in self.acceleration:
            self.add_new_transient(dt)
        self.acceleration[dt][eid] = [ax, ay, az]
        self.velocity[dt][eid] = [vx, vy, vz]
        self.pressure[dt][eid] = pressure


class ComplexCBushForceArray(ScalarObject):
    def get_headers(self):
        headers = ['fx', 'fy', 'fz', 'mx', 'my', 'mz']
        return headers

    def __init__(self, data_code, is_sort1, isubcase, dt):
        ScalarObject.__init__(self, data_code, isubcase)

        self.result_flag = 0
        #self.code = [self.format_code, self.sort_code, self.s_code]

        #self.ntimes = 0  # or frequency/mode
        #self.ntotal = 0
        self.itime = 0
        self.nelements = 0  # result specific
        self.element_type = 'CBUSH'

    def _reset_indices(self):
        self.itotal = 0
        self.ielement = 0

    def build(self):
        #print('ntimes=%s nelements=%s ntotal=%s subtitle=%s' % (self.ntimes, self.nelements, self.ntotal, self.subtitle))
        if self.is_built:
            return
        nnodes = 1

        #self.names = []
        #self.nelements //= nnodes
        self.nelements /= self.ntimes
        #self.ntotal //= self.ntimes
        self.itime = 0
        self.ielement = 0
        self.itotal = 0
        self.is_built = True
        #print('ntotal=%s ntimes=%s nelements=%s' % (self.ntotal, self.ntimes, self.nelements))

        #print("ntimes=%s nelements=%s ntotal=%s" % (self.ntimes, self.nelements, self.ntotal))
        self._times = zeros(self.ntimes, 'float32')
        self.element = zeros(self.ntotal, 'int32')

        # the number is messed up because of the offset for the element's properties

        if not self.nelements * nnodes == self.ntotal:
            msg = 'ntimes=%s nelements=%s nnodes=%s ne*nn=%s ntotal=%s' % (self.ntimes,
                                                                           self.nelements, nnodes,
                                                                           self.nelements * nnodes,
                                                                           self.ntotal)
            raise RuntimeError(msg)
        #[fx, fy, fz, mx, my, mz]
        self.data = zeros((self.ntimes, self.ntotal, 6), 'complex64')

    def add_sort1(self, dt, eid, fx, fy, fz, mx, my, mz):
        #[fx, fy, fz, mx, my, mz]
        self._times[self.itime] = dt
        self.data[self.itime, self.itotal, :] = [fx, fy, fz, mx, my, mz]
        self.element[self.itotal] = eid
        self.itotal += 1

    def get_stats(self):
        if not self.is_built:
            return [
                '<%s>\n' % self.__class__.__name__,
                '  ntimes: %i\n' % self.ntimes,
                '  ntotal: %i\n' % self.ntotal,
            ]

        nelements = self.nelements
        ntimes = self.ntimes
        #ntotal = self.ntotal
        msg = []

        if self.nonlinear_factor is not None:  # transient
            msg.append('  type=%s ntimes=%i nelements=%i\n'
                       % (self.__class__.__name__, ntimes, nelements))
        else:
            msg.append('  type=%s nelements=%i\n' % (self.__class__.__name__, nelements))
        msg.append('  eType, cid\n')
        msg.append('  data: [ntimes, nelements, 6] where 6=[%s]\n' % str(', '.join(self.get_headers())))
        msg.append('  data.shape = %s\n' % str(self.data.shape).replace('L', ''))
        # msg.append('  is_sort1=%s is_sort2=%s\n' % (self.is_sort1(), self.is_sort2()))
        msg.append('  CBUSH\n  ')
        msg += self.get_data_code()
        return msg

    def write_f06(self, header, page_stamp, page_num=1, f=None, is_mag_phase=False, is_sort1=True):
        #msg_temp, nnodes = get_f06_header(self, is_mag_phase, is_sort1)

        # write the f06

        #is_sort1 = False
        if is_mag_phase:
            mag_phase = '                                                          (MAGNITUDE/PHASE)\n\n'
        else:
            mag_phase = '                                                          (REAL/IMAGINARY)\n\n'


        name = self.data_code['name']
        if name == 'freq':
            name = 'FREQUENCY'
        else:
            raise RuntimeError(name)

        # is_sort1 = True
        if is_sort1:
            line2 = '       ID.     FORCE-X       FORCE-Y       FORCE-Z      MOMENT-X      MOMENT-Y      MOMENT-Z  \n'
        else:
            line2 = '   %26s        FORCE-X       FORCE-Y       FORCE-Z      MOMENT-X      MOMENT-Y      MOMENT-Z  \n' % name

        # force
        msg_temp = header + [
            '                         C O M P L E X   F O R C E S   I N   B U S H   E L E M E N T S   ( C B U S H ) \n',
            mag_phase,
            ' ',
            # line1,
            line2,
        ]
        if self.is_sort1():
            if is_sort1:
                page_num = self._write_sort1_as_sort1(f, page_num, page_stamp, header, msg_temp, is_mag_phase)
            else:
                page_num = self._write_sort1_as_sort2(f, page_num, page_stamp, header, msg_temp, is_mag_phase)
        else:
            assert self.is_sort1() == True, str(self)
        return page_num - 1

    def _write_sort1_as_sort1(self, f, page_num, page_stamp, header, msg_temp, is_mag_phase):
        ntimes = self.data.shape[0]
        eids = self.element
        for itime in range(ntimes):
            dt = self._times[itime]
            dt_line = ' %14s = %12.5E\n' % (self.data_code['name'], dt)
            header[1] = dt_line
            msg = header + msg_temp
            f.write(''.join(msg))

            # TODO: can I get this without a reshape?
            #fx, fy, fz, mx, my, mz
            if self.is_sort1():
                fx = self.data[itime, :, 0]
                fy = self.data[itime, :, 1]
                fz = self.data[itime, :, 2]
                mx = self.data[itime, :, 3]
                my = self.data[itime, :, 4]
                mz = self.data[itime, :, 5]
            else:
                fx = self.data[:, itime, 0]
                fy = self.data[:, itime, 1]
                fz = self.data[:, itime, 2]
                mx = self.data[:, itime, 3]
                my = self.data[:, itime, 4]
                mz = self.data[:, itime, 5]
                af = self.data[:, itime, 6]


            for eid, fxi, fyi, fzi, mxi, myi, mzi in zip(eids, fx, fy, fz, mx, my, mz):
                vals = (fxi, fyi, fzi, mxi, myi, mzi)
                (vals2, is_all_zeros) = writeImagFloats13E(vals, is_mag_phase)
                (fxir, fyir, fzir, mxir, myir, mzir,
                 fxii, fyii, fzii, mxii, myii, mzii) = vals2

                f.write('0%26i   %-13s  %-13s  %-13s  %-13s  %-13s  %s\n'
                        ' %26s   %-13s  %-13s  %-13s  %-13s  %-13s  %s\n' % (
                            eid, fxir, fyir, fzir, mxir, myir, mzir,
                            '', fxii, fyii, fzii, mxii, myii, mzii))
            f.write(page_stamp % page_num)
            page_num += 1
        return page_num

    def _write_sort1_as_sort2(self, f, page_num, page_stamp, header, msg_temp, is_mag_phase):
        eids = self.element
        times = self._times
        for ieid, eid in enumerate(eids):
            eid_line = ' ELEMENT-ID = %s' % (eid)
            header[1] = eid_line
            msg = header + msg_temp
            f.write(''.join(msg))

            # TODO: can I get this without a reshape?
            if self.is_sort1():
                fx = self.data[:, ieid, 0]
                fy = self.data[:, ieid, 1]
                fz = self.data[:, ieid, 2]
                mx = self.data[:, ieid, 3]
                my = self.data[:, ieid, 4]
                mz = self.data[:, ieid, 5]
            else:
                raise RuntimeError()

            for dt, fxi, fyi, fzi, mxi, myi, mzi in zip(times, fx, fy, fz, mx, my, mz):
                vals = (fxi, fyi, fzi, mxi, myi, mzi)
                (vals2, is_all_zeros) = writeImagFloats13E(vals, is_mag_phase)
                (fxir, fyir, fzir, mxir, myir, mzir,
                 fxii, fyii, fzii, mxii, myii, mzii) = vals2

                f.write('0%26s   %-13s  %-13s  %-13s  %-13s  %-13s  %s\n'
                        ' %26s   %-13s  %-13s  %-13s  %-13s  %-13s  %s\n' % (
                            write_float_12E(dt),
                            fxir, fyir, fzir, mxir, myir, mzir,
                            '', fxii, fyii, fzii, mxii, myii, mzii))
            f.write(page_stamp % page_num)
            page_num += 1
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


class ComplexForce_VU(ScalarObject):  # 191-VUBEAM
    def __init__(self, data_code, is_sort1, isubcase, dt):
        ScalarObject.__init__(self, data_code, isubcase)
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
            self.add = self.add_sort2

    def get_stats(self):
        msg = ['  '] + self.get_data_code()
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

    def add(self, nnodes, dt, data):
        [eid, parent, coord, icord, forces] = data
        self.parent[eid] = parent
        self.coord[eid] = coord
        self.icord[eid] = icord

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

    def add_sort1(self, nnodes, dt, data):
        [eid, parent, coord, icord, forces] = data
        if dt not in self.forceX:
            self.add_new_transient(dt)
        self.parent[eid] = parent
        self.coord[eid] = coord
        self.icord[eid] = icord

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

    def add_sort2(self, nnodes, eid, data):
        [dt, parent, coord, icord, forces] = data
        if dt not in self.forceX:
            self.add_new_transient(dt)
        self.parent[eid] = parent
        self.coord[eid] = coord
        self.icord[eid] = icord

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


class ComplexForce_VU_2D(ScalarObject):  # 189-VUQUAD,190-VUTRIA
    def __init__(self, data_code, is_sort1, isubcase, dt):
        ScalarObject.__init__(self, data_code, isubcase)
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
            self.add = self.add_sort2

    def get_stats(self):
        msg = ['  '] + self.get_data_code()
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

    def add(self, nnodes, dt, data):
        [eid, parent, coord, icord, theta, forces] = data
        self.parent[eid] = parent
        self.coord[eid] = coord
        self.icord[eid] = icord
        self.theta[eid] = theta

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

    def add_sort1(self, nnodes, dt, data):
        [eid, parent, coord, icord, theta, forces] = data
        self._fillObject(dt, eid, parent, coord, icord, theta, forces)

    def add_sort2(self, nnodes, eid, data):
        [dt, parent, coord, icord, theta, forces] = data
        self._fillObject(dt, eid, parent, coord, icord, theta, forces)

    def _fillObject(self, dt, eid, parent, coord, icord, theta, forces):
        if dt not in self.membraneX:
            self.add_new_transient(dt)
        self.parent[eid] = parent
        self.coord[eid] = coord
        self.icord[eid] = icord
        self.theta[eid] = theta

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
