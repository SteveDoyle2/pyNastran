#pylint disable=C0301
from __future__ import (nested_scopes, generators, division, absolute_import,
                        print_function, unicode_literals)
from six import iteritems
from six.moves import zip, range
from numpy import array, zeros, searchsorted, array_equal
from itertools import cycle

from pyNastran.op2.resultObjects.op2_Objects import ScalarObject
from pyNastran.op2.tables.oes_stressStrain.real.oes_springs import _write_f06_springs
from pyNastran.f06.f06_formatting import writeFloats13E, writeFloats12E, _eigenvalue_header, get_key0

class RealRodForceArray(ScalarObject):
    def __init__(self, data_code, is_sort1, isubcase, dt):
        self.element_type = None
        self.element_name = None
        ScalarObject.__init__(self, data_code, isubcase)
        #self.code = [self.format_code, self.sort_code, self.s_code]

        #self.ntimes = 0  # or frequency/mode
        #self.ntotal = 0
        self.nelements = 0  # result specific

        if is_sort1:
            self.add = self.add_sort1
        else:
            raise NotImplementedError('SORT2')

    def _reset_indices(self):
        self.itotal = 0
        self.ielement = 0

    def get_headers(self):
        headers = ['axial', 'torsion']
        return headers

    #def get_headers(self):
        #headers = ['axial', 'torque']
        #return headers

    def _get_msgs(self):
        base_msg = ['       ELEMENT       AXIAL       TORSIONAL     ELEMENT       AXIAL       TORSIONAL\n',
                    '         ID.         FORCE        MOMENT        ID.          FORCE        MOMENT\n']
        crod_msg = ['                                     F O R C E S   I N   R O D   E L E M E N T S      ( C R O D )\n', ]
        conrod_msg = ['                                     F O R C E S   I N   R O D   E L E M E N T S      ( C O N R O D )\n', ]
        ctube_msg = ['                                     F O R C E S   I N   R O D   E L E M E N T S      ( C T U B E )\n', ]
        crod_msg += base_msg
        conrod_msg += base_msg
        ctube_msg += base_msg
        return crod_msg, conrod_msg, ctube_msg

    def build(self):
        #print('ntimes=%s nelements=%s ntotal=%s' % (self.ntimes, self.nelements, self.ntotal))
        if self.is_built:
            return

        assert self.ntimes > 0, 'ntimes=%s' % self.ntimes
        assert self.nelements > 0, 'nelements=%s' % self.nelements
        assert self.ntotal > 0, 'ntotal=%s' % self.ntotal
        #self.names = []
        self.nelements //= self.ntimes
        self.itime = 0
        self.ielement = 0
        self.itotal = 0
        #self.ntimes = 0
        #self.nelements = 0
        self.is_built = True

        #print("ntimes=%s nelements=%s ntotal=%s" % (self.ntimes, self.nelements, self.ntotal))
        dtype = 'float32'
        if isinstance(self.nonlinear_factor, int):
            dtype = 'int32'
        self._times = zeros(self.ntimes, dtype=dtype)
        self.element = zeros(self.nelements, dtype='int32')

        #[axial_force, torque]
        self.data = zeros((self.ntimes, self.ntotal, 2), dtype='float32')

    def add(self, dt, eid, axial, torque):
        self.add_sort1(dt, eid, axial, torque)

    def add_sort1(self, dt, eid, axial, torque):
        self._times[self.itime] = dt
        self.element[self.ielement] = eid
        self.data[self.itime, self.ielement, :] = [axial, torque]
        self.ielement += 1

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
            ntimes_word = 'ntimes'
        else:
            msg.append('  type=%s nelements=%i\n'
                       % (self.__class__.__name__, nelements))
            ntimes_word = 1
        msg.append('  eType\n')
        headers = self.get_headers()
        n = len(headers)
        msg.append('  data: [%s, nnodes, %i] where %i=[%s]\n' % (ntimes_word, n, n, str(', '.join(headers))))
        msg.append('  data.shape = %s\n' % str(self.data.shape).replace('L', ''))
        #msg.append('  element type: %s\n' % self.element_type)
        msg.append('  element name: %s\n  ' % self.element_name)
        msg += self.get_data_code()
        return msg

    def get_f06_header(self, is_mag_phase=True):
        crod_msg, conrod_msg, ctube_msg = self._get_msgs()
        if 'CROD' in self.element_name:
            msg = crod_msg
        elif 'CONROD' in self.element_name:
            msg = conrod_msg
        elif 'CTUBE' in self.element_name:
            msg = ctube_msg
        return self.element_name, msg

    def get_element_index(self, eids):
        # elements are always sorted; nodes are not
        itot = searchsorted(eids, self.element)  #[0]
        return itot

    def eid_to_element_node_index(self, eids):
        #ind = ravel([searchsorted(self.element == eid) for eid in eids])
        ind = searchsorted(eids, self.element)
        #ind = ind.reshape(ind.size)
        #ind.sort()
        return ind

    def write_f06(self, header, page_stamp, page_num=1, f=None, is_mag_phase=False, is_sort1=True):
        (elem_name, msg_temp) = self.get_f06_header(is_mag_phase)

        # write the f06
        #(ntimes, ntotal, two) = self.data.shape
        ntimes = self.data.shape[0]

        eids = self.element
        is_odd = False
        nwrite = len(eids)
        if len(eids) % 2 == 1:
            nwrite -= 1
            is_odd = True

        #print('len(eids)=%s nwrite=%s is_odd=%s' % (len(eids), nwrite, is_odd))
        for itime in range(ntimes):
            dt = self._times[itime]  # TODO: rename this...
            header = _eigenvalue_header(self, header, itime, ntimes, dt)
            f.write(''.join(header + msg_temp))

            # TODO: can I get this without a reshape?
            #print("self.data.shape=%s itime=%s ieids=%s" % (str(self.data.shape), itime, str(ieids)))
            axial = self.data[itime, :, 0]
            torsion = self.data[itime, :, 1]

            # loop over all the elements
            out = []
            for eid, axiali, torsioni in zip(eids, axial, torsion):
                ([axiali, torsioni], is_all_zeros) = writeFloats13E([axiali, torsioni])
                out.append([eid, axiali, torsioni])

            for i in range(0, nwrite, 2):
                out_line = '      %8i   %-13s  %-13s  %8i   %-13s  %s\n' % tuple(out[i] + out[i + 1])
                f.write(out_line)
            if is_odd:
                out_line = '      %8i   %-13s  %s\n' % tuple(out[-1])
                f.write(out_line)
            f.write(page_stamp % page_num)
            page_num += 1
        return page_num - 1


class RealCBeamForce(ScalarObject):  # 2-CBEAM
    def __init__(self, data_code, is_sort1, isubcase, dt):
        self.element_type = None
        self.element_name = None
        ScalarObject.__init__(self, data_code, isubcase)
        self.nodes = {}
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
        msg = ['  '] + self.get_data_code()
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

    def add_f06_data(self, data, dt=None):
        if dt:
            raise NotImplementedError(dt)

        for d in data:
            (eid, nid, sd, bm1, bm2, ts1, ts2, af, ttrq, wtrq) = d
            #print('eid, nid, sd', eid, nid, sd)
            if eid in self.nodes:
                #if sd in self.nodes[eid]:
                self.nodes[eid][sd] = nid
                self.bendingMoment[eid][sd] = [bm1, bm2]
                self.shear[eid][sd] = [ts1, ts2]
                self.axial[eid][sd] = af
                self.totalTorque[eid][sd] = ttrq
                self.warpingTorque[eid][sd] = wtrq
            else:
                self.nodes[eid] = {sd: [nid]}
                self.bendingMoment[eid] = {sd: [bm1, bm2]}
                self.shear[eid] = {sd: [ts1, ts2]}
                self.axial[eid] = {sd: af}
                self.totalTorque[eid] = {sd: ttrq}
                self.warpingTorque[eid] = {sd: wtrq}
            #print('nodes', self.nodes)

    def add_new_element(self, dt, data):
        [eid, nid, sd, bm1, bm2, ts1, ts2, af, ttrq, wtrq] = data
        self.nodes[eid] = {sd: nid}
        self.bendingMoment[eid] = {sd: [bm1, bm2]}
        self.shear[eid] = {sd: [ts1, ts2]}
        self.axial[eid] = {sd: af}
        self.totalTorque[eid] = {sd: ttrq}
        self.warpingTorque[eid] = {sd: wtrq}

    def add(self, dt, data):
        [eid, nid, sd, bm1, bm2, ts1, ts2, af, ttrq, wtrq] = data
        self.nodes[eid][sd] = nid
        self.bendingMoment[eid][sd] = [bm1, bm2]
        self.shear[eid][sd] = [ts1, ts2]
        self.axial[eid][sd] = af
        self.totalTorque[eid][sd] = ttrq
        self.warpingTorque[eid][sd] = wtrq

    def addNewElementSort1(self, dt, data):
        [eid, nid, sd, bm1, bm2, ts1, ts2, af, ttrq, wtrq] = data
        self._fillObjectNew(dt, eid, nid, sd, bm1, bm2, ts1, ts2, af, ttrq, wtrq)

    def add_sort1(self, dt, data):
        [eid, nid, sd, bm1, bm2, ts1, ts2, af, ttrq, wtrq] = data
        self._fillObject(dt, eid, nid, sd, bm1, bm2, ts1, ts2, af, ttrq, wtrq)

    def addNewElementSort2(self, eid, data):
        [dt, nid, sd, bm1, bm2, ts1, ts2, af, ttrq, wtrq] = data
        self._fillObjectNew(
            dt, eid, nid, sd, bm1, bm2, ts1, ts2, af, ttrq, wtrq)

    def add_sort2(self, eid, data):
        [dt, nid, sd, bm1, bm2, ts1, ts2, af, ttrq, wtrq] = data
        self._fillObject(dt, eid, nid, sd, bm1, bm2, ts1, ts2, af, ttrq, wtrq)

    def _fillObject(self, dt, eid, nid, sd, bm1, bm2, ts1, ts2, af, ttrq, wtrq):
        #if dt not in self.axial:
            #self.add_new_transient(dt)
        self.nodes[eid][sd] = nid
        self.bendingMoment[dt][eid][sd] = [bm1, bm2]
        self.shear[dt][eid][sd] = [ts1, ts2]
        self.axial[dt][eid][sd] = af
        self.totalTorque[dt][eid][sd] = ttrq
        self.warpingTorque[dt][eid][sd] = wtrq

    def _fillObjectNew(self, dt, eid, nid, sd, bm1, bm2, ts1, ts2, af, ttrq, wtrq):
        if dt not in self.axial:
            self.add_new_transient(dt)
        self.nodes[eid] = {sd: nid}
        self.bendingMoment[dt][eid] = {sd: [bm1, bm2]}
        self.shear[dt][eid] = {sd: [ts1, ts2]}
        self.axial[dt][eid] = {sd: af}
        self.totalTorque[dt][eid] = {sd: ttrq}
        self.warpingTorque[dt][eid] = {sd: wtrq}

    def _write_f06_transient(self, header, page_stamp, page_num=1, f=None, is_mag_phase=False, is_sort1=True):
        words = ['                                 F O R C E S   I N   B E A M   E L E M E N T S        ( C B E A M )\n',
                 '                    STAT DIST/   - BENDING MOMENTS -            - WEB  SHEARS -           AXIAL          TOTAL          WARPING\n',
                 '   ELEMENT-ID  GRID   LENGTH    PLANE 1       PLANE 2        PLANE 1       PLANE 2        FORCE          TORQUE         TORQUE\n']

        msg = []
        for dt, bms in sorted(iteritems(self.bendingMoment)):
            header[1] = ' %s = %10.4E\n' % (self.data_code['name'], dt)
            msg += header + words
            for eid, bm in sorted(iteritems(bms)):
                for sd in sorted(bm):
                    nid = self.nodes[eid][sd]
                    bm1, bm2 = self.bendingMoment[dt][eid][sd]
                    ts1, ts2 = self.shear[dt][eid][sd]
                    af = self.axial[dt][eid][sd]
                    ttrq = self.totalTorque[dt][eid][sd]
                    wtrq = self.warpingTorque[dt][eid][sd]
                    (vals2, is_all_zeros) = writeFloats13E([bm1, bm2, ts1, ts2, af, ttrq, wtrq])
                    [bm1, bm2, ts1, ts2, af, ttrq, wtrq] = vals2

                    if sd == 0.:
                        msg.append('0  %8i\n' % (eid))

                    # TODO store grid ID
                    msg.append('           %8i   %.3f   %-13s %-13s  %-13s %-13s  %-13s  %-13s  %s\n' % (nid, sd, bm1, bm2, ts1, ts2, af, ttrq, wtrq))

            msg.append(page_stamp % page_num)
            f.write(''.join(msg))
            msg = ['']
            page_num += 1
        return page_num - 1

    def write_f06(self, header, page_stamp, page_num=1, f=None, is_mag_phase=False, is_sort1=True):
        if self.nonlinear_factor is not None:
            return self._write_f06_transient(header, page_stamp, page_num, f, is_mag_phase=is_mag_phase, is_sort1=is_sort1)
        msg = header + ['                                 F O R C E S   I N   B E A M   E L E M E N T S        ( C B E A M )\n',
                        '                    STAT DIST/   - BENDING MOMENTS -            - WEB  SHEARS -           AXIAL          TOTAL          WARPING\n',
                        '   ELEMENT-ID  GRID   LENGTH    PLANE 1       PLANE 2        PLANE 1       PLANE 2        FORCE          TORQUE         TORQUE\n']
        for eid, bm in sorted(iteritems(self.bendingMoment)):
            msg.append('0  %8i\n' % eid)
            for sd in sorted(bm):
                nid = self.nodes[eid][sd]
                bm1, bm2 = self.bendingMoment[eid][sd]
                ts1, ts2 = self.shear[eid][sd]
                af = self.axial[eid][sd]
                ttrq = self.totalTorque[eid][sd]
                wtrq = self.warpingTorque[eid][sd]
                (vals2, is_all_zeros) = writeFloats13E([bm1, bm2, ts1, ts2, af, ttrq, wtrq])
                [bm1, bm2, ts1, ts2, af, ttrq, wtrq] = vals2
                msg.append('           %8i   %.3f   %-13s %-13s  %-13s %-13s  %-13s  %-13s  %s\n' % (nid, sd, bm1, bm2, ts1, ts2, af, ttrq, wtrq))

        msg.append(page_stamp % page_num)
        f.write(''.join(msg))
        return page_num


class RealCShearForce(ScalarObject):  # 4-CSHEAR
    def __init__(self, data_code, is_sort1, isubcase, dt):
        self.element_type = None
        self.element_name = None
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

    def add_f06_data(self, data, dt=None):
        if dt:
            raise NotImplementedError(dt)

        for d in data:
            [
                eid,
                f41, f21, tau12, kick1,
                f12, f32, tau23, kick2,
                f23, f43, tau34, kick3,
                f34, f14, tau41, kick4,
            ] = d
            #print('eid, nid, sd', eid, nid, sd)
            self.force41[eid] = f41
            self.force14[eid] = f14

            self.force21[eid] = f21
            self.force12[eid] = f12

            self.force32[eid] = f32
            self.force23[eid] = f23

            self.force43[eid] = f43
            self.force34[eid] = f34

            self.kickForce1[eid] = kick1
            self.kickForce2[eid] = kick2
            self.kickForce3[eid] = kick3
            self.kickForce4[eid] = kick4

            self.shear12[eid] = tau12
            self.shear23[eid] = tau23
            self.shear34[eid] = tau34
            self.shear41[eid] = tau41

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

    def write_f06(self, header, page_stamp, page_num=1, f=None, is_mag_phase=False, is_sort1=True):
        if self.nonlinear_factor is not None:
            return self._write_f06_transient(header, page_stamp, page_num, f, is_mag_phase=is_mag_phase, is_sort1=is_sort1)

        words = [
            '                           F O R C E S   A C T I N G   O N   S H E A R   P A N E L   E L E M E N T S   (CSHEAR)\n'
            ' \n'
            '                  ====== POINT  1 ======      ====== POINT  2 ======      ====== POINT  3 ======      ====== POINT  4 ======\n'
            '   ELEMENT        F-FROM-4      F-FROM-2      F-FROM-1      F-FROM-3      F-FROM-2      F-FROM-4      F-FROM-3      F-FROM-1\n'
            '         ID               KICK-1       SHEAR-12       KICK-2       SHEAR-23       KICK-3       SHEAR-34       KICK-4       SHEAR-41\n'
        ]
        msg = header + words
        for eid, forcei in sorted(iteritems(self.force14)):
            f41 = self.force41[eid]
            f14 = self.force14[eid]
            f21 = self.force21[eid]
            f12 = self.force12[eid]
            f32 = self.force32[eid]
            f23 = self.force23[eid]
            f43 = self.force43[eid]
            f34 = self.force34[eid]
            kick1 = self.kickForce1[eid]
            kick2 = self.kickForce2[eid]
            kick3 = self.kickForce3[eid]
            kick4 = self.kickForce4[eid]
            tau12 = self.shear12[eid]
            tau23 = self.shear23[eid]
            tau34 = self.shear34[eid]
            tau41 = self.shear41[eid]

            vals = [
                f14, f12,
                f21, f23,
                f32, f34,
                f43, f41,
                kick1, tau12,
                kick2, tau23,
                kick3, tau34,
                kick4, tau41,
            ]
            (vals2, is_all_zeros) = writeFloats12E(vals)
            [f14, f12,
             f21, f23,
             f32, f34,
             f43, f41,
             kick1, tau12, kick2, tau23, kick3, tau34, kick4, tau41,
            ] = vals2
            msg.append('0%13i%-13s %-13s %-13s %-13s %-13s %-13s %-13s %s\n' % (eid, f14, f12, f21, f23, f32, f34, f43, f41))
            msg.append('                     %-13s %-13s %-13s %-13s %-13s %-13s %-13s %s\n' % (kick1, tau12, kick2, tau23, kick3, tau34, kick4, tau41))

        msg.append(page_stamp % page_num)
        f.write(''.join(msg))
        return page_num


class RealSpringForce(ScalarObject):  # 11-CELAS1,12-CELAS2,13-CELAS3, 14-CELAS4
    def __init__(self, data_code, is_sort1, isubcase, dt):
        self.element_type = None
        self.element_name = None
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

    def add_f06_data(self, data, dt):
        if dt is not None:
            for datai in data:
                (eid, forcei) = datai
                self.force[dt][eid] = forcei
            return

        for datai in data:
            (eid, forcei) = datai
            self.force[eid] = forcei

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

    def _write_f06_transient(self, header, page_stamp, page_num=1, f=None, is_mag_phase=False, is_sort1=True):
        words = ['                              F O R C E S   I N   S C A L A R   S P R I N G S        ( C E L A S 2 )\n',
                 '      ELEMENT         FORCE            ELEMENT         FORCE            ELEMENT         FORCE            ELEMENT         FORCE\n',
                 '        ID.                              ID.                              ID.                              ID.\n',
                 ]
        msg = []
        for dt, Force in sorted(iteritems(self.force)):
            header[1] = ' %s = %10.4E\n' % (self.data_code['name'], dt)
            msg += header + words

            forces = []
            #elements = []
            line = '   '
            for eid, force in sorted(iteritems(Force)):
                #elements.append(eid)
                forces.append(force)
                line += '%13s  %13s     ' % (eid, f)
                if len(forces) == 3:
                    forces = []
                    msg.append(line.rstrip() + '\n')
                    line = '   '

            if forces:
                msg.append(line.rstrip() + '\n')

            msg.append(page_stamp % page_num)
            f.write(''.join(msg))
            msg = ['']
            page_num += 1

        return page_num - 1

    def write_f06(self, header, page_stamp, page_num=1, f=None, is_mag_phase=False, is_sort1=True):
        if self.nonlinear_factor is not None:
            return self._write_f06_transient(header, page_stamp, page_num, f, is_sort1=is_sort1)
        msg = header + ['                              F O R C E S   I N   S C A L A R   S P R I N G S        ( C E L A S 2 )\n',
                        '      ELEMENT         FORCE            ELEMENT         FORCE            ELEMENT         FORCE            ELEMENT         FORCE\n',
                        '        ID.                              ID.                              ID.                              ID.\n',
                        ]
        f.write(''.join(msg))
        _write_f06_springs(f, self.force)
        f.write(page_stamp % page_num)

        return page_num


class RealDamperForce(ScalarObject):  # 20-CDAMP1,21-CDAMP2,22-CDAMP3,23-CDAMP4
    def __init__(self, data_code, is_sort1, isubcase, dt):
        self.element_type = None
        self.element_name = None
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

    def _write_f06_transient(self, header, page_stamp, page_num=1, f=None, is_mag_phase=False, is_sort1=True):
        words = ['                              F O R C E S   I N   S C A L A R   S P R I N G S        ( C E L A S 2 )\n',
                 ' \n',
                 '        TIME          FORCE              TIME          FORCE              TIME          FORCE              TIME          FORCE\n']
        msg = []
        for dt, Force in sorted(self.force.items()):
            header[1] = ' %s = %10.4E\n' % (self.data_code['name'], dt)
            msg += header + words

            #packs = []
            forces = []
            elements = []
            line = '   '
            for eid, force in sorted(Force.items()):
                elements.append(eid)
                forces.append(force)
                #pack.append(eid)
                #pack.append(f)
                line += '%13s  %13s     ' % (eid, f)
                if len(forces) == 3:
                    msg.append(line.rstrip() + '\n')

            if forces:
                msg.append(line.rstrip() + '\n')
            msg.append(page_stamp % page_num)
            f.write(''.join(msg))
            msg = ['']
            page_num += 1

        return page_num - 1

    def write_f06(self, header, page_stamp, page_num=1, f=None, is_mag_phase=False, is_sort1=True):
        if self.nonlinear_factor is not None:
            return self._write_f06_transient(header, page_stamp, page_num, f, is_sort1=is_sort1)
        msg = header + ['                              F O R C E S   I N   S C A L A R   S P R I N G S        ( C E L A S 2 )\n',
                        ' \n',
                        '        TIME          FORCE              TIME          FORCE              TIME          FORCE              TIME          FORCE\n']
        #packs = []
        forces = []
        elements = []
        line = '   '
        for eid, force in sorted(self.force.items()):
            elements.append(eid)
            forces.append(force)
            #pack.append(eid)
            #pack.append(f)
            line += '%13s  %13s     ' % (eid, force)
            if len(forces) == 3:
                msg.append(line.rstrip() + '\n')

        if forces:
            msg.append(line.rstrip() + '\n')

        msg.append(page_stamp % page_num)
        f.write(''.join(msg))
        return page_num


class RealViscForce(ScalarObject):  # 24-CVISC
    def __init__(self, data_code, is_sort1, isubcase, dt):
        self.element_type = None
        self.element_name = None
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


class RealPlateForceArray(ScalarObject):  # 33-CQUAD4, 74-CTRIA3
    def __init__(self, data_code, is_sort1, isubcase, dt):
        self.element_type = None
        self.element_name = None
        ScalarObject.__init__(self, data_code, isubcase)

        self.dt = dt
        self.nelements = 0

        if is_sort1:
            if dt is not None:
                self.add = self.add_sort1
        else:
            assert dt is not None
            self.add = self.add_sort2

    def _reset_indices(self):
        self.itotal = 0
        self.ielement = 0

    def _get_msgs(self):
        raise NotImplementedError()

    def get_headers(self):
        return ['mx', 'my', 'mxy', 'bmx', 'bmy', 'bmxy', 'tx', 'ty']

    def build(self):
        #print('ntimes=%s nelements=%s ntotal=%s' % (self.ntimes, self.nelements, self.ntotal))
        if self.is_built:
            return

        assert self.ntimes > 0, 'ntimes=%s' % self.ntimes
        assert self.nelements > 0, 'nelements=%s' % self.nelements
        assert self.ntotal > 0, 'ntotal=%s' % self.ntotal
        #self.names = []
        #self.nelements //= self.ntimes
        self.itime = 0
        self.ielement = 0
        self.itotal = 0
        #self.ntimes = 0
        #self.nelements = 0
        self.is_built = True

        #print("ntimes=%s nelements=%s ntotal=%s" % (self.ntimes, self.nelements, self.ntotal))
        dtype = 'float32'
        if isinstance(self.nonlinear_factor, int):
            dtype = 'int32'
        self._times = zeros(self.ntimes, dtype=dtype)
        self.element = zeros(self.nelements, dtype='int32')

        #[mx, my, mxy, bmx, bmy, bmxy, tx, ty]
        self.data = zeros((self.ntimes, self.ntotal, 8), dtype='float32')

    #def add_new_eid_sort1(self, dt, eid, axial, SMa, torsion, SMt):
        #self._times[self.itime] = dt
        #self.element[self.ielement] = eid
        #self.data[self.itime, self.ielement, :] = [axial, SMa, torsion, SMt]
        #self.ielement += 1

    def add(self, dt, eid, mx, my, mxy, bmx, bmy, bmxy, tx, ty):
        self.add_sort1(dt, eid, mx, my, mxy, bmx, bmy, bmxy, tx, ty)

    def add_sort1(self, dt, eid, mx, my, mxy, bmx, bmy, bmxy, tx, ty):
        self._times[self.itime] = dt
        self.element[self.itotal] = eid
        self.data[self.itime, self.itotal, :] = [mx, my, mxy, bmx, bmy, bmxy, tx, ty]
        self.itotal += 1

    def add_sort2(self, eid, mx, my, mxy, bmx, bmy, bmxy, tx, ty):
        raise NotImplementedError('SORT2')
        #if dt not in self.mx:
            #self.add_new_transient(dt)
        self.data[self.itime, self.itotal, :] = [mx, my, mxy, bmx, bmy, bmxy, tx, ty]

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
            ntimes_word = 'ntimes'
        else:
            msg.append('  type=%s nelements=%i\n'
                       % (self.__class__.__name__, nelements))
            ntimes_word = 1
        #msg.append('  eType\n')
        headers = self.get_headers()
        n = len(headers)
        msg.append('  data: [%s, nnodes, %i] where %i=[%s]\n' % (ntimes_word, n, n, str(', '.join(headers))))
        msg.append('  data.shape = %s\n' % str(self.data.shape).replace('L', ''))
        msg.append('  element type: %s\n  ' % self.element_name)
        msg += self.get_data_code()
        return msg

    def get_f06_header(self, is_mag_phase=True):
        if 'CTRIA3' in self.element_name:
            msg = [
                '                             F O R C E S   I N   T R I A N G U L A R   E L E M E N T S   ( T R I A 3 )\n'
                ' \n'
                '    ELEMENT                - MEMBRANE  FORCES -                        - BENDING  MOMENTS -              - TRANSVERSE SHEAR FORCES -\n'
                '      ID              FX            FY            FXY             MX            MY            MXY             QX            QY\n'
            ]
            nnodes = 3
        elif 'CQUAD4' in self.element_name:
            msg = [
                '                          F O R C E S   I N   Q U A D R I L A T E R A L   E L E M E N T S   ( Q U A D 4 )\n'
                ' \n'
                '    ELEMENT                    - MEMBRANE  FORCES -                      - BENDING   MOMENTS -            - TRANSVERSE SHEAR FORCES -\n'
                '      ID       GRID-ID     FX            FY            FXY           MX            MY            MXY           QX            QY\n'
            ]
            nnodes = 4
        else:
            raise NotImplementedError(self.element_name)
        return self.element_name, nnodes, msg

    def get_element_index(self, eids):
        # elements are always sorted; nodes are not
        itot = searchsorted(eids, self.element)  #[0]
        return itot

    def eid_to_element_node_index(self, eids):
        #ind = ravel([searchsorted(self.element_node[:, 0] == eid) for eid in eids])
        ind = searchsorted(eids, self.element)
        #ind = ind.reshape(ind.size)
        #ind.sort()
        return ind

    def write_f06(self, header, page_stamp, page_num=1, f=None, is_mag_phase=False, is_sort1=True):
        (elem_name, nnodes, msg_temp) = self.get_f06_header(is_mag_phase)

        # write the f06
        ntimes = self.data.shape[0]

        eids = self.element
        cen_word = 'CEN/%i' % nnodes
        for itime in range(ntimes):
            dt = self._times[itime]  # TODO: rename this...
            header = _eigenvalue_header(self, header, itime, ntimes, dt)
            f.write(''.join(header + msg_temp))

            # TODO: can I get this without a reshape?
            #print("self.data.shape=%s itime=%s ieids=%s" % (str(self.data.shape), itime, str(ieids)))
            #[mx, my, mxy, bmx, bmy, bmxy, tx, ty]
            mx = self.data[itime, :, 0]
            my = self.data[itime, :, 1]
            mxy = self.data[itime, :, 2]
            bmx = self.data[itime, :, 3]
            bmy = self.data[itime, :, 4]
            bmxy = self.data[itime, :, 5]
            tx = self.data[itime, :, 6]
            ty = self.data[itime, :, 7]

            # loop over all the elements
            if self.element_type == 74:
                for eid, mxi, myi, mxyi, bmxi, bmyi, bmxyi, txi, tyi in zip(eids, mx, my, mxy, bmx, bmy, bmxy, tx, ty):
                    ([mxi, myi, mxyi, bmxi, bmyi, bmxyi, txi, tyi],
                    is_all_zeros) = writeFloats13E([mxi, myi, mxyi, bmxi, bmyi, bmxyi, txi, tyi])
                    # ctria3
                    #          8      -7.954568E+01  2.560061E+03 -4.476376E+01    1.925648E+00  1.914048E+00  3.593237E-01    8.491534E+00  5.596094E-01  #
                    f.write('   %8i %18s %13s %13s   %13s %13s %13s   %13s %s\n' % (eid, mxi, myi, mxyi, bmxi, bmyi, bmxyi, txi, tyi))

            elif self.element_type == 33:
                for eid, mxi, myi, mxyi, bmxi, bmyi, bmxyi, txi, tyi in zip(eids, mx, my, mxy, bmx, bmy, bmxy, tx, ty):
                    ([mxi, myi, mxyi, bmxi, bmyi, bmxyi, txi, tyi],
                    is_all_zeros) = writeFloats13E([mxi, myi, mxyi, bmxi, bmyi, bmxyi, txi, tyi])
                    # cquad4
                    #0         6    CEN/4  1.072685E+01  2.504399E+03 -2.455727E+01 -5.017930E+00 -2.081427E+01 -5.902618E-01 -9.126162E+00  4.194400E+01#
                    #Fmt = '% 8i   ' + '%27.20E   ' * 8 + '\n'
                    #f.write(Fmt % (eid, mxi, myi, mxyi, bmxi, bmyi, bmxyi, txi, tyi))
                    #
                    f.write('0 %8i %8s %13s %13s %13s %13s %13s %13s %13s %s\n' % (eid, cen_word, mxi, myi, mxyi, bmxi, bmyi, bmxyi, txi, tyi))
            else:
                raise NotImplementedError(self.element_type)
            f.write(page_stamp % page_num)
            page_num += 1
        return page_num - 1


class RealPlateBilinearForceArray(ScalarObject):  # 144-CQUAD4
    def __init__(self, data_code, is_sort1, isubcase, dt):
        self.element_type = None
        self.element_name = None
        ScalarObject.__init__(self, data_code, isubcase)

        self.dt = dt
        self.nelements = 0

        if is_sort1:
            if dt is not None:
                self.add = self.add_sort1
        else:
            assert dt is not None
            self.add = self.add_sort2

    def _reset_indices(self):
        self.itotal = 0
        self.ielement = 0

    def _get_msgs(self):
        raise NotImplementedError()

    def get_headers(self):
        return ['mx', 'my', 'mxy', 'bmx', 'bmy', 'bmxy', 'tx', 'ty']

    def build(self):
        # print('ntimes=%s nelements=%s ntotal=%s' % (self.ntimes, self.nelements, self.ntotal))
        if self.is_built:
            return

        assert self.ntimes > 0, 'ntimes=%s' % self.ntimes
        assert self.nelements > 0, 'nelements=%s' % self.nelements
        assert self.ntotal > 0, 'ntotal=%s' % self.ntotal
        #self.names = []
        #self.nelements //= self.ntimes
        self.itime = 0
        self.ielement = 0
        self.itotal = 0
        #self.ntimes = 0
        #self.nelements = 0
        self.is_built = True

        #print("ntimes=%s nelements=%s ntotal=%s" % (self.ntimes, self.nelements, self.ntotal))
        dtype = 'float32'
        if isinstance(self.nonlinear_factor, int):
            dtype = 'int32'
        self._times = zeros(self.ntimes, dtype=dtype)
        self.element_node = zeros((self.ntotal, 2), dtype='int32')

        #[mx, my, mxy, bmx, bmy, bmxy, tx, ty]
        self.data = zeros((self.ntimes, self.ntotal, 8), dtype='float32')

    def __eq__(self, table):
        assert self.is_sort1() == table.is_sort1()
        assert self.nonlinear_factor == table.nonlinear_factor
        assert self.ntotal == table.ntotal
        assert self.table_name == table.table_name, 'table_name=%r table.table_name=%r' % (self.table_name, table.table_name)
        assert self.approach_code == table.approach_code
        if not array_equal(self.element_node, table.element_node):
            assert self.element_node.shape == table.element_node.shape, 'element_node shape=%s table.shape=%s' % (self.element_node.shape, table.element_node.shape)
            msg = 'table_name=%r class_name=%s\n' % (self.table_name, self.__class__.__name__)
            msg += '%s\n' % str(self.code_information())
            msg += '(Eid, Nid)\n'
            for (eid, nid), (eid2, nid2) in zip(self.element_node, table.element_node):
                msg += '(%s, %s)    (%s, %s)\n' % (eid, nid, eid2, nid2)
            print(msg)
            raise ValueError(msg)
        if not array_equal(self.data, table.data):
            msg = 'table_name=%r class_name=%s\n' % (self.table_name, self.__class__.__name__)
            msg += '%s\n' % str(self.code_information())
            i = 0
            for itime in range(self.ntimes):
                for ie, e in enumerate(self.element_node):
                    (eid, nid) = e
                    t1 = self.data[itime, ie, :]
                    t2 = table.data[itime, ie, :]
                    (mx1, my1, mxy1, bmx1, bmy1, bmxy1, tx1, ty1) = t1
                    (mx2, my2, mxy2, bmx2, bmy2, bmxy2, tx2, ty2) = t2

                    if not array_equal(t1, t2):
                        msg += '(%s, %s)    (%s, %s, %s, %s, %s, %s, %s, %s)  (%s, %s, %s, %s, %s, %s, %s, %s)\n' % (
                            eid, nid,
                            mx1, my1, mxy1, bmx1, bmy1, bmxy1, tx1, ty1,
                            mx2, my2, mxy2, bmx2, bmy2, bmxy2, tx2, ty2)
                        i += 1
                        if i > 10:
                            print(msg)
                            raise ValueError(msg)
                #print(msg)
                if i > 0:
                    raise ValueError(msg)
        return True

    def add(self, dt, eid, term, nid, mx, my, mxy, bmx, bmy, bmxy, tx, ty):
        self.add_sort1(dt, eid, term, nid, mx, my, mxy, bmx, bmy, bmxy, tx, ty)

    def add_sort1(self, dt, eid, term, nid, mx, my, mxy, bmx, bmy, bmxy, tx, ty):
        self._times[self.itime] = dt
        self.element_node[self.itotal] = [eid, nid]
        self.data[self.itime, self.itotal, :] = [mx, my, mxy, bmx, bmy, bmxy, tx, ty]
        self.itotal += 1

    def add_sort2(self, eid, mx, my, mxy, bmx, bmy, bmxy, tx, ty):
        raise NotImplementedError('SORT2')
        #if dt not in self.mx:
            #self.add_new_transient(dt)
        self.data[self.itime, self.itotal, :] = [mx, my, mxy, bmx, bmy, bmxy, tx, ty]

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
            ntimes_word = 'ntimes'
        else:
            msg.append('  type=%s nelements=%i\n'
                       % (self.__class__.__name__, nelements))
            ntimes_word = 1
        #msg.append('  eType\n')
        headers = self.get_headers()
        n = len(headers)
        msg.append('  data: [%s, nnodes, %i] where %i=[%s]\n' % (ntimes_word, n, n, str(', '.join(headers))))
        msg.append('  data.shape = %s\n' % str(self.data.shape).replace('L', ''))
        msg.append('  element type: %s\n  ' % self.element_name)
        msg += self.get_data_code()
        return msg

    def get_f06_header(self, is_mag_phase=True):
        # if 'CTRIA3' in self.element_name:
            # msg = [
                # '                             F O R C E S   I N   T R I A N G U L A R   E L E M E N T S   ( T R I A 3 )\n'
                # ' \n'
                # '    ELEMENT                - MEMBRANE  FORCES -                        - BENDING  MOMENTS -              - TRANSVERSE SHEAR FORCES -\n'
                # '      ID              FX            FY            FXY             MX            MY            MXY             QX            QY\n'
            # ]
            # nnodes = 3

        if self.element_type == 70:
            # CQUAD4
            element_name = 'CTRIAR'
            msg = [
                '                             F O R C E S   I N   T R I A N G U L A R   E L E M E N T S   ( T R I A R )\n'
                ' \n'
                '    ELEMENT                    - MEMBRANE  FORCES -                      - BENDING   MOMENTS -            - TRANSVERSE SHEAR FORCES -\n'
                '      ID       GRID-ID     FX            FY            FXY           MX            MY            MXY           QX            QY\n'
            ]
            nnodes = 6
        elif self.element_type == 75:
            # CQUAD4
            element_name = 'CTRIA6'
            msg = [
                '                             F O R C E S   I N   T R I A N G U L A R   E L E M E N T S   ( T R I A 6 )\n'
                ' \n'
                '    ELEMENT                    - MEMBRANE  FORCES -                      - BENDING   MOMENTS -            - TRANSVERSE SHEAR FORCES -\n'
                '      ID       GRID-ID     FX            FY            FXY           MX            MY            MXY           QX            QY\n'
            ]
            nnodes = 6
        elif self.element_type == 64:
            # CQUAD4
            element_name = 'CQUAD8'
            msg = [
                '                          F O R C E S   I N   Q U A D R I L A T E R A L   E L E M E N T S   ( Q U A D 8 )\n'
                ' \n'
                '    ELEMENT                    - MEMBRANE  FORCES -                      - BENDING   MOMENTS -            - TRANSVERSE SHEAR FORCES -\n'
                '      ID       GRID-ID     FX            FY            FXY           MX            MY            MXY           QX            QY\n'
            ]
            nnodes = 8
        elif self.element_type == 82:
            # CQUAD4
            element_name = 'CQUADR'
            msg = [
                '                          F O R C E S   I N   Q U A D R I L A T E R A L   E L E M E N T S   ( Q U A D R )\n'
                ' \n'
                '    ELEMENT                    - MEMBRANE  FORCES -                      - BENDING   MOMENTS -            - TRANSVERSE SHEAR FORCES -\n'
                '      ID       GRID-ID     FX            FY            FXY           MX            MY            MXY           QX            QY\n'
            ]
            nnodes = 4
        elif self.element_type == 144:
            # CQUAD4
            element_name = 'CQUAD4'
            msg = [
                '                          F O R C E S   I N   Q U A D R I L A T E R A L   E L E M E N T S   ( Q U A D 4 )        OPTION = BILIN  \n'
                ' \n'
                '    ELEMENT                    - MEMBRANE  FORCES -                      - BENDING   MOMENTS -            - TRANSVERSE SHEAR FORCES -\n'
                '      ID       GRID-ID     FX            FY            FXY           MX            MY            MXY           QX            QY\n'
            ]
            nnodes = 4
        else:
            raise NotImplementedError('element_name=%s element_type=%s' % (self.element_name, self.element_type))
        return element_name, nnodes, msg

    def get_element_index(self, eids):
        # elements are always sorted; nodes are not
        itot = searchsorted(eids, self.element)  #[0]
        return itot

    def eid_to_element_node_index(self, eids):
        #ind = ravel([searchsorted(self.element_node[:, 0] == eid) for eid in eids])
        ind = searchsorted(eids, self.element)
        #ind = ind.reshape(ind.size)
        #ind.sort()
        return ind

    def write_f06(self, header, page_stamp, page_num=1, f=None, is_mag_phase=False, is_sort1=True):
        (elem_name, nnodes, msg_temp) = self.get_f06_header(is_mag_phase)

        # write the f06
        ntimes = self.data.shape[0]

        eids = self.element_node[:, 0]
        nids = self.element_node[:, 1]
        cen_word = 'CEN/%i' % nnodes
        if self.element_type  in [64, 82, 144]: # CQUAD8, CQUADR, CQUAD4
            cyc = cycle([0, 1, 2, 3, 4])
        elif self.element_type  in [70, 75]: # CTRIAR, CTRIA6
            cyc = cycle([0, 1, 2, 3])
        else:
            raise NotImplementedError(self.element_type)

        for itime in range(ntimes):
            dt = self._times[itime]  # TODO: rename this...
            header = _eigenvalue_header(self, header, itime, ntimes, dt)
            f.write(''.join(header + msg_temp))

            # TODO: can I get this without a reshape?
            #print("self.data.shape=%s itime=%s ieids=%s" % (str(self.data.shape), itime, str(ieids)))
            #[mx, my, mxy, bmx, bmy, bmxy, tx, ty]
            mx = self.data[itime, :, 0]
            my = self.data[itime, :, 1]
            mxy = self.data[itime, :, 2]
            bmx = self.data[itime, :, 3]
            bmy = self.data[itime, :, 4]
            bmxy = self.data[itime, :, 5]
            tx = self.data[itime, :, 6]
            ty = self.data[itime, :, 7]

            # loop over all the elements
            #if self.element_type in 144:
            for i, eid, nid, mxi, myi, mxyi, bmxi, bmyi, bmxyi, txi, tyi in zip(cyc, eids, nids, mx, my, mxy, bmx, bmy, bmxy, tx, ty):
                ([mxi, myi, mxyi, bmxi, bmyi, bmxyi, txi, tyi],
                is_all_zeros) = writeFloats13E([mxi, myi, mxyi, bmxi, bmyi, bmxyi, txi, tyi])
                # ctria3
                #          8      -7.954568E+01  2.560061E+03 -4.476376E+01    1.925648E+00  1.914048E+00  3.593237E-01    8.491534E+00  5.596094E-01  #
                if i == 0:
                    f.write('0  %8i    CEN/4 %-13s %-13s %-13s %-13s %-13s %-13s %-13s %s\n' % (eid, mxi, myi, mxyi, bmxi, bmyi, bmxyi, txi, tyi))
                else:
                    f.write('            %8i %-13s %-13s %-13s %-13s %-13s %-13s %-13s %s\n' % (nid, mxi, myi, mxyi, bmxi, bmyi, bmxyi, txi, tyi))
            # else:
                # raise NotImplementedError(self.element_type)
            f.write(page_stamp % page_num)
            page_num += 1
        return page_num - 1

# class RealCBarForce(ScalarObject):  # 34-CBAR
    # def __init__(self, data_code, is_sort1, isubcase, dt):
        # self.element_type = None
        # self.element_name = None
        # ScalarObject.__init__(self, data_code, isubcase)
        # self.bendingMomentA = {}
        # self.bendingMomentB = {}
        # self.shear = {}
        # self.axial = {}
        # self.torque = {}

        # self.dt = dt
        # if is_sort1:
            # if dt is not None:
                # self.add = self.add_sort1
        # else:
            # assert dt is not None
            # self.add = self.add_sort2

    # def get_stats(self):
        # msg = ['  '] + self.get_data_code()
        # if self.dt is not None:  # transient
            # ntimes = len(self.torque)
            # time0 = get_key0(self.torque)
            # nelements = len(self.torque[time0])
            # msg.append('  type=%s ntimes=%s nelements=%s\n'
                       # % (self.__class__.__name__, ntimes, nelements))
        # else:
            # nelements = len(self.torque)
            # msg.append('  type=%s nelements=%s\n' % (self.__class__.__name__,
                                                     # nelements))
        # msg.append('  bendingMomentA, bendingMomentB, shear, axial, torque\n')
        # return msg

    # def add_new_transient(self, dt):
        # self.dt = dt
        # self.bendingMomentA[dt] = {}
        # self.bendingMomentB[dt] = {}
        # self.shear[dt] = {}
        # self.axial[dt] = {}
        # self.torque[dt] = {}

    # def add(self, dt, data):
        # [eid, bm1a, bm2a, bm1b, bm2b, ts1, ts2, af, trq] = data
        # self.bendingMomentA[eid] = [bm1a, bm2a]
        # self.bendingMomentB[eid] = [bm1b, bm2b]
        # self.shear[eid] = [ts1, ts2]
        # self.axial[eid] = af
        # self.torque[eid] = trq

    # def add_sort1(self, dt, data):
        # [eid, bm1a, bm2a, bm1b, bm2b, ts1, ts2, af, trq] = data
        # if dt not in self.axial:
            # self.add_new_transient(dt)

        # self.bendingMomentA[dt][eid] = [bm1a, bm2a]
        # self.bendingMomentB[dt][eid] = [bm1b, bm2b]
        # self.shear[dt][eid] = [ts1, ts2]
        # self.axial[dt][eid] = af
        # self.torque[dt][eid] = trq

    # def add_sort2(self, eid, data):
        # [dt, bm1a, bm2a, bm1b, bm2b, ts1, ts2, af, trq] = data
        # if dt not in self.axial:
            # self.add_new_transient(dt)

        # self.bendingMomentA[dt][eid] = [bm1a, bm2a]
        # self.bendingMomentB[dt][eid] = [bm1b, bm2b]
        # self.shear[dt][eid] = [ts1, ts2]
        # self.axial[dt][eid] = af
        # self.torque[dt][eid] = trq


    # def write_f06(self, header, page_stamp, page_num=1, f=None, is_mag_phase=False, is_sort1=True):
        # if self.nonlinear_factor is not None:
            # return self._write_f06_transient(header, page_stamp, page_num, f, is_mag_phase=False, is_sort1=is_sort1)

        # words = ['                                 F O R C E S   I N   B A R   E L E M E N T S         ( C B A R )\n',
                 # '0    ELEMENT         BEND-MOMENT END-A            BEND-MOMENT END-B                - SHEAR -               AXIAL\n',
                 # '       ID.         PLANE 1       PLANE 2        PLANE 1       PLANE 2        PLANE 1       PLANE 2         FORCE         TORQUE\n']
        # msg = []
        # #header[1] = ' %s = %10.4E\n' % (self.data_code['name'], dt)
        # f.write(''.join(words))
        # for eid in sorted(self.bendingMomentA):
            # bm1a, bm2a = self.bendingMomentA[eid]
            # bm1b, bm2b = self.bendingMomentB[eid]
            # ts1, ts2 = self.shear[eid]
            # af = self.axial[eid]
            # trq = self.torque[eid]
            # (vals2, is_all_zeros) = writeFloats13E([bm1a, bm2a, bm1b, bm2b, ts1, ts2, af, trq])
            # [bm1a, bm2a, bm1b, bm2b, ts1, ts2, af, trq] = vals2
            # f.write('     %8i    %-13s %-13s  %-13s %-13s  %-13s %-13s  %-13s  %s\n' % (eid, bm1a, bm2a, bm1b, bm2b, ts1, ts2, af, trq))
# #            1     2.504029E+06  9.728743E+06   5.088001E+05  1.976808E+06   1.995229E+06  7.751935E+06  -3.684978E-07  -1.180941E-07

            # f.write(page_stamp % page_num)
        # return page_num

    # def _write_f06_transient(self, header, page_stamp, page_num, f, is_mag_phase=False, is_sort1=True):
        # assert f is not None
        # words = ['                                 F O R C E S   I N   B A R   E L E M E N T S         ( C B A R )\n',
                 # '0    ELEMENT         BEND-MOMENT END-A            BEND-MOMENT END-B                - SHEAR -               AXIAL\n',
                 # '       ID.         PLANE 1       PLANE 2        PLANE 1       PLANE 2        PLANE 1       PLANE 2         FORCE         TORQUE\n']

        # for dt, bm in sorted(iteritems(self.bendingMomentA)):
            # header[1] = ' %s = %10.4E\n' % (self.data_code['name'], dt)
            # f.write(''.join(header + words))
            # for eid in sorted(bm):
                # bm1a, bm2a = self.bendingMomentA[dt][eid]
                # bm1b, bm2b = self.bendingMomentB[dt][eid]
                # ts1, ts2 = self.shear[dt][eid]
                # af = self.axial[dt][eid]
                # trq = self.torque[dt][eid]
                # (vals2, is_all_zeros) = writeFloats13E([bm1a, bm2a, bm1b, bm2b, ts1, ts2, af, trq])
                # [bm1a, bm2a, bm1b, bm2b, ts1, ts2, af, trq] = vals2
                # f.write('     %8i    %-13s %-13s  %-13s %-13s  %-13s %-13s  %-13s  %s\n' % (eid, bm1a, bm2a, bm1b, bm2b, ts1, ts2, af, trq))
# #            1     2.504029E+06  9.728743E+06   5.088001E+05  1.976808E+06   1.995229E+06  7.751935E+06  -3.684978E-07  -1.180941E-07

            # f.write(page_stamp % page_num)
            # page_num += 1
        # return page_num - 1

class RealCBarForceArray(ScalarObject):  # 34-CBAR
    def __init__(self, data_code, is_sort1, isubcase, dt):
        self.element_type = None
        self.element_name = None
        ScalarObject.__init__(self, data_code, isubcase)
        #self.code = [self.format_code, self.sort_code, self.s_code]

        #self.ntimes = 0  # or frequency/mode
        #self.ntotal = 0
        self.nelements = 0  # result specific

        if is_sort1:
            self.add = self.add_sort1
        else:
            raise NotImplementedError('SORT2')

    def _reset_indices(self):
        self.itotal = 0
        self.ielement = 0

    def get_headers(self):
        #self.bendingMomentA = {}
        #self.bendingMomentB = {}
        #self.shear = {}
        #self.axial = {}
        #self.torque = {}
        headers = [
            'bending_moment_a1', 'bending_moment_a2',
            'bending_moment_b1', 'bending_moment_b2',
            'shear1', 'shear2',
            'axial', 'torque']
        return headers

    def build(self):
        #print('ntimes=%s nelements=%s ntotal=%s' % (self.ntimes, self.nelements, self.ntotal))
        if self.is_built:
            return

        assert self.ntimes > 0, 'ntimes=%s' % self.ntimes
        assert self.nelements > 0, 'nelements=%s' % self.nelements
        assert self.ntotal > 0, 'ntotal=%s' % self.ntotal
        #self.names = []
        self.nelements //= self.ntimes
        self.itime = 0
        self.ielement = 0
        self.itotal = 0
        #self.ntimes = 0
        #self.nelements = 0
        self.is_built = True

        #print("ntimes=%s nelements=%s ntotal=%s" % (self.ntimes, self.nelements, self.ntotal))
        dtype = 'float32'
        if isinstance(self.nonlinear_factor, int):
            dtype = 'int32'
        self._times = zeros(self.ntimes, dtype=dtype)
        self.element = zeros(self.nelements, dtype='int32')

        #[bending_moment_a1, bending_moment_a2, bending_moment_b1, bending_moment_b2, shear1, shear2, axial, torque]
        self.data = zeros((self.ntimes, self.ntotal, 8), dtype='float32')

    def add(self, dt, data):
        self.add_sort1(dt, data)

    def add_sort1(self, dt, data):
        data = [eid, bending_moment_a1, bending_moment_a2,
                bending_moment_b1, bending_moment_b2, shear1, shear2, axial, torque] = data
        self._times[self.itime] = dt
        self.element[self.ielement] = eid
        self.data[self.itime, self.ielement, :] = [
            bending_moment_a1, bending_moment_a2,
            bending_moment_b1, bending_moment_b2,
            shear1, shear2, axial, torque]
        self.ielement += 1

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
            ntimes_word = 'ntimes'
        else:
            msg.append('  type=%s nelements=%i\n'
                       % (self.__class__.__name__, nelements))
            ntimes_word = 1
        msg.append('  eType\n')
        headers = self.get_headers()
        n = len(headers)
        msg.append('  data: [%s, nnodes, %i] where %i=[%s]\n' % (ntimes_word, n, n, str(', '.join(headers))))
        msg.append('  data.shape = %s\n' % str(self.data.shape).replace('L', ''))
        #msg.append('  element type: %s\n' % self.element_type)
        msg.append('  element name: %s\n  ' % self.element_name)
        msg += self.get_data_code()
        return msg

    def eid_to_element_node_index(self, eids):
        ind = searchsorted(eids, self.element)
        return ind

    #def _get_msgs(self):
        #base_msg = ['       ELEMENT       AXIAL       TORSIONAL     ELEMENT       AXIAL       TORSIONAL\n',
                    #'         ID.         FORCE        MOMENT        ID.          FORCE        MOMENT\n']
        #crod_msg   = ['                                     F O R C E S   I N   R O D   E L E M E N T S      ( C R O D )\n', ]
        #conrod_msg = ['                                     F O R C E S   I N   R O D   E L E M E N T S      ( C O N R O D )\n', ]
        #ctube_msg  = ['                                     F O R C E S   I N   R O D   E L E M E N T S      ( C T U B E )\n', ]
        #crod_msg += base_msg
        #conrod_msg += base_msg
        #ctube_msg += base_msg
        #return crod_msg, conrod_msg, ctube_msg

    #def get_f06_header(self, is_mag_phase=True):
        #crod_msg, conrod_msg, ctube_msg = self._get_msgs()
        #if 'CROD' in self.element_name:
            #msg = crod_msg
        #elif 'CONROD' in self.element_name:
            #msg = conrod_msg
        #elif 'CTUBE' in self.element_name:
            #msg = ctube_msg
        #return self.element_name, msg

    def write_f06(self, header, page_stamp, page_num=1, f=None, is_mag_phase=False, is_sort1=True):
        words = ['                                 F O R C E S   I N   B A R   E L E M E N T S         ( C B A R )\n',
                 '0    ELEMENT         BEND-MOMENT END-A            BEND-MOMENT END-B                - SHEAR -               AXIAL\n',
                 '       ID.         PLANE 1       PLANE 2        PLANE 1       PLANE 2        PLANE 1       PLANE 2         FORCE         TORQUE\n']
        msg = []
        #header[1] = ' %s = %10.4E\n' % (self.data_code['name'], dt)
        eids = self.element
        #f.write(''.join(words))

        ntimes = self.data.shape[0]
        for itime in range(ntimes):
            dt = self._times[itime]  # TODO: rename this...
            header = _eigenvalue_header(self, header, itime, ntimes, dt)
            f.write(''.join(header + words))

            # loop over all the elements
            bm1a = self.data[itime, :, 0]
            bm2a = self.data[itime, :, 1]
            bm1b = self.data[itime, :, 2]
            bm2b = self.data[itime, :, 3]
            ts1 = self.data[itime, :, 4]
            ts2 = self.data[itime, :, 5]
            af = self.data[itime, :, 6]
            trq = self.data[itime, :, 7]
            for eid, bm1ai, bm2ai, bm1bi, bm2bi, ts1i, ts2i, afi, trqi in zip(eids, bm1a, bm2a, bm1b, bm2b, ts1, ts2, af, trq):
                ([bm1ai, bm2ai, bm1bi, bm2bi, ts1i, ts2i, afi, trqi], is_all_zeros) = writeFloats13E([
                  bm1ai, bm2ai, bm1bi, bm2bi, ts1i, ts2i, afi, trqi])
            f.write('     %8i    %-13s %-13s  %-13s %-13s  %-13s %-13s  %-13s  %s\n' % (eid, bm1a, bm2a, bm1b, bm2b, ts1, ts2, af, trq))
            f.write(page_stamp % page_num)
        return page_num

    def __eq__(self, table):
        assert self.is_sort1() == table.is_sort1()
        assert self.nonlinear_factor == table.nonlinear_factor
        assert self.ntotal == table.ntotal
        assert self.table_name == table.table_name, 'table_name=%r table.table_name=%r' % (self.table_name, table.table_name)
        assert self.approach_code == table.approach_code
        if not array_equal(self.element, table.element):
            assert self.element.shape == table.element.shape, 'element shape=%s table.shape=%s' % (self.element.shape, table.element.shape)
            msg = 'table_name=%r class_name=%s\n' % (self.table_name, self.__class__.__name__)
            msg += '%s\n' % str(self.code_information())
            msg += 'Eid\n'
            for eid, eid2 in zip(self.element, table.element):
                msg += '%s, %s\n' % (eid, eid2)
            print(msg)
            raise ValueError(msg)
        if not array_equal(self.data, table.data):
            msg = 'table_name=%r class_name=%s\n' % (self.table_name, self.__class__.__name__)
            msg += '%s\n' % str(self.code_information())
            i = 0
            for itime in range(self.ntimes):
                for ie, e in enumerate(self.element):
                    t1 = self.data[itime, ie, :]
                    t2 = table.data[itime, ie, :]
                    (bm1a1, bm2a1, bm1b1, bm2b1, ts11, ts21, af1, trq1) = t1
                    (bm1a2, bm2a2, bm1b2, bm2b2, ts12, ts22, af2, trq2) = t2

                    if not array_equal(t1, t2):
                        msg += '(%s)    (%s, %s, %s, %s, %s, %s, %s, %s)  (%s, %s, %s, %s, %s, %s, %s, %s)\n' % (
                            eid,
                            bm1a1, bm2a1, bm1b1, bm2b1, ts11, ts21, af1, trq1,
                            bm1a2, bm2a2, bm1b2, bm2b2, ts12, ts22, af2, trq2)
                        i += 1
                        if i > 10:
                            print(msg)
                            raise ValueError(msg)
                #print(msg)
                if i > 0:
                    raise ValueError(msg)
        return True


# class RealCBar100Force(ScalarObject):  # 100-CBAR
    # def __init__(self, data_code, is_sort1, isubcase, dt):
        # self.element_type = None
        # self.element_name = None
        # ScalarObject.__init__(self, data_code, isubcase)
        # self.bendingMoment = {}
        # self.shear = {}
        # self.axial = {}
        # self.torque = {}

        # self.dt = dt
        # if is_sort1:
            # if dt is not None:
                # self.add = self.add_sort1
        # else:
            # assert dt is not None
            # self.add = self.add_sort2

    # def get_stats(self):
        # msg = ['  '] + self.get_data_code()
        # if self.dt is not None:  # transient
            # ntimes = len(self.torque)
            # time0 = get_key0(self.torque)
            # nelements = len(self.torque[time0])
            # msg.append('  type=%s ntimes=%s nelements=%s\n'
                       # % (self.__class__.__name__, ntimes, nelements))
        # else:
            # nelements = len(self.torque)
            # msg.append('  type=%s nelements=%s\n' % (self.__class__.__name__,
                                                     # nelements))
        # msg.append('  bendingMoment, shear, axial, torque\n')
        # return msg

    # def add_new_transient(self, dt):
        # self.dt = dt
        # self.bendingMoment[dt] = {}
        # self.shear[dt] = {}
        # self.axial[dt] = {}
        # self.torque[dt] = {}

    # def add(self, dt, data):
        # [eid, sd, bm1, bm2, ts1, ts2, af, trq] = data
        # self.bendingMoment[eid] = [bm1, bm2]
        # self.shear[eid] = [ts1, ts2]
        # self.axial[eid] = af
        # self.torque[eid] = trq

    # def add_sort1(self, dt, data):
        # [eid, sd, bm1, bm2, ts1, ts2, af, trq] = data
        # if dt not in self.axial:
            # self.add_new_transient(dt)
        # self.bendingMoment[dt][eid] = [bm1, bm2]
        # self.shear[dt][eid] = [ts1, ts2]
        # self.axial[dt][eid] = af
        # self.torque[dt][eid] = trq

    # def add_sort2(self, eid, data):
        # [dt, sd, bm1, bm2, ts1, ts2, af, trq] = data
        # if dt not in self.axial:
            # self.add_new_transient(dt)
        # self.bendingMoment[dt][eid] = [bm1, bm2]
        # self.shear[dt][eid] = [ts1, ts2]
        # self.axial[dt][eid] = af
        # self.torque[dt][eid] = trq


class RealConeAxForce(ScalarObject):  # 35-CCONEAX
    def __init__(self, data_code, is_sort1, isubcase, dt):
        self.element_type = None
        self.element_name = None
        ScalarObject.__init__(self, data_code, isubcase)
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
            self.add = self.add_sort2

    def get_stats(self):
        msg = ['  '] + self.get_data_code()
        if self.dt is not None:  # transient
            ntimes = len(self.hopa)
            time0 = get_key0(self.hopa)
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
        self.hopa[dt][eid] = hopa
        self.bmu[dt][eid] = bmu
        self.bmv[dt][eid] = bmv
        self.tm[dt][eid] = tm
        self.su[dt][eid] = su
        self.sv[dt][eid] = sv

    def add_sort2(self, eid, data):
        [dt, hopa, bmu, bmv, tm, su, sv] = data
        if dt not in self.hopa:
            self.add_new_transient(dt)
        self.hopa[dt][eid] = hopa
        self.bmu[dt][eid] = bmu
        self.bmv[dt][eid] = bmv
        self.tm[dt][eid] = tm
        self.su[dt][eid] = su
        self.sv[dt][eid] = sv


class RealCBar100ForceArray(ScalarObject):  # 100-CBAR
    """
    CBAR-34s are converted to CBAR-100s when you have PLOAD1s
    (distributed bar loads).  The number of stations by default is 2,
    but with a CBARAO, you can change this (max of 8 points; 6 internal
    points).

    If you use a CBARO without PLOAD1s, you wil turn CBAR-34s into
    CBAR-100s as well.
    """
    def __init__(self, data_code, is_sort1, isubcase, dt):
        self.element_type = None
        self.element_name = None
        ScalarObject.__init__(self, data_code, isubcase)
        #self.code = [self.format_code, self.sort_code, self.s_code]

        #self.ntimes = 0  # or frequency/mode
        #self.ntotal = 0
        self.nelements = 0  # result specific

        if is_sort1:
            self.add = self.add_sort1
        else:
            raise NotImplementedError('SORT2')

    def _reset_indices(self):
        self.itotal = 0
        self.ielement = 0

    def get_headers(self):
        headers = [
            'station', 'bending_moment1', 'bending_moment2', 'shear1', 'shear2', 'axial', 'torque'
        ]
        return headers

    def build(self):
        #print('ntimes=%s nelements=%s ntotal=%s' % (self.ntimes, self.nelements, self.ntotal))
        if self.is_built:
            return

        assert self.ntimes > 0, 'ntimes=%s' % self.ntimes
        assert self.nelements > 0, 'nelements=%s' % self.nelements
        assert self.ntotal > 0, 'ntotal=%s' % self.ntotal
        #self.names = []
        self.nelements //= self.ntimes
        self.itime = 0
        self.ielement = 0
        self.itotal = 0
        #self.ntimes = 0
        #self.nelements = 0
        self.is_built = True

        #print("ntimes=%s nelements=%s ntotal=%s" % (self.ntimes, self.nelements, self.ntotal))
        dtype = 'float32'
        if isinstance(self.nonlinear_factor, int):
            dtype = 'int32'
        self._times = zeros(self.ntimes, dtype=dtype)
        self.element = zeros(self.nelements, dtype='int32')

        # [station, bending_moment1, bending_moment2, shear1, shear2, axial, torque]
        self.data = zeros((self.ntimes, self.ntotal, 7), dtype='float32')

    def add(self, dt, data):
        self.add_sort1(dt, data)

    def add_sort1(self, dt, data):
        data = [eid, sd, bm1, bm2, ts1, ts2, af, trq] = data
        self._times[self.itime] = dt
        self.element[self.ielement] = eid

        # station, bending_moment1, bending_moment2, shear1, shear2, axial, torque
        self.data[self.itime, self.ielement, :] = [sd, bm1, bm2, ts1, ts2, af, trq]
        self.ielement += 1

    def get_stats(self):
        if not self.is_built:
            return ['<%s>\n' % self.__class__.__name__,
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
            ntimes_word = 'ntimes'
        else:
            msg.append('  type=%s nelements=%i\n'
                       % (self.__class__.__name__, nelements))
            ntimes_word = 1
        msg.append('  eType\n')
        headers = self.get_headers()
        n = len(headers)
        msg.append('  data: [%s, nnodes, %i] where %i=[%s]\n' % (ntimes_word, n, n, str(', '.join(headers))))
        msg.append('  data.shape = %s\n' % str(self.data.shape).replace('L', ''))
        #msg.append('  element type: %s\n' % self.element_type)
        msg.append('  element name: %s\n  ' % self.element_name)
        msg += self.get_data_code()
        return msg

    def eid_to_element_node_index(self, eids):
        ind = searchsorted(eids, self.element)
        return ind

    def _write_f06(self, header, page_stamp, page_num=1, f=None, is_mag_phase=False, is_sort1=True):
        words = [
            '                         F O R C E   D I S T R I B U T I O N   I N   B A R   E L E M E N T S          ( C B A R )\n'
            '0    ELEMENT  STATION         BEND-MOMENT                      SHEAR FORCE                     AXIAL\n'
            '       ID.     (PCT)     PLANE 1        PLANE 2           PLANE 1        PLANE 2               FORCE              TORQUE\n']
            # '        15893   0.000   1.998833E+02   9.004551E+01      2.316835E+00   1.461960E+00         -2.662207E+03       9.795244E-02'

        # words = ['                                 F O R C E S   I N   B A R   E L E M E N T S         ( C B A R )\n',
                 # '0    ELEMENT         BEND-MOMENT END-A            BEND-MOMENT END-B                - SHEAR -               AXIAL\n',
                 # '       ID.         PLANE 1       PLANE 2        PLANE 1       PLANE 2        PLANE 1       PLANE 2         FORCE         TORQUE\n']
        msg = []
        #header[1] = ' %s = %10.4E\n' % (self.data_code['name'], dt)
        eids = self.element
        #f.write(''.join(words))

        ntimes = self.data.shape[0]
        for itime in range(ntimes):
            dt = self._times[itime]  # TODO: rename this...
            header = _eigenvalue_header(self, header, itime, ntimes, dt)
            f.write(''.join(header + words))

            # loop over all the elements
            # sd, bm1, bm2, ts1, ts2, af, trq
            sd = self.data[itime, :, 0]
            bm1 = self.data[itime, :, 1]
            bm2 = self.data[itime, :, 2]
            ts1 = self.data[itime, :, 4]
            ts2 = self.data[itime, :, 5]
            af = self.data[itime, :, 6]
            trq = self.data[itime, :, 7]
            for eid, bm1i, bm2i, ts1i, ts2i, afi, trqi in zip(eids, sd, bm1, bm2, ts1, ts2, af, trq):
                ([bm1i, bm2i, ts1i, ts2i, afi, trqi], is_all_zeros) = writeFloats13E([
                  bm1i, bm2i, ts1i, ts2i, afi, trqi])
            f.write('     %8i  %4.3f  %-13s  %-13s     %-13s  %-13s     %-13s     %s\n' % (eid, sd, bm1i, bm2i, ts1i, ts2i, afi, trqi))
            f.write(page_stamp % page_num)
        return page_num

    def __eq__(self, table):
        assert self.is_sort1() == table.is_sort1()
        assert self.nonlinear_factor == table.nonlinear_factor
        assert self.ntotal == table.ntotal
        assert self.table_name == table.table_name, 'table_name=%r table.table_name=%r' % (self.table_name, table.table_name)
        assert self.approach_code == table.approach_code
        if not array_equal(self.element, table.element):
            assert self.element.shape == table.element.shape, 'element shape=%s table.shape=%s' % (self.element.shape, table.element.shape)
            msg = 'table_name=%r class_name=%s\n' % (self.table_name, self.__class__.__name__)
            msg += '%s\n' % str(self.code_information())
            msg += 'Eid\n'
            for eid, eid2 in zip(self.element, table.element):
                msg += '%s, %s\n' % (eid, eid2)
            print(msg)
            raise ValueError(msg)
        if not array_equal(self.data, table.data):
            msg = 'table_name=%r class_name=%s\n' % (self.table_name, self.__class__.__name__)
            msg += '%s\n' % str(self.code_information())
            i = 0
            for itime in range(self.ntimes):
                for ie, e in enumerate(self.element):
                    t1 = self.data[itime, ie, :]
                    t2 = table.data[itime, ie, :]
                    (sd1, bm11, bm21, ts11, ts21, af1, trq1) = t1
                    (sd2, bm12, bm22, ts12, ts22, af2, trq2) = t2

                    if not array_equal(t1, t2):
                        msg += '(%s)    (%s, %s, %s, %s, %s, %s, %s)  (%s, %s, %s, %s, %s, %s, %s)\n' % (
                            eid,
                            sd1, bm11, bm21, ts11, ts21, af1, trq1,
                            sd2, bm12, bm22, ts12, ts22, af2, trq2)
                        i += 1
                        if i > 10:
                            print(msg)
                            raise ValueError(msg)
                #print(msg)
                if i > 0:
                    raise ValueError(msg)
        return True


class RealCGapForce(ScalarObject):  # 38-CGAP
    def __init__(self, data_code, is_sort1, isubcase, dt):
        self.element_type = None
        self.element_name = None
        ScalarObject.__init__(self, data_code, isubcase)
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
            self.add = self.add_sort2

    def get_stats(self):
        msg = ['  '] + self.get_data_code()
        if self.dt is not None:  # transient
            ntimes = len(self.fx)
            time0 = get_key0(self.fx)
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
        self.fx[dt][eid] = fx
        self.sfy[dt][eid] = sfy
        self.sfz[dt][eid] = sfz
        self.u[dt][eid] = u
        self.v[dt][eid] = v
        self.w[dt][eid] = w
        self.sv[dt][eid] = sv
        self.sw[dt][eid] = sw

    def add_sort2(self, eid, data):
        [dt, fx, sfy, sfz, u, v, w, sv, sw] = data
        if dt not in self.fx:
            self.add_new_transient(dt)
        self.fx[dt][eid] = fx
        self.sfy[dt][eid] = sfy
        self.sfz[dt][eid] = sfz
        self.u[dt][eid] = u
        self.v[dt][eid] = v
        self.w[dt][eid] = w
        self.sv[dt][eid] = sv
        self.sw[dt][eid] = sw


class RealBendForce(ScalarObject):  # 69-CBEND
    def __init__(self, data_code, is_sort1, isubcase, dt):
        self.element_type = None
        self.element_name = None
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
        self._fillObject(dt, eid, nidA, bm1A, bm2A, sp1A, sp2A, axialA, torqueA,
                                  nidB, bm1B, bm2B, sp1B, sp2B, axialB, torqueB)

    def _fillObject(self, dt, eid, nidA, bm1A, bm2A, sp1A, sp2A, axialA, torqueA,
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


class RealPentaPressureForce(ScalarObject):  # 77-PENTA_PR,78-TETRA_PR
    def __init__(self, data_code, is_sort1, isubcase, dt):
        self.element_type = None
        self.element_name = None
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

    def write_f06(self, header, page_stamp, page_num=1, f=None, is_mag_phase=False, is_sort1=True):
        #words = ['                                   P E A K   A C C E L E R A T I O N S   A N D   P R E S S U R E S\n',
        #         ' \n',
        #         '    TIME         EL-TYPE             X-ACCELERATION            Y-ACCELERATION            Z-ACCELERATION            PRESSURE (DB)\n']
        if self.nonlinear_factor is not None:
            return self._write_f06_transient(header, page_stamp, page_num, f, is_sort1=is_sort1)
        f.write('%s write_f06 not implemented...\n' % self.__class__.__name__)
        #raise NotImplementedError()
        return page_num

    def _write_f06_transient(self, header, page_stamp, page_num=1, f=None, is_mag_phase=False, is_sort1=True):
        words = ['                                   P E A K   A C C E L E R A T I O N S   A N D   P R E S S U R E S\n',
                 ' \n',
                 '    TIME         EL-TYPE             X-ACCELERATION            Y-ACCELERATION            Z-ACCELERATION            PRESSURE (DB)\n']
        msg = []
        for dt, acc in sorted(self.acceleration.items()):
            header[1] = ' %s = %10.4E\n' % (self.data_code['name'], dt)
            msg += header + words
            for eid in sorted(acc):
                ax, ay, az = self.acceleration[dt][eid]
                vx, vy, vz = self.velocity[dt][eid]
                pressure = self.pressure[dt][eid]
                vals = [ax, ay, az, pressure]
                (vals2, is_all_zeros) = writeFloats13E(vals)
                [ax, ay, az, pressure] = vals2
                eType = 'PENPR'
                msg.append('0%13s    %5s               %-13s             %-13s             %-13s             %s\n' % (eid, eType, ax, ay, az, pressure))
            msg.append(page_stamp % page_num)
            f.write(''.join(msg))
            msg = ['']
            page_num += 1
        return page_num - 1


class RealCBushForceArray(ScalarObject):
    def __init__(self, data_code, is_sort1, isubcase, dt):
        self.element_type = None
        self.element_name = None
        ScalarObject.__init__(self, data_code, isubcase)
        #self.code = [self.format_code, self.sort_code, self.s_code]

        #self.ntimes = 0  # or frequency/mode
        #self.ntotal = 0
        self.nelements = 0  # result specific

        if is_sort1:
            self.add = self.add_sort1
        else:
            raise NotImplementedError('SORT2')

    def _reset_indices(self):
        self.itotal = 0
        self.ielement = 0

    def get_headers(self):
        headers = ['fx', 'fy', 'fz', 'mx', 'my', 'mz']
        return headers

    def build(self):
        #print('ntimes=%s nelements=%s ntotal=%s' % (self.ntimes, self.nelements, self.ntotal))
        if self.is_built:
            return

        assert self.ntimes > 0, 'ntimes=%s' % self.ntimes
        assert self.nelements > 0, 'nelements=%s' % self.nelements
        assert self.ntotal > 0, 'ntotal=%s' % self.ntotal
        #self.names = []
        self.nelements //= self.ntimes
        self.itime = 0
        self.ielement = 0
        self.itotal = 0
        #self.ntimes = 0
        #self.nelements = 0
        self.is_built = True

        #print("ntimes=%s nelements=%s ntotal=%s" % (self.ntimes, self.nelements, self.ntotal))
        dtype = 'float32'
        if isinstance(self.nonlinear_factor, int):
            dtype = 'int32'
        self._times = zeros(self.ntimes, dtype=dtype)
        self.element = zeros(self.nelements, dtype='int32')

        #[fx, fy, fz, mx, my, mz]
        self.data = zeros((self.ntimes, self.ntotal, 6), dtype='float32')

    def __eq__(self, table):
        assert self.is_sort1() == table.is_sort1()
        assert self.nonlinear_factor == table.nonlinear_factor
        assert self.ntotal == table.ntotal
        assert self.table_name == table.table_name, 'table_name=%r table.table_name=%r' % (self.table_name, table.table_name)
        assert self.approach_code == table.approach_code
        if not array_equal(self.element, table.element):
            assert self.element.shape == table.element.shape, 'shape=%s element.shape=%s' % (self.element.shape, table.element.shape)
            msg = 'table_name=%r class_name=%s\n' % (self.table_name, self.__class__.__name__)
            msg += '%s\n' % str(self.code_information())
            for eid, eid2 in zip(self.element, table.element):
                msg += '%s, %s\n' % (eid, eid2)
            print(msg)
            raise ValueError(msg)
        if not array_equal(self.data, table.data):
            msg = 'table_name=%r class_name=%s\n' % (self.table_name, self.__class__.__name__)
            msg += '%s\n' % str(self.code_information())
            ntimes = self.data.shape[0]

            i = 0
            if self.is_sort1():
                for itime in range(ntimes):
                    for ieid, eid, in enumerate(self.element):
                        t1 = self.data[itime, ieid, :]
                        t2 = table.data[itime, ieid, :]
                        (fx1, fy1, fz1, mx1, my1, mz1) = t1
                        (fx2, fy2, fz2, mx2, my2, mz2) = t2
                        if not allclose(t1, t2):
                        #if not array_equal(t1, t2):
                            msg += '%s\n  (%s, %s, %s)\n  (%s, %s, %s)\n' % (
                                eid,
                                fx1, fy1, fz1, mx1, my1, mz1,
                                fx2, fy2, fz2, mx2, my2, mz2)
                            i += 1
                        if i > 10:
                            print(msg)
                            raise ValueError(msg)
            else:
                raise NotImplementedError(self.is_sort2())
            if i > 0:
                print(msg)
                raise ValueError(msg)
        return True

    def add(self, dt, eid, fx, fy, fz, mx, my, mz):
        self.add_sort1(dt, eid, fx, fy, fz, mx, my, mz)

    def add_sort1(self, dt, eid, fx, fy, fz, mx, my, mz):
        self._times[self.itime] = dt
        self.element[self.ielement] = eid
        self.data[self.itime, self.ielement, :] = [fx, fy, fz, mx, my, mz]
        self.ielement += 1

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
            ntimes_word = 'ntimes'
        else:
            msg.append('  type=%s nelements=%i\n'
                       % (self.__class__.__name__, nelements))
            ntimes_word = 1
        msg.append('  eType\n')
        headers = self.get_headers()
        n = len(headers)
        msg.append('  data: [%s, nnodes, %i] where %i=[%s]\n' % (ntimes_word, n, n, str(', '.join(headers))))
        msg.append('  data.shape = %s\n' % str(self.data.shape).replace('L', ''))
        #msg.append('  element type: %s\n' % self.element_type)
        msg.append('  element name: %s\n  ' % self.element_name)
        msg += self.get_data_code()
        return msg

    def get_f06_header(self, is_mag_phase=True):
        # TODO: not up to date...
        msg = ['                                     F O R C E S   I N   B U S H   E L E M E N T S      ( C B U S H)\n', ]
        return msg

    def get_element_index(self, eids):
        # elements are always sorted; nodes are not
        itot = searchsorted(eids, self.element)  #[0]
        return itot

    def eid_to_element_node_index(self, eids):
        #ind = ravel([searchsorted(self.element == eid) for eid in eids])
        ind = searchsorted(eids, self.element)
        #ind = ind.reshape(ind.size)
        #ind.sort()
        return ind

    def write_f06(self, header, page_stamp, page_num=1, f=None, is_mag_phase=False, is_sort1=True):
        msg_temp = self.get_f06_header(is_mag_phase)

        # write the f06
        #(ntimes, ntotal, two) = self.data.shape
        ntimes = self.data.shape[0]

        eids = self.element
        is_odd = False
        nwrite = len(eids)
        if len(eids) % 2 == 1:
            nwrite -= 1
            is_odd = True

        #print('len(eids)=%s nwrite=%s is_odd=%s' % (len(eids), nwrite, is_odd))
        for itime in range(ntimes):
            dt = self._times[itime]  # TODO: rename this...
            header = _eigenvalue_header(self, header, itime, ntimes, dt)
            f.write(''.join(header + msg_temp))

            # TODO: can I get this without a reshape?
            #print("self.data.shape=%s itime=%s ieids=%s" % (str(self.data.shape), itime, str(ieids)))
            fx = self.data[itime, :, 0]
            fy = self.data[itime, :, 1]
            fz = self.data[itime, :, 2]
            mx = self.data[itime, :, 3]
            my = self.data[itime, :, 4]
            mz = self.data[itime, :, 5]

            # loop over all the elements
            out = []
            for eid, fxi, fyi, fzi, mxi, myi, mzi in zip(eids, fx, fy, fz, mx, my, mz):
                ([fxi, fyi, fzi, mxi, myi, mzi], is_all_zeros) = writeFloats13E([fxi, fyi, fzi, mxi, myi, mzi])
                f.write('      %8i   %-13s  %-13s  %-13s  %-13s  %-13s  %s\n' % tuple(eid, fxi, fyi, fzi, mxi, myi, mzi))
            f.write(page_stamp % page_num)
            page_num += 1
        return page_num - 1


class RealCBushForce(ScalarObject):  # 102-CBUSH
    def __init__(self, data_code, is_sort1, isubcase, dt):
        self.element_type = None
        self.element_name = None
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

    def add(self, dt, eid, fx, fy, fz, mx, my, mz):
        self.force[eid] = array([fx, fy, fz], dtype='float32')
        self.moment[eid] = array([mx, my, mz], dtype='float32')

    def add_sort1(self, dt, eid, fx, fy, fz, mx, my, mz):
        if dt not in self.force:
            self.add_new_transient(dt)
        self.force[dt][eid] = array([fx, fy, fz], dtype='float32')
        self.moment[dt][eid] = array([mx, my, mz], dtype='float32')

    def add_sort2(self, eid, dt, fx, fy, fz, mx, my, mz):
        if dt not in self.force:
            self.add_new_transient(dt)
        self.force[dt][eid] = array([fx, fy, fz], dtype='float32')
        self.moment[dt][eid] = array([mx, my, mz], dtype='float32')


class RealForce_VU(ScalarObject):  # 191-VUBEAM
    def __init__(self, data_code, is_sort1, isubcase, dt):
        self.element_type = None
        self.element_name = None
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

        # handle SORT1 case
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

    def add(self, nNodes, dt, data):
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

    def add_sort1(self, nNodes, dt, data):
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

    def add_sort2(self, nNodes, eid, data):
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


class RealForce_VU_2D(ScalarObject):  # 190-VUTRIA # 189-VUQUAD
    def __init__(self, data_code, is_sort1, isubcase, dt):
        self.element_type = None
        self.element_name = None
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

        # handle SORT1 case
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

    def add(self, nNodes, dt, data):
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

    def add_sort1(self, nNodes, dt, data):
        [eid, parent, coord, icord, theta, forces] = data
        self._fillObject(dt, eid, parent, coord, icord, theta, forces)

    def add_sort2(self, nNodes, eid, data):
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
