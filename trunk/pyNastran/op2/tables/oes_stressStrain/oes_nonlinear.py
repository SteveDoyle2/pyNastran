from __future__ import (nested_scopes, generators, division, absolute_import,
                        print_function, unicode_literals)
from six import iteritems
from six.moves import range
from math import isnan

from pyNastran.op2.tables.oes_stressStrain.real.oes_objects import StressObject, StrainObject
from pyNastran.f06.f06_formatting import writeFloats13E


class NonlinearQuad(StressObject):

    def __init__(self, data_code, is_sort1, isubcase, dt):
        StressObject.__init__(self, data_code, isubcase)
        #self.eType = 'QUAD4FD' # or CTRIA3

        self.code = [self.format_code, self.sort_code, self.s_code]
        self.eType = {}
        self.fiberDistance = {}
        self.oxx = {}
        self.oyy = {}
        self.ozz = {}
        self.txy = {}

        self.exx = {}
        self.eyy = {}
        self.ezz = {}
        self.exy = {}

        self.es = {}
        self.eps = {}
        self.ecs = {}

        self.dt = dt
        if is_sort1:
            if dt is not None:
                self.add = self.add_sort1
                self.add_new_eid = self.add_new_eid_sort1
        else:
            assert dt is not None
            self.add = self.addSort2
            self.add_new_eid = self.add_new_eid_sort2

    def get_stats(self):
        nelements = len(self.eType)

        msg = self.get_data_code()
        if self.nonlinear_factor is not None:  # transient
            ntimes = len(self.oxx)
            msg.append('  type=%s ntimes=%s nelements=%s\n'
                       % (self.__class__.__name__, ntimes, nelements))
        else:
            msg.append('  type=%s nelements=%s\n' % (self.__class__.__name__,
                                                     nelements))
        msg.append('  eType, fiberDistance, oxx, oyy, ozz, txy, '
                   'exx, eyy, ezz, exy, es, eps, ecs\n')
        return msg

    def delete_transient(self, dt):
        del self.fiberDistance[dt]
        del self.oxx[dt]
        del self.oyy[dt]
        del self.ozz[dt]
        del self.txy[dt]

        del self.exx[dt]
        del self.eyy[dt]
        del self.ezz[dt]
        del self.exy[dt]

        del self.es[dt]
        del self.eps[dt]
        del self.ecs[dt]

    def get_transients(self):
        k = self.oxx.keys()
        k.sort()
        return k

    def add_new_transient(self, dt):
        self.fiberDistance[dt] = {}
        self.oxx[dt] = {}
        self.oyy[dt] = {}
        self.ozz[dt] = {}
        self.txy[dt] = {}

        self.exx[dt] = {}
        self.eyy[dt] = {}
        self.ezz[dt] = {}
        self.exy[dt] = {}

        self.es[dt] = {}
        self.eps[dt] = {}
        self.ecs[dt] = {}

    def add_new_eid_sort1(self, eType, dt, data):
        if dt not in self.oxx:
            self.add_new_transient(dt)
        (eid, fd, sx, sy, sz, txy, es, eps, ecs, ex, ey, ez, exy) = data
        self.fiberDistance[dt][eid] = [fd]
        if isnan(sz):
            sz = 0.
        if isnan(ez):
            ez = 0.
        self.eType[eid] = eType
        self.oxx[dt][eid] = [sx]
        self.oyy[dt][eid] = [sy]
        self.ozz[dt][eid] = [sz]
        self.txy[dt][eid] = [txy]

        self.exx[dt][eid] = [ex]
        self.eyy[dt][eid] = [ey]
        self.ezz[dt][eid] = [ez]
        self.exy[dt][eid] = [exy]

        self.es[dt][eid] = [es]
        self.eps[dt][eid] = [eps]
        self.ecs[dt][eid] = [ecs]

    def add_sort1(self, dt, data):
        (eid, fd, sx, sy, sz, txy, es, eps, ecs, ex, ey, ez, exy) = data
        self.fiberDistance[dt][eid].append(fd)
        if isnan(sz):
            sz = 0.
        if isnan(ez):
            ez = 0.

        self.oxx[dt][eid].append(sx)
        self.oyy[dt][eid].append(sy)
        self.ozz[dt][eid].append(sz)
        self.txy[dt][eid].append(txy)

        self.exx[dt][eid].append(ex)
        self.eyy[dt][eid].append(ey)
        self.ezz[dt][eid].append(ez)
        self.exy[dt][eid].append(exy)

        self.es[dt][eid].append(es)
        self.eps[dt][eid].append(eps)
        self.ecs[dt][eid].append(ecs)

    def write_f06(self, header, pageStamp, page_num=1, f=None, is_mag_phase=None):
        msgStart = ['      ELEMENT-ID =     129\n'
                    '               N O N L I N E A R   S T R E S S E S   I N   Q U A D R I L A T E R A L   E L E M E N T S    ( Q U A D 4 )\n'
                    ' \n',
                    '    TIME         FIBER                        STRESSES/ TOTAL STRAINS                     EQUIVALENT    EFF. STRAIN     EFF. CREEP\n'
                    '               DISTANCE           X              Y             Z               XY           STRESS    PLASTIC/NLELAST     STRAIN\n']
#0 5.000E-05  -5.000000E-01  -4.484895E+01  -1.561594E+02                 -2.008336E-02   1.392609E+02   0.0            0.0
        msgE = {}
        msgT = {}
        for (dt, Oxxs) in sorted(iteritems(self.oxx)):
            header[1] = ' %s = %10.4E\n' % (self.data_code['name'], dt)

            for (eid, oxxs) in sorted(iteritems(Oxxs)):
                msgE[eid] = header + ['      ELEMENT-ID = %8i\n' % (eid)]
                if eid not in msgT:
                    msgT[eid] = []
                for i, oxx in enumerate(oxxs):
                    fd = self.fiberDistance[dt][eid][i]
                    oxx = self.oxx[dt][eid][i]
                    oyy = self.oyy[dt][eid][i]
                    ozz = self.ozz[dt][eid][i]
                    txy = self.txy[dt][eid][i]

                    exx = self.exx[dt][eid][i]
                    eyy = self.eyy[dt][eid][i]
                    ezz = self.ezz[dt][eid][i]
                    exy = self.exy[dt][eid][i]

                    es = self.es[dt][eid][i]
                    eps = self.eps[dt][eid][i]
                    ecs = self.ecs[dt][eid][i]
                    ([oxx, oyy, ozz, txy, exx, eyy, es, eps, ecs, exx, eyy, ezz, exy], is_all_zeros) = writeFloats13E([oxx, oyy, ozz, txy, exx, eyy, es, eps, ecs, exx, eyy, ezz, exy])
                    if i == 0:
                        msgT[eid].append('0 %9.3E %-13s  %-13s  %-13s  %-13s  %-13s  %-13s  %-13s  %s\n' % (dt, fd, oxx, oyy, ozz, txy, es, eps, ecs))
                    else:
                        msgT[eid].append('     %9s %-13s  %-13s  %-13s  %-13s  %-13s\n' % ('', '', exx, eyy, ezz, exy))

        msg = []
        for eid, e in sorted(iteritems(msgE)):
            msg += header + e + msgStart + msgT[eid]
            msg.append(pageStamp % page_num)
            page_num += 1
        f.write(''.join(msg))
        return page_num - 1


class HyperelasticQuad(StressObject):
    def __init__(self, data_code, is_sort1, isubcase, dt):
        StressObject.__init__(self, data_code, isubcase)
        self.eType = 'QUAD4FD'

        self.code = [self.format_code, self.sort_code, self.s_code]
        self.Type = {}
        self.IDs = {}
        self.oxx = {}
        self.oyy = {}
        self.txy = {}
        self.angle = {}
        self.majorP = {}
        self.minorP = {}

        self.dt = dt
        if is_sort1:
            if dt is not None:
                self.add = self.add_sort1
                self.add_new_eid = self.add_new_eid_sort1
        else:
            assert dt is not None
            self.add = self.addSort2
            self.add_new_eid = self.add_new_eid_sort2

    def get_stats(self):
        nelements = len(self.eType)

        msg = self.get_data_code()
        if self.nonlinear_factor is not None:  # transient
            ntimes = len(self.oxx)
            msg.append('  type=%s ntimes=%s nelements=%s\n'
                       % (self.__class__.__name__, ntimes, nelements))
        else:
            msg.append('  type=%s nelements=%s\n' % (self.__class__.__name__,
                                                     nelements))
        msg.append('  Type, oxx, oyy, txy, angle, majorP, minorP\n')
        return msg

    def delete_transient(self, dt):
        del self.fiberDistance[dt]
        del self.oxx[dt]
        del self.oyy[dt]
        del self.txy[dt]

        del self.angle[dt]
        del self.majorP[dt]
        del self.minorP[dt]

    def get_transients(self):
        k = self.oxx.keys()
        k.sort()
        return k

    def add_new_transient(self, dt):
        self.oxx[dt] = {}
        self.oyy[dt] = {}
        self.txy[dt] = {}
        self.angle[dt] = {}
        self.majorP[dt] = {}
        self.minorP[dt] = {}

    def add_new_eid_sort1(self, dt, data):
        if dt not in self.oxx:
            self.add_new_transient(dt)
        (eid, Type, oxx, oyy, txy, angle, majorP, minorP) = data
        self.Type[eid] = Type
        self.oxx[dt] = {eid: [oxx]}
        self.oyy[dt] = {eid: [oyy]}
        self.txy[dt] = {eid: [txy]}
        self.angle[dt] = {eid: [angle]}
        self.majorP[dt] = {eid: [majorP]}
        self.minorP[dt] = {eid: [minorP]}

    def add_sort1(self, dt, eid, data):
        (ID, oxx, oyy, txy, angle, majorP, minorP) = data
        self.oxx[dt][eid].append(oxx)
        self.oyy[dt][eid].append(oyy)
        self.txy[dt][eid].append(txy)
        self.angle[dt][eid].append(angle)
        self.majorP[dt][eid].append(majorP)
        self.minorP[dt][eid].append(minorP)

    def write_f06(self, header, pageStamp, page_num=1, f=None, is_mag_phase=False):  # .. todo:: doesnt support CTRIA3NL (calls them CQUAD4s)
        msg = ['           S T R E S S E S   I N   H Y P E R E L A S T I C   Q U A D R I L A T E R A L   E L E M E N T S  ( QUAD4FD )\n',
               '  ELEMENT     GRID/    POINT       ---------CAUCHY STRESSES--------             PRINCIPAL STRESSES (ZERO SHEAR)\n',
               '     ID       GAUSS      ID      NORMAL-X       NORMAL-Y      SHEAR-XY       ANGLE         MAJOR           MINOR\n', ]
               #0       1     GAUS         1   7.318995E+00   6.367099E-01  -6.551054E+00   -31.4888    1.133173E+01   -3.376026E+00
               #                           2   1.097933E+01   4.149028E+00   6.278160E+00    30.7275    1.471111E+01    4.172537E-01

        for dt, Oxxs in sorted(iteritems(self.oxx)):
            #header[-1] = '     LOAD STEP = %12.5E' %(dt)
            msg += header
            for eid, oxxs in sorted(iteritems(Oxxs)):
                gauss = self.Type[eid]
                oxx = self.oxx[dt][eid]
                oyy = self.oyy[dt][eid]
                txy = self.txy[dt][eid]
                angle = self.angle[dt][eid]
                majorP = self.majorP[dt][eid]
                minorP = self.minorP[dt][eid]

                for i in range(4):  # 1,2,3,4
                    if i == 0:
                        msg.append('0%8i %8s  %8i  %13E.6  %13E.6  %13E.6  %13E.6  %13E.6  %13E.6\n' % (eid, gauss, i + 1, oxx[i], oyy[i], txy[i], angle[i], majorP[i], minorP[i]))
                    else:
                        msg.append(' %8s %8s  %8i  %13E.6  %13E.6  %13E.6  %13E.6  %13E.6  %13E.6\n' % ('', '', i + 1, oxx[i], oyy[i], txy[i], angle[i], majorP[i], minorP[i]))
        f.write(''.join(msg))
        return page_num


class NonlinearRod(StressObject):
    def __init__(self, data_code, is_sort1, isubcase, dt):
        StressObject.__init__(self, data_code, isubcase)
        #self.eType = 'CROD'
        self.eTypeMap = {89: 'CRODNL', 92: 'CONRODNL'}
        self.code = [self.format_code, self.sort_code, self.s_code]

        self.eType = {}
        self.axialStress = {}
        self.equivStress = {}
        self.totalStrain = {}
        self.effectivePlasticCreepStrain = {}
        self.effectiveCreepStrain = {}
        self.linearTorsionalStress = {}

        self.dt = dt
        if is_sort1:
            if dt is not None:
                self.add = self.add_sort1
                #self.add_new_eid = self.add_new_eid_sort1
        else:
            assert dt is not None
            self.add = self.addSort2
            #self.add_new_eid = self.add_new_eid_sort2

    def get_stats(self):
        nelements = len(self.eType)
        msg = self.get_data_code()
        if self.nonlinear_factor is not None:  # transient
            ntimes = len(self.axialStress)
            msg.append('  type=%s ntimes=%s nelements=%s\n'
                       % (self.__class__.__name__, ntimes, nelements))
        else:
            msg.append('  type=%s nelements=%s\n' % (self.__class__.__name__,
                                                     nelements))
        msg.append('  eType, axialStress, equivStress, totalStrain, '
                   'effectivePlasticCreepStrain, effectiveCreepStrain, '
                   'linearTorsionalStress\n')
        return msg

    def delete_transient(self, dt):
        del self.axialStress[dt]
        del self.equivStress[dt]
        del self.totalStrain[dt]
        del self.effectivePlasticCreepStrain[dt]

        del self.effectiveCreepStrain[dt]
        del self.linearTorsionalStress[dt]

    def get_transients(self):
        k = self.axialStress.keys()
        k.sort()
        return k

    def add_new_transient(self, dt):
        self.axialStress[dt] = {}
        self.equivStress[dt] = {}
        self.totalStrain[dt] = {}
        self.effectivePlasticCreepStrain[dt] = {}
        self.effectiveCreepStrain[dt] = {}
        self.linearTorsionalStress[dt] = {}

    def add_sort1(self, eType, dt, data):
        if dt not in self.axialStress:
            self.add_new_transient(dt)
        eid = data[0]
        self.eType[eid] = eType
        self.axialStress[dt][eid] = data[1]
        self.equivStress[dt][eid] = data[2]
        self.totalStrain[dt][eid] = data[3]
        self.effectivePlasticCreepStrain[dt][eid] = data[4]
        self.effectiveCreepStrain[dt][eid] = data[5]
        self.linearTorsionalStress[dt][eid] = data[6]
        #print data

    def write_f06(self, header, pageStamp, page_num=1, f=None, is_mag_phase=False):  # .. todo:: doesnt support CONROD/CTUBE (calls them CRODs)
        """
        ::

          ELEMENT-ID =     102
                                   N O N L I N E A R   S T R E S S E S   I N   R O D   E L E M E N T S      ( C R O D )
            TIME          AXIAL STRESS         EQUIVALENT         TOTAL STRAIN       EFF. STRAIN          EFF. CREEP        LIN. TORSIONAL
                                                 STRESS                             PLASTIC/NLELAST          STRAIN              STRESS
          2.000E-02        1.941367E+01        1.941367E+01        1.941367E-04        0.0                 0.0                 0.0
          3.000E-02        1.941367E+01        1.941367E+01        1.941367E-04        0.0                 0.0                 0.0
        """
        msg = []
        msgStart = ['                         N O N L I N E A R   S T R E S S E S   I N   R O D   E L E M E N T S      ( C R O D )\n',
                    ' \n',
                    '    TIME          AXIAL STRESS         EQUIVALENT         TOTAL STRAIN       EFF. STRAIN          EFF. CREEP        LIN. TORSIONAL\n',
                    '                                         STRESS                             PLASTIC/NLELAST          STRAIN              STRESS\n']
        msgE = {}
        msgT = {}
        for dt, axials in sorted(iteritems(self.axialStress)):
            for eid, axial in sorted(iteritems(axials)):
                eqs = self.equivStress[dt][eid]
                ts = self.totalStrain[dt][eid]
                epcs = self.effectivePlasticCreepStrain[dt][eid]
                ecs = self.effectiveCreepStrain[dt][eid]
                lts = self.linearTorsionalStress[dt][eid]
                #print "dt=%s axials=%s eqs=%s ts=%s epcs=%s ecs=%s lts=%s" %(dt,axial,eqs,ts,epcs,ecs,lts)
                msgE[eid] = '      ELEMENT-ID = %8i\n' % (eid)
                if eid not in msgT:
                    msgT[eid] = []
                msgT[eid].append('  %9.3E       %13.6E       %13.6E       %13.6E       %13.6E       %13.6E       %13.6E\n' % (dt, axial, eqs, ts, epcs, ecs, lts))

        for eid, e in sorted(iteritems(msgE)):
            msg += header + [e] + msgStart + msgT[eid]
            msg.append(pageStamp % page_num)
        f.write(''.join(msg))
        return page_num
