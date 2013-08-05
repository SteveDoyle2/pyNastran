from __future__ import (nested_scopes, generators, division, absolute_import,
                        print_function, unicode_literals)

from .oes_objects import StressObject, StrainObject
from pyNastran.f06.f06_formatting import writeFloats12E


class CompositePlateStressObject(StressObject):
    """
    ::

      # s_code = 0
                      S T R E S S E S   I N   L A Y E R E D   C O M P O S I T E   E L E M E N T S   ( Q U A D 4 )
      ELEMENT  PLY  STRESSES IN FIBER AND MATRIX DIRECTIONS    INTER-LAMINAR  STRESSES  PRINCIPAL STRESSES (ZERO SHEAR)      MAX
        ID      ID    NORMAL-1     NORMAL-2     SHEAR-12     SHEAR XZ-MAT  SHEAR YZ-MAT  ANGLE    MAJOR        MINOR        SHEAR
    """
    def __init__(self, data_code, is_sort1, isubcase, dt, read_mode):
        StressObject.__init__(self, data_code, isubcase, read_mode)
        self.eType = {}

        self.code = [self.format_code, self.sort_code, self.s_code]
        self.o11 = {}
        self.o22 = {}
        self.t12 = {}
        self.t1z = {}
        self.t2z = {}
        self.angle = {}
        self.majorP = {}
        self.minorP = {}
        #self.fiberCurvature = {}
        self.ovmShear = {}
        #print "self.data_code = ",self.data_code
        #if self.isVonMisesStress():
        #if self.code == [1, 0, 0]:
            #assert not self.isVonMises, 'isVonMises=%s' % (self.isVonMises())
            ##self.isVonMises = True
        #else:
            #raise RuntimeError("Invalid Code: compositePlateStress - get the format/sort/stressCode=%s" % (self.code))

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
            ntimes = len(self.o11)
            msg.append('  type=%s ntimes=%s nelements=%s\n'
                       % (self.__class__.__name__, ntimes, nelements))
        else:
            msg.append('  type=%s nelements=%s\n' % (self.__class__.__name__,
                                                     nelements))
        msg.append('  eType, fiberCurvature, o11, o22, t12, t1z, t2z, angle, '
                   'majorP, minorP, ovmShear\n')
        return msg

    def add_f06_data(self, data, transient, eType):
        """
        ::
                         S T R E S S E S   I N   L A Y E R E D   C O M P O S I T E   E L E M E N T S   ( Q U A D 4 )
          ELEMENT  PLY  STRESSES IN FIBER AND MATRIX DIRECTIONS    INTER-LAMINAR  STRESSES  PRINCIPAL STRESSES (ZERO SHEAR)      MAX
            ID      ID    NORMAL-1     NORMAL-2     SHEAR-12     SHEAR XZ-MAT  SHEAR YZ-MAT  ANGLE    MAJOR        MINOR        SHEAR
              151    1  -1.02406E+04  4.18348E+05  4.14359E+02   -8.62021E+00  1.86352E+04   89.94  4.18348E+05 -1.02410E+04  2.14295E+05
        """
        if transient is None:
            #for line in data:
                #print line
            #sys.exit()
            for line in data:
                #print line
                if eType == 'CQUAD4':
                    (eid, iLayer, o11, o22, t12, t1z, t2z, angle,
                        majorP, minorP, ovmShear) = line
                    #if nid=='CEN/4': nid='C'
                    if eid not in self.eType:
                        self.eType[eid] = 'CQUAD4'
                        self.o11[eid] = [o11]
                        self.o22[eid] = [o22]
                        self.t12[eid] = [t12]
                        self.t1z[eid] = [t1z]
                        self.t2z[eid] = [t2z]
                        self.angle[eid] = [angle]
                        self.majorP[eid] = [majorP]
                        self.minorP[eid] = [minorP]
                        self.ovmShear[eid] = [ovmShear]
                    else:
                        self.o11[eid].append(o11)
                        self.o22[eid].append(o22)
                        self.t12[eid].append(t12)
                        self.t1z[eid].append(t1z)
                        self.t2z[eid].append(t2z)
                        self.angle[eid].append(angle)
                        self.majorP[eid].append(majorP)
                        self.minorP[eid].append(minorP)
                        self.ovmShear[eid].append(ovmShear)
                else:
                    msg = 'line=%s not supported...' % (line)
                    raise NotImplementedError(msg)
            return
        #for line in data:
            #print line
        raise NotImplementedError('transient results not supported')

    def get_transients(self):
        k = self.o11.keys()
        k.sort()
        return k

    def add_new_transient(self, dt):
        """
        initializes the transient variables
        """
        #self.fiberDistance[dt] = {}
        self.o11[dt] = {}
        self.o22[dt] = {}
        self.t12[dt] = {}
        self.t1z[dt] = {}
        self.t2z[dt] = {}
        self.angle[dt] = {}
        self.majorP[dt] = {}
        self.minorP[dt] = {}
        self.ovmShear[dt] = {}

    def add_new_eid(self, eType, dt, eid, o11, o22, t12, t1z, t2z, angle,
                  majorP, minorP, ovm):
        """all points are located at the centroid"""
        #print "Composite Plate Strain add..."
        #msg = "eid=%s eType=%s o11=%g o22=%g t12=%g t1z=%g t2z=%g \nangle=%g major=%g minor=%g vm=%g" % (eid, eType, o11, o22, t12, t1z, t2z, angle, majorP, minorP, ovm)
        if eid in self.o11:
            return self.add(dt, eid, o11, o22, t12, t1z, t2z, angle, majorP, minorP, ovm)
        assert eid not in self.o11, msg + '\n  o11=%s eType=%s code=%s' % (
            self.o11[eid], self.eType[eid], self.data_code)

        self.eType[eid] = eType
        self.o11[eid] = [o11]
        self.o22[eid] = [o22]
        self.t12[eid] = [t12]
        self.t1z[eid] = [t1z]
        self.t2z[eid] = [t2z]
        self.angle[eid] = [angle]
        self.majorP[eid] = [majorP]
        self.minorP[eid] = [minorP]
        self.ovmShear[eid] = [ovm]
        #print msg
        #if nodeID==0: raise Exception(msg)

    def add_new_eid_sort1(self, eType, dt, eid, o11, o22, t12, t1z, t2z, angle,
                       majorP, minorP, ovm):
        """all points are located at the centroid"""
        #print "Composite Plate Strain add..."

        if dt not in self.o11:
            self.add_new_transient(dt)
        if eid in self.o11[dt]:
            return self.add(dt, eid, o11, o22, t12, t1z, t2z, angle,
                            majorP, minorP, ovm)

        assert eid not in self.o11[dt]
        self.eType[eid] = eType
        self.o11[dt][eid] = [o11]
        self.o22[dt][eid] = [o22]
        self.t12[dt][eid] = [t12]
        self.t1z[dt][eid] = [t1z]
        self.t2z[dt][eid] = [t2z]
        self.angle[dt][eid] = [angle]
        self.majorP[dt][eid] = [majorP]
        self.minorP[dt][eid] = [minorP]
        self.ovmShear[dt][eid] = [ovm]
        #msg = "eid=%s o11=%g o22=%g t12=%g t1z=%g t2z=%g \nangle=%g major=%g minor=%g vm=%g" % (eid, o11, o22, t12, t1z, t2z, angle, majorP, minorP, ovm)
        #print msg
        #if nodeID==0: raise Exception(msg)

    def add(self, dt, eid, o11, o22, t12, t1z, t2z, angle,
            majorP, minorP, ovm):
        #print "***add"
        #msg = "eid=%s o11=%g o22=%g t12=%g t1z=%g t2z=%g \nangle=%g major=%g minor=%g vm=%g" % (eid, o11, o22, t12, t1z, t2z, angle, majorP, minorP, ovm)
        #print msg
        #print self.o11
        self.o11[eid].append(o11)
        self.o22[eid].append(o22)
        self.t12[eid].append(t12)
        self.t1z[eid].append(t1z)
        self.t2z[eid].append(t2z)
        self.angle[eid].append(angle)
        self.majorP[eid].append(majorP)
        self.minorP[eid].append(minorP)
        self.ovmShear[eid].append(ovm)
        #if nodeID==0: raise Exception(msg)

    def add_sort1(self, dt, eid, o11, o22, t12, t1z, t2z, angle,
                 majorP, minorP, ovm):
        #print "***add"
        msg = "eid=%s o11=%g o22=%g t12=%g t1z=%g t2z=%g \nangle=%g major=%g minor=%g vm=%g" % (eid, o11, o22, t12, t1z, t2z, angle, majorP, minorP, ovm)
        #print msg
        #print self.o11

        self.o11[dt][eid].append(o11)
        self.o22[dt][eid].append(o22)
        self.t12[dt][eid].append(t12)
        self.t1z[dt][eid].append(t1z)
        self.t2z[dt][eid].append(t2z)
        self.angle[dt][eid].append(angle)
        self.majorP[dt][eid].append(majorP)
        self.minorP[dt][eid].append(minorP)
        self.ovmShear[dt][eid].append(ovm)
        #if nodeID==0: raise Exception(msg)

    def getHeaders(self):
        headers = ['o11', 'o22', 't12', 't1z', 't2z']
        if self.isVonMises:
            headers.append('oVonMises')
        else:
            headers.append('maxShear')
        return headers

    def write_f06(self, header, pageStamp, f, pageNum=1, is_mag_phase=False):
        if self.nonlinear_factor is not None:
            return self._write_f06_transient(header, pageStamp, f, pageNum)

        if self.isVonMises():
            von = 'VON'
            mises = 'MISES'
        else:
            von = 'MAX'
            mises = 'SHEAR'

        words = ['   ELEMENT  PLY  STRESSES IN FIBER AND MATRIX DIRECTIONS    INTER-LAMINAR  STRESSES  PRINCIPAL STRESSES (ZERO SHEAR)      %s\n' % (von),
                 '     ID      ID    NORMAL-1     NORMAL-2     SHEAR-12     SHEAR XZ-MAT  SHEAR YZ-MAT  ANGLE    MAJOR        MINOR        %s\n' % (mises)]

        eTypes = self.eType.values()
        if 'CQUAD4' in eTypes or 'QUAD4LC' in eTypes:
            quadMsg = header + ['                   S T R E S S E S   I N   L A Y E R E D   C O M P O S I T E   E L E M E N T S   ( Q U A D 4 )\n'] + words
            isQuad = True
        else:
            quadMsg = []
            isQuad = False

        if 'CTRIA3' in eTypes or 'TRIA3LC' in eTypes:
            isTri = True
            triMsg = header + ['                   S T R E S S E S   I N   L A Y E R E D   C O M P O S I T E   E L E M E N T S   ( T R I A 3 )\n'] + words
        else:
            isTri = False
            triMsg = []

        for eid, o11s in sorted(self.o11.iteritems()):
            out = ''
            eType = self.eType[eid]
            for iLayer in xrange(len(o11s)):
                o11 = self.o11[eid][iLayer]
                o22 = self.o22[eid][iLayer]
                t12 = self.t12[eid][iLayer]
                t1z = self.t1z[eid][iLayer]
                t2z = self.t2z[eid][iLayer]

                angle = self.angle[eid][iLayer]
                major = self.majorP[eid][iLayer]
                minor = self.minorP[eid][iLayer]
                ovm = self.ovmShear[eid][iLayer]
                (vals2, isAllZeros) = writeFloats12E([o11, o22, t12, t1z, t2z, major, minor, ovm])
                [o11, o22, t12, t1z, t2z, major, minor, ovm] = vals2
                out += '0 %8s %4s  %12s %12s %12s   %12s %12s  %6.2F %12s %12s %-s\n' % (eid, iLayer + 1, o11, o22, t12, t1z, t2z, angle, major, minor, ovm)

            if eType in ['CQUAD4', 'QUAD4LC']:
                quadMsg.append(out)
            elif eType in ['CTRIA3', 'TRIA3LC']:
                triMsg.append(out)
            #else:
            #    raise NotImplementedError('eType = |%r|' %(eType)) # CQUAD8LC

        if isQuad:
            quadMsg.append(pageStamp + str(pageNum) + '\n')
            pageNum += 1
        if isTri:
            triMsg.append(pageStamp + str(pageNum) + '\n')
            pageNum += 1

        msg = ''.join(quadMsg + triMsg)
        f.write(''.join(msg))
        return pageNum

    def _write_f06_transient(self, header, pageStamp, f,
                             pageNum=1, is_mag_phase=False):
        if self.isVonMises():
            von = 'VON'
            mises = 'MISES'
        else:
            von = 'MAX'
            mises = 'SHEAR'

        words = ['   ELEMENT  PLY  STRESSES IN FIBER AND MATRIX DIRECTIONS    INTER-LAMINAR  STRESSES  PRINCIPAL STRESSES (ZERO SHEAR)      %s\n' % (von),
                 '     ID      ID    NORMAL-1     NORMAL-2     SHEAR-12     SHEAR XZ-MAT  SHEAR YZ-MAT  ANGLE    MAJOR        MINOR        %s\n' % (mises)]

        eTypes = self.eType.values()
        if 'CQUAD4' in eTypes or 'QUAD4LC' in eTypes:
            quadWords = ['                   S T R E S S E S   I N   L A Y E R E D   C O M P O S I T E   E L E M E N T S   ( Q U A D 4 )\n'] + words
            isQuad = True
        else:
            quadWords = []
            isQuad = False

        if 'CTRIA3' in eTypes or 'TRIA3LC' in eTypes:
            isTri = True
            triWords = ['                   S T R E S S E S   I N   L A Y E R E D   C O M P O S I T E   E L E M E N T S   ( T R I A 3 )\n'] + words
        else:
            isTri = False
            triWords = []

        msg = []
        for dt, O11s in sorted(self.o11.iteritems()):
            quadMsg = []
            triMsg = []
            header[1] = ' %s = %10.4E\n' % (self.data_code['name'], dt)
            if isQuad:
                quadMsg = header + quadWords
            if isTri:
                triMsg = header + triWords

            for eid, o11s in sorted(O11s.iteritems()):
                out = ''
                eType = self.eType[eid]
                for iLayer in xrange(len(o11s)):
                    o11 = self.o11[dt][eid][iLayer]
                    o22 = self.o22[dt][eid][iLayer]
                    t12 = self.t12[dt][eid][iLayer]
                    t1z = self.t1z[dt][eid][iLayer]
                    t2z = self.t2z[dt][eid][iLayer]

                    angle = self.angle[dt][eid][iLayer]
                    major = self.majorP[dt][eid][iLayer]
                    minor = self.minorP[dt][eid][iLayer]
                    ovm = self.ovmShear[dt][eid][iLayer]
                    (vals2, isAllZeros) = writeFloats12E([o11, o22,
                                                          t12, t1z, t2z,
                                                          major, minor, ovm])
                    [o11, o22, t12, t1z, t2z, major, minor, ovm] = vals2
                    out += '0 %8s %4s  %12s %12s %12s   %12s %12s  %6.2F %12s %12s %-s\n' % (eid, iLayer + 1, o11, o22, t12, t1z, t2z, angle, major, minor, ovm)

                if eType in ['CQUAD4', 'QUAD4LC']:
                    quadMsg.append(out)
                elif eType in ['CTRIA3', 'TRIA3LC']:
                    triMsg.append(out)
                #else:
                #    raise NotImplementedError('eType = |%r|' %(eType)) # CQUAD8LC

            if isQuad:
                quadMsg.append(pageStamp + str(pageNum) + '\n')
                pageNum += 1
            if isTri:
                triMsg.append(pageStamp + str(pageNum) + '\n')
                pageNum += 1
            msg += quadMsg + triMsg

        f.write(''.join(msg))
        return pageNum - 1

    def __repr__(self):
        return self.get_stats()


class CompositePlateStrainObject(StrainObject):
    """
    ::
      ???
      ELEMENT  PLY  STRESSES IN FIBER AND MATRIX DIRECTIONS    INTER-LAMINAR  STRESSES  PRINCIPAL STRESSES (ZERO SHEAR)      MAX
        ID      ID    NORMAL-1     NORMAL-2     SHEAR-12     SHEAR XZ-MAT  SHEAR YZ-MAT  ANGLE    MAJOR        MINOR        SHEAR
    """
    def __init__(self, data_code, is_sort1, isubcase, dt=None):
        StrainObject.__init__(self, data_code, isubcase)

        self.eType = {}
        self.code = [self.format_code, self.sort_code, self.s_code]
        self.e11 = {}
        self.e22 = {}
        self.e12 = {}
        self.e1z = {}
        self.e2z = {}
        self.angle = {}
        self.majorP = {}
        self.minorP = {}

        if self.code == [1, 0, 14]:
            self.evmShear = {}
            assert self.isVonMises() == False
        else:
            raise RuntimeError("Invalid Code: compositePlateStrain - get the format/sort/stressCode=%s" % (self.code))

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
            ntimes = len(self.e11)
            msg.append('  type=%s ntimes=%s nelements=%s\n'
                       % (self.__class__.__name__, ntimes, nelements))
        else:
            msg.append('  type=%s nelements=%s\n' % (self.__class__.__name__,
                                                     nelements))
        msg.append('  eType, e11, e22, e12, e1z, e2z, angle, majorP, minorP\n')
        return msg

    def get_transients(self):
        k = self.e11.keys()
        k.sort()
        return k

    def add_new_transient(self, dt):
        """
        initializes the transient variables
        """
        #self.fiberDistance[dt] = {}
        self.e11[dt] = {}
        self.e22[dt] = {}
        self.e12[dt] = {}
        self.e1z[dt] = {}
        self.e2z[dt] = {}
        self.angle[dt] = {}
        self.majorP[dt] = {}
        self.minorP[dt] = {}
        self.evmShear[dt] = {}

    def add_new_eid(self, eType, dt, eid, e11, e22, e12, e1z, e2z, angle, majorP, minorP, evm):
        """all points are located at the centroid"""
        #print "Composite Plate Strain add..."
        if eid in self.e11:
            return self.add(dt, eid, e11, e22, e12, e1z, e2z, angle, majorP, minorP, evm)
        assert eid not in self.e11
        assert isinstance(eid, int)
        self.eType[eid] = eType
        self.e11[eid] = [e11]
        self.e22[eid] = [e22]
        self.e12[eid] = [e12]
        self.e1z[eid] = [e1z]
        self.e2z[eid] = [e2z]
        self.angle[eid] = [angle]
        self.majorP[eid] = [majorP]
        self.minorP[eid] = [minorP]
        self.evmShear[eid] = [evm]
        #msg = "eid=%s e11=%g e22=%g e12=%g e1z=%g e2z=%g \nangle=%g major=%g minor=%g vm=%g" % (eid, e11, e22, e12, e1z, e2z, angle, majorP, minorP, evm)
        #print msg
        #if nodeID==0: raise Exception(msg)

    def add_new_eid_sort1(self, eType, dt, eid, e11, e22, e12, e1z, e2z, angle, majorP, minorP, evm):
        """all points are located at the centroid"""
        #print "Composite Plate Strain add..."

        if dt not in self.e11:
            self.add_new_transient(dt)
        assert eid not in self.e11[dt]
        assert isinstance(eid, int)

        self.eType[eid] = eType
        self.e11[dt][eid] = [e11]
        self.e22[dt][eid] = [e22]
        self.e12[dt][eid] = [e12]
        self.e1z[dt][eid] = [e1z]
        self.e2z[dt][eid] = [e2z]
        self.angle[dt][eid] = [angle]
        self.majorP[dt][eid] = [majorP]
        self.minorP[dt][eid] = [minorP]
        self.evmShear[dt][eid] = [evm]
        msg = "eid=%s e11=%g e22=%g e12=%g e1z=%g e2z=%g \nangle=%g major=%g minor=%g vm=%g" % (eid, e11, e22, e12, e1z, e2z, angle, majorP, minorP, evm)
        #print msg
        #if nodeID==0: raise Exception(msg)

    def add(self, dt, eid, e11, e22, e12, e1z, e2z, angle, majorP, minorP, evm):
        #print "***add"
        msg = "eid=%s e11=%g e22=%g e12=%g e1z=%g e2z=%g \nangle=%g major=%g minor=%g vm=%g" % (eid, e11, e22, e12, e1z, e2z, angle, majorP, minorP, evm)
        #print msg
        #print self.o11
        self.e11[eid].append(e11)
        self.e22[eid].append(e22)
        self.e12[eid].append(e12)
        self.e1z[eid].append(e1z)
        self.e2z[eid].append(e2z)
        self.angle[eid].append(angle)
        self.majorP[eid].append(majorP)
        self.minorP[eid].append(minorP)
        self.evmShear[eid].append(evm)
        #if nodeID==0: raise Exception(msg)

    def add_sort1(self, dt, eid, e11, e22, e12, e1z, e2z, angle, majorP, minorP, evm):
        #print "***add"
        #msg = "eid=%s e11=%g e22=%g e12=%g e1z=%g e2z=%g \nangle=%g major=%g minor=%g vm=%g" % (eid, e11, e22, e12, e1z, e2z, angle, majorP, minorP, evm)
        #print msg
        #print self.o11
        self.e11[dt][eid].append(e11)
        self.e22[dt][eid].append(e22)
        self.e12[dt][eid].append(e12)
        self.e1z[dt][eid].append(e1z)
        self.e2z[dt][eid].append(e2z)
        self.angle[dt][eid].append(angle)
        self.majorP[dt][eid].append(majorP)
        self.minorP[dt][eid].append(minorP)
        self.evmShear[dt][eid].append(evm)
        #if nodeID==0: raise Exception(msg)

    def getHeaders(self):
        headers = ['e11', 'e22', 'e12', 'e1z', 'e2z']
        if self.isVonMises:
            headers.append('eVonMises')
        else:
            headers.append('maxShear')
        return headers

    def write_f06(self, header, pageStamp, f, pageNum=1, is_mag_phase=False):
        if self.nonlinear_factor is not None:
            return self._write_f06_transient(header, pageStamp, f, pageNum)

        if self.isVonMises():
            von = 'VON'
            mises = 'MISES'
        else:
            von = 'MAX'
            mises = 'SHEAR'

        words = ['   ELEMENT  PLY   STRAINS IN FIBER AND MATRIX DIRECTIONS    INTER-LAMINAR   STRAINS  PRINCIPAL  STRAINS (ZERO SHEAR)      %s\n' % (von),
                 '     ID      ID    NORMAL-1     NORMAL-2     SHEAR-12     SHEAR XZ-MAT  SHEAR YZ-MAT  ANGLE    MAJOR        MINOR        %s\n' % (mises)]

        eTypes = self.eType.values()
        if 'CQUAD4' in eTypes or 'QUAD4LC' in eTypes:
            quadMsg = header + ['                     S T R A I N S   I N   L A Y E R E D   C O M P O S I T E   E L E M E N T S   ( Q U A D 4 )\n'] + words
            isQuad = True
        else:
            quadMsg = []
            isQuad = False

        if 'CTRIA3' in eTypes or 'TRIA3LC' in eTypes:
            isTri = True
            triMsg = header + ['                     S T R A I N S   I N   L A Y E R E D   C O M P O S I T E   E L E M E N T S   ( T R I A 3 )\n'] + words
        else:
            isTri = False
            triMsg = []

        for eid, e11s in sorted(self.e11.iteritems()):
            out = ''
            eType = self.eType[eid]
            for iLayer in xrange(len(e11s)):
                e11 = self.e11[eid][iLayer]
                e22 = self.e22[eid][iLayer]
                e12 = self.e12[eid][iLayer]
                e1z = self.e1z[eid][iLayer]
                e2z = self.e2z[eid][iLayer]

                angle = self.angle[eid][iLayer]
                major = self.majorP[eid][iLayer]
                minor = self.minorP[eid][iLayer]
                evm = self.evmShear[eid][iLayer]

                (vals2, isAllZeros) = writeFloats12E([e11,
                                                      e22, e12, e1z, e2z, major, minor, evm])
                [e11, e22, e12, e1z, e2z, major, minor, evm] = vals2
                out += '0 %8s %4s  %12s %12s %12s   %12s %12s  %6.2F %12s %12s %-s\n' % (eid, iLayer + 1, e11, e22, e12, e1z, e2z, angle, major, minor, evm)

            if eType in ['CQUAD4', 'QUAD4LC']:
                quadMsg.append(out)
            elif eType in ['CTRIA3', 'TRIA3LC']:
                triMsg.append(out)
            #else:
            #    raise NotImplementedError('eType = |%r|' %(eType)) # CQUAD8LC

        if isQuad:
            quadMsg.append(pageStamp + str(pageNum) + '\n')
            pageNum += 1
        if isTri:
            triMsg.append(pageStamp + str(pageNum) + '\n')
            pageNum += 1

        msg = ''.join(quadMsg + triMsg)
        f.write(''.join(msg))
        return pageNum

    def _write_f06_transient(self, header, pageStamp, f, pageNum=1, is_mag_phase=False):
        if self.isVonMises():
            von = 'VON'
            mises = 'MISES'
        else:
            von = 'MAX'
            mises = 'SHEAR'

        words = ['   ELEMENT  PLY  STRESSES IN FIBER AND MATRIX DIRECTIONS    INTER-LAMINAR  STRESSES  PRINCIPAL STRESSES (ZERO SHEAR)      %s\n' % (von),
                 '     ID      ID    NORMAL-1     NORMAL-2     SHEAR-12     SHEAR XZ-MAT  SHEAR YZ-MAT  ANGLE    MAJOR        MINOR        %s\n' % (mises)]

        eTypes = self.eType.values()
        if 'CQUAD4' in eTypes or 'QUAD4LC' in eTypes:
            quadWords = ['                     S T R A I N S   I N   L A Y E R E D   C O M P O S I T E   E L E M E N T S   ( Q U A D 4 )\n'] + words
            isQuad = True
        else:
            quadWords = []
            isQuad = False

        if 'CTRIA3' in eTypes or 'TRIA3LC' in eTypes:
            isTri = True
            triWords = ['                     S T R A I N S   I N   L A Y E R E D   C O M P O S I T E   E L E M E N T S   ( T R I A 3 )\n'] + words
        else:
            isTri = False
            triWords = []

        msg = []
        for dt, e11s in sorted(self.e11.iteritems()):
            quadMsg = []
            triMsg = []
            header[1] = ' %s = %10.4E\n' % (self.data_code['name'], dt)
            if isQuad:
                quadMsg = header + quadWords
            if isTri:
                triMsg = header + triWords

            for eid, e11s in sorted(e11s.iteritems()):
                out = ''
                eType = self.eType[eid]
                for iLayer in xrange(len(e11s)):
                    e11 = self.e11[dt][eid][iLayer]
                    e22 = self.e22[dt][eid][iLayer]
                    e12 = self.e12[dt][eid][iLayer]
                    e1z = self.e1z[dt][eid][iLayer]
                    e2z = self.e2z[dt][eid][iLayer]

                    angle = self.angle[dt][eid][iLayer]
                    major = self.majorP[dt][eid][iLayer]
                    minor = self.minorP[dt][eid][iLayer]
                    evm = self.evmShear[dt][eid][iLayer]
                    (vals2, isAllZeros) = writeFloats12E([e11, e22,
                                                          e12, e1z, e2z, major, minor, evm])
                    [e11, e22, e12, e1z, e2z, major, minor, evm] = vals2
                    out += '0 %8s %4s  %12s %12s %12s   %12s %12s  %6.2F %12s %12s %-s\n' % (eid, iLayer + 1, e11, e22, e12, e1z, e2z, angle, major, minor, evm)

                if eType in ['CQUAD4', 'QUAD4LC']:
                    quadMsg.append(out)
                elif eType in ['CTRIA3', 'TRIA3LC']:
                    triMsg.append(out)
                else:
                    raise NotImplementedError('eType = |%r|' % (eType))

            if isQuad:
                quadMsg.append(pageStamp + str(pageNum) + '\n')
                pageNum += 1
            if isTri:
                triMsg.append(pageStamp + str(pageNum) + '\n')
                pageNum += 1
            msg += quadMsg + triMsg

        f.write(''.join(msg))
        return pageNum - 1

    def __repr__(self):
        return self.get_stats()