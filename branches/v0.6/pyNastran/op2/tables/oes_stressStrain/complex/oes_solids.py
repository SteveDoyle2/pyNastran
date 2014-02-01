from __future__ import (nested_scopes, generators, division, absolute_import,
                        print_function, unicode_literals)
from numpy import array, sqrt
from numpy.linalg import eigh

from ..real.oes_objects import StressObject, StrainObject
from pyNastran.f06.f06_formatting import writeFloats13E
from pyNastran.f06.f06_formatting import writeImagFloats13E


class ComplexSolidStressObject(StressObject):
    def __init__(self, data_code, is_sort1, isubcase, dt):
        StressObject.__init__(self, data_code, isubcase)

        self.eType = {}
        self.code = [self.format_code, self.sort_code, self.s_code]

        self.cid = {}  # gridGauss
        self.oxx = {}
        self.oyy = {}
        self.ozz = {}
        self.txy = {}
        self.tyz = {}
        self.txz = {}

        #self.dt = dt
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
            msg.append('  type=%s nelements=%s\n' % (self.__class__.__name__, nelements))
        msg.append('  eType, cid, oxx, oyy, ozz, txy, tyz, txz\n  ')
        msg.append(', '.join(set(self.eType.values())))
        msg.append('\n')
        return msg

    def add_f06_data(self, data, transient):
        if transient is None:
            if not hasattr(self, 'data'):
                self.data = []
            self.data += data
        else:
            dt = transient[1]
            if not hasattr(self, 'data'):
                self.data = {}
            #print(self.data)
            if dt not in self.data:
                self.data[dt] = []
            for line in data:
                self.data[dt] += data

    def processF06Data(self):
        raise NotImplementedError()

    def delete_transient(self, dt):
        del self.oxx[dt]
        del self.oyy[dt]
        del self.ozz[dt]
        del self.txy[dt]
        del self.tyz[dt]
        del self.txz[dt]

    def get_transients(self):
        k = self.oxx.keys()
        k.sort()
        return k

    def add_new_transient(self, dt):
        """
        initializes the transient variables
        """
        self.oxx[dt] = {}
        self.oyy[dt] = {}
        self.ozz[dt] = {}
        self.txy[dt] = {}
        self.tyz[dt] = {}
        self.txz[dt] = {}

    def add_new_eid_sort1(self, eType, dt, eid, cid, ctype, nodef):
        #print "Solid Stress add transient..."
        assert cid >= 0
        assert eid >= 0

        #print "dt=%s eid=%s eType=%s" %(dt,eid,eType)
        if dt not in self.oxx:
            self.add_new_transient(dt)
        assert eid not in self.oxx[dt], self.oxx[dt]

        self.eType[eid] = eType
        self.cid[eid] = cid
        self.oxx[dt][eid] = {}
        self.oyy[dt][eid] = {}
        self.ozz[dt][eid] = {}
        self.txy[dt][eid] = {}
        self.tyz[dt][eid] = {}
        self.txz[dt][eid] = {}

        #msg = "*eid=%s nodeID=%s vm=%g" % (eid, nodeID, ovm)
        #print msg
        #if nodeID == 0:
            #raise ValueError(msg)

    def add_sort1(self, dt, eid, nodeID, oxx, oyy, ozz, txy, tyz, tzx):
        #print "***add"
        #msg = "eid=%s nodeID=%s vm=%g" % (eid, nodeID, ovm)
        #print msg
        #print self.oxx
        #print self.fiberDistance

        #self.eType[eid] = eType
        assert eid != 'C'
        #print("eid=%s nid=%s oyy=%s" %(eid,nodeID, oyy))
        self.oxx[dt][eid][nodeID] = oxx
        self.oyy[dt][eid][nodeID] = oyy
        self.ozz[dt][eid][nodeID] = ozz
        #print self.oxx
        self.txy[dt][eid][nodeID] = txy
        self.tyz[dt][eid][nodeID] = tyz
        self.txz[dt][eid][nodeID] = tzx

    def getHeaders(self):
        headers = ['oxx', 'oyy', 'ozz', 'txy', 'tyz', 'txz']
        return headers

    def directionalVectors(self, oxx, oyy, ozz, txy, tyz, txz):
        A = [[oxx, txy, txz],
             [txy, oyy, tyz],
             [txz, tyz, ozz]]
        (Lambda, v) = eigh(A)  # a hermitian matrix is a symmetric-real matrix
        return v

    def getF06_Header(self):
        raise NotImplementedError()

    def write_f06(self, header, pageStamp, pageNum, f, is_mag_phase=False):
        """
        Not perfect, but it's not bad...
        """
        if is_mag_phase:
            asdf
        else:
            stamp = [
                '                 C O M P L E X     S T R A I N S   I N   T E T R A H E D R O N   E L E M E N T S   ( C T E T R A )',
                '                                                          (REAL/IMAGINARY)',
                '0                   CORNER      --------------------------CENTER AND CORNER POINT  STRAINS---------------------------',
                '     ELEMENT-ID    GRID-ID      NORMAL-X       NORMAL-Y       NORMAL-Z         SHEAR-XY       SHEAR-YZ       SHEAR-ZX',
            ]
        for dt, oxx_dt in sorted(self.oxx.iteritems()):
            f.write('\n'.join(stamp))
            f.write('FREQ = %r\n' % dt)
            for elem, oxx_elem in sorted(oxx_dt.iteritems()):
                node_ids = oxx_elem.keys()
                node_ids.remove('CENTER')
                f.write('%8s\n' % elem)
                for inode in ['CENTER'] + sorted(node_ids):
                    # cid
                    cid = 1
                    oxx = self.oxx[dt][elem][inode]
                    oyy = self.oyy[dt][elem][inode]
                    ozz = self.ozz[dt][elem][inode]
                    txy = self.txy[dt][elem][inode]
                    tyz = self.tyz[dt][elem][inode]
                    txz = self.txz[dt][elem][inode]
                    ([oxxr, oyyr, ozzr, txyr, tyzr, txzr,
                      oxxi, oyyi, ozzi, txyi, tyzi, txzi,], isAllZeros) = writeImagFloats13E([oxx, oyy, ozz,
                                                                                              txy, tyz, txz], is_mag_phase)

                    if inode == 'CENTER':
                        f.write('0   %8s %8s %8s %8s %8s %8s %8s %8s\n' % (inode, cid, oxxr, oyyr, ozzr, txyr, tyzr, txzr))
                        f.write('    %8s %8s %8s %8s %8s %8s %8s %8s\n' % ('', '',     oxxi, oyyi, ozzi, txyi, tyzi, txzi))
                    else:
                        f.write('    %8s %8s %8s %8s %8s %8s %8s %8s\n' % (inode, cid, oxxr, oyyr, ozzr, txyr, tyzr, txzr))
                        f.write('    %8s %8s %8s %8s %8s %8s %8s %8s\n' % ('', '',     oxxi, oyyi, ozzi, txyi, tyzi, txzi))
            f.write(pageStamp % pageNum)
            pageNum += 1
        return pageNum - 1

class ComplexSolidStrainObject(StrainObject):
    def __init__(self, data_code, is_sort1, isubcase, dt):
        StrainObject.__init__(self, data_code, isubcase)

        self.eType = {}
        self.code = [self.format_code, self.sort_code, self.s_code]

        self.cid = {}  # gridGauss
        self.exx = {}
        self.eyy = {}
        self.ezz = {}
        self.exy = {}
        self.eyz = {}
        self.exz = {}

        #self.dt = dt
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
            ntimes = len(self.exx)
            msg.append('  type=%s ntimes=%s nelements=%s\n'
                       % (self.__class__.__name__, ntimes, nelements))
        else:
            msg.append('  type=%s nelements=%s\n' % (self.__class__.__name__, nelements))
        msg.append('  eType, cid, exx, eyy, ezz, exy, eyz, exz\n  ')
        msg.append(', '.join(set(self.eType.values())))
        msg.append('\n')
        return msg

    def add_f06_data(self, data, transient):
        if transient is None:
            if not hasattr(self, 'data'):
                self.data = []
            self.data += data
        else:
            dt = transient[1]
            if not hasattr(self, 'data'):
                self.data = {}
            #print(self.data)
            if dt not in self.data:
                self.data[dt] = []
            for line in data:
                self.data[dt] += data

    def processF06Data(self):
        raise NotImplementedError()

    def delete_transient(self, dt):
        del self.exx[dt]
        del self.eyy[dt]
        del self.ezz[dt]
        del self.exy[dt]
        del self.eyz[dt]
        del self.exz[dt]

    def get_transients(self):
        k = self.exx.keys()
        k.sort()
        return k

    def add_new_transient(self, dt):
        """
        initializes the transient variables
        """
        self.exx[dt] = {}
        self.eyy[dt] = {}
        self.ezz[dt] = {}
        self.exy[dt] = {}
        self.eyz[dt] = {}
        self.exz[dt] = {}

    def add_new_eid_sort1(self, eType, dt, eid, cid, ctype, nodef):
        #print "Solid Stress add transient..."
        assert cid >= 0
        assert eid >= 0

        #print "dt=%s eid=%s eType=%s" %(dt,eid,eType)
        if dt not in self.exx:
            self.add_new_transient(dt)
        assert eid not in self.exx[dt], self.exx[dt]

        self.eType[eid] = eType
        self.cid[eid] = cid
        self.exx[dt][eid] = {}
        self.eyy[dt][eid] = {}
        self.ezz[dt][eid] = {}
        self.exy[dt][eid] = {}
        self.eyz[dt][eid] = {}
        self.exz[dt][eid] = {}

        #msg = "*eid=%s nodeID=%s vm=%g" % (eid, nodeID, ovm)
        #print msg
        #if nodeID == 0:
            #raise ValueError(msg)

    def add_sort1(self, dt, eid, nodeID, ex, ey, ez, etxy, etyz, etzx):
        #print "***add"
        #msg = "eid=%s nodeID=%s vm=%g" % (eid, nodeID, evm)
        #print msg
        #print self.oxx
        #print self.fiberDistance

        #self.eType[eid] = eType
        assert eid != 'C'
        #print "eid=%s nid=%s exx=%s" %(eid, nodeID, exx)
        self.exx[dt][eid][nodeID] = ex
        self.eyy[dt][eid][nodeID] = ey
        self.ezz[dt][eid][nodeID] = ez
        #print self.oxx
        self.exy[dt][eid][nodeID] = etxy
        self.eyz[dt][eid][nodeID] = etyz
        self.exz[dt][eid][nodeID] = etzx

    def getHeaders(self):
        headers = ['exx', 'eyy', 'ezz', 'exy', 'eyz', 'exz']
        return headers

    def directionalVectors(self, exx, eyy, ezz, exy, eyz, exz):
        A = [[exx, exy, exz],
             [exy, eyy, eyz],
             [exz, eyz, ezz]]
        (Lambda, v) = eigh(A)  # a hermitian matrix is a symmetric-real matrix
        return v

    def getF06_Header(self):
        raise NotImplementedError()
