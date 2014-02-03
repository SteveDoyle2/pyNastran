#pylint: disable=C0301,C0111,R0921
from __future__ import (nested_scopes, generators, division, absolute_import,
                        print_function, unicode_literals)
from numpy.linalg import eig

from ..real.oes_objects import StressObject, StrainObject
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

        if is_sort1:
            if dt is not None:
                self.add = self.add_sort1
                self.add_new_eid = self.add_new_eid_sort1
        else:
            assert dt is not None
            raise NotImplementedError('SORT2')
            #self.add = self.addSort2
            #self.add_new_eid = self.add_new_eid_sort2

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

        assert eid != 'C'
        #print("eid=%s nid=%s oyy=%s" %(eid,nodeID, oyy))
        self.oxx[dt][eid][nodeID] = oxx
        self.oyy[dt][eid][nodeID] = oyy
        self.ozz[dt][eid][nodeID] = ozz

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
        (Lambda, v) = eig(A)  # we can't use a hermitian matrix
        return v

    def getF06_Header(self, is_mag_phase):
        if is_mag_phase:
            tetra_msg = [
                '                 C O M P L E X   S T R E S S E S   I N   T E T R A H E D R O N   E L E M E N T S   ( C T E T R A )',
                '                                                          (MAGNITUDE/PHASE)',
                '0                   CORNER      --------------------------CENTER AND CORNER POINT STRESSES---------------------------',
                '     ELEMENT-ID    GRID-ID      NORMAL-X       NORMAL-Y       NORMAL-Z         SHEAR-XY       SHEAR-YZ       SHEAR-ZX',
                '',
            ]
        else:
            tetra_msg = [
                '                 C O M P L E X   S T R E S S E S   I N   T E T R A H E D R O N   E L E M E N T S   ( C T E T R A )',
                '                                                          (REAL/IMAGINARY)',
                '0                   CORNER      --------------------------CENTER AND CORNER POINT STRESSES---------------------------',
                '     ELEMENT-ID    GRID-ID      NORMAL-X       NORMAL-Y       NORMAL-Z         SHEAR-XY       SHEAR-YZ       SHEAR-ZX',
                '',
            ]
        penta_msg = tetra_msg
        hexa_msg = tetra_msg

        tetra_eids = []
        hexa_eids = []
        penta_eids = []

        tetra10_eids = []
        hexa20_eids = []
        penta15_eids = []
        for eid, eType in sorted(self.eType.iteritems()):
            if eType in ['CTETRA4']:
                tetra_eids.append(eid)
            elif eType in ['CTETRA10']:
                tetra10_eids.append(eid)

            elif eType in ['CPENTA']:
                penta_eids.append(eid)
            elif eType in ['CPENTA15']:
                penta_eids.append(eid)

            elif eType in ['CHEXA8']:
                hexa_eids.append(eid)
            elif eType in ['CHEXA20']:
                hexa20_eids.append(eid)

            else:
                raise NotImplementedError('eType=%r' % eType)
        return (tetra_msg, hexa_msg, penta_msg,
                tetra_eids, hexa_eids, penta_eids,
                tetra10_eids, hexa20_eids, penta15_eids)

    def write_f06(self, header, page_stamp, pageNum=1, f=None, is_mag_phase=False):
        (tetra_msg, hexa_msg, penta_msg,
         tetra_eids, hexa_eids, penta_eids,
         tetra10_eids, hexa20_eids, penta15_eids) = self.getF06_Header(is_mag_phase)
        dts = self.oxx.keys()
        for dt in sorted(dts):
            if tetra_eids:
                self.write_element_transient('CTETRA', 4, tetra_eids, dt, header, tetra_msg, f, is_mag_phase)
                f.write(page_stamp % pageNum)
                pageNum += 1
            if tetra10_eids:
                self.write_element_transient('CTETRA', 10, tetra10_eids, dt, header, tetra_msg, f, is_mag_phase)
                f.write(page_stamp % pageNum)
                pageNum += 1

            if hexa_eids:
                self.write_element_transient('CHEXA',  8,  hexa_eids, dt, header, hexa_msg, f, is_mag_phase)
                f.write(page_stamp % pageNum)
                pageNum += 1
            if hexa20_eids:
                self.write_element_transient('CHEXA',  20,  hexa20_eids, dt, header, hexa_msg, f,is_mag_phase)
                f.write(page_stamp % pageNum)
                pageNum += 1

            if penta_eids:
                self.write_element_transient('CPENTA', 6, penta_eids, dt, header, penta_msg, f, is_mag_phase)
                f.write(page_stamp % pageNum)
                pageNum += 1
            if penta15_eids:
                self.write_element_transient('CPENTA', 15, penta15_eids, dt, header, penta_msg, f,is_mag_phase)
                f.write(page_stamp % pageNum)
                pageNum += 1

        return pageNum - 1

    def write_element_transient(self, element_name, nnodes, eids, dt, header, msg, f, is_mag_phase):
        dtLine = ' %14s = %12.5E\n' % (self.data_code['name'], dt)
        header[1] = dtLine
        msg = header + msg

        f.write('\n'.join(msg))
        for eid in eids:
            node_ids = self.oxx[dt][eid].keys()
            node_ids.remove('CENTER')
            cid = 10
            f.write('0 %12i %11sGRID CS %2i GP\n' % (eid, 0, nnodes))
            for inode in ['CENTER'] + sorted(node_ids):
                # cid
                oxx = self.oxx[dt][eid][inode]
                oyy = self.oyy[dt][eid][inode]
                ozz = self.ozz[dt][eid][inode]
                txy = self.txy[dt][eid][inode]
                tyz = self.tyz[dt][eid][inode]
                txz = self.txz[dt][eid][inode]
                ([oxxr, oyyr, ozzr, txyr, tyzr, txzr,
                  oxxi, oyyi, ozzi, txyi, tyzi, txzi,], isAllZeros) = writeImagFloats13E([oxx, oyy, ozz,
                                                                                          txy, tyz, txz], is_mag_phase)

                f.write('0   %22s    %-13s  %-13s  %-13s    %-13s  %-13s  %s\n' % (inode, oxxr, oyyr, ozzr, txyr, tyzr, txzr))
                f.write('    %22s    %-13s  %-13s  %-13s    %-13s  %-13s  %s\n' % ('',    oxxi, oyyi, ozzi, txyi, tyzi, txzi))


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
            raise NotImplementedError('SORT2')
            #self.add = self.add_sort2
            #self.add_new_eid = self.add_new_eid_sort2

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
        assert cid >= 0
        assert eid >= 0

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
        #msg = "eid=%s nodeID=%s vm=%g" % (eid, nodeID, evm)
        assert eid != 'C'
        #print "eid=%s nid=%s exx=%s" %(eid, nodeID, exx)
        self.exx[dt][eid][nodeID] = ex
        self.eyy[dt][eid][nodeID] = ey
        self.ezz[dt][eid][nodeID] = ez

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
        (Lambda, v) = eig(A)  # we can't use a hermitian matrix
        return v

    def getF06_Header(self, is_mag_phase):
        if is_mag_phase:
            tetra_msg = [
                '                 C O M P L E X     S T R A I N S   I N   T E T R A H E D R O N   E L E M E N T S   ( C T E T R A )',
                '                                                          (MAGNITUDE/PHASE)',
                '0                   CORNER      --------------------------CENTER AND CORNER POINT  STRAINS---------------------------',
                '     ELEMENT-ID    GRID-ID      NORMAL-X       NORMAL-Y       NORMAL-Z         SHEAR-XY       SHEAR-YZ       SHEAR-ZX',
                '',
            ]
        else:
            tetra_msg = [
                '                 C O M P L E X     S T R A I N S   I N   T E T R A H E D R O N   E L E M E N T S   ( C T E T R A )',
                '                                                          (REAL/IMAGINARY)',
                '0                   CORNER      --------------------------CENTER AND CORNER POINT  STRAINS---------------------------',
                '     ELEMENT-ID    GRID-ID      NORMAL-X       NORMAL-Y       NORMAL-Z         SHEAR-XY       SHEAR-YZ       SHEAR-ZX',
                '',
            ]
        penta_msg = tetra_msg
        hexa_msg = tetra_msg

        tetra_eids = []
        hexa_eids = []
        penta_eids = []

        tetra10_eids = []
        hexa20_eids = []
        penta15_eids = []
        for eid, eType in sorted(self.eType.iteritems()):
            if eType in ['CTETRA4']:
                tetra_eids.append(eid)
            elif eType in ['CTETRA10']:
                tetra10_eids.append(eid)

            elif eType in ['CPENTA']:
                penta_eids.append(eid)
            elif eType in ['CPENTA15']:
                penta_eids.append(eid)

            elif eType in ['CHEXA8']:
                hexa_eids.append(eid)
            elif eType in ['CHEXA20']:
                hexa20_eids.append(eid)

            else:
                raise NotImplementedError('eType=%r' % eType)
        return (tetra_msg, hexa_msg, penta_msg,
                tetra_eids, hexa_eids, penta_eids,
                tetra10_eids, hexa20_eids, penta15_eids)

    def write_f06(self, header, page_stamp, pageNum=1, f=None, is_mag_phase=False):
        (tetra_msg, hexa_msg, penta_msg,
         tetra_eids, hexa_eids, penta_eids,
         tetra10_eids, hexa20_eids, penta15_eids) = self.getF06_Header(is_mag_phase)
        dts = self.exx.keys()
        for dt in sorted(dts):
            if tetra_eids:
                self.write_element_transient('CTETRA', 4, tetra_eids, dt, header, tetra_msg, f, is_mag_phase)
                f.write(page_stamp % pageNum)
                pageNum += 1
            if tetra10_eids:
                self.write_element_transient('CTETRA', 10, tetra10_eids, dt, header, tetra_msg, f, is_mag_phase)
                f.write(page_stamp % pageNum)
                pageNum += 1

            if hexa_eids:
                self.write_element_transient('CHEXA',  8,  hexa_eids, dt, header, hexa_msg, f, is_mag_phase)
                f.write(page_stamp % pageNum)
                pageNum += 1
            if hexa20_eids:
                self.write_element_transient('CHEXA',  20,  hexa20_eids, dt, header, hexa_msg, f,is_mag_phase)
                f.write(page_stamp % pageNum)
                pageNum += 1

            if penta_eids:
                self.write_element_transient('CPENTA', 6, penta_eids, dt, header, penta_msg, f, is_mag_phase)
                f.write(page_stamp % pageNum)
                pageNum += 1
            if penta15_eids:
                self.write_element_transient('CPENTA', 15, penta15_eids, dt, header, penta_msg, f,is_mag_phase)
                f.write(page_stamp % pageNum)
                pageNum += 1

        return pageNum - 1

    def write_element_transient(self, element_name, nnodes, eids, dt, header, msg, f, is_mag_phase):
        dtLine = ' %14s = %12.5E\n' % (self.data_code['name'], dt)
        header[1] = dtLine
        msg = header + msg

        f.write('\n'.join(msg))
        for eid in eids:
            node_ids = self.exx[dt][eid].keys()
            node_ids.remove('CENTER')
            cid = 10
            f.write('0 %12i %11sGRID CS %2i GP\n' % (eid, 0, nnodes))
            for inode in ['CENTER'] + sorted(node_ids):
                # cid
                oxx = self.exx[dt][eid][inode]
                oyy = self.eyy[dt][eid][inode]
                ozz = self.ezz[dt][eid][inode]
                txy = self.exy[dt][eid][inode]
                tyz = self.eyz[dt][eid][inode]
                txz = self.exz[dt][eid][inode]
                ([oxxr, oyyr, ozzr, txyr, tyzr, txzr,
                  oxxi, oyyi, ozzi, txyi, tyzi, txzi,], isAllZeros) = writeImagFloats13E([oxx, oyy, ozz,
                                                                                          txy, tyz, txz], is_mag_phase)

                f.write('0   %22s    %-13s  %-13s  %-13s    %-13s  %-13s  %s\n' % (inode, oxxr, oyyr, ozzr, txyr, tyzr, txzr))
                f.write('    %22s    %-13s  %-13s  %-13s    %-13s  %-13s  %s\n' % ('',    oxxi, oyyi, ozzi, txyi, tyzi, txzi))
