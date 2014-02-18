# pylint: disable=C0301,W0613,C0103,R0913,R0914,R0904,C0111,R0201,R0902
from __future__ import (nested_scopes, generators, division, absolute_import,
                        print_function, unicode_literals)
from collections import defaultdict
from itertools import izip, count
from numpy import array, sqrt, zeros, where, searchsorted, ravel, argwhere, unique, arange
from numpy.linalg import eigh

from .oes_objects import StressObject, StrainObject, OES_Object
from pyNastran.f06.f06_formatting import writeFloats13E


class RealSolidVector(OES_Object):
    def __init__(self, data_code, is_sort1, isubcase, dt):
        OES_Object.__init__(self, data_code, isubcase, apply_data_code=False)
        self.eType = {}
        #self.code = [self.format_code, self.sort_code, self.s_code]

        #self.ntimes = 0  # or frequency/mode
        #self.ntotal = 0
        self.nelements = 0  # result specific

        if is_sort1:
            #sort1
            self.add_node = self.add_node_sort1
            self.add_eid = self.add_eid_sort1
        else:
            raise NotImplementedError('SORT2')

    def _get_msgs(self):
        raise NotImplementedError()

    def get_headers(self):
        raise NotImplementedError()

    def _reset_indices(self):
        self.itotal = 0
        self.ielement = 0

    def build(self):
        #print('ntimes=%s nelements=%s ntotal=%s' % (self.ntimes, self.nelements, self.ntotal))
        if self.is_built:
            return

        #self.names = []
        self.nelements //= self.ntimes
        self.itime = 0
        self.ielement = 0
        self.itotal = 0
        self.is_built = True

        #print("ntimes=%s nelements=%s ntotal=%s" % (self.ntimes, self.nelements, self.ntotal))
        self.times = zeros(self.ntimes, 'float32')
        #self.element_types2 = array(self.nelements, dtype='|S8')
        self.element_types3 = zeros((self.nelements, 2), dtype='int32')

        # TODO: could be more efficient by using nelements for cid
        self.element_node = zeros((self.ntotal, 2), 'int32')
        self.element_cid = zeros((self.nelements, 2), 'int32')

        #if self.element_name == 'CTETRA':
            #nnodes = 4
        #elif self.element_name == 'CPENTA':
            #nnodes = 6
        #elif self.element_name == 'CHEXA':
            #nnodes = 8
        #self.element_node = zeros((self.ntotal, nnodes, 2), 'int32')

        #[oxx, oyy, ozz, txy, tyz, txz, o1, o2, o3, ovmShear]
        self.data = zeros((self.ntimes, self.ntotal, 10), 'float32')
        #self.data = zeros((self.ntimes, self.nelements, nnodes+1, 10), 'float32')

    def _add_eid_sort1(self, element_num, element_type, dt, eid, cid, ctype, nodef):
        self.times[self.itime] = dt
        #print(self.element_types2, element_type, self.element_types2.dtype)
        #self.element_types2[self.ielement] = string_(element_type)   # TODO: save this...
        #self.element_types2[self.ielement] = element_type

        #try:
        if self.ielement < self.nelements:
            self.element_cid[self.ielement] = [eid, cid]
            self.element_types3[self.ielement, :] = [element_num, nodef]
        #except IndexError:
            #pass
            #print('element_types3', self.element_types3)

        self.element_node[self.itotal, :] = [eid, 0]  # 0 is center
        # no data
        #self.element_node[self.ielement, 0, :] = [eid, 0]  # 0 is center
        self.ielement += 1
        #self.itotal += 1
        #self.data

    def add_eid_sort1(self, eType, cid, dt, eid, node_id, oxx, oyy, ozz, txy, tyz, txz, o1, o2, o3, aCos, bCos, cCos, pressure, ovm):
        #assert cid >= 0
        assert eid >= 0

        #print "dt=%s eid=%s eType=%s" %(dt,eid,eType)
        self.times[self.itime] = dt

        self.eType[self.ielement] = eType
        self.element_node[self.itotal, :] = [eid, 0]  # 0 is center
        self.data[self.itime, self.itotal, :] = [oxx, oyy, ozz, txy, tyz, txz, o1, o2, o3, ovm]
        #self.data[self.itime, self.ielement, 0, :] = [oxx, oyy, ozz, txy, tyz, txz, o1, o2, o3, ovm]

        #print('element_cid[%i, :] = [%s, %s]' % (self.ielement, eid, cid))
        if self.ielement == self.nelements:
            self.ielement = 0
        self.element_cid[self.ielement, :] = [eid, cid]
        self.itotal += 1
        self.ielement += 1

    def add_node_sort1(self, dt, eid, inode, node_id, oxx, oyy, ozz, txy, tyz, txz, o1, o2, o3, aCos, bCos, cCos, pressure, ovm):
        # skipping aCos, bCos, cCos, pressure
        self.data[self.itime, self.itotal, :] = [oxx, oyy, ozz, txy, tyz, txz, o1, o2, o3, ovm]
        #print('data[%s, %s, :] = %s' % (self.itime, self.itotal, str(self.data[self.itime, self.itotal, :])))

        #self.data[self.itime, self.ielement-1, self.inode, :] = [oxx, oyy, ozz, txy, tyz, txz, o1, o2, o3, ovm]

        #print('eid=%i node_id=%i exx=%s' % (eid, node_id, str(ex)))
        self.element_node[self.itotal, :] = [eid, node_id]
        #self.element_node[self.ielement-1, self.inode-1, :] = [eid, node_id]
        self.itotal += 1

    def get_stats(self):
        if not self.is_built:
            return ['<%s>\n' % self.__class__.__name__,
                    '  ntimes: %i\n' % self.ntimes,
                    '  ntotal: %i\n' % self.ntotal,
                    ]

        nelements = self.nelements
        ntimes = self.ntimes
        #ntotal = self.ntotal
        nnodes_per_element = self.element_node.shape[0] // nelements
        nnodes = self.element_node.shape[0]

        msg = []

        if self.nonlinear_factor is not None:  # transient
            msg.append('  type=%s ntimes=%i nelements=%i nnodes=%i\n  nnodes_per_element=%i (including centroid)\n'
                       % (self.__class__.__name__, ntimes, nelements, nnodes, nnodes_per_element))
            ntimes_word = 'ntimes'
        else:
            msg.append('  type=%s nelements=%i nnodes=%i\n  nodes_per_element=%i (including centroid)\n'
                       % (self.__class__.__name__, nelements, nnodes, nnodes_per_element))
            ntimes_word = 1
        msg.append('  eType, cid\n')
        #msg.append('  data.shape=%s' % str(self.data.shape))
        headers = self.get_headers()
        n = len(headers)
        msg.append('  data: [%s, nnodes, %i] where %i=[%s]\n' % (ntimes_word, n, n, str(', '.join(headers))))
        msg.append('  element types: %s\n  ' % ', '.join(self.element_names))
        msg += self.get_data_code()
        return msg

    def get_f06_header(self, is_mag_phase=True):
        tetra_msg, penta_msg, hexa_msg = self._get_msgs()
        if 'CTETRA' in self.element_name:
            msg = tetra_msg
            nnodes = 4
        elif 'CHEXA' in self.element_name:
            msg = hexa_msg
            nnodes = 8
        elif 'CPENTA' in self.element_name:
            msg = penta_msg
            nnodes = 6
        return self.element_name, nnodes, msg

    def get_element_index(self, eids):
        # elements are always sorted; nodes are not
        itot = searchsorted(eids, self.element_node[:, 0])  #[0]
        return itot

    def eid_to_element_node_index(self, eids):
        #ind = ravel([searchsortd(self.element_node[:, 0] == eid) for eid in eids])
        ind = searchsorted(eids, self.element_node[:, 0])
        #ind = ind.reshape(ind.size)
        #ind.sort()
        return ind

    def write_f06(self, header, page_stamp, pageNum=1, f=None, is_mag_phase=False):
        (elem_name, nnodes, msg_temp) = self.get_f06_header(is_mag_phase)

        # write the f06
        (ntimes, ntotal, six) = self.data.shape

        eids2 = self.element_node[:, 0]
        nodes = self.element_node[:, 1]
        for itime in xrange(ntimes):
            dt = self.times[itime]  # TODO: rename this...
            if self.nonlinear_factor is not None:
                dtLine = ' %14s = %12.5E\n' % (self.data_code['name'], dt)
                header[1] = dtLine
            f.write(''.join(header + msg_temp))

            # TODO: can I get this without a reshape?
            #print("self.data.shape=%s itime=%s ieids=%s" % (str(self.data.shape), itime, str(ieids)))
            oxx = self.data[itime, :, 0]
            oyy = self.data[itime, :, 1]
            ozz = self.data[itime, :, 2]
            txy = self.data[itime, :, 3]
            tyz = self.data[itime, :, 4]
            txz = self.data[itime, :, 5]
            o1 = self.data[itime, :, 6]
            o2 = self.data[itime, :, 7]
            o3 = self.data[itime, :, 8]
            ovm = self.data[itime, :, 9]
            p = (o1 + o2 + o3) / -3.

            # loop over all the elements and nodes
            cnnodes = nnodes + 1
            for i, deid, node_id, doxx, doyy, dozz, dtxy, dtyz, dtxz, do1, do2, do3, dp, dovm in izip(
                count(), eids2, nodes, oxx, oyy, ozz, txy, tyz, txz, o1, o2, o3, p, ovm):

                # TODO: cid not supported
                A = [[doxx, dtxy, dtxz],
                     [dtxy, doyy, dtyz],
                     [dtxz, dtyz, dozz]]
                (Lambda, v) = eigh(A)  # a hermitian matrix is a symmetric-real matrix

                ([oxxi, oyyi, ozzi, txyi, tyzi, txzi, o1i, o2i, o3i, pi, ovmi],
                 isAllZeros) = writeFloats13E([doxx, doyy, dozz, dtxy, dtyz, dtxz,
                                               do1, do2, do3, dp, dovm])

                if i % cnnodes == 0:
                    f.write('0  %8s           0GRID CS  %i GP\n' % (deid, nnodes))
                    #msg.append('0              %8s  X  %13s  XY  %13s   A  %13s  LX%5.2f%5.2f%5.2f  %-13s   %s\n'
                    f.write('0              %8s  X  %13s  XY  %13s   A  %13s  LX%5.2f%5.2f%5.2f  %-13s   %s\n'
                            '               %8s  Y  %13s  YZ  %13s   B  %13s  LY%5.2f%5.2f%5.2f\n'
                            '               %8s  Z  %13s  ZX  %13s   C  %13s  LZ%5.2f%5.2f%5.2f\n'
                            % ('CENTER', oxxi, txyi, o1i, v[0, 1], v[0, 2], v[0, 0], pi, ovmi,
                                     '', oyyi, tyzi, o2i, v[1, 1], v[1, 2], v[1, 0],
                                     '', ozzi, txzi, o3i, v[2, 1], v[2, 2], v[2, 0]))
                else:
                    f.write('0              %8s  X  %13s  XY  %13s   A  %13s  LX%5.2f%5.2f%5.2f  %-13s   %s\n'
                            '               %8s  Y  %13s  YZ  %13s   B  %13s  LY%5.2f%5.2f%5.2f\n'
                            '               %8s  Z  %13s  ZX  %13s   C  %13s  LZ%5.2f%5.2f%5.2f\n'
                            % (node_id, oxxi, txyi, o1i, v[0, 1], v[0, 2], v[0, 0], pi, ovmi,
                                    '', oyyi, tyzi, o2i, v[1, 1], v[1, 2], v[1, 0],
                                    '', ozzi, txzi, o3i, v[2, 1], v[2, 2], v[2, 0]))
                i += 1
            f.write(page_stamp % pageNum)
            pageNum += 1
        return pageNum - 1


class RealSolidStressVector(RealSolidVector, StressObject):
    def __init__(self, data_code, is_sort1, isubcase, dt):
        RealSolidVector.__init__(self, data_code, is_sort1, isubcase, dt)
        StressObject.__init__(self, data_code, isubcase)

    def get_headers(self):
        if self.isVonMises():
            von_mises = 'von_mises'
        else:
            von_mises = 'max_shear'
        headers = ['oxx', 'oyy', 'ozz', 'txy', 'tyz', 'txz', 'o1', 'o2', 'o3', von_mises]
        return headers

    def _get_msgs(self):
        if self.isVonMises():
            von_mises = 'VON MISES'
        else:
            von_mises = 'MAX SHEAR'

        base_msg = [
            '0                CORNER        ------CENTER AND CORNER POINT STRESSES---------       DIR.  COSINES       MEAN                   \n',
            '  ELEMENT-ID    GRID-ID        NORMAL              SHEAR             PRINCIPAL       -A-  -B-  -C-     PRESSURE       %s \n' % von_mises,
        ]

        tetra_msg = ['                   S T R E S S E S   I N    T E T R A H E D R O N   S O L I D   E L E M E N T S   ( C T E T R A )\n', ]
        penta_msg = ['                    S T R E S S E S   I N   P E N T A H E D R O N   S O L I D   E L E M E N T S   ( P E N T A )\n', ]
        hexa_msg = ['                      S T R E S S E S   I N   H E X A H E D R O N   S O L I D   E L E M E N T S   ( H E X A )\n', ]
        tetra_msg += base_msg
        penta_msg += base_msg
        hexa_msg += base_msg
        return tetra_msg, penta_msg, hexa_msg

class RealSolidStrainVector(RealSolidVector, StrainObject):
    def __init__(self, data_code, is_sort1, isubcase, dt):
        RealSolidVector.__init__(self, data_code, is_sort1, isubcase, dt)
        StrainObject.__init__(self, data_code, isubcase)

    def get_headers(self):
        if self.isVonMises():
            von_mises = 'von_mises'
        else:
            von_mises = 'max_shear'
        headers = ['exx', 'eyy', 'ezz', 'exy', 'eyz', 'exz', 'e1', 'e2', 'e3', von_mises]
        return headers

    def _get_msgs(self):
        if self.isVonMises():
            von_mises = 'VON MISES'
        else:
            von_mises = 'MAX SHEAR'

        base_msg = [
            '0                CORNER        ------CENTER AND CORNER POINT STRAINS---------       DIR.  COSINES       MEAN\n',
            '  ELEMENT-ID    GRID-ID        NORMAL              SHEAR             PRINCIPAL       -A-  -B-  -C-     PRESSURE       %s \n' % von_mises]
        tetra_msg = ['                     S T R A I N S   I N    T E T R A H E D R O N   S O L I D   E L E M E N T S   ( C T E T R A )\n', ]
        penta_msg = ['                      S T R A I N S   I N   P E N T A H E D R O N   S O L I D   E L E M E N T S   ( P E N T A )\n', ]
        hexa_msg = ['                      S T R A I N S   I N   H E X A H E D R O N   S O L I D   E L E M E N T S   ( H E X A )\n', ]
        tetra_msg += base_msg
        penta_msg += base_msg
        hexa_msg += base_msg
        return tetra_msg, penta_msg, hexa_msg


class RealSolidStressObject(StressObject):
    """
    ::

      # s_code=0
              N O N L I N E A R   S T R E S S E S   I N   T E T R A H E D R O N   S O L I D   E L E M E N T S   ( T E T R A )
      ELEMENT GRID/   POINT                         STRESSES/ TOTAL STRAINS                          EQUIVALENT EFF. STRAIN  EFF. CREEP
         ID   GAUSS     ID       X           Y           Z           XY          YZ          ZX        STRESS   PLAS/NLELAS   STRAIN
          3    GRID   CENTER  6.6667E+02  2.6667E+03  2.6667E+03 -1.3333E+03  2.6667E+03 -1.3333E+03  6.0000E+03  1.5000E-04   0.0
                              1.6667E-05  6.6667E-05  6.6667E-05 -6.6667E-05  1.3333E-04 -6.6667E-05

    ::

      # s_code=1
                            S T R E S S E S   I N   H E X A H E D R O N   S O L I D   E L E M E N T S   ( H E X A )
                     CORNER        ------CENTER AND CORNER POINT STRESSES---------       DIR.  COSINES       MEAN
      ELEMENT-ID    GRID-ID        NORMAL              SHEAR             PRINCIPAL       -A-  -B-  -C-     PRESSURE       VON MISES
          1               0GRID CS  8 GP
                     CENTER  X   4.499200E+02  XY  -5.544791E+02   A   1.000000E+04  LX 0.00 0.69-0.72  -3.619779E+03    9.618462E+03
                             Y   4.094179E+02  YZ   5.456968E-12   B  -1.251798E+02  LY 0.00 0.72 0.69
                             Z   1.000000E+04  ZX  -4.547474E-13   C   9.845177E+02  LZ 1.00 0.00 0.00
    """
    def __init__(self, data_code, is_sort1, isubcase, dt=None):
        StressObject.__init__(self, data_code, isubcase)

        self.eType = {}
        self.data = None
        self._f06_data = None
        self.code = [self.format_code, self.sort_code, self.s_code]

        self.cid = {}  # gridGauss
        self.oxx = {}
        self.oyy = {}
        self.ozz = {}
        self.txy = {}
        self.tyz = {}
        self.txz = {}
        self.o1 = {}
        self.o2 = {}
        self.o3 = {}
        self.ovmShear = {}

        self.dt = dt
        if is_sort1:
            if dt is not None:
                self.add_node = self.add_node_sort1
                self.add_eid = self.add_eid_sort1
        else:
            assert dt is not None
            raise NotImplementedError('SORT2')
            #self.add = self.addSort2
            #self.add_eid = self.add_eid_sort2

    def get_stats(self):
        nelements = len(self.eType)

        msg = []
        if self.nonlinear_factor is not None:  # transient
            ntimes = len(self.oxx)
            msg.append('  type=%s ntimes=%s nelements=%s\n'
                       % (self.__class__.__name__, ntimes, nelements))
        else:
            msg.append('  type=%s nelements=%s\n' % (self.__class__.__name__,
                                                     nelements))
        msg.append('  eType, cid, oxx, oyy, ozz, txy, tyz, txz, '
                   'o1, o2, o3, ovmShear\n  ')
        msg.append('elementTypes: %s\n  ' % ', '.join(set(self.eType.values())))
        msg += self.get_data_code()
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

            if dt not in self.data:
                self.data[dt] = []
            for line in data:
                self.data[dt] += data

    def processF06Data(self):
        """
        ::

                            S T R E S S E S   I N   P E N T A H E D R O N   S O L I D   E L E M E N T S   ( P E N T A )
                         CORNER        ------CENTER AND CORNER POINT STRESSES---------       DIR.  COSINES       MEAN
          ELEMENT-ID    GRID-ID        NORMAL              SHEAR             PRINCIPAL       -A-  -B-  -C-     PRESSURE       VON MISES
                  2           0GRID CS  6 GP
                         CENTER  X  -1.829319E+03  XY   7.883865E+01   A   1.033182E+04  LX-0.12 0.71 0.70  -2.115135E+03    1.232595E+04
                                 Y  -1.825509E+03  YZ  -1.415218E+03   B  -2.080181E+03  LY-0.12 0.69-0.71
                                 Z   1.000023E+04  ZX  -1.415218E+03   C  -1.906232E+03  LZ 0.99 0.16 0.00
        """
        eMap = {'CTETRA4': 5, 'CPENTA6': 7, 'CHEXA8': 9, 'HEXA20': 9,
                'PENTA15': 7, 'TETRA10': 5, }   # +1 for the centroid
        if self.nonlinear_factor is None:
            ipack = []
            i = 0
            n = 0
            while n < len(self.data):
                line = self.data[n]

                eType = line[0]
                eid = int(line[1])
                nNodes = eMap[eType]
                self.eType[eid] = eType
                self.oxx[eid] = {}
                self.oyy[eid] = {}
                self.ozz[eid] = {}
                self.txy[eid] = {}
                self.tyz[eid] = {}
                self.txz[eid] = {}
                self.o1[eid] = {}
                self.o2[eid] = {}
                self.o3[eid] = {}
                self.ovmShear[eid] = {}
                n += 1
                for j in xrange(nNodes):
                    (blank, nodeID, x, oxx, xy, txy, a, o1, lx, d1, d2, d3, pressure, ovmShear) = self.data[n]
                    (blank, blank, y, oyy, yz, tyz, b, o2, ly, d1, d2, d3, blank, blank) = self.data[n+1]
                    (blank, blank, z, ozz, zx, txz, c, o3, lz, d1, d2, d3, blank, blank) = self.data[n+2]
                    if    nodeID.strip() == 'CENTER':
                        nodeID = 'CENTER'
                    else:
                        nodeID = int(nodeID)

                    self.oxx[eid][nodeID] = float(oxx)
                    self.oyy[eid][nodeID] = float(oyy)
                    self.ozz[eid][nodeID] = float(ozz)
                    self.txy[eid][nodeID] = float(txy)
                    self.tyz[eid][nodeID] = float(tyz)
                    self.txz[eid][nodeID] = float(txz)
                    self.o1[eid][nodeID] = float(o1)
                    self.o2[eid][nodeID] = float(o2)
                    self.o3[eid][nodeID] = float(o3)
                    self.ovmShear[eid][nodeID] = float(ovmShear)
                    n += 3
            return
        return
        raise NotImplementedError()
        #(dtName, dt) = transient
        #self.data_code['name'] = dtName
        #if dt not in self.gridTypes:
            #self.add_new_transient()

        #for line in data:
            #(nodeID, grid_type, t1, t2, t3, r1, r2, r3) = line
            #self.gridTypes[dt][nodeID] = array([t1, t2, t3])
            #self.translations[dt][nodeID] = array([t1, t2, t3])
            #self.rotations[dt][nodeID] = array([r1, r2, r3])
        #del self.data

    def delete_transient(self, dt):
        del self.oxx[dt]
        del self.oyy[dt]
        del self.ozz[dt]
        del self.txy[dt]
        del self.tyz[dt]
        del self.txz[dt]
        del self.o1[dt]
        del self.o2[dt]
        del self.o3[dt]

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
        self.o1[dt] = {}
        self.o2[dt] = {}
        self.o3[dt] = {}

        #self.aCos[dt] = {}
        #self.bCos[dt] = {}
        #self.cCos[dt] = {}
        #self.pressure[dt] = {}
        self.ovmShear[dt] = {}

    def add_eid(self, eType, cid, dt, eid, nodeID, oxx, oyy, ozz, txy, tyz, txz, o1, o2, o3, aCos, bCos, cCos, pressure, ovm):
        #assert eid not in self.oxx
        assert cid >= 0
        assert eid >= 0
        self.eType[eid] = eType
        self.cid[eid] = cid
        self.oxx[eid] = {nodeID: oxx}
        self.oyy[eid] = {nodeID: oyy}
        self.ozz[eid] = {nodeID: ozz}
        self.txy[eid] = {nodeID: txy}
        self.tyz[eid] = {nodeID: tyz}
        self.txz[eid] = {nodeID: txz}
        self.o1[eid] = {nodeID: o1}
        self.o2[eid] = {nodeID: o2}
        self.o3[eid] = {nodeID: o3}
        #self.aCos[eid] = {nodeID: aCos}
        #self.bCos[eid] = {nodeID: bCos}
        #self.cCos[eid] = {nodeID: cCos}
        #self.pressure[eid] = {nodeID: pressure}
        self.ovmShear[eid] = {nodeID: ovm}
        msg = "*eid=%s nodeID=%s vm=%g" % (eid, nodeID, ovm)
        #print msg
        if nodeID == 0:
            raise ValueError(msg)

    def add_eid_sort1(self, eType, cid, dt, eid, nodeID, oxx, oyy, ozz, txy, tyz, txz, o1, o2, o3, aCos, bCos, cCos, pressure, ovm):
        assert cid >= 0
        assert eid >= 0

        if dt not in self.oxx:
            self.add_new_transient(dt)
        assert eid not in self.oxx[dt], self.oxx[dt]

        self.eType[eid] = eType
        self.cid[eid] = cid
        self.oxx[dt][eid] = {nodeID: oxx}
        self.oyy[dt][eid] = {nodeID: oyy}
        self.ozz[dt][eid] = {nodeID: ozz}
        self.txy[dt][eid] = {nodeID: txy}
        self.tyz[dt][eid] = {nodeID: tyz}
        self.txz[dt][eid] = {nodeID: txz}

        self.o1[dt][eid] = {nodeID: o1}
        self.o2[dt][eid] = {nodeID: o2}
        self.o3[dt][eid] = {nodeID: o3}

        #self.aCos[dt][eid] = {nodeID: aCos}
        #self.bCos[dt][eid] = {nodeID: bCos}
        #self.cCos[dt][eid] = {nodeID: cCos}
        #self.pressure[dt][eid] = {nodeID: pressure}
        self.ovmShear[dt][eid] = {nodeID: ovm}
        if nodeID == 0:
            msg = "*eid=%s nodeID=%s vm=%g" % (eid, nodeID, ovm)
            #print msg
            raise ValueError(msg)

    def add_node(self, dt, eid, inode, nodeID, oxx, oyy, ozz, txy, tyz, txz, o1, o2, o3, aCos, bCos, cCos, pressure, ovm):
        #msg = "eid=%s nodeID=%s vm=%g" % (eid, nodeID, ovm)
        #print msg
        self.oxx[eid][nodeID] = oxx
        self.oyy[eid][nodeID] = oyy
        self.ozz[eid][nodeID] = ozz

        self.txy[eid][nodeID] = txy
        self.tyz[eid][nodeID] = tyz
        self.txz[eid][nodeID] = txz

        self.o1[eid][nodeID] = o1
        self.o2[eid][nodeID] = o2
        self.o3[eid][nodeID] = o3

        #self.aCos[eid][nodeID] = aCos
        #self.bCos[eid][nodeID] = bCos
        #self.cCos[eid][nodeID] = cCos
        #self.pressure[eid][nodeID] = pressure
        self.ovmShear[eid][nodeID] = ovm
        if nodeID == 0:
            raise ValueError(msg)

    def add_node_sort1(self, dt, eid, inode, nodeID, oxx, oyy, ozz, txy, tyz, txz, o1, o2, o3, aCos, bCos, cCos, pressure, ovm):
        #self.eType[eid] = eType
        #print "eid=%s nid=%s oxx=%s" %(eid,nodeID,oxx)
        self.oxx[dt][eid][nodeID] = oxx
        self.oyy[dt][eid][nodeID] = oyy
        self.ozz[dt][eid][nodeID] = ozz

        self.txy[dt][eid][nodeID] = txy
        self.tyz[dt][eid][nodeID] = tyz
        self.txz[dt][eid][nodeID] = txz

        self.o1[dt][eid][nodeID] = o1
        self.o2[dt][eid][nodeID] = o2
        self.o3[dt][eid][nodeID] = o3

        #self.aCos[dt][eid][nodeID] = aCos
        #self.bCos[dt][eid][nodeID] = bCos
        #self.cCos[dt][eid][nodeID] = cCos
        #self.pressure[dt][eid][nodeID] = pressure
        self.ovmShear[dt][eid][nodeID] = ovm
        if nodeID == 0:
            msg = "eid=%s nodeID=%s vm=%g" % (eid, nodeID, ovm)
            #print msg
            raise ValueError(msg)

    def getHeaders(self):
        headers = ['oxx', 'oyy', 'ozz', 'txy', 'tyz', 'txz']
        if self.isVonMises():
            headers.append('oVonMises')
        else:
            headers.append('oMaxShear')
        return headers

    def pressure(self, o1, o2, o3):
        """
        returns the hydrostatic pressure
        (o1+o2+o3)/-3.
        http://en.wikipedia.org/wiki/Stress_%28mechanics%29
        """
        return (o1 + o2 + o3) / -3.

    def directionalVectors(self, oxx, oyy, ozz, txy, tyz, txz):
        A = [[oxx, txy, txz],
             [txy, oyy, tyz],
             [txz, tyz, ozz]]
        (Lambda, v) = eigh(A)  # a hermitian matrix is a symmetric-real matrix
        return v

    def getF06_Header(self):
        if self.isVonMises():
            vonMises = 'VON MISES'
        else:
            vonMises = 'MAX SHEAR'

        tetraMsg = ['                   S T R E S S E S   I N    T E T R A H E D R O N   S O L I D   E L E M E N T S   ( C T E T R A )\n',
                    '0                CORNER        ------CENTER AND CORNER POINT STRESSES---------       DIR.  COSINES       MEAN                   \n',
                    '  ELEMENT-ID    GRID-ID        NORMAL              SHEAR             PRINCIPAL       -A-  -B-  -C-     PRESSURE       %s \n' % (vonMises)]

        pentaMsg = ['                    S T R E S S E S   I N   P E N T A H E D R O N   S O L I D   E L E M E N T S   ( P E N T A )\n',
                    '0                CORNER        ------CENTER AND CORNER POINT STRESSES---------       DIR.  COSINES       MEAN                   \n',
                    '  ELEMENT-ID    GRID-ID        NORMAL              SHEAR             PRINCIPAL       -A-  -B-  -C-     PRESSURE       %s \n' % (vonMises)]

        hexaMsg = ['                      S T R E S S E S   I N   H E X A H E D R O N   S O L I D   E L E M E N T S   ( H E X A )\n',
                   '0                CORNER        ------CENTER AND CORNER POINT STRESSES---------       DIR.  COSINES       MEAN                   \n',
                   '  ELEMENT-ID    GRID-ID        NORMAL              SHEAR             PRINCIPAL       -A-  -B-  -C-     PRESSURE       %s \n' % (vonMises)]

        #nNodes = {'CTETRA':4,'CPENTA':6,'CHEXA':8,'HEXA':8,'PENTA':6,'TETRA':4,}
        TETRA = defaultdict(list)
        HEXA = defaultdict(list)
        PENTA = defaultdict(list)

        for eid, eType in sorted(self.eType.iteritems()):
            if 'CTETRA' in eType:
                TETRA[eType[6:]].append(eid)
            elif 'CHEXA' in eType:
                HEXA[eType[5:]].append(eid)
            elif 'CPENTA' in eType:
                PENTA[eType[6:]].append(eid)
            else:
                raise NotImplementedError('eType=%r' % eType)

        return (tetraMsg, pentaMsg, hexaMsg,
                TETRA, PENTA, HEXA)

    def write_f06(self, header, pageStamp, pageNum=1, f=None, is_mag_phase=False):
        assert f is not None
        if self.nonlinear_factor is not None:
            return self._write_f06_transient(header, pageStamp, pageNum, f)

        (tetraMsg, pentaMsg, hexaMsg,
         TETRA, PENTA, HEXA) = self.getF06_Header()
        #nNodes = {'CTETRA':4,'CPENTA':6,'CHEXA':8,'HEXA':8,'PENTA':6,'TETRA':4,}
        keys = [int(key) for key in TETRA.keys()]
        for key in sorted(keys):
            eids = TETRA[str(key)]
            self._write_element('CTETRA'+str(key), key, eids, header, tetraMsg, f)
            f.write(pageStamp % pageNum)
            pageNum += 1

        keys = [int(key) for key in PENTA.keys()]
        for key in sorted(keys):
            eids = PENTA[str(key)]
            self._write_element('CPENTA'+str(key), key, eids, header, pentaMsg, f)
            f.write(pageStamp % pageNum)
            pageNum += 1

        keys = [int(key) for key in HEXA.keys()]
        for key in sorted(keys):
            eids = HEXA[str(key)]
            self._write_element('CHEXA'+str(key), key, eids, header, hexaMsg, f)
            f.write(pageStamp % pageNum)
            pageNum += 1
        return pageNum - 1

    def _write_f06_transient(self, header, pageStamp, pageNum, f):
        (tetraMsg, pentaMsg, hexaMsg,
         TETRA, HEXA, PENTA) = self.getF06_Header()
        dts = self.oxx.keys()
        for dt in dts:
            keys = [int(key) for key in TETRA.keys()]
            for key in sorted(keys):
                eids = TETRA[str(key)]
                self._write_element_transient('CTETRA'+str(key), key, eids, dt, header, tetraMsg, f)
                f.write(pageStamp % pageNum)
                pageNum += 1

            keys = [int(key) for key in PENTA.keys()]
            for key in sorted(keys):
                eids = PENTA[str(key)]
                self._write_element_transient('CPENTA'+str(key), key, eids, dt, header, pentaMsg, f)
                f.write(pageStamp % pageNum)
                pageNum += 1

            keys = [int(key) for key in HEXA.keys()]
            for key in sorted(keys):
                eids = HEXA[str(key)]
                self._write_element_transient('CHEXA'+str(key), key, eids, dt, header, hexaMsg, f)
                f.write(pageStamp % pageNum)
                pageNum += 1

        return pageNum - 1

    def _write_element(self, eType, nNodes, eids, header, tetraMsg, f):
        f.write(''.join(header + tetraMsg))
        for eid in eids:
            #eType = self.eType[eid]

            k = self.oxx[eid].keys()
            cen = 'CENTER'
            k.remove(cen)
            k.sort()
            f.write('0  %8s           0GRID CS  %i GP\n' % (eid, nNodes))
            for nid in [cen] + k:
                oxx = self.oxx[eid][nid]
                oyy = self.oyy[eid][nid]
                ozz = self.ozz[eid][nid]
                txy = self.txy[eid][nid]
                tyz = self.tyz[eid][nid]
                txz = self.txz[eid][nid]

                o1 = self.o1[eid][nid]
                o2 = self.o2[eid][nid]
                o3 = self.o3[eid][nid]
                ovm = self.ovmShear[eid][nid]
                p = (o1 + o2 + o3) / -3.

                A = [[oxx, txy, txz],
                     [txy, oyy, tyz],
                     [txz, tyz, ozz]]
                (Lambda, v) = eigh(A)  # a hermitian matrix is a symmetric-real matrix

                ([oxx, oyy, ozz, txy, tyz, txz, o1, o2, o3, p, ovm], isAllZeros) = writeFloats13E([
                  oxx, oyy, ozz, txy, tyz, txz, o1, o2, o3, p, ovm])
                f.write('0              %8s  X  %13s  XY  %13s   A  %13s  LX%5.2f%5.2f%5.2f  %13s   %-s\n'
                        '               %8s  Y  %13s  YZ  %13s   B  %13s  LY%5.2f%5.2f%5.2f\n'
                        '               %8s  Z  %13s  ZX  %13s   C  %13s  LZ%5.2f%5.2f%5.2f\n' % (
                            nid, oxx, txy, o1, v[0, 1], v[0, 2], v[0, 0], p, ovm,
                            '', oyy, tyz, o2, v[1, 1], v[1, 2], v[1, 0],
                            '', ozz, txz, o3, v[2, 1], v[2, 2], v[2, 0]))

    def _write_element_transient(self, eType, nNodes, eids, dt, header, tetraMsg, f):
        dtLine = '%14s = %12.5E\n' % (self.data_code['name'], dt)
        header[1] = dtLine
        f.write(''.join(header + tetraMsg))
        cen = 'CENTER'
        for eid in eids:
            k = self.oxx[dt][eid].keys()
            k.remove(cen)
            k.sort()
            f.write('0  %8s           0GRID CS  %i GP\n' % (eid, nNodes))
            for nid in [cen] + k:
                oxx = self.oxx[dt][eid][nid]
                oyy = self.oyy[dt][eid][nid]
                ozz = self.ozz[dt][eid][nid]
                txy = self.txy[dt][eid][nid]
                tyz = self.tyz[dt][eid][nid]
                txz = self.txz[dt][eid][nid]

                o1 = self.o1[dt][eid][nid]
                o2 = self.o2[dt][eid][nid]
                o3 = self.o3[dt][eid][nid]
                ovm = self.ovmShear[dt][eid][nid]
                p = (o1 + o2 + o3) / -3.

                #print "nid = |%r|" %(nid)
                A = [[oxx, txy, txz],
                     [txy, oyy, tyz],
                     [txz, tyz, ozz]]
                (Lambda, v) = eigh(A)  # a hermitian matrix is a symmetric-real matrix

                ([oxx, oyy, ozz, txy, tyz, txz, o1, o2, o3, p, ovm], isAllZeros) = writeFloats13E([
                  oxx, oyy, ozz, txy, tyz, txz, o1, o2, o3, p, ovm])
                f.write('0              %8s  X  %13s  XY  %13s   A  %13s  LX%5.2f%5.2f%5.2f  %13s   %s\n'
                        '               %8s  Y  %13s  YZ  %13s   B  %13s  LY%5.2f%5.2f%5.2f\n'
                        '               %8s  Z  %13s  ZX  %13s   C  %13s  LZ%5.2f%5.2f%5.2f\n'
                        % (nid, oxx, txy, o1, v[0, 1], v[0, 2], v[0, 0], p, ovm,
                        '', oyy, tyz, o2, v[1, 1], v[1, 2], v[1, 0],
                        '', ozz, txz, o3, v[2, 1], v[2, 2], v[2, 0]))


class RealSolidStrainObject(StrainObject):
    """
    ::

      # code=[1,0,11]
                            S T R A I N S   I N   H E X A H E D R O N   S O L I D   E L E M E N T S   ( H E X A )
                     CORNER        ------CENTER AND CORNER POINT STRESSES---------       DIR.  COSINES       MEAN
      ELEMENT-ID    GRID-ID        NORMAL              SHEAR             PRINCIPAL       -A-  -B-  -C-     PRESSURE       VON MISES
              1           0GRID CS  8 GP
                     CENTER  X   4.499200E+02  XY  -5.544791E+02   A   1.000000E+04  LX 0.00 0.69-0.72  -3.619779E+03    9.618462E+03
                             Y   4.094179E+02  YZ   5.456968E-12   B  -1.251798E+02  LY 0.00 0.72 0.69
                             Z   1.000000E+04  ZX  -4.547474E-13   C   9.845177E+02  LZ 1.00 0.00 0.00

      # code=[1,0,10]
                         S T R A I N S   I N    T E T R A H E D R O N   S O L I D   E L E M E N T S   ( C T E T R A )
                     CORNER        ------CENTER AND CORNER POINT  STRAINS---------       DIR.  COSINES       MEAN         OCTAHEDRAL
      ELEMENT-ID    GRID-ID        NORMAL              SHEAR             PRINCIPAL       -A-  -B-  -C-     PRESSURE         SHEAR
              4           0GRID CS  4 GP
                     CENTER  X  -2.288232E-04  XY   1.240506E-04   A   9.631978E-04  LX-0.10-0.71-0.70  -1.601805E-04    5.692614E-04
                             Y  -2.289814E-04  YZ  -2.369997E-04   B  -2.909276E-04  LY-0.10 0.71-0.70
                             Z   9.383460E-04  ZX  -2.369997E-04   C  -1.917288E-04  LZ 0.99 0.00-0.15
    """
    def __init__(self, data_code, is_sort1, isubcase, dt):
        StrainObject.__init__(self, data_code, isubcase)
        self.eType = {}
        self.data = None
        self._f06_data = None
        self.code = [self.format_code, self.sort_code, self.s_code]

        self.cid = {}
        self.exx = {}
        self.eyy = {}
        self.ezz = {}
        self.exy = {}
        self.eyz = {}
        self.exz = {}
        self.e1 = {}
        self.e2 = {}
        self.e3 = {}
        #self.aCos = {}
        #self.bCos = {}
        #self.cCos = {}
        #self.pressure = {}
        self.evmShear = {}

        self.nonlinear_factor = dt
        if is_sort1:
            if dt is not None:
                print("dt =", dt)
                self.add_node = self.add_node_sort1
                self.add_eid = self.add_eid_sort1
        else:
            assert dt is not None
            self.add_node = self.add_node_sort2
            self.add_eid = self.add_eid_sort2

    def get_stats(self):
        nelements = len(self.eType)

        msg = []
        if self.nonlinear_factor is not None:  # transient
            ntimes = len(self.exx)
            msg.append('  type=%s ntimes=%s nelements=%s\n'
                       % (self.__class__.__name__, ntimes, nelements))
        else:
            msg.append('  type=%s nelements=%s\n' % (self.__class__.__name__,
                                                     nelements))
        msg.append('  eType, cid, exx, eyy, ezz, exy, eyz, exz, '
                   'e1, e2, e3, evmShear\n')
        msg.append('  elementTypes: %s\n  ' % ', '.join(set(self.eType.values())))
        msg += self.get_data_code()
        return msg

    def add_f06_data(self, data, transient):
        if transient is None:
            if not hasattr(self, '_f06_data'):
                self._f06_data = []
            self._f06_data += data
        else:
            if not hasattr(self, '_f06_data'):
                self._f06_data = {}
            dt = transient[1]
            if dt not in self.data:
                self._f06_data[dt] = []
            for line in data:
                self._f06_data[dt] += data

    def processF06Data(self):
        """
        ::

                            S T R E S S E S   I N   P E N T A H E D R O N   S O L I D   E L E M E N T S   ( P E N T A )
                         CORNER        ------CENTER AND CORNER POINT STRESSES---------       DIR.  COSINES       MEAN
          ELEMENT-ID    GRID-ID        NORMAL              SHEAR             PRINCIPAL       -A-  -B-  -C-     PRESSURE       VON MISES
                  2           0GRID CS  6 GP
                         CENTER  X  -1.829319E+03  XY   7.883865E+01   A   1.033182E+04  LX-0.12 0.71 0.70  -2.115135E+03    1.232595E+04
                                 Y  -1.825509E+03  YZ  -1.415218E+03   B  -2.080181E+03  LY-0.12 0.69-0.71
                                 Z   1.000023E+04  ZX  -1.415218E+03   C  -1.906232E+03  LZ 0.99 0.16 0.00
        """
        eMap = {'CTETRA': 5, 'CPENTA': 7, 'CHEXA': 9, 'HEXA': 9,
                'PENTA': 7, 'TETRA': 5, }   # +1 for the centroid
        if self.nonlinear_factor is None:
            pack = []
            i = 0
            n = 0
            while n < len(self.data):
                line = self.data[n]

                eType = line[0]
                eid = int(line[1])
                nNodes = eMap[eType]
                self.eType[eid] = eType
                self.exx[eid] = {}
                self.eyy[eid] = {}
                self.ezz[eid] = {}
                self.exy[eid] = {}
                self.eyz[eid] = {}
                self.exz[eid] = {}
                self.e1[eid] = {}
                self.e2[eid] = {}
                self.e3[eid] = {}
                self.evmShear[eid] = {}
                n += 1
                for j in xrange(nNodes):
                    (blank, nodeID, x, exx, xy, exy, a, e1, lx, d1, d2, d3,
                        pressure, evmShear) = self.data[n]
                    (blank, blank, y, eyy, yz, eyz, b, e2, ly, d1, d2, d3,
                        blank, blank) = self.data[n + 1]
                    (blank, blank, z, ezz, zx, exz, c, e3, lz, d1, d2, d3,
                        blank, blank) = self.data[n + 2]
                    if    nodeID.strip() == 'CENTER':
                        nodeID = 'CENTER'
                    else:
                        nodeID = int(nodeID)
                    self.exx[eid][nodeID] = float(exx)
                    self.eyy[eid][nodeID] = float(eyy)
                    self.ezz[eid][nodeID] = float(ezz)
                    self.exy[eid][nodeID] = float(exy)
                    self.eyz[eid][nodeID] = float(eyz)
                    self.exz[eid][nodeID] = float(exz)
                    self.e1[eid][nodeID] = float(e1)
                    self.e2[eid][nodeID] = float(e2)
                    self.e3[eid][nodeID] = float(e3)
                    self.evmShear[eid][nodeID] = float(evmShear)
                    n += 3
            return
        raise NotImplementedError()
        (dtName, dt) = transient
        self.data_code['name'] = dtName
        if dt not in self.gridTypes:
            self.add_new_transient()

        for line in data:
            (nodeID, grid_type, t1, t2, t3, r1, r2, r3) = line
            self.gridTypes[dt][nodeID] = array([t1, t2, t3])
            self.translations[dt][nodeID] = array([t1, t2, t3])
            self.rotations[dt][nodeID] = array([r1, r2, r3])
        del self._f06_data

    def delete_transient(self, dt):
        del self.exx[dt]
        del self.eyy[dt]
        del self.ezz[dt]
        del self.exy[dt]
        del self.eyz[dt]
        del self.exz[dt]
        del self.e1[dt]
        del self.e2[dt]
        del self.e3[dt]

    def get_transients(self):
        k = self.exx.keys()
        k.sort()
        return k

    def add_new_transient(self, dt):
        """
        initializes the transient variables
        """
        self.nonlinear_factor = dt
        self.exx[dt] = {}
        self.eyy[dt] = {}
        self.ezz[dt] = {}
        self.exy[dt] = {}
        self.eyz[dt] = {}
        self.exz[dt] = {}
        self.e1[dt] = {}
        self.e2[dt] = {}
        self.e3[dt] = {}
        #self.aCos[dt] = {}
        #self.bCos[dt] = {}
        #self.cCos[dt] = {}
        #self.pressure[dt] = {}
        self.evmShear[dt] = {}

    def add_eid(self, eType, cid, dt, eid, nodeID, exx, eyy, ezz, exy, eyz, exz, e1, e2, e3, aCos, bCos, cCos, pressure, evm):
        assert eid not in self.exx
        assert cid >= 0
        assert eid >= 0
        self.eType[eid] = eType
        self.cid[eid] = cid
        self.exx[eid] = {nodeID: exx}
        self.eyy[eid] = {nodeID: eyy}
        self.ezz[eid] = {nodeID: ezz}
        self.exy[eid] = {nodeID: exy}
        self.eyz[eid] = {nodeID: eyz}
        self.exz[eid] = {nodeID: exz}
        self.e1[eid] = {nodeID: e1}
        self.e2[eid] = {nodeID: e2}
        self.e3[eid] = {nodeID: e3}
        #self.aCos[eid] = {nodeID: aCos}
        #self.bCos[eid] = {nodeID: bCos}
        #self.cCos[eid] = {nodeID: cCos}
        #self.pressure[eid] = {nodeID: pressure}
        self.evmShear[eid] = {nodeID: evm}
        if nodeID == 0:
            msg = "*eid=%s nodeID=%s evmShear=%g" % (eid, nodeID, evm)
            #print msg
            raise RuntimeError(msg)

    def add_eid_sort1(self, eType, cid, dt, eid, nodeID, exx, eyy, ezz, exy, eyz, exz, e1, e2, e3, aCos, bCos, cCos, pressure, evm):
        assert cid >= 0
        assert eid >= 0
        if dt not in self.exx:
            self.add_new_transient(dt)

        assert eid not in self.exx[dt], self.exx[dt]
        self.eType[eid] = eType
        self.cid[eid] = cid
        self.exx[dt][eid] = {nodeID: exx}
        self.eyy[dt][eid] = {nodeID: eyy}
        self.ezz[dt][eid] = {nodeID: ezz}
        self.exy[dt][eid] = {nodeID: exy}
        self.eyz[dt][eid] = {nodeID: eyz}
        self.exz[dt][eid] = {nodeID: exz}
        self.e1[dt][eid] = {nodeID: e1}
        self.e2[dt][eid] = {nodeID: e2}
        self.e3[dt][eid] = {nodeID: e3}
        #self.aCos[dt][eid] = {nodeID: aCos}
        #self.bCos[dt][eid] = {nodeID: bCos}
        #self.cCos[dt][eid] = {nodeID: cCos}
        #self.pressure[dt][eid] = {nodeID: pressure}
        self.evmShear[dt][eid] = {nodeID: evm}
        if nodeID == 0:
            msg = "*eid=%s nodeID=%s evmShear=%g" % (eid, nodeID, evm)
            #print msg
            raise Exception(msg)

    def add_node(self, dt, eid, inode, nodeID, exx, eyy, ezz, exy, eyz, exz, e1, e2, e3, aCos, bCos, cCos, pressure, evm):
        self.exx[eid][nodeID] = exx
        self.eyy[eid][nodeID] = eyy
        self.ezz[eid][nodeID] = ezz

        self.exy[eid][nodeID] = exy
        self.eyz[eid][nodeID] = eyz
        self.exz[eid][nodeID] = exz
        self.e1[eid][nodeID] = e1
        self.e2[eid][nodeID] = e2
        self.e3[eid][nodeID] = e3

        #self.aCos[eid][nodeID] = aCos
        #self.bCos[eid][nodeID] = bCos
        #self.cCos[eid][nodeID] = cCos
        #self.pressure[eid][nodeID] = pressure
        self.evmShear[eid][nodeID] = evm
        if nodeID == 0:
            msg = "eid=%s nodeID=%s vm=%g" % (eid, nodeID, evm)
            #print msg
            raise RuntimeError(msg)

    def add_node_sort1(self, dt, eid, inode, nodeID, exx, eyy, ezz, exy, eyz, exz, e1, e2, e3, aCos, bCos, cCos, pressure, evm):
        self.exx[dt][eid][nodeID] = exx
        self.eyy[dt][eid][nodeID] = eyy
        self.ezz[dt][eid][nodeID] = ezz

        self.exy[dt][eid][nodeID] = exy
        self.eyz[dt][eid][nodeID] = eyz
        self.exz[dt][eid][nodeID] = exz

        self.e1[dt][eid][nodeID] = e1
        self.e2[dt][eid][nodeID] = e2
        self.e3[dt][eid][nodeID] = e3

        #self.aCos[dt][eid][nodeID] = aCos
        #self.bCos[dt][eid][nodeID] = bCos
        #self.cCos[dt][eid][nodeID] = cCos
        #self.pressure[dt][eid][nodeID] = pressure
        self.evmShear[dt][eid][nodeID] = evm

        if nodeID == 0:
            msg = "eid=%s nodeID=%s vm=%g" % (eid, nodeID, evm)
            #print msg
            raise Exception(msg)

    def getHeaders(self):
        headers = ['exx', 'eyy', 'ezz', 'exy', 'eyz', 'exz']
        if self.isVonMises():
            headers.append('evm')
        else:
            headers.append('maxShear')
        return headers

    def ovm(self, o11, o22, o33, o12, o13, o23):
        """http://en.wikipedia.org/wiki/Von_Mises_yield_criterion"""
        ovm = 0.5 * ((o11 - o22) ** 2 +
                     (o22 - o33) ** 2 +
                     (o11 - o33) ** 2 +
                     6 * (o23 ** 2 + o13 ** 2 + o12 ** 2))
        return ovm

    #def ovmPlane(self,o11,o22,o12):
        #"""http://en.wikipedia.org/wiki/Von_Mises_yield_criterion"""
        #ovm = sqrt(o11**2+o22**2-o1*o2+3*o12**2)
        #return ovm

    def octahedral(self, o11, o22, o33, o12, o13, o23):
        """http://en.wikipedia.org/wiki/Von_Mises_yield_criterion"""
        ovm = self.ovm(o11, o22, o33, o12, o13, o23)
        return ovm * sqrt(2) / 3.

    def pressure(self, e1, e2, e3):
        """
        returns the hydrostatic pressure
        (e1+e2+e3)/-3.
        http://en.wikipedia.org/wiki/Stress_%28mechanics%29
        """
        return (e1 + e2 + e3) / -3.

    def directionalVectors(self, exx, eyy, ezz, exy, eyz, exz):
        A = [[exx, exy, exz],
             [exy, eyy, eyz],
             [exz, eyz, ezz]]
        (Lambda, v) = eigh(A)  # a hermitian matrix is a symmetric-real matrix
        return v

    def getF06_Header(self):
        if self.isVonMises():
            vonMises = 'VON MISES'
        else:
            vonMises = 'MAX SHEAR'

        tetraMsg = ['                     S T R A I N S   I N    T E T R A H E D R O N   S O L I D   E L E M E N T S   ( C T E T R A )\n',
                    '0                CORNER        ------CENTER AND CORNER POINT STRAINS---------       DIR.  COSINES       MEAN\n',
                    '  ELEMENT-ID    GRID-ID        NORMAL              SHEAR             PRINCIPAL       -A-  -B-  -C-     PRESSURE       %s \n' % (vonMises)]

        pentaMsg = ['                      S T R A I N S   I N   P E N T A H E D R O N   S O L I D   E L E M E N T S   ( P E N T A )\n',
                    '0                CORNER        ------CENTER AND CORNER POINT STRAINS---------       DIR.  COSINES       MEAN\n',
                    '  ELEMENT-ID    GRID-ID        NORMAL              SHEAR             PRINCIPAL       -A-  -B-  -C-     PRESSURE       %s \n' % (vonMises)]

        hexaMsg = ['                      S T R A I N S   I N   H E X A H E D R O N   S O L I D   E L E M E N T S   ( H E X A )\n',
                   '0                CORNER        ------CENTER AND CORNER POINT STRAINS---------       DIR.  COSINES       MEAN\n',
                   '  ELEMENT-ID    GRID-ID        NORMAL              SHEAR             PRINCIPAL       -A-  -B-  -C-     PRESSURE       %s \n' % (vonMises)]

        #nNodes = {'CTETRA':4,'CPENTA':6,'CHEXA':8,'HEXA':8,'PENTA':6,'TETRA':4,}
        TETRA = defaultdict(list)
        HEXA = defaultdict(list)
        PENTA = defaultdict(list)

        for eid, eType in sorted(self.eType.iteritems()):
            if 'CTETRA' in eType:
                TETRA[eType[6:]].append(eid)
            elif 'CHEXA' in eType:
                HEXA[eType[5:]].append(eid)
            elif 'CPENTA' in eType:
                PENTA[eType[6:]].append(eid)
            else:
                raise NotImplementedError('eType=%r' % eType)

        return (tetraMsg, pentaMsg, hexaMsg,
                TETRA, PENTA, HEXA)

    def write_f06(self, header, pageStamp, pageNum=1, f=None, is_mag_phase=False):
        assert f is not None
        if self.nonlinear_factor is not None:
            return self._write_f06_transient(header, pageStamp, pageNum, f)

        (tetraMsg, pentaMsg, hexaMsg,
         TETRA, PENTA, HEXA) = self.getF06_Header()
        #nNodes = {'CTETRA':4,'CPENTA':6,'CHEXA':8,'HEXA':8,'PENTA':6,'TETRA':4,}
        keys = [int(key) for key in TETRA.keys()]
        for key in sorted(keys):
            eids = TETRA[str(key)]
            self._write_element('CTETRA'+str(key), key, eids, header, tetraMsg, f)
            f.write(pageStamp % pageNum)
            pageNum += 1

        keys = [int(key) for key in PENTA.keys()]
        for key in sorted(keys):
            eids = PENTA[str(key)]
            self._write_element('CPENTA'+str(key), key, eids, header, pentaMsg, f)
            f.write(pageStamp % pageNum)
            pageNum += 1

        keys = [int(key) for key in HEXA.keys()]
        for key in sorted(keys):
            eids = HEXA[str(key)]
            self._write_element('CHEXA'+str(key), key, eids, header, hexaMsg, f)
            f.write(pageStamp % pageNum)
            pageNum += 1
        return pageNum - 1

    def _write_f06_transient(self, header, pageStamp, pageNum=1, f=None, is_mag_phase=False):
        msg = []
        (tetraMsg, pentaMsg, hexaMsg,
         TETRA, PENTA, HEXA) = self.getF06_Header()
        dts = self.exx.keys()
        for dt in dts:
            keys = [int(key) for key in TETRA.keys()]
            for key in sorted(keys):
                eids = TETRA[str(key)]
                self._write_element_transient('CTETRA'+str(key), key, eids, dt, header, tetraMsg, f)
                f.write(pageStamp % pageNum)
                pageNum += 1

            keys = [int(key) for key in PENTA.keys()]
            for key in sorted(keys):
                eids = PENTA[str(key)]
                self._write_element_transient('CPENTA'+str(key), key, eids, dt, header, pentaMsg, f)
                f.write(pageStamp % pageNum)
                pageNum += 1

            keys = [int(key) for key in HEXA.keys()]
            for key in sorted(keys):
                eids = HEXA[str(key)]
                self._write_element_transient('CHEXA'+str(key), key, eids, dt, header, hexaMsg, f)
                f.write(pageStamp % pageNum)
                pageNum += 1

        return pageNum - 1

    def _write_element(self, eType, nNodes, eids, header, tetraMsg, f):
        f.write(''.join(header + tetraMsg))
        cen = 'CENTER'
        for eid in eids:
            k = self.exx[eid].keys()
            k.remove(cen)
            k.sort()
            f.write('0  %8s           0GRID CS  %i GP\n' % (eid, nNodes))
            for nid in [cen] + k:
                exx = self.exx[eid][nid]
                eyy = self.eyy[eid][nid]
                ezz = self.ezz[eid][nid]
                exy = self.exy[eid][nid]
                eyz = self.eyz[eid][nid]
                exz = self.exz[eid][nid]

                e1 = self.e1[eid][nid]
                e2 = self.e2[eid][nid]
                e3 = self.e3[eid][nid]
                evm = self.evmShear[eid][nid]
                p = (e1 + e2 + e3) / -3.

                A = [[exx, exy, exz],
                     [exy, eyy, eyz],
                     [exz, eyz, ezz]]
                (Lambda, v) = eigh(A)  # a hermitian matrix is a symmetric-real matrix

                ([exx, eyy, ezz, exy, eyz, exz, e1, e2, e3, p, evm], isAllZeros) = writeFloats13E([exx, eyy, ezz, exy, eyz, exz, e1, e2, e3, p, evm])
                f.write('0              %8s  X  %13s  XY  %13s   A  %13s  LX%5.2f%5.2f%5.2f  %13s   %-s\n'
                           '               %8s  Y  %13s  YZ  %13s   B  %13s  LY%5.2f%5.2f%5.2f\n'
                           '               %8s  Z  %13s  ZX  %13s   C  %13s  LZ%5.2f%5.2f%5.2f\n'
                           % (nid, exx, exy, e1, v[0, 1], v[0, 2], v[0, 0], p, evm,
                              '', eyy, eyz, e2, v[1, 1], v[1, 2], v[1, 0],
                              '', ezz, exz, e3, v[2, 1], v[2, 2], v[2, 0]), )

    def _write_element_transient(self, eType, nNodes, eids, dt, header, tetraMsg, f):
        dtLine = '%14s = %12.5E\n' % (self.data_code['name'], dt)
        header[1] = dtLine
        f.write(''.join(header + tetraMsg))
        cen = 'CENTER'
        for eid in eids:
            k = self.exx[dt][eid].keys()
            k.remove(cen)
            k.sort()
            f.write('0  %8s           0GRID CS  %i GP\n' % (eid, nNodes))
            for nid in [cen] + k:
                exx = self.exx[dt][eid][nid]
                eyy = self.eyy[dt][eid][nid]
                ezz = self.ezz[dt][eid][nid]
                exy = self.exy[dt][eid][nid]
                eyz = self.eyz[dt][eid][nid]
                exz = self.exz[dt][eid][nid]

                e1 = self.e1[dt][eid][nid]
                e2 = self.e2[dt][eid][nid]
                e3 = self.e3[dt][eid][nid]
                evm = self.evmShear[dt][eid][nid]
                p = (e1 + e2 + e3) / -3.

                A = [[exx, exy, exz],
                     [exy, eyy, eyz],
                     [exz, eyz, ezz]]
                (Lambda, v) = eigh(A)  # a hermitian matrix is a symmetric-real matrix

                ([exx, eyy, ezz, exy, eyz, exz, e1, e2, e3, p, evm], isAllZeros) = writeFloats13E([exx, eyy, ezz, exy, eyz, exz, e1, e2, e3, p, evm])
                f.write('0              %8s  X  %13s  XY  %13s   A  %13s  LX%5.2f%5.2f%5.2f  %13s   %-s\n'
                         '               %8s  Y  %13s  YZ  %13s   B  %13s  LY%5.2f%5.2f%5.2f\n'
                         '               %8s  Z  %13s  ZX  %13s   C  %13s  LZ%5.2f%5.2f%5.2f\n' % (
                             nid, exx, exy, e1, v[0, 1], v[0, 2], v[0, 0], p, evm,
                             '', eyy, eyz, e2, v[1, 1], v[1, 2], v[1, 0],
                             '', ezz, exz, e3, v[2, 1], v[2, 2], v[2, 0]))
