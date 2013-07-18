from __future__ import (nested_scopes, generators, division, absolute_import,
                        print_function, unicode_literals)
from numpy import array, sqrt, ones, zeros
from numpy.linalg import eigh
import pandas as pd

from .oes_objects import StressObject, StrainObject
from pyNastran.f06.f06_formatting import writeFloats13E


class SolidStressObject(StressObject):
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
    def __init__(self, data_code, is_sort1, isubcase, dt, read_mode):
        StressObject.__init__(self, data_code, isubcase, read_mode)

        self.shape = {}
        self._ncount = 0
        self._inode_start = None
        self._inode_end = None
        self._ielement_start = None
        self._ielement_end = None

        if read_mode == 0:
            return
        self.code = [self.format_code, self.sort_code, self.s_code]

        self.dt = dt
        #if is_sort1:
            #if dt is not None:
                #self.add = self.add_sort1
                #self.add_new_eid = self.add_new_eid_sort1
        #else:
            #assert dt is not None
            #self.add = self.addSort2
            #self.add_new_eid = self.add_new_eid_sort2

    def add_f06_data(self, data, transient):
        if transient is None:
            if not hasattr(self, 'data'):
                self.data = []
            self.data += data
        else:
            if not hasattr(self, 'data'):
                self.data = {}
            if dt not in self.data:
                self.data[dt] = []
            for line in data:
                self.data[dt] += data

    def process_f06_data(self):
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
            ipack = []
            i = 0
            n = 0
            while n < len(self.data):
                line = self.data[n]
                #print n,line

                eType = line[0]
                eid = int(line[1])
                #print "eType = ",eType
                nNodes = eMap[eType]
                #nodeID = None
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
                    #print self.data[n]
                    (blank,nodeID,x,oxx,xy,txy,a,o1,lx,d1,d2,d3,pressure,ovmShear) = self.data[n]
                    (blank,blank, y,oyy,yz,tyz,b,o2,ly,d1,d2,d3,blank,blank) = self.data[n+1]
                    (blank,blank, z,ozz,zx,txz,c,o3,lz,d1,d2,d3,blank,blank) = self.data[n+2]
                    if    nodeID.strip() == 'CENTER':
                        nodeID = 'C'
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
        raise NotImplementedError()

        (dtName, dt) = transient
        self.data_code['name'] = dtName
        if dt not in self.gridTypes:
            self.add_new_transient()

        for line in data:
            (nodeID, gridType, t1, t2, t3, r1, r2, r3) = line
            self.gridTypes[dt][nodeID] = array([t1, t2, t3])
            self.translations[dt][nodeID] = array([t1, t2, t3])
            self.rotations[dt][nodeID] = array([r1, r2, r3])
        del self.data

    def _increase_size(self, dt, nelements, nnodes):
        if dt in self.shape:  # default dictionary
            self.shape[dt][0] += nelements
            self.shape[dt][1] += nnodes
        else:
            self.shape[dt] = [nelements, nnodes]
        print("shape =", self.shape)

    def _get_shape(self):
        ndt = len(self.shape)
        dts = self.shape.keys()
        shape0 = dts[0]
        nelements = self.shape[shape0][0]
        nnodes = self.shape[shape0][1]
        print("ndt=%s nnodes=%s dts=%s" % (str(ndt), str(nnodes), str(dts)))
        return ndt, nelements, nnodes, dts

    def _increment(self, nnodes, nelements):
        self._inode_start += nnodes
        self._inode_end += nnodes
        self._ielement_start += nelements
        self._ielement_end += nelements
        return self._inode_start, self._inode_end, self._ielement_start, self._ielement_end

    def _preallocate(self, dt, nnodes, nelements):
        ndt, nelements_size, nnodes_size, dts = self._get_shape()
        print("ndt=%s nelements_size=%s nnodes_size=%s dts=%s" % (ndt, nelements_size, nnodes_size, str(dts)))

        if self._inode_start is not None:
            return (self._inode_start, self._inode_start + nnodes,
                    self._ielement_start, self._ielement_start + nelements)
        print('----definition----')
        n = ndt * nnodes_size
        if self._ncount != 0:
            asfd
        self._ncount += 1
        self._inode_start = 0
        self._inode_end = nnodes

        self._ielement_start = 0
        self._ielement_end = nelements

        data = {}
        element_data = {}
        columns = []
        if dts[0] is not None:
            if isinstance(dt, int):
                data['dt'] = pd.Series(zeros((n), dtype='int32'))
            else:
                data['dt'] = pd.Series(zeros((n), dtype='float32'))
            columns.append('dt')

        element_data['element_type'] = pd.Series(zeros(nelements_size, dtype='str'))
        element_data['element_id'] = pd.Series(zeros((nelements_size), dtype='int32'))

        data['element_id'] = pd.Series(zeros((n), dtype='int32'))
        data['node_id'] = pd.Series(zeros((n), dtype='int32'))

        #columns.append('element_type')
        columns.append('element_id')
        columns.append('node_id')

        #data['grid_type'] = pd.Series(zeros(ndt), dtype='int32'))
        #data['grid_type_str'] = pd.Series(zeros(nnodes), dtype='str'))
        print('n =', n)
        data['oxx'] = pd.Series(zeros((n), dtype='float32'))
        data['oyy'] = pd.Series(zeros((n), dtype='float32'))
        data['ozz'] = pd.Series(zeros((n), dtype='float32'))

        data['txy'] = pd.Series(zeros((n), dtype='float32'))
        data['txz'] = pd.Series(zeros((n), dtype='float32'))
        data['tyz'] = pd.Series(zeros((n), dtype='float32'))

        data['o1'] = pd.Series(zeros((n), dtype='float32'))
        data['o2'] = pd.Series(zeros((n), dtype='float32'))
        data['o3'] = pd.Series(zeros((n), dtype='float32'))

        if self.isVonMises():
            data['ovm'] = pd.Series(zeros((n), dtype='float32'))
        else:
            data['max_shear'] = pd.Series(zeros((n), dtype='float32'))

        #pressure
        #aCos
        #bCos
        #cCos
        columns += ['oxx', 'oyy', 'ozz', 'txy', 'txz', 'tyz', 'o1', 'o2', 'o3']

        if self.isVonMises():
            key = 'ovm'
        else:
            key = 'max_shear'
        data[key] = pd.Series(zeros((ndt * nnodes), dtype='float32'))
        columns.append(key)

        self.data = pd.DataFrame(data, columns=columns)
        self.element_data = pd.DataFrame(element_data, columns=['element_id', 'element_type', 'cid'])
        return (self._inode_start, self._inode_end, self._ielement_start, self._ielement_end)

    def _finalize(self, dt):
        ndt, nelements, nnodes, dts = self._get_shape()

        if dt != max(dts):
            return
        #print(self.data.to_string())
        #print("----finalize----")

        #grid_type_str = []
        #for grid_type in self.grid_type:
            #grid_type_str.append('C' if grid_type==0 else grid_type)
        #self.grid_type_str = pd.Series(grid_type_str, dtype='str')

        if dts[0] is not None:
            self.data = self.data.set_index(['dt', 'element_id', 'node_id'])
        else:
            self.data = self.data.set_index(['element_id', 'node_id'])
        #print "final\n", self.data
        del self._inode_start
        del self._inode_end

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
        tetraEids = []
        hexaEids = []
        pentaEids = []
        for eid, eType in sorted(self.eType.iteritems()):
            if eType == 'CTETRA' or eType == 'TETRA':
                tetraEids.append(eid)
            elif eType == 'CPENTA' or eType == 'PENTA':
                pentaEids.append(eid)
            elif eType == 'CHEXA' or eType == 'HEXA':
                hexaEids.append(eid)
            else:
                raise NotImplementedError('eType=|%r|' % (eType))
        return (tetraMsg, pentaMsg, hexaMsg, tetraEids, hexaEids, pentaEids)

    def write_f06(self, header, pageStamp, pageNum=1, f=None, is_mag_phase=False):
        if self.nonlinear_factor is not None:
            return self._write_f06_transient(header, pageStamp, pageNum, f)
        msg = []

        (tetraMsg, pentaMsg, hexaMsg, tetraEids, hexaEids,
            pentaEids) = self.getF06_Header()
        #nNodes = {'CTETRA':4,'CPENTA':6,'CHEXA':8,'HEXA':8,'PENTA':6,'TETRA':4,}
        if tetraEids:
            msg += self.writeElement(
                'CTETRA', 4, tetraEids, header, tetraMsg, f)
            msg.append(pageStamp + str(pageNum) + '\n')
            pageNum += 1
        if hexaEids:
            msg += self.writeElement('CHEXA', 8, hexaEids, header, hexaMsg, f)
            msg.append(pageStamp + str(pageNum) + '\n')
            pageNum += 1
        if pentaEids:
            msg += self.writeElement(
                'CPENTA', 6, pentaEids, header, pentaMsg, f)
            msg.append(pageStamp + str(pageNum) + '\n')
            pageNum += 1

        return (''.join(msg), pageNum - 1)

    def _write_f06_transient(self, header, pageStamp, pageNum=1, f=None, is_mag_phase=False):
        msg = []
        (tetraMsg, pentaMsg, hexaMsg, tetraEids, hexaEids,
            pentaEids) = self.getF06_Header()
        dts = self.oxx.keys()
        for dt in dts:
            if tetraEids:
                msg += self.writeElementTransient('CTETRA',
                                                  4, tetraEids, dt, header, tetraMsg, f)
                msg.append(pageStamp + str(pageNum) + '\n')
                pageNum += 1
            if hexaEids:
                msg += self.writeElementTransient('CHEXA',
                                                  8, hexaEids, dt, header, hexaMsg, f)
                msg.append(pageStamp + str(pageNum) + '\n')
                pageNum += 1
            if pentaEids:
                msg += self.writeElementTransient('CPENTA',
                                                  6, pentaEids, dt, header, pentaMsg, f)
                msg.append(pageStamp + str(pageNum) + '\n')
                pageNum += 1
        return (''.join(msg), pageNum - 1)

    def writeElement(self, eType, nNodes, eids, header, tetraMsg, f=None):
        msg = header + tetraMsg
        for eid in eids:
            #eType = self.eType[eid]

            k = self.oxx[eid].keys()
            k.remove('C')
            k.sort()
            msg.append('0  %8s           0GRID CS  %i GP\n' % (eid, nNodes))
            for nid in ['C'] + k:
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

                if nid == 'C':
                    nid = 'CENTER'
                #print "nid = |%r|" %(nid)
                A = [[oxx, txy, txz],
                     [txy, oyy, tyz],
                     [txz, tyz, ozz]]
                (Lambda, v) = eigh(A)  # a hermitian matrix is a symmetric-real matrix

                ([oxx, oyy, ozz, txy, tyz, txz, o1, o2, o3, p, ovm], isAllZeros) = writeFloats13E([oxx, oyy, ozz, txy, tyz, txz, o1, o2, o3, p, ovm])
                msg.append('0              %8s  X  %13s  XY  %13s   A  %13s  LX%5.2f%5.2f%5.2f  %13s   %-s\n' % (nid, oxx, txy, o1, v[0, 1], v[0, 2], v[0, 0], p, ovm.strip()))
                msg.append('               %8s  Y  %13s  YZ  %13s   B  %13s  LY%5.2f%5.2f%5.2f\n' % ('', oyy, tyz, o2, v[1, 1], v[1, 2], v[1, 0]))
                msg.append('               %8s  Z  %13s  ZX  %13s   C  %13s  LZ%5.2f%5.2f%5.2f\n' % ('', ozz, txz, o3, v[2, 1], v[2, 2], v[2, 0]))

        if f is not None:
            f.write(''.join(msg))
        return ''.join(msg)

    def writeElementTransient(self, eType, nNodes, eids, dt, header, tetraMsg, f=None):
        dtLine = '%14s = %12.5E\n' % (self.data_code['name'], dt)
        header[1] = dtLine
        msg = header + tetraMsg
        for eid in eids:
            #eType = self.eType[eid]

            k = self.oxx[dt][eid].keys()
            k.remove('C')
            k.sort()
            msg.append('0  %8s           0GRID CS  %i GP\n' % (eid, nNodes))
            for nid in ['C'] + k:
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

                if nid == 'C':
                    nid = 'CENTER'
                #print "nid = |%r|" %(nid)
                A = [[oxx, txy, txz],
                     [txy, oyy, tyz],
                     [txz, tyz, ozz]]
                (Lambda, v) = eigh(A)  # a hermitian matrix is a symmetric-real matrix

                ([oxx, oyy, ozz, txy, tyz, txz, o1, o2, o3, p, ovm], isAllZeros) = writeFloats13E([oxx, oyy, ozz, txy, tyz, txz, o1, o2, o3, p, ovm])
                msg.append('0              %8s  X  %13s  XY  %13s   A  %13s  LX%5.2f%5.2f%5.2f  %13s   %-s\n' % (nid, oxx, txy, o1, v[0, 1], v[0, 2], v[0, 0], p, ovm.strip()))
                msg.append('               %8s  Y  %13s  YZ  %13s   B  %13s  LY%5.2f%5.2f%5.2f\n' % ('', oyy, tyz, o2, v[1, 1], v[1, 2], v[1, 0]))
                msg.append('               %8s  Z  %13s  ZX  %13s   C  %13s  LZ%5.2f%5.2f%5.2f\n' % ('', ozz, txz, o3, v[2, 1], v[2, 2], v[2, 0]))

        if f is not None:
            f.write(''.join(msg))
            msg = ['']
        return ''.join(msg)

    def get_stats(self):
        ndt, nelements, nnodes, dts = self._get_shape()
        nelements = len(self.element_data['element_id'])

        msg = self.get_data_code()
        if self.nonlinear_factor is not None:  # transient
            ntimes = len(self.oxx)  # bug
            dtstring = self.data_code['name'] + ', '
            msg.append('  type=%s ntimes=%s nelements=%s\n'
                       % (self.__class__.__name__, ntimes, nelements))
        else:
            dtstring = ''
            msg.append('  type=%s nelements=%s\n' % (self.__class__.__name__,
                                                     nelements))

        if self.isVonMises():
            ovmstr = 'ovm'
        else:
            ovmstr = 'max_shear'

        etypes = self.element_data['element_type']
        msg.append('  element_data: index        : element_id\n')
        msg.append('                results      : element_type, cid\n')
        msg.append('  data:         index        : %selement_id, node_id\n' % dtstring)
        msg.append('                results      : oxx, oyy, ozz, txy, tyz, txz, o1, o2, o3, %s\n' % ovmstr)
        msg.append('                element_types: %s' %(', '.join(set(etypes))))

        #print("self.element_data\n", self.element_data.to_string())
        #print("self.data\n", self.data.to_string())
        msg.append('\n')
        return msg

    def __repr__(self):
        return self.get_stats()


class SolidStrainObject(StrainObject):
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
    def __init__(self, data_code, is_sort1, isubcase, dt, read_mode):
        StrainObject.__init__(self, data_code, isubcase, read_mode)
        self.eType = {}

        if read_mode == 0:
            return
        self.code = [self.format_code, self.sort_code, self.s_code]

        #self.aCos = {}
        #self.bCos = {}
        #self.cCos = {}
        #self.pressure = {}

        self.nonlinear_factor = dt
        #if is_sort1:
            #if dt is not None:
                #self.add = self.add_sort1
                #self.add_new_eid = self.add_new_eid_sort1
        #else:
            #assert dt is not None
            #self.add = self.addSort2
            #self.add_new_eid = self.add_new_eid_sort2

    def add_f06_data(self, data, transient):
        if transient is None:
            if not hasattr(self, 'data'):
                self.data = []
            self.data += data
        else:
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
        eMap = {'CTETRA': 5, 'CPENTA': 7, 'CHEXA': 9, 'HEXA': 9,
                'PENTA': 7, 'TETRA': 5, }   # +1 for the centroid
        if self.nonlinear_factor is None:
            pack = []
            i = 0
            n = 0
            while n < len(self.data):
                line = self.data[n]
                #print(n,line)

                eType = line[0]
                eid = int(line[1])
                #print("eType = ",eType)
                nNodes = eMap[eType]
                #nodeID = None
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
                    #print self.data[n]
                    (blank, nodeID, x, exx, xy, exy, a, e1, lx, d1, d2, d3,
                        pressure, evmShear) = self.data[n]
                    (blank, blank, y, eyy, yz, eyz, b, e2, ly, d1, d2, d3,
                        blank, blank) = self.data[n + 1]
                    (blank, blank, z, ezz, zx, exz, c, e3, lz, d1, d2, d3,
                        blank, blank) = self.data[n + 2]
                    if    nodeID.strip() == 'CENTER':
                        nodeID = 'C'
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
            (nodeID, gridType, t1, t2, t3, r1, r2, r3) = line
            self.gridTypes[dt][nodeID] = array([t1, t2, t3])
            self.translations[dt][nodeID] = array([t1, t2, t3])
            self.rotations[dt][nodeID] = array([r1, r2, r3])
        del self.data

    def _get_shape(self):
        ndt = len(self.shape)
        dts = self.shape.keys()
        shape0 = dts[0]
        nnodes = self.shape[shape0]
        print("ndt=%s nnodes=%s dts=%s" % (ndt, nnodes, str(dts)))
        return ndt, nnodes, dts

    def _increase_size(self, dt, nelements, nnodes):
        #self.shape += 1
        if dt in self.shape:  # default dictionary
            self.shape[dt][0] += nelements
            self.shape[dt][1] += nnodes
        else:
            self.shape[dt] = [nnodes]

    def _preallocate(self, dt, nelements_size, nnodes_size):
        if self.shape is None:
            self._inode_start += nnodes
            self._inode_end += nnodes
        else:
            ndt, nelements, nnodes, dts = self._get_shape()
            #print("ndt=%s nnodes=%s dts=%s" % (ndt, nnodes, str(dts)))

            data = {}
            columns = []
            if dts[0] is not None:
                if isinstance(dt, int):
                    data['dt'] = pd.Series(zeros((ndt * nnodes), dtype='int32'))
                else:
                    data['dt'] = pd.Series(zeros((ndt * nnodes), dtype='float32'))
                columns.append('dt')

            data['element_id'] = pd.Series(zeros((ndt * nnodes), dtype='int32'))
            data['node_id'] = pd.Series(zeros((ndt * nnodes), dtype='int32'))

            columns.append('element_id')
            columns.append('node_id')

            #data['grid_type'] = pd.Series(zeros(ndt), dtype='int32'))
            #data['grid_type_str'] = pd.Series(zeros(nnodes), dtype='str'))
            data['exx'] = pd.Series(zeros((ndt * nnodes_size), dtype='float32'))
            data['eyy'] = pd.Series(zeros((ndt * nnodes_size), dtype='float32'))
            data['ezz'] = pd.Series(zeros((ndt * nnodes_size), dtype='float32'))
            data['exy'] = pd.Series(zeros((ndt * nnodes_size), dtype='float32'))
            data['exz'] = pd.Series(zeros((ndt * nnodes_size), dtype='float32'))
            data['eyz'] = pd.Series(zeros((ndt * nnodes_size), dtype='float32'))
            data['e1'] = pd.Series(zeros((ndt * nnodes_size), dtype='float32'))
            data['e2'] = pd.Series(zeros((ndt * nnodes_size), dtype='float32'))
            data['e3'] = pd.Series(zeros((ndt * nnodes_size), dtype='float32'))

            if self.isVonMises():
                data['evm'] = pd.Series(zeros((ndt * nnodes_size), dtype='float32'))
            else:
                data['max_shear'] = pd.Series(zeros((ndt * nnodes_size), dtype='float32'))
            #pressure
            #aCos
            #bCos
            #cCos
            columns += ['exx', 'eyy', 'ezz', 'exy', 'exz', 'e1', 'e2', 'e3']

            self.data = pd.DataFrame(data, columns=columns)
            self.element_data = pd.DataFrame(element_data, columns=['cid'])
            self._inode_start = 0
            self._inode_end = nnodes
        return self._inode_start, self._inode_end

    def _finalize(self):
        ndt, nelements, nnodes, dts = self._get_shape()

        #grid_type_str = []
        #for grid_type in self.grid_type:
            #grid_type_str.append('C' if grid_type==0 else grid_type)
        #self.grid_type_str = pd.Series(grid_type_str, dtype='str')

        if dts[0] is not None:
            self.data = self.data.set_index(['dt', 'element_id', 'node_id'])
        else:
            self.data = self.data.set_index(['element_id', 'node_id'])
        self.element_data = self.element_data.set_index(['element_id'])
        #print "final\n", self.data
        del self._inode_start
        del self._inode_end

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
        tetraEids = []
        hexaEids = []
        pentaEids = []
        for eid, eType in sorted(self.eType.iteritems()):
            if eType == 'CTETRA' or eType == 'TETRA':
                tetraEids.append(eid)
            elif eType == 'CPENTA' or eType == 'PENTA':
                pentaEids.append(eid)
            elif eType == 'CHEXA' or eType == 'HEXA':
                hexaEids.append(eid)
            else:
                raise NotImplementedError('eType=|%r|' % (eType))

        return (tetraMsg, pentaMsg, hexaMsg, tetraEids, hexaEids, pentaEids)

    def write_f06(self, header, pageStamp, pageNum=1, f=None, is_mag_phase=False):
        if self.nonlinear_factor is not None:
            return self._write_f06_transient(header, pageStamp, pageNum, f)
        msg = []

        (tetraMsg, pentaMsg, hexaMsg, tetraEids, hexaEids,
            pentaEids) = self.getF06_Header()
        #nNodes = {'CTETRA':4,'CPENTA':6,'CHEXA':8,'HEXA':8,'PENTA':6,'TETRA':4,}
        if tetraEids:
            msg += self.writeElement(
                'CTETRA', 4, tetraEids, header, tetraMsg, f)
            msg.append(pageStamp + str(pageNum) + '\n')
            pageNum += 1
        if hexaEids:
            msg += self.writeElement('CHEXA', 8, hexaEids, header, hexaMsg, f)
            msg.append(pageStamp + str(pageNum) + '\n')
            pageNum += 1
        if pentaEids:
            msg += self.writeElement(
                'CPENTA', 6, pentaEids, header, pentaMsg, f)
            msg.append(pageStamp + str(pageNum) + '\n')
            pageNum += 1

        return (''.join(msg), pageNum - 1)

    def _write_f06_transient(self, header, pageStamp, pageNum=1, f=None, is_mag_phase=False):
        msg = []
        (tetraMsg, pentaMsg, hexaMsg, tetraEids, hexaEids,
            pentaEids) = self.getF06_Header()
        dts = self.exx.keys()
        for dt in dts:
            if tetraEids:
                msg += self.writeElementTransient('CTETRA',
                                                  4, tetraEids, dt, header, tetraMsg, f)
                msg.append(pageStamp + str(pageNum) + '\n')
                pageNum += 1
            if hexaEids:
                msg += self.writeElementTransient('CHEXA',
                                                  8, hexaEids, dt, header, hexaMsg, f)
                msg.append(pageStamp + str(pageNum) + '\n')
                pageNum += 1
            if pentaEids:
                msg += self.writeElementTransient('CPENTA',
                                                  6, pentaEids, dt, header, pentaMsg, f)
                msg.append(pageStamp + str(pageNum) + '\n')
                pageNum += 1

        return (''.join(msg), pageNum - 1)

    def writeElement(self, eType, nNodes, eids, header, tetraMsg, f=None):
        msg = header + tetraMsg
        for eid in eids:
            #eType = self.eType[eid]

            k = self.exx[eid].keys()
            k.remove('C')
            k.sort()
            msg.append('0  %8s           0GRID CS  %i GP\n' % (eid, nNodes))
            for nid in ['C'] + k:
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

                if nid == 'C':
                    nid = 'CENTER'
                #print "nid = |%r|" %(nid)
                A = [[exx, exy, exz],
                     [exy, eyy, eyz],
                     [exz, eyz, ezz]]
                (Lambda, v) = eigh(A)  # a hermitian matrix is a symmetric-real matrix

                ([exx, eyy, ezz, exy, eyz, exz, e1, e2, e3, p, evm], isAllZeros) = writeFloats13E([exx, eyy, ezz, exy, eyz, exz, e1, e2, e3, p, evm])
                msg.append('0              %8s  X  %13s  XY  %13s   A  %13s  LX%5.2f%5.2f%5.2f  %13s   %-s\n' % (nid, exx, exy, e1, v[0, 1], v[0, 2], v[0, 0], p, evm.strip()))
                msg.append('               %8s  Y  %13s  YZ  %13s   B  %13s  LY%5.2f%5.2f%5.2f\n' % ('', eyy, eyz, e2, v[1, 1], v[1, 2], v[1, 0]))
                msg.append('               %8s  Z  %13s  ZX  %13s   C  %13s  LZ%5.2f%5.2f%5.2f\n' % ('', ezz, exz, e3, v[2, 1], v[2, 2], v[2, 0]))

        if f is not None:
            f.write(''.join(msg))
            msg = ['']
        return ''.join(msg)

    def writeElementTransient(self, eType, nNodes, eids, dt, header, tetraMsg, pageStamp=1, f=None):
        dtLine = '%14s = %12.5E\n' % (self.data_code['name'], dt)
        header[1] = dtLine
        msg = header + tetraMsg
        for eid in eids:
            #eType = self.eType[eid]

            k = self.exx[dt][eid].keys()
            k.remove('C')
            k.sort()
            msgA = '0  %8s           0GRID CS  %i GP\n' % (eid, nNodes)
            for nid in ['C'] + k:
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

                if nid == 'C':
                    nid = 'CENTER'
                #print "nid = |%r|" %(nid)
                A = [[exx, exy, exz],
                     [exy, eyy, eyz],
                     [exz, eyz, ezz]]
                (Lambda, v) = eigh(A)  # a hermitian matrix is a symmetric-real matrix

                ([exx, eyy, ezz, exy, eyz, exz, e1, e2, e3, p, evm], isAllZeros) = writeFloats13E([exx, eyy, ezz, exy, eyz, exz, e1, e2, e3, p, evm])
                msgA += '0              %8s  X  %13s  XY  %13s   A  %13s  LX%5.2f%5.2f%5.2f  %13s   %-s\n' % (nid, exx, exy, e1, v[0, 1], v[0, 2], v[0, 0], p, evm.strip())
                msgA += '               %8s  Y  %13s  YZ  %13s   B  %13s  LY%5.2f%5.2f%5.2f\n' % ('', eyy, eyz, e2, v[1, 1], v[1, 2], v[1, 0])
                msgA += '               %8s  Z  %13s  ZX  %13s   C  %13s  LZ%5.2f%5.2f%5.2f\n' % ('', ezz, exz, e3, v[2, 1], v[2, 2], v[2, 0])

        if f is not None:
            f.write(''.join(msg))
            msg = ['']
        return ''.join(msg)

    def get_stats():
        ndt, nelements, nnodes, dts = self._get_shape()
        nelements = len(self.element_data['element_id'])

        msg = self.get_data_code()
        if self.nonlinear_factor is not None:  # transient
            ntimes = len(self.exx)  # bug
            dtstring = self.data_code['name'] + ', '
            msg.append('  type=%s ntimes=%s nelements=%s\n'
                       % (self.__class__.__name__, ntimes, nelements))
        else:
            dtstring = ''
            msg.append('  type=%s nelements=%s\n' % (self.__class__.__name__,
                                                     nelements))

        if self.isVonMises():
            ovmstr = 'evm'
        else:
            ovmstr = 'max_shear'
        msg.append('  element_data: index  : element_id\n')
        msg.append('                results: element_type, cid\n')
        msg.append('  data:         index  : %selement_id, node_id\n' % dtstring)
        msg.append('                results: exx, eyy, ezz, exy, eyz, exz, '
                   'e1, e2, e3, %s\n' % ovmstr)
        msg.append(', '.join(set(self.eType.values())))
        msg.append('\n')
        return msg

    def __repr__(self):
        return self.get_stats()