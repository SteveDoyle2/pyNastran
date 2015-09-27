# pylint: disable=C0301,W0613,C0103,R0913,R0914,R0904,C0111,R0201,R0902
from __future__ import (nested_scopes, generators, division, absolute_import,
                        print_function, unicode_literals)
from six import iteritems
from six.moves import zip, range
from itertools import count
from struct import Struct, pack

from numpy import sqrt, zeros, where, searchsorted
from numpy.linalg import eigh

from pyNastran.op2.tables.oes_stressStrain.real.oes_objects import StressObject, StrainObject, OES_Object
from pyNastran.f06.f06_formatting import writeFloats13E, _eigenvalue_header


class RealSolidArray(OES_Object):
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

    def is_real(self):
        return True

    def is_complex(self):
        return False

    def get_headers(self):
        raise NotImplementedError()

    def _reset_indices(self):
        self.itotal = 0
        self.ielement = 0

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

        # TODO: could be more efficient by using nelements for cid
        self.element_node = zeros((self.ntotal, 2), dtype='int32')
        self.element_cid = zeros((self.nelements, 2), dtype='int32')

        #if self.element_name == 'CTETRA':
            #nnodes = 4
        #elif self.element_name == 'CPENTA':
            #nnodes = 6
        #elif self.element_name == 'CHEXA':
            #nnodes = 8
        #self.element_node = zeros((self.ntotal, nnodes, 2), 'int32')

        #[oxx, oyy, ozz, txy, tyz, txz, o1, o2, o3, ovmShear]
        self.data = zeros((self.ntimes, self.ntotal, 10), 'float32')
        self.nnodes = self.element_node.shape[0] // self.nelements
        #self.data = zeros((self.ntimes, self.nelements, nnodes+1, 10), 'float32')

    def add_eid_sort1(self, eType, cid, dt, eid, node_id, oxx, oyy, ozz, txy, tyz, txz, o1, o2, o3, aCos, bCos, cCos, pressure, ovm):
        assert cid >= -1, cid
        assert eid >= 0, eid

        #print "dt=%s eid=%s eType=%s" %(dt,eid,eType)
        self._times[self.itime] = dt

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

        #print('eid=%i node_id=%i exx=%s' % (eid, node_id, str(oxx)))
        self.element_node[self.itotal, :] = [eid, node_id]
        #self.element_node[self.ielement-1, self.inode-1, :] = [eid, node_id]
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
        try:
            nnodes_per_element = self.element_node.shape[0] // nelements
        except ZeroDivisionError:
            nnodes_per_element = '???'
        nnodes = self.element_node.shape[0]

        msg = []

        if self.nonlinear_factor is not None:  # transient
            msg.append('  type=%s ntimes=%i nelements=%i nnodes=%i\n  nnodes_per_element=%s (including centroid)\n'
                       % (self.__class__.__name__, ntimes, nelements, nnodes, nnodes_per_element))
            ntimes_word = 'ntimes'
        else:
            msg.append('  type=%s nelements=%i nnodes=%i\n  nodes_per_element=%i (including centroid)\n'
                       % (self.__class__.__name__, nelements, nnodes, nnodes_per_element))
            ntimes_word = 1
        msg.append('  eType, cid\n')
        headers = self.get_headers()
        n = len(headers)
        msg.append('  data: [%s, nnodes, %i] where %i=[%s]\n' % (ntimes_word, n, n, str(', '.join(headers))))
        msg.append('  data.shape = %s\n' % str(self.data.shape).replace('L', ''))
        msg.append('  element name: %s\n  ' % self.element_name)
        msg += self.get_data_code()
        return msg


    def get_element_index(self, eids):
        # elements are always sorted; nodes are not
        itot = searchsorted(eids, self.element_node[:, 0])  #[0]
        return itot

    def eid_to_element_node_index(self, eids):
        #ind = ravel([searchsorted(self.element_node[:, 0] == eid) for eid in eids])
        ind = searchsorted(eids, self.element_node[:, 0])
        #ind = ind.reshape(ind.size)
        #ind.sort()
        return ind

    def write_f06(self, header, page_stamp, page_num=1, f=None,
                  is_mag_phase=False, is_sort1=True):
        nnodes, msg_temp = _get_f06_header_nnodes(self, is_mag_phase)

        # write the f06
        ntimes = self.data.shape[0]

        eids2 = self.element_node[:, 0]
        nodes = self.element_node[:, 1]

        eids3 = self.element_cid[:, 0]
        cids3 = self.element_cid[:, 1]

        for itime in range(ntimes):
            dt = self._times[itime]
            header = _eigenvalue_header(self, header, itime, ntimes, dt)
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
            for i, deid, node_id, doxx, doyy, dozz, dtxy, dtyz, dtxz, do1, do2, do3, dp, dovm in zip(
                count(), eids2, nodes, oxx, oyy, ozz, txy, tyz, txz, o1, o2, o3, p, ovm):

                j = where(eids3 == deid)[0]
                cid = cids3[j]
                A = [[doxx, dtxy, dtxz],
                     [dtxy, doyy, dtyz],
                     [dtxz, dtyz, dozz]]
                (Lambda, v) = eigh(A)  # a hermitian matrix is a symmetric-real matrix

                ([oxxi, oyyi, ozzi, txyi, tyzi, txzi, o1i, o2i, o3i, pi, ovmi],
                 is_all_zeros) = writeFloats13E([doxx, doyy, dozz, dtxy, dtyz, dtxz,
                                                 do1, do2, do3, dp, dovm])

                if i % cnnodes == 0:
                    f.write('0  %8s    %8iGRID CS  %i GP\n' % (deid, cid, nnodes))
                    f.write('0              %8s  X  %-13s  XY  %-13s   A  %-13s  LX%5.2f%5.2f%5.2f  %-13s   %s\n'
                            '               %8s  Y  %-13s  YZ  %-13s   B  %-13s  LY%5.2f%5.2f%5.2f\n'
                            '               %8s  Z  %-13s  ZX  %-13s   C  %-13s  LZ%5.2f%5.2f%5.2f\n'
                            % ('CENTER', oxxi, txyi, o1i, v[0, 1], v[0, 2], v[0, 0], pi, ovmi,
                               '', oyyi, tyzi, o2i, v[1, 1], v[1, 2], v[1, 0],
                               '', ozzi, txzi, o3i, v[2, 1], v[2, 2], v[2, 0]))
                else:
                    f.write('0              %8s  X  %-13s  XY  %-13s   A  %-13s  LX%5.2f%5.2f%5.2f  %-13s   %s\n'
                            '               %8s  Y  %-13s  YZ  %-13s   B  %-13s  LY%5.2f%5.2f%5.2f\n'
                            '               %8s  Z  %-13s  ZX  %-13s   C  %-13s  LZ%5.2f%5.2f%5.2f\n'
                            % (node_id, oxxi, txyi, o1i, v[0, 1], v[0, 2], v[0, 0], pi, ovmi,
                               '', oyyi, tyzi, o2i, v[1, 1], v[1, 2], v[1, 0],
                               '', ozzi, txzi, o3i, v[2, 1], v[2, 2], v[2, 0]))
                i += 1
            f.write(page_stamp % page_num)
            page_num += 1
        return page_num - 1

    def _write_table_3(self, f, fascii, itable=-3, itime=0):
        import inspect
        frame = inspect.currentframe()
        call_frame = inspect.getouterframes(frame, 2)
        fascii.write('%s.write_table_3: %s\n' % (self.__class__.__name__, call_frame[1][3]))

        f.write(pack('12i', *[4, itable, 4,
                             4, 1, 4,
                             4, 0, 4,
                             4, 146, 4,
                             ]))
        approach_code = self.approach_code
        table_code = self.table_code
        isubcase = self.isubcase
        element_type = self.element_type
        #[
            #'aCode', 'tCode', 'element_type', 'isubcase',
            #'???', '???', '???', 'load_set'
            #'format_code', 'num_wide', 's_code', '???',
            #'???', '???', '???', '???',
            #'???', '???', '???', '???',
            #'???', '???', '???', '???',
            #'???', 'Title', 'subtitle', 'label']
        #random_code = self.random_code
        format_code = 1
        s_code = self.s_code
        num_wide = self.num_wide
        acoustic_flag = 0
        thermal = 0
        Title = '%-128s' % self.Title
        subtitle = '%-128s' % self.subtitle
        label = '%-128s' % self.label
        ftable3 = '50i 128s 128s 128s'
        oCode = 0
        if self.analysis_code == 1:
            lsdvmn = self.lsdvmn
        else:
            raise NotImplementedError(self.analysis_code)

        table3 = [
            approach_code, table_code, element_type, isubcase, lsdvmn,
            0, 0, self.load_set, format_code, num_wide,
            s_code, acoustic_flag, 0, 0, 0,
            0, 0, 0, 0, 0,
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            0, thermal, thermal, 0,
            Title.encode('ascii'), subtitle.encode('ascii'), label.encode('ascii'),
        ]

        n = 0
        for v in table3:
            if isinstance(v, int) or isinstance(v, float):
                n += 4
            else:
                n += len(v)
        assert n == 584, n
        data = [584] + table3 + [584]
        fmt = 'i' + ftable3 + 'i'
        print(fmt)
        print(data)
        #f.write(pack(fascii, '%s header 3c' % self.table_name, fmt, data))
        fascii.write('%s header 3c = %s\n' % (self.table_name, data))
        f.write(pack(fmt, *data))

    def write_op2(self, f, fascii, itable, date, is_mag_phase=False):
        if self.nnodes != 9:
            return
        import inspect
        frame = inspect.currentframe()
        call_frame = inspect.getouterframes(frame, 2)
        fascii.write('%s.write_op2: %s\n' % (self.__class__.__name__, call_frame[1][3]))

        if itable == -1:
            self._write_table_header(f, fascii, date)
            itable = -3

        if isinstance(self.nonlinear_factor, float):
            op2_format = '%sif' % (7 * self.ntimes)
            raise NotImplementedError()
        else:
            op2_format = 'i21f'
        s = Struct(op2_format)
        nnodes_expected = self.nnodes

        eids2 = self.element_node[:, 0]
        nodes = self.element_node[:, 1]

        eids3 = self.element_cid[:, 0]
        cids3 = self.element_cid[:, 1]

        # table 4 info
        #ntimes = self.data.shape[0]
        nnodes = self.data.shape[1]
        nelements = len(eids2)

        # 21 = 1 node, 3 principal, 6 components, 9 vectors, 2 p/ovm
        #ntotal = ((nnodes * 21) + 1) + (nelements * 4)
        ntotal = 4 + 21 * nnodes_expected

        #print('shape = %s' % str(self.data.shape))
        assert nnodes > 1, nnodes
        assert self.ntimes == 1, self.ntimes

        device_code = self.device_code
        fascii.write('  ntimes = %s\n' % self.ntimes)

        fmt = '%2i %6f'
        #print('ntotal=%s' % (ntotal))
        #assert ntotal == 193, ntotal

        struct1 = Struct(self._endian + b'ii4si')
        struct2 = Struct(self._endian + b'i20f')

        cen = b'GRID'
        for itime in range(self.ntimes):
            self._write_table_3(f, fascii, itable, itime)

            # record 4
            header = [4, -4, 4,
                      4, 1, 4,
                      4, 0, 4,
                      4, ntotal, 4,
                      4 * ntotal]
            f.write(pack('%ii' % len(header), *header))
            fascii.write('r4 [4, 0, 4]\n')
            fascii.write('r4 [4, %s, 4]\n' % (itable - 1))
            fascii.write('r4 [4, %i, 4]\n' % (4 * ntotal))

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

            print('eids3', eids3)
            cnnodes = nnodes_expected + 1
            for i, deid, node_id, doxx, doyy, dozz, dtxy, dtyz, dtxz, do1, do2, do3, dp, dovm in zip(
                count(), eids2, nodes, oxx, oyy, ozz, txy, tyz, txz, o1, o2, o3, p, ovm):
                print('  eid =', deid, node_id)

                j = where(eids3 == deid)[0]
                assert len(j) > 0, j
                cid = cids3[j][0]
                A = [[doxx, dtxy, dtxz],
                     [dtxy, doyy, dtyz],
                     [dtxz, dtyz, dozz]]
                (Lambda, v) = eigh(A)  # a hermitian matrix is a symmetric-real matrix

                #node_id, oxxi, txyi, o1i, v[0, 1], v[0, 2], v[0, 0], pi, ovmi,
                    #'', oyyi, tyzi, o2i, v[1, 1], v[1, 2], v[1, 0],
                    #'', ozzi, txzi, o3i, v[2, 1], v[2, 2], v[2, 0]
                    #(grid_device, sxx, sxy, s1, a1, a2, a3, pressure, svm,
                     #syy, syz, s2, b1, b2, b3,
                     #szz, sxz, s3, c1, c2, c3)

                if i % cnnodes == 0:
                    data = [deid * 10 + self.device_code, cid, cen, nnodes_expected]
                    fascii.write('  eid=%s cid=%s cen=%s nnodes = %s\n' % tuple(data))
                    f.write(
                        #struct1.pack(*data)
                        pack(b'2i 4s i', *data)
                    )
                #else:
                fascii.write('  nid=%i\n' % data[0])

                data = [node_id,
                        doxx, dtxy, do1, v[0, 1], v[0, 2], v[0, 0], dp, dovm,
                        doyy, dtyz, do2, v[1, 1], v[1, 2], v[1, 0],
                        dozz, dtxz, do3, v[2, 1], v[2, 2], v[2, 0], ]
                fascii.write('    oxx, txy, o1, v01, v02, v00, p, ovm = %s\n' % data[1:8])
                fascii.write('    oyy, tyz, o2, v11, v12, v10         = %s\n' % data[8:14])
                fascii.write('    ozz, txz, o3, v21, v22, v20         = %s\n' % data[14:])
                f.write(struct2.pack(*data))
                i += 1

            itable -= 2
            header = [4 * ntotal,]
            f.write(pack('i', *header))
            fascii.write('footer = %s' % header)
        header = [
            4, itable, 4,
            4, 1, 4,
            4, 0, 4,
        ]
        f.write(pack('%ii' % len(header), *header))


class RealSolidStressArray(RealSolidArray, StressObject):
    def __init__(self, data_code, is_sort1, isubcase, dt):
        RealSolidArray.__init__(self, data_code, is_sort1, isubcase, dt)
        StressObject.__init__(self, data_code, isubcase)

    def get_headers(self):
        if self.is_von_mises():
            von_mises = 'von_mises'
        else:
            von_mises = 'max_shear'
        headers = ['oxx', 'oyy', 'ozz', 'txy', 'tyz', 'txz', 'o1', 'o2', 'o3', von_mises]
        return headers


class RealSolidStrainArray(RealSolidArray, StrainObject):
    def __init__(self, data_code, is_sort1, isubcase, dt):
        RealSolidArray.__init__(self, data_code, is_sort1, isubcase, dt)
        StrainObject.__init__(self, data_code, isubcase)

    def get_headers(self):
        if self.is_von_mises():
            von_mises = 'von_mises'
        else:
            von_mises = 'max_shear'
        headers = ['exx', 'eyy', 'ezz', 'exy', 'eyz', 'exz', 'e1', 'e2', 'e3', von_mises]
        return headers

def _get_solid_msgs(self):
    if self.is_von_mises():
        von_mises = 'VON MISES'
    else:
        von_mises = 'MAX SHEAR'

    if self.is_stress():

        base_msg = [
            '0                CORNER        ------CENTER AND CORNER POINT STRESSES---------       DIR.  COSINES       MEAN                   \n',
            '  ELEMENT-ID    GRID-ID        NORMAL              SHEAR             PRINCIPAL       -A-  -B-  -C-     PRESSURE       %s \n' % von_mises]
        tetra_msg = ['                   S T R E S S E S   I N    T E T R A H E D R O N   S O L I D   E L E M E N T S   ( C T E T R A )\n', ]
        penta_msg = ['                    S T R E S S E S   I N   P E N T A H E D R O N   S O L I D   E L E M E N T S   ( P E N T A )\n', ]
        hexa_msg = ['                      S T R E S S E S   I N   H E X A H E D R O N   S O L I D   E L E M E N T S   ( H E X A )\n', ]
    else:
        base_msg = [
            '0                CORNER        ------CENTER AND CORNER POINT  STRAINS---------       DIR.  COSINES       MEAN                   \n',
            '  ELEMENT-ID    GRID-ID        NORMAL              SHEAR             PRINCIPAL       -A-  -B-  -C-     PRESSURE       %s \n' % von_mises]
        tetra_msg = ['                     S T R A I N S   I N    T E T R A H E D R O N   S O L I D   E L E M E N T S   ( C T E T R A )\n', ]
        penta_msg = ['                      S T R A I N S   I N   P E N T A H E D R O N   S O L I D   E L E M E N T S   ( P E N T A )\n', ]
        hexa_msg = ['                        S T R A I N S   I N   H E X A H E D R O N   S O L I D   E L E M E N T S   ( H E X A )\n', ]
    tetra_msg += base_msg
    penta_msg += base_msg
    hexa_msg += base_msg
    return tetra_msg, penta_msg, hexa_msg

def _get_f06_header_nnodes(self, is_mag_phase=True):
    tetra_msg, penta_msg, hexa_msg = _get_solid_msgs(self)
    if self.element_type == 39: # CTETRA
        msg = tetra_msg
        nnodes = 4
    elif self.element_type == 67: # CHEXA
        msg = hexa_msg
        nnodes = 8
    elif self.element_type == 68: # CPENTA
        msg = penta_msg
        nnodes = 6
    else:
        raise NotImplementedError('element_name=%s self.element_type=%s' % (self.element_name, self.element_type))
    return nnodes, msg

class RealSolidStress(StressObject):
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
    def __init__(self, data_code, is_sort1, isubcase, dt):
        StressObject.__init__(self, data_code, isubcase)

        self.eType = {}
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
            assert dt is not None, dt
            raise NotImplementedError('SORT2')

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
        msg.append('  scode=%s stress_bits=%s is_stress=%s dist=%s vm=%s\n' % (
            self.s_code, self.stress_bits, self.is_stress(), self.is_fiber_distance(), self.is_von_mises()))
        msg.append('  %s\n  ' % self.element_name)
        msg += self.get_data_code()
        return msg

    def add_f06_data(self, _f06_data, transient):
        if transient is None:
            if self._f06_data is None:
                self._f06_data = []
            self._f06_data += _f06_data
        else:
            dt = transient[1]
            if self._f06_data is None:
                self._f06_data = {}
            if dt not in self._f06_data:
                self._f06_data[dt] = []
            self._f06_data[dt] += _f06_data
            assert 'name' in self.data_code, self.data_code

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
        # +1 for the centroid
        if self.element_type == 39: # CTETRA
            nnodes = 5
        elif self.element_type == 67:  # CHEXA
            nnodes = 9
        elif self.element_type == 68:  # PENTA
            nnodes = 7
        else:
            raise RuntimeError('element_name=%s element_type=%s' % (self.element_name, self.element_type))
        if self.nonlinear_factor is None:
            ipack = []
            i = 0
            n = 0
            while n < len(self._f06_data):
                line = self._f06_data[n]

                etype = line[0]
                eid = int(line[1])
                self.eType[eid] = etype
                self.cid[eid] = 0  ## TODO: set this properly
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
                for j in range(nnodes):
                    (blank, node_id, x, oxx, xy, txy, a, o1, lx, d1, d2, d3, pressure, ovmShear) = self._f06_data[n]
                    (blank, blank, y, oyy, yz, tyz, b, o2, ly, d1, d2, d3, blank, blank) = self._f06_data[n+1]
                    (blank, blank, z, ozz, xz, txz, c, o3, lz, d1, d2, d3, blank, blank) = self._f06_data[n+2]
                    if node_id.strip() == 'CENTER':
                        node_id = 0
                    else:
                        node_id = int(node_id)

                    self.oxx[eid][node_id] = float(oxx)
                    self.oyy[eid][node_id] = float(oyy)
                    self.ozz[eid][node_id] = float(ozz)
                    self.txy[eid][node_id] = float(txy)
                    self.tyz[eid][node_id] = float(tyz)
                    self.txz[eid][node_id] = float(txz)
                    self.o1[eid][node_id] = float(o1)
                    self.o2[eid][node_id] = float(o2)
                    self.o3[eid][node_id] = float(o3)
                    self.ovmShear[eid][node_id] = float(ovmShear)
                    n += 3
            return

        for dt, f06_data in sorted(iteritems(self._f06_data)):
            n = 0
            if dt not in self.oxx:
                self.oxx[dt] = {}
                self.oyy[dt] = {}
                self.ozz[dt] = {}
                self.txy[dt] = {}
                self.tyz[dt] = {}
                self.txz[dt] = {}
                self.o1[dt] = {}
                self.o2[dt] = {}
                self.o3[dt] = {}
                self.ovmShear[dt] = {}

            while n < len(f06_data):
                line = f06_data[n]

                etype = line[0]
                eid = int(line[1])
                #nnodes = eMap[etype]
                self.eType[eid] = etype
                self.oxx[dt][eid] = {}
                self.oyy[dt][eid] = {}
                self.ozz[dt][eid] = {}
                self.txy[dt][eid] = {}
                self.tyz[dt][eid] = {}
                self.txz[dt][eid] = {}
                self.o1[dt][eid] = {}
                self.o2[dt][eid] = {}
                self.o3[dt][eid] = {}
                self.ovmShear[dt][eid] = {}
                n += 1
                for j in range(nnodes):
                    try:
                        (blank, node_id, x, oxx, txy, txy, a, o1, lx, d1, d2, d3, pressure, ovmShear) = f06_data[n]
                    except:
                        print('stress; error reading %s; nnodes=%s' % (self.element_name, nnodes))
                        raise
                    #(blank, node_id, x, oxx, xy, txy, a, o1, lx, d1, d2, d3, pressure, ovmShear) = f06_data[n]
                    (blank, blank, y, oyy, yz, tyz, b, o2, ly, d1, d2, d3, blank, blank) = f06_data[n + 1]
                    (blank, blank, z, ozz, zx, txz, c, o3, lz, d1, d2, d3, blank, blank) = f06_data[n + 2]
                    if node_id.strip() == 'CENTER':
                        node_id = 'CENTER'
                    else:
                        node_id = int(node_id)

                    self.oxx[dt][eid][node_id] = float(oxx)
                    self.oyy[dt][eid][node_id] = float(oyy)
                    self.ozz[dt][eid][node_id] = float(ozz)
                    self.txy[dt][eid][node_id] = float(txy)
                    self.tyz[dt][eid][node_id] = float(tyz)
                    self.txz[dt][eid][node_id] = float(txz)
                    self.o1[dt][eid][node_id] = float(o1)
                    self.o2[dt][eid][node_id] = float(o2)
                    self.o3[dt][eid][node_id] = float(o3)
                    self.ovmShear[dt][eid][node_id] = float(ovmShear)
                    n += 3

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

    def add_eid(self, eType, cid, dt, eid, node_id, oxx, oyy, ozz, txy, tyz, txz, o1, o2, o3, aCos, bCos, cCos, pressure, ovm):
        #assert eid not in self.oxx
        assert cid >= -1, cid
        assert eid >= 0, eid
        assert isinstance(node_id, int), node_id
        self.eType[eid] = eType
        self.cid[eid] = cid
        self.oxx[eid] = {node_id : oxx}
        self.oyy[eid] = {node_id : oyy}
        self.ozz[eid] = {node_id : ozz}
        self.txy[eid] = {node_id : txy}
        self.tyz[eid] = {node_id : tyz}
        self.txz[eid] = {node_id : txz}
        self.o1[eid] = {node_id : o1}
        self.o2[eid] = {node_id : o2}
        self.o3[eid] = {node_id : o3}
        #self.aCos[eid] = {node_id : aCos}
        #self.bCos[eid] = {node_id : bCos}
        #self.cCos[eid] = {node_id : cCos}
        #self.pressure[eid] = {node_id : pressure}
        self.ovmShear[eid] = {node_id : ovm}
        #msg = "*eid=%s node_id=%s vm=%g" % (eid, node_id, ovm)
        #print msg

    def add_eid_sort1(self, eType, cid, dt, eid, node_id, oxx, oyy, ozz, txy, tyz, txz, o1, o2, o3, aCos, bCos, cCos, pressure, ovm):
        assert cid >= -1, cid
        assert eid >= 0, eid
        assert isinstance(node_id, int), node_id

        if dt not in self.oxx:
            self.add_new_transient(dt)
        assert eid not in self.oxx[dt], self.oxx[dt]

        self.eType[eid] = eType
        self.cid[eid] = cid
        self.oxx[dt][eid] = {node_id : oxx}
        self.oyy[dt][eid] = {node_id : oyy}
        self.ozz[dt][eid] = {node_id : ozz}
        self.txy[dt][eid] = {node_id : txy}
        self.tyz[dt][eid] = {node_id : tyz}
        self.txz[dt][eid] = {node_id : txz}

        self.o1[dt][eid] = {node_id : o1}
        self.o2[dt][eid] = {node_id : o2}
        self.o3[dt][eid] = {node_id : o3}

        #self.aCos[dt][eid] = {node_id : aCos}
        #self.bCos[dt][eid] = {node_id : bCos}
        #self.cCos[dt][eid] = {node_id : cCos}
        #self.pressure[dt][eid] = {node_id : pressure}
        self.ovmShear[dt][eid] = {node_id : ovm}

    def add_node(self, dt, eid, inode, node_id, oxx, oyy, ozz, txy, tyz, txz, o1, o2, o3, aCos, bCos, cCos, pressure, ovm):
        assert isinstance(node_id, int), node_id
        #msg = "eid=%s node_id=%s vm=%g" % (eid, node_id, ovm)
        self.oxx[eid][node_id] = oxx
        self.oyy[eid][node_id] = oyy
        self.ozz[eid][node_id] = ozz

        self.txy[eid][node_id] = txy
        self.tyz[eid][node_id] = tyz
        self.txz[eid][node_id] = txz

        self.o1[eid][node_id] = o1
        self.o2[eid][node_id] = o2
        self.o3[eid][node_id] = o3

        #self.aCos[eid][node_id] = aCos
        #self.bCos[eid][node_id] = bCos
        #self.cCos[eid][node_id] = cCos
        #self.pressure[eid][node_id] = pressure
        self.ovmShear[eid][node_id] = ovm

    def add_node_sort1(self, dt, eid, inode, node_id, oxx, oyy, ozz, txy, tyz, txz, o1, o2, o3, aCos, bCos, cCos, pressure, ovm):
        assert isinstance(node_id, int), node_id
        #self.eType[eid] = eType
        #print "eid=%s nid=%s oxx=%s" %(eid,node_id,oxx)
        self.oxx[dt][eid][node_id] = oxx
        self.oyy[dt][eid][node_id] = oyy
        self.ozz[dt][eid][node_id] = ozz

        self.txy[dt][eid][node_id] = txy
        self.tyz[dt][eid][node_id] = tyz
        self.txz[dt][eid][node_id] = txz

        self.o1[dt][eid][node_id] = o1
        self.o2[dt][eid][node_id] = o2
        self.o3[dt][eid][node_id] = o3

        #self.aCos[dt][eid][node_id] = aCos
        #self.bCos[dt][eid][node_id] = bCos
        #self.cCos[dt][eid][node_id] = cCos
        #self.pressure[dt][eid][node_id] = pressure
        self.ovmShear[dt][eid][node_id] = ovm

    def get_headers(self):
        headers = ['oxx', 'oyy', 'ozz', 'txy', 'tyz', 'txz']
        if self.is_von_mises():
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

    def directional_vectors(self, oxx, oyy, ozz, txy, tyz, txz):
        A = [[oxx, txy, txz],
             [txy, oyy, tyz],
             [txz, tyz, ozz]]
        (Lambda, v) = eigh(A)  # a hermitian matrix is a symmetric-real matrix
        return v

    def write_f06(self, header, page_stamp, page_num=1, f=None,
                  is_mag_phase=False, is_sort1=True):
        assert f is not None
        if self.nonlinear_factor is not None:
            return self._write_f06_transient(header, page_stamp, page_num, f,
                                             is_mag_phase=is_mag_phase, is_sort1=is_sort1)

        tetraMsg, pentaMsg, hexaMsg = _get_solid_msgs(self)
        eids = self.oxx.keys()
        if self.element_type == 39: # CTETRA
            nnodes = 4
            self._write_element('CTETRA4', nnodes, eids, header, tetraMsg, f)
        elif self.element_type == 68: # CPENTA
            nnodes = 6
            self._write_element('CPENTA6', nnodes, eids, header, pentaMsg, f)
        elif self.element_type == 67: # CHEXA
            nnodes = 8
            self._write_element('CHEXA8', nnodes, eids, header, hexaMsg, f)
        else:
            raise NotImplementedError('element_name=%r type=%s' % (self.element_name, self.element_type))
        f.write(page_stamp % page_num)
        return page_num

    def _write_f06_transient(self, header, page_stamp, page_num, f,
                             is_mag_phase=False, is_sort1=True):
        tetraMsg, pentaMsg, hexaMsg = _get_solid_msgs(self)

        itime = 0
        ntimes = len(self.oxx)
        if self.element_type == 39: # CTETRA
            nnodes = 4
            for dt, oxx in sorted(iteritems(self.oxx)):
                eids = sorted(oxx.keys())
                header = _eigenvalue_header(self, header, itime, ntimes, dt)
                self._write_element_transient('CTETRA4', nnodes, eids, dt, header, tetraMsg, f)
                f.write(page_stamp % page_num)
                page_num += 1
                itime += 1
        elif self.element_type == 68: # CPENTA
            nnodes = 6
            for dt, oxx in sorted(iteritems(self.oxx)):
                eids = sorted(oxx.keys())
                header = _eigenvalue_header(self, header, itime, ntimes, dt)
                self._write_element_transient('CPENTA6', nnodes, eids, dt, header, pentaMsg, f)
                f.write(page_stamp % page_num)
                page_num += 1
                itime += 1
        elif self.element_type == 67: # CHEXA
            nnodes = 8
            for dt, oxx in sorted(iteritems(self.oxx)):
                eids = sorted(oxx.keys())
                header = _eigenvalue_header(self, header, itime, ntimes, dt)
                self._write_element_transient('CHEXA8', nnodes, eids, dt, header, hexaMsg, f)
                f.write(page_stamp % page_num)
                page_num += 1
                itime += 1
        else:
            raise NotImplementedError('element_name=%r type=%s' % (self.element_name, self.element_type))
        return page_num - 1

    def _write_element(self, eType, nnodes, eids, header, tetraMsg, f):
        f.write(''.join(header + tetraMsg))
        for eid in eids:
            nids = sorted(self.oxx[eid].keys())
            cid = self.cid[eid]

            f.write('0  %8s    %8iGRID CS  %i GP\n' % (eid, cid, nnodes))
            for nid in nids:
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

                ([oxx, oyy, ozz, txy, tyz, txz, o1, o2, o3, p, ovm], is_all_zeros) = writeFloats13E([
                  oxx, oyy, ozz, txy, tyz, txz, o1, o2, o3, p, ovm])
                if nid == 0:
                    nid = 'CENTER'
                f.write('0              %8s  X  %-13s  XY  %-13s   A  %-13s  LX%5.2f%5.2f%5.2f  %-13s   %s\n'
                        '               %8s  Y  %-13s  YZ  %-13s   B  %-13s  LY%5.2f%5.2f%5.2f\n'
                        '               %8s  Z  %-13s  ZX  %-13s   C  %-13s  LZ%5.2f%5.2f%5.2f\n' % (
                            nid, oxx, txy, o1, v[0, 1], v[0, 2], v[0, 0], p, ovm,
                            '', oyy, tyz, o2, v[1, 1], v[1, 2], v[1, 0],
                            '', ozz, txz, o3, v[2, 1], v[2, 2], v[2, 0]))

    def _write_element_transient(self, eType, nnodes, eids, dt, header, tetraMsg, f):
        #dtLine = '%14s = %12.5E\n' % (self.data_code['name'], dt)
        #header[1] = dtLine
        f.write(''.join(header + tetraMsg))
        for eid in eids:
            cid = self.cid[eid]
            nids = sorted(self.oxx[dt][eid].keys())
            f.write('0  %8s    %8iGRID CS  %i GP\n' % (eid, cid, nnodes))
            for nid in nids:
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

                A = [[oxx, txy, txz],
                     [txy, oyy, tyz],
                     [txz, tyz, ozz]]
                (Lambda, v) = eigh(A)  # a hermitian matrix is a symmetric-real matrix

                ([oxx, oyy, ozz, txy, tyz, txz, o1, o2, o3, p, ovm], is_all_zeros) = writeFloats13E([
                  oxx, oyy, ozz, txy, tyz, txz, o1, o2, o3, p, ovm])
                if nid == 0:
                    nid = 'CENTER'
                f.write('0              %8s  X  %-13s  XY  %-13s   A  %-13s  LX%5.2f%5.2f%5.2f  %-13s   %s\n'
                        '               %8s  Y  %-13s  YZ  %-13s   B  %-13s  LY%5.2f%5.2f%5.2f\n'
                        '               %8s  Z  %-13s  ZX  %-13s   C  %-13s  LZ%5.2f%5.2f%5.2f\n'
                        % (nid, oxx, txy, o1, v[0, 1], v[0, 2], v[0, 0], p, ovm,
                           '', oyy, tyz, o2, v[1, 1], v[1, 2], v[1, 0],
                           '', ozz, txz, o3, v[2, 1], v[2, 2], v[2, 0]))


class RealSolidStrain(StrainObject):
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
        self._data = None
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
                self.add_node = self.add_node_sort1
                self.add_eid = self.add_eid_sort1
        else:
            assert dt is not None

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
        msg.append('  %s\n' % self.element_name)
        msg.append('  scode=%s stress_bits=%s is_stress=%s dist=%s vm=%s\n' % (
            self.s_code, self.stress_bits, self.is_stress(), self.is_fiber_distance(), self.is_von_mises()))
        msg += self.get_data_code()
        return msg

    def add_f06_data(self, _f06_data, transient):
        if transient is None:
            if self._f06_data is None:
                self._f06_data = []
            self._f06_data += _f06_data
        else:
            dt = transient[1]
            if self._f06_data is None:
                self._f06_data = {}
            if dt not in self._f06_data:
                self._f06_data[dt] = []
            self._f06_data[dt] += _f06_data
            assert 'name' in self.data_code, self.data_code

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
        # +1 for the centroid
        if self.element_type == 39: # CTETRA
            nnodes = 5
        elif self.element_type == 67:  # CHEXA
            nnodes = 9
        elif self.element_type == 68:  # PENTA
            nnodes = 7
        else:
            raise RuntimeError('element_name=%s element_type=%s' % (self.element_name, self.element_type))

        if self._f06_data is None:
            return
        if self.nonlinear_factor is None:
            pack = []
            i = 0
            n = 0
            while n < len(self._f06_data):
                line = self._f06_data[n]

                etype = line[0]
                eid = int(line[1])
                self.cid[eid] = 0  ## TODO: set this properly
                self.eType[eid] = etype
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
                for j in range(nnodes):
                    (blank, node_id, x, exx, xy, exy, a, e1, lx, d1, d2, d3, pressure, evmShear) = self._f06_data[n]
                    (blank, blank, y, eyy, yz, eyz, b, e2, ly, d1, d2, d3, blank, blank) = self._f06_data[n + 1]
                    (blank, blank, z, ezz, zx, exz, c, e3, lz, d1, d2, d3, blank, blank) = self._f06_data[n + 2]
                    if node_id.strip() == 'CENTER':
                        node_id = 0
                    else:
                        node_id = int(node_id)
                    self.exx[eid][node_id] = float(exx)
                    self.eyy[eid][node_id] = float(eyy)
                    self.ezz[eid][node_id] = float(ezz)
                    self.exy[eid][node_id] = float(exy)
                    self.eyz[eid][node_id] = float(eyz)
                    self.exz[eid][node_id] = float(exz)
                    self.e1[eid][node_id] = float(e1)
                    self.e2[eid][node_id] = float(e2)
                    self.e3[eid][node_id] = float(e3)
                    self.evmShear[eid][node_id] = float(evmShear)
                    n += 3
            return

        #(dtName, dt) = transient
        #self.data_code['name'] = dtName
        for dt, f06_data in sorted(iteritems(self._f06_data)):
            n = 0
            if dt not in self.exx:
                self.exx[dt] = {}
                self.eyy[dt] = {}
                self.ezz[dt] = {}
                self.exy[dt] = {}
                self.eyz[dt] = {}
                self.exz[dt] = {}
                self.e1[dt] = {}
                self.e2[dt] = {}
                self.e3[dt] = {}
                self.evmShear[dt] = {}

            while n < len(f06_data):
                line = f06_data[n]

                etype = line[0]
                eid = int(line[1])
                self.eType[eid] = etype
                self.exx[dt][eid] = {}
                self.eyy[dt][eid] = {}
                self.ezz[dt][eid] = {}
                self.exy[dt][eid] = {}
                self.eyz[dt][eid] = {}
                self.exz[dt][eid] = {}
                self.e1[dt][eid] = {}
                self.e2[dt][eid] = {}
                self.e3[dt][eid] = {}
                self.evmShear[dt][eid] = {}
                n += 1
                for j in range(nnodes):
                    try:
                        (blank, node_id, x, oxx, txy, txy, a, o1, lx, d1, d2, d3, pressure, ovmShear) = f06_data[n]
                    except:
                        print('strain; error reading %s; nnodes=%s\n%r' % (self.element_name, nnodes, f06_data[n]))
                        raise

                    (blank, blank, y, oyy, tyz, tyz, b, o2, ly, d1, d2, d3, blank, blank) = f06_data[n+1]
                    (blank, blank, z, ozz, tzx, txz, c, o3, lz, d1, d2, d3, blank, blank) = f06_data[n+2]
                    if j == 0:
                        assert node_id == 'CENTER', 'j=%s node_id=%s' % (j, node_id)
                        node_id = 0
                    else:
                        node_id = int(node_id)

                    self.exx[dt][eid][node_id] = float(oxx)
                    self.eyy[dt][eid][node_id] = float(oyy)
                    self.ezz[dt][eid][node_id] = float(ozz)
                    self.exy[dt][eid][node_id] = float(txy)
                    self.eyz[dt][eid][node_id] = float(tyz)
                    self.exz[dt][eid][node_id] = float(txz)
                    self.e1[dt][eid][node_id] = float(o1)
                    self.e2[dt][eid][node_id] = float(o2)
                    self.e3[dt][eid][node_id] = float(o3)
                    self.evmShear[dt][eid][node_id] = float(ovmShear)
                    n += 3
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

    def add_eid(self, eType, cid, dt, eid, node_id, exx, eyy, ezz, exy, eyz, exz, e1, e2, e3, aCos, bCos, cCos, pressure, evm):
        assert eid not in self.exx
        assert cid >= -1, cid
        assert eid >= 0, eid
        assert isinstance(node_id, int), node_id
        self.eType[eid] = eType
        self.cid[eid] = cid
        self.exx[eid] = {node_id : exx}
        self.eyy[eid] = {node_id : eyy}
        self.ezz[eid] = {node_id : ezz}
        self.exy[eid] = {node_id : exy}
        self.eyz[eid] = {node_id : eyz}
        self.exz[eid] = {node_id : exz}
        self.e1[eid] = {node_id : e1}
        self.e2[eid] = {node_id : e2}
        self.e3[eid] = {node_id : e3}
        #self.aCos[eid] = {node_id : aCos}
        #self.bCos[eid] = {node_id : bCos}
        #self.cCos[eid] = {node_id : cCos}
        #self.pressure[eid] = {node_id : pressure}
        self.evmShear[eid] = {node_id : evm}

    def add_eid_sort1(self, eType, cid, dt, eid, node_id, exx, eyy, ezz, exy, eyz, exz, e1, e2, e3, aCos, bCos, cCos, pressure, evm):
        assert cid >= -1, cid
        assert eid >= 0, eid
        if dt not in self.exx:
            self.add_new_transient(dt)

        assert eid not in self.exx[dt], self.exx[dt]
        assert isinstance(node_id, int), node_id
        self.eType[eid] = eType
        self.cid[eid] = cid
        self.exx[dt][eid] = {node_id : exx}
        self.eyy[dt][eid] = {node_id : eyy}
        self.ezz[dt][eid] = {node_id : ezz}
        self.exy[dt][eid] = {node_id : exy}
        self.eyz[dt][eid] = {node_id : eyz}
        self.exz[dt][eid] = {node_id : exz}
        self.e1[dt][eid] = {node_id : e1}
        self.e2[dt][eid] = {node_id : e2}
        self.e3[dt][eid] = {node_id : e3}
        #self.aCos[dt][eid] = {node_id : aCos}
        #self.bCos[dt][eid] = {node_id : bCos}
        #self.cCos[dt][eid] = {node_id : cCos}
        #self.pressure[dt][eid] = {node_id : pressure}
        self.evmShear[dt][eid] = {node_id : evm}

    def add_node(self, dt, eid, inode, node_id, exx, eyy, ezz, exy, eyz, exz, e1, e2, e3, aCos, bCos, cCos, pressure, evm):
        assert isinstance(node_id, int), node_id
        self.exx[eid][node_id] = exx
        self.eyy[eid][node_id] = eyy
        self.ezz[eid][node_id] = ezz

        self.exy[eid][node_id] = exy
        self.eyz[eid][node_id] = eyz
        self.exz[eid][node_id] = exz
        self.e1[eid][node_id] = e1
        self.e2[eid][node_id] = e2
        self.e3[eid][node_id] = e3

        #self.aCos[eid][node_id] = aCos
        #self.bCos[eid][node_id] = bCos
        #self.cCos[eid][node_id] = cCos
        #self.pressure[eid][node_id] = pressure
        self.evmShear[eid][node_id] = evm

    def add_node_sort1(self, dt, eid, inode, node_id, exx, eyy, ezz, exy, eyz, exz, e1, e2, e3, aCos, bCos, cCos, pressure, evm):
        assert isinstance(node_id, int), node_id
        self.exx[dt][eid][node_id] = exx
        self.eyy[dt][eid][node_id] = eyy
        self.ezz[dt][eid][node_id] = ezz

        self.exy[dt][eid][node_id] = exy
        self.eyz[dt][eid][node_id] = eyz
        self.exz[dt][eid][node_id] = exz

        self.e1[dt][eid][node_id] = e1
        self.e2[dt][eid][node_id] = e2
        self.e3[dt][eid][node_id] = e3

        #self.aCos[dt][eid][node_id] = aCos
        #self.bCos[dt][eid][node_id] = bCos
        #self.cCos[dt][eid][node_id] = cCos
        #self.pressure[dt][eid][node_id] = pressure
        self.evmShear[dt][eid][node_id] = evm

    def get_headers(self):
        headers = ['exx', 'eyy', 'ezz', 'exy', 'eyz', 'exz']
        if self.is_von_mises():
            headers.append('evm')
        else:
            headers.append('maxShear')
        return headers

    def ovm(self, o11, o22, o33, o12, o13, o23):
        """http://en.wikipedia.org/wiki/Von_Mises_yield_criterion"""
        ovm = sqrt(0.5 * ((o11 - o22) ** 2 +
                          (o22 - o33) ** 2 +
                          (o11 - o33) ** 2 +
                     6 * (o23 ** 2 + o13 ** 2 + o12 ** 2)))
        return ovm

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

    def write_f06(self, header, page_stamp, page_num=1, f=None, is_mag_phase=False, is_sort1=True):
        assert f is not None
        if self.nonlinear_factor is not None:
            return self._write_f06_transient(header, page_stamp, page_num, f,
                                             is_mag_phase=is_mag_phase, is_sort1=is_sort1)

        tetraMsg, pentaMsg, hexaMsg = _get_solid_msgs(self)
        eids = sorted(self.exx.keys())
        if self.element_type == 39: # CTETRA
            nnodes = 4
            self._write_element('CTETRA4', nnodes, eids, header, tetraMsg, f)
        elif self.element_type == 68: # CPENTA
            nnodes = 6
            self._write_element('CPENTA6', nnodes, eids, header, pentaMsg, f)
        elif self.element_type == 67: # CHEXA
            nnodes = 8
            self._write_element('CHEXA8', nnodes, eids, header, hexaMsg, f)
        else:
            raise NotImplementedError('element_name=%r type=%s' % (self.element_name, self.element_type))
        f.write(page_stamp % page_num)
        return page_num

    def _write_f06_transient(self, header, page_stamp, page_num=1, f=None,
                             is_mag_phase=False, is_sort1=True):
        tetraMsg, pentaMsg, hexaMsg = _get_solid_msgs(self)
        itime = 0
        ntimes = len(self.exx)
        if self.element_type == 39: # CTETRA
            nnodes = 4
            for dt, oxx in sorted(iteritems(self.exx)):
                eids = sorted(oxx.keys())
                header = _eigenvalue_header(self, header, itime, ntimes, dt)
                self._write_element_transient('CTETRA4', nnodes, eids, dt, header, tetraMsg, f)
                f.write(page_stamp % page_num)
                page_num += 1
                itime += 1
        elif self.element_type == 68: # CPENTA
            nnodes = 6
            for dt, oxx in sorted(iteritems(self.exx)):
                eids = sorted(oxx.keys())
                header = _eigenvalue_header(self, header, itime, ntimes, dt)
                self._write_element_transient('CPENTA6', nnodes, eids, dt, header, pentaMsg, f)
                f.write(page_stamp % page_num)
                page_num += 1
                itime += 1
        elif self.element_type == 67: # CHEXA
            nnodes = 8
            for dt, oxx in sorted(iteritems(self.exx)):
                eids = sorted(oxx.keys())
                header = _eigenvalue_header(self, header, itime, ntimes, dt)
                self._write_element_transient('CHEXA8', nnodes, eids, dt, header, hexaMsg, f)
                f.write(page_stamp % page_num)
                itime += 1
                page_num += 1
        else:
            raise NotImplementedError('element_name=%r type=%s' % (self.element_name, self.element_type))
        return page_num - 1

    def _write_element(self, eType, nnodes, eids, header, tetraMsg, f):
        f.write(''.join(header + tetraMsg))
        for eid in eids:
            nids = sorted(self.exx[eid].keys())
            cid = self.cid[eid]
            f.write('0  %8s    %8iGRID CS  %i GP\n' % (eid, cid, nnodes))
            for nid in nids:
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

                ([exx, eyy, ezz, exy, eyz, exz, e1, e2, e3, p, evm], is_all_zeros) = writeFloats13E([
                  exx, eyy, ezz, exy, eyz, exz, e1, e2, e3, p, evm])
                if nid == 0:
                    nid = 'CENTER'
                f.write('0              %8s  X  %-13s  XY  %-13s   A  %-13s  LX%5.2f%5.2f%5.2f  %-13s   %s\n'
                        '               %8s  Y  %-13s  YZ  %-13s   B  %-13s  LY%5.2f%5.2f%5.2f\n'
                        '               %8s  Z  %-13s  ZX  %-13s   C  %-13s  LZ%5.2f%5.2f%5.2f\n'
                        % (nid, exx, exy, e1, v[0, 1], v[0, 2], v[0, 0], p, evm,
                           '', eyy, eyz, e2, v[1, 1], v[1, 2], v[1, 0],
                           '', ezz, exz, e3, v[2, 1], v[2, 2], v[2, 0]), )

    def _write_element_transient(self, eType, nnodes, eids, dt, header, tetraMsg, f):
        #dtLine = '%14s = %12.5E\n' % (self.data_code['name'], dt)
        #header[1] = dtLine
        f.write(''.join(header + tetraMsg))
        for eid in eids:
            cid = self.cid[eid]
            f.write('0  %8s    %8iGRID CS  %i GP\n' % (eid, cid, nnodes))
            nids = sorted(self.exx[dt][eid].keys())
            for nid in nids:
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

                ([exx, eyy, ezz, exy, eyz, exz, e1, e2, e3, p, evm], is_all_zeros) = writeFloats13E([
                  exx, eyy, ezz, exy, eyz, exz, e1, e2, e3, p, evm])
                if nid == 0:
                    nid = 'CENTER'
                f.write('0              %8s  X  %-13s  XY  %-13s   A  %-13s  LX%5.2f%5.2f%5.2f  %-13s   %s\n'
                        '               %8s  Y  %-13s  YZ  %-13s   B  %-13s  LY%5.2f%5.2f%5.2f\n'
                        '               %8s  Z  %-13s  ZX  %-13s   C  %-13s  LZ%5.2f%5.2f%5.2f\n' % (
                            nid, exx, exy, e1, v[0, 1], v[0, 2], v[0, 0], p, evm,
                            '', eyy, eyz, e2, v[1, 1], v[1, 2], v[1, 0],
                            '', ezz, exz, e3, v[2, 1], v[2, 2], v[2, 0]))
