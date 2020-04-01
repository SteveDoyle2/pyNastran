# pylint: disable=C0301,C0103,R0913,R0914,R0904,C0111,R0201,R0902
from itertools import count
from struct import Struct, pack
from typing import List

import numpy as np
from numpy import zeros, where, searchsorted
from numpy.linalg import eigh  # type: ignore

from pyNastran.utils.numpy_utils import integer_types, float_types
from pyNastran.f06.f06_formatting import write_floats_13e, _eigenvalue_header
from pyNastran.op2.result_objects.op2_objects import get_times_dtype
from pyNastran.op2.tables.oes_stressStrain.real.oes_objects import StressObject, StrainObject, OES_Object


class RealSolidArray(OES_Object):
    def __init__(self, data_code, is_sort1, isubcase, dt):
        OES_Object.__init__(self, data_code, isubcase, apply_data_code=False)
        #self.code = [self.format_code, self.sort_code, self.s_code]

        #self.ntimes = 0  # or frequency/mode
        #self.ntotal = 0
        self.nelements = 0  # result specific

        #if is_sort1:
            ##sort1
            #self.add_node = self.add_node_sort1
            #self.add_eid = self.add_eid_sort1
        #else:
            #raise NotImplementedError('SORT2')

    @property
    def is_real(self):
        return True

    @property
    def is_complex(self):
        return False

    def get_headers(self):
        raise NotImplementedError()

    def _reset_indices(self):
        self.itotal = 0
        self.ielement = 0

    def update_data_components(self):
        # vm
        oxx = self.data[:, :, 0]
        oyy = self.data[:, :, 1]
        ozz = self.data[:, :, 2]
        txy = self.data[:, :, 3]
        tyz = self.data[:, :, 4]
        txz = self.data[:, :, 5]

        #I1 = oxx + oyy + ozz
        #txyz = txy**2 + tyz**2 + txz ** 2
        #I2 = oxx * oyy + oyy * ozz + ozz * oxx - txyz
        #I3 = oxx * oyy * ozz + 2 * txy * tyz * txz + oxx * tyz**2 - oyy * txz**2 - ozz * txy

        # (n_subarrays, nrows, ncols)

        ovm = np.sqrt((oxx - oyy)**2 + (oyy - ozz)**2 + (oxx - ozz)**2 +
                      3. * (txy**2 + tyz**2 + txz ** 2))
        self.data[:, :, 9] = ovm
        #A = [[doxx, dtxy, dtxz],
             #[dtxy, doyy, dtyz],
             #[dtxz, dtyz, dozz]]
        #(_lambda, v) = eigh(A)  # a hermitian matrix is a symmetric-real matrix

    def __iadd__(self, factor):
        """[A] += b"""
        #[oxx, oyy, ozz, txy, tyz, txz, o1, o2, o3, ovmShear]
        if isinstance(factor, float_types):
            self.data[:, :, :6] += factor
        else:
            # TODO: should support arrays
            raise TypeError('factor=%s and must be a float' % (factor))
        self.update_data_components()

    def __isub__(self, factor):
        """[A] -= b"""
        if isinstance(factor, float_types):
            self.data[:, :, :6] -= factor
        else:
            # TODO: should support arrays
            raise TypeError('factor=%s and must be a float' % (factor))
        self.update_data_components()

    def __imul__(self, factor):
        """[A] *= b"""
        assert isinstance(factor, float_types), 'factor=%s and must be a float' % (factor)
        self.data[:, :, :6] *= factor
        self.update_data_components()

    def __idiv__(self, factor):
        """[A] *= b"""
        assert isinstance(factor, float_types), 'factor=%s and must be a float' % (factor)
        self.data[:, :, :6] *= 1. / factor
        self.update_data_components()

    def build(self):
        """sizes the vectorized attributes of the RealSolidArray"""
        #print('ntimes=%s nelements=%s ntotal=%s' % (self.ntimes, self.nelements, self.ntotal))
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
        dtype, idtype, fdtype = get_times_dtype(self.nonlinear_factor, self.size)
        _times = zeros(self.ntimes, dtype=dtype)

        # TODO: could be more efficient by using nelements for cid
        element_node = zeros((self.ntotal, 2), dtype=idtype)
        element_cid = zeros((self.nelements, 2), dtype=idtype)

        #if self.element_name == 'CTETRA':
            #nnodes = 4
        #elif self.element_name == 'CPENTA':
            #nnodes = 6
        #elif self.element_name == 'CHEXA':
            #nnodes = 8
        #self.element_node = zeros((self.ntotal, nnodes, 2), 'int32')

        #[oxx, oyy, ozz, txy, tyz, txz, o1, o2, o3, ovmShear]
        data = zeros((self.ntimes, self.ntotal, 10), fdtype)
        self.nnodes = element_node.shape[0] // self.nelements
        #self.data = zeros((self.ntimes, self.nelements, nnodes+1, 10), 'float32')

        if self.load_as_h5:
            #for key, value in sorted(self.data_code.items()):
                #print(key, value)
            group = self._get_result_group()
            self._times = group.create_dataset('_times', data=_times)
            self.element_node = group.create_dataset('element_node', data=element_node)
            self.element_cid = group.create_dataset('element_cid', data=element_cid)
            self.data = group.create_dataset('data', data=data)
        else:
            self._times = _times
            self.element_node = element_node
            self.element_cid = element_cid
            self.data = data

    def build_dataframe(self):
        """creates a pandas dataframe"""
        import pandas as pd

        headers = self.get_headers()
        # TODO: cid?
        #element_node = [self.element_node[:, 0], self.element_node[:, 1]]
        if self.nonlinear_factor not in (None, np.nan):
            column_names, column_values = self._build_dataframe_transient_header()
            data_frame = self._build_pandas_transient_element_node(
                column_values, column_names,
                headers, self.element_node, self.data)
            #self.data_frame = pd.Panel(self.data, items=column_values, major_axis=element_node, minor_axis=headers).to_frame()
            #self.data_frame.columns.names = column_names
            #self.data_frame.index.names = ['ElementID', 'NodeID', 'Item']
        else:
            # Static            sxc  sxd  sxe  sxf  smax  smin    MS_tension  MS_compression
            # ElementID NodeID
            # 12        22      0.0  0.0  0.0  0.0   0.0   0.0  1.401298e-45    1.401298e-45
            #           26      0.0  0.0  0.0  0.0   0.0   0.0  1.401298e-45    1.401298e-45
            index = pd.MultiIndex.from_arrays(self.element_node.T, names=['ElementID', 'NodeID'])
            data_frame = pd.DataFrame(self.data[0], columns=headers, index=index)
            data_frame.columns.names = ['Static']
        self.data_frame = data_frame

    def add_eid_sort1(self, unused_etype, cid, dt, eid, unused_node_id,
                      oxx, oyy, ozz, txy, tyz, txz, o1, o2, o3,
                      unused_acos, unused_bcos, unused_ccos, unused_pressure, ovm):
        assert cid >= -1, cid
        assert eid >= 0, eid

        #print "dt=%s eid=%s eType=%s" %(dt,eid,eType)
        self._times[self.itime] = dt
        self.element_node[self.itotal, :] = [eid, 0]  # 0 is center

        omax_mid_min = [o1, o2, o3]
        omin = min(omax_mid_min)
        omax_mid_min.remove(omin)

        omax = max(omax_mid_min)
        omax_mid_min.remove(omax)

        omid = omax_mid_min[0]
        self.data[self.itime, self.itotal, :] = [oxx, oyy, ozz, txy, tyz, txz, omax, omid, omin, ovm]
        #self.data[self.itime, self.ielement, 0, :] = [oxx, oyy, ozz, txy, tyz, txz, o1, o2, o3, ovm]

        #print('element_cid[%i, :] = [%s, %s]' % (self.ielement, eid, cid))
        if self.ielement == self.nelements:
            self.ielement = 0
        self.element_cid[self.ielement, :] = [eid, cid]
        self.itotal += 1
        self.ielement += 1

    def __eq__(self, table):  # pragma: no cover
        assert self.is_sort1 == table.is_sort1
        self._eq_header(table)
        if not np.array_equal(self.data, table.data):
            msg = 'table_name=%r class_name=%s\n' % (self.table_name, self.__class__.__name__)
            msg += '%s\n' % str(self.code_information())
            i = 0
            for itime in range(self.ntimes):
                for ieid, eid_nid in enumerate(self.element_node):
                    (eid, nid) = eid_nid
                    t1 = self.data[itime, ieid, :]
                    t2 = table.data[itime, ieid, :]
                    (oxx1, oyy1, ozz1, txy1, tyz1, txz1, o11, o21, o31, ovm1) = t1
                    (oxx2, oyy2, ozz2, txy2, tyz2, txz2, o12, o22, o32, ovm2) = t2

                    if not np.array_equal(t1, t2):
                        msg += (
                            '(%s, %s)    (%s, %s, %s, %s, %s, %s, %s, %s, %s, %s)\n'
                            '%s      (%s, %s, %s, %s, %s, %s, %s, %s, %s, %s)\n' % (
                                eid, nid,
                                oxx1, oyy1, ozz1, txy1, tyz1, txz1, o11, o21, o31, ovm1,
                                ' ' * (len(str(eid)) + len(str(nid)) + 2),
                                oxx2, oyy2, ozz2, txy2, tyz2, txz2, o12, o22, o32, ovm2))
                        i += 1
                        if i > 10:
                            print(msg)
                            raise ValueError(msg)
                #print(msg)
                if i > 0:
                    raise ValueError(msg)
        return True

    def add_node_sort1(self, dt, eid, unused_inode, node_id,
                       oxx, oyy, ozz, txy, tyz, txz, o1, o2, o3,
                       unused_acos, unused_bcos, unused_ccos, unused_pressure, ovm):
        # skipping aCos, bCos, cCos, pressure
        omax_mid_min = [o1, o2, o3]
        omin = min(omax_mid_min)
        omax_mid_min.remove(omin)

        omax = max(omax_mid_min)
        omax_mid_min.remove(omax)

        omid = omax_mid_min[0]
        self.data[self.itime, self.itotal, :] = [oxx, oyy, ozz, txy, tyz, txz, omax, omid, omin, ovm]
        #print('data[%s, %s, :] = %s' % (self.itime, self.itotal, str(self.data[self.itime, self.itotal, :])))

        #self.data[self.itime, self.ielement-1, self.inode, :] = [oxx, oyy, ozz, txy, tyz, txz, o1, o2, o3, ovm]

        #print('eid=%i node_id=%i exx=%s' % (eid, node_id, str(oxx)))
        self.element_node[self.itotal, :] = [eid, node_id]
        #self.element_node[self.ielement-1, self.inode-1, :] = [eid, node_id]
        self.itotal += 1

    @property
    def nnodes_per_element(self) -> int:
        return self.nnodes_per_element_no_centroid + 1

    @property
    def nnodes_per_element_no_centroid(self) -> int:
        if self.element_type == 39: # CTETRA
            nnodes = 4
        elif self.element_type == 67: # CHEXA
            nnodes = 8
        elif self.element_type == 68: # CPENTA
            nnodes = 6
        elif self.element_type == 255: # CPYRAM
            nnodes = 5
        else:
            raise NotImplementedError(f'element_name={self.element_name} self.element_type={self.element_type}')
        return nnodes

    def get_stats(self, short=False) -> List[str]:
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

        if self.nonlinear_factor not in (None, np.nan):  # transient
            msg.append('  type=%s ntimes=%i nelements=%i nnodes=%i\n  nnodes_per_element=%s (including centroid)\n'
                       % (self.__class__.__name__, ntimes, nelements, nnodes, nnodes_per_element))
            ntimes_word = 'ntimes'
        else:
            msg.append('  type=%s nelements=%i nnodes=%i\n  nodes_per_element=%i (including centroid)\n'
                       % (self.__class__.__name__, nelements, nnodes, nnodes_per_element))
            ntimes_word = '1'
        msg.append('  eType, cid\n')
        headers = self.get_headers()
        n = len(headers)
        msg.append('  data: [%s, nnodes, %i] where %i=[%s]\n' % (ntimes_word, n, n, str(', '.join(headers))))
        msg.append('  element_node.shape = %s\n' % str(self.element_node.shape).replace('L', ''))
        msg.append('  element_cid.shape = %s\n' % str(self.element_cid.shape).replace('L', ''))
        msg.append('  data.shape = %s\n' % str(self.data.shape).replace('L', ''))
        msg.append('  element name: %s\n' % self.element_name)
        msg += self.get_data_code()
        #print(''.join(msg))
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

    def write_f06(self, f06_file, header=None, page_stamp='PAGE %s', page_num=1, is_mag_phase=False, is_sort1=True):
        if header is None:
            header = []
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
            f06_file.write(''.join(header + msg_temp))

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

            cnnodes = nnodes + 1
            for i, deid, node_id, doxx, doyy, dozz, dtxy, dtyz, dtxz, do1, do2, do3, dp, dovm in zip(
                    count(), eids2, nodes, oxx, oyy, ozz, txy, tyz, txz, o1, o2, o3, p, ovm):

                j = where(eids3 == deid)[0][0]
                cid = cids3[j]
                A = [[doxx, dtxy, dtxz],
                     [dtxy, doyy, dtyz],
                     [dtxz, dtyz, dozz]]
                (_lambda, v) = eigh(A)  # a hermitian matrix is a symmetric-real matrix

                # o1-max
                # o2-mid
                # o3-min
                assert do1 >= do2 >= do3, 'o1 >= o2 >= o3; eid=%s o1=%e o2=%e o3=%e' % (deid, do1, do2, do3)
                [oxxi, oyyi, ozzi, txyi, tyzi, txzi, o1i, o2i, o3i, pi, ovmi] = write_floats_13e(
                    [doxx, doyy, dozz, dtxy, dtyz, dtxz, do1, do2, do3, dp, dovm])

                if i % cnnodes == 0:
                    f06_file.write('0  %8s    %8iGRID CS  %i GP\n' % (deid, cid, nnodes))
                    f06_file.write(
                        '0              %8s  X  %-13s  XY  %-13s   A  %-13s  LX%5.2f%5.2f%5.2f  %-13s   %s\n'
                        '               %8s  Y  %-13s  YZ  %-13s   B  %-13s  LY%5.2f%5.2f%5.2f\n'
                        '               %8s  Z  %-13s  ZX  %-13s   C  %-13s  LZ%5.2f%5.2f%5.2f\n'
                        % ('CENTER', oxxi, txyi, o1i, v[0, 1], v[0, 2], v[0, 0], pi, ovmi,
                           '', oyyi, tyzi, o2i, v[1, 1], v[1, 2], v[1, 0],
                           '', ozzi, txzi, o3i, v[2, 1], v[2, 2], v[2, 0]))
                else:
                    f06_file.write(
                        '0              %8s  X  %-13s  XY  %-13s   A  %-13s  LX%5.2f%5.2f%5.2f  %-13s   %s\n'
                        '               %8s  Y  %-13s  YZ  %-13s   B  %-13s  LY%5.2f%5.2f%5.2f\n'
                        '               %8s  Z  %-13s  ZX  %-13s   C  %-13s  LZ%5.2f%5.2f%5.2f\n'
                        % (node_id, oxxi, txyi, o1i, v[0, 1], v[0, 2], v[0, 0], pi, ovmi,
                           '', oyyi, tyzi, o2i, v[1, 1], v[1, 2], v[1, 0],
                           '', ozzi, txzi, o3i, v[2, 1], v[2, 2], v[2, 0]))
                i += 1
            f06_file.write(page_stamp % page_num)
            page_num += 1
        return page_num - 1

    def write_op2(self, op2, op2_ascii, itable, new_result,
                  date, is_mag_phase=False, endian='>'):
        """writes an OP2"""
        import inspect
        frame = inspect.currentframe()
        call_frame = inspect.getouterframes(frame, 2)
        op2_ascii.write('%s.write_op2: %s\n' % (self.__class__.__name__, call_frame[1][3]))

        if itable == -1:
            #print('***************', itable)
            self._write_table_header(op2, op2_ascii, date)
            itable = -3

        #if isinstance(self.nonlinear_factor, float):
            #op2_format = '%sif' % (7 * self.ntimes)
            #raise NotImplementedError()
        #else:
            #op2_format = 'i21f'
        #s = Struct(op2_format)
        nnodes_expected = self.nnodes

        eids2 = self.element_node[:, 0]
        nodes = self.element_node[:, 1]

        eids3 = self.element_cid[:, 0]
        cids3 = self.element_cid[:, 1]

        # table 4 info
        #ntimes = self.data.shape[0]
        nnodes = self.data.shape[1]
        nelements = len(np.unique(eids2))

        # 21 = 1 node, 3 principal, 6 components, 9 vectors, 2 p/ovm
        #ntotal = ((nnodes * 21) + 1) + (nelements * 4)

        nnodes_expected = self.nnodes_per_element_no_centroid
        ntotali = 4 + 21 * nnodes_expected
        ntotali = self.num_wide
        ntotal = ntotali * nelements

        #print('shape = %s' % str(self.data.shape))
        assert nnodes > 1, nnodes
        #assert self.ntimes == 1, self.ntimes

        #device_code = self.device_code
        op2_ascii.write('  ntimes = %s\n' % self.ntimes)

        #fmt = '%2i %6f'
        #print('ntotal=%s' % (ntotal))
        #assert ntotal == 193, ntotal

        if self.is_sort1:
            #op2_format = endian + b'2i6f'
            struct1 = Struct(endian + b'ii4si')
            struct2 = Struct(endian + b'i20f')
        else:
            raise NotImplementedError('SORT2')

        cen = b'GRID'
        op2_ascii.write('nelements=%i\n' % nelements)
        for itime in range(self.ntimes):
            self._write_table_3(op2, op2_ascii, new_result, itable, itime)

            # record 4
            #print('stress itable = %s' % itable)
            itable -= 1
            header = [4, itable, 4,
                      4, 1, 4,
                      4, 0, 4,
                      4, ntotal, 4,
                      4 * ntotal]
            op2.write(pack('%ii' % len(header), *header))
            op2_ascii.write('r4 [4, 0, 4]\n')
            op2_ascii.write('r4 [4, %s, 4]\n' % (itable))
            op2_ascii.write('r4 [4, %i, 4]\n' % (4 * ntotal))

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

            #print('eids3', eids3)
            cnnodes = nnodes_expected + 1
            for i, deid, node_id, doxx, doyy, dozz, dtxy, dtyz, dtxz, do1, do2, do3, dp, dovm in zip(
                    count(), eids2, nodes, oxx, oyy, ozz, txy, tyz, txz, o1, o2, o3, p, ovm):
                #print('  eid =', deid, node_id)

                j = where(eids3 == deid)[0]
                assert len(j) > 0, j
                cid = cids3[j][0]
                A = [[doxx, dtxy, dtxz],
                     [dtxy, doyy, dtyz],
                     [dtxz, dtyz, dozz]]
                (_lambda, v) = eigh(A)  # a hermitian matrix is a symmetric-real matrix

                #node_id, oxxi, txyi, o1i, v[0, 1], v[0, 2], v[0, 0], pi, ovmi,
                    #'', oyyi, tyzi, o2i, v[1, 1], v[1, 2], v[1, 0],
                    #'', ozzi, txzi, o3i, v[2, 1], v[2, 2], v[2, 0]
                    #(grid_device, sxx, sxy, s1, a1, a2, a3, pressure, svm,
                     #syy, syz, s2, b1, b2, b3,
                     #szz, sxz, s3, c1, c2, c3)

                if i % cnnodes == 0:
                    data = [deid * 10 + self.device_code, cid, cen, nnodes_expected]
                    op2_ascii.write('  eid=%s cid=%s cen=%s nnodes = %s\n' % tuple(data))
                    op2.write(
                        struct1.pack(*data)
                        #pack(b'2i 4s i', *data)
                    )
                #else:
                op2_ascii.write('    nid=%i\n' % node_id)

                data = [node_id,
                        doxx, dtxy, do1, v[0, 1], v[0, 2], v[0, 0], dp, dovm,
                        doyy, dtyz, do2, v[1, 1], v[1, 2], v[1, 0],
                        dozz, dtxz, do3, v[2, 1], v[2, 2], v[2, 0], ]
                op2_ascii.write('      oxx, txy, o1, v01, v02, v00, p, ovm = %s\n' % data[1:8])
                op2_ascii.write('      oyy, tyz, o2, v11, v12, v10         = %s\n' % data[8:14])
                op2_ascii.write('      ozz, txz, o3, v21, v22, v20         = %s\n' % data[14:])
                op2.write(struct2.pack(*data))
                i += 1

            itable -= 1
            header = [4 * ntotal,]
            op2.write(pack('i', *header))
            op2_ascii.write('footer = %s\n' % header)
            new_result = False
        return itable


class RealSolidStressArray(RealSolidArray, StressObject):
    def __init__(self, data_code, is_sort1, isubcase, dt):
        RealSolidArray.__init__(self, data_code, is_sort1, isubcase, dt)
        StressObject.__init__(self, data_code, isubcase)

    def get_headers(self) -> List[str]:
        if self.is_von_mises:
            von_mises = 'von_mises'
        else:
            von_mises = 'max_shear'
        headers = ['oxx', 'oyy', 'ozz', 'txy', 'tyz', 'txz', 'omax', 'omid', 'omin', von_mises]
        return headers


class RealSolidStrainArray(RealSolidArray, StrainObject):
    def __init__(self, data_code, is_sort1, isubcase, dt):
        RealSolidArray.__init__(self, data_code, is_sort1, isubcase, dt)
        StrainObject.__init__(self, data_code, isubcase)

    def get_headers(self) -> List[str]:
        if self.is_von_mises:
            von_mises = 'von_mises'
        else:
            von_mises = 'max_shear'
        headers = ['exx', 'eyy', 'ezz', 'exy', 'eyz', 'exz', 'emax', 'emid', 'emin', von_mises]
        return headers

def _get_solid_msgs(self):
    if self.is_von_mises:
        von_mises = 'VON MISES'
    else:
        von_mises = 'MAX SHEAR'

    if self.is_stress:

        base_msg = [
            '0                CORNER        ------CENTER AND CORNER POINT STRESSES---------       DIR.  COSINES       MEAN                   \n',
            '  ELEMENT-ID    GRID-ID        NORMAL              SHEAR             PRINCIPAL       -A-  -B-  -C-     PRESSURE       %s \n' % von_mises]
        tetra_msg = ['                   S T R E S S E S   I N    T E T R A H E D R O N   S O L I D   E L E M E N T S   ( C T E T R A )\n', ]
        penta_msg = ['                    S T R E S S E S   I N   P E N T A H E D R O N   S O L I D   E L E M E N T S   ( P E N T A )\n', ]
        hexa_msg = ['                      S T R E S S E S   I N   H E X A H E D R O N   S O L I D   E L E M E N T S   ( H E X A )\n', ]
        pyram_msg = ['                      S T R E S S E S   I N   P Y R A M I D   S O L I D   E L E M E N T S   ( P Y R A M )\n', ]
    else:
        base_msg = [
            '0                CORNER        ------CENTER AND CORNER POINT  STRAINS---------       DIR.  COSINES       MEAN                   \n',
            '  ELEMENT-ID    GRID-ID        NORMAL              SHEAR             PRINCIPAL       -A-  -B-  -C-     PRESSURE       %s \n' % von_mises]
        tetra_msg = ['                     S T R A I N S   I N    T E T R A H E D R O N   S O L I D   E L E M E N T S   ( C T E T R A )\n', ]
        penta_msg = ['                      S T R A I N S   I N   P E N T A H E D R O N   S O L I D   E L E M E N T S   ( P E N T A )\n', ]
        hexa_msg = ['                        S T R A I N S   I N   H E X A H E D R O N   S O L I D   E L E M E N T S   ( H E X A )\n', ]
        pyram_msg = ['                        S T R A I N S   I N   P Y R A M I D   S O L I D   E L E M E N T S   ( P Y R A M )\n', ]

    tetra_msg += base_msg
    penta_msg += base_msg
    hexa_msg += base_msg
    return tetra_msg, penta_msg, hexa_msg, pyram_msg

def _get_f06_header_nnodes(self, is_mag_phase=True):
    tetra_msg, penta_msg, hexa_msg, pyram_msg = _get_solid_msgs(self)
    if self.element_type == 39: # CTETRA
        msg = tetra_msg
        nnodes = 4
    elif self.element_type == 67: # CHEXA
        msg = hexa_msg
        nnodes = 8
    elif self.element_type == 68: # CPENTA
        msg = penta_msg
        nnodes = 6
    elif self.element_type == 255: # CPYRAM
        msg = pyram_msg
        nnodes = 5
    else:  # pragma: no cover
        msg = f'element_name={self.element_name} self.element_type={self.element_type}'
        raise NotImplementedError(msg)
    return nnodes, msg
