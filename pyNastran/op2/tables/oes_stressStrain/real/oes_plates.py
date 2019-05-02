# coding: utf-8
#pylint disable=C0103
from __future__ import (nested_scopes, generators, division, absolute_import,
                        print_function, unicode_literals)
from itertools import count
import numpy as np

from pyNastran.utils.numpy_utils import integer_types
from pyNastran.op2.tables.oes_stressStrain.real.oes_objects import (
    StressObject, StrainObject, OES_Object)
from pyNastran.f06.f06_formatting import write_floats_13e, _eigenvalue_header


class RealPlateArray(OES_Object):
    def __init__(self, data_code, is_sort1, isubcase, dt):
        OES_Object.__init__(self, data_code, isubcase, apply_data_code=False)
        #self.code = [self.format_code, self.sort_code, self.s_code]

        #self.ntimes = 0  # or frequency/mode
        #self.ntotal = 0
        self.ielement = 0
        self.nelements = 0  # result specific
        self.nnodes = None

        #if is_sort1:
            #pass
        #else:
            #raise NotImplementedError('SORT2')

    @property
    def is_real(self):
        return True

    @property
    def is_complex(self):
        return False

    def _reset_indices(self):
        self.itotal = 0
        self.ielement = 0

    def get_headers(self):
        raise NotImplementedError('%s needs to implement get_headers' % self.__class__.__name__)

    def is_bilinear(self):
        if self.element_type in [33, 74]:  # CQUAD4, CTRIA3
            return False
        elif self.element_type in [144, 64, 82, 70, 75]:  # CQUAD4
            return True
        else:
            raise NotImplementedError('name=%s type=%s' % (self.element_name, self.element_type))

    def build(self):
        """sizes the vectorized attributes of the RealPlateArray"""
        #print("self.ielement = %s" % self.ielement)
        #print('ntimes=%s nelements=%s ntotal=%s' % (self.ntimes, self.nelements, self.ntotal))
        assert self.ntimes > 0, 'ntimes=%s' % self.ntimes
        assert self.nelements > 0, 'nelements=%s' % self.nelements
        assert self.ntotal > 0, 'ntotal=%s' % self.ntotal
        #self.names = []
        if self.element_type in [33, 74]:
            nnodes_per_element = 1
        elif self.element_type == 144:
            nnodes_per_element = 5
        elif self.element_type == 64:  # CQUAD8
            nnodes_per_element = 5
        elif self.element_type == 82:  # CQUADR
            nnodes_per_element = 5
        elif self.element_type == 70:  # CTRIAR
            nnodes_per_element = 4
        elif self.element_type == 75:  # CTRIA6
            nnodes_per_element = 4
        else:
            raise NotImplementedError('name=%r type=%s' % (self.element_name, self.element_type))

        #print('nnodes_per_element[%s, %s] = %s' % (
            #self.isubcase, self.element_type, nnodes_per_element))
        self.nnodes = nnodes_per_element
        #self.nelements //= nnodes_per_element
        self.nelements //= self.ntimes
        self.itime = 0
        self.ielement = 0
        self.itotal = 0
        #self.ntimes = 0
        #self.nelements = 0
        self.is_built = True

        #print("***name=%s type=%s nnodes_per_element=%s ntimes=%s nelements=%s ntotal=%s" % (
            #self.element_name, self.element_type, nnodes_per_element, self.ntimes,
            #self.nelements, self.ntotal))
        dtype = 'float32'
        if isinstance(self.nonlinear_factor, integer_types):
            dtype = 'int32'

        _times = np.zeros(self.ntimes, dtype=dtype)
        element_node = np.zeros((self.ntotal, 2), dtype='int32')

        #[fiber_dist, oxx, oyy, txy, angle, majorP, minorP, ovm]
        data = np.zeros((self.ntimes, self.ntotal, 8), dtype='float32')
        if self.load_as_h5:
            #for key, value in sorted(self.data_code.items()):
                #print(key, value)
            group = self._get_result_group()
            self._times = group.create_dataset('_times', data=_times)
            self.element_node = group.create_dataset('element_node', data=element_node)
            self.data = group.create_dataset('data', data=data)
        else:
            self._times = _times
            self.element_node = element_node
            self.data = data

    def build_dataframe(self):
        """creates a pandas dataframe"""
        import pandas as pd
        headers = self.get_headers()

        nelements = self.element_node.shape[0] // 2
        if self.is_fiber_distance:
            fiber_distance = ['Top', 'Bottom'] * nelements
        else:
            fiber_distance = ['Mean', 'Curvature'] * nelements
        fd = np.array(fiber_distance, dtype='unicode')
        element_node = [self.element_node[:, 0], self.element_node[:, 1], fd]

        if self.nonlinear_factor not in (None, np.nan):
            column_names, column_values = self._build_dataframe_transient_header()
            self.data_frame = pd.Panel(self.data, items=column_values,
                                       major_axis=element_node, minor_axis=headers).to_frame()
            self.data_frame.columns.names = column_names
            self.data_frame.index.names = ['ElementID', 'NodeID', 'Location', 'Item']
        else:
            # option B - nice!
            df1 = pd.DataFrame(element_node).T
            df1.columns = ['ElementID', 'NodeID', 'Location']
            df2 = pd.DataFrame(self.data[0])
            df2.columns = headers
            self.data_frame = df1.join(df2)
        self.data_frame = self.data_frame.reset_index().replace(
            {'NodeID': {0:'CEN'}}).set_index(['ElementID', 'NodeID', 'Location'])
        #print(self.data_frame)

    def __eq__(self, table):
        assert self.is_sort1 == table.is_sort1
        self._eq_header(table)
        if not np.array_equal(self.data, table.data):
            msg = 'table_name=%r class_name=%s\n' % (self.table_name, self.__class__.__name__)
            msg += '%s\n' % str(self.code_information())
            i = 0
            for itime in range(self.ntimes):
                for ie, element_nodei in enumerate(self.element_node):
                    (eid, nid) = element_nodei
                    t1 = self.data[itime, ie, :]
                    t2 = table.data[itime, ie, :]
                    (fiber_dist1, oxx1, oyy1, txy1, angle1, major_p1, minor_p1, ovm1) = t1
                    (fiber_dist2, oxx2, oyy2, txy2, angle2, major_p2, minor_p2, ovm2) = t2

                    # vm stress can be NaN for some reason...
                    if not np.array_equal(t1[:-1], t2[:-1]):
                        msg += '(%s, %s)    (%s, %s, %s, %s, %s, %s, %s, %s)  (%s, %s, %s, %s, %s, %s, %s, %s)\n' % (
                            eid, nid,
                            fiber_dist1, oxx1, oyy1, txy1, angle1, major_p1, minor_p1, ovm1,
                            fiber_dist2, oxx2, oyy2, txy2, angle2, major_p2, minor_p2, ovm2)
                        i += 1
                        if i > 10:
                            print(msg)
                            raise ValueError(msg)
                #print(msg)
                if i > 0:
                    raise ValueError(msg)
        return True

    def add_new_eid_sort1(self, dt, eid, node_id, fiber_dist, oxx, oyy, txy, angle,
                          major_principal, minor_principal, ovm):
        assert isinstance(eid, integer_types), eid
        assert isinstance(node_id, integer_types), node_id
        self._times[self.itime] = dt
        #assert self.itotal == 0, oxx
        self.element_node[self.itotal, :] = [eid, node_id]
        self.data[self.itime, self.itotal, :] = [fiber_dist, oxx, oyy, txy, angle,
                                                 major_principal, minor_principal, ovm]
        self.itotal += 1
        self.ielement += 1

    def add_new_node_sort1(self, dt, eid, node_id, fiber_dist, oxx, oyy, txy, angle,
                           major_principal, minor_principal, ovm):
        self.add_sort1(dt, eid, node_id, fiber_dist, oxx, oyy, txy, angle,
                       major_principal, minor_principal, ovm)

    def add_sort1(self, dt, eid, node_id, fiber_dist, oxx, oyy, txy, angle,
                  major_principal, minor_principal, ovm):
        assert eid is not None, eid
        assert isinstance(eid, integer_types) and eid > 0, 'dt=%s eid=%s' % (dt, eid)
        assert isinstance(node_id, integer_types), node_id
        self.element_node[self.itotal, :] = [eid, node_id]
        self.data[self.itime, self.itotal, :] = [fiber_dist, oxx, oyy, txy, angle,
                                                 major_principal, minor_principal, ovm]
        self.itotal += 1

    def get_stats(self, short=False):
        if not self.is_built:
            return [
                '<%s>\n' % self.__class__.__name__,
                '  ntimes: %i\n' % self.ntimes,
                '  ntotal: %i\n' % self.ntotal,
            ]

        nelements = self.nelements
        ntimes = self.ntimes
        nnodes = self.nnodes
        ntotal = self.ntotal
        nlayers = 2
        nelements = self.ntotal // self.nnodes // 2

        msg = []
        if self.nonlinear_factor not in (None, np.nan):  # transient
            msgi = '  type=%s ntimes=%i nelements=%i nnodes_per_element=%i nlayers=%i ntotal=%i\n' % (
                self.__class__.__name__, ntimes, nelements, nnodes, nlayers, ntotal)
            ntimes_word = 'ntimes'
        else:
            msgi = '  type=%s nelements=%i nnodes_per_element=%i nlayers=%i ntotal=%i\n' % (
                self.__class__.__name__, nelements, nnodes, nlayers, ntotal)
            ntimes_word = '1'
        msg.append(msgi)
        headers = self.get_headers()
        n = len(headers)
        msg.append('  data: [%s, ntotal, %i] where %i=[%s]\n' % (ntimes_word, n, n,
                                                                 str(', '.join(headers))))
        msg.append('  element_node.shape = %s\n' % str(self.element_node.shape).replace('L', ''))
        msg.append('  data.shape=%s\n' % str(self.data.shape).replace('L', ''))
        msg.append('  element type: %s\n' % self.element_name)
        msg.append('  s_code: %s\n' % self.s_code)
        msg += self.get_data_code()
        return msg

    def get_element_index(self, eids):
        # elements are always sorted; nodes are not
        itot = np.searchsorted(eids, self.element_node[:, 0])  #[0]
        return itot

    def eid_to_element_node_index(self, eids):
        ind = np.ravel([np.searchsorted(self.element_node[:, 0] == eid) for eid in eids])
        #ind = searchsorted(eids, self.element)
        #ind = ind.reshape(ind.size)
        #ind.sort()
        return ind

    def write_f06(self, f06_file, header=None, page_stamp='PAGE %s',
                  page_num=1, is_mag_phase=False, is_sort1=True):
        if header is None:
            header = []
        msg, nnodes, cen = _get_plate_msg(self)

        # write the f06
        ntimes = self.data.shape[0]

        eids = self.element_node[:, 0]
        nids = self.element_node[:, 1]

        #cen_word = 'CEN/%i' % nnodes
        cen_word = cen
        for itime in range(ntimes):
            dt = self._times[itime]
            header = _eigenvalue_header(self, header, itime, ntimes, dt)
            f06_file.write(''.join(header + msg))

            #print("self.data.shape=%s itime=%s ieids=%s" % (str(self.data.shape), itime, str(ieids)))

            #[fiber_dist, oxx, oyy, txy, angle, majorP, minorP, ovm]
            fiber_dist = self.data[itime, :, 0]
            oxx = self.data[itime, :, 1]
            oyy = self.data[itime, :, 2]
            txy = self.data[itime, :, 3]
            angle = self.data[itime, :, 4]
            major_principal = self.data[itime, :, 5]
            minor_principal = self.data[itime, :, 6]
            ovm = self.data[itime, :, 7]

            for (i, eid, nid, fdi, oxxi, oyyi, txyi, anglei, major, minor, ovmi) in zip(
                 count(), eids, nids, fiber_dist, oxx, oyy, txy, angle, major_principal, minor_principal, ovm):
                [fdi, oxxi, oyyi, txyi, major, minor, ovmi] = write_floats_13e(
                    [fdi, oxxi, oyyi, txyi, major, minor, ovmi])
                ilayer = i % 2
                # tria3
                if self.element_type in [33, 74]:  # CQUAD4, CTRIA3
                    if ilayer == 0:
                        f06_file.write('0  %6i   %-13s     %-13s  %-13s  %-13s   %8.4f   %-13s   %-13s  %s\n' % (
                            eid, fdi, oxxi, oyyi, txyi, anglei, major, minor, ovmi))
                    else:
                        f06_file.write('   %6s   %-13s     %-13s  %-13s  %-13s   %8.4f   %-13s   %-13s  %s\n' % (
                            '', fdi, oxxi, oyyi, txyi, anglei, major, minor, ovmi))

                elif self.element_type in [64, 70, 75, 82, 144]:  # CQUAD8, CTRIAR, CTRIA6, CQUADR, CQUAD4
                    # bilinear
                    if nid == 0 and ilayer == 0:  # CEN
                        f06_file.write('0  %8i %8s  %-13s  %-13s %-13s %-13s   %8.4f  %-13s %-13s %s\n' % (
                            eid, cen_word, fdi, oxxi, oyyi, txyi, anglei, major, minor, ovmi))
                    elif ilayer == 0:
                        f06_file.write('   %8s %8i  %-13s  %-13s %-13s %-13s   %8.4f  %-13s %-13s %s\n' % (
                            '', nid, fdi, oxxi, oyyi, txyi, anglei, major, minor, ovmi))
                    elif ilayer == 1:
                        f06_file.write('   %8s %8s  %-13s  %-13s %-13s %-13s   %8.4f  %-13s %-13s %s\n\n' % (
                            '', '', fdi, oxxi, oyyi, txyi, anglei, major, minor, ovmi))
                else:
                    msg = 'element_name=%s self.element_type=%s' % (
                        self.element_name, self.element_type)
                    raise NotImplementedError(msg)

            f06_file.write(page_stamp % page_num)
            page_num += 1
        return page_num - 1

    def get_nnodes_bilinear(self):
        """gets the number of nodes and whether or not the element has bilinear results"""
        is_bilinear = False
        if self.element_type == 74:
            nnodes = 3
        elif self.element_type == 33:
            nnodes = 4
        elif self.element_type == 144:
            nnodes = 4
            is_bilinear = True
        elif self.element_type == 82:  # CQUADR
            nnodes = 4
            is_bilinear = True
        elif self.element_type == 64:  # CQUAD8
            nnodes = 4
            is_bilinear = True
        elif self.element_type == 75:  # CTRIA6
            nnodes = 3
            is_bilinear = True
        elif self.element_type == 70:  # CTRIAR
            nnodes = 3
            is_bilinear = True
        else:
            raise NotImplementedError('name=%s type=%s' % (self.element_name, self.element_type))
        return nnodes, is_bilinear

    def write_op2(self, op2, op2_ascii, itable, new_result,
                  date, is_mag_phase=False, endian='>'):
        """writes an OP2"""
        import inspect
        from struct import Struct, pack
        frame = inspect.currentframe()
        call_frame = inspect.getouterframes(frame, 2)
        op2_ascii.write('%s.write_op2: %s\n' % (self.__class__.__name__, call_frame[1][3]))

        if itable == -1:
            self._write_table_header(op2, op2_ascii, date)
            itable = -3

        nnodes, is_bilinear = self.get_nnodes_bilinear()
        if is_bilinear:
            nnodes_all = nnodes + 1
            ntotal = 2 + 17 * nnodes_all
        else:
            nnodes_all = nnodes
        #print("nnodes_all =", nnodes_all)
        cen_word_ascii = 'CEN/%i' % nnodes
        cen_word = b'CEN/%i' % nnodes

        #msg.append('  element_node.shape = %s\n' % str(self.element_node.shape).replace('L', ''))
        #msg.append('  data.shape=%s\n' % str(self.data.shape).replace('L', ''))

        eids = self.element_node[:, 0]
        nids = self.element_node[:, 1]

        eids_device = eids * 10 + self.device_code

        nelements = len(np.unique(eids))
        #print('nelements =', nelements)
        # 21 = 1 node, 3 principal, 6 components, 9 vectors, 2 p/ovm
        #ntotal = ((nnodes * 21) + 1) + (nelements * 4)

        ntotali = self.num_wide
        ntotal = ntotali * nelements

        #print('shape = %s' % str(self.data.shape))
        assert nnodes > 1, nnodes
        #assert self.ntimes == 1, self.ntimes

        device_code = self.device_code
        op2_ascii.write('  ntimes = %s\n' % self.ntimes)

        #fmt = '%2i %6f'
        #print('ntotal=%s' % (ntotal))
        #assert ntotal == 193, ntotal

        #[fiber_dist, oxx, oyy, txy, angle, majorP, minorP, ovm]
        op2_ascii.write('  #elementi = [eid_device, fd1, sx1, sy1, txy1, angle1, major1, minor1, vm1,\n')
        op2_ascii.write('  #                        fd2, sx2, sy2, txy2, angle2, major2, minor2, vm2,]\n')

        if self.is_sort1:
            struct1 = Struct(endian + b'i16f')
        else:
            raise NotImplementedError('SORT2')

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

            fiber_dist = self.data[itime, :, 0]
            oxx = self.data[itime, :, 1]
            oyy = self.data[itime, :, 2]
            txy = self.data[itime, :, 3]
            angle = self.data[itime, :, 4]
            major_principal = self.data[itime, :, 5]
            minor_principal = self.data[itime, :, 6]
            ovm = self.data[itime, :, 7]

            nwide = 0
            for (i, eid_device, eid, nid, fdi, oxxi, oyyi, txyi, anglei, major, minor, ovmi) in zip(
                 count(), eids_device, eids, nids, fiber_dist, oxx, oyy, txy, angle, major_principal, minor_principal, ovm):
                #[fdi, oxxi, oyyi, txyi, major, minor, ovmi] = write_floats_13e(
                    #[fdi, oxxi, oyyi, txyi, major, minor, ovmi])
                ilayer = i % 2
                # tria3
                if self.element_type in [33, 74]:  # CQUAD4, CTRIA3
                    if ilayer == 0:
                        #print([eid, fdi, oxxi, oyyi, txyi, anglei, major, minor, ovmi])
                        data = [eid_device, fdi, oxxi, oyyi, txyi, anglei, major, minor, ovmi]
                        op2.write(pack('i8f', *data))
                        op2_ascii.write('eid=%s ilayer=0 data=%s' % (eid, str(data[1:])))

                    else:
                        data = [fdi, oxxi, oyyi, txyi, anglei, major, minor, ovmi]
                        op2.write(pack('8f', *data))
                        op2_ascii.write('eid=%s ilayer=1 data=%s' % (eid, str(data[1:])))
                    #print('eid=%-2s ilayer=%s data=%s' % (eid_device, ilayer, str(data[1:])))

                elif self.element_type in [64, 70, 75, 82, 144]:  # CQUAD8, CTRIAR, CTRIA6, CQUADR, CQUAD4
                    # bilinear
                    if nid == 0 and ilayer == 0:  # CEN
                        #print()
                        data = [eid_device, cen_word, nid,
                                fdi, oxxi, oyyi, txyi, anglei, major, minor, ovmi]
                        op2.write(pack('i 4s i 8f', *data))
                        op2_ascii.write('0  %8i %8s  %-13s  %-13s %-13s %-13s   %8.4f  %-13s %-13s %s\n' % (
                            eid, cen_word_ascii, fdi, oxxi, oyyi, txyi, anglei, major, minor, ovmi))
                    elif ilayer == 0:
                        data = [nid, fdi, oxxi, oyyi, txyi, anglei, major, minor, ovmi]
                        op2.write(pack('i 8f', *data))
                        op2_ascii.write('   %8s %8i  %-13s  %-13s %-13s %-13s   %8.4f  %-13s %-13s %s\n' % (
                            '', nid, fdi, oxxi, oyyi, txyi, anglei, major, minor, ovmi))
                    elif ilayer == 1:
                        data = [fdi, oxxi, oyyi, txyi, anglei, major, minor, ovmi]
                        op2.write(pack('8f', *data))
                        op2_ascii.write('   %8s %8s  %-13s  %-13s %-13s %-13s   %8.4f  %-13s %-13s %s\n\n' % (
                            '', '', fdi, oxxi, oyyi, txyi, anglei, major, minor, ovmi))
                    else:  # pragma: no cover
                        raise RuntimeError()
                    #print('eid=%-2s nid=%s ilayer=%s data=%s' % (eid, nid, ilayer, str(data[1:])))
                else:
                    msg = 'element_name=%s self.element_type=%s' % (
                        self.element_name, self.element_type)
                    raise NotImplementedError(msg)
                nwide += len(data)

            assert nwide == ntotal, "nwide=%s ntotal=%s" % (nwide, ntotal)
            itable -= 1
            header = [4 * ntotal,]
            op2.write(pack('i', *header))
            op2_ascii.write('footer = %s\n' % header)
            new_result = False
        return itable


class RealPlateStressArray(RealPlateArray, StressObject):
    def __init__(self, data_code, is_sort1, isubcase, dt):
        RealPlateArray.__init__(self, data_code, is_sort1, isubcase, dt)
        StressObject.__init__(self, data_code, isubcase)

    def get_headers(self):
        if self.is_fiber_distance:
            fiber_dist = 'fiber_distance'
        else:
            fiber_dist = 'fiber_curvature'

        if self.is_von_mises:
            ovm = 'von_mises'
        else:
            ovm = 'max_shear'
        headers = [fiber_dist, 'oxx', 'oyy', 'txy', 'angle', 'omax', 'omin', ovm]
        return headers


class RealPlateStrainArray(RealPlateArray, StrainObject):
    """
    used for:
     - RealPlateStressArray
     - RealPlateStrainArray
    """
    def __init__(self, data_code, is_sort1, isubcase, dt):
        RealPlateArray.__init__(self, data_code, is_sort1, isubcase, dt)
        StrainObject.__init__(self, data_code, isubcase)

    def get_headers(self):
        if self.is_fiber_distance:
            fiber_dist = 'fiber_distance'
        else:
            fiber_dist = 'fiber_curvature'

        if self.is_von_mises:
            ovm = 'von_mises'
        else:
            ovm = 'max_shear'
        headers = [fiber_dist, 'exx', 'eyy', 'exy', 'angle', 'emax', 'emin', ovm]
        return headers


def _get_plate_msg(self):
    if self.is_von_mises:
        von_mises = 'VON MISES'
    else:
        von_mises = 'MAX SHEAR'

    if self.is_stress:
        if self.is_fiber_distance:
            quad_msg_temp = ['    ELEMENT              FIBER            STRESSES IN ELEMENT COORD SYSTEM         PRINCIPAL STRESSES (ZERO SHEAR)               \n',
                             '      ID      GRID-ID   DISTANCE        NORMAL-X      NORMAL-Y      SHEAR-XY      ANGLE        MAJOR         MINOR       %s \n' % von_mises]
            tri_msg_temp = ['  ELEMENT      FIBER               STRESSES IN ELEMENT COORD SYSTEM             PRINCIPAL STRESSES (ZERO SHEAR)                 \n',
                            '    ID.       DISTANCE           NORMAL-X       NORMAL-Y      SHEAR-XY       ANGLE         MAJOR           MINOR        %s\n' % von_mises]
        else:
            quad_msg_temp = ['    ELEMENT              FIBER            STRESSES IN ELEMENT COORD SYSTEM         PRINCIPAL STRESSES (ZERO SHEAR)               \n',
                             '      ID      GRID-ID  CURVATURE        NORMAL-X      NORMAL-Y      SHEAR-XY      ANGLE        MAJOR         MINOR       %s \n' % von_mises]
            tri_msg_temp = ['  ELEMENT      FIBER               STRESSES IN ELEMENT COORD SYSTEM             PRINCIPAL STRESSES (ZERO SHEAR)                 \n',
                            '    ID.      CURVATURE           NORMAL-X       NORMAL-Y      SHEAR-XY       ANGLE         MAJOR           MINOR        %s\n' % von_mises]

        #is_bilinear = False
        cquad4_msg = ['                         S T R E S S E S   I N   Q U A D R I L A T E R A L   E L E M E N T S   ( Q U A D 4 )\n'] + tri_msg_temp
        cquad8_msg = ['                         S T R E S S E S   I N   Q U A D R I L A T E R A L   E L E M E N T S   ( Q U A D 8 )\n'] + tri_msg_temp
        cquadr_msg = ['                        S T R E S S E S   I N   Q U A D R I L A T E R A L   E L E M E N T S   ( Q U A D R )\n'] + tri_msg_temp

        ## TODO: can cquad8s be bilinear???
        #is_bilinear = True
        #cquadr_bilinear_msg = ['                         S T R E S S E S   I N   Q U A D R I L A T E R A L   E L E M E N T S   ( Q U A D R )        OPTION = BILIN  \n \n'] + quad_msg_temp
        cquad4_bilinear_msg = ['                         S T R E S S E S   I N   Q U A D R I L A T E R A L   E L E M E N T S   ( Q U A D 4 )        OPTION = BILIN  \n \n'] + quad_msg_temp

        #is_bilinear = False
        ctria3_msg = ['                           S T R E S S E S   I N   T R I A N G U L A R   E L E M E N T S   ( T R I A 3 )\n'] + tri_msg_temp
        ctria6_msg = ['                           S T R E S S E S   I N   T R I A N G U L A R   E L E M E N T S   ( T R I A 6 )\n'] + tri_msg_temp
        ctriar_msg = ['                           S T R E S S E S   I N   T R I A N G U L A R   E L E M E N T S   ( T R I A R )\n'] + tri_msg_temp
    else:
        if self.is_fiber_distance:
            quad_msg_temp = ['    ELEMENT              STRAIN            STRAINS IN ELEMENT COORD SYSTEM         PRINCIPAL  STRAINS (ZERO SHEAR)               \n',
                             '      ID      GRID-ID   DISTANCE        NORMAL-X      NORMAL-Y      SHEAR-XY      ANGLE        MAJOR         MINOR       %s \n' % von_mises]
            tri_msg_temp = ['  ELEMENT      FIBER                STRAINS IN ELEMENT COORD SYSTEM             PRINCIPAL  STRAINS (ZERO SHEAR)                 \n',
                            '    ID.       DISTANCE           NORMAL-X       NORMAL-Y      SHEAR-XY       ANGLE         MAJOR           MINOR        %s\n' % von_mises]
        else:
            quad_msg_temp = ['    ELEMENT              STRAIN            STRAINS IN ELEMENT COORD SYSTEM         PRINCIPAL  STRAINS (ZERO SHEAR)               \n',
                             '      ID      GRID-ID   CURVATURE       NORMAL-X      NORMAL-Y      SHEAR-XY      ANGLE        MAJOR         MINOR       %s \n' % von_mises]
            tri_msg_temp = ['  ELEMENT      STRAIN               STRAINS IN ELEMENT COORD SYSTEM             PRINCIPAL  STRAINS (ZERO SHEAR)                 \n',
                            '    ID.       CURVATURE          NORMAL-X       NORMAL-Y      SHEAR-XY       ANGLE         MAJOR           MINOR        %s\n' % von_mises]

        #is_bilinear = False
        cquad4_msg = ['                         S T R A I N S   I N   Q U A D R I L A T E R A L   E L E M E N T S   ( Q U A D 4 )\n'] + tri_msg_temp
        cquad8_msg = ['                         S T R A I N S   I N   Q U A D R I L A T E R A L   E L E M E N T S   ( Q U A D 8 )\n'] + tri_msg_temp
        cquadr_msg = ['                         S T R A I N S   I N   Q U A D R I L A T E R A L   E L E M E N T S   ( Q U A D R )\n'] + tri_msg_temp

        ## TODO: can cquad8s be bilinear???
        #cquadr_bilinear_msg = ['                           S T R A I N S   I N   Q U A D R I L A T E R A L   E L E M E N T S   ( Q U A D R )        OPTION = BILIN  \n \n'] + quad_msg_temp
        cquad4_bilinear_msg = ['                           S T R A I N S   I N   Q U A D R I L A T E R A L   E L E M E N T S   ( Q U A D 4 )        OPTION = BILIN  \n \n'] + quad_msg_temp

        cquadr_msg = ['                         S T R A I N S   I N   Q U A D R I L A T E R A L   E L E M E N T S   ( Q U A D R )\n'] + tri_msg_temp
        ctria3_msg = ['                             S T R A I N S   I N   T R I A N G U L A R   E L E M E N T S   ( T R I A 3 )\n'] + tri_msg_temp
        ctria6_msg = ['                             S T R A I N S   I N   T R I A N G U L A R   E L E M E N T S   ( T R I A 6 )\n'] + tri_msg_temp
        ctriar_msg = ['                             S T R A I N S   I N   T R I A N G U L A R   E L E M E N T S   ( T R I A R )\n'] + tri_msg_temp


    #is_bilinear = False
    if self.element_type == 74:
        msg = ctria3_msg
        nnodes = 3
        cen = 'CEN/3'
    elif self.element_type == 33:
        msg = cquad4_msg
        nnodes = 4
        cen = 'CEN/4'
    elif self.element_type == 144:
        msg = cquad4_bilinear_msg
        nnodes = 4
        #is_bilinear = True
        cen = 'CEN/4'
    elif self.element_type == 82:  # CQUADR
        msg = cquadr_msg
        nnodes = 4
        #is_bilinear = True
        cen = 'CEN/4'
    elif self.element_type == 64:  # CQUAD8
        msg = cquad8_msg
        nnodes = 4
        #is_bilinear = True
        cen = 'CEN/8'
    elif self.element_type == 75:  # CTRIA6
        msg = ctria6_msg
        nnodes = 3
        #is_bilinear = True
        cen = 'CEN/6'
    elif self.element_type == 70:  # CTRIAR
        msg = ctriar_msg
        nnodes = 3
        #is_bilinear = True
        cen = 'CEN/3'
    else:
        raise NotImplementedError('name=%s type=%s' % (self.element_name, self.element_type))
    return msg, nnodes, cen
