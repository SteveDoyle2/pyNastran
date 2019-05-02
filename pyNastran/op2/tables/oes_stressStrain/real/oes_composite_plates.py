from __future__ import (nested_scopes, generators, division, absolute_import,
                        print_function, unicode_literals)
import numpy as np
from numpy import zeros, searchsorted, unique, ravel

from pyNastran.utils.numpy_utils import integer_types
from pyNastran.op2.tables.oes_stressStrain.real.oes_objects import (
    StressObject, StrainObject, OES_Object)
from pyNastran.f06.f06_formatting import write_floats_12e, _eigenvalue_header


class RealCompositePlateArray(OES_Object):
    def __init__(self, data_code, is_sort1, isubcase, dt):
        OES_Object.__init__(self, data_code, isubcase, apply_data_code=False)
        #self.code = [self.format_code, self.sort_code, self.s_code]

        #self.ntimes = 0  # or frequency/mode
        #self.ntotal = 0
        self.ielement = 0
        self.nelements = 0  # result specific
        self.nnodes = None

        #if is_sort1:
            #if dt is not None:
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

    def _get_msgs(self):
        raise NotImplementedError('%s needs to implement _get_msgs' % self.__class__.__name__)

    def get_headers(self):
        raise NotImplementedError('%s needs to implement get_headers' % self.__class__.__name__)

    def build(self):
        """sizes the vectorized attributes of the RealCompositePlateArray"""
        assert self.ntimes > 0, 'ntimes=%s' % self.ntimes
        assert self.nelements > 0, 'nelements=%s' % self.nelements
        assert self.ntotal > 0, 'ntotal=%s' % self.ntotal

        if self.element_type == 95:  # CQUAD4
            nnodes_per_element = 1
        elif self.element_type == 96:  # CQUAD8
            nnodes_per_element = 1
        elif self.element_type == 97:  # CTRIA3
            nnodes_per_element = 1
        elif self.element_type == 98:  # CTRIA6
            nnodes_per_element = 1
        else:
            msg = 'element_name=%s element_type=%s' %(self.element_name, self.element_type)
            raise NotImplementedError(msg)

        self.nnodes = nnodes_per_element
        self.itime = 0
        self.ielement = 0
        self.itotal = 0
        self.is_built = True

        dtype = 'float32'
        if isinstance(self.nonlinear_factor, integer_types):
            dtype = 'int32'

        _times = zeros(self.ntimes, dtype=dtype)
        element_layer = zeros((self.ntotal, 2), dtype='int32')

        #[o11, o22, t12, t1z, t2z, angle, major, minor, ovm]
        data = zeros((self.ntimes, self.ntotal, 9), dtype='float32')

        if self.load_as_h5:
            #for key, value in sorted(self.data_code.items()):
                #print(key, value)
            group = self._get_result_group()
            self._times = group.create_dataset('_times', data=_times)
            self.element_layer = group.create_dataset('element_layer', data=element_layer)
            self.data = group.create_dataset('data', data=data)
        else:
            self._times = _times
            self.element_layer = element_layer
            self.data = data

    def build_dataframe(self):
        """
        major-axis - the axis

        mode   1     2   3
        freq  1.0   2.0 3.0
        T1
        T2
        T3
        R1
        R2
        R3

        major_axis / top = [
            [1, 2, 3],
            [1.0, 2.0, 3.0]
        ]
        minor_axis / headers = [T1, T2, T3, R1, R2, R3]
        name = mode
        """
        import pandas as pd
        headers = self.get_headers()
        element_layer = [self.element_layer[:, 0], self.element_layer[:, 1]]
        if self.nonlinear_factor not in (None, np.nan):
            column_names, column_values = self._build_dataframe_transient_header()
            self.data_frame = pd.Panel(self.data, items=column_values,
                                       major_axis=element_layer, minor_axis=headers).to_frame()
            self.data_frame.columns.names = column_names
            self.data_frame.index.names = ['ElementID', 'Layer', 'Item']
        else:
            self.data_frame = pd.Panel(self.data,
                                       major_axis=element_layer, minor_axis=headers).to_frame()
            self.data_frame.columns.names = ['Static']
            self.data_frame.index.names = ['ElementID', 'Layer', 'Item']

    def __eq__(self, table):
        assert self.is_sort1 == table.is_sort1
        self._eq_header(table)
        if not np.array_equal(self.element_layer, table.element_layer):
            assert self.element_node.shape == table.element_layer.shape, 'element_layer shape=%s table.shape=%s' % (
                self.element_layer.shape, table.element_layer.shape)
            msg = 'table_name=%r class_name=%s\n' % (self.table_name, self.__class__.__name__)
            msg += '%s\n' % str(self.code_information())
            msg += '(Eid, Layer)\n'
            for (eid, layer1), (eid2, layer2) in zip(self.element_layer, table.element_layer):
                msg += '(%s, %s)    (%s, %s)\n' % (eid, layer1, eid2, layer2)
            print(msg)
            raise ValueError(msg)
        if not np.array_equal(self.data, table.data):
            msg = 'table_name=%r class_name=%s\n' % (self.table_name, self.__class__.__name__)
            msg += '%s\n' % str(self.code_information())
            i = 0
            for itime in range(self.ntimes):
                for ie, e in enumerate(self.element_layer):
                    (eid, layer) = e
                    t1 = self.data[itime, ie, :]
                    t2 = table.data[itime, ie, :]
                    (o11, o22, t12, t1z, t2z, angle, major, minor, ovm) = t1
                    (o112, o222, t122, t1z2, t2z2, angle2, major2, minor2, ovm2) = t2

                    # vm stress can be NaN for some reason...
                    if not np.array_equal(t1[:-1], t2[:-1]):
                        msg += (
                            '(%s, %s)    (%s, %s, %s, %s, %s, %s, %s, %s, %s)'
                            '  (%s, %s, %s, %s, %s, %s, %s, %s, %s)\n' % (
                                eid, layer,
                                o11, o22, t12, t1z, t2z, angle, major, minor, ovm,
                                o112, o222, t122, t1z2, t2z2, angle2, major2, minor2, ovm2))
                        i += 1
                        if i > 10:
                            print(msg)
                            raise ValueError(msg)
                #print(msg)
                if i > 0:
                    raise ValueError(msg)
        return True

    def add_new_eid_sort1(self, etype, dt, eid, layer, o11, o22, t12, t1z, t2z,
                          angle, major, minor, ovm):
        self._times[self.itime] = dt
        self.element_layer[self.itotal, :] = [eid, layer]
        self.data[self.itime, self.itotal, :] = [o11, o22, t12, t1z, t2z, angle, major, minor, ovm]
        self.itotal += 1
        self.ielement += 1

    def add_sort1(self, dt, eid, layer, o11, o22, t12, t1z, t2z, angle,
                  major, minor, ovm):
        """unvectorized method for adding SORT1 transient data"""
        assert eid is not None
        assert isinstance(eid, integer_types) and eid > 0, 'dt=%s eid=%s' % (dt, eid)
        self.element_layer[self.itotal, :] = [eid, layer]
        self.data[self.itime, self.itotal, :] = [o11, o22, t12, t1z, t2z, angle, major, minor, ovm]
        self.itotal += 1

    def get_stats(self, short=False):
        if not self.is_built:
            msg = [
                '<%s>\n' % self.__class__.__name__,
                '  ntimes: %i\n' % self.ntimes,
                '  ntotal: %i\n' % self.ntotal,
            ]
            return msg

        nelements = self.nelements
        ntimes = self.ntimes
        #nnodes = self.nnodes
        ntotal = self.ntotal
        nelements = len(unique(self.element_layer[:, 0]))

        msg = []
        if self.nonlinear_factor not in (None, np.nan):  # transient
            msg.append('  type=%s ntimes=%i nelements=%i ntotal=%i\n'
                       % (self.__class__.__name__, ntimes, nelements, ntotal))
            ntimes_word = 'ntimes'
        else:
            msg.append('  type=%s nelements=%i ntotal=%i\n'
                       % (self.__class__.__name__, nelements, ntotal))
            ntimes_word = '1'
        headers = self.get_headers()
        n = len(headers)
        msg.append('  data: [%s, ntotal, %i] where %i=[%s]\n' % (ntimes_word, n, n, str(', '.join(headers))))
        msg.append('  element_layer.shape = %s\n' % str(self.element_layer.shape).replace('L', ''))
        msg.append('  data.shape = %s\n' % str(self.data.shape).replace('L', ''))
        msg.append('  element type: %s\n' % self.element_name)
        msg += self.get_data_code()
        return msg

    def get_f06_header(self, is_mag_phase=True):
        ctria3_msg, ctria6_msg, cquad4_msg, cquad8_msg = self._get_msgs()

        if self.element_type == 95:  # CQUAD4
            msg = cquad4_msg
            nnodes = 4
        elif self.element_type == 96:  # CQUAD8
            msg = cquad8_msg
            nnodes = 8
        elif self.element_type == 97:  # CTRIA3
            msg = ctria3_msg
            nnodes = 3
        elif self.element_type == 98:  # CTRIA6
            msg = ctria6_msg
            nnodes = 6
        else:
            raise NotImplementedError(self.element_name)

        return self.element_name, nnodes, msg

    def get_element_index(self, eids):
        # elements are always sorted; nodes are not
        itot = searchsorted(eids, self.element_layer[:, 0])  #[0]
        return itot

    def eid_to_element_node_index(self, eids):
        ind = ravel([searchsorted(self.element_layer[:, 0] == eid) for eid in eids])
        #ind = searchsorted(eids, self.element)
        #ind = ind.reshape(ind.size)
        #ind.sort()
        return ind

    def write_f06(self, f06_file, header=None, page_stamp='PAGE %s',
                  page_num=1, is_mag_phase=False, is_sort1=True):
        if header is None:
            header = []
        #msg, nnodes, is_bilinear = self._get_msgs()
        if self.is_von_mises:
            von = 'VON'
            mises = 'MISES'
        else:
            von = 'MAX'
            mises = 'SHEAR'

        if self.is_strain:
            words = ['   ELEMENT  PLY   STRAINS IN FIBER AND MATRIX DIRECTIONS    INTER-LAMINAR   STRAINS  PRINCIPAL  STRAINS (ZERO SHEAR)      %s\n' % von,
                     '     ID      ID    NORMAL-1     NORMAL-2     SHEAR-12     SHEAR XZ-MAT  SHEAR YZ-MAT  ANGLE    MAJOR        MINOR        %s\n' % mises]
        else:
            words = ['   ELEMENT  PLY  STRESSES IN FIBER AND MATRIX DIRECTIONS    INTER-LAMINAR  STRESSES  PRINCIPAL STRESSES (ZERO SHEAR)      %s\n' % von,
                     '     ID      ID    NORMAL-1     NORMAL-2     SHEAR-12     SHEAR XZ-MAT  SHEAR YZ-MAT  ANGLE    MAJOR        MINOR        %s\n' % mises]

        if self.element_type == 95:  # CQUAD4
            if self.is_strain:
                msg = ['                     S T R A I N S   I N   L A Y E R E D   C O M P O S I T E   E L E M E N T S   ( Q U A D 4 )\n'] + words
            else:
                msg = ['                   S T R E S S E S   I N   L A Y E R E D   C O M P O S I T E   E L E M E N T S   ( Q U A D 4 )\n'] + words
        #elif self.element_type == 96:  # CQUAD8
            #nnodes_per_element = 1
        elif self.element_type == 97:  # CTRIA3
            if self.is_strain:
                msg = ['                     S T R A I N S   I N   L A Y E R E D   C O M P O S I T E   E L E M E N T S   ( T R I A 3 )\n'] + words
            else:
                msg = ['                   S T R E S S E S   I N   L A Y E R E D   C O M P O S I T E   E L E M E N T S   ( T R I A 3 )\n'] + words
        elif self.element_type == 96:  # QUAD8
            # good
            if self.is_strain:
                msg = ['                     S T R A I N S   I N   L A Y E R E D   C O M P O S I T E   E L E M E N T S   ( Q U A D 8 )\n'] + words
            else:
                msg = ['                   S T R E S S E S   I N   L A Y E R E D   C O M P O S I T E   E L E M E N T S   ( Q U A D 8 )\n'] + words

        elif self.element_type == 98:  # CTRIA6
            # good
            if self.is_strain:
                msg = ['                     S T R A I N S   I N   L A Y E R E D   C O M P O S I T E   E L E M E N T S   ( T R I A 6 )\n'] + words
            else:
                msg = ['                   S T R E S S E S   I N   L A Y E R E D   C O M P O S I T E   E L E M E N T S   ( T R I A 6 )\n'] + words
        else:
            msg = 'element_name=%s element_type=%s' % (self.element_name, self.element_type)
            raise NotImplementedError(msg)

        # write the f06
        ntimes = self.data.shape[0]

        eids = self.element_layer[:, 0]
        layers = self.element_layer[:, 1]

        for itime in range(ntimes):
            dt = self._times[itime]
            header = _eigenvalue_header(self, header, itime, ntimes, dt)
            f06_file.write(''.join(header + msg))

            #print("self.data.shape=%s itime=%s ieids=%s" % (str(self.data.shape), itime, str(ieids)))

            #[o11, o22, t12, t1z, t2z, angle, major, minor, ovm]
            o11 = self.data[itime, :, 0]
            o22 = self.data[itime, :, 1]
            t12 = self.data[itime, :, 2]
            t1z = self.data[itime, :, 3]
            t2z = self.data[itime, :, 4]
            angle = self.data[itime, :, 5]
            major = self.data[itime, :, 6]
            minor = self.data[itime, :, 7]
            ovm = self.data[itime, :, 8]

            for eid, layer, o11i, o22i, t12i, t1zi, t2zi, anglei, majori, minori, ovmi in zip(
                eids, layers, o11, o22, t12, t1z, t2z, angle, major, minor, ovm):

                [o11i, o22i, t12i, t1zi, t2zi, majori, minori, ovmi] = write_floats_12e([
                 o11i, o22i, t12i, t1zi, t2zi, majori, minori, ovmi])
                f06_file.write('0 %8s %4s  %12s %12s %12s   %12s %12s  %6.2F %12s %12s %s\n'
                               % (eid, layer, o11i, o22i, t12i, t1zi, t2zi, anglei, majori, minori, ovmi))
            f06_file.write(page_stamp % page_num)
            page_num += 1
        return page_num - 1

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

        #print("nnodes_all =", nnodes_all)
        #msg.append('  element_node.shape = %s\n' % str(self.element_node.shape).replace('L', ''))
        #msg.append('  data.shape=%s\n' % str(self.data.shape).replace('L', ''))

        eids = self.element_layer[:, 0]
        layers = self.element_layer[:, 1]
        eids_device = eids * 10 + self.device_code

        nelements = len(np.unique(eids))
        #print('nelements =', nelements)
        # 21 = 1 node, 3 principal, 6 components, 9 vectors, 2 p/ovm
        #ntotal = ((nnodes * 21) + 1) + (nelements * 4)

        ntotali = self.num_wide
        nlayers = self.data.shape[1]
        ntotal = ntotali * nlayers

        #print('shape = %s' % str(self.data.shape))
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

        ntimes = self.data.shape[0]

        for itime in range(ntimes):
            nwide = 0
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

            #dt = self._times[itime]
            #header = _eigenvalue_header(self, header, itime, ntimes, dt)
            #f06_file.write(''.join(header + msg))

            #[o11, o22, t12, t1z, t2z, angle, major, minor, ovm]
            o11 = self.data[itime, :, 0]
            o22 = self.data[itime, :, 1]
            t12 = self.data[itime, :, 2]
            t1z = self.data[itime, :, 3]
            t2z = self.data[itime, :, 4]
            angle = self.data[itime, :, 5]
            major = self.data[itime, :, 6]
            minor = self.data[itime, :, 7]
            ovm = self.data[itime, :, 8]

            for eid_device, eid, layer, o11i, o22i, t12i, t1zi, t2zi, anglei, majori, minori, ovmi in zip(
                eids_device, eids, layers, o11, o22, t12, t1z, t2z, angle, major, minor, ovm):

                data = [eid_device, layer, o11i, o22i, t12i, t1zi, t2zi, anglei, majori, minori, ovmi]
                op2.write(pack('2i 9f', *data))

                [o11i, o22i, t12i, t1zi, t2zi, majori, minori, ovmi] = write_floats_12e([
                 o11i, o22i, t12i, t1zi, t2zi, majori, minori, ovmi])
                op2_ascii.write('0 %8s %4s  %12s %12s %12s   %12s %12s  %6.2F %12s %12s %s\n'
                                % (eid, layer, o11i, o22i, t12i, t1zi, t2zi, anglei, majori, minori, ovmi))

                nwide += len(data)

            assert nwide == ntotal, "nwide=%s ntotal=%s" % (nwide, ntotal)
            itable -= 1
            header = [4 * ntotal,]
            op2.write(pack('i', *header))
            op2_ascii.write('footer = %s\n' % header)
            new_result = False
        return itable


class RealCompositePlateStressArray(RealCompositePlateArray, StressObject):
    def __init__(self, data_code, is_sort1, isubcase, dt):
        RealCompositePlateArray.__init__(self, data_code, is_sort1, isubcase, dt)
        StressObject.__init__(self, data_code, isubcase)

    @property
    def is_stress(self):
        return True

    @property
    def is_strain(self):
        return False

    def get_headers(self):
        if self.is_von_mises:
            ovm = 'von_mises'
        else:
            ovm = 'max_shear'
        headers = ['o11', 'o22', 't12', 't1z', 't2z', 'angle', 'major', 'minor', ovm]
        return headers


class RealCompositePlateStrainArray(RealCompositePlateArray, StrainObject):
    def __init__(self, data_code, is_sort1, isubcase, dt):
        RealCompositePlateArray.__init__(self, data_code, is_sort1, isubcase, dt)
        StrainObject.__init__(self, data_code, isubcase)

    @property
    def is_stress(self):
        return False

    @property
    def is_strain(self):
        return True

    def get_headers(self):
        if self.is_von_mises:
            ovm = 'von_mises'
        else:
            ovm = 'max_shear'
        headers = ['e11', 'e22', 'e12', 'e1z', 'e2z', 'angle', 'major', 'minor', ovm]
        return headers
