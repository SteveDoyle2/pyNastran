from itertools import cycle
from abc import abstractmethod
from typing import List

import numpy as np
from numpy import zeros, searchsorted, allclose

from pyNastran.utils.numpy_utils import integer_types
from pyNastran.op2.result_objects.op2_objects import BaseElement, get_complex_times_dtype
from pyNastran.op2.tables.oef_forces.oef_force_objects import ForceObject
from pyNastran.f06.f06_formatting import write_imag_floats_13e, write_float_12e # get_key0,
from pyNastran.f06.f06_formatting import _eigenvalue_header


class ComplexForceObject(ForceObject):
    def __init__(self, data_code, isubcase, apply_data_code=True):
        ForceObject.__init__(self, data_code, isubcase, apply_data_code=apply_data_code)

    @property
    def is_real(self):
        return False

    @property
    def is_complex(self):
        return True

    @property
    def nnodes_per_element(self):
        return 1


class ComplexRodForceArray(ComplexForceObject):
    def __init__(self, data_code, is_sort1, isubcase, dt):
        self.element_type = None
        self.element_name = None
        ComplexForceObject.__init__(self, data_code, isubcase)
        #self.code = [self.format_code, self.sort_code, self.s_code]

        #self.ntimes = 0  # or frequency/mode
        #self.ntotal = 0
        self.nelements = 0  # result specific

        if is_sort1:
            pass
        else:
            raise NotImplementedError('SORT2')

    def get_headers(self) -> List[str]:
        headers = ['axial_force', 'torque']
        return headers

    #def get_headers(self):
        #headers = ['axial', 'torque']
        #return headers

    def build(self):
        """sizes the vectorized attributes of the ComplexRodForceArray"""
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
        dtype, idtype, cfdtype = get_complex_times_dtype(self.nonlinear_factor, self.size)

        self._times = zeros(self.ntimes, dtype=dtype)
        self.element = zeros(self.nelements, dtype=idtype)

        #[axial_force, torque]
        self.data = zeros((self.ntimes, self.ntotal, 2), dtype=cfdtype)

    def build_dataframe(self):
        """creates a pandas dataframe"""
        headers = self.get_headers()
        column_names, column_values = self._build_dataframe_transient_header()
        data_frame = self._build_pandas_transient_elements(column_values, column_names,
                                                           headers, self.element, self.data)
        #data_frame = pd.Panel(self.data, items=column_values,
                              #major_axis=self.element, minor_axis=headers).to_frame()
        #data_frame.columns.names = column_names
        #data_frame.index.names = ['ElementID', 'Item']
        self.data_frame = data_frame

    def __eq__(self, table):  # pragma: no cover
        self._eq_header(table)
        assert self.is_sort1 == table.is_sort1
        if not np.array_equal(self.data, table.data):
            msg = 'table_name=%r class_name=%s\n' % (self.table_name, self.__class__.__name__)
            msg += '%s\n' % str(self.code_information())
            i = 0
            for itime in range(self.ntimes):
                for ie, eid in enumerate(self.element):
                    t1 = self.data[itime, ie, :]
                    t2 = table.data[itime, ie, :]
                    (axial1, torque1) = t1
                    (axial2, torque2) = t2

                    if not allclose(t1, t2):
                        msg += '(%s)    (%s, %s)  (%s, %s)\n' % (
                            eid,
                            axial1, torque1,
                            axial2, torque2)
                        i += 1
                        if i > 10:
                            print(msg)
                            raise ValueError(msg)
                #print(msg)
                if i > 0:
                    raise ValueError(msg)
        return True

    def add_sort1(self, dt, eid, axial, torque):
        """unvectorized method for adding SORT1 transient data"""
        assert isinstance(eid, integer_types) and eid > 0, 'dt=%s eid=%s' % (dt, eid)
        self._times[self.itime] = dt
        self.element[self.ielement] = eid
        self.data[self.itime, self.ielement, :] = [axial, torque]
        self.ielement += 1

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

        msg = []
        if self.nonlinear_factor not in (None, np.nan):  # transient
            msg.append('  type=%s ntimes=%i nelements=%i; table_name=%r\n'
                       % (self.__class__.__name__, ntimes, nelements, self.table_name))
            ntimes_word = 'ntimes'
        else:
            msg.append('  type=%s nelements=%i; table_name=%r\n'
                       % (self.__class__.__name__, nelements, self.table_name))
            ntimes_word = '1'
        msg.append('  eType\n')
        headers = self.get_headers()
        n = len(headers)
        msg.append('  data: [%s, nnodes, %i] where %i=[%s]\n' % (
            ntimes_word, n, n, str(', '.join(headers))))
        msg.append('  data.shape = %s\n' % str(self.data.shape).replace('L', ''))
        #msg.append('  element type: %s\n' % self.element_type)
        msg.append('  element name: %s\n' % self.element_name)
        msg += self.get_data_code()
        return msg

    def get_f06_header(self, is_mag_phase=True, is_sort1=True):
        if self.element_type == 1: # CROD
            msg = ['                             C O M P L E X   F O R C E S   I N   R O D   E L E M E N T S   ( C R O D )\n']
        elif self.element_type == 10: # CONROD
            msg = ['                           C O M P L E X   F O R C E S   I N   R O D   E L E M E N T S   ( C O N R O D )\n']
        elif self.element_type == 3: # CTUBE
            msg = ['                            C O M P L E X   F O R C E S   I N   R O D   E L E M E N T S   ( C T U B E )\n']
            #pass
        else:
            raise NotImplementedError('element_name=%s element_type=%s' % (self.element_name, self.element_type))

        if is_mag_phase:
            msg += ['                                                          (MAGNITUDE/PHASE)\n']
        else:
            msg += ['                                                          (REAL/IMAGINARY)\n']

        if is_sort1:
            msg += [
                ' \n'
                '                 ELEMENT                             AXIAL                                       TORSIONAL\n'
                '                   ID.                              STRAIN                                         STRAIN\n'
            ]
            #'                      14                  0.0          /  0.0                           0.0          /  0.0'
        else:
            raise NotImplementedError('sort2')


        return self.element_name, msg

    def get_element_index(self, eids):
        # elements are always sorted; nodes are not
        itot = searchsorted(eids, self.element)  #[0]
        return itot

    def eid_to_element_node_index(self, eids):
        #ind = ravel([searchsorted(self.element == eid) for eid in eids])
        ind = searchsorted(eids, self.element)
        #ind = ind.reshape(ind.size)
        #ind.sort()
        return ind

    def write_f06(self, f06_file, header=None, page_stamp='PAGE %s', page_num=1, is_mag_phase=False, is_sort1=True):
        if header is None:
            header = []
        (elem_name, msg_temp) = self.get_f06_header(is_mag_phase=is_mag_phase, is_sort1=is_sort1)

        # write the f06
        #(ntimes, ntotal, two) = self.data.shape
        ntimes = self.data.shape[0]

        eids = self.element
        #is_odd = False
        #nwrite = len(eids)
        #if len(eids) % 2 == 1:
            #nwrite -= 1
            #is_odd = True

        #print('len(eids)=%s nwrite=%s is_odd=%s' % (len(eids), nwrite, is_odd))
        for itime in range(ntimes):
            dt = self._times[itime]  # TODO: rename this...
            header = _eigenvalue_header(self, header, itime, ntimes, dt)
            f06_file.write(''.join(header + msg_temp))

            #print("self.data.shape=%s itime=%s ieids=%s" % (str(self.data.shape), itime, str(ieids)))
            axial = self.data[itime, :, 0]
            torsion = self.data[itime, :, 1]

            for eid, axiali, torsioni in zip(eids, axial, torsion):
                out = write_imag_floats_13e([axiali, torsioni], is_mag_phase)
                [raxial, rtorsion, iaxial, itorsion] = out
                #ELEMENT                             AXIAL                                       TORSIONAL
                    #ID.                              STRESS                                         STRESS
                    #14                  0.0          /  0.0                           0.0          /  0.0

                f06_file.write('      %8i   %-13s / %-13s   %-13s / %s\n' % (eid, raxial, iaxial, rtorsion, itorsion))
            f06_file.write(page_stamp % page_num)
            page_num += 1
        return page_num - 1

    def write_op2(self, op2, op2_ascii, itable, new_result, date,
                  is_mag_phase=False, endian='>'):
        """writes an OP2"""
        import inspect
        from struct import Struct, pack
        frame = inspect.currentframe()
        call_frame = inspect.getouterframes(frame, 2)
        op2_ascii.write('%s.write_op2: %s\n' % (self.__class__.__name__, call_frame[1][3]))

        if itable == -1:
            self._write_table_header(op2, op2_ascii, date)
            itable = -3

        #eids = self.element

        # table 4 info
        #ntimes = self.data.shape[0]
        #nnodes = self.data.shape[1]
        nelements = self.data.shape[1]

        # 21 = 1 node, 3 principal, 6 components, 9 vectors, 2 p/ovm
        #ntotal = ((nnodes * 21) + 1) + (nelements * 4)

        ntotali = self.num_wide
        ntotal = ntotali * nelements

        #device_code = self.device_code
        op2_ascii.write('  ntimes = %s\n' % self.ntimes)

        eids_device = self.element * 10 + self.device_code

        if self.is_sort1:
            struct1 = Struct(endian + b'i4f')
        else:
            raise NotImplementedError('SORT2')

        op2_ascii.write('nelements=%i\n' % nelements)

        for itime in range(self.ntimes):
            self._write_table_3(op2, op2_ascii, new_result, itable, itime)

            # record 4
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

            axial = self.data[itime, :, 0]
            torsion = self.data[itime, :, 1]

            for eid_device, axiali, torsioni in zip(eids_device, axial, torsion):
                data = [eid_device, axiali.real, torsioni.real, axiali.imag, torsioni.imag]
                op2_ascii.write('  eid_device=%s data=%s\n' % (eid_device, tuple(data)))
                op2.write(struct1.pack(*data))

            itable -= 1
            header = [4 * ntotal,]
            op2.write(pack('i', *header))
            op2_ascii.write('footer = %s\n' % header)
            new_result = False
        return itable


class ComplexCShearForceArray(BaseElement):
    def __init__(self, data_code, is_sort1, isubcase, dt):
        self.element_type = None
        self.element_name = None
        BaseElement.__init__(self, data_code, isubcase)
        #self.code = [self.format_code, self.sort_code, self.s_code]

        #self.ntimes = 0  # or frequency/mode
        #self.ntotal = 0
        self.nelements = 0  # result specific

        #if is_sort1:
            #pass
        #else:
            #raise NotImplementedError('SORT2')

    @property
    def is_real(self) -> bool:
        return False

    @property
    def is_complex(self) -> bool:
        return True

    @property
    def nnodes_per_element(self) -> int:
        return 1

    def _reset_indices(self):
        self.itotal = 0
        self.ielement = 0

    def get_headers(self) -> List[str]:
        headers = [
            'force41', 'force14', 'force21', 'force12', 'force32', 'force23',
            'force43', 'force34', 'kickForce1', 'kickForce2', 'kickForce3',
            'kickForce4', 'shear12', 'shear23', 'shear34', 'shear41'
        ]
        return headers

    #def get_headers(self):
        #headers = ['axial', 'torque']
        #return headers

    def build(self):
        """sizes the vectorized attributes of the ComplexCShearForceArray"""
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
        dtype, idtype, cfdtype = get_complex_times_dtype(self.nonlinear_factor, self.size)
        self._times = zeros(self.ntimes, dtype=dtype)
        self.element = zeros(self.nelements, dtype=idtype)

        #[force41, force14, force21, force12, force32, force23, force43, force34,
        #kick_force1, kick_force2, kick_force3, kick_force4,
        #shear12, shear23, shear34, shear41]
        self.data = zeros((self.ntimes, self.ntotal, 16), dtype=cfdtype)

    def build_dataframe(self):
        """creates a pandas dataframe"""
        #Mode                                           1                   2
        #EigenvalueReal                              -0.0                -0.0
        #EigenvalueImag                              -0.0                -0.0
        #Damping                                      0.0                 0.0
        #ElementID Item
        #22        force41     2.927977e-10+0.000000e+00j  0.000000+0.000000j
        #          force14     2.927977e-10+5.855954e-10j  0.000000+0.000000j
        #          force21    -2.927977e-10+0.000000e+00j  0.000000+0.000000j
        #          force12    -2.927977e-10+5.855954e-10j  0.000000+0.000000j
        #          force32     2.927977e-10+0.000000e+00j  0.000000+0.000000j
        #          force23     2.927977e-10+5.855954e-10j  0.000000+0.000000j
        #          force43    -2.927977e-10+0.000000e+00j  0.000000+0.000000j
        #          force34    -2.927977e-10+5.855954e-10j  0.000000+0.000000j
        #          kickForce1  0.000000e+00+0.000000e+00j  0.000000+0.000000j
        #          kickForce2  0.000000e+00+0.000000e+00j  0.000000+0.000000j
        #          kickForce3  0.000000e+00+0.000000e+00j  0.000000+0.000000j
        #          kickForce4  0.000000e+00+0.000000e+00j  0.000000+0.000000j
        #          shear12     0.000000e+00+0.000000e+00j  0.000000+0.000000j
        #          shear23     0.000000e+00+0.000000e+00j  0.000000+0.000000j
        #          shear34     0.000000e+00+0.000000e+00j  0.000000+0.000000j
        #          shear41     0.000000e+00+0.000000e+00j  0.000000+0.000000j
        headers = self.get_headers()
        column_names, column_values = self._build_dataframe_transient_header()
        self.data_frame = self._build_pandas_transient_elements(
            column_values, column_names,
            headers, self.element, self.data)

    def __eq__(self, table):  # pragma: no cover
        self._eq_header(table)
        assert self.is_sort1 == table.is_sort1
        if not np.array_equal(self.data, table.data):
            msg = 'table_name=%r class_name=%s\n' % (self.table_name, self.__class__.__name__)
            msg += '%s\n' % str(self.code_information())
            i = 0
            for itime in range(self.ntimes):
                for ie, eid in enumerate(self.element):
                    t1 = self.data[itime, ie, :]
                    t2 = table.data[itime, ie, :]
                    (force41a, force14a, force21a, force12a, force32a, force23a, force43a, force34a,
                     kick_force1a, kick_force2a, kick_force3a, kick_force4a,
                     shear12a, shear23a, shear34a, shear41a) = t1

                    (force41b, force14b, force21b, force12b, force32b, force23b, force43b, force34b,
                     kick_force1b, kick_force2b, kick_force3b, kick_force4b,
                     shear12b, shear23b, shear34b, shear41b) = t2

                    if not allclose(t1, t2):
                        msg += (
                            '%s   (%s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s)\n'
                            '     (%s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s)\n' % (
                                eid,
                                force41a, force14a, force21a, force12a, force32a, force23a,
                                force43a, force34a, kick_force1a, kick_force2a, kick_force3a,
                                kick_force4a, shear12a, shear23a, shear34a, shear41a,

                                force41b, force14b, force21b, force12b, force32b, force23b,
                                force43b, force34b, kick_force1b, kick_force2b, kick_force3b,
                                kick_force4b, shear12b, shear23b, shear34b, shear41b
                            ))
                        i += 1
                        if i > 10:
                            print(msg)
                            raise ValueError(msg)
                #print(msg)
                if i > 0:
                    raise ValueError(msg)
        return True

    def add_sort1(self, dt, eid,
                  force41, force14, force21, force12, force32, force23, force43, force34,
                  kick_force1, kick_force2, kick_force3, kick_force4,
                  shear12, shear23, shear34, shear41):
        """unvectorized method for adding SORT1 transient data"""
        assert isinstance(eid, integer_types) and eid > 0, 'dt=%s eid=%s' % (dt, eid)
        self._times[self.itime] = dt
        self.element[self.ielement] = eid
        self.data[self.itime, self.ielement, :] = [
            force41, force14, force21, force12, force32, force23, force43, force34,
            kick_force1, kick_force2, kick_force3, kick_force4,
            shear12, shear23, shear34, shear41]
        self.ielement += 1

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

        msg = []
        if self.nonlinear_factor not in (None, np.nan):  # transient
            msg.append('  type=%s ntimes=%i nelements=%i; table_name=%r\n'
                       % (self.__class__.__name__, ntimes, nelements, self.table_name))
            ntimes_word = 'ntimes'
        else:
            msg.append('  type=%s nelements=%i; table_name=%r\n'
                       % (self.__class__.__name__, nelements, self.table_name))
            ntimes_word = '1'
        msg.append('  eType\n')
        headers = self.get_headers()
        n = len(headers)
        msg.append('  data: [%s, nnodes, %i] where %i=[%s]\n' % (
            ntimes_word, n, n, str(', '.join(headers))))
        msg.append('  data.shape = %s\n' % str(self.data.shape).replace('L', ''))
        #msg.append('  element type: %s\n' % self.element_type)
        msg.append('  element name: %s\n' % self.element_name)
        msg += self.get_data_code()
        return msg

    def get_f06_header(self, is_mag_phase=True, is_sort1=True):
        msg = ['                C O M P L E X   F O R C E S   A C T I N G   O N   S H E A R   P A N E L   E L E M E N T S   (CSHEAR)\n']

        if is_mag_phase:
            msg += ['                                                          (MAGNITUDE/PHASE)\n \n']
        else:
            msg += ['                                                          (REAL/IMAGINARY)\n \n']

        if is_sort1:
            msg += [
                '                  ====== POINT  1 ======      ====== POINT  2 ======      ====== POINT  3 ======      ====== POINT  4 ======\n'
                ' ELEMENT          F-FROM-4      F-FROM-2      F-FROM-1      F-FROM-3      F-FROM-2      F-FROM-4      F-FROM-3      F-FROM-1\n'
                '         ID               KICK-1       SHEAR-12       KICK-2       SHEAR-23       KICK-3       SHEAR-34       KICK-4       SHEAR-41\n'
            ]
        else:
            raise NotImplementedError('sort2')
        return msg

    def write_f06(self, f06_file, header=None, page_stamp='PAGE %s',
                  page_num=1, is_mag_phase=False, is_sort1=True):
        if header is None:
            header = []
        msg_temp = self.get_f06_header(is_mag_phase=is_mag_phase, is_sort1=is_sort1)

        # write the f06
        #(ntimes, ntotal, two) = self.data.shape
        ntimes = self.data.shape[0]

        eids = self.element
        for itime in range(ntimes):
            dt = self._times[itime]  # TODO: rename this...
            header = _eigenvalue_header(self, header, itime, ntimes, dt)
            f06_file.write(''.join(header + msg_temp))

            #print("self.data.shape=%s itime=%s ieids=%s" % (str(self.data.shape), itime, str(ieids)))

            ## TODO: I'm sure this ordering is wrong...
            force41 = self.data[itime, :, 0]
            force14 = self.data[itime, :, 1]
            force21 = self.data[itime, :, 2]  # TODO: this is wrong...
            force12 = self.data[itime, :, 3]
            force32 = self.data[itime, :, 4]
            force23 = self.data[itime, :, 5]
            force43 = self.data[itime, :, 6]
            force34 = self.data[itime, :, 7]
            kick_force1 = self.data[itime, :, 8]
            kick_force2 = self.data[itime, :, 9]
            kick_force3 = self.data[itime, :, 10]
            kick_force4 = self.data[itime, :, 11]
            shear12 = self.data[itime, :, 12]
            shear23 = self.data[itime, :, 13]
            shear34 = self.data[itime, :, 14]
            shear41 = self.data[itime, :, 15]
            assert len(force12) > 0, force12

            for (eid, iforce41, force14i, iforce21, iforce12, iforce32, iforce23, iforce43, iforce34,
                 ikick_force1, ikick_force2, ikick_force3, ikick_force4,
                 ishear12, ishear23, ishear34, ishear41) in zip(
                     eids, force41, force14, force21, force12, force32, force23, force43, force34,
                     kick_force1, kick_force2, kick_force3, kick_force4,
                     shear12, shear23, shear34, shear41):

                vals2 = write_imag_floats_13e([
                    iforce41, force14i, iforce21, iforce12, iforce32, iforce23, iforce43, iforce34,
                    ikick_force1, ikick_force2, ikick_force3, ikick_force4,
                    ishear12, ishear23, ishear34, ishear41], is_mag_phase)

                [
                    force41r, force14r, force21i, force12r, force32r, force23r, force43r, force34r,
                    kick_force1r, kick_force2r, kick_force3r, kick_force4r,
                    shear12r, shear23r, shear34r, shear41r,

                    force41i, force14i, force21i, force12i, force32i, force23i, force43i, force34i,
                    kick_force1i, kick_force2i, kick_force3i, kick_force4i,
                    shear12i, shear23i, shear34i, shear41i
                ] = vals2

                #complex_cshear_force_f06
                #'                  ====== POINT  1 ======      ====== POINT  2 ======      ====== POINT  3 ======      ====== POINT  4 ======'
                #' ELEMENT          F-FROM-4      F-FROM-2      F-FROM-1      F-FROM-3      F-FROM-2      F-FROM-4      F-FROM-3      F-FROM-1'
                #'         ID               KICK-1       SHEAR-12       KICK-2       SHEAR-23       KICK-3       SHEAR-34       KICK-4       SHEAR-41'
                #'            25  0.0           0.0           0.0           0.0           0.0           0.0           0.0           0.0'
                #'                0.0           0.0           0.0           0.0           0.0           0.0           0.0           0.0'

                f06_file.write(
                    '      %8i %-13s %-13s %-13s %-13s %-13s %-13s %-13s  %s\n'
                    '               %-13s %-13s %-13s %-13s %-13s %-13s %-13s  %s\n'
                    '                      %-13s %-13s %-13s %-13s %-13s %-13s %-13s  %s\n'
                    '                      %-13s %-13s %-13s %-13s %-13s %-13s %-13s  %s\n'% (
                        eid,
                        force41r, force14r, force21i, force12r, force32r, force23r, force43r, force34r,
                        kick_force1r, kick_force2r, kick_force3r, kick_force4r,
                        shear12r, shear23r, shear34r, shear41r,

                        force41i, force14i, force21i, force12i, force32i, force23i, force43i, force34i,
                        kick_force1i, kick_force2i, kick_force3i, kick_force4i,
                        shear12i, shear23i, shear34i, shear41i
                ))
            f06_file.write(page_stamp % page_num)
            page_num += 1
        return page_num - 1


class ComplexSpringDamperForceArray(ComplexForceObject):
    def __init__(self, data_code, is_sort1, isubcase, dt):
        self.element_type = None
        self.element_name = None
        ComplexForceObject.__init__(self, data_code, isubcase)
        #self.code = [self.format_code, self.sort_code, self.s_code]

        #self.ntimes = 0  # or frequency/mode
        #self.ntotal = 0
        self.nelements = 0  # result specific

        #if is_sort1:
            #pass
        #else:
            #raise NotImplementedError('SORT2')

    def get_headers(self) -> List[str]:
        headers = ['spring_force']
        return headers

    #def get_headers(self):
        #headers = ['axial', 'torque']
        #return headers

    def build(self):
        """sizes the vectorized attributes of the ComplexSpringDamperForceArray"""
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
        dtype, idtype, cfdtype = get_complex_times_dtype(self.nonlinear_factor, self.size)
        self._times = zeros(self.ntimes, dtype=dtype)
        self.element = zeros(self.nelements, dtype=idtype)

        #[axial_force, torque]
        self.data = zeros((self.ntimes, self.ntotal, 1), dtype=cfdtype)

    def build_dataframe(self):
        """creates a pandas dataframe"""
        #Mode                                     1                   2
        #EigenvalueReal                        -0.0                -0.0
        #EigenvalueImag                        -0.0                -0.0
        #Damping                                0.0                 0.0
        #ElementID Item
        #30        spring_force  0.000000+0.000000j  0.000000+0.000000j
        #31        spring_force  0.000000+0.000000j  0.000000+0.000000j
        #32        spring_force  0.000000+0.000000j  0.000000+0.000000j
        #33        spring_force  0.000000+0.000000j  0.000000+0.000000j
        headers = self.get_headers()
        column_names, column_values = self._build_dataframe_transient_header()
        self.data_frame = self._build_pandas_transient_elements(
            column_values, column_names,
            headers, self.element, self.data)


    def __eq__(self, table):  # pragma: no cover
        self._eq_header(table)
        assert self.is_sort1 == table.is_sort1
        if not np.array_equal(self.data, table.data):
            msg = 'table_name=%r class_name=%s\n' % (self.table_name, self.__class__.__name__)
            msg += '%s\n' % str(self.code_information())
            i = 0
            for itime in range(self.ntimes):
                for ie, eid in enumerate(self.element):
                    t1 = self.data[itime, ie, 0]
                    t2 = table.data[itime, ie, 0]

                    if not allclose([t1.real, t1.imag], [t2.real, t2.imag], atol=0.0001):
                        msg += '%s    (%s, %s)  (%s, %s)\n' % (
                            eid,
                            t1.real, t1.imag,
                            t2.real, t2.imag)
                        i += 1
                        if i > 10:
                            print(msg)
                            raise ValueError(msg)
                #print(msg)
                if i > 0:
                    raise ValueError(msg)
        return True

    def add_sort1(self, dt, eid, force):
        """unvectorized method for adding SORT1 transient data"""
        assert isinstance(eid, integer_types) and eid > 0, 'dt=%s eid=%s' % (dt, eid)
        self._times[self.itime] = dt
        self.element[self.ielement] = eid
        self.data[self.itime, self.ielement, 0] = force
        self.ielement += 1

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

        msg = []
        if self.nonlinear_factor not in (None, np.nan):  # transient
            msg.append('  type=%s ntimes=%i nelements=%i; table_name=%r\n'
                       % (self.__class__.__name__, ntimes, nelements, self.table_name))
            ntimes_word = 'ntimes'
        else:
            msg.append('  type=%s nelements=%i; table_name=%r\n'
                       % (self.__class__.__name__, nelements, self.table_name))
            ntimes_word = '1'
        msg.append('  eType\n')
        headers = self.get_headers()
        n = len(headers)
        msg.append('  data: [%s, nelements, %i] where %i=[%s]\n' % (ntimes_word, n, n, str(', '.join(headers))))
        msg.append('  data.shape = %s\n' % str(self.data.shape).replace('L', ''))
        #msg.append('  element type: %s\n' % self.element_type)
        msg.append('  element name: %s\n' % self.element_name)
        msg += self.get_data_code()
        return msg

    def get_f06_header(self, is_mag_phase=True, is_sort1=True):
        # 11-CELAS1, 12-CELAS2, 13-CELAS3, 14-CELAS4
        if self.element_type == 11:
            msg = ['                         C O M P L E X   F O R C E S   I N   S C A L A R   S P R I N G S   ( C E L A S 1 )\n']
        elif self.element_type == 12:
            msg = ['                         C O M P L E X   F O R C E S   I N   S C A L A R   S P R I N G S   ( C E L A S 2 )\n']
        elif self.element_type == 13:
            msg = ['                         C O M P L E X   F O R C E S   I N   S C A L A R   S P R I N G S   ( C E L A S 3 )\n']
        elif self.element_type == 14:
            msg = ['                         C O M P L E X   F O R C E S   I N   S C A L A R   S P R I N G S   ( C E L A S 4 )\n']
        elif self.element_type == 20: # CDAMP1
            msg = ['                         C O M P L E X   F O R C E S   I N   S C A L A R   D A M P E R S   ( C D A M P 1 )\n']
        elif self.element_type == 21: # CDAMP2
            msg = ['                         C O M P L E X   F O R C E S   I N   S C A L A R   D A M P E R S   ( C D A M P 2 )\n']
        elif self.element_type == 22: # CDAMP3
            msg = ['                         C O M P L E X   F O R C E S   I N   S C A L A R   D A M P E R S   ( C D A M P 3 )\n']
        elif self.element_type == 23: # CDAMP4
            msg = ['                         C O M P L E X   F O R C E S   I N   S C A L A R   D A M P E R S   ( C D A M P 4 )\n']
        else:
            raise NotImplementedError('element_name=%s element_type=%s' % (self.element_name, self.element_type))

        if is_mag_phase:
            msg += ['                                                          (MAGNITUDE/PHASE)\n \n']
        else:
            msg += ['                                                          (REAL/IMAGINARY)\n \n']

        if is_sort1:
            msg += [
                '                ELEMENT                                                   ELEMENT\n'
                '                  ID.                    FORCE                              ID.                    FORCE\n'
            ]
            #'                      14                  0.0          /  0.0                           0.0          /  0.0'
        else:
            msg += ['            FREQUENCY                    FORCE                        FREQUENCY                    FORCE\n']


        return msg

    #def get_element_index(self, eids):
        ## elements are always sorted; nodes are not
        #itot = searchsorted(eids, self.element)  #[0]
        #return itot

    #def eid_to_element_node_index(self, eids):
        ##ind = ravel([searchsorted(self.element == eid) for eid in eids])
        #ind = searchsorted(eids, self.element)
        ##ind = ind.reshape(ind.size)
        ##ind.sort()
        #return ind

    def write_f06(self, f06_file, header=None, page_stamp='PAGE %s', page_num=1, is_mag_phase=False, is_sort1=True):
        if header is None:
            header = []
        msg_temp = self.get_f06_header(is_mag_phase=is_mag_phase, is_sort1=is_sort1)

        # write the f06
        #(ntimes, ntotal, two) = self.data.shape
        ntimes = self.data.shape[0]

        eids = self.element
        #is_odd = False
        #nwrite = len(eids)
        #if len(eids) % 2 == 1:
            #nwrite -= 1
            #is_odd = True

        #print('len(eids)=%s nwrite=%s is_odd=%s' % (len(eids), nwrite, is_odd))
        for itime in range(ntimes):
            dt = self._times[itime]  # TODO: rename this...
            header = _eigenvalue_header(self, header, itime, ntimes, dt)
            f06_file.write(''.join(header + msg_temp))

            #print("self.data.shape=%s itime=%s ieids=%s" % (str(self.data.shape), itime, str(ieids)))
            spring_force = self.data[itime, :, 0]

            for eid, spring_forcei in zip(eids, spring_force):
                [rspring, ispring] = write_imag_floats_13e([spring_forcei], is_mag_phase)
                #ELEMENT                             AXIAL                                       TORSIONAL
                    #ID.                              STRESS                                         STRESS
                    #14                  0.0          /  0.0                           0.0          /  0.0

                f06_file.write('      %8i   %-13s / %-13s\n' % (eid, rspring, ispring))
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

        #eids = self.element

        # table 4 info
        #ntimes = self.data.shape[0]
        #nnodes = self.data.shape[1]
        nelements = self.data.shape[1]

        # 21 = 1 node, 3 principal, 6 components, 9 vectors, 2 p/ovm
        #ntotal = ((nnodes * 21) + 1) + (nelements * 4)

        ntotali = self.num_wide
        ntotal = ntotali * nelements

        #print('shape = %s' % str(self.data.shape))
        #assert self.ntimes == 1, self.ntimes

        #device_code = self.device_code
        op2_ascii.write('  ntimes = %s\n' % self.ntimes)

        eids_device = self.element * 10 + self.device_code

        #print('ntotal=%s' % (ntotal))
        #assert ntotal == 193, ntotal

        if self.is_sort1:
            struct1 = Struct(endian + b'i2f')
        else:
            raise NotImplementedError('SORT2')

        op2_ascii.write('%s-nelements=%i\n' % (self.element_name, nelements))
        for itime in range(self.ntimes):
            self._write_table_3(op2, op2_ascii, new_result, itable, itime)

            # record 4
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

            force = self.data[itime, :, 0]

            for eid, forcei in zip(eids_device, force):
                data = [eid, forcei.real, forcei.imag]
                op2_ascii.write('  eid=%s force=%s\n' % (eid, forcei))
                op2.write(struct1.pack(*data))

            itable -= 1
            header = [4 * ntotal,]
            op2.write(pack('i', *header))
            op2_ascii.write('footer = %s\n' % header)
            new_result = False
        return itable


class ComplexSpringForceArray(ComplexSpringDamperForceArray):  # 11-CELAS1,12-CELAS2,13-CELAS3, 14-CELAS4
    def __init__(self, data_code, is_sort1, isubcase, dt):
        ComplexSpringDamperForceArray.__init__(self, data_code, is_sort1, isubcase, dt)

class ComplexDamperForceArray(ComplexSpringDamperForceArray):
    def __init__(self, data_code, is_sort1, isubcase, dt):
        ComplexSpringDamperForceArray.__init__(self, data_code, is_sort1, isubcase, dt)


class ComplexViscForceArray(BaseElement):
    def __init__(self, data_code, is_sort1, isubcase, dt):
        self.element_type = None
        self.element_name = None
        BaseElement.__init__(self, data_code, isubcase)
        #self.code = [self.format_code, self.sort_code, self.s_code]

        #self.ntimes = 0  # or frequency/mode
        #self.ntotal = 0
        self.nelements = 0  # result specific

        if is_sort1:
            pass
        else:
            raise NotImplementedError('SORT2')

    @property
    def is_real(self):
        return False

    @property
    def is_complex(self):
        return True

    def _reset_indices(self):
        self.itotal = 0
        self.ielement = 0

    def get_headers(self) -> List[str]:
        headers = ['axial_force', 'torque']
        return headers

    #def get_headers(self):
        #headers = ['axial', 'torque']
        #return headers

    def build(self):
        """sizes the vectorized attributes of the ComplexViscForceArray"""
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
        dtype = 'float32'
        if isinstance(self.nonlinear_factor, integer_types):
            dtype = 'int32'
        self._times = zeros(self.ntimes, dtype=dtype)
        self.element = zeros(self.nelements, dtype='int32')

        #[axial_force, torque]
        self.data = zeros((self.ntimes, self.ntotal, 2), dtype='complex64')

    def build_dataframe(self):
        """creates a pandas dataframe"""
        #Mode                         1        2        3        4
        #EigenvalueReal            -0.0     -0.0     -0.0     -0.0
        #EigenvalueImag            -0.0     -0.0     -0.0     -0.0
        #Damping                    0.0      0.0      0.0      0.0
        #ElementID Item
        #50        axial_force  (-0+0j)  (-0+0j)  (-0+0j)  (-0+0j)
        #          torque       (-0+0j)  (-0+0j)  (-0+0j)  (-0+0j)
        #51        axial_force  (-0+0j)  (-0+0j)  (-0+0j)  (-0+0j)
        #          torque            0j  (-0+0j)  (-0+0j)  (-0+0j)
        headers = self.get_headers()
        column_names, column_values = self._build_dataframe_transient_header()
        data_frame = self._build_pandas_transient_elements(column_values, column_names,
                                                           headers, self.element, self.data)
        self.data_frame = data_frame

    def __eq__(self, table):  # pragma: no cover
        self._eq_header(table)
        assert self.is_sort1 == table.is_sort1
        if not np.array_equal(self.data, table.data):
            msg = 'table_name=%r class_name=%s\n' % (self.table_name, self.__class__.__name__)
            msg += '%s\n' % str(self.code_information())
            i = 0
            for itime in range(self.ntimes):
                for ie, eid in enumerate(self.element):
                    t1 = self.data[itime, ie, :]
                    t2 = table.data[itime, ie, :]
                    (axial1, torque1) = t1
                    (axial2, torque2) = t2

                    if not allclose(t1, t2):
                        msg += '(%s)    (%s, %s)  (%s, %s)\n' % (
                            eid,
                            axial1, torque1,
                            axial2, torque2)
                        i += 1
                        if i > 10:
                            print(msg)
                            raise ValueError(msg)
                #print(msg)
                if i > 0:
                    raise ValueError(msg)
        return True

    def add_sort1(self, dt, eid, axial, torque):
        """unvectorized method for adding SORT1 transient data"""
        assert isinstance(eid, integer_types) and eid > 0, 'dt=%s eid=%s' % (dt, eid)
        self._times[self.itime] = dt
        self.element[self.ielement] = eid
        self.data[self.itime, self.ielement, :] = [axial, torque]
        self.ielement += 1

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

        msg = []
        if self.nonlinear_factor not in (None, np.nan):  # transient
            msg.append('  type=%s ntimes=%i nelements=%i; table_name=%r\n'
                       % (self.__class__.__name__, ntimes, nelements, self.table_name))
            ntimes_word = 'ntimes'
        else:
            msg.append('  type=%s nelements=%i; table_name=%r\n'
                       % (self.__class__.__name__, nelements, self.table_name))
            ntimes_word = '1'
        msg.append('  eType\n')
        headers = self.get_headers()
        n = len(headers)
        msg.append('  data: [%s, nnodes, %i] where %i=[%s]\n' % (ntimes_word, n, n, str(', '.join(headers))))
        msg.append('  data.shape = %s\n' % str(self.data.shape).replace('L', ''))
        #msg.append('  element type: %s\n' % self.element_type)
        msg.append('  element name: %s\n' % self.element_name)
        msg += self.get_data_code()
        return msg

    def get_f06_header(self, is_mag_phase=True, is_sort1=True):
        #if self.element_type == 1: # CROD
            #msg = ['                             C O M P L E X   F O R C E S   I N   R O D   E L E M E N T S   ( C R O D )\n']
        #elif self.element_type == 10: # CONROD
            #msg = ['                           C O M P L E X   F O R C E S   I N   R O D   E L E M E N T S   ( C O N R O D )\n']
        #elif self.element_type == 3: # CTUBE
            #msg = ['                            C O M P L E X   F O R C E S   I N   R O D   E L E M E N T S   ( C T U B E )\n']
            ##pass
        if self.element_type == 24:
            msg = ['                           C O M P L E X   F O R C E S   I N   V I S C   E L E M E N T S   ( C V I S C )\n']
        else:
            raise NotImplementedError('element_name=%s element_type=%s' % (self.element_name, self.element_type))
        if is_mag_phase:
            msg += ['                                                          (MAGNITUDE/PHASE)\n']
        else:
            msg += ['                                                          (REAL/IMAGINARY)\n']

        if is_sort1:
            msg += [
                ' \n'
                '                 ELEMENT                             AXIAL                                       TORSIONAL\n'
                '                   ID.                              STRAIN                                         STRAIN\n'
            ]
            #'                      14                  0.0          /  0.0                           0.0          /  0.0'
        else:
            raise NotImplementedError('sort2')


        return self.element_name, msg

    def get_element_index(self, eids):
        # elements are always sorted; nodes are not
        itot = searchsorted(eids, self.element)  #[0]
        return itot

    def eid_to_element_node_index(self, eids):
        #ind = ravel([searchsorted(self.element == eid) for eid in eids])
        ind = searchsorted(eids, self.element)
        #ind = ind.reshape(ind.size)
        #ind.sort()
        return ind

    def write_f06(self, f06_file, header=None, page_stamp='PAGE %s',
                  page_num=1, is_mag_phase=False, is_sort1=True):
        if header is None:
            header = []
        (elem_name, msg_temp) = self.get_f06_header(is_mag_phase=is_mag_phase, is_sort1=is_sort1)

        # write the f06
        #(ntimes, ntotal, two) = self.data.shape
        ntimes = self.data.shape[0]

        eids = self.element
        #is_odd = False
        #nwrite = len(eids)
        #if len(eids) % 2 == 1:
            #nwrite -= 1
            #is_odd = True

        #print('len(eids)=%s nwrite=%s is_odd=%s' % (len(eids), nwrite, is_odd))
        for itime in range(ntimes):
            dt = self._times[itime]  # TODO: rename this...
            header = _eigenvalue_header(self, header, itime, ntimes, dt)
            f06_file.write(''.join(header + msg_temp))

            #print("self.data.shape=%s itime=%s ieids=%s" % (str(self.data.shape), itime, str(ieids)))
            axial = self.data[itime, :, 0]
            torsion = self.data[itime, :, 1]

            for eid, axiali, torsioni in zip(eids, axial, torsion):
                out = write_imag_floats_13e([axiali, torsioni], is_mag_phase)
                [raxial, rtorsion, iaxial, itorsion] = out
                #ELEMENT                             AXIAL                                       TORSIONAL
                    #ID.                              STRESS                                         STRESS
                    #14                  0.0          /  0.0                           0.0          /  0.0
                f06_file.write('      %8i   %-13s / %-13s   %-13s / %s\n' %
                               (eid, raxial, iaxial, rtorsion, itorsion))
            f06_file.write(page_stamp % page_num)
            page_num += 1
        return page_num - 1


class ComplexPlateForceArray(ComplexForceObject):
    def __init__(self, data_code, is_sort1, isubcase, dt):
        self.element_type = None
        self.element_name = None
        ComplexForceObject.__init__(self, data_code, isubcase)
        #self.code = [self.format_code, self.sort_code, self.s_code]

        #self.ntimes = 0  # or frequency/mode
        #self.ntotal = 0
        self.nelements = 0  # result specific

        #if is_sort1:
            #pass
        #else:
            #raise NotImplementedError('SORT2')

    def get_headers(self) -> List[str]:
        headers = ['mx', 'my', 'mxy', 'bmx', 'bmy', 'bmxy', 'tx', 'ty']
        return headers

    def build(self):
        """sizes the vectorized attributes of the ComplexPlateForceArray"""
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
        dtype, idtype, cfdtype = get_complex_times_dtype(self.nonlinear_factor, self.size)
        self._times = zeros(self.ntimes, dtype=dtype)
        self.element = zeros(self.nelements, dtype=idtype)

        #[mx, my, mxy, bmx, bmy, bmxy, tx, ty]
        self.data = zeros((self.ntimes, self.ntotal, 8), dtype=cfdtype)

    def build_dataframe(self):
        """creates a pandas dataframe"""
        # Freq           0.00001  10.00000 20.00000 30.00000                 40.00000 50.00000 60.00000
        # ElementID Item
        #8         mx         0j       0j       0j       0j   (-361.6303-680.04156j)       0j       0j
        #          my         0j       0j       0j       0j  (-7884.6196-14826.936j)       0j       0j
        #          mxy        0j       0j       0j       0j    (-237.5723-446.7519j)       0j       0j
        #          bmx        0j       0j       0j       0j   (5.514431+10.3698225j)       0j       0j
        #          bmy        0j       0j       0j       0j    (10.107019+19.00613j)       0j       0j
        #          bmxy       0j       0j       0j       0j  (-16.361727-30.768036j)       0j       0j
        #          tx         0j       0j       0j       0j     (18.819313+35.3895j)       0j       0j
        #          ty         0j       0j       0j       0j   (-61.55238-115.74853j)       0j       0j
        #9         mx         0j       0j       0j       0j   (1086.9078+2043.9175j)       0j       0j
        #          my         0j       0j       0j       0j    (8089.895+15212.953j)       0j       0j
        #          mxy        0j       0j       0j       0j   (-4725.3286-8885.925j)       0j       0j
        #          bmx        0j       0j       0j       0j   (-3.9810739-7.486363j)       0j       0j
        #          bmy        0j       0j       0j       0j  (-10.283798-19.338562j)       0j       0j
        #          bmxy       0j       0j       0j       0j   (-8.663734-16.292051j)       0j       0j
        #          tx         0j       0j       0j       0j    (54.14508+101.81919j)       0j       0j
        #          ty         0j       0j       0j       0j   (-61.92162-116.44288j)       0j       0j
        headers = self.get_headers()
        column_names, column_values = self._build_dataframe_transient_header()
        data_frame = self._build_pandas_transient_elements(column_values, column_names,
                                                           headers, self.element, self.data)
        #data_frame = pd.Panel(self.data, items=column_values,
                                   #major_axis=self.element, minor_axis=headers).to_frame()
        #data_frame.columns.names = column_names
        #data_frame.index.names = ['ElementID', 'Item']
        self.data_frame = data_frame

    def __eq__(self, table):  # pragma: no cover
        assert self.is_sort1 == table.is_sort1
        self._eq_header(table)
        if not np.array_equal(self.data, table.data):
            msg = 'table_name=%r class_name=%s\n' % (self.table_name, self.__class__.__name__)
            msg += '%s\n' % str(self.code_information())
            i = 0
            for itime in range(self.ntimes):
                for ie, eid in enumerate(self.element):
                    t1 = self.data[itime, ie, :]
                    t2 = table.data[itime, ie, :]
                    (mx1, my1, mxy1, bmx1, bmy1, bmxy1, tx1, ty1) = t1
                    (mx2, my2, mxy2, bmx2, bmy2, bmxy2, tx2, ty2) = t2

                    if not allclose(t1, t2):
                    #if not np.array_equal(t1.real, t2.real):
                        msg += ('%-8s (%s, %s, %s, %s, %s, %s, %s, %s)\n'
                                '%-8s (%s, %s, %s, %s, %s, %s, %s, %s)\n' % (
                                    eid,
                                    #mx1.real, my1.real, mxy1.real, bmx1.real, bmy1.real,
                                    #bmxy1.real, tx1.real, ty1.real,
                                    mx1, my1, mxy1, bmx1, bmy1, bmxy1, tx1, ty1,
                                    '',
                                    mx2, my2, mxy2, bmx2, bmy2, bmxy2, tx2, ty2,
                                    #mx2.real, my2.real, mxy2.real, bmx2.real, bmy2.real,
                                    #bmxy2.real, tx2.real, ty2.real,
                        ))
                        i += 1
                        if i > 10:
                            print(msg)
                            raise ValueError(msg)
                #print(msg)
                if i > 0:
                    raise ValueError(msg)
        return True

    def add_sort1(self, dt, eid, mx, my, mxy, bmx, bmy, bmxy, tx, ty):
        """unvectorized method for adding SORT1 transient data"""
        assert isinstance(eid, integer_types) and eid > 0, 'dt=%s eid=%s' % (dt, eid)
        self._times[self.itime] = dt
        self.element[self.ielement] = eid
        self.data[self.itime, self.ielement, :] = [mx, my, mxy, bmx, bmy, bmxy, tx, ty]
        self.ielement += 1
        self.itotal += 1

    #@property
    #def nnodes_per_element(self):
        #return 1

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

        msg = []
        if self.nonlinear_factor not in (None, np.nan):  # transient
            msg.append('  type=%s ntimes=%i nelements=%i; table_name=%r\n'
                       % (self.__class__.__name__, ntimes, nelements, self.table_name))
            ntimes_word = 'ntimes'
        else:
            msg.append('  type=%s nelements=%i; table_name=%r\n'
                       % (self.__class__.__name__, nelements, self.table_name))
            ntimes_word = '1'
        msg.append('  eType\n')
        headers = self.get_headers()
        n = len(headers)
        msg.append('  data: [%s, nnodes, %i] where %i=[%s]\n' % (ntimes_word, n, n, str(', '.join(headers))))
        msg.append('  data.shape = %s\n' % str(self.data.shape).replace('L', ''))
        #msg.append('  element type: %s\n' % self.element_type)
        msg.append('  element name: %s\n' % self.element_name)
        msg += self.get_data_code()
        return msg

    def get_f06_header(self, is_mag_phase=True, is_sort1=True):
        loads = ['    ELEMENT                - MEMBRANE  FORCES -                        - BENDING MOMENTS -               - TRANSVERSE SHEAR FORCES -\n'
                 '      ID              FX            FY            FXY             MX            MY            MXY             QX            QY\n',]
        if is_mag_phase:
            mag_real = ['                                                         (MAGNITUDE/PHASE)\n \n']
        else:
            mag_real = ['                                                          (REAL/IMAGINARY)\n \n']

        cquad4_bilinear = ['                  C O M P L E X   F O R C E S   I N   Q U A D R I L A T E R A L   E L E M E N T S   ( Q U A D 4 )\n']  # good
        cquad4_linear = ['                  C O M P L E X   F O R C E S   I N   Q U A D R I L A T E R A L   E L E M E N T S   ( Q U A D 4 )\n']  # good
        ctria3 = ['                     C O M P L E X   F O R C E S   I N   T R I A N G U L A R   E L E M E N T S   ( T R I A 3 )\n']  # good
        cquad8 = ['                  C O M P L E X   F O R C E S   I N   Q U A D R I L A T E R A L   E L E M E N T S   ( Q U A D 8 )\n']
        cquadr = ['                  C O M P L E X   F O R C E S   I N   Q U A D R I L A T E R A L   E L E M E N T S   ( Q U A D R )\n']
        ctria6 = ['                     C O M P L E X   F O R C E S   I N   T R I A N G U L A R   E L E M E N T S   ( T R I A 6 )\n']
        ctriar = ['                     C O M P L E X   F O R C E S   I N   T R I A N G U L A R   E L E M E N T S   ( T R I A R )\n']

        #is_bilinear = False
        if self.element_type == 144: # CQUAD4
            msg = cquad4_linear + mag_real + loads
        elif self.element_type == 33: # CQUAD4
            msg = cquad4_bilinear + mag_real + loads
        elif self.element_type == 64:  #CQUAD8
            msg = cquad8 + mag_real + loads
        elif self.element_type in [82, 228]:  # CQUADR, CQUADR linear
            msg = cquadr + mag_real + loads
        elif self.element_type == 74: # CTRIA3
            msg = ctria3 + mag_real + loads
        elif self.element_type == 75:  # CTRIA6
            msg = ctria6 + mag_real + loads
        elif self.element_type in [70, 227]:  # CTRIAR, CTRIAR linear
            msg = ctriar + mag_real + loads
        else:
            raise NotImplementedError('name=%r type=%s' % (self.element_name, self.element_type))
        return msg

    def get_element_index(self, eids):
        # elements are always sorted; nodes are not
        itot = searchsorted(eids, self.element)  #[0]
        return itot

    def eid_to_element_node_index(self, eids):
        #ind = ravel([searchsorted(self.element == eid) for eid in eids])
        ind = searchsorted(eids, self.element)
        #ind = ind.reshape(ind.size)
        #ind.sort()
        return ind

    def write_f06(self, f06_file, header=None, page_stamp='PAGE %s', page_num=1, is_mag_phase=False, is_sort1=True):
        if header is None:
            header = []
        msg_temp = self.get_f06_header(is_mag_phase=is_mag_phase, is_sort1=is_sort1)

        # write the f06
        #(ntimes, ntotal, two) = self.data.shape
        ntimes = self.data.shape[0]

        eids = self.element
        for itime in range(ntimes):
            dt = self._times[itime]  # TODO: rename this...
            header = _eigenvalue_header(self, header, itime, ntimes, dt)
            f06_file.write(''.join(header + msg_temp))

            #print("self.data.shape=%s itime=%s ieids=%s" % (str(self.data.shape), itime, str(ieids)))
            mx = self.data[itime, :, 0]
            my = self.data[itime, :, 1]
            mxy = self.data[itime, :, 2]
            bmx = self.data[itime, :, 3]
            bmy = self.data[itime, :, 4]
            bmxy = self.data[itime, :, 5]
            tx = self.data[itime, :, 6]
            ty = self.data[itime, :, 7]

            for eid, mxi, myi, mxyi, bmxi, bmyi, bmxyi, txi, tyi in zip(eids, mx, my, mxy, bmx, bmy, bmxy, tx, ty):
                out = write_imag_floats_13e([mxi, myi, mxyi, bmxi, bmyi, bmxyi, txi, tyi], is_mag_phase)
                [smxr, smyr, smxyr, sbmxr, sbmyr, sbmxyr, stxr, styr,
                 smxi, smyi, smxyi, sbmxi, sbmyi, sbmxyi, stxi, styi] = out
                #"""
                    #ELEMENT                - MEMBRANE  FORCES -                        - BENDING MOMENTS -               - TRANSVERSE SHEAR FORCES -
                      #ID              FX            FY            FXY             MX            MY            MXY             QX            QY
                #0       564       1.543439E+03  7.311177E+02  1.322702E+02    1.080178E+00  1.699104E+00  2.618547E-01    3.877034E+01  4.518554E+00
                                  #358.3129      358.0245      177.5593        177.5292      178.2112        0.0907        358.1465      179.4567
                #"""
                #                fx     fy     fxy     mx     my     mxy    qx      qy
                f06_file.write(
                    '0  %8i   %-13s  %-13s  %-13s   %-13s  %-13s  %-13s   %-13s  %s\n'
                    '   %8s   %-13s  %-13s  %-13s   %-13s  %-13s  %-13s   %-13s  %s\n' % (
                        eid, smxr, smyr, smxyr, sbmxr, sbmyr, sbmxyr, stxr, styr,
                        '', smxi, smyi, smxyi, sbmxi, sbmyi, sbmxyi, stxi, styi))
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

        #if isinstance(self.nonlinear_factor, float):
            #op2_format = '%sif' % (7 * self.ntimes)
            #raise NotImplementedError()
        #else:
            #op2_format = 'i21f'
        #s = Struct(op2_format)

        eids = self.element
        eids_device = eids * 10 + self.device_code

        # table 4 info
        #ntimes = self.data.shape[0]
        #nnodes = self.data.shape[1]
        nelements = self.data.shape[1]

        # 21 = 1 node, 3 principal, 6 components, 9 vectors, 2 p/ovm
        #ntotal = ((nnodes * 21) + 1) + (nelements * 4)

        ntotali = self.num_wide
        ntotal = ntotali * nelements

        #print('shape = %s' % str(self.data.shape))
        #assert self.ntimes == 1, self.ntimes

        op2_ascii.write('  ntimes = %s\n' % self.ntimes)

        #fmt = '%2i %6f'
        #print('ntotal=%s' % (ntotal))
        #assert ntotal == 193, ntotal

        if self.is_sort1:
            struct1 = Struct(endian + b'i 16f')
        else:
            raise NotImplementedError('SORT2')

        op2_ascii.write('%s-nelements=%i\n' % (self.element_name, nelements))
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

            mx = self.data[itime, :, 0]
            my = self.data[itime, :, 1]
            mxy = self.data[itime, :, 2]
            bmx = self.data[itime, :, 3]
            bmy = self.data[itime, :, 4]
            bmxy = self.data[itime, :, 5]
            tx = self.data[itime, :, 6]
            ty = self.data[itime, :, 7]

            for eid_device, mxi, myi, mxyi, bmxi, bmyi, bmxyi, txi, tyi in zip(eids_device, mx, my, mxy, bmx, bmy, bmxy, tx, ty):
                data = [eid_device,
                        mxi.real, myi.real, mxyi.real, bmxi.real, bmyi.real, bmxyi.real, txi.real, tyi.real,
                        mxi.imag, myi.imag, mxyi.imag, bmxi.imag, bmyi.imag, bmxyi.imag, txi.imag, tyi.imag]
                op2_ascii.write('  eid_device=%s data=%s\n' % (eid_device, str(data)))
                op2.write(struct1.pack(*data))

            itable -= 1
            header = [4 * ntotal,]
            op2.write(pack('i', *header))
            op2_ascii.write('footer = %s\n' % header)
            new_result = False
        return itable


class ComplexPlate2ForceArray(ComplexForceObject):
    def __init__(self, data_code, is_sort1, isubcase, dt):
        self.element_type = None
        self.element_name = None
        ComplexForceObject.__init__(self, data_code, isubcase)
        #self.code = [self.format_code, self.sort_code, self.s_code]

        #self.ntimes = 0  # or frequency/mode
        #self.ntotal = 0
        self.nelements = 0  # result specific

        #if is_sort1:
            #pass
        #else:
            #raise NotImplementedError('SORT2')

    def get_headers(self) -> List[str]:
        headers = ['mx', 'my', 'mxy', 'bmx', 'bmy', 'bmxy', 'tx', 'ty']
        return headers

    def build(self):
        """sizes the vectorized attributes of the ComplexPlate2ForceArray"""
        if self.is_built:
            return

        assert self.ntimes > 0, 'ntimes=%s' % self.ntimes
        assert self.nelements > 0, 'nelements=%s' % self.nelements
        assert self.ntotal > 0, 'ntotal=%s' % self.ntotal
        #self.names = []
        self.nelements //= self.ntimes
        #print('ntimes=%s nelements=%s ntotal=%s' % (self.ntimes, self.nelements, self.ntotal))
        self.itime = 0
        self.ielement = 0
        self.itotal = 0
        #self.ntimes = 0
        #self.nelements = 0
        self.is_built = True

        #print("ntimes=%s nelements=%s ntotal=%s" % (self.ntimes, self.nelements, self.ntotal))
        dtype, idtype, cfdtype = get_complex_times_dtype(self.nonlinear_factor, self.size)
        self._times = zeros(self.ntimes, dtype=dtype)

        self.element = zeros(self.nelements, dtype=idtype)
        self.element_node = zeros((self.ntotal, 2), dtype=idtype)

        #[mx, my, mxy, bmx, bmy, bmxy, tx, ty]
        self.data = zeros((self.ntimes, self.ntotal, 8), dtype=cfdtype)

    def build_dataframe(self):
        """creates a pandas dataframe"""
        import pandas as pd
        headers = self.get_headers()
        assert 0 not in self.element
        #print(self.element_node)
        element_node = [self.element_node[:, 0], self.element_node[:, 1]]
        assert 0 not in self.element_node[:, 0]
        if self.nonlinear_factor not in (None, np.nan):
            # Freq                  0.00001  10.00000 20.00000 30.00000                  40.00000 50.00000 60.00000
            # ElementID NodeID Item
            # 6         0      mx         0j       0j       0j       0j    (-705.7376-1327.1312j)       0j       0j
            #                  my         0j       0j       0j       0j      (7404.8853+13924.8j)       0j       0j
            #                  mxy        0j       0j       0j       0j  (-101.319756-190.53061j)       0j       0j
            #                  bmx        0j       0j       0j       0j    (3.0701134+5.7733126j)       0j       0j
            #                  bmy        0j       0j       0j       0j     (98.75731+185.71196j)       0j       0j
            #                  bmxy       0j       0j       0j       0j   (0.25202343+0.4739271j)       0j       0j
            #                  tx         0j       0j       0j       0j    (14.426779+27.129389j)       0j       0j
            #                  ty         0j       0j       0j       0j     (-199.6823-375.5002j)       0j       0j
            #           4      mx         0j       0j       0j       0j    (-2934.639-5518.5537j)       0j       0j
            #                  my         0j       0j       0j       0j    (7516.2485+14134.217j)       0j       0j
            #                  mxy        0j       0j       0j       0j  (-101.319756-190.53061j)       0j       0j
            #                  bmx        0j       0j       0j       0j    (-19.69526-37.036705j)       0j       0j
            #                  bmy        0j       0j       0j       0j     (100.64615+189.2639j)       0j       0j
            #                  bmxy       0j       0j       0j       0j   (0.25202343+0.4739271j)       0j       0j
            #                  tx         0j       0j       0j       0j    (14.426779+27.129389j)       0j       0j
            #                  ty         0j       0j       0j       0j     (-199.6823-375.5002j)       0j       0j
            column_names, column_values = self._build_dataframe_transient_header()
            data_frame = self._build_pandas_transient_element_node(
                column_values, column_names,
                headers, self.element_node, self.data)

            #data_frame = pd.Panel(self.data, items=column_values,
                                  #major_axis=element_node, minor_axis=headers).to_frame()
            #data_frame.columns.names = column_names
            #data_frame.index.names = ['ElementID', 'NodeID', 'Item']
        else:
            data_frame = pd.Panel(self.data,
                                  major_axis=element_node, minor_axis=headers).to_frame()
            data_frame.columns.names = ['Static']
            data_frame.index.names = ['ElementID', 'NodeID', 'Item']
        self.data_frame = data_frame

    def __eq__(self, table):  # pragma: no cover
        self._eq_header(table)
        assert self.is_sort1 == table.is_sort1
        if not np.array_equal(self.data, table.data):
            msg = 'table_name=%r class_name=%s\n' % (self.table_name, self.__class__.__name__)
            msg += '%s\n' % str(self.code_information())
            i = 0
            for itime in range(self.ntimes):
                for ie, e in enumerate(self.element_node):
                    (eid, nid) = e
                    t1 = self.data[itime, ie, :]
                    t2 = table.data[itime, ie, :]
                    (mx1, my1, mxy1, bmx1, bmy1, bmxy1, tx1, ty1) = t1
                    (mx2, my2, mxy2, bmx2, bmy2, bmxy2, tx2, ty2) = t2

                    if not allclose(t1, t2):
                        base1 = '(%s, %s)   ' % (eid, nid)
                        base2 = ' ' * len(base1)
                        msg += (
                            '%s (%s, %s, %s, %s, %s, %s, %s, %s)\n'
                            '%s(%s, %s, %s, %s, %s, %s, %s, %s)\n' % (
                                base1,
                                mx1, my1, mxy1, bmx1, bmy1, bmxy1, tx1, ty1,
                                base2,
                                mx2, my2, mxy2, bmx2, bmy2, bmxy2, tx2, ty2))
                        i += 1
                        if i > 10:
                            print(msg)
                            raise ValueError(msg)
                if i > 0:
                    raise ValueError(msg)
        return True

    def add_new_element_sort1(self, dt, eid, term, nid, mx, my, mxy, bmx, bmy, bmxy, tx, ty):
        self._times[self.itime] = dt
        self.element[self.ielement] = eid
        self.element_node[self.itotal, :] = [eid, nid]
        self.data[self.itime, self.itotal, :] = [mx, my, mxy, bmx, bmy, bmxy, tx, ty]
        self.itotal += 1
        self.ielement += 1

    def add_sort1(self, dt, eid, nid, mx, my, mxy, bmx, bmy, bmxy, tx, ty):
        """unvectorized method for adding SORT1 transient data"""
        assert isinstance(eid, integer_types) and eid > 0, 'dt=%s eid=%s' % (dt, eid)
        self._times[self.itime] = dt
        #assert self.element[self.ielement - 1] == eid, eid
        self.element_node[self.itotal, :] = [eid, nid]
        self.data[self.itime, self.itotal, :] = [mx, my, mxy, bmx, bmy, bmxy, tx, ty]
        self.itotal += 1

    @property
    def nnodes_per_element(self):
        if self.element_type == 144:  # CQUAD4
            nnodes_element = 5
        elif self.element_type == 64:  # CQUAD8
            nnodes_element = 5
        elif self.element_type == 82:  # CQUADR
            nnodes_element = 5
        elif self.element_type == 75:  # CTRIA6
            nnodes_element = 4
        elif self.element_type == 70:  # CTRIAR
            nnodes_element = 4
        else:  # pragma: no cover
            raise NotImplementedError('element_type=%s element_name=%s' % (self.element_type, self.element_name))
        return nnodes_element

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

        msg = []
        if self.nonlinear_factor not in (None, np.nan):  # transient
            msg.append('  type=%s ntimes=%i nelements=%i; table_name=%r\n'
                       % (self.__class__.__name__, ntimes, nelements, self.table_name))
            ntimes_word = 'ntimes'
        else:
            msg.append('  type=%s nelements=%i; table_name=%r\n'
                       % (self.__class__.__name__, nelements, self.table_name))
            ntimes_word = '1'
        msg.append('  eType\n')
        headers = self.get_headers()
        n = len(headers)
        msg.append('  data: [%s, nnodes, %i] where %i=[%s]\n' % (
            ntimes_word, n, n, str(', '.join(headers))))
        msg.append('  data.shape = %s\n' % str(self.data.shape).replace('L', ''))
        #msg.append('  element type: %s\n' % self.element_type)
        msg.append('  element name: %s\n' % self.element_name)
        msg += self.get_data_code()
        return msg

    def get_f06_header(self, is_mag_phase=True, is_sort1=True):
        loads = [
            '    ELEMENT                - MEMBRANE  FORCES -                        - BENDING MOMENTS -               - TRANSVERSE SHEAR FORCES -\n'
            '      ID              FX            FY            FXY             MX            MY            MXY             QX            QY\n',]
        if is_mag_phase:
            mag_real = ['                                                         (MAGNITUDE/PHASE)\n \n']
        else:
            mag_real = ['                                                          (REAL/IMAGINARY)\n \n']

        cquad4_bilinear = ['                  C O M P L E X   F O R C E S   I N   Q U A D R I L A T E R A L   E L E M E N T S   ( Q U A D 4 )\n']  # good
        cquad4_linear = ['                  C O M P L E X   F O R C E S   I N   Q U A D R I L A T E R A L   E L E M E N T S   ( Q U A D 4 )\n']  # good
        ctria3 = ['                     C O M P L E X   F O R C E S   I N   T R I A N G U L A R   E L E M E N T S   ( T R I A 3 )\n']  # good
        cquad8 = ['                  C O M P L E X   F O R C E S   I N   Q U A D R I L A T E R A L   E L E M E N T S   ( Q U A D 8 )\n']
        cquadr = ['                  C O M P L E X   F O R C E S   I N   Q U A D R I L A T E R A L   E L E M E N T S   ( Q U A D R )\n']
        ctria6 = ['                     C O M P L E X   F O R C E S   I N   T R I A N G U L A R   E L E M E N T S   ( T R I A 6 )\n']
        ctriar = ['                     C O M P L E X   F O R C E S   I N   T R I A N G U L A R   E L E M E N T S   ( T R I A R )\n']

        #is_bilinear = False
        if self.element_type == 144: # CQUAD4
            msg = cquad4_linear + mag_real + loads
        elif self.element_type == 33: # CQUAD4
            msg = cquad4_bilinear + mag_real + loads
        elif self.element_type == 64:  #CQUAD8
            msg = cquad8 + mag_real + loads
        elif self.element_type == 82:  # CQUADR
            msg = cquadr + mag_real + loads
        elif self.element_type == 74: # CTRIA3
            msg = ctria3 + mag_real + loads
        elif self.element_type == 75:  # CTRIA6
            msg = ctria6 + mag_real + loads
        elif self.element_type == 70:  # CTRIAR
            msg = ctriar + mag_real + loads
        else:
            raise NotImplementedError('name=%r type=%s' % (self.element_name, self.element_type))
        return msg

    def get_element_index(self, eids):
        # elements are always sorted; nodes are not
        itot = searchsorted(eids, self.element)  #[0]
        return itot

    def eid_to_element_node_index(self, eids):
        #ind = ravel([searchsorted(self.element == eid) for eid in eids])
        ind = searchsorted(eids, self.element)
        #ind = ind.reshape(ind.size)
        #ind.sort()
        return ind

    def write_f06(self, f06_file, header=None, page_stamp='PAGE %s', page_num=1, is_mag_phase=False, is_sort1=True):
        if header is None:
            header = []
        msg_temp = self.get_f06_header(is_mag_phase=is_mag_phase, is_sort1=is_sort1)

        # write the f06
        #(ntimes, ntotal, two) = self.data.shape
        ntimes = self.data.shape[0]

        eids = self.element
        for itime in range(ntimes):
            dt = self._times[itime]  # TODO: rename this...
            header = _eigenvalue_header(self, header, itime, ntimes, dt)
            f06_file.write(''.join(header + msg_temp))

            #print("self.data.shape=%s itime=%s ieids=%s" % (str(self.data.shape), itime, str(ieids)))
            mx = self.data[itime, :, 0]
            my = self.data[itime, :, 1]
            mxy = self.data[itime, :, 2]
            bmx = self.data[itime, :, 3]
            bmy = self.data[itime, :, 4]
            bmxy = self.data[itime, :, 5]
            tx = self.data[itime, :, 6]
            ty = self.data[itime, :, 7]

            for eid, mxi, myi, mxyi, bmxi, bmyi, bmxyi, txi, tyi in zip(eids, mx, my, mxy, bmx, bmy, bmxy, tx, ty):
                out = write_imag_floats_13e([mxi, myi, mxyi, bmxi, bmyi, bmxyi, txi, tyi], is_mag_phase)
                [smxr, smyr, smxyr, sbmxr, sbmyr, sbmxyr, stxr, styr,
                 smxi, smyi, smxyi, sbmxi, sbmyi, sbmxyi, stxi, styi] = out
                #"""
                    #ELEMENT                - MEMBRANE  FORCES -                        - BENDING MOMENTS -               - TRANSVERSE SHEAR FORCES -
                      #ID              FX            FY            FXY             MX            MY            MXY             QX            QY
                #0       564       1.543439E+03  7.311177E+02  1.322702E+02    1.080178E+00  1.699104E+00  2.618547E-01    3.877034E+01  4.518554E+00
                                  #358.3129      358.0245      177.5593        177.5292      178.2112        0.0907        358.1465      179.4567
                #"""
                #                fx     fy     fxy     mx     my     mxy    qx      qy
                f06_file.write(
                    '0  %8i   %-13s  %-13s  %-13s   %-13s  %-13s  %-13s   %-13s  %s\n'
                    '   %8s   %-13s  %-13s  %-13s   %-13s  %-13s  %-13s   %-13s  %s\n' % (
                        eid, smxr, smyr, smxyr, sbmxr, sbmyr, sbmxyr, stxr, styr,
                        '', smxi, smyi, smxyi, sbmxi, sbmyi, sbmxyi, stxi, styi))
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

        #if isinstance(self.nonlinear_factor, float):
            #op2_format = '%sif' % (7 * self.ntimes)
            #raise NotImplementedError()
        #else:
            #op2_format = 'i21f'
        #s = Struct(op2_format)

        #eids = self.element
        eids_device = self.element * 10 + self.device_code
        eids = self.element_node[:, 0]
        nids = self.element_node[:, 1]

        # table 4 info
        #ntimes = self.data.shape[0]
        #nnodes = self.data.shape[1]
        nelements = len(self.element)

        # 21 = 1 node, 3 principal, 6 components, 9 vectors, 2 p/ovm
        #ntotal = ((nnodes * 21) + 1) + (nelements * 4)

        ntotali = self.num_wide

        nnodes_all = 5
        numwide_imag = 2 + nnodes_all * 17
        assert ntotali == numwide_imag

        ntotal = ntotali * nelements

        #print('shape = %s' % str(self.data.shape))
        #assert self.ntimes == 1, self.ntimes

        op2_ascii.write('  ntimes = %s\n' % self.ntimes)

        #fmt = '%2i %6f'
        #print('ntotal=%s' % (ntotal))
        #assert ntotal == 193, ntotal

        if self.is_sort1:
            struct1 = Struct(endian + b'i 4s i 16f')
            struct2 = Struct(endian + b'i 16f')
        else:
            raise NotImplementedError('SORT2')

        op2_ascii.write('%s-nelements=%i\n' % (self.element_name, nelements))
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

            mx = self.data[itime, :, 0]
            my = self.data[itime, :, 1]
            mxy = self.data[itime, :, 2]
            bmx = self.data[itime, :, 3]
            bmy = self.data[itime, :, 4]
            bmxy = self.data[itime, :, 5]
            tx = self.data[itime, :, 6]
            ty = self.data[itime, :, 7]

            nwide = 0
            ielement = -1
            for eid, nid, mxi, myi, mxyi, bmxi, bmyi, bmxyi, txi, tyi in zip(eids, nids, mx, my, mxy, bmx, bmy, bmxy, tx, ty):
                if nid == 0:
                    ielement += 1
                    eid_device = eids_device[ielement]
                    data = [eid_device, b'CEN/', nid,
                            mxi.real, myi.real, mxyi.real, bmxi.real, bmyi.real, bmxyi.real, txi.real, tyi.real,
                            mxi.imag, myi.imag, mxyi.imag, bmxi.imag, bmyi.imag, bmxyi.imag, txi.imag, tyi.imag]
                    op2_ascii.write('  eid_device=%s data=%s\n' % (eid_device, str(data)))
                    op2.write(struct1.pack(*data))
                else:
                    data = [nid,
                            mxi.real, myi.real, mxyi.real, bmxi.real, bmyi.real, bmxyi.real, txi.real, tyi.real,
                            mxi.imag, myi.imag, mxyi.imag, bmxi.imag, bmyi.imag, bmxyi.imag, txi.imag, tyi.imag]
                    op2_ascii.write('    data=%s\n' % (str(data)))
                    op2.write(struct2.pack(*data))
                nwide += len(data)
            assert nwide == ntotal, 'nwide=%s ntotal=%s' % (nwide, ntotal)
            itable -= 1
            header = [4 * ntotal,]
            op2.write(pack('i', *header))
            op2_ascii.write('footer = %s\n' % header)
            new_result = False
        return itable



class ComplexCBarWeldForceArray(ComplexForceObject):
    def __init__(self, data_code, is_sort1, isubcase, dt):
        self.element_type = None
        self.element_name = None
        ComplexForceObject.__init__(self, data_code, isubcase)

        self.result_flag = 0
        #self.code = [self.format_code, self.sort_code, self.s_code]

        #self.ntimes = 0  # or frequency/mode
        #self.ntotal = 0
        self.itime = 0
        self.nelements = 0  # result specific
        #self.element_type = 'CBAR'
        #self.cid = {}  # gridGauss

        #if is_sort1:
            ##sort1
            #pass
        #else:
            #raise NotImplementedError('SORT2')

    def get_headers(self) -> List[str]:
        headers = ['bending_moment_1a', 'bending_moment_2a',
                   'bending_moment_1b', 'bending_moment_2b',
                   'shear1', 'shear2', 'axial', 'torque', ]
        return headers

    def build(self):
        """sizes the vectorized attributes of the ComplexCBarForceArray"""
        #print('ntimes=%s nelements=%s ntotal=%s subtitle=%s' % (
            #self.ntimes, self.nelements, self.ntotal, self.subtitle))
        if self.is_built:
            return
        nnodes = 1

        #self.names = []
        #self.nelements //= nnodes
        self.nelements //= self.ntimes
        #self.ntotal //= self.ntimes
        self.itime = 0
        self.ielement = 0
        self.itotal = 0
        self.is_built = True
        #print('ntotal=%s ntimes=%s nelements=%s' % (self.ntotal, self.ntimes, self.nelements))

        #print("ntimes=%s nelements=%s ntotal=%s" % (self.ntimes, self.nelements, self.ntotal))
        dtype, idtype, cfdtype = get_complex_times_dtype(self.nonlinear_factor, self.size)
        self._times = zeros(self.ntimes, dtype=dtype)
        self.element = zeros(self.ntotal, dtype=idtype)

        # the number is messed up because of the offset for the element's properties

        #if not self.nelements * nnodes == self.ntotal:
            #msg = 'ntimes=%s nelements=%s nnodes=%s ne*nn=%s ntotal=%s' % (
                #self.ntimes, self.nelements, nnodes, self.nelements * nnodes, self.ntotal)
            #raise RuntimeError(msg)
        #[bm1a, bm2a, bm1b, bm2b, ts1, ts2, af, trq]
        self.data = zeros((self.ntimes, self.ntotal, 8), dtype=cfdtype)


    def build_dataframe(self):
        """creates a pandas dataframe"""
        headers = self.get_headers()
        column_names, column_values = self._build_dataframe_transient_header()
        data_frame = self._build_pandas_transient_elements(column_values, column_names,
                                                           headers, self.element, self.data)
        #self.data_frame = pd.Panel(self.data, items=column_values,
                                   #major_axis=self.element, minor_axis=headers).to_frame()
        #self.data_frame.columns.names = column_names
        #self.data_frame.index.names = ['ElementID', 'Item']
        self.data_frame = data_frame

    def __eq__(self, table):  # pragma: no cover
        assert self.is_sort1 == table.is_sort1
        self._eq_header(table)
        if not np.array_equal(self.data, table.data):
            msg = 'table_name=%r class_name=%s\n' % (self.table_name, self.__class__.__name__)
            msg += '%s\n' % str(self.code_information())
            ntimes = self.data.shape[0]

            i = 0
            if self.is_sort1:
                for itime in range(ntimes):
                    for ieid, eid in enumerate(self.element):
                        t1 = self.data[itime, ieid, :]
                        t2 = table.data[itime, ieid, :]
                        (s1a1, s2a1, s3a1, s4a1, axial1, s2a1, s2b1, s2c1, s2d1) = t1
                        (s1a2, s2a2, s3a2, s4a2, axial2, s2a2, s2b2, s2c2, s2d2) = t2
                        #d = t1 - t2
                        if not allclose([s1a1.real, s2a1.real, s3a1.real, s4a1.real, axial1.real, s2a1.real, s2b1.real, s2c1.real, s2d1.real],
                                        [s1a2.real, s2a2.real, s3a2.real, s4a2.real, axial2.real, s2a2.real, s2b2.real, s2c2.real, s2d2.real], atol=0.0001):
                        #if not np.array_equal(t1, t2):
                            msg += '%-4s  (%s, %s, %s, %s, %s, %s, %s, %s, %s)\n      (%s, %s, %s, %s, %s, %s, %s, %s, %s)\n' % (
                                eid,
                                s1a1.real, s2a1.real, s3a1.real, s4a1.real, axial1.real, s2a1.real, s2b1.real, s2c1.real, s2d1.real,
                                s1a2.real, s2a2.real, s3a2.real, s4a2.real, axial2.real, s2a2.real, s2b2.real, s2c2.real, s2d2.real,
                                )
                            i += 1
                        if i > 10:
                            print(msg)
                            raise ValueError(msg)
            else:
                raise NotImplementedError(self.is_sort2)
            if i > 0:
                print(msg)
                raise ValueError(msg)
        return True

    def add_sort1(self, dt, eid, bm1a, bm2a, bm1b, bm2b, ts1, ts2, af, trq):
        """unvectorized method for adding SORT1 transient data"""
        assert isinstance(eid, integer_types) and eid > 0, 'dt=%s eid=%s' % (dt, eid)
        self._times[self.itime] = dt
        self.data[self.itime, self.itotal, :] = [bm1a, bm2a, bm1b, bm2b, ts1, ts2, af, trq]
        self.element[self.itotal] = eid
        self.itotal += 1

    def get_stats(self, short=False) -> List[str]:
        if not self.is_built:
            return [
                '<%s>\n' % self.__class__.__name__,
                '  ntimes: %i\n' % self.ntimes,
                '  ntotal: %i\n' % self.ntotal,
            ]

        nelements = self.nelements
        ntimes = self.ntimes
        msg = []

        if self.nonlinear_factor not in (None, np.nan):  # transient
            msg.append('  type=%s ntimes=%i nelements=%i; table_name=%r\n'
                       % (self.__class__.__name__, ntimes, nelements, self.table_name))
        else:
            msg.append('  type=%s nelements=%i; table_name=%r\n' % (
                self.__class__.__name__, nelements, self.table_name))
        msg.append('  eType, cid\n')
        msg.append('  data: [ntimes, nelements, 8] where 8=[%s]\n' % str(', '.join(self.get_headers())))
        msg.append('  data.shape = %s\n' % str(self.data.shape).replace('L', ''))
        msg.append('  is_sort1=%s is_sort2=%s\n' % (self.is_sort1, self.is_sort2))
        msg.append(f'  {self.element_name}\n')
        msg += self.get_data_code()
        return msg

    def write_f06(self, f06_file, header=None, page_stamp='PAGE %s',
                  page_num=1, is_mag_phase=False, is_sort1=True):
        if header is None:
            header = []
        #msg_temp, nnodes = get_f06_header(self, is_mag_phase, is_sort1)


        #is_sort1 = False
        if is_mag_phase:
            mag_phase = '                                                          (MAGNITUDE/PHASE)\n \n'
        else:
            mag_phase = '                                                          (REAL/IMAGINARY)\n \n'


        name = self.data_code['name']
        if name == 'freq':
            name = 'FREQUENCY'
        #else: # mode
            #raise RuntimeError(name)

        if is_sort1:
            line1 = '0    ELEMENT         BEND-MOMENT-END-A            BEND-MOMENT-END-B                  SHEAR\n'
            line2 = '       ID.         PLANE 1       PLANE 2        PLANE 1       PLANE 2        PLANE 1       PLANE 2        FORCE          TORQUE\n'
        else:
            line1 = '                    BEND-MOMENT-END-A            BEND-MOMENT-END-B                  SHEAR\n'
            line2 = '   %16s       PLANE 1       PLANE 2        PLANE 1       PLANE 2        PLANE 1       PLANE 2        FORCE          TORQUE\n' % name

        # force
        msg_temp = header + [
            '                             C O M P L E X   F O R C E S   I N   B A R   E L E M E N T S   ( C B A R )\n',
            mag_phase,
            ' ',
            line1,
            line2,
        ]
        if self.is_sort1:
            assert self.is_sort1 == True, str(self)
            if is_sort1:
                page_num = self._write_sort1_as_sort1(f06_file, page_num, page_stamp, header, msg_temp, is_mag_phase)
            else:
                self._write_sort1_as_sort2(f06_file, page_num, page_stamp, header, msg_temp, is_mag_phase)
        else:
            assert self.is_sort1 == True, str(self)
        return page_num - 1

    def _write_sort1_as_sort1(self, f06_file, page_num, page_stamp, header, msg_temp, is_mag_phase):
        eids = self.element
        #times = self._times
        ntimes = self.data.shape[0]
        for itime in range(ntimes):
            dt = self._times[itime]
            dt_line = ' %14s = %12.5E\n' % (self.data_code['name'], dt)
            header[1] = dt_line
            msg = header + msg_temp
            f06_file.write(''.join(msg))

            #bm1a, bm2a, bm1b, bm2b, ts1, ts2, af, trq
            assert self.is_sort1 == True, str(self)
            bm1a = self.data[itime, :, 0]
            bm2a = self.data[itime, :, 1]
            bm1b = self.data[itime, :, 2]
            bm2b = self.data[itime, :, 3]
            ts1 = self.data[itime, :, 4]
            ts2 = self.data[itime, :, 5]
            af = self.data[itime, :, 6]
            trq = self.data[itime, :, 7]

            for eid, bm1ai, bm2ai, bm1bi, bm2bi, ts1i, ts2i, afi, trqi in zip(eids, bm1a, bm2a, bm1b, bm2b, ts1, ts2, af, trq):
                vals = (bm1ai, bm2ai, bm1bi, bm2bi, ts1i, ts2i, afi, trqi)
                vals2 = write_imag_floats_13e(vals, is_mag_phase)
                (bm1air, bm2air, bm1bir, bm2bir, ts1ir, ts2ir, afir, trqir,
                 bm1aii, bm2aii, bm1bii, bm2bii, ts1ii, ts2ii, afii, trqii) = vals2

                f06_file.write('0%16i   %-13s  %-13s  %-13s  %-13s  %-13s  %-13s  %-13s  %s\n'
                               ' %14s   %-13s  %-13s  %-13s  %-13s  %-13s  %-13s  %-13s  %s\n' % (
                                   eid, bm1air, bm2air, bm1bir, bm2bir, ts1ir, ts2ir, afir, trqir,
                                   '', bm1aii, bm2aii, bm1bii, bm2bii, ts1ii, ts2ii, afii, trqii))
            f06_file.write(page_stamp % page_num)
            page_num += 1
        return page_num

    def _write_sort1_as_sort2(self, f06_file, page_num, page_stamp, header, msg_temp, is_mag_phase):
        eids = self.element
        times = self._times
        #ntimes = self.data.shape[0]
        for ieid, eid in enumerate(eids):
            eid_line = ' ELEMENT-ID = %s' % (eid)
            header[1] = eid_line
            msg = header + msg_temp
            f06_file.write(''.join(msg))

            #bm1a, bm2a, bm1b, bm2b, ts1, ts2, af, trq
            bm1a = self.data[:, ieid, 0]
            bm2a = self.data[:, ieid, 1]
            bm1b = self.data[:, ieid, 2]
            bm2b = self.data[:, ieid, 3]
            ts1 = self.data[:, ieid, 4]
            ts2 = self.data[:, ieid, 5]
            af = self.data[:, ieid, 6]
            trq = self.data[:, ieid, 7]

            for dt, bm1ai, bm2ai, bm1bi, bm2bi, ts1i, ts2i, afi, trqi in zip(times, bm1a, bm2a, bm1b, bm2b, ts1, ts2, af, trq):
                vals = (bm1ai, bm2ai, bm1bi, bm2bi, ts1i, ts2i, afi, trqi)
                vals2 = write_imag_floats_13e(vals, is_mag_phase)
                (bm1air, bm2air, bm1bir, bm2bir, ts1ir, ts2ir, afir, trqir,
                 bm1aii, bm2aii, bm1bii, bm2bii, ts1ii, ts2ii, afii, trqii) = vals2

                f06_file.write('0%16s  %-13s  %-13s  %-13s  %-13s  %-13s  %-13s  %-13s  %s\n'
                               ' %15s  %-13s  %-13s  %-13s  %-13s  %-13s  %-13s  %-13s  %s\n' % (
                                   write_float_12e(dt),
                                   bm1air, bm2air, bm1bir, bm2bir, ts1ir, ts2ir, afir, trqir,
                                   '', bm1aii, bm2aii, bm1bii, bm2bii, ts1ii, ts2ii, afii, trqii))
            f06_file.write(page_stamp % page_num)
            page_num += 1
        return page_num

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

        #if isinstance(self.nonlinear_factor, float):
            #op2_format = '%sif' % (7 * self.ntimes)
            #raise NotImplementedError()
        #else:
            #op2_format = 'i21f'
        #s = Struct(op2_format)

        eids = self.element
        eids_device = eids * 10 + self.device_code

        # table 4 info
        #ntimes = self.data.shape[0]
        #nnodes = self.data.shape[1]
        nelements = self.data.shape[1]

        # 21 = 1 node, 3 principal, 6 components, 9 vectors, 2 p/ovm
        #ntotal = ((nnodes * 21) + 1) + (nelements * 4)

        ntotali = self.num_wide
        ntotal = ntotali * nelements

        #print('shape = %s' % str(self.data.shape))
        #assert self.ntimes == 1, self.ntimes

        op2_ascii.write('  ntimes = %s\n' % self.ntimes)

        #fmt = '%2i %6f'
        #print('ntotal=%s' % (ntotal))
        #assert ntotal == 193, ntotal

        if self.is_sort1:
            struct1 = Struct(endian + b'i 16f')
        else:
            raise NotImplementedError('SORT2')

        op2_ascii.write('%s-nelements=%i\n' % (self.element_name, nelements))
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

            #bm1a, bm2a, bm1b, bm2b, ts1, ts2, af, trq
            bm1a = self.data[itime, :, 0]
            bm2a = self.data[itime, :, 1]
            bm1b = self.data[itime, :, 2]
            bm2b = self.data[itime, :, 3]
            ts1 = self.data[itime, :, 4]
            ts2 = self.data[itime, :, 5]
            af = self.data[itime, :, 6]
            trq = self.data[itime, :, 7]
            assert len(eids_device) == len(bm1a.real)

            for eid_device, bm1ai, bm2ai, bm1bi, bm2bi, ts1i, ts2i, afi, trqi in zip(
                eids_device, bm1a, bm2a, bm1b, bm2b, ts1, ts2, af, trq):

                data = [eid_device,
                        bm1ai.real, bm2ai.real, bm1bi.real, bm2bi.real, ts1i.real, ts2i.real, afi.real, trqi.real,
                        bm1ai.imag, bm2ai.imag, bm1bi.imag, bm2bi.imag, ts1i.imag, ts2i.imag, afi.imag, trqi.imag]
                op2_ascii.write('  eid_device=%s data=%s\n' % (eid_device, str(data)))
                op2.write(struct1.pack(*data))

            itable -= 1
            header = [4 * ntotal,]
            op2.write(pack('i', *header))
            op2_ascii.write('footer = %s\n' % header)
            new_result = False
        return itable

class ComplexCBarForceArray(ComplexCBarWeldForceArray):
    def __init__(self, data_code, is_sort1, isubcase, dt):
        ComplexCBarWeldForceArray.__init__(self, data_code, is_sort1, isubcase, dt)

class ComplexCWeldForceArray(ComplexCBarWeldForceArray):
    def __init__(self, data_code, is_sort1, isubcase, dt):
        ComplexCBarWeldForceArray.__init__(self, data_code, is_sort1, isubcase, dt)



class ComplexCBeamForceArray(ComplexForceObject):
    def __init__(self, data_code, is_sort1, isubcase, dt):
        self.element_type = None
        self.element_name = None
        ComplexForceObject.__init__(self, data_code, isubcase)

        self.result_flag = 0
        self.itime = 0
        self.nelements = 0  # result specific
        #self.element_type = 'CBEAM'

        #if is_sort1:
            ##sort1
            #pass
        #else:
            #raise NotImplementedError('SORT2')

    def get_headers(self) -> List[str]:
        headers = [
            'sd', 'bending_moment1', 'bending_moment2', 'shear1', 'shear2',
            'axial_force', 'total_torque', 'warping_torque', ]
        return headers

    def build(self):
        """sizes the vectorized attributes of the ComplexCBeamForceArray"""
        #print('ntimes=%s nelements=%s ntotal=%s subtitle=%s' % (
            #self.ntimes, self.nelements, self.ntotal, self.subtitle))
        nnodes = 11

        #self.names = []
        #self.nelements //= nnodes
        self.nelements //= self.ntimes
        #self.ntotal //= self.ntimes
        self.itime = 0
        self.ielement = 0
        self.itotal = 0
        self.is_built = True
        #print('ntotal=%s ntimes=%s nelements=%s' % (self.ntotal, self.ntimes, self.nelements))

        #print("ntimes=%s nelements=%s ntotal=%s" % (self.ntimes, self.nelements, self.ntotal))
        dtype, idtype, cfdtype = get_complex_times_dtype(self.nonlinear_factor, self.size)
        self._times = zeros(self.ntimes, dtype)
        self.element = zeros(self.ntotal, idtype)
        self.element_node = zeros((self.ntotal, 2), idtype)

        # the number is messed up because of the offset for the element's properties

        #if not self.nelements * nnodes == self.ntotal:
            #msg = 'ntimes=%s nelements=%s nnodes=%s ne*nn=%s ntotal=%s' % (
                #self.ntimes, self.nelements, nnodes, self.nelements * nnodes, self.ntotal)
            #raise RuntimeError(msg)
        #[sd, bm1, bm2, ts1, ts2, af, ttrq, wtrq]
        self.data = zeros((self.ntimes, self.ntotal, 8), cfdtype)

    def finalize(self):
        sd = self.data[0, :, 0].real
        i_sd_zero = np.where(sd != 0.0)[0]
        i_node_zero = np.where(self.element_node[:, 1] != 0)[0]
        assert i_node_zero.max() > 0, 'CBEAM element_node hasnt been filled'
        i = np.union1d(i_sd_zero, i_node_zero)
        self.element = self.element[i]
        self.element_node = self.element_node[i, :]
        self.data = self.data[:, i, :]

    def build_dataframe(self):
        """creates a pandas dataframe"""
        # Freq                                          0.00001             10.00000  ...            50.00000            60.00000
        # ElementID Location Item                                                     ...
        # 12.0      12.0     bending_moment1  0.000000+0.000000j  0.000000+0.000000j  ...  0.000000+0.000000j  0.000000+0.000000j
        #                    bending_moment2  0.000000+0.000000j  0.000000+0.000000j  ...  0.000000+0.000000j  0.000000+0.000000j
        #                    shear1           0.000000+0.000000j  0.000000+0.000000j  ...  0.000000+0.000000j  0.000000+0.000000j
        #                    shear2           0.000000+0.000000j  0.000000+0.000000j  ...  0.000000+0.000000j  0.000000+0.000000j
        #                    axial_force      0.000000+0.000000j  0.000000+0.000000j  ...  0.000000+0.000000j  0.000000+0.000000j
        #                    total_torque     0.000000+0.000000j  0.000000+0.000000j  ...  0.000000+0.000000j  0.000000+0.000000j
        #                    warping_torque   0.000000+0.000000j  0.000000+0.000000j  ...  0.000000+0.000000j  0.000000+0.000000j
        # 0.0       1.0      bending_moment1  0.000000+0.000000j  0.000000+0.000000j  ...  0.000000+0.000000j  0.000000+0.000000j
        #                    bending_moment2  0.000000+0.000000j  0.000000+0.000000j  ...  0.000000+0.000000j  0.000000+0.000000j
        #                    shear1           0.000000+0.000000j  0.000000+0.000000j  ...  0.000000+0.000000j  0.000000+0.000000j
        #                    shear2           0.000000+0.000000j  0.000000+0.000000j  ...  0.000000+0.000000j  0.000000+0.000000j
        #                    axial_force      0.000000+0.000000j  0.000000+0.000000j  ...  0.000000+0.000000j  0.000000+0.000000j
        #                    total_torque     0.000000+0.000000j  0.000000+0.000000j  ...  0.000000+0.000000j  0.000000+0.000000j
        #                    warping_torque   0.000000+0.000000j  0.000000+0.000000j  ...  0.000000+0.000000j  0.000000+0.000000j
        import pandas as pd
        headers = self.get_headers()[1:]
        column_names, column_values = self._build_dataframe_transient_header()
        element_location = [
            self.element_node[:, 0],
            self.data[0, :, 0].real,
        ]
        is_v25 = pd.__version__ >= '0.25'
        if is_v25:
            print(f'skipping pandas {self.class_name}')
            return
        # wrong type for ElementID
        #data_frame = self._build_pandas_transient_element_node(
            #column_values, column_names,
            #headers, element_location, self.data[:, :, 1:])
        #data_frame.index.names = ['ElementID', 'Location', 'Item']
        #data_frame.index['ElementID', :]# .astype('int32')
        #print(data_frame)

        data_frame = pd.Panel(self.data[:, :, 1:], items=column_values,
                              major_axis=element_location, minor_axis=headers).to_frame()
        data_frame.columns.names = column_names
        data_frame.index.names = ['ElementID', 'Location', 'Item']
        #print(data_frame)
        self.data_frame = data_frame

    def __eq__(self, table):  # pragma: no cover
        return self.assert_equal(table)

    def assert_equal(self, table, rtol=1.e-5, atol=1.e-8):
        assert self.is_sort1 == table.is_sort1
        self._eq_header(table)
        if not np.allclose(self.data, table.data, atol=atol):
            msg = 'table_name=%r class_name=%s\n' % (self.table_name, self.__class__.__name__)
            msg += '%s\n' % str(self.code_information())
            ntimes = self.data.shape[0]

            i = 0
            if self.is_sort1:
                for itime in range(ntimes):
                    for ieid, eid in enumerate(self.element):
                        t1 = self.data[itime, ieid, :]
                        t2 = table.data[itime, ieid, :]
                        #print(t1)
                        #'sd', 'bending_moment1', 'bending_moment2', 'shear1', 'shear2',
                        #'axial_force', 'total_torque', 'warping_torque', ]

                        (sd1, bm11, bm21, shear11, shear21, axial1, total_torque1, warp_torque1) = t1
                        (sd2, bm12, bm22, shear12, shear22, axial2, total_torque2, warp_torque2) = t2
                        d = t1 - t2
                        if not allclose(t1, t2, atol=atol):
                            msg += (
                                '%-4s  (%s, %sj, %s, %sj, %s, %sj, %s, %sj, %s, %sj, %s, %sj, %s, %sj)\n'
                                '      (%s, %sj, %s, %sj, %s, %sj, %s, %sj, %s, %sj, %s, %sj, %s, %sj)\n'
                                '  dt12=(%s, %sj, %s, %sj, %s, %sj, %s, %sj, %s, %sj, %s, %sj, %s, %sj)\n' % (
                                    eid,
                                    bm11.real, bm11.imag,
                                    bm21.real, bm21.imag,
                                    shear11.real, shear11.imag,
                                    shear21.real, shear21.imag,
                                    axial1.real, axial1.imag,
                                    total_torque1.real, total_torque1.imag,
                                    warp_torque1.real, warp_torque1.imag,

                                    bm12.real, bm12.imag,
                                    bm22.real, bm22.imag,
                                    shear12.real, shear12.imag,
                                    shear22.real, shear22.imag,
                                    axial2.real, axial2.imag,
                                    total_torque2.real, total_torque2.imag,
                                    warp_torque2.real, warp_torque2.imag,

                                    d[0].real, d[0].imag,
                                    d[1].real, d[1].imag,
                                    d[2].real, d[2].imag,
                                    d[3].real, d[3].imag,
                                    d[4].real, d[4].imag,
                                    d[5].real, d[5].imag,
                                    d[6].real, d[6].imag,
                                ))
                        i += 1
                    if i > 10:
                        print(msg)
                        raise ValueError(msg)
            else:
                raise NotImplementedError(self.is_sort2)
            if i > 0:
                print(msg)
                raise ValueError(msg)
        return True

    #def add_new_element_sort1(self, dt, eid, nid, sd, bm1, bm2, ts1, ts2, af, ttrq, wtrq):
        #return self.add_sort1(dt, eid, nid, sd, bm1, bm2, ts1, ts2, af, ttrq, wtrq)

    def add_sort1(self, dt, eid, nid, sd, bm1, bm2, ts1, ts2, af, ttrq, wtrq):
        """unvectorized method for adding SORT1 transient data"""
        assert isinstance(eid, integer_types) and eid > 0, 'dt=%s eid=%s type=%s' % (dt, eid, type(eid))
        self._times[self.itime] = dt
        self.data[self.itime, self.itotal, :] = [sd, bm1, bm2, ts1, ts2, af, ttrq, wtrq]
        self.element[self.itotal] = eid
        self.element_node[self.itotal, :] = [eid, nid]
        self.itotal += 1

    def get_stats(self, short=False) -> List[str]:
        if not self.is_built:
            return [
                '<%s>\n' % self.__class__.__name__,
                '  ntimes: %i\n' % self.ntimes,
                '  ntotal: %i\n' % self.ntotal,
            ]

        nelements = self.nelements
        ntimes = self.ntimes
        msg = []

        if self.nonlinear_factor not in (None, np.nan):  # transient
            msg.append('  type=%s ntimes=%i nelements=%i\n'
                       % (self.__class__.__name__, ntimes, nelements))
        else:
            msg.append('  type=%s nelements=%i\n' % (self.__class__.__name__, nelements))
        #msg.append('  eType, cid\n')
        msg.append('  data: [ntimes, nelements, 8] where 8=[%s]\n' % str(', '.join(self.get_headers())))
        msg.append('  data.shape = %s\n' % str(self.data.shape).replace('L', ''))
        msg.append('  is_sort1=%s is_sort2=%s\n' % (self.is_sort1, self.is_sort2))
        msg.append('  CBEAM\n')
        msg += self.get_data_code()
        return msg

    def write_f06(self, f06_file, header=None, page_stamp='PAGE %s', page_num=1, is_mag_phase=False, is_sort1=True):
        if header is None:
            header = []
        # option B
        #'                           C O M P L E X   F O R C E S   I N   B E A M   E L E M E N T S   ( C B E A M ) '
        #'                                                          (REAL/IMAGINARY)'
        #'                    STAT DIST/   - BENDING MOMENTS -            - WEB  SHEARS -           AXIAL          TOTAL          WARPING'
        #'   ELEMENT-ID  GRID   LENGTH    PLANE 1       PLANE 2        PLANE 1       PLANE 2        FORCE          TORQUE         TORQUE'
        #'0        20'
        #'0                11   0.000    0.0           0.0            0.0           0.0            0.0            0.0            0.0'
        #'                               0.0           0.0            0.0           0.0            0.0            0.0            0.0'
        #'0                12   1.000    0.0           0.0            0.0           0.0            0.0            0.0            0.0'
        #'                               0.0           0.0            0.0           0.0            0.0            0.0            0.0'

        #msg_temp, nnodes = get_f06_header(self, is_mag_phase, is_sort1)
        #print('write_f06 not implemented for ComplexCBeamForceArray')
        #return page_num

        #is_sort1 = False
        if is_mag_phase:
            mag_phase = '                                                          (MAGNITUDE/PHASE)\n \n'
        else:
            mag_phase = '                                                          (REAL/IMAGINARY)\n \n'


        if is_sort1:
            line1 = '0    ELEMENT         BEND-MOMENT-END-A            BEND-MOMENT-END-B                  SHEAR\n'
            line2 = '       ID.         PLANE 1       PLANE 2        PLANE 1       PLANE 2        PLANE 1       PLANE 2        FORCE          TORQUE\n'
        else:
            name = self.data_code['name']
            if name == 'freq':
                name = 'FREQUENCY'
            else: # mode
                raise RuntimeError(name)
            line1 = '                    BEND-MOMENT-END-A            BEND-MOMENT-END-B                  SHEAR\n'
            line2 = '   %16s       PLANE 1       PLANE 2        PLANE 1       PLANE 2        PLANE 1       PLANE 2        FORCE          TORQUE\n' % name

        # force
        msg_temp = header + [
            '                             C O M P L E X   F O R C E S   I N   B A R   E L E M E N T S   ( C B E A M )\n',
            mag_phase,
            ' ',
            line1,
            line2,
        ]
        if self.is_sort1:
            assert self.is_sort1 == True, str(self)
            #if is_sort1:
            page_num = self._write_sort1_as_sort1(f06_file, page_num, page_stamp, header, msg_temp, is_mag_phase)
            #else:
                #self._write_sort1_as_sort2(f06_file, page_num, page_stamp, header, msg_temp, is_mag_phase)
        else:
            assert self.is_sort1 == True, str(self)
        return page_num - 1

    def _write_sort1_as_sort1(self, f06_file, page_num, page_stamp, header, msg_temp, is_mag_phase):
        eids = self.element
        #times = self._times
        ntimes = self.data.shape[0]
        for itime in range(ntimes):
            dt = self._times[itime]
            dt_line = ' %14s = %12.5E\n' % (self.data_code['name'], dt)
            header[1] = dt_line
            msg = header + msg_temp
            f06_file.write(''.join(msg))

            #bm1a, bm2a, bm1b, bm2b, ts1, ts2, af, trq
            assert self.is_sort1 == True, str(self)
            #sd, bm1, bm2, ts1, ts2, af, ttrq, wtrq
            sd = self.data[itime, :, 0]
            bm1 = self.data[itime, :, 1]
            bm2 = self.data[itime, :, 2]
            ts1 = self.data[itime, :, 3]
            ts2 = self.data[itime, :, 4]
            af = self.data[itime, :, 5]
            ttrq = self.data[itime, :, 6]
            wtrq = self.data[itime, :, 7]

            for eid, sdi, bm1i, bm2i, ts1i, ts2i, afi, ttrqi, wtrqi in zip(eids, sd, bm1, bm2, ts1, ts2, af, ttrq, wtrq):
                vals = (sdi, bm1i, bm2i, ts1i, ts2i, afi, ttrqi, wtrqi)
                vals2 = write_imag_floats_13e(vals, is_mag_phase)
                (sdir, bm1ir, bm2ir, ts1ir, ts2ir, afir, ttrqir, wtrqir,
                 sdii, bm1ii, bm2ii, ts1ii, ts2ii, afii, ttrqii, wtrqii) = vals2

                f06_file.write('0%16i   %-13s  %-13s  %-13s  %-13s  %-13s  %-13s  %-13s  %s\n'
                               ' %14s   %-13s  %-13s  %-13s  %-13s  %-13s  %-13s  %-13s  %s\n' % (
                                   eid, sdir, bm1ir, bm2ir, ts1ir, ts2ir, afir, ttrqir, wtrqir,
                                   '', sdii, bm1ii, bm2ii, ts1ii, ts2ii, afii, ttrqii, wtrqii))
            f06_file.write(page_stamp % page_num)
            page_num += 1
        return page_num

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

        eids = self.element_node[:, 0]
        nids = self.element_node[:, 1]
        #long_form = False
        #if nids.min() == 0:
            #long_form = True

        eids_device = eids * 10 + self.device_code
        ueids = np.unique(eids)
        #ieid = np.searchsorted(eids, ueids)
        # table 4 info
        #ntimes = self.data.shape[0]
        #nnodes = self.data.shape[1]
        nelements = len(ueids)

        # 21 = 1 node, 3 principal, 6 components, 9 vectors, 2 p/ovm
        #ntotal = ((nnodes * 21) + 1) + (nelements * 4)

        ntotali = self.num_wide
        ntotal = ntotali * nelements

        op2_ascii.write('  ntimes = %s\n' % self.ntimes)

        #print('ntotal=%s' % (ntotal))
        #assert ntotal == 193, ntotal

        if self.is_sort1:
            struct1 = Struct(endian + b'2i 15f')
            struct2 = Struct(endian + b'i 15f')
        else:
            raise NotImplementedError('SORT2')

        op2_ascii.write('nelements=%i\n' % nelements)
        for itime in range(self.ntimes):
            self._write_table_3(op2, op2_ascii, new_result, itable, itime)

            # record 4
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

            sd = self.data[itime, :, 0]
            bm1 = self.data[itime, :, 1]
            bm2 = self.data[itime, :, 2]
            ts1 = self.data[itime, :, 3]
            ts2 = self.data[itime, :, 4]
            af = self.data[itime, :, 5]
            ttrq = self.data[itime, :, 6]
            wtrq = self.data[itime, :, 7]

            icount = 0
            nwide = 0
            ielement = 0
            assert len(eids) == len(sd)
            for eid, sdi, bm1i, bm2i, ts1i, ts2i, afi, ttrqi, wtrqi in zip(eids, sd, bm1, bm2, ts1, ts2, af, ttrq, wtrq):
                if icount == 0:
                    eid_device = eids_device[ielement]
                    nid = nids[ielement]
                    data = [eid_device, nid, sdi.real,
                            bm1i.real, bm2i.real, ts1i.real, ts2i.real, afi.real, ttrqi.real, wtrqi.real,
                            bm1i.imag, bm2i.imag, ts1i.imag, ts2i.imag, afi.imag, ttrqi.imag, wtrqi.imag] # 17
                    op2.write(struct1.pack(*data))
                    ielement += 1
                    icount = 1
                elif nid > 0 and icount > 0:
                    # 11 total nodes, with 1, 11 getting an nid; the other 9 being
                    # xxb sections
                    data = [0, 0.,
                            0., 0., 0., 0., 0., 0., 0.,
                            0., 0., 0., 0., 0., 0., 0.]
                    #print('***adding %s\n' % (10-icount))
                    for unused_i in range(10 - icount):
                        op2.write(struct2.pack(*data))
                        nwide += len(data)

                    eid_device2 = eids_device[ielement]
                    #print(eids_device)
                    assert eid_device == eid_device2, 'eid_device=%s eid_device2=%s' % (eid_device, eid_device2)
                    nid = nids[ielement]
                    data = [nid, sdi.real,
                            bm1i.real, bm2i.real, ts1i.real, ts2i.real, afi.real, ttrqi.real, wtrqi.real,
                            bm1i.imag, bm2i.imag, ts1i.imag, ts2i.imag, afi.imag, ttrqi.imag, wtrqi.imag] # 16
                    op2.write(struct2.pack(*data))
                    ielement += 1
                    icount = 0
                else:
                    raise RuntimeError('CBEAM OEF op2 writer')
                    #data = [0, xxb, sxc, sxd, sxe, sxf, smax, smin, smt, smc]  # 10
                    #op2.write(struct2.pack(*data))
                    #icount += 1

                op2_ascii.write('  eid_device=%s data=%s\n' % (eid_device, str(data)))
                nwide += len(data)

            assert ntotal == nwide, 'ntotal=%s nwide=%s' % (ntotal, nwide)

            itable -= 1
            header = [4 * ntotal,]
            op2.write(pack('i', *header))
            op2_ascii.write('footer = %s\n' % header)
            new_result = False
        return itable


class ComplexCBendForceArray(BaseElement):  # 69-CBEND
    def __init__(self, data_code, is_sort1, isubcase, dt):
        self.element_type = None
        self.element_name = None
        BaseElement.__init__(self, data_code, isubcase)
        #self.code = [self.format_code, self.sort_code, self.s_code]

        #self.ntimes = 0  # or frequency/mode
        #self.ntotal = 0
        self.nelements = 0  # result specific

        #if is_sort1:
            #pass
        #else:
            #raise NotImplementedError('SORT2')

    @property
    def is_real(self):
        """is the result real?"""
        return False

    @property
    def is_complex(self):
        """is the result complex?"""
        return True

    def _reset_indices(self):
        self.itotal = 0
        self.ielement = 0

    def get_headers(self) -> List[str]:
        headers = [
            'bending_moment_1a', 'bending_moment_2a', 'shear_1a', 'shear_2a', 'axial_a', 'torque_a',
            'bending_moment_1b', 'bending_moment_2b', 'shear_1b', 'shear_2b', 'axial_b', 'torque_b',
        ]
        return headers

    def build(self):
        """sizes the vectorized attributes of the ComplexCBendForceArray"""
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
        dtype = 'float32'
        if isinstance(self.nonlinear_factor, integer_types):
            dtype = 'int32'
        self._times = zeros(self.ntimes, dtype=dtype)
        self.element_node = zeros((self.nelements, 3), dtype='int32')

        #[bending_moment_1a, bending_moment_2a, shear_1a, shear_2a, axial_a, torque_a
        # bending_moment_1b, bending_moment_2b, shear_1b, shear_2b, axial_b, torque_b]
        self.data = zeros((self.ntimes, self.nelements, 12), dtype='complex64')

    def build_dataframe(self):
        """creates a pandas dataframe"""
        #Freq                                       0.0                 2.5
        #ElementID Item
        #6901      bending_moment_1a  1.066567-0.035549j  1.066996-0.035577j
        #          bending_moment_2a  1.101375-0.036709j  1.102188-0.036763j
        #          shear_1a           0.516478-0.017214j  0.516842-0.017239j
        #          shear_2a           0.859292-0.028640j  0.860111-0.028695j
        #          axial_a            0.834822-0.027825j  0.834982-0.027835j
        #          torque_a           0.953420-0.031777j  0.953947-0.031813j
        #          bending_moment_1b -0.284733+0.009490j -0.284828+0.009497j
        #          bending_moment_2b  0.094127-0.003137j  0.093836-0.003118j
        #          shear_1b           0.834822-0.027825j  0.834982-0.027835j
        #          shear_2b           0.859292-0.028640j  0.860111-0.028695j
        #          axial_b           -0.516478+0.017214j -0.516842+0.017239j
        #          torque_b          -0.242082+0.008069j -0.242077+0.008068j
        #6902      bending_moment_1a -0.931214+0.031037j -0.931519+0.031058j
        headers = self.get_headers()
        column_names, column_values = self._build_dataframe_transient_header()

        # element_node is (nelements, 3)
        element = self.element_node[:, 0]
        self.data_frame = self._build_pandas_transient_elements(
            column_values, column_names,
            headers, element, self.data)

    def __eq__(self, table):  # pragma: no cover
        assert self.is_sort1 == table.is_sort1
        if not np.array_equal(self.element_node, table.element_node):
            assert self.element_node.shape == table.element_node.shape, 'element_node shape=%s table.shape=%s' % (self.element_node.shape, table.element_nodes.shape)
            msg = 'table_name=%r class_name=%s\n' % (self.table_name, self.__class__.__name__)
            msg += '%s\n' % str(self.code_information())
            msg += 'Eid, Nid_A, Nid_B\n'
            for (eid1, nida1, nidb1), (eid2, nida2, nidb2) in zip(self.element_node, table.element_node):
                msg += '(%s, %s, %s), (%s, %s, %s)\n' % (eid1, nida1, nidb1, eid2, nida2, nidb2)
            print(msg)
            raise ValueError(msg)
        if not np.array_equal(self.data, table.data):
            msg = 'table_name=%r class_name=%s\n' % (self.table_name, self.__class__.__name__)
            msg += '%s\n' % str(self.code_information())
            i = 0
            eids = self.element_node[:, 0]
            for itime in range(self.ntimes):
                for ie, eid in enumerate(eids):
                    t1 = self.data[itime, ie, :]
                    t2 = table.data[itime, ie, :]
                    (bending_moment_1a1, bending_moment_2a1, shear_1a1, shear_2a1, axial_a1, torque_a1,
                     bending_moment_1b1, bending_moment_2b1, shear_1b1, shear_2b1, axial_b1, torque_b1) = t1
                    (bending_moment_1a2, bending_moment_2a2, shear_1a2, shear_2a2, axial_a2, torque_a2,
                     bending_moment_1b2, bending_moment_2b2, shear_1b2, shear_2b2, axial_b2, torque_b2) = t2

                    if not allclose(t1, t2):
                        msg += '(%s)    (%s, %s)  (%s, %s)\n' % (
                            eid,
                            bending_moment_1a1.real,
                            bending_moment_1b1.real,

                            bending_moment_1a2.real,
                            bending_moment_1b2.real, )
                        i += 1
                        if i > 10:
                            print(msg)
                            raise ValueError(msg)

                    #if not allclose(t1, t2):
                        #msg += '(%s)    (%s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s)  (%s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s)\n' % (
                            #eid,
                            #bending_moment_1a1, bending_moment_2a1, shear_1a1, shear_2a1, axial_a1, torque_a1,
                            #bending_moment_1b1, bending_moment_2b1, shear_1b1, shear_2b1, axial_b1, torque_b1,

                            #bending_moment_1a2, bending_moment_2a2, shear_1a2, shear_2a2, axial_a2, torque_a2,
                            #bending_moment_1b2, bending_moment_2b2, shear_1b2, shear_2b2, axial_b2, torque_b2)
                        #i += 1
                        #if i > 10:
                            #print(msg)
                            #raise ValueError(msg)
                #print(msg)
                if i > 0:
                    raise ValueError(msg)
        return True

    def add_sort1(self, dt, eid,
                  nid_a, bending_moment_1a, bending_moment_2a, shear_1a, shear_2a, axial_a, torque_a,
                  nid_b, bending_moment_1b, bending_moment_2b, shear_1b, shear_2b, axial_b, torque_b):
        """unvectorized method for adding SORT1 transient data"""
        assert isinstance(eid, integer_types) and eid > 0, 'dt=%s eid=%s' % (dt, eid)
        #bending_moment_1a, bending_moment_2a, shear_1a, shear_2a, axial_a, torque_a,
        #bending_moment_1b, bending_moment_2b, shear_1b, shear_2b, axial_b, torque_b

        self._times[self.itime] = dt
        self.element_node[self.ielement] = [eid, nid_a, nid_b]
        self.data[self.itime, self.ielement, :] = [
            bending_moment_1a, bending_moment_2a, shear_1a, shear_2a, axial_a, torque_a,
            bending_moment_1b, bending_moment_2b, shear_1b, shear_2b, axial_b, torque_b
        ]
        self.ielement += 1
        if self.ielement == self.nelements:
            self.ielement = 0

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

        msg = []
        if self.nonlinear_factor not in (None, np.nan):  # transient
            msg.append('  type=%s ntimes=%i nelements=%i; table_name=%r\n'
                       % (self.__class__.__name__, ntimes, nelements, self.table_name))
            ntimes_word = 'ntimes'
        else:
            msg.append('  type=%s nelements=%i; table_name=%r\n'
                       % (self.__class__.__name__, nelements, self.table_name))
            ntimes_word = '1'
        msg.append('  eType\n')
        headers = self.get_headers()
        n = len(headers)
        msg.append('  data: [%s, nnodes, %i] where %i=[%s]\n' % (
            ntimes_word, n, n, str(', '.join(headers))))
        msg.append('  data.shape = %s\n' % str(self.data.shape).replace('L', ''))
        #msg.append('  element type: %s\n' % self.element_type)
        msg.append('  element name: %s\n' % self.element_name)
        msg += self.get_data_code()
        return msg

    def get_f06_header(self, is_mag_phase=True, is_sort1=True):
        msg = ['                           C O M P L E X   F O R C E S   I N   B E N D    E L E M E N T S   ( C B E N D )\n']

        if is_mag_phase:
            msg += ['                                                          (MAGNITUDE/PHASE)\n']
        else:
            msg += ['                                                          (REAL/IMAGINARY)\n']

        if is_sort1:
            msg += [
                '                                 - BENDING MOMENTS -            -   SHEARS   -            AXIAL'
                '   ELEMENT-ID  GRID    END      PLANE 1       PLANE 2        PLANE 1       PLANE 2        FORCE          TORQUE'
            ]
        else:
            raise NotImplementedError('sort2')
        return msg

    def write_f06(self, f06_file, header=None, page_stamp='PAGE %s',
                  page_num=1, is_mag_phase=False, is_sort1=True):
        if header is None:
            header = []
        #'                           C O M P L E X   F O R C E S   I N   B E N D    E L E M E N T S   ( C B E N D )'
        #'                                                          (REAL/IMAGINARY)'
        #'                                 - BENDING MOMENTS -            -   SHEARS   -            AXIAL'
        #'   ELEMENT-ID  GRID    END      PLANE 1       PLANE 2        PLANE 1       PLANE 2        FORCE          TORQUE'
        #'0        27      21      A     0.0           0.0            0.0           0.0            0.0            0.0'
        #'                               0.0           0.0            0.0           0.0            0.0            0.0'
        #'0                22      B     0.0           0.0            0.0           0.0            0.0            0.0'
        #'                               0.0           0.0            0.0           0.0            0.0            0.0'
        msg_temp = self.get_f06_header(is_mag_phase=is_mag_phase, is_sort1=is_sort1)

        # write the f06
        #(ntimes, ntotal, two) = self.data.shape
        ntimes = self.data.shape[0]

        eids = self.element_node[:, 0]
        nid_a = self.element_node[:, 1]
        nid_b = self.element_node[:, 2]
        for itime in range(ntimes):
            dt = self._times[itime]  # TODO: rename this...
            header = _eigenvalue_header(self, header, itime, ntimes, dt)
            f06_file.write(''.join(header + msg_temp))

            #print("self.data.shape=%s itime=%s ieids=%s" % (str(self.data.shape), itime, str(ieids)))
            bending_moment_1a = self.data[itime, :, 0]
            bending_moment_2a = self.data[itime, :, 1]
            shear_1a = self.data[itime, :, 2]
            shear_2a = self.data[itime, :, 3]
            axial_a = self.data[itime, :, 4]
            torque_a = self.data[itime, :, 5]
            bending_moment_1b = self.data[itime, :, 6]
            bending_moment_2b = self.data[itime, :, 7]
            shear_1b = self.data[itime, :, 8]
            shear_2b = self.data[itime, :, 9]
            axial_b = self.data[itime, :, 10]
            torque_b = self.data[itime, :, 11]

            for (eid,
                 nid_ai, bending_moment_1ai, bending_moment_2ai, shear_1ai, shear_2ai, axial_ai, torque_ai,
                 nid_bi, bending_moment_1bi, bending_moment_2bi, shear_1bi, shear_2bi, axial_bi, torque_bi) in zip(eids,
                                                                                                                   nid_a, bending_moment_1a, bending_moment_2a, shear_1a, shear_2a, axial_a, torque_a,
                                                                                                                   nid_b, bending_moment_1b, bending_moment_2b, shear_1b, shear_2b, axial_b, torque_b):
                [bending_moment_1ari, bending_moment_2ari, shear_1ari, shear_2ari, axial_ari, torque_ari,
                 bending_moment_1bri, bending_moment_2bri, shear_1bri, shear_2bri, axial_bri, torque_bri,
                 bending_moment_1aii, bending_moment_2aii, shear_1aii, shear_2aii, axial_aii, torque_aii,
                 bending_moment_1bii, bending_moment_2bii, shear_1bii, shear_2bii, axial_bii, torque_bii,
                 ] = write_imag_floats_13e(
                     [bending_moment_1ai, bending_moment_2ai, shear_1ai, shear_2ai, axial_ai, torque_ai,
                      bending_moment_1bi, bending_moment_2bi, shear_1bi, shear_2bi, axial_bi, torque_bi],
                     is_mag_phase)
                f06_file.write(
                    '0  %8s%8s      A    %13s %13s  %13s %13s  %13s  %s\n'
                    '                              %13s %13s  %13s %13s  %13s  %s\n'
                    '0  %8s%8s      B    %13s %13s  %13s %13s  %13s  %s\n'
                    '                              %13s %13s  %13s %13s  %13s  %s\n'
                    % (
                        eid, nid_ai,
                        bending_moment_1ari, bending_moment_2ari, shear_1ari, shear_2ari, axial_ari, torque_ari,
                        bending_moment_1aii, bending_moment_2aii, shear_1aii, shear_2aii, axial_aii, torque_aii,
                        '', nid_bi,
                        bending_moment_1bri, bending_moment_2bri, shear_1bri, shear_2bri, axial_bri, torque_bri,
                        bending_moment_1bii, bending_moment_2bii, shear_1bii, shear_2bii, axial_bii, torque_bii,))
            f06_file.write(page_stamp % page_num)
            page_num += 1
        return page_num - 1


class ComplexSolidPressureForceArray(ComplexForceObject):
    def __init__(self, data_code, is_sort1, isubcase, dt):
        self.element_type = None
        self.element_name = None
        ComplexForceObject.__init__(self, data_code, isubcase)
        #self.code = [self.format_code, self.sort_code, self.s_code]

        #self.ntimes = 0  # or frequency/mode
        #self.ntotal = 0
        self.nelements = 0  # result specific

        if is_sort1:
            pass
        else:
            raise NotImplementedError('SORT2')

    def _reset_indices(self):
        self.itotal = 0
        self.ielement = 0

    def get_headers(self) -> List[str]:
        headers = ['ax', 'ay', 'az', 'vx', 'vy', 'vz', 'pressure']
        return headers

    #def get_headers(self):
        #headers = ['axial', 'torque']
        #return headers

    def build(self):
        """sizes the vectorized attributes of the ComplexSolidPressureForceArray"""
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
        dtype, idtype, cfdtype = get_complex_times_dtype(self.nonlinear_factor, self.size)
        self._times = zeros(self.ntimes, dtype=dtype)
        self.element = zeros(self.nelements, dtype=idtype)

        #[ax, ay, az, vx, vy, vz, pressure]
        self.data = zeros((self.ntimes, self.ntotal, 7), dtype=cfdtype)

    def build_dataframe(self):
        """creates a pandas dataframe"""
        #Mode                                        1                           2
        #EigenvalueReal                    -0.000000                   -0.000000
        #EigenvalueImag                    -0.000000                   -0.000000
        #Damping                               0.000000                    0.000000
        #ElementID Item
        #1000      ax       -1.887379e-13+2.791559e-13j -1.901257e-13+2.789015e-13j
        #          ay        3.330669e-14-7.316397e-14j  1.776357e-14-7.368508e-14j
        #          az       -1.360023e-13-9.545406e-14j -1.432188e-13-8.333307e-14j
        #          vx        0.000000e+00+0.000000e+00j  0.000000e+00+0.000000e+00j
        #          vy        0.000000e+00+0.000000e+00j  0.000000e+00+0.000000e+00j
        #          vz        0.000000e+00+0.000000e+00j  0.000000e+00+0.000000e+00j
        #          pressure  0.000000e+00+0.000000e+00j  0.000000e+00+0.000000e+00j
        headers = self.get_headers()
        column_names, column_values = self._build_dataframe_transient_header()
        self.data_frame = self._build_pandas_transient_elements(
            column_values, column_names,
            headers, self.element, self.data)

    def __eq__(self, table):  # pragma: no cover
        self._eq_header(table)
        assert self.is_sort1 == table.is_sort1
        if not np.array_equal(self.data, table.data):
            msg = 'table_name=%r class_name=%s\n' % (self.table_name, self.__class__.__name__)
            msg += '%s\n' % str(self.code_information())
            i = 0
            for itime in range(self.ntimes):
                for ie, eid in enumerate(self.element):
                    t1 = self.data[itime, ie, :]
                    t2 = table.data[itime, ie, :]
                    (ax1, ay1, az1, vx1, vy1, vz1, pressure1) = t1
                    (ax2, ay2, az2, vx2, vy2, vz2, pressure2) = t2
                    #rpressure1 = pressure1.real
                    #rpressure2 = pressure2.real
                    if not allclose([ax1, ay1, az1, vx1, vy1, vz1],
                                    [ax2, ay2, az2, vx2, vy2, vz2]):
                        msg += '%s    (%s, %s)  (%s, %s)\n' % (
                            eid,
                            ax1.real, t1.imag,
                            ax2.real, t2.imag)
                        i += 1
                        if i > 10:
                            print(msg)
                            raise ValueError(msg)
                #print(msg)
                if i > 0:
                    raise ValueError(msg)
        return True

    def add_sort1(self, dt, eid, ename, ax, ay, az, vx, vy, vz, pressure):
        """unvectorized method for adding SORT1 transient data"""
        assert isinstance(eid, integer_types) and eid > 0, 'dt=%s eid=%s' % (dt, eid)
        self._times[self.itime] = dt
        self.element[self.ielement] = eid
        self.data[self.itime, self.ielement, :] = [ax, ay, az, vx, vy, vz, pressure]
        self.ielement += 1

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

        msg = []
        if self.nonlinear_factor not in (None, np.nan):  # transient
            msg.append('  type=%s ntimes=%i nelements=%i; table_name=%r\n'
                       % (self.__class__.__name__, ntimes, nelements, self.table_name))
            ntimes_word = 'ntimes'
        else:
            msg.append('  type=%s nelements=%i; table_name=%r\n'
                       % (self.__class__.__name__, nelements, self.table_name))
            ntimes_word = '1'
        msg.append('  eType\n')
        headers = self.get_headers()
        n = len(headers)
        msg.append('  data: [%s, nelements, %i] where %i=[%s]\n' % (
            ntimes_word, n, n, str(', '.join(headers))))
        msg.append('  data.shape = %s\n' % str(self.data.shape).replace('L', ''))
        #msg.append('  element type: %s\n' % self.element_type)
        msg.append('  element name: %s\n' % self.element_name)
        msg += self.get_data_code()
        return msg

    def get_f06_header(self, is_mag_phase=True, is_sort1=True):
        msg = [
            '                                                 ( R O O T   M E A N   S Q U A R E )      \n'
            '              C O M P L E X   A C C E L E R A T I O N S   V E L O C I T I E S   A N D   P R E S S U R E   L E V E L S\n'
            #'                                                          (REAL/IMAGINARY)'
            #'     ELE-ID   EL-TYPE    X-ACCELERATION  Y-ACCELERATION  Z-ACCELERATION     X-VELOCITY     Y-VELOCITY     Z-VELOCITY   PRESSURE (DB)'
            #'       2000    PENPR      6.883253E+06    1.066544E+07   -6.883253E+06     7.288279E+05  -3.134843E+04  -7.288279E+05   1.162309E+02'
            #'                          1.831744E+07   -7.878719E+05   -1.831744E+07    -2.738759E+05  -4.243642E+05   2.738759E+05'
            #''
        ]
        #msg = ['                                   C O M P L E X   A C O U S T I C   P R E S S U R E   R E S U L T S']
        #'                                   C O M P L E X   A C O U S T I C   P R E S S U R E   R E S U L T S'
        #'                                                         (MAGNITUDE/PHASE)'
        #' '
        #'      POINT ID.               TYPE                    P                   P(RMS)              DB                  DB(A)'
        #'0           57                 S                   7.339671E+05        5.189931E+05        1.173135E+02        3.011353E+01'
        #'                                                       249.9102            249.9102            249.9102            249.9102'

        if is_mag_phase:
            msg += ['                                                          (MAGNITUDE/PHASE)\n \n']
        else:
            msg += ['                                                          (REAL/IMAGINARY)\n \n']

        if is_sort1:
            msg += ['     ELE-ID   EL-TYPE    X-ACCELERATION  Y-ACCELERATION  Z-ACCELERATION     X-VELOCITY     Y-VELOCITY     Z-VELOCITY   PRESSURE (DB)\n']
            #msg += [
                #'      POINT ID.               TYPE                    P                   P(RMS)              DB                  DB(A)\n'
            #]
            #'                      14                  0.0          /  0.0                           0.0          /  0.0'
        else:
            raise NotImplementedError('sort2')


        return msg

    #def get_element_index(self, eids):
        ## elements are always sorted; nodes are not
        #itot = searchsorted(eids, self.element)  #[0]
        #return itot

    #def eid_to_element_node_index(self, eids):
        ##ind = ravel([searchsorted(self.element == eid) for eid in eids])
        #ind = searchsorted(eids, self.element)
        ##ind = ind.reshape(ind.size)
        ##ind.sort()
        #return ind

    def write_f06(self, f06_file, header=None, page_stamp='PAGE %s',
                  page_num=1, is_mag_phase=False, is_sort1=True):
        if header is None:
            header = []
        msg_temp = self.get_f06_header(is_mag_phase=is_mag_phase, is_sort1=is_sort1)

        # write the f06
        #(ntimes, ntotal, two) = self.data.shape
        ntimes = self.data.shape[0]

        eids = self.element

        #print('len(eids)=%s nwrite=%s is_odd=%s' % (len(eids), nwrite, is_odd))
        etypei = self.element_type
        for itime in range(ntimes):
            dt = self._times[itime]  # TODO: rename this...
            header = _eigenvalue_header(self, header, itime, ntimes, dt)
            f06_file.write(''.join(header + msg_temp))

            #print("self.data.shape=%s itime=%s ieids=%s" % (str(self.data.shape), itime, str(ieids)))
            ax = self.data[itime, :, 0]
            ay = self.data[itime, :, 0]
            az = self.data[itime, :, 0]
            vx = self.data[itime, :, 0]
            vy = self.data[itime, :, 0]
            vz = self.data[itime, :, 0]
            pressure = self.data[itime, :, 0]

            for eid, axi, ayi, azi, vxi, vyi, vzi, pressurei in zip(eids, ax, ay, az, vx, vy, vz, pressure):
                out = write_imag_floats_13e([axi, ayi, azi, vxi, vyi, vzi, pressurei], is_mag_phase)
                [saxr, sayr, sazr, svxr, svyr, svzr, spressurer,
                 saxi, sayi, sazi, svxi, svyi, svzi, spressurei] = out
                #'       1000    HEXPR      1.582050E-08    5.505425E+06    2.598164E-09    -8.884337E-10  -4.806934E+04   1.046571E-10   9.968034E+01'
                #'                         -1.116439E-08   -6.040572E+05    1.315160E-09    -1.258955E-09  -4.381078E+05  -2.067553E-10'

                f06_file.write('      %8i %8s %-13s %-13s %-13s %-13s %-13s %-13s %s\n'
                               '      %8s %8s %-13s %-13s %-13s %-13s %-13s %s\n\n'
                               % (eid, etypei, saxr, sayr, sazr, svxr, svyr, svzr, spressurer,
                                  '', '',      saxi, sayi, sazi, svxi, svyi, svzi))
            f06_file.write(page_stamp % page_num)
            page_num += 1
        return page_num - 1

    def write_op2(self, op2, op2_ascii, itable, new_result, date,
                  is_mag_phase=False, endian='>'):
        """writes an OP2"""
        import inspect
        from struct import Struct, pack
        frame = inspect.currentframe()
        call_frame = inspect.getouterframes(frame, 2)
        op2_ascii.write('%s.write_op2: %s\n' % (self.__class__.__name__, call_frame[1][3]))

        if itable == -1:
            self._write_table_header(op2, op2_ascii, date)
            itable = -3

        eids = self.element

        # table 4 info
        #ntimes = self.data.shape[0]
        #nnodes = self.data.shape[1]
        nelements = self.data.shape[1]

        # 21 = 1 node, 3 principal, 6 components, 9 vectors, 2 p/ovm
        #ntotal = ((nnodes * 21) + 1) + (nelements * 4)

        ntotali = self.num_wide
        ntotal = ntotali * nelements

        #device_code = self.device_code
        op2_ascii.write('  ntimes = %s\n' % self.ntimes)

        eids_device = self.element * 10 + self.device_code

        if self.is_sort1:
            struct1 = Struct(endian + b'i 8s13f')
        else:
            raise NotImplementedError('SORT2')

        op2_ascii.write('nelements=%i\n' % nelements)

        etypei = self.element_type
        if etypei == 76:
            ename = b'HEXPR'
        elif etypei == 77:
            ename = b'PENPR'
        elif etypei == 78:
            ename = b'TETPR'
        else:
            raise NotImplementedError(self)
        #etypeb = self.element_type#.encode('ascii')
        for itime in range(self.ntimes):
            self._write_table_3(op2, op2_ascii, new_result, itable, itime)

            # record 4
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

            ax = self.data[itime, :, 0]
            ay = self.data[itime, :, 0]
            az = self.data[itime, :, 0]
            vx = self.data[itime, :, 0]
            vy = self.data[itime, :, 0]
            vz = self.data[itime, :, 0]
            pressure = self.data[itime, :, 0]

            for eid, eid_device, axi, ayi, azi, vxi, vyi, vzi, pressurei in zip(
                    eids, eids_device, ax, ay, az, vx, vy, vz, pressure):
                out = write_imag_floats_13e([axi, ayi, azi, vxi, vyi, vzi, pressurei], is_mag_phase)
                [saxr, sayr, sazr, svxr, svyr, svzr, spressurer,
                 saxi, sayi, sazi, svxi, svyi, svzi, spressurei] = out
                #'       1000    HEXPR      1.582050E-08    5.505425E+06    2.598164E-09    -8.884337E-10  -4.806934E+04   1.046571E-10   9.968034E+01'
                #'                         -1.116439E-08   -6.040572E+05    1.315160E-09    -1.258955E-09  -4.381078E+05  -2.067553E-10'
                data = [
                    eid_device, ename,
                    axi.real, ayi.real, azi.real, vxi.real, vyi.real, vzi.real, pressurei.real,
                    axi.imag, ayi.imag, azi.imag, vxi.imag, vyi.imag, vzi.imag,
                ]
                op2_ascii.write('      %8i %8s %-13s %-13s %-13s %-13s %-13s %-13s %s\n'
                                '      %8s %8s %-13s %-13s %-13s %-13s %-13s %s\n\n'
                                % (eid, etypei, saxr, sayr, sazr, svxr, svyr, svzr, spressurer,
                                   '', '',      saxi, sayi, sazi, svxi, svyi, svzi))
                op2.write(struct1.pack(*data))

            #for eid, eid_device, fxi, fyi, fzi, mxi, myi, mzi in zip(eids, eids_device, fx, fy, fz, mx, my, mz):
                #data = [
                    #eid_device,
                    #fxi.real, fyi.real, fzi.real, mxi.real, myi.real, mzi.real,
                    #fxi.imag, fyi.imag, fzi.imag, mxi.imag, myi.imag, mzi.imag,
                #]

                #vals = (fxi, fyi, fzi, mxi, myi, mzi)
                #vals2 = write_imag_floats_13e(vals, is_mag_phase)
                #(fxir, fyir, fzir, mxir, myir, mzir,
                 #fxii, fyii, fzii, mxii, myii, mzii) = vals2
                #op2_ascii.write('0%26i   %-13s  %-13s  %-13s  %-13s  %-13s  %s\n'
                               #' %26s   %-13s  %-13s  %-13s  %-13s  %-13s  %s\n' % (
                                   #eid, fxir, fyir, fzir, mxir, myir, mzir,
                                   #'', fxii, fyii, fzii, mxii, myii, mzii))
                #op2.write(struct1.pack(*data))

            itable -= 1
            header = [4 * ntotal,]
            op2.write(pack('i', *header))
            op2_ascii.write('footer = %s\n' % header)
            new_result = False
        return itable


class ComplexForceMomentArray(ComplexForceObject):

    def __init__(self, data_code, is_sort1, isubcase, dt):
        ComplexForceObject.__init__(self, data_code, isubcase)

        self.result_flag = 0
        #self.code = [self.format_code, self.sort_code, self.s_code]

        #self.ntimes = 0  # or frequency/mode
        #self.ntotal = 0
        self.itime = 0
        self.nelements = 0  # result specific

    def get_headers(self) -> List[str]:
        headers = ['fx', 'fy', 'fz', 'mx', 'my', 'mz']
        return headers

    @property
    def is_real(self):
        return False

    @property
    def is_complex(self):
        return True

    def _reset_indices(self):
        self.itotal = 0
        self.ielement = 0

    def build(self):
        """sizes the vectorized attributes of the ComplexCBushForceArray"""
        #print('ntimes=%s nelements=%s ntotal=%s subtitle=%s' % (
            #self.ntimes, self.nelements, self.ntotal, self.subtitle))
        if self.is_built:
            return
        nnodes = 1

        #self.names = []
        #self.nelements //= nnodes
        self.nelements /= self.ntimes
        #self.ntotal //= self.ntimes
        self.itime = 0
        self.ielement = 0
        self.itotal = 0
        self.is_built = True
        #print('ntotal=%s ntimes=%s nelements=%s' % (self.ntotal, self.ntimes, self.nelements))

        #print("ntimes=%s nelements=%s ntotal=%s" % (self.ntimes, self.nelements, self.ntotal))
        self._times = zeros(self.ntimes, 'float32')
        self.element = zeros(self.ntotal, 'int32')

        # the number is messed up because of the offset for the element's properties

        if self.nelements * nnodes != self.ntotal:
            msg = 'ntimes=%s nelements=%s nnodes=%s ne*nn=%s ntotal=%s' % (
                self.ntimes, self.nelements, nnodes, self.nelements * nnodes, self.ntotal)
            raise RuntimeError(msg)
        #[fx, fy, fz, mx, my, mz]
        self.data = zeros((self.ntimes, self.ntotal, 6), 'complex64')

    def build_dataframe(self):
        """creates a pandas dataframe"""
        #Freq                              10.0
        #ElementID Item
        #123       fx    10000.000000+0.000021j
        #          fy     1000.000000+0.000002j
        #          fz      100.000000+0.000000j
        #          mx     7000.000000+0.000000j
        #          my      700.000000+0.000000j
        #          mz       70.000000+0.000000j
        headers = self.get_headers()
        column_names, column_values = self._build_dataframe_transient_header()
        self.data_frame = self._build_pandas_transient_elements(
            column_values, column_names,
            headers, self.element, self.data)

    def __eq__(self, table):  # pragma: no cover
        self._eq_header(table)
        assert self.is_sort1 == table.is_sort1
        if not np.array_equal(self.data, table.data):
            msg = 'table_name=%r class_name=%s\n' % (self.table_name, self.__class__.__name__)
            msg += '%s\n' % str(self.code_information())
            ntimes = self.data.shape[0]

            i = 0
            if self.is_sort1:
                for itime in range(ntimes):
                    for ieid, eid in enumerate(self.element):
                        t1 = self.data[itime, ieid, :]
                        t2 = table.data[itime, ieid, :]
                        (tx1, ty1, tz1, rx1, ry1, rz1) = t1
                        (tx2, ty2, tz2, rx2, ry2, rz2) = t2
                        d = t1 - t2
                        if not allclose([tx1.real, tx1.imag, ty1.real, ty1.imag],
                                        [tx2.real, tx2.imag, ty2.real, ty2.imag], atol=0.0001):
                        #if not np.array_equal(t1, t2):
                            msg += '%-4s  (%s, %sj, %s, %sj)\n      (%s, %sj, %s, %sj)\n  dt12=(%s, %sj, %s, %sj)\n' % (
                                eid,
                                tx1.real, tx1.imag, ty1.real, ty1.imag,
                                tx2.real, tx2.imag, ty2.real, ty2.imag,
                                d[0].real, d[0].imag, d[1].real, d[1].imag,)
                            i += 1
                        if i > 10:
                            print(msg)
                            raise ValueError(msg)
            else:
                raise NotImplementedError(self.is_sort2)
            if i > 0:
                print(msg)
                raise ValueError(msg)
        return True

    def add_sort1(self, dt, eid, fx, fy, fz, mx, my, mz):
        """unvectorized method for adding SORT1 transient data"""
        assert isinstance(eid, integer_types) and eid > 0, 'dt=%s eid=%s' % (dt, eid)
        #[fx, fy, fz, mx, my, mz]
        self._times[self.itime] = dt
        self.data[self.itime, self.itotal, :] = [fx, fy, fz, mx, my, mz]
        self.element[self.itotal] = eid
        self.itotal += 1

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
        msg = []

        if self.nonlinear_factor not in (None, np.nan):  # transient
            msg.append('  type=%s ntimes=%i nelements=%i; table_name=%r\n'
                       % (self.__class__.__name__, ntimes, nelements, self.table_name))
        else:
            msg.append('  type=%s nelements=%i; table_name=%r\n' % (
                self.__class__.__name__, nelements, self.table_name))
        msg.append('  eType, cid\n')
        msg.append('  data: [ntimes, nelements, 6] where 6=[%s]\n' % str(', '.join(self.get_headers())))
        msg.append('  data.shape = %s\n' % str(self.data.shape).replace('L', ''))
        # msg.append('  is_sort1=%s is_sort2=%s\n' % (self.is_sort1, self.is_sort2))
        msg.append(f'  {self.element_name}\n')
        msg += self.get_data_code()
        return msg

    def write_f06(self, f06_file, header=None, page_stamp='PAGE %s',
                  page_num=1, is_mag_phase=False, is_sort1=True):
        if header is None:
            header = []

        if is_mag_phase:
            mag_phase = '                                                          (MAGNITUDE/PHASE)\n\n'
        else:
            mag_phase = '                                                          (REAL/IMAGINARY)\n\n'


        name = self.data_code['name']
        if name == 'freq':
            name = 'FREQUENCY'
        elif name == 'mode':
            name = 'MODE'
        else:
            raise RuntimeError(name)

        # is_sort1 = True
        if is_sort1:
            line2 = '       ID.     FORCE-X       FORCE-Y       FORCE-Z      MOMENT-X      MOMENT-Y      MOMENT-Z  \n'
        else:
            line2 = '   %26s        FORCE-X       FORCE-Y       FORCE-Z      MOMENT-X      MOMENT-Y      MOMENT-Z  \n' % name

        # force
        words = self._words()
        msg_temp = header + [
            words,
            mag_phase,
            ' ',
            # line1,
            line2,
        ]
        if self.is_sort1:
            if is_sort1:
                page_num = self._write_sort1_as_sort1(f06_file, page_num, page_stamp, header, msg_temp, is_mag_phase)
            else:
                page_num = self._write_sort1_as_sort2(f06_file, page_num, page_stamp, header, msg_temp, is_mag_phase)
        else:
            assert self.is_sort1 == True, str(self)
        return page_num - 1

    def _write_sort1_as_sort1(self, f06_file, page_num, page_stamp, header, msg_temp, is_mag_phase):
        ntimes = self.data.shape[0]
        eids = self.element
        for itime in range(ntimes):
            dt = self._times[itime]
            dt_line = ' %14s = %12.5E\n' % (self.data_code['name'], dt)
            header[1] = dt_line
            msg = header + msg_temp
            f06_file.write(''.join(msg))

            #fx, fy, fz, mx, my, mz
            if self.is_sort1:
                fx = self.data[itime, :, 0]
                fy = self.data[itime, :, 1]
                fz = self.data[itime, :, 2]
                mx = self.data[itime, :, 3]
                my = self.data[itime, :, 4]
                mz = self.data[itime, :, 5]
            else:
                fx = self.data[:, itime, 0]
                fy = self.data[:, itime, 1]
                fz = self.data[:, itime, 2]
                mx = self.data[:, itime, 3]
                my = self.data[:, itime, 4]
                mz = self.data[:, itime, 5]

            for eid, fxi, fyi, fzi, mxi, myi, mzi in zip(eids, fx, fy, fz, mx, my, mz):
                vals = (fxi, fyi, fzi, mxi, myi, mzi)
                vals2 = write_imag_floats_13e(vals, is_mag_phase)
                (fxir, fyir, fzir, mxir, myir, mzir,
                 fxii, fyii, fzii, mxii, myii, mzii) = vals2

                f06_file.write('0%26i   %-13s  %-13s  %-13s  %-13s  %-13s  %s\n'
                               ' %26s   %-13s  %-13s  %-13s  %-13s  %-13s  %s\n' % (
                                   eid, fxir, fyir, fzir, mxir, myir, mzir,
                                   '', fxii, fyii, fzii, mxii, myii, mzii))
            f06_file.write(page_stamp % page_num)
            page_num += 1
        return page_num

    def _write_sort1_as_sort2(self, f06_file, page_num, page_stamp, header, msg_temp, is_mag_phase):
        eids = self.element
        times = self._times
        for ieid, eid in enumerate(eids):
            eid_line = ' ELEMENT-ID = %s' % (eid)
            header[1] = eid_line
            msg = header + msg_temp
            f06_file.write(''.join(msg))

            if self.is_sort1:
                fx = self.data[:, ieid, 0]
                fy = self.data[:, ieid, 1]
                fz = self.data[:, ieid, 2]
                mx = self.data[:, ieid, 3]
                my = self.data[:, ieid, 4]
                mz = self.data[:, ieid, 5]
            else:
                raise RuntimeError()

            for dt, fxi, fyi, fzi, mxi, myi, mzi in zip(times, fx, fy, fz, mx, my, mz):
                vals = (fxi, fyi, fzi, mxi, myi, mzi)
                vals2 = write_imag_floats_13e(vals, is_mag_phase)
                (fxir, fyir, fzir, mxir, myir, mzir,
                 fxii, fyii, fzii, mxii, myii, mzii) = vals2

                f06_file.write('0%26s   %-13s  %-13s  %-13s  %-13s  %-13s  %s\n'
                               ' %26s   %-13s  %-13s  %-13s  %-13s  %-13s  %s\n' % (
                                   write_float_12e(dt),
                                   fxir, fyir, fzir, mxir, myir, mzir,
                                   '', fxii, fyii, fzii, mxii, myii, mzii))
            f06_file.write(page_stamp % page_num)
            page_num += 1
        return page_num

    def write_op2(self, op2, op2_ascii, itable, new_result, date,
                  is_mag_phase=False, endian='>'):
        """writes an OP2"""
        import inspect
        from struct import Struct, pack
        frame = inspect.currentframe()
        call_frame = inspect.getouterframes(frame, 2)
        op2_ascii.write('%s.write_op2: %s\n' % (self.__class__.__name__, call_frame[1][3]))

        if itable == -1:
            self._write_table_header(op2, op2_ascii, date)
            itable = -3

        eids = self.element

        # table 4 info
        #ntimes = self.data.shape[0]
        #nnodes = self.data.shape[1]
        nelements = self.data.shape[1]

        # 21 = 1 node, 3 principal, 6 components, 9 vectors, 2 p/ovm
        #ntotal = ((nnodes * 21) + 1) + (nelements * 4)

        ntotali = self.num_wide
        ntotal = ntotali * nelements

        #device_code = self.device_code
        op2_ascii.write('  ntimes = %s\n' % self.ntimes)

        eids_device = self.element * 10 + self.device_code

        if self.is_sort1:
            struct1 = Struct(endian + b'i 12f')
        else:
            raise NotImplementedError('SORT2')

        op2_ascii.write('nelements=%i\n' % nelements)

        for itime in range(self.ntimes):
            self._write_table_3(op2, op2_ascii, new_result, itable, itime)

            # record 4
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

            fx = self.data[itime, :, 0]
            fy = self.data[itime, :, 1]
            fz = self.data[itime, :, 2]
            mx = self.data[itime, :, 3]
            my = self.data[itime, :, 4]
            mz = self.data[itime, :, 5]

            for eid, eid_device, fxi, fyi, fzi, mxi, myi, mzi in zip(eids, eids_device, fx, fy, fz, mx, my, mz):
                data = [
                    eid_device,
                    fxi.real, fyi.real, fzi.real, mxi.real, myi.real, mzi.real,
                    fxi.imag, fyi.imag, fzi.imag, mxi.imag, myi.imag, mzi.imag,
                ]

                vals = (fxi, fyi, fzi, mxi, myi, mzi)
                vals2 = write_imag_floats_13e(vals, is_mag_phase)
                (fxir, fyir, fzir, mxir, myir, mzir,
                 fxii, fyii, fzii, mxii, myii, mzii) = vals2
                op2_ascii.write('0%26i   %-13s  %-13s  %-13s  %-13s  %-13s  %s\n'
                               ' %26s   %-13s  %-13s  %-13s  %-13s  %-13s  %s\n' % (
                                   eid, fxir, fyir, fzir, mxir, myir, mzir,
                                   '', fxii, fyii, fzii, mxii, myii, mzii))
                op2.write(struct1.pack(*data))

            itable -= 1
            header = [4 * ntotal,]
            op2.write(pack('i', *header))
            op2_ascii.write('footer = %s\n' % header)
            new_result = False
        return itable

    @abstractmethod
    def _words(self) -> str:
        return ''

class ComplexCBushForceArray(ComplexForceMomentArray):
    def __init__(self, data_code, is_sort1, isubcase, dt):
        ComplexForceMomentArray.__init__(self, data_code, is_sort1, isubcase, dt)
        self.element_type = 'CBUSH'

    def _words(self) -> str:
        words = '                         C O M P L E X   F O R C E S   I N   B U S H   E L E M E N T S   ( C B U S H ) \n'
        return words

class ComplexCBearForceArray(ComplexForceMomentArray):
    def __init__(self, data_code, is_sort1, isubcase, dt):
        ComplexForceMomentArray.__init__(self, data_code, is_sort1, isubcase, dt)
        self.element_type = 'CBEAR'

    #def _words(self) -> str:
        #words = '                        C O M P L E X   F O R C E S   I N   B U S H   E L E M E N T S   ( C B E A R ) \n'
        #return words


class ComplexCBeamForceVUArray(BaseElement):  # 191-VUBEAM
    """
    **ELTYPE = 191 Beam view element (VUBEAM)**

    2 PARENT I     Parent p-element identification number
    3 COORD  I     CID coordinate system identification number
    4 ICORD  CHAR4 ICORD flat/curved and so on TCODE,7 =0 Real
    5 VUGRID I     VU grid ID for output grid
    6 POSIT  RS    x/L position of VU grid identification number
    7 POS(3) RS    Y, Z, W coordinate of output point
    10 NX    RS    Normal x
    11 TXY   RS    Shear xy
    12 TZX   RS    Shear zx

    **ELTYPE = 191 Beam view element (VUBEAM)**

    TCODE,7 = 1 Real/imaginary or magnitude/phase
    5 VUGRID   I  VU grid identification number for output grid
    6 POSIT    RS x/L position of VU grid identification number

    7 FORCEXR  RS Force x real/mag.
    8 SHEARYR  RS Shear force y real/mag.
    9 SHEARZR  RS Shear force z real/mag.
    10 TORSINR RS Torsional moment x real/mag.
    11 BENDYR  RS Bending moment y real/mag.
    12 BENDZR  RS Bending moment z real/mag.

    13 FORCEXI RS Force x imag./phase
    14 SHEARYI RS Shear force y imag./phase
    15 SHEARZI RS Shear force z imag./phase
    16 TORSINI RS Torsional moment x imag./phase
    17 BENDYI  RS Bending moment y imag./phase
    18 BENDZI  RS Bending moment z imag./phase
    Words 5 through max repeat 2 times
    """
    def __init__(self, data_code, is_sort1, isubcase, dt):
        BaseElement.__init__(self, data_code, isubcase, apply_data_code=True)
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
        return False

    @property
    def is_complex(self):
        return True

    @property
    def nnodes_per_element(self):
        if self.element_type in [191]:  # VUBEAM
            nnodes_per_element = 2
        else:
            raise NotImplementedError('name=%r type=%s' % (self.element_name, self.element_type))
        return nnodes_per_element

    def _reset_indices(self):
        self.itotal = 0
        self.ielement = 0

    def get_headers(self) -> List[str]:
        return ['xxb', 'force_x', 'shear_y', 'shear_z', 'torsion', 'bending_y', 'bending_z']

    def build(self):
        """sizes the vectorized attributes of the ComplexCBendForceVUArray"""
        #print("self.ielement = %s" % self.ielement)
        #print('ntimes=%s nelements=%s ntotal=%s' % (self.ntimes, self.nelements, self.ntotal))
        assert self.ntimes > 0, 'ntimes=%s' % self.ntimes
        assert self.nelements > 0, 'nelements=%s' % self.nelements
        assert self.ntotal > 0, 'ntotal=%s' % self.ntotal

        nnodes_per_element = self.nnodes_per_element

        #print('nnodes_per_element[%s, %s] = %s' % (self.isubcase, self.element_type, nnodes_per_element))
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
            #self.element_name, self.element_type, nnodes_per_element, self.ntimes, self.nelements, self.ntotal))
        dtype, idtype, cfdtype = get_complex_times_dtype(self.nonlinear_factor, self.size)
        self._times = np.zeros(self.ntimes, dtype=dtype)
        self.element_node = np.zeros((self.ntotal, 2), dtype=idtype)
        self.parent_coord = np.zeros((self.ntotal, 2), dtype=idtype)

        #[xxb, force_x, shear_y, shear_z, torsion, bending_y, bending_z]
        self.data = np.zeros((self.ntimes, self.ntotal, 7), dtype=cfdtype)

    #def build_dataframe(self):
        #"""creates a pandas dataframe"""
        #import pandas as pd
        #headers = self.get_headers()

        #nelements = self.element_node.shape[0] // 2
        #if self.is_fiber_distance:
            #fiber_distance = ['Top', 'Bottom'] * nelements
        #else:
            #fiber_distance = ['Mean', 'Curvature'] * nelements
        #fd = np.array(fiber_distance, dtype='unicode')
        #element_node = [self.element_node[:, 0], self.element_node[:, 1], fd]

        #if self.nonlinear_factor not in (None, np.nan):
            #column_names, column_values = self._build_dataframe_transient_header()
            #self.data_frame = pd.Panel(self.data, items=column_values, major_axis=element_node, minor_axis=headers).to_frame()
            #self.data_frame.columns.names = column_names
            #self.data_frame.index.names = ['ElementID', 'NodeID', 'Location', 'Item']
        #else:
            ## option B - nice!
            #df1 = pd.DataFrame(element_node).T
            #df1.columns = ['ElementID', 'NodeID', 'Location']
            #df2 = pd.DataFrame(self.data[0])
            #df2.columns = headers
            #self.data_frame = df1.join(df2)
        #self.data_frame = self.data_frame.reset_index().replace({'NodeID': {0:'CEN'}}).set_index(['ElementID', 'NodeID', 'Location'])
        #print(self.data_frame)

    def __eq__(self, table):  # pragma: no cover
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
                    (xxb1, fx1, fy1, fz1, mx1, my1, mz1) = t1
                    (xxb2, fx2, fy2, fz2, mx2, my2, mz2) = t2

                    if not np.array_equal(t1, t2):
                        eid_nid1 = '(%s, %s)  ' % (eid, nid)
                        eid_nid2 = ' ' * len(eid_nid1)
                        msg += ('%s(%s, %s, %s, %s, %s, %s, %s)\n%s(%s, %s, %s, %s, %s, %s, %s)\n' % (
                            eid_nid1,
                            xxb1, fx1, fy1, fz1, mx1, my1, mz1,
                            eid_nid2,
                            xxb2, fx2, fy2, fz2, mx2, my2, mz2))
                        i += 1
                        if i > 10:
                            #print(msg.replace('+0j,', '0,'))
                            raise ValueError(msg.replace('0j,', '0,').replace('+0j)', ')'))
                #print(msg)
                if i > 0:
                    raise ValueError(msg.replace('0j,', '0,').replace('+0j)', ')'))
        return True

    def _add_sort1(self, dt, eid, parent, coord, icord,
                   node_id, xxb, force_x, shear_y, shear_z, torsion, bending_y, bending_z):
        assert eid is not None, eid
        assert isinstance(node_id, int), node_id
        self.element_node[self.itotal, :] = [eid, node_id]
        self.parent_coord[self.itotal, :] = [parent, coord]
        # TODO: save ICORD
        #print('parent=%r, coord=%r, icord=%r' % (parent, coord, icord))
        self.data[self.itime, self.itotal, :] = [xxb, force_x, shear_y, shear_z, torsion, bending_y, bending_z]
        self.itotal += 1

    def get_stats(self, short=False) -> List[str]:
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


    def write_f06(self, f06_file, header=None, page_stamp='PAGE %s', page_num=1, is_mag_phase=False, is_sort1=True):
        """
                 C O M P L E X   F O R C E S   I N   P - V E R S I O N   B E A M   E L E M E N T S   ( B E A M )
                                                                (REAL/IMAGINARY)
                         VU-ELEMENT ID=  100001001, P-ELEMENT ID =       1, OUTPUT COORD. ID=       0, P OF EDGES =  3

           VUGRID VUGRID DIST/    - BENDING MOMENTS -            -WEB  SHEARS -             AXIAL         TOTAL
             ID.     LENGTH      PLANE 1       PLANE 2        PLANE 1       PLANE 2          FORCE         TORQUE
        111001001     0.000    0.000000E+00 -1.598690E+05   0.000000E+00 -1.040952E+06   0.000000E+00   0.000000E+00
                               0.000000E+00  0.000000E+00   0.000000E+00  0.000000E+00   0.000000E+00   0.000000E+00
        111001002     0.333    0.000000E+00  5.328967E+04   0.000000E+00  1.872484E+05   0.000000E+00   0.000000E+00
                               0.000000E+00  0.000000E+00   0.000000E+00  0.000000E+00   0.000000E+00   0.000000E+00

                       C O M P L E X    S T R A I N S    I N   P - V E R S I O N   B E A M   E L E M E N T S   ( B E A M )
                                                                 (REAL/IMAGINARY)
                          VU-ELEMENT ID=  100001003, P-ELEMENT ID =       1, OUTPUT COORD. ID=       0, P OF EDGES =  3

            VUGRID VUGRID DIST/     LOCATION         LOCATION         LOCATION         LOCATION
              ID.     LENGTH           C                D  E                F
         111001003     0.667    -2.557904E+00    -2.557904E+00     2.557904E+00     2.557904E+00
                                 0.000000E+00     0.000000E+00     0.000000E+00     0.000000E+00
         111001004     1.000     7.673713E+00     7.673713E+00    -7.673713E+00    -7.673713E+00
                                 0.000000E+00     0.000000E+00     0.000000E+00     0.000000E+00
        """

        msg = [
            '                   C O M P L E X   F O R C E S   I N   P - V E R S I O N   B E A M   E L E M E N T S   ( B E A M )\n'
            '                                                           (REAL/IMAGINARY)\n'
            '                    VU-ELEMENT ID=  %9i, P-ELEMENT ID =%8i, OUTPUT COORD. ID=%8i, P OF EDGES =  3\n'
            '\n'
            '      VUGRID VUGRID DIST/     - BENDING MOMENTS -              -WEB  SHEARS -               AXIAL           TOTAL                    \n'
            '        ID.     LENGTH       PLANE 1       PLANE 2          PLANE 1       PLANE 2            FORCE           TORQUE   \n'
            #'   111001003     0.667     0.000000E+00  5.328967E+04     0.000000E+00 -1.872484E+05     0.000000E+00     0.000000E+00'
            #'                           0.000000E+00  0.000000E+00     0.000000E+00  0.000000E+00     0.000000E+00     0.000000E+00'
            #'   111001004     1.000     0.000000E+00 -1.598690E+05     0.000000E+00  1.040952E+06     0.000000E+00     0.000000E+00'
            #'                           0.000000E+00  0.000000E+00     0.000000E+00  0.000000E+00     0.000000E+00     0.000000E+00'

            #'                 C O M P L E X    S T R A I N S    I N   P - V E R S I O N   B E A M   E L E M E N T S   ( B E A M )\n'
            #'                                                           (REAL/IMAGINARY)\n'
            #'                    VU-ELEMENT ID=  %9i, P-ELEMENT ID =       1, OUTPUT COORD. ID=       0, P OF EDGES =  3\n'
            #'\n'
            #'      VUGRID VUGRID DIST/     LOCATION         LOCATION         LOCATION         LOCATION                                            \n'
            #'        ID.     LENGTH           C                D  E                F                                                              \n'
            #'   111001003     0.667    -2.557904E+00    -2.557904E+00     2.557904E+00     2.557904E+00'
            #'                           0.000000E+00     0.000000E+00     0.000000E+00     0.000000E+00'
            #'   111001004     1.000     7.673713E+00     7.673713E+00    -7.673713E+00    -7.673713E+00'
            #'                           0.000000E+00     0.000000E+00     0.000000E+00     0.000000E+00'
        ]
        if header is None:
            header = []
        #msg, nnodes, cen = _get_plate_msg(self)

        # write the f06
        ntimes = self.data.shape[0]

        eids = self.element_node[:, 0]
        nids = self.element_node[:, 1]
        parent = self.parent_coord[:, 0]
        coord = self.parent_coord[:, 1]

        for itime in range(ntimes):
            dt = self._times[itime]
            header = _eigenvalue_header(self, header, itime, ntimes, dt)

            #[xxb, force_x, shear_y, shear_z, torsion, bending_y, bending_z]
            xxb = self.data[itime, :, 0]
            fx = self.data[itime, :, 1]
            fy = self.data[itime, :, 2]
            fz = self.data[itime, :, 3]
            mx = self.data[itime, :, 4]
            my = self.data[itime, :, 5]
            mz = self.data[itime, :, 6]

            for (i, eid, parenti, coordi, nid, xxbi, fxi, fyi, fzi, mxi, myi, mzi) in zip(
                 cycle(range(2)), eids, parent, coord, nids, xxb, fx, fy, fz, mx, my, mz):
                if i == 0:
                    f06_file.write(''.join(header + msg) % (eid, parenti, coordi))

                #out = write_imag_floats_13e([fxi, fyi, fzi, mxi, myi, mzi], is_mag_phase=is_mag_phase)
                #[fxri, fyri, fzri, mxri, myri, mzri,
                 #fxii, fyii, fzii, mxii, myii, mzii] = out

                        #   nid   xxb
                f06_file.write(
                    '   %9i     %.3f    %13.6E %13.6E    %13.6E %13.6E    %13.6E    %13.6E\n'
                    '                          %13.6E %13.6E    %13.6E %13.6E    %13.6E    %13.6E\n' % (
                        nid, xxbi.real,
                        myi.real, mzi.real, fyi.real, fzi.real, fxi.real, mxi.real,
                        myi.imag, mzi.imag, fyi.imag, fzi.imag, fxi.imag, mxi.imag,
                ))

                # stress/strain
                #f06_file.write(
                    #'   %9i     %.3s      %13.6E  %13.6E  %13.6E  %13.6E  %13.6E  %13.6E\n'
                    #'                          %13.6E  %13.6E  %13.6E  %13.6E  %13.6E  %13.6E\n' % (
                        #nid, xxbi.real,
                        #fxi.real, fyi.real, fzi.real, mxi.real, myi.real, mzi.real,
                        #fxi.imag, fyi.imag, fzi.imag, mxi.imag, myi.imag, mzi.imag,
                #))

                if i == 1:
                    f06_file.write(page_stamp % page_num + '\n')
                    page_num += 1
        return page_num - 1


class ComplexForceVU_2DArray(BaseElement):  # 189-VUQUAD,190-VUTRIA
    def __init__(self, data_code, is_sort1, isubcase, dt):
        BaseElement.__init__(self, data_code, isubcase)

        #self.parent = {}
        #self.coord = {}
        #self.icord = {}
        #self.theta = {}

        #self.ntimes = 0  # or frequency/mode
        #self.ntotal = 0
        self.nelements = 0  # result specific
        self.ntimes = 0

        # TODO if dt=None, handle SORT1 case
        self.dt = dt
        #if is_sort1:
            #if dt is not None:
                #self.add = self.add_sort1
        #else:
            #assert dt is not None
            #self.add = self.add_sort2

    @property
    def is_real(self) -> bool:
        return False
    @property
    def is_complex(self) -> bool:
        return True

    def _reset_indices(self):
        self.itotal = 0
        self.ielement = 0

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

        msg = []
        if self.nonlinear_factor not in (None, np.nan):  # transient
            msg.append('  type=%s ntimes=%i nelements=%i; table_name=%r\n'
                       % (self.__class__.__name__, ntimes, nelements, self.table_name))
            ntimes_word = 'ntimes'
        else:
            msg.append('  type=%s nelements=%i; table_name=%r\n'
                       % (self.__class__.__name__, nelements, self.table_name))
            ntimes_word = '1'
        msg.append('  eType\n')
        headers = self.get_headers()
        n = len(headers)
        msg.append('  data: [%s, nnodes, %i] where %i=[%s]\n' % (
            ntimes_word, n, n, str(', '.join(headers))))
        msg.append('  data.shape = %s\n' % str(self.data.shape).replace('L', ''))
        #msg.append('  element type: %s\n' % self.element_type)
        msg.append('  element name: %s\n' % self.element_name)
        msg += self.get_data_code()
        return msg

    def build(self):
        """sizes the vectorized attributes of the ComplexCShearForceArray"""
        #print('%s ntimes=%s nelements=%s ntotal=%s' % (
            #self.element_type, self.ntimes, self.nelements, self.ntotal))
        assert self.ntimes > 0, 'ntimes=%s' % self.ntimes
        assert self.nelements > 0, 'nelements=%s' % self.nelements
        assert self.ntotal > 0, 'ntotal=%s' % self.ntotal
        #self.names = []
        self.nelements //= self.ntimes
        self.ntotal = self.nelements
        self.itime = 0
        self.ielement = 0
        self.itotal = 0
        #self.ntimes = 0
        #self.nelements = 0
        self.is_built = True

        #print("ntimes=%s nelements=%s ntotal=%s" % (self.ntimes, self.nelements, self.ntotal))
        dtype, idtype, cfdtype = get_complex_times_dtype(self.nonlinear_factor, self.size)
        self._times = zeros(self.ntimes, dtype=dtype)
        self.element_node = zeros((self.nelements, 2), dtype=idtype)

        #[membrane_x, membrane_y, membrane_xy, bending_x, bending_y, bending_xy,
        # shear_yz, shear_xz]
        self.data = zeros((self.ntimes, self.nelements, 8), dtype=cfdtype)

    def get_headers(self) -> List[str]:
        headers = [
            'membrane_x', 'membrane_y', 'membrane_xy',
            'bending_x', 'bending_y', 'bending_xy',
            'shear_yz', 'shear_xz']
        return headers

    def add_sort1(self, nnodes, dt, eid, parent, coord, icord, theta, vugrids, forces):
        """unvectorized method for adding SORT1 transient data"""
        assert isinstance(eid, integer_types) and eid > 0, 'dt=%s eid=%s' % (dt, eid)
        self._times[self.itime] = dt

        #self.parent[eid] = parent
        #self.coord[eid] = coord
        #self.icord[eid] = icord
        #self.theta[eid] = theta
        for vugrid, force in zip(vugrids, forces):
            self.element_node[self.ielement, :] = [eid, vugrid]
            self.data[self.itime, self.ielement, :] = force
            self.ielement += 1
#'                  C O M P L E X   F O R C E S   I N   Q U A D R I L A T E R A L   E L E M E N T S   ( Q U A D 8 )'
#'                                                          (REAL/IMAGINARY)'
#' '
#'    ELEMENT                    - MEMBRANE  FORCES -                      - BENDING   MOMENTS -            - TRANSVERSE SHEAR FORCES -'
#'      ID       GRID-ID     FX            FY            FXY           MX            MY            MXY           QX            QY'
#'0       100    CEN/8  0.0           0.0           0.0           0.0           0.0           0.0          -3.492460E-10 -1.368206E-09'
#'                      0.0           0.0           0.0           0.0           0.0           0.0           2.910383E-11  5.088840E-10'
#''
