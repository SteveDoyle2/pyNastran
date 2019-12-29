from typing import List
import numpy as np
from numpy import zeros, searchsorted, allclose

from pyNastran.utils.numpy_utils import integer_types
from pyNastran.op2.tables.oes_stressStrain.real.oes_objects import (
    StressObject, StrainObject, OES_Object, SORT2_TABLE_NAME_MAP)
from pyNastran.f06.f06_formatting import write_floats_13e, _eigenvalue_header #, get_key0


class RandomRodArray(OES_Object):
    def __init__(self, data_code, is_sort1, isubcase, dt):
        OES_Object.__init__(self, data_code, isubcase, apply_data_code=False)
        #self.code = [self.format_code, self.sort_code, self.s_code]

        self.nelements = 0  # result specific

    @property
    def is_real(self) -> bool:
        return True

    @property
    def is_complex(self) -> bool:
        return False

    @property
    def nnodes_per_element(self) -> int:
        return 1

    def _reset_indices(self):
        self.itotal = 0
        self.ielement = 0

    def _get_msgs(self):
        raise NotImplementedError()

    def get_headers(self):
        raise NotImplementedError()

    def build(self):
        """sizes the vectorized attributes of the RealRodArray"""
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
        self.build_data(self.ntimes, self.nelements, dtype)

    def build_data(self, ntimes, nelements, dtype):
        """actually performs the build step"""
        self.ntimes = ntimes
        self.nelements = nelements
        self._times = zeros(ntimes, dtype=dtype)
        self.element = zeros(nelements, dtype='int32')

        #[axial, torsion]
        self.data = zeros((ntimes, nelements, 2), dtype='float32')

    def build_dataframe(self):
        """creates a pandas dataframe"""
        import pandas as pd
        headers = self.get_headers()
        if self.nonlinear_factor not in (None, np.nan):
            #ElementID                  1101          1102
            #ElementID Item
            #1102      axial    1.079339e+02  4.400043e-03
            #          torsion  2.439666e-04  0.000000e+00
            #          axial    3.145797e+01  1.282420e-03
            #          torsion  7.110654e-05  0.000000e+00
            #          axial    1.572898e+01  6.412102e-04
            #          torsion  3.555327e-05  0.000000e+00
            #          axial    3.495329e+00  1.424911e-04
            #          torsion  7.900725e-06  0.000000e+00
            #          axial    2.928574e-07  1.193838e-11
            #          torsion  6.618416e-13  0.000000e+00
            column_names, column_values = self._build_dataframe_transient_header()
            #if is_v25:
                #print(f'skipping pandas {self.class_name}')
                #return
            data_frame = self._build_pandas_transient_elements(
                column_values, column_names,
                headers, self.element, self.data)
            #column_names, column_values = self._build_dataframe_transient_header()
            #self.data_frame = pd.Panel(self.data, items=column_values,
                                       #major_axis=self.element, minor_axis=headers).to_frame()
            #self.data_frame.columns.names = column_names
            #self.data_frame.index.names = ['ElementID', 'Item']
        else:
            data_frame = pd.Panel(self.data, major_axis=self.element, minor_axis=headers).to_frame()
            data_frame.columns.names = ['Static']
            data_frame.index.names = ['ElementID', 'Item']
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
                        (axial1, torsion1) = t1
                        (axial2, torsion2) = t2
                        if not allclose(t1, t2):
                        #if not np.array_equal(t1, t2):
                            msg += '%s\n  (%s, %s)\n  (%s, %s)\n' % (
                                eid,
                                axial1, torsion1,
                                axial2, torsion2)
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

    def add_sort1(self, dt, eid, axial, torsion):
        self._times[self.itime] = dt
        #if self.itime == 0:
        #print('itime=%s eid=%s' % (self.itime, eid))
        self.element[self.ielement] = eid
        self.data[self.itime, self.ielement, :] = [axial, torsion]
        self.ielement += 1

    def get_stats(self, short=False) -> List[str]:
        if not self.is_built:
            return [
                '<%s>\n' % self.__class__.__name__,
                '  ntimes: %i\n' % self.ntimes,
                '  ntotal: %i\n' % self.ntotal,
            ]

        ntimes, nelements, _ = self.data.shape
        assert self.ntimes == ntimes, 'ntimes=%s expected=%s' % (self.ntimes, ntimes)
        assert self.nelements == nelements, 'nelements=%s expected=%s' % (self.nelements, nelements)

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
        msg.append('  element.shape = %s\n' % str(self.element.shape).replace('L', ''))
        msg.append('  data.shape = %s\n' % str(self.data.shape).replace('L', ''))
        msg.append('  element type: %s\n' % self.element_name)
        msg += self.get_data_code()
        return msg

    def get_f06_header(self, is_mag_phase=True):
        crod_msg, conrod_msg, ctube_msg = self._get_msgs()
        if 'CROD' in self.element_name:
            msg = crod_msg
        elif 'CONROD' in self.element_name:
            msg = conrod_msg
        elif 'CTUBE' in self.element_name:
            msg = ctube_msg
        else:
            raise NotImplementedError(self.element_name)
        return self.element_name, msg

    def get_element_index(self, eids):
        # elements are always sorted; nodes are not
        itot = searchsorted(eids, self.element)  #[0]
        return itot

    def eid_to_element_node_index(self, eids):
        #ind = ravel([searchsorted(self.element_node[:, 0] == eid) for eid in eids])
        ind = searchsorted(eids, self.element)
        #ind = ind.reshape(ind.size)
        #ind.sort()
        return ind

    def write_f06(self, f06_file, header=None, page_stamp='PAGE %s', page_num=1, is_mag_phase=False, is_sort1=True):
        if header is None:
            header = []
        (elem_name, msg_temp) = self.get_f06_header(is_mag_phase)
        if self.is_sort1:
            page_num = self._write_sort1_as_sort1(header, page_stamp, page_num, f06_file, msg_temp)
        return page_num

    def _write_sort1_as_sort1(self, header, page_stamp, page_num, f06_file, msg_temp):
        print('update the RandomRodArray header')
        ntimes = self.data.shape[0]

        eids = self.element
        is_odd = False
        nwrite = len(eids)
        if len(eids) % 2 == 1:
            nwrite -= 1
            is_odd = True

        for itime in range(ntimes):
            dt = self._times[itime]
            header = _eigenvalue_header(self, header, itime, ntimes, dt)
            f06_file.write(''.join(header + msg_temp))

            #print("self.data.shape=%s itime=%s ieids=%s" % (str(self.data.shape), itime, str(ieids)))
            axial = self.data[itime, :, 0]
            torsion = self.data[itime, :, 1]

            out = []
            for eid, axiali, torsioni in zip(eids, axial, torsion):
                [axiali, torsioni] = write_floats_13e([axiali, torsioni])
                out.append([eid, axiali, torsioni])

            for i in range(0, nwrite, 2):
                f06_file.write(
                    '      %8i %-13s  %-13s %-8i   %-13s  %-s\n' % (
                        tuple(out[i] + out[i + 1])))
            if is_odd:
                f06_file.write('      %8i %-13s  %13s\n' % (tuple(out[-1])))
            f06_file.write(page_stamp % page_num)
            page_num += 1
        return page_num - 1


class RandomBushStressArray(RandomRodArray, StressObject):
    def __init__(self, data_code, is_sort1, isubcase, dt):
        RandomRodArray.__init__(self, data_code, is_sort1, isubcase, dt)
        StressObject.__init__(self, data_code, isubcase)

    def get_headers(self) -> List[str]:
        headers = ['axial', 'SMa', 'torsion', 'SMt']
        return headers

    def _get_msgs(self):
        raise NotImplementedError()
        #base_msg = ['       ELEMENT       AXIAL       SAFETY      TORSIONAL     SAFETY       ELEMENT       AXIAL       SAFETY      TORSIONAL     SAFETY\n',
                    #'         ID.        STRESS       MARGIN        STRESS      MARGIN         ID.        STRESS       MARGIN        STRESS      MARGIN\n']
        #crod_msg   = ['                                     S T R E S S E S   I N   R O D   E L E M E N T S      ( C R O D )\n', ]
        #conrod_msg = ['                                     S T R E S S E S   I N   R O D   E L E M E N T S      ( C O N R O D )\n', ]
        #ctube_msg  = ['                                     S T R E S S E S   I N   R O D   E L E M E N T S      ( C T U B E )\n', ]
        ##cbush_msg  = ['                                     S T R E S S E S   I N   R O D   E L E M E N T S      ( C B U S H )\n', ]
        #crod_msg += base_msg
        #conrod_msg += base_msg
        #ctube_msg += base_msg
        ##cbush_msg += base_msg
        #return crod_msg, conrod_msg, ctube_msg


class RandomRodStressArray(RandomRodArray, StressObject):
    def __init__(self, data_code, is_sort1, isubcase, dt):
        RandomRodArray.__init__(self, data_code, is_sort1, isubcase, dt)
        StressObject.__init__(self, data_code, isubcase)

    def get_headers(self) -> List[str]:
        headers = ['axial', 'torsion']
        return headers

    def _get_msgs(self):
        base_msg = ['       ELEMENT       AXIAL       TORSIONAL     SAFETY       ELEMENT       AXIAL       TORSIONAL     SAFETY\n',
                    '         ID.        STRESS         STRESS      MARGIN         ID.        STRESS         STRESS      MARGIN\n']
        crod_msg = ['                                     S T R E S S E S   I N   R O D   E L E M E N T S      ( C R O D )\n', ]
        conrod_msg = ['                                     S T R E S S E S   I N   R O D   E L E M E N T S      ( C O N R O D )\n', ]
        ctube_msg = ['                                     S T R E S S E S   I N   R O D   E L E M E N T S      ( C T U B E )\n', ]
        crod_msg += base_msg
        conrod_msg += base_msg
        ctube_msg += base_msg
        return crod_msg, conrod_msg, ctube_msg

class RandomRodStrainArray(RandomRodArray, StrainObject):
    def __init__(self, data_code, is_sort1, isubcase, dt):
        RandomRodArray.__init__(self, data_code, is_sort1, isubcase, dt)
        StrainObject.__init__(self, data_code, isubcase)

    def get_headers(self) -> List[str]:
        headers = ['axial', 'torsion']
        return headers

    def _get_msgs(self):
        # TODO: update this...
        #
        #                                     S T R E S S E S   I N   R O D   E L E M E N T S      ( C T U B E )
        #                                                     ( ROOT MEAN SQUARE )
        #
        #     ELEMENT         AXIAL                   TORSIONAL                ELEMENT         AXIAL                   TORSIONAL
        #       ID.          STRESS                     STRESS                   ID.          STRESS                     STRESS
        #        3306      1.038912E+01              7.809715E-02                 3307      6.633282E-02              0.0
        base_msg = ['       ELEMENT       AXIAL       TORSIONAL     SAFETY       ELEMENT       AXIAL       TORSIONAL     SAFETY\n',
                    '         ID.        STRAIN         STRAIN      MARGIN         ID.        STRAIN         STRAIN      MARGIN\n']
        crod_msg = ['                                       S T R A I N S   I N   R O D   E L E M E N T S      ( C R O D )\n', ]
        conrod_msg = ['                                       S T R A I N S   I N   R O D   E L E M E N T S      ( C O N R O D )\n', ]
        ctube_msg = ['                                       S T R A I N S   I N   R O D   E L E M E N T S      ( C T U B E )\n', ]
        crod_msg += base_msg
        conrod_msg += base_msg
        ctube_msg += base_msg
        return crod_msg, conrod_msg, ctube_msg
