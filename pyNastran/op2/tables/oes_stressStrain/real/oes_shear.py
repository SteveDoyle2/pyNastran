from typing import List

import numpy as np
from numpy import zeros, allclose

from pyNastran.utils.numpy_utils import integer_types
from pyNastran.op2.tables.oes_stressStrain.real.oes_objects import StressObject, StrainObject, OES_Object
from pyNastran.f06.f06_formatting import _eigenvalue_header #, get_key0


class RealShearArray(OES_Object):
    def __init__(self, data_code, is_sort1, isubcase, dt):
        OES_Object.__init__(self, data_code, isubcase, apply_data_code=False)
        #self.code = [self.format_code, self.sort_code, self.s_code]
        self.nelements = 0  # result specific

        #if is_sort1:
            #self.add_new_eid = self.add_new_eid_sort1
        #else:
            #raise NotImplementedError('SORT2')

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
        """sizes the vectorized attributes of the RealShearArray"""
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
        _times = zeros(self.ntimes, dtype=dtype)
        element = zeros(self.nelements, dtype='int32')

        # [max_shear, avg_shear, margin]
        data = zeros((self.ntimes, self.ntotal, 3), dtype='float32')

        if self.load_as_h5:
            #for key, value in sorted(self.data_code.items()):
                #print(key, value)
            group = self._get_result_group()
            self._times = group.create_dataset('_times', data=_times)
            self.element = group.create_dataset('element', data=element)
            self.data = group.create_dataset('data', data=data)
        else:
            self._times = _times
            self.element = element
            self.data = data

    def build_dataframe(self):
        """creates a pandas dataframe"""
        import pandas as pd

        headers = self.get_headers()
        if self.nonlinear_factor not in (None, np.nan):
            #Mode                            1             2             3
            #Freq                 1.482246e-10  3.353940e-09  1.482246e-10
            #Eigenvalue          -8.673617e-19  4.440892e-16  8.673617e-19
            #Radians              9.313226e-10  2.107342e-08  9.313226e-10
            #ElementID Item
            #22        max_shear  8.050749e-13  5.871460e-07  2.035239e-12
            #         avg_shear -8.050749e-13  5.871460e-07  2.035239e-12
            #         margin     1.401298e-45  1.401298e-45  1.401298e-45
            column_names, column_values = self._build_dataframe_transient_header()
            data_frame = self._build_pandas_transient_elements(
                column_values, column_names,
                headers, self.element, self.data)
        else:
            #Static     axial           SMa  torsion           SMt
            #ElementID
            #14           0.0  1.401298e-45      0.0  1.401298e-45
            #15           0.0  1.401298e-45      0.0  1.401298e-45
            data_frame = pd.DataFrame(self.data[0], columns=headers, index=self.element)
            data_frame.index.name = 'ElementID'
            data_frame.columns.names = ['Static']
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
                    for ieid, eid, in enumerate(self.element):
                        t1 = self.data[itime, ieid, :]
                        t2 = table.data[itime, ieid, :]
                        (max_shear, avg_shear, margin) = t1
                        (max_shear2, avg_shear2, margin2) = t2
                        if not allclose(t1, t2):
                        #if not np.array_equal(t1, t2):
                            msg += '%s\n  (%s, %s, %s)\n  (%s, %s, %s)\n' % (
                                eid,
                                max_shear, avg_shear, margin,
                                max_shear2, avg_shear2, margin2)
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

    def add_sort1(self, dt, eid, max_shear, avg_shear, margin):
        """
        ELEMENT            MAX            AVG        SAFETY         ELEMENT            MAX            AVG        SAFETY
          ID.             SHEAR          SHEAR       MARGIN           ID.             SHEAR          SHEAR       MARGIN
            328        1.721350E+03   1.570314E+03   7.2E+01
        """
        self._times[self.itime] = dt
        self.element[self.ielement] = eid
        self.data[self.itime, self.ielement, :] = [max_shear, avg_shear, margin]
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
            msg.append('  type=%s ntimes=%i nelements=%i\n'
                       % (self.__class__.__name__, ntimes, nelements))
            ntimes_word = 'ntimes'
        else:
            msg.append('  type=%s nelements=%i\n'
                       % (self.__class__.__name__, nelements))
            ntimes_word = '1'
        msg.append('  eType\n')
        headers = self.get_headers()
        n = len(headers)
        msg.append('  data: [%s, nnodes, %i] where %i=[%s]\n' % (ntimes_word, n, n, str(', '.join(headers))))
        msg.append('  element.shape = %s\n' % str(self.element.shape).replace('L', ''))
        msg.append('  data.shape = %s\n' % str(self.data.shape).replace('L', ''))
        msg.append('  element type: %s\n' % self.element_name)
        msg += self.get_data_code()
        return msg

    def get_f06_header(self):
        raise NotImplementedError('CSHEAR...')

    def write_f06(self, f06_file, header=None, page_stamp='PAGE %s',
                  page_num=1, is_mag_phase=False, is_sort1=True):
        if header is None:
            header = []
        msg_temp = self.get_f06_header()

        # write the f06
        ntimes = self.data.shape[0]

        eids = self.element
        is_odd = False
        nwrite = len(eids)
        if len(eids) % 2 == 1:
            nwrite -= 1
            is_odd = True

        for itime in range(ntimes):
            dt = self._times[itime]  # TODO: rename this...
            header = _eigenvalue_header(self, header, itime, ntimes, dt)
            f06_file.write(''.join(header + msg_temp))

            #print("self.data.shape=%s itime=%s ieids=%s" % (str(self.data.shape), itime, str(ieids)))
            max_shear = self.data[itime, :, 0]
            avg_shear = self.data[itime, :, 1]
            margin = self.data[itime, :, 2]

            out = []
            for eid, max_sheari, avg_sheari, margini in zip(eids, max_shear, avg_shear, margin):
                #[max_sheari, avg_sheari, margini] = write_floats_13e([max_sheari, avg_sheari, margini])
                out.append([eid, max_sheari, avg_sheari, margini])

            for i in range(0, nwrite, 2):
                f06_file.write('      %8i   %13s  %10.4E %13s  %8i   %13s  %10.4E %s\n' % (
                    tuple(out[i] + out[i + 1])))
            if is_odd:
                f06_file.write('      %8i   %13s  %10.4E %s\n' % tuple(out[-1]))
            f06_file.write(page_stamp % page_num)
            page_num += 1
        return page_num - 1


class RealShearStressArray(RealShearArray, StressObject):
    def __init__(self, data_code, is_sort1, isubcase, dt):
        RealShearArray.__init__(self, data_code, is_sort1, isubcase, dt)
        StressObject.__init__(self, data_code, isubcase)

    def get_headers(self) -> List[str]:
        headers = ['max_shear', 'avg_shear', 'margin']
        return headers

    def get_f06_header(self):
        msg = [
            '                                     S T R E S S E S   I N   S H E A R   P A N E L S      ( C S H E A R )\n'
            '      ELEMENT            MAX            AVG        SAFETY         ELEMENT            MAX            AVG        SAFETY\n'
            '        ID.             SHEAR          SHEAR       MARGIN           ID.             SHEAR          SHEAR       MARGIN\n'
            #'          328        1.721350E+03   1.570314E+03   7.2E+01'
        ]
        return msg


class RealShearStrainArray(RealShearArray, StrainObject):
    def __init__(self, data_code, is_sort1, isubcase, dt):
        RealShearArray.__init__(self, data_code, is_sort1, isubcase, dt)
        StrainObject.__init__(self, data_code, isubcase)

    def get_headers(self) -> List[str]:
        headers = ['max_shear', 'avg_shear', 'margin']
        return headers

    def get_f06_header(self):
        msg = [
            '                                     S T R A I N S   I N   S H E A R   P A N E L S      ( C S H E A R )\n'
            '      ELEMENT            MAX            AVG        SAFETY         ELEMENT            MAX            AVG        SAFETY\n'
            '        ID.             SHEAR          SHEAR       MARGIN           ID.             SHEAR          SHEAR       MARGIN\n'
            #'          328        1.721350E+03   1.570314E+03   7.2E+01'
        ]
        return msg
