import numpy as np
from numpy import zeros, allclose

from pyNastran.op2.tables.oes_stressStrain.real.oes_objects import StressObject, StrainObject, OES_Object
from pyNastran.op2.result_objects.op2_objects import get_times_dtype
from pyNastran.f06.f06_formatting import _eigenvalue_header #, get_key0


class RandomShearArray(OES_Object):
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

    def _reset_indices(self) -> None:
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

        #print("ntimes=%s nelements=%s ntotal=%s" % (self.ntimes, self.nelements, self.ntotal))
        if self.is_sort1:
            ntimes = self.ntimes
            nelements = self.nelements
        else:
            ntimes = self.nelements
            nelements = self.ntimes

        dtype, idtype, fdtype = get_times_dtype(self.nonlinear_factor, self.size, self.analysis_fmt)
        self._times = zeros(ntimes, dtype=dtype)
        self.element = zeros(nelements, dtype='int32')

        # [max_shear, avg_shear]
        self.data = zeros((ntimes, nelements, 2), dtype='float32')

    def build_dataframe(self):
        """creates a pandas dataframe"""
        import pandas as pd
        headers = self.get_headers()
        if self.nonlinear_factor not in (None, np.nan):
            column_names, column_values = self._build_dataframe_transient_header()
            self.data_frame = pd.Panel(self.data, items=column_values,
                                       major_axis=self.element, minor_axis=headers).to_frame()
            self.data_frame.columns.names = column_names
            self.data_frame.index.names = ['ElementID', 'Item']
        else:
            self.data_frame = pd.Panel(self.data, major_axis=self.element, minor_axis=headers).to_frame()
            self.data_frame.columns.names = ['Static']
            self.data_frame.index.names = ['ElementID', 'Item']

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
                        (max_shear1, avg_shear1) = t1
                        (max_shear2, avg_shear2) = t2
                        if not allclose(t1, t2):
                        #if not np.array_equal(t1, t2):
                            msg += '%s\n  (%s, %s)\n  (%s, %s)\n' % (
                                eid,
                                max_shear1, avg_shear1,
                                max_shear2, avg_shear2)
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

    def add_sort1(self, dt, eid, max_shear, avg_shear):
        """
        ELEMENT            MAX            AVG       ELEMENT            MAX            AVG
          ID.             SHEAR          SHEAR        ID.             SHEAR          SHEAR
            328        1.721350E+03   1.570314E+03
        """
        assert self.sort_method == 1, self
        self._times[self.itime] = dt
        self.element[self.ielement] = eid
        self.data[self.itime, self.ielement, :] = [max_shear, avg_shear]
        self.ielement += 1

    def add_sort2(self, dt, eid, max_shear, avg_shear):
        """
        ELEMENT            MAX            AVG       ELEMENT            MAX            AVG
          ID.             SHEAR          SHEAR        ID.             SHEAR          SHEAR
            328        1.721350E+03   1.570314E+03
        """
        assert self.is_sort2, self
        itime = self.itotal
        ielement = self.itime
        #print(f'dt={dt:g} eid={eid} itime={itime} ielement={ielement}')
        self._times[itime] = dt
        self.element[ielement] = eid
        self.data[itime, ielement, :] = [max_shear, avg_shear]
        self.itotal += 1

    def get_stats(self, short: bool=False) -> list[str]:
        if not self.is_built:
            return [
                f'<{self.__class__.__name__}>; table_name={self.table_name!r}\n',
                f'  ntimes: {self.ntimes:d}\n',
                f'  ntotal: {self.ntotal:d}\n',
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
        msg.append(f'  element.shape = {self.element.shape}\n')
        msg.append(f'  data.shape = {self.data.shape}\n')
        msg.append(f'  element type: {self.element_name}-{self.element_type}\n')
        msg += self.get_data_code()
        return msg

    def get_f06_header(self):
        raise NotImplementedError('CSHEAR...')

    def write_f06(self, f06_file, header=None, page_stamp='PAGE %s',
                  page_num: int=1, is_mag_phase: bool=False, is_sort1: bool=True):
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

            out = []
            for eid, max_sheari, avg_sheari in zip(eids, max_shear, avg_shear):
                #[max_sheari, avg_sheari] = write_floats_13e([max_sheari, avg_sheari])
                out.append([eid, max_sheari, avg_sheari])

            for i in range(0, nwrite, 2):
                f06_file.write('      %8i   %13s  %13s  %8i   %13s  %s\n' % (
                    tuple(out[i] + out[i + 1])))
            if is_odd:
                f06_file.write('      %8i   %13s  %s\n' % tuple(out[-1]))
            f06_file.write(page_stamp % page_num)
            page_num += 1
        return page_num - 1


class RandomShearStressArray(RandomShearArray, StressObject):
    def __init__(self, data_code, is_sort1, isubcase, dt):
        RandomShearArray.__init__(self, data_code, is_sort1, isubcase, dt)
        StressObject.__init__(self, data_code, isubcase)

    def get_headers(self) -> list[str]:
        headers = ['max_shear', 'avg_shear']
        return headers

    def get_f06_header(self):
        msg = [
            '                                     S T R E S S E S   I N   S H E A R   P A N E L S      ( C S H E A R )\n'
            '      ELEMENT            MAX            AVG        ELEMENT            MAX            AVG \n'
            '        ID.             SHEAR          SHEAR         ID.             SHEAR          SHEAR\n'
            #'          328        1.721350E+03   1.570314E+03
        ]
        return msg


class RandomShearStrainArray(RandomShearArray, StrainObject):
    def __init__(self, data_code, is_sort1, isubcase, dt):
        RandomShearArray.__init__(self, data_code, is_sort1, isubcase, dt)
        StrainObject.__init__(self, data_code, isubcase)

    def get_headers(self) -> list[str]:
        headers = ['max_shear', 'avg_shear']
        return headers

    def get_f06_header(self):
        assert self.table_name in [], self.table_name
        msg = [
            '                                     S T R A I N S   I N   S H E A R   P A N E L S      ( C S H E A R )\n'
            '      ELEMENT            MAX            AVG        ELEMENT            MAX            AVG \n'
            '        ID.             SHEAR          SHEAR         ID.             SHEAR          SHEAR\n'
            #'          328        1.721350E+03   1.570314E+03   7.2E+01'
        ]
        return msg
