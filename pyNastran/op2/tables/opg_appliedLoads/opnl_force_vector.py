import numpy as np
from pyNastran.op2.result_objects.table_object import RealTableArray, ComplexTableArray


class RealForceVectorArray(RealTableArray):  # table_code=2, sort_code=0, thermal=0

    def __init__(self, data_code, is_sort1, isubcase, dt):
        RealTableArray.__init__(self, data_code, is_sort1, isubcase, dt)

    def write_f06(self, f06_file, header=None, page_stamp='PAGE %s', page_num=1, is_mag_phase=False, is_sort1=True):
        if header is None:
            header = []
        words = ['                                         N O N - L I N E A R - F O R C E   V E C T O R\n']
        #words += self.get_table_marker()
        if self.nonlinear_factor not in (None, np.nan):
            return self._write_f06_transient_block(words, header, page_stamp, page_num, f06_file, write_words=True,
                                                   is_mag_phase=is_mag_phase, is_sort1=is_sort1)
        return self._write_f06_block(words, header, page_stamp, page_num, f06_file, write_words=True)


class ComplexForceVectorArray(ComplexTableArray):
    def __init__(self, data_code, is_sort1, isubcase, dt):
        ComplexTableArray.__init__(self, data_code, is_sort1, isubcase, dt)

    def write_f06(self, f06_file, header=None, page_stamp='PAGE %s',
                  page_num=1, is_mag_phase=False, is_sort1=True):
        if header is None:
            header = []
        words = ['                                       C O M P L E X   F O R C E   V E C T O R\n',]
        return self._write_f06_transient_block(
            words, header, page_stamp, page_num, f06_file, is_mag_phase=is_mag_phase, is_sort1=is_sort1)
