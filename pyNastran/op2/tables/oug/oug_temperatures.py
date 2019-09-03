# pylint: disable=E1101
import numpy as np
from pyNastran.op2.result_objects.scalar_table_object import RealScalarTableArray

class RealTemperatureArray(RealScalarTableArray):
    def __init__(self, data_code, is_sort1, isubcase, dt):
        RealScalarTableArray.__init__(self, data_code, is_sort1, isubcase, dt)

    def write_f06(self, f06_file, header=None, page_stamp='PAGE %s',
                  page_num=1, is_mag_phase=False, is_sort1=True):
        if header is None:
            header = []
        if self.nonlinear_factor not in (None, np.nan):
            return self._write_f06_transient(header, page_stamp, page_num, f06_file,
                                             is_mag_phase=is_mag_phase, is_sort1=is_sort1)
        words = ['                                              T E M P E R A T U R E   V E C T O R\n',
                 ' \n',
                 '      POINT ID.   TYPE      ID   VALUE     ID+1 VALUE     ID+2 VALUE     ID+3 VALUE     ID+4 VALUE     ID+5 VALUE\n']
        return self._write_f06_block(words, header, page_stamp, page_num, f06_file, write_words=False)

    def _write_f06_transient(self, header, page_stamp, page_num=1, f06_file=None,
                             is_mag_phase=False, is_sort1=True):
        words = ['                                              T E M P E R A T U R E   V E C T O R\n',
                 ' \n',
                 '      POINT ID.   TYPE      ID   VALUE     ID+1 VALUE     ID+2 VALUE     ID+3 VALUE     ID+4 VALUE     ID+5 VALUE\n']
        return self._write_f06_transient_block(words, header, page_stamp, page_num, f06_file, write_words=False)
