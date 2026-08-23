import numpy as np
from pyNastran.op2.result_objects.table_object import RealTableArray, ComplexTableArray
from pyNastran.op2.result_objects.scalar6_table_object import RealScalarTableArray


class RealLoadVectorArray(RealTableArray):  # table_code=2, sort_code=0, thermal=0

    def __init__(self, data_code, is_sort1, isubcase, dt):
        RealTableArray.__init__(self, data_code, is_sort1, isubcase, dt)

    def write_f06(self, f06_file, header=None, page_stamp='PAGE %s',
                  page_num: int=1, is_mag_phase: bool=False, is_sort1: bool=True):
        if header is None:
            header = []
        words = ['                                                     L O A D   V E C T O R\n', ]
        #words += self.get_table_marker()
        write_words = True
        if self.nonlinear_factor not in (None, np.nan):
            return self._write_f06_transient_block(
                words, header, page_stamp, page_num, f06_file, write_words,
                is_mag_phase=is_mag_phase, is_sort1=is_sort1)
        return self._write_f06_block(
            words, header, page_stamp, page_num, f06_file, write_words,
            is_mag_phase=False, is_sort1=True
        )


class ComplexLoadVectorArray(ComplexTableArray):
    def __init__(self, data_code, is_sort1, isubcase, dt):
        ComplexTableArray.__init__(self, data_code, is_sort1, isubcase, dt)

    def write_f06(self, f06_file, header=None, page_stamp='PAGE %s',
                  page_num: int=1, is_mag_phase: bool=False, is_sort1: bool=True):
        if header is None:
            header = []
        words = ['                                               C O M P L E X   L O A D   V E C T O R\n', ]
        return self._write_f06_transient_block(
            words, header, page_stamp, page_num, f06_file, is_mag_phase, is_sort1)


class RealTemperatureVectorArray(RealScalarTableArray):
    def __init__(self, data_code, is_sort1, isubcase, dt):
        RealScalarTableArray.__init__(self, data_code, is_sort1, isubcase, dt)

    def write_f06(self, f06_file, header=None, page_stamp='PAGE %s',
                  page_num: int=1, is_mag_phase: bool=False, is_sort1: bool=True):
        if header is None:
            header = []
        words = [
            '                                              T E M P E R A T U R E   V E C T O R\n',
            ' \n',
            '      POINT ID.   TYPE      ID   VALUE     ID+1 VALUE     ID+2 VALUE     ID+3 VALUE     ID+4 VALUE     ID+5 VALUE\n'
        ]
        #words += self.get_table_marker()
        write_words = False
        if self.nonlinear_factor not in (None, np.nan):
            return self._write_f06_transient_block(
                words, header, page_stamp, page_num, f06_file, write_words,
                is_mag_phase=is_mag_phase, is_sort1=is_sort1)
        return self._write_f06_block(words, header, page_stamp, page_num, f06_file, write_words,
                                     is_mag_phase=is_mag_phase, is_sort1=is_sort1)

class RealThermalVelocityVectorArray(RealScalarTableArray):
    def __init__(self, data_code, is_sort1, isubcase, dt):
        RealScalarTableArray.__init__(self, data_code, is_sort1, isubcase, dt)

    def h5_table_dict(self) -> dict:
        from tables import Int64Col, Float64Col
        h5_table_dict = {
            'ID': Int64Col(pos=0),
            'X': Float64Col(pos=1),
            'Y': Float64Col(pos=2),
            'Z': Float64Col(pos=3),
            'RX': Float64Col(pos=4),
            'RY': Float64Col(pos=5),
            'RZ': Float64Col(pos=6),
            'DOMAIN_ID': Int64Col(pos=7),
        }
        return h5_table_dict

    def add_to_h5_array(self, arr,
                        ntime_nnode0: int, ntime_nnode1: int,
                        itime: int):
        arr["ID"][ntime_nnode0:ntime_nnode1] = self.node_gridtype[:, 0]
        arr["X"][ntime_nnode0:ntime_nnode1] = 0.
        arr["Y"][ntime_nnode0:ntime_nnode1] = 0.
        arr["Z"][ntime_nnode0:ntime_nnode1] = 0.
        arr["RX"][ntime_nnode0:ntime_nnode1] = 0.
        arr["RY"][ntime_nnode0:ntime_nnode1] = 0.
        arr["RZ"][ntime_nnode0:ntime_nnode1] = 0.

    def write_f06(self, f06_file, header=None, page_stamp='PAGE %s',
                  page_num: int=1, is_mag_phase: bool=False, is_sort1: bool=True):
        if header is None:
            header = []
        words = [
            '                                              THERMAL VELOCITY   V E C T O R\n',
            ' \n',
            '      POINT ID.   TYPE      ID   VALUE     ID+1 VALUE     ID+2 VALUE     ID+3 VALUE     ID+4 VALUE     ID+5 VALUE\n'
        ]
        #words += self.get_table_marker()
        write_words = False
        if self.nonlinear_factor not in (None, np.nan):
            return self._write_f06_transient_block(
                words, header, page_stamp, page_num, f06_file, write_words,
                is_mag_phase=is_mag_phase, is_sort1=is_sort1)
        return self._write_f06_block(
            words, header, page_stamp, page_num, f06_file, write_words,
            is_mag_phase=is_mag_phase, is_sort1=is_sort1)
