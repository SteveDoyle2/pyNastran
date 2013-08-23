from numpy import abs, angle
from pyNastran.op2.resultObjects.tableObject import TableObject, ComplexTableObject
from pyNastran.f06.f06_formatting import writeFloats13E


class LoadVectorObject(TableObject):  # table_code=2, sort_code=0, thermal=0

    def __init__(self, data_code, is_sort1, isubcase, dt, read_mode):
        TableObject.__init__(self, data_code, is_sort1, isubcase, dt, read_mode)

    def write_matlab(self, isubcase, f, is_mag_phase=False):
        name = 'loadVector'
        if self.nonlinear_factor is None:
            return self._write_matlab(name, isubcase, f)
        else:
            return self._write_matlab_transient(name, isubcase, f)

    def write_f06(self, header, pageStamp, f, pageNum=1, is_mag_phase=False):
        if self.nonlinear_factor is not None:
            return self._write_f06_transient(header, pageStamp, f, pageNum)

        words = header + [
            '                                                     L O A D   V E C T O R\n',
            ' \n',
            '      POINT ID.   TYPE          T1             T2             T3             R1             R2             R3\n']
        return self._write_f06_block(words, header, pageStamp, f, pageNum)

    def _write_f06_transient(self, header, pageStamp, f, pageNum=1):
        msg = []
        words = ['                                                     L O A D   V E C T O R\n',
                 ' \n',
                 '      POINT ID.   TYPE          T1             T2             T3             R1             R2             R3\n']
        return self._write_f06_transient_block(words, header, pageStamp, f, pageNum)

class ComplexLoadVectorObject(ComplexTableObject):  # table_code=11, approach_code=???
    def __init__(self, data_code, is_sort1, isubcase, dt, read_mode):
        ComplexTableObject.__init__(self, data_code, is_sort1, isubcase, dt, read_mode)

    def write_matlab(self, isubcase, f, is_mag_phase=False):
        name = 'loadVector'
        if self.nonlinear_factor is None:
            return self._write_matlab(name, isubcase, f, is_mag_phase)
        else:
            return self._write_matlab_transient(name, isubcase, f, is_mag_phase)

    def write_f06(self, header, pageStamp, f, pageNum=1, is_mag_phase=False):
        if self.nonlinear_factor is not None:
            return self._write_f06_transient(header, pageStamp, f, pageNum, is_mag_phase)
        words = header + [
            '                                               C O M P L E X   L O A D   V E C T O R\n',
            '                                                          (REAL/IMAGINARY)\n',
            ' \n',
            '      POINT ID.   TYPE          T1             T2             T3             R1             R2             R3\n']
        return self._write_f06_block(words, header, pageStamp, f, pageNum)

    def _write_f06_transient(self, header, pageStamp, f, pageNum=1, is_mag_phase=False):
        words = [
            '                                               C O M P L E X   L O A D   V E C T O R\n',
            '                                                          (REAL/IMAGINARY)\n',
            ' \n',
            '      POINT ID.   TYPE          T1             T2             T3             R1             R2             R3\n']
        return self._write_f06_transient_block(words, header, pageStamp, f, pageNum, is_mag_phase)


class ThermalVector(TableObject):
    def __init__(self, data_code, is_sort1, isubcase, dt, read_mode):
        TableObject.__init__(self, data_code, is_sort1, isubcase, dt, read_mode)

    def write_f06(self, header, pageStamp, f, pageNum=1, is_mag_phase=False):
        if self.nonlinear_factor is not None:
            return self._write_f06_transient(header, pageStamp, f, pageNum)

        msg = header + ['                                              T E M P E R A T U R E   V E C T O R\n',
                        ' \n',
                        '      POINT ID.   TYPE      ID   VALUE     ID+1 VALUE     ID+2 VALUE     ID+3 VALUE     ID+4 VALUE     ID+5 VALUE\n']
        return self._write_f06_block(words, header, pageStamp, f, pageNum)

    def _write_f06_transient(self, header, pageStamp, f, pageNum=1, is_mag_phase=False):
        msg = []
        words = ['                                              T E M P E R A T U R E   V E C T O R\n',
                 ' \n',
                 '      POINT ID.   TYPE      ID   VALUE     ID+1 VALUE     ID+2 VALUE     ID+3 VALUE     ID+4 VALUE     ID+5 VALUE\n']
        return self._write_f06_transient_block(words, header, pageStamp, f, pageNum)


class ThermalLoadVectorObject(ThermalVector):     # table_code=2, thermal=1
    def __init__(self, data_code, is_sort1, isubcase, dt, read_mode):
        ThermalVector.__init__(self, data_code, is_sort1, isubcase, dt, read_mode)


class ThermalVelocityVectorObject(ThermalVector):  # table_code=10, thermal=1
    def __init__(self, data_code, is_sort1, isubcase, dt, read_mode):
        ThermalVector.__init__(self, data_code, is_sort1, isubcase, dt, read_mode)
