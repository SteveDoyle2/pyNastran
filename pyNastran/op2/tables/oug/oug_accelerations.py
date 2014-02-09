from pyNastran.op2.resultObjects.tableObject import TableObject, ComplexTableObject


class AccelerationObject(TableObject):  # approach_code=11, thermal=0
    def __init__(self, data_code, is_sort1, isubcase, dt=None):
        TableObject.__init__(self, data_code, is_sort1, isubcase, dt)

    def write_matlab(self, isubcase, f=None, is_mag_phase=False):
        name = 'accelerations'
        if self.nonlinear_factor is None:
            return self._write_matlab(name, isubcase, f)
        else:
            return self._write_matlab_transient(name, isubcase, f)

    def write_f06(self, header, pageStamp, pageNum=1, f=None, is_mag_phase=False):
        if self.nonlinear_factor is not None:
            return self._write_f06_transient(header, pageStamp, pageNum, f)
        words = ['                                             A C C E L E R A T I O N   V E C T O R\n',
                 ' \n',
                 '      POINT ID.   TYPE          T1             T2             T3             R1             R2             R3\n']
        words += self.get_table_marker()
        return self._write_f06_block(words, header, pageStamp, pageNum, f)

    def _write_f06_transient(self, header, pageStamp, pageNum=1, f=None, is_mag_phase=False):
        words = ['                                             A C C E L E R A T I O N   V E C T O R\n',
                 ' \n',
                 '      POINT ID.   TYPE          T1             T2             T3             R1             R2             R3\n']
        words += self.get_table_marker()
        return self._write_f06_transient_block(words, header, pageStamp, pageNum, f)


class ComplexAccelerationObject(ComplexTableObject):  # table_code=11, approach_code=???
    def __init__(self, data_code, is_sort1, isubcase, dt=None):
        ComplexTableObject.__init__(self, data_code, is_sort1, isubcase, dt)

    def write_matlab(self, isubcase, f=None, is_mag_phase=False):
        name = 'accelerations'
        if self.nonlinear_factor is None:
            return self._write_matlab(name, isubcase, f, is_mag_phase)
        else:
            return self._write_matlab_transient(name, isubcase, f, is_mag_phase)

    def write_f06(self, header, pageStamp, pageNum=1, f=None, is_mag_phase=False):
        if self.nonlinear_factor is not None:
            return self._write_f06_transient(header, pageStamp, pageNum, f, is_mag_phase)

        words = ['                                       C O M P L E X   A C C E L E R A T I O N   V E C T O R\n']
        return self._write_f06_block(words, header, pageStamp, pageNum, f, is_mag_phase)

    def _write_f06_transient(self, header, pageStamp, pageNum=1, f=None, is_mag_phase=False):
        words = ['                                       C O M P L E X   A C C E L E R A T I O N   V E C T O R\n']
        return self._write_f06_transient_block(words, header, pageStamp, pageNum, f, is_mag_phase)