from pyNastran.op2.resultObjects.tableObject import TableObject, ComplexTableObject


class AccelerationObject(TableObject):  # approach_code=11, thermal=0
    def __init__(self, data_code, is_sort1, isubcase, dt=None, read_mode=0):
        TableObject.__init__(self, data_code, is_sort1, isubcase, dt, read_mode)

    def write_f06(self, header, pageStamp, f, pageNum=1, is_mag_phase=False):
        if self.nonlinear_factor is not None:
            return self._write_f06_transient(header, pageStamp, f, pageNum)
        words = ['                                             A C C E L E R A T I O N   V E C T O R\n',
                 ' \n',
                 '      POINT ID.   TYPE          T1             T2             T3             R1             R2             R3\n']
        words += self.get_table_marker()
        return self._write_f06_block(words, header, pageStamp, f, pageNum)

    def _write_f06_transient(self, header, pageStamp, f, pageNum=1, is_mag_phase=False):
        words = ['                                             A C C E L E R A T I O N   V E C T O R\n',
                 ' \n',
                 '      POINT ID.   TYPE          T1             T2             T3             R1             R2             R3\n']
        words += self.get_table_marker()
        return self._write_f06_transient_block(words, header, pageStamp, f, pageNum)


class ComplexAccelerationObject(ComplexTableObject):  # table_code=11, approach_code=???
    def __init__(self, data_code, is_sort1, isubcase, dt=None, read_mode=0):
        ComplexTableObject.__init__(self, data_code, is_sort1, isubcase, dt, read_mode)

    def write_f06(self, header, pageStamp, f, pageNum=1, is_mag_phase=False):
        if self.nonlinear_factor is not None:
            return self._write_f06_transient(header, pageStamp, f, pageNum, is_mag_phase)

        words = ['                                       C O M P L E X   A C C E L E R A T I O N   V E C T O R\n']
        return self._write_f06_block(words, header, pageStamp, f, pageNum, is_mag_phase)

    def _write_f06_transient(self, header, pageStamp, f, pageNum=1, is_mag_phase=False):
        words = ['                                       C O M P L E X   A C C E L E R A T I O N   V E C T O R\n']
        return self._write_f06_transient_block(words, header, pageStamp, f, pageNum, is_mag_phase)