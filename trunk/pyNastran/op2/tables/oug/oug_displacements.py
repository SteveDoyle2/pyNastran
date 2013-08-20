from pyNastran.op2.resultObjects.tableObject import TableObject, ComplexTableObject


class DisplacementObject(TableObject):  # approach_code=1, thermal=0
    def __init__(self, data_code, is_sort1, isubcase, dt=None, read_mode=0):
        TableObject.__init__(self, data_code, is_sort1, isubcase, dt, read_mode)

    def write_matlab(self, isubcase, f, is_mag_phase=False):
        name = 'displacements'
        if self.nonlinear_factor is None:
            return self._write_matlab(name, isubcase, f)
        else:
            return self._write_matlab_transient(name, isubcase, f)

    def write_f06(self, header, pageStamp, f, pageNum=1, is_mag_phase=False):
        if self.nonlinear_factor is not None:
            return self._write_f06_transient(header, pageStamp, f, pageNum)
        words = ['                                             D I S P L A C E M E N T   V E C T O R\n',
                 ' \n',
                 '      POINT ID.   TYPE          T1             T2             T3             R1             R2             R3\n']
        words += self.get_table_marker()
        return self._write_f06_block(words, header, pageStamp, f, pageNum)

    def _write_f06_transient(self, header, pageStamp, f, pageNum=1, is_mag_phase=False):
        words = ['                                             D I S P L A C E M E N T   V E C T O R\n',
                 ' \n',
                 '      POINT ID.   TYPE          T1             T2             T3             R1             R2             R3\n']
        words += self.get_table_marker()
        return self._write_f06_transient_block(words, header, pageStamp, f, pageNum)


class ComplexDisplacementObject(ComplexTableObject):  # approach_code=1, sort_code=0, thermal=0
    def __init__(self, data_code, is_sort1, isubcase, dt=None, read_mode=0):
        ComplexTableObject.__init__(self, data_code, is_sort1, isubcase, dt, read_mode)

    def write_matlab(self, isubcase, f, is_mag_phase=False):
        name = 'displacements'
        if self.nonlinear_factor is None:
            return self._write_matlab(name, isubcase, f, is_mag_phase)
        else:
            return self._write_matlab_transient(name, isubcase, f, is_mag_phase)

    def write_f06(self, header, pageStamp, f, pageNum=1, is_mag_phase=False):
        if self.nonlinear_factor is not None:
            return self._write_f06_transient(header, pageStamp, f, pageNum, is_mag_phase)

        words = ['                                       C O M P L E X   D I S P L A C E M E N T   V E C T O R\n']
        return self._write_f06_block(words, header, pageStamp, f, pageNum, is_mag_phase)

    def _write_f06_transient(self, header, pageStamp, f, pageNum=1, is_mag_phase=False):
        words = ['                                       C O M P L E X   D I S P L A C E M E N T   V E C T O R\n']
        return self._write_f06_transient_block(words, header, pageStamp, f, pageNum, is_mag_phase)