import numpy as np
from pyNastran.op2.result_objects.table_object import RealTableArray, ComplexTableArray


class RealContactForcesArray(RealTableArray):
    def __init__(self, data_code, is_sort1, isubcase, dt):
        RealTableArray.__init__(self, data_code, is_sort1, isubcase, dt)

    def write_f06(self, f06_file, header=None, page_stamp='PAGE %s',
                  page_num: int=1, is_mag_phase: bool=False, is_sort1: bool=True):
        if header is None:
            header = []
        #' \n',
        #'      POINT ID.   TYPE          T1             T2             T3             R1             R2             R3\n']
        if self.table_name in ['OQGCF1']:
            # contact
            words = ['                               C O N T A C T  F O R C E S\n', ]
            pass
        elif self.table_name in ['OQGGF1']:
            # glue
            words = ['                               G L U E  F O R C E S\n', ]
            pass
        #elif self.table_name in ['OQGATO1', 'OQGATO2']:
            #words += ['                                                 ( AUTO-CORRELATION FUNCTION )']
        #elif self.table_name in ['OQGPSD1', 'OQGPSD2']:
            #words += ['                                             ( POWER SPECTRAL DENSITY FUNCTION )']
        #elif self.table_name in ['OQGRMS1', 'OQGRMS2']:
            #words += ['                                                     ( ROOT MEAN SQUARE )']
        #elif self.table_name in ['OQGCRM1', 'OQGCRM2']:
            #words += ['                                               ( CUMULATIVE ROOT MEAN SQUARE )']
        #elif self.table_name in ['OQGNO1', 'OQGNO2']:
            #words += ['                                                 ( NUMBER OF ZERO CROSSINGS )']
        else:  # pragma: no cover
            raise NotImplementedError(self.table_name)

        #words += self.get_table_marker()
        write_words = True
        if self.nonlinear_factor not in (None, np.nan):
            return self._write_f06_transient_block(
                words, header, page_stamp, page_num, f06_file, write_words,
                is_mag_phase=is_mag_phase, is_sort1=is_sort1)
        return self._write_f06_block(
            words, header, page_stamp, page_num, f06_file, write_words,
            is_mag_phase=is_mag_phase, is_sort1=is_sort1)


class ComplexSPCForcesArray(ComplexTableArray):
    def __init__(self, data_code, is_sort1, isubcase, dt):
        ComplexTableArray.__init__(self, data_code, is_sort1, isubcase, dt)

    def write_f06(self, f06_file, header=None, page_stamp='PAGE %s',
                  page_num: int=1, is_mag_phase: bool=False, is_sort1: bool=True):
        if header is None:
            header = []
        words = ['                         C O M P L E X   F O R C E S   O F   S I N G L E   P O I N T   C O N S T R A I N T\n']
        return self._write_f06_transient_block(
            words, header, page_stamp, page_num, f06_file, is_mag_phase, is_sort1)
