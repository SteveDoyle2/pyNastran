import numpy as np
from pyNastran.op2.result_objects.table_object import RealTableArray, ComplexTableArray


class RealAccelerationArray(RealTableArray):
    def __init__(self, data_code, is_sort1, isubcase, dt):
        RealTableArray.__init__(self, data_code, is_sort1, isubcase, dt)

    def write_f06(self, f06_file, header=None, page_stamp='PAGE %s', page_num=1, is_mag_phase=False, is_sort1=True):
        if header is None:
            header = []
        words = ['                                             A C C E L E R A T I O N   V E C T O R\n', ]
        #' \n',
        #'      POINT ID.   TYPE          T1             T2             T3             R1             R2             R3\n']
        if self.table_name in ['OUGV1', 'OUGV2', 'OUPV1', 'BOUGV1']:
            pass
        elif self.table_name in ['OUXY1', 'OUXY2']:
            words = ['                                       A C C E L E R A T I O N   V E C T O R   (SOLUTION SET)']
        elif self.table_name in ['ROUGV1', 'ROUGV2']:
            words += ['                                                (RELATIVE TO ENFORCED MOTION INPUT)']
        elif self.table_name in ['OAGATO1', 'OAGATO2']:
            words += ['                                                 ( AUTO-CORRELATION FUNCTION )']
        elif self.table_name in ['OUGPSD1', 'OUGPSD2', 'OAGPSD1', 'OAGPSD2']:
            words += ['                                             ( POWER SPECTRAL DENSITY FUNCTION )']
        elif self.table_name in ['OUGRMS1', 'OUGRMS2', 'OAGRMS1', 'OAGRMS2']:
            words += ['                                                     ( ROOT MEAN SQUARE )']
        elif self.table_name in ['OUGCRM1', 'OUGCRM2', 'OAGCRM1', 'OAGCRM2']:
            words += ['                                               ( CUMULATIVE ROOT MEAN SQUARE )']
        elif self.table_name in ['OUGNO1', 'OUGNO2', 'OAGNO1', 'OAGNO2']:
            words += ['                                                 ( NUMBER OF ZERO CROSSINGS )']
        else:
            raise NotImplementedError(self.table_name)

                 #' \n',
                 #'      POINT ID.   TYPE          T1             T2             T3             R1             R2             R3\n']
        #words += self.get_table_marker()
        if self.nonlinear_factor not in (None, np.nan):
            return self._write_f06_transient_block(words, header, page_stamp, page_num, f06_file, write_words=True,
                                                   is_mag_phase=is_mag_phase, is_sort1=is_sort1)
        return self._write_f06_block(words, header, page_stamp, page_num, f06_file, write_words=True)


class ComplexAccelerationArray(ComplexTableArray):
    def __init__(self, data_code, is_sort1, isubcase, dt):
        ComplexTableArray.__init__(self, data_code, is_sort1, isubcase, dt)

    def write_f06(self, f06_file, header=None, page_stamp='PAGE %s',
                  page_num=1, is_mag_phase=False, is_sort1=True):
        if header is None:
            header = []
        words = ['                                       C O M P L E X   A C C E L E R A T I O N   V E C T O R\n']
        #words += self.get_table_marker()
        return self._write_f06_transient_block(words, header, page_stamp, page_num, f06_file,
                                               is_mag_phase=is_mag_phase, is_sort1=is_sort1)
