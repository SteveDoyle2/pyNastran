from __future__ import print_function
from six import string_types
import numpy as np


from pyNastran.op2.result_objects.table_object import RealTableArray, ComplexTableArray


#def write_block(f, fascii)
def make_pack_form(data):
    N = 0
    n = 0
    #Form = ''
    fold = None
    old = None
    for d in data:
        if isinstance(d, string_types):
            n = len(d)
            f = 's'
        elif isinstance(d, int):
            n = 4
            f = 'i'
        elif isinstance(d, float):
            n = 4
            f = 'f'
        else:
            raise NotImplementedError(type(d))
        if old and f != fold:
            form = str(N) + fold
            #Form += form
            N = n
        else:
            N += n
        old = d
        fold = f
    if N:
        form = str(N) + f
        #Form += form
    return form


class RealDisplacementArray(RealTableArray):
    def __init__(self, data_code, is_sort1, isubcase, dt):
        RealTableArray.__init__(self, data_code, is_sort1, isubcase, dt)

    def write_f06(self, f06_file, header=None, page_stamp='PAGE %s',
                  page_num=1, is_mag_phase=False, is_sort1=True):
        if header is None:
            header = []
        words = ['                                             D I S P L A C E M E N T   V E C T O R\n', ]
        #' \n',
        #'      POINT ID.   TYPE          T1             T2             T3             R1             R2             R3\n']
        if self.table_name in ['OUGV1', 'OUGV2', 'BOUGV1', 'OUPV1']:
            pass
        elif self.table_name in ['ROUGV1', 'ROUGV2']:
            words += ['                                                (RELATIVE TO ENFORCED MOTION INPUT)']
        elif self.table_name in ['OUGATO1', 'OUGATO2']:
            words += ['                                                 ( AUTO-CORRELATION FUNCTION )']
        elif self.table_name in ['OUGPSD1', 'OUGPSD2']:
            words += ['                                             ( POWER SPECTRAL DENSITY FUNCTION )']
        elif self.table_name in ['OUGRMS1', 'OUGRMS2']:
            words += ['                                                     ( ROOT MEAN SQUARE )']
        elif self.table_name in ['OUGCRM1', 'OUGCRM2']:
            words += ['                                               ( CUMULATIVE ROOT MEAN SQUARE )']
        elif self.table_name in ['OUGNO1', 'OUGNO2']:
            words += ['                                                 ( NUMBER OF ZERO CROSSINGS )']
        elif self.table_name in ['OCRUG']:
            words += ['                                                 ( OCRUG??? )']
        else:
            raise NotImplementedError('table_name=%r' % self.table_name)
        #words += self.get_table_marker()
        write_words = True
        if self.nonlinear_factor not in (None, np.nan):
            return self._write_f06_transient_block(words, header, page_stamp, page_num, f06_file, write_words,
                                                   is_mag_phase=is_mag_phase, is_sort1=is_sort1)
        return self._write_f06_block(words, header, page_stamp, page_num, f06_file, write_words,
                                     is_mag_phase=is_mag_phase, is_sort1=is_sort1)


class ComplexDisplacementArray(ComplexTableArray):
    def __init__(self, data_code, is_sort1, isubcase, dt):
        ComplexTableArray.__init__(self, data_code, is_sort1, isubcase, dt)

    def write_f06(self, f06_file, header=None, page_stamp='PAGE %s', page_num=1, is_mag_phase=False, is_sort1=True):
        if header is None:
            header = []
        words = ['                                       C O M P L E X   D I S P L A C E M E N T   V E C T O R\n']
        return self._write_f06_transient_block(words, header, page_stamp, page_num, f06_file,
                                               is_mag_phase=is_mag_phase, is_sort1=is_sort1)
