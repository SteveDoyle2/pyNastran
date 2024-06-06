
from __future__ import annotations
#import copy
from struct import pack
#from itertools import count
import warnings
#from typing import TextIO

import numpy as np

from pyNastran.bdf import MAX_32_BIT_INT
#from pyNastran.op2.result_objects.op2_objects import (
    #ScalarObject, get_sort_node_sizes, set_as_sort1,
    #NULL_GRIDTYPE, SORT1_TABLES, SORT2_TABLES)

#from pyNastran.f06.f06_formatting import (
    #write_floats_13e, write_floats_13e_long,
    #write_imag_floats_13e, write_float_12e)
from pyNastran.op2.errors import SixtyFourBitError
from pyNastran.op2.op2_interface.write_utils import view_dtype, view_idtype_as_fdtype
#from pyNastran.utils.numpy_utils import integer_types, float_types
#from pyNastran.op2.writer.utils import fix_table3_types
#from pyNastran.op2.tables.oes_stressStrain.real.oes_objects import (
    #set_static_case, set_modal_case, set_transient_case,
    #set_freq_case, set_complex_modes_case)

from pyNastran.op2.result_objects.table_object import RealTableArray # , ComplexTableArray


class RealBoltArray(RealTableArray):
    def __init__(self, data_code, is_sort1, isubcase, dt):
        RealTableArray.__init__(self, data_code, is_sort1, isubcase, dt)

    def write_f06(self, f06_file, header=None, page_stamp='PAGE %s',
                  page_num: int=1, is_mag_phase: bool=False, is_sort1: bool=True):
        if header is None:
            header = []

        if self.table_name in ['OBOLT1']:
            words = [
                '                                                         B O L T   R E S U L T S\n',
                '        ELEMENT ID    AXIAL FORCE    SHEAR FORCE-1    SHEAR FORCE-2    BENDING MOMENT-1    BENDING MOMENT-2    AXIAL PRELOAD STRAIN\n'
            ]
            write_words = False
        else:  # pragma: no cover
            raise NotImplementedError(self.table_name)

        #words += self.get_table_marker()
        write_words = False
        if self.nonlinear_factor not in (None, np.nan):
            return self._write_f06_transient_block(words, header, page_stamp, page_num, f06_file, write_words,
                                                   is_mag_phase=is_mag_phase, is_sort1=is_sort1)
        return self._write_f06_block(words, header, page_stamp, page_num, f06_file, write_words,
                                     is_mag_phase=is_mag_phase, is_sort1=is_sort1)

    def write_op2(self, op2_file, fascii, itable: int, new_result,
                  date, is_mag_phase: bool=False, endian: str='>'):
        """writes an OP2"""
        import inspect
        allowed_tables = [
            'OBOLT1',
        ]

        assert self.table_name in allowed_tables, self.table_name

        frame = inspect.currentframe()
        call_frame = inspect.getouterframes(frame, 2)
        fascii.write('%s.write_op2: %s\n' % (self.__class__.__name__, call_frame[1][3]))

        if itable == -1:
            self._write_table_header(op2_file, fascii, date)
            itable = -3

        #print('nonlinear_factor =', self.nonlinear_factor)
        if not self.is_sort1:
            raise NotImplementedError('SORT2')
        #op2_format = endian + b'2i6f'

        node = self.node_gridtype[:, 0]
        #gridtype = self.node_gridtype[:, 1]

        # table 4 info
        #ntimes = self.data.shape[0]
        nnodes = self.data.shape[1]
        nnodes_device = self.node_gridtype[:, 0] * 10 + self.device_code

        nnodes = len(node)
        max_id = node.max()
        if max_id > MAX_32_BIT_INT:
            raise SixtyFourBitError(f'64-bit OP2 writing is not supported; max id={max_id}')

        fdtype = self.data.dtype
        if self.size == fdtype.itemsize:
            pass
        else:
            warnings.warn(f'downcasting {self.class_name}...this is buggy')
            #idtype = np.int32(1)
            fdtype = np.float32(1.0)

        node_floats = view_idtype_as_fdtype(nnodes_device, fdtype)

        #print(node_gridtype_floats)
        #node_gridtype_floats = nodedevice_gridtype.view(fdtype) # .reshape(nnodes, 2)

        #format_table4_1 = Struct(self._endian + b'15i')
        #format_table4_2 = Struct(self._endian + b'3i')

        #(2+6) => (node_id, t1i, t2i, t3i, r1i, r2i, r3i)
        ntotal = nnodes * 7

        #print('shape = %s' % str(self.data.shape))
        #assert nnodes > 1, nnodes
        assert ntotal > 1, ntotal

        unused_device_code = self.device_code
        fascii.write('  ntimes = %s\n' % self.ntimes)

        #fmt = '%2i %6f'
        #print('ntotal=%s' % (ntotal))
        for itime in range(self.ntimes):
            self._write_table_3(op2_file, fascii, new_result, itable, itime)

            # record 4
            itable -= 1
            header = [4, itable, 4,
                      4, 1, 4,
                      4, 0, 4,
                      4, ntotal, 4,
                      4*ntotal]
            op2_file.write(pack(b'%ii' % len(header), *header))
            fascii.write('r4 [4, 0, 4]\n')
            fascii.write(f'r4 [4, {itable:d}, 4]\n')
            fascii.write(f'r4 [4, {4*ntotal:d}, 4]\n')

            datai = view_dtype(self.data[itime, :, :], fdtype)
            node_gridtype_data = np.column_stack([node_floats, datai])
            op2_file.write(node_gridtype_data)
            assert ntotal == node_gridtype_data.size
            #print(self.data.shape, ntotal)

            itable -= 1
            header = [4 * ntotal,]
            op2_file.write(pack(b'i', *header))
            fascii.write('footer = %s\n' % header)
            new_result = False
        return itable
