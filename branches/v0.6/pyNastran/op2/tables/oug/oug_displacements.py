## GNU Lesser General Public License
##
## Program pyNastran - a python interface to NASTRAN files
## Copyright (C) 2011-2012  Steven Doyle, Al Danial
##
## Authors and copyright holders of pyNastran
## Steven Doyle <mesheb82@gmail.com>
## Al Danial    <al.danial@gmail.com>
##
## This file is part of pyNastran.
##
## pyNastran is free software: you can redistribute it and/or modify
## it under the terms of the GNU Lesser General Public License as published by
## the Free Software Foundation, either version 3 of the License, or
## (at your option) any later version.
##
## pyNastran is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
##
## You should have received a copy of the GNU Lesser General Public License
## along with pyNastran.  If not, see <http://www.gnu.org/licenses/>.
##
from struct import pack
from pyNastran.op2.resultObjects.tableObject import (TableObject,
                                                     ComplexTableObject)



def make_pack_form(data):
    N = 0
    n = 0
    Form = ''
    old = None
    for d in data:
        if isinstance(d, str):
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
            Form += form
            N = n
        else:
            N += n
        old = d
        fold = f
    if N:
        form = str(N) + f
        Form += form
    return form

class DisplacementObject(TableObject):  # approach_code=1, thermal=0
    def __init__(self, data_code, is_sort1, isubcase, dt=None):
        TableObject.__init__(self, data_code, is_sort1, isubcase, dt)

    def write_matlab(self, isubcase, f, is_mag_phase=False):
        name = 'displacements'
        if self.nonlinear_factor is None:
            return self._write_matlab(name, isubcase, f)
        else:
            return self._write_matlab_transient(name, isubcase, f)

    def write_table2(self, f, packing=False):
        i = -2
        marker1 = [4,  0, 4]
        marker2 = [4,  i, 4]
        marker3 = [4,  0, 4]
        marker = marker1 + marker2 + marker3
        if packing:
            npack = len(marker)
            p = pack('%ii' % npack, *marker)
            f.write(p)
        else:
            f.write(str(marker)+'\n')

        device_code = 1
        t_code = 1
        subcase = 1
        approach_code = 1
        analysis_code = 10 * approach_code
        plot_code = 1
        table_code = 1  # OUG

        random_code = 0    # 8
        num_wide = 4       # 10
        acoustic_flag = 1  # 13
        thermal = 0        # 23

        #               # 0      1                 2
        format_code = 0 # real, real/imag-complex, mag/phase

        is_sort2   = 0 # False
        is_complex = 0 # False
        is_random  = 0 # False

        bit = 4 * is_sort2 + 2 * is_complex + is_random
        #0 SORT1 Real No
        #1 SORT1 Complex No
        #2 SORT2 Real No
        #3 SORT2 Complex No
        #4 SORT1 Real Yes
        #5 SORT2 Real Yes

        data_form = ['i'] * 28
        data = [0] * 28

        data[0] = 'OUGV1   '
        data_form[0] = '8s'

        data[1] = analysis_code # a_code
        data[2] = t_code # table_code???
        data[4] = subcase

        lsdvmn = 1
        if t_code == 1:
            if approach_code == 1:
                data[5] = lsdvmn
            else:
                raise NotImplementedError(approach_code)
        else:
            raise NotImplementedError(t_code)

        #data[8] = random_code
        #data[10] = num_wide
        #data[13] = acoustic_flag
        data[23] = thermal

        #data[24] = 0
        self.title = ''
        self.subtitle = ''
        self.label = ''
        data_form[25] = '32s'
        data_form[26] = '32s'
        data_form[27] = '32s'
        data[25] = '%-32s' % self.title
        data[26] = '%-32s' % self.subtitle
        data[27] = '%-32s' % self.label
        data = data[:27]
        data_form = data_form[:27]

        form = '8s'

        stress = 0
        if 0:
            ## random code
            self.add_data_parameter(data, 'randomCode', 'i', 8, False)
            ## format code
            self.add_data_parameter(data, 'format_code', 'i', 9, False)
            ## number of words per entry in record; .. note:: is this needed for this table ???
            self.add_data_parameter(data, 'num_wide', 'i', 10, False)
            ## acoustic pressure flag
            self.add_data_parameter(data, 'acousticFlag', 'f', 13, False)
            ## thermal flag; 1 for heat transfer, 0 otherwise
            self.add_data_parameter(data, 'thermal', 'i', 23, False)
            self.isFlipped = False

            eidDevice = self.get_values(data, 'i', 5)
            floatVal = self.get_values(data, 'f', 5)


        #data = ['OUG     ', analysis_code, plot_code, heat, stress]
        form = '8s4i'
        if packing:
            form = make_pack_form(data)
            p = pack(''.join(data_form), *data)
            f.write(p)
        else:
            f.write(str(data)+'\n')

        n = 50
        i = [4, 0, 4, n, 4, 0, 4]

        ni = len(i)
        form = '%ii' % ni
        if packing:
            p = pack(form, *i)
            f.write(p)
        else:
            f.write(str(i)+'\n')
        return device_code


    def write_op2(self, header, pageStamp, f, is_mag_phase=False, packing=False):
        if self.nonlinear_factor is not None:
            return self._write_op2_transient(header, pageStamp, pageNum, f)

        i = -1
        marker1 = [4,  0, 4]
        marker2 = [4,  i, 4]
        marker3 = [4,  0, 4]
        marker = marker1 + marker2 + marker3
        if packing:
            npack = len(marker)
            p = pack('%ii' % npack, *marker)
            f.write(p)
        else:
            f.write(str(marker)+'\n')

        device_code = self.write_table2(f, packing=packing)
        return self._write_op2_block(f, header, device_code, packing=packing)


    def write_f06(self, header, pageStamp, pageNum=1, f=None, is_mag_phase=False):
        if self.nonlinear_factor is not None:
            return self._write_f06_transient(header, pageStamp, pageNum, f)
        words = ['                                             D I S P L A C E M E N T   V E C T O R\n',
                 ' \n',
                 '      POINT ID.   TYPE          T1             T2             T3             R1             R2             R3\n']
        words += self.get_table_marker()
        return self._write_f06_block(words, header, pageStamp, pageNum, f)

    def _write_f06_transient(self, header, pageStamp, pageNum=1, f=None, is_mag_phase=False):
        words = ['                                             D I S P L A C E M E N T   V E C T O R\n',
                 ' \n',
                 '      POINT ID.   TYPE          T1             T2             T3             R1             R2             R3\n']
        words += self.get_table_marker()
        return self._write_f06_transient_block(words, header, pageStamp, pageNum, f)


class ComplexDisplacementObject(ComplexTableObject):  # approach_code=1, sort_code=0, thermal=0
    def __init__(self, data_code, is_sort1, isubcase, dt=None):
        ComplexTableObject.__init__(self, data_code, is_sort1, isubcase, dt)

    def write_matlab(self, isubcase, f=None, is_mag_phase=False):
        name = 'displacements'
        if self.nonlinear_factor is None:
            return self._write_matlab(name, isubcase, f, is_mag_phase)
        else:
            return self._write_matlab_transient(name, isubcase, f, is_mag_phase)

    def write_f06(self, header, pageStamp, pageNum=1, f=None, is_mag_phase=False):
        if self.nonlinear_factor is not None:
            return self._write_f06_transient(header, pageStamp, pageNum, f, is_mag_phase)

        words = ['                                       C O M P L E X   D I S P L A C E M E N T   V E C T O R\n']
        return self._write_f06_block(words, header, pageStamp, pageNum, f, is_mag_phase)

    def _write_f06_transient(self, header, pageStamp, pageNum=1, f=None, is_mag_phase=False):
        words = ['                                       C O M P L E X   D I S P L A C E M E N T   V E C T O R\n']
        return self._write_f06_transient_block(words, header, pageStamp, pageNum, f, is_mag_phase)