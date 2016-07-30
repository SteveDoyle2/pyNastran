from __future__ import print_function
from six import string_types
#from struct import pack


from pyNastran.op2.result_objects.table_object import RealTableArray, ComplexTableArray
#RealTableObject, ComplexTableObject


#def write_block(f, fascii)
def make_pack_form(data):
    N = 0
    n = 0
    #Form = ''
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


class RealDisplacementArray(RealTableArray):
    def __init__(self, data_code, is_sort1, isubcase, dt):
        RealTableArray.__init__(self, data_code, is_sort1, isubcase, dt)

    def write_f06(self, f, header=None, page_stamp='PAGE %s', page_num=1, is_mag_phase=False, is_sort1=True):
        if header is None:
            header = []
        words = ['                                             D I S P L A C E M E N T   V E C T O R\n', ]
                 #' \n',
                 #'      POINT ID.   TYPE          T1             T2             T3             R1             R2             R3\n']
        #words += self.get_table_marker()
        write_words = True
        if self.nonlinear_factor is not None:
            return self._write_f06_transient_block(words, header, page_stamp, page_num, f, write_words,
                                                   is_mag_phase=is_mag_phase, is_sort1=is_sort1)
        return self._write_f06_block(words, header, page_stamp, page_num, f, write_words,
                                         is_mag_phase=is_mag_phase, is_sort1=is_sort1)


class ComplexDisplacementArray(ComplexTableArray):
    def __init__(self, data_code, is_sort1, isubcase, dt):
        ComplexTableArray.__init__(self, data_code, is_sort1, isubcase, dt)

    def write_f06(self, f, header=None, page_stamp='PAGE %s', page_num=1, is_mag_phase=False, is_sort1=True):
        if header is None:
            header = []
        words = ['                                       C O M P L E X   D I S P L A C E M E N T   V E C T O R\n']
        return self._write_f06_transient_block(words, header, page_stamp, page_num, f,
                                               is_mag_phase=is_mag_phase, is_sort1=is_sort1)


#class RealDisplacement(RealTableObject):  # approach_code=1, thermal=0
    #def __init__(self, data_code, is_sort1, isubcase, dt):
        #RealTableObject.__init__(self, data_code, is_sort1, isubcase, dt)

    #def write_f06(self, f, header=None, page_stamp='PAGE %s', page_num=1, is_mag_phase=False, is_sort1=True):
        #if header is None:
            #header = []
        #words = ['                                             D I S P L A C E M E N T   V E C T O R\n',
                 #' \n',
                 #'      POINT ID.   TYPE          T1             T2             T3             R1             R2             R3\n']
        #words += self.get_table_marker()

        #if self.nonlinear_factor is not None:
            #return self._write_f06_transient_block(words, header, page_stamp, page_num, f,
                                                   #is_mag_phase=is_mag_phase, is_sort1=is_sort1)
        #return self._write_f06_block(words, header, page_stamp, page_num, f,
                                     #is_mag_phase=is_mag_phase, is_sort1=is_sort1)

    #def _write_table_3(self, f, fascii, itable):
        #aCode = 1
        #tCode = 1
        #approach_code = 0 #self.approach_code
        #table_code = self.table_code
        #isubcase = self.isubcase
        #lsdvmn = 1
        #random_code = 0
        #format_code = 1
        #num_wide = self.num_wide
        #acoustic_flag = 0
        #thermal = 0
        #title = '%-128s' % self.title
        #subtitle = '%-128s' % self.subtitle
        #label = '%-128s' % self.label
        #ftable3 = '24i 128s 128s 128s'
        #oCode = 0

        #table3 = [
            #aCode, tCode, 0, isubcase, lsdvmn,
            #0, 0, random_code, format_code, num_wide,
            #oCode, acoustic_flag, 0, 0, 0,
            #0, 0, 0, 0, 0,
            #0, thermal, thermal, 0, Title,
            #subtitle, label, ]

        #write_markers(f, fascii, '%s header3a' % self.table_name, [-3, 1, 0])
        #write_markers(f, fascii, '%s header3b' % self.table_name, [146])

        #data = [584] + table3 + [584]
        #fmt = 'i' + ftable3 + 'i'
        #f.write(pack(fascii, '%s header 3c' % self.table_name, fmt, data))

    #def write_op2(self, f, fascii, is_mag_phase=False):
        #fascii.write('%s.write_op2\n' % self.__class__.__name__)
        #self._write_table_header(f, fascii)
        ##recordi = ['OUG1    ', 2, 26, 14, 0, 1]
        ##[subtable_name, month=2, day=26, year=2014, zero=0, one=1]
        #subtable_name = 'OUG1    '
        #self._write_table3(f, fascii, -3)

        #if self.nonlinear_factor is not None:
            #return self._write_op2_transient_block(f, fascii)
        #return self._write_op2_block(f, fascii)

    #def _write_op2_block(self, f, fascii):
        #nnodes = len(self.translations)
        #nwords = nnodes * self.num_wide
        #nbytes = nwords * 4
        #print('nnodes=%s nbytes=%s' % (nnodes, nbytes))
        ##nbytes = 16386
        ##assert nbytes >= 16386, nbytes
        #print('nbytes=%s' % nbytes)
        #header = [
            #4, 0, 4,
            #4, nwords, 4,
            #nbytes #recordi, nbytes*4,
        #]
        #f.write(pack(fascii, 'header4', '7i', header))

        #data = []
        #fmt = '2i 6f'
        #device_code = self.device_code
        #for nodeID, translation in sorted(iteritems(self.translations)):
            #rotation = self.rotations[nodeID]
            #grid_type = self.gridTypes[nodeID]
            #if grid_type == 'G':
                #grid_type = 1
            #else:
                #raise RuntimeError(grid_type)

            #(dx, dy, dz) = translation
            #(rx, ry, rz) = rotation
            #vals = [dx, dy, dz, rx, ry, rz]
            #data = [device_code + 10*nodeID, grid_type, dx, dy, dz, rx, ry, rz]
            #f.write(pack(fascii, 'nid, gridType, dx, dy, dz, rx, ry, rz', fmt, data))
        #f.write(pack(fascii, 'closer', 'i', nbytes))
