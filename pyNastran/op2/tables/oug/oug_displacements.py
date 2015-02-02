from struct import pack as spack, Struct as sStruct

def pack(fascii, msg, fmt, data):
    s = Struct(fascii, fmt)
    return s.pack(msg, data)

class Struct(object):
    def __init__(self, fascii, fmt):
        self.fascii = fascii
        self.s = sStruct(fmt)

    def pack(self, msg, data):
        self.fascii.write('%s = %s\n' % (msg, data))
        if isinstance(data, list):
            return self.s.pack(*data)
        return self.s.pack(data)

from pyNastran.op2.resultObjects.tableObject import RealTableArray, ComplexTableArray, RealTableObject, ComplexTableObject


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


class RealDisplacementArray(RealTableArray):
    def __init__(self, data_code, is_sort1, isubcase, dt):
        RealTableArray.__init__(self, data_code, is_sort1, isubcase, dt)

    def write_f06(self, header, pageStamp, page_num=1, f=None, is_mag_phase=False):
        words = ['                                             D I S P L A C E M E N T   V E C T O R\n', ]
                 #' \n',
                 #'      POINT ID.   TYPE          T1             T2             T3             R1             R2             R3\n']
        #words += self.get_table_marker()
        if self.nonlinear_factor is not None:
            return self._write_f06_transient_block(words, header, pageStamp, page_num, f)
        return self._write_f06_block(words, header, pageStamp, page_num, f)


class ComplexDisplacementArray(ComplexTableArray):
    def __init__(self, data_code, is_sort1, isubcase, dt):
        ComplexTableArray.__init__(self, data_code, is_sort1, isubcase, dt)

    def write_f06(self, header, pageStamp, page_num=1, f=None, is_mag_phase=False):
        words = ['                                       C O M P L E X   D I S P L A C E M E N T   V E C T O R\n']
        return self._write_f06_transient_block(words, header, pageStamp, page_num, f, is_mag_phase)


class RealDisplacement(RealTableObject):  # approach_code=1, thermal=0
    def __init__(self, data_code, is_sort1, isubcase, dt):
        RealTableObject.__init__(self, data_code, is_sort1, isubcase, dt)

    def write_f06(self, header, pageStamp, page_num=1, f=None, is_mag_phase=False):
        words = ['                                             D I S P L A C E M E N T   V E C T O R\n',
                 ' \n',
                 '      POINT ID.   TYPE          T1             T2             T3             R1             R2             R3\n']
        words += self.get_table_marker()

        if self.nonlinear_factor is not None:
            return self._write_f06_transient_block(words, header, pageStamp, page_num, f)
        return self._write_f06_block(words, header, pageStamp, page_num, f)

    def _write_table_header(self, f):
        header = [
            # table 1
            4, 0, 4,
            4, -1, 4,
            4, 0, 4,  # 9i

            4, 2, 4,  # 3i
            8, self.table_name, 8,  # 1i, 8s, 1i
            4, 2, 4,  #3i

            #===============
            # table 2
            4, 0, 4,
            4, -2, 4,
            4, 0, 4,  # 9i

            4, 7, 4,
            28,  # 4i -> 13i
            'OUG1    ', 3, 6, 14, 0, 1,   # subtable,todays date 3/6/2014, 0, 1
            28,
            ]
        #header_format = '13i 4s 4i' + '13i 8s 6i'
        header_format = '13i4s4i' + '13i8s6i'
        f.write(Struct(fascii, 'OUG header', header_format).pack(header))

    def write_op2(self, f, fascii, is_mag_phase=False):
        if 0:
            #recordi = ['OUG1    ', 2, 26, 14, 0, 1]
            #[subtable_name, month=2, day=26, year=2014, zero=0, one=1]
            subtable_name = 'OUG1    '
            month = 1
            day = 30
            year = 2014
            recordi = [subtable_name, month, day, year, 0, 1]

            approach_code = 1
            table_code = 1
            isubcase = 1
            lsdvmn = 1
            random_code = 0
            format_code = 1
            num_wide = 8
            acoustic_flag = 0
            thermal = 0
            Title = ' ' * 128
            subtitle = ' ' * 128
            label = ' ' * 128
            ftable3 = '24i 128s 128s 128s'
            oCode = 0
            table3 = [
                aCode, tCode, 0, isubcase, lsdvmn,
                0, 0, random_code, format_code, num_wide,
                oCode, acoustic_flag, 0, 0, 0,
                0, 0, 0, 0, 0,
                0, thermal, thermal, 0, Title,
                subtitle, label, ]
            #table3 = [approach_code, table_code, 0, isubcase, lsdvmn,
                         #0, 0, random_code, format_code, num_wide,
                         #0, acoustic_flag, acoustic_flag, 0, 0,
                         #0, 0, 0, 0, 0,
                         #0, thermal, thermal, 0,
                         #Title, subtitle, label,
            #]
            pack(ftable3, *table3)
            fdata = '7i 8s 8i 8s 6i'

            fmt_start = '6i'
            fmt_end = '3i'
            table_num = 3
            start = [
                4, 146, 4,
                4, 584, 4,
                ]
            end = [4, 584, 4]

            data = start + [
                4,0,4,
                4,2,4,
                8, 'OUGV1   ', 8,

                4, 0, 4,
                4, 7, 4,
                28] + recordi + [28,
            ] + table3 + end

            out = pack(fascii, 'table3', fmt_start + fdata + ftable3 + fmt_end, *data)
            f.write(out)

        if self.nonlinear_factor is not None:
            return self._write_op2_transient_block(f, fascii)
        return self._write_op2_block(f, fascii)

    def _write_op2_block(self, f, fascii):
        from six import iteritems
        nnodes = len(self.translations)
        nwords = nnodes * self.num_wide
        nbytes = nwords * 4
        print('nnodes=%s nbytes=%s' % (nnodes, nbytes))
        #nbytes = 16386
        #assert nbytes >= 16386, nbytes
        print('nbytes=%s' % nbytes)
        header = [
            4, 0, 4,
            4, nwords, 4,
            nbytes #recordi, nbytes*4,
        ]
        f.write(pack(fascii, 'header4', '7i', header))

        data = []
        fmt = '2i 6f'
        device_code = self.device_code
        for nodeID, translation in sorted(iteritems(self.translations)):
            rotation = self.rotations[nodeID]
            grid_type = self.gridTypes[nodeID]
            grid_type = 1
            #print('grid_type=%s' % grid_type)

            (dx, dy, dz) = translation
            (rx, ry, rz) = rotation
            vals = [dx, dy, dz, rx, ry, rz]
            data = [device_code + 10*nodeID, grid_type, dx, dy, dz, rx, ry, rz]
            f.write(pack(fascii, 'nid, gridType, dx, dy, dz, rx, ry, rz', fmt, data))
        f.write(pack(fascii, 'closer', 'i', nbytes))

class ComplexDisplacement(ComplexTableObject):  # approach_code=1, sort_code=0, thermal=0
    def __init__(self, data_code, is_sort1, isubcase, dt):
        ComplexTableObject.__init__(self, data_code, is_sort1, isubcase, dt)

    def write_f06(self, header, pageStamp, page_num=1, f=None, is_mag_phase=False):
        words = ['                                       C O M P L E X   D I S P L A C E M E N T   V E C T O R\n']
        return self._write_f06_transient_block(words, header, pageStamp, page_num, f, is_mag_phase)
