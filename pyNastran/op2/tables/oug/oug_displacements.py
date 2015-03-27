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

def write_markers(f, fascii, msg, markers):
    s = sStruct('3i')
    for marker in markers:
        fascii.write('%s=%s ' % (msg, [4, marker, 4]))
        data = [4, marker, 4]
        f.write(s.pack(*data))

from pyNastran.op2.resultObjects.tableObject import RealTableArray, ComplexTableArray, RealTableObject, ComplexTableObject

def write_table_header(f, fascii, table_name):
    table0 = [
        4, 2, 4,
        8, table_name, 8,
    ]
    table0_format = '4i 4s i'
    f.write(pack(fascii, 'OUG header0', table0_format, table0))


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

    def write_f06(self, header, page_stamp, page_num=1, f=None, is_mag_phase=False):
        words = ['                                             D I S P L A C E M E N T   V E C T O R\n', ]
                 #' \n',
                 #'      POINT ID.   TYPE          T1             T2             T3             R1             R2             R3\n']
        #words += self.get_table_marker()
        if self.nonlinear_factor is not None:
            return self._write_f06_transient_block(words, header, page_stamp, page_num, f)
        return self._write_f06_block(words, header, page_stamp, page_num, f)


class ComplexDisplacementArray(ComplexTableArray):
    def __init__(self, data_code, is_sort1, isubcase, dt):
        ComplexTableArray.__init__(self, data_code, is_sort1, isubcase, dt)

    def write_f06(self, header, page_stamp, page_num=1, f=None, is_mag_phase=False):
        words = ['                                       C O M P L E X   D I S P L A C E M E N T   V E C T O R\n']
        return self._write_f06_transient_block(words, header, page_stamp, page_num, f, is_mag_phase)


class RealDisplacement(RealTableObject):  # approach_code=1, thermal=0
    def __init__(self, data_code, is_sort1, isubcase, dt):
        RealTableObject.__init__(self, data_code, is_sort1, isubcase, dt)

    def write_f06(self, header, page_stamp, page_num=1, f=None, is_mag_phase=False):
        words = ['                                             D I S P L A C E M E N T   V E C T O R\n',
                 ' \n',
                 '      POINT ID.   TYPE          T1             T2             T3             R1             R2             R3\n']
        words += self.get_table_marker()

        if self.nonlinear_factor is not None:
            return self._write_f06_transient_block(words, header, page_stamp, page_num, f)
        return self._write_f06_block(words, header, page_stamp, page_num, f)

    def _write_table_header(self, f, fascii):
        #get_nmarkers- [4, 0, 4]
        #marker = [4, 2, 4]
        #table_header = [8, 'BOUGV1  ', 8]
        write_table_header(f, fascii, self.table_name)


        #read_markers -> [4, -1, 4]
        #get_nmarkers- [4, 0, 4]
        #read_record - marker = [4, 7, 4]
        #read_record - record = [28, recordi, 28]

        write_markers(f, fascii, 'OUG header1a', [-1, 0, 7])

        table1_fmt = '9i'
        table1 = [
            28,
            1, 2, 3, 4, 5, 6, 7,
            28,
        ]
        f.write(pack(fascii, 'OUG header1b', table1_fmt, table1))

        import datetime
        today = datetime.datetime.today()
        month = today.month
        day = today.day
        year = today.year - 2000
        #recordi = [subtable_name, month, day, year, 0, 1]

        write_markers(f, fascii, 'OUG header2a', [0, -2, 0, 7])
        table2 = [
            28,  # 4i -> 13i
            'OUG1    ', month, day, year, 0, 1,   # subtable,todays date 3/6/2014, 0, 1
            28,
            ]
        table2_format = 'i8s6i'
        f.write(pack(fascii, 'OUG header2b', table2_format, table2))

    def write_op2(self, f, fascii, is_mag_phase=False):
        self._write_table_header(f, fascii)
        #recordi = ['OUG1    ', 2, 26, 14, 0, 1]
        #[subtable_name, month=2, day=26, year=2014, zero=0, one=1]
        subtable_name = 'OUG1    '

        aCode = 1
        tCode = 1
        approach_code = 0 #self.approach_code
        table_code = self.table_code
        isubcase = self.isubcase
        lsdvmn = 1
        random_code = 0
        format_code = 1
        num_wide = self.num_wide
        acoustic_flag = 0
        thermal = 0
        Title = '%-128s' % self.Title
        subtitle = '%-128s' % self.subtitle
        label = '%-128s' % self.label
        ftable3 = '24i 128s 128s 128s'
        oCode = 0
        table3 = [
            aCode, tCode, 0, isubcase, lsdvmn,
            0, 0, random_code, format_code, num_wide,
            oCode, acoustic_flag, 0, 0, 0,
            0, 0, 0, 0, 0,
            0, thermal, thermal, 0, Title,
            subtitle, label, ]

        write_markers(f, fascii, 'OUG header3a', [-3, 1, 0])
        write_markers(f, fascii, 'OUG header3b', [146])

        data = [584] + table3 + [584]
        fmt = 'i' + ftable3 + 'i'
        f.write(pack(fascii, 'OUG header 3c', fmt, data))

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
            if grid_type == 'G':
                grid_type = 1
            else:
                raise RuntimeError(gridType)

            (dx, dy, dz) = translation
            (rx, ry, rz) = rotation
            vals = [dx, dy, dz, rx, ry, rz]
            data = [device_code + 10*nodeID, grid_type, dx, dy, dz, rx, ry, rz]
            f.write(pack(fascii, 'nid, gridType, dx, dy, dz, rx, ry, rz', fmt, data))
        f.write(pack(fascii, 'closer', 'i', nbytes))

class ComplexDisplacement(ComplexTableObject):  # approach_code=1, sort_code=0, thermal=0
    def __init__(self, data_code, is_sort1, isubcase, dt):
        ComplexTableObject.__init__(self, data_code, is_sort1, isubcase, dt)

    def write_f06(self, header, page_stamp, page_num=1, f=None, is_mag_phase=False):
        words = ['                                       C O M P L E X   D I S P L A C E M E N T   V E C T O R\n']
        return self._write_f06_transient_block(words, header, page_stamp, page_num, f, is_mag_phase)
