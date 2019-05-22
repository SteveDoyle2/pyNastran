class MONPNT1:
    """MONPNT1 table"""
    def __init__(self, frequencies, matrices, comp_matrices):
        self.frequencies = frequencies

        # [inertial, external, flexible_increment, gust, total_aero, total]
        self.comp_matrices = comp_matrices
        self.cx_matrix = matrices[comp_matrices[0]]
        self.cy_matrix = matrices[comp_matrices[1]]
        self.cz_matrix = matrices[comp_matrices[2]]
        self.rx_matrix = matrices[comp_matrices[3]]
        self.ry_matrix = matrices[comp_matrices[4]]
        self.rz_matrix = matrices[comp_matrices[5]]

    def write(self, f06_file, page_stamp='', page_num=1):
        comps = ['CX', 'CY', 'CZ', 'CMX', 'CMY', 'CMZ']
        dok1 = self.cx_matrix.data.todok()
        dok2 = self.cy_matrix.data.todok()
        dok3 = self.cz_matrix.data.todok()
        dok4 = self.rx_matrix.data.todok()
        dok5 = self.ry_matrix.data.todok()
        dok6 = self.rz_matrix.data.todok()

        for icomp, comp in enumerate(comps):
            lines = [
                '                              S T R U C T U R A L   M O N I T O R   P O I N T   I N T E G R A T E D   L O A D S (MONPNT1)\n'
                '                                                       (REAL/IMAGINARY)\n'
                '\n'
                '        MONITOR POINT NAME = AEROSG2D          COMPONENT =  %3s          GENERAL            SUBCASE NO.        1\n' % comp,
                '        LABEL = Full Vehicle Integrated Loads                           \n'
                '          CP =        2          X = 0.000000E+00          Y = 0.000000E+00          Z = 0.000000E+00          CD =       2\n'
                '\n'
                '          FREQUENCY      INERTIAL       EXTERNAL       FLEXIBLE       GUST           TOTAL          TOTAL       \n'
                '                                                       INCREMENT                     AERO                       \n'
                '          ------------   ------------   ------------   ------------   ------------   ------------   ------------\n'
                #'          0.000000E+00   0.000000E+00  -1.000000E+00   0.000000E+00   0.000000E+00   0.000000E+00  -1.000000E+00\n'
                #'                         0.000000E+00   0.000000E+00   0.000000E+00   0.000000E+00   0.000000E+00   0.000000E+00\n'
            ]
            f06_file.write(''.join(lines))
            for ifreq, freq in enumerate(self.frequencies):
                #if icomp == 4:
                v1 = dok1[(icomp, ifreq)]
                v2 = dok2[(icomp, ifreq)]
                v3 = dok3[(icomp, ifreq)]

                v4 = dok4[(icomp, ifreq)]
                v5 = 0.
                v6 = 0.
                #v5 = dok5[(icomp, ifreq)]
                #v6 = dok6[(icomp, ifreq)]

                #v4 = v3
                #v6 = v1 + v2 + v3 + v4 # wrong
                #v6 = v2

                #v2 = v6
                #v4 = v3

                f06_file.write('         %13E  %13E  %13E  %13E  %13E  %13E  %13E\n' % (
                    freq, v1.real, v2.real, v3.real, v4.real, v5.real, v6.real))
                f06_file.write('                        %13E  %13E  %13E  %13E  %13E  %13E\n' % (
                    v1.imag, v2.imag, v3.imag, v4.imag, v5.imag, v6.imag))
                page_num += 1
            f06_file.write(page_stamp % page_num)
        page_num -= 1
        return page_num


class MONPNT3:
    """MONPNT3 table"""
    def __init__(self, frequencies, matrix):
        self.name = matrix.name
        self.frequencies = frequencies
        self.data = matrix.data

    def write(self, f06_file, page_stamp='', page_num=1):
        comps = ['CX', 'CY', 'CZ', 'CMX', 'CMY', 'CMZ']
        matrix = self.data
        dok = matrix.todok()
        assert self.frequencies is not None, self.frequencies

        for icomp, comp in enumerate(comps):
            lines = [
                '                S T R U C T U R A L   I N T E G R A T E D   F R E E   B O D Y   M O N I T O R   P O I N T   L O A D S (MONPNT3)\n'
                '                                                       (REAL/IMAGINARY)\n'
                '\n'
                '        MONITOR POINT NAME = M3_B              COMPONENT =  %3s                                     SUBCASE NO.        1\n' % comp,
                '        LABEL = NODE 2 AND ELEMENT 102                                  \n'
                '          CP =        0          X = 1.000000E+00          Y = 0.000000E+00          Z = 0.000000E+00          CD =       0\n'
                '\n'
                '          FREQUENCY      RESULTANT   \n'
                '                                     \n'
                '          ------------   ------------\n'
            ]
            f06_file.write(''.join(lines))

            for ifreq, freq in enumerate(self.frequencies):
                #print('icomp=%r ifreq=%r' % (icomp, ifreq))
                val = dok[(icomp, ifreq)]
                f06_file.write('          %13E %13E\n                        %13E\n' % (
                    freq, val.real, val.imag))
            page_num += 1
            f06_file.write(page_stamp % page_num)
        page_num -= 1
        return page_num

    def __repr__(self):
        msg = 'MONPNT3; nfreqs=%s\n' % len(self.frequencies)
        return msg
