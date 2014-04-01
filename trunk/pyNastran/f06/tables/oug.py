#pylint: disable=C0301,C0103,W0612
from itertools import izip

from pyNastran.op2.tables.oug.oug_displacements import RealDisplacement, ComplexDisplacement
from pyNastran.op2.tables.oug.oug_eigenvectors import EigenvectorObject  # ,ComplexEigenVectorObject
from pyNastran.op2.tables.oug.oug_temperatures import RealTemperature


class OUG(object):
    def _read_f06_table(self, data_types, debug=False):
        pass
    def __init__(self):
        self.displacements = {}
        self.temperatures = {}
        self.eigenvectors = {}
        self.accelerations = {}
        self.velocities = {}
        self.iSubcases = []
        self.i = 0

    def _real_eigenvectors(self, marker):
        """
        Reads real eigenvector table accounting for blank entries

        :param self:   the object pointer

        ::
                                                                                                                 SUBCASE 1
          EIGENVALUE =  6.158494E+07
              CYCLES =  1.248985E+03         R E A L   E I G E N V E C T O R   N O .          1

          POINT ID.   TYPE          T1             T2             T3             R1             R2             R3
                 1      G      2.547245E-17  -6.388945E-16   2.292728E+00  -1.076928E-15   2.579163E-17   0.0
              2002      G     -6.382321E-17  -1.556607E-15   3.242408E+00  -6.530917E-16   1.747180E-17   0.0
              2003      G     -6.382321E-17  -1.556607E-15   3.242408E+00
              2004      S     -6.382321E-17  -1.556607E-15   3.242408E+00

        * analysis_code = 2 (Normal modes)
        * table_code    = 7 (Eigenvector)
        * device_code   = 1 (Print)
        * sort_code     = 0 (Sort2,Real,Sorted Results) => sort_bits = [0,0,0]
        * format_code   = 1 (Real)
        * #s_code        = 0 (Stress)
        * num_wide      = 8 (???)
        """
        cycle, iMode = marker.strip().split('R E A L   E I G E N V E C T O R   N O .')
        iMode = int(iMode)

        cycles = cycle.strip().split('=')
        #print smarker
        assert 'CYCLES' == cycles[0].strip(), 'marker=%s' % marker
        cycle = float(cycles[1])

        #print "marker = |%s|" %(marker)
        #subcaseName = '???'
        #print self.storedLines
        #isubcase = self.storedLines[-2].strip()[1:]
        #isubcase = int(isubcase.strip('SUBCASE '))
        #print "subcaseName=%s isubcase=%s" %(subcaseName,isubcase)
        (subcaseName, isubcase, transient, dt, analysis_code, is_sort1) = self._read_f06_subcase_header()
        eigenvalue_real = transient[1]
        headers = self.skip(2)

        data_code = {
            'log': self.log, 'analysis_code': analysis_code,
            'device_code': 1, 'table_code': 7, 'sort_code': 0,
            'sort_bits': [0, 0, 0], 'num_wide': 8, 'format_code': 1,
            'mode': iMode, 'eigr': eigenvalue_real, 'mode_cycle': cycle,
            'dataNames': ['mode', 'eigr', 'mode_cycle'],
            'name': 'mode', 'table_name': 'OUGV1',
            'nonlinear_factor': iMode,
            #'s_code':0,
            #'element_name':'CBAR','element_type':34,'stress_bits':stress_bits,
        }
        data = self._real_f06_table_data(allow_blanks=True)

        #print("cycle=%-8s eigen=%s" % (cycle, eigenvalue_real))
        #print "isubcase = %s" % isubcase
        if isubcase in self.eigenvectors:
            self.eigenvectors[isubcase].read_f06_data(data_code, data)
        else:
            is_sort1 = True
            self.eigenvectors[isubcase] = EigenvectorObject(data_code, is_sort1, isubcase, iMode)
            self.eigenvectors[isubcase].read_f06_data(data_code, data)

    def _real_f06_table_data(self, allow_blanks=False):
        """
        Reads real displacement/velocity/spc forces/mpc forces
        Handles GRIDs and SPOINTs.

        :param self:         the object pointer
        :param allow_blanks: Accounting for blank entries (e.g. on eigenvector)
                             default=False
        :returns data:       the parsed data

        ..todo:: support L, H, and R points
        """
        field_length = 15  # width of each eigenvector field
        num_fields = 6     # the number of fields (T1, T2, T3, R1, R2, R3)
        data = []

        n = 0
        while 1:
            line = self.infile.readline()[1:].rstrip()

            # TODO: add catch for FATAL
            if 'PAGE' in line:
                break

            #: point ID (int)
            node_id = int(line[:14].strip())

            #: TYPE (str)
            grid_type = line[14:24].strip()

            if grid_type == 'G':
                sline = [node_id, grid_type]
                fields = [line[24:39], line[39:54], line[54:69], line[69:84], line[84:99], line[99:114]]
                if allow_blanks:
                    sline += [float(val) if val.strip() != '' else 0.0 for val in fields]
                else:
                    sline += [float(val) for val in fields]
                data.append(sline)
            elif grid_type == 'S':
                fields = [line[24:39], line[39:54], line[54:69], line[69:84], line[84:99],  line[99:114]]
                fields = [float(val) if val.strip() != '' else None for val in fields]
                for ioffset, field in enumerate(fields):
                    if field is None:
                        # incomplete spoint line
                        # ID, S, 1.0, 2.0, 3.0, 4.0, 5.0, None
                        break
                    sline = [node_id + ioffset,
                             grid_type, field, 0.0, 0.0, 0.0, 0.0, 0.0]
                    data.append(sline)
            else:
                raise NotImplementedError('grid_type = %r' % grid_type)
            self.i += 1
            n += 1
        return data

    def _displacement_vector(self):
        """
        ::
                                               D I S P L A C E M E N T   V E C T O R

          POINT ID.   TYPE          T1             T2             T3             R1             R2             R3
                 1      G      9.663032E-05   0.0           -2.199001E-04   0.0           -9.121119E-05   0.0
                 2      G      0.0            0.0            0.0            0.0            0.0            0.0
                 3      G      0.0            0.0            0.0            0.0            0.0            0.0

        * analysis_code = 1 (Statics)
        * device_code   = 1 (Print)
        * table_code    = 1 (Displacement)
        * sort_code     = 0 (Sort2,Real,Sorted Results) => sort_bits = [0,0,0]
        * num_wide      = 8 (???)
        """
        (subcaseName, isubcase, transient, dt, analysis_code, is_sort1) = self._read_f06_subcase_header()
        #print "subcaseName=%r isubcase=%s"  % (subcaseName, isubcase)
        headers = self.skip(2)
        #raise RuntimeError(headers)

        data_code = {'log': self.log, 'analysis_code': analysis_code,
                    'device_code': 1, 'table_code': 1,
                    'sort_code': 0, 'sort_bits': [0, 0, 0], 'num_wide': 8,
                    'table_name': 'OUG', 'nonlinear_factor': dt,
                    'lsdvmn': 1, 'format_code': 3,
                    'dataNames':['lsdvmn']}
        #print "headers = %s" %(headers)
        #print "transient =", transient

        data = self._real_f06_table_data(allow_blanks=False)

        if isubcase in self.displacements:
            self.displacements[isubcase].add_f06_data(data, transient)
        else:
            is_sort1 = True
            disp = RealDisplacement(data_code, is_sort1, isubcase, dt)
            disp.add_f06_data(data, transient)
            self.displacements[isubcase] = disp
        self.iSubcases.append(isubcase)

    def _complex_displacement_vector(self):
        """
        ::

            BACKWARD WHIRL
                                                                                                                   SUBCASE 2
            POINT-ID =       101
                                             C O M P L E X   D I S P L A C E M E N T   V E C T O R
                                                               (MAGNITUDE/PHASE)

            FREQUENCY   TYPE          T1             T2             T3             R1             R2             R3
          2.000000E+01     G      3.242295E-16   1.630439E-01   1.630439E-01   1.691497E-17   1.362718E-01   1.362718E-01
                                  196.0668        90.0000       180.0000        63.4349       180.0000       270.0000

        * table_code    = 1 (Displacement)
        * format_code   = 3 (Magnitude/Phase)
        * sort_bits     = [0,1,1]  (Sort1,Real/Imaginary,RandomResponse)
        * analysis_code = 5 (Frequency)
        * sort_code     = 2 (Random Response)
        """
        (subcaseName, isubcase, transient, dt, analysis_code, is_sort1) = self._read_f06_subcase_header()
        #print("transient =", transient)
        #print("dt =", dt)
        name = transient[0]
        data_names = [name]
        headers = self.skip(3)
        data = []

        data_code = {'log': self.log, 'analysis_code': 5, 'device_code': 1,
                    'table_code': 1, 'sort_code': 2, 'sort_bits': [0, 1, 1],
                    'num_wide': 14, 'format_code': 3, 'table_name': 'OUGV1',
                    'nonlinear_factor': dt,
                    #'mode':iMode,'eigr':transient[1], 'mode_cycle':cycle,
                    'dataNames': data_names,
                    'name': name,
                    #'s_code':0,
                    }
        while 1:
            line = self.infile.readline()[1:].rstrip('\r\n ')
            if 'PAGE' in line:
                break
            sline = line.strip().split()
            line2 = self.infile.readline()[1:].rstrip('\r\n ')
            sline += line2.strip().split()
            #print("sline = ", sline)

            freq = float(sline[0])
            grid_type = sline[1].strip()
            #dxyz_ryz = sline[2:]
            if grid_type == 'G':
                dx = float(sline[2]) + float(sline[8])*1j
                dy = float(sline[3]) + float(sline[9])*1j
                dz = float(sline[4]) + float(sline[10])*1j
                rx = float(sline[5]) + float(sline[11])*1j
                ry = float(sline[6]) + float(sline[12])*1j
                rz = float(sline[7]) + float(sline[13])*1j
                out = [freq, grid_type, dx, dy, dz, rx, ry, rz]
                data.append(out)
            elif grid_type == 'S':
                out = line[2:]
            else:
                raise NotImplementedError(grid_type)
            self.i += 2
        self.i += 1

        if isubcase in self.displacements:
            self.displacements[isubcase].add_f06_data(data, transient)
        else:
            is_sort1 = True
            disp = ComplexDisplacement(data_code, is_sort1, isubcase, dt)
            disp.add_f06_data(data, transient)
            self.displacements[isubcase] = disp
        self.iSubcases.append(isubcase)

    def _temperature_vector(self):
        """
        ::

          LOAD STEP =  1.00000E+00
                                                T E M P E R A T U R E   V E C T O R

          POINT ID.   TYPE      ID   VALUE     ID+1 VALUE     ID+2 VALUE     ID+3 VALUE     ID+4 VALUE     ID+5 VALUE
                 1      S      1.300000E+03   1.300000E+03   1.300000E+03   1.300000E+03   1.300000E+03   1.300000E+03
                 7      S      1.300000E+03   1.300000E+03   1.300000E+03   1.300000E+03
          analysis_code = 1 (Statics)
          device_code   = 1 (Print)
          table_code    = 1 (Displacement/Temperature)
          sort_code     = 0 (Sort2,Real,Sorted Results) => sort_bits = [0,0,0]
          format_code   = 1 (Real)
          s_code        = 0 (Stress)
          num_wide      = 8 (???)
        """
        (subcaseName, isubcase, transient, dt, analysis_code, is_sort1) = self._read_f06_subcase_header()
        #print transient

        headers = self.skip(2)
        #print "headers = %s" % headers
        data = self._read_temperature_table()

        data_code = {
            'log': self.log, 'analysis_code': 1, 'device_code': 1,
            'table_code': 1, 'sort_code': 0, 'sort_bits': [0, 0, 0],
            'num_wide': 8, 'table_name': 'OUG', 'nonlinear_factor': dt,
            'dataNames':['lsdvmn'],
            'thermal':1,
            'format_code':1,
            #'element_name':eType,'s_code':0,'stress_bits':stress_bits
        }

        if transient:
            name = transient[0]
            data_code['dataNames'] = [name + 's']
            data_code['name'] = name

        if isubcase in self.temperatures:
            self.temperatures[isubcase].add_f06_data(data, transient)
        else:
            is_sort1 = True
            temp = RealTemperature(data_code, is_sort1, isubcase, dt)
            temp.add_f06_data(data, transient)
            self.temperatures[isubcase] = temp
        self.iSubcases.append(isubcase)

    def _read_temperature_table(self):
        data = []
        Format = [int, str, float, float, float, float, float, float]
        while 1:
            line = self.infile.readline()[1:].rstrip('\r\n ')
            if 'PAGE' in line:
                return data
            sline = [line[0:15], line[15:22].strip(), line[22:40], line[40:55], line[55:70], line[70:85], line[85:100], line[100:115]]
            sline = self._parse_line_temperature(sline, Format)
            data.append(sline)
        return data

    def _parse_line_temperature(self, sline, Format):
        out = []
        for (entry, iFormat) in izip(sline, Format):
            if entry is '':
                return out
            #print "sline=|%r|\n entry=|%r| format=%r" %(sline, entry, iFormat)
            entry2 = iFormat(entry)
            out.append(entry2)
        return out
