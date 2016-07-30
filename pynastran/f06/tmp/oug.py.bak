#pylint: disable=C0301
from __future__ import print_function
from six.moves import zip

from pyNastran.op2.tables.oug.oug_displacements import RealDisplacement, ComplexDisplacement
from pyNastran.op2.tables.oug.oug_eigenvectors import (Eigenvector, ComplexEigenvector,
                                                       RealEigenvectorArray, ComplexEigenvectorArray)
from pyNastran.op2.tables.oug.oug_temperatures import RealTemperature


class OUG(object):
    def _read_f06_table(self, data_types, debug=False):
        pass
    def __init__(self):
        self.iSubcases = []
        self.i = 0

    def _real_eigenvectors(self, marker):
        """
        Reads real eigenvector table accounting for blank entries

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
        cycle, imode = marker.strip().split('R E A L   E I G E N V E C T O R   N O .')
        imode = int(imode)

        cycles = cycle.strip().split('=')
        assert 'CYCLES' == cycles[0].strip(), 'marker=%s' % marker
        cycle = float(cycles[1])

        (subcase_name, isubcase, transient, dt, analysis_code, is_sort1) = self._read_f06_subcase_header()
        eigenvalue_real = transient[1]
        headers = self.skip(2)

        data_code = {
            'log': self.log, 'analysis_code': analysis_code,
            'device_code': 1, 'table_code': 7, 'sort_code': 0,
            'sort_bits': [0, 0, 0], 'num_wide': 8, 'format_code': 1,
            'mode': imode, 'eigr': eigenvalue_real, 'mode_cycle': cycle,
            'data_names': ['mode', 'eigr', 'mode_cycle'],
            'name': 'mode', 'table_name': 'OUGV1',
            'nonlinear_factor': imode,
        }
        data = self._real_f06_real_table_data(allow_blanks=True)

        #print("cycle=%-8s eigen=%s" % (cycle, eigenvalue_real))
        if isubcase in self.eigenvectors:
            self.eigenvectors[isubcase].read_f06_data(data_code, data)
        else:
            is_sort1 = True
            if self.is_vectorized:
                vector = RealEigenvectorArray(data_code, is_sort1, isubcase, imode, f06_flag=True)
            else:
                vector = Eigenvector(data_code, is_sort1, isubcase, imode)
            vector.read_f06_data(data_code, data)
            self.eigenvectors[isubcase] = vector

    def _complex_eigenvectors(self, marker):
        """
        Reads real eigenvector table accounting for blank entries

        ::
                COMPLEX EIGENVALUE =  0.000000E+00,  3.805272E+02
                                                 C O M P L E X   E I G E N V E C T O R   NO.          5
                                                                    (REAL/IMAGINARY)

                POINT ID.   TYPE          T1             T2             T3             R1             R2             R3
          0            1      G      0.0            0.0            0.0           -3.041936E-02  -1.922321E-01   0.0
                                     0.0            0.0            0.0            0.0            0.0            0.0
          0          227      S      1.418276E-03  -2.095675E-04  -1.663478E-03   2.633889E-03   4.171373E-14  -2.633889E-03
                                     0.0            0.0            0.0            0.0            0.0            0.0
          0          233      S      1.663478E-03   2.095675E-04  -1.418276E-03
                                     0.0            0.0            0.0

        * analysis_code = 2 (Normal modes)
        * table_code    = 7 (Eigenvector)
        * device_code   = 1 (Print)
        * sort_code     = ? (Sort2,Real/Imag or Mag/Phase,Sorted Results) => sort_bits = [0,0,0]
        * format_code   = 1 (???)
        * num_wide      = 16 (8*2)
        """

        #self._read_table_dummy()
        #return

        blank_imode = marker.strip().split('C O M P L E X   E I G E N V E C T O R   NO.')
        blank, imode = blank_imode
        assert blank == '', blank_imode
        imode = int(imode)

        (subcase_name, isubcase, transient, _dt, analysis_code, _is_sort1) = self._read_f06_subcase_header()
        eigenvalue_real, eigenvalue_imag = transient[1]

        real_imag_mag_phase = self._get_next_line().strip()
        headers = self.skip(2).split()

        if '(REAL/IMAGINARY)' in real_imag_mag_phase:
            is_mag_phase = False
        else:
            raise RuntimeError(real_imag_mag_phase)

        if 'POINT' in headers and 'ID.' in headers:
            is_sort1 = True
        else:
            msg = 'Only SORT1 is supported; headers=%s' % headers
            raise NotImplementedError(msg)

        data_code = {
            'log': self.log, 'analysis_code': analysis_code,
            'device_code': 1, 'table_code': 7, 'sort_code': 0,
            'sort_bits': [0, 0, 0], 'num_wide': 16, 'format_code': 1,
            'mode': imode, 'eigr': eigenvalue_real, 'eigi' : eigenvalue_imag,
            'data_names': ['mode', 'eigr', 'eigi'],
            'name': 'mode', 'table_name': 'OUGV1',
            'nonlinear_factor': imode,
        }
        data = self._real_f06_complex_table_data(allow_blanks=True, is_mag_phase=is_mag_phase)

        if isubcase in self.eigenvectors:
            self.eigenvectors[isubcase].read_f06_data(data_code, data)
        else:
            is_sort1 = True
            if self.is_vectorized:
                vector = ComplexEigenvectorArray(data_code, is_sort1, isubcase, imode, f06_flag=True)  ## TODO: f06_flag isn't defined...
            else:
                vector = ComplexEigenvector(data_code, is_sort1, isubcase, imode, apply_data_code=False)
            vector.read_f06_data(data_code, data)
            self.eigenvectors[isubcase] = vector

    def _real_f06_real_table_data(self, allow_blanks=False):
        """
        Reads real displacement/velocity/spc forces/mpc forces
        Handles GRIDs and SPOINTs.

        Parameters
        ----------
        allow_blanks : bool; default=False
            Accounting for blank entries (e.g. on eigenvector)

        Returns
        -------
        data : List[entries]
            the parsed data

        .. todo:: support L, H, and R points
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

            fields = [line[24:39], line[39:54], line[54:69], line[69:84], line[84:99], line[99:114]]
            if grid_type in ['G', 'L']:
                # .. todo:: are L points single DOFs?
                # we do know they're greater than the max value...
                sline = [node_id, grid_type]
                if allow_blanks:
                    sline += [float(val) if val.strip() != '' else 0.0 for val in fields]
                else:
                    sline += [float(val) for val in fields]
                data.append(sline)
            elif grid_type == 'S':
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

    def _real_f06_complex_table_data(self, allow_blanks=False, is_mag_phase=False):
        """
        Reads complex displacement/velocity/spc forces/mpc forces
        Handles GRIDs and SPOINTs.

        Parameters
        ----------
        allow_blanks : bool; default=False
            Accounting for blank entries (e.g. on eigenvector)

        Returns
        -------
        data : list[entries]
            the parsed data

        .. todo:: support L, H, and R points
        """
        field_length = 15  # width of each eigenvector field
        num_fields = 6     # the number of fields (T1, T2, T3, R1, R2, R3)
        data = []

        n = 0
        while 1:
            line1 = self.infile.readline()[1:].rstrip()

            # TODO: add catch for FATAL
            if 'PAGE' in line1:
                break
            line2 = self.infile.readline()[1:].rstrip()

            #: point ID (int)
            node_id = int(line1[:14].strip())

            #: TYPE (str)
            grid_type = line1[14:24].strip()

            fields = [
                line1[24:39], line2[24:39],
                line1[39:54], line2[39:54],
                line1[54:69], line2[54:69],
                line1[69:84], line2[69:84],
                line1[84:99], line2[84:99],
                line1[99:114],line2[99:114],
            ]
            if grid_type in ['G', 'L']:
                # .. todo:: are L points single DOFs?
                # we do know they're greater than the max value...
                sline = [node_id, grid_type]
                fields2 = [float(val) if val.strip() != '' else None for val in fields]

                for i in range(0, len(fields), 2):
                    value = complex(float(fields[i]), float(fields[i+1]))
                    sline.append(value)
                assert len(sline) == 8, sline
                data.append(sline)
            elif grid_type == 'S':
                ioffset = 0
                for i in range(0, len(fields), 2):
                    if fields[i].strip() == '':
                        # incomplete spoint line
                        # ID, S, 1.0, 2.0, 3.0, 4.0, 5.0, None
                        break
                    value = complex(float(fields[i]), float(fields[i+1]))
                    sline = [node_id + ioffset,
                             grid_type, value, 0.0, 0.0, 0.0, 0.0, 0.0]
                    data.append(sline)
                    #print(sline)
                    ioffset += 1
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
        (subcase_name, isubcase, transient, dt, analysis_code, is_sort1) = self._read_f06_subcase_header()
        headers = self.skip(2)
        #raise RuntimeError(headers)

        data_code = {
            'analysis_code': analysis_code,
            'device_code': 1, 'table_code': 1,
            'sort_code': 0, 'sort_bits': [0, 0, 0], 'num_wide': 8,
            'table_name': 'OUG', 'nonlinear_factor': dt,
            'lsdvmn': 1, 'format_code': 3,
            'data_names':['lsdvmn']
        }
        #print("headers = %s" % (headers))
        #print("transient =", transient)

        data = self._real_f06_real_table_data(allow_blanks=False)

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
        (subcase_name, isubcase, transient, dt, analysis_code, is_sort1) = self._read_f06_subcase_header()
        #print("transient =", transient)
        #print("dt =", dt)
        name = transient[0]
        data_names = [name]
        headers = self.skip(3)
        data = []

        data_code = {
            'analysis_code': 5, 'device_code': 1,
            'table_code': 1, 'sort_code': 2, 'sort_bits': [0, 1, 1],
            'num_wide': 14, 'format_code': 3, 'table_name': 'OUGV1',
            'nonlinear_factor': dt,
            #'mode':iMode,'eigr':transient[1], 'mode_cycle':cycle,
            'data_names': data_names,
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
        (subcase_name, isubcase, transient, dt, analysis_code, is_sort1) = self._read_f06_subcase_header()

        headers = self.skip(2)
        #print "headers = %s" % headers
        data = self._read_temperature_table()

        data_code = {
            'log': self.log, 'analysis_code': 1, 'device_code': 1,
            'table_code': 1, 'sort_code': 0, 'sort_bits': [0, 0, 0],
            'num_wide': 8, 'table_name': 'OUG', 'nonlinear_factor': dt,
            'data_names':['lsdvmn'],
            'thermal':1,
            'format_code':1,
            #'element_name':eType,'s_code':0,'stress_bits':stress_bits
        }

        if transient:
            name = transient[0]
            data_code['data_names'] = [name + 's']
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

    def _parse_line_temperature(self, sline, formats):
        out = []
        for (entry, iformat) in zip(sline, formats):
            if entry is '':
                return out
            #print "sline=%r\n entry=%r iformat=%r" % (sline, entry, iformat)
            entry2 = iformat(entry)
            out.append(entry2)
        return out
