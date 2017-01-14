from six.moves import zip
#from pyNastran.op2.tables.oef_forces.oef_forceObjects import RealPlateBilinearForce, RealRodForce # RealPlateForce

class OEF(object):
    def __init__(self):
        pass

    def _temperature_gradients_and_fluxes(self):
        (subcase_name, isubcase, transient, dt, analysis_code, is_sort1) = self._read_f06_subcase_header()
        headers = self.skip(2)
        #print "headers = %s" % (headers)
        data = self._read_gradient_fluxes_table()
        return
        #if isubcase in self.temperatureGrad:
            #self.temperatureGrad[isubcase].addData(data)
        #else:
            #self.temperatureGrad[isubcase] = TemperatureGradientObject(isubcase, data)
        #self.iSubcases.append(isubcase)

    def _read_gradient_fluxes_table(self):
        data = []
        Format = [int, str, float, float, float, float, float, float]
        while 1:
            line = self.infile.readline()[1:].rstrip('\r\n ')
            self.i += 1
            if 'PAGE' in line:
                return data
            sline = [line[0:15], line[15:24].strip(), line[24:44], line[44:61], line[61:78], line[78:95], line[95:112], line[112:129]]
            sline = self._parse_line_gradients_fluxes(sline, Format)
            data.append(sline)
        return data

    def _parse_line_gradients_fluxes(self, sline, Format):
        out = []
        for entry, iFormat in zip(sline, Format):
            if entry.strip() is '':
                out.append(0.0)
            else:
                #print "sline=|%r|\n entry=|%r| format=%r" % (sline, entry, iFormat)
                entry2 = iFormat(entry)
                out.append(entry2)
        return out


    def _forces_in_crod_elements(self):
        self._forces_in_rod_elements('CROD', 1, 'crod_force')

    def _forces_in_ctube_elements(self):
        self._forces_in_rod_elements('CTUBE', 3, 'ctube_force')

    def _forces_in_conrod_elements(self):
        self._forces_in_rod_elements('CONROD', 10, 'conrod_force')

    def _forces_in_rod_elements(self, element_name, element_type, result_name):
        """
            # 1-CROD
            # 3-CTUBE
            # 10-CONROD
        ::
                                                     F O R C E S   I N   R O D   E L E M E N T S     ( C R O D )
         ELEMENT           AXIAL                                     ELEMENT           AXIAL
           ID.             FORCE          TORQUE                       ID.             FORCE          TORQUE
             1       -7.007184E+02   0.0                                 2       -4.900904E+04   0.0
             3       -7.141140E+04   0.0
        """
        slot = getattr(self, result_name)
        #print(self.stored_lines)
        (subcase_name, isubcase, transient, dt, analysis_code, is_sort1) = self._read_f06_subcase_header()
        headers = self.skip(3)

        lines = []
        while 1:
            line = self.infile.readline().rstrip(' \n\r')
            self.i += 1
            if 'PAGE' in line:
                break
            lines.append(line)
            self.fatal_check(line)

        data = []
        for line in lines:
            eid, axial, torque = line[1:15], line[15:35], line[35:41]
            try:
                eid = int(eid)
                axial = float(axial)
                torque = float(torque)
            except ValueError:
                raise ValueError('line=%r' % line)
            data.append([eid, axial, torque])
            if line[41:].strip():
                eid, axial, torque = line[41:75], line[75:95], line[95:101]
                try:
                    eid = int(eid)
                    axial = float(axial)
                    torque = float(torque)
                except ValueError:
                    raise ValueError('line=%r' % line[35:])
                data.append([eid, axial, torque])
            #if line[101:]:
                #asdf

        data_code = {
            'analysis_code': analysis_code,
            'device_code': 1,
            'sort_code': 0,
            'sort_bits': [0, 0, 0],

            'table_name': 'OEF1X', # probably wrong
            'table_code': 5, # wrong
            'num_wide': 3,

            'format_code': 1,
            'element_name': element_name, 'element_type': element_type,
            'nonlinear_factor': dt,
            'data_names':['lsdvmn'],
            'lsdvmn': 1,
            }

        is_sort1 = True
        #ngrids = 4
        #print('isubcase =', isubcase)
        if isubcase not in slot:
            assert 'nonlinear_factor' in data_code
            slot[isubcase] = RealRodForce(data_code, is_sort1, isubcase, transient)
        slot[isubcase].add_f06_data(data, transient)
        self.iSubcases.append(isubcase)

    def _forces_in_cquad4s_bilinear(self):
        """
        ::

                                  F O R C E S   I N   Q U A D R I L A T E R A L   E L E M E N T S   ( Q U A D 4 )        OPTION = BILIN

            ELEMENT                    - MEMBRANE  FORCES -                      - BENDING   MOMENTS -            - TRANSVERSE SHEAR FORCES -
              ID       GRID-ID     FX            FY            FXY           MX            MY            MXY           QX            QY
                  1    CEN/4  0.0           0.0           0.0          -7.371223E+01 -4.023861E+02 -2.679984E+01  1.315875E+01 -7.356985E+01
                           1  0.0           0.0           0.0          -1.043592E+02 -3.888291E+02 -2.698050E+01  1.315875E+01 -7.356985E+01
                           2  0.0           0.0           0.0          -1.036512E+02 -4.152917E+02 -2.731157E+01  1.315875E+01 -7.356985E+01
                           8  0.0           0.0           0.0          -4.306526E+01 -4.159432E+02 -2.661917E+01  1.315875E+01 -7.356985E+01
                           7  0.0           0.0           0.0          -4.377329E+01 -3.894806E+02 -2.628810E+01  1.315875E+01 -7.356985E+01

        element_type = 33 b/c not bilinear
        """
        # composite
        element_name = 'CQUAD4'
        element_type = 144   # was listed as 95; has to be 144...i think...
        #print(self.stored_lines)
        (subcase_name, isubcase, transient, dt, analysis_code, is_sort1) = self._read_f06_subcase_header()
        headers = self.skip(3)

        lines = []
        while 1:
            line = self.infile.readline().rstrip('\n\r')
            self.i += 1
            if 'PAGE' in line:
                break
            lines.append(line)
            self.fatal_check(line)

        data = []
        for line in lines:
            eid, grid, fx, fy, fxy, mx, my, mxy, qx, qy = line[1:15], line[15:20], line[20:35], line[35:49], line[49:63], line[63:77], line[77:91], line[91:105], line[105:119], line[119:140]
            eid = eid.strip()
            grid = grid.strip()
            if eid:
                eid = int(eid)
            if 'C' not in grid:
                grid = int(grid)
            fx = float(fx)
            fy = float(fy)
            fxy = float(fxy)
            mx = float(mx)
            my = float(my)
            mxy = float(mxy)
            qx = float(qx)
            qy = float(qy)
            data.append([eid, grid, fx, fy, fxy, mx, my, mxy, qx, qy])

        data_code = {'analysis_code': analysis_code,
                    'device_code': 1,
                    'sort_code': 0,
                    'sort_bits': [0, 0, 0],

                    'table_name': 'OEF1X', # probably wrong
                    'table_code': 5, # wrong
                    'num_wide': 10,

                    'format_code': 1,
                    'element_name': element_name, 'element_type': element_type,
                    'nonlinear_factor': dt,
                    'data_names':['lsdvmn'],
                    'lsdvmn': 1,
                    }

        is_sort1 = True
        ngrids = 4
        #result_name = 'cquad4_composite_plate_force'
        #slot = self.cquad4_composite_plate_force
        result_name = 'cquad4_force'
        slot = self.cquad4_force


        if isubcase not in slot:
            assert 'nonlinear_factor' in data_code
            slot[isubcase] = RealPlateBilinearForce(data_code, is_sort1, isubcase, transient)
        slot[isubcase].add_f06_data(transient, data, element_name, ngrids)
        self.iSubcases.append(isubcase)

